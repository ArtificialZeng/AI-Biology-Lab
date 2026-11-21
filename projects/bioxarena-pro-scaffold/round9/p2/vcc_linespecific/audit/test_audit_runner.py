#!/usr/bin/env python3
"""End-to-end tests for audit_runner.py (Rule 3 gate + synthetic shared dir audit)."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from audit_runner import audit_loo_suite, check_prereg_gate


class TestAuditRunner(unittest.TestCase):
    """Verify PREREG_COMMIT gate and end-to-end synthetic shared-dir audit."""

    def test_prereg_gate_blocks_without_commit(self):
        with tempfile.TemporaryDirectory() as td:
            sdir = Path(td)
            with self.assertRaises(PermissionError):
                check_prereg_gate(sdir, ["K562", "RPE1", "HepG2"], dev_only=False)
            # --dev-only rejects non-Replogle lines
            with self.assertRaises(PermissionError):
                check_prereg_gate(sdir, ["K562", "HepG2"], dev_only=True)
            # --dev-only permits Replogle K562/RPE1
            check_prereg_gate(sdir, ["K562", "RPE1"], dev_only=True)
            # With PREREG_COMMIT.txt present, holdout lines allowed
            (sdir / "PREREG_COMMIT.txt").write_text("abc1234\n")
            check_prereg_gate(sdir, ["K562", "RPE1", "HepG2"], dev_only=False)

    def test_synthetic_shared_dir_audit_confirmed_and_null_not_confirmed(self):
        rng = np.random.default_rng(77)
        lines = ["L1", "L2", "L3", "L4", "L5"]
        n_genes = 70
        n_targets = 35
        genes = [f"G{i}" for i in range(n_genes)]
        targets = genes[:n_targets]
        shared_prog = rng.normal(0.0, 0.8, size=(n_targets, n_genes))

        with tempfile.TemporaryDirectory() as td:
            shared_dir = Path(td)
            (shared_dir / "COMMON_GENES.txt").write_text("\n".join(genes) + "\n")
            for mode in ("signal", "null"):
                sub_dir = shared_dir / mode
                sub_dir.mkdir(parents=True, exist_ok=True)
                (sub_dir / "COMMON_GENES.txt").write_text("\n".join(genes) + "\n")
                for line in lines:
                    ldir = sub_dir / "lines" / line
                    ldir.mkdir(parents=True, exist_ok=True)
                    (ldir / "genes.txt").write_text("\n".join(genes) + "\n")
                    (ldir / "targets.txt").write_text("\n".join(targets) + "\n")
                    load = rng.normal(0.0, 1.0, size=(n_genes, 10))
                    cov = load @ load.T + 0.25 * np.eye(n_genes)
                    ctrl = rng.multivariate_normal(np.zeros(n_genes), cov, size=400)
                    np.savez(ldir / "ctrl_cells.npz", X=ctrl.astype(np.float32))
                    ha = np.zeros((n_targets, n_genes), dtype=np.float32)
                    hb = np.zeros((n_targets, n_genes), dtype=np.float32)
                    for i in range(n_targets):
                        ls = 0.5 * (cov[:, i] / cov[i, i]) if mode == "signal" else np.zeros(n_genes)
                        ha[i] = shared_prog[i] + ls + rng.normal(0.0, 0.12, size=n_genes)
                        hb[i] = shared_prog[i] + ls + rng.normal(0.0, 0.12, size=n_genes)
                    delta = 0.5 * (ha + hb)
                    np.savez(ldir / "pseudobulk_delta.npz", delta_zscore=delta, half_a=ha, half_b=hb)

            res_sig = audit_loo_suite(shared_dir / "signal", lines, n_boot=300, seed=10)
            self.assertEqual(res_sig["verdict"]["overall_verdict"], "CONFIRMED")
            row0 = res_sig["per_line_df"].iloc[0]
            self.assertIn("max_mismatch_ci95_hi", res_sig["per_line_df"].columns)
            self.assertIn("max_mismatch_ci975_hi", res_sig["per_line_df"].columns)
            self.assertIn("sctv_slope", res_sig["per_line_df"].columns)
            self.assertIn("tmvs_cal_r", res_sig["per_line_df"].columns)
            self.assertIn("tmvs_max_ctrl_cal_r_ci975_hi", res_sig["per_line_df"].columns)
            self.assertIn("r14_hotmvs_minus_tmvs_cal_r", res_sig["per_line_df"].columns)
            self.assertIn("s5_tmvs_cal_r", res_sig["per_line_df"].columns)
            self.assertIn("rho_ab_all_nontarget", res_sig["per_line_df"].columns)
            self.assertIn("r_cal_sq", res_sig["per_line_df"].columns)
            self.assertIn("slope_deming", res_sig["per_line_df"].columns)
            self.assertIn("ctrl_disc_ownmasked_hybrid", res_sig["per_line_df"].columns)
            self.assertEqual(len(res_sig["detailed"][row0["line"]]["mismatch_by_source"]), 4)
            # 4 source mismatches + 1 shuffled_var control = 5 Group B controls
            self.assertEqual(len(res_sig["detailed"][row0["line"]]["tmvs_ctrls"]), 5)
            self.assertIn("r4_wording_governance", res_sig["verdict"])
            self.assertIn("r14_mechanism_governance", res_sig["verdict"])
            self.assertIn("Group_B", res_sig["verdict"])

            # Confirmatory subset (L1..L3) excludes null lines cleanly under R2/R10
            res_conf = audit_loo_suite(
                shared_dir / "signal", lines, confirmatory_lines=["L1", "L2", "L3"], n_boot=300, seed=10
            )
            self.assertEqual(res_conf["verdict"]["n_lines"], 3)
            self.assertEqual(res_conf["verdict"]["overall_verdict"], "CONFIRMED")

            res_null = audit_loo_suite(shared_dir / "null", lines, n_boot=300, seed=10)
            self.assertEqual(res_null["verdict"]["overall_verdict"], "未确认 (NOT CONFIRMED)")


if __name__ == "__main__":
    unittest.main()
