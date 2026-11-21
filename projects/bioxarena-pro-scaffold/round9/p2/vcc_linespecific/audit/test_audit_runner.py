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

            res_null = audit_loo_suite(shared_dir / "null", lines, n_boot=300, seed=10)
            self.assertEqual(res_null["verdict"]["overall_verdict"], "未确认 (NOT CONFIRMED)")


if __name__ == "__main__":
    unittest.main()
