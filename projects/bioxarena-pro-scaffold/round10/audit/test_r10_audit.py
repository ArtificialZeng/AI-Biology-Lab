#!/usr/bin/env python3
"""Unit & planted-truth verification suite for Round 10 Group C independent audit (A-T3, PREREG_C)."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from independent_r10_loader import (
    AuditLineData,
    enforce_r10_gate,
    parse_confirmatory_lines,
)
from independent_r10_metrics import (
    adjudicate_round10_group_c,
    evaluate_fct_perm_null,
    evaluate_hc3_noninferiority,
    evaluate_joint_residuals,
    evaluate_s7_confounder_controls,
)
from independent_r10_predictors import (
    align_source_and_target,
    build_predictions,
    fit_3d_shared_subspace,
)


def _make_planted_line(
    name: str,
    genes: list[str],
    targets: list[str],
    mu_lib: np.ndarray,
    fc_true: np.ndarray,
    n_cells: np.ndarray,
    seed: int,
    noise_sd: float = 0.015,
) -> AuditLineData:
    rng = np.random.default_rng(seed)
    n_t, n_g = len(targets), len(genes)
    mu_safe = np.maximum(mu_lib, 1e-8)
    clean_delta = np.log1p(mu_safe[None, :] * fc_true) - np.log1p(mu_safe)[None, :]
    ha = clean_delta + rng.normal(0.0, noise_sd, size=(n_t, n_g))
    hb = clean_delta + rng.normal(0.0, noise_sd, size=(n_t, n_g))
    deltas = 0.5 * (ha + hb)

    g_map = {g: i for i, g in enumerate(genes)}
    t_idx = np.asarray([g_map.get(t, -1) for t in targets], dtype=np.int64)
    kd = np.asarray([abs(float(deltas[i, t_idx[i]])) if t_idx[i] >= 0 else np.nan for i in range(n_t)], dtype=np.float64)
    ctrl_cells = rng.normal(0.0, 0.4, size=(120, n_g)) + np.log1p(mu_safe)[None, :]
    return AuditLineData(
        name=name,
        genes=genes,
        target_names=targets,
        target_indices=t_idx,
        deltas=deltas,
        ctrl_cells=ctrl_cells,
        ctrl_mean_libnorm=mu_lib.copy(),
        ctrl_mean_libnorm_full=mu_lib.copy(),
        repr_mask=np.ones(n_t, dtype=bool),
        half_a_deltas=ha,
        half_b_deltas=hb,
        half_valid=np.ones(n_t, dtype=bool),
        n_cells=n_cells.astype(np.float64),
        kd_depth_log1p=kd,
        metadata={"planted": True},
    )


class TestRound10Audit(unittest.TestCase):
    def test_blinding_gate_refuses_unfrozen_r10(self):
        with tempfile.TemporaryDirectory() as tmp:
            sdir = Path(tmp)
            (sdir / "PREREG_C_COMMIT.txt").write_text("dc4b4a5\n", encoding="utf-8")
            with self.assertRaises(RuntimeError):
                enforce_r10_gate(["dld1"], allow_r9_only=False, status_dir=sdir)
            with self.assertRaises(RuntimeError):
                enforce_r10_gate(["dld1"], allow_r9_only=True, status_dir=sdir)
            (sdir / "CODE_FREEZE_R10.txt").write_text("deadbeef\n", encoding="utf-8")
            prereg_h, freeze_h = enforce_r10_gate(["dld1"], allow_r9_only=False, status_dir=sdir)
            self.assertEqual(prereg_h, "dc4b4a5")
            self.assertEqual(freeze_h, "deadbeef")

    def test_p1_kconf_lt2_no_adjudication_and_unanimous_kconf_ge2(self):
        with tempfile.TemporaryDirectory() as tmp:
            cf = Path(tmp) / "CONFIRMATORY_LINES.txt"
            cf.write_text("# comment\nCONFIRMATORY_LINES=dld1\n", encoding="utf-8")
            self.assertEqual(parse_confirmatory_lines(cf), ["dld1"])
        v1 = adjudicate_round10_group_c({"dld1": {"hc1_pass": True, "hc2_pass": True, "hc3_pass": True}}, ["dld1"])
        self.assertEqual(v1["status"], "NO_CONFIRMATORY_ADJUDICATION")
        self.assertIn("确认性线不足", v1["primary_verdict"])

        v2_fail = adjudicate_round10_group_c(
            {
                "dld1": {"hc1_pass": True, "hc2_pass": True, "hc3_pass": True},
                "cd4t_rest": {"hc1_pass": True, "hc2_pass": False, "hc3_pass": True},
            },
            ["dld1", "cd4t_rest"],
        )
        self.assertEqual(v2_fail["status"], "ADJUDICATED")
        self.assertIn("NOT CONFIRMED", v2_fail["primary_verdict"])
        self.assertIn("PASS", v2_fail["hc3_verdict"])

    def test_planted_fct_beats_mismatch_and_perm100(self):
        rng = np.random.default_rng(77)
        n_g, n_t = 90, 45
        genes = [f"GENE{i}" for i in range(n_g)]
        targets = genes[:n_t]
        # Distinct source and target basal expression profiles create genuine FCT residual signal
        fc_true = np.exp(rng.normal(0.0, 0.35, size=(n_t, n_g)))
        mu_s1 = rng.uniform(0.5, 40.0, size=n_g)
        mu_s2 = rng.uniform(0.5, 40.0, size=n_g)
        mu_tgt = rng.uniform(0.5, 40.0, size=n_g)
        nc = rng.integers(30, 300, size=n_t)

        src1 = _make_planted_line("nadig_hepg2", genes, targets, mu_s1, fc_true, nc, seed=1)
        src2 = _make_planted_line("nadig_jurkat", genes, targets, mu_s2, fc_true, nc, seed=2)
        tgt = _make_planted_line("dld1", genes, targets, mu_tgt, fc_true, nc, seed=3)

        cohort = align_source_and_target([src1, src2], tgt, only_reproducible=True)
        preds, fc_bar, grand_mean = build_predictions(cohort, shrink=1e-3, include_sa3_anchors=True)
        res_fct = evaluate_joint_residuals(
            cohort.delta_obs,
            preds["fct"],
            preds["transferred_mean"],
            cohort.target_cols,
            mode="resp_top50",
            n_boot=250,
            seed=20261122,
            half_a_obs=cohort.half_a_obs,
            half_b_obs=cohort.half_b_obs,
            half_valid=cohort.half_valid,
        )
        res_mis1 = evaluate_joint_residuals(
            cohort.delta_obs,
            preds["fct_mismatch::nadig_hepg2"],
            preds["transferred_mean"],
            cohort.target_cols,
            mode="resp_top50",
            n_boot=250,
            seed=20261122,
            half_a_obs=cohort.half_a_obs,
            half_b_obs=cohort.half_b_obs,
            half_valid=cohort.half_valid,
        )
        perm = evaluate_fct_perm_null(
            cohort.delta_obs,
            preds["transferred_mean"],
            fc_bar,
            cohort.tgt_ctrl_mean_libnorm,
            cohort.target_cols,
            ceil_raw=res_fct["noise_ceiling_r"],
            ceil_orth=res_fct["noise_ceiling_r_orth"],
            n_perm=40,
        )
        self.assertGreater(res_fct["r_orth_cal_ci975_lo"], 0.6)
        self.assertGreater(res_fct["r_orth_cal_ci975_lo"], res_mis1["r_orth_cal_ci975_hi"])
        self.assertGreater(res_fct["r_orth_cal_ci975_lo"], perm["r_orth_cal"]["p975"])

        # R18 GlobalScale must collapse identically to 0 under r_orth_cal
        res_gs = evaluate_joint_residuals(
            cohort.delta_obs,
            preds["transfer_mean_global_scale"],
            preds["transferred_mean"],
            cohort.target_cols,
            mode="resp_top50",
            n_boot=50,
            seed=20261122,
            half_a_obs=cohort.half_a_obs,
            half_b_obs=cohort.half_b_obs,
            half_valid=cohort.half_valid,
        )
        self.assertLess(abs(res_gs["r_orth_cal_est"]), 1e-10)

        # H-C3 should pass & find strict superiority on FCT vs TMean
        hc3 = evaluate_hc3_noninferiority(
            cohort.delta_obs,
            preds["fct"],
            preds["transferred_mean"],
            grand_mean,
            cohort.target_cols,
            n_boot=150,
            seed=20261122,
        )
        self.assertTrue(hc3["hc3_pass_line"])
        self.assertEqual(hc3["n_inferior_conventions"], 0)
        self.assertGreater(hc3["ctrl_de50_r"]["diff_mean"], 0.002)
        self.assertGreater(hc3["ctrl_de50_r"]["diff_ci975_lo"], 0.0)

        # S7 confounder checks + centered OLS invariance + 3 tertiles
        s7 = evaluate_s7_confounder_controls(
            cohort.eval_targets,
            cohort.delta_obs,
            preds["fct"],
            preds["transferred_mean"],
            fc_bar,
            cohort.tgt_ctrl_mean_libnorm,
            cohort.source_names,
            cohort.source_ctrl_mean_libnorm,
            cohort.target_cols,
            cohort.n_cells,
            cohort.kd_depth_log1p,
            half_a_obs=cohort.half_a_obs,
            half_b_obs=cohort.half_b_obs,
            half_valid=cohort.half_valid,
            n_boot=120,
            n_perm=20,
            seed=20261122,
        )
        self.assertLess(s7["ols_centered"]["uncentered_eval_invariance_err"], 1e-10)
        self.assertEqual(set(s7["tertiles"].keys()), {"T1_low_ncells", "T2_mid_ncells", "T3_high_ncells"})
        self.assertTrue(s7["tertiles"]["T1_low_ncells"]["passes_hc1"])
        self.assertTrue(s7["tertiles"]["T1_low_ncells"]["passes_hc2"])

    def test_fit_3d_shared_subspace_rank(self):
        rng = np.random.default_rng(11)
        mat = rng.normal(size=(25, 40)) + 0.5
        q3 = fit_3d_shared_subspace({"s1": {f"t{i}": mat[i] for i in range(mat.shape[0])}}, n_genes=40)
        self.assertEqual(q3.shape, (40, 3))
        np.testing.assert_allclose(q3.T @ q3, np.eye(3), atol=1e-10)

    def test_precomputed_cov_subset_equivalence(self):
        from independent_r10_predictors import sample_covariance
        rng = np.random.default_rng(19)
        n_g, n_t = 40, 20
        genes = [f"G{i}" for i in range(n_g)]
        targets = genes[:n_t]
        fc_true = np.exp(rng.normal(0.0, 0.25, size=(n_t, n_g)))
        mu_s = rng.uniform(1.0, 20.0, size=n_g)
        mu_t = rng.uniform(1.0, 20.0, size=n_g)
        nc = np.full(n_t, 50)
        s1 = _make_planted_line("s1", genes, targets, mu_s, fc_true, nc, seed=1)
        # Drop 5 trailing genes in target to verify G_eval(t) covariance column slicing matches
        sub_genes = genes[:35]
        t_sub = _make_planted_line("t1", sub_genes, sub_genes[:15], mu_t[:35], fc_true[:15, :35], nc[:15], seed=2)
        cohort = align_source_and_target([s1], t_sub, only_reproducible=True)
        p_direct, _, _ = build_predictions(cohort, shrink=0.0, include_sa3_anchors=True)
        p_cached, _, _ = build_predictions(
            cohort,
            shrink=0.0,
            include_sa3_anchors=True,
            precomputed_covs={"s1": sample_covariance(s1.ctrl_cells, 0.0), "t1": sample_covariance(t_sub.ctrl_cells, 0.0)},
            precomputed_gene_orders={"s1": genes, "t1": sub_genes},
        )
        for k in p_direct:
            np.testing.assert_allclose(p_direct[k], p_cached[k], atol=1e-12)


if __name__ == "__main__":
    unittest.main(verbosity=2)
