#!/usr/bin/env python3
"""Unit tests for independent Round 10-E post-hoc baseline ablation audit (A-T3)."""

import sys
import unittest
from pathlib import Path

import numpy as np

THIS_DIR = Path(__file__).resolve().parent
AUDIT_DIR = THIS_DIR.parent
for p in (str(THIS_DIR), str(AUDIT_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from independent_r10_metrics import evaluate_joint_residuals, evaluate_paired_cal_diff  # noqa: E402
from independent_r10_predictors import AlignedCohort, fct_from_fc_bar  # noqa: E402
from independent_r10e_ablation import (  # noqa: E402
    build_e1_constant_median_pred,
    build_e6_log_fc_pred,
    build_e7_linear_fc_pred,
    evaluate_line_ablation,
)


class TestR10EAblationUnit(unittest.TestCase):
    def test_e1_e6_e7_mathematical_identities_and_limits(self) -> None:
        rng = np.random.default_rng(7)
        fc_bar = np.exp(rng.normal(0.0, 0.5, size=(12, 80)))
        mu_t = rng.uniform(0.1, 10.0, size=80)

        # E1 depends strictly on the median of max(mu_t, 0) and is invariant to gene permutation of mu_t
        pred_e1, c_med = build_e1_constant_median_pred(fc_bar, mu_t)
        self.assertAlmostEqual(c_med, float(np.median(mu_t)), places=14)
        pred_e1_perm, _ = build_e1_constant_median_pred(fc_bar, rng.permutation(mu_t))
        self.assertLess(float(np.max(np.abs(pred_e1 - pred_e1_perm))), 1e-15)

        # E6 = ln(FC_bar) equals the large-c asymptotic limit of FCT(c)
        pred_e6 = build_e6_log_fc_pred(fc_bar)
        fct_huge_c = fct_from_fc_bar(fc_bar, np.full_like(mu_t, 1e10))
        self.assertLess(float(np.max(np.abs(pred_e6 - fct_huge_c))), 1e-9)

        # E7 = FC_bar - 1 equals small-c FCT(c)/c limit and preserves exact r_orth_cal invariance under scaling
        pred_e7 = build_e7_linear_fc_pred(fc_bar)
        tiny_c = 1e-6
        fct_tiny_c = fct_from_fc_bar(fc_bar, np.full_like(mu_t, tiny_c)) / tiny_c
        self.assertLess(float(np.max(np.abs(pred_e7 - fct_tiny_c))), 1e-5)

        t_mean = rng.normal(0.0, 0.4, size=(12, 80))
        delta_obs = t_mean + 0.5 * (pred_e7 - t_mean) + rng.normal(0.0, 0.05, size=(12, 80))
        tcols = np.arange(12, dtype=np.int64)
        ha = delta_obs + rng.normal(0.0, 0.02, size=(12, 80))
        hb = delta_obs + rng.normal(0.0, 0.02, size=(12, 80))
        res_e7 = evaluate_joint_residuals(delta_obs, pred_e7, t_mean, tcols, mode="resp_top50", n_boot=50, seed=11, half_a_obs=ha, half_b_obs=hb)
        res_fct_tiny = evaluate_joint_residuals(delta_obs, tiny_c * pred_e7, t_mean, tcols, mode="resp_top50", n_boot=50, seed=11, half_a_obs=ha, half_b_obs=hb)
        self.assertAlmostEqual(res_e7["r_orth_cal_est"], res_fct_tiny["r_orth_cal_est"], places=4)

    def test_paired_fct_minus_e1_exact_zero_on_constant_ctrl_mean(self) -> None:
        rng = np.random.default_rng(19)
        n_t, n_g = 20, 100
        mu_const = np.full(n_g, 2.5, dtype=np.float64)
        src_mu = {"s1": rng.uniform(0.2, 5.0, size=n_g), "s2": rng.uniform(0.2, 5.0, size=n_g)}
        src_d = {
            "s1": {f"G{i}": rng.normal(0.0, 0.3, size=n_g) for i in range(n_t)},
            "s2": {f"G{i}": rng.normal(0.0, 0.3, size=n_g) for i in range(n_t)},
        }
        delta_obs = (
            np.stack([0.5 * (src_d["s1"][f"G{i}"] + src_d["s2"][f"G{i}"]) for i in range(n_t)], axis=0)
            + rng.normal(0.0, 0.25, size=(n_t, n_g))
        )
        cohort = AlignedCohort(
            target_line="synthetic_const",
            eval_genes=[f"G{i}" for i in range(n_g)],
            eval_targets=[f"G{i}" for i in range(n_t)],
            target_cols=np.arange(n_t, dtype=np.int64),
            delta_obs=delta_obs,
            n_cells=np.full(n_t, 50.0),
            kd_depth_log1p=np.ones(n_t),
            half_a_obs=delta_obs + rng.normal(0.0, 0.05, size=(n_t, n_g)),
            half_b_obs=delta_obs + rng.normal(0.0, 0.05, size=(n_t, n_g)),
            half_valid=np.ones(n_t, dtype=bool),
            tgt_ctrl_cells=np.ones((10, n_g)),
            tgt_ctrl_mean_libnorm=mu_const,
            tgt_ctrl_mean_libnorm_full=mu_const,
            source_names=["s1", "s2"],
            source_deltas_by_line=src_d,
            source_ctrl_cells={"s1": np.ones((10, n_g)), "s2": np.ones((10, n_g))},
            source_ctrl_mean_libnorm=src_mu,
        )
        m_rows, p_rows = evaluate_line_ablation(cohort, masks=("resp_top50", "all_nontarget"), n_boot=100, seed=20261122)
        self.assertEqual(len(m_rows), 8)
        self.assertEqual(len(p_rows), 2)
        for pr in p_rows:
            self.assertAlmostEqual(float(pr["diff_est"]), 0.0, places=15)
            self.assertAlmostEqual(float(pr["diff_ci975_lo"]), 0.0, places=15)
            self.assertAlmostEqual(float(pr["diff_ci975_hi"]), 0.0, places=15)
            self.assertAlmostEqual(float(pr["e1_pct_of_fct"]), 100.0, places=12)

    def test_paired_diff_equals_difference_of_independent_point_estimates(self) -> None:
        rng = np.random.default_rng(123)
        n_t, n_g = 30, 120
        mu_var = rng.lognormal(0.0, 1.0, size=n_g)
        src_mu = {"s1": rng.lognormal(0.0, 1.0, size=n_g), "s2": rng.lognormal(0.0, 1.0, size=n_g)}
        src_d = {
            "s1": {f"G{i}": rng.normal(0.0, 0.4, size=n_g) for i in range(n_t)},
            "s2": {f"G{i}": rng.normal(0.0, 0.4, size=n_g) for i in range(n_t)},
        }
        fc = np.exp(0.5 * rng.normal(0.0, 0.4, size=(n_t, n_g)))
        delta_obs = np.log1p(mu_var[None, :] * fc) - np.log1p(mu_var)[None, :] + rng.normal(0.0, 0.05, size=(n_t, n_g))
        cohort = AlignedCohort(
            target_line="synthetic_var",
            eval_genes=[f"G{i}" for i in range(n_g)],
            eval_targets=[f"G{i}" for i in range(n_t)],
            target_cols=np.arange(n_t, dtype=np.int64),
            delta_obs=delta_obs,
            n_cells=np.full(n_t, 60.0),
            kd_depth_log1p=np.ones(n_t),
            half_a_obs=delta_obs + rng.normal(0.0, 0.05, size=(n_t, n_g)),
            half_b_obs=delta_obs + rng.normal(0.0, 0.05, size=(n_t, n_g)),
            half_valid=np.ones(n_t, dtype=bool),
            tgt_ctrl_cells=np.ones((10, n_g)),
            tgt_ctrl_mean_libnorm=mu_var,
            tgt_ctrl_mean_libnorm_full=mu_var,
            source_names=["s1", "s2"],
            source_deltas_by_line=src_d,
            source_ctrl_cells={"s1": np.ones((10, n_g)), "s2": np.ones((10, n_g))},
            source_ctrl_mean_libnorm=src_mu,
        )
        m_rows, p_rows = evaluate_line_ablation(cohort, masks=("resp_top50",), n_boot=200, seed=20261122)
        fct_row = next(r for r in m_rows if r["predictor"] == "FCT")
        e1_row = next(r for r in m_rows if r["predictor"] == "E1")
        pr = p_rows[0]
        self.assertAlmostEqual(float(pr["diff_est"]), float(fct_row["r_orth_cal"]) - float(e1_row["r_orth_cal"]), places=14)
        self.assertGreater(float(pr["diff_est"]), 0.0)


if __name__ == "__main__":
    unittest.main()
