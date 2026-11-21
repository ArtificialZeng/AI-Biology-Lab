#!/usr/bin/env python3
"""Unit and planted-truth end-to-end tests for independent_predictors + independent_stats."""

from __future__ import annotations

import unittest

import numpy as np

from independent_predictors import (
    compute_control_covariance,
    predict_zero_shot_suite,
    src_corr_target_var_covariance,
)
from independent_stats import (
    adjudicate_team_a,
    bootstrap_paired_mean_diff,
    bootstrap_perturbation_stat,
    compute_de50_r_rows,
    compute_disc_raw_rows,
    disattenuated_correlation,
    noise_ceiling_from_splithalf,
    ols_slope,
    pearson_1d,
)


def _simulate_line_data(seed: int, n_genes: int = 80, n_ctrl: int = 400):
    """Generate synthetic 5-line Perturb-seq data with line-specific covariance signal."""
    rng = np.random.default_rng(seed)
    lines = ["L1", "L2", "L3", "L4", "L5"]
    targets = [f"G{i}" for i in range(40)]
    target_to_col = {g: i for i, g in enumerate(targets)}

    # Shared perturbation program across all lines
    shared_deltas = rng.normal(0.0, 0.8, size=(len(targets), n_genes))

    ctrl_covs = {}
    deltas_by_line = {}
    half_a_by_line = {}
    half_b_by_line = {}

    for li, line in enumerate(lines):
        # Distinct block-structured latent regulatory loadings per cell line
        loadings = rng.normal(0.0, 1.0, size=(n_genes, 12))
        true_cov = loadings @ loadings.T + 0.3 * np.eye(n_genes)
        ctrl_cells = rng.multivariate_normal(np.zeros(n_genes), true_cov, size=n_ctrl)
        ctrl_covs[line] = compute_control_covariance(ctrl_cells, shrink=1e-3)

        # True perturbation deltas = shared + 0.5 * true_cov[:, g]/true_cov[g,g] + noise
        line_deltas = {}
        ha_line = {}
        hb_line = {}
        for ti, tgt in enumerate(targets):
            g_col = target_to_col[tgt]
            ls_true = 0.5 * (true_cov[:, g_col] / true_cov[g_col, g_col])
            clean = shared_deltas[ti] + ls_true
            noise_a = rng.normal(0.0, 0.15, size=n_genes)
            noise_b = rng.normal(0.0, 0.15, size=n_genes)
            ha = clean + noise_a
            hb = clean + noise_b
            ha_line[tgt] = ha
            hb_line[tgt] = hb
            line_deltas[tgt] = 0.5 * (ha + hb)
        deltas_by_line[line] = line_deltas
        half_a_by_line[line] = ha_line
        half_b_by_line[line] = hb_line

    return lines, targets, target_to_col, ctrl_covs, deltas_by_line, half_a_by_line, half_b_by_line


class TestIndependentPredictors(unittest.TestCase):
    """Test zero-leakage invariance and planted-truth end-to-end audit pipeline."""

    def test_zero_leakage_invariance(self):
        lines, targets, t2c, covs, deltas, _, _ = _simulate_line_data(101)
        tgt_line = "L1"
        src_covs = {l: covs[l] for l in lines if l != tgt_line}
        src_deltas = {l: deltas[l] for l in lines if l != tgt_line}
        mu_lib = {l: np.full(covs[tgt_line].shape[0], 2.0 + i) for i, l in enumerate(lines)}
        dlib = {l: {t: np.ones(covs[tgt_line].shape[0]) * 0.5 for t in targets} for l in src_covs}
        preds1 = predict_zero_shot_suite(
            targets,
            t2c,
            covs[tgt_line],
            src_covs,
            src_deltas,
            mismatch_cov=covs["L2"],
            target_ctrl_mean_libnorm=mu_lib[tgt_line],
            source_ctrl_mean_libnorm={l: mu_lib[l] for l in src_covs},
            source_deltas_libnorm_by_target=dlib,
        )
        # Predictors never even take target deltas as an argument -> structural zero leakage.
        # Verify determinism with identical inputs:
        preds2 = predict_zero_shot_suite(
            targets,
            t2c,
            covs[tgt_line],
            src_covs,
            src_deltas,
            mismatch_cov=covs["L2"],
            target_ctrl_mean_libnorm=mu_lib[tgt_line],
            source_ctrl_mean_libnorm={l: mu_lib[l] for l in src_covs},
            source_deltas_libnorm_by_target=dlib,
        )
        for k in preds1:
            np.testing.assert_array_equal(preds1[k], preds2[k])
        # R3/R4/R12/R14/R15 predictors emitted cleanly
        for sl in src_covs:
            self.assertIn(f"mismatch_cov_hybrid::{sl}", preds1)
            self.assertIn(f"tmvs_mismatch::{sl}", preds1)
        self.assertIn("src_corr_tgt_var_hybrid", preds1)
        self.assertIn("transfer_mean_var_scaled", preds1)
        self.assertIn("transfer_mean_global_scale", preds1)
        self.assertIn("tmvs_shuffled_var", preds1)
        self.assertIn("hybrid_on_tmvs", preds1)
        self.assertIn("fct", preds1)
        self.assertTrue(np.all(np.isfinite(preds1["fct"])))

        # R4(a) matrix preserves target marginal variances exactly while averaging source correlations
        sctv = src_corr_target_var_covariance(covs[tgt_line], list(src_covs.values()))
        np.testing.assert_allclose(np.diag(sctv), np.diag(covs[tgt_line]), atol=1e-10)

    def test_planted_signal_passes_all_hypotheses(self):
        lines, targets, t2c, covs, deltas, ha, hb = _simulate_line_data(2026)
        n_genes = len(next(iter(deltas["L1"].values())))
        t_cols = [t2c[t] for t in targets]
        line_audit_rows = []

        for li, tgt_line in enumerate(lines):
            src_lines = [l for l in lines if l != tgt_line]
            src_covs = {l: covs[l] for l in src_lines}
            src_deltas = {l: deltas[l] for l in src_lines}
            mismatch = covs[src_lines[0]]

            preds = predict_zero_shot_suite(
                targets, t2c, covs[tgt_line], src_covs, src_deltas, mismatch_cov=mismatch
            )
            y_true = np.stack([deltas[tgt_line][t] for t in targets])
            y_ha = np.stack([ha[tgt_line][t] for t in targets])
            y_hb = np.stack([hb[tgt_line][t] for t in targets])
            t_mean = preds["transferred_mean"]

            # Line-specific residuals
            true_ls = y_true - t_mean
            hyb_ls = preds["hybrid"] - t_mean
            mis_ls = preds["mismatch_cov_hybrid"] - t_mean
            shuf_ls = preds["shuffled_gene_hybrid"] - t_mean
            ha_ls = y_ha - t_mean
            hb_ls = y_hb - t_mean

            # Mask each perturbation's own target gene for evaluation
            mask_mat = np.ones((len(targets), n_genes), dtype=bool)
            for i, gc in enumerate(t_cols):
                mask_mat[i, gc] = False

            def _slope_ci(pred_ls_mat):
                return bootstrap_perturbation_stat(
                    len(targets),
                    lambda idx: ols_slope(
                        true_ls[idx][mask_mat[idx]],
                        pred_ls_mat[idx][mask_mat[idx]],
                        fit_intercept=True,
                    ),
                    n_boot=500,
                    seed=42 + li,
                )

            ci_hyb_slope = _slope_ci(hyb_ls)
            ci_mis_slope = _slope_ci(mis_ls)
            ci_shuf_slope = _slope_ci(shuf_ls)

            # Split-half ceiling and noise-calibrated correlation on line-specific residual
            r_ab = pearson_1d(ha_ls[mask_mat], hb_ls[mask_mat])
            ceil = noise_ceiling_from_splithalf(r_ab)

            ci_cal_corr = bootstrap_perturbation_stat(
                len(targets),
                lambda idx: disattenuated_correlation(
                    pearson_1d(true_ls[idx][mask_mat[idx]], hyb_ls[idx][mask_mat[idx]]),
                    ceil,
                ),
                n_boot=500,
                seed=99 + li,
            )

            # H-A4 metrics
            disc_hyb = compute_disc_raw_rows(y_true, preds["hybrid"], masked_cols=t_cols)
            disc_tm = compute_disc_raw_rows(y_true, t_mean, masked_cols=t_cols)
            de50_hyb = compute_de50_r_rows(y_true, preds["hybrid"], target_col_idx=t_cols)
            de50_tm = compute_de50_r_rows(y_true, t_mean, target_col_idx=t_cols)
            ci_disc = bootstrap_paired_mean_diff(disc_hyb, disc_tm, n_boot=500, seed=7 + li)
            ci_de50 = bootstrap_paired_mean_diff(de50_hyb, de50_tm, n_boot=500, seed=17 + li)

            # Shuffled covariance on hybrid still must lie far below matched hybrid (H-A3 specificity),
            # while pure CIPHER with shuffled covariance must have ~0 slope.
            self.assertLess(ci_shuf_slope.ci95_hi, ci_hyb_slope.ci95_lo)
            self.assertLess(ci_mis_slope.ci95_hi, ci_hyb_slope.ci95_lo)
            line_audit_rows.append(
                {
                    "line": tgt_line,
                    "primary_slope_ci95_lo": ci_hyb_slope.ci95_lo,
                    "noise_cal_corr_ci95_lo": ci_cal_corr.ci95_lo,
                    "mismatch_slope_ci95_hi": ci_mis_slope.ci95_hi,
                    "disc_raw_diff_ci95_hi": ci_disc.ci95_hi,
                    "de50_r_diff_ci95_hi": ci_de50.ci95_hi,
                }
            )

        verdict = adjudicate_team_a(line_audit_rows)
        self.assertEqual(verdict["overall_verdict"], "CONFIRMED")


if __name__ == "__main__":
    unittest.main()
