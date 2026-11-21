#!/usr/bin/env python3
"""Unit tests for independent_stats.py (A-T3 statistical audit primitives)."""

from __future__ import annotations

import math
import unittest

import numpy as np

from independent_stats import (
    adjudicate_team_a,
    bootstrap_paired_mean_diff,
    bootstrap_perturbation_stat,
    ci_strictly_below,
    compute_de50_r_rows,
    compute_disc_raw_rows,
    disattenuated_correlation,
    noise_ceiling_from_splithalf,
    ols_slope,
    pearson_1d,
    spearman_brown_reliability,
)


class TestIndependentStats(unittest.TestCase):
    """Verify exactness of independent statistical and metric computations."""

    def test_ols_slope_and_bootstrap_recovery(self):
        rng = np.random.default_rng(7)
        x = rng.normal(0.0, 1.0, size=(120, 200))
        y = 0.25 * x + rng.normal(0.0, 0.2, size=(120, 200))

        def pooled_slope(idx: np.ndarray) -> float:
            return ols_slope(x[idx].ravel(), y[idx].ravel(), fit_intercept=True)

        ci = bootstrap_perturbation_stat(120, pooled_slope, n_boot=1000, seed=11)
        self.assertGreater(ci.estimate, 0.22)
        self.assertLess(ci.estimate, 0.28)
        self.assertLess(ci.ci95_lo, 0.25)
        self.assertGreater(ci.ci95_hi, 0.25)

    def test_origin_slope_and_null_ci(self):
        rng = np.random.default_rng(9)
        x = rng.normal(0.0, 1.0, size=(80, 100))
        y_null = rng.normal(0.0, 1.0, size=(80, 100))
        ci_null = bootstrap_perturbation_stat(
            80,
            lambda idx: ols_slope(x[idx].ravel(), y_null[idx].ravel(), fit_intercept=False),
            n_boot=1000,
            seed=12,
        )
        self.assertLess(ci_null.ci95_lo, 0.0)
        self.assertGreater(ci_null.ci95_hi, 0.0)

    def test_disc_raw_exact_and_ties(self):
        rng = np.random.default_rng(19)
        t_mat = rng.normal(0.0, 1.0, size=(20, 150))
        # Exact match -> 1.0 everywhere
        d_exact = compute_disc_raw_rows(t_mat, t_mat, masked_cols=[0, 1])
        np.testing.assert_allclose(d_exact, 1.0, atol=1e-12)

        # Constant zero predictions -> 0.5 everywhere
        d_zero = compute_disc_raw_rows(t_mat, np.zeros_like(t_mat))
        np.testing.assert_allclose(d_zero, 0.5, atol=1e-12)

        # Constant identical nonzero predictions -> permutation of ranks -> mean == 0.5
        const_row = rng.normal(0.0, 1.0, size=(1, 150))
        p_const = np.repeat(const_row, 20, axis=0)
        d_const = compute_disc_raw_rows(t_mat, p_const)
        self.assertAlmostEqual(float(np.mean(d_const)), 0.5, places=12)

        # Identical true vectors across two rows -> exact within-row tie -> 0.5
        t_tied = np.repeat(t_mat[:1], 5, axis=0)
        p_tied = np.repeat(t_mat[:1], 5, axis=0)
        d_tied = compute_disc_raw_rows(t_tied, p_tied)
        np.testing.assert_allclose(d_tied, 0.5, atol=1e-12)

    def test_de50_r_masks_target_gene(self):
        rng = np.random.default_rng(23)
        t_mat = rng.normal(0.0, 0.5, size=(5, 100))
        p_mat = t_mat.copy()
        # Make target column huge in true and inverted in pred; masking must ignore it
        target_cols = [10, 20, 30, 40, 50]
        for i, tc in enumerate(target_cols):
            t_mat[i, tc] = -100.0
            p_mat[i, tc] = +100.0
        scores = compute_de50_r_rows(t_mat, p_mat, target_col_idx=target_cols, top_k=50)
        np.testing.assert_allclose(scores, 1.0, atol=1e-12)

    def test_spearman_brown_noise_ceiling(self):
        rng = np.random.default_rng(31)
        signal = rng.normal(0.0, 1.0, size=5000)
        # Two halves with unit noise variance -> full sample has 0.5 noise variance -> SNR_full = 2
        half_a = signal + rng.normal(0.0, 1.0, size=5000)
        half_b = signal + rng.normal(0.0, 1.0, size=5000)
        full_obs = 0.5 * (half_a + half_b)
        r_ab = pearson_1d(half_a, half_b)
        ceil = noise_ceiling_from_splithalf(r_ab)
        r_obs = pearson_1d(signal, full_obs)
        r_cal = disattenuated_correlation(r_obs, ceil)
        # Expected ceiling = sqrt(2/3) ~ 0.8165, disattenuated ~ 1.0
        self.assertAlmostEqual(spearman_brown_reliability(0.5), 2.0 / 3.0, places=6)
        self.assertGreater(r_cal, 0.95)
        self.assertLess(r_cal, 1.05)

        self.assertTrue(math.isnan(disattenuated_correlation(0.2, 0.0)))

    def test_paired_bootstrap_and_adjudication(self):
        a = np.linspace(0.6, 0.8, 50)
        b = np.linspace(0.5, 0.7, 50)
        ci = bootstrap_paired_mean_diff(a, b, n_boot=1000, seed=5)
        self.assertAlmostEqual(ci.estimate, 0.1, places=6)
        self.assertGreater(ci.ci95_lo, 0.09)

        # Non-overlapping CIs
        c_low = bootstrap_paired_mean_diff(np.zeros(20), np.zeros(20) - 0.02, n_boot=500)
        c_high = bootstrap_paired_mean_diff(np.ones(20) * 0.25, np.zeros(20), n_boot=500)
        self.assertTrue(ci_strictly_below(c_low, c_high))

        # Confirmed case (4/5 pass H1-H3, 1/5 inferior on H4)
        rows_pass = [
            {"line": "L1", "primary_slope_ci95_lo": 0.15, "noise_cal_corr_ci95_lo": 0.10, "mismatch_slope_ci95_hi": 0.05, "disc_raw_diff_ci95_hi": 0.01, "de50_r_diff_ci95_hi": 0.01},
            {"line": "L2", "primary_slope_ci95_lo": 0.18, "noise_cal_corr_ci95_lo": 0.12, "mismatch_slope_ci95_hi": 0.02, "disc_raw_diff_ci95_hi": 0.02, "de50_r_diff_ci95_hi": 0.02},
            {"line": "L3", "primary_slope_ci95_lo": 0.14, "noise_cal_corr_ci95_lo": 0.08, "mismatch_slope_ci95_hi": 0.01, "disc_raw_diff_ci95_hi": 0.03, "de50_r_diff_ci95_hi": -0.01},
            {"line": "L4", "primary_slope_ci95_lo": 0.20, "noise_cal_corr_ci95_lo": 0.05, "mismatch_slope_ci95_hi": 0.08, "disc_raw_diff_ci95_hi": 0.01, "de50_r_diff_ci95_hi": 0.02},
            {"line": "L5", "primary_slope_ci95_lo": 0.09, "noise_cal_corr_ci95_lo": -0.01, "mismatch_slope_ci95_hi": 0.11, "disc_raw_diff_ci95_hi": -0.01, "de50_r_diff_ci95_hi": 0.01},
        ]
        verdict = adjudicate_team_a(rows_pass)
        self.assertEqual(verdict["overall_verdict"], "CONFIRMED")

        # Failure on H-A1 (only 3/5 > 0.13) -> 未确认
        rows_fail = [dict(r) for r in rows_pass]
        rows_fail[2]["primary_slope_ci95_lo"] = 0.13  # boundary: must be strictly > 0.13
        verdict_fail = adjudicate_team_a(rows_fail)
        self.assertEqual(verdict_fail["overall_verdict"], "未确认 (NOT CONFIRMED)")


if __name__ == "__main__":
    unittest.main()
