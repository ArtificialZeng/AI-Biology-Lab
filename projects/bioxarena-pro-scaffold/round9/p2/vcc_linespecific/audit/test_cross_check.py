#!/usr/bin/env python3
"""Unit tests for cross_check.py."""

from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from cross_check import (
    compare_prediction_arrays,
    compare_tables,
    normalize_t2_long_table,
    render_audit_markdown,
)


class TestCrossCheck(unittest.TestCase):
    """Verify table comparison, array comparison, and Markdown audit rendering."""

    def test_compare_tables_match_and_flag(self):
        df_a = pd.DataFrame(
            [
                {"line": "K562", "hybrid_slope": 0.2100, "disc_raw_hybrid": 0.6500},
                {"line": "RPE1", "hybrid_slope": 0.1400, "disc_raw_hybrid": 0.6100},
            ]
        )
        df_b = pd.DataFrame(
            [
                {"line": "K562", "hybrid_slope": 0.21002, "disc_raw_hybrid": 0.6500},
                {"line": "RPE1", "hybrid_slope": 0.14500, "disc_raw_hybrid": 0.6100},
            ]
        )
        res = compare_tables(df_a, df_b, tol=1e-4)
        self.assertFalse(res["all_match"])
        self.assertEqual(res["n_flagged"], 1)
        self.assertAlmostEqual(res["max_abs_diff"], 0.005, places=5)

    def test_normalize_t2_long_table(self):
        long_df = pd.DataFrame(
            [
                {
                    "target_line": "L1",
                    "method": "Hybrid_ProjShared",
                    "n_eval_targets": 20,
                    "ls_slope_ols_top50": 0.25,
                    "ls_slope_ols_top50_ci95_lo": 0.18,
                    "ls_slope_ols_top50_ci95_hi": 0.31,
                    "ls_resid_corr_top50": 0.22,
                    "ls_resid_corr_cal_top50": 0.35,
                    "ls_resid_corr_cal_top50_sq": 0.1225,
                    "ls_resid_corr_cal_top50_sq_resamp_ci95_lo": 0.04,
                    "ls_resid_corr_cal_top50_sq_resamp_ci95_hi": 0.25,
                    "ls_sd_ratio_top50": 0.11,
                    "ls_slope_obs_on_pred_top50": 1.85,
                    "ls_slope_deming_top50": 0.26,
                    "ls_resid_corr_cal_top50_resamp_ci95_lo": 0.20,
                    "ls_resid_corr_cal_top50_resamp_ci95_hi": 0.50,
                    "ctrl_disc_raw": 0.72,
                    "ctrl_de50_r": 0.45,
                    "cent_disc_raw": 0.73,
                    "cent_de50_r": 0.46,
                },
                {
                    "target_line": "L1",
                    "method": "Transfer_Mean",
                    "n_eval_targets": 20,
                    "ls_slope_ols_top50": 0.0,
                    "ls_slope_ols_top50_ci95_lo": 0.0,
                    "ls_slope_ols_top50_ci95_hi": 0.0,
                    "ls_resid_corr_cal_top50": 0.0,
                    "ctrl_disc_raw": 0.70,
                    "ctrl_de50_r": 0.44,
                    "cent_disc_raw": 0.71,
                    "cent_de50_r": 0.45,
                },
                {
                    "target_line": "L1",
                    "method": "NegCtrl_HybridProj_MismatchedCov::L2",
                    "n_eval_targets": 20,
                    "ls_slope_ols_top50": 0.01,
                    "ls_slope_ols_top50_ci95_lo": -0.04,
                    "ls_slope_ols_top50_ci95_hi": 0.05,
                    "ls_resid_corr_cal_top50": 0.01,
                    "ctrl_disc_raw": 0.70,
                    "ctrl_de50_r": 0.44,
                    "cent_disc_raw": 0.71,
                    "cent_de50_r": 0.45,
                },
                {
                    "target_line": "L1",
                    "method": "NegCtrl_HybridProj_MismatchedCov(L3)",
                    "n_eval_targets": 20,
                    "ls_slope_ols_top50": 0.02,
                    "ls_slope_ols_top50_ci95_lo": -0.02,
                    "ls_slope_ols_top50_ci95_hi": 0.07,
                    "ls_resid_corr_cal_top50": 0.02,
                    "ctrl_disc_raw": 0.70,
                    "ctrl_de50_r": 0.44,
                    "cent_disc_raw": 0.71,
                    "cent_de50_r": 0.45,
                },
                {
                    "target_line": "L1",
                    "method": "NegCtrl_HybridProj_SrcCorrTgtVar",
                    "n_eval_targets": 20,
                    "ls_slope_ols_top50": 0.03,
                    "ls_slope_ols_top50_ci95_lo": 0.00,
                    "ls_slope_ols_top50_ci95_hi": 0.06,
                    "ls_resid_corr_cal_top50": 0.04,
                    "ctrl_disc_raw": 0.70,
                    "ctrl_de50_r": 0.44,
                    "cent_disc_raw": 0.71,
                    "cent_de50_r": 0.45,
                },
                {
                    "target_line": "L1",
                    "method": "Transfer_Mean_VarScaled",
                    "n_eval_targets": 20,
                    "ls_slope_ols_top50": 0.01,
                    "ls_slope_ols_top50_ci95_lo": -0.01,
                    "ls_slope_ols_top50_ci95_hi": 0.03,
                    "ls_resid_corr_cal_top50": 0.02,
                    "ctrl_disc_raw": 0.70,
                    "ctrl_de50_r": 0.44,
                    "cent_disc_raw": 0.71,
                    "cent_de50_r": 0.45,
                },
            ]
        )
        wide = normalize_t2_long_table(long_df, primary_method="Hybrid_ProjShared")
        self.assertEqual(list(wide["line"]), ["L1"])
        self.assertAlmostEqual(wide.loc[0, "hybrid_slope"], 0.25)
        self.assertAlmostEqual(wide.loc[0, "r_cal_sq"], 0.1225)
        self.assertAlmostEqual(wide.loc[0, "sd_pred_over_sd_obs"], 0.11)
        self.assertAlmostEqual(wide.loc[0, "slope_deming"], 0.26)
        self.assertAlmostEqual(wide.loc[0, "max_mismatch_ci95_hi"], 0.07)
        self.assertAlmostEqual(wide.loc[0, "sctv_slope"], 0.03)
        self.assertAlmostEqual(wide.loc[0, "tmvs_cal_r"], 0.02)
        self.assertAlmostEqual(wide.loc[0, "ctrl_disc_ownmasked_tmean"], 0.70)

    def test_compare_prediction_arrays(self):
        a = np.ones((5, 10), dtype=np.float64)
        b = a + 1e-7
        res = compare_prediction_arrays(a, b, tol=1e-5)
        self.assertTrue(res["within_tol"])

    def test_render_audit_markdown(self):
        summary = {
            "lines": ["L1", "L2"],
            "n_common_genes": 100,
            "verdict": {
                "n_lines": 2,
                "H-A1": {"pass": False, "count": 1, "lines": ["L1"]},
                "H-A2": {"pass": True, "count": 2, "lines": ["L1", "L2"]},
                "H-A3": {"pass": False, "count": 0, "lines": []},
                "H-A4": {
                    "pass": True,
                    "disc_raw_inferior_count": 0,
                    "disc_raw_inferior_lines": [],
                    "de50_r_inferior_count": 0,
                    "de50_r_inferior_lines": [],
                },
                "overall_verdict": "未确认 (NOT CONFIRMED)",
            },
        }
        df = pd.DataFrame([{"line": "L1", "hybrid_slope": 0.15}, {"line": "L2", "hybrid_slope": 0.05}])
        md = render_audit_markdown(summary, df)
        self.assertIn("未确认 (NOT CONFIRMED)", md)
        self.assertIn("L1", md)


if __name__ == "__main__":
    unittest.main()
