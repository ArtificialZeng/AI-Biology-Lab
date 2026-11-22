#!/usr/bin/env python3
"""Independent A-T3 vs A-T2 R10-E cross-check script (post-hoc exploratory)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd

THIS_DIR = Path(__file__).resolve().parent
A2_EXP_DIR = THIS_DIR.parents[1] / "exploratory"

METHOD_MAP = {
    "FCT": "FCT",
    "E1": "E1_ConstMed",
    "E6": "E6_LogFCBar",
    "E7": "E7_FCBarMinus1",
}


def run_r10e_cross_check(tol_pt: float = 1e-6, tol_ci: float = 1e-6) -> Dict[str, Any]:
    a3_m = pd.read_csv(THIS_DIR / "audit_r10e_metrics.csv")
    a3_p = pd.read_csv(THIS_DIR / "audit_r10e_paired_fct_minus_e1.csv")
    a2_m = pd.read_csv(A2_EXP_DIR / "baseline_ablation_metrics.csv")
    a2_p = pd.read_csv(A2_EXP_DIR / "baseline_ablation_paired_fct_vs_e1.csv")

    a2_m_idx = {
        (r["target_line"], r["mask"], r["method"]): r
        for _, r in a2_m.iterrows()
    }
    a2_p_idx = {
        (r["target_line"], r["mask"]): r
        for _, r in a2_p.iterrows()
    }

    checks: List[Dict[str, Any]] = []
    max_pt_diff = 0.0
    max_ci_diff = 0.0

    metric_cols = [
        ("r_orth_cal", "r_orth_cal", False),
        ("r_orth_cal_ci95_lo", "r_orth_cal_ci95_lo", True),
        ("r_orth_cal_ci95_hi", "r_orth_cal_ci95_hi", True),
        ("r_orth_cal_ci975_lo", "r_orth_cal_ci975_lo", True),
        ("r_orth_cal_ci975_hi", "r_orth_cal_ci975_hi", True),
        ("r_orth_raw", "r_orth_raw", False),
        ("r_cal", "r_cal", False),
        ("slope_ols", "slope_ols", False),
        ("rho_ab_orth", "rho_half_ab_orth", False),
        ("noise_ceiling_r_orth", "ceiling_r_orth", False),
        ("e1_const_median", "mu_t_median", False),
    ]

    for _, r3 in a3_m.iterrows():
        key2 = (r3["target_line"], r3["mask"], METHOD_MAP[r3["predictor"]])
        if key2 not in a2_m_idx:
            raise KeyError(f"Missing A-T2 metric row for {key2}")
        r2 = a2_m_idx[key2]
        for c3, c2, is_ci in metric_cols:
            v3, v2 = float(r3[c3]), float(r2[c2])
            diff = abs(v3 - v2)
            if is_ci:
                max_ci_diff = max(max_ci_diff, diff)
            else:
                max_pt_diff = max(max_pt_diff, diff)
            checks.append(
                {
                    "table": "metrics",
                    "line": r3["target_line"],
                    "mask": r3["mask"],
                    "predictor": r3["predictor"],
                    "field": c3,
                    "a3": v3,
                    "a2": v2,
                    "abs_diff": diff,
                    "is_ci": is_ci,
                }
            )

    pair_cols = [
        ("fct_r_orth_cal", "fct_r_orth_cal", False),
        ("e1_r_orth_cal", "e1_r_orth_cal", False),
        ("diff_est", "diff_fct_minus_e1", False),
        ("diff_ci95_lo", "diff_ci95_lo", True),
        ("diff_ci95_hi", "diff_ci95_hi", True),
        ("diff_ci975_lo", "diff_ci975_lo", True),
        ("diff_ci975_hi", "diff_ci975_hi", True),
    ]
    bool_mismatches = []
    for _, r3 in a3_p.iterrows():
        key2 = (r3["target_line"], r3["mask"])
        if key2 not in a2_p_idx:
            raise KeyError(f"Missing A-T2 paired row for {key2}")
        r2 = a2_p_idx[key2]
        if bool(r3["diff_ci975_excludes_zero"]) != bool(r2["ci975_excludes_zero"]):
            bool_mismatches.append({"key": key2, "a3": bool(r3["diff_ci975_excludes_zero"]), "a2": bool(r2["ci975_excludes_zero"])})
        for c3, c2, is_ci in pair_cols:
            v3, v2 = float(r3[c3]), float(r2[c2])
            diff = abs(v3 - v2)
            if is_ci:
                max_ci_diff = max(max_ci_diff, diff)
            else:
                max_pt_diff = max(max_pt_diff, diff)
            checks.append(
                {
                    "table": "paired_fct_minus_e1",
                    "line": r3["target_line"],
                    "mask": r3["mask"],
                    "predictor": "FCT_minus_E1",
                    "field": c3,
                    "a3": v3,
                    "a2": v2,
                    "abs_diff": diff,
                    "is_ci": is_ci,
                }
            )

    checks_sorted = sorted(checks, key=lambda d: d["abs_diff"], reverse=True)
    result = {
        "analysis_type": "post_hoc_exploratory",
        "n_metric_rows_compared": len(a3_m),
        "n_paired_rows_compared": len(a3_p),
        "n_total_scalar_checks": len(checks),
        "n_bool_mismatches": len(bool_mismatches),
        "bool_mismatches": bool_mismatches,
        "max_point_abs_diff": max_pt_diff,
        "max_ci_abs_diff": max_ci_diff,
        "passes_point_tolerance_1e_6": bool(max_pt_diff < tol_pt and len(bool_mismatches) == 0),
        "passes_ci_tolerance_1e_6": bool(max_ci_diff < tol_ci),
        "top_10_diffs": checks_sorted[:10],
    }
    (THIS_DIR / "cross_check_r10e.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


if __name__ == "__main__":
    out = run_r10e_cross_check()
    print(json.dumps(
        {k: v for k, v in out.items() if k != "top_10_diffs"} | {"top_3_diffs": out["top_10_diffs"][:3]},
        indent=2,
    ))
