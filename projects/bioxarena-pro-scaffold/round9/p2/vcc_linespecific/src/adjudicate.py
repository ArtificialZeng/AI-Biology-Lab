#!/usr/bin/env python3
"""
Strict pre-registered hypothesis adjudication & report generator (H-A1..H-A4).

Reads outputs produced by `run_eval.py` (`line_metrics_<tag>.csv`, `bootstrap_ci_<tag>.csv`,
`summary_<tag>.json`) and checks every hypothesis objectively. If any condition fails,
explicitly writes `未确认` (Rule 8: unfavorable results reported faithfully).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd


def adjudicate_results(
    df_lines: pd.DataFrame,
    df_boot: pd.DataFrame,
    summary: Dict[str, Any],
    primary_method: str = "Hybrid_ProjShared",
    min_lines_pass: int = 4,
    max_lines_fail_ha4: int = 1,
    slope_threshold: float = 0.13,
) -> Dict[str, Any]:
    """Evaluate H-A1 through H-A4 strictly on held-out line outputs."""
    lines = sorted(df_lines["target_line"].unique().tolist())
    n_lines = len(lines)

    sub_prim = df_lines[df_lines["method"] == primary_method].set_index("target_line")
    sub_cipher = df_lines[df_lines["method"] == "CIPHER_Linear"].set_index("target_line")

    # H-A1: primary or CIPHER line-specific slope 95% CI lower bound > 0.13 on >= min_lines_pass lines
    ha1_per_line = {}
    ha1_pass_prim = 0
    ha1_pass_either = 0
    for ln in lines:
        lo_p = float(sub_prim.loc[ln, "ls_slope_ols_top50_ci95_lo"])
        est_p = float(sub_prim.loc[ln, "ls_slope_ols_top50"])
        lo_c = float(sub_cipher.loc[ln, "ls_slope_ols_top50_ci95_lo"])
        est_c = float(sub_cipher.loc[ln, "ls_slope_ols_top50"])
        ok_p = bool(np.isfinite(lo_p) and lo_p > slope_threshold)
        ok_c = bool(np.isfinite(lo_c) and lo_c > slope_threshold)
        ha1_pass_prim += int(ok_p)
        ha1_pass_either += int(ok_p or ok_c)
        ha1_per_line[ln] = {
            f"{primary_method}_slope": est_p,
            f"{primary_method}_ci95_lo": lo_p,
            "CIPHER_Linear_slope": est_c,
            "CIPHER_Linear_ci95_lo": lo_c,
            "pass_primary": ok_p,
            "pass_either": bool(ok_p or ok_c),
        }

    # H-A2: noise-calibrated line-specific correlation 95% CI > 0 on >= min_lines_pass lines
    ha2_per_line = {}
    ha2_pass = 0
    ha2_pass_resamp = 0
    for ln in lines:
        cal_r = float(sub_prim.loc[ln, "ls_resid_corr_cal_top50"])
        lo_fixed = float(sub_prim.loc[ln, "ls_resid_corr_cal_top50_ci95_lo"])
        lo_resamp = float(sub_prim.loc[ln, "ls_resid_corr_cal_top50_resamp_ci95_lo"])
        ok_f = bool(np.isfinite(lo_fixed) and lo_fixed > 0.0)
        ok_r = bool(np.isfinite(lo_resamp) and lo_resamp > 0.0)
        ha2_pass += int(ok_f)
        ha2_pass_resamp += int(ok_r)
        ha2_per_line[ln] = {
            "cal_r_top50": cal_r,
            "ci95_lo_fixed_ceil": lo_fixed,
            "ci95_lo_resamp_ceil": lo_resamp,
            "pass_fixed_ceil": ok_f,
            "pass_resamp_ceil": ok_r,
        }

    # H-A3: specificity - mismatched cov control slope significantly below matched cov (disjoint CIs)
    mis_map = {
        "Hybrid_ProjShared": "NegCtrl_HybridProj_MismatchedCov",
        "Hybrid_Raw": "NegCtrl_Hybrid_MismatchedCov",
        "CIPHER_Linear": "NegCtrl_CIPHER_MismatchedCov",
    }
    shuf_map = {
        "Hybrid_ProjShared": "NegCtrl_HybridProj_ShuffledGenes",
        "Hybrid_Raw": "NegCtrl_Hybrid_ShuffledGenes",
        "CIPHER_Linear": "NegCtrl_CIPHER_ShuffledGenes",
    }
    mis_name = mis_map.get(primary_method, "NegCtrl_HybridProj_MismatchedCov")
    shuf_name = shuf_map.get(primary_method, "NegCtrl_HybridProj_ShuffledGenes")
    b_ha3 = df_boot[
        (df_boot["method_a"] == primary_method)
        & (df_boot["method_b"] == mis_name)
        & (df_boot["metric"] == "ls_slope_ols_top50")
    ].set_index("target_line")
    ha3_per_line = {}
    ha3_pass = 0
    for ln in lines:
        if ln in b_ha3.index:
            row = b_ha3.loc[ln]
            ok = bool(row["disjoint_ci95"] and float(row["delta"]) > 0.0)
            ha3_per_line[ln] = {
                "matched_slope": float(row["mean_a"]),
                "mismatched_slope": float(row["mean_b"]),
                "disjoint_ci95": bool(row["disjoint_ci95"]),
                "pass": ok,
            }
            ha3_pass += int(ok)

    b_ha3_shuf = df_boot[
        (df_boot["method_a"] == primary_method)
        & (df_boot["method_b"] == shuf_name)
        & (df_boot["metric"] == "ls_slope_ols_top50")
    ].set_index("target_line")
    ha3_shuf_pass = 0
    for ln in lines:
        if ln in b_ha3_shuf.index:
            r_s = b_ha3_shuf.loc[ln]
            ha3_shuf_pass += int(bool(r_s["disjoint_ci95"] and float(r_s["delta"]) > 0.0))

    # H-A4: non-inferiority on cent_disc_raw, ctrl_disc_raw, cent_de50_r, ctrl_de50_r vs Transfer_Mean
    ha4_metrics = ["cent_disc_raw", "cent_de50_r", "ctrl_disc_raw", "ctrl_de50_r"]
    ha4_details = {}
    ha4_verdicts = {}
    for met in ha4_metrics:
        b_sub = df_boot[
            (df_boot["method_a"] == primary_method)
            & (df_boot["method_b"] == "Transfer_Mean")
            & (df_boot["metric"] == met)
        ].set_index("target_line")
        n_inferior = int(sum(bool(b_sub.loc[ln, "ci_lt_zero"]) for ln in lines if ln in b_sub.index))
        ha4_details[met] = {
            ln: {
                "delta": float(b_sub.loc[ln, "delta"]),
                "ci95_lo": float(b_sub.loc[ln, "ci95_lo"]),
                "ci95_hi": float(b_sub.loc[ln, "ci95_hi"]),
                "inferior_ci_lt_zero": bool(b_sub.loc[ln, "ci_lt_zero"]),
            }
            for ln in lines
            if ln in b_sub.index
        }
        ha4_verdicts[met] = {
            "n_inferior": n_inferior,
            "verdict": "CONFIRMED" if n_inferior <= max_lines_fail_ha4 else "未确认",
        }

    ha4_overall = all(v["verdict"] == "CONFIRMED" for v in ha4_verdicts.values())
    return {
        "n_lines": n_lines,
        "primary_method": primary_method,
        "H-A1": {
            "pass_count_primary": ha1_pass_prim,
            "pass_count_either": ha1_pass_either,
            "required": min_lines_pass,
            "verdict": "CONFIRMED" if ha1_pass_prim >= min_lines_pass else "未确认",
            "per_line": ha1_per_line,
        },
        "H-A2": {
            "pass_count_fixed_ceil": ha2_pass,
            "pass_count_resamp_ceil": ha2_pass_resamp,
            "required": min_lines_pass,
            "verdict": "CONFIRMED" if ha2_pass >= min_lines_pass else "未确认",
            "per_line": ha2_per_line,
        },
        "H-A3": {
            "pass_count": ha3_pass,
            "pass_count_shuffled_genes": ha3_shuf_pass,
            "required": min_lines_pass,
            "verdict": "CONFIRMED" if ha3_pass >= min_lines_pass else "未确认",
            "per_line": ha3_per_line,
        },
        "H-A4": {
            "max_inferior_allowed": max_lines_fail_ha4,
            "verdict": "CONFIRMED" if ha4_overall else "未确认",
            "by_metric": ha4_verdicts,
            "details": ha4_details,
        },
    }


def render_markdown_report(adj: Dict[str, Any], df_lines: pd.DataFrame, summary: Dict[str, Any], tag: str) -> str:
    lines_md: List[str] = [
        f"# Line-Specific Zero-Shot Perturbation Evaluation Report (`{tag}`)",
        "",
        f"- **Evaluated Lines ({adj['n_lines']})**: {', '.join(summary['lines'])}",
        f"- **Primary Method**: `{adj['primary_method']}`",
        "",
        "## Hypothesis Adjudication Summary",
        "",
        "| Hypothesis | Pass / Total | Threshold | Verdict |",
        "|---|---:|---|---|",
        f"| **H-A1** (Line-specific slope 95% CI > 0.13) | {adj['H-A1']['pass_count_primary']}/{adj['n_lines']} | >= {adj['H-A1']['required']} | **{adj['H-A1']['verdict']}** |",
        f"| **H-A2** (Noise-calibrated resid corr CI > 0) | {adj['H-A2']['pass_count_fixed_ceil']}/{adj['n_lines']} | >= {adj['H-A2']['required']} | **{adj['H-A2']['verdict']}** |",
        f"| **H-A3** (Matched > Mismatched cov, disjoint CI) | {adj['H-A3']['pass_count']}/{adj['n_lines']} | >= {adj['H-A3']['required']} | **{adj['H-A3']['verdict']}** |",
        f"| **H-A4** (Non-inferiority on disc_raw & de50_r) | {adj['H-A4']['by_metric']['cent_disc_raw']['n_inferior']}/{adj['n_lines']} inf | <= {adj['H-A4']['max_inferior_allowed']} inf | **{adj['H-A4']['verdict']}** |",
        "",
        "## Per-Line Summary by Method",
        "",
        "| Target Line | Method | N | Slope OLS (Top50) | 95% CI | Slope NoInt | Slope PerPert | Calib Resid r | Cent disc_raw | Cent de50_r | Ctrl disc_raw | Ctrl de50_r | MSE |",
        "|---|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for _, r in df_lines.iterrows():
        ci_str = f"[{r['ls_slope_ols_top50_ci95_lo']:+.4f}, {r['ls_slope_ols_top50_ci95_hi']:+.4f}]"
        lines_md.append(
            f"| {r['target_line']} | `{r['method']}` | {int(r['n_eval_targets'])} | "
            f"{r['ls_slope_ols_top50']:+.4f} | {ci_str} | {r['ls_slope_noint_top50']:+.4f} | "
            f"{r['ls_slope_per_pert_top50']:+.4f} | {r['ls_resid_corr_cal_top50']:+.4f} | "
            f"{r['cent_disc_raw']:.4f} | {r['cent_de50_r']:.4f} | "
            f"{r['ctrl_disc_raw']:.4f} | {r['ctrl_de50_r']:.4f} | {r['mse']:.4f} |"
        )
    return "\n".join(lines_md) + "\n"


def main() -> None:
    p = argparse.ArgumentParser(description="Adjudicate line-specific evaluation results.")
    p.add_argument("--results-dir", type=Path, default=Path(__file__).resolve().parents[1] / "results")
    p.add_argument("--tag", required=True, help="e.g. dev_log1p or holdout_log1p")
    p.add_argument("--primary-method", default="Hybrid_ProjShared")
    p.add_argument("--min-lines-pass", type=int, default=4)
    p.add_argument("--max-lines-fail-ha4", type=int, default=1)
    args = p.parse_args()

    df_lines = pd.read_csv(args.results_dir / f"line_metrics_{args.tag}.csv")
    df_boot = pd.read_csv(args.results_dir / f"bootstrap_ci_{args.tag}.csv")
    summary = json.loads((args.results_dir / f"summary_{args.tag}.json").read_text(encoding="utf-8"))

    adj = adjudicate_results(
        df_lines=df_lines,
        df_boot=df_boot,
        summary=summary,
        primary_method=args.primary_method,
        min_lines_pass=args.min_lines_pass,
        max_lines_fail_ha4=args.max_lines_fail_ha4,
    )
    (args.results_dir / f"adjudication_{args.tag}.json").write_text(
        json.dumps(adj, indent=2), encoding="utf-8"
    )
    md = render_markdown_report(adj, df_lines, summary, args.tag)
    (args.results_dir / f"REPORT_{args.tag}.md").write_text(md, encoding="utf-8")
    print(f"Wrote adjudication_{args.tag}.json and REPORT_{args.tag}.md", flush=True)


if __name__ == "__main__":
    main()
