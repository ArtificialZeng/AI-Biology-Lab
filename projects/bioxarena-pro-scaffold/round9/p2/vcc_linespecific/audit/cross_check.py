#!/usr/bin/env python3
"""Independent cross-check & discrepancy report generator for Team A (A-T3).

Compares A-T3's independently recomputed audit metrics (`audit_per_line.csv` /
`audit_summary.json`) against A-T2's reported outputs in `p2/vcc_linespecific/results/`,
and writes a concise, traceable Markdown audit report (`AUDIT_REPORT.md`).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


DEFAULT_TOL = 1e-4


def normalize_t2_long_table(
    t2_df: pd.DataFrame,
    primary_method: str = "Hybrid_ProjShared",
) -> pd.DataFrame:
    """Pivot A-T2's long `line_metrics_<tag>.csv` into wide format aligned with `audit_per_line.csv`."""
    if "target_line" not in t2_df.columns or "method" not in t2_df.columns:
        return t2_df

    mismatch_name = (
        "NegCtrl_HybridProj_MismatchedCov"
        if primary_method == "Hybrid_ProjShared"
        else "NegCtrl_Hybrid_MismatchedCov"
    )
    shuf_name = (
        "NegCtrl_HybridProj_ShuffledGenes"
        if primary_method == "Hybrid_ProjShared"
        else "NegCtrl_Hybrid_ShuffledGenes"
    )

    rows = []
    for ln, grp in t2_df.groupby("target_line", sort=True):
        by_m = grp.set_index("method")
        rec: Dict[str, object] = {"line": str(ln)}
        if primary_method in by_m.index:
            r = by_m.loc[primary_method]
            rec.update(
                {
                    "n_eval_targets": int(r["n_eval_targets"]),
                    "hybrid_slope": float(r["ls_slope_ols_top50"]),
                    "primary_slope_ci95_lo": float(r["ls_slope_ols_top50_ci95_lo"]),
                    "hybrid_slope_ci95_hi": float(r["ls_slope_ols_top50_ci95_hi"]),
                    "noise_cal_corr": float(r["ls_resid_corr_cal_top50"]),
                    "disc_raw_hybrid": float(r["ctrl_disc_raw"]),
                    "de50_r_hybrid": float(r["ctrl_de50_r"]),
                    "cent_disc_raw_hybrid": float(r["cent_disc_raw"]),
                    "cent_de50_r_hybrid": float(r["cent_de50_r"]),
                }
            )
            if "pds_cosine_panel" in r.index and np.isfinite(r["pds_cosine_panel"]):
                rec["pds_cosine_panel_hybrid"] = float(r["pds_cosine_panel"])
            if "ls_resid_corr_cal_top50_resamp_ci95_lo" in r.index and np.isfinite(
                r["ls_resid_corr_cal_top50_resamp_ci95_lo"]
            ):
                rec["noise_cal_corr_ci95_lo"] = float(r["ls_resid_corr_cal_top50_resamp_ci95_lo"])
                rec["noise_cal_corr_ci95_hi"] = float(r["ls_resid_corr_cal_top50_resamp_ci95_hi"])
        if "Transfer_Mean" in by_m.index:
            rt = by_m.loc["Transfer_Mean"]
            rec["disc_raw_tmean"] = float(rt["ctrl_disc_raw"])
            rec["de50_r_tmean"] = float(rt["ctrl_de50_r"])
            rec["cent_disc_raw_tmean"] = float(rt["cent_disc_raw"])
            rec["cent_de50_r_tmean"] = float(rt["cent_de50_r"])
        if "CIPHER_Linear" in by_m.index:
            rec["cipher_linear_slope"] = float(by_m.loc["CIPHER_Linear", "ls_slope_ols_top50"])
        if mismatch_name in by_m.index:
            rm = by_m.loc[mismatch_name]
            rec["mismatch_slope"] = float(rm["ls_slope_ols_top50"])
            rec["mismatch_slope_ci95_lo"] = float(rm["ls_slope_ols_top50_ci95_lo"])
            rec["mismatch_slope_ci95_hi"] = float(rm["ls_slope_ols_top50_ci95_hi"])
        if shuf_name in by_m.index:
            rec["shuf_target_slope"] = float(by_m.loc[shuf_name, "ls_slope_ols_top50"])
        rows.append(rec)
    return pd.DataFrame(rows)


def compare_tables(
    audit_df: pd.DataFrame,
    t2_df: pd.DataFrame,
    key_col: str = "line",
    column_pairs: Optional[List[Tuple[str, str]]] = None,
    tol: float = DEFAULT_TOL,
) -> dict:
    """Compare aligned numeric columns between A-T3 audit table and A-T2 results table."""
    if column_pairs is None:
        shared_cols = [c for c in audit_df.columns if c in t2_df.columns and c != key_col]
        column_pairs = [(c, c) for c in shared_cols if pd.api.types.is_numeric_dtype(audit_df[c])]

    merged = pd.merge(audit_df, t2_df, on=key_col, suffixes=("_audit", "_t2"), how="inner")
    diffs = []
    max_abs_diff = 0.0
    n_flagged = 0

    for _, row in merged.iterrows():
        line = row[key_col]
        for col_a, col_b in column_pairs:
            ca = f"{col_a}_audit" if col_a == col_b and f"{col_a}_audit" in merged.columns else col_a
            cb = f"{col_b}_t2" if col_a == col_b and f"{col_b}_t2" in merged.columns else col_b
            if ca not in row or cb not in row:
                continue
            va = float(row[ca])
            vb = float(row[cb])
            if np.isnan(va) and np.isnan(vb):
                abs_d = 0.0
                flagged = False
            elif np.isnan(va) or np.isnan(vb):
                abs_d = float("inf")
                flagged = True
            else:
                abs_d = abs(va - vb)
                flagged = abs_d > tol
                if abs_d > max_abs_diff:
                    max_abs_diff = abs_d
            if flagged:
                n_flagged += 1
            diffs.append(
                {
                    "line": str(line),
                    "metric_audit": col_a,
                    "metric_t2": col_b,
                    "val_audit": va,
                    "val_t2": vb,
                    "abs_diff": abs_d,
                    "within_tol": not flagged,
                }
            )

    return {
        "n_lines_matched": int(len(merged)),
        "n_comparisons": int(len(diffs)),
        "n_flagged": int(n_flagged),
        "max_abs_diff": float(max_abs_diff),
        "tol": float(tol),
        "all_match": (n_flagged == 0 and len(merged) == len(audit_df)),
        "rows": diffs,
    }


def compare_prediction_arrays(
    pred_a: np.ndarray,
    pred_b: np.ndarray,
    tol: float = 1e-5,
) -> dict:
    """Verify element-wise numerical agreement between two prediction matrices."""
    a = np.asarray(pred_a, dtype=np.float64)
    b = np.asarray(pred_b, dtype=np.float64)
    if a.shape != b.shape:
        return {"shape_match": False, "shape_a": list(a.shape), "shape_b": list(b.shape), "within_tol": False}
    diff = np.abs(a - b)
    max_d = float(np.nanmax(diff)) if diff.size > 0 else 0.0
    mean_d = float(np.nanmean(diff)) if diff.size > 0 else 0.0
    return {
        "shape_match": True,
        "max_abs_diff": max_d,
        "mean_abs_diff": mean_d,
        "within_tol": bool(max_d <= tol),
    }


def render_audit_markdown(
    audit_summary: dict,
    audit_df: pd.DataFrame,
    comparison: Optional[dict] = None,
    t2_verdict: Optional[str] = None,
) -> str:
    """Render an unvarnished Markdown audit report from independent audit artifacts."""
    v = audit_summary["verdict"]
    lines = [
        "# Team A Independent Statistical Audit Report (A-T3)",
        "",
        f"- **Evaluated Cell Lines ({len(audit_summary['lines'])})**: {', '.join(audit_summary['lines'])}",
        f"- **Common Aligned Genes**: {audit_summary['n_common_genes']}",
        f"- **Independent Overall Verdict**: **{v['overall_verdict']}**",
        "",
        "## 1. Pre-Registered Hypothesis Adjudication",
        "",
        f"- **H-A1** (Line-specific slope 95% CI lower bound > 0.13): {'PASS' if v['H-A1']['pass'] else 'FAIL'} ({v['H-A1']['count']}/{v['n_lines']} lines: {v['H-A1']['lines']})",
        f"- **H-A2** (Noise-calibrated correlation 95% CI lower bound > 0): {'PASS' if v['H-A2']['pass'] else 'FAIL'} ({v['H-A2']['count']}/{v['n_lines']} lines: {v['H-A2']['lines']})",
        f"- **H-A3** (Mismatched covariance slope CI strictly below matched CI): {'PASS' if v['H-A3']['pass'] else 'FAIL'} ({v['H-A3']['count']}/{v['n_lines']} lines: {v['H-A3']['lines']})",
        f"- **H-A4** (Hybrid non-inferior to Transferred Mean on disc_raw and de50_r): {'PASS' if v['H-A4']['pass'] else 'FAIL'} (disc_raw inferior: {v['H-A4']['disc_raw_inferior_count']}/{v['n_lines']}, de50_r inferior: {v['H-A4']['de50_r_inferior_count']}/{v['n_lines']})",
        "",
        "## 2. Per-Line Independent Recomputation Table",
        "",
        audit_df.to_markdown(index=False),
        "",
    ]

    if comparison is not None:
        lines.extend(
            [
                "## 3. Cross-Check Against A-T2 Results",
                "",
                f"- **Matched Lines**: {comparison['n_lines_matched']}",
                f"- **Numeric Comparisons**: {comparison['n_comparisons']}",
                f"- **Flagged Discrepancies (> {comparison['tol']})**: {comparison['n_flagged']}",
                f"- **Max Absolute Difference**: {comparison['max_abs_diff']:.6g}",
                f"- **T2 Verdict vs A-T3 Verdict**: `{t2_verdict}` vs `{v['overall_verdict']}`",
                "",
            ]
        )
        flagged = [r for r in comparison["rows"] if not r["within_tol"]]
        if flagged:
            lines.append("### Flagged Metric Discrepancies")
            lines.append("")
            lines.append(pd.DataFrame(flagged).to_markdown(index=False))
            lines.append("")

    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--audit-dir", type=Path, required=True)
    ap.add_argument("--t2-csv", type=Path, required=False)
    ap.add_argument("--t2-summary-json", type=Path, required=False)
    ap.add_argument("--t2-adj-json", type=Path, required=False)
    ap.add_argument("--t2-primary-method", default="Hybrid_ProjShared")
    ap.add_argument("--tol", type=float, default=DEFAULT_TOL)
    args = ap.parse_args()

    audit_df = pd.read_csv(args.audit_dir / "audit_per_line.csv")
    audit_summary = json.loads((args.audit_dir / "audit_summary.json").read_text())

    comp = None
    t2_verdict = None
    if args.t2_csv and args.t2_csv.exists():
        t2_raw = pd.read_csv(args.t2_csv)
        t2_df = normalize_t2_long_table(t2_raw, primary_method=args.t2_primary_method)
        comp = compare_tables(audit_df, t2_df, tol=args.tol)
        (args.audit_dir / "cross_check_diffs.json").write_text(json.dumps(comp, indent=2))
    if args.t2_adj_json and args.t2_adj_json.exists():
        adj = json.loads(args.t2_adj_json.read_text())
        hyp_verdicts = [adj.get(h, {}).get("verdict") for h in ("H-A1", "H-A2", "H-A3", "H-A4")]
        if all(v is not None for v in hyp_verdicts):
            t2_verdict = "CONFIRMED" if all(v == "CONFIRMED" for v in hyp_verdicts) else "未确认"
    elif args.t2_summary_json and args.t2_summary_json.exists():
        t2_sum = json.loads(args.t2_summary_json.read_text())
        t2_verdict = t2_sum.get("overall_verdict") or t2_sum.get("verdict", {}).get("overall_verdict")

    md = render_audit_markdown(audit_summary, audit_df, comparison=comp, t2_verdict=t2_verdict)
    (args.audit_dir / "AUDIT_REPORT.md").write_text(md)
    print(f"Wrote {args.audit_dir / 'AUDIT_REPORT.md'}")


if __name__ == "__main__":
    main()
