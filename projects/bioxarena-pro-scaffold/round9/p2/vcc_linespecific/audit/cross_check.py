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
                    "ctrl_disc_ownmasked_hybrid": float(r["ctrl_disc_raw"]),
                    "de50_r_hybrid": float(r["ctrl_de50_r"]),
                    "cent_disc_ownmasked_hybrid": float(r["cent_disc_raw"]),
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
            col_aliases = (
                ("ls_resid_corr_top50", "resid_corr"),
                ("ls_resid_corr_top50_ci95_lo", "resid_corr_ci95_lo"),
                ("ls_resid_corr_top50_ci95_hi", "resid_corr_ci95_hi"),
                ("ls_resid_corr_orth_top50", "r_orth_resid"),
                ("ls_resid_corr_orth_cal_top50", "r_orth_cal"),
                ("ls_resid_corr_orth_cal_top50_sq", "r_orth_cal_sq"),
                ("ls_resid_corr_orth_cal_all", "hybrid_r_orth_cal_all_nontarget"),
                ("ls_resid_corr_cal_all", "hybrid_cal_r_all_nontarget"),
                ("ls_resid_corr_cal_top50_sq", "r_cal_sq"),
                ("ls_r_cal_sq_top50", "r_cal_sq"),
                ("ls_resid_corr_cal_top50_sq_resamp_ci95_lo", "r_cal_sq_ci95_lo"),
                ("ls_resid_corr_cal_top50_sq_resamp_ci95_hi", "r_cal_sq_ci95_hi"),
                ("ls_sd_ratio_top50", "sd_pred_over_sd_obs"),
                ("sd_pred_over_sd_obs", "sd_pred_over_sd_obs"),
                ("ls_slope_obs_on_pred_top50", "slope_obs_on_pred"),
                ("slope_obs_on_pred", "slope_obs_on_pred"),
                ("ls_slope_deming_top50", "slope_deming"),
                ("slope_deming", "slope_deming"),
            )
            for col_t2, col_au in col_aliases:
                if col_au not in rec and col_t2 in r.index and np.isfinite(r[col_t2]):
                    rec[col_au] = float(r[col_t2])
        if "Transfer_Mean" in by_m.index:
            rt = by_m.loc["Transfer_Mean"]
            rec["ctrl_disc_ownmasked_tmean"] = float(rt["ctrl_disc_raw"])
            rec["de50_r_tmean"] = float(rt["ctrl_de50_r"])
            rec["cent_disc_ownmasked_tmean"] = float(rt["cent_disc_raw"])
            rec["cent_de50_r_tmean"] = float(rt["cent_de50_r"])
        if "Transfer_Mean_GlobalScale" in by_m.index:
            rg = by_m.loc["Transfer_Mean_GlobalScale"]
            for ck_src, ck_dst in (
                ("ls_slope_ols_top50", "global_scale_slope"),
                ("ls_resid_corr_top50", "global_scale_resid_corr"),
                ("ls_resid_corr_cal_top50", "global_scale_cal_r"),
                ("ls_resid_corr_orth_cal_top50", "global_scale_r_orth_cal"),
                ("ls_resid_corr_orth_cal_all", "global_scale_r_orth_cal_all_nontarget"),
            ):
                if ck_src in rg.index and np.isfinite(rg[ck_src]):
                    rec[ck_dst] = float(rg[ck_src])
        if "Transfer_Mean_VarScaled" in by_m.index:
            rv = by_m.loc["Transfer_Mean_VarScaled"]
            rec["tmvs_slope"] = float(rv["ls_slope_ols_top50"])
            rec["tmvs_cal_r"] = float(rv["ls_resid_corr_cal_top50"])
            for ck_src, ck_dst in (
                ("ls_slope_ols_top50_ci95_lo", "tmvs_slope_ci95_lo"),
                ("ls_slope_ols_top50_ci95_hi", "tmvs_slope_ci95_hi"),
                ("ls_slope_ols_top50_ci975_lo", "tmvs_slope_ci975_lo"),
                ("ls_slope_ols_top50_ci975_hi", "tmvs_slope_ci975_hi"),
                ("ls_resid_corr_cal_top50_resamp_ci95_lo", "tmvs_cal_r_ci95_lo"),
                ("ls_resid_corr_cal_top50_resamp_ci95_hi", "tmvs_cal_r_ci95_hi"),
                ("ls_resid_corr_cal_top50_resamp_ci975_lo", "tmvs_cal_r_ci975_lo"),
                ("ls_resid_corr_cal_top50_resamp_ci975_hi", "tmvs_cal_r_ci975_hi"),
                ("ls_resid_corr_orth_top50", "tmvs_r_orth_resid"),
                ("ls_resid_corr_orth_cal_top50", "tmvs_r_orth_cal"),
                ("ls_resid_corr_orth_cal_top50_sq", "tmvs_r_orth_cal_sq"),
                ("ls_resid_corr_cal_top50_sq", "tmvs_r_cal_sq"),
                ("ls_sd_ratio_top50", "tmvs_sd_pred_over_sd_obs"),
                ("ls_slope_deming_top50", "tmvs_slope_deming"),
                ("ls_slope_obs_on_pred_top50", "tmvs_slope_obs_on_pred"),
                ("ls_resid_corr_cal_all", "tmvs_cal_r_all_nontarget"),
                ("ls_resid_corr_orth_cal_all", "tmvs_r_orth_cal_all_nontarget"),
                ("ls_slope_ols_all", "tmvs_slope_all_nontarget"),
            ):
                if ck_src in rv.index and np.isfinite(rv[ck_src]):
                    rec[ck_dst] = float(rv[ck_src])
        if "Hybrid_on_TMVS" in by_m.index:
            rh = by_m.loc["Hybrid_on_TMVS"]
            rec["hotmvs_slope"] = float(rh["ls_slope_ols_top50"])
            rec["hotmvs_cal_r"] = float(rh["ls_resid_corr_cal_top50"])
            if "ls_resid_corr_orth_cal_top50" in rh.index and np.isfinite(rh["ls_resid_corr_orth_cal_top50"]):
                rec["hotmvs_r_orth_cal"] = float(rh["ls_resid_corr_orth_cal_top50"])
        if "FCT" in by_m.index:
            rf = by_m.loc["FCT"]
            rec["fct_slope"] = float(rf["ls_slope_ols_top50"])
            rec["fct_cal_r"] = float(rf["ls_resid_corr_cal_top50"])
            if "ls_resid_corr_orth_cal_top50" in rf.index and np.isfinite(rf["ls_resid_corr_orth_cal_top50"]):
                rec["fct_r_orth_cal"] = float(rf["ls_resid_corr_orth_cal_top50"])
        if "NegCtrl_TMVS_ShuffledVar" in by_m.index:
            rsv = by_m.loc["NegCtrl_TMVS_ShuffledVar"]
            rec["tmvs_shuf_var_slope"] = float(rsv["ls_slope_ols_top50"])
            rec["tmvs_shuf_var_cal_r"] = float(rsv["ls_resid_corr_cal_top50"])
        m_str = grp["method"].astype(str)
        tmvs_ctrl_rows = grp[
            m_str.str.startswith("NegCtrl_TMVS_MismatchedVar(") | (m_str == "NegCtrl_TMVS_ShuffledVar")
        ]
        if len(tmvs_ctrl_rows) > 0:
            for col_c, dst_c in (
                ("ls_slope_ols_top50_ci95_hi", "tmvs_max_ctrl_slope_ci95_hi"),
                ("ls_slope_ols_top50_ci975_hi", "tmvs_max_ctrl_slope_ci975_hi"),
                ("ls_resid_corr_cal_top50_resamp_ci95_hi", "tmvs_max_ctrl_cal_r_ci95_hi"),
                ("ls_resid_corr_cal_top50_resamp_ci975_hi", "tmvs_max_ctrl_cal_r_ci975_hi"),
            ):
                if col_c in tmvs_ctrl_rows.columns:
                    vals = tmvs_ctrl_rows[col_c].astype(float)
                    vals = vals[np.isfinite(vals)]
                    if len(vals) > 0:
                        rec[dst_c] = float(vals.max())
        if "NegCtrl_HybridProj_SrcCorrTgtVar" in by_m.index:
            rs_tv = by_m.loc["NegCtrl_HybridProj_SrcCorrTgtVar"]
            rec["sctv_slope"] = float(rs_tv["ls_slope_ols_top50"])
            rec["sctv_slope_ci95_hi"] = float(rs_tv["ls_slope_ols_top50_ci95_hi"])
            rec["sctv_cal_r"] = float(rs_tv["ls_resid_corr_cal_top50"])
            if "ls_resid_corr_orth_cal_top50" in rs_tv.index and np.isfinite(rs_tv["ls_resid_corr_orth_cal_top50"]):
                rec["sctv_r_orth_cal"] = float(rs_tv["ls_resid_corr_orth_cal_top50"])
            if "ls_resid_corr_cal_top50_resamp_ci95_hi" in rs_tv.index and np.isfinite(
                rs_tv["ls_resid_corr_cal_top50_resamp_ci95_hi"]
            ):
                rec["sctv_cal_r_ci95_hi"] = float(rs_tv["ls_resid_corr_cal_top50_resamp_ci95_hi"])
        if "CIPHER_Linear" in by_m.index:
            rec["cipher_linear_slope"] = float(by_m.loc["CIPHER_Linear", "ls_slope_ols_top50"])
        if mismatch_name in by_m.index:
            rm = by_m.loc[mismatch_name]
            rec["mismatch_slope"] = float(rm["ls_slope_ols_top50"])
            rec["mismatch_slope_ci95_lo"] = float(rm["ls_slope_ols_top50_ci95_lo"])
            rec["mismatch_slope_ci95_hi"] = float(rm["ls_slope_ols_top50_ci95_hi"])
        # Support R3 multi-source mismatch rows (`NegCtrl_HybridProj_MismatchedCov(<line>)` or `::` or `_`)
        mis_rows = grp[
            m_str.str.startswith(f"{mismatch_name}(")
            | m_str.str.startswith(f"{mismatch_name}::")
            | m_str.str.startswith(f"{mismatch_name}_")
            | (m_str == mismatch_name)
        ]
        if len(mis_rows) > 0:
            idx_max = mis_rows["ls_slope_ols_top50_ci95_hi"].astype(float).idxmax()
            rec["max_mismatch_slope"] = float(mis_rows.loc[idx_max, "ls_slope_ols_top50"])
            rec["max_mismatch_ci95_hi"] = float(mis_rows.loc[idx_max, "ls_slope_ols_top50_ci95_hi"])
            if "ls_slope_ols_top50_ci975_hi" in mis_rows.columns:
                rec["max_mismatch_ci975_hi"] = float(mis_rows["ls_slope_ols_top50_ci975_hi"].astype(float).max())
            if "ls_resid_corr_orth_cal_top50_resamp_ci975_hi" in mis_rows.columns:
                ov = mis_rows["ls_resid_corr_orth_cal_top50_resamp_ci975_hi"].astype(float)
                ov = ov[np.isfinite(ov)]
                rec["max_mismatch_r_orth_cal_ci975_hi"] = float(ov.max()) if len(ov) > 0 else 0.0
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
    conf_lines = v.get("confirmatory_lines", audit_summary["lines"])
    gov = v.get("r4_wording_governance", {})
    gov_r14 = v.get("r14_mechanism_governance", {})
    vb = v.get("Group_B")
    ci_pct = f"{100.0 * v.get('ci_level', 0.975):.1f}%"
    lines = [
        "# Team A Independent Statistical Audit Report (A-T3)",
        "",
        f"- **Evaluated Cell Lines ({len(audit_summary['lines'])})**: {', '.join(audit_summary['lines'])}",
        f"- **Confirmatory Cell Lines ({len(conf_lines)})**: {', '.join(conf_lines)}",
        f"- **Common Aligned Genes**: {audit_summary['n_common_genes']}",
        f"- **Confirmatory Two-Sided CI Level (R13 Bonferroni)**: **{ci_pct}**",
        f"- **Independent Group A Verdict (`Hybrid_ProjShared`)**: **{v['overall_verdict']}**",
        f"- **Independent Group B Verdict (`Transfer_Mean_VarScaled`)**: **{vb['overall_verdict']}**" if vb else "- **Independent Group B Verdict**: NA",
        "",
        "## 1. Pre-Registered Hypothesis Adjudication",
        "",
        f"### Group A (`Hybrid_ProjShared`, {ci_pct} CI)",
        "",
        f"- **H-A1** (Line-specific OLS slope {ci_pct} CI lower bound > 0.13): {'PASS' if v['H-A1']['pass'] else 'FAIL'} ({v['H-A1']['count']}/{v['n_lines']} lines: {v['H-A1']['lines']})",
        f"- **H-A2** (De-amplituded `r_orth_cal` {ci_pct} CI lower bound > 0, R17): {'PASS' if v['H-A2']['pass'] else 'FAIL'} ({v['H-A2']['count']}/{v['n_lines']} lines: {v['H-A2']['lines']})",
        f"- **H-A3** (Matched > ALL K-1 mismatched covs on BOTH `r_orth_cal` AND slope, R18): {'PASS' if v['H-A3']['pass'] else 'FAIL'} ({v['H-A3']['count']}/{v['n_lines']} lines: {v['H-A3']['lines']})",
        f"- **H-A4** (Hybrid non-inferior to Transferred Mean across ALL 6 conventions, {ci_pct} CI): {'PASS' if v['H-A4']['pass'] else 'FAIL'} "
        f"(ctrl_disc_panel inf: {v['H-A4']['disc_raw_inferior_count']}/{v['n_lines']}, "
        f"ctrl_disc_ownmasked inf: {v['H-A4'].get('ctrl_disc_ownmasked_inferior_count', 'NA')}/{v['n_lines']}, "
        f"ctrl_de50_r inf: {v['H-A4']['de50_r_inferior_count']}/{v['n_lines']}, "
        f"cent_disc_panel inf: {v['H-A4'].get('cent_disc_raw_inferior_count', 'NA')}/{v['n_lines']}, "
        f"cent_disc_ownmasked inf: {v['H-A4'].get('cent_disc_ownmasked_inferior_count', 'NA')}/{v['n_lines']}, "
        f"cent_de50_r inf: {v['H-A4'].get('cent_de50_r_inferior_count', 'NA')}/{v['n_lines']})",
        "",
    ]
    if vb:
        hb4 = vb.get("H-B4", {})
        vcc_rng = vb.get("vcc2026_top10_range_check", {})
        lines.extend(
            [
                f"### Group B (`Transfer_Mean_VarScaled`, {ci_pct} CI, R12/R17/R18)",
                "",
                f"- **H-B1** (TMVS OLS slope {ci_pct} CI lower bound > 0.13): {'PASS' if vb['H-B1']['pass'] else 'FAIL'} ({vb['H-B1']['count']}/{vb['n_lines']} lines: {vb['H-B1']['lines']})",
                f"- **H-B2** (TMVS de-amplituded `r_orth_cal` {ci_pct} CI lower bound > 0, R17): {'PASS' if vb['H-B2']['pass'] else 'FAIL'} ({vb['H-B2']['count']}/{vb['n_lines']} lines: {vb['H-B2']['lines']})",
                f"- **H-B3** (TMVS `r_orth_cal` CI_lo > max(all K-1 mismatched var CI_hi, 100-perm 97.5th pct), R18): {'PASS' if vb['H-B3']['pass'] else 'FAIL'} ({vb['H-B3']['count']}/{vb['n_lines']} lines: {vb['H-B3']['lines']})",
                f"- **H-B4** (TMVS non-inferior across ALL 6 conventions, {ci_pct} CI): {'PASS' if hb4.get('pass') else 'FAIL'} "
                f"(ctrl_disc_panel={hb4.get('ctrl_disc_panel_inferior_count', 'NA')}, ctrl_disc_ownmasked={hb4.get('ctrl_disc_ownmasked_inferior_count', 'NA')}, "
                f"ctrl_de50_r={hb4.get('ctrl_de50_r_inferior_count', 'NA')}, cent_disc_panel={hb4.get('cent_disc_panel_inferior_count', 'NA')}, "
                f"cent_disc_ownmasked={hb4.get('cent_disc_ownmasked_inferior_count', 'NA')}, cent_de50_r={hb4.get('cent_de50_r_inferior_count', 'NA')})",
                f"- **Descriptive VCC 2026 Top-10 Range [0.02, 0.13]**: {vcc_rng.get('count', 0)}/{vb['n_lines']} confirmatory lines ({vcc_rng.get('lines', [])}) — *Caveat*: {vcc_rng.get('r6_caveat', '')}",
                "",
            ]
        )
    if gov or gov_r14:
        r14_key = (
            "hotmvs_beats_tmvs_r_orth_cal_all_conf_ci95_lo_gt_0"
            if "hotmvs_beats_tmvs_r_orth_cal_all_conf_ci95_lo_gt_0" in gov_r14
            else "hotmvs_beats_tmvs_cal_r_all_conf_ci95_lo_gt_0"
        )
        lines.extend(
            [
                "### R4 & R14 Mechanism Wording Governance (R17 Updated)",
                "",
                f"- **Beats `SrcCorrTgtVar` on both slope & cal_r (all confirmatory lines, disjoint CI)**: `{gov.get('beats_src_corr_tgt_var_all_conf_slope_and_cal')}`",
                f"- **Beats `Transfer_Mean_VarScaled` on `r_orth_cal` (all confirmatory lines, disjoint CI)**: `{gov.get('beats_transfer_mean_var_scaled_all_conf_cal')}`",
                f"- **Permitted R4 Wording**: `{gov.get('allowed_mechanism_phrase')}`",
                f"- **Variance Caveat**: `{gov.get('marginal_variance_caveat')}`",
                f"- **R14 `Hybrid_on_TMVS > TMVS` on `r_orth_cal` (all confirmatory lines, paired 95% CI lo > 0)**: `{gov_r14.get(r14_key)}`",
                f"- **Permitted R14 Statement**: `{gov_r14.get('allowed_cofluctuation_increment_phrase')}`",
                "",
            ]
        )
    if "global_scale_resid_corr" in audit_df.columns:
        gs_vals = [f"{r['line']}: raw `r_resid={r['global_scale_resid_corr']:+.4f}` (`r_orth_cal={r['global_scale_r_orth_cal']:+.4f}`)" for _, r in audit_df.iterrows()]
        lines.extend(
            [
                f"> **R18 Global-Amplitude Artifact Callout**: Using **a single scalar per cell line (`Transfer_Mean_GlobalScale`)** already achieves raw residual correlation **{'; '.join(gs_vals)}**, collapsing near zero once per-perturbation amplitude is removed (`r_orth_cal`).",
                "",
            ]
        )
    lines.extend(
        [
            "## 2. Per-Line Independent Recomputation Table",
            "",
            audit_df.to_markdown(index=False),
            "",
        ]
    )
    orth_cols = [
        c
        for c in (
            "line",
            "rho_ab_orth_top50",
            "noise_ceiling_r_orth",
            "r_orth_cal",
            "r_orth_cal_ci975_lo",
            "r_orth_cal_ci975_hi",
            "tmvs_r_orth_cal",
            "tmvs_r_orth_cal_ci975_lo",
            "tmvs_r_orth_cal_ci975_hi",
            "global_scale_resid_corr",
            "global_scale_r_orth_cal",
            "tmvs_perm100_r_orth_cal_p025",
            "tmvs_perm100_r_orth_cal_median",
            "tmvs_perm100_r_orth_cal_p975",
            "tmvs_max_ctrl_r_orth_cal_ci975_hi",
            "fct_r_orth_cal",
            "hotmvs_r_orth_cal",
            "r14_hotmvs_minus_tmvs_r_orth_cal",
            "r14_hotmvs_minus_tmvs_r_orth_cal_ci95_lo",
            "rho_ab_orth_all_nontarget",
            "tmvs_r_orth_cal_all_nontarget",
            "hybrid_r_orth_cal_all_nontarget",
            "s5_tmvs_r_orth_cal",
            "s5_hybrid_r_orth_cal",
        )
        if c in audit_df.columns
    ]
    if len(orth_cols) > 1:
        lines.extend(
            [
                "### R17/R18 Orthogonalized (`r_orth_cal`), GlobalScale, & 100-Perm Shuffled-Var Table",
                "",
                audit_df[orth_cols].to_markdown(index=False),
                "",
            ]
        )
    # R14/R15 clean breakdown table
    s5_s6_cols = [
        c
        for c in (
            "line",
            "rho_ab_top50",
            "noise_ceiling_r",
            "rho_ab_all_nontarget",
            "noise_ceiling_r_all_nontarget",
            "hybrid_cal_r_all_nontarget",
            "tmvs_cal_r_all_nontarget",
            "s5_n_expressed_genes",
            "s5_hybrid_cal_r",
            "s5_tmvs_cal_r",
            "s5_tmvs_slope",
            "fct_slope",
            "fct_cal_r",
            "hotmvs_slope",
            "hotmvs_cal_r",
            "r14_hotmvs_minus_tmvs_cal_r",
            "r14_hotmvs_minus_tmvs_cal_r_ci95_lo",
            "r14_hotmvs_minus_tmvs_cal_r_ci95_hi",
        )
        if c in audit_df.columns
    ]
    if len(s5_s6_cols) > 1:
        lines.extend(
            [
                "### R14 Mechanism & R15 Secondary Baselines / S5 Expressed / S6 Raw Ceilings",
                "",
                audit_df[s5_s6_cols].to_markdown(index=False),
                "",
            ]
        )

    detailed = audit_summary.get("detailed", {})
    r19_rows = []
    for ln in audit_summary["lines"]:
        r19d = detailed.get(ln, {}).get("r19_gain_decomposition", {})
        for mask_name in ("resp_top50", "all_nontarget"):
            if mask_name in r19d:
                d_m = r19d[mask_name]
                r19_rows.append(
                    {
                        "target_line": ln,
                        "mask": mask_name,
                        "beta_t": d_m.get("beta_t"),
                        "beta_t_ci95": f"[{d_m.get('beta_t_ci95_lo', float('nan')):.3f}, {d_m.get('beta_t_ci95_hi', float('nan')):.3f}]",
                        "frac_whole_line": d_m.get("frac_whole_line"),
                        "whole_ci95": f"[{d_m.get('frac_whole_line_ci95_lo', float('nan')):.3f}, {d_m.get('frac_whole_line_ci95_hi', float('nan')):.3f}]",
                        "frac_per_gene_inc": d_m.get("frac_per_gene_inc"),
                        "gene_ci95": f"[{d_m.get('frac_per_gene_inc_ci95_lo', float('nan')):.3f}, {d_m.get('frac_per_gene_inc_ci95_hi', float('nan')):.3f}]",
                        "frac_remaining": d_m.get("frac_remaining"),
                        "rem_ci95": f"[{d_m.get('frac_remaining_ci95_lo', float('nan')):.3f}, {d_m.get('frac_remaining_ci95_hi', float('nan')):.3f}]",
                    }
                )
    if r19_rows:
        lines.extend(
            [
                "### R19 Descriptive Multiplicative Gain Variance Decomposition",
                "",
                pd.DataFrame(r19_rows).to_markdown(index=False),
                "",
            ]
        )

    ctrl_rows = []
    for ln in audit_summary["lines"]:
        dl = detailed.get(ln, {})
        for src, vals in dl.get("mismatch_by_source", {}).items():
            ctrl_rows.append(
                {
                    "target_line": ln,
                    "group": "Group_A_MismatchedCov",
                    "control": src,
                    "r_orth_cal": vals.get("r_orth_cal"),
                    "r_orth_cal_ci975_hi": vals.get("r_orth_cal_ci975_hi"),
                    "slope_ols": vals.get("slope_ols"),
                    "slope_ci975_hi": vals.get("slope_ols_ci975_hi"),
                    "noise_cal_corr": vals.get("noise_cal_corr"),
                }
            )
        for ck, vals in dl.get("tmvs_ctrls", {}).items():
            ctrl_rows.append(
                {
                    "target_line": ln,
                    "group": "Group_B_TMVS_Ctrl",
                    "control": ck,
                    "r_orth_cal": vals.get("r_orth_cal"),
                    "r_orth_cal_ci975_hi": vals.get("r_orth_cal_ci975_hi"),
                    "slope_ols": vals.get("slope_ols"),
                    "slope_ci975_hi": vals.get("slope_ols_ci975_hi"),
                    "noise_cal_corr": vals.get("noise_cal_corr"),
                }
            )
    if ctrl_rows:
        lines.extend(
            [
                "### R3 & R12 All Mismatched / Shuffled Control Specificity Details",
                "",
                pd.DataFrame(ctrl_rows).to_markdown(index=False),
                "",
            ]
        )

    if comparison is not None:
        our_dual = (
            f"Group A: {v['overall_verdict']} | Group B: {vb['overall_verdict']}"
            if vb
            else v["overall_verdict"]
        )
        lines.extend(
            [
                "## 3. Cross-Check Against A-T2 Results",
                "",
                f"- **Matched Lines**: {comparison['n_lines_matched']}",
                f"- **Numeric Comparisons**: {comparison['n_comparisons']}",
                f"- **Flagged Discrepancies (> {comparison['tol']})**: {comparison['n_flagged']}",
                f"- **Max Absolute Difference**: {comparison['max_abs_diff']:.6g}",
                f"- **T2 Verdict vs A-T3 Verdict**: `{t2_verdict}` vs `{our_dual}`",
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
        if "R13_Dual_Summary" in adj:
            t2_verdict = adj["R13_Dual_Summary"]
        else:
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
