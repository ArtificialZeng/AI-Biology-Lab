#!/usr/bin/env python3
"""
Strict pre-registered hypothesis adjudication & report generator (Group A: H-A1..A4; Group B: H-B1..B4).

Implements PREREG_A v4 / PREREG_AUDIT_R1..R3 rules:
  - R2/R10: Confirmatory cohort excludes Replogle dev lines and lines with < 50 targets.
  - R3/R12/R17/R18: Primary correlation is amplitude-orthogonalized r_orth_cal; H-A3 requires slope AND r_orth_cal > ALL mismatched covs; H-B3 requires r_orth_cal > max(mismatched vars, 100-perm shuffled var p97.5).
  - R4/R14/R20: Strict wording locks vs SrcCorrTgtVar, TMVS, Hybrid_on_TMVS on r_orth_cal, and PRIOR_ART_B novelty bounds.
  - R5: Strict AND across all 6 H-A4/H-B4 non-inferiority conventions.
  - R13: 97.5% two-sided Bonferroni CIs ([1.25, 98.75] percentiles) across the 2 confirmatory families.
  - R18/R19: Reports Transfer_Mean_GlobalScale scalar baseline and 3-layer residual variance decomposition.
  - Rule 8: Any failed criterion explicitly outputs `未确认`.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


def _is_dev_line(name: str) -> bool:
    low = name.lower()
    return ("replogle" in low) and ("k562" in low or "rpe1" in low)


def _fmt_ci(lo: float, hi: float, fmt: str = "+.4f") -> str:
    if not (np.isfinite(lo) and np.isfinite(hi)):
        return "[nan, nan]"
    return f"[{lo:{fmt}}, {hi:{fmt}}]"


def _resolve_ci_tag(df_lines: pd.DataFrame, df_boot: pd.DataFrame, ci_level: float) -> Tuple[str, float]:
    """Use ci975 when ci_level >= 0.97 and columns exist; fallback to ci95 only for legacy tests."""
    want_975 = ci_level >= 0.97
    has_975 = (
        "ls_slope_ols_top50_ci975_lo" in df_lines.columns
        and "ci975_lo" in df_boot.columns
        and "a_ci975_lo" in df_boot.columns
    )
    if want_975 and has_975:
        return "ci975", 0.975
    return "ci95", 0.95


def _adjudicate_slope_and_calr(
    sub_df: pd.DataFrame,
    all_lines: Sequence[str],
    eval_lines: Sequence[str],
    ci_tag: str,
    slope_threshold: float,
    min_lines_pass: int,
    ceilings: Optional[Dict[str, Any]] = None,
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Evaluate H-X1 (slope > 0.13) and H-X2 (joint-boot r_orth_cal > 0, valid >= 95%, rho_AB_orth > 0)."""
    use_orth = "ls_resid_corr_orth_cal_top50" in sub_df.columns
    cal_col = "ls_resid_corr_orth_cal_top50" if use_orth else "ls_resid_corr_cal_top50"
    sq_col = f"{cal_col}_sq"
    slope_lo_col = f"ls_slope_ols_top50_{ci_tag}_lo"
    resamp_lo_col = (
        f"{cal_col}_resamp_ci975_lo"
        if ci_tag == "ci975" and f"{cal_col}_resamp_ci975_lo" in sub_df.columns
        else f"{cal_col}_resamp_ci95_lo"
    )
    fixed_lo_col = f"{cal_col}_ci95_lo"
    vfrac_col = f"{cal_col}_resamp_valid_frac"
    has_vfrac = vfrac_col in sub_df.columns

    h1_per: Dict[str, Any] = {}
    h2_per: Dict[str, Any] = {}
    h1_pass = 0
    h2_pass_resamp = 0
    h2_pass_fixed = 0

    for ln in all_lines:
        is_conf = ln in eval_lines
        est_s = float(sub_df.loc[ln, "ls_slope_ols_top50"])
        lo_s = float(sub_df.loc[ln, slope_lo_col])
        ok_s = bool(np.isfinite(lo_s) and lo_s > slope_threshold)
        if is_conf:
            h1_pass += int(ok_s)
        h1_per[ln] = {
            "confirmatory": is_conf,
            "slope_ols_top50": est_s,
            f"slope_ols_top50_{ci_tag}_lo": lo_s,
            "pass": ok_s,
        }

        cal_r = float(sub_df.loc[ln, cal_col])
        cal_sq = (
            float(sub_df.loc[ln, sq_col])
            if sq_col in sub_df.columns
            else float(cal_r ** 2 if np.isfinite(cal_r) else np.nan)
        )
        raw_cal_r = float(sub_df.loc[ln, "ls_resid_corr_cal_top50"]) if "ls_resid_corr_cal_top50" in sub_df.columns else float("nan")
        lo_f = float(sub_df.loc[ln, fixed_lo_col])
        lo_r = float(sub_df.loc[ln, resamp_lo_col])
        vfrac = (
            float(sub_df.loc[ln, vfrac_col])
            if has_vfrac
            else (1.0 if np.isfinite(lo_r) else 0.0)
        )
        ceil_ln = (ceilings or {}).get(ln, {}).get("resp_top50", {})
        rho_ab_orth = float(ceil_ln.get("rho_half_ab_orth", float("nan"))) if use_orth else float(ceil_ln.get("rho_half_ab", float("nan")))
        rho_ok = (not ceilings or ln not in ceilings) or (np.isfinite(rho_ab_orth) and rho_ab_orth > 0.0)

        ok_f = bool(rho_ok and np.isfinite(cal_r) and cal_r > 0.0 and np.isfinite(lo_f) and lo_f > 0.0)
        ok_r = bool(
            rho_ok
            and np.isfinite(cal_r)
            and cal_r > 0.0
            and np.isfinite(vfrac)
            and vfrac >= 0.95
            and np.isfinite(lo_r)
            and lo_r > 0.0
        )
        if is_conf:
            h2_pass_fixed += int(ok_f)
            h2_pass_resamp += int(ok_r)
        h2_per[ln] = {
            "confirmatory": is_conf,
            "metric_used": cal_col,
            "rho_half_ab_orth": rho_ab_orth if use_orth else float("nan"),
            "cal_r_orth_top50": cal_r if use_orth else float("nan"),
            "cal_r_orth_sq_top50": cal_sq if use_orth else float("nan"),
            "cal_r_top50": raw_cal_r if use_orth else cal_r,
            "cal_r_sq_top50": float(raw_cal_r ** 2 if np.isfinite(raw_cal_r) else np.nan) if use_orth else cal_sq,
            "resamp_valid_frac": vfrac,
            "ci95_lo_fixed_ceil": lo_f,
            f"{ci_tag}_lo_resamp_ceil": lo_r,
            "pass_fixed_ceil": ok_f,
            "pass_resamp_ceil": ok_r,
        }

    h1_dict = {
        "pass_count": h1_pass,
        "pass_count_primary": h1_pass,
        "required": min_lines_pass,
        "verdict": "CONFIRMED" if h1_pass >= min_lines_pass else "未确认",
        "per_line": h1_per,
    }
    h2_dict = {
        "metric_used": cal_col,
        "pass_count_fixed_ceil": h2_pass_fixed,
        "pass_count_resamp_ceil": h2_pass_resamp,
        "required": min_lines_pass,
        "verdict": "CONFIRMED" if h2_pass_resamp >= min_lines_pass else "未确认",
        "per_line": h2_per,
    }
    return h1_dict, h2_dict


def _adjudicate_noninf_family(
    df_boot: pd.DataFrame,
    method_name: str,
    all_lines: Sequence[str],
    eval_lines: Sequence[str],
    ci_tag: str,
    max_lines_fail: int,
) -> Dict[str, Any]:
    """Evaluate H-X4 strict AND non-inferiority across all 6 disc_raw & de50_r conventions."""
    candidate_metrics = [
        "ctrl_disc_raw",
        "ctrl_de50_r",
        "ctrl_disc_raw_panel",
        "cent_disc_raw",
        "cent_de50_r",
        "cent_disc_raw_panel",
    ]
    avail_metrics = set(df_boot["metric"].unique().tolist())
    metrics = [m for m in candidate_metrics if m in avail_metrics]
    lt_col = "ci975_lt_zero" if (ci_tag == "ci975" and "ci975_lt_zero" in df_boot.columns) else "ci_lt_zero"
    lo_col = f"{ci_tag}_lo" if f"{ci_tag}_lo" in df_boot.columns else "ci95_lo"
    hi_col = f"{ci_tag}_hi" if f"{ci_tag}_hi" in df_boot.columns else "ci95_hi"

    details: Dict[str, Any] = {}
    verdicts: Dict[str, Any] = {}
    inf_union = set()
    inf_union_ci95 = set()

    for met in metrics:
        b_sub = df_boot[
            (df_boot["method_a"] == method_name)
            & (df_boot["method_b"] == "Transfer_Mean")
            & (df_boot["metric"] == met)
        ].set_index("target_line")
        inf_lines = [ln for ln in eval_lines if ln in b_sub.index and bool(b_sub.loc[ln, lt_col])]
        inf_lines_95 = [ln for ln in eval_lines if ln in b_sub.index and bool(b_sub.loc[ln, "ci_lt_zero"])]
        inf_union.update(inf_lines)
        inf_union_ci95.update(inf_lines_95)
        n_inf = len(inf_lines)
        details[met] = {
            ln: {
                "confirmatory": ln in eval_lines,
                "delta": float(b_sub.loc[ln, "delta"]),
                f"{ci_tag}_lo": float(b_sub.loc[ln, lo_col]),
                f"{ci_tag}_hi": float(b_sub.loc[ln, hi_col]),
                "inferior_ci_lt_zero": bool(b_sub.loc[ln, lt_col]),
                "inferior_ci95_lt_zero": bool(b_sub.loc[ln, "ci_lt_zero"]),
            }
            for ln in all_lines
            if ln in b_sub.index
        }
        verdicts[met] = {
            "n_inferior": n_inf,
            "n_inferior_ci95": len(inf_lines_95),
            "inferior_lines": inf_lines,
            "verdict": "CONFIRMED" if n_inf <= max_lines_fail else "未确认",
        }

    overall = bool(verdicts) and all(v["verdict"] == "CONFIRMED" for v in verdicts.values())
    max_inf = max((int(v["n_inferior"]) for v in verdicts.values()), default=0)
    return {
        "max_inferior_allowed": max_lines_fail,
        "max_n_inferior": max_inf,
        "n_inferior_union": len(inf_union),
        "n_inferior_union_ci95": len(inf_union_ci95),
        "inferior_union_lines": sorted(inf_union),
        "verdict": "CONFIRMED" if overall else "未确认",
        "by_metric": verdicts,
        "details": details,
    }


def adjudicate_results(
    df_lines: pd.DataFrame,
    df_boot: pd.DataFrame,
    summary: Dict[str, Any],
    primary_method: str = "Hybrid_ProjShared",
    group_b_method: str = "Transfer_Mean_VarScaled",
    min_lines_pass: Optional[int] = None,
    max_lines_fail_ha4: Optional[int] = None,
    slope_threshold: float = 0.13,
    confirmatory_lines: Optional[List[str]] = None,
    min_confirmatory_targets: int = 50,
    ci_level: float = 0.975,
) -> Dict[str, Any]:
    """Evaluate Group A (H-A1..A4), Group B (H-B1..B4), R3/R4/R5/R10/R12/R13/R14 rules strictly."""
    all_lines = sorted(df_lines["target_line"].unique().tolist())
    sub_prim = df_lines[df_lines["method"] == primary_method].set_index("target_line")
    sub_cipher = df_lines[df_lines["method"] == "CIPHER_Linear"].set_index("target_line")
    has_group_b = group_b_method in set(df_lines["method"].unique().tolist())
    sub_gb = (
        df_lines[df_lines["method"] == group_b_method].set_index("target_line")
        if has_group_b
        else None
    )

    if confirmatory_lines is not None:
        eval_lines = [ln for ln in all_lines if ln in set(confirmatory_lines)]
    elif summary.get("phase") == "holdout":
        eval_lines = [
            ln
            for ln in all_lines
            if (not _is_dev_line(ln))
            and int(sub_prim.loc[ln, "n_eval_targets"]) >= min_confirmatory_targets
        ]
    else:
        eval_lines = list(all_lines)
    if not eval_lines:
        eval_lines = list(all_lines)
    secondary_lines = [ln for ln in all_lines if ln not in set(eval_lines)]
    n_eval = len(eval_lines)

    if min_lines_pass is None:
        if len(secondary_lines) > 0 or len(all_lines) >= 5:
            min_lines_pass = int(math.ceil(0.8 * n_eval))
        else:
            min_lines_pass = max(4, int(math.ceil(0.8 * n_eval)))
    if max_lines_fail_ha4 is None:
        if len(secondary_lines) > 0 or len(all_lines) >= 5:
            max_lines_fail_ha4 = int(math.floor(0.2 * n_eval))
        else:
            max_lines_fail_ha4 = 1

    ci_tag, resolved_ci_level = _resolve_ci_tag(df_lines, df_boot, ci_level)
    a_lo_col = f"a_{ci_tag}_lo" if f"a_{ci_tag}_lo" in df_boot.columns else "a_ci95_lo"
    b_lo_col = f"b_{ci_tag}_lo" if f"b_{ci_tag}_lo" in df_boot.columns else "b_ci95_lo"
    b_hi_col = f"b_{ci_tag}_hi" if f"b_{ci_tag}_hi" in df_boot.columns else "b_ci95_hi"
    lo_col = f"{ci_tag}_lo" if f"{ci_tag}_lo" in df_boot.columns else "ci95_lo"
    gt_col = f"{ci_tag}_gt_zero" if f"{ci_tag}_gt_zero" in df_boot.columns else "ci_gt_zero"
    disj_col = f"disjoint_{ci_tag}" if f"disjoint_{ci_tag}" in df_boot.columns else "disjoint_ci95"

    ceilings = summary.get("pert_noise_ceilings", {})
    boot_metrics = set(df_boot["metric"].unique().tolist()) if not df_boot.empty else set()
    use_orth_lines = "ls_resid_corr_orth_cal_top50" in df_lines.columns
    use_orth_boot = use_orth_lines and ("ls_resid_corr_orth_cal_top50" in boot_metrics)
    orth_metric_boot: Optional[str] = (
        "ls_resid_corr_orth_cal_top50"
        if use_orth_boot
        else ("ls_resid_corr_cal_top50" if "ls_resid_corr_cal_top50" in boot_metrics else None)
    )
    prim_orth_col = "ls_resid_corr_orth_cal_top50" if use_orth_lines else "ls_resid_corr_cal_top50"
    prim_orth_lo_col = (
        f"{prim_orth_col}_resamp_{ci_tag}_lo"
        if f"{prim_orth_col}_resamp_{ci_tag}_lo" in sub_prim.columns
        else f"{prim_orth_col}_resamp_ci95_lo"
    )

    # ---- Group A: H-A1, H-A2, H-A3, H-A4 ----
    ha1, ha2 = _adjudicate_slope_and_calr(
        sub_prim, all_lines, eval_lines, ci_tag, slope_threshold, min_lines_pass, ceilings=ceilings
    )
    # Record CIPHER_Linear reference alongside H-A1 (cannot substitute for primary_method)
    c_lo_col = f"ls_slope_ols_top50_{ci_tag}_lo"
    for ln in all_lines:
        est_c = float(sub_cipher.loc[ln, "ls_slope_ols_top50"])
        lo_c = float(sub_cipher.loc[ln, c_lo_col])
        ha1["per_line"][ln]["CIPHER_Linear_slope"] = est_c
        ha1["per_line"][ln][f"CIPHER_Linear_{ci_tag}_lo"] = lo_c

    # H-A3 (R3/R17/R18): matched_lo > max_{s'} mismatch_hi(s') on BOTH slope AND r_orth_cal (AND)
    ha3_per_line: Dict[str, Any] = {}
    ha3_pass = 0
    for ln in all_lines:
        is_conf = ln in eval_lines
        if primary_method == "Hybrid_ProjShared":
            mis_mask = (
                df_boot["method_b"].str.startswith("NegCtrl_HybridProj_MismatchedCov(")
                | (df_boot["method_b"] == "NegCtrl_HybridProj_MismatchedCov")
            )
        else:
            legacy_mis = (
                "NegCtrl_Hybrid_MismatchedCov"
                if primary_method == "Hybrid_Raw"
                else "NegCtrl_CIPHER_MismatchedCov"
            )
            mis_mask = df_boot["method_b"] == legacy_mis

        b_mis_ln = df_boot[
            (df_boot["target_line"] == ln)
            & (df_boot["method_a"] == primary_method)
            & (df_boot["metric"] == "ls_slope_ols_top50")
            & mis_mask
        ]
        b_mis_orth_map: Dict[str, pd.Series] = {}
        require_ha3_orth = use_orth_boot
        if require_ha3_orth:
            b_mis_o = df_boot[
                (df_boot["target_line"] == ln)
                & (df_boot["method_a"] == primary_method)
                & (df_boot["metric"] == "ls_resid_corr_orth_cal_top50")
                & mis_mask
            ]
            for _, ro in b_mis_o.iterrows():
                b_mis_orth_map[str(ro["method_b"])] = ro

        by_src: Dict[str, Any] = {}
        matched_lo = float(sub_prim.loc[ln, f"ls_slope_ols_top50_{ci_tag}_lo"])
        matched_est = float(sub_prim.loc[ln, "ls_slope_ols_top50"])
        matched_orth_lo = float(sub_prim.loc[ln, prim_orth_lo_col]) if prim_orth_lo_col in sub_prim.columns else float("nan")
        matched_orth_est = float(sub_prim.loc[ln, prim_orth_col]) if prim_orth_col in sub_prim.columns else float("nan")

        max_mis_hi = float("-inf")
        max_mis_orth_hi = float("-inf")
        n_checked = 0
        all_src_pass = True
        for _, r_m in b_mis_ln.iterrows():
            src_label = str(r_m["method_b"])
            b_hi = float(r_m[b_hi_col])
            pass_slope = bool(np.isfinite(matched_lo) and np.isfinite(b_hi) and matched_lo > b_hi)
            if np.isfinite(b_hi):
                max_mis_hi = max(max_mis_hi, b_hi)

            pass_orth = True
            b_o_est = float("nan")
            b_o_lo = float("nan")
            b_o_hi = float("nan")
            if require_ha3_orth:
                if src_label in b_mis_orth_map:
                    ro = b_mis_orth_map[src_label]
                    b_o_est = float(ro["mean_b"])
                    b_o_lo = float(ro[b_lo_col])
                    b_o_hi = float(ro[b_hi_col])
                    zero_res_null = (
                        np.isfinite(float(r_m["mean_b"]))
                        and abs(float(r_m["mean_b"])) <= 1e-12
                        and abs(float(r_m[b_lo_col])) <= 1e-12
                        and np.isfinite(b_hi)
                        and abs(b_hi) <= 1e-12
                        and not np.isfinite(b_o_hi)
                    )
                    eff_o_hi = 0.0 if zero_res_null else b_o_hi
                    pass_orth = bool(np.isfinite(matched_orth_lo) and np.isfinite(eff_o_hi) and matched_orth_lo > eff_o_hi)
                    if np.isfinite(eff_o_hi):
                        max_mis_orth_hi = max(max_mis_orth_hi, eff_o_hi)
                    b_o_hi = eff_o_hi
                else:
                    pass_orth = False

            pass_both = bool(pass_slope and pass_orth)
            all_src_pass = all_src_pass and pass_both
            by_src[src_label] = {
                "mismatched_slope": float(r_m["mean_b"]),
                f"mismatched_{ci_tag}_lo": float(r_m[b_lo_col]),
                f"mismatched_{ci_tag}_hi": b_hi,
                "pass_vs_this_source_slope": pass_slope,
                "mismatched_r_orth_cal": b_o_est,
                f"mismatched_r_orth_cal_{ci_tag}_lo": b_o_lo,
                f"mismatched_r_orth_cal_{ci_tag}_hi": b_o_hi,
                "pass_vs_this_source_r_orth_cal": pass_orth,
                "pass_vs_this_source": pass_both,
            }
            n_checked += 1
        ok_ha3 = bool(n_checked >= 1 and all_src_pass)
        if is_conf:
            ha3_pass += int(ok_ha3)
        ha3_per_line[ln] = {
            "confirmatory": is_conf,
            "matched_slope": matched_est,
            f"matched_{ci_tag}_lo": matched_lo,
            "matched_r_orth_cal": matched_orth_est,
            f"matched_r_orth_cal_{ci_tag}_lo": matched_orth_lo,
            "n_sources_checked": n_checked,
            f"max_mismatch_{ci_tag}_hi": max_mis_hi if np.isfinite(max_mis_hi) else float("nan"),
            f"max_mismatch_r_orth_cal_{ci_tag}_hi": max_mis_orth_hi if np.isfinite(max_mis_orth_hi) else float("nan"),
            "pass": ok_ha3,
            "by_source": by_src,
        }

    shuf_map = {
        "Hybrid_ProjShared": "NegCtrl_HybridProj_ShuffledGenes",
        "Hybrid_Raw": "NegCtrl_Hybrid_ShuffledGenes",
        "CIPHER_Linear": "NegCtrl_CIPHER_ShuffledGenes",
    }
    shuf_name = shuf_map.get(primary_method, "NegCtrl_HybridProj_ShuffledGenes")
    shuf_metric = orth_metric_boot if orth_metric_boot is not None else "ls_slope_ols_top50"
    b_ha3_shuf = df_boot[
        (df_boot["method_a"] == primary_method)
        & (df_boot["method_b"] == shuf_name)
        & (df_boot["metric"] == shuf_metric)
    ].set_index("target_line")
    ha3_shuf_pass = sum(
        1
        for ln in eval_lines
        if ln in b_ha3_shuf.index
        and bool(
            b_ha3_shuf.loc[ln, disj_col]
            and float(b_ha3_shuf.loc[ln, "delta"]) > 0.0
            and float(b_ha3_shuf.loc[ln, b_hi_col]) < float(b_ha3_shuf.loc[ln, a_lo_col])
        )
    )
    ha3 = {
        "pass_count": ha3_pass,
        "pass_count_shuffled_genes": ha3_shuf_pass,
        "required": min_lines_pass,
        "criteria": f"matched > ALL S mismatched cov on BOTH ls_slope_ols_top50 AND {orth_metric_boot or 'slope_only'} ({ci_tag})",
        "verdict": "CONFIRMED" if ha3_pass >= min_lines_pass else "未确认",
        "per_line": ha3_per_line,
    }
    ha4 = _adjudicate_noninf_family(
        df_boot, primary_method, all_lines, eval_lines, ci_tag, max_lines_fail_ha4
    )
    v_group_a = (
        "CONFIRMED"
        if all(h["verdict"] == "CONFIRMED" for h in (ha1, ha2, ha3, ha4))
        else "未确认"
    )

    # ---- Group B: H-B1, H-B2, H-B3, H-B4 (Transfer_Mean_VarScaled) ----
    perm100_all = summary.get("tmvs_shuffled_var_perm100", {})
    if has_group_b and sub_gb is not None:
        hb1, hb2 = _adjudicate_slope_and_calr(
            sub_gb, all_lines, eval_lines, ci_tag, slope_threshold, min_lines_pass, ceilings=ceilings
        )
        hb3_per_line: Dict[str, Any] = {}
        hb3_pass = 0
        gb_cal_col = "ls_resid_corr_orth_cal_top50" if ("ls_resid_corr_orth_cal_top50" in sub_gb.columns and "ls_resid_corr_orth_cal_top50" in boot_metrics) else "ls_resid_corr_cal_top50"
        cal_lo_col = (
            f"{gb_cal_col}_resamp_ci975_lo"
            if ci_tag == "ci975" and f"{gb_cal_col}_resamp_ci975_lo" in sub_gb.columns
            else f"{gb_cal_col}_resamp_ci95_lo"
        )
        perm_sub_key = "r_orth_cal" if gb_cal_col == "ls_resid_corr_orth_cal_top50" else "r_cal"

        for ln in all_lines:
            is_conf = ln in eval_lines
            b_tmvs_s = df_boot[
                (df_boot["target_line"] == ln)
                & (df_boot["method_a"] == group_b_method)
                & (df_boot["metric"] == "ls_slope_ols_top50")
                & df_boot["method_b"].str.startswith("NegCtrl_TMVS_MismatchedVar(")
            ].set_index("method_b")
            b_tmvs_c = df_boot[
                (df_boot["target_line"] == ln)
                & (df_boot["method_a"] == group_b_method)
                & (df_boot["metric"] == gb_cal_col)
                & df_boot["method_b"].str.startswith("NegCtrl_TMVS_MismatchedVar(")
            ].set_index("method_b")

            m_s_lo = float(sub_gb.loc[ln, f"ls_slope_ols_top50_{ci_tag}_lo"])
            m_c_est = float(sub_gb.loc[ln, gb_cal_col])
            m_c_lo = float(sub_gb.loc[ln, cal_lo_col])
            by_ctrl: Dict[str, Any] = {}
            n_mis = 0
            all_mis_pass = True
            max_mis_hi = float("-inf")

            ctrl_names = sorted(set(b_tmvs_s.index.tolist()) | set(b_tmvs_c.index.tolist()))
            for ctrl_name in ctrl_names:
                s_est = float(b_tmvs_s.loc[ctrl_name, "mean_b"]) if ctrl_name in b_tmvs_s.index else float("nan")
                s_lo = float(b_tmvs_s.loc[ctrl_name, b_lo_col]) if ctrl_name in b_tmvs_s.index else float("nan")
                s_hi = float(b_tmvs_s.loc[ctrl_name, b_hi_col]) if ctrl_name in b_tmvs_s.index else float("nan")
                pass_s = bool(np.isfinite(m_s_lo) and np.isfinite(s_hi) and m_s_lo > s_hi)

                c_est = float(b_tmvs_c.loc[ctrl_name, "mean_b"]) if ctrl_name in b_tmvs_c.index else float("nan")
                c_hi = float(b_tmvs_c.loc[ctrl_name, b_hi_col]) if ctrl_name in b_tmvs_c.index else float("nan")
                # On 2-line dev (|S|=1), MismatchedVar(s0) reduces identically to Transfer_Mean (zero residual).
                zero_residual_null = (
                    np.isfinite(s_est)
                    and abs(s_est) <= 1e-12
                    and abs(s_lo) <= 1e-12
                    and abs(s_hi) <= 1e-12
                    and not np.isfinite(c_hi)
                )
                eff_c_hi = 0.0 if zero_residual_null else c_hi
                if np.isfinite(eff_c_hi):
                    max_mis_hi = max(max_mis_hi, eff_c_hi)
                pass_c = bool(np.isfinite(m_c_lo) and np.isfinite(eff_c_hi) and m_c_lo > eff_c_hi)
                by_ctrl[ctrl_name] = {
                    "ctrl_r_orth_cal": 0.0 if zero_residual_null else c_est,
                    f"ctrl_{gb_cal_col}_{ci_tag}_hi": eff_c_hi,
                    "ctrl_slope": s_est,
                    f"ctrl_slope_{ci_tag}_hi": s_hi,
                    "zero_residual_equivalence": zero_residual_null,
                    "pass_slope_info": pass_s,
                    "pass_cal_r": pass_c,
                }
                all_mis_pass = all_mis_pass and pass_c
                n_mis += 1

            # R18: 100-permutation shuffled target variance null (fallback to single df_boot row only in legacy tests)
            perm_entry = perm100_all.get(ln, {}).get(perm_sub_key) if isinstance(perm100_all.get(ln), dict) else None
            if perm_entry is not None:
                perm_p25 = float(perm_entry.get("p2_5", float("nan")))
                perm_med = float(perm_entry.get("median", float("nan")))
                perm_p975 = float(perm_entry.get("p97_5", float("nan")))
                perm_p9875 = float(perm_entry.get("p98_75", float("nan")))
                pass_perm = bool(np.isfinite(m_c_lo) and np.isfinite(perm_p975) and m_c_lo > perm_p975)
                by_ctrl["NegCtrl_TMVS_ShuffledVar_perm100"] = {
                    "n_perm": int(perm100_all[ln].get("n_perm", 100)),
                    "p2_5": perm_p25,
                    "median": perm_med,
                    "p97_5": perm_p975,
                    "p98_75": perm_p9875,
                    "pass_cal_r": pass_perm,
                }
            else:
                b_shuf_c = df_boot[
                    (df_boot["target_line"] == ln)
                    & (df_boot["method_a"] == group_b_method)
                    & (df_boot["metric"] == gb_cal_col)
                    & (df_boot["method_b"] == "NegCtrl_TMVS_ShuffledVar")
                ]
                perm_p975 = float(b_shuf_c.iloc[0][b_hi_col]) if not b_shuf_c.empty else float("nan")
                pass_perm = bool(np.isfinite(m_c_lo) and np.isfinite(perm_p975) and m_c_lo > perm_p975)
                by_ctrl["NegCtrl_TMVS_ShuffledVar_single_legacy"] = {
                    f"ctrl_{gb_cal_col}_{ci_tag}_hi": perm_p975,
                    "pass_cal_r": pass_perm,
                }

            bound_candidates = [x for x in (max_mis_hi, perm_p975) if np.isfinite(x)]
            max_ctrl_bound = max(bound_candidates) if bound_candidates else float("nan")
            ok_hb3 = bool(
                n_mis >= 1
                and all_mis_pass
                and pass_perm
                and np.isfinite(m_c_lo)
                and np.isfinite(max_ctrl_bound)
                and m_c_lo > max_ctrl_bound
            )
            if is_conf:
                hb3_pass += int(ok_hb3)
            hb3_per_line[ln] = {
                "confirmatory": is_conf,
                "metric_used": gb_cal_col,
                "matched_r_orth_cal": m_c_est,
                f"matched_cal_r_{ci_tag}_lo": m_c_lo,
                f"matched_slope_{ci_tag}_lo": m_s_lo,
                "n_mismatch_controls_checked": n_mis,
                "max_mismatch_ci_hi": max_mis_hi if np.isfinite(max_mis_hi) else float("nan"),
                "perm100_p97_5": perm_p975,
                "max_control_bound": max_ctrl_bound,
                "pass": ok_hb3,
                "by_control": by_ctrl,
            }
        hb3 = {
            "pass_count": hb3_pass,
            "required": min_lines_pass,
            "metric_used": gb_cal_col,
            "criteria": f"TMVS {gb_cal_col} {ci_tag} lo > max(ALL mismatched vars {ci_tag} hi, perm100 p97.5)",
            "verdict": "CONFIRMED" if hb3_pass >= min_lines_pass else "未确认",
            "per_line": hb3_per_line,
        }
        hb4 = _adjudicate_noninf_family(
            df_boot, group_b_method, all_lines, eval_lines, ci_tag, max_lines_fail_ha4
        )
        v_group_b = (
            "CONFIRMED"
            if all(h["verdict"] == "CONFIRMED" for h in (hb1, hb2, hb3, hb4))
            else "未确认"
        )
        vcc_band = {
            ln: {
                "slope": float(sub_gb.loc[ln, "ls_slope_ols_top50"]),
                f"{ci_tag}_lo": float(sub_gb.loc[ln, f"ls_slope_ols_top50_{ci_tag}_lo"]),
                f"{ci_tag}_hi": float(sub_gb.loc[ln, f"ls_slope_ols_top50_{ci_tag}_hi"]),
                "point_in_0.02_0.13": bool(0.02 <= float(sub_gb.loc[ln, "ls_slope_ols_top50"]) <= 0.13),
            }
            for ln in eval_lines
        }
    else:
        hb1 = {"verdict": "MISSING"}
        hb2 = {"verdict": "MISSING"}
        hb3 = {"verdict": "MISSING"}
        hb4 = {"verdict": "MISSING"}
        v_group_b = "MISSING"
        vcc_band = {}

    # ---- R4 & R17 Wording Checks (vs SrcCorrTgtVar & Transfer_Mean_VarScaled on r_orth_cal) ----
    r4_cal_metric = "ls_resid_corr_orth_cal_top50" if use_orth_boot else "ls_resid_corr_cal_top50"
    r4_sctv_per_line: Dict[str, Any] = {}
    r4_var_per_line: Dict[str, Any] = {}
    b_sctv_slope = df_boot[
        (df_boot["method_a"] == primary_method)
        & (df_boot["method_b"] == "NegCtrl_HybridProj_SrcCorrTgtVar")
        & (df_boot["metric"] == "ls_slope_ols_top50")
    ].set_index("target_line")
    b_sctv_cal = df_boot[
        (df_boot["method_a"] == primary_method)
        & (df_boot["method_b"] == "NegCtrl_HybridProj_SrcCorrTgtVar")
        & (df_boot["metric"] == r4_cal_metric)
    ].set_index("target_line")
    b_var_cal = df_boot[
        (df_boot["method_a"] == primary_method)
        & (df_boot["method_b"] == group_b_method)
        & (df_boot["metric"] == r4_cal_metric)
    ].set_index("target_line")

    for ln in all_lines:
        s_ok = False
        c_ok = False
        if ln in b_sctv_slope.index:
            rs_r = b_sctv_slope.loc[ln]
            s_ok = bool(
                rs_r[disj_col]
                and float(rs_r["delta"]) > 0.0
                and float(rs_r[b_hi_col]) < float(rs_r[a_lo_col])
            )
        if ln in b_sctv_cal.index:
            rc_r = b_sctv_cal.loc[ln]
            c_ok = bool(
                rc_r[disj_col]
                and float(rc_r["delta"]) > 0.0
                and float(rc_r[b_hi_col]) < float(rc_r[a_lo_col])
            )
        r4_sctv_per_line[ln] = {
            "confirmatory": ln in eval_lines,
            "slope_disjoint_win": s_ok,
            "cal_r_disjoint_win": c_ok,
            "both_disjoint_win": bool(s_ok and c_ok),
        }
        if ln in b_var_cal.index:
            rv_r = b_var_cal.loc[ln]
            r4_var_per_line[ln] = {
                "confirmatory": ln in eval_lines,
                "metric_used": r4_cal_metric,
                "delta_cal_r": float(rv_r["delta"]),
                f"paired_{ci_tag}_lo": float(rv_r[lo_col]),
                "paired_ci_gt_zero": bool(rv_r[gt_col]),
            }

    all_conf_beat_sctv = bool(eval_lines) and all(
        r4_sctv_per_line.get(ln, {}).get("both_disjoint_win", False) for ln in eval_lines
    )
    all_conf_beat_varscaled = bool(eval_lines) and all(
        r4_var_per_line.get(ln, {}).get("paired_ci_gt_zero", False) for ln in eval_lines
    )
    allowed_phrase = "细胞系特异的共波动结构" if all_conf_beat_sctv else "目标线对照统计量"

    # ---- R14 & R20 Mechanism Check (Hybrid_on_TMVS - Transfer_Mean_VarScaled r_orth_cal paired CI > 0) ----
    b_r14 = df_boot[
        (df_boot["method_a"] == "Hybrid_on_TMVS")
        & (df_boot["method_b"] == group_b_method)
        & (df_boot["metric"] == r4_cal_metric)
    ].set_index("target_line")
    r14_per_line: Dict[str, Any] = {}
    r14_pass = 0
    r14_pass_95 = 0
    for ln in all_lines:
        if ln in b_r14.index:
            rr = b_r14.loc[ln]
            ok_r14 = bool(rr[gt_col])
            ok_r14_95 = bool(rr["ci_gt_zero"])
            if ln in eval_lines:
                r14_pass += int(ok_r14)
                r14_pass_95 += int(ok_r14_95)
            r14_per_line[ln] = {
                "confirmatory": ln in eval_lines,
                "metric_used": r4_cal_metric,
                "delta_cal_r": float(rr["delta"]),
                f"{ci_tag}_lo": float(rr[lo_col]),
                "ci95_lo": float(rr["ci95_lo"]),
                f"pass_{ci_tag}": ok_r14,
                "pass_ci95": ok_r14_95,
            }
    r14_confirmed = bool(eval_lines) and (r14_pass_95 == len(eval_lines))
    r14_statement = (
        "共波动结构在边际方差之外提供了额外信息"
        if r14_confirmed
        else "在本数据上，共波动结构没有超出边际方差的可检测增量"
    )

    # ---- R18: Single-scalar GlobalScale summary ----
    gs_sub = df_lines[df_lines["method"] == "Transfer_Mean_GlobalScale"].set_index("target_line")
    gs_per_line: Dict[str, Any] = {}
    for ln in all_lines:
        if ln in gs_sub.index:
            rg = gs_sub.loc[ln]
            gs_per_line[ln] = {
                "confirmatory": ln in eval_lines,
                "ls_resid_corr_top50": float(rg["ls_resid_corr_top50"]),
                "ls_resid_corr_cal_top50": float(rg["ls_resid_corr_cal_top50"]),
                "ls_resid_corr_orth_top50": float(rg.get("ls_resid_corr_orth_top50", float("nan"))),
                "ls_resid_corr_orth_cal_top50": float(rg.get("ls_resid_corr_orth_cal_top50", float("nan"))),
                "ls_slope_ols_top50": float(rg["ls_slope_ols_top50"]),
            }
    conf_gs_raw = [
        gs_per_line[ln]["ls_resid_corr_top50"]
        for ln in eval_lines
        if ln in gs_per_line and np.isfinite(gs_per_line[ln]["ls_resid_corr_top50"])
    ]
    conf_gs_orth = [
        gs_per_line[ln]["ls_resid_corr_orth_cal_top50"]
        for ln in eval_lines
        if ln in gs_per_line and np.isfinite(gs_per_line[ln]["ls_resid_corr_orth_cal_top50"])
    ]
    gs_summary_str = (
        f"单个标量 Transfer_Mean_GlobalScale 即可得到原始残差相关 "
        f"[{min(conf_gs_raw):+.3f}, {max(conf_gs_raw):+.3f}]"
        if conf_gs_raw
        else "Transfer_Mean_GlobalScale 未运行"
    )
    if conf_gs_orth:
        gs_summary_str += f"（去幅度 r_orth_cal 为 [{min(conf_gs_orth):+.3f}, {max(conf_gs_orth):+.3f}]）"

    var_decomp = summary.get("variance_decomposition", {})

    return {
        "code_freeze_hash": summary.get("code_freeze_hash", "UNFROZEN_DEV"),
        "bootstrap_seed": summary.get("bootstrap_seed", 20261121),
        "n_total_lines": len(all_lines),
        "n_lines": n_eval,
        "confirmatory_lines": eval_lines,
        "secondary_lines": secondary_lines,
        "primary_method": primary_method,
        "group_b_method": group_b_method,
        "ci_level_used": resolved_ci_level,
        "ci_tag_used": ci_tag,
        "primary_corr_metric": r4_cal_metric,
        "OVERALL_VERDICT": v_group_a,
        "Group_A_Verdict": v_group_a,
        "Group_B_Verdict": v_group_b,
        "R13_Dual_Summary": f"Group A ({primary_method}): {v_group_a} | Group B ({group_b_method}): {v_group_b} (Bonferroni {resolved_ci_level*100:.1f}% CIs)",
        "H-A1": ha1,
        "H-A2": ha2,
        "H-A3": ha3,
        "H-A4": ha4,
        "H-B1": hb1,
        "H-B2": hb2,
        "H-B3": hb3,
        "H-B4": hb4,
        "R4_Wording_Checks": {
            "all_confirmatory_beat_SrcCorrTgtVar_both_disjoint": all_conf_beat_sctv,
            "allowed_cofluctuation_phrase": allowed_phrase,
            "all_confirmatory_beat_TransferMeanVarScaled_cal_r": all_conf_beat_varscaled,
            "marginal_variance_explanation_note": (
                "NONE (beats Transfer_Mean_VarScaled)"
                if all_conf_beat_varscaled
                else "增益可由边际方差解释"
            ),
            "vs_SrcCorrTgtVar_by_line": r4_sctv_per_line,
            "vs_TransferMeanVarScaled_by_line": r4_var_per_line,
        },
        "R14_Mechanism_Check": {
            "pass_count": r14_pass,
            "pass_count_ci95": r14_pass_95,
            "required": len(eval_lines),
            "confirmed": r14_confirmed,
            "metric_used": r4_cal_metric,
            "allowed_statement": r14_statement,
            "per_line": r14_per_line,
        },
        "R18_GlobalScale_Summary": {
            "summary_statement": gs_summary_str,
            "per_line": gs_per_line,
        },
        "R19_Variance_Decomposition": var_decomp,
        "R12_VCC_Top10_Band_Check": {
            "r6_caveat": "Arc VCC 2026 Wrap-Up (0.02-0.13) slope estimator formula is undisclosed and cell-eval2 does not implement it; OLS slope = r * (sd_pred / sd_obs) is scale-dependent.",
            "per_confirmatory_line": vcc_band,
        },
        "R20_Novelty_Constraint": (
            "FCT is equivalent to VCC 2026 CellularAI (mu_t * 2^log2FC); CASTER performed supervised control-variance scaling. "
            "TMVS is framed strictly as a parameter-free variance-ratio variant of expression-aware transfer (not a new method); "
            "contribution is the pre-registered 3-layer decomposition of recoverable zero-shot line-specific residuals into "
            "(1) line-wide amplitude, (2) gene marginal variance, and (3) off-diagonal co-fluctuations."
        ),
    }


def render_markdown_report(adj: Dict[str, Any], df_lines: pd.DataFrame, summary: Dict[str, Any], tag: str) -> str:
    ci_tag = adj.get("ci_tag_used", "ci975")
    ci_pct = "97.5%" if ci_tag == "ci975" else "95%"
    r4 = adj.get("R4_Wording_Checks", {})
    r14 = adj.get("R14_Mechanism_Check", {})
    r18_gs = adj.get("R18_GlobalScale_Summary", {})
    r19_vd = adj.get("R19_Variance_Decomposition", {})
    conf_set = set(adj["confirmatory_lines"])

    ha4_sub = ", ".join(f"{k}={v['n_inferior']}/{adj['n_lines']}" for k, v in adj["H-A4"]["by_metric"].items())

    lines_md: List[str] = [
        f"# Line-Specific Zero-Shot Perturbation Evaluation Report (`{tag}`)",
        "",
        f"- **Code Freeze Commit (C1)**: `{adj.get('code_freeze_hash', summary.get('code_freeze_hash', 'UNFROZEN_DEV'))}` | **Bootstrap Seed (C2)**: `{adj.get('bootstrap_seed', summary.get('bootstrap_seed', 20261121))}`",
        f"- **Confirmatory Lines ({adj['n_lines']})**: {', '.join(adj['confirmatory_lines'])}",
        f"- **Secondary Lines ({len(adj['secondary_lines'])})**: {', '.join(adj['secondary_lines']) if adj['secondary_lines'] else 'None'}",
        f"- **Dual Confirmatory Summary (R13, {ci_pct} Bonferroni CIs)**: **{adj['R13_Dual_Summary']}**",
        f"- **Group A Primary (`{adj['primary_method']}`)**: **{adj['Group_A_Verdict']}**",
        f"- **Group B Primary (`{adj['group_b_method']}`)**: **{adj['Group_B_Verdict']}**",
        f"- **R17 Primary Correlation Metric**: `{adj.get('primary_corr_metric', 'ls_resid_corr_orth_cal_top50')}` (amplitude-orthogonalized; raw `r_cal` reported as secondary)",
        f"- **R18 Global-Scale Scalar Baseline**: {r18_gs.get('summary_statement', 'N/A')}",
        f"- **R4 Co-fluctuation Structure Claim Allowed**: `{r4.get('allowed_cofluctuation_phrase', '')}` (`all beat SrcCorrTgtVar disjoint={r4.get('all_confirmatory_beat_SrcCorrTgtVar_both_disjoint')}`)",
        f"- **R4 Marginal Variance Note**: `{r4.get('marginal_variance_explanation_note', '')}`",
        f"- **R14 Mechanism Statement (`Hybrid_on_TMVS` vs `TMVS` on `r_orth_cal`, 95% CI)**: `{r14.get('allowed_statement', '')}` (`95%: {r14.get('pass_count_ci95', r14.get('pass_count'))}/{r14.get('required')}`, `97.5%: {r14.get('pass_count')}/{r14.get('required')}`)",
        f"- **R20 Prior-Art Novelty Scope**: {adj.get('R20_Novelty_Constraint', '')}",
        "",
        f"## Hypothesis Adjudication Summary (Confirmatory Cohort, {ci_pct} CIs)",
        "",
        "| Hypothesis | Group | Pass / Confirmatory | Threshold | Verdict |",
        "|---|---|---:|---|---|",
        f"| **H-A1** (Slope {ci_pct} CI lo > 0.13) | A (`Hybrid_ProjShared`) | {adj['H-A1']['pass_count']}/{adj['n_lines']} | >= {adj['H-A1']['required']} | **{adj['H-A1']['verdict']}** |",
        f"| **H-A2** (Joint-boot `r_orth_cal` {ci_pct} CI lo > 0, valid>=95%, rho_AB_orth>0) | A (`Hybrid_ProjShared`) | {adj['H-A2']['pass_count_resamp_ceil']}/{adj['n_lines']} (fixed95: {adj['H-A2']['pass_count_fixed_ceil']}/{adj['n_lines']}) | >= {adj['H-A2']['required']} | **{adj['H-A2']['verdict']}** |",
        f"| **H-A3** (Matched > ALL source mismatched cov on slope AND `r_orth_cal`) | A (`Hybrid_ProjShared`) | {adj['H-A3']['pass_count']}/{adj['n_lines']} | >= {adj['H-A3']['required']} | **{adj['H-A3']['verdict']}** |",
        f"| **H-A4** (Strict AND non-inferiority across 6 conventions) | A (`Hybrid_ProjShared`) | {adj['H-A4']['max_n_inferior']}/{adj['n_lines']} inf ({ha4_sub}; union={adj['H-A4']['n_inferior_union']}) | <= {adj['H-A4']['max_inferior_allowed']} inf | **{adj['H-A4']['verdict']}** |",
    ]
    if adj.get("Group_B_Verdict") != "MISSING":
        hb4_sub = ", ".join(f"{k}={v['n_inferior']}/{adj['n_lines']}" for k, v in adj["H-B4"]["by_metric"].items())
        lines_md.extend([
            f"| **H-B1** (Slope {ci_pct} CI lo > 0.13) | B (`Transfer_Mean_VarScaled`) | {adj['H-B1']['pass_count']}/{adj['n_lines']} | >= {adj['H-B1']['required']} | **{adj['H-B1']['verdict']}** |",
            f"| **H-B2** (Joint-boot `r_orth_cal` {ci_pct} CI lo > 0, valid>=95%, rho_AB_orth>0) | B (`Transfer_Mean_VarScaled`) | {adj['H-B2']['pass_count_resamp_ceil']}/{adj['n_lines']} | >= {adj['H-B2']['required']} | **{adj['H-B2']['verdict']}** |",
            f"| **H-B3** (Matched `r_orth_cal` lo > max(ALL mismatched vars CI hi, 100-perm p97.5)) | B (`Transfer_Mean_VarScaled`) | {adj['H-B3']['pass_count']}/{adj['n_lines']} | >= {adj['H-B3']['required']} | **{adj['H-B3']['verdict']}** |",
            f"| **H-B4** (Strict AND non-inferiority across 6 conventions) | B (`Transfer_Mean_VarScaled`) | {adj['H-B4']['max_n_inferior']}/{adj['n_lines']} inf ({hb4_sub}; union={adj['H-B4']['n_inferior_union']}) | <= {adj['H-B4']['max_inferior_allowed']} inf | **{adj['H-B4']['verdict']}** |",
        ])

    lines_md.extend([
        "",
        f"## Main Table: Primary View (A-T2 Own-Target-Zeroed Control Reference, {ci_pct} CIs, R5/R6/R11/R13/R17/R18)",
        "",
        f"| Cohort | Target Line | Method | N | Slope OLS (Top50) | {ci_pct} CI | r_orth_cal | r_orth_cal^2 | r_orth_cal^2 Joint {ci_pct} CI | Resid r_orth | r_cal (raw) | Resid r (raw) | sd_pred/sd_obs | Deming (delta=1) | Obs~Pred | Ctrl disc_raw | Ctrl de50_r | MSE |",
        "|---|---|---|---:|---:|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ])
    main_prefixes = (
        "Hybrid_ProjShared",
        "Hybrid_Raw",
        "CIPHER_Linear",
        "Transfer_Mean",
        "Transfer_Mean_GlobalScale",
        "Transfer_Mean_VarScaled",
        "Hybrid_on_TMVS",
        "FCT",
        "NegCtrl_HybridProj_SrcCorrTgtVar",
        "NegCtrl_HybridProj_ShuffledGenes",
        "NegCtrl_HybridProj_DiagCov",
        "NegCtrl_TMVS_ShuffledVar",
        "Source_Grand_Mean",
        "NegCtrl_HybridProj_MismatchedCov(",
        "NegCtrl_TMVS_MismatchedVar(",
    )
    slope_lo_k = f"ls_slope_ols_top50_{ci_tag}_lo" if f"ls_slope_ols_top50_{ci_tag}_lo" in df_lines.columns else "ls_slope_ols_top50_ci95_lo"
    slope_hi_k = f"ls_slope_ols_top50_{ci_tag}_hi" if f"ls_slope_ols_top50_{ci_tag}_hi" in df_lines.columns else "ls_slope_ols_top50_ci95_hi"
    sq_orth_lo_k = (
        f"ls_resid_corr_orth_cal_top50_sq_resamp_{ci_tag}_lo"
        if f"ls_resid_corr_orth_cal_top50_sq_resamp_{ci_tag}_lo" in df_lines.columns
        else (
            "ls_resid_corr_orth_cal_top50_sq_resamp_ci95_lo"
            if "ls_resid_corr_orth_cal_top50_sq_resamp_ci95_lo" in df_lines.columns
            else "ls_resid_corr_cal_top50_sq_resamp_ci95_lo"
        )
    )
    sq_orth_hi_k = (
        f"ls_resid_corr_orth_cal_top50_sq_resamp_{ci_tag}_hi"
        if f"ls_resid_corr_orth_cal_top50_sq_resamp_{ci_tag}_hi" in df_lines.columns
        else (
            "ls_resid_corr_orth_cal_top50_sq_resamp_ci95_hi"
            if "ls_resid_corr_orth_cal_top50_sq_resamp_ci95_hi" in df_lines.columns
            else "ls_resid_corr_cal_top50_sq_resamp_ci95_hi"
        )
    )

    for _, r in df_lines.iterrows():
        m_name = str(r["method"])
        if not (m_name in main_prefixes or any(m_name.startswith(p) for p in ("NegCtrl_HybridProj_MismatchedCov(", "NegCtrl_TMVS_MismatchedVar("))):
            continue
        cohort = "Confirmatory" if str(r["target_line"]) in conf_set else "Secondary"
        ci_str = _fmt_ci(r[slope_lo_k], r[slope_hi_k])
        orth_cal_r = r.get("ls_resid_corr_orth_cal_top50", float("nan"))
        orth_sq = r.get("ls_resid_corr_orth_cal_top50_sq", float("nan"))
        orth_sq_ci = _fmt_ci(r.get(sq_orth_lo_k, float("nan")), r.get(sq_orth_hi_k, float("nan")), fmt=".4f")
        orth_r_raw = r.get("ls_resid_corr_orth_top50", float("nan"))
        lines_md.append(
            f"| {cohort} | {r['target_line']} | `{m_name}` | {int(r['n_eval_targets'])} | "
            f"{r['ls_slope_ols_top50']:+.4f} | {ci_str} | {orth_cal_r:+.4f} | "
            f"{orth_sq:.4f} | {orth_sq_ci} | {orth_r_raw:+.4f} | "
            f"{r['ls_resid_corr_cal_top50']:+.4f} | {r['ls_resid_corr_top50']:+.4f} | "
            f"{r.get('ls_sd_ratio_top50', float('nan')):.4f} | {r.get('ls_slope_deming_top50', float('nan')):+.4f} | "
            f"{r.get('ls_slope_obs_on_pred_top50', float('nan')):+.4f} | "
            f"{r['ctrl_disc_raw']:.4f} | {r['ctrl_de50_r']:.4f} | {r['mse']:.4f} |"
        )

    # R18 Table: 100-Permutation TMVS Shuffled-Variance Null vs H-B3 Threshold
    perm100 = summary.get("tmvs_shuffled_var_perm100", {})
    if perm100:
        hb3_pl = adj.get("H-B3", {}).get("per_line", {})
        lines_md.extend([
            "",
            f"## R18: TMVS Shuffled-Variance 100-Permutation Null vs H-B3 ({ci_pct} CI Threshold)",
            "",
            "| Target Line | TMVS r_orth_cal | TMVS CI lo | Perm100 p2.5 | Perm100 median | Perm100 p97.5 | Perm100 max | Max Mismatch CI hi | H-B3 Pass |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---|",
        ])
        for ln in sorted(perm100.keys()):
            po = perm100[ln].get("r_orth_cal", {})
            hln = hb3_pl.get(ln, {})
            lines_md.append(
                f"| {ln} | {hln.get('matched_r_orth_cal', float('nan')):+.4f} | "
                f"{hln.get(f'matched_cal_r_{ci_tag}_lo', float('nan')):+.4f} | "
                f"{po.get('p2_5', float('nan')):+.4f} | {po.get('median', float('nan')):+.4f} | "
                f"{po.get('p97_5', float('nan')):+.4f} | {po.get('max', float('nan')):+.4f} | "
                f"{hln.get('max_mismatch_ci_hi', float('nan')):+.4f} | {hln.get('pass', False)} |"
            )

    # R19 Table: Descriptive 3-Layer Variance Decomposition of Residual Delta_t - Delta_src
    if r19_vd:
        lines_md.extend([
            "",
            "## R19: Descriptive 3-Layer Residual Variance Decomposition (`R = Delta_t - Delta_src`, PREREG_A v4 §5.4)",
            "",
            f"| Target Line | Mask | N | Elig Genes | beta_ols_int [95% CI] | beta_mult | L1 Line Gain [95% CI] | L2 Per-Gene Gain [95% CI] | L3 Remaining [95% CI] |",
            "|---|---|---:|---:|---|---:|---|---|---|",
        ])
        for ln in sorted(r19_vd.keys()):
            for mode in ("resp_top50", "all_nontarget"):
                dm = r19_vd[ln].get(mode)
                if not dm:
                    continue
                bi = dm["beta_ols_int"]
                bm = dm["beta_mult"]
                l1 = dm["frac_line_gain"]
                l2 = dm["frac_gene_gain"]
                l3 = dm["frac_remaining"]
                lines_md.append(
                    f"| {ln} | `{mode}` | {dm.get('n_perts', 0)} | {dm.get('n_eligible_genes', 0)} (C_dist>={dm.get('min_gene_obs_for_gene_layer', 5)}) | "
                    f"{bi['estimate']:.4f} {_fmt_ci(bi['ci95'][0], bi['ci95'][1], fmt='.4f')} | "
                    f"{bm['estimate']:.4f} | "
                    f"{l1['estimate']:.4f} {_fmt_ci(l1['ci95'][0], l1['ci95'][1], fmt='.4f')} | "
                    f"{l2['estimate']:.4f} {_fmt_ci(l2['ci95'][0], l2['ci95'][1], fmt='.4f')} | "
                    f"{l3['estimate']:.4f} {_fmt_ci(l3['ci95'][0], l3['ci95'][1], fmt='.4f')} |"
                )

    # S6 Table: Split-half noise ceiling selection bias (resp_top50 vs all_nontarget, orth + raw)
    ceilings = summary.get("pert_noise_ceilings", {})
    if ceilings:
        lines_md.extend([
            "",
            "## S6: Split-Half Noise Ceiling (`resp_top50` vs Unselected `all_nontarget`, Amplitude-Orthogonalized & Raw)",
            "",
            "| Target Line | N half-valid | rho_AB_orth (top50) | Ceil r_orth (top50) | rho_AB (top50) | Rep Slope (top50) | rho_AB_orth (all) | rho_AB (all) | TMVS r_orth_cal (all) | Hybrid r_orth_cal (all) |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ])
        for ln in sorted(ceilings.keys()):
            c = ceilings[ln]
            t50 = c.get("resp_top50", {})
            allm = c.get("all_nontarget", {})
            sub_ln = df_lines[df_lines["target_line"] == ln].set_index("method")
            tmvs_all_o = float(sub_ln.loc["Transfer_Mean_VarScaled"].get("ls_resid_corr_orth_cal_all", float("nan"))) if "Transfer_Mean_VarScaled" in sub_ln.index else float("nan")
            hyb_all_o = float(sub_ln.loc[adj["primary_method"]].get("ls_resid_corr_orth_cal_all", float("nan"))) if adj["primary_method"] in sub_ln.index else float("nan")
            lines_md.append(
                f"| {ln} | {c.get('n_half_valid_targets', 0)} | "
                f"{t50.get('rho_half_ab_orth', float('nan')):.4f} | {t50.get('ceiling_r_orth', float('nan')):.4f} | "
                f"{t50.get('rho_half_ab', float('nan')):.4f} | {t50.get('slope_rep_mean', float('nan')):.4f} | "
                f"{allm.get('rho_half_ab_orth', float('nan')):.4f} | {allm.get('rho_half_ab', float('nan')):.4f} | "
                f"{tmvs_all_o:+.4f} | {hyb_all_o:+.4f} |"
            )

    # S5 Table: Top-75% expressed genes by target control variance
    s5 = summary.get("s5_expressed_top75_var", {})
    if s5:
        lines_md.extend([
            "",
            "## S5: Expressed Genes Sensitivity (Drop Bottom 25% Target Control Variance Genes)",
            "",
            f"| Target Line | Genes Kept | rho_AB_orth (S5) | Method | Slope OLS (S5) | r_orth_cal (S5) | r_orth_cal^2 (S5) | Joint {ci_pct} CI | r_cal raw (S5) |",
            "|---|---:|---:|---|---:|---:|---:|---|---:|",
        ])
        ci_orth_key = "ls_resid_corr_orth_cal_top50_resamp_ci975" if ci_tag == "ci975" else "ls_resid_corr_orth_cal_top50_resamp_ci95"
        ci_raw_key = "ls_resid_corr_cal_top50_resamp_ci975" if ci_tag == "ci975" else "ls_resid_corr_cal_top50_resamp_ci95"
        for ln in sorted(s5.keys()):
            s5_ln = s5[ln]
            for m_name, mvals in s5_ln.get("methods", {}).items():
                ci_v = mvals.get(ci_orth_key, mvals.get(ci_raw_key, [float("nan"), float("nan")]))
                lines_md.append(
                    f"| {ln} | {s5_ln.get('n_genes_kept')}/{s5_ln.get('n_genes_total')} | "
                    f"{s5_ln.get('rho_half_ab_orth_top50', s5_ln.get('rho_half_ab_top50', float('nan'))):.4f} | `{m_name}` | "
                    f"{mvals.get('ls_slope_ols_top50', float('nan')):+.4f} | "
                    f"{mvals.get('ls_resid_corr_orth_cal_top50', float('nan')):+.4f} | "
                    f"{mvals.get('ls_resid_corr_orth_cal_top50_sq', float('nan')):.4f} | "
                    f"{_fmt_ci(ci_v[0], ci_v[1])} | "
                    f"{mvals.get('ls_resid_corr_cal_top50', float('nan')):+.4f} |"
                )

    lines_md.extend([
        "",
        "## Appendix: Alternative H-A4/H-B4 Discrimination & Centroid Conventions (R5)",
        "",
        "| Target Line | Method | Ctrl disc_raw (A-T2) | Ctrl disc_raw_panel (A-T3) | Cent disc_raw | Cent disc_raw_panel | Cent de50_r | pds_cosine_panel | Slope OLS (All) |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ])
    for _, r in df_lines.iterrows():
        m_name = str(r["method"])
        if m_name.startswith("NegCtrl_HybridProj_MismatchedCov(") or m_name.startswith("NegCtrl_TMVS_MismatchedVar("):
            continue
        if m_name not in main_prefixes:
            continue
        lines_md.append(
            f"| {r['target_line']} | `{m_name}` | "
            f"{r['ctrl_disc_raw']:.4f} | {r.get('ctrl_disc_raw_panel', float('nan')):.4f} | "
            f"{r['cent_disc_raw']:.4f} | {r.get('cent_disc_raw_panel', float('nan')):.4f} | "
            f"{r['cent_de50_r']:.4f} | {r.get('pds_cosine_panel', float('nan')):.4f} | "
            f"{r['ls_slope_ols_all']:+.4f} |"
        )
    return "\n".join(lines_md) + "\n"


def main() -> None:
    p = argparse.ArgumentParser(description="Adjudicate line-specific evaluation results.")
    p.add_argument("--results-dir", type=Path, default=Path(__file__).resolve().parents[1] / "results")
    p.add_argument("--tag", required=True, help="e.g. dev_log1p or holdout_log1p")
    p.add_argument("--primary-method", default="Hybrid_ProjShared")
    p.add_argument("--group-b-method", default="Transfer_Mean_VarScaled")
    p.add_argument("--min-lines-pass", type=int, default=None)
    p.add_argument("--max-lines-fail-ha4", type=int, default=None)
    p.add_argument("--ci-level", type=float, default=0.975)
    p.add_argument("--confirmatory-lines", type=str, default=None, help="Optional comma-separated confirmatory line names.")
    args = p.parse_args()

    df_lines = pd.read_csv(args.results_dir / f"line_metrics_{args.tag}.csv")
    df_boot = pd.read_csv(args.results_dir / f"bootstrap_ci_{args.tag}.csv")
    summary = json.loads((args.results_dir / f"summary_{args.tag}.json").read_text(encoding="utf-8"))

    conf_list = (
        [x.strip() for x in args.confirmatory_lines.split(",") if x.strip()]
        if args.confirmatory_lines
        else None
    )
    adj = adjudicate_results(
        df_lines=df_lines,
        df_boot=df_boot,
        summary=summary,
        primary_method=args.primary_method,
        group_b_method=args.group_b_method,
        min_lines_pass=args.min_lines_pass,
        max_lines_fail_ha4=args.max_lines_fail_ha4,
        ci_level=args.ci_level,
        confirmatory_lines=conf_list,
    )
    (args.results_dir / f"adjudication_{args.tag}.json").write_text(
        json.dumps(adj, indent=2), encoding="utf-8"
    )
    md = render_markdown_report(adj, df_lines, summary, args.tag)
    (args.results_dir / f"REPORT_{args.tag}.md").write_text(md, encoding="utf-8")
    print(f"Wrote adjudication_{args.tag}.json and REPORT_{args.tag}.md", flush=True)


if __name__ == "__main__":
    main()
