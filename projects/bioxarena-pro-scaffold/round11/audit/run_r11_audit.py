#!/usr/bin/env python3
"""
Independent Round 11 Group D Audit Runner, Report Generator, & R10-E Numerical Verifier (A-T3).

Strictly NEVER imports or inspects `round11/src/` or `round11/results/`.
Modes:
  1. `--mode dev_r9r10` (pre-freeze safe testing strictly on R9/R10 unblinded lines):
       Refuses any access to `round11/shared/` paths. Supports optional verification
       against `round10/exploratory/baseline_ablation_*.csv` (`seed=20261122`) as well as
       full `PREREG_D` dry-run (`seed=20261123`) on R9/R10 lines.
  2. `--mode holdout_r11` (post-freeze Round 11 evaluation):
       Enforces `PREREG_D_COMMIT.txt` + `CODE_FREEZE_R11.txt`, reads `round11/shared/lines/<line>`
       + `CONFIRMATORY_LINES_R11.txt`, adjudicates `P1` (`K_conf < 2`), Claim 1 (`H-D1 + A1`, `H-D2`,
       `H-D3`), Claim 2 (`H-D4a/b`), `H-D_util`, `SA-8` (`resp_top50_split`), and `H-D5` (`feng_ipsc`).
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

AUDIT_DIR = Path(__file__).resolve().parent
if str(AUDIT_DIR) not in sys.path:
    sys.path.insert(0, str(AUDIT_DIR))

from independent_r11_loader import (
    R9_LINE_NAMES,
    R9_LINE_SET,
    R10_LINES_DIR,
    R10_SOURCE_EXT_NAMES,
    R11_CANDIDATE_CONFIRMATORY,
    R11_LINES_DIR,
    R11_SECONDARY_ONLY,
    R11_STATUS_DIR,
    R9_LINES_DIR,
    AuditLineData,
    CanonicalHGNC,
    enforce_r11_gate,
    load_line_dataset,
    load_r9_common_genes,
    parse_confirmatory_lines_r11,
)
from independent_r11_metrics import (
    DEFAULT_R11_SEED,
    adjudicate_round11_group_d,
    evaluate_cross_mask_cal_diff,
    evaluate_fct_perm_null,
    evaluate_hd_util_noninferiority,
    evaluate_joint_residuals,
    evaluate_paired_cal_diff,
    evaluate_r19_gain_decomposition,
    evaluate_s7_confounder_controls,
    evaluate_sa8_perm_null,
    evaluate_sa8_resp_top50_split,
)
from independent_r11_predictors import (
    align_source_and_target,
    build_predictions,
    deterministic_subsample,
    sample_covariance,
)


def _fmt(val: Any, digits: int = 4) -> str:
    try:
        fv = float(val)
        return f"{fv:+.{digits}f}" if np.isfinite(fv) else "nan"
    except (TypeError, ValueError):
        return str(val)


def evaluate_target_line_r11(
    target_data: AuditLineData,
    source_list: Sequence[AuditLineData],
    only_reproducible: bool = True,
    shrink: float = 0.0,
    cov_max_cells: Optional[int] = None,
    n_boot: int = 4000,
    n_perm: int = 100,
    seed: int = DEFAULT_R11_SEED,
    include_sa3_anchors: bool = True,
    include_exploratory_e: bool = True,
    include_r19: bool = True,
    include_s7: bool = True,
    precomputed_covs: Optional[Dict[str, np.ndarray]] = None,
    precomputed_gene_orders: Optional[Dict[str, Sequence[str]]] = None,
) -> Dict[str, Any]:
    """Run complete independent Round 11 audit on one target line against its source pool."""
    cohort = align_source_and_target(source_list, target_data, only_reproducible=only_reproducible)
    preds, fc_bar, grand_mean, const_stats = build_predictions(
        cohort=cohort,
        shrink=shrink,
        cov_max_cells=cov_max_cells,
        seed=42,
        include_sa3_anchors=include_sa3_anchors,
        include_exploratory_e=include_exploratory_e,
        precomputed_covs=precomputed_covs,
        precomputed_gene_orders=precomputed_gene_orders,
    )
    t_mean = preds["transferred_mean"]

    # Evaluate primary (`all_nontarget`) & secondary (`resp_top50`) across all predictors & controls
    dual_mask_methods = {
        "fct",
        "e1_const_med",
        "fct_full_nt",
        "e2_const_p10",
        "e2_const_p90",
        "e3_src_mean_ctrl",
        "e5_rank_matched_mu",
        "e6_log_fc_bar",
        "e7_fc_minus_1",
        "transfer_mean_global_scale",
        "transfer_mean_var_scaled",
        "hybrid",
        "hybrid_on_tmvs",
        "cipher_linear",
        "source_grand_mean",
    } | {f"fct_mismatch::{s}" for s in cohort.source_names}

    method_eval: Dict[str, Dict[str, Any]] = {}
    for m_name, p_mat in preds.items():
        m_all = evaluate_joint_residuals(
            delta_obs=cohort.delta_obs,
            pred=p_mat,
            t_mean=t_mean,
            target_cols=cohort.target_cols,
            mode="all_nontarget",
            n_boot=n_boot,
            seed=seed,
            half_a_obs=cohort.half_a_obs,
            half_b_obs=cohort.half_b_obs,
            half_valid=cohort.half_valid,
        )
        m_top50 = (
            evaluate_joint_residuals(
                delta_obs=cohort.delta_obs,
                pred=p_mat,
                t_mean=t_mean,
                target_cols=cohort.target_cols,
                mode="resp_top50",
                n_boot=n_boot,
                seed=seed,
                half_a_obs=cohort.half_a_obs,
                half_b_obs=cohort.half_b_obs,
                half_valid=cohort.half_valid,
            )
            if m_name in dual_mask_methods
            else {}
        )
        method_eval[m_name] = {"all_nontarget": m_all, "resp_top50": m_top50}

    fct_all = method_eval["fct"]["all_nontarget"]
    fct_top50 = method_eval["fct"]["resp_top50"]
    e1_all = method_eval["e1_const_med"]["all_nontarget"]
    e1_top50 = method_eval["e1_const_med"]["resp_top50"]

    # Permutations under BOTH masks (`all_nontarget` for H-D2, `resp_top50` for SA-7 / H-D4 context)
    perm100_all = evaluate_fct_perm_null(
        delta_obs=cohort.delta_obs,
        t_mean=t_mean,
        fc_bar=fc_bar,
        tgt_ctrl_mean_libnorm=cohort.tgt_ctrl_mean_libnorm,
        target_cols=cohort.target_cols,
        ceil_raw=float(fct_all.get("noise_ceiling_r", 0.0)),
        ceil_orth=float(fct_all.get("noise_ceiling_r_orth", 0.0)),
        mode="all_nontarget",
        n_perm=n_perm,
    )
    perm100_top50 = evaluate_fct_perm_null(
        delta_obs=cohort.delta_obs,
        t_mean=t_mean,
        fc_bar=fc_bar,
        tgt_ctrl_mean_libnorm=cohort.tgt_ctrl_mean_libnorm,
        target_cols=cohort.target_cols,
        ceil_raw=float(fct_top50.get("noise_ceiling_r", 0.0)),
        ceil_orth=float(fct_top50.get("noise_ceiling_r_orth", 0.0)),
        mode="resp_top50",
        n_perm=n_perm,
    )

    # Mismatched controls under primary `all_nontarget` & secondary `resp_top50`
    mis_all = {
        s: method_eval[f"fct_mismatch::{s}"]["all_nontarget"]
        for s in cohort.source_names
        if f"fct_mismatch::{s}" in method_eval
    }
    mis_top50 = {
        s: method_eval[f"fct_mismatch::{s}"]["resp_top50"]
        for s in cohort.source_names
        if f"fct_mismatch::{s}" in method_eval
    }
    mis_all_ci975_his = [float(v.get("r_orth_cal_ci975_hi", float("nan"))) for v in mis_all.values()]
    max_mis_all_ci975_hi = float(np.nanmax(mis_all_ci975_his)) if mis_all_ci975_his else float("nan")
    perm_all_p975 = float(perm100_all["r_orth_cal"]["p975"])

    # H-D1 (`all_nontarget` + Addendum A1 validity conditions: valid_frac >= 0.95 & rho_ab_orth > 0)
    fct_all_lo = float(fct_all.get("r_orth_cal_ci975_lo", float("nan")))
    fct_all_vf = float(fct_all.get("r_orth_cal_valid_frac", 0.0))
    fct_all_rho = float(fct_all.get("rho_ab_orth", float("nan")))
    hd1_ci975_gt_0_only = bool(np.isfinite(fct_all_lo) and fct_all_lo > 0.0)
    hd1_valid_frac_ok = bool(np.isfinite(fct_all_vf) and fct_all_vf >= 0.95)
    hd1_rho_ab_ok = bool(np.isfinite(fct_all_rho) and fct_all_rho > 0.0)
    hd1_pass = bool(hd1_ci975_gt_0_only and hd1_valid_frac_ok and hd1_rho_ab_ok)
    hd1_a1_deviation_flag = bool(hd1_ci975_gt_0_only and not (hd1_valid_frac_ok and hd1_rho_ab_ok))

    # H-D2 (`all_nontarget` FCT CI_lo > max(97.5% CI_hi mismatched, Perm100 Q_0.975))
    hd2_pass = bool(
        np.isfinite(fct_all_lo)
        and np.isfinite(max_mis_all_ci975_hi)
        and np.isfinite(perm_all_p975)
        and (fct_all_lo > max(max_mis_all_ci975_hi, perm_all_p975))
    )

    # H-D3 (paired FCT - E1_ConstMed under `all_nontarget`) + R10-E `resp_top50` diagnostic pair
    pair_fct_minus_e1_all = evaluate_paired_cal_diff(
        delta_obs=cohort.delta_obs,
        pred_a=preds["fct"],
        pred_b=preds["e1_const_med"],
        t_mean=t_mean,
        target_cols=cohort.target_cols,
        mode="all_nontarget",
        n_boot=n_boot,
        seed=seed,
        half_a_obs=cohort.half_a_obs,
        half_b_obs=cohort.half_b_obs,
        half_valid=cohort.half_valid,
        ortho=True,
    )
    pair_fct_minus_e1_top50 = evaluate_paired_cal_diff(
        delta_obs=cohort.delta_obs,
        pred_a=preds["fct"],
        pred_b=preds["e1_const_med"],
        t_mean=t_mean,
        target_cols=cohort.target_cols,
        mode="resp_top50",
        n_boot=n_boot,
        seed=seed,
        half_a_obs=cohort.half_a_obs,
        half_b_obs=cohort.half_b_obs,
        half_valid=cohort.half_valid,
        ortho=True,
    )
    hd3_lo = float(pair_fct_minus_e1_all.get("diff_ci975_lo", float("nan")))
    hd3_pass = bool(np.isfinite(hd3_lo) and hd3_lo > 0.0)

    # H-D4a & H-D4b (independent Claim 2: resp_top50 trap)
    e1_top50_lo = float(e1_top50.get("r_orth_cal_ci975_lo", float("nan")))
    hd4a_pass = bool(np.isfinite(e1_top50_lo) and e1_top50_lo > 0.0)
    cross_mask_e1 = evaluate_cross_mask_cal_diff(
        delta_obs=cohort.delta_obs,
        pred=preds["e1_const_med"],
        t_mean=t_mean,
        target_cols=cohort.target_cols,
        mode_a="resp_top50",
        mode_b="all_nontarget",
        n_boot=n_boot,
        seed=seed,
        half_a_obs=cohort.half_a_obs,
        half_b_obs=cohort.half_b_obs,
        half_valid=cohort.half_valid,
        ortho=True,
    )
    hd4b_lo = float(cross_mask_e1.get("diff_ci975_lo", float("nan")))
    hd4b_pass = bool(np.isfinite(hd4b_lo) and hd4b_lo > 0.0)

    # H-D_util (parallel 6-metric non-inferiority)
    hd_util = evaluate_hd_util_noninferiority(
        delta_obs=cohort.delta_obs,
        pred_fct=preds["fct"],
        pred_tmean=t_mean,
        grand_mean=grand_mean,
        target_cols=cohort.target_cols,
        n_boot=n_boot,
        seed=seed,
    )

    # SA-8 (`resp_top50_split`) for FCT, E1_ConstMed, and Perm100
    sa8_fct = evaluate_sa8_resp_top50_split(
        delta_obs=cohort.delta_obs,
        pred=preds["fct"],
        t_mean=t_mean,
        target_cols=cohort.target_cols,
        half_a_obs=cohort.half_a_obs,
        half_b_obs=cohort.half_b_obs,
        half_valid=cohort.half_valid,
        n_boot=n_boot,
        seed=seed,
    )
    sa8_e1 = evaluate_sa8_resp_top50_split(
        delta_obs=cohort.delta_obs,
        pred=preds["e1_const_med"],
        t_mean=t_mean,
        target_cols=cohort.target_cols,
        half_a_obs=cohort.half_a_obs,
        half_b_obs=cohort.half_b_obs,
        half_valid=cohort.half_valid,
        n_boot=n_boot,
        seed=seed,
    )
    sa8_perm100 = evaluate_sa8_perm_null(
        delta_obs=cohort.delta_obs,
        t_mean=t_mean,
        fc_bar=fc_bar,
        tgt_ctrl_mean_libnorm=cohort.tgt_ctrl_mean_libnorm,
        target_cols=cohort.target_cols,
        half_a_obs=cohort.half_a_obs,
        half_b_obs=cohort.half_b_obs,
        half_valid=cohort.half_valid,
        ceil_single_ab=float(sa8_fct.get("ceil_single_half_AtoB", 0.0)),
        ceil_single_ba=float(sa8_fct.get("ceil_single_half_BtoA", 0.0)),
        n_perm=n_perm,
    )

    # H-D5 (`all_nontarget` comparison vs `nourreddine_kolf21j`, primary focus on `feng_ipsc`)
    hd5_block: Optional[Dict[str, Any]] = None
    if "fct_mismatch::nourreddine_kolf21j" in preds:
        mis_kolf = mis_all["nourreddine_kolf21j"]
        pair_fct_minus_kolf = evaluate_paired_cal_diff(
            delta_obs=cohort.delta_obs,
            pred_a=preds["fct"],
            pred_b=preds["fct_mismatch::nourreddine_kolf21j"],
            t_mean=t_mean,
            target_cols=cohort.target_cols,
            mode="all_nontarget",
            n_boot=n_boot,
            seed=seed,
            half_a_obs=cohort.half_a_obs,
            half_b_obs=cohort.half_b_obs,
            half_valid=cohort.half_valid,
            ortho=True,
        )
        mis_kolf_hi = float(mis_kolf.get("r_orth_cal_ci975_hi", float("nan")))
        diff_lo = float(pair_fct_minus_kolf.get("diff_ci975_lo", float("nan")))
        hd5_block = {
            "target_line": cohort.target_line,
            "is_feng_ipsc": cohort.target_line == "feng_ipsc",
            "fct_all": {
                "est": fct_all.get("r_orth_cal_est"),
                "ci975_lo": fct_all.get("r_orth_cal_ci975_lo"),
                "ci975_hi": fct_all.get("r_orth_cal_ci975_hi"),
                "valid_frac": fct_all.get("r_orth_cal_valid_frac"),
            },
            "mis_kolf21j_all": {
                "est": mis_kolf.get("r_orth_cal_est"),
                "ci975_lo": mis_kolf.get("r_orth_cal_ci975_lo"),
                "ci975_hi": mis_kolf.get("r_orth_cal_ci975_hi"),
                "valid_frac": mis_kolf.get("r_orth_cal_valid_frac"),
            },
            "paired_fct_minus_mis_kolf21j": pair_fct_minus_kolf,
            "fct_ci975_lo_gt_kolf21j_ci975_hi": bool(
                np.isfinite(fct_all_lo) and np.isfinite(mis_kolf_hi) and fct_all_lo > mis_kolf_hi
            ),
            "paired_diff_ci975_gt_0": bool(np.isfinite(diff_lo) and diff_lo > 0.0),
        }

    # SA-5 top-75% genes on both masks
    tgt_cov_var = np.var(cohort.tgt_ctrl_cells, axis=0, ddof=1)
    q25 = float(np.percentile(tgt_cov_var, 25.0))
    top75_mask = tgt_cov_var > q25
    if int(np.sum(top75_mask)) < 50:
        top75_mask = tgt_cov_var >= q25
    sa5_top75_all = evaluate_joint_residuals(
        delta_obs=cohort.delta_obs,
        pred=preds["fct"],
        t_mean=t_mean,
        target_cols=cohort.target_cols,
        mode="all_nontarget",
        n_boot=n_boot,
        seed=seed,
        half_a_obs=cohort.half_a_obs,
        half_b_obs=cohort.half_b_obs,
        half_valid=cohort.half_valid,
        gene_keep_mask=top75_mask,
    )

    r19_all = (
        evaluate_r19_gain_decomposition(
            delta_obs=cohort.delta_obs,
            t_mean=t_mean,
            target_cols=cohort.target_cols,
            mode="all_nontarget",
            n_boot=n_boot,
            seed=seed,
        )
        if include_r19
        else {}
    )
    r19_top50 = (
        evaluate_r19_gain_decomposition(
            delta_obs=cohort.delta_obs,
            t_mean=t_mean,
            target_cols=cohort.target_cols,
            mode="resp_top50",
            n_boot=n_boot,
            seed=seed,
        )
        if include_r19
        else {}
    )

    s7_all = (
        evaluate_s7_confounder_controls(
            eval_targets=cohort.eval_targets,
            delta_obs=cohort.delta_obs,
            pred_fct=preds["fct"],
            t_mean=t_mean,
            fc_bar=fc_bar,
            tgt_ctrl_mean_libnorm=cohort.tgt_ctrl_mean_libnorm,
            source_names=cohort.source_names,
            source_ctrl_mean_libnorm=cohort.source_ctrl_mean_libnorm,
            target_cols=cohort.target_cols,
            n_cells=cohort.n_cells,
            kd_depth_log1p=cohort.kd_depth_log1p,
            mode="all_nontarget",
            half_a_obs=cohort.half_a_obs,
            half_b_obs=cohort.half_b_obs,
            half_valid=cohort.half_valid,
            n_boot=n_boot,
            n_perm=n_perm,
            seed=seed,
        )
        if include_s7
        else {}
    )
    s7_top50 = (
        evaluate_s7_confounder_controls(
            eval_targets=cohort.eval_targets,
            delta_obs=cohort.delta_obs,
            pred_fct=preds["fct"],
            t_mean=t_mean,
            fc_bar=fc_bar,
            tgt_ctrl_mean_libnorm=cohort.tgt_ctrl_mean_libnorm,
            source_names=cohort.source_names,
            source_ctrl_mean_libnorm=cohort.source_ctrl_mean_libnorm,
            target_cols=cohort.target_cols,
            n_cells=cohort.n_cells,
            kd_depth_log1p=cohort.kd_depth_log1p,
            mode="resp_top50",
            half_a_obs=cohort.half_a_obs,
            half_b_obs=cohort.half_b_obs,
            half_valid=cohort.half_valid,
            n_boot=n_boot,
            n_perm=n_perm,
            seed=seed,
        )
        if include_s7
        else {}
    )

    addendum_notes = []
    if cohort.target_line == "cd4t_stim48":
        addendum_notes.append(
            "Addendum A2: halves CE0008678 and CE0010866 share the 24 sequencing lanes GSM9393875..98 "
            "(channel-NT subtracted, residual shared channel noise elevates ceiling / lowers r_orth_cal conservatively); "
            "PREREG_D Q2 weaker independence: shares CD4+ alpha-beta T lineage, library, & Marson lab with R10 cd4t_rest."
        )
    if cohort.target_line == "feng_ipsc":
        addendum_notes.append(
            "Addendum A3: sorted 17-line parity halves split clonal replicates across the same HipSci donors; "
            "split-half ceiling captures clonal + channel rather than inter-donor variance. Secondary-only (H-D5)."
        )

    return {
        "target_line": cohort.target_line,
        "n_genes_eval": len(cohort.eval_genes),
        "n_eval_targets": len(cohort.eval_targets),
        "n_sources": len(cohort.source_names),
        "source_names": cohort.source_names,
        "const_stats": const_stats,
        "method_eval": method_eval,
        "mismatched_all": mis_all,
        "mismatched_top50": mis_top50,
        "max_mismatch_all_ci975_hi": max_mis_all_ci975_hi,
        "perm100_all_nontarget": perm100_all,
        "perm100_resp_top50": perm100_top50,
        "pair_fct_minus_e1_all": pair_fct_minus_e1_all,
        "pair_fct_minus_e1_top50": pair_fct_minus_e1_top50,
        "cross_mask_e1_top50_minus_all": cross_mask_e1,
        "hd1_pass": hd1_pass,
        "hd1_ci975_gt_0_only": hd1_ci975_gt_0_only,
        "hd1_valid_frac_ok": hd1_valid_frac_ok,
        "hd1_rho_ab_ok": hd1_rho_ab_ok,
        "hd1_a1_deviation_flag": hd1_a1_deviation_flag,
        "hd2_pass": hd2_pass,
        "hd3_pass": hd3_pass,
        "hd4a_pass": hd4a_pass,
        "hd4b_pass": hd4b_pass,
        "hd4_pass": bool(hd4a_pass and hd4b_pass),
        "hd_util": hd_util,
        "hd_util_pass": bool(hd_util["hd_util_pass_line"]),
        "sa8_fct": sa8_fct,
        "sa8_e1": sa8_e1,
        "sa8_perm100": sa8_perm100,
        "hd5_kolf21j": hd5_block,
        "sa5_top75_all_fct": sa5_top75_all,
        "r19_all": r19_all,
        "r19_top50": r19_top50,
        "s7_all": s7_all,
        "s7_top50": s7_top50,
        "addendum_notes": addendum_notes,
    }


def verify_against_r10e_csvs(
    per_line_results: Dict[str, Dict[str, Any]],
    r10e_dir: Path,
    tol: float = 1e-10,
) -> Dict[str, Any]:
    """
    Numerically verify `round11/audit` outputs (run at `seed=20261122`, `n_boot=4000`, `n_perm=100`,
    `source_mode=r9_6`) against frozen `round10/exploratory/baseline_ablation_*.csv` reference values.
    """
    m_csv = r10e_dir / "baseline_ablation_metrics.csv"
    p_csv = r10e_dir / "baseline_ablation_paired_fct_vs_e1.csv"
    if not m_csv.exists() or not p_csv.exists():
        raise FileNotFoundError(f"Missing R10-E reference CSVs in {r10e_dir}")

    method_map = {
        "FCT": "fct",
        "E1_ConstMed": "e1_const_med",
        "E2_ConstP10": "e2_const_p10",
        "E2_ConstP90": "e2_const_p90",
        "E3_SrcMeanCtrl": "e3_src_mean_ctrl",
        "E5_RankMatchedMu": "e5_rank_matched_mu",
        "E6_LogFCBar": "e6_log_fc_bar",
        "E7_FCBarMinus1": "e7_fc_minus_1",
    }

    checks: List[Dict[str, Any]] = []
    max_diff = 0.0

    with m_csv.open("r", newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            ln = row["target_line"]
            mask = row["mask"]
            meth = row["method"]
            if ln not in per_line_results:
                continue
            info = per_line_results[ln]
            if meth == "Perm100_Median":
                p_blk = info["perm100_all_nontarget"] if mask == "all_nontarget" else info["perm100_resp_top50"]
                val = float(p_blk["r_orth_cal"]["median"])
                ref = float(row["r_orth_cal"])
                d = abs(val - ref)
                max_diff = max(max_diff, d)
                checks.append({"line": ln, "mask": mask, "metric": "Perm100_Median", "audit": val, "ref": ref, "abs_diff": d})
                val_p975 = float(p_blk["r_orth_cal"]["p975"])
                ref_p975 = float(row["r_orth_cal_ci975_hi"])
                d2 = abs(val_p975 - ref_p975)
                max_diff = max(max_diff, d2)
                checks.append({"line": ln, "mask": mask, "metric": "Perm100_p975", "audit": val_p975, "ref": ref_p975, "abs_diff": d2})
                continue
            if meth not in method_map:
                continue
            m_key = method_map[meth]
            res = info["method_eval"][m_key][mask]
            for csv_k, res_k in (
                ("r_orth_cal", "r_orth_cal_est"),
                ("r_orth_cal_ci975_lo", "r_orth_cal_ci975_lo"),
                ("r_orth_cal_ci975_hi", "r_orth_cal_ci975_hi"),
                ("r_orth_raw", "r_orth_raw_est"),
                ("r_cal", "r_cal_est"),
                ("slope_ols", "slope_ols_est"),
                ("rho_half_ab_orth", "rho_ab_orth"),
                ("ceiling_r_orth", "noise_ceiling_r_orth"),
            ):
                val = float(res[res_k])
                ref = float(row[csv_k])
                d = abs(val - ref)
                max_diff = max(max_diff, d)
                checks.append({"line": ln, "mask": mask, "metric": f"{meth}.{csv_k}", "audit": val, "ref": ref, "abs_diff": d})

    with p_csv.open("r", newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            ln = row["target_line"]
            mask = row["mask"]
            if ln not in per_line_results:
                continue
            info = per_line_results[ln]
            p_blk = info["pair_fct_minus_e1_all"] if mask == "all_nontarget" else info["pair_fct_minus_e1_top50"]
            for csv_k, p_k in (
                ("diff_fct_minus_e1", "diff_est"),
                ("diff_ci975_lo", "diff_ci975_lo"),
                ("diff_ci975_hi", "diff_ci975_hi"),
            ):
                val = float(p_blk[p_k])
                ref = float(row[csv_k])
                d = abs(val - ref)
                max_diff = max(max_diff, d)
                checks.append({"line": ln, "mask": mask, "metric": f"FCT-E1.{csv_k}", "audit": val, "ref": ref, "abs_diff": d})

    return {
        "n_checks": len(checks),
        "max_abs_diff": max_diff,
        "passes_tol": bool(len(checks) > 0 and max_diff <= tol),
        "tol": tol,
        "sample_mismatches": [c for c in checks if c["abs_diff"] > tol][:10],
    }


def _write_summary_csvs(per_line: Dict[str, Dict[str, Any]], out_dir: Path) -> None:
    rows: List[Dict[str, Any]] = []
    for ln, info in per_line.items():
        for m_name, masks_d in info["method_eval"].items():
            for mask_name, blk in masks_d.items():
                if not blk:
                    continue
                rows.append(
                    {
                        "target_line": ln,
                        "method": m_name,
                        "mask": mask_name,
                        "n_eval_targets": info["n_eval_targets"],
                        "n_genes_eval": info["n_genes_eval"],
                        "r_orth_cal": blk.get("r_orth_cal_est"),
                        "r_orth_cal_ci975_lo": blk.get("r_orth_cal_ci975_lo"),
                        "r_orth_cal_ci975_hi": blk.get("r_orth_cal_ci975_hi"),
                        "r_orth_cal_valid_frac": blk.get("r_orth_cal_valid_frac"),
                        "r_orth_raw": blk.get("r_orth_raw_est"),
                        "r_cal": blk.get("r_cal_est"),
                        "slope_ols": blk.get("slope_ols_est"),
                        "rho_ab_orth": blk.get("rho_ab_orth"),
                        "noise_ceiling_r_orth": blk.get("noise_ceiling_r_orth"),
                    }
                )
    if rows:
        with (out_dir / "audit_r11_per_line.csv").open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)


def render_audit_markdown(
    summary: Dict[str, Any],
) -> str:
    adj = summary["adjudication"]
    lines: List[str] = [
        f"# Round 11 (`PREREG_D` + Addendum A) Independent Statistical Audit (`{summary['mode']}`)",
        "",
        f"- **Parameters:** `n_boot = {summary['n_boot']}`, `n_perm = {summary['n_perm']}`, `seed = {summary['seed']}`, Bonferroni two-sided `97.5%` CIs (`[1.25%, 98.75%]`), `--source-mode {summary['source_mode']}`.",
        f"- **Adjudication Status:** `{adj['status']}` (`K_conf = {adj['k_conf']}`)",
    ]
    if adj["status"] == "NO_CONFIRMATORY_ADJUDICATION":
        lines.append(f"- **P1 Notice:** `{adj['abstract_notice']}`")
    else:
        lines.extend(
            [
                f"- **Claim 1 (Primary Biological `H-D1 ∩ H-D2 ∩ H-D3` under `all_nontarget`):** **{adj['primary_verdict']}** (`H-D1`: `{adj['hd1_pass_count']}`, `H-D2`: `{adj['hd2_pass_count']}`, `H-D3`: `{adj['hd3_pass_count']}`)",
                f"- **Claim 2 (Methodological `H-D4 = H-D4a ∩ H-D4b`, independent):** **{adj['hd4_verdict']}** (`H-D4a`: `{adj['hd4a_pass_count']}`, `H-D4b`: `{adj['hd4b_pass_count']}`)",
                f"- **Parallel Utility (`H-D_util`):** **{adj['hd_util_verdict']}** (`{adj['hd_util_pass_count']}`)",
            ]
        )
        if adj.get("hd1_a1_deviation_lines"):
            lines.append(
                f"- **Addendum A1 Deviation Flagged (`CI975_lo > 0` while `valid_frac < 0.95` or `rho_ab_orth <= 0`):** `{adj['hd1_a1_deviation_lines']}`"
            )

    lines.extend(
        [
            "",
            "## Table 1: Primary Whole-Transcriptome Endpoints (`all_nontarget`, `97.5% CI`)",
            "",
            "| Line (`N`) | `FCT` `r_orth_cal [97.5% CI]` | `valid_frac` / `rho_AB` | `E1_ConstMed` | `FCT - E1 [97.5% CI]` | `max Mis97.5_hi` | `Perm100 Q975` | `H-D1 / D2 / D3` |",
            "| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |",
        ]
    )
    for ln, info in summary["per_line"].items():
        fa = info["method_eval"]["fct"]["all_nontarget"]
        ea = info["method_eval"]["e1_const_med"]["all_nontarget"]
        p3 = info["pair_fct_minus_e1_all"]
        pm = info["perm100_all_nontarget"]["r_orth_cal"]["p975"]
        lines.append(
            f"| `{ln}` (`N={info['n_eval_targets']}`) "
            f"| `{_fmt(fa['r_orth_cal_est'])} [{_fmt(fa['r_orth_cal_ci975_lo'])}, {_fmt(fa['r_orth_cal_ci975_hi'])}]` "
            f"| `{float(fa['r_orth_cal_valid_frac']):.3f}` / `{_fmt(fa['rho_ab_orth'])}` "
            f"| `{_fmt(ea['r_orth_cal_est'])} [{_fmt(ea['r_orth_cal_ci975_lo'])}, {_fmt(ea['r_orth_cal_ci975_hi'])}]` "
            f"| `{_fmt(p3['diff_est'])} [{_fmt(p3['diff_ci975_lo'])}, {_fmt(p3['diff_ci975_hi'])}]` "
            f"| `{_fmt(info['max_mismatch_all_ci975_hi'])}` | `{_fmt(pm)}` "
            f"| `{info['hd1_pass']} / {info['hd2_pass']} / {info['hd3_pass']}` |"
        )

    lines.extend(
        [
            "",
            "## Table 2: Selection Trap (`H-D4` & `SA-8 resp_top50_split`) and Dual-Mask Permutations",
            "",
            "| Line | `E1(resp_top50) [97.5% CI]` | `E1(top50 - all) [97.5% CI]` | `SA-8 FCT split` | `SA-8 E1 split` | `SA-8 Perm100 Med / Q975` | `Perm100 top50 Q975` | `H-D4a / D4b` |",
            "| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |",
        ]
    )
    for ln, info in summary["per_line"].items():
        et = info["method_eval"]["e1_const_med"]["resp_top50"]
        cm = info["cross_mask_e1_top50_minus_all"]
        sf = info["sa8_fct"]
        se = info["sa8_e1"]
        sp = info["sa8_perm100"]["r_orth_cal_split"]
        pt = info["perm100_resp_top50"]["r_orth_cal"]["p975"]
        lines.append(
            f"| `{ln}` "
            f"| `{_fmt(et['r_orth_cal_est'])} [{_fmt(et['r_orth_cal_ci975_lo'])}, {_fmt(et['r_orth_cal_ci975_hi'])}]` "
            f"| `{_fmt(cm['diff_est'])} [{_fmt(cm['diff_ci975_lo'])}, {_fmt(cm['diff_ci975_hi'])}]` "
            f"| `{_fmt(sf.get('r_orth_cal_split_est'))} [{_fmt(sf.get('r_orth_cal_split_ci975_lo'))}, {_fmt(sf.get('r_orth_cal_split_ci975_hi'))}]` "
            f"| `{_fmt(se.get('r_orth_cal_split_est'))} [{_fmt(se.get('r_orth_cal_split_ci975_lo'))}, {_fmt(se.get('r_orth_cal_split_ci975_hi'))}]` "
            f"| `{_fmt(sp['median'])} / {_fmt(sp['p975'])}` | `{_fmt(pt)}` "
            f"| `{info['hd4a_pass']} / {info['hd4b_pass']}` |"
        )

    lines.extend(
        [
            "",
            "## Table 3: `H-D5` (`all_nontarget` FCT vs `Mis(nourreddine_kolf21j)`) & Addendum Notes",
            "",
            "| Line | `FCT` (`all_nontarget`) | `Mis(kolf21j)` (`all_nontarget`) | `FCT - Mis(kolf21j) [97.5% CI]` | `Diff CI > 0` |",
            "| :--- | :---: | :---: | :---: | :---: |",
        ]
    )
    for ln, info in summary["per_line"].items():
        h5 = info.get("hd5_kolf21j")
        if not h5:
            continue
        f5 = h5["fct_all"]
        m5 = h5["mis_kolf21j_all"]
        p5 = h5["paired_fct_minus_mis_kolf21j"]
        lines.append(
            f"| `{ln}` "
            f"| `{_fmt(f5['est'])} [{_fmt(f5['ci975_lo'])}, {_fmt(f5['ci975_hi'])}]` "
            f"| `{_fmt(m5['est'])} [{_fmt(m5['ci975_lo'])}, {_fmt(m5['ci975_hi'])}]` "
            f"| `{_fmt(p5['diff_est'])} [{_fmt(p5['diff_ci975_lo'])}, {_fmt(p5['diff_ci975_hi'])}]` "
            f"| `{h5['paired_diff_ci975_gt_0']}` |"
        )
        for note in info.get("addendum_notes", []):
            lines.append(f"- **Disclosure (`{ln}`):** {note}")
    return "\n".join(lines) + "\n"


def _resolve_line_dir(line_name: str, mode: str) -> Path:
    if mode == "dev_r9r10":
        if line_name in R9_LINE_SET:
            return R9_LINES_DIR / line_name
        return R10_LINES_DIR / line_name
    if line_name in R9_LINE_SET:
        return R9_LINES_DIR / line_name
    if line_name in ("dld1", "cd4t_rest"):
        return R10_LINES_DIR / line_name
    return R11_LINES_DIR / line_name


def main() -> None:
    parser = argparse.ArgumentParser(description="Independent Round 11 Group D Audit Runner (A-T3)")
    parser.add_argument("--mode", choices=("dev_r9r10", "holdout_r11"), default="dev_r9r10")
    parser.add_argument("--source-mode", choices=("r9_6", "r9_r10_8"), default="r9_6")
    parser.add_argument("--target-lines", nargs="+", default=None)
    parser.add_argument("--n-boot", type=int, default=4000)
    parser.add_argument("--n-perm", type=int, default=100)
    parser.add_argument("--seed", type=int, default=DEFAULT_R11_SEED)
    parser.add_argument("--cov-max-cells", type=int, default=None)
    parser.add_argument("--no-sa3", action="store_true")
    parser.add_argument("--no-r19", action="store_true")
    parser.add_argument("--no-s7", action="store_true")
    parser.add_argument("--no-hgnc", action="store_true")
    parser.add_argument("--verify-r10e-dir", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()

    dev_only = args.mode == "dev_r9r10"
    if dev_only:
        targets = list(args.target_lines or ("dld1", "cd4t_rest"))
        enforce_r11_gate(targets, dev_only=True, extra_lines_dirs=[_resolve_line_dir(t, "dev_r9r10") for t in targets])
        confirmatory_lines: List[str] = []
    else:
        enforce_r11_gate(args.target_lines or [], dev_only=False)
        conf_list = parse_confirmatory_lines_r11()
        confirmatory_lines = list(conf_list)
        if args.target_lines is not None:
            targets = list(args.target_lines)
        else:
            targets = list(conf_list)
            for cand in R11_CANDIDATE_CONFIRMATORY + R11_SECONDARY_ONLY:
                if (R11_LINES_DIR / cand / "pseudobulk_delta.npz").exists() and cand not in targets:
                    targets.append(cand)

    hgnc = None if args.no_hgnc else CanonicalHGNC.from_tsv()
    common_genes = load_r9_common_genes(hgnc=hgnc)

    source_pool_names = list(R9_LINE_NAMES)
    if args.source_mode == "r9_r10_8":
        source_pool_names.extend(R10_SOURCE_EXT_NAMES)

    loaded_cache: Dict[str, AuditLineData] = {}
    for sn in source_pool_names:
        sdir = R9_LINES_DIR / sn if sn in R9_LINE_SET else R10_LINES_DIR / sn
        loaded_cache[sn] = load_line_dataset(sdir, representation="log1p", aligned_genes=common_genes, hgnc=hgnc, min_gene_overlap=3000)

    # Precompute sample covariances for shared sources
    cov_cache: Dict[str, np.ndarray] = {}
    gene_order_cache: Dict[str, Sequence[str]] = {}
    if not args.no_sa3:
        for sn, sdata in loaded_cache.items():
            cov_cache[sn] = sample_covariance(deterministic_subsample(sdata.ctrl_cells, args.cov_max_cells, seed=42), shrink=0.0)
            gene_order_cache[sn] = sdata.genes

    per_line: Dict[str, Dict[str, Any]] = {}
    for tgt in targets:
        if tgt in loaded_cache:
            tdata = loaded_cache[tgt]
        else:
            min_ov = 100 if tgt in R11_SECONDARY_ONLY else 3000
            tdata = load_line_dataset(
                _resolve_line_dir(tgt, args.mode),
                representation="log1p",
                aligned_genes=common_genes,
                hgnc=hgnc,
                min_gene_overlap=min_ov,
            )
            if not args.no_sa3:
                cov_cache[tgt] = sample_covariance(deterministic_subsample(tdata.ctrl_cells, args.cov_max_cells, seed=42), shrink=0.0)
                gene_order_cache[tgt] = tdata.genes

        # Zero self-leakage: exclude `tgt` if it appears in `source_pool_names`
        src_list = [loaded_cache[sn] for sn in source_pool_names if sn != tgt]
        per_line[tgt] = evaluate_target_line_r11(
            target_data=tdata,
            source_list=src_list,
            only_reproducible=True,
            cov_max_cells=args.cov_max_cells,
            n_boot=args.n_boot,
            n_perm=args.n_perm,
            seed=args.seed,
            include_sa3_anchors=not args.no_sa3,
            include_exploratory_e=True,
            include_r19=not args.no_r19,
            include_s7=not args.no_s7,
            precomputed_covs=cov_cache if not args.no_sa3 else None,
            precomputed_gene_orders=gene_order_cache if not args.no_sa3 else None,
        )

    adjudication = adjudicate_round11_group_d(
        per_line_summary=per_line,
        confirmatory_lines=confirmatory_lines if args.mode == "holdout_r11" else [],
    )
    r10e_check = (
        verify_against_r10e_csvs(per_line, args.verify_r10e_dir)
        if args.verify_r10e_dir is not None
        else None
    )

    summary = {
        "mode": args.mode,
        "source_mode": args.source_mode,
        "n_boot": args.n_boot,
        "n_perm": args.n_perm,
        "seed": args.seed,
        "targets": targets,
        "confirmatory_lines": confirmatory_lines,
        "adjudication": adjudication,
        "r10e_verification": r10e_check,
        "per_line": per_line,
    }

    args.out_dir.mkdir(parents=True, exist_ok=True)
    _write_summary_csvs(per_line, args.out_dir)
    (args.out_dir / "audit_r11_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    (args.out_dir / "AUDIT_R11_REPORT.md").write_text(
        render_audit_markdown(summary),
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "status": f"R11_AUDIT_{args.mode.upper()}_COMPLETE",
                "targets": targets,
                "adjudication_status": adjudication["status"],
                "r10e_verification": r10e_check,
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
