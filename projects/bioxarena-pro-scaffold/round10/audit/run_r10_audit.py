#!/usr/bin/env python3
"""
Independent Round 10 Group C Audit Runner & Dual-Team C3 Cross-Checker (A-T3, PREREG_C).

Usage modes:
  1. `--mode r9_dev` (Pre-freeze safe testing strictly on Round 9's 6 lines in `round9/p2/shared/lines/`):
       Runs independent FCT + mismatched controls + 100-perm + SA-1..SA-6 / S7 suite in LOLO on R9 lines,
       refusing any access to Round 10 directories or lines.
  2. `--mode holdout_r10` (Post-freeze Round 10 evaluation):
       Strictly enforces `PREREG_C_COMMIT.txt` + `CODE_FREEZE_R10.txt` gate, reads Round 9 sources +
       `round10/shared/lines/<line>` + `CONFIRMATORY_LINES.txt`, enforces P1 (`K_conf < 2` -> descriptive only),
       runs full primary (`FCT` + 6 source mismatched controls + 100-perm), H-C3, and secondary SA-1..SA-6 / S7.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

from independent_r10_loader import (
    R9_LINE_NAMES,
    R9_LINES_DIR,
    R10_LINES_DIR,
    R10_SECONDARY_ONLY,
    R10_STATUS_DIR,
    CanonicalHGNC,
    enforce_r10_gate,
    load_line_dataset,
    load_r9_common_genes,
    parse_confirmatory_lines,
)
from independent_r10_metrics import (
    adjudicate_round10_group_c,
    evaluate_fct_perm_null,
    evaluate_hc3_noninferiority,
    evaluate_joint_residuals,
    evaluate_paired_cal_diff,
    evaluate_r19_gain_decomposition,
    evaluate_s7_confounder_controls,
)
from independent_r10_predictors import (
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


def evaluate_target_line(
    target_data: Any,
    source_list: Sequence[Any],
    only_reproducible: bool = True,
    shrink: float = 0.0,
    cov_max_cells: Optional[int] = None,
    n_boot: int = 4000,
    n_perm: int = 100,
    seed: int = 20261122,
    include_sa3_anchors: bool = True,
    include_r19: bool = True,
    precomputed_covs: Optional[Dict[str, np.ndarray]] = None,
    precomputed_gene_orders: Optional[Dict[str, Sequence[str]]] = None,
) -> Dict[str, Any]:
    """Run complete independent audit on one target line against its source cohort."""
    cohort = align_source_and_target(source_list, target_data, only_reproducible=only_reproducible)
    preds, fc_bar, grand_mean = build_predictions(
        cohort=cohort,
        shrink=shrink,
        cov_max_cells=cov_max_cells,
        seed=42,
        include_sa3_anchors=include_sa3_anchors,
        precomputed_covs=precomputed_covs,
        precomputed_gene_orders=precomputed_gene_orders,
    )
    t_mean = preds["transferred_mean"]

    # Primary FCT + comparators (resp_top50 on all methods; all_nontarget on primary + SA-2/SA-3/SA-5 anchors)
    all_nt_methods = {"fct", "fct_full_nt", "transfer_mean_global_scale", "hybrid", "transfer_mean_var_scaled", "hybrid_on_tmvs"}
    method_eval: Dict[str, Dict[str, Any]] = {}
    for m_name, p_mat in preds.items():
        m_top50 = evaluate_joint_residuals(
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
        m_all = (
            evaluate_joint_residuals(
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
            if m_name in all_nt_methods
            else {}
        )
        method_eval[m_name] = {"resp_top50": m_top50, "all_nontarget": m_all}

    fct50 = method_eval["fct"]["resp_top50"]
    ceil_raw = float(fct50.get("noise_ceiling_r", 0.0))
    ceil_orth = float(fct50.get("noise_ceiling_r_orth", 0.0))

    # 100-Perm Null (`NegCtrl_FCT_Perm(0..99)`)
    perm100 = evaluate_fct_perm_null(
        delta_obs=cohort.delta_obs,
        t_mean=t_mean,
        fc_bar=fc_bar,
        tgt_ctrl_mean_libnorm=cohort.tgt_ctrl_mean_libnorm,
        target_cols=cohort.target_cols,
        ceil_raw=ceil_raw,
        ceil_orth=ceil_orth,
        mode="resp_top50",
        n_perm=n_perm,
    )

    # Mismatched source control summary
    mis_dict = {
        s_name: method_eval[f"fct_mismatch::{s_name}"]["resp_top50"]
        for s_name in cohort.source_names
        if f"fct_mismatch::{s_name}" in method_eval
    }
    mis_ci975_his = [v.get("r_orth_cal_ci975_hi", float("nan")) for v in mis_dict.values()]
    max_mis_ci975_hi = float(np.nanmax(mis_ci975_his)) if mis_ci975_his else float("nan")
    perm_p975 = float(perm100["r_orth_cal"]["p975"])
    fct_lo = float(fct50.get("r_orth_cal_ci975_lo", float("nan")))

    hc1_pass = bool(np.isfinite(fct_lo) and fct_lo > 0.0)
    hc2_pass = bool(
        np.isfinite(fct_lo)
        and np.isfinite(max_mis_ci975_hi)
        and np.isfinite(perm_p975)
        and (fct_lo > max(max_mis_ci975_hi, perm_p975))
    )

    # H-C3 non-inferiority (6 paired conventions + pds_cosine_panel)
    hc3 = evaluate_hc3_noninferiority(
        delta_obs=cohort.delta_obs,
        pred_fct=preds["fct"],
        pred_tmean=t_mean,
        grand_mean=grand_mean,
        target_cols=cohort.target_cols,
        n_boot=n_boot,
        seed=seed,
    )

    # SA-5 Top-75% expressed/variable genes
    tgt_cov_var = np.var(cohort.tgt_ctrl_cells, axis=0, ddof=1)
    q25 = float(np.percentile(tgt_cov_var, 25.0))
    top75_mask = tgt_cov_var > q25
    if int(np.sum(top75_mask)) < 50:
        top75_mask = tgt_cov_var >= q25
    sa5_top75_fct = evaluate_joint_residuals(
        delta_obs=cohort.delta_obs,
        pred=preds["fct"],
        t_mean=t_mean,
        target_cols=cohort.target_cols,
        mode="resp_top50",
        n_boot=n_boot,
        seed=seed,
        half_a_obs=cohort.half_a_obs,
        half_b_obs=cohort.half_b_obs,
        half_valid=cohort.half_valid,
        gene_keep_mask=top75_mask,
    )

    # R14 paired check if SA-3 anchors computed
    r14_ci = None
    if include_sa3_anchors and "hybrid_on_tmvs" in preds and "transfer_mean_var_scaled" in preds:
        r14_ci = evaluate_paired_cal_diff(
            delta_obs=cohort.delta_obs,
            pred_a=preds["hybrid_on_tmvs"],
            pred_b=preds["transfer_mean_var_scaled"],
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

    # SA-4 R19 3-Tier gain decomposition
    r19 = {}
    if include_r19:
        r19["resp_top50"] = evaluate_r19_gain_decomposition(
            delta_obs=cohort.delta_obs,
            t_mean=t_mean,
            target_cols=cohort.target_cols,
            mode="resp_top50",
            n_boot=n_boot,
            seed=seed,
        )
        r19["all_nontarget"] = evaluate_r19_gain_decomposition(
            delta_obs=cohort.delta_obs,
            t_mean=t_mean,
            target_cols=cohort.target_cols,
            mode="all_nontarget",
            n_boot=n_boot,
            seed=seed,
        )

    # SA-6 / S7 Confounder Controls
    s7 = evaluate_s7_confounder_controls(
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
        half_a_obs=cohort.half_a_obs,
        half_b_obs=cohort.half_b_obs,
        half_valid=cohort.half_valid,
        n_boot=n_boot,
        n_perm=n_perm,
        seed=seed,
    )

    return {
        "line": target_data.name,
        "n_eval_genes": len(cohort.eval_genes),
        "n_eval_targets": len(cohort.eval_targets),
        "only_reproducible": only_reproducible,
        "hc1_pass": hc1_pass,
        "hc2_pass": hc2_pass,
        "hc3_pass": hc3["hc3_pass_line"],
        "max_mismatch_r_orth_cal_ci975_hi": max_mis_ci975_hi,
        "perm100": perm100,
        "methods": method_eval,
        "hc3": hc3,
        "sa5_top75_fct": sa5_top75_fct,
        "r14_hotmvs_minus_tmvs": r14_ci,
        "r19": r19,
        "s7": s7,
    }


def run_cross_check(
    line_results: Dict[str, Dict[str, Any]],
    audit_rows: List[Dict[str, Any]],
    results_csv: Path,
    results_summary_json: Optional[Path] = None,
) -> Dict[str, Any]:
    """Full C3 parity cross-check comparing independent A-T3 audit against A-T2 `round10/results/`."""
    import pandas as pd

    if not results_csv.exists():
        return {"status": "MISSING_RESULTS_CSV", "path": str(results_csv)}
    df_t2 = pd.read_csv(results_csv)

    # Auto-detect sibling summary_*.json if omitted
    if results_summary_json is None and results_csv.name.startswith("line_metrics_"):
        cand = results_csv.with_name(results_csv.name.replace("line_metrics_", "summary_", 1)).with_suffix(".json")
        if cand.exists():
            results_summary_json = cand
    sum_t2 = json.loads(results_summary_json.read_text(encoding="utf-8")) if (results_summary_json and results_summary_json.exists()) else None

    checks: List[Dict[str, Any]] = []
    max_pt_diff = 0.0
    max_ci_diff = 0.0

    def _record(ln: str, key: str, a3_val: Any, a2_val: Any, is_ci: bool = False) -> None:
        nonlocal max_pt_diff, max_ci_diff
        try:
            v3, v2 = float(a3_val), float(a2_val)
        except (TypeError, ValueError):
            return
        if not (np.isfinite(v3) and np.isfinite(v2)):
            return
        diff = abs(v3 - v2)
        if is_ci:
            max_ci_diff = max(max_ci_diff, diff)
        else:
            max_pt_diff = max(max_pt_diff, diff)
        checks.append({"line": ln, "metric": key, "a3": v3, "a2": v2, "abs_diff": diff, "is_ci": is_ci})

    # Detect long format (`target_line`, `method`) from round10/src/
    if "target_line" in df_t2.columns and "method" in df_t2.columns:
        method_map = {
            "Source_Grand_Mean": "source_grand_mean",
            "FCT": "fct",
            "Transfer_Mean_GlobalScale": "transfer_mean_global_scale",
            "Transfer_Mean_VarScaled": "transfer_mean_var_scaled",
            "CIPHER_Linear": "cipher_linear",
            "Hybrid_ProjShared": "hybrid",
            "Hybrid_on_TMVS": "hybrid_on_tmvs",
        }
        for ln, res in line_results.items():
            sub_ln = df_t2[df_t2["target_line"] == ln]
            if sub_ln.empty:
                continue
            # Also add all mismatched mean methods dynamically
            line_method_map = dict(method_map)
            for m_key in res["methods"]:
                if m_key.startswith("fct_mismatch::"):
                    sname = m_key.split("::", 1)[1]
                    line_method_map[f"NegCtrl_FCT_MismatchedMean({sname})"] = m_key

            for a2_meth, a3_meth in line_method_map.items():
                m_row = sub_ln[sub_ln["method"] == a2_meth]
                if len(m_row) != 1 or a3_meth not in res["methods"]:
                    continue
                r2 = m_row.iloc[0]
                m3 = res["methods"][a3_meth]["resp_top50"]
                if a2_meth != "Transfer_Mean_GlobalScale":
                    _record(ln, f"{a2_meth}.r_orth_cal", m3.get("r_orth_cal_est"), r2.get("r_orth_cal"), is_ci=False)
                    _record(ln, f"{a2_meth}.r_orth_cal_ci975_lo", m3.get("r_orth_cal_ci975_lo"), r2.get("r_orth_cal_ci975_lo"), is_ci=True)
                    _record(ln, f"{a2_meth}.r_orth_cal_ci975_hi", m3.get("r_orth_cal_ci975_hi"), r2.get("r_orth_cal_ci975_hi"), is_ci=True)
                _record(ln, f"{a2_meth}.r_cal", m3.get("r_cal_est"), r2.get("r_cal"), is_ci=False)
                _record(ln, f"{a2_meth}.slope_ols", m3.get("slope_ols_est"), r2.get("slope_ols"), is_ci=False)
                m3_all = res["methods"][a3_meth].get("all_nontarget", {})
                if m3_all and "all_nontarget_r_orth_cal" in r2 and a2_meth != "Transfer_Mean_GlobalScale":
                    _record(ln, f"{a2_meth}.all_r_orth_cal", m3_all.get("r_orth_cal_est"), r2.get("all_nontarget_r_orth_cal"), is_ci=False)

            # H-C3 aggregate panel/ownmasked/de50 parity on FCT & Transfer_Mean
            hc3 = res.get("hc3", {})
            for a2_meth, hc3_field in (("FCT", "fct_mean"), ("Transfer_Mean", "tmean_mean")):
                m_row = sub_ln[sub_ln["method"] == a2_meth]
                if len(m_row) != 1:
                    continue
                r2 = m_row.iloc[0]
                for c3_k, c2_k in (
                    ("ctrl_disc_panel", "ctrl_disc_raw_panel"),
                    ("ctrl_disc_ownmasked", "ctrl_disc_raw"),
                    ("ctrl_de50_r", "ctrl_de50_r"),
                    ("cent_disc_panel", "cent_disc_raw_panel"),
                    ("cent_disc_ownmasked", "cent_disc_raw"),
                    ("cent_de50_r", "cent_de50_r"),
                    ("pds_cosine_panel", "pds_cosine_panel"),
                ):
                    if c3_k in hc3 and c2_k in r2:
                        _record(ln, f"{a2_meth}.hc3.{c3_k}", hc3[c3_k].get(hc3_field), r2.get(c2_k), is_ci=False)

            if sum_t2 is not None:
                p2 = sum_t2.get("fct_perm100", {}).get(ln, {}).get("r_orth_cal", {})
                p3 = res["perm100"]["r_orth_cal"]
                _record(ln, "perm100.p975", p3.get("p975"), p2.get("p97_5"), is_ci=False)
                _record(ln, "perm100.median", p3.get("median"), p2.get("median"), is_ci=False)

                s2 = sum_t2.get("s7_cellcount_kd", {}).get(ln, {})
                s3 = res["s7"]
                if s2 and s3:
                    _record(ln, "s7.sp_ln_ncells", s3["spearman_ln_ncells"].get("rho_est"), s2.get("spearman_r_orth_p_vs_ln_ncells", {}).get("estimate"), is_ci=False)
                    _record(ln, "s7.sp_abs_kp", s3["spearman_abs_kp"].get("rho_est"), s2.get("spearman_r_orth_p_vs_abs_kd", {}).get("estimate"), is_ci=False)
                    _record(ln, "s7.ols_intercept", s3.get("ols_centered", {}).get("intercept_adjusted_est"), s2.get("ols_adjusted_both", {}).get("adjusted_intercept", {}).get("estimate"), is_ci=False)
                    for t3_k, t2_k in (("T1_low_ncells", "T1_low"), ("T2_mid_ncells", "T2_mid"), ("T3_high_ncells", "T3_high")):
                        t3_v = s3.get("tertiles", {}).get(t3_k, {})
                        t2_v = s2.get("n_cells_tertiles", {}).get(t2_k, {})
                        if t3_v and t2_v:
                            _record(ln, f"s7.{t2_k}.r_orth_cal", t3_v.get("fct_r_orth_cal"), t2_v.get("fct_r_orth_cal", {}).get("estimate"), is_ci=False)
                            _record(ln, f"s7.{t2_k}.ci975_lo", t3_v.get("fct_r_orth_cal_ci975_lo"), t2_v.get("fct_r_orth_cal", {}).get("ci975_lo"), is_ci=True)

    return {
        "n_checks": len(checks),
        "max_point_abs_diff": max_pt_diff,
        "max_ci_abs_diff": max_ci_diff,
        "passes_c3_tolerance": bool(len(checks) > 0 and max_pt_diff <= 1e-4 and max_ci_diff <= 0.01),
        "top_diffs": sorted(checks, key=lambda x: -x["abs_diff"])[:10],
    }


def _flatten_csv_row(res: Dict[str, Any]) -> Dict[str, Any]:
    f50 = res["methods"]["fct"]["resp_top50"]
    fall = res["methods"]["fct"]["all_nontarget"]
    gs = res["methods"].get("transfer_mean_global_scale", {}).get("resp_top50", {})
    hyb = res["methods"].get("hybrid", {}).get("resp_top50", {})
    tmvs = res["methods"].get("transfer_mean_var_scaled", {}).get("resp_top50", {})
    hc3 = res["hc3"]
    s7 = res["s7"]
    ols = s7.get("ols_centered", {})
    tert = s7.get("tertiles", {})

    row = {
        "line": res["line"],
        "n_eval_genes": res["n_eval_genes"],
        "n_eval_targets": res["n_eval_targets"],
        "fct_r_orth_cal": f50.get("r_orth_cal_est"),
        "fct_r_orth_cal_ci975_lo": f50.get("r_orth_cal_ci975_lo"),
        "fct_r_orth_cal_ci975_hi": f50.get("r_orth_cal_ci975_hi"),
        "fct_r_cal": f50.get("r_cal_est"),
        "fct_slope_ols": f50.get("slope_ols_est"),
        "fct_slope_ols_ci975_lo": f50.get("slope_ols_ci975_lo"),
        "fct_slope_ols_ci975_hi": f50.get("slope_ols_ci975_hi"),
        "fct_sd_ratio": f50.get("sd_pred_over_sd_obs_est"),
        "fct_slope_deming": f50.get("slope_deming_est"),
        "rho_ab_orth_top50": f50.get("rho_ab_orth"),
        "ceil_orth_top50": f50.get("noise_ceiling_r_orth"),
        "fct_all_r_orth_cal": fall.get("r_orth_cal_est"),
        "rho_ab_orth_all": fall.get("rho_ab_orth"),
        "fct_top75_r_orth_cal": res["sa5_top75_fct"].get("r_orth_cal_est"),
        "max_mismatch_r_orth_cal_ci975_hi": res["max_mismatch_r_orth_cal_ci975_hi"],
        "perm100_median": res["perm100"]["r_orth_cal"]["median"],
        "perm100_p975": res["perm100"]["r_orth_cal"]["p975"],
        "hc1_pass": res["hc1_pass"],
        "hc2_pass": res["hc2_pass"],
        "hc3_pass": res["hc3_pass"],
        "hc3_n_inferior": hc3.get("n_inferior_conventions"),
        "hc3_n_superior": hc3.get("n_superior_conventions"),
        "global_scale_r_cal": gs.get("r_cal_est"),
        "global_scale_r_orth_cal": gs.get("r_orth_cal_est"),
        "hybrid3d_r_orth_cal": hyb.get("r_orth_cal_est"),
        "tmvs_r_orth_cal": tmvs.get("r_orth_cal_est"),
        "s7_sp_ln_ncells": s7["spearman_ln_ncells"].get("rho_est"),
        "s7_sp_abs_kp": s7["spearman_abs_kp"].get("rho_est"),
        "s7_ols_intercept": ols.get("intercept_adjusted_est"),
        "s7_ols_intercept_ci975_lo": ols.get("intercept_adjusted_ci975_lo"),
        "s7_ols_intercept_ci975_hi": ols.get("intercept_adjusted_ci975_hi"),
        "s7_t1_r_orth_cal": tert.get("T1_low_ncells", {}).get("fct_r_orth_cal"),
        "s7_t1_ci975_lo": tert.get("T1_low_ncells", {}).get("fct_r_orth_cal_ci975_lo"),
    }
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description="Independent Round 10 Group C Audit Runner (A-T3)")
    parser.add_argument("--mode", choices=["r9_dev", "holdout_r10"], required=True)
    parser.add_argument("--target-lines", nargs="+", default=None)
    parser.add_argument("--only-reproducible", action="store_true", default=True)
    parser.add_argument("--all-targets", action="store_true", default=False)
    parser.add_argument("--shrink", type=float, default=0.0)
    parser.add_argument("--cov-max-cells", type=int, default=20000)
    parser.add_argument("--n-boot", type=int, default=4000)
    parser.add_argument("--n-perm", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20261122)
    parser.add_argument("--no-hgnc", action="store_true", default=False)
    parser.add_argument("--out-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--cross-check-csv", type=Path, default=None)
    args = parser.parse_args()

    only_rep = not args.all_targets
    cov_max = args.cov_max_cells if args.cov_max_cells > 0 else None
    hgnc = None if args.no_hgnc else CanonicalHGNC.from_tsv()
    common_genes = load_r9_common_genes(hgnc=hgnc)

    if args.mode == "r9_dev":
        req_targets = args.target_lines if args.target_lines else list(R9_LINE_NAMES)
        enforce_r10_gate(req_targets, allow_r9_only=True, extra_lines_dir=R9_LINES_DIR)
        r9_loaded = {
            ln: load_line_dataset(R9_LINES_DIR / ln, representation="log1p", aligned_genes=common_genes, hgnc=hgnc)
            for ln in R9_LINE_NAMES
        }
        cov_cache = {
            ln: sample_covariance(deterministic_subsample(ds.ctrl_cells, cov_max, seed=42), shrink=args.shrink)
            for ln, ds in r9_loaded.items()
        }
        gene_cache = {ln: ds.genes for ln, ds in r9_loaded.items()}
        line_results = {}
        csv_rows = []
        for tgt in req_targets:
            src_list = [r9_loaded[s] for s in R9_LINE_NAMES if s != tgt]
            res = evaluate_target_line(
                target_data=r9_loaded[tgt],
                source_list=src_list,
                only_reproducible=only_rep,
                shrink=args.shrink,
                cov_max_cells=cov_max,
                n_boot=args.n_boot,
                n_perm=args.n_perm,
                seed=args.seed,
                precomputed_covs=cov_cache,
                precomputed_gene_orders=gene_cache,
            )
            line_results[tgt] = res
            csv_rows.append(_flatten_csv_row(res))
        prefix = "r9_dev_smoke"
        verdict = {"mode": "r9_dev", "note": "Development verification strictly on Round 9 lines; zero R10 data touched."}
    else:
        conf_lines = parse_confirmatory_lines(R10_STATUS_DIR / "CONFIRMATORY_LINES.txt")
        avail_r10 = sorted(
            d.name
            for d in R10_LINES_DIR.iterdir()
            if d.is_dir() and (d / "pseudobulk_delta.npz").exists()
        ) if R10_LINES_DIR.exists() else []
        req_targets = args.target_lines if args.target_lines else avail_r10
        if not req_targets:
            raise RuntimeError(f"No standardized Round 10 lines available in {R10_LINES_DIR}.")
        prereg_h, freeze_h = enforce_r10_gate(req_targets, allow_r9_only=False)
        r9_sources = [
            load_line_dataset(R9_LINES_DIR / ln, representation="log1p", aligned_genes=common_genes, hgnc=hgnc)
            for ln in R9_LINE_NAMES
        ]
        cov_cache = {
            ds.name: sample_covariance(deterministic_subsample(ds.ctrl_cells, cov_max, seed=42), shrink=args.shrink)
            for ds in r9_sources
        }
        gene_cache = {ds.name: ds.genes for ds in r9_sources}
        line_results = {}
        csv_rows = []
        for tgt in req_targets:
            tgt_data = load_line_dataset(R10_LINES_DIR / tgt, representation="log1p", aligned_genes=common_genes, hgnc=hgnc)
            res = evaluate_target_line(
                target_data=tgt_data,
                source_list=r9_sources,
                only_reproducible=only_rep,
                shrink=args.shrink,
                cov_max_cells=cov_max,
                n_boot=args.n_boot,
                n_perm=args.n_perm,
                seed=args.seed,
                precomputed_covs=cov_cache,
                precomputed_gene_orders=gene_cache,
            )
            line_results[tgt] = res
            csv_rows.append(_flatten_csv_row(res))
        active_conf = [c for c in conf_lines if c in line_results]
        verdict = adjudicate_round10_group_c(line_results, active_conf)
        verdict["prereg_commit"] = prereg_h
        verdict["code_freeze_commit"] = freeze_h
        verdict["secondary_lines_evaluated"] = [l for l in req_targets if l in R10_SECONDARY_ONLY or l not in active_conf]
        prefix = "audit_r10" if only_rep else "audit_r10_s1_alltargets"

    args.out_dir.mkdir(parents=True, exist_ok=True)
    if csv_rows:
        with (args.out_dir / f"{prefix}_per_line.csv").open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(csv_rows[0].keys()))
            writer.writeheader()
            writer.writerows(csv_rows)

    cc = run_cross_check(line_results, csv_rows, args.cross_check_csv) if args.cross_check_csv else None
    payload = {"verdict": verdict, "cross_check": cc, "lines": line_results}
    (args.out_dir / f"{prefix}_summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({"prefix": prefix, "verdict": verdict, "cross_check": cc}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
