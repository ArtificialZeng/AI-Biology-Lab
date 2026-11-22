#!/usr/bin/env python3
"""
Round 11 Evaluation Runner (`PREREG_D` Group D `FCT`, Primary Mask `all_nontarget`).

Execution phases:
  - `--phase dev_r9`: Leave-one-line-out software/development check on Round 9's 6 seen lines.
  - `--phase dev_r10`: Zero-shot check on Round 10's already-unblinded lines against R9 sources.
  - `--phase holdout_r11`: Prospective zero-shot confirmatory evaluation of unseen Round 11 lines,
    strictly blocked unless BOTH `PREREG_D_COMMIT.txt` AND `CODE_FREEZE_R11.txt` exist.

Source modes:
  - `--source-mode r9_6` (primary): 6 Round 9 human CRISPRi lines.
  - `--source-mode r9_r10_8` (sensitivity): 6 Round 9 + 2 Round 10 confirmatory (`dld1`, `cd4t_rest`).
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

R11_SRC_DIR = Path(__file__).resolve().parent
R10_SRC_DIR = R11_SRC_DIR.parents[1] / "round10" / "src"
for p in (str(R10_SRC_DIR), str(R11_SRC_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from adjudicate_r10 import SIX_CONVENTIONS  # noqa: E402
from run_eval_r10 import select_eval_targets_r10  # noqa: E402

from adjudicate_r11 import adjudicate_r11_results, render_r11_markdown_report
from data_loader_r11 import (
    DEFAULT_HGNC_FILE,
    HGNCBridge,
    LineDatasetR10,
    R9_COMMON_GENES_FILE,
    R9_LINES_DIR,
    R10_LINES_DIR,
    R11_CONFIRMATORY_LINES_FILE,
    R11_LINES_DIR,
    R11_STATUS_DIR,
    ROUND11_ROOT,
    check_prereg_d_gate,
    list_standardized_lines,
    load_r9_common_genes,
    load_source_cohort_r11,
    load_standardized_line_r10,
    parse_confirmatory_lines_file,
)
from metrics_r11 import (
    compute_fct_perm100_null_r11,
    compute_fct_perm100_split_r11,
    compute_joint_cross_mask_e1_diff_r11,
    compute_joint_paired_r_orth_diff_r11,
    compute_r19_variance_decomposition_r10,
    compute_residual_suite_r11,
    compute_resp_top50_split_suite_r11,
    compute_row_metrics_r10,
    compute_s7_cellcount_and_kd_r11,
    paired_row_bootstrap_ci_r11,
)
from predictors_r11 import (
    FctSourceState,
    build_source_input_r10,
    compute_control_covariance_r10,
    predict_group_d_suite,
    subsample_ctrl_r10,
)

RESULTS_DIR = ROUND11_ROOT / "results"


def _subset_dataset(ds: LineDatasetR10, aligned_genes: List[str], cols: np.ndarray, g_map: Dict[str, int]) -> LineDatasetR10:
    return LineDatasetR10(
        name=ds.name,
        genes=aligned_genes,
        target_names=ds.target_names,
        target_indices=np.asarray([g_map.get(t, -1) for t in ds.target_names], dtype=np.int64),
        deltas=ds.deltas[:, cols],
        ctrl_cells=ds.ctrl_cells[:, cols],
        ctrl_mean_log1p=ds.ctrl_mean_log1p[cols] if ds.ctrl_mean_log1p is not None else None,
        ctrl_mean_libnorm=ds.ctrl_mean_libnorm[cols] if ds.ctrl_mean_libnorm is not None else None,
        ctrl_mean_libnorm_full=(
            ds.ctrl_mean_libnorm_full[cols] if ds.ctrl_mean_libnorm_full is not None else None
        ),
        repr_mask=ds.repr_mask,
        half_a_deltas=ds.half_a_deltas[:, cols] if ds.half_a_deltas is not None else None,
        half_b_deltas=ds.half_b_deltas[:, cols] if ds.half_b_deltas is not None else None,
        half_valid=ds.half_valid,
        n_cells=ds.n_cells,
        knockdown_depth_log1p=ds.knockdown_depth_log1p,
        knockdown_fraction=ds.knockdown_fraction,
        target_ensembl_ids=ds.target_ensembl_ids,
        representation=ds.representation,
        metadata=ds.metadata,
    )


def evaluate_target_against_sources_r11(
    target_ds: LineDatasetR10,
    source_datasets: Dict[str, LineDatasetR10],
    only_reproducible: bool = True,
    s5_top75_genes: bool = False,
    shrink: float = 0.0,
    cov_max_cells: Optional[int] = None,
    include_r9_anchors: bool = True,
    include_exploratory_e: bool = False,
    min_gene_overlap: int = 3000,
    n_boot: int = 4000,
    seed: int = 20261123,
    n_perm: int = 100,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], Dict[str, Any]]:
    """Evaluate one target line `target_ds` strictly zero-shot against `source_datasets` under `PREREG_D`."""
    if target_ds.name in source_datasets:
        raise ValueError(f"Leakage guard: target line {target_ds.name} present inside source_datasets.")

    src_names = sorted(source_datasets.keys())
    shared_genes_set = set(target_ds.genes).intersection(*(set(source_datasets[s].genes) for s in src_names))
    if len(shared_genes_set) < min_gene_overlap:
        raise RuntimeError(
            f"Target line {target_ds.name} gene overlap with sources {len(shared_genes_set)} < {min_gene_overlap}."
        )
    aligned_genes = [g for g in target_ds.genes if g in shared_genes_set]
    g_map = {g: i for i, g in enumerate(aligned_genes)}
    t_pos = {g: i for i, g in enumerate(target_ds.genes)}
    t_cols = np.asarray([t_pos[g] for g in aligned_genes], dtype=np.int64)
    t_sub_ds = _subset_dataset(target_ds, aligned_genes, t_cols, g_map)

    src_inputs = []
    for sn in src_names:
        s_raw = source_datasets[sn]
        s_pos = {g: i for i, g in enumerate(s_raw.genes)}
        s_cols = np.asarray([s_pos[g] for g in aligned_genes], dtype=np.int64)
        src_inputs.append(
            build_source_input_r10(
                _subset_dataset(s_raw, aligned_genes, s_cols, g_map),
                only_reproducible=only_reproducible,
                shrink=shrink,
                cov_max_cells=cov_max_cells,
                compute_cov=include_r9_anchors,
                seed=42,
            )
        )

    state = FctSourceState.fit(src_inputs, mask_self=True)
    row_idx, q_names, q_gidx = select_eval_targets_r10(t_sub_ds, src_inputs, only_reproducible=only_reproducible)
    if len(row_idx) == 0:
        raise RuntimeError(f"Zero evaluated targets found for {target_ds.name}.")
    if t_sub_ds.ctrl_mean_libnorm is None:
        raise ValueError(f"Target line {target_ds.name} missing `ctrl_mean_libnorm`.")

    t_ctrl_cov = None
    if include_r9_anchors:
        sub_c = subsample_ctrl_r10(t_sub_ds.ctrl_cells, cov_max_cells=cov_max_cells, seed=42)
        t_ctrl_cov = compute_control_covariance_r10(sub_c, shrink=shrink)

    preds, fc_bundle, const_stats = predict_group_d_suite(
        state=state,
        target_ctrl_mean_libnorm=t_sub_ds.ctrl_mean_libnorm,
        query_target_names=q_names,
        query_target_gene_indices=q_gidx,
        target_ctrl_cov=t_ctrl_cov,
        target_ctrl_mean_libnorm_full=t_sub_ds.ctrl_mean_libnorm_full,
        include_r9_anchors=include_r9_anchors,
        include_exploratory_e=include_exploratory_e,
    )

    true_deltas = t_sub_ds.deltas[row_idx]
    transfer_ref = preds["Transfer_Mean"]
    ha_sub = t_sub_ds.half_a_deltas[row_idx] if t_sub_ds.half_a_deltas is not None else None
    hb_sub = t_sub_ds.half_b_deltas[row_idx] if t_sub_ds.half_b_deltas is not None else None
    hv_sub = t_sub_ds.half_valid[row_idx] if t_sub_ds.half_valid is not None else None
    n_cells_sub = t_sub_ds.n_cells[row_idx] if t_sub_ds.n_cells is not None else None
    kd_sub = t_sub_ds.knockdown_depth_log1p[row_idx] if t_sub_ds.knockdown_depth_log1p is not None else None

    gene_keep_mask = None
    if s5_top75_genes:
        t_var = np.var(np.asarray(t_sub_ds.ctrl_cells, dtype=np.float64), axis=0, ddof=1)
        gene_keep_mask = t_var > float(np.percentile(t_var, 25.0))

    line_rows: List[Dict[str, Any]] = []
    boot_rows: List[Dict[str, Any]] = []
    row_metrics_by_method: Dict[str, Dict[str, np.ndarray]] = {}

    for m_name, p_mat in preds.items():
        rm_ctrl = compute_row_metrics_r10(
            p_mat, true_deltas, q_gidx, centroid=None, shared_axis=None, mask_target_in_disc=True
        )
        rm_cent = compute_row_metrics_r10(
            p_mat,
            true_deltas,
            q_gidx,
            centroid=state.grand_mean,
            shared_axis=state.axis_shared,
            mask_target_in_disc=True,
        )
        row_metrics_by_method[m_name] = {
            "ctrl_disc_raw": rm_ctrl["disc_raw"],
            "ctrl_disc_raw_panel": rm_ctrl["disc_raw_panel"],
            "ctrl_de50_r": rm_ctrl["de50_r"],
            "cent_disc_raw": rm_cent["disc_raw"],
            "cent_disc_raw_panel": rm_cent["disc_raw_panel"],
            "cent_de50_r": rm_cent["de50_r"],
            "pds_cosine_panel": rm_ctrl.get("pds_cosine_panel", rm_ctrl["disc_raw_panel"]),
            "mse": rm_ctrl["mse"],
        }
        res_suite = compute_residual_suite_r11(
            delta_obs=true_deltas,
            delta_pred=p_mat,
            delta_transfer=transfer_ref,
            target_gene_idx=q_gidx,
            delta_half_a=ha_sub,
            delta_half_b=hb_sub,
            half_valid=hv_sub,
            gene_keep_mask=gene_keep_mask,
            top_k=50,
            n_boot=n_boot,
            seed=seed,
        )
        # Round 11 PRIMARY mask = all_nontarget; SECONDARY mask = resp_top50
        tall = res_suite["all_nontarget"]
        t50 = res_suite["resp_top50"]
        ro = tall.get("r_orth_cal", {})
        ro2 = tall.get("r_orth_cal_sq", {})
        rc = tall.get("r_cal", {})
        sl = tall.get("slope_ols", {})
        sdr = tall.get("sd_pred_over_sd_obs", {})
        dem = tall.get("slope_deming", {})
        rev = tall.get("slope_obs_on_pred", {})

        t50_ro = t50.get("r_orth_cal", {})
        t50_rc = t50.get("r_cal", {})
        t50_sl = t50.get("slope_ols", {})

        split_res = compute_resp_top50_split_suite_r11(
            delta_pred=p_mat,
            delta_transfer=transfer_ref,
            target_gene_idx=q_gidx,
            delta_half_a=ha_sub,
            delta_half_b=hb_sub,
            half_valid=hv_sub,
            gene_keep_mask=gene_keep_mask,
            top_k=50,
            n_boot=n_boot,
            seed=seed,
        )
        sp_cal = split_res.get("r_orth_cal", {})
        sp_raw = split_res.get("r_orth_raw", {})

        line_rows.append(
            {
                "target_line": target_ds.name,
                "method": m_name,
                "n_genes_eval": len(aligned_genes),
                "n_eval_targets": int(len(q_names)),
                "n_sources": len(src_names),
                "primary_mask": "all_nontarget",
                # Primary (all_nontarget)
                "r_orth_cal": ro.get("estimate", float("nan")),
                "r_orth_cal_ci95_lo": ro.get("ci95_lo", float("nan")),
                "r_orth_cal_ci95_hi": ro.get("ci95_hi", float("nan")),
                "r_orth_cal_ci975_lo": ro.get("ci975_lo", float("nan")),
                "r_orth_cal_ci975_hi": ro.get("ci975_hi", float("nan")),
                "r_orth_cal_ci9833_lo": ro.get("ci9833_lo", float("nan")),
                "r_orth_cal_ci9833_hi": ro.get("ci9833_hi", float("nan")),
                "r_orth_cal_valid_frac": ro.get("valid_frac", float("nan")),
                "r_orth_cal_sq": ro2.get("estimate", float("nan")),
                "r_orth_cal_sq_ci975_lo": ro2.get("ci975_lo", float("nan")),
                "r_orth_cal_sq_ci975_hi": ro2.get("ci975_hi", float("nan")),
                "r_orth_raw": tall.get("r_orth_raw", {}).get("estimate", float("nan")),
                "r_cal": rc.get("estimate", float("nan")),
                "r_cal_ci975_lo": rc.get("ci975_lo", float("nan")),
                "r_cal_ci975_hi": rc.get("ci975_hi", float("nan")),
                "r_raw": tall.get("r_raw", {}).get("estimate", float("nan")),
                "slope_ols": sl.get("estimate", float("nan")),
                "slope_ols_ci975_lo": sl.get("ci975_lo", float("nan")),
                "slope_ols_ci975_hi": sl.get("ci975_hi", float("nan")),
                "sd_pred_over_sd_obs": sdr.get("estimate", float("nan")),
                "slope_deming": dem.get("estimate", float("nan")),
                "slope_obs_on_pred": rev.get("estimate", float("nan")),
                "rho_half_ab_orth": tall.get("rho_half_ab_orth", float("nan")),
                "ceiling_r_orth": tall.get("ceiling_r_orth", float("nan")),
                "rho_half_ab": tall.get("rho_half_ab", float("nan")),
                "ceiling_r": tall.get("ceiling_r", float("nan")),
                # Secondary (resp_top50)
                "resp_top50_r_orth_cal": t50_ro.get("estimate", float("nan")),
                "resp_top50_r_orth_cal_ci95_lo": t50_ro.get("ci95_lo", float("nan")),
                "resp_top50_r_orth_cal_ci95_hi": t50_ro.get("ci95_hi", float("nan")),
                "resp_top50_r_orth_cal_ci975_lo": t50_ro.get("ci975_lo", float("nan")),
                "resp_top50_r_orth_cal_ci975_hi": t50_ro.get("ci975_hi", float("nan")),
                "resp_top50_r_orth_cal_ci9833_lo": t50_ro.get("ci9833_lo", float("nan")),
                "resp_top50_r_orth_cal_ci9833_hi": t50_ro.get("ci9833_hi", float("nan")),
                "resp_top50_r_orth_cal_valid_frac": t50_ro.get("valid_frac", float("nan")),
                "resp_top50_r_orth_raw": t50.get("r_orth_raw", {}).get("estimate", float("nan")),
                "resp_top50_r_cal": t50_rc.get("estimate", float("nan")),
                "resp_top50_slope_ols": t50_sl.get("estimate", float("nan")),
                "resp_top50_rho_half_ab_orth": t50.get("rho_half_ab_orth", float("nan")),
                "resp_top50_ceiling_r_orth": t50.get("ceiling_r_orth", float("nan")),
                "resp_top50_rho_half_ab": t50.get("rho_half_ab", float("nan")),
                "resp_top50_ceiling_r": t50.get("ceiling_r", float("nan")),
                # Secondary SA-8 (`resp_top50_split`, single-half ceiling sqrt(rho_AB))
                "resp_top50_split_n_perts": int(split_res.get("n_perts", 0)),
                "resp_top50_split_r_orth_cal": sp_cal.get("estimate", float("nan")),
                "resp_top50_split_r_orth_cal_ci95_lo": sp_cal.get("ci95_lo", float("nan")),
                "resp_top50_split_r_orth_cal_ci95_hi": sp_cal.get("ci95_hi", float("nan")),
                "resp_top50_split_r_orth_cal_ci975_lo": sp_cal.get("ci975_lo", float("nan")),
                "resp_top50_split_r_orth_cal_ci975_hi": sp_cal.get("ci975_hi", float("nan")),
                "resp_top50_split_r_orth_cal_valid_frac": sp_cal.get("valid_frac", float("nan")),
                "resp_top50_split_r_orth_raw": sp_raw.get("estimate", float("nan")),
                "resp_top50_split_r_orth_raw_ci975_lo": sp_raw.get("ci975_lo", float("nan")),
                "resp_top50_split_r_orth_raw_ci975_hi": sp_raw.get("ci975_hi", float("nan")),
                "resp_top50_split_r_orth_cal_a_to_b": split_res.get("r_orth_cal_a_to_b", float("nan")),
                "resp_top50_split_r_orth_cal_b_to_a": split_res.get("r_orth_cal_b_to_a", float("nan")),
                "resp_top50_split_rho_half_ab_orth_a": split_res.get("rho_half_ab_orth_a", float("nan")),
                "resp_top50_split_rho_half_ab_orth_b": split_res.get("rho_half_ab_orth_b", float("nan")),
                "resp_top50_split_ceiling_single_a": split_res.get("ceiling_single_a", float("nan")),
                "resp_top50_split_ceiling_single_b": split_res.get("ceiling_single_b", float("nan")),
                "cos_agreement_top50": res_suite["cos_agreement_top50"],
                **{c: float(np.mean(row_metrics_by_method[m_name][c])) for c in SIX_CONVENTIONS},
                "pds_cosine_panel": float(np.mean(row_metrics_by_method[m_name]["pds_cosine_panel"])),
                "mse": float(np.mean(row_metrics_by_method[m_name]["mse"])),
            }
        )

    # 1. H-D3: Paired within-mask FCT - E1_ConstMed (all_nontarget primary + resp_top50 diagnostic)
    if ha_sub is not None and hb_sub is not None:
        for mode, comp_name in [
            ("all_nontarget", "H-D3_FCT_minus_E1_all_nontarget"),
            ("resp_top50", "diag_FCT_minus_E1_resp_top50"),
        ]:
            hd3_ci = compute_joint_paired_r_orth_diff_r11(
                delta_obs=true_deltas,
                delta_pred_a=preds["FCT"],
                delta_pred_b=preds["E1_ConstMed"],
                delta_transfer=transfer_ref,
                target_gene_idx=q_gidx,
                delta_half_a=ha_sub,
                delta_half_b=hb_sub,
                half_valid=hv_sub,
                gene_keep_mask=gene_keep_mask,
                mode=mode,
                top_k=50,
                n_boot=n_boot,
                seed=seed,
            )
            fct_val = next(r for r in line_rows if r["method"] == "FCT")
            e1_val = next(r for r in line_rows if r["method"] == "E1_ConstMed")
            ma = fct_val["r_orth_cal"] if mode == "all_nontarget" else fct_val["resp_top50_r_orth_cal"]
            mb = e1_val["r_orth_cal"] if mode == "all_nontarget" else e1_val["resp_top50_r_orth_cal"]
            boot_rows.append(
                {
                    "target_line": target_ds.name,
                    "comparison": comp_name,
                    "metric": f"r_orth_cal_{mode}",
                    "method_a": "FCT",
                    "method_b": "E1_ConstMed",
                    "mean_a": float(ma),
                    "mean_b": float(mb),
                    "estimate": hd3_ci["estimate"],
                    "ci95_lo": hd3_ci["ci95_lo"],
                    "ci95_hi": hd3_ci["ci95_hi"],
                    "ci975_lo": hd3_ci["ci975_lo"],
                    "ci975_hi": hd3_ci["ci975_hi"],
                    "ci9833_lo": hd3_ci["ci9833_lo"],
                    "ci9833_hi": hd3_ci["ci9833_hi"],
                    "valid_frac": hd3_ci["valid_frac"],
                    "n_items": hd3_ci["n_items"],
                    "n_boot": hd3_ci["n_boot"],
                }
            )

        # 2. H-D4: Cross-mask paired E1_ConstMed(resp_top50) - E1_ConstMed(all_nontarget)
        hd4_ci = compute_joint_cross_mask_e1_diff_r11(
            delta_obs=true_deltas,
            delta_pred_e1=preds["E1_ConstMed"],
            delta_transfer=transfer_ref,
            target_gene_idx=q_gidx,
            delta_half_a=ha_sub,
            delta_half_b=hb_sub,
            half_valid=hv_sub,
            gene_keep_mask=gene_keep_mask,
            top_k=50,
            n_boot=n_boot,
            seed=seed,
        )
        boot_rows.append(
            {
                "target_line": target_ds.name,
                "comparison": "H-D4_E1_resp_top50_minus_all_nontarget",
                "metric": "r_orth_cal_cross_mask",
                "method_a": "E1_ConstMed(resp_top50)",
                "method_b": "E1_ConstMed(all_nontarget)",
                "mean_a": float(hd4_ci.get("mean_a", float("nan"))),
                "mean_b": float(hd4_ci.get("mean_b", float("nan"))),
                "estimate": hd4_ci["estimate"],
                "ci95_lo": hd4_ci["ci95_lo"],
                "ci95_hi": hd4_ci["ci95_hi"],
                "ci975_lo": hd4_ci["ci975_lo"],
                "ci975_hi": hd4_ci["ci975_hi"],
                "ci9833_lo": hd4_ci["ci9833_lo"],
                "ci9833_hi": hd4_ci["ci9833_hi"],
                "valid_frac": hd4_ci["valid_frac"],
                "n_items": hd4_ci["n_items"],
                "n_boot": hd4_ci["n_boot"],
            }
        )

        # 2b. Descriptive H-D5: FCT minus NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j) under all_nontarget
        kolf_key = "NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)"
        if kolf_key in preds:
            hd5_ci = compute_joint_paired_r_orth_diff_r11(
                delta_obs=true_deltas,
                delta_pred_a=preds["FCT"],
                delta_pred_b=preds[kolf_key],
                delta_transfer=transfer_ref,
                target_gene_idx=q_gidx,
                delta_half_a=ha_sub,
                delta_half_b=hb_sub,
                half_valid=hv_sub,
                gene_keep_mask=gene_keep_mask,
                mode="all_nontarget",
                top_k=50,
                n_boot=n_boot,
                seed=seed,
            )
            fct_val = next(r for r in line_rows if r["method"] == "FCT")
            kolf_val = next(r for r in line_rows if r["method"] == kolf_key)
            boot_rows.append(
                {
                    "target_line": target_ds.name,
                    "comparison": "H-D5_FCT_minus_Mis_kolf21j_all_nontarget",
                    "metric": "r_orth_cal_all_nontarget",
                    "method_a": "FCT",
                    "method_b": kolf_key,
                    "mean_a": float(fct_val["r_orth_cal"]),
                    "mean_b": float(kolf_val["r_orth_cal"]),
                    "estimate": hd5_ci["estimate"],
                    "ci95_lo": hd5_ci["ci95_lo"],
                    "ci95_hi": hd5_ci["ci95_hi"],
                    "ci975_lo": hd5_ci["ci975_lo"],
                    "ci975_hi": hd5_ci["ci975_hi"],
                    "ci9833_lo": hd5_ci["ci9833_lo"],
                    "ci9833_hi": hd5_ci["ci9833_hi"],
                    "valid_frac": hd5_ci["valid_frac"],
                    "n_items": hd5_ci["n_items"],
                    "n_boot": hd5_ci["n_boot"],
                }
            )

    # 3. H-C3-style Practical Utility Non-Inferiority vs Transfer_Mean across 6 conventions
    pair_list = [("FCT", "Transfer_Mean")]
    if "FCT_FullNTMean" in preds:
        pair_list.append(("FCT_FullNTMean", "FCT"))
    for anc in ("Hybrid_ProjShared", "Transfer_Mean_VarScaled", "Hybrid_on_TMVS"):
        if anc in preds:
            pair_list.append((anc, "Transfer_Mean"))
    for ma, mb in pair_list:
        for conv in SIX_CONVENTIONS + ("pds_cosine_panel",):
            ci = paired_row_bootstrap_ci_r11(
                row_metrics_by_method[ma][conv],
                row_metrics_by_method[mb][conv],
                n_boot=n_boot,
                seed=seed,
            )
            boot_rows.append(
                {
                    "target_line": target_ds.name,
                    "comparison": f"H-C3_{ma}_minus_{mb}",
                    "metric": conv,
                    "method_a": ma,
                    "method_b": mb,
                    "mean_a": ci["mean_a"],
                    "mean_b": ci["mean_b"],
                    "estimate": ci["delta"],
                    "ci95_lo": ci["ci95_lo"],
                    "ci95_hi": ci["ci95_hi"],
                    "ci975_lo": ci["ci975_lo"],
                    "ci975_hi": ci["ci975_hi"],
                    "ci9833_lo": ci["ci9833_lo"],
                    "ci9833_hi": ci["ci9833_hi"],
                    "valid_frac": 1.0,
                    "n_items": ci["n_rows"],
                    "n_boot": int(n_boot),
                }
            )

    # 4. 100-Perm Null under all_nontarget (primary), resp_top50 (secondary), AND SA-8 resp_top50_split
    fct_row_summary = next(r for r in line_rows if r["method"] == "FCT")
    perm100_all = compute_fct_perm100_null_r11(
        state=state,
        target_ctrl_mean_libnorm=t_sub_ds.ctrl_mean_libnorm,
        query_target_names=q_names,
        delta_obs=true_deltas,
        delta_transfer=transfer_ref,
        target_gene_idx=q_gidx,
        ceil_r=float(fct_row_summary["ceiling_r"]),
        ceil_r_orth=float(fct_row_summary["ceiling_r_orth"]),
        fc_bundle=fc_bundle,
        mode="all_nontarget",
        n_perm=n_perm,
        top_k=50,
        gene_keep_mask=gene_keep_mask,
    )
    perm100_top50 = compute_fct_perm100_null_r11(
        state=state,
        target_ctrl_mean_libnorm=t_sub_ds.ctrl_mean_libnorm,
        query_target_names=q_names,
        delta_obs=true_deltas,
        delta_transfer=transfer_ref,
        target_gene_idx=q_gidx,
        ceil_r=float(fct_row_summary["resp_top50_ceiling_r"]),
        ceil_r_orth=float(fct_row_summary["resp_top50_ceiling_r_orth"]),
        fc_bundle=fc_bundle,
        mode="resp_top50",
        n_perm=n_perm,
        top_k=50,
        gene_keep_mask=gene_keep_mask,
    )
    perm100_split = compute_fct_perm100_split_r11(
        state=state,
        target_ctrl_mean_libnorm=t_sub_ds.ctrl_mean_libnorm,
        query_target_names=q_names,
        delta_transfer=transfer_ref,
        target_gene_idx=q_gidx,
        delta_half_a=ha_sub,
        delta_half_b=hb_sub,
        half_valid=hv_sub,
        fc_bundle=fc_bundle,
        gene_keep_mask=gene_keep_mask,
        n_perm=n_perm,
        top_k=50,
    )

    r19_decomp = compute_r19_variance_decomposition_r10(
        delta_obs=true_deltas,
        delta_transfer=transfer_ref,
        target_gene_idx=q_gidx,
        min_perts_per_gene=5,
        top_k=50,
        n_boot=n_boot,
        seed=seed,
    )

    mismatch_preds = {
        k: v for k, v in preds.items() if k.startswith("NegCtrl_FCT_MismatchedMean(")
    }
    s7_by_mode = {
        mode_s7: compute_s7_cellcount_and_kd_r11(
            state=state,
            query_target_names=q_names,
            target_ctrl_mean_libnorm=t_sub_ds.ctrl_mean_libnorm,
            fc_bundle=fc_bundle,
            delta_obs=true_deltas,
            delta_pred_fct=preds["FCT"],
            delta_pred_e1=preds["E1_ConstMed"],
            delta_transfer=transfer_ref,
            target_gene_idx=q_gidx,
            n_cells=n_cells_sub,
            knockdown_depth=kd_sub,
            delta_half_a=ha_sub,
            delta_half_b=hb_sub,
            half_valid=hv_sub,
            mismatch_preds=mismatch_preds,
            gene_keep_mask=gene_keep_mask,
            mode=mode_s7,
            top_k=50,
            n_boot=n_boot,
            seed=seed,
        )
        for mode_s7 in ("all_nontarget", "resp_top50")
    }

    extra_line_info: Dict[str, Any] = {
        "n_genes_eval": len(aligned_genes),
        "n_eval_targets": int(len(q_names)),
        "source_lines": src_names,
        "const_stats": const_stats,
        "perm100": {
            "all_nontarget": perm100_all,
            "resp_top50": perm100_top50,
            "resp_top50_split": perm100_split,
        },
        "r19_variance_decomposition": r19_decomp,
        "s7_cellcount_kd": s7_by_mode["all_nontarget"],
        "s7_cellcount_kd_resp_top50": s7_by_mode["resp_top50"],
    }
    if ha_sub is not None and hb_sub is not None and "Hybrid_on_TMVS" in preds and "Transfer_Mean_VarScaled" in preds:
        extra_line_info["r14_hybrid_on_tmvs_vs_tmvs_all_nontarget"] = compute_joint_paired_r_orth_diff_r11(
            delta_obs=true_deltas,
            delta_pred_a=preds["Hybrid_on_TMVS"],
            delta_pred_b=preds["Transfer_Mean_VarScaled"],
            delta_transfer=transfer_ref,
            target_gene_idx=q_gidx,
            delta_half_a=ha_sub,
            delta_half_b=hb_sub,
            half_valid=hv_sub,
            gene_keep_mask=gene_keep_mask,
            mode="all_nontarget",
            top_k=50,
            n_boot=n_boot,
            seed=seed,
        )

    return line_rows, boot_rows, extra_line_info


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Round 11 Group D (`FCT`, `all_nontarget` primary) runner.")
    p.add_argument("--phase", choices=["dev_r9", "dev_r10", "holdout_r11"], required=True)
    p.add_argument("--source-mode", choices=["r9_6", "r9_r10_8"], default="r9_6")
    p.add_argument("--repr", choices=["log1p", "zscore", "libnorm"], default="log1p")
    p.add_argument("--all-targets", action="store_true", default=False, help="S1 sensitivity: include non-reproducible targets.")
    p.add_argument("--s5-top75", action="store_true", default=False, help="S5 sensitivity: restrict mask to top-75%% ctrl-var genes.")
    p.add_argument("--no-r9-anchors", action="store_true", default=False)
    p.add_argument("--no-exploratory-e", action="store_true", default=False, help="Skip descriptive E2/E3/E5/E6/E7 ablations.")
    p.add_argument("--no-canon-hgnc", action="store_true", default=False)
    p.add_argument("--lines", type=str, default="", help="Optional comma-separated target lines.")
    p.add_argument("--n-boot", type=int, default=4000)
    p.add_argument("--n-perm", type=int, default=100)
    p.add_argument("--seed", type=int, default=20261123)
    p.add_argument("--tag-suffix", type=str, default="")
    p.add_argument("--out-dir", type=Path, default=RESULTS_DIR)
    return p.parse_args()


def main() -> None:
    t0 = time.time()
    args = parse_args()
    only_rep = False if args.all_targets else True
    include_anchors = not args.no_r9_anchors
    hgnc = None if args.no_canon_hgnc else HGNCBridge.from_tsv(DEFAULT_HGNC_FILE)
    common_r9 = load_r9_common_genes(hgnc_bridge=hgnc, common_genes_file=R9_COMMON_GENES_FILE)

    prereg_hash = ""
    freeze_hash = ""
    confirmatory_lines: List[str] = []
    source_pool_names: List[str] = []

    if args.phase == "dev_r9":
        avail_r9 = list_standardized_lines(R9_LINES_DIR)
        if not avail_r9:
            raise RuntimeError(f"Missing Round 9 lines under {R9_LINES_DIR}.")
        all_ds = {
            ln: load_standardized_line_r10(
                R9_LINES_DIR / ln, representation=args.repr, common_genes=common_r9, hgnc_bridge=hgnc
            )
            for ln in avail_r9
        }
        eval_lines = [x.strip() for x in args.lines.split(",") if x.strip()] if args.lines.strip() else avail_r9
        source_pool_names = avail_r9
        target_ds_map = all_ds
        src_ds_base = all_ds
        confirmatory_lines = []
    elif args.phase == "dev_r10":
        src_ds_base = load_source_cohort_r11(
            args.source_mode, representation=args.repr, common_genes=common_r9, hgnc_bridge=hgnc
        )
        avail_r10 = list_standardized_lines(R10_LINES_DIR)
        if not avail_r10:
            raise RuntimeError(f"Missing Round 10 lines under {R10_LINES_DIR}.")
        target_ds_map = {
            ln: load_standardized_line_r10(
                R10_LINES_DIR / ln, representation=args.repr, common_genes=common_r9, hgnc_bridge=hgnc
            )
            for ln in avail_r10
        }
        eval_lines = [x.strip() for x in args.lines.split(",") if x.strip()] if args.lines.strip() else avail_r10
        source_pool_names = sorted(src_ds_base.keys())
        confirmatory_lines = []
    else:
        prereg_hash, freeze_hash = check_prereg_d_gate(R11_STATUS_DIR)
        src_ds_base = load_source_cohort_r11(
            args.source_mode, representation=args.repr, common_genes=common_r9, hgnc_bridge=hgnc
        )
        avail_r11 = list_standardized_lines(R11_LINES_DIR)
        if not avail_r11:
            raise RuntimeError(f"No processed Round 11 target lines found under {R11_LINES_DIR}.")
        target_ds_map = {
            ln: load_standardized_line_r10(
                R11_LINES_DIR / ln,
                representation=args.repr,
                common_genes=common_r9,
                hgnc_bridge=hgnc,
                min_gene_overlap=(50 if ln == "feng_ipsc" else 3000),
            )
            for ln in avail_r11
        }
        eval_lines = [x.strip() for x in args.lines.split(",") if x.strip()] if args.lines.strip() else avail_r11
        source_pool_names = sorted(src_ds_base.keys())
        confirmatory_lines = parse_confirmatory_lines_file(R11_CONFIRMATORY_LINES_FILE)

    all_line_rows: List[Dict[str, Any]] = []
    all_boot_rows: List[Dict[str, Any]] = []
    perm_map: Dict[str, Any] = {}
    r19_map: Dict[str, Any] = {}
    s7_map: Dict[str, Any] = {}
    s7_top50_map: Dict[str, Any] = {}
    r14_map: Dict[str, Any] = {}
    eval_meta: Dict[str, Any] = {}

    for tgt_ln in eval_lines:
        tgt_ds = target_ds_map[tgt_ln]
        src_ds = {k: v for k, v in src_ds_base.items() if k != tgt_ln}

        l_rows, b_rows, extra = evaluate_target_against_sources_r11(
            target_ds=tgt_ds,
            source_datasets=src_ds,
            only_reproducible=only_rep,
            s5_top75_genes=args.s5_top75,
            include_r9_anchors=include_anchors,
            include_exploratory_e=not args.no_exploratory_e,
            min_gene_overlap=(50 if tgt_ln == "feng_ipsc" else 3000),
            n_boot=args.n_boot,
            seed=args.seed,
            n_perm=args.n_perm,
        )
        all_line_rows.extend(l_rows)
        all_boot_rows.extend(b_rows)
        perm_map[tgt_ln] = extra["perm100"]
        r19_map[tgt_ln] = extra["r19_variance_decomposition"]
        s7_map[tgt_ln] = extra["s7_cellcount_kd"]
        s7_top50_map[tgt_ln] = extra["s7_cellcount_kd_resp_top50"]
        if "r14_hybrid_on_tmvs_vs_tmvs_all_nontarget" in extra:
            r14_map[tgt_ln] = extra["r14_hybrid_on_tmvs_vs_tmvs_all_nontarget"]
        eval_meta[tgt_ln] = {
            "n_genes_eval": extra["n_genes_eval"],
            "n_eval_targets": extra["n_eval_targets"],
            "source_lines": extra["source_lines"],
            "const_stats": extra["const_stats"],
        }

    df_lines = pd.DataFrame(all_line_rows)
    df_boot = pd.DataFrame(all_boot_rows)
    summary: Dict[str, Any] = {
        "phase": args.phase,
        "source_mode": args.source_mode,
        "source_lines": source_pool_names,
        "representation": args.repr,
        "primary_mask": "all_nontarget",
        "secondary_mask": "resp_top50",
        "only_reproducible": only_rep,
        "s5_top75": args.s5_top75,
        "canon_hgnc": not args.no_canon_hgnc,
        "prereg_commit_hash": prereg_hash,
        "code_freeze_hash": freeze_hash,
        "confirmatory_lines": confirmatory_lines,
        "bootstrap_seed": int(args.seed),
        "n_boot": int(args.n_boot),
        "n_perm": int(args.n_perm),
        "eval_meta": eval_meta,
        "fct_perm100": perm_map,
        "r19_variance_decomposition": r19_map,
        "s7_cellcount_kd": s7_map,
        "s7_cellcount_kd_resp_top50": s7_top50_map,
        "r14_hybrid_on_tmvs_vs_tmvs_all_nontarget": r14_map,
        "elapsed_sec": round(time.time() - t0, 2),
    }

    args.out_dir.mkdir(parents=True, exist_ok=True)
    tag = f"{args.phase}_{args.source_mode}_{args.repr}{args.tag_suffix}"
    df_lines.to_csv(args.out_dir / f"line_metrics_{tag}.csv", index=False)
    df_boot.to_csv(args.out_dir / f"bootstrap_ci_{tag}.csv", index=False)
    (args.out_dir / f"summary_{tag}.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    adj = adjudicate_r11_results(df_lines, df_boot, summary, confirmatory_lines=confirmatory_lines)
    (args.out_dir / f"adjudication_{tag}.json").write_text(json.dumps(adj, indent=2), encoding="utf-8")
    (args.out_dir / f"REPORT_{tag}.md").write_text(
        render_r11_markdown_report(df_lines, df_boot, summary, adj, tag), encoding="utf-8"
    )
    print(f"Saved Round 11 results (*_{tag}.*) in {summary['elapsed_sec']}s", flush=True)


if __name__ == "__main__":
    main()
