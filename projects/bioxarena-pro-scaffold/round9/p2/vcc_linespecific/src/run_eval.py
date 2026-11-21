#!/usr/bin/env python3
"""
Leave-one-cell-line-out zero-shot evaluation runner for line-specific perturbation response.

Phases:
  --phase dev      : Runs ONLY on Replogle development pair (K562 <-> RPE1).
  --phase holdout  : Runs LOO across all lines in `p2/shared/lines/`.
                     Hard-blocked unless `p2/status/PREREG_COMMIT.txt` exists and is non-empty.
"""

from __future__ import annotations

import argparse
import gc
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

SRC_DIR = Path(__file__).resolve().parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from data_loader import (
    LineDataset,
    check_prereg_gate,
    is_replogle_dev_line,
    list_available_lines,
    load_common_genes,
    load_line_dataset,
)
from metrics import (
    cis_non_overlapping,
    compute_line_specific_metrics_and_boots,
    compute_row_metrics,
    compute_splithalf_noise_ceiling,
    drf_score,
    paired_row_bootstrap_ci,
    pearson_1d,
    single_array_bootstrap_ci,
    spearman_brown_ceiling,
    stratified_row_bootstrap_ci,
)
from predictors import (
    SourceLineInput,
    SourceModelState,
    cipher_unit_response_columns,
    compute_control_covariance,
    predict_zero_shot_suite,
)

RESULTS_DIR = SRC_DIR.parent / "results"
AUDIT_LIT_DIR = SRC_DIR.parent / "audit" / "LIT"


def control_splithalf_covariance_reliability(
    ctrl_cells: np.ndarray,
    target_indices: Sequence[int],
    seed: int = 42,
    shrink: float = 0.0,
    min_var: float = 1e-8,
) -> Dict[str, float]:
    """
    Estimate noise ceiling/reliability of target control covariance vectors C_t[:, g] / C_t[g, g]
    by splitting target control cells 50/50 into Half A and Half B.
    Uses ONLY unperturbed control cells of target line t.
    """
    n = ctrl_cells.shape[0]
    rng = np.random.RandomState(seed)
    perm = rng.permutation(n)
    mid = n // 2
    cov_a = compute_control_covariance(ctrl_cells[perm[:mid]], shrink=shrink)
    cov_b = compute_control_covariance(ctrl_cells[perm[mid:]], shrink=shrink)
    r_a = cipher_unit_response_columns(cov_a, target_indices, min_var=min_var)
    r_b = cipher_unit_response_columns(cov_b, target_indices, min_var=min_var)
    corrs = []
    for i, ti in enumerate(target_indices):
        mask = np.ones(r_a.shape[1], dtype=bool)
        if 0 <= int(ti) < r_a.shape[1]:
            mask[int(ti)] = False
        corrs.append(pearson_1d(r_a[i, mask], r_b[i, mask]))
    arr = np.asarray(corrs, dtype=np.float64)
    rho_ab = float(np.mean(arr))
    rho_sb, ceil_r = spearman_brown_ceiling(rho_ab)
    return {
        "n_ctrl_cells": int(n),
        "ctrl_cov_splithalf_r_mean": rho_ab,
        "ctrl_cov_splithalf_r_median": float(np.median(arr)),
        "ctrl_cov_spearman_brown": rho_sb,
        "ctrl_cov_ceiling_r": ceil_r,
    }


def _subsample_ctrl(ctrl_cells: np.ndarray, cov_max_cells: Optional[int], seed: int = 42) -> np.ndarray:
    if cov_max_cells is None or ctrl_cells.shape[0] <= cov_max_cells:
        return ctrl_cells
    rng = np.random.RandomState(seed)
    sel = np.sort(rng.choice(ctrl_cells.shape[0], size=cov_max_cells, replace=False))
    return ctrl_cells[sel]


def build_source_input(
    ds: LineDataset,
    shrink: float = 0.0,
    only_reproducible: bool = False,
    cov_max_cells: Optional[int] = None,
) -> SourceLineInput:
    """Construct SourceLineInput from a loaded source LineDataset."""
    ctrl_sub = _subsample_ctrl(ds.ctrl_cells, cov_max_cells=cov_max_cells, seed=42)
    cov = compute_control_covariance(ctrl_sub, shrink=shrink)
    if only_reproducible and ds.repr_mask is not None:
        keep = ds.repr_mask & (ds.target_indices >= 0)
    else:
        keep = ds.target_indices >= 0
    return SourceLineInput(
        name=ds.name,
        ctrl_cov=cov,
        deltas=ds.deltas[keep],
        target_indices=ds.target_indices[keep],
        target_names=tuple(ds.target_names[i] for i in np.where(keep)[0]),
        ctrl_mean_libnorm=ds.ctrl_mean_libnorm if ds.representation == "log1p" else None,
    )


def _evaluate_s5_top75_var(
    preds: Dict[str, np.ndarray],
    true_deltas: np.ndarray,
    transfer_ref: np.ndarray,
    target_ctrl_cov: np.ndarray,
    q_gidx: np.ndarray,
    ha_sub: Optional[np.ndarray],
    hb_sub: Optional[np.ndarray],
    obs_sub: Optional[np.ndarray],
    tr_sub: Optional[np.ndarray],
    hv_mask: Optional[np.ndarray],
    n_boot: int,
    seed: int,
) -> Dict[str, Any]:
    """PREREG_AUDIT_R2 S5: Recompute r_cal and slope after dropping bottom 25% target ctrl variance genes."""
    diag_t = np.diag(np.asarray(target_ctrl_cov, dtype=np.float64))
    G_dim = diag_t.shape[0]
    order = np.argsort(diag_t)
    keep_cols = np.zeros(G_dim, dtype=bool)
    keep_cols[order[G_dim // 4 :]] = True
    old_to_new = np.full(G_dim, -1, dtype=np.int64)
    old_to_new[keep_cols] = np.arange(int(np.sum(keep_cols)), dtype=np.int64)
    q_gidx_s5 = np.array(
        [int(old_to_new[ti]) if 0 <= int(ti) < G_dim else -1 for ti in q_gidx],
        dtype=np.int64,
    )

    ceil_s5 = None
    if ha_sub is not None and hb_sub is not None and obs_sub is not None and tr_sub is not None:
        q_hv = q_gidx_s5[hv_mask] if hv_mask is not None else q_gidx_s5
        ceil_s5 = compute_splithalf_noise_ceiling(
            delta_half_a=ha_sub[:, keep_cols],
            delta_half_b=hb_sub[:, keep_cols],
            delta_obs_full=obs_sub[:, keep_cols],
            delta_src_mean=tr_sub[:, keep_cols],
            target_gene_idx=q_hv,
            top_k=50,
            n_boot=n_boot,
            seed=seed,
        )

    s5_methods = [
        "Hybrid_ProjShared",
        "Transfer_Mean_VarScaled",
        "Hybrid_on_TMVS",
        "FCT",
        "NegCtrl_HybridProj_SrcCorrTgtVar",
        "NegCtrl_TMVS_ShuffledVar",
    ]
    out_methods: Dict[str, Any] = {}
    d_obs_s5 = true_deltas[:, keep_cols]
    tr_s5 = transfer_ref[:, keep_cols]
    for m_name in s5_methods:
        if m_name not in preds:
            continue
        ls_s5 = compute_line_specific_metrics_and_boots(
            delta_obs=d_obs_s5,
            delta_pred=preds[m_name][:, keep_cols],
            delta_src_mean=tr_s5,
            target_gene_idx=q_gidx_s5,
            top_k=50,
            n_boot=n_boot,
            seed=seed,
        )["resp_top50"]
        ceil_val = ceil_s5["resp_top50"]["ceiling_r"] if ceil_s5 is not None else float("nan")
        cal_r = (
            ls_s5["resid_corr"] / ceil_val
            if (np.isfinite(ceil_val) and ceil_val > 0)
            else float("nan")
        )
        ci95 = [float("nan"), float("nan")]
        ci975 = [float("nan"), float("nan")]
        vfrac = 0.0
        if ceil_s5 is not None:
            b_cal = ls_s5["boot_resid_corr"] / ceil_s5["resp_top50"]["boot_ceiling_r"]
            vfrac = float(np.mean(np.isfinite(b_cal)))
            if np.isfinite(cal_r) and ceil_val > 0 and vfrac >= 0.95:
                ci95 = np.nanpercentile(b_cal, [2.5, 97.5]).tolist()
                ci975 = np.nanpercentile(b_cal, [1.25, 98.75]).tolist()
        out_methods[m_name] = {
            "ls_slope_ols_top50": ls_s5["slope_ols"],
            "ls_slope_ols_top50_ci95": ls_s5["slope_ols_ci95"],
            "ls_slope_ols_top50_ci975": ls_s5["slope_ols_ci975"],
            "ls_resid_corr_top50": ls_s5["resid_corr"],
            "ls_resid_corr_cal_top50": cal_r,
            "ls_resid_corr_cal_top50_sq": float(cal_r ** 2) if np.isfinite(cal_r) else float("nan"),
            "ls_resid_corr_cal_top50_resamp_ci95": ci95,
            "ls_resid_corr_cal_top50_resamp_ci975": ci975,
            "ls_resid_corr_cal_top50_resamp_valid_frac": vfrac,
        }
    return {
        "n_genes_total": int(G_dim),
        "n_genes_kept": int(np.sum(keep_cols)),
        "rho_half_ab_top50": ceil_s5["resp_top50"]["rho_half_ab"] if ceil_s5 else float("nan"),
        "ceiling_r_top50": ceil_s5["resp_top50"]["ceiling_r"] if ceil_s5 else float("nan"),
        "rho_half_ab_all": ceil_s5["all_nontarget"]["rho_half_ab"] if ceil_s5 else float("nan"),
        "ceiling_r_all": ceil_s5["all_nontarget"]["ceiling_r"] if ceil_s5 else float("nan"),
        "methods": out_methods,
    }


def select_eval_targets(
    target_ds: LineDataset,
    source_inputs: Sequence[SourceLineInput],
    shared_targets_only: bool = True,
    only_reproducible: bool = False,
) -> Tuple[np.ndarray, List[str], np.ndarray]:
    """Select evaluation query target rows for held-out target line."""
    valid = target_ds.target_indices >= 0
    if only_reproducible and target_ds.repr_mask is not None:
        valid = valid & target_ds.repr_mask
    if shared_targets_only:
        src_targets = set()
        for s in source_inputs:
            src_targets.update(s.target_names)
        in_src = np.array([t in src_targets for t in target_ds.target_names], dtype=bool)
        valid = valid & in_src
    row_idx = np.where(valid)[0]
    q_names = [target_ds.target_names[i] for i in row_idx]
    q_gidx = target_ds.target_indices[row_idx]
    return row_idx, q_names, q_gidx


def evaluate_datasets(
    datasets: Dict[str, LineDataset],
    shrink: float = 0.0,
    cov_max_cells: Optional[int] = None,
    mask_self: bool = True,
    shared_targets_only: bool = True,
    only_reproducible: bool = False,
    n_boot: int = 4000,
    seed: int = 12345,
    require_official_spec: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, Any]]:
    """Run full leave-one-cell-line-out evaluation across `datasets`."""
    if require_official_spec:
        lit_metric = AUDIT_LIT_DIR / "LINE_SPECIFIC_METRIC.md"
        if not lit_metric.exists():
            raise RuntimeError(f"Missing official metric spec {lit_metric}.")

    line_names = sorted(datasets.keys())
    all_source_inputs = {
        ln: build_source_input(
            datasets[ln],
            shrink=shrink,
            only_reproducible=only_reproducible,
            cov_max_cells=cov_max_cells,
        )
        for ln in line_names
    }

    line_rows: List[Dict[str, Any]] = []
    boot_rows: List[Dict[str, Any]] = []
    pooled_row_scores: Dict[str, Dict[str, List[float]]] = {}
    pooled_strata: List[str] = []
    ctrl_reliability: Dict[str, Dict[str, float]] = {}
    pert_noise_ceilings: Dict[str, Dict[str, Any]] = {}
    drf_calibration: Dict[str, Dict[str, Any]] = {}
    lam_by_target_line: Dict[str, float] = {}
    mismatched_cov_sources_by_line: Dict[str, List[str]] = {}
    s5_expressed_top75: Dict[str, Dict[str, Any]] = {}
    saved_eval_by_line: Dict[str, Dict[str, Any]] = {}

    comparisons = [
        ("Hybrid_Raw", "Transfer_Mean"),
        ("Hybrid_ProjShared", "Transfer_Mean"),
        ("Hybrid_ProjShared", "Transfer_Mean_VarScaled"),
        ("Hybrid_ProjShared", "NegCtrl_HybridProj_SrcCorrTgtVar"),
        ("Hybrid_ProjShared", "NegCtrl_HybridProj_ShuffledGenes"),
        ("Hybrid_ProjShared", "NegCtrl_HybridProj_DiagCov"),
        ("CIPHER_Linear", "Transfer_Mean"),
        ("CIPHER_Linear", "NegCtrl_CIPHER_ShuffledGenes"),
        ("CIPHER_Linear", "NegCtrl_CIPHER_DiagCov"),
        ("Hybrid_Raw", "NegCtrl_Hybrid_ShuffledGenes"),
        ("NegCtrl_HybridProj_SrcCorrTgtVar", "Transfer_Mean_VarScaled"),
        # PREREG_AUDIT_R2 Group B (TMVS) & R14/R15 comparisons
        ("Transfer_Mean_VarScaled", "Transfer_Mean"),
        ("Transfer_Mean_VarScaled", "NegCtrl_TMVS_ShuffledVar"),
        ("Transfer_Mean_VarScaled", "FCT"),
        ("Hybrid_on_TMVS", "Transfer_Mean_VarScaled"),
        ("Hybrid_on_TMVS", "Transfer_Mean"),
    ]

    for t_idx_line, target_line in enumerate(line_names):
        t_ds = datasets[target_line]
        src_names = [ln for ln in line_names if ln != target_line]
        src_inputs = [all_source_inputs[s] for s in src_names]
        state = SourceModelState.fit(src_inputs, mask_self=mask_self)
        lam_by_target_line[target_line] = float(state.lam)

        row_idx, q_names, q_gidx = select_eval_targets(
            t_ds, src_inputs, shared_targets_only=shared_targets_only, only_reproducible=only_reproducible
        )
        if len(row_idx) == 0:
            raise RuntimeError(f"Zero evaluation targets found for target line {target_line}")

        # Target-line control covariance ONLY (zero perturbation leakage)
        t_ctrl_sub = _subsample_ctrl(t_ds.ctrl_cells, cov_max_cells=cov_max_cells, seed=42)
        t_ctrl_cov = compute_control_covariance(t_ctrl_sub, shrink=shrink)
        mismatched_cov_sources_by_line[target_line] = list(src_names)
        preds = predict_zero_shot_suite(
            state=state,
            target_ctrl_cov=t_ctrl_cov,
            query_target_names=q_names,
            query_target_gene_indices=q_gidx,
            mismatched_ctrl_cov=None,
            shuffle_seed=42,
            target_ctrl_mean_libnorm=(
                t_ds.ctrl_mean_libnorm if t_ds.representation == "log1p" else None
            ),
        )
        true_deltas = t_ds.deltas[row_idx]
        transfer_ref = preds["Transfer_Mean"]
        _cos_keep = {
            "Transfer_Mean",
            "CIPHER_Linear",
            "Hybrid_Raw",
            "Hybrid_ProjShared",
            "Transfer_Mean_VarScaled",
            "Hybrid_on_TMVS",
            "FCT",
            "NegCtrl_HybridProj_SrcCorrTgtVar",
        }
        saved_eval_by_line[target_line] = {
            "target_to_row": {t: idx for idx, t in enumerate(q_names)},
            "q_gidx": q_gidx,
            "true_deltas": true_deltas.astype(np.float32),
            "preds": {k: v.astype(np.float32) for k, v in preds.items() if k in _cos_keep},
        }

        ctrl_reliability[target_line] = control_splithalf_covariance_reliability(
            t_ctrl_sub, q_gidx, seed=42, shrink=shrink
        )

        # Perturbation split-half noise ceiling & DRF positive controls (when half A/B available)
        ceil_info = None
        pos_metrics_cent = None
        pos_metrics_ctrl = None
        ha_sub = hb_sub = obs_sub = tr_sub = hv_mask_eval = None
        if t_ds.half_a_deltas is not None and t_ds.half_b_deltas is not None:
            hv_mask = (
                t_ds.half_valid[row_idx]
                if t_ds.half_valid is not None
                else np.ones(len(row_idx), dtype=bool)
            )
            if np.sum(hv_mask) >= 10:
                hv_mask_eval = hv_mask
                ha_sub = t_ds.half_a_deltas[row_idx][hv_mask]
                hb_sub = t_ds.half_b_deltas[row_idx][hv_mask]
                obs_sub = true_deltas[hv_mask]
                tr_sub = transfer_ref[hv_mask]
                qg_sub = q_gidx[hv_mask]
                ceil_raw = compute_splithalf_noise_ceiling(
                    delta_half_a=ha_sub,
                    delta_half_b=hb_sub,
                    delta_obs_full=obs_sub,
                    delta_src_mean=tr_sub,
                    target_gene_idx=qg_sub,
                    top_k=50,
                    n_boot=n_boot,
                    seed=seed,
                )
                ceil_info = ceil_raw
                pert_noise_ceilings[target_line] = {
                    mode: {
                        k: v
                        for k, v in d.items()
                        if not isinstance(v, np.ndarray)
                    }
                    for mode, d in ceil_raw.items()
                }
                pert_noise_ceilings[target_line]["n_half_valid_targets"] = int(np.sum(hv_mask))

                # DRF positive control: Half A predicting Half B
                pos_metrics_cent = compute_row_metrics(
                    ha_sub,
                    hb_sub,
                    qg_sub,
                    centroid=state.grand_mean,
                    shared_axis=state.axis_shared,
                    mask_target_in_disc=True,
                )
                pos_metrics_ctrl = compute_row_metrics(
                    ha_sub,
                    hb_sub,
                    qg_sub,
                    centroid=None,
                    shared_axis=None,
                    mask_target_in_disc=True,
                )

        row_metrics_by_method: Dict[str, Dict[str, np.ndarray]] = {}
        ls_by_method: Dict[str, Dict[str, Dict[str, Any]]] = {}

        for m_name, p_mat in preds.items():
            rm_cent = compute_row_metrics(
                p_mat,
                true_deltas,
                q_gidx,
                centroid=state.grand_mean,
                shared_axis=state.axis_shared,
                mask_target_in_disc=True,
            )
            rm_ctrl = compute_row_metrics(
                p_mat,
                true_deltas,
                q_gidx,
                centroid=None,
                shared_axis=None,
                mask_target_in_disc=True,
            )
            ls_res = compute_line_specific_metrics_and_boots(
                delta_obs=true_deltas,
                delta_pred=p_mat,
                delta_src_mean=transfer_ref,
                target_gene_idx=q_gidx,
                top_k=50,
                n_boot=n_boot,
                seed=seed,
            )
            ls_by_method[m_name] = ls_res

            merged = {
                "mse": rm_cent["mse"],
                "cent_delta_r": rm_cent["delta_r"],
                "cent_de50_r": rm_cent["de50_r"],
                "cent_de50_sign": rm_cent["de50_sign"],
                "cent_disc_raw": rm_cent["disc_raw"],
                "cent_disc_raw_panel": rm_cent["disc_raw_panel"],
                "cent_v5": rm_cent["v5"],
                "cent_disc_v5": rm_cent["disc_v5"],
                "ctrl_delta_r": rm_ctrl["delta_r"],
                "ctrl_de50_r": rm_ctrl["de50_r"],
                "ctrl_de50_sign": rm_ctrl["de50_sign"],
                "ctrl_disc_raw": rm_ctrl["disc_raw"],
                "ctrl_disc_raw_panel": rm_ctrl["disc_raw_panel"],
                "pds_cosine_panel": rm_ctrl.get("pds_cosine_panel", rm_ctrl["disc_raw"]),
            }
            row_metrics_by_method[m_name] = merged
            pooled_row_scores.setdefault(m_name, {})
            for k, arr in merged.items():
                pooled_row_scores[m_name].setdefault(k, []).extend(arr.tolist())

            top50_ls = ls_res["resp_top50"]
            all_ls = ls_res["all_nontarget"]

            # Noise-calibrated correlation & r_cal^2 (H-A2, R11) if split-half ceiling is positive
            ceil_top50 = ceil_info["resp_top50"]["ceiling_r"] if ceil_info else float("nan")
            ceil_all = ceil_info["all_nontarget"]["ceiling_r"] if ceil_info else float("nan")
            cal_r_top50 = (
                top50_ls["resid_corr"] / ceil_top50
                if (np.isfinite(ceil_top50) and ceil_top50 > 0)
                else float("nan")
            )
            cal_r_sq_top50 = float(cal_r_top50 ** 2) if np.isfinite(cal_r_top50) else float("nan")
            cal_r_top50_ci_fixed = (
                [c / ceil_top50 for c in top50_ls["resid_corr_ci95"]]
                if (np.isfinite(ceil_top50) and ceil_top50 > 0)
                else [float("nan"), float("nan")]
            )
            if np.isfinite(ceil_top50) and ceil_top50 > 0 and np.any(np.isfinite(top50_ls["boot_resid_corr"])):
                cal_r_sq_top50_ci_fixed = np.nanpercentile(
                    np.square(top50_ls["boot_resid_corr"] / ceil_top50), [2.5, 97.5]
                ).tolist()
            else:
                cal_r_sq_top50_ci_fixed = [float("nan"), float("nan")]

            cal_r_top50_ci_resamp = [float("nan"), float("nan")]
            cal_r_top50_ci975_resamp = [float("nan"), float("nan")]
            cal_r_sq_top50_ci_resamp = [float("nan"), float("nan")]
            cal_r_sq_top50_ci975_resamp = [float("nan"), float("nan")]
            cal_r_top50_resamp_valid_frac = 0.0
            boot_cal_top50 = np.full_like(top50_ls["boot_resid_corr"], np.nan)
            if ceil_info is not None:
                b_cal = top50_ls["boot_resid_corr"] / ceil_info["resp_top50"]["boot_ceiling_r"]
                cal_r_top50_resamp_valid_frac = float(np.mean(np.isfinite(b_cal)))
                # PREREG_A §5.2: require >= 95% finite joint-bootstrap draws & positive full-sample ceiling
                if (
                    np.isfinite(cal_r_top50)
                    and ceil_top50 > 0
                    and cal_r_top50_resamp_valid_frac >= 0.95
                ):
                    boot_cal_top50 = b_cal
                    cal_r_top50_ci_resamp = np.nanpercentile(b_cal, [2.5, 97.5]).tolist()
                    cal_r_top50_ci975_resamp = np.nanpercentile(b_cal, [1.25, 98.75]).tolist()
                    cal_r_sq_top50_ci_resamp = np.nanpercentile(np.square(b_cal), [2.5, 97.5]).tolist()
                    cal_r_sq_top50_ci975_resamp = np.nanpercentile(np.square(b_cal), [1.25, 98.75]).tolist()
            top50_ls["resid_corr_cal"] = cal_r_top50
            top50_ls["boot_resid_corr_cal"] = boot_cal_top50

            cal_r_all = (
                all_ls["resid_corr"] / ceil_all
                if (np.isfinite(ceil_all) and ceil_all > 0)
                else float("nan")
            )
            cal_r_sq_all = float(cal_r_all ** 2) if np.isfinite(cal_r_all) else float("nan")
            cal_r_all_ci_fixed = (
                [c / ceil_all for c in all_ls["resid_corr_ci95"]]
                if (np.isfinite(ceil_all) and ceil_all > 0)
                else [float("nan"), float("nan")]
            )
            cal_r_all_ci_resamp = [float("nan"), float("nan")]
            cal_r_all_ci975_resamp = [float("nan"), float("nan")]
            cal_r_sq_all_ci_resamp = [float("nan"), float("nan")]
            cal_r_all_resamp_valid_frac = 0.0
            if ceil_info is not None:
                b_cal_all = all_ls["boot_resid_corr"] / ceil_info["all_nontarget"]["boot_ceiling_r"]
                cal_r_all_resamp_valid_frac = float(np.mean(np.isfinite(b_cal_all)))
                if (
                    np.isfinite(cal_r_all)
                    and ceil_all > 0
                    and cal_r_all_resamp_valid_frac >= 0.95
                ):
                    cal_r_all_ci_resamp = np.nanpercentile(b_cal_all, [2.5, 97.5]).tolist()
                    cal_r_all_ci975_resamp = np.nanpercentile(b_cal_all, [1.25, 98.75]).tolist()
                    cal_r_sq_all_ci_resamp = np.nanpercentile(np.square(b_cal_all), [2.5, 97.5]).tolist()

            rec: Dict[str, Any] = {
                "target_line": target_line,
                "n_eval_targets": int(len(row_idx)),
                "method": m_name,
                "lambda_src": float(state.lam),
                # Primary & sensitivity line-specific metrics (LINE_SPECIFIC_METRIC.md & R6/R11/R13)
                "ls_slope_ols_top50": top50_ls["slope_ols"],
                "ls_slope_ols_top50_ci95_lo": top50_ls["slope_ols_ci95"][0],
                "ls_slope_ols_top50_ci95_hi": top50_ls["slope_ols_ci95"][1],
                "ls_slope_ols_top50_ci975_lo": top50_ls["slope_ols_ci975"][0],
                "ls_slope_ols_top50_ci975_hi": top50_ls["slope_ols_ci975"][1],
                "ls_slope_noint_top50": top50_ls["slope_no_intercept"],
                "ls_slope_per_pert_top50": top50_ls["per_pert_slope_mean"],
                "ls_sd_ratio_top50": top50_ls["sd_ratio"],
                "ls_sd_ratio_top50_ci95_lo": top50_ls["sd_ratio_ci95"][0],
                "ls_sd_ratio_top50_ci95_hi": top50_ls["sd_ratio_ci95"][1],
                "ls_slope_obs_on_pred_top50": top50_ls["slope_obs_on_pred"],
                "ls_slope_obs_on_pred_top50_ci95_lo": top50_ls["slope_obs_on_pred_ci95"][0],
                "ls_slope_obs_on_pred_top50_ci95_hi": top50_ls["slope_obs_on_pred_ci95"][1],
                "ls_slope_deming_top50": top50_ls["slope_deming"],
                "ls_slope_deming_top50_ci95_lo": top50_ls["slope_deming_ci95"][0],
                "ls_slope_deming_top50_ci95_hi": top50_ls["slope_deming_ci95"][1],
                "ls_resid_corr_top50": top50_ls["resid_corr"],
                "ls_resid_corr_top50_ci95_lo": top50_ls["resid_corr_ci95"][0],
                "ls_resid_corr_top50_ci95_hi": top50_ls["resid_corr_ci95"][1],
                "ls_resid_corr_top50_ci975_lo": top50_ls["resid_corr_ci975"][0],
                "ls_resid_corr_top50_ci975_hi": top50_ls["resid_corr_ci975"][1],
                "ls_resid_corr_cal_top50": cal_r_top50,
                "ls_resid_corr_cal_top50_ci95_lo": cal_r_top50_ci_fixed[0],
                "ls_resid_corr_cal_top50_ci95_hi": cal_r_top50_ci_fixed[1],
                "ls_resid_corr_cal_top50_resamp_ci95_lo": cal_r_top50_ci_resamp[0],
                "ls_resid_corr_cal_top50_resamp_ci95_hi": cal_r_top50_ci_resamp[1],
                "ls_resid_corr_cal_top50_resamp_ci975_lo": cal_r_top50_ci975_resamp[0],
                "ls_resid_corr_cal_top50_resamp_ci975_hi": cal_r_top50_ci975_resamp[1],
                "ls_resid_corr_cal_top50_resamp_valid_frac": cal_r_top50_resamp_valid_frac,
                "ls_resid_corr_cal_top50_sq": cal_r_sq_top50,
                "ls_resid_corr_cal_top50_sq_ci95_lo": cal_r_sq_top50_ci_fixed[0],
                "ls_resid_corr_cal_top50_sq_ci95_hi": cal_r_sq_top50_ci_fixed[1],
                "ls_resid_corr_cal_top50_sq_resamp_ci95_lo": cal_r_sq_top50_ci_resamp[0],
                "ls_resid_corr_cal_top50_sq_resamp_ci95_hi": cal_r_sq_top50_ci_resamp[1],
                "ls_resid_corr_cal_top50_sq_resamp_ci975_lo": cal_r_sq_top50_ci975_resamp[0],
                "ls_resid_corr_cal_top50_sq_resamp_ci975_hi": cal_r_sq_top50_ci975_resamp[1],
                "ls_resid_corr_per_pert_top50": top50_ls["per_pert_corr_mean"],
                "pred_vs_obs_cosine_top50": top50_ls["pred_vs_obs_cosine_mean"],
                "ls_slope_ols_all": all_ls["slope_ols"],
                "ls_slope_ols_all_ci95_lo": all_ls["slope_ols_ci95"][0],
                "ls_slope_ols_all_ci95_hi": all_ls["slope_ols_ci95"][1],
                "ls_slope_noint_all": all_ls["slope_no_intercept"],
                "ls_slope_per_pert_all": all_ls["per_pert_slope_mean"],
                "ls_sd_ratio_all": all_ls["sd_ratio"],
                "ls_slope_obs_on_pred_all": all_ls["slope_obs_on_pred"],
                "ls_slope_deming_all": all_ls["slope_deming"],
                "ls_resid_corr_all": all_ls["resid_corr"],
                "ls_resid_corr_all_ci95_lo": all_ls["resid_corr_ci95"][0],
                "ls_resid_corr_all_ci95_hi": all_ls["resid_corr_ci95"][1],
                "ls_resid_corr_cal_all": cal_r_all,
                "ls_resid_corr_cal_all_ci95_lo": cal_r_all_ci_fixed[0],
                "ls_resid_corr_cal_all_ci95_hi": cal_r_all_ci_fixed[1],
                "ls_resid_corr_cal_all_resamp_ci95_lo": cal_r_all_ci_resamp[0],
                "ls_resid_corr_cal_all_resamp_ci95_hi": cal_r_all_ci_resamp[1],
                "ls_resid_corr_cal_all_resamp_ci975_lo": cal_r_all_ci975_resamp[0],
                "ls_resid_corr_cal_all_resamp_ci975_hi": cal_r_all_ci975_resamp[1],
                "ls_resid_corr_cal_all_resamp_valid_frac": cal_r_all_resamp_valid_frac,
                "ls_resid_corr_cal_all_sq": cal_r_sq_all,
                "ls_resid_corr_cal_all_sq_resamp_ci95_lo": cal_r_sq_all_ci_resamp[0],
                "ls_resid_corr_cal_all_sq_resamp_ci95_hi": cal_r_sq_all_ci_resamp[1],
            }
            for k, arr in merged.items():
                rec[k] = float(np.mean(arr))
            line_rows.append(rec)

        if pos_metrics_cent is not None and "Source_Grand_Mean" in row_metrics_by_method:
            neg = row_metrics_by_method["Source_Grand_Mean"]
            drf_calibration[target_line] = {
                "cent_disc_raw_drf": drf_score(float(np.mean(pos_metrics_cent["disc_raw"])), float(np.mean(neg["cent_disc_raw"]))),
                "cent_de50_r_drf": drf_score(float(np.mean(pos_metrics_cent["de50_r"])), float(np.mean(neg["cent_de50_r"]))),
                "ctrl_disc_raw_drf": drf_score(float(np.mean(pos_metrics_ctrl["disc_raw"])), float(np.mean(neg["ctrl_disc_raw"]))),
                "ctrl_de50_r_drf": drf_score(float(np.mean(pos_metrics_ctrl["de50_r"])), float(np.mean(neg["ctrl_de50_r"]))),
                "mse_drf": drf_score(float(np.mean(pos_metrics_cent["mse"])), float(np.mean(neg["mse"])), lower_is_better=True),
                "pos_ctrl_cent_disc_raw": float(np.mean(pos_metrics_cent["disc_raw"])),
                "pos_ctrl_cent_de50_r": float(np.mean(pos_metrics_cent["de50_r"])),
            }

        pooled_strata.extend([target_line] * len(row_idx))

        # Per-line paired bootstrap CIs for row-averaged metrics (including both A-T2 & A-T3 panel disc)
        row_met_keys = [
            "cent_disc_raw",
            "cent_disc_raw_panel",
            "cent_de50_r",
            "ctrl_disc_raw",
            "ctrl_disc_raw_panel",
            "ctrl_de50_r",
            "pds_cosine_panel",
            "mse",
        ]
        # Per-line paired bootstrap CIs for pooled line-specific metrics (plus joint-boot cal_r)
        ls_boot_keys = [
            ("resp_top50", "slope_ols", "ls_slope_ols_top50"),
            ("resp_top50", "slope_no_intercept", "ls_slope_noint_top50"),
            ("resp_top50", "per_pert_slope", "ls_slope_per_pert_top50"),
            ("resp_top50", "resid_corr", "ls_resid_corr_top50"),
            ("resp_top50", "resid_corr_cal", "ls_resid_corr_cal_top50"),
            ("all_nontarget", "slope_ols", "ls_slope_ols_all"),
            ("all_nontarget", "resid_corr", "ls_resid_corr_all"),
        ]
        line_comparisons = (
            list(comparisons)
            + [("Hybrid_ProjShared", f"NegCtrl_HybridProj_MismatchedCov({s_name})") for s_name in src_names]
            + [("Transfer_Mean_VarScaled", f"NegCtrl_TMVS_MismatchedVar({s_name})") for s_name in src_names]
        )
        for ma, mb in line_comparisons:
            if ma not in row_metrics_by_method or mb not in row_metrics_by_method:
                continue
            for met in row_met_keys:
                ci_diff = paired_row_bootstrap_ci(
                    row_metrics_by_method[ma][met],
                    row_metrics_by_method[mb][met],
                    n_boot=n_boot,
                    seed=seed,
                )
                ci_a = single_array_bootstrap_ci(row_metrics_by_method[ma][met], n_boot=n_boot, seed=seed)
                ci_b = single_array_bootstrap_ci(row_metrics_by_method[mb][met], n_boot=n_boot, seed=seed)
                boot_rows.append(
                    {
                        "target_line": target_line,
                        "metric": met,
                        "method_a": ma,
                        "method_b": mb,
                        "disjoint_ci95": cis_non_overlapping(
                            (ci_a["ci95_lo"], ci_a["ci95_hi"]),
                            (ci_b["ci95_lo"], ci_b["ci95_hi"]),
                        ),
                        "disjoint_ci975": cis_non_overlapping(
                            (ci_a["ci975_lo"], ci_a["ci975_hi"]),
                            (ci_b["ci975_lo"], ci_b["ci975_hi"]),
                        ),
                        "a_ci95_lo": ci_a["ci95_lo"],
                        "a_ci95_hi": ci_a["ci95_hi"],
                        "b_ci95_lo": ci_b["ci95_lo"],
                        "b_ci95_hi": ci_b["ci95_hi"],
                        "a_ci975_lo": ci_a["ci975_lo"],
                        "a_ci975_hi": ci_a["ci975_hi"],
                        "b_ci975_lo": ci_b["ci975_lo"],
                        "b_ci975_hi": ci_b["ci975_hi"],
                        **ci_diff,
                    }
                )
            for mask_mode, b_key, met_label in ls_boot_keys:
                ba = ls_by_method[ma][mask_mode][f"boot_{b_key}"]
                bb = ls_by_method[mb][mask_mode][f"boot_{b_key}"]
                pt_key = "per_pert_slope_mean" if b_key == "per_pert_slope" else b_key
                pt_a = float(ls_by_method[ma][mask_mode][pt_key])
                pt_b = float(ls_by_method[mb][mask_mode][pt_key])
                diff_b = ba - bb
                ci_d = (
                    np.nanpercentile(diff_b, [2.5, 97.5])
                    if np.any(np.isfinite(diff_b))
                    else np.array([np.nan, np.nan])
                )
                ci_d975 = (
                    np.nanpercentile(diff_b, [1.25, 98.75])
                    if np.any(np.isfinite(diff_b))
                    else np.array([np.nan, np.nan])
                )
                ci_a = (
                    np.nanpercentile(ba, [2.5, 97.5])
                    if np.any(np.isfinite(ba))
                    else np.array([np.nan, np.nan])
                )
                ci_a975 = (
                    np.nanpercentile(ba, [1.25, 98.75])
                    if np.any(np.isfinite(ba))
                    else np.array([np.nan, np.nan])
                )
                ci_b = (
                    np.nanpercentile(bb, [2.5, 97.5])
                    if np.any(np.isfinite(bb))
                    else np.array([np.nan, np.nan])
                )
                ci_b975 = (
                    np.nanpercentile(bb, [1.25, 98.75])
                    if np.any(np.isfinite(bb))
                    else np.array([np.nan, np.nan])
                )
                boot_se = float(np.nanstd(diff_b)) if np.any(np.isfinite(diff_b)) else float("nan")
                disjoint = (
                    cis_non_overlapping(
                        (float(ci_a[0]), float(ci_a[1])),
                        (float(ci_b[0]), float(ci_b[1])),
                    )
                    if np.all(np.isfinite(ci_a)) and np.all(np.isfinite(ci_b))
                    else False
                )
                disjoint975 = (
                    cis_non_overlapping(
                        (float(ci_a975[0]), float(ci_a975[1])),
                        (float(ci_b975[0]), float(ci_b975[1])),
                    )
                    if np.all(np.isfinite(ci_a975)) and np.all(np.isfinite(ci_b975))
                    else False
                )
                boot_rows.append(
                    {
                        "target_line": target_line,
                        "metric": met_label,
                        "method_a": ma,
                        "method_b": mb,
                        "disjoint_ci95": disjoint,
                        "disjoint_ci975": disjoint975,
                        "a_ci95_lo": float(ci_a[0]),
                        "a_ci95_hi": float(ci_a[1]),
                        "b_ci95_lo": float(ci_b[0]),
                        "b_ci95_hi": float(ci_b[1]),
                        "a_ci975_lo": float(ci_a975[0]),
                        "a_ci975_hi": float(ci_a975[1]),
                        "b_ci975_lo": float(ci_b975[0]),
                        "b_ci975_hi": float(ci_b975[1]),
                        "mean_a": pt_a,
                        "mean_b": pt_b,
                        "delta": pt_a - pt_b,
                        "ci95_lo": float(ci_d[0]),
                        "ci95_hi": float(ci_d[1]),
                        "ci975_lo": float(ci_d975[0]),
                        "ci975_hi": float(ci_d975[1]),
                        "boot_se": boot_se,
                        "ci_gt_zero": bool(np.isfinite(ci_d[0]) and ci_d[0] > 0.0),
                        "ci_lt_zero": bool(np.isfinite(ci_d[1]) and ci_d[1] < 0.0),
                        "ci975_gt_zero": bool(np.isfinite(ci_d975[0]) and ci_d975[0] > 0.0),
                        "ci975_lt_zero": bool(np.isfinite(ci_d975[1]) and ci_d975[1] < 0.0),
                    }
                )
        s5_expressed_top75[target_line] = _evaluate_s5_top75_var(
            preds=preds,
            true_deltas=true_deltas,
            transfer_ref=transfer_ref,
            target_ctrl_cov=t_ctrl_cov,
            q_gidx=q_gidx,
            ha_sub=ha_sub,
            hb_sub=hb_sub,
            obs_sub=obs_sub,
            tr_sub=tr_sub,
            hv_mask=hv_mask_eval,
            n_boot=n_boot,
            seed=seed,
        )
        del preds, row_metrics_by_method, ls_by_method
        gc.collect()

    df_lines = pd.DataFrame(line_rows)
    df_boot = pd.DataFrame(boot_rows)

    strat_summary: Dict[str, Dict[str, Any]] = {}
    for ma, mb in comparisons:
        if ma not in pooled_row_scores or mb not in pooled_row_scores:
            continue
        comp_key = f"{ma}_vs_{mb}"
        strat_summary[comp_key] = {}
        for met in [
            "cent_disc_raw",
            "cent_disc_raw_panel",
            "cent_de50_r",
            "ctrl_disc_raw",
            "ctrl_disc_raw_panel",
            "ctrl_de50_r",
            "pds_cosine_panel",
            "mse",
        ]:
            if len(pooled_row_scores[ma][met]) != len(pooled_strata) or len(pooled_row_scores[mb][met]) != len(pooled_strata):
                continue
            strat_summary[comp_key][met] = stratified_row_bootstrap_ci(
                np.asarray(pooled_row_scores[ma][met]),
                np.asarray(pooled_row_scores[mb][met]),
                pooled_strata,
                n_boot=n_boot,
                seed=seed,
            )

    # Cross-line cosine similarity diagnostic (Wrap-Up §2.6: predicted vs measured similarity)
    cross_line_cosine: Dict[str, Dict[str, float]] = {}
    for i in range(len(line_names)):
        for j in range(i + 1, len(line_names)):
            la, lb = line_names[i], line_names[j]
            da, db = saved_eval_by_line[la], saved_eval_by_line[lb]
            shared_t = sorted(set(da["target_to_row"].keys()) & set(db["target_to_row"].keys()))
            if not shared_t:
                continue
            G_dim = da["true_deltas"].shape[1]
            obs_cos_all, obs_cos_top50 = [], []
            pred_cos: Dict[str, Tuple[List[float], List[float]]] = {
                m: ([], []) for m in da["preds"].keys() if m in db["preds"]
            }
            for t_name in shared_t:
                ia = da["target_to_row"][t_name]
                ib = db["target_to_row"][t_name]
                ga = int(da["q_gidx"][ia])
                va_obs = da["true_deltas"][ia]
                vb_obs = db["true_deltas"][ib]
                mask = np.isfinite(va_obs) & np.isfinite(vb_obs)
                if 0 <= ga < G_dim:
                    mask[ga] = False
                idx_all = np.flatnonzero(mask)
                if idx_all.size < 2:
                    continue
                # Responsive union top-50 across the pair
                mag = np.maximum(np.abs(va_obs[idx_all]), np.abs(vb_obs[idx_all]))
                idx_top = idx_all[np.argsort(mag)[-min(50, idx_all.size) :]]

                def _cos(u: np.ndarray, v: np.ndarray) -> float:
                    den = float(np.linalg.norm(u) * np.linalg.norm(v))
                    return float(np.dot(u, v) / den) if den > 1e-12 else 0.0

                obs_cos_all.append(_cos(va_obs[idx_all], vb_obs[idx_all]))
                obs_cos_top50.append(_cos(va_obs[idx_top], vb_obs[idx_top]))
                for m_name in pred_cos:
                    pa = da["preds"][m_name][ia]
                    pb = db["preds"][m_name][ib]
                    pred_cos[m_name][0].append(_cos(pa[idx_all], pb[idx_all]))
                    pred_cos[m_name][1].append(_cos(pa[idx_top], pb[idx_top]))

            pair_dict: Dict[str, float] = {
                "n_shared_targets": float(len(obs_cos_all)),
                "measured_cross_line_cos_all": float(np.mean(obs_cos_all)),
                "measured_cross_line_cos_top50": float(np.mean(obs_cos_top50)),
            }
            for m_name, (c_all, c_top) in pred_cos.items():
                pair_dict[f"{m_name}_cos_all"] = float(np.mean(c_all))
                pair_dict[f"{m_name}_cos_top50"] = float(np.mean(c_top))
            cross_line_cosine[f"{la}__vs__{lb}"] = pair_dict

    summary = {
        "n_lines": len(line_names),
        "lines": line_names,
        "lambda_by_target_line": lam_by_target_line,
        "mismatched_cov_sources_by_line": mismatched_cov_sources_by_line,
        "ctrl_cov_reliability": ctrl_reliability,
        "pert_noise_ceilings": pert_noise_ceilings,
        "s5_expressed_top75_var": s5_expressed_top75,
        "drf_calibration": drf_calibration,
        "cross_line_cosine_diagnostics": cross_line_cosine,
        "macro_means_by_method": (
            df_lines.groupby("method").mean(numeric_only=True).to_dict(orient="index")
        ),
        "stratified_bootstrap_comparisons": strat_summary,
    }
    return df_lines, df_boot, summary


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Zero-shot line-specific perturbation evaluation.")
    p.add_argument("--phase", choices=["dev", "holdout"], required=True)
    p.add_argument("--repr", choices=["log1p", "zscore", "libnorm"], default="log1p")
    p.add_argument("--shrink", type=float, default=0.0)
    p.add_argument("--cov-max-cells", type=int, default=None)
    p.add_argument("--no-mask-self", action="store_true", default=False)
    p.add_argument("--tag-suffix", type=str, default="")
    p.add_argument("--all-targets", action="store_true", help="Evaluate all targets instead of shared.")
    p.add_argument("--only-reproducible", action="store_true", default=False)
    p.add_argument("--n-boot", type=int, default=4000)
    p.add_argument("--seed", type=int, default=12345)
    return p.parse_args()


def main() -> None:
    t0 = time.time()
    args = parse_args()
    avail = list_available_lines()
    if not avail:
        raise RuntimeError("No datasets found in p2/shared/lines/. Waiting for A-T1.")

    if args.phase == "dev":
        selected_lines = [ln for ln in avail if is_replogle_dev_line(ln)]
        if len(selected_lines) < 2:
            raise RuntimeError(
                f"Dev phase requires >= 2 Replogle lines (K562 & RPE1), found: {selected_lines}"
            )
        # Intersect strictly within Replogle dev lines so held-out lines never influence dev
        g_lists = []
        for ln in selected_lines:
            gf = SRC_DIR.parents[1] / "shared" / "lines" / ln / "genes.txt"
            g_lists.append([x.strip() for x in gf.read_text(encoding="utf-8").splitlines() if x.strip()])
        common_genes = sorted(set(g_lists[0]).intersection(*g_lists[1:]))
    else:
        selected_lines = avail
        check_prereg_gate(selected_lines)
        common_genes = load_common_genes()

    datasets: Dict[str, LineDataset] = {
        ln: load_line_dataset(ln, representation=args.repr, common_genes=common_genes)
        for ln in selected_lines
    }
    df_lines, df_boot, summary = evaluate_datasets(
        datasets=datasets,
        shrink=args.shrink,
        cov_max_cells=args.cov_max_cells,
        mask_self=not args.no_mask_self,
        shared_targets_only=not args.all_targets,
        only_reproducible=args.only_reproducible,
        n_boot=args.n_boot,
        seed=args.seed,
        require_official_spec=True,
    )
    summary["phase"] = args.phase
    summary["representation"] = args.repr
    summary["shrink"] = args.shrink
    summary["cov_max_cells"] = args.cov_max_cells
    summary["mask_self"] = not args.no_mask_self
    summary["only_reproducible"] = args.only_reproducible
    summary["shared_targets_only"] = not args.all_targets
    summary["elapsed_sec"] = round(time.time() - t0, 2)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    tag = f"{args.phase}_{args.repr}{args.tag_suffix}"
    df_lines.to_csv(RESULTS_DIR / f"line_metrics_{tag}.csv", index=False)
    df_boot.to_csv(RESULTS_DIR / f"bootstrap_ci_{tag}.csv", index=False)
    (RESULTS_DIR / f"summary_{tag}.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"Saved results to {RESULTS_DIR} (*_{tag}.csv / .json) in {summary['elapsed_sec']}s", flush=True)


if __name__ == "__main__":
    main()
