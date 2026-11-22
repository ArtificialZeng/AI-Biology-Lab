#!/usr/bin/env python3
"""
Round 12 Metrics Engine (`PREREG_E`).

Reuses frozen Round 10 and Round 11 residual/bootstrap primitives directly without
modifying any Round 10 or Round 11 file, and implements the Round 12 prospective additions:
  1. `DEFAULT_R12_SEED = 20261201`, Bonferroni 2-claim + IUT two-sided `97.5%` CIs.
  2. `compute_paired_fct_minus_all_mismatched_r12`: Joint perturbation bootstrap of
     paired `FCT - NegCtrl_FCT_MismatchedMean(s)` against EVERY source `s` in the pool
     (`H-E1b` individual IUT criteria across all sources).
  3. `compute_joint_sa8_split_minus_all_r12`: Joint perturbation cross-mask bootstrap of
     `r_orth_cal(pred, resp_top50_split) - r_orth_cal(pred, all_nontarget)` evaluated on
     the exact common half-valid perturbation subset (`rows_common`) using single-half
     disattenuation `sqrt(rho_AB)` for `resp_top50_split` and Spearman-Brown for
     `all_nontarget` (`H-E2b`).
  4. Canonical `N_5` Null-Cone semi-partial correlation suite (`compute_n5_suite_r12`,
     `compute_n5_perm100_r12`) for reporting-only structural baseline analysis with explicit
     marking of by-construction zeroed / cone-suppressed controls.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

ROUND12_SRC_DIR = Path(__file__).resolve().parent
REPO_ROOT = ROUND12_SRC_DIR.parents[1]
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
R11_SRC_DIR = REPO_ROOT / "round11" / "src"
for p in (str(R10_SRC_DIR), str(R11_SRC_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from fct_pipeline import FctSourceState, compute_source_transfer_and_fc_bar  # noqa: E402
from metrics_r10 import (  # noqa: E402
    EPS,
    _extract_masked_idx_per_pert,
    _moments_to_corr_and_slope,
    disattenuated_r,
)
from metrics_r11 import (  # noqa: E402
    _boot_disattenuated_orth_from_moments,
    _collect_mode_orth_moments,
    _prepare_split_direction_indices,
    compute_fct_perm100_null_r11,
    compute_fct_perm100_split_r11,
    compute_joint_cross_mask_e1_diff_r11,
    compute_joint_paired_r_orth_diff_r11,
    compute_residual_suite_r11,
    compute_resp_top50_split_suite_r11,
    compute_s7_cellcount_and_kd_r11,
    paired_row_bootstrap_ci_r11,
    single_half_ceiling,
    summarize_boot_r11,
)

DEFAULT_R12_SEED = 20261201
N5_CONE_NODES: Tuple[float, ...] = (0.05, 1.0, 20.0)
N5_EXACT_ZERO_METHODS = frozenset({"Transfer_Mean", "Transfer_Mean_GlobalScale", "Source_Grand_Mean"})
N5_CONE_SUPPRESSED_METHODS = frozenset(
    {"E1_ConstMed", "E2_ConstP10", "E2_ConstP90", "E6_LogFCBar", "E7_FCBarMinus1"}
)


def _eval_single_half_boot(
    pred_mom: Tuple[np.ndarray, ...],
    half_mom: Tuple[np.ndarray, ...],
    draws: np.ndarray,
) -> Tuple[float, np.ndarray]:
    """Evaluate point estimate and joint bootstrap draws under single-half ceiling `sqrt(rho_AB)`."""
    nv, sx, sy, sxx, syy, sxy = pred_mom
    _, pt_r, _, _, _ = _moments_to_corr_and_slope(
        float(nv.sum()), float(sx.sum()), float(sy.sum()),
        float(sxx.sum()), float(syy.sum()), float(sxy.sum()),
    )
    _, b_r, _, _, _ = _moments_to_corr_and_slope(
        nv[draws].sum(axis=1), sx[draws].sum(axis=1), sy[draws].sum(axis=1),
        sxx[draws].sum(axis=1), syy[draws].sum(axis=1), sxy[draws].sum(axis=1),
    )

    hnv, hsx, hsy, hsxx, hsyy, hsxy = half_mom
    _, pt_rho, _, _, _ = _moments_to_corr_and_slope(
        float(hnv.sum()), float(hsx.sum()), float(hsy.sum()),
        float(hsxx.sum()), float(hsyy.sum()), float(hsxy.sum()),
    )
    pt_cal = disattenuated_r(float(pt_r), single_half_ceiling(float(pt_rho)))
    _, b_rho, _, _, _ = _moments_to_corr_and_slope(
        hnv[draws].sum(axis=1), hsx[draws].sum(axis=1), hsy[draws].sum(axis=1),
        hsxx[draws].sum(axis=1), hsyy[draws].sum(axis=1), hsxy[draws].sum(axis=1),
    )
    b_cal = np.array(
        [disattenuated_r(float(r), single_half_ceiling(float(rh))) for r, rh in zip(b_r, b_rho)],
        dtype=np.float64,
    )
    return float(pt_cal), b_cal


def compute_joint_sa8_split_minus_all_r12(
    delta_obs: np.ndarray,
    delta_pred: np.ndarray,
    delta_transfer: np.ndarray,
    target_gene_idx: Sequence[int],
    delta_half_a: Optional[np.ndarray],
    delta_half_b: Optional[np.ndarray],
    half_valid: Optional[np.ndarray] = None,
    gene_keep_mask: Optional[np.ndarray] = None,
    top_k: int = 50,
    n_boot: int = 4000,
    seed: int = DEFAULT_R12_SEED,
) -> Dict[str, Any]:
    """
    `PREREG_E H-E2b`: Joint perturbation bootstrap of the cross-mask contrast:
        `r_orth_cal(delta_pred, resp_top50_split) - r_orth_cal(delta_pred, all_nontarget)`
    restricted strictly to the common perturbations valid across both halves (`rows_common`)
    and `all_nontarget`, using identical bootstrap draw indices `draws`.
    """
    empty = summarize_boot_r11(float("nan"), np.array([float("nan")]), 0, n_boot)
    if delta_half_a is None or delta_half_b is None:
        return {**empty, "n_common_perts": 0, "split_summary": empty, "all_common_summary": empty}

    d_obs = np.asarray(delta_obs, dtype=np.float64)
    p_mat = np.asarray(delta_pred, dtype=np.float64)
    tm = np.asarray(delta_transfer, dtype=np.float64)
    ha = np.asarray(delta_half_a, dtype=np.float64)
    hb = np.asarray(delta_half_b, dtype=np.float64)

    rows_ab, idx_ab = _prepare_split_direction_indices(
        ha, hb, ha, tm, target_gene_idx, half_valid,
        extra_valid=p_mat, gene_keep_mask=gene_keep_mask, top_k=top_k,
    )
    rows_ba, idx_ba = _prepare_split_direction_indices(
        hb, ha, hb, tm, target_gene_idx, half_valid,
        extra_valid=p_mat, gene_keep_mask=gene_keep_mask, top_k=top_k,
    )
    rows_all, idx_all = _extract_masked_idx_per_pert(
        d_obs, target_gene_idx, gene_keep_mask=gene_keep_mask, mode="all_nontarget", top_k=top_k,
    )
    map_ab = dict(zip(rows_ab, idx_ab))
    map_ba = dict(zip(rows_ba, idx_ba))
    map_all = dict(zip(rows_all, idx_all))

    rows_common = [p for p in rows_ab if p in map_ba and p in map_all]
    m = len(rows_common)
    if m == 0:
        return {**empty, "n_common_perts": 0, "split_summary": empty, "all_common_summary": empty}

    idx_ab_c = [map_ab[p] for p in rows_common]
    idx_ba_c = [map_ba[p] for p in rows_common]
    idx_all_c = [map_all[p] for p in rows_common]

    pred_res = p_mat - tm
    r_obs_full = d_obs - tm
    ha_res = ha - tm
    hb_res = hb - tm

    [pm_ab], hm_ab = _collect_mode_orth_moments(
        hb_res, [pred_res], tm, rows_common, idx_ab_c, ha_res, hb_res, half_valid
    )
    [pm_ba], hm_ba = _collect_mode_orth_moments(
        ha_res, [pred_res], tm, rows_common, idx_ba_c, hb_res, ha_res, half_valid
    )
    [pm_all], hm_all = _collect_mode_orth_moments(
        r_obs_full, [pred_res], tm, rows_common, idx_all_c, ha_res, hb_res, half_valid
    )
    if hm_ab is None or hm_ba is None or hm_all is None:
        return {**empty, "n_common_perts": m, "split_summary": empty, "all_common_summary": empty}

    draws = np.random.default_rng(seed).integers(0, m, size=(n_boot, m))
    pt_cal_ab, b_cal_ab = _eval_single_half_boot(pm_ab, hm_ab, draws)
    pt_cal_ba, b_cal_ba = _eval_single_half_boot(pm_ba, hm_ba, draws)
    pt_split = 0.5 * (pt_cal_ab + pt_cal_ba)
    b_split = 0.5 * (b_cal_ab + b_cal_ba)

    pt_all, b_all, _, _ = _boot_disattenuated_orth_from_moments(pm_all, hm_all, draws)
    diff_pt = float(pt_split - pt_all)
    b_diff = b_split - b_all
    diff_sum = summarize_boot_r11(diff_pt, b_diff, m, n_boot)
    return {
        **diff_sum,
        "n_common_perts": int(m),
        "split_summary": summarize_boot_r11(pt_split, b_split, m, n_boot),
        "all_common_summary": summarize_boot_r11(pt_all, b_all, m, n_boot),
    }


def compute_paired_fct_minus_all_mismatched_r12(
    delta_obs: np.ndarray,
    preds: Dict[str, np.ndarray],
    delta_transfer: np.ndarray,
    target_gene_idx: Sequence[int],
    delta_half_a: Optional[np.ndarray],
    delta_half_b: Optional[np.ndarray],
    half_valid: Optional[np.ndarray] = None,
    gene_keep_mask: Optional[np.ndarray] = None,
    n_boot: int = 4000,
    seed: int = DEFAULT_R12_SEED,
) -> Dict[str, Dict[str, Any]]:
    """
    `PREREG_E H-E1b`: Joint perturbation bootstrap of `FCT - NegCtrl_FCT_MismatchedMean(s)`
    for EVERY source line `s` present in `preds` under primary `all_nontarget` mask.
    """
    mis_keys = sorted(k for k in preds if k.startswith("NegCtrl_FCT_MismatchedMean("))
    out: Dict[str, Dict[str, Any]] = {}
    for mk in mis_keys:
        src_name = mk[len("NegCtrl_FCT_MismatchedMean(") : -1]
        out[src_name] = compute_joint_paired_r_orth_diff_r11(
            delta_obs=delta_obs,
            delta_pred_a=preds["FCT"],
            delta_pred_b=preds[mk],
            delta_transfer=delta_transfer,
            target_gene_idx=target_gene_idx,
            delta_half_a=delta_half_a,
            delta_half_b=delta_half_b,
            half_valid=half_valid,
            gene_keep_mask=gene_keep_mask,
            mode="all_nontarget",
            top_k=50,
            n_boot=n_boot,
            seed=seed,
        )
    return out


# ==============================================================================
# Canonical `N_5` Null-Cone Basis & Semi-Partial Residual Correlation Suite
# ==============================================================================


def g_const_cone(fc: np.ndarray, c: float) -> np.ndarray:
    """Normalized constant-baseline response basis function `log1p(c * FC) - log1p(c)`."""
    return np.log1p(float(c) * np.asarray(fc, dtype=np.float64)) - np.log1p(float(c))


def build_n5_basis(
    t_mean_p: np.ndarray,
    grand_mean_p: np.ndarray,
    fc_bar_p: np.ndarray,
    rcond: float = 1e-12,
) -> np.ndarray:
    """
    Orthonormal SVD basis `Q_p` (`m_p x r_p`) spanning:
        `N_5(p) = span{1, t_mean_p, grand_mean_p, g_{0.05}(fc_bar_p), g_{1.0}(fc_bar_p), g_{20.0}(fc_bar_p)}`.
    """
    m = t_mean_p.shape[0]
    ones = np.full(m, 1.0 / np.sqrt(float(m)), dtype=np.float64)
    cols_raw = [
        t_mean_p.astype(np.float64),
        grand_mean_p.astype(np.float64),
        g_const_cone(fc_bar_p, N5_CONE_NODES[0]),
        g_const_cone(fc_bar_p, N5_CONE_NODES[1]),
        g_const_cone(fc_bar_p, N5_CONE_NODES[2]),
    ]
    orth_cols = []
    for v in cols_raw:
        vp = v - float(np.dot(ones, v)) * ones
        nrm = float(np.linalg.norm(vp))
        if nrm > EPS:
            orth_cols.append(vp / nrm)
    if not orth_cols:
        return ones[:, None]
    mat = np.column_stack(orth_cols)
    u, sv, _ = np.linalg.svd(mat, full_matrices=False)
    keep = sv > max(EPS, rcond * float(sv[0]))
    if not np.any(keep):
        return ones[:, None]
    return np.column_stack([ones, u[:, keep]])


def prepare_n5_per_pert_data(
    delta_obs: np.ndarray,
    t_mean: np.ndarray,
    grand_mean: np.ndarray,
    fc_bar: np.ndarray,
    target_gene_idx: Sequence[int],
    half_a_obs: Optional[np.ndarray],
    half_b_obs: Optional[np.ndarray],
    half_valid: Optional[np.ndarray],
    gene_keep_mask: Optional[np.ndarray] = None,
    mode: str = "all_nontarget",
    top_k: int = 50,
    precomputed_items: Optional[Sequence[Tuple[int, np.ndarray]]] = None,
    eval_delta_obs: Optional[np.ndarray] = None,
) -> Tuple[List[Tuple[int, np.ndarray, np.ndarray]], np.ndarray, np.ndarray, np.ndarray, np.ndarray, Dict[str, Any]]:
    """Precompute per-perturbation `N_5` bases `Q_p` and observed/split-half orthogonalized moments."""
    obs_eval = delta_obs if eval_delta_obs is None else eval_delta_obs
    n_t, n_g = obs_eval.shape
    per_p: List[Tuple[int, np.ndarray, np.ndarray]] = []
    ro_ss: List[float] = []
    haa_ss: List[float] = []
    hbb_ss: List[float] = []
    hab_ss: List[float] = []

    has_halves = half_a_obs is not None and half_b_obs is not None
    ranks: List[int] = []

    if precomputed_items is not None:
        item_iter: List[Tuple[int, np.ndarray]] = [(int(p), np.asarray(idx, dtype=np.int64)) for p, idx in precomputed_items if len(idx) > 0]
    else:
        item_iter = []
        for p in range(n_t):
            valid = np.isfinite(obs_eval[p]) & np.isfinite(t_mean[p]) & np.isfinite(fc_bar[p]) & np.isfinite(grand_mean)
            if gene_keep_mask is not None:
                valid = valid & gene_keep_mask
            tc = int(target_gene_idx[p])
            if 0 <= tc < n_g:
                valid[tc] = False
            idx = np.flatnonzero(valid)
            if idx.size == 0:
                continue
            if mode == "resp_top50" and idx.size > top_k:
                idx = idx[np.argsort(np.abs(obs_eval[p, idx]), kind="mergesort")[-top_k:]]
            item_iter.append((p, idx))

    for p, idx in item_iter:
        tm_p = t_mean[p, idx].astype(np.float64)
        gm_p = grand_mean[idx].astype(np.float64)
        fcb_p = fc_bar[p, idx].astype(np.float64)
        q_p = build_n5_basis(tm_p, gm_p, fcb_p)
        ranks.append(q_p.shape[1])

        r_obs = obs_eval[p, idx].astype(np.float64) - tm_p
        r_obs_orth = r_obs - q_p @ (q_p.T @ r_obs)
        per_p.append((p, idx, q_p))
        ro_ss.append(float(np.dot(r_obs_orth, r_obs_orth)))

        if has_halves and (half_valid is None or bool(half_valid[p])):
            ha = half_a_obs[p, idx].astype(np.float64)
            hb = half_b_obs[p, idx].astype(np.float64)
            if np.all(np.isfinite(ha) & np.isfinite(hb)):
                ra = ha - tm_p
                rb = hb - tm_p
                ra_o = ra - q_p @ (q_p.T @ ra)
                rb_o = rb - q_p @ (q_p.T @ rb)
                haa_ss.append(float(np.dot(ra_o, ra_o)))
                hbb_ss.append(float(np.dot(rb_o, rb_o)))
                hab_ss.append(float(np.dot(ra_o, rb_o)))
                continue
        haa_ss.append(0.0)
        hbb_ss.append(0.0)
        hab_ss.append(0.0)

    meta = {
        "n_eval_perts": len(per_p),
        "rank_min": min(ranks) if ranks else 0,
        "rank_max": max(ranks) if ranks else 0,
        "rank_mean": float(np.mean(ranks)) if ranks else 0.0,
    }
    return (
        per_p,
        np.asarray(ro_ss, dtype=np.float64),
        np.asarray(haa_ss, dtype=np.float64),
        np.asarray(hbb_ss, dtype=np.float64),
        np.asarray(hab_ss, dtype=np.float64),
        meta,
    )


def _pred_n5_moments(
    per_p: Sequence[Tuple[int, np.ndarray, np.ndarray]],
    delta_obs: np.ndarray,
    pred: np.ndarray,
    t_mean: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-perturbation `<r_obs_orth, r_pred_orth>`, `||r_pred_orth||^2`, `||r_pred_cent||^2`."""
    m = len(per_p)
    xy = np.empty(m, dtype=np.float64)
    yy_o = np.empty(m, dtype=np.float64)
    yy_c = np.empty(m, dtype=np.float64)
    for k, (p, idx, q_p) in enumerate(per_p):
        tm_p = t_mean[p, idx].astype(np.float64)
        r_obs = delta_obs[p, idx].astype(np.float64) - tm_p
        r_obs_o = r_obs - q_p @ (q_p.T @ r_obs)
        r_pr = pred[p, idx].astype(np.float64) - tm_p
        r_pr_c = r_pr - float(np.mean(r_pr))
        r_pr_o = r_pr - q_p @ (q_p.T @ r_pr)
        xy[k] = float(np.dot(r_obs_o, r_pr_o))
        yy_o[k] = float(np.dot(r_pr_o, r_pr_o))
        yy_c[k] = float(np.dot(r_pr_c, r_pr_c))
    return xy, yy_o, yy_c


def _calc_n5_vectorized(
    s_xy: Any, s_xx: Any, s_yyc: Any, s_yyo: Any,
    s_ha: Any, s_hb: Any, s_hab: Any,
    single_half: bool = False,
) -> Tuple[Any, Any, Any, Any, Any, Any]:
    eta = np.where(s_yyc > EPS, np.sqrt(np.clip(s_yyo / np.maximum(s_yyc, EPS), 0.0, 1.0)), 0.0)
    denom_sp = np.sqrt(np.maximum(s_xx, 0.0) * np.maximum(s_yyc, 0.0))
    r_sp_raw = np.where((denom_sp > EPS) & (s_yyo > 1e-24), np.clip(s_xy / np.maximum(denom_sp, EPS), -1.0, 1.0), 0.0)
    denom_orth = np.sqrt(np.maximum(s_xx, 0.0) * np.maximum(s_yyo, 0.0))
    r_orth_raw = np.where((denom_orth > EPS) & (s_yyo > 1e-24), np.clip(s_xy / np.maximum(denom_orth, EPS), -1.0, 1.0), 0.0)
    hden = np.sqrt(np.maximum(s_ha, 0.0) * np.maximum(s_hb, 0.0))
    rho = np.where(hden > EPS, np.clip(s_hab / np.maximum(hden, EPS), -1.0, 1.0), 0.0)
    rho_pos = np.maximum(rho, 0.0)
    ceil_pos = np.sqrt(rho_pos) if single_half else np.sqrt((2.0 * rho_pos) / (1.0 + rho_pos))
    ceil = np.where(rho > 0.0, ceil_pos, 0.0)
    r_sp_cal = np.where(ceil > EPS, r_sp_raw / np.maximum(ceil, EPS), float("nan"))
    return r_sp_cal, r_sp_raw, r_orth_raw, eta, rho, ceil


def _summarize_paired_mis_n5(
    preds: Dict[str, np.ndarray],
    boot_draws_cal: Dict[str, Tuple[float, np.ndarray]],
    m: int,
    n_boot: int,
) -> Dict[str, Dict[str, Any]]:
    paired_mis: Dict[str, Dict[str, Any]] = {}
    if "FCT" in boot_draws_cal and m > 0:
        fct_pt, fct_b = boot_draws_cal["FCT"]
        for mk in sorted(k for k in preds if k.startswith("NegCtrl_FCT_MismatchedMean(")):
            src_name = mk[len("NegCtrl_FCT_MismatchedMean(") : -1]
            mis_pt, mis_b = boot_draws_cal[mk]
            paired_mis[src_name] = summarize_boot_r11(float(fct_pt - mis_pt), fct_b - mis_b, m, n_boot)
    return paired_mis


def _n5_construction_role(method_name: str) -> str:
    if method_name in N5_EXACT_ZERO_METHODS:
        return "ZERO_BY_EXACT_CONSTRUCTION"
    if method_name in N5_CONE_SUPPRESSED_METHODS:
        return "SUPPRESSED_BY_CONE_CONSTRUCTION"
    return "INFORMATIVE"


def compute_n5_suite_r12(
    delta_obs: np.ndarray,
    preds: Dict[str, np.ndarray],
    delta_transfer: np.ndarray,
    grand_mean: np.ndarray,
    fc_bar: np.ndarray,
    target_gene_idx: Sequence[int],
    delta_half_a: Optional[np.ndarray],
    delta_half_b: Optional[np.ndarray],
    half_valid: Optional[np.ndarray] = None,
    gene_keep_mask: Optional[np.ndarray] = None,
    n_boot: int = 4000,
    seed: int = DEFAULT_R12_SEED,
) -> Dict[str, Any]:
    """
    Evaluate `N_5` Null-Cone semi-partial `r_cone_sp` (disattenuated by `N_5` split-half
    ceiling) across all predictors and paired `FCT - Mis_s` contrasts using identical
    perturbation bootstrap draws (`PREREG_E` reporting-only analysis).
    """
    per_p, xx_o, haa_o, hbb_o, hab_o, meta = prepare_n5_per_pert_data(
        delta_obs=delta_obs,
        t_mean=delta_transfer,
        grand_mean=grand_mean,
        fc_bar=fc_bar,
        target_gene_idx=target_gene_idx,
        half_a_obs=delta_half_a,
        half_b_obs=delta_half_b,
        half_valid=half_valid,
        gene_keep_mask=gene_keep_mask,
        mode="all_nontarget",
    )
    m = len(per_p)
    draws = np.random.default_rng(seed).integers(0, m, size=(n_boot, m)) if m > 0 else None

    methods_out: Dict[str, Dict[str, Any]] = {}
    boot_draws_cal: Dict[str, Tuple[float, np.ndarray]] = {}

    for m_name, p_mat in preds.items():
        if m == 0:
            empty = summarize_boot_r11(float("nan"), np.array([float("nan")]), 0, n_boot)
            methods_out[m_name] = {
                "r_cone_sp": empty,
                "r_cone_sp_raw": float("nan"),
                "r_cone_orth_raw": float("nan"),
                "eta": float("nan"),
                "r2_n5": float("nan"),
                "rho_half_ab_cone": float("nan"),
                "ceiling_cone": float("nan"),
                "construction_status": _n5_construction_role(m_name),
            }
            continue

        pt_cal, b_cal, pt_raw, pt_orth_raw, pt_eta, pt_rho, pt_ceil, pt_r2 = _eval_n5_direction(
            per_p, delta_obs, p_mat, delta_transfer, xx_o, haa_o, hbb_o, hab_o, draws, single_half=False
        )
        boot_draws_cal[m_name] = (pt_cal, b_cal)
        methods_out[m_name] = {
            "r_cone_sp": summarize_boot_r11(pt_cal, b_cal, m, n_boot),
            "r_cone_sp_raw": pt_raw,
            "r_cone_orth_raw": pt_orth_raw,
            "eta": pt_eta,
            "r2_n5": pt_r2,
            "rho_half_ab_cone": pt_rho,
            "ceiling_cone": pt_ceil,
            "construction_status": _n5_construction_role(m_name),
        }

    paired_mis = _summarize_paired_mis_n5(preds, boot_draws_cal, m, n_boot)
    return {
        "meta": meta,
        "methods": methods_out,
        "paired_fct_minus_mis": paired_mis,
    }


def _eval_n5_direction(
    per_p: Sequence[Tuple[int, np.ndarray, np.ndarray]],
    eval_obs: np.ndarray,
    p_mat: np.ndarray,
    delta_transfer: np.ndarray,
    xx_o: np.ndarray,
    haa_o: np.ndarray,
    hbb_o: np.ndarray,
    hab_o: np.ndarray,
    draws: Optional[np.ndarray],
    single_half: bool = False,
) -> Tuple[float, np.ndarray, float, float, float, float, float, float]:
    """Evaluate point estimate and bootstrap draws for one `N_5` direction."""
    s_xx = float(np.sum(xx_o))
    s_ha = float(np.sum(haa_o))
    s_hb = float(np.sum(hbb_o))
    s_hab = float(np.sum(hab_o))
    xy, yy_o, yy_c = _pred_n5_moments(per_p, eval_obs, p_mat, delta_transfer)
    s_yyc = float(np.sum(yy_c))
    pt_cal, pt_raw, pt_orth_raw, pt_eta, pt_rho, pt_ceil = _calc_n5_vectorized(
        float(np.sum(xy)), s_xx, s_yyc, float(np.sum(yy_o)), s_ha, s_hb, s_hab,
        single_half=single_half,
    )
    pt_r2 = float(max(0.0, 1.0 - float(pt_eta) ** 2)) if s_yyc > EPS else 1.0
    if draws is not None:
        b_cal, _, _, _, _, _ = _calc_n5_vectorized(
            np.sum(xy[draws], axis=1), np.sum(xx_o[draws], axis=1),
            np.sum(yy_c[draws], axis=1), np.sum(yy_o[draws], axis=1),
            np.sum(haa_o[draws], axis=1), np.sum(hbb_o[draws], axis=1), np.sum(hab_o[draws], axis=1),
            single_half=single_half,
        )
    else:
        b_cal = np.array([float(pt_cal)], dtype=np.float64)
    return (
        float(pt_cal),
        np.asarray(b_cal, dtype=np.float64),
        float(pt_raw),
        float(pt_orth_raw),
        float(pt_eta),
        float(pt_rho),
        float(pt_ceil),
        float(pt_r2),
    )


def compute_n5_resp_top50_split_r12(
    delta_obs: np.ndarray,
    preds: Dict[str, np.ndarray],
    delta_transfer: np.ndarray,
    grand_mean: np.ndarray,
    fc_bar: np.ndarray,
    target_gene_idx: Sequence[int],
    delta_half_a: Optional[np.ndarray],
    delta_half_b: Optional[np.ndarray],
    half_valid: Optional[np.ndarray] = None,
    gene_keep_mask: Optional[np.ndarray] = None,
    top_k: int = 50,
    n_boot: int = 4000,
    seed: int = DEFAULT_R12_SEED,
) -> Dict[str, Any]:
    """
    Evaluate `PREREG_E §5.3` `N_5` Null-Cone semi-partial `r_cone_sp` under split-half
    out-of-sample `resp_top50_split` (`A->B` and `B->A`, disattenuated by single-half
    cone ceiling `sqrt(rho_AB_cone)`) across all predictors & paired `FCT - Mis_s`.
    """
    if delta_half_a is None or delta_half_b is None:
        return {"meta": {"n_eval_perts": 0, "status": "missing_split_halves"}, "methods": {}, "paired_fct_minus_mis": {}}

    ha = np.asarray(delta_half_a, dtype=np.float64)
    hb = np.asarray(delta_half_b, dtype=np.float64)
    tm = np.asarray(delta_transfer, dtype=np.float64)
    fcb = np.asarray(fc_bar, dtype=np.float64)
    gm = np.asarray(grand_mean, dtype=np.float64)
    gmask = np.isfinite(gm) if gene_keep_mask is None else (np.asarray(gene_keep_mask, dtype=bool) & np.isfinite(gm))

    rows_ab, idx_ab = _prepare_split_direction_indices(
        ha, hb, ha, tm, target_gene_idx, half_valid,
        extra_valid=fcb, gene_keep_mask=gmask, top_k=top_k,
    )
    rows_ba, idx_ba = _prepare_split_direction_indices(
        hb, ha, hb, tm, target_gene_idx, half_valid,
        extra_valid=fcb, gene_keep_mask=gmask, top_k=top_k,
    )
    map_ab = dict(zip(rows_ab, idx_ab))
    map_ba = dict(zip(rows_ba, idx_ba))
    rows_common = [p for p in rows_ab if p in map_ba]
    m = len(rows_common)
    if m == 0:
        return {"meta": {"n_eval_perts": 0, "status": "no_valid_split_perturbations"}, "methods": {}, "paired_fct_minus_mis": {}}

    items_ab = [(p, map_ab[p]) for p in rows_common]
    items_ba = [(p, map_ba[p]) for p in rows_common]

    per_ab, xx_ab, haa_ab, hbb_ab, hab_ab, meta_ab = prepare_n5_per_pert_data(
        delta_obs=delta_obs, t_mean=tm, grand_mean=gm, fc_bar=fcb,
        target_gene_idx=target_gene_idx, half_a_obs=ha, half_b_obs=hb,
        half_valid=half_valid, precomputed_items=items_ab, eval_delta_obs=hb,
    )
    per_ba, xx_ba, haa_ba, hbb_ba, hab_ba, meta_ba = prepare_n5_per_pert_data(
        delta_obs=delta_obs, t_mean=tm, grand_mean=gm, fc_bar=fcb,
        target_gene_idx=target_gene_idx, half_a_obs=ha, half_b_obs=hb,
        half_valid=half_valid, precomputed_items=items_ba, eval_delta_obs=ha,
    )

    draws = np.random.default_rng(seed).integers(0, m, size=(n_boot, m))
    methods_out: Dict[str, Dict[str, Any]] = {}
    boot_draws_cal: Dict[str, Tuple[float, np.ndarray]] = {}

    for m_name, p_mat in preds.items():
        ca, ba_arr, ra, oa, ea, rho_a, ceil_a, r2a = _eval_n5_direction(
            per_ab, hb, p_mat, tm, xx_ab, haa_ab, hbb_ab, hab_ab, draws, single_half=True
        )
        cb, bb_arr, rb, ob, eb, rho_b, ceil_b, r2b = _eval_n5_direction(
            per_ba, ha, p_mat, tm, xx_ba, haa_ba, hbb_ba, hab_ba, draws, single_half=True
        )
        pt_cal = 0.5 * (ca + cb)
        b_cal = 0.5 * (ba_arr + bb_arr)
        boot_draws_cal[m_name] = (float(pt_cal), np.asarray(b_cal, dtype=np.float64))
        methods_out[m_name] = {
            "r_cone_sp": summarize_boot_r11(float(pt_cal), np.asarray(b_cal, dtype=np.float64), m, n_boot),
            "r_cone_sp_raw": float(0.5 * (ra + rb)),
            "r_cone_orth_raw": float(0.5 * (oa + ob)),
            "eta": float(0.5 * (ea + eb)),
            "r2_n5": float(0.5 * (r2a + r2b)),
            "rho_half_ab_cone_a": float(rho_a),
            "rho_half_ab_cone_b": float(rho_b),
            "ceiling_cone_single_a": float(ceil_a),
            "ceiling_cone_single_b": float(ceil_b),
            "construction_status": _n5_construction_role(m_name),
        }

    return {
        "meta": {"n_eval_perts": int(m), "mode": "N5_resp_top50_split", "ab": meta_ab, "ba": meta_ba},
        "methods": methods_out,
        "paired_fct_minus_mis": _summarize_paired_mis_n5(preds, boot_draws_cal, m, n_boot),
    }


def compute_n5_perm100_r12(
    state: FctSourceState,
    target_ctrl_mean_libnorm: np.ndarray,
    query_target_names: Sequence[str],
    delta_obs: np.ndarray,
    delta_transfer: np.ndarray,
    grand_mean: np.ndarray,
    target_gene_idx: Sequence[int],
    delta_half_a: Optional[np.ndarray],
    delta_half_b: Optional[np.ndarray],
    half_valid: Optional[np.ndarray] = None,
    fc_bundle: Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]] = None,
    gene_keep_mask: Optional[np.ndarray] = None,
    n_perm: int = 100,
) -> Dict[str, Any]:
    """
    Evaluate single-draw 100-permutation null (`k = 0..n_perm-1`) under canonical `N_5`.
    (Single-draw empirical null preserves `eta ~ 0.63-0.71` non-cone residual energy.)
    """
    mu_t = np.maximum(np.asarray(target_ctrl_mean_libnorm, dtype=np.float64), 0.0)
    G = mu_t.shape[0]
    if fc_bundle is None:
        fc_bundle = compute_source_transfer_and_fc_bar(state, query_target_names)
    pred_transfer, fc_bar, has_hits = fc_bundle

    per_p, xx_o, haa_o, hbb_o, hab_o, meta = prepare_n5_per_pert_data(
        delta_obs=delta_obs,
        t_mean=delta_transfer,
        grand_mean=grand_mean,
        fc_bar=fc_bar,
        target_gene_idx=target_gene_idx,
        half_a_obs=delta_half_a,
        half_b_obs=delta_half_b,
        half_valid=half_valid,
        gene_keep_mask=gene_keep_mask,
        mode="all_nontarget",
    )
    if not per_p:
        return {"n_perm": int(n_perm), "mode": "N5_all_nontarget", "r_cone_sp": {}, "med_eta": float("nan")}

    s_xx = float(np.sum(xx_o))
    s_ha = float(np.sum(haa_o))
    s_hb = float(np.sum(hbb_o))
    s_hab = float(np.sum(hab_o))

    r_sp_cal_arr = np.empty(n_perm, dtype=np.float64)
    r_sp_raw_arr = np.empty(n_perm, dtype=np.float64)
    eta_arr = np.empty(n_perm, dtype=np.float64)

    for k in range(n_perm):
        perm = np.random.default_rng(k).permutation(G)
        mu_p = mu_t[perm]
        s_xy = 0.0
        s_yyo = 0.0
        s_yyc = 0.0
        for p, idx, q_p in per_p:
            tm_p = delta_transfer[p, idx].astype(np.float64)
            r_obs = delta_obs[p, idx].astype(np.float64) - tm_p
            r_obs_o = r_obs - q_p @ (q_p.T @ r_obs)
            if bool(has_hits[p]):
                m_sub = mu_p[idx]
                pred_v = np.log1p(m_sub * fc_bar[p, idx]) - np.log1p(m_sub)
            else:
                pred_v = pred_transfer[p, idx]
            r_pr = pred_v - tm_p
            r_pr_c = r_pr - float(np.mean(r_pr))
            r_pr_o = r_pr - q_p @ (q_p.T @ r_pr)
            s_xy += float(np.dot(r_obs_o, r_pr_o))
            s_yyo += float(np.dot(r_pr_o, r_pr_o))
            s_yyc += float(np.dot(r_pr_c, r_pr_c))
        c_v, raw_v, _, eta_v, _, _ = _calc_n5_vectorized(s_xy, s_xx, s_yyc, s_yyo, s_ha, s_hb, s_hab)
        r_sp_cal_arr[k] = float(c_v)
        r_sp_raw_arr[k] = float(raw_v)
        eta_arr[k] = float(eta_v)

    fin = r_sp_cal_arr[np.isfinite(r_sp_cal_arr)]
    if fin.size > 0:
        p025, p25, p50, p975, p9875 = np.percentile(fin, [1.25, 2.5, 50.0, 97.5, 98.75])
        pcts = {
            "min": float(np.min(fin)),
            "p1_25": float(p025),
            "p2_5": float(p25),
            "median": float(p50),
            "median_abs": float(np.median(np.abs(fin))),
            "p97_5": float(p975),
            "p98_75": float(p9875),
            "max": float(np.max(fin)),
            "mean": float(np.mean(fin)),
        }
    else:
        pcts = {
            k: float("nan")
            for k in ("min", "p1_25", "p2_5", "median", "median_abs", "p97_5", "p98_75", "max", "mean")
        }

    return {
        "n_perm": int(n_perm),
        "mode": "N5_all_nontarget",
        "meta": meta,
        "r_cone_sp": pcts,
        "perm0_r_cone_sp": float(r_sp_cal_arr[0]),
        "perm0_r_cone_sp_raw": float(r_sp_raw_arr[0]),
        "med_eta": float(np.nanmedian(eta_arr)),
    }


__all__ = [
    "DEFAULT_R12_SEED",
    "EPS",
    "N5_CONE_NODES",
    "N5_CONE_SUPPRESSED_METHODS",
    "N5_EXACT_ZERO_METHODS",
    "build_n5_basis",
    "compute_fct_perm100_null_r11",
    "compute_fct_perm100_split_r11",
    "compute_joint_cross_mask_e1_diff_r11",
    "compute_joint_paired_r_orth_diff_r11",
    "compute_joint_sa8_split_minus_all_r12",
    "compute_n5_perm100_r12",
    "compute_n5_resp_top50_split_r12",
    "compute_n5_suite_r12",
    "compute_paired_fct_minus_all_mismatched_r12",
    "compute_residual_suite_r11",
    "compute_resp_top50_split_suite_r11",
    "compute_s7_cellcount_and_kd_r11",
    "g_const_cone",
    "paired_row_bootstrap_ci_r11",
    "prepare_n5_per_pert_data",
    "single_half_ceiling",
    "summarize_boot_r11",
]
