#!/usr/bin/env python3
"""
Independent Round 12 Confirmatory Endpoints (`H-E1a/b/c`, `H-E2a/b`), IUT + 2-Claim
Bonferroni Adjudication, and Full §5 Non-Confirmatory Suite (A-T3).

Strictly follows `PREREG_E.md` (`v2`, commit `4c606f4`):
  - Reuses frozen `round11/audit/` core primitives (`evaluate_joint_residuals`,
    `evaluate_paired_cal_diff`, `evaluate_fct_perm_null`, `evaluate_sa8_resp_top50_split`,
    `evaluate_sa8_perm_null`, `evaluate_hd_util_noninferiority`, `evaluate_r19_gain_decomposition`,
    `evaluate_s7_confounder_controls`) and frozen `exploratory_check/exploratory_core.py`
    (`evaluate_cross_mask_split_minus_all`).
  - Defaults: `seed = 20261201`, `n_boot = 4,000`, `n_perm = 100`, two-sided Bonferroni
    `97.5%` percentile CIs (`[1.25%, 98.75%]`, one-sided alpha = 0.0125 per tail).
  - Never reads, imports, or references `round12/src/` or `round12/results/`.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

AUDIT_R12_DIR = Path(__file__).resolve().parent
REPO_ROOT = AUDIT_R12_DIR.parents[1]
R11_AUDIT_DIR = REPO_ROOT / "round11" / "audit"
R11_EXP_DIR = R11_AUDIT_DIR / "exploratory_check"
for _p in (str(R11_AUDIT_DIR), str(R11_EXP_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from exploratory_core import evaluate_cross_mask_split_minus_all  # noqa: E402
from independent_r11_loader import AuditLineData  # noqa: E402
from independent_r11_metrics import (  # noqa: E402
    EPS,
    _extract_pert_gene_indices,
    _extract_sa8_masks,
    ceiling_from_halves,
    ceiling_single_half,
    disattenuate,
    evaluate_cross_mask_cal_diff,
    evaluate_fct_perm_null,
    evaluate_hd_util_noninferiority,
    evaluate_joint_residuals,
    evaluate_paired_cal_diff,
    evaluate_r19_gain_decomposition,
    evaluate_s7_confounder_controls,
    evaluate_sa8_perm_null,
    evaluate_sa8_resp_top50_split,
    summarize_boot,
)
from independent_r11_predictors import (  # noqa: E402
    AlignedCohort,
    _align_line_cols,
    build_predictions,
    fct_from_fc_bar,
)

DEFAULT_R12_SEED = 20261201

N5_NON_EVIDENTIARY_TAGS: Dict[str, str] = {
    "transferred_mean": "ZERO_BY_EXACT_CONSTRUCTION (residual identity in N_5)",
    "transfer_mean_global_scale": "ZERO_BY_EXACT_CONSTRUCTION (r_pred in span{t_mean} subset N_5; eta=0)",
    "source_grand_mean": "ZERO_BY_EXACT_CONSTRUCTION (r_pred = b_S - t_mean in N_5; eta=0)",
    "e1_const_med": "SUPPRESSED_BY_CONE_CONSTRUCTION (g_{c_med} in conv{g_0.05, g_1.0, g_20.0}; non-evidentiary)",
    "e2_const_p10": "SUPPRESSED_BY_CONE_CONSTRUCTION (g_{p10} in conv{g_0.05, g_1.0, g_20.0}; non-evidentiary)",
    "e2_const_p90": "SUPPRESSED_BY_CONE_CONSTRUCTION (g_{p90} in conv{g_0.05, g_1.0, g_20.0}; non-evidentiary)",
    "e6_log_fc_bar": "SUPPRESSED_BY_CONE_CONSTRUCTION (c -> inf boundary ray of N_5 cone; non-evidentiary)",
    "e7_fc_minus_1": "SUPPRESSED_BY_CONE_CONSTRUCTION (c -> 0+ boundary ray of N_5 cone; non-evidentiary)",
    "e3_src_mean_ctrl": "INFORMATIVE_SOURCE_BASELINE_REFERENCE",
    "e5_rank_matched_mu": "INFORMATIVE_TARGET_RANK_BASELINE",
    "fct": "PRIMARY_MODEL_DESCRIPTIVE_N5",
    "transfer_mean_var_scaled": "INFORMATIVE_NON_CONE_METHOD",
    "hybrid": "INFORMATIVE_NON_CONE_METHOD",
}


def _n5_g_const(fc_vec: np.ndarray, c: float) -> np.ndarray:
    """Centered constant-baseline transform g_c(fc_bar) = log1p(c*fc_bar) - log1p(c)."""
    fc_safe = np.maximum(np.asarray(fc_vec, dtype=np.float64), 0.0)
    return np.log1p(c * fc_safe) - math.log1p(c)


def _build_n5_orthonormal_basis(
    t_mean_row: np.ndarray,
    grand_mean_row: np.ndarray,
    fc_bar_row: np.ndarray,
    rcond: float = 1e-12,
) -> np.ndarray:
    """
    Construct canonical N_5(p) orthonormal basis spanning
      {1, t_mean[p], grand_mean, g_0.05(fc_bar[p]), g_1.0(fc_bar[p]), g_20.0(fc_bar[p])}
    by explicitly separating the unit intercept vector, centering/unit-normalizing the remaining
    columns, and applying SVD so projected residuals are exact zero-mean per perturbation.
    """
    m = t_mean_row.shape[0]
    ones = np.full(m, 1.0 / math.sqrt(float(m)), dtype=np.float64)
    cols = [
        np.asarray(t_mean_row, dtype=np.float64),
        np.asarray(grand_mean_row, dtype=np.float64),
        _n5_g_const(fc_bar_row, 0.05),
        _n5_g_const(fc_bar_row, 1.0),
        _n5_g_const(fc_bar_row, 20.0),
    ]
    orth_cols: List[np.ndarray] = []
    for v in cols:
        vp = v - float(np.dot(ones, v)) * ones
        nrm = float(np.linalg.norm(vp))
        if nrm > EPS:
            orth_cols.append(vp / nrm)
    if not orth_cols:
        return ones[:, None]
    mat = np.column_stack(orth_cols)
    u, s, _ = np.linalg.svd(mat, full_matrices=False)
    keep = s > max(EPS, rcond * float(s[0]))
    if not np.any(keep):
        return ones[:, None]
    return np.column_stack([ones, u[:, keep]])


def _project_perp(vec: np.ndarray, q: np.ndarray) -> np.ndarray:
    if q.shape[1] == 0:
        return vec.copy()
    return vec - q @ (q.T @ vec)


def _eval_n5_single_stream_moments(
    obs_list: Sequence[np.ndarray],
    pred_cone_list: Sequence[np.ndarray],
    pred_cent_list: Sequence[np.ndarray],
    ha_list: Sequence[np.ndarray],
    hb_list: Sequence[np.ndarray],
    draws: Optional[np.ndarray],
    n_boot: int,
    single_half_ceil: bool,
) -> Dict[str, Any]:
    """
    Compute canonical N_5 semi-partial residual correlation:
      r_cone_sp_raw = sum <r_obs_cone, r_pred_cone> / sqrt(sum ||r_obs_cone||^2 * sum ||r_pred_cent||^2),
      eta           = sqrt(sum ||r_pred_cone||^2 / sum ||r_pred_cent||^2),
      R2_in_N5      = 1 - eta^2,
    disattenuated by the identically projected split-half noise ceiling.
    """
    m = len(obs_list)
    if m == 0:
        empty = summarize_boot(float("nan"), np.array([float("nan")]), 0, n_boot).to_dict("r_cone_sp")
        empty["_b_arr"] = np.full(n_boot, float("nan"), dtype=np.float64)
        empty["eta_est"] = float("nan")
        empty["r2_in_n5"] = float("nan")
        empty["rho_ab_n5"] = float("nan")
        return empty

    xy = np.asarray([float(np.dot(x, y)) for x, y in zip(obs_list, pred_cone_list)], dtype=np.float64)
    xx = np.asarray([float(np.dot(x, x)) for x in obs_list], dtype=np.float64)
    yy_cone = np.asarray([float(np.dot(y, y)) for y in pred_cone_list], dtype=np.float64)
    yy_cent = np.asarray([float(np.dot(y, y)) for y in pred_cent_list], dtype=np.float64)
    haa = np.asarray([float(np.dot(a, a)) for a in ha_list], dtype=np.float64)
    hbb = np.asarray([float(np.dot(b, b)) for b in hb_list], dtype=np.float64)
    hab = np.asarray([float(np.dot(a, b)) for a, b in zip(ha_list, hb_list)], dtype=np.float64)
    has_halves = bool(np.any(haa > 0.0) and np.any(hbb > 0.0))

    def _calc(idx_mat: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        if idx_mat is None:
            s_xy = np.sum(xy)[None]
            s_xx = np.sum(xx)[None]
            s_yy_cone = np.sum(yy_cone)[None]
            s_yy_cent = np.sum(yy_cent)[None]
            s_haa = np.sum(haa)[None]
            s_hbb = np.sum(hbb)[None]
            s_hab = np.sum(hab)[None]
        else:
            s_xy = xy[idx_mat].sum(axis=1)
            s_xx = xx[idx_mat].sum(axis=1)
            s_yy_cone = yy_cone[idx_mat].sum(axis=1)
            s_yy_cent = yy_cent[idx_mat].sum(axis=1)
            s_haa = haa[idx_mat].sum(axis=1)
            s_hbb = hbb[idx_mat].sum(axis=1)
            s_hab = hab[idx_mat].sum(axis=1)

        denom_sp = np.sqrt(np.maximum(s_xx, 0.0) * np.maximum(s_yy_cent, 0.0))
        r_sp_raw = np.where((denom_sp > EPS) & (s_yy_cone > 1e-24), np.clip(s_xy / np.maximum(denom_sp, EPS), -1.0, 1.0), 0.0)
        eta = np.where(s_yy_cent > EPS, np.sqrt(np.clip(s_yy_cone / np.maximum(s_yy_cent, EPS), 0.0, 1.0)), 0.0)

        if has_halves:
            denom_h = np.sqrt(np.maximum(s_haa * s_hbb, 0.0))
            rho_h = np.where(denom_h > EPS, np.clip(s_hab / np.maximum(denom_h, EPS), -1.0, 1.0), np.nan)
            fn_ceil = ceiling_single_half if single_half_ceil else ceiling_from_halves
            ceil_arr = np.asarray([fn_ceil(float(v)) for v in rho_h], dtype=np.float64)
            r_sp_cal = np.where(
                s_yy_cone <= 1e-24,
                0.0,
                np.asarray([disattenuate(float(r), float(c)) for r, c in zip(r_sp_raw, ceil_arr)], dtype=np.float64),
            )
        else:
            rho_h = np.full_like(r_sp_raw, np.nan)
            r_sp_cal = np.full_like(r_sp_raw, np.nan)
        return r_sp_cal, r_sp_raw, eta, rho_h

    pt_cal, pt_raw, pt_eta, pt_rho = _calc(None)
    if draws is not None:
        b_cal, b_raw, _, _ = _calc(draws)
    else:
        b_cal = np.array([float("nan")])
        b_raw = np.array([float("nan")])

    out = summarize_boot(float(pt_cal[0]), b_cal, m, n_boot).to_dict("r_cone_sp")
    out.update(summarize_boot(float(pt_raw[0]), b_raw, m, n_boot).to_dict("r_cone_sp_raw"))
    out["eta_est"] = float(pt_eta[0])
    out["r2_in_n5"] = float(max(0.0, 1.0 - float(pt_eta[0]) ** 2))
    out["rho_ab_n5"] = float(pt_rho[0])
    out["n_eval_perts"] = m
    out["_b_arr"] = b_cal
    return out


def evaluate_n5_suite_r12(
    delta_obs: np.ndarray,
    preds: Dict[str, np.ndarray],
    t_mean: np.ndarray,
    grand_mean: np.ndarray,
    fc_bar: np.ndarray,
    tgt_ctrl_mean_libnorm: np.ndarray,
    target_cols: Sequence[int],
    source_names: Sequence[str],
    half_a_obs: Optional[np.ndarray],
    half_b_obs: Optional[np.ndarray],
    half_valid: Optional[np.ndarray],
    mode: str = "all_nontarget",
    n_boot: int = 4000,
    n_perm: int = 100,
    seed: int = DEFAULT_R12_SEED,
    top_k: int = 50,
) -> Dict[str, Any]:
    """
    Independent `PREREG_E §5` Item 3 `N_5` semi-partial analysis (`r_cone_sp`, `eta`, `R2_in_N5`)
    supporting BOTH `mode = "all_nontarget"` and symmetric `mode = "resp_top50_split"` (`(A->B + B->A)/2`):
      - Evaluates `FCT`, `6` mismatched `Mis_s` (unpaired & paired `FCT - Mis_s`), `Perm0..99` single-draw
        empirical distribution (never expected), `TMVS`, `Hybrid`, `E3`, `E5`, and construction-nulled
        controls (`Transfer_Mean`, `GlobalScale`, `SGM`, `E1`, `E2`, `E6`, `E7`) explicitly tagged
        as non-evidentiary baselines.
    """
    if mode not in ("all_nontarget", "resp_top50_split"):
        raise ValueError(f"Unsupported N_5 mode: {mode!r}")

    core_keys = [
        k for k in (
            "fct",
            "transferred_mean",
            "transfer_mean_global_scale",
            "source_grand_mean",
            "e1_const_med",
            "e2_const_p10",
            "e2_const_p90",
            "e6_log_fc_bar",
            "e7_fc_minus_1",
            "e3_src_mean_ctrl",
            "e5_rank_matched_mu",
            "transfer_mean_var_scaled",
            "hybrid",
        )
        if k in preds
    ]
    eval_keys = core_keys + [f"fct_mismatch::{s}" for s in source_names if f"fct_mismatch::{s}" in preds]
    mu_safe = np.maximum(np.asarray(tgt_ctrl_mean_libnorm, dtype=np.float64), 0.0)
    n_g = delta_obs.shape[1]

    if mode == "all_nontarget":
        kept = _extract_pert_gene_indices(delta_obs, t_mean, preds["fct"], target_cols, mode="all_nontarget")
        m = len(kept)
        if m == 0:
            return {"mode": mode, "n_eval_perts": 0, "methods": {}, "paired_fct_minus_mis": {}, "perm100": {}}
        rng = np.random.default_rng(seed)
        draws = rng.integers(0, m, size=(n_boot, m))

        bases, obs_c, ha_c, hb_c = [], [], [], []
        for p, idx in kept:
            ref = t_mean[p, idx].astype(np.float64)
            q_n5 = _build_n5_orthonormal_basis(ref, grand_mean[idx].astype(np.float64), fc_bar[p, idx].astype(np.float64))
            ro = delta_obs[p, idx].astype(np.float64) - ref
            bases.append(q_n5)
            obs_c.append(_project_perp(ro, q_n5))
            if (
                half_a_obs is not None
                and half_b_obs is not None
                and (half_valid is None or bool(half_valid[p]))
                and np.all(np.isfinite(half_a_obs[p, idx]) & np.isfinite(half_b_obs[p, idx]))
            ):
                ha_c.append(_project_perp(half_a_obs[p, idx].astype(np.float64) - ref, q_n5))
                hb_c.append(_project_perp(half_b_obs[p, idx].astype(np.float64) - ref, q_n5))
            else:
                ha_c.append(np.zeros(0, dtype=np.float64))
                hb_c.append(np.zeros(0, dtype=np.float64))

        methods_out: Dict[str, Any] = {}
        boot_map: Dict[str, np.ndarray] = {}
        for k in eval_keys:
            pmat = preds[k]
            yp_cone, yp_cent = [], []
            for j, (p, idx) in enumerate(kept):
                rp = pmat[p, idx].astype(np.float64) - t_mean[p, idx].astype(np.float64)
                yp_cone.append(_project_perp(rp, bases[j]))
                yp_cent.append(rp - float(np.mean(rp)))
            res = _eval_n5_single_stream_moments(obs_c, yp_cone, yp_cent, ha_c, hb_c, draws, n_boot, single_half_ceil=False)
            boot_map[k] = res.pop("_b_arr")
            res["audit_classification"] = N5_NON_EVIDENTIARY_TAGS.get(k, "INFORMATIVE_MISMATCH_BASELINE")
            methods_out[k] = res

        # Paired N_5 FCT - Mis_s
        paired_mis: Dict[str, Any] = {}
        for s in source_names:
            mk = f"fct_mismatch::{s}"
            if mk not in methods_out:
                continue
            pt_d = float(methods_out["fct"]["r_cone_sp_est"]) - float(methods_out[mk]["r_cone_sp_est"])
            b_d = boot_map["fct"] - boot_map[mk]
            ds = summarize_boot(pt_d, b_d, m, n_boot).to_dict("diff")
            ds["paired_ci975_gt_0"] = bool(np.isfinite(ds["diff_ci975_lo"]) and ds["diff_ci975_lo"] > 0.0)
            paired_mis[s] = ds

        # Single-draw Perm0..99 (non-expected)
        p_cal = np.empty(n_perm, dtype=np.float64)
        p_raw = np.empty(n_perm, dtype=np.float64)
        p_eta = np.empty(n_perm, dtype=np.float64)
        for pk in range(n_perm):
            perm = np.random.default_rng(pk).permutation(n_g)
            pred_pk = fct_from_fc_bar(fc_bar, mu_safe[perm])
            yp_cone, yp_cent = [], []
            for j, (p, idx) in enumerate(kept):
                rp = pred_pk[p, idx] - t_mean[p, idx].astype(np.float64)
                yp_cone.append(_project_perp(rp, bases[j]))
                yp_cent.append(rp - float(np.mean(rp)))
            rk = _eval_n5_single_stream_moments(obs_c, yp_cone, yp_cent, ha_c, hb_c, None, 1, single_half_ceil=False)
            p_cal[pk] = rk["r_cone_sp_est"]
            p_raw[pk] = rk["r_cone_sp_raw_est"]
            p_eta[pk] = rk["eta_est"]

        return {
            "mode": mode,
            "n_eval_perts": m,
            "methods": methods_out,
            "paired_fct_minus_mis": paired_mis,
            "perm100": _summarize_n5_perm_draws(p_cal, p_raw, p_eta),
        }

    # Mode == "resp_top50_split": symmetric A->B and B->A single-half ceiling disattenuation
    masks = _extract_sa8_masks(
        delta_obs=delta_obs,
        t_mean=t_mean,
        pred=preds["fct"],
        target_cols=target_cols,
        half_a_obs=half_a_obs,
        half_b_obs=half_b_obs,
        half_valid=half_valid,
        top_k=top_k,
    )
    m = len(masks)
    if m == 0 or half_a_obs is None or half_b_obs is None:
        return {"mode": mode, "n_eval_perts": 0, "methods": {}, "paired_fct_minus_mis": {}, "perm100": {}}
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, m, size=(n_boot, m))

    q_n5_a, obs_b_on_a, ha_a, hb_a = [], [], [], []
    q_n5_b, obs_a_on_b, ha_b, hb_b = [], [], [], []
    for p, ia, ib in masks:
        ref_a = t_mean[p, ia].astype(np.float64)
        qa = _build_n5_orthonormal_basis(ref_a, grand_mean[ia].astype(np.float64), fc_bar[p, ia].astype(np.float64))
        q_n5_a.append(qa)
        ha_a.append(_project_perp(half_a_obs[p, ia].astype(np.float64) - ref_a, qa))
        hb_va = _project_perp(half_b_obs[p, ia].astype(np.float64) - ref_a, qa)
        hb_a.append(hb_va)
        obs_b_on_a.append(hb_va)

        ref_b = t_mean[p, ib].astype(np.float64)
        qb = _build_n5_orthonormal_basis(ref_b, grand_mean[ib].astype(np.float64), fc_bar[p, ib].astype(np.float64))
        q_n5_b.append(qb)
        ha_vb = _project_perp(half_a_obs[p, ib].astype(np.float64) - ref_b, qb)
        hb_b.append(_project_perp(half_b_obs[p, ib].astype(np.float64) - ref_b, qb))
        ha_b.append(ha_vb)
        obs_a_on_b.append(ha_vb)

    methods_out = {}
    boot_map = {}
    for k in eval_keys:
        pmat = preds[k]
        yp_ca, yp_cent_a, yp_cb, yp_cent_b = [], [], [], []
        for j, (p, ia, ib) in enumerate(masks):
            rpa = pmat[p, ia].astype(np.float64) - t_mean[p, ia].astype(np.float64)
            yp_ca.append(_project_perp(rpa, q_n5_a[j]))
            yp_cent_a.append(rpa - float(np.mean(rpa)))
            rpb = pmat[p, ib].astype(np.float64) - t_mean[p, ib].astype(np.float64)
            yp_cb.append(_project_perp(rpb, q_n5_b[j]))
            yp_cent_b.append(rpb - float(np.mean(rpb)))

        res_ab = _eval_n5_single_stream_moments(obs_b_on_a, yp_ca, yp_cent_a, ha_a, hb_a, draws, n_boot, single_half_ceil=True)
        res_ba = _eval_n5_single_stream_moments(obs_a_on_b, yp_cb, yp_cent_b, ha_b, hb_b, draws, n_boot, single_half_ceil=True)
        b_ab = res_ab.pop("_b_arr")
        b_ba = res_ba.pop("_b_arr")
        pt_avg = 0.5 * (float(res_ab["r_cone_sp_est"]) + float(res_ba["r_cone_sp_est"]))
        b_avg = 0.5 * (b_ab + b_ba)
        out_k = summarize_boot(pt_avg, b_avg, m, n_boot).to_dict("r_cone_sp")
        out_k["r_cone_sp_raw_est"] = 0.5 * (float(res_ab["r_cone_sp_raw_est"]) + float(res_ba["r_cone_sp_raw_est"]))
        out_k["eta_est"] = 0.5 * (float(res_ab["eta_est"]) + float(res_ba["eta_est"]))
        out_k["r2_in_n5"] = 0.5 * (float(res_ab["r2_in_n5"]) + float(res_ba["r2_in_n5"]))
        out_k["rho_ab_n5"] = 0.5 * (float(res_ab["rho_ab_n5"]) + float(res_ba["rho_ab_n5"]))
        out_k["n_eval_perts"] = m
        out_k["audit_classification"] = N5_NON_EVIDENTIARY_TAGS.get(k, "INFORMATIVE_MISMATCH_BASELINE")
        boot_map[k] = b_avg
        methods_out[k] = out_k

    paired_mis = {}
    for s in source_names:
        mk = f"fct_mismatch::{s}"
        if mk not in methods_out:
            continue
        pt_d = float(methods_out["fct"]["r_cone_sp_est"]) - float(methods_out[mk]["r_cone_sp_est"])
        b_d = boot_map["fct"] - boot_map[mk]
        ds = summarize_boot(pt_d, b_d, m, n_boot).to_dict("diff")
        ds["paired_ci975_gt_0"] = bool(np.isfinite(ds["diff_ci975_lo"]) and ds["diff_ci975_lo"] > 0.0)
        paired_mis[s] = ds

    p_cal = np.empty(n_perm, dtype=np.float64)
    p_raw = np.empty(n_perm, dtype=np.float64)
    p_eta = np.empty(n_perm, dtype=np.float64)
    for pk in range(n_perm):
        perm = np.random.default_rng(pk).permutation(n_g)
        pred_pk = fct_from_fc_bar(fc_bar, mu_safe[perm])
        yp_ca, yp_cent_a, yp_cb, yp_cent_b = [], [], [], []
        for j, (p, ia, ib) in enumerate(masks):
            rpa = pred_pk[p, ia] - t_mean[p, ia].astype(np.float64)
            yp_ca.append(_project_perp(rpa, q_n5_a[j]))
            yp_cent_a.append(rpa - float(np.mean(rpa)))
            rpb = pred_pk[p, ib] - t_mean[p, ib].astype(np.float64)
            yp_cb.append(_project_perp(rpb, q_n5_b[j]))
            yp_cent_b.append(rpb - float(np.mean(rpb)))
        ra = _eval_n5_single_stream_moments(obs_b_on_a, yp_ca, yp_cent_a, ha_a, hb_a, None, 1, single_half_ceil=True)
        rb = _eval_n5_single_stream_moments(obs_a_on_b, yp_cb, yp_cent_b, ha_b, hb_b, None, 1, single_half_ceil=True)
        p_cal[pk] = 0.5 * (float(ra["r_cone_sp_est"]) + float(rb["r_cone_sp_est"]))
        p_raw[pk] = 0.5 * (float(ra["r_cone_sp_raw_est"]) + float(rb["r_cone_sp_raw_est"]))
        p_eta[pk] = 0.5 * (float(ra["eta_est"]) + float(rb["eta_est"]))

    return {
        "mode": mode,
        "n_eval_perts": m,
        "methods": methods_out,
        "paired_fct_minus_mis": paired_mis,
        "perm100": _summarize_n5_perm_draws(p_cal, p_raw, p_eta),
    }


def _summarize_n5_perm_draws(p_cal: np.ndarray, p_raw: np.ndarray, p_eta: np.ndarray) -> Dict[str, Any]:
    fin = p_cal[np.isfinite(p_cal)]
    return {
        "n_perm": int(p_cal.size),
        "min": float(np.min(fin)) if fin.size else float("nan"),
        "p025": float(np.percentile(fin, 2.5)) if fin.size else float("nan"),
        "median": float(np.median(fin)) if fin.size else float("nan"),
        "mean": float(np.mean(fin)) if fin.size else float("nan"),
        "p975": float(np.percentile(fin, 97.5)) if fin.size else float("nan"),
        "max": float(np.max(fin)) if fin.size else float("nan"),
        "median_abs": float(np.median(np.abs(fin))) if fin.size else float("nan"),
        "median_raw": float(np.nanmedian(p_raw)),
        "median_eta": float(np.nanmedian(p_eta)),
    }


def align_source_and_target_r12(
    source_lines: Sequence[AuditLineData],
    target_data: AuditLineData,
    only_reproducible: bool = True,
) -> AlignedCohort:
    """
    Construct joint `G_eval(t)` intersecting target genes with ALL source lines in the pool
    (necessary when `r9_r10_8` extension sources miss a few `R9_COMMON` genes) and apply `P5`
    target reproducibility + self-gene presence inside `G_eval(t)`.
    """
    src_gene_sets = [set(s.genes) for s in source_lines]
    tgt_keep_cols = [i for i, g in enumerate(target_data.genes) if all(g in ss for ss in src_gene_sets)]
    eval_genes = [target_data.genes[i] for i in tgt_keep_cols]
    tgt_col_slice = np.asarray(tgt_keep_cols, dtype=np.int64)
    eval_col_map = {g: j for j, g in enumerate(eval_genes)}

    src_deltas: Dict[str, Dict[str, np.ndarray]] = {}
    src_ctrl: Dict[str, np.ndarray] = {}
    src_mu_lib: Dict[str, np.ndarray] = {}
    src_avail_targets = set()

    for s in source_lines:
        cols = _align_line_cols(s.genes, eval_genes)
        src_ctrl[s.name] = np.asarray(s.ctrl_cells[:, cols], dtype=np.float64)
        src_mu_lib[s.name] = np.asarray(s.ctrl_mean_libnorm[cols], dtype=np.float64)
        s_map: Dict[str, np.ndarray] = {}
        for i, tgt in enumerate(s.target_names):
            if only_reproducible and s.repr_mask is not None and not bool(s.repr_mask[i]):
                continue
            s_map[tgt] = np.asarray(s.deltas[i, cols], dtype=np.float64)
            src_avail_targets.add(tgt)
        src_deltas[s.name] = s_map

    kept_rows: List[int] = []
    kept_target_cols: List[int] = []
    for i, tgt in enumerate(target_data.target_names):
        tc = eval_col_map.get(tgt, -1)
        if tc < 0:
            continue
        if only_reproducible and target_data.repr_mask is not None and not bool(target_data.repr_mask[i]):
            continue
        if tgt not in src_avail_targets:
            continue
        kept_rows.append(i)
        kept_target_cols.append(tc)

    if not kept_rows:
        raise RuntimeError(f"Zero shared eligible targets for {target_data.name} (only_reproducible={only_reproducible}).")

    idx_arr = np.asarray(kept_rows, dtype=np.int64)
    eval_targets = [target_data.target_names[i] for i in kept_rows]
    grid = np.ix_(idx_arr, tgt_col_slice)

    return AlignedCohort(
        target_line=target_data.name,
        eval_genes=eval_genes,
        eval_targets=eval_targets,
        target_cols=np.asarray(kept_target_cols, dtype=np.int64),
        delta_obs=np.asarray(target_data.deltas[grid], dtype=np.float64),
        n_cells=np.asarray(target_data.n_cells[idx_arr], dtype=np.float64) if target_data.n_cells is not None else None,
        kd_depth_log1p=np.asarray(target_data.kd_depth_log1p[idx_arr], dtype=np.float64),
        half_a_obs=np.asarray(target_data.half_a_deltas[grid], dtype=np.float64) if target_data.half_a_deltas is not None else None,
        half_b_obs=np.asarray(target_data.half_b_deltas[grid], dtype=np.float64) if target_data.half_b_deltas is not None else None,
        half_valid=np.asarray(target_data.half_valid[idx_arr], dtype=bool) if target_data.half_valid is not None else None,
        tgt_ctrl_cells=np.asarray(target_data.ctrl_cells[:, tgt_col_slice], dtype=np.float64),
        tgt_ctrl_mean_libnorm=np.asarray(target_data.ctrl_mean_libnorm[tgt_col_slice], dtype=np.float64),
        tgt_ctrl_mean_libnorm_full=(
            np.asarray(target_data.ctrl_mean_libnorm_full[tgt_col_slice], dtype=np.float64)
            if target_data.ctrl_mean_libnorm_full is not None
            else None
        ),
        source_names=[s.name for s in source_lines],
        source_deltas_by_line=src_deltas,
        source_ctrl_cells=src_ctrl,
        source_ctrl_mean_libnorm=src_mu_lib,
    )


def evaluate_target_line_r12(
    target_data: AuditLineData,
    source_list: Sequence[AuditLineData],
    only_reproducible: bool = True,
    shrink: float = 0.0,
    cov_max_cells: Optional[int] = None,
    n_boot: int = 4000,
    n_perm: int = 100,
    seed: int = DEFAULT_R12_SEED,
    include_sa3_anchors: bool = True,
    include_exploratory_e: bool = True,
    include_n5: bool = True,
    include_n5_split: bool = True,
    include_r19: bool = True,
    include_s7: bool = True,
    precomputed_covs: Optional[Dict[str, np.ndarray]] = None,
    precomputed_gene_orders: Optional[Dict[str, Sequence[str]]] = None,
) -> Dict[str, Any]:
    """
    Evaluate one target cell line cleanly on BOTH `PREREG_E.md` §4 confirmatory endpoints
    (`H-E1a`, all individual `H-E1b` paired `FCT - Mis_s`, `H-E1c`, `H-E2a`, `H-E2b`) AND
    all §5 pre-specified non-confirmatory descriptive analyses (`只报告、不判定`).
    """
    cohort = align_source_and_target_r12(source_list, target_data, only_reproducible=only_reproducible)
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

    # Evaluate primary `all_nontarget` and secondary in-sample `resp_top50`
    dual_methods = {"fct", "e1_const_med", "fct_full_nt", "transfer_mean_var_scaled", "hybrid"} | {
        f"fct_mismatch::{s}" for s in cohort.source_names
    }
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
            if m_name in dual_methods
            else {}
        )
        method_eval[m_name] = {"all_nontarget": m_all, "resp_top50": m_top50}

    fct_all = method_eval["fct"]["all_nontarget"]
    fct_top50 = method_eval["fct"]["resp_top50"]
    e1_all = method_eval["e1_const_med"]["all_nontarget"]
    e1_top50 = method_eval["e1_const_med"]["resp_top50"]

    # Perm100 on `all_nontarget` (confirmatory H-E1c) and `resp_top50` (descriptive §5.6)
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

    # ---- Confirmatory Claim 1: H-E1a, H-E1b (every source paired diff), H-E1c ----
    fct_all_lo = float(fct_all.get("r_orth_cal_ci975_lo", float("nan")))
    fct_all_vf = float(fct_all.get("r_orth_cal_valid_frac", 0.0))
    fct_all_rho = float(fct_all.get("rho_ab_orth", float("nan")))
    he1a_ci_ok = bool(np.isfinite(fct_all_lo) and fct_all_lo > 0.0)
    he1a_vf_ok = bool(np.isfinite(fct_all_vf) and fct_all_vf >= 0.95)
    he1a_rho_ok = bool(np.isfinite(fct_all_rho) and fct_all_rho > 0.0)
    he1a_pass = bool(he1a_ci_ok and he1a_vf_ok and he1a_rho_ok)

    mis_all: Dict[str, Dict[str, Any]] = {}
    paired_he1b: Dict[str, Dict[str, Any]] = {}
    he1b_flags: List[bool] = []
    mis_all_his: List[float] = []
    for s in cohort.source_names:
        mk = f"fct_mismatch::{s}"
        if mk not in preds:
            continue
        mis_res = method_eval[mk]["all_nontarget"]
        mis_all[s] = mis_res
        mis_all_his.append(float(mis_res.get("r_orth_cal_ci975_hi", float("nan"))))

        p_res = evaluate_paired_cal_diff(
            delta_obs=cohort.delta_obs,
            pred_a=preds["fct"],
            pred_b=preds[mk],
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
        d_lo = float(p_res.get("diff_ci975_lo", float("nan")))
        d_vf = float(p_res.get("diff_valid_frac", 0.0))
        d_pass = bool(np.isfinite(d_lo) and d_lo > 0.0 and d_vf > 0.0)
        p_res["passes_he1b_source"] = d_pass
        paired_he1b[s] = p_res
        he1b_flags.append(d_pass)

    he1b_pass = bool(len(he1b_flags) > 0 and all(he1b_flags))
    perm_all_p975 = float(perm100_all["r_orth_cal"]["p975"])
    he1c_pass = bool(np.isfinite(fct_all_lo) and np.isfinite(perm_all_p975) and fct_all_lo > perm_all_p975)
    he1_pass = bool(he1a_pass and he1b_pass and he1c_pass)

    # ---- Confirmatory Claim 2: H-E2a & H-E2b (`SA-8` resp_top50_split) ----
    sa8_fct = evaluate_sa8_resp_top50_split(
        delta_obs=cohort.delta_obs,
        pred=preds["fct"],
        t_mean=t_mean,
        target_cols=cohort.target_cols,
        half_a_obs=cohort.half_a_obs,
        half_b_obs=cohort.half_b_obs,
        half_valid=cohort.half_valid,
        top_k=50,
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
        top_k=50,
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
    sa8_e1_split_minus_all = evaluate_cross_mask_split_minus_all(
        delta_obs=cohort.delta_obs,
        pred=preds["e1_const_med"],
        t_mean=t_mean,
        target_cols=cohort.target_cols,
        half_a_obs=cohort.half_a_obs,
        half_b_obs=cohort.half_b_obs,
        half_valid=cohort.half_valid,
        top_k=50,
        n_boot=n_boot,
        seed=seed,
    )
    he2a_lo = float(sa8_e1.get("r_orth_cal_split_ci975_lo", float("nan")))
    he2a_vf = float(sa8_e1.get("r_orth_cal_split_valid_frac", 0.0))
    he2a_pass = bool(np.isfinite(he2a_lo) and he2a_lo > 0.0 and he2a_vf > 0.0)
    he2b_lo = float(sa8_e1_split_minus_all.get("diff_ci975_lo", float("nan")))
    he2b_vf = float(sa8_e1_split_minus_all.get("diff_valid_frac", 0.0))
    he2b_pass = bool(np.isfinite(he2b_lo) and he2b_lo > 0.0 and he2b_vf > 0.0)
    he2_pass = bool(he2a_pass and he2b_pass)

    # ---- §5 Non-Confirmatory Descriptive Analyses (`只报告、不判定`) ----
    # 5.1 Legacy unpaired H-D2
    max_mis_hi = float(np.nanmax(mis_all_his)) if mis_all_his else float("nan")
    legacy_hd2_pass = bool(
        np.isfinite(fct_all_lo)
        and np.isfinite(max_mis_hi)
        and np.isfinite(perm_all_p975)
        and (fct_all_lo > max(max_mis_hi, perm_all_p975))
    )

    # 5.2 Paired FCT - E1_ConstMed (`all_nontarget` + `resp_top50`)
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

    # 5.3 N_5 Semi-partial analysis on all_nontarget & resp_top50_split
    n5_all = (
        evaluate_n5_suite_r12(
            delta_obs=cohort.delta_obs,
            preds=preds,
            t_mean=t_mean,
            grand_mean=grand_mean,
            fc_bar=fc_bar,
            tgt_ctrl_mean_libnorm=cohort.tgt_ctrl_mean_libnorm,
            target_cols=cohort.target_cols,
            source_names=cohort.source_names,
            half_a_obs=cohort.half_a_obs,
            half_b_obs=cohort.half_b_obs,
            half_valid=cohort.half_valid,
            mode="all_nontarget",
            n_boot=n_boot,
            n_perm=n_perm,
            seed=seed,
        )
        if include_n5
        else {}
    )
    n5_split = (
        evaluate_n5_suite_r12(
            delta_obs=cohort.delta_obs,
            preds=preds,
            t_mean=t_mean,
            grand_mean=grand_mean,
            fc_bar=fc_bar,
            tgt_ctrl_mean_libnorm=cohort.tgt_ctrl_mean_libnorm,
            target_cols=cohort.target_cols,
            source_names=cohort.source_names,
            half_a_obs=cohort.half_a_obs,
            half_b_obs=cohort.half_b_obs,
            half_valid=cohort.half_valid,
            mode="resp_top50_split",
            n_boot=n_boot,
            n_perm=n_perm,
            seed=seed,
        )
        if (include_n5 and include_n5_split)
        else {}
    )

    # 5.6 In-sample resp_top50 cross-mask E1(resp_top50) - E1(all_nontarget)
    cross_mask_e1_insample = evaluate_cross_mask_cal_diff(
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

    # 5.7 H-D_util (6+1 conventions), SA-5 top-75%, R19, S7
    hd_util = evaluate_hd_util_noninferiority(
        delta_obs=cohort.delta_obs,
        pred_fct=preds["fct"],
        pred_tmean=t_mean,
        grand_mean=grand_mean,
        target_cols=cohort.target_cols,
        n_boot=n_boot,
        seed=seed,
    )
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

    return {
        "target_line": cohort.target_line,
        "n_genes_eval": len(cohort.eval_genes),
        "n_eval_targets": len(cohort.eval_targets),
        "n_sources": len(cohort.source_names),
        "source_names": list(cohort.source_names),
        "const_stats": const_stats,
        "method_eval": method_eval,
        "mismatched_all": mis_all,
        "paired_he1b": paired_he1b,
        "perm100_all_nontarget": perm100_all,
        "perm100_resp_top50": perm100_top50,
        "sa8_fct": sa8_fct,
        "sa8_e1": sa8_e1,
        "sa8_perm100": sa8_perm100,
        "sa8_e1_split_minus_all": sa8_e1_split_minus_all,
        "he1a_pass": he1a_pass,
        "he1a_details": {"ci975_gt_0": he1a_ci_ok, "valid_frac_ge_95": he1a_vf_ok, "rho_ab_orth_gt_0": he1a_rho_ok},
        "he1b_pass": he1b_pass,
        "he1b_n_pass": int(sum(he1b_flags)),
        "he1b_n_total": len(he1b_flags),
        "he1c_pass": he1c_pass,
        "he1_pass": he1_pass,
        "he2a_pass": he2a_pass,
        "he2b_pass": he2b_pass,
        "he2_pass": he2_pass,
        "legacy_hd2_pass": legacy_hd2_pass,
        "max_mismatch_all_ci975_hi": max_mis_hi,
        "pair_fct_minus_e1_all": pair_fct_minus_e1_all,
        "pair_fct_minus_e1_top50": pair_fct_minus_e1_top50,
        "cross_mask_e1_insample": cross_mask_e1_insample,
        "n5_all": n5_all,
        "n5_split": n5_split,
        "hd_util": hd_util,
        "sa5_top75_all_fct": sa5_top75_all,
        "r19_all": r19_all,
        "s7_all": s7_all,
        "e1_all_summary": e1_all,
        "e1_top50_summary": e1_top50,
    }


def adjudicate_round12(
    per_line_summary: Dict[str, Dict[str, Any]],
    gate_info: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Apply PREREG_E.md (§0.2, §2.2, §4.2, §4.3) IUT + 2-Claim Bonferroni adjudication:
      - If `a549_ko` fails mandatory count/ingestion gate (`k_conf == 0`): abort confirmatory
        adjudication and emit descriptive results only.
      - If `k_conf == 1` (`a549_ko` only): execute confirmatory adjudication with explicit
        mandatory limited-generalizability (`推广性有限`) scope disclosure.
      - If `k_conf == 2` (`a549_ko` + `hek293_xucao`): require `2/2` (`100%`) IUT passage per claim.
      - Asymmetric modality interpretation rule on `a549_ko` (`CRISPRi` sources -> `CRISPRko` A549).
    """
    conf_lines = list(gate_info.get("confirmatory_lines", []))
    k_conf = len(conf_lines)
    hek_role = str(gate_info.get("hek293_role", "SECONDARY_ONLY"))

    if k_conf == 0:
        return {
            "status": "CONFIRMATORY_ABORTED_K0",
            "k_conf": 0,
            "confirmatory_lines": [],
            "a549_role": gate_info.get("a549_role"),
            "hek293_role": hek_role,
            "claim1_he1_verdict": "ABORTED_MANDATORY_A549_FAILED",
            "claim2_he2_verdict": "ABORTED_MANDATORY_A549_FAILED",
            "claim1_confirmed": False,
            "claim2_confirmed": False,
            "limited_generalizability_notice": True,
            "modality_interpretation": "CONFIRMATORY_ABORTED_NO_VALID_A549",
            "mandatory_disclosure": (
                "Round 12 confirmatory adjudication was aborted (K_conf = 0) because mandatory line "
                "a549_ko (GSE345058) did not satisfy ingestion or >=50 targets / >=3,000 genes."
            ),
        }

    missing = [ln for ln in conf_lines if ln not in per_line_summary]
    if missing:
        raise RuntimeError(f"Missing confirmatory line summaries in adjudicate_round12: {missing}")

    c1_flags = {ln: bool(per_line_summary[ln]["he1_pass"]) for ln in conf_lines}
    c2_flags = {ln: bool(per_line_summary[ln]["he2_pass"]) for ln in conf_lines}
    c1_ok = all(c1_flags.values())
    c2_ok = all(c2_flags.values())

    a549_res = per_line_summary.get("a549_ko")
    if a549_res is None:
        mod_interp = "A549_NOT_EVALUATED"
    elif bool(a549_res["he1_pass"]) and bool(a549_res["he2_pass"]):
        mod_interp = "STRENGTHENS_CROSS_LINE_AND_CROSS_MODALITY_GENERALIZATION (CRISPRi -> CRISPRko A549)"
    else:
        mod_interp = (
            "STRICTLY_CONFOUNDED_CELL_LINE_VS_MODALITY: A549 (CRISPRko) failed Claim 1 and/or Claim 2; "
            "cannot distinguish lineage biology from CRISPRko vs CRISPRi perturbation modality."
        )

    if k_conf == 1:
        status = "CONFIRMATORY_SINGLE_LINE_K1"
        gen_notice = True
        scope_msg = (
            f"MANDATORY SCOPE NOTICE (K_conf = 1, 推广性有限): Confirmatory adjudication relies on a single line "
            f"(a549_ko, CRISPRko); hek293_xucao status = {hek_role}. Cross-line generalizability is limited."
        )
    else:
        status = "CONFIRMATORY_TWO_LINE_K2"
        gen_notice = False
        scope_msg = "Confirmatory adjudication evaluated unanimously across K_conf = 2 lines (a549_ko & hek293_xucao)."

    return {
        "status": status,
        "k_conf": k_conf,
        "confirmatory_lines": conf_lines,
        "a549_role": gate_info.get("a549_role", "CONFIRMATORY"),
        "hek293_role": hek_role,
        "claim1_he1_per_line": c1_flags,
        "claim2_he2_per_line": c2_flags,
        "claim1_confirmed": c1_ok,
        "claim2_confirmed": c2_ok,
        "claim1_he1_verdict": "CONFIRMED" if c1_ok else "NOT_CONFIRMED",
        "claim2_he2_verdict": "CONFIRMED" if c2_ok else "NOT_CONFIRMED",
        "limited_generalizability_notice": gen_notice,
        "modality_interpretation": mod_interp,
        "mandatory_disclosure": scope_msg,
    }
