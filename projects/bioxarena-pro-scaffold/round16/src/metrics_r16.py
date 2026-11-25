#!/usr/bin/env python3
"""
Round 16 Metrics Engine (`PREREG_I v5 + Addendum A`).

Strictly imports frozen mathematical primitives from `round10/src/` and `round11/src/`
without modifying any file in Rounds 9–15:
  1. Primary endpoint: `r_orth_cal` on `all_nontarget` (bit-compatible with frozen
     `round10/src/metrics_r10.py::_orth_pair` and `round12/src/metrics_r12.py`), with
     2-claim Bonferroni (`alpha_claim = 0.025`, two-sided `97.5%` CI `[1.25%, 98.75%]`),
     TOST equivalence (`delta = 0.04`, two-sided `95.0%` CI `[2.5%, 97.5%]`), and
     Addendum A Item 10 exact-invariance flag (`max_abs_diff <= 1e-12`).
  2. Secondary endpoint (`§4.1`, Systema-inspired shared-axis removal): `r_dorth_cal`,
     projecting residuals orthogonal to `span(1, delta_transfer_p, source_grand_mean)`
     via deterministic float64 two-pass modified Gram-Schmidt with rank guard (`1e-10`)
     and zero-energy guard (`1e-12`). By construction, all target-agnostic affine
     predictors in `span(1, delta_transfer_p, source_grand_mean)` evaluate to 0.0.
  3. Secondary diagnostics:
       - `§4.2` Mechanism Hypothesis M2: perturbation-level saturation-weight excess shift
         `compute_winner_curse_weight_shift_r16` and joint cross-mask shift
         `compute_joint_cross_mask_r16`.
       - Addendum A Item 4 (`ST-HVG` `Perm100`): `compute_gene_label_perm_null_r16`
         (`k = 0..99`, `np.random.default_rng(k).permutation(G)` on `G_eval`).
       - Addendum A Item 8 (within-line disjoint resample reference):
         `compute_paired_with_within_line_resample_r16` computing `delta_pair(resample)`
         and paired `delta_pair(s') - delta_pair(resample)` on identical bootstrap draws.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

ROUND16_SRC_DIR = Path(__file__).resolve().parent
REPO_ROOT = ROUND16_SRC_DIR.parents[1]
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
R11_SRC_DIR = REPO_ROOT / "round11" / "src"
for _p in (str(R10_SRC_DIR), str(R11_SRC_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from metrics_r10 import (  # noqa: E402
    _extract_masked_idx_per_pert,
    _moments,
    _moments_to_corr_and_slope,
    _orth_pair,
    disattenuated_r,
    spearman_brown_ceiling,
)
from metrics_r11 import (  # noqa: E402
    _prepare_split_direction_indices,
    single_half_ceiling,
)

DEFAULT_R16_SEED = 20261601
DEFAULT_RESAMPLE_SEED = 20261602
DEFAULT_TOST_DELTA = 0.04
EXACT_INVARIANCE_TOL = 1e-12
RANK_TOL = 1e-10
ZERO_ENERGY_TOL = 1e-12


def _array_max_abs_diff_and_exact_inv(
    a: np.ndarray,
    b: np.ndarray,
    tol: float = EXACT_INVARIANCE_TOL,
) -> Tuple[float, bool]:
    """Compute max absolute difference between two predicted delta matrices and Addendum A Item 10 exact-invariance flag."""
    arr_a = np.asarray(a, dtype=np.float64)
    arr_b = np.asarray(b, dtype=np.float64)
    if arr_a.shape != arr_b.shape or arr_a.size == 0:
        return float("nan"), False
    if not (np.all(np.isfinite(arr_a)) and np.all(np.isfinite(arr_b))):
        return float("nan"), False
    mad = float(np.max(np.abs(arr_a - arr_b)))
    return mad, bool(np.isfinite(mad) and mad <= float(tol))


def summarize_boot_r16(
    point: float,
    boots: np.ndarray,
    n_items: int,
    n_boot: int,
) -> Dict[str, Any]:
    """
    Summarize bootstrap replicates with:
      - `97.5%` CI (`[1.25%, 98.75%]`, 2-claim Bonferroni alpha_claim = 0.025 two-sided)
      - `95.0%` CI (`[2.5%, 97.5%]`, exact TOST interval at one-sided alpha_claim = 0.025)
    """
    fin = np.asarray(boots, dtype=np.float64)
    valid_mask = np.isfinite(fin)
    valid_frac = float(np.mean(valid_mask)) if fin.size > 0 else 0.0
    vals = fin[valid_mask]
    if vals.size == 0:
        return {
            "estimate": float(point),
            "ci95_lo": float("nan"),
            "ci95_hi": float("nan"),
            "ci975_lo": float("nan"),
            "ci975_hi": float("nan"),
            "boot_se": float("nan"),
            "valid_frac": valid_frac,
            "n_items": int(n_items),
            "n_boot": int(n_boot),
        }
    q0125, q025, q975, q9875 = np.percentile(vals, [1.25, 2.5, 97.5, 98.75])
    return {
        "estimate": float(point),
        "ci95_lo": float(q025),
        "ci95_hi": float(q975),
        "ci975_lo": float(q0125),
        "ci975_hi": float(q9875),
        "boot_se": float(np.std(vals, ddof=1)) if vals.size > 1 else 0.0,
        "valid_frac": valid_frac,
        "n_items": int(n_items),
        "n_boot": int(n_boot),
    }


def _build_centered_basis(
    ref_vec: np.ndarray,
    extra_basis_vec: np.ndarray,
    rank_tol: float = RANK_TOL,
) -> np.ndarray:
    """
    Build deterministic float64 orthonormal columns spanning centered `{ref_vec, extra_basis_vec}`
    via two-pass modified Gram-Schmidt. Returns shape `(m, r)` with `r in {0, 1, 2}`.
    """
    cols: List[np.ndarray] = []
    u1 = np.asarray(ref_vec, dtype=np.float64)
    u1 = u1 - float(np.mean(u1))
    n1 = float(np.linalg.norm(u1))
    if n1 > rank_tol:
        cols.append(u1 / n1)

    u2 = np.asarray(extra_basis_vec, dtype=np.float64)
    u2 = u2 - float(np.mean(u2))
    if cols:
        u2 = u2 - float(np.dot(cols[0], u2)) * cols[0]
        u2 = u2 - float(np.dot(cols[0], u2)) * cols[0]
    n2 = float(np.linalg.norm(u2))
    if n2 > rank_tol:
        cols.append(u2 / n2)

    if not cols:
        return np.zeros((u1.size, 0), dtype=np.float64)
    return np.column_stack(cols)


def _orth_or_dorth_pair(
    x: np.ndarray,
    y: np.ndarray,
    ref: np.ndarray,
    extra_basis: Optional[np.ndarray] = None,
    rank_tol: float = RANK_TOL,
    zero_tol: float = ZERO_ENERGY_TOL,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    When `extra_basis is None` (PRIMARY `r_orth_cal`): delegate verbatim to frozen
    `round10/src/metrics_r10.py::_orth_pair(x, y, ref)` so Rounds 10-16 share exact
    bit-level parity (`<= 1e-15`).

    When `extra_basis is not None` (SECONDARY `r_dorth_cal`, Systema-inspired shared-axis
    removal): project centered `(x, y)` orthogonal to `span(1, ref, extra_basis)` with
    rank and zero-energy guards so null-space predictors yield exact `0.0`.
    """
    xv = np.asarray(x, dtype=np.float64)
    yv = np.asarray(y, dtype=np.float64)
    if extra_basis is None:
        return _orth_pair(xv, yv, np.asarray(ref, dtype=np.float64))

    xc = xv - float(np.mean(xv))
    yc = yv - float(np.mean(yv))
    nx_raw = float(np.linalg.norm(xc))
    ny_raw = float(np.linalg.norm(yc))
    Q = _build_centered_basis(ref, extra_basis, rank_tol=rank_tol)
    if Q.shape[1] > 0:
        xc = xc - Q @ (Q.T @ xc)
        xc = xc - Q @ (Q.T @ xc)
        yc = yc - Q @ (Q.T @ yc)
        yc = yc - Q @ (Q.T @ yc)
    if float(np.linalg.norm(xc)) <= zero_tol * max(1.0, nx_raw):
        xc = np.zeros_like(xc)
    if float(np.linalg.norm(yc)) <= zero_tol * max(1.0, ny_raw):
        yc = np.zeros_like(yc)
    return xc, yc


def _collect_mode_moments(
    r_obs_full: np.ndarray,
    r_pred_list: Sequence[np.ndarray],
    t_mean_full: np.ndarray,
    rows: Sequence[int],
    idx_list: Sequence[np.ndarray],
    extra_basis_gene: Optional[np.ndarray],
    ha_res_full: Optional[np.ndarray],
    hb_res_full: Optional[np.ndarray],
    half_valid: Optional[np.ndarray],
) -> Tuple[List[Tuple[np.ndarray, ...]], Optional[Tuple[np.ndarray, ...]]]:
    """Collect per-perturbation residual moments for `r_orth_cal` or secondary `r_dorth_cal`."""
    k_preds = len(r_pred_list)
    xs_o: List[List[np.ndarray]] = [[] for _ in range(k_preds)]
    ys_o: List[List[np.ndarray]] = [[] for _ in range(k_preds)]
    ha_o: List[np.ndarray] = []
    hb_o: List[np.ndarray] = []
    has_halves = ha_res_full is not None and hb_res_full is not None
    any_hv = False

    for p, idx in zip(rows, idx_list):
        xv = r_obs_full[p, idx]
        rv = t_mean_full[p, idx]
        bv = extra_basis_gene[idx] if extra_basis_gene is not None else None
        for k, rp in enumerate(r_pred_list):
            xo, yo = _orth_or_dorth_pair(xv, rp[p, idx], rv, extra_basis=bv)
            xs_o[k].append(xo)
            ys_o[k].append(yo)

        if has_halves and ha_res_full is not None and hb_res_full is not None:
            hv_ok = half_valid is None or bool(half_valid[p])
            if hv_ok and np.all(np.isfinite(ha_res_full[p, idx])) and np.all(np.isfinite(hb_res_full[p, idx])):
                hao, hbo = _orth_or_dorth_pair(ha_res_full[p, idx], hb_res_full[p, idx], rv, extra_basis=bv)
                ha_o.append(hao)
                hb_o.append(hbo)
                any_hv = True
            else:
                z = np.zeros(0, dtype=np.float64)
                ha_o.append(z)
                hb_o.append(z)

    pred_moments = [_moments(xs_o[k], ys_o[k]) for k in range(k_preds)]
    half_moments = _moments(ha_o, hb_o) if any_hv else None
    return pred_moments, half_moments


def _moments_to_safe_corr(
    nv: Any,
    sx: Any,
    sy: Any,
    sxx: Any,
    syy: Any,
    sxy: Any,
    zero_null_to_zero: bool = False,
    var_floor: float = 1e-24,
) -> Any:
    """Compute pooled Pearson correlation; optionally map exact zero-variance null predictors to `0.0`."""
    _, r_raw, _, _, _ = _moments_to_corr_and_slope(nv, sx, sy, sxx, syy, sxy)
    if not zero_null_to_zero:
        return float(r_raw) if np.ndim(nv) == 0 else np.asarray(r_raw, dtype=np.float64)
    if np.ndim(nv) == 0:
        n_f = float(nv)
        if n_f <= 1.0:
            return 0.0
        vy = float(syy) - (float(sy) ** 2) / n_f
        vx = float(sxx) - (float(sx) ** 2) / n_f
        if vy <= var_floor or vx <= var_floor:
            return 0.0
        return float(r_raw)
    nv_a = np.asarray(nv, dtype=np.float64)
    sy_a = np.asarray(sy, dtype=np.float64)
    sx_a = np.asarray(sx, dtype=np.float64)
    syy_a = np.asarray(syy, dtype=np.float64)
    sxx_a = np.asarray(sxx, dtype=np.float64)
    r_arr = np.asarray(r_raw, dtype=np.float64).copy()
    with np.errstate(divide="ignore", invalid="ignore"):
        vy = np.where(nv_a > 1.0, syy_a - (sy_a ** 2) / np.maximum(nv_a, 1.0), 0.0)
        vx = np.where(nv_a > 1.0, sxx_a - (sx_a ** 2) / np.maximum(nv_a, 1.0), 0.0)
    r_arr[(nv_a <= 1.0) | (vy <= var_floor) | (vx <= var_floor)] = 0.0
    return r_arr


def _eval_boot_cal(
    pred_mom: Tuple[np.ndarray, ...],
    half_mom: Optional[Tuple[np.ndarray, ...]],
    draws: np.ndarray,
    single_half: bool = False,
    zero_null_to_zero: bool = False,
) -> Tuple[float, np.ndarray, float, float, float]:
    """Evaluate raw and disattenuated pooled correlation across joint bootstrap `draws`."""
    nv, sx, sy, sxx, syy, sxy = pred_mom
    pt_r = float(_moments_to_safe_corr(
        float(nv.sum()), float(sx.sum()), float(sy.sum()),
        float(sxx.sum()), float(syy.sum()), float(sxy.sum()),
        zero_null_to_zero=zero_null_to_zero,
    ))
    b_r = _moments_to_safe_corr(
        nv[draws].sum(axis=1), sx[draws].sum(axis=1), sy[draws].sum(axis=1),
        sxx[draws].sum(axis=1), syy[draws].sum(axis=1), sxy[draws].sum(axis=1),
        zero_null_to_zero=zero_null_to_zero,
    )

    n_boot = draws.shape[0]
    ceil_fn = single_half_ceiling if single_half else spearman_brown_ceiling
    if half_mom is not None:
        hnv, hsx, hsy, hsxx, hsyy, hsxy = half_mom
        _, pt_rho, _, _, _ = _moments_to_corr_and_slope(
            float(hnv.sum()), float(hsx.sum()), float(hsy.sum()),
            float(hsxx.sum()), float(hsyy.sum()), float(hsxy.sum()),
        )
        pt_ceil = float(ceil_fn(float(pt_rho)))
        _, b_rho, _, _, _ = _moments_to_corr_and_slope(
            hnv[draws].sum(axis=1), hsx[draws].sum(axis=1), hsy[draws].sum(axis=1),
            hsxx[draws].sum(axis=1), hsyy[draws].sum(axis=1), hsxy[draws].sum(axis=1),
        )
        b_ceil = np.array([ceil_fn(float(x)) for x in b_rho], dtype=np.float64)
    else:
        pt_rho = float("nan")
        pt_ceil = float("nan")
        b_ceil = np.full(n_boot, float("nan"), dtype=np.float64)

    def _safe_disatt(r_val: float, c_val: float) -> float:
        if zero_null_to_zero and r_val == 0.0 and np.isfinite(c_val) and c_val > 0.0:
            return 0.0
        return float(disattenuated_r(r_val, c_val))

    pt_cal = _safe_disatt(pt_r, pt_ceil)
    b_cal = np.array([_safe_disatt(float(r), float(c)) for r, c in zip(b_r, b_ceil)], dtype=np.float64)
    return pt_cal, b_cal, pt_r, float(pt_rho), pt_ceil


def compute_orth_suite_r16(
    delta_obs: np.ndarray,
    delta_pred: np.ndarray,
    delta_transfer: np.ndarray,
    target_gene_idx: Sequence[int],
    source_grand_mean: Optional[np.ndarray] = None,
    delta_half_a: Optional[np.ndarray] = None,
    delta_half_b: Optional[np.ndarray] = None,
    half_valid: Optional[np.ndarray] = None,
    gene_keep_mask: Optional[np.ndarray] = None,
    mode: str = "all_nontarget",
    top_k: int = 50,
    n_boot: int = 4000,
    seed: int = DEFAULT_R16_SEED,
) -> Dict[str, Any]:
    """
    Compute PRIMARY `r_orth_cal` (`source_grand_mean=None`) or SECONDARY Systema-inspired
    `r_dorth_cal` (`source_grand_mean` provided) on `all_nontarget` or `resp_top50`.
    """
    d_obs = np.asarray(delta_obs, dtype=np.float64)
    p_mat = np.asarray(delta_pred, dtype=np.float64)
    tm = np.asarray(delta_transfer, dtype=np.float64)
    b_s = np.asarray(source_grand_mean, dtype=np.float64) if source_grand_mean is not None else None

    rows, idx_list = _extract_masked_idx_per_pert(
        d_obs, target_gene_idx, gene_keep_mask=gene_keep_mask, mode=mode, top_k=top_k
    )
    m = len(rows)
    if m == 0:
        empty = summarize_boot_r16(float("nan"), np.array([float("nan")]), 0, n_boot)
        return {**empty, "raw_r": float("nan"), "rho_ab_orth": float("nan"), "ceiling": float("nan"), "mode": mode}

    r_obs = d_obs - tm
    r_pred = p_mat - tm
    ha_res = (np.asarray(delta_half_a, dtype=np.float64) - tm) if delta_half_a is not None else None
    hb_res = (np.asarray(delta_half_b, dtype=np.float64) - tm) if delta_half_b is not None else None

    [pm], hm = _collect_mode_moments(
        r_obs, [r_pred], tm, rows, idx_list, b_s, ha_res, hb_res, half_valid
    )
    draws = np.random.default_rng(seed).integers(0, m, size=(n_boot, m))
    pt_cal, b_cal, pt_raw, pt_rho, pt_ceil = _eval_boot_cal(
        pm, hm, draws, single_half=False, zero_null_to_zero=(b_s is not None)
    )
    out = summarize_boot_r16(pt_cal, b_cal, m, n_boot)
    out.update({
        "raw_r": pt_raw,
        "rho_ab_orth": pt_rho,
        "ceiling": pt_ceil,
        "mode": mode,
        "double_orth": b_s is not None,
    })
    return out


def compute_paired_contrast_r16(
    delta_obs: np.ndarray,
    delta_pred_a: np.ndarray,
    delta_pred_b: np.ndarray,
    delta_transfer: np.ndarray,
    target_gene_idx: Sequence[int],
    source_grand_mean: Optional[np.ndarray] = None,
    delta_half_a: Optional[np.ndarray] = None,
    delta_half_b: Optional[np.ndarray] = None,
    half_valid: Optional[np.ndarray] = None,
    gene_keep_mask: Optional[np.ndarray] = None,
    mode: str = "all_nontarget",
    top_k: int = 50,
    n_boot: int = 4000,
    seed: int = DEFAULT_R16_SEED,
) -> Dict[str, Any]:
    """
    Joint perturbation bootstrap of `r(pred_a) - r(pred_b)` sharing identical draws and
    noise ceiling (`H-I1` matched vs mismatched, `H-I2` `STATE_SE` vs `FCT`).
    """
    d_obs = np.asarray(delta_obs, dtype=np.float64)
    pa = np.asarray(delta_pred_a, dtype=np.float64)
    pb = np.asarray(delta_pred_b, dtype=np.float64)
    tm = np.asarray(delta_transfer, dtype=np.float64)
    b_s = np.asarray(source_grand_mean, dtype=np.float64) if source_grand_mean is not None else None

    rows, idx_list = _extract_masked_idx_per_pert(
        d_obs, target_gene_idx, gene_keep_mask=gene_keep_mask, mode=mode, top_k=top_k
    )
    m = len(rows)
    if m == 0:
        empty = summarize_boot_r16(float("nan"), np.array([float("nan")]), 0, n_boot)
        return {**empty, "a_summary": empty, "b_summary": empty}

    r_obs = d_obs - tm
    ha_res = (np.asarray(delta_half_a, dtype=np.float64) - tm) if delta_half_a is not None else None
    hb_res = (np.asarray(delta_half_b, dtype=np.float64) - tm) if delta_half_b is not None else None

    [pm_a, pm_b], hm = _collect_mode_moments(
        r_obs, [pa - tm, pb - tm], tm, rows, idx_list, b_s, ha_res, hb_res, half_valid
    )
    draws = np.random.default_rng(seed).integers(0, m, size=(n_boot, m))
    pt_a, b_a, _, pt_rho, pt_ceil = _eval_boot_cal(
        pm_a, hm, draws, single_half=False, zero_null_to_zero=(b_s is not None)
    )
    pt_b, b_b, _, _, _ = _eval_boot_cal(
        pm_b, hm, draws, single_half=False, zero_null_to_zero=(b_s is not None)
    )

    mad, exact_inv = _array_max_abs_diff_and_exact_inv(pa, pb)
    diff_sum = summarize_boot_r16(float(pt_a - pt_b), b_a - b_b, m, n_boot)
    diff_sum["a_summary"] = summarize_boot_r16(pt_a, b_a, m, n_boot)
    diff_sum["b_summary"] = summarize_boot_r16(pt_b, b_b, m, n_boot)
    diff_sum["rho_ab_orth"] = pt_rho
    diff_sum["ceiling"] = pt_ceil
    diff_sum["double_orth"] = b_s is not None
    diff_sum["max_abs_diff"] = mad
    diff_sum["exact_invariance"] = exact_inv
    return diff_sum


def compute_paired_model_vs_all_mismatched_r16(
    delta_obs: np.ndarray,
    delta_pred_matched: np.ndarray,
    mismatched_preds: Dict[str, np.ndarray],
    delta_transfer: np.ndarray,
    target_gene_idx: Sequence[int],
    source_grand_mean: Optional[np.ndarray] = None,
    delta_half_a: Optional[np.ndarray] = None,
    delta_half_b: Optional[np.ndarray] = None,
    half_valid: Optional[np.ndarray] = None,
    gene_keep_mask: Optional[np.ndarray] = None,
    n_boot: int = 4000,
    seed: int = DEFAULT_R16_SEED,
) -> Dict[str, Any]:
    """
    `PREREG_I H-I1`: Compute joint perturbation bootstrap paired differences
    `r(matched) - r(mismatched[s'])` across EVERY mismatched control line `s'`
    sharing identical bootstrap draws, returning `97.5%` and `95%` (TOST) CIs.
    """
    d_obs = np.asarray(delta_obs, dtype=np.float64)
    pm_mat = np.asarray(delta_pred_matched, dtype=np.float64)
    tm = np.asarray(delta_transfer, dtype=np.float64)
    b_s = np.asarray(source_grand_mean, dtype=np.float64) if source_grand_mean is not None else None
    src_keys = sorted(mismatched_preds.keys())

    rows, idx_list = _extract_masked_idx_per_pert(
        d_obs, target_gene_idx, gene_keep_mask=gene_keep_mask, mode="all_nontarget", top_k=50
    )
    m = len(rows)
    if m == 0 or not src_keys:
        return {"n_rows": m, "per_source": {}, "double_orth": b_s is not None}

    r_obs = d_obs - tm
    pred_res_list = [pm_mat - tm] + [np.asarray(mismatched_preds[k], dtype=np.float64) - tm for k in src_keys]
    ha_res = (np.asarray(delta_half_a, dtype=np.float64) - tm) if delta_half_a is not None else None
    hb_res = (np.asarray(delta_half_b, dtype=np.float64) - tm) if delta_half_b is not None else None

    moms, hm = _collect_mode_moments(
        r_obs, pred_res_list, tm, rows, idx_list, b_s, ha_res, hb_res, half_valid
    )
    draws = np.random.default_rng(seed).integers(0, m, size=(n_boot, m))
    pt_m, b_m, _, pt_rho, pt_ceil = _eval_boot_cal(
        moms[0], hm, draws, single_half=False, zero_null_to_zero=(b_s is not None)
    )

    per_source: Dict[str, Dict[str, Any]] = {}
    for j, sk in enumerate(src_keys, start=1):
        pt_s, b_s_boot, _, _, _ = _eval_boot_cal(
            moms[j], hm, draws, single_half=False, zero_null_to_zero=(b_s is not None)
        )
        mad, exact_inv = _array_max_abs_diff_and_exact_inv(pm_mat, mismatched_preds[sk])
        d_sum = summarize_boot_r16(float(pt_m - pt_s), b_m - b_s_boot, m, n_boot)
        d_sum["mismatched_estimate"] = pt_s
        d_sum["max_abs_diff"] = mad
        d_sum["exact_invariance"] = exact_inv
        per_source[sk] = d_sum

    return {
        "n_rows": m,
        "matched_summary": summarize_boot_r16(pt_m, b_m, m, n_boot),
        "rho_ab_orth": pt_rho,
        "ceiling": pt_ceil,
        "double_orth": b_s is not None,
        "per_source": per_source,
    }


def compute_paired_with_within_line_resample_r16(
    delta_obs: np.ndarray,
    delta_pred_matched: np.ndarray,
    delta_pred_resample: np.ndarray,
    mismatched_preds: Dict[str, np.ndarray],
    delta_transfer: np.ndarray,
    target_gene_idx: Sequence[int],
    source_grand_mean: Optional[np.ndarray] = None,
    delta_half_a: Optional[np.ndarray] = None,
    delta_half_b: Optional[np.ndarray] = None,
    half_valid: Optional[np.ndarray] = None,
    gene_keep_mask: Optional[np.ndarray] = None,
    n_boot: int = 4000,
    seed: int = DEFAULT_R16_SEED,
) -> Dict[str, Any]:
    """
    `PREREG_I Addendum A Item 8` (secondary descriptive within-line resample reference):
    Using the exact same joint perturbation bootstrap draws (`seed`, `n_boot`) and noise ceiling
    as primary `H-I1`, compute:
      - `delta_pair_resample = r(matched) - r(resample)` with `97.5%` and `95.0%` CIs;
      - For every mismatched source `s'`:
          `delta_pair_mis = r(matched) - r(mis[s'])` and
          `delta_pair_excess_over_resample = delta_pair(s') - delta_pair(resample)`
        computed on the exact same bootstrap draws.
    """
    d_obs = np.asarray(delta_obs, dtype=np.float64)
    pm_mat = np.asarray(delta_pred_matched, dtype=np.float64)
    pr_mat = np.asarray(delta_pred_resample, dtype=np.float64)
    tm = np.asarray(delta_transfer, dtype=np.float64)
    b_s = np.asarray(source_grand_mean, dtype=np.float64) if source_grand_mean is not None else None
    src_keys = sorted(mismatched_preds.keys())

    rows, idx_list = _extract_masked_idx_per_pert(
        d_obs, target_gene_idx, gene_keep_mask=gene_keep_mask, mode="all_nontarget", top_k=50
    )
    m = len(rows)
    if m == 0:
        empty = summarize_boot_r16(float("nan"), np.array([float("nan")]), 0, n_boot)
        return {
            "n_rows": 0,
            "matched_summary": empty,
            "resample_summary": empty,
            "delta_pair_resample": empty,
            "rho_ab_orth": float("nan"),
            "ceiling": float("nan"),
            "double_orth": b_s is not None,
            "per_source": {},
        }

    r_obs = d_obs - tm
    pred_res_list = [
        pm_mat - tm,
        pr_mat - tm,
    ] + [np.asarray(mismatched_preds[k], dtype=np.float64) - tm for k in src_keys]
    ha_res = (np.asarray(delta_half_a, dtype=np.float64) - tm) if delta_half_a is not None else None
    hb_res = (np.asarray(delta_half_b, dtype=np.float64) - tm) if delta_half_b is not None else None

    moms, hm = _collect_mode_moments(
        r_obs, pred_res_list, tm, rows, idx_list, b_s, ha_res, hb_res, half_valid
    )
    draws = np.random.default_rng(seed).integers(0, m, size=(n_boot, m))
    pt_m, b_m, _, pt_rho, pt_ceil = _eval_boot_cal(
        moms[0], hm, draws, single_half=False, zero_null_to_zero=(b_s is not None)
    )
    pt_r, b_r, _, _, _ = _eval_boot_cal(
        moms[1], hm, draws, single_half=False, zero_null_to_zero=(b_s is not None)
    )

    mad_r, exact_inv_r = _array_max_abs_diff_and_exact_inv(pm_mat, pr_mat)
    pt_dp_resamp = float(pt_m - pt_r)
    b_dp_resamp = b_m - b_r
    dp_resamp_sum = summarize_boot_r16(pt_dp_resamp, b_dp_resamp, m, n_boot)
    dp_resamp_sum["resample_estimate"] = float(pt_r)
    dp_resamp_sum["max_abs_diff"] = mad_r
    dp_resamp_sum["exact_invariance"] = exact_inv_r

    per_source: Dict[str, Dict[str, Any]] = {}
    for j, sk in enumerate(src_keys, start=2):
        pt_s, b_s_boot, _, _, _ = _eval_boot_cal(
            moms[j], hm, draws, single_half=False, zero_null_to_zero=(b_s is not None)
        )
        mad_s, exact_inv_s = _array_max_abs_diff_and_exact_inv(pm_mat, mismatched_preds[sk])
        pt_dp_mis = float(pt_m - pt_s)
        b_dp_mis = b_m - b_s_boot
        dp_mis_sum = summarize_boot_r16(pt_dp_mis, b_dp_mis, m, n_boot)
        dp_mis_sum["mismatched_estimate"] = float(pt_s)
        dp_mis_sum["max_abs_diff"] = mad_s
        dp_mis_sum["exact_invariance"] = exact_inv_s

        pt_excess = float(pt_dp_mis - pt_dp_resamp)
        b_excess = b_dp_mis - b_dp_resamp
        excess_sum = summarize_boot_r16(pt_excess, b_excess, m, n_boot)
        per_source[sk] = {
            "delta_pair_mis": dp_mis_sum,
            "delta_pair_excess_over_resample": excess_sum,
        }

    return {
        "n_rows": m,
        "matched_summary": summarize_boot_r16(pt_m, b_m, m, n_boot),
        "resample_summary": summarize_boot_r16(pt_r, b_r, m, n_boot),
        "delta_pair_resample": dp_resamp_sum,
        "rho_ab_orth": pt_rho,
        "ceiling": pt_ceil,
        "double_orth": b_s is not None,
        "per_source": per_source,
    }


def compute_gene_label_perm_null_r16(
    delta_obs: np.ndarray,
    delta_pred: np.ndarray,
    delta_transfer: np.ndarray,
    target_gene_idx: Sequence[int],
    source_grand_mean: Optional[np.ndarray] = None,
    delta_half_a: Optional[np.ndarray] = None,
    delta_half_b: Optional[np.ndarray] = None,
    half_valid: Optional[np.ndarray] = None,
    gene_keep_mask: Optional[np.ndarray] = None,
    mode: str = "all_nontarget",
    top_k: int = 50,
    n_perm: int = 100,
    perm_seed_offset: int = 0,
) -> Dict[str, Any]:
    """
    `PREREG_I Addendum A Item 4` (`ST-HVG` `Perm100` null, following frozen `PREREG_E §4.1`):
    Permutes the gene columns of `delta_pred` within `G_eval` (`G = delta_pred.shape[1]`)
    using `np.random.default_rng(perm_seed_offset + k).permutation(G)` for `k = 0..n_perm-1`
    (default `perm_seed_offset = 0`, i.e., `seed = k` verbatim per `PREREG_E §4.1`), holding
    `delta_obs`, `delta_transfer`, `target_gene_idx`, and split halves fixed.
    """
    d_obs = np.asarray(delta_obs, dtype=np.float64)
    p_mat = np.asarray(delta_pred, dtype=np.float64)
    tm = np.asarray(delta_transfer, dtype=np.float64)
    b_s = np.asarray(source_grand_mean, dtype=np.float64) if source_grand_mean is not None else None

    if d_obs.shape != p_mat.shape or d_obs.ndim != 2:
        raise ValueError(f"Shape mismatch: delta_obs={d_obs.shape}, delta_pred={p_mat.shape}")
    G = int(d_obs.shape[1])

    rows, idx_list = _extract_masked_idx_per_pert(
        d_obs, target_gene_idx, gene_keep_mask=gene_keep_mask, mode=mode, top_k=top_k
    )
    m = len(rows)
    if m == 0 or G == 0 or int(n_perm) <= 0:
        return {
            "n_perm": int(n_perm),
            "perm_seed_offset": int(perm_seed_offset),
            "perm_seed_rule": f"np.random.default_rng({int(perm_seed_offset)} + k).permutation(G) for k=0..{max(int(n_perm) - 1, 0)}",
            "n_genes_permuted": G,
            "double_orth": b_s is not None,
            "observed_estimate": float("nan"),
            "rho_ab_orth": float("nan"),
            "ceiling": float("nan"),
            "perm_estimates": [],
            "n_finite_perms": 0,
            "median": float("nan"),
            "mean": float("nan"),
            "q025": float("nan"),
            "q975": float("nan"),
            "max": float("nan"),
            "empirical_p_ge": float("nan"),
        }

    r_obs = d_obs - tm
    ha_res = (np.asarray(delta_half_a, dtype=np.float64) - tm) if delta_half_a is not None else None
    hb_res = (np.asarray(delta_half_b, dtype=np.float64) - tm) if delta_half_b is not None else None
    zero_null = b_s is not None

    # Compute observed pred moments + split-half ceiling once (shared across all k permutations)
    [pm_obs], hm = _collect_mode_moments(
        r_obs, [p_mat - tm], tm, rows, idx_list, b_s, ha_res, hb_res, half_valid
    )
    if hm is not None:
        hnv, hsx, hsy, hsxx, hsyy, hsxy = hm
        _, pt_rho, _, _, _ = _moments_to_corr_and_slope(
            float(hnv.sum()), float(hsx.sum()), float(hsy.sum()),
            float(hsxx.sum()), float(hsyy.sum()), float(hsxy.sum()),
        )
        pt_rho = float(pt_rho)
        pt_ceil = float(spearman_brown_ceiling(pt_rho))
    else:
        pt_rho = float("nan")
        pt_ceil = float("nan")

    def _eval_point_from_mom(pm: Tuple[np.ndarray, ...]) -> float:
        nv, sx, sy, sxx, syy, sxy = pm
        pt_r = float(
            _moments_to_safe_corr(
                float(nv.sum()), float(sx.sum()), float(sy.sum()),
                float(sxx.sum()), float(syy.sum()), float(sxy.sum()),
                zero_null_to_zero=zero_null,
            )
        )
        if zero_null and pt_r == 0.0 and np.isfinite(pt_ceil) and pt_ceil > 0.0:
            return 0.0
        return float(disattenuated_r(pt_r, pt_ceil))

    obs_est = _eval_point_from_mom(pm_obs)
    perm_vals = np.empty(int(n_perm), dtype=np.float64)
    for k in range(int(n_perm)):
        perm = np.random.default_rng(int(perm_seed_offset) + k).permutation(G)
        p_perm = p_mat[:, perm]
        [pm_k], _ = _collect_mode_moments(
            r_obs, [p_perm - tm], tm, rows, idx_list, b_s, None, None, None
        )
        perm_vals[k] = _eval_point_from_mom(pm_k)

    fin_vals = perm_vals[np.isfinite(perm_vals)]
    emp_p = (
        float((1 + int(np.sum(fin_vals >= obs_est))) / (int(fin_vals.size) + 1))
        if (fin_vals.size > 0 and np.isfinite(obs_est))
        else float("nan")
    )
    return {
        "n_perm": int(n_perm),
        "perm_seed_offset": int(perm_seed_offset),
        "perm_seed_rule": f"np.random.default_rng({int(perm_seed_offset)} + k).permutation(G) for k=0..{int(n_perm) - 1}",
        "n_genes_permuted": G,
        "double_orth": zero_null,
        "observed_estimate": float(obs_est),
        "rho_ab_orth": float(pt_rho),
        "ceiling": float(pt_ceil),
        "perm_estimates": [float(x) for x in perm_vals],
        "n_finite_perms": int(fin_vals.size),
        "median": float(np.median(fin_vals)) if fin_vals.size > 0 else float("nan"),
        "mean": float(np.mean(fin_vals)) if fin_vals.size > 0 else float("nan"),
        "q025": float(np.percentile(fin_vals, 2.5)) if fin_vals.size > 0 else float("nan"),
        "q975": float(np.percentile(fin_vals, 97.5)) if fin_vals.size > 0 else float("nan"),
        "max": float(np.max(fin_vals)) if fin_vals.size > 0 else float("nan"),
        "empirical_p_ge": emp_p,
    }


def compute_joint_cross_mask_r16(
    delta_obs: np.ndarray,
    delta_pred: np.ndarray,
    delta_transfer: np.ndarray,
    target_gene_idx: Sequence[int],
    source_grand_mean: Optional[np.ndarray] = None,
    delta_half_a: Optional[np.ndarray] = None,
    delta_half_b: Optional[np.ndarray] = None,
    half_valid: Optional[np.ndarray] = None,
    gene_keep_mask: Optional[np.ndarray] = None,
    top_k: int = 50,
    n_boot: int = 4000,
    seed: int = DEFAULT_R16_SEED,
) -> Dict[str, Any]:
    """
    Secondary `§4.2` cross-mask bootstrap:
        `r(delta_pred, resp_top50_split) - r(delta_pred, all_nontarget)`
    restricted to `rows_common`. When `source_grand_mean is None`, matches frozen
    `round12/src/metrics_r12.py::compute_joint_sa8_split_minus_all_r12` to `<= 1e-15`.
    """
    empty = summarize_boot_r16(float("nan"), np.array([float("nan")]), 0, n_boot)
    if delta_half_a is None or delta_half_b is None:
        return {**empty, "n_common_perts": 0, "split_summary": empty, "all_common_summary": empty}

    d_obs = np.asarray(delta_obs, dtype=np.float64)
    p_mat = np.asarray(delta_pred, dtype=np.float64)
    tm = np.asarray(delta_transfer, dtype=np.float64)
    b_s = np.asarray(source_grand_mean, dtype=np.float64) if source_grand_mean is not None else None
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

    [pm_ab], hm_ab = _collect_mode_moments(
        hb_res, [pred_res], tm, rows_common, idx_ab_c, b_s, ha_res, hb_res, half_valid
    )
    [pm_ba], hm_ba = _collect_mode_moments(
        ha_res, [pred_res], tm, rows_common, idx_ba_c, b_s, hb_res, ha_res, half_valid
    )
    [pm_all], hm_all = _collect_mode_moments(
        r_obs_full, [pred_res], tm, rows_common, idx_all_c, b_s, ha_res, hb_res, half_valid
    )
    if hm_ab is None or hm_ba is None or hm_all is None:
        return {**empty, "n_common_perts": m, "split_summary": empty, "all_common_summary": empty}

    draws = np.random.default_rng(seed).integers(0, m, size=(n_boot, m))
    pt_ab, b_ab, _, _, _ = _eval_boot_cal(
        pm_ab, hm_ab, draws, single_half=True, zero_null_to_zero=(b_s is not None)
    )
    pt_ba, b_ba, _, _, _ = _eval_boot_cal(
        pm_ba, hm_ba, draws, single_half=True, zero_null_to_zero=(b_s is not None)
    )
    pt_split = 0.5 * (pt_ab + pt_ba)
    b_split = 0.5 * (b_ab + b_ba)

    pt_all, b_all, _, _, _ = _eval_boot_cal(
        pm_all, hm_all, draws, single_half=False, zero_null_to_zero=(b_s is not None)
    )
    diff_sum = summarize_boot_r16(float(pt_split - pt_all), b_split - b_all, m, n_boot)
    diff_sum["n_common_perts"] = int(m)
    diff_sum["double_orth"] = b_s is not None
    diff_sum["split_summary"] = summarize_boot_r16(pt_split, b_split, m, n_boot)
    diff_sum["all_common_summary"] = summarize_boot_r16(pt_all, b_all, m, n_boot)
    return diff_sum


def compute_winner_curse_weight_shift_r16(
    mu_target: np.ndarray,
    mu_sources: np.ndarray,
    delta_half_a: np.ndarray,
    delta_half_b: np.ndarray,
    delta_transfer: np.ndarray,
    target_gene_idx: Sequence[int],
    half_valid: Optional[np.ndarray] = None,
    gene_keep_mask: Optional[np.ndarray] = None,
    top_k: int = 50,
    n_boot: int = 4000,
    seed: int = DEFAULT_R16_SEED,
) -> Dict[str, Any]:
    """
    Secondary descriptive diagnostic (`§4.2`, Mechanism Hypothesis M2): perturbation-level
    target-over-source saturation weight shift `w(mu_t) - w(mu_s)` (`w(mu) = mu/(1+mu)`)
    on independent Half A / Half B top-50 genes minus `all_nontarget`.
    """
    mu_t = np.maximum(np.asarray(mu_target, dtype=np.float64), 0.0)
    mu_s = np.maximum(np.asarray(mu_sources, dtype=np.float64), 0.0)
    w_s = np.mean(mu_s / (1.0 + mu_s), axis=0) if mu_s.ndim == 2 else (mu_s / (1.0 + mu_s))
    w_t = mu_t / (1.0 + mu_t)
    dw = w_t - w_s

    ha = np.asarray(delta_half_a, dtype=np.float64)
    hb = np.asarray(delta_half_b, dtype=np.float64)
    tm = np.asarray(delta_transfer, dtype=np.float64)

    rows_ab, idx_ab = _prepare_split_direction_indices(
        ha, hb, ha, tm, target_gene_idx, half_valid, gene_keep_mask=gene_keep_mask, top_k=top_k
    )
    rows_ba, idx_ba = _prepare_split_direction_indices(
        hb, ha, hb, tm, target_gene_idx, half_valid, gene_keep_mask=gene_keep_mask, top_k=top_k
    )
    rows_all, idx_all = _extract_masked_idx_per_pert(
        0.5 * (ha + hb), target_gene_idx, gene_keep_mask=gene_keep_mask, mode="all_nontarget", top_k=top_k
    )
    map_ab = dict(zip(rows_ab, idx_ab))
    map_ba = dict(zip(rows_ba, idx_ba))
    map_all = dict(zip(rows_all, idx_all))
    rows_common = [p for p in rows_ab if p in map_ba and p in map_all]
    m = len(rows_common)
    if m == 0:
        empty = summarize_boot_r16(float("nan"), np.array([float("nan")]), 0, n_boot)
        return {**empty, "w_t_all_mean": float("nan"), "w_t_top50_mean": float("nan")}

    dw_shifts = np.empty(m, dtype=np.float64)
    wt_top = np.empty(m, dtype=np.float64)
    wt_all = np.empty(m, dtype=np.float64)
    for i, p in enumerate(rows_common):
        ia, ib, il = map_ab[p], map_ba[p], map_all[p]
        dw_shifts[i] = 0.5 * (float(np.mean(dw[ia])) + float(np.mean(dw[ib]))) - float(np.mean(dw[il]))
        wt_top[i] = 0.5 * (float(np.mean(w_t[ia])) + float(np.mean(w_t[ib])))
        wt_all[i] = float(np.mean(w_t[il]))

    draws = np.random.default_rng(seed).integers(0, m, size=(n_boot, m))
    boot_means = dw_shifts[draws].mean(axis=1)
    out = summarize_boot_r16(float(np.mean(dw_shifts)), boot_means, m, n_boot)
    out["w_t_all_mean"] = float(np.mean(wt_all))
    out["w_t_top50_mean"] = float(np.mean(wt_top))
    out["dw_all_mean"] = float(np.mean(dw))
    return out
