#!/usr/bin/env python3
"""
Clean-Room Round 16 Endpoint Evaluation & Joint Perturbation Bootstrap
(`A-T3`, strictly following `PREREG_I.md v4 §1.1, §3, §4`, `seed = 20261601`).

Implements:
  - Primary `r_orth_cal` (`all_nontarget` default, plus `resp_top50` & `resp_top50_split`)
  - Secondary `r_dorth_cal` (Systema-inspired projection onto `span{1, t_mean_p, b_S}^perp`)
    preserving exact Property P1 null invariance (`0.0 <= 1e-14`) by construction
  - Orthogonalized split-half correlation explicitly keyed `rho_ab_orth` (`M1`)
  - Joint perturbation-level bootstrap (`n_boot = 4,000`, `97.5%` & `95.0%` CIs, `valid_frac`)
  - Paired contrast bootstrap `(pred_a - pred_b)` sharing identical draws & noise ceiling
  - `M2` descriptive saturation-weight shift diagnostic

Never reads or imports `round16/src/`.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

EPS = 1e-12
NULL_NORM_TOL = 1e-12
DEFAULT_R16_SEED = 20261601


def summarize_boot_r16(
    point: float,
    boot_arr: Sequence[float] | np.ndarray,
    n_items: int,
    n_boot: int,
) -> Dict[str, Any]:
    """Summarize bootstrap draws into two-sided 95.0% (`H-I1` TOST) and 97.5% (`H-I1/I2` sup.) CIs."""
    arr = np.asarray(boot_arr, dtype=np.float64)
    fin_mask = np.isfinite(arr)
    valid_frac = float(np.mean(fin_mask)) if arr.size > 0 else 0.0
    fin = arr[fin_mask]
    if fin.size == 0 or not np.isfinite(point):
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
    q0125, q025, q975, q9875 = np.percentile(fin, [1.25, 2.5, 97.5, 98.75])
    return {
        "estimate": float(point),
        "ci95_lo": float(q025),
        "ci95_hi": float(q975),
        "ci975_lo": float(q0125),
        "ci975_hi": float(q9875),
        "boot_se": float(np.std(fin, ddof=1)) if fin.size > 1 else 0.0,
        "valid_frac": valid_frac,
        "n_items": int(n_items),
        "n_boot": int(n_boot),
    }


def spearman_brown_ceiling(rho_ab: float) -> float:
    """Spearman-Brown split-half ceiling `sqrt(2*rho / (1 + rho))`; `NaN` if `rho <= 0`."""
    if not np.isfinite(rho_ab) or rho_ab <= 0.0:
        return float("nan")
    sb = (2.0 * float(rho_ab)) / (1.0 + float(rho_ab))
    return float(np.sqrt(max(sb, 0.0)))


def disattenuate(r_raw: float, ceil: float) -> float:
    if not (np.isfinite(r_raw) and np.isfinite(ceil)) or ceil <= 1e-8:
        return float("nan")
    return float(np.clip(float(r_raw) / float(ceil), -1.0, 1.0))


def _orth_pair(x: np.ndarray, y: np.ndarray, ref: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Uncentered orthogonalization perpendicular to `ref` matching frozen `metrics_r10::_orth_pair`."""
    nrm = float(np.linalg.norm(ref))
    if nrm <= EPS:
        return x.copy(), y.copy()
    b = ref / nrm
    return x - float(np.dot(x, b)) * b, y - float(np.dot(y, b)) * b


def _dorth_basis(t_mean_row: np.ndarray, grand_mean_row: np.ndarray, rcond: float = 1e-12) -> np.ndarray:
    """Orthonormal basis spanning `span{1, t_mean_p, b_S}` over masked gene indices."""
    m = t_mean_row.shape[0]
    ones = np.full(m, 1.0 / math.sqrt(float(m)), dtype=np.float64)
    extra_cols: List[np.ndarray] = []
    for raw_v in (t_mean_row, grand_mean_row):
        v = np.asarray(raw_v, dtype=np.float64)
        vp = v - float(np.dot(ones, v)) * ones
        nrm = float(np.linalg.norm(vp))
        if nrm > EPS:
            extra_cols.append(vp / nrm)
    if not extra_cols:
        return ones[:, None]
    mat = np.column_stack(extra_cols)
    u, s, _ = np.linalg.svd(mat, full_matrices=False)
    keep = s > max(EPS, rcond * float(s[0]))
    if not np.any(keep):
        return ones[:, None]
    return np.column_stack([ones, u[:, keep]])


def _project_perp(vec: np.ndarray, q: np.ndarray, ref_scale: float = 1.0) -> np.ndarray:
    proj = vec - q @ (q.T @ vec)
    if float(np.linalg.norm(proj)) <= NULL_NORM_TOL * max(1.0, ref_scale):
        return np.zeros_like(proj)
    return proj


def _extract_pert_gene_indices(
    delta_obs: np.ndarray,
    t_mean: np.ndarray,
    pred: np.ndarray,
    target_cols: Sequence[int],
    mode: str = "all_nontarget",
    top_k: int = 50,
    extra_pred: Optional[np.ndarray] = None,
    source_grand_mean: Optional[np.ndarray] = None,
    gene_keep_mask: Optional[np.ndarray] = None,
) -> List[Tuple[int, np.ndarray]]:
    n_t, n_g = delta_obs.shape
    kept: List[Tuple[int, np.ndarray]] = []
    gm_valid = np.isfinite(source_grand_mean) if source_grand_mean is not None else None
    for p in range(n_t):
        valid = np.isfinite(delta_obs[p]) & np.isfinite(t_mean[p]) & np.isfinite(pred[p])
        if extra_pred is not None:
            valid &= np.isfinite(extra_pred[p])
        if gm_valid is not None:
            valid &= gm_valid
        if gene_keep_mask is not None:
            valid &= gene_keep_mask
        tc = int(target_cols[p])
        if 0 <= tc < n_g:
            valid[tc] = False
        idx = np.flatnonzero(valid)
        if idx.size < 2:
            continue
        if mode == "resp_top50" and idx.size > top_k:
            idx = idx[np.argsort(np.abs(delta_obs[p, idx]), kind="mergesort")[-top_k:]]
        kept.append((p, idx))
    return kept


def _moments(xs: Sequence[np.ndarray], ys: Sequence[np.ndarray]) -> Tuple[np.ndarray, ...]:
    nv = np.asarray([float(x.size) for x in xs], dtype=np.float64)
    sx = np.asarray([float(np.sum(x)) for x in xs], dtype=np.float64)
    sy = np.asarray([float(np.sum(y)) for y in ys], dtype=np.float64)
    sxx = np.asarray([float(np.dot(x, x)) for x in xs], dtype=np.float64)
    syy = np.asarray([float(np.dot(y, y)) for y in ys], dtype=np.float64)
    sxy = np.asarray([float(np.dot(x, y)) for x, y in zip(xs, ys)], dtype=np.float64)
    return nv, sx, sy, sxx, syy, sxy


def _moments_to_corr_and_slope(
    nv: np.ndarray | float,
    sx: np.ndarray | float,
    sy: np.ndarray | float,
    sxx: np.ndarray | float,
    syy: np.ndarray | float,
    sxy: np.ndarray | float,
    draws: Optional[np.ndarray] = None,
    zero_degenerate_y_as_zero: bool = False,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    if draws is None:
        N = np.asarray(float(np.sum(nv)))
        Sx = np.asarray(float(np.sum(sx)))
        Sy = np.asarray(float(np.sum(sy)))
        Sxx = np.asarray(float(np.sum(sxx)))
        Syy = np.asarray(float(np.sum(syy)))
        Sxy = np.asarray(float(np.sum(sxy)))
    else:
        N = np.sum(nv[draws], axis=1)
        Sx = np.sum(sx[draws], axis=1)
        Sy = np.sum(sy[draws], axis=1)
        Sxx = np.sum(sxx[draws], axis=1)
        Syy = np.sum(syy[draws], axis=1)
        Sxy = np.sum(sxy[draws], axis=1)

    n_safe = np.maximum(N, 1.0)
    cov = Sxy - (Sx * Sy) / n_safe
    vx = Sxx - (Sx * Sx) / n_safe
    vy = Syy - (Sy * Sy) / n_safe
    denom = np.sqrt(np.maximum(vx, 0.0) * np.maximum(vy, 0.0))
    fallback = 0.0 if zero_degenerate_y_as_zero else np.nan
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = np.where(
            denom > EPS,
            np.clip(cov / denom, -1.0, 1.0),
            np.where(vx > EPS, fallback, np.nan),
        )
        slope_ols = np.where(vx > EPS, cov / vx, np.nan)
        sd_ratio = np.where(vx > EPS, np.sqrt(np.maximum(vy, 0.0) / np.maximum(vx, EPS)), np.nan)
    return corr, slope_ols, sd_ratio


def _project_streams(
    kept: Sequence[Tuple[int, np.ndarray]],
    delta_obs: np.ndarray,
    pred: np.ndarray,
    t_mean: np.ndarray,
    source_grand_mean: Optional[np.ndarray],
    half_a_obs: Optional[np.ndarray],
    half_b_obs: Optional[np.ndarray],
    half_valid: Optional[np.ndarray],
) -> Dict[str, Any]:
    r_obs_full = np.asarray(delta_obs, dtype=np.float64) - np.asarray(t_mean, dtype=np.float64)
    r_pred_full = np.asarray(pred, dtype=np.float64) - np.asarray(t_mean, dtype=np.float64)
    t_full = np.asarray(t_mean, dtype=np.float64)
    has_halves = half_a_obs is not None and half_b_obs is not None
    if has_halves:
        ha_full = np.asarray(half_a_obs, dtype=np.float64) - t_full
        hb_full = np.asarray(half_b_obs, dtype=np.float64) - t_full
    else:
        ha_full = hb_full = None

    xo_raw, yo_raw, ha_raw, hb_raw = [], [], [], []
    xo_proj, yo_proj, ha_proj, hb_proj = [], [], [], []
    cos_vals: List[float] = []
    de50_vals: List[float] = []
    any_hv = False

    for p, idx in kept:
        obs_p = delta_obs[p, idx].astype(np.float64)
        pred_p = pred[p, idx].astype(np.float64)
        den_cos = float(np.linalg.norm(obs_p) * np.linalg.norm(pred_p))
        if den_cos > EPS:
            cos_vals.append(float(np.clip(np.dot(obs_p, pred_p) / den_cos, -1.0, 1.0)))

        k_de = min(50, idx.size)
        top_sub = np.argsort(np.abs(obs_p - t_full[p, idx]), kind="mergesort")[-k_de:]
        dt_de = obs_p[top_sub] - np.mean(obs_p[top_sub])
        dp_de = pred_p[top_sub] - np.mean(pred_p[top_sub])
        den_de = float(np.linalg.norm(dt_de) * np.linalg.norm(dp_de))
        if den_de > EPS:
            de50_vals.append(float(np.clip(np.dot(dt_de, dp_de) / den_de, -1.0, 1.0)))

        xv = r_obs_full[p, idx]
        yv = r_pred_full[p, idx]
        rv = t_full[p, idx]
        xo_raw.append(xv)
        yo_raw.append(yv)

        if source_grand_mean is None:
            xp, yp = _orth_pair(xv, yv, rv)
            q_mat = None
        else:
            gm = np.asarray(source_grand_mean[idx], dtype=np.float64)
            q_mat = _dorth_basis(rv, gm)
            xp = _project_perp(xv, q_mat, ref_scale=float(np.linalg.norm(xv)))
            yp = _project_perp(yv, q_mat, ref_scale=float(np.linalg.norm(yv)))
        xo_proj.append(xp)
        yo_proj.append(yp)

        if has_halves and ha_full is not None and hb_full is not None:
            hv_ok = half_valid is None or bool(half_valid[p])
            if hv_ok and np.all(np.isfinite(ha_full[p, idx])) and np.all(np.isfinite(hb_full[p, idx])):
                hav = ha_full[p, idx]
                hbv = hb_full[p, idx]
                ha_raw.append(hav)
                hb_raw.append(hbv)
                if source_grand_mean is None:
                    hap, hbp = _orth_pair(hav, hbv, rv)
                else:
                    assert q_mat is not None
                    hap = _project_perp(hav, q_mat, ref_scale=float(np.linalg.norm(hav)))
                    hbp = _project_perp(hbv, q_mat, ref_scale=float(np.linalg.norm(hbv)))
                ha_proj.append(hap)
                hb_proj.append(hbp)
                any_hv = True
                continue
        z = np.zeros(0, dtype=np.float64)
        ha_raw.append(z)
        hb_raw.append(z)
        ha_proj.append(z)
        hb_proj.append(z)

    return {
        "xo_raw": xo_raw,
        "yo_raw": yo_raw,
        "xo_proj": xo_proj,
        "yo_proj": yo_proj,
        "ha_raw": ha_raw,
        "hb_raw": hb_raw,
        "ha_proj": ha_proj,
        "hb_proj": hb_proj,
        "any_hv": any_hv,
        "cos_mean": float(np.mean(cos_vals)) if cos_vals else float("nan"),
        "de50_mean": float(np.mean(de50_vals)) if de50_vals else float("nan"),
    }


def compute_orth_suite_r16(
    delta_obs: np.ndarray,
    delta_pred: np.ndarray,
    delta_transfer: np.ndarray,
    target_gene_idx: Sequence[int],
    source_grand_mean: Optional[np.ndarray] = None,
    delta_half_a: Optional[np.ndarray] = None,
    delta_half_b: Optional[np.ndarray] = None,
    half_valid: Optional[np.ndarray] = None,
    mode: str = "all_nontarget",
    top_k: int = 50,
    n_boot: int = 4000,
    seed: int = DEFAULT_R16_SEED,
    gene_keep_mask: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """
    Compute `r_orth_cal` (when `source_grand_mean is None`) or `r_dorth_cal` (when provided)
    together with joint perturbation bootstrap `97.5%` / `95.0%` CIs, `valid_frac`,
    and explicit `rho_ab_orth` (`PREREG_I v4 §1.1, §3.1, §4.1`).
    """
    if mode not in ("all_nontarget", "resp_top50"):
        raise ValueError(f"Unsupported mode: {mode!r}")

    kept = _extract_pert_gene_indices(
        delta_obs=delta_obs,
        t_mean=delta_transfer,
        pred=delta_pred,
        target_cols=target_gene_idx,
        mode=mode,
        top_k=top_k,
        source_grand_mean=source_grand_mean,
        gene_keep_mask=gene_keep_mask,
    )
    m = len(kept)
    endpoint = "r_orth_cal" if source_grand_mean is None else "r_dorth_cal"
    if m == 0:
        out = summarize_boot_r16(float("nan"), np.zeros(0), 0, n_boot)
        out.update(
            {
                "endpoint": endpoint,
                "mode": mode,
                "n_eval_perts": 0,
                "rho_ab_orth": float("nan"),
                "noise_ceiling_r_orth": float("nan"),
                "rho_ab_raw": float("nan"),
                "noise_ceiling_r": float("nan"),
                "r_cal_est": float("nan"),
                "r_proj_raw_est": float("nan"),
                "slope_ols_est": float("nan"),
                "cos_agreement": float("nan"),
                "de50_r_mean": float("nan"),
            }
        )
        return out

    streams = _project_streams(
        kept=kept,
        delta_obs=delta_obs,
        pred=delta_pred,
        t_mean=delta_transfer,
        source_grand_mean=source_grand_mean,
        half_a_obs=delta_half_a,
        half_b_obs=delta_half_b,
        half_valid=half_valid,
    )

    is_dorth = source_grand_mean is not None
    m_raw = _moments(streams["xo_raw"], streams["yo_raw"])
    m_proj = _moments(streams["xo_proj"], streams["yo_proj"])

    pt_r_raw, pt_slope_raw, pt_sd_raw = _moments_to_corr_and_slope(*m_raw)
    pt_r_proj, _, _ = _moments_to_corr_and_slope(
        *m_proj, zero_degenerate_y_as_zero=is_dorth
    )

    rng = np.random.default_rng(int(seed))
    draws = rng.integers(0, m, size=(int(n_boot), m))
    b_r_raw, b_slope_raw, _ = _moments_to_corr_and_slope(*m_raw, draws=draws)
    b_r_proj, _, _ = _moments_to_corr_and_slope(
        *m_proj, draws=draws, zero_degenerate_y_as_zero=is_dorth
    )

    if streams["any_hv"]:
        mh_raw = _moments(streams["ha_raw"], streams["hb_raw"])
        mh_proj = _moments(streams["ha_proj"], streams["hb_proj"])

        pt_rho_raw = float(_moments_to_corr_and_slope(*mh_raw)[0])
        pt_ceil_raw = spearman_brown_ceiling(pt_rho_raw)
        b_rho_raw = _moments_to_corr_and_slope(*mh_raw, draws=draws)[0]
        b_ceil_raw = np.asarray([spearman_brown_ceiling(float(v)) for v in b_rho_raw], dtype=np.float64)

        pt_rho_orth = float(_moments_to_corr_and_slope(*mh_proj)[0])
        pt_ceil_orth = spearman_brown_ceiling(pt_rho_orth)
        b_rho_orth = _moments_to_corr_and_slope(*mh_proj, draws=draws)[0]
        b_ceil_orth = np.asarray([spearman_brown_ceiling(float(v)) for v in b_rho_orth], dtype=np.float64)
    else:
        pt_rho_raw = pt_ceil_raw = float("nan")
        pt_rho_orth = pt_ceil_orth = float("nan")
        b_ceil_raw = np.full(int(n_boot), float("nan"), dtype=np.float64)
        b_ceil_orth = np.full(int(n_boot), float("nan"), dtype=np.float64)

    def _safe_disatten(r_val: float, c_val: float) -> float:
        if is_dorth and np.isfinite(r_val) and abs(r_val) <= NULL_NORM_TOL and np.isfinite(c_val) and c_val > 1e-8:
            return 0.0
        return disattenuate(r_val, c_val)

    pt_cal_raw = disattenuate(float(pt_r_raw), pt_ceil_raw)
    pt_cal_proj = _safe_disatten(float(pt_r_proj), pt_ceil_orth)
    b_cal_proj = np.asarray(
        [_safe_disatten(float(r), float(c)) for r, c in zip(b_r_proj, b_ceil_orth)],
        dtype=np.float64,
    )

    boot_sum = summarize_boot_r16(pt_cal_proj, b_cal_proj, m, int(n_boot))
    boot_sum.update(
        {
            "endpoint": endpoint,
            "mode": mode,
            "n_eval_perts": int(m),
            "rho_ab_orth": float(pt_rho_orth),
            "noise_ceiling_r_orth": float(pt_ceil_orth),
            "rho_ab_raw": float(pt_rho_raw),
            "noise_ceiling_r": float(pt_ceil_raw),
            "r_cal_est": float(pt_cal_raw),
            "r_raw_est": float(pt_r_raw),
            "r_proj_raw_est": float(pt_r_proj),
            "slope_ols_est": float(pt_slope_raw),
            "sd_pred_over_sd_obs": float(pt_sd_raw),
            "cos_agreement": float(streams["cos_mean"]),
            "de50_r_mean": float(streams["de50_mean"]),
        }
    )
    return boot_sum


def compute_paired_contrast_r16(
    delta_obs: np.ndarray,
    pred_a: np.ndarray,
    pred_b: np.ndarray,
    delta_transfer: np.ndarray,
    target_gene_idx: Sequence[int],
    source_grand_mean: Optional[np.ndarray] = None,
    delta_half_a: Optional[np.ndarray] = None,
    delta_half_b: Optional[np.ndarray] = None,
    half_valid: Optional[np.ndarray] = None,
    mode: str = "all_nontarget",
    top_k: int = 50,
    n_boot: int = 4000,
    seed: int = DEFAULT_R16_SEED,
    gene_keep_mask: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """
    Joint perturbation bootstrap paired contrast `r(pred_a) - r(pred_b)` under identical
    perturbation draws and identical projected split-half noise ceiling (`PREREG_I v4 §3.1–3.2`).
    """
    if mode not in ("all_nontarget", "resp_top50"):
        raise ValueError(f"Unsupported mode: {mode!r}")

    kept = _extract_pert_gene_indices(
        delta_obs=delta_obs,
        t_mean=delta_transfer,
        pred=pred_a,
        target_cols=target_gene_idx,
        mode=mode,
        top_k=top_k,
        extra_pred=pred_b,
        source_grand_mean=source_grand_mean,
        gene_keep_mask=gene_keep_mask,
    )
    m = len(kept)
    if m == 0:
        out = summarize_boot_r16(float("nan"), np.zeros(0), 0, n_boot)
        out.update({"n_eval_perts": 0, "rho_ab_orth": float("nan")})
        return out

    sa = _project_streams(
        kept=kept,
        delta_obs=delta_obs,
        pred=pred_a,
        t_mean=delta_transfer,
        source_grand_mean=source_grand_mean,
        half_a_obs=delta_half_a,
        half_b_obs=delta_half_b,
        half_valid=half_valid,
    )
    sb = _project_streams(
        kept=kept,
        delta_obs=delta_obs,
        pred=pred_b,
        t_mean=delta_transfer,
        source_grand_mean=source_grand_mean,
        half_a_obs=None,
        half_b_obs=None,
        half_valid=None,
    )

    is_dorth = source_grand_mean is not None
    ma = _moments(sa["xo_proj"], sa["yo_proj"])
    mb = _moments(sb["xo_proj"], sb["yo_proj"])

    pt_ra = float(_moments_to_corr_and_slope(*ma, zero_degenerate_y_as_zero=is_dorth)[0])
    pt_rb = float(_moments_to_corr_and_slope(*mb, zero_degenerate_y_as_zero=is_dorth)[0])

    rng = np.random.default_rng(int(seed))
    draws = rng.integers(0, m, size=(int(n_boot), m))
    b_ra = _moments_to_corr_and_slope(*ma, draws=draws, zero_degenerate_y_as_zero=is_dorth)[0]
    b_rb = _moments_to_corr_and_slope(*mb, draws=draws, zero_degenerate_y_as_zero=is_dorth)[0]

    if sa["any_hv"]:
        mh = _moments(sa["ha_proj"], sa["hb_proj"])
        pt_rho_orth = float(_moments_to_corr_and_slope(*mh)[0])
        pt_ceil = spearman_brown_ceiling(pt_rho_orth)
        b_rho = _moments_to_corr_and_slope(*mh, draws=draws)[0]
        b_ceil = np.asarray([spearman_brown_ceiling(float(v)) for v in b_rho], dtype=np.float64)
    else:
        pt_rho_orth = float("nan")
        pt_ceil = float("nan")
        b_ceil = np.full(int(n_boot), float("nan"), dtype=np.float64)

    def _d(r_val: float, c_val: float) -> float:
        if is_dorth and np.isfinite(r_val) and abs(r_val) <= NULL_NORM_TOL and np.isfinite(c_val) and c_val > 1e-8:
            return 0.0
        return disattenuate(r_val, c_val)

    pt_diff = _d(pt_ra, pt_ceil) - _d(pt_rb, pt_ceil)
    b_diff = np.asarray(
        [_d(float(ra), float(c)) - _d(float(rb), float(c)) for ra, rb, c in zip(b_ra, b_rb, b_ceil)],
        dtype=np.float64,
    )
    out = summarize_boot_r16(pt_diff, b_diff, m, int(n_boot))
    out.update(
        {
            "n_eval_perts": int(m),
            "rho_ab_orth": float(pt_rho_orth),
            "a_est": _d(pt_ra, pt_ceil),
            "b_est": _d(pt_rb, pt_ceil),
        }
    )
    return out


def compute_paired_model_vs_all_mismatched_r16(
    delta_obs: np.ndarray,
    pred_matched: np.ndarray,
    mismatched_preds: Dict[str, np.ndarray],
    delta_transfer: np.ndarray,
    target_gene_idx: Sequence[int],
    source_grand_mean: Optional[np.ndarray] = None,
    delta_half_a: Optional[np.ndarray] = None,
    delta_half_b: Optional[np.ndarray] = None,
    half_valid: Optional[np.ndarray] = None,
    mode: str = "all_nontarget",
    top_k: int = 50,
    n_boot: int = 4000,
    seed: int = DEFAULT_R16_SEED,
) -> Dict[str, Any]:
    """Compute paired `(matched - mismatched_s')` bootstrap summaries across every source `s'`."""
    per_source: Dict[str, Dict[str, Any]] = {}
    for sk in sorted(mismatched_preds.keys()):
        per_source[sk] = compute_paired_contrast_r16(
            delta_obs=delta_obs,
            pred_a=pred_matched,
            pred_b=mismatched_preds[sk],
            delta_transfer=delta_transfer,
            target_gene_idx=target_gene_idx,
            source_grand_mean=source_grand_mean,
            delta_half_a=delta_half_a,
            delta_half_b=delta_half_b,
            half_valid=half_valid,
            mode=mode,
            top_k=top_k,
            n_boot=n_boot,
            seed=seed,
        )
    return {"per_source": per_source, "n_sources": len(per_source)}


def compute_winner_curse_weight_shift_r16(
    ctrl_mean_libnorm: np.ndarray,
    delta_obs: np.ndarray,
    target_gene_idx: Sequence[int],
    top_k: int = 50,
    n_boot: int = 4000,
    seed: int = DEFAULT_R16_SEED,
) -> Dict[str, Any]:
    """
    Secondary Mechanism Hypothesis `M2` diagnostic (`PREREG_I v4 §4.2`):
    mean saturation-weight `w(mu) = mu / (1 + mu)` shift comparing top-k `|delta_obs|` genes
    against `all_nontarget`.
    """
    mu = np.maximum(np.asarray(ctrl_mean_libnorm, dtype=np.float64), 0.0)
    w = mu / (1.0 + mu)
    n_t, n_g = delta_obs.shape
    shifts: List[float] = []
    for p in range(n_t):
        valid = np.isfinite(delta_obs[p]) & np.isfinite(w)
        tc = int(target_gene_idx[p])
        if 0 <= tc < n_g:
            valid[tc] = False
        idx = np.flatnonzero(valid)
        if idx.size < 2:
            continue
        k = min(int(top_k), idx.size)
        top_idx = idx[np.argsort(np.abs(delta_obs[p, idx]), kind="mergesort")[-k:]]
        shifts.append(float(np.mean(w[top_idx]) - np.mean(w[idx])))
    if not shifts:
        return summarize_boot_r16(float("nan"), np.zeros(0), 0, n_boot)
    arr = np.asarray(shifts, dtype=np.float64)
    rng = np.random.default_rng(int(seed))
    draws = rng.integers(0, arr.size, size=(int(n_boot), arr.size))
    return summarize_boot_r16(float(np.mean(arr)), arr[draws].mean(axis=1), arr.size, int(n_boot))
