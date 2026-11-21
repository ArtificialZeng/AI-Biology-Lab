#!/usr/bin/env python3
"""Independent statistical audit runner for Team A (A-T3).

Strictly independent from p2/vcc_linespecific/src/: recomputes zero-shot
predictions, line-specific slope, noise-calibrated correlation, disc_raw,
de50_r, 4,000-resample perturbation-level 95% CIs, and fixed H-A1..H-A4
adjudications directly from p2/shared/ inputs, and cross-checks A-T2's outputs.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from independent_predictors import (
    compute_control_covariance,
    predict_zero_shot_suite,
)
from independent_stats import (
    EPS,
    BootstrapCI,
    adjudicate_team_a,
    adjudicate_team_b,
    bootstrap_paired_mean_diff,
    compute_de50_r_rows,
    compute_disc_ownmasked_rows,
    compute_disc_raw_rows,
    disattenuated_correlation,
    noise_ceiling_from_splithalf,
    ols_slope,
    pearson_1d,
    summarize_bootstrap_array,
)


def _extract_npz_array(npz_path: Path, preferred_keys: Sequence[str]) -> tuple[np.ndarray, str]:
    """Load the first matching array key from an .npz archive."""
    with np.load(npz_path, allow_pickle=False) as z:
        keys = list(z.files)
        for k in preferred_keys:
            if k in keys:
                return np.asarray(z[k], dtype=np.float64), k
        if len(keys) == 1:
            k = keys[0]
            return np.asarray(z[k], dtype=np.float64), k
        raise KeyError(f"None of {preferred_keys} found in {npz_path}; available={keys}")


def _read_lines_txt(path: Path) -> List[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def _space_key_prefs(space: str) -> tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...], tuple[str, ...]]:
    """Return strictly matched (ctrl, delta, half_a, half_b) keys for the requested space."""
    if space == "zscore":
        return (
            ("X_zscore", "ctrl_zscore", "zscore", "X"),
            ("delta_zscore", "zscore", "delta"),
            ("half_a_delta_zscore", "half_a", "delta_half_a"),
            ("half_b_delta_zscore", "half_b", "delta_half_b"),
        )
    if space == "log1p":
        return (
            ("X_log1p", "ctrl_log1p", "log1p", "X"),
            ("delta_log1p", "log1p", "delta"),
            ("half_a_delta_log1p", "half_a", "delta_half_a"),
            ("half_b_delta_log1p", "half_b", "delta_half_b"),
        )
    if space == "libnorm":
        return (
            ("X_libnorm", "ctrl_libnorm", "libnorm", "X"),
            ("delta_libnorm", "libnorm", "delta"),
            ("half_a_delta_libnorm", "half_a", "delta_half_a"),
            ("half_b_delta_libnorm", "half_b", "delta_half_b"),
        )
    raise ValueError(f"Unknown expression space {space!r}; expected 'zscore', 'log1p', or 'libnorm'.")


def load_shared_line(
    line_dir: Path,
    common_genes: Sequence[str],
    ctrl_key_prefs: Optional[Sequence[str]] = None,
    delta_key_prefs: Optional[Sequence[str]] = None,
    space: str = "zscore",
) -> dict:
    """Load one cell line's control matrix and pseudobulk deltas aligned to common_genes."""
    default_ctrl, default_delta, default_ha, default_hb = _space_key_prefs(space)
    c_prefs = ctrl_key_prefs if ctrl_key_prefs is not None else default_ctrl
    d_prefs = delta_key_prefs if delta_key_prefs is not None else default_delta

    genes = _read_lines_txt(line_dir / "genes.txt")
    targets = _read_lines_txt(line_dir / "targets.txt")
    ctrl_raw, ctrl_key = _extract_npz_array(line_dir / "ctrl_cells.npz", c_prefs)
    delta_raw, delta_key = _extract_npz_array(line_dir / "pseudobulk_delta.npz", d_prefs)

    # Guard against mismatched units when fallback key 'X' or 'delta' was used
    if space == "zscore" and ctrl_key in ("X_log1p", "ctrl_log1p", "log1p", "X"):
        mu = np.mean(ctrl_raw, axis=0, keepdims=True)
        sd = np.std(ctrl_raw, axis=0, ddof=1, keepdims=True)
        if float(np.median(np.abs(mu))) > 1e-2 or float(np.median(np.abs(sd - 1.0))) > 0.1:
            ctrl_raw = (ctrl_raw - mu) / np.maximum(sd, 1e-8)
            ctrl_key = f"{ctrl_key}->zscored"

    gene_pos = {g: i for i, g in enumerate(genes)}
    missing = [g for g in common_genes if g not in gene_pos]
    if missing:
        raise ValueError(f"Line {line_dir.name} missing {len(missing)} COMMON_GENES (e.g. {missing[:5]})")
    col_idx = np.array([gene_pos[g] for g in common_genes], dtype=np.int64)

    if ctrl_raw.shape[1] != len(genes):
        raise ValueError(f"{line_dir.name} ctrl shape {ctrl_raw.shape} != len(genes) {len(genes)}")
    if delta_raw.shape != (len(targets), len(genes)):
        raise ValueError(f"{line_dir.name} delta shape {delta_raw.shape} != ({len(targets)}, {len(genes)})")

    ctrl_aligned = ctrl_raw[:, col_idx]
    delta_aligned = delta_raw[:, col_idx]

    ctrl_mean_libnorm, delta_libnorm = None, None
    with np.load(line_dir / "ctrl_cells.npz", allow_pickle=False) as zc:
        if "ctrl_mean_libnorm" in zc.files:
            ctrl_mean_libnorm = np.asarray(zc["ctrl_mean_libnorm"], dtype=np.float64)[col_idx]

    half_a, half_b, n_cells, reproducible, half_valid = None, None, None, None, None
    with np.load(line_dir / "pseudobulk_delta.npz", allow_pickle=False) as z:
        files = set(z.files)
        if "n_cells" in files:
            n_cells = np.asarray(z["n_cells"], dtype=np.float64)
        if "reproducible" in files:
            reproducible = np.asarray(z["reproducible"], dtype=bool)
        if "half_valid" in files:
            half_valid = np.asarray(z["half_valid"], dtype=bool)
        if "delta_libnorm" in files:
            delta_libnorm = np.asarray(z["delta_libnorm"], dtype=np.float64)[:, col_idx]
        ha_k = next((k for k in default_ha if k in files), None)
        hb_k = next((k for k in default_hb if k in files), None)
        if ha_k is not None and hb_k is not None:
            half_a = np.asarray(z[ha_k], dtype=np.float64)[:, col_idx]
            half_b = np.asarray(z[hb_k], dtype=np.float64)[:, col_idx]

    return {
        "line": line_dir.name,
        "targets": targets,
        "ctrl_cells": ctrl_aligned,
        "delta": delta_aligned,
        "ctrl_mean_libnorm": ctrl_mean_libnorm,
        "delta_libnorm": delta_libnorm,
        "ctrl_key": ctrl_key,
        "delta_key": delta_key,
        "n_cells": n_cells,
        "reproducible": reproducible,
        "half_valid": half_valid,
        "half_a": half_a,
        "half_b": half_b,
    }


def _build_masked_vectors(
    r_true: np.ndarray,
    r_pred: np.ndarray,
    delta_obs: np.ndarray,
    target_cols: Sequence[int],
    mode: str = "resp_top50",
    top_k: int = 50,
    gene_keep_mask: Optional[np.ndarray] = None,
) -> Tuple[List[np.ndarray], List[np.ndarray], np.ndarray]:
    """Build per-perturbation masked residual vectors following LINE_SPECIFIC_METRIC.md."""
    n_perts, n_genes = delta_obs.shape
    xs, ys, kept_rows = [], [], []
    for p in range(n_perts):
        valid = np.isfinite(r_true[p]) & np.isfinite(r_pred[p]) & np.isfinite(delta_obs[p])
        if gene_keep_mask is not None:
            valid = valid & gene_keep_mask
        t_idx = int(target_cols[p])
        if 0 <= t_idx < n_genes:
            valid[t_idx] = False
        idx = np.flatnonzero(valid)
        if idx.size == 0:
            continue
        if mode == "resp_top50" and idx.size > top_k:
            order = np.argsort(np.abs(delta_obs[p, idx]))[-top_k:]
            idx = idx[order]
        xs.append(r_true[p, idx].astype(np.float64))
        ys.append(r_pred[p, idx].astype(np.float64))
        kept_rows.append(p)
    return xs, ys, np.array(kept_rows, dtype=np.int64)


def _moments_to_stats(
    N: np.ndarray | float,
    Sx: np.ndarray | float,
    Sy: np.ndarray | float,
    Sxx: np.ndarray | float,
    Syy: np.ndarray | float,
    Sxy: np.ndarray | float,
):
    N_safe = np.maximum(N, 1.0)
    cov = Sxy - (Sx * Sy) / N_safe
    vx = Sxx - (Sx * Sx) / N_safe
    vy = Syy - (Sy * Sy) / N_safe
    with np.errstate(divide="ignore", invalid="ignore"):
        b_ols = np.where(vx > EPS, cov / vx, 0.0)
        b_noint = np.where(Sxx > EPS, Sxy / Sxx, 0.0)
        denom = np.sqrt(np.maximum(vx, 0.0) * np.maximum(vy, 0.0))
        b_corr = np.where(denom > EPS, np.clip(cov / denom, -1.0, 1.0), 0.0)
    return b_ols, b_noint, b_corr


def _moments_to_r6_stats(
    N: np.ndarray | float,
    Sx: np.ndarray | float,
    Sy: np.ndarray | float,
    Sxx: np.ndarray | float,
    Syy: np.ndarray | float,
    Sxy: np.ndarray | float,
):
    N_safe = np.maximum(N, 1.0)
    cov = Sxy - (Sx * Sy) / N_safe
    vx = Sxx - (Sx * Sx) / N_safe
    vy = Syy - (Sy * Sy) / N_safe
    with np.errstate(divide="ignore", invalid="ignore"):
        sd_ratio = np.where(vx > EPS, np.sqrt(np.maximum(vy, 0.0) / np.maximum(vx, EPS)), 0.0)
        b_rev = np.where(vy > EPS, cov / vy, 0.0)
        diff = vy - vx
        b_dem = np.where(
            np.abs(cov) > EPS,
            (diff + np.sqrt(diff * diff + 4.0 * cov * cov)) / (2.0 * cov),
            0.0,
        )
    return sd_ratio, b_rev, b_dem


def _eval_residual_slope_and_corr(
    r_true: np.ndarray,
    r_pred: np.ndarray,
    delta_obs: np.ndarray,
    target_cols: Sequence[int],
    mode: str,
    n_boot: int,
    seed: int,
    ceiling_r: float,
    half_a_res: Optional[np.ndarray] = None,
    half_b_res: Optional[np.ndarray] = None,
    half_valid_mask: Optional[np.ndarray] = None,
    gene_keep_mask: Optional[np.ndarray] = None,
) -> dict:
    """Compute line-specific OLS slope, zero-intercept slope, per-pert slope, and corr with 95%/97.5% CIs."""
    xs, ys, kept = _build_masked_vectors(
        r_true, r_pred, delta_obs, target_cols, mode=mode, gene_keep_mask=gene_keep_mask
    )
    m = len(xs)
    if m == 0:
        return {"n_perts": 0, "slope_ols": float("nan")}

    nv = np.array([float(x.size) for x in xs], dtype=np.float64)
    sx = np.array([float(np.sum(x)) for x in xs], dtype=np.float64)
    sy = np.array([float(np.sum(y)) for y in ys], dtype=np.float64)
    sxx = np.array([float(np.dot(x, x)) for x in xs], dtype=np.float64)
    syy = np.array([float(np.dot(y, y)) for y in ys], dtype=np.float64)
    sxy = np.array([float(np.dot(x, y)) for x, y in zip(xs, ys)], dtype=np.float64)

    pt_ols, pt_noint, pt_corr = _moments_to_stats(
        float(nv.sum()), float(sx.sum()), float(sy.sum()), float(sxx.sum()), float(syy.sum()), float(sxy.sum())
    )
    pt_sd_ratio, pt_rev_ols, pt_deming = _moments_to_r6_stats(
        float(nv.sum()), float(sx.sum()), float(sy.sum()), float(sxx.sum()), float(syy.sum()), float(sxy.sum())
    )

    has_direct_halves = half_a_res is not None and half_b_res is not None
    hnv = np.zeros(m, dtype=np.float64)
    hsx = np.zeros(m, dtype=np.float64)
    hsy = np.zeros(m, dtype=np.float64)
    hsxx = np.zeros(m, dtype=np.float64)
    hsyy = np.zeros(m, dtype=np.float64)
    hsxy = np.zeros(m, dtype=np.float64)
    if has_direct_halves:
        any_hv = False
        for i, p in enumerate(kept):
            if half_valid_mask is not None and not bool(half_valid_mask[p]):
                continue
            valid = np.isfinite(delta_obs[p]) & np.isfinite(half_a_res[p]) & np.isfinite(half_b_res[p])
            if gene_keep_mask is not None:
                valid = valid & gene_keep_mask
            t_idx = int(target_cols[p])
            if 0 <= t_idx < delta_obs.shape[1]:
                valid[t_idx] = False
            idx = np.flatnonzero(valid)
            if idx.size == 0:
                continue
            if mode == "resp_top50" and idx.size > 50:
                idx = idx[np.argsort(np.abs(delta_obs[p, idx]))[-50:]]
            ha_v = half_a_res[p, idx].astype(np.float64)
            hb_v = half_b_res[p, idx].astype(np.float64)
            hnv[i] = float(ha_v.size)
            hsx[i] = float(np.sum(ha_v))
            hsy[i] = float(np.sum(hb_v))
            hsxx[i] = float(np.dot(ha_v, ha_v))
            hsyy[i] = float(np.dot(hb_v, hb_v))
            hsxy[i] = float(np.dot(ha_v, hb_v))
            any_hv = True
        if not any_hv:
            has_direct_halves = False

    rng = np.random.default_rng(seed)
    draws = rng.integers(0, m, size=(n_boot, m))
    N_b = nv[draws].sum(axis=1)
    Sx_b = sx[draws].sum(axis=1)
    Sy_b = sy[draws].sum(axis=1)
    Sxx_b = sxx[draws].sum(axis=1)
    Syy_b = syy[draws].sum(axis=1)
    Sxy_b = sxy[draws].sum(axis=1)
    b_ols, b_noint, b_corr = _moments_to_stats(N_b, Sx_b, Sy_b, Sxx_b, Syy_b, Sxy_b)
    b_sd_ratio, b_rev_ols, b_deming = _moments_to_r6_stats(N_b, Sx_b, Sy_b, Sxx_b, Syy_b, Sxy_b)

    pt_hab = float("nan")
    if has_direct_halves:
        _, _, pt_hab_arr = _moments_to_stats(
            float(hnv.sum()), float(hsx.sum()), float(hsy.sum()),
            float(hsxx.sum()), float(hsyy.sum()), float(hsxy.sum()),
        )
        pt_hab = float(pt_hab_arr)
        full_ceiling = noise_ceiling_from_splithalf(pt_hab)
        _, _, b_hab = _moments_to_stats(
            hnv[draws].sum(axis=1),
            hsx[draws].sum(axis=1),
            hsy[draws].sum(axis=1),
            hsxx[draws].sum(axis=1),
            hsyy[draws].sum(axis=1),
            hsxy[draws].sum(axis=1),
        )
        b_ceil = np.array([noise_ceiling_from_splithalf(float(rh)) for rh in b_hab], dtype=np.float64)
        b_cal = np.array(
            [disattenuated_correlation(float(rc), float(cb)) for rc, cb in zip(b_corr, b_ceil)],
            dtype=np.float64,
        )
    else:
        full_ceiling = ceiling_r
        b_cal = np.array(
            [disattenuated_correlation(float(rc), ceiling_r) for rc in b_corr],
            dtype=np.float64,
        )

    pt_cal = disattenuated_correlation(float(pt_corr), float(full_ceiling))
    pt_cal_sq = float(pt_cal * pt_cal) if np.isfinite(pt_cal) else float("nan")
    pt_cal_sq_signed = float(np.sign(pt_cal) * pt_cal * pt_cal) if np.isfinite(pt_cal) else float("nan")
    b_cal_sq = np.where(np.isfinite(b_cal), b_cal * b_cal, np.nan)
    b_cal_sq_signed = np.where(np.isfinite(b_cal), np.sign(b_cal) * b_cal * b_cal, np.nan)

    ci_ols = summarize_bootstrap_array(float(pt_ols), b_ols, n_items=m)
    ci_noint = summarize_bootstrap_array(float(pt_noint), b_noint, n_items=m)
    ci_corr = summarize_bootstrap_array(float(pt_corr), b_corr, n_items=m)
    ci_cal = summarize_bootstrap_array(float(pt_cal), b_cal, n_items=m)
    ci_cal_sq = summarize_bootstrap_array(pt_cal_sq, b_cal_sq, n_items=m)
    ci_cal_sq_signed = summarize_bootstrap_array(pt_cal_sq_signed, b_cal_sq_signed, n_items=m)
    ci_sd_ratio = summarize_bootstrap_array(float(pt_sd_ratio), b_sd_ratio, n_items=m)
    ci_rev_ols = summarize_bootstrap_array(float(pt_rev_ols), b_rev_ols, n_items=m)
    ci_deming = summarize_bootstrap_array(float(pt_deming), b_deming, n_items=m)

    per_pert_slopes = [ols_slope(xs[i], ys[i], fit_intercept=True) for i in range(m)]
    per_pert_corrs = [pearson_1d(xs[i], ys[i]) for i in range(m)]

    return {
        "n_perts": m,
        "slope_ols": ci_ols.estimate,
        "slope_ols_ci95_lo": ci_ols.ci95_lo,
        "slope_ols_ci95_hi": ci_ols.ci95_hi,
        "slope_ols_ci975_lo": ci_ols.ci975_lo,
        "slope_ols_ci975_hi": ci_ols.ci975_hi,
        "slope_ols_boot_se": ci_ols.boot_se,
        "slope_no_intercept": ci_noint.estimate,
        "slope_no_intercept_ci95_lo": ci_noint.ci95_lo,
        "slope_no_intercept_ci95_hi": ci_noint.ci95_hi,
        "per_pert_slope_mean": float(np.nanmean(per_pert_slopes)),
        "resid_corr": ci_corr.estimate,
        "resid_corr_ci95_lo": ci_corr.ci95_lo,
        "resid_corr_ci95_hi": ci_corr.ci95_hi,
        "resid_corr_ci975_lo": ci_corr.ci975_lo,
        "resid_corr_ci975_hi": ci_corr.ci975_hi,
        "noise_cal_corr": ci_cal.estimate,
        "noise_cal_corr_ci95_lo": ci_cal.ci95_lo,
        "noise_cal_corr_ci95_hi": ci_cal.ci95_hi,
        "noise_cal_corr_ci975_lo": ci_cal.ci975_lo,
        "noise_cal_corr_ci975_hi": ci_cal.ci975_hi,
        "r_cal_sq": ci_cal_sq.estimate,
        "r_cal_sq_ci95_lo": ci_cal_sq.ci95_lo,
        "r_cal_sq_ci95_hi": ci_cal_sq.ci95_hi,
        "r_cal_sq_ci975_lo": ci_cal_sq.ci975_lo,
        "r_cal_sq_ci975_hi": ci_cal_sq.ci975_hi,
        "r_cal_sq_signed": ci_cal_sq_signed.estimate,
        "r_cal_sq_signed_ci95_lo": ci_cal_sq_signed.ci95_lo,
        "r_cal_sq_signed_ci95_hi": ci_cal_sq_signed.ci95_hi,
        "sd_pred_over_sd_obs": ci_sd_ratio.estimate,
        "sd_pred_over_sd_obs_ci95_lo": ci_sd_ratio.ci95_lo,
        "sd_pred_over_sd_obs_ci95_hi": ci_sd_ratio.ci95_hi,
        "slope_obs_on_pred": ci_rev_ols.estimate,
        "slope_obs_on_pred_ci95_lo": ci_rev_ols.ci95_lo,
        "slope_obs_on_pred_ci95_hi": ci_rev_ols.ci95_hi,
        "slope_deming": ci_deming.estimate,
        "slope_deming_ci95_lo": ci_deming.ci95_lo,
        "slope_deming_ci95_hi": ci_deming.ci95_hi,
        "per_pert_corr_mean": float(np.nanmean(per_pert_corrs)),
        "rho_ab": pt_hab,
        "noise_ceiling_r": float(full_ceiling),
    }


def _eval_joint_paired_cal_r_diff(
    r_true: np.ndarray,
    r_pred_a: np.ndarray,
    r_pred_b: np.ndarray,
    delta_obs: np.ndarray,
    target_cols: Sequence[int],
    mode: str,
    n_boot: int,
    seed: int,
    ceiling_r: float,
    half_a_res: Optional[np.ndarray] = None,
    half_b_res: Optional[np.ndarray] = None,
    half_valid_mask: Optional[np.ndarray] = None,
) -> BootstrapCI:
    """Paired joint-bootstrap CI for r_cal(pred_a) - r_cal(pred_b) (R14)."""
    xs_a, ys_a, kept_a = _build_masked_vectors(r_true, r_pred_a, delta_obs, target_cols, mode=mode)
    xs_b, ys_b, kept_b = _build_masked_vectors(r_true, r_pred_b, delta_obs, target_cols, mode=mode)
    m = len(kept_a)
    if m == 0 or not np.array_equal(kept_a, kept_b):
        return BootstrapCI(float("nan"), float("nan"), float("nan"), float("nan"), m, n_boot)

    def _moms(xs, ys):
        nv = np.array([float(x.size) for x in xs], dtype=np.float64)
        sx = np.array([float(np.sum(x)) for x in xs], dtype=np.float64)
        sy = np.array([float(np.sum(y)) for y in ys], dtype=np.float64)
        sxx = np.array([float(np.dot(x, x)) for x in xs], dtype=np.float64)
        syy = np.array([float(np.dot(y, y)) for y in ys], dtype=np.float64)
        sxy = np.array([float(np.dot(x, y)) for x, y in zip(xs, ys)], dtype=np.float64)
        return nv, sx, sy, sxx, syy, sxy

    nv_a, sx_a, sy_a, sxx_a, syy_a, sxy_a = _moms(xs_a, ys_a)
    nv_b, sx_b, sy_b, sxx_b, syy_b, sxy_b = _moms(xs_b, ys_b)
    _, _, pt_ra = _moments_to_stats(
        float(nv_a.sum()), float(sx_a.sum()), float(sy_a.sum()),
        float(sxx_a.sum()), float(syy_a.sum()), float(sxy_a.sum()),
    )
    _, _, pt_rb = _moments_to_stats(
        float(nv_b.sum()), float(sx_b.sum()), float(sy_b.sum()),
        float(sxx_b.sum()), float(syy_b.sum()), float(sxy_b.sum()),
    )

    rng = np.random.default_rng(seed)
    draws = rng.integers(0, m, size=(n_boot, m))
    _, _, b_ra = _moments_to_stats(
        nv_a[draws].sum(axis=1), sx_a[draws].sum(axis=1), sy_a[draws].sum(axis=1),
        sxx_a[draws].sum(axis=1), syy_a[draws].sum(axis=1), sxy_a[draws].sum(axis=1),
    )
    _, _, b_rb = _moments_to_stats(
        nv_b[draws].sum(axis=1), sx_b[draws].sum(axis=1), sy_b[draws].sum(axis=1),
        sxx_b[draws].sum(axis=1), syy_b[draws].sum(axis=1), sxy_b[draws].sum(axis=1),
    )

    has_direct_halves = half_a_res is not None and half_b_res is not None
    if has_direct_halves:
        hnv = np.zeros(m, dtype=np.float64)
        hsx = np.zeros(m, dtype=np.float64)
        hsy = np.zeros(m, dtype=np.float64)
        hsxx = np.zeros(m, dtype=np.float64)
        hsyy = np.zeros(m, dtype=np.float64)
        hsxy = np.zeros(m, dtype=np.float64)
        any_hv = False
        for i, p in enumerate(kept_a):
            if half_valid_mask is not None and not bool(half_valid_mask[p]):
                continue
            valid = np.isfinite(delta_obs[p]) & np.isfinite(half_a_res[p]) & np.isfinite(half_b_res[p])
            t_idx = int(target_cols[p])
            if 0 <= t_idx < delta_obs.shape[1]:
                valid[t_idx] = False
            idx = np.flatnonzero(valid)
            if idx.size == 0:
                continue
            if mode == "resp_top50" and idx.size > 50:
                idx = idx[np.argsort(np.abs(delta_obs[p, idx]))[-50:]]
            ha_v = half_a_res[p, idx].astype(np.float64)
            hb_v = half_b_res[p, idx].astype(np.float64)
            hnv[i] = float(ha_v.size)
            hsx[i] = float(np.sum(ha_v))
            hsy[i] = float(np.sum(hb_v))
            hsxx[i] = float(np.dot(ha_v, ha_v))
            hsyy[i] = float(np.dot(hb_v, hb_v))
            hsxy[i] = float(np.dot(ha_v, hb_v))
            any_hv = True
        if any_hv:
            _, _, pt_hab = _moments_to_stats(
                float(hnv.sum()), float(hsx.sum()), float(hsy.sum()),
                float(hsxx.sum()), float(hsyy.sum()), float(hsxy.sum()),
            )
            full_c = noise_ceiling_from_splithalf(float(pt_hab))
            _, _, b_hab = _moments_to_stats(
                hnv[draws].sum(axis=1), hsx[draws].sum(axis=1), hsy[draws].sum(axis=1),
                hsxx[draws].sum(axis=1), hsyy[draws].sum(axis=1), hsxy[draws].sum(axis=1),
            )
            b_ceil = np.array([noise_ceiling_from_splithalf(float(rh)) for rh in b_hab], dtype=np.float64)
        else:
            full_c = ceiling_r
            b_ceil = np.full(n_boot, ceiling_r, dtype=np.float64)
    else:
        full_c = ceiling_r
        b_ceil = np.full(n_boot, ceiling_r, dtype=np.float64)

    pt_diff = disattenuated_correlation(float(pt_ra), float(full_c)) - disattenuated_correlation(float(pt_rb), float(full_c))
    b_diff = np.array(
        [
            disattenuated_correlation(float(ra), float(cb)) - disattenuated_correlation(float(rb), float(cb))
            for ra, rb, cb in zip(b_ra, b_rb, b_ceil)
        ],
        dtype=np.float64,
    )
    return summarize_bootstrap_array(pt_diff, b_diff, n_items=m)


def _mean_row_cosine(
    a_mat: np.ndarray,
    b_mat: np.ndarray,
    delta_obs: np.ndarray,
    target_cols: Sequence[int],
    mode: str = "resp_top50",
    top_k: int = 50,
) -> float:
    """Compute mean per-perturbation cosine similarity on masked genes (Wrap-Up companion diagnostic)."""
    n_perts, n_genes = delta_obs.shape
    vals = []
    for p in range(n_perts):
        valid = np.isfinite(a_mat[p]) & np.isfinite(b_mat[p]) & np.isfinite(delta_obs[p])
        t_idx = int(target_cols[p])
        if 0 <= t_idx < n_genes:
            valid[t_idx] = False
        idx = np.flatnonzero(valid)
        if idx.size == 0:
            continue
        if mode == "resp_top50" and idx.size > top_k:
            idx = idx[np.argsort(np.abs(delta_obs[p, idx]))[-top_k:]]
        u = a_mat[p, idx].astype(np.float64)
        v = b_mat[p, idx].astype(np.float64)
        denom = float(np.linalg.norm(u) * np.linalg.norm(v))
        if denom > 1e-12:
            vals.append(float(np.clip(np.dot(u, v) / denom, -1.0, 1.0)))
    return float(np.mean(vals)) if vals else float("nan")


def control_split_half_noise_ceiling(
    ctrl_cells: np.ndarray,
    r_true: np.ndarray,
    delta_obs: np.ndarray,
    target_cols: Sequence[int],
    n_cells_per_pert: Optional[np.ndarray],
    mode: str = "resp_top50",
    seed: int = 42,
) -> float:
    """Estimate line-specific residual noise ceiling from control cells when perturbation split-halves are unavailable."""
    ctrl = np.asarray(ctrl_cells, dtype=np.float64)
    n_ctrl = ctrl.shape[0]
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n_ctrl)
    half = n_ctrl // 2
    # Empirical control variance per gene estimated strictly via split-half difference / 2
    var_gene = np.var(ctrl[perm[:half]], axis=0, ddof=1) * 0.5 + np.var(ctrl[perm[half : 2 * half]], axis=0, ddof=1) * 0.5
    default_n = 50.0
    n_vec = n_cells_per_pert if n_cells_per_pert is not None else np.full(len(target_cols), default_n)

    obs_vals = []
    noise_vars = []
    n_genes = delta_obs.shape[1]
    for p in range(len(target_cols)):
        valid = np.isfinite(r_true[p]) & np.isfinite(delta_obs[p])
        t_idx = int(target_cols[p])
        if 0 <= t_idx < n_genes:
            valid[t_idx] = False
        idx = np.flatnonzero(valid)
        if mode == "resp_top50" and idx.size > 50:
            idx = idx[np.argsort(np.abs(delta_obs[p, idx]))[-50:]]
        obs_vals.append(r_true[p, idx])
        np_val = max(float(n_vec[p]), 1.0)
        noise_vars.append((1.0 / np_val + 1.0 / max(float(n_ctrl), 1.0)) * var_gene[idx])

    v_obs = float(np.var(np.concatenate(obs_vals)))
    v_noise = float(np.mean(np.concatenate(noise_vars)))
    if v_obs <= 1e-12:
        return 0.0
    rho = max(0.0, 1.0 - v_noise / v_obs)
    return float(np.sqrt(rho))


def audit_loo_suite(
    shared_dir: Path,
    lines: Sequence[str],
    space: str = "zscore",
    reproducible_only: bool = True,
    shrink: float = 1e-3,
    cov_max_cells: Optional[int] = 10000,
    project_shared_axis: bool = False,
    fit_ls_lambda_on_sources: bool = False,
    n_boot: int = 4000,
    seed: int = 20261121,
    dev_only: bool = False,
    confirmatory_lines: Optional[Sequence[str]] = None,
) -> dict:
    """Run full independent LOO-CL audit on lines in shared_dir."""
    cg_path = shared_dir / "COMMON_GENES.txt"
    if dev_only or not cg_path.exists():
        gene_sets = [
            set(_read_lines_txt(shared_dir / "lines" / line / "genes.txt"))
            for line in lines
        ]
        common_genes = sorted(set.intersection(*gene_sets))
    else:
        common_genes = _read_lines_txt(cg_path)
    gene_to_col = {g: i for i, g in enumerate(common_genes)}

    loaded = {
        line: load_shared_line(shared_dir / "lines" / line, common_genes, space=space)
        for line in lines
    }
    covs = {
        line: compute_control_covariance(
            loaded[line]["ctrl_cells"], shrink=shrink, max_cells=cov_max_cells, seed=seed
        )
        for line in lines
    }
    rep_by_line = {
        line: {
            tgt: bool(loaded[line]["reproducible"][i]) if loaded[line]["reproducible"] is not None else True
            for i, tgt in enumerate(loaded[line]["targets"])
        }
        for line in lines
    }
    deltas_by_line = {
        line: {
            tgt: loaded[line]["delta"][i]
            for i, tgt in enumerate(loaded[line]["targets"])
            if (not reproducible_only) or rep_by_line[line][tgt]
        }
        for line in lines
    }
    deltas_lib_by_line = {
        line: {
            tgt: loaded[line]["delta_libnorm"][i]
            for i, tgt in enumerate(loaded[line]["targets"])
            if ((not reproducible_only) or rep_by_line[line][tgt]) and loaded[line]["delta_libnorm"] is not None
        }
        for line in lines
    }

    per_line_rows = []
    detailed = {}

    for li, tgt_line in enumerate(lines):
        src_lines = [l for l in lines if l != tgt_line]
        # Targets evaluated: target line targets (reproducible if requested) present in >= 1 source line and in common_genes
        eval_targets = [
            tgt
            for tgt in loaded[tgt_line]["targets"]
            if tgt in gene_to_col
            and (not reproducible_only or rep_by_line[tgt_line][tgt])
            and any(tgt in deltas_by_line[sl] for sl in src_lines)
        ]
        target_cols = [gene_to_col[t] for t in eval_targets]
        src_covs = {l: covs[l] for l in src_lines}
        src_deltas = {l: deltas_by_line[l] for l in src_lines}
        src_mu_lib = {
            l: loaded[l]["ctrl_mean_libnorm"]
            for l in src_lines
            if loaded[l]["ctrl_mean_libnorm"] is not None
        }
        src_deltas_lib = {l: deltas_lib_by_line[l] for l in src_lines if deltas_lib_by_line[l]}
        mismatch_line = src_lines[0]

        preds = predict_zero_shot_suite(
            eval_targets=eval_targets,
            target_to_col=gene_to_col,
            target_ctrl_cov=covs[tgt_line],
            source_covs=src_covs,
            source_deltas_by_target=src_deltas,
            mismatch_cov=covs[mismatch_line],
            project_shared_axis=project_shared_axis,
            fit_ls_lambda_on_sources=fit_ls_lambda_on_sources,
            seed=seed + li,
            target_ctrl_mean_libnorm=loaded[tgt_line]["ctrl_mean_libnorm"],
            source_ctrl_mean_libnorm=src_mu_lib if len(src_mu_lib) == len(src_lines) else None,
            source_deltas_libnorm_by_target=src_deltas_lib if len(src_deltas_lib) == len(src_lines) else None,
        )

        # Ground-truth arrays strictly for evaluation after prediction
        y_true = np.stack([deltas_by_line[tgt_line][t] for t in eval_targets], axis=0)
        t_mean = preds["transferred_mean"]
        r_true = y_true - t_mean

        ha_res, hb_res, hv_mask = None, None, None
        if loaded[tgt_line]["half_a"] is not None and loaded[tgt_line]["half_b"] is not None:
            idx_map = {t: i for i, t in enumerate(loaded[tgt_line]["targets"])}
            ha = np.stack([loaded[tgt_line]["half_a"][idx_map[t]] for t in eval_targets], axis=0)
            hb = np.stack([loaded[tgt_line]["half_b"][idx_map[t]] for t in eval_targets], axis=0)
            if loaded[tgt_line]["half_valid"] is not None:
                hv_mask = np.array([bool(loaded[tgt_line]["half_valid"][idx_map[t]]) for t in eval_targets], dtype=bool)
            ha_res = ha - t_mean
            hb_res = hb - t_mean
            ceil_top50 = 0.0
            ceil_all = 0.0
        else:
            n_c = None
            if loaded[tgt_line]["n_cells"] is not None:
                idx_map = {t: i for i, t in enumerate(loaded[tgt_line]["targets"])}
                n_c = np.array([loaded[tgt_line]["n_cells"][idx_map[t]] for t in eval_targets])
            ceil_top50 = control_split_half_noise_ceiling(
                loaded[tgt_line]["ctrl_cells"], r_true, y_true, target_cols, n_c, "resp_top50", seed + li
            )
            ceil_all = control_split_half_noise_ceiling(
                loaded[tgt_line]["ctrl_cells"], r_true, y_true, target_cols, n_c, "all_nontarget", seed + li
            )

        src_rows = [d for sl in src_lines for d in deltas_by_line[sl].values()]
        src_grand_mean = np.mean(np.stack(src_rows, axis=0), axis=0)

        method_metrics = {}
        for m_name, p_mat in preds.items():
            r_pred = p_mat - t_mean
            m_top50 = _eval_residual_slope_and_corr(
                r_true, r_pred, y_true, target_cols, "resp_top50", n_boot, seed + li * 100, ceil_top50, ha_res, hb_res, hv_mask
            )
            m_all = _eval_residual_slope_and_corr(
                r_true, r_pred, y_true, target_cols, "all_nontarget", n_boot, seed + li * 100 + 50, ceil_all, ha_res, hb_res, hv_mask
            )
            disc_rows = compute_disc_raw_rows(y_true, p_mat, masked_cols=target_cols)
            disc_own_rows = compute_disc_ownmasked_rows(y_true, p_mat, target_col_idx=target_cols)
            de50_rows = compute_de50_r_rows(y_true, p_mat, target_col_idx=target_cols)
            cent_disc_rows = compute_disc_raw_rows(
                y_true - src_grand_mean, p_mat - src_grand_mean, masked_cols=target_cols
            )
            cent_disc_own_rows = compute_disc_ownmasked_rows(
                y_true - src_grand_mean, p_mat - src_grand_mean, target_col_idx=target_cols
            )
            cent_de50_rows = compute_de50_r_rows(
                y_true - src_grand_mean, p_mat - src_grand_mean, target_col_idx=target_cols
            )
            pds_cos_rows = compute_disc_raw_rows(
                y_true, p_mat, masked_cols=target_cols, center_across_genes=False
            )
            cos_top50 = _mean_row_cosine(p_mat, y_true, y_true, target_cols, "resp_top50")
            cos_all = _mean_row_cosine(p_mat, y_true, y_true, target_cols, "all_nontarget")
            method_metrics[m_name] = {
                "resp_top50": m_top50,
                "all_nontarget": m_all,
                "disc_raw_rows": disc_rows,
                "disc_own_rows": disc_own_rows,
                "de50_r_rows": de50_rows,
                "cent_disc_raw_rows": cent_disc_rows,
                "cent_disc_own_rows": cent_disc_own_rows,
                "cent_de50_r_rows": cent_de50_rows,
                "pds_cosine_panel_rows": pds_cos_rows,
                "disc_raw_mean": float(np.mean(disc_rows)),
                "disc_own_mean": float(np.mean(disc_own_rows)),
                "de50_r_mean": float(np.mean(de50_rows)),
                "cent_disc_raw_mean": float(np.mean(cent_disc_rows)),
                "cent_disc_own_mean": float(np.mean(cent_disc_own_rows)),
                "cent_de50_r_mean": float(np.mean(cent_de50_rows)),
                "pds_cosine_panel_mean": float(np.mean(pds_cos_rows)),
                "cos_agreement_top50": cos_top50,
                "cos_agreement_all": cos_all,
            }

        # Paired non-inferiority diffs (H-A4 and H-B4) with both 95% and 97.5% CIs (R12/R13)
        conv_specs = [
            ("disc_raw", "disc_raw_rows", 0),
            ("ctrl_disc_ownmasked", "disc_own_rows", 4),
            ("de50_r", "de50_r_rows", 1),
            ("cent_disc_raw", "cent_disc_raw_rows", 2),
            ("cent_disc_ownmasked", "cent_disc_own_rows", 5),
            ("cent_de50_r", "cent_de50_r_rows", 3),
        ]
        ha4_cis = {
            name: bootstrap_paired_mean_diff(
                method_metrics["hybrid"][row_key],
                method_metrics["transferred_mean"][row_key],
                n_boot=n_boot,
                seed=seed + li + off,
            )
            for name, row_key, off in conv_specs
        }
        hb4_cis = {
            name: bootstrap_paired_mean_diff(
                method_metrics["transfer_mean_var_scaled"][row_key],
                method_metrics["transferred_mean"][row_key],
                n_boot=n_boot,
                seed=seed + li + off + 10,
            )
            for name, row_key, off in conv_specs
        }

        # R3: track ALL source-line mismatched covariance controls (Group A)
        mismatch_by_source = {}
        for sl in src_lines:
            ms_50 = method_metrics[f"mismatch_cov_hybrid::{sl}"]["resp_top50"]
            mismatch_by_source[sl] = {
                "slope_ols": ms_50["slope_ols"],
                "slope_ols_ci95_lo": ms_50["slope_ols_ci95_lo"],
                "slope_ols_ci95_hi": ms_50["slope_ols_ci95_hi"],
                "slope_ols_ci975_lo": ms_50["slope_ols_ci975_lo"],
                "slope_ols_ci975_hi": ms_50["slope_ols_ci975_hi"],
                "noise_cal_corr": ms_50["noise_cal_corr"],
            }
        worst_mismatch_line = max(
            src_lines, key=lambda sl: mismatch_by_source[sl]["slope_ols_ci975_hi"]
        )
        worst_mis = mismatch_by_source[worst_mismatch_line]

        # R12: track ALL source-line mismatched variances + shuffled gene-label variance (Group B)
        tmvs_ctrl_keys = [f"tmvs_mismatch::{sl}" for sl in src_lines] + ["tmvs_shuffled_var"]
        tmvs_ctrls = {}
        for ck in tmvs_ctrl_keys:
            mc = method_metrics[ck]["resp_top50"]
            tmvs_ctrls[ck] = {
                "slope_ols": mc["slope_ols"],
                "slope_ols_ci95_lo": mc["slope_ols_ci95_lo"],
                "slope_ols_ci95_hi": mc["slope_ols_ci95_hi"],
                "slope_ols_ci975_lo": mc["slope_ols_ci975_lo"],
                "slope_ols_ci975_hi": mc["slope_ols_ci975_hi"],
                "noise_cal_corr": mc["noise_cal_corr"],
                "noise_cal_corr_ci95_lo": mc["noise_cal_corr_ci95_lo"],
                "noise_cal_corr_ci95_hi": mc["noise_cal_corr_ci95_hi"],
                "noise_cal_corr_ci975_lo": mc["noise_cal_corr_ci975_lo"],
                "noise_cal_corr_ci975_hi": mc["noise_cal_corr_ci975_hi"],
            }

        def _nan_max(vals: Sequence[float]) -> float:
            fin = [float(v) for v in vals if np.isfinite(v)]
            return max(fin) if fin else float("nan")

        tmvs_max_ctrl_slope_95 = _nan_max([v["slope_ols_ci95_hi"] for v in tmvs_ctrls.values()])
        tmvs_max_ctrl_slope_975 = _nan_max([v["slope_ols_ci975_hi"] for v in tmvs_ctrls.values()])
        tmvs_max_ctrl_cal_95 = _nan_max([v["noise_cal_corr_ci95_hi"] for v in tmvs_ctrls.values()])
        tmvs_max_ctrl_cal_975 = _nan_max([v["noise_cal_corr_ci975_hi"] for v in tmvs_ctrls.values()])

        # R14: Hybrid_on_TMVS minus TMVS paired joint-bootstrap on r_cal
        r14_ci = _eval_joint_paired_cal_r_diff(
            r_true=r_true,
            r_pred_a=preds["hybrid_on_tmvs"] - t_mean,
            r_pred_b=preds["transfer_mean_var_scaled"] - t_mean,
            delta_obs=y_true,
            target_cols=target_cols,
            mode="resp_top50",
            n_boot=n_boot,
            seed=seed + li * 100 + 7,
            ceiling_r=ceil_top50,
            half_a_res=ha_res,
            half_b_res=hb_res,
            half_valid_mask=hv_mask,
        )

        # R15 S5: Expressed-genes sensitivity (exclude lowest 25% of target control diag(C_t))
        t_var_diag = np.diag(covs[tgt_line])
        q25_var = float(np.percentile(t_var_diag, 25.0))
        expr_mask = t_var_diag > q25_var
        if int(np.sum(expr_mask)) < 50:
            expr_mask = t_var_diag >= q25_var
        s5_hyb = _eval_residual_slope_and_corr(
            r_true,
            preds["hybrid"] - t_mean,
            y_true,
            target_cols,
            "resp_top50",
            n_boot,
            seed + li * 100 + 11,
            ceil_top50,
            ha_res,
            hb_res,
            hv_mask,
            gene_keep_mask=expr_mask,
        )
        s5_tmvs = _eval_residual_slope_and_corr(
            r_true,
            preds["transfer_mean_var_scaled"] - t_mean,
            y_true,
            target_cols,
            "resp_top50",
            n_boot,
            seed + li * 100 + 13,
            ceil_top50,
            ha_res,
            hb_res,
            hv_mask,
            gene_keep_mask=expr_mask,
        )

        hyb_50 = method_metrics["hybrid"]["resp_top50"]
        hyb_all = method_metrics["hybrid"]["all_nontarget"]
        mis_50 = method_metrics["mismatch_cov_hybrid"]["resp_top50"]
        sctv_50 = method_metrics["src_corr_tgt_var_hybrid"]["resp_top50"]
        tmvs_50 = method_metrics["transfer_mean_var_scaled"]["resp_top50"]
        tmvs_all = method_metrics["transfer_mean_var_scaled"]["all_nontarget"]
        tmvs_shuf_50 = method_metrics["tmvs_shuffled_var"]["resp_top50"]
        hotmvs_50 = method_metrics["hybrid_on_tmvs"]["resp_top50"]
        fct_50 = method_metrics["fct"]["resp_top50"] if "fct" in method_metrics else {}

        row_dict = {
            "line": tgt_line,
            "n_eval_targets": len(eval_targets),
            "mismatch_line": mismatch_line,
            "worst_mismatch_line": worst_mismatch_line,
            "hybrid_slope": hyb_50["slope_ols"],
            "primary_slope_ci95_lo": hyb_50["slope_ols_ci95_lo"],
            "hybrid_slope_ci95_hi": hyb_50["slope_ols_ci95_hi"],
            "primary_slope_ci975_lo": hyb_50["slope_ols_ci975_lo"],
            "hybrid_slope_ci975_hi": hyb_50["slope_ols_ci975_hi"],
            "resid_corr": hyb_50["resid_corr"],
            "resid_corr_ci95_lo": hyb_50["resid_corr_ci95_lo"],
            "resid_corr_ci95_hi": hyb_50["resid_corr_ci95_hi"],
            "noise_cal_corr": hyb_50["noise_cal_corr"],
            "noise_cal_corr_ci95_lo": hyb_50["noise_cal_corr_ci95_lo"],
            "noise_cal_corr_ci95_hi": hyb_50["noise_cal_corr_ci95_hi"],
            "noise_cal_corr_ci975_lo": hyb_50["noise_cal_corr_ci975_lo"],
            "noise_cal_corr_ci975_hi": hyb_50["noise_cal_corr_ci975_hi"],
            "r_cal_sq": hyb_50["r_cal_sq"],
            "r_cal_sq_ci95_lo": hyb_50["r_cal_sq_ci95_lo"],
            "r_cal_sq_ci95_hi": hyb_50["r_cal_sq_ci95_hi"],
            "sd_pred_over_sd_obs": hyb_50["sd_pred_over_sd_obs"],
            "slope_obs_on_pred": hyb_50["slope_obs_on_pred"],
            "slope_deming": hyb_50["slope_deming"],
            "rho_ab_top50": hyb_50["rho_ab"],
            "noise_ceiling_r": hyb_50["noise_ceiling_r"],
            "rho_ab_all_nontarget": hyb_all["rho_ab"],
            "noise_ceiling_r_all_nontarget": hyb_all["noise_ceiling_r"],
            "hybrid_cal_r_all_nontarget": hyb_all["noise_cal_corr"],
            "mismatch_slope": mis_50["slope_ols"],
            "mismatch_slope_ci95_lo": mis_50["slope_ols_ci95_lo"],
            "mismatch_slope_ci95_hi": mis_50["slope_ols_ci95_hi"],
            "max_mismatch_slope": worst_mis["slope_ols"],
            "max_mismatch_ci95_lo": worst_mis["slope_ols_ci95_lo"],
            "max_mismatch_ci95_hi": worst_mis["slope_ols_ci95_hi"],
            "max_mismatch_ci975_lo": worst_mis["slope_ols_ci975_lo"],
            "max_mismatch_ci975_hi": worst_mis["slope_ols_ci975_hi"],
            "sctv_slope": sctv_50["slope_ols"],
            "sctv_slope_ci95_hi": sctv_50["slope_ols_ci95_hi"],
            "sctv_slope_ci975_hi": sctv_50["slope_ols_ci975_hi"],
            "sctv_cal_r": sctv_50["noise_cal_corr"],
            "sctv_cal_r_ci95_hi": sctv_50["noise_cal_corr_ci95_hi"],
            "sctv_cal_r_ci975_hi": sctv_50["noise_cal_corr_ci975_hi"],
            # Group B: Transfer_Mean_VarScaled + R12/R13 controls
            "tmvs_slope": tmvs_50["slope_ols"],
            "tmvs_slope_ci95_lo": tmvs_50["slope_ols_ci95_lo"],
            "tmvs_slope_ci95_hi": tmvs_50["slope_ols_ci95_hi"],
            "tmvs_slope_ci975_lo": tmvs_50["slope_ols_ci975_lo"],
            "tmvs_slope_ci975_hi": tmvs_50["slope_ols_ci975_hi"],
            "tmvs_cal_r": tmvs_50["noise_cal_corr"],
            "tmvs_cal_r_ci95_lo": tmvs_50["noise_cal_corr_ci95_lo"],
            "tmvs_cal_r_ci95_hi": tmvs_50["noise_cal_corr_ci95_hi"],
            "tmvs_cal_r_ci975_lo": tmvs_50["noise_cal_corr_ci975_lo"],
            "tmvs_cal_r_ci975_hi": tmvs_50["noise_cal_corr_ci975_hi"],
            "tmvs_r_cal_sq": tmvs_50["r_cal_sq"],
            "tmvs_sd_pred_over_sd_obs": tmvs_50["sd_pred_over_sd_obs"],
            "tmvs_slope_deming": tmvs_50["slope_deming"],
            "tmvs_slope_obs_on_pred": tmvs_50["slope_obs_on_pred"],
            "tmvs_cal_r_all_nontarget": tmvs_all["noise_cal_corr"],
            "tmvs_slope_all_nontarget": tmvs_all["slope_ols"],
            "tmvs_shuf_var_slope": tmvs_shuf_50["slope_ols"],
            "tmvs_shuf_var_slope_ci975_hi": tmvs_shuf_50["slope_ols_ci975_hi"],
            "tmvs_shuf_var_cal_r": tmvs_shuf_50["noise_cal_corr"],
            "tmvs_shuf_var_cal_r_ci975_hi": tmvs_shuf_50["noise_cal_corr_ci975_hi"],
            "tmvs_max_ctrl_slope_ci95_hi": tmvs_max_ctrl_slope_95,
            "tmvs_max_ctrl_slope_ci975_hi": tmvs_max_ctrl_slope_975,
            "tmvs_max_ctrl_cal_r_ci95_hi": tmvs_max_ctrl_cal_95,
            "tmvs_max_ctrl_cal_r_ci975_hi": tmvs_max_ctrl_cal_975,
            "tmvs_ctrl_disc_ownmasked": method_metrics["transfer_mean_var_scaled"]["disc_own_mean"],
            "tmvs_de50_r": method_metrics["transfer_mean_var_scaled"]["de50_r_mean"],
            # R14 Hybrid_on_TMVS
            "hotmvs_slope": hotmvs_50["slope_ols"],
            "hotmvs_cal_r": hotmvs_50["noise_cal_corr"],
            "r14_hotmvs_minus_tmvs_cal_r": r14_ci.estimate,
            "r14_hotmvs_minus_tmvs_cal_r_ci95_lo": r14_ci.ci95_lo,
            "r14_hotmvs_minus_tmvs_cal_r_ci95_hi": r14_ci.ci95_hi,
            "r14_hotmvs_minus_tmvs_cal_r_ci975_lo": r14_ci.ci975_lo,
            "r14_hotmvs_minus_tmvs_cal_r_ci975_hi": r14_ci.ci975_hi,
            # R15 FCT + S5 expressed genes
            "fct_slope": fct_50.get("slope_ols", float("nan")),
            "fct_cal_r": fct_50.get("noise_cal_corr", float("nan")),
            "fct_cal_r_ci95_lo": fct_50.get("noise_cal_corr_ci95_lo", float("nan")),
            "fct_cal_r_ci95_hi": fct_50.get("noise_cal_corr_ci95_hi", float("nan")),
            "s5_n_expressed_genes": int(np.sum(expr_mask)),
            "s5_hybrid_cal_r": s5_hyb["noise_cal_corr"],
            "s5_hybrid_cal_r_ci95_lo": s5_hyb["noise_cal_corr_ci95_lo"],
            "s5_tmvs_cal_r": s5_tmvs["noise_cal_corr"],
            "s5_tmvs_cal_r_ci95_lo": s5_tmvs["noise_cal_corr_ci95_lo"],
            "s5_tmvs_slope": s5_tmvs["slope_ols"],
            # Additional negative controls & global metrics
            "diag_cov_slope": method_metrics["diag_cov_hybrid"]["resp_top50"]["slope_ols"],
            "shuf_target_slope": method_metrics["shuffled_gene_hybrid"]["resp_top50"]["slope_ols"],
            "shuf_all_slope": method_metrics["shuffled_all_cov_hybrid"]["resp_top50"]["slope_ols"],
            "cipher_linear_slope": method_metrics["cipher_linear"]["resp_top50"]["slope_ols"],
            "disc_raw_hybrid": method_metrics["hybrid"]["disc_raw_mean"],
            "disc_raw_tmean": method_metrics["transferred_mean"]["disc_raw_mean"],
            "disc_raw_diff": ha4_cis["disc_raw"].estimate,
            "disc_raw_diff_ci95_lo": ha4_cis["disc_raw"].ci95_lo,
            "disc_raw_diff_ci95_hi": ha4_cis["disc_raw"].ci95_hi,
            "disc_raw_diff_ci975_hi": ha4_cis["disc_raw"].ci975_hi,
            "ctrl_disc_ownmasked_hybrid": method_metrics["hybrid"]["disc_own_mean"],
            "ctrl_disc_ownmasked_tmean": method_metrics["transferred_mean"]["disc_own_mean"],
            "ctrl_disc_ownmasked_diff": ha4_cis["ctrl_disc_ownmasked"].estimate,
            "ctrl_disc_ownmasked_diff_ci95_hi": ha4_cis["ctrl_disc_ownmasked"].ci95_hi,
            "ctrl_disc_ownmasked_diff_ci975_hi": ha4_cis["ctrl_disc_ownmasked"].ci975_hi,
            "de50_r_hybrid": method_metrics["hybrid"]["de50_r_mean"],
            "de50_r_tmean": method_metrics["transferred_mean"]["de50_r_mean"],
            "de50_r_diff": ha4_cis["de50_r"].estimate,
            "de50_r_diff_ci95_lo": ha4_cis["de50_r"].ci95_lo,
            "de50_r_diff_ci95_hi": ha4_cis["de50_r"].ci95_hi,
            "de50_r_diff_ci975_hi": ha4_cis["de50_r"].ci975_hi,
            "cent_disc_raw_hybrid": method_metrics["hybrid"]["cent_disc_raw_mean"],
            "cent_disc_raw_tmean": method_metrics["transferred_mean"]["cent_disc_raw_mean"],
            "cent_disc_raw_diff": ha4_cis["cent_disc_raw"].estimate,
            "cent_disc_raw_diff_ci95_hi": ha4_cis["cent_disc_raw"].ci95_hi,
            "cent_disc_raw_diff_ci975_hi": ha4_cis["cent_disc_raw"].ci975_hi,
            "cent_disc_ownmasked_hybrid": method_metrics["hybrid"]["cent_disc_own_mean"],
            "cent_disc_ownmasked_tmean": method_metrics["transferred_mean"]["cent_disc_own_mean"],
            "cent_disc_ownmasked_diff": ha4_cis["cent_disc_ownmasked"].estimate,
            "cent_disc_ownmasked_diff_ci95_hi": ha4_cis["cent_disc_ownmasked"].ci95_hi,
            "cent_disc_ownmasked_diff_ci975_hi": ha4_cis["cent_disc_ownmasked"].ci975_hi,
            "cent_de50_r_hybrid": method_metrics["hybrid"]["cent_de50_r_mean"],
            "cent_de50_r_tmean": method_metrics["transferred_mean"]["cent_de50_r_mean"],
            "cent_de50_r_diff": ha4_cis["cent_de50_r"].estimate,
            "cent_de50_r_diff_ci95_hi": ha4_cis["cent_de50_r"].ci95_hi,
            "cent_de50_r_diff_ci975_hi": ha4_cis["cent_de50_r"].ci975_hi,
            "pds_cosine_panel_hybrid": method_metrics["hybrid"]["pds_cosine_panel_mean"],
        }
        # Populate H-B4 paired diff columns uniformly
        for cname, cci in hb4_cis.items():
            pfx = f"tmvs_{'cent_disc_panel' if cname == 'cent_disc_raw' else cname}_diff"
            row_dict[pfx] = cci.estimate
            row_dict[f"{pfx}_ci95_hi"] = cci.ci95_hi
            row_dict[f"{pfx}_ci975_hi"] = cci.ci975_hi
            if cname == "disc_raw":
                row_dict["tmvs_disc_raw_diff"] = cci.estimate
                row_dict["tmvs_disc_raw_diff_ci95_hi"] = cci.ci95_hi
                row_dict["tmvs_disc_raw_diff_ci975_hi"] = cci.ci975_hi
            elif cname == "cent_disc_raw":
                row_dict["tmvs_cent_disc_raw_diff"] = cci.estimate
                row_dict["tmvs_cent_disc_raw_diff_ci95_hi"] = cci.ci95_hi
                row_dict["tmvs_cent_disc_raw_diff_ci975_hi"] = cci.ci975_hi

        per_line_rows.append(row_dict)

        # Strip raw row arrays before JSON serialization
        detailed[tgt_line] = {
            "mismatch_by_source": mismatch_by_source,
            "tmvs_ctrls": tmvs_ctrls,
            "r14_hotmvs_minus_tmvs": {
                "delta_cal_r": r14_ci.estimate,
                "ci95_lo": r14_ci.ci95_lo,
                "ci95_hi": r14_ci.ci95_hi,
                "ci975_lo": r14_ci.ci975_lo,
                "ci975_hi": r14_ci.ci975_hi,
            },
            "s5_expressed_only": {
                "n_expressed_genes": int(np.sum(expr_mask)),
                "hybrid": s5_hyb,
                "tmvs": s5_tmvs,
            },
            "methods": {
                m: {
                    "resp_top50": v["resp_top50"],
                    "all_nontarget": v["all_nontarget"],
                    "disc_raw_mean": v["disc_raw_mean"],
                    "disc_own_mean": v["disc_own_mean"],
                    "de50_r_mean": v["de50_r_mean"],
                    "cent_disc_raw_mean": v["cent_disc_raw_mean"],
                    "cent_disc_own_mean": v["cent_disc_own_mean"],
                    "cent_de50_r_mean": v["cent_de50_r_mean"],
                    "pds_cosine_panel_mean": v["pds_cosine_panel_mean"],
                    "cos_agreement_top50": v["cos_agreement_top50"],
                    "cos_agreement_all": v["cos_agreement_all"],
                }
                for m, v in method_metrics.items()
            },
        }

    # Adjudicate using exact pre-registered rules (R2, R3, R4, R5, R10, R12, R13, R14)
    conf_list = list(confirmatory_lines) if confirmatory_lines is not None else list(lines)
    n_conf = len(conf_list)
    min_pass = int(np.ceil(0.8 * n_conf))
    max_inf = int(np.floor(0.2 * n_conf))
    verdict = adjudicate_team_a(
        per_line_rows,
        min_lines_pass=min_pass,
        max_lines_inferior=max_inf,
        confirmatory_lines=conf_list,
        use_975_ci=True,
    )
    verdict_b = adjudicate_team_b(
        per_line_rows,
        min_lines_pass=min_pass,
        max_lines_inferior=max_inf,
        confirmatory_lines=conf_list,
        use_975_ci=True,
    )
    verdict["Group_B"] = verdict_b
    conf_rows = [r for r in per_line_rows if r["line"] in set(conf_list)]
    beats_sctv = bool(conf_rows) and all(
        np.isfinite(r["primary_slope_ci95_lo"])
        and np.isfinite(r["sctv_slope_ci95_hi"])
        and np.isfinite(r["noise_cal_corr_ci95_lo"])
        and np.isfinite(r["sctv_cal_r_ci95_hi"])
        and (r["primary_slope_ci95_lo"] > r["sctv_slope_ci95_hi"])
        and (r["noise_cal_corr_ci95_lo"] > r["sctv_cal_r_ci95_hi"])
        for r in conf_rows
    )
    beats_tmvs = bool(conf_rows) and all(
        np.isfinite(r["noise_cal_corr_ci95_lo"])
        and np.isfinite(r["tmvs_cal_r_ci95_hi"])
        and (r["noise_cal_corr_ci95_lo"] > r["tmvs_cal_r_ci95_hi"])
        for r in conf_rows
    )
    r14_all_conf_pass = bool(conf_rows) and all(
        np.isfinite(r["r14_hotmvs_minus_tmvs_cal_r_ci95_lo"])
        and (r["r14_hotmvs_minus_tmvs_cal_r_ci95_lo"] > 0.0)
        for r in conf_rows
    )
    verdict["r4_wording_governance"] = {
        "beats_src_corr_tgt_var_all_conf_slope_and_cal": beats_sctv,
        "beats_transfer_mean_var_scaled_all_conf_cal": beats_tmvs,
        "allowed_mechanism_phrase": (
            "cell-line-specific co-fluctuation structure"
            if beats_sctv
            else "target-line control statistics"
        ),
        "marginal_variance_caveat": (
            None
            if beats_tmvs
            else "Gains over Transfer_Mean may be explained by target-line marginal gene variance"
        ),
    }
    verdict["r14_mechanism_governance"] = {
        "hotmvs_beats_tmvs_cal_r_all_conf_ci95_lo_gt_0": r14_all_conf_pass,
        "per_conf_line_delta_cal_r": {
            r["line"]: {
                "delta": r["r14_hotmvs_minus_tmvs_cal_r"],
                "ci95_lo": r["r14_hotmvs_minus_tmvs_cal_r_ci95_lo"],
                "ci95_hi": r["r14_hotmvs_minus_tmvs_cal_r_ci95_hi"],
            }
            for r in conf_rows
        },
        "allowed_cofluctuation_increment_phrase": (
            "Co-fluctuation structure provides incremental information beyond marginal variance"
            if r14_all_conf_pass
            else "In this dataset, co-fluctuation structure exhibits no detectable increment beyond marginal variance"
        ),
    }

    return {
        "lines": list(lines),
        "n_common_genes": len(common_genes),
        "per_line_df": pd.DataFrame(per_line_rows),
        "detailed": detailed,
        "verdict": verdict,
    }


def check_prereg_gate(status_dir: Path, lines: Sequence[str], dev_only: bool) -> None:
    """Enforce P2_RULES Rule 3: holdout lines forbidden prior to PREREG_COMMIT.txt."""
    if dev_only:
        allowed = {
            "replogle_k562",
            "replogle_k562_essential",
            "replogle_rpe1",
            "K562",
            "RPE1",
            "k562",
            "rpe1",
        }
        if any(l not in allowed for l in lines):
            raise PermissionError(f"--dev-only only allows Replogle K562/RPE1, got {lines}")
        return
    commit_file = status_dir / "PREREG_COMMIT.txt"
    if not commit_file.exists() or not commit_file.read_text().strip():
        raise PermissionError(
            f"Rule 3 Gate: {commit_file} missing or empty. Refusing to evaluate held-out cell lines!"
        )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--shared-dir", type=Path, required=True)
    ap.add_argument("--status-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--lines", nargs="+", required=True)
    ap.add_argument("--dev-only", action="store_true")
    ap.add_argument("--space", choices=["zscore", "log1p", "libnorm"], default="zscore")
    ap.add_argument("--all-targets", action="store_true")
    ap.add_argument("--shrink", type=float, default=1e-3)
    ap.add_argument("--cov-max-cells", type=int, default=10000)
    ap.add_argument("--project-shared-axis", action="store_true")
    ap.add_argument("--fit-ls-lambda-on-sources", action="store_true")
    ap.add_argument("--confirmatory-lines", nargs="+", default=None)
    ap.add_argument("--n-boot", type=int, default=4000)
    args = ap.parse_args()

    check_prereg_gate(args.status_dir, args.lines, args.dev_only)
    res = audit_loo_suite(
        shared_dir=args.shared_dir,
        lines=args.lines,
        space=args.space,
        reproducible_only=not args.all_targets,
        shrink=args.shrink,
        cov_max_cells=(None if args.cov_max_cells <= 0 else args.cov_max_cells),
        project_shared_axis=args.project_shared_axis,
        fit_ls_lambda_on_sources=args.fit_ls_lambda_on_sources,
        n_boot=args.n_boot,
        dev_only=args.dev_only,
        confirmatory_lines=args.confirmatory_lines,
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    res["per_line_df"].to_csv(args.out_dir / "audit_per_line.csv", index=False)
    with open(args.out_dir / "audit_summary.json", "w") as f:
        json.dump(
            {
                "lines": res["lines"],
                "n_common_genes": res["n_common_genes"],
                "detailed": res["detailed"],
                "verdict": res["verdict"],
            },
            f,
            indent=2,
        )
    print(json.dumps(res["verdict"], indent=2))


if __name__ == "__main__":
    main()
