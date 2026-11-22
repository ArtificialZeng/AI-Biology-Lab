#!/usr/bin/env python3
"""
Round 11 Evaluation Metrics (`PREREG_D`).

Strictly imports frozen mathematical primitives and R19/row-metric evaluators from
`round10/src/metrics_r10.py` and `round10/src/fct_pipeline.py` without mutating any
Round 10 file.

Round 11 extensions (`PREREG_D`):
  1. `summarize_boot_r11`: Two-sided `98.33%` Bonferroni CIs (`[0.8333%, 99.1667%]`,
     exact alpha/3 across co-primary H-D1, H-D2, H-D3) alongside `97.5%` and `95%` CIs.
  2. `compute_fct_perm100_null_r11`: 100-permutation null (`default_rng(k)`) computed cleanly
     under BOTH `all_nontarget` (primary `H-D2`) and `resp_top50` (secondary diagnostic).
  3. `compute_joint_paired_r_orth_diff_r11`: Within-mask paired joint bootstrap (`H-D3`:
     `FCT - E1_ConstMed` under `all_nontarget`, plus `resp_top50` diagnostic) with optional
     `gene_keep_mask` support.
  4. `compute_joint_cross_mask_e1_diff_r11`: Cross-mask paired joint bootstrap (`H-D4`:
     `E1_ConstMed(resp_top50) - E1_ConstMed(all_nontarget)`) recomputing each mask's
     orthogonalized split-half noise ceiling jointly inside each bootstrap draw.
  5. `compute_s7_cellcount_and_kd_r11`: `S7` cell-count & knockdown-depth control analysis
     parameterized cleanly for `all_nontarget` primary mask (and optional `resp_top50`).
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.stats import spearmanr

R10_SRC_DIR = Path(__file__).resolve().parents[2] / "round10" / "src"
if str(R10_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(R10_SRC_DIR))

from fct_pipeline import FctSourceState, compute_source_transfer_and_fc_bar, fct_from_fc_bar  # noqa: E402
from metrics_r10 import (  # noqa: E402
    EPS,
    _extract_masked_idx_per_pert,
    _moments,
    _moments_to_corr_and_slope,
    _orth_pair,
    compute_r19_variance_decomposition_r10,
    compute_row_metrics_r10,
    disattenuated_r,
    ols_slope_1d,
    pearson_1d_strict,
    spearman_brown_ceiling,
)

# Exact Bonferroni two-sided percentiles for 3 co-primary claims (alpha = 0.05 / 3)
PCT_LO_9833 = 100.0 * (0.05 / 6.0)
PCT_HI_9833 = 100.0 - PCT_LO_9833


def summarize_boot_r11(
    point: float,
    boots: np.ndarray,
    n_items: int,
    n_boot: int,
) -> Dict[str, Any]:
    """Summarize bootstrap replicates with `95%`, `97.5%`, and Bonferroni `98.33%` CIs."""
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
            "ci9833_lo": float("nan"),
            "ci9833_hi": float("nan"),
            "boot_se": float("nan"),
            "valid_frac": valid_frac,
            "n_items": int(n_items),
            "n_boot": int(n_boot),
        }
    q00833, q0125, q025, q975, q9875, q99167 = np.percentile(
        vals, [PCT_LO_9833, 1.25, 2.5, 97.5, 98.75, PCT_HI_9833]
    )
    return {
        "estimate": float(point),
        "ci95_lo": float(q025),
        "ci95_hi": float(q975),
        "ci975_lo": float(q0125),
        "ci975_hi": float(q9875),
        "ci9833_lo": float(q00833),
        "ci9833_hi": float(q99167),
        "boot_se": float(np.std(vals, ddof=1)) if vals.size > 1 else 0.0,
        "valid_frac": valid_frac,
        "n_items": int(n_items),
        "n_boot": int(n_boot),
    }


def paired_row_bootstrap_ci_r11(
    a: np.ndarray,
    b: np.ndarray,
    n_boot: int = 4000,
    seed: int = 20261123,
) -> Dict[str, Any]:
    """Perturbation-level paired bootstrap CI reporting `95%`, `97.5%`, and `98.33%` CIs."""
    xa = np.asarray(a, dtype=np.float64)
    xb = np.asarray(b, dtype=np.float64)
    n = xa.size
    if n == 0:
        return {"mean_a": float("nan"), "mean_b": float("nan"), "delta": float("nan")}
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(n_boot, n))
    boot_a = xa[idx].mean(axis=1)
    boot_b = xb[idx].mean(axis=1)
    boot_d = boot_a - boot_b

    pcts = [PCT_LO_9833, 1.25, 2.5, 97.5, 98.75, PCT_HI_9833]
    a_p = np.percentile(boot_a, pcts)
    b_p = np.percentile(boot_b, pcts)
    d_p = np.percentile(boot_d, pcts)
    mean_a = float(np.mean(xa))
    mean_b = float(np.mean(xb))
    return {
        "mean_a": mean_a,
        "mean_b": mean_b,
        "delta": float(mean_a - mean_b),
        "a_ci95_lo": float(a_p[2]),
        "a_ci95_hi": float(a_p[3]),
        "a_ci975_lo": float(a_p[1]),
        "a_ci975_hi": float(a_p[4]),
        "a_ci9833_lo": float(a_p[0]),
        "a_ci9833_hi": float(a_p[5]),
        "b_ci95_lo": float(b_p[2]),
        "b_ci95_hi": float(b_p[3]),
        "b_ci975_lo": float(b_p[1]),
        "b_ci975_hi": float(b_p[4]),
        "b_ci9833_lo": float(b_p[0]),
        "b_ci9833_hi": float(b_p[5]),
        "ci95_lo": float(d_p[2]),
        "ci95_hi": float(d_p[3]),
        "ci975_lo": float(d_p[1]),
        "ci975_hi": float(d_p[4]),
        "ci9833_lo": float(d_p[0]),
        "ci9833_hi": float(d_p[5]),
        "ci95_gt_zero": bool(d_p[2] > 0.0),
        "ci975_gt_zero": bool(d_p[1] > 0.0),
        "ci9833_gt_zero": bool(d_p[0] > 0.0),
        "n_rows": int(n),
    }


def _collect_mode_orth_moments(
    r_obs_full: np.ndarray,
    r_pred_list: Sequence[np.ndarray],
    t_mean_full: np.ndarray,
    rows: Sequence[int],
    idx_list: Sequence[np.ndarray],
    ha_res_full: Optional[np.ndarray],
    hb_res_full: Optional[np.ndarray],
    half_valid: Optional[np.ndarray],
) -> Tuple[List[Tuple[np.ndarray, ...]], Optional[Tuple[np.ndarray, ...]]]:
    """
    Collect per-perturbation orthogonalized residual moments for 1 or more predictors
    sharing identical `(rows, idx_list)`, plus orthogonalized split-half moments.
    """
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
        for k, rp in enumerate(r_pred_list):
            xo, yo = _orth_pair(xv, rp[p, idx], rv)
            xs_o[k].append(xo)
            ys_o[k].append(yo)

        if has_halves and ha_res_full is not None and hb_res_full is not None:
            hv_ok = half_valid is None or bool(half_valid[p])
            if hv_ok and np.all(np.isfinite(ha_res_full[p, idx])) and np.all(np.isfinite(hb_res_full[p, idx])):
                hao, hbo = _orth_pair(ha_res_full[p, idx], hb_res_full[p, idx], rv)
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


def _boot_disattenuated_orth_from_moments(
    pred_mom: Tuple[np.ndarray, ...],
    half_mom: Optional[Tuple[np.ndarray, ...]],
    draws: np.ndarray,
) -> Tuple[float, np.ndarray, float, float]:
    """Evaluate point estimate and joint bootstrap draws of `r_orth_cal` from moments."""
    nv, sx, sy, sxx, syy, sxy = pred_mom
    _, pt_r, _, _, _ = _moments_to_corr_and_slope(
        float(nv.sum()), float(sx.sum()), float(sy.sum()),
        float(sxx.sum()), float(syy.sum()), float(sxy.sum()),
    )
    _, b_r, _, _, _ = _moments_to_corr_and_slope(
        nv[draws].sum(axis=1), sx[draws].sum(axis=1), sy[draws].sum(axis=1),
        sxx[draws].sum(axis=1), syy[draws].sum(axis=1), sxy[draws].sum(axis=1),
    )

    n_boot = draws.shape[0]
    if half_mom is not None:
        hnv, hsx, hsy, hsxx, hsyy, hsxy = half_mom
        _, pt_rho, _, _, _ = _moments_to_corr_and_slope(
            float(hnv.sum()), float(hsx.sum()), float(hsy.sum()),
            float(hsxx.sum()), float(hsyy.sum()), float(hsxy.sum()),
        )
        pt_ceil = spearman_brown_ceiling(float(pt_rho))
        _, b_rho, _, _, _ = _moments_to_corr_and_slope(
            hnv[draws].sum(axis=1), hsx[draws].sum(axis=1), hsy[draws].sum(axis=1),
            hsxx[draws].sum(axis=1), hsyy[draws].sum(axis=1), hsxy[draws].sum(axis=1),
        )
        b_ceil = np.array([spearman_brown_ceiling(float(x)) for x in b_rho], dtype=np.float64)
    else:
        pt_rho = float("nan")
        pt_ceil = float("nan")
        b_ceil = np.full(n_boot, float("nan"), dtype=np.float64)

    pt_cal = disattenuated_r(float(pt_r), float(pt_ceil))
    b_cal = np.array([disattenuated_r(float(r), float(c)) for r, c in zip(b_r, b_ceil)], dtype=np.float64)
    return float(pt_cal), b_cal, float(pt_rho), float(pt_ceil)


def compute_residual_suite_r11(
    delta_obs: np.ndarray,
    delta_pred: np.ndarray,
    delta_transfer: np.ndarray,
    target_gene_idx: Sequence[int],
    delta_half_a: Optional[np.ndarray] = None,
    delta_half_b: Optional[np.ndarray] = None,
    half_valid: Optional[np.ndarray] = None,
    gene_keep_mask: Optional[np.ndarray] = None,
    top_k: int = 50,
    n_boot: int = 4000,
    seed: int = 20261123,
) -> Dict[str, Any]:
    """
    Compute Group D residual metrics (`r_orth_cal` + `95%`/`97.5%`/`98.33%` joint bootstrap CIs)
    under BOTH `all_nontarget` (primary) and `resp_top50` (secondary), preserving exact
    element-level parity with `compute_residual_suite_r10`.
    """
    r_obs_full = np.asarray(delta_obs, dtype=np.float64) - np.asarray(delta_transfer, dtype=np.float64)
    r_pred_full = np.asarray(delta_pred, dtype=np.float64) - np.asarray(delta_transfer, dtype=np.float64)
    t_mean_full = np.asarray(delta_transfer, dtype=np.float64)
    has_halves = delta_half_a is not None and delta_half_b is not None
    if has_halves:
        ha_res_full = np.asarray(delta_half_a, dtype=np.float64) - t_mean_full
        hb_res_full = np.asarray(delta_half_b, dtype=np.float64) - t_mean_full
    else:
        ha_res_full = hb_res_full = None

    rows_top50, idx_top50 = _extract_masked_idx_per_pert(
        delta_obs, target_gene_idx, gene_keep_mask=gene_keep_mask, mode="resp_top50", top_k=top_k
    )
    cos_vals: List[float] = []
    for p, idx in zip(rows_top50, idx_top50):
        u = delta_obs[p, idx].astype(np.float64)
        v = delta_pred[p, idx].astype(np.float64)
        den = float(np.linalg.norm(u) * np.linalg.norm(v))
        if den > EPS:
            cos_vals.append(float(np.clip(np.dot(u, v) / den, -1.0, 1.0)))

    out: Dict[str, Any] = {
        "cos_agreement_top50": float(np.mean(cos_vals)) if cos_vals else float("nan"),
    }

    for mode in ("resp_top50", "all_nontarget"):
        rows, idx_list = _extract_masked_idx_per_pert(
            delta_obs,
            target_gene_idx,
            extra_valid=r_pred_full,
            gene_keep_mask=gene_keep_mask,
            mode=mode,
            top_k=top_k,
        )
        m = len(rows)
        if m == 0:
            out[mode] = {"n_perts": 0}
            continue

        xs_raw, ys_raw, xs_orth, ys_orth = [], [], [], []
        ha_raw_l, hb_raw_l, ha_orth_l, hb_orth_l = [], [], [], []
        any_hv = False
        per_pert_slopes: List[float] = []
        per_pert_corrs: List[float] = []

        for p, idx in zip(rows, idx_list):
            xv = r_obs_full[p, idx]
            yv = r_pred_full[p, idx]
            rv = t_mean_full[p, idx]
            xov, yov = _orth_pair(xv, yv, rv)
            xs_raw.append(xv)
            ys_raw.append(yv)
            xs_orth.append(xov)
            ys_orth.append(yov)
            per_pert_slopes.append(ols_slope_1d(xv, yv))
            per_pert_corrs.append(pearson_1d_strict(xv, yv))

            if has_halves and ha_res_full is not None and hb_res_full is not None:
                hv_ok = half_valid is None or bool(half_valid[p])
                if hv_ok and np.all(np.isfinite(ha_res_full[p, idx])) and np.all(np.isfinite(hb_res_full[p, idx])):
                    hav = ha_res_full[p, idx]
                    hbv = hb_res_full[p, idx]
                    hao, hbo = _orth_pair(hav, hbv, rv)
                    ha_raw_l.append(hav)
                    hb_raw_l.append(hbv)
                    ha_orth_l.append(hao)
                    hb_orth_l.append(hbo)
                    any_hv = True
                else:
                    z = np.zeros(0, dtype=np.float64)
                    ha_raw_l.append(z)
                    hb_raw_l.append(z)
                    ha_orth_l.append(z)
                    hb_orth_l.append(z)

        nv, sx, sy, sxx, syy, sxy = _moments(xs_raw, ys_raw)
        pt_slope, pt_r_raw, pt_sd_ratio, pt_rev, pt_dem = _moments_to_corr_and_slope(
            float(nv.sum()), float(sx.sum()), float(sy.sum()), float(sxx.sum()), float(syy.sum()), float(sxy.sum())
        )
        onv, osx, osy, osxx, osyy, osxy = _moments(xs_orth, ys_orth)
        _, pt_r_orth, _, _, _ = _moments_to_corr_and_slope(
            float(onv.sum()), float(osx.sum()), float(osy.sum()), float(osxx.sum()), float(osyy.sum()), float(osxy.sum())
        )

        rng = np.random.default_rng(seed)
        draws = rng.integers(0, m, size=(n_boot, m))
        b_slope, b_r_raw, b_sd_ratio, b_rev, b_dem = _moments_to_corr_and_slope(
            nv[draws].sum(axis=1), sx[draws].sum(axis=1), sy[draws].sum(axis=1),
            sxx[draws].sum(axis=1), syy[draws].sum(axis=1), sxy[draws].sum(axis=1),
        )
        _, b_r_orth, _, _, _ = _moments_to_corr_and_slope(
            onv[draws].sum(axis=1), osx[draws].sum(axis=1), osy[draws].sum(axis=1),
            osxx[draws].sum(axis=1), osyy[draws].sum(axis=1), osxy[draws].sum(axis=1),
        )

        if any_hv:
            hnv, hsx, hsy, hsxx, hsyy, hsxy = _moments(ha_raw_l, hb_raw_l)
            _, pt_rho_ab, _, _, _ = _moments_to_corr_and_slope(
                float(hnv.sum()), float(hsx.sum()), float(hsy.sum()),
                float(hsxx.sum()), float(hsyy.sum()), float(hsxy.sum()),
            )
            pt_ceil = spearman_brown_ceiling(float(pt_rho_ab))
            _, b_rho_ab, _, _, _ = _moments_to_corr_and_slope(
                hnv[draws].sum(axis=1), hsx[draws].sum(axis=1), hsy[draws].sum(axis=1),
                hsxx[draws].sum(axis=1), hsyy[draws].sum(axis=1), hsxy[draws].sum(axis=1),
            )
            b_ceil = np.array([spearman_brown_ceiling(float(x)) for x in b_rho_ab], dtype=np.float64)

            ohnv, ohsx, ohsy, ohsxx, ohsyy, ohsxy = _moments(ha_orth_l, hb_orth_l)
            _, pt_rho_ab_orth, _, _, _ = _moments_to_corr_and_slope(
                float(ohnv.sum()), float(ohsx.sum()), float(ohsy.sum()),
                float(ohsxx.sum()), float(ohsyy.sum()), float(ohsxy.sum()),
            )
            pt_ceil_orth = spearman_brown_ceiling(float(pt_rho_ab_orth))
            _, b_rho_ab_orth, _, _, _ = _moments_to_corr_and_slope(
                ohnv[draws].sum(axis=1), ohsx[draws].sum(axis=1), ohsy[draws].sum(axis=1),
                ohsxx[draws].sum(axis=1), ohsyy[draws].sum(axis=1), ohsxy[draws].sum(axis=1),
            )
            b_ceil_orth = np.array([spearman_brown_ceiling(float(x)) for x in b_rho_ab_orth], dtype=np.float64)
        else:
            pt_rho_ab = pt_ceil = float("nan")
            pt_rho_ab_orth = pt_ceil_orth = float("nan")
            b_ceil = np.full(n_boot, float("nan"), dtype=np.float64)
            b_ceil_orth = np.full(n_boot, float("nan"), dtype=np.float64)

        pt_cal = disattenuated_r(float(pt_r_raw), float(pt_ceil))
        b_cal = np.array([disattenuated_r(float(r), float(c)) for r, c in zip(b_r_raw, b_ceil)], dtype=np.float64)
        pt_orth_cal = disattenuated_r(float(pt_r_orth), float(pt_ceil_orth))
        b_orth_cal = np.array(
            [disattenuated_r(float(r), float(c)) for r, c in zip(b_r_orth, b_ceil_orth)],
            dtype=np.float64,
        )

        pt_orth_cal_sq = float(pt_orth_cal * pt_orth_cal) if np.isfinite(pt_orth_cal) else float("nan")
        b_orth_cal_sq = np.where(np.isfinite(b_orth_cal), b_orth_cal * b_orth_cal, np.nan)
        pt_cal_sq = float(pt_cal * pt_cal) if np.isfinite(pt_cal) else float("nan")
        b_cal_sq = np.where(np.isfinite(b_cal), b_cal * b_cal, np.nan)

        finite_slopes = [x for x in per_pert_slopes if np.isfinite(x)]
        finite_corrs = [x for x in per_pert_corrs if np.isfinite(x)]
        out[mode] = {
            "n_perts": int(m),
            "r_orth_cal": summarize_boot_r11(pt_orth_cal, b_orth_cal, m, n_boot),
            "r_orth_cal_sq": summarize_boot_r11(pt_orth_cal_sq, b_orth_cal_sq, m, n_boot),
            "r_orth_raw": summarize_boot_r11(float(pt_r_orth), b_r_orth, m, n_boot),
            "r_cal": summarize_boot_r11(pt_cal, b_cal, m, n_boot),
            "r_cal_sq": summarize_boot_r11(pt_cal_sq, b_cal_sq, m, n_boot),
            "r_raw": summarize_boot_r11(float(pt_r_raw), b_r_raw, m, n_boot),
            "slope_ols": summarize_boot_r11(float(pt_slope), b_slope, m, n_boot),
            "sd_pred_over_sd_obs": summarize_boot_r11(float(pt_sd_ratio), b_sd_ratio, m, n_boot),
            "slope_deming": summarize_boot_r11(float(pt_dem), b_dem, m, n_boot),
            "slope_obs_on_pred": summarize_boot_r11(float(pt_rev), b_rev, m, n_boot),
            "per_pert_slope_mean": float(np.mean(finite_slopes)) if finite_slopes else float("nan"),
            "per_pert_corr_mean": float(np.mean(finite_corrs)) if finite_corrs else float("nan"),
            "rho_half_ab": float(pt_rho_ab),
            "ceiling_r": float(pt_ceil),
            "rho_half_ab_orth": float(pt_rho_ab_orth),
            "ceiling_r_orth": float(pt_ceil_orth),
        }

    return out


def compute_joint_paired_r_orth_diff_r11(
    delta_obs: np.ndarray,
    delta_pred_a: np.ndarray,
    delta_pred_b: np.ndarray,
    delta_transfer: np.ndarray,
    target_gene_idx: Sequence[int],
    delta_half_a: Optional[np.ndarray],
    delta_half_b: Optional[np.ndarray],
    half_valid: Optional[np.ndarray] = None,
    gene_keep_mask: Optional[np.ndarray] = None,
    mode: str = "all_nontarget",
    top_k: int = 50,
    n_boot: int = 4000,
    seed: int = 20261123,
) -> Dict[str, Any]:
    """
    Within-mask paired joint-bootstrap CI for `r_orth_cal(A) - r_orth_cal(B)`
    (primary `H-D3`: `FCT - E1_ConstMed` under `all_nontarget`).
    """
    t_mean = np.asarray(delta_transfer, dtype=np.float64)
    r_obs = np.asarray(delta_obs, dtype=np.float64) - t_mean
    r_a = np.asarray(delta_pred_a, dtype=np.float64) - t_mean
    r_b = np.asarray(delta_pred_b, dtype=np.float64) - t_mean
    ha_res = (np.asarray(delta_half_a, dtype=np.float64) - t_mean) if delta_half_a is not None else None
    hb_res = (np.asarray(delta_half_b, dtype=np.float64) - t_mean) if delta_half_b is not None else None

    rows, idx_list = _extract_masked_idx_per_pert(
        delta_obs,
        target_gene_idx,
        extra_valid=r_a + r_b,
        gene_keep_mask=gene_keep_mask,
        mode=mode,
        top_k=top_k,
    )
    m = len(rows)
    if m == 0:
        return summarize_boot_r11(float("nan"), np.zeros(0), 0, n_boot)

    pred_moms, half_mom = _collect_mode_orth_moments(
        r_obs, [r_a, r_b], t_mean, rows, idx_list, ha_res, hb_res, half_valid
    )
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, m, size=(n_boot, m))
    pt_a, b_a, _, _ = _boot_disattenuated_orth_from_moments(pred_moms[0], half_mom, draws)
    pt_b, b_b, _, _ = _boot_disattenuated_orth_from_moments(pred_moms[1], half_mom, draws)
    return summarize_boot_r11(pt_a - pt_b, b_a - b_b, m, n_boot)


def compute_joint_cross_mask_e1_diff_r11(
    delta_obs: np.ndarray,
    delta_pred_e1: np.ndarray,
    delta_transfer: np.ndarray,
    target_gene_idx: Sequence[int],
    delta_half_a: Optional[np.ndarray],
    delta_half_b: Optional[np.ndarray],
    half_valid: Optional[np.ndarray] = None,
    gene_keep_mask: Optional[np.ndarray] = None,
    top_k: int = 50,
    n_boot: int = 4000,
    seed: int = 20261123,
) -> Dict[str, Any]:
    """
    Cross-mask paired joint-bootstrap CI for `H-D4`:
      `r_orth_cal(E1_ConstMed, mask='resp_top50') - r_orth_cal(E1_ConstMed, mask='all_nontarget')`.

    Recomputes each mask's orthogonalized residuals AND each mask's split-half noise
    ceiling jointly on the exact same perturbation resample indices `draws`.
    """
    t_mean = np.asarray(delta_transfer, dtype=np.float64)
    r_obs = np.asarray(delta_obs, dtype=np.float64) - t_mean
    r_e1 = np.asarray(delta_pred_e1, dtype=np.float64) - t_mean
    ha_res = (np.asarray(delta_half_a, dtype=np.float64) - t_mean) if delta_half_a is not None else None
    hb_res = (np.asarray(delta_half_b, dtype=np.float64) - t_mean) if delta_half_b is not None else None

    rows_top50, idx_top50_raw = _extract_masked_idx_per_pert(
        delta_obs, target_gene_idx, extra_valid=r_e1, gene_keep_mask=gene_keep_mask, mode="resp_top50", top_k=top_k
    )
    rows_all, idx_all_raw = _extract_masked_idx_per_pert(
        delta_obs, target_gene_idx, extra_valid=r_e1, gene_keep_mask=gene_keep_mask, mode="all_nontarget", top_k=top_k
    )

    map_top50 = dict(zip(rows_top50, idx_top50_raw))
    map_all = dict(zip(rows_all, idx_all_raw))
    rows_common = [p for p in rows_all if p in map_top50]
    m = len(rows_common)
    if m == 0:
        return summarize_boot_r11(float("nan"), np.zeros(0), 0, n_boot)

    idx_top50 = [map_top50[p] for p in rows_common]
    idx_all = [map_all[p] for p in rows_common]

    pred_moms_a, half_mom_a = _collect_mode_orth_moments(
        r_obs, [r_e1], t_mean, rows_common, idx_top50, ha_res, hb_res, half_valid
    )
    pred_moms_b, half_mom_b = _collect_mode_orth_moments(
        r_obs, [r_e1], t_mean, rows_common, idx_all, ha_res, hb_res, half_valid
    )

    rng = np.random.default_rng(seed)
    draws = rng.integers(0, m, size=(n_boot, m))
    pt_a, b_a, _, _ = _boot_disattenuated_orth_from_moments(pred_moms_a[0], half_mom_a, draws)
    pt_b, b_b, _, _ = _boot_disattenuated_orth_from_moments(pred_moms_b[0], half_mom_b, draws)

    res = summarize_boot_r11(pt_a - pt_b, b_a - b_b, m, n_boot)
    res["mean_a"] = pt_a
    res["mean_b"] = pt_b
    res["a_summary"] = summarize_boot_r11(pt_a, b_a, m, n_boot)
    res["b_summary"] = summarize_boot_r11(pt_b, b_b, m, n_boot)
    return res


def compute_fct_perm100_null_r11(
    state: FctSourceState,
    target_ctrl_mean_libnorm: np.ndarray,
    query_target_names: Sequence[str],
    delta_obs: np.ndarray,
    delta_transfer: np.ndarray,
    target_gene_idx: Sequence[int],
    ceil_r: float,
    ceil_r_orth: float,
    fc_bundle: Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]] = None,
    mode: str = "all_nontarget",
    n_perm: int = 100,
    top_k: int = 50,
    row_subset_mask: Optional[np.ndarray] = None,
    gene_keep_mask: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """
    Vectorized 100-permutation null distribution (`k = 0..n_perm-1`, `default_rng(k).permutation(G)`)
    supporting BOTH `all_nontarget` (primary `H-D2`) and `resp_top50` (secondary), plus `row_subset_mask`
    (`S7` tertiles) and `gene_keep_mask` (`S5`).
    """
    mu_t = np.maximum(np.asarray(target_ctrl_mean_libnorm, dtype=np.float64), 0.0)
    G = mu_t.shape[0]
    if fc_bundle is None:
        fc_bundle = compute_source_transfer_and_fc_bar(state, query_target_names)
    pred_transfer, fc_bar, has_hits = fc_bundle

    if row_subset_mask is not None:
        sub_idx = np.flatnonzero(row_subset_mask)
        d_obs_s = delta_obs[sub_idx]
        t_mean_s = delta_transfer[sub_idx]
        t_gidx_s = np.asarray(target_gene_idx, dtype=np.int64)[sub_idx]
        fc_bar_s = fc_bar[sub_idx]
        p_trans_s = pred_transfer[sub_idx]
        hits_s = has_hits[sub_idx]
    else:
        d_obs_s = delta_obs
        t_mean_s = delta_transfer
        t_gidx_s = np.asarray(target_gene_idx, dtype=np.int64)
        fc_bar_s = fc_bar
        p_trans_s = pred_transfer
        hits_s = has_hits

    r_obs = np.asarray(d_obs_s, dtype=np.float64) - np.asarray(t_mean_s, dtype=np.float64)
    rows, idx_list = _extract_masked_idx_per_pert(
        d_obs_s,
        t_gidx_s,
        extra_valid=np.asarray(t_mean_s, dtype=np.float64),
        gene_keep_mask=gene_keep_mask,
        mode=mode,
        top_k=top_k,
    )
    if not rows:
        return {"n_perm": int(n_perm), "mode": mode, "r_orth_cal": {}, "r_cal": {}, "slope_ols": {}}

    rows_arr = np.asarray(rows, dtype=np.int64)
    m = len(rows)
    M = np.zeros((m, G), dtype=np.float64)
    B = np.zeros((m, G), dtype=np.float64)
    for ii, (p, idx) in enumerate(zip(rows, idx_list)):
        M[ii, idx] = 1.0
        rv = np.asarray(t_mean_s[p, idx], dtype=np.float64)
        nrm = float(np.linalg.norm(rv))
        if nrm > EPS:
            B[ii, idx] = rv / nrm

    x_mat_raw = r_obs[rows_arr] * M
    proj_x = np.sum(x_mat_raw * B, axis=1, keepdims=True)
    x_mat_orth = x_mat_raw - proj_x * B

    mask_bool = M > 0.5
    x_raw = x_mat_raw[mask_bool]
    x_orth = x_mat_orth[mask_bool]

    fc_bar_eval = fc_bar_s[rows_arr]
    p_trans_eval = p_trans_s[rows_arr]
    hits_eval = hits_s[rows_arr]
    t_mean_eval = np.asarray(t_mean_s[rows_arr], dtype=np.float64)

    r_orth_cal_vals = np.empty(n_perm, dtype=np.float64)
    r_cal_vals = np.empty(n_perm, dtype=np.float64)
    slope_vals = np.empty(n_perm, dtype=np.float64)

    for k in range(n_perm):
        perm = np.random.default_rng(k).permutation(G)
        pred_k = fct_from_fc_bar(mu_t[perm], fc_bar_eval, p_trans_eval, hits_eval)
        res_k = (pred_k - t_mean_eval) * M
        proj_y = np.sum(res_k * B, axis=1, keepdims=True)
        res_k_orth = res_k - proj_y * B

        y_raw = res_k[mask_bool]
        y_orth = res_k_orth[mask_bool]
        slope_vals[k] = ols_slope_1d(x_raw, y_raw)
        r_cal_vals[k] = disattenuated_r(pearson_1d_strict(x_raw, y_raw), ceil_r)
        r_orth_cal_vals[k] = disattenuated_r(pearson_1d_strict(x_orth, y_orth), ceil_r_orth)

    def _pct(arr: np.ndarray) -> Dict[str, float]:
        fin = arr[np.isfinite(arr)]
        if fin.size == 0:
            return {
                k: float("nan")
                for k in ("p0_833", "p2_5", "median", "p97_5", "p98_75", "p99_167", "max", "mean")
            }
        p00833, p25, p50, p975, p9875, p99167 = np.percentile(
            fin, [PCT_LO_9833, 2.5, 50.0, 97.5, 98.75, PCT_HI_9833]
        )
        return {
            "p0_833": float(p00833),
            "p2_5": float(p25),
            "median": float(p50),
            "p97_5": float(p975),
            "p98_75": float(p9875),
            "p99_167": float(p99167),
            "max": float(np.max(fin)),
            "mean": float(np.mean(fin)),
        }

    return {
        "n_perm": int(n_perm),
        "mode": mode,
        "r_orth_cal": _pct(r_orth_cal_vals),
        "r_cal": _pct(r_cal_vals),
        "slope_ols": _pct(slope_vals),
        "perm0_r_orth_cal": float(r_orth_cal_vals[0]),
        "perm0_r_cal": float(r_cal_vals[0]),
        "perm0_slope_ols": float(slope_vals[0]),
    }


def compute_s7_cellcount_and_kd_r11(
    state: FctSourceState,
    query_target_names: Sequence[str],
    target_ctrl_mean_libnorm: np.ndarray,
    fc_bundle: Tuple[np.ndarray, np.ndarray, np.ndarray],
    delta_obs: np.ndarray,
    delta_pred_fct: np.ndarray,
    delta_pred_e1: np.ndarray,
    delta_transfer: np.ndarray,
    target_gene_idx: Sequence[int],
    n_cells: Optional[np.ndarray],
    knockdown_depth: Optional[np.ndarray],
    delta_half_a: Optional[np.ndarray] = None,
    delta_half_b: Optional[np.ndarray] = None,
    half_valid: Optional[np.ndarray] = None,
    mismatch_preds: Optional[Dict[str, np.ndarray]] = None,
    gene_keep_mask: Optional[np.ndarray] = None,
    mode: str = "all_nontarget",
    top_k: int = 50,
    n_boot: int = 4000,
    seed: int = 20261123,
) -> Dict[str, Any]:
    """
    Round 11 `S7` Cell-Count & Knockdown-Depth Confounding Controls under `mode` (`all_nontarget` primary):
      (i)   Per-pert `r_orth_p` across masked genes
      (ii)  Spearman correlation with `ln(n_cells)` and `|k_p|` (+ 97.5% / 95% / 98.33% CIs)
      (iii) Mean-centered OLS `[1, ln(n_cells) - mean, |k_p| - mean]`
      (iv)  Balanced `n_cells` tertiles (`T1_low`, `T2_mid`, `T3_high`) recomputing pooled FCT, E1,
            mismatched controls, perm100, and paired `FCT - E1` under `mode` (`97.5%` CIs per PREREG_D v2).
    """
    N, G = delta_obs.shape
    t_mean = np.asarray(delta_transfer, dtype=np.float64)
    r_obs = np.asarray(delta_obs, dtype=np.float64) - t_mean
    r_pred = np.asarray(delta_pred_fct, dtype=np.float64) - t_mean

    r_orth_p = np.full(N, np.nan, dtype=np.float64)
    kd_arr = (
        np.abs(np.asarray(knockdown_depth, dtype=np.float64))
        if knockdown_depth is not None
        else np.full(N, np.nan, dtype=np.float64)
    )
    nc_arr = (
        np.asarray(n_cells, dtype=np.float64)
        if n_cells is not None
        else np.full(N, np.nan, dtype=np.float64)
    )
    ln_nc_arr = np.where(np.isfinite(nc_arr) & (nc_arr > 0.0), np.log(nc_arr), np.nan)

    rows, idx_list = _extract_masked_idx_per_pert(
        delta_obs,
        target_gene_idx,
        extra_valid=r_pred,
        gene_keep_mask=gene_keep_mask,
        mode=mode,
        top_k=top_k,
    )
    for p, idx in zip(rows, idx_list):
        xo, yo = _orth_pair(r_obs[p, idx], r_pred[p, idx], t_mean[p, idx])
        r_orth_p[p] = pearson_1d_strict(xo, yo)

    def _boot_spearman(x_vec: np.ndarray, y_vec: np.ndarray) -> Dict[str, Any]:
        ok = np.isfinite(x_vec) & np.isfinite(y_vec)
        xv = x_vec[ok]
        yv = y_vec[ok]
        m = xv.size
        if m < 5 or np.std(xv) <= EPS or np.std(yv) <= EPS:
            return summarize_boot_r11(float("nan"), np.zeros(0), int(m), n_boot)
        pt = float(spearmanr(xv, yv).statistic)
        rng = np.random.default_rng(seed)
        b_vals = np.empty(n_boot, dtype=np.float64)
        for b in range(n_boot):
            idx_b = rng.integers(0, m, size=m)
            xb, yb = xv[idx_b], yv[idx_b]
            if np.std(xb) <= EPS or np.std(yb) <= EPS:
                b_vals[b] = float("nan")
            else:
                b_vals[b] = float(spearmanr(xb, yb).statistic)
        return summarize_boot_r11(pt, b_vals, int(m), n_boot)

    spearman_ln_nc = _boot_spearman(r_orth_p, ln_nc_arr)
    spearman_kd = _boot_spearman(r_orth_p, kd_arr)

    def _boot_centered_ols(y_vec: np.ndarray, cols: List[np.ndarray], col_names: List[str]) -> Dict[str, Any]:
        ok = np.isfinite(y_vec)
        for c in cols:
            ok = ok & np.isfinite(c)
        yv = y_vec[ok]
        X_raw = [c[ok] for c in cols]
        m = yv.size
        if m < len(cols) + 3:
            return {"n_valid": int(m), "status": "insufficient_valid_rows"}
        means = [float(np.mean(c)) for c in X_raw]
        X_centered = np.column_stack([np.ones(m, dtype=np.float64)] + [c - mu for c, mu in zip(X_raw, means)])
        if np.linalg.matrix_rank(X_centered) < X_centered.shape[1]:
            return {"n_valid": int(m), "status": "rank_deficient"}
        beta_pt, _, _, _ = np.linalg.lstsq(X_centered, yv, rcond=None)
        rng = np.random.default_rng(seed)
        draws = rng.integers(0, m, size=(n_boot, m))
        n_par = X_centered.shape[1]
        b_betas = np.full((n_boot, n_par), np.nan, dtype=np.float64)
        for b in range(n_boot):
            di = draws[b]
            try:
                b_betas[b], _, _, _ = np.linalg.lstsq(X_centered[di], yv[di], rcond=None)
            except np.linalg.LinAlgError:
                pass
        return {
            "n_valid": int(m),
            "covariate_means": {name: mu for name, mu in zip(col_names, means)},
            "unadjusted_mean_r_orth_p": summarize_boot_r11(float(np.mean(yv)), yv[draws].mean(axis=1), int(m), n_boot),
            "adjusted_intercept": summarize_boot_r11(float(beta_pt[0]), b_betas[:, 0], int(m), n_boot),
            "slopes": {
                name: summarize_boot_r11(float(beta_pt[j + 1]), b_betas[:, j + 1], int(m), n_boot)
                for j, name in enumerate(col_names)
            },
        }

    ols_both = _boot_centered_ols(r_orth_p, [ln_nc_arr, kd_arr], ["ln_n_cells", "abs_kd_depth"])
    ols_nc_only = _boot_centered_ols(r_orth_p, [ln_nc_arr], ["ln_n_cells"])

    tertiles_out: Dict[str, Any] = {}
    valid_nc = np.isfinite(nc_arr) & (nc_arr > 0.0)
    if int(np.sum(valid_nc)) >= 15:
        v_indices = np.flatnonzero(valid_nc)
        order = sorted(v_indices.tolist(), key=lambda i: (float(nc_arr[i]), int(i)))
        chunks = np.array_split(np.asarray(order, dtype=np.int64), 3)
        labels = ("T1_low", "T2_mid", "T3_high")
        tg_arr = np.asarray(target_gene_idx, dtype=np.int64)
        for lbl, c_idx in zip(labels, chunks):
            if c_idx.size < 5:
                continue
            mask_t = np.zeros(N, dtype=bool)
            mask_t[c_idx] = True
            fct_res = compute_residual_suite_r11(
                delta_obs=delta_obs[c_idx],
                delta_pred=delta_pred_fct[c_idx],
                delta_transfer=delta_transfer[c_idx],
                target_gene_idx=tg_arr[c_idx],
                delta_half_a=delta_half_a[c_idx] if delta_half_a is not None else None,
                delta_half_b=delta_half_b[c_idx] if delta_half_b is not None else None,
                half_valid=half_valid[c_idx] if half_valid is not None else None,
                gene_keep_mask=gene_keep_mask,
                top_k=top_k,
                n_boot=n_boot,
                seed=seed,
            )
            md = fct_res.get(mode, {})
            max_mis_hi = float("nan")
            mis_detail: Dict[str, float] = {}
            if mismatch_preds:
                hi_list = []
                for m_name, m_mat in mismatch_preds.items():
                    m_res = compute_residual_suite_r11(
                        delta_obs=delta_obs[c_idx],
                        delta_pred=m_mat[c_idx],
                        delta_transfer=delta_transfer[c_idx],
                        target_gene_idx=tg_arr[c_idx],
                        delta_half_a=delta_half_a[c_idx] if delta_half_a is not None else None,
                        delta_half_b=delta_half_b[c_idx] if delta_half_b is not None else None,
                        half_valid=half_valid[c_idx] if half_valid is not None else None,
                        gene_keep_mask=gene_keep_mask,
                        top_k=top_k,
                        n_boot=n_boot,
                        seed=seed,
                    )
                    hi_val = float(m_res.get(mode, {}).get("r_orth_cal", {}).get("ci975_hi", float("nan")))
                    mis_detail[m_name] = hi_val
                    if np.isfinite(hi_val):
                        hi_list.append(hi_val)
                if hi_list:
                    max_mis_hi = float(max(hi_list))

            ceil_raw = float(md.get("ceiling_r", float("nan")))
            ceil_orth = float(md.get("ceiling_r_orth", float("nan")))
            perm100 = compute_fct_perm100_null_r11(
                state=state,
                target_ctrl_mean_libnorm=target_ctrl_mean_libnorm,
                query_target_names=query_target_names,
                delta_obs=delta_obs,
                delta_transfer=delta_transfer,
                target_gene_idx=target_gene_idx,
                ceil_r=ceil_raw,
                ceil_r_orth=ceil_orth,
                fc_bundle=fc_bundle,
                mode=mode,
                n_perm=100,
                top_k=top_k,
                row_subset_mask=mask_t,
                gene_keep_mask=gene_keep_mask,
            )
            diff_hd3 = (
                compute_joint_paired_r_orth_diff_r11(
                    delta_obs=delta_obs[c_idx],
                    delta_pred_a=delta_pred_fct[c_idx],
                    delta_pred_b=delta_pred_e1[c_idx],
                    delta_transfer=delta_transfer[c_idx],
                    target_gene_idx=tg_arr[c_idx],
                    delta_half_a=delta_half_a[c_idx],
                    delta_half_b=delta_half_b[c_idx],
                    half_valid=half_valid[c_idx] if half_valid is not None else None,
                    gene_keep_mask=gene_keep_mask,
                    mode=mode,
                    top_k=top_k,
                    n_boot=n_boot,
                    seed=seed,
                )
                if (delta_half_a is not None and delta_half_b is not None)
                else {}
            )

            fct_lo975 = float(md.get("r_orth_cal", {}).get("ci975_lo", float("nan")))
            perm_p975 = float(perm100.get("r_orth_cal", {}).get("p97_5", float("nan")))
            hd3_lo975 = float(diff_hd3.get("ci975_lo", float("nan")))
            beats_mis = bool(np.isfinite(fct_lo975) and np.isfinite(max_mis_hi) and fct_lo975 > max_mis_hi)
            beats_perm = bool(np.isfinite(fct_lo975) and np.isfinite(perm_p975) and fct_lo975 > perm_p975)
            passes_hd3 = bool(np.isfinite(hd3_lo975) and hd3_lo975 > 0.0)
            tertiles_out[lbl] = {
                "n_targets": int(c_idx.size),
                "n_cells_min": float(np.min(nc_arr[c_idx])),
                "n_cells_median": float(np.median(nc_arr[c_idx])),
                "n_cells_max": float(np.max(nc_arr[c_idx])),
                "fct_r_orth_cal": md.get("r_orth_cal", {}),
                "fct_slope_ols": md.get("slope_ols", {}),
                "rho_half_ab_orth": float(md.get("rho_half_ab_orth", float("nan"))),
                "ceiling_r_orth": ceil_orth,
                "max_mismatched_ci975_hi": max_mis_hi,
                "perm100_p975": perm_p975,
                "perm100_median": float(perm100.get("r_orth_cal", {}).get("median", float("nan"))),
                "fct_minus_e1": diff_hd3,
                "beats_mismatched_ci975": beats_mis,
                "beats_perm100_p975": beats_perm,
                "passes_hd1": bool(np.isfinite(fct_lo975) and fct_lo975 > 0.0),
                "passes_hd2": bool(beats_mis and beats_perm),
                "passes_hd3": passes_hd3,
                "mismatched_ci975_hi_by_source": mis_detail,
            }

    return {
        "mode": mode,
        "n_evaluated": int(N),
        "n_finite_r_orth_p": int(np.sum(np.isfinite(r_orth_p))),
        "n_finite_n_cells": int(np.sum(np.isfinite(ln_nc_arr))),
        "n_finite_abs_kd_depth": int(np.sum(np.isfinite(kd_arr))),
        "per_pert_summary": {
            "mean": float(np.nanmean(r_orth_p)),
            "median": float(np.nanmedian(r_orth_p)),
            "sd": float(np.nanstd(r_orth_p, ddof=1)) if int(np.sum(np.isfinite(r_orth_p))) > 1 else float("nan"),
        },
        "spearman_r_orth_p_vs_ln_ncells": spearman_ln_nc,
        "spearman_r_orth_p_vs_abs_kd": spearman_kd,
        "ols_adjusted_both": ols_both,
        "ols_adjusted_ln_ncells_only": ols_nc_only,
        "n_cells_tertiles": tertiles_out,
    }


def single_half_ceiling(rho_ab: float) -> float:
    """Single-half noise ceiling `sqrt(rho_AB)` (`PREREG_D §5 SA-8`); NaN if `rho_ab <= 0`."""
    if not np.isfinite(rho_ab) or rho_ab <= 0.0:
        return float("nan")
    return float(np.sqrt(float(rho_ab)))


def _prepare_split_direction_indices(
    delta_sel: np.ndarray,
    delta_eval_obs: np.ndarray,
    delta_other_obs: np.ndarray,
    delta_transfer: np.ndarray,
    target_gene_idx: Sequence[int],
    half_valid: Optional[np.ndarray],
    extra_valid: Optional[np.ndarray] = None,
    gene_keep_mask: Optional[np.ndarray] = None,
    top_k: int = 50,
) -> Tuple[List[int], List[np.ndarray]]:
    """
    Extract per-perturbation top-k gene indices ranked by `|delta_sel[p]|` (non-self, finite
    in selection half, evaluation half, other half for noise ceiling, transfer mean, and optional predictor).
    """
    N, G = delta_sel.shape
    rows: List[int] = []
    idx_list: List[np.ndarray] = []
    for p in range(N):
        if half_valid is not None and not bool(half_valid[p]):
            continue
        valid = (
            np.isfinite(delta_sel[p])
            & np.isfinite(delta_eval_obs[p])
            & np.isfinite(delta_other_obs[p])
            & np.isfinite(delta_transfer[p])
        )
        if extra_valid is not None:
            valid = valid & np.isfinite(extra_valid[p])
        if gene_keep_mask is not None:
            valid = valid & gene_keep_mask
        ti = int(target_gene_idx[p])
        if 0 <= ti < G:
            valid[ti] = False
        idx = np.flatnonzero(valid)
        if idx.size == 0:
            continue
        if idx.size > top_k:
            idx = idx[np.argsort(np.abs(delta_sel[p, idx]), kind="mergesort")[-top_k:]]
        rows.append(p)
        idx_list.append(idx)
    return rows, idx_list


def compute_resp_top50_split_suite_r11(
    delta_pred: np.ndarray,
    delta_transfer: np.ndarray,
    target_gene_idx: Sequence[int],
    delta_half_a: Optional[np.ndarray],
    delta_half_b: Optional[np.ndarray],
    half_valid: Optional[np.ndarray] = None,
    gene_keep_mask: Optional[np.ndarray] = None,
    top_k: int = 50,
    n_boot: int = 4000,
    seed: int = 20261123,
) -> Dict[str, Any]:
    """
    `PREREG_D §5 SA-8` Split-Half Out-of-Sample Top-50 Mask (`resp_top50_split`):
      - Direction `A -> B`: rank top-50 genes by `|HalfA_delta|`, evaluate `pred` residual
        against `HalfB_delta` residual (orthogonalized to `Transfer_Mean` on `M^{(A)}`),
        disattenuated by single-half noise ceiling `sqrt(rho_AB_orth(M^{(A)}))`.
      - Direction `B -> A`: symmetric reverse.
      - Symmetric averaged `r_orth_raw_split` and `r_orth_cal_split` with joint perturbation
        bootstrap (`default_rng(seed)` x `n_boot`, `97.5%` CIs).
    """
    if delta_half_a is None or delta_half_b is None:
        return {"n_perts": 0, "status": "missing_split_halves"}

    ha = np.asarray(delta_half_a, dtype=np.float64)
    hb = np.asarray(delta_half_b, dtype=np.float64)
    tm = np.asarray(delta_transfer, dtype=np.float64)
    rp = np.asarray(delta_pred, dtype=np.float64) - tm
    ha_res = ha - tm
    hb_res = hb - tm

    rows_ab, idx_ab = _prepare_split_direction_indices(
        ha, hb, ha, tm, target_gene_idx, half_valid, extra_valid=rp, gene_keep_mask=gene_keep_mask, top_k=top_k
    )
    rows_ba, idx_ba = _prepare_split_direction_indices(
        hb, ha, hb, tm, target_gene_idx, half_valid, extra_valid=rp, gene_keep_mask=gene_keep_mask, top_k=top_k
    )
    map_ab = dict(zip(rows_ab, idx_ab))
    map_ba = dict(zip(rows_ba, idx_ba))
    rows_common = [p for p in rows_ab if p in map_ba]
    m = len(rows_common)
    if m == 0:
        return {"n_perts": 0, "status": "no_valid_split_perturbations"}

    idx_a_c = [map_ab[p] for p in rows_common]
    idx_b_c = [map_ba[p] for p in rows_common]

    def _collect_dir(eval_obs_res: np.ndarray, idx_c: Sequence[np.ndarray]) -> Tuple[Tuple[np.ndarray, ...], Tuple[np.ndarray, ...]]:
        xs_o, ys_o, ha_o, hb_o = [], [], [], []
        for p, idx in zip(rows_common, idx_c):
            rv = tm[p, idx]
            xo, yo = _orth_pair(eval_obs_res[p, idx], rp[p, idx], rv)
            hao, hbo = _orth_pair(ha_res[p, idx], hb_res[p, idx], rv)
            xs_o.append(xo)
            ys_o.append(yo)
            ha_o.append(hao)
            hb_o.append(hbo)
        return _moments(xs_o, ys_o), _moments(ha_o, hb_o)

    mom_ab, hmom_ab = _collect_dir(hb_res, idx_a_c)
    mom_ba, hmom_ba = _collect_dir(ha_res, idx_b_c)

    def _eval_dir(
        pm: Tuple[np.ndarray, ...],
        hm: Tuple[np.ndarray, ...],
        draws: np.ndarray,
    ) -> Tuple[float, np.ndarray, float, np.ndarray, float, float]:
        nv, sx, sy, sxx, syy, sxy = pm
        hnv, hsx, hsy, hsxx, hsyy, hsxy = hm
        _, pt_r, _, _, _ = _moments_to_corr_and_slope(
            float(nv.sum()), float(sx.sum()), float(sy.sum()), float(sxx.sum()), float(syy.sum()), float(sxy.sum())
        )
        _, pt_rho, _, _, _ = _moments_to_corr_and_slope(
            float(hnv.sum()), float(hsx.sum()), float(hsy.sum()), float(hsxx.sum()), float(hsyy.sum()), float(hsxy.sum())
        )
        pt_ceil = single_half_ceiling(float(pt_rho))
        pt_cal = disattenuated_r(float(pt_r), pt_ceil)

        _, b_r, _, _, _ = _moments_to_corr_and_slope(
            nv[draws].sum(axis=1), sx[draws].sum(axis=1), sy[draws].sum(axis=1),
            sxx[draws].sum(axis=1), syy[draws].sum(axis=1), sxy[draws].sum(axis=1),
        )
        _, b_rho, _, _, _ = _moments_to_corr_and_slope(
            hnv[draws].sum(axis=1), hsx[draws].sum(axis=1), hsy[draws].sum(axis=1),
            hsxx[draws].sum(axis=1), hsyy[draws].sum(axis=1), hsxy[draws].sum(axis=1),
        )
        b_cal = np.array(
            [disattenuated_r(float(r), single_half_ceiling(float(rh))) for r, rh in zip(b_r, b_rho)],
            dtype=np.float64,
        )
        return float(pt_r), b_r, float(pt_cal), b_cal, float(pt_rho), float(pt_ceil)

    rng = np.random.default_rng(seed)
    draws = rng.integers(0, m, size=(n_boot, m))
    pt_raw_a, b_raw_a, pt_cal_a, b_cal_a, rho_a, ceil_a = _eval_dir(mom_ab, hmom_ab, draws)
    pt_raw_b, b_raw_b, pt_cal_b, b_cal_b, rho_b, ceil_b = _eval_dir(mom_ba, hmom_ba, draws)

    pt_raw_split = 0.5 * (pt_raw_a + pt_raw_b)
    b_raw_split = 0.5 * (b_raw_a + b_raw_b)
    pt_cal_split = 0.5 * (pt_cal_a + pt_cal_b)
    b_cal_split = 0.5 * (b_cal_a + b_cal_b)

    return {
        "n_perts": int(m),
        "r_orth_cal": summarize_boot_r11(pt_cal_split, b_cal_split, m, n_boot),
        "r_orth_raw": summarize_boot_r11(pt_raw_split, b_raw_split, m, n_boot),
        "r_orth_cal_a_to_b": pt_cal_a,
        "r_orth_cal_b_to_a": pt_cal_b,
        "r_orth_raw_a_to_b": pt_raw_a,
        "r_orth_raw_b_to_a": pt_raw_b,
        "rho_half_ab_orth_a": rho_a,
        "rho_half_ab_orth_b": rho_b,
        "ceiling_single_a": ceil_a,
        "ceiling_single_b": ceil_b,
    }


def compute_fct_perm100_split_r11(
    state: FctSourceState,
    target_ctrl_mean_libnorm: np.ndarray,
    query_target_names: Sequence[str],
    delta_transfer: np.ndarray,
    target_gene_idx: Sequence[int],
    delta_half_a: Optional[np.ndarray],
    delta_half_b: Optional[np.ndarray],
    half_valid: Optional[np.ndarray] = None,
    fc_bundle: Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]] = None,
    gene_keep_mask: Optional[np.ndarray] = None,
    n_perm: int = 100,
    top_k: int = 50,
) -> Dict[str, Any]:
    """
    100-permutation null (`k=0..n_perm-1`) under `SA-8` (`resp_top50_split`), computing
    symmetric average `0.5 * (r_cal_{A->B}[k] + r_cal_{B->A}[k])` with single-half ceilings
    `sqrt(rho_AB_orth(M^{(A)}))` and `sqrt(rho_AB_orth(M^{(B)}))`.
    """
    if delta_half_a is None or delta_half_b is None:
        return {"n_perm": int(n_perm), "mode": "resp_top50_split", "r_orth_cal": {}, "r_orth_raw": {}}

    mu_t = np.maximum(np.asarray(target_ctrl_mean_libnorm, dtype=np.float64), 0.0)
    G = mu_t.shape[0]
    if fc_bundle is None:
        fc_bundle = compute_source_transfer_and_fc_bar(state, query_target_names)
    pred_transfer, fc_bar, has_hits = fc_bundle

    ha = np.asarray(delta_half_a, dtype=np.float64)
    hb = np.asarray(delta_half_b, dtype=np.float64)
    tm = np.asarray(delta_transfer, dtype=np.float64)
    ha_res = ha - tm
    hb_res = hb - tm

    rows_ab, idx_ab = _prepare_split_direction_indices(
        ha, hb, ha, tm, target_gene_idx, half_valid, gene_keep_mask=gene_keep_mask, top_k=top_k
    )
    rows_ba, idx_ba = _prepare_split_direction_indices(
        hb, ha, hb, tm, target_gene_idx, half_valid, gene_keep_mask=gene_keep_mask, top_k=top_k
    )
    map_ab = dict(zip(rows_ab, idx_ab))
    map_ba = dict(zip(rows_ba, idx_ba))
    rows_common = [p for p in rows_ab if p in map_ba]
    if not rows_common:
        return {"n_perm": int(n_perm), "mode": "resp_top50_split", "r_orth_cal": {}, "r_orth_raw": {}}

    def _setup_dir(eval_obs_res: np.ndarray, idx_c: Sequence[np.ndarray]):
        xo_l, bv_l, hao_l, hbo_l = [], [], [], []
        for p, idx in zip(rows_common, idx_c):
            rv = tm[p, idx]
            nrm = float(np.linalg.norm(rv))
            bv = (rv / nrm) if nrm > EPS else np.zeros_like(rv)
            xv = eval_obs_res[p, idx]
            xo_l.append(xv - float(np.dot(xv, bv)) * bv)
            ha_o, hb_o = _orth_pair(ha_res[p, idx], hb_res[p, idx], rv)
            hao_l.append(ha_o)
            hbo_l.append(hb_o)
            bv_l.append(bv)
        x_pool = np.concatenate(xo_l)
        rho_ab = pearson_1d_strict(np.concatenate(hao_l), np.concatenate(hbo_l))
        ceil_s = single_half_ceiling(rho_ab)
        return x_pool, bv_l, ceil_s

    idx_a_c = [map_ab[p] for p in rows_common]
    idx_b_c = [map_ba[p] for p in rows_common]
    x_a_pool, bv_a, ceil_a = _setup_dir(hb_res, idx_a_c)
    x_b_pool, bv_b, ceil_b = _setup_dir(ha_res, idx_b_c)

    raw_vals = np.empty(n_perm, dtype=np.float64)
    cal_vals = np.empty(n_perm, dtype=np.float64)

    for k in range(n_perm):
        perm = np.random.default_rng(k).permutation(G)
        mu_p = mu_t[perm]
        yp_a_l, yp_b_l = [], []
        for p, ia, ba_v, ib, bb_v in zip(rows_common, idx_a_c, bv_a, idx_b_c, bv_b):
            hit_p = bool(has_hits[p])
            for idx_dir, bv_dir, out_list in ((ia, ba_v, yp_a_l), (ib, bb_v, yp_b_l)):
                if hit_p:
                    m_sub = mu_p[idx_dir]
                    pred_v = np.log1p(m_sub * fc_bar[p, idx_dir]) - np.log1p(m_sub)
                else:
                    pred_v = pred_transfer[p, idx_dir]
                yv = pred_v - tm[p, idx_dir]
                out_list.append(yv - float(np.dot(yv, bv_dir)) * bv_dir)

        r_raw_a = pearson_1d_strict(x_a_pool, np.concatenate(yp_a_l))
        r_raw_b = pearson_1d_strict(x_b_pool, np.concatenate(yp_b_l))
        r_cal_a = disattenuated_r(r_raw_a, ceil_a)
        r_cal_b = disattenuated_r(r_raw_b, ceil_b)
        raw_vals[k] = 0.5 * (r_raw_a + r_raw_b)
        cal_vals[k] = 0.5 * (r_cal_a + r_cal_b)

    def _pct(arr: np.ndarray) -> Dict[str, float]:
        fin = arr[np.isfinite(arr)]
        if fin.size == 0:
            return {k: float("nan") for k in ("p2_5", "median", "p97_5", "p98_75", "max", "mean")}
        p25, p50, p975, p9875 = np.percentile(fin, [2.5, 50.0, 97.5, 98.75])
        return {
            "p2_5": float(p25),
            "median": float(p50),
            "p97_5": float(p975),
            "p98_75": float(p9875),
            "max": float(np.max(fin)),
            "mean": float(np.mean(fin)),
        }

    return {
        "n_perm": int(n_perm),
        "mode": "resp_top50_split",
        "r_orth_cal": _pct(cal_vals),
        "r_orth_raw": _pct(raw_vals),
        "ceiling_single_a": ceil_a,
        "ceiling_single_b": ceil_b,
        "perm0_r_orth_cal": float(cal_vals[0]),
        "perm0_r_orth_raw": float(raw_vals[0]),
    }
