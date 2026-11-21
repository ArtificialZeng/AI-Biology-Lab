#!/usr/bin/env python3
"""Independent statistical and metric audit primitives for Team A (A-T3).

This module is strictly standalone (zero imports from A-T2's codebase or earlier
project scripts) so that A-T3 can independently audit and recompute all main
evaluation metrics, noise ceilings, bootstrap confidence intervals, and
pre-registered hypothesis verdicts.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable, Sequence

import numpy as np


EPS = 1e-12


@dataclass(frozen=True)
class BootstrapCI:
    """Summary of a perturbation-level bootstrap confidence interval."""

    estimate: float
    ci95_lo: float
    ci95_hi: float
    boot_se: float
    n_items: int
    n_boot: int


def pearson_1d(x: np.ndarray, y: np.ndarray) -> float:
    """Compute finite-safe 1D Pearson correlation between two equal-length vectors."""
    x_arr = np.asarray(x, dtype=np.float64)
    y_arr = np.asarray(y, dtype=np.float64)
    if x_arr.size < 2 or y_arr.size != x_arr.size:
        return 0.0
    x_c = x_arr - np.mean(x_arr)
    y_c = y_arr - np.mean(y_arr)
    denom = float(np.linalg.norm(x_c) * np.linalg.norm(y_c))
    if denom <= EPS:
        return 0.0
    val = float(np.dot(x_c, y_c) / denom)
    return float(np.clip(val, -1.0, 1.0)) if np.isfinite(val) else 0.0


def ols_slope(x: np.ndarray, y: np.ndarray, fit_intercept: bool = True) -> float:
    """Compute OLS slope regressing y on x, with or without an intercept."""
    x_arr = np.asarray(x, dtype=np.float64)
    y_arr = np.asarray(y, dtype=np.float64)
    if x_arr.size == 0 or y_arr.size != x_arr.size:
        return 0.0
    if fit_intercept:
        x_c = x_arr - np.mean(x_arr)
        y_c = y_arr - np.mean(y_arr)
        denom = float(np.dot(x_c, x_c))
        return float(np.dot(x_c, y_c) / denom) if denom > EPS else 0.0
    denom = float(np.dot(x_arr, x_arr))
    return float(np.dot(x_arr, y_arr) / denom) if denom > EPS else 0.0


def row_mask_array(
    n_genes: int,
    exclude_indices: Iterable[int] | None = None,
) -> np.ndarray:
    """Build a boolean mask over genes excluding specified column indices."""
    mask = np.ones(n_genes, dtype=bool)
    if exclude_indices is not None:
        for idx in exclude_indices:
            if idx is not None and 0 <= int(idx) < n_genes:
                mask[int(idx)] = False
    return mask


def compute_de50_r_rows(
    true_delta: np.ndarray,
    pred_delta: np.ndarray,
    target_col_idx: Sequence[int | None] | None = None,
    top_k: int = 50,
) -> np.ndarray:
    """Compute per-perturbation top-k |true_delta| Pearson correlation, masking own target."""
    t_mat = np.asarray(true_delta, dtype=np.float64)
    p_mat = np.asarray(pred_delta, dtype=np.float64)
    n_perts, n_genes = t_mat.shape
    out = np.empty(n_perts, dtype=np.float64)
    for i in range(n_perts):
        t_idx = target_col_idx[i] if target_col_idx is not None else None
        mask = row_mask_array(n_genes, [t_idx] if t_idx is not None else None)
        dt = t_mat[i, mask]
        dp = p_mat[i, mask]
        k = min(top_k, dt.size)
        if k < 2:
            out[i] = 0.0
            continue
        top_idx = np.argpartition(np.abs(dt), -k)[-k:]
        out[i] = pearson_1d(dt[top_idx], dp[top_idx])
    return out


def compute_disc_raw_rows(
    true_delta: np.ndarray,
    pred_delta: np.ndarray,
    masked_cols: Sequence[int] | None = None,
    center_across_genes: bool = True,
) -> np.ndarray:
    """Compute tie-aware mid-rank perturbation discrimination percentile per row.

    Matches the exact tie-aware mid-rank formula:
      disc(i) = 1 - [#{j: c_ij > c_ii} + 0.5 * (#{j != i: c_ij == c_ii})] / (M - 1)
    where c_ij is Pearson correlation (if center_across_genes=True) or cosine
    similarity (if center_across_genes=False) between pred_i and true_j.
    """
    t_mat = np.asarray(true_delta, dtype=np.float64)
    p_mat = np.asarray(pred_delta, dtype=np.float64)
    n_perts, n_genes = t_mat.shape
    if n_perts <= 1:
        return np.full(n_perts, 0.5, dtype=np.float64)

    mask = row_mask_array(n_genes, masked_cols)
    t_sub = t_mat[:, mask]
    p_sub = p_mat[:, mask]

    if center_across_genes:
        t_sub = t_sub - np.mean(t_sub, axis=1, keepdims=True)
        p_sub = p_sub - np.mean(p_sub, axis=1, keepdims=True)

    t_norms = np.linalg.norm(t_sub, axis=1, keepdims=True)
    p_norms = np.linalg.norm(p_sub, axis=1, keepdims=True)
    t_unit = t_sub / np.maximum(t_norms, EPS)
    p_unit = p_sub / np.maximum(p_norms, EPS)

    sim = p_unit @ t_unit.T
    # Zero-norm rows carry no directional information -> exact random tie (0.5)
    p_zero = p_norms[:, 0] <= 1e-10
    t_zero = t_norms[:, 0] <= 1e-10
    if np.any(p_zero):
        sim[p_zero, :] = 0.0
    if np.any(t_zero):
        sim[:, t_zero] = 0.0

    diag = np.diag(sim)[:, None]
    strictly_better = np.sum(sim > diag + 1e-12, axis=1)
    tied_others = np.sum(np.abs(sim - diag) <= 1e-12, axis=1) - 1
    eff_better = strictly_better + 0.5 * np.maximum(tied_others, 0)
    out = 1.0 - eff_better / float(n_perts - 1)
    out[p_zero] = 0.5
    return out


def spearman_brown_reliability(r_half_ab: float) -> float:
    """Convert split-half correlation r_AB into full-sample reliability rho_full."""
    if not np.isfinite(r_half_ab) or r_half_ab <= 0.0:
        return 0.0
    return float(np.clip((2.0 * r_half_ab) / (1.0 + r_half_ab), 0.0, 1.0))


def noise_ceiling_from_splithalf(r_half_ab: float) -> float:
    """Compute correlation noise ceiling sqrt(rho_full) from split-half correlation."""
    rho = spearman_brown_reliability(r_half_ab)
    return float(np.sqrt(rho)) if rho > 0.0 else 0.0


def disattenuated_correlation(r_obs: float, r_ceiling: float) -> float:
    """Divide observed correlation by noise ceiling; returns NaN if ceiling <= 0."""
    if not np.isfinite(r_obs) or not np.isfinite(r_ceiling) or r_ceiling <= EPS:
        return float("nan")
    return float(r_obs / r_ceiling)


def bootstrap_paired_mean_diff(
    vals_a: np.ndarray,
    vals_b: np.ndarray,
    n_boot: int = 4000,
    seed: int = 20261121,
) -> BootstrapCI:
    """Bootstrap 95% percentile CI for mean(vals_a - vals_b) across rows."""
    a = np.asarray(vals_a, dtype=np.float64)
    b = np.asarray(vals_b, dtype=np.float64)
    if a.shape != b.shape or a.ndim != 1 or a.size == 0:
        raise ValueError("vals_a and vals_b must be non-empty 1D arrays of equal shape.")
    diff = a - b
    n = diff.size
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(n_boot, n))
    boots = diff[idx].mean(axis=1)
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return BootstrapCI(
        estimate=float(np.mean(diff)),
        ci95_lo=float(lo),
        ci95_hi=float(hi),
        boot_se=float(np.std(boots, ddof=1)),
        n_items=n,
        n_boot=n_boot,
    )


def summarize_bootstrap_array(
    estimate: float,
    boots: np.ndarray,
    n_items: int,
    min_finite_frac: float = 0.95,
) -> BootstrapCI:
    """Summarize 1D bootstrap draw array, rejecting draws if > (1 - min_finite_frac) are non-finite."""
    arr = np.asarray(boots, dtype=np.float64)
    n_boot = int(arr.size)
    finite_boots = arr[np.isfinite(arr)]
    if n_boot == 0 or finite_boots.size < max(1, int(np.ceil(min_finite_frac * n_boot))):
        return BootstrapCI(estimate, float("nan"), float("nan"), float("nan"), n_items, n_boot)
    lo, hi = np.percentile(finite_boots, [2.5, 97.5])
    return BootstrapCI(
        estimate=float(estimate),
        ci95_lo=float(lo),
        ci95_hi=float(hi),
        boot_se=float(np.std(finite_boots, ddof=1)) if finite_boots.size > 1 else 0.0,
        n_items=n_items,
        n_boot=n_boot,
    )


def bootstrap_perturbation_stat(
    n_items: int,
    stat_fn: Callable[[np.ndarray], float],
    n_boot: int = 4000,
    seed: int = 20261121,
) -> BootstrapCI:
    """Bootstrap 95% percentile CI by resampling perturbation indices with replacement."""
    if n_items <= 0:
        raise ValueError("n_items must be positive.")
    full_idx = np.arange(n_items, dtype=np.int64)
    est = float(stat_fn(full_idx))
    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot, dtype=np.float64)
    for b in range(n_boot):
        boot_idx = rng.integers(0, n_items, size=n_items)
        boots[b] = float(stat_fn(boot_idx))
    return summarize_bootstrap_array(est, boots, n_items=n_items)


def ci_strictly_below(ci_low: BootstrapCI, ci_high: BootstrapCI) -> bool:
    """Check non-overlapping 95% CI ordering: ci_low.ci95_hi < ci_high.ci95_lo."""
    if not (np.isfinite(ci_low.ci95_hi) and np.isfinite(ci_high.ci95_lo)):
        return False
    return bool(ci_low.ci95_hi < ci_high.ci95_lo)


def adjudicate_team_a(
    line_results: Sequence[dict],
    min_lines_pass: int = 4,
    max_lines_inferior: int = 1,
    slope_threshold: float = 0.13,
) -> dict:
    """Apply pre-registered fixed adjudication rules for H-A1 through H-A4.

    Each element of `line_results` must contain:
      - line: str
      - primary_slope_ci95_lo: float
      - noise_cal_corr_ci95_lo: float
      - mismatch_slope_ci95_hi: float
      - disc_raw_diff_ci95_hi: float (hybrid - transferred_mean)
      - de50_r_diff_ci95_hi: float   (hybrid - transferred_mean)
    """
    n_lines = len(line_results)
    ha1_pass_lines = [
        r["line"] for r in line_results if np.isfinite(r["primary_slope_ci95_lo"]) and r["primary_slope_ci95_lo"] > slope_threshold
    ]
    ha2_pass_lines = [
        r["line"] for r in line_results if np.isfinite(r["noise_cal_corr_ci95_lo"]) and r["noise_cal_corr_ci95_lo"] > 0.0
    ]
    ha3_pass_lines = [
        r["line"]
        for r in line_results
        if np.isfinite(r["mismatch_slope_ci95_hi"])
        and np.isfinite(r["primary_slope_ci95_lo"])
        and r["mismatch_slope_ci95_hi"] < r["primary_slope_ci95_lo"]
    ]
    ha4_disc_inf_lines = [
        r["line"] for r in line_results if np.isfinite(r["disc_raw_diff_ci95_hi"]) and r["disc_raw_diff_ci95_hi"] < 0.0
    ]
    ha4_de50_inf_lines = [
        r["line"] for r in line_results if np.isfinite(r["de50_r_diff_ci95_hi"]) and r["de50_r_diff_ci95_hi"] < 0.0
    ]
    has_cent = bool(line_results) and ("cent_disc_raw_diff_ci95_hi" in line_results[0]) and ("cent_de50_r_diff_ci95_hi" in line_results[0])
    ha4_cent_disc_inf = (
        [r["line"] for r in line_results if np.isfinite(r["cent_disc_raw_diff_ci95_hi"]) and r["cent_disc_raw_diff_ci95_hi"] < 0.0]
        if has_cent
        else []
    )
    ha4_cent_de50_inf = (
        [r["line"] for r in line_results if np.isfinite(r["cent_de50_r_diff_ci95_hi"]) and r["cent_de50_r_diff_ci95_hi"] < 0.0]
        if has_cent
        else []
    )

    ha1_ok = len(ha1_pass_lines) >= min_lines_pass
    ha2_ok = len(ha2_pass_lines) >= min_lines_pass
    ha3_ok = len(ha3_pass_lines) >= min_lines_pass
    ha4_ok = (
        (len(ha4_disc_inf_lines) <= max_lines_inferior)
        and (len(ha4_de50_inf_lines) <= max_lines_inferior)
        and ((not has_cent) or (len(ha4_cent_disc_inf) <= max_lines_inferior and len(ha4_cent_de50_inf) <= max_lines_inferior))
    )

    overall = "CONFIRMED" if (ha1_ok and ha2_ok and ha3_ok and ha4_ok) else "未确认 (NOT CONFIRMED)"
    ha4_dict = {
        "pass": ha4_ok,
        "disc_raw_inferior_count": len(ha4_disc_inf_lines),
        "disc_raw_inferior_lines": ha4_disc_inf_lines,
        "de50_r_inferior_count": len(ha4_de50_inf_lines),
        "de50_r_inferior_lines": ha4_de50_inf_lines,
    }
    if has_cent:
        ha4_dict["cent_disc_raw_inferior_count"] = len(ha4_cent_disc_inf)
        ha4_dict["cent_disc_raw_inferior_lines"] = ha4_cent_disc_inf
        ha4_dict["cent_de50_r_inferior_count"] = len(ha4_cent_de50_inf)
        ha4_dict["cent_de50_r_inferior_lines"] = ha4_cent_de50_inf
    return {
        "n_lines": n_lines,
        "H-A1": {"pass": ha1_ok, "count": len(ha1_pass_lines), "lines": ha1_pass_lines},
        "H-A2": {"pass": ha2_ok, "count": len(ha2_pass_lines), "lines": ha2_pass_lines},
        "H-A3": {"pass": ha3_ok, "count": len(ha3_pass_lines), "lines": ha3_pass_lines},
        "H-A4": ha4_dict,
        "overall_verdict": overall,
    }
