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
    """Summary of a perturbation-level bootstrap confidence interval (95% and 97.5% two-sided)."""

    estimate: float
    ci95_lo: float
    ci95_hi: float
    boot_se: float
    n_items: int
    n_boot: int
    ci975_lo: float = float("nan")
    ci975_hi: float = float("nan")


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


def deming_slope(x: np.ndarray, y: np.ndarray) -> float:
    """Compute orthogonal (equal error-variance ratio delta=1) Deming regression slope of y on x."""
    x_arr = np.asarray(x, dtype=np.float64)
    y_arr = np.asarray(y, dtype=np.float64)
    if x_arr.size < 2 or y_arr.size != x_arr.size:
        return 0.0
    xc = x_arr - np.mean(x_arr)
    yc = y_arr - np.mean(y_arr)
    sxx = float(np.dot(xc, xc))
    syy = float(np.dot(yc, yc))
    sxy = float(np.dot(xc, yc))
    if abs(sxy) <= EPS:
        return 0.0
    diff = syy - sxx
    return float((diff + np.sqrt(diff * diff + 4.0 * sxy * sxy)) / (2.0 * sxy))


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


def compute_disc_ownmasked_rows(
    true_delta: np.ndarray,
    pred_delta: np.ndarray,
    target_col_idx: Sequence[int | None],
    center_across_genes: bool = True,
) -> np.ndarray:
    """Compute tie-aware perturbation discrimination zeroing only each row's own target gene (A-T2 convention)."""
    t_mat = np.asarray(true_delta, dtype=np.float64).copy()
    p_mat = np.asarray(pred_delta, dtype=np.float64).copy()
    n_perts, n_genes = t_mat.shape
    if n_perts <= 1:
        return np.full(n_perts, 0.5, dtype=np.float64)

    for i, t_idx in enumerate(target_col_idx):
        if t_idx is not None and 0 <= int(t_idx) < n_genes:
            t_mat[i, int(t_idx)] = 0.0
            p_mat[i, int(t_idx)] = 0.0

    if center_across_genes:
        t_mat = t_mat - np.mean(t_mat, axis=1, keepdims=True)
        p_mat = p_mat - np.mean(p_mat, axis=1, keepdims=True)

    t_norms = np.linalg.norm(t_mat, axis=1, keepdims=True)
    p_norms = np.linalg.norm(p_mat, axis=1, keepdims=True)
    t_unit = t_mat / np.maximum(t_norms, EPS)
    p_unit = p_mat / np.maximum(p_norms, EPS)
    sim = p_unit @ t_unit.T

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
    lo975, hi975 = np.percentile(boots, [1.25, 98.75])
    return BootstrapCI(
        estimate=float(np.mean(diff)),
        ci95_lo=float(lo),
        ci95_hi=float(hi),
        boot_se=float(np.std(boots, ddof=1)),
        n_items=n,
        n_boot=n_boot,
        ci975_lo=float(lo975),
        ci975_hi=float(hi975),
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
    lo975, hi975 = np.percentile(finite_boots, [1.25, 98.75])
    return BootstrapCI(
        estimate=float(estimate),
        ci95_lo=float(lo),
        ci95_hi=float(hi),
        boot_se=float(np.std(finite_boots, ddof=1)) if finite_boots.size > 1 else 0.0,
        n_items=n_items,
        n_boot=n_boot,
        ci975_lo=float(lo975),
        ci975_hi=float(hi975),
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


def _pick_bound(row: dict, prefer_keys: Sequence[str]) -> float:
    for k in prefer_keys:
        if k in row:
            return float(row[k])
    return float("nan")


def adjudicate_team_a(
    line_results: Sequence[dict],
    min_lines_pass: int = 4,
    max_lines_inferior: int = 1,
    slope_threshold: float = 0.13,
    confirmatory_lines: Sequence[str] | None = None,
    use_975_ci: bool = True,
) -> dict:
    """Apply pre-registered fixed adjudication rules for Group A H-A1..H-A4 (R2/R3/R5/R10/R13)."""
    eval_rows = (
        [r for r in line_results if r["line"] in set(confirmatory_lines)]
        if confirmatory_lines is not None
        else list(line_results)
    )
    n_lines = len(eval_rows)
    ci_tag = "975" if use_975_ci else "95"

    def _slope_lo(r: dict) -> float:
        return _pick_bound(r, (f"primary_slope_ci{ci_tag}_lo", "primary_slope_ci95_lo"))

    def _orth_cal_lo(r: dict) -> float:
        return _pick_bound(
            r,
            (
                f"r_orth_cal_ci{ci_tag}_lo",
                "r_orth_cal_ci95_lo",
                f"noise_cal_corr_ci{ci_tag}_lo",
                "noise_cal_corr_ci95_lo",
            ),
        )

    def _mis_slope_hi(r: dict) -> float:
        return _pick_bound(
            r,
            (
                f"max_mismatch_ci{ci_tag}_hi",
                f"mismatch_slope_ci{ci_tag}_hi",
                "max_mismatch_ci95_hi",
                "mismatch_slope_ci95_hi",
            ),
        )

    def _mis_orth_hi(r: dict) -> float:
        return _pick_bound(
            r,
            (
                f"max_mismatch_r_orth_cal_ci{ci_tag}_hi",
                "max_mismatch_r_orth_cal_ci95_hi",
            ),
        )

    ha1_pass_lines = [r["line"] for r in eval_rows if np.isfinite(_slope_lo(r)) and _slope_lo(r) > slope_threshold]
    ha2_pass_lines = [r["line"] for r in eval_rows if np.isfinite(_orth_cal_lo(r)) and _orth_cal_lo(r) > 0.0]
    has_r17_orth_a = bool(eval_rows) and any(
        k in eval_rows[0] for k in (f"r_orth_cal_ci{ci_tag}_lo", "r_orth_cal_ci95_lo")
    ) and any(
        k in eval_rows[0] for k in (f"max_mismatch_r_orth_cal_ci{ci_tag}_hi", "max_mismatch_r_orth_cal_ci95_hi")
    )
    ha3_pass_lines = []
    for r in eval_rows:
        slope_ok = np.isfinite(_mis_slope_hi(r)) and np.isfinite(_slope_lo(r)) and (_mis_slope_hi(r) < _slope_lo(r))
        if has_r17_orth_a:
            orth_ok = np.isfinite(_mis_orth_hi(r)) and np.isfinite(_orth_cal_lo(r)) and (_mis_orth_hi(r) < _orth_cal_lo(r))
            if slope_ok and orth_ok:
                ha3_pass_lines.append(r["line"])
        elif slope_ok:
            ha3_pass_lines.append(r["line"])

    ha4_key_groups = [
        ("ctrl_disc_ownmasked", (f"ctrl_disc_ownmasked_diff_ci{ci_tag}_hi", "ctrl_disc_ownmasked_diff_ci95_hi")),
        ("cent_disc_raw", (f"cent_disc_raw_diff_ci{ci_tag}_hi", "cent_disc_raw_diff_ci95_hi")),
        ("cent_disc_ownmasked", (f"cent_disc_ownmasked_diff_ci{ci_tag}_hi", "cent_disc_ownmasked_diff_ci95_hi")),
        ("cent_de50_r", (f"cent_de50_r_diff_ci{ci_tag}_hi", "cent_de50_r_diff_ci95_hi")),
    ]
    ha4_disc_inf_lines = [
        r["line"]
        for r in eval_rows
        if np.isfinite(_pick_bound(r, (f"disc_raw_diff_ci{ci_tag}_hi", "disc_raw_diff_ci95_hi")))
        and _pick_bound(r, (f"disc_raw_diff_ci{ci_tag}_hi", "disc_raw_diff_ci95_hi")) < 0.0
    ]
    ha4_de50_inf_lines = [
        r["line"]
        for r in eval_rows
        if np.isfinite(_pick_bound(r, (f"de50_r_diff_ci{ci_tag}_hi", "de50_r_diff_ci95_hi")))
        and _pick_bound(r, (f"de50_r_diff_ci{ci_tag}_hi", "de50_r_diff_ci95_hi")) < 0.0
    ]

    extra_inf: dict[str, list[str]] = {}
    for label, key_prefs in ha4_key_groups:
        if eval_rows and any(k in eval_rows[0] for k in key_prefs):
            extra_inf[label] = [
                r["line"]
                for r in eval_rows
                if np.isfinite(_pick_bound(r, key_prefs)) and _pick_bound(r, key_prefs) < 0.0
            ]

    ha1_ok = len(ha1_pass_lines) >= min_lines_pass
    ha2_ok = len(ha2_pass_lines) >= min_lines_pass
    ha3_ok = len(ha3_pass_lines) >= min_lines_pass
    ha4_ok = (
        (len(ha4_disc_inf_lines) <= max_lines_inferior)
        and (len(ha4_de50_inf_lines) <= max_lines_inferior)
        and all(len(v) <= max_lines_inferior for v in extra_inf.values())
    )

    overall = "CONFIRMED" if (ha1_ok and ha2_ok and ha3_ok and ha4_ok) else "未确认 (NOT CONFIRMED)"
    ha4_dict = {
        "pass": ha4_ok,
        "disc_raw_inferior_count": len(ha4_disc_inf_lines),
        "disc_raw_inferior_lines": ha4_disc_inf_lines,
        "de50_r_inferior_count": len(ha4_de50_inf_lines),
        "de50_r_inferior_lines": ha4_de50_inf_lines,
    }
    for label, lines_inf in extra_inf.items():
        ha4_dict[f"{label}_inferior_count"] = len(lines_inf)
        ha4_dict[f"{label}_inferior_lines"] = lines_inf
    return {
        "ci_level": 0.975 if use_975_ci else 0.95,
        "n_lines": n_lines,
        "confirmatory_lines": [r["line"] for r in eval_rows],
        "H-A1": {"pass": ha1_ok, "count": len(ha1_pass_lines), "lines": ha1_pass_lines},
        "H-A2": {"pass": ha2_ok, "count": len(ha2_pass_lines), "lines": ha2_pass_lines},
        "H-A3": {"pass": ha3_ok, "count": len(ha3_pass_lines), "lines": ha3_pass_lines},
        "H-A4": ha4_dict,
        "overall_verdict": overall,
    }


def adjudicate_team_b(
    line_results: Sequence[dict],
    min_lines_pass: int = 3,
    max_lines_inferior: int = 0,
    slope_threshold: float = 0.13,
    confirmatory_lines: Sequence[str] | None = None,
    use_975_ci: bool = True,
) -> dict:
    """Apply pre-registered Group B (`Transfer_Mean_VarScaled`) adjudication (R12 & R13)."""
    eval_rows = (
        [r for r in line_results if r["line"] in set(confirmatory_lines)]
        if confirmatory_lines is not None
        else list(line_results)
    )
    n_lines = len(eval_rows)
    ci_tag = "975" if use_975_ci else "95"

    def _slope_lo(r: dict) -> float:
        return _pick_bound(r, (f"tmvs_slope_ci{ci_tag}_lo", "tmvs_slope_ci95_lo"))

    has_r17_orth_b = bool(eval_rows) and any(
        k in eval_rows[0] for k in (f"tmvs_r_orth_cal_ci{ci_tag}_lo", "tmvs_r_orth_cal_ci95_lo")
    )

    def _orth_cal_lo(r: dict) -> float:
        return _pick_bound(
            r,
            (
                f"tmvs_r_orth_cal_ci{ci_tag}_lo",
                "tmvs_r_orth_cal_ci95_lo",
                f"tmvs_cal_r_ci{ci_tag}_lo",
                "tmvs_cal_r_ci95_lo",
            ),
        )

    def _ctrl_orth_hi(r: dict) -> float:
        return _pick_bound(
            r,
            (
                f"tmvs_max_ctrl_r_orth_cal_ci{ci_tag}_hi",
                "tmvs_max_ctrl_r_orth_cal_ci95_hi",
                f"tmvs_max_ctrl_cal_r_ci{ci_tag}_hi",
                "tmvs_max_ctrl_cal_r_ci95_hi",
            ),
        )

    def _ctrl_slope_hi(r: dict) -> float:
        return _pick_bound(r, (f"tmvs_max_ctrl_slope_ci{ci_tag}_hi", "tmvs_max_ctrl_slope_ci95_hi"))

    hb1_pass_lines = [r["line"] for r in eval_rows if np.isfinite(_slope_lo(r)) and _slope_lo(r) > slope_threshold]
    hb2_pass_lines = [r["line"] for r in eval_rows if np.isfinite(_orth_cal_lo(r)) and _orth_cal_lo(r) > 0.0]
    if has_r17_orth_b:
        # R17 + R18: compare r_orth_cal CI lower bound > max(all 5 mismatched r_orth_cal CI_hi, 100-perm 97.5th pct)
        hb3_pass_lines = [
            r["line"]
            for r in eval_rows
            if np.isfinite(_orth_cal_lo(r))
            and np.isfinite(_ctrl_orth_hi(r))
            and (_orth_cal_lo(r) > _ctrl_orth_hi(r))
        ]
    else:
        hb3_pass_lines = [
            r["line"]
            for r in eval_rows
            if np.isfinite(_slope_lo(r))
            and np.isfinite(_ctrl_slope_hi(r))
            and np.isfinite(_orth_cal_lo(r))
            and np.isfinite(_ctrl_orth_hi(r))
            and (_slope_lo(r) > _ctrl_slope_hi(r))
            and (_orth_cal_lo(r) > _ctrl_orth_hi(r))
        ]

    hb4_key_groups = [
        ("ctrl_disc_panel", (f"tmvs_disc_raw_diff_ci{ci_tag}_hi", "tmvs_disc_raw_diff_ci95_hi")),
        ("ctrl_disc_ownmasked", (f"tmvs_ctrl_disc_ownmasked_diff_ci{ci_tag}_hi", "tmvs_ctrl_disc_ownmasked_diff_ci95_hi")),
        ("ctrl_de50_r", (f"tmvs_de50_r_diff_ci{ci_tag}_hi", "tmvs_de50_r_diff_ci95_hi")),
        ("cent_disc_panel", (f"tmvs_cent_disc_raw_diff_ci{ci_tag}_hi", "tmvs_cent_disc_raw_diff_ci95_hi")),
        ("cent_disc_ownmasked", (f"tmvs_cent_disc_ownmasked_diff_ci{ci_tag}_hi", "tmvs_cent_disc_ownmasked_diff_ci95_hi")),
        ("cent_de50_r", (f"tmvs_cent_de50_r_diff_ci{ci_tag}_hi", "tmvs_cent_de50_r_diff_ci95_hi")),
    ]
    hb4_inf: dict[str, list[str]] = {}
    for label, key_prefs in hb4_key_groups:
        if eval_rows and any(k in eval_rows[0] for k in key_prefs):
            hb4_inf[label] = [
                r["line"]
                for r in eval_rows
                if np.isfinite(_pick_bound(r, key_prefs)) and _pick_bound(r, key_prefs) < 0.0
            ]

    hb1_ok = len(hb1_pass_lines) >= min_lines_pass
    hb2_ok = len(hb2_pass_lines) >= min_lines_pass
    hb3_ok = len(hb3_pass_lines) >= min_lines_pass
    hb4_ok = bool(hb4_inf) and all(len(v) <= max_lines_inferior for v in hb4_inf.values())

    in_vcc_top10_lines = [
        r["line"]
        for r in eval_rows
        if "tmvs_slope" in r and np.isfinite(r["tmvs_slope"]) and 0.02 <= float(r["tmvs_slope"]) <= 0.13
    ]

    overall = "CONFIRMED" if (hb1_ok and hb2_ok and hb3_ok and hb4_ok) else "未确认 (NOT CONFIRMED)"
    hb4_dict: dict[str, object] = {"pass": hb4_ok}
    for label, lines_inf in hb4_inf.items():
        hb4_dict[f"{label}_inferior_count"] = len(lines_inf)
        hb4_dict[f"{label}_inferior_lines"] = lines_inf

    return {
        "ci_level": 0.975 if use_975_ci else 0.95,
        "n_lines": n_lines,
        "confirmatory_lines": [r["line"] for r in eval_rows],
        "H-B1": {"pass": hb1_ok, "count": len(hb1_pass_lines), "lines": hb1_pass_lines},
        "H-B2": {"pass": hb2_ok, "count": len(hb2_pass_lines), "lines": hb2_pass_lines},
        "H-B3": {"pass": hb3_ok, "count": len(hb3_pass_lines), "lines": hb3_pass_lines},
        "H-B4": hb4_dict,
        "vcc2026_top10_range_check": {
            "range": [0.02, 0.13],
            "count": len(in_vcc_top10_lines),
            "lines": in_vcc_top10_lines,
            "r6_caveat": "External Arc VCC 2026 Wrap-Up slope estimator is unpublished in cell-eval2; OLS slopes depend on prediction SD ratio.",
        },
        "overall_verdict": overall,
    }
