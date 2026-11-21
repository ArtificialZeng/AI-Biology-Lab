#!/usr/bin/env python3
"""
Evaluation metrics, Miller DRF calibration, and row-level paired bootstrap CIs.

Includes:
  1. BioXArena / Miller DRF metrics (`disc_raw`, `de50_r`, `de50_sign`, `delta_r`, `mse`,
     `v5`, `disc_v5`) in both centroid-referenced (matching `src/round8_r86_miller_drf.py`)
     and control-referenced conventions.
  2. Miller Dynamic Range Fraction `drf_score(m_pos, m_neg, lower_is_better)`.
  3. Row-level paired bootstrap 95% CI (`paired_row_bootstrap_ci` & `stratified_row_bootstrap_ci`)
     matching `src/compute_row_bootstrap_ci.py` (4,000 resamples, RandomState(12345)).
  4. Line-specific residual decomposition metrics (placeholder interface aligned with
     Arc VCC 2026 Wrap-Up; exact official formula to be locked from A-T5's
     `audit/LIT/LINE_SPECIFIC_METRIC.md`).
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


def pearson_1d(a: np.ndarray, b: np.ndarray) -> float:
    """Safe 1D Pearson correlation clamped to [-1, 1]; returns 0.0 on zero variance."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.size < 2:
        return 0.0
    ac = a - a.mean()
    bc = b - b.mean()
    denom = float(np.linalg.norm(ac) * np.linalg.norm(bc))
    if denom <= 1e-12 or not np.isfinite(denom):
        return 0.0
    v = float(np.dot(ac, bc) / denom)
    return float(np.clip(v, -1.0, 1.0)) if np.isfinite(v) else 0.0


def tie_aware_disc_rows(
    P_mat: np.ndarray,
    T_mat: np.ndarray,
    zero_tol: float = 1e-10,
) -> np.ndarray:
    """
    Per-row perturbation discrimination rank percentile (tie-aware mid-rank),
    matching `eval_subset_metrics` in `src/round8_r86_miller_drf.py` and
    `pds_cosine` in `cell-eval2`:
      disc_i = 1.0 - (strictly_better + 0.5 * tied_others) / (N - 1)
    Rows where prediction vector has zero variance/norm get 0.5 if constant/degenerate,
    or exact tie-aware rank when all rows match.
    """
    n = T_mat.shape[0]
    if n <= 1:
        return np.full(n, 0.5, dtype=np.float64)
    A0 = T_mat - T_mat.mean(axis=1, keepdims=True)
    B0 = P_mat - P_mat.mean(axis=1, keepdims=True)
    An = A0 / np.maximum(np.linalg.norm(A0, axis=1, keepdims=True), 1e-12)
    Bn = B0 / np.maximum(np.linalg.norm(B0, axis=1, keepdims=True), 1e-12)
    C = Bn @ An.T  # C[i, j] = corr(pred_i, true_j)
    diag = np.diag(C)[:, None]
    strictly_better = np.sum(C > diag + 1e-12, axis=1)
    tied_others = np.sum(np.abs(C - diag) <= 1e-12, axis=1) - 1
    eff_better = strictly_better + 0.5 * np.maximum(tied_others, 0)
    pcts = 1.0 - eff_better / max(n - 1, 1)
    # Zero-norm predictions have C[i, :] == 0 everywhere -> all tied -> pcts == 0.5 automatically!
    # Verify explicitly for safety when norm(B0) < zero_tol:
    b_zero = np.linalg.norm(B0, axis=1) < zero_tol
    pcts = np.where(b_zero, 0.5, pcts)
    return pcts.astype(np.float64)


def compute_row_metrics(
    pred_deltas: np.ndarray,
    true_deltas: np.ndarray,
    target_indices: Sequence[int],
    centroid: Optional[np.ndarray] = None,
    shared_axis: Optional[np.ndarray] = None,
    top_k_de: int = 50,
    mask_target_in_disc: bool = True,
) -> Dict[str, np.ndarray]:
    """
    Compute per-perturbation row-level metric arrays for one cell line.

    Parameters
    ----------
    pred_deltas : (N, G) predicted delta expression relative to control.
    true_deltas : (N, G) true delta expression relative to control.
    target_indices : length N, index of each perturbation's target gene (-1 if absent).
    centroid : optional (G,) reference centroid (e.g. `Source_Grand_Mean` for centroid-referenced
               DRF metrics matching `round8_r86_miller_drf.py`). If None, uses zeros (control-referenced).
    shared_axis : optional (G,) PC1 axis for `v5` / `disc_v5` orthogonal residual evaluation.
    """
    P = np.asarray(pred_deltas, dtype=np.float64)
    T = np.asarray(true_deltas, dtype=np.float64)
    n, G = T.shape
    t_idx = np.asarray(target_indices, dtype=np.int64)
    ref = np.zeros(G, dtype=np.float64) if centroid is None else np.asarray(centroid, dtype=np.float64)

    Q = None
    if shared_axis is not None:
        M = np.column_stack([np.ones(G, dtype=np.float64), ref, np.asarray(shared_axis, dtype=np.float64)])
        Q_full, R_full = np.linalg.qr(M, mode="reduced")
        Q = Q_full[:, np.abs(np.diag(R_full)) > 1e-10 * max(1.0, np.abs(R_full).max())]

    mse_rows = np.empty(n, dtype=np.float64)
    delta_r_rows = np.empty(n, dtype=np.float64)
    de50_r_rows = np.empty(n, dtype=np.float64)
    de50_sign_rows = np.empty(n, dtype=np.float64)
    v5_rows = np.empty(n, dtype=np.float64) if Q is not None else None

    T_disc = np.empty((n, G), dtype=np.float64)
    P_disc = np.empty((n, G), dtype=np.float64)
    T_v5_mat = np.empty((n, G), dtype=np.float64) if Q is not None else None
    P_v5_mat = np.empty((n, G), dtype=np.float64) if Q is not None else None

    k_de = min(int(top_k_de), G - 1) if G > 1 else 1
    all_target_set = set(int(i) for i in t_idx if 0 <= i < G)

    for i in range(n):
        ti = int(t_idx[i])
        mask = np.ones(G, dtype=bool)
        if 0 <= ti < G:
            mask[ti] = False

        yt = T[i]
        yp = P[i]
        dt = yt - ref
        dp = yp - ref

        mse_rows[i] = float(np.mean((yt[mask] - yp[mask]) ** 2))
        dt_m = dt[mask]
        dp_m = dp[mask]
        delta_r_rows[i] = pearson_1d(dt_m, dp_m)

        top_idx = np.argsort(np.abs(dt_m))[-k_de:]
        de50_r_rows[i] = pearson_1d(dt_m[top_idx], dp_m[top_idx])
        s_dt = np.sign(dt_m[top_idx])
        s_dp = np.sign(dp_m[top_idx])
        de50_sign_rows[i] = float(
            np.mean(np.where(s_dp == 0, 0.5, (s_dt == s_dp).astype(np.float64)))
        )

        # Discrimination vector preparation
        dt_d = dt.copy()
        dp_d = dp.copy()
        if mask_target_in_disc and 0 <= ti < G:
            dt_d[ti] = 0.0
            dp_d[ti] = 0.0
        T_disc[i] = dt_d
        P_disc[i] = dp_d

        if Q is not None:
            vt = dt - Q @ (Q.T @ dt)
            vp = dp - Q @ (Q.T @ dp)
            v5_rows[i] = pearson_1d(vt[mask], vp[mask])
            if mask_target_in_disc and 0 <= ti < G:
                vt_d = vt.copy()
                vp_d = vp.copy()
                vt_d[ti] = 0.0
                vp_d[ti] = 0.0
                T_v5_mat[i] = vt_d
                P_v5_mat[i] = vp_d
            else:
                T_v5_mat[i] = vt
                P_v5_mat[i] = vp

    out: Dict[str, np.ndarray] = {
        "mse": mse_rows,
        "delta_r": delta_r_rows,
        "de50_r": de50_r_rows,
        "de50_sign": de50_sign_rows,
        "disc_raw": tie_aware_disc_rows(P_disc, T_disc),
    }

    # A-T3 disc_raw_panel (all evaluated panel targets masked, centered Pearson) & cell-eval2 pds_cosine_panel
    if len(all_target_set) > 0 and len(all_target_set) < G:
        keep_cols = np.ones(G, dtype=bool)
        for j in all_target_set:
            keep_cols[j] = False
        T_panel_ref = (T - ref[None, :])[:, keep_cols]
        P_panel_ref = (P - ref[None, :])[:, keep_cols]
        out["disc_raw_panel"] = tie_aware_disc_rows(P_panel_ref, T_panel_ref)

        # cell-eval2 pds_cosine uses uncentered cosine against control!
        T_panel = T[:, keep_cols]
        P_panel = P[:, keep_cols]
        Tn = T_panel / np.maximum(np.linalg.norm(T_panel, axis=1, keepdims=True), 1e-12)
        Pn = P_panel / np.maximum(np.linalg.norm(P_panel, axis=1, keepdims=True), 1e-12)
        S_cos = Pn @ Tn.T
        diag_cos = np.diag(S_cos)[:, None]
        sb = np.sum(S_cos > diag_cos + 1e-12, axis=1)
        to = np.sum(np.abs(S_cos - diag_cos) <= 1e-12, axis=1) - 1
        pds = 1.0 - (sb + 0.5 * np.maximum(to, 0)) / max(n - 1, 1)
        p_zero = np.linalg.norm(P_panel, axis=1) < 1e-10
        out["pds_cosine_panel"] = np.where(p_zero, 0.5, pds).astype(np.float64)
    else:
        out["disc_raw_panel"] = out["disc_raw"].copy()

    if Q is not None:
        out["v5"] = v5_rows
        out["disc_v5"] = tie_aware_disc_rows(P_v5_mat, T_v5_mat)
    return out


def drf_score(m_pos: float, m_neg: float, lower_is_better: bool = False) -> float:
    """Miller et al. (Nat. Biotechnol. 2026) Dynamic Range Fraction."""
    if lower_is_better:
        return float((m_neg - m_pos) / max(abs(m_neg), 1e-12))
    return float((m_pos - m_neg) / max(1.0 - m_neg, 1e-12))


def paired_row_bootstrap_ci(
    scores_a: np.ndarray,
    scores_b: np.ndarray,
    n_boot: int = 4000,
    seed: int = 12345,
) -> Dict[str, float | bool]:
    """Single-line paired row-level bootstrap 95% CI for mean(scores_a - scores_b)."""
    a = np.asarray(scores_a, dtype=np.float64)
    b = np.asarray(scores_b, dtype=np.float64)
    if a.shape != b.shape or a.ndim != 1:
        raise ValueError(f"Mismatched score shapes: {a.shape} vs {b.shape}")
    n = len(a)
    d = a - b
    rs = np.random.RandomState(seed)
    idx = rs.choice(n, size=(n_boot, n), replace=True)
    boots = d[idx].mean(axis=1)
    ci_lo, ci_hi = np.percentile(boots, [2.5, 97.5])
    ci975_lo, ci975_hi = np.percentile(boots, [1.25, 98.75])
    return {
        "mean_a": float(a.mean()),
        "mean_b": float(b.mean()),
        "delta": float(d.mean()),
        "ci95_lo": float(ci_lo),
        "ci95_hi": float(ci_hi),
        "ci975_lo": float(ci975_lo),
        "ci975_hi": float(ci975_hi),
        "boot_se": float(boots.std()),
        "ci_gt_zero": bool(ci_lo > 0.0),
        "ci_lt_zero": bool(ci_hi < 0.0),
        "ci975_gt_zero": bool(ci975_lo > 0.0),
        "ci975_lt_zero": bool(ci975_hi < 0.0),
    }


def single_array_bootstrap_ci(
    scores: np.ndarray,
    n_boot: int = 4000,
    seed: int = 12345,
) -> Dict[str, float | bool]:
    """Row-level bootstrap 95% and 97.5% CIs for a single metric array."""
    a = np.asarray(scores, dtype=np.float64)
    n = len(a)
    rs = np.random.RandomState(seed)
    idx = rs.choice(n, size=(n_boot, n), replace=True)
    boots = a[idx].mean(axis=1)
    ci_lo, ci_hi = np.percentile(boots, [2.5, 97.5])
    ci975_lo, ci975_hi = np.percentile(boots, [1.25, 98.75])
    return {
        "mean": float(a.mean()),
        "ci95_lo": float(ci_lo),
        "ci95_hi": float(ci_hi),
        "ci975_lo": float(ci975_lo),
        "ci975_hi": float(ci975_hi),
        "boot_se": float(boots.std()),
        "ci_gt_zero": bool(ci_lo > 0.0),
        "ci_lt_zero": bool(ci_hi < 0.0),
        "ci975_gt_zero": bool(ci975_lo > 0.0),
        "ci975_lt_zero": bool(ci975_hi < 0.0),
    }


def stratified_row_bootstrap_ci(
    scores_a: np.ndarray,
    scores_b: np.ndarray,
    strata: Sequence[str],
    n_boot: int = 4000,
    seed: int = 12345,
) -> Dict[str, float | bool]:
    """Equal-weight macro-averaged stratified paired bootstrap 95% & 97.5% CIs across lines."""
    a = np.asarray(scores_a, dtype=np.float64)
    b = np.asarray(scores_b, dtype=np.float64)
    st = np.asarray(strata)
    uniq = sorted(set(st.tolist()))
    rs = np.random.RandomState(seed)
    d = a - b
    boot_sum = np.zeros(n_boot, dtype=np.float64)
    pt_a = 0.0
    pt_b = 0.0
    pt_d = 0.0
    for u in uniq:
        idx_u = np.where(st == u)[0]
        pt_a += float(a[idx_u].mean())
        pt_b += float(b[idx_u].mean())
        pt_d += float(d[idx_u].mean())
        boot_sum += d[rs.choice(idx_u, size=(n_boot, len(idx_u)), replace=True)].mean(axis=1)
    k = float(len(uniq))
    boots = boot_sum / k
    ci_lo, ci_hi = np.percentile(boots, [2.5, 97.5])
    ci975_lo, ci975_hi = np.percentile(boots, [1.25, 98.75])
    return {
        "mean_a": pt_a / k,
        "mean_b": pt_b / k,
        "delta": pt_d / k,
        "ci95_lo": float(ci_lo),
        "ci95_hi": float(ci_hi),
        "ci975_lo": float(ci975_lo),
        "ci975_hi": float(ci975_hi),
        "boot_se": float(boots.std()),
        "ci_gt_zero": bool(ci_lo > 0.0),
        "ci_lt_zero": bool(ci_hi < 0.0),
        "ci975_gt_zero": bool(ci975_lo > 0.0),
        "ci975_lt_zero": bool(ci975_hi < 0.0),
    }


def cis_non_overlapping(ci1: Tuple[float, float], ci2: Tuple[float, float]) -> bool:
    """Return True iff closed 95% confidence intervals [lo1, hi1] and [lo2, hi2] are disjoint."""
    return bool(ci1[1] < ci2[0] or ci2[1] < ci1[0])


def _mask_indices_by_mode(
    delta_obs: np.ndarray,
    r_true: np.ndarray,
    r_pred: np.ndarray,
    target_gene_idx: Sequence[int],
    mode: str = "resp_top50",
    top_k: int = 50,
) -> Tuple[List[np.ndarray], List[np.ndarray], List[np.ndarray], List[np.ndarray]]:
    """
    Build per-perturbation masked residual arrays (xp, yp) and raw deltas (d_obs_p, d_pred_p)
    matching A-T5 `audit/LIT/LINE_SPECIFIC_METRIC.md` §2.2 & §4.
    """
    P, G = delta_obs.shape
    t_idx = np.asarray(target_gene_idx, dtype=np.int64)
    xs_list: List[np.ndarray] = []
    ys_list: List[np.ndarray] = []
    obs_list: List[np.ndarray] = []
    # Also retain per-perturbation OLS slopes and Pearson correlations
    per_pert_slopes: List[float] = []
    per_pert_corrs: List[float] = []

    for p in range(P):
        valid = np.isfinite(r_true[p]) & np.isfinite(r_pred[p]) & np.isfinite(delta_obs[p])
        ti = int(t_idx[p])
        if 0 <= ti < G:
            valid[ti] = False
        idx = np.flatnonzero(valid)
        if idx.size == 0:
            continue
        if mode == "resp_top50" and idx.size > top_k:
            order = np.argsort(np.abs(delta_obs[p, idx]))[-top_k:]
            idx = idx[order]
        elif mode not in ("resp_top50", "all_nontarget"):
            raise ValueError(f"Unknown line-specific mask mode: {mode}")

        xp = r_true[p, idx].astype(np.float64)
        yp = r_pred[p, idx].astype(np.float64)
        xs_list.append(xp)
        ys_list.append(yp)
        obs_list.append(idx)

        xv = xp - xp.mean()
        yv = yp - yp.mean()
        vx = float(np.dot(xv, xv))
        vy = float(np.dot(yv, yv))
        per_pert_slopes.append(float(np.dot(xv, yv) / vx) if vx > 1e-12 else float("nan"))
        per_pert_corrs.append(
            float(np.dot(xv, yv) / np.sqrt(vx * vy)) if (vx > 1e-12 and vy > 1e-12) else float("nan")
        )
    return xs_list, ys_list, per_pert_slopes, per_pert_corrs


def _sufficient_stats_per_pert(
    xs_list: Sequence[np.ndarray],
    ys_list: Sequence[np.ndarray],
) -> Dict[str, np.ndarray]:
    """Precompute exact per-perturbation sufficient statistics for fast row-bootstrap."""
    m = len(xs_list)
    n = np.empty(m, dtype=np.float64)
    sx = np.empty(m, dtype=np.float64)
    sy = np.empty(m, dtype=np.float64)
    sxx = np.empty(m, dtype=np.float64)
    syy = np.empty(m, dtype=np.float64)
    sxy = np.empty(m, dtype=np.float64)
    for i, (xp, yp) in enumerate(zip(xs_list, ys_list)):
        n[i] = float(xp.size)
        sx[i] = float(np.sum(xp))
        sy[i] = float(np.sum(yp))
        sxx[i] = float(np.dot(xp, xp))
        syy[i] = float(np.dot(yp, yp))
        sxy[i] = float(np.dot(xp, yp))
    return {"n": n, "sx": sx, "sy": sy, "sxx": sxx, "syy": syy, "sxy": sxy}


def _stats_from_sums(
    N: np.ndarray | float,
    Sx: np.ndarray | float,
    Sy: np.ndarray | float,
    Sxx: np.ndarray | float,
    Syy: np.ndarray | float,
    Sxy: np.ndarray | float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute (slope_ols, slope_no_intercept, resid_corr) from sufficient-statistic sums."""
    cov_xy = Sxy - (Sx * Sy) / np.maximum(N, 1.0)
    var_x = Sxx - (Sx * Sx) / np.maximum(N, 1.0)
    var_y = Syy - (Sy * Sy) / np.maximum(N, 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        slope_ols = np.where(var_x > 1e-12, cov_xy / var_x, np.nan)
        slope_noint = np.where(Sxx > 1e-12, Sxy / Sxx, np.nan)
        denom_r = np.sqrt(np.maximum(var_x, 0.0) * np.maximum(var_y, 0.0))
        corr = np.where((var_x > 1e-12) & (var_y > 1e-12), cov_xy / denom_r, np.nan)
    return slope_ols, slope_noint, corr


def _extended_scale_stats_from_sums(
    N: np.ndarray | float,
    Sx: np.ndarray | float,
    Sy: np.ndarray | float,
    Sxx: np.ndarray | float,
    Syy: np.ndarray | float,
    Sxy: np.ndarray | float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute PREREG_AUDIT_R1 R6 scale diagnostics from sufficient statistics:
      - sd_ratio = sd(r_pred) / sd(r_obs)
      - slope_obs_on_pred = cov(r_obs, r_pred) / var(r_pred)
      - slope_deming = orthogonal Deming regression slope (error variance ratio delta = 1)
    """
    cov_xy = np.asarray(Sxy - (Sx * Sy) / np.maximum(N, 1.0), dtype=np.float64)
    var_x = np.maximum(np.asarray(Sxx - (Sx * Sx) / np.maximum(N, 1.0), dtype=np.float64), 0.0)
    var_y = np.maximum(np.asarray(Syy - (Sy * Sy) / np.maximum(N, 1.0), dtype=np.float64), 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        sd_ratio = np.where(var_x > 1e-12, np.sqrt(var_y / var_x), np.nan)
        slope_obs_on_pred = np.where(var_y > 1e-12, cov_xy / var_y, np.nan)

        diff_xy = var_x - var_y
        discrim = np.sqrt(diff_xy * diff_xy + 4.0 * cov_xy * cov_xy)
        deming_lt = (2.0 * cov_xy) / (diff_xy + discrim)  # numerically stable when var_y < var_x
        deming_ge = (-diff_xy + discrim) / np.where(np.abs(cov_xy) > 1e-12, 2.0 * cov_xy, np.nan)
        both_zero = (var_x <= 1e-12) & (var_y <= 1e-12)
        slope_deming = np.where(both_zero, np.nan, np.where(var_y < var_x, deming_lt, deming_ge))
    return sd_ratio, slope_obs_on_pred, slope_deming


def compute_line_specific_metrics_and_boots(
    delta_obs: np.ndarray,
    delta_pred: np.ndarray,
    delta_src_mean: np.ndarray,
    target_gene_idx: Sequence[int],
    top_k: int = 50,
    n_boot: int = 4000,
    seed: int = 12345,
) -> Dict[str, Dict[str, object]]:
    """
    Compute official Arc VCC 2026 Wrap-Up line-specific slope & residual correlation metrics
    matching A-T5 `audit/LIT/LINE_SPECIFIC_METRIC.md` §2 & §4.

    Returns dict keyed by mask mode ("resp_top50", "all_nontarget") containing point
    estimates, per-perturbation arrays, and exact 4,000-resample perturbation-level
    bootstrap vectors (enabling paired differences and disjoint-CI checks).
    """
    r_true = np.asarray(delta_obs, dtype=np.float64) - np.asarray(delta_src_mean, dtype=np.float64)
    r_pred = np.asarray(delta_pred, dtype=np.float64) - np.asarray(delta_src_mean, dtype=np.float64)
    d_obs = np.asarray(delta_obs, dtype=np.float64)
    d_pred = np.asarray(delta_pred, dtype=np.float64)
    t_idx = np.asarray(target_gene_idx, dtype=np.int64)
    P, G = d_obs.shape

    # Wrap-Up companion diagnostic: predicted vs measured cosine agreement on resp_top50
    cos_agree = []
    for p in range(P):
        valid = np.isfinite(d_obs[p]) & np.isfinite(d_pred[p])
        ti = int(t_idx[p])
        if 0 <= ti < G:
            valid[ti] = False
        idx = np.flatnonzero(valid)
        if idx.size == 0:
            continue
        if idx.size > top_k:
            idx = idx[np.argsort(np.abs(d_obs[p, idx]))[-top_k:]]
        a = d_obs[p, idx]
        b = d_pred[p, idx]
        den = float(np.linalg.norm(a) * np.linalg.norm(b))
        cos_agree.append(float(np.dot(a, b) / den) if den > 1e-12 else 0.0)

    out: Dict[str, Dict[str, object]] = {}
    for mode in ("resp_top50", "all_nontarget"):
        xs_list, ys_list, p_slopes, p_corrs = _mask_indices_by_mode(
            d_obs, r_true, r_pred, t_idx, mode=mode, top_k=top_k
        )
        m = len(xs_list)
        if m == 0:
            raise ValueError(f"No valid perturbations for mode={mode}")
        ss = _sufficient_stats_per_pert(xs_list, ys_list)
        N_sum = float(ss["n"].sum())
        sx_sum = float(ss["sx"].sum())
        sy_sum = float(ss["sy"].sum())
        sxx_sum = float(ss["sxx"].sum())
        syy_sum = float(ss["syy"].sum())
        sxy_sum = float(ss["sxy"].sum())
        pt_ols, pt_noint, pt_corr = _stats_from_sums(
            N_sum, sx_sum, sy_sum, sxx_sum, syy_sum, sxy_sum
        )
        pt_sd_ratio, pt_obs_on_pred, pt_deming = _extended_scale_stats_from_sums(
            N_sum, sx_sum, sy_sum, sxx_sum, syy_sum, sxy_sum
        )

        # Deterministic perturbation-level bootstrap with identical seed across models for pairing!
        rs = np.random.RandomState(seed)
        draws = rs.choice(m, size=(n_boot, m), replace=True)
        bN = ss["n"][draws].sum(axis=1)
        bsx = ss["sx"][draws].sum(axis=1)
        bsy = ss["sy"][draws].sum(axis=1)
        bsxx = ss["sxx"][draws].sum(axis=1)
        bsyy = ss["syy"][draws].sum(axis=1)
        bsxy = ss["sxy"][draws].sum(axis=1)
        b_ols, b_noint, b_corr = _stats_from_sums(bN, bsx, bsy, bsxx, bsyy, bsxy)
        b_sd_ratio, b_obs_on_pred, b_deming = _extended_scale_stats_from_sums(
            bN, bsx, bsy, bsxx, bsyy, bsxy
        )
        p_slopes_arr = np.asarray(p_slopes, dtype=np.float64)
        p_corrs_arr = np.asarray(p_corrs, dtype=np.float64)
        finite_s = np.isfinite(p_slopes_arr)
        finite_c = np.isfinite(p_corrs_arr)
        if np.any(finite_s):
            b_pslope = np.nanmean(p_slopes_arr[draws], axis=1)
            ps_mean = float(np.nanmean(p_slopes_arr))
            ps_ci = np.nanpercentile(b_pslope, [2.5, 97.5]).tolist()
        else:
            b_pslope = np.full(n_boot, np.nan, dtype=np.float64)
            ps_mean = float("nan")
            ps_ci = [float("nan"), float("nan")]
        if np.any(finite_c):
            b_pcorr = np.nanmean(p_corrs_arr[draws], axis=1)
            pc_mean = float(np.nanmean(p_corrs_arr))
            pc_ci = np.nanpercentile(b_pcorr, [2.5, 97.5]).tolist()
        else:
            b_pcorr = np.full(n_boot, np.nan, dtype=np.float64)
            pc_mean = float("nan")
            pc_ci = [float("nan"), float("nan")]

        def _safe_ci(arr: np.ndarray, qs: Tuple[float, float] = (2.5, 97.5)) -> List[float]:
            if np.any(np.isfinite(arr)):
                return np.nanpercentile(arr, list(qs)).tolist()
            return [float("nan"), float("nan")]

        out[mode] = {
            "n_perts": int(m),
            "slope_ols": float(pt_ols),
            "slope_ols_ci95": _safe_ci(b_ols, (2.5, 97.5)),
            "slope_ols_ci975": _safe_ci(b_ols, (1.25, 98.75)),
            "slope_no_intercept": float(pt_noint),
            "slope_no_intercept_ci95": _safe_ci(b_noint, (2.5, 97.5)),
            "slope_no_intercept_ci975": _safe_ci(b_noint, (1.25, 98.75)),
            "per_pert_slope_mean": ps_mean,
            "per_pert_slope_ci95": ps_ci,
            "per_pert_slope_ci975": _safe_ci(b_pslope, (1.25, 98.75)),
            "resid_corr": float(pt_corr),
            "resid_corr_ci95": _safe_ci(b_corr, (2.5, 97.5)),
            "resid_corr_ci975": _safe_ci(b_corr, (1.25, 98.75)),
            "sd_ratio": float(pt_sd_ratio),
            "sd_ratio_ci95": _safe_ci(b_sd_ratio, (2.5, 97.5)),
            "sd_ratio_ci975": _safe_ci(b_sd_ratio, (1.25, 98.75)),
            "slope_obs_on_pred": float(pt_obs_on_pred),
            "slope_obs_on_pred_ci95": _safe_ci(b_obs_on_pred, (2.5, 97.5)),
            "slope_obs_on_pred_ci975": _safe_ci(b_obs_on_pred, (1.25, 98.75)),
            "slope_deming": float(pt_deming),
            "slope_deming_ci95": _safe_ci(b_deming, (2.5, 97.5)),
            "slope_deming_ci975": _safe_ci(b_deming, (1.25, 98.75)),
            "per_pert_corr_mean": pc_mean,
            "per_pert_corr_ci95": pc_ci,
            "pred_vs_obs_cosine_mean": float(np.mean(cos_agree)) if cos_agree else float("nan"),
            "boot_slope_ols": b_ols,
            "boot_slope_no_intercept": b_noint,
            "boot_per_pert_slope": b_pslope,
            "boot_resid_corr": b_corr,
            "boot_sd_ratio": b_sd_ratio,
            "boot_slope_obs_on_pred": b_obs_on_pred,
            "boot_slope_deming": b_deming,
        }
    return out


def spearman_brown_ceiling(rho_half: float) -> Tuple[float, float]:
    """
    Compute Spearman-Brown full-depth reliability rho_SB and correlation ceiling sqrt(rho_SB).
    Returns (nan, nan) if rho_half <= 0 or non-finite, matching LINE_SPECIFIC_METRIC.md §2.4.
    """
    if not np.isfinite(rho_half) or rho_half <= 0.0:
        return float("nan"), float("nan")
    rho_sb = (2.0 * rho_half) / (1.0 + rho_half)
    return float(rho_sb), float(np.sqrt(rho_sb))


def compute_splithalf_noise_ceiling(
    delta_half_a: np.ndarray,
    delta_half_b: np.ndarray,
    delta_obs_full: np.ndarray,
    delta_src_mean: np.ndarray,
    target_gene_idx: Sequence[int],
    top_k: int = 50,
    n_boot: int = 4000,
    seed: int = 12345,
) -> Dict[str, Dict[str, object]]:
    """
    Compute split-half line-specific reliability (rho_AB, rho_SB, ceiling sqrt(rho_SB),
    and replicate OLS slope) plus bootstrap ceiling array for both `resp_top50`
    and `all_nontarget`, using identical perturbation-level bootstrap draws (`seed`).
    """
    r_a = np.asarray(delta_half_a, dtype=np.float64) - np.asarray(delta_src_mean, dtype=np.float64)
    r_b = np.asarray(delta_half_b, dtype=np.float64) - np.asarray(delta_src_mean, dtype=np.float64)
    d_obs = np.asarray(delta_obs_full, dtype=np.float64)
    out: Dict[str, Dict[str, object]] = {}
    for mode in ("resp_top50", "all_nontarget"):
        xs_list, ys_list, _, _ = _mask_indices_by_mode(
            d_obs, r_a, r_b, target_gene_idx, mode=mode, top_k=top_k
        )
        m = len(xs_list)
        ss = _sufficient_stats_per_pert(xs_list, ys_list)
        slope_ab, _, rho_ab = _stats_from_sums(
            float(ss["n"].sum()),
            float(ss["sx"].sum()),
            float(ss["sy"].sum()),
            float(ss["sxx"].sum()),
            float(ss["syy"].sum()),
            float(ss["sxy"].sum()),
        )
        # Swap A/B for symmetric replicate slope
        slope_ba, _, _ = _stats_from_sums(
            float(ss["n"].sum()),
            float(ss["sy"].sum()),
            float(ss["sx"].sum()),
            float(ss["syy"].sum()),
            float(ss["sxx"].sum()),
            float(ss["sxy"].sum()),
        )
        rho_ab = float(rho_ab)
        slope_ab = float(slope_ab)
        slope_ba = float(slope_ba)
        rho_sb, ceil_r = spearman_brown_ceiling(rho_ab)

        rs = np.random.RandomState(seed)
        draws = rs.choice(m, size=(n_boot, m), replace=True)
        b_slope_ab, _, b_rho_ab = _stats_from_sums(
            ss["n"][draws].sum(axis=1),
            ss["sx"][draws].sum(axis=1),
            ss["sy"][draws].sum(axis=1),
            ss["sxx"][draws].sum(axis=1),
            ss["syy"][draws].sum(axis=1),
            ss["sxy"][draws].sum(axis=1),
        )
        with np.errstate(divide="ignore", invalid="ignore"):
            b_rho_sb = np.where(b_rho_ab > 0.0, (2.0 * b_rho_ab) / (1.0 + b_rho_ab), np.nan)
            b_ceil_r = np.sqrt(b_rho_sb)

        out[mode] = {
            "rho_half_ab": rho_ab,
            "rho_sb": rho_sb,
            "ceiling_r": ceil_r,
            "slope_ab": slope_ab,
            "slope_ba": slope_ba,
            "slope_rep_mean": float(0.5 * (slope_ab + slope_ba))
            if (np.isfinite(slope_ab) and np.isfinite(slope_ba))
            else float("nan"),
            "boot_ceiling_r": b_ceil_r,
            "boot_slope_ab": b_slope_ab,
        }
    return out

