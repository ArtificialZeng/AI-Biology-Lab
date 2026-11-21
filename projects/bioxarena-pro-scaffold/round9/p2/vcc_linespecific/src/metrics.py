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
    seed: int = 20261121,
) -> Dict[str, float | bool]:
    """Single-line paired row-level bootstrap 95% CI for mean(scores_a - scores_b)."""
    a = np.asarray(scores_a, dtype=np.float64)
    b = np.asarray(scores_b, dtype=np.float64)
    if a.shape != b.shape or a.ndim != 1:
        raise ValueError(f"Mismatched score shapes: {a.shape} vs {b.shape}")
    n = len(a)
    d = a - b
    rs = np.random.default_rng(seed)
    idx = rs.integers(0, n, size=(n_boot, n))
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
    seed: int = 20261121,
) -> Dict[str, float | bool]:
    """Row-level bootstrap 95% and 97.5% CIs for a single metric array."""
    a = np.asarray(scores, dtype=np.float64)
    n = len(a)
    rs = np.random.default_rng(seed)
    idx = rs.integers(0, n, size=(n_boot, n))
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
    seed: int = 20261121,
) -> Dict[str, float | bool]:
    """Equal-weight macro-averaged stratified paired bootstrap 95% & 97.5% CIs across lines."""
    a = np.asarray(scores_a, dtype=np.float64)
    b = np.asarray(scores_b, dtype=np.float64)
    st = np.asarray(strata)
    uniq = sorted(set(st.tolist()))
    rs = np.random.default_rng(seed)
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
        sel = idx_u[rs.integers(0, len(idx_u), size=(n_boot, len(idx_u)))]
        boot_sum += d[sel].mean(axis=1)
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


def _orth_project(vec: np.ndarray, base: np.ndarray) -> np.ndarray:
    """Remove per-perturbation component along `base` on the same masked gene set (R17)."""
    norm = float(np.linalg.norm(base))
    if norm <= 1e-12 or not np.isfinite(norm):
        return vec.copy()
    b = base / norm
    return vec - float(np.dot(vec, b)) * b


def _mask_indices_by_mode(
    delta_obs: np.ndarray,
    r_true: np.ndarray,
    r_pred: np.ndarray,
    target_gene_idx: Sequence[int],
    mode: str = "resp_top50",
    top_k: int = 50,
    delta_src_mean: Optional[np.ndarray] = None,
) -> Tuple[
    List[np.ndarray],
    List[np.ndarray],
    List[float],
    List[float],
    List[np.ndarray],
    List[np.ndarray],
]:
    """
    Build per-perturbation masked residual arrays (xp, yp) and amplitude-orthogonalized
    arrays (xp_orth, yp_orth) matching A-T5 `LINE_SPECIFIC_METRIC.md` §2.2/§4 and R17.
    """
    P, G = delta_obs.shape
    t_idx = np.asarray(target_gene_idx, dtype=np.int64)
    xs_list: List[np.ndarray] = []
    ys_list: List[np.ndarray] = []
    xs_orth_list: List[np.ndarray] = []
    ys_orth_list: List[np.ndarray] = []
    per_pert_slopes: List[float] = []
    per_pert_corrs: List[float] = []

    for p in range(P):
        valid = np.isfinite(r_true[p]) & np.isfinite(r_pred[p]) & np.isfinite(delta_obs[p])
        if delta_src_mean is not None:
            valid = valid & np.isfinite(delta_src_mean[p])
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

        if delta_src_mean is not None:
            bp = delta_src_mean[p, idx].astype(np.float64)
            xs_orth_list.append(_orth_project(xp, bp))
            ys_orth_list.append(_orth_project(yp, bp))
        else:
            xs_orth_list.append(xp)
            ys_orth_list.append(yp)

        xv = xp - xp.mean()
        yv = yp - yp.mean()
        vx = float(np.dot(xv, xv))
        vy = float(np.dot(yv, yv))
        per_pert_slopes.append(float(np.dot(xv, yv) / vx) if vx > 1e-12 else float("nan"))
        per_pert_corrs.append(
            float(np.dot(xv, yv) / np.sqrt(vx * vy)) if (vx > 1e-12 and vy > 1e-12) else float("nan")
        )
    return xs_list, ys_list, per_pert_slopes, per_pert_corrs, xs_orth_list, ys_orth_list


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
    seed: int = 20261121,
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
    d_src_arr = np.asarray(delta_src_mean, dtype=np.float64)
    for mode in ("resp_top50", "all_nontarget"):
        xs_list, ys_list, p_slopes, p_corrs, xs_orth, ys_orth = _mask_indices_by_mode(
            d_obs, r_true, r_pred, t_idx, mode=mode, top_k=top_k, delta_src_mean=d_src_arr
        )
        m = len(xs_list)
        if m == 0:
            raise ValueError(f"No valid perturbations for mode={mode}")
        ss = _sufficient_stats_per_pert(xs_list, ys_list)
        ss_orth = _sufficient_stats_per_pert(xs_orth, ys_orth)
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
        _, _, pt_corr_orth = _stats_from_sums(
            float(ss_orth["n"].sum()),
            float(ss_orth["sx"].sum()),
            float(ss_orth["sy"].sum()),
            float(ss_orth["sxx"].sum()),
            float(ss_orth["syy"].sum()),
            float(ss_orth["sxy"].sum()),
        )

        # Deterministic perturbation-level bootstrap with identical seed across models for pairing!
        rs = np.random.default_rng(seed)
        draws = rs.integers(0, m, size=(n_boot, m))
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
        _, _, b_corr_orth = _stats_from_sums(
            ss_orth["n"][draws].sum(axis=1),
            ss_orth["sx"][draws].sum(axis=1),
            ss_orth["sy"][draws].sum(axis=1),
            ss_orth["sxx"][draws].sum(axis=1),
            ss_orth["syy"][draws].sum(axis=1),
            ss_orth["sxy"][draws].sum(axis=1),
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
            "resid_corr_orth": float(pt_corr_orth),
            "resid_corr_orth_ci95": _safe_ci(b_corr_orth, (2.5, 97.5)),
            "resid_corr_orth_ci975": _safe_ci(b_corr_orth, (1.25, 98.75)),
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
            "boot_resid_corr_orth": b_corr_orth,
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
    seed: int = 20261121,
) -> Dict[str, Dict[str, object]]:
    """
    Compute split-half line-specific reliability (raw and R17 amplitude-orthogonalized)
    for both `resp_top50` and `all_nontarget`, using identical perturbation draws (`seed`).
    """
    d_src_arr = np.asarray(delta_src_mean, dtype=np.float64)
    r_a = np.asarray(delta_half_a, dtype=np.float64) - d_src_arr
    r_b = np.asarray(delta_half_b, dtype=np.float64) - d_src_arr
    d_obs = np.asarray(delta_obs_full, dtype=np.float64)
    out: Dict[str, Dict[str, object]] = {}
    for mode in ("resp_top50", "all_nontarget"):
        xs_list, ys_list, _, _, xs_orth, ys_orth = _mask_indices_by_mode(
            d_obs, r_a, r_b, target_gene_idx, mode=mode, top_k=top_k, delta_src_mean=d_src_arr
        )
        m = len(xs_list)
        ss = _sufficient_stats_per_pert(xs_list, ys_list)
        ss_orth = _sufficient_stats_per_pert(xs_orth, ys_orth)
        slope_ab, _, rho_ab = _stats_from_sums(
            float(ss["n"].sum()),
            float(ss["sx"].sum()),
            float(ss["sy"].sum()),
            float(ss["sxx"].sum()),
            float(ss["syy"].sum()),
            float(ss["sxy"].sum()),
        )
        slope_ba, _, _ = _stats_from_sums(
            float(ss["n"].sum()),
            float(ss["sy"].sum()),
            float(ss["sx"].sum()),
            float(ss["syy"].sum()),
            float(ss["sxx"].sum()),
            float(ss["sxy"].sum()),
        )
        _, _, rho_ab_orth = _stats_from_sums(
            float(ss_orth["n"].sum()),
            float(ss_orth["sx"].sum()),
            float(ss_orth["sy"].sum()),
            float(ss_orth["sxx"].sum()),
            float(ss_orth["syy"].sum()),
            float(ss_orth["sxy"].sum()),
        )
        rho_ab = float(rho_ab)
        rho_ab_orth = float(rho_ab_orth)
        slope_ab = float(slope_ab)
        slope_ba = float(slope_ba)
        rho_sb, ceil_r = spearman_brown_ceiling(rho_ab)
        rho_sb_orth, ceil_r_orth = spearman_brown_ceiling(rho_ab_orth)

        rs = np.random.default_rng(seed)
        draws = rs.integers(0, m, size=(n_boot, m))
        b_slope_ab, _, b_rho_ab = _stats_from_sums(
            ss["n"][draws].sum(axis=1),
            ss["sx"][draws].sum(axis=1),
            ss["sy"][draws].sum(axis=1),
            ss["sxx"][draws].sum(axis=1),
            ss["syy"][draws].sum(axis=1),
            ss["sxy"][draws].sum(axis=1),
        )
        _, _, b_rho_ab_orth = _stats_from_sums(
            ss_orth["n"][draws].sum(axis=1),
            ss_orth["sx"][draws].sum(axis=1),
            ss_orth["sy"][draws].sum(axis=1),
            ss_orth["sxx"][draws].sum(axis=1),
            ss_orth["syy"][draws].sum(axis=1),
            ss_orth["sxy"][draws].sum(axis=1),
        )
        with np.errstate(divide="ignore", invalid="ignore"):
            b_rho_sb = np.where(b_rho_ab > 0.0, (2.0 * b_rho_ab) / (1.0 + b_rho_ab), np.nan)
            b_ceil_r = np.sqrt(b_rho_sb)
            b_rho_sb_orth = np.where(
                b_rho_ab_orth > 0.0, (2.0 * b_rho_ab_orth) / (1.0 + b_rho_ab_orth), np.nan
            )
            b_ceil_r_orth = np.sqrt(b_rho_sb_orth)

        out[mode] = {
            "rho_half_ab": rho_ab,
            "rho_sb": rho_sb,
            "ceiling_r": ceil_r,
            "rho_half_ab_orth": rho_ab_orth,
            "rho_sb_orth": rho_sb_orth,
            "ceiling_r_orth": ceil_r_orth,
            "slope_ab": slope_ab,
            "slope_ba": slope_ba,
            "slope_rep_mean": float(0.5 * (slope_ab + slope_ba))
            if (np.isfinite(slope_ab) and np.isfinite(slope_ba))
            else float("nan"),
            "boot_ceiling_r": b_ceil_r,
            "boot_ceiling_r_orth": b_ceil_r_orth,
            "boot_slope_ab": b_slope_ab,
        }
    return out


def compute_tmvs_shuffled_var_perm100(
    delta_obs: np.ndarray,
    delta_src_mean: np.ndarray,
    target_gene_idx: Sequence[int],
    target_diag_var: np.ndarray,
    mean_src_var: np.ndarray,
    ceil_r_top50: float,
    ceil_r_orth_top50: float,
    n_perm: int = 100,
    top_k: int = 50,
    seed: int = 42,
    min_var: float = 1e-8,
) -> Dict[str, object]:
    """
    PREREG_AUDIT_R3 R18: Evaluate 100 random gene permutations of target marginal variances
    to obtain null distributions (median, p2.5, p97.5, p98.75) for `r_orth_cal`, `r_cal`, and `slope_ols`.
    """
    d_obs = np.asarray(delta_obs, dtype=np.float64)
    d_src = np.asarray(delta_src_mean, dtype=np.float64)
    t_idx = np.asarray(target_gene_idx, dtype=np.int64)
    diag_t = np.maximum(np.asarray(target_diag_var, dtype=np.float64), min_var)
    denom_var = np.maximum(np.asarray(mean_src_var, dtype=np.float64), min_var)
    P, G = d_obs.shape

    row_masks: List[np.ndarray] = []
    s_list: List[np.ndarray] = []
    r_obs_list: List[np.ndarray] = []
    b_list: List[np.ndarray] = []
    r_obs_orth_list: List[np.ndarray] = []

    for p in range(P):
        valid = np.isfinite(d_obs[p]) & np.isfinite(d_src[p])
        ti = int(t_idx[p])
        if 0 <= ti < G:
            valid[ti] = False
        idx = np.flatnonzero(valid)
        if idx.size == 0:
            continue
        if idx.size > top_k:
            idx = idx[np.argsort(np.abs(d_obs[p, idx]))[-top_k:]]
        row_masks.append(idx)
        sp = d_src[p, idx]
        rp = d_obs[p, idx] - sp
        s_list.append(sp)
        r_obs_list.append(rp)
        norm = float(np.linalg.norm(sp))
        bp = sp / norm if norm > 1e-12 else np.zeros_like(sp)
        b_list.append(bp)
        r_obs_orth_list.append(rp - float(np.dot(rp, bp)) * bp)

    ro_flat = np.concatenate(r_obs_list)
    ro_orth_flat = np.concatenate(r_obs_orth_list)
    ro_c = ro_flat - ro_flat.mean()
    ro_var = float(np.dot(ro_c, ro_c))
    ro_orth_c = ro_orth_flat - ro_orth_flat.mean()
    ro_orth_var = float(np.dot(ro_orth_c, ro_orth_c))

    uniform_len = len({len(idx) for idx in row_masks}) == 1
    if uniform_len:
        idx_mat = np.vstack(row_masks)
        S_mat = np.vstack(s_list)
        B_mat = np.vstack(b_list)

    r_orth_cal_vals = np.full(n_perm, np.nan, dtype=np.float64)
    r_cal_vals = np.full(n_perm, np.nan, dtype=np.float64)
    slope_vals = np.full(n_perm, np.nan, dtype=np.float64)

    for k in range(n_perm):
        perm = np.random.default_rng(k).permutation(G)
        scale_perm = np.sqrt(diag_t[perm] / denom_var) - 1.0
        if uniform_len:
            rp_mat = S_mat * scale_perm[idx_mat]
            dot_pb = np.sum(rp_mat * B_mat, axis=1, keepdims=True)
            rp_orth_flat = (rp_mat - dot_pb * B_mat).ravel()
            rp_flat = rp_mat.ravel()
        else:
            rp_parts = []
            rp_orth_parts = []
            for idx, sp, bp in zip(row_masks, s_list, b_list):
                rp = sp * scale_perm[idx]
                rp_parts.append(rp)
                rp_orth_parts.append(rp - float(np.dot(rp, bp)) * bp)
            rp_flat = np.concatenate(rp_parts)
            rp_orth_flat = np.concatenate(rp_orth_parts)

        rp_c = rp_flat - rp_flat.mean()
        rp_var = float(np.dot(rp_c, rp_c))
        if ro_var > 1e-12:
            cov = float(np.dot(ro_c, rp_c))
            slope_vals[k] = cov / ro_var
            if rp_var > 1e-12 and np.isfinite(ceil_r_top50) and ceil_r_top50 > 0:
                r_cal_vals[k] = (cov / np.sqrt(ro_var * rp_var)) / ceil_r_top50

        rp_oc = rp_orth_flat - rp_orth_flat.mean()
        rp_ovar = float(np.dot(rp_oc, rp_oc))
        if ro_orth_var > 1e-12 and rp_ovar > 1e-12 and np.isfinite(ceil_r_orth_top50) and ceil_r_orth_top50 > 0:
            cov_o = float(np.dot(ro_orth_c, rp_oc))
            r_orth_cal_vals[k] = (cov_o / np.sqrt(ro_orth_var * rp_ovar)) / ceil_r_orth_top50

    def _summarize(arr: np.ndarray) -> Dict[str, object]:
        fin = arr[np.isfinite(arr)]
        if fin.size == 0:
            return {
                "p2_5": float("nan"),
                "median": float("nan"),
                "p97_5": float("nan"),
                "p98_75": float("nan"),
                "min": float("nan"),
                "max": float("nan"),
            }
        return {
            "p2_5": float(np.percentile(fin, 2.5)),
            "median": float(np.percentile(fin, 50.0)),
            "p97_5": float(np.percentile(fin, 97.5)),
            "p98_75": float(np.percentile(fin, 98.75)),
            "min": float(np.min(fin)),
            "max": float(np.max(fin)),
        }

    return {
        "n_perm": int(n_perm),
        "r_orth_cal": _summarize(r_orth_cal_vals),
        "r_cal": _summarize(r_cal_vals),
        "slope_ols": _summarize(slope_vals),
    }


def _solve_r19_batch(
    W_dense: np.ndarray,
    nv: np.ndarray,
    ss_sum: np.ndarray,
    yy_sum: np.ndarray,
    s_sq_sum: np.ndarray,
    sy_sum: np.ndarray,
    y_sq_sum: np.ndarray,
    rr_sum: np.ndarray,
    r_sq_sum: np.ndarray,
    cnt_m: np.ndarray,
    S_m: np.ndarray,
    S2_m: np.ndarray,
    SY_m: np.ndarray,
    S2Y_m: np.ndarray,
    S2Y2_m: np.ndarray,
    S3_m: np.ndarray,
    S3Y_m: np.ndarray,
    S4_m: np.ndarray,
    min_perts_per_gene: int = 5,
    eps: float = 1e-12,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    PREREG_A v4 §5.4 vectorized 3-tier SST(r_obs) variance decomposition with distinct-row
    gating (C_dist >= min_perts_per_gene) and cluster-W^2 finite-df noise-overfit subtraction.
    """
    N = np.maximum(W_dense @ nv, 1.0)
    Ss = W_dense @ ss_sum
    Sy = W_dense @ yy_sum
    Sss = W_dense @ s_sq_sum
    Ssy = W_dense @ sy_sum
    Syy = W_dense @ y_sq_sum
    Sr = W_dense @ rr_sum
    Srr = W_dense @ r_sq_sum

    cov_sy = Ssy - (Ss * Sy) / N
    vx = Sss - (Ss * Ss) / N
    beta = np.where(vx > eps, cov_sy / np.maximum(vx, eps), 1.0)
    alpha = (Sy - beta * Ss) / N
    sst_r = np.maximum(Srr - (Sr * Sr) / N, eps)
    sse_1 = np.maximum(
        Syy - 2.0 * alpha * Sy + N * alpha * alpha - 2.0 * beta * (Ssy - alpha * Ss) + (beta * beta) * Sss,
        0.0,
    )
    f1 = np.clip(1.0 - sse_1 / sst_r, 0.0, 1.0)

    C_dist = (W_dense > 0.0).astype(np.float64) @ cnt_m
    Ss_g = W_dense @ S_m
    Sss_g = W_dense @ S2_m
    Ssy_g = W_dense @ SY_m

    a_col = alpha[:, None]
    b_col = beta[:, None]
    ok_g = (C_dist >= min_perts_per_gene) & (Sss_g > eps)
    mg = np.where(ok_g, (Ssy_g - a_col * Ss_g) / np.maximum(Sss_g, eps), b_col)
    imp_raw = np.where(ok_g, ((mg - b_col) ** 2) * Sss_g, 0.0)

    W2 = W_dense * W_dense
    S2 = W2 @ S2_m
    S2Y = W2 @ S2Y_m
    S2Y2 = W2 @ S2Y2_m
    S3 = W2 @ S3_m
    S3Y = W2 @ S3Y_m
    S4 = W2 @ S4_m
    num_w2 = (
        S2Y2
        - 2.0 * a_col * S2Y
        + (a_col * a_col) * S2
        - 2.0 * mg * (S3Y - a_col * S3)
        + (mg * mg) * S4
    )
    neff_g = np.where(S4 > eps, (Sss_g * Sss_g) / np.maximum(S4, eps), C_dist)
    df_adj = np.maximum(1.0 - 1.0 / np.maximum(neff_g, 2.0), 0.25)
    bias_g = np.where(ok_g, np.maximum(num_w2, 0.0) / (np.maximum(Sss_g, eps) * df_adj), 0.0)
    imp_unbiased = np.maximum(np.sum(imp_raw - bias_g, axis=1), 0.0)
    sse_2 = np.maximum(sse_1 - imp_unbiased, 0.0)

    f12 = np.clip(1.0 - sse_2 / sst_r, f1, 1.0)
    f2 = np.clip(f12 - f1, 0.0, 1.0)
    f3 = np.clip(1.0 - f12, 0.0, 1.0)
    return beta, f1, f2, f3


def compute_residual_variance_decomposition(
    delta_obs: np.ndarray,
    delta_src_mean: np.ndarray,
    target_gene_idx: Sequence[int],
    top_k: int = 50,
    n_boot: int = 4000,
    seed: int = 20261121,
    min_perts_per_gene: int = 5,
) -> Dict[str, Dict[str, object]]:
    """
    PREREG_AUDIT_R3 R19 (PREREG_A v4 §5.4): Descriptive 3-tier SST(r_obs) variance decomposition
    of R = Y - X for both `resp_top50` and `all_nontarget`:
      1. Tier 1 (`frac_line_gain`): Whole-line pooled OLS alpha + beta_t * X
      2. Tier 2 (`frac_gene_gain`): Incremental per-gene slopes m_{t,g} beyond beta_t
         (`C_distinct >= min_perts_per_gene` with cluster-W^2 finite-df noise-overfit subtraction)
      3. Tier 3 (`frac_remaining`): Remaining residual variance fraction
    """
    Y = np.asarray(delta_obs, dtype=np.float64)
    X = np.asarray(delta_src_mean, dtype=np.float64)
    t_idx = np.asarray(target_gene_idx, dtype=np.int64)
    P, G = Y.shape

    rs = np.random.default_rng(seed)
    draws = rs.integers(0, P, size=(n_boot, P))

    out: Dict[str, Dict[str, object]] = {}
    for mode in ("resp_top50", "all_nontarget"):
        M = np.isfinite(Y) & np.isfinite(X)
        for p in range(P):
            ti = int(t_idx[p])
            if 0 <= ti < G:
                M[p, ti] = False
            if mode == "resp_top50":
                idx = np.flatnonzero(M[p])
                if idx.size > top_k:
                    keep = idx[np.argsort(np.abs(Y[p, idx]))[-top_k:]]
                    M[p] = False
                    M[p, keep] = True

        S_m = np.where(M, X, 0.0)
        Y_m = np.where(M, Y, 0.0)
        R_m = Y_m - S_m
        cnt_m = M.astype(np.float64)

        S2_m = S_m * S_m
        SY_m = S_m * Y_m
        Y2_m = Y_m * Y_m
        S2Y_m = S2_m * Y_m
        S2Y2_m = S2_m * Y2_m
        S3_m = S2_m * S_m
        S3Y_m = S3_m * Y_m
        S4_m = S2_m * S2_m

        nv = cnt_m.sum(axis=1)
        ss_sum = S_m.sum(axis=1)
        yy_sum = Y_m.sum(axis=1)
        s_sq_sum = S2_m.sum(axis=1)
        sy_sum = SY_m.sum(axis=1)
        y_sq_sum = Y2_m.sum(axis=1)
        rr_sum = R_m.sum(axis=1)
        r_sq_sum = np.einsum("pg,pg->p", R_m, R_m)

        pt_beta, pt_f1, pt_f2, pt_f3 = _solve_r19_batch(
            np.ones((1, P), dtype=np.float64),
            nv, ss_sum, yy_sum, s_sq_sum, sy_sum, y_sq_sum, rr_sum, r_sq_sum,
            cnt_m, S_m, S2_m, SY_m, S2Y_m, S2Y2_m, S3_m, S3Y_m, S4_m,
            min_perts_per_gene=min_perts_per_gene,
        )

        # Also compute no-intercept multiplicative ratio beta_mult = sum(X*Y)/sum(X^2) for diagnostics
        tot_xx_pt = float(s_sq_sum.sum())
        pt_beta_mult = float(sy_sum.sum() / tot_xx_pt) if tot_xx_pt > 1e-12 else float("nan")

        b_beta = np.empty(n_boot, dtype=np.float64)
        b_beta_mult = np.empty(n_boot, dtype=np.float64)
        b_f1 = np.empty(n_boot, dtype=np.float64)
        b_f2 = np.empty(n_boot, dtype=np.float64)
        b_f3 = np.empty(n_boot, dtype=np.float64)

        chunk = 400
        for s0 in range(0, n_boot, chunk):
            s1 = min(s0 + chunk, n_boot)
            d_chunk = draws[s0:s1]
            W = np.zeros((s1 - s0, P), dtype=np.float64)
            rows_idx = np.repeat(np.arange(s1 - s0), P)
            np.add.at(W, (rows_idx, d_chunk.ravel()), 1.0)

            bb, bf1, bf2, bf3 = _solve_r19_batch(
                W,
                nv, ss_sum, yy_sum, s_sq_sum, sy_sum, y_sq_sum, rr_sum, r_sq_sum,
                cnt_m, S_m, S2_m, SY_m, S2Y_m, S2Y2_m, S3_m, S3Y_m, S4_m,
                min_perts_per_gene=min_perts_per_gene,
            )
            b_beta[s0:s1] = bb
            b_f1[s0:s1] = bf1
            b_f2[s0:s1] = bf2
            b_f3[s0:s1] = bf3
            wxx = W @ s_sq_sum
            wxy = W @ sy_sum
            b_beta_mult[s0:s1] = np.where(wxx > 1e-12, wxy / np.maximum(wxx, 1e-12), np.nan)

        def _ci_pair(val: float, arr: np.ndarray) -> Dict[str, object]:
            return {
                "estimate": float(val),
                "ci95": np.nanpercentile(arr, [2.5, 97.5]).tolist()
                if np.any(np.isfinite(arr))
                else [float("nan"), float("nan")],
                "ci975": np.nanpercentile(arr, [1.25, 98.75]).tolist()
                if np.any(np.isfinite(arr))
                else [float("nan"), float("nan")],
            }

        out[mode] = {
            "n_perts": int(P),
            "min_gene_obs_for_gene_layer": int(min_perts_per_gene),
            "n_eligible_genes": int(np.sum((cnt_m.sum(axis=0) >= min_perts_per_gene) & (S2_m.sum(axis=0) > 1e-12))),
            "beta_ols_int": _ci_pair(float(pt_beta[0]), b_beta),
            "beta_mult": _ci_pair(pt_beta_mult, b_beta_mult),
            "frac_line_gain": _ci_pair(float(pt_f1[0]), b_f1),
            "frac_gene_gain": _ci_pair(float(pt_f2[0]), b_f2),
            "frac_remaining": _ci_pair(float(pt_f3[0]), b_f3),
        }
    return out


