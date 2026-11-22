#!/usr/bin/env python3
"""
Round 10 Evaluation Metrics, Joint Perturbation Bootstrap (RNG=20261122),
100-Perm FCT Null, R19 3-Tier Variance Decomposition, and S7 Confounder Controls.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.stats import spearmanr

from fct_pipeline import FctSourceState, fct_from_fc_bar

EPS = 1e-12


def pearson_1d(a: np.ndarray, b: np.ndarray) -> float:
    """Safe 1D Pearson correlation clamped to [-1, 1]; returns 0.0 on zero variance."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.size < 2:
        return 0.0
    ac = a - a.mean()
    bc = b - b.mean()
    denom = float(np.linalg.norm(ac) * np.linalg.norm(bc))
    if denom <= EPS or not np.isfinite(denom):
        return 0.0
    v = float(np.dot(ac, bc) / denom)
    return float(np.clip(v, -1.0, 1.0)) if np.isfinite(v) else 0.0


def pearson_1d_strict(a: np.ndarray, b: np.ndarray) -> float:
    """1D Pearson correlation returning NaN on degenerate/zero variance."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.size < 2:
        return float("nan")
    ac = a - a.mean()
    bc = b - b.mean()
    denom = float(np.linalg.norm(ac) * np.linalg.norm(bc))
    if denom <= EPS or not np.isfinite(denom):
        return float("nan")
    return float(np.clip(np.dot(ac, bc) / denom, -1.0, 1.0))


def ols_slope_1d(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if x.size < 2:
        return float("nan")
    xc = x - x.mean()
    vx = float(np.dot(xc, xc))
    if vx <= EPS:
        return float("nan")
    return float(np.dot(xc, y - y.mean()) / vx)


def spearman_brown_ceiling(rho_ab: float) -> float:
    """Spearman-Brown noise ceiling sqrt(2*rho / (1 + rho)); NaN when rho <= 0."""
    if not np.isfinite(rho_ab) or rho_ab <= 0.0:
        return float("nan")
    sb = (2.0 * float(rho_ab)) / (1.0 + float(rho_ab))
    return float(np.sqrt(max(sb, 0.0)))


def disattenuated_r(r_raw: float, ceiling_r: float) -> float:
    if not (np.isfinite(r_raw) and np.isfinite(ceiling_r)) or ceiling_r <= 1e-8:
        return float("nan")
    return float(np.clip(float(r_raw) / float(ceiling_r), -1.0, 1.0))


def summarize_boot(
    point: float,
    boots: np.ndarray,
    n_items: int,
    n_boot: int,
) -> Dict[str, Any]:
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


def tie_aware_disc_rows(P_mat: np.ndarray, T_mat: np.ndarray, zero_tol: float = 1e-10) -> np.ndarray:
    """Tie-aware mid-rank discrimination percentile matching Round 9 `metrics.py`."""
    n = T_mat.shape[0]
    if n <= 1:
        return np.full(n, 0.5, dtype=np.float64)
    A0 = T_mat - T_mat.mean(axis=1, keepdims=True)
    B0 = P_mat - P_mat.mean(axis=1, keepdims=True)
    An = A0 / np.maximum(np.linalg.norm(A0, axis=1, keepdims=True), EPS)
    Bn = B0 / np.maximum(np.linalg.norm(B0, axis=1, keepdims=True), EPS)
    C = Bn @ An.T
    diag = np.diag(C)[:, None]
    strictly_better = np.sum(C > diag + EPS, axis=1)
    tied_others = np.sum(np.abs(C - diag) <= EPS, axis=1) - 1
    eff_better = strictly_better + 0.5 * np.maximum(tied_others, 0)
    pcts = 1.0 - eff_better / max(n - 1, 1)
    b_zero = np.linalg.norm(B0, axis=1) < zero_tol
    return np.where(b_zero, 0.5, pcts).astype(np.float64)


def compute_row_metrics_r10(
    pred_deltas: np.ndarray,
    true_deltas: np.ndarray,
    target_indices: Sequence[int],
    centroid: Optional[np.ndarray] = None,
    shared_axis: Optional[np.ndarray] = None,
    top_k_de: int = 50,
    mask_target_in_disc: bool = True,
) -> Dict[str, np.ndarray]:
    """Compute per-perturbation aggregate conventions (`H-C3`: disc_raw, disc_raw_panel, de50_r, pds_cosine_panel)."""
    P = np.asarray(pred_deltas, dtype=np.float64)
    T = np.asarray(true_deltas, dtype=np.float64)
    n, G = T.shape
    t_idx = np.asarray(target_indices, dtype=np.int64)
    ref = np.zeros(G, dtype=np.float64) if centroid is None else np.asarray(centroid, dtype=np.float64)

    Q = None
    if shared_axis is not None:
        M = np.column_stack([np.ones(G, dtype=np.float64), ref, np.asarray(shared_axis, dtype=np.float64)])
        Q_full, R_full = np.linalg.qr(M, mode="reduced")
        Q = Q_full[:, np.abs(np.diag(R_full)) > 1e-10 * max(1.0, float(np.abs(R_full).max()))]

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
    all_target_set = {int(i) for i in t_idx if 0 <= i < G}

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
        de50_sign_rows[i] = float(np.mean(np.where(s_dp == 0, 0.5, (s_dt == s_dp).astype(np.float64))))

        dt_d = dt.copy()
        dp_d = dp.copy()
        if mask_target_in_disc and 0 <= ti < G:
            dt_d[ti] = 0.0
            dp_d[ti] = 0.0
        T_disc[i] = dt_d
        P_disc[i] = dp_d

        if Q is not None and v5_rows is not None and T_v5_mat is not None and P_v5_mat is not None:
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

    if 0 < len(all_target_set) < G:
        keep_cols = np.ones(G, dtype=bool)
        for j in all_target_set:
            keep_cols[j] = False
        T_panel_ref = (T - ref[None, :])[:, keep_cols]
        P_panel_ref = (P - ref[None, :])[:, keep_cols]
        out["disc_raw_panel"] = tie_aware_disc_rows(P_panel_ref, T_panel_ref)

        T_panel = T[:, keep_cols]
        P_panel = P[:, keep_cols]
        Tn = T_panel / np.maximum(np.linalg.norm(T_panel, axis=1, keepdims=True), EPS)
        Pn = P_panel / np.maximum(np.linalg.norm(P_panel, axis=1, keepdims=True), EPS)
        S_cos = Pn @ Tn.T
        diag_cos = np.diag(S_cos)[:, None]
        sb = np.sum(S_cos > diag_cos + EPS, axis=1)
        to = np.sum(np.abs(S_cos - diag_cos) <= EPS, axis=1) - 1
        pds = 1.0 - (sb + 0.5 * np.maximum(to, 0)) / max(n - 1, 1)
        p_zero = np.linalg.norm(P_panel, axis=1) < 1e-10
        out["pds_cosine_panel"] = np.where(p_zero, 0.5, pds).astype(np.float64)
    else:
        out["disc_raw_panel"] = out["disc_raw"].copy()

    if Q is not None and v5_rows is not None and T_v5_mat is not None and P_v5_mat is not None:
        out["v5"] = v5_rows
        out["disc_v5"] = tie_aware_disc_rows(P_v5_mat, T_v5_mat)
    return out


def paired_row_bootstrap_ci_r10(
    a: np.ndarray,
    b: np.ndarray,
    n_boot: int = 4000,
    seed: int = 20261122,
) -> Dict[str, Any]:
    """Perturbation-level paired bootstrap CI using `default_rng(seed)` (`C2`)."""
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

    a_q0125, a_q025, a_q975, a_q9875 = np.percentile(boot_a, [1.25, 2.5, 97.5, 98.75])
    b_q0125, b_q025, b_q975, b_q9875 = np.percentile(boot_b, [1.25, 2.5, 97.5, 98.75])
    d_q0125, d_q025, d_q975, d_q9875 = np.percentile(boot_d, [1.25, 2.5, 97.5, 98.75])
    mean_a = float(np.mean(xa))
    mean_b = float(np.mean(xb))
    return {
        "mean_a": mean_a,
        "mean_b": mean_b,
        "delta": float(mean_a - mean_b),
        "a_ci95_lo": float(a_q025),
        "a_ci95_hi": float(a_q975),
        "a_ci975_lo": float(a_q0125),
        "a_ci975_hi": float(a_q9875),
        "b_ci95_lo": float(b_q025),
        "b_ci95_hi": float(b_q975),
        "b_ci975_lo": float(b_q0125),
        "b_ci975_hi": float(b_q9875),
        "ci95_lo": float(d_q025),
        "ci95_hi": float(d_q975),
        "ci975_lo": float(d_q0125),
        "ci975_hi": float(d_q9875),
        "disjoint_ci95": bool(a_q025 > b_q975),
        "disjoint_ci975": bool(a_q0125 > b_q9875),
        "ci95_gt_zero": bool(d_q025 > 0.0),
        "ci95_lt_zero": bool(d_q975 < 0.0),
        "ci975_gt_zero": bool(d_q0125 > 0.0),
        "ci975_lt_zero": bool(d_q9875 < 0.0),
        "n_rows": int(n),
    }


def _orth_pair(x: np.ndarray, y: np.ndarray, ref: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    nrm = float(np.linalg.norm(ref))
    if nrm <= EPS:
        return x.copy(), y.copy()
    b = ref / nrm
    return x - float(np.dot(x, b)) * b, y - float(np.dot(y, b)) * b


def _moments(xs: Sequence[np.ndarray], ys: Sequence[np.ndarray]) -> Tuple[np.ndarray, ...]:
    nv = np.array([float(x.size) for x in xs], dtype=np.float64)
    sx = np.array([float(np.sum(x)) for x in xs], dtype=np.float64)
    sy = np.array([float(np.sum(y)) for y in ys], dtype=np.float64)
    sxx = np.array([float(np.dot(x, x)) for x in xs], dtype=np.float64)
    syy = np.array([float(np.dot(y, y)) for y in ys], dtype=np.float64)
    sxy = np.array([float(np.dot(x, y)) for x, y in zip(xs, ys)], dtype=np.float64)
    return nv, sx, sy, sxx, syy, sxy


def _moments_to_corr_and_slope(
    nv: np.ndarray | float,
    sx: np.ndarray | float,
    sy: np.ndarray | float,
    sxx: np.ndarray | float,
    syy: np.ndarray | float,
    sxy: np.ndarray | float,
    zero_degenerate_nan: bool = True,
) -> Tuple[np.ndarray | float, np.ndarray | float, np.ndarray | float, np.ndarray | float, np.ndarray | float]:
    """Compute (slope_ols, r, sd_ratio, rev_ols, deming_slope) from aggregated moments."""
    n_safe = np.maximum(nv, 1.0)
    cov = sxy - (sx * sy) / n_safe
    vx = sxx - (sx * sx) / n_safe
    vy = syy - (sy * sy) / n_safe
    with np.errstate(divide="ignore", invalid="ignore"):
        slope_ols = np.where(vx > EPS, cov / vx, np.nan if zero_degenerate_nan else 0.0)
        denom = np.sqrt(np.maximum(vx, 0.0) * np.maximum(vy, 0.0))
        corr = np.where(denom > EPS, np.clip(cov / denom, -1.0, 1.0), np.nan)
        sd_ratio = np.where(vx > EPS, np.sqrt(np.maximum(vy, 0.0) / np.maximum(vx, EPS)), np.nan)
        rev_ols = np.where(vy > EPS, cov / vy, np.nan)
        diff = vy - vx
        deming = np.where(
            np.abs(cov) > EPS,
            (diff + np.sqrt(diff * diff + 4.0 * cov * cov)) / (2.0 * cov),
            np.where(vx > EPS, 0.0, np.nan),
        )
    return slope_ols, corr, sd_ratio, rev_ols, deming


def _extract_masked_idx_per_pert(
    delta_obs: np.ndarray,
    target_gene_idx: Sequence[int],
    extra_valid: Optional[np.ndarray] = None,
    gene_keep_mask: Optional[np.ndarray] = None,
    mode: str = "resp_top50",
    top_k: int = 50,
) -> Tuple[List[int], List[np.ndarray]]:
    N, G = delta_obs.shape
    rows: List[int] = []
    cols_list: List[np.ndarray] = []
    for p in range(N):
        valid = np.isfinite(delta_obs[p])
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
        if mode == "resp_top50" and idx.size > top_k:
            idx = idx[np.argsort(np.abs(delta_obs[p, idx]), kind="mergesort")[-top_k:]]
        rows.append(p)
        cols_list.append(idx)
    return rows, cols_list


def compute_residual_suite_r10(
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
    seed: int = 20261122,
) -> Dict[str, Any]:
    """
    Compute Group C primary (`r_orth_cal` + joint bootstrap 97.5%/95% CIs across
    `resp_top50` & `all_nontarget`) and secondary slope / raw `r_cal` / cosine metrics.
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

    # Mean per-perturbation cosine similarity on resp_top50 (Wrap-Up diagnostic)
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
            "r_orth_cal": summarize_boot(pt_orth_cal, b_orth_cal, m, n_boot),
            "r_orth_cal_sq": summarize_boot(pt_orth_cal_sq, b_orth_cal_sq, m, n_boot),
            "r_orth_raw": summarize_boot(float(pt_r_orth), b_r_orth, m, n_boot),
            "r_cal": summarize_boot(pt_cal, b_cal, m, n_boot),
            "r_cal_sq": summarize_boot(pt_cal_sq, b_cal_sq, m, n_boot),
            "r_raw": summarize_boot(float(pt_r_raw), b_r_raw, m, n_boot),
            "slope_ols": summarize_boot(float(pt_slope), b_slope, m, n_boot),
            "sd_pred_over_sd_obs": summarize_boot(float(pt_sd_ratio), b_sd_ratio, m, n_boot),
            "slope_deming": summarize_boot(float(pt_dem), b_dem, m, n_boot),
            "slope_obs_on_pred": summarize_boot(float(pt_rev), b_rev, m, n_boot),
            "per_pert_slope_mean": float(np.mean(finite_slopes)) if finite_slopes else float("nan"),
            "per_pert_corr_mean": float(np.mean(finite_corrs)) if finite_corrs else float("nan"),
            "rho_half_ab": float(pt_rho_ab),
            "ceiling_r": float(pt_ceil),
            "rho_half_ab_orth": float(pt_rho_ab_orth),
            "ceiling_r_orth": float(pt_ceil_orth),
        }

    return out


def compute_joint_paired_r_orth_diff_r10(
    delta_obs: np.ndarray,
    delta_pred_a: np.ndarray,
    delta_pred_b: np.ndarray,
    delta_transfer: np.ndarray,
    target_gene_idx: Sequence[int],
    delta_half_a: np.ndarray,
    delta_half_b: np.ndarray,
    half_valid: Optional[np.ndarray] = None,
    mode: str = "resp_top50",
    top_k: int = 50,
    n_boot: int = 4000,
    seed: int = 20261122,
) -> Dict[str, Any]:
    """Paired joint-bootstrap CI for `r_orth_cal(A) - r_orth_cal(B)` (e.g., `Hybrid_on_TMVS` - `TMVS`)."""
    r_obs = np.asarray(delta_obs, dtype=np.float64) - np.asarray(delta_transfer, dtype=np.float64)
    r_a = np.asarray(delta_pred_a, dtype=np.float64) - np.asarray(delta_transfer, dtype=np.float64)
    r_b = np.asarray(delta_pred_b, dtype=np.float64) - np.asarray(delta_transfer, dtype=np.float64)
    t_mean = np.asarray(delta_transfer, dtype=np.float64)
    ha_res = np.asarray(delta_half_a, dtype=np.float64) - t_mean
    hb_res = np.asarray(delta_half_b, dtype=np.float64) - t_mean

    rows, idx_list = _extract_masked_idx_per_pert(
        delta_obs, target_gene_idx, extra_valid=r_a + r_b, mode=mode, top_k=top_k
    )
    m = len(rows)
    if m == 0:
        return summarize_boot(float("nan"), np.zeros(0), 0, n_boot)

    xa_o, ya_o, xb_o, yb_o, ha_o, hb_o = [], [], [], [], [], []
    for p, idx in zip(rows, idx_list):
        xv = r_obs[p, idx]
        rv = t_mean[p, idx]
        xa, ya = _orth_pair(xv, r_a[p, idx], rv)
        xb, yb = _orth_pair(xv, r_b[p, idx], rv)
        xa_o.append(xa); ya_o.append(ya)
        xb_o.append(xb); yb_o.append(yb)
        hv_ok = half_valid is None or bool(half_valid[p])
        if hv_ok and np.all(np.isfinite(ha_res[p, idx])) and np.all(np.isfinite(hb_res[p, idx])):
            hao, hbo = _orth_pair(ha_res[p, idx], hb_res[p, idx], rv)
            ha_o.append(hao); hb_o.append(hbo)
        else:
            z = np.zeros(0, dtype=np.float64)
            ha_o.append(z); hb_o.append(z)

    nva, sxa, sya, sxxa, syya, sxya = _moments(xa_o, ya_o)
    nvb, sxb, syb, sxxb, syyb, sxyb = _moments(xb_o, yb_o)
    hnv, hsx, hsy, hsxx, hsyy, hsxy = _moments(ha_o, hb_o)

    _, pt_ra, _, _, _ = _moments_to_corr_and_slope(
        float(nva.sum()), float(sxa.sum()), float(sya.sum()), float(sxxa.sum()), float(syya.sum()), float(sxya.sum())
    )
    _, pt_rb, _, _, _ = _moments_to_corr_and_slope(
        float(nvb.sum()), float(sxb.sum()), float(syb.sum()), float(sxxb.sum()), float(syyb.sum()), float(sxyb.sum())
    )
    _, pt_hab, _, _, _ = _moments_to_corr_and_slope(
        float(hnv.sum()), float(hsx.sum()), float(hsy.sum()), float(hsxx.sum()), float(hsyy.sum()), float(hsxy.sum())
    )
    pt_ceil = spearman_brown_ceiling(float(pt_hab))
    pt_diff = disattenuated_r(float(pt_ra), pt_ceil) - disattenuated_r(float(pt_rb), pt_ceil)

    rng = np.random.default_rng(seed)
    draws = rng.integers(0, m, size=(n_boot, m))
    _, b_ra, _, _, _ = _moments_to_corr_and_slope(
        nva[draws].sum(axis=1), sxa[draws].sum(axis=1), sya[draws].sum(axis=1),
        sxxa[draws].sum(axis=1), syya[draws].sum(axis=1), sxya[draws].sum(axis=1),
    )
    _, b_rb, _, _, _ = _moments_to_corr_and_slope(
        nvb[draws].sum(axis=1), sxb[draws].sum(axis=1), syb[draws].sum(axis=1),
        sxxb[draws].sum(axis=1), syyb[draws].sum(axis=1), sxyb[draws].sum(axis=1),
    )
    _, b_hab, _, _, _ = _moments_to_corr_and_slope(
        hnv[draws].sum(axis=1), hsx[draws].sum(axis=1), hsy[draws].sum(axis=1),
        hsxx[draws].sum(axis=1), hsyy[draws].sum(axis=1), hsxy[draws].sum(axis=1),
    )
    b_diff = np.array(
        [
            disattenuated_r(float(ra), spearman_brown_ceiling(float(hab)))
            - disattenuated_r(float(rb), spearman_brown_ceiling(float(hab)))
            for ra, rb, hab in zip(b_ra, b_rb, b_hab)
        ],
        dtype=np.float64,
    )
    return summarize_boot(pt_diff, b_diff, m, n_boot)


def compute_fct_perm100_null(
    state: FctSourceState,
    target_ctrl_mean_libnorm: np.ndarray,
    query_target_names: Sequence[str],
    delta_obs: np.ndarray,
    delta_transfer: np.ndarray,
    target_gene_idx: Sequence[int],
    ceil_r_top50: float,
    ceil_r_orth_top50: float,
    fc_bundle: Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]] = None,
    n_perm: int = 100,
    top_k: int = 50,
    row_subset_mask: Optional[np.ndarray] = None,
    gene_keep_mask: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """
    100-permutation null distribution for `FCT` (`NegCtrl_FCT_Perm(k)`, `k = 0..n_perm-1`):
    permuting `target_ctrl_mean_libnorm` across `G` genes with `np.random.default_rng(k)`
    while keeping source `fc_bar` fixed (`PREREG_C §3.2`, `H-C2`, `C2`).
    """
    mu_t = np.maximum(np.asarray(target_ctrl_mean_libnorm, dtype=np.float64), 0.0)
    G = mu_t.shape[0]
    if fc_bundle is None:
        from fct_pipeline import compute_source_transfer_and_fc_bar
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
        d_obs_s, t_gidx_s, gene_keep_mask=gene_keep_mask, mode="resp_top50", top_k=top_k
    )
    if not rows:
        return {"n_perm": n_perm, "r_orth_cal": {}, "r_cal": {}, "slope_ols": {}}

    xo_raw_l, xo_orth_l, unit_ref_l = [], [], []
    for p, idx in zip(rows, idx_list):
        xv = r_obs[p, idx]
        rv = t_mean_s[p, idx].astype(np.float64)
        nrm = float(np.linalg.norm(rv))
        bv = (rv / nrm) if nrm > EPS else np.zeros_like(rv)
        xo_raw_l.append(xv)
        xo_orth_l.append(xv - float(np.dot(xv, bv)) * bv)
        unit_ref_l.append(bv)

    x_raw = np.concatenate(xo_raw_l)
    x_orth = np.concatenate(xo_orth_l)
    r_orth_cal_vals = np.empty(n_perm, dtype=np.float64)
    r_cal_vals = np.empty(n_perm, dtype=np.float64)
    slope_vals = np.empty(n_perm, dtype=np.float64)

    for k in range(n_perm):
        perm = np.random.default_rng(k).permutation(G)
        pred_k = fct_from_fc_bar(mu_t[perm], fc_bar_s, p_trans_s, hits_s)
        res_k = pred_k - t_mean_s
        yp_raw_l, yp_orth_l = [], []
        for p, idx, bv in zip(rows, idx_list, unit_ref_l):
            yv = res_k[p, idx]
            yp_raw_l.append(yv)
            yp_orth_l.append(yv - float(np.dot(yv, bv)) * bv)
        y_raw = np.concatenate(yp_raw_l)
        y_orth = np.concatenate(yp_orth_l)
        slope_vals[k] = ols_slope_1d(x_raw, y_raw)
        r_cal_vals[k] = disattenuated_r(pearson_1d_strict(x_raw, y_raw), ceil_r_top50)
        r_orth_cal_vals[k] = disattenuated_r(pearson_1d_strict(x_orth, y_orth), ceil_r_orth_top50)

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
        "r_orth_cal": _pct(r_orth_cal_vals),
        "r_cal": _pct(r_cal_vals),
        "slope_ols": _pct(slope_vals),
        "perm0_r_orth_cal": float(r_orth_cal_vals[0]),
        "perm0_r_cal": float(r_cal_vals[0]),
        "perm0_slope_ols": float(slope_vals[0]),
    }


def compute_r19_variance_decomposition_r10(
    delta_obs: np.ndarray,
    delta_transfer: np.ndarray,
    target_gene_idx: Sequence[int],
    min_perts_per_gene: int = 5,
    top_k: int = 50,
    n_boot: int = 4000,
    seed: int = 20261122,
    boot_chunk: int = 250,
) -> Dict[str, Any]:
    """
    `PREREG_A v4 §5.4` / `PREREG_C SA-4`: Vectorized 3-tier multiplicative gain decomposition
    of `SST(r_obs)` into L1 (whole-line gain), L2 (per-gene incremental gain for genes in
    >=5 distinct perturbations with exact cluster-df W^2 bias correction), and L3 (remaining residual).
    """
    Y_obs = np.asarray(delta_obs, dtype=np.float64)
    X_tm = np.asarray(delta_transfer, dtype=np.float64)
    N_all, G = Y_obs.shape
    out: Dict[str, Any] = {}

    for mode in ("resp_top50", "all_nontarget"):
        rows, idx_list = _extract_masked_idx_per_pert(Y_obs, target_gene_idx, mode=mode, top_k=top_k)
        m = len(rows)
        if m == 0:
            out[mode] = {}
            continue

        rows_arr = np.asarray(rows, dtype=np.int64)
        M = np.zeros((m, G), dtype=np.float64)
        for ii, cols in enumerate(idx_list):
            M[ii, cols] = 1.0
        Xm = np.nan_to_num(X_tm[rows_arr], nan=0.0) * M
        Ym = np.nan_to_num(Y_obs[rows_arr], nan=0.0) * M
        Rm = Ym - Xm

        # Precompute matrices for BLAS chunked aggregation across (n_boot_chunk, m) @ (m, G)
        M_X = Xm
        M_Y = Ym
        M_XX = Xm * Xm
        M_XY = Xm * Ym
        M_YY = Ym * Ym
        M_RR = Rm * Rm
        M_R = Rm

        def _eval_weight_batch(W: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
            # W shape (B, m) integer/float multiplicities
            N_g = W @ M               # (B, G)
            Sx_g = W @ M_X
            Sy_g = W @ M_Y
            Sxx_g = W @ M_XX
            Sxy_g = W @ M_XY
            Syy_g = W @ M_YY
            C_g = (W > 0).astype(np.float64) @ M  # distinct perturbations per gene

            W2 = W * W
            W2_XX = W2 @ M_XX

            # Tier 1 totals (pooled across all masked entries)
            N_tot = N_g.sum(axis=1)
            Sx_tot = Sx_g.sum(axis=1)
            Sy_tot = Sy_g.sum(axis=1)
            Sxx_tot = Sxx_g.sum(axis=1)
            Sxy_tot = Sxy_g.sum(axis=1)
            Syy_tot = Syy_g.sum(axis=1)
            Sr_tot = (W @ M_R).sum(axis=1)
            Srr_tot = (W @ M_RR).sum(axis=1)

            sst = np.maximum(Srr_tot - (Sr_tot * Sr_tot) / np.maximum(N_tot, 1.0), 0.0)
            mx = Sx_tot / np.maximum(N_tot, 1.0)
            my = Sy_tot / np.maximum(N_tot, 1.0)
            vx_tot = Sxx_tot - (Sx_tot * Sx_tot) / np.maximum(N_tot, 1.0)
            cxy_tot = Sxy_tot - (Sx_tot * Sy_tot) / np.maximum(N_tot, 1.0)
            beta1 = np.where(vx_tot > EPS, cxy_tot / np.maximum(vx_tot, EPS), 1.0)
            a1 = my - beta1 * mx

            # SSE1 = sum((Y - a1 - beta1*X)^2)
            sse1 = (
                Syy_tot
                + (a1 * a1) * N_tot
                + (beta1 * beta1) * Sxx_tot
                - 2.0 * a1 * Sy_tot
                - 2.0 * beta1 * Sxy_tot
                + 2.0 * a1 * beta1 * Sx_tot
            )
            sse1 = np.maximum(sse1, 0.0)

            # Tier 2 within-gene ZERO-INTERCEPT multiplicative incremental slope fit:
            #   r_t1 = Y - (a1 + b1 * X) ~ dm_g * X   (PREREG_A v4 §5.4 / PREREG_C SA-4)
            # Only 1 parameter dm_g per eligible gene (C_distinct >= min_perts_per_gene).
            a1_c = a1[:, None]
            b1_c = beta1[:, None]
            Sxr_g = Sxy_g - a1_c * Sx_g - b1_c * Sxx_g
            Srr_g = (
                Syy_g
                - 2.0 * a1_c * Sy_g
                - 2.0 * b1_c * Sxy_g
                + (a1_c * a1_c) * N_g
                + 2.0 * a1_c * b1_c * Sx_g
                + (b1_c * b1_c) * Sxx_g
            )
            Srr_g = np.maximum(Srr_g, 0.0)

            elig = (C_g >= float(min_perts_per_gene)) & (N_g > 1.0) & (Sxx_g > EPS)
            dmg = np.where(elig, Sxr_g / np.maximum(Sxx_g, EPS), 0.0)
            gain_raw_g = np.where(elig, (dmg * dmg) * Sxx_g, 0.0)
            sse_new_g = np.maximum(Srr_g - gain_raw_g, 0.0)

            # Exact 1-df cluster-correlated overfit bias subtraction using W^2
            df_cluster_g = np.where(
                elig,
                W2_XX / np.maximum(Sxx_g, EPS),
                0.0,
            )
            denom_df = np.maximum(N_g - df_cluster_g, 1.0)
            sig2_g = np.where(elig, sse_new_g / denom_df, 0.0)
            bias_g = sig2_g * df_cluster_g
            gain_net_g = np.where(elig, gain_raw_g - bias_g, 0.0)

            l1 = np.where(sst > EPS, np.clip(1.0 - sse1 / np.maximum(sst, EPS), 0.0, 1.0), np.nan)
            total_gain_adj = np.maximum(0.0, gain_net_g.sum(axis=1))
            l2_raw = np.where(sst > EPS, total_gain_adj / np.maximum(sst, EPS), np.nan)
            l2 = np.where(np.isfinite(l1), np.clip(l2_raw, 0.0, 1.0 - l1), np.nan)
            l3 = np.where(np.isfinite(l1), np.clip(1.0 - l1 - l2, 0.0, 1.0), np.nan)
            return l1, l2, l3

        W_pt = np.ones((1, m), dtype=np.float64)
        l1_pt, l2_pt, l3_pt = _eval_weight_batch(W_pt)

        rng = np.random.default_rng(seed)
        draws = rng.integers(0, m, size=(n_boot, m))
        b_l1 = np.empty(n_boot, dtype=np.float64)
        b_l2 = np.empty(n_boot, dtype=np.float64)
        b_l3 = np.empty(n_boot, dtype=np.float64)
        for s in range(0, n_boot, boot_chunk):
            e = min(s + boot_chunk, n_boot)
            W_chunk = np.zeros((e - s, m), dtype=np.float64)
            row_ids = np.repeat(np.arange(e - s), m)
            np.add.at(W_chunk, (row_ids, draws[s:e].ravel()), 1.0)
            c1, c2, c3 = _eval_weight_batch(W_chunk)
            b_l1[s:e] = c1
            b_l2[s:e] = c2
            b_l3[s:e] = c3

        out[mode] = {
            "n_perts": int(m),
            "L1_whole_line_frac": summarize_boot(float(l1_pt[0]), b_l1, m, n_boot),
            "L2_gene_gain_frac": summarize_boot(float(l2_pt[0]), b_l2, m, n_boot),
            "L3_remaining_resid_frac": summarize_boot(float(l3_pt[0]), b_l3, m, n_boot),
        }
    return out


def compute_s7_cellcount_and_kd_analysis(
    state: FctSourceState,
    query_target_names: Sequence[str],
    target_ctrl_mean_libnorm: np.ndarray,
    fc_bundle: Tuple[np.ndarray, np.ndarray, np.ndarray],
    delta_obs: np.ndarray,
    delta_pred_fct: np.ndarray,
    delta_transfer: np.ndarray,
    target_gene_idx: Sequence[int],
    n_cells: Optional[np.ndarray],
    knockdown_depth: Optional[np.ndarray],
    delta_half_a: Optional[np.ndarray] = None,
    delta_half_b: Optional[np.ndarray] = None,
    half_valid: Optional[np.ndarray] = None,
    mismatch_preds: Optional[Dict[str, np.ndarray]] = None,
    top_k: int = 50,
    n_boot: int = 4000,
    seed: int = 20261122,
) -> Dict[str, Any]:
    """
    Pre-registered `PREREG_C §5 SA-6 / S7` (`P6` Exact Specification):
      (i)   Per-perturbation orthogonalized residual correlation `r_orth_p` across masked top-50 genes.
      (ii)  Spearman correlation of `r_orth_p` with `ln(n_cells_p)` and `|k_p|` (+ 97.5% boot CIs).
      (iii) Mean-centered OLS `r_orth_p ~ [1, ln(n_cells_p) - mean, |k_p| - mean]` reporting
            adjusted intercept (`97.5%` boot CI at observed cohort mean covariates) & slopes,
            plus 1-covariate `ln(n_cells)` companion if any perturbations lack self-gene KD.
      (iv)  Balanced tertile stratification (`T1_low`, `T2_mid`, `T3_high`) by `n_cells`
            recomputing pooled `FCT` `r_orth_cal [97.5% CI]`, mismatched controls & perm100.
    """
    N = delta_obs.shape[0]
    r_obs = np.asarray(delta_obs, dtype=np.float64) - np.asarray(delta_transfer, dtype=np.float64)
    r_pred = np.asarray(delta_pred_fct, dtype=np.float64) - np.asarray(delta_transfer, dtype=np.float64)
    t_mean = np.asarray(delta_transfer, dtype=np.float64)

    # (i) Compute per-perturbation r_orth_p on resp_top50
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

    for p in range(N):
        valid = np.isfinite(delta_obs[p]) & np.isfinite(r_pred[p]) & np.isfinite(t_mean[p])
        ti = int(target_gene_idx[p])
        if 0 <= ti < delta_obs.shape[1]:
            valid[ti] = False
        idx = np.flatnonzero(valid)
        if idx.size == 0:
            continue
        if idx.size > top_k:
            idx = idx[np.argsort(np.abs(delta_obs[p, idx]), kind="mergesort")[-top_k:]]
        xo, yo = _orth_pair(r_obs[p, idx], r_pred[p, idx], t_mean[p, idx])
        r_orth_p[p] = pearson_1d_strict(xo, yo)

    def _boot_spearman(x_vec: np.ndarray, y_vec: np.ndarray) -> Dict[str, Any]:
        ok = np.isfinite(x_vec) & np.isfinite(y_vec)
        xv = x_vec[ok]
        yv = y_vec[ok]
        m = xv.size
        if m < 5 or np.std(xv) <= EPS or np.std(yv) <= EPS:
            return summarize_boot(float("nan"), np.zeros(0), int(m), n_boot)
        pt = float(spearmanr(xv, yv).statistic)
        rng = np.random.default_rng(seed)
        b_vals = np.empty(n_boot, dtype=np.float64)
        for b in range(n_boot):
            idx = rng.integers(0, m, size=m)
            xb, yb = xv[idx], yv[idx]
            if np.std(xb) <= EPS or np.std(yb) <= EPS:
                b_vals[b] = float("nan")
            else:
                b_vals[b] = float(spearmanr(xb, yb).statistic)
        return summarize_boot(pt, b_vals, int(m), n_boot)

    # (ii) Spearman correlations
    spearman_ln_nc = _boot_spearman(r_orth_p, ln_nc_arr)
    spearman_kd = _boot_spearman(r_orth_p, kd_arr)

    # (iii) Covariate-adjusted OLS intercept [1, ln(n_cells) - mean_full, |k_p| - mean_full]
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
        # Center at observed full-sample means across all bootstrap draws so the target adjustment
        # point E[r_orth_p | mean_ln_ncells, mean_|k_p|] stays fixed across replicates
        n_par = X_centered.shape[1]
        b_betas = np.full((n_boot, n_par), np.nan, dtype=np.float64)
        for b in range(n_boot):
            di = draws[b]
            Xb = X_centered[di]
            yb = yv[di]
            try:
                b_betas[b], _, _, _ = np.linalg.lstsq(Xb, yb, rcond=None)
            except np.linalg.LinAlgError:
                pass
        res: Dict[str, Any] = {
            "n_valid": int(m),
            "covariate_means": {name: mu for name, mu in zip(col_names, means)},
            "unadjusted_mean_r_orth_p": summarize_boot(float(np.mean(yv)), yv[draws].mean(axis=1), int(m), n_boot),
            "adjusted_intercept": summarize_boot(float(beta_pt[0]), b_betas[:, 0], int(m), n_boot),
            "slopes": {
                name: summarize_boot(float(beta_pt[j + 1]), b_betas[:, j + 1], int(m), n_boot)
                for j, name in enumerate(col_names)
            },
        }
        return res

    ols_both = _boot_centered_ols(r_orth_p, [ln_nc_arr, kd_arr], ["ln_n_cells", "abs_kd_depth"])
    ols_nc_only = _boot_centered_ols(r_orth_p, [ln_nc_arr], ["ln_n_cells"])

    # (iv) Balanced tertile stratification (`T1_low`, `T2_mid`, `T3_high`) by `n_cells`
    tertiles_out: Dict[str, Any] = {}
    valid_nc = np.isfinite(nc_arr) & (nc_arr > 0.0)
    if int(np.sum(valid_nc)) >= 15:
        v_indices = np.flatnonzero(valid_nc)
        # Deterministic sort by (n_cells, orig_idx)
        order = sorted(v_indices.tolist(), key=lambda i: (float(nc_arr[i]), int(i)))
        chunks = np.array_split(np.asarray(order, dtype=np.int64), 3)
        labels = ("T1_low", "T2_mid", "T3_high")
        for lbl, c_idx in zip(labels, chunks):
            if c_idx.size < 5:
                continue
            mask_t = np.zeros(N, dtype=bool)
            mask_t[c_idx] = True
            fct_res = compute_residual_suite_r10(
                delta_obs=delta_obs[c_idx],
                delta_pred=delta_pred_fct[c_idx],
                delta_transfer=delta_transfer[c_idx],
                target_gene_idx=np.asarray(target_gene_idx, dtype=np.int64)[c_idx],
                delta_half_a=delta_half_a[c_idx] if delta_half_a is not None else None,
                delta_half_b=delta_half_b[c_idx] if delta_half_b is not None else None,
                half_valid=half_valid[c_idx] if half_valid is not None else None,
                top_k=top_k,
                n_boot=n_boot,
                seed=seed,
            )
            t50 = fct_res.get("resp_top50", {})
            max_mis_hi = float("nan")
            mis_detail: Dict[str, float] = {}
            if mismatch_preds:
                hi_list = []
                for m_name, m_mat in mismatch_preds.items():
                    m_res = compute_residual_suite_r10(
                        delta_obs=delta_obs[c_idx],
                        delta_pred=m_mat[c_idx],
                        delta_transfer=delta_transfer[c_idx],
                        target_gene_idx=np.asarray(target_gene_idx, dtype=np.int64)[c_idx],
                        delta_half_a=delta_half_a[c_idx] if delta_half_a is not None else None,
                        delta_half_b=delta_half_b[c_idx] if delta_half_b is not None else None,
                        half_valid=half_valid[c_idx] if half_valid is not None else None,
                        top_k=top_k,
                        n_boot=n_boot,
                        seed=seed,
                    )
                    hi_val = float(m_res.get("resp_top50", {}).get("r_orth_cal", {}).get("ci975_hi", float("nan")))
                    mis_detail[m_name] = hi_val
                    if np.isfinite(hi_val):
                        hi_list.append(hi_val)
                if hi_list:
                    max_mis_hi = float(max(hi_list))

            ceil_raw = float(t50.get("ceiling_r", float("nan")))
            ceil_orth = float(t50.get("ceiling_r_orth", float("nan")))
            perm100 = compute_fct_perm100_null(
                state=state,
                target_ctrl_mean_libnorm=target_ctrl_mean_libnorm,
                query_target_names=query_target_names,
                delta_obs=delta_obs,
                delta_transfer=delta_transfer,
                target_gene_idx=target_gene_idx,
                ceil_r_top50=ceil_raw,
                ceil_r_orth_top50=ceil_orth,
                fc_bundle=fc_bundle,
                n_perm=100,
                top_k=top_k,
                row_subset_mask=mask_t,
            )
            fct_lo975 = float(t50.get("r_orth_cal", {}).get("ci975_lo", float("nan")))
            perm_p975 = float(perm100.get("r_orth_cal", {}).get("p97_5", float("nan")))
            beats_mis = bool(np.isfinite(fct_lo975) and np.isfinite(max_mis_hi) and fct_lo975 > max_mis_hi)
            beats_perm = bool(np.isfinite(fct_lo975) and np.isfinite(perm_p975) and fct_lo975 > perm_p975)
            tertiles_out[lbl] = {
                "n_targets": int(c_idx.size),
                "n_cells_min": float(np.min(nc_arr[c_idx])),
                "n_cells_median": float(np.median(nc_arr[c_idx])),
                "n_cells_max": float(np.max(nc_arr[c_idx])),
                "fct_r_orth_cal": t50.get("r_orth_cal", {}),
                "fct_slope_ols": t50.get("slope_ols", {}),
                "rho_half_ab_orth": float(t50.get("rho_half_ab_orth", float("nan"))),
                "ceiling_r_orth": ceil_orth,
                "max_mismatched_ci975_hi": max_mis_hi,
                "perm100_p975": perm_p975,
                "perm100_median": float(perm100.get("r_orth_cal", {}).get("median", float("nan"))),
                "beats_mismatched_ci975": beats_mis,
                "beats_perm100_p975": beats_perm,
                "passes_hc1": bool(np.isfinite(fct_lo975) and fct_lo975 > 0.0),
                "passes_hc2": bool(beats_mis and beats_perm),
                "mismatched_ci975_hi_by_source": mis_detail,
            }

    return {
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
