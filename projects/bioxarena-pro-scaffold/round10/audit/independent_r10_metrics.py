#!/usr/bin/env python3
"""
Independent Round 10 Evaluation Metrics, Joint Bootstrap, 100-Perm Null, H-C3, R19, & S7 (A-T3).

Written independently from PREREG_C §4 & §5:
  - Joint perturbation bootstrap (`n_boot = 4,000`, `np.random.default_rng(20261122)`):
    recomputes BOTH r_orth_raw and orthogonalized split-half Spearman-Brown noise ceiling jointly.
  - Bonferroni two-sided `97.5%` CIs (`[1.25%, 98.75%]`).
  - `NegCtrl_FCT_Perm(k)` (`k = 0..99`, `np.random.default_rng(k).permutation(G)`).
  - `H-C3`: 6 paired aggregate conventions (`FCT - Transfer_Mean`) + `pds_cosine_panel`.
  - `SA-6` / `S7` (`P6`): per-pert `r_orth_p`, Spearman CIs (`ln(n_cells)`, `|k_p|`),
    mean-centered OLS `[1, ln(n_cells), |k_p|]`, and `n_cells` tertiles (`T1`, `T2`, `T3`).
  - `K_conf < 2` (`P1`): strict no-adjudication guard.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.stats import spearmanr

from independent_r10_predictors import fct_from_fc_bar

EPS = 1e-12


@dataclass
class BootstrapSummary:
    estimate: float
    ci95_lo: float
    ci95_hi: float
    ci975_lo: float
    ci975_hi: float
    n_items: int
    n_boot: int

    def to_dict(self, prefix: str = "") -> Dict[str, Any]:
        p = f"{prefix}_" if prefix else ""
        return {
            f"{p}est": self.estimate,
            f"{p}ci95_lo": self.ci95_lo,
            f"{p}ci95_hi": self.ci95_hi,
            f"{p}ci975_lo": self.ci975_lo,
            f"{p}ci975_hi": self.ci975_hi,
        }


def summarize_boot(point: float, boot_arr: np.ndarray, n_items: int, n_boot: int) -> BootstrapSummary:
    fin = boot_arr[np.isfinite(boot_arr)]
    if fin.size == 0 or not np.isfinite(point):
        return BootstrapSummary(float(point), float("nan"), float("nan"), float("nan"), float("nan"), n_items, n_boot)
    lo95, hi95 = np.percentile(fin, [2.5, 97.5])
    lo975, hi975 = np.percentile(fin, [1.25, 98.75])
    return BootstrapSummary(float(point), float(lo95), float(hi95), float(lo975), float(hi975), n_items, n_boot)


def pearson_vec(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.size < 2:
        return float("nan")
    ac = a - np.mean(a)
    bc = b - np.mean(b)
    denom = float(np.linalg.norm(ac) * np.linalg.norm(bc))
    if denom <= EPS or not np.isfinite(denom):
        return float("nan")
    return float(np.clip(np.dot(ac, bc) / denom, -1.0, 1.0))


def spearman_safe(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.size < 3 or float(np.std(a)) <= EPS or float(np.std(b)) <= EPS:
        return float("nan")
    res = spearmanr(a, b)
    val = float(res.statistic) if hasattr(res, "statistic") else float(res[0])
    return val if np.isfinite(val) else float("nan")


def ceiling_from_halves(r_ab: float) -> float:
    if not np.isfinite(r_ab) or r_ab <= 0.0:
        return 0.0
    sb = (2.0 * r_ab) / (1.0 + r_ab)
    return float(np.sqrt(max(sb, 0.0)))


def disattenuate(r_raw: float, ceil: float) -> float:
    if not np.isfinite(r_raw) or not np.isfinite(ceil) or ceil <= EPS:
        return float("nan")
    return float(np.clip(r_raw / ceil, -1.0, 1.0))


def orthogonalize_to_ref(vec: np.ndarray, ref: np.ndarray) -> np.ndarray:
    denom = float(np.dot(ref, ref))
    if denom <= EPS:
        return vec.copy()
    return vec - (float(np.dot(vec, ref)) / denom) * ref


def _extract_pert_gene_indices(
    delta_obs: np.ndarray,
    t_mean: np.ndarray,
    pred: np.ndarray,
    target_cols: Sequence[int],
    mode: str = "resp_top50",
    top_k: int = 50,
    gene_keep_mask: Optional[np.ndarray] = None,
    row_subset: Optional[Sequence[int]] = None,
) -> List[Tuple[int, np.ndarray]]:
    n_t, n_g = delta_obs.shape
    rows = range(n_t) if row_subset is None else [int(r) for r in row_subset]
    kept: List[Tuple[int, np.ndarray]] = []
    for p in rows:
        valid = np.isfinite(delta_obs[p]) & np.isfinite(t_mean[p]) & np.isfinite(pred[p])
        if gene_keep_mask is not None:
            valid &= gene_keep_mask
        tc = int(target_cols[p])
        if 0 <= tc < n_g:
            valid[tc] = False
        idx = np.flatnonzero(valid)
        if idx.size == 0:
            continue
        if mode == "resp_top50" and idx.size > top_k:
            idx = idx[np.argsort(np.abs(delta_obs[p, idx]))[-top_k:]]
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


def _moments_stats(
    nv: np.ndarray,
    sx: np.ndarray,
    sy: np.ndarray,
    sxx: np.ndarray,
    syy: np.ndarray,
    sxy: np.ndarray,
    draws: Optional[np.ndarray] = None,
) -> Tuple[Any, ...]:
    if draws is None:
        N = float(np.sum(nv))
        Sx, Sy = float(np.sum(sx)), float(np.sum(sy))
        Sxx, Syy, Sxy = float(np.sum(sxx)), float(np.sum(syy)), float(np.sum(sxy))
    else:
        N = np.sum(nv[draws], axis=1)
        Sx = np.sum(sx[draws], axis=1)
        Sy = np.sum(sy[draws], axis=1)
        Sxx = np.sum(sxx[draws], axis=1)
        Syy = np.sum(syy[draws], axis=1)
        Sxy = np.sum(sxy[draws], axis=1)

    N_s = np.maximum(N, 1.0)
    cov = Sxy - (Sx * Sy) / N_s
    vx = Sxx - (Sx * Sx) / N_s
    vy = Syy - (Sy * Sy) / N_s
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = np.where(np.sqrt(np.maximum(vx, 0.0) * np.maximum(vy, 0.0)) > EPS,
                        np.clip(cov / np.sqrt(np.maximum(vx, 0.0) * np.maximum(vy, 0.0)), -1.0, 1.0),
                        0.0)
        slope_ols = np.where(vx > EPS, cov / vx, 0.0)
        sd_ratio = np.where(vx > EPS, np.sqrt(np.maximum(vy, 0.0) / np.maximum(vx, EPS)), 0.0)
        diff = vy - vx
        slope_deming = np.where(
            np.abs(cov) > EPS,
            (diff + np.sqrt(diff * diff + 4.0 * cov * cov)) / (2.0 * cov),
            0.0,
        )
    return corr, slope_ols, sd_ratio, slope_deming


def evaluate_joint_residuals(
    delta_obs: np.ndarray,
    pred: np.ndarray,
    t_mean: np.ndarray,
    target_cols: Sequence[int],
    mode: str = "resp_top50",
    n_boot: int = 4000,
    seed: int = 20261122,
    half_a_obs: Optional[np.ndarray] = None,
    half_b_obs: Optional[np.ndarray] = None,
    half_valid: Optional[np.ndarray] = None,
    gene_keep_mask: Optional[np.ndarray] = None,
    row_subset: Optional[Sequence[int]] = None,
) -> Dict[str, Any]:
    """
    Independently compute raw and `Transfer_Mean`-orthogonalized residual endpoints
    with joint perturbation-level bootstrap (`rng = default_rng(20261122)`).
    """
    kept = _extract_pert_gene_indices(
        delta_obs=delta_obs,
        t_mean=t_mean,
        pred=pred,
        target_cols=target_cols,
        mode=mode,
        gene_keep_mask=gene_keep_mask,
        row_subset=row_subset,
    )
    m = len(kept)
    if m == 0:
        return {"n_eval_perts": 0}

    xr_raw, yr_raw = [], []
    xr_orth, yr_orth = [], []
    ha_raw, hb_raw = [], []
    ha_orth, hb_orth = [], []
    any_half = False

    has_halves = half_a_obs is not None and half_b_obs is not None
    for p, idx in kept:
        ref = t_mean[p, idx].astype(np.float64)
        ro = delta_obs[p, idx].astype(np.float64) - ref
        rp = pred[p, idx].astype(np.float64) - ref
        xr_raw.append(ro)
        yr_raw.append(rp)
        xr_orth.append(orthogonalize_to_ref(ro, ref))
        yr_orth.append(orthogonalize_to_ref(rp, ref))

        hv_ok = has_halves and (half_valid is None or bool(half_valid[p]))
        if hv_ok:
            ha_p = half_a_obs[p, idx].astype(np.float64)
            hb_p = half_b_obs[p, idx].astype(np.float64)
            fin = np.isfinite(ha_p) & np.isfinite(hb_p)
            if np.all(fin):
                hra = ha_p - ref
                hrb = hb_p - ref
                ha_raw.append(hra)
                hb_raw.append(hrb)
                ha_orth.append(orthogonalize_to_ref(hra, ref))
                hb_orth.append(orthogonalize_to_ref(hrb, ref))
                any_half = True
                continue
        ha_raw.append(np.zeros(0, dtype=np.float64))
        hb_raw.append(np.zeros(0, dtype=np.float64))
        ha_orth.append(np.zeros(0, dtype=np.float64))
        hb_orth.append(np.zeros(0, dtype=np.float64))

    rng = np.random.default_rng(seed)
    draws = rng.integers(0, m, size=(n_boot, m))

    # Raw & orthogonalized predictor moments
    m_raw = _moments(xr_raw, yr_raw)
    m_orth = _moments(xr_orth, yr_orth)
    pt_r_raw, pt_slope_raw, pt_sd_raw, pt_dem_raw = _moments_stats(*m_raw)
    pt_r_orth, pt_slope_orth, pt_sd_orth, pt_dem_orth = _moments_stats(*m_orth)
    b_r_raw, b_slope_raw, b_sd_raw, b_dem_raw = _moments_stats(*m_raw, draws=draws)
    b_r_orth, b_slope_orth, b_sd_orth, b_dem_orth = _moments_stats(*m_orth, draws=draws)

    if any_half:
        mh_raw = _moments(ha_raw, hb_raw)
        mh_orth = _moments(ha_orth, hb_orth)
        pt_hab_raw = float(_moments_stats(*mh_raw)[0])
        pt_hab_orth = float(_moments_stats(*mh_orth)[0])
        pt_ceil_raw = ceiling_from_halves(pt_hab_raw)
        pt_ceil_orth = ceiling_from_halves(pt_hab_orth)
        b_hab_raw = _moments_stats(*mh_raw, draws=draws)[0]
        b_hab_orth = _moments_stats(*mh_orth, draws=draws)[0]
        b_ceil_raw = np.asarray([ceiling_from_halves(float(x)) for x in b_hab_raw], dtype=np.float64)
        b_ceil_orth = np.asarray([ceiling_from_halves(float(x)) for x in b_hab_orth], dtype=np.float64)
    else:
        pt_hab_raw, pt_hab_orth = float("nan"), float("nan")
        pt_ceil_raw, pt_ceil_orth = 0.0, 0.0
        b_ceil_raw = np.zeros(n_boot, dtype=np.float64)
        b_ceil_orth = np.zeros(n_boot, dtype=np.float64)

    pt_cal = disattenuate(float(pt_r_raw), pt_ceil_raw)
    b_cal = np.asarray([disattenuate(float(r), float(c)) for r, c in zip(b_r_raw, b_ceil_raw)], dtype=np.float64)

    pt_orth_cal = disattenuate(float(pt_r_orth), pt_ceil_orth)
    b_orth_cal = np.asarray(
        [disattenuate(float(r), float(c)) for r, c in zip(b_r_orth, b_ceil_orth)], dtype=np.float64
    )

    out: Dict[str, Any] = {
        "n_eval_perts": m,
        "rho_ab": pt_hab_raw,
        "noise_ceiling_r": pt_ceil_raw,
        "rho_ab_orth": pt_hab_orth,
        "noise_ceiling_r_orth": pt_ceil_orth,
    }
    out.update(summarize_boot(pt_orth_cal, b_orth_cal, m, n_boot).to_dict("r_orth_cal"))
    out.update(summarize_boot(float(pt_r_orth), b_r_orth, m, n_boot).to_dict("r_orth_raw"))
    out.update(summarize_boot(pt_cal, b_cal, m, n_boot).to_dict("r_cal"))
    out.update(summarize_boot(float(pt_slope_raw), b_slope_raw, m, n_boot).to_dict("slope_ols"))
    out.update(summarize_boot(float(pt_sd_raw), b_sd_raw, m, n_boot).to_dict("sd_pred_over_sd_obs"))
    out.update(summarize_boot(float(pt_dem_raw), b_dem_raw, m, n_boot).to_dict("slope_deming"))
    out.update(summarize_boot(float(pt_slope_orth), b_slope_orth, m, n_boot).to_dict("slope_orth_ols"))
    return out


def evaluate_paired_cal_diff(
    delta_obs: np.ndarray,
    pred_a: np.ndarray,
    pred_b: np.ndarray,
    t_mean: np.ndarray,
    target_cols: Sequence[int],
    mode: str = "resp_top50",
    n_boot: int = 4000,
    seed: int = 20261122,
    half_a_obs: Optional[np.ndarray] = None,
    half_b_obs: Optional[np.ndarray] = None,
    half_valid: Optional[np.ndarray] = None,
    ortho: bool = True,
) -> Dict[str, Any]:
    """Paired joint-bootstrap CI for `pred_a - pred_b` (R14 Hybrid_on_TMVS vs TMVS)."""
    kept_a = _extract_pert_gene_indices(delta_obs, t_mean, pred_a, target_cols, mode=mode)
    kept_b = _extract_pert_gene_indices(delta_obs, t_mean, pred_b, target_cols, mode=mode)
    if len(kept_a) == 0 or [p for p, _ in kept_a] != [p for p, _ in kept_b]:
        return summarize_boot(float("nan"), np.array([float("nan")]), 0, n_boot).to_dict("diff")

    m = len(kept_a)
    xa, ya, xb, yb, ha_l, hb_l = [], [], [], [], [], []
    any_half = False
    has_halves = half_a_obs is not None and half_b_obs is not None
    for (p, idx) in kept_a:
        ref = t_mean[p, idx].astype(np.float64)
        ro = delta_obs[p, idx].astype(np.float64) - ref
        rpa = pred_a[p, idx].astype(np.float64) - ref
        rpb = pred_b[p, idx].astype(np.float64) - ref
        if ortho:
            xa.append(orthogonalize_to_ref(ro, ref))
            ya.append(orthogonalize_to_ref(rpa, ref))
            xb.append(orthogonalize_to_ref(ro, ref))
            yb.append(orthogonalize_to_ref(rpb, ref))
        else:
            xa.append(ro); ya.append(rpa); xb.append(ro); yb.append(rpb)
        if has_halves and (half_valid is None or bool(half_valid[p])):
            hap = half_a_obs[p, idx].astype(np.float64)
            hbp = half_b_obs[p, idx].astype(np.float64)
            if np.all(np.isfinite(hap) & np.isfinite(hbp)):
                hra, hrb = hap - ref, hbp - ref
                if ortho:
                    hra = orthogonalize_to_ref(hra, ref)
                    hrb = orthogonalize_to_ref(hrb, ref)
                ha_l.append(hra); hb_l.append(hrb)
                any_half = True
                continue
        ha_l.append(np.zeros(0, dtype=np.float64))
        hb_l.append(np.zeros(0, dtype=np.float64))

    rng = np.random.default_rng(seed)
    draws = rng.integers(0, m, size=(n_boot, m))
    ma, mb = _moments(xa, ya), _moments(xb, yb)
    pt_ra = float(_moments_stats(*ma)[0])
    pt_rb = float(_moments_stats(*mb)[0])
    b_ra = _moments_stats(*ma, draws=draws)[0]
    b_rb = _moments_stats(*mb, draws=draws)[0]

    if any_half:
        mh = _moments(ha_l, hb_l)
        pt_ceil = ceiling_from_halves(float(_moments_stats(*mh)[0]))
        b_ceil = np.asarray([ceiling_from_halves(float(v)) for v in _moments_stats(*mh, draws=draws)[0]], dtype=np.float64)
    else:
        pt_ceil = 0.0
        b_ceil = np.zeros(n_boot, dtype=np.float64)

    pt_diff = disattenuate(pt_ra, pt_ceil) - disattenuate(pt_rb, pt_ceil)
    b_diff = np.asarray(
        [disattenuate(float(ra), float(cb)) - disattenuate(float(rb), float(cb)) for ra, rb, cb in zip(b_ra, b_rb, b_ceil)],
        dtype=np.float64,
    )
    return summarize_boot(pt_diff, b_diff, m, n_boot).to_dict("diff")


def evaluate_fct_perm_null(
    delta_obs: np.ndarray,
    t_mean: np.ndarray,
    fc_bar: np.ndarray,
    tgt_ctrl_mean_libnorm: np.ndarray,
    target_cols: Sequence[int],
    ceil_raw: float,
    ceil_orth: float,
    mode: str = "resp_top50",
    n_perm: int = 100,
    row_subset: Optional[Sequence[int]] = None,
) -> Dict[str, Any]:
    """
    PREREG_C §3.2 & §6: `NegCtrl_FCT_Perm(k)` (`k = 0..99`).
    Permutes target control mean `mu_t` across `G` genes using `np.random.default_rng(k).permutation(G)`
    while keeping `FC_bar` unpermuted.
    """
    n_g = delta_obs.shape[1]
    kept = _extract_pert_gene_indices(delta_obs, t_mean, fc_bar, target_cols, mode=mode, row_subset=row_subset)
    if not kept:
        nan_pct = {"p025": float("nan"), "median": float("nan"), "p975": float("nan")}
        return {"n_perm": n_perm, "r_orth_cal": nan_pct, "r_cal": nan_pct, "slope_ols": nan_pct}

    # Pre-cache r_obs_raw, r_obs_orth, t_mean slices, fc_bar slices, and orthogonalization projectors
    ro_raw_list, ro_orth_list, tm_list, fc_list, unit_ref_list = [], [], [], [], []
    for p, idx in kept:
        ref = t_mean[p, idx].astype(np.float64)
        ro = delta_obs[p, idx].astype(np.float64) - ref
        nrm = float(np.linalg.norm(ref))
        u_ref = (ref / nrm) if nrm > EPS else np.zeros_like(ref)
        ro_raw_list.append(ro)
        ro_orth_list.append(ro - float(np.dot(ro, u_ref)) * u_ref)
        tm_list.append(ref)
        fc_list.append(fc_bar[p, idx].astype(np.float64))
        unit_ref_list.append((p, idx, u_ref))

    x_raw = np.concatenate(ro_raw_list)
    x_orth = np.concatenate(ro_orth_list)
    mu_safe = np.maximum(np.asarray(tgt_ctrl_mean_libnorm, dtype=np.float64), 0.0)

    r_orth_cal_vals = np.empty(n_perm, dtype=np.float64)
    r_cal_vals = np.empty(n_perm, dtype=np.float64)
    slope_vals = np.empty(n_perm, dtype=np.float64)

    for k in range(n_perm):
        perm = np.random.default_rng(k).permutation(n_g)
        mu_p = mu_safe[perm]
        log1p_mu_p = np.log1p(mu_p)
        yp_raw_l, yp_orth_l = [], []
        for idx_i, (p, col_idx, u_ref) in enumerate(unit_ref_list):
            mp = mu_p[col_idx]
            pred_slice = np.log1p(mp * fc_list[idx_i]) - log1p_mu_p[col_idx]
            rp = pred_slice - tm_list[idx_i]
            yp_raw_l.append(rp)
            yp_orth_l.append(rp - float(np.dot(rp, u_ref)) * u_ref)
        y_raw = np.concatenate(yp_raw_l)
        y_orth = np.concatenate(yp_orth_l)

        xc = x_raw - np.mean(x_raw)
        yc = y_raw - np.mean(y_raw)
        vx = float(np.dot(xc, xc))
        slope_vals[k] = float(np.dot(xc, yc) / vx) if vx > EPS else 0.0
        r_cal_vals[k] = disattenuate(pearson_vec(x_raw, y_raw), ceil_raw)
        r_orth_cal_vals[k] = disattenuate(pearson_vec(x_orth, y_orth), ceil_orth)

    def _pcts(arr: np.ndarray) -> Dict[str, float]:
        fin = arr[np.isfinite(arr)]
        if fin.size == 0:
            return {"p025": float("nan"), "median": float("nan"), "p975": float("nan")}
        p025, p50, p975 = np.percentile(fin, [2.5, 50.0, 97.5])
        return {"p025": float(p025), "median": float(p50), "p975": float(p975)}

    return {
        "n_perm": n_perm,
        "r_orth_cal": _pcts(r_orth_cal_vals),
        "r_cal": _pcts(r_cal_vals),
        "slope_ols": _pcts(slope_vals),
    }


def _row_mask(n_g: int, cols: Sequence[int] | None) -> np.ndarray:
    m = np.ones(n_g, dtype=bool)
    if cols is not None:
        for c in cols:
            if c is not None and 0 <= int(c) < n_g:
                m[int(c)] = False
    return m


def _midrank_disc(sim: np.ndarray, p_zero: np.ndarray) -> np.ndarray:
    n_p = sim.shape[0]
    if n_p <= 1:
        return np.full(n_p, 0.5, dtype=np.float64)
    diag = np.diag(sim)[:, None]
    better = np.sum(sim > diag + 1e-12, axis=1)
    ties = np.sum(np.abs(sim - diag) <= 1e-12, axis=1) - 1
    out = 1.0 - (better + 0.5 * np.maximum(ties, 0)) / float(n_p - 1)
    out[p_zero] = 0.5
    return out


def compute_disc_panel_rows(
    true_mat: np.ndarray,
    pred_mat: np.ndarray,
    masked_cols: Sequence[int] | None,
    center: bool = True,
) -> np.ndarray:
    n_p, n_g = true_mat.shape
    mask = _row_mask(n_g, masked_cols)
    t_s = np.asarray(true_mat[:, mask], dtype=np.float64)
    p_s = np.asarray(pred_mat[:, mask], dtype=np.float64)
    if center:
        t_s = t_s - np.mean(t_s, axis=1, keepdims=True)
        p_s = p_s - np.mean(p_s, axis=1, keepdims=True)
    tn = np.linalg.norm(t_s, axis=1, keepdims=True)
    pn = np.linalg.norm(p_s, axis=1, keepdims=True)
    sim = (p_s / np.maximum(pn, EPS)) @ (t_s / np.maximum(tn, EPS)).T
    pz, tz = pn[:, 0] <= 1e-10, tn[:, 0] <= 1e-10
    sim[pz, :] = 0.0
    sim[:, tz] = 0.0
    return _midrank_disc(sim, pz)


def compute_disc_ownmasked_rows(
    true_mat: np.ndarray,
    pred_mat: np.ndarray,
    target_cols: Sequence[int],
    center: bool = True,
) -> np.ndarray:
    n_p, n_g = true_mat.shape
    t_m = np.asarray(true_mat, dtype=np.float64).copy()
    p_m = np.asarray(pred_mat, dtype=np.float64).copy()
    for i, tc in enumerate(target_cols):
        if tc is not None and 0 <= int(tc) < n_g:
            t_m[i, int(tc)] = 0.0
            p_m[i, int(tc)] = 0.0
    if center:
        t_m = t_m - np.mean(t_m, axis=1, keepdims=True)
        p_m = p_m - np.mean(p_m, axis=1, keepdims=True)
    tn = np.linalg.norm(t_m, axis=1, keepdims=True)
    pn = np.linalg.norm(p_m, axis=1, keepdims=True)
    sim = (p_m / np.maximum(pn, EPS)) @ (t_m / np.maximum(tn, EPS)).T
    pz, tz = pn[:, 0] <= 1e-10, tn[:, 0] <= 1e-10
    sim[pz, :] = 0.0
    sim[:, tz] = 0.0
    return _midrank_disc(sim, pz)


def compute_de50_r_rows(
    true_mat: np.ndarray,
    pred_mat: np.ndarray,
    target_cols: Sequence[int],
    top_k: int = 50,
) -> np.ndarray:
    n_p, n_g = true_mat.shape
    out = np.empty(n_p, dtype=np.float64)
    for i in range(n_p):
        mask = _row_mask(n_g, [target_cols[i]] if target_cols[i] >= 0 else None)
        dt = true_mat[i, mask]
        dp = pred_mat[i, mask]
        k = min(top_k, dt.size)
        if k < 2:
            out[i] = 0.0
            continue
        top = np.argpartition(np.abs(dt), -k)[-k:]
        val = pearson_vec(dt[top], dp[top])
        out[i] = val if np.isfinite(val) else 0.0
    return out


def evaluate_hc3_noninferiority(
    delta_obs: np.ndarray,
    pred_fct: np.ndarray,
    pred_tmean: np.ndarray,
    grand_mean: np.ndarray,
    target_cols: Sequence[int],
    n_boot: int = 4000,
    seed: int = 20261122,
) -> Dict[str, Any]:
    """PREREG_C §4.2 `H-C3`: 6 paired conventions (`FCT - Transfer_Mean`) + descriptive `pds_cosine_panel`."""
    specs = (
        ("ctrl_disc_panel", "ctrl", "panel"),
        ("ctrl_disc_ownmasked", "ctrl", "ownmasked"),
        ("ctrl_de50_r", "ctrl", "de50_r"),
        ("cent_disc_panel", "cent", "panel"),
        ("cent_disc_ownmasked", "cent", "ownmasked"),
        ("cent_de50_r", "cent", "de50_r"),
    )
    rng = np.random.default_rng(seed)
    n_p = delta_obs.shape[0]
    draws = rng.integers(0, n_p, size=(n_boot, n_p)) if n_p > 0 else np.zeros((0, 0), dtype=np.int64)

    true_cent = delta_obs - grand_mean[None, :]
    fct_cent = pred_fct - grand_mean[None, :]
    tm_cent = pred_tmean - grand_mean[None, :]

    results: Dict[str, Any] = {}
    n_inferior = 0
    n_superior = 0
    for name, space_mode, agg_mode in specs:
        t_mat = delta_obs if space_mode == "ctrl" else true_cent
        f_mat = pred_fct if space_mode == "ctrl" else fct_cent
        b_mat = pred_tmean if space_mode == "ctrl" else tm_cent
        if agg_mode == "panel":
            rf = compute_disc_panel_rows(t_mat, f_mat, target_cols, center=True)
            rb = compute_disc_panel_rows(t_mat, b_mat, target_cols, center=True)
        elif agg_mode == "ownmasked":
            rf = compute_disc_ownmasked_rows(t_mat, f_mat, target_cols, center=True)
            rb = compute_disc_ownmasked_rows(t_mat, b_mat, target_cols, center=True)
        else:
            rf = compute_de50_r_rows(t_mat, f_mat, target_cols)
            rb = compute_de50_r_rows(t_mat, b_mat, target_cols)

        diffs = rf - rb
        pt_diff = float(np.mean(diffs))
        b_diff = np.mean(diffs[draws], axis=1) if n_p > 0 else np.array([pt_diff])
        ci = summarize_boot(pt_diff, b_diff, n_p, n_boot)
        non_inf = bool(ci.ci975_hi >= 0.0) if np.isfinite(ci.ci975_hi) else False
        sup = bool(ci.ci975_lo > 0.0) if np.isfinite(ci.ci975_lo) else False
        if not non_inf:
            n_inferior += 1
        if sup:
            n_superior += 1
        results[name] = {
            "fct_mean": float(np.mean(rf)),
            "tmean_mean": float(np.mean(rb)),
            "diff_mean": ci.estimate,
            "diff_ci975_lo": ci.ci975_lo,
            "diff_ci975_hi": ci.ci975_hi,
            "non_inferior": non_inf,
            "strictly_superior": sup,
        }

    # Companion uncentered cosine diagnostic (`pds_cosine_panel`)
    cos_f = compute_disc_panel_rows(delta_obs, pred_fct, target_cols, center=False)
    cos_b = compute_disc_panel_rows(delta_obs, pred_tmean, target_cols, center=False)
    cos_diff = cos_f - cos_b
    b_cos = np.mean(cos_diff[draws], axis=1) if n_p > 0 else np.array([float("nan")])
    ci_cos = summarize_boot(float(np.mean(cos_diff)), b_cos, n_p, n_boot)
    results["pds_cosine_panel"] = {
        "fct_mean": float(np.mean(cos_f)),
        "tmean_mean": float(np.mean(cos_b)),
        "diff_mean": ci_cos.estimate,
        "diff_ci975_lo": ci_cos.ci975_lo,
        "diff_ci975_hi": ci_cos.ci975_hi,
    }
    results["hc3_pass_line"] = bool(n_inferior == 0 and n_p > 0)
    results["n_inferior_conventions"] = int(n_inferior)
    results["n_superior_conventions"] = int(n_superior)
    return results


def evaluate_r19_gain_decomposition(
    delta_obs: np.ndarray,
    t_mean: np.ndarray,
    target_cols: Sequence[int],
    mode: str = "resp_top50",
    n_boot: int = 4000,
    seed: int = 20261122,
    min_perts_per_gene: int = 5,
) -> Dict[str, Any]:
    """R19 three-tier multiplicative SST(r_obs) decomposition (L1 whole-line beta_t, L2 gene m_tg, L3 rem)."""
    n_t, n_g = delta_obs.shape
    rows_p, xs_l, ys_l, g_idx_l = [], [], [], []
    for p in range(n_t):
        valid = np.isfinite(delta_obs[p]) & np.isfinite(t_mean[p])
        tc = int(target_cols[p])
        if 0 <= tc < n_g:
            valid[tc] = False
        idx = np.flatnonzero(valid)
        if idx.size == 0:
            continue
        if mode == "resp_top50" and idx.size > 50:
            idx = idx[np.argsort(np.abs(delta_obs[p, idx]))[-50:]]
        rows_p.append(p)
        xs_l.append(t_mean[p, idx].astype(np.float64))
        ys_l.append(delta_obs[p, idx].astype(np.float64))
        g_idx_l.append(idx)

    P = len(rows_p)
    if P == 0:
        return {}

    m_N = np.asarray([float(x.size) for x in xs_l], dtype=np.float64)
    m_S = np.asarray([float(np.sum(x)) for x in xs_l], dtype=np.float64)
    m_Y = np.asarray([float(np.sum(y)) for y in ys_l], dtype=np.float64)
    m_SS = np.asarray([float(np.dot(x, x)) for x in xs_l], dtype=np.float64)
    m_SY = np.asarray([float(np.dot(x, y)) for x, y in zip(xs_l, ys_l)], dtype=np.float64)
    m_YY = np.asarray([float(np.dot(y, y)) for y in ys_l], dtype=np.float64)

    m_SY_r = m_SY - m_SS
    m_YY_r = m_YY - 2.0 * m_SY + m_SS
    m_Y_r = m_Y - m_S

    # Sparse Pert x Gene incidence accumulators
    import scipy.sparse as sp
    flat_rows = np.concatenate([np.full(g.size, i, dtype=np.int64) for i, g in enumerate(g_idx_l)])
    flat_cols = np.concatenate(g_idx_l)
    flat_x = np.concatenate(xs_l)
    flat_yr = np.concatenate([y - x for x, y in zip(xs_l, ys_l)])

    ones_m = sp.csr_matrix((np.ones(flat_rows.size), (flat_rows, flat_cols)), shape=(P, n_g))
    s_m = sp.csr_matrix((flat_x, (flat_rows, flat_cols)), shape=(P, n_g))
    y_m = sp.csr_matrix((flat_yr, (flat_rows, flat_cols)), shape=(P, n_g))
    ss_m = sp.csr_matrix((flat_x * flat_x, (flat_rows, flat_cols)), shape=(P, n_g))
    sy_m = sp.csr_matrix((flat_x * flat_yr, (flat_rows, flat_cols)), shape=(P, n_g))
    s2_m = ss_m
    s2y_m = sp.csr_matrix(((flat_x ** 2) * flat_yr, (flat_rows, flat_cols)), shape=(P, n_g))
    s2y2_m = sp.csr_matrix(((flat_x ** 2) * (flat_yr ** 2), (flat_rows, flat_cols)), shape=(P, n_g))
    s3_m = sp.csr_matrix((flat_x ** 3, (flat_rows, flat_cols)), shape=(P, n_g))
    s3y_m = sp.csr_matrix(((flat_x ** 3) * flat_yr, (flat_rows, flat_cols)), shape=(P, n_g))
    s4_m = sp.csr_matrix((flat_x ** 4, (flat_rows, flat_cols)), shape=(P, n_g))

    def _solve(W: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        Nt = np.maximum(W @ m_N, 1.0)
        St, Yt = W @ m_S, W @ m_Y_r
        SSt, SYt, YYt = W @ m_SS, W @ m_SY_r, W @ m_YY_r
        sst_r = np.maximum(YYt - (Yt * Yt) / Nt, EPS)
        cov_syr = SYt - (St * Yt) / Nt
        var_sr = np.maximum(SSt - (St * St) / Nt, EPS)
        beta = cov_syr / var_sr
        alpha = (Yt - beta * St) / Nt
        sse1 = np.maximum(
            YYt - 2.0 * alpha * Yt - 2.0 * beta * SYt + (alpha * alpha) * Nt + 2.0 * alpha * beta * St + (beta * beta) * SSt,
            0.0,
        )
        f1 = np.clip(1.0 - sse1 / sst_r, 0.0, 1.0)

        W_bin = (W > 0).astype(np.float64)
        c_dist = W_bin @ ones_m
        s_g = W @ s_m
        y_g = W @ y_m
        ss_g = W @ ss_m
        sy_g = W @ sy_m
        ac, bc = alpha[:, None], beta[:, None]
        ok = (c_dist >= min_perts_per_gene) & (ss_g > EPS)
        mg = np.where(ok, (sy_g - ac * s_g) / np.maximum(ss_g, EPS), bc)
        imp_raw = np.where(ok, ((mg - bc) ** 2) * ss_g, 0.0)

        W2 = W * W
        s2 = W2 @ s2_m
        s2y = W2 @ s2y_m
        s2y2 = W2 @ s2y2_m
        s3 = W2 @ s3_m
        s3y = W2 @ s3y_m
        s4 = W2 @ s4_m
        num_w2 = s2y2 - 2.0 * ac * s2y + (ac * ac) * s2 - 2.0 * mg * (s3y - ac * s3) + (mg * mg) * s4
        neff = np.where(s4 > EPS, (ss_g * ss_g) / np.maximum(s4, EPS), c_dist)
        df_adj = np.maximum(1.0 - 1.0 / np.maximum(neff, 2.0), 0.25)
        bias = np.where(ok, np.maximum(num_w2, 0.0) / (np.maximum(ss_g, EPS) * df_adj), 0.0)
        imp_unb = np.maximum(np.sum(imp_raw - bias, axis=1), 0.0)
        sse2 = np.maximum(sse1 - imp_unb, 0.0)
        f12 = np.clip(1.0 - sse2 / sst_r, f1, 1.0)
        f2 = np.clip(f12 - f1, 0.0, 1.0)
        f3 = np.clip(1.0 - f12, 0.0, 1.0)
        return beta, f1, f2, f3

    pt_b, pt_f1, pt_f2, pt_f3 = _solve(np.ones((1, P), dtype=np.float64))
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, P, size=(n_boot, P))
    bb, bf1, bf2, bf3 = [np.empty(n_boot, dtype=np.float64) for _ in range(4)]
    for s0 in range(0, n_boot, 500):
        s1 = min(s0 + 500, n_boot)
        W = np.zeros((s1 - s0, P), dtype=np.float64)
        for bi in range(s1 - s0):
            W[bi] = np.bincount(draws[s0 + bi], minlength=P)
        rb, rf1, rf2, rf3 = _solve(W)
        bb[s0:s1], bf1[s0:s1], bf2[s0:s1], bf3[s0:s1] = rb, rf1, rf2, rf3

    out = {"n_eval_perts": P}
    out.update(summarize_boot(float(pt_b[0]), bb, P, n_boot).to_dict("beta_t"))
    out.update(summarize_boot(float(pt_f1[0]), bf1, P, n_boot).to_dict("frac_whole_line"))
    out.update(summarize_boot(float(pt_f2[0]), bf2, P, n_boot).to_dict("frac_per_gene_inc"))
    out.update(summarize_boot(float(pt_f3[0]), bf3, P, n_boot).to_dict("frac_remaining"))
    return out


def _fit_ols(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    coef, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    return np.asarray(coef, dtype=np.float64)


def evaluate_s7_confounder_controls(
    eval_targets: Sequence[str],
    delta_obs: np.ndarray,
    pred_fct: np.ndarray,
    t_mean: np.ndarray,
    fc_bar: np.ndarray,
    tgt_ctrl_mean_libnorm: np.ndarray,
    source_names: Sequence[str],
    source_ctrl_mean_libnorm: Dict[str, np.ndarray],
    target_cols: Sequence[int],
    n_cells: Optional[np.ndarray],
    kd_depth_log1p: np.ndarray,
    half_a_obs: Optional[np.ndarray] = None,
    half_b_obs: Optional[np.ndarray] = None,
    half_valid: Optional[np.ndarray] = None,
    n_boot: int = 4000,
    n_perm: int = 100,
    seed: int = 20261122,
) -> Dict[str, Any]:
    """
    Independently implement PREREG_C §5 (`SA-6` / `S7` / `P6`) exact 4-part confounder controls:
      (i)   Per-perturbation orthogonalized Pearson `r_orth_p` (50 genes)
      (ii)  Spearman(r_orth_p, ln(n_cells)) and Spearman(r_orth_p, |k_p|) + 97.5% bootstrap CIs
      (iii) Mean-centered OLS r_orth_p ~ [1, ln(n_cells) - mean, |k_p| - mean] + 97.5% CIs
      (iv)  Balanced `n_cells` tertiles T1/T2/T3 recomputing pooled FCT, 6 mismatched means, & 100 perms.
    """
    n_t = len(eval_targets)
    r_orth_p = np.full(n_t, np.nan, dtype=np.float64)
    kept_list = _extract_pert_gene_indices(delta_obs, t_mean, pred_fct, target_cols, mode="resp_top50")
    for p, idx in kept_list:
        ref = t_mean[p, idx].astype(np.float64)
        ro = orthogonalize_to_ref(delta_obs[p, idx].astype(np.float64) - ref, ref)
        rp = orthogonalize_to_ref(pred_fct[p, idx].astype(np.float64) - ref, ref)
        r_orth_p[p] = pearson_vec(ro, rp)

    ln_nc = (
        np.where(np.isfinite(n_cells) & (n_cells > 0), np.log(n_cells), np.nan)
        if n_cells is not None
        else np.full(n_t, np.nan, dtype=np.float64)
    )
    abs_kp = np.asarray(kd_depth_log1p, dtype=np.float64)

    # (ii) Spearman correlations + 97.5% bootstrap CIs
    rng = np.random.default_rng(seed)

    def _boot_spearman(x: np.ndarray, y: np.ndarray) -> BootstrapSummary:
        valid = np.isfinite(x) & np.isfinite(y)
        xv, yv = x[valid], y[valid]
        k = int(xv.size)
        if k < 5:
            return BootstrapSummary(float("nan"), float("nan"), float("nan"), float("nan"), float("nan"), k, n_boot)
        pt = spearman_safe(xv, yv)
        draws = rng.integers(0, k, size=(n_boot, k))
        vals = np.asarray([spearman_safe(xv[d], yv[d]) for d in draws], dtype=np.float64)
        return summarize_boot(pt, vals, k, n_boot)

    sp_nc = _boot_spearman(r_orth_p, ln_nc)
    sp_kp = _boot_spearman(r_orth_p, abs_kp)

    # (iii) Mean-centered OLS r_orth_p ~ [1, ln(n_cells) - c1, |k_p| - c2]
    cc_mask = np.isfinite(r_orth_p) & np.isfinite(ln_nc) & np.isfinite(abs_kp)
    n_cc = int(np.sum(cc_mask))
    ols_out: Dict[str, Any] = {"n_complete_cases": n_cc}
    if n_cc >= 6:
        y_cc = r_orth_p[cc_mask]
        x1_cc = ln_nc[cc_mask]
        x2_cc = abs_kp[cc_mask]
        c1 = float(np.mean(x1_cc))
        c2 = float(np.mean(x2_cc))
        Xc = np.column_stack([np.ones(n_cc, dtype=np.float64), x1_cc - c1, x2_cc - c2])
        Xu = np.column_stack([np.ones(n_cc, dtype=np.float64), x1_cc, x2_cc])
        pt_coef = _fit_ols(Xc, y_cc)
        pt_uncent = _fit_ols(Xu, y_cc)
        inv_err = abs(float(pt_coef[0]) - float(pt_uncent[0] + pt_uncent[1] * c1 + pt_uncent[2] * c2))

        y_hat = Xc @ pt_coef
        sst = float(np.sum((y_cc - np.mean(y_cc)) ** 2))
        sse = float(np.sum((y_cc - y_hat) ** 2))
        r2 = float(1.0 - sse / sst) if sst > EPS else 0.0

        draws_cc = rng.integers(0, n_cc, size=(n_boot, n_cc))
        b_int_fixed = np.empty(n_boot, dtype=np.float64)
        b_int_recent = np.empty(n_boot, dtype=np.float64)
        b_s1 = np.empty(n_boot, dtype=np.float64)
        b_s2 = np.empty(n_boot, dtype=np.float64)
        for bi, d in enumerate(draws_cc):
            coef_b = _fit_ols(Xc[d], y_cc[d])
            b_int_fixed[bi] = coef_b[0]
            b_int_recent[bi] = float(np.mean(y_cc[d]))
            b_s1[bi] = coef_b[1]
            b_s2[bi] = coef_b[2]

        ols_out["mean_ln_ncells"] = c1
        ols_out["mean_abs_kp"] = c2
        ols_out["uncentered_eval_invariance_err"] = inv_err
        ols_out["r_squared"] = r2
        ols_out.update(summarize_boot(float(pt_coef[0]), b_int_fixed, n_cc, n_boot).to_dict("intercept_adjusted"))
        ols_out.update(summarize_boot(float(pt_coef[0]), b_int_recent, n_cc, n_boot).to_dict("intercept_recentered"))
        ols_out.update(summarize_boot(float(pt_coef[1]), b_s1, n_cc, n_boot).to_dict("slope_ln_ncells"))
        ols_out.update(summarize_boot(float(pt_coef[2]), b_s2, n_cc, n_boot).to_dict("slope_abs_kp"))
    else:
        ols_out["status"] = "INSUFFICIENT_COMPLETE_CASES"

    # (iv) Balanced n_cells tertiles T1, T2, T3
    tertiles_out: Dict[str, Any] = {}
    if n_cells is not None and int(np.sum(np.isfinite(n_cells) & (n_cells > 0))) >= 9:
        valid_nc_rows = [i for i in range(n_t) if np.isfinite(n_cells[i]) and n_cells[i] > 0]
        valid_nc_rows.sort(key=lambda i: (float(n_cells[i]), str(eval_targets[i])))
        splits = np.array_split(np.asarray(valid_nc_rows, dtype=np.int64), 3)
        for label, sub in zip(("T1_low_ncells", "T2_mid_ncells", "T3_high_ncells"), splits):
            sub_list = sub.tolist()
            if len(sub_list) < 3:
                continue
            fct_t = evaluate_joint_residuals(
                delta_obs=delta_obs,
                pred=pred_fct,
                t_mean=t_mean,
                target_cols=target_cols,
                mode="resp_top50",
                n_boot=n_boot,
                seed=seed,
                half_a_obs=half_a_obs,
                half_b_obs=half_b_obs,
                half_valid=half_valid,
                row_subset=sub_list,
            )
            mis_hi_list = []
            for s_name in source_names:
                pred_mis = fct_from_fc_bar(fc_bar, source_ctrl_mean_libnorm[s_name])
                m_res = evaluate_joint_residuals(
                    delta_obs=delta_obs,
                    pred=pred_mis,
                    t_mean=t_mean,
                    target_cols=target_cols,
                    mode="resp_top50",
                    n_boot=n_boot,
                    seed=seed,
                    half_a_obs=half_a_obs,
                    half_b_obs=half_b_obs,
                    half_valid=half_valid,
                    row_subset=sub_list,
                )
                mis_hi_list.append(float(m_res.get("r_orth_cal_ci975_hi", float("nan"))))
            perm_t = evaluate_fct_perm_null(
                delta_obs=delta_obs,
                t_mean=t_mean,
                fc_bar=fc_bar,
                tgt_ctrl_mean_libnorm=tgt_ctrl_mean_libnorm,
                target_cols=target_cols,
                ceil_raw=float(fct_t.get("noise_ceiling_r", 0.0)),
                ceil_orth=float(fct_t.get("noise_ceiling_r_orth", 0.0)),
                mode="resp_top50",
                n_perm=n_perm,
                row_subset=sub_list,
            )
            max_mis = float(np.nanmax(mis_hi_list)) if mis_hi_list else float("nan")
            perm975 = float(perm_t["r_orth_cal"]["p975"])
            fct_lo = float(fct_t.get("r_orth_cal_ci975_lo", float("nan")))
            tertiles_out[label] = {
                "n_perts": len(sub_list),
                "n_cells_min": float(np.min(n_cells[sub_list])),
                "n_cells_max": float(np.max(n_cells[sub_list])),
                "n_cells_median": float(np.median(n_cells[sub_list])),
                "fct_r_orth_cal": fct_t.get("r_orth_cal_est"),
                "fct_r_orth_cal_ci975_lo": fct_lo,
                "fct_r_orth_cal_ci975_hi": fct_t.get("r_orth_cal_ci975_hi"),
                "max_mismatch_ci975_hi": max_mis,
                "perm100_p975": perm975,
                "perm100_median": float(perm_t["r_orth_cal"]["median"]),
                "passes_hc1": bool(np.isfinite(fct_lo) and fct_lo > 0.0),
                "passes_hc2": bool(
                    np.isfinite(fct_lo) and np.isfinite(max_mis) and np.isfinite(perm975)
                    and fct_lo > max(max_mis, perm975)
                ),
            }

    return {
        "n_eval": n_t,
        "r_orth_p_summary": {
            "n_finite": int(np.sum(np.isfinite(r_orth_p))),
            "mean": float(np.nanmean(r_orth_p)),
            "median": float(np.nanmedian(r_orth_p)),
        },
        "spearman_ln_ncells": sp_nc.to_dict("rho"),
        "spearman_abs_kp": sp_kp.to_dict("rho"),
        "ols_centered": ols_out,
        "tertiles": tertiles_out,
    }


def adjudicate_round10_group_c(
    per_line_summary: Dict[str, Dict[str, Any]],
    confirmatory_lines: Sequence[str],
) -> Dict[str, Any]:
    """
    PREREG_C §4.2 & P1 Confirmatory Decision Rule:
      - If K_conf < 2: NO confirmatory hypothesis testing (descriptive only).
      - If K_conf >= 2: 100% unanimous pass across all K_conf confirmatory lines for (H-C1 AND H-C2),
        and independently 100% unanimous pass across all K_conf lines & 6 metrics for H-C3.
    """
    k_conf = len(confirmatory_lines)
    if k_conf < 2:
        return {
            "k_conf": k_conf,
            "confirmatory_lines": list(confirmatory_lines),
            "status": "NO_CONFIRMATORY_ADJUDICATION",
            "primary_verdict": "确认性线不足，未做确认性检验",
            "hc3_verdict": "确认性线不足，未做确认性检验",
            "abstract_notice": "Insufficient confirmatory cell lines (K_conf < 2); confirmatory hypothesis testing was not performed (确认性线不足，未做确认性检验).",
        }

    hc1_pass, hc2_pass, hc3_pass = 0, 0, 0
    line_verdicts: Dict[str, Dict[str, bool]] = {}
    for line in confirmatory_lines:
        info = per_line_summary[line]
        c1 = bool(info["hc1_pass"])
        c2 = bool(info["hc2_pass"])
        c3 = bool(info["hc3_pass"])
        hc1_pass += int(c1)
        hc2_pass += int(c2)
        hc3_pass += int(c3)
        line_verdicts[line] = {"hc1": c1, "hc2": c2, "hc3": c3}

    confirmed = (hc1_pass == k_conf) and (hc2_pass == k_conf)
    hc3_confirmed = hc3_pass == k_conf
    return {
        "k_conf": k_conf,
        "confirmatory_lines": list(confirmatory_lines),
        "status": "ADJUDICATED",
        "hc1_pass_count": f"{hc1_pass}/{k_conf}",
        "hc2_pass_count": f"{hc2_pass}/{k_conf}",
        "hc3_pass_count": f"{hc3_pass}/{k_conf}",
        "primary_verdict": "已确认 (CONFIRMED)" if confirmed else "未确认 (NOT CONFIRMED)",
        "hc3_verdict": "非劣效通过 (H-C3 PASS)" if hc3_confirmed else "未通过非劣效 (H-C3 FAIL)",
        "line_verdicts": line_verdicts,
    }
