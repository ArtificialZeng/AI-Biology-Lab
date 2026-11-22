#!/usr/bin/env python3
"""
Universal Target-Free Null-Cone Evaluation Metrics for Zero-Shot Cross-Line Perturbation Prediction.

Implements:
  - N_1: Legacy 1D residual projection `span{1, Delta_TMean}` (exact parity with `round10/src/metrics_r10.py`).
  - N_4: Universal 4D nonlinear baseline-free cone
         `span{1, Delta_TMean, g_{0.05}(FC_bar), g_{1.0}(FC_bar), g_{20.0}(FC_bar)}`.
  - N_5: Universal 5D nonlinear baseline-free + perturbation-invariant shared-template cone
         `span{1, Delta_TMean, b_S, g_{0.05}(FC_bar), g_{1.0}(FC_bar), g_{20.0}(FC_bar)}`.

Every vector in N_1, N_4, N_5 is computed STRICTLY from source cell lines (zero target-line information).
Uses BLAS-batched QR projection + per-perturbation sufficient statistics for exact, millisecond-speed
joint perturbation bootstrap CIs (95% & Bonferroni 97.5%) and 100-permutation null evaluation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

EPS = 1e-12
UNIVERSAL_LOG1P_SCALES: Tuple[float, float, float] = (0.05, 1.0, 20.0)


def log1p_cone_basis(fc_bar: np.ndarray, scales: Sequence[float] = UNIVERSAL_LOG1P_SCALES) -> List[np.ndarray]:
    """Compute target-free universal `log1p(c * FC_bar) - log1p(c)` basis arrays for fixed constants `c`."""
    fc_safe = np.maximum(np.asarray(fc_bar, dtype=np.float64), 1e-12)
    return [np.log1p(float(c) * fc_safe) - np.log1p(float(c)) for c in scales]


def spearman_brown_ceiling(rho_ab: np.ndarray | float) -> np.ndarray | float:
    arr = np.asarray(rho_ab, dtype=np.float64)
    out = np.full_like(arr, np.nan, dtype=np.float64)
    ok = np.isfinite(arr) & (arr > 0.0)
    out[ok] = np.sqrt(np.maximum((2.0 * arr[ok]) / (1.0 + arr[ok]), 0.0))
    return float(out) if np.ndim(rho_ab) == 0 else out


def single_half_ceiling(rho_ab: np.ndarray | float) -> np.ndarray | float:
    arr = np.asarray(rho_ab, dtype=np.float64)
    out = np.full_like(arr, np.nan, dtype=np.float64)
    ok = np.isfinite(arr) & (arr > 0.0)
    out[ok] = np.sqrt(arr[ok])
    return float(out) if np.ndim(rho_ab) == 0 else out


def orthonormal_centered_basis(cols: Sequence[np.ndarray]) -> np.ndarray:
    """
    Build orthonormal basis Q (m, r) for centered columns (orthogonal to 1_m).
    Returns shape (m, r) with r <= len(cols).
    """
    if len(cols) == 0:
        return np.zeros((0, 0), dtype=np.float64)
    M = np.column_stack([np.asarray(c, dtype=np.float64) for c in cols])
    M = M - M.mean(axis=0, keepdims=True)
    Q, R = np.linalg.qr(M, mode="reduced")
    diag = np.abs(np.diag(R))
    if diag.size == 0:
        return np.zeros((M.shape[0], 0), dtype=np.float64)
    scale = max(1.0, float(np.max(diag)))
    keep = diag > 1e-10 * scale
    return Q[:, keep]


def project_out_centered(x: np.ndarray, Y: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Center x (m,) and Y (K, m), then project out orthonormal centered basis Q (m, r) via BLAS.
    Returns: (x_orth (m,), Y_orth (K, m), Y_centered_norm2 (K,)).
    """
    xc = x - float(np.mean(x))
    Yc = Y - np.mean(Y, axis=1, keepdims=True)
    y_tot2 = np.einsum("ij,ij->i", Yc, Yc)
    if Q.shape[1] > 0:
        xo = xc - Q @ (Q.T @ xc)
        Yo = Yc - (Yc @ Q) @ Q.T
    else:
        xo = xc
        Yo = Yc
    return xo, Yo, y_tot2


@dataclass
class PerturbationMoments:
    """Per-perturbation sufficient statistics enabling exact O(1)-per-boot joint bootstrap."""

    predictor_names: List[str]
    sxx: np.ndarray          # (N,)   ||x_p^perp||^2
    syy_orth: np.ndarray     # (N, K) ||y_{p,k}^perp||^2
    syy_tot: np.ndarray      # (N, K) ||y_{p,k,centered}||^2
    sxy: np.ndarray          # (N, K) <x_p^perp, y_{p,k}^perp>
    hsxx: np.ndarray         # (N,)   ||ha_p^perp||^2  (0 if half invalid)
    hsyy: np.ndarray         # (N,)   ||hb_p^perp||^2
    hsxy: np.ndarray         # (N,)   <ha_p^perp, hb_p^perp>
    is_split_half: bool


def compute_perturbation_moments(
    r_obs: np.ndarray,               # (N, G) observed delta - TMean
    Y_preds: np.ndarray,             # (K, N, G) predictor residual matrices (delta_pred - TMean)
    predictor_names: List[str],
    t_mean: np.ndarray,              # (N, G) source Transfer_Mean
    fc_bar: np.ndarray,              # (N, G) source geometric-mean fold-change
    source_grand_mean: Optional[np.ndarray],  # (G,) source grand mean b_S (required for N_5)
    ha_res: Optional[np.ndarray],
    hb_res: Optional[np.ndarray],
    half_valid: Optional[np.ndarray],
    idx_list: List[np.ndarray],
    null_mode: str = "N5",           # "N1", "N4", "N5"
    is_split_half: bool = False,
    split_idx_ab: Optional[Tuple[List[np.ndarray], List[np.ndarray]]] = None,
) -> PerturbationMoments:
    """
    Compute exact per-perturbation orthogonalized sufficient statistics for K predictors simultaneously.
    Because all vectors are centered per perturbation before orthogonalization, sum_p sxx[p], sum_p sxy[p, k]
    yield exact pooled per-perturbation-centered orthogonal correlations matching `metrics_r10`.
    """
    if null_mode not in ("N1", "N4", "N5"):
        raise ValueError(f"Unsupported null_mode: {null_mode}")
    if null_mode == "N5" and source_grand_mean is None:
        raise ValueError("source_grand_mean (b_S) is required for null_mode='N5'.")

    K, N, _ = Y_preds.shape
    g_basis = log1p_cone_basis(fc_bar) if null_mode in ("N4", "N5") else []

    sxx = np.zeros(N, dtype=np.float64)
    syy_orth = np.zeros((N, K), dtype=np.float64)
    syy_tot = np.zeros((N, K), dtype=np.float64)
    sxy = np.zeros((N, K), dtype=np.float64)
    hsxx = np.zeros(N, dtype=np.float64)
    hsyy = np.zeros(N, dtype=np.float64)
    hsxy = np.zeros(N, dtype=np.float64)

    min_genes = {"N1": 4, "N4": 7, "N5": 8}[null_mode]

    def build_Q(p: int, idx: np.ndarray) -> np.ndarray:
        cols = [t_mean[p, idx]]
        if null_mode == "N5":
            assert source_grand_mean is not None
            cols.append(source_grand_mean[idx])
        if null_mode in ("N4", "N5"):
            for gb in g_basis:
                cols.append(gb[p, idx])
        return orthonormal_centered_basis(cols)

    if not is_split_half:
        for p in range(N):
            idx = idx_list[p]
            if idx.size < min_genes:
                continue
            Q = build_Q(p, idx)
            xo, Yo, y_tot2 = project_out_centered(r_obs[p, idx], Y_preds[:, p, idx], Q)
            sxx[p] = float(np.dot(xo, xo))
            syy_orth[p] = np.einsum("ij,ij->i", Yo, Yo)
            syy_tot[p] = y_tot2
            sxy[p] = Yo @ xo

            if ha_res is not None and hb_res is not None and (half_valid is None or bool(half_valid[p])):
                av = ha_res[p, idx]
                bv = hb_res[p, idx]
                if np.all(np.isfinite(av)) and np.all(np.isfinite(bv)):
                    ao, bo_m, _ = project_out_centered(av, bv[None, :], Q)
                    bo = bo_m[0]
                    hsxx[p] = float(np.dot(ao, ao))
                    hsyy[p] = float(np.dot(bo, bo))
                    hsxy[p] = float(np.dot(ao, bo))
    else:
        assert split_idx_ab is not None and ha_res is not None and hb_res is not None
        idx_a, idx_b = split_idx_ab
        for p in range(N):
            if half_valid is not None and not bool(half_valid[p]):
                continue
            for obs_m, idx in ((hb_res, idx_a[p]), (ha_res, idx_b[p])):
                if idx.size < min_genes:
                    continue
                ov = obs_m[p, idx]
                av = ha_res[p, idx]
                bv = hb_res[p, idx]
                if not (np.all(np.isfinite(ov)) and np.all(np.isfinite(av)) and np.all(np.isfinite(bv))):
                    continue
                Q = build_Q(p, idx)
                xo, Yo, y_tot2 = project_out_centered(ov, Y_preds[:, p, idx], Q)
                sxx[p] += float(np.dot(xo, xo))
                syy_orth[p] += np.einsum("ij,ij->i", Yo, Yo)
                syy_tot[p] += y_tot2
                sxy[p] += Yo @ xo

                ao, bo_m, _ = project_out_centered(av, bv[None, :], Q)
                bo = bo_m[0]
                hsxx[p] += float(np.dot(ao, ao))
                hsyy[p] += float(np.dot(bo, bo))
                hsxy[p] += float(np.dot(ao, bo))

    return PerturbationMoments(
        predictor_names=list(predictor_names),
        sxx=sxx,
        syy_orth=syy_orth,
        syy_tot=syy_tot,
        sxy=sxy,
        hsxx=hsxx,
        hsyy=hsyy,
        hsxy=hsxy,
        is_split_half=is_split_half,
    )


def _moments_to_scores(
    sxx: np.ndarray,       # (B,) or scalar
    syy_orth: np.ndarray,  # (B, K) or (K,)
    syy_tot: np.ndarray,   # (B, K) or (K,)
    sxy: np.ndarray,       # (B, K) or (K,)
    hsxx: np.ndarray,      # (B,) or scalar
    hsyy: np.ndarray,      # (B,) or scalar
    hsxy: np.ndarray,      # (B,) or scalar
    is_split_half: bool,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (r_orth_cal, r_semi_cal, eta, rho_ab, ceil) for point or batched bootstrap arrays."""
    h_den = np.sqrt(np.maximum(hsxx * hsyy, 0.0))
    rho = np.where(h_den > EPS, np.clip(hsxy / np.maximum(h_den, EPS), -1.0, 1.0), np.nan)
    ceil = single_half_ceiling(rho) if is_split_half else spearman_brown_ceiling(rho)

    sxx_exp = np.expand_dims(sxx, axis=-1) if np.ndim(sxx) > 0 else sxx
    ceil_exp = np.expand_dims(ceil, axis=-1) if np.ndim(ceil) > 0 else ceil

    den_o = np.sqrt(np.maximum(sxx_exp * syy_orth, 0.0))
    den_t = np.sqrt(np.maximum(sxx_exp * syy_tot, 0.0))
    raw_o = np.where(den_o > EPS, np.clip(sxy / np.maximum(den_o, EPS), -1.0, 1.0), np.nan)
    raw_s = np.where(den_t > EPS, np.clip(sxy / np.maximum(den_t, EPS), -1.0, 1.0), np.nan)
    eta = np.where(syy_tot > EPS, np.sqrt(np.maximum(syy_orth / np.maximum(syy_tot, EPS), 0.0)), np.nan)

    valid_c = np.isfinite(ceil_exp) & (ceil_exp > 1e-8)
    cal_o = np.where(valid_c & np.isfinite(raw_o), np.clip(raw_o / np.maximum(ceil_exp, 1e-8), -1.0, 1.0), np.nan)
    cal_s = np.where(valid_c & np.isfinite(raw_s), np.clip(raw_s / np.maximum(ceil_exp, 1e-8), -1.0, 1.0), np.nan)
    return cal_o, cal_s, eta, rho, ceil


def joint_bootstrap_summary(
    mom: PerturbationMoments,
    n_boot: int = 4000,
    seed: int = 20261124,
    paired_specs: Sequence[Tuple[str, str, str]] = (),
) -> Dict[str, object]:
    """
    Exact joint perturbation bootstrap computing 95% and Bonferroni 97.5% CIs for every predictor
    (both `orth` and `semi` calibrated scores, plus `eta`) and optional paired differences
    `(label, pred_a, pred_b)` using identical bootstrap draws!
    """
    N = mom.sxx.shape[0]
    # Point estimates
    pt_o, pt_s, pt_eta, pt_rho, pt_ceil = _moments_to_scores(
        float(np.sum(mom.sxx)),
        np.sum(mom.syy_orth, axis=0),
        np.sum(mom.syy_tot, axis=0),
        np.sum(mom.sxy, axis=0),
        float(np.sum(mom.hsxx)),
        float(np.sum(mom.hsyy)),
        float(np.sum(mom.hsxy)),
        mom.is_split_half,
    )

    # Multinomial bootstrap weights W: (n_boot, N) float64 -> exact O(n_boot * N * K) BLAS matmul!
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, N, size=(n_boot, N))
    W = np.zeros((n_boot, N), dtype=np.float64)
    rows = np.repeat(np.arange(n_boot, dtype=np.int64), N)
    np.add.at(W, (rows, idx.ravel()), 1.0)

    b_sxx = W @ mom.sxx
    b_syy_o = W @ mom.syy_orth
    b_syy_t = W @ mom.syy_tot
    b_sxy = W @ mom.sxy
    b_hsxx = W @ mom.hsxx
    b_hsyy = W @ mom.hsyy
    b_hsxy = W @ mom.hsxy

    bo, bs, _, _, _ = _moments_to_scores(b_sxx, b_syy_o, b_syy_t, b_sxy, b_hsxx, b_hsyy, b_hsxy, mom.is_split_half)

    def summarize_vec(pt: float, vec: np.ndarray) -> Dict[str, float]:
        fin = vec[np.isfinite(vec)]
        vf = float(fin.size / max(1, vec.size))
        if fin.size == 0:
            return {"est": float(pt), "ci95_lo": float("nan"), "ci95_hi": float("nan"),
                    "ci975_lo": float("nan"), "ci975_hi": float("nan"), "valid_frac": vf}
        q0125, q025, q975, q9875 = np.percentile(fin, [1.25, 2.5, 97.5, 98.75])
        return {
            "est": float(pt),
            "ci95_lo": float(q025),
            "ci95_hi": float(q975),
            "ci975_lo": float(q0125),
            "ci975_hi": float(q9875),
            "valid_frac": vf,
        }

    name_to_k = {nm: k for k, nm in enumerate(mom.predictor_names)}
    per_pred: Dict[str, Dict[str, object]] = {}
    for k, nm in enumerate(mom.predictor_names):
        per_pred[nm] = {
            "orth": summarize_vec(float(pt_o[k]), bo[:, k]),
            "semi": summarize_vec(float(pt_s[k]), bs[:, k]),
            "eta": float(pt_eta[k]),
        }

    paired_out: Dict[str, Dict[str, object]] = {}
    for label, pa, pb in paired_specs:
        if pa in name_to_k and pb in name_to_k:
            ka, kb = name_to_k[pa], name_to_k[pb]
            paired_out[label] = {
                "orth": summarize_vec(float(pt_o[ka] - pt_o[kb]), bo[:, ka] - bo[:, kb]),
                "semi": summarize_vec(float(pt_s[ka] - pt_s[kb]), bs[:, ka] - bs[:, kb]),
            }

    return {
        "n_perturbations": int(N),
        "n_boot": int(n_boot),
        "seed": int(seed),
        "rho_ab": float(pt_rho),
        "ceiling": float(pt_ceil),
        "predictors": per_pred,
        "paired": paired_out,
        "_boot_semi": bs,
        "_boot_orth": bo,
        "_pt_semi": pt_s,
        "_pt_orth": pt_o,
        "_name_to_k": name_to_k,
    }
