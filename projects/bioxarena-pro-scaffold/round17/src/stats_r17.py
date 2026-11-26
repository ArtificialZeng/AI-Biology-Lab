#!/usr/bin/env python3
"""
Round 17 Statistical, Cross-Fitting, Variance-Component, Label-Free Selection,
Mechanism (`J5`), & Holm-Bonferroni Adjudication Engine (`PREREG_J` v2).

Strictly imports frozen mathematical primitives from `round16/src/` and `round10/src/`
without modifying any file under `round9/`–`round16/` or any locked `round17/` file:
  - Primary endpoint (`r_orth_cal`) & secondary (`r_dorth_cal`) bit-compatible with
    frozen `round16/src/metrics_r16.py`.
  - Selection-free 56-index sweep (`b in 0..55`, never dropped post-inference per `J4`),
    shared `n_boot x 56` joint perturbation bootstrap (`seed = 20261701`, `n_boot = 4000`).
  - `H-J1`: Cross-fitted perturbation half-split spread (`s = 0` locked primary, `s = 0..19`
    20-split stability diagnostic), exact within-half bootstrap (`20261701` on P2,
    `20261701 + 10000` on P1), noise-debiased between-index SD (`tau_batch`), centered-bootstrap
    homogeneity test (`p_var`), `J2` floor-effect detectable-signal guard, and 6-branch rule.
  - `H-J2`: Per-index `PREREG_I §3.2` `97.5%` CI `STATE - FCT` calls across all 56 indices,
    flip fraction `f_flip`, sign-reversal detection, `J2` qualifier on `CALL_STABLE`, AND
    companion Bonferroni-56 hardened check under both normal-calibrated `z = norm.ppf(1 - 0.025/112)`
    and `n_boot = 20,000` empirical quantiles side-by-side.
  - `H-J3`: Label-free `L1_ens56` (`H-J3a`) and `L2_ntc_recon_match` (`H-J3b`, Ruling 3 Pearson
    primary with mandatory `top1 - top2` margin, negative-L2 distance secondary, per-fold
    selection secondary), paired vs `b = 0`, out-of-sample vs cross-fitted oracle on `s = 0`,
    vs `FCT`, plus secondary neutral-batch comparators (`bypass_zero`, `mean_batch_emb`,
    `multi_batch_ens8`, `medoid_idx`).
  - `J5a`: Cross-fitted `G4` 2x2 (`X_ctrl x batch`) decomposition at `(b*_1, b_dag_1)` (and
    reversed `(b_dag_1, b*_1)`) evaluated strictly on disjoint `P2`, auxiliary `(b*_2, b_dag_2)`
    on `P1`, and pre-fixed continuity pairs `(0, 17)` / `(17, 0)` & `(0, 35)` / `(35, 0)` on all `N_t`.
  - `J5b`: Cross-line & cross-model 56-index Spearman rank correlation with exact `10,000`-perm
    index-label null (`seed = 20261701`), BH FDR, and top-5 / bottom-5 shared index overlaps.
  - Multiplicity: Holm-Bonferroni FWER <= 0.05 across `M = 3 * K_conf` per model (`STATE_SE`,
    `STATE_HVG` separately) plus pooled-model `6 * K_conf` sensitivity, and locked `k / K_conf`
    cohort claim rule (`general claim iff K_conf >= 2 and k >= 2`, zero pooling across lines).
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
from scipy.stats import chi2, norm, rankdata, spearmanr

ROUND17_SRC_DIR = Path(__file__).resolve().parent
REPO_ROOT = ROUND17_SRC_DIR.parents[1]
R16_SRC_DIR = REPO_ROOT / "round16" / "src"
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
for _p in (str(R16_SRC_DIR), str(R10_SRC_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from adjudicate_r16 import MIN_BOOT_VALID_FRAC, adjudicate_claim_i2_line  # noqa: E402
from metrics_r10 import (  # noqa: E402
    _extract_masked_idx_per_pert,
    _moments_to_corr_and_slope,
    disattenuated_r,
    spearman_brown_ceiling,
)
from metrics_r16 import (  # noqa: E402
    _collect_mode_moments,
    _moments_to_safe_corr,
    compute_paired_contrast_r16,
)
from run_eval_r16 import (  # noqa: E402
    compute_g4_factorial_decomposition_r16,
    compute_practical_accuracy_metrics_r16,
)

DEFAULT_R17_SEED = 20261701
DEFAULT_R17_RESAMPLE_SEED = 20261702
DEFAULT_R17_N_BOOT = 4000
DEFAULT_R17_N_BOOT_BONF56 = 20000
DEFAULT_R17_N_PERM = 10000
DEFAULT_R17_N_SPLITS_DIAG = 20
N_BATCH_INDICES_STATE = 56
DEFAULT_PRACTICAL_DELTA = 0.04
DEFAULT_TAU_THRESHOLD = 0.02  # 2 * tau_batch > 0.04 <=> tau_batch > 0.02
ALPHA_FAMILY = 0.05
ALPHA_LINE_BASE = ALPHA_FAMILY / 3.0  # 0.016666666666666666 (98.333% two-sided CI)
ALPHA_HI2 = 0.025  # 97.5% two-sided CI (PREREG_I §3.2)
ALPHA_BONF56 = ALPHA_HI2 / float(N_BATCH_INDICES_STATE)  # 0.025 / 56
Z_BONF56 = float(norm.ppf(1.0 - ALPHA_BONF56 / 2.0))  # 3.511653...
HALF_P1_BOOT_SEED_OFFSET = 10000

# R16 §4.3 parametric precision constants
R16_SE_SCALE_HI = 0.1502  # a549_ko-like high rho_ab_orth anchor
R16_SE_SCALE_LO = 0.1750  # cd4t_rest-like conservative low rho_ab_orth anchor


# ==============================================================================
# Core Bootstrap & CI Utilities
# ==============================================================================


def _two_sided_ci_from_vals(vals: np.ndarray, alpha: float) -> Tuple[float, float]:
    """Two-sided `(1 - alpha)` percentile interval from finite bootstrap draws."""
    if vals.size == 0 or not np.isfinite(alpha) or alpha <= 0.0 or alpha >= 1.0:
        return float("nan"), float("nan")
    lo_q = 100.0 * (alpha / 2.0)
    hi_q = 100.0 * (1.0 - alpha / 2.0)
    lo, hi = np.percentile(vals, [lo_q, hi_q])
    return float(lo), float(hi)


def _empirical_two_sided_p(vals: np.ndarray) -> float:
    """Two-sided percentile-bootstrap p-value vs 0.0 with (+1)/(B+1) finite-sample floor."""
    n_v = int(vals.size)
    if n_v == 0:
        return float("nan")
    n_le0 = int(np.sum(vals <= 0.0))
    n_ge0 = int(np.sum(vals >= 0.0))
    return float(min(1.0, 2.0 * (1 + min(n_le0, n_ge0)) / float(n_v + 1)))


def _empirical_positive_two_sided_p(vals: np.ndarray) -> float:
    """
    Two-sided-equivalent p-value for a positive-direction contrast (`ci_{1-alpha}_lo > 0`):
    `min(1.0, 2 * (1 + #{vals <= 0}) / (B_valid + 1))`. Equals `_empirical_two_sided_p` when
    median > 0, and maps non-positive distributions conservatively to 1.0 so negative spreads
    cannot consume Holm-Bonferroni ranks for `H-J1`.
    """
    n_v = int(vals.size)
    if n_v == 0:
        return float("nan")
    n_le0 = int(np.sum(vals <= 0.0))
    return float(min(1.0, 2.0 * (1 + n_le0) / float(n_v + 1)))


def summarize_boot_r17(
    point: float,
    boots: np.ndarray,
    n_items: int,
    n_boot: int,
    extra_alphas: Optional[Mapping[str, float]] = None,
) -> Dict[str, Any]:
    """
    Summarize bootstrap draws with standard Round 17 intervals:
      - `ci95`: `[2.5%, 97.5%]` (`alpha = 0.05`)
      - `ci975`: `[1.25%, 98.75%]` (`alpha = 0.025`, `PREREG_I` `H-I2`)
      - `ci9833`: `[0.8333%, 99.1667%]` (`alpha = 0.05/3`, base 3-test line Bonferroni & `J2`)
      - `p_two_sided`: empirical two-sided bootstrap p-value vs 0.0
      - Optional `extra_alphas` (e.g. Holm step-down adjusted CI).
    """
    fin = np.asarray(boots, dtype=np.float64)
    valid_mask = np.isfinite(fin)
    valid_frac = float(np.mean(valid_mask)) if fin.size > 0 else 0.0
    vals = fin[valid_mask]

    if vals.size == 0:
        out: Dict[str, Any] = {
            "estimate": float(point),
            "ci95_lo": float("nan"),
            "ci95_hi": float("nan"),
            "ci975_lo": float("nan"),
            "ci975_hi": float("nan"),
            "ci9833_lo": float("nan"),
            "ci9833_hi": float("nan"),
            "boot_se": float("nan"),
            "p_two_sided": float("nan"),
            "valid_frac": valid_frac,
            "n_valid_boot": 0,
            "n_items": int(n_items),
            "n_boot": int(n_boot),
        }
        if extra_alphas:
            for k in extra_alphas:
                out[f"{k}_lo"] = float("nan")
                out[f"{k}_hi"] = float("nan")
        return out

    lo95, hi95 = _two_sided_ci_from_vals(vals, 0.05)
    lo975, hi975 = _two_sided_ci_from_vals(vals, ALPHA_HI2)
    lo9833, hi9833 = _two_sided_ci_from_vals(vals, ALPHA_LINE_BASE)
    out = {
        "estimate": float(point),
        "ci95_lo": lo95,
        "ci95_hi": hi95,
        "ci975_lo": lo975,
        "ci975_hi": hi975,
        "ci9833_lo": lo9833,
        "ci9833_hi": hi9833,
        "boot_se": float(np.std(vals, ddof=1)) if vals.size > 1 else 0.0,
        "p_two_sided": _empirical_two_sided_p(vals),
        "valid_frac": valid_frac,
        "n_valid_boot": int(vals.size),
        "n_items": int(n_items),
        "n_boot": int(n_boot),
    }
    if extra_alphas:
        for k, a_val in extra_alphas.items():
            lo_e, hi_e = _two_sided_ci_from_vals(vals, float(a_val))
            out[f"{k}_lo"] = lo_e
            out[f"{k}_hi"] = hi_e
    return out


def _compute_half_ceiling_draws(
    half_mom: Optional[Tuple[np.ndarray, ...]],
    draws: np.ndarray,
) -> Tuple[float, float, np.ndarray]:
    """Evaluate observed `(rho_ab_orth, ceiling)` and bootstrap `b_ceil` ONCE for shared draws."""
    n_boot = draws.shape[0]
    if half_mom is None:
        return float("nan"), float("nan"), np.full(n_boot, float("nan"), dtype=np.float64)
    hnv, hsx, hsy, hsxx, hsyy, hsxy = half_mom
    _, pt_rho, _, _, _ = _moments_to_corr_and_slope(
        float(hnv.sum()), float(hsx.sum()), float(hsy.sum()),
        float(hsxx.sum()), float(hsyy.sum()), float(hsxy.sum()),
    )
    pt_rho_f = float(pt_rho)
    pt_ceil = float(spearman_brown_ceiling(pt_rho_f))

    _, b_rho, _, _, _ = _moments_to_corr_and_slope(
        hnv[draws].sum(axis=1), hsx[draws].sum(axis=1), hsy[draws].sum(axis=1),
        hsxx[draws].sum(axis=1), hsyy[draws].sum(axis=1), hsxy[draws].sum(axis=1),
    )
    b_rho_arr = np.asarray(b_rho, dtype=np.float64)
    b_ceil = np.full(n_boot, float("nan"), dtype=np.float64)
    pos = np.isfinite(b_rho_arr) & (b_rho_arr > 0.0)
    if np.any(pos):
        sb = (2.0 * b_rho_arr[pos]) / (1.0 + b_rho_arr[pos])
        b_ceil[pos] = np.sqrt(np.maximum(sb, 0.0))
    return pt_rho_f, pt_ceil, b_ceil


def _eval_pred_mom_with_ceiling(
    pred_mom: Tuple[np.ndarray, ...],
    pt_ceil: float,
    b_ceil: np.ndarray,
    draws: np.ndarray,
    zero_null_to_zero: bool = False,
) -> Tuple[float, np.ndarray, float]:
    """Bit-compatible fast evaluation of `(pt_cal, b_cal, pt_raw)` reusing shared `b_ceil`."""
    nv, sx, sy, sxx, syy, sxy = pred_mom
    pt_r = float(
        _moments_to_safe_corr(
            float(nv.sum()), float(sx.sum()), float(sy.sum()),
            float(sxx.sum()), float(syy.sum()), float(sxy.sum()),
            zero_null_to_zero=zero_null_to_zero,
        )
    )
    b_r = np.asarray(
        _moments_to_safe_corr(
            nv[draws].sum(axis=1), sx[draws].sum(axis=1), sy[draws].sum(axis=1),
            sxx[draws].sum(axis=1), syy[draws].sum(axis=1), sxy[draws].sum(axis=1),
            zero_null_to_zero=zero_null_to_zero,
        ),
        dtype=np.float64,
    )

    if zero_null_to_zero and pt_r == 0.0 and np.isfinite(pt_ceil) and pt_ceil > 0.0:
        pt_cal = 0.0
    else:
        pt_cal = float(disattenuated_r(pt_r, pt_ceil))

    b_cal = np.full(b_r.shape[0], float("nan"), dtype=np.float64)
    ok_c = np.isfinite(b_r) & np.isfinite(b_ceil) & (b_ceil > 1e-8)
    if np.any(ok_c):
        b_cal[ok_c] = np.clip(b_r[ok_c] / b_ceil[ok_c], -1.0, 1.0)
    if zero_null_to_zero:
        z_ok = (b_r == 0.0) & np.isfinite(b_ceil) & (b_ceil > 0.0)
        b_cal[z_ok] = 0.0
    return pt_cal, b_cal, pt_r


def _eval_point_from_moments_subset(
    pred_mom: Tuple[np.ndarray, ...],
    half_mom: Optional[Tuple[np.ndarray, ...]],
    subset_idx: np.ndarray,
    zero_null_to_zero: bool = False,
) -> Tuple[float, float, float, float]:
    """Evaluate point `(r_cal, r_raw, rho_ab_orth, ceiling)` on a subset of row indices."""
    if subset_idx.size == 0:
        return float("nan"), float("nan"), float("nan"), float("nan")
    nv, sx, sy, sxx, syy, sxy = pred_mom
    pt_r = float(
        _moments_to_safe_corr(
            float(nv[subset_idx].sum()), float(sx[subset_idx].sum()), float(sy[subset_idx].sum()),
            float(sxx[subset_idx].sum()), float(syy[subset_idx].sum()), float(sxy[subset_idx].sum()),
            zero_null_to_zero=zero_null_to_zero,
        )
    )
    if half_mom is None:
        return float("nan"), pt_r, float("nan"), float("nan")
    hnv, hsx, hsy, hsxx, hsyy, hsxy = half_mom
    _, pt_rho, _, _, _ = _moments_to_corr_and_slope(
        float(hnv[subset_idx].sum()), float(hsx[subset_idx].sum()), float(hsy[subset_idx].sum()),
        float(hsxx[subset_idx].sum()), float(hsyy[subset_idx].sum()), float(hsxy[subset_idx].sum()),
    )
    pt_rho_f = float(pt_rho)
    pt_ceil = float(spearman_brown_ceiling(pt_rho_f))
    if zero_null_to_zero and pt_r == 0.0 and np.isfinite(pt_ceil) and pt_ceil > 0.0:
        pt_cal = 0.0
    else:
        pt_cal = float(disattenuated_r(pt_r, pt_ceil))
    return pt_cal, pt_r, pt_rho_f, pt_ceil


def _subset_moments(
    mom: Tuple[np.ndarray, ...],
    subset_idx: np.ndarray,
) -> Tuple[np.ndarray, ...]:
    return tuple(arr[subset_idx] for arr in mom)


def _argmax_lowest_tie(vals: np.ndarray) -> int:
    """Argmax with explicit lowest-integer-index tie-break; non-finite mapped to -inf."""
    arr = np.asarray(vals, dtype=np.float64)
    safe = np.where(np.isfinite(arr), arr, -np.inf)
    return int(np.argmax(safe))


def _argmin_lowest_tie(vals: np.ndarray) -> int:
    """Argmin with explicit lowest-integer-index tie-break; non-finite mapped to +inf."""
    arr = np.asarray(vals, dtype=np.float64)
    safe = np.where(np.isfinite(arr), arr, np.inf)
    return int(np.argmin(safe))


def _rank_1based_desc_lowest_tie(vals: np.ndarray, idx: int) -> int:
    """Exact 1-based descending rank (1 = highest score, lowest-index tie-break)."""
    arr = np.asarray(vals, dtype=np.float64)
    safe = np.where(np.isfinite(arr), arr, -np.inf)
    order = sorted(range(arr.size), key=lambda b: (-safe[b], b))
    return int(order.index(int(idx)) + 1)


# ==============================================================================
# 1. Perturbation Stratified Half-Split (`PREREG_J §3.2`)
# ==============================================================================


def make_stratified_half_split_r17(
    target_names: Sequence[str],
    half_valid: np.ndarray,
    seed: int = DEFAULT_R17_SEED,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    `PREREG_J §3.2` deterministic perturbation half-split (`s = 0` uses `seed = 20261701`,
    multi-split diagnostic `s = 0..19` uses `seed = 20261701 + s`):
      1. Sort evaluated perturbations canonically by target symbol.
      2. Permute via `np.random.default_rng(seed)` stratified by `half_valid` status
         whenever any `half_valid == False` perturbations exist (within-stratum permutation
         with `half_valid == True` first, then `half_valid == False`, concatenated).
      3. Exact alternating allocation producing `|P1| = floor(N/2)` and `|P2| = ceil(N/2)`:
         odd positions (`1::2`) -> `P1`, even positions (`0::2`) -> `P2`.
      4. Sort row indices canonically within `P1` and `P2` (by canonical target order).
    Returns `(p1_indices, p2_indices)` indexing into the input array order.
    """
    n = len(target_names)
    if n < 2:
        raise ValueError(f"Half-split requires at least 2 perturbations, got N={n}.")
    hv = np.asarray(half_valid, dtype=bool)
    if hv.shape != (n,):
        raise ValueError(f"half_valid shape {hv.shape} != ({n},).")

    # Canonical order by target symbol (stable tie-break by original index)
    canon_order = np.asarray(
        sorted(range(n), key=lambda i: (str(target_names[i]), i)),
        dtype=np.int64,
    )
    hv_canon = hv[canon_order]
    rng = np.random.default_rng(int(seed))

    if bool(np.any(~hv_canon)) and bool(np.any(hv_canon)):
        idx_true = canon_order[hv_canon]
        idx_false = canon_order[~hv_canon]
        perm_true = rng.permutation(idx_true)
        perm_false = rng.permutation(idx_false)
        combined = np.concatenate([perm_true, perm_false], axis=0)
    else:
        combined = rng.permutation(canon_order)

    # Exact alternating allocation: odd positions -> P1 (floor(N/2)), even -> P2 (ceil(N/2))
    p1_raw = combined[1::2]
    p2_raw = combined[0::2]
    assert len(p1_raw) == n // 2 and len(p2_raw) == (n + 1) // 2

    # Sort each half canonically by target symbol order so bootstrap indexing is order-canonical
    pos_in_canon = {int(orig_i): c_pos for c_pos, orig_i in enumerate(canon_order)}
    p1_sorted = np.asarray(sorted(p1_raw.tolist(), key=lambda i: pos_in_canon[i]), dtype=np.int64)
    p2_sorted = np.asarray(sorted(p2_raw.tolist(), key=lambda i: pos_in_canon[i]), dtype=np.int64)
    return p1_sorted, p2_sorted


# ==============================================================================
# 2. Selection-Free 56-Index Sweep + Joint Bootstrap (`§3.1`, `§3.2`, `§3.3`, `§3.4`)
# ==============================================================================


def compute_sweep56_moments_r17(
    delta_obs: np.ndarray,
    delta_pred_56: np.ndarray,
    fct_delta_pred: np.ndarray,
    delta_transfer: np.ndarray,
    target_gene_idx: Sequence[int],
    source_grand_mean: Optional[np.ndarray] = None,
    delta_half_a: Optional[np.ndarray] = None,
    delta_half_b: Optional[np.ndarray] = None,
    half_valid: Optional[np.ndarray] = None,
    extra_preds: Optional[Mapping[str, np.ndarray]] = None,
) -> Dict[str, Any]:
    """
    Collect per-perturbation residual moments ONCE across:
      - all 56 indices `b = 0..55` (never dropped per `J4`)
      - `FCT`
      - `L1_ens56` (exact uniform mean of `delta_pred_56` across axis 0)
      - any `extra_preds` (e.g. `bypass_zero`, `mean_batch_emb`, `multi_batch_ens8`, `medoid_idx`)
    using frozen `round16/src/metrics_r16.py::_collect_mode_moments`.
    """
    d_obs = np.asarray(delta_obs, dtype=np.float64)
    p56 = np.asarray(delta_pred_56, dtype=np.float64)
    fct = np.asarray(fct_delta_pred, dtype=np.float64)
    tm = np.asarray(delta_transfer, dtype=np.float64)
    b_s = np.asarray(source_grand_mean, dtype=np.float64) if source_grand_mean is not None else None

    if p56.ndim != 3 or p56.shape[0] != N_BATCH_INDICES_STATE:
        raise ValueError(f"delta_pred_56 must have shape (56, P, G), got {p56.shape}.")
    if d_obs.shape != p56.shape[1:] or fct.shape != d_obs.shape or tm.shape != d_obs.shape:
        raise ValueError(
            f"Shape mismatch: obs={d_obs.shape}, p56={p56.shape}, fct={fct.shape}, tm={tm.shape}."
        )

    rows, idx_list = _extract_masked_idx_per_pert(
        d_obs, target_gene_idx, gene_keep_mask=None, mode="all_nontarget", top_k=50
    )
    m = len(rows)
    if m == 0:
        raise ValueError("Zero valid evaluation rows after non-target masking.")

    r_obs = d_obs - tm
    ha_res = (np.asarray(delta_half_a, dtype=np.float64) - tm) if delta_half_a is not None else None
    hb_res = (np.asarray(delta_half_b, dtype=np.float64) - tm) if delta_half_b is not None else None

    l1_pred = np.mean(p56, axis=0)
    extra_keys = sorted(extra_preds.keys()) if extra_preds else []

    pred_res_list: List[np.ndarray] = [p56[b] - tm for b in range(N_BATCH_INDICES_STATE)]
    pred_res_list.append(fct - tm)
    pred_res_list.append(l1_pred - tm)
    for ek in extra_keys:
        assert extra_preds is not None
        pred_res_list.append(np.asarray(extra_preds[ek], dtype=np.float64) - tm)

    all_moms, half_mom = _collect_mode_moments(
        r_obs, pred_res_list, tm, rows, idx_list, b_s, ha_res, hb_res, half_valid
    )
    moms_56 = all_moms[:N_BATCH_INDICES_STATE]
    mom_fct = all_moms[N_BATCH_INDICES_STATE]
    mom_l1 = all_moms[N_BATCH_INDICES_STATE + 1]
    moms_extra = {ek: all_moms[N_BATCH_INDICES_STATE + 2 + j] for j, ek in enumerate(extra_keys)}

    return {
        "rows": np.asarray(rows, dtype=np.int64),
        "m": int(m),
        "double_orth": b_s is not None,
        "moms_56": moms_56,
        "mom_fct": mom_fct,
        "mom_l1": mom_l1,
        "moms_extra": moms_extra,
        "half_mom": half_mom,
    }


def compute_variance_component_r17(
    r_hat_56: np.ndarray,
    r_boot_56: np.ndarray,
) -> Dict[str, Any]:
    """
    `PREREG_J §3.2` Item 3 — Noise-Debiased Between-Index SD (`tau_batch`) and
    Centered-Bootstrap Homogeneity Test (`p_var`):
      - `S_obs_sq = (1/55) * ||C_56 r_hat||_2^2`
      - `V_noise = (1/55) * sum_{b=0}^{55} Var_boot(R_boot[u, b] - mean_b R_boot[u, :])` (ddof=1)
      - `tau_batch = sqrt(max(0, S_obs_sq - V_noise))`
      - `Z[u, :] = R_boot[u, :] - mean_u(R_boot[:, :])`
      - `p_var = (1 + sum_u I((1/55) ||C_56 Z[u]||_2^2 >= S_obs_sq)) / (n_valid_rows + 1)`
    """
    r_hat = np.asarray(r_hat_56, dtype=np.float64)
    r_boot = np.asarray(r_boot_56, dtype=np.float64)
    if r_hat.shape != (N_BATCH_INDICES_STATE,) or r_boot.ndim != 2 or r_boot.shape[1] != N_BATCH_INDICES_STATE:
        raise ValueError(f"Expected (56,) and (B, 56), got {r_hat.shape} and {r_boot.shape}.")

    df = float(N_BATCH_INDICES_STATE - 1)  # 55.0
    n_boot = int(r_boot.shape[0])
    row_valid = np.all(np.isfinite(r_boot), axis=1)
    joint_valid_frac = float(np.mean(row_valid)) if n_boot > 0 else 0.0
    per_idx_valid_frac = np.mean(np.isfinite(r_boot), axis=0) if n_boot > 0 else np.zeros(N_BATCH_INDICES_STATE)
    min_idx_valid_frac = float(np.min(per_idx_valid_frac))

    if not np.all(np.isfinite(r_hat)) or int(np.sum(row_valid)) < 2:
        return {
            "s_obs_sq": float("nan"),
            "s_obs": float("nan"),
            "v_noise": float("nan"),
            "tau_batch": float("nan"),
            "two_tau_batch": float("nan"),
            "p_var": float("nan"),
            "n_valid_joint_rows": int(np.sum(row_valid)),
            "joint_valid_frac": joint_valid_frac,
            "min_per_index_valid_frac": min_idx_valid_frac,
            "n_boot": n_boot,
        }

    r_centered = r_hat - float(np.mean(r_hat))
    s_obs_sq = float(np.dot(r_centered, r_centered) / df)
    s_obs = float(np.sqrt(max(s_obs_sq, 0.0)))

    rv = r_boot[row_valid]
    n_v = int(rv.shape[0])
    # D[u, b] = R_boot[u, b] - mean_b(R_boot[u, :])
    d_mat = rv - np.mean(rv, axis=1, keepdims=True)
    v_noise = float(np.sum(np.var(d_mat, axis=0, ddof=1)) / df)
    tau_batch = float(np.sqrt(max(0.0, s_obs_sq - v_noise)))

    # Centered bootstrap null: Z[u, :] = R_boot[u, :] - mean_u(R_boot[:, :])
    z_mat = rv - np.mean(rv, axis=0, keepdims=True)
    cz_mat = z_mat - np.mean(z_mat, axis=1, keepdims=True)
    t_boot = np.sum(cz_mat * cz_mat, axis=1) / df
    p_var = float((1 + int(np.sum(t_boot >= s_obs_sq))) / float(n_v + 1))

    return {
        "s_obs_sq": s_obs_sq,
        "s_obs": s_obs,
        "v_noise": v_noise,
        "tau_batch": tau_batch,
        "two_tau_batch": 2.0 * tau_batch,
        "p_var": p_var,
        "t_boot_mean": float(np.mean(t_boot)),
        "t_boot_q95": float(np.percentile(t_boot, 95.0)),
        "n_valid_joint_rows": n_v,
        "joint_valid_frac": joint_valid_frac,
        "min_per_index_valid_frac": min_idx_valid_frac,
        "n_boot": n_boot,
    }


def compute_cross_fitted_spread_r17(
    moments_bundle: Dict[str, Any],
    target_names: Sequence[str],
    half_valid: np.ndarray,
    l2_selected_idx: Optional[int] = None,
    n_boot: int = DEFAULT_R17_N_BOOT,
    seed: int = DEFAULT_R17_SEED,
    n_splits_diag: int = DEFAULT_R17_N_SPLITS_DIAG,
) -> Dict[str, Any]:
    """
    `PREREG_J §3.2` Item 2 & `§3.4` Out-of-Sample Oracle Contrast:
      - Locked primary split `s = 0` (`seed = 20261701`):
          Selects `(b*_1, b_dag_1)` strictly on `P1`, evaluates on `P2`.
          Selects `(b*_2, b_dag_2)` strictly on `P2`, evaluates on `P1`.
          Forms `Delta_CF,0 = 0.5 * (Delta_{1->2} + Delta_{2->1})` and
          `r_CF_Oracle,0 = 0.5 * (r_P2(b*_1) + r_P1(b*_2))`.
      - Exact within-half bootstrap (`n_boot = 4,000`):
          `P2` draws with `seed`, `P1` draws with `seed + 10000` (independent).
          Also computes paired out-of-sample contrasts `L1_ens56 - CF_Oracle` and
          `L2_ntc_recon_match - CF_Oracle` on those exact same within-half draws!
      - Multi-split stability check (`s = 0..n_splits_diag-1`, `seed + s`).
    """
    rows = moments_bundle["rows"]
    m = int(moments_bundle["m"])
    if m != len(target_names):
        raise ValueError(f"Expected all {len(target_names)} perturbations valid, got m={m}.")
    zero_null = bool(moments_bundle["double_orth"])
    moms_56 = moments_bundle["moms_56"]
    mom_l1 = moments_bundle["mom_l1"]
    half_mom = moments_bundle["half_mom"]

    # Map original perturbation row -> position in `rows` (usually identity 0..m-1)
    row_to_pos = {int(r): pos for pos, r in enumerate(rows)}

    def _eval_split_point(split_seed: int) -> Dict[str, Any]:
        p1_orig, p2_orig = make_stratified_half_split_r17(target_names, half_valid, seed=split_seed)
        p1 = np.asarray([row_to_pos[int(i)] for i in p1_orig], dtype=np.int64)
        p2 = np.asarray([row_to_pos[int(i)] for i in p2_orig], dtype=np.int64)

        r_p1 = np.empty(N_BATCH_INDICES_STATE, dtype=np.float64)
        r_p2 = np.empty(N_BATCH_INDICES_STATE, dtype=np.float64)
        rho_p1 = float("nan")
        rho_p2 = float("nan")
        for b in range(N_BATCH_INDICES_STATE):
            v1, _, r1, _ = _eval_point_from_moments_subset(moms_56[b], half_mom, p1, zero_null_to_zero=zero_null)
            v2, _, r2, _ = _eval_point_from_moments_subset(moms_56[b], half_mom, p2, zero_null_to_zero=zero_null)
            r_p1[b] = v1
            r_p2[b] = v2
            if b == 0:
                rho_p1, rho_p2 = r1, r2

        b_star_1 = _argmax_lowest_tie(r_p1)
        b_dag_1 = _argmin_lowest_tie(r_p1)
        b_star_2 = _argmax_lowest_tie(r_p2)
        b_dag_2 = _argmin_lowest_tie(r_p2)

        d_1to2 = float(r_p2[b_star_1] - r_p2[b_dag_1])
        d_2to1 = float(r_p1[b_star_2] - r_p1[b_dag_2])
        delta_cf = 0.5 * (d_1to2 + d_2to1)
        r_oracle = 0.5 * (float(r_p2[b_star_1]) + float(r_p1[b_star_2]))

        return {
            "p1_orig": p1_orig,
            "p2_orig": p2_orig,
            "p1_pos": p1,
            "p2_pos": p2,
            "b_star_1": b_star_1,
            "b_dag_1": b_dag_1,
            "b_star_2": b_star_2,
            "b_dag_2": b_dag_2,
            "r_p1_56": r_p1,
            "r_p2_56": r_p2,
            "rho_ab_orth_p1": rho_p1,
            "rho_ab_orth_p2": rho_p2,
            "delta_1to2": d_1to2,
            "delta_2to1": d_2to1,
            "delta_cf": delta_cf,
            "r_cf_oracle": r_oracle,
        }

    s0 = _eval_split_point(int(seed))
    p1_pos = s0["p1_pos"]
    p2_pos = s0["p2_pos"]
    m1 = int(p1_pos.size)
    m2 = int(p2_pos.size)

    # Exact within-half bootstraps: P2 with `seed`, P1 with `seed + 10000`
    draws_p2 = np.random.default_rng(int(seed)).integers(0, m2, size=(int(n_boot), m2))
    draws_p1 = np.random.default_rng(int(seed) + HALF_P1_BOOT_SEED_OFFSET).integers(0, m1, size=(int(n_boot), m1))

    hm_p2 = _subset_moments(half_mom, p2_pos) if half_mom is not None else None
    hm_p1 = _subset_moments(half_mom, p1_pos) if half_mom is not None else None
    rho_p2, ceil_p2, b_ceil_p2 = _compute_half_ceiling_draws(hm_p2, draws_p2)
    rho_p1, ceil_p1, b_ceil_p1 = _compute_half_ceiling_draws(hm_p1, draws_p1)

    b_star_1 = int(s0["b_star_1"])
    b_dag_1 = int(s0["b_dag_1"])
    b_star_2 = int(s0["b_star_2"])
    b_dag_2 = int(s0["b_dag_2"])

    pt_p2_best, boot_p2_best, _ = _eval_pred_mom_with_ceiling(
        _subset_moments(moms_56[b_star_1], p2_pos), ceil_p2, b_ceil_p2, draws_p2, zero_null_to_zero=zero_null
    )
    pt_p2_worst, boot_p2_worst, _ = _eval_pred_mom_with_ceiling(
        _subset_moments(moms_56[b_dag_1], p2_pos), ceil_p2, b_ceil_p2, draws_p2, zero_null_to_zero=zero_null
    )
    pt_p1_best, boot_p1_best, _ = _eval_pred_mom_with_ceiling(
        _subset_moments(moms_56[b_star_2], p1_pos), ceil_p1, b_ceil_p1, draws_p1, zero_null_to_zero=zero_null
    )
    pt_p1_worst, boot_p1_worst, _ = _eval_pred_mom_with_ceiling(
        _subset_moments(moms_56[b_dag_2], p1_pos), ceil_p1, b_ceil_p1, draws_p1, zero_null_to_zero=zero_null
    )

    boot_d_1to2 = boot_p2_best - boot_p2_worst
    boot_d_2to1 = boot_p1_best - boot_p1_worst
    boot_delta_cf = 0.5 * (boot_d_1to2 + boot_d_2to1)
    boot_r_oracle = 0.5 * (boot_p2_best + boot_p1_best)

    pt_d_1to2 = float(pt_p2_best - pt_p2_worst)
    pt_d_2to1 = float(pt_p1_best - pt_p1_worst)
    pt_delta_cf = 0.5 * (pt_d_1to2 + pt_d_2to1)
    pt_r_oracle = 0.5 * (float(pt_p2_best) + float(pt_p1_best))

    sum_d_1to2 = summarize_boot_r17(pt_d_1to2, boot_d_1to2, m2, n_boot)
    sum_d_2to1 = summarize_boot_r17(pt_d_2to1, boot_d_2to1, m1, n_boot)
    sum_delta_cf = summarize_boot_r17(pt_delta_cf, boot_delta_cf, m, n_boot)
    sum_r_oracle = summarize_boot_r17(pt_r_oracle, boot_r_oracle, m, n_boot)

    # Positive-direction two-sided p-value for H-J1 Holm ordering
    fin_dcf = boot_delta_cf[np.isfinite(boot_delta_cf)]
    sum_delta_cf["p_pos_two_sided"] = _empirical_positive_two_sided_p(fin_dcf)

    # Out-of-sample contrasts vs Cross-Fitted Oracle on locked s=0 (`§3.4` line 249)
    pt_p2_l1, boot_p2_l1, _ = _eval_pred_mom_with_ceiling(
        _subset_moments(mom_l1, p2_pos), ceil_p2, b_ceil_p2, draws_p2, zero_null_to_zero=zero_null
    )
    pt_p1_l1, boot_p1_l1, _ = _eval_pred_mom_with_ceiling(
        _subset_moments(mom_l1, p1_pos), ceil_p1, b_ceil_p1, draws_p1, zero_null_to_zero=zero_null
    )
    pt_l1_vs_oracle = 0.5 * (float(pt_p2_l1) + float(pt_p1_l1)) - pt_r_oracle
    boot_l1_vs_oracle = 0.5 * (boot_p2_l1 + boot_p1_l1) - boot_r_oracle
    l1_vs_oracle_sum = summarize_boot_r17(pt_l1_vs_oracle, boot_l1_vs_oracle, m, n_boot)

    l2_vs_oracle_sum: Optional[Dict[str, Any]] = None
    if l2_selected_idx is not None:
        b_l2 = int(l2_selected_idx)
        pt_p2_l2, boot_p2_l2, _ = _eval_pred_mom_with_ceiling(
            _subset_moments(moms_56[b_l2], p2_pos), ceil_p2, b_ceil_p2, draws_p2, zero_null_to_zero=zero_null
        )
        pt_p1_l2, boot_p1_l2, _ = _eval_pred_mom_with_ceiling(
            _subset_moments(moms_56[b_l2], p1_pos), ceil_p1, b_ceil_p1, draws_p1, zero_null_to_zero=zero_null
        )
        pt_l2_vs_oracle = 0.5 * (float(pt_p2_l2) + float(pt_p1_l2)) - pt_r_oracle
        boot_l2_vs_oracle = 0.5 * (boot_p2_l2 + boot_p1_l2) - boot_r_oracle
        l2_vs_oracle_sum = summarize_boot_r17(pt_l2_vs_oracle, boot_l2_vs_oracle, m, n_boot)

    # Multi-split stability diagnostic (`s = 0..n_splits_diag-1`, point estimates only)
    diag_deltas: List[float] = [pt_delta_cf]
    diag_oracles: List[float] = [pt_r_oracle]
    diag_records: List[Dict[str, Any]] = [
        {
            "s": 0,
            "seed": int(seed),
            "b_star_1": b_star_1,
            "b_dag_1": b_dag_1,
            "b_star_2": b_star_2,
            "b_dag_2": b_dag_2,
            "delta_cf": pt_delta_cf,
            "r_cf_oracle": pt_r_oracle,
        }
    ]
    for s_idx in range(1, max(1, int(n_splits_diag))):
        sp_res = _eval_split_point(int(seed) + s_idx)
        diag_deltas.append(float(sp_res["delta_cf"]))
        diag_oracles.append(float(sp_res["r_cf_oracle"]))
        diag_records.append(
            {
                "s": int(s_idx),
                "seed": int(seed) + s_idx,
                "b_star_1": int(sp_res["b_star_1"]),
                "b_dag_1": int(sp_res["b_dag_1"]),
                "b_star_2": int(sp_res["b_star_2"]),
                "b_dag_2": int(sp_res["b_dag_2"]),
                "delta_cf": float(sp_res["delta_cf"]),
                "r_cf_oracle": float(sp_res["r_cf_oracle"]),
            }
        )

    fin_diag = np.asarray([x for x in diag_deltas if np.isfinite(x)], dtype=np.float64)
    fin_ora = np.asarray([x for x in diag_oracles if np.isfinite(x)], dtype=np.float64)
    multi_split_summary = {
        "n_splits": int(len(diag_records)),
        "n_finite_splits": int(fin_diag.size),
        "delta_cf_mean": float(np.mean(fin_diag)) if fin_diag.size else float("nan"),
        "delta_cf_median": float(np.median(fin_diag)) if fin_diag.size else float("nan"),
        "delta_cf_min": float(np.min(fin_diag)) if fin_diag.size else float("nan"),
        "delta_cf_max": float(np.max(fin_diag)) if fin_diag.size else float("nan"),
        "delta_cf_iqr": (
            float(np.percentile(fin_diag, 75.0) - np.percentile(fin_diag, 25.0)) if fin_diag.size else float("nan")
        ),
        "r_cf_oracle_mean": float(np.mean(fin_ora)) if fin_ora.size else float("nan"),
        "r_cf_oracle_median": float(np.median(fin_ora)) if fin_ora.size else float("nan"),
        "per_split": diag_records,
    }

    return {
        "split_s0": {
            "seed": int(seed),
            "p1_boot_seed": int(seed) + HALF_P1_BOOT_SEED_OFFSET,
            "p2_boot_seed": int(seed),
            "n_p1": m1,
            "n_p2": m2,
            "p1_indices": [int(i) for i in s0["p1_orig"]],
            "p2_indices": [int(i) for i in s0["p2_orig"]],
            "b_star_1": b_star_1,
            "b_dag_1": b_dag_1,
            "b_star_2": b_star_2,
            "b_dag_2": b_dag_2,
            "rho_ab_orth_p1": float(rho_p1),
            "rho_ab_orth_p2": float(rho_p2),
            "ceiling_p1": float(ceil_p1),
            "ceiling_p2": float(ceil_p2),
            "delta_1to2": sum_d_1to2,
            "delta_2to1": sum_d_2to1,
            "delta_cf_0": sum_delta_cf,
            "r_cf_oracle_0": sum_r_oracle,
            "l1_vs_cf_oracle": l1_vs_oracle_sum,
            "l2_vs_cf_oracle": l2_vs_oracle_sum,
            "half_valid_both_halves": bool(
                sum_d_1to2["valid_frac"] >= MIN_BOOT_VALID_FRAC
                and sum_d_2to1["valid_frac"] >= MIN_BOOT_VALID_FRAC
                and np.isfinite(rho_p1) and rho_p1 > 0.0
                and np.isfinite(rho_p2) and rho_p2 > 0.0
            ),
        },
        "multi_split_diag": multi_split_summary,
        # Raw bootstrap arrays (stripped before JSON serialization, used for exact Holm CI check)
        "_boot_delta_cf_0": boot_delta_cf,
        "_boot_r_cf_oracle_0": boot_r_oracle,
    }


# ==============================================================================
# 3. `L2_ntc_recon_match` Label-Free Selection (`PREREG_J §3.4`, Ruling 3)
# ==============================================================================


def _pearson_safe(a: np.ndarray, b: np.ndarray) -> float:
    x = np.asarray(a, dtype=np.float64)
    y = np.asarray(b, dtype=np.float64)
    if x.size < 2 or y.size != x.size:
        return float("nan")
    if not (np.all(np.isfinite(x)) and np.all(np.isfinite(y))):
        return float("nan")
    xc = x - float(np.mean(x))
    yc = y - float(np.mean(y))
    denom = float(np.linalg.norm(xc) * np.linalg.norm(yc))
    if denom <= 1e-12 or not np.isfinite(denom):
        return float("nan")
    return float(np.clip(float(np.dot(xc, yc)) / denom, -1.0, 1.0))


def _summarize_criterion_56(scores: np.ndarray) -> Dict[str, Any]:
    """Summarize a 56-index selection criterion with mandatory Ruling 3 top-1 / top-2 margin."""
    arr = np.asarray(scores, dtype=np.float64)
    if arr.shape != (N_BATCH_INDICES_STATE,):
        raise ValueError(f"Expected 56 criterion values, got {arr.shape}.")
    safe = np.where(np.isfinite(arr), arr, -np.inf)
    order = sorted(range(N_BATCH_INDICES_STATE), key=lambda b: (-safe[b], b))
    top1_idx = int(order[0])
    top2_idx = int(order[1])
    c1 = float(arr[top1_idx])
    c2 = float(arr[top2_idx])
    margin = float(c1 - c2) if (np.isfinite(c1) and np.isfinite(c2)) else float("nan")

    fin = arr[np.isfinite(arr)]
    q25 = float(np.percentile(fin, 25.0)) if fin.size else float("nan")
    q75 = float(np.percentile(fin, 75.0)) if fin.size else float("nan")
    return {
        "top1_idx": top1_idx,
        "criterion_top1": c1,
        "top2_idx": top2_idx,
        "criterion_top2": c2,
        "criterion_margin_top1_minus_top2": margin,
        "b0_criterion": float(arr[0]),
        "b0_rank_1based": int(order.index(0) + 1),
        "min": float(np.min(fin)) if fin.size else float("nan"),
        "q25": q25,
        "median": float(np.median(fin)) if fin.size else float("nan"),
        "q75": q75,
        "iqr": float(q75 - q25) if (np.isfinite(q25) and np.isfinite(q75)) else float("nan"),
        "max": float(np.max(fin)) if fin.size else float("nan"),
        "criterion_values_56": [float(x) for x in arr],
    }


def select_l2_ntc_recon_match_r17(
    ntc_pb_ensemble_56: np.ndarray,
    obs_ntc_mean_eval: np.ndarray,
    ntc_pb_per_fold_56: Optional[np.ndarray] = None,
    fold_names: Optional[Sequence[str]] = None,
) -> Dict[str, Any]:
    """
    `PREREG_J §3.4` (Ruling 3 Locked) Label-Free NTC Reconstruction Match:
      - Primary: `C_NTC(b) = Pearson(Y_hat_ntc_4fold[b], Y_obs_ctrl_official_median_log1p)` on `G_eval`.
      - Mandatory disclosures: `top1_idx`, `criterion_top1`, `top2_idx`, `criterion_top2`,
        `criterion_margin_top1_minus_top2`, 56-index summary (`min/median/IQR/max`), rank of `b=0`.
      - Secondary 1: Negative L2 distance `-||Y_hat_ntc_4fold[b] - Y_obs_ctrl||_2` with its own margin.
      - Secondary 2: Per-fold Pearson selection (`f = 1..4`) and cross-fold agreement.
    """
    pb_ens = np.asarray(ntc_pb_ensemble_56, dtype=np.float64)
    y_obs = np.asarray(obs_ntc_mean_eval, dtype=np.float64)
    if pb_ens.ndim != 2 or pb_ens.shape[0] != N_BATCH_INDICES_STATE or pb_ens.shape[1] != y_obs.size:
        raise ValueError(
            f"Incompatible NTC profile shapes: pb_ens={pb_ens.shape}, y_obs={y_obs.shape}."
        )

    pearson_scores = np.array([_pearson_safe(pb_ens[b], y_obs) for b in range(N_BATCH_INDICES_STATE)], dtype=np.float64)
    primary = _summarize_criterion_56(pearson_scores)

    neg_l2_scores = np.array(
        [
            -float(np.linalg.norm(pb_ens[b] - y_obs))
            if (np.all(np.isfinite(pb_ens[b])) and np.all(np.isfinite(y_obs)))
            else float("nan")
            for b in range(N_BATCH_INDICES_STATE)
        ],
        dtype=np.float64,
    )
    secondary_neg_l2 = _summarize_criterion_56(neg_l2_scores)

    per_fold_out: Dict[str, Any] = {}
    fold_agreement: Dict[str, Any] = {}
    if ntc_pb_per_fold_56 is not None:
        pb_pf = np.asarray(ntc_pb_per_fold_56, dtype=np.float64)
        n_f = int(pb_pf.shape[0])
        fnames = list(fold_names) if fold_names is not None else [f"fold_{i}" for i in range(n_f)]
        selected_indices: List[int] = []
        fold_score_vecs: List[np.ndarray] = []
        for fi in range(n_f):
            f_scores = np.array(
                [_pearson_safe(pb_pf[fi, b], y_obs) for b in range(N_BATCH_INDICES_STATE)],
                dtype=np.float64,
            )
            f_sum = _summarize_criterion_56(f_scores)
            per_fold_out[str(fnames[fi])] = f_sum
            selected_indices.append(int(f_sum["top1_idx"]))
            fold_score_vecs.append(f_scores)

        ens_top1 = int(primary["top1_idx"])
        pair_rhos: Dict[str, float] = {}
        for i in range(n_f):
            for j in range(i + 1, n_f):
                r_ij = spearmanr(fold_score_vecs[i], fold_score_vecs[j]).statistic
                pair_rhos[f"{fnames[i]}__vs__{fnames[j]}"] = float(r_ij) if np.isfinite(r_ij) else float("nan")
        fold_agreement = {
            "per_fold_selected_idx": {str(fn): int(si) for fn, si in zip(fnames, selected_indices)},
            "n_unique_selected": int(len(set(selected_indices))),
            "all_folds_agree": bool(len(set(selected_indices)) == 1),
            "n_folds_agree_with_ensemble": int(sum(1 for si in selected_indices if si == ens_top1)),
            "pairwise_fold_criterion_spearman": pair_rhos,
        }

    return {
        "selected_batch_idx": int(primary["top1_idx"]),
        "primary_pearson": primary,
        "secondary_neg_l2": secondary_neg_l2,
        "secondary_per_fold": per_fold_out,
        "secondary_per_fold_agreement": fold_agreement,
    }


# ==============================================================================
# 4. `H-J1`, `H-J2`, `H-J3`, `J2` Floor Guard, & Holm-Bonferroni (`§3.1`–`§3.4`)
# ==============================================================================


def evaluate_detectable_signal_r17(
    r_cf_oracle_0_sum: Mapping[str, Any],
    l1_ens56_sum: Mapping[str, Any],
) -> Dict[str, Any]:
    """
    `PREREG_J §3.2` Item 4 (`J2` Floor-Effect Guard):
      `has_detectable_signal <=> (ci9833_lo(r_CF_Oracle,0) > 0 AND vf >= 0.95)
                                 OR (ci9833_lo(r(L1_ens56)) > 0 AND vf >= 0.95)`.
    """
    ora_vf = float(r_cf_oracle_0_sum.get("valid_frac", 0.0))
    ora_lo = float(r_cf_oracle_0_sum.get("ci9833_lo", float("nan")))
    ora_pos = bool(ora_vf >= MIN_BOOT_VALID_FRAC and np.isfinite(ora_lo) and ora_lo > 0.0)

    l1_vf = float(l1_ens56_sum.get("valid_frac", 0.0))
    l1_lo = float(l1_ens56_sum.get("ci9833_lo", float("nan")))
    l1_pos = bool(l1_vf >= MIN_BOOT_VALID_FRAC and np.isfinite(l1_lo) and l1_lo > 0.0)

    return {
        "has_detectable_signal": bool(ora_pos or l1_pos),
        "cf_oracle_positive_ci9833": ora_pos,
        "l1_ens56_positive_ci9833": l1_pos,
        "cf_oracle_ci9833_lo": ora_lo,
        "l1_ens56_ci9833_lo": l1_lo,
    }


def apply_holm_bonferroni_r17(
    p_values: Mapping[str, float],
    alpha_family: float = ALPHA_FAMILY,
) -> Dict[str, Dict[str, Any]]:
    """
    Exact step-down Holm-Bonferroni FWER <= `alpha_family` procedure over `M = len(p_values)`.
    Deterministic ordering by `(p_value, test_key)`; NaN p-values treated as `1.0`.
    Returns per-key `{raw_p, rank_1based, m_tests, step_alpha, holm_adj_p, rejected, effective_alpha}`.
    """
    keys = sorted(p_values.keys())
    m = len(keys)
    if m == 0:
        return {}

    clean_p = {
        k: (float(p_values[k]) if (p_values[k] is not None and np.isfinite(float(p_values[k]))) else 1.0)
        for k in keys
    }
    ordered = sorted(keys, key=lambda k: (clean_p[k], k))

    out: Dict[str, Dict[str, Any]] = {}
    still_rejecting = True
    running_adj_p = 0.0
    stop_alpha = float(alpha_family) / float(m)

    for rank_0, k in enumerate(ordered):
        denom = float(m - rank_0)
        step_alpha = float(alpha_family) / denom
        p_raw = clean_p[k]
        adj_p = min(1.0, max(running_adj_p, p_raw * denom))
        running_adj_p = adj_p

        if still_rejecting and p_raw < step_alpha:
            rej = True
            eff_alpha = step_alpha
        else:
            if still_rejecting:
                stop_alpha = step_alpha
            still_rejecting = False
            rej = False
            eff_alpha = stop_alpha

        out[k] = {
            "raw_p": float(p_raw),
            "rank_1based": int(rank_0 + 1),
            "m_tests": int(m),
            "step_alpha": float(step_alpha),
            "effective_alpha": float(eff_alpha),
            "holm_adj_p": float(adj_p),
            "rejected": bool(rej),
        }
    return out


def adjudicate_h_j1_r17(
    var_comp: Mapping[str, Any],
    cf_s0: Mapping[str, Any],
    boot_delta_cf: np.ndarray,
    detectable_signal: Mapping[str, Any],
    landscape_summary: Mapping[str, Any],
    l1_ens56_sum: Mapping[str, Any],
    rho_ab_orth_full: float,
    alpha_adj: float,
    holm_rejected: Optional[bool] = None,
    delta_practical: float = DEFAULT_PRACTICAL_DELTA,
    tau_threshold: float = DEFAULT_TAU_THRESHOLD,
) -> Dict[str, Any]:
    """
    `PREREG_J §3.2` Item 5 Locked Decision Rule for `H-J1` (with `J2` Floor-Effect Guard):
      1. `INCONCLUSIVE` (`VALIDITY_FAIL`): `valid_frac < 0.95` on any of the 56 indices or
         either cross-fit half, OR `rho_ab_orth <= 0`, OR non-finite required estimates/CIs.
      2. Statistically detected: `p_var < alpha_adj` AND `ci_{1 - alpha_adj}_lo(delta_cf_0) > 0.0`
         (and `holm_rejected is not False` when inside Holm family):
         - `BATCH_DEPENDENCE_PRACTICALLY_MEANINGFUL` iff `delta_cf_0 > delta` AND `2 * tau_batch > delta` (`tau_batch > 0.02`);
         - `BATCH_DEPENDENCE_DETECTED_SUB_DELTA` otherwise.
      3. `INCONCLUSIVE_LOW_SIGNAL` (`J2`): NOT detected AND `has_detectable_signal == False`.
      4. `BATCH_INVARIANT_WITHIN_DELTA`: NOT detected AND `has_detectable_signal == True`
         AND `ci95_hi(delta_cf_0) < delta` AND `tau_batch <= 0.02`.
      5. `INCONCLUSIVE` (`INSUFFICIENT_PRECISION`): all remaining cases.
    """
    dcf_sum = cf_s0["delta_cf_0"]
    ora_sum = cf_s0["r_cf_oracle_0"]

    p_var = float(var_comp.get("p_var", float("nan")))
    tau_batch = float(var_comp.get("tau_batch", float("nan")))
    min_vf56 = float(var_comp.get("min_per_index_valid_frac", 0.0))
    joint_vf56 = float(var_comp.get("joint_valid_frac", 0.0))

    dcf_est = float(dcf_sum.get("estimate", float("nan")))
    dcf_vf = float(dcf_sum.get("valid_frac", 0.0))
    dcf_ci95_hi = float(dcf_sum.get("ci95_hi", float("nan")))
    halves_ok = bool(cf_s0.get("half_valid_both_halves", False))

    fin_dcf = np.asarray(boot_delta_cf, dtype=np.float64)
    fin_dcf = fin_dcf[np.isfinite(fin_dcf)]
    ci_adj_lo, ci_adj_hi = _two_sided_ci_from_vals(fin_dcf, float(alpha_adj))

    validity_ok = bool(
        min_vf56 >= MIN_BOOT_VALID_FRAC
        and joint_vf56 >= MIN_BOOT_VALID_FRAC
        and dcf_vf >= MIN_BOOT_VALID_FRAC
        and halves_ok
        and np.isfinite(rho_ab_orth_full) and rho_ab_orth_full > 0.0
        and np.isfinite(p_var)
        and np.isfinite(tau_batch)
        and np.isfinite(dcf_est)
        and np.isfinite(ci_adj_lo) and np.isfinite(ci_adj_hi)
        and np.isfinite(dcf_ci95_hi)
    )

    p_dcf_pos = float(dcf_sum.get("p_pos_two_sided", float("nan")))
    p_iut = float(max(p_var, p_dcf_pos)) if (np.isfinite(p_var) and np.isfinite(p_dcf_pos)) else 1.0

    ci_and_p_detected = bool(validity_ok and p_var < float(alpha_adj) and ci_adj_lo > 0.0)
    if holm_rejected is None:
        statistically_detected = ci_and_p_detected
    else:
        # Require BOTH Holm step-down rejection AND exact two-sided (1 - alpha_adj) CI + p_var verification
        statistically_detected = bool(holm_rejected and ci_and_p_detected)

    has_sig = bool(detectable_signal.get("has_detectable_signal", False))

    subreason = ""
    if not validity_ok:
        verdict = "INCONCLUSIVE"
        subreason = "VALIDITY_FAIL"
    elif statistically_detected:
        if dcf_est > float(delta_practical) and tau_batch > float(tau_threshold):
            verdict = "BATCH_DEPENDENCE_PRACTICALLY_MEANINGFUL"
        else:
            verdict = "BATCH_DEPENDENCE_DETECTED_SUB_DELTA"
    elif not has_sig:
        verdict = "INCONCLUSIVE_LOW_SIGNAL"
        subreason = "FLOOR_EFFECT_NO_DETECTABLE_SIGNAL"
    elif dcf_ci95_hi < float(delta_practical) and tau_batch <= float(tau_threshold):
        verdict = "BATCH_INVARIANT_WITHIN_DELTA"
    else:
        verdict = "INCONCLUSIVE"
        subreason = "INSUFFICIENT_PRECISION"

    return {
        "verdict": verdict,
        "subreason": subreason,
        "validity_ok": validity_ok,
        "alpha_adj": float(alpha_adj),
        "holm_rejected": holm_rejected,
        "statistically_detected": statistically_detected,
        "p_var": p_var,
        "p_delta_cf_pos_two_sided": p_dcf_pos,
        "p_h_j1_iut": p_iut,
        "delta_cf_0_estimate": dcf_est,
        "delta_cf_0_ci_adj_lo": ci_adj_lo,
        "delta_cf_0_ci_adj_hi": ci_adj_hi,
        "delta_cf_0_ci95_hi": dcf_ci95_hi,
        "tau_batch": tau_batch,
        "two_tau_batch": 2.0 * tau_batch if np.isfinite(tau_batch) else float("nan"),
        "has_detectable_signal": has_sig,
        # Mandatory J2 reporting alongside every H-J1 call
        "mandatory_j2_context": {
            "max_b_r_orth_cal": float(landscape_summary["max"]),
            "argmax_b_idx": int(landscape_summary["argmax_idx"]),
            "l1_ens56_r_orth_cal": float(l1_ens56_sum["estimate"]),
            "l1_ens56_ci9833_lo": float(l1_ens56_sum["ci9833_lo"]),
            "l1_ens56_ci9833_hi": float(l1_ens56_sum["ci9833_hi"]),
            "cf_oracle_0_r_orth_cal": float(ora_sum["estimate"]),
            "cf_oracle_0_ci9833_lo": float(ora_sum["ci9833_lo"]),
            "cf_oracle_0_ci9833_hi": float(ora_sum["ci9833_hi"]),
        },
    }


def _classify_vs_fct_from_ci(est: float, lo: float, hi: float, vf: float) -> str:
    if vf < MIN_BOOT_VALID_FRAC or not (np.isfinite(est) and np.isfinite(lo) and np.isfinite(hi)):
        return "INDISTINGUISHABLE"
    if lo > 0.0:
        return "MODEL_SUPERIOR_TO_FCT"
    if hi < 0.0:
        return "FCT_SUPERIOR_TO_MODEL"
    return "INDISTINGUISHABLE"


def _summarize_flip_calls(calls_56: Sequence[str], has_detectable_signal: bool) -> Dict[str, Any]:
    c0 = str(calls_56[0])
    idx_mod = [b for b, c in enumerate(calls_56) if c == "MODEL_SUPERIOR_TO_FCT"]
    idx_ind = [b for b, c in enumerate(calls_56) if c == "INDISTINGUISHABLE"]
    idx_fct = [b for b, c in enumerate(calls_56) if c == "FCT_SUPERIOR_TO_MODEL"]
    flip_idx = [b for b, c in enumerate(calls_56) if c != c0]
    n_flip = len(flip_idx)
    f_flip = float(n_flip) / float(N_BATCH_INDICES_STATE)
    has_full_rev = bool(len(idx_mod) >= 1 and len(idx_fct) >= 1)

    if has_full_rev:
        cls = "CALL_UNSTABLE_SIGN_REVERSAL"
    elif n_flip >= 5:
        cls = "CALL_UNSTABLE_SUBSTANTIAL"
    elif n_flip >= 1:
        cls = "CALL_UNSTABLE_MINOR"
    elif has_detectable_signal:
        cls = "CALL_STABLE"
    else:
        cls = "INCONCLUSIVE_LOW_SIGNAL"

    return {
        "classification": cls,
        "baseline_b0_call": c0,
        "n_flip": int(n_flip),
        "f_flip": f_flip,
        "has_full_reversal": has_full_rev,
        "counts": {
            "MODEL_SUPERIOR_TO_FCT": len(idx_mod),
            "INDISTINGUISHABLE": len(idx_ind),
            "FCT_SUPERIOR_TO_MODEL": len(idx_fct),
        },
        "indices_by_call": {
            "MODEL_SUPERIOR_TO_FCT": idx_mod,
            "INDISTINGUISHABLE": idx_ind,
            "FCT_SUPERIOR_TO_MODEL": idx_fct,
        },
        "flip_indices": flip_idx,
        "calls_56": list(calls_56),
    }


def adjudicate_h_j2_r17(
    r_hat_56: np.ndarray,
    r_boot_56: np.ndarray,
    fct_pt: float,
    fct_boot: np.ndarray,
    fct_summary: Mapping[str, Any],
    moments_bundle: Mapping[str, Any],
    detectable_signal: Mapping[str, Any],
    landscape_summary: Mapping[str, Any],
    l1_ens56_sum: Mapping[str, Any],
    n_boot: int = DEFAULT_R17_N_BOOT,
    n_boot_bonf56: int = DEFAULT_R17_N_BOOT_BONF56,
    seed: int = DEFAULT_R17_SEED,
) -> Dict[str, Any]:
    """
    `PREREG_J §3.3` Hypothesis `H-J2` — `STATE - FCT` Call Instability Across All 56 Indices:
      - Primary: exact `PREREG_I §3.2` `97.5%` CI call at every `b in 0..55` (via frozen
        `adjudicate_claim_i2_line`), `f_flip`, `has_full_reversal`, 5-way audit class with `J2`.
      - Companion Bonferroni-56 Hardened Stability Check (`alpha = 0.025 / 56`) side-by-side:
          (a) Normal-calibrated `z = norm.ppf(1 - 0.025 / 112) = 3.5117` using 4,000-boot `boot_se`;
          (b) `n_boot_bonf56 = 20,000` empirical quantile CIs (`[0.02232%, 99.97768%]`).
    """
    has_sig = bool(detectable_signal.get("has_detectable_signal", False))
    m = int(moments_bundle["m"])
    zero_null = bool(moments_bundle["double_orth"])

    calls_primary: List[str] = []
    calls_bonf56_norm: List[str] = []
    per_index_vs_fct: List[Dict[str, Any]] = []

    for b in range(N_BATCH_INDICES_STATE):
        diff_pt = float(r_hat_56[b] - fct_pt)
        diff_boot = r_boot_56[:, b] - fct_boot
        d_sum = summarize_boot_r17(diff_pt, diff_boot, m, n_boot)
        adj_i2 = adjudicate_claim_i2_line(d_sum)
        calls_primary.append(str(adj_i2["classification"]))

        # Companion (a): normal-calibrated Bonferroni-56 (`z = 3.5117`)
        se_b = float(d_sum["boot_se"])
        lo_n = float(diff_pt - Z_BONF56 * se_b) if (np.isfinite(diff_pt) and np.isfinite(se_b)) else float("nan")
        hi_n = float(diff_pt + Z_BONF56 * se_b) if (np.isfinite(diff_pt) and np.isfinite(se_b)) else float("nan")
        cls_n = _classify_vs_fct_from_ci(diff_pt, lo_n, hi_n, float(d_sum["valid_frac"]))
        calls_bonf56_norm.append(cls_n)

        per_index_vs_fct.append(
            {
                "batch_idx": int(b),
                "estimate": diff_pt,
                "ci975_lo": float(d_sum["ci975_lo"]),
                "ci975_hi": float(d_sum["ci975_hi"]),
                "boot_se": se_b,
                "valid_frac": float(d_sum["valid_frac"]),
                "call_975": str(adj_i2["classification"]),
                "bonf56_norm_lo": lo_n,
                "bonf56_norm_hi": hi_n,
                "call_bonf56_norm": cls_n,
            }
        )

    # Companion (b): `n_boot = 20,000` empirical quantile Bonferroni-56
    draws_20k = np.random.default_rng(int(seed)).integers(0, m, size=(int(n_boot_bonf56), m))
    _, ceil_20k, b_ceil_20k = _compute_half_ceiling_draws(moments_bundle["half_mom"], draws_20k)
    _, fct_boot_20k, _ = _eval_pred_mom_with_ceiling(
        moments_bundle["mom_fct"], ceil_20k, b_ceil_20k, draws_20k, zero_null_to_zero=zero_null
    )
    calls_bonf56_emp: List[str] = []
    for b in range(N_BATCH_INDICES_STATE):
        pt_b, boot_b_20k, _ = _eval_pred_mom_with_ceiling(
            moments_bundle["moms_56"][b], ceil_20k, b_ceil_20k, draws_20k, zero_null_to_zero=zero_null
        )
        diff_20k = boot_b_20k - fct_boot_20k
        vm = np.isfinite(diff_20k)
        vf_20k = float(np.mean(vm)) if diff_20k.size else 0.0
        lo_e, hi_e = _two_sided_ci_from_vals(diff_20k[vm], ALPHA_BONF56)
        cls_e = _classify_vs_fct_from_ci(float(pt_b - fct_pt), lo_e, hi_e, vf_20k)
        calls_bonf56_emp.append(cls_e)
        per_index_vs_fct[b]["bonf56_emp20k_lo"] = lo_e
        per_index_vs_fct[b]["bonf56_emp20k_hi"] = hi_e
        per_index_vs_fct[b]["bonf56_emp20k_valid_frac"] = vf_20k
        per_index_vs_fct[b]["call_bonf56_emp20k"] = cls_e

    prim_sum = _summarize_flip_calls(calls_primary, has_sig)
    norm_sum = _summarize_flip_calls(calls_bonf56_norm, has_sig)
    emp_sum = _summarize_flip_calls(calls_bonf56_emp, has_sig)

    return {
        **prim_sum,
        "has_detectable_signal": has_sig,
        "all_56_valid_boot": bool(all(x["valid_frac"] >= MIN_BOOT_VALID_FRAC for x in per_index_vs_fct)),
        "companion_bonf56_normal": {
            "alpha": float(ALPHA_BONF56),
            "z_critical": float(Z_BONF56),
            **norm_sum,
        },
        "companion_bonf56_emp20k": {
            "alpha": float(ALPHA_BONF56),
            "n_boot": int(n_boot_bonf56),
            "seed": int(seed),
            **emp_sum,
        },
        "mandatory_context": {
            "fct_r_orth_cal": float(fct_summary["estimate"]),
            "fct_ci975_lo": float(fct_summary["ci975_lo"]),
            "fct_ci975_hi": float(fct_summary["ci975_hi"]),
            "max_b_r_orth_cal": float(landscape_summary["max"]),
            "l1_ens56_r_orth_cal": float(l1_ens56_sum["estimate"]),
        },
        "per_index_vs_fct": per_index_vs_fct,
    }


def adjudicate_h_j3_rule_r17(
    rule_name: str,
    diff_vs_b0_pt: float,
    diff_vs_b0_boot: np.ndarray,
    rule_vs_fct_sum: Mapping[str, Any],
    rule_vs_cf_oracle_sum: Optional[Mapping[str, Any]],
    alpha_adj: float,
    m_items: int,
    n_boot: int,
    holm_rejected: Optional[bool] = None,
    exact_tie_selected_zero: bool = False,
    delta_practical: float = DEFAULT_PRACTICAL_DELTA,
) -> Dict[str, Any]:
    """
    `PREREG_J §3.4` Item 3 Decision Rule for `L in {L1_ens56, L2_ntc_recon_match}` vs `b = 0`:
      - `EXACT_TIE_SELECTED_ZERO` if `L == L2_ntc_recon_match` and `b_ntc == 0`;
      - `RULE_SUPERIOR_TO_ZERO` iff `valid_frac >= 0.95` AND `ci_adj_lo > 0.0` (and Holm rejected);
        annotated `PRACTICAL_GT_DELTA` (`> 0.04`) vs `SUB_DELTA_LE_DELTA` (`<= 0.04`);
      - `ZERO_SUPERIOR_TO_RULE` iff `valid_frac >= 0.95` AND `ci_adj_hi < 0.0` (and Holm rejected);
      - `INDISTINGUISHABLE` otherwise.
    """
    base_sum = summarize_boot_r17(
        diff_vs_b0_pt,
        diff_vs_b0_boot,
        m_items,
        n_boot,
        extra_alphas={"ci_adj": float(alpha_adj)},
    )
    vf = float(base_sum["valid_frac"])
    lo_adj = float(base_sum["ci_adj_lo"])
    hi_adj = float(base_sum["ci_adj_hi"])
    est = float(base_sum["estimate"])
    valid_boot = bool(vf >= MIN_BOOT_VALID_FRAC and np.isfinite(est) and np.isfinite(lo_adj) and np.isfinite(hi_adj))
    allow_reject = True if holm_rejected is None else bool(holm_rejected)

    mag_note = ""
    if exact_tie_selected_zero:
        verdict = "EXACT_TIE_SELECTED_ZERO"
    elif valid_boot and allow_reject and lo_adj > 0.0:
        verdict = "RULE_SUPERIOR_TO_ZERO"
        mag_note = "PRACTICAL_GT_DELTA" if est > float(delta_practical) else "SUB_DELTA_LE_DELTA"
    elif valid_boot and allow_reject and hi_adj < 0.0:
        verdict = "ZERO_SUPERIOR_TO_RULE"
        mag_note = "PRACTICAL_GT_DELTA" if (-est) > float(delta_practical) else "SUB_DELTA_LE_DELTA"
    else:
        verdict = "INDISTINGUISHABLE"

    return {
        "rule_name": str(rule_name),
        "verdict": verdict,
        "magnitude_annotation": mag_note,
        "alpha_adj": float(alpha_adj),
        "holm_rejected": holm_rejected,
        "vs_b0": base_sum,
        "vs_cf_oracle_s0": dict(rule_vs_cf_oracle_sum) if rule_vs_cf_oracle_sum is not None else None,
        "vs_fct": dict(rule_vs_fct_sum),
        "vs_fct_h_i2_call": adjudicate_claim_i2_line(dict(rule_vs_fct_sum)),
    }


# ==============================================================================
# 5. `J5a` Cross-Fitted `G4` 2x2 Decomposition (`PREREG_J §3.5` Item 1)
# ==============================================================================


def _summarize_g4_strata(
    g4_pairs: Mapping[str, Any],
    prefix: str,
    seen_sources: Sequence[str],
    unseen_sources: Sequence[str],
    same_lineage_sources: Sequence[str] = (),
) -> Dict[str, Any]:
    strata_defs = [
        ("seen_by_checkpoint", sorted(seen_sources)),
        ("unseen_by_checkpoint", sorted(unseen_sources)),
    ]
    if same_lineage_sources:
        strata_defs.append(("same_lineage_resting_comparator", sorted(same_lineage_sources)))

    out: Dict[str, Any] = {}
    for st_name, st_list in strata_defs:
        matched_sks = [sk for sk in st_list if f"{prefix}{sk}" in g4_pairs]
        if not matched_sks:
            continue
        st_metrics: Dict[str, Any] = {"n_sources": len(matched_sks), "sources": matched_sks}
        for mkey in ("r_orth_cal", "r_dorth_cal"):
            m_out: Dict[str, Any] = {}
            for fld in ("delta_batch", "delta_ctrl", "delta_both", "delta_interaction_point", "batch_share_point"):
                vals: List[float] = []
                for sk in matched_sks:
                    entry = g4_pairs[f"{prefix}{sk}"][mkey][fld]
                    v = float(entry["estimate"]) if isinstance(entry, dict) else float(entry)
                    if np.isfinite(v):
                        vals.append(v)
                m_out[fld] = {
                    "n_finite": len(vals),
                    "median": float(np.median(vals)) if vals else float("nan"),
                    "mean": float(np.mean(vals)) if vals else float("nan"),
                    "min": float(min(vals)) if vals else float("nan"),
                    "max": float(max(vals)) if vals else float("nan"),
                }
            st_metrics[mkey] = m_out
        out[st_name] = st_metrics
    return out


def compute_j5a_g4_r17(
    delta_obs: np.ndarray,
    delta_half_a: np.ndarray,
    delta_half_b: np.ndarray,
    half_valid: np.ndarray,
    delta_transfer: np.ndarray,
    source_grand_mean: np.ndarray,
    target_indices: np.ndarray,
    delta_pred_56_matched: np.ndarray,
    mismatched_by_source_and_idx: Mapping[str, Mapping[int, np.ndarray]],
    seen_sources: Sequence[str],
    unseen_sources: Sequence[str],
    same_lineage_sources: Sequence[str],
    p1_indices: Sequence[int],
    p2_indices: Sequence[int],
    b_star_1: int,
    b_dag_1: int,
    b_star_2: Optional[int] = None,
    b_dag_2: Optional[int] = None,
    n_boot: int = DEFAULT_R17_N_BOOT,
    seed: int = DEFAULT_R17_SEED,
) -> Dict[str, Any]:
    """
    `PREREG_J §3.5` Item 1 (`J5a`) — Cross-Fitted `G4` 2x2 (`X_ctrl x batch`) Decomposition
    plus Pre-Fixed Continuity Pairs `(0, 17)` and `(0, 35)`:
      - Primary `J5a` readout: `(B_best = b*_1, B_worst = b_dag_1)` (and reversed `(b_dag_1, b*_1)`)
        evaluated STRICTLY on disjoint `P2` only (`seed = 20261701`).
      - Optional auxiliary diagnostic: `(b*_2, b_dag_2)` (and reversed) evaluated strictly on `P1`
        (`seed = 20261701 + 10000`).
      - Pre-fixed continuity pairs on ALL `N_t` perturbations (`seed = 20261701`):
        `(0, 17)` + `(17, 0)` and `(0, 35)` + `(35, 0)`.
    """
    p1 = np.asarray(p1_indices, dtype=np.int64)
    p2 = np.asarray(p2_indices, dtype=np.int64)
    all_sources = sorted(set(seen_sources) | set(unseen_sources) | set(same_lineage_sources))

    def _run_pair_set(
        b_m: int,
        b_s: int,
        row_sel: Optional[np.ndarray],
        pair_tag: str,
        boot_seed: int,
    ) -> Tuple[Dict[str, Any], Dict[str, Any]]:
        if row_sel is not None:
            obs_s = delta_obs[row_sel]
            ha_s = delta_half_a[row_sel]
            hb_s = delta_half_b[row_sel]
            hv_s = half_valid[row_sel]
            tm_s = delta_transfer[row_sel]
            ti_s = target_indices[row_sel]
            m_bm = delta_pred_56_matched[b_m][row_sel]
            m_bs = delta_pred_56_matched[b_s][row_sel]
        else:
            obs_s, ha_s, hb_s, hv_s, tm_s, ti_s = (
                delta_obs, delta_half_a, delta_half_b, half_valid, delta_transfer, target_indices
            )
            m_bm = delta_pred_56_matched[b_m]
            m_bs = delta_pred_56_matched[b_s]

        fwd_prefix = f"{pair_tag}_fwd_mis_"
        rev_prefix = f"{pair_tag}_rev_mis_"
        pairs_out: Dict[str, Any] = {}
        for sk in all_sources:
            s_map = mismatched_by_source_and_idx.get(sk, {})
            if b_m not in s_map or b_s not in s_map:
                continue
            s_bm = s_map[b_m][row_sel] if row_sel is not None else s_map[b_m]
            s_bs = s_map[b_s][row_sel] if row_sel is not None else s_map[b_s]

            pairs_out[f"{fwd_prefix}{sk}"] = compute_g4_factorial_decomposition_r16(
                delta_obs=obs_s,
                pred_m_bm=m_bm,
                pred_m_bs=m_bs,
                pred_s_bm=s_bm,
                pred_s_bs=s_bs,
                delta_transfer=tm_s,
                target_indices=ti_s,
                source_grand_mean=source_grand_mean,
                delta_half_a=ha_s,
                delta_half_b=hb_s,
                half_valid=hv_s,
                n_boot=n_boot,
                seed=boot_seed,
            )
            pairs_out[f"{rev_prefix}{sk}"] = compute_g4_factorial_decomposition_r16(
                delta_obs=obs_s,
                pred_m_bm=m_bs,
                pred_m_bs=m_bm,
                pred_s_bm=s_bs,
                pred_s_bs=s_bm,
                delta_transfer=tm_s,
                target_indices=ti_s,
                source_grand_mean=source_grand_mean,
                delta_half_a=ha_s,
                delta_half_b=hb_s,
                half_valid=hv_s,
                n_boot=n_boot,
                seed=boot_seed,
            )
        strata_out = {
            "forward": _summarize_g4_strata(pairs_out, fwd_prefix, seen_sources, unseen_sources, same_lineage_sources),
            "reversed": _summarize_g4_strata(pairs_out, rev_prefix, seen_sources, unseen_sources, same_lineage_sources),
        }
        return pairs_out, strata_out

    # 1. Locked primary cross-fitted J5a readout: (b*_1, b_dag_1) on disjoint P2 ONLY
    cf_p2_pairs, cf_p2_strata = _run_pair_set(int(b_star_1), int(b_dag_1), p2, "cf_p1_to_p2", int(seed))
    out: Dict[str, Any] = {
        "cross_fitted_primary_p1_to_p2": {
            "b_best_from_p1": int(b_star_1),
            "b_worst_from_p1": int(b_dag_1),
            "evaluated_on": "P2_only",
            "n_eval_perts": int(p2.size),
            "boot_seed": int(seed),
            "strata_summary": cf_p2_strata,
            "per_source": cf_p2_pairs,
        }
    }

    # 2. Optional auxiliary symmetric P2 -> P1 diagnostic
    if b_star_2 is not None and b_dag_2 is not None:
        has_all_p1 = any(
            (int(b_star_2) in mismatched_by_source_and_idx.get(sk, {}) and int(b_dag_2) in mismatched_by_source_and_idx.get(sk, {}))
            for sk in all_sources
        )
        if has_all_p1:
            cf_p1_pairs, cf_p1_strata = _run_pair_set(
                int(b_star_2), int(b_dag_2), p1, "cf_p2_to_p1", int(seed) + HALF_P1_BOOT_SEED_OFFSET
            )
            out["cross_fitted_auxiliary_p2_to_p1"] = {
                "b_best_from_p2": int(b_star_2),
                "b_worst_from_p2": int(b_dag_2),
                "evaluated_on": "P1_only",
                "n_eval_perts": int(p1.size),
                "boot_seed": int(seed) + HALF_P1_BOOT_SEED_OFFSET,
                "strata_summary": cf_p1_strata,
                "per_source": cf_p1_pairs,
            }

    # 3. Pre-fixed continuity pairs (0, 17) and (0, 35) on ALL N_t perturbations
    c017_pairs, c017_strata = _run_pair_set(0, 17, None, "continuity_b0_vs_b17", int(seed))
    c035_pairs, c035_strata = _run_pair_set(0, 35, None, "continuity_b0_vs_b35", int(seed))
    out["continuity_b0_vs_b17"] = {
        "b_m": 0,
        "b_s": 17,
        "evaluated_on": "all_N_eval",
        "n_eval_perts": int(delta_obs.shape[0]),
        "boot_seed": int(seed),
        "strata_summary": c017_strata,
        "per_source": c017_pairs,
    }
    out["continuity_b0_vs_b35"] = {
        "b_m": 0,
        "b_s": 35,
        "evaluated_on": "all_N_eval",
        "n_eval_perts": int(delta_obs.shape[0]),
        "boot_seed": int(seed),
        "strata_summary": c035_strata,
        "per_source": c035_pairs,
    }
    return out


# ==============================================================================
# 6. `J5b` Cross-Line & Cross-Model 56-Index Spearman Rank Correlation (`§3.5`)
# ==============================================================================


def _bh_fdr(p_vals: Sequence[float]) -> List[float]:
    """Standard Benjamini-Hochberg FDR q-values preserving input order."""
    m = len(p_vals)
    if m == 0:
        return []
    clean = [float(p) if (p is not None and np.isfinite(float(p))) else 1.0 for p in p_vals]
    order = sorted(range(m), key=lambda i: (clean[i], i))
    q_ordered = [0.0] * m
    running_min = 1.0
    for rev_i in range(m - 1, -1, -1):
        rank_1 = rev_i + 1
        orig_i = order[rev_i]
        q = min(1.0, clean[orig_i] * float(m) / float(rank_1))
        running_min = min(running_min, q)
        q_ordered[orig_i] = running_min
    return q_ordered


def _build_index_perm_table_r17(n_perm: int = DEFAULT_R17_N_PERM, seed: int = DEFAULT_R17_SEED) -> np.ndarray:
    """Shared `(n_perm, 56)` permutation table so every pair's p_perm is invariant to pair ordering."""
    rng = np.random.default_rng(int(seed))
    perms = np.empty((int(n_perm), N_BATCH_INDICES_STATE), dtype=np.int64)
    for k in range(int(n_perm)):
        perms[k] = rng.permutation(N_BATCH_INDICES_STATE)
    return perms


def compute_spearman_perm_pair_r17(
    vec_a: np.ndarray,
    vec_b: np.ndarray,
    perm_table: Optional[np.ndarray] = None,
    n_perm: int = DEFAULT_R17_N_PERM,
    seed: int = DEFAULT_R17_SEED,
) -> Dict[str, Any]:
    """
    `PREREG_J §3.5` Item 2 (`J5b`) exact 56-index Spearman rank correlation +
    `10,000`-permutation index-label null + top-5 / bottom-5 shared index overlaps.
    """
    a = np.asarray(vec_a, dtype=np.float64)
    b = np.asarray(vec_b, dtype=np.float64)
    if a.shape != (N_BATCH_INDICES_STATE,) or b.shape != (N_BATCH_INDICES_STATE,):
        raise ValueError(f"Expected 56-length vectors, got {a.shape} and {b.shape}.")

    if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))):
        return {
            "spearman_rho": float("nan"),
            "p_perm_two_sided": float("nan"),
            "n_perm": int(n_perm if perm_table is None else perm_table.shape[0]),
            "top5_a": [],
            "top5_b": [],
            "top5_overlap_count": 0,
            "top5_overlap_indices": [],
            "bottom5_a": [],
            "bottom5_b": [],
            "bottom5_overlap_count": 0,
            "bottom5_overlap_indices": [],
        }

    # Exact deterministic top-5 / bottom-5 (lowest-index tie-break)
    ord_a_desc = sorted(range(N_BATCH_INDICES_STATE), key=lambda i: (-a[i], i))
    ord_b_desc = sorted(range(N_BATCH_INDICES_STATE), key=lambda i: (-b[i], i))
    ord_a_asc = sorted(range(N_BATCH_INDICES_STATE), key=lambda i: (a[i], i))
    ord_b_asc = sorted(range(N_BATCH_INDICES_STATE), key=lambda i: (b[i], i))
    top5_a, top5_b = ord_a_desc[:5], ord_b_desc[:5]
    bot5_a, bot5_b = ord_a_asc[:5], ord_b_asc[:5]
    top5_ov = sorted(set(top5_a) & set(top5_b))
    bot5_ov = sorted(set(bot5_a) & set(bot5_b))

    ra = rankdata(a, method="average")
    rb = rankdata(b, method="average")
    ra_c = ra - float(np.mean(ra))
    rb_c = rb - float(np.mean(rb))
    na = float(np.linalg.norm(ra_c))
    nb = float(np.linalg.norm(rb_c))

    perms = perm_table if perm_table is not None else _build_index_perm_table_r17(n_perm=n_perm, seed=seed)
    n_p = int(perms.shape[0])

    if na <= 1e-12 or nb <= 1e-12:
        return {
            "spearman_rho": float("nan"),
            "p_perm_two_sided": 1.0,
            "n_perm": n_p,
            "top5_a": top5_a,
            "top5_b": top5_b,
            "top5_overlap_count": len(top5_ov),
            "top5_overlap_indices": top5_ov,
            "bottom5_a": bot5_a,
            "bottom5_b": bot5_b,
            "bottom5_overlap_count": len(bot5_ov),
            "bottom5_overlap_indices": bot5_ov,
        }

    za = ra_c / na
    zb = rb_c / nb
    rho_obs = float(np.clip(float(np.dot(za, zb)), -1.0, 1.0))
    # Exact vectorized permutation null: dot(za, zb[perms[k]])
    rho_perm = zb[perms] @ za
    p_perm = float((1 + int(np.sum(np.abs(rho_perm) >= abs(rho_obs)))) / float(n_p + 1))

    return {
        "spearman_rho": rho_obs,
        "p_perm_two_sided": p_perm,
        "n_perm": n_p,
        "perm_abs_q95": float(np.percentile(np.abs(rho_perm), 95.0)),
        "top5_a": top5_a,
        "top5_b": top5_b,
        "top5_overlap_count": len(top5_ov),
        "top5_overlap_indices": top5_ov,
        "bottom5_a": bot5_a,
        "bottom5_b": bot5_b,
        "bottom5_overlap_count": len(bot5_ov),
        "bottom5_overlap_indices": bot5_ov,
    }


def compute_j5b_rank_correlations_r17(
    scores_by_model_and_line: Mapping[str, Mapping[str, np.ndarray]],
    line_roles: Mapping[str, str],
    line_n_eval: Mapping[str, int],
    n_perm: int = DEFAULT_R17_N_PERM,
    seed: int = DEFAULT_R17_SEED,
) -> Dict[str, Any]:
    """
    `PREREG_J §3.5` Item 2 (`J5b`) across all 4 locked pair families:
      1. `confirmatory_pairs`: all unordered pairs among confirmatory units `C_conf` (per model).
      2. `pfizer_within_family_pair`: `(pfizer_telohaec_basal, pfizer_telohaec_cytokine)` (per model).
      3. `all_fresh_n_ge_16_pairs`: all unordered pairs among fresh lines with `N_eval >= 16`
         (excluding R16 subset sensitivity `a549_ko_cas9pos`; per model).
      4. `cross_model_within_line`: `STATE_SE` vs `STATE_HVG` on the identical line `t`.
    Also logs BH FDR within each pair set.
    """
    perm_table = _build_index_perm_table_r17(n_perm=n_perm, seed=seed)

    def _eval_pair_list(
        pair_list: Sequence[Tuple[str, str]],
        vec_lookup_a: Mapping[str, np.ndarray],
        vec_lookup_b: Mapping[str, np.ndarray],
        label_fmt: str = "{a}__vs__{b}",
    ) -> Dict[str, Any]:
        entries: Dict[str, Dict[str, Any]] = {}
        keys: List[str] = []
        pvals: List[float] = []
        for la, lb in pair_list:
            if la not in vec_lookup_a or lb not in vec_lookup_b:
                continue
            pk = label_fmt.format(a=la, b=lb)
            res = compute_spearman_perm_pair_r17(vec_lookup_a[la], vec_lookup_b[lb], perm_table=perm_table)
            res["unit_a"] = la
            res["unit_b"] = lb
            entries[pk] = res
            keys.append(pk)
            pvals.append(float(res["p_perm_two_sided"]))
        qvals = _bh_fdr(pvals)
        for pk, q in zip(keys, qvals):
            entries[pk]["bh_fdr_q"] = float(q)
        rhos = [entries[pk]["spearman_rho"] for pk in keys if np.isfinite(entries[pk]["spearman_rho"])]
        return {
            "n_pairs": len(keys),
            "median_spearman_rho": float(np.median(rhos)) if rhos else float("nan"),
            "mean_spearman_rho": float(np.mean(rhos)) if rhos else float("nan"),
            "min_spearman_rho": float(min(rhos)) if rhos else float("nan"),
            "max_spearman_rho": float(max(rhos)) if rhos else float("nan"),
            "pairs": entries,
        }

    models = sorted(scores_by_model_and_line.keys())
    per_model: Dict[str, Any] = {}
    for mkey in models:
        lmap = scores_by_model_and_line[mkey]
        avail = sorted(lmap.keys())
        conf_lines = sorted(ln for ln in avail if line_roles.get(ln) == "confirmatory")
        conf_pairs = [(conf_lines[i], conf_lines[j]) for i in range(len(conf_lines)) for j in range(i + 1, len(conf_lines))]

        pfizer_pair = (
            [("pfizer_telohaec_basal", "pfizer_telohaec_cytokine")]
            if ("pfizer_telohaec_basal" in lmap and "pfizer_telohaec_cytokine" in lmap)
            else []
        )

        fresh_ge16 = sorted(
            ln
            for ln in avail
            if line_roles.get(ln) != "r16_subset_sensitivity" and int(line_n_eval.get(ln, 0)) >= 16
        )
        fresh_pairs = [(fresh_ge16[i], fresh_ge16[j]) for i in range(len(fresh_ge16)) for j in range(i + 1, len(fresh_ge16))]

        per_model[mkey] = {
            "confirmatory_pairs": _eval_pair_list(conf_pairs, lmap, lmap),
            "pfizer_within_family_pair": _eval_pair_list(pfizer_pair, lmap, lmap),
            "all_fresh_n_ge_16_pairs": _eval_pair_list(fresh_pairs, lmap, lmap),
        }

    cross_model: Dict[str, Any] = {}
    if "STATE_SE" in scores_by_model_and_line and "STATE_HVG" in scores_by_model_and_line:
        se_map = scores_by_model_and_line["STATE_SE"]
        hvg_map = scores_by_model_and_line["STATE_HVG"]
        common_lines = sorted(set(se_map.keys()) & set(hvg_map.keys()))
        cm_pairs = [(ln, ln) for ln in common_lines]
        cross_model = _eval_pair_list(cm_pairs, se_map, hvg_map, label_fmt="{a}__STATE_SE_vs_STATE_HVG")

    return {
        "n_perm": int(n_perm),
        "seed": int(seed),
        "shared_perm_table_order_invariant": True,
        "by_model": per_model,
        "cross_model_within_line": cross_model,
    }


# ==============================================================================
# 7. Parametric Pre-Inference Power Fill-In (`PREREG_J §4.3`)
# ==============================================================================


def compute_parametric_power_r17(n_eval: int, k_conf: int = 3) -> Dict[str, Any]:
    """Closed-form `PREREG_J §4.3` pre-inference precision, MDE, and conditional power."""
    n = max(int(n_eval), 1)
    k = max(int(k_conf), 1)
    alpha_step1 = ALPHA_FAMILY / float(3 * k)
    z_line = float(norm.ppf(1.0 - ALPHA_LINE_BASE / 2.0))
    z_step1 = float(norm.ppf(1.0 - alpha_step1 / 2.0))
    z_80 = float(norm.ppf(0.80))

    se_hi = R16_SE_SCALE_HI / math.sqrt(float(n))
    se_lo = R16_SE_SCALE_LO / math.sqrt(float(n))

    def _cond_pow(delta: float, se: float, z_crit: float) -> float:
        return float(1.0 - norm.cdf(z_crit - delta / se))

    # Variance-component power at tau = 0.02 (df = 55, sigma_e = se_lo / sqrt(2))
    sigma_e = se_lo / math.sqrt(2.0)
    crit_chi2 = float(chi2.ppf(1.0 - alpha_step1, df=55))
    scale_factor = 1.0 + (0.02 ** 2) / (sigma_e ** 2)
    var_pow_lo = float(1.0 - chi2.cdf(crit_chi2 / scale_factor, df=55))

    return {
        "n_eval": int(n_eval),
        "k_conf": int(k_conf),
        "se_pair_hi": se_hi,
        "se_pair_lo": se_lo,
        "hw_9833_hi": z_line * se_hi,
        "hw_9833_lo": z_line * se_lo,
        "hw_holm_step1_hi": z_step1 * se_hi,
        "hw_holm_step1_lo": z_step1 * se_lo,
        "mde_80_holm_step1_hi": (z_step1 + z_80) * se_hi,
        "mde_80_holm_step1_lo": (z_step1 + z_80) * se_lo,
        "cond_power_delta004_holm_step1_hi": _cond_pow(0.04, se_hi, z_step1),
        "cond_power_delta004_holm_step1_lo": _cond_pow(0.04, se_lo, z_step1),
        "cond_power_delta006_holm_step1_hi": _cond_pow(0.06, se_hi, z_step1),
        "cond_power_delta006_holm_step1_lo": _cond_pow(0.06, se_lo, z_step1),
        "var_test_power_tau002_holm_step1_lo": var_pow_lo,
    }


# ==============================================================================
# 8. Full Single-Line Single-Model Evaluator & Cohort Holm Adjudicator
# ==============================================================================


def evaluate_model_on_line_r17(
    line_name: str,
    model_key: str,
    delta_obs: np.ndarray,
    delta_half_a: np.ndarray,
    delta_half_b: np.ndarray,
    half_valid: np.ndarray,
    delta_transfer: np.ndarray,
    source_grand_mean: np.ndarray,
    target_indices: np.ndarray,
    target_names: Sequence[str],
    fct_delta_pred: np.ndarray,
    delta_pred_56: np.ndarray,
    ntc_pb_ensemble_56: np.ndarray,
    obs_ntc_mean_eval: np.ndarray,
    ntc_pb_per_fold_56: Optional[np.ndarray] = None,
    fold_names: Optional[Sequence[str]] = None,
    per_fold_r_orth_cal_56: Optional[np.ndarray] = None,
    secondary_comparators: Optional[Mapping[str, np.ndarray]] = None,
    medoid_batch_idx: Optional[int] = None,
    mismatched_by_source_and_idx: Optional[Mapping[str, Mapping[int, np.ndarray]]] = None,
    seen_sources: Sequence[str] = (),
    unseen_sources: Sequence[str] = (),
    same_lineage_sources: Sequence[str] = (),
    extra_fct_sensitivity: Optional[Mapping[str, Mapping[str, np.ndarray]]] = None,
    n_boot: int = DEFAULT_R17_N_BOOT,
    n_boot_bonf56: int = DEFAULT_R17_N_BOOT_BONF56,
    seed: int = DEFAULT_R17_SEED,
    n_splits_diag: int = DEFAULT_R17_N_SPLITS_DIAG,
) -> Dict[str, Any]:
    """
    Evaluate ALL Round 17 single-line statistics for `(line_name, model_key)` from arrays alone
    (base line-level `alpha = 0.05/3`; confirmatory Holm-Bonferroni overrides are applied in
    `adjudicate_cohort_r17` without recomputing bootstraps).
    """
    comp_map: Dict[str, np.ndarray] = dict(secondary_comparators or {})
    if "multi_batch_ens8" not in comp_map:
        ens8_idx = (0, 7, 14, 21, 28, 35, 42, 49)
        comp_map["multi_batch_ens8"] = np.mean(delta_pred_56[list(ens8_idx)], axis=0)
    if medoid_batch_idx is not None and "medoid_idx" not in comp_map:
        comp_map["medoid_idx"] = delta_pred_56[int(medoid_batch_idx)]

    # 1. Label-free L2 selection (uses strictly unlabeled NTC profiles; zero perturbation labels)
    l2_sel = select_l2_ntc_recon_match_r17(
        ntc_pb_ensemble_56=ntc_pb_ensemble_56,
        obs_ntc_mean_eval=obs_ntc_mean_eval,
        ntc_pb_per_fold_56=ntc_pb_per_fold_56,
        fold_names=fold_names,
    )
    b_ntc = int(l2_sel["selected_batch_idx"])

    # 2. Primary `r_orth_cal` moments + shared 4000 x 56 joint bootstrap
    moms_orth = compute_sweep56_moments_r17(
        delta_obs=delta_obs,
        delta_pred_56=delta_pred_56,
        fct_delta_pred=fct_delta_pred,
        delta_transfer=delta_transfer,
        target_gene_idx=target_indices,
        source_grand_mean=None,
        delta_half_a=delta_half_a,
        delta_half_b=delta_half_b,
        half_valid=half_valid,
        extra_preds=comp_map,
    )
    m = int(moms_orth["m"])
    draws = np.random.default_rng(int(seed)).integers(0, m, size=(int(n_boot), m))
    rho_orth, ceil_orth, b_ceil_orth = _compute_half_ceiling_draws(moms_orth["half_mom"], draws)

    r_hat_56 = np.empty(N_BATCH_INDICES_STATE, dtype=np.float64)
    r_raw_56 = np.empty(N_BATCH_INDICES_STATE, dtype=np.float64)
    r_boot_56 = np.empty((int(n_boot), N_BATCH_INDICES_STATE), dtype=np.float64)
    per_index_orth: List[Dict[str, Any]] = []
    for b in range(N_BATCH_INDICES_STATE):
        pt_b, boot_b, raw_b = _eval_pred_mom_with_ceiling(
            moms_orth["moms_56"][b], ceil_orth, b_ceil_orth, draws, zero_null_to_zero=False
        )
        r_hat_56[b] = pt_b
        r_raw_56[b] = raw_b
        r_boot_56[:, b] = boot_b
        s_b = summarize_boot_r17(pt_b, boot_b, m, n_boot)
        s_b["batch_idx"] = int(b)
        s_b["raw_r"] = float(raw_b)
        s_b["rho_ab_orth"] = float(rho_orth)
        s_b["ceiling"] = float(ceil_orth)
        per_index_orth.append(s_b)

    fct_pt, fct_boot, fct_raw = _eval_pred_mom_with_ceiling(
        moms_orth["mom_fct"], ceil_orth, b_ceil_orth, draws, zero_null_to_zero=False
    )
    fct_sum = summarize_boot_r17(fct_pt, fct_boot, m, n_boot)
    fct_sum["raw_r"] = float(fct_raw)
    fct_sum["rho_ab_orth"] = float(rho_orth)
    fct_sum["ceiling"] = float(ceil_orth)

    l1_pt, l1_boot, l1_raw = _eval_pred_mom_with_ceiling(
        moms_orth["mom_l1"], ceil_orth, b_ceil_orth, draws, zero_null_to_zero=False
    )
    l1_sum = summarize_boot_r17(l1_pt, l1_boot, m, n_boot)
    l1_sum["raw_r"] = float(l1_raw)
    l1_sum["rho_ab_orth"] = float(rho_orth)
    l1_sum["ceiling"] = float(ceil_orth)

    # Full 56-index descriptive landscape (`§3.2` Item 1)
    fin_56 = r_hat_56[np.isfinite(r_hat_56)]
    b0_rank = _rank_1based_desc_lowest_tie(r_hat_56, 0)
    b0_pct = float(100.0 * (N_BATCH_INDICES_STATE - b0_rank) / float(N_BATCH_INDICES_STATE - 1))
    q10, q25, q50, q75, q90 = (
        [float(x) for x in np.percentile(fin_56, [10.0, 25.0, 50.0, 75.0, 90.0])]
        if fin_56.size
        else [float("nan")] * 5
    )
    landscape_orth = {
        "n_indices": N_BATCH_INDICES_STATE,
        "n_finite": int(fin_56.size),
        "min": float(np.min(fin_56)) if fin_56.size else float("nan"),
        "q10": q10,
        "q25": q25,
        "median": q50,
        "mean": float(np.mean(fin_56)) if fin_56.size else float("nan"),
        "q75": q75,
        "q90": q90,
        "max": float(np.max(fin_56)) if fin_56.size else float("nan"),
        "iqr": float(q75 - q25) if (np.isfinite(q25) and np.isfinite(q75)) else float("nan"),
        "raw_range": float(np.max(fin_56) - np.min(fin_56)) if fin_56.size else float("nan"),
        "observed_sd": float(np.std(fin_56, ddof=1)) if fin_56.size > 1 else 0.0,
        "argmin_idx": _argmin_lowest_tie(r_hat_56),
        "argmax_idx": _argmax_lowest_tie(r_hat_56),
        "b0_score": float(r_hat_56[0]),
        "b0_rank_1based_desc": b0_rank,
        "b0_percentile": b0_pct,
        "rho_ab_orth": float(rho_orth),
        "ceiling": float(ceil_orth),
        "scores_56": [float(x) for x in r_hat_56],
        "per_index": per_index_orth,
    }

    # 3. Variance-component estimator & cross-fitted spread (`H-J1`)
    var_comp = compute_variance_component_r17(r_hat_56, r_boot_56)
    cf_res = compute_cross_fitted_spread_r17(
        moments_bundle=moms_orth,
        target_names=target_names,
        half_valid=half_valid,
        l2_selected_idx=b_ntc,
        n_boot=n_boot,
        seed=seed,
        n_splits_diag=n_splits_diag,
    )
    detectable_sig = evaluate_detectable_signal_r17(cf_res["split_s0"]["r_cf_oracle_0"], l1_sum)

    h_j1_base = adjudicate_h_j1_r17(
        var_comp=var_comp,
        cf_s0=cf_res["split_s0"],
        boot_delta_cf=cf_res["_boot_delta_cf_0"],
        detectable_signal=detectable_sig,
        landscape_summary=landscape_orth,
        l1_ens56_sum=l1_sum,
        rho_ab_orth_full=rho_orth,
        alpha_adj=ALPHA_LINE_BASE,
        holm_rejected=None,
    )

    # 4. `H-J2` Call Instability + Bonferroni-56 companion
    h_j2_res = adjudicate_h_j2_r17(
        r_hat_56=r_hat_56,
        r_boot_56=r_boot_56,
        fct_pt=fct_pt,
        fct_boot=fct_boot,
        fct_summary=fct_sum,
        moments_bundle=moms_orth,
        detectable_signal=detectable_sig,
        landscape_summary=landscape_orth,
        l1_ens56_sum=l1_sum,
        n_boot=n_boot,
        n_boot_bonf56=n_boot_bonf56,
        seed=seed,
    )

    # 5. `H-J3`: L1_ens56 (`H-J3a`), L2_ntc_recon_match (`H-J3b`), & secondary comparators
    diff_l1_vs_b0_pt = float(l1_pt - r_hat_56[0])
    diff_l1_vs_b0_boot = l1_boot - r_boot_56[:, 0]
    l1_vs_fct_sum = summarize_boot_r17(float(l1_pt - fct_pt), l1_boot - fct_boot, m, n_boot)
    h_j3a_base = adjudicate_h_j3_rule_r17(
        rule_name="L1_ens56",
        diff_vs_b0_pt=diff_l1_vs_b0_pt,
        diff_vs_b0_boot=diff_l1_vs_b0_boot,
        rule_vs_fct_sum=l1_vs_fct_sum,
        rule_vs_cf_oracle_sum=cf_res["split_s0"]["l1_vs_cf_oracle"],
        alpha_adj=ALPHA_LINE_BASE,
        m_items=m,
        n_boot=n_boot,
        holm_rejected=None,
        exact_tie_selected_zero=False,
    )

    diff_l2_vs_b0_pt = float(r_hat_56[b_ntc] - r_hat_56[0])
    diff_l2_vs_b0_boot = r_boot_56[:, b_ntc] - r_boot_56[:, 0]
    l2_vs_fct_sum = summarize_boot_r17(float(r_hat_56[b_ntc] - fct_pt), r_boot_56[:, b_ntc] - fct_boot, m, n_boot)
    h_j3b_base = adjudicate_h_j3_rule_r17(
        rule_name="L2_ntc_recon_match",
        diff_vs_b0_pt=diff_l2_vs_b0_pt,
        diff_vs_b0_boot=diff_l2_vs_b0_boot,
        rule_vs_fct_sum=l2_vs_fct_sum,
        rule_vs_cf_oracle_sum=cf_res["split_s0"]["l2_vs_cf_oracle"],
        alpha_adj=ALPHA_LINE_BASE,
        m_items=m,
        n_boot=n_boot,
        holm_rejected=None,
        exact_tie_selected_zero=(b_ntc == 0),
    )
    # Also report negative-L2 secondary index's r_orth_cal score
    b_negl2 = int(l2_sel["secondary_neg_l2"]["top1_idx"])
    l2_sel["secondary_neg_l2"]["selected_r_orth_cal"] = float(r_hat_56[b_negl2])
    l2_sel["secondary_neg_l2"]["vs_b0"] = summarize_boot_r17(
        float(r_hat_56[b_negl2] - r_hat_56[0]), r_boot_56[:, b_negl2] - r_boot_56[:, 0], m, n_boot
    )

    # Secondary neutral-batch comparators (`§3.4` Item 4)
    sec_comp_orth: Dict[str, Any] = {}
    for ck in sorted(moms_orth["moms_extra"].keys()):
        pt_c, boot_c, raw_c = _eval_pred_mom_with_ceiling(
            moms_orth["moms_extra"][ck], ceil_orth, b_ceil_orth, draws, zero_null_to_zero=False
        )
        c_sum = summarize_boot_r17(pt_c, boot_c, m, n_boot)
        c_sum["raw_r"] = float(raw_c)
        c_vs_b0 = summarize_boot_r17(float(pt_c - r_hat_56[0]), boot_c - r_boot_56[:, 0], m, n_boot)
        c_vs_fct = summarize_boot_r17(float(pt_c - fct_pt), boot_c - fct_boot, m, n_boot)
        sec_comp_orth[ck] = {
            "r_orth_cal": c_sum,
            "vs_b0": c_vs_b0,
            "vs_fct": c_vs_fct,
            "vs_fct_h_i2_call": adjudicate_claim_i2_line(c_vs_fct),
        }
        if ck == "medoid_idx" and medoid_batch_idx is not None:
            sec_comp_orth[ck]["medoid_batch_idx"] = int(medoid_batch_idx)

    # 6. Secondary `r_dorth_cal` mirror (`§1.2`)
    moms_dorth = compute_sweep56_moments_r17(
        delta_obs=delta_obs,
        delta_pred_56=delta_pred_56,
        fct_delta_pred=fct_delta_pred,
        delta_transfer=delta_transfer,
        target_gene_idx=target_indices,
        source_grand_mean=source_grand_mean,
        delta_half_a=delta_half_a,
        delta_half_b=delta_half_b,
        half_valid=half_valid,
        extra_preds=comp_map,
    )
    rho_dorth, ceil_dorth, b_ceil_dorth = _compute_half_ceiling_draws(moms_dorth["half_mom"], draws)
    r_dorth_56 = np.empty(N_BATCH_INDICES_STATE, dtype=np.float64)
    r_dorth_boot_56 = np.empty((int(n_boot), N_BATCH_INDICES_STATE), dtype=np.float64)
    for b in range(N_BATCH_INDICES_STATE):
        pt_db, boot_db, _ = _eval_pred_mom_with_ceiling(
            moms_dorth["moms_56"][b], ceil_dorth, b_ceil_dorth, draws, zero_null_to_zero=True
        )
        r_dorth_56[b] = pt_db
        r_dorth_boot_56[:, b] = boot_db
    fct_d_pt, fct_d_boot, _ = _eval_pred_mom_with_ceiling(
        moms_dorth["mom_fct"], ceil_dorth, b_ceil_dorth, draws, zero_null_to_zero=True
    )
    l1_d_pt, l1_d_boot, _ = _eval_pred_mom_with_ceiling(
        moms_dorth["mom_l1"], ceil_dorth, b_ceil_dorth, draws, zero_null_to_zero=True
    )
    var_comp_dorth = compute_variance_component_r17(r_dorth_56, r_dorth_boot_56)
    cf_dorth = compute_cross_fitted_spread_r17(
        moments_bundle=moms_dorth,
        target_names=target_names,
        half_valid=half_valid,
        l2_selected_idx=b_ntc,
        n_boot=n_boot,
        seed=seed,
        n_splits_diag=n_splits_diag,
    )
    dorth_fin = r_dorth_56[np.isfinite(r_dorth_56)]
    dorth_summary = {
        "rho_ab_orth": float(rho_dorth),
        "ceiling": float(ceil_dorth),
        "scores_56": [float(x) for x in r_dorth_56],
        "min": float(np.min(dorth_fin)) if dorth_fin.size else float("nan"),
        "median": float(np.median(dorth_fin)) if dorth_fin.size else float("nan"),
        "max": float(np.max(dorth_fin)) if dorth_fin.size else float("nan"),
        "raw_range": float(np.max(dorth_fin) - np.min(dorth_fin)) if dorth_fin.size else float("nan"),
        "b0_score": float(r_dorth_56[0]),
        "fct": summarize_boot_r17(fct_d_pt, fct_d_boot, m, n_boot),
        "l1_ens56": summarize_boot_r17(l1_d_pt, l1_d_boot, m, n_boot),
        "l2_ntc_recon_match": summarize_boot_r17(r_dorth_56[b_ntc], r_dorth_boot_56[:, b_ntc], m, n_boot),
        "variance_component": var_comp_dorth,
        "cross_fitted_s0": cf_dorth["split_s0"],
    }

    # 7. `J5a` Cross-Fitted `G4` 2x2 Decomposition (when mismatched arrays provided)
    j5a_res: Optional[Dict[str, Any]] = None
    if mismatched_by_source_and_idx:
        j5a_res = compute_j5a_g4_r17(
            delta_obs=delta_obs,
            delta_half_a=delta_half_a,
            delta_half_b=delta_half_b,
            half_valid=half_valid,
            delta_transfer=delta_transfer,
            source_grand_mean=source_grand_mean,
            target_indices=target_indices,
            delta_pred_56_matched=delta_pred_56,
            mismatched_by_source_and_idx=mismatched_by_source_and_idx,
            seen_sources=seen_sources,
            unseen_sources=unseen_sources,
            same_lineage_sources=same_lineage_sources,
            p1_indices=cf_res["split_s0"]["p1_indices"],
            p2_indices=cf_res["split_s0"]["p2_indices"],
            b_star_1=int(cf_res["split_s0"]["b_star_1"]),
            b_dag_1=int(cf_res["split_s0"]["b_dag_1"]),
            b_star_2=int(cf_res["split_s0"]["b_star_2"]),
            b_dag_2=int(cf_res["split_s0"]["b_dag_2"]),
            n_boot=n_boot,
            seed=seed,
        )

    # 8. Optional `gse295214_h1` full `r9_6` FCT sensitivity (fixed on primary `(T_eval, G_eval)`)
    fct_sens_out: Dict[str, Any] = {}
    if extra_fct_sensitivity:
        for skey, s_pay in sorted(extra_fct_sensitivity.items()):
            s_fct = np.asarray(s_pay["fct_delta_pred"], dtype=np.float64)
            s_tm = np.asarray(s_pay["delta_transfer"], dtype=np.float64)
            moms_s = compute_sweep56_moments_r17(
                delta_obs=delta_obs,
                delta_pred_56=delta_pred_56,
                fct_delta_pred=s_fct,
                delta_transfer=s_tm,
                target_gene_idx=target_indices,
                source_grand_mean=None,
                delta_half_a=delta_half_a,
                delta_half_b=delta_half_b,
                half_valid=half_valid,
            )
            _, c_s, bc_s = _compute_half_ceiling_draws(moms_s["half_mom"], draws)
            r56_s = np.empty(N_BATCH_INDICES_STATE, dtype=np.float64)
            rb56_s = np.empty((int(n_boot), N_BATCH_INDICES_STATE), dtype=np.float64)
            for b in range(N_BATCH_INDICES_STATE):
                r56_s[b], rb56_s[:, b], _ = _eval_pred_mom_with_ceiling(
                    moms_s["moms_56"][b], c_s, bc_s, draws, zero_null_to_zero=False
                )
            fct_s_pt, fct_s_boot, _ = _eval_pred_mom_with_ceiling(
                moms_s["mom_fct"], c_s, bc_s, draws, zero_null_to_zero=False
            )
            l1_s_pt, l1_s_boot, _ = _eval_pred_mom_with_ceiling(
                moms_s["mom_l1"], c_s, bc_s, draws, zero_null_to_zero=False
            )
            fct_s_sum = summarize_boot_r17(fct_s_pt, fct_s_boot, m, n_boot)
            l1_s_sum = summarize_boot_r17(l1_s_pt, l1_s_boot, m, n_boot)
            fin_s = r56_s[np.isfinite(r56_s)]
            land_s = {"max": float(np.max(fin_s)) if fin_s.size else float("nan")}
            h_j2_s = adjudicate_h_j2_r17(
                r_hat_56=r56_s,
                r_boot_56=rb56_s,
                fct_pt=fct_s_pt,
                fct_boot=fct_s_boot,
                fct_summary=fct_s_sum,
                moments_bundle=moms_s,
                detectable_signal=detectable_sig,
                landscape_summary=land_s,
                l1_ens56_sum=l1_s_sum,
                n_boot=n_boot,
                n_boot_bonf56=n_boot_bonf56,
                seed=seed,
            )
            fct_sens_out[skey] = {
                "fct_r_orth_cal": fct_s_sum,
                "b0_vs_fct": compute_paired_contrast_r16(
                    delta_obs, delta_pred_56[0], s_fct, s_tm, target_indices,
                    source_grand_mean=None, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
                    half_valid=half_valid, n_boot=n_boot, seed=seed,
                ),
                "l1_ens56_vs_fct": summarize_boot_r17(float(l1_s_pt - fct_s_pt), l1_s_boot - fct_s_boot, m, n_boot),
                "h_j2_sensitivity": h_j2_s,
            }

    # 9. Per-fold 56-index secondary summary (`§1.2`)
    per_fold_diag: Optional[Dict[str, Any]] = None
    if per_fold_r_orth_cal_56 is not None:
        pf_arr = np.asarray(per_fold_r_orth_cal_56, dtype=np.float64)
        fnames = list(fold_names) if fold_names is not None else [f"fold_{i}" for i in range(pf_arr.shape[0])]
        per_fold_diag = {}
        for fi, fn in enumerate(fnames):
            fv = pf_arr[fi]
            ff = fv[np.isfinite(fv)]
            per_fold_diag[str(fn)] = {
                "min": float(np.min(ff)) if ff.size else float("nan"),
                "median": float(np.median(ff)) if ff.size else float("nan"),
                "max": float(np.max(ff)) if ff.size else float("nan"),
                "raw_range": float(np.max(ff) - np.min(ff)) if ff.size else float("nan"),
                "b0_score": float(fv[0]),
                "argmax_idx": _argmax_lowest_tie(fv),
                "argmin_idx": _argmin_lowest_tie(fv),
                "scores_56": [float(x) for x in fv],
            }

    return {
        "line_name": line_name,
        "model_key": model_key,
        "n_eval_targets": m,
        "n_eval_genes": int(delta_obs.shape[1]),
        "landscape_r_orth_cal": landscape_orth,
        "fct_r_orth_cal": fct_sum,
        "l1_ens56_r_orth_cal": l1_sum,
        "l2_ntc_selection": l2_sel,
        "variance_component": var_comp,
        "cross_fitting": {
            "split_s0": cf_res["split_s0"],
            "multi_split_diag": cf_res["multi_split_diag"],
        },
        "detectable_signal_j2": detectable_sig,
        "H_J1_line_base": h_j1_base,
        "H_J2": h_j2_res,
        "H_J3_line_base": {
            "H_J3a_L1_ens56": h_j3a_base,
            "H_J3b_L2_ntc_recon_match": h_j3b_base,
        },
        "secondary_neutral_comparators": sec_comp_orth,
        "secondary_r_dorth_cal": dorth_summary,
        "secondary_per_fold_56": per_fold_diag,
        "secondary_j5a_g4": j5a_res,
        "secondary_extra_fct_sensitivity": fct_sens_out,
        "practical_accuracy": {
            "b0": compute_practical_accuracy_metrics_r16(delta_obs, delta_pred_56[0], target_indices),
            "L1_ens56": compute_practical_accuracy_metrics_r16(delta_obs, np.mean(delta_pred_56, axis=0), target_indices),
            "L2_ntc_recon_match": compute_practical_accuracy_metrics_r16(delta_obs, delta_pred_56[b_ntc], target_indices),
            "FCT": compute_practical_accuracy_metrics_r16(delta_obs, fct_delta_pred, target_indices),
        },
        # Internal bootstrap vectors for exact Holm-Bonferroni CI recomputation at cohort level
        "_internal": {
            "boot_delta_cf_0": cf_res["_boot_delta_cf_0"],
            "diff_l1_vs_b0_pt": diff_l1_vs_b0_pt,
            "diff_l1_vs_b0_boot": diff_l1_vs_b0_boot,
            "diff_l2_vs_b0_pt": diff_l2_vs_b0_pt,
            "diff_l2_vs_b0_boot": diff_l2_vs_b0_boot,
            "l1_vs_fct_sum": l1_vs_fct_sum,
            "l2_vs_fct_sum": l2_vs_fct_sum,
            "b_ntc": b_ntc,
            "rho_ab_orth_full": float(rho_orth),
            "m": m,
            "n_boot": int(n_boot),
        },
    }


def _strip_internal_fields(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {k: _strip_internal_fields(v) for k, v in obj.items() if not str(k).startswith("_")}
    if isinstance(obj, list):
        return [_strip_internal_fields(x) for x in obj]
    return obj


def adjudicate_cohort_r17(
    per_line_results: Mapping[str, Mapping[str, Dict[str, Any]]],
    line_roles: Mapping[str, str],
    line_n_eval: Mapping[str, int],
    n_perm: int = DEFAULT_R17_N_PERM,
    seed: int = DEFAULT_R17_SEED,
) -> Dict[str, Any]:
    """
    Apply `PREREG_J §2.1, §3.1, §3.5`:
      1. Separate Holm-Bonferroni FWER <= 0.05 across `M = 3 * K_conf` confirmatory tests for
         `STATE_SE` (primary family) and `STATE_HVG` (co-evaluated second model family).
      2. Mandatory pooled-model Holm-Bonferroni sensitivity across all `6 * K_conf` confirmatory
         tests (`STATE_SE + STATE_HVG` jointly).
      3. Locked cohort-level `k / K_conf` claim rule per model (zero pooling; general cross-context
         claim permitted iff `K_conf >= 2` and `k >= 2` independent confirmatory units reach
         `BATCH_DEPENDENCE_PRACTICALLY_MEANINGFUL`).
      4. `J5b` cross-line & cross-model 56-index Spearman rank correlations + `10,000`-perm null.
    """
    conf_lines = sorted(ln for ln, role in line_roles.items() if role == "confirmatory" and ln in per_line_results)
    k_conf = len(conf_lines)
    models_present = sorted({mkey for lres in per_line_results.values() for mkey in lres.keys()})

    def _collect_family_pvals(models_subset: Sequence[str]) -> Dict[str, float]:
        pvals: Dict[str, float] = {}
        for mkey in models_subset:
            for ln in conf_lines:
                if mkey not in per_line_results[ln]:
                    continue
                res = per_line_results[ln][mkey]
                # H-J1 IUT p-value (if invalid, force 1.0)
                hj1_ok = bool(res["H_J1_line_base"]["validity_ok"])
                p_hj1 = float(res["H_J1_line_base"]["p_h_j1_iut"]) if hj1_ok else 1.0
                pvals[f"{mkey}::{ln}::H-J1"] = p_hj1

                # H-J3a L1_ens56 vs b=0
                l1_s = res["H_J3_line_base"]["H_J3a_L1_ens56"]["vs_b0"]
                p_l1 = (
                    float(l1_s["p_two_sided"])
                    if float(l1_s["valid_frac"]) >= MIN_BOOT_VALID_FRAC and np.isfinite(float(l1_s["p_two_sided"]))
                    else 1.0
                )
                pvals[f"{mkey}::{ln}::H-J3a"] = p_l1

                # H-J3b L2_ntc_recon_match vs b=0 (exact tie -> p = 1.0)
                l2_entry = res["H_J3_line_base"]["H_J3b_L2_ntc_recon_match"]
                l2_s = l2_entry["vs_b0"]
                if l2_entry["verdict"] == "EXACT_TIE_SELECTED_ZERO":
                    p_l2 = 1.0
                elif float(l2_s["valid_frac"]) >= MIN_BOOT_VALID_FRAC and np.isfinite(float(l2_s["p_two_sided"])):
                    p_l2 = float(l2_s["p_two_sided"])
                else:
                    p_l2 = 1.0
                pvals[f"{mkey}::{ln}::H-J3b"] = p_l2
        return pvals

    def _apply_holm_to_confirmatory(
        holm_table: Mapping[str, Dict[str, Any]],
        models_subset: Sequence[str],
    ) -> Dict[str, Dict[str, Any]]:
        adj_by_model: Dict[str, Dict[str, Any]] = {}
        for mkey in models_subset:
            m_lines: Dict[str, Any] = {}
            k_meaningful = 0
            for ln in conf_lines:
                if mkey not in per_line_results[ln]:
                    continue
                res = per_line_results[ln][mkey]
                intr = res["_internal"]
                h1 = holm_table[f"{mkey}::{ln}::H-J1"]
                h3a = holm_table[f"{mkey}::{ln}::H-J3a"]
                h3b = holm_table[f"{mkey}::{ln}::H-J3b"]

                hj1_holm = adjudicate_h_j1_r17(
                    var_comp=res["variance_component"],
                    cf_s0=res["cross_fitting"]["split_s0"],
                    boot_delta_cf=intr["boot_delta_cf_0"],
                    detectable_signal=res["detectable_signal_j2"],
                    landscape_summary=res["landscape_r_orth_cal"],
                    l1_ens56_sum=res["l1_ens56_r_orth_cal"],
                    rho_ab_orth_full=intr["rho_ab_orth_full"],
                    alpha_adj=float(h1["effective_alpha"]),
                    holm_rejected=bool(h1["rejected"]),
                )
                hj1_holm["holm_step"] = h1

                hj3a_holm = adjudicate_h_j3_rule_r17(
                    rule_name="L1_ens56",
                    diff_vs_b0_pt=intr["diff_l1_vs_b0_pt"],
                    diff_vs_b0_boot=intr["diff_l1_vs_b0_boot"],
                    rule_vs_fct_sum=intr["l1_vs_fct_sum"],
                    rule_vs_cf_oracle_sum=res["cross_fitting"]["split_s0"]["l1_vs_cf_oracle"],
                    alpha_adj=float(h3a["effective_alpha"]),
                    m_items=intr["m"],
                    n_boot=intr["n_boot"],
                    holm_rejected=bool(h3a["rejected"]),
                    exact_tie_selected_zero=False,
                )
                hj3a_holm["holm_step"] = h3a

                hj3b_holm = adjudicate_h_j3_rule_r17(
                    rule_name="L2_ntc_recon_match",
                    diff_vs_b0_pt=intr["diff_l2_vs_b0_pt"],
                    diff_vs_b0_boot=intr["diff_l2_vs_b0_boot"],
                    rule_vs_fct_sum=intr["l2_vs_fct_sum"],
                    rule_vs_cf_oracle_sum=res["cross_fitting"]["split_s0"]["l2_vs_cf_oracle"],
                    alpha_adj=float(h3b["effective_alpha"]),
                    m_items=intr["m"],
                    n_boot=intr["n_boot"],
                    holm_rejected=bool(h3b["rejected"]),
                    exact_tie_selected_zero=(int(intr["b_ntc"]) == 0),
                )
                hj3b_holm["holm_step"] = h3b

                if hj1_holm["verdict"] == "BATCH_DEPENDENCE_PRACTICALLY_MEANINGFUL":
                    k_meaningful += 1
                m_lines[ln] = {
                    "H_J1_holm": hj1_holm,
                    "H_J2": res["H_J2"]["classification"],
                    "H_J3a_L1_ens56_holm": hj3a_holm,
                    "H_J3b_L2_ntc_recon_match_holm": hj3b_holm,
                }

            adj_by_model[mkey] = {
                "K_conf": int(k_conf),
                "confirmatory_units": list(conf_lines),
                "k_batch_dependence_practically_meaningful": int(k_meaningful),
                "k_over_K": f"{k_meaningful}/{k_conf}",
                "general_cross_context_claim_permitted": bool(k_conf >= 2 and k_meaningful >= 2),
                "claim_rule_note": (
                    "GENERAL_CROSS_CONTEXT_CLAIM_PERMITTED (k >= 2 independent confirmatory units)"
                    if (k_conf >= 2 and k_meaningful >= 2)
                    else "PER_LINE_STATEMENTS_ONLY (k < 2 or K_conf < 2; no general cross-context claim permitted)"
                ),
                "per_confirmatory_line": m_lines,
            }
        return adj_by_model

    # Primary separate Holm families (`M = 3 * K_conf` per model)
    separate_holm_tables: Dict[str, Any] = {}
    separate_cohort_calls: Dict[str, Any] = {}
    for mkey in models_present:
        pvals_m = _collect_family_pvals([mkey])
        ht_m = apply_holm_bonferroni_r17(pvals_m, alpha_family=ALPHA_FAMILY)
        separate_holm_tables[mkey] = ht_m
        separate_cohort_calls[mkey] = _apply_holm_to_confirmatory(ht_m, [mkey])[mkey]

    # Mandatory pooled-model Holm sensitivity (`6 * K_conf` across STATE_SE + STATE_HVG jointly)
    pooled_pvals = _collect_family_pvals(models_present)
    pooled_holm_table = apply_holm_bonferroni_r17(pooled_pvals, alpha_family=ALPHA_FAMILY)
    pooled_cohort_calls = _apply_holm_to_confirmatory(pooled_holm_table, models_present)

    # J5b Cross-line & cross-model 56-index Spearman rank correlations
    scores_by_model_and_line: Dict[str, Dict[str, np.ndarray]] = {mkey: {} for mkey in models_present}
    for ln, lres in sorted(per_line_results.items()):
        for mkey, mres in lres.items():
            scores_by_model_and_line[mkey][ln] = np.asarray(
                mres["landscape_r_orth_cal"]["scores_56"], dtype=np.float64
            )
    j5b_res = compute_j5b_rank_correlations_r17(
        scores_by_model_and_line=scores_by_model_and_line,
        line_roles=line_roles,
        line_n_eval=line_n_eval,
        n_perm=n_perm,
        seed=seed,
    )

    # Parametric pre-inference power fill-in (`§4.3`) for every line
    power_fill_in = {
        ln: compute_parametric_power_r17(int(line_n_eval[ln]), k_conf=max(k_conf, 1))
        for ln in sorted(line_n_eval.keys())
    }

    return {
        "K_conf": int(k_conf),
        "confirmatory_units": list(conf_lines),
        "separate_model_holm_fwer05": {
            "holm_tables": separate_holm_tables,
            "cohort_adjudications": separate_cohort_calls,
        },
        "pooled_model_holm_sensitivity_fwer05": {
            "m_total_tests": len(pooled_pvals),
            "holm_table": pooled_holm_table,
            "cohort_adjudications": pooled_cohort_calls,
        },
        "j5b_cross_line_rank_correlations": j5b_res,
        "parametric_power_fill_in_s43": power_fill_in,
        "lines_clean": {ln: _strip_internal_fields(lres) for ln, lres in sorted(per_line_results.items())},
    }
