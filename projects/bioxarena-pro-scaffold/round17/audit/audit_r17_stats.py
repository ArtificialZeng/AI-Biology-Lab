#!/usr/bin/env python3
"""
Clean-Room Round 17 Statistical Engine, Recomputation-From-Arrays, & Mechanical Adjudicator
(`A-T3`, strictly following `PREREG_J_DRAFT_v2.md` & `PREREG_J_CLARIFICATIONS_C1.md`).

Implements from exported `.npz` arrays alone (`allow_pickle=False`):
  1. Vectorized `4,000`-draw joint perturbation bootstrap across all 56 batch indices (`0..55`),
     exact `r_orth_cal` (`all_nontarget`) + secondary `r_dorth_cal`.
  2. `H-J1`: Winner's-curse-free half-split cross-fit `Delta_CF,0` / `r_CF-Oracle,0` (`C1.4`),
     multi-split diagnostic `s = 0..19`, noise-debiased `tau_batch` (`ddof=1` `V_noise`),
     centered-bootstrap `p_var` (`C1.6`), intersection-union `p_H-J1 = max(p_var, p_Delta)` (`C1.5`),
     floor-effect predicate `has_detectable_signal` (`J2`), and 6-way classification (`C1.7`).
  3. `H-J2`: `56`-index `STATE - FCT` call instability (`f_flip`) + companion Bonferroni-56
     (`z = Phi^-1(1 - 0.025/112)` normal AND `20,000`-draw empirical quantile CIs per `C1.8`).
  4. `H-J3`: Label-free `L1_ens56` & `L2_ntc_recon_match` (`Pearson` primary + `-L2` + per-fold
     + mandatory `top1 - top2` margin per Ruling 3) vs `b=0` and `CF-Oracle`, plus neutral-batch
     comparators (`bypass_zero`, `mean_batch_emb`, `multi_batch_ens8`, Euclidean/cosine `medoid`).
  5. `J5a`: Cross-fitted `(b*_1, b_dag_1)` on `P2` + auxiliary `(b*_2, b_dag_2)` on `P1` (`C1.11`)
     + pre-fixed continuity `(0, 17)` & `(0, 35)` `G4` 2x2 decompositions.
  6. `J5b`: Cross-line 56-index Spearman rank correlations using ONE shared `10,000 x 56`
     permutation table (`seed=20261701`, `C1.10`) + BH FDR.
  7. Holm-Bonferroni FWER (`M = 3 * K_conf`) separately per model (`STATE_SE`, `STATE_HVG`)
     and pooled-model sensitivity (`6 * K_conf`).
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.stats import norm, rankdata

AUDIT_R17_DIR = Path(__file__).resolve().parent
R16_AUDIT_DIR = AUDIT_R17_DIR.parents[1] / "round16" / "audit"
for _p in (str(AUDIT_R17_DIR), str(R16_AUDIT_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from audit_r16_metrics import (  # noqa: E402
    EPS,
    NULL_NORM_TOL,
    _extract_pert_gene_indices,
    _moments,
    _moments_to_corr_and_slope,
    _project_streams,
    compute_orth_suite_r16,
    compute_paired_contrast_r16,
    disattenuate,
    spearman_brown_ceiling,
)
from audit_r17_controls import (  # noqa: E402
    DEFAULT_R17_SEED,
    R17_ALL_EVAL_LINES,
    R17_CONFIRMATORY_LINES,
    stratified_split_halves_c1,
)

MIN_BOOT_VALID_FRAC = 0.95
TOST_DELTA_R17 = 0.04
TAU_BATCH_THRESH = 0.02  # 2 * tau_batch > 0.04


def _percentile_ci(arr: np.ndarray, alpha: float) -> Tuple[float, float, float]:
    """Return `(lo, hi, boot_se)` at two-sided level `alpha` over finite bootstrap draws."""
    fin = arr[np.isfinite(arr)]
    if fin.size == 0:
        return float("nan"), float("nan"), float("nan")
    q_lo = 100.0 * (float(alpha) / 2.0)
    q_hi = 100.0 * (1.0 - float(alpha) / 2.0)
    lo, hi = np.percentile(fin, [q_lo, q_hi])
    se = float(np.std(fin, ddof=1)) if fin.size > 1 else 0.0
    return float(lo), float(hi), se


def _bootstrap_two_sided_p(draws_fin: np.ndarray) -> float:
    """Two-sided empirical bootstrap p-value (`C1.5`): `min(1, 2*(1 + min(#{d<=0}, #{d>=0}))/(B+1))`."""
    b_val = int(draws_fin.size)
    if b_val == 0:
        return 1.0
    n_le = int(np.sum(draws_fin <= 0.0))
    n_ge = int(np.sum(draws_fin >= 0.0))
    return float(min(1.0, 2.0 * (1.0 + min(n_le, n_ge)) / float(b_val + 1)))


def _compute_orth_vectors_and_boot(
    delta_obs: np.ndarray,
    delta_pred: np.ndarray,
    delta_transfer: np.ndarray,
    target_indices: Sequence[int],
    delta_half_a: np.ndarray,
    delta_half_b: np.ndarray,
    half_valid: np.ndarray,
    draws: Optional[np.ndarray],
    source_grand_mean: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """
    Compute point `r_orth_cal` (or `r_dorth_cal` if `source_grand_mean is not None`) AND exact
    bootstrap vector `b_cal` reusing caller-supplied `draws` `(n_boot, m)`.
    """
    if len(target_indices) != delta_obs.shape[0]:
        raise ValueError(
            f"target_indices length ({len(target_indices)}) != delta_obs rows ({delta_obs.shape[0]})"
        )
    kept = _extract_pert_gene_indices(
        delta_obs=delta_obs,
        t_mean=delta_transfer,
        pred=delta_pred,
        target_cols=target_indices,
        mode="all_nontarget",
        source_grand_mean=source_grand_mean,
    )
    m = len(kept)
    n_boot = 0 if draws is None else int(draws.shape[0])
    if m == 0:
        return {
            "estimate": float("nan"),
            "boot": np.full(n_boot, float("nan"), dtype=np.float64),
            "rho_ab_orth": float("nan"),
            "valid_frac": 0.0,
            "n_eval_perts": 0,
        }

    streams = _project_streams(
        kept=kept,
        delta_obs=delta_obs,
        pred=delta_pred,
        t_mean=delta_transfer,
        source_grand_mean=source_grand_mean,
        half_a_obs=delta_half_a,
        half_b_obs=delta_half_b,
        half_valid=half_valid,
    )
    is_dorth = source_grand_mean is not None
    m_proj = _moments(streams["xo_proj"], streams["yo_proj"])
    pt_r_proj = float(_moments_to_corr_and_slope(*m_proj, zero_degenerate_y_as_zero=is_dorth)[0])

    if streams["any_hv"]:
        mh_proj = _moments(streams["ha_proj"], streams["hb_proj"])
        pt_rho_orth = float(_moments_to_corr_and_slope(*mh_proj)[0])
        pt_ceil_orth = spearman_brown_ceiling(pt_rho_orth)
    else:
        pt_rho_orth = float("nan")
        pt_ceil_orth = float("nan")

    def _safe_dis(rv: float, cv: float) -> float:
        if is_dorth and np.isfinite(rv) and abs(rv) <= NULL_NORM_TOL and np.isfinite(cv) and cv > 1e-8:
            return 0.0
        return disattenuate(rv, cv)

    pt_est = _safe_dis(pt_r_proj, pt_ceil_orth)

    if draws is None or draws.size == 0:
        b_cal = np.zeros(0, dtype=np.float64)
        valid_frac = 1.0 if np.isfinite(pt_est) else 0.0
    else:
        b_r_proj = _moments_to_corr_and_slope(
            *m_proj, draws=draws, zero_degenerate_y_as_zero=is_dorth
        )[0]
        if streams["any_hv"]:
            b_rho = _moments_to_corr_and_slope(*mh_proj, draws=draws)[0]
            b_ceil = np.asarray([spearman_brown_ceiling(float(v)) for v in b_rho], dtype=np.float64)
        else:
            b_ceil = np.full(n_boot, float("nan"), dtype=np.float64)
        b_cal = np.asarray(
            [_safe_dis(float(r), float(c)) for r, c in zip(b_r_proj, b_ceil)],
            dtype=np.float64,
        )
        valid_frac = float(np.mean(np.isfinite(b_cal)))

    return {
        "estimate": float(pt_est),
        "boot": b_cal,
        "rho_ab_orth": float(pt_rho_orth),
        "valid_frac": float(valid_frac),
        "n_eval_perts": int(m),
    }


def _compute_56_sweep_and_fct_boot(
    delta_obs: np.ndarray,
    delta_pred_56: np.ndarray,
    fct_delta_pred: np.ndarray,
    delta_transfer: np.ndarray,
    target_indices: Sequence[int],
    delta_half_a: np.ndarray,
    delta_half_b: np.ndarray,
    half_valid: np.ndarray,
    n_boot: int,
    seed: int,
    source_grand_mean: Optional[np.ndarray] = None,
    draws_override: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """
    Fast vectorized evaluation across all 56 indices AND `FCT` using ONE shared `(n_boot, m)`
    perturbation bootstrap matrix (identical draws and noise ceiling across indices and `FCT`).
    """
    p_cnt = delta_obs.shape[0]
    if draws_override is not None:
        draws = draws_override
    else:
        rng = np.random.default_rng(int(seed))
        draws = rng.integers(0, p_cnt, size=(int(n_boot), p_cnt))

    # Compute FCT first (provides shared kept mask/noise ceiling baseline check cleanly)
    fct_res = _compute_orth_vectors_and_boot(
        delta_obs=delta_obs,
        delta_pred=fct_delta_pred,
        delta_transfer=delta_transfer,
        target_indices=target_indices,
        delta_half_a=delta_half_a,
        delta_half_b=delta_half_b,
        half_valid=half_valid,
        draws=draws,
        source_grand_mean=source_grand_mean,
    )

    r_obs = np.empty(56, dtype=np.float64)
    r_boot = np.empty((int(draws.shape[0]), 56), dtype=np.float64)
    valid_fracs = np.empty(56, dtype=np.float64)
    rho_orth = float(fct_res["rho_ab_orth"])

    # Pre-check per-index
    for b in range(56):
        res_b = _compute_orth_vectors_and_boot(
            delta_obs=delta_obs,
            delta_pred=delta_pred_56[b],
            delta_transfer=delta_transfer,
            target_indices=target_indices,
            delta_half_a=delta_half_a,
            delta_half_b=delta_half_b,
            half_valid=half_valid,
            draws=draws,
            source_grand_mean=source_grand_mean,
        )
        r_obs[b] = float(res_b["estimate"])
        r_boot[:, b] = res_b["boot"]
        valid_fracs[b] = float(res_b["valid_frac"])
        if np.isfinite(res_b["rho_ab_orth"]):
            rho_orth = float(res_b["rho_ab_orth"])

    return {
        "r_obs_56": r_obs,
        "r_boot_56": r_boot,
        "valid_fracs_56": valid_fracs,
        "fct_est": float(fct_res["estimate"]),
        "fct_boot": fct_res["boot"],
        "fct_valid_frac": float(fct_res["valid_frac"]),
        "rho_ab_orth": rho_orth,
        "draws": draws,
    }


def compute_h_j1_for_subset(
    delta_obs: np.ndarray,
    delta_pred_56: np.ndarray,
    delta_pred_L1: np.ndarray,
    delta_transfer: np.ndarray,
    target_indices: Sequence[int],
    delta_half_a: np.ndarray,
    delta_half_b: np.ndarray,
    half_valid: np.ndarray,
    n_boot: int = 4000,
    seed: int = DEFAULT_R17_SEED,
    n_splits_diag: int = 20,
    source_grand_mean: Optional[np.ndarray] = None,
    alpha_adj: float = 0.05 / 3.0,
) -> Dict[str, Any]:
    """
    Compute full `H-J1` (`PREREG_J §3.2`, `C1.4`–`C1.7`):
      - `56`-index landscape statistics (`min`, quantiles, `median`, `mean`, `max`, `IQR`, `s_obs`, `R_obs`, `rank_b0`)
      - Variance-component estimators: `V_noise` (`ddof=1`), `tau_batch`, centered-bootstrap `p_var`
      - Cross-fitted spread `Delta_CF,0` & `r_CF-Oracle,0` (`s = 0`, exact within-half paired bootstrap)
      - Multi-split stability (`s = 0..n_splits_diag-1`)
      - Floor-effect `has_detectable_signal` (`J2`)
      - Mechanical 6-way call under `alpha_adj`.
    """
    p_cnt = delta_obs.shape[0]
    rng_all = np.random.default_rng(int(seed))
    draws_all = rng_all.integers(0, p_cnt, size=(int(n_boot), p_cnt))

    r_obs = np.empty(56, dtype=np.float64)
    r_boot = np.empty((int(n_boot), 56), dtype=np.float64)
    vfracs = np.empty(56, dtype=np.float64)
    rho_orth = float("nan")
    for b in range(56):
        rb = _compute_orth_vectors_and_boot(
            delta_obs=delta_obs,
            delta_pred=delta_pred_56[b],
            delta_transfer=delta_transfer,
            target_indices=target_indices,
            delta_half_a=delta_half_a,
            delta_half_b=delta_half_b,
            half_valid=half_valid,
            draws=draws_all,
            source_grand_mean=source_grand_mean,
        )
        r_obs[b] = float(rb["estimate"])
        r_boot[:, b] = rb["boot"]
        vfracs[b] = float(rb["valid_frac"])
        if np.isfinite(rb["rho_ab_orth"]):
            rho_orth = float(rb["rho_ab_orth"])

    l1_res = _compute_orth_vectors_and_boot(
        delta_obs=delta_obs,
        delta_pred=delta_pred_L1,
        delta_transfer=delta_transfer,
        target_indices=target_indices,
        delta_half_a=delta_half_a,
        delta_half_b=delta_half_b,
        half_valid=half_valid,
        draws=draws_all,
        source_grand_mean=source_grand_mean,
    )

    # Landscape summaries
    fin56 = np.all(np.isfinite(r_obs))
    if fin56:
        q10, q25, q50, q75, q90 = [float(v) for v in np.percentile(r_obs, [10.0, 25.0, 50.0, 75.0, 90.0])]
        r_min = float(np.min(r_obs))
        r_max = float(np.max(r_obs))
        r_mean = float(np.mean(r_obs))
        s_obs = float(np.std(r_obs, ddof=1))
        s2_obs = float(np.sum((r_obs - r_mean) ** 2) / 55.0)
        order_desc = np.argsort(-r_obs, kind="mergesort")
        rank_b0 = int(np.flatnonzero(order_desc == 0)[0]) + 1
        pct_b0 = float(np.mean(r_obs <= r_obs[0]))
    else:
        q10 = q25 = q50 = q75 = q90 = r_min = r_max = r_mean = s_obs = s2_obs = pct_b0 = float("nan")
        rank_b0 = -1

    # C1.6 V_noise & p_var
    joint_mask = np.all(np.isfinite(r_boot), axis=1)
    joint_valid_frac = float(np.mean(joint_mask)) if r_boot.shape[0] > 0 else 0.0
    per_index_valid_ok = bool(np.all(np.isfinite(vfracs) & (vfracs >= MIN_BOOT_VALID_FRAC)))
    r_fin = r_boot[joint_mask]
    if r_fin.shape[0] >= 2 and fin56:
        r_row_centered = r_fin - r_fin.mean(axis=1, keepdims=True)
        v_noise = float(np.sum(np.var(r_row_centered, axis=0, ddof=1)) / 55.0)
        tau_batch = float(math.sqrt(max(0.0, s2_obs - v_noise)))
        z_null = r_fin - r_fin.mean(axis=0, keepdims=True)
        z_row_c = z_null - z_null.mean(axis=1, keepdims=True)
        s2_null = np.sum(z_row_c ** 2, axis=1) / 55.0
        p_var = float(1 + int(np.sum(s2_null >= s2_obs))) / float(r_fin.shape[0] + 1)
    else:
        v_noise = tau_batch = float("nan")
        p_var = 1.0

    # Cross-fitted s=0 (C1.4 & PREREG_J §3.2.2)
    def _eval_split(s_idx: int, with_boot: bool) -> Dict[str, Any]:
        p1_idx, p2_idx = stratified_split_halves_c1(half_valid, split_idx=s_idx, base_seed=seed)
        if p1_idx.size < 2 or p2_idx.size < 2:
            return {"valid": False, "reason": "HALF_SIZE_LT_2", "delta_cf": float("nan"), "r_cf_oracle": float("nan")}

        def _half_pt(h_idx: np.ndarray) -> np.ndarray:
            out_h = np.empty(56, dtype=np.float64)
            for b in range(56):
                out_h[b] = _compute_orth_vectors_and_boot(
                    delta_obs[h_idx], delta_pred_56[b, h_idx], delta_transfer[h_idx],
                    target_indices[h_idx], delta_half_a[h_idx], delta_half_b[h_idx],
                    half_valid[h_idx], draws=None, source_grand_mean=source_grand_mean,
                )["estimate"]
            return out_h

        rp1 = _half_pt(p1_idx)
        rp2 = _half_pt(p2_idx)
        if not (np.all(np.isfinite(rp1)) and np.all(np.isfinite(rp2))):
            return {"valid": False, "delta_cf": float("nan"), "r_cf_oracle": float("nan")}

        b_star_1 = int(np.argmax(rp1))
        b_dag_1 = int(np.argmin(rp1))
        b_star_2 = int(np.argmax(rp2))
        b_dag_2 = int(np.argmin(rp2))

        d_12 = float(rp2[b_star_1] - rp2[b_dag_1])
        d_21 = float(rp1[b_star_2] - rp1[b_dag_2])
        delta_cf = 0.5 * (d_12 + d_21)
        r_oracle = 0.5 * (float(rp2[b_star_1]) + float(rp1[b_star_2]))

        res_sp: Dict[str, Any] = {
            "valid": True,
            "b_star_1": b_star_1,
            "b_dag_1": b_dag_1,
            "b_star_2": b_star_2,
            "b_dag_2": b_dag_2,
            "delta_1_to_2": d_12,
            "delta_2_to_1": d_21,
            "delta_cf": float(delta_cf),
            "r_cf_oracle": float(r_oracle),
            "p1_indices": p1_idx,
            "p2_indices": p2_idx,
        }
        if not with_boot:
            return res_sp

        dr2 = np.random.default_rng(int(seed)).integers(0, p2_idx.size, size=(int(n_boot), p2_idx.size))
        dr1 = np.random.default_rng(int(seed) + 10000).integers(0, p1_idx.size, size=(int(n_boot), p1_idx.size))

        s2_star1 = _compute_orth_vectors_and_boot(
            delta_obs[p2_idx], delta_pred_56[b_star_1, p2_idx], delta_transfer[p2_idx],
            target_indices[p2_idx], delta_half_a[p2_idx], delta_half_b[p2_idx],
            half_valid[p2_idx], draws=dr2, source_grand_mean=source_grand_mean,
        )
        s2_dag1 = _compute_orth_vectors_and_boot(
            delta_obs[p2_idx], delta_pred_56[b_dag_1, p2_idx], delta_transfer[p2_idx],
            target_indices[p2_idx], delta_half_a[p2_idx], delta_half_b[p2_idx],
            half_valid[p2_idx], draws=dr2, source_grand_mean=source_grand_mean,
        )
        s1_star2 = _compute_orth_vectors_and_boot(
            delta_obs[p1_idx], delta_pred_56[b_star_2, p1_idx], delta_transfer[p1_idx],
            target_indices[p1_idx], delta_half_a[p1_idx], delta_half_b[p1_idx],
            half_valid[p1_idx], draws=dr1, source_grand_mean=source_grand_mean,
        )
        s1_dag2 = _compute_orth_vectors_and_boot(
            delta_obs[p1_idx], delta_pred_56[b_dag_2, p1_idx], delta_transfer[p1_idx],
            target_indices[p1_idx], delta_half_a[p1_idx], delta_half_b[p1_idx],
            half_valid[p1_idx], draws=dr1, source_grand_mean=source_grand_mean,
        )

        delta_boot = 0.5 * ((s2_star1["boot"] - s2_dag1["boot"]) + (s1_star2["boot"] - s1_dag2["boot"]))
        oracle_boot = 0.5 * (s2_star1["boot"] + s1_star2["boot"])
        fin_d = delta_boot[np.isfinite(delta_boot)]
        fin_o = oracle_boot[np.isfinite(oracle_boot)]
        vf_cf = float(np.mean(np.isfinite(delta_boot))) if delta_boot.size > 0 else 0.0
        vf_or = float(np.mean(np.isfinite(oracle_boot))) if oracle_boot.size > 0 else 0.0
        half_rho_ok = (
            np.isfinite(s2_star1["rho_ab_orth"]) and s2_star1["rho_ab_orth"] > 0.0
            and np.isfinite(s1_star2["rho_ab_orth"]) and s1_star2["rho_ab_orth"] > 0.0
        )

        res_sp.update(
            {
                "delta_cf_boot": delta_boot,
                "oracle_boot": oracle_boot,
                "delta_cf_valid_frac": vf_cf,
                "oracle_valid_frac": vf_or,
                "half_rho_ok": bool(half_rho_ok),
                "p_delta": _bootstrap_two_sided_p(fin_d),
                "p1_boot_draws": dr1,
                "p2_boot_draws": dr2,
            }
        )
        return res_sp

    s0 = _eval_split(0, with_boot=True)
    # Multi-split diagnostic (s = 0..n_splits_diag-1)
    multi_deltas: List[float] = []
    for s_i in range(int(n_splits_diag)):
        r_si = s0 if s_i == 0 else _eval_split(s_i, with_boot=False)
        if r_si.get("valid") and np.isfinite(r_si["delta_cf"]):
            multi_deltas.append(float(r_si["delta_cf"]))
    if multi_deltas:
        mq25, mq50, mq75 = [float(v) for v in np.percentile(multi_deltas, [25.0, 50.0, 75.0])]
        multi_summary = {
            "n_splits": int(n_splits_diag),
            "n_valid": int(len(multi_deltas)),
            "mean": float(np.mean(multi_deltas)),
            "median": mq50,
            "min": float(np.min(multi_deltas)),
            "max": float(np.max(multi_deltas)),
            "iqr": mq75 - mq25,
            "values": multi_deltas,
        }
    else:
        multi_summary = {"n_splits": int(n_splits_diag), "n_valid": 0, "values": []}

    # Floor-effect detectable-signal predicate (`J2`, fixed 98.333% level = 0.05/3)
    base_alpha = 0.05 / 3.0
    l1_lo9833, l1_hi9833, l1_se = _percentile_ci(l1_res["boot"], alpha=base_alpha)
    l1_valid = bool(l1_res["valid_frac"] >= MIN_BOOT_VALID_FRAC and np.isfinite(l1_lo9833))

    if s0.get("valid"):
        or_lo9833, or_hi9833, or_se = _percentile_ci(s0["oracle_boot"], alpha=base_alpha)
        or_valid = bool(s0.get("oracle_valid_frac", 0.0) >= MIN_BOOT_VALID_FRAC and np.isfinite(or_lo9833))
    else:
        or_lo9833 = or_hi9833 = or_se = float("nan")
        or_valid = False

    has_detectable_signal = bool(
        (or_valid and or_lo9833 > 0.0) or (l1_valid and l1_lo9833 > 0.0)
    )

    # Validity gate (`C1.6`, `C1.7`)
    valid_all = bool(
        fin56
        and np.isfinite(rho_orth)
        and rho_orth > 0.0
        and per_index_valid_ok
        and joint_valid_frac >= MIN_BOOT_VALID_FRAC
        and s0.get("valid", False)
        and s0.get("delta_cf_valid_frac", 0.0) >= MIN_BOOT_VALID_FRAC
        and s0.get("half_rho_ok", False)
        and np.isfinite(v_noise)
        and np.isfinite(tau_batch)
    )

    p_delta = float(s0.get("p_delta", 1.0)) if s0.get("valid") else 1.0
    p_hj1 = max(float(p_var), p_delta) if valid_all else 1.0

    verdict_adj = classify_h_j1_under_alpha(
        valid_all=valid_all,
        p_var=float(p_var),
        delta_cf_est=float(s0.get("delta_cf", float("nan"))),
        delta_cf_boot=s0.get("delta_cf_boot", np.zeros(0)),
        tau_batch=float(tau_batch),
        has_detectable_signal=has_detectable_signal,
        alpha_adj=alpha_adj,
    )

    return {
        "r_obs_56": r_obs.tolist(),
        "r_boot_se_56": np.asarray([float(np.std(r_boot[np.isfinite(r_boot[:, b]), b], ddof=1)) if np.sum(np.isfinite(r_boot[:, b])) > 1 else float("nan") for b in range(56)]).tolist(),
        "valid_fracs_56": vfracs.tolist(),
        "joint_valid_frac": joint_valid_frac,
        "rho_ab_orth": rho_orth,
        "landscape": {
            "min": r_min,
            "q10": q10,
            "q25": q25,
            "median": q50,
            "mean": r_mean,
            "q75": q75,
            "q90": q90,
            "max": r_max,
            "iqr": q75 - q25 if np.isfinite(q75) and np.isfinite(q25) else float("nan"),
            "s_obs": s_obs,
            "range_obs": r_max - r_min if np.isfinite(r_max) and np.isfinite(r_min) else float("nan"),
            "b0_score": float(r_obs[0]) if fin56 else float("nan"),
            "b0_rank_desc": rank_b0,
            "b0_percentile": pct_b0,
            "l1_ens56": {
                "estimate": float(l1_res["estimate"]),
                "ci9833_lo": l1_lo9833,
                "ci9833_hi": l1_hi9833,
                "boot_se": l1_se,
                "valid_frac": float(l1_res["valid_frac"]),
            },
            "cf_oracle_s0": {
                "estimate": float(s0.get("r_cf_oracle", float("nan"))),
                "ci9833_lo": or_lo9833,
                "ci9833_hi": or_hi9833,
                "boot_se": or_se,
                "valid_frac": float(s0.get("oracle_valid_frac", 0.0)) if s0.get("valid") else 0.0,
            },
        },
        "v_noise": v_noise,
        "s2_obs": s2_obs,
        "tau_batch": tau_batch,
        "p_var": float(p_var),
        "p_delta": p_delta,
        "p_hj1_iu": p_hj1,
        "cross_fit_s0": {
            "b_star_1": s0.get("b_star_1"),
            "b_dag_1": s0.get("b_dag_1"),
            "b_star_2": s0.get("b_star_2"),
            "b_dag_2": s0.get("b_dag_2"),
            "delta_1_to_2": s0.get("delta_1_to_2"),
            "delta_2_to_1": s0.get("delta_2_to_1"),
            "delta_cf_0": float(s0.get("delta_cf", float("nan"))),
            "delta_cf_valid_frac": float(s0.get("delta_cf_valid_frac", 0.0)) if s0.get("valid") else 0.0,
            "r_cf_oracle_0": float(s0.get("r_cf_oracle", float("nan"))),
        },
        "multi_split_diag": multi_summary,
        "has_detectable_signal": has_detectable_signal,
        "valid_all": valid_all,
        **verdict_adj,
        "_internal": {
            "r_boot_56": r_boot,
            "draws_all": draws_all,
            "s0": s0,
            "l1_boot": l1_res["boot"],
            "l1_est": float(l1_res["estimate"]),
            "l1_valid_frac": float(l1_res["valid_frac"]),
        },
    }


def classify_h_j1_under_alpha(
    valid_all: bool,
    p_var: float,
    delta_cf_est: float,
    delta_cf_boot: np.ndarray,
    tau_batch: float,
    has_detectable_signal: bool,
    alpha_adj: float,
) -> Dict[str, Any]:
    """
    Exact 6-way `H-J1` classifier (`PREREG_J §3.2.5`, `C1.5`, `C1.7`) under any `alpha_adj`:
      - Returns `classification` + explicit `subreason` (`VALIDITY_FAIL` / `INSUFFICIENT_PRECISION` when `INCONCLUSIVE`).
    """
    ci_adj_lo, ci_adj_hi, dcf_se = _percentile_ci(delta_cf_boot, alpha=alpha_adj)
    ci95_lo, ci95_hi, _ = _percentile_ci(delta_cf_boot, alpha=0.05)

    if not valid_all or not (np.isfinite(ci_adj_lo) and np.isfinite(ci95_hi)):
        return {
            "classification": "INCONCLUSIVE",
            "subreason": "VALIDITY_FAIL",
            "alpha_adj": float(alpha_adj),
            "ci_adj_lo": ci_adj_lo,
            "ci_adj_hi": ci_adj_hi,
            "ci95_lo": ci95_lo,
            "ci95_hi": ci95_hi,
            "delta_cf_boot_se": dcf_se,
        }

    detected = bool(np.isfinite(p_var) and p_var < float(alpha_adj) and ci_adj_lo > 0.0)
    if detected:
        if float(delta_cf_est) > TOST_DELTA_R17 and float(tau_batch) > TAU_BATCH_THRESH:
            cls = "BATCH_DEPENDENCE_PRACTICALLY_MEANINGFUL"
        else:
            cls = "BATCH_DEPENDENCE_DETECTED_SUB_DELTA"
        sub = None
    elif not has_detectable_signal:
        cls = "INCONCLUSIVE_LOW_SIGNAL"
        sub = None
    elif ci95_hi < TOST_DELTA_R17 and float(tau_batch) <= TAU_BATCH_THRESH:
        cls = "BATCH_INVARIANT_WITHIN_DELTA"
        sub = None
    else:
        cls = "INCONCLUSIVE"
        sub = "INSUFFICIENT_PRECISION"

    return {
        "classification": cls,
        "subreason": sub,
        "alpha_adj": float(alpha_adj),
        "ci_adj_lo": ci_adj_lo,
        "ci_adj_hi": ci_adj_hi,
        "ci95_lo": ci95_lo,
        "ci95_hi": ci95_hi,
        "delta_cf_boot_se": dcf_se,
    }


def compute_h_j2_for_subset(
    delta_obs: np.ndarray,
    delta_pred_56: np.ndarray,
    fct_delta_pred: np.ndarray,
    delta_transfer: np.ndarray,
    target_indices: Sequence[int],
    delta_half_a: np.ndarray,
    delta_half_b: np.ndarray,
    half_valid: np.ndarray,
    has_detectable_signal: bool,
    fct_est: Optional[float] = None,
    r_obs_56: Optional[np.ndarray] = None,
    r_boot_56: Optional[np.ndarray] = None,
    fct_boot: Optional[np.ndarray] = None,
    n_boot: int = 4000,
    n_boot_bonf56: int = 20000,
    seed: int = DEFAULT_R17_SEED,
    source_grand_mean: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """
    Compute `H-J2` (`PREREG_J §3.3` + `C1.8`):
      - Primary per-index `97.5%` CI (`PREREG_I §3.2`) calls across `b = 0..55`
      - `f_flip`, counts, `has_full_reversal`, and 5-way `J2` classification (`CALL_STABLE` vs `INCONCLUSIVE_LOW_SIGNAL`)
      - Companion Bonferroni-56 hardened flip check (`C1.8`: normal `z = norm.ppf(1 - 0.025/112)` AND
        empirical `20,000`-draw `np.random.default_rng(20261701)` quantiles `0.025/112` & `1 - 0.025/112`).
    """
    p_cnt = delta_obs.shape[0]
    if r_obs_56 is None or r_boot_56 is None or fct_boot is None or fct_est is None:
        sweep = _compute_56_sweep_and_fct_boot(
            delta_obs, delta_pred_56, fct_delta_pred, delta_transfer, target_indices,
            delta_half_a, delta_half_b, half_valid, n_boot=n_boot, seed=seed,
            source_grand_mean=source_grand_mean,
        )
        r_obs_56 = sweep["r_obs_56"]
        r_boot_56 = sweep["r_boot_56"]
        fct_est = sweep["fct_est"]
        fct_boot = sweep["fct_boot"]

    def _call_vec(diff_est: float, lo: float, hi: float, vf: float) -> str:
        if not (np.isfinite(vf) and vf >= MIN_BOOT_VALID_FRAC and np.isfinite(lo) and np.isfinite(hi)):
            return "INDISTINGUISHABLE"
        if lo > 0.0:
            return "MODEL_SUPERIOR_TO_FCT"
        if hi < 0.0:
            return "FCT_SUPERIOR_TO_MODEL"
        return "INDISTINGUISHABLE"

    z_bonf = float(norm.ppf(1.0 - 0.025 / 112.0))
    calls_975: List[str] = []
    calls_bonf_norm: List[str] = []
    diffs_56: List[float] = []
    ci975_56: List[Tuple[float, float]] = []

    for b in range(56):
        d_est = float(r_obs_56[b] - fct_est)
        d_boot = r_boot_56[:, b] - fct_boot
        vf = float(np.mean(np.isfinite(d_boot))) if d_boot.size > 0 else 0.0
        lo975, hi975, se = _percentile_ci(d_boot, alpha=0.025)
        calls_975.append(_call_vec(d_est, lo975, hi975, vf))
        lo_n = d_est - z_bonf * se if np.isfinite(se) else float("nan")
        hi_n = d_est + z_bonf * se if np.isfinite(se) else float("nan")
        calls_bonf_norm.append(_call_vec(d_est, lo_n, hi_n, vf))
        diffs_56.append(d_est)
        ci975_56.append((lo975, hi975))

    # Empirical 20,000-draw Bonferroni-56 (`C1.8`)
    rng_20k = np.random.default_rng(int(seed))
    draws_20k = rng_20k.integers(0, p_cnt, size=(int(n_boot_bonf56), p_cnt))
    sw_20k = _compute_56_sweep_and_fct_boot(
        delta_obs, delta_pred_56, fct_delta_pred, delta_transfer, target_indices,
        delta_half_a, delta_half_b, half_valid, n_boot=n_boot_bonf56, seed=seed,
        source_grand_mean=source_grand_mean, draws_override=draws_20k,
    )
    q_bonf_lo = 100.0 * (0.025 / 112.0)
    q_bonf_hi = 100.0 * (1.0 - 0.025 / 112.0)
    calls_bonf_emp: List[str] = []
    for b in range(56):
        d_est_20k = float(sw_20k["r_obs_56"][b] - sw_20k["fct_est"])
        d_20k = sw_20k["r_boot_56"][:, b] - sw_20k["fct_boot"]
        fin_20k = d_20k[np.isfinite(d_20k)]
        vf_20k = float(np.mean(np.isfinite(d_20k))) if d_20k.size > 0 else 0.0
        if fin_20k.size > 0:
            lo_e, hi_e = np.percentile(fin_20k, [q_bonf_lo, q_bonf_hi])
        else:
            lo_e = hi_e = float("nan")
        calls_bonf_emp.append(_call_vec(d_est_20k, float(lo_e), float(hi_e), vf_20k))

    def _summarize_calls(calls: List[str]) -> Dict[str, Any]:
        c0 = calls[0]
        n_flip = sum(1 for c in calls if c != c0)
        idx_mod = [i for i, c in enumerate(calls) if c == "MODEL_SUPERIOR_TO_FCT"]
        idx_ind = [i for i, c in enumerate(calls) if c == "INDISTINGUISHABLE"]
        idx_fct = [i for i, c in enumerate(calls) if c == "FCT_SUPERIOR_TO_MODEL"]
        return {
            "b0_call": c0,
            "n_flip": int(n_flip),
            "f_flip": float(n_flip) / 56.0,
            "n_MODEL_SUPERIOR": len(idx_mod),
            "n_INDISTINGUISHABLE": len(idx_ind),
            "n_FCT_SUPERIOR": len(idx_fct),
            "has_full_reversal": bool(len(idx_mod) >= 1 and len(idx_fct) >= 1),
            "indices_MODEL_SUPERIOR": idx_mod,
            "indices_INDISTINGUISHABLE": idx_ind,
            "indices_FCT_SUPERIOR": idx_fct,
            "calls_56": calls,
        }

    s_975 = _summarize_calls(calls_975)
    s_bn = _summarize_calls(calls_bonf_norm)
    s_be = _summarize_calls(calls_bonf_emp)

    if s_975["has_full_reversal"]:
        cls = "CALL_UNSTABLE_SIGN_REVERSAL"
    elif s_975["n_flip"] >= 5:
        cls = "CALL_UNSTABLE_SUBSTANTIAL"
    elif s_975["n_flip"] >= 1:
        cls = "CALL_UNSTABLE_MINOR"
    elif has_detectable_signal:
        cls = "CALL_STABLE"
    else:
        cls = "INCONCLUSIVE_LOW_SIGNAL"

    return {
        "classification": cls,
        "has_detectable_signal": bool(has_detectable_signal),
        "fct_r_orth_cal": float(fct_est),
        "delta_vs_fct_56": diffs_56,
        "ci975_56": ci975_56,
        "primary_975": s_975,
        "bonf56_normal": {"z_critical": z_bonf, **s_bn},
        "bonf56_empirical_20k": {"n_boot": int(n_boot_bonf56), **s_be},
    }


def select_l2_ntc_index_and_diagnostics(
    ntc_pb_ensemble_56: np.ndarray,
    obs_ntc_mean: np.ndarray,
    ntc_pb_per_fold_56: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """
    Select label-free `b_NTC` via primary Pearson correlation `C_NTC(b)` (`PREREG_J §3.4` + Ruling 3)
    with negative L2 criterion, `top1 - top2` margins, `b=0` rank, and per-fold diagnostics.
    """
    obs_v = np.asarray(obs_ntc_mean, dtype=np.float64)
    ens = np.asarray(ntc_pb_ensemble_56, dtype=np.float64)

    def _score_criteria(mat_56: np.ndarray) -> Dict[str, Any]:
        c_p = np.empty(56, dtype=np.float64)
        c_neg_l2 = np.empty(56, dtype=np.float64)
        oc = obs_v - obs_v.mean()
        o_norm = float(np.linalg.norm(oc))
        for b in range(56):
            v = mat_56[b]
            vc = v - v.mean()
            den = float(np.linalg.norm(vc)) * o_norm
            c_p[b] = float(np.dot(vc, oc) / den) if den > EPS else 0.0
            c_neg_l2[b] = -float(np.linalg.norm(v - obs_v))

        ord_p = np.argsort(-c_p, kind="mergesort")
        t1, t2 = int(ord_p[0]), int(ord_p[1])
        ord_l2 = np.argsort(-c_neg_l2, kind="mergesort")
        l1, l2 = int(ord_l2[0]), int(ord_l2[1])
        q25, q50, q75 = [float(v) for v in np.percentile(c_p, [25.0, 50.0, 75.0])]
        return {
            "top1_idx": t1,
            "criterion_top1": float(c_p[t1]),
            "top2_idx": t2,
            "criterion_top2": float(c_p[t2]),
            "criterion_margin_top1_minus_top2": float(c_p[t1] - c_p[t2]),
            "b0_rank_desc": int(np.flatnonzero(ord_p == 0)[0]) + 1,
            "b0_criterion": float(c_p[0]),
            "c_ntc_min": float(np.min(c_p)),
            "c_ntc_median": q50,
            "c_ntc_iqr": q75 - q25,
            "c_ntc_max": float(np.max(c_p)),
            "c_ntc_all_56": c_p.tolist(),
            "neg_l2_diag": {
                "top1_idx": l1,
                "criterion_top1": float(c_neg_l2[l1]),
                "top2_idx": l2,
                "criterion_top2": float(c_neg_l2[l2]),
                "criterion_margin_top1_minus_top2": float(c_neg_l2[l1] - c_neg_l2[l2]),
                "c_neg_l2_all_56": c_neg_l2.tolist(),
            },
        }

    primary = _score_criteria(ens)
    per_fold_diag: Dict[str, Any] = {}
    if ntc_pb_per_fold_56 is not None and ntc_pb_per_fold_56.shape[0] > 0:
        pf = np.asarray(ntc_pb_per_fold_56, dtype=np.float64)
        fold_idxs: List[int] = []
        for f_i in range(pf.shape[0]):
            s_f = _score_criteria(pf[f_i])
            per_fold_diag[f"fold_{f_i}"] = s_f
            fold_idxs.append(int(s_f["top1_idx"]))
        primary["per_fold_top1_indices"] = fold_idxs
        primary["per_fold_all_agree"] = bool(len(set(fold_idxs)) == 1)
        primary["per_fold_agree_with_ensemble"] = int(sum(1 for idx in fold_idxs if idx == primary["top1_idx"]))
        primary["per_fold_details"] = per_fold_diag
    return primary


def classify_h_j3_rule(
    rule_name: str,
    delta_est: float,
    delta_boot: np.ndarray,
    alpha_adj: float,
    selected_idx: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Classify `H-J3a` (`L1_ens56`) or `H-J3b` (`L2_ntc_recon_match`) vs `b=0` (`PREREG_J §3.4.3`, `C1.5`).
    """
    is_exact_zero_tie = bool(rule_name == "L2_ntc_recon_match" and selected_idx == 0)
    vf = float(np.mean(np.isfinite(delta_boot))) if delta_boot.size > 0 else 0.0
    fin = delta_boot[np.isfinite(delta_boot)]
    lo, hi, se = _percentile_ci(delta_boot, alpha=alpha_adj)
    p_val = 1.0 if is_exact_zero_tie else _bootstrap_two_sided_p(fin)

    if is_exact_zero_tie:
        cls = "EXACT_TIE_SELECTED_ZERO"
    elif not (vf >= MIN_BOOT_VALID_FRAC and np.isfinite(lo) and np.isfinite(hi)):
        cls = "INDISTINGUISHABLE"
    elif lo > 0.0:
        cls = "RULE_SUPERIOR_TO_ZERO"
    elif hi < 0.0:
        cls = "ZERO_SUPERIOR_TO_RULE"
    else:
        cls = "INDISTINGUISHABLE"

    return {
        "rule": rule_name,
        "selected_idx": selected_idx,
        "classification": cls,
        "practical_gt_0_04": bool(cls == "RULE_SUPERIOR_TO_ZERO" and float(delta_est) > TOST_DELTA_R17),
        "estimate": float(delta_est),
        "ci_adj_lo": lo,
        "ci_adj_hi": hi,
        "boot_se": se,
        "valid_frac": vf,
        "p_value": float(p_val),
        "alpha_adj": float(alpha_adj),
    }


def apply_holm_fwer_r17(
    p_entries: List[Dict[str, Any]],
    alpha_family: float = 0.05,
) -> Dict[str, Any]:
    """
    Apply Holm-Bonferroni step-down across `M` confirmatory tests (`C1.5`, `C2.1`):
      - Step-down runs on raw p-values only (`p_raw < alpha_family / (M - rank)`).
      - Each test gets `alpha_eff`: its step `alpha` if Holm rejects it, else `stopping_alpha`
        (the `alpha_step` at which the p-value step-down first failed).
      - A test is declared `rejected` iff Holm rejects it AND its interval gate passes at
        `alpha_eff`. Interval-gate failure removes only that test without halting step-down.
    """
    m = len(p_entries)
    if m == 0:
        return {
            "m_tests": 0,
            "alpha_family": float(alpha_family),
            "stopping_alpha": None,
            "steps": [],
            "by_key": {},
        }
    ordered = sorted(enumerate(p_entries), key=lambda kv: (float(kv[1]["p_value"]), str(kv[1]["key"])))

    # Pass 1: Holm step-down on raw p-values alone
    holm_p_still_rejecting = True
    stopping_alpha: Optional[float] = None
    running_max_p_holm = 0.0
    step_alphas: List[float] = []
    p_holms: List[float] = []
    p_holm_rejected: List[bool] = []

    for rank, (_, item) in enumerate(ordered):
        alpha_step = float(alpha_family) / float(m - rank)
        p_raw = float(item["p_value"])
        p_h = min(1.0, p_raw * float(m - rank))
        running_max_p_holm = max(running_max_p_holm, p_h)
        p_hol_rej = bool(holm_p_still_rejecting and np.isfinite(p_raw) and p_raw < alpha_step)
        if not p_hol_rej and holm_p_still_rejecting:
            holm_p_still_rejecting = False
            stopping_alpha = alpha_step
        step_alphas.append(alpha_step)
        p_holms.append(float(running_max_p_holm))
        p_holm_rejected.append(p_hol_rej)

    # Pass 2: Evaluate each test's own interval gate at its effective alpha
    steps: List[Dict[str, Any]] = []
    by_key: Dict[str, Dict[str, Any]] = {}
    for rank, (_, item) in enumerate(ordered):
        p_hol_rej = p_holm_rejected[rank]
        alpha_eff = step_alphas[rank] if p_hol_rej else float(stopping_alpha)  # type: ignore[arg-type]
        ci_ok = (
            True
            if "ci_check_fn" not in item or item["ci_check_fn"] is None
            else bool(item["ci_check_fn"](alpha_eff))
        )
        rej = bool(p_hol_rej and ci_ok)
        rec = {
            "key": item["key"],
            "rank": rank + 1,
            "m_tests": m,
            "p_raw": float(item["p_value"]),
            "p_holm": p_holms[rank],
            "alpha_step": step_alphas[rank],
            "holm_p_rejected": bool(p_hol_rej),
            "alpha_eff": float(alpha_eff),
            "ci_gate_passed": ci_ok,
            "rejected": rej,
        }
        steps.append(rec)
        by_key[item["key"]] = rec

    return {
        "m_tests": m,
        "alpha_family": float(alpha_family),
        "stopping_alpha": stopping_alpha,
        "steps": steps,
        "by_key": by_key,
    }


def _decompose_g4_cell_quad(
    obs: np.ndarray,
    mm: np.ndarray,
    ms: np.ndarray,
    sm: np.ndarray,
    ss: np.ndarray,
    t_mean: np.ndarray,
    t_gidx: np.ndarray,
    grand_mean: np.ndarray,
    ha: np.ndarray,
    hb: np.ndarray,
    hv: np.ndarray,
    n_boot: int,
    seed: int,
) -> Dict[str, Any]:
    """Compute `(Delta_batch, Delta_ctrl, Delta_both, Delta_interaction, batch_share)` under `r_orth_cal` & `r_dorth_cal`."""
    d_batch = compute_paired_contrast_r16(obs, mm, ms, t_mean, t_gidx, None, ha, hb, hv, n_boot=n_boot, seed=seed)
    d_ctrl = compute_paired_contrast_r16(obs, mm, sm, t_mean, t_gidx, None, ha, hb, hv, n_boot=n_boot, seed=seed)
    d_both = compute_paired_contrast_r16(obs, mm, ss, t_mean, t_gidx, None, ha, hb, hv, n_boot=n_boot, seed=seed)
    b_est = float(d_batch["estimate"])
    c_est = float(d_ctrl["estimate"])
    both_est = float(d_both["estimate"])
    den = abs(c_est) + abs(b_est)
    d_batch_d = compute_paired_contrast_r16(obs, mm, ms, t_mean, t_gidx, grand_mean, ha, hb, hv, n_boot=n_boot, seed=seed)
    d_ctrl_d = compute_paired_contrast_r16(obs, mm, sm, t_mean, t_gidx, grand_mean, ha, hb, hv, n_boot=n_boot, seed=seed)
    return {
        "delta_batch_orth": d_batch,
        "delta_ctrl_orth": d_ctrl,
        "delta_both_orth": d_both,
        "delta_interaction_orth": float(both_est - c_est - b_est),
        "batch_share_orth": float(abs(b_est) / den) if den > 1e-12 else float("nan"),
        "delta_batch_dorth": d_batch_d,
        "delta_ctrl_dorth": d_ctrl_d,
    }


def compute_j5a_g4_for_line(
    obs_npz: Dict[str, np.ndarray],
    sweep56_npz: Dict[str, np.ndarray],
    g4_mis_npz: Dict[str, np.ndarray],
    h_j1_s0: Dict[str, Any],
    n_boot: int = 4000,
    seed: int = DEFAULT_R17_SEED,
) -> Dict[str, Any]:
    """
    Compute `J5a` (`PREREG_J §3.5.1` + `C1.11`):
      1. Cross-fitted `(b*_1, b_dag_1)` evaluated strictly on `P2` ONLY (`seed`), + reversed symmetry;
      2. Auxiliary `(b*_2, b_dag_2)` evaluated strictly on `P1` ONLY (`seed + 10000`), + reversed symmetry;
      3. Continuity `(0, 17)` & `(0, 35)` (with reversed pairs) evaluated on ALL `N_eval` perturbations (`seed`).
    """
    obs = obs_npz["delta_obs"]
    ha = obs_npz["delta_half_a"]
    hb = obs_npz["delta_half_b"]
    hv = obs_npz["half_valid"]
    tm = obs_npz["delta_transfer"]
    gm = obs_npz["source_grand_mean"]
    tg = obs_npz["target_indices"]
    d56 = sweep56_npz["delta_pred_56"]

    mis_names = [str(x) for x in g4_mis_npz["mismatched_sources"].tolist()]
    g4_b_arr = [int(x) for x in g4_mis_npz["g4_batch_indices"].tolist()]
    g4_b_pos = {b: i for i, b in enumerate(g4_b_arr)}
    dmis = g4_mis_npz["delta_pred_mismatched"]  # (S_mis, B_g4, P, G)

    seen_s = [str(x) for x in g4_mis_npz["seen_sources"].tolist()]
    unseen_s = [str(x) for x in g4_mis_npz["unseen_sources"].tolist()]
    same_s = [str(x) for x in g4_mis_npz["same_lineage_sources"].tolist()]

    def _run_pair_set(bm: int, bs: int, row_idx: np.ndarray, pair_seed: int) -> Dict[str, Any]:
        bi_m = g4_b_pos[int(bm)]
        bi_s = g4_b_pos[int(bs)]
        mm = d56[int(bm), row_idx]
        ms = d56[int(bs), row_idx]
        per_src: Dict[str, Any] = {}
        rev_src: Dict[str, Any] = {}
        for si, sname in enumerate(mis_names):
            sm = dmis[si, bi_m, row_idx]
            ss = dmis[si, bi_s, row_idx]
            per_src[sname] = _decompose_g4_cell_quad(
                obs[row_idx], mm, ms, sm, ss, tm[row_idx], tg[row_idx], gm, ha[row_idx], hb[row_idx], hv[row_idx],
                n_boot=n_boot, seed=pair_seed,
            )
            rev_src[sname] = _decompose_g4_cell_quad(
                obs[row_idx], ms, mm, ss, sm, tm[row_idx], tg[row_idx], gm, ha[row_idx], hb[row_idx], hv[row_idx],
                n_boot=n_boot, seed=pair_seed,
            )

        def _stratum_sum(sl: Sequence[str]) -> Dict[str, Any]:
            vals = [per_src[s] for s in sl if s in per_src]
            if not vals:
                return {"n_sources": 0}
            shares = [float(v["batch_share_orth"]) for v in vals if np.isfinite(v["batch_share_orth"])]
            batches = [float(v["delta_batch_orth"]["estimate"]) for v in vals]
            ctrls = [float(v["delta_ctrl_orth"]["estimate"]) for v in vals]
            return {
                "n_sources": len(vals),
                "median_batch_share_orth": float(np.median(shares)) if shares else float("nan"),
                "median_delta_batch_orth": float(np.median(batches)),
                "median_delta_ctrl_orth": float(np.median(ctrls)),
            }

        return {
            "b_m": int(bm),
            "b_s": int(bs),
            "n_eval_rows": int(row_idx.size),
            "seed": int(pair_seed),
            "per_source": per_src,
            "reversed_symmetry": rev_src,
            "stratum_summaries": {
                "seen_by_checkpoint": _stratum_sum(seen_s),
                "unseen_by_checkpoint": _stratum_sum(unseen_s),
                "same_lineage_sources": _stratum_sum(same_s),
            },
        }

    all_rows = np.arange(obs.shape[0], dtype=np.int64)
    out: Dict[str, Any] = {
        "continuity_0_17": _run_pair_set(0, 17, all_rows, seed),
        "continuity_0_35": _run_pair_set(0, 35, all_rows, seed),
    }

    p1_idx = obs_npz["split_s0_p1_indices"]
    p2_idx = obs_npz["split_s0_p2_indices"]
    if h_j1_s0.get("valid") and p1_idx.size >= 2 and p2_idx.size >= 2:
        b1_s, b1_d = int(h_j1_s0["b_star_1"]), int(h_j1_s0["b_dag_1"])
        b2_s, b2_d = int(h_j1_s0["b_star_2"]), int(h_j1_s0["b_dag_2"])
        out["cross_fit_primary_on_p2"] = _run_pair_set(b1_s, b1_d, p2_idx, seed)
        out["cross_fit_auxiliary_on_p1"] = _run_pair_set(b2_s, b2_d, p1_idx, int(seed) + 10000)
    else:
        out["cross_fit_primary_on_p2"] = None
        out["cross_fit_auxiliary_on_p1"] = None
    return out


def _spearman_vectorized_against_perms(
    vec_a: np.ndarray,
    vec_b: np.ndarray,
    perm_table: np.ndarray,
) -> Dict[str, Any]:
    """Compute exact Spearman correlation & two-sided `p_perm` using shared `perm_table` (`(N_perm, 56)`)."""
    ra = rankdata(np.asarray(vec_a, dtype=np.float64))
    rb = rankdata(np.asarray(vec_b, dtype=np.float64))
    rac = ra - ra.mean()
    rbc = rb - rb.mean()
    den = math.sqrt(float(np.dot(rac, rac)) * float(np.dot(rbc, rbc)))
    if den <= EPS:
        rho_obs = 0.0
        p_perm = 1.0
    else:
        rho_obs = float(np.dot(rac, rbc) / den)
        perm_rbc = rbc[perm_table]  # (N_perm, 56)
        rho_null = (perm_rbc @ rac) / den
        n_perm = int(perm_table.shape[0])
        p_perm = float(1 + int(np.sum(np.abs(rho_null) >= abs(rho_obs)))) / float(n_perm + 1)

    top5_a = set(np.argsort(-np.asarray(vec_a, dtype=np.float64), kind="mergesort")[:5].tolist())
    top5_b = set(np.argsort(-np.asarray(vec_b, dtype=np.float64), kind="mergesort")[:5].tolist())
    bot5_a = set(np.argsort(np.asarray(vec_a, dtype=np.float64), kind="mergesort")[:5].tolist())
    bot5_b = set(np.argsort(np.asarray(vec_b, dtype=np.float64), kind="mergesort")[:5].tolist())
    return {
        "spearman_rho": float(rho_obs),
        "p_perm": float(p_perm),
        "top5_shared": int(len(top5_a & top5_b)),
        "bot5_shared": int(len(bot5_a & bot5_b)),
    }


def compute_bh_fdr(p_vals: Sequence[float]) -> List[float]:
    """Standard Benjamini-Hochberg step-up adjusted p-values."""
    m = len(p_vals)
    if m == 0:
        return []
    order = np.argsort(np.asarray(p_vals, dtype=np.float64), kind="mergesort")
    adj = np.empty(m, dtype=np.float64)
    prev = 1.0
    for rev_i in range(m - 1, -1, -1):
        idx = int(order[rev_i])
        rank = rev_i + 1
        val = min(prev, float(p_vals[idx]) * float(m) / float(rank))
        adj[idx] = val
        prev = val
    return [float(x) for x in np.clip(adj, 0.0, 1.0)]


def compute_j5b_rank_correlations(
    line_scores_by_model: Dict[str, Dict[str, np.ndarray]],
    n_eval_by_line: Dict[str, int],
    n_perm: int = 10000,
    seed: int = DEFAULT_R17_SEED,
) -> Dict[str, Any]:
    """
    Compute `J5b` (`PREREG_J §3.5.2` + `C1.10`) using ONE shared `(n_perm, 56)` permutation table.
    """
    rng = np.random.default_rng(int(seed))
    perm_table = np.stack([rng.permutation(56) for _ in range(int(n_perm))], axis=0)

    def _run_pair_group(scores_map: Dict[str, np.ndarray], lines: Sequence[str]) -> List[Dict[str, Any]]:
        avail = sorted(ln for ln in lines if ln in scores_map)
        rows: List[Dict[str, Any]] = []
        for i in range(len(avail)):
            for j in range(i + 1, len(avail)):
                l1, l2 = avail[i], avail[j]
                st = _spearman_vectorized_against_perms(scores_map[l1], scores_map[l2], perm_table)
                rows.append({"line_1": l1, "line_2": l2, **st})
        if rows:
            fdrs = compute_bh_fdr([r["p_perm"] for r in rows])
            for r, fv in zip(rows, fdrs):
                r["bh_fdr"] = fv
        return rows

    by_model: Dict[str, Any] = {}
    for mkey, smap in sorted(line_scores_by_model.items()):
        conf_pairs = _run_pair_group(smap, R17_CONFIRMATORY_LINES)
        pfizer_pair = _run_pair_group(smap, ("pfizer_telohaec_basal", "pfizer_telohaec_cytokine"))
        # C1.10: exclude a549_ko_cas9pos from fresh lines with N >= 16
        fresh_ge16 = [
            ln for ln, n_val in sorted(n_eval_by_line.items())
            if ln in smap and n_val >= 16 and ln != "a549_ko_cas9pos"
        ]
        fresh_ge16_indep = [ln for ln in fresh_ge16 if ln != "pfizer_telohaec_cytokine"]
        fresh_pairs = _run_pair_group(smap, fresh_ge16)
        fresh_indep_pairs = _run_pair_group(smap, fresh_ge16_indep)
        by_model[mkey] = {
            "confirmatory_pairs": conf_pairs,
            "pfizer_within_family": pfizer_pair,
            "fresh_n_ge_16_lines": fresh_ge16,
            "fresh_n_ge_16_pairs": fresh_pairs,
            "fresh_n_ge_16_excl_cytokine_pairs": fresh_indep_pairs,
        }

    cross_model: List[Dict[str, Any]] = []
    if "STATE_SE" in line_scores_by_model and "STATE_HVG" in line_scores_by_model:
        common_lines = sorted(
            set(line_scores_by_model["STATE_SE"].keys()) & set(line_scores_by_model["STATE_HVG"].keys())
        )
        for ln in common_lines:
            st = _spearman_vectorized_against_perms(
                line_scores_by_model["STATE_SE"][ln],
                line_scores_by_model["STATE_HVG"][ln],
                perm_table,
            )
            cross_model.append({"line": ln, **st})
        if cross_model:
            fdrs = compute_bh_fdr([r["p_perm"] for r in cross_model])
            for r, fv in zip(cross_model, fdrs):
                r["bh_fdr"] = fv

    return {
        "n_perm": int(n_perm),
        "seed": int(seed),
        "by_model": by_model,
        "cross_model_within_line": cross_model,
    }


def evaluate_line_model_from_arrays(
    obs_npz: Dict[str, np.ndarray],
    sweep56_npz: Dict[str, np.ndarray],
    g4_mis_npz: Optional[Dict[str, np.ndarray]] = None,
    n_boot: int = 4000,
    n_boot_bonf56: int = 20000,
    n_splits_diag: int = 20,
    seed: int = DEFAULT_R17_SEED,
    alpha_base: float = 0.05 / 3.0,
) -> Dict[str, Any]:
    """
    Evaluate one `(line, model)` purely from `observed.npz`, `preds_<M>_sweep56.npz`,
    and optional `preds_<M>_g4_mismatched.npz` (`Stage 4` zero-leakage recomputation).
    """
    d_obs = obs_npz["delta_obs"]
    ha = obs_npz["delta_half_a"]
    hb = obs_npz["delta_half_b"]
    hv = obs_npz["half_valid"]
    tm = obs_npz["delta_transfer"]
    gm = obs_npz["source_grand_mean"]
    tg = obs_npz["target_indices"]
    fct_d = obs_npz["fct_delta_pred"]
    obs_ntc_mean = obs_npz["obs_ntc_mean_official_median_log1p"]

    d56 = sweep56_npz["delta_pred_56"]
    d_l1 = sweep56_npz["delta_pred_L1_ens56"]
    ntc_56 = sweep56_npz["ntc_pb_ensemble_56"]
    ntc_pf = sweep56_npz.get("ntc_pb_per_fold_56")

    h_j1 = compute_h_j1_for_subset(
        delta_obs=d_obs,
        delta_pred_56=d56,
        delta_pred_L1=d_l1,
        delta_transfer=tm,
        target_indices=tg,
        delta_half_a=ha,
        delta_half_b=hb,
        half_valid=hv,
        n_boot=n_boot,
        seed=seed,
        n_splits_diag=n_splits_diag,
        source_grand_mean=None,
        alpha_adj=alpha_base,
    )

    internal = h_j1.pop("_internal")
    r_obs_56 = np.asarray(h_j1["r_obs_56"], dtype=np.float64)
    r_boot_56 = internal["r_boot_56"]
    draws_all = internal["draws_all"]
    s0 = internal["s0"]

    sw_fct = _compute_orth_vectors_and_boot(
        d_obs, fct_d, tm, tg, ha, hb, hv, draws=draws_all, source_grand_mean=None
    )
    h_j2 = compute_h_j2_for_subset(
        delta_obs=d_obs,
        delta_pred_56=d56,
        fct_delta_pred=fct_d,
        delta_transfer=tm,
        target_indices=tg,
        delta_half_a=ha,
        delta_half_b=hb,
        half_valid=hv,
        has_detectable_signal=bool(h_j1["has_detectable_signal"]),
        fct_est=float(sw_fct["estimate"]),
        r_obs_56=r_obs_56,
        r_boot_56=r_boot_56,
        fct_boot=sw_fct["boot"],
        n_boot=n_boot,
        n_boot_bonf56=n_boot_bonf56,
        seed=seed,
        source_grand_mean=None,
    )

    # H-J3a & H-J3b
    l2_sel = select_l2_ntc_index_and_diagnostics(ntc_56, obs_ntc_mean, ntc_pb_per_fold_56=ntc_pf)
    b_ntc = int(l2_sel["top1_idx"])
    # Verify schema consistency if array carried pre-computed l2_selected_batch_idx
    if "l2_selected_batch_idx" in sweep56_npz:
        arr_b = int(sweep56_npz["l2_selected_batch_idx"])
        if arr_b != b_ntc:
            raise ValueError(f"L2 index mismatch: recomputed {b_ntc} != exported {arr_b}")

    r_b0_est = float(r_obs_56[0])
    r_b0_boot = r_boot_56[:, 0]
    d_l1_vs_0_est = float(internal["l1_est"] - r_b0_est)
    d_l1_vs_0_boot = internal["l1_boot"] - r_b0_boot
    h_j3a = classify_h_j3_rule(
        "L1_ens56", d_l1_vs_0_est, d_l1_vs_0_boot, alpha_adj=alpha_base, selected_idx=None
    )

    d_l2_vs_0_est = float(r_obs_56[b_ntc] - r_b0_est)
    d_l2_vs_0_boot = r_boot_56[:, b_ntc] - r_b0_boot
    h_j3b = classify_h_j3_rule(
        "L2_ntc_recon_match", d_l2_vs_0_est, d_l2_vs_0_boot, alpha_adj=alpha_base, selected_idx=b_ntc
    )

    # Rule vs FCT (97.5% CI matching H-I2) and out-of-sample vs CF-Oracle on s=0
    def _rule_extra(d_rule: np.ndarray, rule_boot: np.ndarray, rule_est: float) -> Dict[str, Any]:
        vs_fct_boot = rule_boot - sw_fct["boot"]
        lo975, hi975, _ = _percentile_ci(vs_fct_boot, 0.025)
        vf = float(np.mean(np.isfinite(vs_fct_boot))) if vs_fct_boot.size > 0 else 0.0
        if vf >= MIN_BOOT_VALID_FRAC and lo975 > 0.0:
            c_fct = "MODEL_SUPERIOR_TO_FCT"
        elif vf >= MIN_BOOT_VALID_FRAC and hi975 < 0.0:
            c_fct = "FCT_SUPERIOR_TO_MODEL"
        else:
            c_fct = "INDISTINGUISHABLE"

        vs_cf = None
        if s0.get("valid"):
            p1, p2 = s0["p1_indices"], s0["p2_indices"]
            dr1, dr2 = s0["p1_boot_draws"], s0["p2_boot_draws"]
            r_r_p2 = _compute_orth_vectors_and_boot(
                d_obs[p2], d_rule[p2], tm[p2], tg[p2], ha[p2], hb[p2], hv[p2], draws=dr2
            )
            r_r_p1 = _compute_orth_vectors_and_boot(
                d_obs[p1], d_rule[p1], tm[p1], tg[p1], ha[p1], hb[p1], hv[p1], draws=dr1
            )
            oos_est = 0.5 * (float(r_r_p2["estimate"]) + float(r_r_p1["estimate"])) - float(s0["r_cf_oracle"])
            oos_boot = 0.5 * (r_r_p2["boot"] + r_r_p1["boot"]) - s0["oracle_boot"]
            lo_cf, hi_cf, _ = _percentile_ci(oos_boot, alpha_base)
            vs_cf = {"estimate": float(oos_est), "ci9833_lo": lo_cf, "ci9833_hi": hi_cf}

        return {
            "r_orth_cal": float(rule_est),
            "vs_fct": {"estimate": float(rule_est - sw_fct["estimate"]), "ci975_lo": lo975, "ci975_hi": hi975, "call": c_fct},
            "vs_cf_oracle_s0": vs_cf,
        }

    h_j3a.update(_rule_extra(d_l1, internal["l1_boot"], internal["l1_est"]))
    h_j3b.update(_rule_extra(d56[b_ntc], r_boot_56[:, b_ntc], float(r_obs_56[b_ntc])))

    # Secondary neutral-batch comparators
    comp_res: Dict[str, Any] = {}
    for cname, ckey in (
        ("bypass_zero", "delta_pred_bypass_zero"),
        ("mean_batch_emb", "delta_pred_mean_batch_emb"),
        ("multi_batch_ens8", "delta_pred_multi_batch_ens8"),
        ("medoid_idx", "delta_pred_medoid_idx"),
    ):
        if ckey in sweep56_npz:
            rc = _compute_orth_vectors_and_boot(
                d_obs, sweep56_npz[ckey], tm, tg, ha, hb, hv, draws=draws_all
            )
            lo975, hi975, se_c = _percentile_ci(rc["boot"], 0.025)
            comp_res[cname] = {
                "estimate": float(rc["estimate"]),
                "ci975_lo": lo975,
                "ci975_hi": hi975,
                "boot_se": se_c,
                "valid_frac": float(rc["valid_frac"]),
            }

    # Secondary r_dorth_cal summary on b=0, L1, and 56-index point estimates
    dorth_r_obs = np.asarray([
        _compute_orth_vectors_and_boot(
            d_obs, d56[b], tm, tg, ha, hb, hv, draws=None, source_grand_mean=gm
        )["estimate"]
        for b in range(56)
    ], dtype=np.float64)
    dorth_l1 = compute_orth_suite_r16(
        d_obs, d_l1, tm, tg, source_grand_mean=gm, delta_half_a=ha, delta_half_b=hb, half_valid=hv, n_boot=n_boot, seed=seed
    )

    j5a_res = (
        compute_j5a_g4_for_line(obs_npz, sweep56_npz, g4_mis_npz, s0, n_boot=n_boot, seed=seed)
        if g4_mis_npz is not None
        else None
    )

    return {
        "line_name": str(obs_npz["line_name"]),
        "role": str(obs_npz["role"]),
        "modality": str(obs_npz["modality"]),
        "n_eval": int(d_obs.shape[0]),
        "g_eval": int(d_obs.shape[1]),
        "h_j1": h_j1,
        "h_j2": h_j2,
        "h_j3": {
            "h_j3a_L1_ens56": h_j3a,
            "h_j3b_L2_ntc_recon_match": h_j3b,
            "l2_selection_diagnostics": l2_sel,
        },
        "secondary_comparators": comp_res,
        "secondary_dorth": {
            "r_dorth_56_point": dorth_r_obs.tolist(),
            "dorth_L1_ens56": dorth_l1,
        },
        "j5a_g4": j5a_res,
        "_holm_inputs": {
            "delta_cf_boot": s0.get("delta_cf_boot", np.zeros(0)),
            "delta_cf_est": float(s0.get("delta_cf", float("nan"))),
            "valid_all": bool(h_j1["valid_all"]),
            "p_var": float(h_j1["p_var"]),
            "tau_batch": float(h_j1["tau_batch"]),
            "has_detectable_signal": bool(h_j1["has_detectable_signal"]),
            "l1_vs_0_est": d_l1_vs_0_est,
            "l1_vs_0_boot": d_l1_vs_0_boot,
            "l2_vs_0_est": d_l2_vs_0_est,
            "l2_vs_0_boot": d_l2_vs_0_boot,
            "l2_selected_idx": b_ntc,
        },
    }


def recompute_cohort_from_exported_arrays_r17(
    arrays_dir: Path,
    models: Sequence[str] = ("STATE_SE", "STATE_HVG"),
    n_boot: int = 4000,
    n_boot_bonf56: int = 20000,
    n_perm: int = 10000,
    n_splits_diag: int = 20,
    seed: int = DEFAULT_R17_SEED,
) -> Dict[str, Any]:
    """
    Recompute the entire Round 17 summary (`H-J1`, `H-J2`, `H-J3`, `J5a`, `J5b`, separate-model
    & pooled Holm-Bonferroni FWER, and `k/K_conf` general claim rules) strictly from exported
    `arrays/<line>/*.npz` with `allow_pickle=False`.
    """
    adir = Path(arrays_dir)
    if not adir.is_dir():
        raise FileNotFoundError(f"Arrays directory not found: {adir}")

    def _load_npz(p: Path) -> Dict[str, Any]:
        with np.load(p, allow_pickle=False) as z:
            out: Dict[str, Any] = {}
            for k in z.files:
                val = z[k]
                out[k] = val.item() if val.ndim == 0 else val
            return out

    line_dirs = {d.name: d for d in sorted(adir.iterdir()) if d.is_dir() and (d / "observed.npz").is_file()}
    ordered_lines = [ln for ln in R17_ALL_EVAL_LINES if ln in line_dirs] + [
        ln for ln in sorted(line_dirs.keys()) if ln not in R17_ALL_EVAL_LINES
    ]

    per_model: Dict[str, Dict[str, Any]] = {m: {} for m in models}
    n_eval_by_line: Dict[str, int] = {}
    scores_56_by_model: Dict[str, Dict[str, np.ndarray]] = {m: {} for m in models}
    sens_g3: Dict[str, Any] = {}

    for ln in ordered_lines:
        ld = line_dirs[ln]
        obs_npz = _load_npz(ld / "observed.npz")
        n_eval_by_line[ln] = int(obs_npz["delta_obs"].shape[0])

        for mkey in models:
            sw_path = ld / f"preds_{mkey}_sweep56.npz"
            if not sw_path.is_file():
                continue
            sw_npz = _load_npz(sw_path)
            g4_path = ld / f"preds_{mkey}_g4_mismatched.npz"
            g4_npz = _load_npz(g4_path) if g4_path.is_file() else None
            ev = evaluate_line_model_from_arrays(
                obs_npz=obs_npz,
                sweep56_npz=sw_npz,
                g4_mis_npz=g4_npz,
                n_boot=n_boot,
                n_boot_bonf56=n_boot_bonf56,
                n_splits_diag=n_splits_diag,
                seed=seed,
                alpha_base=0.05 / 3.0,
            )
            per_model[mkey][ln] = ev
            scores_56_by_model[mkey][ln] = np.asarray(ev["h_j1"]["r_obs_56"], dtype=np.float64)

        sens_path = ld / "SENS_FCT_full_r9_6_incl_vcc2025_h1.npz"
        if sens_path.is_file():
            snz = _load_npz(sens_path)
            obs_s = dict(obs_npz)
            obs_s["fct_delta_pred"] = snz["fct_delta_pred"]
            obs_s["delta_transfer"] = snz["delta_transfer"]
            obs_s["source_grand_mean"] = snz["source_grand_mean"]
            for mkey in models:
                sw_path = ld / f"preds_{mkey}_sweep56.npz"
                if sw_path.is_file():
                    sw_npz = _load_npz(sw_path)
                    sens_ev = evaluate_line_model_from_arrays(
                        obs_npz=obs_s, sweep56_npz=sw_npz, g4_mis_npz=None,
                        n_boot=n_boot, n_boot_bonf56=n_boot_bonf56, n_splits_diag=n_splits_diag, seed=seed,
                    )
                    sens_ev.pop("_internal", None)
                    sens_ev.pop("_holm_inputs", None)
                    sens_g3[mkey] = sens_ev

    # Build Holm-Bonferroni FWER across confirmatory lines
    def _build_holm_entries(mkey: str) -> List[Dict[str, Any]]:
        entries: List[Dict[str, Any]] = []
        for ln in R17_CONFIRMATORY_LINES:
            ev = per_model.get(mkey, {}).get(ln)
            if ev is None:
                continue
            if int(ev["n_eval"]) < 50:
                continue
            hi = ev["_holm_inputs"]

            def _hj1_ci(alpha_s: float, hi_ref=hi) -> bool:
                if not hi_ref["valid_all"]:
                    return False
                c_res = classify_h_j1_under_alpha(
                    hi_ref["valid_all"], hi_ref["p_var"], hi_ref["delta_cf_est"],
                    hi_ref["delta_cf_boot"], hi_ref["tau_batch"], hi_ref["has_detectable_signal"],
                    alpha_adj=alpha_s,
                )
                return bool(np.isfinite(c_res["ci_adj_lo"]) and c_res["ci_adj_lo"] > 0.0)

            def _rule_ci(alpha_s: float, boot_arr=None) -> bool:
                lo, hi_c, _ = _percentile_ci(boot_arr, alpha=alpha_s)
                vf = float(np.mean(np.isfinite(boot_arr))) if boot_arr.size > 0 else 0.0
                return bool(vf >= MIN_BOOT_VALID_FRAC and np.isfinite(lo) and np.isfinite(hi_c) and (lo > 0.0 or hi_c < 0.0))

            entries.append({
                "key": f"{mkey}:{ln}:H-J1",
                "model_key": mkey,
                "line": ln,
                "test": "H-J1",
                "p_value": float(ev["h_j1"]["p_hj1_iu"]),
                "ci_check_fn": _hj1_ci,
            })
            entries.append({
                "key": f"{mkey}:{ln}:H-J3a",
                "model_key": mkey,
                "line": ln,
                "test": "H-J3a",
                "p_value": float(ev["h_j3"]["h_j3a_L1_ens56"]["p_value"]),
                "ci_check_fn": lambda a_s, ba=hi["l1_vs_0_boot"]: _rule_ci(a_s, ba),
            })
            entries.append({
                "key": f"{mkey}:{ln}:H-J3b",
                "model_key": mkey,
                "line": ln,
                "test": "H-J3b",
                "p_value": float(ev["h_j3"]["h_j3b_L2_ntc_recon_match"]["p_value"]),
                "ci_check_fn": (
                    (lambda a_s: False)
                    if hi["l2_selected_idx"] == 0
                    else (lambda a_s, ba=hi["l2_vs_0_boot"]: _rule_ci(a_s, ba))
                ),
            })
        return entries

    holm_by_model: Dict[str, Any] = {}
    pooled_entries: List[Dict[str, Any]] = []
    for mkey in models:
        m_entries = _build_holm_entries(mkey)
        pooled_entries.extend(m_entries)
        holm_res = apply_holm_fwer_r17(m_entries, alpha_family=0.05)
        holm_by_model[mkey] = holm_res

        # Attach Holm-adjusted calls cleanly to confirmatory lines
        for ln in R17_CONFIRMATORY_LINES:
            ev = per_model.get(mkey, {}).get(ln)
            if ev is None:
                continue
            hi = ev["_holm_inputs"]
            h1_step = holm_res["by_key"].get(f"{mkey}:{ln}:H-J1")
            if h1_step is not None:
                c_holm = classify_h_j1_under_alpha(
                    hi["valid_all"], hi["p_var"], hi["delta_cf_est"], hi["delta_cf_boot"],
                    hi["tau_batch"], hi["has_detectable_signal"], alpha_adj=float(h1_step["alpha_eff"]),
                )
                # Ensure Holm classification stays synchronized with p-value step-down & CI gate
                if not h1_step["rejected"] and c_holm["classification"] in (
                    "BATCH_DEPENDENCE_PRACTICALLY_MEANINGFUL",
                    "BATCH_DEPENDENCE_DETECTED_SUB_DELTA",
                ):
                    c_holm["classification"] = "INCONCLUSIVE" if hi["has_detectable_signal"] else "INCONCLUSIVE_LOW_SIGNAL"
                    c_holm["subreason"] = "INSUFFICIENT_PRECISION" if hi["has_detectable_signal"] else None
                ev["h_j1_holm"] = {**c_holm, **h1_step}

    holm_pooled = apply_holm_fwer_r17(pooled_entries, alpha_family=0.05)

    # Compute k / K_conf claim summaries per model
    claim_rollups: Dict[str, Any] = {}
    for mkey in models:
        conf_evaluated = [ln for ln in R17_CONFIRMATORY_LINES if ln in per_model.get(mkey, {})]
        k_conf = len(conf_evaluated)
        k_meaningful = sum(
            1 for ln in conf_evaluated
            if per_model[mkey][ln].get("h_j1_holm", per_model[mkey][ln]["h_j1"])["classification"]
            == "BATCH_DEPENDENCE_PRACTICALLY_MEANINGFUL"
        )
        # Clean out internal bootstrap tensors before JSON serialization
        for ln, ev in per_model[mkey].items():
            ev.pop("_holm_inputs", None)
        claim_rollups[mkey] = {
            "confirmatory_lines": conf_evaluated,
            "K_conf": int(k_conf),
            "k_batch_dependence_practically_meaningful": int(k_meaningful),
            "general_cross_context_claim_permitted": bool(k_conf >= 2 and k_meaningful >= 2),
        }

    j5b = compute_j5b_rank_correlations(
        line_scores_by_model=scores_56_by_model,
        n_eval_by_line=n_eval_by_line,
        n_perm=n_perm,
        seed=seed,
    )

    return {
        "n_boot": int(n_boot),
        "n_boot_bonf56": int(n_boot_bonf56),
        "n_perm": int(n_perm),
        "n_splits_diag": int(n_splits_diag),
        "seed": int(seed),
        "models": per_model,
        "holm_by_model": holm_by_model,
        "holm_pooled_sensitivity": holm_pooled,
        "claim_rollups": claim_rollups,
        "j5b_rank_correlations": j5b,
        "sens_gse295214_h1_full_r9_6": sens_g3,
    }
