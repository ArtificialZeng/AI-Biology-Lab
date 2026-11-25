#!/usr/bin/env python3
"""
Pure-Synthetic Numerical Checks & Supporting Illustrations for `PREREG_I v2 §4` (`theory_sim_r16.py`).

Touches ZERO real perturbation datasets.
  1. Property P1 (By Construction, `I2`, `I3`): Verifies numerically that the secondary
     Systema-inspired double-orthogonalized metric `r_dorth_cal` maps `Transfer_Mean`,
     `Transfer_Mean_GlobalScale`, `Source_Grand_Mean`, and 100 random perturbation-varying
     affine combinations in `span(1, delta_transfer_p, source_grand_mean)` identically to
     0.0 (`<= 1e-14`), whereas legacy `r_orth_cal` leaks non-zero correlation on
     `Source_Grand_Mean` when global stress amplitude differs between target and source.
  2. Mechanism Hypothesis M2 (Supporting Synthetic Illustration Only, `I4`): Illustrates
     how finite-source `FC_bar` estimation with low-expression zero-floor genes + cross-line
     baseline divergence yields low `E1_ConstMed(all_nontarget)`, elevated
     `E1_ConstMed(resp_top50_split)` under independent Half A / Half B selection, and
     positive saturation-weight shift `Delta_W > 0`. This is a supporting synthetic
     illustration, NOT an empirical proof or mathematical theorem.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict

import numpy as np

_THIS_DIR = str(Path(__file__).resolve().parent)
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)

from metrics_r16 import (  # noqa: E402
    DEFAULT_R16_SEED,
    compute_joint_cross_mask_r16,
    compute_orth_suite_r16,
    compute_winner_curse_weight_shift_r16,
)
from model_adapters_r16 import build_closed_form_reference_preds  # noqa: E402


def _fct(mu: np.ndarray, fc: np.ndarray) -> np.ndarray:
    return np.log1p(mu * fc) - np.log1p(mu)


def verify_property_p1_synthetic(
    seed: int = DEFAULT_R16_SEED,
    n_random_affine: int = 100,
    n_boot: int = 200,
) -> Dict[str, Any]:
    """Numerical check of Property P1 (null-space invariance by construction) for secondary `r_dorth_cal`."""
    rng = np.random.default_rng(seed)
    G, P, S = 1000, 120, 6
    z_shared = rng.normal(0.0, 1.8, size=G)
    mu_src = np.exp(z_shared[None, :] + rng.normal(0.0, 0.7, size=(S, G)) - 2.0)
    mu_t = np.exp(z_shared + rng.normal(0.0, 0.8, size=G) - 2.0)

    b_stress = rng.normal(0.0, 0.25, size=G)
    delta_true = rng.normal(0.0, 0.35, size=(P, G)) * (rng.random((P, G)) < 0.20) + 0.45 * b_stress[None, :]
    fc_true = np.exp(delta_true)

    src_deltas = np.stack([_fct(mu_src[s][None, :], fc_true) for s in range(S)], axis=0)
    tm = src_deltas.mean(axis=0)
    b_s = tm.mean(axis=0)

    gamma_t = 1.7
    target_true = (
        _fct(mu_t[None, :], fc_true)
        + (gamma_t - 1.0) * 0.45 * b_stress[None, :]
        + rng.normal(0.0, 0.06, size=(P, G))
    )
    ha = target_true + rng.normal(0.0, 0.04, size=(P, G))
    hb = target_true + rng.normal(0.0, 0.04, size=(P, G))
    obs = 0.5 * (ha + hb)
    t_idx = np.arange(P, dtype=np.int64)
    hv = np.ones(P, dtype=bool)

    refs = build_closed_form_reference_preds(mu_t, fc_true, tm, b_s, global_scale_alpha=1.4)
    summary: Dict[str, Any] = {"named": {}}
    for name, pred in refs.items():
        r_orth = compute_orth_suite_r16(
            obs, pred, tm, t_idx, source_grand_mean=None,
            delta_half_a=ha, delta_half_b=hb, half_valid=hv, n_boot=n_boot, seed=seed
        )
        r_dorth = compute_orth_suite_r16(
            obs, pred, tm, t_idx, source_grand_mean=b_s,
            delta_half_a=ha, delta_half_b=hb, half_valid=hv, n_boot=n_boot, seed=seed
        )
        summary["named"][name] = {
            "r_orth_cal": r_orth["estimate"],
            "r_dorth_cal": r_dorth["estimate"],
            "r_dorth_ci975_lo": r_dorth["ci975_lo"],
            "r_dorth_ci975_hi": r_dorth["ci975_hi"],
        }

    max_abs_affine = 0.0
    for k in range(n_random_affine):
        a0 = rng.normal(0.0, 2.0, size=(P, 1))
        a1 = rng.normal(0.0, 2.0, size=(P, 1))
        a2 = rng.normal(0.0, 2.0, size=(P, 1))
        pred_aff = a0 + a1 * tm + a2 * b_s[None, :]
        res = compute_orth_suite_r16(
            obs, pred_aff, tm, t_idx, source_grand_mean=b_s,
            delta_half_a=ha, delta_half_b=hb, half_valid=hv, n_boot=20, seed=seed + k + 1
        )
        max_abs_affine = max(max_abs_affine, abs(float(res["estimate"])))

    summary["n_random_affine_tested"] = n_random_affine
    summary["max_abs_random_affine_r_dorth_cal"] = max_abs_affine
    summary["property_p1_pass"] = bool(
        abs(summary["named"]["Source_Grand_Mean"]["r_orth_cal"]) > 0.15
        and abs(summary["named"]["Source_Grand_Mean"]["r_dorth_cal"]) <= 1e-14
        and abs(summary["named"]["Transfer_Mean_GlobalScale"]["r_dorth_cal"]) <= 1e-14
        and abs(summary["named"]["Transfer_Mean"]["r_dorth_cal"]) <= 1e-14
        and max_abs_affine <= 1e-14
        and summary["named"]["FCT"]["r_dorth_ci975_lo"] > 0.0
    )
    return summary


def illustrate_mechanism_m2_synthetic(
    seed: int = DEFAULT_R16_SEED,
    n_boot: int = 400,
) -> Dict[str, Any]:
    """Supporting synthetic illustration for Mechanism Hypothesis M2 (`PREREG_I v2 §4.2`)."""
    rng = np.random.default_rng(seed + 100)
    G, P, S = 2000, 150, 6

    z_shared = rng.normal(-2.2, 2.0, size=G)
    mu_src = np.exp(z_shared[None, :] + rng.normal(0.0, 0.75, size=(S, G)))
    mu_t = np.exp(z_shared + rng.normal(0.0, 0.85, size=G))
    c_med = float(np.median(mu_t))

    delta_true = rng.normal(0.0, 0.55, size=(P, G)) * (rng.random((P, G)) < 0.12)
    fc_true = np.exp(delta_true)

    src_deltas_obs = np.stack([
        _fct(mu_src[s][None, :], fc_true) + rng.normal(0.0, 0.025, size=(P, G))
        for s in range(S)
    ], axis=0)
    tm = src_deltas_obs.mean(axis=0)
    b_s = tm.mean(axis=0)

    mu_safe = np.maximum(mu_src[:, None, :], 1e-3)
    fc_est_s = np.clip((np.expm1(src_deltas_obs) * (1.0 + mu_safe) + mu_safe) / mu_safe, 0.05, 20.0)
    fc_bar = np.exp(np.mean(np.log(fc_est_s), axis=0))

    target_true = 0.75 * _fct(mu_t[None, :], fc_true) + rng.normal(0.0, 0.03, size=(P, G))
    ha = target_true + rng.normal(0.0, 0.02, size=(P, G))
    hb = target_true + rng.normal(0.0, 0.02, size=(P, G))
    obs = 0.5 * (ha + hb)
    t_idx = np.arange(P, dtype=np.int64)
    hv = np.ones(P, dtype=bool)

    pred_fct = _fct(mu_t[None, :], fc_bar)
    pred_e1 = _fct(np.full_like(mu_t, c_med)[None, :], fc_bar)

    e1_cross_orth = compute_joint_cross_mask_r16(
        obs, pred_e1, tm, t_idx, source_grand_mean=None,
        delta_half_a=ha, delta_half_b=hb, half_valid=hv, top_k=50, n_boot=n_boot, seed=seed
    )
    e1_cross_dorth = compute_joint_cross_mask_r16(
        obs, pred_e1, tm, t_idx, source_grand_mean=b_s,
        delta_half_a=ha, delta_half_b=hb, half_valid=hv, top_k=50, n_boot=n_boot, seed=seed
    )
    fct_orth = compute_orth_suite_r16(
        obs, pred_fct, tm, t_idx, source_grand_mean=None,
        delta_half_a=ha, delta_half_b=hb, half_valid=hv, n_boot=n_boot, seed=seed
    )
    w_shift = compute_winner_curse_weight_shift_r16(
        mu_t, mu_src, ha, hb, tm, t_idx, hv, top_k=50, n_boot=n_boot, seed=seed
    )

    e1_all_orth = float(e1_cross_orth["all_common_summary"]["estimate"])
    e1_split_orth = float(e1_cross_orth["split_summary"]["estimate"])
    return {
        "note": "Supporting synthetic illustration only (not empirical proof)",
        "e1_orth_all": e1_all_orth,
        "e1_orth_split": e1_split_orth,
        "e1_orth_jump": float(e1_cross_orth["estimate"]),
        "e1_orth_jump_ci975_lo": float(e1_cross_orth["ci975_lo"]),
        "e1_dorth_jump": float(e1_cross_dorth["estimate"]),
        "e1_dorth_jump_ci975_lo": float(e1_cross_dorth["ci975_lo"]),
        "fct_orth_all": float(fct_orth["estimate"]),
        "weight_shift_delta_w": float(w_shift["estimate"]),
        "weight_shift_ci975_lo": float(w_shift["ci975_lo"]),
        "w_t_all_mean": float(w_shift["w_t_all_mean"]),
        "w_t_top50_mean": float(w_shift["w_t_top50_mean"]),
        "illustration_consistent": bool(
            e1_all_orth < 0.05
            and e1_split_orth > 0.15
            and float(e1_cross_orth["ci975_lo"]) > 0.0
            and float(e1_cross_dorth["ci975_lo"]) > 0.0
            and float(w_shift["ci975_lo"]) > 0.0
        ),
    }


def main() -> None:
    out = {
        "property_p1_by_construction": verify_property_p1_synthetic(),
        "mechanism_m2_synthetic_illustration": illustrate_mechanism_m2_synthetic(),
    }
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
