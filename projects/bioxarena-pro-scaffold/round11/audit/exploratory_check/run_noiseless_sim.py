#!/usr/bin/env python3
"""
Exploratory Check 4: Noiseless Real-Distribution Simulation for Theorem 1.

Uses:
  - REAL target unperturbed baseline expression mu_{t,g} (`ctrl_mean_libnorm`)
  - REAL R9 source baseline expression mu_{s,g} (`ctrl_mean_libnorm`)
  - REAL perturbation-level shared multiplicative fold-changes delta_{p,g} = log(FC_bar_{p,g})
with ZERO sampling noise and ZERO biological rewiring:
    Delta_{c,p,g}^{sim} = log(1 + mu_{c,g} * exp(delta_{p,g})) - log(1 + mu_{c,g})

Evaluates BOTH:
  1. Unfiltered real `delta_{p,g} = log(FC_bar_{p,g})` AND source-expressed filtered
     `delta_{p,g}` (clipping extreme low-mu source division artifacts only in sensitivity),
     reporting simulated `E1_ConstMed` under BOTH `all_nontarget` and `resp_top50`
     (which in zero-noise equals `resp_top50_split` since Half A == Half B), plus
     Theorem 1's analytical `a_{M, delta^2} = E_{M,delta^2}[w_t] - beta_1 E_{M,delta^2}[w_S]`
     and baseline means `E[w_t]`, `E[w_S]`, `Corr(w_t, w_S)`.
  2. Direct comparison against ACTUAL real observed `E1(all_nontarget)` and
     `E1(resp_top50_split)` across all 12 lines (highlighting `cd4t_stim48` real
     `-0.0203` vs the toy `+0.227` positive-drift claim).
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

EXPLORATORY_CHECK_DIR = Path(__file__).resolve().parent
AUDIT_DIR = EXPLORATORY_CHECK_DIR.parent
for p in (str(EXPLORATORY_CHECK_DIR), str(AUDIT_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from exploratory_core import (  # noqa: E402
    DEFAULT_AUDIT_SEED,
    LINES_12,
    build_r9_lolo_cohort,
    load_12_lines,
)
from independent_r11_metrics import (  # noqa: E402
    EPS,
    _moments,
    _moments_stats,
    evaluate_joint_residuals,
    evaluate_sa8_resp_top50_split,
    orthogonalize_to_ref,
)
from independent_r11_predictors import (  # noqa: E402
    build_predictions,
    compute_source_fc_bar,
    fct_from_fc_bar,
)


def _eval_noiseless_line(
    fc_mat: np.ndarray,
    mu_t: np.ndarray,
    source_names: List[str],
    source_ctrl_mu: Dict[str, np.ndarray],
    source_deltas: Dict[str, Dict[str, np.ndarray]],
    eval_targets: List[str],
    target_cols: np.ndarray,
    mask_finite: np.ndarray,
    delta_clip_abs: Optional[float] = None,
) -> Dict[str, Any]:
    """
    Simulate zero-noise, zero-rewiring response using real (mu_t, {mu_s}, delta_p)
    and compute E1 (`all_nontarget`, `resp_top50`), plus analytical Theorem 1 diagnostics.
    """
    n_t, n_g = fc_mat.shape
    mu_t_safe = np.maximum(np.asarray(mu_t, dtype=np.float64), 0.0)
    c_med = float(np.median(mu_t_safe[np.isfinite(mu_t_safe)]))
    mu_s_safe = {s: np.maximum(np.asarray(source_ctrl_mu[s], dtype=np.float64), 0.0) for s in source_names}

    # Precompute noiseless Delta matrices using exact availability map per target p
    fc_use = fc_mat.copy()
    if delta_clip_abs is not None:
        log_fc = np.clip(np.log(np.maximum(fc_use, 1e-300)), -delta_clip_abs, delta_clip_abs)
        fc_use = np.exp(log_fc)

    delta_sim_t = fct_from_fc_bar(fc_use, mu_t_safe)
    pred_e1_sim = fct_from_fc_bar(fc_use, np.full(n_g, c_med, dtype=np.float64))

    # Build source TMean per pert over the exact source subset avail(p)
    sim_by_src = {s: fct_from_fc_bar(fc_use, mu_s_safe[s]) for s in source_names}
    tmean_sim = np.zeros((n_t, n_g), dtype=np.float64)
    w_sbar_rows = np.zeros((n_t, n_g), dtype=np.float64)
    w_t = mu_t_safe / (1.0 + mu_t_safe)
    w_by_s = {s: mu_s_safe[s] / (1.0 + mu_s_safe[s]) for s in source_names}

    for i, tgt in enumerate(eval_targets):
        avail = [s for s in source_names if tgt in source_deltas[s]]
        tmean_sim[i] = np.mean([sim_by_src[s][i] for s in avail], axis=0)
        w_sbar_rows[i] = np.mean([w_by_s[s] for s in avail], axis=0)

    def _pooled_e1_and_thm1(mode: str) -> Dict[str, Any]:
        xr_o, yr_o = [], []
        num_a, den_a, num_b, den_b = 0.0, 0.0, 0.0, 0.0
        sum_d2_wt, sum_d2_ws, sum_d2 = 0.0, 0.0, 0.0
        for p in range(n_t):
            valid = mask_finite[p].copy()
            tc = int(target_cols[p])
            if 0 <= tc < n_g:
                valid[tc] = False
            idx = np.flatnonzero(valid)
            if idx.size == 0:
                continue
            if mode == "resp_top50" and idx.size > 50:
                idx = idx[np.argsort(np.abs(delta_sim_t[p, idx]))[-50:]]

            ref = tmean_sim[p, idx]
            ro = delta_sim_t[p, idx] - ref
            rp = pred_e1_sim[p, idx] - ref
            xr_o.append(orthogonalize_to_ref(ro, ref))
            yr_o.append(orthogonalize_to_ref(rp, ref))

            # Theorem 1 first-order diagnostics
            d = np.log(np.maximum(fc_use[p, idx], 1e-300))
            d2 = d * d
            wt_p = w_t[idx]
            ws_p = w_sbar_rows[p, idx]
            num_b += float(np.sum(d2 * wt_p * ws_p))
            den_b += float(np.sum(d2 * ws_p * ws_p))
            sum_d2_wt += float(np.sum(d2 * wt_p))
            sum_d2_ws += float(np.sum(d2 * ws_p))
            sum_d2 += float(np.sum(d2))

        r_orth = float(_moments_stats(*_moments(xr_o, yr_o))[0])
        beta1 = num_b / max(den_b, EPS)
        ewt_d2 = sum_d2_wt / max(sum_d2, EPS)
        ews_d2 = sum_d2_ws / max(sum_d2, EPS)
        a_thm1 = ewt_d2 - beta1 * ews_d2
        return {
            "e1_r_orth": r_orth,
            "beta_1_delta2": beta1,
            "E_delta2_wt": ewt_d2,
            "E_delta2_ws": ews_d2,
            "a_thm1_delta2": a_thm1,
        }

    # Global mu/w statistics on valid non-NaN genes
    g_ok = np.isfinite(mu_t)
    w_s_mean_global = np.mean([w_by_s[s] for s in source_names], axis=0)
    wt_g = w_t[g_ok]
    ws_g = w_s_mean_global[g_ok]
    corr_w = float(np.corrcoef(wt_g, ws_g)[0, 1]) if wt_g.size > 2 else float("nan")
    mu_s_mean_global = np.mean([mu_s_safe[s] for s in source_names], axis=0)

    return {
        "c_med": c_med,
        "mean_mu_t": float(np.mean(mu_t_safe[g_ok])),
        "median_mu_t": float(np.median(mu_t_safe[g_ok])),
        "mean_mu_src": float(np.mean(mu_s_mean_global[g_ok])),
        "median_mu_src": float(np.median(mu_s_mean_global[g_ok])),
        "mean_w_t": float(np.mean(wt_g)),
        "mean_w_src": float(np.mean(ws_g)),
        "corr_w_t_w_src": corr_w,
        "all_nontarget": _pooled_e1_and_thm1("all_nontarget"),
        "resp_top50_noiseless_split": _pooled_e1_and_thm1("resp_top50"),
    }


def run_check4(
    seed: int = DEFAULT_AUDIT_SEED,
    out_dir: Path = EXPLORATORY_CHECK_DIR / "outputs",
) -> Dict[str, Any]:
    t0 = time.time()
    out_dir.mkdir(parents=True, exist_ok=True)
    lines = load_12_lines()

    rows: List[Dict[str, Any]] = []
    per_line: Dict[str, Any] = {}

    for ln, grp in LINES_12:
        cohort = build_r9_lolo_cohort(lines, ln)
        fc_bar, t_mean, _ = compute_source_fc_bar(
            eval_targets=cohort.eval_targets,
            source_names=cohort.source_names,
            source_deltas_by_line=cohort.source_deltas_by_line,
            source_ctrl_mean_libnorm=cohort.source_ctrl_mean_libnorm,
        )
        preds, _, _, _ = build_predictions(cohort, include_sa3_anchors=False, include_exploratory_e=False)

        # Real observed E1 (all_nontarget and resp_top50_split) for direct sign/magnitude audit
        real_e1_all = evaluate_joint_residuals(
            delta_obs=cohort.delta_obs,
            pred=preds["e1_const_med"],
            t_mean=t_mean,
            target_cols=cohort.target_cols,
            mode="all_nontarget",
            n_boot=200,
            seed=seed,
            half_a_obs=cohort.half_a_obs,
            half_b_obs=cohort.half_b_obs,
            half_valid=cohort.half_valid,
        )
        real_e1_split = evaluate_sa8_resp_top50_split(
            delta_obs=cohort.delta_obs,
            pred=preds["e1_const_med"],
            t_mean=t_mean,
            target_cols=cohort.target_cols,
            half_a_obs=cohort.half_a_obs,
            half_b_obs=cohort.half_b_obs,
            half_valid=cohort.half_valid,
            top_k=50,
            n_boot=200,
            seed=seed,
        )

        # Mask matching real finite delta_obs
        mask_fin = np.isfinite(cohort.delta_obs) & np.isfinite(t_mean) & np.isfinite(fc_bar)
        sim_unclipped = _eval_noiseless_line(
            fc_mat=fc_bar,
            mu_t=cohort.tgt_ctrl_mean_libnorm,
            source_names=cohort.source_names,
            source_ctrl_mu=cohort.source_ctrl_mean_libnorm,
            source_deltas=cohort.source_deltas_by_line,
            eval_targets=cohort.eval_targets,
            target_cols=cohort.target_cols,
            mask_finite=mask_fin,
            delta_clip_abs=None,
        )
        # Sensitivity check: realistic biological |logFC| <= 2.0 (excluding zero-mu extreme ratios)
        sim_clipped2 = _eval_noiseless_line(
            fc_mat=fc_bar,
            mu_t=cohort.tgt_ctrl_mean_libnorm,
            source_names=cohort.source_names,
            source_ctrl_mu=cohort.source_ctrl_mean_libnorm,
            source_deltas=cohort.source_deltas_by_line,
            eval_targets=cohort.eval_targets,
            target_cols=cohort.target_cols,
            mask_finite=mask_fin,
            delta_clip_abs=2.0,
        )

        e1_obs_all = float(real_e1_all["r_orth_cal_est"])
        e1_obs_split = float(real_e1_split["r_orth_cal_split_est"])
        e1_sim_all = float(sim_unclipped["all_nontarget"]["e1_r_orth"])
        e1_sim_split = float(sim_unclipped["resp_top50_noiseless_split"]["e1_r_orth"])
        e1_sim2_all = float(sim_clipped2["all_nontarget"]["e1_r_orth"])
        e1_sim2_split = float(sim_clipped2["resp_top50_noiseless_split"]["e1_r_orth"])

        per_line[ln] = {
            "group": grp,
            "n_eval": len(cohort.eval_targets),
            "real_observed": {
                "e1_all_cal": e1_obs_all,
                "e1_all_raw": float(real_e1_all["r_orth_raw_est"]),
                "e1_split_cal": e1_obs_split,
                "e1_split_raw": float(real_e1_split["r_orth_raw_split_est"]),
            },
            "noiseless_unclipped_real_logfc": sim_unclipped,
            "noiseless_clipped2_real_logfc": sim_clipped2,
            "sign_match_all_unclipped": bool(np.sign(e1_obs_all) == np.sign(e1_sim_all)),
            "sign_match_all_clipped2": bool(np.sign(e1_obs_all) == np.sign(e1_sim2_all)),
            "sign_match_split_unclipped": bool(np.sign(e1_obs_split) == np.sign(e1_sim_split)),
        }

        row = {
            "line": ln,
            "group": grp,
            "n_eval": len(cohort.eval_targets),
            "mean_wt": sim_unclipped["mean_w_t"],
            "mean_ws": sim_unclipped["mean_w_src"],
            "corr_wt_ws": sim_unclipped["corr_w_t_w_src"],
            "real_e1_all_cal": e1_obs_all,
            "sim_e1_all_unclip": e1_sim_all,
            "sim_e1_all_clip2": e1_sim2_all,
            "thm1_a_all": sim_unclipped["all_nontarget"]["a_thm1_delta2"],
            "sign_match_all": bool(np.sign(e1_obs_all) == np.sign(e1_sim_all)),
            "real_e1_split_cal": e1_obs_split,
            "sim_e1_split_unclip": e1_sim_split,
            "sim_e1_split_clip2": e1_sim2_split,
            "thm1_a_split": sim_unclipped["resp_top50_noiseless_split"]["a_thm1_delta2"],
            "sign_match_split": bool(np.sign(e1_obs_split) == np.sign(e1_sim_split)),
        }
        rows.append(row)
        print(
            f"[Check4 Noiseless] {ln:24s} | wt={row['mean_wt']:.3f} vs ws={row['mean_ws']:.3f} "
            f"| ALL real={row['real_e1_all_cal']:+.4f} vs sim={row['sim_e1_all_unclip']:+.4f} "
            f"(clip2={row['sim_e1_all_clip2']:+.4f}, a={row['thm1_a_all']:+.4f}) match={row['sign_match_all']} "
            f"| SPLIT real={row['real_e1_split_cal']:+.4f} vs sim={row['sim_e1_split_unclip']:+.4f} "
            f"(clip2={row['sim_e1_split_clip2']:+.4f}, a={row['thm1_a_split']:+.4f}) match={row['sign_match_split']}"
        )

    payload = {
        "check": "Noiseless_Theorem1_RealDist_Audit",
        "elapsed_sec": round(time.time() - t0, 2),
        "n_sign_match_all_unclipped": int(sum(1 for r in rows if r["sign_match_all"])),
        "n_sign_match_all_clipped2": int(sum(1 for r in rows if per_line[r["line"]]["sign_match_all_clipped2"])),
        "n_sign_match_split_unclipped": int(sum(1 for r in rows if r["sign_match_split"])),
        "per_line": per_line,
    }
    with open(out_dir / "noiseless_sim.json", "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    with open(out_dir / "noiseless_sim.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    return payload


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=DEFAULT_AUDIT_SEED)
    parser.add_argument("--out-dir", type=Path, default=EXPLORATORY_CHECK_DIR / "outputs")
    args = parser.parse_args()
    run_check4(seed=args.seed, out_dir=args.out_dir)
