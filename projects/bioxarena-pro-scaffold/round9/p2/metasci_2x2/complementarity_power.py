#!/usr/bin/env python3
"""
Fast Monte Carlo power calculation for matched-FPR 2-call complementarity
(Jaccard difference: J_within - J_cross at equal marginal TPR/FPR=0.10,
and 2-call union/mean TPR gain at equal total FPR=0.10).
"""

from __future__ import annotations

import json
import math
from pathlib import Path
import numpy as np
from scipy.stats import norm

OUT_DIR = Path(__file__).resolve().parent / "outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def run_comp_sim() -> None:
    rng = np.random.default_rng(20261122)
    configs = [
        # (M_tasks, D_defects_per_task, rho_within_vendor, rho_cross_vendor)
        (60, 2, 0.65, 0.50),
        (60, 2, 0.85, 0.45),
        (80, 2, 0.65, 0.50),
        (80, 2, 0.65, 0.45),
        (80, 2, 0.85, 0.45),
        (80, 2, 0.85, 0.35),
        (120, 2, 0.65, 0.50),
        (120, 2, 0.85, 0.45),
        (120, 2, 0.85, 0.35),
    ]
    n_mc = 1000
    base_tpr = 0.35
    thr = norm.ppf(1.0 - base_tpr)

    rows = []
    for m_tasks, d_per, rho_w, rho_c in configs:
        # Simulate per-task correlated bivariate normal latent detection across
        # (same_run1, same_run2, cross_run1) with true equal marginal detection rate base_tpr
        # Cov matrix among (S1, S2, X1):
        # corr(S1, S2) = rho_w (within-vendor repeatability), corr(S1, X1) = corr(S2, X1) = rho_c
        cov = np.array([
            [1.0, rho_w, rho_c],
            [rho_w, 1.0, rho_c],
            [rho_c, rho_c, 1.0],
        ])
        L = np.linalg.cholesky(cov)

        d_jacc_list = []
        p_signflip_list = []
        union_gain_list = []

        n_items = m_tasks * d_per
        for _ in range(n_mc):
            # Task cluster effect included implicitly via total corr; decompose task vs item
            z = rng.standard_normal((n_items, 3)) @ L.T
            det = (z > thr).astype(float)  # shape (n_items, 3): S1, S2, X1
            s1, s2, x1 = det[:, 0], det[:, 1], det[:, 2]

            # Global Jaccard difference
            inter_w = np.sum(s1 * s2)
            union_w = np.sum((s1 + s2) > 0)
            j_w = inter_w / max(union_w, 1.0)

            inter_c = np.sum(s1 * x1)
            union_c = np.sum((s1 + x1) > 0)
            j_c = inter_c / max(union_c, 1.0)
            d_jacc_list.append(j_w - j_c)

            # 2-call union gain: (S1 OR X1) - (S1 OR S2) per task cluster
            u_mix = ((s1 + x1) > 0).reshape(m_tasks, d_per).mean(axis=1)
            u_same = ((s1 + s2) > 0).reshape(m_tasks, d_per).mean(axis=1)
            diff_m = u_mix - u_same
            union_gain_list.append(float(np.mean(diff_m)))

            # Normal approximation / paired t-test at task cluster level
            se = float(np.std(diff_m, ddof=1) / math.sqrt(m_tasks))
            t_stat = float(np.mean(diff_m)) / max(se, 1e-12)
            p_signflip_list.append(abs(t_stat) > 1.96)

        rows.append({
            "n_tasks": m_tasks,
            "n_defects": n_items,
            "rho_within": rho_w,
            "rho_cross": rho_c,
            "mean_delta_jaccard_within_minus_cross": round(float(np.mean(d_jacc_list)), 4),
            "sd_delta_jaccard": round(float(np.std(d_jacc_list, ddof=1)), 4),
            "power_delta_jaccard_norm": round(float(np.mean(np.abs(d_jacc_list) > 1.96 * np.std(d_jacc_list, ddof=1) * 0 + (np.array(d_jacc_list) / np.std(d_jacc_list, ddof=1) > 1.96) if False else 0.0)), 4),
            "mean_2call_union_gain": round(float(np.mean(union_gain_list)), 4),
            "sd_2call_union_gain": round(float(np.std(union_gain_list, ddof=1)), 4),
            "power_2call_union_gain": round(float(np.mean(p_signflip_list)), 4),
        })
        # Fix clean analytic/MC normal test for delta_jaccard using per-replicate task bootstrap SE or signal/SD
        z_j = float(np.mean(d_jacc_list)) / float(np.std(d_jacc_list, ddof=1))
        rows[-1]["snr_delta_jaccard"] = round(z_j, 3)
        # Power of two-sided z-test with this empirical SNR:
        rows[-1]["power_delta_jaccard_approx"] = round(float(1.0 - norm.cdf(1.96 - z_j) + norm.cdf(-1.96 - z_j)), 4)
        del rows[-1]["power_delta_jaccard_norm"]
        print(rows[-1], flush=True)

    with open(OUT_DIR / "complementarity_power.json", "w", encoding="utf-8") as f:
        json.dump(rows, f, indent=2)


if __name__ == "__main__":
    run_comp_sim()
