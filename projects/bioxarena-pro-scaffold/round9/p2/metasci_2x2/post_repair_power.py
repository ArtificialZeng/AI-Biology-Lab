#!/usr/bin/env python3
"""
Fast Monte Carlo power calculation for Endpoint E4 (Post-Repair Deterministic Pass Rate)
and Endpoint E3 (True-Localization Rate at matched 10% FPR) across M in {60, 80, 120, 160}.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
import numpy as np

OUT_DIR = Path(__file__).resolve().parent / "outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def run_repair_power() -> None:
    rng = np.random.default_rng(20261123)
    # In the 2x2 post-repair contrast, each defective item starts at 0% pass rate;
    # after 1 bounded repair round triggered ONLY if blocked at matched 10% FPR,
    # an item passes iff (blocked AND repair_succeeds).
    # Let p_pass_same = tpr_same * fix_rate_same, p_pass_cross = tpr_cross * fix_rate_cross.
    # Evaluate paired task-cluster sign-flip / bootstrap over M tasks with D=2 defects/task.
    configs = [
        # (M, D, p_pass_same, p_pass_cross, corr_across_cells)
        (60, 2, 0.21, 0.21, 0.30),  # Null: 0.35*0.60 vs 0.35*0.60
        (60, 2, 0.21, 0.31, 0.30),  # +10 pp repair pass gain
        (80, 2, 0.21, 0.21, 0.30),  # Null
        (80, 2, 0.21, 0.29, 0.30),  # +8 pp
        (80, 2, 0.21, 0.31, 0.30),  # +10 pp
        (80, 2, 0.21, 0.33, 0.30),  # +12 pp
        (120, 2, 0.21, 0.21, 0.30), # Null
        (120, 2, 0.21, 0.29, 0.30), # +8 pp
        (120, 2, 0.21, 0.31, 0.30), # +10 pp
        (160, 2, 0.21, 0.29, 0.30), # +8 pp
    ]
    n_mc = 1000
    rows = []
    for m_tasks, d_per, p_s, p_c, corr in configs:
        # Task-clustered correlated Bernoulli approximation via Gaussian copula
        from scipy.stats import norm
        thr_s = norm.ppf(1.0 - p_s)
        thr_c = norm.ppf(1.0 - p_c)

        sig_count = 0
        diffs = []
        for _ in range(n_mc):
            # Generate 4 cells per task-defect: (00: same, 01: cross, 10: cross, 11: same)
            # Task random effect (var=0.20) + shared item effect (var=0.25) + cell noise (var=0.55)
            u_task = rng.normal(0.0, math.sqrt(0.20), size=(m_tasks, 1))
            v_item = u_task + rng.normal(0.0, math.sqrt(0.25), size=(m_tasks, d_per))
            sd_e = math.sqrt(0.55)
            y00 = (v_item + rng.normal(0.0, sd_e, size=(m_tasks, d_per)) > thr_s).astype(float)
            y01 = (v_item + rng.normal(0.0, sd_e, size=(m_tasks, d_per)) > thr_c).astype(float)
            y10 = (v_item + rng.normal(0.0, sd_e, size=(m_tasks, d_per)) > thr_c).astype(float)
            y11 = (v_item + rng.normal(0.0, sd_e, size=(m_tasks, d_per)) > thr_s).astype(float)

            # Task-level interaction contrast
            d_m = 0.5 * ((y01.mean(axis=1) - y00.mean(axis=1)) + (y10.mean(axis=1) - y11.mean(axis=1)))
            mean_d = float(np.mean(d_m))
            se_d = float(np.std(d_m, ddof=1) / math.sqrt(m_tasks))
            diffs.append(mean_d)
            if abs(mean_d) > 1.96 * max(se_d, 1e-12):
                sig_count += 1

        rows.append({
            "n_tasks": m_tasks,
            "n_defects_per_gen": m_tasks * d_per,
            "p_repair_same": p_s,
            "p_repair_cross": p_c,
            "true_delta_repair": round(p_c - p_s, 4),
            "est_sd_delta": round(float(np.std(diffs, ddof=1)), 4),
            "power_or_t1e": round(sig_count / n_mc, 4),
        })
        print(rows[-1], flush=True)

    with open(OUT_DIR / "post_repair_power.json", "w", encoding="utf-8") as f:
        json.dump(rows, f, indent=2)


if __name__ == "__main__":
    run_repair_power()
