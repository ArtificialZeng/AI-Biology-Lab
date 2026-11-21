#!/usr/bin/env python3
"""
Power, sample-size, bootstrap-coverage, Type-I-error, complementarity-power,
and API-call/cost simulation for the symmetric 2x2 (Generator Vendor x Auditor Vendor)
matched-FPR audit experiment (PREREG_2x2.md).

Uses O(N log N) fast ROC computation and Threshold-Re-Estimated Task-Cluster Bootstrap.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from scipy.stats import norm


OUT_DIR = Path(__file__).resolve().parent / "outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def neyman_pearson_matched_tpr(
    clean_scores: np.ndarray,
    defect_scores: np.ndarray,
    alpha: float,
) -> Tuple[float, np.ndarray]:
    """
    Exact deterministic Neyman-Pearson matched-FPR TPR and per-defect flags z in [0, 1].
    Solves for (tau, theta) on clean_scores such that
      mean(clean > tau) + theta * mean(clean == tau) == alpha
    exactly, then evaluates defect_scores with the same (tau, theta).
    """
    n0 = len(clean_scores)
    target_flags = alpha * n0
    sorted_clean = np.sort(clean_scores)[::-1]
    k_floor = int(math.floor(target_flags))

    if k_floor <= 0:
        tau = float(sorted_clean[0])
    elif k_floor >= n0:
        tau = float(sorted_clean[-1]) - 1e-12
    elif abs(target_flags - k_floor) < 1e-12:
        tau = float(sorted_clean[k_floor - 1])
    else:
        tau = float(sorted_clean[k_floor])

    n_gt = int(np.sum(clean_scores > tau))
    n_eq = int(np.sum(clean_scores == tau))
    theta = (target_flags - n_gt) / n_eq if n_eq > 0 else 0.0
    theta = float(min(1.0, max(0.0, theta)))

    z_defect = (defect_scores > tau) + theta * (defect_scores == tau)
    return float(np.mean(z_defect)), z_defect


def roc_metrics_fast(
    clean_scores: np.ndarray,
    defect_scores: np.ndarray,
    alpha: float = 0.10,
    max_fpr: float = 0.20,
) -> Tuple[float, float, float]:
    """
    Compute (TPR@alpha_FPR, pAUROC[0, max_fpr], full_AUROC) in one sorted pass.
    """
    c_asc = np.sort(clean_scores)
    d_asc = np.sort(defect_scores)
    n0 = len(c_asc)
    n1 = len(d_asc)

    # 1. Exact Neyman-Pearson TPR@alpha
    target_flags = alpha * n0
    k_floor = int(math.floor(target_flags))
    if k_floor <= 0:
        tau = c_asc[-1]
    elif k_floor >= n0:
        tau = c_asc[0] - 1e-12
    elif abs(target_flags - k_floor) < 1e-12:
        tau = c_asc[n0 - k_floor]
    else:
        tau = c_asc[n0 - 1 - k_floor]

    idx_hi_c = int(np.searchsorted(c_asc, tau, side="right"))
    idx_lo_c = int(np.searchsorted(c_asc, tau, side="left"))
    n_gt_c = n0 - idx_hi_c
    n_eq_c = idx_hi_c - idx_lo_c
    theta = (target_flags - n_gt_c) / n_eq_c if n_eq_c > 0 else 0.0
    if theta < 0.0:
        theta = 0.0
    elif theta > 1.0:
        theta = 1.0

    idx_hi_d = int(np.searchsorted(d_asc, tau, side="right"))
    idx_lo_d = int(np.searchsorted(d_asc, tau, side="left"))
    tpr_at_alpha = ((n1 - idx_hi_d) + theta * (idx_hi_d - idx_lo_d)) / n1

    # 2. Exact ROC vertices via adjacent difference mask (avoids np.unique sort overhead)
    mask = np.empty(n0, dtype=bool)
    mask[:-1] = c_asc[:-1] != c_asc[1:]
    mask[-1] = True
    uniq_desc = c_asc[mask][::-1]

    c_gt = 1.0 - np.searchsorted(c_asc, uniq_desc, side="right") / n0
    c_ge = 1.0 - np.searchsorted(c_asc, uniq_desc, side="left") / n0
    d_gt = 1.0 - np.searchsorted(d_asc, uniq_desc, side="right") / n1
    d_ge = 1.0 - np.searchsorted(d_asc, uniq_desc, side="left") / n1

    m = len(uniq_desc)
    fprs = np.empty(2 * m + 2, dtype=float)
    tprs = np.empty(2 * m + 2, dtype=float)
    fprs[0], tprs[0] = 0.0, d_gt[0]
    fprs[1:-1:2] = c_gt
    tprs[1:-1:2] = d_gt
    fprs[2:-1:2] = c_ge
    tprs[2:-1:2] = d_ge
    fprs[-1], tprs[-1] = 1.0, 1.0

    full_auroc = float(np.sum(0.5 * (tprs[1:] + tprs[:-1]) * (fprs[1:] - fprs[:-1])))

    idx = int(np.searchsorted(fprs, max_fpr, side="right"))
    xs = fprs[:idx]
    ys = tprs[:idx]
    if xs[-1] < max_fpr:
        x0, x1 = fprs[idx - 1], fprs[idx]
        y0, y1 = tprs[idx - 1], tprs[idx]
        y_interp = y0 + (y1 - y0) * ((max_fpr - x0) / max(x1 - x0, 1e-12))
        xs = np.append(xs, max_fpr)
        ys = np.append(ys, y_interp)
    pauc = float(np.sum(0.5 * (ys[1:] + ys[:-1]) * (xs[1:] - xs[:-1])) / max_fpr)

    return float(tpr_at_alpha), pauc, full_auroc


def simulate_2x2_study(
    rng: np.random.Generator,
    n_tasks: int,
    n_clean_per_task: int,
    n_defects_per_task: int,
    base_tpr_at_10: float,
    delta_cross_tpr: float,
    auditor_asym_logit: float,
    gen_asym_logit: float,
    icc_task: float,
    corr_auditor_pair: float,
    score_levels: int = 101,
) -> Dict[str, np.ndarray]:
    z_thr = norm.ppf(0.90)
    mu_same = z_thr - norm.ppf(1.0 - base_tpr_at_10)
    mu_cross = z_thr - norm.ppf(1.0 - min(0.99, max(0.01, base_tpr_at_10 + delta_cross_tpr)))

    var_task = icc_task
    var_shared_art = (1.0 - var_task) * corr_auditor_pair
    var_indep = max(1.0 - var_task - var_shared_art, 1e-6)

    sd_task = math.sqrt(var_task)
    sd_art = math.sqrt(var_shared_art)
    sd_ind = math.sqrt(var_indep)

    task_u = rng.normal(0.0, sd_task, size=(n_tasks,))
    clean_scores = np.zeros((2, 2, n_tasks, n_clean_per_task), dtype=float)
    defect_scores = np.zeros((2, 2, n_tasks, n_defects_per_task), dtype=float)
    step = 100.0 / ((score_levels - 1) * 4)

    for g in (0, 1):
        g_shift = gen_asym_logit * (0.5 if g == 1 else -0.5)
        art_clean = task_u[:, None] + rng.normal(0.0, sd_art, size=(n_tasks, n_clean_per_task))
        art_def = task_u[:, None] + rng.normal(0.0, sd_art, size=(n_tasks, n_defects_per_task))
        for a in (0, 1):
            a_shift = auditor_asym_logit * (0.5 if a == 1 else -0.5)
            sep = mu_cross if (g != a) else mu_same
            lat_c = g_shift + a_shift + art_clean + rng.normal(0.0, sd_ind, size=(n_tasks, n_clean_per_task))
            lat_d = g_shift + a_shift + sep + art_def + rng.normal(0.0, sd_ind, size=(n_tasks, n_defects_per_task))
            clean_scores[g, a] = np.round((100.0 * norm.cdf(0.75 * lat_c)) / step) * step
            defect_scores[g, a] = np.round((100.0 * norm.cdf(0.75 * lat_d)) / step) * step

    return {"clean": clean_scores, "defect": defect_scores}


def evaluate_interaction(
    clean_scores: np.ndarray,
    defect_scores: np.ndarray,
    alpha: float = 0.10,
    max_fpr_pauc: float = 0.20,
    n_boot: int = 200,
    rng: np.random.Generator | None = None,
) -> Dict[str, float]:
    n_tasks = clean_scores.shape[2]

    t00, p00, a00 = roc_metrics_fast(clean_scores[0, 0].ravel(), defect_scores[0, 0].ravel(), alpha, max_fpr_pauc)
    t01, p01, a01 = roc_metrics_fast(clean_scores[0, 1].ravel(), defect_scores[0, 1].ravel(), alpha, max_fpr_pauc)
    t10, p10, a10 = roc_metrics_fast(clean_scores[1, 0].ravel(), defect_scores[1, 0].ravel(), alpha, max_fpr_pauc)
    t11, p11, a11 = roc_metrics_fast(clean_scores[1, 1].ravel(), defect_scores[1, 1].ravel(), alpha, max_fpr_pauc)

    delta_tpr = 0.5 * ((t01 - t00) + (t10 - t11))
    delta_pauc = 0.5 * ((p01 - p00) + (p10 - p11))
    delta_auroc = 0.5 * ((a01 - a00) + (a10 - a11))

    assert rng is not None
    boot_tpr = np.empty(n_boot)
    boot_pauc = np.empty(n_boot)
    boot_auroc = np.empty(n_boot)
    boot_idx_mat = rng.integers(0, n_tasks, size=(n_boot, n_tasks))

    c00, d00 = clean_scores[0, 0], defect_scores[0, 0]
    c01, d01 = clean_scores[0, 1], defect_scores[0, 1]
    c10, d10 = clean_scores[1, 0], defect_scores[1, 0]
    c11, d11 = clean_scores[1, 1], defect_scores[1, 1]

    for b in range(n_boot):
        idx = boot_idx_mat[b]
        bt00, bp00, ba00 = roc_metrics_fast(c00[idx].ravel(), d00[idx].ravel(), alpha, max_fpr_pauc)
        bt01, bp01, ba01 = roc_metrics_fast(c01[idx].ravel(), d01[idx].ravel(), alpha, max_fpr_pauc)
        bt10, bp10, ba10 = roc_metrics_fast(c10[idx].ravel(), d10[idx].ravel(), alpha, max_fpr_pauc)
        bt11, bp11, ba11 = roc_metrics_fast(c11[idx].ravel(), d11[idx].ravel(), alpha, max_fpr_pauc)
        boot_tpr[b] = 0.5 * ((bt01 - bt00) + (bt10 - bt11))
        boot_pauc[b] = 0.5 * ((bp01 - bp00) + (bp10 - bp11))
        boot_auroc[b] = 0.5 * ((ba01 - ba00) + (ba10 - ba11))

    ci_t_lo, ci_t_hi = np.percentile(boot_tpr, [2.5, 97.5])
    ci_p_lo, ci_p_hi = np.percentile(boot_pauc, [2.5, 97.5])
    ci_a_lo, ci_a_hi = np.percentile(boot_auroc, [2.5, 97.5])

    return {
        "delta_tpr": float(delta_tpr),
        "delta_pauc": float(delta_pauc),
        "delta_auroc": float(delta_auroc),
        "ci_t_lo": float(ci_t_lo),
        "ci_t_hi": float(ci_t_hi),
        "ci_p_lo": float(ci_p_lo),
        "ci_p_hi": float(ci_p_hi),
        "ci_a_lo": float(ci_a_lo),
        "ci_a_hi": float(ci_a_hi),
        "sig_tpr": float(ci_t_lo > 0.0 or ci_t_hi < 0.0),
        "sig_pauc": float(ci_p_lo > 0.0 or ci_p_hi < 0.0),
        "sig_auroc": float(ci_a_lo > 0.0 or ci_a_hi < 0.0),
    }


def run_power_grid() -> List[Dict[str, float]]:
    rng = np.random.default_rng(20261121)
    configs = [
        (40, 1, 2, 0.00),
        (40, 1, 2, 0.10),
        (40, 1, 2, 0.15),
        (60, 1, 2, 0.00),
        (60, 1, 2, 0.08),
        (60, 1, 2, 0.10),
        (60, 1, 2, 0.15),
        (80, 2, 2, 0.00),
        (80, 2, 2, 0.08),
        (80, 2, 2, 0.10),
        (80, 2, 2, 0.12),
        (80, 2, 2, 0.15),
        (120, 2, 2, 0.00),
        (120, 2, 2, 0.08),
        (120, 2, 2, 0.10),
        (120, 2, 2, 0.12),
        (160, 2, 2, 0.00),
        (160, 2, 2, 0.08),
        (160, 2, 2, 0.10),
    ]

    n_mc = 250
    n_boot = 200
    results = []

    pop_cache: Dict[float, Tuple[float, float]] = {}
    for _, _, _, d_val in configs:
        if d_val not in pop_cache:
            sim_pop = simulate_2x2_study(
                rng=rng,
                n_tasks=8000,
                n_clean_per_task=1,
                n_defects_per_task=2,
                base_tpr_at_10=0.35,
                delta_cross_tpr=d_val,
                auditor_asym_logit=0.80,
                gen_asym_logit=0.30,
                icc_task=0.20,
                corr_auditor_pair=0.35,
            )
            _, p00, a00 = roc_metrics_fast(sim_pop["clean"][0, 0].ravel(), sim_pop["defect"][0, 0].ravel(), 0.10, 0.20)
            _, p01, a01 = roc_metrics_fast(sim_pop["clean"][0, 1].ravel(), sim_pop["defect"][0, 1].ravel(), 0.10, 0.20)
            _, p10, a10 = roc_metrics_fast(sim_pop["clean"][1, 0].ravel(), sim_pop["defect"][1, 0].ravel(), 0.10, 0.20)
            _, p11, a11 = roc_metrics_fast(sim_pop["clean"][1, 1].ravel(), sim_pop["defect"][1, 1].ravel(), 0.10, 0.20)
            pop_cache[d_val] = (
                0.5 * ((p01 - p00) + (p10 - p11)),
                0.5 * ((a01 - a00) + (a10 - a11)),
            )

    for n_tasks, n_clean, n_def, delta in configs:
        recs = []
        for _ in range(n_mc):
            sim = simulate_2x2_study(
                rng=rng,
                n_tasks=n_tasks,
                n_clean_per_task=n_clean,
                n_defects_per_task=n_def,
                base_tpr_at_10=0.35,
                delta_cross_tpr=delta,
                auditor_asym_logit=0.80,
                gen_asym_logit=0.30,
                icc_task=0.20,
                corr_auditor_pair=0.35,
            )
            res = evaluate_interaction(
                sim["clean"],
                sim["defect"],
                alpha=0.10,
                max_fpr_pauc=0.20,
                n_boot=n_boot,
                rng=rng,
            )
            recs.append(res)

        true_dpauc, true_dauroc = pop_cache[delta]
        mean_dtpr = float(np.mean([r["delta_tpr"] for r in recs]))
        mean_dpauc = float(np.mean([r["delta_pauc"] for r in recs]))
        mean_dauroc = float(np.mean([r["delta_auroc"] for r in recs]))
        pw_tpr = float(np.mean([r["sig_tpr"] for r in recs]))
        pw_pauc = float(np.mean([r["sig_pauc"] for r in recs]))
        pw_auroc = float(np.mean([r["sig_auroc"] for r in recs]))
        cov_tpr = float(np.mean([r["ci_t_lo"] <= delta <= r["ci_t_hi"] for r in recs]))
        cov_pauc = float(np.mean([r["ci_p_lo"] <= true_dpauc <= r["ci_p_hi"] for r in recs]))
        cov_auroc = float(np.mean([r["ci_a_lo"] <= true_dauroc <= r["ci_a_hi"] for r in recs]))
        sd_tpr = float(np.std([r["delta_tpr"] for r in recs], ddof=1))
        sd_pauc = float(np.std([r["delta_pauc"] for r in recs], ddof=1))
        sd_auroc = float(np.std([r["delta_auroc"] for r in recs], ddof=1))

        row = {
            "n_tasks": n_tasks,
            "n_clean_per_task": n_clean,
            "n_defects_per_task": n_def,
            "total_items_per_gen": n_tasks * (n_clean + n_def),
            "true_delta_tpr10": delta,
            "est_mean_delta_tpr10": round(mean_dtpr, 4),
            "est_sd_delta_tpr10": round(sd_tpr, 4),
            "power_or_t1e_tpr10": round(pw_tpr, 4),
            "ci95_cov_tpr10": round(cov_tpr, 4),
            "true_delta_pauc20": round(true_dpauc, 4),
            "est_mean_delta_pauc20": round(mean_dpauc, 4),
            "est_sd_delta_pauc20": round(sd_pauc, 4),
            "power_or_t1e_pauc20": round(pw_pauc, 4),
            "ci95_cov_pauc20": round(cov_pauc, 4),
            "true_delta_auroc": round(true_dauroc, 4),
            "est_mean_delta_auroc": round(mean_dauroc, 4),
            "est_sd_delta_auroc": round(sd_auroc, 4),
            "power_or_t1e_auroc": round(pw_auroc, 4),
            "ci95_cov_auroc": round(cov_auroc, 4),
            "n_mc": n_mc,
        }
        results.append(row)
        print(
            f"M={n_tasks:3d} C={n_clean} D={n_def} dTPR={delta:.2f} | "
            f"TPR10 pw={pw_tpr:.3f} cov={cov_tpr:.3f} sd={sd_tpr:.3f} | "
            f"pAUC20 pw={pw_pauc:.3f} cov={cov_pauc:.3f} sd={sd_pauc:.3f} | "
            f"AUROC pw={pw_auroc:.3f} cov={cov_auroc:.3f} sd={sd_auroc:.3f}",
            flush=True,
        )
    return results


def compute_cost_table() -> List[Dict[str, object]]:
    scenarios = [
        {"name": "2x2_Minimal_M40_C1_D2_K3", "V": 2, "M_target": 40, "M_gen_pool": 55, "C": 1, "D": 2, "K": 3, "pilot_M": 8},
        {"name": "2x2_Standard_M60_C1_D2_K4", "V": 2, "M_target": 60, "M_gen_pool": 80, "C": 1, "D": 2, "K": 4, "pilot_M": 10},
        {"name": "2x2_Recommended_M80_C2_D2_K4", "V": 2, "M_target": 80, "M_gen_pool": 105, "C": 2, "D": 2, "K": 4, "pilot_M": 10},
        {"name": "2x2_HighPower_M120_C2_D2_K4", "V": 2, "M_target": 120, "M_gen_pool": 155, "C": 2, "D": 2, "K": 4, "pilot_M": 12},
        {"name": "2x2_UltraPower_M160_C2_D2_K3", "V": 2, "M_target": 160, "M_gen_pool": 205, "C": 2, "D": 2, "K": 3, "pilot_M": 12},
        {"name": "3x3_Extension_M80_C2_D2_K3", "V": 3, "M_target": 80, "M_gen_pool": 115, "C": 2, "D": 2, "K": 3, "pilot_M": 10},
        {"name": "4x4_AllLocalVendors_M80_C2_D2_K3", "V": 4, "M_target": 80, "M_gen_pool": 130, "C": 2, "D": 2, "K": 3, "pilot_M": 10},
    ]

    gen_in_k, gen_out_k = 2.5, 1.5
    aud_in_k, aud_out_k = 5.0, 1.0
    rep_in_k, rep_out_k = 6.0, 1.5

    rows = []
    for s in scenarios:
        V = s["V"]
        M = s["M_target"]
        M_pool = s["M_gen_pool"]
        C = s["C"]
        D = s["D"]
        K = s["K"]
        pilot_M = s["pilot_M"]

        pilot_gen_calls = int(round(V * pilot_M * 1.5))
        pilot_aud_calls = V * V * pilot_M * 2 * 2
        pilot_calls = pilot_gen_calls + pilot_aud_calls

        gen_calls = int(round(V * M_pool * C * 1.5))
        items_per_gen = M * (C + D)
        total_artifacts = V * items_per_gen
        audit_calls = V * V * items_per_gen * K
        expected_repair_calls = int(round(V * V * M * (0.10 * C + 0.40 * D)))
        total_llm_calls = pilot_calls + gen_calls + audit_calls + expected_repair_calls

        tot_gen = pilot_gen_calls + gen_calls
        tot_aud = pilot_aud_calls + audit_calls
        tot_rep = expected_repair_calls
        in_mtok = (tot_gen * gen_in_k + tot_aud * aud_in_k + tot_rep * rep_in_k) / 1000.0
        out_mtok = (tot_gen * gen_out_k + tot_aud * aud_out_k + tot_rep * rep_out_k) / 1000.0

        usd_mid = in_mtok * 3.0 + out_mtok * 15.0
        usd_high = in_mtok * 5.0 + out_mtok * 25.0
        usd_low = in_mtok * 0.5 + out_mtok * 2.5

        rows.append({
            **s,
            "total_artifacts": total_artifacts,
            "pilot_calls": pilot_calls,
            "phase1_gen_calls": gen_calls,
            "phase2_audit_calls": audit_calls,
            "phase3_repair_calls": expected_repair_calls,
            "total_llm_calls": total_llm_calls,
            "est_input_mtok": round(in_mtok, 2),
            "est_output_mtok": round(out_mtok, 2),
            "est_usd_tier_low_0.5_2.5": round(usd_low, 2),
            "est_usd_tier_mid_3_15": round(usd_mid, 2),
            "est_usd_tier_high_5_25": round(usd_high, 2),
        })
    return rows


def main() -> None:
    print("Running streamlined O(N log N) 2x2 matched-FPR power and coverage simulation...", flush=True)
    power_rows = run_power_grid()
    cost_rows = compute_cost_table()

    with open(OUT_DIR / "power_simulation.json", "w", encoding="utf-8") as f:
        json.dump(power_rows, f, indent=2)
    with open(OUT_DIR / "cost_estimates.json", "w", encoding="utf-8") as f:
        json.dump(cost_rows, f, indent=2)

    md_lines = [
        "# Power, Coverage, and Call-Count Simulation Results\n",
        "Generated by `power_and_cost_sim.py` (seed `20261121`, `n_mc=250`, `n_boot=200`).\n",
        "Assumptions: baseline TPR@10%FPR = 0.35, auditor asymmetry logit = 0.80 (~5% vs ~25% uncalibrated FPR), "
        "generator asymmetry logit = 0.30, task ICC = 0.20, cross-auditor correlation = 0.35.\n\n",
        "## 1. Power, Type-I Error, and 95% Cluster-Bootstrap Coverage\n",
        "| M (tasks) | C | D | True ΔTPR@10% | Mean ΔTPR | SD(TPR10) | Pw/T1E TPR10 | Cov TPR10 | True ΔpAUC20 | Mean ΔpAUC20 | SD(pAUC20) | Pw/T1E pAUC20 | Cov pAUC20 | True ΔAUROC | SD(AUROC) | Pw/T1E AUROC | Cov AUROC |\n",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n",
    ]
    for r in power_rows:
        md_lines.append(
            f"| {r['n_tasks']} | {r['n_clean_per_task']} | {r['n_defects_per_task']} | "
            f"{r['true_delta_tpr10']:.2f} | {r['est_mean_delta_tpr10']:.4f} | {r['est_sd_delta_tpr10']:.4f} | "
            f"{r['power_or_t1e_tpr10']:.4f} | {r['ci95_cov_tpr10']:.4f} | "
            f"{r['true_delta_pauc20']:.4f} | {r['est_mean_delta_pauc20']:.4f} | {r['est_sd_delta_pauc20']:.4f} | "
            f"{r['power_or_t1e_pauc20']:.4f} | {r['ci95_cov_pauc20']:.4f} | "
            f"{r['true_delta_auroc']:.4f} | {r['est_sd_delta_auroc']:.4f} | {r['power_or_t1e_auroc']:.4f} | {r['ci95_cov_auroc']:.4f} |\n"
        )

    md_lines.extend([
        "\n## 2. Call-Count and Indicative API Cost Estimates\n",
        "| Scenario | V | M_target | C | D | K | Total Artifacts | Pilot Calls | Gen Calls | Audit Calls | Repair Calls | Total Calls | In MTok | Out MTok | Low ($0.5/$2.5) | Mid ($3/$15) | High ($5/$25) |\n",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n",
    ])
    for c in cost_rows:
        md_lines.append(
            f"| {c['name']} | {c['V']} | {c['M_target']} | {c['C']} | {c['D']} | {c['K']} | "
            f"{c['total_artifacts']} | {c['pilot_calls']} | {c['phase1_gen_calls']} | {c['phase2_audit_calls']} | "
            f"{c['phase3_repair_calls']} | {c['total_llm_calls']} | {c['est_input_mtok']:.2f} | {c['est_output_mtok']:.2f} | "
            f"${c['est_usd_tier_low_0.5_2.5']:.2f} | ${c['est_usd_tier_mid_3_15']:.2f} | ${c['est_usd_tier_high_5_25']:.2f} |\n"
        )

    with open(OUT_DIR / "simulation_summary.md", "w", encoding="utf-8") as f:
        f.writelines(md_lines)
    print("Saved outputs to", OUT_DIR, flush=True)


if __name__ == "__main__":
    main()
