#!/usr/bin/env python3
"""
Frozen statistical evaluation and JSON-schema verification module for PREREG_2x2.md.
Implements:
1. Strict schema & Constitution citation validator (CrossAudit I3).
2. Deterministic Neyman-Pearson matched-FPR TPR(alpha), pAUROC[0,0.20], and full AUROC.
3. Threshold-re-estimated task-cluster bootstrap 95% CIs & two-sided bootstrap p-values
   for the 2x2 cross-vs-same interaction contrast, simple effects, and Holm-Bonferroni adjustment.
4. 5-fold out-of-fold (OOF) threshold calibration check.
5. Equal-budget 2-call complementarity & matched-FPR Jaccard / independence reference.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np


RULE_PATTERN = re.compile(r"^CA-(SCI|DATA|STAT|SPEC|META|DCL)-[0-9]{3}$")


def validate_audit_response(
    payload: Dict[str, Any],
    expected_constitution_commit: str,
) -> Tuple[bool, List[str]]:
    """
    Validate an auditor JSON response against audit_response_schema.json and I3.
    Returns (is_valid, error_reasons). Invalid responses map to ESCALATE / audit_failure.
    """
    errors: List[str] = []
    req_top = {"constitution_commit", "verdict", "overall_block_score", "findings"}
    missing = req_top - set(payload.keys())
    if missing:
        errors.append(f"Missing top-level keys: {sorted(missing)}")
        return False, errors

    if set(payload.keys()) - req_top:
        errors.append(f"Unexpected top-level keys: {sorted(set(payload.keys()) - req_top)}")

    if str(payload["constitution_commit"]) != expected_constitution_commit:
        errors.append(
            f"Constitution commit mismatch: got {payload['constitution_commit']!r}, "
            f"expected {expected_constitution_commit!r}"
        )

    if payload["verdict"] not in ("PASS", "BLOCKED"):
        errors.append(f"Invalid verdict: {payload['verdict']!r}")

    score = payload["overall_block_score"]
    if not isinstance(score, (int, float)) or math.isnan(score) or score < 0 or score > 100:
        errors.append(f"Invalid overall_block_score: {score!r}")

    findings = payload["findings"]
    if not isinstance(findings, list):
        errors.append("findings must be a list")
        return False, errors

    req_finding = {
        "rule_id", "severity", "confidence", "defect_family",
        "file_path", "symbol_or_line", "evidence_quote", "mechanism",
    }
    has_blocker = False
    for idx, f in enumerate(findings):
        if not isinstance(f, dict):
            errors.append(f"finding[{idx}] is not an object")
            continue
        f_miss = req_finding - set(f.keys())
        if f_miss:
            errors.append(f"finding[{idx}] missing keys: {sorted(f_miss)}")
            continue
        if not RULE_PATTERN.match(str(f["rule_id"])):
            errors.append(f"finding[{idx}] invalid rule_id: {f['rule_id']!r}")
        if f["severity"] not in ("BLOCKER", "ADVISORY"):
            errors.append(f"finding[{idx}] invalid severity: {f['severity']!r}")
        elif f["severity"] == "BLOCKER":
            has_blocker = True
        conf = f["confidence"]
        if not isinstance(conf, (int, float)) or conf < 0 or conf > 100:
            errors.append(f"finding[{idx}] invalid confidence: {conf!r}")
        if not str(f["evidence_quote"]).strip():
            errors.append(f"finding[{idx}] empty evidence_quote")

    if payload["verdict"] == "BLOCKED" and not has_blocker:
        errors.append("verdict=BLOCKED without any BLOCKER finding (I3/I6 violation)")
    if payload["verdict"] == "PASS" and has_blocker:
        errors.append("verdict=PASS despite BLOCKER finding present")

    return len(errors) == 0, errors


def neyman_pearson_threshold(
    clean_scores: np.ndarray,
    alpha: float,
) -> Tuple[float, float]:
    """
    Solve deterministically on clean_scores for (tau, theta) with theta in [0, 1] such that
      mean(clean > tau) + theta * mean(clean == tau) == alpha.
    """
    n0 = len(clean_scores)
    if n0 == 0:
        raise ValueError("clean_scores cannot be empty")
    target_flags = alpha * n0
    c_desc = np.sort(clean_scores)[::-1]
    k_floor = int(math.floor(target_flags))

    if k_floor <= 0:
        tau = float(c_desc[0])
    elif k_floor >= n0:
        tau = float(c_desc[-1]) - 1e-12
    elif abs(target_flags - k_floor) < 1e-12:
        tau = float(c_desc[k_floor - 1])
    else:
        tau = float(c_desc[k_floor])

    n_gt = int(np.sum(clean_scores > tau))
    n_eq = int(np.sum(clean_scores == tau))
    theta = (target_flags - n_gt) / n_eq if n_eq > 0 else 0.0
    theta = float(min(1.0, max(0.0, theta)))
    return tau, theta


def apply_threshold(scores: np.ndarray, tau: float, theta: float) -> np.ndarray:
    """Return fractional indicator z in [0, 1] for matched operating point (tau, theta)."""
    return (scores > tau).astype(float) + theta * (scores == tau).astype(float)


def roc_summary(
    clean_scores: np.ndarray,
    defect_scores: np.ndarray,
    alpha: float = 0.10,
    max_fpr: float = 0.20,
) -> Dict[str, float]:
    """Compute exact TPR@alpha, normalized pAUROC[0, max_fpr], and full AUROC."""
    tau, theta = neyman_pearson_threshold(clean_scores, alpha)
    tpr_alpha = float(np.mean(apply_threshold(defect_scores, tau, theta)))
    realized_fpr = float(np.mean(apply_threshold(clean_scores, tau, theta)))

    c_asc = np.sort(clean_scores)
    d_asc = np.sort(defect_scores)
    n0, n1 = len(c_asc), len(d_asc)
    mask = np.empty(n0, dtype=bool)
    mask[:-1] = c_asc[:-1] != c_asc[1:]
    mask[-1] = True
    uniq_desc = c_asc[mask][::-1]

    c_gt = 1.0 - np.searchsorted(c_asc, uniq_desc, side="right") / n0
    c_ge = 1.0 - np.searchsorted(c_asc, uniq_desc, side="left") / n0
    d_gt = 1.0 - np.searchsorted(d_asc, uniq_desc, side="right") / n1
    d_ge = 1.0 - np.searchsorted(d_asc, uniq_desc, side="left") / n1

    fprs = np.empty(2 * len(uniq_desc) + 2, dtype=float)
    tprs = np.empty(2 * len(uniq_desc) + 2, dtype=float)
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

    return {
        "tau": tau,
        "theta": theta,
        "realized_fpr": realized_fpr,
        "tpr_at_alpha": tpr_alpha,
        "pauroc": pauc,
        "auroc": full_auroc,
    }


def holm_adjust(p_values: Sequence[float]) -> List[float]:
    """Holm-Bonferroni step-down adjusted p-values preserving original order."""
    m = len(p_values)
    order = sorted(range(m), key=lambda i: p_values[i])
    adj = [0.0] * m
    running_max = 0.0
    for rank, idx in enumerate(order):
        val = min(1.0, (m - rank) * float(p_values[idx]))
        running_max = max(running_max, val)
        adj[idx] = running_max
    return adj


def oof_threshold_calibration_check(
    clean_scores_by_task: np.ndarray,
    defect_scores_by_task: np.ndarray,
    alpha: float = 0.10,
    n_folds: int = 5,
    seed: int = 20261121,
) -> Dict[str, float]:
    """
    5-fold task-cluster cross-validated threshold calibration check for a single (G, A) cell.
    Fits (tau, theta) on 4 training folds of clean items, applies out-of-fold to held-out
    clean and defect items, and reports realized OOF FPR & OOF TPR.
    """
    n_tasks = clean_scores_by_task.shape[0]
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n_tasks)
    folds = np.array_split(perm, n_folds)

    oof_clean_flags: List[np.ndarray] = []
    oof_defect_flags: List[np.ndarray] = []
    for k in range(n_folds):
        test_idx = folds[k]
        train_idx = np.concatenate([folds[j] for j in range(n_folds) if j != k])
        tau_tr, theta_tr = neyman_pearson_threshold(clean_scores_by_task[train_idx].ravel(), alpha)
        oof_clean_flags.append(apply_threshold(clean_scores_by_task[test_idx].ravel(), tau_tr, theta_tr))
        oof_defect_flags.append(apply_threshold(defect_scores_by_task[test_idx].ravel(), tau_tr, theta_tr))

    return {
        "target_fpr": alpha,
        "oof_realized_fpr": float(np.mean(np.concatenate(oof_clean_flags))),
        "oof_tpr": float(np.mean(np.concatenate(oof_defect_flags))),
    }


def cluster_bootstrap_2x2(
    clean_scores: np.ndarray,
    defect_scores: np.ndarray,
    alpha: float = 0.10,
    max_fpr_pauc: float = 0.20,
    n_boot: int = 4000,
    seed: int = 20261121,
) -> Dict[str, Any]:
    """
    Threshold-Re-Estimated Task-Cluster Bootstrap over M tasks for all 2x2 contrasts.
    clean_scores shape: (2, 2, M, C)
    defect_scores shape: (2, 2, M, D)
    """
    n_tasks = clean_scores.shape[2]
    rng = np.random.default_rng(seed)

    cell_obs = {}
    for g in (0, 1):
        for a in (0, 1):
            cell_obs[(g, a)] = roc_summary(clean_scores[g, a].ravel(), defect_scores[g, a].ravel(), alpha, max_fpr_pauc)

    def contrasts(m_dict: Dict[Tuple[int, int], Dict[str, float]], key: str) -> Tuple[float, float, float]:
        sg0 = m_dict[(0, 1)][key] - m_dict[(0, 0)][key]
        sg1 = m_dict[(1, 0)][key] - m_dict[(1, 1)][key]
        return 0.5 * (sg0 + sg1), sg0, sg1

    metrics = ["tpr_at_alpha", "pauroc", "auroc"]
    obs = {k: contrasts(cell_obs, k) for k in metrics}
    boots = {k: np.empty((n_boot, 3), dtype=float) for k in metrics}

    idx_mat = rng.integers(0, n_tasks, size=(n_boot, n_tasks))
    for b in range(n_boot):
        idx = idx_mat[b]
        cell_b = {}
        for g in (0, 1):
            for a in (0, 1):
                cell_b[(g, a)] = roc_summary(
                    clean_scores[g, a, idx].ravel(),
                    defect_scores[g, a, idx].ravel(),
                    alpha,
                    max_fpr_pauc,
                )
        for k in metrics:
            boots[k][b] = contrasts(cell_b, k)

    out: Dict[str, Any] = {"cells": {f"{g}_{a}": cell_obs[(g, a)] for g in (0, 1) for a in (0, 1)}, "contrasts": {}}
    for k in metrics:
        for col, label in enumerate(["delta_interaction", "simple_gen0", "simple_gen1"]):
            arr = boots[k][:, col]
            lo, hi = np.percentile(arr, [2.5, 97.5])
            # Two-sided bootstrap p-value (continuity-adjusted)
            p_boot = 2.0 * min(
                (1.0 + np.sum(arr <= 0.0)) / (n_boot + 1.0),
                (1.0 + np.sum(arr >= 0.0)) / (n_boot + 1.0),
            )
            out["contrasts"][f"{k}__{label}"] = {
                "estimate": float(obs[k][col]),
                "ci95": [float(lo), float(hi)],
                "boot_p": float(min(1.0, p_boot)),
            }
    return out


def matched_fpr_two_call_complementarity(
    clean_a1: np.ndarray,
    defect_a1: np.ndarray,
    clean_a2: np.ndarray,
    defect_a2: np.ndarray,
    clean_b1: np.ndarray,
    defect_b1: np.ndarray,
    alpha: float = 0.10,
) -> Dict[str, float]:
    """
    Compare 2-call same-vendor (A1 + A2) vs 2-call cross-vendor (A1 + B1)
    at EQUAL total FPR = alpha (via mean-score fusion), AND compute defect-level
    Jaccard overlap when EACH individual call is thresholded at marginal FPR = alpha.
    """
    # 1. Equal-marginal-FPR Jaccard (fixes Song 2026 §5.3 marginal-rate confound)
    tau_a1, th_a1 = neyman_pearson_threshold(clean_a1, alpha)
    tau_a2, th_a2 = neyman_pearson_threshold(clean_a2, alpha)
    tau_b1, th_b1 = neyman_pearson_threshold(clean_b1, alpha)

    z_a1 = apply_threshold(defect_a1, tau_a1, th_a1)
    z_a2 = apply_threshold(defect_a2, tau_a2, th_a2)
    z_b1 = apply_threshold(defect_b1, tau_b1, th_b1)

    # Fractional soft Jaccard for tie-consistency: sum(min(u,v)) / sum(max(u,v))
    j_within = float(np.sum(np.minimum(z_a1, z_a2)) / max(np.sum(np.maximum(z_a1, z_a2)), 1e-12))
    j_cross = float(np.sum(np.minimum(z_a1, z_b1)) / max(np.sum(np.maximum(z_a1, z_b1)), 1e-12))

    r_a1 = float(np.mean(z_a1))
    r_b1 = float(np.mean(z_b1))
    indep_ref = r_a1 + r_b1 - r_a1 * r_b1
    obs_union_cross = float(np.mean(np.maximum(z_a1, z_b1)))

    # 2. Equal-total-FPR (alpha) 2-call mean-score fusion
    same_pair_clean = 0.5 * (clean_a1 + clean_a2)
    same_pair_def = 0.5 * (defect_a1 + defect_a2)
    cross_pair_clean = 0.5 * (clean_a1 + clean_b1)
    cross_pair_def = 0.5 * (defect_a1 + defect_b1)

    tau_sp, th_sp = neyman_pearson_threshold(same_pair_clean, alpha)
    tau_cp, th_cp = neyman_pearson_threshold(cross_pair_clean, alpha)
    tpr_same_2call = float(np.mean(apply_threshold(same_pair_def, tau_sp, th_sp)))
    tpr_cross_2call = float(np.mean(apply_threshold(cross_pair_def, tau_cp, th_cp)))

    return {
        "jaccard_within_at_alpha": j_within,
        "jaccard_cross_at_alpha": j_cross,
        "delta_jaccard_within_minus_cross": j_within - j_cross,
        "obs_union_cross_at_marginal_alpha": obs_union_cross,
        "independence_ref_union": indep_ref,
        "tpr_2call_same_mean_at_total_alpha": tpr_same_2call,
        "tpr_2call_cross_mean_at_total_alpha": tpr_cross_2call,
        "delta_2call_cross_minus_same_at_total_alpha": tpr_cross_2call - tpr_same_2call,
    }


def _self_test() -> None:
    # 1. Schema validation test
    ok, errs = validate_audit_response({
        "constitution_commit": "abc123",
        "verdict": "BLOCKED",
        "overall_block_score": 85.0,
        "findings": [{
            "rule_id": "CA-SCI-001",
            "severity": "BLOCKER",
            "confidence": 90,
            "defect_family": "ALGO_MATH_LOGIC",
            "file_path": "solver.py",
            "symbol_or_line": "L42",
            "evidence_quote": "norm = np.sum(np.abs(v))",
            "mechanism": "Spec requires L2 norm, code uses L1 norm.",
        }],
    }, "abc123")
    assert ok and not errs, errs

    ok_bad, errs_bad = validate_audit_response({
        "constitution_commit": "abc123",
        "verdict": "BLOCKED",
        "overall_block_score": 85.0,
        "findings": [],
    }, "abc123")
    assert not ok_bad and "BLOCKED without any BLOCKER" in errs_bad[0]

    # 2. Matched-FPR & 2x2 bootstrap smoke test
    rng = np.random.default_rng(42)
    clean = rng.normal(20, 10, size=(2, 2, 30, 2)).clip(0, 100)
    defect = rng.normal(55, 15, size=(2, 2, 30, 2)).clip(0, 100)
    res = cluster_bootstrap_2x2(clean, defect, n_boot=100)
    assert abs(res["cells"]["0_0"]["realized_fpr"] - 0.10) < 1e-9
    oof = oof_threshold_calibration_check(clean[0, 0], defect[0, 0])
    assert 0.0 <= oof["oof_realized_fpr"] <= 0.30
    assert holm_adjust([0.01, 0.04, 0.03]) == [0.03, 0.06, 0.06]
    comp = matched_fpr_two_call_complementarity(
        clean[0, 0].ravel(), defect[0, 0].ravel(),
        clean[0, 0].ravel(), defect[0, 0].ravel(),
        clean[0, 1].ravel(), defect[0, 1].ravel(),
        alpha=0.10,
    )
    assert abs(comp["jaccard_within_at_alpha"] - 1.0) < 1e-9
    print("matched_fpr_eval.py self-test PASSED!")


if __name__ == "__main__":
    _self_test()
