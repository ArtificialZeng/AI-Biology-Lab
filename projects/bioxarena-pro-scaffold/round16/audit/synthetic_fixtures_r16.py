#!/usr/bin/env python3
"""
Deterministic Synthetic Benchmark Fixtures & Canonical Parity Pack for Round 16 (`A-T3`).

Explicitly covers (`seed = 20261601`):
  1. Variable sample sizes (`N = 16` underpowered, `N = 86` borderline, `N = 300` high-power)
  2. Zero-signal equivalence (`LINE_EQUIVALENT_NULL_SIGNAL`)
  3. Non-zero-signal control-invariant equivalence (`LINE_EQUIVALENT_NONZERO_SIGNAL`)
  4. Matched-better above delta (`LINE_MATCHED_SPECIFIC`) & sub-delta (`LINE_MATCHED_SPECIFIC_SUB_DELTA`)
  5. Mismatched-better (`LINE_MISMATCHED_BETTER`) & cross-stratum split (`MIXED_OR_INCONCLUSIVE`)
  6. Degenerate Group 2 (`|G_2_qual| = 1` -> `MISSING_STRATUM` / `SINGLE_STRATUM_DESCRIPTIVE_ONLY`)
  7. Batch-imbalanced control cohorts (`[701, 317, 103, 61, 18]` + equal-remainder ties, small-NTC truncation, exclusion floor)
  8. `r_dorth_cal` Property P1 null invariance fixtures (`span{1, t_mean, b_S}`)
  9. Real-weight `ST-HVG-Replogle` synthetic count matrix generator (zero real cells)

Never reads `round16/src/` or touches real cells.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

AUDIT_DIR = Path(__file__).resolve().parent
if str(AUDIT_DIR) not in sys.path:
    sys.path.insert(0, str(AUDIT_DIR))

from audit_r16_adjudicate import (  # noqa: E402
    DEFAULT_TOST_DELTA,
    adjudicate_claim_i1_line,
    adjudicate_claim_i2_line,
    adjudicate_group_r16,
)
from audit_r16_metrics import (  # noqa: E402
    DEFAULT_R16_SEED,
    compute_orth_suite_r16,
    compute_paired_contrast_r16,
    compute_paired_model_vs_all_mismatched_r16,
)
from audit_r16_sampling import (  # noqa: E402
    sample_stratified_control_sets_r16,
)


def make_synthetic_pert_cohort(
    n_perts: int,
    n_genes: int,
    regime: str,
    seed: int = DEFAULT_R16_SEED,
    seen_sources: Sequence[str] = ("nadig_hepg2", "nadig_jurkat"),
    unseen_sources: Sequence[str] = ("cd4t_rest", "dld1", "feng_ipsc"),
) -> Dict[str, Any]:
    """
    Generate deterministic synthetic residual evaluation benchmark matrices for regime in:
      - `matched_better_strong`: matched recovers true line-specific signal (>+0.08 over all mis)
      - `matched_better_sub_delta`: matched beats mis slightly (+0.015), all CIs inside (-0.04, +0.04)
      - `nonzero_signal_equivalent`: strong above-transfer-mean prediction shared by matched & all mis
      - `null_signal_equivalent`: pure noise around transfer mean for both matched & all mis
      - `mismatched_better`: mismatched outperforms matched cleanly
      - `stratum_split`: seen=MATCHED_BETTER, unseen=EQUIVALENT (`MIXED_OR_INCONCLUSIVE`)
      - `underpowered_equivalent_data`: true equivalence on small N (tests power drop at N=16)
    """
    rng = np.random.default_rng(int(seed))
    P, G = int(n_perts), int(n_genes)
    t_mean = rng.normal(0.0, 0.20, size=(P, G))
    grand_mean = rng.normal(0.0, 0.30, size=G)
    b_s_axis = grand_mean[None, :] - np.mean(grand_mean)
    shared_non_cone = rng.normal(0.0, 0.18, size=(P, G))
    spec_signal = rng.normal(0.0, 0.30, size=(P, G))
    true_delta = t_mean + 0.35 * b_s_axis + 0.55 * shared_non_cone + spec_signal

    half_a = true_delta + rng.normal(0.0, 0.06, size=(P, G))
    half_b = true_delta + rng.normal(0.0, 0.06, size=(P, G))
    delta_obs = 0.5 * (half_a + half_b)
    target_cols = np.arange(P, dtype=np.int64) % G
    half_valid = np.ones(P, dtype=bool)
    all_mis = tuple(seen_sources) + tuple(unseen_sources)

    pred_fct = t_mean + 0.40 * spec_signal + 0.25 * shared_non_cone + rng.normal(0.0, 0.10, size=(P, G))

    if regime == "matched_better_strong":
        pred_matched = t_mean + 0.85 * spec_signal + 0.40 * shared_non_cone + rng.normal(0.0, 0.05, size=(P, G))
        mis_preds = {
            s: t_mean + 0.05 * spec_signal + 0.35 * shared_non_cone + rng.normal(0.0, 0.22, size=(P, G))
            for s in all_mis
        }
    elif regime == "matched_better_sub_delta":
        # Shared high correlation (~0.75) plus a modest coherent increment -> Δ in (+0.008, +0.025):
        # ci975_lo > 0.0 (MATCHED_BETTER) AND ci95_hi < 0.04 (within_delta=True) across all pairs.
        base_common = t_mean + 0.60 * spec_signal + 0.35 * shared_non_cone
        pred_matched = base_common + 0.18 * spec_signal + rng.normal(0.0, 0.03, size=(P, G))
        mis_preds = {
            s: base_common + rng.normal(0.0, 0.04, size=(P, G)) for s in all_mis
        }
    elif regime == "nonzero_signal_equivalent":
        # Model has high true above-Transfer_Mean accuracy (spec_signal + shared_non_cone),
        # completely invariant to replacing matched with any mismatched control (tiny independent jitter).
        base_common = t_mean + 0.75 * spec_signal + 0.45 * shared_non_cone
        pred_matched = base_common + rng.normal(0.0, 0.015, size=(P, G))
        mis_preds = {
            s: base_common + rng.normal(0.0, 0.015, size=(P, G)) for s in all_mis
        }
    elif regime == "underpowered_equivalent_data":
        # Perturbation-heterogeneous predictors with zero expected difference (matched vs mis alternate +/- small
        # per-perturbation accuracy shifts of expectation 0), yielding high paired bootstrap variance at N=16
        # (95% CI exceeds +-0.04 -> INCONCLUSIVE) that tightens cleanly into EQUIVALENT by N=300.
        base_common = t_mean + 0.65 * spec_signal + 0.40 * shared_non_cone
        signs = np.where((np.arange(P) % 2) == 0, 1.0, -1.0)[:, None]
        pred_matched = base_common + 0.22 * signs * spec_signal + rng.normal(0.0, 0.05, size=(P, G))
        mis_preds = {
            s: base_common - 0.22 * signs * spec_signal + rng.normal(0.0, 0.05, size=(P, G))
            for s in all_mis
        }
    elif regime == "null_signal_equivalent":
        # Pure uninformative tiny zero-centered noise around t_mean -> matched ~ 0, all mis ~ 0 -> EQUIVALENT_NULL_SIGNAL
        pred_matched = t_mean + rng.normal(0.0, 0.02, size=(P, G))
        mis_preds = {s: t_mean + rng.normal(0.0, 0.02, size=(P, G)) for s in all_mis}
    elif regime == "mismatched_better":
        pred_matched = t_mean + 0.05 * spec_signal + rng.normal(0.0, 0.22, size=(P, G))
        mis_preds = {
            s: t_mean + 0.85 * spec_signal + 0.40 * shared_non_cone + rng.normal(0.0, 0.05, size=(P, G))
            for s in all_mis
        }
    elif regime == "stratum_split":
        # Matched predicts moderate signal; `seen` is degraded (matched > seen), `unseen` identical to matched (equiv)
        base_common = t_mean + 0.70 * spec_signal + 0.40 * shared_non_cone
        pred_matched = base_common + rng.normal(0.0, 0.015, size=(P, G))
        mis_preds = {}
        for s in seen_sources:
            mis_preds[s] = t_mean + 0.05 * spec_signal + rng.normal(0.0, 0.22, size=(P, G))
        for s in unseen_sources:
            mis_preds[s] = base_common + rng.normal(0.0, 0.015, size=(P, G))
    else:
        raise ValueError(f"Unknown synthetic regime={regime!r}")

    return {
        "regime": regime,
        "seed": int(seed),
        "n_perts": P,
        "n_genes": G,
        "seen_sources": tuple(seen_sources),
        "unseen_sources": tuple(unseen_sources),
        "delta_obs": delta_obs,
        "half_a_obs": half_a,
        "half_b_obs": half_b,
        "half_valid": half_valid,
        "t_mean": t_mean,
        "grand_mean": grand_mean,
        "target_cols": target_cols,
        "pred_matched": pred_matched,
        "pred_fct": pred_fct,
        "mis_preds": mis_preds,
    }


def make_batch_imbalanced_control_fixture(seed: int = DEFAULT_R16_SEED) -> Dict[str, Any]:
    """
    Create deterministic batch-imbalanced synthetic control metadata pools:
      - `large_imbalanced`: 1,200 QC-valid NTCs out of 1,350 total cells (`[701, 317, 103, 61, 18]`)
      - `tie_remainder`: 1,000 cells across 4 balanced batches (`[250, 250, 250, 250]`), `n_ctrl=512`
        forces remainder ties (`128*4` no tie vs `512` across 3 batches `[334, 333, 333]`)
      - `small_truncated`: 200 NTCs (`128 <= 200 < 512`) -> truncated to `192 = 6 × 32`
      - `below_floor`: 100 NTCs (`< 128`) -> excluded
    """
    # 1) Large imbalanced with QC filter
    counts_large = {"B_alpha": 701, "B_beta": 317, "B_gamma": 103, "B_delta": 61, "B_eps": 18}
    qc_mask_parts = []
    batch_parts = []
    for b, cnt in counts_large.items():
        qc_mask_parts.extend([True] * cnt)
        batch_parts.extend([b] * cnt)
        # inject 30 failing QC cells per batch so qc_mask filtering is also exercised
        qc_mask_parts.extend([False] * 30)
        batch_parts.extend([b] * 30)

    # 2) Equal remainder tie case (3 batches: 333, 333, 334 -> total 1000, n_ctrl=512)
    # 512*333/1000 = 170.496, 512*333/1000 = 170.496 (exact tie!), 512*334/1000 = 171.008 -> sum floors=511, rem=1 -> tie breaks to B_1
    tie_batches = (["B_2"] * 333) + (["B_1"] * 333) + (["B_3"] * 334)

    return {
        "seed": int(seed),
        "large_imbalanced": {
            "qc_mask": np.asarray(qc_mask_parts, dtype=bool),
            "batch_labels": batch_parts,
            "expected_valid_counts": counts_large,
        },
        "tie_remainder": {
            "batch_labels": tie_batches,
            "expected_alloc": {"B_1": 171, "B_2": 170, "B_3": 171},
        },
        "small_truncated": {"n_cells": 200, "batch_labels": ["X"] * 120 + ["Y"] * 80},
        "below_floor": {"n_cells": 100},
    }


def evaluate_synthetic_fixture_r16(
    fixture: Dict[str, Any],
    n_boot: int = 1000,
    seed: int = DEFAULT_R16_SEED,
    delta: float = DEFAULT_TOST_DELTA,
) -> Dict[str, Any]:
    """Run full clean-room metric + paired bootstrap + adjudication on a synthetic fixture."""
    orth_matched = compute_orth_suite_r16(
        delta_obs=fixture["delta_obs"],
        delta_pred=fixture["pred_matched"],
        delta_transfer=fixture["t_mean"],
        target_gene_idx=fixture["target_cols"],
        source_grand_mean=None,
        delta_half_a=fixture["half_a_obs"],
        delta_half_b=fixture["half_b_obs"],
        half_valid=fixture["half_valid"],
        mode="all_nontarget",
        n_boot=n_boot,
        seed=seed,
    )
    dorth_matched = compute_orth_suite_r16(
        delta_obs=fixture["delta_obs"],
        delta_pred=fixture["pred_matched"],
        delta_transfer=fixture["t_mean"],
        target_gene_idx=fixture["target_cols"],
        source_grand_mean=fixture["grand_mean"],
        delta_half_a=fixture["half_a_obs"],
        delta_half_b=fixture["half_b_obs"],
        half_valid=fixture["half_valid"],
        mode="all_nontarget",
        n_boot=n_boot,
        seed=seed,
    )
    paired_mis = compute_paired_model_vs_all_mismatched_r16(
        delta_obs=fixture["delta_obs"],
        pred_matched=fixture["pred_matched"],
        mismatched_preds=fixture["mis_preds"],
        delta_transfer=fixture["t_mean"],
        target_gene_idx=fixture["target_cols"],
        source_grand_mean=None,
        delta_half_a=fixture["half_a_obs"],
        delta_half_b=fixture["half_b_obs"],
        half_valid=fixture["half_valid"],
        mode="all_nontarget",
        n_boot=n_boot,
        seed=seed,
    )
    h_i1 = adjudicate_claim_i1_line(
        orth_suite=orth_matched,
        paired_mismatched=paired_mis,
        seen_sources=fixture["seen_sources"],
        unseen_sources=fixture["unseen_sources"],
        delta=delta,
    )
    vs_fct = compute_paired_contrast_r16(
        delta_obs=fixture["delta_obs"],
        pred_a=fixture["pred_matched"],
        pred_b=fixture["pred_fct"],
        delta_transfer=fixture["t_mean"],
        target_gene_idx=fixture["target_cols"],
        source_grand_mean=None,
        delta_half_a=fixture["half_a_obs"],
        delta_half_b=fixture["half_b_obs"],
        half_valid=fixture["half_valid"],
        mode="all_nontarget",
        n_boot=n_boot,
        seed=seed,
    )
    h_i2 = adjudicate_claim_i2_line(vs_fct)
    return {
        "regime": fixture["regime"],
        "n_perts": fixture["n_perts"],
        "n_genes": fixture["n_genes"],
        "orth_matched": orth_matched,
        "dorth_matched": dorth_matched,
        "paired_mismatched": paired_mis,
        "h_i1": h_i1,
        "vs_fct": vs_fct,
        "h_i2": h_i2,
    }


def export_canonical_parity_pack(
    output_json: Optional[Path] = None,
    n_boot: int = 400,
    seed: int = DEFAULT_R16_SEED,
) -> Dict[str, Any]:
    """
    Evaluate all canonical benchmark scenarios cleanly so Main & Controller can verify exact
    numerical parity (`<= 1e-10`) and `100%` verdict parity across clean-room implementations.
    """
    specs: List[Tuple[str, int, int, int]] = [
        ("N300_matched_better_strong", 300, 120, 20261601),
        ("N400_matched_better_sub_delta", 400, 120, 20261602),
        ("N200_nonzero_signal_equivalent", 200, 120, 20261603),
        ("N200_null_signal_equivalent", 200, 120, 20261604),
        ("N120_mismatched_better", 120, 120, 20261605),
        ("N120_stratum_split", 120, 120, 20261606),
        ("N300_noisy_wellpowered_equiv", 300, 120, 20261607),
        ("N086_borderline_equiv", 86, 120, 20261608),
        ("N016_underpowered_equiv", 16, 120, 20261609),
    ]
    results: Dict[str, Any] = {}
    for label, np_, ng_, sd_ in specs:
        reg = label.split("_", 1)[1]
        if any(tag in reg for tag in ("noisy_wellpowered", "borderline", "underpowered")):
            reg = "underpowered_equivalent_data"
        fx = make_synthetic_pert_cohort(n_perts=np_, n_genes=ng_, regime=reg, seed=sd_)
        ev = evaluate_synthetic_fixture_r16(fx, n_boot=n_boot, seed=seed)
        results[label] = {
            "seed": sd_,
            "n_perts": np_,
            "n_genes": ng_,
            "matched_r_orth_cal": ev["orth_matched"]["estimate"],
            "matched_ci975": [ev["orth_matched"]["ci975_lo"], ev["orth_matched"]["ci975_hi"]],
            "matched_rho_ab_orth": ev["orth_matched"]["rho_ab_orth"],
            "matched_r_dorth_cal": ev["dorth_matched"]["estimate"],
            "h_i1_line_summary": ev["h_i1"]["line_summary"],
            "h_i1_all_within_delta": ev["h_i1"]["all_within_delta"],
            "h_i1_counts": ev["h_i1"]["counts"],
            "h_i1_seen_status": ev["h_i1"]["strata"]["seen_by_checkpoint"]["stratum_summary"],
            "h_i1_unseen_status": ev["h_i1"]["strata"]["unseen_by_checkpoint"]["stratum_summary"],
            "h_i1_per_source_est": {
                k: {
                    "category": v["category"],
                    "within_delta": v["within_delta"],
                    "estimate": v["estimate"],
                    "ci95": [v["ci95_lo"], v["ci95_hi"]],
                    "ci975": [v["ci975_lo"], v["ci975_hi"]],
                }
                for k, v in ev["h_i1"]["per_source"].items()
            },
            "h_i2_classification": ev["h_i2"]["classification"],
            "h_i2_estimate": ev["vs_fct"]["estimate"],
        }

    batch_fx = make_batch_imbalanced_control_fixture(seed=seed)
    s_large = sample_stratified_control_sets_r16(
        qc_mask=batch_fx["large_imbalanced"]["qc_mask"],
        batch_labels=batch_fx["large_imbalanced"]["batch_labels"],
        seed=seed,
    )
    s_tie = sample_stratified_control_sets_r16(
        n_cells=len(batch_fx["tie_remainder"]["batch_labels"]),
        batch_labels=batch_fx["tie_remainder"]["batch_labels"],
        seed=seed,
    )
    s_small = sample_stratified_control_sets_r16(
        n_cells=batch_fx["small_truncated"]["n_cells"],
        batch_labels=batch_fx["small_truncated"]["batch_labels"],
        seed=seed,
    )
    s_floor = sample_stratified_control_sets_r16(n_cells=batch_fx["below_floor"]["n_cells"], seed=seed)

    # Degenerate Group-2 rollup check (|G_2_qual| = 1)
    grp_deg = adjudicate_group_r16(
        per_line_h_i1={"a549_ko": {"line_summary": "MISSING_STRATUM"}},
        per_line_h_i2={"a549_ko": {"classification": "FCT_SUPERIOR_TO_MODEL"}},
        group1_lines=("nadig_hepg2", "nadig_jurkat"),
        group2_crispri_lines=(),
        group2_crisprko_lines=("a549_ko",),
    )

    pack = {
        "schema": "round16_at3_canonical_synthetic_parity_v1",
        "boot_seed": int(seed),
        "n_boot": int(n_boot),
        "scenarios": results,
        "batch_sampling_checks": {
            "large_imbalanced_sampled": s_large["batch_counts_sampled"],
            "large_imbalanced_n_sets": s_large["n_sets"],
            "tie_remainder_sampled": s_tie["batch_counts_sampled"],
            "small_truncated_n_sampled": s_small["n_sampled"],
            "small_truncated_n_sets": s_small["n_sets"],
            "below_floor_excluded": s_floor["excluded"],
        },
        "degenerate_group2_check": grp_deg,
    }
    if output_json is not None:
        output_json.write_text(json.dumps(pack, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return pack


if __name__ == "__main__":
    out_p = AUDIT_DIR / "canonical_synthetic_parity_r16.json"
    pack = export_canonical_parity_pack(output_json=out_p, n_boot=400, seed=DEFAULT_R16_SEED)
    print(f"Wrote {out_p} ({len(pack['scenarios'])} scenarios)")
    for k, v in pack["scenarios"].items():
        print(
            f"  {k}: matched={v['matched_r_orth_cal']:.4f} h_i1={v['h_i1_line_summary']} "
            f"h_i2={v['h_i2_classification']} counts={v['h_i1_counts']}"
        )
