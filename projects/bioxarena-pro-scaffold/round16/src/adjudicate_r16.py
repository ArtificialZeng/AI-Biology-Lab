#!/usr/bin/env python3
"""
Round 16 Mechanical Adjudicator (`PREREG_I v4 §3`, incorporating `I14`, `M1–M6`, `Q4`).

Applies exact 2-claim Bonferroni (`alpha_claim = 0.05 / 2 = 0.025`) adjudication for:
  - `H-I1` (Primary model `STATE_SE`, matched vs. mismatched target controls under `r_orth_cal`
    on `all_nontarget`, stratified into `seen_by_checkpoint` and `unseen_by_checkpoint` local
    raw-count mismatched pools per `I14` & `Q4`):
      * Orthogonal split-half denominator strictly keyed `rho_ab_orth` (`M1`).
      * Two-sided `97.5%` CI (`[1.25%, 98.75%]`) for superiority/inferiority.
      * Two-sided `95.0%` CI (`[2.5%, 97.5%]`) inside `(-delta, +delta)` for TOST equivalence
        (`delta = 0.04`), with `within_delta` gated on `valid_frac >= 0.95` and finite CIs (`M2`).
      * 4 mutually exclusive per-pair categories:
        `MATCHED_BETTER` / `EQUIVALENT` / `INCONCLUSIVE` / `MISMATCHED_BETTER`.
      * Mandatory two-stratum reporting and line-level verdicts (`M2`, `M3`):
        `MISSING_STRATUM` | `INCONCLUSIVE` | `LINE_MATCHED_SPECIFIC` |
        `LINE_MATCHED_SPECIFIC_SUB_DELTA` | `LINE_EQUIVALENT_NONZERO_SIGNAL` |
        `LINE_EQUIVALENT_NULL_SIGNAL` | `LINE_MISMATCHED_BETTER` | `MIXED_OR_INCONCLUSIVE`.
      * Group-level `k / K` reporting stratified by `(group, modality)` without data pooling,
        and pre-locked Group 2 `< 2` single-stratum descriptive safeguard (`M3`).
  - `H-I2` (`STATE_SE` vs. `FCT` under `r_orth_cal` on `all_nontarget`):
      * Two-sided `97.5%` CI symmetric 3-way classification:
        `MODEL_SUPERIOR_TO_FCT` / `FCT_SUPERIOR_TO_MODEL` / `INDISTINGUISHABLE`.

Secondary analyses (`r_dorth_cal`, Mechanism Hypothesis `M2`, `STATE_HVG`,
`CellularAI_OpenCode_PublicRefit`) are strictly descriptive and are never adjudicated here.
"""

from __future__ import annotations

from typing import Any, Dict, Mapping, Sequence

import numpy as np

DEFAULT_TOST_DELTA = 0.04
MIN_BOOT_VALID_FRAC = 0.95
MIN_GROUP2_PRIMARY_LINES = 2


def classify_pair_tost_r16(
    pair_summary: Dict[str, Any],
    delta: float = DEFAULT_TOST_DELTA,
) -> Dict[str, Any]:
    """
    Classify a single `(line t, mismatched source s')` paired contrast into exactly one
    of 4 mutually exclusive categories (`PREREG_I v5 §3.1`, audit items `I6`, `M2`,
    and conservative `ci975 == 0.0` tie-precedence ruling `2026-11-25`):
      - `within_delta`: True iff `valid_frac >= 0.95`, finite `97.5%` and `95.0%` CIs,
        AND strict two-sided `ci95_lo > -delta and ci95_hi < +delta`.
      - `INCONCLUSIVE` (highest precedence): `valid_frac < 0.95` or non-finite `97.5%` CI
        (with `within_delta=False`), OR either endpoint of the `97.5%` CI equals `0.0`
        (`ci975_lo == 0.0` or `ci975_hi == 0.0`, taking precedence over `EQUIVALENT` even
        when `within_delta == True`).
      - `MATCHED_BETTER`: `ci975_lo > 0.0` (accompanied by gated two-sided `within_delta`)
      - `MISMATCHED_BETTER`: `ci975_hi < 0.0` (accompanied by gated two-sided `within_delta`)
      - `EQUIVALENT`: NOT matched/mismatched better, neither `97.5%` CI endpoint equals `0.0`,
        AND `within_delta == True`
      - `INCONCLUSIVE`: all remaining cases (including exact `±delta` open-interval ties)
    """
    lo975 = float(pair_summary.get("ci975_lo", float("nan")))
    hi975 = float(pair_summary.get("ci975_hi", float("nan")))
    lo95 = float(pair_summary.get("ci95_lo", float("nan")))
    hi95 = float(pair_summary.get("ci95_hi", float("nan")))
    est = float(pair_summary.get("estimate", float("nan")))
    vf = float(pair_summary.get("valid_frac", 0.0))

    finite_975 = bool(np.isfinite(lo975) and np.isfinite(hi975))
    finite_95 = bool(np.isfinite(lo95) and np.isfinite(hi95))
    valid_boot = bool(vf >= MIN_BOOT_VALID_FRAC and finite_975)
    within_delta = bool(
        valid_boot and finite_95 and lo95 > -float(delta) and hi95 < float(delta)
    )

    if not valid_boot or lo975 == 0.0 or hi975 == 0.0:
        cat = "INCONCLUSIVE"
    elif lo975 > 0.0:
        cat = "MATCHED_BETTER"
    elif hi975 < 0.0:
        cat = "MISMATCHED_BETTER"
    elif within_delta:
        cat = "EQUIVALENT"
    else:
        cat = "INCONCLUSIVE"

    return {
        "category": cat,
        "within_delta": within_delta,
        "delta": float(delta),
        "estimate": est,
        "ci975_lo": lo975,
        "ci975_hi": hi975,
        "ci95_lo": lo95,
        "ci95_hi": hi95,
        "valid_frac": vf,
    }


def _summarize_stratum_r16(
    stratum_sources: Sequence[str],
    per_source_all: Dict[str, Dict[str, Any]],
    matched_arm_valid: bool,
    matched_positive: bool,
    delta: float,
) -> Dict[str, Any]:
    """Summarize one mismatched-control stratum (`seen_by_checkpoint` or `unseen_by_checkpoint`)."""
    counts = {
        "MATCHED_BETTER": 0,
        "EQUIVALENT": 0,
        "INCONCLUSIVE": 0,
        "MISMATCHED_BETTER": 0,
    }
    pair_calls: Dict[str, Dict[str, Any]] = {}
    for sk in sorted(stratum_sources):
        c = classify_pair_tost_r16(per_source_all[sk], delta=delta)
        pair_calls[sk] = c
        counts[c["category"]] += 1

    n_pairs = len(pair_calls)
    all_within = bool(n_pairs > 0 and all(c["within_delta"] for c in pair_calls.values()))
    any_within = bool(any(c["within_delta"] for c in pair_calls.values()))

    if n_pairs == 0:
        status = "STRATUM_EMPTY"
    elif not matched_arm_valid:
        status = "STRATUM_INCONCLUSIVE"
    elif matched_positive and counts["MATCHED_BETTER"] == n_pairs:
        status = "STRATUM_MATCHED_SPECIFIC_SUB_DELTA" if all_within else "STRATUM_MATCHED_SPECIFIC"
    elif counts["EQUIVALENT"] == n_pairs:
        status = "STRATUM_EQUIVALENT_NONZERO_SIGNAL" if matched_positive else "STRATUM_EQUIVALENT_NULL_SIGNAL"
    elif counts["MISMATCHED_BETTER"] == n_pairs:
        status = "STRATUM_MISMATCHED_BETTER"
    else:
        status = "STRATUM_MIXED_OR_INCONCLUSIVE"

    return {
        "n_pairs": n_pairs,
        "sources": sorted(stratum_sources),
        "counts": counts,
        "all_within_delta": all_within,
        "any_within_delta": any_within,
        "per_source": pair_calls,
        "stratum_summary": status,
    }


def adjudicate_claim_i1_line(
    orth_suite: Dict[str, Any],
    paired_mismatched: Dict[str, Any],
    seen_sources: Sequence[str],
    unseen_sources: Sequence[str],
    delta: float = DEFAULT_TOST_DELTA,
) -> Dict[str, Any]:
    """
    Adjudicate `H-I1` for primary `STATE_SE` on a single qualifying line under `r_orth_cal`
    (`all_nontarget`) using 2-claim Bonferroni (`97.5%` superiority + `95%` TOST equivalence),
    orthogonal split-half key `rho_ab_orth` (`M1`), gated two-sided `within_delta` (`M2`),
    matched-arm validity & equivalent-signal splitting (`M3`), and two-stratum separation (`I14`).

    Exact priority order:
      1. `MISSING_STRATUM`: either `seen_by_checkpoint` or `unseen_by_checkpoint` is empty.
      2. `INCONCLUSIVE`: matched arm `valid_frac < 0.95` or `rho_ab_orth <= 0` or non-finite CI.
      3. `LINE_MATCHED_SPECIFIC_SUB_DELTA`: `matched_positive_ci975` AND both strata all
         `MATCHED_BETTER` AND `all_within_delta == True` (all pairs strictly inside `±delta`).
      4. `LINE_MATCHED_SPECIFIC`: `matched_positive_ci975` AND both strata all `MATCHED_BETTER`
         AND `all_within_delta == False`.
      5. `LINE_EQUIVALENT_NONZERO_SIGNAL`: `matched_positive_ci975` AND both strata all `EQUIVALENT`
         (ONLY this label supports "non-zero above-transfer-mean signal ignoring target controls").
      6. `LINE_EQUIVALENT_NULL_SIGNAL`: `NOT matched_positive_ci975` AND both strata all `EQUIVALENT`
         (trivial null-signal equivalence; matched arm indistinguishable from 0 or negative).
      7. `LINE_MISMATCHED_BETTER`: both strata all `MISMATCHED_BETTER`.
      8. `MIXED_OR_INCONCLUSIVE`: all remaining cases (including cross-stratum disagreement).
    """
    seen_set = set(seen_sources)
    unseen_set = set(unseen_sources)
    overlap = seen_set & unseen_set
    if overlap:
        raise ValueError(f"Mismatched strata must be disjoint; overlap={sorted(overlap)}")

    per_src = paired_mismatched.get("per_source", {})
    observed_keys = set(per_src.keys())
    expected_keys = seen_set | unseen_set
    if observed_keys != expected_keys:
        missing = sorted(expected_keys - observed_keys)
        extra = sorted(observed_keys - expected_keys)
        raise ValueError(
            "paired_mismatched['per_source'] keys must match seen_sources | unseen_sources exactly "
            f"(missing={missing}, extra={extra})."
        )

    vf = float(orth_suite.get("valid_frac", 0.0))
    # M1: strictly require orthogonal split-half key `rho_ab_orth` (no fallback to legacy `rho_ab`).
    rho = float(orth_suite.get("rho_ab_orth", float("nan")))
    lo975 = float(orth_suite.get("ci975_lo", float("nan")))
    hi975 = float(orth_suite.get("ci975_hi", float("nan")))
    matched_arm_valid = bool(
        vf >= MIN_BOOT_VALID_FRAC
        and np.isfinite(rho)
        and rho > 0.0
        and np.isfinite(lo975)
        and np.isfinite(hi975)
    )
    matched_positive = bool(matched_arm_valid and lo975 > 0.0)

    seen_res = _summarize_stratum_r16(
        seen_sources, per_src, matched_arm_valid, matched_positive, delta=delta
    )
    unseen_res = _summarize_stratum_r16(
        unseen_sources, per_src, matched_arm_valid, matched_positive, delta=delta
    )

    combined_counts = {
        k: seen_res["counts"][k] + unseen_res["counts"][k]
        for k in ("MATCHED_BETTER", "EQUIVALENT", "INCONCLUSIVE", "MISMATCHED_BETTER")
    }
    combined_per_source: Dict[str, Dict[str, Any]] = {}
    for sk, val in seen_res["per_source"].items():
        combined_per_source[sk] = {**val, "stratum": "seen_by_checkpoint"}
    for sk, val in unseen_res["per_source"].items():
        combined_per_source[sk] = {**val, "stratum": "unseen_by_checkpoint"}

    n_total = seen_res["n_pairs"] + unseen_res["n_pairs"]
    all_within_delta = bool(
        n_total > 0 and all(v["within_delta"] for v in combined_per_source.values())
    )
    any_within_delta = bool(any(v["within_delta"] for v in combined_per_source.values()))

    s_stat = seen_res["stratum_summary"]
    u_stat = unseen_res["stratum_summary"]
    mb_states = ("STRATUM_MATCHED_SPECIFIC", "STRATUM_MATCHED_SPECIFIC_SUB_DELTA")

    if s_stat == "STRATUM_EMPTY" or u_stat == "STRATUM_EMPTY":
        line_summary = "MISSING_STRATUM"
    elif not matched_arm_valid:
        line_summary = "INCONCLUSIVE"
    elif s_stat in mb_states and u_stat in mb_states:
        line_summary = (
            "LINE_MATCHED_SPECIFIC_SUB_DELTA" if all_within_delta else "LINE_MATCHED_SPECIFIC"
        )
    elif s_stat == "STRATUM_EQUIVALENT_NONZERO_SIGNAL" and u_stat == "STRATUM_EQUIVALENT_NONZERO_SIGNAL":
        line_summary = "LINE_EQUIVALENT_NONZERO_SIGNAL"
    elif s_stat == "STRATUM_EQUIVALENT_NULL_SIGNAL" and u_stat == "STRATUM_EQUIVALENT_NULL_SIGNAL":
        line_summary = "LINE_EQUIVALENT_NULL_SIGNAL"
    elif s_stat == "STRATUM_MISMATCHED_BETTER" and u_stat == "STRATUM_MISMATCHED_BETTER":
        line_summary = "LINE_MISMATCHED_BETTER"
    else:
        line_summary = "MIXED_OR_INCONCLUSIVE"

    return {
        "matched_arm_valid": matched_arm_valid,
        "matched_positive_ci975": matched_positive,
        "rho_ab_orth": rho,
        "n_pairs": n_total,
        "delta": float(delta),
        "all_within_delta": all_within_delta,
        "any_within_delta": any_within_delta,
        "strata": {
            "seen_by_checkpoint": seen_res,
            "unseen_by_checkpoint": unseen_res,
        },
        "counts": combined_counts,
        "per_source": combined_per_source,
        "line_summary": line_summary,
    }


def adjudicate_group_r16(
    group_key: str,
    line_verdicts: Mapping[str, str],
    line_modalities: Mapping[str, str],
    min_group2_lines: int = MIN_GROUP2_PRIMARY_LINES,
) -> Dict[str, Any]:
    """
    Group-level `H-I1` summary (`PREREG_I v4 §3.1`, audit item `M3`):
      - Reports `k / K` counts separately by `(group_key, modality)` (`CRISPRi` vs `CRISPRko`),
        strictly without pooling perturbations across lines or merging modalities.
      - Pre-locked Group 2 degenerate safeguard: if `group_key == "group2"` and fewer than
        `min_group2_lines` (`< 2`) primary target lines qualify, Group 2 is downgraded to
        `SINGLE_STRATUM_DESCRIPTIVE_ONLY` and no line-level verdicts are awarded.
    """
    if set(line_verdicts.keys()) != set(line_modalities.keys()):
        raise ValueError("line_verdicts and line_modalities must share identical line keys.")

    n_lines = len(line_verdicts)
    downgraded = bool(group_key == "group2" and n_lines < int(min_group2_lines))

    by_modality: Dict[str, Dict[str, Any]] = {}
    for line, verdict in sorted(line_verdicts.items()):
        mod = str(line_modalities[line])
        entry = by_modality.setdefault(mod, {"K": 0, "lines": [], "verdict_counts": {}, "k_over_K": {}})
        entry["K"] += 1
        entry["lines"].append(line)
        v_key = "DESCRIPTIVE_ONLY_GROUP2_LT2" if downgraded else str(verdict)
        entry["verdict_counts"][v_key] = entry["verdict_counts"].get(v_key, 0) + 1

    for entry in by_modality.values():
        k_tot = entry["K"]
        entry["k_over_K"] = {v: f"{cnt}/{k_tot}" for v, cnt in sorted(entry["verdict_counts"].items())}

    return {
        "group_key": str(group_key),
        "K_total": n_lines,
        "line_level_verdicts_awarded": not downgraded,
        "group_status": (
            "SINGLE_STRATUM_DESCRIPTIVE_ONLY" if downgraded else "ADJUDICATED_BY_MODALITY"
        ),
        "by_modality": by_modality,
    }


def adjudicate_claim_i2_line(paired_vs_fct: Dict[str, Any]) -> Dict[str, Any]:
    """
    Adjudicate `H-I2` (`STATE_SE - FCT` under `r_orth_cal`, `all_nontarget`) symmetrically
    at 2-claim Bonferroni two-sided `97.5%` CI.
    """
    lo975 = float(paired_vs_fct.get("ci975_lo", float("nan")))
    hi975 = float(paired_vs_fct.get("ci975_hi", float("nan")))
    vf = float(paired_vs_fct.get("valid_frac", 0.0))
    if vf < MIN_BOOT_VALID_FRAC or not (np.isfinite(lo975) and np.isfinite(hi975)):
        cls = "INDISTINGUISHABLE"
    elif lo975 > 0.0:
        cls = "MODEL_SUPERIOR_TO_FCT"
    elif hi975 < 0.0:
        cls = "FCT_SUPERIOR_TO_MODEL"
    else:
        cls = "INDISTINGUISHABLE"
    return {
        "classification": cls,
        "estimate": float(paired_vs_fct.get("estimate", float("nan"))),
        "ci975_lo": lo975,
        "ci975_hi": hi975,
        "valid_frac": vf,
    }
