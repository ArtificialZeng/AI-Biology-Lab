#!/usr/bin/env python3
"""
Clean-Room Round 16 Mechanical Adjudicator (`A-T3`, strictly `PREREG_I.md v5 §3` & Addendum A).

Implements:
  - Two co-primary claims (`H-I1`, `H-I2`) under Bonferroni `alpha_claim = 0.05 / 2 = 0.025`
    (`97.5%` CI for superiority/inferiority, `95.0%` CI for TOST equivalence at `delta = 0.04`).
  - Mutually exclusive & exhaustive 4-way per-pair TOST classification with `valid_frac >= 0.95`
    gated `within_delta` (`M2`) and explicit conservative zero-boundary tie precedence
    (`ci975_lo == 0.0 or ci975_hi == 0.0 -> INCONCLUSIVE`).
  - Mandatory two-stratum separation (`seen_by_checkpoint` vs `unseen_by_checkpoint`, `I14`),
    explicit `rho_ab_orth > 0` check (`M1`), `LINE_MATCHED_SPECIFIC` vs `LINE_MATCHED_SPECIFIC_SUB_DELTA`
    (`M2`), and `LINE_EQUIVALENT_NONZERO_SIGNAL` vs `LINE_EQUIVALENT_NULL_SIGNAL` (`M3`).
  - Symmetric 3-way `H-I2` classification, degenerate Group-2 safeguard (`SINGLE_STRATUM_DESCRIPTIVE_ONLY`),
    and strictly separate prospective Group 3 (`gse295214_h1`) reporting with `H1`/Replogle-Weissman disclosures.

Never reads or imports `round16/src/`.
"""

from __future__ import annotations

from typing import Any, Dict, Sequence

import numpy as np

DEFAULT_TOST_DELTA = 0.04
MIN_BOOT_VALID_FRAC = 0.95


def classify_pair_tost_r16(
    pair_summary: Dict[str, Any],
    delta: float = DEFAULT_TOST_DELTA,
) -> Dict[str, Any]:
    """
    Classify one `(target_line, mismatched_source)` contrast into mutually exclusive classes:
      - `INCONCLUSIVE` if `valid_frac < 0.95` or `97.5%` CI non-finite
      - `EXACT_INVARIANCE` (`Addendum A Item 10`) if `max_pred_abs_diff <= 1e-12` (treated as
        `EQUIVALENT` in line-level adjudication so exact `[0, 0]` invariance is not masked)
      - `INCONCLUSIVE` (`PREREG_I v5 §3.1`) if `ci975_lo == 0.0` or `ci975_hi == 0.0`
      - `MATCHED_BETTER` if `ci975_lo > 0.0`
      - `MISMATCHED_BETTER` if `ci975_hi < 0.0`
      - `EQUIVALENT` if `within_delta` (`95%` CI strictly inside `(-delta, +delta)`)
      - `INCONCLUSIVE` otherwise
    """
    lo975 = float(pair_summary.get("ci975_lo", float("nan")))
    hi975 = float(pair_summary.get("ci975_hi", float("nan")))
    lo95 = float(pair_summary.get("ci95_lo", float("nan")))
    hi95 = float(pair_summary.get("ci95_hi", float("nan")))
    est = float(pair_summary.get("estimate", float("nan")))
    vf = float(pair_summary.get("valid_frac", 0.0))
    max_abs = float(pair_summary.get("max_pred_abs_diff", float("nan")))
    exact_inv = bool(
        pair_summary.get("exact_invariance", False)
        or (np.isfinite(max_abs) and max_abs <= 1e-12)
    )

    d = float(delta)
    valid_ok = bool(np.isfinite(vf) and vf >= MIN_BOOT_VALID_FRAC)
    finite_975 = bool(np.isfinite(lo975) and np.isfinite(hi975))
    finite_95 = bool(np.isfinite(lo95) and np.isfinite(hi95))
    within_delta = bool(valid_ok and finite_975 and finite_95 and lo95 > -d and hi95 < d)
    # Conservative tie rule (PREREG_I §3.1): either endpoint of the 97.5% CI equal to 0.0
    # evaluates to INCONCLUSIVE prior to MATCHED_BETTER / MISMATCHED_BETTER / EQUIVALENT,
    # except exact invariance (max_abs_diff <= 1e-12, Addendum A Item 10).
    zero_tie_975 = bool(finite_975 and (lo975 == 0.0 or hi975 == 0.0))

    if not valid_ok or not finite_975:
        cat = "INCONCLUSIVE"
    elif exact_inv:
        cat = "EXACT_INVARIANCE"
    elif zero_tie_975:
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
        "exact_invariance": exact_inv,
        "max_pred_abs_diff": max_abs,
        "within_delta": within_delta,
        "delta": d,
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
    matched_positive_ci975: bool,
    delta: float,
) -> Dict[str, Any]:
    counts = {
        "MATCHED_BETTER": 0,
        "EQUIVALENT": 0,
        "EXACT_INVARIANCE": 0,
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
    any_within = bool(n_pairs > 0 and any(c["within_delta"] for c in pair_calls.values()))
    equiv_like = counts["EQUIVALENT"] + counts["EXACT_INVARIANCE"]

    if n_pairs == 0:
        status = "STRATUM_EMPTY"
    elif not matched_arm_valid:
        status = "STRATUM_INCONCLUSIVE"
    elif matched_positive_ci975 and counts["MATCHED_BETTER"] == n_pairs:
        status = "STRATUM_MATCHED_SPECIFIC_SUB_DELTA" if all_within else "STRATUM_MATCHED_SPECIFIC"
    elif equiv_like == n_pairs:
        status = (
            "STRATUM_EQUIVALENT_NONZERO_SIGNAL"
            if matched_positive_ci975
            else "STRATUM_EQUIVALENT_NULL_SIGNAL"
        )
    elif counts["MISMATCHED_BETTER"] == n_pairs:
        status = "STRATUM_MISMATCHED_BETTER"
    else:
        status = "STRATUM_MIXED_OR_INCONCLUSIVE"

    return {
        "n_pairs": n_pairs,
        "sources": sorted(stratum_sources),
        "counts": counts,
        "equiv_like_count": int(equiv_like),
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
    resample_comparison: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Adjudicate `H-I1` per `PREREG_I.md v5 §3.1` + `Addendum A Items 8–10` exact priority order:
      1. `MISSING_STRATUM`
      2. `INCONCLUSIVE` (`not matched_arm_valid`: `valid_frac < 0.95` or `rho_ab_orth <= 0` or non-finite CI)
      3. `LINE_MATCHED_SPECIFIC` (`matched_positive_ci975` & all `MATCHED_BETTER` & not `all_within_delta`)
      4. `LINE_MATCHED_SPECIFIC_SUB_DELTA` (`matched_positive_ci975` & all `MATCHED_BETTER` & `all_within_delta`)
      5. `LINE_EQUIVALENT_NONZERO_SIGNAL` (`matched_positive_ci975` & all `EQUIVALENT`/`EXACT_INVARIANCE` in both strata)
      6. `LINE_EQUIVALENT_NULL_SIGNAL` (`not matched_positive_ci975` & all `EQUIVALENT`/`EXACT_INVARIANCE` in both strata)
      7. `LINE_MISMATCHED_BETTER` (all `MISMATCHED_BETTER` in both strata)
      8. `MIXED_OR_INCONCLUSIVE`
    """
    seen_set = set(seen_sources)
    unseen_set = set(unseen_sources)
    overlap = seen_set & unseen_set
    if overlap:
        raise ValueError(f"Strata must be disjoint; overlap={sorted(overlap)}")

    per_src = paired_mismatched.get("per_source", {})
    expected_keys = seen_set | unseen_set
    observed_keys = set(per_src.keys())
    if observed_keys != expected_keys:
        raise ValueError(
            "paired_mismatched['per_source'] keys must match seen | unseen exactly "
            f"(missing={sorted(expected_keys - observed_keys)}, extra={sorted(observed_keys - expected_keys)})"
        )

    vf = float(orth_suite.get("valid_frac", 0.0))
    rho_orth = float(orth_suite.get("rho_ab_orth", float("nan")))
    lo975 = float(orth_suite.get("ci975_lo", float("nan")))
    hi975 = float(orth_suite.get("ci975_hi", float("nan")))
    matched_arm_valid = bool(
        np.isfinite(vf)
        and vf >= MIN_BOOT_VALID_FRAC
        and np.isfinite(rho_orth)
        and rho_orth > 0.0
        and np.isfinite(lo975)
        and np.isfinite(hi975)
    )
    matched_positive_ci975 = bool(matched_arm_valid and lo975 > 0.0)

    seen_res = _summarize_stratum_r16(
        seen_sources, per_src, matched_arm_valid, matched_positive_ci975, delta=delta
    )
    unseen_res = _summarize_stratum_r16(
        unseen_sources, per_src, matched_arm_valid, matched_positive_ci975, delta=delta
    )

    combined_counts = {
        k: seen_res["counts"][k] + unseen_res["counts"][k]
        for k in ("MATCHED_BETTER", "EQUIVALENT", "EXACT_INVARIANCE", "INCONCLUSIVE", "MISMATCHED_BETTER")
    }
    combined_per_source: Dict[str, Dict[str, Any]] = {}
    for sk, val in seen_res["per_source"].items():
        combined_per_source[sk] = {**val, "stratum": "seen_by_checkpoint"}
    for sk, val in unseen_res["per_source"].items():
        combined_per_source[sk] = {**val, "stratum": "unseen_by_checkpoint"}

    n_seen = seen_res["n_pairs"]
    n_unseen = unseen_res["n_pairs"]
    n_total = n_seen + n_unseen
    all_within_delta = bool(
        n_total > 0 and all(v["within_delta"] for v in combined_per_source.values())
    )
    any_within_delta = bool(
        n_total > 0 and any(v["within_delta"] for v in combined_per_source.values())
    )

    if n_seen == 0 or n_unseen == 0:
        line_summary = "MISSING_STRATUM"
    elif not matched_arm_valid:
        line_summary = "INCONCLUSIVE"
    elif (
        matched_positive_ci975
        and seen_res["counts"]["MATCHED_BETTER"] == n_seen
        and unseen_res["counts"]["MATCHED_BETTER"] == n_unseen
    ):
        line_summary = (
            "LINE_MATCHED_SPECIFIC_SUB_DELTA"
            if all_within_delta
            else "LINE_MATCHED_SPECIFIC"
        )
    elif (
        seen_res["equiv_like_count"] == n_seen
        and unseen_res["equiv_like_count"] == n_unseen
    ):
        line_summary = (
            "LINE_EQUIVALENT_NONZERO_SIGNAL"
            if matched_positive_ci975
            else "LINE_EQUIVALENT_NULL_SIGNAL"
        )
    elif (
        seen_res["counts"]["MISMATCHED_BETTER"] == n_seen
        and unseen_res["counts"]["MISMATCHED_BETTER"] == n_unseen
    ):
        line_summary = "LINE_MISMATCHED_BETTER"
    else:
        line_summary = "MIXED_OR_INCONCLUSIVE"

    interpretation_note: Optional[str] = None
    if line_summary == "LINE_MATCHED_SPECIFIC_SUB_DELTA":
        interpretation_note = (
            "能检测到，但小于事先规定的最小有意义效应 δ=0.04（须与同线重抽参照臂结果合并报告，不可单独写成‘模型用到了对照’）"
        )

    out: Dict[str, Any] = {
        "matched_arm_valid": matched_arm_valid,
        "matched_positive_ci975": matched_positive_ci975,
        "matched_rho_ab_orth": rho_orth,
        "matched_valid_frac": vf,
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
        "interpretation_note": interpretation_note,
    }
    if resample_comparison is not None:
        rc = dict(resample_comparison)
        if not matched_positive_ci975 and rc.get("overall_status") == "ALL_EXCEED_WITH_POSITIVE_SIGNAL":
            rc["overall_status"] = "ALL_EXCEED_WITHOUT_POSITIVE_SIGNAL"
            rc["overall_interpretation"] = "换成其他细胞系的对照，影响不超过重抽本细胞系的对照（或无正向匹配信号）"
        out["resample_reference"] = rc
    return out


def adjudicate_claim_i2_line(paired_vs_fct: Dict[str, Any]) -> Dict[str, Any]:
    """Symmetric `H-I2` (`STATE_SE - FCT`) classification at 97.5% Bonferroni CI (`PREREG_I v4 §3.2`)."""
    lo975 = float(paired_vs_fct.get("ci975_lo", float("nan")))
    hi975 = float(paired_vs_fct.get("ci975_hi", float("nan")))
    vf = float(paired_vs_fct.get("valid_frac", 0.0))
    if not (np.isfinite(vf) and vf >= MIN_BOOT_VALID_FRAC and np.isfinite(lo975) and np.isfinite(hi975)):
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


def adjudicate_group_r16(
    per_line_h_i1: Dict[str, Dict[str, Any]],
    per_line_h_i2: Dict[str, Dict[str, Any]],
    group1_lines: Sequence[str] = ("nadig_hepg2", "nadig_jurkat"),
    group2_crispri_lines: Sequence[str] = ("dld1", "cd4t_rest"),
    group2_crisprko_lines: Sequence[str] = ("a549_ko",),
    group3_lines: Sequence[str] = ("gse295214_h1",),
    group3_qualified: bool = True,
) -> Dict[str, Any]:
    """
    Group-level rollup strictly separated across Group 1 (`K1`), Group 2 CRISPRi,
    Group 2 CRISPRko, and Prospective Group 3 (`gse295214_h1`), enforcing the pre-registered
    degenerate Group-2 safeguard (`SINGLE_STRATUM_DESCRIPTIVE_ONLY` if `|G_2_qual| < 2`)
    per `PREREG_I.md v5 §2.4.5, §3.1`.
    """
    def _tally(lines: Sequence[str]) -> Dict[str, Any]:
        verdicts: Dict[str, str] = {}
        i2_verdicts: Dict[str, str] = {}
        counts: Dict[str, int] = {}
        for ln in lines:
            if ln in per_line_h_i1:
                v = str(per_line_h_i1[ln]["line_summary"])
                verdicts[ln] = v
                counts[v] = counts.get(v, 0) + 1
            if ln in per_line_h_i2:
                i2_verdicts[ln] = str(per_line_h_i2[ln]["classification"])
        return {
            "lines": [ln for ln in lines if ln in per_line_h_i1 or ln in per_line_h_i2],
            "K": len([ln for ln in lines if ln in per_line_h_i1]),
            "h_i1_verdicts": verdicts,
            "h_i1_counts": counts,
            "h_i2_verdicts": i2_verdicts,
        }

    g2_total_qual = sum(
        1 for ln in (list(group2_crispri_lines) + list(group2_crisprko_lines)) if ln in per_line_h_i1
    )
    g2_downgraded = g2_total_qual < 2

    g2_i_tally = _tally(group2_crispri_lines)
    g2_ko_tally = _tally(group2_crisprko_lines)
    if g2_downgraded:
        g2_i_tally["adjudication_status"] = "SINGLE_STRATUM_DESCRIPTIVE_ONLY"
        g2_ko_tally["adjudication_status"] = "SINGLE_STRATUM_DESCRIPTIVE_ONLY"
    else:
        g2_i_tally["adjudication_status"] = "RETROSPECTIVE_ADJUDICATED"
        g2_ko_tally["adjudication_status"] = "RETROSPECTIVE_ADJUDICATED"

    g3_tally = _tally(group3_lines)
    if g3_tally["K"] == 0:
        g3_status = "NOT_EVALUATED_OR_OMITTED"
    elif not group3_qualified:
        g3_status = "SENSITIVITY_ONLY"
    else:
        g3_status = "PROSPECTIVE_HOLDOUT_SEPARATE"
    g3_tally["adjudication_status"] = g3_status
    g3_tally["merged_with_group1_or_group2"] = False
    g3_tally["mandatory_disclosures"] = [
        "shared H1 hESC lineage with Round 9 vcc2025_h1",
        "shared Replogle/Weissman laboratory provenance with Replogle-4",
    ]

    return {
        "group1_replogle4_loo": {**_tally(group1_lines), "note_n_seen": 1, "adjudication_status": "RETROSPECTIVE_ADJUDICATED"},
        "group2_crispri": g2_i_tally,
        "group2_crisprko": g2_ko_tally,
        "group3_prospective": g3_tally,
        "group2_qualifying_count": g2_total_qual,
        "group2_downgraded_single_stratum": g2_downgraded,
    }
