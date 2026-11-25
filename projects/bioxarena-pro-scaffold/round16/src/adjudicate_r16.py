#!/usr/bin/env python3
"""
Round 16 Mechanical Adjudicator (`PREREG_I v5 §3 + Addendum A Items 5, 8, 9, 10`).

Applies exact 2-claim Bonferroni (`alpha_claim = 0.05 / 2 = 0.025`) adjudication for:
  - `H-I1` (Primary model `STATE_SE`, matched vs. mismatched target controls under `r_orth_cal`
    on `all_nontarget`, stratified into `seen_by_checkpoint` and `unseen_by_checkpoint` local
    raw-count mismatched pools per `I14` & `Q4`):
      * Orthogonal split-half denominator strictly keyed `rho_ab_orth` (`M1`).
      * Two-sided `97.5%` CI (`[1.25%, 98.75%]`) for superiority/inferiority.
      * Two-sided `95.0%` CI (`[2.5%, 97.5%]`) inside `(-delta, +delta)` for TOST equivalence
        (`delta = 0.04`), with `within_delta` gated on `valid_frac >= 0.95` and finite CIs (`M2`).
      * Per-pair mutually exclusive categories (Addendum A Item 10):
        `MATCHED_BETTER` / `EQUIVALENT` / `EXACT_INVARIANCE` / `INCONCLUSIVE` / `MISMATCHED_BETTER`,
        where `EXACT_INVARIANCE` (`max_abs_diff <= 1e-12` between matched and mismatched predicted
        delta arrays under valid bootstrap) counts identically as `EQUIVALENT` in stratum and line
        rules, while non-identical zero-boundary ties (`ci975_lo == 0.0` or `ci975_hi == 0.0`) remain
        conservatively `INCONCLUSIVE`.
      * Mandatory two-stratum reporting and line-level verdicts (`M2`, `M3`, Addendum A Item 9):
        `MISSING_STRATUM` | `INCONCLUSIVE` | `LINE_MATCHED_SPECIFIC` |
        `LINE_MATCHED_SPECIFIC_SUB_DELTA` | `LINE_EQUIVALENT_NONZERO_SIGNAL` |
        `LINE_EQUIVALENT_NULL_SIGNAL` | `LINE_MISMATCHED_BETTER` | `MIXED_OR_INCONCLUSIVE`.
      * Group-level `k / K` reporting stratified by `(group, modality)` without data pooling,
        and pre-locked Group 2 `< 2` single-stratum descriptive safeguard (`M3`).
  - `H-I2` (`STATE_SE` vs. `FCT` under `r_orth_cal` on `all_nontarget`):
      * Two-sided `97.5%` CI symmetric 3-way classification:
        `MODEL_SUPERIOR_TO_FCT` / `FCT_SUPERIOR_TO_MODEL` / `INDISTINGUISHABLE`.

Secondary analyses (`r_dorth_cal`, Mechanism Hypothesis `M2`, `STATE_HVG`, `G4` 2x2 decomposition,
and Addendum A Item 8 within-line resample reference) are strictly descriptive and never alter
primary `H-I1`/`H-I2` adjudication. `CellularAI_OpenCode_PublicRefit` is dropped pre-data
(Addendum A Item 5).
"""

from __future__ import annotations

from typing import Any, Dict, Mapping, Sequence

import numpy as np

DEFAULT_TOST_DELTA = 0.04
MIN_BOOT_VALID_FRAC = 0.95
MIN_GROUP2_PRIMARY_LINES = 2
EXACT_INVARIANCE_TOL = 1e-12

SUB_DELTA_INTERPRETATION_NOTE = (
    "能检测到，但小于事先规定的最小有意义效应 δ=0.04"
    "（须结合同线重抽参照臂结果一并报告，不能单独写成「模型用到了对照」）。"
)
RESAMPLE_EXCEEDS_WORDING = "模型用到了目标线特有的对照信息"
RESAMPLE_NOT_EXCEEDS_WORDING = "换成其他细胞系的对照，影响不超过重抽本细胞系的对照"


def classify_pair_tost_r16(
    pair_summary: Dict[str, Any],
    delta: float = DEFAULT_TOST_DELTA,
    exact_invariance_tol: float = EXACT_INVARIANCE_TOL,
) -> Dict[str, Any]:
    """
    Classify a single `(line t, mismatched source s')` paired contrast into mutually exclusive
    categories (`PREREG_I v5 §3.1` + `Addendum A Item 10`):
      - `within_delta`: True iff `valid_frac >= 0.95`, finite `97.5%` and `95.0%` CIs,
        AND strict two-sided `ci95_lo > -delta and ci95_hi < +delta`.
      - Priority order:
        1. `INCONCLUSIVE`: `valid_frac < 0.95` or non-finite `97.5%` CI (`within_delta = False`).
        2. `EXACT_INVARIANCE` (Addendum A Item 10): finite matched and mismatched predicted delta
           arrays are identical (`max_abs_diff <= 1e-12`, or `exact_invariance == True` when
           `max_abs_diff` is omitted); counts identically as `EQUIVALENT` at stratum/line level.
        3. `INCONCLUSIVE` (conservative zero-boundary tie precedence for non-identical arrays):
           `ci975_lo == 0.0` or `ci975_hi == 0.0` (takes precedence over `EQUIVALENT` even when
           `within_delta == True`).
        4. `MATCHED_BETTER`: `ci975_lo > 0.0` (accompanied by gated two-sided `within_delta`).
        5. `MISMATCHED_BETTER`: `ci975_hi < 0.0` (accompanied by gated two-sided `within_delta`).
        6. `EQUIVALENT`: NOT matched/mismatched better, not exact zero-tie, AND `within_delta == True`.
        7. `INCONCLUSIVE`: all remaining cases (including exact `±delta` open-interval ties).
    """
    lo975 = float(pair_summary.get("ci975_lo", float("nan")))
    hi975 = float(pair_summary.get("ci975_hi", float("nan")))
    lo95 = float(pair_summary.get("ci95_lo", float("nan")))
    hi95 = float(pair_summary.get("ci95_hi", float("nan")))
    est = float(pair_summary.get("estimate", float("nan")))
    vf = float(pair_summary.get("valid_frac", 0.0))

    max_abs_diff_raw = pair_summary.get("max_abs_diff", None)
    exact_inv_raw = pair_summary.get("exact_invariance", None)
    if max_abs_diff_raw is not None:
        mad_f = float(max_abs_diff_raw)
        is_exact_inv = bool(np.isfinite(mad_f) and 0.0 <= mad_f <= float(exact_invariance_tol))
    elif exact_inv_raw is not None:
        mad_f = float("nan")
        is_exact_inv = bool(exact_inv_raw)
    else:
        mad_f = float("nan")
        is_exact_inv = False

    finite_975 = bool(np.isfinite(lo975) and np.isfinite(hi975))
    finite_95 = bool(np.isfinite(lo95) and np.isfinite(hi95))
    valid_boot = bool(vf >= MIN_BOOT_VALID_FRAC and finite_975)
    within_delta = bool(
        valid_boot and finite_95 and lo95 > -float(delta) and hi95 < float(delta)
    )

    if not valid_boot:
        cat = "INCONCLUSIVE"
    elif is_exact_inv:
        cat = "EXACT_INVARIANCE"
    elif lo975 == 0.0 or hi975 == 0.0:
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
        "exact_invariance": is_exact_inv,
        "max_abs_diff": mad_f,
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
    any_within = bool(any(c["within_delta"] for c in pair_calls.values()))
    # Addendum A Item 10: EXACT_INVARIANCE counts identically as EQUIVALENT in stratum/line rules.
    n_equiv_effective = counts["EQUIVALENT"] + counts["EXACT_INVARIANCE"]

    if n_pairs == 0:
        status = "STRATUM_EMPTY"
    elif not matched_arm_valid:
        status = "STRATUM_INCONCLUSIVE"
    elif matched_positive and counts["MATCHED_BETTER"] == n_pairs:
        status = "STRATUM_MATCHED_SPECIFIC_SUB_DELTA" if all_within else "STRATUM_MATCHED_SPECIFIC"
    elif n_equiv_effective == n_pairs:
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
        for k in ("MATCHED_BETTER", "EQUIVALENT", "EXACT_INVARIANCE", "INCONCLUSIVE", "MISMATCHED_BETTER")
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

    if line_summary == "LINE_MATCHED_SPECIFIC_SUB_DELTA":
        interp_note = SUB_DELTA_INTERPRETATION_NOTE
    elif line_summary == "LINE_EQUIVALENT_NONZERO_SIGNAL":
        interp_note = (
            "Non-zero matched signal with matched-vs-mismatched equivalence within ±delta "
            "(EXACT_INVARIANCE pairs counted as EQUIVALENT per Addendum A Item 10)."
        )
    elif line_summary == "LINE_EQUIVALENT_NULL_SIGNAL":
        interp_note = (
            "Trivial null-signal equivalence (matched arm not strictly positive; NEVER report as "
            "non-zero signal ignoring target controls)."
        )
    elif line_summary in ("INCONCLUSIVE", "MIXED_OR_INCONCLUSIVE"):
        interp_note = (
            "Any INCONCLUSIVE pair reflects insufficient statistical power, NEVER equivalence ruled out "
            "(PREREG_I §3.3)."
        )
    else:
        interp_note = ""

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
        "interpretation_note": interp_note,
    }


def summarize_within_line_resample_r16(
    resample_result: Dict[str, Any],
    seen_sources: Sequence[str],
    unseen_sources: Sequence[str],
    matched_positive_ci975: bool = True,
) -> Dict[str, Any]:
    """
    Secondary descriptive interpretation helper for `PREREG_I Addendum A Items 8 & 9`.
    Strictly descriptive — NEVER modifies primary `H-I1` or `H-I2` adjudication.
    """
    seen_set = set(seen_sources)
    unseen_set = set(unseen_sources)
    overlap = seen_set & unseen_set
    if overlap:
        raise ValueError(f"Mismatched strata must be disjoint; overlap={sorted(overlap)}")

    per_src_in = resample_result.get("per_source", {})
    observed_keys = set(per_src_in.keys())
    expected_keys = seen_set | unseen_set
    if observed_keys != expected_keys:
        raise ValueError(
            "resample_result['per_source'] keys must match seen_sources | unseen_sources exactly "
            f"(missing={sorted(expected_keys - observed_keys)}, extra={sorted(observed_keys - expected_keys)})."
        )

    per_source_out: Dict[str, Dict[str, Any]] = {}
    for stratum_name, s_list in (
        ("seen_by_checkpoint", sorted(seen_sources)),
        ("unseen_by_checkpoint", sorted(unseen_sources)),
    ):
        for sk in s_list:
            exc = per_src_in[sk]["delta_pair_excess_over_resample"]
            vf = float(exc.get("valid_frac", 0.0))
            lo975 = float(exc.get("ci975_lo", float("nan")))
            hi975 = float(exc.get("ci975_hi", float("nan")))
            lo95 = float(exc.get("ci95_lo", float("nan")))
            hi95 = float(exc.get("ci95_hi", float("nan")))
            est = float(exc.get("estimate", float("nan")))
            valid_boot_exc = bool(
                vf >= MIN_BOOT_VALID_FRAC
                and np.isfinite(lo975)
                and np.isfinite(hi975)
            )
            exceeds = bool(valid_boot_exc and lo975 > 0.0)
            if not valid_boot_exc:
                pair_interp_zh = "同线重抽比较无效（bootstrap 有效率不足或非有限值）"
            elif exceeds:
                pair_interp_zh = RESAMPLE_EXCEEDS_WORDING
            else:
                pair_interp_zh = RESAMPLE_NOT_EXCEEDS_WORDING
            per_source_out[sk] = {
                "stratum": stratum_name,
                "excess_estimate": est,
                "excess_ci975_lo": lo975,
                "excess_ci975_hi": hi975,
                "excess_ci95_lo": lo95,
                "excess_ci95_hi": hi95,
                "valid_frac": vf,
                "valid_boot": valid_boot_exc,
                "exceeds_resample_loss": exceeds,
                "pair_interpretation_zh": pair_interp_zh,
            }

    n_total = len(per_source_out)
    n_valid = sum(1 for v in per_source_out.values() if v["valid_boot"])
    n_exceed = sum(1 for v in per_source_out.values() if v["exceeds_resample_loss"])
    all_valid = bool(n_total > 0 and n_valid == n_total)
    all_exceed = bool(all_valid and n_exceed == n_total)
    none_exceed = bool(all_valid and n_exceed == 0)

    if n_total == 0:
        line_interp_zh = "无错配来源可与同线重抽参照臂比较"
        line_status = "NO_MISMATCHED_SOURCES"
    elif not all_valid:
        line_interp_zh = "部分或全部错配来源的同线重抽比较无效（bootstrap 有效率不足或非有限值）"
        line_status = "INCONCLUSIVE_VS_WITHIN_LINE_RESAMPLE"
    elif all_exceed and bool(matched_positive_ci975):
        line_interp_zh = RESAMPLE_EXCEEDS_WORDING
        line_status = "EXCEEDS_WITHIN_LINE_RESAMPLE"
    elif all_exceed:
        line_interp_zh = "所有错配来源损失超过同线重抽，但匹配臂未显著为正（不作目标线特异对照结论）"
        line_status = "EXCEEDS_RESAMPLE_BUT_MATCHED_NOT_POSITIVE"
    elif none_exceed:
        line_interp_zh = RESAMPLE_NOT_EXCEEDS_WORDING
        line_status = "NOT_EXCEEDING_WITHIN_LINE_RESAMPLE"
    else:
        line_interp_zh = "部分错配来源超过同线重抽损失、部分未超过（按来源逐条报告，不作全线目标特异对照结论）"
        line_status = "MIXED_VS_WITHIN_LINE_RESAMPLE"

    return {
        "role": "secondary_descriptive_only_does_not_alter_primary_H_I1_H_I2",
        "n_pairs": n_total,
        "n_valid_pairs": n_valid,
        "n_exceeds_resample_loss": n_exceed,
        "all_exceeds_resample_loss": all_exceed,
        "none_exceeds_resample_loss": none_exceed,
        "resample_comparison_status": line_status,
        "interpretation_zh": line_interp_zh,
        "delta_pair_resample": resample_result.get("delta_pair_resample", {}),
        "per_source": per_source_out,
    }


def adjudicate_group_r16(
    group_key: str,
    line_verdicts: Mapping[str, str],
    line_modalities: Mapping[str, str],
    min_group2_lines: int = MIN_GROUP2_PRIMARY_LINES,
) -> Dict[str, Any]:
    """
    Group-level `H-I1` summary (`PREREG_I v5 §3.1 + Addendum A`, audit item `M3`):
      - Reports `k / K` counts separately by `(group_key, modality)` (`CRISPRi` vs `CRISPRko`),
        strictly without pooling perturbations across lines or merging modalities or merging
        Group 3 (`gse295214_h1`) with Groups 1/2.
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
