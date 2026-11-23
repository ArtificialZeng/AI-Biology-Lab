#!/usr/bin/env python3
"""
Round 15 (`PREREG_H` v2) Confirmatory Adjudicator, Multi-Wave Aggregator, & Markdown Reporter.

Enforces `PREREG_H §0, §3, §8, §9, §10` using frozen `round12` IUT + 2-Claim Bonferroni
(`97.5%` CIs, `seed = 20261231`), supporting:
  - Per-wave evaluation (`§8.4`: each release wave adjudicated independently and labelled
    `"interim, more lines pending"` until all 5 candidate lines are resolved or sunset
    `2027-06-30` is reached).
  - Cross-wave IUT aggregation (`aggregate_r15_wave_adjudications`): overall confirmation
    requires unanimous `K_conf_total / K_conf_total` (`1 <= K_conf_total <= 5`) passage
    across ALL qualifying lines of ALL waves; a later wave can only keep or revoke an
    earlier confirmation, never rescue an earlier wave's failure.
  - Mandatory `§8.2` scope disclosures (`K_conf=0` aborted, `K_conf=1` single-line limited
    generalizability, `K_conf>=2` single-lab/platform `Arc VCC 2026, 10x Flex CRISPRi`).
  - Strictly descriptive `§9` reporting (`SA-H1` leave-`vcc2025_h1`-out & `SA-Panel` `129`
    val-test overlap vs remainder).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

ROUND15_SRC_DIR = Path(__file__).resolve().parent
REPO_ROOT = ROUND15_SRC_DIR.parents[1]
R12_SRC_DIR = REPO_ROOT / "round12" / "src"
for _p in (str(ROUND15_SRC_DIR), str(R12_SRC_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from adjudicate_r12 import (  # noqa: E402
    _eval_claim1_he1,
    _eval_claim2_he2,
    _evaluate_util_noninferiority,
    render_r12_markdown_report,
)
from data_loader_r15 import (  # noqa: E402
    ALLOWED_R15_LINES,
    ALLOWED_R15_LINES_SET,
    EXCLUDED_R15_LINES,
    normalize_vcc_line_token,
    validate_r15_candidate_line,
)

MANDATORY_SINGLE_LAB_DISCLOSURE = (
    "Arc Institute VCC 2026, 10x Flex CRISPRi (single laboratory and platform; "
    "lines within the same panel share target sets and are not independent target samples)"
)
MANDATORY_SINGLE_LINE_DISCLOSURE = (
    "K_conf = 1 single-line confirmation only; cross-line generalizability is limited."
)


def _normalize_excluded_list(excluded_lines: Sequence[str]) -> List[str]:
    out: List[str] = []
    seen = set()
    for raw in excluded_lines:
        tok = normalize_vcc_line_token(raw)
        if tok in ALLOWED_R15_LINES_SET and tok not in seen:
            seen.add(tok)
            out.append(tok)
    return out


def _norm_verdict(v: Any) -> str:
    return "CONFIRMED" if str(v).strip() == "CONFIRMED" else "NOT_CONFIRMED"


def _eval_conf_set(
    df_lines: pd.DataFrame,
    df_boot: pd.DataFrame,
    summary: Dict[str, Any],
    eval_set: Sequence[str],
    min_valid_boot_frac: float,
) -> Dict[str, Any]:
    expected_sources = sorted(str(x) for x in summary.get("source_lines", []))
    he1 = dict(
        _eval_claim1_he1(
            df_lines=df_lines,
            df_boot=df_boot,
            summary=summary,
            eval_set=list(eval_set),
            expected_sources=expected_sources,
            min_valid_boot_frac=min_valid_boot_frac,
        )
    )
    he2 = dict(_eval_claim2_he2(df_lines=df_lines, df_boot=df_boot, eval_set=list(eval_set)))
    he1["verdict"] = _norm_verdict(he1.get("verdict"))
    he2["verdict"] = _norm_verdict(he2.get("verdict"))
    for hd in (he1, he2):
        subs = dict(hd.get("sub_hypotheses", {}))
        hd["sub_hypotheses"] = {k: {**dict(v), "verdict": _norm_verdict(v.get("verdict"))} for k, v in subs.items()}
    util = _evaluate_util_noninferiority(df_boot=df_boot, eval_set=list(eval_set), noninf_ci_level="ci975")
    return {
        "Group_E_Claim1_Verdict": he1["verdict"],
        "Group_E_Claim2_Verdict": he2["verdict"],
        "H-E1": he1,
        "H-E2": he2,
        "utility_noninferiority": util,
    }


def _extract_single_line_verdict(
    df_lines: pd.DataFrame,
    df_boot: pd.DataFrame,
    summary: Dict[str, Any],
    line: str,
    min_valid_boot_frac: float,
) -> Dict[str, Any]:
    sub_adj = _eval_conf_set(
        df_lines=df_lines,
        df_boot=df_boot,
        summary=summary,
        eval_set=[line],
        min_valid_boot_frac=min_valid_boot_frac,
    )
    s1 = sub_adj["H-E1"]["sub_hypotheses"]
    s2 = sub_adj["H-E2"]["sub_hypotheses"]
    return {
        "line": line,
        "Claim1_H_E1": str(sub_adj["Group_E_Claim1_Verdict"]),
        "Claim2_H_E2": str(sub_adj["Group_E_Claim2_Verdict"]),
        "sub": {
            "H-E1a": str(s1["H-E1a"]["verdict"]),
            "H-E1b": str(s1["H-E1b"]["verdict"]),
            "H-E1c": str(s1["H-E1c"]["verdict"]),
            "H-E2a": str(s2["H-E2a"]["verdict"]),
            "H-E2b": str(s2["H-E2b"]["verdict"]),
        },
    }


def _scope_meta_for_k(k_conf: int, study_complete: bool) -> Dict[str, Any]:
    if k_conf == 0:
        status = "NO_CONFIRMATORY_LINES_DESCRIPTIVE_ONLY"
        phrase = "K_conf = 0 (CONFIRMATORY_ABORTED_K0): zero confirmatory claims tested; descriptive reporting only."
    elif k_conf == 1:
        status = "CONFIRMATORY_SINGLE_LINE_K1"
        phrase = MANDATORY_SINGLE_LINE_DISCLOSURE
    else:
        status = "CONFIRMATORY_MULTI_LINE"
        phrase = f"Evaluated across K_conf = {k_conf} unseen human cell lines from a single laboratory and platform ({MANDATORY_SINGLE_LAB_DISCLOSURE})."
    return {
        "confirmatory_status": status,
        "single_line_limited_generalizability_notice": (k_conf == 1),
        "single_lab_platform_notice": (k_conf >= 2),
        "scope_disclosure": phrase,
        "study_complete": bool(study_complete),
        "statement_label": "final" if study_complete else "interim, more lines pending",
    }


def adjudicate_r15_results(
    df_lines: pd.DataFrame,
    df_boot: pd.DataFrame,
    summary: Dict[str, Any],
    confirmatory_lines: Sequence[str],
    excluded_lines: Sequence[str] = (),
    wave_id: str = "wave1",
    sunset_reached: bool = False,
    min_valid_boot_frac: float = 0.95,
    allow_legacy_lines: bool = False,
) -> Dict[str, Any]:
    """
    Adjudicate one Round 15 release wave (`PREREG_H` v2) with `K_conf_wave in 0..5`.
    """
    conf_norm: List[str] = []
    seen = set()
    for raw in confirmatory_lines:
        if allow_legacy_lines:
            tok = normalize_vcc_line_token(raw)
            if tok in EXCLUDED_R15_LINES:
                raise ValueError("EXCLUDED_LINE: 'jurkat' is strictly excluded.")
            cl = str(raw).strip()
        else:
            cl = validate_r15_candidate_line(raw)
        if cl and cl not in seen:
            seen.add(cl)
            conf_norm.append(cl)
    if not allow_legacy_lines and len(conf_norm) > len(ALLOWED_R15_LINES):
        raise ValueError(f"PREREG_H allows at most {len(ALLOWED_R15_LINES)} confirmatory lines; got {conf_norm}")

    if not df_lines.empty and "target_line" in df_lines.columns:
        for raw_eval in df_lines["target_line"].astype(str).unique():
            tok = normalize_vcc_line_token(raw_eval)
            if tok in EXCLUDED_R15_LINES:
                raise ValueError("EXCLUDED_LINE: 'jurkat' appeared in Round 15 evaluation rows.")
        eval_lines = [str(x) for x in sorted(df_lines["target_line"].astype(str).unique())]
    else:
        eval_lines = []

    secondary_lines = [ln for ln in eval_lines if ln not in seen]
    excl_norm = [ln for ln in _normalize_excluded_list(excluded_lines) if ln not in seen and ln not in set(secondary_lines)]
    resolved_in_wave = list(dict.fromkeys([*conf_norm, *secondary_lines, *excl_norm]))
    pending_lines = [ln for ln in ALLOWED_R15_LINES if ln not in set(resolved_in_wave)]
    study_complete = (len(pending_lines) == 0) or bool(sunset_reached)
    scope = _scope_meta_for_k(len(conf_norm), study_complete=study_complete)
    all_lines_desc = (
        _eval_conf_set(
            df_lines=df_lines,
            df_boot=df_boot,
            summary=summary,
            eval_set=eval_lines,
            min_valid_boot_frac=min_valid_boot_frac,
        )
        if eval_lines
        else {}
    )

    if len(conf_norm) == 0:
        return {
            "protocol": "PREREG_H_v2",
            "wave_id": str(wave_id),
            "k_conf": 0,
            "confirmatory_lines": [],
            "secondary_lines": secondary_lines,
            "excluded_lines": excl_norm,
            "jurkat_excluded": True,
            "resolved_lines": resolved_in_wave,
            "pending_lines": pending_lines,
            "sunset_reached": bool(sunset_reached),
            "Group_E_Claim1_Verdict": "CONFIRMATORY_ABORTED_K0",
            "Group_E_Claim2_Verdict": "CONFIRMATORY_ABORTED_K0",
            "H-E1": {"verdict": "CONFIRMATORY_ABORTED_K0", "pass_count": 0, "required": 0, "sub_hypotheses": {}},
            "H-E2": {"verdict": "CONFIRMATORY_ABORTED_K0", "pass_count": 0, "required": 0, "sub_hypotheses": {}},
            "all_lines_descriptive": all_lines_desc,
            "per_line_verdicts": {},
            "note": "PREREG_H §8.2: K_conf = 0 confirmatory lines; zero confirmatory inference executed.",
            **scope,
        }

    base = _eval_conf_set(
        df_lines=df_lines,
        df_boot=df_boot,
        summary=summary,
        eval_set=conf_norm,
        min_valid_boot_frac=min_valid_boot_frac,
    )
    per_line = {
        ln: _extract_single_line_verdict(df_lines, df_boot, summary, ln, min_valid_boot_frac=min_valid_boot_frac)
        for ln in conf_norm
    }
    return {
        **base,
        "all_lines_descriptive": all_lines_desc,
        "protocol": "PREREG_H_v2",
        "wave_id": str(wave_id),
        "k_conf": len(conf_norm),
        "confirmatory_lines": conf_norm,
        "secondary_lines": secondary_lines,
        "excluded_lines": excl_norm,
        "jurkat_excluded": True,
        "resolved_lines": resolved_in_wave,
        "pending_lines": pending_lines,
        "sunset_reached": bool(sunset_reached),
        "per_line_verdicts": per_line,
        **scope,
    }


def aggregate_r15_wave_adjudications(
    wave_adjs: Sequence[Dict[str, Any]],
    extra_excluded_lines: Sequence[str] = (),
    sunset_reached: bool = False,
) -> Dict[str, Any]:
    """
    Aggregate sequential `PREREG_H §8.4` release waves via strict cross-wave Intersection-Union Test.

    Every qualifying line from every wave enters `confirmatory_lines_total` (`0 <= K_conf_total <= 5`).
    Overall Claim 1 (`H-E1`) and Claim 2 (`H-E2`) pass CONFIRMED iff ALL qualifying lines across
    ALL waves pass unanimously (`K_conf_total / K_conf_total`). Later waves can ONLY maintain or
    revoke an earlier confirmation — NEVER rescue a failed earlier wave.
    """
    if not wave_adjs:
        raise ValueError("Cannot aggregate zero Round 15 wave adjudications.")

    all_conf: List[str] = []
    seen_conf = set()
    all_sec: List[str] = []
    all_excl: List[str] = []
    combined_per_line: Dict[str, Dict[str, Any]] = {}
    any_sunset = bool(sunset_reached)
    trajectory: List[Dict[str, Any]] = []

    for idx, w in enumerate(wave_adjs):
        wid = str(w.get("wave_id") or f"wave_{idx + 1}")
        any_sunset = any_sunset or bool(w.get("sunset_reached", False))
        w_conf = [validate_r15_candidate_line(x) for x in w.get("confirmatory_lines", [])]
        w_plv = w.get("per_line_verdicts", {})
        for ln in w_conf:
            if ln in seen_conf:
                raise ValueError(f"Duplicate confirmatory line '{ln}' across Round 15 release waves.")
            if ln not in w_plv:
                raise ValueError(f"Wave '{wid}' missing per_line_verdicts['{ln}'] required for §8.4 IUT aggregation.")
            seen_conf.add(ln)
            all_conf.append(ln)
            combined_per_line[ln] = dict(w_plv[ln])

        for s in w.get("secondary_lines", []):
            tok = normalize_vcc_line_token(s)
            if tok and tok not in seen_conf and tok not in all_sec:
                all_sec.append(tok)
        for ex in w.get("excluded_lines", []):
            tok = normalize_vcc_line_token(ex)
            if tok in ALLOWED_R15_LINES_SET and tok not in seen_conf and tok not in all_excl:
                all_excl.append(tok)

        cum_v1, cum_v2 = _compute_cross_line_verdict(all_conf, combined_per_line)
        trajectory.append(
            {
                "wave_index": idx + 1,
                "wave_id": wid,
                "wave_confirmatory_lines": w_conf,
                "wave_Claim1": str(w.get("Group_E_Claim1_Verdict")),
                "wave_Claim2": str(w.get("Group_E_Claim2_Verdict")),
                "cumulative_k_conf": len(all_conf),
                "cumulative_Claim1": cum_v1,
                "cumulative_Claim2": cum_v2,
            }
        )

    for ex in _normalize_excluded_list(extra_excluded_lines):
        if ex not in seen_conf and ex not in all_sec and ex not in all_excl:
            all_excl.append(ex)

    all_sec = [ln for ln in all_sec if ln not in seen_conf]
    resolved_set = seen_conf | set(all_sec) | set(all_excl)
    pending_lines = [ln for ln in ALLOWED_R15_LINES if ln not in resolved_set]
    study_complete = (len(pending_lines) == 0) or any_sunset
    k_tot = len(all_conf)
    scope = _scope_meta_for_k(k_tot, study_complete=study_complete)

    if k_tot == 0:
        return {
            "protocol": "PREREG_H_v2",
            "mode": "cross_wave_iut",
            "n_waves": len(wave_adjs),
            "k_conf_total": 0,
            "confirmatory_lines_total": [],
            "secondary_lines_total": all_sec,
            "excluded_lines_total": all_excl,
            "resolved_lines": [ln for ln in ALLOWED_R15_LINES if ln in resolved_set],
            "pending_lines": pending_lines,
            "sunset_reached": any_sunset,
            "Overall_Claim1_Verdict": "CONFIRMATORY_ABORTED_K0",
            "Overall_Claim2_Verdict": "CONFIRMATORY_ABORTED_K0",
            "H-E1": {"verdict": "CONFIRMATORY_ABORTED_K0", "pass_count": 0, "required": 0, "sub_hypotheses": {}},
            "H-E2": {"verdict": "CONFIRMATORY_ABORTED_K0", "pass_count": 0, "required": 0, "sub_hypotheses": {}},
            "per_line_verdicts": {},
            "wave_trajectory": trajectory,
            **scope,
        }

    sub_counts: Dict[str, Dict[str, Any]] = {}
    for sub_k in ("H-E1a", "H-E1b", "H-E1c", "H-E2a", "H-E2b"):
        pc = sum(1 for ln in all_conf if combined_per_line[ln]["sub"][sub_k] in {"CONFIRMED", "PASS"})
        sub_counts[sub_k] = {
            "pass_count": pc,
            "required": k_tot,
            "verdict": "CONFIRMED" if pc == k_tot else "NOT_CONFIRMED",
        }

    c1_pass_lines = sum(1 for ln in all_conf if combined_per_line[ln]["Claim1_H_E1"] == "CONFIRMED")
    c2_pass_lines = sum(1 for ln in all_conf if combined_per_line[ln]["Claim2_H_E2"] == "CONFIRMED")
    c1_ok = all(sub_counts[s]["verdict"] == "CONFIRMED" for s in ("H-E1a", "H-E1b", "H-E1c")) and c1_pass_lines == k_tot
    c2_ok = all(sub_counts[s]["verdict"] == "CONFIRMED" for s in ("H-E2a", "H-E2b")) and c2_pass_lines == k_tot
    overall_c1 = "CONFIRMED" if c1_ok else "NOT_CONFIRMED"
    overall_c2 = "CONFIRMED" if c2_ok else "NOT_CONFIRMED"

    return {
        "protocol": "PREREG_H_v2",
        "mode": "cross_wave_iut",
        "n_waves": len(wave_adjs),
        "k_conf_total": k_tot,
        "confirmatory_lines_total": all_conf,
        "secondary_lines_total": all_sec,
        "excluded_lines_total": all_excl,
        "resolved_lines": [ln for ln in ALLOWED_R15_LINES if ln in resolved_set],
        "pending_lines": pending_lines,
        "sunset_reached": any_sunset,
        "Overall_Claim1_Verdict": overall_c1,
        "Overall_Claim2_Verdict": overall_c2,
        "H-E1": {
            "verdict": overall_c1,
            "pass_count": c1_pass_lines,
            "required": k_tot,
            "sub_hypotheses": {k: sub_counts[k] for k in ("H-E1a", "H-E1b", "H-E1c")},
        },
        "H-E2": {
            "verdict": overall_c2,
            "pass_count": c2_pass_lines,
            "required": k_tot,
            "sub_hypotheses": {k: sub_counts[k] for k in ("H-E2a", "H-E2b")},
        },
        "per_line_verdicts": combined_per_line,
        "wave_trajectory": trajectory,
        **scope,
    }


def _compute_cross_line_verdict(conf_lines: Sequence[str], per_line: Dict[str, Dict[str, Any]]) -> tuple[str, str]:
    if not conf_lines:
        return "CONFIRMATORY_ABORTED_K0", "CONFIRMATORY_ABORTED_K0"
    c1 = "CONFIRMED" if all(per_line[ln]["Claim1_H_E1"] == "CONFIRMED" for ln in conf_lines) else "NOT_CONFIRMED"
    c2 = "CONFIRMED" if all(per_line[ln]["Claim2_H_E2"] == "CONFIRMED" for ln in conf_lines) else "NOT_CONFIRMED"
    return c1, c2


def _fmt_f(v: Any) -> str:
    try:
        fv = float(v)
        return f"{fv:+.4f}" if np.isfinite(fv) else "nan"
    except Exception:
        return "nan"


def render_r15_markdown_report(
    df_lines: pd.DataFrame,
    df_boot: pd.DataFrame,
    summary: Dict[str, Any],
    adj: Dict[str, Any],
    run_tag: str,
) -> str:
    """Render per-wave Round 15 Markdown report with mandatory §8 & §9 disclosures."""
    k_conf = int(adj.get("k_conf", 0))
    if k_conf >= 1:
        base_md = render_r12_markdown_report(df_lines=df_lines, df_boot=df_boot, summary=summary, adj=adj, run_tag=run_tag)
    else:
        base_md = (
            f"# Round 15 `PREREG_H` Report (`{run_tag}`)\n\n"
            f"- **Wave ID**: `{adj.get('wave_id')}` | **Status Label**: `{adj.get('statement_label')}`\n"
            f"- **Confirmatory Claim 1 (`H-E1`)**: **{adj.get('Group_E_Claim1_Verdict')}**\n"
            f"- **Confirmatory Claim 2 (`H-E2`)**: **{adj.get('Group_E_Claim2_Verdict')}**\n"
            f"- **Disclosure**: {adj.get('scope_disclosure')}\n"
        )

    lines_md: List[str] = [
        base_md.rstrip(),
        "",
        "---",
        "",
        "## Round 15 (`PREREG_H` v2) Governance, Multi-Line, & Wave Disclosures",
        "",
        f"- **Protocol**: `PREREG_H` v2 (`seed = {summary.get('bootstrap_seed')}`, `n_boot = {summary.get('n_boot')}`, `n_perm = {summary.get('n_perm')}`)",
        f"- **Release Wave**: `{adj.get('wave_id')}` (`{adj.get('statement_label')}`) | **Pending Candidate Lines**: `{adj.get('pending_lines', [])}` | **Excluded `jurkat`**: `True`",
        f"- **Mandatory Scope Notice (`§8.2`)**: {adj.get('scope_disclosure')}",
        "- **Arc Challenge Curation Disclosure (`§4.2` / `§10.3`)**: Competition challenge panels filter targets by `>= 80%` median knockdown, select the stronger of 2 sgRNAs per gene, and uniformly downsample to `400` cells/perturbation.",
        "- **Prior-Art Notice (`§10.6`)**: Inspect NeurIPS 2026 `OpenReview lnpzb1HvwX` and `CanPert/CASTER` boundaries (`PRIOR_ART_H.md` / `PRIOR_ART_H_CONTROLLER_FAN2026.md`) cleanly.",
    ]

    # SA-H1 (Leave-vcc2025_h1-Out Descriptive Sensitivity)
    sa_h1 = summary.get("sa_h1", {})
    if sa_h1:
        lines_md.extend([
            "",
            "## `SA-H1` (Leave-`vcc2025_h1`-Out Source Sensitivity — `5` Remaining Sources, Strictly Descriptive)",
            "",
            "| Target Line | `n_targets` | `FCT r_orth_cal` [97.5% CI] | `rho_half_ab_orth` | Paired `FCT - Mis_s` Ranges (Min/Max 97.5% CI lo) |",
            "|---|---:|---|---:|---|",
        ])
        for ln in sorted(sa_h1.keys()):
            d = sa_h1[ln]
            fct_d = d.get("fct", {})
            paired = d.get("paired_fct_minus_mismatched", {})
            lo_vals = [float(v.get("ci975_lo", np.nan)) for v in paired.values()] if paired else []
            lo_rng = f"[{min(lo_vals):+.4f}, {max(lo_vals):+.4f}]" if lo_vals and all(np.isfinite(x) for x in lo_vals) else "n/a"
            lines_md.append(
                f"| `{ln}` | {d.get('n_eval_targets', 0)} | "
                f"{_fmt_f(fct_d.get('r_orth_cal'))} [{_fmt_f(fct_d.get('ci975_lo'))}, {_fmt_f(fct_d.get('ci975_hi'))}] | "
                f"{_fmt_f(fct_d.get('rho_ab_orth'))} | {lo_rng} (`{len(paired)}` sources) |"
            )

    # SA-Panel (129 Shared Val-Test Targets vs Remainder)
    panel_meta = summary.get("sa_panel_meta", {})
    sa_panel = summary.get("sa_panel", {})
    lines_md.extend(["", "## `SA-Panel` (`129` Shared Val-Test Targets vs Remainder — Strictly Descriptive)", ""])
    if not panel_meta.get("available", False) or not sa_panel:
        lines_md.append(f"- **Status**: `{panel_meta.get('status', 'PANEL_OVERLAP_LIST_NOT_AVAILABLE')}` (`path={panel_meta.get('path', '')}`).")
    else:
        lines_md.extend([
            "| Target Line | Subgroup | `n_targets` | `rho_half_ab_orth` | `FCT r_orth_cal` [97.5% CI] |",
            "|---|---|---:|---:|---|",
        ])
        for ln in sorted(sa_panel.keys()):
            for sg_name in ("in_val_test_overlap_panel", "panel_unique_remainder"):
                sg = sa_panel[ln].get(sg_name, {})
                fct_s = sg.get("fct", {})
                lines_md.append(
                    f"| `{ln}` | `{sg_name}` | {sg.get('n_targets', 0)} | {_fmt_f(sg.get('rho_half_ab_orth'))} | "
                    f"{_fmt_f(fct_s.get('r_orth_cal'))} [{_fmt_f(fct_s.get('ci975_lo'))}, {_fmt_f(fct_s.get('ci975_hi'))}] |"
                )

    lines_md.append("")
    return "\n".join(lines_md)


def render_r15_overall_markdown_report(overall_adj: Dict[str, Any], report_tag: str = "overall_r15") -> str:
    """Render cross-wave overall `PREREG_H §8.4` IUT Markdown summary report."""
    lines_md: List[str] = [
        f"# Round 15 Cross-Wave Overall IUT Adjudication (`{report_tag}`)",
        "",
        f"- **Statement Status**: **`{overall_adj.get('statement_label')}`** (`study_complete={overall_adj.get('study_complete')}`, `sunset_reached={overall_adj.get('sunset_reached')}`)",
        f"- **Total Qualifying Lines (`K_conf_total`)**: **`{overall_adj.get('k_conf_total')}`** (`{overall_adj.get('confirmatory_lines_total')}`)",
        f"- **Resolved Secondary / Excluded**: secondary=`{overall_adj.get('secondary_lines_total')}`, excluded=`{overall_adj.get('excluded_lines_total')}`",
        f"- **Pending Candidate Lines**: `{overall_adj.get('pending_lines')}`",
        f"- **Overall Claim 1 (`H-E1`)**: **{overall_adj.get('Overall_Claim1_Verdict')}**",
        f"- **Overall Claim 2 (`H-E2`)**: **{overall_adj.get('Overall_Claim2_Verdict')}**",
        f"- **Mandatory Scope Notice (`§8.2`)**: {overall_adj.get('scope_disclosure')}",
        "",
        "## Per-Wave Progression (`§8.4` Monotonic Cross-Wave IUT)",
        "",
        "| Wave # | Wave ID | Wave Qualifying Lines | Wave `H-E1` | Wave `H-E2` | Cumulative `K_conf` | Cumulative `H-E1` | Cumulative `H-E2` |",
        "|---:|---|---|---|---|---:|---|---|",
    ]
    for wt in overall_adj.get("wave_trajectory", []):
        lines_md.append(
            f"| {wt['wave_index']} | `{wt['wave_id']}` | `{wt['wave_confirmatory_lines']}` | "
            f"**{wt['wave_Claim1']}** | **{wt['wave_Claim2']}** | {wt['cumulative_k_conf']} | "
            f"**{wt['cumulative_Claim1']}** | **{wt['cumulative_Claim2']}** |"
        )
    lines_md.append("")
    return "\n".join(lines_md)


def _parse_aggregate_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Round 15 §8.4 Cross-Wave IUT Adjudication Aggregator")
    p.add_argument("--wave-jsons", type=Path, nargs="+", required=True, help="Ordered per-wave adjudicate JSON files.")
    p.add_argument("--excluded-lines", type=str, default="", help="Optional comma-separated extra excluded lines.")
    p.add_argument("--sunset-reached", action="store_true", default=False)
    p.add_argument("--out-json", type=Path, required=True)
    p.add_argument("--out-md", type=Path, default=None)
    p.add_argument("--tag", type=str, default="overall_r15")
    return p.parse_args(argv)


def main() -> None:
    args = _parse_aggregate_args()
    waves = [json.loads(Path(p).read_text(encoding="utf-8")) for p in args.wave_jsons]
    extra_excl = [x.strip() for x in args.excluded_lines.split(",") if x.strip()]
    overall = aggregate_r15_wave_adjudications(
        waves,
        extra_excluded_lines=extra_excl,
        sunset_reached=bool(args.sunset_reached),
    )
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(overall, indent=2) + "\n", encoding="utf-8")
    if args.out_md is not None:
        args.out_md.parent.mkdir(parents=True, exist_ok=True)
        args.out_md.write_text(render_r15_overall_markdown_report(overall, report_tag=args.tag), encoding="utf-8")
    print(json.dumps(overall, indent=2))


if __name__ == "__main__":
    main()
