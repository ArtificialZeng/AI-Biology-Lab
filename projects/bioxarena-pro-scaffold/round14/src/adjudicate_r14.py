#!/usr/bin/env python3
"""
Round 14 (`PREREG_G` / `GSE337804`) Confirmatory Adjudicator & Report Renderer.

Reuses frozen Round 12 (`PREREG_E`) claim adjudication machinery cleanly without modifying
any file under `round10/`, `round11/`, `round12/`, or `round13/`:
  - Enforces Round 14 `K_conf in {0, 1}` scope (`a549_myllia` ONLY).
  - Interprets results strictly as an independent same-cell-line replication on `A549`
    (independent lab/chemistry/guides), NEVER cross-cell-line generalization and NEVER
    pooled with Round 12 `a549_ko`.
  - Renders the pre-registered Item 7 descriptive subgroup breakdown (`FCT` and all 6
    paired `FCT - NegCtrl_FCT_MismatchedMean(s')` contrasts stratified by Round 12
    `a549_ko` overlap vs novel targets, reporting-only).
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Dict, Sequence

import pandas as pd

ROUND14_SRC_DIR = Path(__file__).resolve().parent
REPO_ROOT = ROUND14_SRC_DIR.parents[1]
R12_SRC_DIR = REPO_ROOT / "round12" / "src"
for p in (str(R12_SRC_DIR), str(ROUND14_SRC_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from adjudicate_r12 import (  # noqa: E402
    _eval_claim1_he1,
    _eval_claim2_he2,
    _evaluate_util_noninferiority,
    render_r12_markdown_report,
)


def adjudicate_r14_results(
    df_lines: pd.DataFrame,
    df_boot: pd.DataFrame,
    summary: Dict[str, Any],
    confirmatory_lines: Sequence[str],
    min_valid_boot_frac: float = 0.95,
) -> Dict[str, Any]:
    """
    Execute `PREREG_G` confirmatory adjudication (`K_conf in {0, 1}`; `a549_myllia` ONLY;
    no cross-round pooling with Round 12 or Round 13).
    """
    conf_lines = list(confirmatory_lines)
    if any(str(c).lower() != "a549_myllia" for c in conf_lines):
        raise ValueError(
            f"PREREG_G §2 violation: ONLY 'a549_myllia' can enter Round 14 confirmatory cohort; got {conf_lines}"
        )
    if len(conf_lines) > 1:
        raise ValueError(f"PREREG_G §2 violation: K_conf <= 1 in Round 14; got {conf_lines}")

    k_conf = len(conf_lines)
    all_lines = sorted(df_lines["target_line"].unique().tolist()) if not df_lines.empty else []
    secondary_lines = [ln for ln in all_lines if ln not in set(conf_lines)]
    expected_sources = list(summary.get("source_lines", []))

    desc_he1 = _eval_claim1_he1(df_lines, df_boot, summary, all_lines, expected_sources, min_valid_boot_frac)
    desc_he2 = _eval_claim2_he2(df_lines, df_boot, all_lines)
    desc_util = _evaluate_util_noninferiority(df_boot, all_lines, noninf_ci_level="ci975")

    out: Dict[str, Any] = {
        "schema_version": "r14_prereg_g_v1",
        "parent_protocol": "r12_prereg_e_v2",
        "primary_mask": "all_nontarget",
        "sa8_mask": "resp_top50_split",
        "primary_ci_level": "ci975",
        "k_conf": int(k_conf),
        "confirmatory_lines": conf_lines,
        "secondary_lines": secondary_lines,
        "same_cell_line_independent_dataset_only_notice": True,
        "no_cross_round_pooling_notice": True,
        "all_lines_descriptive": {
            "H-E1": desc_he1,
            "H-E2": desc_he2,
            "H-D_util_reporting": desc_util,
        },
    }

    if k_conf == 0:
        out["confirmatory_status"] = "NO_CONFIRMATORY_LINES_DESCRIPTIVE_ONLY"
        out["Group_E_Claim1_Verdict"] = "NOT_EVALUATED_NO_CONFIRMATORY_LINES"
        out["Group_E_Claim2_Verdict"] = "NOT_EVALUATED_NO_CONFIRMATORY_LINES"
        return out

    conf_he1 = _eval_claim1_he1(df_lines, df_boot, summary, conf_lines, expected_sources, min_valid_boot_frac)
    conf_he2 = _eval_claim2_he2(df_lines, df_boot, conf_lines)
    conf_util = _evaluate_util_noninferiority(df_boot, conf_lines, noninf_ci_level="ci975")

    out["confirmatory_status"] = "CONFIRMATORY_SINGLE_LINE_K1"
    out["H-E1"] = conf_he1
    out["H-E2"] = conf_he2
    out["H-D_util_reporting"] = conf_util
    out["Group_E_Claim1_Verdict"] = conf_he1["verdict"]
    out["Group_E_Claim2_Verdict"] = conf_he2["verdict"]
    return out


def _fmt_est_ci975(d: Dict[str, Any]) -> str:
    est = float(d.get("estimate", float("nan")))
    lo = float(d.get("ci975_lo", float("nan")))
    hi = float(d.get("ci975_hi", float("nan")))
    return f"{est:+.4f} [{lo:+.4f}, {hi:+.4f}]"


def render_r14_markdown_report(
    df_lines: pd.DataFrame,
    df_boot: pd.DataFrame,
    summary: Dict[str, Any],
    adj: Dict[str, Any],
    run_tag: str,
) -> str:
    """
    Render Round 14 (`PREREG_G`) Markdown report reusing `render_r12_markdown_report`
    and appending Item 7 Round 12 `a549_ko` overlap vs novel subgroup descriptive tables.
    """
    base_md = render_r12_markdown_report(df_lines, df_boot, summary, adj, run_tag)
    base_md = base_md.replace(
        f"# Round 12 Prospective Group E Report (`{run_tag}`)",
        f"# Round 14 Prospective Same-Line Independent Replication Report (`PREREG_G` — `{run_tag}`)",
        1,
    )

    lines_extra = [
        "",
        "## Pre-Registered Descriptive Subgroup Analysis (`PREREG_G` Item 7 — Reporting-Only)",
        "",
        "*Descriptive-only stratification of evaluated targets by overlap with Round 12 `a549_ko` (`in_r12_a549_overlap` vs `novel_not_in_r12_a549`); zero confirmatory adjudication weight.*",
        "",
        "| Target Line | Subgroup | N_targets | `FCT` (`all_nontarget` `r_orth_cal`) [97.5% CI] | `rho_AB_orth` | `valid_frac` | Min Paired `FCT − Mis(s')` CI97.5_lo |",
        "|---|---|---:|---|---:|---:|---:|",
    ]
    subgroups_by_line = summary.get("overlap_subgroups", {})
    for tgt_ln in sorted(subgroups_by_line.keys()):
        line_sg = subgroups_by_line[tgt_ln]
        for sg_key, sg_label in (
            ("in_r12_a549_overlap", "In Round 12 `a549_ko` Overlap"),
            ("novel_not_in_r12_a549", "Novel (Not Tested in Round 12 `a549_ko`)"),
        ):
            sg = line_sg.get(sg_key, {})
            n_t = int(sg.get("n_targets", 0))
            if n_t == 0:
                lines_extra.append(f"| {tgt_ln} | {sg_label} | 0 | N/A | N/A | N/A | N/A |")
                continue
            fct_d = sg.get("fct", {})
            paired = sg.get("paired_fct_minus_mismatched", {})
            min_lo = min((float(v.get("ci975_lo", float("nan"))) for v in paired.values()), default=float("nan"))
            rho_v = float(sg.get("rho_half_ab_orth", float("nan")))
            vf = float(fct_d.get("valid_frac", float("nan")))
            lines_extra.append(
                f"| {tgt_ln} | {sg_label} | {n_t} | {_fmt_est_ci975(fct_d)} | {rho_v:+.4f} | {vf:.4f} | {min_lo:+.4f} |"
            )

    lines_extra.extend(
        [
            "",
            "| Target Line | Subgroup | N_targets | Mismatched Source | Paired `FCT − Mis(s')` (`all_nontarget`) [97.5% CI] | `valid_frac` |",
            "|---|---|---:|---|---|---:|",
        ]
    )
    for tgt_ln in sorted(subgroups_by_line.keys()):
        line_sg = subgroups_by_line[tgt_ln]
        for sg_key, sg_label in (
            ("in_r12_a549_overlap", "Overlap (`R12 a549_ko`)"),
            ("novel_not_in_r12_a549", "Novel (Not in `R12 a549_ko`)"),
        ):
            sg = line_sg.get(sg_key, {})
            n_t = int(sg.get("n_targets", 0))
            paired = sg.get("paired_fct_minus_mismatched", {})
            for src_name in sorted(paired.keys()):
                pd_row = paired[src_name]
                vf = float(pd_row.get("valid_frac", float("nan")))
                lines_extra.append(
                    f"| {tgt_ln} | {sg_label} | {n_t} | `{src_name}` | {_fmt_est_ci975(pd_row)} | {vf:.4f} |"
                )

    lines_extra.extend(
        [
            "",
            "## Round 14 (`PREREG_G`) Governance & Interpretation Rules",
            "",
            "- **Protocol Inheritance**: Claims 1–2 (`H-E1`, `H-E2`), `r_orth_cal`, `all_nontarget` / `resp_top50_split`, `r9_6`, 2-claim Bonferroni + IUT (`97.5%` CIs), `4,000` bootstraps, and `100` permutations are inherited identically from locked `PREREG_E`; seed=`20261224`.",
            "- **Same-Cell-Line Independent-Dataset Scope ONLY**: Round 14 tests `a549_myllia` (`GSE337804`, Adduri et al. *Cell* 2026 / Myllia) independently under `K_conf <= 1` (`ch_1 -> Half A`, `ch_2 -> Half B`). Confirmation constitutes replication across two independent A549 datasets in two separately pre-registered prospective rounds (NEVER pooled with Round 12 `a549_ko`, and NEVER claimed as cross-cell-line generalization).",
            "- **Chemistry & Technical Context**: 10x `5' HT v2` + Jumpcode `DepleteX`, `SpCas9` CRISPRn (`3` guides/gene), `20` `non-targeting` NTCs (`ctrl_T` `gene_desert` excluded), RNA QC `umi >= 300 & ngenes >= 200`, `max_pct_mt=None`.",
        ]
    )
    return base_md.rstrip("\n") + "\n" + "\n".join(lines_extra) + "\n"


__all__ = ["adjudicate_r14_results", "render_r14_markdown_report"]
