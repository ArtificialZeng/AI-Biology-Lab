#!/usr/bin/env python3
"""
Round 13 (`PREREG_F`) Automated Confirmatory Adjudicator & Report Generator.

Identical hypotheses, thresholds, masks, IUT, 2-claim Bonferroni (`97.5%` CIs), and
reporting-only analyses as `PREREG_E` (`round12/src/adjudicate_r12.py`), adapted solely
to enforce Round 13 §2 confirmatory cohort scoping (`mac_ctrl` ONLY; `mac_lps` strictly
descriptive) without modifying any frozen `round12/` file.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Dict, Sequence

import pandas as pd

ROUND13_SRC_DIR = Path(__file__).resolve().parent
REPO_ROOT = ROUND13_SRC_DIR.parents[1]
R12_SRC_DIR = REPO_ROOT / "round12" / "src"
for p in (str(R12_SRC_DIR), str(ROUND13_SRC_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from adjudicate_r12 import (  # noqa: E402
    _eval_claim1_he1,
    _eval_claim2_he2,
    _evaluate_util_noninferiority,
    render_r12_markdown_report,
)


def adjudicate_r13_results(
    df_lines: pd.DataFrame,
    df_boot: pd.DataFrame,
    summary: Dict[str, Any],
    confirmatory_lines: Sequence[str],
    min_valid_boot_frac: float = 0.95,
) -> Dict[str, Any]:
    """
    Execute `PREREG_F` confirmatory adjudication (`K_conf in {0, 1}`; `mac_ctrl` ONLY;
    `mac_lps` descriptive only; no cross-round pooling with Round 12).
    """
    conf_lines = list(confirmatory_lines)
    if any(str(c).lower() != "mac_ctrl" for c in conf_lines):
        raise ValueError(
            f"PREREG_F §2 violation: ONLY 'mac_ctrl' can enter Round 13 confirmatory cohort; got {conf_lines}"
        )
    if len(conf_lines) > 1:
        raise ValueError(f"PREREG_F §2 violation: K_conf <= 1 in Round 13; got {conf_lines}")

    k_conf = len(conf_lines)
    all_lines = sorted(df_lines["target_line"].unique().tolist()) if not df_lines.empty else []
    secondary_lines = [l for l in all_lines if l not in set(conf_lines)]
    expected_sources = list(summary.get("source_lines", []))

    desc_he1 = _eval_claim1_he1(df_lines, df_boot, summary, all_lines, expected_sources, min_valid_boot_frac)
    desc_he2 = _eval_claim2_he2(df_lines, df_boot, all_lines)
    desc_util = _evaluate_util_noninferiority(df_boot, all_lines, noninf_ci_level="ci975")

    out: Dict[str, Any] = {
        "schema_version": "r13_prereg_f_v1",
        "parent_protocol": "r12_prereg_e_v2",
        "primary_mask": "all_nontarget",
        "sa8_mask": "resp_top50_split",
        "primary_ci_level": "ci975",
        "k_conf": int(k_conf),
        "confirmatory_lines": conf_lines,
        "secondary_lines": secondary_lines,
        "limited_single_line_generalizability_notice": bool(k_conf == 1),
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


def render_r13_markdown_report(
    df_lines: pd.DataFrame,
    df_boot: pd.DataFrame,
    summary: Dict[str, Any],
    adj: Dict[str, Any],
    run_tag: str,
) -> str:
    """Render Round 13 (`PREREG_F`) Markdown report reusing `render_r12_markdown_report`."""
    base_md = render_r12_markdown_report(df_lines, df_boot, summary, adj, run_tag)
    base_md = base_md.replace(
        f"# Round 12 Prospective Group E Report (`{run_tag}`)",
        f"# Round 13 Prospective Replication Report (`PREREG_F` — `{run_tag}`)",
        1,
    )
    addendum = [
        "",
        "## Round 13 (`PREREG_F`) Governance & Interpretation Rules",
        "",
        "- **Protocol Inheritance**: Claims 1–2 (`H-E1`, `H-E2`), `r_orth_cal`, `all_nontarget` / `resp_top50_split`, `r9_6`, 2-claim Bonferroni + IUT (`97.5%` CIs), `4,000` bootstraps, and `100` permutations are inherited identically from locked `PREREG_E`; seed=`20261223`.",
        "- **Separate-Round Independent Scope**: Round 13 has `K_conf <= 1` (`mac_ctrl` ONLY; `mac_lps` descriptive only). Confirmation in Round 13 MUST be reported strictly as independent prospective replication across two separately pre-registered rounds/lines without pooling Round 12 and Round 13 data; failure on any sub-criterion in `mac_ctrl` MUST be reported verbatim as non-replication on the second line.",
        "- **Technical Boundaries**: Single donor (`RG5581`), 10x GEM-X Flex probe panel (`~18k` targeted probes, no mitochondrial cutoff), channel parity (`L01/L03` vs `L02/L04`).",
    ]
    return base_md.rstrip("\n") + "\n" + "\n".join(addendum) + "\n"


__all__ = ["adjudicate_r13_results", "render_r13_markdown_report"]
