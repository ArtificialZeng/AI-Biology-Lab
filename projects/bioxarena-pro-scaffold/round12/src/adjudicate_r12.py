#!/usr/bin/env python3
"""
Round 12 (`PREREG_E`) Automated Confirmatory Adjudicator & Report Generator.

Implements the pre-registered Round 12 rules (`R12_第十二轮行动指令.md` §1–§2):
  - Bonferroni 2-claim + Intersection–Union Test (IUT) using two-sided `97.5%` CIs
    (`[1.25%, 98.75%]`, `ci975`) across all `K_conf` confirmatory lines:
      * Claim 1 (`H-E1` = Target-Line Specificity): `H-E1a ∩ H-E1b ∩ H-E1c`:
          - `H-E1a`: `FCT(all_nontarget)` `ci975_lo > 0` AND `valid_frac >= 0.95`
            AND `rho_half_ab_orth > 0` (`Addendum A1`).
          - `H-E1b`: Paired `FCT - NegCtrl_FCT_MismatchedMean(s)` (`all_nontarget`
            joint bootstrap) `ci975_lo > 0` for EVERY source line `s` in pool.
          - `H-E1c`: `FCT(all_nontarget)` `ci975_lo > Perm100(all_nontarget) Q0.975`.
      * Claim 2 (`H-E2` = Split-Half Top-50 Still Inflates Constant Baseline): `H-E2a ∩ H-E2b`:
          - `H-E2a`: `E1_ConstMed(resp_top50_split)` `ci975_lo > 0`.
          - `H-E2b`: Joint cross-mask `E1_ConstMed(resp_top50_split) - E1_ConstMed(all_nontarget)`
            over common half-valid targets `ci975_lo > 0`.
  - Confirmatory cohort size rule (§1): `K_conf == 1` (`a549` only) IS allowed for
    confirmatory adjudication with mandatory limited-generalizability disclosure;
    `K_conf == 0` aborts confirmatory adjudication (`NOT_EVALUATED_NO_CONFIRMATORY_LINES`).
  - Parallel reporting-only tables: unpaired `H-D2`, `FCT - E1`, `N_5` null cone
    (explicitly flagging by-construction zero/suppressed controls), `resp_top50`,
    `H-D_util`, and non-confirmatory/secondary lines.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Sequence

import numpy as np
import pandas as pd

ROUND12_SRC_DIR = Path(__file__).resolve().parent
REPO_ROOT = ROUND12_SRC_DIR.parents[1]
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
R11_SRC_DIR = REPO_ROOT / "round11" / "src"
for p in (str(R10_SRC_DIR), str(R11_SRC_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from adjudicate_r10 import SIX_CONVENTIONS, _fval  # noqa: E402


def _perm_q(perm_entry: Dict[str, Any], mask_mode: str = "all_nontarget", q_key: str = "p97_5") -> float:
    return _fval(perm_entry.get(mask_mode, {}).get("r_orth_cal", {}).get(q_key))


def _evaluate_util_noninferiority(
    df_boot: pd.DataFrame,
    eval_set: Sequence[str],
    noninf_ci_level: str = "ci975",
) -> Dict[str, Any]:
    """Reporting-only `H-D_util` / `H-C3`-style non-inferiority across `SIX_CONVENTIONS`."""
    ni_by_metric: Dict[str, Any] = {}
    inferior_pairs: List[str] = []
    superior_count = 0
    total_checks = 0
    lo_col = f"{noninf_ci_level}_lo"
    hi_col = f"{noninf_ci_level}_hi"
    for conv in SIX_CONVENTIONS:
        n_inf = 0
        n_sup = 0
        line_det: Dict[str, Any] = {}
        for ln in eval_set:
            total_checks += 1
            m = (
                df_boot[
                    (df_boot["target_line"] == ln)
                    & (df_boot["metric"] == conv)
                    & (df_boot["method_a"] == "FCT")
                    & (df_boot["method_b"] == "Transfer_Mean")
                ]
                if not df_boot.empty
                else pd.DataFrame()
            )
            if m.empty:
                n_inf += 1
                inferior_pairs.append(f"{ln}:{conv}(missing)")
                line_det[ln] = {"non_inferior": False, "superior": False, "reason": "missing"}
                continue
            r0 = m.iloc[0]
            lo_v = _fval(r0.get(lo_col))
            hi_v = _fval(r0.get(hi_col))
            non_inf = bool(np.isfinite(hi_v) and hi_v >= 0.0)
            sup = bool(np.isfinite(lo_v) and lo_v > 0.0)
            if not non_inf:
                n_inf += 1
                inferior_pairs.append(f"{ln}:{conv}")
            if sup:
                n_sup += 1
                superior_count += 1
            line_det[ln] = {
                "mean_fct": _fval(r0.get("mean_a")),
                "mean_tmean": _fval(r0.get("mean_b")),
                "delta": _fval(r0.get("estimate")),
                "ci_lo": lo_v,
                "ci_hi": hi_v,
                "non_inferior": non_inf,
                "superior": sup,
            }
        ni_by_metric[conv] = {
            "n_inferior": int(n_inf),
            "n_superior": int(n_sup),
            "verdict": "CONFIRMED" if eval_set and n_inf == 0 else "未确认",
            "per_line": line_det,
        }
    return {
        "ci_level": noninf_ci_level,
        "verdict": "CONFIRMED" if eval_set and len(inferior_pairs) == 0 else "未确认",
        "n_inferior_pairs": len(inferior_pairs),
        "inferior_pairs": inferior_pairs,
        "strict_superiority_count": int(superior_count),
        "total_checks": int(total_checks),
        "by_metric": ni_by_metric,
    }


def _eval_claim1_he1(
    df_lines: pd.DataFrame,
    df_boot: pd.DataFrame,
    summary: Dict[str, Any],
    eval_set: Sequence[str],
    expected_sources: Sequence[str],
    min_valid_boot_frac: float = 0.95,
) -> Dict[str, Any]:
    """Evaluate `H-E1a`, `H-E1b` (all paired `FCT - Mis_s`), and `H-E1c` across `eval_set`."""
    per_line: Dict[str, Any] = {}
    pa = pb = pc = p_all = 0
    perm_map = summary.get("fct_perm100", {})

    for ln in eval_set:
        sub_l = df_lines[df_lines["target_line"] == ln]
        fct_df = sub_l[sub_l["method"] == "FCT"]
        fr = fct_df.iloc[0] if not fct_df.empty else None

        # H-E1a
        if fr is None:
            est = lo975 = hi975 = vfrac = rho_orth = float("nan")
            ci_pos = vfrac_ok = rho_ok = ok_a = False
        else:
            est = _fval(fr.get("r_orth_cal"))
            lo975 = _fval(fr.get("r_orth_cal_ci975_lo"))
            hi975 = _fval(fr.get("r_orth_cal_ci975_hi"))
            vfrac = _fval(fr.get("r_orth_cal_valid_frac"))
            rho_orth = _fval(fr.get("rho_half_ab_orth"))
            ci_pos = bool(np.isfinite(lo975) and lo975 > 0.0)
            vfrac_ok = bool(np.isfinite(vfrac) and vfrac >= min_valid_boot_frac)
            rho_ok = bool(np.isfinite(rho_orth) and rho_orth > 0.0)
            ok_a = ci_pos and vfrac_ok and rho_ok
        if ok_a:
            pa += 1

        # H-E1b: Paired FCT - Mis_s > 0 for EVERY mismatched source in pool
        sub_b = df_boot[df_boot["target_line"] == ln] if not df_boot.empty else pd.DataFrame()
        paired_rows = (
            sub_b[sub_b["comparison"].str.startswith("H-E1b_FCT_minus_Mis_")]
            if not sub_b.empty
            else pd.DataFrame()
        )
        per_src_paired: Dict[str, Any] = {}
        min_pair_lo = float("inf")
        n_pair_pass = 0
        for _, pr in paired_rows.iterrows():
            mb = str(pr.get("method_b", ""))
            sname = mb[len("NegCtrl_FCT_MismatchedMean(") : -1] if mb.startswith("NegCtrl_FCT_MismatchedMean(") else mb
            plo = _fval(pr.get("ci975_lo"))
            phi = _fval(pr.get("ci975_hi"))
            pest = _fval(pr.get("estimate"))
            psrc_ok = bool(np.isfinite(plo) and plo > 0.0)
            if psrc_ok:
                n_pair_pass += 1
            if np.isfinite(plo) and plo < min_pair_lo:
                min_pair_lo = plo
            per_src_paired[sname] = {
                "estimate": pest,
                "ci975_lo": plo,
                "ci975_hi": phi,
                "valid_frac": _fval(pr.get("valid_frac")),
                "pass": psrc_ok,
            }
        expected_others = [s for s in expected_sources if s != ln]
        ok_b = bool(len(expected_others) > 0 and all(per_src_paired.get(s, {}).get("pass", False) for s in expected_others))
        if ok_b:
            pb += 1

        # Descriptive unpaired H-D2 reference
        mis_rows = sub_l[sub_l["method"].str.startswith("NegCtrl_FCT_MismatchedMean(")]
        max_mis_hi = float(np.nanmax(mis_rows["r_orth_cal_ci975_hi"].to_numpy(dtype=np.float64))) if not mis_rows.empty else float("nan")
        unpaired_ok = bool(np.isfinite(lo975) and np.isfinite(max_mis_hi) and lo975 > max_mis_hi)

        # H-E1c: FCT ci975_lo > perm100(all_nontarget) p97.5
        pq975 = _perm_q(perm_map.get(ln, {}), mask_mode="all_nontarget", q_key="p97_5")
        ok_c = bool(np.isfinite(lo975) and np.isfinite(pq975) and lo975 > pq975)
        if ok_c:
            pc += 1

        line_ok = ok_a and ok_b and ok_c
        if line_ok:
            p_all += 1

        per_line[ln] = {
            "fct_r_orth_cal": est,
            "fct_ci975_lo": lo975,
            "fct_ci975_hi": hi975,
            "valid_frac": vfrac,
            "rho_half_ab_orth": rho_orth,
            "addendum_a1_ci_gt0_but_low_valid_frac": bool(ci_pos and not vfrac_ok),
            "H-E1a_pass": bool(ok_a),
            "H-E1b_pass": bool(ok_b),
            "H-E1b_n_sources_checked": int(len(expected_others)),
            "H-E1b_n_sources_pass": int(n_pair_pass),
            "H-E1b_min_paired_ci975_lo": float(min_pair_lo) if math.isfinite(min_pair_lo) else float("nan"),
            "H-E1b_per_source": per_src_paired,
            "H-E1c_pass": bool(ok_c),
            "perm100_all_nontarget_q975": pq975,
            "descriptive_unpaired_max_mis_ci975_hi": max_mis_hi,
            "descriptive_unpaired_pass": unpaired_ok,
            "pass": bool(line_ok),
        }

    req = len(eval_set)
    return {
        "verdict": ("CONFIRMED" if req > 0 and p_all == req else "未确认"),
        "pass_count": int(p_all),
        "required": int(req),
        "sub_hypotheses": {
            "H-E1a": {"pass_count": int(pa), "required": int(req), "verdict": "CONFIRMED" if req > 0 and pa == req else "未确认"},
            "H-E1b": {"pass_count": int(pb), "required": int(req), "verdict": "CONFIRMED" if req > 0 and pb == req else "未确认"},
            "H-E1c": {"pass_count": int(pc), "required": int(req), "verdict": "CONFIRMED" if req > 0 and pc == req else "未确认"},
        },
        "per_line": per_line,
    }


def _eval_claim2_he2(
    df_lines: pd.DataFrame,
    df_boot: pd.DataFrame,
    eval_set: Sequence[str],
) -> Dict[str, Any]:
    """Evaluate `H-E2a` (`E1(resp_top50_split) > 0`) and `H-E2b` (joint `E1(split) - E1(all) > 0`)."""
    per_line: Dict[str, Any] = {}
    pa = pb = p_all = 0

    for ln in eval_set:
        sub_l = df_lines[df_lines["target_line"] == ln]
        e1_df = sub_l[sub_l["method"] == "E1_ConstMed"]
        er = e1_df.iloc[0] if not e1_df.empty else None

        if er is None:
            e1_split = e1_split_lo = e1_split_hi = split_vfrac = float("nan")
            ok_a = False
        else:
            e1_split = _fval(er.get("resp_top50_split_r_orth_cal"))
            e1_split_lo = _fval(er.get("resp_top50_split_r_orth_cal_ci975_lo"))
            e1_split_hi = _fval(er.get("resp_top50_split_r_orth_cal_ci975_hi"))
            split_vfrac = _fval(er.get("resp_top50_split_valid_frac"))
            ok_a = bool(np.isfinite(e1_split_lo) and e1_split_lo > 0.0)
        if ok_a:
            pa += 1

        sub_b = df_boot[df_boot["target_line"] == ln] if not df_boot.empty else pd.DataFrame()
        heb_df = sub_b[sub_b["comparison"] == "H-E2b_E1_split_minus_all"] if not sub_b.empty else pd.DataFrame()
        br = heb_df.iloc[0] if not heb_df.empty else None
        if br is None:
            diff_est = diff_lo = diff_hi = float("nan")
            n_common = 0
            ok_b = False
        else:
            diff_est = _fval(br.get("estimate"))
            diff_lo = _fval(br.get("ci975_lo"))
            diff_hi = _fval(br.get("ci975_hi"))
            n_common = int(br.get("n_items", 0))
            ok_b = bool(np.isfinite(diff_lo) and diff_lo > 0.0)
        if ok_b:
            pb += 1

        line_ok = ok_a and ok_b
        if line_ok:
            p_all += 1

        per_line[ln] = {
            "e1_split_r_orth_cal": e1_split,
            "e1_split_ci975_lo": e1_split_lo,
            "e1_split_ci975_hi": e1_split_hi,
            "e1_split_valid_frac": split_vfrac,
            "H-E2a_pass": bool(ok_a),
            "e1_split_minus_all_est": diff_est,
            "e1_split_minus_all_ci975_lo": diff_lo,
            "e1_split_minus_all_ci975_hi": diff_hi,
            "n_common_perts": n_common,
            "H-E2b_pass": bool(ok_b),
            "pass": bool(line_ok),
        }

    req = len(eval_set)
    return {
        "verdict": ("CONFIRMED" if req > 0 and p_all == req else "未确认"),
        "pass_count": int(p_all),
        "required": int(req),
        "sub_hypotheses": {
            "H-E2a": {"pass_count": int(pa), "required": int(req), "verdict": "CONFIRMED" if req > 0 and pa == req else "未确认"},
            "H-E2b": {"pass_count": int(pb), "required": int(req), "verdict": "CONFIRMED" if req > 0 and pb == req else "未确认"},
        },
        "per_line": per_line,
    }


def adjudicate_r12_results(
    df_lines: pd.DataFrame,
    df_boot: pd.DataFrame,
    summary: Dict[str, Any],
    confirmatory_lines: Sequence[str],
    min_valid_boot_frac: float = 0.95,
) -> Dict[str, Any]:
    """
    Execute strict `PREREG_E` confirmatory & reporting-only adjudication.
    """
    conf_lines = list(confirmatory_lines)
    k_conf = len(conf_lines)
    all_lines = sorted(df_lines["target_line"].unique().tolist()) if not df_lines.empty else []
    secondary_lines = [l for l in all_lines if l not in set(conf_lines)]
    expected_sources = list(summary.get("source_lines", []))

    desc_he1 = _eval_claim1_he1(df_lines, df_boot, summary, all_lines, expected_sources, min_valid_boot_frac)
    desc_he2 = _eval_claim2_he2(df_lines, df_boot, all_lines)
    desc_util = _evaluate_util_noninferiority(df_boot, all_lines, noninf_ci_level="ci975")

    out: Dict[str, Any] = {
        "schema_version": "r12_prereg_e_v1",
        "primary_mask": "all_nontarget",
        "sa8_mask": "resp_top50_split",
        "primary_ci_level": "ci975",
        "k_conf": int(k_conf),
        "confirmatory_lines": conf_lines,
        "secondary_lines": secondary_lines,
        "limited_single_line_generalizability_notice": bool(k_conf == 1),
        "all_lines_descriptive": {
            "H-E1": desc_he1,
            "H-E2": desc_he2,
            "H-D_util_reporting": desc_util,
        },
    }

    if k_conf == 0:
        out["confirmatory_status"] = "NO_CONFIRMATORY_LINES"
        out["Group_E_Claim1_Verdict"] = "NOT_EVALUATED_NO_CONFIRMATORY_LINES"
        out["Group_E_Claim2_Verdict"] = "NOT_EVALUATED_NO_CONFIRMATORY_LINES"
        return out

    conf_he1 = _eval_claim1_he1(df_lines, df_boot, summary, conf_lines, expected_sources, min_valid_boot_frac)
    conf_he2 = _eval_claim2_he2(df_lines, df_boot, conf_lines)
    conf_util = _evaluate_util_noninferiority(df_boot, conf_lines, noninf_ci_level="ci975")

    out["confirmatory_status"] = (
        "CONFIRMATORY_SINGLE_LINE_K1" if k_conf == 1 else f"CONFIRMATORY_K{k_conf}"
    )
    out["H-E1"] = conf_he1
    out["H-E2"] = conf_he2
    out["H-D_util_reporting"] = conf_util
    out["Group_E_Claim1_Verdict"] = conf_he1["verdict"]
    out["Group_E_Claim2_Verdict"] = conf_he2["verdict"]
    return out


def render_r12_markdown_report(
    df_lines: pd.DataFrame,
    df_boot: pd.DataFrame,
    summary: Dict[str, Any],
    adj: Dict[str, Any],
    run_tag: str,
) -> str:
    """Render comprehensive Round 12 Markdown report separating Confirmatory & Reporting-Only."""
    conf_set = set(adj.get("confirmatory_lines", []))
    k_conf = int(adj.get("k_conf", 0))
    lines_md: List[str] = [
        f"# Round 12 Prospective Group E Report (`{run_tag}`)",
        "",
        (
            f"- **Phase**: `{summary.get('phase')}` | **Source Mode**: `{summary.get('source_mode')}` "
            f"(`{len(summary.get('source_lines', []))}` sources) | "
            f"**PREREG_E Commit**: `{summary.get('prereg_commit_hash', '')}` | "
            f"**Code Freeze**: `{summary.get('code_freeze_hash', '')}`"
        ),
        (
            f"- **Bootstrap**: `np.random.default_rng({summary.get('bootstrap_seed')})` x "
            f"`{summary.get('n_boot')}` (Bonferroni 2-claim + IUT two-sided `97.5%` CIs `[1.25%, 98.75%]`)"
        ),
    ]
    if k_conf == 0:
        lines_md.append("- **Cohort Disclosure**: **K_conf = 0 (count gate failed); descriptive reporting only.**")
    elif k_conf == 1:
        lines_md.append(
            f"- **Cohort Disclosure**: **K_conf = 1 confirmatory line** (`{adj['confirmatory_lines']}`) — "
            "confirmatory inference permitted under `PREREG_E §1`; **mandatory disclosure: single-line generalizability is limited.**"
        )
    else:
        lines_md.append(f"- **Cohort Disclosure**: **K_conf = {k_conf} confirmatory lines** (`{adj['confirmatory_lines']}`)")

    if adj.get("secondary_lines"):
        lines_md.append(f"- **Secondary / Non-Confirmatory Lines (Reporting-Only)**: `{adj['secondary_lines']}`")
    lines_md.append(f"- **Confirmatory Claim 1 (`H-E1` Target-Line Specificity)**: **{adj['Group_E_Claim1_Verdict']}**")
    lines_md.append(f"- **Confirmatory Claim 2 (`H-E2` Persistent `SA-8` Constant-Baseline Inflation)**: **{adj['Group_E_Claim2_Verdict']}**")

    if k_conf >= 1:
        he1 = adj["H-E1"]
        he2 = adj["H-E2"]
        s1 = he1["sub_hypotheses"]
        s2 = he2["sub_hypotheses"]
        lines_md.extend([
            "",
            "## Confirmatory Hypothesis Adjudication (`97.5%` Two-Sided CIs via IUT + 2-Claim Bonferroni)",
            "",
            "| Claim / Hypothesis | Pass / Total | Criterion | Verdict |",
            "|---|---:|---|---|",
            f"| **Claim 1 (`H-E1`)** | {he1['pass_count']}/{he1['required']} | `H-E1a ∩ H-E1b ∩ H-E1c` (Unanimous) | **{he1['verdict']}** |",
            f"| ├─ `H-E1a` (`FCT` all_nt 97.5% lo > 0, valid >= 0.95, rho_orth > 0) | {s1['H-E1a']['pass_count']}/{s1['H-E1a']['required']} | Addendum A1 | **{s1['H-E1a']['verdict']}** |",
            f"| ├─ `H-E1b` (Paired `FCT − Mis_s` 97.5% lo > 0, ALL sources) | {s1['H-E1b']['pass_count']}/{s1['H-E1b']['required']} | Source-wise IUT | **{s1['H-E1b']['verdict']}** |",
            f"| └─ `H-E1c` (`FCT` 97.5% lo > Perm100 `all_nontarget` Q0.975) | {s1['H-E1c']['pass_count']}/{s1['H-E1c']['required']} | Perm null | **{s1['H-E1c']['verdict']}** |",
            f"| **Claim 2 (`H-E2`)** | {he2['pass_count']}/{he2['required']} | `H-E2a ∩ H-E2b` (Unanimous) | **{he2['verdict']}** |",
            f"| ├─ `H-E2a` (`E1_ConstMed` under `resp_top50_split` 97.5% lo > 0) | {s2['H-E2a']['pass_count']}/{s2['H-E2a']['required']} | SA-8 inflation | **{s2['H-E2a']['verdict']}** |",
            f"| └─ `H-E2b` (Joint `E1(resp_top50_split) − E1(all_nontarget)` 97.5% lo > 0) | {s2['H-E2b']['pass_count']}/{s2['H-E2b']['required']} | Common-pert IUT | **{s2['H-E2b']['verdict']}** |",
        ])

    # Paired FCT - Mis_s table (Claim 1 H-E1b)
    paired_sub = df_boot[df_boot["comparison"].str.startswith("H-E1b_FCT_minus_Mis_")] if not df_boot.empty else pd.DataFrame()
    if not paired_sub.empty:
        lines_md.extend([
            "",
            "## `H-E1b`: Paired Specificity Contrasts `FCT − NegCtrl_FCT_MismatchedMean(s)` (`all_nontarget`, 97.5% CI)",
            "",
            "| Target Line | Role | Mismatched Source | Estimate | 97.5% CI | `lo > 0` |",
            "|---|---|---|---:|---|---|",
        ])
        for _, r in paired_sub.iterrows():
            role = "Confirmatory" if r["target_line"] in conf_set else "Secondary"
            mb = str(r.get("method_b", ""))
            src = mb[len("NegCtrl_FCT_MismatchedMean(") : -1] if mb.startswith("NegCtrl_FCT_MismatchedMean(") else mb
            lo = _fval(r.get("ci975_lo"))
            hi = _fval(r.get("ci975_hi"))
            ok = "PASS" if (np.isfinite(lo) and lo > 0.0) else "FAIL"
            lines_md.append(
                f"| {r['target_line']} | {role} | `{src}` | {_fval(r.get('estimate')):+.4f} | [{lo:+.4f}, {hi:+.4f}] | **{ok}** |"
            )

    # SA-8 & H-E2 table
    lines_md.extend([
        "",
        "## `H-E2` / `SA-8`: Split-Half Top-50 (`resp_top50_split`) & Joint Cross-Mask `E1(split) − E1(all)`",
        "",
        "| Target Line | Role | N_common | `E1(all_common)` | `E1(split)` [97.5% CI] | Joint `E1(split) − E1(all)` [97.5% CI] | `FCT(split)` [97.5% CI] | `Perm100 split` Q97.5 |",
        "|---|---|---:|---:|---|---|---|---:|",
    ])
    perm_map = summary.get("fct_perm100", {})
    for ln in sorted(df_lines["target_line"].unique()):
        role = "Confirmatory" if ln in conf_set else "Secondary"
        sub_l = df_lines[df_lines["target_line"] == ln]
        e1r = sub_l[sub_l["method"] == "E1_ConstMed"].iloc[0]
        fctr = sub_l[sub_l["method"] == "FCT"].iloc[0]
        h2b = (
            df_boot[(df_boot["target_line"] == ln) & (df_boot["comparison"] == "H-E2b_E1_split_minus_all")]
            if not df_boot.empty
            else pd.DataFrame()
        )
        r2b = h2b.iloc[0] if not h2b.empty else None
        nc = int(r2b["n_items"]) if r2b is not None else 0
        all_c = _fval(r2b["mean_b"]) if r2b is not None else _fval(e1r.get("r_orth_cal"))
        d_est = _fval(r2b["estimate"]) if r2b is not None else float("nan")
        d_lo = _fval(r2b["ci975_lo"]) if r2b is not None else float("nan")
        d_hi = _fval(r2b["ci975_hi"]) if r2b is not None else float("nan")
        sp_q975 = _perm_q(perm_map.get(ln, {}), mask_mode="resp_top50_split", q_key="p97_5")
        lines_md.append(
            f"| {ln} | {role} | {nc} | {all_c:+.4f} | "
            f"{_fval(e1r.get('resp_top50_split_r_orth_cal')):+.4f} [{_fval(e1r.get('resp_top50_split_r_orth_cal_ci975_lo')):+.4f}, {_fval(e1r.get('resp_top50_split_r_orth_cal_ci975_hi')):+.4f}] | "
            f"{d_est:+.4f} [{d_lo:+.4f}, {d_hi:+.4f}] | "
            f"{_fval(fctr.get('resp_top50_split_r_orth_cal')):+.4f} [{_fval(fctr.get('resp_top50_split_r_orth_cal_ci975_lo')):+.4f}, {_fval(fctr.get('resp_top50_split_r_orth_cal_ci975_hi')):+.4f}] | "
            f"{sp_q975:+.4f} |"
        )

    # Reporting-only Legacy Contrasts: Unpaired H-D2 check & Paired FCT - E1_ConstMed
    desc_he1_pl = adj.get("all_lines_descriptive", {}).get("H-E1", {}).get("per_line", {})
    lines_md.extend([
        "",
        "## Reporting-Only Legacy Contrasts (Unpaired `H-D2` & Paired `FCT − E1_ConstMed`)",
        "",
        "| Target Line | Role | Unpaired `FCT_lo − max(Mis_hi)` | `FCT − E1 (all_nt)` [97.5% CI] | `FCT − E1 (resp_top50)` [97.5% CI] |",
        "|---|---|---:|---|---|",
    ])
    for ln in sorted(df_lines["target_line"].unique()):
        role = "Confirmatory" if ln in conf_set else "Secondary"
        d1 = desc_he1_pl.get(ln, {})
        margin_unpaired = _fval(d1.get("fct_ci975_lo")) - _fval(d1.get("descriptive_unpaired_max_mis_ci975_hi"))
        u_tag = "PASS" if d1.get("descriptive_unpaired_pass") else "FAIL"
        r_all_df = df_boot[(df_boot["target_line"] == ln) & (df_boot["comparison"] == "H-D3_FCT_minus_E1_all_nontarget")] if not df_boot.empty else pd.DataFrame()
        r_t50_df = df_boot[(df_boot["target_line"] == ln) & (df_boot["comparison"] == "diag_FCT_minus_E1_resp_top50")] if not df_boot.empty else pd.DataFrame()
        ra = r_all_df.iloc[0] if not r_all_df.empty else {}
        rt = r_t50_df.iloc[0] if not r_t50_df.empty else {}
        lines_md.append(
            f"| {ln} | {role} | {margin_unpaired:+.4f} ({u_tag}) | "
            f"{_fval(ra.get('estimate')):+.4f} [{_fval(ra.get('ci975_lo')):+.4f}, {_fval(ra.get('ci975_hi')):+.4f}] | "
            f"{_fval(rt.get('estimate')):+.4f} [{_fval(rt.get('ci975_lo')):+.4f}, {_fval(rt.get('ci975_hi')):+.4f}] |"
        )

    # Full Per-Line Mask Table (reporting-only methods incl. naive resp_top50)
    lines_md.extend([
        "",
        "## Reporting-Only Full Line Metrics (`all_nontarget` & Naive `resp_top50`, 97.5% CIs)",
        "",
        "| Target Line | Role | Method | N | `all_nontarget r_orth_cal` | 97.5% CI | Naive `resp_top50 r_orth_cal` | `rho_AB_orth` | `valid_frac` |",
        "|---|---|---|---:|---:|---|---:|---:|---:|",
    ])
    for _, r in df_lines.iterrows():
        role = "Confirmatory" if r["target_line"] in conf_set else "Secondary"
        lines_md.append(
            f"| {r['target_line']} | {role} | `{r['method']}` | {int(r['n_eval_targets'])} | "
            f"{_fval(r['r_orth_cal']):+.4f} | [{_fval(r['r_orth_cal_ci975_lo']):+.4f}, {_fval(r['r_orth_cal_ci975_hi']):+.4f}] | "
            f"{_fval(r['resp_top50_r_orth_cal']):+.4f} | {_fval(r['rho_half_ab_orth']):+.4f} | {_fval(r['r_orth_cal_valid_frac']):.4f} |"
        )

    # Reporting-only FCT - E1 & N_5 Null Cone section
    n5_map = summary.get("n5_cone_metrics", {})
    if n5_map:
        lines_md.extend([
            "",
            "## Reporting-Only `N_5` Null-Cone Analysis (`r_cone_sp`, 97.5% CIs)",
            "",
            "*Note: `Transfer_Mean_GlobalScale` and `Source_Grand_Mean` are **zero by exact mathematical construction** under `N_5`, and `E1_ConstMed`, `E2_ConstP10`, `E2_ConstP90`, `E6_LogFCBar`, `E7_FCBarMinus1` are **suppressed by cone construction** (`R2_N5 >= 0.99`); only `FCT`, `Mis_s`, single-draw `Perm0..99`, `E3`/`E5`, and non-cone target models are empirically informative.*",
            "",
            "| Target Line | Method | `N_5 r_cone_sp` | 97.5% CI | `eta` | `R2_N5` | Structural Role |",
            "|---|---|---:|---|---:|---:|---|",
        ])
        for ln in sorted(n5_map.keys()):
            methods_d = n5_map[ln].get("methods", {})
            for mn, mv in methods_d.items():
                rc = mv.get("r_cone_sp", {})
                lines_md.append(
                    f"| {ln} | `{mn}` | {_fval(rc.get('estimate')):+.4f} | "
                    f"[{_fval(rc.get('ci975_lo')):+.4f}, {_fval(rc.get('ci975_hi')):+.4f}] | "
                    f"{_fval(mv.get('eta')):.4f} | {_fval(mv.get('r2_n5')):.6f} | `{mv.get('construction_status')}` |"
                )

    lines_md.extend(["", "## 100-Perm Empirical Null (`all_nontarget`, `SA-8`, & Single-Draw `N_5`)", ""])
    n5_perm_map = summary.get("n5_perm100", {})
    for ln in sorted(perm_map.keys()):
        pa = perm_map[ln].get("all_nontarget", {}).get("r_orth_cal", {})
        ps = perm_map[ln].get("resp_top50_split", {}).get("r_orth_cal", {})
        pn5 = n5_perm_map.get(ln, {}).get("r_cone_sp", {})
        lines_md.append(
            f"- **{ln}**: `all_nontarget` med={_fval(pa.get('median')):+.4f}, Q97.5={_fval(pa.get('p97_5')):+.4f}; "
            f"`SA-8 split` med={_fval(ps.get('median')):+.4f}, Q97.5={_fval(ps.get('p97_5')):+.4f}; "
            f"`N_5` single-draw med={_fval(pn5.get('median')):+.4f}, med|r|={_fval(pn5.get('median_abs')):+.4f}, "
            f"Q97.5={_fval(pn5.get('p97_5')):+.4f} (med_eta={_fval(n5_perm_map.get(ln, {}).get('med_eta')):.4f})"
        )

    return "\n".join(lines_md) + "\n"
