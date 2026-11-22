#!/usr/bin/env python3
"""
Round 11 (`PREREG_D v2` + `Addendum A`) Confirmatory Adjudicator & Report Generator.

Primary mask: `all_nontarget`. Secondary masks: `resp_top50`, `SA-8` `resp_top50_split`.

Multiplicity (`PREREG_D v2 §4.1`, `Q1`, `Q5`, `Q8`):
  Two confirmatory claims evaluated with Bonferroni `alpha_claim = 0.05 / 2 = 0.025` ->
  two-sided `97.5%` CIs (`[1.25%, 98.75%]`), with conjunctive IUT inside each claim:
  - Claim 1 (`Group_D_Verdict`, unanimous `K_conf / K_conf`, `K_conf < 2` aborts):
    * `H-D1` (`Addendum A1`): `FCT r_orth_cal(all_nontarget)` `97.5% CI_lo > 0`, `valid_frac >= 0.95`, `rho_AB_orth > 0`
    * `H-D2` (`Q5`): `FCT r_orth_cal(all_nontarget)` `97.5% CI_lo > max(mismatched 97.5% CI_hi)` AND `> perm100 p97.5`
    * `H-D3`: Paired `FCT - E1_ConstMed` (`all_nontarget`) `97.5% CI_lo > 0` (+ `valid_frac >= 0.95`)
  - Claim 2 (`H-D4`, independent methodological claim, unanimous `K_conf / K_conf`):
    * `H-D4a`: `E1_ConstMed(resp_top50)` `r_orth_cal` `97.5% CI_lo > 0` (+ `valid_frac >= 0.95`)
    * `H-D4b`: Paired cross-mask `E1(resp_top50) - E1(all_nontarget)` `97.5% CI_lo > 0` (+ `valid_frac >= 0.95`)

Parallel practical utility non-inferiority (`H-D_util`, `H-C3` continuation, `97.5%` CI):
  - Paired `FCT - Transfer_Mean` across all 6 conventions (`0` inferior pairs).
Descriptive `H-D5` (`feng_ipsc`): `FCT` vs `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` under `all_nontarget`.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Sequence

import numpy as np
import pandas as pd

R11_SRC_DIR = Path(__file__).resolve().parent
R10_SRC_DIR = R11_SRC_DIR.parents[1] / "round10" / "src"
for p in (str(R10_SRC_DIR), str(R11_SRC_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from adjudicate_r10 import SIX_CONVENTIONS, _fval  # noqa: E402


def _evaluate_hd4_lines(
    target_lines: Sequence[str],
    df_lines: pd.DataFrame,
    df_boot: pd.DataFrame,
    min_valid_boot_frac: float = 0.95,
    primary_ci_key: str = "ci975",
) -> Dict[str, Any]:
    """
    Evaluate `H-D4` (part A: `E1_ConstMed(resp_top50) > 0`; part B: paired cross-mask
    `E1(resp_top50) - E1(all_nontarget) > 0`) across `target_lines`, recording `ci95`,
    `ci975`, and `ci9833` checks explicitly.
    """
    per_line: Dict[str, Any] = {}
    pass_counts = {"ci95": 0, "ci975": 0, "ci9833": 0}
    e1_rows = (
        {r["target_line"]: r for _, r in df_lines[df_lines["method"] == "E1_ConstMed"].iterrows()}
        if not df_lines.empty
        else {}
    )
    for ln in target_lines:
        row = e1_rows.get(ln)
        cross = (
            df_boot[
                (df_boot["target_line"] == ln)
                & (df_boot["comparison"] == "H-D4_E1_resp_top50_minus_all_nontarget")
            ]
            if not df_boot.empty
            else pd.DataFrame()
        )
        if row is None or cross.empty:
            per_line[ln] = {"pass": False, "reason": "missing_e1_or_cross_mask_row"}
            continue
        cr = cross.iloc[0]
        vfrac_t50 = _fval(row.get("resp_top50_r_orth_cal_valid_frac"))
        vfrac_cross = _fval(cr.get("valid_frac"))
        vfrac_ok = (
            np.isfinite(vfrac_t50)
            and vfrac_t50 >= min_valid_boot_frac
            and (not np.isfinite(vfrac_cross) or vfrac_cross >= min_valid_boot_frac)
        )

        rec: Dict[str, Any] = {
            "e1_resp_top50": _fval(row.get("resp_top50_r_orth_cal")),
            "e1_all_nontarget": _fval(row.get("r_orth_cal")),
            "cross_diff": _fval(cr.get("estimate")),
            "valid_frac_resp_top50": vfrac_t50,
            "valid_frac_cross": vfrac_cross,
        }
        for ck in ("ci95", "ci975", "ci9833"):
            t50_lo = _fval(row.get(f"resp_top50_r_orth_cal_{ck}_lo"))
            t50_hi = _fval(row.get(f"resp_top50_r_orth_cal_{ck}_hi"))
            diff_lo = _fval(cr.get(f"{ck}_lo"))
            diff_hi = _fval(cr.get(f"{ck}_hi"))
            ok_ck = bool(vfrac_ok and np.isfinite(t50_lo) and t50_lo > 0.0 and np.isfinite(diff_lo) and diff_lo > 0.0)
            if ok_ck:
                pass_counts[ck] += 1
            rec[ck] = {
                "e1_resp_top50_lo": t50_lo,
                "e1_resp_top50_hi": t50_hi,
                "cross_diff_lo": diff_lo,
                "cross_diff_hi": diff_hi,
                "part_a_e1_top50_gt_0": bool(np.isfinite(t50_lo) and t50_lo > 0.0),
                "part_b_cross_diff_gt_0": bool(np.isfinite(diff_lo) and diff_lo > 0.0),
                "pass": ok_ck,
            }
        rec["pass"] = rec[primary_ci_key]["pass"]
        per_line[ln] = rec

    n_req = len(target_lines)
    verdicts = {ck: ("CONFIRMED" if (n_req > 0 and pass_counts[ck] == n_req) else "未确认") for ck in pass_counts}
    return {
        "primary_ci_level": primary_ci_key,
        "verdict": verdicts[primary_ci_key],
        "pass_count": int(pass_counts[primary_ci_key]),
        "required": int(n_req),
        "verdicts_by_ci_level": verdicts,
        "pass_counts_by_ci_level": pass_counts,
        "per_line": per_line,
    }


def _evaluate_hd5_feng_ipsc(df_lines: pd.DataFrame, df_boot: pd.DataFrame) -> Dict[str, Any]:
    """
    `PREREG_D §5 H-D5`: Descriptive non-confirmatory comparison on `feng_ipsc`
    (`all_nontarget`) between `FCT` (using non-KOLF own mu_t) and
    `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)`.
    """
    if df_lines.empty or "feng_ipsc" not in set(df_lines["target_line"].tolist()):
        return {"present": False}
    f_sub = df_lines[df_lines["target_line"] == "feng_ipsc"]
    fct_r = f_sub[f_sub["method"] == "FCT"]
    kolf_r = f_sub[f_sub["method"] == "NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)"]
    if fct_r.empty or kolf_r.empty:
        return {"present": True, "status": "missing_fct_or_kolf21j_row"}

    fr = fct_r.iloc[0]
    kr = kolf_r.iloc[0]
    paired_sub = (
        df_boot[
            (df_boot["target_line"] == "feng_ipsc")
            & (df_boot["comparison"] == "H-D5_FCT_minus_Mis_kolf21j_all_nontarget")
        ]
        if not df_boot.empty
        else pd.DataFrame()
    )
    pr = paired_sub.iloc[0] if not paired_sub.empty else None
    fct_lo = _fval(fr.get("r_orth_cal_ci975_lo"))
    kolf_hi = _fval(kr.get("r_orth_cal_ci975_hi"))
    diff_lo = _fval(pr.get("ci975_lo")) if pr is not None else float("nan")
    return {
        "present": True,
        "role": "SECONDARY_NON_CONFIRMATORY",
        "n_eval_targets": int(fr.get("n_eval_targets", 0)),
        "n_genes_eval": int(fr.get("n_genes_eval", 0)),
        "fct_r_orth_cal": _fval(fr.get("r_orth_cal")),
        "fct_ci975_lo": fct_lo,
        "fct_ci975_hi": _fval(fr.get("r_orth_cal_ci975_hi")),
        "mis_kolf21j_r_orth_cal": _fval(kr.get("r_orth_cal")),
        "mis_kolf21j_ci975_lo": _fval(kr.get("r_orth_cal_ci975_lo")),
        "mis_kolf21j_ci975_hi": kolf_hi,
        "paired_fct_minus_mis_kolf21j": _fval(pr.get("estimate")) if pr is not None else float("nan"),
        "paired_ci975_lo": diff_lo,
        "paired_ci975_hi": _fval(pr.get("ci975_hi")) if pr is not None else float("nan"),
        "fct_strictly_exceeds_kolf21j_paired_975": bool(np.isfinite(diff_lo) and diff_lo > 0.0),
        "fct_ci_disjoint_above_kolf21j_975": bool(np.isfinite(fct_lo) and np.isfinite(kolf_hi) and fct_lo > kolf_hi),
        "addendum_a3_note": "Split halves share genetic donors across sibling iPSC clones; noise ceiling reflects clone/channel variability rather than inter-donor variation.",
    }


def adjudicate_r11_results(
    df_lines: pd.DataFrame,
    df_boot: pd.DataFrame,
    summary: Dict[str, Any],
    confirmatory_lines: Sequence[str],
    min_valid_boot_frac: float = 0.95,
    hd4_ci_level: str = "ci975",
    noninf_ci_level: str = "ci975",
) -> Dict[str, Any]:
    conf_lines = list(confirmatory_lines)
    k_conf = len(conf_lines)
    all_lines = sorted(df_lines["target_line"].unique().tolist()) if not df_lines.empty else []
    secondary_lines = [l for l in all_lines if l not in set(conf_lines)]

    out: Dict[str, Any] = {
        "prereg_commit_hash": summary.get("prereg_commit_hash", ""),
        "code_freeze_hash": summary.get("code_freeze_hash", ""),
        "source_mode": summary.get("source_mode", "r9_6"),
        "primary_mask": "all_nontarget",
        "secondary_masks": ["resp_top50", "resp_top50_split"],
        "ci_level_primary": "ci975",
        "bootstrap_seed": summary.get("bootstrap_seed", 20261123),
        "n_boot": summary.get("n_boot", 4000),
        "phase": summary.get("phase", ""),
        "confirmatory_lines": conf_lines,
        "secondary_lines": secondary_lines,
        "k_conf": int(k_conf),
        "k_conf_disclosure": (
            f"K_conf = {k_conf} confirmatory cell lines"
            if k_conf >= 2
            else "Insufficient confirmatory cell lines (K_conf < 2); confirmatory hypothesis testing was not performed (确认性线不足，未做确认性检验)."
        ),
        "cd4t_stim48_in_k_conf": bool("cd4t_stim48" in set(conf_lines)),
        "descriptive_H_D4_all_lines": _evaluate_hd4_lines(
            all_lines, df_lines, df_boot, min_valid_boot_frac=min_valid_boot_frac, primary_ci_key=hd4_ci_level
        ),
        "descriptive_H_D5_feng_ipsc": _evaluate_hd5_feng_ipsc(df_lines, df_boot),
    }

    if k_conf < 2:
        out["confirmatory_status"] = "INSUFFICIENT_LINES"
        out["Group_D_Verdict"] = "NOT_EVALUATED_INSUFFICIENT_LINES"
        out["H-D1"] = {"verdict": "NOT_EVALUATED_INSUFFICIENT_LINES", "k_conf": int(k_conf)}
        out["H-D2"] = {"verdict": "NOT_EVALUATED_INSUFFICIENT_LINES", "k_conf": int(k_conf)}
        out["H-D3"] = {"verdict": "NOT_EVALUATED_INSUFFICIENT_LINES", "k_conf": int(k_conf)}
        out["H-D4"] = {"verdict": "NOT_EVALUATED_INSUFFICIENT_LINES", "k_conf": int(k_conf)}
        out["H-D_util"] = {"verdict": "NOT_EVALUATED_INSUFFICIENT_LINES", "k_conf": int(k_conf)}
        return out

    out["confirmatory_status"] = "EVALUATED"
    sub_df = df_lines[df_lines["target_line"].isin(conf_lines)]
    fct_rows = {r["target_line"]: r for _, r in sub_df[sub_df["method"] == "FCT"].iterrows()}
    perm_map = summary.get("fct_perm100", {})

    # H-D1 (PREREG_D v2 §4.2 & Addendum A1):
    #   FCT r_orth_cal(all_nontarget) 97.5% CI_lo > 0 AND valid_frac >= 0.95 AND rho_ab_orth > 0
    hd1_lines: Dict[str, Any] = {}
    hd1_pass = 0
    for ln in conf_lines:
        row = fct_rows.get(ln)
        if row is None:
            hd1_lines[ln] = {"pass": False, "reason": "missing_fct_row"}
            continue
        est = _fval(row.get("r_orth_cal"))
        lo975 = _fval(row.get("r_orth_cal_ci975_lo"))
        hi975 = _fval(row.get("r_orth_cal_ci975_hi"))
        vfrac = _fval(row.get("r_orth_cal_valid_frac"))
        rho_orth = _fval(row.get("rho_half_ab_orth"))
        ci_pos = bool(np.isfinite(lo975) and lo975 > 0.0)
        vfrac_ok = bool(np.isfinite(vfrac) and vfrac >= min_valid_boot_frac)
        rho_ok = bool(np.isfinite(rho_orth) and rho_orth > 0.0)
        ok = ci_pos and vfrac_ok and rho_ok
        if ok:
            hd1_pass += 1
        hd1_lines[ln] = {
            "estimate": est,
            "ci975_lo": lo975,
            "ci975_hi": hi975,
            "ci9833_lo": _fval(row.get("r_orth_cal_ci9833_lo")),
            "ci9833_hi": _fval(row.get("r_orth_cal_ci9833_hi")),
            "valid_frac": vfrac,
            "rho_half_ab_orth": rho_orth,
            "ci975_lo_gt_0": ci_pos,
            "valid_frac_ge_095": vfrac_ok,
            "rho_ab_orth_gt_0": rho_ok,
            "addendum_a1_ci_gt0_but_low_valid_frac": bool(ci_pos and not vfrac_ok),
            "pass": bool(ok),
        }
    hd1_verdict = "CONFIRMED" if hd1_pass == k_conf else "未确认"
    out["H-D1"] = {
        "verdict": hd1_verdict,
        "pass_count": int(hd1_pass),
        "required": int(k_conf),
        "per_line": hd1_lines,
    }

    # H-D2 (PREREG_D v2 §4.2 & Q5):
    #   FCT(all_nontarget) ci975_lo > max(all mismatched ci975_hi) AND > perm100(all_nontarget) p97.5
    hd2_lines: Dict[str, Any] = {}
    hd2_pass = 0
    for ln in conf_lines:
        row = fct_rows.get(ln)
        if row is None:
            hd2_lines[ln] = {"pass": False, "reason": "missing_fct_row"}
            continue
        fct_lo = _fval(row.get("r_orth_cal_ci975_lo"))
        mis_rows = sub_df[
            (sub_df["target_line"] == ln)
            & (sub_df["method"].str.startswith("NegCtrl_FCT_MismatchedMean("))
        ]
        mis_hi_by_src: Dict[str, float] = {}
        for _, mr in mis_rows.iterrows():
            mis_hi_by_src[str(mr["method"])] = _fval(mr.get("r_orth_cal_ci975_hi"))
        mis_vals = list(mis_hi_by_src.values())
        mis_finite = bool(mis_vals) and all(np.isfinite(v) for v in mis_vals)
        max_mis_hi = float(max(mis_vals)) if mis_finite else float("nan")
        perm_info = perm_map.get(ln, {}).get("all_nontarget", {}).get("r_orth_cal", {})
        perm_p975 = _fval(perm_info.get("p97_5"))
        perm_med = _fval(perm_info.get("median"))
        beats_mis = bool(np.isfinite(fct_lo) and mis_finite and fct_lo > max_mis_hi)
        beats_perm = bool(np.isfinite(fct_lo) and np.isfinite(perm_p975) and fct_lo > perm_p975)
        ok = beats_mis and beats_perm
        if ok:
            hd2_pass += 1
        hd2_lines[ln] = {
            "fct_ci975_lo": fct_lo,
            "n_mismatched_sources": len(mis_hi_by_src),
            "max_mismatched_ci975_hi": max_mis_hi,
            "mismatched_ci975_hi_by_source": mis_hi_by_src,
            "perm100_all_nontarget_p975": perm_p975,
            "perm100_all_nontarget_median": perm_med,
            "beats_mismatched": beats_mis,
            "beats_perm100": beats_perm,
            "pass": bool(ok),
        }
    hd2_verdict = "CONFIRMED" if hd2_pass == k_conf else "未确认"
    out["H-D2"] = {
        "verdict": hd2_verdict,
        "pass_count": int(hd2_pass),
        "required": int(k_conf),
        "per_line": hd2_lines,
    }

    # H-D3 (PREREG_D v2 §4.2): Paired FCT - E1_ConstMed under all_nontarget (`97.5% CI_lo > 0`)
    hd3_lines: Dict[str, Any] = {}
    hd3_pass = 0
    for ln in conf_lines:
        m = (
            df_boot[
                (df_boot["target_line"] == ln)
                & (df_boot["comparison"] == "H-D3_FCT_minus_E1_all_nontarget")
            ]
            if not df_boot.empty
            else pd.DataFrame()
        )
        if m.empty:
            hd3_lines[ln] = {"pass": False, "reason": "missing_hd3_boot_row"}
            continue
        r0 = m.iloc[0]
        est = _fval(r0.get("estimate"))
        lo975 = _fval(r0.get("ci975_lo"))
        hi975 = _fval(r0.get("ci975_hi"))
        vfrac = _fval(r0.get("valid_frac"))
        ok = bool(
            np.isfinite(lo975)
            and lo975 > 0.0
            and (not np.isfinite(vfrac) or vfrac >= min_valid_boot_frac)
        )
        if ok:
            hd3_pass += 1
        hd3_lines[ln] = {
            "estimate": est,
            "ci975_lo": lo975,
            "ci975_hi": hi975,
            "ci9833_lo": _fval(r0.get("ci9833_lo")),
            "ci9833_hi": _fval(r0.get("ci9833_hi")),
            "valid_frac": vfrac,
            "pass": ok,
        }
    hd3_verdict = "CONFIRMED" if hd3_pass == k_conf else "未确认"
    out["H-D3"] = {
        "verdict": hd3_verdict,
        "pass_count": int(hd3_pass),
        "required": int(k_conf),
        "per_line": hd3_lines,
    }

    out["Group_D_Verdict"] = (
        "CONFIRMED"
        if (hd1_verdict == "CONFIRMED" and hd2_verdict == "CONFIRMED" and hd3_verdict == "CONFIRMED")
        else "未确认"
    )

    # Confirmatory Claim 2 (`H-D4`): Unanimous K_conf / K_conf at `97.5%` CI
    out["H-D4"] = _evaluate_hd4_lines(
        conf_lines, df_lines, df_boot, min_valid_boot_frac=min_valid_boot_frac, primary_ci_key=hd4_ci_level
    )

    # H-C3-style Parallel Practical Utility Non-Inferiority vs Transfer_Mean across ALL 6 conventions
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
        for ln in conf_lines:
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
            delta = _fval(r0.get("estimate"))
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
                "delta": delta,
                "ci_lo": lo_v,
                "ci_hi": hi_v,
                "non_inferior": non_inf,
                "superior": sup,
            }
        ni_by_metric[conv] = {
            "n_inferior": int(n_inf),
            "n_superior": int(n_sup),
            "verdict": "CONFIRMED" if n_inf == 0 else "未确认",
            "per_line": line_det,
        }
    out["H-D_util"] = {
        "ci_level": noninf_ci_level,
        "verdict": "CONFIRMED" if len(inferior_pairs) == 0 else "未确认",
        "n_inferior_pairs": len(inferior_pairs),
        "inferior_pairs": inferior_pairs,
        "strict_superiority_count": int(superior_count),
        "total_checks": int(total_checks),
        "by_metric": ni_by_metric,
    }
    return out


def render_r11_markdown_report(
    df_lines: pd.DataFrame,
    df_boot: pd.DataFrame,
    summary: Dict[str, Any],
    adj: Dict[str, Any],
    tag: str,
) -> str:
    lines_md: List[str] = [
        f"# Round 11 Confirmatory Group D (`FCT`, Primary Mask `all_nontarget`) Report (`{tag}`)",
        "",
        f"- **Phase**: `{summary.get('phase')}` | **Source Mode**: `{summary.get('source_mode')}` (`{len(summary.get('source_lines', []))}` sources) | **PREREG_D Commit**: `{adj.get('prereg_commit_hash')}` | **Code Freeze**: `{adj.get('code_freeze_hash')}`",
        f"- **Bootstrap**: `np.random.default_rng({adj.get('bootstrap_seed')})` x `{adj.get('n_boot')}` (Bonferroni 2-claim + IUT two-sided `97.5%` CIs `[1.25%, 98.75%]`)",
        f"- **Cohort Disclosure**: **{adj.get('k_conf_disclosure')}** (`{adj.get('confirmatory_lines')}`)",
        f"- **Secondary Lines**: `{adj.get('secondary_lines')}`",
        f"- **Confirmatory Claim 1 — Primary Group D (`FCT` 全基因线特异信号) Verdict**: **{adj.get('Group_D_Verdict')}**",
        f"- **Confirmatory Claim 2 — Methodological Control (`H-D4` `resp_top50` Selection Artifact)**: **{adj['H-D4'].get('verdict')}** (all lines descriptive: `{adj['descriptive_H_D4_all_lines'].get('verdict')}`)",
        f"- **Parallel Practical Utility (`H-D_util` / `H-C3`-style Non-Inferiority)**: **{adj['H-D_util'].get('verdict')}**",
        "",
    ]

    if adj.get("cd4t_stim48_in_k_conf") or "cd4t_stim48" in set(df_lines["target_line"].tolist() if not df_lines.empty else []):
        lines_md.extend([
            "- **Mandatory `cd4t_stim48` Disclosures (`Q7` & `Addendum A2`)**: `cd4t_stim48` shares primary human CD4+ alpha-beta T lineage, sgRNA library, and Marson lab with Round 10 `cd4t_rest` (weaker biological/laboratory independence). Furthermore, both donors (`CE0008678`, `CE0010866`) share 24 processing channels (`CD4i_R2L25..R2L48`), so residual channel noise slightly elevates split-half noise ceilings (`r_orth_cal` biased conservatively downward).",
            "",
        ])

    if adj.get("confirmatory_status") == "EVALUATED":
        hd1, hd2, hd3, hd4, hni = adj["H-D1"], adj["H-D2"], adj["H-D3"], adj["H-D4"], adj["H-D_util"]
        lines_md.extend([
            "## Confirmatory Hypothesis Adjudication (`97.5%` Two-Sided CIs via IUT + 2-Claim Bonferroni)",
            "",
            "| Hypothesis | Pass / Total | Criterion | Verdict |",
            "|---|---:|---|---|",
            f"| **H-D1** (`FCT` `all_nontarget` 97.5% CI lo > 0, valid >= 0.95, rho_AB_orth > 0) | {hd1['pass_count']}/{hd1['required']} | Unanimous K_conf | **{hd1['verdict']}** |",
            f"| **H-D2** (`FCT` `all_nontarget` 97.5% lo > max mismatched 97.5% hi & perm100 p97.5) | {hd2['pass_count']}/{hd2['required']} | Unanimous K_conf | **{hd2['verdict']}** |",
            f"| **H-D3** (Paired `FCT − E1_ConstMed` under `all_nontarget` 97.5% CI lo > 0) | {hd3['pass_count']}/{hd3['required']} | Unanimous K_conf | **{hd3['verdict']}** |",
            f"| **H-D4** *(Claim 2)* (`E1` `resp_top50 > 0` & `E1(resp_top50) − E1(all_nt) > 0`, `{hd4['primary_ci_level']}`) | {hd4['pass_count']}/{hd4['required']} | Unanimous K_conf | **{hd4['verdict']}** |",
            f"| **H-D_util** (`FCT` vs `Transfer_Mean` across 6 conventions, `{hni['ci_level']}`) | {hni['total_checks'] - hni['n_inferior_pairs']}/{hni['total_checks']} non-inf (`{hni['strict_superiority_count']}` strict sup) | 0 inferior pairs | **{hni['verdict']}** |",
            "",
        ])

    lines_md.extend([
        "## Per-Line Primary Mask (`all_nontarget`, Joint-Bootstrap 97.5% CIs) & Secondary (`resp_top50`)",
        "",
        "| Target Line | Role | Method | N | Primary `r_orth_cal` | Primary 97.5% CI | `resp_top50` `r_orth_cal` | `resp_top50` 97.5% CI | `rho_AB_orth` |",
        "|---|---|---|---:|---:|---|---:|---|---:|",
    ])
    conf_set = set(adj.get("confirmatory_lines", []))
    for _, r in df_lines.iterrows():
        role = "Confirmatory" if r["target_line"] in conf_set else "Secondary"
        lines_md.append(
            f"| {r['target_line']} | {role} | `{r['method']}` | {int(r['n_eval_targets'])} | "
            f"{_fval(r['r_orth_cal']):+.4f} | [{_fval(r['r_orth_cal_ci975_lo']):+.4f}, {_fval(r['r_orth_cal_ci975_hi']):+.4f}] | "
            f"{_fval(r['resp_top50_r_orth_cal']):+.4f} | [{_fval(r['resp_top50_r_orth_cal_ci975_lo']):+.4f}, {_fval(r['resp_top50_r_orth_cal_ci975_hi']):+.4f}] | "
            f"{_fval(r['rho_half_ab_orth']):+.4f} |"
        )

    lines_md.extend(["", "## `H-D3` (`FCT − E1`) & `H-D4` (`E1(resp_top50) − E1(all_nontarget)`) Paired Comparisons", ""])
    if not df_boot.empty:
        sub_b = df_boot[df_boot["comparison"].isin([
            "H-D3_FCT_minus_E1_all_nontarget",
            "diag_FCT_minus_E1_resp_top50",
            "H-D4_E1_resp_top50_minus_all_nontarget",
            "H-D5_FCT_minus_Mis_kolf21j_all_nontarget",
        ])]
        lines_md.append("| Target Line | Comparison | Estimate | 95% CI | 97.5% CI |")
        lines_md.append("|---|---|---:|---|---|")
        for _, br in sub_b.iterrows():
            lines_md.append(
                f"| {br['target_line']} | `{br['comparison']}` | {_fval(br['estimate']):+.4f} | "
                f"[{_fval(br['ci95_lo']):+.4f}, {_fval(br['ci95_hi']):+.4f}] | "
                f"[{_fval(br['ci975_lo']):+.4f}, {_fval(br['ci975_hi']):+.4f}] |"
            )

    lines_md.extend(["", "## `SA-8` Split-Half Out-of-Sample Top-50 (`resp_top50_split`, Single-Half Ceiling `sqrt(rho_AB)`)", ""])
    if not df_lines.empty and "resp_top50_split_r_orth_cal" in df_lines.columns:
        sa8_sub = df_lines[df_lines["method"].isin(["FCT", "E1_ConstMed"])]
        if not sa8_sub.empty:
            lines_md.append("| Target Line | Method | N | `r_orth_raw_split` | `r_orth_cal_split` | 97.5% CI | `Perm100 split` Med / p97.5 |")
            lines_md.append("|---|---|---:|---:|---:|---|---|")
            perm_map = summary.get("fct_perm100", {})
            for _, sr in sa8_sub.iterrows():
                ps = perm_map.get(sr["target_line"], {}).get("resp_top50_split", {}).get("r_orth_cal", {})
                p_str = (
                    f"{_fval(ps.get('median')):+.4f} / {_fval(ps.get('p97_5')):+.4f}"
                    if str(sr["method"]) == "FCT" and ps
                    else "—"
                )
                lines_md.append(
                    f"| {sr['target_line']} | `{sr['method']}` | {int(sr['n_eval_targets'])} | "
                    f"{_fval(sr.get('resp_top50_split_r_orth_raw')):+.4f} | "
                    f"{_fval(sr.get('resp_top50_split_r_orth_cal')):+.4f} | "
                    f"[{_fval(sr.get('resp_top50_split_r_orth_cal_ci975_lo')):+.4f}, {_fval(sr.get('resp_top50_split_r_orth_cal_ci975_hi')):+.4f}] | "
                    f"{p_str} |"
                )

    hd5 = adj.get("descriptive_H_D5_feng_ipsc", {})
    if hd5.get("present"):
        lines_md.extend([
            "",
            "## Descriptive `H-D5` (`feng_ipsc` Non-KOLF vs `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)`, Non-Confirmatory)",
            "",
            f"- **`feng_ipsc`** (N={hd5.get('n_eval_targets')}, G_eval={hd5.get('n_genes_eval')}): "
            f"`FCT` = {_fval(hd5.get('fct_r_orth_cal')):+.4f} [{_fval(hd5.get('fct_ci975_lo')):+.4f}, {_fval(hd5.get('fct_ci975_hi')):+.4f}] vs "
            f"`Mis(kolf21j)` = {_fval(hd5.get('mis_kolf21j_r_orth_cal')):+.4f} [{_fval(hd5.get('mis_kolf21j_ci975_lo')):+.4f}, {_fval(hd5.get('mis_kolf21j_ci975_hi')):+.4f}]; "
            f"paired `FCT − Mis(kolf21j)` = {_fval(hd5.get('paired_fct_minus_mis_kolf21j')):+.4f} [{_fval(hd5.get('paired_ci975_lo')):+.4f}, {_fval(hd5.get('paired_ci975_hi')):+.4f}].",
            f"- *Addendum A3 Notice:* {hd5.get('addendum_a3_note')}",
        ])

    lines_md.extend(["", "## 100-Perm Null (`all_nontarget` Primary, `resp_top50`, & `SA-8 resp_top50_split`)", ""])
    for ln, pinfo in summary.get("fct_perm100", {}).items():
        pa = pinfo.get("all_nontarget", {}).get("r_orth_cal", {})
        pt = pinfo.get("resp_top50", {}).get("r_orth_cal", {})
        ps = pinfo.get("resp_top50_split", {}).get("r_orth_cal", {})
        lines_md.append(
            f"- **{ln}**: `all_nontarget` median={_fval(pa.get('median')):+.4f}, p97.5={_fval(pa.get('p97_5')):+.4f}; "
            f"`resp_top50` median={_fval(pt.get('median')):+.4f}, p97.5={_fval(pt.get('p97_5')):+.4f}; "
            f"`resp_top50_split` median={_fval(ps.get('median')):+.4f}, p97.5={_fval(ps.get('p97_5')):+.4f}"
        )

    return "\n".join(lines_md) + "\n"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Adjudicate Round 11 Group D evaluation outputs.")
    p.add_argument("--tag", required=True)
    p.add_argument("--results-dir", type=Path, default=Path(__file__).resolve().parents[1] / "results")
    p.add_argument("--confirmatory-lines", type=str, default="")
    p.add_argument("--hd4-ci", choices=["ci95", "ci975", "ci9833"], default="ci975")
    p.add_argument("--noninf-ci", choices=["ci95", "ci975", "ci9833"], default="ci975")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    df_lines = pd.read_csv(args.results_dir / f"line_metrics_{args.tag}.csv")
    df_boot = pd.read_csv(args.results_dir / f"bootstrap_ci_{args.tag}.csv")
    summary = json.loads((args.results_dir / f"summary_{args.tag}.json").read_text(encoding="utf-8"))
    if args.confirmatory_lines.strip():
        conf = [x.strip() for x in args.confirmatory_lines.split(",") if x.strip()]
    else:
        conf = list(summary.get("confirmatory_lines", []))

    adj = adjudicate_r11_results(
        df_lines,
        df_boot,
        summary,
        confirmatory_lines=conf,
        hd4_ci_level=args.hd4_ci,
        noninf_ci_level=args.noninf_ci,
    )
    (args.results_dir / f"adjudication_{args.tag}.json").write_text(json.dumps(adj, indent=2), encoding="utf-8")
    (args.results_dir / f"REPORT_{args.tag}.md").write_text(
        render_r11_markdown_report(df_lines, df_boot, summary, adj, args.tag), encoding="utf-8"
    )
    print(f"Wrote adjudication_{args.tag}.json and REPORT_{args.tag}.md", flush=True)


if __name__ == "__main__":
    main()
