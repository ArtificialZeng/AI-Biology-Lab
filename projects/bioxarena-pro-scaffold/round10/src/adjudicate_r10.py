#!/usr/bin/env python3
"""
Round 10 (`PREREG_C`) Confirmatory Adjudicator & Report Generator.

Implements:
  - `P1` K_conf < 2 guard ("Insufficient confirmatory cell lines (K_conf < 2); ... 确认性线不足，未做确认性检验")
  - `K_conf = 2` mandatory explicit cohort-size disclosure
  - Confirmatory Group C rules (`H-C1`, `H-C2`, Bonferroni 97.5% CIs, unanimous K_conf / K_conf)
  - Parallel Practical Utility (`H-C3`, strict AND across all 6 conventions vs `Transfer_Mean`)
  - Secondary reports: `SA-1` (slopes), `SA-2` (GlobalScale), `SA-3` (R9 anchors w/ P7 3D note),
    `SA-4` (R19 3-tier variance decomp), `SA-5` (S6 all_nontarget, S5, SA-5b full-NT), `SA-6 / S7` (`P6`).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List, Sequence

import numpy as np
import pandas as pd

SIX_CONVENTIONS: tuple[str, ...] = (
    "ctrl_disc_raw",
    "ctrl_disc_raw_panel",
    "ctrl_de50_r",
    "cent_disc_raw",
    "cent_disc_raw_panel",
    "cent_de50_r",
)


def _fval(x: Any) -> float:
    try:
        return float(x)
    except Exception:
        return float("nan")


def adjudicate_r10_results(
    df_lines: pd.DataFrame,
    df_boot: pd.DataFrame,
    summary: Dict[str, Any],
    confirmatory_lines: Sequence[str],
    min_valid_boot_frac: float = 0.95,
) -> Dict[str, Any]:
    conf_lines = list(confirmatory_lines)
    k_conf = len(conf_lines)
    all_lines = sorted(df_lines["target_line"].unique().tolist()) if not df_lines.empty else []
    secondary_lines = [l for l in all_lines if l not in set(conf_lines)]

    out: Dict[str, Any] = {
        "prereg_commit_hash": summary.get("prereg_commit_hash", ""),
        "code_freeze_hash": summary.get("code_freeze_hash", ""),
        "bootstrap_seed": summary.get("bootstrap_seed", 20261122),
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
    }

    if k_conf < 2:
        out["confirmatory_status"] = "INSUFFICIENT_LINES"
        out["Group_C_Verdict"] = "NOT_EVALUATED_INSUFFICIENT_LINES"
        out["H-C1"] = {"verdict": "NOT_EVALUATED_INSUFFICIENT_LINES", "k_conf": int(k_conf)}
        out["H-C2"] = {"verdict": "NOT_EVALUATED_INSUFFICIENT_LINES", "k_conf": int(k_conf)}
        out["H-C3"] = {"verdict": "NOT_EVALUATED_INSUFFICIENT_LINES", "k_conf": int(k_conf)}
        return out

    out["confirmatory_status"] = "EVALUATED"
    sub_df = df_lines[df_lines["target_line"].isin(conf_lines)]
    fct_rows = {
        r["target_line"]: r
        for _, r in sub_df[sub_df["method"] == "FCT"].iterrows()
    }
    perm_map = summary.get("fct_perm100", {})

    # H-C1: FCT r_orth_cal 97.5% CI_lo > 0, valid_frac >= 0.95, rho_ab_orth > 0 on K_conf / K_conf lines
    hc1_lines: Dict[str, Any] = {}
    hc1_pass = 0
    for ln in conf_lines:
        row = fct_rows.get(ln)
        if row is None:
            hc1_lines[ln] = {"pass": False, "reason": "missing_fct_row"}
            continue
        est = _fval(row.get("r_orth_cal"))
        lo975 = _fval(row.get("r_orth_cal_ci975_lo"))
        hi975 = _fval(row.get("r_orth_cal_ci975_hi"))
        vfrac = _fval(row.get("r_orth_cal_valid_frac"))
        rho_orth = _fval(row.get("rho_half_ab_orth"))
        ok = (
            np.isfinite(lo975)
            and lo975 > 0.0
            and np.isfinite(vfrac)
            and vfrac >= min_valid_boot_frac
            and np.isfinite(rho_orth)
            and rho_orth > 0.0
        )
        if ok:
            hc1_pass += 1
        hc1_lines[ln] = {
            "estimate": est,
            "ci975_lo": lo975,
            "ci975_hi": hi975,
            "valid_frac": vfrac,
            "rho_half_ab_orth": rho_orth,
            "pass": bool(ok),
        }
    hc1_verdict = "CONFIRMED" if hc1_pass == k_conf else "未确认"
    out["H-C1"] = {
        "verdict": hc1_verdict,
        "pass_count": int(hc1_pass),
        "required": int(k_conf),
        "per_line": hc1_lines,
    }

    # H-C2: FCT r_orth_cal ci975_lo > max(all mismatched mean ci975_hi) AND > perm100 p97.5 on K_conf / K_conf
    hc2_lines: Dict[str, Any] = {}
    hc2_pass = 0
    for ln in conf_lines:
        row = fct_rows.get(ln)
        if row is None:
            hc2_lines[ln] = {"pass": False, "reason": "missing_fct_row"}
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
        perm_info = perm_map.get(ln, {}).get("r_orth_cal", {})
        perm_p975 = _fval(perm_info.get("p97_5"))
        perm_med = _fval(perm_info.get("median"))
        beats_mis = bool(np.isfinite(fct_lo) and mis_finite and fct_lo > max_mis_hi)
        beats_perm = bool(np.isfinite(fct_lo) and np.isfinite(perm_p975) and fct_lo > perm_p975)
        ok = beats_mis and beats_perm
        if ok:
            hc2_pass += 1
        hc2_lines[ln] = {
            "fct_ci975_lo": fct_lo,
            "n_mismatched_sources": len(mis_hi_by_src),
            "max_mismatched_ci975_hi": max_mis_hi,
            "mismatched_ci975_hi_by_source": mis_hi_by_src,
            "perm100_p975": perm_p975,
            "perm100_median": perm_med,
            "beats_mismatched": beats_mis,
            "beats_perm100": beats_perm,
            "pass": bool(ok),
        }
    hc2_verdict = "CONFIRMED" if hc2_pass == k_conf else "未确认"
    out["H-C2"] = {
        "verdict": hc2_verdict,
        "pass_count": int(hc2_pass),
        "required": int(k_conf),
        "per_line": hc2_lines,
    }

    out["Group_C_Verdict"] = "CONFIRMED" if (hc1_verdict == "CONFIRMED" and hc2_verdict == "CONFIRMED") else "未确认"

    # H-C3: Parallel Practical Utility Non-Inferiority vs Transfer_Mean across ALL 6 conventions
    hc3_by_metric: Dict[str, Any] = {}
    inferior_pairs: List[str] = []
    superior_count = 0
    total_checks = 0
    for conv in SIX_CONVENTIONS:
        n_inf = 0
        n_sup = 0
        line_det: Dict[str, Any] = {}
        for ln in conf_lines:
            total_checks += 1
            m = df_boot[
                (df_boot["target_line"] == ln)
                & (df_boot["metric"] == conv)
                & (df_boot["method_a"] == "FCT")
                & (df_boot["method_b"] == "Transfer_Mean")
            ]
            if m.empty:
                n_inf += 1
                inferior_pairs.append(f"{ln}:{conv}(missing)")
                line_det[ln] = {"non_inferior": False, "superior": False, "reason": "missing"}
                continue
            r0 = m.iloc[0]
            delta = _fval(r0.get("delta"))
            lo975 = _fval(r0.get("ci975_lo"))
            hi975 = _fval(r0.get("ci975_hi"))
            non_inf = bool(np.isfinite(hi975) and hi975 >= 0.0)
            sup = bool(np.isfinite(lo975) and lo975 > 0.0)
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
                "ci975_lo": lo975,
                "ci975_hi": hi975,
                "non_inferior": non_inf,
                "superior": sup,
            }
        hc3_by_metric[conv] = {
            "n_inferior": int(n_inf),
            "n_superior": int(n_sup),
            "verdict": "CONFIRMED" if n_inf == 0 else "未确认",
            "per_line": line_det,
        }
    out["H-C3"] = {
        "verdict": "CONFIRMED" if len(inferior_pairs) == 0 else "未确认",
        "n_inferior_pairs": len(inferior_pairs),
        "inferior_pairs": inferior_pairs,
        "strict_superiority_count": int(superior_count),
        "total_checks": int(total_checks),
        "by_metric": hc3_by_metric,
    }
    return out


def render_r10_markdown_report(
    df_lines: pd.DataFrame,
    df_boot: pd.DataFrame,
    summary: Dict[str, Any],
    adj: Dict[str, Any],
    tag: str,
) -> str:
    lines_md: List[str] = [
        f"# Round 10 Confirmatory Group C (`FCT`) Report (`{tag}`)",
        "",
        f"- **Phase**: `{summary.get('phase')}` | **PREREG_C Commit**: `{adj.get('prereg_commit_hash')}` | **Code Freeze**: `{adj.get('code_freeze_hash')}`",
        f"- **Bootstrap**: `np.random.default_rng({adj.get('bootstrap_seed')})` x `{adj.get('n_boot')}` (Bonferroni two-sided `97.5%` CIs)",
        f"- **Cohort Disclosure**: **{adj.get('k_conf_disclosure')}** (`{adj.get('confirmatory_lines')}`)",
        f"- **Secondary Lines**: `{adj.get('secondary_lines')}`",
        f"- **Primary Group C (`FCT`) Verdict**: **{adj.get('Group_C_Verdict')}**",
        f"- **Parallel Practical Utility (`H-C3` Non-Inferiority)**: **{adj['H-C3'].get('verdict')}**",
        "",
    ]

    if adj.get("confirmatory_status") == "EVALUATED":
        hc1, hc2, hc3 = adj["H-C1"], adj["H-C2"], adj["H-C3"]
        lines_md.extend([
            "## Confirmatory Hypothesis Adjudication (97.5% Bonferroni CIs)",
            "",
            "| Hypothesis | Pass / Total | Criterion | Verdict |",
            "|---|---:|---|---|",
            f"| **H-C1** (`FCT` `r_orth_cal` 97.5% CI lo > 0) | {hc1['pass_count']}/{hc1['required']} | Unanimous K_conf | **{hc1['verdict']}** |",
            f"| **H-C2** (`FCT` `r_orth_cal` lo > max 6 mismatched hi & perm100 p97.5) | {hc2['pass_count']}/{hc2['required']} | Unanimous K_conf | **{hc2['verdict']}** |",
            f"| **H-C3** (Non-inferiority vs `Transfer_Mean` across 6 conventions) | {hc3['total_checks'] - hc3['n_inferior_pairs']}/{hc3['total_checks']} non-inf (`{hc3['strict_superiority_count']}` strict sup) | 0 inferior pairs | **{hc3['verdict']}** |",
            "",
        ])

    lines_md.extend([
        "## Per-Line Results (`resp_top50`, Joint-Bootstrap 97.5% CIs)",
        "",
        "| Target Line | Role | Method | N | r_orth_cal | 97.5% CI | Slope OLS | 97.5% CI | r_cal (raw) | rho_AB_orth | pds_cos |",
        "|---|---|---|---:|---:|---|---:|---|---:|---:|---:|",
    ])
    conf_set = set(adj.get("confirmatory_lines", []))
    for _, r in df_lines.iterrows():
        role = "Confirmatory" if r["target_line"] in conf_set else "Secondary"
        lines_md.append(
            f"| {r['target_line']} | {role} | `{r['method']}` | {int(r['n_eval_targets'])} | "
            f"{_fval(r['r_orth_cal']):+.4f} | [{_fval(r['r_orth_cal_ci975_lo']):+.4f}, {_fval(r['r_orth_cal_ci975_hi']):+.4f}] | "
            f"{_fval(r['slope_ols']):+.4f} | [{_fval(r['slope_ols_ci975_lo']):+.4f}, {_fval(r['slope_ols_ci975_hi']):+.4f}] | "
            f"{_fval(r['r_cal']):+.4f} | {_fval(r['rho_half_ab_orth']):+.4f} | {_fval(r['pds_cosine_panel']):.4f} |"
        )

    lines_md.extend(["", "## 100-Perm `NegCtrl_FCT_Perm` (`k=0..99`) & Mismatched Mean Summary", ""])
    for ln, pinfo in summary.get("fct_perm100", {}).items():
        ro = pinfo.get("r_orth_cal", {})
        hc2_l = adj.get("H-C2", {}).get("per_line", {}).get(ln, {})
        max_mis_hi = _fval(hc2_l.get("max_mismatched_ci975_hi"))
        if not np.isfinite(max_mis_hi) and not df_lines.empty:
            mis_sub = df_lines[
                (df_lines["target_line"] == ln)
                & (df_lines["method"].str.startswith("NegCtrl_FCT_MismatchedMean("))
            ]
            if not mis_sub.empty:
                vals = pd.to_numeric(mis_sub["r_orth_cal_ci975_hi"], errors="coerce").dropna()
                if not vals.empty:
                    max_mis_hi = float(vals.max())
        lines_md.append(
            f"- **{ln}**: perm100 `r_orth_cal` median={_fval(ro.get('median')):+.4f}, "
            f"p97.5={_fval(ro.get('p97_5')):+.4f}, p98.75={_fval(ro.get('p98_75')):+.4f}; "
            f"max mismatched hi={max_mis_hi:+.4f}"
        )

    lines_md.extend([
        "",
        "## `SA-3` Round 9 Anchors Notice (`P7`) & `SA-6 / S7` Cell-Count Controls (`P6`)",
        "",
        "- *Benchmark Anchor Note (`P7`):* `Hybrid_ProjShared` follows frozen Round 9 commit `1ca9f60f` 3D `span{1, b_S, a_S}` projection (whereas Round 9 `PREREG_A v4` prose described rank-1 `q_1`).",
    ])
    for ln, s7 in summary.get("s7_cellcount_kd", {}).items():
        sp_nc = s7.get("spearman_r_orth_p_vs_ln_ncells", {})
        sp_kd = s7.get("spearman_r_orth_p_vs_abs_kd", {})
        ols_b = s7.get("ols_adjusted_both", {}).get("adjusted_intercept", {})
        lines_md.append(
            f"- **{ln} S7**: Spearman(`r_orth_p`, ln n_cells)={_fval(sp_nc.get('estimate')):+.4f} "
            f"[{_fval(sp_nc.get('ci975_lo')):+.4f}, {_fval(sp_nc.get('ci975_hi')):+.4f}], "
            f"Spearman(`r_orth_p`, |k_p|)={_fval(sp_kd.get('estimate')):+.4f} "
            f"[{_fval(sp_kd.get('ci975_lo')):+.4f}, {_fval(sp_kd.get('ci975_hi')):+.4f}], "
            f"OLS-adjusted intercept={_fval(ols_b.get('estimate')):+.4f} "
            f"[{_fval(ols_b.get('ci975_lo')):+.4f}, {_fval(ols_b.get('ci975_hi')):+.4f}]"
        )
        for tname, td in s7.get("n_cells_tertiles", {}).items():
            fr = td.get("fct_r_orth_cal", {})
            lines_md.append(
                f"  * `{tname}` (N={td.get('n_targets')}, n_cells=[{td.get('n_cells_min'):.0f}..{td.get('n_cells_max'):.0f}]): "
                f"FCT `r_orth_cal`={_fval(fr.get('estimate')):+.4f} [{_fval(fr.get('ci975_lo')):+.4f}, {_fval(fr.get('ci975_hi')):+.4f}], "
                f"max_mis_hi={_fval(td.get('max_mismatched_ci975_hi')):+.4f}, perm100_p975={_fval(td.get('perm100_p975')):+.4f}, "
                f"HC1={td.get('passes_hc1')}, HC2={td.get('passes_hc2')}"
            )

    return "\n".join(lines_md) + "\n"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Adjudicate Round 10 Group C evaluation outputs.")
    p.add_argument("--tag", required=True, help="Result file tag under round10/results/")
    p.add_argument("--results-dir", type=Path, default=Path(__file__).resolve().parents[1] / "results")
    p.add_argument("--confirmatory-lines", type=str, default="")
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

    adj = adjudicate_r10_results(df_lines, df_boot, summary, confirmatory_lines=conf)
    (args.results_dir / f"adjudication_{args.tag}.json").write_text(
        json.dumps(adj, indent=2), encoding="utf-8"
    )
    (args.results_dir / f"REPORT_{args.tag}.md").write_text(
        render_r10_markdown_report(df_lines, df_boot, summary, adj, args.tag), encoding="utf-8"
    )
    print(f"Wrote adjudication_{args.tag}.json and REPORT_{args.tag}.md", flush=True)


if __name__ == "__main__":
    main()
