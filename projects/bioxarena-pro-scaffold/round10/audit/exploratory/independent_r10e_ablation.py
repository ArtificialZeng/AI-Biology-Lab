#!/usr/bin/env python3
"""
Independent Round 10-E Post-Hoc Baseline Ablation Audit (A-T3).

Post-hoc exploratory check (NOT pre-registered confirmatory; every report table
MUST carry explicit "事后探索性（post hoc）" labeling):
  - Independently implements:
      * FCT: log1p(max(mu_t, 0) * FC_bar) - log1p(max(mu_t, 0))
      * E1:  Constant gene-invariant baseline c = median(max(mu_t, 0)) inside FCT
      * E6:  ln(FC_bar)                         (constant baseline c -> inf limit)
      * E7:  FC_bar - 1                         (constant baseline c -> 0+ limit)
      * FCT - E1: paired difference on the exact same 4,000 joint bootstrap draws
  - Evaluates 5 target lines (`dld1`, `cd4t_rest` against 6 R9 sources;
    `nadig_hepg2`, `vcc2025_h1`, `nadig_jurkat` LOLO against the other 5 R9 sources)
    under both `resp_top50` and `all_nontarget` (S6) masks (`n_boot=4000`, `seed=20261122`,
    two-sided Bonferroni 97.5% CIs `[1.25%, 98.75%]`).
  - Strictly imports frozen functions from `round10/audit/` without modifying any
    file listed in `FREEZE_MANIFEST_R10.txt`.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

AUDIT_DIR = Path(__file__).resolve().parents[1]
if str(AUDIT_DIR) not in sys.path:
    sys.path.insert(0, str(AUDIT_DIR))

from independent_r10_loader import (  # noqa: E402
    R9_LINE_NAMES,
    R9_LINE_SET,
    R9_LINES_DIR,
    R10_LINES_DIR,
    CanonicalHGNC,
    enforce_r10_gate,
    load_line_dataset,
    load_r9_common_genes,
)
from independent_r10_metrics import (  # noqa: E402
    evaluate_joint_residuals,
    evaluate_paired_cal_diff,
)
from independent_r10_predictors import (  # noqa: E402
    AlignedCohort,
    align_source_and_target,
    compute_source_fc_bar,
    fct_from_fc_bar,
)

DEFAULT_AUDIT_LINES: Tuple[str, ...] = (
    "dld1",
    "cd4t_rest",
    "nadig_hepg2",
    "vcc2025_h1",
    "nadig_jurkat",
)
EVAL_MASKS: Tuple[str, ...] = ("resp_top50", "all_nontarget")
PREDICTOR_ORDER: Tuple[str, ...] = ("FCT", "E1", "E6", "E7")


def build_e1_constant_median_pred(
    fc_bar: np.ndarray,
    tgt_ctrl_mean_libnorm: np.ndarray,
) -> Tuple[np.ndarray, float]:
    """
    E1: Replace gene-specific target control mean mu_t with a single scalar median constant
    c_med = median_g(max(mu_{t,g}, 0)) across all genes, eliminating gene-specific baseline info.
    """
    mu_safe = np.maximum(np.asarray(tgt_ctrl_mean_libnorm, dtype=np.float64), 0.0)
    c_med = float(np.median(mu_safe))
    const_mu = np.full(mu_safe.shape, c_med, dtype=np.float64)
    return fct_from_fc_bar(fc_bar, const_mu), c_med


def build_e6_log_fc_pred(fc_bar: np.ndarray, floor: float = 1e-300) -> np.ndarray:
    """
    E6: Limit of FCT(c) as constant baseline c -> +infinity:
      lim_{c->inf} [log1p(c * FC_bar) - log1p(c)] = ln(FC_bar).
    """
    return np.log(np.maximum(np.asarray(fc_bar, dtype=np.float64), floor))


def build_e7_linear_fc_pred(fc_bar: np.ndarray) -> np.ndarray:
    """
    E7: Zero-baseline first-order limit of FCT(c) as c -> 0+ (scale-invariant in correlation):
      lim_{c->0+} [log1p(c * FC_bar) - log1p(c)] / c = FC_bar - 1.
    """
    return np.asarray(fc_bar, dtype=np.float64) - 1.0


def build_ablation_predictors(
    cohort: AlignedCohort,
) -> Tuple[Dict[str, np.ndarray], np.ndarray, float]:
    """Independently compute FCT, E1, E6, E7, and Transfer_Mean for one aligned cohort."""
    fc_bar, t_mean, _ = compute_source_fc_bar(
        eval_targets=cohort.eval_targets,
        source_names=cohort.source_names,
        source_deltas_by_line=cohort.source_deltas_by_line,
        source_ctrl_mean_libnorm=cohort.source_ctrl_mean_libnorm,
    )
    pred_fct = fct_from_fc_bar(fc_bar, cohort.tgt_ctrl_mean_libnorm)
    pred_e1, c_med = build_e1_constant_median_pred(fc_bar, cohort.tgt_ctrl_mean_libnorm)
    pred_e6 = build_e6_log_fc_pred(fc_bar)
    pred_e7 = build_e7_linear_fc_pred(fc_bar)
    return (
        {
            "FCT": pred_fct,
            "E1": pred_e1,
            "E6": pred_e6,
            "E7": pred_e7,
        },
        t_mean,
        c_med,
    )


def evaluate_line_ablation(
    cohort: AlignedCohort,
    masks: Sequence[str] = EVAL_MASKS,
    n_boot: int = 4000,
    seed: int = 20261122,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """
    Run joint bootstrap evaluation (`FCT`, `E1`, `E6`, `E7`) + paired `FCT - E1`
    difference for each requested mask on one target line.
    """
    preds, t_mean, c_med = build_ablation_predictors(cohort)
    metric_rows: List[Dict[str, Any]] = []
    pair_rows: List[Dict[str, Any]] = []

    cohort_tag = "R9_LOLO" if cohort.target_line in R9_LINE_SET else "R10_Holdout"

    for mask in masks:
        fct_est = float("nan")
        per_pred: Dict[str, Dict[str, Any]] = {}
        for p_name in PREDICTOR_ORDER:
            res = evaluate_joint_residuals(
                delta_obs=cohort.delta_obs,
                pred=preds[p_name],
                t_mean=t_mean,
                target_cols=cohort.target_cols,
                mode=mask,
                n_boot=n_boot,
                seed=seed,
                half_a_obs=cohort.half_a_obs,
                half_b_obs=cohort.half_b_obs,
                half_valid=cohort.half_valid,
            )
            per_pred[p_name] = res
            if p_name == "FCT":
                fct_est = float(res.get("r_orth_cal_est", float("nan")))

        for p_name in PREDICTOR_ORDER:
            res = per_pred[p_name]
            r_est = float(res.get("r_orth_cal_est", float("nan")))
            ratio_vs_fct = (
                r_est / fct_est
                if np.isfinite(r_est) and np.isfinite(fct_est) and abs(fct_est) > 1e-12
                else float("nan")
            )
            metric_rows.append(
                {
                    "analysis_type": "post_hoc_exploratory",
                    "cohort": cohort_tag,
                    "target_line": cohort.target_line,
                    "mask": mask,
                    "predictor": p_name,
                    "n_targets": int(res.get("n_eval_perts", 0)),
                    "n_eval_genes": len(cohort.eval_genes),
                    "e1_const_median": c_med,
                    "r_orth_cal": r_est,
                    "r_orth_cal_ci95_lo": res.get("r_orth_cal_ci95_lo"),
                    "r_orth_cal_ci95_hi": res.get("r_orth_cal_ci95_hi"),
                    "r_orth_cal_ci975_lo": res.get("r_orth_cal_ci975_lo"),
                    "r_orth_cal_ci975_hi": res.get("r_orth_cal_ci975_hi"),
                    "r_orth_raw": res.get("r_orth_raw_est"),
                    "r_cal": res.get("r_cal_est"),
                    "slope_ols": res.get("slope_ols_est"),
                    "rho_ab_orth": res.get("rho_ab_orth"),
                    "noise_ceiling_r_orth": res.get("noise_ceiling_r_orth"),
                    "ratio_vs_fct": ratio_vs_fct,
                }
            )

        paired = evaluate_paired_cal_diff(
            delta_obs=cohort.delta_obs,
            pred_a=preds["FCT"],
            pred_b=preds["E1"],
            t_mean=t_mean,
            target_cols=cohort.target_cols,
            mode=mask,
            n_boot=n_boot,
            seed=seed,
            half_a_obs=cohort.half_a_obs,
            half_b_obs=cohort.half_b_obs,
            half_valid=cohort.half_valid,
            ortho=True,
        )
        pair_rows.append(
            {
                "analysis_type": "post_hoc_exploratory",
                "cohort": cohort_tag,
                "target_line": cohort.target_line,
                "mask": mask,
                "comparison": "FCT_minus_E1",
                "n_targets": len(cohort.eval_targets),
                "e1_const_median": c_med,
                "fct_r_orth_cal": per_pred["FCT"].get("r_orth_cal_est"),
                "e1_r_orth_cal": per_pred["E1"].get("r_orth_cal_est"),
                "e1_pct_of_fct": (
                    100.0
                    * float(per_pred["E1"].get("r_orth_cal_est", float("nan")))
                    / float(per_pred["FCT"].get("r_orth_cal_est", float("nan")))
                    if abs(float(per_pred["FCT"].get("r_orth_cal_est", 0.0))) > 1e-12
                    else float("nan")
                ),
                "diff_est": paired.get("diff_est"),
                "diff_ci95_lo": paired.get("diff_ci95_lo"),
                "diff_ci95_hi": paired.get("diff_ci95_hi"),
                "diff_ci975_lo": paired.get("diff_ci975_lo"),
                "diff_ci975_hi": paired.get("diff_ci975_hi"),
                "diff_ci975_excludes_zero": bool(
                    np.isfinite(float(paired.get("diff_ci975_lo", float("nan"))))
                    and np.isfinite(float(paired.get("diff_ci975_hi", float("nan"))))
                    and (
                        float(paired["diff_ci975_lo"]) > 0.0
                        or float(paired["diff_ci975_hi"]) < 0.0
                    )
                ),
            }
        )

    return metric_rows, pair_rows


def _fmt(val: Any, digits: int = 4) -> str:
    try:
        fv = float(val)
        return f"{fv:+.{digits}f}" if np.isfinite(fv) else "nan"
    except (TypeError, ValueError):
        return str(val)


def render_posthoc_markdown(
    metric_rows: Sequence[Dict[str, Any]],
    pair_rows: Sequence[Dict[str, Any]],
    n_boot: int,
    seed: int,
) -> str:
    """Format post-hoc report with explicit '事后探索性（post hoc）' labels on every table."""
    lines: List[str] = [
        "# Round 10-E Independent Post-Hoc Baseline Ablation Audit (`A-T3`)",
        "",
        "> **重要声明（事后探索性分析 / Post Hoc Exploratory）：** 本文全部结果均为事后探索性消融审计（`R10-E`），非第十轮预注册确认性检验；主要科学判定仍以 `PREREG_C` 确认性结论为准（`H-C1 2/2`, `H-C2 0/2 -> 未确认`, `H-C3 2/2 PASS`）。",
        "",
        f"- **Parameters:** `n_boot = {n_boot}`, `seed = {seed}`, two-sided Bonferroni `97.5%` CIs (`[1.25%, 98.75%]`), zero frozen manifest edits.",
        "- **Predictors Verified:** `FCT` (real target baseline $\\mu_t$), `E1` (constant gene-invariant median $c = \\operatorname{median}(\\mu_t)$), `E6` ($\\ln(\\overline{\\mathrm{FC}})$, $c\\to\\infty$ limit), `E7` ($\\overline{\\mathrm{FC}} - 1$, $c\\to 0^+$ limit), plus paired `FCT - E1` on identical bootstrap draws.",
        "",
        "## 表 1：事后探索性（post hoc）预测器去幅度表现（`r_orth_cal [97.5% CI]`）",
        "",
        "| 性质 | 细胞系 (`N`) | 掩码 | `FCT` | `E1` (中位数常数) | `E1/FCT` | `E6` (`ln FC_bar`) | `E7` (`FC_bar - 1`) | 常数 $c_{\\mathrm{med}}$ |",
        "| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |",
    ]

    by_key = {
        (r["target_line"], r["mask"], r["predictor"]): r
        for r in metric_rows
    }
    seen_pairs: List[Tuple[str, str]] = []
    for r in metric_rows:
        key = (r["target_line"], r["mask"])
        if key not in seen_pairs:
            seen_pairs.append(key)

    for ln, msk in seen_pairs:
        rf = by_key[(ln, msk, "FCT")]
        re1 = by_key[(ln, msk, "E1")]
        re6 = by_key[(ln, msk, "E6")]
        re7 = by_key[(ln, msk, "E7")]
        pct = 100.0 * float(re1["ratio_vs_fct"]) if np.isfinite(float(re1["ratio_vs_fct"])) else float("nan")
        lines.append(
            f"| 事后探索性（post hoc） | `{ln}` (`N={rf['n_targets']}`) | `{msk}` "
            f"| `{_fmt(rf['r_orth_cal'])} [{_fmt(rf['r_orth_cal_ci975_lo'])}, {_fmt(rf['r_orth_cal_ci975_hi'])}]` "
            f"| `{_fmt(re1['r_orth_cal'])} [{_fmt(re1['r_orth_cal_ci975_lo'])}, {_fmt(re1['r_orth_cal_ci975_hi'])}]` "
            f"| `{pct:.1f}%` "
            f"| `{_fmt(re6['r_orth_cal'])} [{_fmt(re6['r_orth_cal_ci975_lo'])}, {_fmt(re6['r_orth_cal_ci975_hi'])}]` "
            f"| `{_fmt(re7['r_orth_cal'])} [{_fmt(re7['r_orth_cal_ci975_lo'])}, {_fmt(re7['r_orth_cal_ci975_hi'])}]` "
            f"| `{float(rf['e1_const_median']):.4f}` |"
        )

    lines.extend(
        [
            "",
            "## 表 2：事后探索性（post hoc）配对差 `FCT - E1`（同一组 4,000 次联合 bootstrap 抽样）",
            "",
            "| 性质 | 细胞系 (`N`) | 掩码 | `FCT` | `E1` | `E1 / FCT` | `FCT - E1 [97.5% CI]` | `97.5% CI` 不含 `0` |",
            "| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: |",
        ]
    )
    for pr in pair_rows:
        pct = float(pr["e1_pct_of_fct"])
        lines.append(
            f"| 事后探索性（post hoc） | `{pr['target_line']}` (`N={pr['n_targets']}`) | `{pr['mask']}` "
            f"| `{_fmt(pr['fct_r_orth_cal'])}` | `{_fmt(pr['e1_r_orth_cal'])}` | `{pct:.1f}%` "
            f"| `{_fmt(pr['diff_est'])} [{_fmt(pr['diff_ci975_lo'])}, {_fmt(pr['diff_ci975_hi'])}]` "
            f"| `{pr['diff_ci975_excludes_zero']}` |"
        )
    return "\n".join(lines) + "\n"


def _write_csv(path: Path, rows: Sequence[Dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Independent R10-E Post-Hoc Baseline Ablation Audit (A-T3)")
    parser.add_argument("--lines", nargs="+", default=list(DEFAULT_AUDIT_LINES))
    parser.add_argument("--masks", nargs="+", default=list(EVAL_MASKS))
    parser.add_argument("--n-boot", type=int, default=4000)
    parser.add_argument("--seed", type=int, default=20261122)
    parser.add_argument("--out-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()

    need_r10 = [ln for ln in args.lines if ln not in R9_LINE_SET]
    if need_r10:
        enforce_r10_gate(need_r10, allow_r9_only=False)
    else:
        enforce_r10_gate(args.lines, allow_r9_only=True, extra_lines_dir=R9_LINES_DIR)

    hgnc = CanonicalHGNC.from_tsv()
    common_genes = load_r9_common_genes(hgnc=hgnc)

    r9_loaded = {
        ln: load_line_dataset(R9_LINES_DIR / ln, representation="log1p", aligned_genes=common_genes, hgnc=hgnc)
        for ln in R9_LINE_NAMES
    }

    all_metrics: List[Dict[str, Any]] = []
    all_pairs: List[Dict[str, Any]] = []
    for tgt in args.lines:
        if tgt in R9_LINE_SET:
            tgt_data = r9_loaded[tgt]
            sources = [r9_loaded[s] for s in R9_LINE_NAMES if s != tgt]
        else:
            tgt_data = load_line_dataset(R10_LINES_DIR / tgt, representation="log1p", aligned_genes=common_genes, hgnc=hgnc)
            sources = [r9_loaded[s] for s in R9_LINE_NAMES]
        cohort = align_source_and_target(sources, tgt_data, only_reproducible=True)
        m_rows, p_rows = evaluate_line_ablation(cohort, masks=args.masks, n_boot=args.n_boot, seed=args.seed)
        all_metrics.extend(m_rows)
        all_pairs.extend(p_rows)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(args.out_dir / "audit_r10e_metrics.csv", all_metrics)
    _write_csv(args.out_dir / "audit_r10e_paired_fct_minus_e1.csv", all_pairs)
    (args.out_dir / "audit_r10e_summary.json").write_text(
        json.dumps(
            {
                "analysis_type": "post_hoc_exploratory",
                "n_boot": args.n_boot,
                "seed": args.seed,
                "lines": list(args.lines),
                "masks": list(args.masks),
                "metrics": all_metrics,
                "paired_fct_minus_e1": all_pairs,
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    (args.out_dir / "AUDIT_R10E_REPORT.md").write_text(
        render_posthoc_markdown(all_metrics, all_pairs, n_boot=args.n_boot, seed=args.seed),
        encoding="utf-8",
    )
    print(json.dumps({"status": "R10E_INDEPENDENT_COMPLETE", "n_metric_rows": len(all_metrics), "n_pair_rows": len(all_pairs)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
