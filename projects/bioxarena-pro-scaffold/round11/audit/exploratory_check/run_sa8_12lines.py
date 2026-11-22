#!/usr/bin/env python3
"""
Exploratory Check 1 (`SA-8` on 12 Cell Lines):
Re-runs `E1_ConstMed` and `FCT` under `all_nontarget`, `resp_top50`, and `resp_top50_split`,
plus `Perm100` (`resp_top50_split` & `all_nontarget`) and joint cross-mask contrasts
(`E1(resp_top50) - E1(all)` and `E1(resp_top50_split) - E1(all)`) across all 12 human
Perturb-seq lines using frozen `round11/audit/` `r_orth_cal` with:
    seed = 20261124, n_boot = 4000, n_perm = 100.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

EXPLORATORY_CHECK_DIR = Path(__file__).resolve().parent
AUDIT_DIR = EXPLORATORY_CHECK_DIR.parent
for p in (str(EXPLORATORY_CHECK_DIR), str(AUDIT_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from exploratory_core import (  # noqa: E402
    DEFAULT_AUDIT_SEED,
    LINES_12,
    build_r9_lolo_cohort,
    evaluate_cross_mask_split_minus_all,
    load_12_lines,
)
from independent_r11_metrics import (  # noqa: E402
    evaluate_cross_mask_cal_diff,
    evaluate_fct_perm_null,
    evaluate_joint_residuals,
    evaluate_sa8_perm_null,
    evaluate_sa8_resp_top50_split,
)
from independent_r11_predictors import (  # noqa: E402
    build_predictions,
    compute_source_fc_bar,
)


def run_check1(
    seed: int = DEFAULT_AUDIT_SEED,
    n_boot: int = 4000,
    n_perm: int = 100,
    out_dir: Path = EXPLORATORY_CHECK_DIR / "outputs",
) -> Dict[str, Any]:
    t0 = time.time()
    out_dir.mkdir(parents=True, exist_ok=True)
    lines = load_12_lines()

    rows: List[Dict[str, Any]] = []
    per_line: Dict[str, Any] = {}

    for ln, grp in LINES_12:
        cohort = build_r9_lolo_cohort(lines, ln)
        fc_bar, t_mean, _ = compute_source_fc_bar(
            eval_targets=cohort.eval_targets,
            source_names=cohort.source_names,
            source_deltas_by_line=cohort.source_deltas_by_line,
            source_ctrl_mean_libnorm=cohort.source_ctrl_mean_libnorm,
        )
        preds, _, _, const_stats = build_predictions(
            cohort,
            include_sa3_anchors=False,
            include_exploratory_e=False,
        )

        # Evaluate E1 and FCT across all_nontarget, resp_top50, and resp_top50_split
        e1_all = evaluate_joint_residuals(
            delta_obs=cohort.delta_obs,
            pred=preds["e1_const_med"],
            t_mean=t_mean,
            target_cols=cohort.target_cols,
            mode="all_nontarget",
            n_boot=n_boot,
            seed=seed,
            half_a_obs=cohort.half_a_obs,
            half_b_obs=cohort.half_b_obs,
            half_valid=cohort.half_valid,
        )
        e1_top50 = evaluate_joint_residuals(
            delta_obs=cohort.delta_obs,
            pred=preds["e1_const_med"],
            t_mean=t_mean,
            target_cols=cohort.target_cols,
            mode="resp_top50",
            n_boot=n_boot,
            seed=seed,
            half_a_obs=cohort.half_a_obs,
            half_b_obs=cohort.half_b_obs,
            half_valid=cohort.half_valid,
        )
        e1_split = evaluate_sa8_resp_top50_split(
            delta_obs=cohort.delta_obs,
            pred=preds["e1_const_med"],
            t_mean=t_mean,
            target_cols=cohort.target_cols,
            half_a_obs=cohort.half_a_obs,
            half_b_obs=cohort.half_b_obs,
            half_valid=cohort.half_valid,
            top_k=50,
            n_boot=n_boot,
            seed=seed,
        )
        e1_top50_minus_all = evaluate_cross_mask_cal_diff(
            delta_obs=cohort.delta_obs,
            pred=preds["e1_const_med"],
            t_mean=t_mean,
            target_cols=cohort.target_cols,
            mode_a="resp_top50",
            mode_b="all_nontarget",
            n_boot=n_boot,
            seed=seed,
            half_a_obs=cohort.half_a_obs,
            half_b_obs=cohort.half_b_obs,
            half_valid=cohort.half_valid,
        )
        e1_split_minus_all = evaluate_cross_mask_split_minus_all(
            delta_obs=cohort.delta_obs,
            pred=preds["e1_const_med"],
            t_mean=t_mean,
            target_cols=cohort.target_cols,
            half_a_obs=cohort.half_a_obs,
            half_b_obs=cohort.half_b_obs,
            half_valid=cohort.half_valid,
            top_k=50,
            n_boot=n_boot,
            seed=seed,
        )

        fct_all = evaluate_joint_residuals(
            delta_obs=cohort.delta_obs,
            pred=preds["fct"],
            t_mean=t_mean,
            target_cols=cohort.target_cols,
            mode="all_nontarget",
            n_boot=n_boot,
            seed=seed,
            half_a_obs=cohort.half_a_obs,
            half_b_obs=cohort.half_b_obs,
            half_valid=cohort.half_valid,
        )
        fct_top50 = evaluate_joint_residuals(
            delta_obs=cohort.delta_obs,
            pred=preds["fct"],
            t_mean=t_mean,
            target_cols=cohort.target_cols,
            mode="resp_top50",
            n_boot=n_boot,
            seed=seed,
            half_a_obs=cohort.half_a_obs,
            half_b_obs=cohort.half_b_obs,
            half_valid=cohort.half_valid,
        )
        fct_split = evaluate_sa8_resp_top50_split(
            delta_obs=cohort.delta_obs,
            pred=preds["fct"],
            t_mean=t_mean,
            target_cols=cohort.target_cols,
            half_a_obs=cohort.half_a_obs,
            half_b_obs=cohort.half_b_obs,
            half_valid=cohort.half_valid,
            top_k=50,
            n_boot=n_boot,
            seed=seed,
        )

        # Single-half ceilings for Perm100 SA-8 (matching frozen SA-8 ceiling definition)
        perm_all = evaluate_fct_perm_null(
            delta_obs=cohort.delta_obs,
            t_mean=t_mean,
            fc_bar=fc_bar,
            tgt_ctrl_mean_libnorm=cohort.tgt_ctrl_mean_libnorm,
            target_cols=cohort.target_cols,
            ceil_raw=float(fct_all["noise_ceiling_r"]),
            ceil_orth=float(fct_all["noise_ceiling_r_orth"]),
            mode="all_nontarget",
            n_perm=n_perm,
        )
        perm_top50 = evaluate_fct_perm_null(
            delta_obs=cohort.delta_obs,
            t_mean=t_mean,
            fc_bar=fc_bar,
            tgt_ctrl_mean_libnorm=cohort.tgt_ctrl_mean_libnorm,
            target_cols=cohort.target_cols,
            ceil_raw=float(fct_top50["noise_ceiling_r"]),
            ceil_orth=float(fct_top50["noise_ceiling_r_orth"]),
            mode="resp_top50",
            n_perm=n_perm,
        )
        perm_split = evaluate_sa8_perm_null(
            delta_obs=cohort.delta_obs,
            t_mean=t_mean,
            fc_bar=fc_bar,
            tgt_ctrl_mean_libnorm=cohort.tgt_ctrl_mean_libnorm,
            target_cols=cohort.target_cols,
            half_a_obs=cohort.half_a_obs,
            half_b_obs=cohort.half_b_obs,
            half_valid=cohort.half_valid,
            ceil_single_ab=float(fct_split["ceil_single_half_AtoB"]),
            ceil_single_ba=float(fct_split["ceil_single_half_BtoA"]),
            top_k=50,
            n_perm=n_perm,
        )

        per_line[ln] = {
            "group": grp,
            "n_genes_eval": len(cohort.eval_genes),
            "n_eval_targets": len(cohort.eval_targets),
            "n_sources": len(cohort.source_names),
            "const_stats": const_stats,
            "e1_all": e1_all,
            "e1_resp_top50": e1_top50,
            "e1_resp_top50_split": e1_split,
            "e1_top50_minus_all": e1_top50_minus_all,
            "e1_split_minus_all": e1_split_minus_all,
            "fct_all": fct_all,
            "fct_resp_top50": fct_top50,
            "fct_resp_top50_split": fct_split,
            "perm100_all": perm_all,
            "perm100_resp_top50": perm_top50,
            "perm100_resp_top50_split": perm_split,
        }

        row = {
            "line": ln,
            "group": grp,
            "n_genes_eval": len(cohort.eval_genes),
            "n_eval_targets": len(cohort.eval_targets),
            "e1_all_est": e1_all["r_orth_cal_est"],
            "e1_all_ci975_lo": e1_all["r_orth_cal_ci975_lo"],
            "e1_all_ci975_hi": e1_all["r_orth_cal_ci975_hi"],
            "e1_top50_est": e1_top50["r_orth_cal_est"],
            "e1_top50_ci975_lo": e1_top50["r_orth_cal_ci975_lo"],
            "e1_top50_ci975_hi": e1_top50["r_orth_cal_ci975_hi"],
            "e1_split_est": e1_split["r_orth_cal_split_est"],
            "e1_split_ci975_lo": e1_split["r_orth_cal_split_ci975_lo"],
            "e1_split_ci975_hi": e1_split["r_orth_cal_split_ci975_hi"],
            "e1_split_valid_frac": e1_split["r_orth_cal_split_valid_frac"],
            "e1_split_minus_all_est": e1_split_minus_all["diff_est"],
            "e1_split_minus_all_ci975_lo": e1_split_minus_all["diff_ci975_lo"],
            "e1_split_minus_all_ci975_hi": e1_split_minus_all["diff_ci975_hi"],
            "e1_split_minus_all_valid_frac": e1_split_minus_all["diff_valid_frac"],
            "e1_top50_minus_all_est": e1_top50_minus_all["diff_est"],
            "e1_top50_minus_all_ci975_lo": e1_top50_minus_all["diff_ci975_lo"],
            "e1_top50_minus_all_ci975_hi": e1_top50_minus_all["diff_ci975_hi"],
            "fct_all_est": fct_all["r_orth_cal_est"],
            "fct_top50_est": fct_top50["r_orth_cal_est"],
            "fct_split_est": fct_split["r_orth_cal_split_est"],
            "fct_split_ci975_lo": fct_split["r_orth_cal_split_ci975_lo"],
            "fct_split_ci975_hi": fct_split["r_orth_cal_split_ci975_hi"],
            "perm100_split_med": perm_split["r_orth_cal_split"]["median"],
            "perm100_split_p975": perm_split["r_orth_cal_split"]["p975"],
        }
        rows.append(row)
        print(
            f"[Check1 SA-8] {ln:24s} N={row['n_eval_targets']:4d} | "
            f"E1(all)={row['e1_all_est']:+.4f} E1(top50)={row['e1_top50_est']:+.4f} "
            f"E1(split)={row['e1_split_est']:+.4f} [{row['e1_split_ci975_lo']:+.4f},{row['e1_split_ci975_hi']:+.4f}] | "
            f"E1(split-all)={row['e1_split_minus_all_est']:+.4f} [{row['e1_split_minus_all_ci975_lo']:+.4f},{row['e1_split_minus_all_ci975_hi']:+.4f}] | "
            f"FCT(split)={row['fct_split_est']:+.4f} Perm(split)={row['perm100_split_med']:+.4f}/{row['perm100_split_p975']:+.4f}"
        )

    payload = {
        "check": "SA-8_12_lines",
        "seed": seed,
        "n_boot": n_boot,
        "n_perm": n_perm,
        "elapsed_sec": round(time.time() - t0, 2),
        "per_line": per_line,
    }
    with open(out_dir / "sa8_12lines.json", "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    with open(out_dir / "sa8_12lines.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return payload


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=DEFAULT_AUDIT_SEED)
    parser.add_argument("--n-boot", type=int, default=4000)
    parser.add_argument("--n-perm", type=int, default=100)
    parser.add_argument("--out-dir", type=Path, default=EXPLORATORY_CHECK_DIR / "outputs")
    args = parser.parse_args()
    run_check1(seed=args.seed, n_boot=args.n_boot, n_perm=args.n_perm, out_dir=args.out_dir)
