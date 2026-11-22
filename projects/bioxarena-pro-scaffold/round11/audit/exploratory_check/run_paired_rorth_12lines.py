#!/usr/bin/env python3
"""
Paired `FCT - Mis_s` specificity check across all 12 cell lines using frozen
`round11/audit/` `r_orth_cal` implementation (`evaluate_paired_cal_diff`,
`evaluate_joint_residuals`, `evaluate_fct_perm_null`) under `all_nontarget`:
    locked P5 cohort (`only_reproducible=True`), R9 leave-one-out (`5` sources
    for R9 lines; `6` sources for R10/R11 lines), `seed = 20261124`,
    `n_boot = 4,000`, `n_perm = 100`, two-sided Bonferroni `97.5%` CIs.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

EXPLORATORY_CHECK_DIR = Path(__file__).resolve().parent
AUDIT_DIR = EXPLORATORY_CHECK_DIR.parent
for p in (str(EXPLORATORY_CHECK_DIR), str(AUDIT_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from exploratory_core import (  # noqa: E402
    DEFAULT_AUDIT_SEED,
    LINES_12,
    build_r9_lolo_cohort,
    load_12_lines,
)
from independent_r11_loader import R9_LINE_NAMES  # noqa: E402
from independent_r11_metrics import (  # noqa: E402
    evaluate_fct_perm_null,
    evaluate_joint_residuals,
    evaluate_paired_cal_diff,
)
from independent_r11_predictors import (  # noqa: E402
    build_predictions,
    compute_source_fc_bar,
)


def run_paired_rorth_check(
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
        preds, _, _, _ = build_predictions(
            cohort,
            include_sa3_anchors=False,
            include_exploratory_e=False,
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
        perm100_all = evaluate_fct_perm_null(
            delta_obs=cohort.delta_obs,
            t_mean=t_mean,
            fc_bar=fc_bar,
            tgt_ctrl_mean_libnorm=cohort.tgt_ctrl_mean_libnorm,
            target_cols=cohort.target_cols,
            ceil_raw=float(fct_all.get("noise_ceiling_r", 0.0)),
            ceil_orth=float(fct_all.get("noise_ceiling_r_orth", 0.0)),
            mode="all_nontarget",
            n_perm=n_perm,
        )

        line_pair_results: Dict[str, Any] = {}
        line_mis_results: Dict[str, Any] = {}
        paired_pass_flags: List[bool] = []
        unpaired_mis_his: List[float] = []

        for src in R9_LINE_NAMES:
            if src == ln:
                continue
            mis_res = evaluate_joint_residuals(
                delta_obs=cohort.delta_obs,
                pred=preds[f"fct_mismatch::{src}"],
                t_mean=t_mean,
                target_cols=cohort.target_cols,
                mode="all_nontarget",
                n_boot=n_boot,
                seed=seed,
                half_a_obs=cohort.half_a_obs,
                half_b_obs=cohort.half_b_obs,
                half_valid=cohort.half_valid,
            )
            pair_res = evaluate_paired_cal_diff(
                delta_obs=cohort.delta_obs,
                pred_a=preds["fct"],
                pred_b=preds[f"fct_mismatch::{src}"],
                t_mean=t_mean,
                target_cols=cohort.target_cols,
                mode="all_nontarget",
                n_boot=n_boot,
                seed=seed,
                half_a_obs=cohort.half_a_obs,
                half_b_obs=cohort.half_b_obs,
                half_valid=cohort.half_valid,
                ortho=True,
            )
            line_mis_results[src] = mis_res
            line_pair_results[src] = pair_res

            diff_lo = float(pair_res.get("diff_ci975_lo", float("nan")))
            diff_vf = float(pair_res.get("diff_valid_frac", 0.0))
            mis_hi = float(mis_res.get("r_orth_cal_ci975_hi", float("nan")))
            paired_pass_flags.append(bool(np.isfinite(diff_lo) and diff_lo > 0.0 and diff_vf > 0.0))
            unpaired_mis_his.append(mis_hi)

        n_mis = len(cohort.source_names)
        n_paired_pass = int(sum(paired_pass_flags))
        all_paired_pass = bool(n_paired_pass == n_mis and n_mis > 0)

        fct_est = float(fct_all.get("r_orth_cal_est", float("nan")))
        fct_lo = float(fct_all.get("r_orth_cal_ci975_lo", float("nan")))
        fct_hi = float(fct_all.get("r_orth_cal_ci975_hi", float("nan")))
        fct_vf = float(fct_all.get("r_orth_cal_valid_frac", 0.0))
        fct_rho = float(fct_all.get("rho_ab_orth", float("nan")))
        n_eval = int(fct_all.get("n_eval_perts", 0))
        perm_p975 = float(perm100_all["r_orth_cal"]["p975"])
        max_mis_hi = float(np.nanmax(unpaired_mis_his)) if unpaired_mis_his else float("nan")
        max_mis_src = max(
            cohort.source_names,
            key=lambda s: float(line_mis_results[s].get("r_orth_cal_ci975_hi", float("-inf"))),
        )
        unpaired_hd2_pass = bool(
            np.isfinite(fct_lo)
            and np.isfinite(max_mis_hi)
            and np.isfinite(perm_p975)
            and fct_lo > max(max_mis_hi, perm_p975)
        )
        hd1_pass = bool(
            np.isfinite(fct_lo)
            and fct_lo > 0.0
            and np.isfinite(fct_vf)
            and fct_vf >= 0.95
            and np.isfinite(fct_rho)
            and fct_rho > 0.0
        )
        he1c_pass = bool(np.isfinite(fct_lo) and np.isfinite(perm_p975) and fct_lo > perm_p975)
        he1_full_pass = bool(hd1_pass and all_paired_pass and he1c_pass)

        # Track bottleneck paired contrast (minimum CI975_lo)
        min_pair_src = min(
            cohort.source_names,
            key=lambda s: float(line_pair_results[s].get("diff_ci975_lo", float("inf"))),
        )
        min_pair_res = line_pair_results[min_pair_src]

        per_line[ln] = {
            "group": grp,
            "n_eval": n_eval,
            "n_mis_sources": n_mis,
            "fct_r_orth_cal": {
                "est": fct_est,
                "ci975_lo": fct_lo,
                "ci975_hi": fct_hi,
                "valid_frac": fct_vf,
                "rho_ab_orth": fct_rho,
            },
            "perm100_p975": perm_p975,
            "perm100_median": float(perm100_all["r_orth_cal"]["median"]),
            "max_mis_source": max_mis_src,
            "max_mis_ci975_hi": max_mis_hi,
            "min_paired_diff_source": min_pair_src,
            "min_paired_diff_est": float(min_pair_res.get("diff_est", float("nan"))),
            "min_paired_diff_ci975_lo": float(min_pair_res.get("diff_ci975_lo", float("nan"))),
            "min_paired_diff_ci975_hi": float(min_pair_res.get("diff_ci975_hi", float("nan"))),
            "n_paired_pass": n_paired_pass,
            "all_mis_paired_ci975_gt_0": all_paired_pass,
            "unpaired_hd2_pass": unpaired_hd2_pass,
            "hd1_pass": hd1_pass,
            "he1c_pass": he1c_pass,
            "he1_full_pass": he1_full_pass,
            "by_source": {},
        }

        for src in R9_LINE_NAMES:
            if src == ln:
                per_line[ln]["by_source"][src] = {"status": "SELF_EXCLUDED_LOO"}
                rows.append({
                    "line": ln,
                    "group": grp,
                    "n_eval": n_eval,
                    "n_mis_sources": n_mis,
                    "source_line": src,
                    "status": "SELF_EXCLUDED_LOO",
                    "fct_est": fct_est,
                    "fct_ci975_lo": fct_lo,
                    "fct_ci975_hi": fct_hi,
                    "fct_valid_frac": fct_vf,
                    "fct_rho_ab_orth": fct_rho,
                    "perm100_p975": perm_p975,
                    "mis_est": float("nan"),
                    "mis_ci975_lo": float("nan"),
                    "mis_ci975_hi": float("nan"),
                    "mis_valid_frac": float("nan"),
                    "unpaired_fct_lo_gt_mis_hi": "",
                    "paired_diff_est": float("nan"),
                    "paired_diff_ci975_lo": float("nan"),
                    "paired_diff_ci975_hi": float("nan"),
                    "paired_diff_valid_frac": float("nan"),
                    "paired_ci975_gt_0": "",
                    "line_n_paired_pass": n_paired_pass,
                    "line_all_mis_paired_pass": all_paired_pass,
                    "line_unpaired_hd2_pass": unpaired_hd2_pass,
                    "line_he1_full_pass": he1_full_pass,
                })
                continue

            mr = line_mis_results[src]
            pr = line_pair_results[src]
            mis_est = float(mr.get("r_orth_cal_est", float("nan")))
            mis_lo = float(mr.get("r_orth_cal_ci975_lo", float("nan")))
            mis_hi = float(mr.get("r_orth_cal_ci975_hi", float("nan")))
            mis_vf = float(mr.get("r_orth_cal_valid_frac", 0.0))
            diff_est = float(pr.get("diff_est", float("nan")))
            diff_lo = float(pr.get("diff_ci975_lo", float("nan")))
            diff_hi = float(pr.get("diff_ci975_hi", float("nan")))
            diff_vf = float(pr.get("diff_valid_frac", 0.0))
            unpaired_ok = bool(np.isfinite(fct_lo) and np.isfinite(mis_hi) and fct_lo > mis_hi)
            paired_ok = bool(np.isfinite(diff_lo) and diff_lo > 0.0 and diff_vf > 0.0)

            per_line[ln]["by_source"][src] = {
                "status": "EVALUATED",
                "mis_r_orth_cal": {
                    "est": mis_est,
                    "ci975_lo": mis_lo,
                    "ci975_hi": mis_hi,
                    "valid_frac": mis_vf,
                },
                "paired_fct_minus_mis": pr,
                "unpaired_fct_lo_gt_mis_hi": unpaired_ok,
                "paired_ci975_gt_0": paired_ok,
            }
            rows.append({
                "line": ln,
                "group": grp,
                "n_eval": n_eval,
                "n_mis_sources": n_mis,
                "source_line": src,
                "status": "EVALUATED",
                "fct_est": fct_est,
                "fct_ci975_lo": fct_lo,
                "fct_ci975_hi": fct_hi,
                "fct_valid_frac": fct_vf,
                "fct_rho_ab_orth": fct_rho,
                "perm100_p975": perm_p975,
                "mis_est": mis_est,
                "mis_ci975_lo": mis_lo,
                "mis_ci975_hi": mis_hi,
                "mis_valid_frac": mis_vf,
                "unpaired_fct_lo_gt_mis_hi": unpaired_ok,
                "paired_diff_est": diff_est,
                "paired_diff_ci975_lo": diff_lo,
                "paired_diff_ci975_hi": diff_hi,
                "paired_diff_valid_frac": diff_vf,
                "paired_ci975_gt_0": paired_ok,
                "line_n_paired_pass": n_paired_pass,
                "line_all_mis_paired_pass": all_paired_pass,
                "line_unpaired_hd2_pass": unpaired_hd2_pass,
                "line_he1_full_pass": he1_full_pass,
            })

        print(
            f"[Paired r_orth_cal] {ln:24s} ({grp}) N={n_eval:4d} | "
            f"FCT={fct_est:+.4f} [{fct_lo:+.4f},{fct_hi:+.4f}] (vf={fct_vf:.3f}) | "
            f"unpaired_HD2={'PASS' if unpaired_hd2_pass else 'FAIL'} | "
            f"paired={n_paired_pass}/{n_mis} ({'ALL_PASS' if all_paired_pass else 'FAIL'}) | "
            f"min_pair={min_pair_src}: {float(min_pair_res['diff_est']):+.4f} "
            f"[{float(min_pair_res['diff_ci975_lo']):+.4f}, {float(min_pair_res['diff_ci975_hi']):+.4f}]"
        )

    all12_pass = sum(1 for ln, _ in LINES_12 if per_line[ln]["all_mis_paired_ci975_gt_0"])
    r9_pass = sum(1 for ln, grp in LINES_12 if grp == "R9" and per_line[ln]["all_mis_paired_ci975_gt_0"])
    r10_r11_pass = sum(1 for ln, grp in LINES_12 if grp in ("R10", "R11") and per_line[ln]["all_mis_paired_ci975_gt_0"])
    unpaired_hd2_cnt = sum(1 for ln, _ in LINES_12 if per_line[ln]["unpaired_hd2_pass"])
    he1_full_cnt = sum(1 for ln, _ in LINES_12 if per_line[ln]["he1_full_pass"])

    print(
        f"\nSUMMARY: all_mis_paired_ci975_gt_0 = {all12_pass}/12 "
        f"(R9 LOO 5/5 = {r9_pass}/6; R10+R11 unseen 6/6 = {r10_r11_pass}/6) | "
        f"unpaired_HD2 = {unpaired_hd2_cnt}/12 | full_HE1 = {he1_full_cnt}/12"
    )

    csv_path = out_dir / "paired_rorth_12lines.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    payload = {
        "check": "paired_fct_minus_mis_r_orth_cal_12lines",
        "seed": seed,
        "n_boot": n_boot,
        "n_perm": n_perm,
        "elapsed_sec": round(time.time() - t0, 2),
        "summary_counts": {
            "all_mis_paired_ci975_gt_0_of_12": all12_pass,
            "r9_loo_all5_paired_pass_of_6": r9_pass,
            "r10_r11_all6_paired_pass_of_6": r10_r11_pass,
            "unpaired_hd2_pass_of_12": unpaired_hd2_cnt,
            "full_he1_pass_of_12": he1_full_cnt,
        },
        "per_line": per_line,
    }
    json_path = out_dir / "paired_rorth_12lines.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    return payload


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=DEFAULT_AUDIT_SEED)
    parser.add_argument("--n-boot", type=int, default=4000)
    parser.add_argument("--n-perm", type=int, default=100)
    parser.add_argument("--out-dir", type=Path, default=EXPLORATORY_CHECK_DIR / "outputs")
    args = parser.parse_args()
    run_paired_rorth_check(seed=args.seed, n_boot=args.n_boot, n_perm=args.n_perm, out_dir=args.out_dir)
