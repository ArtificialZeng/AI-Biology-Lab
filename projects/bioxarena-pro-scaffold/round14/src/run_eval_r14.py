#!/usr/bin/env python3
"""
Round 14 (`PREREG_G`, `GSE337804`) Prospective Same-Cell-Line Replication Runner.

Reuses frozen Round 10 / 11 / 12 evaluation machinery cleanly without modifying any
file under `round10/`, `round11/`, `round12/`, or `round13/`:
  - Default seed `DEFAULT_R14_SEED = 20261224` (`n_boot = 4000`, `n_perm = 100`).
  - Hard blinding gate (`PREREG_G_COMMIT.txt` + `CODE_FREEZE_R14.txt`) on `--phase holdout_r14`.
  - Pre-registered Item 7 descriptive stratification of evaluated targets into
    `in_r12_a549_overlap` (also tested in Round 12 `a549_ko`) vs `novel_not_in_r12_a549`,
    reporting `FCT` (`all_nontarget` `r_orth_cal`) and all 6 paired contrasts `FCT - Mis(s')`.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple

import numpy as np
import pandas as pd

ROUND14_SRC_DIR = Path(__file__).resolve().parent
ROUND14_ROOT = ROUND14_SRC_DIR.parent
REPO_ROOT = ROUND14_ROOT.parent
R11_SRC_DIR = REPO_ROOT / "round11" / "src"
R12_SRC_DIR = REPO_ROOT / "round12" / "src"
for p in (str(R11_SRC_DIR), str(R12_SRC_DIR), str(ROUND14_SRC_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from adjudicate_r14 import adjudicate_r14_results, render_r14_markdown_report  # noqa: E402
from data_loader_r14 import (  # noqa: E402
    DEFAULT_HGNC_FILE,
    DEFAULT_R14_SEED,
    HGNCBridge,
    LineDatasetR10,
    R9_COMMON_GENES_FILE,
    R9_LINES_DIR,
    R10_LINES_DIR,
    R11_LINES_DIR,
    R12_A549_LINE_DIR,
    R14_CONFIRMATORY_LINES_FILE,
    R14_LINES_DIR,
    R14_STATUS_DIR,
    check_prereg_g_gate,
    list_standardized_lines,
    load_r9_common_genes,
    load_r12_a549_tested_targets,
    load_source_cohort_r14,
    load_standardized_line_r10,
    parse_confirmatory_lines_r14,
)
from metrics_r11 import compute_residual_suite_r11  # noqa: E402
from metrics_r12 import compute_paired_fct_minus_all_mismatched_r12  # noqa: E402
from run_eval_r12 import (  # noqa: E402
    _build_aligned_predictions_r12,
    _dev_min_overlap,
    evaluate_target_against_sources_r12,
)

RESULTS_DIR_R14 = ROUND14_ROOT / "results"


def compute_overlap_subgroup_descriptive_r14(
    target_ds: LineDatasetR10,
    source_datasets: Dict[str, LineDatasetR10],
    r12_overlap_targets: Set[str],
    only_reproducible: bool = True,
    s5_top75_genes: bool = False,
    include_r9_anchors: bool = True,
    include_exploratory_e: bool = True,
    n_boot: int = 4000,
    seed: int = DEFAULT_R14_SEED,
) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    """
    Compute `PREREG_G` Item 7 reporting-only subgroup metrics (`FCT` and all 6 paired
    `FCT - NegCtrl_FCT_MismatchedMean(s')` contrasts under `all_nontarget`) separately for:
      - `in_r12_a549_overlap`: evaluated targets also present in Round 12 `a549_ko`
      - `novel_not_in_r12_a549`: evaluated targets never tested in Round 12 `a549_ko`
    """
    al = _build_aligned_predictions_r12(
        target_ds=target_ds,
        source_datasets=source_datasets,
        only_reproducible=only_reproducible,
        s5_top75_genes=s5_top75_genes,
        include_r9_anchors=include_r9_anchors,
        include_exploratory_e=include_exploratory_e,
    )
    q_names: List[str] = list(al["q_names"])
    n_q = len(q_names)
    preds: Dict[str, np.ndarray] = al["preds"]
    true_deltas: np.ndarray = al["true_deltas"]
    transfer_ref: np.ndarray = preds["Transfer_Mean"]
    q_gidx = list(al["q_gidx"])
    ha_sub = al["ha_sub"]
    hb_sub = al["hb_sub"]
    hv_sub = al["hv_sub"]
    gene_keep_mask = al["gene_keep_mask"]

    overlap_idx = [i for i, sym in enumerate(q_names) if sym in r12_overlap_targets]
    novel_idx = [i for i, sym in enumerate(q_names) if sym not in r12_overlap_targets]

    sg_dict: Dict[str, Any] = {}
    flat_rows: List[Dict[str, Any]] = []

    for sg_key, idx_list in (
        ("in_r12_a549_overlap", overlap_idx),
        ("novel_not_in_r12_a549", novel_idx),
    ):
        sub_names = [q_names[i] for i in idx_list]
        if not idx_list:
            sg_dict[sg_key] = {"n_targets": 0, "targets": []}
            continue

        arr = np.asarray(idx_list, dtype=np.int64)
        sub_true = true_deltas[arr]
        sub_tm = transfer_ref[arr]
        sub_gidx = [q_gidx[i] for i in idx_list]
        sub_ha = ha_sub[arr] if ha_sub is not None else None
        sub_hb = hb_sub[arr] if hb_sub is not None else None
        sub_hv = hv_sub[arr] if hv_sub is not None else None
        sub_preds = {
            k: v[arr]
            for k, v in preds.items()
            if isinstance(v, np.ndarray) and v.ndim >= 1 and v.shape[0] == n_q
        }

        fct_suite = compute_residual_suite_r11(
            delta_obs=sub_true,
            delta_pred=sub_preds["FCT"],
            delta_transfer=sub_tm,
            target_gene_idx=sub_gidx,
            delta_half_a=sub_ha,
            delta_half_b=sub_hb,
            half_valid=sub_hv,
            gene_keep_mask=gene_keep_mask,
            n_boot=n_boot,
            seed=seed,
        )
        fct_all = fct_suite.get("all_nontarget", {})
        fct_stat = fct_all.get("r_orth_cal", {})
        rho_ab_orth = float(fct_all.get("rho_half_ab_orth", float("nan")))

        paired_map = compute_paired_fct_minus_all_mismatched_r12(
            delta_obs=sub_true,
            preds=sub_preds,
            delta_transfer=sub_tm,
            target_gene_idx=sub_gidx,
            delta_half_a=sub_ha,
            delta_half_b=sub_hb,
            half_valid=sub_hv,
            gene_keep_mask=gene_keep_mask,
            n_boot=n_boot,
            seed=seed,
        )

        sg_dict[sg_key] = {
            "n_targets": len(sub_names),
            "targets": sub_names,
            "rho_half_ab_orth": rho_ab_orth,
            "fct": fct_stat,
            "paired_fct_minus_mismatched": paired_map,
        }

        flat_rows.append(
            {
                "target_line": target_ds.name,
                "subgroup": sg_key,
                "n_targets": len(sub_names),
                "comparison": "FCT_all_nontarget",
                "mismatched_source": "",
                "estimate": fct_stat.get("estimate", float("nan")),
                "ci95_lo": fct_stat.get("ci95_lo", float("nan")),
                "ci95_hi": fct_stat.get("ci95_hi", float("nan")),
                "ci975_lo": fct_stat.get("ci975_lo", float("nan")),
                "ci975_hi": fct_stat.get("ci975_hi", float("nan")),
                "valid_frac": fct_stat.get("valid_frac", float("nan")),
                "rho_half_ab_orth": rho_ab_orth,
                "n_boot": int(n_boot),
            }
        )
        for src_name in sorted(paired_map.keys()):
            pm = paired_map[src_name]
            flat_rows.append(
                {
                    "target_line": target_ds.name,
                    "subgroup": sg_key,
                    "n_targets": len(sub_names),
                    "comparison": f"FCT_minus_Mis({src_name})",
                    "mismatched_source": src_name,
                    "estimate": pm.get("estimate", float("nan")),
                    "ci95_lo": pm.get("ci95_lo", float("nan")),
                    "ci95_hi": pm.get("ci95_hi", float("nan")),
                    "ci975_lo": pm.get("ci975_lo", float("nan")),
                    "ci975_hi": pm.get("ci975_hi", float("nan")),
                    "valid_frac": pm.get("valid_frac", float("nan")),
                    "rho_half_ab_orth": rho_ab_orth,
                    "n_boot": int(n_boot),
                }
            )

    return sg_dict, flat_rows


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Round 14 (`PREREG_G`) prospective runner.")
    p.add_argument("--phase", choices=["dev_r9", "dev_r10", "dev_r11", "holdout_r14"], required=True)
    p.add_argument("--source-mode", choices=["r9_6", "r9_r10_8"], default="r9_6")
    p.add_argument("--repr", choices=["log1p", "zscore", "libnorm"], default="log1p")
    p.add_argument("--all-targets", action="store_true", default=False)
    p.add_argument("--s5-top75", action="store_true", default=False)
    p.add_argument("--no-r9-anchors", action="store_true", default=False)
    p.add_argument("--no-exploratory-e", action="store_true", default=False)
    p.add_argument("--no-canon-hgnc", action="store_true", default=False)
    p.add_argument("--lines", type=str, default="", help="Optional comma-separated target lines.")
    p.add_argument("--r12-a549-dir", type=Path, default=R12_A549_LINE_DIR)
    p.add_argument("--n-boot", type=int, default=4000)
    p.add_argument("--n-perm", type=int, default=100)
    p.add_argument("--seed", type=int, default=DEFAULT_R14_SEED)
    p.add_argument("--tag-suffix", type=str, default="")
    p.add_argument("--out-dir", type=Path, default=RESULTS_DIR_R14)
    return p.parse_args()


def main() -> None:
    t0 = time.time()
    args = parse_args()
    only_rep = not args.all_targets
    include_anchors = not args.no_r9_anchors
    hgnc = None if args.no_canon_hgnc else HGNCBridge.from_tsv(DEFAULT_HGNC_FILE)
    common_r9 = load_r9_common_genes(hgnc_bridge=hgnc, common_genes_file=R9_COMMON_GENES_FILE)
    r12_overlap_targets = load_r12_a549_tested_targets(args.r12_a549_dir, hgnc_bridge=hgnc)

    prereg_hash = freeze_hash = ""
    if args.phase in ("dev_r9", "dev_r10", "dev_r11"):
        src_ds_base = load_source_cohort_r14("r9_6" if args.phase == "dev_r9" else args.source_mode, args.repr, common_r9, hgnc)
        lines_root = {"dev_r9": R9_LINES_DIR, "dev_r10": R10_LINES_DIR, "dev_r11": R11_LINES_DIR}[args.phase]
        avail = list_standardized_lines(lines_root)
        eval_lines = [x.strip() for x in args.lines.split(",") if x.strip()] if args.lines.strip() else avail
        target_ds_map = {
            ln: load_standardized_line_r10(
                lines_root / ln,
                representation=args.repr,
                common_genes=common_r9,
                hgnc_bridge=hgnc,
                min_gene_overlap=_dev_min_overlap(args.phase, ln),
            )
            for ln in eval_lines
        }
        confirmatory_lines: List[str] = []
    else:
        prereg_hash, freeze_hash = check_prereg_g_gate(R14_STATUS_DIR)
        src_ds_base = load_source_cohort_r14(args.source_mode, args.repr, common_r9, hgnc)
        avail = list_standardized_lines(R14_LINES_DIR)
        if not avail:
            raise RuntimeError(f"No processed Round 14 target lines found under {R14_LINES_DIR}.")
        confirmatory_lines = parse_confirmatory_lines_r14(R14_CONFIRMATORY_LINES_FILE)
        conf_set = set(confirmatory_lines)
        eval_lines = [x.strip() for x in args.lines.split(",") if x.strip()] if args.lines.strip() else avail
        target_ds_map = {
            ln: load_standardized_line_r10(
                R14_LINES_DIR / ln,
                representation=args.repr,
                common_genes=common_r9,
                hgnc_bridge=hgnc,
                min_gene_overlap=(3000 if ln in conf_set else 50),
            )
            for ln in eval_lines
        }

    source_pool_names = sorted(src_ds_base.keys())
    all_line_rows: List[Dict[str, Any]] = []
    all_boot_rows: List[Dict[str, Any]] = []
    all_subgroup_rows: List[Dict[str, Any]] = []
    overlap_subgroups_map: Dict[str, Any] = {}
    perm_map: Dict[str, Any] = {}
    n5_map: Dict[str, Any] = {}
    n5_split_map: Dict[str, Any] = {}
    n5_perm_map: Dict[str, Any] = {}
    r19_map: Dict[str, Any] = {}
    s7_map: Dict[str, Any] = {}
    s7_top50_map: Dict[str, Any] = {}
    eval_meta: Dict[str, Any] = {}

    conf_set = set(confirmatory_lines)
    for tgt_ln in eval_lines:
        tgt_ds = target_ds_map[tgt_ln]
        src_ds = {k: v for k, v in src_ds_base.items() if k != tgt_ln}
        min_ov = (
            _dev_min_overlap(args.phase, tgt_ln)
            if args.phase != "holdout_r14"
            else (3000 if tgt_ln in conf_set else 50)
        )
        l_rows, b_rows, extra = evaluate_target_against_sources_r12(
            target_ds=tgt_ds,
            source_datasets=src_ds,
            only_reproducible=only_rep,
            s5_top75_genes=args.s5_top75,
            include_r9_anchors=include_anchors,
            include_exploratory_e=not args.no_exploratory_e,
            min_gene_overlap=min_ov,
            n_boot=args.n_boot,
            seed=args.seed,
            n_perm=args.n_perm,
        )
        sg_dict, sg_rows = compute_overlap_subgroup_descriptive_r14(
            target_ds=tgt_ds,
            source_datasets=src_ds,
            r12_overlap_targets=r12_overlap_targets,
            only_reproducible=only_rep,
            s5_top75_genes=args.s5_top75,
            include_r9_anchors=include_anchors,
            include_exploratory_e=not args.no_exploratory_e,
            n_boot=args.n_boot,
            seed=args.seed,
        )
        overlap_subgroups_map[tgt_ln] = sg_dict
        all_subgroup_rows.extend(sg_rows)

        all_line_rows.extend(l_rows)
        all_boot_rows.extend(b_rows)
        perm_map[tgt_ln] = extra["perm100"]
        n5_map[tgt_ln] = extra["n5_cone_metrics"]
        n5_split_map[tgt_ln] = extra["n5_cone_split_metrics"]
        n5_perm_map[tgt_ln] = extra["n5_perm100"]
        r19_map[tgt_ln] = extra["r19_variance_decomposition"]
        s7_map[tgt_ln] = extra["s7_cellcount_kd"]
        s7_top50_map[tgt_ln] = extra["s7_cellcount_kd_resp_top50"]
        eval_meta[tgt_ln] = {
            "n_genes_eval": extra["n_genes_eval"],
            "n_eval_targets": extra["n_eval_targets"],
            "source_lines": extra["source_lines"],
            "const_stats": extra["const_stats"],
            "n_r12_a549_overlap": int(sg_dict.get("in_r12_a549_overlap", {}).get("n_targets", 0)),
            "n_novel_not_in_r12_a549": int(sg_dict.get("novel_not_in_r12_a549", {}).get("n_targets", 0)),
        }

    df_lines = pd.DataFrame(all_line_rows)
    df_boot = pd.DataFrame(all_boot_rows)
    df_subgroups = pd.DataFrame(all_subgroup_rows)
    summary: Dict[str, Any] = {
        "phase": args.phase,
        "source_mode": args.source_mode,
        "source_lines": source_pool_names,
        "representation": args.repr,
        "primary_mask": "all_nontarget",
        "sa8_mask": "resp_top50_split",
        "only_reproducible": only_rep,
        "s5_top75": args.s5_top75,
        "canon_hgnc": not args.no_canon_hgnc,
        "prereg_commit_hash": prereg_hash,
        "code_freeze_hash": freeze_hash,
        "confirmatory_lines": confirmatory_lines,
        "bootstrap_seed": int(args.seed),
        "n_boot": int(args.n_boot),
        "n_perm": int(args.n_perm),
        "r12_a549_overlap_target_pool_size": len(r12_overlap_targets),
        "eval_meta": eval_meta,
        "overlap_subgroups": overlap_subgroups_map,
        "fct_perm100": perm_map,
        "n5_cone_metrics": n5_map,
        "n5_cone_split_metrics": n5_split_map,
        "n5_perm100": n5_perm_map,
        "r19_variance_decomposition": r19_map,
        "s7_cellcount_kd": s7_map,
        "s7_cellcount_kd_resp_top50": s7_top50_map,
        "elapsed_sec": round(time.time() - t0, 2),
    }

    args.out_dir.mkdir(parents=True, exist_ok=True)
    run_tag = f"{args.phase}_{args.source_mode}_{args.repr}{args.tag_suffix}"
    adj = adjudicate_r14_results(df_lines, df_boot, summary, confirmatory_lines)
    md = render_r14_markdown_report(df_lines, df_boot, summary, adj, run_tag)

    df_lines.to_csv(args.out_dir / f"line_metrics_{run_tag}.csv", index=False)
    df_boot.to_csv(args.out_dir / f"bootstrap_ci_{run_tag}.csv", index=False)
    df_subgroups.to_csv(args.out_dir / f"overlap_subgroups_{run_tag}.csv", index=False)
    (args.out_dir / f"adjudication_{run_tag}.json").write_text(json.dumps(adj, indent=2), encoding="utf-8")
    (args.out_dir / f"summary_{run_tag}.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (args.out_dir / f"REPORT_{run_tag}.md").write_text(md, encoding="utf-8")
    print(f"Saved Round 14 results (*_{run_tag}.*) in {summary['elapsed_sec']}s", flush=True)


if __name__ == "__main__":
    main()
