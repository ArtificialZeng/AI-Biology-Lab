#!/usr/bin/env python3
"""
Round 13 Prospective Group E Replication Runner (`PREREG_F`, `GSE327644`).

Reuses frozen Round 10 / 11 / 12 evaluation machinery cleanly without modifying any
file under `round10/`, `round11/`, or `round12/`, using `DEFAULT_R13_SEED = 20261223`
and hard `PREREG_F_COMMIT.txt` + `CODE_FREEZE_R13.txt` blinding gates for `--phase holdout_r13`.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd

ROUND13_SRC_DIR = Path(__file__).resolve().parent
ROUND13_ROOT = ROUND13_SRC_DIR.parent
REPO_ROOT = ROUND13_ROOT.parent
R12_SRC_DIR = REPO_ROOT / "round12" / "src"
for p in (str(R12_SRC_DIR), str(ROUND13_SRC_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from adjudicate_r13 import adjudicate_r13_results, render_r13_markdown_report  # noqa: E402
from data_loader_r13 import (  # noqa: E402
    DEFAULT_HGNC_FILE,
    DEFAULT_R13_SEED,
    HGNCBridge,
    R9_COMMON_GENES_FILE,
    R9_LINES_DIR,
    R10_LINES_DIR,
    R11_LINES_DIR,
    R13_CONFIRMATORY_LINES_FILE,
    R13_LINES_DIR,
    R13_STATUS_DIR,
    check_prereg_f_gate,
    list_standardized_lines,
    load_r9_common_genes,
    load_source_cohort_r13,
    load_standardized_line_r10,
    parse_confirmatory_lines_r13,
)
from run_eval_r12 import _dev_min_overlap, evaluate_target_against_sources_r12  # noqa: E402

RESULTS_DIR_R13 = ROUND13_ROOT / "results"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Round 13 (`PREREG_F`) prospective runner.")
    p.add_argument("--phase", choices=["dev_r9", "dev_r10", "dev_r11", "holdout_r13"], required=True)
    p.add_argument("--source-mode", choices=["r9_6", "r9_r10_8"], default="r9_6")
    p.add_argument("--repr", choices=["log1p", "zscore", "libnorm"], default="log1p")
    p.add_argument("--all-targets", action="store_true", default=False)
    p.add_argument("--s5-top75", action="store_true", default=False)
    p.add_argument("--no-r9-anchors", action="store_true", default=False)
    p.add_argument("--no-exploratory-e", action="store_true", default=False)
    p.add_argument("--no-canon-hgnc", action="store_true", default=False)
    p.add_argument("--lines", type=str, default="", help="Optional comma-separated target lines.")
    p.add_argument("--n-boot", type=int, default=4000)
    p.add_argument("--n-perm", type=int, default=100)
    p.add_argument("--seed", type=int, default=DEFAULT_R13_SEED)
    p.add_argument("--tag-suffix", type=str, default="")
    p.add_argument("--out-dir", type=Path, default=RESULTS_DIR_R13)
    return p.parse_args()


def main() -> None:
    t0 = time.time()
    args = parse_args()
    only_rep = not args.all_targets
    include_anchors = not args.no_r9_anchors
    hgnc = None if args.no_canon_hgnc else HGNCBridge.from_tsv(DEFAULT_HGNC_FILE)
    common_r9 = load_r9_common_genes(hgnc_bridge=hgnc, common_genes_file=R9_COMMON_GENES_FILE)

    prereg_hash = freeze_hash = ""
    if args.phase in ("dev_r9", "dev_r10", "dev_r11"):
        src_ds_base = load_source_cohort_r13("r9_6" if args.phase == "dev_r9" else args.source_mode, args.repr, common_r9, hgnc)
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
        prereg_hash, freeze_hash = check_prereg_f_gate(R13_STATUS_DIR)
        src_ds_base = load_source_cohort_r13(args.source_mode, args.repr, common_r9, hgnc)
        avail = list_standardized_lines(R13_LINES_DIR)
        if not avail:
            raise RuntimeError(f"No processed Round 13 target lines found under {R13_LINES_DIR}.")
        confirmatory_lines = parse_confirmatory_lines_r13(R13_CONFIRMATORY_LINES_FILE)
        conf_set = set(confirmatory_lines)
        eval_lines = [x.strip() for x in args.lines.split(",") if x.strip()] if args.lines.strip() else avail
        target_ds_map = {
            ln: load_standardized_line_r10(
                R13_LINES_DIR / ln,
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
            if args.phase != "holdout_r13"
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
        }

    df_lines = pd.DataFrame(all_line_rows)
    df_boot = pd.DataFrame(all_boot_rows)
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
        "eval_meta": eval_meta,
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
    adj = adjudicate_r13_results(df_lines, df_boot, summary, confirmatory_lines)
    md = render_r13_markdown_report(df_lines, df_boot, summary, adj, run_tag)

    df_lines.to_csv(args.out_dir / f"line_metrics_{run_tag}.csv", index=False)
    df_boot.to_csv(args.out_dir / f"bootstrap_ci_{run_tag}.csv", index=False)
    (args.out_dir / f"adjudication_{run_tag}.json").write_text(json.dumps(adj, indent=2), encoding="utf-8")
    (args.out_dir / f"summary_{run_tag}.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (args.out_dir / f"REPORT_{run_tag}.md").write_text(md, encoding="utf-8")
    print(f"Saved Round 13 results (*_{run_tag}.*) in {summary['elapsed_sec']}s", flush=True)


if __name__ == "__main__":
    main()
