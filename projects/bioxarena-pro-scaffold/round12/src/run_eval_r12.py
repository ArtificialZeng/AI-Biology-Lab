#!/usr/bin/env python3
"""
Round 12 Prospective Group E Evaluation Runner (`PREREG_E`).

Reuses frozen Round 10 and Round 11 evaluation machinery cleanly without modifying
any file under `round10/` or `round11/`, adding the Round 12 prospective analyses:
  - `H-E1b`: Paired joint bootstrap `FCT - NegCtrl_FCT_MismatchedMean(s)` against EVERY
    mismatched source `s` individually (`97.5%` CIs, `seed=20261201`, `n_boot=4000`).
  - `H-E2b`: Joint cross-mask bootstrap `E1_ConstMed(resp_top50_split) - E1_ConstMed(all_nontarget)`
    over common half-valid perturbations (`rows_common`).
  - `N_5`: Null-cone semi-partial `r_cone_sp` + paired `FCT - Mis_s` + 100-perm
    empirical null (`reporting-only`, explicitly marking by-construction zero controls).
  - Hard blinding gate via `check_prereg_e_gate()` blocking `--phase holdout_r12` prior
    to `PREREG_E_COMMIT.txt` & `CODE_FREEZE_R12.txt`.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

ROUND12_SRC_DIR = Path(__file__).resolve().parent
ROUND12_ROOT = ROUND12_SRC_DIR.parent
REPO_ROOT = ROUND12_ROOT.parent
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
R11_SRC_DIR = REPO_ROOT / "round11" / "src"
for p in (str(R10_SRC_DIR), str(R11_SRC_DIR), str(ROUND12_SRC_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from adjudicate_r12 import adjudicate_r12_results, render_r12_markdown_report  # noqa: E402
from data_loader_r12 import (  # noqa: E402
    DEFAULT_HGNC_FILE,
    HGNCBridge,
    R9_COMMON_GENES_FILE,
    R9_LINES_DIR,
    R10_LINES_DIR,
    R11_LINES_DIR,
    R12_CONFIRMATORY_LINES_FILE,
    R12_LINES_DIR,
    R12_STATUS_DIR,
    LineDatasetR10,
    check_prereg_e_gate,
    list_standardized_lines,
    load_r9_common_genes,
    load_source_cohort_r12,
    load_standardized_line_r10,
    parse_confirmatory_lines_r12,
)
from fct_pipeline import (  # noqa: E402
    FctSourceState,
    build_source_input_r10,
    compute_control_covariance_r10,
    subsample_ctrl_r10,
)
from run_eval_r10 import select_eval_targets_r10  # noqa: E402
from metrics_r12 import (  # noqa: E402
    DEFAULT_R12_SEED,
    compute_joint_sa8_split_minus_all_r12,
    compute_n5_perm100_r12,
    compute_n5_resp_top50_split_r12,
    compute_n5_suite_r12,
    compute_paired_fct_minus_all_mismatched_r12,
)
from predictors_r11 import predict_group_d_suite  # noqa: E402
from run_eval_r11 import _subset_dataset, evaluate_target_against_sources_r11  # noqa: E402

RESULTS_DIR_R12 = ROUND12_ROOT / "results"


def _build_aligned_predictions_r12(
    target_ds: LineDatasetR10,
    source_datasets: Dict[str, LineDatasetR10],
    only_reproducible: bool = True,
    s5_top75_genes: bool = False,
    shrink: float = 0.0,
    cov_max_cells: Optional[int] = None,
    include_r9_anchors: bool = True,
    include_exploratory_e: bool = True,
) -> Dict[str, Any]:
    """Reproduce exact aligned state & predictions needed by R12 paired/SA-8/N_5 modules."""
    src_names = sorted(source_datasets.keys())
    shared_genes_set = set(target_ds.genes).intersection(*(set(source_datasets[s].genes) for s in src_names))
    aligned_genes = [g for g in target_ds.genes if g in shared_genes_set]
    g_map = {g: i for i, g in enumerate(aligned_genes)}
    t_pos = {g: i for i, g in enumerate(target_ds.genes)}
    t_cols = np.asarray([t_pos[g] for g in aligned_genes], dtype=np.int64)
    t_sub_ds = _subset_dataset(target_ds, aligned_genes, t_cols, g_map)

    src_inputs = []
    for sn in src_names:
        s_raw = source_datasets[sn]
        s_pos = {g: i for i, g in enumerate(s_raw.genes)}
        s_cols = np.asarray([s_pos[g] for g in aligned_genes], dtype=np.int64)
        src_inputs.append(
            build_source_input_r10(
                _subset_dataset(s_raw, aligned_genes, s_cols, g_map),
                only_reproducible=only_reproducible,
                shrink=shrink,
                cov_max_cells=cov_max_cells,
                compute_cov=include_r9_anchors,
                seed=42,
            )
        )

    state = FctSourceState.fit(src_inputs, mask_self=True)
    row_idx, q_names, q_gidx = select_eval_targets_r10(t_sub_ds, src_inputs, only_reproducible=only_reproducible)
    t_ctrl_cov = None
    if include_r9_anchors:
        sub_c = subsample_ctrl_r10(t_sub_ds.ctrl_cells, cov_max_cells=cov_max_cells, seed=42)
        t_ctrl_cov = compute_control_covariance_r10(sub_c, shrink=shrink)

    preds, fc_bundle, _ = predict_group_d_suite(
        state=state,
        target_ctrl_mean_libnorm=t_sub_ds.ctrl_mean_libnorm,
        query_target_names=q_names,
        query_target_gene_indices=q_gidx,
        target_ctrl_cov=t_ctrl_cov,
        target_ctrl_mean_libnorm_full=t_sub_ds.ctrl_mean_libnorm_full,
        include_r9_anchors=include_r9_anchors,
        include_exploratory_e=include_exploratory_e,
    )
    gene_keep_mask = None
    if s5_top75_genes:
        t_var = np.var(np.asarray(t_sub_ds.ctrl_cells, dtype=np.float64), axis=0, ddof=1)
        gene_keep_mask = t_var > float(np.percentile(t_var, 25.0))

    return {
        "state": state,
        "preds": preds,
        "fc_bundle": fc_bundle,
        "t_sub_ds": t_sub_ds,
        "true_deltas": t_sub_ds.deltas[row_idx],
        "q_names": q_names,
        "q_gidx": q_gidx,
        "ha_sub": t_sub_ds.half_a_deltas[row_idx] if t_sub_ds.half_a_deltas is not None else None,
        "hb_sub": t_sub_ds.half_b_deltas[row_idx] if t_sub_ds.half_b_deltas is not None else None,
        "hv_sub": t_sub_ds.half_valid[row_idx] if t_sub_ds.half_valid is not None else None,
        "gene_keep_mask": gene_keep_mask,
    }


def evaluate_target_against_sources_r12(
    target_ds: LineDatasetR10,
    source_datasets: Dict[str, LineDatasetR10],
    only_reproducible: bool = True,
    s5_top75_genes: bool = False,
    shrink: float = 0.0,
    cov_max_cells: Optional[int] = None,
    include_r9_anchors: bool = True,
    include_exploratory_e: bool = True,
    min_gene_overlap: int = 3000,
    n_boot: int = 4000,
    seed: int = DEFAULT_R12_SEED,
    n_perm: int = 100,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], Dict[str, Any]]:
    """
    Evaluate one target line `target_ds` zero-shot against `source_datasets` under `PREREG_E`:
      1. Runs frozen `evaluate_target_against_sources_r11` for all shared Group D metrics,
         bootstrap CIs, dual/SA-8 `perm100`, `S7`, and `R19` using `seed`.
      2. Appends `H-E1b` paired `FCT - Mis_s` rows for EVERY source line `s`.
      3. Appends `H-E2b` joint `E1(resp_top50_split) - E1(all_nontarget)` (plus diagnostic `FCT`).
      4. Computes reporting-only `N_5` null-cone suite (`r_cone_sp`, paired `FCT - Mis_s`,
         and single-draw `Perm0..n_perm-1`).
    """
    line_rows, boot_rows, extra = evaluate_target_against_sources_r11(
        target_ds=target_ds,
        source_datasets=source_datasets,
        only_reproducible=only_reproducible,
        s5_top75_genes=s5_top75_genes,
        shrink=shrink,
        cov_max_cells=cov_max_cells,
        include_r9_anchors=include_r9_anchors,
        include_exploratory_e=include_exploratory_e,
        min_gene_overlap=min_gene_overlap,
        n_boot=n_boot,
        seed=seed,
        n_perm=n_perm,
    )

    al = _build_aligned_predictions_r12(
        target_ds=target_ds,
        source_datasets=source_datasets,
        only_reproducible=only_reproducible,
        s5_top75_genes=s5_top75_genes,
        shrink=shrink,
        cov_max_cells=cov_max_cells,
        include_r9_anchors=include_r9_anchors,
        include_exploratory_e=include_exploratory_e,
    )
    preds = al["preds"]
    state = al["state"]
    _, fc_bar, _ = al["fc_bundle"]
    t_sub_ds = al["t_sub_ds"]
    true_deltas = al["true_deltas"]
    transfer_ref = preds["Transfer_Mean"]
    q_gidx = al["q_gidx"]
    ha_sub = al["ha_sub"]
    hb_sub = al["hb_sub"]
    hv_sub = al["hv_sub"]
    gene_keep_mask = al["gene_keep_mask"]

    # 1. H-E1b: Paired FCT - Mis_s across EVERY mismatched source
    paired_mis_map = compute_paired_fct_minus_all_mismatched_r12(
        delta_obs=true_deltas,
        preds=preds,
        delta_transfer=transfer_ref,
        target_gene_idx=q_gidx,
        delta_half_a=ha_sub,
        delta_half_b=hb_sub,
        half_valid=hv_sub,
        gene_keep_mask=gene_keep_mask,
        n_boot=n_boot,
        seed=seed,
    )
    fct_row = next(r for r in line_rows if r["method"] == "FCT")
    row_by_method = {r["method"]: r for r in line_rows}
    for src_name, pci in paired_mis_map.items():
        mk = f"NegCtrl_FCT_MismatchedMean({src_name})"
        boot_rows.append(
            {
                "target_line": target_ds.name,
                "comparison": f"H-E1b_FCT_minus_Mis_{src_name}_all_nontarget",
                "metric": "r_orth_cal_all_nontarget",
                "method_a": "FCT",
                "method_b": mk,
                "mean_a": float(fct_row["r_orth_cal"]),
                "mean_b": float(row_by_method[mk]["r_orth_cal"]),
                "estimate": pci["estimate"],
                "ci95_lo": pci["ci95_lo"],
                "ci95_hi": pci["ci95_hi"],
                "ci975_lo": pci["ci975_lo"],
                "ci975_hi": pci["ci975_hi"],
                "ci9833_lo": pci["ci9833_lo"],
                "ci9833_hi": pci["ci9833_hi"],
                "valid_frac": pci["valid_frac"],
                "n_items": pci["n_items"],
                "n_boot": pci["n_boot"],
            }
        )

    # 2. H-E2b: Joint cross-mask E1_ConstMed(resp_top50_split) - E1_ConstMed(all_nontarget)
    #    plus diagnostic FCT(resp_top50_split) - FCT(all_nontarget)
    for m_key, comp_name in (
        ("E1_ConstMed", "H-E2b_E1_split_minus_all"),
        ("FCT", "diag_FCT_split_minus_all"),
    ):
        jdiff = compute_joint_sa8_split_minus_all_r12(
            delta_obs=true_deltas,
            delta_pred=preds[m_key],
            delta_transfer=transfer_ref,
            target_gene_idx=q_gidx,
            delta_half_a=ha_sub,
            delta_half_b=hb_sub,
            half_valid=hv_sub,
            gene_keep_mask=gene_keep_mask,
            top_k=50,
            n_boot=n_boot,
            seed=seed,
        )
        boot_rows.append(
            {
                "target_line": target_ds.name,
                "comparison": comp_name,
                "metric": "r_orth_cal_sa8_minus_all",
                "method_a": f"{m_key}(resp_top50_split)",
                "method_b": f"{m_key}(all_nontarget_common)",
                "mean_a": float(jdiff["split_summary"]["estimate"]),
                "mean_b": float(jdiff["all_common_summary"]["estimate"]),
                "estimate": jdiff["estimate"],
                "ci95_lo": jdiff["ci95_lo"],
                "ci95_hi": jdiff["ci95_hi"],
                "ci975_lo": jdiff["ci975_lo"],
                "ci975_hi": jdiff["ci975_hi"],
                "ci9833_lo": jdiff["ci9833_lo"],
                "ci9833_hi": jdiff["ci9833_hi"],
                "valid_frac": jdiff["valid_frac"],
                "n_items": jdiff["n_common_perts"],
                "n_boot": jdiff["n_boot"],
            }
        )

    # 3. Reporting-only N_5 Null-Cone metrics (all_nontarget & resp_top50_split) + 100-perm null
    n5_suite = compute_n5_suite_r12(
        delta_obs=true_deltas,
        preds=preds,
        delta_transfer=transfer_ref,
        grand_mean=state.grand_mean,
        fc_bar=fc_bar,
        target_gene_idx=q_gidx,
        delta_half_a=ha_sub,
        delta_half_b=hb_sub,
        half_valid=hv_sub,
        gene_keep_mask=gene_keep_mask,
        n_boot=n_boot,
        seed=seed,
    )
    n5_split_suite = compute_n5_resp_top50_split_r12(
        delta_obs=true_deltas,
        preds=preds,
        delta_transfer=transfer_ref,
        grand_mean=state.grand_mean,
        fc_bar=fc_bar,
        target_gene_idx=q_gidx,
        delta_half_a=ha_sub,
        delta_half_b=hb_sub,
        half_valid=hv_sub,
        gene_keep_mask=gene_keep_mask,
        top_k=50,
        n_boot=n_boot,
        seed=seed,
    )
    n5_perm = compute_n5_perm100_r12(
        state=state,
        target_ctrl_mean_libnorm=t_sub_ds.ctrl_mean_libnorm,
        query_target_names=al["q_names"],
        delta_obs=true_deltas,
        delta_transfer=transfer_ref,
        grand_mean=state.grand_mean,
        target_gene_idx=q_gidx,
        delta_half_a=ha_sub,
        delta_half_b=hb_sub,
        half_valid=hv_sub,
        fc_bundle=al["fc_bundle"],
        gene_keep_mask=gene_keep_mask,
        n_perm=n_perm,
    )
    extra["n5_cone_metrics"] = n5_suite
    extra["n5_cone_split_metrics"] = n5_split_suite
    extra["n5_perm100"] = n5_perm
    return line_rows, boot_rows, extra


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Round 12 Group E (`PREREG_E`) prospective runner.")
    p.add_argument("--phase", choices=["dev_r9", "dev_r10", "dev_r11", "holdout_r12"], required=True)
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
    p.add_argument("--seed", type=int, default=DEFAULT_R12_SEED)
    p.add_argument("--tag-suffix", type=str, default="")
    p.add_argument("--out-dir", type=Path, default=RESULTS_DIR_R12)
    return p.parse_args()


def _dev_min_overlap(phase: str, ln: str) -> int:
    if phase == "dev_r9":
        return 3000
    if phase == "dev_r10":
        return 3000 if ln in ("dld1", "cd4t_rest") else 50
    if phase == "dev_r11":
        return 50 if ln == "feng_ipsc" else 3000
    return 3000


def main() -> None:
    t0 = time.time()
    args = parse_args()
    only_rep = not args.all_targets
    include_anchors = not args.no_r9_anchors
    hgnc = None if args.no_canon_hgnc else HGNCBridge.from_tsv(DEFAULT_HGNC_FILE)
    common_r9 = load_r9_common_genes(hgnc_bridge=hgnc, common_genes_file=R9_COMMON_GENES_FILE)

    prereg_hash = freeze_hash = ""
    if args.phase == "dev_r9":
        src_ds_base = load_source_cohort_r12("r9_6", representation=args.repr, common_genes=common_r9, hgnc_bridge=hgnc)
        avail = list_standardized_lines(R9_LINES_DIR)
        eval_lines = [x.strip() for x in args.lines.split(",") if x.strip()] if args.lines.strip() else avail
        target_ds_map = {
            ln: load_standardized_line_r10(R9_LINES_DIR / ln, representation=args.repr, common_genes=common_r9, hgnc_bridge=hgnc)
            for ln in eval_lines
        }
        confirmatory_lines: List[str] = []
    elif args.phase == "dev_r10":
        src_ds_base = load_source_cohort_r12(args.source_mode, representation=args.repr, common_genes=common_r9, hgnc_bridge=hgnc)
        avail = list_standardized_lines(R10_LINES_DIR)
        eval_lines = [x.strip() for x in args.lines.split(",") if x.strip()] if args.lines.strip() else avail
        target_ds_map = {
            ln: load_standardized_line_r10(
                R10_LINES_DIR / ln,
                representation=args.repr,
                common_genes=common_r9,
                hgnc_bridge=hgnc,
                min_gene_overlap=_dev_min_overlap("dev_r10", ln),
            )
            for ln in eval_lines
        }
        confirmatory_lines = []
    elif args.phase == "dev_r11":
        src_ds_base = load_source_cohort_r12(args.source_mode, representation=args.repr, common_genes=common_r9, hgnc_bridge=hgnc)
        avail = list_standardized_lines(R11_LINES_DIR)
        eval_lines = [x.strip() for x in args.lines.split(",") if x.strip()] if args.lines.strip() else avail
        target_ds_map = {
            ln: load_standardized_line_r10(
                R11_LINES_DIR / ln,
                representation=args.repr,
                common_genes=common_r9,
                hgnc_bridge=hgnc,
                min_gene_overlap=_dev_min_overlap("dev_r11", ln),
            )
            for ln in eval_lines
        }
        confirmatory_lines = []
    else:
        prereg_hash, freeze_hash = check_prereg_e_gate(R12_STATUS_DIR)
        src_ds_base = load_source_cohort_r12(args.source_mode, representation=args.repr, common_genes=common_r9, hgnc_bridge=hgnc)
        avail = list_standardized_lines(R12_LINES_DIR)
        if not avail:
            raise RuntimeError(f"No processed Round 12 target lines found under {R12_LINES_DIR}.")
        confirmatory_lines = parse_confirmatory_lines_r12(R12_CONFIRMATORY_LINES_FILE)
        conf_set = set(confirmatory_lines)
        eval_lines = [x.strip() for x in args.lines.split(",") if x.strip()] if args.lines.strip() else avail
        target_ds_map = {
            ln: load_standardized_line_r10(
                R12_LINES_DIR / ln,
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
            if args.phase != "holdout_r12"
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
    adj = adjudicate_r12_results(df_lines, df_boot, summary, confirmatory_lines)
    md = render_r12_markdown_report(df_lines, df_boot, summary, adj, run_tag)

    df_lines.to_csv(args.out_dir / f"line_metrics_{run_tag}.csv", index=False)
    df_boot.to_csv(args.out_dir / f"bootstrap_ci_{run_tag}.csv", index=False)
    (args.out_dir / f"adjudication_{run_tag}.json").write_text(json.dumps(adj, indent=2), encoding="utf-8")
    (args.out_dir / f"summary_{run_tag}.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (args.out_dir / f"REPORT_{run_tag}.md").write_text(md, encoding="utf-8")
    print(f"Saved Round 12 results (*_{run_tag}.*) in {summary['elapsed_sec']}s", flush=True)


if __name__ == "__main__":
    main()
