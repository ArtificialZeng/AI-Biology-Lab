#!/usr/bin/env python3
"""
Round 10 Evaluation Runner (`PREREG_C` Group C `FCT`).

Execution modes:
  - `--phase dev_r9`: Leave-one-line-out software/development verification strictly on
    Round 9's existing 6 lines (`round9/p2/shared/lines/`). Never accesses any Round 10 data.
  - `--phase holdout_r10`: Zero-shot confirmatory evaluation of Round 10 target lines
    using all 6 Round 9 lines as source cohort (`round9/p2/shared/lines/`), strictly
    blocked unless `PREREG_C_COMMIT.txt` AND `CODE_FREEZE_R10.txt` exist (`check_prereg_c_gate`).
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from adjudicate_r10 import SIX_CONVENTIONS, adjudicate_r10_results, render_r10_markdown_report
from data_loader_r10 import (
    CONFIRMATORY_LINES_FILE,
    DEFAULT_HGNC_FILE,
    HGNCBridge,
    LineDatasetR10,
    R10_LINES_DIR,
    R10_STATUS_DIR,
    R9_COMMON_GENES_FILE,
    R9_LINES_DIR,
    check_prereg_c_gate,
    list_standardized_lines,
    load_r9_common_genes,
    load_standardized_line_r10,
)
from fct_pipeline import (
    FctSourceState,
    build_source_input_r10,
    compute_control_covariance_r10,
    predict_group_c_suite,
    subsample_ctrl_r10,
)
from metrics_r10 import (
    compute_fct_perm100_null,
    compute_joint_paired_r_orth_diff_r10,
    compute_r19_variance_decomposition_r10,
    compute_residual_suite_r10,
    compute_row_metrics_r10,
    compute_s7_cellcount_and_kd_analysis,
    paired_row_bootstrap_ci_r10,
)

ROUND10_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROUND10_ROOT / "results"


def select_eval_targets_r10(
    target_ds: LineDatasetR10,
    source_inputs: Sequence[Any],
    only_reproducible: bool = True,
) -> Tuple[np.ndarray, List[str], np.ndarray]:
    """
    Exact `P5` target eligibility:
      1) target symbol in `G_eval(t)` (`target_indices >= 0`)
      2) reproducible on target line (`repr_mask`) when `only_reproducible=True`
      3) available (and reproducible, since sources were filtered) in >=1 source line.
    """
    valid = target_ds.target_indices >= 0
    if only_reproducible and target_ds.repr_mask is not None:
        valid = valid & target_ds.repr_mask
    src_targets = set()
    for s in source_inputs:
        src_targets.update(s.target_names)
    in_src = np.array([t in src_targets for t in target_ds.target_names], dtype=bool)
    valid = valid & in_src
    row_idx = np.where(valid)[0]
    q_names = [target_ds.target_names[int(i)] for i in row_idx]
    q_gidx = target_ds.target_indices[row_idx]
    return row_idx, q_names, q_gidx


def evaluate_target_against_sources_r10(
    target_ds: LineDatasetR10,
    source_datasets: Dict[str, LineDatasetR10],
    only_reproducible: bool = True,
    s5_top75_genes: bool = False,
    shrink: float = 0.0,
    cov_max_cells: Optional[int] = None,
    include_r9_anchors: bool = True,
    n_boot: int = 4000,
    seed: int = 20261122,
    n_perm: int = 100,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], Dict[str, Any]]:
    """Evaluate one target line `target_ds` strictly zero-shot from `source_datasets`."""
    # Intersect gene space strictly via existing target finite aligned genes & all source genes
    src_names = sorted(source_datasets.keys())
    src_gene_sets = [set(source_datasets[s].genes) for s in src_names]
    shared_genes_set = set(target_ds.genes).intersection(*src_gene_sets)
    if len(shared_genes_set) < 3000:
        raise RuntimeError(
            f"Target line {target_ds.name} gene overlap with sources {len(shared_genes_set)} < 3000."
        )
    aligned_genes = [g for g in target_ds.genes if g in shared_genes_set]
    t_pos = {g: i for i, g in enumerate(target_ds.genes)}
    t_cols = np.asarray([t_pos[g] for g in aligned_genes], dtype=np.int64)
    g_map = {g: i for i, g in enumerate(aligned_genes)}

    def _subset_ds(ds: LineDatasetR10, cols: np.ndarray) -> LineDatasetR10:
        return LineDatasetR10(
            name=ds.name,
            genes=aligned_genes,
            target_names=ds.target_names,
            target_indices=np.asarray([g_map.get(t, -1) for t in ds.target_names], dtype=np.int64),
            deltas=ds.deltas[:, cols],
            ctrl_cells=ds.ctrl_cells[:, cols],
            ctrl_mean_log1p=ds.ctrl_mean_log1p[cols] if ds.ctrl_mean_log1p is not None else None,
            ctrl_mean_libnorm=ds.ctrl_mean_libnorm[cols] if ds.ctrl_mean_libnorm is not None else None,
            ctrl_mean_libnorm_full=(
                ds.ctrl_mean_libnorm_full[cols] if ds.ctrl_mean_libnorm_full is not None else None
            ),
            repr_mask=ds.repr_mask,
            half_a_deltas=ds.half_a_deltas[:, cols] if ds.half_a_deltas is not None else None,
            half_b_deltas=ds.half_b_deltas[:, cols] if ds.half_b_deltas is not None else None,
            half_valid=ds.half_valid,
            n_cells=ds.n_cells,
            knockdown_depth_log1p=ds.knockdown_depth_log1p,
            knockdown_fraction=ds.knockdown_fraction,
            target_ensembl_ids=ds.target_ensembl_ids,
            representation=ds.representation,
            metadata=ds.metadata,
        )

    t_sub_ds = _subset_ds(target_ds, t_cols)
    src_inputs = []
    for sn in src_names:
        s_raw = source_datasets[sn]
        s_pos = {g: i for i, g in enumerate(s_raw.genes)}
        s_cols = np.asarray([s_pos[g] for g in aligned_genes], dtype=np.int64)
        src_inputs.append(
            build_source_input_r10(
                _subset_ds(s_raw, s_cols),
                only_reproducible=only_reproducible,
                shrink=shrink,
                cov_max_cells=cov_max_cells,
                compute_cov=include_r9_anchors,
                seed=42,
            )
        )

    state = FctSourceState.fit(src_inputs, mask_self=True)
    row_idx, q_names, q_gidx = select_eval_targets_r10(
        t_sub_ds, src_inputs, only_reproducible=only_reproducible
    )
    if len(row_idx) == 0:
        raise RuntimeError(f"Zero evaluated targets found for {target_ds.name}.")

    if t_sub_ds.ctrl_mean_libnorm is None:
        raise ValueError(f"Target line {target_ds.name} missing `ctrl_mean_libnorm`.")

    t_ctrl_cov = None
    if include_r9_anchors:
        sub_c = subsample_ctrl_r10(t_sub_ds.ctrl_cells, cov_max_cells=cov_max_cells, seed=42)
        t_ctrl_cov = compute_control_covariance_r10(sub_c, shrink=shrink)

    preds, fc_bundle = predict_group_c_suite(
        state=state,
        target_ctrl_mean_libnorm=t_sub_ds.ctrl_mean_libnorm,
        query_target_names=q_names,
        query_target_gene_indices=q_gidx,
        target_ctrl_cov=t_ctrl_cov,
        target_ctrl_mean_libnorm_full=t_sub_ds.ctrl_mean_libnorm_full,
        include_r9_anchors=include_r9_anchors,
    )

    true_deltas = t_sub_ds.deltas[row_idx]
    transfer_ref = preds["Transfer_Mean"]
    ha_sub = t_sub_ds.half_a_deltas[row_idx] if t_sub_ds.half_a_deltas is not None else None
    hb_sub = t_sub_ds.half_b_deltas[row_idx] if t_sub_ds.half_b_deltas is not None else None
    hv_sub = t_sub_ds.half_valid[row_idx] if t_sub_ds.half_valid is not None else None
    n_cells_sub = t_sub_ds.n_cells[row_idx] if t_sub_ds.n_cells is not None else None
    kd_sub = t_sub_ds.knockdown_depth_log1p[row_idx] if t_sub_ds.knockdown_depth_log1p is not None else None

    # Optional S5 top-75% genes ranked by unperturbed target control variance (`var > Q_0.25`)
    gene_keep_mask = None
    if s5_top75_genes:
        t_var = np.var(np.asarray(t_sub_ds.ctrl_cells, dtype=np.float64), axis=0, ddof=1)
        var_thresh = float(np.percentile(t_var, 25.0))
        gene_keep_mask = t_var > var_thresh

    line_rows: List[Dict[str, Any]] = []
    boot_rows: List[Dict[str, Any]] = []
    row_metrics_by_method: Dict[str, Dict[str, np.ndarray]] = {}

    for m_name, p_mat in preds.items():
        rm_ctrl = compute_row_metrics_r10(
            p_mat, true_deltas, q_gidx, centroid=None, shared_axis=None, mask_target_in_disc=True
        )
        rm_cent = compute_row_metrics_r10(
            p_mat,
            true_deltas,
            q_gidx,
            centroid=state.grand_mean,
            shared_axis=state.axis_shared,
            mask_target_in_disc=True,
        )
        row_metrics_by_method[m_name] = {
            "ctrl_disc_raw": rm_ctrl["disc_raw"],
            "ctrl_disc_raw_panel": rm_ctrl["disc_raw_panel"],
            "ctrl_de50_r": rm_ctrl["de50_r"],
            "cent_disc_raw": rm_cent["disc_raw"],
            "cent_disc_raw_panel": rm_cent["disc_raw_panel"],
            "cent_de50_r": rm_cent["de50_r"],
            "pds_cosine_panel": rm_ctrl.get("pds_cosine_panel", rm_ctrl["disc_raw_panel"]),
            "mse": rm_ctrl["mse"],
        }
        res_suite = compute_residual_suite_r10(
            delta_obs=true_deltas,
            delta_pred=p_mat,
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
        t50 = res_suite["resp_top50"]
        tall = res_suite["all_nontarget"]
        ro = t50.get("r_orth_cal", {})
        ro2 = t50.get("r_orth_cal_sq", {})
        rc = t50.get("r_cal", {})
        sl = t50.get("slope_ols", {})
        sdr = t50.get("sd_pred_over_sd_obs", {})
        dem = t50.get("slope_deming", {})
        rev = t50.get("slope_obs_on_pred", {})
        aro = tall.get("r_orth_cal", {})

        line_rows.append(
            {
                "target_line": target_ds.name,
                "method": m_name,
                "n_genes_eval": len(aligned_genes),
                "n_eval_targets": int(len(q_names)),
                "r_orth_cal": ro.get("estimate", float("nan")),
                "r_orth_cal_ci95_lo": ro.get("ci95_lo", float("nan")),
                "r_orth_cal_ci95_hi": ro.get("ci95_hi", float("nan")),
                "r_orth_cal_ci975_lo": ro.get("ci975_lo", float("nan")),
                "r_orth_cal_ci975_hi": ro.get("ci975_hi", float("nan")),
                "r_orth_cal_valid_frac": ro.get("valid_frac", float("nan")),
                "r_orth_cal_sq": ro2.get("estimate", float("nan")),
                "r_orth_cal_sq_ci975_lo": ro2.get("ci975_lo", float("nan")),
                "r_orth_cal_sq_ci975_hi": ro2.get("ci975_hi", float("nan")),
                "r_orth_raw": t50.get("r_orth_raw", {}).get("estimate", float("nan")),
                "r_cal": rc.get("estimate", float("nan")),
                "r_cal_ci975_lo": rc.get("ci975_lo", float("nan")),
                "r_cal_ci975_hi": rc.get("ci975_hi", float("nan")),
                "r_raw": t50.get("r_raw", {}).get("estimate", float("nan")),
                "slope_ols": sl.get("estimate", float("nan")),
                "slope_ols_ci975_lo": sl.get("ci975_lo", float("nan")),
                "slope_ols_ci975_hi": sl.get("ci975_hi", float("nan")),
                "sd_pred_over_sd_obs": sdr.get("estimate", float("nan")),
                "slope_deming": dem.get("estimate", float("nan")),
                "slope_obs_on_pred": rev.get("estimate", float("nan")),
                "cos_agreement_top50": res_suite["cos_agreement_top50"],
                "rho_half_ab_orth": t50.get("rho_half_ab_orth", float("nan")),
                "ceiling_r_orth": t50.get("ceiling_r_orth", float("nan")),
                "rho_half_ab": t50.get("rho_half_ab", float("nan")),
                "ceiling_r": t50.get("ceiling_r", float("nan")),
                "all_nontarget_r_orth_cal": aro.get("estimate", float("nan")),
                "all_nontarget_r_orth_cal_ci975_lo": aro.get("ci975_lo", float("nan")),
                "all_nontarget_r_orth_cal_ci975_hi": aro.get("ci975_hi", float("nan")),
                "all_nontarget_rho_half_ab_orth": tall.get("rho_half_ab_orth", float("nan")),
                "all_nontarget_ceiling_r_orth": tall.get("ceiling_r_orth", float("nan")),
                **{c: float(np.mean(row_metrics_by_method[m_name][c])) for c in SIX_CONVENTIONS},
                "pds_cosine_panel": float(np.mean(row_metrics_by_method[m_name]["pds_cosine_panel"])),
                "mse": float(np.mean(row_metrics_by_method[m_name]["mse"])),
            }
        )

    # Paired bootstrap comparisons (`H-C3` FCT vs Transfer_Mean across 6 conventions + SA-3 anchors)
    pair_list = [("FCT", "Transfer_Mean")]
    if "FCT_FullNTMean" in preds:
        pair_list.append(("FCT_FullNTMean", "FCT"))
    for anc in ("Hybrid_ProjShared", "Transfer_Mean_VarScaled", "Hybrid_on_TMVS"):
        if anc in preds:
            pair_list.append((anc, "Transfer_Mean"))
    for ma, mb in pair_list:
        for conv in SIX_CONVENTIONS + ("pds_cosine_panel",):
            ci = paired_row_bootstrap_ci_r10(
                row_metrics_by_method[ma][conv],
                row_metrics_by_method[mb][conv],
                n_boot=n_boot,
                seed=seed,
            )
            boot_rows.append(
                {
                    "target_line": target_ds.name,
                    "metric": conv,
                    "method_a": ma,
                    "method_b": mb,
                    **ci,
                }
            )

    fct_row_summary = next(r for r in line_rows if r["method"] == "FCT")
    perm100 = compute_fct_perm100_null(
        state=state,
        target_ctrl_mean_libnorm=t_sub_ds.ctrl_mean_libnorm,
        query_target_names=q_names,
        delta_obs=true_deltas,
        delta_transfer=transfer_ref,
        target_gene_idx=q_gidx,
        ceil_r_top50=float(fct_row_summary["ceiling_r"]),
        ceil_r_orth_top50=float(fct_row_summary["ceiling_r_orth"]),
        fc_bundle=fc_bundle,
        n_perm=n_perm,
        top_k=50,
        gene_keep_mask=gene_keep_mask,
    )

    r19_decomp = compute_r19_variance_decomposition_r10(
        delta_obs=true_deltas,
        delta_transfer=transfer_ref,
        target_gene_idx=q_gidx,
        min_perts_per_gene=5,
        top_k=50,
        n_boot=n_boot,
        seed=seed,
    )

    mismatch_preds = {
        k: v for k, v in preds.items() if k.startswith("NegCtrl_FCT_MismatchedMean(")
    }
    s7_analysis = compute_s7_cellcount_and_kd_analysis(
        state=state,
        query_target_names=q_names,
        target_ctrl_mean_libnorm=t_sub_ds.ctrl_mean_libnorm,
        fc_bundle=fc_bundle,
        delta_obs=true_deltas,
        delta_pred_fct=preds["FCT"],
        delta_transfer=transfer_ref,
        target_gene_idx=q_gidx,
        n_cells=n_cells_sub,
        knockdown_depth=kd_sub,
        delta_half_a=ha_sub,
        delta_half_b=hb_sub,
        half_valid=hv_sub,
        mismatch_preds=mismatch_preds,
        top_k=50,
        n_boot=n_boot,
        seed=seed,
    )

    extra_line_info: Dict[str, Any] = {
        "n_genes_eval": len(aligned_genes),
        "n_eval_targets": int(len(q_names)),
        "source_lines": src_names,
        "perm100": perm100,
        "r19_variance_decomposition": r19_decomp,
        "s7_cellcount_kd": s7_analysis,
    }
    if ha_sub is not None and hb_sub is not None and "Hybrid_on_TMVS" in preds and "Transfer_Mean_VarScaled" in preds:
        extra_line_info["r14_hybrid_on_tmvs_vs_tmvs"] = compute_joint_paired_r_orth_diff_r10(
            delta_obs=true_deltas,
            delta_pred_a=preds["Hybrid_on_TMVS"],
            delta_pred_b=preds["Transfer_Mean_VarScaled"],
            delta_transfer=transfer_ref,
            target_gene_idx=q_gidx,
            delta_half_a=ha_sub,
            delta_half_b=hb_sub,
            half_valid=hv_sub,
            n_boot=n_boot,
            seed=seed,
        )

    return line_rows, boot_rows, extra_line_info


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Round 10 Group C (`FCT`) evaluation runner.")
    p.add_argument("--phase", choices=["dev_r9", "holdout_r10"], required=True)
    p.add_argument("--repr", choices=["log1p", "zscore", "libnorm"], default="log1p")
    p.add_argument("--only-reproducible", action="store_true", default=True)
    p.add_argument("--all-targets", action="store_true", default=False, help="S1 sensitivity: include non-reproducible targets.")
    p.add_argument("--s5-top75", action="store_true", default=False, help="S5 sensitivity: restrict residual mask to top-75%% expr/var genes.")
    p.add_argument("--no-r9-anchors", action="store_true", default=False, help="Skip SA-3 R9 anchors covariance computations.")
    p.add_argument("--no-canon-hgnc", action="store_true", default=False, help="Disable HGNC canonicalization (for raw R9 parity check).")
    p.add_argument("--lines", type=str, default="", help="Optional comma-separated target lines to evaluate.")
    p.add_argument("--n-boot", type=int, default=4000)
    p.add_argument("--n-perm", type=int, default=100)
    p.add_argument("--seed", type=int, default=20261122)
    p.add_argument("--tag-suffix", type=str, default="")
    return p.parse_args()


def main() -> None:
    t0 = time.time()
    args = parse_args()
    only_rep = False if args.all_targets else True
    include_anchors = not args.no_r9_anchors
    hgnc = None if args.no_canon_hgnc else HGNCBridge.from_tsv(DEFAULT_HGNC_FILE)
    common_r9 = load_r9_common_genes(hgnc_bridge=hgnc, common_genes_file=R9_COMMON_GENES_FILE)

    prereg_hash = ""
    freeze_hash = ""
    confirmatory_lines: List[str] = []

    if args.phase == "dev_r9":
        avail_r9 = list_standardized_lines(R9_LINES_DIR)
        if not avail_r9:
            raise RuntimeError(f"Missing Round 9 standardized lines under {R9_LINES_DIR}.")
        all_ds = {
            ln: load_standardized_line_r10(
                R9_LINES_DIR / ln, representation=args.repr, common_genes=common_r9, hgnc_bridge=hgnc
            )
            for ln in avail_r9
        }
        if args.lines.strip():
            eval_lines = [x.strip() for x in args.lines.split(",") if x.strip()]
        else:
            eval_lines = avail_r9
        confirmatory_lines = []  # dev_r9 is strictly non-confirmatory software development check
    else:
        prereg_hash, freeze_hash = check_prereg_c_gate(R10_STATUS_DIR)
        avail_r9 = list_standardized_lines(R9_LINES_DIR)
        avail_r10 = list_standardized_lines(R10_LINES_DIR)
        if not avail_r10:
            raise RuntimeError(f"No processed Round 10 target lines found under {R10_LINES_DIR}.")
        r9_ds = {
            ln: load_standardized_line_r10(
                R9_LINES_DIR / ln, representation=args.repr, common_genes=common_r9, hgnc_bridge=hgnc
            )
            for ln in avail_r9
        }
        r10_ds = {
            ln: load_standardized_line_r10(
                R10_LINES_DIR / ln, representation=args.repr, common_genes=common_r9, hgnc_bridge=hgnc
            )
            for ln in avail_r10
        }
        eval_lines = [x.strip() for x in args.lines.split(",") if x.strip()] if args.lines.strip() else avail_r10
        if CONFIRMATORY_LINES_FILE.exists() and CONFIRMATORY_LINES_FILE.read_text(encoding="utf-8").strip():
            confirmatory_lines = [
                x.strip()
                for x in CONFIRMATORY_LINES_FILE.read_text(encoding="utf-8").splitlines()
                if x.strip() and not x.strip().startswith("#")
            ]

    all_line_rows: List[Dict[str, Any]] = []
    all_boot_rows: List[Dict[str, Any]] = []
    perm_map: Dict[str, Any] = {}
    r19_map: Dict[str, Any] = {}
    s7_map: Dict[str, Any] = {}
    r14_map: Dict[str, Any] = {}
    eval_meta: Dict[str, Any] = {}

    for tgt_ln in eval_lines:
        if args.phase == "dev_r9":
            tgt_ds = all_ds[tgt_ln]
            src_ds = {k: v for k, v in all_ds.items() if k != tgt_ln}
        else:
            tgt_ds = r10_ds[tgt_ln]
            src_ds = r9_ds

        l_rows, b_rows, extra = evaluate_target_against_sources_r10(
            target_ds=tgt_ds,
            source_datasets=src_ds,
            only_reproducible=only_rep,
            s5_top75_genes=args.s5_top75,
            include_r9_anchors=include_anchors,
            n_boot=args.n_boot,
            seed=args.seed,
            n_perm=args.n_perm,
        )
        all_line_rows.extend(l_rows)
        all_boot_rows.extend(b_rows)
        perm_map[tgt_ln] = extra["perm100"]
        r19_map[tgt_ln] = extra["r19_variance_decomposition"]
        s7_map[tgt_ln] = extra["s7_cellcount_kd"]
        if "r14_hybrid_on_tmvs_vs_tmvs" in extra:
            r14_map[tgt_ln] = extra["r14_hybrid_on_tmvs_vs_tmvs"]
        eval_meta[tgt_ln] = {
            "n_genes_eval": extra["n_genes_eval"],
            "n_eval_targets": extra["n_eval_targets"],
            "source_lines": extra["source_lines"],
        }

    df_lines = pd.DataFrame(all_line_rows)
    df_boot = pd.DataFrame(all_boot_rows)
    summary: Dict[str, Any] = {
        "phase": args.phase,
        "representation": args.repr,
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
        "r19_variance_decomposition": r19_map,
        "s7_cellcount_kd": s7_map,
        "r14_hybrid_on_tmvs_vs_tmvs": r14_map,
        "elapsed_sec": round(time.time() - t0, 2),
    }

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    tag = f"{args.phase}_{args.repr}{args.tag_suffix}"
    df_lines.to_csv(RESULTS_DIR / f"line_metrics_{tag}.csv", index=False)
    df_boot.to_csv(RESULTS_DIR / f"bootstrap_ci_{tag}.csv", index=False)
    (RESULTS_DIR / f"summary_{tag}.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    adj = adjudicate_r10_results(df_lines, df_boot, summary, confirmatory_lines=confirmatory_lines)
    (RESULTS_DIR / f"adjudication_{tag}.json").write_text(json.dumps(adj, indent=2), encoding="utf-8")
    (RESULTS_DIR / f"REPORT_{tag}.md").write_text(
        render_r10_markdown_report(df_lines, df_boot, summary, adj, tag), encoding="utf-8"
    )
    print(f"Saved Round 10 results (*_{tag}.*) in {summary['elapsed_sec']}s", flush=True)


if __name__ == "__main__":
    main()
