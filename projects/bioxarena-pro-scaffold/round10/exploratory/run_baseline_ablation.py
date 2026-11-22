#!/usr/bin/env python3
"""
R10-E Post-Hoc Exploratory Baseline Ablation for Fold-Change Transfer (FCT).

Strictly POST-HOC EXPLORATORY (事后探索性、非预注册).
Imports ONLY frozen modules from `round10/src/` without modifying any frozen file.

Evaluates FCT and baseline-ablated variants (`E1`, `E2_p10`, `E2_p90`, `E3`, `E5`,
`E6=log(FC_bar)`, `E7=FC_bar - 1`, and 100-perm null) across:
  - 6 Round 9 lines (strict LOLO, 5 sources per target line)
  - 4 Round 10 lines (`dld1`, `cd4t_rest` confirmatory; `calu3`, `wtc11_i3neuron` secondary; 6 R9 sources)
Under BOTH gene masks:
  - `resp_top50`
  - `all_nontarget` (S6 selection-free full non-target panel)
With joint perturbation-level bootstrap (`seed=20261122`, `n_boot=4000`, two-sided `97.5%` CIs)
plus paired difference `FCT - E1` (`97.5%` CI on the exact same bootstrap draws).
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

ROUND10_ROOT = Path(__file__).resolve().parents[1]
R10_SRC_DIR = ROUND10_ROOT / "src"
if str(R10_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(R10_SRC_DIR))

from data_loader_r10 import (
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
    compute_source_transfer_and_fc_bar,
    fct_from_fc_bar,
)
from metrics_r10 import (
    EPS,
    _extract_masked_idx_per_pert,
    compute_fct_perm100_null,
    compute_joint_paired_r_orth_diff_r10,
    compute_residual_suite_r10,
    disattenuated_r,
    ols_slope_1d,
    pearson_1d_strict,
)
from run_eval_r10 import select_eval_targets_r10

OUT_DIR = Path(__file__).resolve().parent


def build_ablation_predictors(
    state: FctSourceState,
    target_ctrl_mean_libnorm: np.ndarray,
    query_target_names: Sequence[str],
    min_var: float = 1e-8,
) -> Tuple[Dict[str, np.ndarray], Tuple[np.ndarray, np.ndarray, np.ndarray], Dict[str, float]]:
    """
    Construct FCT and post-hoc baseline ablation predictors strictly from `fc_bar`
    and unperturbed control summaries (`mu_t`, `mu_S`).
    """
    mu_t = np.maximum(np.asarray(target_ctrl_mean_libnorm, dtype=np.float64), 0.0)
    G = mu_t.shape[0]
    fc_bundle = compute_source_transfer_and_fc_bar(state, query_target_names, min_var=min_var)
    pred_transfer, fc_bar, has_hits = fc_bundle

    # Scalar constant baselines from target control mean distribution
    c_med = float(np.median(mu_t))
    c_p10 = float(np.percentile(mu_t, 10.0))
    c_p90 = float(np.percentile(mu_t, 90.0))
    mu_e1 = np.full(G, c_med, dtype=np.float64)
    mu_e2_p10 = np.full(G, c_p10, dtype=np.float64)
    mu_e2_p90 = np.full(G, c_p90, dtype=np.float64)

    # E3: arithmetic mean of source lines' libnorm control means (5 sources in R9 LOLO, 6 in R10)
    mu_e3 = np.mean(
        [np.maximum(np.asarray(s.ctrl_mean_libnorm, dtype=np.float64), 0.0) for s in state.source_lines],
        axis=0,
    )

    # E5: quantile-reassigned mu_t ordered by source mean rank (preserves mu_t marginal distribution, strips gene correspondence)
    mu_t_sorted = np.sort(mu_t, kind="mergesort")
    rank_src = np.argsort(np.argsort(mu_e3, kind="mergesort"), kind="mergesort")
    mu_e5 = mu_t_sorted[rank_src]

    # E6: large-baseline limit mu -> +inf  => log(FC_bar)
    fc_safe = np.maximum(fc_bar, 1e-12)
    pred_e6 = np.log(fc_safe)
    if not np.all(has_hits):
        pred_e6[~has_hits] = pred_transfer[~has_hits]

    # E7: small-baseline limit mu -> 0+   => FC_bar - 1
    pred_e7 = fc_bar - 1.0
    if not np.all(has_hits):
        pred_e7[~has_hits] = pred_transfer[~has_hits]

    preds: Dict[str, np.ndarray] = {
        "FCT": fct_from_fc_bar(mu_t, fc_bar, pred_transfer, has_hits),
        "E1_ConstMed": fct_from_fc_bar(mu_e1, fc_bar, pred_transfer, has_hits),
        "E2_ConstP10": fct_from_fc_bar(mu_e2_p10, fc_bar, pred_transfer, has_hits),
        "E2_ConstP90": fct_from_fc_bar(mu_e2_p90, fc_bar, pred_transfer, has_hits),
        "E3_SrcMeanCtrl": fct_from_fc_bar(mu_e3, fc_bar, pred_transfer, has_hits),
        "E5_RankMatchedMu": fct_from_fc_bar(mu_e5, fc_bar, pred_transfer, has_hits),
        "E6_LogFCBar": pred_e6,
        "E7_FCBarMinus1": pred_e7,
    }
    const_info = {
        "mu_t_p10": c_p10,
        "mu_t_median": c_med,
        "mu_t_p90": c_p90,
        "mu_e3_median": float(np.median(mu_e3)),
    }
    return preds, fc_bundle, const_info


def compute_fct_perm100_for_mode(
    state: FctSourceState,
    target_ctrl_mean_libnorm: np.ndarray,
    query_target_names: Sequence[str],
    delta_obs: np.ndarray,
    delta_transfer: np.ndarray,
    target_gene_idx: Sequence[int],
    ceil_r: float,
    ceil_r_orth: float,
    fc_bundle: Tuple[np.ndarray, np.ndarray, np.ndarray],
    mode: str,
    n_perm: int = 100,
    top_k: int = 50,
) -> Dict[str, Any]:
    """
    Fast vectorized 100-perm FCT null distribution for both `resp_top50` and `all_nontarget`,
    strictly matching `compute_fct_perm100_null` (`default_rng(k).permutation(G)`).
    """
    if mode == "resp_top50":
        return compute_fct_perm100_null(
            state=state,
            target_ctrl_mean_libnorm=target_ctrl_mean_libnorm,
            query_target_names=query_target_names,
            delta_obs=delta_obs,
            delta_transfer=delta_transfer,
            target_gene_idx=target_gene_idx,
            ceil_r_top50=ceil_r,
            ceil_r_orth_top50=ceil_r_orth,
            fc_bundle=fc_bundle,
            n_perm=n_perm,
            top_k=top_k,
        )

    mu_t = np.maximum(np.asarray(target_ctrl_mean_libnorm, dtype=np.float64), 0.0)
    G = mu_t.shape[0]
    pred_transfer, fc_bar, has_hits = fc_bundle
    r_obs = np.asarray(delta_obs, dtype=np.float64) - np.asarray(delta_transfer, dtype=np.float64)
    rows, idx_list = _extract_masked_idx_per_pert(
        delta_obs, target_gene_idx, mode=mode, top_k=top_k
    )
    if not rows:
        return {"n_perm": n_perm, "r_orth_cal": {}, "r_cal": {}, "slope_ols": {}}

    rows_arr = np.asarray(rows, dtype=np.int64)
    m = len(rows)
    M = np.zeros((m, G), dtype=np.float64)
    B = np.zeros((m, G), dtype=np.float64)
    for ii, (p, idx) in enumerate(zip(rows, idx_list)):
        M[ii, idx] = 1.0
        rv = delta_transfer[p, idx].astype(np.float64)
        nrm = float(np.linalg.norm(rv))
        if nrm > EPS:
            B[ii, idx] = rv / nrm

    x_mat_raw = r_obs[rows_arr] * M
    proj_x = np.sum(x_mat_raw * B, axis=1, keepdims=True)
    x_mat_orth = x_mat_raw - proj_x * B

    mask_bool = M > 0.5
    x_raw = x_mat_raw[mask_bool]
    x_orth = x_mat_orth[mask_bool]
    n_tot = float(x_raw.size)

    r_orth_cal_vals = np.empty(n_perm, dtype=np.float64)
    r_cal_vals = np.empty(n_perm, dtype=np.float64)
    slope_vals = np.empty(n_perm, dtype=np.float64)

    fc_bar_s = fc_bar[rows_arr]
    p_trans_s = pred_transfer[rows_arr]
    hits_s = has_hits[rows_arr]
    t_mean_s = delta_transfer[rows_arr]

    for k in range(n_perm):
        perm = np.random.default_rng(k).permutation(G)
        pred_k = fct_from_fc_bar(mu_t[perm], fc_bar_s, p_trans_s, hits_s)
        res_k = (pred_k - t_mean_s) * M
        proj_y = np.sum(res_k * B, axis=1, keepdims=True)
        res_k_orth = res_k - proj_y * B

        y_raw = res_k[mask_bool]
        y_orth = res_k_orth[mask_bool]
        slope_vals[k] = ols_slope_1d(x_raw, y_raw)
        r_cal_vals[k] = disattenuated_r(pearson_1d_strict(x_raw, y_raw), ceil_r)
        r_orth_cal_vals[k] = disattenuated_r(pearson_1d_strict(x_orth, y_orth), ceil_r_orth)

    def _pct(arr: np.ndarray) -> Dict[str, float]:
        fin = arr[np.isfinite(arr)]
        if fin.size == 0:
            return {k: float("nan") for k in ("p2_5", "median", "p97_5", "p98_75", "max", "mean")}
        p25, p50, p975, p9875 = np.percentile(fin, [2.5, 50.0, 97.5, 98.75])
        return {
            "p2_5": float(p25),
            "median": float(p50),
            "p97_5": float(p975),
            "p98_75": float(p9875),
            "max": float(np.max(fin)),
            "mean": float(np.mean(fin)),
        }

    return {
        "n_perm": int(n_perm),
        "r_orth_cal": _pct(r_orth_cal_vals),
        "r_cal": _pct(r_cal_vals),
        "slope_ols": _pct(slope_vals),
    }


def evaluate_line_ablation(
    target_ds: LineDatasetR10,
    source_datasets: Dict[str, LineDatasetR10],
    cohort_label: str,
    n_boot: int = 4000,
    n_perm: int = 100,
    seed: int = 20261122,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """Run FCT baseline ablation (E1..E7 + perm100 + paired diff FCT - E1) on one target line."""
    src_names = sorted(source_datasets.keys())
    src_gene_sets = [set(source_datasets[s].genes) for s in src_names]
    shared_genes_set = set(target_ds.genes).intersection(*src_gene_sets)
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
                only_reproducible=True,
                compute_cov=False,
                seed=42,
            )
        )

    state = FctSourceState.fit(src_inputs, mask_self=True)
    row_idx, q_names, q_gidx = select_eval_targets_r10(t_sub_ds, src_inputs, only_reproducible=True)

    preds, fc_bundle, const_info = build_ablation_predictors(
        state=state,
        target_ctrl_mean_libnorm=t_sub_ds.ctrl_mean_libnorm,
        query_target_names=q_names,
    )
    true_deltas = t_sub_ds.deltas[row_idx]
    transfer_ref = fc_bundle[0]
    ha_sub = t_sub_ds.half_a_deltas[row_idx]
    hb_sub = t_sub_ds.half_b_deltas[row_idx]
    hv_sub = t_sub_ds.half_valid[row_idx]

    ablation_rows: List[Dict[str, Any]] = []
    diff_rows: List[Dict[str, Any]] = []

    ceil_by_mode: Dict[str, Tuple[float, float]] = {}
    for m_name, p_mat in preds.items():
        res_suite = compute_residual_suite_r10(
            delta_obs=true_deltas,
            delta_pred=p_mat,
            delta_transfer=transfer_ref,
            target_gene_idx=q_gidx,
            delta_half_a=ha_sub,
            delta_half_b=hb_sub,
            half_valid=hv_sub,
            top_k=50,
            n_boot=n_boot,
            seed=seed,
        )
        for mode in ("resp_top50", "all_nontarget"):
            md = res_suite[mode]
            ro = md.get("r_orth_cal", {})
            rc = md.get("r_cal", {})
            sl = md.get("slope_ols", {})
            ceil_by_mode[mode] = (
                float(md.get("ceiling_r", float("nan"))),
                float(md.get("ceiling_r_orth", float("nan"))),
            )
            ablation_rows.append(
                {
                    "cohort": cohort_label,
                    "target_line": target_ds.name,
                    "mask": mode,
                    "method": m_name,
                    "n_genes_eval": len(aligned_genes),
                    "n_eval_targets": int(len(q_names)),
                    "n_sources": len(src_names),
                    "r_orth_cal": ro.get("estimate", float("nan")),
                    "r_orth_cal_ci95_lo": ro.get("ci95_lo", float("nan")),
                    "r_orth_cal_ci95_hi": ro.get("ci95_hi", float("nan")),
                    "r_orth_cal_ci975_lo": ro.get("ci975_lo", float("nan")),
                    "r_orth_cal_ci975_hi": ro.get("ci975_hi", float("nan")),
                    "r_orth_raw": md.get("r_orth_raw", {}).get("estimate", float("nan")),
                    "r_cal": rc.get("estimate", float("nan")),
                    "slope_ols": sl.get("estimate", float("nan")),
                    "slope_ols_ci975_lo": sl.get("ci975_lo", float("nan")),
                    "slope_ols_ci975_hi": sl.get("ci975_hi", float("nan")),
                    "rho_half_ab_orth": md.get("rho_half_ab_orth", float("nan")),
                    "ceiling_r_orth": md.get("ceiling_r_orth", float("nan")),
                    **const_info,
                }
            )

    # 100-perm null median/p97.5 & paired difference FCT - E1 under BOTH masks
    for mode in ("resp_top50", "all_nontarget"):
        cr, cro = ceil_by_mode[mode]
        p100 = compute_fct_perm100_for_mode(
            state=state,
            target_ctrl_mean_libnorm=t_sub_ds.ctrl_mean_libnorm,
            query_target_names=q_names,
            delta_obs=true_deltas,
            delta_transfer=transfer_ref,
            target_gene_idx=q_gidx,
            ceil_r=cr,
            ceil_r_orth=cro,
            fc_bundle=fc_bundle,
            mode=mode,
            n_perm=n_perm,
            top_k=50,
        )
        pro = p100.get("r_orth_cal", {})
        psl = p100.get("slope_ols", {})
        ablation_rows.append(
            {
                "cohort": cohort_label,
                "target_line": target_ds.name,
                "mask": mode,
                "method": "Perm100_Median",
                "n_genes_eval": len(aligned_genes),
                "n_eval_targets": int(len(q_names)),
                "n_sources": len(src_names),
                "r_orth_cal": pro.get("median", float("nan")),
                "r_orth_cal_ci95_lo": pro.get("p2_5", float("nan")),
                "r_orth_cal_ci95_hi": pro.get("p97_5", float("nan")),
                "r_orth_cal_ci975_lo": pro.get("p2_5", float("nan")),
                "r_orth_cal_ci975_hi": pro.get("p97_5", float("nan")),
                "r_orth_raw": float("nan"),
                "r_cal": p100.get("r_cal", {}).get("median", float("nan")),
                "slope_ols": psl.get("median", float("nan")),
                "slope_ols_ci975_lo": psl.get("p2_5", float("nan")),
                "slope_ols_ci975_hi": psl.get("p97_5", float("nan")),
                "rho_half_ab_orth": float("nan"),
                "ceiling_r_orth": cro,
                **const_info,
            }
        )

        # Joint perturbation-level paired diff: FCT - E1_ConstMed
        d_fct_e1 = compute_joint_paired_r_orth_diff_r10(
            delta_obs=true_deltas,
            delta_pred_a=preds["FCT"],
            delta_pred_b=preds["E1_ConstMed"],
            delta_transfer=transfer_ref,
            target_gene_idx=q_gidx,
            delta_half_a=ha_sub,
            delta_half_b=hb_sub,
            half_valid=hv_sub,
            mode=mode,
            top_k=50,
            n_boot=n_boot,
            seed=seed,
        )
        fct_est = next(
            r["r_orth_cal"]
            for r in ablation_rows
            if r["target_line"] == target_ds.name and r["mask"] == mode and r["method"] == "FCT"
        )
        e1_est = next(
            r["r_orth_cal"]
            for r in ablation_rows
            if r["target_line"] == target_ds.name and r["mask"] == mode and r["method"] == "E1_ConstMed"
        )
        diff_rows.append(
            {
                "cohort": cohort_label,
                "target_line": target_ds.name,
                "mask": mode,
                "n_eval_targets": int(len(q_names)),
                "fct_r_orth_cal": fct_est,
                "e1_r_orth_cal": e1_est,
                "e1_over_fct_ratio": (e1_est / fct_est) if (np.isfinite(fct_est) and abs(fct_est) > 1e-12) else float("nan"),
                "diff_fct_minus_e1": d_fct_e1["estimate"],
                "diff_ci95_lo": d_fct_e1["ci95_lo"],
                "diff_ci95_hi": d_fct_e1["ci95_hi"],
                "diff_ci975_lo": d_fct_e1["ci975_lo"],
                "diff_ci975_hi": d_fct_e1["ci975_hi"],
                "ci975_excludes_zero": bool(
                    np.isfinite(d_fct_e1["ci975_lo"])
                    and np.isfinite(d_fct_e1["ci975_hi"])
                    and (d_fct_e1["ci975_lo"] > 0.0 or d_fct_e1["ci975_hi"] < 0.0)
                ),
                "fct_gt_e1_ci975": bool(np.isfinite(d_fct_e1["ci975_lo"]) and d_fct_e1["ci975_lo"] > 0.0),
                "perm100_median_r_orth_cal": pro.get("median", float("nan")),
                "perm100_p975_r_orth_cal": pro.get("p97_5", float("nan")),
                "e1_gt_perm100_p975": bool(
                    np.isfinite(e1_est) and np.isfinite(pro.get("p97_5", float("nan"))) and e1_est > pro["p97_5"]
                ),
            }
        )

    return ablation_rows, diff_rows


def render_markdown_report(df_ab: pd.DataFrame, df_diff: pd.DataFrame, elapsed_sec: float) -> str:
    md: List[str] = [
        "# Round 10-E Post-Hoc Exploratory Baseline Ablation Report (`round10/exploratory`)",
        "",
        "> **重要声明：事后探索性（post hoc，非预注册）分析。** 本报告仅作为下一轮假设生成依据，不改变 [`R10_VERDICT.md`](../status/R10_VERDICT.md) 的确认性判定（`H-C1: 2/2`, `H-C2: 0/2`, 主要科学论断 `未确认`, `H-C3: 12/12 非劣`）。",
        "",
        f"- **冻结代码基线**：复用 `round10/src/`（`f8d2ebfc`），未改动任何冻结清单文件；联合 perturbation bootstrap `default_rng(20261122)` x `4,000`，双侧 `97.5%` CI；置换 `default_rng(k)` `k=0..99`；总运行耗时 `{elapsed_sec:.2f}s`。",
        "- **对比模型（事后探索性 post hoc）**：",
        "  - `FCT`（真实目标线对照 libnorm 均值 $\\mu_t$，对照）",
        "  - `E1_ConstMed`：全基因统一常数基线 $c = \\operatorname{median}(\\mu_t)$（完全不含目标基因差异信息）",
        "  - `E2_ConstP10` / `E2_ConstP90`：全基因常数取 $\\mu_t$ 的 10% / 90% 分位数",
        "  - `E3_SrcMeanCtrl`：源线对照 libnorm 均值平均（R10 取 6 条源线平均，R9 LOLO 取 5 条源线平均）",
        "  - `E5_RankMatchedMu`：保留 $\\mu_t$ 经验分布、但按源线平均排序重分配（剥离目标线基因特异排布）",
        "  - `E6_LogFCBar`：$\\mu_t \\to +\\infty$ 极限 $\\ln(\\overline{\\mathrm{FC}})$",
        "  - `E7_FCBarMinus1`：$\\mu_t \\to 0^+$ 极限 $\\overline{\\mathrm{FC}} - 1$",
        "  - `Perm100_Median`：100 次基因标签置换目标均值 $\\mu_t$ 的中位数（方括号列为 `[p2.5, p97.5]`）",
        "",
    ]

    for mode, title in [
        ("resp_top50", "## 表 1（事后探索性 post hoc）：`resp_top50` 掩码下基线消融 `r_orth_cal [97.5% CI]`"),
        ("all_nontarget", "## 表 2（事后探索性 post hoc）：`all_nontarget` (`S6` 无选基因偏倚全掩码) 基线消融 `r_orth_cal [97.5% CI]`"),
    ]:
        md.extend([title, ""])
        md.append("| Cohort | Target Line | N | FCT | E1 (Med) | E2 (P10) | E2 (P90) | E3 (SrcMean) | E5 (Rank) | E6 log(FC̄) | E7 (FC̄−1) | Perm100 Med |")
        md.append("|---|---|---:|---|---|---|---|---|---|---|---|---|")
        sub = df_ab[df_ab["mask"] == mode]
        line_order = [
            "dld1", "cd4t_rest", "calu3", "wtc11_i3neuron",
            "nadig_hepg2", "nadig_jurkat", "vcc2025_h1",
            "replogle_k562_essential", "replogle_rpe1", "nourreddine_kolf21j",
        ]
        cols = [
            "FCT", "E1_ConstMed", "E2_ConstP10", "E2_ConstP90",
            "E3_SrcMeanCtrl", "E5_RankMatchedMu", "E6_LogFCBar", "E7_FCBarMinus1", "Perm100_Median",
        ]
        for ln in line_order:
            ls = sub[sub["target_line"] == ln]
            if ls.empty:
                continue
            cohort = str(ls.iloc[0]["cohort"])
            n_t = int(ls.iloc[0]["n_eval_targets"])
            by_m = {r["method"]: r for _, r in ls.iterrows()}
            cells = []
            for m in cols:
                r = by_m.get(m)
                if r is None:
                    cells.append("—")
                elif m == "Perm100_Median":
                    cells.append(f"`{r['r_orth_cal']:+.4f}` (`p97.5={r['r_orth_cal_ci975_hi']:+.4f}`)")
                else:
                    cells.append(f"`{r['r_orth_cal']:+.4f}` [`{r['r_orth_cal_ci975_lo']:+.4f}`, `{r['r_orth_cal_ci975_hi']:+.4f}`]")
            md.append(f"| {cohort} | **{ln}** | {n_t} | " + " | ".join(cells) + " |")
        md.append("")

    md.extend([
        "## 表 3（事后探索性 post hoc）：同一组联合 Bootstrap 下配对差 `FCT − E1` (`97.5% CI`)",
        "",
        "| Cohort | Target Line | Mask | N | FCT | E1 (ConstMed) | E1 / FCT | FCT − E1 | 97.5% CI | FCT > E1 (`ci975_lo > 0`)? |",
        "|---|---|---|---:|---:|---:|---:|---:|---|---|",
    ])
    for _, r in df_diff.iterrows():
        ratio_s = f"{r['e1_over_fct_ratio']*100:.1f}%" if np.isfinite(r["e1_over_fct_ratio"]) else "nan"
        md.append(
            f"| {r['cohort']} | **{r['target_line']}** | `{r['mask']}` | {int(r['n_eval_targets'])} | "
            f"`{r['fct_r_orth_cal']:+.4f}` | `{r['e1_r_orth_cal']:+.4f}` | {ratio_s} | "
            f"`{r['diff_fct_minus_e1']:+.4f}` | [`{r['diff_ci975_lo']:+.4f}`, `{r['diff_ci975_hi']:+.4f}`] | "
            f"{'True' if r['fct_gt_e1_ci975'] else 'False'} |"
        )
    md.append("")
    return "\n".join(md) + "\n"


def main() -> None:
    t0 = time.time()
    p = argparse.ArgumentParser(description="R10-E Post-Hoc Baseline Ablation.")
    p.add_argument("--n-boot", type=int, default=4000)
    p.add_argument("--n-perm", type=int, default=100)
    p.add_argument("--seed", type=int, default=20261122)
    args = p.parse_args()

    check_prereg_c_gate(R10_STATUS_DIR)
    hgnc = HGNCBridge.from_tsv(DEFAULT_HGNC_FILE)
    common_r9 = load_r9_common_genes(hgnc_bridge=hgnc, common_genes_file=R9_COMMON_GENES_FILE)

    r9_names = list_standardized_lines(R9_LINES_DIR)
    r10_names = list_standardized_lines(R10_LINES_DIR)
    r9_ds = {
        ln: load_standardized_line_r10(R9_LINES_DIR / ln, "log1p", common_genes=common_r9, hgnc_bridge=hgnc)
        for ln in r9_names
    }
    r10_ds = {
        ln: load_standardized_line_r10(R10_LINES_DIR / ln, "log1p", common_genes=common_r9, hgnc_bridge=hgnc)
        for ln in r10_names
    }

    r10_conf_set = {"dld1", "cd4t_rest"}
    all_ab: List[Dict[str, Any]] = []
    all_diff: List[Dict[str, Any]] = []

    for ln in r10_names:
        lbl = "R10_Confirmatory" if ln in r10_conf_set else "R10_Secondary"
        ab, df_r = evaluate_line_ablation(
            r10_ds[ln], r9_ds, cohort_label=lbl, n_boot=args.n_boot, n_perm=args.n_perm, seed=args.seed
        )
        all_ab.extend(ab)
        all_diff.extend(df_r)

    for ln in r9_names:
        src_lolo = {k: v for k, v in r9_ds.items() if k != ln}
        ab, df_r = evaluate_line_ablation(
            r9_ds[ln], src_lolo, cohort_label="R9_LOLO", n_boot=args.n_boot, n_perm=args.n_perm, seed=args.seed
        )
        all_ab.extend(ab)
        all_diff.extend(df_r)

    elapsed = round(time.time() - t0, 2)
    df_ab = pd.DataFrame(all_ab)
    df_diff = pd.DataFrame(all_diff)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df_ab.to_csv(OUT_DIR / "baseline_ablation_metrics.csv", index=False)
    df_diff.to_csv(OUT_DIR / "baseline_ablation_paired_fct_vs_e1.csv", index=False)
    (OUT_DIR / "BASELINE_ABLATION.md").write_text(
        render_markdown_report(df_ab, df_diff, elapsed), encoding="utf-8"
    )
    print(f"Saved R10-E baseline ablation outputs to {OUT_DIR} in {elapsed}s", flush=True)


if __name__ == "__main__":
    main()
