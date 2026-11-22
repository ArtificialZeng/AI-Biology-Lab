#!/usr/bin/env python3
"""
Round 11 Group D Predictor Suite (`FCT` + Target-Agnostic Constant Baseline `E1_ConstMed`
+ Mismatched Source Controls + R9 Anchors).

Strictly imports frozen predictor implementations from `round10/src/fct_pipeline.py`
without modifying any Round 10 file. Zero target perturbation leakage.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple

import numpy as np

R10_SRC_DIR = Path(__file__).resolve().parents[2] / "round10" / "src"
if str(R10_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(R10_SRC_DIR))

from fct_pipeline import (  # noqa: E402
    FctSourceState,
    SourceLineInputR10,
    build_source_input_r10,
    compute_control_covariance_r10,
    compute_source_transfer_and_fc_bar,
    fct_from_fc_bar,
    predict_group_c_suite,
    subsample_ctrl_r10,
)


def predict_group_d_suite(
    state: FctSourceState,
    target_ctrl_mean_libnorm: np.ndarray,
    query_target_names: Sequence[str],
    query_target_gene_indices: Optional[Sequence[int]] = None,
    target_ctrl_cov: Optional[np.ndarray] = None,
    target_ctrl_mean_libnorm_full: Optional[np.ndarray] = None,
    min_var: float = 1e-8,
    include_r9_anchors: bool = True,
    include_exploratory_e: bool = False,
) -> Tuple[Dict[str, np.ndarray], Tuple[np.ndarray, np.ndarray, np.ndarray], Dict[str, float]]:
    """
    Generate Round 11 predictors reusing frozen `predict_group_c_suite`:
      - Primary: `FCT` (`PREREG_A v4 §4.6` / `PREREG_C §3.1`, unchanged)
      - Preregistered Constant-Baseline Null (`H-D3` & `H-D4`): `E1_ConstMed` using scalar
        `c = median(max(mu_t, 0))` identically across all `G` genes inside `fct_from_fc_bar`
      - Mismatched Source Controls (`NegCtrl_FCT_MismatchedMean(s')`, 6 R9 or 8 R9+R10)
      - Anchors (`Transfer_Mean`, `Source_Grand_Mean`, `GlobalScale`, `TMVS`, `Hybrid_ProjShared`, etc.)
      - Optional `E2/E3/E5/E6/E7` diagnostic ablations when `include_exploratory_e=True`.
    """
    preds, fc_bundle = predict_group_c_suite(
        state=state,
        target_ctrl_mean_libnorm=target_ctrl_mean_libnorm,
        query_target_names=query_target_names,
        query_target_gene_indices=query_target_gene_indices,
        target_ctrl_cov=target_ctrl_cov,
        target_ctrl_mean_libnorm_full=target_ctrl_mean_libnorm_full,
        min_var=min_var,
        include_r9_anchors=include_r9_anchors,
    )
    pred_transfer, fc_bar, has_hits = fc_bundle
    mu_t = np.maximum(np.asarray(target_ctrl_mean_libnorm, dtype=np.float64), 0.0)
    G = mu_t.shape[0]

    c_med = float(np.median(mu_t))
    c_p10 = float(np.percentile(mu_t, 10.0))
    c_p90 = float(np.percentile(mu_t, 90.0))
    mu_e1 = np.full(G, c_med, dtype=np.float64)
    preds["E1_ConstMed"] = fct_from_fc_bar(mu_e1, fc_bar, pred_transfer, has_hits)

    const_stats: Dict[str, float] = {
        "mu_t_median": c_med,
        "mu_t_p10": c_p10,
        "mu_t_p90": c_p90,
    }

    if include_exploratory_e:
        mu_e3 = np.mean(
            [np.maximum(np.asarray(s.ctrl_mean_libnorm, dtype=np.float64), 0.0) for s in state.source_lines],
            axis=0,
        )
        mu_t_sorted = np.sort(mu_t, kind="mergesort")
        rank_src = np.argsort(np.argsort(mu_e3, kind="mergesort"), kind="mergesort")
        mu_e5 = mu_t_sorted[rank_src]

        fc_safe = np.maximum(fc_bar, 1e-12)
        pred_e6 = np.log(fc_safe)
        pred_e7 = fc_bar - 1.0
        if not np.all(has_hits):
            pred_e6[~has_hits] = pred_transfer[~has_hits]
            pred_e7[~has_hits] = pred_transfer[~has_hits]

        preds["E2_ConstP10"] = fct_from_fc_bar(np.full(G, c_p10, dtype=np.float64), fc_bar, pred_transfer, has_hits)
        preds["E2_ConstP90"] = fct_from_fc_bar(np.full(G, c_p90, dtype=np.float64), fc_bar, pred_transfer, has_hits)
        preds["E3_SrcMeanCtrl"] = fct_from_fc_bar(mu_e3, fc_bar, pred_transfer, has_hits)
        preds["E5_RankMatchedMu"] = fct_from_fc_bar(mu_e5, fc_bar, pred_transfer, has_hits)
        preds["E6_LogFCBar"] = pred_e6
        preds["E7_FCBarMinus1"] = pred_e7
        const_stats["mu_e3_median"] = float(np.median(mu_e3))

    return preds, fc_bundle, const_stats
