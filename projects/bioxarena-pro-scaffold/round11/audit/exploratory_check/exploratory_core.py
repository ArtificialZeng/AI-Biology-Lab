#!/usr/bin/env python3
"""
Shared read-only data loaders and exact cross-mask bootstrap helpers for the
Round 11 Exploratory Review Audit (`round11/audit/exploratory_check/`).

Strictly imports ONLY frozen Round 11 audit modules (`round11/audit/independent_r11_*.py`),
reads ONLY standardized expression caches (`round9/p2/shared/lines`, `round10/shared/lines`,
`round11/shared/lines`), and NEVER touches `exploratory/null_cone_theorem/*.py` or `round11/src/`.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

EXPLORATORY_CHECK_DIR = Path(__file__).resolve().parent
AUDIT_DIR = EXPLORATORY_CHECK_DIR.parent
if str(AUDIT_DIR) not in sys.path:
    sys.path.insert(0, str(AUDIT_DIR))

from independent_r11_loader import (  # noqa: E402
    R9_LINES_DIR,
    R9_LINE_NAMES,
    R10_LINES_DIR,
    R11_LINES_DIR,
    AuditLineData,
    CanonicalHGNC,
    enforce_r11_gate,
    load_line_dataset,
    load_r9_common_genes,
)
from independent_r11_metrics import (  # noqa: E402
    EPS,
    _extract_pert_gene_indices,
    _extract_sa8_masks,
    _moments,
    _moments_stats,
    _prepare_residual_streams,
    ceiling_from_halves,
    ceiling_single_half,
    disattenuate,
    orthogonalize_to_ref,
    summarize_boot,
)
from independent_r11_predictors import (  # noqa: E402
    AlignedCohort,
    align_source_and_target,
)

LINES_12: Tuple[Tuple[str, str], ...] = (
    ("nadig_hepg2", "R9"),
    ("nadig_jurkat", "R9"),
    ("nourreddine_kolf21j", "R9"),
    ("replogle_k562_essential", "R9"),
    ("replogle_rpe1", "R9"),
    ("vcc2025_h1", "R9"),
    ("dld1", "R10"),
    ("cd4t_rest", "R10"),
    ("calu3", "R10"),
    ("cd4t_stim48", "R11"),
    ("telohaec", "R11"),
    ("feng_ipsc", "R11"),
)

DEFAULT_AUDIT_SEED = 20261124


def load_12_lines() -> Dict[str, AuditLineData]:
    """Load all 12 standardized cell lines (excluding N=4 wtc11_i3neuron) after freeze check."""
    enforce_r11_gate(["cd4t_stim48", "telohaec", "feng_ipsc"], dev_only=False)
    hgnc = CanonicalHGNC.from_tsv()
    common_genes = load_r9_common_genes(hgnc=hgnc)
    lines: Dict[str, AuditLineData] = {}
    for ln, grp in LINES_12:
        if grp == "R9":
            p = R9_LINES_DIR / ln
            min_ov = 3000
        elif grp == "R10":
            p = R10_LINES_DIR / ln
            min_ov = 50
        elif grp == "R11":
            p = R11_LINES_DIR / ln
            min_ov = 50
        else:
            raise ValueError(f"Unknown group {grp!r}")
        lines[ln] = load_line_dataset(
            p,
            representation="log1p",
            aligned_genes=common_genes,
            hgnc=hgnc,
            min_gene_overlap=min_ov,
        )
    return lines


def build_r9_lolo_cohort(lines: Dict[str, AuditLineData], target_name: str) -> AlignedCohort:
    """Align `target_name` against R9 source lines (excluding `target_name` if it is in R9)."""
    src_names = [s for s in R9_LINE_NAMES if s != target_name]
    return align_source_and_target(
        source_lines=[lines[s] for s in src_names],
        target_data=lines[target_name],
        only_reproducible=True,
    )


def evaluate_cross_mask_split_minus_all(
    delta_obs: np.ndarray,
    pred: np.ndarray,
    t_mean: np.ndarray,
    target_cols: Sequence[int],
    half_a_obs: Optional[np.ndarray],
    half_b_obs: Optional[np.ndarray],
    half_valid: Optional[np.ndarray] = None,
    top_k: int = 50,
    n_boot: int = 4000,
    seed: int = DEFAULT_AUDIT_SEED,
) -> Dict[str, Any]:
    """
    Joint perturbation-level bootstrap contrast:
        `r_orth_cal(pred, resp_top50_split) - r_orth_cal(pred, all_nontarget)`
    with `resp_top50_split` using symmetric `(A->B + B->A) / 2` and single-half ceiling
    `sqrt(rho_AB)`, and `all_nontarget` using Spearman-Brown ceiling `sqrt(2*rho/(1+rho))`
    inside identical perturbation resample indices.
    """
    split_masks = _extract_sa8_masks(
        delta_obs=delta_obs,
        t_mean=t_mean,
        pred=pred,
        target_cols=target_cols,
        half_a_obs=half_a_obs,
        half_b_obs=half_b_obs,
        half_valid=half_valid,
        top_k=top_k,
    )
    all_masks = _extract_pert_gene_indices(
        delta_obs=delta_obs,
        t_mean=t_mean,
        pred=pred,
        target_cols=target_cols,
        mode="all_nontarget",
    )
    if not split_masks or not all_masks or half_a_obs is None or half_b_obs is None:
        out = summarize_boot(float("nan"), np.array([float("nan")]), 0, n_boot).to_dict("diff")
        out["n_eval_perts"] = 0
        return out

    map_split = {p: (ia, ib) for p, ia, ib in split_masks}
    map_all = {p: idx for p, idx in all_masks}
    common_p = [p for p, _, _ in split_masks if p in map_all]
    m = len(common_p)
    if m == 0:
        out = summarize_boot(float("nan"), np.array([float("nan")]), 0, n_boot).to_dict("diff")
        out["n_eval_perts"] = 0
        return out

    # Prepare split streams
    xb_ab, yp_ab, ha_ab, hb_ab = [], [], [], []
    xa_ba, yp_ba, ha_ba, hb_ba = [], [], [], []
    for p in common_p:
        ia, ib = map_split[p]
        ref_a = t_mean[p, ia].astype(np.float64)
        ha_a = orthogonalize_to_ref(half_a_obs[p, ia].astype(np.float64) - ref_a, ref_a)
        hb_a = orthogonalize_to_ref(half_b_obs[p, ia].astype(np.float64) - ref_a, ref_a)
        rp_a = orthogonalize_to_ref(pred[p, ia].astype(np.float64) - ref_a, ref_a)
        xb_ab.append(hb_a)
        yp_ab.append(rp_a)
        ha_ab.append(ha_a)
        hb_ab.append(hb_a)

        ref_b = t_mean[p, ib].astype(np.float64)
        ha_b = orthogonalize_to_ref(half_a_obs[p, ib].astype(np.float64) - ref_b, ref_b)
        hb_b = orthogonalize_to_ref(half_b_obs[p, ib].astype(np.float64) - ref_b, ref_b)
        rp_b = orthogonalize_to_ref(pred[p, ib].astype(np.float64) - ref_b, ref_b)
        xa_ba.append(ha_b)
        yp_ba.append(rp_b)
        ha_ba.append(ha_b)
        hb_ba.append(hb_b)

    # Prepare all_nontarget stream on identical common_p
    kept_all = [(p, map_all[p]) for p in common_p]
    x_all, y_all, ha_all, hb_all, any_half_all = _prepare_residual_streams(
        kept_all, delta_obs, pred, t_mean, half_a_obs, half_b_obs, half_valid, ortho=True
    )
    if not any_half_all:
        out = summarize_boot(float("nan"), np.array([float("nan")]), 0, n_boot).to_dict("diff")
        out["n_eval_perts"] = 0
        return out

    rng = np.random.default_rng(seed)
    draws = rng.integers(0, m, size=(n_boot, m))

    m_ab, mh_ab = _moments(xb_ab, yp_ab), _moments(ha_ab, hb_ab)
    m_ba, mh_ba = _moments(xa_ba, yp_ba), _moments(ha_ba, hb_ba)
    m_all, mh_all = _moments(x_all, y_all), _moments(ha_all, hb_all)

    pt_cal_ab = disattenuate(float(_moments_stats(*m_ab)[0]), ceiling_single_half(float(_moments_stats(*mh_ab)[0])))
    pt_cal_ba = disattenuate(float(_moments_stats(*m_ba)[0]), ceiling_single_half(float(_moments_stats(*mh_ba)[0])))
    pt_split = 0.5 * (pt_cal_ab + pt_cal_ba) if (np.isfinite(pt_cal_ab) and np.isfinite(pt_cal_ba)) else float("nan")
    pt_all = disattenuate(float(_moments_stats(*m_all)[0]), ceiling_from_halves(float(_moments_stats(*mh_all)[0])))
    pt_diff = pt_split - pt_all

    b_r_ab = _moments_stats(*m_ab, draws=draws)[0]
    b_rho_ab = _moments_stats(*mh_ab, draws=draws)[0]
    b_r_ba = _moments_stats(*m_ba, draws=draws)[0]
    b_rho_ba = _moments_stats(*mh_ba, draws=draws)[0]
    b_r_all = _moments_stats(*m_all, draws=draws)[0]
    b_rho_all = _moments_stats(*mh_all, draws=draws)[0]

    b_diff = np.empty(n_boot, dtype=np.float64)
    for bi in range(n_boot):
        c_ab = disattenuate(float(b_r_ab[bi]), ceiling_single_half(float(b_rho_ab[bi])))
        c_ba = disattenuate(float(b_r_ba[bi]), ceiling_single_half(float(b_rho_ba[bi])))
        c_all = disattenuate(float(b_r_all[bi]), ceiling_from_halves(float(b_rho_all[bi])))
        if np.isfinite(c_ab) and np.isfinite(c_ba) and np.isfinite(c_all):
            b_diff[bi] = 0.5 * (c_ab + c_ba) - c_all
        else:
            b_diff[bi] = float("nan")

    out = summarize_boot(pt_diff, b_diff, m, n_boot).to_dict("diff")
    out["n_eval_perts"] = m
    out["split_est"] = pt_split
    out["all_est"] = pt_all
    return out
