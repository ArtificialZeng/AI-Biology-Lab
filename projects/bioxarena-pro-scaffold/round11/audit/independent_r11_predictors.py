#!/usr/bin/env python3
"""
Independent Round 11 Group D Predictors, Constant/Rank/Limit Nulls (`E1..E7`),
Mismatched Controls, & SA-3 R9 Anchors (A-T3, PREREG_D v2 §3 & §5).

Strict zero target-perturbation leakage:
  Predictors receive ONLY:
    - source line perturbations & control covariances / control means
    - unperturbed target control mean (`ctrl_mean_libnorm`, optional `_full`),
      unperturbed target control covariance (`C_t`), and target gene symbols.
  - Primary: `fct` (`PREREG_D §3.1`, unchanged from `PREREG_A v4 §4.6` / `PREREG_C §3.1`)
  - Constant target-baseline-free null (`H-D3`, `H-D4`): `e1_const_med` (`c = median(max(mu_t, 0))`)
  - Mismatched controls (`H-D2`, `H-D5`): `fct_mismatch::<s'>` for each source `s'`
  - Exploratory ablations (`SA-5`): `e2_const_p10`, `e2_const_p90`, `e3_src_mean_ctrl`,
    `e5_rank_matched_mu`, `e6_log_fc_bar`, `e7_fc_minus_1`
  - SA-2 & SA-3 anchors (`1ca9f60f`): `transfer_mean_global_scale`, `transfer_mean_var_scaled`,
    `hybrid`, `hybrid_on_tmvs`, `cipher_linear`, `source_grand_mean`, `fct_full_nt`.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

AUDIT_DIR = Path(__file__).resolve().parent
if str(AUDIT_DIR) not in sys.path:
    sys.path.insert(0, str(AUDIT_DIR))

from independent_r11_loader import AuditLineData

EPS = 1e-12


@dataclass
class AlignedCohort:
    """Aligned source pool + target metadata ready for zero-shot Round 11 evaluation."""

    target_line: str
    eval_genes: List[str]
    eval_targets: List[str]
    target_cols: np.ndarray
    delta_obs: np.ndarray
    n_cells: Optional[np.ndarray]
    kd_depth_log1p: np.ndarray
    half_a_obs: Optional[np.ndarray]
    half_b_obs: Optional[np.ndarray]
    half_valid: Optional[np.ndarray]
    tgt_ctrl_cells: np.ndarray
    tgt_ctrl_mean_libnorm: np.ndarray
    tgt_ctrl_mean_libnorm_full: Optional[np.ndarray]
    source_names: List[str]
    source_deltas_by_line: Dict[str, Dict[str, np.ndarray]]
    source_ctrl_cells: Dict[str, np.ndarray]
    source_ctrl_mean_libnorm: Dict[str, np.ndarray]


def sample_covariance(ctrl_cells: np.ndarray, shrink: float = 0.0) -> np.ndarray:
    """Unbiased sample covariance over control cells with optional diagonal shrinkage."""
    x = np.asarray(ctrl_cells, dtype=np.float64)
    n, g = x.shape
    if n < 2:
        return np.eye(g, dtype=np.float64)
    xc = x - np.mean(x, axis=0, keepdims=True)
    cov = (xc.T @ xc) / float(n - 1)
    cov = 0.5 * (cov + cov.T)
    if shrink > 0.0:
        cov = (1.0 - shrink) * cov + shrink * np.diag(np.diag(cov))
    return cov


def deterministic_subsample(ctrl_cells: np.ndarray, max_cells: Optional[int], seed: int = 42) -> np.ndarray:
    """Subsample rows deterministically (`seed=42`, sorted indices)."""
    if max_cells is None or max_cells <= 0 or ctrl_cells.shape[0] <= max_cells:
        return ctrl_cells
    rng = np.random.default_rng(seed)
    chosen = np.sort(rng.choice(ctrl_cells.shape[0], size=max_cells, replace=False))
    return ctrl_cells[chosen]


def _align_line_cols(line_genes: Sequence[str], target_genes: Sequence[str]) -> np.ndarray:
    pos = {g: i for i, g in enumerate(line_genes)}
    missing = [g for g in target_genes if g not in pos]
    if missing:
        raise ValueError(f"Missing {len(missing)} aligned target genes (e.g. {missing[:5]})")
    return np.asarray([pos[g] for g in target_genes], dtype=np.int64)


def align_source_and_target(
    source_lines: Sequence[AuditLineData],
    target_data: AuditLineData,
    only_reproducible: bool = True,
) -> AlignedCohort:
    """
    Construct G_eval(t) and P5-eligible shared perturbation targets:
      - Target perturbations must satisfy `reproducible` (when `only_reproducible=True`)
      - At least 1 source line must have that perturbation (and satisfy `reproducible` when `only_reproducible=True`).
      - Target self-gene g_p must resolve inside G_eval(t).
    """
    eval_genes = list(target_data.genes)
    src_deltas: Dict[str, Dict[str, np.ndarray]] = {}
    src_ctrl: Dict[str, np.ndarray] = {}
    src_mu_lib: Dict[str, np.ndarray] = {}
    src_avail_targets = set()

    for s in source_lines:
        cols = _align_line_cols(s.genes, eval_genes)
        src_ctrl[s.name] = np.asarray(s.ctrl_cells[:, cols], dtype=np.float64)
        src_mu_lib[s.name] = np.asarray(s.ctrl_mean_libnorm[cols], dtype=np.float64)
        s_map: Dict[str, np.ndarray] = {}
        for i, tgt in enumerate(s.target_names):
            if only_reproducible and s.repr_mask is not None and not bool(s.repr_mask[i]):
                continue
            s_map[tgt] = np.asarray(s.deltas[i, cols], dtype=np.float64)
            src_avail_targets.add(tgt)
        src_deltas[s.name] = s_map

    kept_rows: List[int] = []
    for i, tgt in enumerate(target_data.target_names):
        if int(target_data.target_indices[i]) < 0:
            continue
        if only_reproducible and target_data.repr_mask is not None and not bool(target_data.repr_mask[i]):
            continue
        if tgt not in src_avail_targets:
            continue
        kept_rows.append(i)

    if not kept_rows:
        raise RuntimeError(f"Zero shared eligible targets for {target_data.name} (only_reproducible={only_reproducible}).")

    idx_arr = np.asarray(kept_rows, dtype=np.int64)
    eval_targets = [target_data.target_names[i] for i in kept_rows]

    return AlignedCohort(
        target_line=target_data.name,
        eval_genes=eval_genes,
        eval_targets=eval_targets,
        target_cols=np.asarray(target_data.target_indices[idx_arr], dtype=np.int64),
        delta_obs=np.asarray(target_data.deltas[idx_arr], dtype=np.float64),
        n_cells=np.asarray(target_data.n_cells[idx_arr], dtype=np.float64) if target_data.n_cells is not None else None,
        kd_depth_log1p=np.asarray(target_data.kd_depth_log1p[idx_arr], dtype=np.float64),
        half_a_obs=np.asarray(target_data.half_a_deltas[idx_arr], dtype=np.float64) if target_data.half_a_deltas is not None else None,
        half_b_obs=np.asarray(target_data.half_b_deltas[idx_arr], dtype=np.float64) if target_data.half_b_deltas is not None else None,
        half_valid=np.asarray(target_data.half_valid[idx_arr], dtype=bool) if target_data.half_valid is not None else None,
        tgt_ctrl_cells=np.asarray(target_data.ctrl_cells, dtype=np.float64),
        tgt_ctrl_mean_libnorm=np.asarray(target_data.ctrl_mean_libnorm, dtype=np.float64),
        tgt_ctrl_mean_libnorm_full=(
            np.asarray(target_data.ctrl_mean_libnorm_full, dtype=np.float64)
            if target_data.ctrl_mean_libnorm_full is not None
            else None
        ),
        source_names=[s.name for s in source_lines],
        source_deltas_by_line=src_deltas,
        source_ctrl_cells=src_ctrl,
        source_ctrl_mean_libnorm=src_mu_lib,
    )


def compute_source_fc_bar(
    eval_targets: Sequence[str],
    source_names: Sequence[str],
    source_deltas_by_line: Dict[str, Dict[str, np.ndarray]],
    source_ctrl_mean_libnorm: Dict[str, np.ndarray],
    floor: float = 1e-8,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Independently implement PREREG_D §3.1 geometric-mean fold change + Transfer_Mean + Source_Grand_Mean:
      mu_{s,g} = max(ctrl_mean_libnorm_s, 1e-8)
      p_{s,g}  = max(expm1(log1p(mu_{s,g}) + delta_{s,p,g}), 1e-8)
      FC_bar_{p,g} = exp(mean_{s in S_p} ln(p_{s,g} / mu_{s,g}))
    """
    n_t = len(eval_targets)
    any_s = next(iter(source_ctrl_mean_libnorm.values()))
    n_g = any_s.shape[0]
    safe_mu = {s: np.maximum(np.asarray(source_ctrl_mean_libnorm[s], dtype=np.float64), floor) for s in source_names}
    safe_log1p_mu = {s: np.log1p(safe_mu[s]) for s in source_names}

    fc_bar = np.ones((n_t, n_g), dtype=np.float64)
    t_mean = np.zeros((n_t, n_g), dtype=np.float64)

    all_src_rows: List[np.ndarray] = []
    for s in source_names:
        for d in source_deltas_by_line[s].values():
            all_src_rows.append(np.asarray(d, dtype=np.float64))
    grand_mean = np.mean(np.stack(all_src_rows, axis=0), axis=0) if all_src_rows else np.zeros(n_g, dtype=np.float64)

    for i, tgt in enumerate(eval_targets):
        avail = [s for s in source_names if tgt in source_deltas_by_line[s]]
        if not avail:
            raise ValueError(f"Evaluated target {tgt!r} missing from all sources.")
        t_mean[i] = np.mean([source_deltas_by_line[s][tgt] for s in avail], axis=0)
        log_fcs = []
        for s in avail:
            d = np.asarray(source_deltas_by_line[s][tgt], dtype=np.float64)
            pert = np.maximum(np.expm1(safe_log1p_mu[s] + d), floor)
            log_fcs.append(np.log(pert / safe_mu[s]))
        fc_bar[i] = np.exp(np.mean(log_fcs, axis=0))

    return fc_bar, t_mean, grand_mean


def fct_from_fc_bar(fc_bar: np.ndarray, ctrl_mean_libnorm: np.ndarray) -> np.ndarray:
    """PREREG_D §3.1: delta_FCT = log1p(max(mu, 0) * FC_bar) - log1p(max(mu, 0))."""
    mu = np.maximum(np.asarray(ctrl_mean_libnorm, dtype=np.float64), 0.0)
    return np.log1p(mu[None, :] * fc_bar) - np.log1p(mu)[None, :]


def _unit_col_vec(cov: np.ndarray, g_col: int, min_var: float = 1e-8) -> np.ndarray:
    n_g = cov.shape[0]
    if g_col < 0 or g_col >= n_g:
        return np.zeros(n_g, dtype=np.float64)
    denom = max(float(cov[g_col, g_col]), min_var)
    return cov[:, g_col] / denom


def fit_3d_shared_subspace(
    source_deltas_by_line: Dict[str, Dict[str, np.ndarray]],
    n_genes: int,
) -> np.ndarray:
    """Frozen Round 9 (`1ca9f60f`) 3D shared subspace basis span{1, b_S, a_S}."""
    rows = [np.asarray(d, dtype=np.float64) for s_map in source_deltas_by_line.values() for d in s_map.values()]
    if not rows:
        return np.zeros((n_genes, 0), dtype=np.float64)
    mat = np.stack(rows, axis=0)
    ones = np.ones(n_genes, dtype=np.float64)
    grand_mean = np.mean(mat, axis=0)

    basis_cols: List[np.ndarray] = []
    for c in (ones, grand_mean):
        v = c.astype(np.float64).copy()
        for q in basis_cols:
            v = v - float(np.dot(q, v)) * q
        nrm = float(np.linalg.norm(v))
        if nrm > 1e-10:
            basis_cols.append(v / nrm)

    q12 = np.column_stack(basis_cols) if basis_cols else np.zeros((n_genes, 0), dtype=np.float64)
    res = mat - (mat @ q12) @ q12.T if q12.shape[1] > 0 else mat.copy()
    row_norms = np.linalg.norm(res, axis=1)
    keep = row_norms > 1e-10
    if int(np.sum(keep)) >= 1:
        res_normed = res[keep] / row_norms[keep, None]
        _, _, vt = np.linalg.svd(res_normed, full_matrices=False)
        a_s = vt[0].astype(np.float64).copy()
        for q in basis_cols:
            a_s = a_s - float(np.dot(q, a_s)) * q
        a_nrm = float(np.linalg.norm(a_s))
        if a_nrm > 1e-10:
            basis_cols.append(a_s / a_nrm)

    return np.column_stack(basis_cols) if basis_cols else np.zeros((n_genes, 0), dtype=np.float64)


def fit_source_cipher_lambdas(
    source_names: Sequence[str],
    source_covs: Dict[str, np.ndarray],
    source_deltas_by_line: Dict[str, Dict[str, np.ndarray]],
    gene_to_col: Dict[str, int],
    q_shared: np.ndarray,
) -> Tuple[float, float]:
    """Fit global CIPHER slope and LOO cross-source residual scaling lambda_ls on source lines only."""
    n_g = next(iter(source_covs.values())).shape[0]
    g_num, g_den = 0.0, 0.0
    ls_num, ls_den = 0.0, 0.0

    for s_hold in source_names:
        others = [s for s in source_names if s != s_hold]
        for tgt, d_true in source_deltas_by_line[s_hold].items():
            g_col = gene_to_col.get(tgt, -1)
            if g_col < 0:
                continue
            v_hold = _unit_col_vec(source_covs[s_hold], g_col)
            mask = np.ones(n_g, dtype=bool)
            mask[g_col] = False
            g_num += float(np.dot(v_hold[mask], d_true[mask]))
            g_den += float(np.dot(v_hold[mask], v_hold[mask]))

            if not others:
                continue
            avail_other = [o for o in others if tgt in source_deltas_by_line[o]]
            if not avail_other:
                continue
            res_true = d_true - np.mean([source_deltas_by_line[o][tgt] for o in avail_other], axis=0)
            res_v = v_hold - np.mean([_unit_col_vec(source_covs[o], g_col) for o in avail_other], axis=0)
            if q_shared.shape[1] > 0:
                res_true = res_true - q_shared @ (q_shared.T @ res_true)
                res_v = res_v - q_shared @ (q_shared.T @ res_v)
            ls_num += float(np.dot(res_v[mask], res_true[mask]))
            ls_den += float(np.dot(res_v[mask], res_v[mask]))

    lam_global = float(g_num / g_den) if g_den > EPS else 0.0
    lam_ls = float(ls_num / ls_den) if ls_den > EPS else 0.0
    return lam_global, lam_ls


def build_predictions(
    cohort: AlignedCohort,
    shrink: float = 0.0,
    cov_max_cells: Optional[int] = None,
    seed: int = 42,
    include_sa3_anchors: bool = True,
    include_exploratory_e: bool = True,
    precomputed_covs: Optional[Dict[str, np.ndarray]] = None,
    precomputed_gene_orders: Optional[Dict[str, Sequence[str]]] = None,
) -> Tuple[Dict[str, np.ndarray], np.ndarray, np.ndarray, Dict[str, float]]:
    """
    Build all Round 11 zero-shot prediction matrices + reusable `fc_bar`, `grand_mean`, and baseline quantiles.
    """
    n_t = len(cohort.eval_targets)
    n_g = len(cohort.eval_genes)

    fc_bar, t_mean, grand_mean = compute_source_fc_bar(
        eval_targets=cohort.eval_targets,
        source_names=cohort.source_names,
        source_deltas_by_line=cohort.source_deltas_by_line,
        source_ctrl_mean_libnorm=cohort.source_ctrl_mean_libnorm,
    )

    mu_t = np.maximum(np.asarray(cohort.tgt_ctrl_mean_libnorm, dtype=np.float64), 0.0)
    c_med = float(np.median(mu_t))
    c_p10 = float(np.percentile(mu_t, 10.0))
    c_p90 = float(np.percentile(mu_t, 90.0))

    preds: Dict[str, np.ndarray] = {
        "transferred_mean": t_mean,
        "source_grand_mean": np.broadcast_to(grand_mean[None, :], (n_t, n_g)).copy(),
        "fct": fct_from_fc_bar(fc_bar, mu_t),
        "e1_const_med": fct_from_fc_bar(fc_bar, np.full(n_g, c_med, dtype=np.float64)),
    }
    const_stats: Dict[str, float] = {
        "mu_t_p10": c_p10,
        "mu_t_median": c_med,
        "mu_t_p90": c_p90,
    }

    if cohort.tgt_ctrl_mean_libnorm_full is not None:
        preds["fct_full_nt"] = fct_from_fc_bar(fc_bar, cohort.tgt_ctrl_mean_libnorm_full)

    for s_name in cohort.source_names:
        preds[f"fct_mismatch::{s_name}"] = fct_from_fc_bar(fc_bar, cohort.source_ctrl_mean_libnorm[s_name])

    if include_exploratory_e:
        mu_src_mean = np.mean(
            [np.maximum(np.asarray(cohort.source_ctrl_mean_libnorm[s], dtype=np.float64), 0.0) for s in cohort.source_names],
            axis=0,
        )
        mu_t_sorted = np.sort(mu_t, kind="mergesort")
        rank_src = np.argsort(np.argsort(mu_src_mean, kind="mergesort"), kind="mergesort")
        mu_rank_matched = mu_t_sorted[rank_src]

        preds["e2_const_p10"] = fct_from_fc_bar(fc_bar, np.full(n_g, c_p10, dtype=np.float64))
        preds["e2_const_p90"] = fct_from_fc_bar(fc_bar, np.full(n_g, c_p90, dtype=np.float64))
        preds["e3_src_mean_ctrl"] = fct_from_fc_bar(fc_bar, mu_src_mean)
        preds["e5_rank_matched_mu"] = fct_from_fc_bar(fc_bar, mu_rank_matched)
        preds["e6_log_fc_bar"] = np.log(np.maximum(fc_bar, 1e-300))
        preds["e7_fc_minus_1"] = fc_bar - 1.0
        const_stats["mu_e3_median"] = float(np.median(mu_src_mean))

    if include_sa3_anchors:
        def _get_cov(line_name: str, cells: np.ndarray) -> np.ndarray:
            if precomputed_covs is not None and line_name in precomputed_covs:
                base_cov = precomputed_covs[line_name]
                if precomputed_gene_orders is not None and line_name in precomputed_gene_orders:
                    src_genes = precomputed_gene_orders[line_name]
                    if list(src_genes) != list(cohort.eval_genes):
                        idx = _align_line_cols(src_genes, cohort.eval_genes)
                        return np.asarray(base_cov[np.ix_(idx, idx)], dtype=np.float64)
                if base_cov.shape == (n_g, n_g):
                    return base_cov
            return sample_covariance(deterministic_subsample(cells, cov_max_cells, seed=seed), shrink=shrink)

        tgt_cov = _get_cov(cohort.target_line, cohort.tgt_ctrl_cells)
        src_covs = {s: _get_cov(s, cohort.source_ctrl_cells[s]) for s in cohort.source_names}
        tgt_diag = np.maximum(np.diag(tgt_cov), 1e-8)
        src_diag_mean = np.maximum(
            np.mean([np.maximum(np.diag(src_covs[s]), 1e-8) for s in cohort.source_names], axis=0),
            1e-8,
        )
        c_global = float(np.sqrt(np.sum(tgt_diag) / max(float(np.sum(src_diag_mean)), EPS)))
        preds["transfer_mean_global_scale"] = c_global * t_mean

        var_scale = np.sqrt(tgt_diag / src_diag_mean)
        tmvs = t_mean * var_scale[None, :]
        preds["transfer_mean_var_scaled"] = tmvs

        q_3d = fit_3d_shared_subspace(cohort.source_deltas_by_line, n_genes=n_g)
        g_to_col = {g: idx for idx, g in enumerate(cohort.eval_genes)}
        lam_global, _ = fit_source_cipher_lambdas(
            source_names=cohort.source_names,
            source_covs=src_covs,
            source_deltas_by_line=cohort.source_deltas_by_line,
            gene_to_col=g_to_col,
            q_shared=q_3d,
        )

        cipher = np.zeros((n_t, n_g), dtype=np.float64)
        hyb = np.zeros((n_t, n_g), dtype=np.float64)
        hyb_tmvs = np.zeros((n_t, n_g), dtype=np.float64)
        for i, tgt in enumerate(cohort.eval_targets):
            g_col = int(cohort.target_cols[i])
            avail = [s for s in cohort.source_names if tgt in cohort.source_deltas_by_line[s]]
            v_t = _unit_col_vec(tgt_cov, g_col)
            v_s_mean = np.mean([_unit_col_vec(src_covs[s], g_col) for s in avail], axis=0)
            diff = v_t - v_s_mean
            if q_3d.shape[1] > 0:
                diff = diff - q_3d @ (q_3d.T @ diff)
            ls_term = lam_global * diff
            cipher[i] = lam_global * v_t
            hyb[i] = t_mean[i] + ls_term
            hyb_tmvs[i] = tmvs[i] + ls_term

        preds["cipher_linear"] = cipher
        preds["hybrid"] = hyb
        preds["hybrid_on_tmvs"] = hyb_tmvs

    return preds, fc_bar, grand_mean, const_stats
