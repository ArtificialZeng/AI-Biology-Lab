#!/usr/bin/env python3
"""Independent zero-shot cross-line predictor implementations for Team A audit (A-T3).

Written strictly from first principles without importing any code from
p2/vcc_linespecific/src/. Every predictor takes ONLY:
  - Source cell lines' perturbation deltas and control covariances
  - Target cell line's unperturbed control cells / covariance and target gene list
Zero target perturbation responses are ever read or accepted during fitting or prediction.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


EPS = 1e-12


def compute_control_covariance(
    ctrl_cells: np.ndarray,
    shrink: float = 1e-3,
    max_cells: Optional[int] = None,
    seed: int = 42,
) -> np.ndarray:
    """Estimate (G x G) gene-gene covariance from unperturbed target/source control cells."""
    x = np.asarray(ctrl_cells, dtype=np.float64)
    if x.ndim != 2 or x.shape[0] < 2:
        raise ValueError("ctrl_cells must be a 2D array with at least 2 cells.")
    if max_cells is not None and x.shape[0] > max_cells:
        rng = np.random.default_rng(seed)
        idx = np.sort(rng.choice(x.shape[0], size=max_cells, replace=False))
        x = x[idx]
    cov = np.cov(x, rowvar=False)
    cov = np.atleast_2d(cov).astype(np.float64)
    if shrink > 0.0:
        diag = np.diag(cov)
        good = np.isfinite(diag) & (diag > 0.0)
        scale = float(np.mean(diag[good])) if np.any(good) else 1.0
        cov = (1.0 - shrink) * cov + shrink * scale * np.eye(cov.shape[0])
    return cov


def shuffle_gene_labels_covariance(cov: np.ndarray, seed: int = 42) -> np.ndarray:
    """Negative control: symmetrically permute gene labels of covariance matrix C."""
    c = np.asarray(cov, dtype=np.float64)
    n_genes = c.shape[0]
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n_genes)
    # Ensure non-trivial permutation for small synthetic checks
    if n_genes > 1 and np.array_equal(perm, np.arange(n_genes)):
        perm = np.roll(perm, 1)
    return c[np.ix_(perm, perm)]


def meanfield_null_covariance(ctrl_cells: np.ndarray, shrink: float = 1e-3, seed: int = 42) -> np.ndarray:
    """Negative control: permute each gene independently across control cells."""
    x = np.asarray(ctrl_cells, dtype=np.float64).copy()
    rng = np.random.default_rng(seed)
    for g in range(x.shape[1]):
        rng.shuffle(x[:, g])
    return compute_control_covariance(x, shrink=shrink, max_cells=None, seed=seed)


def cipher_unit_vector(cov: np.ndarray, gene_col: int, min_var: float = 1e-8) -> np.ndarray:
    """Return CIPHER linear-response unit column C[:, g] / C[g, g]."""
    var_g = float(cov[gene_col, gene_col])
    if not np.isfinite(var_g) or var_g <= min_var:
        return np.zeros(cov.shape[0], dtype=np.float64)
    return np.asarray(cov[:, gene_col], dtype=np.float64) / var_g


def fit_source_cipher_lambda(
    source_covs: Dict[str, np.ndarray],
    source_deltas: Dict[str, np.ndarray],
    source_target_cols: Dict[str, Sequence[int]],
    exclude_self: bool = True,
    min_var: float = 1e-8,
) -> float:
    """Fit global scalar lambda strictly across source lines via OLS through origin."""
    num = 0.0
    den = 0.0
    for line, deltas in source_deltas.items():
        cov = source_covs[line]
        cols = source_target_cols[line]
        n_genes = cov.shape[0]
        for row_i, g_col in enumerate(cols):
            if g_col < 0 or g_col >= n_genes:
                continue
            v = cipher_unit_vector(cov, g_col, min_var=min_var)
            d = np.asarray(deltas[row_i], dtype=np.float64)
            if exclude_self:
                mask = np.ones(n_genes, dtype=bool)
                mask[g_col] = False
                v_m = v[mask]
                d_m = d[mask]
            else:
                v_m = v
                d_m = d
            num += float(np.dot(v_m, d_m))
            den += float(np.dot(v_m, v_m))
    return float(num / den) if den > EPS else 0.0


def fit_source_shared_axis_basis(source_deltas: Dict[str, np.ndarray]) -> np.ndarray:
    """Fit orthonormal basis Q for span{1, centroid, PC1} strictly on source deltas."""
    rows = [np.asarray(d, dtype=np.float64) for d in source_deltas.values() if len(d) > 0]
    if not rows:
        raise ValueError("Cannot fit shared basis with empty source deltas.")
    mat = np.vstack(rows)
    n_genes = mat.shape[1]
    centroid = np.mean(mat, axis=0)
    ones = np.ones(n_genes, dtype=np.float64)
    q0, r0 = np.linalg.qr(np.column_stack([ones, centroid]), mode="reduced")
    q0 = q0[:, np.abs(np.diag(r0)) > 1e-10 * max(1.0, float(np.max(np.abs(r0))))]

    centered = mat - centroid[None, :]
    proj_c = centered - (centered @ q0) @ q0.T
    norms = np.linalg.norm(proj_c, axis=1, keepdims=True)
    unit_rows = proj_c / np.maximum(norms, EPS)
    if unit_rows.shape[0] >= 2:
        _, _, vt = np.linalg.svd(unit_rows, full_matrices=False)
        pc1 = vt[0]
        m_full = np.column_stack([ones, centroid, pc1])
    else:
        m_full = np.column_stack([ones, centroid])
    q_full, r_full = np.linalg.qr(m_full, mode="reduced")
    keep = np.abs(np.diag(r_full)) > 1e-10 * max(1.0, float(np.max(np.abs(r_full))))
    return q_full[:, keep]


def fit_source_ls_lambda(
    source_covs: Dict[str, np.ndarray],
    source_deltas_by_target: Dict[str, Dict[str, np.ndarray]],
    target_to_col: Dict[str, int],
    exclude_self: bool = True,
    project_q: Optional[np.ndarray] = None,
    min_var: float = 1e-8,
) -> float:
    """Fit line-specific residual scaling lambda_ls strictly among source lines (LOO within sources)."""
    lines = sorted(source_covs.keys())
    if len(lines) < 2:
        return 1.0
    num = 0.0
    den = 0.0
    n_genes = next(iter(source_covs.values())).shape[0]
    for s_hold in lines:
        other_lines = [l for l in lines if l != s_hold]
        for tgt, d_true in source_deltas_by_target[s_hold].items():
            g_col = target_to_col.get(tgt, -1)
            if g_col < 0:
                continue
            other_ds = [source_deltas_by_target[ol][tgt] for ol in other_lines if tgt in source_deltas_by_target[ol]]
            if not other_ds:
                continue
            other_vs = [cipher_unit_vector(source_covs[ol], g_col, min_var=min_var) for ol in other_lines if tgt in source_deltas_by_target[ol]]
            res_true = np.asarray(d_true, dtype=np.float64) - np.mean(other_ds, axis=0)
            res_pred_unit = cipher_unit_vector(source_covs[s_hold], g_col, min_var=min_var) - np.mean(other_vs, axis=0)
            if project_q is not None:
                res_true = res_true - project_q @ (project_q.T @ res_true)
                res_pred_unit = res_pred_unit - project_q @ (project_q.T @ res_pred_unit)
            if exclude_self:
                mask = np.ones(n_genes, dtype=bool)
                mask[g_col] = False
                rt, rp = res_true[mask], res_pred_unit[mask]
            else:
                rt, rp = res_true, res_pred_unit
            num += float(np.dot(rp, rt))
            den += float(np.dot(rp, rp))
    return float(num / den) if den > EPS else 0.0


def predict_zero_shot_suite(
    eval_targets: Sequence[str],
    target_to_col: Dict[str, int],
    target_ctrl_cov: np.ndarray,
    source_covs: Dict[str, np.ndarray],
    source_deltas_by_target: Dict[str, Dict[str, np.ndarray]],
    mismatch_cov: np.ndarray,
    project_shared_axis: bool = False,
    fit_ls_lambda_on_sources: bool = False,
    seed: int = 42,
) -> Dict[str, np.ndarray]:
    """Produce all zero-shot predictors strictly without target perturbation data.

    Returns dict with keys:
      - transferred_mean
      - cipher_linear
      - hybrid
      - mismatch_cov_hybrid
      - shuffled_gene_hybrid
    """
    n_genes = target_ctrl_cov.shape[0]
    source_lines = sorted(source_covs.keys())

    # Assemble source deltas per line for lambda fitting
    s_deltas_arr: Dict[str, np.ndarray] = {}
    s_cols: Dict[str, List[int]] = {}
    for line in source_lines:
        d_list, c_list = [], []
        for tgt, vec in source_deltas_by_target[line].items():
            if tgt in target_to_col:
                d_list.append(vec)
                c_list.append(target_to_col[tgt])
        s_deltas_arr[line] = np.asarray(d_list, dtype=np.float64) if d_list else np.zeros((0, n_genes))
        s_cols[line] = c_list

    lam_global = fit_source_cipher_lambda(source_covs, s_deltas_arr, s_cols, exclude_self=True)
    q_shared = fit_source_shared_axis_basis(s_deltas_arr) if project_shared_axis else None
    if fit_ls_lambda_on_sources:
        lam_ls = fit_source_ls_lambda(
            source_covs, source_deltas_by_target, target_to_col, exclude_self=True, project_q=q_shared
        )
    else:
        lam_ls = lam_global

    shuf_cov = shuffle_gene_labels_covariance(target_ctrl_cov, seed=seed)
    shuf_src_covs = {
        l: shuffle_gene_labels_covariance(source_covs[l], seed=seed + 100 + idx)
        for idx, l in enumerate(source_lines)
    }

    n_eval = len(eval_targets)
    out = {
        "transferred_mean": np.zeros((n_eval, n_genes), dtype=np.float64),
        "cipher_linear": np.zeros((n_eval, n_genes), dtype=np.float64),
        "hybrid": np.zeros((n_eval, n_genes), dtype=np.float64),
        "mismatch_cov_hybrid": np.zeros((n_eval, n_genes), dtype=np.float64),
        "shuffled_gene_hybrid": np.zeros((n_eval, n_genes), dtype=np.float64),
        "shuffled_all_cov_hybrid": np.zeros((n_eval, n_genes), dtype=np.float64),
    }

    for i, tgt in enumerate(eval_targets):
        g_col = target_to_col[tgt]
        avail_lines = [l for l in source_lines if tgt in source_deltas_by_target[l]]
        if not avail_lines:
            raise ValueError(f"Target {tgt!r} missing from all source lines.")
        t_mean = np.mean([source_deltas_by_target[l][tgt] for l in avail_lines], axis=0).astype(np.float64)
        v_src_mean = np.mean([cipher_unit_vector(source_covs[l], g_col) for l in avail_lines], axis=0)
        v_src_shuf_mean = np.mean([cipher_unit_vector(shuf_src_covs[l], g_col) for l in avail_lines], axis=0)

        v_tgt = cipher_unit_vector(target_ctrl_cov, g_col)
        v_mis = cipher_unit_vector(mismatch_cov, g_col)
        v_shuf = cipher_unit_vector(shuf_cov, g_col)

        def _ls_term(v_cand: np.ndarray, v_ref: np.ndarray = v_src_mean) -> np.ndarray:
            diff = v_cand - v_ref
            if q_shared is not None:
                diff = diff - q_shared @ (q_shared.T @ diff)
            return lam_ls * diff

        out["transferred_mean"][i] = t_mean
        out["cipher_linear"][i] = lam_global * v_tgt
        out["hybrid"][i] = t_mean + _ls_term(v_tgt)
        out["mismatch_cov_hybrid"][i] = t_mean + _ls_term(v_mis)
        out["shuffled_gene_hybrid"][i] = t_mean + _ls_term(v_shuf)
        out["shuffled_all_cov_hybrid"][i] = t_mean + _ls_term(v_shuf, v_src_shuf_mean)

    return out
