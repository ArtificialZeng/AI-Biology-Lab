#!/usr/bin/env python3
"""
Zero-shot cross-cell-line perturbation response predictors.

Strict zero-leakage contract:
  For any held-out target cell line `t`, predictors may ONLY access:
    1) Target control cells `X_ctrl_t` (or its covariance `C_t`), and
    2) The query target gene list `target_genes`.
  No perturbation responses from target line `t` may ever be passed to `fit_*` or `predict_*`.

Predictors implemented:
  1. Transfer_Mean                 : Average observed response across source lines for target g
                                     (falls back to Source_Grand_Mean if g unseen in sources).
  2. Source_Grand_Mean             : Constant mean across all source perturbations (DRF neg ctrl).
  3. CIPHER_Linear                 : lambda * C_t[:, g] / C_t[g, g], with single scalar lambda
                                     fitted strictly on source lines (Kuznets-Speck et al. 2026).
  4. Hybrid_Raw                    : Transfer_Mean + (pred_CIPHER(t, g) - mean_{s} pred_CIPHER(s, g)).
  5. Hybrid_ProjShared             : Transfer_Mean + P_perp(pred_CIPHER(t, g) - mean_{s} pred_CIPHER(s, g)),
                                     projecting out source-derived shared subspace span{1, b_S, a_S}.
  6. NegCtrl_MismatchedCov         : Replaces C_t with a mismatched cell line covariance C_{t'}.
  7. NegCtrl_ShuffledGenes         : Replaces C_t with P_pi C_t P_pi^T (gene labels permuted).
  8. NegCtrl_DiagonalCov           : Keeps only diag(C_t) (mean-field null destroying co-fluctuations).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np


@dataclass(frozen=True)
class SourceLineInput:
    """Source-line data permitted for training/fitting (never contains target line)."""

    name: str
    ctrl_cov: np.ndarray            # (G, G) float64 control covariance C_s
    deltas: np.ndarray              # (N_s, G) float64 perturbation responses Delta_s
    target_indices: np.ndarray      # (N_s,) int64 index of perturbed target gene in 0..G-1
    target_names: Tuple[str, ...]   # length N_s canonical target identifiers
    ctrl_mean_libnorm: Optional[np.ndarray] = None  # (G,) float64 unperturbed libnorm mean (for FCT)


def compute_control_covariance(
    ctrl_cells: np.ndarray,
    shrink: float = 0.0,
    min_var: float = 0.0,
) -> np.ndarray:
    """
    Compute (G, G) float64 gene-gene sample covariance from unperturbed control cells (N_c, G).

    Matches CIPHER's `compute_covariance(..., shrink=shrink)` using unbiased (N_c - 1) denominator
    and optional convex shrinkage `(1 - shrink) * Sigma + shrink * mean_diag * I`.
    """
    X = np.asarray(ctrl_cells, dtype=np.float64)
    if X.ndim != 2 or X.shape[0] < 2:
        raise ValueError(f"ctrl_cells must be 2D with >= 2 cells, got shape {X.shape}")
    Xc = X - X.mean(axis=0, keepdims=True)
    Sigma = (Xc.T @ Xc) / float(X.shape[0] - 1)
    Sigma = 0.5 * (Sigma + Sigma.T)  # enforce exact symmetry
    if shrink > 0.0:
        diag = np.diag(Sigma)
        good = np.isfinite(diag) & (diag > 0)
        mean_d = float(np.mean(diag[good])) if np.any(good) else 1.0
        Sigma = (1.0 - float(shrink)) * Sigma + float(shrink) * mean_d * np.eye(Sigma.shape[0])
    if min_var > 0.0:
        np.fill_diagonal(Sigma, np.maximum(np.diag(Sigma), float(min_var)))
    return Sigma


def shuffle_covariance_genes(cov: np.ndarray, seed: int = 42) -> np.ndarray:
    """Permute gene labels symmetrically: P_pi @ cov @ P_pi^T."""
    rng = np.random.RandomState(seed)
    G = cov.shape[0]
    perm = rng.permutation(G)
    return cov[np.ix_(perm, perm)].copy()


def diagonal_covariance(cov: np.ndarray) -> np.ndarray:
    """Zero out all off-diagonal gene-gene covariances while keeping marginal variances."""
    return np.diag(np.diag(cov).copy())


def cov_to_corr(cov: np.ndarray, min_var: float = 1e-8) -> np.ndarray:
    """Convert covariance matrix C to symmetric clipped Pearson correlation R = D^{-1/2} C D^{-1/2}."""
    C = np.asarray(cov, dtype=np.float64)
    diag = np.maximum(np.diag(C), float(min_var))
    inv_sd = 1.0 / np.sqrt(diag)
    corr = C * (inv_sd[:, None] * inv_sd[None, :])
    corr = np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)
    corr = np.clip(0.5 * (corr + corr.T), -1.0, 1.0)
    np.fill_diagonal(corr, 1.0)
    return corr


def cipher_unit_response_columns(
    cov: np.ndarray,
    target_indices: Sequence[int],
    min_var: float = 1e-8,
) -> np.ndarray:
    """
    Compute unscaled CIPHER unit-target-shift linear response vectors:
        r(g) = C[:, g] / C[g, g]
    Returns shape (len(target_indices), G). Columns with C[g, g] < min_var or non-finite return zeros.
    """
    G = cov.shape[0]
    idx = np.asarray(target_indices, dtype=np.int64)
    out = np.zeros((len(idx), G), dtype=np.float64)
    valid = (idx >= 0) & (idx < G)
    if not np.any(valid):
        return out
    v_idx = idx[valid]
    diag = np.diag(cov)[v_idx]
    ok = np.isfinite(diag) & (diag >= min_var)
    cols = cov[:, v_idx[ok]].T  # (n_ok, G)
    scaled = cols / diag[ok, None]
    scaled = np.nan_to_num(scaled, nan=0.0, posinf=0.0, neginf=0.0)
    ok_positions = np.where(valid)[0][ok]
    out[ok_positions] = scaled
    return out


def fit_shared_subspace(source_deltas: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Fit source-only shared perturbation centroid `b_S` and PC1 residual axis `a_S`,
    and return orthonormal QR basis `Q_S` for span{1, b_S, a_S}.
    Matches `fit_baseline_and_axis` in `src/xatlas_panhuman_solver.py`.
    """
    Y = np.asarray(source_deltas, dtype=np.float64)
    G = Y.shape[1]
    base = Y.mean(axis=0)
    Ys = Y - base[None, :]
    M0 = np.column_stack([np.ones(G, dtype=np.float64), base])
    Q0, R0 = np.linalg.qr(M0, mode="reduced")
    Q0 = Q0[:, np.abs(np.diag(R0)) > 1e-10 * max(1.0, np.abs(R0).max())]
    Yp = Ys - (Ys @ Q0) @ Q0.T
    norms = np.linalg.norm(Yp, axis=1, keepdims=True)
    keep = norms.ravel() > 1e-12
    if np.sum(keep) >= 2:
        Yn = Yp[keep] / norms[keep]
        _, _, Vt = np.linalg.svd(Yn, full_matrices=False)
        axis = Vt[0]
        if axis[np.argmax(np.abs(axis))] < 0:
            axis = -axis
    else:
        axis = np.zeros(G, dtype=np.float64)

    M_full = np.column_stack([np.ones(G, dtype=np.float64), base, axis])
    Q_full, R_full = np.linalg.qr(M_full, mode="reduced")
    Q_full = Q_full[:, np.abs(np.diag(R_full)) > 1e-10 * max(1.0, np.abs(R_full).max())]
    return base, axis, Q_full


def project_out_basis(mat: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Project columns/rows of `mat` (..., G) orthogonal to orthonormal basis `Q` (G, k)."""
    if Q.shape[1] == 0:
        return mat.copy()
    return mat - (mat @ Q) @ Q.T


def fit_cipher_lambda(
    source_lines: Sequence[SourceLineInput],
    mask_self: bool = True,
    min_var: float = 1e-8,
) -> float:
    """
    Fit the single global amplitude scalar `lambda` strictly on source lines by least squares:
        min_lambda sum_{s, g} || Delta_s(g) - lambda * r_s(g) ||_2^2
    where r_s(g) = C_s[:, g] / C_s[g, g].
    When `mask_self=True`, excludes each perturbation's own target gene g from the dot products.
    """
    num = 0.0
    den = 0.0
    for src in source_lines:
        R = cipher_unit_response_columns(src.ctrl_cov, src.target_indices, min_var=min_var)
        D = np.asarray(src.deltas, dtype=np.float64)
        if mask_self:
            G = D.shape[1]
            row_mask = np.ones_like(D, dtype=np.float64)
            valid_t = (src.target_indices >= 0) & (src.target_indices < G)
            row_mask[np.where(valid_t)[0], src.target_indices[valid_t]] = 0.0
            R = R * row_mask
            D = D * row_mask
        num += float(np.sum(R * D))
        den += float(np.sum(R * R))
    if den <= 1e-12 or not np.isfinite(num / den):
        return 0.0
    return float(num / den)


@dataclass
class SourceModelState:
    """Pre-fitted source-only state used for zero-shot prediction on a held-out target line."""

    source_lines: Tuple[SourceLineInput, ...]
    lam: float
    grand_mean: np.ndarray                # (G,)
    base_shared: np.ndarray               # (G,)
    axis_shared: np.ndarray               # (G,)
    Q_shared: np.ndarray                  # (G, k)
    target_to_source_rows: Dict[str, List[Tuple[int, int]]]  # target -> [(src_idx, row_idx)]

    @classmethod
    def fit(
        cls,
        source_lines: Sequence[SourceLineInput],
        mask_self: bool = True,
        min_var: float = 1e-8,
    ) -> "SourceModelState":
        if len(source_lines) == 0:
            raise ValueError("At least one source line is required.")
        lam = fit_cipher_lambda(source_lines, mask_self=mask_self, min_var=min_var)
        all_deltas = np.vstack([np.asarray(s.deltas, dtype=np.float64) for s in source_lines])
        grand_mean = all_deltas.mean(axis=0)
        base_shared, axis_shared, Q_shared = fit_shared_subspace(all_deltas)

        target_map: Dict[str, List[Tuple[int, int]]] = {}
        for s_idx, src in enumerate(source_lines):
            for r_idx, t_name in enumerate(src.target_names):
                target_map.setdefault(t_name, []).append((s_idx, r_idx))

        return cls(
            source_lines=tuple(source_lines),
            lam=lam,
            grand_mean=grand_mean,
            base_shared=base_shared,
            axis_shared=axis_shared,
            Q_shared=Q_shared,
            target_to_source_rows=target_map,
        )


def predict_zero_shot_suite(
    state: SourceModelState,
    target_ctrl_cov: np.ndarray,
    query_target_names: Sequence[str],
    query_target_gene_indices: Sequence[int],
    mismatched_ctrl_cov: Optional[np.ndarray] = None,
    shuffle_seed: int = 42,
    min_var: float = 1e-8,
    target_ctrl_mean_libnorm: Optional[np.ndarray] = None,
) -> Dict[str, np.ndarray]:
    """
    Generate predictions (N_q, G) for all required methods & negative controls on a target line.

    Target line inputs are strictly limited to `target_ctrl_cov` and `query_target_*`.
    """
    N_q = len(query_target_names)
    G = state.grand_mean.shape[0]
    q_idx = np.asarray(query_target_gene_indices, dtype=np.int64)

    # 1. Source_Grand_Mean and Transfer_Mean
    pred_grand_mean = np.broadcast_to(state.grand_mean[None, :], (N_q, G)).copy()
    pred_transfer = np.empty((N_q, G), dtype=np.float64)

    # Precompute unit responses on every source line for all query target genes
    src_cipher_preds = [
        state.lam * cipher_unit_response_columns(src.ctrl_cov, q_idx, min_var=min_var)
        for src in state.source_lines
    ]
    mean_src_cipher_all = np.mean(src_cipher_preds, axis=0)  # (N_q, G)
    mean_src_cipher_matched = np.empty((N_q, G), dtype=np.float64)

    for i, t_name in enumerate(query_target_names):
        hits = state.target_to_source_rows.get(t_name, [])
        if hits:
            pred_transfer[i] = np.mean(
                [state.source_lines[s_i].deltas[r_i] for s_i, r_i in hits], axis=0
            )
            src_ids = sorted(set(s_i for s_i, _ in hits))
            mean_src_cipher_matched[i] = np.mean(
                [src_cipher_preds[s_i][i] for s_i in src_ids], axis=0
            )
        else:
            pred_transfer[i] = state.grand_mean
            mean_src_cipher_matched[i] = mean_src_cipher_all[i]

    # 2. Target CIPHER linear response
    pred_cipher = state.lam * cipher_unit_response_columns(
        target_ctrl_cov, q_idx, min_var=min_var
    )

    # 3. Hybrid (raw differential and shared-subspace-projected differential)
    diff_matched = pred_cipher - mean_src_cipher_matched
    diff_matched_proj = project_out_basis(diff_matched, state.Q_shared)
    pred_hybrid_raw = pred_transfer + diff_matched
    pred_hybrid_proj = pred_transfer + diff_matched_proj

    preds: Dict[str, np.ndarray] = {
        "Source_Grand_Mean": pred_grand_mean,
        "Transfer_Mean": pred_transfer,
        "CIPHER_Linear": pred_cipher,
        "Hybrid_Raw": pred_hybrid_raw,
        "Hybrid_ProjShared": pred_hybrid_proj,
    }

    # 4a. Negative control: Mismatched cell line covariance (ALL source lines + optional explicit)
    for s_i, src in enumerate(state.source_lines):
        diff_mis_s = src_cipher_preds[s_i] - mean_src_cipher_matched
        preds[f"NegCtrl_HybridProj_MismatchedCov({src.name})"] = (
            pred_transfer + project_out_basis(diff_mis_s, state.Q_shared)
        )
    if mismatched_ctrl_cov is not None:
        pred_cipher_mis = state.lam * cipher_unit_response_columns(
            mismatched_ctrl_cov, q_idx, min_var=min_var
        )
        diff_mis = pred_cipher_mis - mean_src_cipher_matched
        preds["NegCtrl_CIPHER_MismatchedCov"] = pred_cipher_mis
        preds["NegCtrl_Hybrid_MismatchedCov"] = pred_transfer + diff_mis
        preds["NegCtrl_HybridProj_MismatchedCov"] = pred_transfer + project_out_basis(
            diff_mis, state.Q_shared
        )

    # 4b. Negative control: Shuffled gene labels covariance
    cov_shuf = shuffle_covariance_genes(target_ctrl_cov, seed=shuffle_seed)
    pred_cipher_shuf = state.lam * cipher_unit_response_columns(
        cov_shuf, q_idx, min_var=min_var
    )
    diff_shuf = pred_cipher_shuf - mean_src_cipher_matched
    preds["NegCtrl_CIPHER_ShuffledGenes"] = pred_cipher_shuf
    preds["NegCtrl_Hybrid_ShuffledGenes"] = pred_transfer + diff_shuf
    preds["NegCtrl_HybridProj_ShuffledGenes"] = pred_transfer + project_out_basis(
        diff_shuf, state.Q_shared
    )

    # 4c. Negative control: Diagonal covariance (mean-field null)
    cov_diag = diagonal_covariance(target_ctrl_cov)
    pred_cipher_diag = state.lam * cipher_unit_response_columns(
        cov_diag, q_idx, min_var=min_var
    )
    diff_diag = pred_cipher_diag - mean_src_cipher_matched
    preds["NegCtrl_CIPHER_DiagCov"] = pred_cipher_diag
    preds["NegCtrl_Hybrid_DiagCov"] = pred_transfer + diff_diag
    preds["NegCtrl_HybridProj_DiagCov"] = pred_transfer + project_out_basis(
        diff_diag, state.Q_shared
    )

    # 4d. PREREG_AUDIT_R1 R4(a): Source mean correlation matrix + Target marginal variances
    #     C = D_t^{1/2} * R_bar_S * D_t^{1/2}
    mean_src_corr = np.mean([cov_to_corr(s.ctrl_cov, min_var=min_var) for s in state.source_lines], axis=0)
    diag_t = np.maximum(np.diag(np.asarray(target_ctrl_cov, dtype=np.float64)), min_var)
    std_t = np.sqrt(diag_t)
    cov_src_corr_tgt_var = (std_t[:, None] * mean_src_corr) * std_t[None, :]
    pred_cipher_sctv = state.lam * cipher_unit_response_columns(
        cov_src_corr_tgt_var, q_idx, min_var=min_var
    )
    diff_sctv = pred_cipher_sctv - mean_src_cipher_matched
    preds["NegCtrl_HybridProj_SrcCorrTgtVar"] = pred_transfer + project_out_basis(
        diff_sctv, state.Q_shared
    )

    # 4e. PREREG_AUDIT_R1 R4(b) & PREREG_AUDIT_R2 R12: Marginal-variance-scaled Transfer_Mean (TMVS)
    #     Delta_hat_{p,t,g} = Delta_src_{p,g} * sqrt(C_{t,gg} / mean_s C_{s,gg})
    src_vars_by_line = [
        np.maximum(np.diag(np.asarray(s.ctrl_cov, dtype=np.float64)), min_var)
        for s in state.source_lines
    ]
    mean_src_var = np.maximum(np.mean(src_vars_by_line, axis=0), min_var)
    var_scale = np.sqrt(diag_t / mean_src_var)
    pred_tmvs = pred_transfer * var_scale[None, :]
    preds["Transfer_Mean_VarScaled"] = pred_tmvs

    # 4f. PREREG_AUDIT_R2 R12 (H-B3): TMVS MismatchedVar (all source lines) & ShuffledVar controls
    for s_i, src in enumerate(state.source_lines):
        mis_var_scale = np.sqrt(src_vars_by_line[s_i] / mean_src_var)
        preds[f"NegCtrl_TMVS_MismatchedVar({src.name})"] = pred_transfer * mis_var_scale[None, :]
    diag_shuf = np.maximum(np.diag(np.asarray(cov_shuf, dtype=np.float64)), min_var)
    shuf_var_scale = np.sqrt(diag_shuf / mean_src_var)
    preds["NegCtrl_TMVS_ShuffledVar"] = pred_transfer * shuf_var_scale[None, :]

    # 4g. PREREG_AUDIT_R2 R14: Hybrid_on_TMVS = TMVS + lambda * P_perp(v_t - v_bar_S)
    preds["Hybrid_on_TMVS"] = pred_tmvs + diff_matched_proj

    # 4h. PREREG_AUDIT_R2 R15: Fold-Change Transfer (FCT)
    #     Delta_hat = log1p(mu_t * FC_bar_src) - log1p(mu_t), geometric mean FC across source lines
    if target_ctrl_mean_libnorm is not None and all(
        s.ctrl_mean_libnorm is not None for s in state.source_lines
    ):
        mu_t = np.maximum(np.asarray(target_ctrl_mean_libnorm, dtype=np.float64), 0.0)
        mu_s_list = [
            np.maximum(np.asarray(s.ctrl_mean_libnorm, dtype=np.float64), min_var)
            for s in state.source_lines
        ]
        pred_fct = np.empty((N_q, G), dtype=np.float64)
        for i, t_name in enumerate(query_target_names):
            hits = state.target_to_source_rows.get(t_name, [])
            if hits:
                log_fcs = []
                for s_i, r_i in hits:
                    mu_s = mu_s_list[s_i]
                    pert_s = np.maximum(
                        np.expm1(np.log1p(mu_s) + state.source_lines[s_i].deltas[r_i]),
                        min_var,
                    )
                    log_fcs.append(np.log(pert_s / mu_s))
                fc_bar = np.exp(np.mean(log_fcs, axis=0))
                pred_fct[i] = np.log1p(mu_t * fc_bar) - np.log1p(mu_t)
            else:
                pred_fct[i] = pred_transfer[i]
        preds["FCT"] = pred_fct

    return preds
