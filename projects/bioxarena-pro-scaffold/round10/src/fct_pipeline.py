#!/usr/bin/env python3
"""
Round 10 Group C (`FCT`) Predictor Suite, Mismatched Controls, 100-Perm Null, and SA-3 R9 Anchors.

Zero target perturbation leakage: Target inputs are strictly limited to
unperturbed target control mean (`ctrl_mean_libnorm`), unperturbed control
covariance (`target_ctrl_cov`, only for SA-2/SA-3 comparators), and target gene indices.

FCT (`PREREG_A v4 §4.6` / `PREREG_C §3.1`): zero tuning, verbatim unchanged from Round 9.
Round 9 Anchor `Hybrid_ProjShared` (`SA-3` / `P7`): uses frozen `1ca9f60f` 3D `span{1, b_S, a_S}` projection.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

import numpy as np

from data_loader_r10 import LineDatasetR10


@dataclass
class SourceLineInputR10:
    """Source-line summary fitted strictly on source lines (zero target leakage)."""

    name: str
    target_names: List[str]
    target_indices: np.ndarray          # (N_t,) index in aligned genes (-1 if absent)
    deltas: np.ndarray                  # (N_t, G) float64
    ctrl_cov: Optional[np.ndarray]      # (G, G) float64
    ctrl_mean_libnorm: np.ndarray       # (G,) float64


def compute_control_covariance_r10(ctrl_cells: np.ndarray, shrink: float = 0.0) -> np.ndarray:
    """Unbiased sample covariance on unperturbed control cells with optional shrinkage toward diag."""
    X = np.asarray(ctrl_cells, dtype=np.float64)
    if X.ndim != 2 or X.shape[0] < 2:
        G = X.shape[1] if X.ndim == 2 else 0
        return np.eye(G, dtype=np.float64)
    Xc = X - X.mean(axis=0, keepdims=True)
    cov = (Xc.T @ Xc) / float(X.shape[0] - 1)
    cov = 0.5 * (cov + cov.T)
    if shrink > 0.0:
        cov = (1.0 - shrink) * cov + shrink * np.diag(np.diag(cov))
    return cov


def subsample_ctrl_r10(ctrl_cells: np.ndarray, cov_max_cells: Optional[int], seed: int = 42) -> np.ndarray:
    """Deterministic control subsampling matching Round 9 (`seed=42`)."""
    if cov_max_cells is None or ctrl_cells.shape[0] <= cov_max_cells:
        return ctrl_cells
    rng = np.random.default_rng(seed)
    idx = np.sort(rng.choice(ctrl_cells.shape[0], size=cov_max_cells, replace=False))
    return ctrl_cells[idx]


def build_source_input_r10(
    ds: LineDatasetR10,
    only_reproducible: bool = True,
    shrink: float = 0.0,
    cov_max_cells: Optional[int] = None,
    compute_cov: bool = True,
    seed: int = 42,
) -> SourceLineInputR10:
    """Construct source-only input structure filtering reproducible targets when `only_reproducible=True`."""
    if ds.ctrl_mean_libnorm is None:
        raise ValueError(f"Source dataset {ds.name} missing `ctrl_mean_libnorm` required for FCT.")
    if only_reproducible and ds.repr_mask is not None:
        keep = np.asarray(ds.repr_mask, dtype=bool)
    else:
        keep = np.ones(len(ds.target_names), dtype=bool)
    k_idx = np.where(keep)[0]
    t_names = [ds.target_names[int(i)] for i in k_idx]
    t_indices = ds.target_indices[k_idx]
    deltas = np.asarray(ds.deltas[k_idx], dtype=np.float64)

    cov = None
    if compute_cov:
        sub_c = subsample_ctrl_r10(ds.ctrl_cells, cov_max_cells=cov_max_cells, seed=seed)
        cov = compute_control_covariance_r10(sub_c, shrink=shrink)

    return SourceLineInputR10(
        name=ds.name,
        target_names=t_names,
        target_indices=t_indices,
        deltas=deltas,
        ctrl_cov=cov,
        ctrl_mean_libnorm=np.asarray(ds.ctrl_mean_libnorm, dtype=np.float64),
    )


def _cipher_unit_response_columns(
    cov: np.ndarray,
    target_indices: np.ndarray,
    min_var: float = 1e-8,
) -> np.ndarray:
    """Compute r(g) = C[:, g] / max(C[g, g], min_var) returning shape (len(target_indices), G)."""
    G = cov.shape[0]
    idx = np.asarray(target_indices, dtype=np.int64)
    out = np.zeros((len(idx), G), dtype=np.float64)
    valid = (idx >= 0) & (idx < G)
    if not np.any(valid):
        return out
    v_idx = idx[valid]
    diag = np.diag(cov)[v_idx]
    ok = np.isfinite(diag) & (diag >= min_var)
    cols = cov[:, v_idx[ok]].T
    scaled = np.nan_to_num(cols / diag[ok, None], nan=0.0, posinf=0.0, neginf=0.0)
    out[np.where(valid)[0][ok]] = scaled
    return out


def _fit_cipher_lambda(
    source_lines: Sequence[SourceLineInputR10],
    mask_self: bool = True,
    min_var: float = 1e-8,
) -> float:
    num = 0.0
    den = 0.0
    for src in source_lines:
        if src.ctrl_cov is None:
            continue
        R = _cipher_unit_response_columns(src.ctrl_cov, src.target_indices, min_var=min_var)
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


def _fit_shared_subspace_span3(source_deltas: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Frozen Round 9 (`1ca9f60f`) 3D orthonormal subspace `span{1, b_S, a_S}` per `PREREG_C §5 SA-3 P7`.
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


def _project_out_basis(mat: np.ndarray, Q: np.ndarray) -> np.ndarray:
    if Q.shape[1] == 0:
        return mat.copy()
    return mat - (mat @ Q) @ Q.T


@dataclass
class FctSourceState:
    """Source-only pre-fitted state for Group C (`FCT`), mismatched nulls, and SA-3 anchors."""

    source_lines: Tuple[SourceLineInputR10, ...]
    grand_mean: np.ndarray
    base_shared: np.ndarray
    axis_shared: np.ndarray
    Q_shared: np.ndarray
    lam: float
    target_to_source_rows: Dict[str, List[Tuple[int, int]]]
    has_covs: bool

    @classmethod
    def fit(
        cls,
        source_lines: Sequence[SourceLineInputR10],
        mask_self: bool = True,
        min_var: float = 1e-8,
    ) -> "FctSourceState":
        if len(source_lines) == 0:
            raise ValueError("At least one source line is required.")
        all_deltas = np.vstack([np.asarray(s.deltas, dtype=np.float64) for s in source_lines])
        grand_mean = all_deltas.mean(axis=0)
        base_shared, axis_shared, Q_shared = _fit_shared_subspace_span3(all_deltas)
        has_covs = all(s.ctrl_cov is not None for s in source_lines)
        lam = _fit_cipher_lambda(source_lines, mask_self=mask_self, min_var=min_var) if has_covs else 0.0

        target_map: Dict[str, List[Tuple[int, int]]] = {}
        for s_idx, src in enumerate(source_lines):
            for r_idx, t_name in enumerate(src.target_names):
                target_map.setdefault(t_name, []).append((s_idx, r_idx))

        return cls(
            source_lines=tuple(source_lines),
            grand_mean=grand_mean,
            base_shared=base_shared,
            axis_shared=axis_shared,
            Q_shared=Q_shared,
            lam=lam,
            target_to_source_rows=target_map,
            has_covs=has_covs,
        )


def compute_source_transfer_and_fc_bar(
    state: FctSourceState,
    query_target_names: Sequence[str],
    min_var: float = 1e-8,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Precompute `Transfer_Mean` (N_q, G), geometric-mean source fold-change `fc_bar` (N_q, G),
    and `has_hits` boolean mask (N_q,) strictly from source lines.
    Matches `PREREG_A v4 §4.6` / `PREREG_C §3.1` verbatim.
    """
    N_q = len(query_target_names)
    G = state.grand_mean.shape[0]
    mu_s_list = [
        np.maximum(np.asarray(s.ctrl_mean_libnorm, dtype=np.float64), min_var)
        for s in state.source_lines
    ]
    pred_transfer = np.empty((N_q, G), dtype=np.float64)
    fc_bar = np.ones((N_q, G), dtype=np.float64)
    has_hits = np.zeros(N_q, dtype=bool)

    for i, t_name in enumerate(query_target_names):
        hits = state.target_to_source_rows.get(t_name, [])
        if hits:
            has_hits[i] = True
            pred_transfer[i] = np.mean(
                [state.source_lines[s_i].deltas[r_i] for s_i, r_i in hits], axis=0
            )
            log_fcs = []
            for s_i, r_i in hits:
                mu_s = mu_s_list[s_i]
                pert_s = np.maximum(
                    np.expm1(np.log1p(mu_s) + state.source_lines[s_i].deltas[r_i]),
                    min_var,
                )
                log_fcs.append(np.log(pert_s / mu_s))
            fc_bar[i] = np.exp(np.mean(log_fcs, axis=0))
        else:
            has_hits[i] = False
            pred_transfer[i] = state.grand_mean
            fc_bar[i] = 1.0

    return pred_transfer, fc_bar, has_hits


def fct_from_fc_bar(
    mu_vec: np.ndarray,
    fc_bar: np.ndarray,
    pred_transfer: np.ndarray,
    has_hits: np.ndarray,
) -> np.ndarray:
    """Vectorized Fold-Change Transfer prediction given any 1D control mean vector `mu_vec`."""
    mu_safe = np.maximum(np.asarray(mu_vec, dtype=np.float64), 0.0)
    pred = np.log1p(mu_safe[None, :] * fc_bar) - np.log1p(mu_safe[None, :])
    if not np.all(has_hits):
        pred[~has_hits] = pred_transfer[~has_hits]
    return pred


def iter_fct_perm_predictions(
    state: FctSourceState,
    target_ctrl_mean_libnorm: np.ndarray,
    query_target_names: Sequence[str],
    n_perm: int = 100,
    min_var: float = 1e-8,
    precomputed: Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]] = None,
) -> Iterator[Tuple[int, np.ndarray]]:
    """
    Yield `(k, pred_perm_k)` for `k = 0..n_perm-1` permuting target control mean `mu_t` across
    `G` genes using `np.random.default_rng(k).permutation(G)` while keeping `fc_bar` unpermuted (`C2`).
    """
    mu_t = np.maximum(np.asarray(target_ctrl_mean_libnorm, dtype=np.float64), 0.0)
    G = mu_t.shape[0]
    if precomputed is None:
        pred_transfer, fc_bar, has_hits = compute_source_transfer_and_fc_bar(
            state, query_target_names, min_var=min_var
        )
    else:
        pred_transfer, fc_bar, has_hits = precomputed

    for k in range(n_perm):
        perm = np.random.default_rng(k).permutation(G)
        yield k, fct_from_fc_bar(mu_t[perm], fc_bar, pred_transfer, has_hits)


def predict_group_c_suite(
    state: FctSourceState,
    target_ctrl_mean_libnorm: np.ndarray,
    query_target_names: Sequence[str],
    query_target_gene_indices: Optional[Sequence[int]] = None,
    target_ctrl_cov: Optional[np.ndarray] = None,
    target_ctrl_mean_libnorm_full: Optional[np.ndarray] = None,
    min_var: float = 1e-8,
    include_r9_anchors: bool = True,
) -> Tuple[Dict[str, np.ndarray], Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Generate all zero-shot predictors for Round 10 Group C + negative controls + SA-2/SA-3 anchors.

    Returns:
      preds: dict of predictor name -> (N_q, G) matrix
      fc_bundle: `(pred_transfer, fc_bar, has_hits)` for streaming 100-perm computation
    """
    N_q = len(query_target_names)
    G = state.grand_mean.shape[0]
    mu_t = np.maximum(np.asarray(target_ctrl_mean_libnorm, dtype=np.float64), 0.0)
    fc_bundle = compute_source_transfer_and_fc_bar(state, query_target_names, min_var=min_var)
    pred_transfer, fc_bar, has_hits = fc_bundle

    preds: Dict[str, np.ndarray] = {
        "Source_Grand_Mean": np.broadcast_to(state.grand_mean[None, :], (N_q, G)).copy(),
        "Transfer_Mean": pred_transfer,
        "FCT": fct_from_fc_bar(mu_t, fc_bar, pred_transfer, has_hits),
    }

    # All source mismatched control means (6 lines in R10 holdout; 5 lines in R9 dev LOLO)
    for src in state.source_lines:
        mu_mis = np.maximum(np.asarray(src.ctrl_mean_libnorm, dtype=np.float64), 0.0)
        preds[f"NegCtrl_FCT_MismatchedMean({src.name})"] = fct_from_fc_bar(
            mu_mis, fc_bar, pred_transfer, has_hits
        )

    # Perm k=0 diagnostic row (full k=0..99 distribution computed via iter_fct_perm_predictions)
    perm0 = np.random.default_rng(0).permutation(G)
    preds["NegCtrl_FCT_Perm0"] = fct_from_fc_bar(mu_t[perm0], fc_bar, pred_transfer, has_hits)

    # SA-5b full-NT mean sensitivity when distinct full-cohort NT mean exists
    if target_ctrl_mean_libnorm_full is not None:
        mu_full = np.maximum(np.asarray(target_ctrl_mean_libnorm_full, dtype=np.float64), 0.0)
        preds["FCT_FullNTMean"] = fct_from_fc_bar(mu_full, fc_bar, pred_transfer, has_hits)

    # SA-2 GlobalScale & SA-3 Round 9 anchors (`P7`: Hybrid_ProjShared uses 3D span{1, b_S, a_S})
    if (
        include_r9_anchors
        and state.has_covs
        and target_ctrl_cov is not None
        and query_target_gene_indices is not None
    ):
        q_idx = np.asarray(query_target_gene_indices, dtype=np.int64)
        diag_t = np.maximum(np.diag(np.asarray(target_ctrl_cov, dtype=np.float64)), min_var)
        src_vars = [
            np.maximum(np.diag(np.asarray(s.ctrl_cov, dtype=np.float64)), min_var)
            for s in state.source_lines
            if s.ctrl_cov is not None
        ]
        mean_src_var = np.maximum(np.mean(src_vars, axis=0), min_var)

        # SA-2 R18 scalar global scaling control
        c_global = float(
            np.sqrt(max(float(np.mean(diag_t)), min_var) / max(float(np.mean(mean_src_var)), min_var))
        )
        preds["Transfer_Mean_GlobalScale"] = pred_transfer * c_global

        # Group B TMVS anchor
        var_scale = np.sqrt(diag_t / mean_src_var)
        pred_tmvs = pred_transfer * var_scale[None, :]
        preds["Transfer_Mean_VarScaled"] = pred_tmvs

        # CIPHER_Linear & Group A Hybrid_ProjShared (frozen 1ca9f60f 3D span{1, b_S, a_S})
        src_cipher_preds = [
            state.lam * _cipher_unit_response_columns(s.ctrl_cov, q_idx, min_var=min_var)
            for s in state.source_lines
            if s.ctrl_cov is not None
        ]
        mean_src_cipher_all = np.mean(src_cipher_preds, axis=0)
        mean_src_cipher_matched = np.empty((N_q, G), dtype=np.float64)
        for i, t_name in enumerate(query_target_names):
            hits = state.target_to_source_rows.get(t_name, [])
            if hits:
                src_ids = sorted(set(s_i for s_i, _ in hits))
                mean_src_cipher_matched[i] = np.mean(
                    [src_cipher_preds[s_i][i] for s_i in src_ids], axis=0
                )
            else:
                mean_src_cipher_matched[i] = mean_src_cipher_all[i]

        pred_cipher = state.lam * _cipher_unit_response_columns(
            target_ctrl_cov, q_idx, min_var=min_var
        )
        diff_matched_proj = _project_out_basis(pred_cipher - mean_src_cipher_matched, state.Q_shared)
        preds["CIPHER_Linear"] = pred_cipher
        preds["Hybrid_ProjShared"] = pred_transfer + diff_matched_proj
        preds["Hybrid_on_TMVS"] = pred_tmvs + diff_matched_proj

    return preds, fc_bundle
