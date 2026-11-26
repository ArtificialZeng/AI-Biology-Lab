#!/usr/bin/env python3
"""
Round 17 Cached-Embedding 56-Index Sweep & Phase-2 `G4` Transition Inference Engine (`PREREG_J` v2).

Reuses frozen Round 16 primitives (`encode_controls_se600m_r16`, `preprocess_hvg_counts_r16`,
`assert_zero_training_leakage`, `STATE_CHECKPOINT_FOLDS`, `MULTI_BATCH_ENS8_INDICES`,
`PRED_CLIP_MAX`) without modifying any file under `round16/`.

Key architectural properties:
  1. Caches `SE-600M` `X_state` (`2,058`-dim) and `ST-HVG` `X_hvg` (`2,000`-dim) control
     representations ONCE per pool and loops all 56 batch indices strictly inside
     `StateTransitionPerturbationModel` (loading each zero-shot fold checkpoint ONCE per phase).
  2. Preserves exact bit-level decoding parity with frozen `ArcStateAdapter._infer_single_fold`:
     `padded=True` (`16 x 32` sets), `clip(gene_pred, 0, 14)` (`gene_decoder` for `STATE_SE`,
     direct `preds` for `STATE_HVG`), equal-weight 4-fold mean (`hepg2, jurkat, k562, rpe1`),
     subtracting `'non-targeting'` under the identical batch index `b`.
  3. Saves 4-fold average and per-fold `'non-targeting'` pseudobulks on `G_eval` for Ruling 3
     `L2_ntc_recon_match`, plus secondary neutral-batch comparators (`bypass_zero`,
     `mean_batch_emb`, `multi_batch_ens8`, weight-only `medoid_idx`) and per-fold 56-index
     point `r_orth_cal` summaries without retaining `(4, 56, P, G)` in RAM.
  4. Strictly enforces `assert_all_execution_gates_r17` (`R17_CODE_FREEZE.txt` + `R17_GO.txt`)
     before any real-cell forward pass.
"""

from __future__ import annotations

import gc
import os
import pickle
import random
import sys
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

ROUND17_SRC_DIR = Path(__file__).resolve().parent
REPO_ROOT = ROUND17_SRC_DIR.parents[1]
R16_SRC_DIR = REPO_ROOT / "round16" / "src"
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
for _p in (str(ROUND17_SRC_DIR), str(R16_SRC_DIR), str(R10_SRC_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from metrics_r10 import (  # noqa: E402
    _extract_masked_idx_per_pert,
    _moments_to_corr_and_slope,
    disattenuated_r,
    spearman_brown_ceiling,
)
from metrics_r16 import _collect_mode_moments, _moments_to_safe_corr  # noqa: E402
from model_adapters_r16 import (  # noqa: E402
    DEFAULT_MODELS_DIR,
    MULTI_BATCH_ENS8_INDICES,
    PRED_CLIP_MAX,
    STATE_CHECKPOINT_FOLDS,
    assert_zero_training_leakage,
)
from ntc_extractors_r17 import (  # noqa: E402
    DEFAULT_R17_PREREG_LOCK_PATH,
    DEFAULT_R17_STATUS_DIR,
    assert_all_execution_gates_r17,
)
from stats_r17 import DEFAULT_R17_SEED, N_BATCH_INDICES_STATE  # noqa: E402


def _resolve_fold_dirs_r17(variant: str, models_dir: Path = DEFAULT_MODELS_DIR) -> List[Path]:
    if variant not in ("STATE_SE", "STATE_HVG"):
        raise ValueError(f"Unsupported STATE variant: {variant!r}")
    repo_name = "ST-SE-Replogle" if variant == "STATE_SE" else "ST-HVG-Replogle"
    dirs = [Path(models_dir) / repo_name / "zeroshot" / f for f in STATE_CHECKPOINT_FOLDS]
    for d in dirs:
        ckpt = d / "checkpoints" / "final.ckpt"
        if not ckpt.is_file():
            raise FileNotFoundError(f"STATE checkpoint not found: {ckpt}")
    return dirs


def compute_batch_encoder_medoid_from_weights_r17(
    fold_weights: Sequence[np.ndarray],
) -> Dict[str, Any]:
    """
    Compute the weight-only geometric medoid index `b_medoid in {0..55}` from shipped
    `batch_encoder.weight` matrices (`(56, 128)` per fold):
      - Primary (locked): Euclidean medoid on the 4-fold mean weight matrix `W_mean`,
        minimizing `sum_{b'=0}^{55} ||W_mean[b] - W_mean[b']||_2` (lowest-index tie-break).
      - Companion diagnostic: cosine medoid maximizing mean cosine similarity to all 56 rows
        (reported transparently alongside Euclidean medoid).
    """
    w_stack = np.stack([np.asarray(w, dtype=np.float64) for w in fold_weights], axis=0)
    if w_stack.ndim != 3 or w_stack.shape[1] != N_BATCH_INDICES_STATE:
        raise ValueError(f"Expected (n_folds, 56, D) weights, got {w_stack.shape}.")
    w_mean = np.mean(w_stack, axis=0)  # (56, D)

    # Pairwise Euclidean distance matrix
    diff = w_mean[:, None, :] - w_mean[None, :, :]
    dist_mat = np.linalg.norm(diff, axis=-1)  # (56, 56)
    sum_l2 = np.sum(dist_mat, axis=1)
    euc_medoid = int(np.argmin(sum_l2))

    # Companion cosine medoid
    norms = np.maximum(np.linalg.norm(w_mean, axis=1, keepdims=True), 1e-12)
    w_unit = w_mean / norms
    cos_mat = w_unit @ w_unit.T
    mean_cos = np.mean(cos_mat, axis=1)
    cos_medoid = int(np.argmax(mean_cos))

    return {
        "medoid_idx": euc_medoid,
        "medoid_rule": "euclidean_geometric_medoid_on_4fold_mean_batch_encoder_weight",
        "sum_l2_min": float(sum_l2[euc_medoid]),
        "sum_l2_all_56": [float(x) for x in sum_l2],
        "cosine_medoid_idx_diag": cos_medoid,
        "row_l2_norms_mean_weight": [float(x) for x in norms.ravel()],
    }


def compute_batch_encoder_medoid_idx_r17(
    variant: str,
    models_dir: Path = DEFAULT_MODELS_DIR,
) -> Dict[str, Any]:
    """Read-only extraction of `batch_encoder.weight` from the 4 zero-shot folds to get `medoid_idx`."""
    import torch

    fold_dirs = _resolve_fold_dirs_r17(variant, models_dir=models_dir)
    weights: List[np.ndarray] = []
    for fdir in fold_dirs:
        ckpt = torch.load(fdir / "checkpoints" / "final.ckpt", map_location="cpu", weights_only=False)
        sd = ckpt["state_dict"] if isinstance(ckpt, dict) and "state_dict" in ckpt else ckpt
        w_key = next((k for k in sd.keys() if k.endswith("batch_encoder.weight")), None)
        if w_key is None:
            raise KeyError(f"batch_encoder.weight missing from {fdir}/checkpoints/final.ckpt")
        weights.append(sd[w_key].detach().cpu().float().numpy())
        del ckpt, sd
    res = compute_batch_encoder_medoid_from_weights_r17(weights)
    res["variant"] = variant
    return res


def _eval_fold_point_r_orth_56(
    delta_fold_56: np.ndarray,
    delta_obs: np.ndarray,
    delta_half_a: np.ndarray,
    delta_half_b: np.ndarray,
    half_valid: np.ndarray,
    delta_transfer: np.ndarray,
    target_indices: Sequence[int],
) -> np.ndarray:
    """Compute exact 56-index point `r_orth_cal` for a single zero-shot fold (secondary diagnostic)."""
    d_obs = np.asarray(delta_obs, dtype=np.float64)
    tm = np.asarray(delta_transfer, dtype=np.float64)
    rows, idx_list = _extract_masked_idx_per_pert(
        d_obs, target_indices, gene_keep_mask=None, mode="all_nontarget", top_k=50
    )
    if not rows:
        return np.full(N_BATCH_INDICES_STATE, float("nan"), dtype=np.float64)
    r_obs = d_obs - tm
    ha_res = np.asarray(delta_half_a, dtype=np.float64) - tm
    hb_res = np.asarray(delta_half_b, dtype=np.float64) - tm
    pred_res = [delta_fold_56[b] - tm for b in range(N_BATCH_INDICES_STATE)]
    moms, hm = _collect_mode_moments(r_obs, pred_res, tm, rows, idx_list, None, ha_res, hb_res, half_valid)
    if hm is None:
        return np.full(N_BATCH_INDICES_STATE, float("nan"), dtype=np.float64)
    hnv, hsx, hsy, hsxx, hsyy, hsxy = hm
    _, pt_rho, _, _, _ = _moments_to_corr_and_slope(
        float(hnv.sum()), float(hsx.sum()), float(hsy.sum()),
        float(hsxx.sum()), float(hsyy.sum()), float(hsxy.sum()),
    )
    pt_ceil = float(spearman_brown_ceiling(float(pt_rho)))
    out = np.empty(N_BATCH_INDICES_STATE, dtype=np.float64)
    for b in range(N_BATCH_INDICES_STATE):
        nv, sx, sy, sxx, syy, sxy = moms[b]
        pt_r = float(
            _moments_to_safe_corr(
                float(nv.sum()), float(sx.sum()), float(sy.sum()),
                float(sxx.sum()), float(syy.sum()), float(sxy.sum()),
                zero_null_to_zero=False,
            )
        )
        out[b] = float(disattenuated_r(pt_r, pt_ceil))
    return out


def _map_eval_gene_cols_r17(
    g_model: Sequence[str],
    gene_names: Sequence[str],
    gene_canon_map: Optional[Mapping[str, str]] = None,
) -> np.ndarray:
    g_model_canon_to_idx: Dict[str, int] = {}
    for idx, gm in enumerate(g_model):
        cgm = gene_canon_map.get(gm, gm) if gene_canon_map is not None else gm
        g_model_canon_to_idx.setdefault(cgm, idx)
    col_sel: List[int] = []
    for g in gene_names:
        cg = gene_canon_map.get(g, g) if gene_canon_map is not None else g
        if cg not in g_model_canon_to_idx:
            raise KeyError(f"Requested eval gene {g!r} (canon={cg!r}) absent from G_model.")
        col_sel.append(g_model_canon_to_idx[cg])
    return np.asarray(col_sel, dtype=np.int64)


def _eval_batch_tensor_r17(
    model: Any,
    ctrl_dev: Any,
    pert_dev_list: Sequence[Any],
    b_dev: Any,
    col_sel: np.ndarray,
    variant: str,
    n_perts: int,
    n_genes: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """Run one batch tensor over `['non-targeting'] + pert_names` and return `(pb_ntc, delta_perts)`."""
    pb_ntc = np.empty(n_genes, dtype=np.float64)
    pb_perts = np.empty((n_perts, n_genes), dtype=np.float64)
    for qi, p_dev in enumerate(pert_dev_list):
        pred_raw = model.forward(
            {"ctrl_cell_emb": ctrl_dev, "pert_emb": p_dev, "batch": b_dev},
            padded=True,
        )
        gene_pred = model.gene_decoder(pred_raw) if variant == "STATE_SE" else pred_raw
        g_np = np.clip(gene_pred.detach().cpu().numpy().astype(np.float64), 0.0, PRED_CLIP_MAX)
        mean_vec = g_np.mean(axis=0)[col_sel]
        if qi == 0:
            pb_ntc[:] = mean_vec
        else:
            pb_perts[qi - 1] = mean_vec
    return pb_ntc, pb_perts - pb_ntc[None, :]


def run_state_sweep56_r17(
    variant: str,
    target_line: str,
    ctrl_tensor_np: np.ndarray,
    pert_names: Sequence[str],
    gene_names: Sequence[str],
    pert_canon_to_model: Mapping[str, str],
    gene_canon_map: Mapping[str, str],
    models_dir: Path = DEFAULT_MODELS_DIR,
    status_dir: Path = DEFAULT_R17_STATUS_DIR,
    prereg_lock_path: Optional[Path] = DEFAULT_R17_PREREG_LOCK_PATH,
    device: str = "mps",
    seed: int = DEFAULT_R17_SEED,
    obs_for_per_fold_diag: Optional[Mapping[str, Any]] = None,
    allow_synthetic_without_freeze: bool = False,
) -> Dict[str, Any]:
    """
    Phase 1 (`PREREG_J §1.2, §3.4`): Run all 56 trained batch indices (`b = 0..55`) plus
    `bypass_zero` and `mean_batch_emb` on cached `ctrl_tensor_np` (`X_state` or `X_hvg`),
    loading each of the 4 zero-shot folds ONCE.
    Returns exact float64 4-fold averaged arrays, per-fold NTC profiles, per-fold 56-index
    point `r_orth_cal` summaries, `L1_ens56`, `multi_batch_ens8`, and weight-only `medoid_idx`.
    """
    assert_zero_training_leakage(variant, target_line, checkpoint_id="zeroshot/ensemble4")
    if not allow_synthetic_without_freeze:
        assert_all_execution_gates_r17(status_dir=status_dir, prereg_lock_path=prereg_lock_path)

    os.environ.setdefault("WANDB_MODE", "disabled")
    os.environ.setdefault("HF_HUB_OFFLINE", "1")
    os.environ.setdefault("TRANSFORMERS_OFFLINE", "1")

    import torch
    from state.tx.models.state_transition import StateTransitionPerturbationModel

    ctrl_np = np.asarray(ctrl_tensor_np, dtype=np.float32)
    n_cells = int(ctrl_np.shape[0])
    if n_cells <= 0 or n_cells % 32 != 0:
        raise ValueError(f"ctrl_tensor_np rows ({n_cells}) must be a positive multiple of 32.")

    fold_dirs = _resolve_fold_dirs_r17(variant, models_dir=models_dir)
    n_folds = len(fold_dirs)
    fold_names = [f"zeroshot/{d.name}" for d in fold_dirs]
    n_perts = len(pert_names)
    n_genes = len(gene_names)

    pert_keys_model = [pert_canon_to_model.get(p, p) for p in pert_names]
    all_query = ["non-targeting"] + pert_keys_model

    dev_use = device if (device != "mps" or torch.backends.mps.is_available()) else "cpu"
    dev = torch.device(dev_use)

    delta_ens_56 = np.zeros((N_BATCH_INDICES_STATE, n_perts, n_genes), dtype=np.float64)
    ntc_pb_per_fold = np.zeros((n_folds, N_BATCH_INDICES_STATE, n_genes), dtype=np.float64)
    delta_bypass_zero = np.zeros((n_perts, n_genes), dtype=np.float64)
    delta_mean_batch_emb = np.zeros((n_perts, n_genes), dtype=np.float64)
    per_fold_r_orth_56: Optional[np.ndarray] = (
        np.full((n_folds, N_BATCH_INDICES_STATE), float("nan"), dtype=np.float64)
        if obs_for_per_fold_diag is not None
        else None
    )
    fold_batch_weights: List[np.ndarray] = []

    for fi, fdir in enumerate(fold_dirs):
        with open(fdir / "var_dims.pkl", "rb") as fh:
            g_model = list(pickle.load(fh)["gene_names"])
        col_sel = _map_eval_gene_cols_r17(g_model, gene_names, gene_canon_map)
        pert_map = torch.load(fdir / "pert_onehot_map.pt", map_location="cpu", weights_only=False)

        random.seed(int(seed))
        np.random.seed(int(seed))
        torch.manual_seed(int(seed))
        if dev_use == "mps" and hasattr(torch, "mps") and hasattr(torch.mps, "manual_seed"):
            torch.mps.manual_seed(int(seed))

        ckpt_file = fdir / "checkpoints" / "final.ckpt"
        model = StateTransitionPerturbationModel.load_from_checkpoint(str(ckpt_file), map_location="cpu")
        model.eval().to(dev)

        orig_w = model.batch_encoder.weight.detach().clone()
        fold_batch_weights.append(orig_w.detach().cpu().float().numpy())
        ctrl_dev = torch.from_numpy(ctrl_np).to(dev)

        # Pre-stage perturbation embeddings once per fold (exact same tensor values, avoids 58x host->device copies)
        pert_dev_list: List[torch.Tensor] = []
        for pk in all_query:
            if pk not in pert_map:
                raise KeyError(f"Perturbation key {pk!r} missing from {fdir}/pert_onehot_map.pt")
            pert_dev_list.append(pert_map[pk].to(device=dev, dtype=torch.float32).unsqueeze(0).expand(n_cells, -1))

        delta_fold_56 = np.empty((N_BATCH_INDICES_STATE, n_perts, n_genes), dtype=np.float64)

        with torch.no_grad():
            model.batch_encoder.weight.copy_(orig_w)
            for b in range(N_BATCH_INDICES_STATE):
                b_dev = torch.full((n_cells,), int(b), dtype=torch.long, device=dev)
                ntc_v, d_b = _eval_batch_tensor_r17(
                    model, ctrl_dev, pert_dev_list, b_dev, col_sel, variant, n_perts, n_genes
                )
                ntc_pb_per_fold[fi, b] = ntc_v
                delta_fold_56[b] = d_b
                delta_ens_56[b] += d_b / float(n_folds)

            # Secondary neutral-batch comparators (`PREREG_J §3.4` Item 4): bypass_zero & mean_batch_emb
            b_zero = torch.zeros(n_cells, dtype=torch.long, device=dev)
            model.batch_encoder.weight.copy_(orig_w)
            model.batch_encoder.weight[0].zero_()
            _, d_bz = _eval_batch_tensor_r17(
                model, ctrl_dev, pert_dev_list, b_zero, col_sel, variant, n_perts, n_genes
            )
            delta_bypass_zero += d_bz / float(n_folds)

            model.batch_encoder.weight.copy_(orig_w)
            model.batch_encoder.weight[0].copy_(orig_w.mean(dim=0))
            _, d_mb = _eval_batch_tensor_r17(
                model, ctrl_dev, pert_dev_list, b_zero, col_sel, variant, n_perts, n_genes
            )
            delta_mean_batch_emb += d_mb / float(n_folds)

            model.batch_encoder.weight.copy_(orig_w)
            if dev_use == "mps":
                torch.mps.synchronize()

        if obs_for_per_fold_diag is not None and per_fold_r_orth_56 is not None:
            per_fold_r_orth_56[fi] = _eval_fold_point_r_orth_56(
                delta_fold_56=delta_fold_56,
                delta_obs=obs_for_per_fold_diag["delta_obs"],
                delta_half_a=obs_for_per_fold_diag["delta_half_a"],
                delta_half_b=obs_for_per_fold_diag["delta_half_b"],
                half_valid=obs_for_per_fold_diag["half_valid"],
                delta_transfer=obs_for_per_fold_diag["delta_transfer"],
                target_indices=obs_for_per_fold_diag["target_indices"],
            )

        del delta_fold_56, pert_dev_list, ctrl_dev, orig_w, model
        gc.collect()
        if dev_use == "mps":
            torch.mps.empty_cache()

    medoid_info = compute_batch_encoder_medoid_from_weights_r17(fold_batch_weights)
    medoid_idx = int(medoid_info["medoid_idx"])

    ntc_pb_ensemble_56 = np.mean(ntc_pb_per_fold, axis=0)
    delta_l1_ens56 = np.mean(delta_ens_56, axis=0)
    delta_ens8 = np.mean(delta_ens_56[list(MULTI_BATCH_ENS8_INDICES)], axis=0)
    delta_medoid = delta_ens_56[medoid_idx].copy()

    return {
        "delta_pred_56": delta_ens_56,
        "ntc_pb_ensemble_56": ntc_pb_ensemble_56,
        "ntc_pb_per_fold_56": ntc_pb_per_fold,
        "fold_names": fold_names,
        "per_fold_r_orth_cal_56": per_fold_r_orth_56,
        "delta_pred_L1_ens56": delta_l1_ens56,
        "delta_pred_bypass_zero": delta_bypass_zero,
        "delta_pred_mean_batch_emb": delta_mean_batch_emb,
        "delta_pred_multi_batch_ens8": delta_ens8,
        "delta_pred_medoid_idx": delta_medoid,
        "medoid_info": medoid_info,
        "device_used": dev_use,
    }


def run_state_g4_mismatched_r17(
    variant: str,
    target_line: str,
    mismatched_ctrl_tensors: Mapping[str, np.ndarray],
    g4_batch_indices: Sequence[int],
    pert_names: Sequence[str],
    gene_names: Sequence[str],
    pert_canon_to_model: Mapping[str, str],
    gene_canon_map: Mapping[str, str],
    models_dir: Path = DEFAULT_MODELS_DIR,
    status_dir: Path = DEFAULT_R17_STATUS_DIR,
    prereg_lock_path: Optional[Path] = DEFAULT_R17_PREREG_LOCK_PATH,
    device: str = "mps",
    seed: int = DEFAULT_R17_SEED,
    allow_synthetic_without_freeze: bool = False,
) -> Dict[str, Any]:
    """
    Phase 2 (`PREREG_J §3.5` `J5a`): Run ONLY the required unique `g4_batch_indices`
    (`sorted({0, 17, 35, b*_1, b_dag_1, [b*_2, b_dag_2]})`) across all mismatched control
    pools `S`, loading each of the 4 zero-shot folds ONCE.
    Returns `{mismatched_sources, g4_batch_indices, delta_pred_mismatched [S, B_g4, P, G], by_source_and_idx}`.
    """
    assert_zero_training_leakage(variant, target_line, checkpoint_id="zeroshot/ensemble4")
    if not allow_synthetic_without_freeze:
        assert_all_execution_gates_r17(status_dir=status_dir, prereg_lock_path=prereg_lock_path)

    os.environ.setdefault("WANDB_MODE", "disabled")
    os.environ.setdefault("HF_HUB_OFFLINE", "1")
    os.environ.setdefault("TRANSFORMERS_OFFLINE", "1")

    import torch
    from state.tx.models.state_transition import StateTransitionPerturbationModel

    src_keys = sorted(mismatched_ctrl_tensors.keys())
    b_indices = sorted({int(b) for b in g4_batch_indices})
    for b in b_indices:
        if not (0 <= b < N_BATCH_INDICES_STATE):
            raise ValueError(f"Batch index out of range [0..55]: {b}")

    fold_dirs = _resolve_fold_dirs_r17(variant, models_dir=models_dir)
    n_folds = len(fold_dirs)
    n_src = len(src_keys)
    n_b = len(b_indices)
    n_perts = len(pert_names)
    n_genes = len(gene_names)

    pert_keys_model = [pert_canon_to_model.get(p, p) for p in pert_names]
    all_query = ["non-targeting"] + pert_keys_model

    dev_use = device if (device != "mps" or torch.backends.mps.is_available()) else "cpu"
    dev = torch.device(dev_use)

    out_arr = np.zeros((n_src, n_b, n_perts, n_genes), dtype=np.float64)

    for fdir in fold_dirs:
        with open(fdir / "var_dims.pkl", "rb") as fh:
            g_model = list(pickle.load(fh)["gene_names"])
        col_sel = _map_eval_gene_cols_r17(g_model, gene_names, gene_canon_map)
        pert_map = torch.load(fdir / "pert_onehot_map.pt", map_location="cpu", weights_only=False)

        random.seed(int(seed))
        np.random.seed(int(seed))
        torch.manual_seed(int(seed))
        if dev_use == "mps" and hasattr(torch, "mps") and hasattr(torch.mps, "manual_seed"):
            torch.mps.manual_seed(int(seed))

        ckpt_file = fdir / "checkpoints" / "final.ckpt"
        model = StateTransitionPerturbationModel.load_from_checkpoint(str(ckpt_file), map_location="cpu")
        model.eval().to(dev)

        with torch.no_grad():
            for si, sk in enumerate(src_keys):
                ctrl_np = np.asarray(mismatched_ctrl_tensors[sk], dtype=np.float32)
                n_cells = int(ctrl_np.shape[0])
                if n_cells <= 0 or n_cells % 32 != 0:
                    raise ValueError(f"Mismatched pool {sk!r} rows ({n_cells}) must be a multiple of 32.")
                ctrl_dev = torch.from_numpy(ctrl_np).to(dev)
                pert_dev_list = [
                    pert_map[pk].to(device=dev, dtype=torch.float32).unsqueeze(0).expand(n_cells, -1)
                    for pk in all_query
                ]
                for bi_pos, b_val in enumerate(b_indices):
                    b_dev = torch.full((n_cells,), int(b_val), dtype=torch.long, device=dev)
                    pb_ntc = np.empty(n_genes, dtype=np.float64)
                    pb_perts = np.empty((n_perts, n_genes), dtype=np.float64)
                    for qi, p_dev in enumerate(pert_dev_list):
                        pred_raw = model.forward(
                            {"ctrl_cell_emb": ctrl_dev, "pert_emb": p_dev, "batch": b_dev},
                            padded=True,
                        )
                        gene_pred = model.gene_decoder(pred_raw) if variant == "STATE_SE" else pred_raw
                        g_np = np.clip(gene_pred.detach().cpu().numpy().astype(np.float64), 0.0, PRED_CLIP_MAX)
                        mv = g_np.mean(axis=0)[col_sel]
                        if qi == 0:
                            pb_ntc[:] = mv
                        else:
                            pb_perts[qi - 1] = mv
                    out_arr[si, bi_pos] += (pb_perts - pb_ntc[None, :]) / float(n_folds)
                del ctrl_dev, pert_dev_list

            if dev_use == "mps":
                torch.mps.synchronize()

        del model
        gc.collect()
        if dev_use == "mps":
            torch.mps.empty_cache()

    by_source_and_idx: Dict[str, Dict[int, np.ndarray]] = {}
    for si, sk in enumerate(src_keys):
        by_source_and_idx[sk] = {b_val: out_arr[si, bi_pos] for bi_pos, b_val in enumerate(b_indices)}

    return {
        "mismatched_sources": src_keys,
        "g4_batch_indices": b_indices,
        "delta_pred_mismatched": out_arr,
        "by_source_and_idx": by_source_and_idx,
        "device_used": dev_use,
    }
