#!/usr/bin/env python3
"""
Clean-Room Round 17 Runner & Three-Party Tolerance Audit CLI (`round17/audit/run_audit_r17.py`, `A-T3`).

Strictly governed by:
  - `round17/prereg/PREREG_J_DRAFT_v2.md` (`PREREG_J_LOCK.md`)
  - `round17/prereg/PREREG_J_CLARIFICATIONS_C1.md` (`C1.1–C1.13`)
  - `round17/src/R17_ARRAY_SCHEMA.md`

HARD SECURITY & BLINDING BOUNDARIES:
  - NEVER reads `round17/src/*.py` EXCEPT locked `pfizer_telohaec_reader_r17.py`.
  - NEVER reads `round17/results/`, `round16/results/`, `round16/src/`, `round16_addB/src/`,
    `round16_addB/results/`, or any `MAIN.txt` file.
  - `--phase real_r17` REFUSES unless `PREREG_J_LOCK.md`, `R17_CODE_FREEZE_AT3.txt`,
    and `R17_GO_AT3.txt` all exist and are non-empty.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import torch

AUDIT_R17_DIR = Path(__file__).resolve().parent
R17_ROOT = AUDIT_R17_DIR.parent
REPO_ROOT = R17_ROOT.parent
R16_AUDIT_DIR = REPO_ROOT / "round16" / "audit"
R10_SRC_DIR = REPO_ROOT / "round10" / "src"

for _p in (str(AUDIT_R17_DIR), str(R16_AUDIT_DIR), str(R10_SRC_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from audit_r16_raw_controls import SampledRawControlPool  # noqa: E402
from audit_r16_sampling import (  # noqa: E402
    CELL_SET_LEN_DEFAULT,
    MIN_CTRL_FLOOR_DEFAULT,
    N_CTRL_DEFAULT,
)
from audit_r16_state_adapter import (  # noqa: E402
    MODELS_ROOT_DEFAULT,
    CleanRoomArcStateAdapter,
    CleanRoomSE600MEncoder,
    align_hvg_expression_to_model,
    get_hgnc_bridge,
    preprocess_hvg_expression,
    resolve_compute_device,
    resolve_pert_lookup_keys,
)
from audit_r17_controls import (  # noqa: E402
    DEFAULT_GSE295214_RAW_DIR,
    DEFAULT_NADIG_RAW_DIR,
    DEFAULT_PFIZER_H5AD,
    DEFAULT_R17_RESAMPLE_SEED,
    DEFAULT_R17_SEED,
    R16_MISMATCHED_POOLS,
    R17_ALL_EVAL_LINES,
    R17_LINE_MODALITIES,
    R17_LINE_ROLES,
    assert_freeze_gates_r17,
    extract_raw_control_pool_r17,
    load_line_eval_reference_r17,
    resolve_r17_mismatched_strata,
    sha256_file_r17,
    stratified_split_halves_c1,
)
from audit_r17_stats import (  # noqa: E402
    _compute_orth_vectors_and_boot,
    recompute_cohort_from_exported_arrays_r17,
    select_l2_ntc_index_and_diagnostics,
)
from run_audit_r16 import project_model_deltas_to_eval_genes  # noqa: E402

FOLD_CHECKPOINTS_R17: Tuple[str, ...] = (
    "zeroshot/hepg2",
    "zeroshot/jurkat",
    "zeroshot/k562",
    "zeroshot/rpe1",
)
FOLD_SHORT_NAMES_R17: Tuple[str, ...] = ("hepg2", "jurkat", "k562", "rpe1")
ENS8_BATCH_INDICES: Tuple[int, ...] = (0, 7, 14, 21, 28, 35, 42, 49)


def _describe_barcode_prefix_r17(pool_name: str) -> str:
    """Document party barcode prefix convention cleanly per `PREREG_J C2.4`."""
    mapping = {
        "cd4t_rest": "donor_lane:barcode (strip before ':')",
        "cd4t_stim48": "donor_lane:barcode (strip before ':')",
        "mac_lps": "lane_tag:barcode (strip before ':')",
        "mac_ctrl": "lane_tag:barcode (strip before ':')",
        "a549_myllia": "channel:barcode (strip before ':')",
        "calu3": "gem_group:barcode (strip before ':')",
        "a549_ko": "cell_id (exact .obs string; no synthetic prefix)",
        "a549_ko_cas9pos": "cell_id (exact .obs string; no synthetic prefix)",
        "pfizer_telohaec_basal": "exact .obs barcode (no synthetic prefix)",
        "pfizer_telohaec_cytokine": "exact .obs barcode (no synthetic prefix)",
        "hek293_xucao": "exact .obs barcode (no synthetic prefix)",
        "telohaec": "exact .obs barcode (no synthetic prefix)",
        "gse295214_h1": "donor_lane:barcode (strip before ':') or exact .obs string",
        "dld1": "lane_tag:barcode (strip before ':') or exact .obs string",
        "nadig_hepg2": "exact .obs barcode (no synthetic prefix)",
        "nadig_jurkat": "exact .obs barcode (no synthetic prefix)",
        "feng_ipsc": "exact .obs barcode (no synthetic prefix)",
    }
    return mapping.get(pool_name, "pool_prefix:barcode (strip substring through first ':')")


def _json_default_r17(obj: Any) -> Any:
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        v = float(obj)
        return None if (math.isnan(v) or math.isinf(v)) else v
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)):
        return None
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


def _dump_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=False, default=_json_default_r17)
        fh.write("\n")


def _as_strict_str_array(items: Sequence[Any]) -> np.ndarray:
    """Convert string sequences into fixed-width Unicode `dtype=str` arrays (`allow_pickle=False`)."""
    strs = [str(x) for x in items]
    if not strs:
        return np.asarray([], dtype="<U1")
    return np.asarray(strs, dtype=str)


def write_npz_strict_r17(path: Path, arrays: Dict[str, Any]) -> None:
    """
    Save `.npz` strictly via `np.savez_compressed` with zero `object`-dtype arrays
    and verify immediate readability with `np.load(..., allow_pickle=False)`.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    clean: Dict[str, np.ndarray] = {}
    for k, v in arrays.items():
        if isinstance(v, np.ndarray):
            if v.dtype == object:
                clean[k] = _as_strict_str_array(v.tolist()) if v.ndim == 1 else v.astype(str)
            else:
                clean[k] = np.ascontiguousarray(v)
        elif isinstance(v, (str, np.str_)):
            clean[k] = np.asarray(str(v), dtype=str)
        elif isinstance(v, (bool, np.bool_)):
            clean[k] = np.asarray(bool(v), dtype=bool)
        elif isinstance(v, (int, np.integer)):
            clean[k] = np.asarray(int(v), dtype=np.int64)
        elif isinstance(v, (float, np.floating)):
            clean[k] = np.asarray(float(v), dtype=np.float64)
        elif isinstance(v, (list, tuple)):
            if len(v) == 0:
                clean[k] = np.asarray([], dtype="<U1")
            elif all(isinstance(x, (str, np.str_)) for x in v):
                clean[k] = _as_strict_str_array(v)
            elif all(isinstance(x, (bool, np.bool_)) for x in v):
                clean[k] = np.asarray(v, dtype=bool)
            elif all(isinstance(x, (int, np.integer)) and not isinstance(x, bool) for x in v):
                clean[k] = np.asarray(v, dtype=np.int64)
            else:
                clean[k] = np.asarray(v, dtype=np.float64)
        else:
            raise TypeError(f"Unsupported type for strict NPZ key {k!r}: {type(v).__name__}")
        if clean[k].dtype == object:
            raise ValueError(f"Key {k!r} resolved to object dtype, violating allow_pickle=False")

    np.savez_compressed(path, **clean)
    with np.load(path, allow_pickle=False) as check:
        _ = check.files


def build_manifest_sha256(directory: Path, pattern: str = "*.npz") -> Dict[str, str]:
    files = sorted(p for p in directory.glob(pattern) if p.is_file())
    return {p.name: sha256_file_r17(p) for p in files}


def compute_batch_encoder_medoids(weights_4_folds: np.ndarray) -> Dict[str, Any]:
    """
    Compute Euclidean geometric medoid (`C1.9` primary) and Cosine geometric medoid
    (descriptive diagnostic) across the 4-fold mean `batch_encoder.weight` (`(56, D)`).
    Ties broken strictly by lowest index (`np.argmin`).
    """
    arr = np.asarray(weights_4_folds, dtype=np.float64)
    if arr.ndim != 3 or arr.shape[1] != 56:
        raise ValueError(f"Expected weights_4_folds shape (4, 56, D), got {arr.shape}")
    w_mean = arr.mean(axis=0)  # (56, D)
    diff = w_mean[:, None, :] - w_mean[None, :, :]
    d_euc = np.linalg.norm(diff, axis=2)
    sum_euc = d_euc.sum(axis=1)
    euc_idx = int(np.argmin(sum_euc))

    norms = np.maximum(np.linalg.norm(w_mean, axis=1, keepdims=True), 1e-12)
    u = w_mean / norms
    cos_sim = np.clip(u @ u.T, -1.0, 1.0)
    d_cos = np.clip(1.0 - cos_sim, 0.0, 2.0)
    sum_cos = d_cos.sum(axis=1)
    cos_idx = int(np.argmin(sum_cos))

    per_fold_euc: List[int] = []
    for f in range(arr.shape[0]):
        df = np.linalg.norm(arr[f][:, None, :] - arr[f][None, :, :], axis=2)
        per_fold_euc.append(int(np.argmin(df.sum(axis=1))))

    return {
        "euclid_medoid_idx": euc_idx,
        "cosine_medoid_idx": cos_idx,
        "euclid_sum_distances_56": sum_euc.tolist(),
        "cosine_sum_distances_56": sum_cos.tolist(),
        "per_fold_euclid_medoid_idx": per_fold_euc,
    }


def save_x_state_embedding_npz(
    path: Path,
    x_state: np.ndarray,
    control_barcodes: Sequence[str],
    seed: int,
    dataset_name: str,
    se_audit: Dict[str, Any],
) -> None:
    arr = np.ascontiguousarray(np.asarray(x_state, dtype=np.float32))
    if arr.ndim != 2 or arr.shape[1] != 2058:
        raise ValueError(f"Expected X_state shape (n_ctrl, 2058), got {arr.shape}")
    write_npz_strict_r17(
        path,
        {
            "X_state": arr,
            "control_barcodes": _as_strict_str_array(control_barcodes),
            "seed": int(seed),
            "dataset_name": str(dataset_name),
            "se_audit_json": json.dumps(se_audit, sort_keys=True, default=_json_default_r17),
        },
    )


def build_observed_npz_dict(
    ref_bundle: Dict[str, Any],
    raw_pool: SampledRawControlPool,
) -> Dict[str, Any]:
    n_total = int(raw_pool.sampling_audit.get("n_valid_pool", len(raw_pool.batch_labels)))
    n_sampled = int(raw_pool.raw_counts.shape[0])
    frac = float(n_sampled) / float(n_total) if n_total > 0 else float("nan")
    return {
        "line_name": str(ref_bundle["line_name"]),
        "role": str(ref_bundle["role"]),
        "modality": str(ref_bundle["modality"]),
        "checkpoint_id": str(ref_bundle["checkpoint_id"]),
        "delta_obs": np.asarray(ref_bundle["delta_obs"], dtype=np.float64),
        "delta_half_a": np.asarray(ref_bundle["delta_half_a"], dtype=np.float64),
        "delta_half_b": np.asarray(ref_bundle["delta_half_b"], dtype=np.float64),
        "half_valid": np.asarray(ref_bundle["half_valid"], dtype=bool),
        "delta_transfer": np.asarray(ref_bundle["delta_transfer"], dtype=np.float64),
        "source_grand_mean": np.asarray(ref_bundle["source_grand_mean"], dtype=np.float64),
        "fct_delta_pred": np.asarray(ref_bundle["fct_delta_pred"], dtype=np.float64),
        "e1_constmed_delta_pred": np.asarray(ref_bundle["e1_constmed_delta_pred"], dtype=np.float64),
        "target_indices": np.asarray(ref_bundle["target_indices"], dtype=np.int64),
        "target_names": _as_strict_str_array(ref_bundle["eval_perts"]),
        "gene_names": _as_strict_str_array(ref_bundle["eval_genes"]),
        "mu_target": np.asarray(ref_bundle["mu_target"], dtype=np.float64),
        "mu_sources": np.asarray(ref_bundle["mu_sources"], dtype=np.float64),
        "obs_ntc_mean_official_median_log1p": np.asarray(
            ref_bundle["obs_ntc_mean_official_median_log1p"], dtype=np.float64
        ),
        "split_s0_p1_indices": np.asarray(ref_bundle["split_s0_p1_indices"], dtype=np.int64),
        "split_s0_p2_indices": np.asarray(ref_bundle["split_s0_p2_indices"], dtype=np.int64),
        "seen_sources": _as_strict_str_array(ref_bundle["seen_sources"]),
        "unseen_sources": _as_strict_str_array(ref_bundle["unseen_sources"]),
        "same_lineage_sources": _as_strict_str_array(ref_bundle["same_lineage_sources"]),
        "matched_ctrl_barcodes": _as_strict_str_array(raw_pool.cell_barcodes.tolist()),
        "matched_ctrl_batch_labels": _as_strict_str_array(raw_pool.batch_labels.tolist()),
        "n_ctrl_total": int(n_total),
        "n_ctrl_sampled": int(n_sampled),
        "ctrl_sampling_fraction": float(frac),
    }


class RealStateSweepBackend:
    """
    Wraps clean-room `CleanRoomSE600MEncoder` & `CleanRoomArcStateAdapter` (`zeroshot/ensemble4`)
    for full 56-index + bypass/mean_emb/mismatched inference without reading `round17/src/`.
    """

    def __init__(
        self,
        models_root: Path = MODELS_ROOT_DEFAULT,
        device: str = "auto",
        hgnc_bridge: Optional[Any] = None,
        status_dir: Optional[Path] = None,
    ) -> None:
        self.models_root = Path(models_root)
        self.device = resolve_compute_device(device)
        self.hgnc_bridge = get_hgnc_bridge(hgnc_bridge)
        self.status_dir = Path(status_dir) if status_dir is not None else (REPO_ROOT / "round16" / "status")
        self._se_encoder: Optional[CleanRoomSE600MEncoder] = None
        self._adapters: Dict[str, CleanRoomArcStateAdapter] = {}
        self._hvg_ctrl_cache: Dict[Tuple[str, int], Tuple[np.ndarray, Dict[str, Any]]] = {}

    def get_se_encoder(self) -> CleanRoomSE600MEncoder:
        if self._se_encoder is None:
            self._se_encoder = CleanRoomSE600MEncoder(
                se600m_dir=self.models_root / "SE-600M",
                device=self.device,
                batch_size=16,
                hgnc_bridge=self.hgnc_bridge,
            )
            self._se_encoder.load()
        return self._se_encoder

    def get_adapter(self, model_key: str) -> CleanRoomArcStateAdapter:
        if model_key not in self._adapters:
            ad = CleanRoomArcStateAdapter(
                model_key=model_key,
                checkpoint_id="zeroshot/ensemble4",
                models_root=self.models_root,
                status_dir=self.status_dir,
                device=self.device,
                neutral_batch_mode="index_zero",
                delta_mode="minus_pred_ntc",
                hvg_preproc_mode="official_median_log1p",
                allow_synthetic_smoke=False,
                hgnc_bridge=self.hgnc_bridge,
            )
            ad.load_weights()
            self._adapters[model_key] = ad
        return self._adapters[model_key]

    def get_model_genes_and_perts(self, model_key: str) -> Tuple[List[str], List[str]]:
        ad = self.get_adapter(model_key)
        common_genes: Optional[set[str]] = None
        common_perts: Optional[set[str]] = None
        first_genes: List[str] = []
        for ck in FOLD_SHORT_NAMES_R17:
            b = ad._bundles[ck]
            gs = [self.hgnc_bridge.canonicalize(str(g)) or str(g) for g in b["gene_names"]]
            ps = [self.hgnc_bridge.canonicalize(str(p)) or str(p) for p in b["pert_lookup"].keys() if p != b["control_pert"]]
            if common_genes is None:
                first_genes = list(b["gene_names"])
                common_genes = set(gs)
                common_perts = set(ps)
            else:
                common_genes &= set(gs)
                common_perts &= set(ps)
        ordered_genes = [g for g in first_genes if (self.hgnc_bridge.canonicalize(str(g)) or str(g)) in (common_genes or set())]
        return ordered_genes, sorted(common_perts or set())

    def encode_pool_se600m(
        self,
        raw_pool: SampledRawControlPool,
        seed: int,
    ) -> Tuple[np.ndarray, Dict[str, Any]]:
        enc = self.get_se_encoder()
        x_state, map_stats = enc.encode_cells(
            ctrl_counts=raw_pool.raw_counts,
            source_genes=raw_pool.gene_names,
            source_ensembl_ids=raw_pool.gene_ids,
            seed=int(seed),
            use_hgnc_rescue=False,
            dataset_name=raw_pool.line_name,
        )
        audit = dict(map_stats)
        audit["x_state_sha256"] = hashlib.sha256(np.ascontiguousarray(x_state, dtype=np.float32).tobytes()).hexdigest()
        audit["device"] = self.device
        audit["dtype"] = "float32"
        return x_state, audit

    def _prepare_fold_ctrl_tensor(
        self,
        model_key: str,
        fold_short: str,
        raw_pool: SampledRawControlPool,
        x_state_pool: Optional[np.ndarray],
        n_sets: int,
        cell_set_len: int,
    ) -> torch.Tensor:
        ad = self.get_adapter(model_key)
        bundle = ad._bundles[fold_short]
        dev = bundle["device"]
        if model_key == "STATE_SE":
            if x_state_pool is None:
                raise ValueError(f"Missing cached X_state for pool {raw_pool.line_name!r}")
            arr = np.asarray(x_state_pool, dtype=np.float32).reshape(n_sets, cell_set_len, 2058)
            return torch.tensor(arr, dtype=torch.float32, device=dev)

        fold_genes = list(bundle["gene_names"])
        cache_key = (raw_pool.line_name, hash(tuple(fold_genes)))
        if cache_key not in self._hvg_ctrl_cache:
            raw_mat = np.asarray(raw_pool.raw_counts, dtype=np.float32)
            raw_col_means = np.mean(raw_mat.astype(np.float64), axis=0) if raw_mat.shape[0] > 0 else None
            norm_mat = preprocess_hvg_expression(
                raw_mat,
                mode="official_median_log1p",
                cell_total_umis=raw_pool.cell_total_umis,
                target_median_umi=raw_pool.median_ntc_umi,
            )
            aligned, ainfo = align_hvg_expression_to_model(
                norm_mat,
                raw_pool.gene_names,
                fold_genes,
                source_ensembl_ids=raw_pool.gene_ids,
                hgnc_bridge=self.hgnc_bridge,
                col_mean_scores=raw_col_means,
            )
            self._hvg_ctrl_cache[cache_key] = (np.asarray(aligned, dtype=np.float32), ainfo)
        aligned_arr, _ = self._hvg_ctrl_cache[cache_key]
        arr = aligned_arr.reshape(n_sets, cell_set_len, len(fold_genes))
        return torch.tensor(arr, dtype=torch.float32, device=dev)

    def get_batch_encoder_weights(self, model_key: str) -> np.ndarray:
        ad = self.get_adapter(model_key)
        ws: List[np.ndarray] = []
        for ck in FOLD_SHORT_NAMES_R17:
            b = ad._bundles[ck]
            ws.append(b["model"].batch_encoder.weight.detach().cpu().numpy().astype(np.float64))
        return np.stack(ws, axis=0)  # (4, 56, D)

    def predict_batch_indices_for_pool(
        self,
        model_key: str,
        raw_pool: SampledRawControlPool,
        x_state_pool: Optional[np.ndarray],
        eval_perts: Sequence[str],
        eval_genes: Sequence[str],
        batch_schemes: Sequence[Tuple[str, int]],
        cell_set_len: int = CELL_SET_LEN_DEFAULT,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Evaluate exact `(B_schemes,)` batch schemes across all 4 ensemble folds on `raw_pool`.
        Returns:
          - `delta_per_fold`: `(4, B_schemes, P, G_eval)` `float64`
          - `ntc_pb_per_fold`: `(4, B_schemes, G_eval)` `float64`
        """
        ad = self.get_adapter(model_key)
        n_ctrl = int(raw_pool.raw_counts.shape[0])
        if n_ctrl % int(cell_set_len) != 0:
            raise ValueError(f"Pool {raw_pool.line_name} n_ctrl={n_ctrl} not divisible by cell_set_len={cell_set_len}")
        n_sets = n_ctrl // int(cell_set_len)
        n_b = len(batch_schemes)
        p_cnt = len(eval_perts)
        g_cnt = len(eval_genes)

        delta_pf = np.empty((len(FOLD_SHORT_NAMES_R17), n_b, p_cnt, g_cnt), dtype=np.float64)
        ntc_pf = np.empty((len(FOLD_SHORT_NAMES_R17), n_b, g_cnt), dtype=np.float64)

        for f_i, ck in enumerate(FOLD_SHORT_NAMES_R17):
            bundle = ad._bundles[ck]
            model = bundle["model"]
            dev = bundle["device"]
            pert_lookup = bundle["pert_lookup"]
            ntc_key = bundle["control_pert"]
            fold_genes = list(bundle["gene_names"])
            n_genes = len(fold_genes)
            p_map = resolve_pert_lookup_keys(
                eval_perts, pert_lookup, control_pert=ntc_key, hgnc_bridge=self.hgnc_bridge
            )
            query_raw_keys = [p_map[str(p).strip()] for p in eval_perts] + [ntc_key]
            ctrl_tensor = self._prepare_fold_ctrl_tensor(
                model_key=model_key,
                fold_short=ck,
                raw_pool=raw_pool,
                x_state_pool=x_state_pool,
                n_sets=n_sets,
                cell_set_len=int(cell_set_len),
            )
            for bi, (bmode, bidx) in enumerate(batch_schemes):
                preds_q_cells = ad._predict_single_batch_index(
                    model=model,
                    dev=dev,
                    pert_lookup=pert_lookup,
                    query_raw_keys=query_raw_keys,
                    ctrl_tensor=ctrl_tensor,
                    batch_mode=bmode,
                    fixed_batch_idx=int(bidx),
                    n_genes=n_genes,
                    per_cell_batch_idx=None,
                )  # (P + 1, n_cells, n_genes) float32 clipped [0, 14]
                p_pb = preds_q_cells[:-1].astype(np.float64).mean(axis=1)  # (P, n_genes)
                ntc_pb = preds_q_cells[-1].astype(np.float64).mean(axis=0)  # (n_genes,)
                d_fold_full = p_pb - ntc_pb[None, :]
                delta_pf[f_i, bi] = project_model_deltas_to_eval_genes(
                    d_fold_full, fold_genes, eval_genes, hgnc_bridge=self.hgnc_bridge
                )
                ntc_proj = project_model_deltas_to_eval_genes(
                    ntc_pb[None, :], fold_genes, eval_genes, hgnc_bridge=self.hgnc_bridge
                )[0]
                ntc_pf[f_i, bi] = ntc_proj

        return delta_pf, ntc_pf


class SyntheticStateSweepBackend:
    """
    Pure-CPU deterministic synthetic State backend implementing the exact 4-fold
    `zeroshot/{hepg2,jurkat,k562,rpe1}` + `56`-index `batch_encoder.weight` interface
    for `--phase synthetic_e2e` and CPU unit verification (ZERO real data or weight I/O).
    """

    def __init__(
        self,
        model_genes: Sequence[str],
        eval_perts: Sequence[str],
        seed: int = DEFAULT_R17_SEED,
    ) -> None:
        self.model_genes = list(model_genes)
        self.eval_perts = list(eval_perts)
        self.seed = int(seed)
        rng = np.random.default_rng(self.seed + 77)
        g_cnt = len(self.model_genes)
        p_cnt = len(self.eval_perts)

        # Deterministic perturbation effects & batch vectors per model and fold
        self.pert_effects: Dict[str, np.ndarray] = {
            mk: rng.normal(0.0, 0.25, size=(4, p_cnt, g_cnt)).astype(np.float64)
            for mk in ("STATE_SE", "STATE_HVG")
        }
        self.batch_weights: Dict[str, np.ndarray] = {}
        self.batch_gene_proj: Dict[str, np.ndarray] = {}
        for mk_i, mk in enumerate(("STATE_SE", "STATE_HVG")):
            mrng = np.random.default_rng(self.seed + 500 + mk_i * 101)
            w = mrng.normal(0.0, 0.4, size=(4, 56, 128)).astype(np.float64)
            # Anchor clusters in 128-dim weight space cleanly so Euclidean medoid is unambiguous
            w[:, 10, :] *= 0.15
            w[:, 17, :] += 1.4
            w[:, 35, :] -= 1.1
            self.batch_weights[mk] = w
            self.batch_gene_proj[mk] = mrng.normal(0.0, 0.08, size=(128, g_cnt)).astype(np.float64)

        self._gene_to_idx = {g: i for i, g in enumerate(self.model_genes)}
        self._pert_to_idx = {p: i for i, p in enumerate(self.eval_perts)}

    def get_model_genes_and_perts(self, model_key: str) -> Tuple[List[str], List[str]]:
        return list(self.model_genes), list(self.eval_perts)

    def encode_pool_se600m(
        self,
        raw_pool: SampledRawControlPool,
        seed: int,
    ) -> Tuple[np.ndarray, Dict[str, Any]]:
        prng = np.random.default_rng(int(seed) + sum(ord(c) for c in raw_pool.line_name) * 13)
        base = np.mean(np.log1p(raw_pool.raw_counts.astype(np.float32)), axis=1, keepdims=True)
        n_ctrl = raw_pool.raw_counts.shape[0]
        x_state = (base + prng.normal(0.0, 0.25, size=(n_ctrl, 2058))).astype(np.float32)
        audit = {
            "source": "synthetic_se600m",
            "device": "cpu",
            "dtype": "float32",
            "n_mapped_to_pe19790": int(raw_pool.raw_counts.shape[1]),
            "x_state_sha256": hashlib.sha256(np.ascontiguousarray(x_state, dtype=np.float32).tobytes()).hexdigest(),
        }
        return x_state, audit

    def get_batch_encoder_weights(self, model_key: str) -> np.ndarray:
        return self.batch_weights[model_key].copy()

    def predict_batch_indices_for_pool(
        self,
        model_key: str,
        raw_pool: SampledRawControlPool,
        x_state_pool: Optional[np.ndarray],
        eval_perts: Sequence[str],
        eval_genes: Sequence[str],
        batch_schemes: Sequence[Tuple[str, int]],
        cell_set_len: int = CELL_SET_LEN_DEFAULT,
    ) -> Tuple[np.ndarray, np.ndarray]:
        norm_c = preprocess_hvg_expression(
            np.asarray(raw_pool.raw_counts, dtype=np.float32),
            mode="official_median_log1p",
            cell_total_umis=raw_pool.cell_total_umis,
            target_median_umi=raw_pool.median_ntc_umi,
        ).astype(np.float64)
        ctrl_mean_in = norm_c.mean(axis=0)  # matches model_genes order in synthetic mode
        se_shift = float(np.mean(x_state_pool)) * 0.05 if x_state_pool is not None else 0.0

        p_idx = np.asarray([self._pert_to_idx[p] for p in eval_perts], dtype=np.int64)
        g_idx = np.asarray([self._gene_to_idx[g] for g in eval_genes], dtype=np.int64)
        w4 = self.batch_weights[model_key]
        proj = self.batch_gene_proj[model_key][:, g_idx]  # (128, G_eval)

        n_b = len(batch_schemes)
        p_cnt = len(eval_perts)
        g_cnt = len(eval_genes)
        delta_pf = np.empty((4, n_b, p_cnt, g_cnt), dtype=np.float64)
        ntc_pf = np.empty((4, n_b, g_cnt), dtype=np.float64)

        # Control scale sensitivity makes cd4t_stim48 / pfizer batch-sensitive while keeping other pools distinct
        ctrl_scale = float(np.std(ctrl_mean_in[g_idx])) + 0.35
        for f_i in range(4):
            p_base = self.pert_effects[model_key][f_i][np.ix_(p_idx, g_idx)]
            w_f = w4[f_i]
            for bi, (bmode, bidx) in enumerate(batch_schemes):
                if bmode == "bypass_zero":
                    b_vec = np.zeros(128, dtype=np.float64)
                elif bmode == "mean_batch_emb":
                    b_vec = w_f.mean(axis=0)
                else:
                    b_vec = w_f[int(bidx)]

                b_gene = b_vec @ proj  # (G_eval,)
                # Perturbation * batch interaction drives true r_orth_cal variance across b = 0..55
                b_amp = 1.0 + 0.45 * float(np.tanh(b_vec[:8].mean())) * ctrl_scale
                if model_key == "STATE_HVG":
                    b_amp = 1.0 + 0.12 * float(np.tanh(b_vec[:8].mean()))
                ntc_pb = np.clip(ctrl_mean_in[g_idx] + 0.25 * b_gene + se_shift, 0.0, 14.0)
                p_pb = np.clip(ntc_pb[None, :] + p_base * b_amp + 0.10 * b_gene[None, :], 0.0, 14.0)
                ntc_pf[f_i, bi] = ntc_pb
                delta_pf[f_i, bi] = p_pb - ntc_pb[None, :]

        return delta_pf, ntc_pf


def export_line_model_arrays_r17(
    line_dir: Path,
    model_key: str,
    backend: Any,
    ref_bundle: Dict[str, Any],
    obs_npz: Dict[str, Any],
    raw_pool_matched: SampledRawControlPool,
    x_state_matched: Optional[np.ndarray],
    mismatched_pools: Dict[str, SampledRawControlPool],
    x_state_mismatched: Dict[str, Optional[np.ndarray]],
    cell_set_len: int = CELL_SET_LEN_DEFAULT,
    seed: int = DEFAULT_R17_SEED,
) -> Dict[str, Any]:
    """
    Run full 56-index sweep + comparator arms + `G4` mismatched evaluation for one
    `(target_line, model_key)` and export `preds_<M>_sweep56.npz` & `preds_<M>_g4_mismatched.npz`.
    """
    eval_perts = ref_bundle["eval_perts"]
    eval_genes = ref_bundle["eval_genes"]
    d_obs = obs_npz["delta_obs"]
    ha = obs_npz["delta_half_a"]
    hb = obs_npz["delta_half_b"]
    hv = obs_npz["half_valid"]
    tm = obs_npz["delta_transfer"]
    tg = obs_npz["target_indices"]
    obs_ntc_mean = obs_npz["obs_ntc_mean_official_median_log1p"]

    sweep_schemes: List[Tuple[str, int]] = [("index_zero", b) for b in range(56)] + [
        ("bypass_zero", 0),
        ("mean_batch_emb", 0),
    ]
    d_pf, ntc_pf = backend.predict_batch_indices_for_pool(
        model_key=model_key,
        raw_pool=raw_pool_matched,
        x_state_pool=x_state_matched,
        eval_perts=eval_perts,
        eval_genes=eval_genes,
        batch_schemes=sweep_schemes,
        cell_set_len=cell_set_len,
    )

    d_pf_56 = d_pf[:, :56, :, :]  # (4, 56, P, G)
    ntc_pf_56 = ntc_pf[:, :56, :]  # (4, 56, G)
    delta_pred_56 = d_pf_56.mean(axis=0)  # (56, P, G) float64
    ntc_pb_ensemble_56 = ntc_pf_56.mean(axis=0)  # (56, G) float64
    delta_bypass_zero = d_pf[:, 56, :, :].mean(axis=0)
    delta_mean_batch_emb = d_pf[:, 57, :, :].mean(axis=0)
    delta_multi_batch_ens8 = delta_pred_56[list(ENS8_BATCH_INDICES)].mean(axis=0)
    delta_l1 = delta_pred_56.mean(axis=0)

    # Per-fold 56-index point r_orth_cal diagnostic (R17_ARRAY_SCHEMA.md §3.2)
    per_fold_r56 = np.empty((4, 56), dtype=np.float64)
    for f_i in range(4):
        for b in range(56):
            per_fold_r56[f_i, b] = _compute_orth_vectors_and_boot(
                d_obs, d_pf_56[f_i, b], tm, tg, ha, hb, hv, draws=None, source_grand_mean=None
            )["estimate"]

    # L2 selection via C_NTC(b)
    l2_diag = select_l2_ntc_index_and_diagnostics(ntc_pb_ensemble_56, obs_ntc_mean, ntc_pb_per_fold_56=ntc_pf_56)
    l2_idx = int(l2_diag["top1_idx"])
    delta_l2 = delta_pred_56[l2_idx]

    # Weight-only Euclidean & Cosine geometric medoids (C1.9)
    w4 = backend.get_batch_encoder_weights(model_key)
    medoid_info = compute_batch_encoder_medoids(w4)
    medoid_idx = int(medoid_info["euclid_medoid_idx"])
    delta_medoid = delta_pred_56[medoid_idx]

    sweep_path = line_dir / f"preds_{model_key}_sweep56.npz"
    write_npz_strict_r17(
        sweep_path,
        {
            "delta_pred_56": np.asarray(delta_pred_56, dtype=np.float64),
            "batch_indices": np.arange(56, dtype=np.int64),
            "ntc_pb_ensemble_56": np.asarray(ntc_pb_ensemble_56, dtype=np.float64),
            "ntc_pb_per_fold_56": np.asarray(ntc_pf_56, dtype=np.float64),
            "fold_names": _as_strict_str_array(FOLD_CHECKPOINTS_R17),
            "per_fold_r_orth_cal_56": np.asarray(per_fold_r56, dtype=np.float64),
            "delta_pred_L1_ens56": np.asarray(delta_l1, dtype=np.float64),
            "delta_pred_L2_ntc_recon_match": np.asarray(delta_l2, dtype=np.float64),
            "l2_selected_batch_idx": int(l2_idx),
            "delta_pred_bypass_zero": np.asarray(delta_bypass_zero, dtype=np.float64),
            "delta_pred_mean_batch_emb": np.asarray(delta_mean_batch_emb, dtype=np.float64),
            "delta_pred_multi_batch_ens8": np.asarray(delta_multi_batch_ens8, dtype=np.float64),
            "delta_pred_medoid_idx": np.asarray(delta_medoid, dtype=np.float64),
            "medoid_batch_idx": int(medoid_idx),
            "target_names": _as_strict_str_array(eval_perts),
            "gene_names": _as_strict_str_array(eval_genes),
            "model_key": str(model_key),
            "checkpoint_id": "zeroshot/ensemble4",
            "hvg_norm_mode": "official_median_log1p" if model_key == "STATE_HVG" else "",
            "control_source": str(ref_bundle["line_name"]),
            "control_barcodes": _as_strict_str_array(raw_pool_matched.cell_barcodes.tolist()),
        },
    )

    # Stage 3 G4 Mismatched Control Evaluation (`s = 0` cross-fitted pairs + continuity (0,17), (0,35))
    p1_idx = np.asarray(obs_npz["split_s0_p1_indices"], dtype=np.int64)
    p2_idx = np.asarray(obs_npz["split_s0_p2_indices"], dtype=np.int64)
    g4_set = {0, 17, 35}
    if p1_idx.size >= 2 and p2_idx.size >= 2:
        rp1 = np.asarray(
            [
                _compute_orth_vectors_and_boot(
                    d_obs[p1_idx], delta_pred_56[b, p1_idx], tm[p1_idx], tg[p1_idx], ha[p1_idx], hb[p1_idx], hv[p1_idx],
                    draws=None, source_grand_mean=None,
                )["estimate"]
                for b in range(56)
            ],
            dtype=np.float64,
        )
        rp2 = np.asarray(
            [
                _compute_orth_vectors_and_boot(
                    d_obs[p2_idx], delta_pred_56[b, p2_idx], tm[p2_idx], tg[p2_idx], ha[p2_idx], hb[p2_idx], hv[p2_idx],
                    draws=None, source_grand_mean=None,
                )["estimate"]
                for b in range(56)
            ],
            dtype=np.float64,
        )
        b_star_1 = int(np.argmax(rp1))
        b_dag_1 = int(np.argmin(rp1))
        b_star_2 = int(np.argmax(rp2))
        b_dag_2 = int(np.argmin(rp2))
        g4_set.update({b_star_1, b_dag_1, b_star_2, b_dag_2})
    else:
        b_star_1 = b_dag_1 = b_star_2 = b_dag_2 = -1

    g4_batch_indices = sorted(g4_set)
    seen_s = list(ref_bundle["seen_sources"])
    unseen_s = list(ref_bundle["unseen_sources"])
    same_s = list(ref_bundle["same_lineage_sources"])
    mismatched_sources = sorted(set(seen_s) | set(unseen_s) | set(same_s))
    g4_schemes = [("index_zero", b) for b in g4_batch_indices]

    dmis = np.empty((len(mismatched_sources), len(g4_batch_indices), len(eval_perts), len(eval_genes)), dtype=np.float64)
    for si, sname in enumerate(mismatched_sources):
        mpool = mismatched_pools[sname]
        mx = x_state_mismatched.get(sname)
        d_s_pf, _ = backend.predict_batch_indices_for_pool(
            model_key=model_key,
            raw_pool=mpool,
            x_state_pool=mx,
            eval_perts=eval_perts,
            eval_genes=eval_genes,
            batch_schemes=g4_schemes,
            cell_set_len=cell_set_len,
        )
        dmis[si] = d_s_pf.mean(axis=0)  # (B_g4, P, G)

    g4_path = line_dir / f"preds_{model_key}_g4_mismatched.npz"
    write_npz_strict_r17(
        g4_path,
        {
            "mismatched_sources": _as_strict_str_array(mismatched_sources),
            "seen_sources": _as_strict_str_array(seen_s),
            "unseen_sources": _as_strict_str_array(unseen_s),
            "same_lineage_sources": _as_strict_str_array(same_s),
            "g4_batch_indices": np.asarray(g4_batch_indices, dtype=np.int64),
            "delta_pred_mismatched": np.asarray(dmis, dtype=np.float64),
            "b_star_1": int(b_star_1),
            "b_dag_1": int(b_dag_1),
            "b_star_2": int(b_star_2),
            "b_dag_2": int(b_dag_2),
            "target_names": _as_strict_str_array(eval_perts),
            "gene_names": _as_strict_str_array(eval_genes),
            "model_key": str(model_key),
            "checkpoint_id": "zeroshot/ensemble4",
        },
    )
    return {
        "medoid_diagnostics": medoid_info,
        "l2_selected_batch_idx": int(l2_idx),
        "g4_batch_indices": g4_batch_indices,
    }


def _collect_numeric_leaves(obj: Any, prefix: str = "") -> Dict[str, float]:
    out: Dict[str, float] = {}
    if isinstance(obj, dict):
        for k, v in obj.items():
            p = f"{prefix}.{k}" if prefix else str(k)
            out.update(_collect_numeric_leaves(v, p))
    elif isinstance(obj, (list, tuple)):
        for i, v in enumerate(obj):
            out.update(_collect_numeric_leaves(v, f"{prefix}[{i}]"))
    elif isinstance(obj, (int, float, np.integer, np.floating)) and not isinstance(obj, bool):
        out[prefix] = float(obj)
    return out


def compare_summaries_stage4(summary_a: Dict[str, Any], summary_b: Dict[str, Any], tol: float = 1e-10) -> Dict[str, Any]:
    """Compare intersecting numeric leaves of `recompute_cohort_from_exported_arrays_r17` vs saved summary."""
    shared_roots = (
        "n_boot",
        "n_boot_bonf56",
        "n_perm",
        "n_splits_diag",
        "seed",
        "models",
        "holm_by_model",
        "holm_pooled_sensitivity",
        "claim_rollups",
        "j5b_rank_correlations",
        "sens_gse295214_h1_full_r9_6",
    )
    la = _collect_numeric_leaves({k: summary_a[k] for k in shared_roots if k in summary_a})
    lb = _collect_numeric_leaves({k: summary_b[k] for k in shared_roots if k in summary_b})
    if set(la.keys()) != set(lb.keys()):
        missing_a = sorted(set(lb.keys()) - set(la.keys()))[:5]
        missing_b = sorted(set(la.keys()) - set(lb.keys()))[:5]
        raise ValueError(f"Summary key mismatch: missing_in_a={missing_a}, missing_in_b={missing_b}")

    max_abs = 0.0
    worst_key = ""
    for k, va in la.items():
        vb = lb[k]
        if math.isnan(va) and math.isnan(vb):
            continue
        if math.isnan(va) != math.isnan(vb):
            raise ValueError(f"NaN mismatch at leaf {k}: {va} vs {vb}")
        diff = abs(va - vb)
        if diff > max_abs:
            max_abs = diff
            worst_key = k
    if max_abs > float(tol):
        raise ValueError(f"Stage 4 parity violation at {worst_key}: max_abs={max_abs:.3e} > tol={tol:.3e}")
    return {"n_leaves": len(la), "max_abs_diff": float(max_abs), "worst_key": worst_key, "passed": True}


def _resolve_target_lines(cohort: str, lines_arg: Optional[Sequence[str]]) -> List[str]:
    if lines_arg:
        for ln in lines_arg:
            if ln not in R17_ALL_EVAL_LINES:
                raise ValueError(f"Unknown Round 17 line {ln!r}; allowed: {R17_ALL_EVAL_LINES}")
        return list(lines_arg)
    if cohort == "all":
        return list(R17_ALL_EVAL_LINES)
    if cohort == "confirmatory":
        return [ln for ln in R17_ALL_EVAL_LINES if R17_LINE_ROLES[ln] == "confirmatory"]
    raise ValueError(f"Unsupported --cohort {cohort!r}")


def run_real_r17(
    target_lines: Sequence[str],
    models: Sequence[str],
    output_dir: Path,
    models_root: Path,
    status_dir: Optional[Path],
    prereg_dir: Optional[Path],
    raw_nadig_dir: Path,
    raw_gse295214_dir: Path,
    pfizer_h5ad_path: Path,
    device: str = "auto",
    n_ctrl: int = N_CTRL_DEFAULT,
    cell_set_len: int = CELL_SET_LEN_DEFAULT,
    ctrl_seed: int = DEFAULT_R17_SEED,
    ctrl_resample_seed: int = DEFAULT_R17_RESAMPLE_SEED,
    n_boot: int = 4000,
    n_boot_bonf56: int = 20000,
    n_perm: int = 10000,
    n_splits_diag: int = 20,
    seed: int = DEFAULT_R17_SEED,
) -> Dict[str, Any]:
    """
    Execute real-cell Round 17 confirmatory/replication/descriptive inference on MPS/CPU.
    STRICTLY refuses unless `PREREG_J_LOCK.md`, `R17_CODE_FREEZE_AT3.txt`, and `R17_GO_AT3.txt` exist.
    """
    gate_receipts = assert_freeze_gates_r17(status_dir=status_dir, prereg_dir=prereg_dir)

    out_dir = Path(output_dir)
    emb_dir = out_dir / "embeddings"
    arrays_dir = out_dir / "arrays"
    out_dir.mkdir(parents=True, exist_ok=True)
    emb_dir.mkdir(parents=True, exist_ok=True)
    arrays_dir.mkdir(parents=True, exist_ok=True)

    hb = get_hgnc_bridge()
    backend = RealStateSweepBackend(
        models_root=models_root,
        device=device,
        hgnc_bridge=hb,
        status_dir=status_dir,
    )
    for mkey in models:
        backend.get_adapter(mkey)
    model_genes, model_perts = backend.get_model_genes_and_perts("STATE_SE" if "STATE_SE" in models else models[0])
    # Use intersection of available genes/perts across requested models so all models evaluate aligned R17 sets
    if len(models) > 1:
        g_sets = []
        p_sets = []
        for mkey in models:
            mg, mp = backend.get_model_genes_and_perts(mkey)
            g_sets.append(set(hb.canonicalize(x) or x for x in mg))
            p_sets.append(set(mp))
        inter_g = set.intersection(*g_sets)
        inter_p = set.intersection(*p_sets)
        model_genes = [g for g in model_genes if (hb.canonicalize(g) or g) in inter_g]
        model_perts = sorted(inter_p)

    # Identify all required pools (target lines + mismatched pools referenced by target_lines)
    needed_pools: List[str] = list(target_lines)
    for ln in target_lines:
        st = resolve_r17_mismatched_strata(ln)
        for p in st["seen_sources"] + st["unseen_sources"] + st["same_lineage_sources"]:
            if p not in needed_pools:
                needed_pools.append(p)

    # Stage 1 & Stage 2: Extract fresh control pools and encode SE-600M if STATE_SE requested
    raw_pools: Dict[str, SampledRawControlPool] = {}
    ctrl_barcodes_json: Dict[str, Any] = {}
    x_state_map: Dict[str, Optional[np.ndarray]] = {}
    need_se = "STATE_SE" in models

    for pname in needed_pools:
        pool = extract_raw_control_pool_r17(
            line_name=pname,
            seed=ctrl_seed,
            n_ctrl=n_ctrl,
            cell_set_len=cell_set_len,
            min_ctrl_floor=MIN_CTRL_FLOOR_DEFAULT,
            custom_nadig_dir=raw_nadig_dir,
            custom_gse295214_dir=raw_gse295214_dir,
            pfizer_h5ad_path=pfizer_h5ad_path,
        )
        raw_pools[pname] = pool
        n_total = int(pool.sampling_audit.get("n_valid_pool", len(pool.batch_labels)))
        ctrl_barcodes_json[pname] = {
            "ctrl_seed": int(ctrl_seed),
            "ctrl_resample_seed": int(ctrl_resample_seed),
            "n_ctrl_total": int(n_total),
            "n_ctrl_sampled": int(pool.raw_counts.shape[0]),
            "barcodes": pool.cell_barcodes.tolist(),
            "barcode_prefix_convention": _describe_barcode_prefix_r17(pname),
            "cell_total_umis": [float(u) for u in pool.cell_total_umis.tolist()],
            "batch_labels": pool.batch_labels.tolist(),
            "sampling_mode": pool.sampling_audit.get("sampling_mode"),
        }
        if need_se:
            xs, se_audit = backend.encode_pool_se600m(pool, seed=ctrl_seed)
            x_state_map[pname] = xs
            safe_pool = pname.replace("/", "_")
            save_x_state_embedding_npz(
                emb_dir / f"x_state_{safe_pool}.npz",
                x_state=xs,
                control_barcodes=pool.cell_barcodes.tolist(),
                seed=ctrl_seed,
                dataset_name=pname,
                se_audit=se_audit,
            )
        else:
            x_state_map[pname] = None

    _dump_json(out_dir / "control_barcodes_r17.json", ctrl_barcodes_json)
    if need_se:
        _dump_json(emb_dir / "embeddings_manifest_sha256.json", build_manifest_sha256(emb_dir))

    # Stage 3: Export observed.npz, preds_<M>_sweep56.npz, preds_<M>_g4_mismatched.npz, SENS, manifests
    line_run_meta: Dict[str, Any] = {}
    for ln in target_lines:
        ld = arrays_dir / ln
        ld.mkdir(parents=True, exist_ok=True)
        pool = raw_pools[ln]
        ref = load_line_eval_reference_r17(
            target_line=ln,
            model_genes=model_genes,
            model_pert_names=model_perts,
            raw_ctrl_pool=pool,
            hgnc_bridge=hb,
            fct_pool_mode="primary",
            split_seed=seed,
        )
        obs_npz = build_observed_npz_dict(ref, pool)
        write_npz_strict_r17(ld / "observed.npz", obs_npz)

        if "sens_full_r9_6" in ref:
            s_ref = ref["sens_full_r9_6"]
            write_npz_strict_r17(
                ld / "SENS_FCT_full_r9_6_incl_vcc2025_h1.npz",
                {
                    "fct_delta_pred": np.asarray(s_ref["fct_delta_pred"], dtype=np.float64),
                    "delta_transfer": np.asarray(s_ref["delta_transfer"], dtype=np.float64),
                    "source_grand_mean": np.asarray(s_ref["source_grand_mean"], dtype=np.float64),
                    "target_names": _as_strict_str_array(s_ref["target_names"]),
                    "gene_names": _as_strict_str_array(s_ref["gene_names"]),
                },
            )

        per_m_meta: Dict[str, Any] = {}
        for mkey in models:
            per_m_meta[mkey] = export_line_model_arrays_r17(
                line_dir=ld,
                model_key=mkey,
                backend=backend,
                ref_bundle=ref,
                obs_npz=obs_npz,
                raw_pool_matched=pool,
                x_state_matched=x_state_map[ln],
                mismatched_pools=raw_pools,
                x_state_mismatched=x_state_map,
                cell_set_len=cell_set_len,
                seed=seed,
            )
        _dump_json(ld / "array_manifest_sha256.json", build_manifest_sha256(ld))
        line_run_meta[ln] = per_m_meta

    # Stage 4: Exact zero-leakage computation directly from exported arrays
    summary = recompute_cohort_from_exported_arrays_r17(
        arrays_dir=arrays_dir,
        models=models,
        n_boot=n_boot,
        n_boot_bonf56=n_boot_bonf56,
        n_perm=n_perm,
        n_splits_diag=n_splits_diag,
        seed=seed,
    )
    summary["run_metadata"] = {
        "phase": "real_r17",
        "device": resolve_compute_device(device),
        "ctrl_seed": int(ctrl_seed),
        "ctrl_resample_seed": int(ctrl_resample_seed),
        "n_ctrl": int(n_ctrl),
        "cell_set_len": int(cell_set_len),
        "gate_receipts": gate_receipts,
        "line_model_diagnostics": line_run_meta,
    }
    _dump_json(out_dir / "summary_r17.json", summary)
    return summary


def build_synthetic_r17_cohort(
    out_dir: Path,
    models: Sequence[str] = ("STATE_SE", "STATE_HVG"),
    target_lines: Sequence[str] = R17_ALL_EVAL_LINES,
    g_cnt: int = 120,
    n_ctrl: int = 64,
    cell_set_len: int = CELL_SET_LEN_DEFAULT,
    n_boot: int = 400,
    n_boot_bonf56: int = 1000,
    n_perm: int = 500,
    n_splits_diag: int = 10,
    seed: int = DEFAULT_R17_SEED,
    ctrl_seed: int = DEFAULT_R17_SEED,
    ctrl_resample_seed: int = DEFAULT_R17_RESAMPLE_SEED,
) -> Dict[str, Any]:
    """
    Build a complete pure-CPU synthetic Round 17 `R17_ARRAY_SCHEMA.md` cohort spanning ALL
    11 evaluation lines (including `N_eval >= 52` confirmatory lines, `gse295214_h1` with
    `SENS_FCT_full_r9_6_incl_vcc2025_h1.npz`, `mac_ctrl` `N_eval=3` edge case, `a549_ko_cas9pos`)
    plus all 6 mismatched control pools, export `.npz` strictly with `allow_pickle=False`,
    compute `summary_r17.json`, AND run full Stage 4 verification (`<= 1e-12`).
    """
    out_dir = Path(out_dir)
    emb_dir = out_dir / "embeddings"
    arrays_dir = out_dir / "arrays"
    out_dir.mkdir(parents=True, exist_ok=True)
    emb_dir.mkdir(parents=True, exist_ok=True)
    arrays_dir.mkdir(parents=True, exist_ok=True)

    eval_genes = [f"GENE_{i:04d}" for i in range(int(g_cnt))]
    max_p = 64
    master_perts = sorted(f"PERT_{i:03d}" for i in range(max_p))

    backend = SyntheticStateSweepBackend(model_genes=eval_genes, eval_perts=master_perts, seed=seed)

    # Synthetic N_eval sizes mirroring real PREREG_J structure cleanly
    syn_n_eval = {
        "cd4t_stim48": 54,
        "gse295214_h1": 62,
        "pfizer_telohaec_basal": 56,
        "pfizer_telohaec_cytokine": 56,
        "a549_myllia": 20,
        "hek293_xucao": 16,
        "telohaec": 9,
        "calu3": 8,
        "mac_lps": 4,
        "mac_ctrl": 3,
        "a549_ko_cas9pos": 56,
    }

    rng = np.random.default_rng(int(seed) + 999)
    n_sets = int(n_ctrl) // int(cell_set_len)
    all_pools: List[str] = list(target_lines) + [p for p in R16_MISMATCHED_POOLS if p not in target_lines]

    raw_pools: Dict[str, SampledRawControlPool] = {}
    x_state_map: Dict[str, Optional[np.ndarray]] = {}
    ctrl_barcodes_json: Dict[str, Any] = {}
    need_se = "STATE_SE" in models

    for pi, pname in enumerate(all_pools):
        lam = 1.5 + 0.5 * (pi % 5)
        raw_c = rng.poisson(lam, size=(n_ctrl, g_cnt)).astype(np.float32)
        b_lbls = np.asarray([f"{pname}_b{i % 3}" for i in range(n_ctrl)], dtype=str)
        bcs = np.asarray([f"{pname}_CELL_{i:04d}" for i in range(n_ctrl)], dtype=str)
        umis = np.sum(raw_c, axis=1, dtype=np.float64)
        pool = SampledRawControlPool(
            line_name=pname,
            raw_counts=raw_c,
            gene_names=eval_genes,
            gene_ids=eval_genes,
            batch_labels=b_lbls,
            cell_barcodes=bcs,
            cell_total_umis=umis,
            median_ntc_umi=float(np.median(umis)),
            sampling_audit={
                "n_valid_pool": int(n_ctrl * 4),
                "n_sampled": int(n_ctrl),
                "n_sets": int(n_sets),
                "sampling_mode": "synthetic_stratified_hamilton",
            },
            is_log1p=False,
            median_ntc_umi_all_valid=float(np.median(umis)),
        )
        raw_pools[pname] = pool
        ctrl_barcodes_json[pname] = {
            "ctrl_seed": int(ctrl_seed),
            "ctrl_resample_seed": int(ctrl_resample_seed),
            "n_ctrl_total": int(n_ctrl * 4),
            "n_ctrl_sampled": int(n_ctrl),
            "barcodes": bcs.tolist(),
            "barcode_prefix_convention": _describe_barcode_prefix_r17(pname),
            "cell_total_umis": [float(u) for u in umis.tolist()],
            "batch_labels": b_lbls.tolist(),
            "sampling_mode": "synthetic_stratified_hamilton",
        }
        if need_se:
            xs, se_audit = backend.encode_pool_se600m(pool, seed=ctrl_seed)
            x_state_map[pname] = xs
            save_x_state_embedding_npz(
                emb_dir / f"x_state_{pname}.npz",
                x_state=xs,
                control_barcodes=bcs.tolist(),
                seed=ctrl_seed,
                dataset_name=pname,
                se_audit=se_audit,
            )
        else:
            x_state_map[pname] = None

    _dump_json(out_dir / "control_barcodes_r17.json", ctrl_barcodes_json)
    if need_se:
        _dump_json(emb_dir / "embeddings_manifest_sha256.json", build_manifest_sha256(emb_dir))

    line_run_meta: Dict[str, Any] = {}
    for ln in target_lines:
        ld = arrays_dir / ln
        ld.mkdir(parents=True, exist_ok=True)
        pool = raw_pools[ln]
        p_cnt = syn_n_eval[ln]
        eval_perts = master_perts[:p_cnt]

        # Build coherent synthetic ground truth where batch effect from STATE_SE index 17/35 creates signal
        d_se_true_b0 = backend.pert_effects["STATE_SE"].mean(axis=0)[:p_cnt]
        d_obs = 0.85 * d_se_true_b0 + rng.normal(0.0, 0.08, size=(p_cnt, g_cnt))
        ha = d_obs + rng.normal(0.0, 0.04, size=(p_cnt, g_cnt))
        hb = d_obs + rng.normal(0.0, 0.04, size=(p_cnt, g_cnt))
        hv = np.ones(p_cnt, dtype=bool)
        if p_cnt >= 16:
            hv[-1] = False  # exercise half_valid stratification cleanly
        p1_idx, p2_idx = stratified_split_halves_c1(hv, split_idx=0, base_seed=seed)

        tm = 0.25 * d_obs + rng.normal(0.0, 0.05, size=(p_cnt, g_cnt))
        fct_d = 0.65 * d_obs + rng.normal(0.0, 0.06, size=(p_cnt, g_cnt))
        e1_d = 0.55 * d_obs + rng.normal(0.0, 0.07, size=(p_cnt, g_cnt))
        tg = np.asarray([(i * 2) % g_cnt for i in range(p_cnt)], dtype=np.int64)
        mu_tgt = np.mean(
            preprocess_hvg_expression(pool.raw_counts, "official_median_log1p", pool.cell_total_umis, pool.median_ntc_umi),
            axis=0,
        ).astype(np.float64)
        n_src = 5 if ln == "gse295214_h1" else 6
        mu_srcs = np.tile(mu_tgt[None, :], (n_src, 1)) + rng.normal(0.0, 0.05, size=(n_src, g_cnt))
        strata = resolve_r17_mismatched_strata(ln)

        ref = {
            "line_name": ln,
            "role": R17_LINE_ROLES[ln],
            "modality": R17_LINE_MODALITIES[ln],
            "checkpoint_id": "zeroshot/ensemble4",
            "eval_genes": eval_genes,
            "eval_perts": eval_perts,
            "target_indices": tg,
            "delta_obs": d_obs,
            "delta_half_a": ha,
            "delta_half_b": hb,
            "half_valid": hv,
            "delta_transfer": tm,
            "source_grand_mean": mu_srcs.mean(axis=0),
            "fct_delta_pred": fct_d,
            "e1_constmed_delta_pred": e1_d,
            "mu_target": mu_tgt,
            "mu_sources": mu_srcs,
            "obs_ntc_mean_official_median_log1p": mu_tgt,
            "split_s0_p1_indices": p1_idx,
            "split_s0_p2_indices": p2_idx,
            "seen_sources": list(strata["seen_sources"]),
            "unseen_sources": list(strata["unseen_sources"]),
            "same_lineage_sources": list(strata["same_lineage_sources"]),
        }
        obs_npz = build_observed_npz_dict(ref, pool)
        write_npz_strict_r17(ld / "observed.npz", obs_npz)

        if ln == "gse295214_h1":
            write_npz_strict_r17(
                ld / "SENS_FCT_full_r9_6_incl_vcc2025_h1.npz",
                {
                    "fct_delta_pred": fct_d + 0.01,
                    "delta_transfer": tm + 0.005,
                    "source_grand_mean": mu_srcs.mean(axis=0) + 0.01,
                    "target_names": _as_strict_str_array(eval_perts),
                    "gene_names": _as_strict_str_array(eval_genes),
                },
            )

        per_m_meta: Dict[str, Any] = {}
        for mkey in models:
            per_m_meta[mkey] = export_line_model_arrays_r17(
                line_dir=ld,
                model_key=mkey,
                backend=backend,
                ref_bundle=ref,
                obs_npz=obs_npz,
                raw_pool_matched=pool,
                x_state_matched=x_state_map[ln],
                mismatched_pools=raw_pools,
                x_state_mismatched=x_state_map,
                cell_set_len=cell_set_len,
                seed=seed,
            )
        _dump_json(ld / "array_manifest_sha256.json", build_manifest_sha256(ld))
        line_run_meta[ln] = per_m_meta

    summary = recompute_cohort_from_exported_arrays_r17(
        arrays_dir=arrays_dir,
        models=models,
        n_boot=n_boot,
        n_boot_bonf56=n_boot_bonf56,
        n_perm=n_perm,
        n_splits_diag=n_splits_diag,
        seed=seed,
    )
    summary["run_metadata"] = {
        "phase": "synthetic_e2e",
        "device": "cpu",
        "ctrl_seed": int(ctrl_seed),
        "ctrl_resample_seed": int(ctrl_resample_seed),
        "n_ctrl": int(n_ctrl),
        "cell_set_len": int(cell_set_len),
        "line_model_diagnostics": line_run_meta,
    }
    summary_path = out_dir / "summary_r17.json"
    _dump_json(summary_path, summary)

    # Mandatory Stage 4 round-trip audit check (`<= 1e-12`)
    with open(summary_path, "r", encoding="utf-8") as fh:
        saved_summary = json.load(fh)
    recomputed = recompute_cohort_from_exported_arrays_r17(
        arrays_dir=arrays_dir,
        models=models,
        n_boot=n_boot,
        n_boot_bonf56=n_boot_bonf56,
        n_perm=n_perm,
        n_splits_diag=n_splits_diag,
        seed=seed,
    )
    parity = compare_summaries_stage4(recomputed, saved_summary, tol=1e-12)
    summary["stage4_self_audit"] = parity
    _dump_json(summary_path, summary)
    return summary


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Clean-room Round 17 runner and Stage 4 recomputation CLI (`A-T3`)")
    p.add_argument("--phase", choices=("real_r17", "synthetic_e2e"), default=None)
    p.add_argument("--recompute-from-arrays", type=Path, default=None, help="Path to exported arrays/ directory")
    p.add_argument("--summary-json", type=Path, default=None, help="Optional summary_r17.json to verify parity")
    p.add_argument("--out-json", type=Path, default=None, help="Optional output path for Stage 4 recomputed JSON")
    p.add_argument("--cohort", default="all", choices=("all", "confirmatory"))
    p.add_argument("--lines", nargs="*", default=None)
    p.add_argument("--models", nargs="+", default=["STATE_SE", "STATE_HVG"], choices=("STATE_SE", "STATE_HVG"))
    p.add_argument("--output-dir", type=Path, default=AUDIT_R17_DIR / "results")
    p.add_argument("--models-root", type=Path, default=MODELS_ROOT_DEFAULT)
    p.add_argument("--status-dir", type=Path, default=R17_ROOT / "status")
    p.add_argument("--prereg-dir", type=Path, default=R17_ROOT / "prereg")
    p.add_argument("--raw-nadig-dir", type=Path, default=DEFAULT_NADIG_RAW_DIR)
    p.add_argument("--raw-gse295214-dir", type=Path, default=DEFAULT_GSE295214_RAW_DIR)
    p.add_argument("--pfizer-h5ad", type=Path, default=DEFAULT_PFIZER_H5AD)
    p.add_argument("--device", default="auto", choices=("auto", "cpu", "mps"))
    p.add_argument("--n-ctrl", type=int, default=N_CTRL_DEFAULT)
    p.add_argument("--cell-set-len", type=int, default=CELL_SET_LEN_DEFAULT)
    p.add_argument("--n-boot", type=int, default=4000)
    p.add_argument("--n-boot-bonf56", type=int, default=20000)
    p.add_argument("--n-perm", type=int, default=10000)
    p.add_argument("--n-splits-diag", type=int, default=20)
    p.add_argument("--seed", type=int, default=DEFAULT_R17_SEED)
    p.add_argument("--ctrl-seed", type=int, default=DEFAULT_R17_SEED)
    p.add_argument("--ctrl-resample-seed", type=int, default=DEFAULT_R17_RESAMPLE_SEED)
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = parse_args(argv)
    if args.recompute_from_arrays is not None:
        adir = Path(args.recompute_from_arrays)
        recomputed = recompute_cohort_from_exported_arrays_r17(
            arrays_dir=adir,
            models=args.models,
            n_boot=args.n_boot,
            n_boot_bonf56=args.n_boot_bonf56,
            n_perm=args.n_perm,
            n_splits_diag=args.n_splits_diag,
            seed=args.seed,
        )
        sum_path = args.summary_json if args.summary_json is not None else (adir.parent / "summary_r17.json")
        parity = None
        if sum_path is not None and sum_path.is_file():
            with open(sum_path, "r", encoding="utf-8") as fh:
                saved = json.load(fh)
            parity = compare_summaries_stage4(recomputed, saved, tol=1e-10)
            recomputed["stage4_parity_vs_saved"] = parity
        if args.out_json is not None:
            _dump_json(args.out_json, recomputed)
        msg = {"status": "RECOMPUTE_OK", "arrays_dir": str(adir), "parity": parity}
        print(json.dumps(msg, sort_keys=True))
        return

    lines = _resolve_target_lines(args.cohort, args.lines)
    if args.phase == "real_r17":
        summary = run_real_r17(
            target_lines=lines,
            models=args.models,
            output_dir=args.output_dir,
            models_root=args.models_root,
            status_dir=args.status_dir,
            prereg_dir=args.prereg_dir,
            raw_nadig_dir=args.raw_nadig_dir,
            raw_gse295214_dir=args.raw_gse295214_dir,
            pfizer_h5ad_path=args.pfizer_h5ad,
            device=args.device,
            n_ctrl=args.n_ctrl,
            cell_set_len=args.cell_set_len,
            ctrl_seed=args.ctrl_seed,
            ctrl_resample_seed=args.ctrl_resample_seed,
            n_boot=args.n_boot,
            n_boot_bonf56=args.n_boot_bonf56,
            n_perm=args.n_perm,
            n_splits_diag=args.n_splits_diag,
            seed=args.seed,
        )
        print(json.dumps({"status": "REAL_R17_DONE", "claim_rollups": summary["claim_rollups"]}, sort_keys=True))
        return

    if args.phase == "synthetic_e2e":
        # Guard CPU smoke vs real weights
        syn_n_ctrl = int(args.n_ctrl) if args.n_ctrl != N_CTRL_DEFAULT else 64
        summary = build_synthetic_r17_cohort(
            out_dir=args.output_dir,
            models=args.models,
            target_lines=lines,
            n_ctrl=syn_n_ctrl,
            cell_set_len=args.cell_set_len,
            n_boot=args.n_boot,
            n_boot_bonf56=args.n_boot_bonf56,
            n_perm=args.n_perm,
            n_splits_diag=args.n_splits_diag,
            seed=args.seed,
            ctrl_seed=args.ctrl_seed,
            ctrl_resample_seed=args.ctrl_resample_seed,
        )
        print(
            json.dumps(
                {
                    "status": "SYNTHETIC_E2E_OK",
                    "claim_rollups": summary["claim_rollups"],
                    "stage4_self_audit": summary["stage4_self_audit"],
                },
                sort_keys=True,
            )
        )
        return

    raise SystemExit("Must pass --phase {real_r17,synthetic_e2e} or --recompute-from-arrays <path>")


if __name__ == "__main__":
    main()
