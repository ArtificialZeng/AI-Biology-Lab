#!/usr/bin/env python3
"""
Round 16 Model Adapters, Zero-Leakage Guards, and Triple Hard Execution Gate (`PREREG_I v5 + Addendum A`).

Strict Safety & Licensing Guarantees (`I7`, `I8`, `I9`, `I10`, Gate 2 Verified, Addendum A Items 1–11):
  1. `assert_zero_training_leakage(...)`: Refuses any combination where `target_line`
     (or a homologous lineage such as `jurkat` / `nadig_jurkat` or `vcc2025_h1` / `gse295214_h1`)
     overlaps the checkpoint's training set or the public source summary catalog.
  2. `assert_all_execution_gates_r16(...)`: Triple hard gate refusing real-cell execution
     unless ALL THREE non-empty files exist in `round16/status/`:
       - `PREREG_I_COMMIT.txt`
       - `CODE_FREEZE_R16.txt`
       - `R16_GO.txt`
  3. Gate-2 Verified `ArcStateAdapter` Runtime Hooks:
       - `SE-600M`: loads `model.safetensors` (`715M` params, `0` missing keys) + `protein_embeddings.pt`
         (`19,790 x 5,120`, exact gene-symbol lookup per `PREREG_I §2.1.4`) via `StateEmbeddingModel`
         in pure `float32` (`mps`/`cpu`; `bf16`/`fp16` disallowed on MPS due to Metal assertion),
         using official `arc-state 0.11.1` `FilteredGenesCounts` + `VCIDatasetSentenceCollator`
         (`training=False`, `batch_size=16`, `shuffle=False`, `num_workers=0`) with per-pool RNG
         seeding immediately before encoding (`20261601` primary, `20261602` resample; Addendum A Item 11).
       - `ST-SE`: `padded=True` `16 x 32` sets -> trained `LatentToGeneDecoder` (`2058 -> 2000`) ->
         clip `[0.0, 14.0]` -> `mean(Y(p)) - mean(Y('non-targeting'))` on identical input controls.
       - `ST-HVG`: `padded=True` `16 x 32` sets -> direct `preds` (`2000`, NEVER untrained `gene_decoder`) ->
         clip `[0.0, 14.0]` -> `mean(Y(p)) - mean(Y('non-targeting'))`.
         Primary input normalization: official `arc-state` `normalize_total(median) -> log1p`
         (`state/_cli/_tx/_preprocess_train.py:45-55`) with `G_model` mapped via `HGNCBridge`
         (Addendum A Item 4); secondary sensitivity: `libnorm(1e4) -> log1p`.
       - Batch modes: `fixed_index_zero` (locked primary, `batch_idx=0`), `bypass_zero` (`0 in R^128`),
         `mean_batch_emb` (`mean_b W[b]`), `multi_batch_ens8` (`{0,7,14,21,28,35,42,49}`),
         `fixed_index:<k>`, and `native_gem_group` (for `G4` 2x2 decomposition on Group 1 and
         `nadig_hepg2` Arc archived `adata_pred.h5ad` check per Addendum A Item 3).
  4. Addendum A Item 5: `CellularAI_OpenCode_PublicRefit` is explicitly dropped pre-data (`不做`).
  5. Addendum A Item 8: `sample_disjoint_resample_control_sets_r16` draws a disjoint within-line
     second control set (`seed = 20261602`) without altering the primary `seed = 20261601` sample.
"""

from __future__ import annotations

import gc
import hashlib
import os
import pickle
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, FrozenSet, List, Optional, Sequence, Tuple

import numpy as np

from metrics_r16 import DEFAULT_R16_SEED, DEFAULT_RESAMPLE_SEED

ROUND16_DIR = Path(__file__).resolve().parent.parent
DEFAULT_STATUS_DIR = ROUND16_DIR / "status"
DEFAULT_MODELS_DIR = ROUND16_DIR / "shared" / "models"

STATE_CHECKPOINT_FOLDS: Tuple[str, ...] = ("hepg2", "jurkat", "k562", "rpe1")
MULTI_BATCH_ENS8_INDICES: Tuple[int, ...] = (0, 7, 14, 21, 28, 35, 42, 49)
PRED_CLIP_MAX: float = 14.0

LINE_CANONICAL_FAMILY: Dict[str, str] = {
    "hepg2": "hepg2",
    "nadig_hepg2": "hepg2",
    "jurkat": "jurkat",
    "nadig_jurkat": "jurkat",
    "k562": "k562",
    "replogle_k562_essential": "k562",
    "replogle_k562_gwps": "k562",
    "rpe1": "rpe1",
    "replogle_rpe1": "rpe1",
    "vcc2025_h1": "h1_hesc",
    "gse295214_h1": "h1_hesc",
    "nourreddine_kolf21j": "kolf21j",
    "a549_ko": "a549_gse345058",
    "a549_ko_cas9pos": "a549_gse345058",
    "a549_myllia": "a549_gse337804",
    "dld1": "dld1",
    "cd4t_rest": "cd4t_gse314342",
    "cd4t_stim48": "cd4t_gse314342",
    "feng_ipsc": "feng_ipsc",
    "hek293_xucao": "hek293",
    "mac_ctrl": "macrophage_gse327644",
    "mac_lps": "macrophage_gse327644",
    "cal33": "cal33",
    "hela": "hela",
    "ac16": "ac16",
    "molm13": "molm13",
    "u87mg": "u87mg",
}

STATE_REPLOGLE_TRAIN_FAMILIES: Dict[str, FrozenSet[str]] = {
    "zeroshot/hepg2": frozenset({"jurkat", "k562", "rpe1"}),
    "zeroshot/jurkat": frozenset({"hepg2", "k562", "rpe1"}),
    "zeroshot/k562": frozenset({"hepg2", "jurkat", "rpe1"}),
    "zeroshot/rpe1": frozenset({"hepg2", "jurkat", "k562"}),
    "zeroshot/ensemble4": frozenset({"hepg2", "jurkat", "k562", "rpe1"}),
}

REPLOGLE4_LOO_PRIMARY_LINES: Tuple[str, ...] = (
    "nadig_hepg2",
    "nadig_jurkat",
)
# A-T5 Q4: `feng_ipsc` (ENA 2024-11-08 < SE-600M 2025-06) demoted to SENSITIVITY_ONLY as a target
# line, while RETAINED as an `unseen_by_checkpoint` mismatched-control source.
GROUP2_UNSEEN_PRIMARY_LINES: Tuple[str, ...] = (
    "a549_ko",
    "dld1",
    "cd4t_rest",
)
# Group 3 prospective holdout (`PREREG_I v5 §2.4`): pre-download locked; reported strictly
# separately from Groups 1 & 2, with `vcc2025_h1` excluded from primary `FCT` source pool.
GROUP3_PROSPECTIVE_LINES: Tuple[str, ...] = (
    "gse295214_h1",
)
GROUP3_PRIMARY_FCT_SOURCES: Tuple[str, ...] = (
    "nadig_hepg2",
    "nadig_jurkat",
    "nourreddine_kolf21j",
    "replogle_k562_essential",
    "replogle_rpe1",
)
RETAINED_UNSEEN_MISMATCHED_ONLY_LINES: Tuple[str, ...] = (
    "feng_ipsc",
)
PRIMARY_LOCAL_RAW_LINES: Tuple[str, ...] = (
    REPLOGLE4_LOO_PRIMARY_LINES
    + GROUP2_UNSEEN_PRIMARY_LINES
    + RETAINED_UNSEEN_MISMATCHED_ONLY_LINES
)
OPTIONAL_DOWNLOAD_EXTENSION_LINES: Tuple[str, ...] = (
    "replogle_k562_essential",
    "replogle_rpe1",
    "vcc2025_h1",
    "nourreddine_kolf21j",
)
EXCLUDED_NO_LOCAL_RAW_FAMILIES: FrozenSet[str] = frozenset(
    LINE_CANONICAL_FAMILY[x] for x in OPTIONAL_DOWNLOAD_EXTENSION_LINES
)
SENSITIVITY_ONLY_LINES: Tuple[str, ...] = (
    "feng_ipsc",
    "cd4t_stim48",
    "a549_ko_cas9pos",
    "hek293_xucao",
    "a549_myllia",
)


def build_mismatched_strata_r16(
    target_line: str,
    checkpoint_id: str,
    group2_qualifying_lines: Sequence[str] = GROUP2_UNSEEN_PRIMARY_LINES,
    extra_unseen_mismatched_sources: Sequence[str] = RETAINED_UNSEEN_MISMATCHED_ONLY_LINES,
    seen_local_pool: Sequence[str] = REPLOGLE4_LOO_PRIMARY_LINES,
) -> Dict[str, Tuple[str, ...]]:
    """
    Build the exact `PREREG_I v5` (`I14`, `Q4`) two-stratum local mismatched-control pool:
      - `seen_by_checkpoint`: local raw-count lines in `seen_local_pool` present in the
        checkpoint's training set (excluding `target_line`'s canonical family).
      - `unseen_by_checkpoint`: qualifying local raw-count lines in
        `group2_qualifying_lines | extra_unseen_mismatched_sources` (retaining `feng_ipsc`
        as an unseen mismatched control per `Q4`) unseen by the checkpoint (excluding
        `target_line`'s canonical family).

    Strictly rejects any line lacking local raw single-cell counts (`nourreddine_kolf21j`,
    `replogle_k562_essential`, `replogle_rpe1`, `vcc2025_h1`) from the primary mismatched pool.
    """
    assert_zero_training_leakage("STATE_SE", target_line, checkpoint_id=checkpoint_id)
    train_fams = STATE_REPLOGLE_TRAIN_FAMILIES[checkpoint_id]
    t_fam = canonical_family(target_line)

    unseen_candidates = tuple(
        dict.fromkeys(tuple(group2_qualifying_lines) + tuple(extra_unseen_mismatched_sources))
    )
    for src in tuple(seen_local_pool) + unseen_candidates:
        if canonical_family(src) in EXCLUDED_NO_LOCAL_RAW_FAMILIES:
            raise ValueError(
                f"Primary mismatched pool cannot include non-local raw-count line {src!r} "
                "(requires separate human download approval per I14)."
            )

    seen_out = []
    for src in seen_local_pool:
        sf = canonical_family(src)
        if sf == t_fam:
            continue
        if sf not in train_fams:
            raise ValueError(
                f"Line {src!r} (family={sf!r}) in seen_local_pool is not in checkpoint "
                f"{checkpoint_id!r} training families {sorted(train_fams)}."
            )
        seen_out.append(src)

    unseen_out = []
    for src in unseen_candidates:
        sf = canonical_family(src)
        if sf == t_fam:
            continue
        if sf in train_fams:
            raise ValueError(
                f"Line {src!r} (family={sf!r}) in unseen mismatched pool overlaps checkpoint "
                f"{checkpoint_id!r} training families {sorted(train_fams)}."
            )
        unseen_out.append(src)

    return {
        "seen_by_checkpoint": tuple(sorted(seen_out)),
        "unseen_by_checkpoint": tuple(sorted(unseen_out)),
    }


class BlindingGateError(RuntimeError):
    """Raised when real-cell inference or archived `adata_pred.h5ad` access is attempted pre-gate."""


class TrainingLeakageError(ValueError):
    """Raised when a target cell line overlaps a checkpoint training set or public source pool."""


def canonical_family(line_name: str) -> str:
    key = line_name.strip().lower().replace("-", "").replace(" ", "")
    return LINE_CANONICAL_FAMILY.get(key, key)


def assert_prereg_i_locked(status_dir: Path = DEFAULT_STATUS_DIR) -> Path:
    """Gate required before opening Arc's archived `adata_pred.h5ad` or running Gate 5 (`I9`)."""
    lock_path = Path(status_dir) / "PREREG_I_COMMIT.txt"
    if not lock_path.is_file() or lock_path.stat().st_size == 0:
        raise BlindingGateError(
            f"Refusing execution/access before PREREG_I lock ({lock_path})."
        )
    return lock_path


def assert_code_freeze_r16(status_dir: Path = DEFAULT_STATUS_DIR) -> Path:
    """Hard gate: refuse any real-cell inference unless `CODE_FREEZE_R16.txt` exists and is non-empty."""
    freeze_path = Path(status_dir) / "CODE_FREEZE_R16.txt"
    if not freeze_path.is_file() or freeze_path.stat().st_size == 0:
        raise BlindingGateError(
            f"Refusing real-cell execution: missing or empty freeze gate {freeze_path}. "
            "Master Controller must lock PREREG_I and write CODE_FREEZE_R16.txt first."
        )
    return freeze_path


def assert_r16_go(status_dir: Path = DEFAULT_STATUS_DIR) -> Path:
    """Hard gate: refuse any real-cell inference unless `R16_GO.txt` exists and is non-empty."""
    go_path = Path(status_dir) / "R16_GO.txt"
    if not go_path.is_file() or go_path.stat().st_size == 0:
        raise BlindingGateError(
            f"Refusing real-cell execution: missing or empty release gate {go_path}. "
            "Master Controller must authorize Gate 5 by writing R16_GO.txt."
        )
    return go_path


def assert_all_execution_gates_r16(status_dir: Path = DEFAULT_STATUS_DIR) -> Dict[str, Path]:
    """
    Triple hard execution gate (`PREREG_I v5 §5` & `2026-11-25 23:52` Ruling 7):
    ALL THREE non-empty files (`PREREG_I_COMMIT.txt`, `CODE_FREEZE_R16.txt`, `R16_GO.txt`)
    must exist in `status_dir` before any real-cell data read or inference.
    """
    return {
        "prereg_lock": assert_prereg_i_locked(status_dir),
        "code_freeze": assert_code_freeze_r16(status_dir),
        "r16_go": assert_r16_go(status_dir),
    }


def assert_zero_training_leakage(
    model_key: str,
    target_line: str,
    checkpoint_id: Optional[str] = None,
    source_pool_lines: Optional[Sequence[str]] = None,
) -> None:
    """Verify zero training-set leakage for `STATE_SE`, `STATE_HVG`, or public source refits."""
    t_fam = canonical_family(target_line)
    if model_key in ("STATE_SE", "STATE_HVG"):
        if checkpoint_id is None or checkpoint_id not in STATE_REPLOGLE_TRAIN_FAMILIES:
            raise TrainingLeakageError(
                f"Unknown or missing STATE checkpoint_id={checkpoint_id!r} for model={model_key}"
            )
        train_fams = STATE_REPLOGLE_TRAIN_FAMILIES[checkpoint_id]
        if t_fam in train_fams:
            raise TrainingLeakageError(
                f"Training leakage blocked: target_line={target_line!r} (family={t_fam!r}) "
                f"is inside {model_key} checkpoint {checkpoint_id!r} training set {sorted(train_fams)}."
            )

    if source_pool_lines is not None:
        src_fams = {canonical_family(s) for s in source_pool_lines}
        if t_fam in src_fams:
            raise TrainingLeakageError(
                f"Source-pool leakage blocked: target_line={target_line!r} (family={t_fam!r}) "
                f"appears inside source_pool_lines={list(source_pool_lines)}."
            )


def verify_file_sha256(path: Path, expected_sha256: str) -> str:
    """Compute and verify SHA-256 digest of a downloaded checkpoint or catalog file."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    digest = h.hexdigest()
    if digest.lower() != expected_sha256.strip().lower():
        raise ValueError(f"SHA-256 mismatch for {path}: got {digest}, expected {expected_sha256}")
    return digest


class InsufficientControlCellsError(ValueError):
    """Raised when a control source has fewer than `min_ctrl` (`128`) QC-passing NTC cells (`M5.1`)."""


def sample_control_indices(
    n_total: int,
    n_sample: int,
    seed: int = DEFAULT_R16_SEED,
) -> np.ndarray:
    """Deterministically select `n_sample` control cell indices WITHOUT replacement (`M5.1`)."""
    if n_total <= 0 or n_sample <= 0:
        raise ValueError(f"Invalid control counts: n_total={n_total}, n_sample={n_sample}")
    if n_sample > n_total:
        raise InsufficientControlCellsError(
            f"Cannot sample {n_sample} control cells without replacement from n_total={n_total}."
        )
    rng = np.random.default_rng(seed)
    idx = rng.choice(n_total, size=n_sample, replace=False)
    return np.sort(idx.astype(np.int64))


def sample_stratified_control_sets_r16(
    n_total: int,
    batch_labels: Optional[Sequence[object]] = None,
    target_n: int = 512,
    cell_set_len: int = 32,
    min_ctrl: int = 128,
    seed: int = DEFAULT_R16_SEED,
) -> np.ndarray:
    """
    Locked `PREREG_I v5 §2.2` (`M5.1`) deterministic NTC control-cell sampler & shared partitioner:
      - If `n_total < min_ctrl` (`< 128`): raises `InsufficientControlCellsError` (source excluded).
      - If `min_ctrl <= n_total < target_n`: uses all cells WITHOUT replacement, shuffled once
        with `seed`, truncated to `(n_total // cell_set_len) * cell_set_len`.
      - If `n_total >= target_n` (`>= 512`): samples `target_n` cells WITHOUT replacement,
        proportionally stratified by `batch_labels` (Hamilton / largest-remainder allocation with
        lexicographic tie-break on batch label), shuffled once with `seed` into a single shared
        `(n_sets, cell_set_len)` (`16 x 32`) partition reused identically across all perturbations.
    """
    if cell_set_len <= 0 or target_n <= 0 or target_n % cell_set_len != 0:
        raise ValueError(
            f"target_n={target_n} must be a positive integer multiple of cell_set_len={cell_set_len}."
        )
    if n_total < min_ctrl:
        raise InsufficientControlCellsError(
            f"Control source has n_total={n_total} < min_ctrl={min_ctrl} NTC cells; exclude per M5.1."
        )
    if batch_labels is not None and len(batch_labels) != n_total:
        raise ValueError(
            f"batch_labels length ({len(batch_labels)}) must equal n_total ({n_total})."
        )

    rng = np.random.default_rng(seed)
    if n_total < target_n:
        n_use = (n_total // cell_set_len) * cell_set_len
        chosen = rng.permutation(n_total)[:n_use].astype(np.int64)
        return chosen.reshape(-1, cell_set_len)

    if batch_labels is None:
        chosen = rng.choice(n_total, size=target_n, replace=False).astype(np.int64)
        rng.shuffle(chosen)
        return chosen.reshape(-1, cell_set_len)

    groups: Dict[str, np.ndarray] = {}
    labels_str = [str(b) for b in batch_labels]
    unique_keys = sorted(set(labels_str))
    for k in unique_keys:
        groups[k] = np.array([i for i, lbl in enumerate(labels_str) if lbl == k], dtype=np.int64)

    # Exact integer Hamilton / largest-remainder allocation (zero float64 rounding ties)
    quotas = {k: min(len(groups[k]), (target_n * len(groups[k])) // n_total) for k in unique_keys}
    rem = target_n - sum(quotas.values())
    if rem > 0:
        order = sorted(
            unique_keys,
            key=lambda k: (-(target_n * len(groups[k]) - quotas[k] * n_total), k),
        )
        while rem > 0:
            progressed = False
            for k in order:
                if rem == 0:
                    break
                if quotas[k] < len(groups[k]):
                    quotas[k] += 1
                    rem -= 1
                    progressed = True
            if not progressed:
                break

    parts = []
    for k in unique_keys:
        q = quotas[k]
        if q > 0:
            parts.append(rng.choice(groups[k], size=q, replace=False).astype(np.int64))
    chosen = np.concatenate(parts, axis=0)
    rng.shuffle(chosen)
    return chosen.reshape(-1, cell_set_len)


def sample_disjoint_resample_control_sets_r16(
    n_total: int,
    primary_sets: np.ndarray,
    batch_labels: Optional[Sequence[object]] = None,
    target_n: int = 512,
    cell_set_len: int = 32,
    min_ctrl: int = 128,
    seed: int = DEFAULT_RESAMPLE_SEED,
) -> Tuple[Optional[np.ndarray], Dict[str, Any]]:
    """
    `PREREG_I Addendum A Item 8` deterministic disjoint within-line resample control sampler
    (`seed = 20261602`), leaving the locked `PREREG_I §2.2` primary sample (`primary_sets`,
    `seed = 20261601`) completely unchanged:
      - Excludes `primary_sets` from `0..n_total-1`, preserving ascending original cell order.
      - If `n_remaining >= target_n` (i.e., `n_total >= 2 * target_n = 1024`): draws a second
        batch-stratified `target_n` (`512`, `16 x 32`) set via `sample_stratified_control_sets_r16`
        with `seed = 20261602` (`overlap_with_primary == 0`).
      - If `min_ctrl <= n_remaining < target_n` (`640 <= n_total < 1024`): uses all disjoint
        remaining cells without replacement, shuffled with `seed = 20261602`, truncated to
        `(n_remaining // cell_set_len) * cell_set_len` (`>= 128`, `overlap_with_primary == 0`).
      - If `n_remaining < min_ctrl` (`< 128` disjoint cells remain): returns `(None, audit)`
        with `eligible = False`.
    """
    if n_total <= 0:
        raise ValueError(f"n_total must be positive, got {n_total}.")
    if batch_labels is not None and len(batch_labels) != n_total:
        raise ValueError(
            f"batch_labels length ({len(batch_labels)}) must equal n_total ({n_total})."
        )
    primary_flat = np.asarray(primary_sets, dtype=np.int64).ravel()
    if primary_flat.size == 0 or np.any(primary_flat < 0) or np.any(primary_flat >= n_total):
        raise ValueError(f"primary_sets contains out-of-range indices for n_total={n_total}.")
    if len(set(primary_flat.tolist())) != int(primary_flat.size):
        raise ValueError("primary_sets must not contain duplicate cell indices.")

    in_primary = np.zeros(n_total, dtype=bool)
    in_primary[primary_flat] = True
    rem_orig_idx = np.flatnonzero(~in_primary).astype(np.int64)
    n_rem = int(rem_orig_idx.size)

    if n_rem < int(min_ctrl):
        return None, {
            "eligible": False,
            "reason": f"remaining_disjoint_ntc={n_rem} < min_ctrl={int(min_ctrl)}",
            "n_total": int(n_total),
            "n_primary": int(primary_flat.size),
            "n_remaining_disjoint": n_rem,
            "n_resample": 0,
            "n_sets": 0,
            "cell_set_len": int(cell_set_len),
            "seed": int(seed),
            "overlap_with_primary": 0,
        }

    rem_batches = [batch_labels[int(i)] for i in rem_orig_idx] if batch_labels is not None else None
    local_sets = sample_stratified_control_sets_r16(
        n_total=n_rem,
        batch_labels=rem_batches,
        target_n=target_n,
        cell_set_len=cell_set_len,
        min_ctrl=min_ctrl,
        seed=seed,
    )
    resample_sets = rem_orig_idx[local_sets].astype(np.int64)
    overlap = len(set(primary_flat.tolist()) & set(resample_sets.ravel().tolist()))
    if overlap != 0:
        raise RuntimeError(f"Internal error: disjoint resample overlapped primary by {overlap} cells.")

    regime = f"full_{int(target_n)}_disjoint" if n_rem >= int(target_n) else "truncated_disjoint_complement"
    return resample_sets, {
        "eligible": True,
        "regime": regime,
        "n_total": int(n_total),
        "n_primary": int(primary_flat.size),
        "n_remaining_disjoint": n_rem,
        "n_resample": int(resample_sets.size),
        "n_sets": int(resample_sets.shape[0]),
        "cell_set_len": int(resample_sets.shape[1]),
        "seed": int(seed),
        "overlap_with_primary": 0,
    }


def preprocess_hvg_counts_r16(
    raw_counts: np.ndarray,
    raw_gene_names: Sequence[str],
    g_model_names: Sequence[str],
    norm_mode: str = "official_median_log1p",
    canon_map: Optional[Dict[str, str]] = None,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    Preprocess full-transcriptome raw UMI counts into `ST-HVG`'s ordered `2,000`-gene `X_hvg` input:
      - `official_median_log1p` (locked PRIMARY per `PREREG_I v5 §2.1`, matching
        `state/_cli/_tx/_preprocess_train.py:45-55`): scales each cell by
        `median_i(total_umi_i > 0) / total_umi_i` (exact `sc.pp.normalize_total(target_sum=None)`)
        on the full measured transcriptome, applies `log1p`, then slices ordered `g_model_names`
        (zero-filling missing genes).
      - `libnorm1e4_log1p` (locked SECONDARY sensitivity): scales each cell by `1e4 / total_umi_i`
        on the full measured transcriptome, applies `log1p`, then slices ordered `g_model_names`.
    """
    X = np.asarray(raw_counts, dtype=np.float64)
    if X.ndim != 2 or X.shape[1] != len(raw_gene_names):
        raise ValueError(
            f"raw_counts shape {X.shape} incompatible with len(raw_gene_names)={len(raw_gene_names)}."
        )
    X = np.maximum(X, 0.0)
    totals = X.sum(axis=1)
    pos = totals > 0.0
    if norm_mode == "official_median_log1p":
        target_sum = float(np.median(totals[pos])) if np.any(pos) else 1.0
    elif norm_mode == "libnorm1e4_log1p":
        target_sum = 1e4
    else:
        raise ValueError(f"Unsupported HVG norm_mode: {norm_mode!r}")

    scale = np.where(pos, target_sum / np.maximum(totals, 1e-12), 0.0)

    # Map canonicalized gene symbol -> column index in raw_counts (keep highest mean count on dups)
    raw_mean = X.mean(axis=0)
    sym_to_col: Dict[str, Tuple[float, int]] = {}
    for j, g in enumerate(raw_gene_names):
        cg = canon_map.get(g, g) if canon_map is not None else g
        prev = sym_to_col.get(cg)
        if prev is None or raw_mean[j] > prev[0]:
            sym_to_col[cg] = (float(raw_mean[j]), j)

    n_cells = X.shape[0]
    n_hvg = len(g_model_names)
    X_hvg = np.zeros((n_cells, n_hvg), dtype=np.float32)
    n_mapped = 0
    for k, gm in enumerate(g_model_names):
        cgm = canon_map.get(gm, gm) if canon_map is not None else gm
        hit = sym_to_col.get(cgm)
        if hit is not None:
            col_j = hit[1]
            X_hvg[:, k] = np.log1p(X[:, col_j] * scale).astype(np.float32)
            n_mapped += 1

    audit = {
        "norm_mode": norm_mode,
        "target_sum": float(target_sum),
        "n_cells": int(n_cells),
        "n_raw_genes": int(len(raw_gene_names)),
        "n_g_model": int(n_hvg),
        "n_g_model_mapped": int(n_mapped),
    }
    return X_hvg, audit


def compute_embedding_separation_diagnostics(
    emb_a: np.ndarray,
    emb_b: np.ndarray,
) -> Dict[str, float]:
    """
    Compute `A-T3 §2.3` `SE-600M` embedding non-degeneracy diagnostics between two control pools:
    centroid cosine similarity and Euclidean distance in both `2,048`-dim cell space and full `2,058`-dim space.
    """
    A = np.asarray(emb_a, dtype=np.float64)
    B = np.asarray(emb_b, dtype=np.float64)
    if A.ndim != 2 or B.ndim != 2 or A.shape[1] != B.shape[1]:
        raise ValueError(f"Embedding shape mismatch: A={A.shape}, B={B.shape}")
    mu_a = A.mean(axis=0)
    mu_b = B.mean(axis=0)
    d_full = float(np.linalg.norm(mu_a - mu_b))
    n_a = float(np.linalg.norm(mu_a))
    n_b = float(np.linalg.norm(mu_b))
    cos_full = float(np.dot(mu_a, mu_b) / (n_a * n_b)) if (n_a > 0 and n_b > 0) else float("nan")

    d_cell = min(2048, A.shape[1])
    mu_a_c = mu_a[:d_cell]
    mu_b_c = mu_b[:d_cell]
    d_2048 = float(np.linalg.norm(mu_a_c - mu_b_c))
    na_c = float(np.linalg.norm(mu_a_c))
    nb_c = float(np.linalg.norm(mu_b_c))
    cos_2048 = float(np.dot(mu_a_c, mu_b_c) / (na_c * nb_c)) if (na_c > 0 and nb_c > 0) else float("nan")

    return {
        "centroid_cos_2048": cos_2048,
        "centroid_l2_2048": d_2048,
        "centroid_cos_full": cos_full,
        "centroid_l2_full": d_full,
    }


def compute_pseudobulk_log1p_shift(
    pred_cells: np.ndarray,
    ctrl_cells: np.ndarray,
    input_is_log1p: bool = False,
    libnorm_scale: float = 1e4,
) -> np.ndarray:
    """
    Compute pseudobulk `log1p` perturbation shift `mean(log1p(X_pred)) - mean(log1p(X_ctrl))`.
    When `input_is_log1p=False`, applies standard `10^4` library-size normalization first.
    """
    xp = np.asarray(pred_cells, dtype=np.float64)
    xc = np.asarray(ctrl_cells, dtype=np.float64)
    if xp.ndim != 2 or xc.ndim != 2 or xp.shape[1] != xc.shape[1]:
        raise ValueError(f"Shape mismatch: pred={xp.shape}, ctrl={xc.shape}")
    if input_is_log1p:
        return np.mean(xp, axis=0) - np.mean(xc, axis=0)

    xp_pos = np.maximum(xp, 0.0)
    xc_pos = np.maximum(xc, 0.0)
    sp = np.maximum(xp_pos.sum(axis=1, keepdims=True), 1e-12)
    sc = np.maximum(xc_pos.sum(axis=1, keepdims=True), 1e-12)
    lp = np.log1p(xp_pos * (libnorm_scale / sp))
    lc = np.log1p(xc_pos * (libnorm_scale / sc))
    return np.mean(lp, axis=0) - np.mean(lc, axis=0)


def encode_controls_se600m_r16(
    raw_counts: np.ndarray,
    raw_gene_names: Sequence[str],
    se_dir: Path = DEFAULT_MODELS_DIR / "SE-600M",
    device: str = "mps",
    batch_size: int = 16,
    dataset_name: str = "r16_ctrl",
    seed: int = DEFAULT_R16_SEED,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    Encode raw integer UMI control cells into `2,058`-dim `X_state` embeddings using frozen
    `SE-600M/model.safetensors` + `protein_embeddings.pt` + `config.yaml` in pure `float32`
    (`PREREG_I Addendum A Item 11`):
      - Official `arc-state 0.11.1` `FilteredGenesCounts` and `VCIDatasetSentenceCollator`
        (`is_train=False`, `collator.training=False`).
      - `batch_size=16`, `shuffle=False`, `num_workers=0`, preserving sampled cell order and
        raw file gene-column order with exact gene-symbol lookup into `protein_embeddings.pt`.
      - Calls `random.seed(seed)`, `np.random.seed(seed)`, `torch.manual_seed(seed)`, and
        (on MPS) `torch.mps.manual_seed(seed)` ONCE per pool IMMEDIATELY before `for batch in dl:`
        (strictly after `StateEmbeddingModel` construction & weight loading so model init never
        perturbs global torch RNG prior to `loader.py` `torch.randperm`/`torch.randint`).
      - Strictly avoids `MPS bfloat16`/`float16`, which abort inside macOS Metal.
    """
    raw_arr = np.asarray(raw_counts, dtype=np.float32)
    if raw_arr.ndim != 2 or raw_arr.shape[1] != len(raw_gene_names):
        raise ValueError(
            f"raw_counts shape {raw_arr.shape} incompatible with len(raw_gene_names)={len(raw_gene_names)}."
        )
    if raw_arr.shape[0] <= 0:
        raise ValueError(f"raw_counts must contain at least 1 cell, got shape {raw_arr.shape}.")
    if not np.all(np.isfinite(raw_arr)):
        raise ValueError("raw_counts contains non-finite values (NaN/Inf).")
    if isinstance(batch_size, bool) or int(batch_size) <= 0:
        raise ValueError(f"batch_size must be a positive integer, got {batch_size!r}.")
    if isinstance(seed, bool) or seed is None:
        raise ValueError(f"seed must be an integer, got {seed!r}.")
    seed_int = int(seed)
    batch_size_int = int(batch_size)

    os.environ.setdefault("WANDB_MODE", "disabled")
    os.environ.setdefault("HF_HUB_OFFLINE", "1")
    os.environ.setdefault("TRANSFORMERS_OFFLINE", "1")

    import anndata as ad
    import scipy.sparse as sp
    import torch
    import torch.nn as nn
    from omegaconf import OmegaConf
    from safetensors.torch import load_file as load_safetensors
    from torch.utils.data import DataLoader
    from state.emb.data.loader import FilteredGenesCounts, VCIDatasetSentenceCollator
    from state.emb.nn.model import StateEmbeddingModel
    from state.emb.utils import get_embedding_cfg

    se_path = Path(se_dir)
    cfg = OmegaConf.load(se_path / "config.yaml")
    pe = torch.load(se_path / "protein_embeddings.pt", map_location="cpu", weights_only=False)
    pe_set = set(pe.keys())
    mapped_raw_rows = sum(1 for g in raw_gene_names if g in pe_set)
    mapped_unique = len({g for g in raw_gene_names if g in pe_set})
    min_required_genes = int(getattr(cfg.dataset, "cls_token_idx", 3)) + 1
    if mapped_unique < min_required_genes:
        raise ValueError(
            f"SE-600M requires at least {min_required_genes} unique genes mapped to protein_embeddings.pt; "
            f"got {mapped_unique} out of {len(raw_gene_names)} raw genes."
        )

    dev_str = device if (device != "mps" or torch.backends.mps.is_available()) else "cpu"
    dev = torch.device(dev_str)

    se_model = StateEmbeddingModel(
        token_dim=get_embedding_cfg(cfg).size,
        d_model=cfg.model.emsize,
        nhead=cfg.model.nhead,
        d_hid=cfg.model.d_hid,
        nlayers=cfg.model.nlayers,
        output_dim=cfg.model.output_dim,
        dropout=0.0,
        compiled=False,
        emb_size=get_embedding_cfg(cfg).size,
        cfg=cfg,
    )
    sd = load_safetensors(str(se_path / "model.safetensors"), device="cpu")
    incomp = se_model.load_state_dict(sd, strict=False)
    del sd
    if len(incomp.missing_keys) > 0:
        raise RuntimeError(f"SE-600M missing keys on load: {incomp.missing_keys}")

    all_pe = torch.vstack(list(pe.values()))
    se_model.pe_embedding = nn.Embedding.from_pretrained(all_pe)
    se_model.binary_decoder.requires_grad = False
    se_model.to(device=dev, dtype=torch.float32)
    se_model.pe_embedding.to(device=dev, dtype=torch.float32)
    del all_pe
    se_model.eval()

    X_csr = sp.csr_matrix(raw_arr)
    adata = ad.AnnData(X=X_csr)
    adata.var["gene_name"] = list(raw_gene_names)
    adata.var_names = list(raw_gene_names)
    adata.var_names_make_unique()
    adata.obs_names = [f"{dataset_name}_{i}" for i in range(X_csr.shape[0])]

    ds = FilteredGenesCounts(
        cfg,
        adata=adata,
        adata_name=dataset_name,
        shape_dict={dataset_name: adata.shape},
        protein_embeds=pe,
        gene_column="gene_name",
    )
    collator = VCIDatasetSentenceCollator(
        cfg,
        valid_gene_mask=ds.valid_gene_index,
        ds_emb_mapping_inference=ds.ds_emb_map,
        is_train=False,
        precision=torch.float32,
    )
    collator.training = False
    dl = DataLoader(ds, batch_size=batch_size_int, shuffle=False, collate_fn=collator, num_workers=0)

    # Addendum A Item 11: seed ONCE per control pool immediately before `for batch in dl:`
    # (strictly after `StateEmbeddingModel` construction & weight loading).
    random.seed(seed_int)
    np.random.seed(seed_int)
    torch.manual_seed(seed_int)
    if dev_str == "mps" and hasattr(torch, "mps") and hasattr(torch.mps, "manual_seed"):
        torch.mps.manual_seed(seed_int)

    embs, ds_embs = [], []
    with torch.no_grad():
        for batch in dl:
            if isinstance(batch, tuple):
                bdev = tuple(x.to(dev) if isinstance(x, torch.Tensor) else x for x in batch)
            elif isinstance(batch, dict):
                bdev = {k: (v.to(dev) if isinstance(v, torch.Tensor) else v) for k, v in batch.items()}
            else:
                bdev = batch
            _, _, _, emb, ds_emb = se_model._compute_embedding_for_batch(bdev)
            ds_proj = se_model.dataset_embedder(ds_emb)
            embs.append(emb.detach().cpu().float().numpy())
            ds_embs.append(ds_proj.detach().cpu().float().numpy())
    if dev_str == "mps":
        torch.mps.synchronize()

    embs_np = np.concatenate(embs, axis=0).astype(np.float32)
    ds_np = np.concatenate(ds_embs, axis=0).astype(np.float32)
    x_state = np.ascontiguousarray(np.concatenate([embs_np, ds_np], axis=-1), dtype=np.float32)
    x_state_sha256 = hashlib.sha256(x_state.tobytes()).hexdigest()

    del se_model
    gc.collect()
    if dev_str == "mps":
        torch.mps.empty_cache()

    audit = {
        "device": dev_str,
        "dtype": "float32",
        "seed": seed_int,
        "batch_size": batch_size_int,
        "shuffle": False,
        "num_workers": 0,
        "collator_training": False,
        "n_cells": int(x_state.shape[0]),
        "embed_dim": int(x_state.shape[1]),
        "n_raw_genes": int(len(raw_gene_names)),
        "n_mapped_to_pe19790": int(mapped_unique),
        "n_mapped_to_pe19790_raw_rows": int(mapped_raw_rows),
        "pe_total_genes": int(len(pe)),
        "emb2048_norm_mean": float(np.linalg.norm(embs_np, axis=1).mean()),
        "ds10_norm_mean": float(np.linalg.norm(ds_np, axis=1).mean()),
        "all_finite": bool(np.isfinite(x_state).all()),
        "x_state_sha256": x_state_sha256,
    }
    return x_state, audit


@dataclass(frozen=True)
class ModelAuditSpec:
    model_key: str
    role: str
    license_id: str
    checkpoint_id: Optional[str]
    supports_gene_perm: bool
    n_ctrl_infer: int
    cell_set_len: Optional[int] = None


class BasePerturbModelAdapter:
    """Abstract zero-shot perturbation model adapter enforcing leakage and triple execution gates."""

    def __init__(self, spec: ModelAuditSpec, status_dir: Path = DEFAULT_STATUS_DIR):
        self.spec = spec
        self.status_dir = Path(status_dir)

    def predict_deltas(
        self,
        target_line: str,
        control_line_label: str,
        ctrl_cells: np.ndarray,
        pert_names: Sequence[str],
        gene_names: Sequence[str],
        source_pool_lines: Optional[Sequence[str]] = None,
        allow_synthetic_without_freeze: bool = False,
        **kwargs: Any,
    ) -> np.ndarray:
        assert_zero_training_leakage(
            model_key=self.spec.model_key,
            target_line=target_line,
            checkpoint_id=self.spec.checkpoint_id,
            source_pool_lines=source_pool_lines,
        )
        if not allow_synthetic_without_freeze:
            assert_all_execution_gates_r16(self.status_dir)
        return self._predict_impl(
            target_line=target_line,
            control_line_label=control_line_label,
            ctrl_cells=ctrl_cells,
            pert_names=pert_names,
            gene_names=gene_names,
            **kwargs,
        )

    def _predict_impl(
        self,
        target_line: str,
        control_line_label: str,
        ctrl_cells: np.ndarray,
        pert_names: Sequence[str],
        gene_names: Sequence[str],
        **kwargs: Any,
    ) -> np.ndarray:
        raise NotImplementedError


class ArcStateAdapter(BasePerturbModelAdapter):
    """
    Gate-2 Verified Adapter for Arc `STATE` (`STATE_SE` primary; `STATE_HVG` secondary).
    Enforces zero-leakage guard and triple hard gate (`PREREG_I_COMMIT.txt` +
    `CODE_FREEZE_R16.txt` + `R16_GO.txt`) unless `allow_synthetic_without_freeze=True`.
    """

    def __init__(
        self,
        variant: str,
        checkpoint_id: str,
        checkpoint_path: Optional[Path] = None,
        se_weights_path: Optional[Path] = None,
        models_dir: Path = DEFAULT_MODELS_DIR,
        device: str = "mps",
        batch_mode: str = "fixed_index_zero",
        hvg_norm_mode: str = "official_median_log1p",
        status_dir: Path = DEFAULT_STATUS_DIR,
    ):
        if variant not in ("STATE_SE", "STATE_HVG"):
            raise ValueError(f"Unsupported STATE variant: {variant}")
        spec = ModelAuditSpec(
            model_key=variant,
            role="primary" if variant == "STATE_SE" else "secondary",
            license_id="Arc-Research-Institute-State-Non-Commercial-2025-06-23",
            checkpoint_id=checkpoint_id,
            supports_gene_perm=(variant == "STATE_HVG"),
            n_ctrl_infer=512,
            cell_set_len=32,
        )
        super().__init__(spec=spec, status_dir=status_dir)
        self.variant = variant
        self.checkpoint_id = checkpoint_id
        self.models_dir = Path(models_dir)
        self.checkpoint_path = Path(checkpoint_path) if checkpoint_path is not None else None
        self.se_weights_path = Path(se_weights_path) if se_weights_path is not None else None
        self.device = device
        self.batch_mode = batch_mode
        self.hvg_norm_mode = hvg_norm_mode

    def _resolve_fold_dirs(self) -> List[Path]:
        if self.checkpoint_path is not None:
            p = self.checkpoint_path
            if p.is_file():
                # checkpoints/final.ckpt -> fold dir is parent.parent
                return [p.parent.parent if p.parent.name == "checkpoints" else p.parent]
            if not p.exists():
                raise FileNotFoundError(f"STATE checkpoint not present ({p}).")
            return [p]

        repo_name = "ST-SE-Replogle" if self.variant == "STATE_SE" else "ST-HVG-Replogle"
        if self.checkpoint_id == "zeroshot/ensemble4":
            dirs = [self.models_dir / repo_name / "zeroshot" / f for f in STATE_CHECKPOINT_FOLDS]
        elif self.checkpoint_id.startswith("zeroshot/"):
            fold = self.checkpoint_id.split("/", 1)[1]
            dirs = [self.models_dir / repo_name / "zeroshot" / fold]
        else:
            raise ValueError(f"Unsupported checkpoint_id: {self.checkpoint_id!r}")

        for d in dirs:
            ckpt = d / "checkpoints" / "final.ckpt"
            if not ckpt.is_file():
                raise FileNotFoundError(f"STATE checkpoint not present ({ckpt}).")
        return dirs

    @staticmethod
    def _infer_single_fold(
        fold_dir: Path,
        variant: str,
        ctrl_tensor_np: np.ndarray,
        pert_keys_model: Sequence[str],
        batch_mode: str,
        device_str: str,
        native_batch_labels: Optional[Sequence[str]] = None,
    ) -> Tuple[np.ndarray, List[str]]:
        """
        Run one frozen `ST-SE` or `ST-HVG` fold on `ctrl_tensor_np` (shape `(N_ctrl, D_in)`,
        exact multiple of `32`), returning `(delta_2000 [len(pert_keys_model), 2000], g_model_2000)`.
        Uses `padded=True` (`16 x 32` sets) and `mean(clip(Y(p),0,14)) - mean(clip(Y('non-targeting'),0,14))`.
        """
        os.environ.setdefault("WANDB_MODE", "disabled")
        os.environ.setdefault("HF_HUB_OFFLINE", "1")
        os.environ.setdefault("TRANSFORMERS_OFFLINE", "1")

        import torch
        from state.tx.models.state_transition import StateTransitionPerturbationModel

        n_cells = int(ctrl_tensor_np.shape[0])
        if n_cells <= 0 or n_cells % 32 != 0:
            raise ValueError(f"ctrl_tensor_np rows ({n_cells}) must be a positive multiple of 32.")

        with open(fold_dir / "var_dims.pkl", "rb") as f:
            vd = pickle.load(f)
        g_model = list(vd["gene_names"])
        pert_map = torch.load(fold_dir / "pert_onehot_map.pt", map_location="cpu", weights_only=False)
        batch_map = torch.load(fold_dir / "batch_onehot_map.torch", map_location="cpu", weights_only=False)
        b_str_to_idx = {str(k): int(torch.argmax(v).item()) for k, v in batch_map.items()}

        dev_use = device_str if (device_str != "mps" or torch.backends.mps.is_available()) else "cpu"
        dev = torch.device(dev_use)
        random.seed(DEFAULT_R16_SEED)
        np.random.seed(DEFAULT_R16_SEED)
        torch.manual_seed(DEFAULT_R16_SEED)
        if dev_use == "mps" and hasattr(torch, "mps") and hasattr(torch.mps, "manual_seed"):
            torch.mps.manual_seed(DEFAULT_R16_SEED)
        ckpt_file = fold_dir / "checkpoints" / "final.ckpt"
        model = StateTransitionPerturbationModel.load_from_checkpoint(str(ckpt_file), map_location="cpu")
        model.eval().to(dev)

        orig_w = model.batch_encoder.weight.detach().clone()
        ctrl_dev = torch.from_numpy(np.asarray(ctrl_tensor_np, dtype=np.float32)).to(dev)

        # Build list of (override_mode, batch_indices_tensor) runs to execute & average
        runs: List[Tuple[str, torch.Tensor]] = []
        if batch_mode == "fixed_index_zero":
            runs.append(("orig", torch.zeros(n_cells, dtype=torch.long, device=dev)))
        elif batch_mode == "bypass_zero":
            runs.append(("zero", torch.zeros(n_cells, dtype=torch.long, device=dev)))
        elif batch_mode == "mean_batch_emb":
            runs.append(("mean", torch.zeros(n_cells, dtype=torch.long, device=dev)))
        elif batch_mode == "multi_batch_ens8":
            for b_i in MULTI_BATCH_ENS8_INDICES:
                runs.append(("orig", torch.full((n_cells,), int(b_i), dtype=torch.long, device=dev)))
        elif batch_mode.startswith("fixed_index:"):
            b_i = int(batch_mode.split(":", 1)[1])
            runs.append(("orig", torch.full((n_cells,), b_i, dtype=torch.long, device=dev)))
        elif batch_mode == "native_gem_group":
            if native_batch_labels is None or len(native_batch_labels) != n_cells:
                raise ValueError("native_gem_group mode requires native_batch_labels of length n_cells.")
            idxs = []
            for lbl in native_batch_labels:
                s_lbl = str(lbl).strip()
                if s_lbl in b_str_to_idx:
                    idxs.append(b_str_to_idx[s_lbl])
                elif s_lbl.endswith(".0") and s_lbl[:-2] in b_str_to_idx:
                    idxs.append(b_str_to_idx[s_lbl[:-2]])
                else:
                    idxs.append(0)
            runs.append(("orig", torch.tensor(idxs, dtype=torch.long, device=dev)))
        else:
            raise ValueError(f"Unsupported batch_mode: {batch_mode!r}")

        all_query = ["non-targeting"] + list(pert_keys_model)
        pb_accum = np.zeros((len(all_query), len(g_model)), dtype=np.float64)

        with torch.no_grad():
            for override, b_dev in runs:
                model.batch_encoder.weight.copy_(orig_w)
                if override == "zero":
                    model.batch_encoder.weight[0].zero_()
                elif override == "mean":
                    model.batch_encoder.weight[0].copy_(orig_w.mean(dim=0))

                for qi, pk in enumerate(all_query):
                    if pk not in pert_map:
                        raise KeyError(f"Perturbation key {pk!r} missing from pert_onehot_map.pt")
                    p_dev = pert_map[pk].to(device=dev, dtype=torch.float32).unsqueeze(0).expand(n_cells, -1)
                    pred_raw = model.forward(
                        {"ctrl_cell_emb": ctrl_dev, "pert_emb": p_dev, "batch": b_dev},
                        padded=True,
                    )
                    if variant == "STATE_SE":
                        gene_pred = model.gene_decoder(pred_raw)
                    else:
                        # Gate-2 verified: ST-HVG gene_decoder is untrained init; use preds directly.
                        gene_pred = pred_raw
                    g_np = np.clip(gene_pred.detach().cpu().numpy().astype(np.float64), 0.0, PRED_CLIP_MAX)
                    pb_accum[qi] += g_np.mean(axis=0)
            model.batch_encoder.weight.copy_(orig_w)
            if dev_use == "mps":
                torch.mps.synchronize()

        pb_mean = pb_accum / float(len(runs))
        ntc_mean = pb_mean[0:1, :]
        delta_2000 = pb_mean[1:, :] - ntc_mean

        del model
        if dev_use == "mps":
            torch.mps.empty_cache()
        return delta_2000, g_model

    def _predict_impl(
        self,
        target_line: str,
        control_line_label: str,
        ctrl_cells: np.ndarray,
        pert_names: Sequence[str],
        gene_names: Sequence[str],
        ctrl_is_embedded_x_state: bool = False,
        ctrl_raw_gene_names: Optional[Sequence[str]] = None,
        pert_canon_to_model: Optional[Dict[str, str]] = None,
        gene_canon_map: Optional[Dict[str, str]] = None,
        native_batch_labels: Optional[Sequence[str]] = None,
        return_per_fold: bool = False,
        se_encode_seed: int = DEFAULT_R16_SEED,
    ) -> Any:
        """
        Execute frozen `STATE_SE` or `STATE_HVG` inference aligned strictly to `(pert_names, gene_names)`.
          - For `STATE_SE`: pass pre-cached `2,058`-dim `X_state` (`ctrl_is_embedded_x_state=True`),
            or raw integer UMI counts + `ctrl_raw_gene_names` (encoded on the fly via `SE-600M`
            with per-pool seed `se_encode_seed`, default `20261601` per Addendum A Item 11).
          - For `STATE_HVG`: pass `2,000`-gene pre-normalized `X_hvg` (`ctrl_raw_gene_names=None`),
            or full-transcriptome raw integer UMI counts + `ctrl_raw_gene_names` (preprocessed via
            `self.hvg_norm_mode`, default `official_median_log1p`).
        """
        fold_dirs = self._resolve_fold_dirs()
        with open(fold_dirs[0] / "var_dims.pkl", "rb") as f:
            g_model = list(pickle.load(f)["gene_names"])

        if self.variant == "STATE_SE":
            if ctrl_is_embedded_x_state:
                ctrl_in = np.asarray(ctrl_cells, dtype=np.float32)
            else:
                if ctrl_raw_gene_names is None:
                    raise ValueError("STATE_SE requires ctrl_raw_gene_names when ctrl_is_embedded_x_state=False.")
                se_dir = self.se_weights_path.parent if (self.se_weights_path and self.se_weights_path.is_file()) else (self.se_weights_path or (self.models_dir / "SE-600M"))
                ctrl_in, _ = encode_controls_se600m_r16(
                    raw_counts=ctrl_cells,
                    raw_gene_names=ctrl_raw_gene_names,
                    se_dir=se_dir,
                    device=self.device,
                    batch_size=16,
                    dataset_name=control_line_label,
                    seed=int(se_encode_seed),
                )
        else:
            if ctrl_raw_gene_names is not None:
                ctrl_in, _ = preprocess_hvg_counts_r16(
                    raw_counts=ctrl_cells,
                    raw_gene_names=ctrl_raw_gene_names,
                    g_model_names=g_model,
                    norm_mode=self.hvg_norm_mode,
                    canon_map=gene_canon_map,
                )
            else:
                ctrl_in = np.asarray(ctrl_cells, dtype=np.float32)

        pert_keys_model = [
            pert_canon_to_model[p] if (pert_canon_to_model and p in pert_canon_to_model) else p
            for p in pert_names
        ]

        # Build canonical column mapping from g_model -> requested eval gene_names
        g_model_canon_to_idx: Dict[str, int] = {}
        for idx, gm in enumerate(g_model):
            cgm = gene_canon_map.get(gm, gm) if gene_canon_map is not None else gm
            g_model_canon_to_idx.setdefault(cgm, idx)
        col_sel = []
        for g in gene_names:
            cg = gene_canon_map.get(g, g) if gene_canon_map is not None else g
            if cg not in g_model_canon_to_idx:
                raise KeyError(f"Requested eval gene {g!r} (canon={cg!r}) absent from G_model.")
            col_sel.append(g_model_canon_to_idx[cg])
        col_sel_arr = np.asarray(col_sel, dtype=np.int64)

        per_fold_out: Dict[str, np.ndarray] = {}
        for fdir in fold_dirs:
            fold_name = f"zeroshot/{fdir.name}"
            d2000, _ = self._infer_single_fold(
                fold_dir=fdir,
                variant=self.variant,
                ctrl_tensor_np=ctrl_in,
                pert_keys_model=pert_keys_model,
                batch_mode=self.batch_mode,
                device_str=self.device,
                native_batch_labels=native_batch_labels,
            )
            per_fold_out[fold_name] = d2000[:, col_sel_arr]

        ens_mean = np.mean([per_fold_out[k] for k in sorted(per_fold_out.keys())], axis=0)
        if return_per_fold:
            return ens_mean, per_fold_out
        return ens_mean


class CellularAIOpenCodePublicRefitAdapter(BasePerturbModelAdapter):
    """
    DROPPED PRE-DATA (`PREREG_I Addendum A Item 5`):
    `CellularAI_OpenCode_PublicRefit` will NOT be executed in Round 16 to eliminate
    optional-analysis degrees of freedom. Retained only as a disabled stub preserving
    leakage-guard unit coverage.
    """

    def __init__(
        self,
        catalog_json_path: Optional[Path] = None,
        status_dir: Path = DEFAULT_STATUS_DIR,
    ):
        spec = ModelAuditSpec(
            model_key="CellularAI_OpenCode_PublicRefit",
            role="dropped_per_addendum_a",
            license_id="Apache-2.0",
            checkpoint_id=None,
            supports_gene_perm=True,
            n_ctrl_infer=400,
            cell_set_len=None,
        )
        super().__init__(spec=spec, status_dir=status_dir)
        self.catalog_json_path = catalog_json_path

    def _predict_impl(
        self,
        target_line: str,
        control_line_label: str,
        ctrl_cells: np.ndarray,
        pert_names: Sequence[str],
        gene_names: Sequence[str],
        **kwargs: Any,
    ) -> np.ndarray:
        raise NotImplementedError(
            "CellularAI_OpenCode_PublicRefit was explicitly dropped pre-data per "
            "PREREG_I_ADDENDUM_A.md Item 5 and will NOT be executed in Round 16."
        )


def build_closed_form_reference_preds(
    mu_target: np.ndarray,
    fc_bar: np.ndarray,
    delta_transfer: np.ndarray,
    source_grand_mean: np.ndarray,
    global_scale_alpha: float = 1.35,
) -> Dict[str, np.ndarray]:
    """Build exact closed-form reference predictors for parity and unit verification."""
    mu_t = np.maximum(np.asarray(mu_target, dtype=np.float64), 0.0)
    fc = np.maximum(np.asarray(fc_bar, dtype=np.float64), 0.0)
    tm = np.asarray(delta_transfer, dtype=np.float64)
    b_s = np.asarray(source_grand_mean, dtype=np.float64)
    c_med = float(np.median(mu_t))

    return {
        "FCT": np.log1p(mu_t[None, :] * fc) - np.log1p(mu_t[None, :]),
        "E1_ConstMed": np.log1p(c_med * fc) - np.log1p(c_med),
        "Transfer_Mean": tm.copy(),
        "Transfer_Mean_GlobalScale": float(global_scale_alpha) * tm,
        "Source_Grand_Mean": np.broadcast_to(b_s[None, :], tm.shape).copy(),
    }
