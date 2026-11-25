#!/usr/bin/env python3
"""
Clean-Room Round 16 Control-Cell Sampling & Stratified Mismatched Pool Builder
(`A-T3`, strictly following `PREREG_I.md v4 §2.2–§2.3`, `seed = 20261601`).

Never reads or imports `round16/src/`.
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

DEFAULT_R16_SEED = 20261601
N_CTRL_DEFAULT = 512
CELL_SET_LEN_DEFAULT = 32
MIN_CTRL_FLOOR_DEFAULT = 128

REPLOGLE4_ALL_LINES: Tuple[str, ...] = (
    "nadig_hepg2",
    "nadig_jurkat",
    "replogle_k562_essential",
    "replogle_rpe1",
)
REPLOGLE4_LOO_PRIMARY_LINES: Tuple[str, ...] = ("nadig_hepg2", "nadig_jurkat")
GROUP2_PRIMARY_LINES: Tuple[str, ...] = ("a549_ko", "dld1", "cd4t_rest")
GROUP2_UNSEEN_RETAINED_SOURCES: Tuple[str, ...] = ("feng_ipsc",)
GROUP2_CRISPRKO_LINES: Tuple[str, ...] = ("a549_ko",)
GROUP2_CRISPRI_LINES: Tuple[str, ...] = ("dld1", "cd4t_rest")
EXCLUDED_NO_LOCAL_RAW_LINES: Tuple[str, ...] = ("nourreddine_kolf21j",)
OPTIONAL_EXTENSION_LINES: Tuple[str, ...] = (
    "replogle_k562_essential",
    "replogle_rpe1",
    "vcc2025_h1",
)
LOCAL_RAW_CONTROL_SOURCES: Tuple[str, ...] = (
    "nadig_hepg2",
    "nadig_jurkat",
    "a549_ko",
    "dld1",
    "cd4t_rest",
    "feng_ipsc",
)


def libnorm_log1p(counts: np.ndarray, target_sum: float = 1e4) -> np.ndarray:
    """Library-size normalize each row cell to `target_sum` and apply `log1p`."""
    x = np.asarray(counts, dtype=np.float64)
    if x.ndim != 2:
        raise ValueError(f"Expected 2D counts matrix (cells x genes), got shape {x.shape}")
    row_sums = np.sum(x, axis=1, keepdims=True)
    safe_sums = np.where(row_sums > 0.0, row_sums, 1.0)
    normed = np.where(row_sums > 0.0, (x / safe_sums) * float(target_sum), 0.0)
    return np.log1p(np.maximum(normed, 0.0))


def _allocate_largest_remainder_exact(
    batch_sizes: Dict[str, int],
    n_target: int,
) -> Dict[str, int]:
    """
    Exact integer largest-remainder (Hamilton) allocation for proportional batch quotas,
    with deterministic lexicographic ascending batch-label tie-breaking (`PREREG_I v4 §2.2`).
    """
    total = sum(batch_sizes.values())
    if n_target <= 0 or total <= 0:
        return {b: 0 for b in batch_sizes}
    if n_target >= total:
        return dict(batch_sizes)

    ordered = sorted(batch_sizes.keys())
    alloc: Dict[str, int] = {}
    rem_num: Dict[str, int] = {}
    for b in ordered:
        nb = batch_sizes[b]
        exact_num = n_target * nb
        floor_q = exact_num // total
        capped = min(nb, floor_q)
        alloc[b] = capped
        rem_num[b] = (exact_num - capped * total) if capped < nb else -1

    deficit = n_target - sum(alloc.values())
    while deficit > 0:
        eligible = [b for b in ordered if alloc[b] < batch_sizes[b]]
        if not eligible:
            break
        eligible.sort(key=lambda b: (-rem_num[b], b))
        chosen = eligible[0]
        alloc[chosen] += 1
        rem_num[chosen] = -1
        deficit -= 1
    return alloc


def sample_stratified_control_sets_r16(
    n_cells: Optional[int] = None,
    qc_mask: Optional[Sequence[bool]] = None,
    batch_labels: Optional[Sequence[Any]] = None,
    n_ctrl: int = N_CTRL_DEFAULT,
    cell_set_len: int = CELL_SET_LEN_DEFAULT,
    min_ctrl_floor: int = MIN_CTRL_FLOOR_DEFAULT,
    seed: int = DEFAULT_R16_SEED,
) -> Dict[str, Any]:
    """
    Deterministically sample control-cell sets per `PREREG_I.md v4 §2.2`:
      1. Filter to QC-passing NTC indices.
      2. If `N_NTC < 128`, return `excluded=True` (zero sets).
      3. If `128 <= N_NTC < 512`, take all QC-valid NTCs without replacement, shuffle once under
         `seed`, and truncate to `(N_NTC // 32) * 32` (`N_NTC // 32` sets of `32`).
      4. If `N_NTC >= 512`, draw `512` NTCs without replacement, proportionally stratified by batch
         via exact largest-remainder with lexicographic tie-breaking (or unstratified if no batch),
         shuffle once under the same RNG, and partition into `16` disjoint sets of `32`.
    """
    if qc_mask is not None:
        mask_arr = np.asarray(qc_mask, dtype=bool)
        if mask_arr.ndim != 1:
            raise ValueError(f"qc_mask must be 1D, got shape {mask_arr.shape}")
        if n_cells is not None and int(n_cells) != mask_arr.size:
            raise ValueError(f"n_cells ({n_cells}) != len(qc_mask) ({mask_arr.size})")
        valid_idx = np.flatnonzero(mask_arr)
    else:
        if n_cells is None or int(n_cells) < 0:
            raise ValueError("Provide non-negative `n_cells` or `qc_mask`.")
        valid_idx = np.arange(int(n_cells), dtype=np.int64)

    n_valid = int(valid_idx.size)
    if batch_labels is not None:
        bl_arr = np.asarray([str(x) for x in batch_labels], dtype=object)
        expected_len = int(mask_arr.size) if qc_mask is not None else int(n_cells or 0)
        if bl_arr.size != expected_len:
            raise ValueError(f"len(batch_labels)={bl_arr.size} != total cells {expected_len}")
        valid_batches = bl_arr[valid_idx]
        total_batch_counts: Dict[str, int] = {}
        batch_pools: Dict[str, List[int]] = {}
        for g_idx, b_val in zip(valid_idx.tolist(), valid_batches.tolist()):
            b_str = str(b_val)
            batch_pools.setdefault(b_str, []).append(int(g_idx))
            total_batch_counts[b_str] = total_batch_counts.get(b_str, 0) + 1
        batch_pools = {k: sorted(v) for k, v in sorted(batch_pools.items())}
        total_batch_counts = dict(sorted(total_batch_counts.items()))
    else:
        total_batch_counts = {"__ALL__": n_valid}
        batch_pools = {"__ALL__": valid_idx.tolist()}

    if n_valid < int(min_ctrl_floor):
        return {
            "seed": int(seed),
            "n_ntc_valid": n_valid,
            "n_ctrl_requested": int(n_ctrl),
            "cell_set_len": int(cell_set_len),
            "min_ctrl_floor": int(min_ctrl_floor),
            "excluded": True,
            "exclusion_reason": f"N_NTC={n_valid} < floor {min_ctrl_floor}",
            "sampling_mode": "excluded_below_floor",
            "n_sampled": 0,
            "n_sets": 0,
            "sampling_fraction": 0.0,
            "selected_indices": np.zeros(0, dtype=np.int64),
            "sets": np.zeros((0, int(cell_set_len)), dtype=np.int64),
            "batch_counts_total": total_batch_counts,
            "batch_counts_sampled": {k: 0 for k in total_batch_counts},
        }

    rng = np.random.default_rng(int(seed))
    if n_valid < int(n_ctrl):
        n_keep = (n_valid // int(cell_set_len)) * int(cell_set_len)
        shuffled = valid_idx[rng.permutation(n_valid)]
        selected = shuffled[:n_keep].astype(np.int64)
        mode = "small_ntc_truncated"
        if batch_labels is not None:
            sel_set = set(selected.tolist())
            sampled_batch_counts = {
                b: sum(1 for gi in batch_pools[b] if gi in sel_set) for b in total_batch_counts
            }
        else:
            sampled_batch_counts = {"__ALL__": int(n_keep)}
    else:
        n_keep = int(n_ctrl)
        if n_keep % int(cell_set_len) != 0:
            raise ValueError(f"n_ctrl ({n_keep}) must be divisible by cell_set_len ({cell_set_len}).")
        if batch_labels is not None and len(batch_pools) > 1:
            alloc = _allocate_largest_remainder_exact(total_batch_counts, n_keep)
            chosen_parts: List[np.ndarray] = []
            for b in sorted(alloc.keys()):
                kb = alloc[b]
                if kb <= 0:
                    continue
                pool = np.asarray(batch_pools[b], dtype=np.int64)
                picked = rng.choice(pool, size=kb, replace=False)
                chosen_parts.append(picked)
            selected_unshuffled = np.concatenate(chosen_parts)
            selected = selected_unshuffled[rng.permutation(selected_unshuffled.size)].astype(np.int64)
            sampled_batch_counts = dict(sorted(alloc.items()))
            mode = "stratified_by_batch"
        else:
            picked = rng.choice(valid_idx, size=n_keep, replace=False)
            selected = picked[rng.permutation(picked.size)].astype(np.int64)
            key = next(iter(total_batch_counts.keys()))
            sampled_batch_counts = {key: int(n_keep)}
            mode = "unstratified_single_batch"

    n_sets = int(selected.size) // int(cell_set_len)
    sets = selected.reshape(n_sets, int(cell_set_len))
    return {
        "seed": int(seed),
        "n_ntc_valid": n_valid,
        "n_ctrl_requested": int(n_ctrl),
        "cell_set_len": int(cell_set_len),
        "min_ctrl_floor": int(min_ctrl_floor),
        "excluded": False,
        "exclusion_reason": None,
        "sampling_mode": mode,
        "n_sampled": int(selected.size),
        "n_sets": n_sets,
        "sampling_fraction": float(selected.size) / float(n_valid),
        "selected_indices": selected,
        "sets": sets,
        "batch_counts_total": total_batch_counts,
        "batch_counts_sampled": sampled_batch_counts,
    }


def _ckpt_seen_training_lines(checkpoint_id: str) -> Tuple[str, ...]:
    """Return Replogle-4 training lines seen by `checkpoint_id`."""
    ck = str(checkpoint_id).strip().lower()
    if ck.endswith("/ensemble4") or ck == "ensemble4":
        return REPLOGLE4_ALL_LINES
    for short, full in (
        ("hepg2", "nadig_hepg2"),
        ("jurkat", "nadig_jurkat"),
        ("k562", "replogle_k562_essential"),
        ("rpe1", "replogle_rpe1"),
    ):
        if ck.endswith(f"/{short}") or ck == short:
            return tuple(x for x in REPLOGLE4_ALL_LINES if x != full)
    raise ValueError(f"Unrecognized checkpoint_id={checkpoint_id!r}")


def build_mismatched_strata_r16(
    target_line: str,
    checkpoint_id: str,
    group2_qualifying_lines: Sequence[str] = GROUP2_PRIMARY_LINES,
    retained_unseen_sources: Sequence[str] = GROUP2_UNSEEN_RETAINED_SOURCES,
    allow_optional_extensions: Sequence[str] = (),
) -> Dict[str, Tuple[str, ...]]:
    """
    Build disjoint (`seen_by_checkpoint`, `unseen_by_checkpoint`) mismatched control pools
    strictly following `PREREG_I.md v4 §2.3` (`I14`, `Q4`):
      - Rejects `nourreddine_kolf21j` and non-local lines unless explicitly permitted via
        `allow_optional_extensions`.
      - Target line itself is always excluded from mismatched sources.
      - Includes retained unseen source `feng_ipsc` (`Q4`) in `unseen_by_checkpoint`.
    """
    ext_set = set(allow_optional_extensions)
    allowed_local = set(LOCAL_RAW_CONTROL_SOURCES) | ext_set

    for lines_arg in (group2_qualifying_lines, retained_unseen_sources):
        for ln in lines_arg:
            if ln in EXCLUDED_NO_LOCAL_RAW_LINES:
                raise ValueError(f"Excluded non-local line cannot enter mismatched pool: {ln!r}")
            if ln in OPTIONAL_EXTENSION_LINES and ln not in ext_set:
                raise ValueError(
                    f"Optional extension line {ln!r} requires explicit human download approval."
                )
            if ln not in allowed_local:
                raise ValueError(f"Line {ln!r} lacks verified local raw single-cell counts.")

    ckpt_seen = set(_ckpt_seen_training_lines(checkpoint_id))
    local_seen_candidates = tuple(sorted((set(REPLOGLE4_LOO_PRIMARY_LINES) | ext_set) & ckpt_seen))
    seen_sources = tuple(s for s in local_seen_candidates if s != target_line)

    unseen_pool = (set(group2_qualifying_lines) | set(retained_unseen_sources)) - {target_line} - ckpt_seen
    unseen_sources = tuple(sorted(unseen_pool))

    return {
        "seen_by_checkpoint": seen_sources,
        "unseen_by_checkpoint": unseen_sources,
    }
