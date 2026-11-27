#!/usr/bin/env python3
"""
Model-free preprocessing, RHJ axis scatter, `pairwise_common_axis` (`P4`) masks,
`full_44728_api` secondary input scatter, Two-Tier perturbation-key resolution (`§3.5`),
and output gene projection for Round 16 Addendum B (`Occam`).

Zero Occam imports or vendored Occam code.
"""

from __future__ import annotations

import csv
import hashlib
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

_SRC_DIR = str(Path(__file__).resolve().parent)
if _SRC_DIR not in sys.path:
    sys.path.insert(0, _SRC_DIR)

from raw_controls_addb import strip_ensembl_version  # noqa: E402

OCCAM_TARGET_SUM = 2000.0
OCCAM_N_ENCODER_GENES = 44728
RHJ_N_AXIS_GENES = 6983
RHJ_GENE_AXIS_SHA256 = "ef64f2ba7ff7816ed947796066cd8a68014884c6176619d750ef79462c78d4c0"


@dataclass(frozen=True)
class RhjGeneAxis:
    """Frozen 6,983-gene RHJ benchmark axis loaded from `gene_axis.csv`."""

    ensembl_ids: List[str]
    gene_symbols: List[str]
    axis_sha256: str


@dataclass(frozen=True)
class PairwiseCommonAxisMask:
    """Pre-registered `P4` pairwise common input mask between target `t` and source `s'`."""

    target_line: str
    source_line: str
    mask: np.ndarray  # bool array of length 6,983
    n_axis_total: int
    n_common: int
    n_target_present: int
    n_source_present: int
    n_target_only_zeroed: int
    n_source_only_zeroed: int
    n_union_absent: int
    requires_common_axis_arm: bool
    mask_sha256: str


@dataclass(frozen=True)
class PerturbationResolutionResult:
    """Result of Two-Tier (`rhj_fewshot` eligible -> direct `pert2id`) resolution (`§3.5`)."""

    ordered_targets: List[str]
    resolved_keys: List[str]
    resolved_ids: List[int]
    tiers: List[int]  # 1 or 2 per target
    tier1_targets: List[str]
    tier2_targets: List[str]
    bridged_pairs: Dict[str, str]  # canonical target -> distinct legacy Occam key
    unresolved_targets: List[str]


def compute_string_list_sha256(items: Sequence[str]) -> str:
    """Exact newline-joined UTF-8 SHA-256 used by Occam's `gene_axis_sha256` and `cohort_sha256`."""
    return hashlib.sha256("\n".join(str(x) for x in items).encode("utf-8")).hexdigest()


def load_rhj_gene_axis(csv_path: Path, *, verify_sha256: bool = True) -> RhjGeneAxis:
    """Load `benchmarks/rhj_fewshot/gene_axis.csv` and verify count (`6,983`) and SHA-256."""
    ens_ids: List[str] = []
    symbols: List[str] = []
    with open(csv_path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames or "ensembl_gene_id" not in reader.fieldnames:
            raise KeyError(f"{csv_path} missing 'ensembl_gene_id' column: {reader.fieldnames}")
        sym_col = "gene_name" if "gene_name" in reader.fieldnames else None
        for row in reader:
            ens_ids.append(strip_ensembl_version(row["ensembl_gene_id"]))
            symbols.append(str(row[sym_col]).strip() if sym_col else "")

    if len(set(ens_ids)) != len(ens_ids):
        raise ValueError(f"{csv_path} contains duplicate Ensembl IDs.")
    digest = compute_string_list_sha256(ens_ids)
    if verify_sha256:
        if len(ens_ids) != RHJ_N_AXIS_GENES:
            raise ValueError(f"Expected {RHJ_N_AXIS_GENES} RHJ axis genes, got {len(ens_ids)}.")
        if digest != RHJ_GENE_AXIS_SHA256:
            raise ValueError(
                f"RHJ gene_axis_sha256 mismatch: got {digest}, expected {RHJ_GENE_AXIS_SHA256}."
            )
    return RhjGeneAxis(ensembl_ids=ens_ids, gene_symbols=symbols, axis_sha256=digest)


def dedup_leading_ids(raw_ids: Sequence[str]) -> Tuple[np.ndarray, List[str], List[str]]:
    """
    `PREREG_I_ADDENDUM_B_DRAFT_v2.md §3.3.2` (`DEVIATION_R16_D7` mirror):
    Strip `.N` version suffixes and strictly keep the LEADING occurrence of each non-empty ID
    in original raw file column order, dropping subsequent duplicates.
    """
    seen = set()
    keep_indices: List[int] = []
    unique_ids: List[str] = []
    dropped_dups: List[str] = []
    for idx, raw_id in enumerate(raw_ids):
        clean = strip_ensembl_version(raw_id)
        if not clean:
            continue
        if clean in seen:
            dropped_dups.append(clean)
            continue
        seen.add(clean)
        keep_indices.append(idx)
        unique_ids.append(clean)
    return np.asarray(keep_indices, dtype=np.int64), unique_ids, dropped_dups


def normalize_counts_target_sum_log1p(
    raw_counts: np.ndarray,
    *,
    target_sum: float = OCCAM_TARGET_SUM,
    full_row_sums: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, int]:
    """
    Normalize non-negative raw UMI counts to `target_sum` (`2000.0`) using full-panel per-cell
    library sizes (`float64` sum) followed by `log1p` in `float32` (`§3.3.3`).

    Zero-UMI cells (`L_i == 0`) remain `0.0` (never `NaN`/`Inf`). Negative or non-finite raw
    counts raise `ValueError`.
    """
    x = np.asarray(raw_counts, dtype=np.float32)
    if x.ndim != 2:
        raise ValueError(f"Expected 2D raw_counts array, got shape {x.shape}.")
    if not np.all(np.isfinite(x)):
        raise ValueError("raw_counts contains non-finite values.")
    if x.size > 0 and float(x.min()) < 0.0:
        raise ValueError(f"raw_counts contains negative values (min={float(x.min())}).")

    if full_row_sums is None:
        rs = x.sum(axis=1, dtype=np.float64)
    else:
        rs = np.asarray(full_row_sums, dtype=np.float64)
        if rs.shape != (x.shape[0],):
            raise ValueError(f"full_row_sums shape {rs.shape} != ({x.shape[0]},).")
        if not np.all(np.isfinite(rs)) or (rs.size > 0 and float(rs.min()) < 0.0):
            raise ValueError("full_row_sums must be finite and non-negative.")

    n_zero_lib = int(np.sum(rs <= 0.0))
    safe_rs = np.where(rs > 0.0, rs, 1.0).astype(np.float32)
    scaled = (x / safe_rs[:, None]) * np.float32(target_sum)
    if n_zero_lib > 0:
        scaled[rs <= 0.0, :] = 0.0
    normed = np.log1p(scaled).astype(np.float32, copy=False)
    return np.ascontiguousarray(normed, dtype=np.float32), n_zero_lib


def build_pairwise_common_axis_mask(
    rhj_axis_ensembl: Sequence[str],
    target_raw_ensembl: Sequence[str],
    source_raw_ensembl: Sequence[str],
    *,
    target_line: str = "",
    source_line: str = "",
) -> PairwiseCommonAxisMask:
    """
    Construct the pre-registered `P4` pairwise common input mask (`§3.3.4`):
        A_common(t, s') = A_RHJ ∩ G_raw(t) ∩ G_raw(s')
    Symmetric by construction (`mask(t, s') == mask(s', t)`).
    """
    axis = [strip_ensembl_version(g) for g in rhj_axis_ensembl]
    t_set = {strip_ensembl_version(g) for g in target_raw_ensembl if strip_ensembl_version(g)}
    s_set = {strip_ensembl_version(g) for g in source_raw_ensembl if strip_ensembl_version(g)}

    t_in = np.asarray([g in t_set for g in axis], dtype=bool)
    s_in = np.asarray([g in s_set for g in axis], dtype=bool)
    common = t_in & s_in

    n_total = len(axis)
    n_common = int(common.sum())
    n_t = int(t_in.sum())
    n_s = int(s_in.sum())
    n_t_only = int(np.sum(t_in & ~s_in))
    n_s_only = int(np.sum(s_in & ~t_in))
    n_union_absent = int(n_total - n_common)
    mask_sha256 = hashlib.sha256(np.ascontiguousarray(common, dtype=np.uint8).tobytes()).hexdigest()

    return PairwiseCommonAxisMask(
        target_line=target_line,
        source_line=source_line,
        mask=common,
        n_axis_total=n_total,
        n_common=n_common,
        n_target_present=n_t,
        n_source_present=n_s,
        n_target_only_zeroed=n_t_only,
        n_source_only_zeroed=n_s_only,
        n_union_absent=n_union_absent,
        requires_common_axis_arm=(n_union_absent > 0),
        mask_sha256=mask_sha256,
    )


def scatter_rhj_6983_to_encoder(
    raw_counts: np.ndarray,
    raw_ensembl_ids: Sequence[str],
    rhj_axis_ensembl: Sequence[str],
    encoder_gene_names: Sequence[str],
    *,
    target_sum: float = OCCAM_TARGET_SUM,
    active_rhj_mask: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    Locked primary input preprocessing (`rhj_6983_axis`, `§3.3.3`) and optional `P4`
    `pairwise_common_axis` masking (`§3.3.4`):
      1. Compute per-cell library size `L_i` over ALL raw measured `Gene Expression` features;
      2. Keep strictly the leading occurrence of each version-stripped Ensembl ID (`§3.3.2`);
      3. Select RHJ axis columns present in `raw_ensembl_ids` (and True in `active_rhj_mask`,
         when supplied), normalize `log1p(x_raw / L_i * target_sum)` in `float32`;
      4. Scatter into exact columns of the `len(encoder_gene_names)` (`44,728`) encoder tensor,
         zeroing all non-axis columns, absent axis genes, and masked-out axis genes.
    """
    x_raw = np.asarray(raw_counts, dtype=np.float32)
    if x_raw.ndim != 2 or x_raw.shape[1] != len(raw_ensembl_ids):
        raise ValueError(
            f"raw_counts shape {x_raw.shape} incompatible with len(raw_ensembl_ids)={len(raw_ensembl_ids)}."
        )
    if not np.all(np.isfinite(x_raw)):
        raise ValueError("raw_counts contains non-finite values.")
    if x_raw.size > 0 and float(x_raw.min()) < 0.0:
        raise ValueError(f"raw_counts contains negative values (min={float(x_raw.min())}).")

    axis = [strip_ensembl_version(g) for g in rhj_axis_ensembl]
    enc_genes = [strip_ensembl_version(g) for g in encoder_gene_names]
    if len(set(axis)) != len(axis):
        raise ValueError("rhj_axis_ensembl contains duplicate Ensembl IDs.")
    if len(set(enc_genes)) != len(enc_genes):
        raise ValueError("encoder_gene_names contains duplicate Ensembl IDs.")

    enc_pos = {g: i for i, g in enumerate(enc_genes)}
    missing_from_enc = [g for g in axis if g not in enc_pos]
    if missing_from_enc:
        raise ValueError(
            f"{len(missing_from_enc)} RHJ axis gene(s) absent from encoder_gene_names "
            f"(e.g. {missing_from_enc[:5]})."
        )

    if active_rhj_mask is not None:
        mask = np.asarray(active_rhj_mask, dtype=bool)
        if mask.shape != (len(axis),):
            raise ValueError(f"active_rhj_mask shape {mask.shape} != ({len(axis)},).")
    else:
        mask = np.ones(len(axis), dtype=bool)

    # Full-panel library size over ALL raw measured Gene Expression features (before axis restriction)
    full_rs = x_raw.sum(axis=1, dtype=np.float64)

    # Leading-occurrence deduplication on raw Ensembl IDs (D7 mirror)
    keep_raw_idx, uniq_raw_ids, dropped_dups = dedup_leading_ids(raw_ensembl_ids)
    raw_pos = {g: int(keep_raw_idx[k]) for k, g in enumerate(uniq_raw_ids)}

    src_cols: List[int] = []
    dst_cols: List[int] = []
    present_rhj = 0
    absent_rhj: List[str] = []
    zeroed_by_mask = 0

    for ax_i, ens in enumerate(axis):
        r_col = raw_pos.get(ens)
        if r_col is None:
            absent_rhj.append(ens)
            continue
        present_rhj += 1
        if not mask[ax_i]:
            zeroed_by_mask += 1
            continue
        src_cols.append(r_col)
        dst_cols.append(enc_pos[ens])

    n_cells = x_raw.shape[0]
    n_enc = len(enc_genes)
    x_enc = np.zeros((n_cells, n_enc), dtype=np.float32)

    if src_cols:
        sub_raw = x_raw[:, np.asarray(src_cols, dtype=np.int64)]
        sub_norm, n_zero_lib = normalize_counts_target_sum_log1p(
            sub_raw, target_sum=target_sum, full_row_sums=full_rs
        )
        x_enc[:, np.asarray(dst_cols, dtype=np.int64)] = sub_norm
    else:
        n_zero_lib = int(np.sum(full_rs <= 0.0))

    mask_sha256 = hashlib.sha256(np.ascontiguousarray(mask, dtype=np.uint8).tobytes()).hexdigest()
    audit = {
        "input_axis_mode": "pairwise_common_axis" if active_rhj_mask is not None else "rhj_6983_axis",
        "target_sum": float(target_sum),
        "n_cells": int(n_cells),
        "n_raw_genes": int(len(raw_ensembl_ids)),
        "n_encoder_genes": int(n_enc),
        "n_axis_total": int(len(axis)),
        "n_present_rhj": int(present_rhj),
        "n_absent_rhj": int(len(absent_rhj)),
        "n_active_rhj_populated": int(len(src_cols)),
        "n_present_but_zeroed_by_mask": int(zeroed_by_mask),
        "n_dropped_duplicate_raw_ensembl": int(len(dropped_dups)),
        "n_zero_libsize_cells": int(n_zero_lib),
        "active_rhj_mask_sha256": mask_sha256,
        "all_finite": bool(np.all(np.isfinite(x_enc))),
    }
    return np.ascontiguousarray(x_enc, dtype=np.float32), audit


def scatter_full_44728_api(
    raw_counts: np.ndarray,
    raw_ensembl_ids: Sequence[str],
    encoder_gene_names: Sequence[str],
    *,
    target_sum: float = OCCAM_TARGET_SUM,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    Secondary descriptive input-axis sensitivity (`full_44728_api`, `§3.3.5`):
    Matches Occam's `GeneSpace.prepare` convention — populate ALL leading-occurrence Ensembl genes
    matching `encoder_gene_names` (`44,728`), computing per-cell row sums over the `44,728`-matched
    columns, scaling to `target_sum = 2000.0`, and applying `log1p`.
    """
    x_raw = np.asarray(raw_counts, dtype=np.float32)
    if x_raw.ndim != 2 or x_raw.shape[1] != len(raw_ensembl_ids):
        raise ValueError(
            f"raw_counts shape {x_raw.shape} incompatible with len(raw_ensembl_ids)={len(raw_ensembl_ids)}."
        )
    if not np.all(np.isfinite(x_raw)):
        raise ValueError("raw_counts contains non-finite values.")
    if x_raw.size > 0 and float(x_raw.min()) < 0.0:
        raise ValueError(f"raw_counts contains negative values (min={float(x_raw.min())}).")

    enc_genes = [strip_ensembl_version(g) for g in encoder_gene_names]
    if len(set(enc_genes)) != len(enc_genes):
        raise ValueError("encoder_gene_names contains duplicate Ensembl IDs.")
    enc_pos = {g: i for i, g in enumerate(enc_genes)}

    keep_raw_idx, uniq_raw_ids, dropped_dups = dedup_leading_ids(raw_ensembl_ids)
    src_cols: List[int] = []
    dst_cols: List[int] = []
    for k, ens in enumerate(uniq_raw_ids):
        e_col = enc_pos.get(ens)
        if e_col is not None:
            src_cols.append(int(keep_raw_idx[k]))
            dst_cols.append(e_col)

    n_cells = x_raw.shape[0]
    n_enc = len(enc_genes)
    x_enc = np.zeros((n_cells, n_enc), dtype=np.float32)

    if src_cols:
        sub_raw = x_raw[:, np.asarray(src_cols, dtype=np.int64)]
        # Row sum strictly over the 44,728-matched Ensembl genes (GeneSpace.prepare convention)
        matched_rs = sub_raw.sum(axis=1, dtype=np.float64)
        sub_norm, n_zero_lib = normalize_counts_target_sum_log1p(
            sub_raw, target_sum=target_sum, full_row_sums=matched_rs
        )
        x_enc[:, np.asarray(dst_cols, dtype=np.int64)] = sub_norm
    else:
        n_zero_lib = n_cells

    audit = {
        "input_axis_mode": "full_44728_api",
        "target_sum": float(target_sum),
        "n_cells": int(n_cells),
        "n_raw_genes": int(len(raw_ensembl_ids)),
        "n_encoder_genes": int(n_enc),
        "n_matched_encoder_genes": int(len(src_cols)),
        "n_dropped_duplicate_raw_ensembl": int(len(dropped_dups)),
        "n_zero_libsize_cells": int(n_zero_lib),
        "all_finite": bool(np.all(np.isfinite(x_enc))),
    }
    return np.ascontiguousarray(x_enc, dtype=np.float32), audit


def load_rhj_perturbations_table(csv_path: Path) -> Dict[str, Any]:
    """
    Parse `benchmarks/rhj_fewshot/perturbations.csv` returning:
      - `eligible_keys`: sorted `2,381` eligible perturbation keys (`eligible == "1"`)
      - `test_keys`: sorted `621` held-out HepG2 test perturbation keys (`split == "test"`, `eligible == "1"`)
      - `train_keys`: sorted `1,657` HepG2 training perturbation keys (`split == "train"`, `eligible == "1"`)
      - `val_keys`: sorted `103` HepG2 validation perturbation keys (`split == "val"`, `eligible == "1"`)
      - `pert2id`: `{gene_name: int(pert2id)}` for all rows with valid `pert2id`
    """
    eligible_keys: List[str] = []
    test_keys: List[str] = []
    train_keys: List[str] = []
    val_keys: List[str] = []
    pert2id: Dict[str, int] = {}

    with open(csv_path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            g = str(row["gene_name"]).strip()
            split = str(row.get("split", "")).strip()
            elig = str(row.get("eligible", "")).strip() == "1"
            pid_raw = str(row.get("pert2id", "")).strip()
            if pid_raw and pid_raw.lstrip("-").isdigit():
                pert2id[g] = int(pid_raw)
            if elig:
                eligible_keys.append(g)
                if split == "test":
                    test_keys.append(g)
                elif split == "train":
                    train_keys.append(g)
                elif split == "val":
                    val_keys.append(g)

    return {
        "eligible_keys": sorted(eligible_keys),
        "test_keys": sorted(test_keys),
        "train_keys": sorted(train_keys),
        "val_keys": sorted(val_keys),
        "pert2id": pert2id,
    }


def load_master_pert2id_csv(csv_path: Path) -> Dict[str, int]:
    """Load `master_gene_set_sorted.csv` (`gene_id,gene_name`) into `{gene_name: int(gene_id)}`."""
    out: Dict[str, int] = {}
    with open(csv_path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            g = str(row["gene_name"]).strip()
            gid = int(str(row["gene_id"]).strip())
            if g in out:
                raise ValueError(f"Duplicate gene_name {g!r} in {csv_path}")
            out[g] = gid
    return out


def resolve_perturbation_keys_two_tier(
    target_names: Sequence[str],
    rhj_eligible_keys: Sequence[str],
    ckpt_pert2id: Mapping[str, int],
    canonicalize_fn: Callable[[str], Optional[str]],
    *,
    require_all_resolved: bool = True,
) -> PerturbationResolutionResult:
    """
    Locked Two-Tier Perturbation-Key Resolution (`PREREG_I_ADDENDUM_B_DRAFT_v2.md §3.5`):

    Occam's `18,216`-gene `pert2id` contains `317` duplicate HGNC pairs where both a legacy
    symbol (e.g. `AARS`, `DARS`, `TAZ`, `H2AFZ`) and the approved HGNC symbol (`AARS1`, `DARS1`,
    `TAFAZZIN`, `H2AZ1`) exist as separate nodes, while `rhj_fewshot` (`2,381` eligible keys,
    `0` HGNC collisions) trained ONLY the legacy node.

      - **Tier 1 (Trained `rhj_fewshot` priority):** Map canonical target `p` to the unique key
        `k` in `rhj_eligible_keys` with `canon(k) == canon(p)` (verifying `k in ckpt_pert2id`).
      - **Tier 2 (Direct `pert2id` fallback):** If absent from Tier 1, resolve `p` directly in
        `ckpt_pert2id` (exact `p` first, then `canon(p)` if exact string differs).
    """
    def _canon(sym: str) -> str:
        c = canonicalize_fn(sym)
        return str(c).strip() if c else str(sym).strip()

    # Build Tier-1 canonical map and assert ZERO collisions among rhj_eligible_keys
    tier1_canon_to_key: Dict[str, str] = {}
    collisions: Dict[str, List[str]] = {}
    for k in sorted(set(str(x).strip() for x in rhj_eligible_keys)):
        ck = _canon(k)
        if ck in tier1_canon_to_key and tier1_canon_to_key[ck] != k:
            collisions.setdefault(ck, [tier1_canon_to_key[ck]]).append(k)
        else:
            tier1_canon_to_key[ck] = k
    if collisions:
        raise ValueError(
            f"Tier-1 rhj_eligible_keys has {len(collisions)} HGNC collision(s): "
            f"{dict(list(collisions.items())[:5])}"
        )

    # Validate every Tier-1 key exists in ckpt_pert2id
    missing_t1 = [k for k in tier1_canon_to_key.values() if k not in ckpt_pert2id]
    if missing_t1:
        raise KeyError(
            f"{len(missing_t1)} rhj_eligible_keys absent from ckpt_pert2id (e.g. {missing_t1[:5]})."
        )

    ordered_targets: List[str] = []
    resolved_keys: List[str] = []
    resolved_ids: List[int] = []
    tiers: List[int] = []
    tier1_targets: List[str] = []
    tier2_targets: List[str] = []
    bridged_pairs: Dict[str, str] = {}
    unresolved: List[str] = []

    for raw_t in target_names:
        t = str(raw_t).strip()
        ct = _canon(t)
        if ct in tier1_canon_to_key:
            k = tier1_canon_to_key[ct]
            ordered_targets.append(t)
            resolved_keys.append(k)
            resolved_ids.append(int(ckpt_pert2id[k]))
            tiers.append(1)
            tier1_targets.append(t)
            if k != t:
                bridged_pairs[t] = k
            continue

        # Tier 2: direct pert2id fallback (exact string first, then canonical string)
        t2_key: Optional[str] = None
        if t in ckpt_pert2id:
            t2_key = t
        elif ct in ckpt_pert2id:
            t2_key = ct

        if t2_key is not None:
            ordered_targets.append(t)
            resolved_keys.append(t2_key)
            resolved_ids.append(int(ckpt_pert2id[t2_key]))
            tiers.append(2)
            tier2_targets.append(t)
            if t2_key != t:
                bridged_pairs[t] = t2_key
        else:
            unresolved.append(t)

    if require_all_resolved and unresolved:
        raise KeyError(
            f"{len(unresolved)} target(s) unresolved by Two-Tier pert2id rule: {unresolved[:10]}"
        )

    return PerturbationResolutionResult(
        ordered_targets=ordered_targets,
        resolved_keys=resolved_keys,
        resolved_ids=resolved_ids,
        tiers=tiers,
        tier1_targets=tier1_targets,
        tier2_targets=tier2_targets,
        bridged_pairs=bridged_pairs,
        unresolved_targets=unresolved,
    )


def build_output_gene_projection(
    rhj_axis: RhjGeneAxis,
    encoder_gene_names: Sequence[str],
    target_gene_names: Sequence[str],
    canonicalize_fn: Callable[[str], Optional[str]],
    *,
    require_full_coverage: bool = True,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    Map Occam's `6,983` RHJ outputs onto the exact ordered columns of `G_eval^R16(t)` (`§3.5`).

    Canonicalizes `rhj_axis.gene_symbols` via frozen R10 `HGNCBridge` (asserting `0` HGNC
    collisions across the `6,983` RHJ axis genes) and returns the exact `int64` column indices
    into `encoder_gene_names` (`44,728`) corresponding to `target_gene_names`.
    """
    def _canon(sym: str) -> str:
        c = canonicalize_fn(sym)
        return str(c).strip() if c else str(sym).strip()

    enc_pos = {strip_ensembl_version(g): i for i, g in enumerate(encoder_gene_names)}
    canon_to_enc_col: Dict[str, int] = {}
    canon_to_ens: Dict[str, str] = {}
    collisions: Dict[str, List[str]] = {}

    for ens, sym in zip(rhj_axis.ensembl_ids, rhj_axis.gene_symbols):
        if ens not in enc_pos:
            raise KeyError(f"RHJ axis Ensembl ID {ens!r} missing from encoder_gene_names.")
        csym = _canon(sym)
        if csym in canon_to_enc_col and canon_to_ens[csym] != ens:
            collisions.setdefault(csym, [canon_to_ens[csym]]).append(ens)
        else:
            canon_to_enc_col[csym] = enc_pos[ens]
            canon_to_ens[csym] = ens

    if collisions:
        raise ValueError(
            f"RHJ axis canonicalization yielded {len(collisions)} HGNC collision(s): "
            f"{dict(list(collisions.items())[:5])}"
        )

    out_cols: List[int] = []
    missing: List[str] = []
    for g in target_gene_names:
        cg = _canon(str(g).strip())
        if cg in canon_to_enc_col:
            out_cols.append(canon_to_enc_col[cg])
        else:
            missing.append(str(g).strip())

    if require_full_coverage and missing:
        raise KeyError(
            f"{len(missing)}/{len(target_gene_names)} target gene(s) absent from RHJ 6,983 axis: "
            f"{missing[:10]}"
        )

    return np.asarray(out_cols, dtype=np.int64), {
        "n_target_genes": int(len(target_gene_names)),
        "n_mapped": int(len(out_cols)),
        "n_missing": int(len(missing)),
        "missing_genes": missing,
        "rhj_hgnc_collisions": 0,
    }
