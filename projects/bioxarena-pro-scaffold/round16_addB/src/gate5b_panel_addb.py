#!/usr/bin/env python3
"""
Model-free Gate 5B Split-Half Panel Builder for local `nadig_hepg2`
(`PREREG_I_ADDENDUM_B_DRAFT_v2.md §5.2`, `PREREG_I_ADDENDUM_B_LOCK.md L1`).

Strict compliance:
  - Clean-room original implementation (ZERO imports or vendoring of Occam code).
  - Lock `L1`: Decodes `obs["gene"]` and `var["gene_name"]` through `read_h5ad_table_col`
    (`DEVIATION_R16_D3` legacy AnnData categorical branch + modern categorical branch;
    raw integer codes are NEVER used as labels).
  - Hard Halt Guard (`L1`): Enforces exact counts (`611` cohort perturbations, `4,976`
    `non-targeting` controls, `42,402` panel rows) unless overridden for synthetic unit tests;
    halts immediately with `Gate5BPanelCountError` on any mismatch.
  - Deterministic per-`(split, group)` `blake2b` batch-stratified equal halves matching
    `splits_manifest.json` (`SPLIT_SALT = "occam.vcell.panel.splithalf/v1"`, `base_seed = 42`,
    `n_splits = 5`), checking all `20` split SHA-256 digests (`Case A` vs `Case B` Gate 5B
    tolerance selection) and `ctrl_batch_max_abs_diff <= 1`.
  - Model-free experimental replicate agreement scorer (`pearson_delta_mean`, `mae_mean`) on
    the built panel (author reference: `0.372174`).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import h5py
import numpy as np
import scipy.sparse as sp

_SRC_DIR = str(Path(__file__).resolve().parent)
if _SRC_DIR not in sys.path:
    sys.path.insert(0, _SRC_DIR)

from preprocess_addb import (  # noqa: E402
    OCCAM_TARGET_SUM,
    compute_string_list_sha256,
    load_rhj_gene_axis,
    load_rhj_perturbations_table,
    normalize_counts_target_sum_log1p,
)
from raw_controls_addb import (  # noqa: E402
    DEFAULT_LOCAL_NADIG_DIR,
    read_h5ad_table_col,
    strip_ensembl_version,
)

SPLITHALF_CONTROL_LABEL = "non-targeting"
SPLITHALF_SALT = "occam.vcell.panel.splithalf/v1"
SPLITHALF_CONTROL_GROUP = "__controls__"

# Locked Gate 5B expected counts (`PREREG_I_ADDENDUM_B_LOCK.md L1` & `splits_manifest.json`)
GATE5B_EXPECTED_N_COHORT_PERTS = 611
GATE5B_EXPECTED_N_CONTROLS = 4976
GATE5B_EXPECTED_N_PANEL_ROWS = 42402
GATE5B_EXPECTED_COHORT_SHA256 = "86c91f938024b2a25e2aeeece0532a36be16fe83ec9f2e41881d4a523a7b534a"
GATE5B_EXPECTED_TOO_SHALLOW = (
    "CHMP3",
    "COG3",
    "CPSF2",
    "EIF3F",
    "HAUS1",
    "PRPF38B",
    "RFFL",
    "RPF2",
    "TCP1",
    "ZNF676",
)
GATE5B_AUTHOR_EXP_PEARSON_DELTA = 0.372174


class Gate5BPanelCountError(RuntimeError):
    """Raised when Gate 5B fails Lock L1 exact count assertions (611 / 4,976 / 42,402)."""


@dataclass(frozen=True)
class BalancedHalfSplit:
    """One balanced batch-stratified split over panel row indices `0..n_panel_rows-1`."""

    split_index: int
    pert_A: np.ndarray
    pert_B: np.ndarray
    ctrl_A: np.ndarray
    ctrl_B: np.ndarray
    dropped: np.ndarray
    n_dropped_odd: int
    ctrl_batch_max_abs_diff: int
    pert_A_sha256: str
    pert_B_sha256: str
    ctrl_A_sha256: str
    ctrl_B_sha256: str


@dataclass
class Gate5BPanel:
    """Built Gate 5B HepG2 split-half panel (model-free)."""

    X_panel: Optional[np.ndarray]  # float32, shape (n_panel_rows, n_axis_genes) or None if metadata_only
    labels: np.ndarray  # str array, shape (n_panel_rows,)
    batches: np.ndarray  # str array, shape (n_panel_rows,)
    source_rows: np.ndarray  # int64 row indices into source h5ad
    axis_ensembl: List[str]
    cohort_perts: List[str]
    too_shallow: List[str]
    cells_per_pert: Dict[str, int]
    half_size: Dict[str, int]
    n_controls: int
    n_panel_rows: int
    cohort_sha256: str
    splits: List[BalancedHalfSplit]
    manifest_comparison: Dict[str, Any]
    experimental_metrics: Optional[Dict[str, Any]] = None


def _sha256_int64_array(arr: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(arr, dtype=np.int64).tobytes()).hexdigest()


def splithalf_group_seed(base_seed: int, split_index: int, group_name: str) -> int:
    """Deterministic independent 64-bit seed per `(base_seed, split_index, group_name)` via blake2b."""
    payload = f"{SPLITHALF_SALT}|{int(base_seed)}|{int(split_index)}|{group_name}".encode("utf-8")
    return int.from_bytes(hashlib.blake2b(payload, digest_size=8).digest(), "little")


def batch_stratified_equal_halves(
    batches: Sequence[str],
    rng: np.random.Generator,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Partition one group's local indices `0..len(batches)-1` into balanced batch-stratified halves
    `(idx_a, idx_b, idx_dropped)` by permuting unique batches, permuting within each batch,
    dealing alternately across batch boundaries, and dropping one random cell from the larger
    half if total count is odd.
    """
    b_arr = np.asarray(batches, dtype=str)
    n = len(b_arr)
    ordered: List[int] = []
    for batch_val in rng.permutation(np.unique(b_arr)):
        in_batch = np.flatnonzero(b_arr == batch_val)
        ordered.extend(in_batch[rng.permutation(len(in_batch))].tolist())

    half_a = np.asarray(ordered[0::2], dtype=np.int64)
    half_b = np.asarray(ordered[1::2], dtype=np.int64)
    dropped = np.zeros(0, dtype=np.int64)
    if n % 2 == 1:
        if len(half_a) > len(half_b):
            drop_pos = int(rng.integers(0, len(half_a)))
            dropped = np.asarray([half_a[drop_pos]], dtype=np.int64)
            half_a = np.delete(half_a, drop_pos)
        else:
            drop_pos = int(rng.integers(0, len(half_b)))
            dropped = np.asarray([half_b[drop_pos]], dtype=np.int64)
            half_b = np.delete(half_b, drop_pos)

    if len(half_a) != len(half_b):
        raise AssertionError(f"Unbalanced halves: {len(half_a)} vs {len(half_b)}")
    return half_a, half_b, dropped


def freeze_splithalf_cohort(
    labels: Sequence[str],
    eligible_test_perts: Sequence[str],
    *,
    min_cells: int = 10,
    control_label: str = SPLITHALF_CONTROL_LABEL,
) -> Dict[str, Any]:
    """Build deterministic frozen `min_cells` cohort from decoded perturbation labels."""
    lab_arr = np.asarray(labels, dtype=str)
    n_controls = int(np.sum(lab_arr == control_label))
    non_ctrl = lab_arr[lab_arr != control_label]
    counts: Dict[str, int] = {}
    if non_ctrl.size > 0:
        uniq_labs, uniq_cnts = np.unique(non_ctrl, return_counts=True)
        for u, c in zip(uniq_labs.tolist(), uniq_cnts.tolist()):
            counts[str(u)] = int(c)

    elig_set = set(str(x).strip() for x in eligible_test_perts)
    kept: List[str] = []
    too_shallow: List[str] = []
    for p in sorted(elig_set):
        cnt = counts.get(p, 0)
        if cnt >= int(min_cells):
            kept.append(p)
        else:
            too_shallow.append(p)

    return {
        "n_perts": len(kept),
        "perts": kept,
        "too_shallow": too_shallow,
        "cells_per_pert": {p: counts[p] for p in kept},
        "half_size": {p: counts[p] // 2 for p in kept},
        "n_controls": n_controls,
        "min_cells": int(min_cells),
        "cohort_sha256": compute_string_list_sha256(kept),
    }


def draw_splithalf_splits(
    panel_labels: Sequence[str],
    panel_batches: Sequence[str],
    cohort_perts: Sequence[str],
    *,
    n_splits: int = 5,
    base_seed: int = 42,
    control_label: str = SPLITHALF_CONTROL_LABEL,
) -> List[BalancedHalfSplit]:
    """Generate `n_splits` deterministic balanced batch-stratified splits on panel row indices."""
    lab_arr = np.asarray(panel_labels, dtype=str)
    bat_arr = np.asarray(panel_batches, dtype=str)
    perts = list(cohort_perts)

    ctrl_rows = np.flatnonzero(lab_arr == control_label)
    rows_by_pert = {p: np.flatnonzero(lab_arr == p) for p in perts}

    splits: List[BalancedHalfSplit] = []
    for s in range(int(n_splits)):
        pa_parts: List[np.ndarray] = []
        pb_parts: List[np.ndarray] = []
        drop_parts: List[np.ndarray] = []
        for p in perts:
            rows = rows_by_pert[p]
            rng = np.random.default_rng(splithalf_group_seed(base_seed, s, p))
            loc_a, loc_b, loc_d = batch_stratified_equal_halves(bat_arr[rows], rng)
            pa_parts.append(rows[loc_a])
            pb_parts.append(rows[loc_b])
            drop_parts.append(rows[loc_d])

        c_rng = np.random.default_rng(splithalf_group_seed(base_seed, s, SPLITHALF_CONTROL_GROUP))
        ca_loc, cb_loc, cd_loc = batch_stratified_equal_halves(bat_arr[ctrl_rows], c_rng)
        ctrl_A = ctrl_rows[ca_loc]
        ctrl_B = ctrl_rows[cb_loc]
        drop_parts.append(ctrl_rows[cd_loc])

        pert_A = np.concatenate(pa_parts) if pa_parts else np.zeros(0, dtype=np.int64)
        pert_B = np.concatenate(pb_parts) if pb_parts else np.zeros(0, dtype=np.int64)
        dropped = np.concatenate(drop_parts) if drop_parts else np.zeros(0, dtype=np.int64)

        # Disjointness & balance assertions
        if np.intersect1d(pert_A, pert_B).size > 0:
            raise ValueError(f"Split {s}: pert_A and pert_B overlap.")
        if np.intersect1d(ctrl_A, ctrl_B).size > 0:
            raise ValueError(f"Split {s}: ctrl_A and ctrl_B overlap.")
        if len(pert_A) != len(pert_B) or len(ctrl_A) != len(ctrl_B):
            raise ValueError(f"Split {s}: unequal half lengths.")

        # Control batch balance check (must be <= 1 per batch across all batches)
        ua, ca_cnt = np.unique(bat_arr[ctrl_A], return_counts=True)
        ub, cb_cnt = np.unique(bat_arr[ctrl_B], return_counts=True)
        map_a = dict(zip(ua.tolist(), ca_cnt.tolist()))
        map_b = dict(zip(ub.tolist(), cb_cnt.tolist()))
        all_b = set(map_a) | set(map_b)
        max_b_diff = max((abs(map_a.get(b, 0) - map_b.get(b, 0)) for b in all_b), default=0)

        splits.append(
            BalancedHalfSplit(
                split_index=s,
                pert_A=pert_A,
                pert_B=pert_B,
                ctrl_A=ctrl_A,
                ctrl_B=ctrl_B,
                dropped=dropped,
                n_dropped_odd=int(len(dropped)),
                ctrl_batch_max_abs_diff=int(max_b_diff),
                pert_A_sha256=_sha256_int64_array(pert_A),
                pert_B_sha256=_sha256_int64_array(pert_B),
                ctrl_A_sha256=_sha256_int64_array(ctrl_A),
                ctrl_B_sha256=_sha256_int64_array(ctrl_B),
            )
        )
    return splits


def compare_splits_against_manifest(
    splits: Sequence[BalancedHalfSplit],
    manifest: Optional[Mapping[str, Any]],
) -> Dict[str, Any]:
    """
    Compare generated split SHA-256 digests against `splits_manifest.json` (`§5.2` Case A vs Case B).
    """
    if manifest is None:
        return {
            "manifest_provided": False,
            "n_hash_matches": 0,
            "n_hash_total": 4 * len(splits),
            "tolerance_case": "NO_MANIFEST",
        }
    ref_by_idx = {int(rec["split"]): rec for rec in manifest.get("splits", [])}
    n_match = 0
    n_total = 0
    per_split: List[Dict[str, Any]] = []
    for sp_obj in splits:
        ref = ref_by_idx.get(sp_obj.split_index, {})
        flags = {}
        for key in ("pert_A_sha256", "pert_B_sha256", "ctrl_A_sha256", "ctrl_B_sha256"):
            n_total += 1
            ok = getattr(sp_obj, key) == ref.get(key)
            flags[key] = bool(ok)
            if ok:
                n_match += 1
        per_split.append(
            {
                "split": sp_obj.split_index,
                "n_pert_per_half": int(len(sp_obj.pert_A)),
                "n_ctrl_per_half": int(len(sp_obj.ctrl_A)),
                "n_dropped_odd": int(sp_obj.n_dropped_odd),
                "ctrl_batch_max_abs_diff": int(sp_obj.ctrl_batch_max_abs_diff),
                "hashes_match": flags,
                "counts_match_manifest": (
                    int(len(sp_obj.pert_A)) == int(ref.get("n_pert_per_half", -1))
                    and int(len(sp_obj.ctrl_A)) == int(ref.get("n_ctrl_per_half", -1))
                    and int(sp_obj.n_dropped_odd) == int(ref.get("n_dropped_odd", -1))
                ),
            }
        )

    case_a = n_total > 0 and n_match == n_total
    return {
        "manifest_provided": True,
        "n_hash_matches": int(n_match),
        "n_hash_total": int(n_total),
        "all_counts_match_manifest": all(r["counts_match_manifest"] for r in per_split),
        "max_ctrl_batch_abs_diff": max((r["ctrl_batch_max_abs_diff"] for r in per_split), default=0),
        "tolerance_case": "CASE_A_EXACT_SPLIT_MATCH" if case_a else "CASE_B_FRESH_DETERMINISTIC_SPLITS",
        "gate5b_tolerances": (
            {"pearson_delta": 0.005, "mae": 0.002, "experimental_pearson_delta": 0.002}
            if case_a
            else {"pearson_delta": 0.015, "mae": 0.005, "experimental_pearson_delta": 0.015}
        ),
        "per_split": per_split,
    }


def _stream_panel_rows_normalized_on_axis(
    x_node: Any,
    keep_rows: np.ndarray,
    axis_file_cols: np.ndarray,
    *,
    target_sum: float = OCCAM_TARGET_SUM,
    chunk_size: int = 4096,
) -> np.ndarray:
    """
    Stream strictly `keep_rows` from AnnData `/X` (Dense, CSR, or CSC), compute each cell's
    full-panel library size over ALL file genes (`float64`), restrict to `axis_file_cols`,
    and apply `log1p(x / L_i * target_sum)` in `float32`.
    """
    n_keep = int(len(keep_rows))
    n_axis = int(len(axis_file_cols))
    out = np.zeros((n_keep, n_axis), dtype=np.float32)
    if n_keep == 0:
        return out

    if isinstance(x_node, h5py.Dataset):
        # keep_rows is strictly increasing (`np.flatnonzero`)
        for s in range(0, n_keep, int(chunk_size)):
            e = min(s + int(chunk_size), n_keep)
            chunk_rows = keep_rows[s:e]
            block = np.asarray(x_node[chunk_rows, :], dtype=np.float32)
            full_rs = block.sum(axis=1, dtype=np.float64)
            axis_block = block[:, axis_file_cols]
            norm_block, _ = normalize_counts_target_sum_log1p(
                axis_block, target_sum=target_sum, full_row_sums=full_rs
            )
            out[s:e, :] = norm_block
        return out

    if isinstance(x_node, h5py.Group) and "indptr" in x_node and "indices" in x_node and "data" in x_node:
        enc = str(x_node.attrs.get("encoding-type", "")).lower()
        n_vars = int(x_node.attrs["shape"][1]) if "shape" in x_node.attrs else int(axis_file_cols.max()) + 1
        if "csc" in enc:
            shape = tuple(int(v) for v in x_node.attrs["shape"])
            mat_csr = sp.csc_matrix(
                (x_node["data"][:], x_node["indices"][:], x_node["indptr"][:]),
                shape=shape,
            ).tocsr()
            for s in range(0, n_keep, int(chunk_size)):
                e = min(s + int(chunk_size), n_keep)
                sub = mat_csr[keep_rows[s:e]]
                full_rs = np.asarray(sub.sum(axis=1), dtype=np.float64).ravel()
                axis_block = sub[:, axis_file_cols].toarray().astype(np.float32, copy=False)
                norm_block, _ = normalize_counts_target_sum_log1p(
                    axis_block, target_sum=target_sum, full_row_sums=full_rs
                )
                out[s:e, :] = norm_block
            return out

        # CSR sparse `/X`
        indptr = np.asarray(x_node["indptr"][:], dtype=np.int64)
        data_ds = x_node["data"]
        idx_ds = x_node["indices"]
        col_to_axis = np.full(n_vars, -1, dtype=np.int64)
        col_to_axis[axis_file_cols] = np.arange(n_axis, dtype=np.int64)

        for out_i, r_i in enumerate(keep_rows.tolist()):
            p0, p1 = int(indptr[r_i]), int(indptr[r_i + 1])
            if p1 <= p0:
                continue
            vals = np.asarray(data_ds[p0:p1], dtype=np.float32)
            cols = np.asarray(idx_ds[p0:p1], dtype=np.int64)
            rs = float(vals.sum(dtype=np.float64))
            if rs <= 0.0:
                continue
            tgt = col_to_axis[cols]
            keep = tgt >= 0
            if np.any(keep):
                out[out_i, tgt[keep]] = np.log1p(
                    (vals[keep] / np.float32(rs)) * np.float32(target_sum)
                )
        return out

    raise TypeError(f"Unsupported AnnData /X node type: {type(x_node).__name__}")


def compute_splithalf_experimental_metrics(
    X_panel: np.ndarray,
    panel_labels: Sequence[str],
    cohort_perts: Sequence[str],
    splits: Sequence[BalancedHalfSplit],
) -> Dict[str, Any]:
    """
    Model-free split-half experimental replicate agreement across `splits`:
    For each split `s` and perturbation `p` in `cohort_perts`:
      - Reference delta (Half A): `mean(X[pert_A == p]) - mean(X[ctrl_A])`
      - Experimental replicate delta (Half B): `mean(X[pert_B == p]) - mean(X[ctrl_B])`
    Computes mean per-perturbation Pearson correlation and MAE over the 6,983-gene axis
    (matching Occam's `nocal` `experimental` replicate convention; author ref = `0.372174`).
    """
    X = np.asarray(X_panel, dtype=np.float32)
    lab = np.asarray(panel_labels, dtype=str)
    perts = list(cohort_perts)

    split_pearsons: List[float] = []
    split_maes: List[float] = []
    for sp_obj in splits:
        ctrl_a_mean = X[sp_obj.ctrl_A].mean(axis=0, dtype=np.float64)
        ctrl_b_mean = X[sp_obj.ctrl_B].mean(axis=0, dtype=np.float64)
        lab_a = lab[sp_obj.pert_A]
        lab_b = lab[sp_obj.pert_B]

        # Group indices in O(N_pert_half) once per split
        rows_a_by_p: Dict[str, List[int]] = {}
        for row_i, p_name in zip(sp_obj.pert_A.tolist(), lab_a.tolist()):
            rows_a_by_p.setdefault(p_name, []).append(row_i)
        rows_b_by_p: Dict[str, List[int]] = {}
        for row_i, p_name in zip(sp_obj.pert_B.tolist(), lab_b.tolist()):
            rows_b_by_p.setdefault(p_name, []).append(row_i)

        rs_list: List[float] = []
        maes_list: List[float] = []
        for p in perts:
            ra = rows_a_by_p[p]
            rb = rows_b_by_p[p]
            da = X[ra].mean(axis=0, dtype=np.float64) - ctrl_a_mean
            db = X[rb].mean(axis=0, dtype=np.float64) - ctrl_b_mean
            da_c = da - da.mean()
            db_c = db - db.mean()
            denom = float(np.linalg.norm(da_c) * np.linalg.norm(db_c))
            if denom > 1e-12:
                rs_list.append(float(np.dot(da_c, db_c) / denom))
            maes_list.append(float(np.mean(np.abs(db - da))))

        split_pearsons.append(float(np.mean(rs_list)))
        split_maes.append(float(np.mean(maes_list)))

    mean_p = float(np.mean(split_pearsons))
    return {
        "n_splits": len(splits),
        "n_perts": len(perts),
        "experimental_pearson_delta_mean": mean_p,
        "experimental_pearson_delta_per_split": split_pearsons,
        "experimental_mae_mean": float(np.mean(split_maes)),
        "experimental_mae_per_split": split_maes,
        "author_reference_pearson_delta": GATE5B_AUTHOR_EXP_PEARSON_DELTA,
        "abs_diff_vs_author": float(abs(mean_p - GATE5B_AUTHOR_EXP_PEARSON_DELTA)),
    }


def build_gate5b_hepg2_panel(
    h5ad_path: Path,
    axis_ensembl: Sequence[str],
    eligible_test_perts: Sequence[str],
    *,
    splits_manifest: Optional[Mapping[str, Any]] = None,
    min_cells: int = 10,
    n_splits: int = 5,
    base_seed: int = 42,
    target_sum: float = OCCAM_TARGET_SUM,
    load_expression: bool = True,
    compute_experimental: bool = True,
    enforce_lock_l1_counts: bool = True,
    expected_n_perts: int = GATE5B_EXPECTED_N_COHORT_PERTS,
    expected_n_controls: int = GATE5B_EXPECTED_N_CONTROLS,
    expected_n_panel_rows: int = GATE5B_EXPECTED_N_PANEL_ROWS,
    expected_cohort_sha256: Optional[str] = GATE5B_EXPECTED_COHORT_SHA256,
    chunk_size: int = 4096,
) -> Gate5BPanel:
    """
    Build the model-free Gate 5B HepG2 split-half panel from local `nadig_hepg2` `.h5ad`.

    Lock `L1` enforcement:
      - Reads `obs["gene"]` (perturbation) and `var["gene_name"]` using `read_h5ad_table_col`
        with `forbid_raw_int_codes=True` (legacy `__categories` decoded; raw codes rejected).
      - When `enforce_lock_l1_counts=True`, halts with `Gate5BPanelCountError` unless exact
        counts `(611 perturbations, 4,976 controls, 42,402 panel rows)` and `cohort_sha256` match.
    """
    h5ad_path = Path(h5ad_path)
    axis = [strip_ensembl_version(g) for g in axis_ensembl]
    if len(set(axis)) != len(axis):
        raise ValueError("axis_ensembl contains duplicates.")

    with h5py.File(h5ad_path, "r") as f:
        obs = f["obs"]
        obs_keys = set(obs.keys())
        pert_col = next((c for c in ("gene", "gene_name", "target_gene", "perturbation") if c in obs_keys), None)
        if pert_col is None:
            raise KeyError(f"{h5ad_path}: obs missing perturbation column: {sorted(obs_keys)}")
        batch_col = next((c for c in ("gem_group", "dataset_batch_num", "batch") if c in obs_keys), None)
        if batch_col is None:
            raise KeyError(f"{h5ad_path}: obs missing batch column: {sorted(obs_keys)}")

        # Lock L1: strictly decode legacy AnnData categoricals; never treat int codes as labels
        all_labels = [s.strip() for s in read_h5ad_table_col(obs, pert_col, forbid_raw_int_codes=True)]
        all_batches = [s.strip() for s in read_h5ad_table_col(obs, batch_col)]

        var = f["var"]
        var_keys = set(var.keys())
        ens_col = next(
            (c for c in ("gene_id", "ensembl_gene_id", "ensembl_id") if c in var_keys),
            None,
        )
        if ens_col is None:
            raise KeyError(f"{h5ad_path}: var missing Ensembl ID column: {sorted(var_keys)}")
        var_ens = [strip_ensembl_version(x) for x in read_h5ad_table_col(var, ens_col, forbid_raw_int_codes=True)]
        if "gene_name" in var_keys:
            # Exercise Lock L1 legacy categorical verification on var["gene_name"] as well
            _ = read_h5ad_table_col(var, "gene_name", forbid_raw_int_codes=True)

        var_pos: Dict[str, int] = {}
        for i, ens in enumerate(var_ens):
            var_pos.setdefault(ens, i)
        missing_axis = [g for g in axis if g not in var_pos]
        if missing_axis:
            raise ValueError(
                f"{len(missing_axis)} RHJ axis gene(s) absent from {h5ad_path} (e.g. {missing_axis[:5]})."
            )
        axis_file_cols = np.asarray([var_pos[g] for g in axis], dtype=np.int64)

        cohort = freeze_splithalf_cohort(
            all_labels, eligible_test_perts, min_cells=min_cells, control_label=SPLITHALF_CONTROL_LABEL
        )
        labels_arr = np.asarray(all_labels, dtype=str)
        batches_arr = np.asarray(all_batches, dtype=str)
        keep_set = set(cohort["perts"]) | {SPLITHALF_CONTROL_LABEL}
        keep_rows = np.flatnonzero(np.isin(labels_arr, sorted(keep_set))).astype(np.int64)
        panel_labels = labels_arr[keep_rows]
        panel_batches = batches_arr[keep_rows]

        n_perts_found = int(cohort["n_perts"])
        n_ctrl_found = int(np.sum(panel_labels == SPLITHALF_CONTROL_LABEL))
        n_panel_found = int(len(keep_rows))

        if enforce_lock_l1_counts:
            errors: List[str] = []
            if n_perts_found != expected_n_perts:
                errors.append(f"n_cohort_perts={n_perts_found} != expected {expected_n_perts}")
            if n_ctrl_found != expected_n_controls:
                errors.append(f"n_controls={n_ctrl_found} != expected {expected_n_controls}")
            if n_panel_found != expected_n_panel_rows:
                errors.append(f"n_panel_rows={n_panel_found} != expected {expected_n_panel_rows}")
            if expected_cohort_sha256 is not None and cohort["cohort_sha256"] != expected_cohort_sha256:
                errors.append(
                    f"cohort_sha256={cohort['cohort_sha256']} != expected {expected_cohort_sha256}"
                )
            if errors:
                raise Gate5BPanelCountError(
                    "Gate 5B Lock L1 count verification FAILED (halting): " + "; ".join(errors)
                )

        X_panel: Optional[np.ndarray] = None
        if load_expression:
            X_panel = _stream_panel_rows_normalized_on_axis(
                f["X"],
                keep_rows,
                axis_file_cols,
                target_sum=target_sum,
                chunk_size=chunk_size,
            )
            if not np.all(np.isfinite(X_panel)):
                raise ValueError("Built Gate 5B X_panel contains non-finite values.")

    splits = draw_splithalf_splits(
        panel_labels,
        panel_batches,
        cohort["perts"],
        n_splits=n_splits,
        base_seed=base_seed,
        control_label=SPLITHALF_CONTROL_LABEL,
    )
    if any(sp_obj.ctrl_batch_max_abs_diff > 1 for sp_obj in splits):
        bad = [(sp_obj.split_index, sp_obj.ctrl_batch_max_abs_diff) for sp_obj in splits if sp_obj.ctrl_batch_max_abs_diff > 1]
        raise ValueError(f"Gate 5B control batch stratification violated (max_abs_diff > 1): {bad}")

    manifest_cmp = compare_splits_against_manifest(splits, splits_manifest)
    exp_metrics = (
        compute_splithalf_experimental_metrics(X_panel, panel_labels, cohort["perts"], splits)
        if (load_expression and compute_experimental and X_panel is not None)
        else None
    )

    return Gate5BPanel(
        X_panel=X_panel,
        labels=panel_labels,
        batches=panel_batches,
        source_rows=keep_rows,
        axis_ensembl=axis,
        cohort_perts=list(cohort["perts"]),
        too_shallow=list(cohort["too_shallow"]),
        cells_per_pert=dict(cohort["cells_per_pert"]),
        half_size=dict(cohort["half_size"]),
        n_controls=n_ctrl_found,
        n_panel_rows=n_panel_found,
        cohort_sha256=str(cohort["cohort_sha256"]),
        splits=splits,
        manifest_comparison=manifest_cmp,
        experimental_metrics=exp_metrics,
    )


def main() -> int:
    ap = argparse.ArgumentParser(description="Model-free Gate 5B HepG2 panel count & split check.")
    ap.add_argument(
        "--h5ad",
        type=Path,
        default=DEFAULT_LOCAL_NADIG_DIR / "GSE264667_hepg2_raw_singlecell_01.h5ad",
    )
    ap.add_argument("--gene-axis-csv", type=Path, required=True)
    ap.add_argument("--perturbations-csv", type=Path, required=True)
    ap.add_argument("--splits-manifest-json", type=Path, default=None)
    ap.add_argument("--metadata-only", action="store_true", help="Verify counts & splits without reading /X.")
    ap.add_argument("--out-json", type=Path, default=None)
    args = ap.parse_args()

    axis = load_rhj_gene_axis(args.gene_axis_csv, verify_sha256=True)
    ptable = load_rhj_perturbations_table(args.perturbations_csv)
    manifest = (
        json.loads(args.splits_manifest_json.read_text(encoding="utf-8"))
        if args.splits_manifest_json is not None and args.splits_manifest_json.is_file()
        else None
    )

    panel = build_gate5b_hepg2_panel(
        h5ad_path=args.h5ad,
        axis_ensembl=axis.ensembl_ids,
        eligible_test_perts=ptable["test_keys"],
        splits_manifest=manifest,
        load_expression=(not args.metadata_only),
        compute_experimental=(not args.metadata_only),
        enforce_lock_l1_counts=True,
    )

    summary = {
        "h5ad_path": str(args.h5ad),
        "n_axis_genes": len(panel.axis_ensembl),
        "gene_axis_sha256": axis.axis_sha256,
        "n_eligible_test_perts": len(ptable["test_keys"]),
        "n_cohort_perts": len(panel.cohort_perts),
        "n_controls": panel.n_controls,
        "n_panel_rows": panel.n_panel_rows,
        "n_perturbed_rows": int(panel.n_panel_rows - panel.n_controls),
        "cohort_sha256": panel.cohort_sha256,
        "too_shallow": panel.too_shallow,
        "manifest_comparison": panel.manifest_comparison,
        "experimental_metrics": panel.experimental_metrics,
    }
    text = json.dumps(summary, indent=2, sort_keys=True)
    print(text)
    if args.out_json is not None:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(text + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
