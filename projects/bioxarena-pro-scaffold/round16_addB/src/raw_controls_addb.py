#!/usr/bin/env python3
"""
Model-free raw control loaders for Round 16 Addendum B (`Occam`), loading exact
cells by the locked Round 16 barcode lists (`PREREG_I_ADDENDUM_B_DRAFT_v2.md §3.2–§3.3`,
`PREREG_I_ADDENDUM_B_LOCK.md L1`).

Design invariants:
  - Strictly read-only on all raw datasets and `round16/results/` JSON files.
  - Zero Occam imports or vendored Occam code.
  - Implements `DEVIATION_R16_D3` / Lock `L1` legacy AnnData categorical decoding
    (`codes` dataset + sibling `__categories/<col>`), so integer codes are NEVER
    returned as perturbation or gene labels.
  - Implements `DEVIATION_R16_D1` lane gene-order alignment across multi-lane 10x `.h5`
    files (`cd4t_rest`, `dld1`, optional `gse295214_h1`).
  - Preserves exact locked barcode order (`0..N-1`) and extracts version-stripped
    Ensembl gene IDs alongside gene symbols for Occam's `44,728`-gene Ensembl space.
"""

from __future__ import annotations

import csv
import gzip
import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import h5py
import numpy as np
import scipy.sparse as sp

ADDB_SRC_DIR = Path(__file__).resolve().parent
ADDB_ROOT = ADDB_SRC_DIR.parent
REPO_ROOT = ADDB_ROOT.parent

DEFAULT_PRIMARY_BARCODES_JSON = REPO_ROOT / "round16" / "results" / "control_barcodes_r16.json"
DEFAULT_SENSITIVITY_BARCODES_JSON = (
    REPO_ROOT / "round16" / "results" / "sensitivity" / "control_barcodes_r16.json"
)
DEFAULT_LOCAL_NADIG_DIR = Path("/Users/mac/.cache/r9_p2_raw/nadig2025")

LOCKED_PRIMARY_SEED = 20261601
LOCKED_RESAMPLE_SEED = 20261602


@dataclass
class RawControlPool:
    """Standardized model-free raw UMI control matrix loaded by locked R16 barcodes."""

    line_name: str
    pool_kind: str  # "primary" or "resample"
    seed: int
    barcodes: List[str]
    raw_counts: np.ndarray  # float32, shape (n_cells, n_raw_genes)
    raw_ensembl_ids: List[str]  # version-stripped Ensembl IDs (len == n_raw_genes)
    raw_gene_symbols: List[str]  # gene symbols (len == n_raw_genes)
    source_files: List[str]
    audit: Dict[str, Any] = field(default_factory=dict)


def decode_hdf5_str_array(arr: Any) -> List[str]:
    """Decode an HDF5 byte/unicode/numeric 1D array cleanly into Python strings."""
    out: List[str] = []
    for val in np.asarray(arr).tolist():
        if isinstance(val, (bytes, np.bytes_)):
            out.append(val.decode("utf-8", errors="replace"))
        else:
            out.append(str(val))
    return out


def strip_ensembl_version(ens_id: str) -> str:
    """Strip Ensembl `.N` version suffix (e.g., `ENSG00000187634.13` -> `ENSG00000187634`)."""
    s = str(ens_id).strip()
    if not s:
        return ""
    return s.split(".", 1)[0]


def read_h5ad_table_col(
    table_grp: h5py.Group,
    col: str,
    *,
    forbid_raw_int_codes: bool = False,
) -> List[str]:
    """
    Read an AnnData `.obs` or `.var` column from an open `h5py.Group`.

    Supports all three AnnData storage encodings:
      1. Modern categorical (`h5py.Group` with `codes` and `categories`);
      2. Legacy AnnData categorical (`DEVIATION_R16_D3` / Lock `L1`: integer `Dataset`
         with sibling `__categories/<col>`);
      3. Plain `h5py.Dataset` (strings, bytes, or numbers).

    If `forbid_raw_int_codes=True`, raises `ValueError` when the returned column
    still comes from an undecoded integer dataset (guarding Lock `L1`).
    """
    if col not in table_grp:
        raise KeyError(f"Column {col!r} not found in HDF5 group {table_grp.name!r}: {sorted(table_grp.keys())}")
    node = table_grp[col]

    # 1. Modern AnnData categorical: group containing `categories` + `codes`
    if isinstance(node, h5py.Group) and "codes" in node and "categories" in node:
        cats = decode_hdf5_str_array(node["categories"][:])
        codes = np.asarray(node["codes"][:], dtype=np.int64)
        n_cats = len(cats)
        return [cats[c] if 0 <= c < n_cats else "" for c in codes]

    # 2. Legacy AnnData categorical (DEVIATION_R16_D3 / Lock L1): integer dataset + `__categories/<col>`
    if (
        isinstance(node, h5py.Dataset)
        and "__categories" in table_grp
        and isinstance(table_grp["__categories"], h5py.Group)
        and col in table_grp["__categories"]
        and np.issubdtype(node.dtype, np.integer)
    ):
        cats = decode_hdf5_str_array(table_grp["__categories"][col][:])
        codes = np.asarray(node[:], dtype=np.int64)
        n_cats = len(cats)
        return [cats[c] if 0 <= c < n_cats else "" for c in codes]

    # 3. Plain dataset
    if isinstance(node, h5py.Dataset):
        if forbid_raw_int_codes and np.issubdtype(node.dtype, np.integer):
            raise ValueError(
                f"Column {col!r} in {table_grp.name!r} is an undecoded integer dataset "
                f"(dtype={node.dtype}) with no `__categories/{col}` table; refusing to use "
                f"integer codes as labels (Lock L1)."
            )
        return decode_hdf5_str_array(node[:])

    raise TypeError(f"Unsupported HDF5 node type for {table_grp.name}/{col}: {type(node).__name__}")


def resolve_h5ad_index_col(table_grp: h5py.Group) -> Optional[str]:
    """Resolve AnnData `.obs`/`.var` index dataset name (`attrs['_index']` first, per D4)."""
    keys = set(table_grp.keys())
    attr_idx = table_grp.attrs.get("_index")
    if attr_idx is not None:
        raw_val = (
            attr_idx.reshape(-1)[0]
            if (isinstance(attr_idx, np.ndarray) and attr_idx.size >= 1)
            else attr_idx
        )
        attr_s = (
            raw_val.decode("utf-8", errors="replace").strip()
            if isinstance(raw_val, (bytes, np.bytes_))
            else str(raw_val).strip()
        )
        if attr_s in keys:
            return attr_s
    for cand in ("_index", "index", "cell_barcode", "barcode", "cell_id", "gene_id"):
        if cand in keys:
            return cand
    return None


def read_h5ad_x_rows(x_node: Any, row_indices: Sequence[int], n_genes: int) -> np.ndarray:
    """
    Read exact ordered `row_indices` from an AnnData `/X` node (Dense, CSR, or CSC)
    into a contiguous `float32` array of shape `(len(row_indices), n_genes)`.
    """
    sel = np.asarray(row_indices, dtype=np.int64)
    if sel.size == 0:
        return np.zeros((0, int(n_genes)), dtype=np.float32)

    if isinstance(x_node, h5py.Dataset):
        # h5py fancy indexing requires strictly increasing unique indices
        uniq_rows, inv = np.unique(sel, return_inverse=True)
        uniq_mat = np.asarray(x_node[uniq_rows, :], dtype=np.float32)
        return np.ascontiguousarray(uniq_mat[inv], dtype=np.float32)

    if isinstance(x_node, h5py.Group) and "indptr" in x_node and "indices" in x_node and "data" in x_node:
        enc = str(x_node.attrs.get("encoding-type", "")).lower()
        if "csc" in enc:
            shape = tuple(int(v) for v in x_node.attrs["shape"])
            mat_csr = sp.csc_matrix(
                (x_node["data"][:], x_node["indices"][:], x_node["indptr"][:]),
                shape=shape,
            ).tocsr()
            return np.ascontiguousarray(mat_csr[sel].toarray(), dtype=np.float32)
        indptr = np.asarray(x_node["indptr"][:], dtype=np.int64)
        data_ds = x_node["data"]
        idx_ds = x_node["indices"]
        out = np.zeros((sel.size, int(n_genes)), dtype=np.float32)
        for out_i, r_i in enumerate(sel.tolist()):
            p0, p1 = int(indptr[r_i]), int(indptr[r_i + 1])
            if p1 > p0:
                out[out_i, np.asarray(idx_ds[p0:p1], dtype=np.int64)] = np.asarray(
                    data_ds[p0:p1], dtype=np.float32
                )
        return np.ascontiguousarray(out, dtype=np.float32)

    raise TypeError(f"Unsupported AnnData /X node type: {type(x_node).__name__}")


def load_locked_r16_barcodes(
    line_name: str,
    pool: str = "primary",
    *,
    primary_json: Path = DEFAULT_PRIMARY_BARCODES_JSON,
    sensitivity_json: Path = DEFAULT_SENSITIVITY_BARCODES_JSON,
) -> Tuple[List[str], int, Dict[str, Any]]:
    """
    Load exact ordered locked Round 16 control barcodes (`PREREG_I_ADDENDUM_B_DRAFT_v2.md §3.2`).

    - `pool == "primary"`: loads `seed = 20261601` ordered barcodes from `primary_json`
      (falling back to `sensitivity_json` only if `line_name` is absent from `primary_json`).
    - `pool == "resample"`: loads `seed = 20261602` disjoint within-line resample barcodes,
      using `sensitivity_json` for `feng_ipsc` (where `primary_json` has `eligible=False`)
      and `primary_json` for Group 1–2 lines, asserting `overlap_with_primary == 0`.
    """
    if pool not in ("primary", "resample"):
        raise ValueError(f"pool must be 'primary' or 'resample', got {pool!r}")

    prim_data = json.loads(Path(primary_json).read_text(encoding="utf-8")) if Path(primary_json).is_file() else {}
    sens_data = (
        json.loads(Path(sensitivity_json).read_text(encoding="utf-8"))
        if Path(sensitivity_json).is_file()
        else {}
    )

    if pool == "primary":
        entry = prim_data.get(line_name) or sens_data.get(line_name)
        src_path = str(primary_json if line_name in prim_data else sensitivity_json)
        if entry is None:
            raise KeyError(f"Line {line_name!r} not found in locked R16 barcode JSONs.")
        bcs = [str(b) for b in entry.get("barcodes", [])]
        seed = int(entry.get("primary_seed", LOCKED_PRIMARY_SEED))
        if not bcs:
            raise ValueError(f"Empty primary barcode list for {line_name!r} in {src_path}.")
        if len(set(bcs)) != len(bcs):
            raise ValueError(f"Duplicate primary barcodes found for {line_name!r} in {src_path}.")
        return bcs, seed, {
            "source_json": src_path,
            "n_ctrl_total": int(entry.get("n_ctrl_total", 0)),
            "n_ctrl_sampled": len(bcs),
            "pool": "primary",
        }

    # Resample pool (`seed = 20261602`): prefer whichever JSON has eligible non-empty resample barcodes
    chosen_entry: Optional[Dict[str, Any]] = None
    src_path = ""
    for cand_path, cand_data in ((primary_json, prim_data), (sensitivity_json, sens_data)):
        e = cand_data.get(line_name)
        if not e:
            continue
        w = e.get("within_line_resample") or {}
        if w.get("eligible") and w.get("barcodes"):
            chosen_entry = e
            src_path = str(cand_path)
            break
    if chosen_entry is None:
        raise KeyError(f"No eligible within_line_resample barcodes found for {line_name!r}.")

    prim_bcs = [str(b) for b in chosen_entry.get("barcodes", [])]
    w = chosen_entry["within_line_resample"]
    resamp_bcs = [str(b) for b in w["barcodes"]]
    seed = int(w.get("seed", LOCKED_RESAMPLE_SEED))
    if len(set(resamp_bcs)) != len(resamp_bcs):
        raise ValueError(f"Duplicate resample barcodes found for {line_name!r} in {src_path}.")
    overlap = len(set(prim_bcs) & set(resamp_bcs))
    if overlap != 0:
        raise ValueError(
            f"Resample barcodes for {line_name!r} overlap primary by {overlap} cells (must be 0)."
        )
    return resamp_bcs, seed, {
        "source_json": src_path,
        "n_ctrl_total": int(chosen_entry.get("n_ctrl_total", 0)),
        "n_resample": len(resamp_bcs),
        "overlap_with_primary": 0,
        "pool": "resample",
    }


def _validate_requested_barcodes(barcodes: Sequence[str], context: str) -> List[str]:
    req = [str(b).strip() for b in barcodes]
    if not req:
        raise ValueError(f"[{context}] Requested barcode list is empty.")
    if len(set(req)) != len(req):
        raise ValueError(f"[{context}] Requested barcode list contains duplicates.")
    return req


def load_nadig_controls_by_barcodes(
    h5ad_path: Path,
    barcodes: Sequence[str],
    *,
    line_name: str,
    pool_kind: str = "primary",
    seed: int = LOCKED_PRIMARY_SEED,
    verify_ntc: bool = True,
) -> RawControlPool:
    """
    Load exact raw counts for locked barcodes from a Nadig 2025 `.h5ad` file
    (local `/Users/mac/.cache/r9_p2_raw/nadig2025/GSE264667_{hepg2,jurkat}_raw_singlecell_01.h5ad`),
    decoding legacy AnnData categoricals per `DEVIATION_R16_D3` / Lock `L1`.
    """
    req = _validate_requested_barcodes(barcodes, line_name)
    h5ad_path = Path(h5ad_path)
    with h5py.File(h5ad_path, "r") as f:
        obs = f["obs"]
        obs_idx_col = resolve_h5ad_index_col(obs)
        if obs_idx_col is None:
            raise KeyError(f"Nadig h5ad {h5ad_path} missing obs index column.")
        obs_bcs = [s.strip() for s in read_h5ad_table_col(obs, obs_idx_col)]
        bc_to_row: Dict[str, int] = {}
        for idx, bc in enumerate(obs_bcs):
            if bc in bc_to_row:
                raise ValueError(f"Duplicate obs barcode {bc!r} in {h5ad_path}.")
            bc_to_row[bc] = idx

        missing = [bc for bc in req if bc not in bc_to_row]
        if missing:
            raise KeyError(
                f"[{line_name}] {len(missing)}/{len(req)} locked barcode(s) absent from {h5ad_path} "
                f"(e.g. {missing[:5]})."
            )
        sel_rows = np.asarray([bc_to_row[bc] for bc in req], dtype=np.int64)

        obs_keys = set(obs.keys())
        gene_col = next((c for c in ("gene", "target_gene", "perturbation") if c in obs_keys), None)
        if verify_ntc:
            if gene_col is None:
                raise KeyError(f"Nadig h5ad {h5ad_path} missing obs perturbation column for NTC verification.")
            pert_vals = read_h5ad_table_col(obs, gene_col, forbid_raw_int_codes=True)
            bad_ntc = [
                (req[i], pert_vals[int(r)])
                for i, r in enumerate(sel_rows)
                if str(pert_vals[int(r)]).strip().lower() not in ("non-targeting", "nontargeting")
            ]
            if bad_ntc:
                raise ValueError(
                    f"[{line_name}] {len(bad_ntc)} locked barcode(s) are NOT non-targeting in {h5ad_path}: "
                    f"{bad_ntc[:5]}"
                )

        var = f["var"]
        var_keys = set(var.keys())
        ens_col = next(
            (c for c in ("gene_id", "ensembl_gene_id", "ensembl_id", "gene_ids") if c in var_keys),
            None,
        )
        if ens_col is None:
            ens_col = resolve_h5ad_index_col(var)
        if ens_col is None:
            raise KeyError(f"Nadig h5ad {h5ad_path} missing var Ensembl ID column: {sorted(var_keys)}")
        raw_ens = [strip_ensembl_version(x) for x in read_h5ad_table_col(var, ens_col, forbid_raw_int_codes=True)]

        sym_col = next(
            (c for c in ("gene_name", "gene_symbol", "feature_name", "symbol") if c in var_keys),
            None,
        )
        raw_syms = (
            [s.strip() for s in read_h5ad_table_col(var, sym_col, forbid_raw_int_codes=True)]
            if sym_col is not None
            else list(raw_ens)
        )

        raw_counts = read_h5ad_x_rows(f["X"], sel_rows, len(raw_ens))

    return RawControlPool(
        line_name=line_name,
        pool_kind=pool_kind,
        seed=int(seed),
        barcodes=list(req),
        raw_counts=raw_counts,
        raw_ensembl_ids=raw_ens,
        raw_gene_symbols=raw_syms,
        source_files=[str(h5ad_path)],
        audit={
            "loader": "load_nadig_controls_by_barcodes",
            "n_cells": int(raw_counts.shape[0]),
            "n_raw_genes": int(raw_counts.shape[1]),
            "verified_ntc": bool(verify_ntc),
            "obs_index_col": obs_idx_col,
            "var_ensembl_col": ens_col,
            "var_symbol_col": sym_col,
        },
    )


def _read_10x_h5_gex_meta(h5_path: Path) -> Tuple[List[str], List[str], np.ndarray, List[str]]:
    """Read `(gex_ensembl_ids, gex_gene_symbols, gex_row_indices, raw_barcodes)` from a 10x `.h5`."""
    with h5py.File(h5_path, "r") as f:
        m = f["matrix"]
        raw_bcs = [s.strip() for s in decode_hdf5_str_array(m["barcodes"][:])]
        fids = [strip_ensembl_version(x) for x in decode_hdf5_str_array(m["features/id"][:])]
        fnames = [s.strip() for s in decode_hdf5_str_array(m["features/name"][:])]
        if "feature_type" in m["features"]:
            ftypes = [s.strip() for s in decode_hdf5_str_array(m["features/feature_type"][:])]
            gex_rows = np.asarray([i for i, t in enumerate(ftypes) if t == "Gene Expression"], dtype=np.int64)
        else:
            gex_rows = np.arange(len(fids), dtype=np.int64)
    gex_ens = [fids[int(i)] for i in gex_rows]
    gex_sym = [fnames[int(i)] for i in gex_rows]
    return gex_ens, gex_sym, gex_rows, raw_bcs


def _read_10x_h5_selected_columns(
    h5_path: Path,
    col_indices: Sequence[int],
    gex_rows: np.ndarray,
    n_total_features: int,
) -> np.ndarray:
    """Read exact `col_indices` of a 10x CSC matrix restricted to `gex_rows`, returning `(n_sel, len(gex_rows))`."""
    n_sel = len(col_indices)
    n_gex = len(gex_rows)
    out = np.zeros((n_sel, n_gex), dtype=np.float32)
    if n_sel == 0:
        return out
    feat_to_gex = np.full(int(n_total_features), -1, dtype=np.int32)
    feat_to_gex[gex_rows] = np.arange(n_gex, dtype=np.int32)

    with h5py.File(h5_path, "r") as f:
        m = f["matrix"]
        indptr = np.asarray(m["indptr"][:], dtype=np.int64)
        data_ds = m["data"]
        idx_ds = m["indices"]
        for out_i, col_j in enumerate(col_indices):
            p0, p1 = int(indptr[int(col_j)]), int(indptr[int(col_j) + 1])
            if p1 <= p0:
                continue
            r_idx = np.asarray(idx_ds[p0:p1], dtype=np.int64)
            r_val = np.asarray(data_ds[p0:p1], dtype=np.float32)
            mapped = feat_to_gex[r_idx]
            keep = mapped >= 0
            if np.any(keep):
                out[out_i, mapped[keep]] = r_val[keep]
    return out


def load_10x_multi_h5_controls_by_barcodes(
    prefix_to_h5: Dict[str, Path],
    barcodes: Sequence[str],
    *,
    line_name: str,
    canonical_ref_prefix: Optional[str] = None,
    pool_kind: str = "primary",
    seed: int = LOCKED_PRIMARY_SEED,
    allow_duplicate_ensembl: bool = False,
) -> RawControlPool:
    """
    Load exact raw counts from multi-lane 10x `.h5` files (`dld1`, `cd4t_rest`, or `gse295214_h1`)
    by locked `f"{prefix}:{raw_barcode}"` identifiers.

    Implements `DEVIATION_R16_D1`: if a lane contains the exact same unique Ensembl gene set as
    the canonical reference lane (`canonical_ref_prefix`, e.g., `CD4i_R2L01` for `cd4t_rest`)
    in a different column order, its columns are permuted onto the canonical reference order and
    the lane prefix is recorded in `audit["reordered_lane_prefixes"]`. Any gene-set mismatch raises.
    """
    req = _validate_requested_barcodes(barcodes, line_name)
    if not prefix_to_h5:
        raise ValueError(f"[{line_name}] Empty prefix_to_h5 map.")

    ref_prefix = canonical_ref_prefix or sorted(prefix_to_h5.keys())[0]
    if ref_prefix not in prefix_to_h5:
        raise KeyError(f"[{line_name}] canonical_ref_prefix {ref_prefix!r} not in prefix_to_h5.")

    ref_ens, ref_sym, _, _ = _read_10x_h5_gex_meta(prefix_to_h5[ref_prefix])
    n_genes = len(ref_ens)
    ref_unique = len(set(ref_ens)) == n_genes
    if not allow_duplicate_ensembl and not ref_unique:
        raise ValueError(
            f"[{line_name}] Reference lane {ref_prefix!r} contains duplicate Ensembl IDs "
            f"({n_genes - len(set(ref_ens))} duplicates)."
        )

    # Group requested barcodes by prefix while retaining exact target output slot `out_idx`
    by_prefix: Dict[str, List[Tuple[int, str]]] = {}
    for out_idx, full_bc in enumerate(req):
        if ":" not in full_bc:
            raise ValueError(f"[{line_name}] Expected 'prefix:barcode' format, got {full_bc!r}.")
        pref, raw_bc = full_bc.split(":", 1)
        if pref not in prefix_to_h5:
            raise KeyError(f"[{line_name}] Barcode prefix {pref!r} not found in available lanes.")
        by_prefix.setdefault(pref, []).append((out_idx, raw_bc))

    raw_counts = np.zeros((len(req), n_genes), dtype=np.float32)
    reordered_lanes: List[str] = []
    used_files: List[str] = []

    for pref in sorted(by_prefix.keys()):
        h5_path = prefix_to_h5[pref]
        used_files.append(str(h5_path))
        lane_ens, _, gex_rows, lane_bcs = _read_10x_h5_gex_meta(h5_path)
        with h5py.File(h5_path, "r") as f:
            n_total_features = int(f["matrix/shape"][0])

        # D1 gene-order check and permutation onto canonical reference axis
        perm_to_ref: Optional[np.ndarray] = None
        if lane_ens != ref_ens:
            if not ref_unique or len(lane_ens) != n_genes or len(set(lane_ens)) != n_genes:
                raise ValueError(
                    f"[{line_name}] Lane {pref!r} Ensembl gene axis cannot be aligned uniquely to "
                    f"reference {ref_prefix!r} (len {len(lane_ens)} vs {n_genes})."
                )
            if set(lane_ens) != set(ref_ens):
                raise ValueError(
                    f"[{line_name}] Lane {pref!r} has a different Ensembl gene set than reference {ref_prefix!r}."
                )
            lane_pos = {g: i for i, g in enumerate(lane_ens)}
            perm_to_ref = np.asarray([lane_pos[g] for g in ref_ens], dtype=np.int64)
            reordered_lanes.append(pref)

        bc_map = {bc: j for j, bc in enumerate(lane_bcs)}
        req_items = by_prefix[pref]
        missing = [rb for _, rb in req_items if rb not in bc_map]
        if missing:
            raise KeyError(
                f"[{line_name}] {len(missing)} barcode(s) absent from {h5_path} (e.g. {missing[:5]})."
            )
        col_indices = [bc_map[rb] for _, rb in req_items]
        lane_mat = _read_10x_h5_selected_columns(h5_path, col_indices, gex_rows, n_total_features)
        if perm_to_ref is not None:
            lane_mat = lane_mat[:, perm_to_ref]

        out_rows = np.asarray([oi for oi, _ in req_items], dtype=np.int64)
        raw_counts[out_rows, :] = lane_mat

    return RawControlPool(
        line_name=line_name,
        pool_kind=pool_kind,
        seed=int(seed),
        barcodes=list(req),
        raw_counts=np.ascontiguousarray(raw_counts, dtype=np.float32),
        raw_ensembl_ids=ref_ens,
        raw_gene_symbols=ref_sym,
        source_files=used_files,
        audit={
            "loader": "load_10x_multi_h5_controls_by_barcodes",
            "canonical_ref_prefix": ref_prefix,
            "n_lanes_opened": len(used_files),
            "reordered_lane_prefixes": reordered_lanes,
            "n_cells": int(raw_counts.shape[0]),
            "n_raw_genes": int(raw_counts.shape[1]),
        },
    )


def build_dld1_prefix_map(raw_dir: Path) -> Tuple[Dict[str, Path], str]:
    """Map `GSM98642xx` -> filtered 10x `.h5` path for `dld1`."""
    ddir = Path(raw_dir) / "dld1" if (Path(raw_dir) / "dld1").is_dir() else Path(raw_dir)
    rna_files = sorted(ddir.glob("GSM98642*.filtered_feature_bc_matrix.h5"))
    if not rna_files:
        raise FileNotFoundError(f"No DLD-1 RNA .h5 files found in {ddir}")
    pmap = {p.name.split("_")[0]: p for p in rna_files}
    return pmap, rna_files[0].name.split("_")[0]


def build_cd4t_rest_prefix_map(raw_dir: Path) -> Tuple[Dict[str, Path], str]:
    """Map `f'{donor_id}_{lane_tag}'` -> filtered 10x `.h5` path for `cd4t_rest` (canonical ref = D4[0])."""
    cdir = Path(raw_dir) / "cd4t_rest" if (Path(raw_dir) / "cd4t_rest").is_dir() else Path(raw_dir)
    d4_files = sorted(cdir.glob("GSM93938*_D4_Rest_sample_filtered_feature_bc_matrix.h5"))
    d1_files = sorted(cdir.glob("GSM93938*_D1_Rest_sample_filtered_feature_bc_matrix.h5"))
    if not d4_files or not d1_files:
        raise FileNotFoundError(f"Missing cd4t_rest D4/D1 .h5 files in {cdir}")
    pmap: Dict[str, Path] = {}
    canonical_ref: Optional[str] = None
    for donor_id, files in (("CE0006864_D4", d4_files), ("CE0008162_D1", d1_files)):
        for lp in files:
            toks = lp.name.split("_")
            lane_tag = f"{toks[1]}_{toks[2]}"
            pref = f"{donor_id}_{lane_tag}"
            pmap[pref] = lp
            if canonical_ref is None:
                canonical_ref = pref
    assert canonical_ref is not None
    return pmap, canonical_ref


def _find_matching_file(directory: Path, preferred: Sequence[str], glob_pattern: str) -> Path:
    for name in preferred:
        cand = directory / name
        if cand.is_file():
            return cand
    matches = sorted(directory.glob(glob_pattern))
    if not matches:
        raise FileNotFoundError(f"No file matching {glob_pattern!r} in {directory}")
    return matches[0]


def _pick_col_ci(columns: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
    low_map = {str(c).lower(): str(c) for c in columns}
    for cand in candidates:
        if cand in columns:
            return cand
        if cand.lower() in low_map:
            return low_map[cand.lower()]
    return None


def _read_mtx_header(mtx_path: Path) -> Tuple[int, int, int, int, str]:
    opener = gzip.open if str(mtx_path).endswith(".gz") else open
    n_skip = 0
    dims_line: Optional[str] = None
    first_data: Optional[str] = None
    with opener(mtx_path, "rt", encoding="utf-8", errors="replace") as fh:
        for raw_line in fh:
            n_skip += 1
            s = raw_line.strip()
            if not s or s.startswith("%"):
                continue
            if dims_line is None:
                dims_line = s
            else:
                first_data = s
                n_skip -= 1
                break
    if dims_line is None:
        raise RuntimeError(f"Empty or malformed MTX header in {mtx_path}")
    toks = dims_line.split()
    if len(toks) < 3:
        raise RuntimeError(f"Malformed MTX dimension line {dims_line!r} in {mtx_path}")
    sep = "\t" if (first_data is not None and "\t" in first_data) else " "
    return int(toks[0]), int(toks[1]), int(toks[2]), n_skip, sep


def load_a549_mtx_controls_by_barcodes(
    raw_dir: Path,
    barcodes: Sequence[str],
    *,
    line_name: str = "a549_ko",
    pool_kind: str = "primary",
    seed: int = LOCKED_PRIMARY_SEED,
    chunk_lines: int = 1_500_000,
) -> RawControlPool:
    """
    Stream exact locked-barcode raw UMI rows from A549-KO MatrixMarket files (`GSE345058`)
    directly into memory (excluding non-human `CAS9` transgene, yielding `29,891` human RNA genes;
    zero temporary files written).
    """
    import pandas as pd

    req = _validate_requested_barcodes(barcodes, line_name)
    base = Path(raw_dir)
    adir = base / "a549_ko" if (base / "a549_ko").is_dir() else base
    cell_path = _find_matching_file(
        adir, ["GSE345058_cell_annotations.csv.gz", "GSE345058_cell_annotations.csv"], "*cell_annotations*.csv*"
    )
    gene_path = _find_matching_file(
        adir, ["GSE345058_gene_annotations.csv.gz", "GSE345058_gene_annotations.csv"], "*gene_annotations*.csv*"
    )
    mtx_path = _find_matching_file(
        adir, ["GSE345058_raw_counts.mtx.gz", "GSE345058_raw_counts.mtx"], "*raw_counts*.mtx*"
    )

    cell_df = pd.read_csv(cell_path, dtype=str, keep_default_na=False)
    gene_df = pd.read_csv(gene_path, dtype=str, keep_default_na=False)
    n_annot_cells = len(cell_df)
    n_annot_genes = len(gene_df)

    cid_col = _pick_col_ci(cell_df.columns, ["cell_id", "Cell_ID", "barcode", "Barcode", "cell"])
    if cid_col is None:
        raise KeyError(f"[{line_name}] Missing cell_id column in {cell_path}: {list(cell_df.columns)}")
    cids_raw = cell_df[cid_col].astype(str).str.strip().tolist()
    cid_to_raw_row: Dict[str, int] = {}
    for i, cid in enumerate(cids_raw):
        if cid in cid_to_raw_row:
            raise ValueError(f"[{line_name}] Duplicate cell_id {cid!r} in {cell_path}.")
        cid_to_raw_row[cid] = i

    missing = [bc for bc in req if bc not in cid_to_raw_row]
    if missing:
        raise KeyError(f"[{line_name}] {len(missing)}/{len(req)} locked barcode(s) absent from {cell_path}: {missing[:5]}")

    raw_cell_to_out = np.full(n_annot_cells, -1, dtype=np.int32)
    for out_i, bc in enumerate(req):
        raw_cell_to_out[cid_to_raw_row[bc]] = out_i

    g_sym_col = _pick_col_ci(gene_df.columns, ["gene_symbol", "gene_name", "symbol", "name", "gene"]) or gene_df.columns[0]
    g_ens_col = _pick_col_ci(gene_df.columns, ["ensembl_id", "gene_id", "gene_ids", "ensembl"])
    if g_ens_col is None:
        raise KeyError(f"[{line_name}] Missing Ensembl ID column in {gene_path}: {list(gene_df.columns)}")
    g_type_col = _pick_col_ci(gene_df.columns, ["feature_types", "feature_type", "type"])
    g_gen_col = _pick_col_ci(gene_df.columns, ["genome", "species", "reference"])

    raw_g_syms = [str(x).strip() for x in gene_df[g_sym_col].tolist()]
    raw_g_ens = [strip_ensembl_version(x) for x in gene_df[g_ens_col].tolist()]
    raw_g_gen = [str(x).strip().upper() for x in gene_df[g_gen_col].tolist()] if g_gen_col else [""] * n_annot_genes

    def _is_cas9(i_g: int) -> bool:
        return raw_g_syms[i_g].upper() == "CAS9" or raw_g_ens[i_g].upper() == "CAS9" or raw_g_gen[i_g] == "CAS9"

    if g_type_col is not None:
        g_types = [str(x).strip() for x in gene_df[g_type_col].tolist()]
        rna_idx = [i for i, t in enumerate(g_types) if t in ("Gene Expression", "RNA", "") and not _is_cas9(i)]
        if not rna_idx:
            rna_idx = [i for i in range(n_annot_genes) if not _is_cas9(i)]
    else:
        rna_idx = [i for i in range(n_annot_genes) if not _is_cas9(i)]

    raw_gene_to_out = np.full(n_annot_genes, -1, dtype=np.int32)
    out_ens: List[str] = []
    out_syms: List[str] = []
    for k, ri in enumerate(rna_idx):
        raw_gene_to_out[ri] = k
        out_ens.append(raw_g_ens[ri])
        out_syms.append(raw_g_syms[ri])

    n_mtx_rows, n_mtx_cols, _, n_skip, sep = _read_mtx_header(mtx_path)
    if n_mtx_rows == n_annot_cells and n_mtx_cols == n_annot_genes:
        cell_first = True
    elif n_mtx_rows == n_annot_genes and n_mtx_cols == n_annot_cells:
        cell_first = False
    else:
        raise RuntimeError(
            f"[{line_name}] MTX shape ({n_mtx_rows}, {n_mtx_cols}) incompatible with "
            f"(cells={n_annot_cells}, genes={n_annot_genes}) in {mtx_path}."
        )

    raw_counts = np.zeros((len(req), len(out_ens)), dtype=np.float32)
    reader = pd.read_csv(
        mtx_path,
        sep=sep,
        header=None,
        skiprows=n_skip,
        names=["r", "c", "v"],
        dtype={"r": np.int32, "c": np.int32, "v": np.float32},
        chunksize=int(chunk_lines),
        engine="c",
    )
    for chunk in reader:
        rr = chunk["r"].to_numpy(dtype=np.int64) - 1
        cc = chunk["c"].to_numpy(dtype=np.int64) - 1
        vv = chunk["v"].to_numpy(dtype=np.float32)
        c_raw = rr if cell_first else cc
        g_raw = cc if cell_first else rr
        out_c = raw_cell_to_out[c_raw]
        out_g = raw_gene_to_out[g_raw]
        keep = (out_c >= 0) & (out_g >= 0) & (vv != 0.0)
        if np.any(keep):
            np.add.at(raw_counts, (out_c[keep], out_g[keep]), vv[keep])

    return RawControlPool(
        line_name=line_name,
        pool_kind=pool_kind,
        seed=int(seed),
        barcodes=list(req),
        raw_counts=np.ascontiguousarray(raw_counts, dtype=np.float32),
        raw_ensembl_ids=out_ens,
        raw_gene_symbols=out_syms,
        source_files=[str(cell_path), str(gene_path), str(mtx_path)],
        audit={
            "loader": "load_a549_mtx_controls_by_barcodes",
            "n_cells": int(raw_counts.shape[0]),
            "n_raw_genes": int(raw_counts.shape[1]),
            "n_excluded_features": int(n_annot_genes - len(out_ens)),
            "cell_first_mtx": bool(cell_first),
        },
    )


def _parse_feng_row_label(raw_label: str) -> Tuple[str, str, bool]:
    """Parse Feng CSV row label `ENSG...:SYMBOL:Gene-Expression` -> `(ensembl_id, symbol, is_gex)`."""
    lbl = raw_label.strip().strip('"')
    if ":CRISPR-Guide-Capture" in lbl:
        return "", "", False
    parts = lbl.split(":")
    if len(parts) >= 3:
        return strip_ensembl_version(parts[0]), parts[1].strip(), (parts[-1].strip().lower() == "gene-expression")
    if len(parts) == 2:
        if parts[0].strip().startswith("ENSG"):
            return strip_ensembl_version(parts[0]), parts[1].strip(), True
        return "", parts[0].strip(), (parts[1].strip().lower() == "gene-expression")
    return "", lbl, True


def load_feng_csv_controls_by_barcodes(
    raw_dir: Path,
    barcodes: Sequence[str],
    *,
    line_name: str = "feng_ipsc",
    pool_kind: str = "primary",
    seed: int = LOCKED_PRIMARY_SEED,
) -> RawControlPool:
    """
    Stream exact locked-barcode raw UMI columns from Feng iPSC `TargetedScreen_RNA-UMI-Counts.csv.gz`
    (already D2-canonical CSV header barcodes), keeping strictly `Gene-Expression` rows (`36,518` genes).
    """
    req = _validate_requested_barcodes(barcodes, line_name)
    base = Path(raw_dir)
    fdir = base / "feng_ipsc" if (base / "feng_ipsc").is_dir() else base
    csv_path = _find_matching_file(
        fdir,
        ["TargetedScreen_RNA-UMI-Counts.csv.gz", "TargetedScreen_RNA-UMI-Counts.csv"],
        "*RNA-UMI-Counts*.csv*",
    )

    opener = gzip.open if str(csv_path).endswith(".gz") else open
    out_ens: List[str] = []
    out_syms: List[str] = []
    gene_rows: List[np.ndarray] = []

    with opener(csv_path, "rt", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh)
        header = [h.strip() for h in next(reader)]
        first_row = next(reader, None)
        if first_row is None:
            raise RuntimeError(f"[{line_name}] Empty data rows in {csv_path}")

        col_offset = len(first_row) - len(header)
        if col_offset not in (0, 1):
            raise RuntimeError(
                f"[{line_name}] Header/row column mismatch in {csv_path}: header={len(header)}, row1={len(first_row)}"
            )
        hdr_start = 1 if col_offset == 0 else 0
        hdr_to_row_col: Dict[str, int] = {}
        for h_idx in range(hdr_start, len(header)):
            cid = header[h_idx]
            if cid in hdr_to_row_col:
                raise ValueError(f"[{line_name}] Duplicate CSV header barcode {cid!r} in {csv_path}.")
            hdr_to_row_col[cid] = h_idx + col_offset

        missing = [bc for bc in req if bc not in hdr_to_row_col]
        if missing:
            raise KeyError(
                f"[{line_name}] {len(missing)}/{len(req)} locked barcode(s) absent from {csv_path}: {missing[:5]}"
            )
        sel_cols = [hdr_to_row_col[bc] for bc in req]
        n_sel = len(sel_cols)

        def _consume(row_vals: List[str]) -> None:
            if not row_vals:
                return
            ens, sym, is_gex = _parse_feng_row_label(row_vals[0])
            if not is_gex or not sym:
                return
            out_ens.append(ens)
            out_syms.append(sym)
            gene_rows.append(np.fromiter((float(row_vals[j]) for j in sel_cols), dtype=np.float32, count=n_sel))

        _consume(first_row)
        for r_vals in reader:
            _consume(r_vals)

    if not out_ens:
        raise RuntimeError(f"[{line_name}] Zero Gene-Expression rows parsed from {csv_path}.")
    raw_counts = np.ascontiguousarray(np.vstack(gene_rows).T, dtype=np.float32)

    return RawControlPool(
        line_name=line_name,
        pool_kind=pool_kind,
        seed=int(seed),
        barcodes=list(req),
        raw_counts=raw_counts,
        raw_ensembl_ids=out_ens,
        raw_gene_symbols=out_syms,
        source_files=[str(csv_path)],
        audit={
            "loader": "load_feng_csv_controls_by_barcodes",
            "n_cells": int(raw_counts.shape[0]),
            "n_raw_genes": int(raw_counts.shape[1]),
            "csv_col_offset": int(col_offset),
        },
    )


def resolve_local_nadig_h5ad_path(line_name: str, raw_nadig_dir: Optional[Path] = None) -> Path:
    """Resolve local SSD Nadig `.h5ad` path (`PREREG_I_ADDENDUM_B_DRAFT_v2.md §3.2`, never sync share)."""
    if line_name == "nadig_hepg2":
        fname = "GSE264667_hepg2_raw_singlecell_01.h5ad"
    elif line_name == "nadig_jurkat":
        fname = "GSE264667_jurkat_raw_singlecell_01.h5ad"
    else:
        raise ValueError(f"Unsupported Nadig line: {line_name!r}")

    if raw_nadig_dir is not None:
        base = Path(raw_nadig_dir)
    elif os.environ.get("NADIG_DATA_DIR"):
        base = Path(os.environ["NADIG_DATA_DIR"])
    else:
        base = DEFAULT_LOCAL_NADIG_DIR

    direct = base / fname
    if direct.is_file():
        return direct
    nested = base / "nadig2025" / fname
    if nested.is_file():
        return nested
    raise FileNotFoundError(f"Local Nadig file {fname} not found under {base}")


def load_raw_controls_by_locked_barcodes(
    line_name: str,
    barcodes: Sequence[str],
    *,
    pool_kind: str = "primary",
    seed: int = LOCKED_PRIMARY_SEED,
    repo_root: Path = REPO_ROOT,
    raw_nadig_dir: Optional[Path] = None,
    raw_dir_override: Optional[Path] = None,
    verify_nadig_ntc: bool = True,
) -> RawControlPool:
    """Unified read-only dispatcher loading exact raw control rows for any local Addendum B line."""
    if line_name in ("nadig_hepg2", "nadig_jurkat"):
        h5ad = resolve_local_nadig_h5ad_path(line_name, raw_nadig_dir=raw_nadig_dir)
        return load_nadig_controls_by_barcodes(
            h5ad,
            barcodes,
            line_name=line_name,
            pool_kind=pool_kind,
            seed=seed,
            verify_ntc=verify_nadig_ntc,
        )
    if line_name == "dld1":
        raw_d = raw_dir_override or (Path(repo_root) / "round10" / "shared" / "raw")
        pmap, ref_pref = build_dld1_prefix_map(raw_d)
        return load_10x_multi_h5_controls_by_barcodes(
            pmap, barcodes, line_name=line_name, canonical_ref_prefix=ref_pref, pool_kind=pool_kind, seed=seed
        )
    if line_name == "cd4t_rest":
        raw_d = raw_dir_override or (Path(repo_root) / "round10" / "shared" / "raw")
        pmap, ref_pref = build_cd4t_rest_prefix_map(raw_d)
        return load_10x_multi_h5_controls_by_barcodes(
            pmap, barcodes, line_name=line_name, canonical_ref_prefix=ref_pref, pool_kind=pool_kind, seed=seed
        )
    if line_name == "a549_ko":
        raw_d = raw_dir_override or (Path(repo_root) / "round12" / "shared" / "raw")
        return load_a549_mtx_controls_by_barcodes(
            raw_d, barcodes, line_name=line_name, pool_kind=pool_kind, seed=seed
        )
    if line_name == "feng_ipsc":
        raw_d = raw_dir_override or (Path(repo_root) / "round11" / "shared" / "raw")
        return load_feng_csv_controls_by_barcodes(
            raw_d, barcodes, line_name=line_name, pool_kind=pool_kind, seed=seed
        )
    raise ValueError(f"Unsupported Addendum B raw control line: {line_name!r}")
