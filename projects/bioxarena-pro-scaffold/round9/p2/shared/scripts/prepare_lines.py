#!/usr/bin/env python3
"""
Unified Addendum-B Single-Cell-to-Shared-Line Preprocessing Pipeline (A-T1).

Produces standardized per-cell-line directories under `p2/shared/lines/<line>/`:
  - `genes.txt`            : official HGNC gene symbols (aligned, deduplicated, finite across all batches)
  - `targets.txt`          : perturbed target gene symbols (deduplicated 1 row per target)
  - `ctrl_cells.npz`       : randomly subsampled (<= 20,000, batch-stratified) non-targeting control cells
                             with `X_log1p`, `X_zscore`, `X_libnorm`, `ctrl_mean_log1p`, `ctrl_mean_libnorm`,
                             `n_ctrl_total`, `batch_labels`
  - `pseudobulk_delta.npz` : per-target pseudobulk deltas (`delta_log1p`, `delta_zscore`, `delta_libnorm`),
                             disjoint 50/50 split-half deltas (`half_{a,b}_delta_{log1p,zscore,libnorm}`,
                             `half_valid`), `targets`, `target_ensembl_ids`, `reproducible`, `n_cells`,
                             `ad_counts`, `energy_test_p_value`
  - `line_metadata.json`   : provenance, SHA-256, and alignment/filtering loss audit

Reuses verified Addendum-B statistical tests (`ad_two_sample_pvalues_vectorized`,
`energy_distance_perm_pvalue`, `stratified_nt_sample`) from `src/prepare_nadig2025_cache.py`.
"""

from __future__ import annotations

import argparse
import gc
import json
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.decomposition import TruncatedSVD
from statsmodels.stats.multitest import multipletests

PROJECT_ROOT = Path(__file__).resolve().parents[4]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from prepare_nadig2025_cache import (
    ad_two_sample_pvalues_vectorized,
    energy_distance_perm_pvalue,
    sha256_file,
    stratified_nt_sample,
    verify_ad_implementation,
)

SHARED_DIR = Path(__file__).resolve().parents[1]
LINES_DIR = SHARED_DIR / "lines"


def decode_str_array(arr) -> List[str]:
    return [x.decode("utf-8", errors="replace") if isinstance(x, (bytes, np.bytes_)) else str(x) for x in arr]


def read_obs_or_var_col(grp: h5py.Group, col: str) -> List[str]:
    """Read an AnnData obs/var column handling both categorical (legacy & new h5ad) and flat datasets."""
    node = grp[col]
    if isinstance(node, h5py.Group):
        cats = decode_str_array(node["categories"][:])
        codes = node["codes"][:]
        return [cats[int(c)] if int(c) >= 0 else "" for c in codes]
    if "__categories" in grp and col in grp["__categories"] and np.issubdtype(node.dtype, np.integer):
        cats = decode_str_array(grp["__categories"][col][:])
        codes = node[:]
        return [cats[int(c)] if int(c) >= 0 else "" for c in codes]
    return decode_str_array(node[:])


def read_h5ad_chunk_dense(x_node, s: int, e: int, n_genes_raw: int) -> np.ndarray:
    """Read rows [s:e] from a dense HDF5 dataset, CSR HDF5 group, or in-memory scipy CSR matrix."""
    if sp.issparse(x_node):
        return x_node[s:e].toarray().astype(np.float64, copy=False)
    if isinstance(x_node, h5py.Dataset):
        return np.asarray(x_node[s:e, :], dtype=np.float64)
    if isinstance(x_node, h5py.Group):
        enc = str(x_node.attrs.get("encoding-type", ""))
        if "csc" in enc.lower():
            raise NotImplementedError("CSC HDF5 node must be pre-converted to CSR via open_x_node().")
        indptr = x_node["indptr"][s : e + 1]
        p0, p1 = int(indptr[0]), int(indptr[-1])
        data = x_node["data"][p0:p1]
        indices = x_node["indices"][p0:p1]
        rel_ptr = indptr - p0
        return sp.csr_matrix((data, indices, rel_ptr), shape=(e - s, n_genes_raw)).toarray().astype(np.float64, copy=False)
    raise TypeError(f"Unsupported HDF5 X node type: {type(x_node)}")


def open_x_node(f: h5py.File, x_layer: str = "X"):
    """Open HDF5 X/layer node; convert CSC sparse groups to fast in-memory CSR."""
    x_node = f[x_layer] if x_layer == "X" else f[f"layers/{x_layer}"]
    if isinstance(x_node, h5py.Group):
        enc = str(x_node.attrs.get("encoding-type", ""))
        shape = tuple(int(v) for v in x_node.attrs.get("shape"))
        if "csc" in enc.lower():
            csc = sp.csc_matrix((x_node["data"][:], x_node["indices"][:], x_node["indptr"][:]), shape=shape)
            return csc.tocsr(), shape
        return x_node, shape
    if isinstance(x_node, h5py.Dataset):
        return x_node, (int(x_node.shape[0]), int(x_node.shape[1]))
    raise TypeError(f"Unsupported X node: {type(x_node)}")


def compute_row_umi_sums(x_node, n_cells: int, n_genes_raw: int, chunk_size: int = 8192) -> np.ndarray:
    """Compute per-cell total UMI counts directly from X when obs/UMI_count is absent."""
    if sp.issparse(x_node):
        return np.asarray(x_node.sum(axis=1), dtype=np.float64).ravel()
    if isinstance(x_node, h5py.Group):
        indptr = x_node["indptr"][:]
        data = x_node["data"][:]
        umi = np.empty(n_cells, dtype=np.float64)
        for s in range(0, n_cells, chunk_size):
            e = min(s + chunk_size, n_cells)
            ptrs = indptr[s : e + 1]
            p0, p1 = int(ptrs[0]), int(ptrs[-1])
            seg = data[p0:p1].astype(np.float64, copy=False)
            lens = np.diff(ptrs)
            if len(seg) == 0:
                umi[s:e] = 0.0
            elif np.all(lens > 0):
                starts = np.empty(len(lens), dtype=np.int64)
                starts[0] = 0
                np.cumsum(lens[:-1], out=starts[1:])
                umi[s:e] = np.add.reduceat(seg, starts)
            else:
                rel = ptrs - p0
                umi[s:e] = [float(np.sum(seg[rel[i] : rel[i + 1]])) if lens[i] > 0 else 0.0 for i in range(len(lens))]
        return umi
    umi = np.zeros(n_cells, dtype=np.float64)
    for s in range(0, n_cells, chunk_size):
        e = min(s + chunk_size, n_cells)
        umi[s:e] = np.asarray(x_node[s:e, :], dtype=np.float64).sum(axis=1)
    return umi


def collect_existing_line_symbols(exclude_line: str) -> set:
    """Return union of gene & target symbols already present in other prepared lines."""
    syms = set()
    if not LINES_DIR.exists():
        return syms
    for ld in LINES_DIR.iterdir():
        if not ld.is_dir() or ld.name == exclude_line:
            continue
        for fn in ("genes.txt", "targets.txt"):
            p = ld / fn
            if p.exists():
                syms.update(x.strip() for x in p.read_text(encoding="utf-8").splitlines() if x.strip())
    return syms


def build_gene_keep_mask(
    var_gnames: Sequence[str],
    finite_mask: np.ndarray,
    control_mean: np.ndarray,
    extra_keep_symbols: Optional[set] = None,
    max_genes: Optional[int] = None,
) -> Tuple[np.ndarray, Dict[str, object]]:
    """
    Filter to finite non-empty gene symbols, deduplicate symbols deterministically
    (retaining highest control_mean), and optionally cap to top `max_genes` while
    force-retaining `extra_keep_symbols`.
    """
    n_raw = len(var_gnames)
    n_finite = int(finite_mask.sum())
    n_empty = 0
    best_by_sym: Dict[str, Tuple[float, int]] = {}
    dup_counter: Dict[str, int] = defaultdict(int)

    for idx in np.where(finite_mask)[0]:
        sym = str(var_gnames[idx]).strip()
        if not sym or sym.lower() in ("nan", "none", "na"):
            n_empty += 1
            continue
        dup_counter[sym] += 1
        score = float(control_mean[idx])
        if sym not in best_by_sym or score > best_by_sym[sym][0]:
            best_by_sym[sym] = (score, int(idx))

    dup_syms = sorted([s for s, c in dup_counter.items() if c > 1])
    n_dup_dropped = sum(c - 1 for c in dup_counter.values() if c > 1)

    kept_items = list(best_by_sym.items())  # [(sym, (score, idx)), ...]
    n_pre_cap = len(kept_items)
    if max_genes is not None and len(kept_items) > max_genes:
        keep_set = extra_keep_symbols or set()
        forced = [it for it in kept_items if it[0] in keep_set]
        others = sorted([it for it in kept_items if it[0] not in keep_set], key=lambda x: x[1][0], reverse=True)
        n_rem = max(0, max_genes - len(forced))
        kept_items = forced + others[:n_rem]

    selected = sorted([it[1][1] for it in kept_items])
    audit = {
        "n_genes_raw": int(n_raw),
        "n_genes_finite_all_batches": int(n_finite),
        "n_genes_empty_symbol_dropped": int(n_empty),
        "n_genes_duplicate_symbol_dropped": int(n_dup_dropped),
        "n_genes_low_expr_cap_dropped": int(n_pre_cap - len(selected)),
        "duplicate_symbols": dup_syms,
        "n_genes_final": int(len(selected)),
    }
    return np.asarray(selected, dtype=np.int64), audit


def stratified_split_halves(
    indices: np.ndarray,
    batches: np.ndarray,
    seed: int = 42,
) -> Tuple[np.ndarray, np.ndarray]:
    """Deterministic 50/50 stratified split within batches."""
    rng = np.random.RandomState(seed)
    if len(indices) < 2:
        return indices.copy(), np.array([], dtype=indices.dtype)
    b_vals = batches[indices]
    uniq_b = np.unique(b_vals)
    ha, hb = [], []
    toggle = 0
    for b in uniq_b:
        pool = indices[b_vals == b].copy()
        rng.shuffle(pool)
        n = len(pool)
        mid = n // 2
        if n % 2 == 1:
            if toggle % 2 == 0:
                mid += 1
            toggle += 1
        ha.append(pool[:mid])
        hb.append(pool[mid:])
    a_arr = np.sort(np.concatenate(ha)) if ha else np.array([], dtype=indices.dtype)
    b_arr = np.sort(np.concatenate(hb)) if hb else np.array([], dtype=indices.dtype)
    return a_arr, b_arr


def process_singlecell_line(
    line_name: str,
    h5ad_paths: "Path | Sequence[Tuple[Path, str]]",
    adapter_fn: Callable[..., Dict[str, object]],
    max_ctrl_cells: int = 20000,
    max_genes: Optional[int] = None,
    extra_keep_symbols: Optional[set] = None,
    precomputed_stats_fn: Optional[Callable[[], Dict[str, Tuple[int, float, float]]]] = None,
    chunk_size: int = 8192,
    seed: int = 42,
    x_layer: str = "X",
) -> Path:
    """
    Generic chunked Addendum-B single-cell processor for one or more h5ad files.

    `h5ad_paths` can be a single `Path` or a sequence of `(Path, file_tag)` pairs
    concatenated virtually along cells (intersecting gene symbols when panels differ).
    `adapter_fn(f, tag="")` returns dict with:
      - `nt_mask`: bool array (n_cells_part,)
      - `batches`: int/str array (n_cells_part,)
      - `umi`: float64 array (n_cells_part,) or None (computed from X/layer)
      - `pert_keys`: list[str] (n_cells_part,) unique perturbation group key ("" for NT/invalid)
      - `pert_meta`: dict[key -> (ens_id, gene_symbol)]
      - `var_gene_ids`: list[str] (n_genes_part,)
      - `var_gene_names`: list[str] (n_genes_part,)
    """
    from contextlib import ExitStack

    t0 = time.time()
    out_dir = LINES_DIR / line_name
    out_dir.mkdir(parents=True, exist_ok=True)
    specs: List[Tuple[Path, str]] = (
        [(h5ad_paths, "")] if isinstance(h5ad_paths, Path) else [(Path(p), str(tag)) for p, tag in h5ad_paths]
    )
    print(f"\n=== [{line_name}] Reading {[p.name for p, _ in specs]} ===", flush=True)
    input_sha256s = {p.name: sha256_file(p) for p, _ in specs}
    input_sha256 = input_sha256s[specs[0][0].name] if len(specs) == 1 else json.dumps(input_sha256s)

    with ExitStack() as stack:
        files = [stack.enter_context(h5py.File(p, "r")) for p, _ in specs]
        x_nodes, part_shapes, part_metas, part_umis = [], [], [], []
        for f_i, (_, tag) in zip(files, specs):
            xn, shp = open_x_node(f_i, x_layer=x_layer)
            x_nodes.append(xn)
            part_shapes.append(shp)
            try:
                m_i = adapter_fn(f_i, tag)
            except TypeError:
                m_i = adapter_fn(f_i)
            part_metas.append(m_i)
            if m_i["umi"] is None:
                print(f"  Computing per-cell UMI totals from {x_layer} ({tag or specs[0][0].name})...", flush=True)
                part_umis.append(compute_row_umi_sums(xn, shp[0], shp[1], chunk_size=chunk_size))
            else:
                part_umis.append(np.asarray(m_i["umi"], dtype=np.float64))

        # Align var columns across files if multi-file
        if len(specs) == 1:
            col_maps: List[Optional[np.ndarray]] = [None]
            var_gids = list(part_metas[0]["var_gene_ids"])
            var_gnames = list(part_metas[0]["var_gene_names"])
            n_multifile_dropped = 0
        else:
            ref_names = [str(x).strip() for x in part_metas[0]["var_gene_names"]]
            if all([str(x).strip() for x in m["var_gene_names"]] == ref_names for m in part_metas[1:]):
                col_maps = [None] * len(specs)
                var_gids = list(part_metas[0]["var_gene_ids"])
                var_gnames = ref_names
                n_multifile_dropped = 0
            else:
                # Deterministic intersection of unique non-empty symbols in file-0 order
                sym_maps = []
                for m in part_metas:
                    smap: Dict[str, int] = {}
                    for ci, nm in enumerate(m["var_gene_names"]):
                        s = str(nm).strip()
                        if s and s not in smap:
                            smap[s] = ci
                    sym_maps.append(smap)
                common_syms = [
                    s for s in sym_maps[0].keys() if all(s in sm for sm in sym_maps[1:])
                ]
                col_maps = [np.array([sm[s] for s in common_syms], dtype=np.int64) for sm in sym_maps]
                var_gids = [str(part_metas[0]["var_gene_ids"][sym_maps[0][s]]) for s in common_syms]
                var_gnames = common_syms
                n_multifile_dropped = int(max(shp[1] for shp in part_shapes) - len(common_syms))
        n_dup_nt_skipped = int(sum(m.get("n_duplicate_nt_skipped", 0) for m in part_metas))

        n_genes_raw = len(var_gnames)
        offsets = np.cumsum([0] + [shp[0] for shp in part_shapes[:-1]]).tolist()
        n_cells = int(sum(shp[0] for shp in part_shapes))

        def read_virtual_chunk(s: int, e: int) -> np.ndarray:
            parts = []
            for idx_f, (xn, shp, cmap, off) in enumerate(zip(x_nodes, part_shapes, col_maps, offsets)):
                fs, fe = off, off + shp[0]
                lo, hi = max(s, fs), min(e, fe)
                if lo >= hi:
                    continue
                arr = read_h5ad_chunk_dense(xn, lo - fs, hi - fs, shp[1])
                if cmap is not None:
                    arr = arr[:, cmap]
                parts.append(arr)
            return parts[0] if len(parts) == 1 else np.vstack(parts)

        nt_mask = np.concatenate([np.asarray(m["nt_mask"], dtype=bool) for m in part_metas])
        batches_raw = np.concatenate([np.asarray(m["batches"], dtype=str) for m in part_metas])
        umi = np.concatenate(part_umis)
        pert_keys: List[str] = [pk for m in part_metas for pk in m["pert_keys"]]
        pert_meta: Dict[str, Tuple[str, str]] = {}
        for m in part_metas:
            pert_meta.update(m["pert_meta"])
        target_remaps: Dict[str, str] = {}
        for m in part_metas:
            target_remaps.update(m.get("target_remaps", {}))

        uniq_batches, batch_idx_per_cell = np.unique(batches_raw, return_inverse=True)
        n_batches = len(uniq_batches)

        # Split NT cells 50/50 per batch for disjoint split-half noise ceiling
        nt_all_idx = np.where(nt_mask)[0]
        nt_half_a_idx, nt_half_b_idx = stratified_split_halves(nt_all_idx, batch_idx_per_cell, seed=seed)
        nt_half_flag = np.zeros(n_cells, dtype=np.int8)  # 1 = half A, 2 = half B
        nt_half_flag[nt_half_a_idx] = 1
        nt_half_flag[nt_half_b_idx] = 2

        # Pass 1: Sequential scan over chunks for exact per-batch NT statistics (full + half A + half B)
        print(f"  [Pass 1] NT stats across {n_cells} cells x {n_genes_raw} genes ({n_batches} batches)...", flush=True)
        nt_sum_lib = np.zeros((n_batches, n_genes_raw), dtype=np.float64)
        nt_sq_lib = np.zeros((n_batches, n_genes_raw), dtype=np.float64)
        nt_sum_log = np.zeros((n_batches, n_genes_raw), dtype=np.float64)
        nt_cnt = np.zeros(n_batches, dtype=np.int64)

        ha_sum_lib = np.zeros((n_batches, n_genes_raw), dtype=np.float64)
        ha_sq_lib = np.zeros((n_batches, n_genes_raw), dtype=np.float64)
        ha_sum_log = np.zeros((n_batches, n_genes_raw), dtype=np.float64)
        ha_cnt = np.zeros(n_batches, dtype=np.int64)

        hb_sum_lib = np.zeros((n_batches, n_genes_raw), dtype=np.float64)
        hb_sq_lib = np.zeros((n_batches, n_genes_raw), dtype=np.float64)
        hb_sum_log = np.zeros((n_batches, n_genes_raw), dtype=np.float64)
        hb_cnt = np.zeros(n_batches, dtype=np.int64)

        for s in range(0, n_cells, chunk_size):
            e = min(s + chunk_size, n_cells)
            m_loc = nt_mask[s:e]
            if not np.any(m_loc):
                continue
            Xc = read_virtual_chunk(s, e)
            rows = np.where(m_loc)[0]
            X_lib = Xc[rows] * (1e4 / np.maximum(umi[s:e][rows], 1.0))[:, None]
            X_log = np.log1p(X_lib)
            b_nt = batch_idx_per_cell[s:e][rows]
            hf = nt_half_flag[s:e][rows]

            np.add.at(nt_sum_lib, b_nt, X_lib)
            np.add.at(nt_sq_lib, b_nt, X_lib * X_lib)
            np.add.at(nt_sum_log, b_nt, X_log)
            np.add.at(nt_cnt, b_nt, 1)

            ma = hf == 1
            if np.any(ma):
                np.add.at(ha_sum_lib, b_nt[ma], X_lib[ma])
                np.add.at(ha_sq_lib, b_nt[ma], X_lib[ma] * X_lib[ma])
                np.add.at(ha_sum_log, b_nt[ma], X_log[ma])
                np.add.at(ha_cnt, b_nt[ma], 1)

            mb = hf == 2
            if np.any(mb):
                np.add.at(hb_sum_lib, b_nt[mb], X_lib[mb])
                np.add.at(hb_sq_lib, b_nt[mb], X_lib[mb] * X_lib[mb])
                np.add.at(hb_sum_log, b_nt[mb], X_log[mb])
                np.add.at(hb_cnt, b_nt[mb], 1)

        mu_lib = nt_sum_lib / np.maximum(nt_cnt[:, None], 1)
        var_lib = np.maximum((nt_sq_lib - nt_cnt[:, None] * (mu_lib ** 2)) / np.maximum(nt_cnt[:, None] - 1, 1), 0.0)
        sd_lib = np.sqrt(var_lib)
        mu_log = nt_sum_log / np.maximum(nt_cnt[:, None], 1)

        nt_tot = max(int(nt_cnt.sum()), 1)
        ctrl_mean_lib = nt_sum_lib.sum(axis=0) / nt_tot
        ctrl_mean_log = nt_sum_log.sum(axis=0) / nt_tot
        var_global_lib = np.maximum((nt_sq_lib.sum(axis=0) - nt_tot * (ctrl_mean_lib ** 2)) / max(nt_tot - 1, 1), 0.0)
        sd_global_lib = np.sqrt(var_global_lib)

        # Require non-zero std across all batches with >= 15 NT cells (and overall NT std > 1e-12)
        # so 1-2 tiny batches (e.g. HepG2 gem_groups with 5-7 NT cells) do not drop ~2,700 valid genes.
        well_pop = nt_cnt >= 15
        if np.any(well_pop):
            finite_batch_ok = np.all(np.isfinite(mu_lib[well_pop]) & np.isfinite(sd_lib[well_pop]) & (sd_lib[well_pop] > 1e-12), axis=0)
        else:
            finite_batch_ok = np.ones(n_genes_raw, dtype=bool)
        finite_mask = (
            finite_batch_ok
            & np.all(np.isfinite(mu_lib), axis=0)
            & np.isfinite(ctrl_mean_lib)
            & (ctrl_mean_lib > 0.0)
            & np.isfinite(sd_global_lib)
            & (sd_global_lib > 1e-12)
        )
        # Fill sd_lib in small batches (<15 NT cells) with global NT sd if zero
        sd_lib = np.where((nt_cnt[:, None] >= 15) & (sd_lib > 1e-12), sd_lib, sd_global_lib[None, :])

        # Force-include target gene symbols that exist in pert_meta if finite
        target_sym_set = {sym for _, sym in pert_meta.values()}
        keep_union = set(target_sym_set)
        if extra_keep_symbols:
            keep_union |= set(extra_keep_symbols)

        sel_idx, gene_audit = build_gene_keep_mask(
            var_gnames,
            finite_mask=finite_mask,
            control_mean=ctrl_mean_lib,
            extra_keep_symbols=keep_union,
            max_genes=max_genes,
        )
        # Also prepare top-5000 finite index for exact Addendum B AD/Energy test if needed
        ok_raw_idx = np.where(finite_mask)[0]
        top5000_raw_idx = np.sort(ok_raw_idx[np.argsort(ctrl_mean_lib[ok_raw_idx])[-min(5000, len(ok_raw_idx)) :]])

        sel_gnames = [str(var_gnames[i]).strip() for i in sel_idx]
        sel_gids = [str(var_gids[i]).strip() for i in sel_idx]
        G = len(sel_idx)
        print(
            f"  [Pass 1 DONE] finite={gene_audit['n_genes_finite_all_batches']}/{n_genes_raw}, "
            f"final G={G} (dup_dropped={gene_audit['n_genes_duplicate_symbol_dropped']}), "
            f"nt_total={len(nt_all_idx)} ({time.time()-t0:.1f}s)",
            flush=True,
        )

        # Half-A and Half-B batch NT stats on selected genes
        mu_lib_s = mu_lib[:, sel_idx].astype(np.float32)
        sd_lib_s = sd_lib[:, sel_idx].astype(np.float32)
        mu_log_s = mu_log[:, sel_idx].astype(np.float32)

        ha_mu_lib = (ha_sum_lib[:, sel_idx] / np.maximum(ha_cnt[:, None], 1)).astype(np.float32)
        ha_var_lib = np.maximum(
            (ha_sq_lib[:, sel_idx] - ha_cnt[:, None] * (ha_mu_lib.astype(np.float64) ** 2))
            / np.maximum(ha_cnt[:, None] - 1, 1),
            0.0,
        )
        ha_sd_raw = np.sqrt(ha_var_lib).astype(np.float32)
        ha_sd_lib = np.where((ha_cnt[:, None] >= 10) & (ha_sd_raw > 1e-12), ha_sd_raw, sd_lib_s)
        ha_mu_log = (ha_sum_log[:, sel_idx] / np.maximum(ha_cnt[:, None], 1)).astype(np.float32)

        hb_mu_lib = (hb_sum_lib[:, sel_idx] / np.maximum(hb_cnt[:, None], 1)).astype(np.float32)
        hb_var_lib = np.maximum(
            (hb_sq_lib[:, sel_idx] - hb_cnt[:, None] * (hb_mu_lib.astype(np.float64) ** 2))
            / np.maximum(hb_cnt[:, None] - 1, 1),
            0.0,
        )
        hb_sd_raw = np.sqrt(hb_var_lib).astype(np.float32)
        hb_sd_lib = np.where((hb_cnt[:, None] >= 10) & (hb_sd_raw > 1e-12), hb_sd_raw, sd_lib_s)
        hb_mu_log = (hb_sum_log[:, sel_idx] / np.maximum(hb_cnt[:, None], 1)).astype(np.float32)

        # Map non-NT perturbation cell groups
        pkey_to_rows: Dict[str, List[int]] = defaultdict(list)
        for i in range(n_cells):
            if nt_mask[i]:
                continue
            pk = pert_keys[i]
            if pk:
                pkey_to_rows[pk].append(i)

        all_pkeys = sorted(pkey_to_rows.keys())
        pkey_to_idx = {pk: i for i, pk in enumerate(all_pkeys)}
        cell_pkey_idx = np.full(n_cells, -1, dtype=np.int32)
        pert_half_flag = np.zeros(n_cells, dtype=np.int8)  # 1 = half A, 2 = half B
        rng_split = np.random.RandomState(seed)
        for pk, rlist in pkey_to_rows.items():
            pi = pkey_to_idx[pk]
            arr = np.asarray(rlist, dtype=np.int64)
            cell_pkey_idx[arr] = pi
            if len(arr) >= 2:
                ha, hb = stratified_split_halves(arr, batch_idx_per_cell, seed=seed + pi)
                pert_half_flag[ha] = 1
                pert_half_flag[hb] = 2

        n_pk = len(all_pkeys)

        # Subsample <= max_ctrl_cells NT control cells (stratified by batch)
        ctrl_sample_idx = stratified_nt_sample(nt_all_idx, batch_idx_per_cell, n_target=max_ctrl_cells, seed=seed)
        ctrl_pos_map = np.full(n_cells, -1, dtype=np.int32)
        ctrl_pos_map[ctrl_sample_idx] = np.arange(len(ctrl_sample_idx), dtype=np.int32)
        N_ctrl = len(ctrl_sample_idx)

        X_ctrl_lib = np.empty((N_ctrl, G), dtype=np.float32)
        X_ctrl_log = np.empty((N_ctrl, G), dtype=np.float32)
        X_ctrl_z = np.empty((N_ctrl, G), dtype=np.float32)

        # Accumulate full and split-half pseudobulk deltas
        pb_cnt = np.zeros(n_pk, dtype=np.int64)
        pb_sum_lib = np.zeros((n_pk, G), dtype=np.float64)
        pb_sum_log = np.zeros((n_pk, G), dtype=np.float64)
        pb_sum_z = np.zeros((n_pk, G), dtype=np.float64)

        ha_pb_cnt = np.zeros(n_pk, dtype=np.int64)
        ha_pb_sum_lib = np.zeros((n_pk, G), dtype=np.float64)
        ha_pb_sum_log = np.zeros((n_pk, G), dtype=np.float64)
        ha_pb_sum_z = np.zeros((n_pk, G), dtype=np.float64)

        hb_pb_cnt = np.zeros(n_pk, dtype=np.int64)
        hb_sum_pb_lib = np.zeros((n_pk, G), dtype=np.float64)
        hb_sum_pb_log = np.zeros((n_pk, G), dtype=np.float64)
        hb_sum_pb_z = np.zeros((n_pk, G), dtype=np.float64)

        # Check if we need to run AD + Energy tests (no precomputed stats)
        keep_z = nt_mask | (cell_pkey_idx >= 0)
        need_ad_energy = precomputed_stats_fn is None
        if need_ad_energy:
            keep_z_idx = np.where(keep_z)[0]
            z_pos_map = np.full(n_cells, -1, dtype=np.int32)
            z_pos_map[keep_z_idx] = np.arange(len(keep_z_idx), dtype=np.int32)
            G_top = len(top5000_raw_idx)
            Z_5000 = np.empty((len(keep_z_idx), G_top), dtype=np.float32)
            mu_top = mu_lib[:, top5000_raw_idx].astype(np.float32)
            sd_top = sd_lib[:, top5000_raw_idx].astype(np.float32)

        t2 = time.time()
        print(f"  [Pass 2] Building ctrl_cells ({N_ctrl} x {G}) + pseudobulk deltas ({n_pk} groups)...", flush=True)
        for s in range(0, n_cells, chunk_size):
            e = min(s + chunk_size, n_cells)
            if not np.any(keep_z[s:e]):
                continue
            Xc_raw = read_virtual_chunk(s, e)
            scale = (1e4 / np.maximum(umi[s:e], 1.0)).astype(np.float32)
            bi = batch_idx_per_cell[s:e]

            if need_ad_energy:
                zp = z_pos_map[s:e]
                z_loc = np.where(zp >= 0)[0]
                if len(z_loc) > 0:
                    Zt = Xc_raw[z_loc][:, top5000_raw_idx].astype(np.float32)
                    Zt *= scale[z_loc][:, None]
                    Zt -= mu_top[bi[z_loc]]
                    Zt /= sd_top[bi[z_loc]]
                    Z_5000[zp[z_loc]] = Zt

            Xs_lib = Xc_raw[:, sel_idx].astype(np.float32) * scale[:, None]
            Xs_log = np.log1p(Xs_lib)
            Xs_z = (Xs_lib - mu_lib_s[bi]) / sd_lib_s[bi]

            # Store sampled control cells
            cpos = ctrl_pos_map[s:e]
            c_loc = np.where(cpos >= 0)[0]
            if len(c_loc) > 0:
                dst = cpos[c_loc]
                X_ctrl_lib[dst] = Xs_lib[c_loc]
                X_ctrl_log[dst] = Xs_log[c_loc]
                X_ctrl_z[dst] = Xs_z[c_loc]

            # Accumulate perturbations
            pidx = cell_pkey_idx[s:e]
            p_loc = np.where(pidx >= 0)[0]
            if len(p_loc) > 0:
                pi_loc = pidx[p_loc]
                bi_loc = bi[p_loc]
                d_lib = Xs_lib[p_loc] - mu_lib_s[bi_loc]
                d_log = Xs_log[p_loc] - mu_log_s[bi_loc]
                d_z = Xs_z[p_loc]

                np.add.at(pb_sum_lib, pi_loc, d_lib)
                np.add.at(pb_sum_log, pi_loc, d_log)
                np.add.at(pb_sum_z, pi_loc, d_z)
                np.add.at(pb_cnt, pi_loc, 1)

                hf = pert_half_flag[s:e][p_loc]
                ma = hf == 1
                if np.any(ma):
                    pi_a, bi_a, loc_a = pi_loc[ma], bi_loc[ma], p_loc[ma]
                    np.add.at(ha_pb_sum_lib, pi_a, Xs_lib[loc_a] - ha_mu_lib[bi_a])
                    np.add.at(ha_pb_sum_log, pi_a, Xs_log[loc_a] - ha_mu_log[bi_a])
                    np.add.at(ha_pb_sum_z, pi_a, (Xs_lib[loc_a] - ha_mu_lib[bi_a]) / ha_sd_lib[bi_a])
                    np.add.at(ha_pb_cnt, pi_a, 1)

                mb = hf == 2
                if np.any(mb):
                    pi_b, bi_b, loc_b = pi_loc[mb], bi_loc[mb], p_loc[mb]
                    np.add.at(hb_sum_pb_lib, pi_b, Xs_lib[loc_b] - hb_mu_lib[bi_b])
                    np.add.at(hb_sum_pb_log, pi_b, Xs_log[loc_b] - hb_mu_log[bi_b])
                    np.add.at(hb_sum_pb_z, pi_b, (Xs_lib[loc_b] - hb_mu_lib[bi_b]) / hb_sd_lib[bi_b])
                    np.add.at(hb_pb_cnt, pi_b, 1)

            if (e // chunk_size) % 10 == 0 or e == n_cells:
                print(f"    chunk {e}/{n_cells} ({time.time()-t2:.1f}s)", flush=True)

    # Reproducibility stats: either precomputed or computed via exact Addendum B
    if precomputed_stats_fn is not None:
        stats_map = precomputed_stats_fn()
    else:
        print(f"  [Stats] Computing exact Addendum-B AD + 20-PC Energy tests on {n_pk} groups...", flush=True)
        verify_ad_implementation()
        svd = TruncatedSVD(n_components=20, random_state=42)
        Z_pca20 = svd.fit_transform(Z_5000).astype(np.float32)
        nt_ad_idx = stratified_nt_sample(nt_all_idx, batch_idx_per_cell, n_target=1000, seed=42)
        nt_en_idx = stratified_nt_sample(nt_all_idx, batch_idx_per_cell, n_target=500, seed=42)
        X_ctrl_ad = Z_5000[z_pos_map[nt_ad_idx]]
        C_pca_en = Z_pca20[z_pos_map[nt_en_idx]]
        sigma_cache: Dict[Tuple[int, int], float] = {}
        rng_en = np.random.RandomState(42)
        stats_map = {}
        for k_i, pk in enumerate(all_pkeys):
            ridx = np.asarray(pkey_to_rows[pk], dtype=np.int64)
            n_c = len(ridx)
            if n_c < 25:
                stats_map[pk] = (0, 0.0, 1.0)
                continue
            rz = z_pos_map[ridx]
            p_ad = ad_two_sample_pvalues_vectorized(Z_5000[rz], X_ctrl_ad, sigma_cache)
            _, q_ad, _, _ = multipletests(p_ad, method="fdr_bh")
            ad_cnt = float(np.sum(q_ad < 0.05))
            ep_val = float(energy_distance_perm_pvalue(Z_pca20[rz], C_pca_en, rng_en, n_perm=2000, batch_size=100))
            rep = 1 if (ad_cnt >= 50.0 and ep_val < 0.001) else 0
            stats_map[pk] = (rep, ad_cnt, ep_val)
            if (k_i + 1) % 100 == 0 or (k_i + 1) == len(all_pkeys):
                print(f"    [{k_i+1}/{len(all_pkeys)}] rep_so_far={sum(v[0] for v in stats_map.values())}", flush=True)
        del Z_5000, Z_pca20, X_ctrl_ad, C_pca_en
        gc.collect()

    # Deduplicate to 1 row per unique target gene symbol (retain n_cells >= 10)
    best_by_sym: Dict[str, Tuple[Tuple[int, float, float], int, str, str]] = {}
    n_pert_min10 = 0
    for pk in all_pkeys:
        pi = pkey_to_idx[pk]
        n_c = int(pb_cnt[pi])
        if n_c < 10:
            continue
        n_pert_min10 += 1
        ens, sym = pert_meta[pk]
        sym = str(sym).strip()
        if not sym or sym.lower() in ("nan", "none", "na", "non-targeting"):
            continue
        rep, ad_val, ep_val = stats_map.get(pk, (0, 0.0, 1.0))
        rep = 1 if (n_c >= 25 and ad_val >= 50.0 and ep_val < 0.001) else 0
        key_score = (int(rep), float(ad_val), float(n_c))
        if sym not in best_by_sym or key_score > best_by_sym[sym][0]:
            best_by_sym[sym] = (key_score, pi, ens, pk)

    sorted_syms = sorted(best_by_sym.keys())
    P = len(sorted_syms)
    sel_pi = np.array([best_by_sym[s][1] for s in sorted_syms], dtype=np.int64)
    target_ens = np.array([best_by_sym[s][2] for s in sorted_syms], dtype=str)
    reproducible = np.array([best_by_sym[s][0][0] == 1 for s in sorted_syms], dtype=bool)
    ad_counts_arr = np.array([best_by_sym[s][0][1] for s in sorted_syms], dtype=np.float32)
    ep_arr = np.array([stats_map.get(best_by_sym[s][3], (0, 0.0, 1.0))[2] for s in sorted_syms], dtype=np.float32)
    n_cells_arr = pb_cnt[sel_pi].astype(np.float32)

    cnt_f = pb_cnt[sel_pi, None].astype(np.float64)
    delta_lib = (pb_sum_lib[sel_pi] / cnt_f).astype(np.float32)
    delta_log = (pb_sum_log[sel_pi] / cnt_f).astype(np.float32)
    delta_z = (pb_sum_z[sel_pi] / cnt_f).astype(np.float32)

    ha_c = ha_pb_cnt[sel_pi]
    hb_c = hb_pb_cnt[sel_pi]
    half_valid = (ha_c >= 1) & (hb_c >= 1)
    ha_den = np.maximum(ha_c[:, None], 1).astype(np.float64)
    hb_den = np.maximum(hb_c[:, None], 1).astype(np.float64)
    ha_d_lib = (ha_pb_sum_lib[sel_pi] / ha_den).astype(np.float32)
    ha_d_log = (ha_pb_sum_log[sel_pi] / ha_den).astype(np.float32)
    ha_d_z = (ha_pb_sum_z[sel_pi] / ha_den).astype(np.float32)
    hb_d_lib = (hb_sum_pb_lib[sel_pi] / hb_den).astype(np.float32)
    hb_d_log = (hb_sum_pb_log[sel_pi] / hb_den).astype(np.float32)
    hb_d_z = (hb_sum_pb_z[sel_pi] / hb_den).astype(np.float32)

    # Save outputs
    (out_dir / "genes.txt").write_text("\n".join(sel_gnames) + "\n", encoding="utf-8")
    (out_dir / "targets.txt").write_text("\n".join(sorted_syms) + "\n", encoding="utf-8")

    batch_labels = np.array([str(uniq_batches[batch_idx_per_cell[i]]) for i in ctrl_sample_idx])
    np.savez_compressed(
        out_dir / "ctrl_cells.npz",
        X_log1p=X_ctrl_log,
        X_zscore=X_ctrl_z,
        X_libnorm=X_ctrl_lib,
        ctrl_mean_log1p=ctrl_mean_log[sel_idx].astype(np.float32),
        ctrl_mean_libnorm=ctrl_mean_lib[sel_idx].astype(np.float32),
        n_ctrl_total=np.int64(len(nt_all_idx)),
        batch_labels=batch_labels,
    )
    np.savez_compressed(
        out_dir / "pseudobulk_delta.npz",
        delta_log1p=delta_log,
        delta_zscore=delta_z,
        delta_libnorm=delta_lib,
        half_a_delta_log1p=ha_d_log,
        half_b_delta_log1p=hb_d_log,
        half_a_delta_zscore=ha_d_z,
        half_b_delta_zscore=hb_d_z,
        half_a_delta_libnorm=ha_d_lib,
        half_b_delta_libnorm=hb_d_lib,
        half_valid=half_valid,
        targets=np.array(sorted_syms, dtype=str),
        target_ensembl_ids=target_ens,
        reproducible=reproducible,
        n_cells=n_cells_arr,
        ad_counts=ad_counts_arr,
        energy_test_p_value=ep_arr,
    )

    sym_set = set(sel_gnames)
    line_meta = {
        "line": line_name,
        "input_file": ", ".join(p.name for p, _ in specs),
        "input_files": [p.name for p, _ in specs],
        "input_sha256": input_sha256,
        "input_sha256s": input_sha256s,
        "n_genes_multifile_dropped": int(n_multifile_dropped),
        "n_duplicate_nt_skipped": int(n_dup_nt_skipped),
        "target_symbol_remaps": target_remaps,
        "raw_shape": [int(n_cells - n_dup_nt_skipped), int(n_genes_raw)],
        "n_batches": int(n_batches),
        "n_nt_total": int(len(nt_all_idx)),
        "n_nt_sampled": int(N_ctrl),
        **gene_audit,
        "n_raw_pert_groups": int(n_pk),
        "n_pert_groups_min10_cells": int(n_pert_min10),
        "n_targets_deduplicated": int(P),
        "n_targets_reproducible": int(reproducible.sum()),
        "n_targets_in_line_genes": int(sum(1 for t in sorted_syms if t in sym_set)),
        "n_targets_half_valid": int(half_valid.sum()),
        "elapsed_sec": round(float(time.time() - t0), 1),
    }
    (out_dir / "line_metadata.json").write_text(json.dumps(line_meta, indent=2), encoding="utf-8")
    print(
        f"[DONE {line_name}] G={G}, ctrl={N_ctrl}/{len(nt_all_idx)}, "
        f"targets={P} (rep={int(reproducible.sum())}, in_genes={line_meta['n_targets_in_line_genes']}), "
        f"elapsed={line_meta['elapsed_sec']}s",
        flush=True,
    )
    return out_dir


def nadig_adapter(f: h5py.File) -> Dict[str, object]:
    g_list = read_obs_or_var_col(f["obs"], "gene")
    gid_list = read_obs_or_var_col(f["obs"], "gene_id")
    gt_list = read_obs_or_var_col(f["obs"], "gene_transcript")
    umi = f["obs/UMI_count"][:].astype(np.float64)
    gems = f["obs/gem_group"][:]
    nt_mask = np.array([g == "non-targeting" for g in g_list], dtype=bool)
    pert_keys = []
    pert_meta: Dict[str, Tuple[str, str]] = {}
    for i in range(len(g_list)):
        if nt_mask[i]:
            pert_keys.append("")
        else:
            ens, sym = gid_list[i], g_list[i]
            if not ens.startswith("ENSG"):
                pert_keys.append("")
            else:
                pk = gt_list[i]
                pert_keys.append(pk)
                if pk not in pert_meta:
                    pert_meta[pk] = (ens, sym)
    return {
        "nt_mask": nt_mask,
        "batches": gems,
        "umi": umi,
        "pert_keys": pert_keys,
        "pert_meta": pert_meta,
        "var_gene_ids": read_obs_or_var_col(f["var"], "gene_id"),
        "var_gene_names": read_obs_or_var_col(f["var"], "gene_name"),
    }


def nadig_precomputed_stats(h5ad_path: Path, cache_meta_csv: Path, cache_x_npy: Path, cache_genes_csv: Path):
    """Map exact R7 Addendum-B AD/Energy statistics onto Nadig gene_transcript keys."""
    meta_df = pd.read_csv(cache_meta_csv)
    with h5py.File(h5ad_path, "r") as f:
        g_list = read_obs_or_var_col(f["obs"], "gene")
        gid_list = read_obs_or_var_col(f["obs"], "gene_id")
        gt_list = read_obs_or_var_col(f["obs"], "gene_transcript")
        gems = f["obs/gem_group"][:]
        umi = f["obs/UMI_count"][:].astype(np.float64)
        var_gids = read_obs_or_var_col(f["var"], "gene_id")
        ens_to_gts = defaultdict(lambda: defaultdict(list))
        for i in range(len(g_list)):
            if g_list[i] == "non-targeting":
                continue
            if gid_list[i].startswith("ENSG"):
                ens_to_gts[gid_list[i]][gt_list[i]].append(i)

        # For rare ties in n_cells across transcripts of same ENSG, disambiguate via R7 cache_X correlation
        cache_X = np.load(cache_x_npy)
        cache_genes = pd.read_csv(cache_genes_csv)["gene_id"].tolist()
        gid_pos = {g: i for i, g in enumerate(var_gids)}
        c_col_idx = np.array([gid_pos[g] for g in cache_genes], dtype=np.int64)
        nt_mask = np.array([g == "non-targeting" for g in g_list], dtype=bool)

        stats_map: Dict[str, Tuple[int, float, float]] = {}
        ambig_rows = []
        for row_idx, row in meta_df.iterrows():
            ens = row["perturbation"]
            nc_target = int(row["n_cells"])
            cands = ens_to_gts[ens]
            matches = [gt for gt, rlist in cands.items() if len(rlist) == nc_target]
            if len(matches) == 1:
                stats_map[matches[0]] = (int(row["reproducible"]), float(row["ad_counts"]), float(row["energy_test_p_value"]))
            else:
                ambig_rows.append((row_idx, row, matches if matches else list(cands.keys())))

        if ambig_rows:
            # Compute per-gem mean/std on cache_genes sequentially for the rare ambiguous ENSG
            uniq_g, ginv = np.unique(gems, return_inverse=True)
            n_c_tot = len(g_list)
            nt_sum_c = np.zeros((len(uniq_g), len(c_col_idx)), dtype=np.float64)
            nt_sq_c = np.zeros((len(uniq_g), len(c_col_idx)), dtype=np.float64)
            nt_cnt_c = np.zeros(len(uniq_g), dtype=np.int64)
            X_ds = f["X"]
            for s in range(0, n_c_tot, 8192):
                e = min(s + 8192, n_c_tot)
                m_loc = nt_mask[s:e]
                if not np.any(m_loc):
                    continue
                rows = np.where(m_loc)[0]
                X_nt = X_ds[s:e, :][rows][:, c_col_idx].astype(np.float64)
                X_nt *= (1e4 / np.maximum(umi[s:e][rows], 1.0))[:, None]
                b_loc = ginv[s:e][rows]
                np.add.at(nt_sum_c, b_loc, X_nt)
                np.add.at(nt_sq_c, b_loc, X_nt * X_nt)
                np.add.at(nt_cnt_c, b_loc, 1)
            mu_c = nt_sum_c / np.maximum(nt_cnt_c[:, None], 1)
            var_c = np.maximum((nt_sq_c - nt_cnt_c[:, None] * (mu_c ** 2)) / np.maximum(nt_cnt_c[:, None] - 1, 1), 0.0)
            sd_c = np.maximum(np.sqrt(var_c), 1e-12)
            for row_idx, row, cand_gts in ambig_rows:
                ref_vec = cache_X[row_idx].astype(np.float64)
                best_gt, best_err = cand_gts[0], float("inf")
                for gt in cand_gts:
                    ridx = np.sort(np.asarray(ens_to_gts[row["perturbation"]][gt], dtype=np.int64))
                    Xp = X_ds[ridx, :][:, c_col_idx].astype(np.float64)
                    Xp *= (1e4 / np.maximum(umi[ridx], 1.0))[:, None]
                    bi_p = ginv[ridx]
                    Zp = ((Xp - mu_c[bi_p]) / sd_c[bi_p]).mean(axis=0)
                    err = float(np.mean((Zp - ref_vec) ** 2))
                    if err < best_err:
                        best_err = err
                        best_gt = gt
                stats_map[best_gt] = (int(row["reproducible"]), float(row["ad_counts"]), float(row["energy_test_p_value"]))
    return stats_map


def replogle_precomputed_stats(bulk_h5ad_path: Path) -> Dict[str, Tuple[int, float, float]]:
    """Load official Replogle et al. (2022) per-gene_transcript AD counts and energy p-values."""
    stats_map: Dict[str, Tuple[int, float, float]] = {}
    with h5py.File(bulk_h5ad_path, "r") as fb:
        gt = decode_str_array(fb["obs/gene_transcript"][:])
        ctrl = fb["obs/core_control"][:]
        ad = fb["obs/anderson_darling_counts"][:]
        ep = fb["obs/energy_test_p_value"][:]
        nc = fb["obs/num_cells_filtered"][:]
        for i, pk in enumerate(gt):
            if ctrl[i]:
                continue
            ad_val = float(ad[i]) if np.isfinite(ad[i]) else 0.0
            ep_val = float(ep[i]) if np.isfinite(ep[i]) else 1.0
            nc_val = float(nc[i]) if np.isfinite(nc[i]) else 0.0
            rep = 1 if (nc_val >= 25.0 and ad_val >= 50.0 and ep_val < 0.001) else 0
            stats_map[pk] = (rep, ad_val, ep_val)
    return stats_map


def load_hgnc_lookups(hgnc_tsv: Optional[Path] = None) -> Tuple[set, Dict[str, str], Dict[str, str]]:
    """Load official HGNC approved symbols, ENSG->symbol, and symbol->ENSG maps if table exists."""
    p = hgnc_tsv or Path("/Users/mac/.cache/r9_p2_raw/hgnc_complete_set.txt")
    if not p.exists():
        return set(), {}, {}
    df = pd.read_csv(p, sep="\t", low_memory=False)
    approved = set()
    ens_to_sym: Dict[str, str] = {}
    sym_to_ens: Dict[str, str] = {}
    for row in df.itertuples(index=False):
        if row.status != "Approved":
            continue
        sym = str(row.symbol).strip()
        approved.add(sym)
        ens = row.ensembl_gene_id if isinstance(row.ensembl_gene_id, str) else ""
        if ens.startswith("ENSG"):
            ens_to_sym[ens] = sym
            sym_to_ens[sym] = ens
    return approved, ens_to_sym, sym_to_ens


def kolf_adapter(
    f: h5py.File,
    tag: str = "",
    approved: Optional[set] = None,
    ens_to_sym: Optional[Dict[str, str]] = None,
) -> Dict[str, object]:
    """Adapter for Nourreddine et al. (2026) KOLF2.1J Strong Perturbations h5ads (layers/counts)."""
    if approved is None or ens_to_sym is None:
        approved, ens_to_sym, _ = load_hgnc_lookups()
    perturbed = read_obs_or_var_col(f["obs"], "perturbed")
    gt = read_obs_or_var_col(f["obs"], "gene_target")
    ge = read_obs_or_var_col(f["obs"], "gene_target_ensembl_id")
    pt = read_obs_or_var_col(f["obs"], "perturbation")
    ch = read_obs_or_var_col(f["obs"], "channel")
    umi = f["obs/n_UMI_counts"][:].astype(np.float64)
    nt_mask = np.array([p == "False" or g == "NTC" for p, g in zip(perturbed, gt)], dtype=bool)
    batches = [f"{tag}_{c}" if tag else c for c in ch]
    pert_keys: List[str] = []
    pert_meta: Dict[str, Tuple[str, str]] = {}
    target_remaps: Dict[str, str] = {}
    for i in range(len(gt)):
        if nt_mask[i]:
            pert_keys.append("")
            continue
        raw_sym = gt[i].strip()
        ens = ge[i].strip()
        sym = raw_sym
        if approved and raw_sym not in approved and ens in ens_to_sym:
            sym = ens_to_sym[ens]
            target_remaps[raw_sym] = sym
        pk = f"{tag}_{pt[i]}" if tag else pt[i]
        pert_keys.append(pk)
        if pk not in pert_meta:
            pert_meta[pk] = (ens, sym)
    return {
        "nt_mask": nt_mask,
        "batches": batches,
        "umi": umi,
        "pert_keys": pert_keys,
        "pert_meta": pert_meta,
        "target_remaps": target_remaps,
        "var_gene_ids": read_obs_or_var_col(f["var"], "gene_ids"),
        "var_gene_names": read_obs_or_var_col(f["var"], "_index"),
    }


def vcc_adapter(
    f: h5py.File,
    tag: str = "",
    sym_to_ens: Optional[Dict[str, str]] = None,
    seen_nt_barcodes: Optional[set] = None,
) -> Dict[str, object]:
    """
    Adapter for Arc Virtual Cell Challenge 2025 H1 hESC h5ads (Training / Validation / Test).
    Deduplicates identical `non-targeting` barcodes shared across Train/Val/Test splits and
    preserves native shared `Flex_*` batch identifiers.
    """
    if sym_to_ens is None:
        _, _, sym_to_ens = load_hgnc_lookups()
    tg = read_obs_or_var_col(f["obs"], "target_gene")
    gd = read_obs_or_var_col(f["obs"], "guide_id")
    bt = read_obs_or_var_col(f["obs"], "batch")
    barcodes = read_obs_or_var_col(f["obs"], "_index")
    nt_mask = np.zeros(len(tg), dtype=bool)
    pert_keys: List[str] = []
    pert_meta: Dict[str, Tuple[str, str]] = {}
    n_dup_nt_skipped = 0
    for i in range(len(tg)):
        if tg[i] == "non-targeting":
            bc = barcodes[i]
            if seen_nt_barcodes is not None:
                if bc in seen_nt_barcodes:
                    n_dup_nt_skipped += 1
                    pert_keys.append("")
                    continue
                seen_nt_barcodes.add(bc)
            nt_mask[i] = True
            pert_keys.append("")
            continue
        sym = tg[i].strip()
        if not sym or sym.lower() in ("nan", "none", "na", "non-targeting"):
            pert_keys.append("")
            continue
        pk = f"{tag}_{gd[i]}" if tag else gd[i]
        pert_keys.append(pk)
        if pk not in pert_meta:
            pert_meta[pk] = (sym_to_ens.get(sym, ""), sym)
    var_names = read_obs_or_var_col(f["var"], "_index")
    if "gene_id" in f["var"]:
        raw_gids = read_obs_or_var_col(f["var"], "gene_id")
        var_gids = [raw_gids[i] or sym_to_ens.get(s, "") for i, s in enumerate(var_names)]
    else:
        var_gids = [sym_to_ens.get(s, "") for s in var_names]
    return {
        "nt_mask": nt_mask,
        "batches": bt,
        "umi": None,
        "pert_keys": pert_keys,
        "pert_meta": pert_meta,
        "target_remaps": {},
        "n_duplicate_nt_skipped": int(n_dup_nt_skipped),
        "var_gene_ids": var_gids,
        "var_gene_names": var_names,
    }


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Prepare standardized shared cell-line inputs (A-T1).")
    p.add_argument("--line", required=True, choices=["nadig_hepg2", "nadig_jurkat", "replogle_k562_essential", "replogle_rpe1", "vcc2025_h1", "nourreddine_kolf21j"])
    p.add_argument("--input-dir", type=Path, required=True)
    p.add_argument("--bulk-dir", type=Path, default=None, help="Optional directory containing Replogle normalized bulk h5ads.")
    p.add_argument("--max-ctrl", type=int, default=20000)
    p.add_argument("--max-genes", type=int, default=None)
    p.add_argument("--include-test", action="store_true", help="Include adata_Test.h5ad for vcc2025_h1.")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if args.line == "nadig_hepg2":
        fp = args.input_dir / "GSE264667_hepg2_raw_singlecell_01.h5ad"
        cdir = args.input_dir / "cache"
        process_singlecell_line(
            "nadig_hepg2",
            fp,
            nadig_adapter,
            max_ctrl_cells=args.max_ctrl,
            max_genes=args.max_genes,
            precomputed_stats_fn=lambda: nadig_precomputed_stats(
                fp, cdir / "HepG2_meta.csv", cdir / "HepG2_X.npy", cdir / "HepG2_genes.csv"
            ),
        )
    elif args.line == "nadig_jurkat":
        fp = args.input_dir / "GSE264667_jurkat_raw_singlecell_01.h5ad"
        cdir = args.input_dir / "cache"
        process_singlecell_line(
            "nadig_jurkat",
            fp,
            nadig_adapter,
            max_ctrl_cells=args.max_ctrl,
            max_genes=args.max_genes,
            precomputed_stats_fn=lambda: nadig_precomputed_stats(
                fp, cdir / "Jurkat_meta.csv", cdir / "Jurkat_X.npy", cdir / "Jurkat_genes.csv"
            ),
        )
    elif args.line == "replogle_k562_essential":
        fp = args.input_dir / "K562_essential_raw_singlecell_01.h5ad"
        bdir = args.bulk_dir or args.input_dir
        bk = bdir / "K562_essential_normalized_bulk_01.h5ad"
        process_singlecell_line(
            "replogle_k562_essential",
            fp,
            nadig_adapter,
            max_ctrl_cells=args.max_ctrl,
            max_genes=args.max_genes,
            precomputed_stats_fn=(lambda: replogle_precomputed_stats(bk)) if bk.exists() else None,
        )
    elif args.line == "replogle_rpe1":
        fp = args.input_dir / "rpe1_raw_singlecell_01.h5ad"
        bdir = args.bulk_dir or args.input_dir
        bk = bdir / "rpe1_normalized_bulk_01.h5ad"
        process_singlecell_line(
            "replogle_rpe1",
            fp,
            nadig_adapter,
            max_ctrl_cells=args.max_ctrl,
            max_genes=args.max_genes,
            precomputed_stats_fn=(lambda: replogle_precomputed_stats(bk)) if bk.exists() else None,
        )
    elif args.line == "nourreddine_kolf21j":
        approved, ens_to_sym, _ = load_hgnc_lookups()
        specs = [
            (args.input_dir / "KOLF_Chromatin_Modifiers_Strong_Perturbations.h5ad", "CM"),
            (args.input_dir / "KOLF_Metabolic_Enzymes_Strong_Perturbations.h5ad", "ME"),
        ]
        keep_syms = collect_existing_line_symbols("nourreddine_kolf21j")
        max_g = args.max_genes if args.max_genes is not None else 15000
        process_singlecell_line(
            "nourreddine_kolf21j",
            specs,
            lambda f, tag: kolf_adapter(f, tag=tag, approved=approved, ens_to_sym=ens_to_sym),
            max_ctrl_cells=args.max_ctrl,
            max_genes=max_g,
            extra_keep_symbols=keep_syms,
            x_layer="counts",
        )
    elif args.line == "vcc2025_h1":
        _, _, sym_to_ens = load_hgnc_lookups()
        specs = [
            (args.input_dir / "adata_Training.h5ad", "train"),
            (args.input_dir / "adata_Validation.h5ad", "val"),
        ]
        test_fp = args.input_dir / "adata_Test.h5ad"
        if args.include_test:
            if not test_fp.exists():
                raise FileNotFoundError(f"Requested --include-test but missing {test_fp}")
            specs.append((test_fp, "test"))
        keep_syms = collect_existing_line_symbols("vcc2025_h1")
        max_g = args.max_genes if args.max_genes is not None else 15000
        seen_nt: set = set()
        process_singlecell_line(
            "vcc2025_h1",
            specs,
            lambda f, tag: vcc_adapter(f, tag=tag, sym_to_ens=sym_to_ens, seen_nt_barcodes=seen_nt),
            max_ctrl_cells=args.max_ctrl,
            max_genes=max_g,
            extra_keep_symbols=keep_syms,
            x_layer="X",
        )
    else:
        raise ValueError(f"Unsupported line: {args.line}")
