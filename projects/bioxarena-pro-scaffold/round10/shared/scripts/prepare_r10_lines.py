#!/usr/bin/env python3
"""
Unified Round 10 Single-Cell Preprocessing & Addendum-B QC Engine (A-T1).

Implements PREREG_C §2.2-2.4 and R10_DOWNLOAD_PLAN.md:
- Disjoint Half A / Half B batch partitioning (per-half NT control centering/SD,
  requiring `>= 2` cells in both Half A and Half B for `half_valid`).
- Canonical HGNC symbol resolution via frozen `hgnc_complete_set.txt`
  (`approved -> ENSG -> unambiguous prev_symbol -> unambiguous alias_symbol`).
- Strictly single-sgRNA perturbed and NTC cells (`P4`).
- Full Addendum-B reproducibility tests (top-5000 finite genes AD BH<0.05 >= 50
  and 20-PC Energy Distance p < 0.001 over 2,000 perms, n_cells >= 25),
  accelerated via pre-sorted Timsort AD + optional multi-worker sharding and
  out-of-core sequential target-group sharding for multi-million-cell cohorts.
- Saves `S7` covariates (`n_cells`, `knockdown_depth_log1p`, `knockdown_fraction`,
  `half_a_n_cells`, `half_b_n_cells`, `target_raw_symbols`, `target_ensembl_ids`,
  `target_hgnc_status`) and full/half control means.
"""

from __future__ import annotations

import argparse
import gc
import json
import multiprocessing as mp
import os
import re
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from statsmodels.stats.multitest import multipletests

PROJECT_ROOT = Path(__file__).resolve().parents[3]
SRC_DIR = PROJECT_ROOT / "src"
R9_SCRIPTS = PROJECT_ROOT / "round9" / "p2" / "shared" / "scripts"
for p in (SRC_DIR, R9_SCRIPTS):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from prepare_nadig2025_cache import (
    ad_two_sample_pvalues_vectorized,
    energy_distance_perm_pvalue,
    sha256_file,
    stratified_nt_sample,
    verify_ad_implementation,
)
from prepare_lines import (
    decode_str_array,
    open_x_node,
    read_h5ad_chunk_dense,
    read_obs_or_var_col,
)

R10_SHARED_DIR = Path(__file__).resolve().parents[1]
R10_RAW_DIR = R10_SHARED_DIR / "raw"
R10_LINES_DIR = R10_SHARED_DIR / "lines"
DEFAULT_HGNC_PATH = Path("/Users/mac/.cache/r9_p2_raw/hgnc_complete_set.txt")
R9_COMMON_GENES_PATH = PROJECT_ROOT / "round9" / "p2" / "shared" / "COMMON_GENES.txt"

DLD1_TARGET_RE = re.compile(r"^([A-Za-z0-9.-]+)(?:_(ENSG\d+))?_(?:P\d+|ENST\d+)")
DLD1_LOW_MOI_CMOS = {f"A{i}" for i in range(2, 13)} | {f"B{i}" for i in range(1, 8)}
CD4_EXCEL_GUIDE_REMAP = {"1-Jun": "JUN-1", "2-Jun": "JUN-2"}


@dataclass(frozen=True)
class HGNCLookup:
    approved: set
    ens_to_sym: Dict[str, str]
    sym_to_ens: Dict[str, str]
    prev_to_sym: Dict[str, str]
    alias_to_sym: Dict[str, str]

    def resolve(self, raw_sym: str, ens_id: str = "") -> Tuple[str, str]:
        """Resolve `(raw_sym, ens_id)` to `(canonical_symbol, resolution_status)`."""
        s = str(raw_sym).strip()
        e = str(ens_id).strip().split(".")[0]
        if s in self.approved:
            return s, "approved"
        if e.startswith("ENSG") and e in self.ens_to_sym:
            return self.ens_to_sym[e], "ensg"
        if s in self.prev_to_sym:
            return self.prev_to_sym[s], "prev_symbol"
        if s in self.alias_to_sym:
            return self.alias_to_sym[s], "alias_symbol"
        if not s and e.startswith("ENSG"):
            return e, "unmapped"
        return s, "unmapped"


def build_hgnc_lookup(hgnc_tsv: Optional[Path] = None) -> HGNCLookup:
    p = hgnc_tsv or DEFAULT_HGNC_PATH
    if not p.exists():
        raise FileNotFoundError(f"Missing frozen HGNC table: {p}")
    df = pd.read_csv(p, sep="\t", low_memory=False)
    app = df[df["status"] == "Approved"]

    approved: set = set()
    ens_to_sym: Dict[str, str] = {}
    ambig_ens: set = set()
    sym_to_ens: Dict[str, str] = {}
    prev_to_sym: Dict[str, str] = {}
    ambig_prev: set = set()
    alias_to_sym: Dict[str, str] = {}
    ambig_alias: set = set()

    for r in app.itertuples(index=False):
        sym = str(r.symbol).strip()
        approved.add(sym)
        ens = str(r.ensembl_gene_id).strip() if isinstance(r.ensembl_gene_id, str) else ""
        if ens.startswith("ENSG"):
            if ens in ens_to_sym and ens_to_sym[ens] != sym:
                ambig_ens.add(ens)
            else:
                ens_to_sym[ens] = sym
            sym_to_ens[sym] = ens
        if isinstance(r.prev_symbol, str):
            for pv in r.prev_symbol.split("|"):
                pv = pv.strip()
                if pv:
                    if pv in prev_to_sym and prev_to_sym[pv] != sym:
                        ambig_prev.add(pv)
                    else:
                        prev_to_sym[pv] = sym
        if isinstance(r.alias_symbol, str):
            for al in r.alias_symbol.split("|"):
                al = al.strip()
                if al:
                    if al in alias_to_sym and alias_to_sym[al] != sym:
                        ambig_alias.add(al)
                    else:
                        alias_to_sym[al] = sym

    for e in ambig_ens:
        ens_to_sym.pop(e, None)
    for pv in ambig_prev:
        prev_to_sym.pop(pv, None)
    for al in ambig_alias:
        alias_to_sym.pop(al, None)

    return HGNCLookup(approved, ens_to_sym, sym_to_ens, prev_to_sym, alias_to_sym)


def select_canonical_genes(
    var_gene_ids: Sequence[str],
    var_gene_names: Sequence[str],
    finite_mask: np.ndarray,
    control_mean: np.ndarray,
    hgnc: HGNCLookup,
    max_genes: Optional[int] = None,
    extra_keep_symbols: Optional[set] = None,
) -> Tuple[np.ndarray, List[str], List[str], Dict[str, object]]:
    """
    Canonicalize finite genes to HGNC symbols, deduplicate collisions by highest
    control_mean (tie-break smallest index), and apply optional `max_genes` cap
    with force-keep set.
    """
    n_raw = len(var_gene_names)
    n_finite = int(finite_mask.sum())
    n_empty = 0
    n_unmapped = 0
    dup_counter: Dict[str, int] = defaultdict(int)
    best_by_sym: Dict[str, Tuple[float, int, str]] = {}

    for idx in np.where(finite_mask)[0]:
        raw_s = str(var_gene_names[idx]).strip()
        raw_e = str(var_gene_ids[idx]).strip() if idx < len(var_gene_ids) else ""
        if not raw_s or raw_s.lower() in ("nan", "none", "na"):
            if not raw_e.startswith("ENSG"):
                n_empty += 1
                continue
        canon_s, status = hgnc.resolve(raw_s, raw_e)
        if status == "unmapped" and canon_s not in hgnc.approved:
            n_unmapped += 1
            if not canon_s:
                n_empty += 1
                continue
        dup_counter[canon_s] += 1
        score = float(control_mean[idx])
        eid = raw_e.split(".")[0] if raw_e.startswith("ENSG") else hgnc.sym_to_ens.get(canon_s, "")
        if canon_s not in best_by_sym or score > best_by_sym[canon_s][0]:
            best_by_sym[canon_s] = (score, int(idx), eid)

    dup_syms = sorted(s for s, c in dup_counter.items() if c > 1)
    n_dup_dropped = sum(c - 1 for c in dup_counter.values() if c > 1)
    kept = list(best_by_sym.items())
    n_pre_cap = len(kept)
    if max_genes is not None and len(kept) > max_genes:
        keep_set = extra_keep_symbols or set()
        forced = [it for it in kept if it[0] in keep_set]
        others = sorted([it for it in kept if it[0] not in keep_set], key=lambda x: x[1][0], reverse=True)
        n_rem = max(0, max_genes - len(forced))
        kept = forced + others[:n_rem]

    kept.sort(key=lambda x: x[1][1])
    sel_idx = np.asarray([it[1][1] for it in kept], dtype=np.int64)
    sel_names = [it[0] for it in kept]
    sel_eids = [it[1][2] for it in kept]
    audit = {
        "n_genes_raw": int(n_raw),
        "n_genes_finite_all_batches": int(n_finite),
        "n_genes_empty_symbol_dropped": int(n_empty),
        "n_genes_unmapped_to_approved_hgnc": int(n_unmapped),
        "n_genes_duplicate_symbol_dropped": int(n_dup_dropped),
        "n_genes_low_expr_cap_dropped": int(n_pre_cap - len(kept)),
        "duplicate_symbols_sample": dup_syms[:20],
        "n_genes_final": int(len(sel_idx)),
    }
    return sel_idx, sel_names, sel_eids, audit


@dataclass
class PartSpec:
    """Virtual cell-matrix partition (one file, channel, or lane)."""
    tag: str
    path: Path
    n_cells: int
    n_genes_raw: int
    nt_mask: np.ndarray            # bool (n_cells,)
    batches: np.ndarray            # str (n_cells,)
    half_flags: np.ndarray         # int8 (n_cells,) 1=Half A, 2=Half B, 0=excluded
    umi: np.ndarray                # float64 (n_cells,)
    pert_keys: List[str]           # len n_cells ("" if NT or excluded)
    pert_meta: Dict[str, Tuple[str, str, str, str]]  # key -> (ens_id, canon_sym, raw_sym, hgnc_status)
    var_gene_ids: List[str]
    var_gene_names: List[str]
    read_chunk_fn: Callable[[int, int], np.ndarray]  # returns dense float64 (e-s, n_genes_raw)
    read_rows_fn: Optional[Callable[[np.ndarray], np.ndarray]] = None


def _compute_half_stats(
    half_sum_lib: np.ndarray,
    half_sq_lib: np.ndarray,
    half_sum_log: np.ndarray,
    half_cnt: np.ndarray,
    batches_in_half: np.ndarray,
    sel_idx: np.ndarray,
    min_batch_nt: int = 15,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute per-batch (mu_lib, sd_lib, mu_log) and global half pooled `ctrl_mean_lib`
    strictly from NTC cells belonging to `batches_in_half`.
    Falls back to the half's pooled NTC SD (or mean if 0 NTC in batch) per PREREG_C §2.3.
    """
    tot_h = int(half_cnt[batches_in_half].sum())
    if tot_h == 0:
        raise RuntimeError("Half has zero NT control cells!")
    g_sum_lib = half_sum_lib[batches_in_half].sum(axis=0)[sel_idx]
    g_sq_lib = half_sq_lib[batches_in_half].sum(axis=0)[sel_idx]
    g_sum_log = half_sum_log[batches_in_half].sum(axis=0)[sel_idx]
    g_mu_lib = (g_sum_lib / tot_h).astype(np.float32)
    g_var_lib = np.maximum((g_sq_lib - tot_h * (g_mu_lib.astype(np.float64) ** 2)) / max(tot_h - 1, 1), 0.0)
    g_sd_lib = np.maximum(np.sqrt(g_var_lib), 1e-12).astype(np.float32)
    g_mu_log = (g_sum_log / tot_h).astype(np.float32)

    has_cnt = half_cnt[:, None] > 0
    cnt_b = np.maximum(half_cnt[:, None], 1)
    raw_mu_lib = (half_sum_lib[:, sel_idx] / cnt_b).astype(np.float32)
    raw_mu_log = (half_sum_log[:, sel_idx] / cnt_b).astype(np.float32)
    mu_lib_b = np.where(has_cnt, raw_mu_lib, g_mu_lib[None, :]).astype(np.float32)
    mu_log_b = np.where(has_cnt, raw_mu_log, g_mu_log[None, :]).astype(np.float32)

    var_b = np.maximum(
        (half_sq_lib[:, sel_idx] - half_cnt[:, None] * (raw_mu_lib.astype(np.float64) ** 2))
        / np.maximum(half_cnt[:, None] - 1, 1),
        0.0,
    )
    sd_b_raw = np.sqrt(var_b).astype(np.float32)
    valid_b = (half_cnt[:, None] >= min_batch_nt) & (sd_b_raw > 1e-12)
    sd_lib_b = np.where(valid_b, sd_b_raw, g_sd_lib[None, :]).astype(np.float32)
    return mu_lib_b, sd_lib_b, mu_log_b, g_mu_lib


# Global worker state initialized once per pool per shard (safe under both spawn & fork)
_W_Z = None
_W_XC_SORTED = None
_W_C_PCA = None
_W_V20 = None
_W_SIGMA_CACHE: Dict[Tuple[int, int], float] = {}


def _ad_worker_init(
    z_mmap_path: Optional[str],
    z_shape: Tuple[int, int],
    xc_sorted: np.ndarray,
    c_pca_en: np.ndarray,
    v20: np.ndarray,
) -> None:
    global _W_Z, _W_XC_SORTED, _W_C_PCA, _W_V20, _W_SIGMA_CACHE
    _W_Z = np.memmap(z_mmap_path, dtype=np.float32, mode="r", shape=z_shape) if z_mmap_path is not None else None
    _W_XC_SORTED = xc_sorted
    _W_C_PCA = c_pca_en
    _W_V20 = v20
    _W_SIGMA_CACHE = {}


def _ad_worker_task(item: Tuple[int, Optional[np.ndarray], Optional[np.ndarray], int]) -> Tuple[int, int, float, float]:
    pi, local_rows, zp_direct, task_seed = item
    Zp = zp_direct if zp_direct is not None else np.asarray(_W_Z[local_rows], dtype=np.float32)
    P_pca = (Zp @ _W_V20).astype(np.float32)
    Zp_s = np.sort(Zp, axis=0)
    p_ad = ad_two_sample_pvalues_vectorized(Zp_s, _W_XC_SORTED, _W_SIGMA_CACHE)
    _, q_ad, _, _ = multipletests(p_ad, method="fdr_bh")
    ad_cnt = float(np.sum(q_ad < 0.05))
    rng_k = np.random.RandomState(task_seed)
    ep_val = float(energy_distance_perm_pvalue(P_pca, _W_C_PCA, rng_k, n_perm=2000, batch_size=100))
    rep = 1 if (ad_cnt >= 50.0 and ep_val < 0.001) else 0
    return pi, rep, ad_cnt, ep_val


def process_r10_line(
    line_name: str,
    parts: Sequence[PartSpec],
    hgnc: HGNCLookup,
    max_ctrl_cells: int = 20000,
    max_genes: Optional[int] = None,
    chunk_size: int = 8192,
    seed: int = 42,
    out_root: Optional[Path] = None,
    extra_meta: Optional[dict] = None,
    max_ram_z_bytes: int = 2_500_000_000,
    shard_max_bytes: int = 8_000_000_000,
    n_workers: int = 8,
) -> Path:
    """
    Chunked, low-RAM processor for R10 standardized line outputs with disjoint split halves.
    """
    t0 = time.time()
    out_dir = (out_root or R10_LINES_DIR) / line_name
    out_dir.mkdir(parents=True, exist_ok=True)
    if not parts:
        raise ValueError(f"{line_name}: empty parts")

    # 1. Align raw var genes across multi-file parts deterministically
    ref_names = [str(x).strip() for x in parts[0].var_gene_names]
    ref_ids = [str(x).strip() for x in parts[0].var_gene_ids]
    same_var = all(
        [str(x).strip() for x in p.var_gene_names] == ref_names and [str(x).strip() for x in p.var_gene_ids] == ref_ids
        for p in parts[1:]
    )
    if same_var:
        col_maps: List[Optional[np.ndarray]] = [None] * len(parts)
        var_gids = ref_ids
        var_gnames = ref_names
        n_multifile_dropped = 0
    else:
        sym_maps: List[Dict[str, int]] = []
        raw_total = len(set(ref_names))
        for p in parts:
            smap: Dict[str, int] = {}
            for j, (gid, gnm) in enumerate(zip(p.var_gene_ids, p.var_gene_names)):
                key = str(gnm).strip() or str(gid).strip()
                if key and key not in smap:
                    smap[key] = j
            sym_maps.append(smap)
        common_keys = [k for k in sym_maps[0] if all(k in sm for sm in sym_maps[1:])]
        col_maps = [np.array([sm[k] for k in common_keys], dtype=np.int64) for sm in sym_maps]
        var_gnames = common_keys
        var_gids = [parts[0].var_gene_ids[j] for j in col_maps[0]]
        n_multifile_dropped = int(raw_total - len(common_keys))

    n_genes_raw = len(var_gnames)

    # Concatenate lightweight per-cell index vectors across parts
    part_offsets = []
    cur = 0
    all_nt: List[np.ndarray] = []
    all_batches: List[str] = []
    all_half: List[np.ndarray] = []
    all_umi: List[np.ndarray] = []
    all_pk: List[str] = []
    pert_meta: Dict[str, Tuple[str, str, str, str]] = {}
    file_sha256s: Dict[str, str] = {}

    for p in parts:
        part_offsets.append((cur, cur + p.n_cells))
        cur += p.n_cells
        all_nt.append(p.nt_mask)
        all_batches.extend(str(b) for b in p.batches)
        all_half.append(p.half_flags)
        all_umi.append(p.umi.astype(np.float64))
        all_pk.extend(p.pert_keys)
        pert_meta.update(p.pert_meta)
        if p.path.name not in file_sha256s and p.path.exists():
            file_sha256s[p.path.name] = sha256_file(p.path)

    n_cells = cur
    nt_mask = np.concatenate(all_nt)
    half_arr = np.concatenate(all_half)
    umi = np.concatenate(all_umi)
    uniq_batches, batch_idx_per_cell = np.unique(np.array(all_batches, dtype=str), return_inverse=True)
    n_batches = len(uniq_batches)

    half_a_batches = np.where([np.any((batch_idx_per_cell == b) & (half_arr == 1)) for b in range(n_batches)])[0]
    half_b_batches = np.where([np.any((batch_idx_per_cell == b) & (half_arr == 2)) for b in range(n_batches)])[0]

    def iter_chunks():
        for p, cmap, (s_g, e_g) in zip(parts, col_maps, part_offsets):
            for loc_s in range(0, p.n_cells, chunk_size):
                loc_e = min(loc_s + chunk_size, p.n_cells)
                yield s_g + loc_s, s_g + loc_e, p, cmap, loc_s, loc_e

    # 2. Pass 1: NT batch stats (overall + Half A + Half B)
    sum_lib = np.zeros((n_batches, n_genes_raw), dtype=np.float64)
    sq_lib = np.zeros((n_batches, n_genes_raw), dtype=np.float64)
    sum_log = np.zeros((n_batches, n_genes_raw), dtype=np.float64)
    nt_cnt = np.zeros(n_batches, dtype=np.int64)

    ha_sum_lib = np.zeros((n_batches, n_genes_raw), dtype=np.float64)
    ha_sq_lib = np.zeros((n_batches, n_genes_raw), dtype=np.float64)
    ha_sum_log = np.zeros((n_batches, n_genes_raw), dtype=np.float64)
    ha_cnt = np.zeros(n_batches, dtype=np.int64)

    hb_sum_lib = np.zeros((n_batches, n_genes_raw), dtype=np.float64)
    hb_sq_lib = np.zeros((n_batches, n_genes_raw), dtype=np.float64)
    hb_sum_log = np.zeros((n_batches, n_genes_raw), dtype=np.float64)
    hb_cnt = np.zeros(n_batches, dtype=np.int64)

    print(
        f"\n=== [{line_name}] Pass 1 (NTC stats over {int(nt_mask.sum()):,} NTCs / {n_cells:,} cells, {n_batches} batches) ===",
        flush=True,
    )
    for g_s, g_e, p, cmap, loc_s, loc_e in iter_chunks():
        m_nt = nt_mask[g_s:g_e]
        if not np.any(m_nt):
            continue
        rows = np.where(m_nt)[0]
        if p.read_rows_fn is not None:
            X_nt = p.read_rows_fn(loc_s + rows).astype(np.float64, copy=False)
        else:
            X_nt = p.read_chunk_fn(loc_s, loc_e)[rows].astype(np.float64, copy=False)
        if cmap is not None:
            X_nt = X_nt[:, cmap]

        scale = (1e4 / np.maximum(umi[g_s:g_e][rows], 1.0))[:, None]
        X_lib = X_nt * scale
        X_log = np.log1p(X_lib)
        bi = batch_idx_per_cell[g_s:g_e][rows]
        hf = half_arr[g_s:g_e][rows]

        for b in np.unique(bi):
            mb = bi == b
            Xlb = X_lib[mb]
            sum_lib[b] += Xlb.sum(axis=0)
            sq_lib[b] += (Xlb * Xlb).sum(axis=0)
            sum_log[b] += X_log[mb].sum(axis=0)
            nt_cnt[b] += int(mb.sum())

            ma = mb & (hf == 1)
            if np.any(ma):
                Xla = X_lib[ma]
                ha_sum_lib[b] += Xla.sum(axis=0)
                ha_sq_lib[b] += (Xla * Xla).sum(axis=0)
                ha_sum_log[b] += X_log[ma].sum(axis=0)
                ha_cnt[b] += int(ma.sum())

            mhb = mb & (hf == 2)
            if np.any(mhb):
                Xlb2 = X_lib[mhb]
                hb_sum_lib[b] += Xlb2.sum(axis=0)
                hb_sq_lib[b] += (Xlb2 * Xlb2).sum(axis=0)
                hb_sum_log[b] += X_log[mhb].sum(axis=0)
                hb_cnt[b] += int(mhb.sum())

    nt_all_idx = np.where(nt_mask)[0]
    n_nt_tot = int(nt_cnt.sum())
    if n_nt_tot == 0:
        raise RuntimeError(f"{line_name}: 0 NTC cells found!")

    ctrl_mean_lib = sum_lib.sum(axis=0) / n_nt_tot
    ctrl_mean_log = sum_log.sum(axis=0) / n_nt_tot
    var_global = np.maximum((sq_lib.sum(axis=0) - n_nt_tot * (ctrl_mean_lib ** 2)) / max(n_nt_tot - 1, 1), 0.0)
    sd_global = np.sqrt(var_global)

    has_nt = nt_cnt[:, None] > 0
    raw_mu_lib = sum_lib / np.maximum(nt_cnt[:, None], 1)
    mu_lib = np.where(has_nt, raw_mu_lib, ctrl_mean_lib[None, :])
    mu_log = np.where(has_nt, sum_log / np.maximum(nt_cnt[:, None], 1), ctrl_mean_log[None, :])
    var_lib = np.maximum((sq_lib - nt_cnt[:, None] * (raw_mu_lib ** 2)) / np.maximum(nt_cnt[:, None] - 1, 1), 0.0)
    sd_lib = np.sqrt(var_lib)

    well_pop = nt_cnt >= 15
    finite_batch_ok = (
        np.all(np.isfinite(mu_lib[well_pop]) & np.isfinite(sd_lib[well_pop]) & (sd_lib[well_pop] > 1e-12), axis=0)
        if np.any(well_pop)
        else np.ones(n_genes_raw, dtype=bool)
    )
    finite_mask = (
        finite_batch_ok
        & np.all(np.isfinite(mu_lib), axis=0)
        & np.isfinite(ctrl_mean_lib)
        & (ctrl_mean_lib > 0.0)
        & np.isfinite(sd_global)
        & (sd_global > 1e-12)
    )
    sd_lib = np.where((nt_cnt[:, None] >= 15) & (sd_lib > 1e-12), sd_lib, sd_global[None, :])

    # Canonical HGNC selection on finite genes, force-retaining canonical R9 COMMON_GENES & targets
    extra_keep = {sym for _, sym, _, _ in pert_meta.values()}
    if R9_COMMON_GENES_PATH.exists():
        for x in R9_COMMON_GENES_PATH.read_text(encoding="utf-8").splitlines():
            s = x.strip()
            if s:
                cs, _ = hgnc.resolve(s)
                extra_keep.add(cs)

    sel_idx, sel_gnames, sel_gids, gene_audit = select_canonical_genes(
        var_gids, var_gnames, finite_mask, ctrl_mean_lib, hgnc, max_genes=max_genes, extra_keep_symbols=extra_keep
    )
    ok_raw_idx = np.where(finite_mask)[0]
    top5000_raw_idx = np.sort(ok_raw_idx[np.argsort(ctrl_mean_lib[ok_raw_idx])[-min(5000, len(ok_raw_idx)) :]])
    G = len(sel_idx)
    print(
        f"  [Pass 1 DONE] finite={gene_audit['n_genes_finite_all_batches']}/{n_genes_raw}, "
        f"G_final={G}, nt_total={n_nt_tot:,} ({time.time()-t0:.1f}s)",
        flush=True,
    )

    mu_lib_s = mu_lib[:, sel_idx].astype(np.float32)
    sd_lib_s = sd_lib[:, sel_idx].astype(np.float32)
    mu_log_s = mu_log[:, sel_idx].astype(np.float32)

    ha_mu_lib, ha_sd_lib, ha_mu_log, ha_ctrl_mean_lib = _compute_half_stats(
        ha_sum_lib, ha_sq_lib, ha_sum_log, ha_cnt, half_a_batches, sel_idx
    )
    hb_mu_lib, hb_sd_lib, hb_mu_log, hb_ctrl_mean_lib = _compute_half_stats(
        hb_sum_lib, hb_sq_lib, hb_sum_log, hb_cnt, half_b_batches, sel_idx
    )

    # 3. Map perturbations & pre-plan sharded Z_5000 storage for eligible groups (n >= 25)
    pkey_to_rows: Dict[str, List[int]] = defaultdict(list)
    for i in range(n_cells):
        if nt_mask[i] or half_arr[i] not in (1, 2):
            continue
        pk = all_pk[i]
        if pk:
            pkey_to_rows[pk].append(i)

    all_pkeys = sorted(pkey_to_rows.keys())
    pkey_to_idx = {pk: i for i, pk in enumerate(all_pkeys)}
    cell_pkey_idx = np.full(n_cells, -1, dtype=np.int32)
    for pk, rlist in pkey_to_rows.items():
        cell_pkey_idx[np.asarray(rlist, dtype=np.int64)] = pkey_to_idx[pk]
    n_pk = len(all_pkeys)

    ctrl_sample_idx = stratified_nt_sample(nt_all_idx, batch_idx_per_cell, n_target=max_ctrl_cells, seed=seed)
    ctrl_pos_map = np.full(n_cells, -1, dtype=np.int32)
    ctrl_pos_map[ctrl_sample_idx] = np.arange(len(ctrl_sample_idx), dtype=np.int32)
    N_ctrl = len(ctrl_sample_idx)

    nt_ad_idx = stratified_nt_sample(nt_all_idx, batch_idx_per_cell, n_target=1000, seed=42)
    nt_en_idx = stratified_nt_sample(nt_all_idx, batch_idx_per_cell, n_target=500, seed=42)
    nt_ad_pos = np.full(n_cells, -1, dtype=np.int32)
    nt_ad_pos[nt_ad_idx] = np.arange(len(nt_ad_idx), dtype=np.int32)
    nt_en_pos = np.full(n_cells, -1, dtype=np.int32)
    nt_en_pos[nt_en_idx] = np.arange(len(nt_en_idx), dtype=np.int32)

    X_ctrl_lib = np.empty((N_ctrl, G), dtype=np.float32)
    X_ctrl_log = np.empty((N_ctrl, G), dtype=np.float32)
    X_ctrl_z = np.empty((N_ctrl, G), dtype=np.float32)

    G_top = len(top5000_raw_idx)
    X_ctrl_ad = np.empty((len(nt_ad_idx), G_top), dtype=np.float32)
    X_ctrl_en_5000 = np.empty((len(nt_en_idx), G_top), dtype=np.float32)

    pb_cnt = np.zeros(n_pk, dtype=np.int64)
    pb_sum_lib = np.zeros((n_pk, G), dtype=np.float64)
    pb_sum_log = np.zeros((n_pk, G), dtype=np.float64)
    pb_sum_z = np.zeros((n_pk, G), dtype=np.float64)

    ha_pb_cnt = np.zeros(n_pk, dtype=np.int64)
    ha_pb_sum_lib = np.zeros((n_pk, G), dtype=np.float64)
    ha_pb_sum_log = np.zeros((n_pk, G), dtype=np.float64)
    ha_pb_sum_z = np.zeros((n_pk, G), dtype=np.float64)

    hb_pb_cnt = np.zeros(n_pk, dtype=np.int64)
    hb_pb_sum_lib = np.zeros((n_pk, G), dtype=np.float64)
    hb_pb_sum_log = np.zeros((n_pk, G), dtype=np.float64)
    hb_pb_sum_z = np.zeros((n_pk, G), dtype=np.float64)

    elig_pidx = [i for i, pk in enumerate(all_pkeys) if len(pkey_to_rows[pk]) >= 25]
    n_elig_cells = sum(len(pkey_to_rows[all_pkeys[pi]]) for pi in elig_pidx)
    row_bytes = max(G_top, 1) * 4
    use_ram_z = (n_elig_cells * row_bytes) <= max_ram_z_bytes

    z_ram_pos = np.full(n_cells, -1, dtype=np.int32)
    cell_shard = np.full(n_cells, -1, dtype=np.int32)
    shard_groups: List[List[int]] = []
    shard_counts: List[int] = []
    Z_ram: Optional[np.ndarray] = None

    if use_ram_z:
        Z_ram = np.empty((n_elig_cells, G_top), dtype=np.float32)
        pos_cur = 0
        for pi in elig_pidx:
            r_arr = np.asarray(pkey_to_rows[all_pkeys[pi]], dtype=np.int64)
            z_ram_pos[r_arr] = np.arange(pos_cur, pos_cur + len(r_arr), dtype=np.int32)
            pos_cur += len(r_arr)
        if elig_pidx:
            shard_groups = [elig_pidx]
            shard_counts = [n_elig_cells]
    else:
        shard_cap_rows = max(1, int(shard_max_bytes // row_bytes))
        cur_grp: List[int] = []
        cur_rows = 0
        for pi in elig_pidx:
            nr = len(pkey_to_rows[all_pkeys[pi]])
            if cur_grp and cur_rows + nr > shard_cap_rows:
                shard_groups.append(cur_grp)
                shard_counts.append(cur_rows)
                cur_grp = []
                cur_rows = 0
            sid = len(shard_groups)
            r_arr = np.asarray(pkey_to_rows[all_pkeys[pi]], dtype=np.int64)
            cell_shard[r_arr] = sid
            cur_grp.append(pi)
            cur_rows += nr
        if cur_grp:
            shard_groups.append(cur_grp)
            shard_counts.append(cur_rows)

    keep_all = nt_mask | (cell_pkey_idx >= 0)
    mu_top = mu_lib[:, top5000_raw_idx].astype(np.float32)
    sd_top = sd_lib[:, top5000_raw_idx].astype(np.float32)
    ZtZ = np.zeros((G_top, G_top), dtype=np.float64)

    shard_row_paths = [out_dir / f"_tmp_shard_{s}.rows.bin" for s in range(len(shard_groups))]
    shard_pid_paths = [out_dir / f"_tmp_shard_{s}.pidx.bin" for s in range(len(shard_groups))]

    t2 = time.time()
    mode_desc = f"RAM({n_elig_cells:,}x{G_top})" if use_ram_z else f"{len(shard_groups)} SSD shards (~{n_elig_cells*row_bytes/1e9:.1f}GB)"
    print(
        f"  [Pass 2] Populating ctrl_cells ({N_ctrl:,}x{G}), pseudobulk ({n_pk} groups, elig={len(elig_pidx)}), Z_5000={mode_desc}...",
        flush=True,
    )

    try:
        shard_row_fps = [] if use_ram_z else [open(pth, "wb") for pth in shard_row_paths]
        shard_pid_fps = [] if use_ram_z else [open(pth, "wb") for pth in shard_pid_paths]
        try:
            for g_s, g_e, p, cmap, loc_s, loc_e in iter_chunks():
                m_keep = keep_all[g_s:g_e]
                if not np.any(m_keep):
                    continue
                Xc = p.read_chunk_fn(loc_s, loc_e)
                if cmap is not None:
                    Xc = Xc[:, cmap]
                scale = (1e4 / np.maximum(umi[g_s:g_e], 1.0)).astype(np.float32)
                bi = batch_idx_per_cell[g_s:g_e]

                # Top-5000 Z for keep_all -> ZtZ Gram accumulation + NT / eligible cell Z_5000 routing
                k_loc = np.where(m_keep)[0]
                Zt_k = Xc[k_loc][:, top5000_raw_idx].astype(np.float32)
                Zt_k *= scale[k_loc][:, None]
                Zt_k -= mu_top[bi[k_loc]]
                Zt_k /= sd_top[bi[k_loc]]
                Zt_k64 = Zt_k.astype(np.float64)
                ZtZ += Zt_k64.T @ Zt_k64

                # Map global chunk indices -> position in Zt_k
                loc_to_zpos = np.full(g_e - g_s, -1, dtype=np.int32)
                loc_to_zpos[k_loc] = np.arange(len(k_loc), dtype=np.int32)

                ad_loc = np.where(nt_ad_pos[g_s:g_e] >= 0)[0]
                if len(ad_loc) > 0:
                    X_ctrl_ad[nt_ad_pos[g_s:g_e][ad_loc]] = Zt_k[loc_to_zpos[ad_loc]]

                en_loc = np.where(nt_en_pos[g_s:g_e] >= 0)[0]
                if len(en_loc) > 0:
                    X_ctrl_en_5000[nt_en_pos[g_s:g_e][en_loc]] = Zt_k[loc_to_zpos[en_loc]]

                if use_ram_z:
                    zr_loc = np.where(z_ram_pos[g_s:g_e] >= 0)[0]
                    if len(zr_loc) > 0:
                        Z_ram[z_ram_pos[g_s:g_e][zr_loc]] = Zt_k[loc_to_zpos[zr_loc]]
                else:
                    sh_chunk = cell_shard[g_s:g_e]
                    for sid in np.unique(sh_chunk[sh_chunk >= 0]):
                        sl = np.where(sh_chunk == sid)[0]
                        Zt_k[loc_to_zpos[sl]].tofile(shard_row_fps[sid])
                        cell_pkey_idx[g_s:g_e][sl].astype(np.int32).tofile(shard_pid_fps[sid])

                # Canonical G matrix -> ctrl_cells + overall + Half A + Half B pseudobulks
                Xs_lib = Xc[:, sel_idx].astype(np.float32) * scale[:, None]
                Xs_log = np.log1p(Xs_lib)
                Xs_z = (Xs_lib - mu_lib_s[bi]) / sd_lib_s[bi]

                cpos = ctrl_pos_map[g_s:g_e]
                c_loc = np.where(cpos >= 0)[0]
                if len(c_loc) > 0:
                    dst = cpos[c_loc]
                    X_ctrl_lib[dst] = Xs_lib[c_loc]
                    X_ctrl_log[dst] = Xs_log[c_loc]
                    X_ctrl_z[dst] = Xs_z[c_loc]

                pidx = cell_pkey_idx[g_s:g_e]
                p_loc = np.where(pidx >= 0)[0]
                if len(p_loc) > 0:
                    pi_loc = pidx[p_loc]
                    bi_loc = bi[p_loc]
                    np.add.at(pb_sum_lib, pi_loc, Xs_lib[p_loc] - mu_lib_s[bi_loc])
                    np.add.at(pb_sum_log, pi_loc, Xs_log[p_loc] - mu_log_s[bi_loc])
                    np.add.at(pb_sum_z, pi_loc, Xs_z[p_loc])
                    np.add.at(pb_cnt, pi_loc, 1)

                    hf = half_arr[g_s:g_e][p_loc]
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
                        np.add.at(hb_pb_sum_lib, pi_b, Xs_lib[loc_b] - hb_mu_lib[bi_b])
                        np.add.at(hb_pb_sum_log, pi_b, Xs_log[loc_b] - hb_mu_log[bi_b])
                        np.add.at(hb_pb_sum_z, pi_b, (Xs_lib[loc_b] - hb_mu_lib[bi_b]) / hb_sd_lib[bi_b])
                        np.add.at(hb_pb_cnt, pi_b, 1)
        finally:
            for fp in shard_row_fps + shard_pid_fps:
                fp.close()

        print(f"  [Pass 2 DONE] ({time.time()-t2:.1f}s). Computing exact SVD eigenvectors & Addendum-B QC...", flush=True)
        verify_ad_implementation()
        k_pc = min(20, G_top)
        _, evecs = np.linalg.eigh(ZtZ)
        V20 = evecs[:, -k_pc:][:, ::-1].astype(np.float32)
        Xc_sorted = np.sort(X_ctrl_ad, axis=0)
        C_pca_en = (X_ctrl_en_5000 @ V20).astype(np.float32)

        stats_map: Dict[str, Tuple[int, float, float]] = {}
        for k_i, pk in enumerate(all_pkeys):
            if len(pkey_to_rows[pk]) < 25:
                stats_map[pk] = (0, 0.0, 1.0)

        done_elig = 0
        eff_workers = min(max(1, int(n_workers)), os.cpu_count() or 1)
        for sid, grp_pis in enumerate(shard_groups):
            if use_ram_z:
                z_mmap_path = None
                z_shape = (0, 0)
                tasks = [
                    (
                        pi,
                        None,
                        Z_ram[z_ram_pos[np.asarray(pkey_to_rows[all_pkeys[pi]], dtype=np.int64)]],
                        (seed + pi * 1000003) & 0x7FFFFFFF,
                    )
                    for pi in grp_pis
                ]
            else:
                n_s = shard_counts[sid]
                z_mmap_path = str(shard_row_paths[sid])
                z_shape = (n_s, G_top)
                pid_s = np.fromfile(shard_pid_paths[sid], dtype=np.int32)
                if len(pid_s) != n_s:
                    raise RuntimeError(f"Shard {sid} size mismatch: {len(pid_s)} != {n_s}")
                order_s = np.argsort(pid_s, kind="mergesort")
                sorted_pids = pid_s[order_s]
                uniq_p, starts, counts = np.unique(sorted_pids, return_index=True, return_counts=True)
                offset_map = {
                    int(p_val): order_s[int(st) : int(st + cnt)]
                    for p_val, st, cnt in zip(uniq_p, starts, counts)
                }
                tasks = [
                    (pi, offset_map[pi], None, (seed + pi * 1000003) & 0x7FFFFFFF)
                    for pi in grp_pis
                ]

            if eff_workers > 1 and len(tasks) >= 16:
                ctx = mp.get_context("spawn")
                with ProcessPoolExecutor(
                    max_workers=eff_workers,
                    mp_context=ctx,
                    initializer=_ad_worker_init,
                    initargs=(z_mmap_path, z_shape, Xc_sorted, C_pca_en, V20),
                ) as ex:
                    for pi_res, rep, ad_cnt, ep_val in ex.map(_ad_worker_task, tasks, chunksize=16):
                        stats_map[all_pkeys[pi_res]] = (rep, ad_cnt, ep_val)
            else:
                _ad_worker_init(z_mmap_path, z_shape, Xc_sorted, C_pca_en, V20)
                for item in tasks:
                    pi_res, rep, ad_cnt, ep_val = _ad_worker_task(item)
                    stats_map[all_pkeys[pi_res]] = (rep, ad_cnt, ep_val)
                _ad_worker_init(None, (0, 0), np.empty((0, 0), np.float32), np.empty((0, 0), np.float32), np.empty((0, 0), np.float32))

            done_elig += len(grp_pis)
            rep_cnt = sum(v[0] for v in stats_map.values())
            print(
                f"    [Shard {sid+1}/{len(shard_groups)}] evaluated {done_elig}/{len(elig_pidx)}, rep_so_far={rep_cnt} ({time.time()-t2:.1f}s)",
                flush=True,
            )
            if not use_ram_z:
                shard_row_paths[sid].unlink(missing_ok=True)
                shard_pid_paths[sid].unlink(missing_ok=True)

        del Z_ram, ZtZ, V20, X_ctrl_ad, X_ctrl_en_5000, Xc_sorted, C_pca_en
        gc.collect()
    finally:
        for pth in shard_row_paths + shard_pid_paths:
            if pth.exists():
                pth.unlink(missing_ok=True)

    # 4. Deduplicate to 1 row per unique canonical target gene symbol (retain n_cells >= 10)
    best_by_sym: Dict[str, Tuple[Tuple[int, float, float], int, str, str, str]] = {}
    n_pert_min10 = 0
    for pk in all_pkeys:
        pi = pkey_to_idx[pk]
        n_c = int(pb_cnt[pi])
        if n_c < 10:
            continue
        n_pert_min10 += 1
        ens, sym, raw_sym, hstat = pert_meta[pk]
        sym = str(sym).strip()
        if not sym or sym.lower() in ("nan", "none", "na", "non-targeting"):
            continue
        rep, ad_val, ep_val = stats_map[pk]
        rep = 1 if (n_c >= 25 and ad_val >= 50.0 and ep_val < 0.001) else 0
        key_score = (int(rep), float(ad_val), float(n_c))
        if sym not in best_by_sym or key_score > best_by_sym[sym][0]:
            best_by_sym[sym] = (key_score, pi, ens, raw_sym, hstat)

    sorted_syms = sorted(best_by_sym.keys())
    P = len(sorted_syms)
    sel_pi = np.array([best_by_sym[s][1] for s in sorted_syms], dtype=np.int64)
    target_ens = np.array([best_by_sym[s][2] for s in sorted_syms], dtype=str)
    target_raw = np.array([best_by_sym[s][3] for s in sorted_syms], dtype=str)
    target_hstat = np.array([best_by_sym[s][4] for s in sorted_syms], dtype=str)
    reproducible = np.array([best_by_sym[s][0][0] == 1 for s in sorted_syms], dtype=bool)
    ad_counts_arr = np.array([best_by_sym[s][0][1] for s in sorted_syms], dtype=np.float32)
    ep_arr = np.array([stats_map[all_pkeys[pi]][2] for pi in sel_pi], dtype=np.float32)
    n_cells_arr = pb_cnt[sel_pi].astype(np.float32)

    cnt_f = np.maximum(pb_cnt[sel_pi, None], 1).astype(np.float64)
    delta_lib = (pb_sum_lib[sel_pi] / cnt_f).astype(np.float32)
    delta_log = (pb_sum_log[sel_pi] / cnt_f).astype(np.float32)
    delta_z = (pb_sum_z[sel_pi] / cnt_f).astype(np.float32)

    ha_c = ha_pb_cnt[sel_pi]
    hb_c = hb_pb_cnt[sel_pi]
    half_valid = (ha_c >= 2) & (hb_c >= 2)  # PREREG_C §2.3
    ha_den = np.maximum(ha_c[:, None], 1).astype(np.float64)
    hb_den = np.maximum(hb_c[:, None], 1).astype(np.float64)
    ha_d_lib = (ha_pb_sum_lib[sel_pi] / ha_den).astype(np.float32)
    ha_d_log = (ha_pb_sum_log[sel_pi] / ha_den).astype(np.float32)
    ha_d_z = (ha_pb_sum_z[sel_pi] / ha_den).astype(np.float32)
    hb_d_lib = (hb_pb_sum_lib[sel_pi] / hb_den).astype(np.float32)
    hb_d_log = (hb_pb_sum_log[sel_pi] / hb_den).astype(np.float32)
    hb_d_z = (hb_pb_sum_z[sel_pi] / hb_den).astype(np.float32)

    gname_to_pos = {g: j for j, g in enumerate(sel_gnames)}
    kd_log1p = np.full(P, np.nan, dtype=np.float32)
    kd_frac = np.full(P, np.nan, dtype=np.float32)
    cmean_lib_s = ctrl_mean_lib[sel_idx].astype(np.float32)
    for row_i, sym in enumerate(sorted_syms):
        g_j = gname_to_pos.get(sym)
        if g_j is not None:
            kd_log1p[row_i] = delta_log[row_i, g_j]
            base_c = float(cmean_lib_s[g_j])
            kd_frac[row_i] = float(-delta_lib[row_i, g_j] / max(base_c, 1e-8))

    (out_dir / "genes.txt").write_text("\n".join(sel_gnames) + "\n", encoding="utf-8")
    (out_dir / "targets.txt").write_text("\n".join(sorted_syms) + "\n", encoding="utf-8")

    batch_labels = np.array([str(uniq_batches[batch_idx_per_cell[i]]) for i in ctrl_sample_idx])
    ctrl_mean_log_s = ctrl_mean_log[sel_idx].astype(np.float32)
    np.savez_compressed(
        out_dir / "ctrl_cells.npz",
        X_log1p=X_ctrl_log,
        X_zscore=X_ctrl_z,
        X_libnorm=X_ctrl_lib,
        ctrl_mean_log1p=ctrl_mean_log_s,
        ctrl_mean_libnorm=cmean_lib_s,
        ctrl_mean_log1p_full=ctrl_mean_log_s,
        ctrl_mean_libnorm_full=cmean_lib_s,
        ctrl_mean_log1p_sub20k=X_ctrl_log.mean(axis=0).astype(np.float32),
        ctrl_mean_libnorm_sub20k=X_ctrl_lib.mean(axis=0).astype(np.float32),
        half_a_ctrl_mean_libnorm=ha_ctrl_mean_lib,
        half_b_ctrl_mean_libnorm=hb_ctrl_mean_lib,
        n_ctrl_total=np.int64(n_nt_tot),
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
        half_a_n_cells=ha_c.astype(np.float32),
        half_b_n_cells=hb_c.astype(np.float32),
        half_valid=half_valid,
        targets=np.array(sorted_syms, dtype=str),
        target_raw_symbols=target_raw,
        target_ensembl_ids=target_ens,
        target_hgnc_status=target_hstat,
        reproducible=reproducible,
        n_cells=n_cells_arr,
        knockdown_depth_log1p=kd_log1p,
        knockdown_fraction=kd_frac,
        ad_counts=ad_counts_arr,
        energy_test_p_value=ep_arr,
    )

    sym_set = set(sel_gnames)
    line_meta = {
        "line": line_name,
        "n_input_files": len(file_sha256s),
        "input_sha256s": file_sha256s,
        "n_genes_multifile_dropped": int(n_multifile_dropped),
        "raw_shape": [int(n_cells), int(n_genes_raw)],
        "n_batches": int(n_batches),
        "n_half_a_batches": int(len(half_a_batches)),
        "n_half_b_batches": int(len(half_b_batches)),
        "n_nt_total": int(n_nt_tot),
        "n_nt_sampled": int(N_ctrl),
        "n_nt_half_a": int(ha_cnt.sum()),
        "n_nt_half_b": int(hb_cnt.sum()),
        **gene_audit,
        "n_raw_pert_groups": int(n_pk),
        "n_pert_groups_min10_cells": int(n_pert_min10),
        "n_targets_deduplicated": int(P),
        "n_targets_reproducible": int(reproducible.sum()),
        "n_targets_in_line_genes": int(sum(1 for t in sorted_syms if t in sym_set)),
        "n_targets_half_valid": int(half_valid.sum()),
        "elapsed_sec": round(float(time.time() - t0), 1),
        **(extra_meta or {}),
    }
    (out_dir / "line_metadata.json").write_text(json.dumps(line_meta, indent=2), encoding="utf-8")
    print(
        f"[DONE {line_name}] G={G}, ctrl={N_ctrl}/{n_nt_tot}, targets={P} "
        f"(rep={int(reproducible.sum())}, half_valid={int(half_valid.sum())})",
        flush=True,
    )
    return out_dir


# ==============================================================================
# Shared 10x CSC .h5 Utilities
# ==============================================================================


def compute_10x_csc_rna_qc(
    h5_path: Path,
    rna_row_idx: Optional[np.ndarray],
    mt_mask_in_raw: np.ndarray,
    batch_cols: int = 20000,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute `(umi, ngenes, pct_mt)` directly from 10x CSC `.h5` without densifying."""
    with h5py.File(h5_path, "r") as f:
        m = f["matrix"]
        n_feat, n_col = int(m["shape"][0]), int(m["shape"][1])
        indptr = m["indptr"][:]
        umi = np.zeros(n_col, dtype=np.float64)
        ngenes = np.zeros(n_col, dtype=np.int32)
        mt_umi = np.zeros(n_col, dtype=np.float64)

        is_rna_feat = None
        if rna_row_idx is not None and len(rna_row_idx) < n_feat:
            is_rna_feat = np.zeros(n_feat, dtype=bool)
            is_rna_feat[rna_row_idx] = True

        for cs in range(0, n_col, batch_cols):
            ce = min(cs + batch_cols, n_col)
            p0, p1 = int(indptr[cs]), int(indptr[ce])
            if p1 == p0:
                continue
            dat = m["data"][p0:p1].astype(np.float64, copy=False)
            idx = m["indices"][p0:p1]
            col_lens = np.diff(indptr[cs : ce + 1])
            col_ids = np.repeat(np.arange(cs, ce, dtype=np.int32), col_lens)
            if is_rna_feat is not None:
                m_r = is_rna_feat[idx]
                dat = dat[m_r]
                idx = idx[m_r]
                col_ids = col_ids[m_r]
            if len(dat) == 0:
                continue
            umi_seg = np.bincount(col_ids - cs, weights=dat, minlength=ce - cs)
            ng_seg = np.bincount(col_ids - cs, minlength=ce - cs)
            mt_m = mt_mask_in_raw[idx]
            mt_seg = np.bincount(col_ids[mt_m] - cs, weights=dat[mt_m], minlength=ce - cs) if np.any(mt_m) else 0.0
            umi[cs:ce] = umi_seg
            ngenes[cs:ce] = ng_seg
            mt_umi[cs:ce] = mt_seg
        pct_mt = 100.0 * mt_umi / np.maximum(umi, 1.0)
    return umi, ngenes, pct_mt


def make_10x_csc_subset_readers(
    h5_path: Path,
    keep_cols: np.ndarray,
    n_raw_features: int,
    rna_row_idx: Optional[np.ndarray],
) -> Tuple[Callable[[int, int], np.ndarray], Callable[[np.ndarray], np.ndarray]]:
    """
    Create memory-bounded `(read_chunk_fn, read_rows_fn)` over a sorted subset `keep_cols`
    of columns in a 10x CSC `.h5` file `(n_raw_features, n_raw_cells)`.
    """
    n_out_genes = len(rna_row_idx) if rna_row_idx is not None else n_raw_features
    if rna_row_idx is not None and len(rna_row_idx) == n_raw_features and np.array_equal(rna_row_idx, np.arange(n_raw_features)):
        rna_row_idx = None

    def _read_cols(raw_cols: np.ndarray) -> np.ndarray:
        if len(raw_cols) == 0:
            return np.zeros((0, n_out_genes), dtype=np.float64)
        if len(raw_cols) > 1 and not np.all(raw_cols[:-1] <= raw_cols[1:]):
            raise ValueError("raw_cols must be sorted ascending")
        c_min, c_max = int(raw_cols[0]), int(raw_cols[-1])
        span = c_max - c_min + 1
        with h5py.File(h5_path, "r") as f:
            m = f["matrix"]
            indptr = m["indptr"]
            data = m["data"]
            indices = m["indices"]
            if len(raw_cols) / max(span, 1) >= 0.15:
                ptrs = indptr[c_min : c_max + 2]
                p0, p1 = int(ptrs[0]), int(ptrs[-1])
                csc_span = sp.csc_matrix(
                    (data[p0:p1], indices[p0:p1], ptrs - p0),
                    shape=(n_raw_features, span),
                )
                sub = csc_span[:, raw_cols - c_min]
            else:
                # Sparse cell-subset access (e.g., Pass 1 NTC-only reads): pull only target columns
                col_starts = indptr[raw_cols]
                col_ends = indptr[raw_cols + 1]
                all_d, all_i = [], []
                new_ptr = np.zeros(len(raw_cols) + 1, dtype=np.int64)
                cur_p = 0
                for j, (s0, s1) in enumerate(zip(col_starts, col_ends)):
                    if s1 > s0:
                        all_d.append(data[int(s0) : int(s1)])
                        all_i.append(indices[int(s0) : int(s1)])
                        cur_p += int(s1 - s0)
                    new_ptr[j + 1] = cur_p
                d_cat = np.concatenate(all_d) if all_d else np.empty(0, dtype=np.float32)
                i_cat = np.concatenate(all_i) if all_i else np.empty(0, dtype=np.int32)
                sub = sp.csc_matrix((d_cat, i_cat, new_ptr), shape=(n_raw_features, len(raw_cols)))
            if rna_row_idx is not None:
                sub = sub[rna_row_idx, :]
            return sub.T.toarray().astype(np.float64, copy=False)

    return (lambda s, e: _read_cols(keep_cols[s:e])), (lambda loc_idx: _read_cols(keep_cols[loc_idx]))


def extract_top2_csc(csc_mat: sp.csc_matrix) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return `(top1_umi, top2_umi, argmax1)` across columns of CSC matrix."""
    nc = csc_mat.shape[1]
    t1 = np.zeros(nc, dtype=np.float64)
    t2 = np.zeros(nc, dtype=np.float64)
    a1 = np.zeros(nc, dtype=np.int64)
    indptr = csc_mat.indptr
    data = csc_mat.data
    indices = csc_mat.indices
    for j in range(nc):
        s, e = int(indptr[j]), int(indptr[j + 1])
        n = e - s
        if n == 0:
            continue
        if n == 1:
            t1[j] = data[s]
            a1[j] = indices[s]
        else:
            d = data[s:e]
            o = np.argpartition(d, -2)[-2:]
            i1, i2 = (o[1], o[0]) if d[o[1]] >= d[o[0]] else (o[0], o[1])
            t1[j] = d[i1]
            t2[j] = d[i2]
            a1[j] = indices[s + i1]
    return t1, t2, a1


# ==============================================================================
# Dataset Adapters (WTC11 i3Neuron, Calu-3, DLD-1 Low-MOI, CD4+ T Rest)
# ==============================================================================


def build_wtc11_i3neuron_parts(raw_dir: Path, hgnc: HGNCLookup) -> Tuple[List[PartSpec], dict]:
    """Tian & Kampmann 2021 WTC11 iPSC-derived i3Neuron CRISPRi (`TianKampmann2021_CRISPRi.h5ad`)."""
    p = raw_dir / "wtc11_i3neuron" / "TianKampmann2021_CRISPRi.h5ad"
    if not p.exists():
        raise FileNotFoundError(p)
    f = h5py.File(p, "r")
    x_node, shape = open_x_node(f, x_layer="X")
    obs = f["obs"]
    batches_raw = np.array(read_obs_or_var_col(obs, "batch"), dtype=str)
    perts_raw = read_obs_or_var_col(obs, "perturbation")
    guides_raw = read_obs_or_var_col(obs, "guide_id")
    nperts = obs["nperts"][:]
    pm = obs["percent_mito"][:]
    ng = obs["ngenes"][:]
    umi_all = obs["ncounts"][:].astype(np.float64)

    is_nt_all = np.array([x == "control" for x in perts_raw], dtype=bool)
    qc_ok = (pm < 20.0) & (ng >= 500) & (is_nt_all | (nperts == 1))
    keep_idx = np.where(qc_ok)[0]
    nt_mask = is_nt_all[keep_idx]
    batches = np.array([f"lib_{b}" for b in batches_raw[keep_idx]], dtype=str)
    half_flags = np.where(np.isin(batches_raw[keep_idx], ["1", "3"]), 1, 2).astype(np.int8)

    pert_keys: List[str] = []
    pert_meta: Dict[str, Tuple[str, str, str, str]] = {}
    for gi in keep_idx:
        if is_nt_all[gi]:
            pert_keys.append("")
            continue
        raw_sym = str(perts_raw[gi]).strip()
        pk = str(guides_raw[gi]).strip()
        canon_sym, status = hgnc.resolve(raw_sym, "")
        eid = hgnc.sym_to_ens.get(canon_sym, "")
        pert_keys.append(pk)
        if pk not in pert_meta:
            pert_meta[pk] = (eid, canon_sym, raw_sym, status)

    # Remap in-memory CSR directly to kept rows so chunks are strictly contiguous
    x_kept = x_node[keep_idx]
    var_ids = read_obs_or_var_col(f["var"], "ensembl_id")
    var_names = read_obs_or_var_col(f["var"], "_index")
    f.close()
    n_kept, n_genes_raw = x_kept.shape

    part = PartSpec(
        tag="TianKampmann2021",
        path=p,
        n_cells=n_kept,
        n_genes_raw=n_genes_raw,
        nt_mask=nt_mask,
        batches=batches,
        half_flags=half_flags,
        umi=umi_all[keep_idx],
        pert_keys=pert_keys,
        pert_meta=pert_meta,
        var_gene_ids=var_ids,
        var_gene_names=var_names,
        read_chunk_fn=lambda s, e, mat=x_kept: read_h5ad_chunk_dense(mat, s, e, n_genes_raw),
    )
    audit = {
        "n_raw_cells": int(shape[0]),
        "n_qc_kept_cells": int(n_kept),
        "half_a_batches_spec": ["lib_1", "lib_3"],
        "half_b_batches_spec": ["lib_2", "lib_4"],
    }
    return [part], audit


def build_calu3_parts(raw_dir: Path, hgnc: HGNCLookup) -> Tuple[List[PartSpec], dict]:
    """Sunshine & Hein 2023 Calu-3 CRISPRi (`SunshineHein2023.h5ad`), single-gene elements ONLY."""
    p = raw_dir / "calu3" / "SunshineHein2023.h5ad"
    if not p.exists():
        raise FileNotFoundError(p)
    f = h5py.File(p, "r")
    x_node, shape = open_x_node(f, x_layer="X")
    obs = f["obs"]
    match_type = np.array(read_obs_or_var_col(obs, "match_type"), dtype=str)
    mle = read_obs_or_var_col(obs, "matched_library_element")
    pert = read_obs_or_var_col(obs, "perturbation")
    gem = obs["gem_group"][:]
    pm = obs["percent_mito"][:]
    ng = obs["ngenes"][:]
    umi_all = obs["ncounts"][:].astype(np.float64)

    var_ids = read_obs_or_var_col(f["var"], "ensembl_id")
    var_names = read_obs_or_var_col(f["var"], "gene_symbol")

    valid_match = np.isin(match_type, ["exact_match", "likely_match"]) & (pm < 20.0) & (ng >= 500)
    is_ctrl = np.array([x == "control" for x in pert], dtype=bool)
    is_single_gene = np.array([("_" not in str(x)) and (str(x) not in ("-1", "control", "")) for x in pert], dtype=bool)
    n_dual_dropped = int(np.sum(valid_match & ~is_ctrl & ~is_single_gene))
    keep_ok = valid_match & (is_ctrl | is_single_gene)
    keep_idx = np.where(keep_ok)[0]

    nt_mask = is_ctrl[keep_idx]
    gem_k = gem[keep_idx]
    batches = np.array([f"gem_{int(g)}" for g in gem_k], dtype=str)
    half_flags = np.where(np.isin(gem_k, [1, 3, 5, 7]), 1, 2).astype(np.int8)

    pert_keys: List[str] = []
    pert_meta: Dict[str, Tuple[str, str, str, str]] = {}
    for gi in keep_idx:
        if is_ctrl[gi]:
            pert_keys.append("")
            continue
        raw_sym = str(pert[gi]).strip()
        pk = str(mle[gi]).strip()
        canon_sym, status = hgnc.resolve(raw_sym, "")
        eid = hgnc.sym_to_ens.get(canon_sym, "")
        pert_keys.append(pk)
        if pk not in pert_meta:
            pert_meta[pk] = (eid, canon_sym, raw_sym, status)

    x_kept = x_node[keep_idx]
    f.close()
    n_kept, n_genes_raw = x_kept.shape
    part = PartSpec(
        tag="SunshineHein2023",
        path=p,
        n_cells=n_kept,
        n_genes_raw=n_genes_raw,
        nt_mask=nt_mask,
        batches=batches,
        half_flags=half_flags,
        umi=umi_all[keep_idx],
        pert_keys=pert_keys,
        pert_meta=pert_meta,
        var_gene_ids=var_ids,
        var_gene_names=var_names,
        read_chunk_fn=lambda s, e, mat=x_kept: read_h5ad_chunk_dense(mat, s, e, n_genes_raw),
    )
    audit = {
        "n_raw_cells": int(shape[0]),
        "n_dual_gene_dropped": int(n_dual_dropped),
        "n_qc_kept_single_gene_and_ntc_cells": int(n_kept),
        "half_a_batches_spec": ["gem_1", "gem_3", "gem_5", "gem_7"],
        "half_b_batches_spec": ["gem_2", "gem_4", "gem_6", "gem_8"],
    }
    return [part], audit


def build_dld1_parts(raw_dir: Path, hgnc: HGNCLookup) -> Tuple[List[PartSpec], dict]:
    """DLD-1 Low-MOI (`GSE337988`, 24 channels `GSM9864227..GSM9864250`), excluding `A1` & `NTC_MENG*`."""
    ddir = raw_dir / "dld1"
    rna_files = sorted(ddir.glob("GSM98642*.filtered_feature_bc_matrix.h5"))
    if len(rna_files) != 24:
        raise RuntimeError(f"Expected 24 DLD-1 Low-MOI RNA .h5 files, found {len(rna_files)}")

    parts: List[PartSpec] = []
    tot_raw_rna = tot_sing_h = tot_s1 = tot_ntc = tot_meng_excl = 0

    for ch_idx, r_path in enumerate(rna_files, start=1):
        gsm = r_path.name.split("_")[0]
        h_path = sorted(ddir.glob(f"{gsm}_*.hashing.umi_correct.h5"))[0]
        c_path = sorted(ddir.glob(f"{gsm}_*.crispr.umi_correct.h5"))[0]

        with h5py.File(r_path, "r") as fr:
            rm = fr["matrix"]
            n_genes_raw, n_rna_cells = int(rm["shape"][0]), int(rm["shape"][1])
            rna_bcs = [x.decode().split("-")[0] for x in rm["barcodes"][:]]
            var_ids = [x.decode() for x in rm["features/id"][:]]
            var_names = [x.decode() for x in rm["features/name"][:]]
        mt_raw = np.array([g.startswith("MT-") for g in var_names], dtype=bool)
        rna_umi, rna_ng, rna_mt = compute_10x_csc_rna_qc(r_path, None, mt_raw)
        qc_rna = (rna_mt < 20.0) & (rna_ng >= 500)

        with h5py.File(h_path, "r") as fh:
            hm = fh["matrix"]
            hbcs = [x.decode().split("-")[0] for x in hm["barcodes"][:]]
            cmo_names = [x.decode() for x in hm["features/name"][:]]
            H_csc = sp.csc_matrix((hm["data"][:], hm["indices"][:], hm["indptr"][:]), shape=tuple(hm["shape"][:]))
        with h5py.File(c_path, "r") as fc:
            cm = fc["matrix"]
            cbcs = [x.decode().split("-")[0] for x in cm["barcodes"][:]]
            gnames = [x.decode() for x in cm["features/name"][:]]
            C_csc = sp.csc_matrix((cm["data"][:], cm["indices"][:], cm["indptr"][:]), shape=tuple(cm["shape"][:]))

        h_map = {b: j for j, b in enumerate(hbcs)}
        c_map = {b: j for j, b in enumerate(cbcs)}

        matched_r_idx = [i for i, b in enumerate(rna_bcs) if qc_rna[i] and b in h_map and b in c_map]
        matched_h_cols = [h_map[rna_bcs[i]] for i in matched_r_idx]
        matched_c_cols = [c_map[rna_bcs[i]] for i in matched_r_idx]

        H = H_csc[:, matched_h_cols].T.toarray()
        Hs = np.sort(H, axis=1)
        h_arg = H.argmax(axis=1)
        is_cmo_ok = np.isin([cmo_names[a] for a in h_arg], list(DLD1_LOW_MOI_CMOS))
        sing_h = (Hs[:, -1] >= 10) & (Hs[:, -1] / np.maximum(Hs[:, -2], 0.5) >= 3.0) & is_cmo_ok

        C_sub = C_csc[:, np.asarray(matched_c_cols)[sing_h]]
        c_t1, c_t2, c_a1 = extract_top2_csc(C_sub)
        sing_c = (c_t1 >= 5.0) & (c_t2 < 3.0) & (c_t1 / np.maximum(c_t2, 0.5) >= 3.0)

        h_pass_r = np.asarray(matched_r_idx, dtype=np.int64)[sing_h]
        keep_cols: List[int] = []
        nt_flags: List[bool] = []
        pkeys: List[str] = []
        pmeta: Dict[str, Tuple[str, str, str, str]] = {}

        for loc_j, r_i in enumerate(h_pass_r):
            if not sing_c[loc_j]:
                continue
            gname = gnames[int(c_a1[loc_j])]
            if gname.startswith("NTC_MENG"):
                tot_meng_excl += 1
                continue
            if gname.startswith("NTC_"):
                keep_cols.append(int(r_i))
                nt_flags.append(True)
                pkeys.append("")
                tot_ntc += 1
                continue
            m = DLD1_TARGET_RE.match(gname)
            if m is None:
                continue
            raw_sym = m.group(1)
            raw_ens = m.group(2) or ""
            canon_sym, hstat = hgnc.resolve(raw_sym, raw_ens)
            ens = raw_ens or hgnc.sym_to_ens.get(canon_sym, "")
            keep_cols.append(int(r_i))
            nt_flags.append(False)
            pkeys.append(gname)
            if gname not in pmeta:
                pmeta[gname] = (ens, canon_sym, raw_sym, hstat)

        keep_arr = np.asarray(keep_cols, dtype=np.int64)
        ord_k = np.argsort(keep_arr)
        keep_arr = keep_arr[ord_k]
        nt_arr = np.asarray(nt_flags, dtype=bool)[ord_k]
        pkeys = [pkeys[i] for i in ord_k]

        tot_raw_rna += n_rna_cells
        tot_sing_h += int(sing_h.sum())
        tot_s1 += len(keep_arr)

        read_chunk, read_rows = make_10x_csc_subset_readers(r_path, keep_arr, n_genes_raw, None)
        half_flag = np.full(len(keep_arr), 1 if (ch_idx % 2 == 1) else 2, dtype=np.int8)
        parts.append(
            PartSpec(
                tag=gsm,
                path=r_path,
                n_cells=len(keep_arr),
                n_genes_raw=n_genes_raw,
                nt_mask=nt_arr,
                batches=np.full(len(keep_arr), f"ch_{ch_idx:02d}", dtype=str),
                half_flags=half_flag,
                umi=rna_umi[keep_arr],
                pert_keys=pkeys,
                pert_meta=pmeta,
                var_gene_ids=var_ids,
                var_gene_names=var_names,
                read_chunk_fn=read_chunk,
                read_rows_fn=read_rows,
            )
        )

    audit = {
        "n_channels": 24,
        "n_raw_rna_cells": int(tot_raw_rna),
        "n_low_moi_hashing_singlets_excl_A1": int(tot_sing_h),
        "n_single_sgRNA_and_ntc_kept": int(tot_s1),
        "n_ntc_kept": int(tot_ntc),
        "n_ntc_meng_excluded": int(tot_meng_excl),
        "half_rule": "odd_channels_1..23_vs_even_channels_2..24",
    }
    return parts, audit


def fit_poisson_gaussian_thresholds_from_hist(
    hist: np.ndarray,
    min_thr: int = 5,
    max_thr: int = 20,
    default_thr: int = 7,
    min_entries: int = 50,
) -> np.ndarray:
    """
    Fit 2-component bimodal EM per guide on log10 binned counts (`hist` shape `(n_guides, 66)`,
    index 1..65) to separate background ambient UMIs from signal UMIs, floored at `min_thr`.
    """
    n_guides = hist.shape[0]
    thrs = np.full(n_guides, default_thr, dtype=np.int32)
    vals = np.arange(1, hist.shape[1], dtype=np.float64)
    lx = np.log10(vals)
    u_grid = np.arange(min_thr, max_thr + 1, dtype=np.float64)
    lu = np.log10(u_grid)

    for g in range(n_guides):
        w = hist[g, 1:].astype(np.float64)
        tot = float(w.sum())
        if tot < min_entries:
            continue
        # Weighted quantiles for robust initialization
        cw = np.cumsum(w) / tot
        m0 = float(lx[min(len(lx) - 1, int(np.searchsorted(cw, 0.25)))])
        m1 = max(float(lx[min(len(lx) - 1, int(np.searchsorted(cw, 0.85)))]), np.log10(float(min_thr)))
        if m1 <= m0 + 0.15:
            continue
        pi1 = 0.25
        s0 = s1 = 0.25
        for _ in range(15):
            p0 = (1.0 - pi1) * np.exp(-0.5 * ((lx - m0) / s0) ** 2) / s0
            p1 = pi1 * np.exp(-0.5 * ((lx - m1) / s1) ** 2) / s1
            denom = np.maximum(p0 + p1, 1e-300)
            r1 = p1 / denom
            w1 = float(np.sum(w * r1))
            w0 = tot - w1
            if w1 < 2.0 or w0 < 2.0:
                break
            pi1 = float(np.clip(w1 / tot, 0.01, 0.99))
            m0 = float(np.sum(w * (1.0 - r1) * lx) / w0)
            m1 = float(np.sum(w * r1 * lx) / w1)
            s0 = float(np.clip(np.sqrt(np.sum(w * (1.0 - r1) * (lx - m0) ** 2) / w0), 0.10, 1.0))
            s1 = float(np.clip(np.sqrt(np.sum(w * r1 * (lx - m1) ** 2) / w1), 0.10, 1.0))
        if m1 > m0 + 0.25:
            p0g = (1.0 - pi1) * np.exp(-0.5 * ((lu - m0) / s0) ** 2) / s0
            p1g = pi1 * np.exp(-0.5 * ((lu - m1) / s1) ** 2) / s1
            pos = np.where(p1g >= p0g)[0]
            if len(pos) > 0:
                thrs[g] = int(u_grid[pos[0]])
    return thrs


def build_cd4t_rest_parts(raw_dir: Path, hgnc: HGNCLookup) -> Tuple[List[PartSpec], dict]:
    """
    CD4+ T Resting CRISPRi (`GSE314342`), two donors:
      - Half A (`1`): Donor `CE0006864` (`D4_Rest`, `CD4i_R2L01..24`, 24 lanes)
      - Half B (`2`): Donor `CE0008162` (`D1_Rest`, `CD4i_R1L01..23`, 23 lanes)
    Uses pooled per-donor Poisson-Gaussian guide threshold (`>= 5` UMI floor + single-sgRNA
    ratio >= 3) and post-guide paper QC (`<= 5% MT & >= 500 genes`).
    """
    cdir = raw_dir / "cd4t_rest"
    sg_df = pd.read_csv(cdir / "sgrna_library_metadata.suppl_table.csv")
    guide_to_target: Dict[str, Tuple[str, str, str, str]] = {}
    ntc_guide_set: set = set()
    for r in sg_df.itertuples(index=False):
        g = str(r.sgRNA).strip()
        if g.startswith("NTC-"):
            ntc_guide_set.add(g)
        elif pd.notna(r.target_gene_name) and str(r.target_gene_name).strip():
            raw_sym = str(r.target_gene_name).strip()
            raw_ens = str(r.target_gene_id).strip() if pd.notna(r.target_gene_id) else ""
            canon_sym, hstat = hgnc.resolve(raw_sym, raw_ens)
            ens = raw_ens or hgnc.sym_to_ens.get(canon_sym, "")
            guide_to_target[g] = (ens, canon_sym, raw_sym, hstat)

    d4_files = sorted(cdir.glob("GSM93938*_D4_Rest_sample_filtered_feature_bc_matrix.h5"))
    d1_files = sorted(cdir.glob("GSM93938*_D1_Rest_sample_filtered_feature_bc_matrix.h5"))
    if len(d4_files) != 24 or len(d1_files) != 23:
        raise RuntimeError(f"Expected 24 D4 + 23 D1 files, found D4={len(d4_files)}, D1={len(d1_files)}")

    parts: List[PartSpec] = []
    donor_audit: Dict[str, dict] = {}

    for donor_id, half_val, files in [
        ("CE0006864_D4", 1, d4_files),
        ("CE0008162_D1", 2, d1_files),
    ]:
        print(f"  [cd4t_rest] Scanning {donor_id} ({len(files)} lanes) for QC & pooled guide histogram...", flush=True)
        # First pass over donor lanes: extract RNA QC and CRISPR matrix slice, accumulate histogram
        lane_cache = []
        hist = None
        crp_names_ref = None
        for lp in files:
            parts_name = lp.name.split("_")
            lane_tag = f"{parts_name[1]}_{parts_name[2]}"  # CD4i_R2L01 / CD4i_R1L01
            with h5py.File(lp, "r") as f:
                m = f["matrix"]
                n_feat, n_col = int(m["shape"][0]), int(m["shape"][1])
                ftypes = np.array([x.decode() for x in m["features/feature_type"][:]])
                fids = np.array([x.decode() for x in m["features/id"][:]])
                fnames = np.array([x.decode() for x in m["features/name"][:]])
                rna_rows = np.where(ftypes == "Gene Expression")[0]
                crp_rows = np.where(ftypes != "Gene Expression")[0]
                if crp_names_ref is None:
                    crp_names_ref = [CD4_EXCEL_GUIDE_REMAP.get(x, x) for x in fnames[crp_rows]]
                    hist = np.zeros((len(crp_rows), 66), dtype=np.int64)
                mt_raw = np.zeros(n_feat, dtype=bool)
                mt_raw[rna_rows] = [x.startswith("MT-") for x in fnames[rna_rows]]

            umi_l, ng_l, mt_l = compute_10x_csc_rna_qc(lp, rna_rows, mt_raw)
            ok_init = (mt_l < 20.0) & (ng_l >= 200)
            ok_hq = (mt_l <= 5.0) & (ng_l >= 500)

            # Extract CRISPR submatrix only for ok_hq cells & initial cells for histogram
            is_crp = np.zeros(n_feat, dtype=bool)
            is_crp[crp_rows] = True
            crp_remap_idx = np.full(n_feat, -1, dtype=np.int32)
            crp_remap_idx[crp_rows] = np.arange(len(crp_rows), dtype=np.int32)

            with h5py.File(lp, "r") as f:
                m = f["matrix"]
                indptr = m["indptr"][:]
                data = m["data"][:]
                indices = m["indices"][:]
            m_all = is_crp[indices]
            c_idx = crp_remap_idx[indices[m_all]]
            c_dat = data[m_all].astype(np.int64, copy=False)
            col_lens = np.diff(indptr)
            c_col = np.repeat(np.arange(n_col, dtype=np.int32), col_lens)[m_all]

            # Build CSC matrix ONLY for ok_init cells
            ok_init_mask = ok_init[c_col]
            bins = np.clip(c_dat[ok_init_mask], 1, 65)
            np.add.at(hist, (c_idx[ok_init_mask], bins), 1)

            hq_cols = np.where(ok_hq)[0]
            col_map_hq = np.full(n_col, -1, dtype=np.int32)
            col_map_hq[hq_cols] = np.arange(len(hq_cols), dtype=np.int32)
            hq_m = col_map_hq[c_col] >= 0
            # Construct compact CSC across ok_hq cells
            C_hq = sp.csr_matrix(
                (c_dat[hq_m].astype(np.float64), (c_idx[hq_m], col_map_hq[c_col[hq_m]])),
                shape=(len(crp_rows), len(hq_cols)),
            ).tocsc()
            del indptr, data, indices, c_idx, c_dat, c_col, m_all, col_lens

            lane_cache.append(
                (
                    lp,
                    lane_tag,
                    n_feat,
                    rna_rows,
                    [fids[i] for i in rna_rows],
                    [fnames[i] for i in rna_rows],
                    umi_l,
                    hq_cols,
                    C_hq,
                )
            )

        thrs = fit_poisson_gaussian_thresholds_from_hist(hist, min_thr=5, max_thr=20, default_thr=7)
        d_ntc = d_targ = 0
        for lp, lane_tag, n_feat, rna_rows, rids, rnames, umi_l, hq_cols, C_hq in lane_cache:
            t1, t2, a1 = extract_top2_csc(C_hq)
            pos_mask = C_hq.data >= thrs[C_hq.indices]
            pos_csc = sp.csc_matrix((pos_mask.astype(np.int8), C_hq.indices, C_hq.indptr), shape=C_hq.shape)
            n_pos = np.asarray(pos_csc.sum(axis=0)).ravel()
            s1_mask = (n_pos == 1) & (t1 >= 5.0) & ((t2 < 3.0) | (t1 / np.maximum(t2, 0.5) >= 3.0))

            keep_cols: List[int] = []
            nt_flags: List[bool] = []
            pkeys: List[str] = []
            pmeta: Dict[str, Tuple[str, str, str, str]] = {}

            for loc_j in np.where(s1_mask)[0]:
                gname = crp_names_ref[int(a1[loc_j])]
                raw_col = int(hq_cols[loc_j])
                if gname in ntc_guide_set:
                    keep_cols.append(raw_col)
                    nt_flags.append(True)
                    pkeys.append("")
                    d_ntc += 1
                elif gname in guide_to_target:
                    keep_cols.append(raw_col)
                    nt_flags.append(False)
                    pkeys.append(gname)
                    if gname not in pmeta:
                        pmeta[gname] = guide_to_target[gname]
                    d_targ += 1

            keep_arr = np.asarray(keep_cols, dtype=np.int64)
            nt_arr = np.asarray(nt_flags, dtype=bool)
            read_chunk, read_rows = make_10x_csc_subset_readers(lp, keep_arr, n_feat, rna_rows)
            parts.append(
                PartSpec(
                    tag=lane_tag,
                    path=lp,
                    n_cells=len(keep_arr),
                    n_genes_raw=len(rna_rows),
                    nt_mask=nt_arr,
                    batches=np.full(len(keep_arr), f"{donor_id}_{lane_tag}", dtype=str),
                    half_flags=np.full(len(keep_arr), half_val, dtype=np.int8),
                    umi=umi_l[keep_arr],
                    pert_keys=pkeys,
                    pert_meta=pmeta,
                    var_gene_ids=rids,
                    var_gene_names=rnames,
                    read_chunk_fn=read_chunk,
                    read_rows_fn=read_rows,
                )
            )
        donor_audit[donor_id] = {
            "n_lanes": len(files),
            "n_single_ntc": int(d_ntc),
            "n_single_targeting": int(d_targ),
            "median_guide_threshold": float(np.median(thrs)),
        }
        print(f"    -> {donor_id} single-NTC={d_ntc:,}, single-targeting={d_targ:,}, median_thr={np.median(thrs):.1f}", flush=True)
        del lane_cache
        gc.collect()

    return parts, {"donors": donor_audit, "half_rule": "HalfA=CE0006864_D4_Rest, HalfB=CE0008162_D1_Rest"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Prepare R10 standardized shared cell lines (A-T1).")
    p.add_argument(
        "--line",
        required=True,
        choices=["wtc11_i3neuron", "calu3", "dld1", "cd4t_rest", "all"],
        help="Which R10 line(s) to process.",
    )
    p.add_argument("--workers", type=int, default=8, help="Parallel worker processes for AD/Energy tests.")
    p.add_argument("--raw-dir", type=Path, default=R10_RAW_DIR)
    p.add_argument("--out-root", type=Path, default=R10_LINES_DIR)
    return p.parse_args()


def run_line(line: str, raw_dir: Path, out_root: Path, hgnc: HGNCLookup, workers: int) -> Path:
    builders = {
        "wtc11_i3neuron": build_wtc11_i3neuron_parts,
        "calu3": build_calu3_parts,
        "dld1": build_dld1_parts,
        "cd4t_rest": build_cd4t_rest_parts,
    }
    parts, audit = builders[line](raw_dir, hgnc)
    return process_r10_line(
        line_name=line,
        parts=parts,
        hgnc=hgnc,
        out_root=out_root,
        extra_meta={"adapter_audit": audit},
        n_workers=workers,
    )


def main() -> None:
    args = parse_args()
    hgnc = build_hgnc_lookup()
    targets = ["wtc11_i3neuron", "calu3", "dld1", "cd4t_rest"] if args.line == "all" else [args.line]
    for ln in targets:
        run_line(ln, args.raw_dir, args.out_root, hgnc, args.workers)


if __name__ == "__main__":
    main()
