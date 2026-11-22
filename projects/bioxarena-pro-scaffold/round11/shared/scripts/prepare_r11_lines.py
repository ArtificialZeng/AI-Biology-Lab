#!/usr/bin/env python3
"""
Unified Round 11 Single-Cell Preprocessing & Addendum-B QC Engine (A-T1).

Implements PREREG_D §2.2-2.3 & Addendum A strictly reusing Round 10's frozen
standardization, HGNC lookup, split-half centering, and Addendum-B reproducibility
primitives (`round10/shared/scripts/prepare_r10_lines.py`).

BLINDING INVARIANT:
  Computes ONLY counts, QC, split-half deltas, and Addendum-B reproducibility flags.
  Strictly prohibits computing any evaluation predictions, correlations, slopes, or FCT.
"""

from __future__ import annotations

import argparse
import gc
import gzip
import json
import os
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

SCRIPTS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = Path(__file__).resolve().parents[3]
R10_SCRIPTS = PROJECT_ROOT / "round10" / "shared" / "scripts"
if str(R10_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(R10_SCRIPTS))

from prepare_r10_lines import (  # noqa: E402
    CD4_EXCEL_GUIDE_REMAP,
    DEFAULT_HGNC_PATH,
    HGNCLookup,
    PartSpec,
    build_hgnc_lookup,
    compute_10x_csc_rna_qc,
    extract_top2_csc,
    fit_poisson_gaussian_thresholds_from_hist,
    make_10x_csc_subset_readers,
    process_r10_line,
)

R11_SHARED_DIR = PROJECT_ROOT / "round11" / "shared"
R11_RAW_DIR = R11_SHARED_DIR / "raw"
R11_LINES_DIR = R11_SHARED_DIR / "lines"
GO_PATH = PROJECT_ROOT / "round11" / "status" / "R11_GO.txt"

TELOHAEC_CHAN_ORDER: Tuple[str, ...] = tuple(
    [f"scRNAseq_2kG_11AMDox_{i}" for i in range(1, 7)]
    + [f"scRNAseq_2kG_2PMDox1_{i}" for i in range(1, 8)]
    + [f"scRNAseq_2kG_2PMDox2_{i}" for i in range(1, 8)]
)
TELOHAEC_CHAN_TO_LANE: Dict[str, int] = {ch: idx for idx, ch in enumerate(TELOHAEC_CHAN_ORDER, start=1)}

HCASMC_SUBLIBS_SORTED: Tuple[str, ...] = (
    "3A", "3B", "3C", "3D", "3E", "3F", "3G", "3H",
    "4A", "4B", "4C", "4D", "4E", "4F",
)
FENG_KOLF_EXCLUDED = {"kolf_2", "kolf_3"}


def require_r11_go() -> None:
    if not GO_PATH.exists():
        raise RuntimeError(f"Missing {GO_PATH}! Refusing to touch R11 data before R11_GO.txt.")


def is_eligible_telohaec_plain_target(raw_sym: str, hgnc: HGNCLookup) -> Tuple[bool, str, str, str]:
    """
    Filter TeloHAEC targeting symbols strictly per PREREG_D §2.2.1:
      - Reject empty, negative-control, safe-targeting
      - Reject alternative TSS (-TSS*) & enhancers/multi-gene loci (Enhancer*, 9p21*, *-and-*)
      - Require clean resolution via frozen HGNC lookup (or approved symbol)
    """
    s = str(raw_sym).strip()
    if not s or s.lower() in ("nan", "none", "na", "negative-control", "safe-targeting"):
        return False, "", "", "control_or_empty"
    if any(tok in s for tok in ("-TSS", "Enhancer", "9p21", "-and-", "_and_")):
        return False, "", "", "excluded_non_plain_tss_or_cis"
    canon_sym, hstat = hgnc.resolve(s, "")
    if hstat == "unmapped" and canon_sym not in hgnc.approved:
        return False, "", "", "unmapped"
    ens = hgnc.sym_to_ens.get(canon_sym, "")
    return True, ens, canon_sym, hstat


def make_csr_h5_subset_readers(
    h5_path: Path,
    x_group_path: str,
    keep_rows: np.ndarray,
    n_genes_raw: int,
) -> Tuple[
    Callable[[int, int], np.ndarray],
    Callable[[np.ndarray], np.ndarray],
    Callable[[int, int, np.ndarray], np.ndarray],
]:
    """
    Memory-bounded row-subset readers for CSR HDF5 (`cells x genes`) datasets such as
    `.h5ad` `/X` or temporary streamed CSR `.h5`. `keep_rows` must be sorted ascending.
    """
    if len(keep_rows) > 1 and not np.all(keep_rows[:-1] <= keep_rows[1:]):
        raise ValueError("keep_rows must be sorted ascending")

    def _read_rows(raw_rows: np.ndarray, feat_cols: Optional[np.ndarray] = None) -> np.ndarray:
        out_dtype = np.float32 if feat_cols is not None else np.float64
        out_w = len(feat_cols) if feat_cols is not None else n_genes_raw
        if len(raw_rows) == 0:
            return np.zeros((0, out_w), dtype=out_dtype)
        r_min, r_max = int(raw_rows[0]), int(raw_rows[-1])
        span = r_max - r_min + 1
        with h5py.File(h5_path, "r") as f:
            xg = f[x_group_path]
            indptr = xg["indptr"]
            data = xg["data"]
            indices = xg["indices"]
            if len(raw_rows) / max(span, 1) >= 0.20:
                ptrs = indptr[r_min : r_max + 2].astype(np.int64)
                p0, p1 = int(ptrs[0]), int(ptrs[-1])
                csr_span = sp.csr_matrix(
                    (data[p0:p1], indices[p0:p1], ptrs - np.int64(p0)),
                    shape=(span, n_genes_raw),
                )
                sub = csr_span[raw_rows - r_min, :]
            else:
                r_starts = indptr[raw_rows]
                r_ends = indptr[raw_rows + 1]
                all_d, all_i = [], []
                new_ptr = np.zeros(len(raw_rows) + 1, dtype=np.int64)
                cur_p = 0
                for j, (s0, s1) in enumerate(zip(r_starts, r_ends)):
                    if s1 > s0:
                        all_d.append(data[int(s0) : int(s1)])
                        all_i.append(indices[int(s0) : int(s1)])
                        cur_p += int(s1 - s0)
                    new_ptr[j + 1] = cur_p
                d_cat = np.concatenate(all_d) if all_d else np.empty(0, dtype=np.float32)
                i_cat = np.concatenate(all_i) if all_i else np.empty(0, dtype=np.int32)
                sub = sp.csr_matrix((d_cat, i_cat, new_ptr), shape=(len(raw_rows), n_genes_raw))

            if feat_cols is not None:
                sub = sub[:, feat_cols]
            return sub.toarray().astype(out_dtype, copy=False)

    return (
        (lambda s, e: _read_rows(keep_rows[s:e], None)),
        (lambda loc_idx: _read_rows(keep_rows[loc_idx], None)),
        (lambda s, e, fcols: _read_rows(keep_rows[s:e], fcols)),
    )


def apply_stratified_seed42_split(pert_keys: Sequence[str], nt_mask: np.ndarray, seed: int = 42) -> np.ndarray:
    """
    Deterministic stratified 50/50 split (`np.random.default_rng(42)`) across QC-valid cells
    per target group and NTC pool, per PREREG_D §2.2.2 missing-unit rule.
    """
    n = len(nt_mask)
    half_flags = np.zeros(n, dtype=np.int8)
    rng = np.random.default_rng(seed)

    nt_idx = np.where(nt_mask)[0]
    perm_nt = nt_idx[rng.permutation(len(nt_idx))]
    half_flags[perm_nt[0::2]] = 1
    half_flags[perm_nt[1::2]] = 2

    by_pk: Dict[str, List[int]] = defaultdict(list)
    for i in range(n):
        if not nt_mask[i] and pert_keys[i]:
            by_pk[pert_keys[i]].append(i)
    for pk in sorted(by_pk.keys()):
        arr = np.asarray(by_pk[pk], dtype=np.int64)
        perm = arr[rng.permutation(len(arr))]
        half_flags[perm[0::2]] = 1
        half_flags[perm[1::2]] = 2
    return half_flags


def build_telohaec_parts(raw_dir: Path, hgnc: HGNCLookup) -> Tuple[List[PartSpec], dict]:
    """
    TeloHAEC_2285 (`IGVFDS7943IIWZ` / `GSE210681`):
      - `IGVFFI8645YMGX.h5ad` (214,449 QC'd cells x 17,472 genes `SYMBOL:ENSG`)
      - Cross-checked against `GSE210681_Dialout_guides_per_CBC.txt.gz` (`max_umi_ct>=4 & X2nd<4`)
      - Strictly `type == 'non-targeting'` (`400` guides) as NTC; excludes `safe-targeting` (`600`),
        `-TSS*`, `Enhancer*`, `9p21*`, `*-and-*`
      - Pools `P4` single-sgRNA cells per plain canonical target symbol (15 guides/TSS library)
      - Split-half: odd vs even sorted 10x channel (`1..20`: `1,3..19` vs `2,4..20`)
    """
    tdir = raw_dir / "telohaec"
    h5_p = tdir / "IGVFFI8645YMGX.h5ad"
    g_p = tdir / "IGVFFI9754AGFB.tsv.gz"
    d_p = tdir / "GSE210681_Dialout_guides_per_CBC.txt.gz"
    for p in (h5_p, g_p, d_p):
        if not p.exists():
            raise FileNotFoundError(p)

    with h5py.File(h5_p, "r") as f:
        bcs = [x.decode("utf-8", errors="replace") for x in f["obs/barcodes"][:]]
        vidx = [x.decode("utf-8", errors="replace") for x in f["var/_index"][:]]
        umi_all = f["obs/n_counts"][:].astype(np.float64)

    var_names = [v.split(":", 1)[0].strip() for v in vidx]
    var_ids = [v.split(":", 1)[1].strip() if ":" in v else "" for v in vidx]

    gdf = pd.read_csv(g_p, sep="\t", low_memory=False)
    sp_to_type = dict(zip(gdf["spacer"].astype(str), gdf["type"].astype(str)))

    dial = pd.read_csv(d_p, sep="\t", low_memory=False).set_index("X10x_cbc")
    x10_keys = []
    b_fields = []
    for b in bcs:
        raw_sym, spacer, cbc_chan = b.split(":", 2)
        cbc16, chan = cbc_chan.split("-", 1)
        lane = TELOHAEC_CHAN_TO_LANE[chan]
        b_fields.append((raw_sym, spacer, lane, chan))
        x10_keys.append(f"{cbc16}-{lane}")

    d_sub = dial.reindex(x10_keys)
    dial_sp = d_sub["guide_w_max"].astype(str).to_numpy()
    p4_ok = (
        d_sub["max_umi_ct"].notna().to_numpy()
        & (d_sub["max_umi_ct"].fillna(0).to_numpy() >= 4)
        & (d_sub["X2nd_place_umi_ct"].fillna(999).to_numpy() < 4)
        & (dial_sp == np.array([x[1] for x in b_fields], dtype=str))
    )

    keep_cols: List[int] = []
    nt_flags: List[bool] = []
    lanes_k: List[int] = []
    pkeys: List[str] = []
    pmeta: Dict[str, Tuple[str, str, str, str]] = {}

    n_fail_p4 = int((~p4_ok).sum())
    n_safe_excl = n_tss_excl = n_unmapped_excl = 0

    for idx_c, ((raw_sym, spacer, lane, chan), ok) in enumerate(zip(b_fields, p4_ok)):
        if not ok:
            continue
        gtype = sp_to_type.get(spacer, "")
        if gtype == "non-targeting":
            keep_cols.append(idx_c)
            nt_flags.append(True)
            lanes_k.append(lane)
            pkeys.append("")
            continue
        if gtype == "safe-targeting":
            n_safe_excl += 1
            continue
        if gtype != "targeting":
            continue
        elig, ens, canon_sym, hstat = is_eligible_telohaec_plain_target(raw_sym, hgnc)
        if not elig:
            if hstat == "excluded_non_plain_tss_or_cis":
                n_tss_excl += 1
            else:
                n_unmapped_excl += 1
            continue
        # Pool P4 single-sgRNA cells per target gene symbol (15 guides/TSS library)
        keep_cols.append(idx_c)
        nt_flags.append(False)
        lanes_k.append(lane)
        pkeys.append(canon_sym)
        if canon_sym not in pmeta:
            pmeta[canon_sym] = (ens, canon_sym, raw_sym, hstat)

    keep_arr = np.asarray(keep_cols, dtype=np.int64)
    nt_arr = np.asarray(nt_flags, dtype=bool)
    lanes_arr = np.asarray(lanes_k, dtype=np.int32)
    half_flags = np.where((lanes_arr % 2) == 1, 1, 2).astype(np.int8)

    # Verify >= 15 NTC per lane; merge any lane < 15 into HalfA/HalfB pooled
    lane_ntc_counts = Counter(int(l) for l, nt in zip(lanes_arr, nt_arr) if nt)
    low_ntc_lanes = {l for l in range(1, 21) if lane_ntc_counts.get(l, 0) < 15}
    batches = np.array(
        [
            (f"Half{'A' if l % 2 == 1 else 'B'}_pooled" if l in low_ntc_lanes else f"ch_{l:02d}")
            for l in lanes_arr
        ],
        dtype=str,
    )

    read_chunk, read_rows, read_chunk_cols = make_csr_h5_subset_readers(h5_p, "X", keep_arr, len(var_names))
    part = PartSpec(
        tag="TeloHAEC_2285",
        path=h5_p,
        n_cells=len(keep_arr),
        n_genes_raw=len(var_names),
        nt_mask=nt_arr,
        batches=batches,
        half_flags=half_flags,
        umi=umi_all[keep_arr],
        pert_keys=pkeys,
        pert_meta=pmeta,
        var_gene_ids=var_ids,
        var_gene_names=var_names,
        read_chunk_fn=read_chunk,
        read_rows_fn=read_rows,
        read_chunk_cols_fn=read_chunk_cols,
    )
    audit = {
        "n_h5ad_raw_cells": len(bcs),
        "n_failed_dialout_p4": n_fail_p4,
        "n_safe_targeting_excluded": n_safe_excl,
        "n_tss_or_cis_excluded": n_tss_excl,
        "n_unmapped_excluded": n_unmapped_excl,
        "n_kept_cells": len(keep_arr),
        "n_kept_ntc": int(nt_arr.sum()),
        "lane_ntc_counts": {str(k): int(v) for k, v in sorted(lane_ntc_counts.items())},
        "target_pooling_rule": "P4 single-sgRNA cells grouped by plain canonical HGNC symbol (15-guide/TSS library)",
        "half_rule": "odd_10x_lanes_1..19_vs_even_10x_lanes_2..20",
    }
    return [part], audit


def _parse_hcasmc_sublibrary(bc: str) -> str:
    if "__" in bc:
        return bc.split("__")[-1].strip()
    if "-" in bc:
        suf = bc.rsplit("-", 1)[-1].strip()
        if suf in HCASMC_SUBLIBS_SORTED or suf.isdigit():
            return suf
    return ""


def build_hcasmc_parts(raw_dir: Path, hgnc: HGNCLookup, mt_max_pct: float = 10.0) -> Tuple[List[PartSpec], dict]:
    """
    HCASMC-hTERT (`IGVFDS6248NDUW`):
      - `IGVFFI5999NXAF.h5` (`Cell Ranger HDF5` CSC from `parse-biosciences-v1.2.1`)
      - Checks `Q3` feasibility JSON first; uses `IGVFFI4575UMXX.tsv.gz` & `IGVFFI9546UECF.tsv.gz`
      - Strictly `type == 'non-targeting'` (`non-targeting_`) as NTC; excludes `0Safe_` (`safe-targeting`)
      - Split-half: sorted 14 Parse sublibraries (`3A..4F`) parity if present (`>= 2` suffixes),
        else deterministic `seed=42` stratified 50/50 split (`PREREG_D §2.2.2`).
    """
    hdir = raw_dir / "hcasmc"
    h5_p = hdir / "IGVFFI5999NXAF.h5"
    g_p = hdir / "IGVFFI4575UMXX.tsv.gz"
    g_coord_p = hdir / "IGVFFI9546UECF.tsv.gz"
    for p in (h5_p, g_p, g_coord_p):
        if not p.exists():
            raise FileNotFoundError(p)

    gdf = pd.read_csv(g_p, sep="\t", low_memory=False)
    gc_df = pd.read_csv(g_coord_p, sep="\t", low_memory=False)
    ens_map = dict(zip(gc_df["guide_id"].astype(str), gc_df["intended_target_name"].astype(str)))

    ref_ntc = set(gdf.loc[gdf["type"] == "non-targeting", "guide_id"].astype(str))
    ref_safe = set(gdf.loc[gdf["type"] == "safe-targeting", "guide_id"].astype(str))
    guide_to_target: Dict[str, Tuple[str, str, str, str]] = {}
    for gid in gdf.loc[gdf["type"] == "targeting", "guide_id"].astype(str):
        raw_sym = gid.rsplit("_", 1)[0].strip()
        raw_ens = str(ens_map.get(gid, "")).strip()
        if not raw_ens.startswith("ENSG"):
            raw_ens = ""
        canon_sym, hstat = hgnc.resolve(raw_sym, raw_ens)
        ens = raw_ens or hgnc.sym_to_ens.get(canon_sym, "")
        guide_to_target[gid] = (ens, canon_sym, raw_sym, hstat)

    with h5py.File(h5_p, "r") as f:
        m = f["matrix"]
        n_feat, n_col = int(m["shape"][0]), int(m["shape"][1])
        ftypes = np.array([x.decode("utf-8", errors="replace") for x in m["features/feature_type"][:]])
        fids = [x.decode("utf-8", errors="replace") for x in m["features/id"][:]]
        fnames = [x.decode("utf-8", errors="replace") for x in m["features/name"][:]]
        bcs = [x.decode("utf-8", errors="replace") for x in m["barcodes"][:]]
        indptr_all = m["indptr"][:]

        rna_rows = np.where(ftypes == "Gene Expression")[0]
        crp_rows = np.where(ftypes != "Gene Expression")[0]
        if len(crp_rows) == 0:
            raise RuntimeError("HCASMC missing non-GEX guide features inside .h5")

        id_hits = sum((fids[i] in ref_ntc or fids[i] in ref_safe or fids[i] in guide_to_target) for i in crp_rows)
        nm_hits = sum((fnames[i] in ref_ntc or fnames[i] in ref_safe or fnames[i] in guide_to_target) for i in crp_rows)
        crp_keys = [fnames[i] for i in crp_rows] if nm_hits >= id_hits else [fids[i] for i in crp_rows]

        mt_raw = np.zeros(n_feat, dtype=bool)
        mt_raw[rna_rows] = [fnames[i].startswith("MT-") for i in rna_rows]
        is_crp = np.zeros(n_feat, dtype=bool)
        is_crp[crp_rows] = True
        crp_remap_idx = np.full(n_feat, -1, dtype=np.int32)
        crp_remap_idx[crp_rows] = np.arange(len(crp_rows), dtype=np.int32)
        is_rna_feat = ~is_crp

        umi_l = np.zeros(n_col, dtype=np.float64)
        ng_l = np.zeros(n_col, dtype=np.int32)
        mt_umi_l = np.zeros(n_col, dtype=np.float64)
        crp_idx_parts, crp_dat_parts, crp_col_parts = [], [], []
        hist = np.zeros((len(crp_rows), 66), dtype=np.int64)

        data_ds = m["data"]
        indices_ds = m["indices"]
        for cs in range(0, n_col, 25000):
            ce = min(cs + 25000, n_col)
            p0, p1 = int(indptr_all[cs]), int(indptr_all[ce])
            if p1 == p0:
                continue
            dat = data_ds[p0:p1]
            idx = indices_ds[p0:p1]
            col_lens = np.diff(indptr_all[cs : ce + 1])
            col_ids = np.repeat(np.arange(cs, ce, dtype=np.int32), col_lens)

            m_r = is_rna_feat[idx]
            if np.any(m_r):
                r_dat = dat[m_r].astype(np.float64, copy=False)
                r_idx = idx[m_r]
                r_col = col_ids[m_r] - cs
                umi_l[cs:ce] = np.bincount(r_col, weights=r_dat, minlength=ce - cs)
                ng_l[cs:ce] = np.bincount(r_col, minlength=ce - cs)
                mt_m = mt_raw[r_idx]
                if np.any(mt_m):
                    mt_umi_l[cs:ce] = np.bincount(r_col[mt_m], weights=r_dat[mt_m], minlength=ce - cs)

            m_c = ~m_r
            if np.any(m_c):
                crp_idx_parts.append(crp_remap_idx[idx[m_c]])
                crp_dat_parts.append(dat[m_c].astype(np.int64, copy=False))
                crp_col_parts.append(col_ids[m_c])

    mt_l = 100.0 * mt_umi_l / np.maximum(umi_l, 1.0)
    ok_init = (mt_l < 20.0) & (ng_l >= 200)
    ok_hq = (mt_l <= mt_max_pct) & (ng_l >= 500)

    c_idx = np.concatenate(crp_idx_parts) if crp_idx_parts else np.empty(0, dtype=np.int32)
    c_dat = np.concatenate(crp_dat_parts) if crp_dat_parts else np.empty(0, dtype=np.int64)
    c_col = np.concatenate(crp_col_parts) if crp_col_parts else np.empty(0, dtype=np.int32)
    del indptr_all, crp_idx_parts, crp_dat_parts, crp_col_parts

    ok_init_mask = ok_init[c_col]
    bins = np.clip(c_dat[ok_init_mask], 1, 65)
    np.add.at(hist, (c_idx[ok_init_mask], bins), 1)
    thrs = fit_poisson_gaussian_thresholds_from_hist(hist, min_thr=5, max_thr=20, default_thr=7)

    hq_cols = np.where(ok_hq)[0]
    col_map_hq = np.full(n_col, -1, dtype=np.int32)
    col_map_hq[hq_cols] = np.arange(len(hq_cols), dtype=np.int32)
    hq_m = col_map_hq[c_col] >= 0
    C_hq = sp.csr_matrix(
        (c_dat[hq_m].astype(np.float64), (c_idx[hq_m], col_map_hq[c_col[hq_m]])),
        shape=(len(crp_rows), len(hq_cols)),
    ).tocsc()
    del c_idx, c_dat, c_col
    gc.collect()

    t1, t2, a1 = extract_top2_csc(C_hq)
    pos_mask = C_hq.data >= thrs[C_hq.indices]
    pos_csc = sp.csc_matrix((pos_mask.astype(np.int8), C_hq.indices, C_hq.indptr), shape=C_hq.shape)
    n_pos = np.asarray(pos_csc.sum(axis=0)).ravel()
    s1_mask = (n_pos == 1) & (t1 >= 5.0) & ((t2 < 3.0) | (t1 / np.maximum(t2, 0.5) >= 3.0))

    raw_keep: List[int] = []
    nt_flags: List[bool] = []
    g_pkeys: List[str] = []
    sublibs_k: List[str] = []
    pmeta_g: Dict[str, Tuple[str, str, str, str]] = {}
    n_safe_excl = 0

    for loc_j in np.where(s1_mask)[0]:
        gk = crp_keys[int(a1[loc_j])]
        raw_col = int(hq_cols[loc_j])
        if gk in ref_safe:
            n_safe_excl += 1
            continue
        if gk in ref_ntc:
            raw_keep.append(raw_col)
            nt_flags.append(True)
            g_pkeys.append("")
            sublibs_k.append(_parse_hcasmc_sublibrary(bcs[raw_col]))
        elif gk in guide_to_target:
            raw_keep.append(raw_col)
            nt_flags.append(False)
            g_pkeys.append(gk)
            if gk not in pmeta_g:
                pmeta_g[gk] = guide_to_target[gk]
            sublibs_k.append(_parse_hcasmc_sublibrary(bcs[raw_col]))

    # Check whether guide-level vs gene-level pooling is needed (count-only check)
    nt_arr = np.asarray(nt_flags, dtype=bool)
    keep_arr = np.asarray(raw_keep, dtype=np.int64)
    g_counts = Counter(pk for pk, nt in zip(g_pkeys, nt_flags) if not nt and pk)
    n_guides_25 = sum(v >= 25 for v in g_counts.values())
    sym_counts: Dict[str, int] = defaultdict(int)
    for pk, c in g_counts.items():
        sym_counts[pmeta_g[pk][1]] += c
    n_syms_25 = sum(v >= 25 for v in sym_counts.values())

    # If average guide count leaves < 150 guides at >=25 (e.g. 3 guides/gene dilution), pool P4 cells by gene
    if n_guides_25 < 150 and n_syms_25 > n_guides_25:
        pkeys = [pmeta_g[pk][1] if pk else "" for pk in g_pkeys]
        pmeta = {pmeta_g[pk][1]: pmeta_g[pk] for pk in g_pkeys if pk}
        target_key_mode = "gene_pooled_p4_cells"
    else:
        pkeys = g_pkeys
        pmeta = pmeta_g
        target_key_mode = "best_single_sgrna_guide"

    uniq_subs = sorted({s for s in sublibs_k if s})
    if len(uniq_subs) >= 2:
        sub_rank = {s: idx for idx, s in enumerate(uniq_subs, start=1)}
        half_flags = np.array([1 if (sub_rank.get(s, 1) % 2 == 1) else 2 for s in sublibs_k], dtype=np.int8)
        sub_ntc = Counter(s for s, nt in zip(sublibs_k, nt_arr) if nt)
        if sum(sub_ntc[s] for s in uniq_subs[0::2]) > 0 and sum(sub_ntc[s] for s in uniq_subs[1::2]) > 0:
            batches = np.array(
                [
                    (f"Half{'A' if sub_rank[s] % 2 == 1 else 'B'}_pooled" if sub_ntc[s] < 15 else f"sub_{s}")
                    for s in sublibs_k
                ],
                dtype=str,
            )
            split_mode = f"sorted_sublib_parity ({len(uniq_subs)} sublibraries: {uniq_subs})"
        else:
            half_flags = apply_stratified_seed42_split(pkeys, nt_arr, seed=42)
            batches = np.array([f"Half{'A' if h == 1 else 'B'}_strat42" for h in half_flags], dtype=str)
            split_mode = "missing_unit_stratified_seed42 (zero NTC in one sublibrary half)"
    else:
        half_flags = apply_stratified_seed42_split(pkeys, nt_arr, seed=42)
        batches = np.array([f"Half{'A' if h == 1 else 'B'}_strat42" for h in half_flags], dtype=str)
        split_mode = f"missing_unit_stratified_seed42 (parsed sublibrary suffixes={uniq_subs})"

    rna_ids = [fids[i] for i in rna_rows]
    rna_names = [fnames[i] for i in rna_rows]
    read_chunk, read_rows, read_chunk_cols = make_10x_csc_subset_readers(h5_p, keep_arr, n_feat, rna_rows)
    part = PartSpec(
        tag="HCASMC_hTERT",
        path=h5_p,
        n_cells=len(keep_arr),
        n_genes_raw=len(rna_rows),
        nt_mask=nt_arr,
        batches=batches,
        half_flags=half_flags,
        umi=umi_l[keep_arr],
        pert_keys=pkeys,
        pert_meta=pmeta,
        var_gene_ids=rna_ids,
        var_gene_names=rna_names,
        read_chunk_fn=read_chunk,
        read_rows_fn=read_rows,
        read_chunk_cols_fn=read_chunk_cols,
    )
    audit = {
        "n_raw_barcodes": int(n_col),
        "n_hq_cells": int(len(hq_cols)),
        "n_safe_targeting_excluded": int(n_safe_excl),
        "n_kept_cells": int(len(keep_arr)),
        "n_kept_ntc": int(nt_arr.sum()),
        "n_guides_ge25": int(n_guides_25),
        "n_syms_ge25": int(n_syms_25),
        "target_key_mode": target_key_mode,
        "split_half_mode": split_mode,
        "median_guide_threshold": float(np.median(thrs)),
    }
    return [part], audit


def build_cd4t_stim48_parts(raw_dir: Path, hgnc: HGNCLookup) -> Tuple[List[PartSpec], dict]:
    """
    CD4+ T Stim48hr CRISPRi (`GSE314342`), lexicographic donor split (`PREREG_D §2.2.3`):
      - Half A (`1`): Donor `CE0008678` (`D3_Stim48hr`, `CD4i_R2L25..48`, 24 lanes)
      - Half B (`2`): Donor `CE0010866` (`D2_Stim48hr`, `CD4i_R2L25..48`, 24 lanes)
    Uses pooled per-donor Poisson-Gaussian guide threshold (`>= 5` UMI floor + single-sgRNA
    ratio >= 3) and post-guide QC (`<= 5% MT & >= 500 genes`), identical to `cd4t_rest`.
    """
    cdir = raw_dir / "cd4t_stim48"
    # Pre-QC Q6 contamination guard
    forbidden = list(cdir.glob("*_D1_*")) + list(cdir.glob("*_D4_*")) + list(cdir.glob("*_Rest_*")) + list(cdir.glob("*_Stim8hr_*"))
    if forbidden:
        raise RuntimeError(f"Q6 overlap violation in cd4t_stim48 dir: {[p.name for p in forbidden[:5]]}")

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

    d3_files = sorted(cdir.glob("GSM93938*_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5"))
    d2_files = sorted(cdir.glob("GSM93938*_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5"))
    if len(d3_files) != 24 or len(d2_files) != 24:
        raise RuntimeError(f"Expected 24 D3 + 24 D2 Stim48hr files, found D3={len(d3_files)}, D2={len(d2_files)}")

    parts: List[PartSpec] = []
    donor_audit: Dict[str, dict] = {}

    for donor_id, half_val, files in [
        ("CE0008678_D3", 1, d3_files),
        ("CE0010866_D2", 2, d2_files),
    ]:
        print(f"  [cd4t_stim48] Scanning {donor_id} ({len(files)} lanes) for QC & pooled guide histogram...", flush=True)
        lane_cache = []
        hist = None
        crp_names_ref = None
        for lp in files:
            parts_name = lp.name.split("_")
            lane_tag = f"{parts_name[1]}_{parts_name[2]}"  # CD4i_R2L25 .. R2L48
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

            is_crp = np.zeros(n_feat, dtype=bool)
            is_crp[crp_rows] = True
            crp_remap_idx = np.full(n_feat, -1, dtype=np.int32)
            crp_remap_idx[crp_rows] = np.arange(len(crp_rows), dtype=np.int32)
            is_rna_feat = ~is_crp

            umi_l = np.zeros(n_col, dtype=np.float64)
            ng_l = np.zeros(n_col, dtype=np.int32)
            mt_umi_l = np.zeros(n_col, dtype=np.float64)
            crp_idx_parts, crp_dat_parts, crp_col_parts = [], [], []

            with h5py.File(lp, "r") as f:
                m = f["matrix"]
                indptr_all = m["indptr"][:]
                data_ds = m["data"]
                indices_ds = m["indices"]
                for cs in range(0, n_col, 25000):
                    ce = min(cs + 25000, n_col)
                    p0, p1 = int(indptr_all[cs]), int(indptr_all[ce])
                    if p1 == p0:
                        continue
                    dat = data_ds[p0:p1]
                    idx = indices_ds[p0:p1]
                    col_lens = np.diff(indptr_all[cs : ce + 1])
                    col_ids = np.repeat(np.arange(cs, ce, dtype=np.int32), col_lens)

                    m_r = is_rna_feat[idx]
                    if np.any(m_r):
                        r_dat = dat[m_r].astype(np.float64, copy=False)
                        r_idx = idx[m_r]
                        r_col = col_ids[m_r] - cs
                        umi_l[cs:ce] = np.bincount(r_col, weights=r_dat, minlength=ce - cs)
                        ng_l[cs:ce] = np.bincount(r_col, minlength=ce - cs)
                        mt_m = mt_raw[r_idx]
                        if np.any(mt_m):
                            mt_umi_l[cs:ce] = np.bincount(r_col[mt_m], weights=r_dat[mt_m], minlength=ce - cs)

                    m_c = ~m_r
                    if np.any(m_c):
                        crp_idx_parts.append(crp_remap_idx[idx[m_c]])
                        crp_dat_parts.append(dat[m_c].astype(np.int64, copy=False))
                        crp_col_parts.append(col_ids[m_c])

            mt_l = 100.0 * mt_umi_l / np.maximum(umi_l, 1.0)
            ok_init = (mt_l < 20.0) & (ng_l >= 200)
            ok_hq = (mt_l <= 5.0) & (ng_l >= 500)

            c_idx = np.concatenate(crp_idx_parts) if crp_idx_parts else np.empty(0, dtype=np.int32)
            c_dat = np.concatenate(crp_dat_parts) if crp_dat_parts else np.empty(0, dtype=np.int64)
            c_col = np.concatenate(crp_col_parts) if crp_col_parts else np.empty(0, dtype=np.int32)

            ok_init_mask = ok_init[c_col]
            bins = np.clip(c_dat[ok_init_mask], 1, 65)
            np.add.at(hist, (c_idx[ok_init_mask], bins), 1)

            hq_cols = np.where(ok_hq)[0]
            col_map_hq = np.full(n_col, -1, dtype=np.int32)
            col_map_hq[hq_cols] = np.arange(len(hq_cols), dtype=np.int32)
            hq_m = col_map_hq[c_col] >= 0
            C_hq = sp.csr_matrix(
                (c_dat[hq_m].astype(np.float64), (c_idx[hq_m], col_map_hq[c_col[hq_m]])),
                shape=(len(crp_rows), len(hq_cols)),
            ).tocsc()
            del indptr_all, c_idx, c_dat, c_col, crp_idx_parts, crp_dat_parts, crp_col_parts

            lane_cache.append(
                (lp, lane_tag, n_feat, rna_rows, [fids[i] for i in rna_rows], [fnames[i] for i in rna_rows], umi_l, hq_cols, C_hq)
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
            read_chunk, read_rows, read_chunk_cols = make_10x_csc_subset_readers(lp, keep_arr, n_feat, rna_rows)
            parts.append(
                PartSpec(
                    tag=lane_tag,
                    path=lp,
                    n_cells=len(keep_arr),
                    n_genes_raw=len(rna_rows),
                    nt_mask=nt_arr,
                    batches=np.full(len(keep_arr), f"{donor_id}_{lane_tag}"),
                    half_flags=np.full(len(keep_arr), half_val, dtype=np.int8),
                    umi=umi_l[keep_arr],
                    pert_keys=pkeys,
                    pert_meta=pmeta,
                    var_gene_ids=rids,
                    var_gene_names=rnames,
                    read_chunk_fn=read_chunk,
                    read_rows_fn=read_rows,
                    read_chunk_cols_fn=read_chunk_cols,
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

    return parts, {
        "donors": donor_audit,
        "half_rule": "HalfA=CE0008678_D3_Stim48hr, HalfB=CE0010866_D2_Stim48hr (lexicographic)",
        "shared_lineage_lab_disclosure": "Shares CD4+ alpha-beta T cell type, library, and Marson lab with R10 cd4t_rest",
    }


def _stream_feng_csv_to_csr_h5(
    csv_gz_path: Path,
    keep_cell_ids: Sequence[str],
    out_csr_h5: Path,
    cell_shard_size: int = 30000,
) -> Tuple[List[str], List[str], np.ndarray, np.ndarray, np.ndarray]:
    """
    Stream transposed `TargetedScreen_RNA-UMI-Counts.csv.gz` (`rows = genes`, `cols = cells`)
    line-by-line in < 1.5 GB RAM using temporary cell-range binary shards, retaining ONLY
    `keep_cell_ids` columns, and transpose into cell-major `CSR` HDF5 `(n_cells, n_genes)`.
    """
    cell_to_kidx = {cid: j for j, cid in enumerate(keep_cell_ids)}
    n_kcells = len(keep_cell_ids)

    gene_ids: List[str] = []
    gene_names: List[str] = []
    umi = np.zeros(n_kcells, dtype=np.float64)
    ng = np.zeros(n_kcells, dtype=np.int32)
    mt_umi = np.zeros(n_kcells, dtype=np.float64)

    n_shards = (n_kcells + cell_shard_size - 1) // cell_shard_size
    out_csr_h5.parent.mkdir(parents=True, exist_ok=True)
    shard_paths = [out_csr_h5.parent / f"_tmp_feng_shard_{s}.bin" for s in range(n_shards)]
    rec_dtype = np.dtype([("c", "<u2"), ("g", "<u4"), ("v", "<f4")])

    t0 = time.time()
    try:
        shard_fps = [open(sp_p, "wb", buffering=1024 * 1024) for sp_p in shard_paths]
        try:
            with gzip.open(csv_gz_path, "rt", encoding="utf-8") as fh:
                header_line = fh.readline()
                raw_cols = [c.strip().strip('"') for c in header_line.rstrip("\r\n").split(",")]
                hdr_cells = raw_cols[1:] if len(raw_cols) > 1 else raw_cols
                csv_pos_for_k = np.full(len(hdr_cells), -1, dtype=np.int32)
                found_keep = 0
                for ci, cid in enumerate(hdr_cells):
                    k = cell_to_kidx.get(cid, -1)
                    if k >= 0:
                        csv_pos_for_k[ci] = k
                        found_keep += 1
                if found_keep != n_kcells:
                    raise RuntimeError(f"Feng CSV header matched only {found_keep}/{n_kcells} keep_cell_ids!")

                src_idx = np.where(csv_pos_for_k >= 0)[0]
                dst_k = csv_pos_for_k[src_idx]
                shard_of_k = dst_k // cell_shard_size
                loc_of_k = (dst_k % cell_shard_size).astype(np.uint16)
                n_skipped_guide_rows = 0

                for line in fh:
                    comma = line.find(",")
                    if comma < 0:
                        continue
                    feat_label = line[:comma].strip().strip('"')
                    if feat_label.endswith(":CRISPR-Guide-Capture"):
                        n_skipped_guide_rows += 1
                        continue
                    if feat_label.endswith(":Gene-Expression"):
                        feat_core = feat_label[: -len(":Gene-Expression")]
                    else:
                        feat_core = feat_label

                    if ":" in feat_core:
                        ens_id, sym = feat_core.split(":", 1)
                        if not ens_id.startswith("ENSG") and sym.startswith("ENSG"):
                            ens_id, sym = sym, ens_id
                    else:
                        ens_id, sym = ("", feat_core) if not feat_core.startswith("ENSG") else (feat_core, feat_core)

                    arr = np.fromstring(line[comma + 1 :], dtype=np.int32, sep=",")
                    if len(arr) != len(hdr_cells):
                        raise RuntimeError(f"Row {feat_label} len={len(arr)} != {len(hdr_cells)}")
                    sub = arr[src_idx]
                    nz_mask = sub > 0
                    g_idx = len(gene_names)
                    gene_ids.append(ens_id.strip())
                    gene_names.append(sym.strip())

                    if np.any(nz_mask):
                        k_hit = dst_k[nz_mask]
                        v_hit = sub[nz_mask]
                        umi[k_hit] += v_hit
                        ng[k_hit] += 1
                        if sym.strip().startswith("MT-"):
                            mt_umi[k_hit] += v_hit

                        sh_hit = shard_of_k[nz_mask]
                        lc_hit = loc_of_k[nz_mask]
                        for sid in np.unique(sh_hit):
                            m_s = sh_hit == sid
                            rec = np.empty(int(m_s.sum()), dtype=rec_dtype)
                            rec["c"] = lc_hit[m_s]
                            rec["g"] = g_idx
                            rec["v"] = v_hit[m_s].astype(np.float32)
                            rec.tofile(shard_fps[int(sid)])

                    if (g_idx + 1) % 2000 == 0:
                        print(f"    [Feng stream] parsed {g_idx+1} genes ({time.time()-t0:.1f}s)", flush=True)
        finally:
            for fp in shard_fps:
                fp.close()

        n_genes = len(gene_names)
        print(
            f"    [Feng stream] {n_genes} RNA genes read ({n_skipped_guide_rows} guide capture rows skipped) "
            f"in {time.time()-t0:.1f}s; building chunked CSR HDF5 ({n_shards} cell shards)...",
            flush=True,
        )
        with h5py.File(out_csr_h5, "w") as hf:
            g = hf.create_group("X")
            ds_data = g.create_dataset("data", shape=(0,), maxshape=(None,), dtype="float32", chunks=(1 << 20,))
            ds_idx = g.create_dataset("indices", shape=(0,), maxshape=(None,), dtype="int32", chunks=(1 << 20,))
            indptr = np.zeros(n_kcells + 1, dtype=np.int64)
            cur_nnz = 0

            for sid, sp_p in enumerate(shard_paths):
                c0 = sid * cell_shard_size
                c1 = min(c0 + cell_shard_size, n_kcells)
                nc = c1 - c0
                rec = np.fromfile(sp_p, dtype=rec_dtype)
                sp_p.unlink(missing_ok=True)
                if len(rec) > 0:
                    csr_s = sp.csr_matrix(
                        (rec["v"], (rec["c"].astype(np.int32), rec["g"].astype(np.int32))),
                        shape=(nc, n_genes),
                        dtype=np.float32,
                    )
                    del rec
                    csr_s.sum_duplicates()
                    nnz_s = len(csr_s.data)
                    if nnz_s > 0:
                        ds_data.resize((cur_nnz + nnz_s,))
                        ds_idx.resize((cur_nnz + nnz_s,))
                        ds_data[cur_nnz : cur_nnz + nnz_s] = csr_s.data
                        ds_idx[cur_nnz : cur_nnz + nnz_s] = csr_s.indices.astype(np.int32)
                    indptr[c0 + 1 : c1 + 1] = np.int64(cur_nnz) + csr_s.indptr[1:].astype(np.int64)
                    cur_nnz += nnz_s
                    print(f"    [Feng shard] {sid+1}/{n_shards}: nnz_s={nnz_s:,}, cur_nnz={cur_nnz:,}", flush=True)
                    del csr_s
                else:
                    indptr[c0 + 1 : c1 + 1] = cur_nnz
                gc.collect()
            g.create_dataset("indptr", data=indptr, dtype=np.int64)
    finally:
        for sp_p in shard_paths:
            sp_p.unlink(missing_ok=True)

    pct_mt = 100.0 * mt_umi / np.maximum(umi, 1.0)
    return gene_ids, gene_names, umi, ng, pct_mt


def build_feng_ipsc_parts(raw_dir: Path, hgnc: HGNCLookup) -> Tuple[List[PartSpec], dict]:
    """
    Feng 2026 TargetedScreen (`feng_ipsc`, Secondary ONLY — H-D5 / SA-KOLF):
      - Excludes `kolf_2` & `kolf_3` (`Q6`) leaving 17 non-KOLF HipSci lines (`9` donors)
      - Strictly excludes all `"unassigned"` cells (`PREREG_D §2.2.4`)
      - Keeps `NonTarget_*` NTC (`7,266` cells) + single-sgRNA targeting (`550,406` cells)
      - Sorted `17` Cell_Line parity split (`9` odd vs `8` even); batches grouped at
        `Cell_Line x Inlet` if `>= 15` NTC else `Cell_Line_pooled` (`>= 86` NTC/line)
    """
    fdir = raw_dir / "feng_ipsc"
    meta_p = fdir / "TargetedScreen_Cell-Metadata.tsv.gz"
    csv_p = fdir / "TargetedScreen_RNA-UMI-Counts.csv.gz"
    for p in (meta_p, csv_p):
        if not p.exists():
            raise FileNotFoundError(p)

    df = pd.read_csv(meta_p, sep="\t", low_memory=False)
    n_raw_meta = len(df)
    n_kolf = int(df["Cell_Line"].isin(FENG_KOLF_EXCLUDED).sum())
    df = df[~df["Cell_Line"].isin(FENG_KOLF_EXCLUDED)].copy()
    assert not df["Cell_Line"].isin(FENG_KOLF_EXCLUDED).any(), "kolf_2/kolf_3 leak detected!"

    n_unassigned = int((df["Guide_Call"] == "unassigned").sum())
    df = df[df["Guide_Call"] != "unassigned"].copy().reset_index(drop=True)

    bparts = df["Batch"].astype(str).str.split("_", expand=True)
    df["Pool"] = bparts[1]
    df["Inlet"] = bparts[3]
    sorted_lines = sorted(df["Cell_Line"].unique())
    if len(sorted_lines) != 17:
        raise RuntimeError(f"Expected 17 non-KOLF lines in TargetedScreen, got {len(sorted_lines)}: {sorted_lines}")
    line_to_half = {cl: np.int8(1 if idx % 2 == 1 else 2) for idx, cl in enumerate(sorted_lines, start=1)}

    is_ntc = df["Guide_Call"].astype(str).str.startswith("NonTarget_").to_numpy()
    # Count NTC per Cell_Line x Inlet; pool to Cell_Line if < 15 NTC
    cl_inlet = [f"{cl}_{inl}" for cl, inl in zip(df["Cell_Line"], df["Inlet"])]
    unit_ntc_cnt = Counter(u for u, nt in zip(cl_inlet, is_ntc) if nt)
    batch_all = np.array(
        [(u if unit_ntc_cnt.get(u, 0) >= 15 else f"{cl}_pooled") for u, cl in zip(cl_inlet, df["Cell_Line"])],
        dtype=str,
    )

    pkeys: List[str] = []
    pmeta: Dict[str, Tuple[str, str, str, str]] = {}
    tmp_csr_h5 = R11_LINES_DIR / "feng_ipsc" / "_tmp_feng_csr.h5"
    # RNA-UMI-Counts columns omit `PC-` and `-D<day>` (e.g. `PC-P4-D3_I73_...` -> `P4_I73_...`)
    norm_cell_ids = [
        re.sub(r"^PC-(P\d+)-D\d+_(I\d+_.*)$", r"\1_\2", str(cid))
        for cid in df["Cell_ID"].astype(str)
    ]
    var_ids, var_names, umi_all, ng_all, pct_mt_all = _stream_feng_csv_to_csr_h5(
        csv_p, norm_cell_ids, tmp_csr_h5
    )

    # Apply conservative RNA sanity filter matching R10 (`MT < 20%` & `ng >= 500`)
    qc_ok = (pct_mt_all < 20.0) & (ng_all >= 500)
    keep_rows = np.where(qc_ok)[0]
    nt_arr = is_ntc[keep_rows]
    batches = batch_all[keep_rows]
    half_flags = np.array([line_to_half[cl] for cl in df["Cell_Line"].iloc[keep_rows]], dtype=np.int8)

    for row_i in keep_rows:
        gc_str = str(df["Guide_Call"].iloc[row_i]).strip()
        if gc_str.startswith("NonTarget_"):
            pkeys.append("")
            continue
        raw_sym = gc_str.split("_", 1)[0].strip()
        canon_sym, hstat = hgnc.resolve(raw_sym, "")
        ens = hgnc.sym_to_ens.get(canon_sym, "")
        pkeys.append(gc_str)
        if gc_str not in pmeta:
            pmeta[gc_str] = (ens, canon_sym, raw_sym, hstat)

    read_chunk, read_rows, read_chunk_cols = make_csr_h5_subset_readers(tmp_csr_h5, "X", keep_rows, len(var_names))
    part = PartSpec(
        tag="Feng2026_TargetedScreen",
        path=tmp_csr_h5,
        n_cells=len(keep_rows),
        n_genes_raw=len(var_names),
        nt_mask=nt_arr,
        batches=batches,
        half_flags=half_flags,
        umi=umi_all[keep_rows],
        pert_keys=pkeys,
        pert_meta=pmeta,
        var_gene_ids=var_ids,
        var_gene_names=var_names,
        read_chunk_fn=read_chunk,
        read_rows_fn=read_rows,
        read_chunk_cols_fn=read_chunk_cols,
    )
    audit = {
        "n_metadata_total": int(n_raw_meta),
        "n_kolf_excluded": int(n_kolf),
        "n_unassigned_excluded": int(n_unassigned),
        "n_qc_kept_cells": int(len(keep_rows)),
        "n_kept_ntc": int(nt_arr.sum()),
        "n_non_kolf_lines": len(sorted_lines),
        "half_a_lines": sorted_lines[0::2],
        "half_b_lines": sorted_lines[1::2],
        "tmp_csr_h5": str(tmp_csr_h5),
    }
    return [part], audit


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Prepare R11 standardized shared cell lines (A-T1).")
    p.add_argument(
        "--line",
        required=True,
        choices=["telohaec", "hcasmc", "cd4t_stim48", "feng_ipsc", "all"],
        help="Which R11 line(s) to process.",
    )
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--chunk-size", type=int, default=8192)
    p.add_argument("--raw-dir", type=Path, default=R11_RAW_DIR)
    p.add_argument("--out-root", type=Path, default=R11_LINES_DIR)
    return p.parse_args()


def run_line(line: str, raw_dir: Path, out_root: Path, hgnc: HGNCLookup, workers: int, chunk_size: int = 8192) -> Path:
    require_r11_go()
    builders = {
        "telohaec": build_telohaec_parts,
        "hcasmc": build_hcasmc_parts,
        "cd4t_stim48": build_cd4t_stim48_parts,
        "feng_ipsc": build_feng_ipsc_parts,
    }
    parts, audit = builders[line](raw_dir, hgnc)
    tmp_to_clean = audit.pop("tmp_csr_h5", None)
    try:
        return process_r10_line(
            line_name=line,
            parts=parts,
            hgnc=hgnc,
            chunk_size=chunk_size,
            out_root=out_root,
            extra_meta={"adapter_audit": audit},
            n_workers=workers,
        )
    finally:
        if tmp_to_clean:
            Path(tmp_to_clean).unlink(missing_ok=True)


def main() -> None:
    args = parse_args()
    require_r11_go()
    hgnc = build_hgnc_lookup(DEFAULT_HGNC_PATH)
    targets = ["telohaec", "hcasmc", "cd4t_stim48", "feng_ipsc"] if args.line == "all" else [args.line]
    for ln in targets:
        run_line(ln, args.raw_dir, args.out_root, hgnc, args.workers, args.chunk_size)


if __name__ == "__main__":
    main()
