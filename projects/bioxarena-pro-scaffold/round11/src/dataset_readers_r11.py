#!/usr/bin/env python3
"""
Round 11 Standardized Dataset Readers (`PREREG_D v2` & `R11_DOWNLOAD_PLAN.md`).

Implements four bounded-memory, zero-expression-leakage ingestion adapters:
  1. `build_telohaec_parts`: `.h5ad` + IGVF guide TSV + GEO Dialout (`max_umi_ct >= 4 & X2nd_place_umi_ct < 4`),
     strict `type == "non-targeting"` NTC-only (`safe-targeting` excluded), `-TSS*` / dual-locus exclusion.
  2. `check_hcasmc_feasibility` & `build_hcasmc_parts`: Parse v1.2.1 CSC `.h5` chunked streaming reader
     (`<= 300 MB` slice footprint), strict `non-targeting` vs `safe-targeting` (`0Safe_`) exclusion,
     Poisson-Gaussian `P4` guide assignment, `__<sublib>` parity split + `seed=42` backup.
  3. `build_cd4t_stim48_parts`: Reuses Round 10 `cd4t_rest` Flex `.h5` Poisson-Gaussian pipeline on
     donor/timepoint-disjoint `D3_Stim48hr` (`Half A = CE0008678`) vs `D2_Stim48hr` (`Half B = CE0010866`).
  4. `build_feng_ipsc_parts`: Pure row-by-row gzip streaming of `TargetedScreen_RNA-UMI-Counts.csv.gz`
     into cell-bucketed CSC HDF5 cache, `Q6` `kolf_2`/`kolf_3` exclusion, `"unassigned"` exclusion,
     `NonTarget_*` NTC-only, 17-line parity split, and `< 15` NTC `Cell_Line x Inlet` fallback.

Imports shared preprocessing helpers cleanly from frozen `round10/shared/scripts/prepare_r10_lines.py`
and `round10/src/data_loader_r10.py` without modifying any Round 10 file.
"""

from __future__ import annotations

import csv
import gzip
import re
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROUND11_SRC_DIR = Path(__file__).resolve().parent
REPO_ROOT = ROUND11_SRC_DIR.parents[1]
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
R10_SCRIPTS_DIR = REPO_ROOT / "round10" / "shared" / "scripts"
for p in (str(R10_SRC_DIR), str(R10_SCRIPTS_DIR), str(ROUND11_SRC_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from data_loader_r10 import decode_str_array, iter_h5ad_csr_chunks_r10  # noqa: E402
from prepare_r10_lines import (  # noqa: E402
    CD4_EXCEL_GUIDE_REMAP,
    HGNCLookup,
    PartSpec,
    build_hgnc_lookup,
    compute_10x_csc_rna_qc,
    extract_top2_csc,
    fit_poisson_gaussian_thresholds_from_hist,
    make_10x_csc_subset_readers,
    open_x_node,
    process_r10_line,
    read_h5ad_chunk_dense,
    read_obs_or_var_col,
)

_TSS_SUFFIX_RE = re.compile(r"-TSS\d*$", re.IGNORECASE)
_TELOHAEC_LANE_RE = re.compile(r"-(\d+)$")
_PARSE_SUBLIB_RE = re.compile(r"__([A-Za-z0-9]+)$")
_FENG_KOLF_EXCLUDED = {"kolf_2", "kolf_3", "kolf2", "kolf3", "kolf_2_1j"}
_FENG_DELIM_RE = re.compile(r"[,;|+]")


def _find_first_file(directory: Path, preferred: Sequence[str], glob_pattern: str) -> Path:
    for name in preferred:
        cand = directory / name
        if cand.exists():
            return cand
    matches = sorted(directory.glob(glob_pattern))
    if not matches:
        raise FileNotFoundError(f"Missing expected file matching {glob_pattern} in {directory}")
    return matches[0]


def _pick_column(columns: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
    low_map = {str(c).lower(): str(c) for c in columns}
    for cand in candidates:
        if cand in columns:
            return cand
        if cand.lower() in low_map:
            return low_map[cand.lower()]
    return None


def _is_plain_single_gene_symbol(sym: str) -> bool:
    s = str(sym).strip()
    if not s or s.lower() in ("nan", "none", "na", "null"):
        return False
    if _TSS_SUFFIX_RE.search(s):
        return False
    low = s.lower()
    if (
        low.startswith("enhancer")
        or low.startswith("9p21")
        or "-and-" in low
        or "negative" in low
        or "safe" in low
        or "non-targeting" in low
        or "nontargeting" in low
    ):
        return False
    if any(ch in s for ch in ("_", ",", ";", "/", "+", "|", " ")):
        return False
    return True


def _apply_min_ntc_batch_fallback(
    batches: np.ndarray,
    nt_mask: np.ndarray,
    half_flags: np.ndarray,
    group_labels: Optional[np.ndarray] = None,
    min_batch_nt: int = 15,
) -> Tuple[np.ndarray, int]:
    """
    Merge batches with `< min_batch_nt` NTCs into per-group (if group >= min_batch_nt NTCs)
    or half-pooled NTC fallback bucket per `PREREG_D §2.3`.
    """
    out = np.asarray(batches, dtype=object).copy()
    nt_counts: Dict[str, int] = defaultdict(int)
    for b, is_nt in zip(out, nt_mask):
        if is_nt:
            nt_counts[str(b)] += 1

    group_nt_counts: Dict[str, int] = defaultdict(int)
    if group_labels is not None:
        for g, is_nt in zip(group_labels, nt_mask):
            if is_nt:
                group_nt_counts[str(g)] += 1

    n_fallback = 0
    for i in range(len(out)):
        b_str = str(out[i])
        if nt_counts.get(b_str, 0) >= min_batch_nt:
            continue
        n_fallback += 1
        h = int(half_flags[i])
        if group_labels is not None:
            g_str = str(group_labels[i])
            if group_nt_counts.get(g_str, 0) >= min_batch_nt:
                out[i] = f"{g_str}_half{h}_pooled"
                continue
        out[i] = f"half_{h}_pooled"

    return out.astype(str), n_fallback


def _compute_h5ad_feature_subsets(f: h5py.File) -> Tuple[np.ndarray, np.ndarray, List[str], List[str], List[str]]:
    """Inspect `var` group cleanly to separate RNA features and extract IDs/symbols."""
    var = f["var"]
    var_keys = set(var.keys())
    idx_col = "_index" if "_index" in var_keys else ("index" if "index" in var_keys else None)
    raw_idx = read_obs_or_var_col(var, idx_col) if idx_col else [str(i) for i in range(int(f["X"].attrs["shape"][1]) if isinstance(f["X"], h5py.Group) else int(f["X"].shape[1]))]
    n_var = len(raw_idx)

    ft_col = next((c for c in ("feature_types", "feature_type") if c in var_keys), None)
    if ft_col is not None:
        ftypes = [str(x) for x in read_obs_or_var_col(var, ft_col)]
        rna_rows = np.array([i for i, t in enumerate(ftypes) if t == "Gene Expression"], dtype=np.int64)
        if len(rna_rows) == 0:
            rna_rows = np.arange(n_var, dtype=np.int64)
    else:
        ftypes = ["Gene Expression"] * n_var
        rna_rows = np.arange(n_var, dtype=np.int64)

    id_col = next((c for c in ("ensembl_id", "gene_ids", "gene_id", "feature_id") if c in var_keys), idx_col)
    name_col = next((c for c in ("gene_symbol", "gene_symbols", "feature_name", "gene_name") if c in var_keys), idx_col)
    var_ids_all = read_obs_or_var_col(var, id_col) if id_col else raw_idx
    var_names_all = read_obs_or_var_col(var, name_col) if name_col else raw_idx

    r_ids = [str(var_ids_all[i]).strip() for i in rna_rows]
    r_names = [str(var_names_all[i]).strip() for i in rna_rows]
    return rna_rows, np.array([gn.startswith("MT-") for gn in r_names], dtype=bool), r_ids, r_names, ftypes


# ==============================================================================
# 1. TeloHAEC_2285 Reader (.h5ad + IGVF guide TSV + GEO Dialout)
# ==============================================================================


def build_telohaec_parts(
    raw_dir: Path,
    hgnc: HGNCLookup,
    max_pct_mt: float = 20.0,
    min_ngenes: int = 500,
    min_umi_p4: int = 4,
    include_tss_sensitivity: bool = False,
) -> Tuple[List[PartSpec], dict]:
    """
    TeloHAEC_2285 (`IGVFDS7943IIWZ` / `GSE210681`):
      - Single-sgRNA P4 lock: `max_umi_ct >= 4 & X2nd_place_umi_ct < 4`
      - Primary NTC ONLY (`type == "non-targeting"` / `negative_control__*`), strictly
        excluding `safe-targeting` (`safe_targeting__*`) even when `targeting == True`
      - Plain HGNC-approved single-gene targets only (excluding `-TSS*`, `Enhancer-*`, `9p21-*`, `*-and-*`)
      - Odd sorted 10x lanes (`1,3,...,19`) -> Half A (`1`); Even lanes (`2,4,...,20`) -> Half B (`2`).
    """
    tdir = raw_dir / "telohaec"
    h5ad_path = _find_first_file(tdir, ["IGVFFI8645YMGX.h5ad"], "*.h5ad")
    tsv_path = _find_first_file(tdir, ["IGVFFI9754AGFB.tsv.gz", "IGVFFI9754AGFB.tsv"], "*9754*.tsv*")
    dial_path = _find_first_file(
        tdir,
        ["GSE210681_Dialout_guides_per_CBC.txt.gz", "GSE210681_Dialout_guides_per_CBC.txt"],
        "*Dialout*.txt*",
    )

    guide_tsv = pd.read_csv(tsv_path, sep="\t", dtype=str, keep_default_na=False)
    gid_col = _pick_column(guide_tsv.columns, ["guide_id", "feature_id", "sgrna_id", "id"])
    gtype_col = _pick_column(guide_tsv.columns, ["type", "guide_type", "category"])
    gtarget_col = _pick_column(guide_tsv.columns, ["intended_target_name", "gene_symbol", "target_gene_name"])

    # guide_id -> (type_lower, raw_target)
    guide_type_map: Dict[str, Tuple[str, str]] = {}
    if gid_col is not None:
        for row in guide_tsv.itertuples(index=False):
            gid = str(getattr(row, gid_col)).strip()
            gtype = str(getattr(row, gtype_col)).strip().lower() if gtype_col else ""
            gtgt = str(getattr(row, gtarget_col)).strip() if gtarget_col else ""
            if gid:
                guide_type_map[gid] = (gtype, gtgt)

    with (gzip.open(dial_path, "rt") if str(dial_path).endswith(".gz") else open(dial_path, "r")) as f_hdr:
        first_line = f_hdr.readline()
    sep_dial = "\t" if "\t" in first_line else ("," if "," in first_line else r"\s+")
    dial_df = pd.read_csv(dial_path, sep=sep_dial, low_memory=False)

    cbc_col = _pick_column(dial_df.columns, ["CBC", "cell_barcode", "barcode", "cell_id"]) or dial_df.columns[0]
    m1_col = _pick_column(dial_df.columns, ["max_umi_ct", "max_umi", "top1_umi"])
    m2_col = _pick_column(dial_df.columns, ["X2nd_place_umi_ct", "2nd_place_umi_ct", "second_place_umi_ct", "top2_umi"])
    dg_col = _pick_column(dial_df.columns, ["guide_id", "top_guide", "best_guide", "sgrna", "guide"])
    ds_col = _pick_column(dial_df.columns, ["gene_symbol", "target_gene", "intended_target_name", "gene"])
    dlane_col = _pick_column(dial_df.columns, ["lane", "channel", "sample"])
    if m1_col is None or m2_col is None:
        raise KeyError(f"Dialout table missing max/2nd UMI columns in {dial_path}: {list(dial_df.columns)}")

    cbc_list = dial_df[cbc_col].astype(str).str.strip().tolist()
    m1_vals = pd.to_numeric(dial_df[m1_col], errors="coerce").to_numpy(dtype=np.float64)
    m2_vals = pd.to_numeric(dial_df[m2_col], errors="coerce").to_numpy(dtype=np.float64)
    dg_list = dial_df[dg_col].astype(str).str.strip().tolist() if dg_col else [""] * len(dial_df)
    ds_list = dial_df[ds_col].astype(str).str.strip().tolist() if ds_col else [""] * len(dial_df)
    dlane_raw = dial_df[dlane_col].tolist() if dlane_col is not None else [None] * len(dial_df)

    cbc_map: Dict[str, Tuple[str, str, Optional[int]]] = {}
    n_p4_pass = n_safe_excluded = n_tss_excluded = n_nonplain_excluded = n_ntc_dial = n_tgt_dial = 0
    for cbc, m1, m2, guide_str, sym_raw, lane_obj in zip(
        cbc_list, m1_vals, m2_vals, dg_list, ds_list, dlane_raw
    ):
        if not (np.isfinite(m1) and np.isfinite(m2) and m1 >= float(min_umi_p4) and m2 < float(min_umi_p4)):
            continue
        n_p4_pass += 1
        lane_val = None
        if lane_obj is not None:
            try:
                fv = float(lane_obj)
                if np.isfinite(fv):
                    lane_val = int(fv)
            except Exception:
                lane_val = None

        gtype, tsv_tgt = guide_type_map.get(guide_str, ("", ""))
        sym_use = sym_raw or tsv_tgt
        low_gid = guide_str.lower()
        low_sym = sym_use.lower()

        # Safe-targeting MUST be excluded cleanly before NTC check
        if (
            "safe" in gtype
            or low_gid.startswith("safe_targeting")
            or low_sym in ("safe-targeting", "safe_targeting")
        ):
            n_safe_excluded += 1
            continue

        if (
            gtype == "non-targeting"
            or low_gid.startswith("negative_control")
            or low_sym in ("negative-control", "negative_control", "non-targeting", "nontargeting", "ntc")
        ):
            cbc_map[cbc] = ("ntc", "", lane_val)
            n_ntc_dial += 1
            continue

        if _TSS_SUFFIX_RE.search(sym_use):
            if not include_tss_sensitivity:
                n_tss_excluded += 1
                continue
            sym_use = _TSS_SUFFIX_RE.sub("", sym_use)

        if not _is_plain_single_gene_symbol(sym_use):
            n_nonplain_excluded += 1
            continue

        canon_sym, status = hgnc.resolve(sym_use, "")
        if canon_sym not in hgnc.approved:
            n_nonplain_excluded += 1
            continue

        n_tgt_dial += 1
        cbc_map[cbc] = ("tgt", f"{canon_sym}||{sym_use}||{status}", lane_val)

    # Read `.h5ad` metadata & compute RNA QC chunk-by-chunk
    f = h5py.File(h5ad_path, "r")
    try:
        x_node, shape = open_x_node(f, x_layer="X")
        n_raw_cells, n_raw_var = int(shape[0]), int(shape[1])
        obs = f["obs"]
        obs_keys = set(obs.keys())
        obs_idx_col = "_index" if "_index" in obs_keys else "index"
        obs_barcodes = read_obs_or_var_col(obs, obs_idx_col)
        obs_lane_col = next((c for c in ("lane", "channel", "gem_group") if c in obs_keys), None)
        obs_lanes = read_obs_or_var_col(obs, obs_lane_col) if obs_lane_col else None

        rna_rows, mt_mask_rna, rna_ids, rna_names, _ = _compute_h5ad_feature_subsets(f)
        full_rna_subset = (len(rna_rows) == n_raw_var and np.array_equal(rna_rows, np.arange(n_raw_var)))

        # Compute RNA QC streaming chunks of X cleanly
        umi_all = np.zeros(n_raw_cells, dtype=np.float64)
        ng_all = np.zeros(n_raw_cells, dtype=np.int32)
        mt_umi_all = np.zeros(n_raw_cells, dtype=np.float64)

        if isinstance(x_node, h5py.Group) and "indptr" in x_node:
            enc_qc = str(x_node.attrs.get("encoding-type", "")).lower()
            if "csc" in enc_qc:
                is_rna_var = np.zeros(n_raw_var, dtype=bool)
                is_rna_var[rna_rows] = True
                is_mt_var = np.zeros(n_raw_var, dtype=bool)
                if np.any(mt_mask_rna):
                    is_mt_var[rna_rows[mt_mask_rna]] = True
                indptr_full = np.asarray(x_node["indptr"][:], dtype=np.int64)
                for cs in range(0, n_raw_var, 2048):
                    ce = min(cs + 2048, n_raw_var)
                    p0, p1 = int(indptr_full[cs]), int(indptr_full[ce])
                    if p1 <= p0:
                        continue
                    idx_seg = np.asarray(x_node["indices"][p0:p1], dtype=np.int64)
                    dat_seg = np.asarray(x_node["data"][p0:p1], dtype=np.float64)
                    col_off = np.repeat(np.arange(cs, ce, dtype=np.int64), np.diff(indptr_full[cs : ce + 1]))
                    m_rna = is_rna_var[col_off]
                    if np.any(m_rna):
                        r_idx = idx_seg[m_rna]
                        d_val = dat_seg[m_rna]
                        umi_all += np.bincount(r_idx, weights=d_val, minlength=n_raw_cells)
                        pos = d_val > 0
                        if np.any(pos):
                            ng_all += np.bincount(r_idx[pos], minlength=n_raw_cells)
                    m_mt = is_mt_var[col_off]
                    if np.any(m_mt):
                        mt_umi_all += np.bincount(idx_seg[m_mt], weights=dat_seg[m_mt], minlength=n_raw_cells)
            else:
                for cs, ce, blk in iter_h5ad_csr_chunks_r10(h5ad_path, row_chunk=4096):
                    sub = blk if full_rna_subset else blk[:, rna_rows]
                    umi_all[cs:ce] = sub.sum(axis=1)
                    ng_all[cs:ce] = np.sum(sub > 0, axis=1)
                    if np.any(mt_mask_rna):
                        mt_umi_all[cs:ce] = sub[:, mt_mask_rna].sum(axis=1)
        elif sp.issparse(x_node):
            for cs in range(0, n_raw_cells, 4096):
                ce = min(cs + 4096, n_raw_cells)
                blk = x_node[cs:ce, :]
                sub = blk if full_rna_subset else blk[:, rna_rows]
                umi_all[cs:ce] = np.asarray(sub.sum(axis=1), dtype=np.float64).ravel()
                ng_all[cs:ce] = np.asarray((sub > 0).sum(axis=1), dtype=np.int32).ravel()
                if np.any(mt_mask_rna):
                    mt_umi_all[cs:ce] = np.asarray(sub[:, mt_mask_rna].sum(axis=1), dtype=np.float64).ravel()
        else:
            for cs in range(0, n_raw_cells, 4096):
                ce = min(cs + 4096, n_raw_cells)
                blk = np.asarray(x_node[cs:ce, :], dtype=np.float64)
                sub = blk if full_rna_subset else blk[:, rna_rows]
                umi_all[cs:ce] = sub.sum(axis=1)
                ng_all[cs:ce] = np.sum(sub > 0, axis=1)
                if np.any(mt_mask_rna):
                    mt_umi_all[cs:ce] = sub[:, mt_mask_rna].sum(axis=1)
        pct_mt_all = 100.0 * mt_umi_all / np.maximum(umi_all, 1.0)

        keep_rows: List[int] = []
        nt_flags: List[bool] = []
        batches: List[str] = []
        half_flags: List[int] = []
        pkeys: List[str] = []
        pmeta: Dict[str, Tuple[str, str, str, str]] = {}

        for r_i, bc in enumerate(obs_barcodes):
            if not (umi_all[r_i] > 0.0 and ng_all[r_i] >= min_ngenes and pct_mt_all[r_i] < max_pct_mt):
                continue
            call = cbc_map.get(bc)
            if call is None and "-" in bc:
                call = cbc_map.get(bc.split("-")[0])
            if call is None:
                continue
            kind, payload, lane_d = call
            lane_num = lane_d
            if lane_num is None:
                m_l = _TELOHAEC_LANE_RE.search(bc)
                if m_l:
                    lane_num = int(m_l.group(1))
                elif obs_lanes is not None:
                    m_ol = re.search(r"(\d+)$", str(obs_lanes[r_i]))
                    if m_ol:
                        lane_num = int(m_ol.group(1))
            if lane_num is None:
                continue

            keep_rows.append(r_i)
            batches.append(f"lane_{lane_num:02d}")
            # PREREG_D §2.2.1: Half A = odd sorted lanes (1,3,...,19); Half B = even sorted lanes (2,4,...,20)
            half_flags.append(1 if (lane_num % 2 == 1) else 2)
            if kind == "ntc":
                nt_flags.append(True)
                pkeys.append("")
            else:
                nt_flags.append(False)
                canon_sym, raw_sym, status = payload.split("||", 2)
                pkeys.append(canon_sym)
                if canon_sym not in pmeta:
                    pmeta[canon_sym] = (hgnc.sym_to_ens.get(canon_sym, ""), canon_sym, raw_sym, status)

        keep_arr = np.asarray(keep_rows, dtype=np.int64)
        nt_arr = np.asarray(nt_flags, dtype=bool)
        half_arr = np.asarray(half_flags, dtype=np.int8)
        batch_arr, n_fb = _apply_min_ntc_batch_fallback(np.array(batches, dtype=str), nt_arr, half_arr, min_batch_nt=15)

        if len(keep_arr) == 0:
            X_kept = sp.csr_matrix((0, len(rna_rows)), dtype=np.float32)
        elif isinstance(x_node, h5py.Group):
            enc = str(x_node.attrs.get("encoding-type", "")).lower()
            if "csc" in enc:
                x_full = sp.csc_matrix(
                    (x_node["data"][:], x_node["indices"][:], x_node["indptr"][:]),
                    shape=(n_raw_cells, n_raw_var),
                ).tocsr()
            else:
                x_full = sp.csr_matrix(
                    (x_node["data"][:], x_node["indices"][:], x_node["indptr"][:]),
                    shape=(n_raw_cells, n_raw_var),
                )
            X_kept = x_full[keep_arr][:, rna_rows].tocsr() if not full_rna_subset else x_full[keep_arr].tocsr()
        else:
            X_kept = x_node[keep_arr][:, rna_rows] if not full_rna_subset else x_node[keep_arr]
    finally:
        f.close()

    n_genes_rna = len(rna_rows)
    part = PartSpec(
        tag="TeloHAEC_2285",
        path=h5ad_path,
        n_cells=len(keep_arr),
        n_genes_raw=n_genes_rna,
        nt_mask=nt_arr,
        batches=batch_arr,
        half_flags=half_arr,
        umi=umi_all[keep_arr],
        pert_keys=pkeys,
        pert_meta=pmeta,
        var_gene_ids=rna_ids,
        var_gene_names=rna_names,
        read_chunk_fn=lambda s, e, mat=X_kept: read_h5ad_chunk_dense(mat, s, e, n_genes_rna),
    )
    audit = {
        "n_raw_cells_h5ad": n_raw_cells,
        "n_dialout_p4_pass": int(n_p4_pass),
        "n_safe_targeting_excluded": int(n_safe_excluded),
        "n_tss_excluded": int(n_tss_excluded),
        "n_nonplain_or_unapproved_excluded": int(n_nonplain_excluded),
        "n_dialout_p4_ntc": int(n_ntc_dial),
        "n_dialout_p4_targeting": int(n_tgt_dial),
        "n_qc_kept_cells": int(len(keep_arr)),
        "n_kept_ntc": int(nt_arr.sum()),
        "n_kept_half_a_ntc": int(np.sum(nt_arr & (half_arr == 1))),
        "n_kept_half_b_ntc": int(np.sum(nt_arr & (half_arr == 2))),
        "n_low_ntc_batch_fallback_cells": int(n_fb),
        "half_rule": "odd_sorted_10x_lanes_vs_even_sorted_10x_lanes",
    }
    return [part], audit


# ==============================================================================
# 2. HCASMC-hTERT Parse v1.2.1 Chunked CSC HDF5 Reader & Feasibility Gate (Q3)
# ==============================================================================


def _load_hcasmc_guide_metadata(raw_dir: Path, hgnc: HGNCLookup) -> Tuple[Set[str], Set[str], Dict[str, Tuple[str, str, str, str]]]:
    """
    Parse IGVF guide metadata TSVs (`IGVFFI4575UMXX.tsv.gz` / `IGVFFI9546UECF.tsv.gz`).
    CRITICAL (`PREREG_D §2.2.2`): BOTH `non-targeting` AND `safe-targeting` (`0Safe_`) carry
    `targeting == False`; identify NTC strictly via `type == "non-targeting"`.
    """
    hdir = raw_dir / "hcasmc"
    tsv_files = sorted(hdir.glob("*.tsv*"))
    if not tsv_files:
        raise FileNotFoundError(f"No guide metadata .tsv/.tsv.gz found in {hdir}")

    ntc_guides: Set[str] = set()
    safe_guides: Set[str] = set()
    guide_to_target: Dict[str, Tuple[str, str, str, str]] = {}

    for tpath in tsv_files:
        df = pd.read_csv(tpath, sep="\t", dtype=str, keep_default_na=False)
        gid_col = _pick_column(df.columns, ["guide_id", "feature_id", "id", "sgrna_id"]) or df.columns[0]
        gtype_col = _pick_column(df.columns, ["type", "guide_type", "category"])
        gtgt_col = _pick_column(df.columns, ["intended_target_name", "gene_symbol", "target_gene_name"])
        for r in df.itertuples(index=False):
            gid = str(getattr(r, gid_col)).strip()
            if not gid:
                continue
            gtype = str(getattr(r, gtype_col)).strip().lower() if gtype_col else ""
            raw_sym = str(getattr(r, gtgt_col)).strip() if gtgt_col else ""
            low_gid = gid.lower()
            if "safe" in gtype or low_gid.startswith("0safe_") or low_gid.startswith("safe"):
                safe_guides.add(gid)
                continue
            if gtype == "non-targeting" or low_gid.startswith("non-targeting") or low_gid.startswith("ntc"):
                ntc_guides.add(gid)
                continue
            if raw_sym:
                canon_sym, status = hgnc.resolve(raw_sym, "")
                if canon_sym in hgnc.approved:
                    guide_to_target[gid] = (hgnc.sym_to_ens.get(canon_sym, ""), canon_sym, raw_sym, status)
    return ntc_guides, safe_guides, guide_to_target


def check_hcasmc_feasibility(
    raw_dir: Path,
    hgnc: Optional[HGNCLookup] = None,
    batch_cols: int = 20000,
) -> Dict[str, Any]:
    """
    `PREREG_D §2.2.2 Q3` Mandatory Pre-Evaluation Technical Feasibility Gate:
      1) Verify `/matrix` Parse CSC HDF5 layout is chunk-sliceable in `batch_cols` windows;
      2) Verify `feature_type != "Gene Expression"` guide-capture rows exist and overlap TSV metadata.
    Never computes any expression correlation or perturbation summary (`counts/structure ONLY`).
    """
    try:
        hgnc_l = hgnc or build_hgnc_lookup()
        hdir = raw_dir / "hcasmc"
        h5_path = _find_first_file(hdir, ["IGVFFI5999NXAF.h5"], "*.h5*")
        ntc_guides, safe_guides, guide_to_target = _load_hcasmc_guide_metadata(raw_dir, hgnc_l)

        with h5py.File(h5_path, "r") as f:
            if "matrix" not in f:
                return {"feasible": False, "reason": "Missing /matrix group in HDF5"}
            m = f["matrix"]
            for req in ("shape", "barcodes", "data", "indices", "indptr", "features"):
                if req not in m:
                    return {"feasible": False, "reason": f"Missing /matrix/{req} dataset"}
            n_feat, n_col = int(m["shape"][0]), int(m["shape"][1])
            ftypes = np.array(decode_str_array(m["features/feature_type"][:]))
            fnames = np.array(decode_str_array(m["features/name"][:]))
            fids = np.array(decode_str_array(m["features/id"][:]))
            rna_rows = np.where(ftypes == "Gene Expression")[0]
            crp_rows = np.where(ftypes != "Gene Expression")[0]
            if len(rna_rows) == 0 or len(crp_rows) == 0:
                return {
                    "feasible": False,
                    "reason": f"Missing RNA ({len(rna_rows)}) or guide-capture ({len(crp_rows)}) feature_type rows",
                }
            # Test chunk sliceability of indptr/data/indices
            ce = min(batch_cols, n_col)
            ptrs = m["indptr"][0 : ce + 1]
            p0, p1 = int(ptrs[0]), int(ptrs[-1])
            _ = m["data"][p0 : min(p1, p0 + 1000)]
            _ = m["indices"][p0 : min(p1, p0 + 1000)]

            bcs_sample = decode_str_array(m["barcodes"][0 : min(5000, n_col)])
            sublibs = sorted({mm.group(1) for bc in bcs_sample if (mm := _PARSE_SUBLIB_RE.search(bc))})
            crp_overlap_ntc = sum(1 for i in crp_rows if fids[i] in ntc_guides or fnames[i] in ntc_guides)
            crp_overlap_tgt = sum(1 for i in crp_rows if fids[i] in guide_to_target or fnames[i] in guide_to_target)

        return {
            "feasible": bool(len(crp_rows) > 0 and (crp_overlap_ntc > 0 or len(ntc_guides) > 0)),
            "reason": "OK" if crp_overlap_ntc > 0 or len(ntc_guides) > 0 else "No NTC guide features matched",
            "n_features": n_feat,
            "n_cells": n_col,
            "n_rna_features": int(len(rna_rows)),
            "n_guide_features": int(len(crp_rows)),
            "n_sample_sublibs": len(sublibs),
            "sample_sublibs": sublibs[:20],
            "n_matched_ntc_guides": int(crp_overlap_ntc),
            "n_matched_tgt_guides": int(crp_overlap_tgt),
            "n_safe_guides_excluded_meta": int(len(safe_guides)),
        }
    except Exception as exc:
        return {"feasible": False, "reason": f"{type(exc).__name__}: {exc}"}


def build_hcasmc_parts(
    raw_dir: Path,
    hgnc: HGNCLookup,
    max_pct_mt: float = 10.0,
    min_ngenes: int = 500,
    batch_cols: int = 20000,
) -> Tuple[List[PartSpec], dict]:
    """
    HCASMC-hTERT (`IGVFDS6248NDUW`, Parse v1.2.1 CSC HDF5):
      - Strictly chunked `batch_cols` reads (zero full-matrix densification)
      - Strict `type == "non-targeting"` ONLY (excluding `safe-targeting` / `0Safe_`)
      - Poisson-Gaussian bimodal EM (`min_thr=5, max_thr=20, default_thr=7`) + `n_pos == 1 & dominance >= 3`
      - Parse sublibrary ID parity split (`3A..4F`), with locked `seed=42` 50/50 fallback if `< 2` sublibs.
    """
    hdir = raw_dir / "hcasmc"
    h5_path = _find_first_file(hdir, ["IGVFFI5999NXAF.h5"], "*.h5*")
    ntc_guides, safe_guides, guide_to_target = _load_hcasmc_guide_metadata(raw_dir, hgnc)

    with h5py.File(h5_path, "r") as f:
        m = f["matrix"]
        n_feat, n_col = int(m["shape"][0]), int(m["shape"][1])
        barcodes = decode_str_array(m["barcodes"][:])
        ftypes = np.array(decode_str_array(m["features/feature_type"][:]))
        fids = np.array(decode_str_array(m["features/id"][:]))
        fnames = np.array(decode_str_array(m["features/name"][:]))
        rna_rows = np.where(ftypes == "Gene Expression")[0]
        crp_rows = np.where(ftypes != "Gene Expression")[0]
        if len(rna_rows) == 0 or len(crp_rows) == 0:
            raise RuntimeError(f"HCASMC HDF5 {h5_path} missing RNA or guide features.")
        mt_raw = np.zeros(n_feat, dtype=bool)
        mt_raw[rna_rows] = [str(x).startswith("MT-") for x in fnames[rna_rows]]

    crp_keys = [
        fids[i] if (fids[i] in ntc_guides or fids[i] in safe_guides or fids[i] in guide_to_target) else fnames[i]
        for i in crp_rows
    ]

    rna_umi, rna_ng, rna_mt = compute_10x_csc_rna_qc(h5_path, rna_rows, mt_raw, batch_cols=batch_cols)
    ok_init = (rna_mt < 20.0) & (rna_ng >= 200)
    ok_hq = (rna_mt <= max_pct_mt) & (rna_ng >= min_ngenes) & (rna_umi > 0.0)
    hq_cols = np.where(ok_hq)[0]

    # Second compact CSC pass strictly extracting guide-capture entries for histogram & HQ matrix
    is_crp = np.zeros(n_feat, dtype=bool)
    is_crp[crp_rows] = True
    crp_remap = np.full(n_feat, -1, dtype=np.int32)
    crp_remap[crp_rows] = np.arange(len(crp_rows), dtype=np.int32)
    col_map_hq = np.full(n_col, -1, dtype=np.int32)
    col_map_hq[hq_cols] = np.arange(len(hq_cols), dtype=np.int32)

    hist = np.zeros((len(crp_rows), 66), dtype=np.int64)
    hq_row_l, hq_col_l, hq_dat_l = [], [], []
    with h5py.File(h5_path, "r") as f:
        m = f["matrix"]
        indptr = m["indptr"][:]
        for cs in range(0, n_col, batch_cols):
            ce = min(cs + batch_cols, n_col)
            p0, p1 = int(indptr[cs]), int(indptr[ce])
            if p1 == p0:
                continue
            dat = m["data"][p0:p1]
            idx = m["indices"][p0:p1]
            m_c = is_crp[idx]
            if not np.any(m_c):
                continue
            col_lens = np.diff(indptr[cs : ce + 1])
            col_ids = np.repeat(np.arange(cs, ce, dtype=np.int32), col_lens)[m_c]
            c_idx = crp_remap[idx[m_c]]
            c_dat = dat[m_c].astype(np.int64, copy=False)

            m_i = ok_init[col_ids]
            if np.any(m_i):
                np.add.at(hist, (c_idx[m_i], np.clip(c_dat[m_i], 1, 65)), 1)

            hc = col_map_hq[col_ids]
            m_h = hc >= 0
            if np.any(m_h):
                hq_row_l.append(c_idx[m_h])
                hq_col_l.append(hc[m_h])
                hq_dat_l.append(c_dat[m_h].astype(np.float64, copy=False))

    if not hq_dat_l:
        raise RuntimeError("HCASMC yielded zero guide-capture counts on QC-passing cells.")

    C_hq = sp.csr_matrix(
        (np.concatenate(hq_dat_l), (np.concatenate(hq_row_l), np.concatenate(hq_col_l))),
        shape=(len(crp_rows), len(hq_cols)),
    ).tocsc()

    thrs = fit_poisson_gaussian_thresholds_from_hist(hist, min_thr=5, max_thr=20, default_thr=7)
    t1, t2, a1 = extract_top2_csc(C_hq)
    pos_mask = C_hq.data >= thrs[C_hq.indices]
    pos_csc = sp.csc_matrix((pos_mask.astype(np.int8), C_hq.indices, C_hq.indptr), shape=C_hq.shape)
    n_pos = np.asarray(pos_csc.sum(axis=0)).ravel()
    s1_mask = (n_pos == 1) & (t1 >= 5.0) & ((t2 < 3.0) | (t1 / np.maximum(t2, 0.5) >= 3.0))

    raw_keep: List[int] = []
    nt_flags: List[bool] = []
    pkeys: List[str] = []
    pmeta: Dict[str, Tuple[str, str, str, str]] = {}
    sublibs_kept: List[Optional[str]] = []
    n_safe_excl = 0

    for loc_j in np.where(s1_mask)[0]:
        gkey = crp_keys[int(a1[loc_j])]
        if gkey in safe_guides:
            n_safe_excl += 1
            continue
        raw_col = int(hq_cols[loc_j])
        mm = _PARSE_SUBLIB_RE.search(barcodes[raw_col])
        sublib = mm.group(1) if mm else None
        if gkey in ntc_guides:
            raw_keep.append(raw_col)
            nt_flags.append(True)
            pkeys.append("")
            sublibs_kept.append(sublib)
        elif gkey in guide_to_target:
            raw_keep.append(raw_col)
            nt_flags.append(False)
            pkeys.append(gkey)
            pmeta.setdefault(gkey, guide_to_target[gkey])
            sublibs_kept.append(sublib)

    keep_arr = np.asarray(raw_keep, dtype=np.int64)
    nt_arr = np.asarray(nt_flags, dtype=bool)
    distinct_sublibs = sorted({s for s in sublibs_kept if s is not None})
    used_seed42_backup = False

    if len(distinct_sublibs) >= 2 and all(s is not None for s in sublibs_kept):
        sub_to_half = {s: (1 if i % 2 == 0 else 2) for i, s in enumerate(distinct_sublibs)}
        half_arr = np.array([sub_to_half[str(s)] for s in sublibs_kept], dtype=np.int8)
        batch_raw = np.array([f"sublib_{s}" for s in sublibs_kept], dtype=str)
        batch_arr, n_fb = _apply_min_ntc_batch_fallback(batch_raw, nt_arr, half_arr, min_batch_nt=15)
        if np.sum(nt_arr & (half_arr == 1)) == 0 or np.sum(nt_arr & (half_arr == 2)) == 0:
            used_seed42_backup = True
    else:
        used_seed42_backup = True

    if used_seed42_backup:
        rng42 = np.random.default_rng(42)
        half_arr = np.ones(len(keep_arr), dtype=np.int8)
        # Stratified deterministic 50/50 split within NTC and per perturbation target
        grp_to_idx: Dict[str, List[int]] = defaultdict(list)
        for i, (is_nt, pk) in enumerate(zip(nt_arr, pkeys)):
            grp_to_idx["__NTC__" if is_nt else pmeta[pk][1]].append(i)
        for g_key in sorted(grp_to_idx.keys()):
            idxs = np.array(grp_to_idx[g_key], dtype=np.int64)
            perm = rng42.permutation(len(idxs))
            half_arr[idxs[perm[len(idxs) // 2 :]]] = 2
        batch_arr = np.full(len(keep_arr), "pooled_all", dtype=str)
        n_fb = 0

    read_chunk, read_rows, read_chunk_cols = make_10x_csc_subset_readers(h5_path, keep_arr, n_feat, rna_rows)
    rids = [str(fids[i]) for i in rna_rows]
    rnames = [str(fnames[i]) for i in rna_rows]
    part = PartSpec(
        tag="HCASMC_hTERT",
        path=h5_path,
        n_cells=len(keep_arr),
        n_genes_raw=len(rna_rows),
        nt_mask=nt_arr,
        batches=batch_arr,
        half_flags=half_arr,
        umi=rna_umi[keep_arr],
        pert_keys=pkeys,
        pert_meta=pmeta,
        var_gene_ids=rids,
        var_gene_names=rnames,
        read_chunk_fn=read_chunk,
        read_rows_fn=read_rows,
        read_chunk_cols_fn=read_chunk_cols,
    )
    audit = {
        "n_raw_cells": n_col,
        "n_hq_rna_cells": int(len(hq_cols)),
        "n_safe_targeting_excluded": int(n_safe_excl),
        "n_qc_kept_cells": int(len(keep_arr)),
        "n_kept_ntc": int(nt_arr.sum()),
        "n_kept_half_a_ntc": int(np.sum(nt_arr & (half_arr == 1))),
        "n_kept_half_b_ntc": int(np.sum(nt_arr & (half_arr == 2))),
        "distinct_sublibs": distinct_sublibs,
        "used_stratified_seed42_backup": bool(used_seed42_backup),
        "n_low_ntc_sublib_fallback_cells": int(n_fb),
        "median_guide_threshold": float(np.median(thrs)),
    }
    return [part], audit


# ==============================================================================
# 3. CD4+ T Stimulated 48h Reader (Reusing Round 10 CD4+ T Flex Logic + Q6)
# ==============================================================================


def build_cd4t_stim48_parts(
    raw_dir: Path,
    hgnc: HGNCLookup,
    shared_meta_dir: Optional[Path] = None,
    allow_partial_lanes: bool = False,
) -> Tuple[List[PartSpec], dict]:
    """
    CD4+ T Stimulated 48hr CRISPRi (`GSE314342 Stim48hr`):
      - Half A (`1`): Donor `CE0008678` (`D3_Stim48hr`, `CD4i_R2L25..48`, 24 lanes)
      - Half B (`2`): Donor `CE0010866` (`D2_Stim48hr`, `CD4i_R2L25..48`, 24 lanes)
      - Q6 Disjointness: rejects any R10 `Rest` / `D1` / `D4` / `CE0006864` / `CE0008162` file.
      - Identical RNA QC (`MT <= 5.0%`, `ngenes >= 500`), per-donor Poisson-Gaussian `P4` guide
        filtering, `CD4_EXCEL_GUIDE_REMAP`, and strict `NTC-` control selection as Round 10.
    """
    cdir = raw_dir / "cd4t_stim48"
    csv_name = "sgrna_library_metadata.suppl_table.csv"
    if (cdir / csv_name).exists():
        sg_csv = cdir / csv_name
    elif shared_meta_dir is not None and (shared_meta_dir / csv_name).exists():
        sg_csv = shared_meta_dir / csv_name
    else:
        r10_cached = REPO_ROOT / "round10" / "shared" / "raw" / "cd4t_rest" / csv_name
        if not r10_cached.exists():
            raise FileNotFoundError(f"Missing {csv_name} in {cdir} or {r10_cached}")
        sg_csv = r10_cached

    sg_df = pd.read_csv(sg_csv)
    guide_to_target: Dict[str, Tuple[str, str, str, str]] = {}
    ntc_guide_set: Set[str] = set()
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

    # Q6 overlap check
    bad_files = [
        p.name
        for p in cdir.glob("*.h5")
        if any(tok in p.name for tok in ("Rest", "_D1_", "_D4_", "CE0006864", "CE0008162"))
    ]
    if bad_files:
        raise RuntimeError(f"Q6 disjointness violation: found Rest/D1/D4 files in cd4t_stim48: {bad_files}")

    d3_files = sorted(cdir.glob("GSM93938*_D3_Stim48hr_*.h5"))
    d2_files = sorted(cdir.glob("GSM93938*_D2_Stim48hr_*.h5"))
    if not allow_partial_lanes and (len(d3_files) != 24 or len(d2_files) != 24):
        raise RuntimeError(f"Expected 24 D3_Stim48hr + 24 D2_Stim48hr .h5 files, found D3={len(d3_files)}, D2={len(d2_files)}")
    if not d3_files or not d2_files:
        raise RuntimeError(f"Missing D3_Stim48hr ({len(d3_files)}) or D2_Stim48hr ({len(d2_files)}) files in {cdir}")

    parts: List[PartSpec] = []
    donor_audit: Dict[str, dict] = {}

    for donor_id, half_val, files in [
        ("CE0008678_D3_Stim48hr", 1, d3_files),
        ("CE0010866_D2_Stim48hr", 2, d2_files),
    ]:
        lane_cache = []
        hist = None
        crp_names_ref = None
        for lp in files:
            parts_name = lp.name.split("_")
            lane_tag = f"{parts_name[1]}_{parts_name[2]}" if len(parts_name) >= 3 else lp.stem
            with h5py.File(lp, "r") as f:
                m = f["matrix"]
                n_feat, n_col = int(m["shape"][0]), int(m["shape"][1])
                ftypes = np.array(decode_str_array(m["features/feature_type"][:]))
                fids = np.array(decode_str_array(m["features/id"][:]))
                fnames = np.array(decode_str_array(m["features/name"][:]))
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
            np.add.at(hist, (c_idx[ok_init_mask], np.clip(c_dat[ok_init_mask], 1, 65)), 1)

            hq_cols = np.where(ok_hq)[0]
            col_map_hq = np.full(n_col, -1, dtype=np.int32)
            col_map_hq[hq_cols] = np.arange(len(hq_cols), dtype=np.int32)
            hq_m = col_map_hq[c_col] >= 0
            C_hq = sp.csr_matrix(
                (c_dat[hq_m].astype(np.float64), (c_idx[hq_m], col_map_hq[c_col[hq_m]])),
                shape=(len(crp_rows), len(hq_cols)),
            ).tocsc()
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
                    pmeta.setdefault(gname, guide_to_target[gname])
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

    return parts, {
        "donors": donor_audit,
        "q6_disjoint_from_rest": True,
        "half_rule": "HalfA=CE0008678_D3_Stim48hr, HalfB=CE0010866_D2_Stim48hr",
    }


# ==============================================================================
# 4. Feng et al. 2026 TargetedScreen Streaming CSV Reader (PREREG_D §2.2.4, H-D5)
# ==============================================================================


def _parse_feng_gene_row_label(raw_label: str) -> Tuple[str, str, bool]:
    """Parse `ENSG...:SYMBOL:Gene-Expression` -> `(ens_id, symbol, is_gene_expression)`."""
    lbl = raw_label.strip().strip('"')
    if ":CRISPR-Guide-Capture" in lbl:
        return "", "", False
    parts = lbl.split(":")
    if len(parts) >= 3:
        return parts[0].split(".")[0], parts[1], (parts[-1].lower() == "gene-expression")
    if len(parts) == 2:
        if parts[0].startswith("ENSG"):
            return parts[0].split(".")[0], parts[1], True
        return "", parts[0], (parts[1].lower() == "gene-expression")
    return "", lbl, True


def build_feng_ipsc_parts(
    raw_dir: Path,
    hgnc: HGNCLookup,
    cache_dir: Optional[Path] = None,
    bucket_size: int = 40000,
    min_batch_nt: int = 15,
) -> Tuple[List[PartSpec], dict]:
    """
    Feng 2026 `TargetedScreen` (`27989294 v2`):
      - `Q6`: Strictly excludes `kolf_2` & `kolf_3` before CSV ingestion
      - `P4` & NTC-Only: Strictly excludes `"unassigned"` (never pooled w/ NTC); NTC strictly `NonTarget_*`
      - Sorted 17 non-KOLF `Cell_Line` parity split (`odd` 1-based -> Half A, `even` -> Half B)
      - Unit centering: `Cell_Line x Inlet` if `>= 15` NTCs else `Cell_Line` / half pooled NTC
      - Pure row-by-row gzip streaming of `TargetedScreen_RNA-UMI-Counts.csv.gz` into cell-bucketed
        10x CSC HDF5 cache (`< 1 GB` peak RAM regardless of CSV size).
    """
    fdir = raw_dir / "feng_ipsc"
    meta_path = _find_first_file(fdir, ["TargetedScreen_Cell-Metadata.tsv.gz", "TargetedScreen_Cell-Metadata.tsv"], "*Cell-Metadata*.tsv*")
    csv_path = _find_first_file(fdir, ["TargetedScreen_RNA-UMI-Counts.csv.gz", "TargetedScreen_RNA-UMI-Counts.csv"], "*RNA-UMI-Counts*.csv*")

    meta_df = pd.read_csv(meta_path, sep="\t", dtype=str, keep_default_na=False)
    cid_col = _pick_column(meta_df.columns, ["Cell_ID", "cell_id", "barcode", "CBC"])
    cline_col = _pick_column(meta_df.columns, ["Cell_Line", "cell_line", "line"])
    gcall_col = _pick_column(meta_df.columns, ["Guide_Call", "guide_call", "sgrna", "perturbation"])
    tgt_col = _pick_column(meta_df.columns, ["Target_Gene", "target_gene", "gene_symbol"])
    inlet_col = _pick_column(meta_df.columns, ["Inlet", "inlet", "batch", "lane"])
    if cid_col is None or cline_col is None or gcall_col is None:
        raise KeyError(f"Feng metadata missing required columns in {meta_path}: {list(meta_df.columns)}")

    n_kolf_excluded = n_unassigned_excl = n_multiguide_excl = 0
    cell_records: Dict[str, Tuple[bool, str, str, str]] = {}
    kept_lines_set: Set[str] = set()
    pmeta: Dict[str, Tuple[str, str, str, str]] = {}

    for r in meta_df.itertuples(index=False):
        cid = str(getattr(r, cid_col)).strip()
        cline = str(getattr(r, cline_col)).strip()
        gcall = str(getattr(r, gcall_col)).strip()
        if not cid or not cline:
            continue
        if cline.lower() in _FENG_KOLF_EXCLUDED or cline.lower().startswith("kolf_2") or cline.lower().startswith("kolf_3"):
            n_kolf_excluded += 1
            continue
        low_g = gcall.lower()
        if not gcall or low_g in ("unassigned", "na", "nan", "none"):
            n_unassigned_excl += 1
            continue
        if _FENG_DELIM_RE.search(gcall):
            n_multiguide_excl += 1
            continue

        inlet = str(getattr(r, inlet_col)).strip() if inlet_col else "inlet0"
        unit_name = f"{cline}__inlet{inlet}"
        kept_lines_set.add(cline)

        if gcall.startswith("NonTarget_") or low_g.startswith("nontarget"):
            cell_records[cid] = (True, "", cline, unit_name)
        else:
            raw_sym = (str(getattr(r, tgt_col)).strip() if tgt_col else "") or re.split(r"[_-]", gcall, maxsplit=1)[0]
            canon_sym, hstat = hgnc.resolve(raw_sym, "")
            if canon_sym not in hgnc.approved:
                continue
            cell_records[cid] = (False, canon_sym, cline, unit_name)
            if canon_sym not in pmeta:
                pmeta[canon_sym] = (hgnc.sym_to_ens.get(canon_sym, ""), canon_sym, raw_sym, hstat)

    lines_sorted = sorted(kept_lines_set)
    if len(lines_sorted) < 2:
        raise RuntimeError(f"Expected >= 2 non-KOLF cell lines in Feng metadata, got {lines_sorted}")
    line_to_half = {ln: (1 if i % 2 == 0 else 2) for i, ln in enumerate(lines_sorted)}

    cache_root = cache_dir if cache_dir is not None else fdir
    cache_root.mkdir(parents=True, exist_ok=True)
    cache_h5_path = cache_root / "feng_ipsc_kept_csc.h5"

    # Stream CSV row-by-row into cell buckets on temp disk
    opener = gzip.open if str(csv_path).endswith(".gz") else open
    with opener(csv_path, "rt", encoding="utf-8", errors="replace") as f_csv:
        reader = csv.reader(f_csv)
        header = [h.strip() for h in next(reader)]
        first_data_row = next(reader, None)
        if first_data_row is None:
            raise RuntimeError(f"Empty data rows in {csv_path}")

        # Support BOTH standard header (["", c1, c2, ...], diff=0) and R write.csv header ([c1, c2, ...], diff=1)
        col_offset = len(first_data_row) - len(header)
        if col_offset not in (0, 1):
            raise RuntimeError(
                f"Header/row width mismatch in {csv_path}: header={len(header)}, row1={len(first_data_row)}"
            )
        hdr_start = 1 if col_offset == 0 else 0

        keep_csv_cols_l: List[int] = []
        kept_cids: List[str] = []
        nt_list: List[bool] = []
        pkeys: List[str] = []
        unit_list: List[str] = []
        line_list: List[str] = []
        half_list: List[int] = []

        for h_idx in range(hdr_start, len(header)):
            cid = header[h_idx]
            rec = cell_records.get(cid)
            if rec is None:
                continue
            is_nt, pk, cline, uname = rec
            keep_csv_cols_l.append(h_idx + col_offset)
            kept_cids.append(cid)
            nt_list.append(is_nt)
            pkeys.append(pk)
            unit_list.append(uname)
            line_list.append(cline)
            half_list.append(line_to_half[cline])

        n_kept = len(kept_cids)
        if n_kept == 0:
            raise RuntimeError("Zero kept cells matched header of Feng RNA-UMI-Counts CSV.")

        keep_csv_cols = np.asarray(keep_csv_cols_l, dtype=np.int64)
        n_buckets = max(1, (n_kept + bucket_size - 1) // bucket_size)
        bucket_of_col = np.arange(n_kept, dtype=np.int64) // bucket_size
        local_col_in_bucket = (np.arange(n_kept, dtype=np.int64) % bucket_size).astype(np.int32)

        rna_gene_ids: List[str] = []
        rna_gene_names: List[str] = []
        n_guide_rows_skipped = 0

        with tempfile.TemporaryDirectory() as tmpd:
            tmp_path = Path(tmpd)
            g_fps = [open(tmp_path / f"b_{bi}.g.bin", "wb") for bi in range(n_buckets)]
            c_fps = [open(tmp_path / f"b_{bi}.c.bin", "wb") for bi in range(n_buckets)]
            v_fps = [open(tmp_path / f"b_{bi}.v.bin", "wb") for bi in range(n_buckets)]
            try:
                g_idx = 0

                def _consume_row(row_vals: List[str]) -> None:
                    nonlocal g_idx, n_guide_rows_skipped
                    if not row_vals:
                        return
                    eid, sym, is_gex = _parse_feng_gene_row_label(row_vals[0])
                    if not is_gex or not sym:
                        n_guide_rows_skipped += 1
                        return
                    vals = np.fromiter((float(row_vals[j]) for j in keep_csv_cols), dtype=np.float32, count=n_kept)
                    nz = np.flatnonzero(vals > 0.0)
                    if nz.size > 0:
                        b_ids = bucket_of_col[nz]
                        loc_c = local_col_in_bucket[nz]
                        nz_v = vals[nz]
                        g_arr = np.full(nz.size, g_idx, dtype=np.int32)
                        for bi in np.unique(b_ids):
                            mb = b_ids == bi
                            g_arr[mb].tofile(g_fps[int(bi)])
                            loc_c[mb].tofile(c_fps[int(bi)])
                            nz_v[mb].tofile(v_fps[int(bi)])
                    rna_gene_ids.append(eid or sym)
                    rna_gene_names.append(sym)
                    g_idx += 1

                _consume_row(first_data_row)
                for row in reader:
                    _consume_row(row)
            finally:
                for fp in g_fps + c_fps + v_fps:
                    fp.close()

            n_genes_rna = len(rna_gene_names)
            if n_genes_rna == 0:
                raise RuntimeError(f"Zero Gene-Expression rows parsed from {csv_path}.")

            # Assemble bucket-by-bucket directly into chunked 10x CSC HDF5
            umi_kept = np.zeros(n_kept, dtype=np.float64)
            with h5py.File(cache_h5_path, "w") as hf:
                mg = hf.create_group("matrix")
                mg.create_dataset("shape", data=np.array([n_genes_rna, n_kept], dtype=np.int64))
                mg.create_dataset("barcodes", data=np.array(kept_cids, dtype=h5py.string_dtype("utf-8")))
                fg = mg.create_group("features")
                fg.create_dataset("id", data=np.array(rna_gene_ids, dtype=h5py.string_dtype("utf-8")))
                fg.create_dataset("name", data=np.array(rna_gene_names, dtype=h5py.string_dtype("utf-8")))
                fg.create_dataset("feature_type", data=np.array(["Gene Expression"] * n_genes_rna, dtype=h5py.string_dtype("utf-8")))

                indptr_ds = mg.create_dataset("indptr", shape=(n_kept + 1,), dtype=np.int64)
                indptr_ds[0] = 0
                data_ds = mg.create_dataset("data", shape=(0,), maxshape=(None,), dtype=np.float32, chunks=(65536,))
                idx_ds = mg.create_dataset("indices", shape=(0,), maxshape=(None,), dtype=np.int32, chunks=(65536,))

                cur_nnz = 0
                for bi in range(n_buckets):
                    cs = bi * bucket_size
                    ce = min(cs + bucket_size, n_kept)
                    nc_b = ce - cs
                    ga = np.fromfile(tmp_path / f"b_{bi}.g.bin", dtype=np.int32)
                    ca = np.fromfile(tmp_path / f"b_{bi}.c.bin", dtype=np.int32)
                    va = np.fromfile(tmp_path / f"b_{bi}.v.bin", dtype=np.float32)
                    if ga.size > 0:
                        csc_b = sp.coo_matrix((va, (ga, ca)), shape=(n_genes_rna, nc_b)).tocsc()
                        umi_kept[cs:ce] = np.asarray(csc_b.sum(axis=0)).ravel()
                        add_nnz = int(csc_b.nnz)
                        if add_nnz > 0:
                            data_ds.resize((cur_nnz + add_nnz,))
                            idx_ds.resize((cur_nnz + add_nnz,))
                            data_ds[cur_nnz : cur_nnz + add_nnz] = csc_b.data
                            idx_ds[cur_nnz : cur_nnz + add_nnz] = csc_b.indices
                        indptr_ds[cs + 1 : ce + 1] = cur_nnz + csc_b.indptr[1:]
                        cur_nnz += add_nnz
                    else:
                        indptr_ds[cs + 1 : ce + 1] = cur_nnz

    nt_arr = np.asarray(nt_list, dtype=bool)
    half_arr = np.asarray(half_list, dtype=np.int8)
    batch_arr, n_fb = _apply_min_ntc_batch_fallback(
        np.array(unit_list, dtype=str),
        nt_arr,
        half_arr,
        group_labels=np.array(line_list, dtype=str),
        min_batch_nt=min_batch_nt,
    )

    keep_all_idx = np.arange(n_kept, dtype=np.int64)
    read_chunk, read_rows, read_chunk_cols = make_10x_csc_subset_readers(
        cache_h5_path, keep_all_idx, n_genes_rna, None
    )
    part = PartSpec(
        tag="Feng2026_TargetedScreen_NonKOLF",
        path=csv_path,
        n_cells=n_kept,
        n_genes_raw=n_genes_rna,
        nt_mask=nt_arr,
        batches=batch_arr,
        half_flags=half_arr,
        umi=umi_kept,
        pert_keys=pkeys,
        pert_meta=pmeta,
        var_gene_ids=rna_gene_ids,
        var_gene_names=rna_gene_names,
        read_chunk_fn=read_chunk,
        read_rows_fn=read_rows,
        read_chunk_cols_fn=read_chunk_cols,
    )
    audit = {
        "n_kolf_cells_excluded": int(n_kolf_excluded),
        "n_unassigned_excluded": int(n_unassigned_excl),
        "n_multiguide_excluded": int(n_multiguide_excl),
        "n_guide_rows_skipped_in_csv": int(n_guide_rows_skipped),
        "n_non_kolf_lines": len(lines_sorted),
        "non_kolf_lines_sorted": lines_sorted,
        "half_a_lines": [ln for ln in lines_sorted if line_to_half[ln] == 1],
        "half_b_lines": [ln for ln in lines_sorted if line_to_half[ln] == 2],
        "n_qc_kept_cells": int(n_kept),
        "n_kept_ntc": int(nt_arr.sum()),
        "n_kept_half_a_ntc": int(np.sum(nt_arr & (half_arr == 1))),
        "n_kept_half_b_ntc": int(np.sum(nt_arr & (half_arr == 2))),
        "n_low_ntc_unit_fallback_cells": int(n_fb),
        "cache_h5_path": str(cache_h5_path),
    }
    return [part], audit


def prepare_r11_line(
    line: str,
    raw_dir: Path,
    out_root: Path,
    hgnc: Optional[HGNCLookup] = None,
    workers: int = 8,
    chunk_size: int = 8192,
) -> Path:
    """Standardize any Round 11 candidate line cleanly using Round 10 `process_r10_line`."""
    hgnc_l = hgnc or build_hgnc_lookup()
    if line == "telohaec":
        parts, audit = build_telohaec_parts(raw_dir, hgnc_l)
    elif line == "hcasmc":
        parts, audit = build_hcasmc_parts(raw_dir, hgnc_l)
    elif line == "cd4t_stim48":
        parts, audit = build_cd4t_stim48_parts(raw_dir, hgnc_l)
    elif line == "feng_ipsc":
        parts, audit = build_feng_ipsc_parts(raw_dir, hgnc_l)
    else:
        raise ValueError(f"Unsupported Round 11 line: {line}")

    return process_r10_line(
        line_name=line,
        parts=parts,
        hgnc=hgnc_l,
        chunk_size=chunk_size,
        out_root=out_root,
        extra_meta={"adapter_audit": audit},
        n_workers=workers,
    )
