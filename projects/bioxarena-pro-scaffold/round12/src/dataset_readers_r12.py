#!/usr/bin/env python3
"""
Round 12 Dataset Ingestion Adapters (`PREREG_E` §3).

Implements bounded-memory, zero-expression-leakage readers for:
  1. `build_a549_parts`: `GSE345058` (`A549-doxCas9` CRISPRko, `CROP-seq`).
     - Reads `raw_counts.mtx.gz` + `cell_annotations.csv.gz` + `gene_annotations.csv.gz`
       + optional `sgRNA.csv.gz` (NO `.h5ad` required).
     - Streams `.mtx.gz` via fixed-size chunks into cell-bucketed CSC HDF5 cache
       (`< 1.5 GB` peak RAM), auto-detecting matrix orientation (`cells x genes` vs
       `genes x cells`) and deleting temp bucket binaries incrementally.
     - Half A = `batch1 + batch3`, Half B = `batch2 + batch4`.
     - NTC = `Nontargeting` only; single-sgRNA verified; optional `exclude_zero_cas9`
       (`total_counts_CAS9 > 0`) sensitivity toggle.
  2. `check_hek293_feasibility` & `build_hek293_parts`: `XuCao2023.h5ad` (`HEK293` CRISPRi).
     - Supports dense, CSR, and CSC `.h5ad` storage (`open_x_node`).
     - Replicate number odd/even parity split from `.obs` if present and bimodal across
       halves, else deterministic `seed=42` stratified 50/50 split within NTC & each target.
     - Single-sgRNA cells without outcome-based `>40% KD` filtering by default
       (`apply_author_kd40_filter=False`, optional sensitivity toggle).

Strictly reuses frozen `round10` / `round11` helpers cleanly without modifying any file
under `round10/` or `round11/`.
"""

from __future__ import annotations

import gzip
import re
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROUND12_SRC_DIR = Path(__file__).resolve().parent
REPO_ROOT = ROUND12_SRC_DIR.parents[1]
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
R10_SCRIPTS_DIR = REPO_ROOT / "round10" / "shared" / "scripts"
R11_SRC_DIR = REPO_ROOT / "round11" / "src"
for p in (str(R10_SRC_DIR), str(R10_SCRIPTS_DIR), str(R11_SRC_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from data_loader_r10 import iter_h5ad_csr_chunks_r10  # noqa: E402
from dataset_readers_r11 import (  # noqa: E402
    _apply_min_ntc_batch_fallback,
    _compute_h5ad_feature_subsets,
    _find_first_file,
    _is_plain_single_gene_symbol,
    _pick_column,
)
from prepare_r10_lines import (  # noqa: E402
    HGNCLookup,
    PartSpec,
    build_hgnc_lookup,
    make_10x_csc_subset_readers,
    open_x_node,
    process_r10_line,
    read_h5ad_chunk_dense,
    read_obs_or_var_col,
)

_MULTI_GUIDE_RE = re.compile(r"[,;|+]")
_A549_NTC_TOKENS = frozenset({"nontargeting", "non-targeting", "ntc", "nontarget"})
_HEK_NTC_TOKENS = frozenset(
    {"control", "ctrl", "ntc", "nt", "non-targeting", "nontargeting", "nontarget", "scramble"}
)
_UNASSIGNED_TOKENS = frozenset(
    {"", "na", "nan", "none", "null", "unassigned", "multiplet", "doublet", "ambiguous", "unknown"}
)


def _resolve_raw_subdir(raw_dir: Path, candidates: Tuple[str, ...]) -> Path:
    base = Path(raw_dir)
    for name in candidates:
        p = base / name
        if p.exists():
            return p
    return base / candidates[0]


def _extract_trailing_int(val: Any) -> Optional[int]:
    m = re.search(r"(\d+)$", str(val).strip())
    return int(m.group(1)) if m else None


def _read_mtx_header(mtx_path: Path) -> Tuple[int, int, int, int, str]:
    """Return `(n_rows, n_cols, nnz, n_skip_lines, sep)` from `.mtx` / `.mtx.gz`."""
    opener = gzip.open if str(mtx_path).endswith(".gz") else open
    n_skip = 0
    dims_line: Optional[str] = None
    first_data: Optional[str] = None
    with opener(mtx_path, "rt", encoding="utf-8", errors="replace") as f:
        for raw_line in f:
            n_skip += 1
            s = raw_line.strip()
            if not s or s.startswith("%"):
                continue
            if dims_line is None:
                dims_line = s
            else:
                first_data = s
                n_skip -= 1  # rewind count before first triplet row
                break
    if dims_line is None:
        raise RuntimeError(f"Empty or malformed MTX header in {mtx_path}")
    parts = dims_line.split()
    if len(parts) < 3:
        raise RuntimeError(f"Invalid MTX dimension line {dims_line!r} in {mtx_path}")
    n_rows, n_cols, nnz = int(parts[0]), int(parts[1]), int(parts[2])
    sep = "\t" if (first_data is not None and "\t" in first_data) else " "
    return n_rows, n_cols, nnz, n_skip, sep


# ==============================================================================
# 1. A549 (GSE345058) Streaming Matrix Market Reader
# ==============================================================================


def build_a549_parts(
    raw_dir: Path,
    hgnc: HGNCLookup,
    exclude_zero_cas9: bool = False,
    cache_dir: Optional[Path] = None,
    bucket_size: int = 40000,
    chunk_lines: int = 1_500_000,
    min_batch_nt: int = 15,
    max_pct_mt: float = 20.0,
    min_ngenes: int = 500,
) -> Tuple[List[PartSpec], dict]:
    """
    `GSE345058` (`A549-doxCas9` CRISPRko):
      - Streams `raw_counts.mtx.gz` into cell-bucketed 10x CSC HDF5 cache (`< 1.5 GB` RAM).
      - Half A (`1`) = `batch1 + batch3`; Half B (`2`) = `batch2 + batch4`.
      - Strict single-sgRNA & `Nontargeting`-only controls; plain single HGNC genes.
      - Optional `exclude_zero_cas9=True` drops cells with `total_counts_CAS9 == 0`.
    """
    adir = _resolve_raw_subdir(raw_dir, ("a549_ko", "a549"))
    cell_path = _find_first_file(adir, ["GSE345058_cell_annotations.csv.gz", "GSE345058_cell_annotations.csv"], "*cell_annotations*.csv*")
    gene_path = _find_first_file(adir, ["GSE345058_gene_annotations.csv.gz", "GSE345058_gene_annotations.csv"], "*gene_annotations*.csv*")
    mtx_path = _find_first_file(adir, ["GSE345058_raw_counts.mtx.gz", "GSE345058_raw_counts.mtx"], "*raw_counts*.mtx*")
    sgrna_candidates = sorted(adir.glob("*sgRNA*.csv*")) + sorted(adir.glob("*sgrna*.csv*"))

    cell_df = pd.read_csv(cell_path, dtype=str, keep_default_na=False)
    gene_df = pd.read_csv(gene_path, dtype=str, keep_default_na=False)

    # Known sgRNA allowlist check if sgRNA table exists
    known_ntc_sgrnas: Set[str] = set()
    if sgrna_candidates:
        sg_df = pd.read_csv(sgrna_candidates[0], dtype=str, keep_default_na=False)
        sg_id_c = _pick_column(sg_df.columns, ["sgRNA", "guide_id", "sgrna", "sequence", "id"])
        sg_pert_c = _pick_column(sg_df.columns, ["Perturbation", "perturbation", "gene", "target", "type"])
        if sg_id_c and sg_pert_c:
            for gval, pval in zip(sg_df[sg_id_c], sg_df[sg_pert_c]):
                if str(pval).strip().lower() in _A549_NTC_TOKENS:
                    known_ntc_sgrnas.add(str(gval).strip())

    c_cellid_col = _pick_column(cell_df.columns, ["cell_id", "Cell_ID", "barcode", "Barcode", "cell"])
    c_batch_col = _pick_column(cell_df.columns, ["batch", "Batch", "batch_id"])
    c_unit_col = _pick_column(cell_df.columns, ["Lane", "lane", "channel", "sample", "batch", "Batch"])
    c_sgrna_col = _pick_column(cell_df.columns, ["sgRNA", "sgrna", "guide_id", "guide"])
    c_pert_col = _pick_column(cell_df.columns, ["Perturbation", "perturbation", "target_gene", "gene_symbol", "gene"])
    c_cond_col = _pick_column(cell_df.columns, ["Condition", "condition"])
    c_cas9_col = _pick_column(cell_df.columns, ["total_counts_CAS9", "cas9_counts", "CAS9_counts", "CAS9"])
    if c_batch_col is None or c_pert_col is None or c_sgrna_col is None:
        raise KeyError(
            f"A549 cell annotations missing batch/sgRNA/Perturbation columns in {cell_path}: {list(cell_df.columns)}"
        )
    if c_cellid_col is not None:
        cids = [str(x).strip() for x in cell_df[c_cellid_col].tolist()]
        if any(not c for c in cids) or len(set(cids)) != len(cids):
            raise ValueError(
                f"A549 cell annotations contain empty or duplicate {c_cellid_col!r} entries in {cell_path}"
            )

    n_annot_cells = len(cell_df)
    n_annot_genes = len(gene_df)

    # Parse gene annotations
    g_sym_col = _pick_column(gene_df.columns, ["gene_symbol", "gene_name", "symbol", "Symbol", "name", "gene"]) or gene_df.columns[0]
    g_ens_col = _pick_column(gene_df.columns, ["ensembl_id", "gene_id", "gene_ids", "Ensembl", "ensembl"])
    g_type_col = _pick_column(gene_df.columns, ["feature_types", "feature_type", "type"])
    g_gen_col = _pick_column(gene_df.columns, ["genome", "species", "reference"])
    raw_g_syms = [str(x).strip() for x in gene_df[g_sym_col].tolist()]
    raw_g_ens = (
        [str(x).strip().split(".")[0] for x in gene_df[g_ens_col].tolist()]
        if g_ens_col is not None
        else [""] * n_annot_genes
    )
    raw_g_gen = (
        [str(x).strip().upper() for x in gene_df[g_gen_col].tolist()]
        if g_gen_col is not None
        else [""] * n_annot_genes
    )
    def _is_cas9_feat(i_g: int) -> bool:
        return (
            raw_g_syms[i_g].upper() == "CAS9"
            or raw_g_ens[i_g].upper() == "CAS9"
            or raw_g_gen[i_g] == "CAS9"
        )

    if g_type_col is not None:
        g_types = [str(x).strip() for x in gene_df[g_type_col].tolist()]
        rna_gene_raw_idx = [
            i for i, t in enumerate(g_types)
            if t in ("Gene Expression", "RNA", "") and not _is_cas9_feat(i)
        ]
        if not rna_gene_raw_idx:
            rna_gene_raw_idx = [i for i in range(n_annot_genes) if not _is_cas9_feat(i)]
    else:
        rna_gene_raw_idx = [i for i in range(n_annot_genes) if not _is_cas9_feat(i)]

    raw_to_rna = np.full(n_annot_genes, -1, dtype=np.int32)
    rna_ids: List[str] = []
    rna_names: List[str] = []
    for k, ri in enumerate(rna_gene_raw_idx):
        raw_to_rna[ri] = k
        rna_ids.append(raw_g_ens[ri] if raw_g_ens[ri] else raw_g_syms[ri])
        rna_names.append(raw_g_syms[ri])
    n_genes_rna = len(rna_names)
    mt_mask_rna = np.array([gn.startswith("MT-") for gn in rna_names], dtype=bool)

    cas9_raw_idx = next((i for i in range(n_annot_genes) if _is_cas9_feat(i)), -1)
    cas9_vals = (
        pd.to_numeric(cell_df[c_cas9_col], errors="coerce").to_numpy(dtype=np.float64)
        if c_cas9_col is not None
        else None
    )
    if exclude_zero_cas9 and cas9_vals is None and cas9_raw_idx < 0:
        raise KeyError(
            f"exclude_zero_cas9=True requested, but neither CAS9 column nor CAS9 feature found in {adir}"
        )

    batches_raw = cell_df[c_batch_col].astype(str).str.strip().tolist()
    units_raw = cell_df[c_unit_col].astype(str).str.strip().tolist() if c_unit_col else batches_raw
    sgrnas_raw = cell_df[c_sgrna_col].astype(str).str.strip().tolist()
    perts_raw = cell_df[c_pert_col].astype(str).str.strip().tolist()
    conds_raw = cell_df[c_cond_col].astype(str).str.strip().tolist() if c_cond_col is not None else None

    raw_to_kept = np.full(n_annot_cells, -1, dtype=np.int32)
    kept_nt_l: List[bool] = []
    kept_half_l: List[int] = []
    kept_unit_l: List[str] = []
    kept_pkey_l: List[str] = []
    pmeta: Dict[str, Tuple[str, str, str, str]] = {}

    n_multiguide_excl = n_cas9_zero_excl = n_bad_batch_excl = n_nonplain_excl = 0

    for ci in range(n_annot_cells):
        sg = sgrnas_raw[ci]
        pr = perts_raw[ci]
        if not sg or sg.lower() in _UNASSIGNED_TOKENS or _MULTI_GUIDE_RE.search(sg) or " " in sg:
            n_multiguide_excl += 1
            continue
        if _MULTI_GUIDE_RE.search(pr):
            n_multiguide_excl += 1
            continue

        b_num = _extract_trailing_int(batches_raw[ci])
        if b_num not in (1, 2, 3, 4):
            n_bad_batch_excl += 1
            continue
        # PREREG_E §3: Half A = batch1 + batch3; Half B = batch2 + batch4
        half_flag = 1 if b_num in (1, 3) else 2

        if exclude_zero_cas9 and cas9_vals is not None:
            cv = cas9_vals[ci]
            if not (np.isfinite(cv) and cv > 0.0):
                n_cas9_zero_excl += 1
                continue

        pert_is_ntc = pr.lower() in _A549_NTC_TOKENS
        if conds_raw is not None:
            cond_low = conds_raw[ci].lower()
            if cond_low in _UNASSIGNED_TOKENS:
                n_multiguide_excl += 1
                continue
            cond_is_ntc = cond_low in _A549_NTC_TOKENS
            if pert_is_ntc != cond_is_ntc:
                n_nonplain_excl += 1
                continue
            is_ntc = pert_is_ntc and cond_is_ntc
        else:
            is_ntc = pert_is_ntc
        if is_ntc and known_ntc_sgrnas and sg in known_ntc_sgrnas:
            pass  # verified NTC guide identity
        if is_ntc:
            k_idx = len(kept_nt_l)
            raw_to_kept[ci] = k_idx
            kept_nt_l.append(True)
            kept_half_l.append(half_flag)
            kept_unit_l.append(units_raw[ci])
            kept_pkey_l.append("")
        else:
            if not _is_plain_single_gene_symbol(pr):
                n_nonplain_excl += 1
                continue
            canon_sym, hstat = hgnc.resolve(pr, "")
            if canon_sym not in hgnc.approved:
                n_nonplain_excl += 1
                continue
            k_idx = len(kept_nt_l)
            raw_to_kept[ci] = k_idx
            kept_nt_l.append(False)
            kept_half_l.append(half_flag)
            kept_unit_l.append(units_raw[ci])
            kept_pkey_l.append(canon_sym)
            if canon_sym not in pmeta:
                pmeta[canon_sym] = (hgnc.sym_to_ens.get(canon_sym, ""), canon_sym, pr, hstat)

    n_pre_qc = len(kept_nt_l)
    if n_pre_qc == 0:
        raise RuntimeError(f"A549 metadata yielded zero valid cells in {adir}")

    n_mtx_rows, n_mtx_cols, _, n_skip_lines, sep = _read_mtx_header(mtx_path)
    if n_mtx_rows == n_annot_cells and n_mtx_cols == n_annot_genes:
        cell_first = True
    elif n_mtx_rows == n_annot_genes and n_mtx_cols == n_annot_cells:
        cell_first = False
    else:
        raise RuntimeError(
            f"A549 MTX shape ({n_mtx_rows}, {n_mtx_cols}) incompatible with "
            f"(cells={n_annot_cells}, genes={n_annot_genes}) in {mtx_path}"
        )

    cache_root = cache_dir if cache_dir is not None else adir
    cache_root.mkdir(parents=True, exist_ok=True)
    tag_name = "a549_ko_cas9pos" if exclude_zero_cas9 else "a549_ko"
    cache_h5_path = cache_root / f"{tag_name}_kept_csc.h5"

    n_buckets = (n_pre_qc + bucket_size - 1) // bucket_size
    cas9_from_mtx = np.zeros(n_pre_qc, dtype=np.float32) if (exclude_zero_cas9 and cas9_vals is None) else None

    with tempfile.TemporaryDirectory() as tmp_bin_dir:
        tmp_path = Path(tmp_bin_dir)
        f_g = [open(tmp_path / f"b_{bi}.g.bin", "wb") for bi in range(n_buckets)]
        f_c = [open(tmp_path / f"b_{bi}.c.bin", "wb") for bi in range(n_buckets)]
        f_v = [open(tmp_path / f"b_{bi}.v.bin", "wb") for bi in range(n_buckets)]

        try:
            reader = pd.read_csv(
                mtx_path,
                sep=sep,
                header=None,
                skiprows=n_skip_lines,
                names=["r", "c", "v"],
                dtype={"r": np.int32, "c": np.int32, "v": np.float32},
                chunksize=chunk_lines,
                engine="c",
            )
            for chunk in reader:
                rr = chunk["r"].to_numpy(dtype=np.int64) - 1
                cc = chunk["c"].to_numpy(dtype=np.int64) - 1
                vv = chunk["v"].to_numpy(dtype=np.float32)
                c_raw = rr if cell_first else cc
                g_raw = cc if cell_first else rr

                kept_ci = raw_to_kept[c_raw]
                if cas9_from_mtx is not None and cas9_raw_idx >= 0:
                    m_c9 = (kept_ci >= 0) & (g_raw == cas9_raw_idx) & (vv > 0)
                    if np.any(m_c9):
                        np.add.at(cas9_from_mtx, kept_ci[m_c9], vv[m_c9])

                g_rna = raw_to_rna[g_raw]
                m_ok = (kept_ci >= 0) & (g_rna >= 0) & (vv != 0.0)
                if not np.any(m_ok):
                    continue
                k_arr = kept_ci[m_ok].astype(np.int32)
                gr_arr = g_rna[m_ok].astype(np.int32)
                vr_arr = vv[m_ok].astype(np.float32)
                b_ids = k_arr // bucket_size
                for bi in np.unique(b_ids):
                    mb = b_ids == bi
                    gr_arr[mb].tofile(f_g[int(bi)])
                    (k_arr[mb] - int(bi) * bucket_size).astype(np.int32).tofile(f_c[int(bi)])
                    vr_arr[mb].tofile(f_v[int(bi)])
        finally:
            for fh in f_g + f_c + f_v:
                fh.close()

        umi_pre = np.zeros(n_pre_qc, dtype=np.float64)
        ng_pre = np.zeros(n_pre_qc, dtype=np.int32)
        mt_umi_pre = np.zeros(n_pre_qc, dtype=np.float64)

        with h5py.File(cache_h5_path, "w") as hf:
            mg = hf.create_group("matrix")
            mg.create_dataset("shape", data=np.array([n_genes_rna, n_pre_qc], dtype=np.int64))
            data_ds = mg.create_dataset("data", shape=(0,), maxshape=(None,), dtype=np.float32, chunks=(65536,))
            idx_ds = mg.create_dataset("indices", shape=(0,), maxshape=(None,), dtype=np.int32, chunks=(65536,))
            indptr_ds = mg.create_dataset("indptr", shape=(n_pre_qc + 1,), dtype=np.int64)
            indptr_ds[0] = 0
            cur_nnz = 0
            for bi in range(n_buckets):
                cs = bi * bucket_size
                ce = min(cs + bucket_size, n_pre_qc)
                nc_b = ce - cs
                gp = tmp_path / f"b_{bi}.g.bin"
                cp = tmp_path / f"b_{bi}.c.bin"
                vp = tmp_path / f"b_{bi}.v.bin"
                ga = np.fromfile(gp, dtype=np.int32)
                ca = np.fromfile(cp, dtype=np.int32)
                va = np.fromfile(vp, dtype=np.float32)
                gp.unlink(missing_ok=True)
                cp.unlink(missing_ok=True)
                vp.unlink(missing_ok=True)
                if ga.size > 0:
                    csc_b = sp.coo_matrix((va, (ga, ca)), shape=(n_genes_rna, nc_b)).tocsc()
                    umi_pre[cs:ce] = np.asarray(csc_b.sum(axis=0)).ravel()
                    ng_pre[cs:ce] = np.diff(csc_b.indptr).astype(np.int32)
                    if np.any(mt_mask_rna):
                        mt_umi_pre[cs:ce] = np.asarray(csc_b[mt_mask_rna, :].sum(axis=0)).ravel()
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

    pct_mt_pre = 100.0 * mt_umi_pre / np.maximum(umi_pre, 1.0)
    qc_ok = (umi_pre > 0.0) & (ng_pre >= min_ngenes) & (pct_mt_pre <= max_pct_mt)
    if cas9_from_mtx is not None:
        n_cas9_zero_excl += int(np.sum(qc_ok & (cas9_from_mtx <= 0.0)))
        qc_ok &= cas9_from_mtx > 0.0
    keep_idx = np.where(qc_ok)[0].astype(np.int64)
    if keep_idx.size == 0:
        raise RuntimeError(f"A549 yielded zero cells after RNA QC in {adir}")

    nt_arr = np.asarray(kept_nt_l, dtype=bool)[keep_idx]
    half_arr = np.asarray(kept_half_l, dtype=np.int8)[keep_idx]
    unit_arr = np.asarray(kept_unit_l, dtype=str)[keep_idx]
    batch_arr, n_fb = _apply_min_ntc_batch_fallback(unit_arr, nt_arr, half_arr, min_batch_nt=min_batch_nt)
    pkeys_kept = [kept_pkey_l[int(i)] for i in keep_idx]

    read_chunk, read_rows, read_chunk_cols = make_10x_csc_subset_readers(
        cache_h5_path, keep_idx, n_genes_rna, None
    )
    part = PartSpec(
        tag="A549_GSE345058_Cas9Pos" if exclude_zero_cas9 else "A549_GSE345058",
        path=mtx_path,
        n_cells=len(keep_idx),
        n_genes_raw=n_genes_rna,
        nt_mask=nt_arr,
        batches=batch_arr,
        half_flags=half_arr,
        umi=umi_pre[keep_idx],
        pert_keys=pkeys_kept,
        pert_meta=pmeta,
        var_gene_ids=rna_ids,
        var_gene_names=rna_names,
        read_chunk_fn=read_chunk,
        read_rows_fn=read_rows,
        read_chunk_cols_fn=read_chunk_cols,
    )
    audit = {
        "n_raw_annot_cells": int(n_annot_cells),
        "exclude_zero_cas9": bool(exclude_zero_cas9),
        "n_multiguide_or_unassigned_excluded": int(n_multiguide_excl),
        "n_bad_batch_excluded": int(n_bad_batch_excl),
        "n_cas9_zero_excluded": int(n_cas9_zero_excl),
        "n_nonplain_or_unapproved_excluded": int(n_nonplain_excl),
        "n_pre_qc_kept": int(n_pre_qc),
        "n_qc_kept_cells": int(len(keep_idx)),
        "n_kept_ntc": int(nt_arr.sum()),
        "n_kept_half_a_ntc": int(np.sum(nt_arr & (half_arr == 1))),
        "n_kept_half_b_ntc": int(np.sum(nt_arr & (half_arr == 2))),
        "n_low_ntc_unit_fallback_cells": int(n_fb),
        "half_rule": "batch1_plus_batch3_vs_batch2_plus_batch4",
        "cache_h5_path": str(cache_h5_path),
    }
    return [part], audit


# ==============================================================================
# 2. HEK293 (`XuCao2023.h5ad`) Reader & Pre-Eval Feasibility Check
# ==============================================================================


def _sample_node_nonneg_int(node: Any, max_elems: int = 1000) -> bool:
    """Check whether sparse group or dataset sample contains non-negative integer counts."""
    if isinstance(node, h5py.Group) and "data" in node:
        ds = node["data"]
        n = min(int(ds.shape[0]), max_elems)
        arr = np.asarray(ds[:n], dtype=np.float64)
    elif isinstance(node, h5py.Dataset):
        nr = min(int(node.shape[0]), 32)
        nc = min(int(node.shape[1]), max_elems // max(nr, 1)) if len(node.shape) >= 2 else min(int(node.shape[0]), max_elems)
        arr = np.asarray(node[:nr, :nc] if len(node.shape) >= 2 else node[:nc], dtype=np.float64).ravel()
    else:
        return False
    if arr.size == 0:
        return False
    return bool(np.all(np.isfinite(arr) & (arr >= 0.0) & (np.abs(arr - np.round(arr)) <= 1e-6)))


def check_hek293_feasibility(raw_dir: Path, hgnc: Optional[HGNCLookup] = None) -> Dict[str, Any]:
    """
    `PREREG_E §3` metadata-only feasibility gate for `XuCao2023.h5ad` (`counts/structure ONLY`,
    never computes any expression correlation or perturbation response delta).
    """
    try:
        hgnc_l = hgnc or build_hgnc_lookup()
        hdir = _resolve_raw_subdir(raw_dir, ("hek293_xucao", "hek293"))
        h5ad_path = _find_first_file(hdir, ["XuCao2023.h5ad"], "*.h5ad")
        with h5py.File(h5ad_path, "r") as f:
            obs = f["obs"]
            obs_keys = sorted(obs.keys())
            var = f["var"]
            var_keys = sorted(var.keys())
            layer_keys = sorted(f["layers"].keys()) if "layers" in f else []
            int_count_layers: List[str] = []
            if "X" in f and _sample_node_nonneg_int(f["X"]):
                int_count_layers.append("X")
            for lk in layer_keys:
                if _sample_node_nonneg_int(f["layers"][lk]):
                    int_count_layers.append(f"layers/{lk}")
            idx_col = "_index" if "_index" in obs_keys else "index"
            n_cells = len(read_obs_or_var_col(obs, idx_col))
            _, _, r_ids, r_names, _ = _compute_h5ad_feature_subsets(f)

            rep_col = next(
                (c for c in ("replicate", "rep", "sample", "batch", "technical_replicate", "tech_rep", "replicate_id") if c in obs_keys),
                None,
            )
            rep_vals = read_obs_or_var_col(obs, rep_col) if rep_col else []
            rep_labels = sorted({str(x).strip() for x in rep_vals if str(x).strip()})

            pert_col = next((c for c in ("perturbation", "target_gene", "target", "gene", "gene_symbol") if c in obs_keys), None)
            guide_col = next((c for c in ("sgRNA", "sgrna", "guide_id", "top_sgrna", "guide", "perturbation_2") if c in obs_keys), None)
            nperts_col = next((c for c in ("nperts", "n_guides", "num_guides", "n_perts") if c in obs_keys), None)
            n_ntc = n_tgt = 0
            canon_tgts: Set[str] = set()
            if pert_col:
                pvals = read_obs_or_var_col(obs, pert_col)
                gvals = read_obs_or_var_col(obs, guide_col) if guide_col else None
                npvals = (
                    pd.to_numeric(pd.Series(read_obs_or_var_col(obs, nperts_col)), errors="coerce").to_numpy()
                    if nperts_col
                    else None
                )
                for idx_c, p in enumerate(pvals):
                    ps = str(p).strip()
                    gs = str(gvals[idx_c]).strip() if gvals is not None else ""
                    if npvals is not None and np.isfinite(npvals[idx_c]) and int(npvals[idx_c]) > 1:
                        continue
                    if _MULTI_GUIDE_RE.search(ps) or (gs and _MULTI_GUIDE_RE.search(gs)):
                        continue
                    if ps.lower() in _HEK_NTC_TOKENS:
                        n_ntc += 1
                    elif _is_plain_single_gene_symbol(ps):
                        cs, _ = hgnc_l.resolve(ps, "")
                        if cs in hgnc_l.approved:
                            n_tgt += 1
                            canon_tgts.add(cs)

        has_guide_info = bool(guide_col is not None or nperts_col is not None or pert_col is not None)
        return {
            "feasible": bool(
                n_ntc > 0
                and len(canon_tgts) > 0
                and len(r_ids) > 0
                and len(int_count_layers) > 0
                and has_guide_info
            ),
            "h5ad_path": str(h5ad_path),
            "n_cells": int(n_cells),
            "n_rna_features": int(len(r_names)),
            "layers": layer_keys,
            "int_count_sources": int_count_layers,
            "obs_columns": obs_keys,
            "var_columns": var_keys,
            "detected_rep_col": rep_col,
            "detected_rep_labels": rep_labels[:20],
            "detected_pert_col": pert_col,
            "detected_guide_col": guide_col,
            "detected_nperts_col": nperts_col,
            "n_ntc_meta": int(n_ntc),
            "n_tgt_cells_meta": int(n_tgt),
            "n_unique_canon_targets_meta": int(len(canon_tgts)),
        }
    except Exception as exc:
        return {"feasible": False, "reason": f"{type(exc).__name__}: {exc}"}


def _is_nonhuman_or_vector_var(sym: str, gid: str) -> bool:
    su = sym.upper()
    gu = gid.upper()
    if gu.startswith("ENSMUSG") or su in ("CAS9", "DCAS9", "DCAS9-KRAB-MECP2", "KRAB", "MECP2_LENTIVIRUS"):
        return True
    if "DCAS9" in su or "KRAB-MECP2" in su:
        return True
    return False


def _resolve_hek293_integer_x_layer(f: h5py.File, x_layer: Optional[str]) -> str:
    """
    `PREREG_E §3.2 (E3)`: Select raw non-negative integer UMI matrix layer.
    When `x_layer` is `None` or `'auto'`, prefer `'X'` if integer, else inspect
    `layers/*` (`counts`, `raw_counts`, `raw`, etc.), or raise `EXCLUDED_FEASIBILITY_Q3`.
    """
    req = "auto" if x_layer is None else str(x_layer).strip()
    if req.lower() != "auto":
        node = f["X"] if req == "X" else f["layers"][req]
        if not _sample_node_nonneg_int(node):
            raise RuntimeError(
                f"EXCLUDED_FEASIBILITY_Q3: HEK293 layer {req!r} does not contain raw non-negative integer UMI counts."
            )
        return req

    if "X" in f and _sample_node_nonneg_int(f["X"]):
        return "X"
    if "layers" in f:
        avail = sorted(f["layers"].keys())
        preferred = [k for k in ("counts", "raw_counts", "raw", "UMI", "umi") if k in avail]
        ordered = preferred + [k for k in avail if k not in preferred]
        for lk in ordered:
            if _sample_node_nonneg_int(f["layers"][lk]):
                return lk
    raise RuntimeError(
        "EXCLUDED_FEASIBILITY_Q3: HEK293 XuCao2023.h5ad lacks raw non-negative integer UMI counts in X or layers/*."
    )


def build_hek293_parts(
    raw_dir: Path,
    hgnc: HGNCLookup,
    apply_author_kd40_filter: bool = False,
    max_pct_mt: float = 20.0,
    min_ngenes: int = 200,
    min_umi: int = 300,
    min_batch_nt: int = 15,
    force_seed42_split: bool = False,
    x_layer: str = "auto",
) -> Tuple[List[PartSpec], dict]:
    """
    `XuCao2023.h5ad` (`HEK293` PerturbSci-Kinetics CRISPRi):
      - `E3`: verifies raw non-negative integer UMI counts (`X` or `layers/*`); auto-selects
        integer layer when `X` is float-normalized, else raises `EXCLUDED_FEASIBILITY_Q3`.
      - All single-sgRNA cells by default (`apply_author_kd40_filter=False` avoids outcome selection).
      - Excludes synthetic dCas9/lentiviral & non-human spike-in features if present.
      - Stage 1: lexicographic 1-based parity on `.obs` replicate column if `>= 2` distinct
        labels populate both halves; Stage 2 fallback: deterministic `seed=42` stratified
        50/50 split on barcode-sorted indices within NTC and per target gene.
    """
    hdir = _resolve_raw_subdir(raw_dir, ("hek293_xucao", "hek293"))
    h5ad_path = _find_first_file(hdir, ["XuCao2023.h5ad"], "*.h5ad")

    with h5py.File(h5ad_path, "r") as f:
        resolved_layer = _resolve_hek293_integer_x_layer(f, x_layer)
        x_node, shape = open_x_node(f, x_layer=resolved_layer)
        n_raw_cells, n_raw_var = int(shape[0]), int(shape[1])
        obs = f["obs"]
        obs_keys = set(obs.keys())
        obs_idx_col = "_index" if "_index" in obs_keys else ("index" if "index" in obs_keys else None)
        obs_barcodes = read_obs_or_var_col(obs, obs_idx_col) if obs_idx_col else [f"cell_{i}" for i in range(n_raw_cells)]

        raw_rna_rows, _, raw_rids, raw_rnames, _ = _compute_h5ad_feature_subsets(f)
        human_mask = [not _is_nonhuman_or_vector_var(nm, gid) for nm, gid in zip(raw_rnames, raw_rids)]
        rna_rows = raw_rna_rows[np.asarray(human_mask, dtype=bool)]
        rna_ids = [gid for gid, ok in zip(raw_rids, human_mask) if ok]
        rna_names = [nm for nm, ok in zip(raw_rnames, human_mask) if ok]
        mt_mask_rna = np.array([gn.startswith("MT-") for gn in rna_names], dtype=bool)
        full_rna_subset = (len(rna_rows) == n_raw_var and np.array_equal(rna_rows, np.arange(n_raw_var)))

        pert_col = next((c for c in ("perturbation", "target_gene", "target", "gene_symbol", "gene") if c in obs_keys), None)
        if pert_col is None:
            raise KeyError(f"HEK293 .h5ad missing perturbation column in {obs_keys}")
        pert_vals = read_obs_or_var_col(obs, pert_col)

        guide_col = next((c for c in ("sgRNA", "sgrna", "guide_id", "top_sgrna", "guide", "perturbation_2") if c in obs_keys), None)
        guide_vals = read_obs_or_var_col(obs, guide_col) if guide_col else None
        nperts_col = next((c for c in ("nperts", "n_guides", "num_guides", "n_perts") if c in obs_keys), None)
        nperts_vals = pd.to_numeric(pd.Series(read_obs_or_var_col(obs, nperts_col)), errors="coerce").to_numpy() if nperts_col else None

        rep_col = next(
            (c for c in ("replicate", "rep", "sample", "batch", "technical_replicate", "tech_rep", "replicate_id") if c in obs_keys),
            None,
        )
        rep_vals = read_obs_or_var_col(obs, rep_col) if rep_col else None
        batch_col = next((c for c in ("batch", "sample", "lane", "channel", "gem_group") if c in obs_keys), None)
        batch_vals = read_obs_or_var_col(obs, batch_col) if batch_col else rep_vals

        kd_vals = None
        if apply_author_kd40_filter:
            kd_col = next((c for c in ("kd_pass", "pass_kd40", "kd_pct", "knockdown_pct", "knockdown_fraction", "kd_efficiency") if c in obs_keys), None)
            if kd_col is None:
                raise KeyError(f"apply_author_kd40_filter=True requested, but no KD column found in {obs_keys}")
            raw_kd = read_obs_or_var_col(obs, kd_col)
            kd_vals = []
            for v in raw_kd:
                vs = str(v).strip().lower()
                if vs in ("true", "yes", "1", "pass"):
                    kd_vals.append(1.0)
                elif vs in ("false", "no", "0", "fail"):
                    kd_vals.append(0.0)
                else:
                    try:
                        fv = float(vs)
                        kd_vals.append((fv / 100.0) if fv > 1.0 else fv)
                    except Exception:
                        kd_vals.append(float("nan"))

        # Streaming RNA QC computation (CSR Group, CSC-via-open_x_node CSR sparse, or dense Dataset)
        umi_all = np.zeros(n_raw_cells, dtype=np.float64)
        ng_all = np.zeros(n_raw_cells, dtype=np.int32)
        mt_umi_all = np.zeros(n_raw_cells, dtype=np.float64)
        if isinstance(x_node, h5py.Group) and "indptr" in x_node:
            layer_arg = None if resolved_layer == "X" else resolved_layer
            for cs, ce, blk in iter_h5ad_csr_chunks_r10(h5ad_path, row_chunk=4096, layer_name=layer_arg):
                sub = blk if full_rna_subset else blk[:, rna_rows]
                umi_all[cs:ce] = sub.sum(axis=1)
                ng_all[cs:ce] = np.sum(sub > 0, axis=1)
                if np.any(mt_mask_rna):
                    mt_umi_all[cs:ce] = sub[:, mt_mask_rna].sum(axis=1)
        elif sp.issparse(x_node):
            for cs in range(0, n_raw_cells, 4096):
                ce = min(cs + 4096, n_raw_cells)
                sub = x_node[cs:ce, :] if full_rna_subset else x_node[cs:ce, :][:, rna_rows]
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
        batches_l: List[str] = []
        rep_strs_l: List[str] = []
        kept_bcs_l: List[str] = []
        pkeys: List[str] = []
        pmeta: Dict[str, Tuple[str, str, str, str]] = {}
        n_multiguide_excl = n_kd40_excl = n_nonplain_excl = 0

        for r_i in range(n_raw_cells):
            if not (umi_all[r_i] >= float(min_umi) and ng_all[r_i] >= min_ngenes and pct_mt_all[r_i] <= max_pct_mt):
                continue
            ps = str(pert_vals[r_i]).strip()
            plow = ps.lower()
            gs = str(guide_vals[r_i]).strip() if guide_vals is not None else ""
            is_ntc = plow in _HEK_NTC_TOKENS or plow.startswith("nontarget") or plow.startswith("non-targeting")

            if nperts_vals is not None and np.isfinite(nperts_vals[r_i]) and int(nperts_vals[r_i]) > 1:
                n_multiguide_excl += 1
                continue
            if _MULTI_GUIDE_RE.search(ps) or (gs and _MULTI_GUIDE_RE.search(gs)):
                n_multiguide_excl += 1
                continue
            if not is_ntc and plow in _UNASSIGNED_TOKENS:
                n_multiguide_excl += 1
                continue

            b_lbl = str(batch_vals[r_i]).strip() if batch_vals is not None else "batch_0"
            r_str = str(rep_vals[r_i]).strip() if rep_vals is not None else ""

            if is_ntc:
                keep_rows.append(r_i)
                nt_flags.append(True)
                batches_l.append(b_lbl)
                rep_strs_l.append(r_str)
                kept_bcs_l.append(str(obs_barcodes[r_i]))
                pkeys.append("")
            else:
                if kd_vals is not None:
                    kv = kd_vals[r_i]
                    if not (np.isfinite(kv) and kv > 0.40):
                        n_kd40_excl += 1
                        continue
                if not _is_plain_single_gene_symbol(ps):
                    n_nonplain_excl += 1
                    continue
                canon_sym, hstat = hgnc.resolve(ps, "")
                if canon_sym not in hgnc.approved:
                    n_nonplain_excl += 1
                    continue
                keep_rows.append(r_i)
                nt_flags.append(False)
                batches_l.append(b_lbl)
                rep_strs_l.append(r_str)
                kept_bcs_l.append(str(obs_barcodes[r_i]))
                pkeys.append(canon_sym)
                if canon_sym not in pmeta:
                    pmeta[canon_sym] = (hgnc.sym_to_ens.get(canon_sym, ""), canon_sym, ps, hstat)

        keep_arr = np.asarray(keep_rows, dtype=np.int64)
        nt_arr = np.asarray(nt_flags, dtype=bool)
        if len(keep_arr) == 0:
            raise RuntimeError(f"HEK293 {h5ad_path} yielded zero valid single-sgRNA cells.")

        # Stage 1: lexicographic 1-based parity on obs replicate column if >= 2 distinct labels
        distinct_reps = sorted({x for x in rep_strs_l if x and x.lower() not in _UNASSIGNED_TOKENS})
        used_seed42_fallback = False
        half_rule_used = f"obs_replicate_parity:{rep_col}"
        if (not force_seed42_split) and rep_col is not None and len(distinct_reps) >= 2 and all(bool(x) for x in rep_strs_l):
            rep_to_half = {lbl: (1 if i % 2 == 0 else 2) for i, lbl in enumerate(distinct_reps)}
            half_arr = np.array([rep_to_half[x] for x in rep_strs_l], dtype=np.int8)
            batch_arr, n_fb = _apply_min_ntc_batch_fallback(
                np.array(batches_l, dtype=str), nt_arr, half_arr, min_batch_nt=min_batch_nt
            )
            if np.sum(nt_arr & (half_arr == 1)) == 0 or np.sum(nt_arr & (half_arr == 2)) == 0:
                used_seed42_fallback = True
        else:
            used_seed42_fallback = True

        # Stage 2 fallback: stratified 50/50 split on barcode-sorted indices (seed=42)
        if used_seed42_fallback:
            half_rule_used = "seed42_stratified_50_50"
            rng42 = np.random.default_rng(42)
            half_arr = np.ones(len(keep_arr), dtype=np.int8)
            grp_to_idx: Dict[str, List[int]] = defaultdict(list)
            for i, (is_nt, pk) in enumerate(zip(nt_arr, pkeys)):
                grp_to_idx["__NTC__" if is_nt else pk].append(i)
            for g_key in sorted(grp_to_idx.keys()):
                idxs_sorted = sorted(grp_to_idx[g_key], key=lambda idx: (kept_bcs_l[idx], idx))
                idxs = np.array(idxs_sorted, dtype=np.int64)
                perm = rng42.permutation(len(idxs))
                half_arr[idxs[perm[len(idxs) // 2 :]]] = 2
            batch_arr, n_fb = _apply_min_ntc_batch_fallback(
                np.array(batches_l, dtype=str), nt_arr, half_arr, min_batch_nt=min_batch_nt
            )

        if isinstance(x_node, h5py.Group):
            x_full = sp.csr_matrix(
                (x_node["data"][:], x_node["indices"][:], x_node["indptr"][:]),
                shape=(n_raw_cells, n_raw_var),
            )
            x_kept = x_full[keep_arr][:, rna_rows].tocsr() if not full_rna_subset else x_full[keep_arr].tocsr()
        elif sp.issparse(x_node):
            x_kept = x_node[keep_arr][:, rna_rows].tocsr() if not full_rna_subset else x_node[keep_arr].tocsr()
        else:
            x_kept = x_node[keep_arr][:, rna_rows] if not full_rna_subset else x_node[keep_arr]

    n_genes_rna = len(rna_rows)
    part = PartSpec(
        tag="HEK293_XuCao2023",
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
        read_chunk_fn=lambda s, e, mat=x_kept: read_h5ad_chunk_dense(mat, s, e, n_genes_rna),
    )
    audit = {
        "n_raw_cells": int(n_raw_cells),
        "x_layer_used": resolved_layer,
        "apply_author_kd40_filter": bool(apply_author_kd40_filter),
        "n_multiguide_or_unassigned_excluded": int(n_multiguide_excl),
        "n_kd40_excluded": int(n_kd40_excl),
        "n_nonplain_or_unapproved_excluded": int(n_nonplain_excl),
        "n_qc_kept_cells": int(len(keep_arr)),
        "n_kept_ntc": int(nt_arr.sum()),
        "n_kept_half_a_ntc": int(np.sum(nt_arr & (half_arr == 1))),
        "n_kept_half_b_ntc": int(np.sum(nt_arr & (half_arr == 2))),
        "n_low_ntc_batch_fallback_cells": int(n_fb),
        "used_seed42_fallback": bool(used_seed42_fallback),
        "half_rule_used": half_rule_used,
    }
    return [part], audit


def prepare_r12_line(
    line: str,
    raw_dir: Path,
    out_root: Path,
    hgnc: Optional[HGNCLookup] = None,
    workers: int = 8,
    chunk_size: int = 8192,
    exclude_zero_cas9: bool = False,
    apply_author_kd40_filter: bool = False,
) -> Path:
    """Standardize any Round 12 candidate line via `process_r10_line`."""
    hgnc_l = hgnc or build_hgnc_lookup()
    if line in ("a549_ko", "a549_ko_cas9pos", "a549", "a549_cas9pos"):
        use_c9 = exclude_zero_cas9 or (line in ("a549_ko_cas9pos", "a549_cas9pos"))
        parts, audit = build_a549_parts(
            raw_dir,
            hgnc_l,
            exclude_zero_cas9=use_c9,
            cache_dir=Path(out_root) / f".cache_{line}",
        )
    elif line in ("hek293_xucao", "hek293"):
        parts, audit = build_hek293_parts(raw_dir, hgnc_l, apply_author_kd40_filter=apply_author_kd40_filter)
    else:
        raise ValueError(f"Unsupported Round 12 line: {line}")

    return process_r10_line(
        line_name=line,
        parts=parts,
        hgnc=hgnc_l,
        chunk_size=chunk_size,
        out_root=out_root,
        extra_meta={"adapter_audit": audit},
        n_workers=workers,
    )
