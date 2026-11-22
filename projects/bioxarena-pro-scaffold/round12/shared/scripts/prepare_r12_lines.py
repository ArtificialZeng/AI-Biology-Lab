#!/usr/bin/env python3
"""
Round 12 Standardized Line Cache Builder (`PREREG_E`).

Produces standardized caches under `round12/shared/lines/<line>/` strictly
matching frozen Round 10/11 schemas (`process_r10_line`):
  - `a549_ko`           (GSE345058, all QC-passed cells; confirmatory primary)
  - `a549_ko_cas9pos`   (GSE345058, CAS9 > 0 sensitivity cohort per PREREG_E §5.5)
  - `hek293_xucao`      (XuCao2023.h5ad; ONLY if `E3` feasibility gate passes)

Always invoke via `python3 -I`. Cleans temporary `.cache_*` / `*_kept_csc.h5`
intermediates so `round12/shared/raw/` and `round12/shared/lines/` remain clean.

Note on `a549_ko` / `a549_ko_cas9pos` int64 HDF5 indptr fix:
  `round12/src/dataset_readers_r12.py:424` computes `cur_nnz + csc_b.indptr[1:]`
  where `csc_b.indptr` is `np.int32`. On `GSE345058` (`nnz = 3,598,670,133 > 2^31-1`),
  NumPy 2.x raises `OverflowError: Python integer ... out of bounds for int32` once
  `cur_nnz` crosses `2^31-1`. `build_a549_parts_int64` runs the exact same logic
  as `dataset_readers_r12.build_a549_parts` with `.astype(np.int64)` on `csc_b.indptr[1:]`.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
import tempfile
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

SCRIPTS_DIR = Path(__file__).resolve().parent
ROUND12_ROOT = SCRIPTS_DIR.parents[1]
REPO_ROOT = ROUND12_ROOT.parent
R10_SRC = REPO_ROOT / "round10" / "src"
R10_SCRIPTS = REPO_ROOT / "round10" / "shared" / "scripts"
R11_SRC = REPO_ROOT / "round11" / "src"
R12_SRC = ROUND12_ROOT / "src"
for p in (str(R10_SRC), str(R10_SCRIPTS), str(R11_SRC), str(R12_SRC)):
    if p not in sys.path:
        sys.path.insert(0, p)

from prepare_r10_lines import (  # noqa: E402
    DEFAULT_HGNC_PATH,
    HGNCLookup,
    PartSpec,
    build_hgnc_lookup,
    make_10x_csc_subset_readers,
    process_r10_line,
)
from dataset_readers_r11 import (  # noqa: E402
    _apply_min_ntc_batch_fallback,
    _find_first_file,
    _is_plain_single_gene_symbol,
    _pick_column,
)
import dataset_readers_r12 as dr12  # noqa: E402

RAW_ROOT = ROUND12_ROOT / "shared" / "raw"
LINES_ROOT = ROUND12_ROOT / "shared" / "lines"
DEFAULT_FEAS_JSON = ROUND12_ROOT / "shared" / "hek293_feasibility.json"


def build_a549_parts_int64(
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
    """Exact mirror of `dataset_readers_r12.build_a549_parts` with int64 `indptr` cast for NumPy 2.x."""
    adir = dr12._resolve_raw_subdir(raw_dir, ("a549_ko", "a549"))
    cell_path = _find_first_file(adir, ["GSE345058_cell_annotations.csv.gz", "GSE345058_cell_annotations.csv"], "*cell_annotations*.csv*")
    gene_path = _find_first_file(adir, ["GSE345058_gene_annotations.csv.gz", "GSE345058_gene_annotations.csv"], "*gene_annotations*.csv*")
    mtx_path = _find_first_file(adir, ["GSE345058_raw_counts.mtx.gz", "GSE345058_raw_counts.mtx"], "*raw_counts*.mtx*")
    sgrna_candidates = sorted(adir.glob("*sgRNA*.csv*")) + sorted(adir.glob("*sgrna*.csv*"))

    cell_df = pd.read_csv(cell_path, dtype=str, keep_default_na=False)
    gene_df = pd.read_csv(gene_path, dtype=str, keep_default_na=False)

    known_ntc_sgrnas = set()
    if sgrna_candidates:
        sg_df = pd.read_csv(sgrna_candidates[0], dtype=str, keep_default_na=False)
        sg_id_c = _pick_column(sg_df.columns, ["sgRNA", "guide_id", "sgrna", "sequence", "id"])
        sg_pert_c = _pick_column(sg_df.columns, ["Perturbation", "perturbation", "gene", "target", "type"])
        if sg_id_c and sg_pert_c:
            for gval, pval in zip(sg_df[sg_id_c], sg_df[sg_pert_c]):
                if str(pval).strip().lower() in dr12._A549_NTC_TOKENS:
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
        if not sg or sg.lower() in dr12._UNASSIGNED_TOKENS or dr12._MULTI_GUIDE_RE.search(sg) or " " in sg:
            n_multiguide_excl += 1
            continue
        if dr12._MULTI_GUIDE_RE.search(pr):
            n_multiguide_excl += 1
            continue

        b_num = dr12._extract_trailing_int(batches_raw[ci])
        if b_num not in (1, 2, 3, 4):
            n_bad_batch_excl += 1
            continue
        half_flag = 1 if b_num in (1, 3) else 2

        if exclude_zero_cas9 and cas9_vals is not None:
            cv = cas9_vals[ci]
            if not (np.isfinite(cv) and cv > 0.0):
                n_cas9_zero_excl += 1
                continue

        pert_is_ntc = pr.lower() in dr12._A549_NTC_TOKENS
        if conds_raw is not None:
            cond_low = conds_raw[ci].lower()
            if cond_low in dr12._UNASSIGNED_TOKENS:
                n_multiguide_excl += 1
                continue
            cond_is_ntc = cond_low in dr12._A549_NTC_TOKENS
            if pert_is_ntc != cond_is_ntc:
                n_nonplain_excl += 1
                continue
            is_ntc = pert_is_ntc and cond_is_ntc
        else:
            is_ntc = pert_is_ntc
        if is_ntc and known_ntc_sgrnas and sg in known_ntc_sgrnas:
            pass
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

    n_mtx_rows, n_mtx_cols, _, n_skip_lines, sep = dr12._read_mtx_header(mtx_path)
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
                    # Cast indptr to int64 BEFORE addition to avoid NumPy 2.x int32 overflow when cur_nnz >= 2^31-1.
                    indptr_ds[cs + 1 : ce + 1] = cur_nnz + csc_b.indptr[1:].astype(np.int64)
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


def assert_hek293_feasible(feas_json: Path) -> None:
    if not feas_json.exists():
        raise RuntimeError(f"Blocked hek293_xucao prep: missing E3 feasibility JSON {feas_json}")
    data = json.loads(feas_json.read_text(encoding="utf-8"))
    if not data.get("e3_gate_pass", False):
        raise RuntimeError(
            f"Blocked hek293_xucao prep: E3 feasibility gate did not pass ({data.get('e3_verdict', 'UNKNOWN')})"
        )


def main() -> int:
    p = argparse.ArgumentParser(description="Prepare Round 12 standardized line caches")
    p.add_argument(
        "--lines",
        nargs="+",
        required=True,
        choices=["a549_ko", "a549_ko_cas9pos", "hek293_xucao"],
        help="Which line(s) to standardize",
    )
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--chunk-size", type=int, default=8192)
    p.add_argument("--feasibility-json", type=Path, default=DEFAULT_FEAS_JSON)
    p.add_argument("--keep-temp-h5", action="store_true", help="Do not delete .cache_* intermediates")
    args = p.parse_args()

    hgnc = build_hgnc_lookup(DEFAULT_HGNC_PATH)
    LINES_ROOT.mkdir(parents=True, exist_ok=True)

    for ln in args.lines:
        t0 = time.time()
        print(f"[PREP START] {ln} from {RAW_ROOT}...", flush=True)
        cache_dir = LINES_ROOT / f".cache_{ln}"
        try:
            if ln in ("a549_ko", "a549_ko_cas9pos"):
                use_c9 = (ln == "a549_ko_cas9pos")
                parts, audit = build_a549_parts_int64(
                    RAW_ROOT,
                    hgnc,
                    exclude_zero_cas9=use_c9,
                    cache_dir=cache_dir,
                )
                out_dir = process_r10_line(
                    line_name=ln,
                    parts=parts,
                    hgnc=hgnc,
                    chunk_size=args.chunk_size,
                    out_root=LINES_ROOT,
                    extra_meta={"adapter_audit": audit},
                    n_workers=args.workers,
                )
            elif ln == "hek293_xucao":
                assert_hek293_feasible(args.feasibility_json)
                out_dir = dr12.prepare_r12_line(
                    line="hek293_xucao",
                    raw_dir=RAW_ROOT,
                    out_root=LINES_ROOT,
                    hgnc=hgnc,
                    workers=args.workers,
                    chunk_size=args.chunk_size,
                    apply_author_kd40_filter=False,
                )
            else:
                raise ValueError(f"Unsupported line {ln}")
            print(f"[PREP DONE] {ln} -> {out_dir} ({time.time() - t0:.1f}s)", flush=True)
        finally:
            if not args.keep_temp_h5 and cache_dir.exists():
                shutil.rmtree(cache_dir, ignore_errors=True)
                print(f"[CLEANUP] Removed {cache_dir}", flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
