#!/usr/bin/env python3
"""
Clean-Room Round 16 Raw Single-Cell Control Extractor & Standardized Reference Loader (`A-T3`).

Extracts ONLY QC-passing unperturbed `non-targeting` control (`NTC`) cells from raw local
datasets (or explicitly authorized standardized caches if requested), deterministically sampling
the exact `16 × 32` (`n_ctrl = 512`, `seed = 20261601`) stratified control sets BEFORE reading
dense expression rows so peak RAM stays `<< 2 GB` even on multi-million-cell cohorts (`cd4t_rest`,
`a549_ko`, `dld1`, `feng_ipsc`, `nadig_hepg2`, `nadig_jurkat`).

Never reads or imports `round16/src/`, never opens `adata_pred.h5ad`, and never executes
without freeze gates when invoked from `run_audit_r16.py --phase real_r16`.
"""

from __future__ import annotations

import csv
import gzip
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

AUDIT_DIR = Path(__file__).resolve().parent
REPO_ROOT = AUDIT_DIR.parents[1]
R9_SCRIPTS_DIR = REPO_ROOT / "round9" / "p2" / "shared" / "scripts"
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
R10_SCRIPTS_DIR = REPO_ROOT / "round10" / "shared" / "scripts"
R11_SRC_DIR = REPO_ROOT / "round11" / "src"
R12_SRC_DIR = REPO_ROOT / "round12" / "src"

for _p in (
    str(AUDIT_DIR),
    str(R9_SCRIPTS_DIR),
    str(R10_SRC_DIR),
    str(R10_SCRIPTS_DIR),
    str(R11_SRC_DIR),
    str(R12_SRC_DIR),
):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from audit_r16_sampling import (  # noqa: E402
    CELL_SET_LEN_DEFAULT,
    DEFAULT_R16_SEED,
    MIN_CTRL_FLOOR_DEFAULT,
    N_CTRL_DEFAULT,
    sample_stratified_control_sets_r16,
)
from audit_r16_state_adapter import get_hgnc_bridge  # noqa: E402
from data_loader_r10 import (  # noqa: E402
    DEFAULT_HGNC_FILE,
    HGNCBridge,
    R9_LINES_DIR,
    R10_LINES_DIR,
    decode_str_array,
    load_r9_common_genes,
    load_standardized_line_r10,
)
from data_loader_r11 import R11_LINES_DIR  # noqa: E402
from data_loader_r12 import R12_LINES_DIR  # noqa: E402
from dataset_readers_r11 import (  # noqa: E402
    _FENG_DELIM_RE,
    _FENG_KOLF_EXCLUDED,
    _apply_min_ntc_batch_fallback,
    _find_first_file,
    _parse_feng_gene_row_label,
    _pick_column,
)
from dataset_readers_r12 import (  # noqa: E402
    _A549_NTC_TOKENS,
    _MULTI_GUIDE_RE,
    _UNASSIGNED_TOKENS,
    _extract_trailing_int,
    _read_mtx_header,
    _resolve_raw_subdir,
)
from fct_pipeline import (  # noqa: E402
    FctSourceState,
    build_source_input_r10,
    predict_group_c_suite,
)
from prepare_lines import read_obs_or_var_col  # noqa: E402
from prepare_r10_lines import (  # noqa: E402
    CD4_EXCEL_GUIDE_REMAP,
    DLD1_LOW_MOI_CMOS,
    build_cd4t_rest_parts,
    build_hgnc_lookup,
    compute_10x_csc_rna_qc,
    extract_top2_csc,
)

SYNC_NADIG_DEFAULT = Path(
    "/Volumes/雷雳传输/A顶会论文同步盘/A_BIG_MAC顶会论文写稿/AI自主研究/gemini/external_data/nadig2025"
)
STANDARD_LINE_DIRS: Dict[str, Path] = {
    "nadig_hepg2": R9_LINES_DIR / "nadig_hepg2",
    "nadig_jurkat": R9_LINES_DIR / "nadig_jurkat",
    "replogle_k562_essential": R9_LINES_DIR / "replogle_k562_essential",
    "replogle_rpe1": R9_LINES_DIR / "replogle_rpe1",
    "vcc2025_h1": R9_LINES_DIR / "vcc2025_h1",
    "nourreddine_kolf21j": R9_LINES_DIR / "nourreddine_kolf21j",
    "dld1": R10_LINES_DIR / "dld1",
    "cd4t_rest": R10_LINES_DIR / "cd4t_rest",
    "feng_ipsc": R11_LINES_DIR / "feng_ipsc",
    "a549_ko": R12_LINES_DIR / "a549_ko",
}


@dataclass
class SampledRawControlPool:
    """Exact sampled `16 × 32` raw single-cell NTC bundle (`float32`, zero non-NTC cells)."""

    line_name: str
    raw_counts: np.ndarray              # (n_sampled, G_raw) float32 raw integer UMI counts
    gene_names: List[str]               # (G_raw,)
    gene_ids: List[str]                 # (G_raw,)
    batch_labels: np.ndarray            # (n_sampled,) str
    cell_total_umis: np.ndarray         # (n_sampled,) float64
    median_ntc_umi: float               # median total UMI over ALL QC-valid NTCs in line
    sampling_audit: Dict[str, Any]
    is_log1p: bool = False


def resolve_nadig_raw_h5ad(line_name: str, custom_nadig_dir: Optional[Path] = None) -> Path:
    """Locate `GSE264667_{hepg2,jurkat}_raw_singlecell_01.h5ad` across local/sync paths."""
    fname = {
        "nadig_hepg2": "GSE264667_hepg2_raw_singlecell_01.h5ad",
        "nadig_jurkat": "GSE264667_jurkat_raw_singlecell_01.h5ad",
    }[line_name]
    cands: List[Path] = []
    if custom_nadig_dir is not None:
        cands.append(Path(custom_nadig_dir) / fname)
    if os.environ.get("NADIG_DATA_DIR"):
        cands.append(Path(os.environ["NADIG_DATA_DIR"]) / fname)
    cands.extend(
        [
            SYNC_NADIG_DEFAULT / fname,
            Path("/Users/mac/.cache/r9_p2_raw/nadig2025") / fname,
            REPO_ROOT / "external_data" / "nadig2025" / fname,
        ]
    )
    for c in cands:
        if c.exists():
            return c
    raise FileNotFoundError(
        f"Missing raw .h5ad for {line_name} ({fname}). Checked: {[str(c) for c in cands]}. "
        "Mount the Thunderbolt sync drive or pass `--raw-nadig-dir` / `NADIG_DATA_DIR`."
    )


def _extract_nadig_raw_ntc_sample(
    line_name: str,
    custom_nadig_dir: Optional[Path] = None,
    n_ctrl: int = N_CTRL_DEFAULT,
    cell_set_len: int = CELL_SET_LEN_DEFAULT,
    min_ctrl_floor: int = MIN_CTRL_FLOOR_DEFAULT,
    seed: int = DEFAULT_R16_SEED,
) -> SampledRawControlPool:
    """Scan Nadig `.h5ad` `obs` metadata, sample `16 × 32` NTCs stratified by `gem_group`, then slice rows."""
    h5_path = resolve_nadig_raw_h5ad(line_name, custom_nadig_dir=custom_nadig_dir)
    with h5py.File(h5_path, "r") as f:
        obs_g = read_obs_or_var_col(f["obs"], "gene")
        nt_raw_idx = np.flatnonzero([str(g) == "non-targeting" for g in obs_g])
        gems_all = read_obs_or_var_col(f["obs"], "gem_group") if "gem_group" in f["obs"] else ["0"] * len(obs_g)
        nt_batches = np.asarray([str(gems_all[i]) for i in nt_raw_idx], dtype=object)
        if "UMI_count" in f["obs"]:
            umis_all = np.asarray(f["obs/UMI_count"][:], dtype=np.float64)
            nt_umis = umis_all[nt_raw_idx]
        else:
            nt_umis = None

        var_g = f["var"]
        g_names = read_obs_or_var_col(var_g, "gene_name" if "gene_name" in var_g else "_index")
        g_ids = read_obs_or_var_col(var_g, "gene_id") if "gene_id" in var_g else list(g_names)

        samp = sample_stratified_control_sets_r16(
            n_cells=int(nt_raw_idx.size),
            batch_labels=nt_batches,
            n_ctrl=n_ctrl,
            cell_set_len=cell_set_len,
            min_ctrl_floor=min_ctrl_floor,
            seed=seed,
        )
        if samp["excluded"]:
            raise ValueError(f"{line_name} NTC pool excluded: {samp['exclusion_reason']}")

        sel_local = samp["selected_indices"]
        sel_raw_rows = nt_raw_idx[sel_local]

        # Read strictly the selected rows while preserving exact sel_local permutation order
        ord_pos = np.argsort(sel_raw_rows, kind="mergesort")
        sorted_raw = sel_raw_rows[ord_pos]
        x_node = f["X"]
        n_sel = len(sel_raw_rows)
        n_genes = len(g_names)
        mat_sorted = np.empty((n_sel, n_genes), dtype=np.float32)

        if isinstance(x_node, h5py.Dataset):
            for k, r in enumerate(sorted_raw):
                mat_sorted[k] = np.asarray(x_node[int(r), :], dtype=np.float32)
        else:
            enc = str(x_node.attrs.get("encoding-type", "")).lower()
            if "csc" in enc:
                raise NotImplementedError(f"Unexpected CSC layout in Nadig {h5_path}")
            indptr = x_node["indptr"]
            data = x_node["data"]
            indices = x_node["indices"]
            for k, r in enumerate(sorted_raw):
                p0, p1 = int(indptr[int(r)]), int(indptr[int(r) + 1])
                row_dense = np.zeros(n_genes, dtype=np.float32)
                if p1 > p0:
                    row_dense[indices[p0:p1]] = np.asarray(data[p0:p1], dtype=np.float32)
                mat_sorted[k] = row_dense

        raw_counts = np.empty_like(mat_sorted)
        raw_counts[ord_pos] = mat_sorted

    sel_umis = raw_counts.sum(axis=1, dtype=np.float64) if nt_umis is None else nt_umis[sel_local]
    med_umi = float(np.median(nt_umis[nt_umis > 0.0])) if nt_umis is not None else float(np.median(sel_umis))
    return SampledRawControlPool(
        line_name=line_name,
        raw_counts=raw_counts,
        gene_names=[str(x).strip() for x in g_names],
        gene_ids=[str(x).strip() for x in g_ids],
        batch_labels=nt_batches[sel_local].astype(str),
        cell_total_umis=np.asarray(sel_umis, dtype=np.float64),
        median_ntc_umi=med_umi,
        sampling_audit=samp,
        is_log1p=False,
    )


def _extract_a549_raw_ntc_sample(
    raw_dir: Path,
    exclude_zero_cas9: bool = False,
    min_pct_mt: float = 20.0,
    min_ngenes: int = 500,
    n_ctrl: int = N_CTRL_DEFAULT,
    cell_set_len: int = CELL_SET_LEN_DEFAULT,
    min_ctrl_floor: int = MIN_CTRL_FLOOR_DEFAULT,
    seed: int = DEFAULT_R16_SEED,
) -> SampledRawControlPool:
    """Stream `GSE345058_raw_counts.mtx.gz` retaining ONLY pre-QC `Nontargeting` rows (`< 1 GB` RAM)."""
    adir = _resolve_raw_subdir(raw_dir, ("a549_ko", "a549"))
    cell_path = _find_first_file(
        adir, ["GSE345058_cell_annotations.csv.gz", "GSE345058_cell_annotations.csv"], "*cell_annotations*.csv*"
    )
    gene_path = _find_first_file(
        adir, ["GSE345058_gene_annotations.csv.gz", "GSE345058_gene_annotations.csv"], "*gene_annotations*.csv*"
    )
    mtx_path = _find_first_file(
        adir, ["GSE345058_raw_counts.mtx.gz", "GSE345058_raw_counts.mtx"], "*raw_counts*.mtx*"
    )

    cell_df = pd.read_csv(cell_path, dtype=str, keep_default_na=False)
    gene_df = pd.read_csv(gene_path, dtype=str, keep_default_na=False)
    c_batch_col = _pick_column(cell_df.columns, ["batch", "Batch", "batch_id"])
    c_unit_col = _pick_column(cell_df.columns, ["Lane", "lane", "channel", "sample", "batch", "Batch"])
    c_sgrna_col = _pick_column(cell_df.columns, ["sgRNA", "sgrna", "guide_id", "guide"])
    c_pert_col = _pick_column(cell_df.columns, ["Perturbation", "perturbation", "target_gene", "gene_symbol", "gene"])
    c_cond_col = _pick_column(cell_df.columns, ["Condition", "condition"])
    c_cas9_col = _pick_column(cell_df.columns, ["total_counts_CAS9", "cas9_counts", "CAS9_counts", "CAS9"])

    g_sym_col = _pick_column(gene_df.columns, ["gene_symbol", "gene_name", "symbol", "Symbol", "name", "gene"]) or gene_df.columns[0]
    g_ens_col = _pick_column(gene_df.columns, ["ensembl_id", "gene_id", "gene_ids", "Ensembl", "ensembl"])
    g_type_col = _pick_column(gene_df.columns, ["feature_types", "feature_type", "type"])
    g_gen_col = _pick_column(gene_df.columns, ["genome", "species", "reference"])

    raw_g_syms = [str(x).strip() for x in gene_df[g_sym_col].tolist()]
    raw_g_ens = [str(x).strip().split(".")[0] for x in gene_df[g_ens_col].tolist()] if g_ens_col else [""] * len(gene_df)
    raw_g_gen = [str(x).strip().upper() for x in gene_df[g_gen_col].tolist()] if g_gen_col else [""] * len(gene_df)
    g_types = [str(x).strip() for x in gene_df[g_type_col].tolist()] if g_type_col else ["Gene Expression"] * len(gene_df)

    def _is_cas9(ig: int) -> bool:
        return raw_g_syms[ig].upper() == "CAS9" or raw_g_ens[ig].upper() == "CAS9" or raw_g_gen[ig] == "CAS9"

    rna_raw_idx = [i for i, t in enumerate(g_types) if t in ("Gene Expression", "RNA", "") and not _is_cas9(i)]
    if not rna_raw_idx:
        rna_raw_idx = [i for i in range(len(gene_df)) if not _is_cas9(i)]

    raw_to_rna = np.full(len(gene_df), -1, dtype=np.int32)
    rna_names: List[str] = []
    rna_ids: List[str] = []
    for k, ri in enumerate(rna_raw_idx):
        raw_to_rna[ri] = k
        rna_names.append(raw_g_syms[ri])
        rna_ids.append(raw_g_ens[ri] if raw_g_ens[ri] else raw_g_syms[ri])
    mt_mask_rna = np.asarray([g.startswith("MT-") for g in rna_names], dtype=bool)
    cas9_raw_idx = next((i for i in range(len(gene_df)) if _is_cas9(i)), -1)

    batches_raw = cell_df[c_batch_col].astype(str).str.strip().tolist()
    units_raw = cell_df[c_unit_col].astype(str).str.strip().tolist() if c_unit_col else batches_raw
    sgrnas_raw = cell_df[c_sgrna_col].astype(str).str.strip().tolist()
    perts_raw = cell_df[c_pert_col].astype(str).str.strip().tolist()
    conds_raw = cell_df[c_cond_col].astype(str).str.strip().tolist() if c_cond_col else None
    cas9_vals = pd.to_numeric(cell_df[c_cas9_col], errors="coerce").to_numpy(dtype=np.float64) if c_cas9_col else None

    raw_to_ntc = np.full(len(cell_df), -1, dtype=np.int32)
    ntc_units: List[str] = []
    for ci in range(len(cell_df)):
        sg = sgrnas_raw[ci]
        pr = perts_raw[ci]
        if not sg or sg.lower() in _UNASSIGNED_TOKENS or _MULTI_GUIDE_RE.search(sg) or " " in sg or _MULTI_GUIDE_RE.search(pr):
            continue
        if _extract_trailing_int(batches_raw[ci]) not in (1, 2, 3, 4):
            continue
        if exclude_zero_cas9 and cas9_vals is not None:
            cv = cas9_vals[ci]
            if not (np.isfinite(cv) and cv > 0.0):
                continue
        if pr.lower() not in _A549_NTC_TOKENS:
            continue
        if conds_raw is not None:
            cl = conds_raw[ci].lower()
            if cl in _UNASSIGNED_TOKENS or cl not in _A549_NTC_TOKENS:
                continue
        raw_to_ntc[ci] = len(ntc_units)
        ntc_units.append(units_raw[ci])

    n_pre_ntc = len(ntc_units)
    if n_pre_ntc == 0:
        raise RuntimeError(f"Zero A549 pre-QC NTC cells in {adir}")

    n_mtx_rows, n_mtx_cols, _, n_skip_lines, sep = _read_mtx_header(mtx_path)
    cell_first = n_mtx_rows == len(cell_df) and n_mtx_cols == len(gene_df)
    if not cell_first and not (n_mtx_rows == len(gene_df) and n_mtx_cols == len(cell_df)):
        raise RuntimeError(f"MTX shape mismatch {(n_mtx_rows, n_mtx_cols)} in {mtx_path}")

    r_chunks, c_chunks, v_chunks = [], [], []
    cas9_mtx = np.zeros(n_pre_ntc, dtype=np.float32) if (exclude_zero_cas9 and cas9_vals is None) else None
    reader = pd.read_csv(
        mtx_path,
        sep=sep,
        header=None,
        skiprows=n_skip_lines,
        names=["r", "c", "v"],
        dtype={"r": np.int32, "c": np.int32, "v": np.float32},
        chunksize=1_500_000,
        engine="c",
    )
    for ch in reader:
        rr = ch["r"].to_numpy(dtype=np.int64) - 1
        cc = ch["c"].to_numpy(dtype=np.int64) - 1
        vv = ch["v"].to_numpy(dtype=np.float32)
        c_raw = rr if cell_first else cc
        g_raw = cc if cell_first else rr
        nt_i = raw_to_ntc[c_raw]
        if cas9_mtx is not None and cas9_raw_idx >= 0:
            mc9 = (nt_i >= 0) & (g_raw == cas9_raw_idx) & (vv > 0)
            if np.any(mc9):
                np.add.at(cas9_mtx, nt_i[mc9], vv[mc9])
        g_r = raw_to_rna[g_raw]
        m_ok = (nt_i >= 0) & (g_r >= 0) & (vv != 0.0)
        if np.any(m_ok):
            r_chunks.append(nt_i[m_ok].astype(np.int32))
            c_chunks.append(g_r[m_ok].astype(np.int32))
            v_chunks.append(vv[m_ok].astype(np.float32))

    ntc_csr = sp.csr_matrix(
        (
            np.concatenate(v_chunks) if v_chunks else np.empty(0, dtype=np.float32),
            (
                np.concatenate(r_chunks) if r_chunks else np.empty(0, dtype=np.int32),
                np.concatenate(c_chunks) if c_chunks else np.empty(0, dtype=np.int32),
            ),
        ),
        shape=(n_pre_ntc, len(rna_names)),
    )
    umi_pre = np.asarray(ntc_csr.sum(axis=1), dtype=np.float64).ravel()
    ng_pre = np.diff(ntc_csr.indptr)
    mt_umi = np.asarray(ntc_csr[:, mt_mask_rna].sum(axis=1), dtype=np.float64).ravel() if np.any(mt_mask_rna) else np.zeros(n_pre_ntc)
    pct_mt = 100.0 * mt_umi / np.maximum(umi_pre, 1.0)
    qc_ok = (umi_pre > 0.0) & (ng_pre >= min_ngenes) & (pct_mt <= min_pct_mt)
    if cas9_mtx is not None:
        qc_ok &= cas9_mtx > 0.0

    keep_rows = np.flatnonzero(qc_ok)
    ntc_csr_qc = ntc_csr[keep_rows].tocsr()
    umi_qc = umi_pre[keep_rows]
    units_qc = np.asarray([ntc_units[int(i)] for i in keep_rows], dtype=object)

    samp = sample_stratified_control_sets_r16(
        n_cells=int(len(keep_rows)),
        batch_labels=units_qc,
        n_ctrl=n_ctrl,
        cell_set_len=cell_set_len,
        min_ctrl_floor=min_ctrl_floor,
        seed=seed,
    )
    if samp["excluded"]:
        raise ValueError(f"a549_ko NTC pool excluded: {samp['exclusion_reason']}")

    sel = samp["selected_indices"]
    raw_counts = ntc_csr_qc[sel].toarray().astype(np.float32)
    return SampledRawControlPool(
        line_name="a549_ko",
        raw_counts=raw_counts,
        gene_names=rna_names,
        gene_ids=rna_ids,
        batch_labels=units_qc[sel].astype(str),
        cell_total_umis=umi_qc[sel],
        median_ntc_umi=float(np.median(umi_qc[umi_qc > 0.0])),
        sampling_audit=samp,
        is_log1p=False,
    )


def _extract_dld1_raw_ntc_sample(
    raw_dir: Path,
    n_ctrl: int = N_CTRL_DEFAULT,
    cell_set_len: int = CELL_SET_LEN_DEFAULT,
    min_ctrl_floor: int = MIN_CTRL_FLOOR_DEFAULT,
    seed: int = DEFAULT_R16_SEED,
) -> SampledRawControlPool:
    """Scan all 24 DLD-1 channels for QC-valid NTCs, sample `16 × 32`, then slice only selected cols."""
    ddir = raw_dir / "dld1"
    rna_files = sorted(ddir.glob("GSM98642*_LIB*_SAM24482910.filtered_feature_bc_matrix.h5"))
    hash_files = sorted(ddir.glob("GSM98642*_LIB*_SAM24482910.hashing.umi_correct.h5"))
    crp_files = sorted(ddir.glob("GSM98642*_LIB*_SAM24482910.crispr.umi_correct.h5"))
    if len(rna_files) != 24 or len(hash_files) != 24 or len(crp_files) != 24:
        raise RuntimeError(f"Expected 24 DLD-1 HDF5 trios in {ddir}")

    ntc_records: List[Tuple[Path, int, float, str]] = []
    ref_names: Optional[List[str]] = None
    ref_ids: Optional[List[str]] = None

    for ch_idx, (r_path, h_path, c_path) in enumerate(zip(rna_files, hash_files, crp_files), start=1):
        with h5py.File(r_path, "r") as fr:
            rm = fr["matrix"]
            n_genes_raw = int(rm["shape"][0])
            rna_bcs = [x.decode().split("-")[0] for x in rm["barcodes"][:]]
            var_ids = [x.decode() for x in rm["features/id"][:]]
            var_names = [x.decode() for x in rm["features/name"][:]]
        if ref_names is None:
            ref_names = var_names
            ref_ids = var_ids
        mt_raw = np.asarray([g.startswith("MT-") for g in var_names], dtype=bool)
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
        ch_label = f"ch_{ch_idx:02d}"
        for loc_j, r_i in enumerate(h_pass_r):
            if not sing_c[loc_j]:
                continue
            gname = gnames[int(c_a1[loc_j])]
            if gname.startswith("NTC_") and not gname.startswith("NTC_MENG"):
                ntc_records.append((r_path, int(r_i), float(rna_umi[r_i]), ch_label))

    assert ref_names is not None and ref_ids is not None
    all_umis = np.asarray([r[2] for r in ntc_records], dtype=np.float64)
    all_batches = np.asarray([r[3] for r in ntc_records], dtype=object)
    samp = sample_stratified_control_sets_r16(
        n_cells=len(ntc_records),
        batch_labels=all_batches,
        n_ctrl=n_ctrl,
        cell_set_len=cell_set_len,
        min_ctrl_floor=min_ctrl_floor,
        seed=seed,
    )
    if samp["excluded"]:
        raise ValueError(f"dld1 NTC pool excluded: {samp['exclusion_reason']}")

    sel = samp["selected_indices"]
    by_file: Dict[Path, List[Tuple[int, int]]] = defaultdict(list)
    for pos_out, rec_i in enumerate(sel.tolist()):
        r_path, raw_col, _, _ = ntc_records[rec_i]
        by_file[r_path].append((pos_out, raw_col))

    n_genes = len(ref_names)
    raw_counts = np.zeros((len(sel), n_genes), dtype=np.float32)
    for r_path, items in by_file.items():
        with h5py.File(r_path, "r") as fr:
            m = fr["matrix"]
            indptr = m["indptr"]
            data = m["data"]
            indices = m["indices"]
            for pos_out, col_i in items:
                p0, p1 = int(indptr[col_i]), int(indptr[col_i + 1])
                if p1 > p0:
                    raw_counts[pos_out, indices[p0:p1]] = np.asarray(data[p0:p1], dtype=np.float32)

    return SampledRawControlPool(
        line_name="dld1",
        raw_counts=raw_counts,
        gene_names=ref_names,
        gene_ids=ref_ids,
        batch_labels=all_batches[sel].astype(str),
        cell_total_umis=all_umis[sel],
        median_ntc_umi=float(np.median(all_umis[all_umis > 0.0])),
        sampling_audit=samp,
        is_log1p=False,
    )


def _extract_cd4t_rest_raw_ntc_sample(
    raw_dir: Path,
    hgnc_bridge: Optional[HGNCBridge] = None,
    n_ctrl: int = N_CTRL_DEFAULT,
    cell_set_len: int = CELL_SET_LEN_DEFAULT,
    min_ctrl_floor: int = MIN_CTRL_FLOOR_DEFAULT,
    seed: int = DEFAULT_R16_SEED,
) -> SampledRawControlPool:
    """Use frozen `build_cd4t_rest_parts` to index all QC-valid NTCs across 47 lanes, then slice 512."""
    hgnc_l = build_hgnc_lookup()
    parts, _ = build_cd4t_rest_parts(raw_dir, hgnc_l)
    nt_recs: List[Tuple[int, int, str, float]] = []
    for p_idx, p in enumerate(parts):
        nt_local = np.flatnonzero(p.nt_mask)
        for loc_i in nt_local.tolist():
            nt_recs.append((p_idx, int(loc_i), str(p.batches[loc_i]), float(p.umi[loc_i])))

    all_batches = np.asarray([r[2] for r in nt_recs], dtype=object)
    all_umis = np.asarray([r[3] for r in nt_recs], dtype=np.float64)
    samp = sample_stratified_control_sets_r16(
        n_cells=len(nt_recs),
        batch_labels=all_batches,
        n_ctrl=n_ctrl,
        cell_set_len=cell_set_len,
        min_ctrl_floor=min_ctrl_floor,
        seed=seed,
    )
    if samp["excluded"]:
        raise ValueError(f"cd4t_rest NTC pool excluded: {samp['exclusion_reason']}")

    sel = samp["selected_indices"]
    by_part: Dict[int, List[Tuple[int, int]]] = defaultdict(list)
    for pos_out, rec_i in enumerate(sel.tolist()):
        p_idx, loc_i, _, _ = nt_recs[rec_i]
        by_part[p_idx].append((pos_out, loc_i))

    ref_names = [str(x) for x in parts[0].var_gene_names]
    ref_ids = [str(x) for x in parts[0].var_gene_ids]
    raw_counts = np.zeros((len(sel), len(ref_names)), dtype=np.float32)

    for p_idx, items in by_part.items():
        p = parts[p_idx]
        items_sorted = sorted(items, key=lambda kv: kv[1])
        loc_arr = np.asarray([kv[1] for kv in items_sorted], dtype=np.int64)
        out_arr = np.asarray([kv[0] for kv in items_sorted], dtype=np.int64)
        assert p.read_rows_fn is not None
        rows_mat = np.asarray(p.read_rows_fn(loc_arr), dtype=np.float32)
        if [str(x) for x in p.var_gene_names] == ref_names:
            raw_counts[out_arr] = rows_mat
        else:
            pm = {str(g): j for j, g in enumerate(p.var_gene_names)}
            cidx = [pm[g] for g in ref_names if g in pm]
            oidx = [i for i, g in enumerate(ref_names) if g in pm]
            raw_counts[np.ix_(out_arr, oidx)] = rows_mat[:, cidx]

    return SampledRawControlPool(
        line_name="cd4t_rest",
        raw_counts=raw_counts,
        gene_names=ref_names,
        gene_ids=ref_ids,
        batch_labels=all_batches[sel].astype(str),
        cell_total_umis=all_umis[sel],
        median_ntc_umi=float(np.median(all_umis[all_umis > 0.0])),
        sampling_audit=samp,
        is_log1p=False,
    )


def _extract_feng_ipsc_raw_ntc_sample(
    raw_dir: Path,
    n_ctrl: int = N_CTRL_DEFAULT,
    cell_set_len: int = CELL_SET_LEN_DEFAULT,
    min_ctrl_floor: int = MIN_CTRL_FLOOR_DEFAULT,
    seed: int = DEFAULT_R16_SEED,
) -> SampledRawControlPool:
    """
    Stream `TargetedScreen_RNA-UMI-Counts.csv.gz` retaining ONLY non-KOLF `NonTarget_*` NTC columns,
    compute total UMIs & stratified `16 × 32` selection deterministically.
    """
    fdir = raw_dir / "feng_ipsc"
    meta_path = _find_first_file(fdir, ["TargetedScreen_Cell-Metadata.tsv.gz", "TargetedScreen_Cell-Metadata.tsv"], "*Cell-Metadata*.tsv*")
    csv_path = _find_first_file(fdir, ["TargetedScreen_RNA-UMI-Counts.csv.gz", "TargetedScreen_RNA-UMI-Counts.csv"], "*RNA-UMI-Counts*.csv*")

    meta_df = pd.read_csv(meta_path, sep="\t", dtype=str, keep_default_na=False)
    cid_col = _pick_column(meta_df.columns, ["Cell_ID", "cell_id", "barcode", "CBC"])
    cline_col = _pick_column(meta_df.columns, ["Cell_Line", "cell_line", "line"])
    gcall_col = _pick_column(meta_df.columns, ["Guide_Call", "guide_call", "sgrna", "perturbation"])
    inlet_col = _pick_column(meta_df.columns, ["Inlet", "inlet", "batch", "lane"])

    nt_meta: Dict[str, Tuple[str, str]] = {}
    kept_lines: Set[str] = set()
    for r in meta_df.itertuples(index=False):
        cid = str(getattr(r, cid_col)).strip()
        cline = str(getattr(r, cline_col)).strip()
        gcall = str(getattr(r, gcall_col)).strip()
        if not cid or not cline:
            continue
        if cline.lower() in _FENG_KOLF_EXCLUDED or cline.lower().startswith("kolf_2") or cline.lower().startswith("kolf_3"):
            continue
        lg = gcall.lower()
        if not gcall or lg in ("unassigned", "na", "nan", "none") or _FENG_DELIM_RE.search(gcall):
            continue
        if not (gcall.startswith("NonTarget_") or lg.startswith("nontarget")):
            continue
        inlet = str(getattr(r, inlet_col)).strip() if inlet_col else "inlet0"
        kept_lines.add(cline)
        nt_meta[cid] = (cline, f"{cline}__inlet{inlet}")

    lines_sorted = sorted(kept_lines)
    line_to_half = {ln: (1 if i % 2 == 0 else 2) for i, ln in enumerate(lines_sorted)}

    opener = gzip.open if str(csv_path).endswith(".gz") else open
    with opener(csv_path, "rt", encoding="utf-8", errors="replace") as f_csv:
        reader = csv.reader(f_csv)
        header = [h.strip() for h in next(reader)]
        first_row = next(reader, None)
        if first_row is None:
            raise RuntimeError(f"Empty data rows in {csv_path}")
        col_offset = len(first_row) - len(header)
        hdr_start = 1 if col_offset == 0 else 0

        keep_csv_cols: List[int] = []
        unit_list: List[str] = []
        line_list: List[str] = []
        half_list: List[int] = []
        for h_idx in range(hdr_start, len(header)):
            cid = header[h_idx]
            rec = nt_meta.get(cid)
            if rec is None:
                continue
            cline, uname = rec
            keep_csv_cols.append(h_idx + col_offset)
            unit_list.append(uname)
            line_list.append(cline)
            half_list.append(line_to_half[cline])

        n_ntc = len(keep_csv_cols)
        if n_ntc == 0:
            raise RuntimeError(f"Zero kept non-KOLF NTC columns in {csv_path}")

        nt_arr = np.ones(n_ntc, dtype=bool)
        half_arr = np.asarray(half_list, dtype=np.int8)
        batch_arr, _ = _apply_min_ntc_batch_fallback(
            np.asarray(unit_list, dtype=str),
            nt_arr,
            half_arr,
            group_labels=np.asarray(line_list, dtype=str),
            min_batch_nt=15,
        )
        samp = sample_stratified_control_sets_r16(
            n_cells=n_ntc,
            batch_labels=batch_arr,
            n_ctrl=n_ctrl,
            cell_set_len=cell_set_len,
            min_ctrl_floor=min_ctrl_floor,
            seed=seed,
        )
        if samp["excluded"]:
            raise ValueError(f"feng_ipsc NTC pool excluded: {samp['exclusion_reason']}")

        sel_idx = samp["selected_indices"]
        sel_csv_cols = np.asarray([keep_csv_cols[i] for i in sel_idx], dtype=np.int64)
        all_csv_cols = np.asarray(keep_csv_cols, dtype=np.int64)

        umi_all = np.zeros(n_ntc, dtype=np.float64)
        gene_names: List[str] = []
        gene_ids: List[str] = []
        sel_rows: List[np.ndarray] = []

        def _consume(row_vals: List[str]) -> None:
            if not row_vals:
                return
            eid, sym, is_gex = _parse_feng_gene_row_label(row_vals[0])
            if not is_gex or not sym:
                return
            all_v = np.fromiter((float(row_vals[j]) for j in all_csv_cols), dtype=np.float64, count=n_ntc)
            umi_all[:] += all_v
            sel_v = all_v[sel_idx].astype(np.float32)
            sel_rows.append(sel_v)
            gene_names.append(sym)
            gene_ids.append(eid or sym)

        _consume(first_row)
        for r in reader:
            _consume(r)

    raw_counts = np.stack(sel_rows, axis=1).astype(np.float32)  # (n_sampled, G_raw)
    return SampledRawControlPool(
        line_name="feng_ipsc",
        raw_counts=raw_counts,
        gene_names=gene_names,
        gene_ids=gene_ids,
        batch_labels=batch_arr[sel_idx].astype(str),
        cell_total_umis=umi_all[sel_idx],
        median_ntc_umi=float(np.median(umi_all[umi_all > 0.0])),
        sampling_audit=samp,
        is_log1p=False,
    )


def extract_raw_control_pool_r16(
    line_name: str,
    custom_nadig_dir: Optional[Path] = None,
    n_ctrl: int = N_CTRL_DEFAULT,
    cell_set_len: int = CELL_SET_LEN_DEFAULT,
    min_ctrl_floor: int = MIN_CTRL_FLOOR_DEFAULT,
    seed: int = DEFAULT_R16_SEED,
    allow_npz_fallback_for_hvg: bool = False,
) -> SampledRawControlPool:
    """Dispatch exact raw NTC single-cell extraction for all 6 Gate-2 local raw-count cohorts."""
    if line_name in ("nadig_hepg2", "nadig_jurkat"):
        try:
            return _extract_nadig_raw_ntc_sample(
                line_name=line_name,
                custom_nadig_dir=custom_nadig_dir,
                n_ctrl=n_ctrl,
                cell_set_len=cell_set_len,
                min_ctrl_floor=min_ctrl_floor,
                seed=seed,
            )
        except FileNotFoundError:
            if not allow_npz_fallback_for_hvg:
                raise
            ld = STANDARD_LINE_DIRS[line_name]
            genes = [x.strip() for x in (ld / "genes.txt").read_text(encoding="utf-8").splitlines() if x.strip()]
            ctrl_npz = np.load(ld / "ctrl_cells.npz", allow_pickle=False)
            batches = np.asarray([str(x) for x in ctrl_npz["batch_labels"]], dtype=object)
            x_log1p = np.asarray(ctrl_npz["X_log1p"], dtype=np.float32)
            samp = sample_stratified_control_sets_r16(
                n_cells=x_log1p.shape[0],
                batch_labels=batches,
                n_ctrl=n_ctrl,
                cell_set_len=cell_set_len,
                min_ctrl_floor=min_ctrl_floor,
                seed=seed,
            )
            sel = samp["selected_indices"]
            return SampledRawControlPool(
                line_name=line_name,
                raw_counts=x_log1p[sel],
                gene_names=genes,
                gene_ids=genes,
                batch_labels=batches[sel].astype(str),
                cell_total_umis=np.full(len(sel), 1e4, dtype=np.float64),
                median_ntc_umi=1e4,
                sampling_audit=samp,
                is_log1p=True,
            )

    if line_name == "a549_ko":
        return _extract_a549_raw_ntc_sample(
            REPO_ROOT / "round12" / "shared" / "raw",
            n_ctrl=n_ctrl,
            cell_set_len=cell_set_len,
            min_ctrl_floor=min_ctrl_floor,
            seed=seed,
        )
    if line_name == "dld1":
        return _extract_dld1_raw_ntc_sample(
            REPO_ROOT / "round10" / "shared" / "raw",
            n_ctrl=n_ctrl,
            cell_set_len=cell_set_len,
            min_ctrl_floor=min_ctrl_floor,
            seed=seed,
        )
    if line_name == "cd4t_rest":
        return _extract_cd4t_rest_raw_ntc_sample(
            REPO_ROOT / "round10" / "shared" / "raw",
            n_ctrl=n_ctrl,
            cell_set_len=cell_set_len,
            min_ctrl_floor=min_ctrl_floor,
            seed=seed,
        )
    if line_name == "feng_ipsc":
        return _extract_feng_ipsc_raw_ntc_sample(
            REPO_ROOT / "round11" / "shared" / "raw",
            n_ctrl=n_ctrl,
            cell_set_len=cell_set_len,
            min_ctrl_floor=min_ctrl_floor,
            seed=seed,
        )
    raise ValueError(f"Unsupported local raw control source {line_name!r}")


def load_line_eval_reference_r16(
    target_line: str,
    model_genes: Sequence[str],
    model_pert_names: Sequence[str],
    hgnc_bridge: Optional[HGNCBridge] = None,
    require_self_gene_in_eval: bool = False,
) -> Dict[str, Any]:
    """
    Build standardized Round 16 reference evaluation data for `target_line`:
      - Harmonizes `G_R9_6148`, `target_line` genes, source genes, and `G_model` (`2,000` HVGs)
        through frozen R10 `HGNCBridge`.
      - Locks `N_eval^R16(t)`: frozen `P5`-qualified reproducible targets in `canon(P_model)`.
        Per `R16_GATE2_SCHEMA.md §2.3`, does NOT drop targets whose self-gene is outside `G_eval^R16(t)`
        (uses `target_gene_idx = -1` unless `require_self_gene_in_eval=True`).
      - Computes zero-leakage source-cohort `Transfer_Mean`, `source_grand_mean`, and frozen `FCT` (`H-I2`).
    """
    hb = get_hgnc_bridge(hgnc_bridge)
    r9_common_canon = load_r9_common_genes(hgnc_bridge=hb)

    model_g_canon: List[str] = []
    seen_mg: Set[str] = set()
    for g in model_genes:
        cg = hb.canonicalize(str(g))
        if cg and cg not in seen_mg:
            seen_mg.add(cg)
            model_g_canon.append(cg)

    model_p_canon: Set[str] = set()
    for p in model_pert_names:
        if str(p) == "non-targeting":
            continue
        cp = hb.canonicalize(str(p))
        if cp:
            model_p_canon.add(cp)

    src_lines = (
        sorted(set(["nadig_hepg2", "nadig_jurkat", "replogle_k562_essential", "replogle_rpe1"]) - {target_line})
        if target_line in ("nadig_hepg2", "nadig_jurkat")
        else ["nadig_hepg2", "nadig_jurkat", "nourreddine_kolf21j", "replogle_k562_essential", "replogle_rpe1", "vcc2025_h1"]
    )

    tgt_ds = load_standardized_line_r10(
        STANDARD_LINE_DIRS[target_line],
        representation="log1p",
        common_genes=r9_common_canon,
        hgnc_bridge=hb,
        min_gene_overlap=1000,
    )
    src_ds_map = {
        sn: load_standardized_line_r10(
            STANDARD_LINE_DIRS[sn],
            representation="log1p",
            common_genes=r9_common_canon,
            hgnc_bridge=hb,
            min_gene_overlap=1000,
        )
        for sn in src_lines
    }

    # Exact intersection G_eval^R16(t) across R9_COMMON, all source lines, target line, and model output genes
    shared_set = set(tgt_ds.genes) & set(model_g_canon)
    for sd in src_ds_map.values():
        shared_set &= set(sd.genes)
    eval_genes = [g for g in r9_common_canon if g in shared_set]
    eval_g_to_idx = {g: i for i, g in enumerate(eval_genes)}

    def _align_ds_cols(ds: Any) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        pos = {g: i for i, g in enumerate(ds.genes)}
        cols = np.asarray([pos[g] for g in eval_genes], dtype=np.int64)
        t_idx = np.asarray([eval_g_to_idx.get(t, -1) for t in ds.target_names], dtype=np.int64)
        return (
            np.asarray(ds.deltas[:, cols], dtype=np.float64),
            np.asarray(ds.ctrl_cells[:, cols], dtype=np.float32),
            np.asarray(ds.ctrl_mean_libnorm[cols], dtype=np.float64),
            t_idx,
            np.asarray(ds.half_a_deltas[:, cols], dtype=np.float64),
            np.asarray(ds.half_b_deltas[:, cols], dtype=np.float64),
            cols,
        )

    tgt_deltas, _, tgt_cmean_lib, tgt_gidx_full, tgt_ha, tgt_hb, _ = _align_ds_cols(tgt_ds)
    p5_gene_set = set(tgt_ds.genes)
    src_inputs = []
    for sn in src_lines:
        sd = src_ds_map[sn]
        p5_gene_set &= set(sd.genes)
        pos = {g: i for i, g in enumerate(sd.genes)}
        cols = np.asarray([pos[g] for g in eval_genes], dtype=np.int64)
        # Subset LineDatasetR10 fields via shallow copy for build_source_input_r10
        from data_loader_r10 import LineDatasetR10

        sub_sd = LineDatasetR10(
            name=sd.name,
            genes=eval_genes,
            target_names=sd.target_names,
            target_indices=np.asarray([eval_g_to_idx.get(t, -1) for t in sd.target_names], dtype=np.int64),
            deltas=sd.deltas[:, cols],
            ctrl_cells=sd.ctrl_cells[:, cols],
            ctrl_mean_log1p=sd.ctrl_mean_log1p[cols] if sd.ctrl_mean_log1p is not None else None,
            ctrl_mean_libnorm=sd.ctrl_mean_libnorm[cols] if sd.ctrl_mean_libnorm is not None else None,
            ctrl_mean_libnorm_full=sd.ctrl_mean_libnorm_full[cols] if sd.ctrl_mean_libnorm_full is not None else None,
            repr_mask=sd.repr_mask,
            half_a_deltas=sd.half_a_deltas[:, cols] if sd.half_a_deltas is not None else None,
            half_b_deltas=sd.half_b_deltas[:, cols] if sd.half_b_deltas is not None else None,
            half_valid=sd.half_valid,
            n_cells=sd.n_cells,
            knockdown_depth_log1p=sd.knockdown_depth_log1p,
            knockdown_fraction=sd.knockdown_fraction,
            target_ensembl_ids=sd.target_ensembl_ids,
            representation=sd.representation,
            metadata=sd.metadata,
        )
        src_inputs.append(build_source_input_r10(sub_sd, only_reproducible=True, compute_cov=False))

    fct_state = FctSourceState.fit(src_inputs, mask_self=True)
    src_rep_targets = set().union(*(set(si.target_names) for si in src_inputs))

    # Exact R9-R12 P5 self-gene eligibility over G_eval^P5(t) before intersecting w/ canon(P_model)
    row_sel: List[int] = []
    eval_perts: List[str] = []
    eval_target_gidx: List[int] = []
    for ri, t_sym in enumerate(tgt_ds.target_names):
        if tgt_ds.repr_mask is not None and not bool(tgt_ds.repr_mask[ri]):
            continue
        if t_sym not in src_rep_targets:
            continue
        if t_sym not in model_p_canon:
            continue
        if require_self_gene_in_eval and t_sym not in eval_g_to_idx:
            continue
        if not require_self_gene_in_eval and t_sym not in p5_gene_set:
            continue
        row_sel.append(ri)
        eval_perts.append(t_sym)
        eval_target_gidx.append(int(tgt_gidx_full[ri]))

    row_arr = np.asarray(row_sel, dtype=np.int64)
    q_gidx_arr = np.asarray(eval_target_gidx, dtype=np.int64)
    preds_fct, _ = predict_group_c_suite(
        state=fct_state,
        target_ctrl_mean_libnorm=tgt_cmean_lib,
        query_target_names=eval_perts,
        query_target_gene_indices=q_gidx_arr,
        include_r9_anchors=False,
    )

    return {
        "target_line": target_line,
        "eval_genes": eval_genes,
        "eval_perts": eval_perts,
        "target_gene_idx": q_gidx_arr,
        "delta_obs": tgt_deltas[row_arr],
        "delta_half_a": tgt_ha[row_arr],
        "delta_half_b": tgt_hb[row_arr],
        "half_valid": np.asarray(tgt_ds.half_valid[row_arr], dtype=bool),
        "delta_transfer": np.asarray(preds_fct["Transfer_Mean"], dtype=np.float64),
        "source_grand_mean": np.asarray(fct_state.grand_mean, dtype=np.float64),
        "delta_fct": np.asarray(preds_fct["FCT"], dtype=np.float64),
        "ctrl_mean_libnorm": tgt_cmean_lib,
        "source_lines_fct": src_lines,
    }
