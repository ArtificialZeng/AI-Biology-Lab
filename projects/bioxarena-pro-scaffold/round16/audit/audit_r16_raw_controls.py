#!/usr/bin/env python3
"""
Clean-Room Round 16 Raw Single-Cell Control Extractor & Standardized Reference Loader (`A-T3`,
strictly aligned to `PREREG_I.md v5` + `Addendum A`).

Features:
  - Extracts ONLY QC-passing unperturbed `non-targeting` control (`NTC`) cells from raw local
    datasets (`nadig_hepg2`, `nadig_jurkat`, `a549_ko`, `dld1`, `cd4t_rest`, `feng_ipsc`) and
    pre-download locked prospective Group 3 (`gse295214_h1` `85` GSM CellRanger `.h5` files).
  - Deterministically samples the `16 × 32` (`n_ctrl = 512`, `seed = 20261601`) stratified primary
    control sets AND optional strictly non-overlapping `seed = 20261602` same-line resample control
    sets (`Addendum A Item 8`) BEFORE reading dense expression rows.
  - Computes `median_ntc_umi` STRICTLY over the sampled control pool's full measured transcriptome
    (`PREREG_I v5 §2.1.4 L83`) and tracks exact sampled `cell_barcodes` for 1-to-1 verification.
  - Builds frozen zero-leakage reference bundles (`r9_6` for Group 2; `r9_6 \\ {vcc2025_h1}` primary
    for Group 3 + locked-target `full r9_6` secondary sensitivity per `Addendum A`).

Never reads or imports `round16/src/`, never downloads data, and never opens `adata_pred.h5ad`.
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
    RESAMPLE_R16_SEED,
    sample_nonoverlapping_resample_sets_r16,
    sample_stratified_control_sets_r16,
)
from audit_r16_state_adapter import get_hgnc_bridge  # noqa: E402
from data_loader_r10 import (  # noqa: E402
    DEFAULT_HGNC_FILE,
    HGNCBridge,
    LineDatasetR10,
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
    _is_plain_single_gene_symbol,
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
    DLD1_LOW_MOI_CMOS,
    HGNCLookup,
    PartSpec,
    build_cd4t_rest_parts,
    build_hgnc_lookup,
    compute_10x_csc_rna_qc,
    extract_top2_csc,
    fit_poisson_gaussian_thresholds_from_hist,
    make_10x_csc_subset_readers,
    process_r10_line,
)

SYNC_NADIG_DEFAULT = Path(
    "/Volumes/雷雳传输/A顶会论文同步盘/A_BIG_MAC顶会论文写稿/AI自主研究/gemini/external_data/nadig2025"
)
R16_LINES_DIR = REPO_ROOT / "round16" / "shared" / "lines"
R16_RAW_DIR = REPO_ROOT / "round16" / "shared" / "raw"

R9_6_ALL_LINES: Tuple[str, ...] = (
    "nadig_hepg2",
    "nadig_jurkat",
    "nourreddine_kolf21j",
    "replogle_k562_essential",
    "replogle_rpe1",
    "vcc2025_h1",
)
GROUP3_PRIMARY_FCT_SOURCES: Tuple[str, ...] = (
    "nadig_hepg2",
    "nadig_jurkat",
    "nourreddine_kolf21j",
    "replogle_k562_essential",
    "replogle_rpe1",
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
    "cd4t_stim48": R11_LINES_DIR / "cd4t_stim48",
    "feng_ipsc": R11_LINES_DIR / "feng_ipsc",
    "a549_ko": R12_LINES_DIR / "a549_ko",
    "a549_ko_cas9pos": R12_LINES_DIR / "a549_ko_cas9pos",
    "a549_myllia": REPO_ROOT / "round14" / "shared" / "lines" / "a549_myllia",
    "hek293_xucao": R12_LINES_DIR / "hek293_xucao",
    "gse295214_h1": R16_LINES_DIR / "gse295214_h1",
}

_GSE295214_ILLEGAL_TOK_RE = re.compile(r"[,;+\s]")
_GSE295214_ENSG_RE = re.compile(r"(ENSG\d+)", re.IGNORECASE)
_GSE295214_IGNORE_SUFFIX_SEG_RE = re.compile(
    r"^(?:P\d.*|ENST\d+|ENSG\d+|\d+|sg\d+|TSS\d*|[+-]\d+.*|AB|A|B|[ACGT]{15,})$",
    re.IGNORECASE,
)
_GSE295214_HYPHEN_SUFFIX_RE = re.compile(
    r"^(.+?)-(?:\d+|sg\d+|P\d+.*|TSS\d*)$",
    re.IGNORECASE,
)
_GSE295214_BLACKLIST_SUBSTRS: Tuple[str, ...] = (
    "safe",
    "intergenic",
    "unassigned",
    "multiplet",
    "doublet",
    "gene_desert",
)


@dataclass
class SampledRawControlPool:
    """
    Sampled raw single-cell NTC bundle (`float32`, zero non-NTC cells).
    `median_ntc_umi` strictly reflects the sampled pool (`PREREG_I v5 §2.1.4 L83`).
    """

    line_name: str
    raw_counts: np.ndarray              # (n_sampled, G_raw) float32 raw integer UMI counts
    gene_names: List[str]               # (G_raw,)
    gene_ids: List[str]                 # (G_raw,)
    batch_labels: np.ndarray            # (n_sampled,) str
    cell_barcodes: np.ndarray           # (n_sampled,) str unique barcodes in exact sampled order
    cell_total_umis: np.ndarray         # (n_sampled,) float64
    median_ntc_umi: float               # median total UMI over SAMPLED pool (v5 §2.1.4 L83)
    sampling_audit: Dict[str, Any]
    is_log1p: bool = False
    median_ntc_umi_all_valid: Optional[float] = None


def _read_obs_cell_barcodes(obs_grp: h5py.Group, n_expected: int, line_name: str) -> np.ndarray:
    for col in ("_index", "index", "cell_id", "barcode", "Barcode", "Cell_ID"):
        if col in obs_grp:
            vals = read_obs_or_var_col(obs_grp, col)
            if len(vals) == n_expected:
                return np.asarray([str(x).strip() for x in vals], dtype=object)
    return np.asarray([f"{line_name}_row_{i}" for i in range(n_expected)], dtype=object)


def _median_positive(arr: np.ndarray, fallback: float = 1.0) -> float:
    pos = np.asarray(arr, dtype=np.float64)
    pos = pos[np.isfinite(pos) & (pos > 0.0)]
    return float(np.median(pos)) if pos.size > 0 else float(fallback)


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
    return_resample: bool = False,
    resample_seed: int = RESAMPLE_R16_SEED,
) -> SampledRawControlPool | Tuple[SampledRawControlPool, Optional[SampledRawControlPool]]:
    """Scan Nadig `.h5ad` `obs`, sample primary + optional non-overlapping resample `16 × 32`, slice rows."""
    h5_path = resolve_nadig_raw_h5ad(line_name, custom_nadig_dir=custom_nadig_dir)
    with h5py.File(h5_path, "r") as f:
        obs_g = read_obs_or_var_col(f["obs"], "gene")
        all_obs_bcs = _read_obs_cell_barcodes(f["obs"], len(obs_g), line_name)
        nt_raw_idx = np.flatnonzero([str(g) == "non-targeting" for g in obs_g])
        gems_all = read_obs_or_var_col(f["obs"], "gem_group") if "gem_group" in f["obs"] else ["0"] * len(obs_g)
        nt_batches = np.asarray([str(gems_all[i]) for i in nt_raw_idx], dtype=object)
        nt_bcs = all_obs_bcs[nt_raw_idx]
        if "UMI_count" in f["obs"]:
            umis_all = np.asarray(f["obs/UMI_count"][:], dtype=np.float64)
            nt_umis = umis_all[nt_raw_idx]
        else:
            nt_umis = None

        var_g = f["var"]
        g_names = read_obs_or_var_col(var_g, "gene_name" if "gene_name" in var_g else "_index")
        g_ids = read_obs_or_var_col(var_g, "gene_id") if "gene_id" in var_g else list(g_names)
        g_names_clean = [str(x).strip() for x in g_names]
        g_ids_clean = [str(x).strip() for x in g_ids]

        samp = sample_stratified_control_sets_r16(
            n_cells=int(nt_raw_idx.size),
            batch_labels=nt_batches,
            n_ctrl=n_ctrl,
            cell_set_len=cell_set_len,
            min_ctrl_floor=min_ctrl_floor,
            seed=seed,
            cell_barcodes=nt_bcs,
        )
        if samp["excluded"]:
            raise ValueError(f"{line_name} NTC pool excluded: {samp['exclusion_reason']}")

        resamp = (
            sample_nonoverlapping_resample_sets_r16(
                primary_sampling_audit=samp,
                n_cells=int(nt_raw_idx.size),
                batch_labels=nt_batches,
                n_ctrl=n_ctrl,
                cell_set_len=cell_set_len,
                min_ctrl_floor=min_ctrl_floor,
                resample_seed=resample_seed,
                cell_barcodes=nt_bcs,
            )
            if return_resample
            else None
        )

        x_node = f["X"]
        n_genes = len(g_names_clean)

        def _slice_nadig(sel_local: np.ndarray, audit_dict: Dict[str, Any]) -> SampledRawControlPool:
            sel_raw_rows = nt_raw_idx[sel_local]
            ord_pos = np.argsort(sel_raw_rows, kind="mergesort")
            sorted_raw = sel_raw_rows[ord_pos]
            mat_sorted = np.empty((len(sel_raw_rows), n_genes), dtype=np.float32)

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
            row_sum_umis = raw_counts.sum(axis=1, dtype=np.float64)
            sel_umis = row_sum_umis if nt_umis is None else nt_umis[sel_local]
            return SampledRawControlPool(
                line_name=line_name,
                raw_counts=raw_counts,
                gene_names=g_names_clean,
                gene_ids=g_ids_clean,
                batch_labels=nt_batches[sel_local].astype(str),
                cell_barcodes=nt_bcs[sel_local].astype(str),
                cell_total_umis=np.asarray(sel_umis, dtype=np.float64),
                median_ntc_umi=_median_positive(sel_umis),
                sampling_audit=audit_dict,
                is_log1p=False,
                median_ntc_umi_all_valid=_median_positive(nt_umis) if nt_umis is not None else None,
            )

        pri_pool = _slice_nadig(samp["selected_indices"], samp)
        if not return_resample:
            return pri_pool
        sec_pool = None if (resamp is None or resamp["excluded"]) else _slice_nadig(resamp["selected_indices"], resamp)
        return pri_pool, sec_pool


def _extract_a549_raw_ntc_sample(
    raw_dir: Path,
    exclude_zero_cas9: bool = False,
    min_pct_mt: float = 20.0,
    min_ngenes: int = 500,
    n_ctrl: int = N_CTRL_DEFAULT,
    cell_set_len: int = CELL_SET_LEN_DEFAULT,
    min_ctrl_floor: int = MIN_CTRL_FLOOR_DEFAULT,
    seed: int = DEFAULT_R16_SEED,
    return_resample: bool = False,
    resample_seed: int = RESAMPLE_R16_SEED,
) -> SampledRawControlPool | Tuple[SampledRawControlPool, Optional[SampledRawControlPool]]:
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
    c_id_col = _pick_column(cell_df.columns, ["cell_id", "cell_barcode", "barcode", "Barcode", "Cell_ID", "CBC"])
    c_raw_bc_col = _pick_column(cell_df.columns, ["cell_barcode", "barcode", "Barcode", "CBC"])
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

    cids_raw = cell_df[c_id_col].astype(str).str.strip().tolist() if c_id_col else [f"a549_row_{i}" for i in range(len(cell_df))]
    raw_bcs_col = cell_df[c_raw_bc_col].astype(str).str.strip().tolist() if c_raw_bc_col else cids_raw
    batches_raw = cell_df[c_batch_col].astype(str).str.strip().tolist()
    units_raw = cell_df[c_unit_col].astype(str).str.strip().tolist() if c_unit_col else batches_raw
    sgrnas_raw = cell_df[c_sgrna_col].astype(str).str.strip().tolist()
    perts_raw = cell_df[c_pert_col].astype(str).str.strip().tolist()
    conds_raw = cell_df[c_cond_col].astype(str).str.strip().tolist() if c_cond_col else None
    cas9_vals = pd.to_numeric(cell_df[c_cas9_col], errors="coerce").to_numpy(dtype=np.float64) if c_cas9_col else None

    raw_to_ntc = np.full(len(cell_df), -1, dtype=np.int32)
    ntc_units: List[str] = []
    ntc_cids: List[str] = []
    ntc_raw_bcs: List[str] = []
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
        ntc_cids.append(cids_raw[ci])
        ntc_raw_bcs.append(raw_bcs_col[ci])

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
    bcs_qc = np.asarray([ntc_cids[int(i)] for i in keep_rows], dtype=object)
    raw_bcs_qc = [ntc_raw_bcs[int(i)] for i in keep_rows]

    samp = sample_stratified_control_sets_r16(
        n_cells=int(len(keep_rows)),
        batch_labels=units_qc,
        n_ctrl=n_ctrl,
        cell_set_len=cell_set_len,
        min_ctrl_floor=min_ctrl_floor,
        seed=seed,
        cell_barcodes=bcs_qc,
    )
    if samp["excluded"]:
        raise ValueError(f"a549_ko NTC pool excluded: {samp['exclusion_reason']}")
    samp["selected_raw_cell_barcodes"] = [raw_bcs_qc[int(i)] for i in samp["selected_indices"]]

    def _build_a549_pool(audit_dict: Dict[str, Any]) -> SampledRawControlPool:
        sel = audit_dict["selected_indices"]
        raw_counts = ntc_csr_qc[sel].toarray().astype(np.float32)
        sel_umis = umi_qc[sel]
        return SampledRawControlPool(
            line_name="a549_ko",
            raw_counts=raw_counts,
            gene_names=rna_names,
            gene_ids=rna_ids,
            batch_labels=units_qc[sel].astype(str),
            cell_barcodes=bcs_qc[sel].astype(str),
            cell_total_umis=sel_umis,
            median_ntc_umi=_median_positive(sel_umis),
            sampling_audit=audit_dict,
            is_log1p=False,
            median_ntc_umi_all_valid=_median_positive(umi_qc),
        )

    pri_pool = _build_a549_pool(samp)
    if not return_resample:
        return pri_pool
    resamp = sample_nonoverlapping_resample_sets_r16(
        primary_sampling_audit=samp,
        n_cells=int(len(keep_rows)),
        batch_labels=units_qc,
        n_ctrl=n_ctrl,
        cell_set_len=cell_set_len,
        min_ctrl_floor=min_ctrl_floor,
        resample_seed=resample_seed,
        cell_barcodes=bcs_qc,
    )
    if resamp["excluded"]:
        return pri_pool, None
    resamp["selected_raw_cell_barcodes"] = [raw_bcs_qc[int(i)] for i in resamp["selected_indices"]]
    return pri_pool, _build_a549_pool(resamp)


def _extract_dld1_raw_ntc_sample(
    raw_dir: Path,
    n_ctrl: int = N_CTRL_DEFAULT,
    cell_set_len: int = CELL_SET_LEN_DEFAULT,
    min_ctrl_floor: int = MIN_CTRL_FLOOR_DEFAULT,
    seed: int = DEFAULT_R16_SEED,
    return_resample: bool = False,
    resample_seed: int = RESAMPLE_R16_SEED,
) -> SampledRawControlPool | Tuple[SampledRawControlPool, Optional[SampledRawControlPool]]:
    """Scan all 24 DLD-1 channels for QC-valid NTCs, sample `16 × 32`, then slice only selected cols."""
    ddir = raw_dir / "dld1"
    rna_files = sorted(ddir.glob("GSM98642*_LIB*_SAM24482910.filtered_feature_bc_matrix.h5"))
    hash_files = sorted(ddir.glob("GSM98642*_LIB*_SAM24482910.hashing.umi_correct.h5"))
    crp_files = sorted(ddir.glob("GSM98642*_LIB*_SAM24482910.crispr.umi_correct.h5"))
    if len(rna_files) != 24 or len(hash_files) != 24 or len(crp_files) != 24:
        raise RuntimeError(f"Expected 24 DLD-1 HDF5 trios in {ddir}")

    # (r_path, raw_col, umi, ch_label, unique_bc)
    ntc_records: List[Tuple[Path, int, float, str, str]] = []
    ref_names: Optional[List[str]] = None
    ref_ids: Optional[List[str]] = None

    for ch_idx, (r_path, h_path, c_path) in enumerate(zip(rna_files, hash_files, crp_files), start=1):
        gsm = r_path.name.split("_")[0].upper()
        with h5py.File(r_path, "r") as fr:
            rm = fr["matrix"]
            raw_barcodes = [decode_str_array([x])[0] for x in rm["barcodes"][:]]
            rna_bcs = [x.split("-")[0] for x in raw_barcodes]
            var_ids = [x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x) for x in rm["features/id"][:]]
            var_names = [x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x) for x in rm["features/name"][:]]
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
                ntc_records.append((r_path, int(r_i), float(rna_umi[r_i]), ch_label, f"{gsm}:{raw_barcodes[r_i]}"))

    assert ref_names is not None and ref_ids is not None
    all_umis = np.asarray([r[2] for r in ntc_records], dtype=np.float64)
    all_batches = np.asarray([r[3] for r in ntc_records], dtype=object)
    all_bcs = np.asarray([r[4] for r in ntc_records], dtype=object)

    samp = sample_stratified_control_sets_r16(
        n_cells=len(ntc_records),
        batch_labels=all_batches,
        n_ctrl=n_ctrl,
        cell_set_len=cell_set_len,
        min_ctrl_floor=min_ctrl_floor,
        seed=seed,
        cell_barcodes=all_bcs,
    )
    if samp["excluded"]:
        raise ValueError(f"dld1 NTC pool excluded: {samp['exclusion_reason']}")

    n_genes = len(ref_names)

    def _slice_dld1(audit_dict: Dict[str, Any]) -> SampledRawControlPool:
        sel = audit_dict["selected_indices"]
        by_file: Dict[Path, List[Tuple[int, int]]] = defaultdict(list)
        for pos_out, rec_i in enumerate(sel.tolist()):
            r_path, raw_col, _, _, _ = ntc_records[rec_i]
            by_file[r_path].append((pos_out, raw_col))

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
        sel_umis = all_umis[sel]
        return SampledRawControlPool(
            line_name="dld1",
            raw_counts=raw_counts,
            gene_names=ref_names,
            gene_ids=ref_ids,
            batch_labels=all_batches[sel].astype(str),
            cell_barcodes=all_bcs[sel].astype(str),
            cell_total_umis=sel_umis,
            median_ntc_umi=_median_positive(sel_umis),
            sampling_audit=audit_dict,
            is_log1p=False,
            median_ntc_umi_all_valid=_median_positive(all_umis),
        )

    pri_pool = _slice_dld1(samp)
    if not return_resample:
        return pri_pool
    resamp = sample_nonoverlapping_resample_sets_r16(
        primary_sampling_audit=samp,
        n_cells=len(ntc_records),
        batch_labels=all_batches,
        n_ctrl=n_ctrl,
        cell_set_len=cell_set_len,
        min_ctrl_floor=min_ctrl_floor,
        resample_seed=resample_seed,
        cell_barcodes=all_bcs,
    )
    if resamp["excluded"]:
        return pri_pool, None
    return pri_pool, _slice_dld1(resamp)


def _match_10x_part_selected_barcodes(
    part: PartSpec,
    rows_mat: np.ndarray,
    expected_umis: np.ndarray,
) -> List[str]:
    """Recover exact `f"{donor_id}_{lane_tag}:{barcode}"` strings (`Addendum A §7`) from `part.path`."""
    if rows_mat.shape[0] == 0:
        return []
    prefix = str(part.batches[0]).strip() if len(part.batches) > 0 else part.tag
    if not part.path.exists():
        return [f"{prefix}:cell_{k}" for k in range(rows_mat.shape[0])]

    with h5py.File(part.path, "r") as f:
        if "matrix" not in f or "barcodes" not in f["matrix"]:
            return [f"{prefix}:cell_{k}" for k in range(rows_mat.shape[0])]
        m = f["matrix"]
        barcodes = decode_str_array(m["barcodes"][:])
        ftypes = np.asarray(decode_str_array(m["features/feature_type"][:]), dtype=object)
        rna_rows = np.flatnonzero(ftypes == "Gene Expression")
        indptr = m["indptr"][:]
        data = m["data"]
        indices = m["indices"]

        rna_pos = np.full(int(m["shape"][0]), -1, dtype=np.int32)
        rna_pos[rna_rows] = np.arange(len(rna_rows), dtype=np.int32)
        n_col = len(barcodes)

        # Quick scan of column RNA sums to narrow candidate columns cleanly
        col_rna_umi = np.zeros(n_col, dtype=np.float64)
        for c in range(n_col):
            p0, p1 = int(indptr[c]), int(indptr[c + 1])
            if p1 > p0:
                idx_c = indices[p0:p1]
                dat_c = data[p0:p1]
                mr = rna_pos[idx_c] >= 0
                if np.any(mr):
                    col_rna_umi[c] = float(np.sum(dat_c[mr], dtype=np.float64))

        out_bcs: List[str] = []
        used_cols: Set[int] = set()
        for r_k in range(rows_mat.shape[0]):
            u_exp = float(expected_umis[r_k])
            cand_cols = np.flatnonzero(np.abs(col_rna_umi - u_exp) <= 1e-3)
            target_nz = np.flatnonzero(rows_mat[r_k] != 0.0)
            target_vals = rows_mat[r_k, target_nz]
            found_col: Optional[int] = None
            for col_c in cand_cols:
                if int(col_c) in used_cols:
                    continue
                p0, p1 = int(indptr[int(col_c)]), int(indptr[int(col_c) + 1])
                idx_c = indices[p0:p1]
                dat_c = data[p0:p1]
                rp = rna_pos[idx_c]
                mr = rp >= 0
                rp_rna = rp[mr]
                dv_rna = np.asarray(dat_c[mr], dtype=np.float32)
                nz_sub = dv_rna != 0.0
                rp_rna = rp_rna[nz_sub]
                dv_rna = dv_rna[nz_sub]
                if rp_rna.shape == target_nz.shape and np.array_equal(rp_rna, target_nz) and np.allclose(dv_rna, target_vals, atol=1e-5):
                    found_col = int(col_c)
                    break
            if found_col is not None:
                used_cols.add(found_col)
                out_bcs.append(f"{prefix}:{barcodes[found_col]}")
            else:
                out_bcs.append(f"{prefix}:unmatched_{r_k}")
    return out_bcs


def _extract_cd4t_rest_raw_ntc_sample(
    raw_dir: Path,
    hgnc_bridge: Optional[HGNCBridge] = None,
    n_ctrl: int = N_CTRL_DEFAULT,
    cell_set_len: int = CELL_SET_LEN_DEFAULT,
    min_ctrl_floor: int = MIN_CTRL_FLOOR_DEFAULT,
    seed: int = DEFAULT_R16_SEED,
    return_resample: bool = False,
    resample_seed: int = RESAMPLE_R16_SEED,
) -> SampledRawControlPool | Tuple[SampledRawControlPool, Optional[SampledRawControlPool]]:
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

    ref_names = [str(x) for x in parts[0].var_gene_names]
    ref_ids = [str(x) for x in parts[0].var_gene_ids]

    def _slice_cd4t(audit_dict: Dict[str, Any]) -> SampledRawControlPool:
        sel = audit_dict["selected_indices"]
        by_part: Dict[int, List[Tuple[int, int, float]]] = defaultdict(list)
        for pos_out, rec_i in enumerate(sel.tolist()):
            p_idx, loc_i, _, umi_v = nt_recs[rec_i]
            by_part[p_idx].append((pos_out, loc_i, umi_v))

        raw_counts = np.zeros((len(sel), len(ref_names)), dtype=np.float32)
        barcodes_out = [""] * len(sel)

        for p_idx, items in by_part.items():
            p = parts[p_idx]
            items_sorted = sorted(items, key=lambda kv: kv[1])
            loc_arr = np.asarray([kv[1] for kv in items_sorted], dtype=np.int64)
            out_arr = np.asarray([kv[0] for kv in items_sorted], dtype=np.int64)
            umi_arr = np.asarray([kv[2] for kv in items_sorted], dtype=np.float64)
            assert p.read_rows_fn is not None
            rows_mat = np.asarray(p.read_rows_fn(loc_arr), dtype=np.float32)
            matched_bcs = _match_10x_part_selected_barcodes(p, rows_mat, umi_arr)
            for pos_i, bc_str in zip(out_arr.tolist(), matched_bcs):
                barcodes_out[pos_i] = bc_str
            if [str(x) for x in p.var_gene_names] == ref_names:
                raw_counts[out_arr] = rows_mat
            else:
                pm = {str(g): j for j, g in enumerate(p.var_gene_names)}
                cidx = [pm[g] for g in ref_names if g in pm]
                oidx = [i for i, g in enumerate(ref_names) if g in pm]
                raw_counts[np.ix_(out_arr, oidx)] = rows_mat[:, cidx]

        audit_dict["selected_barcodes"] = barcodes_out
        sel_umis = all_umis[sel]
        return SampledRawControlPool(
            line_name="cd4t_rest",
            raw_counts=raw_counts,
            gene_names=ref_names,
            gene_ids=ref_ids,
            batch_labels=all_batches[sel].astype(str),
            cell_barcodes=np.asarray(barcodes_out, dtype=str),
            cell_total_umis=sel_umis,
            median_ntc_umi=_median_positive(sel_umis),
            sampling_audit=audit_dict,
            is_log1p=False,
            median_ntc_umi_all_valid=_median_positive(all_umis),
        )

    pri_pool = _slice_cd4t(samp)
    if not return_resample:
        return pri_pool
    resamp = sample_nonoverlapping_resample_sets_r16(
        primary_sampling_audit=samp,
        n_cells=len(nt_recs),
        batch_labels=all_batches,
        n_ctrl=n_ctrl,
        cell_set_len=cell_set_len,
        min_ctrl_floor=min_ctrl_floor,
        resample_seed=resample_seed,
    )
    if resamp["excluded"]:
        return pri_pool, None
    return pri_pool, _slice_cd4t(resamp)


def _extract_feng_ipsc_raw_ntc_sample(
    raw_dir: Path,
    n_ctrl: int = N_CTRL_DEFAULT,
    cell_set_len: int = CELL_SET_LEN_DEFAULT,
    min_ctrl_floor: int = MIN_CTRL_FLOOR_DEFAULT,
    seed: int = DEFAULT_R16_SEED,
    return_resample: bool = False,
    resample_seed: int = RESAMPLE_R16_SEED,
) -> SampledRawControlPool | Tuple[SampledRawControlPool, Optional[SampledRawControlPool]]:
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
        keep_cids: List[str] = []
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
            keep_cids.append(cid)
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
        cids_arr = np.asarray(keep_cids, dtype=object)
        samp = sample_stratified_control_sets_r16(
            n_cells=n_ntc,
            batch_labels=batch_arr,
            n_ctrl=n_ctrl,
            cell_set_len=cell_set_len,
            min_ctrl_floor=min_ctrl_floor,
            seed=seed,
            cell_barcodes=cids_arr,
        )
        if samp["excluded"]:
            raise ValueError(f"feng_ipsc NTC pool excluded: {samp['exclusion_reason']}")

        resamp = (
            sample_nonoverlapping_resample_sets_r16(
                primary_sampling_audit=samp,
                n_cells=n_ntc,
                batch_labels=batch_arr,
                n_ctrl=n_ctrl,
                cell_set_len=cell_set_len,
                min_ctrl_floor=min_ctrl_floor,
                resample_seed=resample_seed,
                cell_barcodes=cids_arr,
            )
            if return_resample
            else None
        )

        sel_idx = samp["selected_indices"]
        sec_idx = None if (resamp is None or resamp["excluded"]) else resamp["selected_indices"]
        all_csv_cols = np.asarray(keep_csv_cols, dtype=np.int64)

        umi_all = np.zeros(n_ntc, dtype=np.float64)
        gene_names: List[str] = []
        gene_ids: List[str] = []
        sel_rows_pri: List[np.ndarray] = []
        sel_rows_sec: List[np.ndarray] = []

        def _consume(row_vals: List[str]) -> None:
            if not row_vals:
                return
            eid, sym, is_gex = _parse_feng_gene_row_label(row_vals[0])
            if not is_gex or not sym:
                return
            all_v = np.fromiter((float(row_vals[j]) for j in all_csv_cols), dtype=np.float64, count=n_ntc)
            umi_all[:] += all_v
            sel_rows_pri.append(all_v[sel_idx].astype(np.float32))
            if sec_idx is not None:
                sel_rows_sec.append(all_v[sec_idx].astype(np.float32))
            gene_names.append(sym)
            gene_ids.append(eid or sym)

        _consume(first_row)
        for r in reader:
            _consume(r)

    def _build_feng_pool(rows_list: List[np.ndarray], idx_arr: np.ndarray, audit_dict: Dict[str, Any]) -> SampledRawControlPool:
        raw_counts = np.stack(rows_list, axis=1).astype(np.float32)
        sel_umis = umi_all[idx_arr]
        return SampledRawControlPool(
            line_name="feng_ipsc",
            raw_counts=raw_counts,
            gene_names=gene_names,
            gene_ids=gene_ids,
            batch_labels=batch_arr[idx_arr].astype(str),
            cell_barcodes=cids_arr[idx_arr].astype(str),
            cell_total_umis=sel_umis,
            median_ntc_umi=_median_positive(sel_umis),
            sampling_audit=audit_dict,
            is_log1p=False,
            median_ntc_umi_all_valid=_median_positive(umi_all),
        )

    pri_pool = _build_feng_pool(sel_rows_pri, sel_idx, samp)
    if not return_resample:
        return pri_pool
    sec_pool = None if sec_idx is None else _build_feng_pool(sel_rows_sec, sec_idx, resamp)  # type: ignore[arg-type]
    return pri_pool, sec_pool


# ==============================================================================
# Group 3 Prospective Holdout (`GSE295214`, H1 hESC, `PREREG_I v5 §2.4`)
# ==============================================================================


def _extract_gsm_channel_id(path: Path) -> str:
    m = re.match(r"^(GSM\d+)", path.name, flags=re.IGNORECASE)
    return m.group(1).upper() if m else path.stem


def _find_gse295214_h5_files(raw_dir: Path, allow_partial_channels: bool = False) -> List[Tuple[str, Path]]:
    """Locate and sort `GSE295214` 10x CellRanger `.h5` files in lexicographical GSM order."""
    gdir = Path(raw_dir)
    if (gdir / "gse295214_h1").exists():
        gdir = gdir / "gse295214_h1"
    h5s = sorted(p for p in gdir.glob("*.h5") if p.is_file())
    if not h5s:
        raise FileNotFoundError(f"No .h5 files found in GSE295214 dir: {gdir}")
    paired = sorted(((_extract_gsm_channel_id(p), p) for p in h5s), key=lambda kv: (kv[0], kv[1].name))
    if not allow_partial_channels and len(paired) != 85:
        raise RuntimeError(f"PREREG_I §2.4 requires 85 GSM .h5 channels for GSE295214, found {len(paired)} in {gdir}")
    return paired


def _is_gse295214_ntc(guide_name: str, guide_id: str = "") -> bool:
    """Strictly match `non-targeting_*` single or pipe-delimited dual NTCs (`Addendum A §2.1`)."""
    gn = str(guide_name).strip()
    gi = str(guide_id).strip()
    if gn.lower().startswith("non-targeting_") or gi.lower().startswith("non-targeting_"):
        return True
    raw = gn or gi
    if "|" in raw:
        subs = [s.strip() for s in raw.split("|") if s.strip()]
        if subs and all(s.lower().startswith("non-targeting_") for s in subs):
            return True
    return False


def _resolve_single_guide_token(
    token: str,
    hgnc: HGNCLookup,
) -> Optional[Tuple[str, str, str, str]]:
    """
    Exact three-stage deterministic guide-token resolution (`PREREG_I_ADDENDUM_A §2.1`):
      1. Exact full-token `_try_resolve(tok, ens_raw)` (protects `HLA-A`, `MT-CO1`)
      2. Underscore `_` delimited multi-part tokens (`len(parts) >= 2`) with multi-gene guard
      3. Trailing hyphen `-` suffix strip (`-(?:\\d+|sg\\d+|P\\d+.*|TSS\\d*)`)
    Accepts ONLY canonical HGNC symbols present in `hgnc.approved`.
    """
    tok = str(token).strip()
    if not tok or _GSE295214_ILLEGAL_TOK_RE.search(tok):
        return None
    low = tok.lower()
    if low.startswith("non-targeting_") or low == "non-targeting":
        return None
    if any(b in low for b in _GSE295214_BLACKLIST_SUBSTRS) or "-and-" in low:
        return None

    m_ens = _GSE295214_ENSG_RE.search(tok)
    ens_raw = m_ens.group(1).upper() if m_ens else ""

    def _try_resolve(sym_cand: str, ens_cand: str) -> Optional[Tuple[str, str, str, str]]:
        sc = str(sym_cand).strip()
        if not sc:
            return None
        canon_sym, hstat = hgnc.resolve(sc, str(ens_cand).strip())
        if canon_sym not in hgnc.approved:
            return None
        ens_out = str(ens_cand).strip() or hgnc.sym_to_ens.get(canon_sym, "")
        return ens_out, canon_sym, sc, hstat

    # Stage 1: Exact full-token lookup
    hit1 = _try_resolve(tok, ens_raw)
    if hit1 is not None:
        return hit1

    # Stage 2: Underscore `_` delimited tokens (`len(parts) >= 2`)
    if "_" in tok:
        parts = [p for p in tok.split("_") if p]
        if len(parts) >= 2:
            matched_syms: Set[str] = set()
            for seg in parts:
                if _GSE295214_IGNORE_SUFFIX_SEG_RE.match(seg):
                    continue
                seg_hit = _try_resolve(seg, "")
                if seg_hit is not None:
                    matched_syms.add(seg_hit[1])
            if len(matched_syms) > 1:
                return None

            if parts[0].upper().startswith("ENSG"):
                hit2a = _try_resolve(parts[1], parts[0].split(".")[0])
                if hit2a is not None:
                    return hit2a
            hit2b = _try_resolve(parts[0], ens_raw)
            if hit2b is not None:
                return hit2b

    # Stage 3: Trailing hyphen `-` suffix strip ONLY if stages 1–2 did not match
    m_hy = _GSE295214_HYPHEN_SUFFIX_RE.match(tok)
    if m_hy:
        hit3 = _try_resolve(m_hy.group(1), ens_raw)
        if hit3 is not None:
            return hit3

    return None


def parse_gse295214_guide_r16(
    guide_name: str,
    guide_id: str,
    hgnc: HGNCLookup,
) -> Tuple[str, Optional[Tuple[str, str, str, str]]]:
    """
    Classify a GSE295214 non-GEX feature strictly per `PREREG_I_ADDENDUM_A §2.1` into
    `("ntc", None)`, `("target", (ens, canon_sym, raw_sym, hstat))`, or `("exclude", None)`.
    """
    gnm = str(guide_name).strip()
    gid = str(guide_id).strip()
    if gid.lower().startswith("non-targeting_") or gnm.lower().startswith("non-targeting_"):
        return "ntc", None

    raw = gnm or gid
    if not raw:
        return "exclude", None

    if "|" in raw:
        subs = [s.strip() for s in raw.split("|") if s.strip()]
        if not subs:
            return "exclude", None
        if all(s.lower().startswith("non-targeting_") for s in subs):
            return "ntc", None
        resolved_subs: List[Tuple[str, str, str, str]] = []
        for sub in subs:
            r_sub = _resolve_single_guide_token(sub, hgnc)
            if r_sub is None:
                return "exclude", None
            resolved_subs.append(r_sub)
        if len({r[1] for r in resolved_subs}) != 1:
            return "exclude", None
        return "target", resolved_subs[0]

    hit = _resolve_single_guide_token(raw, hgnc)
    if hit is None and gid and gid != raw:
        hit = _resolve_single_guide_token(gid, hgnc)
    if hit is None:
        return "exclude", None
    return "target", hit


def _parse_gse295214_target(
    guide_name: str,
    guide_id: str,
    hgnc: HGNCLookup,
) -> Optional[Tuple[str, str, str, str]]:
    """Backward-compatible helper returning `(ens, canon_sym, raw_sym, hstat)` or `None`."""
    kind, info = parse_gse295214_guide_r16(guide_name, guide_id, hgnc)
    return info if kind == "target" else None


def _validate_10x_rna_integers(h5_path: Path, rna_row_idx: np.ndarray) -> None:
    """Verify 100% of stored non-zero RNA counts are finite, non-negative exact integers (`§2.4.2 L132`)."""
    with h5py.File(h5_path, "r") as f:
        m = f["matrix"]
        n_feat = int(m["shape"][0])
        is_rna = np.zeros(n_feat, dtype=bool)
        is_rna[rna_row_idx] = True
        idx = m["indices"][:]
        dat = m["data"][:]
        r_dat = np.asarray(dat[is_rna[idx]], dtype=np.float64)
        if r_dat.size == 0:
            return
        if not np.all(np.isfinite(r_dat) & (r_dat >= 0.0) & (r_dat == np.floor(r_dat))):
            raise ValueError(f"Non-integer or negative RNA counts detected in {h5_path}")


def _scan_gse295214_channels_common(
    raw_dir: Path,
    hgnc: HGNCLookup,
    allow_partial_channels: bool = False,
) -> Tuple[List[Dict[str, Any]], np.ndarray, Dict[str, Any]]:
    """
    Shared Pass-1 + Pass-2 single-guide QC scanner strictly adhering to `PREREG_I v5 §2.4.2`:
      - Checks exact non-negative integers on `Gene Expression`
      - Pooled `fit_poisson_gaussian_thresholds_from_hist(hist, min_thr=5, max_thr=20, default_thr=7)`
        over cells with `n_genes >= 200`
      - Single-guide `n_pos == 1 & top1 >= 5 & (top2 < 3 or top1 / max(top2, 0.5) >= 3.0)`
      - Cell QC `n_umis >= 300 & n_genes >= 200`, `max_pct_mt = None`
      - Lexicographical odd/even split-half across channels + within-half `< 15` NTC fallback
    """
    ch_files = _find_gse295214_h5_files(raw_dir, allow_partial_channels=allow_partial_channels)
    hist: Optional[np.ndarray] = None
    crp_names_ref: Optional[List[str]] = None
    crp_ids_ref: Optional[List[str]] = None
    channel_cache: List[Dict[str, Any]] = []

    for ch_pos_1based, (gsm, h5_path) in enumerate(ch_files, start=1):
        half_val = 1 if (ch_pos_1based % 2 == 1) else 2
        with h5py.File(h5_path, "r") as f:
            m = f["matrix"]
            n_feat, n_col = int(m["shape"][0]), int(m["shape"][1])
            barcodes = [str(x) for x in decode_str_array(m["barcodes"][:])]
            ftypes = np.asarray(decode_str_array(m["features/feature_type"][:]), dtype=object)
            fids = [str(x) for x in decode_str_array(m["features/id"][:])]
            fnames = [str(x) for x in decode_str_array(m["features/name"][:])]

            rna_rows = np.flatnonzero(ftypes == "Gene Expression")
            crp_rows = np.flatnonzero(ftypes == "CRISPR Guide Capture")
            if crp_rows.size == 0:
                crp_rows = np.flatnonzero(ftypes != "Gene Expression")
            if rna_rows.size == 0 or crp_rows.size == 0:
                raise RuntimeError(f"Expected both Gene Expression & CRISPR Guide Capture in {h5_path}")

            if crp_names_ref is None:
                crp_names_ref = [fnames[int(i)] for i in crp_rows]
                crp_ids_ref = [fids[int(i)] for i in crp_rows]
                hist = np.zeros((len(crp_rows), 66), dtype=np.int64)

            _validate_10x_rna_integers(h5_path, rna_rows)
            mt_raw = np.zeros(n_feat, dtype=bool)
            mt_raw[rna_rows] = [fnames[int(i)].startswith("MT-") for i in rna_rows]
            umi_l, ng_l, _ = compute_10x_csc_rna_qc(h5_path, rna_rows, mt_raw)

            ok_hist = ng_l >= 200
            ok_qc = (umi_l >= 300.0) & (ng_l >= 200)

            indptr = m["indptr"][:]
            dat_all = m["data"][:]
            idx_all = m["indices"][:]
            crp_map = np.full(n_feat, -1, dtype=np.int32)
            crp_map[crp_rows] = np.arange(len(crp_rows), dtype=np.int32)

            # Accumulate histogram strictly on `n_genes >= 200` cells and extract CRISPR matrix for `ok_qc`
            col_lens = np.diff(indptr)
            col_ids = np.repeat(np.arange(n_col, dtype=np.int32), col_lens)
            m_c = crp_map[idx_all] >= 0
            c_idx = crp_map[idx_all[m_c]]
            c_dat = dat_all[m_c].astype(np.int64, copy=False)
            c_col = col_ids[m_c]

            hist_mask = ok_hist[c_col] & (c_dat >= 1)
            if np.any(hist_mask):
                np.add.at(hist, (c_idx[hist_mask], np.clip(c_dat[hist_mask], 1, 65)), 1)

            qc_cols = np.flatnonzero(ok_qc)
            col_map_qc = np.full(n_col, -1, dtype=np.int32)
            col_map_qc[qc_cols] = np.arange(len(qc_cols), dtype=np.int32)
            m_qc_c = col_map_qc[c_col] >= 0
            C_qc = sp.csr_matrix(
                (c_dat[m_qc_c].astype(np.float64), (c_idx[m_qc_c], col_map_qc[c_col[m_qc_c]])),
                shape=(len(crp_rows), len(qc_cols)),
            ).tocsc()

        channel_cache.append(
            {
                "gsm": gsm,
                "half_val": half_val,
                "path": h5_path,
                "n_feat": n_feat,
                "barcodes": barcodes,
                "rna_rows": rna_rows,
                "rids": [fids[int(i)] for i in rna_rows],
                "rnames": [fnames[int(i)] for i in rna_rows],
                "crp_names": [fnames[int(i)] for i in crp_rows],
                "crp_ids": [fids[int(i)] for i in crp_rows],
                "umi_l": umi_l,
                "qc_cols": qc_cols,
                "C_qc": C_qc,
            }
        )

    assert hist is not None and crp_names_ref is not None and crp_ids_ref is not None
    thrs = fit_poisson_gaussian_thresholds_from_hist(hist, min_thr=5, max_thr=20, default_thr=7)

    # Identify NTC guides (`non-targeting_*`) and target guides strictly per Addendum A §2.1
    ntc_guide_indices: Set[int] = set()
    guide_target_map: Dict[int, Tuple[str, str, str, str]] = {}
    for gi, (gn, gid) in enumerate(zip(crp_names_ref, crp_ids_ref)):
        kind, info = parse_gse295214_guide_r16(gn, gid, hgnc)
        if kind == "ntc":
            ntc_guide_indices.add(gi)
        elif kind == "target" and info is not None:
            guide_target_map[gi] = info

    # Per-channel Pass 2: single-guide call + raw column extraction + within-half `< 15` NTC check
    all_nt_flags: List[bool] = []
    all_batches: List[str] = []
    all_halves: List[int] = []
    ch_results: List[Dict[str, Any]] = []
    tot_ntc = tot_targ = 0

    for rec in channel_cache:
        C_qc = rec["C_qc"]
        t1, t2, a1 = extract_top2_csc(C_qc)
        pos_mask = C_qc.data >= thrs[C_qc.indices]
        pos_csc = sp.csc_matrix((pos_mask.astype(np.int8), C_qc.indices, C_qc.indptr), shape=C_qc.shape)
        n_pos = np.asarray(pos_csc.sum(axis=0)).ravel()
        s1_mask = (n_pos == 1) & (t1 >= 5.0) & ((t2 < 3.0) | (t1 / np.maximum(t2, 0.5) >= 3.0))

        keep_cols: List[int] = []
        keep_bcs: List[str] = []
        nt_flags: List[bool] = []
        pkeys: List[str] = []
        pmeta: Dict[str, Tuple[str, str, str, str]] = {}
        gsm = str(rec["gsm"]).upper()
        hv = int(rec["half_val"])
        qc_cols = rec["qc_cols"]
        bcs = rec["barcodes"]
        cn_local = rec["crp_names"]
        cid_local = rec["crp_ids"]

        for loc_j in np.flatnonzero(s1_mask):
            gi = int(a1[loc_j])
            raw_col = int(qc_cols[loc_j])
            if gi in ntc_guide_indices or _is_gse295214_ntc(cn_local[gi], cid_local[gi]):
                keep_cols.append(raw_col)
                keep_bcs.append(f"{gsm}:{bcs[raw_col]}")
                nt_flags.append(True)
                pkeys.append("")
                tot_ntc += 1
            elif gi in guide_target_map:
                gkey = str(cn_local[gi]).strip() or str(cid_local[gi]).strip()
                keep_cols.append(raw_col)
                keep_bcs.append(f"{gsm}:{bcs[raw_col]}")
                nt_flags.append(False)
                pkeys.append(gkey)
                if gkey not in pmeta:
                    pmeta[gkey] = guide_target_map[gi]
                tot_targ += 1

        all_nt_flags.extend(nt_flags)
        all_batches.extend([gsm] * len(keep_cols))
        all_halves.extend([hv] * len(keep_cols))
        ch_results.append(
            {
                **rec,
                "keep_cols": np.asarray(keep_cols, dtype=np.int64),
                "keep_bcs": keep_bcs,
                "nt_flags": np.asarray(nt_flags, dtype=bool),
                "pkeys": pkeys,
                "pmeta": pmeta,
            }
        )

    nt_arr_all = np.asarray(all_nt_flags, dtype=bool)
    half_arr_all = np.asarray(all_halves, dtype=np.int8)
    _, n_fb = _apply_min_ntc_batch_fallback(
        np.asarray(all_batches, dtype=object),
        nt_arr_all,
        half_arr_all,
        group_labels=None,
        min_batch_nt=15,
    )
    nt_half_a = int(np.sum(nt_arr_all & (half_arr_all == 1)))
    nt_half_b = int(np.sum(nt_arr_all & (half_arr_all == 2)))

    audit = {
        "n_channels": len(ch_files),
        "n_crp_guides": len(crp_names_ref),
        "n_ntc_guides": len(ntc_guide_indices),
        "n_single_ntc": int(tot_ntc),
        "n_single_targeting": int(tot_targ),
        "n_ntc_half_a": nt_half_a,
        "n_ntc_half_b": nt_half_b,
        "n_low_ntc_fallback_cells": int(n_fb),
        "median_guide_threshold": float(np.median(thrs)),
        "half_rule": "lexicographic_gsm_odd_HalfA_vs_even_HalfB",
        "integer_rna_verified": True,
    }
    return ch_results, thrs, audit


def build_gse295214_parts(
    raw_dir: Path,
    hgnc: Optional[HGNCLookup] = None,
    allow_partial_channels: bool = False,
) -> Tuple[List[PartSpec], Dict[str, Any]]:
    """
    Build standardized `PartSpec` objects across `GSE295214` sorted GSM channels (`PREREG_I v5 §2.4.2`).
    """
    hgnc_l = hgnc or build_hgnc_lookup()
    ch_results, _, audit = _scan_gse295214_channels_common(
        raw_dir=raw_dir, hgnc=hgnc_l, allow_partial_channels=allow_partial_channels
    )
    parts: List[PartSpec] = []
    for rec in ch_results:
        keep_arr = rec["keep_cols"]
        if keep_arr.size == 0:
            continue
        read_chunk, read_rows, read_chunk_cols = make_10x_csc_subset_readers(
            rec["path"], keep_arr, rec["n_feat"], rec["rna_rows"]
        )
        parts.append(
            PartSpec(
                tag=rec["gsm"],
                path=rec["path"],
                n_cells=int(keep_arr.size),
                n_genes_raw=int(len(rec["rna_rows"])),
                nt_mask=rec["nt_flags"],
                batches=np.full(int(keep_arr.size), rec["gsm"], dtype=object),
                half_flags=np.full(int(keep_arr.size), int(rec["half_val"]), dtype=np.int8),
                umi=rec["umi_l"][keep_arr],
                pert_keys=rec["pkeys"],
                pert_meta=rec["pmeta"],
                var_gene_ids=rec["rids"],
                var_gene_names=rec["rnames"],
                read_chunk_fn=read_chunk,
                read_rows_fn=read_rows,
                read_chunk_cols_fn=read_chunk_cols,
            )
        )
    return parts, audit


def prepare_gse295214_standardized_line(
    raw_dir: Path,
    out_root: Path = R16_LINES_DIR,
    workers: int = 8,
    chunk_size: int = 8192,
    allow_partial_channels: bool = False,
) -> Tuple[Path, Dict[str, Any]]:
    """Run frozen `process_r10_line` standardization on `gse295214_h1` cleanly."""
    hgnc_l = build_hgnc_lookup()
    parts, audit = build_gse295214_parts(raw_dir, hgnc=hgnc_l, allow_partial_channels=allow_partial_channels)
    out_dir = process_r10_line(
        line_name="gse295214_h1",
        parts=parts,
        hgnc=hgnc_l,
        chunk_size=chunk_size,
        out_root=out_root,
        extra_meta={"adapter_audit": audit},
        n_workers=workers,
    )
    return out_dir, audit


def extract_gse295214_ntc_raw_counts(
    raw_dir: Path,
    n_ctrl: int = N_CTRL_DEFAULT,
    cell_set_len: int = CELL_SET_LEN_DEFAULT,
    min_ctrl_floor: int = MIN_CTRL_FLOOR_DEFAULT,
    seed: int = DEFAULT_R16_SEED,
    return_resample: bool = False,
    resample_seed: int = RESAMPLE_R16_SEED,
    allow_partial_channels: bool = False,
) -> SampledRawControlPool | Tuple[SampledRawControlPool, Optional[SampledRawControlPool]]:
    """
    Extract ONLY QC-valid `non-targeting_*` single-guide NTC raw counts from `GSE295214`,
    proportionally stratified across sorted GSM channels (`16 × 32`, `seed = 20261601`).
    """
    hgnc_l = build_hgnc_lookup() if DEFAULT_HGNC_FILE.exists() else HGNCLookup(set(), {}, {}, {}, {})
    ch_results, _, ch_audit = _scan_gse295214_channels_common(
        raw_dir=raw_dir, hgnc=hgnc_l, allow_partial_channels=allow_partial_channels
    )

    nt_recs: List[Tuple[int, int, str, float, str]] = []
    for ch_idx, rec in enumerate(ch_results):
        keep_cols = rec["keep_cols"]
        keep_bcs = rec["keep_bcs"]
        nt_flags = rec["nt_flags"]
        umi_l = rec["umi_l"]
        gsm = rec["gsm"]
        for loc_i in np.flatnonzero(nt_flags):
            raw_col = int(keep_cols[loc_i])
            nt_recs.append((ch_idx, raw_col, gsm, float(umi_l[raw_col]), str(keep_bcs[loc_i])))

    if not nt_recs:
        raise RuntimeError("Zero QC-valid non-targeting_* NTC singlets found in GSE295214.")

    ref_names = ch_results[0]["rnames"]
    ref_ids = ch_results[0]["rids"]
    all_batches = np.asarray([r[2] for r in nt_recs], dtype=object)
    all_umis = np.asarray([r[3] for r in nt_recs], dtype=np.float64)
    all_bcs = np.asarray([r[4] for r in nt_recs], dtype=object)

    samp = sample_stratified_control_sets_r16(
        n_cells=len(nt_recs),
        batch_labels=all_batches,
        n_ctrl=n_ctrl,
        cell_set_len=cell_set_len,
        min_ctrl_floor=min_ctrl_floor,
        seed=seed,
        cell_barcodes=all_bcs,
    )
    if samp["excluded"]:
        raise ValueError(f"gse295214_h1 NTC pool excluded: {samp['exclusion_reason']}")
    samp["gse295214_channel_audit"] = ch_audit

    n_genes = len(ref_names)

    def _slice_gse(audit_dict: Dict[str, Any]) -> SampledRawControlPool:
        sel = audit_dict["selected_indices"]
        by_ch: Dict[int, List[Tuple[int, int]]] = defaultdict(list)
        for pos_out, rec_i in enumerate(sel.tolist()):
            ch_idx, raw_col, _, _, _ = nt_recs[rec_i]
            by_ch[ch_idx].append((pos_out, raw_col))

        raw_counts = np.zeros((len(sel), n_genes), dtype=np.float32)
        for ch_idx, items in by_ch.items():
            rec = ch_results[ch_idx]
            rna_rows = rec["rna_rows"]
            rna_pos = np.full(rec["n_feat"], -1, dtype=np.int32)
            if rec["rnames"] == ref_names:
                rna_pos[rna_rows] = np.arange(n_genes, dtype=np.int32)
            else:
                rmap = {g: j for j, g in enumerate(ref_names)}
                for src_ri, nm in zip(rna_rows, rec["rnames"]):
                    if nm in rmap:
                        rna_pos[int(src_ri)] = rmap[nm]
            with h5py.File(rec["path"], "r") as f:
                m = f["matrix"]
                indptr = m["indptr"]
                data = m["data"]
                indices = m["indices"]
                for pos_out, col_i in items:
                    p0, p1 = int(indptr[col_i]), int(indptr[col_i + 1])
                    if p1 > p0:
                        rp = rna_pos[indices[p0:p1]]
                        mr = rp >= 0
                        if np.any(mr):
                            raw_counts[pos_out, rp[mr]] = np.asarray(data[p0:p1][mr], dtype=np.float32)
        sel_umis = all_umis[sel]
        return SampledRawControlPool(
            line_name="gse295214_h1",
            raw_counts=raw_counts,
            gene_names=list(ref_names),
            gene_ids=list(ref_ids),
            batch_labels=all_batches[sel].astype(str),
            cell_barcodes=all_bcs[sel].astype(str),
            cell_total_umis=sel_umis,
            median_ntc_umi=_median_positive(sel_umis),
            sampling_audit=audit_dict,
            is_log1p=False,
            median_ntc_umi_all_valid=_median_positive(all_umis),
        )

    pri_pool = _slice_gse(samp)
    if not return_resample:
        return pri_pool
    resamp = sample_nonoverlapping_resample_sets_r16(
        primary_sampling_audit=samp,
        n_cells=len(nt_recs),
        batch_labels=all_batches,
        n_ctrl=n_ctrl,
        cell_set_len=cell_set_len,
        min_ctrl_floor=min_ctrl_floor,
        resample_seed=resample_seed,
        cell_barcodes=all_bcs,
    )
    if resamp["excluded"]:
        return pri_pool, None
    return pri_pool, _slice_gse(resamp)


def evaluate_gse295214_counts_gate_r16(
    std_line_dir: Path,
    model_genes: Sequence[str],
    model_pert_names: Sequence[str],
    ntc_half_a: int,
    ntc_half_b: int,
    hgnc_bridge: Optional[HGNCBridge] = None,
    out_md_path: Optional[Path] = None,
    src_line_dirs: Optional[Dict[str, Path]] = None,
) -> Dict[str, Any]:
    """
    Counts-only qualifying gate BEFORE any metric calculation (`PREREG_I v5 §2.4.3 L140–144`):
      - `N_eval^R16(gse295214_h1) >= 50` (using `r9_6 \\ {vcc2025_h1}` primary pool)
      - `|G_eval^R9(gse295214_h1)| = |canon(G_R9_6148) ∩ canon(G_fin)| >= 3,000`
      - `per-half NTC >= 15` (`Half A >= 15` and `Half B >= 15`)
    """
    hb = get_hgnc_bridge(hgnc_bridge)
    r9_common_canon = load_r9_common_genes(hgnc_bridge=hb)
    tgt_ds = load_standardized_line_r10(
        Path(std_line_dir), representation="log1p", common_genes=r9_common_canon, hgnc_bridge=hb, min_gene_overlap=100
    )
    g_eval_r9_n = len(tgt_ds.genes)

    dir_map = src_line_dirs or STANDARD_LINE_DIRS
    src_rep_targets: Set[str] = set()
    for sn in GROUP3_PRIMARY_FCT_SOURCES:
        sd = load_standardized_line_r10(
            dir_map[sn], representation="log1p", common_genes=r9_common_canon, hgnc_bridge=hb, min_gene_overlap=100
        )
        if sd.repr_mask is not None:
            src_rep_targets.update(t for t, ok in zip(sd.target_names, sd.repr_mask) if ok)
        else:
            src_rep_targets.update(sd.target_names)

    model_g_canon = {hb.canonicalize(str(g)) for g in model_genes if hb.canonicalize(str(g))}
    model_p_canon = {hb.canonicalize(str(p)) for p in model_pert_names if str(p) != "non-targeting" and hb.canonicalize(str(p))}
    g_eval_r9_set = set(tgt_ds.genes)
    g_eval_r16_set = g_eval_r9_set & model_g_canon

    n_eval = 0
    n_self_masked = 0
    n_half_valid = 0
    for ri, t_sym in enumerate(tgt_ds.target_names):
        if tgt_ds.repr_mask is not None and not bool(tgt_ds.repr_mask[ri]):
            continue
        if t_sym not in src_rep_targets:
            continue
        if t_sym not in model_p_canon:
            continue
        if t_sym not in g_eval_r9_set:
            continue
        n_eval += 1
        if t_sym in g_eval_r16_set:
            n_self_masked += 1
        if tgt_ds.half_valid is None or bool(tgt_ds.half_valid[ri]):
            n_half_valid += 1

    pass_n = n_eval >= 50
    pass_g = g_eval_r9_n >= 3000
    pass_ntc = (int(ntc_half_a) >= 15) and (int(ntc_half_b) >= 15)
    qualified = bool(pass_n and pass_g and pass_ntc)
    status = "QUALIFIED_PRIMARY" if qualified else "SENSITIVITY_ONLY"

    res = {
        "line_name": "gse295214_h1",
        "status": status,
        "qualified_primary": qualified,
        "n_eval_r16": int(n_eval),
        "n_half_valid": int(n_half_valid),
        "n_self_masked": int(n_self_masked),
        "g_eval_r9": int(g_eval_r9_n),
        "g_eval_r16": int(len(g_eval_r16_set)),
        "ntc_half_a": int(ntc_half_a),
        "ntc_half_b": int(ntc_half_b),
        "pass_n_eval_ge_50": pass_n,
        "pass_g_eval_r9_ge_3000": pass_g,
        "pass_per_half_ntc_ge_15": pass_ntc,
        "fct_primary_pool": list(GROUP3_PRIMARY_FCT_SOURCES),
    }

    if out_md_path is not None:
        out_p = Path(out_md_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        md = "\n".join(
            [
                "# Round 16 Group 3 Counts-Only QC Gate (`R16_GSE295214_QC_COUNTS.md`)",
                "",
                f"- **Qualification Status:** `{status}` (`qualified_primary={qualified}`)",
                f"- **`N_eval^R16` (`>= 50`, excl. `vcc2025_h1`):** `{n_eval}` (`half_valid={n_half_valid}`, `self_masked={n_self_masked}`, pass=`{pass_n}`)",
                f"- **`|G_eval^R9|` (`>= 3,000`):** `{g_eval_r9_n}` (`|G_eval^R16|={len(g_eval_r16_set)}`, pass=`{pass_g}`)",
                f"- **Per-Half NTCs (`>= 15` each):** Half A=`{int(ntc_half_a)}`, Half B=`{int(ntc_half_b)}` (pass=`{pass_ntc}`)",
                f"- **Locked Primary `FCT` Pool:** `{', '.join(GROUP3_PRIMARY_FCT_SOURCES)}`",
                "",
            ]
        )
        out_p.write_text(md, encoding="utf-8")
    return res


def extract_raw_control_pool_r16(
    line_name: str,
    custom_nadig_dir: Optional[Path] = None,
    custom_gse295214_dir: Optional[Path] = None,
    n_ctrl: int = N_CTRL_DEFAULT,
    cell_set_len: int = CELL_SET_LEN_DEFAULT,
    min_ctrl_floor: int = MIN_CTRL_FLOOR_DEFAULT,
    seed: int = DEFAULT_R16_SEED,
    allow_npz_fallback_for_hvg: bool = False,
    return_resample: bool = False,
    resample_seed: int = RESAMPLE_R16_SEED,
    allow_partial_gse295214_channels: bool = False,
) -> SampledRawControlPool | Tuple[SampledRawControlPool, Optional[SampledRawControlPool]]:
    """
    Dispatch exact raw NTC single-cell extraction for all Gate-2 local raw cohorts plus
    prospective Group 3 (`gse295214_h1`, `PREREG_I v5 §2.4`).
    """
    if line_name in ("nadig_hepg2", "nadig_jurkat"):
        try:
            return _extract_nadig_raw_ntc_sample(
                line_name=line_name,
                custom_nadig_dir=custom_nadig_dir,
                n_ctrl=n_ctrl,
                cell_set_len=cell_set_len,
                min_ctrl_floor=min_ctrl_floor,
                seed=seed,
                return_resample=return_resample,
                resample_seed=resample_seed,
            )
        except FileNotFoundError:
            if not allow_npz_fallback_for_hvg:
                raise
            ld = STANDARD_LINE_DIRS[line_name]
            genes = [x.strip() for x in (ld / "genes.txt").read_text(encoding="utf-8").splitlines() if x.strip()]
            ctrl_npz = np.load(ld / "ctrl_cells.npz", allow_pickle=False)
            batches = np.asarray([str(x) for x in ctrl_npz["batch_labels"]], dtype=object)
            x_log1p = np.asarray(ctrl_npz["X_log1p"], dtype=np.float32)
            bcs = (
                np.asarray([str(x) for x in ctrl_npz["cell_barcodes"]], dtype=object)
                if "cell_barcodes" in ctrl_npz
                else np.asarray([f"{line_name}_npz_{i}" for i in range(x_log1p.shape[0])], dtype=object)
            )
            samp = sample_stratified_control_sets_r16(
                n_cells=x_log1p.shape[0],
                batch_labels=batches,
                n_ctrl=n_ctrl,
                cell_set_len=cell_set_len,
                min_ctrl_floor=min_ctrl_floor,
                seed=seed,
                cell_barcodes=bcs,
            )

            def _npz_pool(audit_d: Dict[str, Any]) -> SampledRawControlPool:
                sel = audit_d["selected_indices"]
                return SampledRawControlPool(
                    line_name=line_name,
                    raw_counts=x_log1p[sel],
                    gene_names=genes,
                    gene_ids=genes,
                    batch_labels=batches[sel].astype(str),
                    cell_barcodes=bcs[sel].astype(str),
                    cell_total_umis=np.full(len(sel), 1e4, dtype=np.float64),
                    median_ntc_umi=1e4,
                    sampling_audit=audit_d,
                    is_log1p=True,
                )

            pri_p = _npz_pool(samp)
            if not return_resample:
                return pri_p
            resamp = sample_nonoverlapping_resample_sets_r16(
                primary_sampling_audit=samp,
                n_cells=x_log1p.shape[0],
                batch_labels=batches,
                n_ctrl=n_ctrl,
                cell_set_len=cell_set_len,
                min_ctrl_floor=min_ctrl_floor,
                resample_seed=resample_seed,
                cell_barcodes=bcs,
            )
            return pri_p, (None if resamp["excluded"] else _npz_pool(resamp))

    if line_name == "a549_ko":
        return _extract_a549_raw_ntc_sample(
            REPO_ROOT / "round12" / "shared" / "raw",
            n_ctrl=n_ctrl,
            cell_set_len=cell_set_len,
            min_ctrl_floor=min_ctrl_floor,
            seed=seed,
            return_resample=return_resample,
            resample_seed=resample_seed,
        )
    if line_name == "dld1":
        return _extract_dld1_raw_ntc_sample(
            REPO_ROOT / "round10" / "shared" / "raw",
            n_ctrl=n_ctrl,
            cell_set_len=cell_set_len,
            min_ctrl_floor=min_ctrl_floor,
            seed=seed,
            return_resample=return_resample,
            resample_seed=resample_seed,
        )
    if line_name == "cd4t_rest":
        return _extract_cd4t_rest_raw_ntc_sample(
            REPO_ROOT / "round10" / "shared" / "raw",
            n_ctrl=n_ctrl,
            cell_set_len=cell_set_len,
            min_ctrl_floor=min_ctrl_floor,
            seed=seed,
            return_resample=return_resample,
            resample_seed=resample_seed,
        )
    if line_name == "feng_ipsc":
        return _extract_feng_ipsc_raw_ntc_sample(
            REPO_ROOT / "round11" / "shared" / "raw",
            n_ctrl=n_ctrl,
            cell_set_len=cell_set_len,
            min_ctrl_floor=min_ctrl_floor,
            seed=seed,
            return_resample=return_resample,
            resample_seed=resample_seed,
        )
    if line_name == "gse295214_h1":
        gdir = Path(custom_gse295214_dir) if custom_gse295214_dir is not None else (R16_RAW_DIR / "gse295214_h1")
        return extract_gse295214_ntc_raw_counts(
            raw_dir=gdir,
            n_ctrl=n_ctrl,
            cell_set_len=cell_set_len,
            min_ctrl_floor=min_ctrl_floor,
            seed=seed,
            return_resample=return_resample,
            resample_seed=resample_seed,
            allow_partial_channels=allow_partial_gse295214_channels,
        )
    raise ValueError(f"Unsupported control source {line_name!r}")


def _resolve_fct_source_lines(target_line: str, fct_pool_mode: str = "primary") -> List[str]:
    if target_line in ("nadig_hepg2", "nadig_jurkat"):
        return sorted(set(["nadig_hepg2", "nadig_jurkat", "replogle_k562_essential", "replogle_rpe1"]) - {target_line})
    if target_line == "gse295214_h1":
        if fct_pool_mode == "primary":
            return list(GROUP3_PRIMARY_FCT_SOURCES)  # r9_6 \ {vcc2025_h1} per v5 §2.4.4 L147
        if fct_pool_mode == "full_r9_6_sensitivity":
            return list(R9_6_ALL_LINES)
        raise ValueError(f"Unsupported fct_pool_mode={fct_pool_mode!r}")
    return list(R9_6_ALL_LINES)


def _subset_ds_eval_cols(sd: LineDatasetR10, eval_genes: Sequence[str], eval_g_to_idx: Dict[str, int]) -> Any:
    pos = {g: i for i, g in enumerate(sd.genes)}
    cols = np.asarray([pos[g] for g in eval_genes], dtype=np.int64)
    sub_sd = LineDatasetR10(
        name=sd.name,
        genes=list(eval_genes),
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
    return build_source_input_r10(sub_sd, only_reproducible=True, compute_cov=False)


def load_line_eval_reference_r16(
    target_line: str,
    model_genes: Sequence[str],
    model_pert_names: Sequence[str],
    hgnc_bridge: Optional[HGNCBridge] = None,
    require_self_gene_in_eval: bool = False,
    fct_pool_mode: str = "primary",
    line_dirs_override: Optional[Dict[str, Path]] = None,
) -> Dict[str, Any]:
    """
    Build standardized Round 16 reference evaluation data (`PREREG_I v5 §2.1–§2.4`):
      - Harmonizes `G_R9_6148`, `target_line` genes, source genes, and `G_model` (`2,000` HVGs)
        through frozen R10 `HGNCBridge`.
      - Uses `Replogle-4 \\ {t}` for Group 1, `r9_6` for Group 2, and `r9_6 \\ {vcc2025_h1}` primary
        for Group 3 (`gse295214_h1`), keeping `full_r9_6_sensitivity` strictly evaluated on the
        same primary locked targets (`Addendum A`).
    """
    hb = get_hgnc_bridge(hgnc_bridge)
    r9_common_canon = load_r9_common_genes(hgnc_bridge=hb)
    dir_map = line_dirs_override or STANDARD_LINE_DIRS

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

    src_lines = _resolve_fct_source_lines(target_line, fct_pool_mode=fct_pool_mode)

    tgt_ds = load_standardized_line_r10(
        dir_map[target_line],
        representation="log1p",
        common_genes=r9_common_canon,
        hgnc_bridge=hb,
        min_gene_overlap=100 if target_line == "gse295214_h1" else 1000,
    )
    min_ov_src = 50 if line_dirs_override is not None else 1000
    src_ds_map = {
        sn: load_standardized_line_r10(
            dir_map[sn],
            representation="log1p",
            common_genes=r9_common_canon,
            hgnc_bridge=hb,
            min_gene_overlap=min_ov_src,
        )
        for sn in src_lines
    }

    shared_set = set(tgt_ds.genes) & set(model_g_canon)
    for sd in src_ds_map.values():
        shared_set &= set(sd.genes)
    eval_genes = [g for g in r9_common_canon if g in shared_set]
    eval_g_to_idx = {g: i for i, g in enumerate(eval_genes)}

    pos_t = {g: i for i, g in enumerate(tgt_ds.genes)}
    t_cols = np.asarray([pos_t[g] for g in eval_genes], dtype=np.int64)
    tgt_deltas = np.asarray(tgt_ds.deltas[:, t_cols], dtype=np.float64)
    tgt_cmean_lib = np.asarray(tgt_ds.ctrl_mean_libnorm[t_cols], dtype=np.float64)
    tgt_gidx_full = np.asarray([eval_g_to_idx.get(t, -1) for t in tgt_ds.target_names], dtype=np.int64)
    tgt_ha = np.asarray(tgt_ds.half_a_deltas[:, t_cols], dtype=np.float64)
    tgt_hb = np.asarray(tgt_ds.half_b_deltas[:, t_cols], dtype=np.float64)

    p5_gene_set = set(tgt_ds.genes)
    src_inputs = []
    for sn in src_lines:
        sd = src_ds_map[sn]
        p5_gene_set &= set(sd.genes)
        src_inputs.append(_subset_ds_eval_cols(sd, eval_genes, eval_g_to_idx))

    fct_state = FctSourceState.fit(src_inputs, mask_self=True)
    src_rep_targets = set().union(*(set(si.target_names) for si in src_inputs))

    # Eligibility against P5 reproducible target set & canon(P_model)
    g_self_scope = set(tgt_ds.genes) if target_line == "gse295214_h1" else p5_gene_set
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
        if not require_self_gene_in_eval and t_sym not in g_self_scope:
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

    out: Dict[str, Any] = {
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
        "fct_pool_mode": fct_pool_mode,
        "g_eval_r9": len(tgt_ds.genes),
    }

    # Addendum A: For Group 3 (`gse295214_h1`), compute secondary `full r9_6` strictly on the SAME locked targets & genes
    if target_line == "gse295214_h1" and fct_pool_mode == "primary":
        sens_res = build_fct_sensitivity_on_locked_targets_r16(
            primary_ref_bundle=out,
            full_source_lines=R9_6_ALL_LINES,
            hgnc_bridge=hb,
            line_dirs_override=dir_map,
        )
        out["secondary_full_r9_6_same_targets"] = sens_res
    return out


def build_fct_sensitivity_on_locked_targets_r16(
    primary_ref_bundle: Dict[str, Any],
    full_source_lines: Sequence[str] = R9_6_ALL_LINES,
    hgnc_bridge: Optional[HGNCBridge] = None,
    line_dirs_override: Optional[Dict[str, Path]] = None,
) -> Dict[str, Any]:
    """
    Addendum A (`第 3 组敏感性分析固定在同一靶点集`):
    Refit `FCT` & `Transfer_Mean` using `full_source_lines` (e.g. `r9_6` incl. `vcc2025_h1`),
    evaluating strictly on `primary_ref_bundle["eval_perts"]` and `primary_ref_bundle["eval_genes"]`.
    """
    hb = get_hgnc_bridge(hgnc_bridge)
    r9_common_canon = load_r9_common_genes(hgnc_bridge=hb)
    dir_map = line_dirs_override or STANDARD_LINE_DIRS
    eval_genes = primary_ref_bundle["eval_genes"]
    eval_g_to_idx = {g: i for i, g in enumerate(eval_genes)}

    src_inputs = []
    for sn in full_source_lines:
        sd = load_standardized_line_r10(
            dir_map[sn],
            representation="log1p",
            common_genes=r9_common_canon,
            hgnc_bridge=hb,
            min_gene_overlap=1000,
        )
        src_inputs.append(_subset_ds_eval_cols(sd, eval_genes, eval_g_to_idx))

    fct_state = FctSourceState.fit(src_inputs, mask_self=True)
    preds_fct, _ = predict_group_c_suite(
        state=fct_state,
        target_ctrl_mean_libnorm=primary_ref_bundle["ctrl_mean_libnorm"],
        query_target_names=primary_ref_bundle["eval_perts"],
        query_target_gene_indices=primary_ref_bundle["target_gene_idx"],
        include_r9_anchors=False,
    )
    return {
        "source_lines_fct": list(full_source_lines),
        "eval_genes": eval_genes,
        "eval_perts": primary_ref_bundle["eval_perts"],
        "delta_transfer": np.asarray(preds_fct["Transfer_Mean"], dtype=np.float64),
        "source_grand_mean": np.asarray(fct_state.grand_mean, dtype=np.float64),
        "delta_fct": np.asarray(preds_fct["FCT"], dtype=np.float64),
    }
