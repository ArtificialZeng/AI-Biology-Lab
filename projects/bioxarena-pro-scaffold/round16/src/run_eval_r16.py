#!/usr/bin/env python3
"""
Round 16 Pre-Registered Evaluation Runner & Pre-Download GSE295214 Reader (`PREREG_I v5 + Addendum A`).

Hard Execution & Blinding Guards (`PREREG_I v5 §5`, Addendum A Items 1–11):
  1. Triple Hard Gate (`assert_all_execution_gates_r16`): Refuses ANY real-cell data access or
     model execution unless ALL THREE non-empty files exist in `round16/status/`:
       - `PREREG_I_COMMIT.txt`
       - `CODE_FREEZE_R16.txt`
       - `R16_GO.txt`
  2. Auditable Array, `X_state` Embedding, & Control-Barcode Exports (`round16/results/arrays/<line>/`,
     `round16/results/embeddings/`, `round16/results/control_barcodes_r16.json`):
     Exports per-line `observed.npz`, per-arm `preds_<arm>.npz` (including Item 8 within-line
     disjoint resample arm `seed = 20261602`), per-pool `2,058`-dim `X_state` embeddings
     `x_state_<safe_key>.npz` + `embeddings_manifest_sha256.json` (Item 11, per-pool seeded
     `20261601` primary / `20261602` resample for cross-party `<= 1e-5` verification prior to
     predicted deltas), exact sampled cell barcode lists for every control source (Item 7), and
     `array_manifest_sha256.json` for three-party verification.
  3. Group 3 (`GSE295214`, `H1 hESC`) Pre-Download Reader (`build_gse295214_parts`,
     `extract_gse295214_ntc_raw_counts`, `prepare_gse295214_line_r16`):
     Written and unit-tested strictly on synthetic 10x `.h5` fixtures BEFORE any download;
     implements exact `PREREG_I v5 §2.4` + Addendum A Items 1–2 rules (`full r9_6` sensitivity
     locked to primary `N_eval^R16` updating both `FCT` and `Transfer_Mean`).
  4. Arc Archived `adata_pred.h5ad` Adapter Check (Addendum A Item 3):
     Strictly fixed to `nadig_hepg2` ONLY (`median_per_pert_pearson >= 0.90` AND
     `pooled_pearson >= 0.95`; stops immediately for controller diagnosis on failure).
  5. Secondary Diagnostics: `STATE_HVG` `Perm100` (`k=0..99`, `default_rng(k)`, Item 4),
     `G4` 2x2 positional alignment on Group 1 and `(0, 17)` + `(17, 0)` symmetry check on
     Groups 2–3 (Item 6), within-line disjoint resample contrast (Items 8–9), and
     `EXACT_INVARIANCE` boundary handling (Item 10). `CellularAI_OpenCode_PublicRefit` is
     explicitly dropped pre-data (Item 5).
  6. `--synthetic-smoke`: Runs the complete evaluation, adjudication, `Perm100`, `G4` 2x2 +
     `(17, 0)` symmetry, Item 8 within-line resample, Item 11 `X_state` embedding export,
     barcode export, and `.npz`/`.json` pipeline strictly on 100% synthetic tensors (ZERO real
     cells touched).
"""

from __future__ import annotations

import argparse
import csv
import gc
import gzip
import hashlib
import json
import os
import pickle
import re
import sys
import tempfile
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Set, Tuple

import h5py
import numpy as np
import scipy.sparse as sp

ROUND16_SRC_DIR = Path(__file__).resolve().parent
ROUND16_ROOT = ROUND16_SRC_DIR.parent
REPO_ROOT = ROUND16_ROOT.parent
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
R10_SCRIPTS_DIR = REPO_ROOT / "round10" / "shared" / "scripts"

for _p in (str(ROUND16_SRC_DIR), str(R10_SRC_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from adjudicate_r16 import (  # noqa: E402
    DEFAULT_TOST_DELTA,
    adjudicate_claim_i1_line,
    adjudicate_claim_i2_line,
    adjudicate_group_r16,
    summarize_within_line_resample_r16,
)
from data_loader_r10 import (  # noqa: E402
    DEFAULT_HGNC_FILE,
    HGNCBridge,
    LineDatasetR10,
    R9_COMMON_GENES_FILE,
    R9_LINES_DIR,
    R9_SOURCE_LINE_ORDER,
    R10_LINES_DIR,
    decode_str_array,
    load_r9_common_genes,
    load_standardized_line_r10,
)
from fct_pipeline import (  # noqa: E402
    FctSourceState,
    build_source_input_r10,
    compute_source_transfer_and_fc_bar,
    fct_from_fc_bar,
)
from metrics_r16 import (  # noqa: E402
    DEFAULT_R16_SEED,
    DEFAULT_RESAMPLE_SEED,
    compute_gene_label_perm_null_r16,
    compute_joint_cross_mask_r16,
    compute_orth_suite_r16,
    compute_paired_contrast_r16,
    compute_paired_model_vs_all_mismatched_r16,
    compute_paired_with_within_line_resample_r16,
    compute_winner_curse_weight_shift_r16,
)
from model_adapters_r16 import (  # noqa: E402
    DEFAULT_MODELS_DIR,
    DEFAULT_STATUS_DIR,
    GROUP2_UNSEEN_PRIMARY_LINES,
    GROUP3_PRIMARY_FCT_SOURCES,
    GROUP3_PROSPECTIVE_LINES,
    REPLOGLE4_LOO_PRIMARY_LINES,
    RETAINED_UNSEEN_MISMATCHED_ONLY_LINES,
    ArcStateAdapter,
    assert_all_execution_gates_r16,
    assert_prereg_i_locked,
    assert_zero_training_leakage,
    build_mismatched_strata_r16,
    compute_embedding_separation_diagnostics,
    encode_controls_se600m_r16,
    preprocess_hvg_counts_r16,
    sample_disjoint_resample_control_sets_r16,
    sample_stratified_control_sets_r16,
)

SYNC_EXTERNAL_DATA_DIR = Path(
    "/Volumes/雷雳传输/A顶会论文同步盘/A_BIG_MAC顶会论文写稿/AI自主研究/gemini/external_data"
)
DEFAULT_RESULTS_DIR = ROUND16_ROOT / "results"
DEFAULT_R16_SHARED_DIR = ROUND16_ROOT / "shared"

REPLOGLE4_LINES: Tuple[str, ...] = (
    "nadig_hepg2",
    "nadig_jurkat",
    "replogle_k562_essential",
    "replogle_rpe1",
)

LINE_SPECS_R16: Dict[str, Dict[str, Any]] = {
    "nadig_hepg2": {
        "group": "group1",
        "modality": "CRISPRi",
        "checkpoint_id": "zeroshot/hepg2",
        "fct_sources": ("nadig_jurkat", "replogle_k562_essential", "replogle_rpe1"),
        "line_dir": R9_LINES_DIR / "nadig_hepg2",
    },
    "nadig_jurkat": {
        "group": "group1",
        "modality": "CRISPRi",
        "checkpoint_id": "zeroshot/jurkat",
        "fct_sources": ("nadig_hepg2", "replogle_k562_essential", "replogle_rpe1"),
        "line_dir": R9_LINES_DIR / "nadig_jurkat",
    },
    "a549_ko": {
        "group": "group2",
        "modality": "CRISPRko",
        "checkpoint_id": "zeroshot/ensemble4",
        "fct_sources": R9_SOURCE_LINE_ORDER,
        "line_dir": REPO_ROOT / "round12" / "shared" / "lines" / "a549_ko",
    },
    "dld1": {
        "group": "group2",
        "modality": "CRISPRi",
        "checkpoint_id": "zeroshot/ensemble4",
        "fct_sources": R9_SOURCE_LINE_ORDER,
        "line_dir": R10_LINES_DIR / "dld1",
    },
    "cd4t_rest": {
        "group": "group2",
        "modality": "CRISPRi",
        "checkpoint_id": "zeroshot/ensemble4",
        "fct_sources": R9_SOURCE_LINE_ORDER,
        "line_dir": R10_LINES_DIR / "cd4t_rest",
    },
    "feng_ipsc": {
        "group": "sensitivity",
        "modality": "CRISPRi",
        "checkpoint_id": "zeroshot/ensemble4",
        "fct_sources": R9_SOURCE_LINE_ORDER,
        "line_dir": REPO_ROOT / "round11" / "shared" / "lines" / "feng_ipsc",
    },
    "gse295214_h1": {
        "group": "group3",
        "modality": "CRISPRi",
        "checkpoint_id": "zeroshot/ensemble4",
        "fct_sources": GROUP3_PRIMARY_FCT_SOURCES,
        "fct_sources_sensitivity": R9_SOURCE_LINE_ORDER,
        "line_dir": DEFAULT_R16_SHARED_DIR / "lines" / "gse295214_h1",
    },
}

ARC_ADAPTER_CHECK_LINE = "nadig_hepg2"
ARC_ADAPTER_MEDIAN_PEARSON_MIN = 0.90
ARC_ADAPTER_POOLED_PEARSON_MIN = 0.95
G4_SPREAD_BATCH_INDICES: Tuple[int, ...] = (0, 1, 17, 35)

_GSM_RE = re.compile(r"(GSM\d+)", re.IGNORECASE)
_ENSG_RE = re.compile(r"(ENSG\d+)", re.IGNORECASE)
_MULTI_DELIM_RE = re.compile(r"[,;\s]")
_ILLEGAL_PLUS_RE = re.compile(r"(?<!_)\+|\+(?![_\d])")
_MULTI_GUIDE_DELIM_RE = re.compile(r"[,;|+]")
_SAFE_OR_BAD_GUIDE_TOKENS = ("safe", "intergenic", "unassigned", "multiplet", "doublet", "gene_desert")
_GUIDE_SUFFIX_TOKEN_RE = re.compile(
    r"^(?:P\d+(?:[-]?P\d+)*|ENST\d+|ENSG\d+|\d+(?:\.\d+)?(?:-P\d+.*)?|sg\d+|TSS\d*|[+-]|[+-]\d+.*|AB|A|B|[ACGT]{15,})$",
    re.IGNORECASE,
)
_DLD1_LOW_MOI_CMOS = frozenset({f"A{i}" for i in range(2, 13)} | {f"B{i}" for i in range(1, 8)})
_CD4_EXCEL_GUIDE_REMAP = {"1-Jun": "JUN-1", "2-Jun": "JUN-2"}
_A549_NTC_TOKENS = frozenset({"nontargeting", "non-targeting", "ntc", "nontarget"})
_A549_UNASSIGNED_TOKENS = frozenset(
    {"", "na", "nan", "none", "null", "unassigned", "multiplet", "doublet", "ambiguous", "unknown"}
)
_FENG_KOLF_EXCLUDED = frozenset({"kolf_2", "kolf_3", "kolf2", "kolf3", "kolf_2_1j"})


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _safe_arm_filename(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(name).strip())


# ==============================================================================
# Exact-Parity Pure-NumPy/SciPy Helpers (Work Identically in `.venv` & Anaconda)
# ==============================================================================


@dataclass
class LightweightPartSpec:
    """Drop-in compatible mirror of `prepare_r10_lines.PartSpec` when statsmodels is absent."""

    tag: str
    path: Path
    n_cells: int
    n_genes_raw: int
    nt_mask: np.ndarray
    batches: np.ndarray
    half_flags: np.ndarray
    umi: np.ndarray
    pert_keys: List[str]
    pert_meta: Dict[str, Tuple[str, str, str, str]]
    var_gene_ids: List[str]
    var_gene_names: List[str]
    read_chunk_fn: Callable[[int, int], np.ndarray]
    read_rows_fn: Optional[Callable[[np.ndarray], np.ndarray]] = None
    read_chunk_cols_fn: Optional[Callable[[int, int, np.ndarray], np.ndarray]] = None
    cell_barcodes: Optional[List[str]] = None


def _fit_poisson_gaussian_thresholds_r16(
    hist: np.ndarray,
    min_thr: int = 5,
    max_thr: int = 20,
    default_thr: int = 7,
    min_entries: int = 50,
) -> np.ndarray:
    """Exact bit-for-bit implementation of frozen `prepare_r10_lines.fit_poisson_gaussian_thresholds_from_hist`."""
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


def _extract_top2_csc_r16(csc_mat: sp.csc_matrix) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Exact bit-for-bit implementation of frozen `prepare_r10_lines.extract_top2_csc`."""
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


def _apply_min_ntc_batch_fallback_r16(
    batches: np.ndarray,
    nt_mask: np.ndarray,
    half_flags: np.ndarray,
    min_batch_nt: int = 15,
    group_labels: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, int]:
    """Exact within-half `< min_batch_nt` NTC batch fallback (optional group fallback; zero cross-half sharing)."""
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
        if nt_counts.get(str(out[i]), 0) >= min_batch_nt:
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


def _compute_10x_csc_rna_qc_r16(
    h5_path: Path,
    rna_row_idx: Optional[np.ndarray],
    mt_mask_in_raw: np.ndarray,
    batch_cols: int = 20000,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Exact mirror of frozen `prepare_r10_lines.compute_10x_csc_rna_qc` (no statsmodels dependency)."""
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
            umi[cs:ce] = np.bincount(col_ids - cs, weights=dat, minlength=ce - cs)
            ngenes[cs:ce] = np.bincount(col_ids - cs, minlength=ce - cs)
            mt_m = mt_mask_in_raw[idx]
            mt_umi[cs:ce] = (
                np.bincount(col_ids[mt_m] - cs, weights=dat[mt_m], minlength=ce - cs)
                if np.any(mt_m)
                else 0.0
            )
        pct_mt = 100.0 * mt_umi / np.maximum(umi, 1.0)
    return umi, ngenes, pct_mt


def _make_10x_csc_subset_readers_r16(
    h5_path: Path,
    keep_cols: np.ndarray,
    n_raw_features: int,
    rna_row_idx: Optional[np.ndarray],
) -> Tuple[
    Callable[[int, int], np.ndarray],
    Callable[[np.ndarray], np.ndarray],
    Callable[[int, int, np.ndarray], np.ndarray],
]:
    """Exact mirror of `prepare_r10_lines.make_10x_csc_subset_readers` (memory-bounded CSC reads)."""
    n_out_genes = len(rna_row_idx) if rna_row_idx is not None else n_raw_features
    if (
        rna_row_idx is not None
        and len(rna_row_idx) == n_raw_features
        and np.array_equal(rna_row_idx, np.arange(n_raw_features))
    ):
        rna_row_idx = None

    def _read_cols(raw_cols: np.ndarray, feat_cols: Optional[np.ndarray] = None) -> np.ndarray:
        out_dtype = np.float32 if feat_cols is not None else np.float64
        out_w = len(feat_cols) if feat_cols is not None else n_out_genes
        if len(raw_cols) == 0:
            return np.zeros((0, out_w), dtype=out_dtype)
        order = np.argsort(raw_cols)
        sorted_cols = raw_cols[order]
        c_min, c_max = int(sorted_cols[0]), int(sorted_cols[-1])
        span = c_max - c_min + 1
        with h5py.File(h5_path, "r") as f:
            m = f["matrix"]
            indptr = m["indptr"][:]
            data = m["data"]
            indices = m["indices"]
            if len(sorted_cols) / max(span, 1) >= 0.15:
                ptrs = indptr[c_min : c_max + 2]
                p0, p1 = int(ptrs[0]), int(ptrs[-1])
                csc_span = sp.csc_matrix(
                    (data[p0:p1], indices[p0:p1], ptrs - p0),
                    shape=(n_raw_features, span),
                )
                sub = csc_span[:, sorted_cols - c_min]
            else:
                col_starts = indptr[sorted_cols]
                col_ends = indptr[sorted_cols + 1]
                all_d, all_i = [], []
                new_ptr = np.zeros(len(sorted_cols) + 1, dtype=np.int64)
                cur_p = 0
                for j, (s0, s1) in enumerate(zip(col_starts, col_ends)):
                    if s1 > s0:
                        all_d.append(data[int(s0) : int(s1)])
                        all_i.append(indices[int(s0) : int(s1)])
                        cur_p += int(s1 - s0)
                    new_ptr[j + 1] = cur_p
                d_cat = np.concatenate(all_d) if all_d else np.empty(0, dtype=np.float32)
                i_cat = np.concatenate(all_i) if all_i else np.empty(0, dtype=np.int32)
                sub = sp.csc_matrix((d_cat, i_cat, new_ptr), shape=(n_raw_features, len(sorted_cols)))

            if feat_cols is not None:
                target_rows = rna_row_idx[feat_cols] if rna_row_idx is not None else feat_cols
                sub = sub[target_rows, :]
            elif rna_row_idx is not None:
                sub = sub[rna_row_idx, :]
            dense = sub.T.toarray().astype(out_dtype, copy=False)
            if not np.array_equal(order, np.arange(len(raw_cols))):
                inv = np.empty_like(order)
                inv[order] = np.arange(len(raw_cols))
                dense = dense[inv]
            return dense

    return (
        (lambda s, e: _read_cols(keep_cols[s:e], None)),
        (lambda loc_idx: _read_cols(keep_cols[loc_idx], None)),
        (lambda s, e, fcols: _read_cols(keep_cols[s:e], fcols)),
    )


# ==============================================================================
# Group 3 (`GSE295214`, H1 hESC) Pre-Download 10x `.h5` Reader (`PREREG_I v5 §2.4`)
# ==============================================================================


def _resolve_single_guide_token(
    token: str,
    hgnc: Any,
) -> Optional[Tuple[str, str, str, str]]:
    """Resolve one targeting guide token to `(ens_id, canon_sym, raw_sym, status)` or None."""
    tok = token.strip()
    if not tok or _MULTI_DELIM_RE.search(tok) or _ILLEGAL_PLUS_RE.search(tok):
        return None
    low = tok.lower()
    if low.startswith("non-targeting_") or low == "non-targeting":
        return None
    if any(bad in low for bad in _SAFE_OR_BAD_GUIDE_TOKENS) or "-and-" in low:
        return None

    ens_m = _ENSG_RE.search(tok)
    ens_raw = ens_m.group(1) if ens_m else ""

    def _try_resolve(sym_cand: str, ens_cand: str) -> Optional[Tuple[str, str, str, str]]:
        if hasattr(hgnc, "canonicalize_with_source"):
            cs, st = hgnc.canonicalize_with_source(sym_cand, ensembl_id=ens_cand or None)
            if cs is not None and cs in hgnc.approved:
                ens_out = ens_cand or getattr(hgnc, "sym_to_ens", {}).get(cs, "")
                return (ens_out, cs, sym_cand, st)
            return None
        cs, st = hgnc.resolve(sym_cand, ens_cand)
        if cs in hgnc.approved:
            ens_out = ens_cand or hgnc.sym_to_ens.get(cs, "")
            return (ens_out, cs, sym_cand, st)
        return None

    # 1. Exact match on full token
    direct = _try_resolve(tok, ens_raw)
    if direct is not None:
        return direct

    # 2. Underscore-delimited tokens (e.g., `GENE_1`, `GENE_P1P2_1`, `GENE_ENSG..._P1`, `ENSG..._GENE_1`)
    parts = [p for p in tok.split("_") if p]
    if len(parts) >= 2:
        # Reject multi-gene constructs where two distinct approved gene symbols appear
        distinct_gene_hits: Set[str] = set()
        for p in parts:
            if _GUIDE_SUFFIX_TOKEN_RE.match(p):
                continue
            r_p = _try_resolve(p, "")
            if r_p is not None:
                distinct_gene_hits.add(r_p[1])
        if len(distinct_gene_hits) > 1:
            return None
        first_cand = parts[0]
        if first_cand.upper().startswith("ENSG") and len(parts) >= 2:
            hit = _try_resolve(parts[1], first_cand.split(".")[0])
            if hit is not None:
                return hit
        hit = _try_resolve(first_cand, ens_raw)
        if hit is not None:
            return hit

    # 3. Trailing hyphen guide index (`GENE-1`, `GENE-sg1`, `GENE-P1P2`, `GENE-TSS1`) ONLY when full token did not match
    m_hyph = re.match(r"^(.+?)-(?:\d+|sg\d+|P\d+(?:[-]?P\d+)*|TSS\d*)$", tok, re.IGNORECASE)
    if m_hyph:
        hit = _try_resolve(m_hyph.group(1), ens_raw)
        if hit is not None:
            return hit

    return None


def parse_gse295214_guide_r16(
    guide_id: str,
    guide_name: str,
    hgnc: Any,
) -> Tuple[str, Optional[Tuple[str, str, str, str]]]:
    """
    Classify a `GSE295214` `CRISPR Guide Capture` feature per `PREREG_I v5 §2.4`:
      - Returns `("ntc", None)` strictly for `non-targeting_*` guides (safe-targeting &
        other non-gene controls strictly excluded).
      - Returns `("target", (ens_id, canon_sym, raw_sym, status))` for unambiguous
        single-gene HGNC-approved guides (supporting dual-guide single-gene `A|B` cassettes
        iff both halves resolve to the identical canonical HGNC gene).
      - Returns `("exclude", None)` otherwise.
    """
    gid = str(guide_id).strip()
    gnm = str(guide_name).strip()
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
        resolved = [_resolve_single_guide_token(s, hgnc) for s in subs]
        if any(r is None for r in resolved):
            return "exclude", None
        canon_syms = {r[1] for r in resolved if r is not None}
        if len(canon_syms) != 1:
            return "exclude", None
        return "target", resolved[0]

    res = _resolve_single_guide_token(raw, hgnc)
    if res is None and gid and gid != raw:
        res = _resolve_single_guide_token(gid, hgnc)
    if res is not None:
        return "target", res
    return "exclude", None


def _extract_gsm_channel_id(h5_path: Path) -> str:
    m = _GSM_RE.search(h5_path.name)
    return m.group(1).upper() if m else h5_path.stem


def build_gse295214_parts(
    raw_dir: Path,
    hgnc: Any,
    min_umi: int = 300,
    min_ngenes: int = 200,
    min_thr: int = 5,
    max_thr: int = 20,
    default_thr: int = 7,
    min_batch_nt: int = 15,
    batch_cols: int = 25000,
    use_r10_partspec: bool = False,
) -> Tuple[List[Any], Dict[str, Any]]:
    """
    Pre-download locked `GSE295214` (`H1 hESC`, 85 GSM 10x CellRanger `.h5` channels) reader
    (`PREREG_I v5 §2.4`):
      - Verifies 100% of stored `Gene Expression` non-zero counts are finite non-negative integers.
      - Pools `CRISPR Guide Capture` histograms across cells with `n_genes >= 200` (`max_pct_mt=None`)
        and fits exact `PREREG_F §3.2` bimodal Poisson-Gaussian thresholds (`min_thr=5, max_thr=20,
        default_thr=7`).
      - Retains singlets with `n_umis >= 300`, `n_genes >= 200`, `max_pct_mt=None`, `n_pos == 1`,
        `top1 >= 5`, and `(top2 < 3 or top1 / max(top2, 0.5) >= 3.0)`.
      - Strictly `non-targeting_*` NTCs; unambiguous single-gene HGNC targeting guides.
      - `PREREG_H §6` Tier-1 lexicographical GSM parity: 1-based odd -> Half A (`1`), even -> Half B (`2`),
        with global within-half `< 15` NTC fallback (zero cross-half sharing).
    """
    rdir = Path(raw_dir)
    if (rdir / "gse295214_h1").is_dir():
        rdir = rdir / "gse295214_h1"
    elif (rdir / "GSE295214").is_dir():
        rdir = rdir / "GSE295214"

    h5_files_raw = sorted(rdir.glob("*.h5"))
    if not h5_files_raw:
        raise FileNotFoundError(f"No 10x .h5 channel files found in {rdir}")

    # Sort channels strictly in Python lexicographical order of channel ID (`PREREG_H §6`)
    ch_pairs = sorted(
        [(_extract_gsm_channel_id(p), p) for p in h5_files_raw],
        key=lambda x: (x[0], x[1].name),
    )
    ch_ids = [c for c, _ in ch_pairs]
    if len(set(ch_ids)) != len(ch_ids):
        raise ValueError(f"Duplicate GSM channel IDs detected in {rdir}: {ch_ids[:10]}")

    PartCls: Any = LightweightPartSpec
    if use_r10_partspec:
        if str(R10_SCRIPTS_DIR) not in sys.path:
            sys.path.insert(0, str(R10_SCRIPTS_DIR))
        from prepare_r10_lines import PartSpec  # type: ignore

        PartCls = PartSpec

    lane_cache = []
    hist: Optional[np.ndarray] = None
    crp_ids_ref: Optional[List[str]] = None
    crp_names_ref: Optional[List[str]] = None
    tot_raw_cells = 0
    tot_hq_rna_cells = 0

    for ch_pos_1based, (ch_id, lp) in enumerate(ch_pairs, start=1):
        half_val = 1 if (ch_pos_1based % 2 == 1) else 2
        with h5py.File(lp, "r") as f:
            if "matrix" not in f:
                raise KeyError(f"EXCLUDED_FEASIBILITY: {lp} missing CellRanger /matrix group.")
            m = f["matrix"]
            n_feat, n_col = int(m["shape"][0]), int(m["shape"][1])
            raw_bcs = decode_str_array(m["barcodes"][:]) if "barcodes" in m else [f"cell_{i}" for i in range(n_col)]
            ftypes = np.array(decode_str_array(m["features/feature_type"][:]))
            fids = np.array(decode_str_array(m["features/id"][:]))
            fnames = np.array(decode_str_array(m["features/name"][:]))
            rna_rows = np.where(ftypes == "Gene Expression")[0]
            crp_rows = np.where(ftypes != "Gene Expression")[0]
            if len(rna_rows) == 0 or len(crp_rows) == 0:
                raise ValueError(
                    f"EXCLUDED_FEASIBILITY: {lp} lacks Gene Expression ({len(rna_rows)}) or "
                    f"CRISPR Guide Capture ({len(crp_rows)}) rows."
                )

            cur_crp_ids = [str(fids[i]) for i in crp_rows]
            cur_crp_names = [str(fnames[i]) for i in crp_rows]
            if crp_names_ref is None:
                crp_ids_ref = cur_crp_ids
                crp_names_ref = cur_crp_names
                hist = np.zeros((len(crp_rows), 66), dtype=np.int64)
            elif cur_crp_names != crp_names_ref or cur_crp_ids != crp_ids_ref:
                raise ValueError(f"CRISPR guide feature axis mismatch in {lp} vs {ch_pairs[0][1]}")

            is_crp = np.zeros(n_feat, dtype=bool)
            is_crp[crp_rows] = True
            crp_remap_idx = np.full(n_feat, -1, dtype=np.int32)
            crp_remap_idx[crp_rows] = np.arange(len(crp_rows), dtype=np.int32)
            is_rna_feat = ~is_crp

            umi_l = np.zeros(n_col, dtype=np.float64)
            ng_l = np.zeros(n_col, dtype=np.int32)
            crp_idx_parts, crp_dat_parts, crp_col_parts = [], [], []

            indptr_all = m["indptr"][:]
            data_ds = m["data"]
            indices_ds = m["indices"]
            for cs in range(0, n_col, batch_cols):
                ce = min(cs + batch_cols, n_col)
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
                    # Mechanical 100% integer raw-UMI check (`PREREG_H §5.1` / `PREREG_I §2.4`)
                    if not np.all(
                        np.isfinite(r_dat) & (r_dat >= 0.0) & (np.abs(r_dat - np.round(r_dat)) <= 1e-6)
                    ):
                        raise ValueError(
                            f"EXCLUDED_FEASIBILITY: non-integer or negative Gene Expression counts in {lp}"
                        )
                    r_col = col_ids[m_r] - cs
                    umi_l[cs:ce] = np.bincount(r_col, weights=r_dat, minlength=ce - cs)
                    pos_r = r_dat > 0.0
                    if np.any(pos_r):
                        ng_l[cs:ce] = np.bincount(r_col[pos_r], minlength=ce - cs)

                m_c = ~m_r
                if np.any(m_c):
                    c_d = dat[m_c].astype(np.float64, copy=False)
                    pos_c = c_d > 0.0
                    if np.any(pos_c):
                        crp_idx_parts.append(crp_remap_idx[idx[m_c][pos_c]])
                        crp_dat_parts.append(np.rint(c_d[pos_c]).astype(np.int64))
                        crp_col_parts.append(col_ids[m_c][pos_c])

        tot_raw_cells += n_col
        # PREREG_F §3.2 / PREREG_I §2.4: histogram over cells with >= 200 genes, NO MT% cutoff
        ok_init = ng_l >= int(min_ngenes)
        ok_hq = (umi_l >= float(min_umi)) & (ng_l >= int(min_ngenes))
        tot_hq_rna_cells += int(np.sum(ok_hq))

        c_idx = np.concatenate(crp_idx_parts) if crp_idx_parts else np.empty(0, dtype=np.int32)
        c_dat = np.concatenate(crp_dat_parts) if crp_dat_parts else np.empty(0, dtype=np.int64)
        c_col = np.concatenate(crp_col_parts) if crp_col_parts else np.empty(0, dtype=np.int32)

        if c_col.size > 0 and hist is not None:
            ok_init_mask = ok_init[c_col]
            if np.any(ok_init_mask):
                bins = np.clip(c_dat[ok_init_mask], 1, 65)
                np.add.at(hist, (c_idx[ok_init_mask], bins), 1)

        hq_cols = np.where(ok_hq)[0]
        col_map_hq = np.full(n_col, -1, dtype=np.int32)
        col_map_hq[hq_cols] = np.arange(len(hq_cols), dtype=np.int32)
        if c_col.size > 0 and len(hq_cols) > 0:
            hq_m = col_map_hq[c_col] >= 0
            C_hq = sp.csr_matrix(
                (c_dat[hq_m].astype(np.float64), (c_idx[hq_m], col_map_hq[c_col[hq_m]])),
                shape=(len(crp_rows), len(hq_cols)),
            ).tocsc()
        else:
            C_hq = sp.csc_matrix((len(crp_rows), len(hq_cols)), dtype=np.float64)

        lane_cache.append(
            (
                ch_id,
                half_val,
                lp,
                n_feat,
                rna_rows,
                [str(fids[i]) for i in rna_rows],
                [str(fnames[i]) for i in rna_rows],
                raw_bcs,
                umi_l,
                hq_cols,
                C_hq,
            )
        )

    assert hist is not None and crp_ids_ref is not None and crp_names_ref is not None
    thrs = _fit_poisson_gaussian_thresholds_r16(
        hist, min_thr=min_thr, max_thr=max_thr, default_thr=default_thr
    )

    guide_class = [
        parse_gse295214_guide_r16(gid, gnm, hgnc) for gid, gnm in zip(crp_ids_ref, crp_names_ref)
    ]
    n_ntc_guides = sum(1 for k, _ in guide_class if k == "ntc")
    n_tgt_guides = sum(1 for k, _ in guide_class if k == "target")
    n_excl_guides = sum(1 for k, _ in guide_class if k == "exclude")

    # Stage 1: collect per-channel singlets so we can run exact global `< 15` NTC batch fallback
    staged: List[Dict[str, Any]] = []
    all_batches_l: List[np.ndarray] = []
    all_nt_l: List[np.ndarray] = []
    all_half_l: List[np.ndarray] = []
    tot_ntc = tot_tgt = tot_excluded_guide_cells = 0

    for ch_id, half_val, lp, n_feat, rna_rows, rids, rnames, raw_bcs, umi_l, hq_cols, C_hq in lane_cache:
        t1, t2, a1 = _extract_top2_csc_r16(C_hq)
        pos_mask = C_hq.data >= thrs[C_hq.indices]
        pos_csc = sp.csc_matrix((pos_mask.astype(np.int8), C_hq.indices, C_hq.indptr), shape=C_hq.shape)
        n_pos = np.asarray(pos_csc.sum(axis=0)).ravel()
        s1_mask = (n_pos == 1) & (t1 >= 5.0) & ((t2 < 3.0) | (t1 / np.maximum(t2, 0.5) >= 3.0))

        keep_cols: List[int] = []
        nt_flags: List[bool] = []
        pkeys: List[str] = []
        pmeta: Dict[str, Tuple[str, str, str, str]] = {}

        for loc_j in np.where(s1_mask)[0]:
            g_i = int(a1[loc_j])
            g_kind, g_info = guide_class[g_i]
            raw_col = int(hq_cols[loc_j])
            if g_kind == "ntc":
                keep_cols.append(raw_col)
                nt_flags.append(True)
                pkeys.append("")
                tot_ntc += 1
            elif g_kind == "target" and g_info is not None:
                gkey = crp_names_ref[g_i] or crp_ids_ref[g_i]
                keep_cols.append(raw_col)
                nt_flags.append(False)
                pkeys.append(gkey)
                pmeta.setdefault(gkey, g_info)
                tot_tgt += 1
            else:
                tot_excluded_guide_cells += 1

        keep_arr = np.asarray(keep_cols, dtype=np.int64)
        ord_k = np.argsort(keep_arr)
        keep_arr = keep_arr[ord_k]
        nt_arr = np.asarray(nt_flags, dtype=bool)[ord_k]
        pkeys_ord = [pkeys[int(i)] for i in ord_k]
        bcs_kept = [f"{ch_id}:{raw_bcs[int(c)]}" for c in keep_arr]
        b_arr = np.full(len(keep_arr), ch_id, dtype=object)
        h_arr = np.full(len(keep_arr), half_val, dtype=np.int8)

        staged.append(
            {
                "ch_id": ch_id,
                "half_val": half_val,
                "lp": lp,
                "n_feat": n_feat,
                "rna_rows": rna_rows,
                "rids": rids,
                "rnames": rnames,
                "bcs_kept": bcs_kept,
                "umi_l": umi_l,
                "keep_arr": keep_arr,
                "nt_arr": nt_arr,
                "pkeys": pkeys_ord,
                "pmeta": pmeta,
            }
        )
        all_batches_l.append(b_arr)
        all_nt_l.append(nt_arr)
        all_half_l.append(h_arr)

    if sum(len(x) for x in all_nt_l) == 0:
        raise RuntimeError(f"GSE295214 yielded zero QC-passing singlet cells across {len(ch_pairs)} channels.")

    flat_batches = np.concatenate(all_batches_l)
    flat_nt = np.concatenate(all_nt_l)
    flat_half = np.concatenate(all_half_l)
    fb_batches, n_fb_cells = _apply_min_ntc_batch_fallback_r16(
        flat_batches, flat_nt, flat_half, min_batch_nt=min_batch_nt
    )

    parts: List[Any] = []
    offset = 0
    for item in staged:
        n_c = len(item["keep_arr"])
        b_slice = fb_batches[offset : offset + n_c]
        h_slice = flat_half[offset : offset + n_c]
        offset += n_c
        if n_c == 0:
            continue
        read_chunk, read_rows, read_chunk_cols = _make_10x_csc_subset_readers_r16(
            item["lp"], item["keep_arr"], item["n_feat"], item["rna_rows"]
        )
        part_kwargs: Dict[str, Any] = {
            "tag": item["ch_id"],
            "path": item["lp"],
            "n_cells": n_c,
            "n_genes_raw": len(item["rna_rows"]),
            "nt_mask": item["nt_arr"],
            "batches": b_slice,
            "half_flags": h_slice,
            "umi": item["umi_l"][item["keep_arr"]],
            "pert_keys": item["pkeys"],
            "pert_meta": item["pmeta"],
            "var_gene_ids": item["rids"],
            "var_gene_names": item["rnames"],
            "read_chunk_fn": read_chunk,
            "read_rows_fn": read_rows,
            "read_chunk_cols_fn": read_chunk_cols,
        }
        if not use_r10_partspec:
            part_kwargs["cell_barcodes"] = item["bcs_kept"]
        parts.append(PartCls(**part_kwargs))

    n_ha_ntc = int(np.sum(flat_nt & (flat_half == 1)))
    n_hb_ntc = int(np.sum(flat_nt & (flat_half == 2)))
    audit = {
        "n_channels": len(ch_pairs),
        "channels_sorted_sample": ch_ids[:5] + (["..."] + ch_ids[-3:] if len(ch_ids) > 8 else []),
        "half_a_n_channels": sum(1 for i in range(1, len(ch_pairs) + 1) if i % 2 == 1),
        "half_b_n_channels": sum(1 for i in range(1, len(ch_pairs) + 1) if i % 2 == 0),
        "n_raw_cells": int(tot_raw_cells),
        "n_hq_rna_cells": int(tot_hq_rna_cells),
        "qc_spec": {"min_umi": int(min_umi), "min_ngenes": int(min_ngenes), "max_pct_mt": None},
        "n_guide_features": len(crp_names_ref),
        "n_ntc_guides": int(n_ntc_guides),
        "n_targeting_guides": int(n_tgt_guides),
        "n_excluded_guides": int(n_excl_guides),
        "median_guide_threshold": float(np.median(thrs)),
        "n_kept_ntc": int(tot_ntc),
        "n_kept_targeting": int(tot_tgt),
        "n_excluded_guide_singlet_cells": int(tot_excluded_guide_cells),
        "n_kept_half_a_ntc": n_ha_ntc,
        "n_kept_half_b_ntc": n_hb_ntc,
        "per_half_ntc_ge15": bool(n_ha_ntc >= 15 and n_hb_ntc >= 15),
        "n_low_ntc_channel_fallback_cells": int(n_fb_cells),
        "half_rule": "PREREG_H_Tier1_lexicographical_GSM_1based_odd_HalfA_even_HalfB",
    }
    return parts, audit


def _read_part_rows_by_indices(
    parts: Sequence[Any],
    ntc_refs: Sequence[Tuple[int, int]],
    flat_idx: np.ndarray,
    n_genes: int,
) -> np.ndarray:
    """Read exact `flat_idx` rows across `parts` into a `(len(flat_idx), n_genes)` float32 matrix."""
    raw_out = np.empty((len(flat_idx), n_genes), dtype=np.float32)
    by_part: Dict[int, List[Tuple[int, int]]] = defaultdict(list)
    for out_pos, g_i in enumerate(flat_idx):
        p_idx, loc_i = ntc_refs[int(g_i)]
        by_part[p_idx].append((out_pos, loc_i))

    for p_idx, pairs in sorted(by_part.items()):
        part = parts[p_idx]
        loc_arr = np.asarray([loc_i for _, loc_i in pairs], dtype=np.int64)
        out_arr = np.asarray([out_pos for out_pos, _ in pairs], dtype=np.int64)
        ord_loc = np.argsort(loc_arr)
        loc_sorted = loc_arr[ord_loc]
        if part.read_rows_fn is not None:
            rows_mat = np.asarray(part.read_rows_fn(loc_sorted), dtype=np.float32)
        else:
            rows_mat = np.vstack(
                [part.read_chunk_fn(int(r), int(r) + 1) for r in loc_sorted]
            ).astype(np.float32, copy=False)
        raw_out[out_arr[ord_loc]] = rows_mat
    return raw_out


def extract_ntc_bundle_from_parts_r16(
    parts: Sequence[Any],
    target_n: int = 512,
    cell_set_len: int = 32,
    min_ctrl: int = 128,
    seed: int = DEFAULT_R16_SEED,
    resample_seed: int = DEFAULT_RESAMPLE_SEED,
    include_resample: bool = True,
) -> Dict[str, Any]:
    """
    Sample primary (`seed = 20261601`, `PREREG_I v5 §2.2`) and optional disjoint within-line
    resample (`seed = 20261602`, Addendum A Item 8) QC-passing NTC cells from `PartSpec` /
    `LightweightPartSpec` objects, preserving exact ordered cell barcodes (Addendum A Item 7).
    """
    ref_genes = list(parts[0].var_gene_names)
    ntc_refs: List[Tuple[int, int]] = []
    ntc_batches: List[str] = []
    ntc_barcodes: List[str] = []
    for p_idx, part in enumerate(parts):
        if list(part.var_gene_names) != ref_genes:
            raise ValueError(f"Gene axis mismatch across parts: {part.tag} != {parts[0].tag}")
        part_bcs = getattr(part, "cell_barcodes", None)
        for loc_i in np.where(part.nt_mask)[0]:
            li = int(loc_i)
            ntc_refs.append((p_idx, li))
            ntc_batches.append(str(part.batches[li]))
            if part_bcs is not None:
                ntc_barcodes.append(str(part_bcs[li]))
            else:
                ntc_barcodes.append(f"{part.tag}:row_{li}")

    n_total = len(ntc_refs)
    sets = sample_stratified_control_sets_r16(
        n_total=n_total,
        batch_labels=ntc_batches,
        target_n=target_n,
        cell_set_len=cell_set_len,
        min_ctrl=min_ctrl,
        seed=seed,
    )
    flat_idx = sets.ravel()
    n_sampled = len(flat_idx)
    raw_out = _read_part_rows_by_indices(parts, ntc_refs, flat_idx, len(ref_genes))
    sampled_batches = [ntc_batches[int(i)] for i in flat_idx]
    sampled_barcodes = [ntc_barcodes[int(i)] for i in flat_idx]

    resamp_raw: Optional[np.ndarray] = None
    resamp_batches: Optional[List[str]] = None
    resamp_barcodes: Optional[List[str]] = None
    resamp_audit: Dict[str, Any] = {"eligible": False, "reason": "not_requested"}
    if include_resample:
        resamp_sets, resamp_audit = sample_disjoint_resample_control_sets_r16(
            n_total=n_total,
            primary_sets=sets,
            batch_labels=ntc_batches,
            target_n=target_n,
            cell_set_len=cell_set_len,
            min_ctrl=min_ctrl,
            seed=resample_seed,
        )
        if resamp_sets is not None:
            r_flat = resamp_sets.ravel()
            resamp_raw = _read_part_rows_by_indices(parts, ntc_refs, r_flat, len(ref_genes))
            resamp_batches = [ntc_batches[int(i)] for i in r_flat]
            resamp_barcodes = [ntc_barcodes[int(i)] for i in r_flat]
            resamp_audit["sampled_barcodes"] = resamp_barcodes
            resamp_audit["sampled_batch_labels"] = resamp_batches

    audit = {
        "n_ctrl_total": int(n_total),
        "n_ctrl_sampled": int(n_sampled),
        "n_sets": int(sets.shape[0]),
        "cell_set_len": int(sets.shape[1]),
        "ctrl_sampling_fraction": float(n_sampled / max(n_total, 1)),
        "n_raw_genes": int(len(ref_genes)),
        "seed": int(seed),
        "sampled_barcodes": sampled_barcodes,
        "sampled_batch_labels": sampled_batches,
        "within_line_resample": resamp_audit,
    }
    return {
        "raw_counts": raw_out,
        "gene_names": ref_genes,
        "batch_labels": sampled_batches,
        "barcodes": sampled_barcodes,
        "ctrl_audit": audit,
        "resample_raw_counts": resamp_raw,
        "resample_batch_labels": resamp_batches,
        "resample_barcodes": resamp_barcodes,
        "resample_audit": resamp_audit,
    }


def extract_ntc_raw_counts_from_parts_r16(
    parts: Sequence[Any],
    target_n: int = 512,
    cell_set_len: int = 32,
    min_ctrl: int = 128,
    seed: int = DEFAULT_R16_SEED,
) -> Tuple[np.ndarray, List[str], List[str], Dict[str, Any]]:
    """
    Sample `target_n` (`512`, `16 x 32`) QC-passing NTC cells proportionally stratified across
    batches from a list of `PartSpec` / `LightweightPartSpec` objects using `PREREG_I v5 §2.2`.
    Returns `(raw_counts [n_ctrl, n_genes], gene_names, sampled_batch_labels, ctrl_audit)`
    where `ctrl_audit["sampled_barcodes"]` stores the exact ordered barcode list (Item 7).
    """
    b = extract_ntc_bundle_from_parts_r16(
        parts=parts,
        target_n=target_n,
        cell_set_len=cell_set_len,
        min_ctrl=min_ctrl,
        seed=seed,
        include_resample=False,
    )
    return b["raw_counts"], b["gene_names"], b["batch_labels"], b["ctrl_audit"]


def extract_gse295214_ntc_raw_counts(
    raw_dir: Path,
    hgnc: Any,
    seed: int = DEFAULT_R16_SEED,
    parts_override: Optional[Sequence[Any]] = None,
) -> Tuple[np.ndarray, List[str], List[str], Dict[str, Any]]:
    """Extract stratified `512`-NTC raw integer UMI matrix + barcodes for `GSE295214` (`PREREG_I v5 §2.2, §2.4`)."""
    if parts_override is not None:
        parts = list(parts_override)
        reader_audit: Dict[str, Any] = {}
    else:
        parts, reader_audit = build_gse295214_parts(raw_dir=raw_dir, hgnc=hgnc, use_r10_partspec=False)
    raw_mat, genes, batches, ctrl_audit = extract_ntc_raw_counts_from_parts_r16(parts, seed=seed)
    ctrl_audit["reader_audit"] = reader_audit
    return raw_mat, genes, batches, ctrl_audit


def prepare_gse295214_line_r16(
    raw_dir: Path,
    out_root: Path = DEFAULT_R16_SHARED_DIR / "lines",
    qc_md_path: Path = DEFAULT_R16_SHARED_DIR / "R16_GSE295214_QC_COUNTS.md",
    status_dir: Path = DEFAULT_STATUS_DIR,
    models_dir: Path = DEFAULT_MODELS_DIR,
    workers: int = 8,
    chunk_size: int = 8192,
) -> Dict[str, Any]:
    """
    Post-lock, post-human-approval standardization & COUNTS-ONLY QC gate for Group 3 `GSE295214`
    (`PREREG_I v5 §2.4`). Refuses execution unless `PREREG_I_COMMIT.txt` exists.
    Never computes any model prediction or correlation metric.
    """
    assert_prereg_i_locked(status_dir)
    if str(R10_SCRIPTS_DIR) not in sys.path:
        sys.path.insert(0, str(R10_SCRIPTS_DIR))
    from prepare_r10_lines import build_hgnc_lookup, process_r10_line  # type: ignore

    hgnc_l = build_hgnc_lookup()
    parts, audit = build_gse295214_parts(raw_dir=raw_dir, hgnc=hgnc_l, use_r10_partspec=True)
    line_dir = process_r10_line(
        line_name="gse295214_h1",
        parts=parts,
        hgnc=hgnc_l,
        chunk_size=chunk_size,
        out_root=out_root,
        extra_meta={"adapter_audit": audit},
        n_workers=workers,
    )

    # Counts-only eligibility check (`N_eval^R16 >= 50`, `|G_eval^R9| >= 3000`, `per-half NTC >= 15`)
    bridge = HGNCBridge.from_tsv(DEFAULT_HGNC_FILE)
    common_genes = load_r9_common_genes(hgnc_bridge=bridge, common_genes_file=R9_COMMON_GENES_FILE)
    ds_t = load_standardized_line_r10(line_dir, representation="log1p", common_genes=common_genes, hgnc_bridge=bridge, min_gene_overlap=1)
    n_g_r9 = len(ds_t.genes)

    src_targets: Set[str] = set()
    for sn in GROUP3_PRIMARY_FCT_SOURCES:
        s_ds = load_standardized_line_r10(R9_LINES_DIR / sn, representation="log1p", common_genes=common_genes, hgnc_bridge=bridge)
        if s_ds.repr_mask is not None:
            src_targets.update(t for t, ok in zip(s_ds.target_names, s_ds.repr_mask) if ok)
        else:
            src_targets.update(s_ds.target_names)

    p5_targets = [
        t
        for i, t in enumerate(ds_t.target_names)
        if ds_t.target_indices[i] >= 0
        and (ds_t.repr_mask is None or bool(ds_t.repr_mask[i]))
        and t in src_targets
    ]
    g_model, p_model_canon_map = load_state_vocab_r16(models_dir, bridge)
    g_model_canon = {bridge.canonicalize(g) or g for g in g_model}
    n_g_r16 = sum(1 for g in ds_t.genes if g in g_model_canon)
    n_eval_r16 = sum(1 for t in p5_targets if t in p_model_canon_map)

    passes = bool(n_eval_r16 >= 50 and n_g_r9 >= 3000 and audit["per_half_ntc_ge15"])
    qc_summary = {
        "line": "gse295214_h1",
        "qualifies_group3_primary": passes,
        "n_eval_p5_r9_excl_vcc2025_h1": len(p5_targets),
        "n_eval_r16": int(n_eval_r16),
        "n_genes_g_eval_r9": int(n_g_r9),
        "n_genes_g_eval_r16": int(n_g_r16),
        "n_kept_half_a_ntc": int(audit["n_kept_half_a_ntc"]),
        "n_kept_half_b_ntc": int(audit["n_kept_half_b_ntc"]),
        "per_half_ntc_ge15": bool(audit["per_half_ntc_ge15"]),
        "adapter_audit": audit,
    }
    qc_md_path.parent.mkdir(parents=True, exist_ok=True)
    qc_md_path.write_text(
        "# Round 16 Group 3 (`GSE295214`, `H1 hESC`) Counts-Only QC Gate\n\n"
        f"- **Status:** `{'PASS_GROUP3_PRIMARY' if passes else 'SENSITIVITY_ONLY'}`\n"
        f"- **`N_eval^R16` (`r9_6 \\ {{vcc2025_h1}} ∩ P_model`):** `{n_eval_r16}` (`>= 50` required)\n"
        f"- **`|G_eval^R9|`:** `{n_g_r9}` (`>= 3,000` required); **`|G_eval^R16|`:** `{n_g_r16}`\n"
        f"- **Per-Half NTC (`non-targeting_*`):** Half A = `{audit['n_kept_half_a_ntc']}`, Half B = `{audit['n_kept_half_b_ntc']}` (`>= 15` each required)\n\n"
        "```json\n" + json.dumps(qc_summary, indent=2) + "\n```\n",
        encoding="utf-8",
    )
    return qc_summary


# ==============================================================================
# Local Raw-Count NTC & Barcode Extractors (Read-Only, Zero Writes to R9–R15)
# ==============================================================================


def _read_h5ad_obs_col_r16(obs_grp: h5py.Group, col: str) -> List[str]:
    node = obs_grp[col]
    if isinstance(node, h5py.Group) and "codes" in node and "categories" in node:
        cats = decode_str_array(node["categories"][:])
        codes = np.asarray(node["codes"][:], dtype=np.int64)
        return [cats[c] if 0 <= c < len(cats) else "" for c in codes]
    return decode_str_array(node[:])


def _resolve_h5ad_obs_index_col_r16(group: h5py.Group) -> Optional[str]:
    """Resolve AnnData `.obs`/`.var` index column cleanly (`attrs['_index']` -> `_index` -> `index`)."""
    keys = set(group.keys())
    attr_idx = group.attrs.get("_index")
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
    for cand in ("_index", "index", "cell_barcode", "barcode", "cell_id"):
        if cand in keys:
            return cand
    return None


def _read_h5ad_x_rows_r16(x_node: Any, sel_raw_rows: np.ndarray, n_genes: int) -> np.ndarray:
    """Read exact `sel_raw_rows` from AnnData `/X` (dense, CSC, or CSR) into float32."""
    raw_mat = np.zeros((len(sel_raw_rows), n_genes), dtype=np.float32)
    if len(sel_raw_rows) == 0:
        return raw_mat
    if isinstance(x_node, h5py.Dataset):
        uniq_r, inv_r = np.unique(sel_raw_rows, return_inverse=True)
        uniq_mat = np.asarray(x_node[uniq_r, :], dtype=np.float32)
        raw_mat = uniq_mat[inv_r]
    else:
        enc = str(x_node.attrs.get("encoding-type", "")).lower()
        if "csc" in enc:
            shape = tuple(int(x) for x in x_node.attrs["shape"])
            mat = sp.csc_matrix(
                (x_node["data"][:], x_node["indices"][:], x_node["indptr"][:]), shape=shape
            ).tocsr()
            raw_mat = mat[sel_raw_rows].toarray().astype(np.float32, copy=False)
        else:
            indptr = x_node["indptr"][:]
            data_ds = x_node["data"]
            idx_ds = x_node["indices"]
            for out_i, r_i in enumerate(sel_raw_rows):
                p0, p1 = int(indptr[int(r_i)]), int(indptr[int(r_i) + 1])
                if p1 > p0:
                    raw_mat[out_i, idx_ds[p0:p1]] = np.asarray(data_ds[p0:p1], dtype=np.float32)
    return raw_mat


def extract_nadig_h5ad_ntc_bundle_r16(
    h5ad_path: Path,
    seed: int = DEFAULT_R16_SEED,
    resample_seed: int = DEFAULT_RESAMPLE_SEED,
    target_n: int = 512,
    include_resample: bool = True,
) -> Dict[str, Any]:
    """
    Bounded-memory pure-`h5py` extractor reading ONLY the sampled primary (`seed = 20261601`)
    and optional disjoint within-line resample (`seed = 20261602`) `non-targeting` NTC rows,
    their native `gem_group` labels, and exact AnnData `.obs` cell barcodes from Nadig 2025 `.h5ad`.
    """
    with h5py.File(h5ad_path, "r") as f:
        obs = f["obs"]
        obs_keys = set(obs.keys())
        gene_col = next((c for c in ("gene", "target_gene", "perturbation") if c in obs_keys), None)
        if gene_col is None:
            raise KeyError(f"Nadig h5ad {h5ad_path} missing gene/perturbation obs column: {sorted(obs_keys)}")
        pert_vals = _read_h5ad_obs_col_r16(obs, gene_col)
        gem_col = next((c for c in ("gem_group", "batch", "channel") if c in obs_keys), None)
        gem_vals = _read_h5ad_obs_col_r16(obs, gem_col) if gem_col else ["1"] * len(pert_vals)
        obs_idx_col = _resolve_h5ad_obs_index_col_r16(obs)
        obs_bcs = (
            _read_h5ad_obs_col_r16(obs, obs_idx_col)
            if obs_idx_col is not None
            else [f"obs_{i}" for i in range(len(pert_vals))]
        )

        var = f["var"]
        var_keys = set(var.keys())
        vname_col = next(
            (c for c in ("gene_name", "gene_symbol", "feature_name", "gene_names", "symbol") if c in var_keys),
            None,
        )
        if vname_col is None:
            vname_col = _resolve_h5ad_obs_index_col_r16(var)
        if vname_col is None:
            raise KeyError(f"Nadig h5ad {h5ad_path} missing var gene-name column: {sorted(var_keys)}")
        var_genes = _read_h5ad_obs_col_r16(var, vname_col)

        ntc_rows = np.array(
            [i for i, v in enumerate(pert_vals) if str(v).strip().lower() in ("non-targeting", "nontargeting")],
            dtype=np.int64,
        )
        ntc_gems = [str(gem_vals[int(i)]).strip() for i in ntc_rows]
        ntc_bcs = [str(obs_bcs[int(i)]).strip() for i in ntc_rows]

        sets = sample_stratified_control_sets_r16(
            n_total=len(ntc_rows),
            batch_labels=ntc_gems,
            target_n=target_n,
            cell_set_len=32,
            min_ctrl=128,
            seed=seed,
        )
        flat_loc = sets.ravel()
        sel_raw_rows = ntc_rows[flat_loc]
        sel_gems = [ntc_gems[int(i)] for i in flat_loc]
        sel_barcodes = [ntc_bcs[int(i)] for i in flat_loc]

        x_node = f["X"]
        n_genes = len(var_genes)
        raw_mat = _read_h5ad_x_rows_r16(x_node, sel_raw_rows, n_genes)

        resamp_raw: Optional[np.ndarray] = None
        resamp_gems: Optional[List[str]] = None
        resamp_barcodes: Optional[List[str]] = None
        resamp_audit: Dict[str, Any] = {"eligible": False, "reason": "not_requested"}
        if include_resample:
            resamp_sets, resamp_audit = sample_disjoint_resample_control_sets_r16(
                n_total=len(ntc_rows),
                primary_sets=sets,
                batch_labels=ntc_gems,
                target_n=target_n,
                cell_set_len=32,
                min_ctrl=128,
                seed=resample_seed,
            )
            if resamp_sets is not None:
                r_loc = resamp_sets.ravel()
                r_raw_rows = ntc_rows[r_loc]
                resamp_raw = _read_h5ad_x_rows_r16(x_node, r_raw_rows, n_genes)
                resamp_gems = [ntc_gems[int(i)] for i in r_loc]
                resamp_barcodes = [ntc_bcs[int(i)] for i in r_loc]
                resamp_audit["sampled_barcodes"] = resamp_barcodes
                resamp_audit["sampled_batch_labels"] = resamp_gems

    audit = {
        "source_file": str(h5ad_path),
        "n_ctrl_total": int(len(ntc_rows)),
        "n_ctrl_sampled": int(len(sel_raw_rows)),
        "n_sets": int(sets.shape[0]),
        "cell_set_len": int(sets.shape[1]),
        "ctrl_sampling_fraction": float(len(sel_raw_rows) / max(len(ntc_rows), 1)),
        "n_raw_genes": int(len(var_genes)),
        "n_distinct_gem_groups_sampled": int(len(set(sel_gems))),
        "seed": int(seed),
        "sampled_barcodes": sel_barcodes,
        "sampled_batch_labels": sel_gems,
        "within_line_resample": resamp_audit,
    }
    return {
        "raw_counts": raw_mat,
        "gene_names": var_genes,
        "batch_labels": sel_gems,
        "barcodes": sel_barcodes,
        "ctrl_audit": audit,
        "resample_raw_counts": resamp_raw,
        "resample_batch_labels": resamp_gems,
        "resample_barcodes": resamp_barcodes,
        "resample_audit": resamp_audit,
    }


def extract_nadig_h5ad_ntc_raw_counts(
    h5ad_path: Path,
    seed: int = DEFAULT_R16_SEED,
    target_n: int = 512,
) -> Tuple[np.ndarray, List[str], List[str], Dict[str, Any]]:
    """Backward-compatible wrapper returning `(raw_mat, var_genes, sel_gems, audit)` with barcodes in `audit`."""
    b = extract_nadig_h5ad_ntc_bundle_r16(
        h5ad_path=h5ad_path, seed=seed, target_n=target_n, include_resample=False
    )
    return b["raw_counts"], b["gene_names"], b["batch_labels"], b["ctrl_audit"]


def _pick_col_r16(columns: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
    low_map = {str(c).lower(): str(c) for c in columns}
    for cand in candidates:
        if cand in columns:
            return cand
        if cand.lower() in low_map:
            return low_map[cand.lower()]
    return None


def _find_file_r16(directory: Path, preferred: Sequence[str], glob_pattern: str) -> Path:
    for name in preferred:
        cand = directory / name
        if cand.exists():
            return cand
    matches = sorted(directory.glob(glob_pattern))
    if not matches:
        raise FileNotFoundError(f"Missing expected file matching {glob_pattern} in {directory}")
    return matches[0]


def _read_mtx_header_r16(mtx_path: Path) -> Tuple[int, int, int, int, str]:
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
                n_skip -= 1
                break
    if dims_line is None:
        raise RuntimeError(f"Empty or malformed MTX header in {mtx_path}")
    parts = dims_line.split()
    if len(parts) < 3:
        raise RuntimeError(f"Invalid MTX dimension line {dims_line!r} in {mtx_path}")
    sep = "\t" if (first_data is not None and "\t" in first_data) else " "
    return int(parts[0]), int(parts[1]), int(parts[2]), n_skip, sep


def _parse_feng_gene_row_label_r16(raw_label: str) -> Tuple[str, str, bool]:
    """Exact mirror of `dataset_readers_r11._parse_feng_gene_row_label`."""
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


def _build_dld1_ntc_lightweight_parts_r16(raw_dir: Path) -> List[LightweightPartSpec]:
    """
    Read-only exact mirror of `prepare_r10_lines.build_dld1_parts` for QC-passing NTC cells,
    preserving exact `f"{gsm}:{raw_barcode}"` barcodes with ZERO writes to `round10/`.
    """
    ddir = Path(raw_dir) / "dld1" if (Path(raw_dir) / "dld1").is_dir() else Path(raw_dir)
    rna_files = sorted(ddir.glob("GSM98642*.filtered_feature_bc_matrix.h5"))
    if not rna_files:
        raise FileNotFoundError(f"No DLD-1 Low-MOI RNA .h5 files found in {ddir}")

    parts: List[LightweightPartSpec] = []
    for ch_idx, r_path in enumerate(rna_files, start=1):
        gsm = r_path.name.split("_")[0]
        h_path = sorted(ddir.glob(f"{gsm}_*.hashing.umi_correct.h5"))[0]
        c_path = sorted(ddir.glob(f"{gsm}_*.crispr.umi_correct.h5"))[0]

        with h5py.File(r_path, "r") as fr:
            rm = fr["matrix"]
            n_genes_raw = int(rm["shape"][0])
            raw_bcs_full = decode_str_array(rm["barcodes"][:])
            rna_bcs = [x.split("-")[0] for x in raw_bcs_full]
            var_ids = decode_str_array(rm["features/id"][:])
            var_names = decode_str_array(rm["features/name"][:])
        mt_raw = np.array([g.startswith("MT-") for g in var_names], dtype=bool)
        rna_umi, rna_ng, rna_mt = _compute_10x_csc_rna_qc_r16(r_path, None, mt_raw)
        qc_rna = (rna_mt < 20.0) & (rna_ng >= 500)

        with h5py.File(h_path, "r") as fh:
            hm = fh["matrix"]
            hbcs = [x.split("-")[0] for x in decode_str_array(hm["barcodes"][:])]
            cmo_names = decode_str_array(hm["features/name"][:])
            H_csc = sp.csc_matrix((hm["data"][:], hm["indices"][:], hm["indptr"][:]), shape=tuple(hm["shape"][:]))
        with h5py.File(c_path, "r") as fc:
            cm = fc["matrix"]
            cbcs = [x.split("-")[0] for x in decode_str_array(cm["barcodes"][:])]
            gnames = decode_str_array(cm["features/name"][:])
            C_csc = sp.csc_matrix((cm["data"][:], cm["indices"][:], cm["indptr"][:]), shape=tuple(cm["shape"][:]))

        h_map = {b: j for j, b in enumerate(hbcs)}
        c_map = {b: j for j, b in enumerate(cbcs)}
        matched_r_idx = [i for i, b in enumerate(rna_bcs) if qc_rna[i] and b in h_map and b in c_map]
        matched_h_cols = [h_map[rna_bcs[i]] for i in matched_r_idx]
        matched_c_cols = [c_map[rna_bcs[i]] for i in matched_r_idx]

        H = H_csc[:, matched_h_cols].T.toarray()
        Hs = np.sort(H, axis=1)
        h_arg = H.argmax(axis=1)
        is_cmo_ok = np.isin([cmo_names[a] for a in h_arg], sorted(_DLD1_LOW_MOI_CMOS))
        sing_h = (Hs[:, -1] >= 10) & (Hs[:, -1] / np.maximum(Hs[:, -2], 0.5) >= 3.0) & is_cmo_ok

        C_sub = C_csc[:, np.asarray(matched_c_cols, dtype=np.int64)[sing_h]]
        c_t1, c_t2, c_a1 = _extract_top2_csc_r16(C_sub)
        sing_c = (c_t1 >= 5.0) & (c_t2 < 3.0) & (c_t1 / np.maximum(c_t2, 0.5) >= 3.0)

        h_pass_r = np.asarray(matched_r_idx, dtype=np.int64)[sing_h]
        keep_cols: List[int] = []
        for loc_j, r_i in enumerate(h_pass_r):
            if not sing_c[loc_j]:
                continue
            gname = gnames[int(c_a1[loc_j])]
            if gname.startswith("NTC_MENG"):
                continue
            if gname.startswith("NTC_"):
                keep_cols.append(int(r_i))

        keep_arr = np.sort(np.asarray(keep_cols, dtype=np.int64))
        bcs_kept = [f"{gsm}:{raw_bcs_full[int(c)]}" for c in keep_arr]
        read_chunk, read_rows, read_chunk_cols = _make_10x_csc_subset_readers_r16(
            r_path, keep_arr, n_genes_raw, None
        )
        parts.append(
            LightweightPartSpec(
                tag=gsm,
                path=r_path,
                n_cells=len(keep_arr),
                n_genes_raw=n_genes_raw,
                nt_mask=np.ones(len(keep_arr), dtype=bool),
                batches=np.full(len(keep_arr), f"ch_{ch_idx:02d}"),
                half_flags=np.full(len(keep_arr), 1 if (ch_idx % 2 == 1) else 2, dtype=np.int8),
                umi=rna_umi[keep_arr],
                pert_keys=[""] * len(keep_arr),
                pert_meta={},
                var_gene_ids=var_ids,
                var_gene_names=var_names,
                read_chunk_fn=read_chunk,
                read_rows_fn=read_rows,
                read_chunk_cols_fn=read_chunk_cols,
                cell_barcodes=bcs_kept,
            )
        )
    return parts


def _build_cd4t_rest_ntc_lightweight_parts_r16(raw_dir: Path) -> List[LightweightPartSpec]:
    """
    Read-only exact mirror of `prepare_r10_lines.build_cd4t_rest_parts` for QC-passing NTC cells,
    preserving exact `f"{donor_id}_{lane_tag}:{raw_barcode}"` barcodes with ZERO writes to `round10/`.
    """
    cdir = Path(raw_dir) / "cd4t_rest" if (Path(raw_dir) / "cd4t_rest").is_dir() else Path(raw_dir)
    ntc_guide_set: Set[str] = set()
    with open(cdir / "sgrna_library_metadata.suppl_table.csv", "r", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh):
            g = str(row.get("sgRNA", "")).strip()
            if g.startswith("NTC-"):
                ntc_guide_set.add(g)

    d4_files = sorted(cdir.glob("GSM93938*_D4_Rest_sample_filtered_feature_bc_matrix.h5"))
    d1_files = sorted(cdir.glob("GSM93938*_D1_Rest_sample_filtered_feature_bc_matrix.h5"))
    if not d4_files or not d1_files:
        raise FileNotFoundError(f"Missing cd4t_rest D4/D1 .h5 files in {cdir}")

    parts: List[LightweightPartSpec] = []
    for donor_id, half_val, files in [
        ("CE0006864_D4", 1, d4_files),
        ("CE0008162_D1", 2, d1_files),
    ]:
        lane_cache = []
        hist: Optional[np.ndarray] = None
        crp_names_ref: Optional[List[str]] = None
        for lp in files:
            parts_name = lp.name.split("_")
            lane_tag = f"{parts_name[1]}_{parts_name[2]}"
            with h5py.File(lp, "r") as f:
                m = f["matrix"]
                n_feat, n_col = int(m["shape"][0]), int(m["shape"][1])
                raw_bcs = decode_str_array(m["barcodes"][:])
                ftypes = np.array(decode_str_array(m["features/feature_type"][:]))
                fids = np.array(decode_str_array(m["features/id"][:]))
                fnames = np.array(decode_str_array(m["features/name"][:]))
                rna_rows = np.where(ftypes == "Gene Expression")[0]
                crp_rows = np.where(ftypes != "Gene Expression")[0]
                if crp_names_ref is None:
                    crp_names_ref = [_CD4_EXCEL_GUIDE_REMAP.get(x, x) for x in fnames[crp_rows]]
                    hist = np.zeros((len(crp_rows), 66), dtype=np.int64)
                mt_raw = np.zeros(n_feat, dtype=bool)
                mt_raw[rna_rows] = [str(x).startswith("MT-") for x in fnames[rna_rows]]

                is_crp = np.zeros(n_feat, dtype=bool)
                is_crp[crp_rows] = True
                crp_remap_idx = np.full(n_feat, -1, dtype=np.int32)
                crp_remap_idx[crp_rows] = np.arange(len(crp_rows), dtype=np.int32)
                is_rna_feat = ~is_crp

                umi_l = np.zeros(n_col, dtype=np.float64)
                ng_l = np.zeros(n_col, dtype=np.int32)
                mt_umi_l = np.zeros(n_col, dtype=np.float64)
                crp_idx_parts, crp_dat_parts, crp_col_parts = [], [], []

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
            if c_col.size > 0 and hist is not None:
                ok_init_mask = ok_init[c_col]
                if np.any(ok_init_mask):
                    np.add.at(hist, (c_idx[ok_init_mask], np.clip(c_dat[ok_init_mask], 1, 65)), 1)

            hq_cols = np.where(ok_hq)[0]
            col_map_hq = np.full(n_col, -1, dtype=np.int32)
            col_map_hq[hq_cols] = np.arange(len(hq_cols), dtype=np.int32)
            if c_col.size > 0 and len(hq_cols) > 0:
                hq_m = col_map_hq[c_col] >= 0
                C_hq = sp.csr_matrix(
                    (c_dat[hq_m].astype(np.float64), (c_idx[hq_m], col_map_hq[c_col[hq_m]])),
                    shape=(len(crp_rows), len(hq_cols)),
                ).tocsc()
            else:
                C_hq = sp.csc_matrix((len(crp_rows), len(hq_cols)), dtype=np.float64)
            lane_cache.append(
                (
                    lp,
                    lane_tag,
                    n_feat,
                    rna_rows,
                    [str(fids[i]) for i in rna_rows],
                    [str(fnames[i]) for i in rna_rows],
                    raw_bcs,
                    umi_l,
                    hq_cols,
                    C_hq,
                )
            )

        assert hist is not None and crp_names_ref is not None
        thrs = _fit_poisson_gaussian_thresholds_r16(hist, min_thr=5, max_thr=20, default_thr=7)
        for lp, lane_tag, n_feat, rna_rows, rids, rnames, raw_bcs, umi_l, hq_cols, C_hq in lane_cache:
            t1, t2, a1 = _extract_top2_csc_r16(C_hq)
            pos_mask = C_hq.data >= thrs[C_hq.indices]
            pos_csc = sp.csc_matrix((pos_mask.astype(np.int8), C_hq.indices, C_hq.indptr), shape=C_hq.shape)
            n_pos = np.asarray(pos_csc.sum(axis=0)).ravel()
            s1_mask = (n_pos == 1) & (t1 >= 5.0) & ((t2 < 3.0) | (t1 / np.maximum(t2, 0.5) >= 3.0))

            keep_cols = [
                int(hq_cols[loc_j])
                for loc_j in np.where(s1_mask)[0]
                if crp_names_ref[int(a1[loc_j])] in ntc_guide_set
            ]
            keep_arr = np.asarray(keep_cols, dtype=np.int64)
            bcs_kept = [f"{donor_id}_{lane_tag}:{raw_bcs[int(c)]}" for c in keep_arr]
            read_chunk, read_rows, read_chunk_cols = _make_10x_csc_subset_readers_r16(
                lp, keep_arr, n_feat, rna_rows
            )
            parts.append(
                LightweightPartSpec(
                    tag=lane_tag,
                    path=lp,
                    n_cells=len(keep_arr),
                    n_genes_raw=len(rna_rows),
                    nt_mask=np.ones(len(keep_arr), dtype=bool),
                    batches=np.full(len(keep_arr), f"{donor_id}_{lane_tag}"),
                    half_flags=np.full(len(keep_arr), half_val, dtype=np.int8),
                    umi=umi_l[keep_arr],
                    pert_keys=[""] * len(keep_arr),
                    pert_meta={},
                    var_gene_ids=rids,
                    var_gene_names=rnames,
                    read_chunk_fn=read_chunk,
                    read_rows_fn=read_rows,
                    read_chunk_cols_fn=read_chunk_cols,
                    cell_barcodes=bcs_kept,
                )
            )
    return parts


def _extract_a549_ntc_bundle_r16(
    raw_dir: Path,
    seed: int = DEFAULT_R16_SEED,
    resample_seed: int = DEFAULT_RESAMPLE_SEED,
    target_n: int = 512,
    include_resample: bool = True,
    bucket_size: int = 40000,
    chunk_lines: int = 1_500_000,
    min_batch_nt: int = 15,
    max_pct_mt: float = 20.0,
    min_ngenes: int = 500,
) -> Dict[str, Any]:
    """
    Read-only exact mirror of `dataset_readers_r12.build_a549_parts` strictly on NTC cells,
    streaming only NTC cells into a `TemporaryDirectory` (ZERO writes to `round12/shared/raw/`)
    and preserving exact `cell_id` barcodes (Addendum A Items 7 & 8).
    """
    import pandas as pd

    base = Path(raw_dir)
    adir = base / "a549_ko" if (base / "a549_ko").is_dir() else (base / "a549" if (base / "a549").is_dir() else base)
    cell_path = _find_file_r16(
        adir, ["GSE345058_cell_annotations.csv.gz", "GSE345058_cell_annotations.csv"], "*cell_annotations*.csv*"
    )
    gene_path = _find_file_r16(
        adir, ["GSE345058_gene_annotations.csv.gz", "GSE345058_gene_annotations.csv"], "*gene_annotations*.csv*"
    )
    mtx_path = _find_file_r16(
        adir, ["GSE345058_raw_counts.mtx.gz", "GSE345058_raw_counts.mtx"], "*raw_counts*.mtx*"
    )

    cell_df = pd.read_csv(cell_path, dtype=str, keep_default_na=False)
    gene_df = pd.read_csv(gene_path, dtype=str, keep_default_na=False)
    n_annot_cells = len(cell_df)
    n_annot_genes = len(gene_df)

    c_cellid_col = _pick_col_r16(cell_df.columns, ["cell_id", "Cell_ID", "barcode", "Barcode", "cell"])
    c_batch_col = _pick_col_r16(cell_df.columns, ["batch", "Batch", "batch_id"])
    c_unit_col = _pick_col_r16(cell_df.columns, ["Lane", "lane", "channel", "sample", "batch", "Batch"])
    c_sgrna_col = _pick_col_r16(cell_df.columns, ["sgRNA", "sgrna", "guide_id", "guide"])
    c_pert_col = _pick_col_r16(cell_df.columns, ["Perturbation", "perturbation", "target_gene", "gene_symbol", "gene"])
    c_cond_col = _pick_col_r16(cell_df.columns, ["Condition", "condition"])
    if c_batch_col is None or c_pert_col is None or c_sgrna_col is None:
        raise KeyError(f"A549 cell annotations missing required columns in {cell_path}")

    g_sym_col = _pick_col_r16(gene_df.columns, ["gene_symbol", "gene_name", "symbol", "Symbol", "name", "gene"]) or gene_df.columns[0]
    g_ens_col = _pick_col_r16(gene_df.columns, ["ensembl_id", "gene_id", "gene_ids", "Ensembl", "ensembl"])
    g_type_col = _pick_col_r16(gene_df.columns, ["feature_types", "feature_type", "type"])
    g_gen_col = _pick_col_r16(gene_df.columns, ["genome", "species", "reference"])
    raw_g_syms = [str(x).strip() for x in gene_df[g_sym_col].tolist()]
    raw_g_ens = [str(x).strip().split(".")[0] for x in gene_df[g_ens_col].tolist()] if g_ens_col else [""] * n_annot_genes
    raw_g_gen = [str(x).strip().upper() for x in gene_df[g_gen_col].tolist()] if g_gen_col else [""] * n_annot_genes

    def _is_cas9(i_g: int) -> bool:
        return raw_g_syms[i_g].upper() == "CAS9" or raw_g_ens[i_g].upper() == "CAS9" or raw_g_gen[i_g] == "CAS9"

    if g_type_col is not None:
        g_types = [str(x).strip() for x in gene_df[g_type_col].tolist()]
        rna_gene_raw_idx = [i for i, t in enumerate(g_types) if t in ("Gene Expression", "RNA", "") and not _is_cas9(i)]
        if not rna_gene_raw_idx:
            rna_gene_raw_idx = [i for i in range(n_annot_genes) if not _is_cas9(i)]
    else:
        rna_gene_raw_idx = [i for i in range(n_annot_genes) if not _is_cas9(i)]

    raw_to_rna = np.full(n_annot_genes, -1, dtype=np.int32)
    rna_names: List[str] = []
    for k, ri in enumerate(rna_gene_raw_idx):
        raw_to_rna[ri] = k
        rna_names.append(raw_g_syms[ri])
    n_genes_rna = len(rna_names)
    mt_mask_rna = np.array([gn.startswith("MT-") for gn in rna_names], dtype=bool)

    cids_raw = (
        cell_df[c_cellid_col].astype(str).str.strip().tolist()
        if c_cellid_col is not None
        else [f"a549_row_{i}" for i in range(n_annot_cells)]
    )
    batches_raw = cell_df[c_batch_col].astype(str).str.strip().tolist()
    units_raw = cell_df[c_unit_col].astype(str).str.strip().tolist() if c_unit_col else batches_raw
    sgrnas_raw = cell_df[c_sgrna_col].astype(str).str.strip().tolist()
    perts_raw = cell_df[c_pert_col].astype(str).str.strip().tolist()
    conds_raw = cell_df[c_cond_col].astype(str).str.strip().tolist() if c_cond_col is not None else None

    raw_to_ntc = np.full(n_annot_cells, -1, dtype=np.int32)
    ntc_pre_cids: List[str] = []
    ntc_pre_units: List[str] = []
    ntc_pre_halves: List[int] = []
    for ci in range(n_annot_cells):
        sg = sgrnas_raw[ci]
        pr = perts_raw[ci]
        if not sg or sg.lower() in _A549_UNASSIGNED_TOKENS or _MULTI_GUIDE_DELIM_RE.search(sg) or " " in sg:
            continue
        if _MULTI_GUIDE_DELIM_RE.search(pr):
            continue
        m_b = re.search(r"(\d+)$", batches_raw[ci])
        b_num = int(m_b.group(1)) if m_b else None
        if b_num not in (1, 2, 3, 4):
            continue
        half_flag = 1 if b_num in (1, 3) else 2
        pert_is_ntc = pr.lower() in _A549_NTC_TOKENS
        if conds_raw is not None:
            cond_low = conds_raw[ci].lower()
            if cond_low in _A549_UNASSIGNED_TOKENS:
                continue
            cond_is_ntc = cond_low in _A549_NTC_TOKENS
            if pert_is_ntc != cond_is_ntc:
                continue
            is_ntc = pert_is_ntc and cond_is_ntc
        else:
            is_ntc = pert_is_ntc
        if is_ntc:
            raw_to_ntc[ci] = len(ntc_pre_cids)
            ntc_pre_cids.append(cids_raw[ci])
            ntc_pre_units.append(units_raw[ci])
            ntc_pre_halves.append(half_flag)

    n_pre_ntc = len(ntc_pre_cids)
    if n_pre_ntc == 0:
        raise RuntimeError(f"A549 metadata yielded zero NTC cells in {adir}")

    n_mtx_rows, n_mtx_cols, _, n_skip_lines, sep = _read_mtx_header_r16(mtx_path)
    if n_mtx_rows == n_annot_cells and n_mtx_cols == n_annot_genes:
        cell_first = True
    elif n_mtx_rows == n_annot_genes and n_mtx_cols == n_annot_cells:
        cell_first = False
    else:
        raise RuntimeError(f"A549 MTX shape ({n_mtx_rows}, {n_mtx_cols}) incompatible with annotations in {mtx_path}")

    n_buckets = max(1, (n_pre_ntc + bucket_size - 1) // bucket_size)
    with tempfile.TemporaryDirectory(prefix="r16_a549_ntc_") as tmp_dir:
        tmp_path = Path(tmp_dir)
        tmp_h5 = tmp_path / "a549_ntc_csc.h5"
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
                ntc_ci = raw_to_ntc[c_raw]
                g_rna = raw_to_rna[g_raw]
                m_ok = (ntc_ci >= 0) & (g_rna >= 0) & (vv != 0.0)
                if not np.any(m_ok):
                    continue
                k_arr = ntc_ci[m_ok].astype(np.int32)
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

        umi_pre = np.zeros(n_pre_ntc, dtype=np.float64)
        ng_pre = np.zeros(n_pre_ntc, dtype=np.int32)
        mt_umi_pre = np.zeros(n_pre_ntc, dtype=np.float64)
        with h5py.File(tmp_h5, "w") as hf:
            mg = hf.create_group("matrix")
            mg.create_dataset("shape", data=np.array([n_genes_rna, n_pre_ntc], dtype=np.int64))
            data_ds = mg.create_dataset("data", shape=(0,), maxshape=(None,), dtype=np.float32, chunks=(65536,))
            idx_ds = mg.create_dataset("indices", shape=(0,), maxshape=(None,), dtype=np.int32, chunks=(65536,))
            indptr_ds = mg.create_dataset("indptr", shape=(n_pre_ntc + 1,), dtype=np.int64)
            indptr_ds[0] = 0
            cur_nnz = 0
            for bi in range(n_buckets):
                cs = bi * bucket_size
                ce = min(cs + bucket_size, n_pre_ntc)
                nc_b = ce - cs
                ga = np.fromfile(tmp_path / f"b_{bi}.g.bin", dtype=np.int32)
                ca = np.fromfile(tmp_path / f"b_{bi}.c.bin", dtype=np.int32)
                va = np.fromfile(tmp_path / f"b_{bi}.v.bin", dtype=np.float32)
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
                    indptr_ds[cs + 1 : ce + 1] = cur_nnz + csc_b.indptr[1:].astype(np.int64)
                    cur_nnz += add_nnz
                else:
                    indptr_ds[cs + 1 : ce + 1] = cur_nnz

        pct_mt_pre = 100.0 * mt_umi_pre / np.maximum(umi_pre, 1.0)
        qc_ok = (umi_pre > 0.0) & (ng_pre >= int(min_ngenes)) & (pct_mt_pre <= float(max_pct_mt))
        keep_idx = np.where(qc_ok)[0].astype(np.int64)
        ntc_cids = [ntc_pre_cids[int(i)] for i in keep_idx]
        ntc_units = np.asarray([ntc_pre_units[int(i)] for i in keep_idx], dtype=str)
        ntc_halves = np.asarray([ntc_pre_halves[int(i)] for i in keep_idx], dtype=np.int8)
        ntc_batches_arr, _ = _apply_min_ntc_batch_fallback_r16(
            ntc_units, np.ones(len(keep_idx), dtype=bool), ntc_halves, min_batch_nt=min_batch_nt
        )
        ntc_batches = [str(b) for b in ntc_batches_arr]

        sets = sample_stratified_control_sets_r16(
            n_total=len(keep_idx),
            batch_labels=ntc_batches,
            target_n=target_n,
            cell_set_len=32,
            min_ctrl=128,
            seed=seed,
        )
        flat_loc = sets.ravel()
        _, read_rows, _ = _make_10x_csc_subset_readers_r16(tmp_h5, keep_idx, n_genes_rna, None)
        raw_mat = np.asarray(read_rows(flat_loc), dtype=np.float32)
        sel_batches = [ntc_batches[int(i)] for i in flat_loc]
        sel_barcodes = [ntc_cids[int(i)] for i in flat_loc]

        resamp_raw: Optional[np.ndarray] = None
        resamp_batches: Optional[List[str]] = None
        resamp_barcodes: Optional[List[str]] = None
        resamp_audit: Dict[str, Any] = {"eligible": False, "reason": "not_requested"}
        if include_resample:
            resamp_sets, resamp_audit = sample_disjoint_resample_control_sets_r16(
                n_total=len(keep_idx),
                primary_sets=sets,
                batch_labels=ntc_batches,
                target_n=target_n,
                cell_set_len=32,
                min_ctrl=128,
                seed=resample_seed,
            )
            if resamp_sets is not None:
                r_loc = resamp_sets.ravel()
                resamp_raw = np.asarray(read_rows(r_loc), dtype=np.float32)
                resamp_batches = [ntc_batches[int(i)] for i in r_loc]
                resamp_barcodes = [ntc_cids[int(i)] for i in r_loc]
                resamp_audit["sampled_barcodes"] = resamp_barcodes
                resamp_audit["sampled_batch_labels"] = resamp_batches

    audit = {
        "source_file": str(mtx_path),
        "n_ctrl_total": int(len(keep_idx)),
        "n_ctrl_sampled": int(len(sel_barcodes)),
        "n_sets": int(sets.shape[0]),
        "cell_set_len": int(sets.shape[1]),
        "ctrl_sampling_fraction": float(len(sel_barcodes) / max(len(keep_idx), 1)),
        "n_raw_genes": int(n_genes_rna),
        "seed": int(seed),
        "sampled_barcodes": sel_barcodes,
        "sampled_batch_labels": sel_batches,
        "within_line_resample": resamp_audit,
    }
    return {
        "raw_counts": raw_mat,
        "gene_names": rna_names,
        "batch_labels": sel_batches,
        "barcodes": sel_barcodes,
        "ctrl_audit": audit,
        "resample_raw_counts": resamp_raw,
        "resample_batch_labels": resamp_batches,
        "resample_barcodes": resamp_barcodes,
        "resample_audit": resamp_audit,
    }


def _extract_feng_ipsc_ntc_bundle_r16(
    raw_dir: Path,
    seed: int = DEFAULT_R16_SEED,
    resample_seed: int = DEFAULT_RESAMPLE_SEED,
    target_n: int = 512,
    include_resample: bool = True,
    min_batch_nt: int = 15,
) -> Dict[str, Any]:
    """
    Read-only exact mirror of `dataset_readers_r11.build_feng_ipsc_parts` strictly on NTC cells:
    samples primary + optional disjoint resample columns BEFORE streaming the CSV once
    (ZERO temp files, ZERO writes to `round11/shared/raw/`, exact `Cell_ID` barcodes).
    """
    base = Path(raw_dir)
    fdir = base / "feng_ipsc" if (base / "feng_ipsc").is_dir() else base
    meta_path = _find_file_r16(
        fdir, ["TargetedScreen_Cell-Metadata.tsv.gz", "TargetedScreen_Cell-Metadata.tsv"], "*Cell-Metadata*.tsv*"
    )
    csv_path = _find_file_r16(
        fdir, ["TargetedScreen_RNA-UMI-Counts.csv.gz", "TargetedScreen_RNA-UMI-Counts.csv"], "*RNA-UMI-Counts*.csv*"
    )

    m_open = gzip.open if str(meta_path).endswith(".gz") else open
    ntc_meta: Dict[str, Tuple[str, str]] = {}
    kept_lines_set: Set[str] = set()
    with m_open(meta_path, "rt", encoding="utf-8", errors="replace") as fm:
        m_reader = csv.DictReader(fm, delimiter="\t")
        cols = m_reader.fieldnames or []
        cid_col = _pick_col_r16(cols, ["Cell_ID", "cell_id", "barcode", "CBC"])
        cline_col = _pick_col_r16(cols, ["Cell_Line", "cell_line", "line"])
        gcall_col = _pick_col_r16(cols, ["Guide_Call", "guide_call", "sgrna", "perturbation"])
        inlet_col = _pick_col_r16(cols, ["Inlet", "inlet", "batch", "lane"])
        if cid_col is None or cline_col is None or gcall_col is None:
            raise KeyError(f"Feng metadata missing required columns in {meta_path}: {cols}")

        for row in m_reader:
            cid = str(row.get(cid_col, "")).strip()
            cline = str(row.get(cline_col, "")).strip()
            gcall = str(row.get(gcall_col, "")).strip()
            if not cid or not cline:
                continue
            cl_low = cline.lower()
            if cl_low in _FENG_KOLF_EXCLUDED or cl_low.startswith("kolf_2") or cl_low.startswith("kolf_3"):
                continue
            low_g = gcall.lower()
            if not gcall or low_g in ("unassigned", "na", "nan", "none"):
                continue
            if _MULTI_GUIDE_DELIM_RE.search(gcall):
                continue
            inlet = str(row.get(inlet_col, "inlet0")).strip() if inlet_col else "inlet0"
            unit_name = f"{cline}__inlet{inlet}"
            kept_lines_set.add(cline)
            if gcall.startswith("NonTarget_") or low_g.startswith("nontarget"):
                ntc_meta[cid] = (cline, unit_name)

    lines_sorted = sorted(kept_lines_set)
    if len(lines_sorted) < 2:
        raise RuntimeError(f"Expected >= 2 non-KOLF cell lines in Feng metadata, got {lines_sorted}")
    line_to_half = {ln: (1 if i % 2 == 0 else 2) for i, ln in enumerate(lines_sorted)}

    c_open = gzip.open if str(csv_path).endswith(".gz") else open
    with c_open(csv_path, "rt", encoding="utf-8", errors="replace") as fc:
        reader = csv.reader(fc)
        header = [h.strip() for h in next(reader)]
        first_data_row = next(reader, None)
        if first_data_row is None:
            raise RuntimeError(f"Empty data rows in {csv_path}")

        col_offset = len(first_data_row) - len(header)
        if col_offset not in (0, 1):
            raise RuntimeError(
                f"Header/row width mismatch in {csv_path}: header={len(header)}, row1={len(first_data_row)}"
            )
        hdr_start = 1 if col_offset == 0 else 0

        ntc_csv_cols: List[int] = []
        ntc_cids: List[str] = []
        ntc_units: List[str] = []
        ntc_lines: List[str] = []
        ntc_halves: List[int] = []
        for h_idx in range(hdr_start, len(header)):
            cid = header[h_idx]
            rec = ntc_meta.get(cid)
            if rec is None:
                continue
            cline, uname = rec
            ntc_csv_cols.append(h_idx + col_offset)
            ntc_cids.append(cid)
            ntc_units.append(uname)
            ntc_lines.append(cline)
            ntc_halves.append(line_to_half[cline])

        n_total = len(ntc_cids)
        ntc_batches_arr, _ = _apply_min_ntc_batch_fallback_r16(
            np.asarray(ntc_units, dtype=str),
            np.ones(n_total, dtype=bool),
            np.asarray(ntc_halves, dtype=np.int8),
            min_batch_nt=min_batch_nt,
            group_labels=np.asarray(ntc_lines, dtype=str),
        )
        ntc_batches = [str(b) for b in ntc_batches_arr]

        sets = sample_stratified_control_sets_r16(
            n_total=n_total,
            batch_labels=ntc_batches,
            target_n=target_n,
            cell_set_len=32,
            min_ctrl=128,
            seed=seed,
        )
        flat_prim = sets.ravel()
        prim_csv_cols = [ntc_csv_cols[int(i)] for i in flat_prim]
        sel_batches = [ntc_batches[int(i)] for i in flat_prim]
        sel_barcodes = [ntc_cids[int(i)] for i in flat_prim]

        resamp_csv_cols: Optional[List[int]] = None
        resamp_batches: Optional[List[str]] = None
        resamp_barcodes: Optional[List[str]] = None
        resamp_audit: Dict[str, Any] = {"eligible": False, "reason": "not_requested"}
        if include_resample:
            resamp_sets, resamp_audit = sample_disjoint_resample_control_sets_r16(
                n_total=n_total,
                primary_sets=sets,
                batch_labels=ntc_batches,
                target_n=target_n,
                cell_set_len=32,
                min_ctrl=128,
                seed=resample_seed,
            )
            if resamp_sets is not None:
                r_flat = resamp_sets.ravel()
                resamp_csv_cols = [ntc_csv_cols[int(i)] for i in r_flat]
                resamp_batches = [ntc_batches[int(i)] for i in r_flat]
                resamp_barcodes = [ntc_cids[int(i)] for i in r_flat]
                resamp_audit["sampled_barcodes"] = resamp_barcodes
                resamp_audit["sampled_batch_labels"] = resamp_batches

        rna_names: List[str] = []
        prim_rows: List[np.ndarray] = []
        resamp_rows: List[np.ndarray] = []

        def _consume(r_vals: List[str]) -> None:
            if not r_vals:
                return
            _, sym, is_gex = _parse_feng_gene_row_label_r16(r_vals[0])
            if not is_gex or not sym:
                return
            rna_names.append(sym)
            prim_rows.append(
                np.fromiter((float(r_vals[j]) for j in prim_csv_cols), dtype=np.float32, count=len(prim_csv_cols))
            )
            if resamp_csv_cols is not None:
                resamp_rows.append(
                    np.fromiter(
                        (float(r_vals[j]) for j in resamp_csv_cols), dtype=np.float32, count=len(resamp_csv_cols)
                    )
                )

        _consume(first_data_row)
        for row in reader:
            _consume(row)

    if not rna_names:
        raise RuntimeError(f"Zero Gene-Expression rows parsed from {csv_path}.")
    raw_mat = np.vstack(prim_rows).T.astype(np.float32, copy=True)
    resamp_raw = np.vstack(resamp_rows).T.astype(np.float32, copy=True) if resamp_csv_cols is not None else None

    audit = {
        "source_file": str(csv_path),
        "n_ctrl_total": int(n_total),
        "n_ctrl_sampled": int(len(sel_barcodes)),
        "n_sets": int(sets.shape[0]),
        "cell_set_len": int(sets.shape[1]),
        "ctrl_sampling_fraction": float(len(sel_barcodes) / max(n_total, 1)),
        "n_raw_genes": int(len(rna_names)),
        "seed": int(seed),
        "sampled_barcodes": sel_barcodes,
        "sampled_batch_labels": sel_batches,
        "within_line_resample": resamp_audit,
    }
    return {
        "raw_counts": raw_mat,
        "gene_names": rna_names,
        "batch_labels": sel_batches,
        "barcodes": sel_barcodes,
        "ctrl_audit": audit,
        "resample_raw_counts": resamp_raw,
        "resample_batch_labels": resamp_batches,
        "resample_barcodes": resamp_barcodes,
        "resample_audit": resamp_audit,
    }


def _resolve_nadig_h5ad_path_r16(
    line_name: str,
    raw_nadig_dir: Optional[Path] = None,
) -> Path:
    """Resolve Nadig 2025 `.h5ad` path via `--raw-nadig-dir`, `NADIG_DATA_DIR`, or sync-drive default."""
    if line_name == "nadig_hepg2":
        fname = "GSE264667_hepg2_raw_singlecell_01.h5ad"
    elif line_name == "nadig_jurkat":
        fname = "GSE264667_jurkat_raw_singlecell_01.h5ad"
    else:
        raise ValueError(f"Not a Nadig line: {line_name!r}")

    if raw_nadig_dir is not None:
        base = Path(raw_nadig_dir)
    elif os.environ.get("NADIG_DATA_DIR"):
        base = Path(os.environ["NADIG_DATA_DIR"])
    else:
        base = SYNC_EXTERNAL_DATA_DIR / "nadig2025"

    direct = base / fname
    if direct.is_file():
        return direct
    nested = base / "nadig2025" / fname
    if nested.is_file():
        return nested
    return direct


def extract_local_ntc_bundle_r16(
    line_name: str,
    status_dir: Path = DEFAULT_STATUS_DIR,
    seed: int = DEFAULT_R16_SEED,
    resample_seed: int = DEFAULT_RESAMPLE_SEED,
    include_resample: bool = True,
    gse295214_raw_dir: Optional[Path] = None,
    raw_nadig_dir: Optional[Path] = None,
) -> Dict[str, Any]:
    """
    Triple-gated, read-only local raw-count NTC + barcode extractor for all 6 local lines
    plus optional Group 3 `gse295214_h1` (primary `seed = 20261601` + Item 8 disjoint resample
    `seed = 20261602`; writes ZERO files into Rounds 9–15).
    """
    assert_all_execution_gates_r16(status_dir)

    if line_name in ("nadig_hepg2", "nadig_jurkat"):
        p = _resolve_nadig_h5ad_path_r16(line_name, raw_nadig_dir=raw_nadig_dir)
        return extract_nadig_h5ad_ntc_bundle_r16(
            p, seed=seed, resample_seed=resample_seed, include_resample=include_resample
        )
    if line_name == "gse295214_h1":
        raw_d = gse295214_raw_dir or (DEFAULT_R16_SHARED_DIR / "raw" / "gse295214_h1")
        bridge = HGNCBridge.from_tsv(DEFAULT_HGNC_FILE)
        parts, reader_audit = build_gse295214_parts(raw_dir=raw_d, hgnc=bridge, use_r10_partspec=False)
        b = extract_ntc_bundle_from_parts_r16(
            parts, seed=seed, resample_seed=resample_seed, include_resample=include_resample
        )
        b["ctrl_audit"]["reader_audit"] = reader_audit
        return b
    if line_name == "a549_ko":
        return _extract_a549_ntc_bundle_r16(
            REPO_ROOT / "round12" / "shared" / "raw",
            seed=seed,
            resample_seed=resample_seed,
            include_resample=include_resample,
        )
    if line_name == "dld1":
        parts = _build_dld1_ntc_lightweight_parts_r16(REPO_ROOT / "round10" / "shared" / "raw")
        return extract_ntc_bundle_from_parts_r16(
            parts, seed=seed, resample_seed=resample_seed, include_resample=include_resample
        )
    if line_name == "cd4t_rest":
        parts = _build_cd4t_rest_ntc_lightweight_parts_r16(REPO_ROOT / "round10" / "shared" / "raw")
        return extract_ntc_bundle_from_parts_r16(
            parts, seed=seed, resample_seed=resample_seed, include_resample=include_resample
        )
    if line_name == "feng_ipsc":
        return _extract_feng_ipsc_ntc_bundle_r16(
            REPO_ROOT / "round11" / "shared" / "raw",
            seed=seed,
            resample_seed=resample_seed,
            include_resample=include_resample,
        )
    raise ValueError(f"Unsupported local raw-count line: {line_name!r}")


def extract_local_ntc_raw_counts(
    line_name: str,
    status_dir: Path = DEFAULT_STATUS_DIR,
    seed: int = DEFAULT_R16_SEED,
    gse295214_raw_dir: Optional[Path] = None,
    raw_nadig_dir: Optional[Path] = None,
) -> Tuple[np.ndarray, List[str], List[str], Dict[str, Any]]:
    """Backward-compatible wrapper returning `(raw_counts, gene_names, batch_labels, ctrl_audit)`."""
    b = extract_local_ntc_bundle_r16(
        line_name=line_name,
        status_dir=status_dir,
        seed=seed,
        include_resample=False,
        gse295214_raw_dir=gse295214_raw_dir,
        raw_nadig_dir=raw_nadig_dir,
    )
    return b["raw_counts"], b["gene_names"], b["batch_labels"], b["ctrl_audit"]


# ==============================================================================
# Arc Archived `adata_pred.h5ad` Adapter Fidelity Check (`Addendum A Item 3`)
# ==============================================================================


def compare_predicted_deltas_vs_arc_archive_r16(
    adapter_deltas: np.ndarray,
    adapter_perts: Sequence[str],
    adapter_genes: Sequence[str],
    arc_deltas: np.ndarray,
    arc_perts: Sequence[str],
    arc_genes: Sequence[str],
    median_pearson_min: float = ARC_ADAPTER_MEDIAN_PEARSON_MIN,
    pooled_pearson_min: float = ARC_ADAPTER_POOLED_PEARSON_MIN,
) -> Dict[str, Any]:
    """
    `PREREG_I Addendum A Item 3`: Compare per-perturbation predicted shifts from our
    `ArcStateAdapter` (`nadig_hepg2`, native `gem_group` batch labels) against Arc's archived
    `adata_pred.h5ad` on all common perturbations and genes.
    Pass criteria (both required):
      - `median_per_pert_pearson >= 0.90`
      - `pooled_pearson >= 0.95`
    Failure halts execution for Master Controller adapter diagnosis without changing rules.
    """
    ad_mat = np.asarray(adapter_deltas, dtype=np.float64)
    ar_mat = np.asarray(arc_deltas, dtype=np.float64)
    arc_p_idx: Dict[str, int] = {}
    for i, p in enumerate(arc_perts):
        arc_p_idx.setdefault(str(p), i)
    arc_g_idx: Dict[str, int] = {}
    for j, g in enumerate(arc_genes):
        arc_g_idx.setdefault(str(g), j)

    common_p_pairs = [(i, arc_p_idx[str(p)], str(p)) for i, p in enumerate(adapter_perts) if str(p) in arc_p_idx]
    common_g_pairs = [(j, arc_g_idx[str(g)], str(g)) for j, g in enumerate(adapter_genes) if str(g) in arc_g_idx]

    if not common_p_pairs or len(common_g_pairs) < 3:
        return {
            "line_name": ARC_ADAPTER_CHECK_LINE,
            "n_common_perts": len(common_p_pairs),
            "n_common_genes": len(common_g_pairs),
            "n_finite_perts": 0,
            "median_per_pert_pearson": float("nan"),
            "pooled_pearson": float("nan"),
            "min_per_pert_pearson": float("nan"),
            "q25_per_pert_pearson": float("nan"),
            "q75_per_pert_pearson": float("nan"),
            "median_pearson_min": float(median_pearson_min),
            "pooled_pearson_min": float(pooled_pearson_min),
            "passed": False,
            "stop_for_controller_diagnosis": True,
        }

    ad_rows = np.asarray([x[0] for x in common_p_pairs], dtype=np.int64)
    ar_rows = np.asarray([x[1] for x in common_p_pairs], dtype=np.int64)
    ad_cols = np.asarray([x[0] for x in common_g_pairs], dtype=np.int64)
    ar_cols = np.asarray([x[1] for x in common_g_pairs], dtype=np.int64)

    A = ad_mat[np.ix_(ad_rows, ad_cols)]
    B = ar_mat[np.ix_(ar_rows, ar_cols)]

    per_pert_rs: List[float] = []
    for k in range(A.shape[0]):
        x = A[k]
        y = B[k]
        if np.all(np.isfinite(x)) and np.all(np.isfinite(y)):
            sx = float(np.std(x))
            sy = float(np.std(y))
            if sx > 1e-12 and sy > 1e-12:
                per_pert_rs.append(float(np.corrcoef(x, y)[0, 1]))

    r_arr = np.asarray(per_pert_rs, dtype=np.float64)
    med_r = float(np.median(r_arr)) if r_arr.size > 0 else float("nan")
    a_flat = A.ravel()
    b_flat = B.ravel()
    fin = np.isfinite(a_flat) & np.isfinite(b_flat)
    if int(np.sum(fin)) >= 3 and float(np.std(a_flat[fin])) > 1e-12 and float(np.std(b_flat[fin])) > 1e-12:
        pool_r = float(np.corrcoef(a_flat[fin], b_flat[fin])[0, 1])
    else:
        pool_r = float("nan")

    passed = bool(
        np.isfinite(med_r)
        and np.isfinite(pool_r)
        and med_r >= float(median_pearson_min)
        and pool_r >= float(pooled_pearson_min)
    )
    return {
        "line_name": ARC_ADAPTER_CHECK_LINE,
        "n_common_perts": len(common_p_pairs),
        "n_common_genes": len(common_g_pairs),
        "n_finite_perts": int(r_arr.size),
        "median_per_pert_pearson": med_r,
        "pooled_pearson": pool_r,
        "min_per_pert_pearson": float(np.min(r_arr)) if r_arr.size > 0 else float("nan"),
        "q25_per_pert_pearson": float(np.percentile(r_arr, 25)) if r_arr.size > 0 else float("nan"),
        "q75_per_pert_pearson": float(np.percentile(r_arr, 75)) if r_arr.size > 0 else float("nan"),
        "median_pearson_min": float(median_pearson_min),
        "pooled_pearson_min": float(pooled_pearson_min),
        "passed": passed,
        "stop_for_controller_diagnosis": not passed,
    }


def verify_arc_archived_adata_pred_r16(
    adata_pred_path: Path,
    adapter_deltas: np.ndarray,
    adapter_perts: Sequence[str],
    adapter_genes: Sequence[str],
    line_name: str = ARC_ADAPTER_CHECK_LINE,
    status_dir: Path = DEFAULT_STATUS_DIR,
    hgnc_bridge: Optional[HGNCBridge] = None,
    median_pearson_min: float = ARC_ADAPTER_MEDIAN_PEARSON_MIN,
    pooled_pearson_min: float = ARC_ADAPTER_POOLED_PEARSON_MIN,
) -> Dict[str, Any]:
    """
    Triple-gated Arc archived `adata_pred.h5ad` adapter fidelity check (`Addendum A Item 3`).
    Strictly restricted to `nadig_hepg2` ONLY and strictly post-`R16_GO`.
    """
    assert_all_execution_gates_r16(status_dir)
    if str(line_name) != ARC_ADAPTER_CHECK_LINE:
        raise ValueError(
            f"Addendum A Item 3 locks Arc adata_pred.h5ad verification strictly to "
            f"{ARC_ADAPTER_CHECK_LINE!r} ONLY; got {line_name!r}."
        )

    p = Path(adata_pred_path)
    if not p.is_file():
        raise FileNotFoundError(f"Arc archived adata_pred.h5ad not found: {p}")

    with h5py.File(p, "r") as f:
        obs = f["obs"]
        obs_keys = set(obs.keys())
        pert_col = next(
            (c for c in ("gene", "target_gene", "perturbation", "pert", "pert_name", "target", "condition") if c in obs_keys),
            None,
        )
        if pert_col is None:
            raise KeyError(f"Archived adata_pred {p} missing perturbation column: {sorted(obs_keys)}")
        pert_vals = _read_h5ad_obs_col_r16(obs, pert_col)

        var = f["var"]
        var_keys = set(var.keys())
        vname_col = next(
            (c for c in ("gene_name", "gene_symbol", "feature_name", "gene_names", "symbol") if c in var_keys),
            None,
        )
        if vname_col is None:
            vname_col = _resolve_h5ad_obs_index_col_r16(var)
        if vname_col is None:
            raise KeyError(f"Archived adata_pred {p} missing var gene-name column: {sorted(var_keys)}")
        var_genes_raw = _read_h5ad_obs_col_r16(var, vname_col)

        def _canon(s: str) -> str:
            return (hgnc_bridge.canonicalize(s) or s) if hgnc_bridge is not None else s

        arc_genes = [_canon(g) for g in var_genes_raw]
        adapter_genes_c = [_canon(g) for g in adapter_genes]
        adapter_perts_c = [_canon(pt) for pt in adapter_perts]
        adapter_pert_set = set(adapter_perts_c)

        ntc_rows = np.array(
            [i for i, v in enumerate(pert_vals) if str(v).strip().lower() in ("non-targeting", "nontargeting")],
            dtype=np.int64,
        )
        if ntc_rows.size == 0:
            raise RuntimeError(f"Archived adata_pred {p} contains zero 'non-targeting' rows.")

        pert_to_rows: Dict[str, List[int]] = defaultdict(list)
        for i, v in enumerate(pert_vals):
            sv = str(v).strip()
            if sv.lower() in ("non-targeting", "nontargeting"):
                continue
            cv = _canon(sv)
            if cv in adapter_pert_set:
                pert_to_rows[cv].append(i)

        arc_perts = sorted(pert_to_rows.keys())
        n_cells_total = len(pert_vals)
        n_genes = len(arc_genes)
        n_buckets = 1 + len(arc_perts)
        row_bucket = np.full(n_cells_total, -1, dtype=np.int32)
        row_bucket[ntc_rows] = 0
        for k, pk in enumerate(arc_perts):
            row_bucket[np.asarray(pert_to_rows[pk], dtype=np.int64)] = k + 1

        sums = np.zeros((n_buckets, n_genes), dtype=np.float64)
        counts = np.zeros(n_buckets, dtype=np.int64)
        x_node = f["X"]
        chunk_rows = 4096

        if isinstance(x_node, h5py.Dataset):
            # Single sequential pass over dense /X (bounded ~64 MB RAM on 3.15 GB adata_pred.h5ad)
            for s in range(0, n_cells_total, chunk_rows):
                e = min(s + chunk_rows, n_cells_total)
                b_chunk = row_bucket[s:e]
                keep_local = np.flatnonzero(b_chunk >= 0)
                if keep_local.size == 0:
                    continue
                block = np.asarray(x_node[s:e, :], dtype=np.float64)
                b_keep = b_chunk[keep_local]
                for bi in np.unique(b_keep):
                    sel = keep_local[b_keep == bi]
                    sums[int(bi)] += block[sel].sum(axis=0)
                    counts[int(bi)] += int(sel.size)
        else:
            enc = str(x_node.attrs.get("encoding-type", "")).lower()
            if "csc" in enc:
                shape = tuple(int(x) for x in x_node.attrs["shape"])
                x_csr = sp.csc_matrix(
                    (np.asarray(x_node["data"][:], dtype=np.float64), x_node["indices"][:], x_node["indptr"][:]),
                    shape=shape,
                ).tocsr()
                sums[0] = np.asarray(x_csr[ntc_rows].sum(axis=0), dtype=np.float64).ravel()
                counts[0] = int(ntc_rows.size)
                for k, pk in enumerate(arc_perts):
                    rk = np.asarray(pert_to_rows[pk], dtype=np.int64)
                    sums[k + 1] = np.asarray(x_csr[rk].sum(axis=0), dtype=np.float64).ravel()
                    counts[k + 1] = int(rk.size)
            else:
                indptr = x_node["indptr"][:]
                data_ds = x_node["data"]
                idx_ds = x_node["indices"]
                for s in range(0, n_cells_total, chunk_rows):
                    e = min(s + chunk_rows, n_cells_total)
                    b_chunk = row_bucket[s:e]
                    keep_local = np.flatnonzero(b_chunk >= 0)
                    if keep_local.size == 0:
                        continue
                    b_keep = b_chunk[keep_local]
                    p0, p1 = int(indptr[s]), int(indptr[e])
                    if p1 == p0:
                        for bi in np.unique(b_keep):
                            counts[int(bi)] += int(np.sum(b_keep == bi))
                        continue
                    csr_blk = sp.csr_matrix(
                        (np.asarray(data_ds[p0:p1], dtype=np.float64), idx_ds[p0:p1], indptr[s : e + 1] - p0),
                        shape=(e - s, n_genes),
                    )
                    for bi in np.unique(b_keep):
                        sel = keep_local[b_keep == bi]
                        sums[int(bi)] += np.asarray(csr_blk[sel].sum(axis=0), dtype=np.float64).ravel()
                        counts[int(bi)] += int(sel.size)

        ntc_mean = sums[0] / float(max(int(counts[0]), 1))
        arc_deltas = np.empty((len(arc_perts), n_genes), dtype=np.float64)
        for k in range(len(arc_perts)):
            arc_deltas[k] = (sums[k + 1] / float(max(int(counts[k + 1]), 1))) - ntc_mean

    res = compare_predicted_deltas_vs_arc_archive_r16(
        adapter_deltas=adapter_deltas,
        adapter_perts=adapter_perts_c,
        adapter_genes=adapter_genes_c,
        arc_deltas=arc_deltas,
        arc_perts=arc_perts,
        arc_genes=arc_genes,
        median_pearson_min=median_pearson_min,
        pooled_pearson_min=pooled_pearson_min,
    )
    res["adata_pred_path"] = str(p)
    res["n_arc_ntc_cells"] = int(ntc_rows.size)
    return res


# ==============================================================================
# Model Vocabulary Loading, Practical Metrics, & Auditable Array Engine
# ==============================================================================


def load_state_vocab_r16(
    models_dir: Path = DEFAULT_MODELS_DIR,
    hgnc_bridge: Optional[HGNCBridge] = None,
) -> Tuple[List[str], Dict[str, str]]:
    """
    Load exact Gate-2 verified `G_model` (`2,000` genes) and `P_model` (`2,023` non-control keys
    mapped from canonical HGNC -> exact `pert_onehot_map.pt` key).
    """
    import torch

    fold_dir = Path(models_dir) / "ST-SE-Replogle" / "zeroshot" / "hepg2"
    with open(fold_dir / "var_dims.pkl", "rb") as fh:
        g_model = list(pickle.load(fh)["gene_names"])
    pmap = torch.load(fold_dir / "pert_onehot_map.pt", map_location="cpu", weights_only=False)
    pert_keys = sorted(k for k in pmap.keys() if k != "non-targeting")
    if hgnc_bridge is None:
        return g_model, {k: k for k in pert_keys}
    canon_to_model: Dict[str, str] = {}
    for k in pert_keys:
        ck = hgnc_bridge.canonicalize(k) or k
        canon_to_model.setdefault(ck, k)
    return g_model, canon_to_model


def compute_practical_accuracy_metrics_r16(
    delta_obs: np.ndarray,
    delta_pred: np.ndarray,
    target_indices: Sequence[int],
    top_k: int = 50,
) -> Dict[str, float]:
    """Descriptive `H-D_util` practical accuracy metrics (`raw_pearson`, `centered_de50_pearson`, `cosine`)."""
    obs = np.asarray(delta_obs, dtype=np.float64)
    pred = np.asarray(delta_pred, dtype=np.float64)
    P, G = obs.shape
    raw_rs, de50_rs, cos_vals = [], [], []
    for i in range(P):
        mask = np.ones(G, dtype=bool)
        ti = int(target_indices[i])
        if 0 <= ti < G:
            mask[ti] = False
        x = obs[i, mask]
        y = pred[i, mask]
        sx, sy = float(np.std(x)), float(np.std(y))
        if sx > 1e-12 and sy > 1e-12:
            raw_rs.append(float(np.corrcoef(x, y)[0, 1]))
        nx, ny = float(np.linalg.norm(x)), float(np.linalg.norm(y))
        if nx > 1e-12 and ny > 1e-12:
            cos_vals.append(float(np.dot(x, y) / (nx * ny)))
        k = min(int(top_k), x.size)
        if k >= 3:
            top_idx = np.argpartition(np.abs(x), -k)[-k:]
            xt, yt = x[top_idx], y[top_idx]
            sxt, syt = float(np.std(xt)), float(np.std(yt))
            if sxt > 1e-12 and syt > 1e-12:
                de50_rs.append(float(np.corrcoef(xt, yt)[0, 1]))
    return {
        "mean_raw_pearson": float(np.mean(raw_rs)) if raw_rs else float("nan"),
        "mean_centered_de50_pearson": float(np.mean(de50_rs)) if de50_rs else float("nan"),
        "mean_cosine_similarity": float(np.mean(cos_vals)) if cos_vals else float("nan"),
        "n_perts": int(P),
    }


def compute_g4_factorial_decomposition_r16(
    delta_obs: np.ndarray,
    pred_m_bm: np.ndarray,
    pred_m_bs: np.ndarray,
    pred_s_bm: np.ndarray,
    pred_s_bs: np.ndarray,
    delta_transfer: np.ndarray,
    target_indices: Sequence[int],
    source_grand_mean: Optional[np.ndarray],
    delta_half_a: Optional[np.ndarray],
    delta_half_b: Optional[np.ndarray],
    half_valid: Optional[np.ndarray],
    n_boot: int = 4000,
    seed: int = DEFAULT_R16_SEED,
) -> Dict[str, Any]:
    """
    `PREREG_I v5 §4.4` Ledger `G4` 2x2 factorial decomposition (`X_ctrl` x `batch`):
      - `delta_batch` ("只换 batch 不换对照"): `r(M, B_M) - r(M, B_S)`
      - `delta_ctrl`  ("只换对照不换 batch"): `r(M, B_M) - r(S, B_M)`
      - `delta_both`: `r(M, B_M) - r(S, B_S)`
      - `delta_int`: `delta_both - delta_ctrl - delta_batch`
      - `batch_share`: `|delta_batch| / (|delta_ctrl| + |delta_batch|)`
    """
    out: Dict[str, Any] = {}
    for metric_key, b_s in (("r_orth_cal", None), ("r_dorth_cal", source_grand_mean)):
        c_batch = compute_paired_contrast_r16(
            delta_obs, pred_m_bm, pred_m_bs, delta_transfer, target_indices,
            source_grand_mean=b_s, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
            half_valid=half_valid, n_boot=n_boot, seed=seed,
        )
        c_ctrl = compute_paired_contrast_r16(
            delta_obs, pred_m_bm, pred_s_bm, delta_transfer, target_indices,
            source_grand_mean=b_s, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
            half_valid=half_valid, n_boot=n_boot, seed=seed,
        )
        c_both = compute_paired_contrast_r16(
            delta_obs, pred_m_bm, pred_s_bs, delta_transfer, target_indices,
            source_grand_mean=b_s, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
            half_valid=half_valid, n_boot=n_boot, seed=seed,
        )
        db = float(c_batch["estimate"])
        dc = float(c_ctrl["estimate"])
        dt = float(c_both["estimate"])
        denom = abs(db) + abs(dc)
        out[metric_key] = {
            "r_M_BM": c_batch["a_summary"],
            "r_M_BS": c_batch["b_summary"],
            "r_S_BM": c_ctrl["b_summary"],
            "r_S_BS": c_both["b_summary"],
            "delta_batch": c_batch,
            "delta_ctrl": c_ctrl,
            "delta_both": c_both,
            "delta_interaction_point": float(dt - dc - db),
            "batch_share_point": float(abs(db) / denom) if denom > 1e-12 else float("nan"),
        }
    return out


def save_line_arrays_r16(
    arrays_dir: Path,
    line_name: str,
    observed_payload: Dict[str, Any],
    pred_arms: Dict[str, Dict[str, Any]],
) -> Dict[str, str]:
    """
    Save auditable `observed.npz` and per-arm `preds_<arm>.npz` under `arrays_dir / line_name`,
    returning `{relative_filename: sha256}` manifest for independent three-party verification.
    """
    ldir = Path(arrays_dir) / line_name
    ldir.mkdir(parents=True, exist_ok=True)
    manifest: Dict[str, str] = {}

    obs_path = ldir / "observed.npz"
    np.savez_compressed(obs_path, **observed_payload)
    manifest[obs_path.name] = _sha256_file(obs_path)

    for arm_name, payload in sorted(pred_arms.items()):
        fname = f"preds_{_safe_arm_filename(arm_name)}.npz"
        fpath = ldir / fname
        np.savez_compressed(fpath, **payload)
        manifest[fname] = _sha256_file(fpath)

    manifest_path = ldir / "array_manifest_sha256.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def save_control_embeddings_r16(
    embeddings_dir: Path,
    control_embeddings: Dict[str, Dict[str, Any]],
) -> Dict[str, str]:
    """
    `PREREG_I Addendum A Item 11`: Save each control pool's `2,058`-dim `float32` `X_state`
    embedding array as `x_state_<safe_key>.npz` (`X_state`, `control_barcodes`, `seed`,
    `dataset_name`, `se_audit_json`) plus `embeddings_manifest_sha256.json` (`{filename: sha256}`)
    for Master Controller cross-party comparison (`max_abs <= 1e-5`) prior to predicted deltas.
    """
    edir = Path(embeddings_dir)
    edir.mkdir(parents=True, exist_ok=True)
    manifest: Dict[str, str] = {}
    seen_fnames: Set[str] = set()

    for pool_key, entry in sorted(control_embeddings.items()):
        safe_key = _safe_arm_filename(pool_key)
        if not safe_key or not re.search(r"[A-Za-z0-9]", safe_key):
            raise ValueError(
                f"Control embedding pool_key must contain alphanumeric characters, got {pool_key!r}."
            )
        x_arr = np.ascontiguousarray(entry["X_state"], dtype=np.float32)
        if x_arr.ndim != 2 or x_arr.shape[0] <= 0 or x_arr.shape[1] != 2058:
            raise ValueError(
                f"Control embedding {pool_key!r} must have shape (N > 0, 2058), got {x_arr.shape}."
            )
        if not np.all(np.isfinite(x_arr)):
            raise ValueError(f"Control embedding {pool_key!r} contains non-finite values.")
        raw_barcodes = entry.get("control_barcodes")
        barcodes = [str(b) for b in raw_barcodes] if raw_barcodes is not None else []
        if barcodes and len(barcodes) != x_arr.shape[0]:
            raise ValueError(
                f"Control embedding {pool_key!r}: len(control_barcodes)={len(barcodes)} "
                f"!= X_state rows={x_arr.shape[0]}."
            )
        seed_raw = entry.get("seed")
        if seed_raw is None or isinstance(seed_raw, bool):
            raise ValueError(f"Control embedding {pool_key!r} requires an integer seed, got {seed_raw!r}.")
        seed_val = int(seed_raw)
        ds_name = str(entry.get("dataset_name", pool_key))
        raw_audit = entry.get("se_audit_json")
        if isinstance(raw_audit, str):
            se_audit_json = raw_audit
        else:
            se_audit_json = json.dumps(entry.get("se_audit") or {}, sort_keys=True)

        fname = f"x_state_{safe_key}.npz"
        if fname in seen_fnames:
            raise ValueError(f"Duplicate sanitized embedding filename {fname!r} for key {pool_key!r}.")
        seen_fnames.add(fname)

        fpath = edir / fname
        np.savez_compressed(
            fpath,
            X_state=x_arr,
            control_barcodes=np.asarray(barcodes, dtype=str),
            seed=np.int64(seed_val),
            dataset_name=np.asarray(ds_name),
            se_audit_json=np.asarray(se_audit_json),
        )
        manifest[fname] = _sha256_file(fpath)

    manifest_path = edir / "embeddings_manifest_sha256.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def evaluate_line_from_arrays_r16(
    line_name: str,
    group_key: str,
    modality: str,
    checkpoint_id: str,
    seen_sources: Sequence[str],
    unseen_sources: Sequence[str],
    delta_obs: np.ndarray,
    delta_half_a: np.ndarray,
    delta_half_b: np.ndarray,
    half_valid: np.ndarray,
    delta_transfer: np.ndarray,
    source_grand_mean: np.ndarray,
    target_indices: np.ndarray,
    target_names: Sequence[str],
    gene_names: Sequence[str],
    mu_target: np.ndarray,
    mu_sources: np.ndarray,
    fct_delta_pred: np.ndarray,
    e1_constmed_delta_pred: np.ndarray,
    state_se_matched: np.ndarray,
    state_se_mismatched: Dict[str, np.ndarray],
    state_se_resample: Optional[np.ndarray] = None,
    state_hvg_matched: Optional[np.ndarray] = None,
    state_hvg_mismatched: Optional[Dict[str, np.ndarray]] = None,
    state_hvg_resample: Optional[np.ndarray] = None,
    state_hvg_libnorm1e4_matched: Optional[np.ndarray] = None,
    state_se_neutral_batch_arms: Optional[Dict[str, Tuple[np.ndarray, Dict[str, np.ndarray]]]] = None,
    state_se_per_fold_arms: Optional[Dict[str, Tuple[np.ndarray, Dict[str, np.ndarray]]]] = None,
    g4_factorial_cells: Optional[Dict[str, Dict[str, np.ndarray]]] = None,
    g4_single_index_arms: Optional[Dict[str, np.ndarray]] = None,
    embedding_diagnostics: Optional[Dict[str, Any]] = None,
    ctrl_audit: Optional[Dict[str, Any]] = None,
    control_barcodes_by_source: Optional[Dict[str, List[str]]] = None,
    resample_barcodes: Optional[List[str]] = None,
    extra_fct_sensitivity: Optional[Dict[str, Any]] = None,
    arrays_dir: Optional[Path] = None,
    n_boot: int = 4000,
    n_perm_hvg: int = 100,
    seed: int = DEFAULT_R16_SEED,
    delta_tost: float = DEFAULT_TOST_DELTA,
) -> Dict[str, Any]:
    """
    Pure, deterministic Round 16 line evaluator & auditable `.npz` exporter
    (`PREREG_I v5 + Addendum A Items 1, 4, 6, 7, 8, 9, 10`).
    Can be called identically on Gate 5 real arrays or Gate 4 `--synthetic-smoke` arrays.
    """
    all_mis_expected = set(seen_sources) | set(unseen_sources)
    if set(state_se_mismatched.keys()) != all_mis_expected:
        raise ValueError(
            f"[{line_name}] state_se_mismatched keys {sorted(state_se_mismatched.keys())} "
            f"!= expected {sorted(all_mis_expected)}"
        )

    # 1. PRIMARY `H-I1` & `H-I2` (`STATE_SE`, `r_orth_cal`, `all_nontarget`)
    se_orth = compute_orth_suite_r16(
        delta_obs, state_se_matched, delta_transfer, target_indices,
        source_grand_mean=None, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
        half_valid=half_valid, n_boot=n_boot, seed=seed,
    )
    se_paired_mis = compute_paired_model_vs_all_mismatched_r16(
        delta_obs, state_se_matched, state_se_mismatched, delta_transfer, target_indices,
        source_grand_mean=None, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
        half_valid=half_valid, n_boot=n_boot, seed=seed,
    )
    h_i1_adj = adjudicate_claim_i1_line(
        se_orth, se_paired_mis,
        seen_sources=seen_sources, unseen_sources=unseen_sources, delta=delta_tost,
    )

    fct_orth = compute_orth_suite_r16(
        delta_obs, fct_delta_pred, delta_transfer, target_indices,
        source_grand_mean=None, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
        half_valid=half_valid, n_boot=n_boot, seed=seed,
    )
    se_vs_fct = compute_paired_contrast_r16(
        delta_obs, state_se_matched, fct_delta_pred, delta_transfer, target_indices,
        source_grand_mean=None, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
        half_valid=half_valid, n_boot=n_boot, seed=seed,
    )
    h_i2_adj = adjudicate_claim_i2_line(se_vs_fct)

    # 2. SECONDARY `§4.1` `r_dorth_cal` mirror
    se_dorth = compute_orth_suite_r16(
        delta_obs, state_se_matched, delta_transfer, target_indices,
        source_grand_mean=source_grand_mean, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
        half_valid=half_valid, n_boot=n_boot, seed=seed,
    )
    se_dorth_mis = compute_paired_model_vs_all_mismatched_r16(
        delta_obs, state_se_matched, state_se_mismatched, delta_transfer, target_indices,
        source_grand_mean=source_grand_mean, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
        half_valid=half_valid, n_boot=n_boot, seed=seed,
    )
    fct_dorth = compute_orth_suite_r16(
        delta_obs, fct_delta_pred, delta_transfer, target_indices,
        source_grand_mean=source_grand_mean, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
        half_valid=half_valid, n_boot=n_boot, seed=seed,
    )
    se_vs_fct_dorth = compute_paired_contrast_r16(
        delta_obs, state_se_matched, fct_delta_pred, delta_transfer, target_indices,
        source_grand_mean=source_grand_mean, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
        half_valid=half_valid, n_boot=n_boot, seed=seed,
    )

    # 2b. SECONDARY `Addendum A Items 8 & 9` Within-line disjoint resample reference (`STATE_SE`)
    within_line_resample_summary: Optional[Dict[str, Any]] = None
    if state_se_resample is not None:
        resamp_orth = compute_paired_with_within_line_resample_r16(
            delta_obs, state_se_matched, state_se_resample, state_se_mismatched,
            delta_transfer, target_indices, source_grand_mean=None,
            delta_half_a=delta_half_a, delta_half_b=delta_half_b,
            half_valid=half_valid, n_boot=n_boot, seed=seed,
        )
        resamp_dorth = compute_paired_with_within_line_resample_r16(
            delta_obs, state_se_matched, state_se_resample, state_se_mismatched,
            delta_transfer, target_indices, source_grand_mean=source_grand_mean,
            delta_half_a=delta_half_a, delta_half_b=delta_half_b,
            half_valid=half_valid, n_boot=n_boot, seed=seed,
        )
        within_line_resample_summary = {
            "r_orth_cal": resamp_orth,
            "r_dorth_cal": resamp_dorth,
            "descriptive_interpretation": summarize_within_line_resample_r16(
                resamp_orth,
                seen_sources=seen_sources,
                unseen_sources=unseen_sources,
                matched_positive_ci975=h_i1_adj["matched_positive_ci975"],
            ),
        }

    # 3. SECONDARY `§4.2` Mechanism Hypothesis `M2`
    m2_weight_shift = compute_winner_curse_weight_shift_r16(
        mu_target, mu_sources, delta_half_a, delta_half_b, delta_transfer, target_indices,
        half_valid=half_valid, n_boot=n_boot, seed=seed,
    )
    m2_cross_mask: Dict[str, Any] = {}
    cm_models: Dict[str, np.ndarray] = {
        "E1_ConstMed": e1_constmed_delta_pred,
        "FCT": fct_delta_pred,
        "STATE_SE": state_se_matched,
    }
    if state_hvg_matched is not None:
        cm_models["STATE_HVG"] = state_hvg_matched
    for mk, mp_arr in cm_models.items():
        m2_cross_mask[mk] = {
            "r_orth_cal": compute_joint_cross_mask_r16(
                delta_obs, mp_arr, delta_transfer, target_indices,
                source_grand_mean=None, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
                half_valid=half_valid, n_boot=n_boot, seed=seed,
            ),
            "r_dorth_cal": compute_joint_cross_mask_r16(
                delta_obs, mp_arr, delta_transfer, target_indices,
                source_grand_mean=source_grand_mean, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
                half_valid=half_valid, n_boot=n_boot, seed=seed,
            ),
        }

    # 4. SECONDARY `§4.3` & `Addendum A Item 4` `STATE_HVG` (`official_median_log1p`, `Perm100`, `libnorm1e4`)
    hvg_summary: Optional[Dict[str, Any]] = None
    if state_hvg_matched is not None and state_hvg_mismatched is not None:
        hvg_orth = compute_orth_suite_r16(
            delta_obs, state_hvg_matched, delta_transfer, target_indices,
            source_grand_mean=None, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
            half_valid=half_valid, n_boot=n_boot, seed=seed,
        )
        hvg_paired_mis = compute_paired_model_vs_all_mismatched_r16(
            delta_obs, state_hvg_matched, state_hvg_mismatched, delta_transfer, target_indices,
            source_grand_mean=None, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
            half_valid=half_valid, n_boot=n_boot, seed=seed,
        )
        hvg_desc_hi1 = adjudicate_claim_i1_line(
            hvg_orth, hvg_paired_mis, seen_sources=seen_sources, unseen_sources=unseen_sources, delta=delta_tost
        )
        hvg_vs_fct = compute_paired_contrast_r16(
            delta_obs, state_hvg_matched, fct_delta_pred, delta_transfer, target_indices,
            source_grand_mean=None, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
            half_valid=half_valid, n_boot=n_boot, seed=seed,
        )
        hvg_dorth = compute_orth_suite_r16(
            delta_obs, state_hvg_matched, delta_transfer, target_indices,
            source_grand_mean=source_grand_mean, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
            half_valid=half_valid, n_boot=n_boot, seed=seed,
        )
        hvg_perm_orth = compute_gene_label_perm_null_r16(
            delta_obs, state_hvg_matched, delta_transfer, target_indices,
            source_grand_mean=None, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
            half_valid=half_valid, n_perm=n_perm_hvg, perm_seed_offset=0,
        )
        hvg_perm_dorth = compute_gene_label_perm_null_r16(
            delta_obs, state_hvg_matched, delta_transfer, target_indices,
            source_grand_mean=source_grand_mean, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
            half_valid=half_valid, n_perm=n_perm_hvg, perm_seed_offset=0,
        )
        hvg_summary = {
            "hvg_norm_primary": "official_median_log1p",
            "g_model_mapping": "HGNCBridge",
            "r_orth_cal_matched": hvg_orth,
            "r_orth_cal_paired_vs_mismatched": hvg_paired_mis,
            "r_orth_cal_descriptive_h_i1": hvg_desc_hi1,
            "r_orth_cal_vs_fct": hvg_vs_fct,
            "r_dorth_cal_matched": hvg_dorth,
            "perm100_r_orth_cal": hvg_perm_orth,
            "perm100_r_dorth_cal": hvg_perm_dorth,
            "practical_accuracy": compute_practical_accuracy_metrics_r16(delta_obs, state_hvg_matched, target_indices),
        }
        if state_hvg_resample is not None:
            hvg_resamp_orth = compute_paired_with_within_line_resample_r16(
                delta_obs, state_hvg_matched, state_hvg_resample, state_hvg_mismatched,
                delta_transfer, target_indices, source_grand_mean=None,
                delta_half_a=delta_half_a, delta_half_b=delta_half_b,
                half_valid=half_valid, n_boot=n_boot, seed=seed,
            )
            hvg_summary["within_line_resample_reference"] = {
                "r_orth_cal": hvg_resamp_orth,
                "descriptive_interpretation": summarize_within_line_resample_r16(
                    hvg_resamp_orth,
                    seen_sources=seen_sources,
                    unseen_sources=unseen_sources,
                    matched_positive_ci975=hvg_desc_hi1["matched_positive_ci975"],
                ),
            }
        if state_hvg_libnorm1e4_matched is not None:
            hvg_summary["sensitivity_libnorm1e4_log1p_matched"] = compute_orth_suite_r16(
                delta_obs, state_hvg_libnorm1e4_matched, delta_transfer, target_indices,
                source_grand_mean=None, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
                half_valid=half_valid, n_boot=n_boot, seed=seed,
            )

    # 5. SECONDARY `§4.4` Neutral-batch robustness panel, `G4` 2x2 factorial, & single-index spread
    neutral_panel_summary: Dict[str, Any] = {}
    if state_se_neutral_batch_arms:
        for bmode, (m_arr, mis_dict) in sorted(state_se_neutral_batch_arms.items()):
            o_m = compute_orth_suite_r16(
                delta_obs, m_arr, delta_transfer, target_indices,
                source_grand_mean=None, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
                half_valid=half_valid, n_boot=n_boot, seed=seed,
            )
            p_mis = compute_paired_model_vs_all_mismatched_r16(
                delta_obs, m_arr, mis_dict, delta_transfer, target_indices,
                source_grand_mean=None, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
                half_valid=half_valid, n_boot=n_boot, seed=seed,
            )
            neutral_panel_summary[bmode] = {
                "matched_orth": o_m,
                "paired_vs_mismatched": p_mis,
                "descriptive_h_i1": adjudicate_claim_i1_line(
                    o_m, p_mis, seen_sources=seen_sources, unseen_sources=unseen_sources, delta=delta_tost
                ),
                "vs_fct": compute_paired_contrast_r16(
                    delta_obs, m_arr, fct_delta_pred, delta_transfer, target_indices,
                    source_grand_mean=None, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
                    half_valid=half_valid, n_boot=n_boot, seed=seed,
                ),
            }

    per_fold_summary: Dict[str, Any] = {}
    if state_se_per_fold_arms:
        for fkey, (m_arr, mis_dict) in sorted(state_se_per_fold_arms.items()):
            o_f = compute_orth_suite_r16(
                delta_obs, m_arr, delta_transfer, target_indices,
                source_grand_mean=None, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
                half_valid=half_valid, n_boot=n_boot, seed=seed,
            )
            p_f = compute_paired_model_vs_all_mismatched_r16(
                delta_obs, m_arr, mis_dict, delta_transfer, target_indices,
                source_grand_mean=None, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
                half_valid=half_valid, n_boot=n_boot, seed=seed,
            )
            per_fold_summary[fkey] = {"matched_orth": o_f, "paired_vs_mismatched": p_f}

    g4_summary: Dict[str, Any] = {}
    if g4_factorial_cells:
        for pair_label, cells in sorted(g4_factorial_cells.items()):
            g4_summary[pair_label] = compute_g4_factorial_decomposition_r16(
                delta_obs=delta_obs,
                pred_m_bm=cells["M_BM"],
                pred_m_bs=cells["M_BS"],
                pred_s_bm=cells["S_BM"],
                pred_s_bs=cells["S_BS"],
                delta_transfer=delta_transfer,
                target_indices=target_indices,
                source_grand_mean=source_grand_mean,
                delta_half_a=delta_half_a,
                delta_half_b=delta_half_b,
                half_valid=half_valid,
                n_boot=n_boot,
                seed=seed,
            )

    # Addendum A Item 6: Stratum-level G4 summaries (`seen_by_checkpoint` vs `unseen_by_checkpoint`)
    g4_stratum_summary: Dict[str, Any] = {}
    if g4_summary:
        for prefix, p_name in (
            ("b0_vs_b17_mis_", "b0_vs_b17"),
            ("b17_vs_b0_sym_mis_", "b17_vs_b0_sym"),
            ("native_vs_", "native"),
        ):
            pairing_out: Dict[str, Any] = {}
            for st_name, st_list in (
                ("seen_by_checkpoint", sorted(seen_sources)),
                ("unseen_by_checkpoint", sorted(unseen_sources)),
            ):
                matched_sks = [sk for sk in st_list if f"{prefix}{sk}" in g4_summary]
                if not matched_sks:
                    continue
                st_metrics: Dict[str, Any] = {"n_sources": len(matched_sks), "sources": matched_sks}
                for mkey in ("r_orth_cal", "r_dorth_cal"):
                    m_out: Dict[str, Any] = {}
                    for fld in (
                        "delta_batch",
                        "delta_ctrl",
                        "delta_both",
                        "delta_interaction_point",
                        "batch_share_point",
                    ):
                        vals: List[float] = []
                        for sk in matched_sks:
                            entry = g4_summary[f"{prefix}{sk}"][mkey][fld]
                            v = float(entry["estimate"]) if isinstance(entry, dict) else float(entry)
                            if np.isfinite(v):
                                vals.append(v)
                        m_out[fld] = {
                            "n_finite": len(vals),
                            "median": float(np.median(vals)) if vals else float("nan"),
                            "mean": float(np.mean(vals)) if vals else float("nan"),
                            "min": float(min(vals)) if vals else float("nan"),
                            "max": float(max(vals)) if vals else float("nan"),
                        }
                    st_metrics[mkey] = m_out
                pairing_out[st_name] = st_metrics
            if pairing_out:
                g4_stratum_summary[p_name] = pairing_out

    g4_spread_summary: Dict[str, Any] = {}
    if g4_single_index_arms:
        per_idx_orth: Dict[str, Any] = {}
        per_idx_dorth: Dict[str, Any] = {}
        orth_pts: List[float] = []
        dorth_pts: List[float] = []
        for b_key, b_arr in sorted(g4_single_index_arms.items()):
            ro = compute_orth_suite_r16(
                delta_obs, b_arr, delta_transfer, target_indices,
                source_grand_mean=None, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
                half_valid=half_valid, n_boot=n_boot, seed=seed,
            )
            rd = compute_orth_suite_r16(
                delta_obs, b_arr, delta_transfer, target_indices,
                source_grand_mean=source_grand_mean, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
                half_valid=half_valid, n_boot=n_boot, seed=seed,
            )
            per_idx_orth[b_key] = ro
            per_idx_dorth[b_key] = rd
            if np.isfinite(ro["estimate"]):
                orth_pts.append(float(ro["estimate"]))
            if np.isfinite(rd["estimate"]):
                dorth_pts.append(float(rd["estimate"]))
        g4_spread_summary = {
            "per_index_r_orth_cal": per_idx_orth,
            "per_index_r_dorth_cal": per_idx_dorth,
            "r_orth_cal_min": float(min(orth_pts)) if orth_pts else float("nan"),
            "r_orth_cal_max": float(max(orth_pts)) if orth_pts else float("nan"),
            "r_orth_cal_range": float(max(orth_pts) - min(orth_pts)) if orth_pts else float("nan"),
            "r_orth_cal_std": float(np.std(orth_pts, ddof=1)) if len(orth_pts) > 1 else 0.0,
            "r_dorth_cal_min": float(min(dorth_pts)) if dorth_pts else float("nan"),
            "r_dorth_cal_max": float(max(dorth_pts)) if dorth_pts else float("nan"),
            "r_dorth_cal_range": float(max(dorth_pts) - min(dorth_pts)) if dorth_pts else float("nan"),
        }

    # Addendum A Item 1: Group 3 `full r9_6` sensitivity fixed on primary `N_eval^R16` using its own `FCT` & `Transfer_Mean`
    fct_sens_summary: Dict[str, Any] = {}
    fct_sens_exports: Dict[str, Dict[str, np.ndarray]] = {}
    if extra_fct_sensitivity:
        for skey, s_payload in sorted(extra_fct_sensitivity.items()):
            if isinstance(s_payload, dict):
                s_pred = np.asarray(s_payload["fct_delta_pred"], dtype=np.float64)
                s_tm = np.asarray(s_payload.get("delta_transfer", delta_transfer), dtype=np.float64)
                s_bs = np.asarray(s_payload.get("source_grand_mean", source_grand_mean), dtype=np.float64)
            else:
                s_pred = np.asarray(s_payload, dtype=np.float64)
                s_tm = np.asarray(delta_transfer, dtype=np.float64)
                s_bs = np.asarray(source_grand_mean, dtype=np.float64)
            fct_sens_exports[skey] = {
                "fct_delta_pred": s_pred,
                "delta_transfer": s_tm,
                "source_grand_mean": s_bs,
            }
            fct_sens_summary[skey] = {
                "fixed_to_primary_n_eval": int(len(target_names)),
                "fct_orth": compute_orth_suite_r16(
                    delta_obs, s_pred, s_tm, target_indices,
                    source_grand_mean=None, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
                    half_valid=half_valid, n_boot=n_boot, seed=seed,
                ),
                "state_se_orth": compute_orth_suite_r16(
                    delta_obs, state_se_matched, s_tm, target_indices,
                    source_grand_mean=None, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
                    half_valid=half_valid, n_boot=n_boot, seed=seed,
                ),
                "state_se_vs_sens_fct": compute_paired_contrast_r16(
                    delta_obs, state_se_matched, s_pred, s_tm, target_indices,
                    source_grand_mean=None, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
                    half_valid=half_valid, n_boot=n_boot, seed=seed,
                ),
                "fct_dorth": compute_orth_suite_r16(
                    delta_obs, s_pred, s_tm, target_indices,
                    source_grand_mean=s_bs, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
                    half_valid=half_valid, n_boot=n_boot, seed=seed,
                ),
                "state_se_dorth": compute_orth_suite_r16(
                    delta_obs, state_se_matched, s_tm, target_indices,
                    source_grand_mean=s_bs, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
                    half_valid=half_valid, n_boot=n_boot, seed=seed,
                ),
                "state_se_vs_sens_fct_dorth": compute_paired_contrast_r16(
                    delta_obs, state_se_matched, s_pred, s_tm, target_indices,
                    source_grand_mean=s_bs, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
                    half_valid=half_valid, n_boot=n_boot, seed=seed,
                ),
            }

    # 6. Auditable `.npz` array & barcode exports (`observed.npz` + `preds_*.npz`)
    bc_map = control_barcodes_by_source or {}
    matched_bcs = list(bc_map.get(line_name, (ctrl_audit or {}).get("sampled_barcodes", [])))
    resamp_bcs = list(
        resample_barcodes
        if resample_barcodes is not None
        else (ctrl_audit or {}).get("within_line_resample", {}).get("sampled_barcodes", [])
    )

    array_manifest: Dict[str, str] = {}
    if arrays_dir is not None:
        obs_payload = {
            "delta_obs": np.asarray(delta_obs, dtype=np.float64),
            "delta_half_a": np.asarray(delta_half_a, dtype=np.float64),
            "delta_half_b": np.asarray(delta_half_b, dtype=np.float64),
            "half_valid": np.asarray(half_valid, dtype=bool),
            "delta_transfer": np.asarray(delta_transfer, dtype=np.float64),
            "source_grand_mean": np.asarray(source_grand_mean, dtype=np.float64),
            "fct_delta_pred": np.asarray(fct_delta_pred, dtype=np.float64),
            "e1_constmed_delta_pred": np.asarray(e1_constmed_delta_pred, dtype=np.float64),
            "target_indices": np.asarray(target_indices, dtype=np.int64),
            "target_names": np.asarray(list(target_names), dtype=str),
            "gene_names": np.asarray(list(gene_names), dtype=str),
            "mu_target": np.asarray(mu_target, dtype=np.float64),
            "mu_sources": np.asarray(mu_sources, dtype=np.float64),
            "n_ctrl_total": np.int64((ctrl_audit or {}).get("n_ctrl_total", -1)),
            "n_ctrl_sampled": np.int64((ctrl_audit or {}).get("n_ctrl_sampled", -1)),
            "ctrl_sampling_fraction": np.float64((ctrl_audit or {}).get("ctrl_sampling_fraction", float("nan"))),
            "matched_ctrl_barcodes": np.asarray(matched_bcs, dtype=str),
            "resample_ctrl_barcodes": np.asarray(resamp_bcs, dtype=str),
        }
        pred_arms: Dict[str, Dict[str, Any]] = {}

        def _pack(
            arr: np.ndarray,
            mkey: str,
            bmode: str,
            csrc: str,
            hnorm: str = "",
            barcodes: Optional[Sequence[str]] = None,
        ) -> Dict[str, Any]:
            bcs_use = list(barcodes) if barcodes is not None else list(bc_map.get(csrc, []))
            return {
                "delta_pred": np.asarray(arr, dtype=np.float64),
                "target_names": np.asarray(list(target_names), dtype=str),
                "gene_names": np.asarray(list(gene_names), dtype=str),
                "model_key": np.asarray(mkey),
                "checkpoint_id": np.asarray(checkpoint_id),
                "batch_mode": np.asarray(bmode),
                "control_source": np.asarray(csrc),
                "hvg_norm_mode": np.asarray(hnorm),
                "control_barcodes": np.asarray(bcs_use, dtype=str),
            }

        pred_arms["STATE_SE_fixed_index_zero_matched"] = _pack(
            state_se_matched, "STATE_SE", "fixed_index_zero", line_name, barcodes=matched_bcs
        )
        if state_se_resample is not None:
            pred_arms["STATE_SE_fixed_index_zero_resample"] = _pack(
                state_se_resample,
                "STATE_SE",
                "fixed_index_zero",
                f"{line_name}__within_line_resample_seed{DEFAULT_RESAMPLE_SEED}",
                barcodes=resamp_bcs,
            )
        for sk, arr in sorted(state_se_mismatched.items()):
            pred_arms[f"STATE_SE_fixed_index_zero_mis_{sk}"] = _pack(
                arr, "STATE_SE", "fixed_index_zero", sk
            )
        if state_hvg_matched is not None:
            pred_arms["STATE_HVG_official_median_log1p_matched"] = _pack(
                state_hvg_matched, "STATE_HVG", "fixed_index_zero", line_name, "official_median_log1p", barcodes=matched_bcs
            )
        if state_hvg_resample is not None:
            pred_arms["STATE_HVG_official_median_log1p_resample"] = _pack(
                state_hvg_resample,
                "STATE_HVG",
                "fixed_index_zero",
                f"{line_name}__within_line_resample_seed{DEFAULT_RESAMPLE_SEED}",
                "official_median_log1p",
                barcodes=resamp_bcs,
            )
        if state_hvg_mismatched is not None:
            for sk, arr in sorted(state_hvg_mismatched.items()):
                pred_arms[f"STATE_HVG_official_median_log1p_mis_{sk}"] = _pack(
                    arr, "STATE_HVG", "fixed_index_zero", sk, "official_median_log1p"
                )
        if state_hvg_libnorm1e4_matched is not None:
            pred_arms["STATE_HVG_libnorm1e4_log1p_matched"] = _pack(
                state_hvg_libnorm1e4_matched, "STATE_HVG", "fixed_index_zero", line_name, "libnorm1e4_log1p", barcodes=matched_bcs
            )
        if state_se_neutral_batch_arms:
            for bmode, (m_arr, mis_dict) in sorted(state_se_neutral_batch_arms.items()):
                pred_arms[f"STATE_SE_{bmode}_matched"] = _pack(
                    m_arr, "STATE_SE", bmode, line_name, barcodes=matched_bcs
                )
                for sk, arr in sorted(mis_dict.items()):
                    pred_arms[f"STATE_SE_{bmode}_mis_{sk}"] = _pack(arr, "STATE_SE", bmode, sk)
        if state_se_per_fold_arms:
            for fkey, (m_arr, mis_dict) in sorted(state_se_per_fold_arms.items()):
                pred_arms[f"STATE_SE_fold_{fkey}_matched"] = _pack(
                    m_arr, "STATE_SE", "fixed_index_zero", line_name, barcodes=matched_bcs
                )
                for sk, arr in sorted(mis_dict.items()):
                    pred_arms[f"STATE_SE_fold_{fkey}_mis_{sk}"] = _pack(arr, "STATE_SE", "fixed_index_zero", sk)
        if g4_factorial_cells:
            for plbl, cells in sorted(g4_factorial_cells.items()):
                if plbl.startswith("b0_vs_b17_mis_"):
                    s_src = plbl[len("b0_vs_b17_mis_") :]
                elif plbl.startswith("b17_vs_b0_sym_mis_"):
                    s_src = plbl[len("b17_vs_b0_sym_mis_") :]
                elif plbl.startswith("native_vs_"):
                    s_src = plbl[len("native_vs_") :]
                else:
                    s_src = plbl
                for ckey, carr in sorted(cells.items()):
                    if ckey.startswith("M_"):
                        g4_csrc = line_name
                        g4_bcs = matched_bcs
                    else:
                        g4_csrc = s_src
                        g4_bcs = list(bc_map.get(s_src, []))
                    pred_arms[f"STATE_SE_G4_{plbl}_{ckey}"] = _pack(
                        carr, "STATE_SE", f"G4:{plbl}:{ckey}", g4_csrc, barcodes=g4_bcs
                    )
        if g4_single_index_arms:
            for b_key, b_arr in sorted(g4_single_index_arms.items()):
                if b_key == "b0":
                    continue
                pred_arms[f"STATE_SE_fixed_index_{b_key}_matched"] = _pack(
                    b_arr, "STATE_SE", f"fixed_index:{b_key}", line_name, barcodes=matched_bcs
                )
        for skey, s_exp in sorted(fct_sens_exports.items()):
            pred_arms[f"SENS_{skey}"] = {
                "fct_delta_pred": s_exp["fct_delta_pred"],
                "delta_transfer": s_exp["delta_transfer"],
                "source_grand_mean": s_exp["source_grand_mean"],
                "target_names": np.asarray(list(target_names), dtype=str),
                "gene_names": np.asarray(list(gene_names), dtype=str),
            }

        array_manifest = save_line_arrays_r16(arrays_dir, line_name, obs_payload, pred_arms)

    return {
        "line_name": line_name,
        "group_key": group_key,
        "modality": modality,
        "checkpoint_id": checkpoint_id,
        "n_eval_targets": int(len(target_names)),
        "n_eval_genes": int(len(gene_names)),
        "n_half_valid": int(np.sum(half_valid)),
        "n_self_masked": int(np.sum(target_indices >= 0)),
        "strata_sources": {
            "seen_by_checkpoint": list(seen_sources),
            "unseen_by_checkpoint": list(unseen_sources),
        },
        "primary": {
            "STATE_SE_matched_orth": se_orth,
            "STATE_SE_paired_vs_mismatched": se_paired_mis,
            "H_I1_adjudication": h_i1_adj,
            "FCT_matched_orth": fct_orth,
            "STATE_SE_vs_FCT": se_vs_fct,
            "H_I2_adjudication": h_i2_adj,
        },
        "secondary": {
            "r_dorth_cal": {
                "STATE_SE_matched": se_dorth,
                "STATE_SE_paired_vs_mismatched": se_dorth_mis,
                "FCT_matched": fct_dorth,
                "STATE_SE_vs_FCT": se_vs_fct_dorth,
            },
            "within_line_resample_reference": within_line_resample_summary,
            "mechanism_m2": {
                "saturation_weight_shift": m2_weight_shift,
                "cross_mask_shifts": m2_cross_mask,
            },
            "STATE_HVG": hvg_summary,
            "neutral_batch_robustness_panel": neutral_panel_summary,
            "per_checkpoint_breakdown": per_fold_summary,
            "g4_factorial_decomposition": g4_summary,
            "g4_stratum_summary": g4_stratum_summary,
            "g4_single_index_spread": g4_spread_summary,
            "practical_accuracy_STATE_SE": compute_practical_accuracy_metrics_r16(
                delta_obs, state_se_matched, target_indices
            ),
            "practical_accuracy_FCT": compute_practical_accuracy_metrics_r16(
                delta_obs, fct_delta_pred, target_indices
            ),
            "extra_fct_sensitivity": fct_sens_summary,
        },
        "ctrl_audit": ctrl_audit or {},
        "embedding_diagnostics": embedding_diagnostics or {},
        "array_manifest_sha256": array_manifest,
    }


# ==============================================================================
# Gate 5 Real-Data Orchestrator (Strictly Triple-Gated)
# ==============================================================================


def _align_line_and_sources_r16(
    line_name: str,
    fct_source_names: Sequence[str],
    line_dir: Path,
    g_model_2000: Sequence[str],
    p_model_canon_to_key: Dict[str, str],
    hgnc_bridge: HGNCBridge,
    common_genes_r9: Sequence[str],
    fixed_target_names: Optional[Sequence[str]] = None,
    fixed_gene_names: Optional[Sequence[str]] = None,
) -> Dict[str, Any]:
    """
    Align a target line and its locked `FCT` sources strictly onto `(T_eval^R16, G_eval^R16)`.
    When `fixed_target_names` and/or `fixed_gene_names` are supplied (Addendum A Item 1),
    strictly locks evaluation onto those exact primary targets/genes while letting the
    sensitivity source pool (`full r9_6`) update ONLY `FCT`, `Transfer_Mean`, and `Source_Grand_Mean`.
    """
    ds_t = load_standardized_line_r10(
        line_dir, representation="log1p", common_genes=common_genes_r9, hgnc_bridge=hgnc_bridge
    )
    src_ds_map: Dict[str, LineDatasetR10] = {
        sn: load_standardized_line_r10(
            R9_LINES_DIR / sn, representation="log1p", common_genes=common_genes_r9, hgnc_bridge=hgnc_bridge
        )
        for sn in fct_source_names
    }

    shared_r9_genes_set = set(ds_t.genes).intersection(*(set(src_ds_map[sn].genes) for sn in fct_source_names))
    if fixed_gene_names is not None:
        missing_g = [g for g in fixed_gene_names if g not in shared_r9_genes_set]
        if missing_g:
            raise ValueError(
                f"[{line_name}] fixed_gene_names has {len(missing_g)} genes absent from source pool "
                f"{list(fct_source_names)} (e.g. {missing_g[:5]})."
            )
        g_eval_r16 = list(fixed_gene_names)
        g_r9_set = set(shared_r9_genes_set)
    else:
        g_r9_aligned = [g for g in ds_t.genes if g in shared_r9_genes_set]
        g_r9_set = set(g_r9_aligned)
        g_model_canon_set = {hgnc_bridge.canonicalize(g) or g for g in g_model_2000}
        g_eval_r16 = [g for g in g_r9_aligned if g in g_model_canon_set]

    g_eval_map = {g: i for i, g in enumerate(g_eval_r16)}
    t_pos = {g: i for i, g in enumerate(ds_t.genes)}
    t_cols = np.asarray([t_pos[g] for g in g_eval_r16], dtype=np.int64)

    if fixed_target_names is not None:
        t_row_map = {t: i for i, t in enumerate(ds_t.target_names)}
        missing_t = [t for t in fixed_target_names if t not in t_row_map]
        if missing_t:
            raise ValueError(
                f"[{line_name}] fixed_target_names has {len(missing_t)} targets absent from line "
                f"(e.g. {missing_t[:5]})."
            )
        row_arr = np.asarray([t_row_map[t] for t in fixed_target_names], dtype=np.int64)
        q_names = list(fixed_target_names)
    else:
        src_repr_targets: Set[str] = set()
        for sn in fct_source_names:
            s_ds = src_ds_map[sn]
            if s_ds.repr_mask is not None:
                src_repr_targets.update(t for t, ok in zip(s_ds.target_names, s_ds.repr_mask) if ok)
            else:
                src_repr_targets.update(s_ds.target_names)
        keep_rows = [
            i
            for i, t in enumerate(ds_t.target_names)
            if t in g_r9_set
            and (ds_t.repr_mask is None or bool(ds_t.repr_mask[i]))
            and t in src_repr_targets
            and t in p_model_canon_to_key
        ]
        row_arr = np.asarray(keep_rows, dtype=np.int64)
        q_names = [ds_t.target_names[int(i)] for i in row_arr]
    q_gidx = np.asarray([g_eval_map.get(t, -1) for t in q_names], dtype=np.int64)

    # Fit `FctSourceState` strictly on `G_eval^R16`
    src_inputs = []
    for sn in fct_source_names:
        s_ds = src_ds_map[sn]
        s_pos = {g: i for i, g in enumerate(s_ds.genes)}
        s_cols = np.asarray([s_pos[g] for g in g_eval_r16], dtype=np.int64)
        s_tidx = np.asarray([g_eval_map.get(t, -1) for t in s_ds.target_names], dtype=np.int64)
        s_sub = LineDatasetR10(
            name=s_ds.name,
            genes=list(g_eval_r16),
            target_names=list(s_ds.target_names),
            target_indices=s_tidx,
            deltas=s_ds.deltas[:, s_cols],
            ctrl_cells=s_ds.ctrl_cells[:, s_cols],
            ctrl_mean_log1p=s_ds.ctrl_mean_log1p[s_cols] if s_ds.ctrl_mean_log1p is not None else None,
            ctrl_mean_libnorm=s_ds.ctrl_mean_libnorm[s_cols] if s_ds.ctrl_mean_libnorm is not None else None,
            ctrl_mean_libnorm_full=(
                s_ds.ctrl_mean_libnorm_full[s_cols] if s_ds.ctrl_mean_libnorm_full is not None else None
            ),
            repr_mask=s_ds.repr_mask,
            half_a_deltas=s_ds.half_a_deltas[:, s_cols] if s_ds.half_a_deltas is not None else None,
            half_b_deltas=s_ds.half_b_deltas[:, s_cols] if s_ds.half_b_deltas is not None else None,
            half_valid=s_ds.half_valid,
            n_cells=s_ds.n_cells,
            knockdown_depth_log1p=s_ds.knockdown_depth_log1p,
            knockdown_fraction=s_ds.knockdown_fraction,
            target_ensembl_ids=s_ds.target_ensembl_ids,
            representation=s_ds.representation,
            metadata=s_ds.metadata,
        )
        src_inputs.append(build_source_input_r10(s_sub, only_reproducible=True, compute_cov=False))

    fct_state = FctSourceState.fit(src_inputs, mask_self=True)
    delta_transfer, fc_bar, has_hits = compute_source_transfer_and_fc_bar(fct_state, q_names)
    mu_target = np.maximum(np.asarray(ds_t.ctrl_mean_libnorm[t_cols], dtype=np.float64), 0.0)
    mu_sources = np.vstack([np.maximum(s.ctrl_mean_libnorm, 0.0) for s in src_inputs])
    fct_pred = fct_from_fc_bar(mu_target, fc_bar, delta_transfer, has_hits)
    c_med = float(np.median(mu_target))
    e1_pred = fct_from_fc_bar(np.full_like(mu_target, c_med), fc_bar, delta_transfer, has_hits)

    return {
        "gene_names": g_eval_r16,
        "target_names": q_names,
        "target_indices": q_gidx,
        "delta_obs": ds_t.deltas[row_arr][:, t_cols],
        "delta_half_a": ds_t.half_a_deltas[row_arr][:, t_cols],
        "delta_half_b": ds_t.half_b_deltas[row_arr][:, t_cols],
        "half_valid": ds_t.half_valid[row_arr],
        "delta_transfer": delta_transfer,
        "source_grand_mean": fct_state.grand_mean,
        "mu_target": mu_target,
        "mu_sources": mu_sources,
        "fc_bar": fc_bar,
        "has_hits": has_hits,
        "fct_delta_pred": fct_pred,
        "e1_constmed_delta_pred": e1_pred,
    }


def run_real_evaluation_r16(
    lines_to_run: Sequence[str],
    status_dir: Path = DEFAULT_STATUS_DIR,
    models_dir: Path = DEFAULT_MODELS_DIR,
    results_dir: Path = DEFAULT_RESULTS_DIR,
    device: str = "mps",
    n_boot: int = 4000,
    n_perm_hvg: int = 100,
    seed: int = DEFAULT_R16_SEED,
    resample_seed: int = DEFAULT_RESAMPLE_SEED,
    include_secondary_panels: bool = True,
    gse295214_raw_dir: Optional[Path] = None,
    raw_nadig_dir: Optional[Path] = None,
    verify_arc_adata_pred_path: Optional[Path] = None,
) -> Dict[str, Any]:
    """
    Gate 5 real-cell evaluation orchestrator (`PREREG_I v5 + Addendum A`).
    STRICTLY refuses execution unless `PREREG_I_COMMIT.txt`, `CODE_FREEZE_R16.txt`, and `R16_GO.txt`
    all exist and are non-empty.
    """
    gate_paths = assert_all_execution_gates_r16(status_dir)
    if verify_arc_adata_pred_path is not None and ARC_ADAPTER_CHECK_LINE not in set(lines_to_run):
        raise ValueError(
            f"Addendum A Item 3 restricts --verify-arc-adata-pred strictly to {ARC_ADAPTER_CHECK_LINE!r}, "
            f"which is not in lines_to_run={list(lines_to_run)}."
        )
    os.environ.setdefault("WANDB_MODE", "disabled")
    os.environ.setdefault("HF_HUB_OFFLINE", "1")
    os.environ.setdefault("TRANSFORMERS_OFFLINE", "1")

    t0 = time.time()
    bridge = HGNCBridge.from_tsv(DEFAULT_HGNC_FILE)
    common_genes_r9 = load_r9_common_genes(hgnc_bridge=bridge, common_genes_file=R9_COMMON_GENES_FILE)
    g_model_2000, p_model_canon_to_key = load_state_vocab_r16(models_dir, bridge)
    gene_canon_map = {g: (bridge.canonicalize(g) or g) for g in g_model_2000}

    # Determine all control pools needed across requested lines + their two mismatched strata
    target_line_set = set(lines_to_run)
    ctrl_pool_lines: Set[str] = set()
    strata_by_line: Dict[str, Dict[str, Tuple[str, ...]]] = {}
    for ln in lines_to_run:
        spec = LINE_SPECS_R16[ln]
        assert_zero_training_leakage(
            "STATE_SE", ln, checkpoint_id=spec["checkpoint_id"], source_pool_lines=spec["fct_sources"]
        )
        strata = build_mismatched_strata_r16(ln, spec["checkpoint_id"])
        strata_by_line[ln] = strata
        ctrl_pool_lines.add(ln)
        ctrl_pool_lines.update(strata["seen_by_checkpoint"])
        ctrl_pool_lines.update(strata["unseen_by_checkpoint"])

    # Step 1: Extract raw NTC counts + exact barcodes & cache SE-600M (`X_state`) + ST-HVG (`X_hvg`)
    ctrl_cache: Dict[str, Dict[str, Any]] = {}
    control_source_barcodes: Dict[str, Dict[str, Any]] = {}
    control_embeddings_to_save: Dict[str, Dict[str, Any]] = {}
    for cln in sorted(ctrl_pool_lines):
        want_resamp = cln in target_line_set
        bundle = extract_local_ntc_bundle_r16(
            cln,
            status_dir=status_dir,
            seed=seed,
            resample_seed=resample_seed,
            include_resample=want_resamp,
            gse295214_raw_dir=gse295214_raw_dir,
            raw_nadig_dir=raw_nadig_dir,
        )
        raw_mat = bundle["raw_counts"]
        raw_genes = bundle["gene_names"]
        batch_labels = bundle["batch_labels"]
        barcodes = bundle["barcodes"]
        c_audit = bundle["ctrl_audit"]

        x_state, se_audit = encode_controls_se600m_r16(
            raw_counts=raw_mat,
            raw_gene_names=raw_genes,
            se_dir=Path(models_dir) / "SE-600M",
            device=device,
            batch_size=16,
            dataset_name=cln,
            seed=seed,
        )
        control_embeddings_to_save[cln] = {
            "X_state": x_state,
            "control_barcodes": list(barcodes),
            "seed": int(seed),
            "dataset_name": cln,
            "se_audit": se_audit,
        }
        raw_canon = {g: (bridge.canonicalize(g) or g) for g in raw_genes}
        full_canon = {**raw_canon, **gene_canon_map}
        x_hvg_official, hvg_off_audit = preprocess_hvg_counts_r16(
            raw_mat, raw_genes, g_model_2000, norm_mode="official_median_log1p", canon_map=full_canon
        )
        x_hvg_1e4, hvg_1e4_audit = preprocess_hvg_counts_r16(
            raw_mat, raw_genes, g_model_2000, norm_mode="libnorm1e4_log1p", canon_map=full_canon
        )

        # Addendum A Items 8 & 11: encode disjoint within-line resample control pool with `resample_seed` (20261602)
        x_state_resamp = None
        x_hvg_off_resamp = None
        se_resamp_audit: Optional[Dict[str, Any]] = None
        if bundle["resample_raw_counts"] is not None:
            r_mat = bundle["resample_raw_counts"]
            x_state_resamp, se_resamp_audit = encode_controls_se600m_r16(
                raw_counts=r_mat,
                raw_gene_names=raw_genes,
                se_dir=Path(models_dir) / "SE-600M",
                device=device,
                batch_size=16,
                dataset_name=f"{cln}_resample",
                seed=resample_seed,
            )
            control_embeddings_to_save[f"{cln}_resample"] = {
                "X_state": x_state_resamp,
                "control_barcodes": list(bundle["resample_barcodes"] or []),
                "seed": int(resample_seed),
                "dataset_name": f"{cln}_resample",
                "se_audit": se_resamp_audit or {},
            }
            x_hvg_off_resamp, _ = preprocess_hvg_counts_r16(
                r_mat, raw_genes, g_model_2000, norm_mode="official_median_log1p", canon_map=full_canon
            )
            del r_mat

        ctrl_cache[cln] = {
            "x_state": x_state,
            "x_state_resample": x_state_resamp,
            "x_hvg_official": x_hvg_official,
            "x_hvg_official_resample": x_hvg_off_resamp,
            "x_hvg_1e4": x_hvg_1e4,
            "batch_labels": batch_labels,
            "barcodes": barcodes,
            "resample_barcodes": bundle["resample_barcodes"] or [],
            "resample_batch_labels": bundle["resample_batch_labels"] or [],
            "ctrl_audit": c_audit,
            "se_audit": se_audit,
            "se_resample_audit": se_resamp_audit,
            "hvg_official_audit": hvg_off_audit,
            "hvg_1e4_audit": hvg_1e4_audit,
        }
        control_source_barcodes[cln] = {
            "primary_seed": int(seed),
            "n_ctrl_total": int(c_audit["n_ctrl_total"]),
            "n_ctrl_sampled": int(c_audit["n_ctrl_sampled"]),
            "barcodes": list(barcodes),
            "batch_labels": list(batch_labels),
            "within_line_resample": {
                "seed": int(resample_seed),
                "eligible": bool(bundle["resample_audit"].get("eligible", False)),
                "regime": bundle["resample_audit"].get("regime", ""),
                "n_resample": int(bundle["resample_audit"].get("n_resample", 0)),
                "overlap_with_primary": int(bundle["resample_audit"].get("overlap_with_primary", 0)),
                "barcodes": list(bundle["resample_barcodes"] or []),
                "batch_labels": list(bundle["resample_batch_labels"] or []),
            },
        }
        del raw_mat, bundle
        gc.collect()

    Path(results_dir).mkdir(parents=True, exist_ok=True)
    barcodes_json_path = Path(results_dir) / "control_barcodes_r16.json"
    barcodes_json_path.write_text(
        json.dumps(control_source_barcodes, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    barcodes_by_source = {cln: ctrl_cache[cln]["barcodes"] for cln in sorted(ctrl_cache.keys())}

    embeddings_dir = Path(results_dir) / "embeddings"
    embeddings_manifest = save_control_embeddings_r16(embeddings_dir, control_embeddings_to_save)

    arrays_dir = Path(results_dir) / "arrays"
    line_results: Dict[str, Any] = {}
    arc_adapter_check_res: Optional[Dict[str, Any]] = None

    for ln in lines_to_run:
        spec = LINE_SPECS_R16[ln]
        strata = strata_by_line[ln]
        seen_s = strata["seen_by_checkpoint"]
        unseen_s = strata["unseen_by_checkpoint"]
        all_mis = list(seen_s) + list(unseen_s)

        aligned = _align_line_and_sources_r16(
            line_name=ln,
            fct_source_names=spec["fct_sources"],
            line_dir=spec["line_dir"],
            g_model_2000=g_model_2000,
            p_model_canon_to_key=p_model_canon_to_key,
            hgnc_bridge=bridge,
            common_genes_r9=common_genes_r9,
        )
        q_names = aligned["target_names"]
        g_eval = aligned["gene_names"]

        # Primary `STATE_SE` (`fixed_index_zero`)
        se_adapter = ArcStateAdapter(
            variant="STATE_SE",
            checkpoint_id=spec["checkpoint_id"],
            models_dir=models_dir,
            device=device,
            batch_mode="fixed_index_zero",
            status_dir=status_dir,
        )
        return_folds = spec["checkpoint_id"] == "zeroshot/ensemble4"
        res_m = se_adapter.predict_deltas(
            target_line=ln, control_line_label=ln,
            ctrl_cells=ctrl_cache[ln]["x_state"],
            pert_names=q_names, gene_names=g_eval,
            ctrl_is_embedded_x_state=True,
            pert_canon_to_model=p_model_canon_to_key,
            gene_canon_map=gene_canon_map,
            return_per_fold=return_folds,
        )
        se_matched = res_m[0] if return_folds else res_m
        per_fold_m: Dict[str, np.ndarray] = res_m[1] if return_folds else {}
        per_fold_mis: Dict[str, Dict[str, np.ndarray]] = {fk: {} for fk in per_fold_m}

        # Addendum A Item 8: Within-line disjoint resample arm (`seed = 20261602`)
        se_resample: Optional[np.ndarray] = None
        if ctrl_cache[ln]["x_state_resample"] is not None:
            se_resample = se_adapter.predict_deltas(
                target_line=ln, control_line_label=f"{ln}_resample",
                ctrl_cells=ctrl_cache[ln]["x_state_resample"],
                pert_names=q_names, gene_names=g_eval,
                ctrl_is_embedded_x_state=True,
                pert_canon_to_model=p_model_canon_to_key,
                gene_canon_map=gene_canon_map,
                return_per_fold=False,
            )

        se_mis: Dict[str, np.ndarray] = {}
        emb_diag: Dict[str, Any] = {
            "matched_se_audit": ctrl_cache[ln]["se_audit"],
            "resample_se_audit": ctrl_cache[ln]["se_resample_audit"],
            "vs_mismatched_centroid_separation": {},
        }
        if ctrl_cache[ln]["x_state_resample"] is not None:
            emb_diag["vs_within_line_resample_centroid_separation"] = compute_embedding_separation_diagnostics(
                ctrl_cache[ln]["x_state"], ctrl_cache[ln]["x_state_resample"]
            )
        for sk in all_mis:
            r_s = se_adapter.predict_deltas(
                target_line=ln, control_line_label=sk,
                ctrl_cells=ctrl_cache[sk]["x_state"],
                pert_names=q_names, gene_names=g_eval,
                ctrl_is_embedded_x_state=True,
                pert_canon_to_model=p_model_canon_to_key,
                gene_canon_map=gene_canon_map,
                return_per_fold=return_folds,
            )
            se_mis[sk] = r_s[0] if return_folds else r_s
            if return_folds:
                for fk, farr in r_s[1].items():
                    per_fold_mis[fk][sk] = farr
            emb_diag["vs_mismatched_centroid_separation"][sk] = compute_embedding_separation_diagnostics(
                ctrl_cache[ln]["x_state"], ctrl_cache[sk]["x_state"]
            )

        per_fold_arms = (
            {fk: (per_fold_m[fk], per_fold_mis[fk]) for fk in sorted(per_fold_m.keys())}
            if return_folds
            else None
        )

        # Secondary `STATE_HVG`, neutral-batch panel, & `G4` 2x2 factorial
        hvg_matched = hvg_mis = hvg_resample = hvg_1e4_m = None
        neutral_arms = None
        g4_cells: Optional[Dict[str, Dict[str, np.ndarray]]] = None
        g4_single_idx_arms: Optional[Dict[str, np.ndarray]] = None
        if include_secondary_panels:
            hvg_adapter = ArcStateAdapter(
                variant="STATE_HVG", checkpoint_id=spec["checkpoint_id"],
                models_dir=models_dir, device=device, batch_mode="fixed_index_zero",
                hvg_norm_mode="official_median_log1p", status_dir=status_dir,
            )
            hvg_matched = hvg_adapter.predict_deltas(
                target_line=ln, control_line_label=ln,
                ctrl_cells=ctrl_cache[ln]["x_hvg_official"],
                pert_names=q_names, gene_names=g_eval,
                pert_canon_to_model=p_model_canon_to_key, gene_canon_map=gene_canon_map,
            )
            if ctrl_cache[ln]["x_hvg_official_resample"] is not None:
                hvg_resample = hvg_adapter.predict_deltas(
                    target_line=ln, control_line_label=f"{ln}_resample",
                    ctrl_cells=ctrl_cache[ln]["x_hvg_official_resample"],
                    pert_names=q_names, gene_names=g_eval,
                    pert_canon_to_model=p_model_canon_to_key, gene_canon_map=gene_canon_map,
                )
            hvg_mis = {
                sk: hvg_adapter.predict_deltas(
                    target_line=ln, control_line_label=sk,
                    ctrl_cells=ctrl_cache[sk]["x_hvg_official"],
                    pert_names=q_names, gene_names=g_eval,
                    pert_canon_to_model=p_model_canon_to_key, gene_canon_map=gene_canon_map,
                )
                for sk in all_mis
            }
            hvg_1e4_m = hvg_adapter.predict_deltas(
                target_line=ln, control_line_label=ln,
                ctrl_cells=ctrl_cache[ln]["x_hvg_1e4"],
                pert_names=q_names, gene_names=g_eval,
                pert_canon_to_model=p_model_canon_to_key, gene_canon_map=gene_canon_map,
            )

            # Neutral-batch robustness panel (`§4.4`)
            neutral_arms = {}
            for bmode in ("bypass_zero", "mean_batch_emb", "multi_batch_ens8"):
                nb_ad = ArcStateAdapter(
                    variant="STATE_SE", checkpoint_id=spec["checkpoint_id"],
                    models_dir=models_dir, device=device, batch_mode=bmode, status_dir=status_dir,
                )
                m_nb = nb_ad.predict_deltas(
                    target_line=ln, control_line_label=ln, ctrl_cells=ctrl_cache[ln]["x_state"],
                    pert_names=q_names, gene_names=g_eval, ctrl_is_embedded_x_state=True,
                    pert_canon_to_model=p_model_canon_to_key, gene_canon_map=gene_canon_map,
                )
                mis_nb = {
                    sk: nb_ad.predict_deltas(
                        target_line=ln, control_line_label=sk, ctrl_cells=ctrl_cache[sk]["x_state"],
                        pert_names=q_names, gene_names=g_eval, ctrl_is_embedded_x_state=True,
                        pert_canon_to_model=p_model_canon_to_key, gene_canon_map=gene_canon_map,
                    )
                    for sk in all_mis
                }
                neutral_arms[bmode] = (m_nb, mis_nb)

            # `G4` 2x2 factorial (`PREREG_I §4.4` + Addendum A Item 6)
            g4_cells = {}
            if spec["group"] == "group1":
                counterpart = seen_s[0]
                if len(ctrl_cache[ln]["batch_labels"]) != len(ctrl_cache[counterpart]["batch_labels"]):
                    raise ValueError(
                        f"Group 1 G4 1-to-1 positional alignment requires equal sampled cell counts: "
                        f"{ln}={len(ctrl_cache[ln]['batch_labels'])} vs {counterpart}={len(ctrl_cache[counterpart]['batch_labels'])}."
                    )
                nat_ad = ArcStateAdapter(
                    variant="STATE_SE", checkpoint_id=spec["checkpoint_id"],
                    models_dir=models_dir, device=device, batch_mode="native_gem_group", status_dir=status_dir,
                )
                m_bm = nat_ad.predict_deltas(
                    ln, ln, ctrl_cache[ln]["x_state"], q_names, g_eval,
                    ctrl_is_embedded_x_state=True, pert_canon_to_model=p_model_canon_to_key,
                    gene_canon_map=gene_canon_map, native_batch_labels=ctrl_cache[ln]["batch_labels"],
                )
                g4_cells[f"native_vs_{counterpart}"] = {
                    "M_BM": m_bm,
                    "M_BS": nat_ad.predict_deltas(
                        ln, ln, ctrl_cache[ln]["x_state"], q_names, g_eval,
                        ctrl_is_embedded_x_state=True, pert_canon_to_model=p_model_canon_to_key,
                        gene_canon_map=gene_canon_map, native_batch_labels=ctrl_cache[counterpart]["batch_labels"],
                    ),
                    "S_BM": nat_ad.predict_deltas(
                        ln, counterpart, ctrl_cache[counterpart]["x_state"], q_names, g_eval,
                        ctrl_is_embedded_x_state=True, pert_canon_to_model=p_model_canon_to_key,
                        gene_canon_map=gene_canon_map, native_batch_labels=ctrl_cache[ln]["batch_labels"],
                    ),
                    "S_BS": nat_ad.predict_deltas(
                        ln, counterpart, ctrl_cache[counterpart]["x_state"], q_names, g_eval,
                        ctrl_is_embedded_x_state=True, pert_canon_to_model=p_model_canon_to_key,
                        gene_canon_map=gene_canon_map, native_batch_labels=ctrl_cache[counterpart]["batch_labels"],
                    ),
                }
            else:
                b17_ad = ArcStateAdapter(
                    variant="STATE_SE", checkpoint_id=spec["checkpoint_id"],
                    models_dir=models_dir, device=device, batch_mode="fixed_index:17", status_dir=status_dir,
                )
                m_b17 = b17_ad.predict_deltas(
                    ln, ln, ctrl_cache[ln]["x_state"], q_names, g_eval,
                    ctrl_is_embedded_x_state=True, pert_canon_to_model=p_model_canon_to_key,
                    gene_canon_map=gene_canon_map,
                )
                # Full G4 coverage across ALL mismatched sources (`(0, 17)` + `(17, 0)` symmetry check)
                for sk in all_mis:
                    s_b17 = b17_ad.predict_deltas(
                        ln, sk, ctrl_cache[sk]["x_state"], q_names, g_eval,
                        ctrl_is_embedded_x_state=True, pert_canon_to_model=p_model_canon_to_key,
                        gene_canon_map=gene_canon_map,
                    )
                    g4_cells[f"b0_vs_b17_mis_{sk}"] = {
                        "M_BM": se_matched,
                        "M_BS": m_b17,
                        "S_BM": se_mis[sk],
                        "S_BS": s_b17,
                    }
                    g4_cells[f"b17_vs_b0_sym_mis_{sk}"] = {
                        "M_BM": m_b17,
                        "M_BS": se_matched,
                        "S_BM": s_b17,
                        "S_BS": se_mis[sk],
                    }
                # Single-index spread across `b in {0, 1, 17, 35}` (`PREREG_I §4.4.3`)
                g4_single_idx_arms = {"b0": se_matched, "b17": m_b17}
                for b_idx in G4_SPREAD_BATCH_INDICES:
                    if b_idx in (0, 17):
                        continue
                    bi_ad = ArcStateAdapter(
                        variant="STATE_SE", checkpoint_id=spec["checkpoint_id"],
                        models_dir=models_dir, device=device, batch_mode=f"fixed_index:{b_idx}", status_dir=status_dir,
                    )
                    g4_single_idx_arms[f"b{b_idx}"] = bi_ad.predict_deltas(
                        ln, ln, ctrl_cache[ln]["x_state"], q_names, g_eval,
                        ctrl_is_embedded_x_state=True, pert_canon_to_model=p_model_canon_to_key,
                        gene_canon_map=gene_canon_map,
                    )

        # Addendum A Item 3: Arc archived adata_pred.h5ad adapter check (strictly `nadig_hepg2` ONLY)
        if verify_arc_adata_pred_path is not None and ln == ARC_ADAPTER_CHECK_LINE:
            nat_m_bm = None
            if g4_cells:
                for _, cells in sorted(g4_cells.items()):
                    if "M_BM" in cells:
                        nat_m_bm = cells["M_BM"]
                        break
            if nat_m_bm is None:
                nat_ad_chk = ArcStateAdapter(
                    variant="STATE_SE", checkpoint_id=spec["checkpoint_id"],
                    models_dir=models_dir, device=device, batch_mode="native_gem_group", status_dir=status_dir,
                )
                nat_m_bm = nat_ad_chk.predict_deltas(
                    ln, ln, ctrl_cache[ln]["x_state"], q_names, g_eval,
                    ctrl_is_embedded_x_state=True, pert_canon_to_model=p_model_canon_to_key,
                    gene_canon_map=gene_canon_map, native_batch_labels=ctrl_cache[ln]["batch_labels"],
                )
            arc_adapter_check_res = verify_arc_archived_adata_pred_r16(
                adata_pred_path=verify_arc_adata_pred_path,
                adapter_deltas=nat_m_bm,
                adapter_perts=q_names,
                adapter_genes=g_eval,
                line_name=ln,
                status_dir=status_dir,
                hgnc_bridge=bridge,
            )
            (Path(results_dir) / "arc_adapter_check_nadig_hepg2.json").write_text(
                json.dumps(arc_adapter_check_res, indent=2) + "\n", encoding="utf-8"
            )
            if not arc_adapter_check_res["passed"]:
                raise RuntimeError(
                    "Addendum A Item 3 STOP: nadig_hepg2 Arc archived adata_pred.h5ad adapter check FAILED "
                    f"(median_per_pert_pearson={arc_adapter_check_res['median_per_pert_pearson']:.4f} < 0.90 or "
                    f"pooled_pearson={arc_adapter_check_res['pooled_pearson']:.4f} < 0.95). "
                    "Halting for Master Controller adapter diagnosis; adjudication rules unchanged."
                )

        # Addendum A Item 1: Group 3 secondary `FCT` sensitivity (`full r9_6` fixed on primary `N_eval^R16`)
        extra_fct_sens: Optional[Dict[str, Any]] = None
        if "fct_sources_sensitivity" in spec:
            sens_aligned = _align_line_and_sources_r16(
                line_name=ln,
                fct_source_names=spec["fct_sources_sensitivity"],
                line_dir=spec["line_dir"],
                g_model_2000=g_model_2000,
                p_model_canon_to_key=p_model_canon_to_key,
                hgnc_bridge=bridge,
                common_genes_r9=common_genes_r9,
                fixed_target_names=q_names,
                fixed_gene_names=g_eval,
            )
            extra_fct_sens = {
                "FCT_full_r9_6_incl_vcc2025_h1": {
                    "fct_delta_pred": sens_aligned["fct_delta_pred"],
                    "delta_transfer": sens_aligned["delta_transfer"],
                    "source_grand_mean": sens_aligned["source_grand_mean"],
                }
            }

        line_results[ln] = evaluate_line_from_arrays_r16(
            line_name=ln,
            group_key=spec["group"],
            modality=spec["modality"],
            checkpoint_id=spec["checkpoint_id"],
            seen_sources=seen_s,
            unseen_sources=unseen_s,
            delta_obs=aligned["delta_obs"],
            delta_half_a=aligned["delta_half_a"],
            delta_half_b=aligned["delta_half_b"],
            half_valid=aligned["half_valid"],
            delta_transfer=aligned["delta_transfer"],
            source_grand_mean=aligned["source_grand_mean"],
            target_indices=aligned["target_indices"],
            target_names=q_names,
            gene_names=g_eval,
            mu_target=aligned["mu_target"],
            mu_sources=aligned["mu_sources"],
            fct_delta_pred=aligned["fct_delta_pred"],
            e1_constmed_delta_pred=aligned["e1_constmed_delta_pred"],
            state_se_matched=se_matched,
            state_se_mismatched=se_mis,
            state_se_resample=se_resample,
            state_hvg_matched=hvg_matched,
            state_hvg_mismatched=hvg_mis,
            state_hvg_resample=hvg_resample,
            state_hvg_libnorm1e4_matched=hvg_1e4_m,
            state_se_neutral_batch_arms=neutral_arms,
            state_se_per_fold_arms=per_fold_arms,
            g4_factorial_cells=g4_cells,
            g4_single_index_arms=g4_single_idx_arms,
            embedding_diagnostics=emb_diag,
            ctrl_audit=ctrl_cache[ln]["ctrl_audit"],
            control_barcodes_by_source=barcodes_by_source,
            resample_barcodes=ctrl_cache[ln]["resample_barcodes"],
            extra_fct_sensitivity=extra_fct_sens,
            arrays_dir=arrays_dir,
            n_boot=n_boot,
            n_perm_hvg=n_perm_hvg,
            seed=seed,
        )

    # Group-level adjudications (strictly separate Group 1, Group 2, and Group 3)
    groups_out: Dict[str, Any] = {}
    for gkey in ("group1", "group2", "group3"):
        glines = [ln for ln in lines_to_run if LINE_SPECS_R16[ln]["group"] == gkey]
        if not glines:
            continue
        verdicts = {ln: line_results[ln]["primary"]["H_I1_adjudication"]["line_summary"] for ln in glines}
        mods = {ln: LINE_SPECS_R16[ln]["modality"] for ln in glines}
        g_adj = adjudicate_group_r16(gkey, verdicts, mods)
        if gkey == "group3":
            g_adj["mandatory_disclosures"] = [
                "Reported strictly separately from Group 1 and Group 2 (never pooled or merged).",
                "Shared H1 hESC cell lineage with Round 9 vcc2025_h1 (excluded from primary FCT pool).",
                "Shared Replogle/Weissman laboratory provenance with Replogle-4 training cell lines.",
            ]
        groups_out[gkey] = g_adj

    summary = {
        "protocol": "PREREG_I_v5_plus_Addendum_A",
        "seed": int(seed),
        "resample_seed": int(resample_seed),
        "n_boot": int(n_boot),
        "n_perm_hvg": int(n_perm_hvg),
        "delta_tost": float(DEFAULT_TOST_DELTA),
        "cellularai_open_code_public_refit": "DROPPED_PRE_DATA_PER_ADDENDUM_A_ITEM_5",
        "three_party_tolerances": {
            "control_barcode_exact_match": True,
            "se600m_embedding_max_abs_diff": 1e-5,
            "independent_pred_delta_max_abs_diff": 1e-4,
            "controller_same_array_metric_tol": 1e-10,
            "cross_party_metric_tol": 1e-5,
            "verdict_parity": "100%",
        },
        "gates_verified": {k: str(v) for k, v in gate_paths.items()},
        "elapsed_sec": round(float(time.time() - t0), 2),
        "control_barcodes_file": str(barcodes_json_path),
        "control_source_barcodes": control_source_barcodes,
        "embeddings_dir": str(embeddings_dir),
        "embeddings_manifest_sha256": embeddings_manifest,
        "arc_adapter_check_nadig_hepg2": arc_adapter_check_res,
        "group_adjudications": groups_out,
        "lines": line_results,
    }
    out_json = Path(results_dir) / "summary_r16.json"
    out_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


# ==============================================================================
# Pure-Synthetic End-to-End Smoke Runner (ZERO Real Cells Touched)
# ==============================================================================


def run_synthetic_smoke_r16(
    out_dir: Path,
    n_boot: int = 300,
    seed: int = DEFAULT_R16_SEED,
) -> Dict[str, Any]:
    """
    Execute the complete `run_eval_r16` evaluation, `H-I1`/`H-I2`/`r_dorth_cal`/`M2`/`Perm100`/
    `G4` 2x2 + `(17,0)` symmetry / `{0,1,17,35}` spread / Item 8 within-line resample / barcode
    export pipeline strictly on 100% synthetic tensors (ZERO real cells touched).
    """
    rng = np.random.default_rng(seed)
    P, G = 60, 240
    gene_names = [f"SYN_G{j:04d}" for j in range(G)]
    target_names = [f"SYN_G{i:04d}" for i in range(P)]
    target_indices = np.arange(P, dtype=np.int64)
    half_valid = np.ones(P, dtype=bool)

    b_s = rng.normal(0.0, 0.12, size=G)
    tm = b_s[None, :] + rng.normal(0.0, 0.18, size=(P, G))
    true_spec = rng.normal(0.0, 0.28, size=(P, G))
    ha = tm + true_spec + rng.normal(0.0, 0.04, size=(P, G))
    hb = tm + true_spec + rng.normal(0.0, 0.04, size=(P, G))
    obs = 0.5 * (ha + hb)

    mu_sources = np.exp(rng.normal(0.0, 0.8, size=(5, G)))
    mu_target = np.exp(rng.normal(0.2, 0.8, size=G))
    fc_bar = np.exp(rng.normal(0.0, 0.25, size=(P, G)))
    has_hits = np.ones(P, dtype=bool)
    fct_pred = fct_from_fc_bar(mu_target, fc_bar, tm, has_hits)
    e1_pred = fct_from_fc_bar(np.full(G, float(np.median(mu_target))), fc_bar, tm, has_hits)

    # Synthetic Group 2 primary line (`a549_ko` schema, matched-specific signal)
    strata = build_mismatched_strata_r16("a549_ko", "zeroshot/ensemble4")
    all_mis = strata["seen_by_checkpoint"] + strata["unseen_by_checkpoint"]
    se_m = tm + 0.75 * true_spec + rng.normal(0.0, 0.05, size=(P, G))
    se_resamp = tm + 0.74 * true_spec + rng.normal(0.0, 0.052, size=(P, G))
    se_mis = {sk: tm + 0.05 * true_spec + rng.normal(0.0, 0.20, size=(P, G)) for sk in all_mis}
    hvg_m = tm + 0.55 * true_spec + rng.normal(0.0, 0.08, size=(P, G))
    hvg_resamp = tm + 0.54 * true_spec + rng.normal(0.0, 0.082, size=(P, G))
    hvg_mis = {sk: tm + 0.05 * true_spec + rng.normal(0.0, 0.22, size=(P, G)) for sk in all_mis}
    hvg_1e4 = hvg_m + rng.normal(0.0, 0.01, size=(P, G))

    m_b17 = tm + 0.65 * true_spec + rng.normal(0.0, 0.06, size=(P, G))
    g4_cells: Dict[str, Dict[str, np.ndarray]] = {}
    for sk in all_mis:
        s_b17 = tm + 0.02 * true_spec + rng.normal(0.0, 0.22, size=(P, G))
        g4_cells[f"b0_vs_b17_mis_{sk}"] = {"M_BM": se_m, "M_BS": m_b17, "S_BM": se_mis[sk], "S_BS": s_b17}
        g4_cells[f"b17_vs_b0_sym_mis_{sk}"] = {"M_BM": m_b17, "M_BS": se_m, "S_BM": s_b17, "S_BS": se_mis[sk]}

    g4_single_idx_arms = {
        "b0": se_m,
        "b1": tm + 0.70 * true_spec + rng.normal(0.0, 0.06, size=(P, G)),
        "b17": m_b17,
        "b35": tm + 0.68 * true_spec + rng.normal(0.0, 0.06, size=(P, G)),
    }
    neutral_arms = {
        "bypass_zero": (se_m + rng.normal(0.0, 0.01, size=(P, G)), se_mis),
    }

    # Synthetic barcode lists (primary 512 + disjoint 512 resample, 0 overlap)
    syn_bcs_primary = [f"SYN_A549_P_{i:04d}" for i in range(512)]
    syn_bcs_resample = [f"SYN_A549_R_{i:04d}" for i in range(512)]
    bc_by_src = {"synthetic_a549_ko": syn_bcs_primary}
    for sk in all_mis:
        bc_by_src[sk] = [f"SYN_{sk.upper()}_{i:04d}" for i in range(512)]

    ctrl_audit_syn = {
        "n_ctrl_total": 2048,
        "n_ctrl_sampled": 512,
        "ctrl_sampling_fraction": 0.25,
        "sampled_barcodes": syn_bcs_primary,
        "within_line_resample": {
            "eligible": True,
            "regime": "full_512_disjoint",
            "n_resample": 512,
            "seed": DEFAULT_RESAMPLE_SEED,
            "overlap_with_primary": 0,
            "sampled_barcodes": syn_bcs_resample,
        },
    }

    arrays_dir = Path(out_dir) / "arrays"
    line_res = evaluate_line_from_arrays_r16(
        line_name="synthetic_a549_ko",
        group_key="group2",
        modality="CRISPRko",
        checkpoint_id="zeroshot/ensemble4",
        seen_sources=strata["seen_by_checkpoint"],
        unseen_sources=strata["unseen_by_checkpoint"],
        delta_obs=obs,
        delta_half_a=ha,
        delta_half_b=hb,
        half_valid=half_valid,
        delta_transfer=tm,
        source_grand_mean=b_s,
        target_indices=target_indices,
        target_names=target_names,
        gene_names=gene_names,
        mu_target=mu_target,
        mu_sources=mu_sources,
        fct_delta_pred=fct_pred,
        e1_constmed_delta_pred=e1_pred,
        state_se_matched=se_m,
        state_se_mismatched=se_mis,
        state_se_resample=se_resamp,
        state_hvg_matched=hvg_m,
        state_hvg_mismatched=hvg_mis,
        state_hvg_resample=hvg_resamp,
        state_hvg_libnorm1e4_matched=hvg_1e4,
        state_se_neutral_batch_arms=neutral_arms,
        g4_factorial_cells=g4_cells,
        g4_single_index_arms=g4_single_idx_arms,
        ctrl_audit=ctrl_audit_syn,
        control_barcodes_by_source=bc_by_src,
        resample_barcodes=syn_bcs_resample,
        arrays_dir=arrays_dir,
        n_boot=n_boot,
        n_perm_hvg=100,
        seed=seed,
    )

    # Addendum A Item 11: export synthetic `2,058`-dim `X_state` embeddings (primary + resample + mismatched)
    rng_emb = np.random.default_rng(seed)

    def _syn_x_state(n_rows: int) -> np.ndarray:
        e2048 = rng_emb.normal(0.0, 1.0, size=(n_rows, 2048)).astype(np.float32)
        e2048 /= np.maximum(np.linalg.norm(e2048, axis=1, keepdims=True), 1e-12)
        d10 = rng_emb.normal(0.0, 0.15, size=(n_rows, 10)).astype(np.float32)
        return np.ascontiguousarray(np.concatenate([e2048, d10], axis=1), dtype=np.float32)

    def _syn_emb_entry(ds_key: str, bcs: List[str], pool_seed: int) -> Dict[str, Any]:
        x_arr = _syn_x_state(len(bcs))
        return {
            "X_state": x_arr,
            "control_barcodes": bcs,
            "seed": int(pool_seed),
            "dataset_name": ds_key,
            "se_audit": {
                "device": "synthetic",
                "dtype": "float32",
                "seed": int(pool_seed),
                "batch_size": 16,
                "shuffle": False,
                "num_workers": 0,
                "collator_training": False,
                "n_cells": len(bcs),
                "embed_dim": 2058,
                "all_finite": True,
                "x_state_sha256": hashlib.sha256(x_arr.tobytes()).hexdigest(),
            },
        }

    syn_ctrl_embeddings: Dict[str, Dict[str, Any]] = {
        "synthetic_a549_ko": _syn_emb_entry("synthetic_a549_ko", syn_bcs_primary, int(seed)),
        "synthetic_a549_ko_resample": _syn_emb_entry(
            "synthetic_a549_ko_resample", syn_bcs_resample, int(DEFAULT_RESAMPLE_SEED)
        ),
    }
    for sk in sorted(all_mis):
        syn_ctrl_embeddings[sk] = _syn_emb_entry(sk, bc_by_src[sk], int(seed))
    embeddings_manifest = save_control_embeddings_r16(Path(out_dir) / "embeddings", syn_ctrl_embeddings)

    summary = {
        "mode": "synthetic_smoke",
        "protocol": "PREREG_I_v5_plus_Addendum_A",
        "real_cells_touched": 0,
        "embeddings_manifest_sha256": embeddings_manifest,
        "control_source_barcodes": {
            k: {
                "barcodes": v,
                "within_line_resample_barcodes": syn_bcs_resample if k == "synthetic_a549_ko" else [],
            }
            for k, v in bc_by_src.items()
        },
        "lines": {"synthetic_a549_ko": line_res},
    }
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    (Path(out_dir) / "control_barcodes_synthetic_smoke.json").write_text(
        json.dumps(summary["control_source_barcodes"], indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (Path(out_dir) / "summary_synthetic_smoke.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Round 16 Evaluation Runner (PREREG_I v5 + Addendum A).")
    p.add_argument(
        "--synthetic-smoke",
        action="store_true",
        help="Run end-to-end pipeline + array export strictly on synthetic arrays (zero real cells).",
    )
    p.add_argument(
        "--cohort",
        choices=["groups12", "group1", "group2", "group3", "sensitivity", "all_local"],
        default="groups12",
        help="Which cohort set to evaluate in real-data Gate 5 mode.",
    )
    p.add_argument("--lines", nargs="*", default=None, help="Explicit subset of line keys to run.")
    p.add_argument("--status-dir", type=Path, default=DEFAULT_STATUS_DIR)
    p.add_argument("--models-dir", type=Path, default=DEFAULT_MODELS_DIR)
    p.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    p.add_argument("--gse295214-raw-dir", type=Path, default=None)
    p.add_argument(
        "--raw-nadig-dir",
        type=Path,
        default=None,
        help="Directory holding GSE264667_{hepg2,jurkat}_raw_singlecell_01.h5ad (overrides NADIG_DATA_DIR / sync drive).",
    )
    p.add_argument(
        "--verify-arc-adata-pred",
        type=Path,
        default=None,
        help="Post-R16_GO path to archived nadig_hepg2 adata_pred.h5ad for Addendum A Item 3 adapter check.",
    )
    p.add_argument("--device", default="mps", choices=["mps", "cpu"])
    p.add_argument("--n-boot", type=int, default=4000)
    p.add_argument("--n-perm-hvg", type=int, default=100)
    p.add_argument("--seed", type=int, default=DEFAULT_R16_SEED)
    p.add_argument("--resample-seed", type=int, default=DEFAULT_RESAMPLE_SEED)
    p.add_argument("--skip-secondary-panels", action="store_true")
    p.add_argument("--smoke-out-dir", type=Path, default=None)
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.synthetic_smoke:
        if args.smoke_out_dir is not None:
            res = run_synthetic_smoke_r16(args.smoke_out_dir, n_boot=min(args.n_boot, 300), seed=args.seed)
        else:
            with tempfile.TemporaryDirectory(prefix="r16_synth_smoke_") as tmpd:
                res = run_synthetic_smoke_r16(Path(tmpd), n_boot=min(args.n_boot, 300), seed=args.seed)
        line_r = res["lines"]["synthetic_a549_ko"]
        print(
            json.dumps(
                {
                    "mode": "synthetic_smoke",
                    "protocol": res["protocol"],
                    "real_cells_touched": 0,
                    "h_i1_line_summary": line_r["primary"]["H_I1_adjudication"]["line_summary"],
                    "h_i2_classification": line_r["primary"]["H_I2_adjudication"]["classification"],
                    "within_line_resample_status": (
                        line_r["secondary"]["within_line_resample_reference"]["descriptive_interpretation"][
                            "resample_comparison_status"
                        ]
                    ),
                    "hvg_perm100_n_finite": line_r["secondary"]["STATE_HVG"]["perm100_r_orth_cal"]["n_finite_perms"],
                    "n_exported_arrays": len(line_r["array_manifest_sha256"]),
                    "n_exported_embeddings": len(res["embeddings_manifest_sha256"]),
                },
                indent=2,
            )
        )
        return 0

    # STRICT TRIPLE HARD GATE BEFORE ANY REAL-DATA EXECUTION
    assert_all_execution_gates_r16(args.status_dir)

    if args.lines:
        target_lines = list(args.lines)
    elif args.cohort == "group1":
        target_lines = list(REPLOGLE4_LOO_PRIMARY_LINES)
    elif args.cohort == "group2":
        target_lines = list(GROUP2_UNSEEN_PRIMARY_LINES)
    elif args.cohort == "groups12":
        target_lines = list(REPLOGLE4_LOO_PRIMARY_LINES) + list(GROUP2_UNSEEN_PRIMARY_LINES)
    elif args.cohort == "group3":
        target_lines = list(GROUP3_PROSPECTIVE_LINES)
    elif args.cohort == "sensitivity":
        target_lines = list(RETAINED_UNSEEN_MISMATCHED_ONLY_LINES)
    else:
        target_lines = (
            list(REPLOGLE4_LOO_PRIMARY_LINES)
            + list(GROUP2_UNSEEN_PRIMARY_LINES)
            + list(RETAINED_UNSEEN_MISMATCHED_ONLY_LINES)
        )

    run_real_evaluation_r16(
        lines_to_run=target_lines,
        status_dir=args.status_dir,
        models_dir=args.models_dir,
        results_dir=args.results_dir,
        device=args.device,
        n_boot=args.n_boot,
        n_perm_hvg=args.n_perm_hvg,
        seed=args.seed,
        resample_seed=args.resample_seed,
        include_secondary_panels=not args.skip_secondary_panels,
        gse295214_raw_dir=args.gse295214_raw_dir,
        raw_nadig_dir=args.raw_nadig_dir,
        verify_arc_adata_pred_path=args.verify_arc_adata_pred,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
