#!/usr/bin/env python3
"""
Round 17 Pre-Download Dataset Reader for Pfizer TeloHAEC Perturb-seq (`Zenodo 18792213`,
`imru_full.h5ad`, Wong et al. 2026, `bioRxiv 10.64898/2026.03.04.709607`, CC BY 4.0).

Built strictly BEFORE download and verified on synthetic `.h5ad` fixtures
(`round17/src/test_pfizer_telohaec_reader_r17.py`). Zero model imports (`torch`/`state`),
zero downloads, zero correlation/metric calculations.

Locked ingestion & QC specification (`PREREG_J_DRAFT_v2.md §2.3`):
  1. Strict Pre-Metric Raw Integer UMI Verification (`EXCLUDED_FEASIBILITY` / `SLOT_AMBIGUOUS`):
     - Reuses frozen Round 14 `_resolve_h5ad_integer_slot` (`X`, non-excluded `layers/*`, `raw/X`)
       to verify 100% of stored non-zero entries are finite non-negative exact integers.
     - Streams `CSR`, `CSC`, or `dense` HDF5 layouts in bounded-RAM chunks (`chunk_rows=4096`)
       via `h5py` + `scipy.sparse` (never loads the 53.29 GB matrix densely into RAM).
  2. Condition Separation (Untreated Basal Primary vs. IL-1b+TNF-a Secondary):
     - `condition == "untreated"` (`pfizer_telohaec_basal`): primary confirmatory unit.
     - `condition == "treated"` (`pfizer_telohaec_cytokine`, IL-1b + TNF-a): strictly secondary
       within-family replication (same cells under a second condition; NEVER a second independent line).
  3. Single-sgRNA & Strict NTC vs. Safe-Targeting Separation:
     - Uses `.obs['gene_target']` & `.obs['guide']` (supporting modern & legacy AnnData categoricals).
     - Strict NTC: `NO-TARGET` (and canonical non-targeting aliases) ONLY.
     - Strictly excludes `SAFE_TARGET` (never pooled into NTC or targets), unassigned/NA/multiplet
       cells, and multi-guide delimiter strings (`[,;|+]`).
     - Resolves plain single-gene targets through frozen `HGNCLookup.resolve` (approved HGNC only);
       retains per-guide `pert_key` for frozen `process_r10_line` per-guide Addendum-B QC + dedup.
  4. Cell RNA QC (`PREREG_G` / `PREREG_I §2.4` Parity):
     - Retain cells with `umi >= 300` AND `ngenes >= 200` over human RNA features (`max_pct_mt = None`).
  5. Disjoint Split-Half & Hardware Batch Centering:
     - Primary hardware unit: `flask` (if `>= 2` distinct flasks each have `>= 15` QC-passing NTCs),
       else `welltag` (if `>= 2` distinct welltags each have `>= 15` NTCs), else deterministic
       stratified 50/50 fallback (`seed = 20261701`, explicitly flagged).
     - Lexicographical 1-based odd -> Half A (`1`), even -> Half B (`2`), with within-half
       `< 15` NTC batch fallback (`_apply_min_ntc_batch_fallback(..., min_batch_nt=15)`).
  6. Counts-Only Qualifying Gate (`prepare_pfizer_telohaec_r17`):
     - Applies exact `PREREG_I` floors (`N_eval >= 50` vs. `r9_6 ∩ canon(P_model)`,
       `|G_eval^R9| >= 3,000`, `per-half NTC >= 15`) and logs counts-only QC strictly pre-inference.
"""

from __future__ import annotations

import hashlib
import json
import pickle
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import h5py
import numpy as np

ROUND17_SRC_DIR = Path(__file__).resolve().parent
REPO_ROOT = ROUND17_SRC_DIR.parents[1]
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
R10_SCRIPTS_DIR = REPO_ROOT / "round10" / "shared" / "scripts"
R11_SRC_DIR = REPO_ROOT / "round11" / "src"
R12_SRC_DIR = REPO_ROOT / "round12" / "src"
R14_SRC_DIR = REPO_ROOT / "round14" / "src"
for _p in (str(R10_SRC_DIR), str(R10_SCRIPTS_DIR), str(R11_SRC_DIR), str(R12_SRC_DIR), str(R14_SRC_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from data_loader_r10 import (  # noqa: E402
    DEFAULT_HGNC_FILE,
    HGNCBridge,
    load_r9_common_genes,
    load_standardized_line_r10,
)
from dataset_readers_r11 import (  # noqa: E402
    _apply_min_ntc_batch_fallback,
    _is_plain_single_gene_symbol,
)
from dataset_readers_r12 import (  # noqa: E402
    _is_nonhuman_or_vector_var,
    _resolve_h5ad_index_col,
)
from dataset_readers_r14 import (  # noqa: E402
    _get_h5ad_slot_node,
    _get_var_group_for_slot,
    _inspect_slot_shape_and_layout,
    _resolve_h5ad_integer_slot,
)
from prepare_r10_lines import (  # noqa: E402
    HGNCLookup,
    PartSpec,
    build_hgnc_lookup,
    process_r10_line,
)

# Locked Zenodo 18792213 provenance constants
PFIZER_TELOHAEC_ZENODO_ID = "18792213"
PFIZER_TELOHAEC_FILENAME = "imru_full.h5ad"
PFIZER_TELOHAEC_EXPECTED_BYTES = 53_290_123_849
PFIZER_TELOHAEC_EXPECTED_MD5 = "0594ddd25e224e5e9cc22414e71daf48"
PFIZER_TELOHAEC_BIORXIV_DOI = "10.64898/2026.03.04.709607"
PFIZER_TELOHAEC_LICENSE = "CC-BY-4.0"
PFIZER_TELOHAEC_RELEASE_DATE = "2026-03-04"

DEFAULT_R17_SEED = 20261701
DEFAULT_R17_MIN_UMI = 300.0
DEFAULT_R17_MIN_NGENES = 200
DEFAULT_R17_MIN_BATCH_NT = 15

R9_LINES_DIR = REPO_ROOT / "round9" / "p2" / "shared" / "lines"
R9_COMMON_GENES_FILE = REPO_ROOT / "round9" / "p2" / "shared" / "COMMON_GENES.txt"
R16_VAR_DIMS_PKL = (
    REPO_ROOT
    / "round16"
    / "shared"
    / "models"
    / "ST-SE-Replogle"
    / "zeroshot"
    / "hepg2"
    / "var_dims.pkl"
)
R9_SOURCE_POOL_6: Tuple[str, ...] = (
    "nadig_hepg2",
    "nadig_jurkat",
    "nourreddine_kolf21j",
    "replogle_k562_essential",
    "replogle_rpe1",
    "vcc2025_h1",
)

_MULTI_GUIDE_RE = re.compile(r"[,;|+]")
_UNASSIGNED_TOKENS = frozenset(
    {
        "",
        "na",
        "n/a",
        "<na>",
        "nan",
        "none",
        "null",
        "unassigned",
        "multiplet",
        "doublet",
        "ambiguous",
        "unknown",
        "multiple",
    }
)
_NTC_TOKENS = frozenset(
    {
        "no-target",
        "no_target",
        "notarget",
        "non-targeting",
        "nontargeting",
        "non_targeting",
        "non-target",
        "ntc",
        "ctrl_nt",
    }
)
_SAFE_TOKENS = frozenset(
    {
        "safe_target",
        "safe-target",
        "safetarget",
        "safe_targeting",
        "safe-targeting",
        "safetargeting",
        "safe",
        "intergenic",
        "gene_desert",
        "ctrl_t",
    }
)
_UNTREATED_COND_TOKENS = frozenset({"untreated", "basal", "unstimulated", "ctrl", "control", "vehicle", "dmso"})
_TREATED_COND_TOKENS = frozenset(
    {"treated", "il1b_tnfa", "il-1b+tnf-a", "il1b+tnfa", "il1b_tnf", "cytokine", "stimulated", "inflamed"}
)


def read_h5ad_table_col_r17(grp: h5py.Group, col: str) -> List[str]:
    """
    Decode an `.obs` or `.var` column from an HDF5 AnnData group cleanly across:
      1. Modern categorical group (`grp[col]/{codes, categories}`),
      2. Legacy AnnData categorical (`grp['__categories'][col]` + integer codes in `grp[col]`),
      3. Plain string/bytes/numeric dataset.
    Never returns raw integer codes for categorical columns.
    """
    if col not in grp:
        raise KeyError(f"Column {col!r} not found in HDF5 group {grp.name!r}; keys={sorted(grp.keys())}")
    node = grp[col]
    if isinstance(node, h5py.Group) and "codes" in node and "categories" in node:
        codes = np.asarray(node["codes"][:], dtype=np.int64)
        cats_raw = node["categories"][:]
        cats = [
            c.decode("utf-8", errors="replace") if isinstance(c, (bytes, bytearray, np.bytes_)) else str(c)
            for c in cats_raw
        ]
        return [cats[int(c)] if 0 <= int(c) < len(cats) else "" for c in codes]

    arr = node[:]
    if "__categories" in grp and isinstance(grp["__categories"], h5py.Group) and col in grp["__categories"]:
        if np.issubdtype(arr.dtype, np.integer):
            cats_raw = grp["__categories"][col][:]
            cats = [
                c.decode("utf-8", errors="replace") if isinstance(c, (bytes, bytearray, np.bytes_)) else str(c)
                for c in cats_raw
            ]
            return [cats[int(c)] if 0 <= int(c) < len(cats) else "" for c in arr]

    out: List[str] = []
    for val in arr:
        if isinstance(val, (bytes, bytearray, np.bytes_)):
            out.append(val.decode("utf-8", errors="replace"))
        else:
            out.append(str(val))
    return out


def _pick_obs_col_r17(obs_keys: Sequence[str], preferred: Sequence[str]) -> Optional[str]:
    key_set = set(obs_keys)
    for cand in preferred:
        if cand in key_set:
            return cand
    low_map = {k.lower(): k for k in sorted(obs_keys)}
    for cand in preferred:
        if cand.lower() in low_map:
            return low_map[cand.lower()]
    return None


def classify_pfizer_condition_r17(raw_cond: str) -> str:
    """Normalize `condition` string to `'untreated'`, `'treated'`, or `'other'`."""
    s = str(raw_cond).strip().lower()
    if s in _UNTREATED_COND_TOKENS or s.startswith("untreat") or s.startswith("basal"):
        return "untreated"
    if s in _TREATED_COND_TOKENS or ("il" in s and "tnf" in s) or s.startswith("treat"):
        return "treated"
    return "other"


def classify_pfizer_guide_r17(
    gene_target: str,
    guide_id: str,
    hgnc: HGNCLookup,
) -> Tuple[str, Optional[Tuple[str, str, str, str]], str]:
    """
    Classify one cell's `(gene_target, guide_id)` pair strictly into:
      - `("ntc", None, guide_key)` for strict `NO-TARGET` / non-targeting controls ONLY;
      - `("safe", None, "")` for `SAFE_TARGET` / safe-targeting / intergenic controls (EXCLUDED);
      - `("target", (ens_id, canon_sym, raw_sym, hgnc_status), guide_key)` for valid plain single-gene targets;
      - `("exclude", None, "")` for unassigned, multi-guide, or non-approved symbols.
    """
    gt = str(gene_target).strip()
    gd = str(guide_id).strip()
    gt_l = gt.lower()
    gd_l = gd.lower()

    if gt_l in _UNASSIGNED_TOKENS:
        return "exclude", None, ""
    if _MULTI_GUIDE_RE.search(gt) or (gd and _MULTI_GUIDE_RE.search(gd)):
        return "exclude", None, ""

    # Safe-targeting check BEFORE NTC/target so safe controls never enter NTC
    for tok in (gt_l, gd_l):
        if not tok:
            continue
        if (
            tok in _SAFE_TOKENS
            or tok.startswith(("safe_target", "safe-target", "safetarget", "safe_"))
            or "gene_desert" in tok
            or "intergenic" in tok
        ):
            return "safe", None, ""

    # Strict NTC check
    if (
        gt_l in _NTC_TOKENS
        or gt_l.startswith(("no-target", "no_target", "non-targeting", "nontargeting", "non_targeting"))
        or gd_l.startswith(("no-target", "no_target", "non-targeting", "nontargeting", "non_targeting"))
    ):
        return "ntc", None, (gd or gt)

    # Plain single-gene target resolution via frozen HGNC lookup
    if not _is_plain_single_gene_symbol(gt):
        return "exclude", None, ""
    canon_sym, hstat = hgnc.resolve(gt, "")
    if canon_sym not in hgnc.approved:
        return "exclude", None, ""
    ens_id = hgnc.sym_to_ens.get(canon_sym, "")
    gkey = gd if (gd and gd_l not in _UNASSIGNED_TOKENS) else canon_sym
    return "target", (ens_id, canon_sym, gt, hstat), gkey


def _resolve_h5ad_path_r17(raw_path: Path) -> Path:
    p = Path(raw_path)
    if p.is_file():
        return p
    if p.is_dir():
        pref = p / PFIZER_TELOHAEC_FILENAME
        if pref.is_file():
            return pref
        cands = sorted(q for q in p.rglob("*.h5ad") if q.is_file() and not q.name.startswith("."))
        if len(cands) == 1:
            return cands[0]
        if cands:
            raise RuntimeError(f"Multiple .h5ad files found under {p}: {[c.name for c in cands]}")
    raise FileNotFoundError(f"Pfizer TeloHAEC .h5ad not found at {raw_path}")


def _extract_var_rna_features_r17(var_grp: h5py.Group, n_var_expected: int) -> Dict[str, Any]:
    """Identify human RNA features and exclude non-human/vector/CRISPR-capture features."""
    var_keys = sorted(k for k in var_grp.keys() if k != "__categories")
    idx_col = _resolve_h5ad_index_col(var_grp, group_name="var")
    raw_idx = [x.strip() for x in read_h5ad_table_col_r17(var_grp, idx_col)]
    if len(raw_idx) != n_var_expected:
        raise ValueError(f"var[{idx_col!r}] length {len(raw_idx)} != matrix cols {n_var_expected}")

    id_col = _pick_obs_col_r17(var_keys, ["ensembl_id", "gene_ids", "gene_id", "feature_id", "ensembl"]) or idx_col
    sym_col = _pick_obs_col_r17(var_keys, ["gene_symbol", "gene_name", "symbol", "feature_name", "gene"]) or idx_col
    ft_col = _pick_obs_col_r17(var_keys, ["feature_types", "feature_type", "modality", "type"])

    raw_ids = [x.strip() for x in read_h5ad_table_col_r17(var_grp, id_col)]
    raw_syms = [x.strip() for x in read_h5ad_table_col_r17(var_grp, sym_col)]
    feat_types = [x.strip() for x in read_h5ad_table_col_r17(var_grp, ft_col)] if ft_col else [""] * n_var_expected

    rna_cols: List[int] = []
    rna_ids: List[str] = []
    rna_syms: List[str] = []
    n_non_rna = 0
    for j in range(n_var_expected):
        ft = feat_types[j].lower()
        if ft and ft not in ("gene expression", "rna", "gene"):
            n_non_rna += 1
            continue
        if _is_nonhuman_or_vector_var(raw_ids[j], raw_syms[j]):
            n_non_rna += 1
            continue
        rna_cols.append(j)
        rna_ids.append(raw_ids[j])
        rna_syms.append(raw_syms[j])

    if not rna_cols:
        raise RuntimeError("EXCLUDED_FEASIBILITY: zero human RNA features found in .var")

    return {
        "rna_cols": np.asarray(rna_cols, dtype=np.int64),
        "rna_ids": rna_ids,
        "rna_syms": rna_syms,
        "n_var_total": n_var_expected,
        "n_non_rna_excluded": int(n_non_rna),
        "id_col": id_col,
        "sym_col": sym_col,
    }


def _compute_rna_umi_and_ngenes_r17(
    h5_path: Path,
    slot: str,
    row_indices: np.ndarray,
    rna_cols: np.ndarray,
    chunk_rows: int = 4096,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Bounded-RAM per-cell total UMI and detected gene count strictly over `rna_cols`
    for the subset `row_indices` (ascending int64).
    """
    n_sub = int(row_indices.size)
    umi = np.zeros(n_sub, dtype=np.float64)
    ngenes = np.zeros(n_sub, dtype=np.int64)
    if n_sub == 0:
        return umi, ngenes

    with h5py.File(h5_path, "r") as f:
        obs_grp = f["obs"]
        n_obs = len(read_h5ad_table_col_r17(obs_grp, _resolve_h5ad_index_col(obs_grp, group_name="obs")))
        var_grp = _get_var_group_for_slot(f, slot)
        n_var = len(read_h5ad_table_col_r17(var_grp, _resolve_h5ad_index_col(var_grp, group_name="var")))
        node = _get_h5ad_slot_node(f, slot)
        shape, layout = _inspect_slot_shape_and_layout(node, n_obs, n_var)
        rna_mask = np.zeros(shape[1], dtype=bool)
        rna_mask[rna_cols] = True

        if layout == "csr":
            indptr = np.asarray(node["indptr"][:], dtype=np.int64)
            data_ds = node["data"]
            idx_ds = node["indices"]
            for s in range(0, n_sub, chunk_rows):
                e = min(s + chunk_rows, n_sub)
                sub_rows = row_indices[s:e]
                # Read contiguous span covering sub_rows[0]..sub_rows[-1] when compact, else per-row
                r0, r1 = int(sub_rows[0]), int(sub_rows[-1]) + 1
                if (r1 - r0) <= 4 * (e - s):
                    p0, p1 = int(indptr[r0]), int(indptr[r1])
                    dat = np.asarray(data_ds[p0:p1], dtype=np.float64)
                    idx = np.asarray(idx_ds[p0:p1], dtype=np.int64)
                    for k, r in enumerate(sub_rows):
                        a = int(indptr[r]) - p0
                        b = int(indptr[r + 1]) - p0
                        if b <= a:
                            continue
                        d_r = dat[a:b]
                        i_r = idx[a:b]
                        m = rna_mask[i_r] & (d_r > 0.0)
                        if np.any(m):
                            umi[s + k] = float(np.sum(d_r[m]))
                            ngenes[s + k] = int(np.sum(m))
                else:
                    for k, r in enumerate(sub_rows):
                        a, b = int(indptr[r]), int(indptr[r + 1])
                        if b <= a:
                            continue
                        d_r = np.asarray(data_ds[a:b], dtype=np.float64)
                        i_r = np.asarray(idx_ds[a:b], dtype=np.int64)
                        m = rna_mask[i_r] & (d_r > 0.0)
                        if np.any(m):
                            umi[s + k] = float(np.sum(d_r[m]))
                            ngenes[s + k] = int(np.sum(m))
        elif layout == "csc":
            indptr = np.asarray(node["indptr"][:], dtype=np.int64)
            data_ds = node["data"]
            idx_ds = node["indices"]
            row_to_pos = np.full(shape[0], -1, dtype=np.int64)
            row_to_pos[row_indices] = np.arange(n_sub, dtype=np.int64)
            for c_start in range(0, len(rna_cols), 1024):
                cols_blk = rna_cols[c_start : c_start + 1024]
                for c in cols_blk:
                    p0, p1 = int(indptr[c]), int(indptr[c + 1])
                    if p1 <= p0:
                        continue
                    r_col = np.asarray(idx_ds[p0:p1], dtype=np.int64)
                    pos = row_to_pos[r_col]
                    keep = pos >= 0
                    if not np.any(keep):
                        continue
                    d_col = np.asarray(data_ds[p0:p1], dtype=np.float64)[keep]
                    pos_k = pos[keep]
                    nz = d_col > 0.0
                    if np.any(nz):
                        np.add.at(umi, pos_k[nz], d_col[nz])
                        np.add.at(ngenes, pos_k[nz], 1)
        else:
            for s in range(0, n_sub, 1024):
                e = min(s + 1024, n_sub)
                sub_rows = row_indices[s:e]
                blk = np.asarray(node[sub_rows.tolist(), :], dtype=np.float64)[:, rna_cols]
                umi[s:e] = np.sum(blk, axis=1)
                ngenes[s:e] = np.sum(blk > 0.0, axis=1)

    return umi, ngenes


def _make_h5ad_subset_readers_r17(
    h5_path: Path,
    slot: str,
    keep_rows: np.ndarray,
    rna_cols: np.ndarray,
    n_obs_total: int,
    n_var_total: int,
    layout: str,
) -> Tuple[Any, Any, Any]:
    """Build bounded-RAM `(read_chunk_fn, read_rows_fn, read_chunk_cols_fn)` for `PartSpec`."""
    keep_rows_arr = np.asarray(keep_rows, dtype=np.int64)
    rna_cols_arr = np.asarray(rna_cols, dtype=np.int64)
    n_rna = int(rna_cols_arr.size)

    # Precompute map from raw column index -> dense RNA column position [0..n_rna-1]
    col_to_rna = np.full(n_var_total, -1, dtype=np.int64)
    col_to_rna[rna_cols_arr] = np.arange(n_rna, dtype=np.int64)

    def _read_selected_rows(abs_rows: np.ndarray, out_cols: Optional[np.ndarray] = None, dtype: Any = np.float64) -> np.ndarray:
        abs_rows = np.asarray(abs_rows, dtype=np.int64)
        n_r = int(abs_rows.size)
        sel_rna_cols = rna_cols_arr if out_cols is None else rna_cols_arr[np.asarray(out_cols, dtype=np.int64)]
        n_c = int(sel_rna_cols.size)
        out = np.zeros((n_r, n_c), dtype=dtype)
        if n_r == 0 or n_c == 0:
            return out

        if out_cols is None:
            c_map = col_to_rna
        else:
            c_map = np.full(n_var_total, -1, dtype=np.int64)
            c_map[sel_rna_cols] = np.arange(n_c, dtype=np.int64)

        with h5py.File(h5_path, "r") as f:
            node = _get_h5ad_slot_node(f, slot)
            if layout == "csr":
                indptr = np.asarray(node["indptr"][:], dtype=np.int64)
                data_ds = node["data"]
                idx_ds = node["indices"]
                for k, r in enumerate(abs_rows):
                    a, b = int(indptr[r]), int(indptr[r + 1])
                    if b <= a:
                        continue
                    idx_r = np.asarray(idx_ds[a:b], dtype=np.int64)
                    mapped = c_map[idx_r]
                    m = mapped >= 0
                    if np.any(m):
                        dat_r = np.asarray(data_ds[a:b], dtype=dtype)[m]
                        out[k, mapped[m]] = dat_r
            elif layout == "csc":
                indptr = np.asarray(node["indptr"][:], dtype=np.int64)
                data_ds = node["data"]
                idx_ds = node["indices"]
                r_map = np.full(n_obs_total, -1, dtype=np.int64)
                r_map[abs_rows] = np.arange(n_r, dtype=np.int64)
                for j_out, c_raw in enumerate(sel_rna_cols):
                    a, b = int(indptr[c_raw]), int(indptr[c_raw + 1])
                    if b <= a:
                        continue
                    r_col = np.asarray(idx_ds[a:b], dtype=np.int64)
                    pos = r_map[r_col]
                    m = pos >= 0
                    if np.any(m):
                        out[pos[m], j_out] = np.asarray(data_ds[a:b], dtype=dtype)[m]
            else:
                for k, r in enumerate(abs_rows):
                    row_v = np.asarray(node[int(r), :], dtype=dtype)
                    out[k, :] = row_v[sel_rna_cols]
        return out

    def read_chunk_fn(start: int, end: int) -> np.ndarray:
        return _read_selected_rows(keep_rows_arr[start:end], out_cols=None, dtype=np.float64)

    def read_rows_fn(local_rows: np.ndarray) -> np.ndarray:
        return _read_selected_rows(keep_rows_arr[np.asarray(local_rows, dtype=np.int64)], out_cols=None, dtype=np.float64)

    def read_chunk_cols_fn(start: int, end: int, cols: np.ndarray) -> np.ndarray:
        return _read_selected_rows(keep_rows_arr[start:end], out_cols=cols, dtype=np.float32)

    return read_chunk_fn, read_rows_fn, read_chunk_cols_fn


def _assign_hardware_split_halves_r17(
    nt_mask: np.ndarray,
    flasks: Optional[Sequence[str]],
    welltags: Optional[Sequence[str]],
    seed: int = DEFAULT_R17_SEED,
    min_batch_nt: int = DEFAULT_R17_MIN_BATCH_NT,
) -> Tuple[np.ndarray, np.ndarray, str, List[str], List[str], int]:
    """
    Assign `(batches, half_flags, half_rule_used, half_a_units, half_b_units, n_fallback_cells)`:
      1. Primary unit: `flask` if `>= 2` distinct flasks each have `>= min_batch_nt` NTCs;
      2. Secondary unit: `welltag` if `>= 2` distinct welltags each have `>= min_batch_nt` NTCs;
      3. Fallback: deterministic stratified 50/50 split (`seed = 20261701`, explicitly flagged).
    Within hardware-unit splits, lexicographical 1-based odd -> Half A (`1`), even -> Half B (`2`),
    followed by within-half `< min_batch_nt` fallback (`_apply_min_ntc_batch_fallback`).
    """
    n_cells = int(nt_mask.size)
    for col_label, raw_units in (("flask", flasks), ("welltag", welltags)):
        if raw_units is None or len(raw_units) != n_cells:
            continue
        u_arr = np.asarray([str(x).strip() for x in raw_units], dtype=object)
        valid_u = sorted({u for u in u_arr if u and u.lower() not in _UNASSIGNED_TOKENS})
        if len(valid_u) < 2:
            continue
        # Count NTCs per unit
        nt_counts = {u: int(np.sum(nt_mask & (u_arr == u))) for u in valid_u}
        n_qual_u = sum(1 for c in nt_counts.values() if c >= min_batch_nt)
        if n_qual_u >= 2:
            # 1-based odd -> Half A (1), even -> Half B (2)
            half_a_units = [u for idx, u in enumerate(valid_u, start=1) if idx % 2 == 1]
            half_b_units = [u for idx, u in enumerate(valid_u, start=1) if idx % 2 == 0]
            a_set = set(half_a_units)
            b_set = set(half_b_units)
            half_flags = np.zeros(n_cells, dtype=np.int8)
            batches = np.full(n_cells, "", dtype=object)
            for i in range(n_cells):
                u = str(u_arr[i])
                if u in a_set:
                    half_flags[i] = 1
                    batches[i] = f"{col_label}:{u}"
                elif u in b_set:
                    half_flags[i] = 2
                    batches[i] = f"{col_label}:{u}"
                else:
                    half_flags[i] = 0
                    batches[i] = ""
            n_ha_nt = int(np.sum(nt_mask & (half_flags == 1)))
            n_hb_nt = int(np.sum(nt_mask & (half_flags == 2)))
            if n_ha_nt >= min_batch_nt and n_hb_nt >= min_batch_nt:
                batches_fb, n_fb = _apply_min_ntc_batch_fallback(
                    batches.astype(str), half_flags, nt_mask, min_batch_nt=min_batch_nt
                )
                return (
                    batches_fb,
                    half_flags,
                    f"{col_label}_lexicographic_odd_even",
                    half_a_units,
                    half_b_units,
                    int(n_fb),
                )

    # Fallback: deterministic 50/50 split stratified by nt_mask (`seed = 20261701`)
    rng = np.random.default_rng(int(seed))
    half_flags = np.zeros(n_cells, dtype=np.int8)
    for is_nt in (True, False):
        idx = np.flatnonzero(nt_mask == is_nt)
        perm = rng.permutation(idx)
        half_flags[perm[0::2]] = 1
        half_flags[perm[1::2]] = 2
    batches = np.where(half_flags == 1, "half_1", "half_2").astype(str)
    return batches, half_flags, f"seed{seed}_stratified_50_50_fallback", ["half_1"], ["half_2"], 0


def build_pfizer_telohaec_parts_r17(
    raw_path: Path,
    hgnc: HGNCLookup,
    condition: str = "untreated",
    x_layer: str = "auto",
    min_umi: float = DEFAULT_R17_MIN_UMI,
    min_ngenes: int = DEFAULT_R17_MIN_NGENES,
    min_batch_nt: int = DEFAULT_R17_MIN_BATCH_NT,
    seed: int = DEFAULT_R17_SEED,
    chunk_rows: int = 4096,
) -> Tuple[List[PartSpec], Dict[str, Any]]:
    """
    Build `PartSpec` shard(s) and counts-only reader audit for Pfizer TeloHAEC (`imru_full.h5ad`)
    strictly filtered to `condition` (`'untreated'` primary basal OR `'treated'` IL-1b+TNF-a secondary).
    """
    cond_req = str(condition).strip().lower()
    if cond_req not in ("untreated", "treated"):
        raise ValueError(f"Unsupported Pfizer TeloHAEC condition: {condition!r} (expected 'untreated' or 'treated')")

    h5_path = _resolve_h5ad_path_r17(Path(raw_path))
    with h5py.File(h5_path, "r") as f:
        chosen_slot, int_slots, slot_diag = _resolve_h5ad_integer_slot(
            f, x_layer=x_layer, file_label=h5_path.name
        )
        obs_grp = f["obs"]
        obs_keys = sorted(k for k in obs_grp.keys() if k != "__categories")
        obs_idx_col = _resolve_h5ad_index_col(obs_grp, group_name="obs")
        obs_index = [x.strip() for x in read_h5ad_table_col_r17(obs_grp, obs_idx_col)]
        n_obs = len(obs_index)

        var_grp = _get_var_group_for_slot(f, chosen_slot)
        node = _get_h5ad_slot_node(f, chosen_slot)
        (mat_r, mat_c), layout = _inspect_slot_shape_and_layout(node, n_obs, len(read_h5ad_table_col_r17(var_grp, _resolve_h5ad_index_col(var_grp, group_name="var"))))
        var_meta = _extract_var_rna_features_r17(var_grp, n_var_expected=mat_c)

        c_cond = _pick_obs_col_r17(obs_keys, ["condition", "cell_treat", "treatment", "stim", "stimulation"])
        c_tgt = _pick_obs_col_r17(obs_keys, ["gene_target", "target_gene", "perturbation", "gene", "target"])
        c_guide = _pick_obs_col_r17(obs_keys, ["guide", "sgRNA", "sgrna", "guide_id", "sgrna_id"])
        c_cellid = _pick_obs_col_r17(obs_keys, ["cellID", "cell_id", "barcode"])
        c_flask = _pick_obs_col_r17(obs_keys, ["flask", "replicate", "batch", "gem_group", "channel", "lane"])
        c_well = _pick_obs_col_r17(obs_keys, ["welltag", "well", "plate_well"])
        c_nguide = _pick_obs_col_r17(obs_keys, ["num_features", "n_guides", "guide_count", "n_sgrna"])

        if c_cond is None or c_tgt is None:
            raise KeyError(
                f"Pfizer TeloHAEC .obs missing required condition/gene_target columns "
                f"(c_cond={c_cond}, c_tgt={c_tgt}, obs_keys={obs_keys})"
            )

        cond_vals = read_h5ad_table_col_r17(obs_grp, c_cond)
        tgt_vals = read_h5ad_table_col_r17(obs_grp, c_tgt)
        guide_vals = read_h5ad_table_col_r17(obs_grp, c_guide) if c_guide else [""] * n_obs
        cellid_vals = read_h5ad_table_col_r17(obs_grp, c_cellid) if c_cellid else obs_index
        flask_vals = read_h5ad_table_col_r17(obs_grp, c_flask) if c_flask else None
        well_vals = read_h5ad_table_col_r17(obs_grp, c_well) if c_well else None
        nguide_vals = read_h5ad_table_col_r17(obs_grp, c_nguide) if c_nguide else None

    # First pass over obs metadata: filter condition + guide/target rules
    cond_counts = {"untreated": 0, "treated": 0, "other": 0}
    cand_rows: List[int] = []
    cand_is_nt: List[bool] = []
    cand_pkeys: List[str] = []
    cand_barcodes: List[str] = []
    cand_flasks: List[str] = []
    cand_wells: List[str] = []
    pert_meta: Dict[str, Tuple[str, str, str, str]] = {}

    n_safe_excl = 0
    n_multiguide_or_unassigned_excl = 0
    n_unapproved_excl = 0

    for i in range(n_obs):
        cnorm = classify_pfizer_condition_r17(cond_vals[i])
        cond_counts[cnorm] += 1
        if cnorm != cond_req:
            continue

        if nguide_vals is not None:
            try:
                ng_i = float(nguide_vals[i])
                if np.isfinite(ng_i) and int(round(ng_i)) != 1:
                    n_multiguide_or_unassigned_excl += 1
                    continue
            except ValueError:
                pass

        cls, pinfo, gkey = classify_pfizer_guide_r17(tgt_vals[i], guide_vals[i], hgnc)
        if cls == "safe":
            n_safe_excl += 1
            continue
        if cls == "exclude":
            gt_s = str(tgt_vals[i]).strip().lower()
            if gt_s in _UNASSIGNED_TOKENS or _MULTI_GUIDE_RE.search(str(tgt_vals[i])) or _MULTI_GUIDE_RE.search(str(guide_vals[i])):
                n_multiguide_or_unassigned_excl += 1
            else:
                n_unapproved_excl += 1
            continue

        cand_rows.append(i)
        if cls == "ntc":
            cand_is_nt.append(True)
            cand_pkeys.append("")
        else:
            assert pinfo is not None
            cand_is_nt.append(False)
            cand_pkeys.append(gkey)
            pert_meta.setdefault(gkey, pinfo)
        bc_base = str(cellid_vals[i]).strip() or obs_index[i]
        cand_barcodes.append(f"{cond_req}:{bc_base}")
        cand_flasks.append(flask_vals[i] if flask_vals is not None else "")
        cand_wells.append(well_vals[i] if well_vals is not None else "")

    cand_rows_arr = np.asarray(cand_rows, dtype=np.int64)
    umi_cand, ngenes_cand = _compute_rna_umi_and_ngenes_r17(
        h5_path, chosen_slot, cand_rows_arr, var_meta["rna_cols"], chunk_rows=chunk_rows
    )

    qc_ok = (umi_cand >= float(min_umi)) & (ngenes_cand >= int(min_ngenes))
    n_rna_qc_excl = int(np.sum(~qc_ok))
    keep_local = np.flatnonzero(qc_ok)
    if keep_local.size == 0:
        raise RuntimeError(f"Pfizer TeloHAEC ({cond_req}): zero cells passed RNA QC (umi>={min_umi}, ngenes>={min_ngenes}).")

    keep_rows = cand_rows_arr[keep_local]
    nt_mask = np.asarray([cand_is_nt[int(k)] for k in keep_local], dtype=bool)
    umi_kept = umi_cand[keep_local]
    pkeys_kept = [cand_pkeys[int(k)] for k in keep_local]
    barcodes_kept = [cand_barcodes[int(k)] for k in keep_local]
    flasks_kept = [cand_flasks[int(k)] for k in keep_local] if flask_vals is not None else None
    wells_kept = [cand_wells[int(k)] for k in keep_local] if well_vals is not None else None

    batches, half_flags, half_rule, ha_units, hb_units, n_fb_cells = _assign_hardware_split_halves_r17(
        nt_mask=nt_mask,
        flasks=flasks_kept,
        welltags=wells_kept,
        seed=seed,
        min_batch_nt=min_batch_nt,
    )

    read_chunk_fn, read_rows_fn, read_chunk_cols_fn = _make_h5ad_subset_readers_r17(
        h5_path=h5_path,
        slot=chosen_slot,
        keep_rows=keep_rows,
        rna_cols=var_meta["rna_cols"],
        n_obs_total=n_obs,
        n_var_total=mat_c,
        layout=layout,
    )

    part = PartSpec(
        tag=f"Pfizer_TeloHAEC_{cond_req}",
        path=h5_path,
        n_cells=int(keep_rows.size),
        n_genes_raw=int(var_meta["rna_cols"].size),
        nt_mask=nt_mask,
        batches=batches,
        half_flags=half_flags,
        umi=umi_kept,
        pert_keys=pkeys_kept,
        pert_meta=pert_meta,
        var_gene_ids=var_meta["rna_ids"],
        var_gene_names=var_meta["rna_syms"],
        read_chunk_fn=read_chunk_fn,
        read_rows_fn=read_rows_fn,
        read_chunk_cols_fn=read_chunk_cols_fn,
    )

    n_ha_ntc = int(np.sum(nt_mask & (half_flags == 1)))
    n_hb_ntc = int(np.sum(nt_mask & (half_flags == 2)))
    audit = {
        "source_zenodo_id": PFIZER_TELOHAEC_ZENODO_ID,
        "source_file": h5_path.name,
        "condition_selected": cond_req,
        "condition_cell_counts_raw": cond_counts,
        "x_slot_used": chosen_slot,
        "layout": layout,
        "int_slots_found": int_slots,
        "slot_selection_mode": slot_diag["selection_mode"],
        "obs_columns_resolved": {
            "condition": c_cond,
            "gene_target": c_tgt,
            "guide": c_guide,
            "cellID": c_cellid,
            "flask": c_flask,
            "welltag": c_well,
            "n_guides": c_nguide,
        },
        "n_raw_obs_total": int(n_obs),
        "n_var_total": int(mat_c),
        "n_var_rna_kept": int(var_meta["rna_cols"].size),
        "n_var_non_rna_excluded": int(var_meta["n_non_rna_excluded"]),
        "n_safe_targeting_excluded": int(n_safe_excl),
        "n_multiguide_or_unassigned_excluded": int(n_multiguide_or_unassigned_excl),
        "n_unapproved_symbol_excluded": int(n_unapproved_excl),
        "n_rna_qc_excluded": int(n_rna_qc_excl),
        "min_umi": float(min_umi),
        "min_ngenes": int(min_ngenes),
        "max_pct_mt": None,
        "n_qc_kept_cells": int(keep_rows.size),
        "n_kept_ntc": int(np.sum(nt_mask)),
        "n_kept_targeting_singlets": int(np.sum(~nt_mask)),
        "n_kept_half_a_ntc": n_ha_ntc,
        "n_kept_half_b_ntc": n_hb_ntc,
        "per_half_ntc_ge15": bool(n_ha_ntc >= min_batch_nt and n_hb_ntc >= min_batch_nt),
        "half_rule_used": half_rule,
        "half_a_units": ha_units,
        "half_b_units": hb_units,
        "n_low_ntc_batch_fallback_cells": int(n_fb_cells),
        "n_unique_pert_keys": len(pert_meta),
        "n_unique_canon_target_genes": len({v[1] for v in pert_meta.values()}),
        "kept_barcodes_sha256": hashlib.sha256("\n".join(barcodes_kept).encode("utf-8")).hexdigest(),
    }
    return [part], audit


def extract_pfizer_telohaec_ntc_raw_counts_r17(
    raw_path: Path,
    hgnc: HGNCLookup,
    condition: str = "untreated",
    x_layer: str = "auto",
    min_umi: float = DEFAULT_R17_MIN_UMI,
    min_ngenes: int = DEFAULT_R17_MIN_NGENES,
    seed: int = DEFAULT_R17_SEED,
) -> Dict[str, Any]:
    """
    Model-free extraction of QC-passing `NO-TARGET` NTC raw integer UMI counts, canonical
    barcodes, and batch labels for deterministic control sampling (`PREREG_J §1.2`).
    """
    parts, audit = build_pfizer_telohaec_parts_r17(
        raw_path=raw_path,
        hgnc=hgnc,
        condition=condition,
        x_layer=x_layer,
        min_umi=min_umi,
        min_ngenes=min_ngenes,
        seed=seed,
    )
    part = parts[0]
    nt_idx = np.flatnonzero(part.nt_mask)
    raw_ntc = part.read_rows_fn(nt_idx).astype(np.float32)

    # Reconstruct canonical barcodes identically for NTC rows
    h5_path = _resolve_h5ad_path_r17(Path(raw_path))
    cond_req = str(condition).strip().lower()
    with h5py.File(h5_path, "r") as f:
        obs_grp = f["obs"]
        obs_idx_col = _resolve_h5ad_index_col(obs_grp, group_name="obs")
        obs_index = [x.strip() for x in read_h5ad_table_col_r17(obs_grp, obs_idx_col)]
        c_cond = audit["obs_columns_resolved"]["condition"]
        c_tgt = audit["obs_columns_resolved"]["gene_target"]
        c_guide = audit["obs_columns_resolved"]["guide"]
        c_cellid = audit["obs_columns_resolved"]["cellID"]
        c_nguide = audit["obs_columns_resolved"]["n_guides"]
        cond_vals = read_h5ad_table_col_r17(obs_grp, c_cond)
        tgt_vals = read_h5ad_table_col_r17(obs_grp, c_tgt)
        guide_vals = read_h5ad_table_col_r17(obs_grp, c_guide) if c_guide else [""] * len(obs_index)
        cellid_vals = read_h5ad_table_col_r17(obs_grp, c_cellid) if c_cellid else obs_index
        nguide_vals = read_h5ad_table_col_r17(obs_grp, c_nguide) if c_nguide else None

    cand_rows: List[int] = []
    cand_bcs: List[str] = []
    for i in range(len(obs_index)):
        if classify_pfizer_condition_r17(cond_vals[i]) != cond_req:
            continue
        if nguide_vals is not None:
            try:
                ng_i = float(nguide_vals[i])
                if np.isfinite(ng_i) and int(round(ng_i)) != 1:
                    continue
            except ValueError:
                pass
        cls, _, _ = classify_pfizer_guide_r17(tgt_vals[i], guide_vals[i], hgnc)
        if cls in ("ntc", "target"):
            cand_rows.append(i)
            bc_base = str(cellid_vals[i]).strip() or obs_index[i]
            cand_bcs.append(f"{cond_req}:{bc_base}")

    with h5py.File(h5_path, "r") as f:
        var_grp = _get_var_group_for_slot(f, audit["x_slot_used"])
        var_meta = _extract_var_rna_features_r17(var_grp, n_var_expected=audit["n_var_total"])
    umi_c, ng_c = _compute_rna_umi_and_ngenes_r17(
        h5_path, audit["x_slot_used"], np.asarray(cand_rows, dtype=np.int64), var_meta["rna_cols"]
    )
    keep_local = np.flatnonzero((umi_c >= float(min_umi)) & (ng_c >= int(min_ngenes)))
    kept_bcs = [cand_bcs[int(k)] for k in keep_local]
    ntc_barcodes = [kept_bcs[int(i)] for i in nt_idx]
    ntc_batches = [str(part.batches[int(i)]) for i in nt_idx]

    return {
        "raw_counts": raw_ntc,
        "gene_names": list(part.var_gene_names),
        "barcodes": ntc_barcodes,
        "batch_labels": ntc_batches,
        "audit": audit,
    }


def prepare_pfizer_telohaec_r17(
    raw_path: Path,
    condition: str = "untreated",
    out_root: Path = REPO_ROOT / "round17" / "shared" / "lines",
    qc_md_path: Optional[Path] = None,
    var_dims_pkl: Path = R16_VAR_DIMS_PKL,
    workers: int = 8,
    chunk_size: int = 8192,
) -> Dict[str, Any]:
    """
    Post-human-approval standardization & COUNTS-ONLY QC gate for Pfizer TeloHAEC (`Zenodo 18792213`).
    Zero model inference and zero correlation metrics are computed.
    """
    cond_req = str(condition).strip().lower()
    line_name = "pfizer_telohaec_basal" if cond_req == "untreated" else "pfizer_telohaec_cytokine"
    hgnc_l = build_hgnc_lookup()
    parts, audit = build_pfizer_telohaec_parts_r17(raw_path=raw_path, hgnc=hgnc_l, condition=cond_req)
    line_dir = process_r10_line(
        line_name=line_name,
        parts=parts,
        hgnc=hgnc_l,
        chunk_size=chunk_size,
        out_root=out_root,
        extra_meta={"adapter_audit": audit},
        n_workers=workers,
    )

    bridge = HGNCBridge.from_tsv(DEFAULT_HGNC_FILE)
    common_genes = load_r9_common_genes(hgnc_bridge=bridge, common_genes_file=R9_COMMON_GENES_FILE)
    ds_t = load_standardized_line_r10(
        line_dir, representation="log1p", common_genes=common_genes, hgnc_bridge=bridge, min_gene_overlap=1
    )
    n_g_r9 = len(ds_t.genes)

    src_targets: Set[str] = set()
    for sn in R9_SOURCE_POOL_6:
        s_ds = load_standardized_line_r10(
            R9_LINES_DIR / sn, representation="log1p", common_genes=common_genes, hgnc_bridge=bridge
        )
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
    with open(var_dims_pkl, "rb") as fh:
        vd = pickle.load(fh)
    g_model_canon = {bridge.canonicalize(str(g)) or str(g) for g in vd["gene_names"]}
    p_model_canon = {bridge.canonicalize(str(p)) or str(p) for p in vd["pert_names"] if str(p) != "non-targeting"}

    n_g_state = sum(1 for g in ds_t.genes if g in g_model_canon)
    n_eval_r17 = sum(1 for t in p5_targets if t in p_model_canon)
    passes = bool(n_eval_r17 >= 50 and n_g_r9 >= 3000 and audit["per_half_ntc_ge15"])
    role_out = (
        ("PASS_CONFIRMATORY_PRIMARY" if cond_req == "untreated" else "SECONDARY_WITHIN_FAMILY_REPLICATION")
        if passes
        else "SECONDARY_UNDERPOWERED_OR_GATE_FAIL"
    )

    qc_summary = {
        "line": line_name,
        "condition": cond_req,
        "role": role_out,
        "passes_prereg_i_floors": passes,
        "n_eval_p5_r9": len(p5_targets),
        "n_eval_r17": int(n_eval_r17),
        "n_genes_g_eval_r9": int(n_g_r9),
        "n_genes_g_eval_state": int(n_g_state),
        "n_kept_half_a_ntc": int(audit["n_kept_half_a_ntc"]),
        "n_kept_half_b_ntc": int(audit["n_kept_half_b_ntc"]),
        "per_half_ntc_ge15": bool(audit["per_half_ntc_ge15"]),
        "adapter_audit": audit,
    }
    if qc_md_path is not None:
        qc_md_path.parent.mkdir(parents=True, exist_ok=True)
        qc_md_path.write_text(
            f"# Round 17 Pfizer TeloHAEC (`{line_name}`, Zenodo `{PFIZER_TELOHAEC_ZENODO_ID}`) Counts-Only QC Gate\n\n"
            f"- **Role:** `{role_out}`\n"
            f"- **`N_eval^R17` (`r9_6 ∩ P_model`):** `{n_eval_r17}` (`>= 50` required)\n"
            f"- **`|G_eval^R9|`:** `{n_g_r9}` (`>= 3,000` required); **`|G_eval^STATE|`:** `{n_g_state}`\n"
            f"- **Per-Half NTC (`NO-TARGET` only; `SAFE_TARGET` excluded):** Half A = `{audit['n_kept_half_a_ntc']}`, Half B = `{audit['n_kept_half_b_ntc']}` (`>= 15` required)\n\n"
            "```json\n" + json.dumps(qc_summary, indent=2) + "\n```\n",
            encoding="utf-8",
        )
    return qc_summary
