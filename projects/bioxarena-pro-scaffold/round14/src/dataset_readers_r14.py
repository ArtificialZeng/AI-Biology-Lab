#!/usr/bin/env python3
"""
Round 14 Dataset Ingestion Adapter for `GSE337804` (`PREREG_G` / `a549_myllia`).

Reads the two `GSE337804_A549_ch_*.RDS.h5ad` files (`ch_1`, `ch_2`) and optional
`GSE337804_sgRNA_lib.csv.gz` cleanly, reusing frozen Round 10 / 11 / 12 primitives
without modifying any file under `round10/`, `round11/`, `round12/`, or `round13/`:

  1. Strict Pre-Metric Feasibility & Integer Count Slot Detection (`EXCLUDED_FEASIBILITY`):
     - Automatically scans `X`, `layers/*`, and `raw/X` (plus `raw/var` when `raw/X` is
       selected) in each `.h5ad` file to locate raw non-negative integer UMI counts; logs the
       detected slot per channel before computing any downstream metric. Fails cleanly with
       `EXCLUDED_FEASIBILITY` if no integer count slot exists.
  2. Single-sgRNA Calling & NTC Separation (`PREREG_G` Item 4):
     - Primary path: uses author per-cell `.obs` guide annotations (`sgrna_id`, `sgrna_symbol`,
       optional multiplicity column), strictly keeping 1-guide cells (excluding `[,;|+]`
       multi-guide / unassigned / multiplet cells).
     - Deterministic fallback (`PREREG_F §3.2`): if `.obs` lacks guide calls and `.var`
       contains non-RNA guide count features, fits per-guide Poisson-Gaussian EM thresholds
       on cells with `ngenes >= 200` within `a549_myllia` (`n_pos == 1`, `top1 >= 5`,
       `top2 < 3 or top1 / max(top2, 0.5) >= 3`).
     - Strict NTC vs Safe-Targeting rule: ONLY the `20` `ctrl_NT` (`non-targeting_1..20`)
       guides count as NTC; all `20` `ctrl_T` (`gene_desert_1..20` safe-targeting) cells
       are excluded strictly from both NTC and targets.
     - Target symbols canonicalized strictly via frozen `HGNCLookup.resolve`.
  3. RNA QC (`PREREG_G` Item 4):
     - Retain cells with `umi >= 300` AND `ngenes >= 200` over human RNA genes, with NO
       mitochondrial percentage upper bound (`max_pct_mt = None`, `DepleteX`-safe).
  4. Half Split & Unit Centering (`PREREG_G` Item 5):
     - `ch_1 -> Half A (1)` (`batch = "ch_1"`), `ch_2 -> Half B (2)` (`batch = "ch_2"`),
       with channel-level NTC subtraction (`< 15` NTC falls back to `half_{h}_pooled`).
  5. Integration: passes `PartSpec` shards into frozen `process_r10_line`.
"""

from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROUND14_SRC_DIR = Path(__file__).resolve().parent
REPO_ROOT = ROUND14_SRC_DIR.parents[1]
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
R10_SCRIPTS_DIR = REPO_ROOT / "round10" / "shared" / "scripts"
R11_SRC_DIR = REPO_ROOT / "round11" / "src"
R12_SRC_DIR = REPO_ROOT / "round12" / "src"
for p in (str(R10_SRC_DIR), str(R10_SCRIPTS_DIR), str(R11_SRC_DIR), str(R12_SRC_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from dataset_readers_r11 import (  # noqa: E402
    _apply_min_ntc_batch_fallback,
    _is_plain_single_gene_symbol,
    _pick_column,
)
from dataset_readers_r12 import (  # noqa: E402
    _is_nonhuman_or_vector_var,
    _resolve_h5ad_index_col,
)
from prepare_r10_lines import (  # noqa: E402
    HGNCLookup,
    PartSpec,
    build_hgnc_lookup,
    extract_top2_csc,
    fit_poisson_gaussian_thresholds_from_hist,
    process_r10_line,
    read_h5ad_chunk_dense,
    read_obs_or_var_col,
)

_CH_RE = re.compile(r"(?:^|[_\-/])(?:ch|channel|batch)[_-]?0*(\d+)(?:[_\-.]|$)", re.IGNORECASE)
_MULTI_GUIDE_RE = re.compile(r"[,;|+]")
_UNASSIGNED_TOKENS = frozenset(
    {"", "na", "n/a", "<na>", "nan", "none", "null", "unassigned", "multiplet", "doublet", "ambiguous", "unknown"}
)
_EXPLICIT_NT_TOKENS = frozenset({"non-targeting", "nontargeting", "non_targeting", "non-target", "ntc", "ctrl_nt"})
_SAFE_TOKENS = frozenset({"safe-targeting", "safe_targeting", "safetargeting", "safe", "ctrl_t", "gene_desert"})
_EXCLUDED_LAYER_TOKENS = (
    "nascent",
    "unspliced",
    "spliced",
    "ambiguous",
    "velocity",
    "intron",
    "exon",
    "lognorm",
    "scale",
)


def _resolve_raw_dir_r14(raw_dir: Path, line: str = "a549_myllia") -> Path:
    base = Path(raw_dir)
    for cand in (base / "gse337804", base / line, base):
        if cand.exists():
            return cand
    return base / "gse337804"


def _extract_channel_num(path: Path) -> Optional[int]:
    for text in (path.stem, path.name, path.parent.name):
        m = _CH_RE.search(text)
        if m:
            return int(m.group(1))
    return None


def _discover_gse337804_inputs(search_dir: Path, line: str = "a549_myllia") -> Tuple[List[Tuple[int, str, Path]], Optional[Path]]:
    """Locate `ch_1`/`ch_2` `.h5ad` files and optional `GSE337804_sgRNA_lib.csv.gz`."""
    if line != "a549_myllia":
        raise ValueError(f"Unsupported Round 14 GSE337804 line: {line!r} (expected 'a549_myllia')")

    h5ads = sorted(
        p for p in search_dir.rglob("*.h5ad")
        if p.is_file() and not p.name.startswith("._") and not p.name.startswith(".")
    )
    if not h5ads:
        raise FileNotFoundError(f"No GSE337804 .h5ad files found under {search_dir}")

    channel_entries: List[Tuple[int, str, Path]] = []
    for p in h5ads:
        ch_num = _extract_channel_num(p)
        if ch_num is not None:
            channel_entries.append((ch_num, f"ch_{ch_num}", p))
    if len(channel_entries) != len(h5ads):
        channel_entries = [(i + 1, f"ch_{i + 1}", p) for i, p in enumerate(h5ads)]
    channel_entries.sort(key=lambda t: (t[0], str(t[2])))

    lib_cands: List[Path] = []
    for pat in ("*sgrna*lib*.csv*", "*sgRNA*lib*.csv*", "*sgrna*.csv*", "*guide*lib*.csv*", "*sgrna*.tsv*", "*sgrna*.txt*"):
        for p in sorted(search_dir.rglob(pat)):
            if p.is_file() and not p.name.startswith("._") and p not in lib_cands:
                lib_cands.append(p)
    return channel_entries, (lib_cands[0] if lib_cands else None)


def _read_tabular_auto_r14(path: Path) -> pd.DataFrame:
    opener = gzip.open if str(path).endswith(".gz") else open
    first_data = ""
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if s and not s.startswith("#"):
                first_data = line
                break
    sep = "\t" if "\t" in first_data else ","
    return pd.read_csv(path, sep=sep, comment="#", dtype=str, keep_default_na=False)


def _is_safe_targeting_r14(gid: str, sym: str, cls_str: str = "") -> bool:
    for val in (gid, sym, cls_str):
        v = str(val).strip().lower()
        if not v:
            continue
        if v in _SAFE_TOKENS or "gene_desert" in v or "safe-target" in v or "safe_target" in v or "safetarget" in v:
            return True
    return False


def _is_non_targeting_r14(gid: str, sym: str, cls_str: str = "") -> bool:
    gid_l = str(gid).strip().lower()
    sym_l = str(sym).strip().lower()
    cls_l = str(cls_str).strip().lower()
    if cls_l in _EXPLICIT_NT_TOKENS:
        return True
    if gid_l.startswith(("non-targeting", "nontargeting", "non_targeting", "non-target")):
        return True
    if gid_l in _EXPLICIT_NT_TOKENS or sym_l in _EXPLICIT_NT_TOKENS:
        return True
    return False


def load_gse337804_sgrna_lib(lib_path: Optional[Path], hgnc: HGNCLookup) -> Dict[str, Any]:
    """
    Parse optional `GSE337804_sgRNA_lib.csv.gz` (`id,symbol,rank,class,entrez,guide,count`).
    Separates strictly `ctrl_NT` (`non-targeting_1..20`) as NTC and `ctrl_T`
    (`gene_desert_1..20` safe-targeting) as excluded.
    """
    ntc_keys: Set[str] = set()
    safe_keys: Set[str] = set()
    guide_to_target: Dict[str, Tuple[str, str, str, str]] = {}

    if lib_path is None or not Path(lib_path).exists():
        return {
            "source": "none_obs_only",
            "ntc_keys": ntc_keys,
            "safe_keys": safe_keys,
            "guide_to_target": guide_to_target,
        }

    df = _read_tabular_auto_r14(Path(lib_path))
    c_id = _pick_column(df.columns, ["id", "guide_id", "sgrna_id", "sgRNA", "sgrna", "name"])
    c_sym = _pick_column(df.columns, ["symbol", "gene_symbol", "sgrna_symbol", "target_gene_name", "gene", "target"])
    c_cls = _pick_column(df.columns, ["class", "guide_class", "guide_type", "type", "role", "category"])
    if c_id is None and c_sym is None:
        raise KeyError(f"sgRNA library {lib_path} missing id/symbol columns: {list(df.columns)}")

    for i in range(len(df)):
        gid = str(df[c_id].iloc[i]).strip() if c_id else ""
        sym = str(df[c_sym].iloc[i]).strip() if c_sym else ""
        cls_s = str(df[c_cls].iloc[i]).strip() if c_cls else ""
        key = gid or sym
        if not key:
            continue
        if _is_safe_targeting_r14(gid, sym, cls_s):
            safe_keys.add(key)
            if sym and sym.lower() == "gene_desert":
                safe_keys.add(sym)
            continue
        if _is_non_targeting_r14(gid, sym, cls_s):
            ntc_keys.add(key)
            if sym and sym.lower() in _EXPLICIT_NT_TOKENS:
                ntc_keys.add(sym)
            continue
        if sym and _is_plain_single_gene_symbol(sym):
            canon_sym, hstat = hgnc.resolve(sym, "")
            if canon_sym in hgnc.approved:
                info = (hgnc.sym_to_ens.get(canon_sym, ""), canon_sym, sym, hstat)
                guide_to_target[key] = info

    return {
        "source": f"sgrna_lib:{Path(lib_path).name}",
        "ntc_keys": ntc_keys,
        "safe_keys": safe_keys,
        "guide_to_target": guide_to_target,
    }


def _is_excluded_layer_name(layer_name: str) -> bool:
    low = str(layer_name).strip().lower()
    return any(tok in low for tok in _EXCLUDED_LAYER_TOKENS)


def _normalize_slot_name(f: h5py.File, slot: str) -> str:
    s = str(slot).strip()
    if s in ("X", "raw/X") or s.startswith("layers/"):
        return s
    if s == "raw" and "raw" in f and isinstance(f["raw"], h5py.Group) and "X" in f["raw"] and (
        "layers" not in f or "raw" not in f["layers"]
    ):
        return "raw/X"
    return f"layers/{s}"


def _ordered_h5ad_candidate_slots(f: h5py.File) -> Tuple[List[str], List[str]]:
    """
    Return `(candidate_slots, excluded_layer_slots)`:
      - `candidate_slots` covers `X`, non-excluded `layers/*`, and `raw/X`.
      - Any layer containing subfraction/derived tokens (`nascent`, `unspliced`, `spliced`,
        `ambiguous`, `velocity`, `intron`, `exon`, `lognorm`, `scale`) is excluded.
    """
    candidates: List[str] = []
    excluded: List[str] = []
    if "X" in f:
        candidates.append("X")
    if "layers" in f:
        avail = sorted(f["layers"].keys())
        kept_lk = []
        for lk in avail:
            if _is_excluded_layer_name(lk):
                excluded.append(f"layers/{lk}")
            else:
                kept_lk.append(lk)
        preferred = [k for k in ("counts", "raw_counts", "raw", "UMI", "umi") if k in kept_lk]
        for lk in preferred + [k for k in kept_lk if k not in preferred]:
            candidates.append(f"layers/{lk}")
    if "raw" in f and isinstance(f["raw"], h5py.Group) and "X" in f["raw"]:
        candidates.append("raw/X")
    return candidates, excluded


def _get_h5ad_slot_node(f: h5py.File, slot: str) -> Any:
    s = _normalize_slot_name(f, slot)
    if s == "X":
        return f["X"]
    if s == "raw/X":
        return f["raw"]["X"]
    if s.startswith("layers/"):
        return f["layers"][s.split("/", 1)[1]]
    return f["layers"][s]


def _get_var_group_for_slot(f: h5py.File, slot: str) -> h5py.Group:
    s = _normalize_slot_name(f, slot)
    if s == "raw/X" and "raw" in f and isinstance(f["raw"], h5py.Group) and "var" in f["raw"]:
        return f["raw"]["var"]
    return f["var"]


def _inspect_slot_shape_and_layout(node: Any, n_obs: int, n_var: int) -> Tuple[Tuple[int, int], str]:
    """Return `((n_obs, n_var), layout)` where layout in `{'csr', 'csc', 'dense'}`."""
    if isinstance(node, h5py.Dataset):
        return (int(node.shape[0]), int(node.shape[1])), "dense"
    if isinstance(node, h5py.Group) and "indptr" in node:
        enc = str(node.attrs.get("encoding-type", "")).lower()
        shape_attr = node.attrs.get("shape")
        shape = (int(shape_attr[0]), int(shape_attr[1])) if shape_attr is not None else (n_obs, n_var)
        if "csc" in enc:
            return shape, "csc"
        if "csr" in enc:
            return shape, "csr"
        indptr_len = int(node["indptr"].shape[0])
        if indptr_len == shape[0] + 1 and indptr_len != shape[1] + 1:
            return shape, "csr"
        if indptr_len == shape[1] + 1 and indptr_len != shape[0] + 1:
            return shape, "csc"
        return shape, "csr"
    raise TypeError(f"Unsupported HDF5 matrix node for GSE337804: {type(node)}")


def _validate_slot_full_nonneg_int(
    f: h5py.File,
    slot: str,
    batch_rows: int = 4096,
) -> Dict[str, Any]:
    """
    Mechanically inspect EVERY value / ALL non-zero values of `slot` for finite non-negative
    exact integers and compute per-cell total UMI median (`PREREG_G §4.1`).
    """
    obs_grp = f["obs"]
    obs_idx_col = _resolve_h5ad_index_col(obs_grp, group_name="obs")
    n_obs = len(read_obs_or_var_col(obs_grp, obs_idx_col))
    var_grp = _get_var_group_for_slot(f, slot)
    var_idx_col = _resolve_h5ad_index_col(var_grp, group_name="var")
    n_var = len(read_obs_or_var_col(var_grp, var_idx_col))

    node = _get_h5ad_slot_node(f, slot)
    (mat_r, mat_c), layout = _inspect_slot_shape_and_layout(node, n_obs, n_var)
    if mat_r != n_obs or mat_c != n_var:
        raise ValueError(f"Shape mismatch in slot {slot!r}: matrix {(mat_r, mat_c)} vs obs/var {(n_obs, n_var)}")

    cell_umi = np.zeros(mat_r, dtype=np.float64)
    n_nonzero = 0
    n_negative = 0
    n_nonfinite = 0
    n_nonint = 0

    if layout == "csr":
        indptr_all = np.asarray(node["indptr"][:], dtype=np.int64)
        for cs in range(0, mat_r, batch_rows):
            ce = min(cs + batch_rows, mat_r)
            p0, p1 = int(indptr_all[cs]), int(indptr_all[ce])
            if p1 <= p0:
                continue
            dat = np.asarray(node["data"][p0:p1], dtype=np.float64)
            idx = np.asarray(node["indices"][p0:p1], dtype=np.int64)
            fin = np.isfinite(dat)
            n_nonfinite += int(np.sum(~fin))
            n_negative += int(np.sum(fin & (dat < 0.0)))
            nz = dat[fin & (dat != 0.0)]
            n_nonzero += int(nz.size)
            n_nonint += int(np.sum(np.abs(nz - np.round(nz)) > 1e-6))
            safe_dat = np.where(fin, dat, 0.0)
            blk = sp.csr_matrix((safe_dat, idx, indptr_all[cs : ce + 1] - p0), shape=(ce - cs, mat_c))
            cell_umi[cs:ce] = np.asarray(blk.sum(axis=1), dtype=np.float64).ravel()
    elif layout == "csc":
        indptr_all = np.asarray(node["indptr"][:], dtype=np.int64)
        for cs in range(0, mat_c, 2048):
            ce = min(cs + 2048, mat_c)
            p0, p1 = int(indptr_all[cs]), int(indptr_all[ce])
            if p1 <= p0:
                continue
            r_idx = np.asarray(node["indices"][p0:p1], dtype=np.int64)
            dat = np.asarray(node["data"][p0:p1], dtype=np.float64)
            fin = np.isfinite(dat)
            n_nonfinite += int(np.sum(~fin))
            n_negative += int(np.sum(fin & (dat < 0.0)))
            nz = dat[fin & (dat != 0.0)]
            n_nonzero += int(nz.size)
            n_nonint += int(np.sum(np.abs(nz - np.round(nz)) > 1e-6))
            safe_dat = np.where(fin, dat, 0.0)
            cell_umi += np.bincount(r_idx, weights=safe_dat, minlength=mat_r)
    else:
        for cs in range(0, mat_r, 2048):
            ce = min(cs + 2048, mat_r)
            blk = np.asarray(node[cs:ce, :], dtype=np.float64)
            fin = np.isfinite(blk)
            n_nonfinite += int(np.sum(~fin))
            n_negative += int(np.sum(fin & (blk < 0.0)))
            nz = blk[fin & (blk != 0.0)]
            n_nonzero += int(nz.size)
            n_nonint += int(np.sum(np.abs(nz - np.round(nz)) > 1e-6))
            cell_umi[cs:ce] = np.where(fin, blk, 0.0).sum(axis=1)

    ok_int = bool(n_nonzero > 0 and n_nonfinite == 0 and n_negative == 0 and n_nonint == 0)
    if n_nonfinite > 0 or n_negative > 0 or n_nonint > 0:
        status = f"non_integer_or_negative_{layout}"
    elif n_nonzero == 0:
        status = f"all_zero_{layout}"
    else:
        status = f"{layout}_nonneg_int"

    return {
        "eligible_candidate": True,
        "is_nonneg_int": ok_int,
        "status": status,
        "layout": layout,
        "shape": [int(mat_r), int(mat_c)],
        "n_nonzero": int(n_nonzero),
        "n_negative": int(n_negative),
        "n_nonfinite": int(n_nonfinite),
        "n_nonint": int(n_nonint),
        "median_total_umi": (float(np.median(cell_umi)) if cell_umi.size > 0 else 0.0),
    }


def _read_slot_rows_csr(
    node: Any,
    layout: str,
    cs: int,
    ce: int,
    n_var: int,
    indptr_all: Optional[np.ndarray] = None,
    csc_col_mat: Optional[sp.csr_matrix] = None,
) -> sp.csr_matrix:
    """Return rows `[cs:ce]` of a slot as `float64` CSR."""
    if layout == "csr":
        assert indptr_all is not None
        p0, p1 = int(indptr_all[cs]), int(indptr_all[ce])
        if p1 <= p0:
            return sp.csr_matrix((ce - cs, n_var), dtype=np.float64)
        dat = np.asarray(node["data"][p0:p1], dtype=np.float64)
        idx = np.asarray(node["indices"][p0:p1], dtype=np.int32)
        return sp.csr_matrix((dat, idx, indptr_all[cs : ce + 1] - p0), shape=(ce - cs, n_var))
    if layout == "csc":
        assert csc_col_mat is not None
        return csc_col_mat[cs:ce]
    return sp.csr_matrix(np.asarray(node[cs:ce, :], dtype=np.float64))


def _compare_slots_elementwise(
    f: h5py.File,
    slot_a: str,
    slot_b: str,
    diag_map: Dict[str, Dict[str, Any]],
    batch_rows: int = 4096,
) -> Dict[str, Any]:
    """
    Compare `slot_a` and `slot_b` element-wise after shape & gene alignment.
    """
    var_a = _get_var_group_for_slot(f, slot_a)
    var_b = _get_var_group_for_slot(f, slot_b)
    genes_a = [str(x).strip() for x in read_obs_or_var_col(var_a, _resolve_h5ad_index_col(var_a, group_name="var"))]
    genes_b = [str(x).strip() for x in read_obs_or_var_col(var_b, _resolve_h5ad_index_col(var_b, group_name="var"))]

    info_a = diag_map[slot_a]
    info_b = diag_map[slot_b]
    shape_a = tuple(info_a["shape"])
    shape_b = tuple(info_b["shape"])
    n_obs = shape_a[0]
    med_a = info_a["median_total_umi"]
    med_b = info_b["median_total_umi"]

    if shape_a[0] != shape_b[0]:
        return {
            "identical": False,
            "reason": f"row_count_mismatch({shape_a[0]} vs {shape_b[0]})",
            "shape_a": list(shape_a),
            "shape_b": list(shape_b),
            "median_total_umi_a": med_a,
            "median_total_umi_b": med_b,
        }

    col_perm: Optional[np.ndarray] = None
    if genes_a != genes_b:
        if (
            len(genes_a) == len(genes_b)
            and len(set(genes_a)) == len(genes_a)
            and len(set(genes_b)) == len(genes_b)
            and set(genes_a) == set(genes_b)
        ):
            pos_b = {g: j for j, g in enumerate(genes_b)}
            col_perm = np.asarray([pos_b[g] for g in genes_a], dtype=np.int64)
        else:
            only_a = len(set(genes_a) - set(genes_b))
            only_b = len(set(genes_b) - set(genes_a))
            return {
                "identical": False,
                "reason": (
                    f"shape_or_gene_mismatch(shape_a={shape_a}, shape_b={shape_b}, "
                    f"only_in_a={only_a}, only_in_b={only_b})"
                ),
                "shape_a": list(shape_a),
                "shape_b": list(shape_b),
                "median_total_umi_a": med_a,
                "median_total_umi_b": med_b,
            }

    node_a = _get_h5ad_slot_node(f, slot_a)
    node_b = _get_h5ad_slot_node(f, slot_b)
    lay_a = info_a["layout"]
    lay_b = info_b["layout"]
    ptr_a = np.asarray(node_a["indptr"][:], dtype=np.int64) if lay_a == "csr" else None
    ptr_b = np.asarray(node_b["indptr"][:], dtype=np.int64) if lay_b == "csr" else None
    csc_a = (
        sp.csc_matrix((node_a["data"][:], node_a["indices"][:], node_a["indptr"][:]), shape=shape_a).astype(np.float64).tocsr()
        if lay_a == "csc"
        else None
    )
    csc_b = (
        sp.csc_matrix((node_b["data"][:], node_b["indices"][:], node_b["indptr"][:]), shape=shape_b).astype(np.float64).tocsr()
        if lay_b == "csc"
        else None
    )

    n_diff_entries = 0
    max_abs_diff = 0.0
    for cs in range(0, n_obs, batch_rows):
        ce = min(cs + batch_rows, n_obs)
        blk_a = _read_slot_rows_csr(node_a, lay_a, cs, ce, shape_a[1], ptr_a, csc_a)
        blk_b = _read_slot_rows_csr(node_b, lay_b, cs, ce, shape_b[1], ptr_b, csc_b)
        if col_perm is not None:
            blk_b = blk_b[:, col_perm]
        diff = (blk_a - blk_b).tocsr()
        diff.eliminate_zeros()
        if diff.nnz > 0:
            n_diff_entries += int(diff.nnz)
            mval = float(np.max(np.abs(diff.data)))
            if mval > max_abs_diff:
                max_abs_diff = mval

    if n_diff_entries > 0 or max_abs_diff > 1e-6:
        return {
            "identical": False,
            "reason": f"value_mismatch(n_diff_entries={n_diff_entries}, max_abs_diff={max_abs_diff:.6g})",
            "shape_a": list(shape_a),
            "shape_b": list(shape_b),
            "n_diff_entries": int(n_diff_entries),
            "max_abs_diff": float(max_abs_diff),
            "median_total_umi_a": med_a,
            "median_total_umi_b": med_b,
        }
    return {
        "identical": True,
        "reason": "elementwise_identical",
        "shape_a": list(shape_a),
        "shape_b": list(shape_b),
        "n_diff_entries": 0,
        "max_abs_diff": 0.0,
        "median_total_umi_a": med_a,
        "median_total_umi_b": med_b,
    }


def _resolve_h5ad_integer_slot(
    f: h5py.File,
    x_layer: str = "auto",
    file_label: str = "",
) -> Tuple[str, List[str], Dict[str, Any]]:
    """
    Detect and validate candidate integer UMI slots under `PREREG_G §4.1`:
      - Excludes derived/subfraction layers matching `_EXCLUDED_LAYER_TOKENS`.
      - Validates ALL non-zero values of candidate & selected slots.
      - Under `auto`, if multiple candidate integer slots exist, compares them element-wise
        after shape & gene alignment: selects `X` iff identical; raises `RuntimeError`
        (`SLOT_AMBIGUOUS`) listing each slot's median total UMI and differences otherwise.
    """
    candidates, excluded_layers = _ordered_h5ad_candidate_slots(f)
    diag_map: Dict[str, Dict[str, Any]] = {}
    for ex_slot in excluded_layers:
        diag_map[ex_slot] = {
            "eligible_candidate": False,
            "is_nonneg_int": False,
            "status": "excluded_subfraction_or_derived_layer",
        }

    int_slots: List[str] = []
    for slot in candidates:
        res = _validate_slot_full_nonneg_int(f, slot)
        diag_map[slot] = res
        if res["is_nonneg_int"]:
            int_slots.append(slot)

    req = "auto" if x_layer is None else str(x_layer).strip()
    if req.lower() != "auto":
        norm_req = _normalize_slot_name(f, req)
        if norm_req not in diag_map or not diag_map[norm_req].get("eligible_candidate", False):
            # Explicitly validate caller-specified slot over all non-zero values
            diag_map[norm_req] = _validate_slot_full_nonneg_int(f, norm_req)
        if not diag_map[norm_req]["is_nonneg_int"]:
            raise RuntimeError(
                f"EXCLUDED_FEASIBILITY: GSE337804 slot {norm_req!r} failed full non-negative integer check "
                f"({diag_map[norm_req]['status']})."
            )
        return norm_req, int_slots, {
            "selection_mode": f"explicit:{norm_req}",
            "excluded_layers": excluded_layers,
            "slots": diag_map,
        }

    if not int_slots:
        raise RuntimeError(
            f"EXCLUDED_FEASIBILITY: GSE337804 {file_label or '.h5ad'} lacks verified non-negative integer UMI counts "
            f"across candidate slots {candidates} (diag={diag_map})."
        )

    if len(int_slots) == 1:
        chosen = int_slots[0]
        return chosen, int_slots, {
            "selection_mode": "unique_candidate",
            "excluded_layers": excluded_layers,
            "slots": diag_map,
        }

    ref_slot = "X" if "X" in int_slots else int_slots[0]
    comparisons: Dict[str, Dict[str, Any]] = {}
    mismatches: List[str] = []
    for oth in int_slots:
        if oth == ref_slot:
            continue
        cmp_res = _compare_slots_elementwise(f, ref_slot, oth, diag_map)
        comparisons[f"{ref_slot}_vs_{oth}"] = cmp_res
        if not cmp_res["identical"]:
            mismatches.append(f"{ref_slot} vs {oth}: {cmp_res['reason']}")

    if mismatches:
        slot_summaries = "; ".join(
            f"{s}(shape={diag_map[s]['shape']}, median_total_umi={diag_map[s]['median_total_umi']:.2f}, nnz={diag_map[s]['n_nonzero']})"
            for s in int_slots
        )
        prefix = f"{file_label}: " if file_label else ""
        raise RuntimeError(
            f"SLOT_AMBIGUOUS: {prefix}multiple non-identical integer UMI slots detected [{slot_summaries}]; "
            f"differences=[{'; '.join(mismatches)}]. Specify explicit `x_layer` after controller ruling."
        )

    return ref_slot, int_slots, {
        "selection_mode": "identical_multi_slot_prefer_X",
        "excluded_layers": excluded_layers,
        "comparisons": comparisons,
        "slots": diag_map,
    }


def _extract_var_metadata_r14(var_grp: h5py.Group, n_var_expected: Optional[int] = None) -> Dict[str, Any]:
    """Extract RNA vs CRISPR feature indices, gene IDs, symbols, and MT mask from `var_grp`."""
    var_keys = set(var_grp.keys())
    idx_col = _resolve_h5ad_index_col(var_grp, group_name="var")
    raw_idx = [str(x).strip() for x in read_obs_or_var_col(var_grp, idx_col)]
    n_var = len(raw_idx)
    if n_var_expected is not None and n_var != n_var_expected:
        raise ValueError(f"var[{idx_col!r}] length {n_var} != matrix cols {n_var_expected}")

    id_col = next((c for c in ("ensembl_id", "gene_ids", "gene_id", "feature_id", "ensembl") if c in var_keys), idx_col)
    name_col = next(
        (c for c in ("gene_symbol", "gene_symbols", "feature_name", "gene_name", "features", "symbol", "name") if c in var_keys),
        idx_col,
    )
    var_ids_all = [str(x).strip() for x in read_obs_or_var_col(var_grp, id_col)]
    var_names_all = [str(x).strip() for x in read_obs_or_var_col(var_grp, name_col)]

    ft_col = next((c for c in ("feature_types", "feature_type") if c in var_keys), None)
    if ft_col is not None:
        ftypes = [str(x).strip() for x in read_obs_or_var_col(var_grp, ft_col)]
        rna_raw = [i for i, t in enumerate(ftypes) if t.lower() == "gene expression"]
        crp_raw = [i for i, t in enumerate(ftypes) if "crispr" in t.lower() or "guide" in t.lower()]
        if not rna_raw and not crp_raw:
            rna_raw = list(range(n_var))
    else:
        ftypes = ["Gene Expression"] * n_var
        rna_raw = list(range(n_var))
        crp_raw = []

    rna_human = [i for i in rna_raw if not _is_nonhuman_or_vector_var(var_names_all[i], var_ids_all[i])]
    rna_rows = np.asarray(rna_human, dtype=np.int64)
    crp_rows = np.asarray(crp_raw, dtype=np.int64)
    r_ids = [var_ids_all[i].split(".")[0] or var_names_all[i] for i in rna_human]
    r_names = [var_names_all[i] for i in rna_human]
    c_ids = [var_ids_all[i] for i in crp_raw]
    c_names = [var_names_all[i] for i in crp_raw]
    return {
        "n_var": n_var,
        "var_idx_col": idx_col,
        "rna_rows": rna_rows,
        "crp_rows": crp_rows,
        "mt_mask_rna": np.array([gn.startswith("MT-") for gn in r_names], dtype=bool),
        "r_ids": r_ids,
        "r_names": r_names,
        "c_ids": c_ids,
        "c_names": c_names,
        "ftypes": ftypes,
    }


def _detect_obs_guide_columns(obs_grp: h5py.Group) -> Dict[str, Optional[str]]:
    obs_keys = list(obs_grp.keys())
    return {
        "guide_col": _pick_column(
            obs_keys,
            ["sgrna_id", "sgRNA_id", "sgRNA", "sgrna", "guide_id", "feature_call", "guide_call", "protospacer_call", "guide", "top_sgrna"],
        ),
        "sym_col": _pick_column(
            obs_keys,
            ["sgrna_symbol", "sgRNA_symbol", "perturbation", "target_gene", "target", "gene_symbol", "gene"],
        ),
        "class_col": _pick_column(
            obs_keys,
            ["sgrna_class", "guide_class", "guide_type", "class", "role", "guide_role"],
        ),
        "ncount_col": _pick_column(
            obs_keys,
            ["num_features", "n_features", "num_guides", "n_guides", "n_sgrnas", "n_perts", "nperts"],
        ),
    }


def _stream_h5ad_rna_qc_and_guides(
    x_node: Any,
    layout: str,
    n_obs: int,
    n_var: int,
    rna_rows: np.ndarray,
    mt_mask_rna: np.ndarray,
    crp_rows: Optional[np.ndarray] = None,
    min_ngenes_fit: int = 200,
    batch_size: int = 4096,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Optional[np.ndarray]]:
    """
    Stream per-cell `(umi_all, ng_all, pct_mt_all)` over RNA features strictly in bounded RAM.
    When `crp_rows` is non-empty, also compute `ng_all` in pass 1 and accumulate `fit_hist`
    over `ng_all >= min_ngenes_fit` per `PREREG_F §3.2`.
    """
    umi_all = np.zeros(n_obs, dtype=np.float64)
    ng_all = np.zeros(n_obs, dtype=np.int32)
    mt_umi_all = np.zeros(n_obs, dtype=np.float64)
    full_rna_subset = (len(rna_rows) == n_var and np.array_equal(rna_rows, np.arange(n_var)))
    mt_raw_idx = rna_rows[mt_mask_rna]

    if layout == "csr":
        indptr_all = np.asarray(x_node["indptr"][:], dtype=np.int64)
        for cs in range(0, n_obs, batch_size):
            ce = min(cs + batch_size, n_obs)
            p0, p1 = int(indptr_all[cs]), int(indptr_all[ce])
            if p1 <= p0:
                continue
            dat = np.asarray(x_node["data"][p0:p1], dtype=np.float64)
            idx = np.asarray(x_node["indices"][p0:p1], dtype=np.int64)
            rel_ptr = indptr_all[cs : ce + 1] - p0
            blk = sp.csr_matrix((dat, idx, rel_ptr), shape=(ce - cs, n_var))
            sub = blk if full_rna_subset else blk[:, rna_rows]
            umi_all[cs:ce] = np.asarray(sub.sum(axis=1), dtype=np.float64).ravel()
            ng_all[cs:ce] = np.asarray((sub > 0).sum(axis=1), dtype=np.int32).ravel()
            if len(mt_raw_idx) > 0:
                mt_umi_all[cs:ce] = np.asarray(blk[:, mt_raw_idx].sum(axis=1), dtype=np.float64).ravel()
    elif layout == "csc":
        is_rna_var = np.zeros(n_var, dtype=bool)
        is_rna_var[rna_rows] = True
        is_mt_var = np.zeros(n_var, dtype=bool)
        if len(mt_raw_idx) > 0:
            is_mt_var[mt_raw_idx] = True
        indptr_all = np.asarray(x_node["indptr"][:], dtype=np.int64)
        for cs in range(0, n_var, 2048):
            ce = min(cs + 2048, n_var)
            p0, p1 = int(indptr_all[cs]), int(indptr_all[ce])
            if p1 <= p0:
                continue
            r_idx = np.asarray(x_node["indices"][p0:p1], dtype=np.int64)
            dat = np.asarray(x_node["data"][p0:p1], dtype=np.float64)
            col_off = np.repeat(np.arange(cs, ce, dtype=np.int64), np.diff(indptr_all[cs : ce + 1]))
            m_rna = is_rna_var[col_off]
            if np.any(m_rna):
                rr = r_idx[m_rna]
                dd = dat[m_rna]
                umi_all += np.bincount(rr, weights=dd, minlength=n_obs)
                pos = dd > 0.0
                if np.any(pos):
                    ng_all += np.bincount(rr[pos], minlength=n_obs).astype(np.int32, copy=False)
            m_mt = is_mt_var[col_off]
            if np.any(m_mt):
                mt_umi_all += np.bincount(r_idx[m_mt], weights=dat[m_mt], minlength=n_obs)
    else:
        for cs in range(0, n_obs, batch_size):
            ce = min(cs + batch_size, n_obs)
            blk = np.asarray(x_node[cs:ce, :], dtype=np.float64)
            sub = blk if full_rna_subset else blk[:, rna_rows]
            umi_all[cs:ce] = sub.sum(axis=1)
            ng_all[cs:ce] = np.sum(sub > 0.0, axis=1).astype(np.int32, copy=False)
            if np.any(mt_mask_rna):
                mt_umi_all[cs:ce] = sub[:, mt_mask_rna].sum(axis=1)

    pct_mt_all = 100.0 * mt_umi_all / np.maximum(umi_all, 1.0)

    fit_hist: Optional[np.ndarray] = None
    if crp_rows is not None and len(crp_rows) > 0:
        n_crp = len(crp_rows)
        fit_hist = np.zeros((n_crp, 66), dtype=np.int64)
        ok_fit = ng_all >= int(min_ngenes_fit)
        if np.any(ok_fit):
            if layout == "csr":
                indptr_all = np.asarray(x_node["indptr"][:], dtype=np.int64)
                for cs in range(0, n_obs, batch_size):
                    ce = min(cs + batch_size, n_obs)
                    fit_local = np.where(ok_fit[cs:ce])[0]
                    if len(fit_local) == 0:
                        continue
                    p0, p1 = int(indptr_all[cs]), int(indptr_all[ce])
                    if p1 <= p0:
                        continue
                    dat = np.asarray(x_node["data"][p0:p1], dtype=np.float64)
                    idx = np.asarray(x_node["indices"][p0:p1], dtype=np.int64)
                    blk = sp.csr_matrix((dat, idx, indptr_all[cs : ce + 1] - p0), shape=(ce - cs, n_var))
                    c_fit = blk[fit_local][:, crp_rows].tocoo()
                    if c_fit.nnz > 0:
                        pos = c_fit.data > 0.0
                        if np.any(pos):
                            np.add.at(
                                fit_hist,
                                (c_fit.col[pos], np.clip(c_fit.data[pos].astype(np.int64), 1, 65)),
                                1,
                            )
            elif layout == "csc":
                indptr_all = np.asarray(x_node["indptr"][:], dtype=np.int64)
                for g_pos, raw_col in enumerate(crp_rows):
                    p0, p1 = int(indptr_all[int(raw_col)]), int(indptr_all[int(raw_col) + 1])
                    if p1 <= p0:
                        continue
                    r_idx = np.asarray(x_node["indices"][p0:p1], dtype=np.int64)
                    dat = np.asarray(x_node["data"][p0:p1], dtype=np.float64)
                    m = ok_fit[r_idx] & (dat > 0.0)
                    if np.any(m):
                        np.add.at(fit_hist[g_pos], np.clip(dat[m].astype(np.int64), 1, 65), 1)
            else:
                fit_cols = np.where(ok_fit)[0]
                for fs in range(0, len(fit_cols), batch_size):
                    f_sel = fit_cols[fs : fs + batch_size]
                    sub_c = np.asarray(x_node[f_sel][:, crp_rows], dtype=np.float64)
                    rr, cc = np.where(sub_c > 0.0)
                    if len(rr) > 0:
                        np.add.at(fit_hist, (cc, np.clip(sub_c[rr, cc].astype(np.int64), 1, 65)), 1)

    return umi_all, ng_all, pct_mt_all, fit_hist


def _extract_kept_h5ad_csr(
    x_node: Any,
    layout: str,
    keep_arr: np.ndarray,
    col_rows: np.ndarray,
    n_obs: int,
    n_var: int,
    batch_size: int = 4096,
) -> sp.csr_matrix:
    """Extract sorted `keep_arr` cells x `col_rows` features into compact `float32` CSR."""
    n_keep = len(keep_arr)
    n_cols_out = len(col_rows)
    if n_keep == 0 or n_cols_out == 0:
        return sp.csr_matrix((n_keep, n_cols_out), dtype=np.float32)

    full_cols = (n_cols_out == n_var and np.array_equal(col_rows, np.arange(n_var)))
    if layout == "csr":
        indptr_all = np.asarray(x_node["indptr"][:], dtype=np.int64)
        blocks: List[sp.csr_matrix] = []
        k_ptr = 0
        for cs in range(0, n_obs, batch_size):
            ce = min(cs + batch_size, n_obs)
            k_start = k_ptr
            while k_ptr < n_keep and keep_arr[k_ptr] < ce:
                k_ptr += 1
            if k_ptr == k_start:
                continue
            local_rows = keep_arr[k_start:k_ptr] - cs
            p0, p1 = int(indptr_all[cs]), int(indptr_all[ce])
            if p1 <= p0:
                blocks.append(sp.csr_matrix((len(local_rows), n_cols_out), dtype=np.float32))
                continue
            dat = np.asarray(x_node["data"][p0:p1], dtype=np.float32)
            idx = np.asarray(x_node["indices"][p0:p1], dtype=np.int32)
            blk = sp.csr_matrix((dat, idx, indptr_all[cs : ce + 1] - p0), shape=(ce - cs, n_var))
            sub = blk[local_rows] if full_cols else blk[local_rows][:, col_rows]
            blocks.append(sub.tocsr())
        return sp.vstack(blocks, format="csr") if len(blocks) > 1 else blocks[0]

    if layout == "csc":
        row_remap = np.full(n_obs, -1, dtype=np.int32)
        row_remap[keep_arr] = np.arange(n_keep, dtype=np.int32)
        indptr_all = np.asarray(x_node["indptr"][:], dtype=np.int64)
        r_l, c_l, d_l = [], [], []
        for out_j, raw_c in enumerate(col_rows):
            rc = int(raw_c)
            p0, p1 = int(indptr_all[rc]), int(indptr_all[rc + 1])
            if p1 <= p0:
                continue
            r_idx = np.asarray(x_node["indices"][p0:p1], dtype=np.int64)
            mapped = row_remap[r_idx]
            m = mapped >= 0
            if np.any(m):
                r_l.append(mapped[m])
                c_l.append(np.full(int(np.sum(m)), out_j, dtype=np.int32))
                d_l.append(np.asarray(x_node["data"][p0:p1], dtype=np.float32)[m])
        if not r_l:
            return sp.csr_matrix((n_keep, n_cols_out), dtype=np.float32)
        return sp.csr_matrix(
            (np.concatenate(d_l), (np.concatenate(r_l), np.concatenate(c_l))),
            shape=(n_keep, n_cols_out),
            dtype=np.float32,
        )

    blocks_d: List[sp.csr_matrix] = []
    for ks in range(0, n_keep, batch_size):
        k_sel = keep_arr[ks : ks + batch_size]
        blk = np.asarray(x_node[k_sel, :], dtype=np.float32)
        sub = blk if full_cols else blk[:, col_rows]
        blocks_d.append(sp.csr_matrix(sub))
    return sp.vstack(blocks_d, format="csr") if len(blocks_d) > 1 else blocks_d[0]


def _make_csr_part_callbacks(x_kept: sp.csr_matrix, n_genes_rna: int):
    def read_chunk(s: int, e: int, mat: sp.csr_matrix = x_kept) -> np.ndarray:
        return read_h5ad_chunk_dense(mat, s, e, n_genes_rna)

    def read_rows(rows: np.ndarray, mat: sp.csr_matrix = x_kept) -> np.ndarray:
        return mat[rows].toarray().astype(np.float64, copy=False)

    def read_chunk_cols(s: int, e: int, cols: np.ndarray, mat: sp.csr_matrix = x_kept) -> np.ndarray:
        return mat[s:e][:, cols].toarray().astype(np.float32, copy=False)

    return read_chunk, read_rows, read_chunk_cols


def check_gse337804_feasibility(
    raw_dir: Path,
    line: str = "a549_myllia",
    hgnc: Optional[HGNCLookup] = None,
    x_layer: str = "auto",
) -> Dict[str, Any]:
    """
    `PREREG_G` pre-metric feasibility gate (metadata & raw integer count structure ONLY;
    never computes any expression delta or predictor correlation).
    """
    try:
        hgnc_l = hgnc or build_hgnc_lookup()
        search_dir = _resolve_raw_dir_r14(raw_dir, line)
        ch_entries, lib_path = _discover_gse337804_inputs(search_dir, line)
        lib_meta = load_gse337804_sgrna_lib(lib_path, hgnc_l)

        per_channel_slots: Dict[str, Any] = {}
        n_rna_total = 0
        n_crp_total = 0
        n_cells_total = 0
        n_obs_ntc = n_obs_safe = n_obs_tgt = 0
        canon_tgts: Set[str] = set()
        used_obs_calls = False
        used_h5ad_crp_fallback = False

        for ch_num, ch_tag, hp in ch_entries:
            with h5py.File(hp, "r") as f:
                resolved_slot, int_slots, diag = _resolve_h5ad_integer_slot(f, x_layer=x_layer, file_label=hp.name)
                per_channel_slots[ch_tag] = {
                    "path": str(hp),
                    "x_slot_used": resolved_slot,
                    "int_count_sources": int_slots,
                    "slot_diagnostics": diag,
                    "selected_median_total_umi": float(diag["slots"][resolved_slot]["median_total_umi"]),
                }
                obs_grp = f["obs"]
                obs_idx_col = _resolve_h5ad_index_col(obs_grp, group_name="obs")
                n_obs = len(read_obs_or_var_col(obs_grp, obs_idx_col))
                var_grp = _get_var_group_for_slot(f, resolved_slot)
                vmeta = _extract_var_metadata_r14(var_grp)
                x_node = _get_h5ad_slot_node(f, resolved_slot)
                (mat_r, mat_c), _ = _inspect_slot_shape_and_layout(x_node, n_obs, vmeta["n_var"])
                if mat_r != n_obs or mat_c != vmeta["n_var"]:
                    raise ValueError(
                        f"Shape mismatch in {hp.name} ({resolved_slot}): matrix {(mat_r, mat_c)} vs obs/var {(n_obs, vmeta['n_var'])}"
                    )
                if len(vmeta["rna_rows"]) == 0:
                    return {"feasible": False, "reason": f"Zero human RNA features in {hp.name}"}

                n_cells_total += n_obs
                n_rna_total = len(vmeta["rna_rows"])
                n_crp_total = len(vmeta["crp_rows"])

                obs_cols = _detect_obs_guide_columns(obs_grp)
                has_obs = obs_cols["guide_col"] is not None or obs_cols["sym_col"] is not None
                if has_obs:
                    used_obs_calls = True
                    gvals = read_obs_or_var_col(obs_grp, obs_cols["guide_col"]) if obs_cols["guide_col"] else [""] * n_obs
                    svals = read_obs_or_var_col(obs_grp, obs_cols["sym_col"]) if obs_cols["sym_col"] else [""] * n_obs
                    cvals = read_obs_or_var_col(obs_grp, obs_cols["class_col"]) if obs_cols["class_col"] else [""] * n_obs
                    nvals = (
                        pd.to_numeric(pd.Series(read_obs_or_var_col(obs_grp, obs_cols["ncount_col"])), errors="coerce").to_numpy()
                        if obs_cols["ncount_col"]
                        else None
                    )
                    for i in range(n_obs):
                        if nvals is not None and np.isfinite(nvals[i]) and int(nvals[i]) != 1:
                            continue
                        gid = str(gvals[i]).strip()
                        sym = str(svals[i]).strip()
                        cls_s = str(cvals[i]).strip()
                        if _MULTI_GUIDE_RE.search(gid) or _MULTI_GUIDE_RE.search(sym):
                            continue
                        if gid.lower() in _UNASSIGNED_TOKENS and sym.lower() in _UNASSIGNED_TOKENS:
                            continue
                        if _is_safe_targeting_r14(gid, sym, cls_s) or gid in lib_meta["safe_keys"]:
                            n_obs_safe += 1
                            continue
                        if _is_non_targeting_r14(gid, sym, cls_s) or gid in lib_meta["ntc_keys"]:
                            n_obs_ntc += 1
                            continue
                        info = lib_meta["guide_to_target"].get(gid)
                        if info is not None:
                            n_obs_tgt += 1
                            canon_tgts.add(info[1])
                        elif sym and _is_plain_single_gene_symbol(sym):
                            cs, _ = hgnc_l.resolve(sym, "")
                            if cs in hgnc_l.approved:
                                n_obs_tgt += 1
                                canon_tgts.add(cs)
                elif len(vmeta["crp_rows"]) > 0:
                    used_h5ad_crp_fallback = True
                else:
                    return {
                        "feasible": False,
                        "reason": f"{hp.name} lacks both per-cell .obs guide calls and CRISPR guide count rows",
                    }

        ok_guides = (
            (used_obs_calls and n_obs_ntc > 0 and len(canon_tgts) > 0)
            or (used_h5ad_crp_fallback and n_crp_total > 0 and len(lib_meta["ntc_keys"]) > 0 and len(lib_meta["guide_to_target"]) > 0)
        )
        ok = bool(len(ch_entries) >= 2 and n_rna_total > 0 and ok_guides)
        return {
            "feasible": ok,
            "status": "OK" if ok else "EXCLUDED_FEASIBILITY",
            "reason": "OK" if ok else "Insufficient resolvable NTC/target guides or missing >=2 channels",
            "line": line,
            "n_h5ad_files": len(ch_entries),
            "channels": [t[1] for t in ch_entries],
            "channel_slots": per_channel_slots,
            "sgrna_lib_source": lib_meta["source"],
            "guide_call_mode": "obs_author_single_guide" if used_obs_calls else "h5ad_crispr_em_fallback",
            "n_raw_cells_total": int(n_cells_total),
            "n_rna_features": int(n_rna_total),
            "n_obs_ntc_meta": int(n_obs_ntc),
            "n_obs_safe_meta": int(n_obs_safe),
            "n_obs_tgt_meta": int(n_obs_tgt),
            "n_unique_canon_targets_meta": int(len(canon_tgts) if used_obs_calls else len({v[1] for v in lib_meta["guide_to_target"].values()})),
        }
    except Exception as exc:
        msg = f"{exc}" if isinstance(exc, RuntimeError) and str(exc).startswith(("SLOT_AMBIGUOUS:", "EXCLUDED_FEASIBILITY:")) else f"{type(exc).__name__}: {exc}"
        st = "SLOT_AMBIGUOUS" if "SLOT_AMBIGUOUS" in msg else "EXCLUDED_FEASIBILITY"
        return {"feasible": False, "status": st, "reason": msg}


def build_gse337804_parts(
    raw_dir: Path,
    hgnc: HGNCLookup,
    line: str = "a549_myllia",
    min_umi: int = 300,
    min_ngenes: int = 200,
    # PREREG_G Item 4: NO mitochondrial% limit (`DepleteX` removes MT/ribosomal transcripts).
    max_pct_mt: Optional[float] = None,
    min_batch_nt: int = 15,
    min_guide_thr: int = 5,
    default_guide_thr: int = 7,
    x_layer: str = "auto",
) -> Tuple[List[PartSpec], dict]:
    """
    Build standardized `PartSpec` shards for `GSE337804` (`a549_myllia`).
    """
    feas = check_gse337804_feasibility(raw_dir=raw_dir, line=line, hgnc=hgnc, x_layer=x_layer)
    if not feas.get("feasible", False):
        if feas.get("status") == "SLOT_AMBIGUOUS":
            raise RuntimeError(str(feas.get("reason")))
        raise RuntimeError(f"EXCLUDED_FEASIBILITY: {line} failed GSE337804 feasibility gate: {feas.get('reason')}")

    search_dir = _resolve_raw_dir_r14(raw_dir, line)
    ch_entries, lib_path = _discover_gse337804_inputs(search_dir, line)
    lib_meta = load_gse337804_sgrna_lib(lib_path, hgnc)
    ntc_keys: Set[str] = lib_meta["ntc_keys"]
    safe_keys: Set[str] = lib_meta["safe_keys"]
    guide_to_target: Dict[str, Tuple[str, str, str, str]] = lib_meta["guide_to_target"]

    ch_records: List[Dict[str, Any]] = []
    n_rna_qc_excl = 0
    using_fallback_em = False

    for ch_num, ch_tag, hp in ch_entries:
        with h5py.File(hp, "r") as f:
            resolved_slot, int_slots, slot_diag = _resolve_h5ad_integer_slot(f, x_layer=x_layer, file_label=hp.name)
            obs_grp = f["obs"]
            obs_idx_col = _resolve_h5ad_index_col(obs_grp, group_name="obs")
            barcodes = [str(x).strip() for x in read_obs_or_var_col(obs_grp, obs_idx_col)]
            n_obs = len(barcodes)

            var_grp = _get_var_group_for_slot(f, resolved_slot)
            vmeta = _extract_var_metadata_r14(var_grp)
            x_node = _get_h5ad_slot_node(f, resolved_slot)
            (mat_r, mat_c), layout = _inspect_slot_shape_and_layout(x_node, n_obs, vmeta["n_var"])

            obs_cols = _detect_obs_guide_columns(obs_grp)
            has_obs = obs_cols["guide_col"] is not None or obs_cols["sym_col"] is not None
            crp_rows = np.empty(0, dtype=np.int64) if has_obs else vmeta["crp_rows"]
            if not has_obs:
                using_fallback_em = True

            umi_all, ng_all, pct_mt_all, fit_hist = _stream_h5ad_rna_qc_and_guides(
                x_node=x_node,
                layout=layout,
                n_obs=n_obs,
                n_var=mat_c,
                rna_rows=vmeta["rna_rows"],
                mt_mask_rna=vmeta["mt_mask_rna"],
                crp_rows=crp_rows,
                min_ngenes_fit=min_ngenes,
            )

            ok_qc = (umi_all >= float(min_umi)) & (ng_all >= int(min_ngenes))
            if max_pct_mt is not None:
                ok_qc &= pct_mt_all <= float(max_pct_mt)
            n_rna_qc_excl += int(np.sum(~ok_qc))
            qc_cols = np.where(ok_qc)[0]

            g_vals = read_obs_or_var_col(obs_grp, obs_cols["guide_col"]) if obs_cols["guide_col"] else [""] * n_obs
            s_vals = read_obs_or_var_col(obs_grp, obs_cols["sym_col"]) if obs_cols["sym_col"] else [""] * n_obs
            c_vals = read_obs_or_var_col(obs_grp, obs_cols["class_col"]) if obs_cols["class_col"] else [""] * n_obs
            n_vals = (
                pd.to_numeric(pd.Series(read_obs_or_var_col(obs_grp, obs_cols["ncount_col"])), errors="coerce").to_numpy()
                if obs_cols["ncount_col"]
                else None
            )

            crp_qc_csc: Optional[sp.csc_matrix] = None
            crp_keys: List[str] = []
            if not has_obs and len(crp_rows) > 0:
                crp_qc_csr = _extract_kept_h5ad_csr(
                    x_node=x_node,
                    layout=layout,
                    keep_arr=qc_cols,
                    col_rows=crp_rows,
                    n_obs=n_obs,
                    n_var=mat_c,
                )
                crp_qc_csc = crp_qc_csr.T.tocsc()
                crp_keys = [
                    cid if (cid in ntc_keys or cid in safe_keys or cid in guide_to_target) else cname
                    for cid, cname in zip(vmeta["c_ids"], vmeta["c_names"])
                ]

        ch_records.append(
            {
                "path": hp,
                "ch_num": ch_num,
                "ch_tag": ch_tag,
                "slot": resolved_slot,
                "int_slots": int_slots,
                "slot_diag": slot_diag,
                "layout": layout,
                "n_obs": n_obs,
                "n_var": mat_c,
                "obs_idx_col": obs_idx_col,
                "rna_rows": vmeta["rna_rows"],
                "r_ids": vmeta["r_ids"],
                "r_names": vmeta["r_names"],
                "umi_all": umi_all,
                "qc_cols": qc_cols,
                "has_obs": has_obs,
                "g_vals": g_vals,
                "s_vals": s_vals,
                "c_vals": c_vals,
                "n_vals": n_vals,
                "fit_hist": fit_hist,
                "crp_qc_csc": crp_qc_csc,
                "crp_keys": crp_keys,
            }
        )

    thrs: Optional[np.ndarray] = None
    if using_fallback_em:
        hist_sum: Optional[np.ndarray] = None
        for rec in ch_records:
            fh = rec["fit_hist"]
            if fh is not None:
                hist_sum = fh.copy() if hist_sum is None else (hist_sum + fh)
        if hist_sum is None:
            raise RuntimeError("EXCLUDED_FEASIBILITY: guide EM fallback lacked CRISPR guide rows.")
        thrs = fit_poisson_gaussian_thresholds_from_hist(
            hist_sum,
            min_thr=min_guide_thr,
            max_thr=20,
            default_thr=default_guide_thr,
        )

    all_kept_parts: List[Dict[str, Any]] = []
    n_safe_excl = n_unmapped_excl = n_multi_or_zero_excl = 0

    for rec in ch_records:
        qc_cols = rec["qc_cols"]
        keep_rows: List[int] = []
        nt_flags: List[bool] = []
        pkeys: List[str] = []
        pmeta: Dict[str, Tuple[str, str, str, str]] = {}

        if rec["has_obs"]:
            g_vals, s_vals, c_vals, n_vals = rec["g_vals"], rec["s_vals"], rec["c_vals"], rec["n_vals"]
            for raw_r in qc_cols:
                r_i = int(raw_r)
                if n_vals is not None and np.isfinite(n_vals[r_i]) and int(n_vals[r_i]) != 1:
                    n_multi_or_zero_excl += 1
                    continue
                gid = str(g_vals[r_i]).strip()
                sym = str(s_vals[r_i]).strip()
                cls_s = str(c_vals[r_i]).strip()
                if _MULTI_GUIDE_RE.search(gid) or _MULTI_GUIDE_RE.search(sym):
                    n_multi_or_zero_excl += 1
                    continue
                if gid.lower() in _UNASSIGNED_TOKENS and sym.lower() in _UNASSIGNED_TOKENS:
                    n_multi_or_zero_excl += 1
                    continue
                if _is_safe_targeting_r14(gid, sym, cls_s) or gid in safe_keys:
                    n_safe_excl += 1
                    continue
                if _is_non_targeting_r14(gid, sym, cls_s) or gid in ntc_keys:
                    keep_rows.append(r_i)
                    nt_flags.append(True)
                    pkeys.append("")
                    continue
                info = guide_to_target.get(gid)
                if info is None and sym and _is_plain_single_gene_symbol(sym):
                    canon_sym, hstat = hgnc.resolve(sym, "")
                    if canon_sym in hgnc.approved:
                        info = (hgnc.sym_to_ens.get(canon_sym, ""), canon_sym, sym, hstat)
                if info is not None:
                    keep_rows.append(r_i)
                    nt_flags.append(False)
                    pkeys.append(info[1])
                    pmeta.setdefault(info[1], info)
                else:
                    n_unmapped_excl += 1
        else:
            C_qc = rec["crp_qc_csc"]
            crp_keys = rec["crp_keys"]
            assert C_qc is not None and thrs is not None
            t1, t2, a1 = extract_top2_csc(C_qc)
            pos_mask = C_qc.data >= thrs[C_qc.indices]
            pos_csc = sp.csc_matrix((pos_mask.astype(np.int8), C_qc.indices, C_qc.indptr), shape=C_qc.shape)
            n_pos = np.asarray(pos_csc.sum(axis=0)).ravel()
            s1 = (n_pos == 1) & (t1 >= float(min_guide_thr)) & ((t2 < 3.0) | (t1 / np.maximum(t2, 0.5) >= 3.0))
            for loc_j, raw_r in enumerate(qc_cols):
                if not bool(s1[loc_j]):
                    n_multi_or_zero_excl += 1
                    continue
                gcall = crp_keys[int(a1[loc_j])]
                if gcall in safe_keys or _is_safe_targeting_r14(gcall, gcall, ""):
                    n_safe_excl += 1
                    continue
                if gcall in ntc_keys or _is_non_targeting_r14(gcall, gcall, ""):
                    keep_rows.append(int(raw_r))
                    nt_flags.append(True)
                    pkeys.append("")
                    continue
                info = guide_to_target.get(gcall)
                if info is not None:
                    keep_rows.append(int(raw_r))
                    nt_flags.append(False)
                    pkeys.append(info[1])
                    pmeta.setdefault(info[1], info)
                else:
                    n_unmapped_excl += 1

        if not keep_rows:
            continue

        keep_arr = np.asarray(keep_rows, dtype=np.int64)
        order = np.argsort(keep_arr, kind="mergesort")
        keep_arr = keep_arr[order]
        nt_arr = np.asarray(nt_flags, dtype=bool)[order]
        pkeys_ord = [pkeys[int(i)] for i in order]

        with h5py.File(rec["path"], "r") as f:
            x_node = _get_h5ad_slot_node(f, rec["slot"])
            x_kept = _extract_kept_h5ad_csr(
                x_node=x_node,
                layout=rec["layout"],
                keep_arr=keep_arr,
                col_rows=rec["rna_rows"],
                n_obs=rec["n_obs"],
                n_var=rec["n_var"],
            )

        all_kept_parts.append(
            {
                "rec": rec,
                "keep_arr": keep_arr,
                "nt_arr": nt_arr,
                "pkeys": pkeys_ord,
                "pmeta": pmeta,
                "x_kept": x_kept,
            }
        )

    if len(all_kept_parts) < 2:
        raise RuntimeError(f"EXCLUDED_FEASIBILITY: {line} did not yield kept cells across >=2 channels.")

    tot_kept = sum(len(x["keep_arr"]) for x in all_kept_parts)
    global_nt = np.concatenate([x["nt_arr"] for x in all_kept_parts])
    # PREREG_G Item 5: ch_1 -> Half A (1), ch_2 -> Half B (2); channel-parity generalizable to >=2 channels.
    global_half = np.concatenate(
        [
            np.full(len(x["keep_arr"]), 1 if int(x["rec"]["ch_num"]) % 2 == 1 else 2, dtype=np.int8)
            for x in all_kept_parts
        ]
    )
    if np.sum(global_nt & (global_half == 1)) == 0 or np.sum(global_nt & (global_half == 2)) == 0:
        raise RuntimeError(f"EXCLUDED_FEASIBILITY: {line} has zero NTC cells in Half A or Half B.")

    raw_batches = np.concatenate(
        [np.full(len(x["keep_arr"]), str(x["rec"]["ch_tag"]), dtype=object) for x in all_kept_parts]
    )
    global_batches, n_fb = _apply_min_ntc_batch_fallback(
        raw_batches, global_nt, global_half, min_batch_nt=min_batch_nt
    )

    parts: List[PartSpec] = []
    offset = 0
    for item in all_kept_parts:
        rec = item["rec"]
        klen = len(item["keep_arr"])
        sl = slice(offset, offset + klen)
        offset += klen
        n_genes_rna = len(rec["rna_rows"])
        rc, rr, rcc = _make_csr_part_callbacks(item["x_kept"], n_genes_rna)
        parts.append(
            PartSpec(
                tag=f"{line}_{rec['ch_tag']}",
                path=rec["path"],
                n_cells=klen,
                n_genes_raw=n_genes_rna,
                nt_mask=item["nt_arr"],
                batches=global_batches[sl],
                half_flags=global_half[sl],
                umi=rec["umi_all"][item["keep_arr"]],
                pert_keys=item["pkeys"],
                pert_meta=item["pmeta"],
                var_gene_ids=rec["r_ids"],
                var_gene_names=rec["r_names"],
                read_chunk_fn=rc,
                read_rows_fn=rr,
                read_chunk_cols_fn=rcc,
            )
        )

    audit = {
        "line": line,
        "n_h5ad_files": len(ch_entries),
        "channel_slots": {
            rec["ch_tag"]: {
                "x_slot_used": rec["slot"],
                "selected_median_total_umi": float(rec["slot_diag"]["slots"][rec["slot"]]["median_total_umi"]),
                "int_slots": rec["int_slots"],
                "slot_diag": rec["slot_diag"],
            }
            for rec in ch_records
        },
        "sgrna_lib_source": lib_meta["source"],
        "guide_call_mode": "h5ad_crispr_em_fallback" if using_fallback_em else "obs_author_single_guide",
        "min_umi": int(min_umi),
        "min_ngenes": int(min_ngenes),
        "max_pct_mt": (None if max_pct_mt is None else float(max_pct_mt)),
        "n_rna_qc_excluded": int(n_rna_qc_excl),
        "n_multi_or_zero_guide_excluded": int(n_multi_or_zero_excl),
        "n_safe_targeting_excluded": int(n_safe_excl),
        "n_unmapped_guide_excluded": int(n_unmapped_excl),
        "n_qc_kept_cells": int(tot_kept),
        "n_kept_ntc": int(global_nt.sum()),
        "n_kept_half_a_ntc": int(np.sum(global_nt & (global_half == 1))),
        "n_kept_half_b_ntc": int(np.sum(global_nt & (global_half == 2))),
        "n_low_ntc_channel_fallback_cells": int(n_fb),
        "half_rule_used": "ch_1_HalfA_vs_ch_2_HalfB",
        "channels_used": sorted({str(b) for b in global_batches}),
        "median_guide_threshold": (None if thrs is None else float(np.median(thrs))),
    }
    return parts, audit


def prepare_r14_line(
    line: str,
    raw_dir: Path,
    out_root: Path,
    hgnc: Optional[HGNCLookup] = None,
    workers: int = 8,
    chunk_size: int = 8192,
    x_layer: str = "auto",
) -> Path:
    """Standardize `GSE337804` (`a549_myllia`) cleanly via frozen `process_r10_line`."""
    hgnc_l = hgnc or build_hgnc_lookup()
    parts, audit = build_gse337804_parts(raw_dir=raw_dir, hgnc=hgnc_l, line=line, x_layer=x_layer)
    return process_r10_line(
        line_name=line,
        parts=parts,
        hgnc=hgnc_l,
        chunk_size=chunk_size,
        out_root=out_root,
        extra_meta={"adapter_audit": audit},
        n_workers=workers,
    )


__all__ = [
    "build_gse337804_parts",
    "check_gse337804_feasibility",
    "load_gse337804_sgrna_lib",
    "prepare_r14_line",
]
