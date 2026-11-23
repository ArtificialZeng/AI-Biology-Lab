#!/usr/bin/env python3
"""
Independent Round 14 GSE337804 (`a549_myllia`) Reader, Slot Checker, & Blinding Gate (A-T3).

Strictly implements locked `PREREG_G.md` (`commit d34560d`):
  - Clean-room isolation: NEVER reads, imports, or references `round14/src/`,
    `round14/results/`, or `round14/status/A-T2.txt`.
  - Blinding gates (`enforce_r14_gate`):
      * `synthetic_test`: synthetic `.h5ad` fixtures ONLY; strictly blocks `round14/shared/`.
      * `dev_r9..dev_r13`: historical unblinded caches only; blocks `round14/shared/`.
      * `reader_xcheck_r14`: raw `GSE337804` reading permitted ONLY after BOTH
        `round14/status/PREREG_G_COMMIT.txt` AND `round14/status/R14_LINES_AUTHORIZED.txt`
        exist and are non-empty. Never downloads data.
      * `holdout_r14`: evaluation on `round14/shared/lines/` permitted ONLY after
        `PREREG_G_COMMIT.txt`, `CODE_FREEZE_R14.txt`, and `CONFIRMATORY_LINES_R14.txt`.
  - Pre-metric technical feasibility & slot verification (`§4.1` & Multi-Slot Stop Guard):
      * Checks `/raw/X`, `/X`, and all non-excluded `/layers/*` (excluding subfractions/
        normalized layers: `nascent`, `unspliced`, `spliced`, `ambiguous`, `velocity`,
        `intron`, `exon`, `lognorm`, `log1p`, `scale`, `norm`, `cpm`, `tpm`, `vst`, `sct`, etc.).
      * Verifies **100% of stored non-zero entries** (streamed chunk-by-chunk, zero sampling)
        are finite, non-negative exact integers (`>= 0`, `x == floor(x)`).
      * If `>= 2` non-excluded slots pass integer raw-count verification, checks exact
        shape, `/var` identifier/symbol list order, and bit-exact matrix values across ALL rows;
        if ANY discrepancy exists, **STOPS with `EXCLUDED_FEASIBILITY`**
        (`MULTIPLE_INCONSISTENT_INTEGER_SLOTS`) — never chooses one arbitrarily.
  - Guide / NTC / Target / QC / Split rules (`§4.2–§4.5, §5`):
      * Author `.obs` single-sgRNA calls (`sgrna_id` / `sgrna_symbol`) matched to
        `GSE337804_sgRNA_lib.csv.gz`; excludes multi-guide delimiters, `nperts > 1`,
        unassigned, multiplet, and empty entries. Falls back to `PREREG_F §3.2`
        Poisson-Gaussian histogram thresholds (`n_genes >= 200`, `min_thr=5, max_thr=20,
        default_thr=7`) ONLY if author `.obs` calls are absent and raw guide counts exist.
      * Strictly retains ONLY the `20` `ctrl_NT` (`non-targeting_1..20`, `symbol == "non-targeting"`)
        as `NTC`; strictly EXCLUDES all `20` `ctrl_T` (`gene_desert_1..20`,
        `symbol == "gene_desert"`) safe-targeting cells from both NTC and targets.
      * Resolves targeting genes through frozen `hgnc_complete_set.txt` (`HGNCLookup.resolve`)
        and pools all QC-valid single-sgRNA cells per canonical HGNC gene (`pert_key = canon_sym`).
      * Cell QC: `n_umis >= 300` AND `n_genes >= 200` with **`max_pct_mt = None`** (`DepleteX`).
      * Disjoint split-half & NTC subtraction: `ch_1 -> Half A (1)`, `ch_2 -> Half B (2)`,
        `batches = ["ch_1", "ch_2"]` (`>= 15` NTC/channel, with safe pooled fallback).
"""

from __future__ import annotations

import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Set, Tuple

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROUND14_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = ROUND14_ROOT.parent
R10_SCRIPTS_DIR = REPO_ROOT / "round10" / "shared" / "scripts"
R11_AUDIT_DIR = REPO_ROOT / "round11" / "audit"
R11_SRC_DIR = REPO_ROOT / "round11" / "src"
R12_AUDIT_DIR = REPO_ROOT / "round12" / "audit"
R13_AUDIT_DIR = REPO_ROOT / "round13" / "audit"

for _p in (str(R10_SCRIPTS_DIR), str(R11_AUDIT_DIR), str(R11_SRC_DIR), str(R12_AUDIT_DIR), str(R13_AUDIT_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from dataset_readers_r11 import _apply_min_ntc_batch_fallback  # noqa: E402
from independent_r11_loader import (  # noqa: E402
    DEFAULT_HGNC_TSV,
    R9_LINE_NAMES,
    R9_LINE_SET,
    R9_LINES_DIR,
    R10_LINES_DIR,
    R10_SOURCE_EXT_NAMES,
    R11_LINES_DIR,
    AuditLineData,
    CanonicalHGNC,
    _first_token,
    load_line_dataset,
)
from independent_r12_loader import R12_LINES_DIR, R12_TARGET_KEYS  # noqa: E402
from independent_r13_loader import R13_LINES_DIR, R13_TARGET_KEYS  # noqa: E402
from prepare_r10_lines import (  # noqa: E402
    HGNCLookup,
    PartSpec,
    build_hgnc_lookup,
    extract_top2_csc,
    fit_poisson_gaussian_thresholds_from_hist,
    process_r10_line,
)

R14_SHARED_DIR = ROUND14_ROOT / "shared"
R14_RAW_DIR = R14_SHARED_DIR / "raw"
R14_LINES_DIR = R14_SHARED_DIR / "lines"
R14_STATUS_DIR = ROUND14_ROOT / "status"

R14_TARGET_KEY = "a549_myllia"
R14_TARGET_KEYS: Tuple[str, ...] = (R14_TARGET_KEY,)
R14_TARGET_SET: Set[str] = {R14_TARGET_KEY}

R10_ALL_LINES: Tuple[str, ...] = ("dld1", "cd4t_rest", "wtc11_i3neuron", "calu3")
R11_ALL_LINES: Tuple[str, ...] = ("telohaec", "cd4t_stim48", "feng_ipsc")
R12_ALL_LINES: Tuple[str, ...] = R12_TARGET_KEYS
R13_ALL_LINES: Tuple[str, ...] = R13_TARGET_KEYS

DEV_ALLOWED_LINES_BY_PHASE: Dict[str, Tuple[str, ...]] = {
    "dev_r9": R9_LINE_NAMES,
    "dev_r10": R10_ALL_LINES,
    "dev_r11": R11_ALL_LINES,
    "dev_r12": R12_ALL_LINES,
    "dev_r13": R13_ALL_LINES,
}
ALL_DEV_LINES: Set[str] = (
    set(R9_LINE_NAMES) | set(R10_ALL_LINES) | set(R11_ALL_LINES) | set(R12_ALL_LINES) | set(R13_ALL_LINES)
)

EXCLUDED_LAYER_TOKENS: Tuple[str, ...] = (
    "nascent",
    "unspliced",
    "spliced",
    "ambiguous",
    "velocity",
    "intron",
    "exon",
    "lognorm",
    "log1p",
    "scale",
    "scaled",
    "norm",
    "cpm",
    "tpm",
    "vst",
    "sct",
    "pearson",
    "residual",
    "pca",
    "umap",
    "data",
)

_MULTI_GUIDE_DELIM_RE = re.compile(r"[,;|+\s]")
_BAD_GUIDE_TOKENS: Set[str] = {
    "",
    "na",
    "nan",
    "none",
    "null",
    "unassigned",
    "multiplet",
    "doublet",
    "ambiguous",
    "negative",
    "empty",
    "undetected",
}


def _assert_not_inside_r14_shared(path: Path, label: str) -> None:
    r14_shared = R14_SHARED_DIR.resolve()
    resolved = path.resolve()
    if resolved == r14_shared or r14_shared in resolved.parents:
        raise RuntimeError(f"PREREG_G blinding lock ({label}): refused accessing {resolved} inside {r14_shared}.")


def enforce_r14_gate(
    phase: str,
    target_lines: Sequence[str] = (),
    status_dir: Path = R14_STATUS_DIR,
    extra_dirs: Sequence[Path] = (),
    require_confirmatory_file: bool = True,
) -> Dict[str, str]:
    """
    Enforce `PREREG_G.md §1–§2` blinding & pre-freeze locks before accessing expression matrices:
      - `synthetic_test`: synthetic fixtures ONLY; strictly blocks `round14/shared/`.
      - `dev_r9..dev_r13`: historical unblinded caches only; blocks `round14/shared/` and `a549_myllia`.
      - `reader_xcheck_r14`: raw `GSE337804` reading requires BOTH non-empty
        `PREREG_G_COMMIT.txt` AND `R14_LINES_AUTHORIZED.txt` in `status_dir`.
      - `holdout_r14`: evaluation on `round14/shared/lines/` requires `PREREG_G_COMMIT.txt`,
        `CODE_FREEZE_R14.txt`, and `CONFIRMATORY_LINES_R14.txt`.
    """
    if phase == "synthetic_test":
        for ed in extra_dirs:
            _assert_not_inside_r14_shared(ed, phase)
        return {"phase": phase, "mode": "SYNTHETIC_FIXTURE_ONLY"}

    if phase in DEV_ALLOWED_LINES_BY_PHASE:
        bad_targets = [t for t in target_lines if (t in R14_TARGET_SET or t not in ALL_DEV_LINES)]
        if bad_targets:
            raise RuntimeError(
                f"PREREG_G dev blinding lock ({phase}): refused non-dev / Round 14 targets {bad_targets}."
            )
        for ed in extra_dirs:
            _assert_not_inside_r14_shared(ed, phase)
        return {"phase": phase, "mode": "DEV_HISTORICAL_UNBLINDED"}

    prereg_hash = _first_token(status_dir / "PREREG_G_COMMIT.txt")
    if not prereg_hash:
        raise RuntimeError(
            f"PREREG_G_COMMIT hard gate blocked ({phase}): missing or empty {status_dir / 'PREREG_G_COMMIT.txt'}."
        )

    if phase == "reader_xcheck_r14":
        auth_tok = _first_token(status_dir / "R14_LINES_AUTHORIZED.txt")
        if not auth_tok:
            raise RuntimeError(
                f"R14_LINES_AUTHORIZED hard gate blocked ({phase}): missing or empty "
                f"{status_dir / 'R14_LINES_AUTHORIZED.txt'}."
            )
        return {
            "phase": phase,
            "prereg_commit": prereg_hash,
            "r14_lines_authorized": auth_tok,
            "code_freeze": _first_token(status_dir / "CODE_FREEZE_R14.txt"),
        }

    if phase == "holdout_r14":
        freeze_hash = _first_token(status_dir / "CODE_FREEZE_R14.txt")
        if not freeze_hash:
            raise RuntimeError(
                f"CODE_FREEZE_R14 hard gate blocked: missing or empty {status_dir / 'CODE_FREEZE_R14.txt'}."
            )
        conf_p = status_dir / "CONFIRMATORY_LINES_R14.txt"
        if require_confirmatory_file and (not conf_p.exists() or not conf_p.read_text(encoding="utf-8").strip()):
            raise RuntimeError(f"CONFIRMATORY_LINES_R14 gate blocked: missing or empty {conf_p}.")
        return {
            "phase": phase,
            "prereg_commit": prereg_hash,
            "code_freeze": freeze_hash,
            "r14_lines_authorized": _first_token(status_dir / "R14_LINES_AUTHORIZED.txt"),
        }

    raise ValueError(f"Unsupported Round 14 phase: {phase!r}")


def parse_confirmatory_lines_r14(conf_path: Path) -> Dict[str, Any]:
    """
    Parse `CONFIRMATORY_LINES_R14.txt` deterministically under `PREREG_G.md` (`max K_conf = 1`).
    If `a549_myllia` is marked `EXCLUDED_FEASIBILITY`, `FAILED_GATE`, `FAILED_COUNT`, or
    `INGESTION_FAILED`, confirmatory adjudication aborts (`k_conf = 0`).
    """
    if not conf_path.exists():
        raise FileNotFoundError(f"Missing confirmatory gate file: {conf_path}")
    text = conf_path.read_text(encoding="utf-8")
    if not text.strip():
        raise RuntimeError(f"Empty confirmatory gate file: {conf_path}")

    listed_conf: List[str] = []
    a549_role: Optional[str] = None
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        upper = line.upper()
        if upper.startswith("CONFIRMATORY_LINES"):
            rhs = line.split("=", 1)[-1] if "=" in line else line.split(":", 1)[-1]
            for tok in rhs.replace(",", " ").split():
                t = tok.strip().lower()
                if t == R14_TARGET_KEY and t not in listed_conf:
                    listed_conf.append(t)
        elif upper.startswith("A549_MYLLIA_ROLE") or line.lower().startswith("a549_myllia="):
            rhs = line.split("=", 1)[-1] if "=" in line else line.split(":", 1)[-1]
            a549_role = rhs.strip().split()[0].upper()
        elif line.lower() == R14_TARGET_KEY and R14_TARGET_KEY not in listed_conf:
            listed_conf.append(R14_TARGET_KEY)

    upper_all = text.upper()
    if a549_role is None:
        if "EXCLUDED_FEASIBILITY" in upper_all:
            a549_role = "EXCLUDED_FEASIBILITY"
        elif "FAILED_COUNT" in upper_all or "FAILED_GATE" in upper_all:
            a549_role = "FAILED_GATE"
        elif "INGESTION_FAILED" in upper_all:
            a549_role = "INGESTION_FAILED"
        elif R14_TARGET_KEY in listed_conf:
            a549_role = "CONFIRMATORY"
        else:
            a549_role = "FAILED_GATE"

    confirmatory_lines = [R14_TARGET_KEY] if a549_role == "CONFIRMATORY" else []
    secondary_eval_lines = (
        [R14_TARGET_KEY]
        if (a549_role not in ("CONFIRMATORY", "EXCLUDED_FEASIBILITY", "INGESTION_FAILED"))
        else []
    )
    return {
        "confirmatory_lines": confirmatory_lines,
        "k_conf": len(confirmatory_lines),
        "a549_myllia_role": a549_role,
        "excluded_feasibility": a549_role == "EXCLUDED_FEASIBILITY",
        "aborted_confirmatory": len(confirmatory_lines) == 0,
        "secondary_eval_lines": secondary_eval_lines,
    }


def _decode_str(arr: Any) -> List[str]:
    out: List[str] = []
    for x in arr:
        if isinstance(x, (bytes, np.bytes_)):
            out.append(x.decode("utf-8", errors="replace"))
        else:
            out.append(str(x))
    return out


def _read_obs_or_var_column(group: h5py.Group, col_name: str) -> Optional[np.ndarray]:
    """Read an AnnData HDF5 `/obs` or `/var` column handling categorical and flat encodings."""
    if col_name not in group:
        return None
    node = group[col_name]
    if isinstance(node, h5py.Group) and "categories" in node and "codes" in node:
        cats = np.asarray(_decode_str(node["categories"][()]), dtype=object)
        codes = np.asarray(node["codes"][()], dtype=np.int64)
        res = np.empty(len(codes), dtype=object)
        valid = (codes >= 0) & (codes < len(cats))
        res[valid] = cats[codes[valid]]
        res[~valid] = ""
        return res
    raw = node[()]
    if "__categories" in group and col_name in group["__categories"] and raw.dtype.kind not in ("S", "O", "U"):
        cats = np.asarray(_decode_str(group["__categories"][col_name][()]), dtype=object)
        codes = np.asarray(raw, dtype=np.int64)
        res = np.empty(len(codes), dtype=object)
        valid = (codes >= 0) & (codes < len(cats))
        res[valid] = cats[codes[valid]]
        res[~valid] = ""
        return res
    if raw.dtype.kind in ("S", "O", "U"):
        return np.asarray(_decode_str(raw), dtype=object)
    return np.asarray(raw)


def _list_group_columns(group: h5py.Group) -> List[str]:
    if "column-order" in group.attrs:
        return _decode_str(group.attrs["column-order"])
    return sorted(k for k in group.keys() if k != "__categories")


def _read_anndata_index(group: h5py.Group, group_name: str) -> List[str]:
    """Read AnnData HDF5 `/obs` or `/var` index (`attrs["_index"]` -> `"_index"` -> `"index"`)."""
    candidates: List[str] = []
    if "_index" in group.attrs:
        raw_attr = group.attrs["_index"]
        if isinstance(raw_attr, (bytes, np.bytes_)):
            idx_name = raw_attr.decode("utf-8", errors="replace").strip()
        elif isinstance(raw_attr, np.ndarray):
            if raw_attr.ndim == 0:
                scalar = raw_attr.item()
                idx_name = (
                    scalar.decode("utf-8", errors="replace").strip()
                    if isinstance(scalar, (bytes, np.bytes_))
                    else str(scalar).strip()
                )
            elif raw_attr.size > 0:
                idx_name = _decode_str([raw_attr.flat[0]])[0].strip()
            else:
                idx_name = ""
        else:
            idx_name = str(raw_attr).strip()
        if idx_name:
            candidates.append(idx_name)
    for fallback_name in ("_index", "index"):
        if fallback_name not in candidates:
            candidates.append(fallback_name)
    for cname in candidates:
        arr = _read_obs_or_var_column(group, cname)
        if arr is not None:
            return [str(x).strip() for x in arr]
    raise RuntimeError(
        f"MISSING_ANNDATA_INDEX: /{group_name} lacks valid attrs['_index'], '_index', or 'index' dataset."
    )


def _extract_var_metadata(var: h5py.Group, group_label: str) -> Dict[str, Any]:
    var_index = _read_anndata_index(var, group_label)
    var_cols = _list_group_columns(var)

    v_sym: Optional[List[str]] = None
    for vc in ("gene_symbol", "gene_name", "symbol", "feature_name", "gene_short_name"):
        tmp_v = _read_obs_or_var_column(var, vc)
        if tmp_v is not None:
            v_sym = [str(x).strip() for x in tmp_v]
            break
    if v_sym is None:
        v_sym = list(var_index)

    v_ens: Optional[List[str]] = None
    for vc in ("ensembl_id", "gene_id", "gene_ids", "ensembl_gene_id", "feature_id"):
        tmp_v = _read_obs_or_var_column(var, vc)
        if tmp_v is not None:
            v_ens = [str(x).strip().split(".")[0] for x in tmp_v]
            break
    if v_ens is None:
        v_ens = [x.split(".")[0] if str(x).startswith("ENSG") else "" for x in var_index]

    v_ftype: Optional[List[str]] = None
    for fc in ("feature_types", "feature_type"):
        tmp_f = _read_obs_or_var_column(var, fc)
        if tmp_f is not None:
            v_ftype = [str(x).strip() for x in tmp_f]
            break

    return {
        "var_index": var_index,
        "var_cols": var_cols,
        "var_syms": v_sym,
        "var_ens": v_ens,
        "var_ftypes": v_ftype,
    }


def _find_sgrna_lib_path(raw_dir: Path) -> Path:
    candidates = [
        raw_dir / "a549_myllia" / "GSE337804_sgRNA_lib.csv.gz",
        raw_dir / "a549_myllia" / "GSE337804_sgRNA_lib.csv",
        raw_dir / "GSE337804_sgRNA_lib.csv.gz",
        raw_dir / "GSE337804_sgRNA_lib.csv",
    ]
    for c in candidates:
        if c.exists():
            return c
    for pattern in ("*sgRNA_lib*.csv*", "*sgrna_lib*.csv*"):
        hits = sorted((raw_dir / "a549_myllia").glob(pattern)) if (raw_dir / "a549_myllia").exists() else []
        if not hits:
            hits = sorted(raw_dir.glob(pattern))
        if hits:
            return hits[0]
    raise FileNotFoundError(f"Missing GSE337804_sgRNA_lib.csv(.gz) under {raw_dir}")


def discover_gse337804_h5ad_channels(raw_dir: Path) -> List[Tuple[int, str, Path]]:
    """
    Locate strictly the two `GSE337804` A549 `.h5ad` channels (`ch_1`, `ch_2`).
    Channel mapping per `PREREG_G §5`: `ch_1 -> Half A (1)`, `ch_2 -> Half B (2)`.
    """
    search_dirs = [d for d in (raw_dir / "a549_myllia", raw_dir) if d.exists()]
    ch1_path: Optional[Path] = None
    ch2_path: Optional[Path] = None

    for d in search_dirs:
        for exact1 in ("GSE337804_A549_ch_1.RDS.h5ad", "GSE337804_A549_ch_1.h5ad"):
            if (d / exact1).exists():
                ch1_path = d / exact1
                break
        for exact2 in ("GSE337804_A549_ch_2.RDS.h5ad", "GSE337804_A549_ch_2.h5ad"):
            if (d / exact2).exists():
                ch2_path = d / exact2
                break
        if ch1_path is None:
            hits1 = sorted(d.glob("*ch_1*.h5ad")) + sorted(d.glob("*ch1*.h5ad"))
            if hits1:
                ch1_path = hits1[0]
        if ch2_path is None:
            hits2 = sorted(d.glob("*ch_2*.h5ad")) + sorted(d.glob("*ch2*.h5ad"))
            if hits2:
                ch2_path = hits2[0]
        if ch1_path is not None and ch2_path is not None:
            break

    if ch1_path is None or ch2_path is None:
        raise FileNotFoundError(
            f"Missing GSE337804 ch_1/ch_2 .h5ad files under {raw_dir}: ch_1={ch1_path}, ch_2={ch2_path}"
        )
    if ch1_path.resolve() == ch2_path.resolve():
        raise RuntimeError(f"ch_1 and ch_2 resolved to identical file: {ch1_path}")
    return [(1, "ch_1", ch1_path), (2, "ch_2", ch2_path)]


def load_gse337804_sgrna_lib(lib_path: Path, hgnc: HGNCLookup) -> Dict[str, Any]:
    """
    Parse `GSE337804_sgRNA_lib.csv.gz` strictly per `PREREG_G.md §4.3–§4.4`:
      - `class == "ctrl_NT"` (`non-targeting_1..20`, `symbol == "non-targeting"`): NTC.
      - `class == "ctrl_T"` (`gene_desert_1..20`, `symbol == "gene_desert"`): strictly
        EXCLUDED (neither NTC nor target).
      - `class == "gene"`: resolved through frozen `hgnc.resolve(symbol, "")`; only
        valid HGNC-mapped genes are kept as targeting perturbations (`pert_key = canon_sym`).
    """
    if not lib_path.exists():
        raise FileNotFoundError(f"Missing GSE337804 sgRNA library file: {lib_path}")

    open_fn = gzip.open if lib_path.name.endswith(".gz") else open
    with open_fn(lib_path, "rt", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    required = {"id", "symbol", "class"}
    if not rows or not required.issubset(set(rows[0].keys())):
        raise ValueError(f"Invalid GSE337804 sgRNA library columns in {lib_path}")

    ntc_guide_ids: Set[str] = set()
    ntc_symbols: Set[str] = set()
    excluded_safe_ids: Set[str] = set()
    excluded_safe_symbols: Set[str] = set()
    guide_to_target: Dict[str, Tuple[str, str, str, str]] = {}
    symbol_to_target: Dict[str, Tuple[str, str, str, str]] = {}
    unmapped_guide_ids: Dict[str, str] = {}
    ordered_guide_ids: List[str] = []

    for row in rows:
        gid = str(row["id"]).strip()
        sym = str(row["symbol"]).strip()
        cls = str(row["class"]).strip()
        ordered_guide_ids.append(gid)

        if cls == "ctrl_NT":
            if sym.lower() in ("gene_desert", "safe-targeting", "safe_targeting") or gid.lower().startswith("gene_desert"):
                raise ValueError(f"Library integrity error: ctrl_NT row looks safe-targeting: {row}")
            ntc_guide_ids.add(gid)
            ntc_symbols.add(sym)
        elif cls == "ctrl_T":
            excluded_safe_ids.add(gid)
            excluded_safe_symbols.add(sym)
        elif cls == "gene":
            canon_sym, hstat = hgnc.resolve(sym, "")
            is_mapped = bool(canon_sym) and (hstat != "unmapped" or canon_sym in hgnc.approved)
            if not is_mapped:
                unmapped_guide_ids[gid] = sym
                continue
            ens_id = hgnc.sym_to_ens.get(canon_sym, "")
            meta = (ens_id, canon_sym, sym, hstat)
            guide_to_target[gid] = meta
            symbol_to_target.setdefault(sym, meta)
        else:
            raise ValueError(f"Unexpected GSE337804 guide class {cls!r} in row {row}")

    if not ntc_guide_ids:
        raise ValueError(f"No ctrl_NT non-targeting guides parsed from {lib_path}")
    if not guide_to_target:
        raise ValueError(f"No mapped targeting guides parsed from {lib_path}")

    return {
        "ntc_guide_ids": ntc_guide_ids,
        "ntc_symbols": ntc_symbols,
        "excluded_safe_ids": excluded_safe_ids,
        "excluded_safe_symbols": excluded_safe_symbols,
        "guide_to_target": guide_to_target,
        "symbol_to_target": symbol_to_target,
        "unmapped_guide_ids": unmapped_guide_ids,
        "ordered_guide_ids": ordered_guide_ids,
    }


def _is_excluded_layer_name(layer_name: str) -> bool:
    low = layer_name.strip().lower()
    if low in ("counts", "raw_counts", "raw", "umi", "umi_counts", "total_counts", "x_raw"):
        return False
    return any(tok in low for tok in EXCLUDED_LAYER_TOKENS)


def _detect_matrix_encoding(node: Any, var_len: int, obs_len: int) -> Tuple[str, int, int]:
    enc_attr = node.attrs.get("encoding-type", "")
    if isinstance(enc_attr, bytes):
        enc_attr = enc_attr.decode("utf-8", errors="ignore")
    if isinstance(node, h5py.Dataset):
        if node.ndim != 2:
            raise ValueError(f"Expected 2D matrix dataset, got shape={node.shape}")
        return "dense", int(node.shape[0]), int(node.shape[1])
    if isinstance(node, h5py.Group) and {"data", "indices", "indptr"}.issubset(set(node.keys())):
        sh = node.attrs.get("shape", None)
        if "csc" in str(enc_attr).lower():
            n_r = int(sh[0]) if sh is not None else int(obs_len)
            n_c = int(sh[1]) if sh is not None else int(len(node["indptr"]) - 1)
            return "csc", n_r, n_c
        n_r = int(sh[0]) if sh is not None else int(len(node["indptr"]) - 1)
        n_c = int(sh[1]) if sh is not None else int(var_len)
        return "csr", n_r, n_c
    raise ValueError("unrecognized_matrix_encoding")


def _full_slot_integer_check(
    f: h5py.File,
    slot_path: str,
    obs_len: int,
    chunk_entries: int = 2_000_000,
    chunk_rows: int = 4096,
) -> Dict[str, Any]:
    """
    Verify **100% of stored entries** across an HDF5 matrix slot (`/raw/X`, `/X`, or `/layers/*`)
    are finite, non-negative exact integers (`>= 0`, `x == floor(x)`). Zero sampling!
    """
    h5_key = slot_path.lstrip("/")
    if h5_key not in f:
        return {"slot_path": slot_path, "is_integer": False, "reason": "missing"}

    var_key = "raw/var" if slot_path == "/raw/X" else "var"
    if var_key not in f or not isinstance(f[var_key], h5py.Group):
        return {"slot_path": slot_path, "is_integer": False, "reason": f"missing_{var_key}"}

    var_meta = _extract_var_metadata(f[var_key], var_key)
    node = f[h5_key]
    try:
        encoding, n_rows, n_cols = _detect_matrix_encoding(node, len(var_meta["var_index"]), obs_len)
    except Exception as exc:
        return {"slot_path": slot_path, "is_integer": False, "reason": str(exc)}

    if n_rows != obs_len or n_cols != len(var_meta["var_index"]):
        return {
            "slot_path": slot_path,
            "is_integer": False,
            "reason": (
                f"dimension_mismatch: matrix=({n_rows}, {n_cols}) vs "
                f"(obs={obs_len}, {var_key}={len(var_meta['var_index'])})"
            ),
        }

    nnz_pos = 0
    sum_total = 0
    max_val = 0
    if encoding == "dense":
        dtype = node.dtype
        is_int_dt = bool(np.issubdtype(dtype, np.integer))
        for r0 in range(0, n_rows, chunk_rows):
            r1 = min(n_rows, r0 + chunk_rows)
            vals = np.asarray(node[r0:r1, :])
            if vals.size == 0:
                continue
            if not is_int_dt and np.any(~np.isfinite(vals)):
                return {"slot_path": slot_path, "is_integer": False, "reason": f"nonfinite_in_rows_{r0}_{r1}"}
            if np.any(vals < 0):
                return {"slot_path": slot_path, "is_integer": False, "reason": f"negative_in_rows_{r0}_{r1}"}
            if not is_int_dt:
                if np.any(vals != np.floor(vals)):
                    max_frac = float(np.max(np.abs(vals - np.rint(vals))))
                    return {
                        "slot_path": slot_path,
                        "is_integer": False,
                        "reason": f"non_integer_float_in_rows_{r0}_{r1}_maxdiff={max_frac:.4g}",
                    }
            pos = vals > 0
            if np.any(pos):
                pvals = vals[pos].astype(np.int64, copy=False)
                nnz_pos += int(pvals.size)
                sum_total += int(np.sum(pvals, dtype=np.int64))
                cmax = int(np.max(pvals))
                if cmax > max_val:
                    max_val = cmax
    else:
        data_ds = node["data"]
        dtype = data_ds.dtype
        total_stored = int(data_ds.shape[0])
        if total_stored == 0:
            return {"slot_path": slot_path, "is_integer": False, "reason": "empty_sparse_matrix"}
        is_int_dt = bool(np.issubdtype(dtype, np.integer))
        for s0 in range(0, total_stored, chunk_entries):
            s1 = min(total_stored, s0 + chunk_entries)
            vals = np.asarray(data_ds[s0:s1])
            if not is_int_dt and np.any(~np.isfinite(vals)):
                return {"slot_path": slot_path, "is_integer": False, "reason": f"sparse_nonfinite_{s0}_{s1}"}
            if np.any(vals < 0):
                return {"slot_path": slot_path, "is_integer": False, "reason": f"sparse_negative_{s0}_{s1}"}
            if not is_int_dt:
                if np.any(vals != np.floor(vals)):
                    max_frac = float(np.max(np.abs(vals - np.rint(vals))))
                    return {
                        "slot_path": slot_path,
                        "is_integer": False,
                        "reason": f"sparse_non_integer_{s0}_{s1}_maxdiff={max_frac:.4g}",
                    }
            pos = vals > 0
            if np.any(pos):
                pvals = vals[pos].astype(np.int64, copy=False)
                nnz_pos += int(pvals.size)
                sum_total += int(np.sum(pvals, dtype=np.int64))
                cmax = int(np.max(pvals))
                if cmax > max_val:
                    max_val = cmax

    if nnz_pos == 0 or max_val < 1:
        return {"slot_path": slot_path, "is_integer": False, "reason": "all_zero_matrix"}

    return {
        "slot_path": slot_path,
        "var_group": var_key,
        "is_integer": True,
        "encoding": encoding,
        "shape": (n_rows, n_cols),
        "dtype": str(dtype),
        "nnz_pos": int(nnz_pos),
        "sum_total": int(sum_total),
        "max_val": int(max_val),
        "var_meta": var_meta,
        "reason": f"full_{encoding}_int_{dtype}_nnz={nnz_pos}_max={max_val}",
    }


def _read_slot_rows_as_csr(
    node: Any,
    encoding: str,
    n_rows: int,
    n_cols: int,
    r0: int,
    r1: int,
) -> sp.csr_matrix:
    n_blk = r1 - r0
    if encoding == "dense":
        mat = sp.csr_matrix(np.asarray(node[r0:r1, :], dtype=np.int64))
        mat.eliminate_zeros()
        mat.sum_duplicates()
        mat.sort_indices()
        return mat
    if encoding == "csr":
        ptrs = np.asarray(node["indptr"][r0 : r1 + 1], dtype=np.int64)
        p0, p1 = int(ptrs[0]), int(ptrs[-1])
        if p1 == p0:
            return sp.csr_matrix((n_blk, n_cols), dtype=np.int64)
        mat = sp.csr_matrix(
            (
                np.asarray(node["data"][p0:p1], dtype=np.int64),
                np.asarray(node["indices"][p0:p1], dtype=np.int32),
                ptrs - np.int64(p0),
            ),
            shape=(n_blk, n_cols),
        )
        mat.eliminate_zeros()
        mat.sum_duplicates()
        mat.sort_indices()
        return mat
    if encoding == "csc":
        all_ptrs = np.asarray(node["indptr"][:], dtype=np.int64)
        indices_ds = node["indices"]
        data_ds = node["data"]
        coo_r: List[np.ndarray] = []
        coo_c: List[np.ndarray] = []
        coo_v: List[np.ndarray] = []
        for col_idx in range(n_cols):
            p0, p1 = int(all_ptrs[col_idx]), int(all_ptrs[col_idx + 1])
            if p1 <= p0:
                continue
            r_col = np.asarray(indices_ds[p0:p1], dtype=np.int64)
            mask = (r_col >= r0) & (r_col < r1)
            if np.any(mask):
                coo_r.append((r_col[mask] - r0).astype(np.int32))
                coo_c.append(np.full(int(np.sum(mask)), col_idx, dtype=np.int32))
                coo_v.append(np.asarray(data_ds[p0:p1][mask], dtype=np.int64))
        if not coo_r:
            return sp.csr_matrix((n_blk, n_cols), dtype=np.int64)
        mat = sp.csr_matrix(
            (np.concatenate(coo_v), (np.concatenate(coo_r), np.concatenate(coo_c))),
            shape=(n_blk, n_cols),
            dtype=np.int64,
        )
        mat.eliminate_zeros()
        mat.sum_duplicates()
        mat.sort_indices()
        return mat
    raise ValueError(f"Unsupported matrix encoding: {encoding!r}")


def _compare_two_integer_slots(
    f: h5py.File,
    slot_a: Dict[str, Any],
    slot_b: Dict[str, Any],
    chunk_rows: int = 8192,
) -> Tuple[bool, str]:
    """
    Strictly verify two non-excluded integer candidate slots match in shape, `/var` genes,
    and 100% of element values across all rows. Never selects silently on mismatch.
    """
    pa, pb = slot_a["slot_path"], slot_b["slot_path"]
    if slot_a["shape"] != slot_b["shape"]:
        return False, f"shape_mismatch({pa}={slot_a['shape']} vs {pb}={slot_b['shape']})"

    va, vb = slot_a["var_meta"], slot_b["var_meta"]
    if va["var_index"] != vb["var_index"] or va["var_syms"] != vb["var_syms"]:
        return False, f"var_gene_mismatch({pa}/{slot_a['var_group']} vs {pb}/{slot_b['var_group']})"

    if (
        slot_a["nnz_pos"] != slot_b["nnz_pos"]
        or slot_a["sum_total"] != slot_b["sum_total"]
        or slot_a["max_val"] != slot_b["max_val"]
    ):
        return (
            False,
            f"summary_mismatch({pa}: nnz={slot_a['nnz_pos']},sum={slot_a['sum_total']},max={slot_a['max_val']} "
            f"vs {pb}: nnz={slot_b['nnz_pos']},sum={slot_b['sum_total']},max={slot_b['max_val']})",
        )

    n_rows, n_cols = slot_a["shape"]
    node_a = f[pa.lstrip("/")]
    node_b = f[pb.lstrip("/")]
    for r0 in range(0, n_rows, chunk_rows):
        r1 = min(n_rows, r0 + chunk_rows)
        csr_a = _read_slot_rows_as_csr(node_a, slot_a["encoding"], n_rows, n_cols, r0, r1)
        csr_b = _read_slot_rows_as_csr(node_b, slot_b["encoding"], n_rows, n_cols, r0, r1)
        diff = (csr_a - csr_b).tocsr()
        diff.eliminate_zeros()
        if diff.nnz > 0:
            first_row = int(r0 + int(np.flatnonzero(np.diff(diff.indptr) > 0)[0]))
            return False, f"matrix_value_mismatch({pa} vs {pb}, diff_nnz={diff.nnz}, first_diff_row={first_row})"

    return True, "EXACT_MATCH"


def _inspect_obs_guide_columns(obs: h5py.Group, lib: Dict[str, Any]) -> Dict[str, Any]:
    obs_cols = _list_group_columns(obs)
    low_map = {c.lower(): c for c in obs_cols}

    id_col = next((low_map[k] for k in ("sgrna_id", "guide_id", "sgrna", "guide") if k in low_map), None)
    sym_col = next((low_map[k] for k in ("sgrna_symbol", "guide_symbol", "perturbation", "target", "gene") if k in low_map), None)
    nperts_col = next((low_map[k] for k in ("nperts", "n_guides", "num_guides", "n_sgrna") if k in low_map), None)

    valid_ids = lib["ntc_guide_ids"] | lib["excluded_safe_ids"] | set(lib["guide_to_target"].keys()) | set(lib["unmapped_guide_ids"].keys())
    valid_syms = lib["ntc_symbols"] | lib["excluded_safe_symbols"] | set(lib["symbol_to_target"].keys())

    matched_id = 0
    matched_sym = 0
    if id_col is not None:
        vals_id = _read_obs_or_var_column(obs, id_col)
        if vals_id is not None:
            matched_id = sum(1 for v in vals_id if str(v).strip() in valid_ids)
    if sym_col is not None:
        vals_sym = _read_obs_or_var_column(obs, sym_col)
        if vals_sym is not None:
            matched_sym = sum(1 for v in vals_sym if str(v).strip() in valid_syms)

    has_obs_calls = (matched_id > 0) or (matched_sym > 0)
    return {
        "has_obs_calls": has_obs_calls,
        "id_col": id_col if matched_id > 0 else None,
        "sym_col": sym_col if matched_sym > 0 else None,
        "nperts_col": nperts_col,
        "matched_id_cells": int(matched_id),
        "matched_sym_cells": int(matched_sym),
        "obs_cols": obs_cols,
    }


def _detect_raw_guide_source(
    f: h5py.File,
    var_meta: Dict[str, Any],
    obs_len: int,
    lib: Dict[str, Any],
) -> Optional[Dict[str, Any]]:
    """Detect raw guide counts in `/obsm/*` or `CRISPR Guide Capture` `/var` rows for `PREREG_F §3.2` fallback."""
    all_gids: Set[str] = lib["ntc_guide_ids"] | lib["excluded_safe_ids"] | set(lib["guide_to_target"].keys())

    ftypes = var_meta.get("var_ftypes")
    v_ids = var_meta["var_index"]
    v_syms = var_meta["var_syms"]
    crp_rows: List[int] = []
    if ftypes is not None:
        for i, ft in enumerate(ftypes):
            if "crispr" in ft.lower() or "guide" in ft.lower():
                if v_ids[i] in all_gids or v_syms[i] in all_gids:
                    crp_rows.append(i)
    if crp_rows:
        return {"source": "var_crispr_features", "crp_rows": crp_rows}

    if "obsm" in f and isinstance(f["obsm"], h5py.Group):
        for ok in sorted(f["obsm"].keys()):
            low = ok.lower()
            if not any(tok in low for tok in ("sgrna", "guide", "crispr")):
                continue
            node = f["obsm"][ok]
            g_names: Optional[List[str]] = None
            if "uns" in f and isinstance(f["uns"], h5py.Group):
                for uk in (f"{ok}_names", "sgrna_ids", "guide_ids", "crispr_guide_ids"):
                    if uk in f["uns"]:
                        g_names = _decode_str(f["uns"][uk][()])
                        break
            if isinstance(node, h5py.Dataset) and node.ndim == 2 and node.shape[0] == obs_len:
                if g_names is None and node.shape[1] == len(lib["ordered_guide_ids"]):
                    g_names = list(lib["ordered_guide_ids"])
                if g_names is not None and any(g in all_gids for g in g_names):
                    return {"source": f"/obsm/{ok}", "encoding": "dense", "guide_names": g_names}
            elif isinstance(node, h5py.Group) and {"data", "indices", "indptr"}.issubset(set(node.keys())):
                try:
                    g_enc, n_r_g, n_c_g = _detect_matrix_encoding(
                        node,
                        len(g_names) if g_names is not None else len(lib["ordered_guide_ids"]),
                        obs_len,
                    )
                except Exception:
                    continue
                if n_r_g == obs_len:
                    if g_names is None and n_c_g == len(lib["ordered_guide_ids"]):
                        g_names = list(lib["ordered_guide_ids"])
                    if g_names is not None and any(g in all_gids for g in g_names):
                        return {"source": f"/obsm/{ok}", "encoding": g_enc, "guide_names": g_names}
    return None


def check_gse337804_h5ad_feasibility(
    raw_dir: Path,
    hgnc: Optional[HGNCLookup] = None,
    lib_path: Optional[Path] = None,
) -> Dict[str, Any]:
    """
    `PREREG_G.md §4.1–§4.2` Pre-Metric Technical Feasibility Gate on `GSE337804` (`ch_1` + `ch_2`):
      1. Inspects `/raw/X`, `/X`, and all non-excluded `/layers/*` in BOTH channels, checking
         100% of stored entries for finite non-negative exact integers.
      2. If `>= 2` non-excluded slots pass in a channel, verifies identical shape, `/var` gene
         list/order, and bit-exact matrix values; if ANY mismatch exists, **STOPS immediately**
         with `status = "EXCLUDED_FEASIBILITY"` (`MULTIPLE_INCONSISTENT_INTEGER_SLOTS`).
      3. Verifies author `.obs` single-sgRNA calls (`sgrna_id`/`sgrna_symbol`) OR fallback raw
         guide count matrix (`PREREG_F §3.2`).
    """
    hgnc_obj = hgnc or build_hgnc_lookup()
    try:
        lp = lib_path or _find_sgrna_lib_path(raw_dir)
        lib = load_gse337804_sgrna_lib(lp, hgnc_obj)
        channels = discover_gse337804_h5ad_channels(raw_dir)
    except Exception as exc:
        return {"line": R14_TARGET_KEY, "feasible": False, "status": "EXCLUDED_FEASIBILITY", "reason": str(exc)}

    channel_reports: List[Dict[str, Any]] = []
    obs_modes: List[str] = []

    for ch_num, ch_tag, hp in channels:
        try:
            with h5py.File(hp, "r") as f:
                if "obs" not in f or "var" not in f:
                    return {
                        "line": R14_TARGET_KEY,
                        "feasible": False,
                        "status": "EXCLUDED_FEASIBILITY",
                        "reason": f"MISSING_ANNDATA_GROUP in {hp.name}",
                    }
                barcodes = _read_anndata_index(f["obs"], "obs")
                obs_len = len(barcodes)

                candidate_paths: List[str] = []
                excluded_layers: List[str] = []
                if "raw" in f and isinstance(f["raw"], h5py.Group) and "X" in f["raw"]:
                    candidate_paths.append("/raw/X")
                if "X" in f:
                    candidate_paths.append("/X")
                if "layers" in f and isinstance(f["layers"], h5py.Group):
                    for lk in sorted(f["layers"].keys()):
                        if _is_excluded_layer_name(lk):
                            excluded_layers.append(f"/layers/{lk}")
                        else:
                            candidate_paths.append(f"/layers/{lk}")

                slot_checks: List[Dict[str, Any]] = []
                integer_slots: List[Dict[str, Any]] = []
                for spath in candidate_paths:
                    chk = _full_slot_integer_check(f, spath, obs_len=obs_len)
                    slot_checks.append(
                        {
                            "slot_path": spath,
                            "is_integer": chk["is_integer"],
                            "reason": chk["reason"],
                            "shape": chk.get("shape"),
                            "nnz_pos": chk.get("nnz_pos"),
                            "sum_total": chk.get("sum_total"),
                        }
                    )
                    if chk["is_integer"]:
                        integer_slots.append(chk)

                if not integer_slots:
                    return {
                        "line": R14_TARGET_KEY,
                        "feasible": False,
                        "status": "EXCLUDED_FEASIBILITY",
                        "reason": f"no_integer_raw_umi_slot_in_{ch_tag}: {slot_checks}",
                        "channel": ch_tag,
                        "slot_checks": slot_checks,
                        "excluded_layers": excluded_layers,
                    }

                chosen = integer_slots[0]
                if len(integer_slots) >= 2:
                    for other in integer_slots[1:]:
                        eq_ok, eq_msg = _compare_two_integer_slots(f, chosen, other)
                        if not eq_ok:
                            return {
                                "line": R14_TARGET_KEY,
                                "feasible": False,
                                "status": "EXCLUDED_FEASIBILITY",
                                "reason": (
                                    f"MULTIPLE_INCONSISTENT_INTEGER_SLOTS in {ch_tag} ({hp.name}): {eq_msg}"
                                ),
                                "channel": ch_tag,
                                "slot_checks": slot_checks,
                                "excluded_layers": excluded_layers,
                            }
                    slot_consistency = f"MULTIPLE_IDENTICAL_INTEGER_SLOTS_VERIFIED ({[s['slot_path'] for s in integer_slots]})"
                else:
                    slot_consistency = f"UNIQUE_INTEGER_SLOT ({chosen['slot_path']})"

                obs_info = _inspect_obs_guide_columns(f["obs"], lib)
                raw_guide_info = None
                if not obs_info["has_obs_calls"]:
                    raw_guide_info = _detect_raw_guide_source(f, chosen["var_meta"], obs_len, lib)
                    if raw_guide_info is None:
                        return {
                            "line": R14_TARGET_KEY,
                            "feasible": False,
                            "status": "EXCLUDED_FEASIBILITY",
                            "reason": (
                                f"missing_obs_single_guide_calls_and_raw_guide_counts in {ch_tag}: "
                                f"obs_cols={obs_info['obs_cols']}"
                            ),
                            "channel": ch_tag,
                        }
                    obs_modes.append("PREREG_F_RAW_GUIDE_FALLBACK")
                else:
                    obs_modes.append("OBS_SINGLE_SGRNA_CALLS")

                channel_reports.append(
                    {
                        "channel": ch_tag,
                        "ch_num": ch_num,
                        "file": str(hp),
                        "obs_len": obs_len,
                        "chosen_slot": chosen["slot_path"],
                        "encoding": chosen["encoding"],
                        "shape": chosen["shape"],
                        "slot_consistency": slot_consistency,
                        "slot_checks": slot_checks,
                        "excluded_layers": excluded_layers,
                        "guide_info": obs_info if obs_info["has_obs_calls"] else raw_guide_info,
                    }
                )
        except Exception as exc:
            return {
                "line": R14_TARGET_KEY,
                "feasible": False,
                "status": "EXCLUDED_FEASIBILITY",
                "reason": f"hdf5_read_error in {hp.name}: {exc}",
            }

    if len(set(obs_modes)) != 1:
        return {
            "line": R14_TARGET_KEY,
            "feasible": False,
            "status": "EXCLUDED_FEASIBILITY",
            "reason": f"INCONSISTENT_GUIDE_MODES_ACROSS_CHANNELS: {obs_modes}",
            "channels": channel_reports,
        }

    return {
        "line": R14_TARGET_KEY,
        "feasible": True,
        "status": "FEASIBLE",
        "guide_call_mode": obs_modes[0],
        "lib_path": str(lp),
        "n_ntc_guides": len(lib["ntc_guide_ids"]),
        "n_excluded_safe_guides": len(lib["excluded_safe_ids"]),
        "n_mapped_target_guides": len(lib["guide_to_target"]),
        "channels": channel_reports,
    }


def make_slot_subset_readers(
    h5_path: Path,
    slot_path: str,
    keep_rows: np.ndarray,
    n_rows_raw: int,
    n_cols_raw: int,
    encoding: str,
) -> Tuple[
    Callable[[int, int], np.ndarray],
    Callable[[np.ndarray], np.ndarray],
    Callable[[int, int, np.ndarray], np.ndarray],
]:
    """Memory-bounded CSR/CSC/dense row-subset readers for `/raw/X`, `/X`, or `/layers/*`."""
    h5_key = slot_path.lstrip("/")
    keep_rows = np.asarray(keep_rows, dtype=np.int64)
    if keep_rows.size > 1 and not np.all(keep_rows[:-1] <= keep_rows[1:]):
        raise ValueError("keep_rows must be strictly sorted ascending")

    def _read(raw_rows: np.ndarray, feat_cols: Optional[np.ndarray] = None) -> np.ndarray:
        out_dtype = np.float32 if feat_cols is not None else np.float64
        out_w = len(feat_cols) if feat_cols is not None else n_cols_raw
        if raw_rows.size == 0:
            return np.zeros((0, out_w), dtype=out_dtype)

        with h5py.File(h5_path, "r") as f:
            node = f[h5_key]
            if encoding == "dense":
                arr = np.asarray(node[raw_rows, :], dtype=out_dtype)
                return arr[:, feat_cols] if feat_cols is not None else arr

            indptr = node["indptr"]
            indices = node["indices"]
            data = node["data"]
            if encoding == "csr":
                r_min, r_max = int(raw_rows[0]), int(raw_rows[-1])
                span = r_max - r_min + 1
                if len(raw_rows) / max(span, 1) >= 0.20:
                    ptrs = np.asarray(indptr[r_min : r_max + 2], dtype=np.int64)
                    p0, p1 = int(ptrs[0]), int(ptrs[-1])
                    mat = sp.csr_matrix(
                        (data[p0:p1], indices[p0:p1], ptrs - np.int64(p0)),
                        shape=(span, n_cols_raw),
                    )[raw_rows - r_min, :]
                else:
                    all_d, all_i = [], []
                    new_ptr = np.zeros(len(raw_rows) + 1, dtype=np.int64)
                    cur_p = 0
                    starts = indptr[raw_rows]
                    ends = indptr[raw_rows + 1]
                    for idx_r, (s0, s1) in enumerate(zip(starts, ends)):
                        if s1 > s0:
                            all_d.append(data[int(s0) : int(s1)])
                            all_i.append(indices[int(s0) : int(s1)])
                            cur_p += int(s1 - s0)
                        new_ptr[idx_r + 1] = cur_p
                    d_cat = np.concatenate(all_d) if all_d else np.empty(0, dtype=np.float32)
                    i_cat = np.concatenate(all_i) if all_i else np.empty(0, dtype=np.int32)
                    mat = sp.csr_matrix((d_cat, i_cat, new_ptr), shape=(len(raw_rows), n_cols_raw))
                if feat_cols is not None:
                    mat = mat[:, feat_cols]
                return mat.toarray().astype(out_dtype, copy=False)

            if encoding == "csc":
                col_targets = (
                    np.arange(n_cols_raw, dtype=np.int64)
                    if feat_cols is None
                    else np.asarray(feat_cols, dtype=np.int64)
                )
                pos_map = np.full(n_rows_raw, -1, dtype=np.int64)
                pos_map[raw_rows] = np.arange(len(raw_rows), dtype=np.int64)
                coo_r, coo_c, coo_v = [], [], []
                all_ptrs = np.asarray(indptr[:], dtype=np.int64)
                for j_out, c_idx in enumerate(col_targets):
                    p0, p1 = int(all_ptrs[c_idx]), int(all_ptrs[c_idx + 1])
                    if p1 <= p0:
                        continue
                    rows_c = indices[p0:p1]
                    mapped = pos_map[rows_c]
                    hit = mapped >= 0
                    if np.any(hit):
                        coo_r.append(mapped[hit])
                        coo_c.append(np.full(int(np.sum(hit)), j_out, dtype=np.int32))
                        coo_v.append(data[p0:p1][hit])
                out_arr = np.zeros((len(raw_rows), len(col_targets)), dtype=out_dtype)
                if coo_r:
                    out_arr[np.concatenate(coo_r), np.concatenate(coo_c)] = np.concatenate(coo_v).astype(out_dtype)
                return out_arr

            raise ValueError(f"Unsupported HDF5 matrix encoding: {encoding!r}")

    return (
        (lambda s, e: _read(keep_rows[s:e], None)),
        (lambda loc_idx: _read(keep_rows[loc_idx], None)),
        (lambda s, e, fcols: _read(keep_rows[s:e], fcols)),
    )


def _classify_obs_cell(
    gid_raw: str,
    sym_raw: str,
    nperts_val: float,
    lib: Dict[str, Any],
) -> Tuple[str, str]:
    """
    Classify one cell from `.obs` single-sgRNA metadata strictly per `PREREG_G §4.2–§4.4`:
      Returns `(role, key)`:
        - `("EXCLUDED_MULTIGUIDE", "")`
        - `("EXCLUDED_UNASSIGNED", "")`
        - `("EXCLUDED_SAFE_TARGETING", "")` (`ctrl_T` / `gene_desert_*`)
        - `("EXCLUDED_UNKNOWN_OR_UNMAPPED", "")`
        - `("NTC", "")` (`ctrl_NT` / `non-targeting*`)
        - `("TARGET", canon_sym)`
    """
    if np.isfinite(nperts_val) and nperts_val > 1.0 + 1e-6:
        return "EXCLUDED_MULTIGUIDE", ""

    g_str = gid_raw.strip()
    s_str = sym_raw.strip()
    g_low = g_str.lower()
    s_low = s_str.lower()

    if any(tok in (g_low, s_low) for tok in ("multiplet", "doublet")):
        return "EXCLUDED_MULTIGUIDE", ""
    if _MULTI_GUIDE_DELIM_RE.search(g_str) or _MULTI_GUIDE_DELIM_RE.search(s_str):
        return "EXCLUDED_MULTIGUIDE", ""
    if (g_low in _BAD_GUIDE_TOKENS) and (s_low in _BAD_GUIDE_TOKENS):
        return "EXCLUDED_UNASSIGNED", ""
    if g_low == "unassigned" or s_low == "unassigned":
        return "EXCLUDED_UNASSIGNED", ""

    if g_str in lib["excluded_safe_ids"] or s_str in lib["excluded_safe_symbols"] or "gene_desert" in g_low or "gene_desert" in s_low:
        return "EXCLUDED_SAFE_TARGETING", ""

    if g_str in lib["ntc_guide_ids"] or (not g_str and s_str in lib["ntc_symbols"]):
        return "NTC", ""
    if g_str and g_str in lib["unmapped_guide_ids"]:
        return "EXCLUDED_UNKNOWN_OR_UNMAPPED", ""
    if g_str and g_str in lib["guide_to_target"]:
        canon_g = lib["guide_to_target"][g_str][1]
        if s_str and s_str in lib["symbol_to_target"] and lib["symbol_to_target"][s_str][1] != canon_g:
            return "EXCLUDED_UNKNOWN_OR_UNMAPPED", ""
        return "TARGET", canon_g

    if not g_str and s_str in lib["symbol_to_target"]:
        return "TARGET", lib["symbol_to_target"][s_str][1]

    if g_str in lib["ntc_symbols"] or s_str in lib["ntc_symbols"]:
        if not g_str or g_str in lib["ntc_guide_ids"]:
            return "NTC", ""

    return "EXCLUDED_UNKNOWN_OR_UNMAPPED", ""


def _read_guide_csc_for_qc_cells(
    h5_path: Path,
    slot_path: str,
    encoding: str,
    n_rows_raw: int,
    n_cols_raw: int,
    guide_info: Dict[str, Any],
    qc_cells: np.ndarray,
    init_cells: np.ndarray,
    lib: Dict[str, Any],
) -> Tuple[sp.csc_matrix, np.ndarray]:
    """Read raw guide counts `(n_guides x len(qc_cells))` and fill Poisson-Gaussian histogram on `init_cells`."""
    ordered_g = lib["ordered_guide_ids"]
    gid_to_pos = {g: i for i, g in enumerate(ordered_g)}
    hist = np.zeros((len(ordered_g), 66), dtype=np.int64)

    all_cells = np.unique(np.concatenate([init_cells, qc_cells]))
    qc_map = np.full(n_rows_raw, -1, dtype=np.int64)
    qc_map[qc_cells] = np.arange(len(qc_cells), dtype=np.int64)
    is_init = np.zeros(n_rows_raw, dtype=bool)
    is_init[init_cells] = True

    coo_r: List[np.ndarray] = []
    coo_c: List[np.ndarray] = []
    coo_v: List[np.ndarray] = []

    if guide_info["source"] == "var_crispr_features":
        crp_rows = np.asarray(guide_info["crp_rows"], dtype=np.int64)
        with h5py.File(h5_path, "r") as f:
            var_key = "raw/var" if slot_path == "/raw/X" else "var"
            vm = _extract_var_metadata(f[var_key], var_key)
        col_gpos = []
        valid_crp = []
        for col_i in crp_rows:
            gp = gid_to_pos.get(vm["var_index"][col_i], gid_to_pos.get(vm["var_syms"][col_i], -1))
            if gp >= 0:
                valid_crp.append(int(col_i))
                col_gpos.append(gp)
        valid_crp_arr = np.asarray(valid_crp, dtype=np.int64)
        _, _, rcc = make_slot_subset_readers(h5_path, slot_path, all_cells, n_rows_raw, n_cols_raw, encoding)
        for s0 in range(0, len(all_cells), 8192):
            e0 = min(len(all_cells), s0 + 8192)
            blk = rcc(s0, e0, valid_crp_arr)
            rr_loc, cc_loc = np.nonzero(blk > 0)
            if rr_loc.size == 0:
                continue
            vals = blk[rr_loc, cc_loc].astype(np.int64, copy=False)
            raw_cell_idx = all_cells[s0 + rr_loc]
            g_pos = np.asarray(col_gpos, dtype=np.int32)[cc_loc]
            in_m = is_init[raw_cell_idx]
            if np.any(in_m):
                np.add.at(hist, (g_pos[in_m], np.clip(vals[in_m], 1, 65)), 1)
            q_loc = qc_map[raw_cell_idx]
            qm = q_loc >= 0
            if np.any(qm):
                coo_r.append(g_pos[qm])
                coo_c.append(q_loc[qm].astype(np.int32))
                coo_v.append(vals[qm].astype(np.float64))
    else:
        src_key = str(guide_info["source"]).lstrip("/")
        g_names = guide_info["guide_names"]
        col_gpos = np.asarray([gid_to_pos.get(g, -1) for g in g_names], dtype=np.int32)
        keep_gcols = np.flatnonzero(col_gpos >= 0)
        g_enc = guide_info["encoding"]
        _, _, rcc = make_slot_subset_readers(h5_path, src_key, all_cells, n_rows_raw, len(g_names), g_enc)
        for s0 in range(0, len(all_cells), 8192):
            e0 = min(len(all_cells), s0 + 8192)
            blk = rcc(s0, e0, keep_gcols)
            rr_loc, cc_loc = np.nonzero(blk > 0)
            if rr_loc.size == 0:
                continue
            vals = blk[rr_loc, cc_loc].astype(np.int64, copy=False)
            raw_cell_idx = all_cells[s0 + rr_loc]
            g_pos = col_gpos[keep_gcols[cc_loc]]
            in_m = is_init[raw_cell_idx]
            if np.any(in_m):
                np.add.at(hist, (g_pos[in_m], np.clip(vals[in_m], 1, 65)), 1)
            q_loc = qc_map[raw_cell_idx]
            qm = q_loc >= 0
            if np.any(qm):
                coo_r.append(g_pos[qm])
                coo_c.append(q_loc[qm].astype(np.int32))
                coo_v.append(vals[qm].astype(np.float64))

    if coo_r:
        csc = sp.csc_matrix(
            (np.concatenate(coo_v), (np.concatenate(coo_r), np.concatenate(coo_c))),
            shape=(len(ordered_g), len(qc_cells)),
        )
    else:
        csc = sp.csc_matrix((len(ordered_g), len(qc_cells)), dtype=np.float64)
    return csc, hist


def build_a549_myllia_audit_parts(
    raw_dir: Path,
    hgnc: HGNCLookup,
    lib_path: Optional[Path] = None,
    min_umi: float = 300.0,
    min_genes: int = 200,
    min_batch_nt: int = 15,
) -> Tuple[List[PartSpec], Dict[str, Any]]:
    """
    Independent `GSE337804` (`a549_myllia`) raw `.h5ad` reader (`PREREG_G.md §4–§5`):
      - Runs full pre-metric feasibility + 100%-entry integer raw-count check + multi-slot guard.
      - Removes non-human / Cas9 / CRISPR guide rows before computing cell `umi` and `n_genes`.
      - Enforces `n_umis >= 300`, `n_genes >= 200`, **`max_pct_mt = None` (`DepleteX`)**.
      - Strictly keeps ONLY `ctrl_NT` (`non-targeting_*`) as `NTC`; excludes ALL `ctrl_T` (`gene_desert_*`).
      - Pools all QC-valid single-sgRNA cells per canonical HGNC gene (`pert_key = canon_sym`).
      - Assigns disjoint `Half A (1) = ch_1`, `Half B (2) = ch_2`, `batches = ["ch_1", "ch_2"]`.
    """
    feas = check_gse337804_h5ad_feasibility(raw_dir=raw_dir, hgnc=hgnc, lib_path=lib_path)
    if not feas["feasible"]:
        raise RuntimeError(f"EXCLUDED_FEASIBILITY: {feas['reason']}")

    lib = load_gse337804_sgrna_lib(Path(feas["lib_path"]), hgnc)
    guide_mode = feas["guide_call_mode"]
    ordered_guides = lib["ordered_guide_ids"]

    ch_prepared: List[Dict[str, Any]] = []
    hist_total = np.zeros((len(ordered_guides), 66), dtype=np.int64)

    for ch_rep in feas["channels"]:
        ch_num = int(ch_rep["ch_num"])
        ch_tag = str(ch_rep["channel"])
        hp = Path(ch_rep["file"])
        slot_path = str(ch_rep["chosen_slot"])
        encoding = str(ch_rep["encoding"])
        n_rows_raw, n_cols_raw = int(ch_rep["shape"][0]), int(ch_rep["shape"][1])

        with h5py.File(hp, "r") as f:
            var_key = "raw/var" if slot_path == "/raw/X" else "var"
            vm = _extract_var_metadata(f[var_key], var_key)
            v_syms = vm["var_syms"]
            v_ens = vm["var_ens"]
            v_ftypes = vm["var_ftypes"]
            v_ids = vm["var_index"]

            keep_gene_mask = np.ones(n_cols_raw, dtype=bool)
            for j in range(n_cols_raw):
                s_up = v_syms[j].upper()
                e_up = v_ens[j].upper()
                if ("CAS9" in s_up) or ("CAS9" in e_up) or e_up.startswith("ENSMUSG") or v_syms[j].lower().startswith("mm10"):
                    keep_gene_mask[j] = False
                    continue
                if v_ftypes is not None and ("crispr" in v_ftypes[j].lower() or "guide" in v_ftypes[j].lower()):
                    keep_gene_mask[j] = False
                    continue
                if guide_mode == "PREREG_F_RAW_GUIDE_FALLBACK" and (
                    v_ids[j] in lib["ntc_guide_ids"]
                    or v_ids[j] in lib["excluded_safe_ids"]
                    or v_ids[j] in lib["guide_to_target"]
                ):
                    keep_gene_mask[j] = False

            keep_gene_cols = np.flatnonzero(keep_gene_mask)
            if keep_gene_cols.size == 0:
                raise RuntimeError(f"Zero human RNA features remaining in {ch_tag} ({hp.name})")

            g_info = ch_rep["guide_info"]
            obs_gid = (
                _read_obs_or_var_column(f["obs"], g_info["id_col"])
                if (guide_mode == "OBS_SINGLE_SGRNA_CALLS" and g_info.get("id_col"))
                else None
            )
            obs_sym = (
                _read_obs_or_var_column(f["obs"], g_info["sym_col"])
                if (guide_mode == "OBS_SINGLE_SGRNA_CALLS" and g_info.get("sym_col"))
                else None
            )
            obs_nperts = (
                np.asarray(pd.to_numeric(_read_obs_or_var_column(f["obs"], g_info["nperts_col"]), errors="coerce"), dtype=np.float64)
                if (guide_mode == "OBS_SINGLE_SGRNA_CALLS" and g_info.get("nperts_col"))
                else None
            )

        all_rows = np.arange(n_rows_raw, dtype=np.int64)
        rc_all, _, _ = make_slot_subset_readers(hp, slot_path, all_rows, n_rows_raw, n_cols_raw, encoding)
        umi_raw = np.zeros(n_rows_raw, dtype=np.float64)
        ngenes_raw = np.zeros(n_rows_raw, dtype=np.int64)
        for s0 in range(0, n_rows_raw, 16384):
            e0 = min(n_rows_raw, s0 + 16384)
            blk = rc_all(s0, e0)[:, keep_gene_cols]
            umi_raw[s0:e0] = blk.sum(axis=1)
            ngenes_raw[s0:e0] = (blk > 0).sum(axis=1)

        qc_ok = (umi_raw >= float(min_umi)) & (ngenes_raw >= int(min_genes))
        qc_cells = np.flatnonzero(qc_ok)
        init_cells = np.flatnonzero(ngenes_raw >= int(min_genes))

        c_qc = None
        if guide_mode == "PREREG_F_RAW_GUIDE_FALLBACK":
            c_qc, ch_hist = _read_guide_csc_for_qc_cells(
                h5_path=hp,
                slot_path=slot_path,
                encoding=encoding,
                n_rows_raw=n_rows_raw,
                n_cols_raw=n_cols_raw,
                guide_info=g_info,
                qc_cells=qc_cells,
                init_cells=init_cells,
                lib=lib,
            )
            hist_total += ch_hist

        ch_prepared.append(
            {
                "ch_num": ch_num,
                "ch_tag": ch_tag,
                "file": hp,
                "slot_path": slot_path,
                "encoding": encoding,
                "n_rows_raw": n_rows_raw,
                "n_cols_raw": n_cols_raw,
                "keep_gene_cols": keep_gene_cols,
                "v_ens": [v_ens[j] for j in keep_gene_cols],
                "v_sym": [v_syms[j] for j in keep_gene_cols],
                "umi_raw": umi_raw,
                "qc_cells": qc_cells,
                "n_failed_qc": int(n_rows_raw - len(qc_cells)),
                "obs_gid": obs_gid,
                "obs_sym": obs_sym,
                "obs_nperts": obs_nperts,
                "c_qc": c_qc,
            }
        )

    thrs: Optional[np.ndarray] = None
    if guide_mode == "PREREG_F_RAW_GUIDE_FALLBACK":
        thrs = fit_poisson_gaussian_thresholds_from_hist(hist_total, min_thr=5, max_thr=20, default_thr=7)

    counts_audit: Dict[str, int] = defaultdict(int)
    per_ch_records: List[Dict[str, Any]] = []

    for ch in ch_prepared:
        qc_cells = ch["qc_cells"]
        counts_audit["n_raw_cells"] += int(ch["n_rows_raw"])
        counts_audit["n_failed_qc_umi_or_genes"] += int(ch["n_failed_qc"])

        keep_rows_ch: List[int] = []
        nt_mask_ch: List[bool] = []
        pkeys_ch: List[str] = []

        if guide_mode == "OBS_SINGLE_SGRNA_CALLS":
            obs_gid = ch["obs_gid"]
            obs_sym = ch["obs_sym"]
            obs_nperts = ch["obs_nperts"]
            for raw_r in qc_cells:
                g_val = str(obs_gid[raw_r]) if obs_gid is not None else ""
                s_val = str(obs_sym[raw_r]) if obs_sym is not None else ""
                np_val = float(obs_nperts[raw_r]) if obs_nperts is not None else float("nan")
                role, key = _classify_obs_cell(g_val, s_val, np_val, lib)
                if role == "NTC":
                    keep_rows_ch.append(int(raw_r))
                    nt_mask_ch.append(True)
                    pkeys_ch.append("")
                    counts_audit["n_kept_ntc_ctrl_NT"] += 1
                elif role == "TARGET":
                    keep_rows_ch.append(int(raw_r))
                    nt_mask_ch.append(False)
                    pkeys_ch.append(key)
                    counts_audit["n_kept_targeting_cells"] += 1
                else:
                    counts_audit[f"n_{role.lower()}"] += 1
        else:
            c_qc = ch["c_qc"]
            assert c_qc is not None and thrs is not None
            t1, t2, a1 = extract_top2_csc(c_qc)
            pos_mask = c_qc.data >= thrs[c_qc.indices]
            pos_csc = sp.csc_matrix((pos_mask.astype(np.int8), c_qc.indices, c_qc.indptr), shape=c_qc.shape)
            n_pos = np.asarray(pos_csc.sum(axis=0)).ravel()
            s1_mask = (n_pos == 1) & (t1 >= 5.0) & ((t2 < 3.0) | (t1 / np.maximum(t2, 0.5) >= 3.0))
            for loc_j in range(len(qc_cells)):
                if not s1_mask[loc_j]:
                    counts_audit["n_excluded_multiguide_or_unassigned_fallback"] += 1
                    continue
                gid = ordered_guides[int(a1[loc_j])]
                raw_r = int(qc_cells[loc_j])
                if gid in lib["excluded_safe_ids"]:
                    counts_audit["n_excluded_safe_targeting"] += 1
                elif gid in lib["ntc_guide_ids"]:
                    keep_rows_ch.append(raw_r)
                    nt_mask_ch.append(True)
                    pkeys_ch.append("")
                    counts_audit["n_kept_ntc_ctrl_NT"] += 1
                elif gid in lib["guide_to_target"]:
                    keep_rows_ch.append(raw_r)
                    nt_mask_ch.append(False)
                    pkeys_ch.append(lib["guide_to_target"][gid][1])
                    counts_audit["n_kept_targeting_cells"] += 1
                else:
                    counts_audit["n_excluded_unknown_or_unmapped"] += 1

        keep_arr = np.asarray(keep_rows_ch, dtype=np.int64)
        nt_arr = np.asarray(nt_mask_ch, dtype=bool)
        if keep_arr.size == 0 or int(np.sum(nt_arr)) == 0:
            raise RuntimeError(f"Channel {ch['ch_tag']} yielded zero kept/NTC cells after QC & guide filtering")

        per_ch_records.append(
            {
                **ch,
                "keep_arr": keep_arr,
                "nt_arr": nt_arr,
                "pkeys": pkeys_ch,
            }
        )

    # Apply `_apply_min_ntc_batch_fallback` across channels (`ch_1` / `ch_2`) safely
    all_batches_raw = np.concatenate([np.full(len(r["keep_arr"]), r["ch_tag"], dtype=object) for r in per_ch_records])
    all_nt = np.concatenate([r["nt_arr"] for r in per_ch_records])
    all_halves = np.concatenate(
        [np.full(len(r["keep_arr"]), np.int8(1 if r["ch_num"] == 1 else 2), dtype=np.int8) for r in per_ch_records]
    )
    batches_adj, n_nt_fallback = _apply_min_ntc_batch_fallback(
        all_batches_raw,
        all_nt,
        all_halves,
        group_labels=None,
        min_batch_nt=min_batch_nt,
    )

    # Build global `pert_meta` pooled strictly by canonical HGNC symbol (`pert_key = canon_sym`)
    canon_to_meta: Dict[str, Tuple[str, str, str, str]] = {}
    for _, meta in lib["guide_to_target"].items():
        canon_to_meta.setdefault(meta[1], meta)

    parts: List[PartSpec] = []
    offset = 0
    unique_kept_targets: Set[str] = set()
    channel_summary: Dict[str, Dict[str, int]] = {}

    for rec in per_ch_records:
        n_k = len(rec["keep_arr"])
        ch_batches = batches_adj[offset : offset + n_k]
        offset += n_k
        ch_halves = np.full(n_k, np.int8(1 if rec["ch_num"] == 1 else 2), dtype=np.int8)

        pmeta_ch: Dict[str, Tuple[str, str, str, str]] = {}
        for pk in rec["pkeys"]:
            if pk:
                unique_kept_targets.add(pk)
                pmeta_ch[pk] = canon_to_meta[pk]

        keep_gcols = rec["keep_gene_cols"]
        raw_rc, raw_rr, raw_rcc = make_slot_subset_readers(
            h5_path=rec["file"],
            slot_path=rec["slot_path"],
            keep_rows=rec["keep_arr"],
            n_rows_raw=rec["n_rows_raw"],
            n_cols_raw=rec["n_cols_raw"],
            encoding=rec["encoding"],
        )

        # Bind slice loops cleanly via default args inside lambdas
        def _make_fns(
            _rc: Callable[[int, int], np.ndarray] = raw_rc,
            _rr: Callable[[np.ndarray], np.ndarray] = raw_rr,
            _rcc: Callable[[int, int, np.ndarray], np.ndarray] = raw_rcc,
            _gcols: np.ndarray = keep_gcols,
        ) -> Tuple[
            Callable[[int, int], np.ndarray],
            Callable[[np.ndarray], np.ndarray],
            Callable[[int, int, np.ndarray], np.ndarray],
        ]:
            return (
                (lambda s, e: _rc(s, e)[:, _gcols]),
                (lambda loc: _rr(loc)[:, _gcols]),
                (lambda s, e, fcols: _rcc(s, e, _gcols[fcols])),
            )

        fn_rc, fn_rr, fn_rcc = _make_fns()
        parts.append(
            PartSpec(
                tag=f"GSE337804_A549_{rec['ch_tag']}",
                path=rec["file"],
                n_cells=n_k,
                n_genes_raw=len(keep_gcols),
                nt_mask=rec["nt_arr"],
                batches=ch_batches.astype(str),
                half_flags=ch_halves,
                umi=rec["umi_raw"][rec["keep_arr"]],
                pert_keys=rec["pkeys"],
                pert_meta=pmeta_ch,
                var_gene_ids=rec["v_ens"],
                var_gene_names=rec["v_sym"],
                read_chunk_fn=fn_rc,
                read_rows_fn=fn_rr,
                read_chunk_cols_fn=fn_rcc,
            )
        )
        channel_summary[rec["ch_tag"]] = {
            "half": int(1 if rec["ch_num"] == 1 else 2),
            "n_kept": int(n_k),
            "n_ntc": int(np.sum(rec["nt_arr"])),
            "n_targeting": int(n_k - np.sum(rec["nt_arr"])),
        }

    audit_meta = {
        "line": R14_TARGET_KEY,
        "feasibility": feas,
        "guide_call_mode": guide_mode,
        "pooling_unit": "canonical_hgnc_gene_symbol",
        "max_pct_mt": None,
        "min_umi": float(min_umi),
        "min_genes": int(min_genes),
        "split_rule": "ch_1->HalfA(1),ch_2->HalfB(2)",
        "n_nt_fallback_cells": int(n_nt_fallback),
        "n_unique_targets_kept": int(len(unique_kept_targets)),
        "counts": dict(counts_audit),
        "channels": channel_summary,
    }
    return parts, audit_meta


def standardize_a549_myllia_to_audit_line(
    raw_dir: Path,
    cache_out_root: Path,
    aligned_genes: Optional[Sequence[str]] = None,
    hgnc_tsv: Path = DEFAULT_HGNC_TSV,
    lib_path: Optional[Path] = None,
    min_gene_overlap: int = 100,
    n_workers: int = 2,
    phase: str = "reader_xcheck_r14",
    status_dir: Path = R14_STATUS_DIR,
) -> Tuple[AuditLineData, Dict[str, Any]]:
    """
    End-to-end clean-room ingestion + Addendum-B standardization for `a549_myllia`
    (`GSE337804`), enforcing `enforce_r14_gate(phase, ...)` BEFORE touching `.h5ad` matrices.
    """
    enforce_r14_gate(
        phase=phase,
        target_lines=(R14_TARGET_KEY,),
        status_dir=status_dir,
        extra_dirs=(raw_dir, cache_out_root),
        require_confirmatory_file=False,
    )
    hgnc_prep = build_hgnc_lookup(hgnc_tsv)
    parts, meta = build_a549_myllia_audit_parts(raw_dir=raw_dir, hgnc=hgnc_prep, lib_path=lib_path)
    line_dir = process_r10_line(
        line_name=R14_TARGET_KEY,
        parts=parts,
        hgnc=hgnc_prep,
        chunk_size=8192,
        out_root=cache_out_root,
        extra_meta={"audit_reader": meta},
        n_workers=n_workers,
    )
    hgnc_canon = CanonicalHGNC.from_tsv(hgnc_tsv)
    adata = load_line_dataset(
        line_dir=line_dir,
        representation="log1p",
        aligned_genes=aligned_genes,
        hgnc=hgnc_canon,
        min_gene_overlap=min_gene_overlap,
    )
    return adata, meta


def resolve_r14_line_cache_dir(line_name: str, phase: str) -> Path:
    """Resolve standardized cache directory cleanly under the current Round 14 blinding phase."""
    if phase == "synthetic_test" or phase.startswith("dev_"):
        if line_name in R9_LINE_SET:
            return R9_LINES_DIR / line_name
        if line_name in R10_ALL_LINES:
            return R10_LINES_DIR / line_name
        if line_name in R11_ALL_LINES:
            return R11_LINES_DIR / line_name
        if line_name in R12_ALL_LINES:
            return R12_LINES_DIR / line_name
        if line_name in R13_ALL_LINES:
            return R13_LINES_DIR / line_name
        raise RuntimeError(f"Phase {phase} cannot resolve Round 14 line {line_name!r}.")
    if line_name in R9_LINE_SET:
        return R9_LINES_DIR / line_name
    if line_name in R10_SOURCE_EXT_NAMES:
        return R10_LINES_DIR / line_name
    if line_name in R14_TARGET_SET:
        return R14_LINES_DIR / line_name
    raise ValueError(f"Unknown line for Round 14 holdout evaluation: {line_name!r}")
