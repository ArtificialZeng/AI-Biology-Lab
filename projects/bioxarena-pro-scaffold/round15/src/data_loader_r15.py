#!/usr/bin/env python3
"""
Round 15 (`PREREG_H` v2) Data Loader, Line Validator, & Hard Blinding Gate.

Strictly reuses frozen `round10` / `round11` / `round12` primitives without modifying any
file under `round9/`–`round14/`:
  - Enforces hard blinding gate (`check_prereg_h_gate`): refuses `--phase holdout_r15`
    unless BOTH `PREREG_H_COMMIT.txt` and `CODE_FREEZE_R15.txt` exist and are non-empty.
  - Locks candidate lines to `("cal33", "hela", "ac16", "molm13", "u87mg")` and strictly
    excludes `"jurkat"` (`EXCLUDED_R15_LINES`) under `PREREG_H §3`.
  - Exposes `DEFAULT_R15_SEED = 20261231`, per-wave `CONFIRMATORY_LINES_R15.txt` parser
    (`0 <= K_conf <= 5`), `SA-H1` (`leave-vcc2025_h1-out`) cohort filter, `SA-Panel`
    (`129` shared val-test overlap) loader, and Arc official context-map loader (`§3.1`).
"""

from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

import pandas as pd

ROUND15_SRC_DIR = Path(__file__).resolve().parent
ROUND15_ROOT = ROUND15_SRC_DIR.parent
REPO_ROOT = ROUND15_ROOT.parent
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
R10_SCRIPTS_DIR = REPO_ROOT / "round10" / "shared" / "scripts"
R11_SRC_DIR = REPO_ROOT / "round11" / "src"
R12_SRC_DIR = REPO_ROOT / "round12" / "src"
for _p in (str(R10_SRC_DIR), str(R10_SCRIPTS_DIR), str(R11_SRC_DIR), str(R12_SRC_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from data_loader_r12 import (  # noqa: E402
    DEFAULT_HGNC_FILE,
    HGNCBridge,
    LineDatasetR10,
    R9_COMMON_GENES_FILE,
    R9_LINES_DIR,
    R10_LINES_DIR,
    R11_LINES_DIR,
    list_standardized_lines,
    load_r9_common_genes,
    load_source_cohort_r12,
    load_standardized_line_r10,
    read_gate_hash,
)
from prepare_r10_lines import HGNCLookup, build_hgnc_lookup  # noqa: E402

DEFAULT_R15_SEED = 20261231
ALLOWED_R15_LINES: Tuple[str, ...] = ("cal33", "hela", "ac16", "molm13", "u87mg")
ALLOWED_R15_LINES_SET = frozenset(ALLOWED_R15_LINES)
EXCLUDED_R15_LINES = frozenset({"jurkat"})
SA_H1_EXCLUDED_SOURCE = "vcc2025_h1"

R15_PREREG_DIR = ROUND15_ROOT / "prereg"
R15_SHARED_DIR = ROUND15_ROOT / "shared"
R15_RAW_DIR = R15_SHARED_DIR / "raw" / "vcc2026"
R15_LINES_DIR = R15_SHARED_DIR / "lines"
R15_STATUS_DIR = ROUND15_ROOT / "status"
R15_RESULTS_DIR = ROUND15_ROOT / "results"

PREREG_H_COMMIT_FILE = R15_STATUS_DIR / "PREREG_H_COMMIT.txt"
CODE_FREEZE_R15_FILE = R15_STATUS_DIR / "CODE_FREEZE_R15.txt"
R15_CONFIRMATORY_LINES_FILE = R15_STATUS_DIR / "CONFIRMATORY_LINES_R15.txt"
DEFAULT_PANEL_OVERLAP_129_FILE = R15_SHARED_DIR / "vcc2026_val_test_overlap_129.txt"
DEFAULT_ARC_CONTEXT_MAP_FILE = R15_SHARED_DIR / "vcc2026_context_to_line.json"
DEFAULT_LINE_MANIFEST_FILE = R15_SHARED_DIR / "vcc2026_line_manifest.json"

_CANON_LINE_TOKENS: Dict[str, str] = {
    "cal33": "cal33",
    "hela": "hela",
    "ac16": "ac16",
    "molm13": "molm13",
    "u87mg": "u87mg",
    "jurkat": "jurkat",
}


def normalize_vcc_line_token(raw_label: str) -> Optional[str]:
    """
    Normalize explicit cell-line string case-insensitively after stripping `-`, `_`,
    and whitespace (`PREREG_H §3.1`), returning one of the 6 known VCC tokens or `None`.
    """
    if raw_label is None:
        return None
    s = str(raw_label).strip()
    if not s or s.lower() in {"nan", "none", "na", "null", "unknown", "unassigned", "multiplet", "doublet"}:
        return None
    compact = re.sub(r"[-_\s]+", "", s.lower())
    return _CANON_LINE_TOKENS.get(compact)


def validate_r15_candidate_line(raw_label: str) -> str:
    """Validate and canonicalize a candidate target line; hard-reject `jurkat` or unknown names."""
    canon = normalize_vcc_line_token(raw_label)
    if canon in EXCLUDED_R15_LINES:
        raise ValueError(
            f"EXCLUDED_LINE: '{raw_label}' normalizes to '{canon}', which is strictly excluded "
            "under PREREG_H §3 (homologous to nadig_jurkat in r9_6)."
        )
    if canon not in ALLOWED_R15_LINES_SET:
        raise ValueError(
            f"Unsupported Round 15 line {raw_label!r}; allowed candidate lines are {ALLOWED_R15_LINES}."
        )
    return str(canon)


def check_prereg_h_gate(gate_dir: Optional[Path] = None) -> Tuple[str, str]:
    """
    Hard blinding lock (`PREREG_H`): refuse Round 15 held-out evaluation unless BOTH
    `PREREG_H_COMMIT.txt` and `CODE_FREEZE_R15.txt` exist and contain valid non-empty tokens.

    If `gate_dir` is omitted, checks `R15_STATUS_DIR` (falling back to `R15_PREREG_DIR`
    for `PREREG_H_COMMIT.txt` or `ROUND15_SRC_DIR` for `CODE_FREEZE_R15.txt` if present).
    """
    if gate_dir is not None:
        gd = Path(gate_dir)
        prereg_file = gd / "PREREG_H_COMMIT.txt"
        freeze_file = gd / "CODE_FREEZE_R15.txt"
    else:
        prereg_file = PREREG_H_COMMIT_FILE if PREREG_H_COMMIT_FILE.exists() else (R15_PREREG_DIR / "PREREG_H_COMMIT.txt")
        freeze_file = CODE_FREEZE_R15_FILE if CODE_FREEZE_R15_FILE.exists() else (ROUND15_SRC_DIR / "CODE_FREEZE_R15.txt")

    prereg_hash = read_gate_hash(prereg_file)
    freeze_hash = read_gate_hash(freeze_file)
    if not prereg_hash:
        raise RuntimeError(
            f"PREREG_H hard gate blocked: missing/empty {prereg_file}. "
            "Round 15 held-out lines cannot be evaluated before PREREG_H lock."
        )
    if not freeze_hash:
        raise RuntimeError(
            f"CODE_FREEZE_R15 hard gate blocked: missing/empty {freeze_file}. "
            "Round 15 held-out lines cannot be evaluated before Round 15 code freeze."
        )
    return prereg_hash, freeze_hash


def parse_confirmatory_lines_r15(
    conf_file: Optional[Path] = None,
    raw_lines: Optional[Sequence[str]] = None,
) -> List[str]:
    """
    Parse confirmatory target lines for a wave (`0 <= K_conf <= 5`).

    Supports `"NONE"` or empty file (`K_conf = 0`), hard-rejects `"jurkat"` (`EXCLUDED_LINE`),
    and enforces that every listed line belongs to `ALLOWED_R15_LINES`.
    """
    items: List[str] = []
    if raw_lines is not None:
        for entry in raw_lines:
            for tok in str(entry).split(","):
                s = tok.strip()
                if s and not s.startswith("#"):
                    items.append(s)
    else:
        cf = Path(conf_file) if conf_file is not None else R15_CONFIRMATORY_LINES_FILE
        if not cf.exists():
            raise FileNotFoundError(
                f"Missing Round 15 confirmatory lines manifest: {cf}. "
                "Generate counts-only CONFIRMATORY_LINES_R15.txt before confirmatory adjudication."
            )
        for line in cf.read_text(encoding="utf-8").splitlines():
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            for tok in s.split(","):
                t = tok.strip()
                if t:
                    items.append(t)

    if len(items) == 1 and items[0].upper() in {"NONE", "EMPTY", "K0"}:
        return []
    out: List[str] = []
    seen = set()
    for raw in items:
        if raw.upper() in {"NONE", "EMPTY", "K0"}:
            continue
        canon = validate_r15_candidate_line(raw)
        if canon not in seen:
            seen.add(canon)
            out.append(canon)
    if len(out) > len(ALLOWED_R15_LINES):
        raise ValueError(f"Round 15 allows at most {len(ALLOWED_R15_LINES)} confirmatory lines; got {out}")
    return out


def load_source_cohort_r15(
    source_mode: str,
    representation: str,
    common_genes: Sequence[str],
    hgnc_bridge: Optional[HGNCBridge],
    exclude_sources: Iterable[str] = (),
) -> Dict[str, LineDatasetR10]:
    """
    Load frozen reference source cohort (`r9_6` or `r9_r10_8`) via `load_source_cohort_r12`,
    optionally excluding specified sources (e.g. `vcc2025_h1` for `PREREG_H §9 SA-H1`).
    """
    base = load_source_cohort_r12(
        source_mode=source_mode,
        representation=representation,
        common_genes=common_genes,
        hgnc_bridge=hgnc_bridge,
    )
    excl = {str(s).strip() for s in exclude_sources if str(s).strip()}
    if not excl:
        return base
    return {k: v for k, v in base.items() if k not in excl}


def _resolve_target_symbol(raw_sym: str, hgnc_lookup: Optional[HGNCLookup], hgnc_bridge: Optional[HGNCBridge]) -> str:
    s = str(raw_sym).strip()
    if not s:
        return ""
    if hgnc_lookup is not None:
        canon, _ = hgnc_lookup.resolve(s, "")
        return canon or s
    if hgnc_bridge is not None:
        # HGNCBridge compatibility if resolve/map method differs; fall back to HGNCLookup
        for attr in ("resolve_symbol", "canonicalize", "map_symbol", "resolve"):
            fn = getattr(hgnc_bridge, attr, None)
            if callable(fn):
                res = fn(s)
                if isinstance(res, tuple):
                    return str(res[0] or s)
                if isinstance(res, str) and res:
                    return res
    return s


def load_vcc2026_panel_overlap_targets(
    overlap_file: Optional[Path] = None,
    hgnc_lookup: Optional[HGNCLookup] = None,
    hgnc_bridge: Optional[HGNCBridge] = None,
) -> Tuple[Optional[Set[str]], Dict[str, Any]]:
    """
    Load the optional `129` shared val-test perturbation list for `PREREG_H §9.2 SA-Panel`.

    If `overlap_file` is None and the default file does not exist yet (pre-release state),
    returns `(None, {"available": False, "status": "PANEL_OVERLAP_LIST_NOT_AVAILABLE", ...})`.
    If an explicit `overlap_file` is passed and missing, raises `FileNotFoundError`.
    """
    explicit = overlap_file is not None
    path = Path(overlap_file) if explicit else DEFAULT_PANEL_OVERLAP_129_FILE
    if not path.exists():
        if explicit:
            raise FileNotFoundError(f"Explicit SA-Panel overlap target file not found: {path}")
        return None, {
            "available": False,
            "status": "PANEL_OVERLAP_LIST_NOT_AVAILABLE",
            "path": str(path),
            "n_targets": 0,
        }

    hl = hgnc_lookup if hgnc_lookup is not None else build_hgnc_lookup()
    raw_syms: List[str] = []
    suf = path.suffix.lower()
    if suf in {".csv", ".tsv"}:
        df = pd.read_csv(path, sep="\t" if suf == ".tsv" else ",")
        if df.empty:
            raw_syms = []
        else:
            col = None
            lower_map = {c.lower(): c for c in df.columns}
            for cand in ("target_gene", "canon_sym", "gene", "gene_symbol", "symbol", "perturbation", "target"):
                if cand in lower_map:
                    col = lower_map[cand]
                    break
            if col is None:
                col = str(df.columns[0])
            raw_syms = [str(v).strip() for v in df[col].dropna().tolist()]
    else:
        for line in path.read_text(encoding="utf-8").splitlines():
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            raw_syms.append(s.split(",")[0].split("\t")[0].strip())

    resolved: Set[str] = set()
    for sym in raw_syms:
        if not sym or sym.lower() in {"non-targeting", "ntc", "target_gene", "gene", "symbol"}:
            continue
        cs = _resolve_target_symbol(sym, hl, hgnc_bridge)
        if cs:
            resolved.add(cs)

    return resolved, {
        "available": True,
        "status": "OK",
        "path": str(path),
        "n_targets": len(resolved),
    }


def load_context_to_line_map(
    context_map: Optional[Dict[str, str]] = None,
    context_map_file: Optional[Path] = None,
) -> Tuple[Optional[Dict[str, str]], Dict[str, Any]]:
    """
    Load an official Arc `context -> cell_line` mapping table (`PREREG_H §3.1(b)`).

    Accepts JSON dict `{"Context A": "cal33", ...}`, JSON object `{"mapping": ..., "source_url": ...}`,
    or CSV/TSV with columns `(context, cell_line|line)`. Normalizes every target value and preserves
    `"jurkat"` so downstream ingestion can immediately drop all `jurkat` cells prior to QC.
    """
    raw_dict: Optional[Dict[str, str]] = None
    source_meta: Dict[str, Any] = {"source": None}
    if context_map is not None:
        raw_dict = {str(k).strip(): str(v).strip() for k, v in context_map.items()}
        source_meta = {"source": "in_memory_dict", "n_entries": len(raw_dict)}
    else:
        path = Path(context_map_file) if context_map_file is not None else DEFAULT_ARC_CONTEXT_MAP_FILE
        if not path.exists():
            if context_map_file is not None:
                raise FileNotFoundError(f"Explicit Arc context-to-line map file not found: {path}")
            return None, {"source": None, "status": "CONTEXT_MAP_NOT_PROVIDED", "path": str(path)}
        source_meta["path"] = str(path)
        source_meta["source"] = "file"
        if path.suffix.lower() == ".json":
            obj = json.loads(path.read_text(encoding="utf-8"))
            if isinstance(obj, dict) and isinstance(obj.get("mapping"), dict):
                raw_dict = {str(k).strip(): str(v).strip() for k, v in obj["mapping"].items()}
                source_meta["source_url"] = obj.get("source_url", "")
            elif isinstance(obj, dict):
                raw_dict = {str(k).strip(): str(v).strip() for k, v in obj.items() if not str(k).startswith("_")}
            else:
                raise ValueError(f"Invalid context map JSON format in {path}")
        else:
            sep = "\t" if path.suffix.lower() == ".tsv" else ","
            with path.open("r", encoding="utf-8") as fh:
                reader = csv.DictReader(fh, delimiter=sep)
                if not reader.fieldnames:
                    raise ValueError(f"Empty context map table: {path}")
                lmap = {c.lower(): c for c in reader.fieldnames}
                c_col = lmap.get("context") or lmap.get("context_id")
                l_col = lmap.get("cell_line") or lmap.get("line") or lmap.get("target_line")
                if not c_col or not l_col:
                    raise ValueError(f"Context map {path} must have 'context' and 'cell_line'/'line' columns; got {reader.fieldnames}")
                raw_dict = {str(r[c_col]).strip(): str(r[l_col]).strip() for r in reader if str(r[c_col]).strip()}

    out: Dict[str, str] = {}
    for ctx_key, line_val in raw_dict.items():
        canon = normalize_vcc_line_token(line_val)
        if canon is None:
            raise ValueError(
                f"EXCLUDED_FEASIBILITY: context map entry {ctx_key!r} -> {line_val!r} "
                "does not resolve to a valid PREREG_H line or 'jurkat'."
            )
        out[str(ctx_key).strip()] = canon
    source_meta["n_entries"] = len(out)
    return out, source_meta


def load_line_manifest(
    line_manifest: Optional[Dict[str, Any]] = None,
    line_manifest_file: Optional[Path] = None,
) -> Tuple[Optional[Dict[str, Dict[str, str]]], Dict[str, Any]]:
    """
    Load an official Arc per-file single-line manifest (`PREREG_H Addendum A §3.1(c)`).

    Accepts dict or JSON/CSV/TSV mapping `filename_or_relpath` to either:
      - `{"line": <cell_line>, "url": <release_url>, "quote": <exact_sentence_or_entry>}`
      - `<cell_line>` (optionally paired with top-level `"source_url"` / `"quote"`)
    File names alone without an explicit official manifest NEVER confer a line label.
    """
    raw_entries: Optional[Dict[str, Any]] = None
    default_url = ""
    source_meta: Dict[str, Any] = {"source": None}
    if line_manifest is not None:
        if isinstance(line_manifest, dict) and isinstance(line_manifest.get("files"), dict):
            default_url = str(line_manifest.get("source_url", "")).strip()
            raw_entries = dict(line_manifest["files"])
        else:
            raw_entries = {str(k).strip(): v for k, v in line_manifest.items() if not str(k).startswith("_")}
        source_meta = {"source": "in_memory_dict", "n_entries": len(raw_entries)}
    else:
        path = Path(line_manifest_file) if line_manifest_file is not None else DEFAULT_LINE_MANIFEST_FILE
        if not path.exists():
            if line_manifest_file is not None:
                raise FileNotFoundError(f"Explicit Arc line manifest file not found: {path}")
            return None, {"source": None, "status": "LINE_MANIFEST_NOT_PROVIDED", "path": str(path)}
        source_meta["path"] = str(path)
        source_meta["source"] = "file"
        if path.suffix.lower() == ".json":
            obj = json.loads(path.read_text(encoding="utf-8"))
            if isinstance(obj, dict) and isinstance(obj.get("files"), dict):
                default_url = str(obj.get("source_url", "")).strip()
                raw_entries = dict(obj["files"])
            elif isinstance(obj, dict) and isinstance(obj.get("entries"), dict):
                default_url = str(obj.get("source_url", "")).strip()
                raw_entries = dict(obj["entries"])
            elif isinstance(obj, dict):
                default_url = str(obj.get("source_url", "")).strip()
                raw_entries = {str(k).strip(): v for k, v in obj.items() if str(k) != "source_url" and not str(k).startswith("_")}
            else:
                raise ValueError(f"Invalid line manifest JSON format in {path}")
        else:
            sep = "\t" if path.suffix.lower() == ".tsv" else ","
            with path.open("r", encoding="utf-8") as fh:
                reader = csv.DictReader(fh, delimiter=sep)
                if not reader.fieldnames:
                    raise ValueError(f"Empty line manifest table: {path}")
                lmap = {c.lower(): c for c in reader.fieldnames}
                f_col = lmap.get("file") or lmap.get("filename") or lmap.get("path")
                l_col = lmap.get("cell_line") or lmap.get("line")
                u_col = lmap.get("url") or lmap.get("source_url")
                q_col = lmap.get("quote") or lmap.get("statement") or lmap.get("manifest_entry")
                if not f_col or not l_col:
                    raise ValueError(f"Line manifest {path} requires 'file'/'filename' and 'line'/'cell_line'; got {reader.fieldnames}")
                raw_entries = {}
                for row in reader:
                    fk = str(row[f_col]).strip()
                    if fk:
                        raw_entries[fk] = {
                            "line": str(row[l_col]).strip(),
                            "url": str(row.get(u_col, "")).strip() if u_col else "",
                            "quote": str(row.get(q_col, "")).strip() if q_col else "",
                        }

    out: Dict[str, Dict[str, str]] = {}
    for fk, spec in raw_entries.items():
        if isinstance(spec, dict):
            raw_line = str(spec.get("line") or spec.get("cell_line") or "").strip()
            url = str(spec.get("url") or spec.get("source_url") or default_url).strip()
            quote = str(spec.get("quote") or spec.get("statement") or spec.get("manifest_entry") or "").strip()
        else:
            raw_line = str(spec).strip()
            url = default_url
            quote = ""
        canon = normalize_vcc_line_token(raw_line)
        if canon is None:
            raise ValueError(
                f"EXCLUDED_FEASIBILITY: line_manifest entry {fk!r} -> {raw_line!r} "
                "does not resolve to a valid PREREG_H line or 'jurkat'."
            )
        norm_fk = str(fk).replace("\\", "/").strip()
        out[norm_fk] = {"line": canon, "url": url, "quote": quote}
    source_meta["n_entries"] = len(out)
    return out, source_meta


__all__ = [
    "ALLOWED_R15_LINES",
    "ALLOWED_R15_LINES_SET",
    "CODE_FREEZE_R15_FILE",
    "DEFAULT_ARC_CONTEXT_MAP_FILE",
    "DEFAULT_HGNC_FILE",
    "DEFAULT_LINE_MANIFEST_FILE",
    "DEFAULT_PANEL_OVERLAP_129_FILE",
    "DEFAULT_R15_SEED",
    "EXCLUDED_R15_LINES",
    "HGNCBridge",
    "HGNCLookup",
    "LineDatasetR10",
    "PREREG_H_COMMIT_FILE",
    "R9_COMMON_GENES_FILE",
    "R9_LINES_DIR",
    "R10_LINES_DIR",
    "R11_LINES_DIR",
    "R15_CONFIRMATORY_LINES_FILE",
    "R15_LINES_DIR",
    "R15_PREREG_DIR",
    "R15_RAW_DIR",
    "R15_RESULTS_DIR",
    "R15_SHARED_DIR",
    "R15_STATUS_DIR",
    "SA_H1_EXCLUDED_SOURCE",
    "build_hgnc_lookup",
    "check_prereg_h_gate",
    "list_standardized_lines",
    "load_context_to_line_map",
    "load_line_manifest",
    "load_r9_common_genes",
    "load_source_cohort_r15",
    "load_standardized_line_r10",
    "load_vcc2026_panel_overlap_targets",
    "normalize_vcc_line_token",
    "parse_confirmatory_lines_r15",
    "validate_r15_candidate_line",
]
