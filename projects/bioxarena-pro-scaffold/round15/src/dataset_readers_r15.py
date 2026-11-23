#!/usr/bin/env python3
"""
Round 15 (`PREREG_H` v2) VCC-Format `.h5ad` Feasibility Checker & Standardized Line Builder.

Supports BOTH published/synthetic `.obs` schemas:
  - 2025 Atlas style: `target_gene`, `guide_id` (allowing single-construct dual-promoter
    `GENE_1|GENE_2` / `non-targeting_1|non-targeting_2`), optional `batch` + line column
  - 2026 Control-Pack style: `target_gene`, `context`, `ntc_id`, optional `batch`/`channel`

Implements `PREREG_H` v2 locks strictly prior to metric evaluation:
  - `§3` / `§3.1`: Candidate lines in `{cal33, hela, ac16, molm13, u87mg}`; hard `jurkat`
    rejection; explicit line column (`cell_line`, `cell_type`, `line`, case/punct-insensitive)
    then official Arc `context -> line` mapping table ONLY; immediate `jurkat` drop prior
    to ANY RNA/guide/QC counting; `EXCLUDED_FEASIBILITY` if cells lack unambiguous lines.
  - `§5`: 100%-entry integer raw-UMI verification via `_resolve_h5ad_integer_slot`
    (`SLOT_AMBIGUOUS` on conflicting integer slots, `EXCLUDED_FEASIBILITY` if none);
    NTC strictly `target_gene == "non-targeting"`; single-construct filter (+ `PREREG_F §3.2`
    EM fallback fitted strictly on `target_line` cells if `.obs` assignments omitted);
    frozen HGNC canonicalization; RNA QC `n_umis >= 300`, `n_genes >= 200`, `max_pct_mt = None`.
  - `§6`: Two-tier disjoint split-half (Tier 1 lexicographic 1-based channel/batch parity
    with `>= 15` NTC/half & `_apply_min_ntc_batch_fallback`; Tier 2 stratified
    `np.random.default_rng(20261231)` sorted by `(barcode, raw_row_idx)`).
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROUND15_SRC_DIR = Path(__file__).resolve().parent
REPO_ROOT = ROUND15_SRC_DIR.parents[1]
R10_SCRIPTS_DIR = REPO_ROOT / "round10" / "shared" / "scripts"
R11_SRC_DIR = REPO_ROOT / "round11" / "src"
R12_SRC_DIR = REPO_ROOT / "round12" / "src"
R14_SRC_DIR = REPO_ROOT / "round14" / "src"
for _p in (str(ROUND15_SRC_DIR), str(R10_SCRIPTS_DIR), str(R11_SRC_DIR), str(R12_SRC_DIR), str(R14_SRC_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from data_loader_r15 import (  # noqa: E402
    DEFAULT_R15_SEED,
    R15_LINES_DIR,
    R15_RAW_DIR,
    load_context_to_line_map,
    load_line_manifest,
    normalize_vcc_line_token,
    validate_r15_candidate_line,
)
from dataset_readers_r11 import _apply_min_ntc_batch_fallback, _is_plain_single_gene_symbol  # noqa: E402
from dataset_readers_r12 import _resolve_h5ad_index_col  # noqa: E402
from dataset_readers_r14 import (  # noqa: E402
    _UNASSIGNED_TOKENS,
    _extract_kept_h5ad_csr,
    _extract_var_metadata_r14,
    _get_h5ad_slot_node,
    _get_var_group_for_slot,
    _inspect_slot_shape_and_layout,
    _make_csr_part_callbacks,
    _resolve_h5ad_integer_slot,
    _stream_h5ad_rna_qc_and_guides,
)
from prepare_r10_lines import (  # noqa: E402
    HGNCLookup,
    PartSpec,
    build_hgnc_lookup,
    extract_top2_csc,
    fit_poisson_gaussian_thresholds_from_hist,
    process_r10_line,
    read_obs_or_var_col,
)

_MISSING_STR_TOKENS = frozenset({"", "nan", "none", "na", "null", "<na>", "nat"})
_SAFE_EXACT_TOKENS = frozenset(
    {"safe-targeting", "safe_targeting", "safetargeting", "gene_desert", "gene-desert", "ctrl_t", "intergenic", "safe"}
)
_EXPLICIT_LINE_COL_CANDIDATES = ("cell_line", "cell_type", "line", "cellline", "celltype")
_CONTEXT_COL_CANDIDATES = ("context", "cell_context", "context_id")
_BATCH_COL_CANDIDATES = (
    "batch",
    "flex_channel",
    "channel",
    "lane",
    "replicate",
    "tech_rep",
    "technical_replicate",
    "replicate_id",
    "sample",
)
_TARGET_COL_CANDIDATES = ("target_gene", "gene_target", "perturbation", "gene", "target")
_GUIDE_COL_CANDIDATES = ("guide_id", "sgrna_id", "sgrna", "sg_id", "guide", "feature_call")
_NTC_COL_CANDIDATES = ("ntc_id", "ntc", "control_id")
_CLASS_COL_CANDIDATES = ("perturbation_type", "class", "guide_class", "control_type")
_NCOUNT_COL_CANDIDATES = ("num_features", "n_guides", "num_guides", "n_sgrnas", "nperts")
_MULTI_DELIM_RE = re.compile(r"[,;+]")


def _clean_str(val: Any) -> str:
    if val is None:
        return ""
    s = str(val).strip()
    return "" if s.lower() in _MISSING_STR_TOKENS else s


def _pick_obs_col(obs_keys: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
    low_map = {str(k).lower(): str(k) for k in obs_keys}
    for cand in candidates:
        if cand in obs_keys:
            return cand
        if cand.lower() in low_map:
            return low_map[cand.lower()]
    return None


def _find_vcc_h5ad_files(raw_dir: Path, line: str) -> List[Path]:
    p = Path(raw_dir)
    if p.is_file():
        if p.suffix.lower() != ".h5ad":
            raise RuntimeError(f"EXCLUDED_FEASIBILITY: Expected .h5ad input for line '{line}', got {p}")
        return [p]
    if not p.is_dir():
        raise RuntimeError(f"EXCLUDED_FEASIBILITY: Raw input directory does not exist for '{line}': {p}")
    files = sorted(
        [fp for fp in p.rglob("*.h5ad") if fp.is_file() and not fp.name.startswith(".")],
        key=lambda fp: fp.relative_to(p).as_posix(),
    )
    if not files:
        raise RuntimeError(f"EXCLUDED_FEASIBILITY: Zero .h5ad files found under {p} for line '{line}'.")
    return files


def _lookup_manifest_entry(
    file_name: str,
    rel_path: str,
    line_manifest: Optional[Dict[str, Dict[str, str]]],
) -> Optional[Dict[str, str]]:
    if not line_manifest:
        return None
    rel_norm = str(rel_path).replace("\\", "/").strip()
    if rel_norm in line_manifest:
        return line_manifest[rel_norm]
    if file_name in line_manifest:
        return line_manifest[file_name]
    return None


def _resolve_obs_cell_lines_v2(
    obs_grp: h5py.Group,
    n_obs: int,
    file_name: str,
    rel_path: str,
    context_map: Optional[Dict[str, str]],
    context_map_meta: Dict[str, Any],
    line_manifest: Optional[Dict[str, Dict[str, str]]] = None,
    line_manifest_meta: Optional[Dict[str, Any]] = None,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    Resolve per-cell line identity strictly following `PREREG_H §3.1 (v2)` + `Addendum A`:
      (a) Explicit cell-line column (`cell_line`, `cell_type`, `line`, ...);
      (b) Else `context` column mapped strictly through official Arc table `context_map`;
      (c) Else explicit official Arc `line_manifest` stating this file holds one named line
          (file names alone NEVER count; if `(a)` or mapped `(b)` disagrees with `(c)` for
          any cell, fail with `EXCLUDED_FEASIBILITY`).
    Never infers line identity from expression. Raises `EXCLUDED_FEASIBILITY` if ambiguous.
    """
    obs_keys = [str(k) for k in obs_grp.keys() if not str(k).startswith("_")]
    exp_col = _pick_obs_col(obs_keys, _EXPLICIT_LINE_COL_CANDIDATES)
    ctx_col = _pick_obs_col(obs_keys, _CONTEXT_COL_CANDIDATES)
    man_entry = _lookup_manifest_entry(file_name=file_name, rel_path=rel_path, line_manifest=line_manifest)

    if exp_col is not None:
        raw_vals = read_obs_or_var_col(obs_grp, exp_col)
        if len(raw_vals) != n_obs:
            raise RuntimeError(f"EXCLUDED_FEASIBILITY: Length mismatch in {file_name} obs[{exp_col}]")
        mapped: List[str] = []
        bad: List[str] = []
        for rv in raw_vals:
            tok = normalize_vcc_line_token(_clean_str(rv))
            if tok is None:
                bad.append(str(rv))
            else:
                mapped.append(tok)
        if bad:
            raise RuntimeError(
                f"EXCLUDED_FEASIBILITY: Explicit line column {exp_col!r} in {file_name} has unresolved/ambiguous "
                f"values (examples={sorted(set(bad))[:5]})."
            )
        if man_entry is not None and any(x != man_entry["line"] for x in mapped):
            obs_lines = sorted(set(mapped))
            raise RuntimeError(
                f"EXCLUDED_FEASIBILITY: PREREG_H Addendum A conflict in {file_name}: explicit column "
                f"{exp_col!r} resolves lines {obs_lines}, disagreeing with line_manifest {man_entry['line']!r}."
            )
        return np.asarray(mapped, dtype=object), {
            "rule": f"explicit_col:{exp_col}",
            "column": exp_col,
            "context_map_used": False,
            "manifest_consistent": (man_entry is not None),
        }

    if ctx_col is not None and context_map:
        raw_ctx = read_obs_or_var_col(obs_grp, ctx_col)
        if len(raw_ctx) != n_obs:
            raise RuntimeError(f"EXCLUDED_FEASIBILITY: Length mismatch in {file_name} obs[{ctx_col}]")
        ci_map = {str(k).strip().lower(): v for k, v in context_map.items()}
        mapped = []
        unmapped: Set[str] = set()
        for rv in raw_ctx:
            key = _clean_str(rv)
            val = context_map.get(key)
            if val is None and key:
                val = ci_map.get(key.lower())
            if val is None:
                unmapped.add(str(rv))
            else:
                mapped.append(val)
        if unmapped:
            raise RuntimeError(
                f"EXCLUDED_FEASIBILITY: Unmapped .obs[{ctx_col!r}] values in {file_name} not present in "
                f"official context map: {sorted(unmapped)[:8]}."
            )
        if man_entry is not None and any(x != man_entry["line"] for x in mapped):
            ctx_lines = sorted(set(mapped))
            raise RuntimeError(
                f"EXCLUDED_FEASIBILITY: PREREG_H Addendum A conflict in {file_name}: mapped .obs[{ctx_col!r}] "
                f"resolves lines {ctx_lines}, disagreeing with line_manifest {man_entry['line']!r}."
            )
        return np.asarray(mapped, dtype=object), {
            "rule": f"context_map:{ctx_col}",
            "column": ctx_col,
            "context_map_used": True,
            "context_map_source": context_map_meta,
            "manifest_consistent": (man_entry is not None),
        }

    if man_entry is not None:
        return np.full(n_obs, man_entry["line"], dtype=object), {
            "rule": f"line_manifest:{file_name}",
            "manifest_line": man_entry["line"],
            "manifest_url": man_entry.get("url", ""),
            "manifest_quote": man_entry.get("quote", ""),
            "manifest_source": line_manifest_meta or {},
            "context_map_used": False,
        }

    if ctx_col is not None:
        raise RuntimeError(
            f"EXCLUDED_FEASIBILITY: {file_name} uses .obs[{ctx_col!r}] without an official Arc "
            "context-to-line mapping table or line_manifest entry (PREREG_H §3.1 + Addendum A)."
        )

    raise RuntimeError(
        f"EXCLUDED_FEASIBILITY: {file_name} lacks an explicit cell-line column ({_EXPLICIT_LINE_COL_CANDIDATES}), "
        "a mapped 'context' column, and an official line_manifest entry (PREREG_H §3.1 + Addendum A; filename alone invalid)."
    )


def _detect_vcc_obs_columns(obs_grp: h5py.Group) -> Dict[str, Optional[str]]:
    obs_keys = [str(k) for k in obs_grp.keys() if not str(k).startswith("_")]
    return {
        "target_col": _pick_obs_col(obs_keys, _TARGET_COL_CANDIDATES),
        "guide_col": _pick_obs_col(obs_keys, _GUIDE_COL_CANDIDATES),
        "ntc_col": _pick_obs_col(obs_keys, _NTC_COL_CANDIDATES),
        "class_col": _pick_obs_col(obs_keys, _CLASS_COL_CANDIDATES),
        "ncount_col": _pick_obs_col(obs_keys, _NCOUNT_COL_CANDIDATES),
        "batch_col": _pick_obs_col(obs_keys, _BATCH_COL_CANDIDATES),
    }


def _has_token(text: str, tokens: Sequence[str]) -> bool:
    low = text.lower().strip()
    return low in tokens or any(tok in low for tok in ("safe-targeting", "safe_targeting", "gene_desert", "ctrl_t"))


def _is_ntc_cassette_segment(seg: str) -> bool:
    low = seg.strip().lower()
    if low == "non-targeting":
        return True
    return low.startswith(("non-targeting", "nontargeting", "ntc_", "ntc-", "ntc"))


def _extract_guide_gene_prefix(seg: str) -> str:
    s = seg.strip()
    for sep in ("_", ".", ":"):
        if sep in s:
            return s.split(sep, 1)[0].strip()
    return s.split("-", 1)[0].strip() if "-" in s else s


def _classify_vcc_obs_cell(
    tgt_raw: str,
    guide_raw: str,
    ntc_raw: str,
    cls_raw: str,
    n_feat: Optional[float],
    ncount_is_nperts: bool,
    hgnc: HGNCLookup,
) -> Tuple[str, str, Optional[Tuple[str, str, str, str]]]:
    """
    Classify one QC-eligible cell's `.obs` perturbation metadata under `PREREG_H §5.2–5.4`.

    Returns `(category, pert_key, meta_tuple)` where `category` is one of:
      - `"ntc"`: strict `non-targeting` single construct
      - `"target"`: single HGNC-canonicalized target gene
      - `"safe"`: safe-targeting / gene_desert control (excluded)
      - `"multi_or_zero"`: multiplet / doublet / unassigned / dual-target conflict (excluded)
      - `"unmapped"`: non-HGNC target symbol (excluded)
    """
    t = _clean_str(tgt_raw)
    g = _clean_str(guide_raw)
    ntc_s = _clean_str(ntc_raw)
    cls_s = _clean_str(cls_raw)

    if n_feat is not None and np.isfinite(n_feat):
        nf = int(n_feat)
        if nf > 1 or nf < 0:
            return "multi_or_zero", "", None
        if nf == 0 and not (ncount_is_nperts and t.lower() == "non-targeting"):
            return "multi_or_zero", "", None

    if _MULTI_DELIM_RE.search(t) or _MULTI_DELIM_RE.search(g) or _MULTI_DELIM_RE.search(ntc_s):
        return "multi_or_zero", "", None
    if "|" in t or len(t.split()) > 1:
        return "multi_or_zero", "", None

    t_low, g_low, ntc_low, cls_low = t.lower(), g.lower(), ntc_s.lower(), cls_s.lower()
    if any(x in _UNASSIGNED_TOKENS for x in (t_low, cls_low) if x):
        return "multi_or_zero", "", None
    if t_low in _SAFE_EXACT_TOKENS or _has_token(g_low, _SAFE_EXACT_TOKENS) or _has_token(cls_low, _SAFE_EXACT_TOKENS) or _has_token(ntc_low, _SAFE_EXACT_TOKENS):
        return "safe", "", None
    if any(tok in cls_low for tok in ("multiplet", "doublet", "unassigned", "ambiguous")):
        return "multi_or_zero", "", None

    if not t:
        return "multi_or_zero", "", None

    # Locked NTC rule (`PREREG_H §5.3`): strictly `target_gene == "non-targeting"`.
    if t_low == "non-targeting":
        if "|" in g:
            segs = [s.strip() for s in g.split("|")]
            if not segs or any((not s or not _is_ntc_cassette_segment(s)) for s in segs):
                return "multi_or_zero", "", None
        elif g_low in _UNASSIGNED_TOKENS and ntc_low == "" and g != "":
            return "multi_or_zero", "", None
        if any(tok in g_low for tok in ("multiplet", "doublet", "unassigned")):
            return "multi_or_zero", "", None
        label = g or ntc_s or "non-targeting"
        return "ntc", "non-targeting", ("", "non-targeting", label, "non-targeting")

    # Single targeting construct: verify optional 2026 `ntc_id` does not claim NTC, and verify
    # any 2025 dual-promoter `|` cassette in `guide_id` belongs exclusively to this single target.
    if ntc_s and (ntc_low != "-" and (_is_ntc_cassette_segment(ntc_s) or "non-targeting" in ntc_low or ntc_low.startswith("ntc"))):
        return "multi_or_zero", "", None
    if not _is_plain_single_gene_symbol(t):
        return "unmapped", "", None

    canon, hgnc_status = hgnc.resolve(t, "")
    if not canon or canon not in hgnc.approved:
        return "unmapped", "", None

    if "|" in g:
        segs = [s.strip() for s in g.split("|")]
        if not segs or any(not s for s in segs):
            return "multi_or_zero", "", None
        t_boundary = re.compile(rf"(?:^|[_.:\s/]){re.escape(t_low)}(?:$|[_.:\s/])")
        for seg in segs:
            if _is_ntc_cassette_segment(seg) or _has_token(seg, _SAFE_EXACT_TOKENS):
                return "multi_or_zero", "", None
            pref = _extract_guide_gene_prefix(seg)
            pref_canon, _ = hgnc.resolve(pref, "") if _is_plain_single_gene_symbol(pref) else ("", "")
            if pref_canon in hgnc.approved and pref_canon != canon:
                return "multi_or_zero", "", None
            if pref.lower() != t_low and pref_canon != canon:
                if not t_boundary.search(seg.lower()):
                    return "multi_or_zero", "", None

    ensg = str(hgnc.sym_to_ens.get(canon, ""))
    return "target", canon, (ensg, canon, t, str(hgnc_status))


def _load_guide_library_fallback(guide_lib_file: Optional[Path], hgnc: HGNCLookup) -> Dict[str, Any]:
    ntc_keys: Set[str] = set()
    safe_keys: Set[str] = set()
    guide_to_target: Dict[str, Tuple[str, str, str, str]] = {}
    if guide_lib_file is None:
        return {"ntc_keys": ntc_keys, "safe_keys": safe_keys, "guide_to_target": guide_to_target}
    p = Path(guide_lib_file)
    if not p.exists():
        raise FileNotFoundError(f"Guide library file not found: {p}")
    sep = "\t" if p.suffix.lower().endswith("tsv") or ".tsv" in p.name.lower() else ","
    with p.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter=sep)
        for row in reader:
            gid = _clean_str(row.get("id") or row.get("guide_id") or row.get("sgrna_id") or row.get("name"))
            sym = _clean_str(row.get("symbol") or row.get("target_gene") or row.get("gene"))
            cls_s = _clean_str(row.get("class") or row.get("perturbation_type"))
            if not gid:
                continue
            if _has_token(sym, _SAFE_EXACT_TOKENS) or _has_token(cls_s, _SAFE_EXACT_TOKENS):
                safe_keys.add(gid)
            elif sym.lower() == "non-targeting" or _is_ntc_cassette_segment(sym):
                ntc_keys.add(gid)
            else:
                cs, st = hgnc.resolve(sym, "")
                if cs and cs in hgnc.approved:
                    guide_to_target[gid] = (str(hgnc.sym_to_ens.get(cs, "")), cs, sym, str(st))
    return {"ntc_keys": ntc_keys, "safe_keys": safe_keys, "guide_to_target": guide_to_target}


def _classify_em_guide_key(key: str, lib_meta: Dict[str, Any], hgnc: HGNCLookup) -> Tuple[str, Optional[Tuple[str, str, str, str]]]:
    k = _clean_str(key)
    if not k:
        return "unmapped", None
    if k in lib_meta["safe_keys"] or _has_token(k, _SAFE_EXACT_TOKENS):
        return "safe", None
    if k in lib_meta["ntc_keys"] or _is_ntc_cassette_segment(k):
        return "ntc", ("", "non-targeting", k, "non-targeting")
    if k in lib_meta["guide_to_target"]:
        return "target", lib_meta["guide_to_target"][k]
    pref = _extract_guide_gene_prefix(k)
    if _is_plain_single_gene_symbol(pref):
        cs, st = hgnc.resolve(pref, "")
        if cs and cs in hgnc.approved:
            return "target", (str(hgnc.sym_to_ens.get(cs, "")), cs, pref, str(st))
    return "unmapped", None


def _compute_line_crp_hist(
    x_node: Any,
    layout: str,
    line_row_idx: np.ndarray,
    ng_line: np.ndarray,
    crp_rows: np.ndarray,
    n_obs: int,
    n_var: int,
    min_ngenes: int,
) -> np.ndarray:
    """Compute `(n_guides, 512)` guide UMI histogram strictly over `target_line` cells with `ngenes >= min_ngenes`."""
    n_crp = len(crp_rows)
    hist = np.zeros((n_crp, 512), dtype=np.int64)
    if n_crp == 0:
        return hist
    fit_rows = line_row_idx[ng_line >= int(min_ngenes)]
    if len(fit_rows) == 0:
        return hist
    crp_mat = _extract_kept_h5ad_csr(
        x_node=x_node,
        layout=layout,
        keep_arr=fit_rows,
        col_rows=crp_rows,
        n_obs=n_obs,
        n_var=n_var,
    )
    if crp_mat.nnz == 0:
        return hist
    crp_csc = crp_mat.tocsc()
    for gi in range(n_crp):
        s, e = int(crp_csc.indptr[gi]), int(crp_csc.indptr[gi + 1])
        if e <= s:
            continue
        vals = np.rint(crp_csc.data[s:e]).astype(np.int64, copy=False)
        vals = np.clip(vals[vals > 0], 1, 511)
        if len(vals) > 0:
            np.add.at(hist[gi], vals, 1)
    return hist


def _assign_two_tier_splits(
    kept_meta: List[Dict[str, Any]],
    has_batch_col: bool,
    min_half_nt: int,
    min_batch_nt: int,
    seed: int,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    """
    Assign Half A (`1`) / Half B (`2`) and subtraction batch labels under `PREREG_H §6`.

    Tier 1: explicit `.obs` batch/channel column with `>= 2` valid distinct unit labels
    producing `>= min_half_nt` NTCs in BOTH odd (Half A) and even (Half B) lexicographical
    parity groups, with per-unit `< min_batch_nt` NTC falling back within half via
    `_apply_min_ntc_batch_fallback`.

    Tier 2: deterministic stratified 50/50 split with `np.random.default_rng(seed)` over
    sorted strata (`"__NTC__"` + sorted canonical symbols), sorting each stratum by
    `(barcode, raw_row_idx, file_idx)`.
    """
    n_kept = len(kept_meta)
    nt_flags = np.array([bool(m["is_nt"]) for m in kept_meta], dtype=bool)
    raw_batches = np.array([str(m["raw_batch"]) for m in kept_meta], dtype=object)

    tier1_reason = "missing_batch_column"
    if has_batch_col and n_kept > 0:
        valid_nonempty = all(bool(b) for b in raw_batches)
        uniq_units = sorted(set(raw_batches)) if valid_nonempty else []
        if not valid_nonempty:
            tier1_reason = "empty_or_nan_batch_labels"
        elif len(uniq_units) < 2:
            tier1_reason = f"fewer_than_2_batch_units ({len(uniq_units)})"
        else:
            unit_half = {u: (1 if ((idx + 1) % 2 == 1) else 2) for idx, u in enumerate(uniq_units)}
            half_flags = np.array([unit_half[b] for b in raw_batches], dtype=np.int32)
            nt_a = int(np.sum(nt_flags & (half_flags == 1)))
            nt_b = int(np.sum(nt_flags & (half_flags == 2)))
            if nt_a >= int(min_half_nt) and nt_b >= int(min_half_nt):
                batches_out, n_fb = _apply_min_ntc_batch_fallback(
                    raw_batches, nt_flags, half_flags, min_batch_nt=int(min_batch_nt)
                )
                return half_flags, np.asarray(batches_out, dtype=object), {
                    "split_tier": "tier1_batch_channel_parity",
                    "half_rule_used": "lex_1based_odd_HalfA_vs_even_HalfB",
                    "n_units": len(uniq_units),
                    "odd_half_a_units": [u for u in uniq_units if unit_half[u] == 1],
                    "even_half_b_units": [u for u in uniq_units if unit_half[u] == 2],
                    "n_ntc_half_a": nt_a,
                    "n_ntc_half_b": nt_b,
                    "n_batch_ntc_fallback_cells": int(n_fb),
                    "batches_used": sorted({str(x) for x in batches_out}),
                }
            tier1_reason = f"insufficient_ntc_under_parity (half_a={nt_a}, half_b={nt_b}, min={min_half_nt})"

    # Tier 2: stratified barcode-sorted RNG split (`seed = 20261231`)
    strata: Dict[str, List[int]] = {}
    for idx, m in enumerate(kept_meta):
        skey = "__NTC__" if m["is_nt"] else str(m["canon_sym"])
        strata.setdefault(skey, []).append(idx)

    half_flags = np.zeros(n_kept, dtype=np.int32)
    rng = np.random.default_rng(int(seed))
    for skey in sorted(strata.keys()):
        idxs = strata[skey]
        idxs.sort(key=lambda i: (str(kept_meta[i]["barcode"]), int(kept_meta[i]["raw_row_idx"]), int(kept_meta[i]["file_idx"])))
        arr = np.asarray(idxs, dtype=np.int64)
        perm = arr[rng.permutation(len(arr))]
        n_half_a = len(perm) // 2
        half_flags[perm[:n_half_a]] = 1
        half_flags[perm[n_half_a:]] = 2

    batches_out = np.array([f"half_{int(h)}" for h in half_flags], dtype=object)
    nt_a = int(np.sum(nt_flags & (half_flags == 1)))
    nt_b = int(np.sum(nt_flags & (half_flags == 2)))
    return half_flags, batches_out, {
        "split_tier": "tier2_stratified_seed",
        "tier1_fallback_reason": tier1_reason,
        "half_rule_used": f"stratified_barcode_sorted_seed_{int(seed)}",
        "split_seed": int(seed),
        "n_strata": len(strata),
        "n_ntc_half_a": nt_a,
        "n_ntc_half_b": nt_b,
        "n_batch_ntc_fallback_cells": 0,
        "batches_used": ["half_1", "half_2"],
    }


def _scan_vcc_inputs(
    raw_dir: Path,
    line: str,
    hgnc: HGNCLookup,
    x_layer: str = "auto",
    context_map: Optional[Dict[str, str]] = None,
    context_map_file: Optional[Path] = None,
    line_manifest: Optional[Dict[str, Any]] = None,
    line_manifest_file: Optional[Path] = None,
    guide_lib_file: Optional[Path] = None,
    min_umi: int = 300,
    min_ngenes: int = 200,
    max_pct_mt: Optional[float] = None,
    min_guide_thr: int = 5,
    default_guide_thr: int = 7,
    load_matrices: bool = False,
) -> Dict[str, Any]:
    canon_line = validate_r15_candidate_line(line)
    h5_files = _find_vcc_h5ad_files(raw_dir, canon_line)
    resolved_ctx_map, ctx_meta = load_context_to_line_map(context_map=context_map, context_map_file=context_map_file)
    resolved_manifest, man_meta = load_line_manifest(line_manifest=line_manifest, line_manifest_file=line_manifest_file)
    lib_meta = _load_guide_library_fallback(guide_lib_file, hgnc)

    file_records: List[Dict[str, Any]] = []
    n_jurkat_dropped_total = 0
    n_other_lines_dropped_total = 0
    n_line_raw_cells = 0
    line_identity_rules: Set[str] = set()
    used_obs_calls = False
    used_em_fallback = False
    base_dir = Path(raw_dir) if Path(raw_dir).is_dir() else Path(raw_dir).parent

    for file_idx, hp in enumerate(h5_files):
        rel_path = hp.relative_to(base_dir).as_posix() if hp.is_relative_to(base_dir) else hp.name
        with h5py.File(hp, "r") as f:
            obs_grp = f["obs"]
            obs_idx_col = _resolve_h5ad_index_col(obs_grp, group_name="obs")
            all_barcodes = [_clean_str(x) for x in read_obs_or_var_col(obs_grp, obs_idx_col)]
            n_obs = len(all_barcodes)

            cell_lines, id_meta = _resolve_obs_cell_lines_v2(
                obs_grp=obs_grp,
                n_obs=n_obs,
                file_name=hp.name,
                rel_path=rel_path,
                context_map=resolved_ctx_map,
                context_map_meta=ctx_meta,
                line_manifest=resolved_manifest,
                line_manifest_meta=man_meta,
            )
            line_identity_rules.add(str(id_meta["rule"]))

            # Immediately drop ALL Jurkat cells BEFORE inspecting slots or counting ANY QC metric (§3.1)
            is_jurkat = cell_lines == "jurkat"
            n_jurkat_file = int(np.sum(is_jurkat))
            n_jurkat_dropped_total += n_jurkat_file

            line_mask = cell_lines == canon_line
            line_rows = np.where(line_mask)[0]
            n_other_lines_dropped_total += int(np.sum((~line_mask) & (~is_jurkat)))
            if len(line_rows) == 0:
                continue
            n_line_raw_cells += len(line_rows)

            resolved_slot, int_slots, slot_diag = _resolve_h5ad_integer_slot(f, x_layer=x_layer, file_label=hp.name)
            var_grp = _get_var_group_for_slot(f, resolved_slot)
            vmeta = _extract_var_metadata_r14(var_grp)
            x_node = _get_h5ad_slot_node(f, resolved_slot)
            (mat_r, mat_c), layout = _inspect_slot_shape_and_layout(x_node, n_obs, vmeta["n_var"])
            if mat_r != n_obs or mat_c != vmeta["n_var"]:
                raise RuntimeError(
                    f"EXCLUDED_FEASIBILITY: Shape mismatch in {hp.name} ({resolved_slot}): "
                    f"{(mat_r, mat_c)} vs obs/var {(n_obs, vmeta['n_var'])}"
                )
            if len(vmeta["rna_rows"]) == 0:
                raise RuntimeError(f"EXCLUDED_FEASIBILITY: Zero human RNA features in {hp.name}.")

            obs_cols = _detect_vcc_obs_columns(obs_grp)
            has_obs = obs_cols["target_col"] is not None or obs_cols["guide_col"] is not None
            if has_obs:
                used_obs_calls = True
            elif len(vmeta["crp_rows"]) > 0:
                used_em_fallback = True
            else:
                raise RuntimeError(
                    f"EXCLUDED_FEASIBILITY: {hp.name} lacks both .obs target/guide columns and CRISPR guide capture rows."
                )

            slot_summary = {
                "path": str(hp),
                "x_slot_used": resolved_slot,
                "int_count_sources": int_slots,
                "slot_diagnostics": slot_diag,
                "selected_median_total_umi": float(slot_diag["slots"][resolved_slot]["median_total_umi"]),
                "layout": layout,
                "n_raw_file_cells": int(n_obs),
                "n_jurkat_dropped": n_jurkat_file,
                "n_target_line_raw_cells": int(len(line_rows)),
                "n_rna_features": int(len(vmeta["rna_rows"])),
                "line_identity": id_meta,
            }

            if not load_matrices:
                # Counts-only feasibility inspection strictly on non-Jurkat target-line rows
                tvals = read_obs_or_var_col(obs_grp, obs_cols["target_col"]) if obs_cols["target_col"] else [""] * n_obs
                gvals = read_obs_or_var_col(obs_grp, obs_cols["guide_col"]) if obs_cols["guide_col"] else [""] * n_obs
                nvals = read_obs_or_var_col(obs_grp, obs_cols["ntc_col"]) if obs_cols["ntc_col"] else [""] * n_obs
                cvals = read_obs_or_var_col(obs_grp, obs_cols["class_col"]) if obs_cols["class_col"] else [""] * n_obs
                kvals = (
                    pd.to_numeric(pd.Series(read_obs_or_var_col(obs_grp, obs_cols["ncount_col"])), errors="coerce").to_numpy()
                    if obs_cols["ncount_col"]
                    else None
                )
                nc_is_nperts = str(obs_cols["ncount_col"] or "").lower() == "nperts"
                file_records.append(
                    {
                        "slot_summary": slot_summary,
                        "obs_cols": obs_cols,
                        "has_obs": has_obs,
                        "line_rows": line_rows,
                        "tvals": [tvals[i] for i in line_rows] if has_obs else [],
                        "gvals": [gvals[i] for i in line_rows] if has_obs else [],
                        "nvals": [nvals[i] for i in line_rows] if has_obs else [],
                        "cvals": [cvals[i] for i in line_rows] if has_obs else [],
                        "kvals": kvals[line_rows] if (has_obs and kvals is not None) else None,
                        "nc_is_nperts": nc_is_nperts,
                        "vmeta": vmeta,
                    }
                )
                continue

            # Full ingestion path (`load_matrices=True`): RNA QC only over `target_line` rows (`Jurkat` excluded)
            umi_all, ng_all, pct_mt_all, _ = _stream_h5ad_rna_qc_and_guides(
                x_node=x_node,
                layout=layout,
                n_obs=n_obs,
                n_var=mat_c,
                rna_rows=vmeta["rna_rows"],
                mt_mask_rna=vmeta["mt_mask_rna"],
                crp_rows=np.empty(0, dtype=np.int64),
                min_ngenes_fit=min_ngenes,
            )
            umi_line = umi_all[line_rows]
            ng_line = ng_all[line_rows]
            pct_mt_line = pct_mt_all[line_rows]

            fit_hist = None
            if not has_obs:
                fit_hist = _compute_line_crp_hist(
                    x_node=x_node,
                    layout=layout,
                    line_row_idx=line_rows,
                    ng_line=ng_line,
                    crp_rows=vmeta["crp_rows"],
                    n_obs=n_obs,
                    n_var=mat_c,
                    min_ngenes=min_ngenes,
                )

            ok_qc = (umi_line >= float(min_umi)) & (ng_line >= int(min_ngenes))
            if max_pct_mt is not None:
                ok_qc &= pct_mt_line <= float(max_pct_mt)
            qc_pos = np.where(ok_qc)[0]
            qc_abs_rows = line_rows[qc_pos]
            n_rna_qc_excl_file = int(np.sum(~ok_qc))

            crp_qc_csc: Optional[sp.csc_matrix] = None
            crp_keys: List[str] = []
            if not has_obs and len(qc_abs_rows) > 0:
                crp_qc_csr = _extract_kept_h5ad_csr(
                    x_node=x_node,
                    layout=layout,
                    keep_arr=qc_abs_rows,
                    col_rows=vmeta["crp_rows"],
                    n_obs=n_obs,
                    n_var=mat_c,
                )
                crp_qc_csc = crp_qc_csr.T.tocsc()
                crp_keys = [
                    cid if (cid in lib_meta["ntc_keys"] or cid in lib_meta["safe_keys"] or cid in lib_meta["guide_to_target"]) else cname
                    for cid, cname in zip(vmeta["c_ids"], vmeta["c_names"])
                ]

            tvals = read_obs_or_var_col(obs_grp, obs_cols["target_col"]) if obs_cols["target_col"] else [""] * n_obs
            gvals = read_obs_or_var_col(obs_grp, obs_cols["guide_col"]) if obs_cols["guide_col"] else [""] * n_obs
            nvals = read_obs_or_var_col(obs_grp, obs_cols["ntc_col"]) if obs_cols["ntc_col"] else [""] * n_obs
            cvals = read_obs_or_var_col(obs_grp, obs_cols["class_col"]) if obs_cols["class_col"] else [""] * n_obs
            bvals = read_obs_or_var_col(obs_grp, obs_cols["batch_col"]) if obs_cols["batch_col"] else [""] * n_obs
            kvals = (
                pd.to_numeric(pd.Series(read_obs_or_var_col(obs_grp, obs_cols["ncount_col"])), errors="coerce").to_numpy()
                if obs_cols["ncount_col"]
                else None
            )

            file_records.append(
                {
                    "path": hp,
                    "file_idx": file_idx,
                    "slot_summary": slot_summary,
                    "slot": resolved_slot,
                    "layout": layout,
                    "n_obs": n_obs,
                    "n_var": mat_c,
                    "vmeta": vmeta,
                    "obs_cols": obs_cols,
                    "has_obs": has_obs,
                    "all_barcodes": all_barcodes,
                    "umi_all": umi_all,
                    "qc_abs_rows": qc_abs_rows,
                    "n_rna_qc_excl": n_rna_qc_excl_file,
                    "tvals": tvals,
                    "gvals": gvals,
                    "nvals": nvals,
                    "cvals": cvals,
                    "bvals": bvals,
                    "kvals": kvals,
                    "nc_is_nperts": str(obs_cols["ncount_col"] or "").lower() == "nperts",
                    "fit_hist": fit_hist,
                    "crp_qc_csc": crp_qc_csc,
                    "crp_keys": crp_keys,
                }
            )

    if len(file_records) == 0 or n_line_raw_cells == 0:
        raise RuntimeError(
            f"EXCLUDED_FEASIBILITY: Zero cells belonging to candidate line '{canon_line}' in {raw_dir} "
            f"(jurkat_dropped={n_jurkat_dropped_total})."
        )
    if used_obs_calls and used_em_fallback:
        raise RuntimeError(
            f"EXCLUDED_FEASIBILITY: Mixed .obs guide mode and CRISPR EM mode across files for '{canon_line}'."
        )

    return {
        "line": canon_line,
        "files": h5_files,
        "file_records": file_records,
        "line_identity_rules": sorted(line_identity_rules),
        "n_jurkat_dropped": int(n_jurkat_dropped_total),
        "n_other_lines_dropped": int(n_other_lines_dropped_total),
        "n_line_raw_cells": int(n_line_raw_cells),
        "used_obs_calls": bool(used_obs_calls),
        "used_em_fallback": bool(used_em_fallback),
        "lib_meta": lib_meta,
        "context_map_meta": ctx_meta,
    }


def check_vcc_feasibility(
    raw_dir: Path,
    line: str,
    hgnc: Optional[HGNCLookup] = None,
    x_layer: str = "auto",
    context_map: Optional[Dict[str, str]] = None,
    context_map_file: Optional[Path] = None,
    line_manifest: Optional[Dict[str, Any]] = None,
    line_manifest_file: Optional[Path] = None,
    guide_lib_file: Optional[Path] = None,
) -> Dict[str, Any]:
    """
    `PREREG_H §3.1 + Addendum A / §5` counts-only pre-metric feasibility gate (zero expression deltas or
    correlations computed). Rejects `jurkat`, unmapped `context`-only files, filename-only guesses,
    manifest-vs-obs conflicts, conflicting integer slots (`SLOT_AMBIGUOUS`), and missing integers (`EXCLUDED_FEASIBILITY`).
    """
    try:
        hgnc_l = hgnc or build_hgnc_lookup()
        scan = _scan_vcc_inputs(
            raw_dir=raw_dir,
            line=line,
            hgnc=hgnc_l,
            x_layer=x_layer,
            context_map=context_map,
            context_map_file=context_map_file,
            line_manifest=line_manifest,
            line_manifest_file=line_manifest_file,
            guide_lib_file=guide_lib_file,
            load_matrices=False,
        )
        n_obs_ntc = n_obs_safe = n_obs_tgt = 0
        canon_tgts: Set[str] = set()
        slots_meta: Dict[str, Any] = {}
        n_rna = 0
        has_batch_col = any(bool(rec["obs_cols"]["batch_col"]) for rec in scan["file_records"])
        for idx, rec in enumerate(scan["file_records"]):
            sm = rec["slot_summary"]
            slots_meta[f"file_{idx}"] = sm
            n_rna = int(sm["n_rna_features"])
            if rec["has_obs"]:
                kvals = rec["kvals"]
                for i in range(len(rec["tvals"])):
                    kf = float(kvals[i]) if kvals is not None else None
                    cat, pkey, _ = _classify_vcc_obs_cell(
                        tgt_raw=rec["tvals"][i],
                        guide_raw=rec["gvals"][i],
                        ntc_raw=rec["nvals"][i],
                        cls_raw=rec["cvals"][i],
                        n_feat=kf,
                        ncount_is_nperts=rec["nc_is_nperts"],
                        hgnc=hgnc_l,
                    )
                    if cat == "ntc":
                        n_obs_ntc += 1
                    elif cat == "safe":
                        n_obs_safe += 1
                    elif cat == "target":
                        n_obs_tgt += 1
                        canon_tgts.add(pkey)
            else:
                for cid, cname in zip(rec["vmeta"]["c_ids"], rec["vmeta"]["c_names"]):
                    ck = cid if (cid in scan["lib_meta"]["ntc_keys"] or cid in scan["lib_meta"]["safe_keys"] or cid in scan["lib_meta"]["guide_to_target"]) else cname
                    ck_cat, ck_meta = _classify_em_guide_key(ck, scan["lib_meta"], hgnc_l)
                    if ck_cat == "target" and ck_meta is not None:
                        canon_tgts.add(ck_meta[1])

        ok_guides = (n_obs_ntc > 0 and len(canon_tgts) > 0) if scan["used_obs_calls"] else len(canon_tgts) > 0
        return {
            "feasible": bool(ok_guides and n_rna > 0),
            "status": "OK" if (ok_guides and n_rna > 0) else "EXCLUDED_FEASIBILITY",
            "reason": "OK" if (ok_guides and n_rna > 0) else "Insufficient resolvable NTC/target guides after Jurkat exclusion",
            "line": scan["line"],
            "n_files_contributing": len(scan["file_records"]),
            "line_identity_rules": scan["line_identity_rules"],
            "n_jurkat_dropped": scan["n_jurkat_dropped"],
            "n_other_lines_dropped": scan["n_other_lines_dropped"],
            "n_target_line_raw_cells": scan["n_line_raw_cells"],
            "n_rna_features": n_rna,
            "has_batch_column": has_batch_col,
            "guide_call_mode": "obs_single_construct" if scan["used_obs_calls"] else "h5ad_crispr_em_fallback",
            "n_obs_ntc_meta": n_obs_ntc,
            "n_obs_safe_meta": n_obs_safe,
            "n_obs_tgt_meta": n_obs_tgt,
            "n_unique_canon_targets_meta": len(canon_tgts),
            "file_slots": slots_meta,
        }
    except Exception as exc:
        msg = str(exc)
        if isinstance(exc, ValueError) and msg.startswith("EXCLUDED_LINE:"):
            return {"feasible": False, "status": "EXCLUDED_LINE", "reason": msg}
        st = "SLOT_AMBIGUOUS" if "SLOT_AMBIGUOUS" in msg else "EXCLUDED_FEASIBILITY"
        return {"feasible": False, "status": st, "reason": msg if msg.startswith((st, "EXCLUDED_FEASIBILITY")) else f"{type(exc).__name__}: {msg}"}


def build_vcc_parts(
    raw_dir: Path,
    hgnc: HGNCLookup,
    line: str,
    min_umi: int = 300,
    min_ngenes: int = 200,
    max_pct_mt: Optional[float] = None,
    min_batch_nt: int = 15,
    min_half_nt: int = 15,
    min_guide_thr: int = 5,
    default_guide_thr: int = 7,
    x_layer: str = "auto",
    context_map: Optional[Dict[str, str]] = None,
    context_map_file: Optional[Path] = None,
    line_manifest: Optional[Dict[str, Any]] = None,
    line_manifest_file: Optional[Path] = None,
    guide_lib_file: Optional[Path] = None,
    seed: int = DEFAULT_R15_SEED,
) -> Tuple[List[PartSpec], Dict[str, Any]]:
    """
    Build standardized `PartSpec` shards for a Round 15 VCC candidate line (`PREREG_H v2 + Addendum A`).
    """
    canon_line = validate_r15_candidate_line(line)
    feas = check_vcc_feasibility(
        raw_dir=raw_dir,
        line=canon_line,
        hgnc=hgnc,
        x_layer=x_layer,
        context_map=context_map,
        context_map_file=context_map_file,
        line_manifest=line_manifest,
        line_manifest_file=line_manifest_file,
        guide_lib_file=guide_lib_file,
    )
    if not feas.get("feasible", False):
        if feas.get("status") == "SLOT_AMBIGUOUS":
            raise RuntimeError(str(feas.get("reason")))
        raise RuntimeError(f"EXCLUDED_FEASIBILITY: {canon_line} failed VCC feasibility gate: {feas.get('reason')}")

    scan = _scan_vcc_inputs(
        raw_dir=raw_dir,
        line=canon_line,
        hgnc=hgnc,
        x_layer=x_layer,
        context_map=context_map,
        context_map_file=context_map_file,
        line_manifest=line_manifest,
        line_manifest_file=line_manifest_file,
        guide_lib_file=guide_lib_file,
        min_umi=min_umi,
        min_ngenes=min_ngenes,
        max_pct_mt=max_pct_mt,
        min_guide_thr=min_guide_thr,
        default_guide_thr=default_guide_thr,
        load_matrices=True,
    )

    thrs: Optional[np.ndarray] = None
    if scan["used_em_fallback"]:
        hist_sum: Optional[np.ndarray] = None
        for rec in scan["file_records"]:
            fh = rec["fit_hist"]
            if fh is not None:
                hist_sum = fh.copy() if hist_sum is None else (hist_sum + fh)
        if hist_sum is None:
            raise RuntimeError("EXCLUDED_FEASIBILITY: guide EM fallback lacked CRISPR guide rows.")
        thrs = fit_poisson_gaussian_thresholds_from_hist(
            hist_sum,
            min_thr=int(min_guide_thr),
            max_thr=20,
            default_thr=int(default_guide_thr),
        )

    all_kept_meta: List[Dict[str, Any]] = []
    per_file_kept: Dict[int, Dict[str, Any]] = {}
    n_rna_qc_excl = n_safe_excl = n_unmapped_excl = n_multi_or_zero_excl = 0
    has_batch_col = all(bool(rec["obs_cols"]["batch_col"]) for rec in scan["file_records"])
    lib_meta = scan["lib_meta"]

    for rec in scan["file_records"]:
        fidx = int(rec["file_idx"])
        n_rna_qc_excl += int(rec["n_rna_qc_excl"])
        qc_abs_rows = rec["qc_abs_rows"]
        keep_rows: List[int] = []
        pmeta: Dict[str, Tuple[str, str, str, str]] = {}

        if rec["has_obs"]:
            tvals, gvals, nvals, cvals, bvals, kvals = (
                rec["tvals"],
                rec["gvals"],
                rec["nvals"],
                rec["cvals"],
                rec["bvals"],
                rec["kvals"],
            )
            for abs_r in qc_abs_rows:
                ri = int(abs_r)
                kf = float(kvals[ri]) if kvals is not None else None
                cat, pkey, meta = _classify_vcc_obs_cell(
                    tgt_raw=tvals[ri],
                    guide_raw=gvals[ri],
                    ntc_raw=nvals[ri],
                    cls_raw=cvals[ri],
                    n_feat=kf,
                    ncount_is_nperts=rec["nc_is_nperts"],
                    hgnc=hgnc,
                )
                if cat == "safe":
                    n_safe_excl += 1
                    continue
                if cat == "unmapped":
                    n_unmapped_excl += 1
                    continue
                if cat == "multi_or_zero" or meta is None:
                    n_multi_or_zero_excl += 1
                    continue
                pmeta[pkey] = meta
                keep_rows.append(ri)
                all_kept_meta.append(
                    {
                        "file_idx": fidx,
                        "local_slot_pos": len(keep_rows) - 1,
                        "raw_row_idx": ri,
                        "barcode": rec["all_barcodes"][ri],
                        "is_nt": (cat == "ntc"),
                        "canon_sym": pkey,
                        "raw_batch": _clean_str(bvals[ri]) if has_batch_col else "",
                    }
                )
        else:
            crp_qc_csc = rec["crp_qc_csc"]
            crp_keys = rec["crp_keys"]
            if crp_qc_csc is not None and len(qc_abs_rows) > 0:
                assert thrs is not None
                top1_v, top2_v, top1_row = extract_top2_csc(crp_qc_csc)
                pass_mask = crp_qc_csc.data >= thrs[crp_qc_csc.indices]
                csum = np.concatenate(([0], np.cumsum(pass_mask, dtype=np.int64)))
                n_pos = csum[crp_qc_csc.indptr[1:]] - csum[crp_qc_csc.indptr[:-1]]
                bvals = rec["bvals"]
                for col_i, abs_r in enumerate(qc_abs_rows):
                    ri = int(abs_r)
                    np_i = int(n_pos[col_i])
                    r_top = int(top1_row[col_i])
                    v1 = float(top1_v[col_i])
                    v2 = float(top2_v[col_i])
                    if np_i != 1 or r_top < 0 or v1 < float(thrs[r_top]) or v1 < float(min_guide_thr) or not (v2 < 3.0 or v1 / max(v2, 0.5) >= 3.0):
                        n_multi_or_zero_excl += 1
                        continue
                    g_key = crp_keys[r_top]
                    g_cat, g_meta = _classify_em_guide_key(g_key, lib_meta, hgnc)
                    if g_cat == "safe":
                        n_safe_excl += 1
                        continue
                    if g_cat == "unmapped" or g_meta is None:
                        n_unmapped_excl += 1
                        continue
                    pkey = g_meta[1]
                    pmeta[pkey] = g_meta
                    keep_rows.append(ri)
                    all_kept_meta.append(
                        {
                            "file_idx": fidx,
                            "local_slot_pos": len(keep_rows) - 1,
                            "raw_row_idx": ri,
                            "barcode": rec["all_barcodes"][ri],
                            "is_nt": (g_cat == "ntc"),
                            "canon_sym": pkey,
                            "raw_batch": _clean_str(bvals[ri]) if has_batch_col else "",
                        }
                    )

        if keep_rows:
            keep_arr = np.asarray(keep_rows, dtype=np.int64)
            with h5py.File(rec["path"], "r") as f:
                x_node = _get_h5ad_slot_node(f, rec["slot"])
                x_kept = _extract_kept_h5ad_csr(
                    x_node=x_node,
                    layout=rec["layout"],
                    keep_arr=keep_arr,
                    col_rows=rec["vmeta"]["rna_rows"],
                    n_obs=rec["n_obs"],
                    n_var=rec["n_var"],
                )
            per_file_kept[fidx] = {
                "rec": rec,
                "keep_arr": keep_arr,
                "x_kept": x_kept,
                "pmeta": pmeta,
                "n_kept": len(keep_arr),
            }

    if len(all_kept_meta) == 0:
        raise RuntimeError(f"EXCLUDED_FEASIBILITY: Zero single-construct QC cells retained for line '{canon_line}'.")

    half_flags, batches_final, split_audit = _assign_two_tier_splits(
        kept_meta=all_kept_meta,
        has_batch_col=has_batch_col,
        min_half_nt=int(min_half_nt),
        min_batch_nt=int(min_batch_nt),
        seed=int(seed),
    )

    # Slice assigned half_flags & batches back onto per-file PartSpecs cleanly
    parts: List[PartSpec] = []
    offset = 0
    file_slots_audit: Dict[str, Any] = {}
    for rec in scan["file_records"]:
        fidx = int(rec["file_idx"])
        entry = per_file_kept.get(fidx)
        if entry is None:
            continue
        nk = int(entry["n_kept"])
        sl = slice(offset, offset + nk)
        offset += nk
        h_part = half_flags[sl]
        b_part = batches_final[sl]
        ka = entry["keep_arr"]
        nt_part = np.array([bool(m["is_nt"]) for m in all_kept_meta[sl]], dtype=bool)
        pkeys_part = [str(m["canon_sym"]) for m in all_kept_meta[sl]]
        umi_part = rec["umi_all"][ka].astype(np.float64)
        vmeta = rec["vmeta"]
        read_chunk_fn, read_rows_fn, read_chunk_cols_fn = _make_csr_part_callbacks(entry["x_kept"], len(vmeta["rna_rows"]))
        parts.append(
            PartSpec(
                tag=f"{canon_line}_part{fidx}",
                path=rec["path"],
                n_cells=nk,
                n_genes_raw=len(vmeta["rna_rows"]),
                nt_mask=nt_part,
                batches=b_part,
                half_flags=h_part,
                umi=umi_part,
                pert_keys=pkeys_part,
                pert_meta=entry["pmeta"],
                var_gene_ids=list(vmeta["r_ids"]),
                var_gene_names=list(vmeta["r_names"]),
                read_chunk_fn=read_chunk_fn,
                read_rows_fn=read_rows_fn,
                read_chunk_cols_fn=read_chunk_cols_fn,
            )
        )
        file_slots_audit[f"file_{fidx}"] = rec["slot_summary"]

    audit: Dict[str, Any] = {
        "line": canon_line,
        "prereg_protocol": "PREREG_H_v2",
        "line_identity_rules": scan["line_identity_rules"],
        "n_jurkat_dropped_pre_qc": int(scan["n_jurkat_dropped"]),
        "n_other_lines_dropped_pre_qc": int(scan["n_other_lines_dropped"]),
        "n_target_line_raw_cells": int(scan["n_line_raw_cells"]),
        "min_umi": int(min_umi),
        "min_ngenes": int(min_ngenes),
        "max_pct_mt": max_pct_mt,
        "min_batch_nt": int(min_batch_nt),
        "min_half_nt": int(min_half_nt),
        "guide_call_mode": "obs_single_construct" if scan["used_obs_calls"] else "h5ad_crispr_em_fallback",
        "median_guide_threshold": float(np.median(thrs)) if thrs is not None else None,
        "n_rna_qc_excluded": int(n_rna_qc_excl),
        "n_safe_targeting_excluded": int(n_safe_excl),
        "n_unmapped_target_excluded": int(n_unmapped_excl),
        "n_multi_or_zero_guide_excluded": int(n_multi_or_zero_excl),
        "n_qc_kept_cells": len(all_kept_meta),
        "n_kept_ntc": int(feas_ntc := sum(1 for m in all_kept_meta if m["is_nt"])),
        "n_kept_targeting": int(len(all_kept_meta) - feas_ntc),
        "n_kept_unique_targets": len({m["canon_sym"] for m in all_kept_meta if not m["is_nt"]}),
        "file_slots": file_slots_audit,
        **split_audit,
    }
    return parts, audit


def prepare_r15_line(
    line: str,
    raw_dir: Path = R15_RAW_DIR,
    out_root: Path = R15_LINES_DIR,
    hgnc: Optional[HGNCLookup] = None,
    min_umi: int = 300,
    min_ngenes: int = 200,
    max_pct_mt: Optional[float] = None,
    min_batch_nt: int = 15,
    min_half_nt: int = 15,
    min_guide_thr: int = 5,
    default_guide_thr: int = 7,
    x_layer: str = "auto",
    context_map: Optional[Dict[str, str]] = None,
    context_map_file: Optional[Path] = None,
    line_manifest: Optional[Dict[str, Any]] = None,
    line_manifest_file: Optional[Path] = None,
    guide_lib_file: Optional[Path] = None,
    seed: int = DEFAULT_R15_SEED,
    workers: int = 4,
    out_name: Optional[str] = None,
) -> Path:
    """Standardize a Round 15 VCC candidate line into `out_root/<line>/` and log full audit."""
    canon_line = validate_r15_candidate_line(line)
    target_out_name = validate_r15_candidate_line(out_name) if out_name else canon_line
    hgnc_l = hgnc or build_hgnc_lookup()
    parts, audit = build_vcc_parts(
        raw_dir=raw_dir,
        hgnc=hgnc_l,
        line=canon_line,
        min_umi=min_umi,
        min_ngenes=min_ngenes,
        max_pct_mt=max_pct_mt,
        min_batch_nt=min_batch_nt,
        min_half_nt=min_half_nt,
        min_guide_thr=min_guide_thr,
        default_guide_thr=default_guide_thr,
        x_layer=x_layer,
        context_map=context_map,
        context_map_file=context_map_file,
        line_manifest=line_manifest,
        line_manifest_file=line_manifest_file,
        guide_lib_file=guide_lib_file,
        seed=seed,
    )
    out_dir = process_r10_line(
        line_name=target_out_name,
        parts=parts,
        hgnc=hgnc_l,
        seed=int(seed),
        out_root=Path(out_root),
        extra_meta=audit,
        n_workers=max(1, int(workers)),
    )
    meta_path = out_dir / "line_metadata.json"
    if meta_path.exists():
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        meta["r15_audit"] = audit
        meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    (out_dir / "r15_prep_audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    return out_dir


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Round 15 (`PREREG_H` v2) VCC `.h5ad` feasibility checker & line builder.")
    p.add_argument("--mode", choices=["check", "prepare"], default="check")
    p.add_argument("--line", required=True, help="Candidate line: cal33, hela, ac16, molm13, u87mg (jurkat rejected).")
    p.add_argument("--raw-dir", type=Path, default=R15_RAW_DIR)
    p.add_argument("--out-root", type=Path, default=R15_LINES_DIR)
    p.add_argument("--x-layer", type=str, default="auto")
    p.add_argument("--context-map-file", type=Path, default=None)
    p.add_argument("--line-manifest-file", type=Path, default=None)
    p.add_argument("--guide-lib-file", type=Path, default=None)
    p.add_argument("--seed", type=int, default=DEFAULT_R15_SEED)
    p.add_argument("--workers", type=int, default=4)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    hgnc = build_hgnc_lookup()
    if args.mode == "check":
        res = check_vcc_feasibility(
            raw_dir=args.raw_dir,
            line=args.line,
            hgnc=hgnc,
            x_layer=args.x_layer,
            context_map_file=args.context_map_file,
            line_manifest_file=args.line_manifest_file,
            guide_lib_file=args.guide_lib_file,
        )
        print(json.dumps(res, indent=2))
        if not res.get("feasible", False):
            sys.exit(2)
        return

    out_dir = prepare_r15_line(
        line=args.line,
        raw_dir=args.raw_dir,
        out_root=args.out_root,
        hgnc=hgnc,
        x_layer=args.x_layer,
        context_map_file=args.context_map_file,
        line_manifest_file=args.line_manifest_file,
        guide_lib_file=args.guide_lib_file,
        seed=args.seed,
        workers=args.workers,
    )
    print(json.dumps({"status": "PREPARED", "line": args.line, "out_dir": str(out_dir)}, indent=2))


if __name__ == "__main__":
    main()
