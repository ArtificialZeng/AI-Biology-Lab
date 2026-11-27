#!/usr/bin/env python3
"""
Clean-Room Round 16 Addendum B Occam Adapter (`A-T3`, strictly aligned to
`PREREG_I_ADDENDUM_B_DRAFT_v2.md` and `PREREG_I_ADDENDUM_B_LOCK.md`).

Design boundaries:
  - Imports ONLY Occam's unchanged frozen model / split-half utilities from
    `round16_addB/shared/models/occam_src_9ddb890e/src/Occam-9ddb890e65ca1ec2f81882c78311002a75b52f67/src`.
  - NEVER reads or imports `round16_addB/src/`, `round16_addB/results/`,
    `round16/src/`, `round16/results/group3/`, or any removed archive member.
  - Strictly Python 3.8.19 compatible (`round16_addB/.venv`).
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import re
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Set, Tuple, Union

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import pearsonr
import torch

AUDIT_DIR = Path(__file__).resolve().parent
ADDB_ROOT = AUDIT_DIR.parent
REPO_ROOT = ADDB_ROOT.parent

DEFAULT_OCCAM_SRC_ROOT = (
    ADDB_ROOT
    / "shared"
    / "models"
    / "occam_src_9ddb890e"
    / "src"
    / "Occam-9ddb890e65ca1ec2f81882c78311002a75b52f67"
    / "src"
)
DEFAULT_HF_ROOT = ADDB_ROOT / "shared" / "models" / "occam_hf_acd080b7"
DEFAULT_CKPT_PATH = (
    DEFAULT_HF_ROOT
    / "checkpoints"
    / "occam_seeds"
    / "rhj_fewshot"
    / "all_sources"
    / "seed_42"
    / "best_model.pth"
)
DEFAULT_BENCH_DIR = DEFAULT_HF_ROOT / "benchmarks" / "rhj_fewshot"
DEFAULT_GENE_AXIS_CSV = DEFAULT_BENCH_DIR / "gene_axis.csv"
DEFAULT_PERTURBATIONS_CSV = DEFAULT_BENCH_DIR / "perturbations.csv"
DEFAULT_COHORT_JSON = DEFAULT_BENCH_DIR / "splithalf" / "cohort.json"
DEFAULT_SPLITS_MANIFEST_JSON = DEFAULT_BENCH_DIR / "splithalf" / "splits_manifest.json"
DEFAULT_SCORES_CSV = DEFAULT_HF_ROOT / "official_comparisons" / "selected_scores_by_seed_split.csv"
DEFAULT_CHECKSUMS_FILE = ADDB_ROOT / "R16_ADDB_CHECKSUMS.sha256"

DEFAULT_NADIG_RAW_DIR = Path("/Users/mac/.cache/r9_p2_raw/nadig2025")
DEFAULT_NADIG_HEPG2_H5AD = DEFAULT_NADIG_RAW_DIR / "GSE264667_hepg2_raw_singlecell_01.h5ad"
DEFAULT_HGNC_FILE = Path("/Users/mac/.cache/r9_p2_raw/hgnc_complete_set.txt")

EXPECTED_BEST_MODEL_SHA256 = "11d5ba5515e79710b3a952532c60aa29515046117292528c59d452a0937685ff"
EXPECTED_GENE_AXIS_CONTENT_SHA256 = "ef64f2ba7ff7816ed947796066cd8a68014884c6176619d750ef79462c78d4c0"
EXPECTED_NADIG_HEPG2_LOCAL_SHA256 = "e1ad7c3c5a201c861a207a858aa7e59f5e6ac1955674c415f7de0d1dadadb52e"

OCCAM_TARGET_SUM = 2000.0
OCCAM_NUM_ENCODER_GENES = 44728
OCCAM_NUM_PERT_NODES = 18216
OCCAM_NUM_RHJ_AXIS_GENES = 6983
CONTROL_LABEL = "non-targeting"


def ensure_occam_src_on_path(src_root: Union[str, Path] = DEFAULT_OCCAM_SRC_ROOT) -> Path:
    """Insert unchanged pinned Occam `src/` root onto `sys.path` if missing."""
    p = Path(src_root).resolve()
    if not p.is_dir():
        raise FileNotFoundError(f"Occam source root not found: {p}")
    sp_str = str(p)
    if sp_str not in sys.path:
        sys.path.insert(0, sp_str)
    return p


def load_hgnc_bridge(hgnc_path: Union[str, Path] = DEFAULT_HGNC_FILE):
    """Load frozen Round 10 `HGNCBridge` cleanly (Python 3.8 compatible)."""
    r10_src = REPO_ROOT / "round10" / "src"
    if str(r10_src) not in sys.path:
        sys.path.insert(0, str(r10_src))
    from data_loader_r10 import HGNCBridge  # type: ignore

    return HGNCBridge.from_tsv(Path(hgnc_path))


def sha256_file(path: Union[str, Path], chunk_bytes: int = 1 << 22) -> str:
    """Stream SHA-256 hex digest of a file."""
    h = hashlib.sha256()
    with open(Path(path), "rb") as fh:
        for block in iter(lambda: fh.read(chunk_bytes), b""):
            h.update(block)
    return h.hexdigest()


def sha256_str_lines(items: Sequence[str]) -> str:
    """Compute Occam `gene_axis_sha256` content digest: `sha256('\\n'.join(items))`."""
    return hashlib.sha256("\n".join(str(x) for x in items).encode("utf-8")).hexdigest()


def sha256_array_bytes(arr: np.ndarray) -> str:
    """Compute Occam split-array digest: `sha256(np.ascontiguousarray(arr).tobytes())`."""
    return hashlib.sha256(np.ascontiguousarray(arr).tobytes()).hexdigest()


def verify_addb_checksums(
    checksums_file: Union[str, Path] = DEFAULT_CHECKSUMS_FILE,
    models_base_dir: Union[str, Path] = ADDB_ROOT / "shared" / "models",
) -> List[Dict[str, str]]:
    """Verify all file checksums listed in `R16_ADDB_CHECKSUMS.sha256`."""
    base = Path(models_base_dir)
    records: List[Dict[str, str]] = []
    with open(Path(checksums_file), "r", encoding="utf-8") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line or line.startswith("#") or "(tree," in line:
                continue
            parts = line.split(None, 1)
            if len(parts) != 2:
                continue
            want_sha, rel_path = parts[0].strip(), parts[1].strip()
            fpath = base / rel_path
            if not fpath.is_file():
                raise FileNotFoundError(f"Checksum target missing: {fpath}")
            got_sha = sha256_file(fpath)
            if got_sha != want_sha:
                raise ValueError(
                    f"SHA-256 mismatch for {rel_path}: got {got_sha}, expected {want_sha}"
                )
            records.append({"rel_path": rel_path, "sha256": got_sha, "status": "OK"})
    if len(records) < 8:
        raise RuntimeError(f"Expected >= 8 checksum entries, verified {len(records)}")
    return records


def decode_str_array(arr: Any) -> List[str]:
    """Decode HDF5 byte/string array cleanly into Python strings."""
    return [
        x.decode("utf-8", errors="replace") if isinstance(x, (bytes, np.bytes_)) else str(x)
        for x in arr
    ]


def read_obs_or_var_col_legacy(
    grp: h5py.Group,
    col: str,
    indices: Optional[np.ndarray] = None,
) -> List[str]:
    """
    Decode an AnnData HDF5 `obs`/`var` column supporting:
      1. Modern anndata categorical groups (`grp[col]/categories` + `codes`)
      2. Nullable groups (`grp[col]/values` + optional `mask`)
      3. Legacy anndata `__categories` (`grp['__categories'][col]` + integer dataset `grp[col]`)
      4. Flat datasets (`bytes` / `str` / numeric)

    Negative category codes (`< 0`) and masked entries map to empty string `""`.
    """
    node = grp[col]
    idx = None if indices is None else np.asarray(indices, dtype=np.int64)

    if isinstance(node, h5py.Group):
        enc = node.attrs.get("encoding-type")
        if isinstance(enc, (bytes, np.bytes_)):
            enc = enc.decode("utf-8", errors="replace")
        if "categories" in node and "codes" in node:
            cats = decode_str_array(node["categories"][:])
            codes = node["codes"][:] if idx is None else node["codes"][idx]
            return [cats[int(c)] if int(c) >= 0 else "" for c in codes]
        if "values" in node:
            vals = decode_str_array(node["values"][:] if idx is None else node["values"][idx])
            if "mask" in node:
                masks = node["mask"][:] if idx is None else node["mask"][idx]
                return ["" if bool(m) else v for v, m in zip(vals, masks)]
            return vals
        raise ValueError(f"Unsupported HDF5 group encoding for {col!r} (keys={list(node.keys())})")

    if "__categories" in grp and col in grp["__categories"] and np.issubdtype(node.dtype, np.integer):
        cats = decode_str_array(grp["__categories"][col][:])
        codes = node[:] if idx is None else node[idx]
        return [cats[int(c)] if int(c) >= 0 else "" for c in codes]

    vals_raw = node[:] if idx is None else node[idx]
    return decode_str_array(vals_raw)


def strip_ensembl_version(gid: str) -> str:
    """Strip `.N` Ensembl version suffix (e.g., `ENSG00000123456.12` -> `ENSG00000123456`)."""
    s = str(gid).strip()
    if s.startswith("ENSG") and "." in s:
        return s.split(".", 1)[0]
    return s


def dedup_leading_ensembl_positions(raw_ensembl_ids: Sequence[str]) -> Dict[str, int]:
    """
    Strip Ensembl versions and keep ONLY the leading occurrence in raw file column order
    (`PREREG_I_ADDENDUM_B_DRAFT_v2.md §3.3.2`).
    """
    pos: Dict[str, int] = {}
    for col_idx, raw_id in enumerate(raw_ensembl_ids):
        clean_id = strip_ensembl_version(raw_id)
        if not clean_id or clean_id in pos:
            continue
        pos[clean_id] = col_idx
    return pos


def load_rhj_gene_axis(
    gene_axis_csv: Union[str, Path] = DEFAULT_GENE_AXIS_CSV,
    expected_content_sha256: str = EXPECTED_GENE_AXIS_CONTENT_SHA256,
) -> Tuple[List[str], List[str]]:
    """Load the `6,983`-gene `rhj_fewshot` Ensembl + symbol axis and assert exact content digest."""
    axis_ens: List[str] = []
    axis_sym: List[str] = []
    with open(Path(gene_axis_csv), "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            eid = strip_ensembl_version(row["ensembl_gene_id"])
            sym = str(row["gene_name"]).strip()
            if eid:
                axis_ens.append(eid)
                axis_sym.append(sym)
    if len(axis_ens) != OCCAM_NUM_RHJ_AXIS_GENES:
        raise ValueError(
            f"Expected {OCCAM_NUM_RHJ_AXIS_GENES} RHJ axis genes, got {len(axis_ens)}"
        )
    got_sha = sha256_str_lines(axis_ens)
    if expected_content_sha256 and got_sha != expected_content_sha256:
        raise ValueError(
            f"gene_axis_sha256 mismatch: got {got_sha}, expected {expected_content_sha256}"
        )
    return axis_ens, axis_sym


def load_rhj_eligible_perturbations(
    perturbations_csv: Union[str, Path] = DEFAULT_PERTURBATIONS_CSV,
) -> List[str]:
    """Load sorted list of eligible perturbation keys (`eligible == 1`) from `perturbations.csv`."""
    eligible: List[str] = []
    with open(Path(perturbations_csv), "r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if int(row.get("eligible", "0")) == 1:
                eligible.append(str(row["gene_name"]).strip())
    if len(eligible) != 2381:
        raise ValueError(f"Expected 2,381 eligible rhj_fewshot keys, got {len(eligible)}")
    return eligible


class TwoTierPertResolver:
    """
    Two-Tier `pert2id` resolver (`PREREG_I_ADDENDUM_B_DRAFT_v2.md §3.5`):
      - Tier 1: Canonicalize the `2,381` eligible `rhj_fewshot` trained keys via `HGNCBridge`,
        mapping e.g. `AARS1 -> AARS`, `DARS1 -> DARS`, `TAFAZZIN -> TAZ`, `H2AZ1 -> H2AFZ`.
      - Tier 2: Direct lookup in `ckpt['pert2id']` when absent from Tier 1 (e.g. `RBM14`).
    """

    def __init__(
        self,
        ckpt_pert2id: Mapping[str, int],
        perturbations_csv: Union[str, Path] = DEFAULT_PERTURBATIONS_CSV,
        hgnc_bridge: Optional[Any] = None,
    ) -> None:
        self.ckpt_pert2id: Dict[str, int] = {str(k): int(v) for k, v in ckpt_pert2id.items()}
        self.hgnc = hgnc_bridge if hgnc_bridge is not None else load_hgnc_bridge()
        self.eligible_keys = load_rhj_eligible_perturbations(perturbations_csv)

        self.tier1_canon_to_key: Dict[str, str] = {}
        collisions: List[Tuple[str, str, str]] = []
        for k in self.eligible_keys:
            if k not in self.ckpt_pert2id:
                raise KeyError(f"Eligible perturbation {k!r} missing from ckpt['pert2id']")
            canon = self.hgnc.canonicalize(k) or k
            if canon in self.tier1_canon_to_key and self.tier1_canon_to_key[canon] != k:
                collisions.append((canon, self.tier1_canon_to_key[canon], k))
            else:
                self.tier1_canon_to_key[canon] = k
        if collisions:
            raise ValueError(f"Unexpected canonical collisions in rhj_2381: {collisions[:5]}")

    def resolve(self, target_symbol: str) -> Tuple[str, int, str]:
        """Resolve one evaluation target symbol -> `(resolved_ckpt_key, pert_id, tier)`."""
        raw = str(target_symbol).strip()
        canon = self.hgnc.canonicalize(raw) or raw
        if canon in self.tier1_canon_to_key:
            key = self.tier1_canon_to_key[canon]
            return key, self.ckpt_pert2id[key], "tier1_rhj2381"
        if raw in self.ckpt_pert2id:
            return raw, self.ckpt_pert2id[raw], "tier2_direct"
        if canon in self.ckpt_pert2id:
            return canon, self.ckpt_pert2id[canon], "tier2_direct_canon"
        raise KeyError(f"Target symbol {target_symbol!r} (canon={canon!r}) not resolvable in pert2id")

    def resolve_many(
        self, target_symbols: Sequence[str]
    ) -> Tuple[List[int], List[Dict[str, Any]]]:
        """Resolve an ordered list of target symbols -> `(pert_ids, audit_rows)`."""
        ids: List[int] = []
        rows: List[Dict[str, Any]] = []
        for sym in target_symbols:
            key, pid, tier = self.resolve(sym)
            ids.append(pid)
            rows.append(
                {
                    "eval_target": str(sym),
                    "resolved_key": key,
                    "pert_id": int(pid),
                    "tier": tier,
                    "bridged": bool(key != str(sym)),
                }
            )
        return ids, rows


def resolve_output_projection(
    eval_genes: Sequence[str],
    rhj_axis_ensembl: Sequence[str],
    rhj_axis_symbols: Sequence[str],
    encoder_gene_names: Sequence[str],
    hgnc_bridge: Optional[Any] = None,
) -> Dict[str, Any]:
    """
    Canonicalize `6,983` RHJ axis outputs via `HGNCBridge` (`0` collisions required) and align
    ordered `eval_genes` (`G_eval^R16`) onto `44,728`-dim decoder outputs and `6,983`-dim axis.
    """
    hb = hgnc_bridge if hgnc_bridge is not None else load_hgnc_bridge()
    enc_idx = {str(g): i for i, g in enumerate(encoder_gene_names)}

    canon_to_pos: Dict[str, int] = {}
    collisions: List[Tuple[str, str, str]] = []
    for pos, (ens, sym) in enumerate(zip(rhj_axis_ensembl, rhj_axis_symbols)):
        c_sym = hb.canonicalize(sym) if sym else ""
        c_ens = hb.canonicalize(ens) if (not c_sym or c_sym == sym) else c_sym
        canon = c_sym or c_ens or ens
        if canon in canon_to_pos:
            collisions.append((canon, rhj_axis_symbols[canon_to_pos[canon]], sym))
        else:
            canon_to_pos[canon] = pos
    if collisions:
        raise ValueError(f"HGNC collision on RHJ_6983 outputs: {collisions[:5]}")

    rhj_cols: List[int] = []
    enc_cols: List[int] = []
    missing: List[str] = []
    for g in eval_genes:
        g_clean = hb.canonicalize(str(g).strip()) or str(g).strip()
        pos = canon_to_pos.get(g_clean)
        if pos is None:
            missing.append(str(g))
            continue
        eid = rhj_axis_ensembl[pos]
        ec = enc_idx.get(eid)
        if ec is None:
            missing.append(str(g))
            continue
        rhj_cols.append(pos)
        enc_cols.append(ec)

    return {
        "eval_genes": list(eval_genes),
        "n_eval_genes": len(eval_genes),
        "n_mapped": len(enc_cols),
        "missing_genes": missing,
        "rhj_cols": np.asarray(rhj_cols, dtype=np.int64),
        "enc_cols": np.asarray(enc_cols, dtype=np.int64),
    }


def compute_pairwise_common_rhj_mask(
    raw_ensembl_a: Sequence[str],
    raw_ensembl_b: Sequence[str],
    rhj_axis_ensembl: Sequence[str],
) -> np.ndarray:
    """Boolean `(6983,)` mask retaining RHJ axis genes present in BOTH raw gene lists (`§3.3.4`)."""
    set_a = set(dedup_leading_ensembl_positions(raw_ensembl_a).keys())
    set_b = set(dedup_leading_ensembl_positions(raw_ensembl_b).keys())
    common = set_a & set_b
    return np.asarray([g in common for g in rhj_axis_ensembl], dtype=bool)


def normalize_counts_to_rhj_44728(
    raw_counts: np.ndarray,
    raw_ensembl_ids: Sequence[str],
    rhj_axis_ensembl: Sequence[str],
    encoder_gene_names: Sequence[str],
    target_sum: float = OCCAM_TARGET_SUM,
    active_rhj_mask: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    """
    Preprocess raw single-cell counts (`N_cells, G_raw`) into Occam's `44,728`-gene input space
    under the `rhj_6983_axis` (or `pairwise_common_axis` when `active_rhj_mask` is passed):
      1. Full-panel library size `L_i = sum_g raw_counts[i, g]` in float64 across ALL raw columns
         (never restricting before summing). `L_i == 0` is safe-guarded to `1.0` (`-> 0.0`).
      2. Deduplicate raw Ensembl IDs keeping the leading occurrence (`§3.3.2`).
      3. For matched (and active) RHJ axis genes: `log1p(raw_val / float32(L_i) * float(target_sum))`.
      4. Scatter strictly into `(N_cells, 44728)` float32 tensor (and `(N_cells, 6983)` float32).
    """
    x = np.asarray(raw_counts, dtype=np.float32)
    if x.ndim != 2 or x.shape[1] != len(raw_ensembl_ids):
        raise ValueError(
            f"Shape mismatch: raw_counts={x.shape}, len(raw_ensembl_ids)={len(raw_ensembl_ids)}"
        )
    if active_rhj_mask is not None:
        mask_arr = np.asarray(active_rhj_mask, dtype=bool)
        if mask_arr.shape != (len(rhj_axis_ensembl),):
            raise ValueError(
                f"active_rhj_mask shape {mask_arr.shape} != ({len(rhj_axis_ensembl)},)"
            )
    else:
        mask_arr = np.ones(len(rhj_axis_ensembl), dtype=bool)

    raw_pos = dedup_leading_ensembl_positions(raw_ensembl_ids)
    enc_pos = {str(g): i for i, g in enumerate(encoder_gene_names)}

    raw_cols: List[int] = []
    rhj_cols: List[int] = []
    enc_cols: List[int] = []
    absent_raw: List[str] = []
    masked_out: List[str] = []

    for r_i, eid in enumerate(rhj_axis_ensembl):
        if eid not in raw_pos:
            absent_raw.append(eid)
            continue
        if not mask_arr[r_i]:
            masked_out.append(eid)
            continue
        ec = enc_pos.get(eid)
        if ec is None:
            raise KeyError(f"RHJ axis gene {eid!r} absent from encoder gene_names")
        raw_cols.append(raw_pos[eid])
        rhj_cols.append(r_i)
        enc_cols.append(ec)

    row_sums = x.sum(axis=1, dtype=np.float64)
    row_sums = np.where(row_sums > 0.0, row_sums, 1.0).astype(np.float32)

    X_6983 = np.zeros((x.shape[0], len(rhj_axis_ensembl)), dtype=np.float32)
    if raw_cols:
        sub = x[:, np.asarray(raw_cols, dtype=np.int64)]
        normed = np.log1p(sub / row_sums[:, None] * np.float32(target_sum)).astype(np.float32)
        X_6983[:, np.asarray(rhj_cols, dtype=np.int64)] = normed

    X_44728 = np.zeros((x.shape[0], len(encoder_gene_names)), dtype=np.float32)
    if enc_cols:
        X_44728[:, np.asarray(enc_cols, dtype=np.int64)] = X_6983[
            :, np.asarray(rhj_cols, dtype=np.int64)
        ]

    audit = {
        "n_cells": int(x.shape[0]),
        "n_raw_genes": int(x.shape[1]),
        "n_rhj_present_in_raw": int(len(rhj_axis_ensembl) - len(absent_raw)),
        "n_rhj_absent_in_raw": int(len(absent_raw)),
        "n_rhj_masked_out": int(len(masked_out)),
        "n_active_injected": int(len(enc_cols)),
        "absent_sample": absent_raw[:5],
    }
    return X_44728, X_6983, audit


def normalize_counts_full_44728_api(
    raw_counts: np.ndarray,
    raw_ensembl_ids: Sequence[str],
    encoder_gene_names: Sequence[str],
    target_sum: float = OCCAM_TARGET_SUM,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    Occam `GeneSpace.prepare` API convention (`PREREG_I_ADDENDUM_B_DRAFT_v2.md §3.3.5`):
    scatter all matched `44,728`-vocabulary columns first, then `NormalizePerCell(target_sum) -> Log1p`
    over matched encoder columns.
    """
    x = np.asarray(raw_counts, dtype=np.float32)
    raw_pos = dedup_leading_ensembl_positions(raw_ensembl_ids)
    enc_pos = {str(g): i for i, g in enumerate(encoder_gene_names)}

    raw_cols: List[int] = []
    enc_cols: List[int] = []
    for eid, ec in enc_pos.items():
        if eid in raw_pos:
            raw_cols.append(raw_pos[eid])
            enc_cols.append(ec)

    X_44728 = np.zeros((x.shape[0], len(encoder_gene_names)), dtype=np.float32)
    if raw_cols:
        sub = x[:, np.asarray(raw_cols, dtype=np.int64)]
        rs = sub.sum(axis=1, dtype=np.float64)
        rs = np.where(rs > 0.0, rs, 1.0).astype(np.float32)
        X_44728[:, np.asarray(enc_cols, dtype=np.int64)] = np.log1p(
            sub / rs[:, None] * np.float32(target_sum)
        ).astype(np.float32)

    return X_44728, {"n_cells": int(x.shape[0]), "n_matched_44728": int(len(enc_cols))}


def count_state_dict_params(state_dict: Mapping[str, torch.Tensor]) -> int:
    """Count total scalar parameters across tensors in a state dict."""
    total = 0
    for tensor in state_dict.values():
        n = 1
        for d in tuple(tensor.shape):
            n *= int(d)
        total += n
    return total


def verify_checkpoint_contract(
    ckpt: Mapping[str, Any],
    model: torch.nn.Module,
    ckpt_path: Optional[Union[str, Path]] = None,
    gene_axis_csv: Union[str, Path] = DEFAULT_GENE_AXIS_CSV,
    verify_sha256: bool = True,
) -> Dict[str, Any]:
    """
    Verify all structural assertions in `PREREG_I_ADDENDUM_B_DRAFT_v2.md §5.1.1`.
    Raises `AssertionError` / `ValueError` on any discrepancy.
    """
    sha = None
    if ckpt_path is not None and verify_sha256:
        sha = sha256_file(ckpt_path)
        if sha != EXPECTED_BEST_MODEL_SHA256:
            raise AssertionError(
                f"best_model.pth SHA-256 mismatch: {sha} != {EXPECTED_BEST_MODEL_SHA256}"
            )

    if ckpt.get("is_virtual_cell_checkpoint") is not True:
        raise AssertionError(
            f"is_virtual_cell_checkpoint expected True, got {ckpt.get('is_virtual_cell_checkpoint')!r}"
        )
    target_sum = float(ckpt.get("target_sum", -1.0))
    if target_sum != OCCAM_TARGET_SUM:
        raise AssertionError(f"target_sum expected {OCCAM_TARGET_SUM}, got {target_sum!r}")

    dc = dict(ckpt.get("decoder_config") or {})
    if dc.get("output_activation") != "softplus":
        raise AssertionError(f"output_activation expected 'softplus', got {dc.get('output_activation')!r}")

    ic = dict(ckpt.get("interaction_config") or {})
    checks_ic = [
        ("interaction_type", ic.get("interaction_type"), "film"),
        ("film_mode", ic.get("film_mode"), "multiplicative_only"),
        ("film_mean_free", bool(ic.get("film_mean_free")), True),
        ("film_hidden_dim", int(ic.get("film_hidden_dim", -1)), 256),
    ]
    for name, got, want in checks_ic:
        if got != want:
            raise AssertionError(f"interaction_config.{name}: expected {want!r}, got {got!r}")

    gene_names = list(ckpt.get("gene_names") or [])
    pert2id = dict(ckpt.get("pert2id") or {})
    if len(gene_names) != OCCAM_NUM_ENCODER_GENES:
        raise AssertionError(f"len(gene_names) expected {OCCAM_NUM_ENCODER_GENES}, got {len(gene_names)}")
    if len(pert2id) != OCCAM_NUM_PERT_NODES:
        raise AssertionError(f"len(pert2id) expected {OCCAM_NUM_PERT_NODES}, got {len(pert2id)}")

    gti = ckpt.get("gene_to_index") or {}
    if len(gti) != len(gene_names) or any(gti.get(g) != i for i, g in enumerate(gene_names)):
        raise AssertionError("ckpt['gene_to_index'] is not exact inverse of ckpt['gene_names']")

    axis_ens, _ = load_rhj_gene_axis(
        gene_axis_csv=gene_axis_csv,
        expected_content_sha256=EXPECTED_GENE_AXIS_CONTENT_SHA256,
    )
    missing_axis = [g for g in axis_ens if g not in gti]
    if missing_axis:
        raise AssertionError(f"{len(missing_axis)} RHJ axis genes missing from encoder gene_names")

    if getattr(model, "source_emb", "SENTINEL") is not None:
        raise AssertionError("Expected model.source_emb is None")
    if getattr(model, "modality_emb", "SENTINEL") is not None:
        raise AssertionError("Expected model.modality_emb is None")
    if hasattr(model, "batch_encoder"):
        raise AssertionError("Unexpected batch_encoder attribute on model")
    batch_named = [n for n, _ in model.named_parameters() if "batch_encoder" in n]
    if batch_named:
        raise AssertionError(f"Unexpected batch_encoder parameters: {batch_named}")

    # Strict re-verification of loaded module state_dicts
    esd = ckpt["encoder_state_dict"]
    dsd = ckpt["decoder_state_dict"]
    gsd = ckpt["gnn_state_dict"]
    isd = ckpt["interaction_state_dict"]
    model.decoder.load_state_dict(dsd, strict=True)
    model.gnn.load_state_dict(gsd, strict=True)
    model.interaction.load_state_dict(isd, strict=True)

    return {
        "ckpt_sha256": sha,
        "gene_axis_content_sha256": EXPECTED_GENE_AXIS_CONTENT_SHA256,
        "is_virtual_cell_checkpoint": True,
        "target_sum": target_sum,
        "output_activation": dc.get("output_activation"),
        "interaction_type": ic.get("interaction_type"),
        "film_mode": ic.get("film_mode"),
        "film_mean_free": bool(ic.get("film_mean_free")),
        "film_hidden_dim": int(ic.get("film_hidden_dim")),
        "n_encoder_genes": len(gene_names),
        "n_pert_nodes": len(pert2id),
        "source_emb_is_none": True,
        "modality_emb_is_none": True,
        "n_batch_encoder_params": 0,
        "encoder_params": count_state_dict_params(esd),
        "decoder_params": count_state_dict_params(dsd),
        "gnn_params": count_state_dict_params(gsd),
        "interaction_params": count_state_dict_params(isd),
        "epoch": int(ckpt.get("epoch", -1)),
        "strict_reload_ok": True,
    }


class OccamCleanRoomAdapter:
    """
    Independent clean-room inference wrapper around unchanged Occam `VirtualCellModel`.
    """

    def __init__(
        self,
        ckpt_path: Union[str, Path] = DEFAULT_CKPT_PATH,
        device: str = "cpu",
        num_threads: int = 4,
        occam_src_root: Union[str, Path] = DEFAULT_OCCAM_SRC_ROOT,
        gene_axis_csv: Union[str, Path] = DEFAULT_GENE_AXIS_CSV,
        verify_sha256: bool = True,
    ) -> None:
        self.device = str(device).strip().lower()
        if self.device not in ("cpu", "mps"):
            raise ValueError(f"Unsupported device={device!r}; expected 'cpu' or 'mps'")
        if self.device == "cpu":
            torch.set_num_threads(int(num_threads))
        elif not torch.backends.mps.is_available():
            raise RuntimeError("MPS requested but torch.backends.mps.is_available() is False")

        ensure_occam_src_on_path(occam_src_root)
        from occam.vcell.pipeline import VirtualCellTrainer  # type: ignore

        self.ckpt_path = Path(ckpt_path)
        self.trainer, self.ckpt = VirtualCellTrainer.from_vcell_checkpoint(
            str(self.ckpt_path),
            device=self.device,
        )
        self.model = self.trainer.model
        self.model.eval()
        # Invalidate and prime GNN eval cache cleanly.
        if hasattr(self.model, "invalidate_gnn_cache"):
            self.model.invalidate_gnn_cache()

        self.gene_names: List[str] = list(self.ckpt["gene_names"])
        self.gene_to_index: Dict[str, int] = {str(g): int(i) for i, g in enumerate(self.gene_names)}
        self.pert2id: Dict[str, int] = {str(k): int(v) for k, v in self.ckpt["pert2id"].items()}
        self.target_sum: float = float(self.ckpt["target_sum"])

        self.contract_summary = verify_checkpoint_contract(
            self.ckpt,
            self.model,
            ckpt_path=self.ckpt_path,
            gene_axis_csv=gene_axis_csv,
            verify_sha256=verify_sha256,
        )
        self.model.eval()

    def encode_controls(
        self,
        X_ctrl_44728: np.ndarray,
        batch_size: int = 256,
    ) -> np.ndarray:
        """Encode control cells `(N, 44728) -> (N, 512)` deterministically with `shuffle=False`."""
        x = np.asarray(X_ctrl_44728, dtype=np.float32)
        if x.ndim != 2 or x.shape[1] != OCCAM_NUM_ENCODER_GENES:
            raise ValueError(f"Expected (N, {OCCAM_NUM_ENCODER_GENES}), got {x.shape}")
        n = x.shape[0]
        hidden = int(getattr(self.model, "encoder_hidden_size", 512))
        out = np.empty((n, hidden), dtype=np.float32)
        self.model.eval()
        with torch.no_grad():
            for s in range(0, n, int(batch_size)):
                e = min(s + int(batch_size), n)
                xt = torch.from_numpy(np.ascontiguousarray(x[s:e])).to(self.device)
                zt = self.model.encoder(xt).float().cpu().numpy()
                out[s:e] = zt
        return out

    def predict_per_cell(
        self,
        inputs: np.ndarray,
        pert_indices: Sequence[Sequence[int]],
        batch_size: int = 256,
        out_cols: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """
        Run deterministic per-cell `VirtualCellModel` forward (`shuffle=False`).
        `inputs` may be either `(N, 44728)` expression or `(N, 512)` pre-encoded `z_ctrl`.
        """
        arr = np.asarray(inputs, dtype=np.float32)
        n = arr.shape[0]
        if n != len(pert_indices):
            raise ValueError(f"Row mismatch: inputs={n}, pert_indices={len(pert_indices)}")
        cols = None if out_cols is None else np.asarray(out_cols, dtype=np.int64)
        n_out = OCCAM_NUM_ENCODER_GENES if cols is None else int(cols.size)
        out = np.empty((n, n_out), dtype=np.float32)
        self.model.eval()
        with torch.no_grad():
            for s in range(0, n, int(batch_size)):
                e = min(s + int(batch_size), n)
                xt = torch.from_numpy(np.ascontiguousarray(arr[s:e])).to(self.device)
                pidx = [list(p) for p in pert_indices[s:e]]
                pred = self.model(xt, pidx)["prediction"].float()
                if cols is not None:
                    pred_np = pred.cpu().numpy()[:, cols]
                else:
                    pred_np = pred.cpu().numpy()
                out[s:e] = pred_np
        return out

    def predict_reconstructed_baseline(
        self,
        X_ctrl_44728_or_z: np.ndarray,
        batch_size: int = 256,
        out_cols: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """Predict null/unperturbed reference (`pert_indices = [[-1]]`) per control cell."""
        n = X_ctrl_44728_or_z.shape[0]
        return self.predict_per_cell(
            X_ctrl_44728_or_z,
            [[-1]] * n,
            batch_size=batch_size,
            out_cols=out_cols,
        )

    def predict_deltas(
        self,
        X_ctrl_44728: np.ndarray,
        pert_ids: Sequence[int],
        output_col_idx: np.ndarray,
        baseline: str = "reconstructed",
        batch_size: int = 256,
        z_ctrl: Optional[np.ndarray] = None,
    ) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
        """
        Compute pseudobulk delta matrix `(N_targets, G_out)` over `N_ctrl` control cells:
          `delta[p] = mean_i pred_i(ctrl_i, p) - baseline_mean`
        matching `PREREG_I_ADDENDUM_B_DRAFT_v2.md §3.1 & §3.4`.
        """
        if baseline not in ("reconstructed", "real"):
            raise ValueError(f"Unsupported baseline={baseline!r}; expected 'reconstructed' or 'real'")
        out_cols = np.asarray(output_col_idx, dtype=np.int64)
        x = np.asarray(X_ctrl_44728, dtype=np.float32)
        z = z_ctrl if z_ctrl is not None else self.encode_controls(x, batch_size=batch_size)

        if baseline == "reconstructed":
            base_cells = self.predict_reconstructed_baseline(z, batch_size=batch_size, out_cols=out_cols)
            base_mean = base_cells.mean(axis=0, dtype=np.float64)
        else:
            base_mean = x[:, out_cols].mean(axis=0, dtype=np.float64)

        deltas = np.empty((len(pert_ids), out_cols.size), dtype=np.float32)
        n_cells = x.shape[0]
        nonfinite_count = 0
        for row_i, pid in enumerate(pert_ids):
            pred_cells = self.predict_per_cell(
                z,
                [[int(pid)]] * n_cells,
                batch_size=batch_size,
                out_cols=out_cols,
            )
            bad = ~np.isfinite(pred_cells)
            if bad.any():
                nonfinite_count += int(bad.sum())
                ctrl_sub = x[:, out_cols]
                pred_cells[bad] = ctrl_sub[bad]
            deltas[row_i] = (pred_cells.mean(axis=0, dtype=np.float64) - base_mean).astype(np.float32)

        return deltas, z, {
            "baseline": baseline,
            "n_ctrl_cells": int(n_cells),
            "n_targets": int(len(pert_ids)),
            "n_output_genes": int(out_cols.size),
            "nonfinite_values_replaced": int(nonfinite_count),
        }


def save_delta_bundle_npz(
    out_path: Union[str, Path],
    delta: np.ndarray,
    z_ctrl: np.ndarray,
    target_genes: Sequence[str],
    eval_genes: Sequence[str],
    control_barcodes: Sequence[str],
    meta: Mapping[str, Any],
) -> Path:
    """Export `.npz` audit bundle matching `PREREG_I_ADDENDUM_B_DRAFT_v2.md §6.2`."""
    p = Path(out_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        p,
        delta=np.asarray(delta, dtype=np.float32),
        z_ctrl=np.asarray(z_ctrl, dtype=np.float32),
        target_genes=np.asarray(list(target_genes), dtype=object),
        eval_genes=np.asarray(list(eval_genes), dtype=object),
        control_barcodes=np.asarray(list(control_barcodes), dtype=object),
        meta_json=np.asarray([json.dumps(dict(meta), sort_keys=True)], dtype=object),
    )
    return p


def load_nadig_hepg2_panel_raw(
    h5ad_path: Union[str, Path] = DEFAULT_NADIG_HEPG2_H5AD,
    cohort_perts: Sequence[str] = (),
    expected_sha256: Optional[str] = EXPECTED_NADIG_HEPG2_LOCAL_SHA256,
) -> Dict[str, Any]:
    """
    Read local `GSE264667_hepg2_raw_singlecell_01.h5ad` (`e1ad7c3c...`) cleanly using HDF5:
      - Verifies optional file SHA-256 (`expected_sha256`).
      - Decodes legacy categorical `obs["gene"]`, `obs["gem_group"]`, `var["gene_id"]`.
      - Filters to the union of `cohort_perts` and `'non-targeting'` (`42,402` rows for the `611`-pert cohort).
      - Reads ONLY those selected rows from dense/CSR `/X` in bounded chunks (`float32`),
        computing full-panel row sums over all `9,624` genes (`-> float64`), BEFORE restricting
        to caller axis (`build_splithalf_panel.py:97-154`).
    """
    p = Path(h5ad_path)
    if not p.is_file():
        raise FileNotFoundError(f"Missing Nadig HepG2 raw h5ad: {p}")
    got_sha = sha256_file(p)
    if expected_sha256 and got_sha != expected_sha256:
        raise ValueError(f"Nadig HepG2 SHA-256 mismatch: {got_sha} != {expected_sha256}")

    cohort_set = set(str(x) for x in cohort_perts) | {CONTROL_LABEL}
    with h5py.File(str(p), "r") as f:
        obs_genes = np.asarray(read_obs_or_var_col_legacy(f["obs"], "gene"), dtype=str)
        obs_batches = np.asarray(read_obs_or_var_col_legacy(f["obs"], "gem_group"), dtype=str)
        var_ens = [strip_ensembl_version(g) for g in read_obs_or_var_col_legacy(f["var"], "gene_id")]

        keep_rows = np.flatnonzero(np.isin(obs_genes, list(cohort_set)))
        lp = obs_genes[keep_rows]
        bp = obs_batches[keep_rows]

        X_node = f["X"]
        n_vars = len(var_ens)
        raw_keep = np.empty((keep_rows.size, n_vars), dtype=np.float32)

        if isinstance(X_node, h5py.Dataset):
            # Stream contiguous HDF5 slices and pick matching panel rows without loading full 5.6 GB X.
            total_obs = int(X_node.shape[0])
            chunk = 8192
            write_pos = 0
            for s in range(0, total_obs, chunk):
                e = min(s + chunk, total_obs)
                lo = int(np.searchsorted(keep_rows, s, side="left"))
                hi = int(np.searchsorted(keep_rows, e, side="left"))
                if hi > lo:
                    block = np.asarray(X_node[s:e, :], dtype=np.float32)
                    local_idx = keep_rows[lo:hi] - s
                    raw_keep[write_pos : write_pos + (hi - lo)] = block[local_idx]
                    write_pos += hi - lo
            assert write_pos == keep_rows.size
        else:
            indptr = np.asarray(X_node["indptr"][:], dtype=np.int64)
            data = X_node["data"]
            indices = X_node["indices"]
            for i, r in enumerate(keep_rows):
                p0, p1 = int(indptr[int(r)]), int(indptr[int(r) + 1])
                row_dense = np.zeros(n_vars, dtype=np.float32)
                if p1 > p0:
                    row_dense[indices[p0:p1]] = np.asarray(data[p0:p1], dtype=np.float32)
                raw_keep[i] = row_dense

    return {
        "file_sha256": got_sha,
        "n_total_obs": int(obs_genes.size),
        "n_total_vars": int(len(var_ens)),
        "keep_rows": keep_rows,
        "raw_counts": raw_keep,
        "labels": lp,
        "batches": bp,
        "var_ensembl": var_ens,
    }


def compute_celleval_pseudobulk_metrics(
    real_pert_X: np.ndarray,
    real_pert_labels: np.ndarray,
    real_ctrl_X: np.ndarray,
    pred_pert_X: np.ndarray,
    pred_pert_labels: np.ndarray,
    pred_ctrl_X: np.ndarray,
    cohort_perts: Sequence[str],
) -> Dict[str, Any]:
    """
    Exact numpy/scipy equivalent of `cell-eval`'s pseudobulk `pearson_delta` and `mae`
    (`cell_eval._types._anndata` + `cell_eval.metrics._anndata`):
      - `delta_real(p) = mean(real_pert[p]) - mean(real_ctrl)`
      - `delta_pred(p) = mean(pred_pert[p]) - mean(pred_ctrl)`
      - `pearson_delta = mean_p pearsonr(delta_pred(p), delta_real(p))`
      - `mae = mean_p mean_g |mean(pred_pert[p]) - mean(real_pert[p])|`
    """
    ctrl_real = np.asarray(real_ctrl_X, dtype=np.float64).mean(axis=0)
    ctrl_pred = np.asarray(pred_ctrl_X, dtype=np.float64).mean(axis=0)

    r_map: Dict[str, np.ndarray] = {}
    p_map: Dict[str, np.ndarray] = {}
    for idx_map, arr, lbls in (
        (r_map, np.asarray(real_pert_X, dtype=np.float64), np.asarray(real_pert_labels, dtype=str)),
        (p_map, np.asarray(pred_pert_X, dtype=np.float64), np.asarray(pred_pert_labels, dtype=str)),
    ):
        uniq, inv = np.unique(lbls, return_inverse=True)
        for u_i, p_name in enumerate(uniq):
            if str(p_name) == CONTROL_LABEL:
                continue
            idx_map[str(p_name)] = arr[inv == u_i].mean(axis=0)

    pds: List[float] = []
    maes: List[float] = []
    for p_name in sorted(set(str(x) for x in cohort_perts)):
        mu_r = r_map[p_name]
        mu_p = p_map[p_name]
        d_r = mu_r - ctrl_real
        d_p = mu_p - ctrl_pred
        if d_r.std() == 0.0 or d_p.std() == 0.0:
            corr = 0.0
        else:
            c_val = float(pearsonr(d_p, d_r)[0])
            corr = 0.0 if np.isnan(c_val) else c_val
        pds.append(corr)
        maes.append(float(np.mean(np.abs(mu_p - mu_r))))

    return {
        "pearson_delta_mean": float(np.mean(pds)),
        "mae_mean": float(np.mean(maes)),
        "n_perts": int(len(cohort_perts)),
    }
