#!/usr/bin/env python3
"""
VCC 2026 submission format builder, validator, `.vcc` packager, and SHA-256 timestamp tools.

Mirrors the official validation contract enforced by `vcc-cli==0.3.0` (`vcc/prep.py`) and
`ArcInstitute/cell-eval2==0.22.0` (`docs/vcc2026_metrics/vcc2026-metrics.md` §7):
  - Upload contract (`vcc prep` / `prediction.h5ad`):
      * One AnnData covering all required contexts (`obs["context"]`, e.g. D/E/F).
      * `obs["target_gene"]` contains gene symbols matching `pert_counts.csv` (no extras, none missing).
      * NO `non-targeting` control rows (`reject_controls=True`).
      * Exact `cells_per_pert` (default 400) cells per perturbation per context.
      * Exact gene axis order matching `gene_names.csv` (`18,533` genes for VCC 2026).
      * `.X` is sparse CSR containing finite, non-negative WHOLE INTEGER counts (no normalization/log1p),
        no explicit stored zeros, per-cell count sum <= 1,000,000, total stored nnz <= 4,750,000,000,
        and no all-zero perturbation.
  - Local `cell-eval2` scoring bridge (`attach_real_controls_for_celleval2`):
      * `cell_eval2.io.validate_pair` requires `pred` and `real` to share the exact same
        perturbation label set INCLUDING `non-targeting`. This helper splices real reference
        `non-targeting` cells onto an upload-format prediction for local `cell-eval2` scoring.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import tarfile
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

DEFAULT_PERT_COL = "target_gene"
DEFAULT_CONTEXT_COL = "context"
DEFAULT_NTC_NAME = "non-targeting"

EXPECTED_GENE_DIM = 18533
EXPECTED_CELLS_PER_PERT = 400
MAX_CELL_DIM = 400_000
MAX_NNZ = 4_750_000_000
MAX_COUNTS_PER_CELL = 1_000_000

PRED_MEMBER = "pred.h5ad.zst"
META_MEMBER = "meta.json"
META_SCHEMA = 1

_GENE_LIST_HEADERS = frozenset({"gene_name", "gene", "gene_symbol", "gene_names", "genes"})
_GUIDE_SUFFIX_RE = re.compile(r"-\d$")


class FormatValidationError(ValueError):
    """Raised when a prediction violates the VCC 2026 submission format contract."""


@dataclass
class ValidationReport:
    passed: bool
    n_cells: int
    n_genes: int
    nnz: int
    contexts: Dict[str, int]
    n_perts_per_context: Dict[str, int]
    all_finite: bool
    all_nonnegative: bool
    all_integer: bool
    min_cell_total: float
    max_cell_total: float
    empty_perturbations: List[str]
    errors: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "passed": self.passed,
            "n_cells": self.n_cells,
            "n_genes": self.n_genes,
            "nnz": self.nnz,
            "contexts": self.contexts,
            "n_perts_per_context": self.n_perts_per_context,
            "all_finite": self.all_finite,
            "all_nonnegative": self.all_nonnegative,
            "all_integer": self.all_integer,
            "min_cell_total": self.min_cell_total,
            "max_cell_total": self.max_cell_total,
            "empty_perturbations": self.empty_perturbations,
            "errors": self.errors,
            "notes": self.notes,
        }


def sha256_file(path: str | Path) -> str:
    """Compute hex SHA-256 of a file using 1 MiB streaming chunks."""
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_gene_names(path: str | Path) -> List[str]:
    """
    Read ordered gene symbol list from `gene_names.csv` (tolerating header row like `gene_name`).
    Rejects empty or duplicate symbols.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Gene list not found: {p}")
    series = pd.read_csv(p, header=None).iloc[:, 0]
    if series.empty:
        raise ValueError(f"Gene list is empty: {p}")
    genes = [str(x).strip() for x in series.tolist() if pd.notna(x) and str(x).strip()]
    if genes and genes[0].lower() in _GENE_LIST_HEADERS:
        genes = genes[1:]
    if not genes:
        raise ValueError(f"Gene list contains no gene symbols: {p}")
    seen: set[str] = set()
    dups = sorted({g for g in genes if g in seen or seen.add(g)})
    if dups:
        raise ValueError(f"Duplicate gene symbols in {p}: {dups[:10]}")
    return genes


def load_pert_counts(
    path: str | Path,
    *,
    pert_col: str = DEFAULT_PERT_COL,
    context_col: str = DEFAULT_CONTEXT_COL,
    ntc_name: str = DEFAULT_NTC_NAME,
    default_n_cells: int = EXPECTED_CELLS_PER_PERT,
) -> Dict[Optional[str], Dict[str, int]]:
    """
    Load official perturbation list from `pert_counts.csv` or a plain-text `.txt` list.
    Returns `{context_or_None: {target_gene: n_cells}}`, excluding `non-targeting`.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Perturbation list not found: {p}")

    if p.suffix.lower() == ".txt":
        lines = [
            ln.strip()
            for ln in p.read_text(encoding="utf-8").splitlines()
            if ln.strip() and not ln.strip().startswith("#")
        ]
        if lines and lines[0].lower() in {pert_col.lower(), "target", "perturbation", "gene"}:
            lines = lines[1:]
        rows = {ln: int(default_n_cells) for ln in lines if ln != ntc_name}
        if not rows:
            raise ValueError(f"Perturbation list is empty: {p}")
        return {None: rows}

    df = pd.read_csv(p)
    if pert_col not in df.columns:
        if len(df.columns) == 1:
            # Handle headerless single-column CSV gracefully while warning if first entry looks like target
            df = pd.read_csv(p, header=None, names=[pert_col])
        elif "target" in df.columns:
            pert_col = "target"
        else:
            raise ValueError(f"Column {pert_col!r} not found in {p}; columns={list(df.columns)}")

    has_counts = "n_cells" in df.columns

    def _extract(group: pd.DataFrame) -> Dict[str, int]:
        out: Dict[str, int] = {}
        for _, row in group.iterrows():
            val = row[pert_col]
            if pd.isna(val):
                continue
            name = str(val).strip()
            if not name or name == ntc_name:
                continue
            n = int(default_n_cells)
            if has_counts and pd.notna(row["n_cells"]):
                cand = int(float(row["n_cells"]))
                if cand > 0:
                    n = cand
            out[name] = n
        return out

    expected: Dict[Optional[str], Dict[str, int]] = {}
    if context_col in df.columns:
        for ctx, grp in df.groupby(df[context_col].astype(str), sort=False):
            expected[str(ctx)] = _extract(grp)
    else:
        expected[None] = _extract(df)

    if not any(expected.values()):
        raise ValueError(f"No perturbations found in {p}")
    return expected


def load_manifest(path: str | Path) -> Optional[Dict[str, Any]]:
    """Load `manifest.json` from a VCC controls bundle if present and valid."""
    p = Path(path)
    if not p.is_file():
        return None
    data = json.loads(p.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        return None
    return data


def _has_fractional_values(data: np.ndarray, epsilon: float = 1e-6) -> bool:
    block = 1 << 22
    for start in range(0, data.size, block):
        chunk = data[start : start + block]
        if np.any(np.abs(chunk - np.trunc(chunk)) > epsilon):
            return True
    return False


def _validate_sparse_structure(X: Any) -> Optional[str]:
    if not sp.issparse(X):
        return None
    fmt = getattr(X, "format", None)
    if fmt not in ("csr", "csc"):
        return None
    indptr = np.asarray(X.indptr).astype(np.int64, copy=False)
    indices = np.asarray(X.indices)
    data_len = int(X.data.shape[0])
    n_major, n_minor = (X.shape[0], X.shape[1]) if fmt == "csr" else (X.shape[1], X.shape[0])
    if indptr.shape[0] != n_major + 1:
        return f"indptr length {indptr.shape[0]} != {n_major + 1}"
    if int(indptr[0]) != 0:
        return f"indptr[0] = {int(indptr[0])} != 0"
    if np.any(np.diff(indptr) < 0):
        return "indptr is not monotonic"
    if indices.shape[0] != data_len:
        return f"indices len {indices.shape[0]} != data len {data_len}"
    nnz = int(indptr[-1])
    if nnz > data_len:
        return f"indptr[-1] ({nnz}) > data len ({data_len})"
    live = indices[:nnz]
    if live.size and (int(live.min()) < 0 or int(live.max()) >= n_minor):
        return f"column/row index out of bounds [0, {n_minor})"
    return None


def validate_vcc2026_anndata(
    adata: ad.AnnData,
    *,
    gene_names: Optional[Sequence[str]] = None,
    pert_counts: Optional[Dict[Optional[str], Dict[str, int]]] = None,
    required_contexts: Optional[Sequence[str]] = None,
    cells_per_pert: Optional[int] = EXPECTED_CELLS_PER_PERT,
    pert_col: str = DEFAULT_PERT_COL,
    context_col: str = DEFAULT_CONTEXT_COL,
    ntc_name: str = DEFAULT_NTC_NAME,
    expected_gene_dim: Optional[int] = EXPECTED_GENE_DIM,
    max_cell_dim: Optional[int] = MAX_CELL_DIM,
    max_nnz: Optional[int] = MAX_NNZ,
    max_counts_per_cell: Optional[int] = MAX_COUNTS_PER_CELL,
    reject_controls: bool = True,
    require_sparse: bool = True,
) -> ValidationReport:
    """
    Validate an AnnData prediction against every VCC 2026 submission requirement.
    """
    errors: List[str] = []
    notes: List[str] = []

    if adata.X is None:
        return ValidationReport(
            passed=False, n_cells=0, n_genes=0, nnz=0, contexts={}, n_perts_per_context={},
            all_finite=False, all_nonnegative=False, all_integer=False,
            min_cell_total=0.0, max_cell_total=0.0, empty_perturbations=[],
            errors=["adata.X is None"],
        )

    n_cells, n_genes = int(adata.shape[0]), int(adata.shape[1])
    if n_cells == 0:
        errors.append("AnnData contains 0 cells.")

    if adata.var_names.astype(str).duplicated().any():
        errors.append("Duplicate gene names in adata.var_names.")
    if gene_names is not None:
        want_genes = [str(g) for g in gene_names]
        got_genes = adata.var_names.astype(str).tolist()
        if got_genes != want_genes:
            missing = sorted(set(want_genes) - set(got_genes))
            extra = sorted(set(got_genes) - set(want_genes))
            if not missing and not extra:
                errors.append("Gene names match set but are out of order relative to gene_names.")
            else:
                errors.append(
                    f"Gene axis mismatch: missing={len(missing)} {missing[:5]}, "
                    f"extra={len(extra)} {extra[:5]}"
                )
    if expected_gene_dim is not None and n_genes != expected_gene_dim:
        errors.append(f"Gene dimension {n_genes} != expected {expected_gene_dim}.")
    if max_cell_dim is not None and n_cells > max_cell_dim:
        errors.append(f"Cell count {n_cells} exceeds max_cell_dim={max_cell_dim}.")

    if pert_col not in adata.obs.columns:
        errors.append(f"Missing perturbation column {pert_col!r} in adata.obs.")
        pert_vals = pd.Series([], dtype=str)
    else:
        if adata.obs[pert_col].isna().any():
            errors.append(f"Null values found in adata.obs[{pert_col!r}].")
        pert_vals = adata.obs[pert_col].astype(str)
        if reject_controls:
            n_ctrl = int((pert_vals == ntc_name).sum())
            if n_ctrl > 0:
                errors.append(
                    f"Submission contains {n_ctrl} '{ntc_name}' control cells "
                    "(controls are forbidden on upload)."
                )

    ctx_counts: Dict[str, int] = {}
    perts_per_ctx: Dict[str, int] = {}
    if required_contexts is not None and len(required_contexts) > 0:
        req_tuple = tuple(str(c) for c in required_contexts)
        if context_col not in adata.obs.columns:
            errors.append(f"Missing context column {context_col!r} in adata.obs.")
        else:
            ctx_vals = adata.obs[context_col].astype(str)
            present = set(ctx_vals.unique())
            unknown = sorted(present - set(req_tuple))
            missing = [c for c in req_tuple if c not in present]
            if unknown:
                errors.append(f"Unknown context labels in {context_col!r}: {unknown}")
            if missing:
                errors.append(f"Missing required context labels in {context_col!r}: {missing}")
            for c in req_tuple:
                mask = ctx_vals == c
                ctx_counts[c] = int(mask.sum())
                if pert_col in adata.obs.columns:
                    perts_per_ctx[c] = int(pert_vals[mask].nunique())

            if pert_counts is not None and pert_col in adata.obs.columns and not missing:
                for c in req_tuple:
                    want_map = pert_counts.get(c, pert_counts.get(None))
                    if want_map is None:
                        errors.append(f"No target list available for context {c!r}.")
                        continue
                    mask = ctx_vals == c
                    obs_labels = pert_vals[mask]
                    obs_set = set(obs_labels.unique()) - {ntc_name}
                    want_set = set(want_map.keys())
                    miss_p = sorted(want_set - obs_set)
                    extra_p = sorted(obs_set - want_set)
                    if miss_p or extra_p:
                        msg = f"Context {c}: target mismatch (missing={len(miss_p)} {miss_p[:5]}, unexpected={len(extra_p)} {extra_p[:5]})."
                        if any(_GUIDE_SUFFIX_RE.search(x) for x in extra_p):
                            msg += " Detected guide suffixes like '-1'; use gene symbols."
                        errors.append(msg)
                    vc = obs_labels.value_counts()
                    wrong_counts = []
                    for g, want_n in sorted(want_map.items()):
                        exp_n = want_n if want_n is not None else cells_per_pert
                        if exp_n is None or exp_n < 1:
                            continue
                        got_n = int(vc.get(g, 0))
                        if got_n != exp_n:
                            wrong_counts.append((g, got_n, exp_n))
                    if wrong_counts:
                        sample = ", ".join(f"{g}:{got}!={exp}" for g, got, exp in wrong_counts[:5])
                        errors.append(f"Context {c}: {len(wrong_counts)} targets have wrong cell counts ({sample}).")
    elif pert_counts is not None and pert_col in adata.obs.columns:
        want_map = pert_counts.get(None) or {k: v for m in pert_counts.values() for k, v in m.items()}
        obs_set = set(pert_vals.unique()) - {ntc_name}
        want_set = set(want_map.keys())
        miss_p = sorted(want_set - obs_set)
        extra_p = sorted(obs_set - want_set)
        if miss_p or extra_p:
            errors.append(f"Target mismatch: missing={len(miss_p)} {miss_p[:5]}, unexpected={len(extra_p)} {extra_p[:5]}.")
        vc = pert_vals.value_counts()
        wrong_counts = []
        for g, want_n in sorted(want_map.items()):
            exp_n = want_n if want_n is not None else cells_per_pert
            if exp_n is not None and exp_n >= 1 and int(vc.get(g, 0)) != exp_n:
                wrong_counts.append((g, int(vc.get(g, 0)), exp_n))
        if wrong_counts:
            errors.append(f"{len(wrong_counts)} targets have wrong cell counts ({wrong_counts[:5]}).")

    X = adata.X
    if require_sparse and not sp.issparse(X):
        errors.append("adata.X is dense; must be stored sparse (scipy.sparse.csr_matrix).")

    struct_err = _validate_sparse_structure(X)
    if struct_err is not None:
        errors.append(f"Sparse matrix structure error: {struct_err}")

    if sp.issparse(X):
        nnz = int(X.nnz)
        data = X.data
        if data.size and np.any(data == 0):
            notes.append(f"Sparse X contains {int((data == 0).sum())} explicit zero entries (counts toward nnz).")
    else:
        arr = np.asarray(X)
        nnz = int(np.count_nonzero(arr))
        data = arr.ravel()

    if max_nnz is not None and nnz > max_nnz:
        errors.append(f"Stored nnz={nnz:,} exceeds max_nnz={max_nnz:,}.")

    all_finite = True
    all_nonnegative = True
    all_integer = True
    if data.size > 0:
        if not np.all(np.isfinite(data)):
            all_finite = False
            errors.append("adata.X contains NaN or Inf values.")
        if float(np.min(data)) < 0.0:
            all_nonnegative = False
            errors.append(f"adata.X contains negative values (min={float(np.min(data)):.6g}).")
        if not np.issubdtype(data.dtype, np.integer) and _has_fractional_values(data):
            all_integer = False
            errors.append("adata.X contains fractional values; raw integer counts required.")

    totals = np.asarray(X.sum(axis=1)).ravel() if sp.issparse(X) else np.asarray(X).sum(axis=1)
    min_cell_total = float(totals.min()) if totals.size else 0.0
    max_cell_total = float(totals.max()) if totals.size else 0.0
    if max_counts_per_cell is not None and totals.size and max_cell_total > max_counts_per_cell:
        over = int((totals > max_counts_per_cell).sum())
        errors.append(
            f"{over} cell(s) exceed max_counts_per_cell={max_counts_per_cell:,} "
            f"(max={max_cell_total:,.1f})."
        )

    empty_perts: List[str] = []
    if pert_col in adata.obs.columns and totals.size:
        p_arr = pert_vals.to_numpy()
        for p_name in sorted(set(p_arr.tolist())):
            if not np.any(totals[p_arr == p_name] > 0):
                empty_perts.append(p_name)
        if empty_perts:
            errors.append(f"{len(empty_perts)} all-zero perturbations: {empty_perts[:5]}")

    return ValidationReport(
        passed=len(errors) == 0,
        n_cells=n_cells,
        n_genes=n_genes,
        nnz=nnz,
        contexts=ctx_counts,
        n_perts_per_context=perts_per_ctx,
        all_finite=all_finite,
        all_nonnegative=all_nonnegative,
        all_integer=all_integer,
        min_cell_total=min_cell_total,
        max_cell_total=max_cell_total,
        empty_perturbations=empty_perts,
        errors=errors,
        notes=notes,
    )


def emit_dispersed_counts_for_context(
    ctrl_counts: Any,
    pert_names: Sequence[str],
    gene_ratios: np.ndarray,
    *,
    cells_per_pert: int = EXPECTED_CELLS_PER_PERT,
    max_counts_per_cell: int = MAX_COUNTS_PER_CELL,
    seed: int = 42,
) -> sp.csr_matrix:
    """
    Synthesize realistic sparse integer raw-count cells for a list of perturbations in one context
    by resampling empirical control cells with replacement and multiplicatively scaling by
    predicted gene-wise fold ratios `gene_ratios` (shape P x G), followed by unbiased stochastic
    rounding and strict per-cell library-size capping.

    Preserves single-cell dispersion and dropout sparsity while enforcing exact non-negative
    integers stored in canonical float32 CSR without explicit zeros.
    """
    ctrl_csr = ctrl_counts.tocsr() if sp.issparse(ctrl_counts) else sp.csr_matrix(ctrl_counts)
    ctrl_csr = ctrl_csr.astype(np.float64, copy=False)
    n_ctrl, n_genes = ctrl_csr.shape
    if n_ctrl == 0:
        raise ValueError("Control count matrix has 0 cells.")

    ratios = np.asarray(gene_ratios, dtype=np.float64)
    if ratios.shape != (len(pert_names), n_genes):
        raise ValueError(f"gene_ratios shape {ratios.shape} != ({len(pert_names)}, {n_genes})")
    ratios = np.where(np.isfinite(ratios) & (ratios >= 0.0), ratios, 1.0)

    # Identify non-empty control cells so a drawn cell can never leave an all-zero perturbation
    ctrl_totals = np.asarray(ctrl_csr.sum(axis=1)).ravel()
    nonempty_ctrl = np.flatnonzero(ctrl_totals > 0)
    sample_pool = nonempty_ctrl if nonempty_ctrl.size > 0 else np.arange(n_ctrl)
    top_ctrl_gene = int(np.argmax(np.asarray(ctrl_csr.sum(axis=0)).ravel()))

    rng = np.random.default_rng(seed)
    blocks: List[sp.csr_matrix] = []
    cap_soft = float(max_counts_per_cell) * 0.999 if max_counts_per_cell > 0 else None

    for p_idx in range(len(pert_names)):
        draw = sample_pool[rng.integers(0, sample_pool.size, size=cells_per_pert)]
        blk = ctrl_csr[draw].copy()
        r_vec = ratios[p_idx]
        if blk.data.size > 0:
            blk.data = blk.data * r_vec[blk.indices]

            # Pre-round library cap so stochastic rounding cannot exceed max_counts_per_cell
            if cap_soft is not None:
                row_sums = np.asarray(blk.sum(axis=1)).ravel()
                over = row_sums > cap_soft
                if np.any(over):
                    scale = np.ones_like(row_sums)
                    scale[over] = cap_soft / row_sums[over]
                    row_nnz = np.diff(blk.indptr)
                    blk.data *= np.repeat(scale, row_nnz)

            # Unbiased stochastic rounding to non-negative integers
            floored = np.floor(blk.data)
            frac = blk.data - floored
            blk.data = (floored + (rng.random(blk.data.size) < frac)).astype(np.float32)

            # Hard guard against rare rounding-up overflow past max_counts_per_cell
            if max_counts_per_cell > 0:
                row_sums_int = np.asarray(blk.sum(axis=1)).ravel()
                for row_i in np.flatnonzero(row_sums_int > max_counts_per_cell):
                    s_i, e_i = int(blk.indptr[row_i]), int(blk.indptr[row_i + 1])
                    excess = float(row_sums_int[row_i] - max_counts_per_cell)
                    if e_i > s_i:
                        j_max = s_i + int(np.argmax(blk.data[s_i:e_i]))
                        blk.data[j_max] = max(0.0, float(blk.data[j_max] - excess))

        blk.eliminate_zeros()
        blk.sum_duplicates()
        blk.sort_indices()

        # Guarantee no all-zero perturbation
        if blk.nnz == 0:
            fallback = ctrl_csr[draw].copy()
            fallback.data = np.rint(fallback.data).astype(np.float32)
            fallback.eliminate_zeros()
            if fallback.nnz > 0:
                blk = fallback
            else:
                blk = sp.csr_matrix(
                    ([1.0], ([0], [top_ctrl_gene])),
                    shape=(cells_per_pert, n_genes),
                    dtype=np.float32,
                )
        blocks.append(blk)

    stacked = sp.vstack(blocks, format="csr")
    stacked.data = stacked.data.astype(np.float32, copy=False)
    stacked.eliminate_zeros()
    stacked.sort_indices()
    return stacked


def package_vcc_archive(
    adata: ad.AnnData,
    output_vcc: str | Path,
    *,
    contexts: Sequence[str],
    perts_path: Optional[str | Path] = None,
    panel_id: Optional[str] = None,
    cli_version: str = "prospective-0.1.0",
) -> int:
    """
    Package an already-validated AnnData into a `.vcc` tar archive (`meta.json` + `pred.h5ad.zst`)
    matching `vcc prep`'s deterministic normalized-header tar layout. Requires `zstandard`.
    """
    import zstandard as zstd

    perts_fp = None
    if perts_path is not None and Path(perts_path).is_file():
        perts_dict = load_pert_counts(perts_path)
        names = sorted({k for m in perts_dict.values() for k in m.keys()})
        if names:
            perts_fp = "sha256:" + hashlib.sha256("\n".join(names).encode("utf-8")).hexdigest()

    out_path = Path(output_vcc)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    nnz = int(adata.X.nnz) if sp.issparse(adata.X) else int(np.count_nonzero(adata.X))

    with TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        h5ad_tmp = tmp_dir / "pred.h5ad"
        zst_tmp = tmp_dir / PRED_MEMBER
        meta_tmp = tmp_dir / META_MEMBER

        adata.write_h5ad(h5ad_tmp)
        meta = {
            "schema": META_SCHEMA,
            "nnz": nnz,
            "n_obs": int(adata.n_obs),
            "n_vars": int(adata.n_vars),
            "cli_version": cli_version,
            "contexts": list(contexts),
        }
        if perts_fp is not None:
            meta["perts_fingerprint"] = perts_fp
        if panel_id is not None:
            meta["panel_id"] = panel_id
        meta_tmp.write_text(json.dumps(meta, sort_keys=True), encoding="utf-8")

        cctx = zstd.ZstdCompressor(level=3, threads=min(8, os.cpu_count() or 1))
        with h5ad_tmp.open("rb") as src, zst_tmp.open("wb") as dst:
            cctx.copy_stream(src, dst)
        h5ad_tmp.unlink()

        def _norm_tar(ti: tarfile.TarInfo) -> tarfile.TarInfo:
            ti.uid = ti.gid = 0
            ti.uname = ti.gname = ""
            ti.mtime = 0
            return ti

        with tarfile.open(out_path, "w") as tar:
            tar.add(meta_tmp, arcname=META_MEMBER, filter=_norm_tar)
            tar.add(zst_tmp, arcname=PRED_MEMBER, filter=_norm_tar)

    return nnz


def attach_real_controls_for_celleval2(
    pred_upload_h5ad: str | Path,
    real_h5ad: str | Path,
    out_scoring_h5ad: str | Path,
    *,
    pert_col: str = DEFAULT_PERT_COL,
    control_label: str = DEFAULT_NTC_NAME,
    context_col: Optional[str] = None,
    context_value: Optional[str] = None,
) -> Path:
    """
    Bridge VCC 2026 upload contract (no `non-targeting` rows) to `cell-eval2` scoring contract
    (`io.validate_pair` requires `non-targeting` present in both `pred` and `real` with identical
    gene axes). Extracts real control cells from `real_h5ad` and appends them onto `pred`.
    """
    pred = ad.read_h5ad(pred_upload_h5ad)
    if context_col and context_value is not None and context_col in pred.obs.columns:
        pred = pred[pred.obs[context_col].astype(str) == str(context_value)].copy()

    real = ad.read_h5ad(real_h5ad)
    if pert_col not in real.obs.columns:
        raise ValueError(f"Real h5ad missing pert_col={pert_col!r}")
    real_ctrl = real[real.obs[pert_col].astype(str) == control_label].copy()
    if real_ctrl.n_obs == 0:
        raise ValueError(f"No {control_label!r} cells found in {real_h5ad}")

    # Align gene order to real.var_names
    real_genes = real.var_names.astype(str).tolist()
    pred = pred[:, real_genes].copy()
    real_ctrl = real_ctrl[:, real_genes].copy()

    obs_pred = pd.DataFrame({pert_col: pred.obs[pert_col].astype(str).to_numpy()},
                            index=[f"pred_{i}" for i in range(pred.n_obs)])
    obs_ctrl = pd.DataFrame({pert_col: real_ctrl.obs[pert_col].astype(str).to_numpy()},
                            index=[f"ctrl_{i}" for i in range(real_ctrl.n_obs)])

    X_pred = pred.X.tocsr() if sp.issparse(pred.X) else sp.csr_matrix(pred.X)
    X_ctrl = real_ctrl.X.tocsr() if sp.issparse(real_ctrl.X) else sp.csr_matrix(real_ctrl.X)
    X_comb = sp.vstack([X_pred.astype(np.float32), X_ctrl.astype(np.float32)], format="csr")
    obs_comb = pd.concat([obs_pred, obs_ctrl], axis=0)
    var_comb = pd.DataFrame(index=real_genes)

    out_ad = ad.AnnData(X=X_comb, obs=obs_comb, var=var_comb)
    out_p = Path(out_scoring_h5ad)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    out_ad.write_h5ad(out_p)
    return out_p


def write_timestamp_bundle(
    outdir: str | Path,
    *,
    output_files: Sequence[Path],
    input_files: Dict[str, str],
    validation_report: ValidationReport,
    method_meta: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Compute SHA-256 checksums for all generated prediction artifacts and write:
      - `SHA256SUMS.txt`
      - `prediction_manifest.json`
      - `TIMESTAMP_DECLARATION.md`
    """
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    utc_now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    output_hashes: Dict[str, Dict[str, Any]] = {}
    sum_lines: List[str] = []
    for p in output_files:
        if p.exists() and p.is_file():
            digest = sha256_file(p)
            size_b = int(p.stat().st_size)
            output_hashes[p.name] = {"sha256": digest, "size_bytes": size_b}
            sum_lines.append(f"{digest}  {p.name}")

    input_hashes: Dict[str, Dict[str, Any]] = {}
    for label, path_str in sorted(input_files.items()):
        p = Path(path_str)
        if p.exists() and p.is_file():
            input_hashes[label] = {
                "filename": p.name,
                "sha256": sha256_file(p),
                "size_bytes": int(p.stat().st_size),
            }

    manifest = {
        "schema": "vcc_prospective_prediction_v1",
        "created_utc": utc_now,
        "method": method_meta,
        "inputs": input_hashes,
        "outputs": output_hashes,
        "validation": validation_report.to_dict(),
    }

    (out / "SHA256SUMS.txt").write_text("\n".join(sum_lines) + "\n", encoding="utf-8")
    (out / "prediction_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    decl_lines = [
        "# Prospective Prediction Cryptographic Timestamp Declaration",
        "",
        f"- **Generated (UTC)**: `{utc_now}`",
        f"- **Model / Predictor**: `{method_meta.get('method', 'unknown')}`",
        f"- **Source Lines**: `{method_meta.get('source_lines', [])}`",
        f"- **Format Validation Passed**: `{validation_report.passed}`",
        f"- **Shape**: `{validation_report.n_cells:,}` cells × `{validation_report.n_genes:,}` genes "
        f"(stored `nnz = {validation_report.nnz:,}`)",
        f"- **Contexts**: `{validation_report.contexts}`",
        f"- **Counts Integrity**: `all_finite={validation_report.all_finite}`, "
        f"`all_nonnegative={validation_report.all_nonnegative}`, "
        f"`all_integer={validation_report.all_integer}`, "
        f"`max_cell_total={validation_report.max_cell_total:,.0f}`",
        "",
        "## Output SHA-256 Checksums",
        "",
        "```text",
        *sum_lines,
        "```",
        "",
        "## Input SHA-256 Checksums (Not Redistributed)",
        "",
        "```json",
        json.dumps(input_hashes, indent=2, sort_keys=True),
        "```",
        "",
        "## Pre-Release Verification & Timestamping Instructions",
        "",
        "1. Verify locally: `shasum -a 256 -c SHA256SUMS.txt`.",
        "2. Before ground-truth public release, commit and push `SHA256SUMS.txt`, "
        "`prediction_manifest.json`, `TIMESTAMP_DECLARATION.md`, and the prediction code/weights "
        "to a public GitHub repository (and optionally anchor via OpenTimestamps `ots stamp SHA256SUMS.txt`).",
        "3. Do **not** commit or redistribute restricted competition input files "
        "(`context_*.h5ad`, `pert_counts.csv`, `gene_names.csv`) unless authorized or released under CC0.",
        "",
    ]
    (out / "TIMESTAMP_DECLARATION.md").write_text("\n".join(decl_lines), encoding="utf-8")
    return manifest
