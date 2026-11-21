#!/usr/bin/env python3
"""
End-to-end prospective VCC 2026 prediction pipeline:
  unperturbed control .h5ad + target perturbation list
    -> zero-shot cross-line prediction (using A-T2 predictors.py or standalone source .h5ad)
    -> dispersed sparse raw-integer count synthesis (400 cells / perturbation / context)
    -> strict VCC 2026 format validation (matching `vcc prep`)
    -> optional `.vcc` archive packaging
    -> cryptographic SHA-256 checksums + timestamp declaration

Usage examples:
  # 1. Official VCC 2026 final bundle (when provided by human or public VCA release):
  python3 prospective_predict.py \
      --bundle-dir ~/vcc-final \
      --source-lines-dir ../../shared/lines \
      --method Hybrid_ProjShared \
      --outdir ../../results/prospective_vcc2026

  # 2. Simulation on Replogle/Nadig or VCC 2025 control data:
  python3 prospective_predict.py \
      --controls SIM=ctrl_only.h5ad \
      --perts targets.txt \
      --source-lines-dir ../../shared/lines \
      --method Hybrid_ProjShared \
      --simulation \
      --outdir /tmp/sim_prospective_out
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

PROSPECTIVE_DIR = Path(__file__).resolve().parent
SRC_DIR = PROSPECTIVE_DIR.parent
P2_ROOT = PROSPECTIVE_DIR.parents[2]
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))
if str(PROSPECTIVE_DIR) not in sys.path:
    sys.path.insert(0, str(PROSPECTIVE_DIR))

from data_loader import check_prereg_gate, list_available_lines, load_line_dataset
from predictors import (
    SourceLineInput,
    SourceModelState,
    cipher_unit_response_columns,
    compute_control_covariance,
    predict_zero_shot_suite,
)
from vcc2026_format import (
    DEFAULT_CONTEXT_COL,
    DEFAULT_NTC_NAME,
    DEFAULT_PERT_COL,
    EXPECTED_CELLS_PER_PERT,
    EXPECTED_GENE_DIM,
    MAX_CELL_DIM,
    MAX_COUNTS_PER_CELL,
    MAX_NNZ,
    emit_dispersed_counts_for_context,
    load_gene_names,
    load_manifest,
    load_pert_counts,
    package_vcc_archive,
    validate_vcc2026_anndata,
    write_timestamp_bundle,
)


def _normalize_counts_log1p(
    X_raw: Any,
    target_sum: Optional[float] = None,
) -> Tuple[np.ndarray, float]:
    """
    Addendum-B library-size normalization + log1p on raw count matrix (N_cells, G).
    If target_sum is None, uses median non-zero library size across the control cells.
    Returns `(X_log1p_dense_float64, resolved_target_sum)`.
    """
    if sp.issparse(X_raw):
        X = X_raw.toarray().astype(np.float64, copy=False)
    else:
        X = np.asarray(X_raw, dtype=np.float64)

    totals = X.sum(axis=1)
    nz = totals[totals > 0]
    if nz.size == 0:
        raise ValueError("All control cells have zero total counts.")
    ts = float(np.median(nz)) if target_sum is None else float(target_sum)
    safe = np.where(totals > 0, totals, 1.0)[:, None]
    X_norm = X * (ts / safe)
    X_norm[totals <= 0] = 0.0
    return np.log1p(np.maximum(X_norm, 0.0)), ts


def _extract_control_matrix(
    h5ad_path: Path,
    panel_genes: Sequence[str],
    *,
    pert_col: str = DEFAULT_PERT_COL,
    ntc_name: str = DEFAULT_NTC_NAME,
    max_cells: Optional[int] = None,
    seed: int = 42,
) -> sp.csr_matrix:
    """
    Load raw counts for unperturbed control cells from `h5ad_path` aligned to `panel_genes`.
    Accepts either a control-only h5ad or a mixed h5ad containing `pert_col == ntc_name`.
    """
    adata = ad.read_h5ad(h5ad_path)
    if pert_col in adata.obs.columns:
        vals = adata.obs[pert_col].astype(str)
        if (vals == ntc_name).any():
            adata = adata[vals == ntc_name].copy()
    if adata.n_obs == 0:
        raise ValueError(f"No control cells found in {h5ad_path}")

    if max_cells is not None and adata.n_obs > max_cells:
        rng = np.random.default_rng(seed)
        keep_idx = np.sort(rng.choice(adata.n_obs, size=max_cells, replace=False))
        adata = adata[keep_idx].copy()

    adata.var_names = adata.var_names.astype(str)
    if adata.var_names.duplicated().any():
        raise ValueError(f"Duplicate gene names in control h5ad: {h5ad_path}")

    src_genes = adata.var_names.tolist()
    if src_genes != list(panel_genes):
        missing = sorted(set(panel_genes) - set(src_genes))
        if missing:
            raise ValueError(
                f"{h5ad_path} is missing {len(missing)} panel genes (first 5: {missing[:5]})"
            )
        adata = adata[:, list(panel_genes)].copy()

    X = adata.X
    if sp.issparse(X):
        X_csr = X.tocsr().astype(np.float32, copy=False)
    else:
        X_csr = sp.csr_matrix(np.asarray(X, dtype=np.float32))
    X_csr.eliminate_zeros()
    return X_csr


def _build_source_from_h5ad(
    name: str,
    h5ad_path: Path,
    shared_genes: Sequence[str],
    *,
    pert_col: str = DEFAULT_PERT_COL,
    ntc_name: str = DEFAULT_NTC_NAME,
    min_cells_per_pert: int = 5,
    max_ctrl_cells: int = 20000,
    seed: int = 42,
) -> SourceLineInput:
    """Build a `SourceLineInput` directly from a raw-count Perturb-seq `.h5ad`."""
    adata = ad.read_h5ad(h5ad_path)
    adata.var_names = adata.var_names.astype(str)
    if pert_col not in adata.obs.columns:
        raise ValueError(f"Source h5ad {h5ad_path} missing obs[{pert_col!r}]")
    adata = adata[:, list(shared_genes)].copy()
    labels = adata.obs[pert_col].astype(str).to_numpy()
    ctrl_idx = np.flatnonzero(labels == ntc_name)
    if ctrl_idx.size < 2:
        raise ValueError(f"Source h5ad {h5ad_path} needs >= 2 control cells ({ntc_name!r})")
    if ctrl_idx.size > max_ctrl_cells:
        rng = np.random.default_rng(seed)
        ctrl_idx = np.sort(rng.choice(ctrl_idx, size=max_ctrl_cells, replace=False))

    X_ctrl_raw = adata.X[ctrl_idx]
    ctrl_log1p, ts = _normalize_counts_log1p(X_ctrl_raw, target_sum=None)
    ctrl_mean = ctrl_log1p.mean(axis=0)
    ctrl_cov = compute_control_covariance(ctrl_log1p)

    g2i = {g: i for i, g in enumerate(shared_genes)}
    deltas_list: List[np.ndarray] = []
    t_idx_list: List[int] = []
    t_names_list: List[str] = []

    for pert in sorted(set(labels.tolist()) - {ntc_name}):
        p_idx = np.flatnonzero(labels == pert)
        if p_idx.size < min_cells_per_pert:
            continue
        p_log1p, _ = _normalize_counts_log1p(adata.X[p_idx], target_sum=ts)
        deltas_list.append(p_log1p.mean(axis=0) - ctrl_mean)
        t_idx_list.append(g2i.get(pert, -1))
        t_names_list.append(pert)

    if not deltas_list:
        raise ValueError(f"No perturbations with >= {min_cells_per_pert} cells in {h5ad_path}")

    return SourceLineInput(
        name=name,
        ctrl_cov=ctrl_cov,
        deltas=np.asarray(deltas_list, dtype=np.float64),
        target_indices=np.asarray(t_idx_list, dtype=np.int64),
        target_names=tuple(t_names_list),
    )


def _load_source_lines_from_shared(
    shared_dir: Path,
    panel_genes: Sequence[str],
    representation: str = "log1p",
    source_line_names: Optional[Sequence[str]] = None,
) -> Tuple[List[SourceLineInput], List[str]]:
    """
    Load pre-processed lines from `p2/shared/lines/` respecting `PREREG_COMMIT.txt` gate,
    restricted to genes shared with `panel_genes`.
    """
    import data_loader as dl

    orig_lines_dir = dl.LINES_DIR
    orig_prereg = dl.PREREG_COMMIT_FILE
    try:
        dl.LINES_DIR = shared_dir
        avail = dl.list_available_lines()
        chosen = list(source_line_names) if source_line_names else avail
        if not chosen:
            raise FileNotFoundError(f"No line directories found in {shared_dir}")
        check_prereg_gate(chosen)

        # Intersect gene sets across chosen source lines and target panel
        common_set = set(panel_genes)
        for ln in chosen:
            g_file = shared_dir / ln / "genes.txt"
            ln_genes = {x.strip() for x in g_file.read_text(encoding="utf-8").splitlines() if x.strip()}
            common_set &= ln_genes
        shared_genes = [g for g in panel_genes if g in common_set]
        if len(shared_genes) < 50:
            raise ValueError(f"Too few shared genes ({len(shared_genes)}) between sources and target panel")

        sources: List[SourceLineInput] = []
        for ln in chosen:
            ds = load_line_dataset(ln, representation=representation, common_genes=shared_genes)
            cov = compute_control_covariance(ds.ctrl_cells)
            # Restrict to reproducible perturbations if mask exists
            if ds.repr_mask is not None and np.any(ds.repr_mask):
                keep = ds.repr_mask
                deltas = ds.deltas[keep]
                t_idx = ds.target_indices[keep]
                t_names = tuple(t for t, k in zip(ds.target_names, keep) if k)
            else:
                deltas = ds.deltas
                t_idx = ds.target_indices
                t_names = ds.target_names
            sources.append(
                SourceLineInput(
                    name=ln,
                    ctrl_cov=cov,
                    deltas=deltas,
                    target_indices=t_idx,
                    target_names=t_names,
                )
            )
        return sources, shared_genes
    finally:
        dl.LINES_DIR = orig_lines_dir
        dl.PREREG_COMMIT_FILE = orig_prereg


def _deltas_to_fold_ratios(
    deltas_shared: np.ndarray,
    ctrl_log1p_full: np.ndarray,
    shared_col_indices: np.ndarray,
    *,
    representation: str = "log1p",
    max_abs_delta: float = 4.0,
    max_fold: float = 50.0,
) -> np.ndarray:
    """
    Convert predicted deltas `(P, G_shared)` into non-negative normalized count multipliers
    `(P, G_panel)` relative to target control baseline expression.
    """
    n_perts = deltas_shared.shape[0]
    n_panel = ctrl_log1p_full.shape[1]
    ratios = np.ones((n_perts, n_panel), dtype=np.float64)

    ctrl_mean_log1p = ctrl_log1p_full[:, shared_col_indices].mean(axis=0)
    if representation == "zscore":
        ctrl_std = ctrl_log1p_full[:, shared_col_indices].std(axis=0)
        deltas_log1p = deltas_shared * ctrl_std[None, :]
    else:
        deltas_log1p = deltas_shared

    deltas_clipped = np.clip(np.nan_to_num(deltas_log1p, nan=0.0), -max_abs_delta, max_abs_delta)
    ctrl_mean_norm = np.expm1(ctrl_mean_log1p)
    pred_mean_norm = np.expm1(np.maximum(0.0, ctrl_mean_log1p[None, :] + deltas_clipped))

    sub_ratios = np.ones_like(pred_mean_norm)
    pos = ctrl_mean_norm > 1e-8
    sub_ratios[:, pos] = pred_mean_norm[:, pos] / ctrl_mean_norm[None, pos]
    sub_ratios = np.clip(sub_ratios, 0.0, max_fold)
    ratios[:, shared_col_indices] = sub_ratios
    return ratios


def _resolve_inputs(args: argparse.Namespace) -> Tuple[
    Dict[str, Path], List[str], Dict[Optional[str], Dict[str, int]], Optional[str], Dict[str, str]
]:
    """Parse `--bundle-dir` or explicit `--controls` / `--genes` / `--perts` / `--contexts`."""
    input_files: Dict[str, str] = {}
    panel_id: Optional[str] = None
    contexts_order: Optional[List[str]] = None

    if args.bundle_dir is not None:
        bdir = Path(args.bundle_dir).expanduser().resolve()
        manifest_path = bdir / "manifest.json"
        manifest = load_manifest(manifest_path)
        if manifest is not None:
            input_files["manifest_json"] = str(manifest_path)
            panel_id = manifest.get("panel_id")
            if isinstance(manifest.get("contexts"), list):
                contexts_order = [str(c) for c in manifest["contexts"]]

        genes_path = args.genes or (bdir / "gene_names.csv")
        perts_path = args.perts or (bdir / "pert_counts.csv")
        if args.contexts:
            contexts_order = [c.strip() for c in args.contexts.split(",") if c.strip()]
        if contexts_order is None:
            # Discover context_*.h5ad
            found = sorted(bdir.glob("context_*.h5ad"))
            contexts_order = [p.stem.replace("context_", "", 1) for p in found]
        if not contexts_order:
            raise FileNotFoundError(f"Could not determine contexts in {bdir}")
        controls_map = {c: bdir / f"context_{c}.h5ad" for c in contexts_order}
    else:
        if not args.controls:
            raise ValueError("Must supply either --bundle-dir or --controls CTX=path.h5ad")
        controls_map = {}
        contexts_order = []
        for item in args.controls:
            if "=" not in item:
                raise ValueError(f"--controls expects CTX=PATH, got {item!r}")
            ctx, p_str = item.split("=", 1)
            ctx = ctx.strip()
            controls_map[ctx] = Path(p_str).expanduser().resolve()
            contexts_order.append(ctx)
        genes_path = args.genes
        perts_path = args.perts

    for ctx, p in controls_map.items():
        if not p.is_file():
            raise FileNotFoundError(f"Control h5ad for context {ctx!r} not found: {p}")
        input_files[f"control_{ctx}"] = str(p)

    if genes_path is not None:
        genes_p = Path(genes_path).expanduser().resolve()
        panel_genes = load_gene_names(genes_p)
        input_files["gene_names"] = str(genes_p)
    else:
        # Infer gene axis from first control h5ad
        first_ad = ad.read_h5ad(next(iter(controls_map.values())), backed="r")
        panel_genes = first_ad.var_names.astype(str).tolist()
        if hasattr(first_ad, "file") and first_ad.file is not None:
            first_ad.file.close()

    if perts_path is None:
        raise ValueError("Must supply --perts (pert_counts.csv or targets.txt)")
    perts_p = Path(perts_path).expanduser().resolve()
    pert_counts = load_pert_counts(perts_p, default_n_cells=args.cells_per_pert)
    input_files["pert_counts"] = str(perts_p)

    return controls_map, panel_genes, pert_counts, panel_id, input_files


def run_prospective(args: argparse.Namespace) -> Dict[str, Any]:
    controls_map, panel_genes, pert_counts, panel_id, input_files = _resolve_inputs(args)
    contexts = list(controls_map.keys())
    panel_gene_to_idx = {g: i for i, g in enumerate(panel_genes)}

    # 1. Load source training lines (or explicit self-control fallback)
    sources: List[SourceLineInput] = []
    source_mode = "unknown"
    if args.source_h5ad:
        # Intersect var_names across all source h5ads and panel_genes
        common_set = set(panel_genes)
        src_paths: List[Tuple[str, Path]] = []
        for item in args.source_h5ad:
            if "=" not in item:
                raise ValueError(f"--source-h5ad expects NAME=PATH, got {item!r}")
            nm, ps = item.split("=", 1)
            p = Path(ps).expanduser().resolve()
            src_ad = ad.read_h5ad(p, backed="r")
            common_set &= set(src_ad.var_names.astype(str).tolist())
            if hasattr(src_ad, "file") and src_ad.file is not None:
                src_ad.file.close()
            src_paths.append((nm.strip(), p))
            input_files[f"source_{nm.strip()}"] = str(p)
        shared_genes = [g for g in panel_genes if g in common_set]
        sources = [
            _build_source_from_h5ad(
                nm, p, shared_genes,
                pert_col=args.source_pert_col,
                ntc_name=args.ntc_name,
                min_cells_per_pert=args.min_source_cells,
                seed=args.seed,
            )
            for nm, p in src_paths
        ]
        source_mode = "source_h5ad"
    elif args.source_lines_dir is not None and Path(args.source_lines_dir).exists():
        wanted = [s.strip() for s in args.source_lines.split(",") if s.strip()] if args.source_lines else None
        sources, shared_genes = _load_source_lines_from_shared(
            Path(args.source_lines_dir).expanduser().resolve(),
            panel_genes,
            representation=args.representation,
            source_line_names=wanted,
        )
        source_mode = "shared_lines"
    elif args.self_control_only:
        shared_genes = list(panel_genes)
        source_mode = "self_control_fallback"
    else:
        raise FileNotFoundError(
            "No source training lines available. Provide --source-lines-dir (e.g. p2/shared/lines), "
            "--source-h5ad NAME=PATH, or pass --self-control-only for control-covariance fallback."
        )

    shared_col_idx = np.array([panel_gene_to_idx[g] for g in shared_genes], dtype=np.int64)
    shared_g2i = {g: i for i, g in enumerate(shared_genes)}

    model_state: Optional[SourceModelState] = None
    if sources:
        model_state = SourceModelState.fit(sources, mask_self=True)

    # 2. Predict and synthesize raw count cells per context
    X_blocks: List[sp.csr_matrix] = []
    obs_perts: List[str] = []
    obs_contexts: List[str] = []
    per_context_summary: Dict[str, Any] = {}

    for ctx_idx, ctx in enumerate(contexts):
        ctrl_csr = _extract_control_matrix(
            controls_map[ctx],
            panel_genes,
            pert_col=args.pert_col,
            ntc_name=args.ntc_name,
            max_cells=args.max_ctrl_cells,
            seed=args.seed + ctx_idx,
        )
        ctrl_log1p_full, lib_ts = _normalize_counts_log1p(ctrl_csr, target_sum=None)
        ctrl_log1p_shared = ctrl_log1p_full[:, shared_col_idx]
        cov_t = compute_control_covariance(ctrl_log1p_shared, shrink=args.cov_shrink)

        ctx_perts_map = pert_counts.get(ctx, pert_counts.get(None))
        if not ctx_perts_map:
            raise ValueError(f"No target perturbations resolved for context {ctx!r}")
        target_names = list(ctx_perts_map.keys())
        target_gene_idx = [shared_g2i.get(t, -1) for t in target_names]

        if model_state is not None:
            mismatched_cov = sources[0].ctrl_cov if len(sources) > 0 else None
            suite = predict_zero_shot_suite(
                model_state,
                cov_t,
                target_names,
                target_gene_idx,
                mismatched_ctrl_cov=mismatched_cov,
                shuffle_seed=args.seed + ctx_idx,
            )
            if args.method not in suite:
                raise ValueError(f"Unknown method {args.method!r}; available: {sorted(suite.keys())}")
            deltas_shared = suite[args.method]
            rep_used = args.representation if source_mode == "shared_lines" else "log1p"
        else:
            # Pure control-covariance unit-response fallback (zero-shot without source deltas)
            deltas_shared = float(args.default_lambda) * cipher_unit_response_columns(
                cov_t, target_gene_idx
            )
            rep_used = "log1p"

        ratios = _deltas_to_fold_ratios(
            deltas_shared,
            ctrl_log1p_full,
            shared_col_idx,
            representation=rep_used,
            max_abs_delta=args.max_abs_delta,
            max_fold=args.max_fold,
        )

        # Support uniform cells_per_pert (400) or per-target n_cells
        unique_counts = set(ctx_perts_map.values())
        if len(unique_counts) == 1:
            n_emit = next(iter(unique_counts))
            X_ctx = emit_dispersed_counts_for_context(
                ctrl_csr,
                target_names,
                ratios,
                cells_per_pert=n_emit,
                max_counts_per_cell=MAX_COUNTS_PER_CELL,
                seed=args.seed + 1000 * (ctx_idx + 1),
            )
            for t_name in target_names:
                obs_perts.extend([t_name] * n_emit)
                obs_contexts.extend([ctx] * n_emit)
        else:
            sub_blocks = []
            for p_i, t_name in enumerate(target_names):
                n_emit = int(ctx_perts_map[t_name])
                sub_blk = emit_dispersed_counts_for_context(
                    ctrl_csr,
                    [t_name],
                    ratios[p_i : p_i + 1],
                    cells_per_pert=n_emit,
                    max_counts_per_cell=MAX_COUNTS_PER_CELL,
                    seed=args.seed + 1000 * (ctx_idx + 1) + p_i,
                )
                sub_blocks.append(sub_blk)
                obs_perts.extend([t_name] * n_emit)
                obs_contexts.extend([ctx] * n_emit)
            X_ctx = sp.vstack(sub_blocks, format="csr")

        X_blocks.append(X_ctx)
        per_context_summary[ctx] = {
            "n_ctrl_cells_used": int(ctrl_csr.shape[0]),
            "ctrl_median_umi": round(lib_ts, 2),
            "n_targets": len(target_names),
            "n_targets_in_shared_genes": int(sum(1 for x in target_gene_idx if x >= 0)),
            "n_emitted_cells": int(X_ctx.shape[0]),
            "emitted_nnz": int(X_ctx.nnz),
        }

    X_all = sp.vstack(X_blocks, format="csr")
    X_all.data = X_all.data.astype(np.float32, copy=False)
    X_all.eliminate_zeros()
    X_all.sort_indices()

    obs_df = pd.DataFrame(
        {
            DEFAULT_PERT_COL: obs_perts,
            DEFAULT_CONTEXT_COL: obs_contexts,
        },
        index=[str(i) for i in range(X_all.shape[0])],
    )
    var_df = pd.DataFrame(index=list(panel_genes))
    pred_ad = ad.AnnData(X=X_all, obs=obs_df, var=var_df)

    # 3. Validate strictly against VCC 2026 rules
    exp_gene_dim = None if args.simulation else EXPECTED_GENE_DIM
    max_cell_dim = None if args.simulation else MAX_CELL_DIM
    report = validate_vcc2026_anndata(
        pred_ad,
        gene_names=panel_genes,
        pert_counts=pert_counts,
        required_contexts=contexts,
        cells_per_pert=args.cells_per_pert,
        expected_gene_dim=exp_gene_dim,
        max_cell_dim=max_cell_dim,
        max_nnz=MAX_NNZ,
        max_counts_per_cell=MAX_COUNTS_PER_CELL,
        reject_controls=True,
        require_sparse=True,
    )
    if not report.passed:
        raise RuntimeError(f"VCC 2026 format validation FAILED: {report.errors}")

    # 4. Save prediction.h5ad, optional .vcc, SHA-256 checksums, and timestamp declaration
    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    h5ad_out = outdir / "prediction.h5ad"
    pred_ad.write_h5ad(h5ad_out)
    output_files = [h5ad_out]

    if args.package_vcc:
        vcc_out = outdir / "prediction.prep.vcc"
        perts_arg = input_files.get("pert_counts")
        package_vcc_archive(
            pred_ad,
            vcc_out,
            contexts=contexts,
            perts_path=perts_arg,
            panel_id=panel_id,
        )
        output_files.append(vcc_out)

    method_meta = {
        "method": args.method if model_state is not None else "Control_Covariance_Fallback",
        "source_mode": source_mode,
        "source_lines": [s.name for s in sources],
        "fitted_cipher_lambda": round(model_state.lam, 6) if model_state is not None else float(args.default_lambda),
        "n_shared_genes": len(shared_genes),
        "n_panel_genes": len(panel_genes),
        "simulation_mode": bool(args.simulation),
        "seed": int(args.seed),
        "per_context": per_context_summary,
    }
    manifest = write_timestamp_bundle(
        outdir,
        output_files=output_files,
        input_files=input_files,
        validation_report=report,
        method_meta=method_meta,
    )
    return manifest


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="End-to-end VCC 2026 prospective prediction & SHA-256 packager")
    p.add_argument("--bundle-dir", type=Path, default=None,
                   help="Directory with manifest.json, gene_names.csv, pert_counts.csv, context_*.h5ad")
    p.add_argument("--controls", action="append", default=None,
                   help="Context control h5ad as CTX=PATH (repeatable, alternative to --bundle-dir)")
    p.add_argument("--genes", type=Path, default=None, help="Path to gene_names.csv")
    p.add_argument("--perts", type=Path, default=None, help="Path to pert_counts.csv or targets.txt")
    p.add_argument("--contexts", type=str, default=None, help="Comma-separated contexts (e.g. D,E,F)")

    p.add_argument("--source-lines-dir", type=Path, default=None,
                   help="Path to p2/shared/lines/ prepared by A-T1")
    p.add_argument("--source-lines", type=str, default=None,
                   help="Comma-separated subset of source lines to load from --source-lines-dir")
    p.add_argument("--source-h5ad", action="append", default=None,
                   help="Standalone source Perturb-seq h5ad as NAME=PATH (repeatable)")
    p.add_argument("--source-pert-col", type=str, default=DEFAULT_PERT_COL)
    p.add_argument("--representation", type=str, default="log1p", choices=["log1p", "zscore"])
    p.add_argument("--self-control-only", action="store_true",
                   help="Fallback: predict via target control covariance only when no source lines exist")
    p.add_argument("--default-lambda", type=float, default=-0.5,
                   help="Fallback CIPHER lambda used only with --self-control-only")

    p.add_argument("--method", type=str, default="Hybrid_ProjShared",
                   choices=["Hybrid_ProjShared", "Hybrid_Raw", "CIPHER_Linear",
                            "Transfer_Mean", "Source_Grand_Mean"])
    p.add_argument("--pert-col", type=str, default=DEFAULT_PERT_COL)
    p.add_argument("--ntc-name", type=str, default=DEFAULT_NTC_NAME)
    p.add_argument("--cells-per-pert", type=int, default=EXPECTED_CELLS_PER_PERT)
    p.add_argument("--max-ctrl-cells", type=int, default=20000)
    p.add_argument("--min-source-cells", type=int, default=5)
    p.add_argument("--cov-shrink", type=float, default=0.0)
    p.add_argument("--max-abs-delta", type=float, default=4.0)
    p.add_argument("--max-fold", type=float, default=50.0)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--simulation", action="store_true",
                   help="Relax fixed 18,533 gene dim / 400k total cell cap for non-VCC-2026 test panels")
    p.add_argument("--package-vcc", action="store_true",
                   help="Also package .vcc tar archive (requires zstandard)")
    p.add_argument("--outdir", type=Path, required=True)
    return p


if __name__ == "__main__":
    parsed = build_parser().parse_args()
    res = run_prospective(parsed)
    print(json.dumps({
        "created_utc": res["created_utc"],
        "method": res["method"]["method"],
        "validation_passed": res["validation"]["passed"],
        "n_cells": res["validation"]["n_cells"],
        "n_genes": res["validation"]["n_genes"],
        "nnz": res["validation"]["nnz"],
        "outputs": res["outputs"],
    }, indent=2))
