#!/usr/bin/env python3
"""
Verify local `ArcInstitute/cell-eval2` (`v0.22.0`) reproduction on public VCC 2025 data
(`H1-VCC-2025-training.h5ad`, or full `adata_Validation.h5ad` if passed via `--real-h5ad`).

Exercises the exact end-to-end competition workflow:
  1. Synthesize two test predictions in VCC 2026 upload format (no `non-targeting` rows):
       - `pred_signal`: dispersed perturbed count profiles
       - `pred_null`: unperturbed control resample (no perturbation specificity)
  2. Validate upload format via `vcc2026_format.validate_vcc2026_anndata`.
  3. Bridge upload -> `cell-eval2` scoring format via `attach_real_controls_for_celleval2`.
  4. Execute `cell-eval2` CLI commands (`baseline`, `run --preset vcc2026 --anchor`,
     `prep-real-bundle`, `score --real-bundle`, and `run --profile vcc`).
  5. Save verified outputs & summary JSON.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

PROSPECTIVE_DIR = Path(__file__).resolve().parent
if str(PROSPECTIVE_DIR) not in sys.path:
    sys.path.insert(0, str(PROSPECTIVE_DIR))

from vcc2026_format import (
    DEFAULT_NTC_NAME,
    DEFAULT_PERT_COL,
    MAX_COUNTS_PER_CELL,
    MAX_NNZ,
    attach_real_controls_for_celleval2,
    emit_dispersed_counts_for_context,
    validate_vcc2026_anndata,
)

DEFAULT_VCC2025_SAMPLE = Path(
    "/tmp/round9_S4/repos/cell-eval2/cell-eval2-main/docs/data/H1-VCC-2025-training.h5ad"
)


def _run_cmd(cmd: List[str]) -> str:
    proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return proc.stdout.strip()


def _find_celleval2_bin(explicit: str | None) -> str:
    if explicit:
        return explicit
    cand = Path(sys.executable).parent / "cell-eval2"
    if cand.exists():
        return str(cand)
    which = shutil.which("cell-eval2")
    if which:
        return which
    raise FileNotFoundError("Could not find `cell-eval2` CLI binary. Pass --cell-eval2-bin.")


def verify_on_vcc2025(
    real_h5ad: Path,
    outdir: Path,
    cell_eval2_bin: str,
    *,
    pert_col: str = DEFAULT_PERT_COL,
    control: str = DEFAULT_NTC_NAME,
    max_perts: int | None = None,
    seed: int = 42,
) -> Dict[str, Any]:
    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    real_ad = ad.read_h5ad(real_h5ad)
    real_ad.var_names = real_ad.var_names.astype(str)
    genes = real_ad.var_names.tolist()

    labels = real_ad.obs[pert_col].astype(str).to_numpy()
    perts = sorted(set(labels.tolist()) - {control})
    if max_perts is not None and len(perts) > max_perts:
        perts = perts[:max_perts]
        keep_mask = np.isin(labels, [control] + perts)
        real_sub = real_ad[keep_mask].copy()
        real_use_path = outdir / "real_subset.h5ad"
        real_sub.write_h5ad(real_use_path)
        real_ad = real_sub
        labels = real_ad.obs[pert_col].astype(str).to_numpy()
    else:
        real_use_path = real_h5ad

    ctrl_mask = labels == control
    ctrl_X = real_ad.X[ctrl_mask]
    ctrl_csr = ctrl_X.tocsr() if sp.issparse(ctrl_X) else sp.csr_matrix(ctrl_X)

    # Normalize control & perturbed pseudobulk means to compute ground-truth fold ratios
    def _norm_mean(mat):
        arr = mat.toarray().astype(np.float64) if sp.issparse(mat) else np.asarray(mat, dtype=np.float64)
        tot = arr.sum(axis=1)
        nz = tot > 0
        arr_n = np.zeros_like(arr)
        arr_n[nz] = arr[nz] * (1e6 / tot[nz, None])
        return arr_n.mean(axis=0)

    ctrl_mean = _norm_mean(ctrl_csr)
    ratios_signal = np.ones((len(perts), len(genes)), dtype=np.float64)
    obs_perts: List[str] = []
    n_cells_map: Dict[str, int] = {}

    for i, p_name in enumerate(perts):
        p_mask = labels == p_name
        n_p = int(p_mask.sum())
        n_cells_map[p_name] = n_p
        obs_perts.extend([p_name] * n_p)
        p_mean = _norm_mean(real_ad.X[p_mask])
        pos = ctrl_mean > 1e-8
        r = np.ones(len(genes), dtype=np.float64)
        r[pos] = np.clip(p_mean[pos] / ctrl_mean[pos], 0.0, 50.0)
        ratios_signal[i] = r

    ratios_null = np.ones_like(ratios_signal)

    # Build upload-format predictions (non-targeting forbidden)
    # Handle variable cell counts per pert in VCC 2025 real data
    sig_blocks, null_blocks = [], []
    for i, p_name in enumerate(perts):
        n_p = n_cells_map[p_name]
        sig_blocks.append(
            emit_dispersed_counts_for_context(
                ctrl_csr, [p_name], ratios_signal[i : i + 1],
                cells_per_pert=n_p, max_counts_per_cell=MAX_COUNTS_PER_CELL, seed=seed + i,
            )
        )
        null_blocks.append(
            emit_dispersed_counts_for_context(
                ctrl_csr, [p_name], ratios_null[i : i + 1],
                cells_per_pert=n_p, max_counts_per_cell=MAX_COUNTS_PER_CELL, seed=seed + 1000 + i,
            )
        )

    X_sig = sp.vstack(sig_blocks, format="csr")
    X_null = sp.vstack(null_blocks, format="csr")
    obs_df = pd.DataFrame(
        {pert_col: obs_perts},
        index=[f"pred_{i}" for i in range(len(obs_perts))],
    )
    var_df = pd.DataFrame(index=genes)

    pert_dict: Dict[str | None, Dict[str, int]] = {None: n_cells_map}

    sig_upload_path = outdir / "pred_signal_upload.h5ad"
    null_upload_path = outdir / "pred_null_upload.h5ad"
    sig_ad = ad.AnnData(X=X_sig, obs=obs_df.copy(), var=var_df.copy())
    null_ad = ad.AnnData(X=X_null, obs=obs_df.copy(), var=var_df.copy())

    rep_sig = validate_vcc2026_anndata(
        sig_ad, gene_names=genes, pert_counts=pert_dict,
        required_contexts=None, expected_gene_dim=None, max_cell_dim=None,
        max_nnz=MAX_NNZ, max_counts_per_cell=MAX_COUNTS_PER_CELL,
        reject_controls=True, require_sparse=True,
    )
    rep_null = validate_vcc2026_anndata(
        null_ad, gene_names=genes, pert_counts=pert_dict,
        required_contexts=None, expected_gene_dim=None, max_cell_dim=None,
        max_nnz=MAX_NNZ, max_counts_per_cell=MAX_COUNTS_PER_CELL,
        reject_controls=True, require_sparse=True,
    )
    if not rep_sig.passed or not rep_null.passed:
        raise RuntimeError(f"Upload validation failed: sig={rep_sig.errors}, null={rep_null.errors}")

    sig_ad.write_h5ad(sig_upload_path)
    null_ad.write_h5ad(null_upload_path)

    # Bridge upload -> scoring format (attach real non-targeting controls)
    sig_score_path = attach_real_controls_for_celleval2(
        sig_upload_path, real_use_path, outdir / "pred_signal_for_celleval2.h5ad",
        pert_col=pert_col, control_label=control,
    )
    null_score_path = attach_real_controls_for_celleval2(
        null_upload_path, real_use_path, outdir / "pred_null_for_celleval2.h5ad",
        pert_col=pert_col, control_label=control,
    )

    # Execute official cell-eval2 CLI pipeline (vcc2026 preset on CPU + pdex)
    common_flags = [
        "--preset", "vcc2026",
        "--pert-col", pert_col,
        "--control", control,
        "--set", "device=cpu",
        "--set", "de.backend=pdex",
    ]
    base_dir = outdir / "baseline"
    sig_dir = outdir / "run_signal_vcc2026"
    null_dir = outdir / "run_null_vcc2026"
    bundle_dir = outdir / "real_bundle"
    vcc2025_dir = outdir / "run_signal_vcc2025_profile"

    base_pred_path = base_dir / "baseline_pred.h5ad"

    _run_cmd([
        cell_eval2_bin, "baseline",
        "-ar", str(real_use_path),
        "-o", str(base_dir),
        "--save-pred", str(base_pred_path),
        *common_flags,
    ])
    _run_cmd([
        cell_eval2_bin, "run", "-ap", str(sig_score_path), "-ar", str(real_use_path),
        "-o", str(sig_dir), "--anchor", "--anchor-splits", "5", "--anchor-base-seed", "0",
        *common_flags,
    ])
    _run_cmd([
        cell_eval2_bin, "run", "-ap", str(null_score_path), "-ar", str(real_use_path),
        "-o", str(null_dir), *common_flags,
    ])
    _run_cmd([
        cell_eval2_bin, "prep-real-bundle",
        "--real", str(real_use_path),
        "--baseline", str(base_pred_path),
        "-o", str(bundle_dir),
        "--id", "vcc2025-verify-bundle",
        "--force",
        *common_flags,
    ])
    _run_cmd([
        cell_eval2_bin, "score",
        "--user-agg", str(sig_dir / "agg_results.csv"),
        "--real-bundle", str(bundle_dir),
        "-o", str(sig_dir / "scored.csv"),
    ])
    _run_cmd([
        cell_eval2_bin, "score",
        "--user-agg", str(null_dir / "agg_results.csv"),
        "--real-bundle", str(bundle_dir),
        "-o", str(null_dir / "scored.csv"),
    ])
    _run_cmd([
        cell_eval2_bin, "run", "-ap", str(sig_score_path), "-ar", str(real_use_path),
        "-o", str(vcc2025_dir),
        "--profile", "vcc", "--pert-col", pert_col, "--control", control,
        "--input-type", "counts", "--set", "device=cpu", "--set", "de.backend=pdex",
    ])

    sig_scored_df = pd.read_csv(sig_dir / "scored.csv")
    null_scored_df = pd.read_csv(null_dir / "scored.csv")
    sig_vcc25_df = pd.read_csv(vcc2025_dir / "agg_results.csv")
    eval_doc = json.loads((sig_dir / "evaluation_results.json").read_text(encoding="utf-8"))

    sig_avg = float(sig_scored_df.loc[sig_scored_df["metric"] == "avg_score", "from_replicate"].iloc[0])
    null_avg = float(null_scored_df.loc[null_scored_df["metric"] == "avg_score", "from_replicate"].iloc[0])
    if not (sig_avg > null_avg):
        raise AssertionError(f"Expected signal avg_score ({sig_avg}) > null avg_score ({null_avg})")

    summary = {
        "real_h5ad": str(real_use_path),
        "real_shape": [int(real_ad.n_obs), int(real_ad.n_vars)],
        "n_perturbations": len(perts),
        "n_control_cells": int(ctrl_mask.sum()),
        "upload_validation_signal": rep_sig.to_dict(),
        "upload_validation_null": rep_null.to_dict(),
        "vcc2026_scored_signal": sig_scored_df.to_dict(orient="records"),
        "vcc2026_scored_null": null_scored_df.to_dict(orient="records"),
        "vcc2025_profile_signal_agg": sig_vcc25_df.to_dict(orient="records"),
        "evaluation_results_suite_settings": eval_doc.get("suite", {}).get("settings"),
        "signal_avg_score_from_replicate": sig_avg,
        "null_avg_score_from_replicate": null_avg,
    }
    (outdir / "celleval2_vcc2025_verify.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    ap = argparse.ArgumentParser(description="Verify cell-eval2 locally on public VCC 2025 data")
    ap.add_argument("--real-h5ad", type=Path, default=DEFAULT_VCC2025_SAMPLE)
    ap.add_argument("--outdir", type=Path, default=Path("/tmp/round9_A_T4/celleval2_vcc2025_verify"))
    ap.add_argument("--cell-eval2-bin", type=str, default=None)
    ap.add_argument("--pert-col", type=str, default=DEFAULT_PERT_COL)
    ap.add_argument("--control", type=str, default=DEFAULT_NTC_NAME)
    ap.add_argument("--max-perts", type=int, default=None)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    ce_bin = _find_celleval2_bin(args.cell_eval2_bin)
    summary = verify_on_vcc2025(
        args.real_h5ad,
        args.outdir,
        ce_bin,
        pert_col=args.pert_col,
        control=args.control,
        max_perts=args.max_perts,
        seed=args.seed,
    )
    print(json.dumps({
        "real_shape": summary["real_shape"],
        "n_perturbations": summary["n_perturbations"],
        "signal_avg_score_from_replicate": summary["signal_avg_score_from_replicate"],
        "null_avg_score_from_replicate": summary["null_avg_score_from_replicate"],
        "vcc2026_scored_signal": summary["vcc2026_scored_signal"],
        "vcc2026_scored_null": summary["vcc2026_scored_null"],
    }, indent=2))


if __name__ == "__main__":
    main()
