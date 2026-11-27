#!/usr/bin/env python3
"""
Clean-Room Round 16 Addendum B (`Occam`) Full Evaluation & Stage-4 Array Recomputation Driver
(`round16_addB/audit/run_audit_addb.py`, `A-T3`).

Implements `ADDB_ARRAY_SCHEMA.md` and `PREREG_I_ADDENDUM_B_DRAFT_v2.md` (`L1`/`L2`, `CPU float32`,
`torch.set_num_threads(4)`) cleanly across three CLI modes:
  1. `--phase real_addb`: Hard-gated by `R16_ADDB_CODE_FREEZE_AT3.txt` + `R16_ADDB_GO_AT3.txt`
     + Gate 5A/5B check; exports exact `control_barcodes_addb.json`, `z_ctrl/*.npz`,
     `arrays/<line>/{observed,preds_*}.npz` (`allow_pickle=False`), `summary_addb.json`,
     and verifies Stage-4 same-array recomputation (`<= 1e-10` contract / `<= 1e-12` target).
  2. `--phase synthetic_e2e`: Hermetic CPU synthetic end-to-end run over all 5 target lines & 7 pools
     exercising exact `ADDB_ARRAY_SCHEMA.md` schema, `pairwise_common_axis` (`n_union_absent > 0`
     vs `reused_primary_rhj_6983`), Tier-2 `RBM14` `402`-target sensitivity, `gse295214_h1` extra
     `FCT` sensitivity, `allow_pickle=False`, and Stage-4 recomputation parity (`0.0 <= 1e-12`).
  3. `--recompute-from-arrays <results_dir>`: Recomputes `summary_addb.json` strictly from exported
     `.npz` arrays + committed R16 `STATE_SE` bundles and verifies `max_abs_diff <= 1e-10`.

Zero reads of `round16_addB/src/*.py`, `round16_addB/results/`, `round16/results/group3/`, or `MAIN.txt`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import torch

AUDIT_ADDB_DIR = Path(__file__).resolve().parent
ADDB_ROOT = AUDIT_ADDB_DIR.parent
REPO_ROOT = ADDB_ROOT.parent
R16_AUDIT_DIR = REPO_ROOT / "round16" / "audit"
R10_SRC_DIR = REPO_ROOT / "round10" / "src"

for _p in (str(AUDIT_ADDB_DIR), str(R16_AUDIT_DIR), str(R10_SRC_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from audit_addb_controls_and_metrics import (  # noqa: E402
    ADDB_ALL_CONTROL_POOLS,
    ADDB_ALL_TARGET_LINES,
    ADDB_LINE_ROLES,
    ADDB_SEEN_SOURCES,
    ADDB_UNSEEN_SOURCES,
    DEFAULT_GSE295214_RAW_DIR,
    DEFAULT_LOCKED_MAIN_R16_DIR,
    DEFAULT_LOCKED_MAIN_SENS_DIR,
    DEFAULT_LOCKED_SE_AUDIT_DIR,
    AddBBlindingGateError,
    CommittedR16ReferenceBundle,
    ExtractedRawControlCells,
    assert_freeze_gates_addb,
    compute_group_adjudications_addb,
    evaluate_addb_line_complete,
    extract_all_required_raw_pools_addb,
    format_control_barcodes_addb_json,
    load_committed_r16_reference_and_se,
    load_locked_control_barcodes_addb,
)
from audit_addb_occam_adapter import (  # noqa: E402
    DEFAULT_CKPT_PATH,
    DEFAULT_GENE_AXIS_CSV,
    DEFAULT_NADIG_RAW_DIR,
    DEFAULT_PERTURBATIONS_CSV,
    EXPECTED_BEST_MODEL_SHA256,
    OCCAM_NUM_ENCODER_GENES,
    OCCAM_NUM_RHJ_AXIS_GENES,
    OCCAM_TARGET_SUM,
    OccamCleanRoomAdapter,
    TwoTierPertResolver,
    compute_pairwise_common_rhj_mask,
    load_hgnc_bridge,
    load_rhj_eligible_perturbations,
    load_rhj_gene_axis,
    normalize_counts_full_44728_api,
    normalize_counts_to_rhj_44728,
    resolve_output_projection,
    sha256_array_bytes,
    sha256_file,
    verify_addb_checksums,
)
from audit_r16_metrics import DEFAULT_R16_SEED  # noqa: E402
from run_addb_gates_cpu import run_stage_5b_cpu  # noqa: E402

LOCKED_DEVICE = "cpu"
LOCKED_NUM_THREADS = 4
OCCAM_MODEL_KEY = "Occam"
OCCAM_CHECKPOINT_ID = "rhj_fewshot/all_sources/seed_42"
RESAMPLE_SEED = 20261602
EXPECTED_AXIS_LEN_REAL = OCCAM_NUM_RHJ_AXIS_GENES


def _log_progress(msg: str) -> None:
    """Emit timestamped line with process CPU seconds per Addendum B §6.1 10-min heartbeat discipline."""
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    t = os.times()
    print(f"[{ts}] [cpu_u={t.user:.1f}s cpu_s={t.system:.1f}s] {msg}", flush=True)


def _to_strict_npz_array(val: Any) -> np.ndarray:
    """Convert array or scalar into strict non-object NumPy array compatible with `allow_pickle=False`."""
    if isinstance(val, np.ndarray):
        if val.dtype == object:
            return np.asarray([str(x) for x in val.tolist()], dtype=str)
        return np.ascontiguousarray(val)
    if isinstance(val, (list, tuple)):
        if len(val) == 0:
            return np.asarray([], dtype="U1")
        if isinstance(val[0], str):
            return np.asarray([str(x) for x in val], dtype=str)
        return np.asarray(val)
    if isinstance(val, str):
        return np.asarray(val, dtype=str)
    if isinstance(val, (bool, np.bool_)):
        return np.asarray(bool(val), dtype=bool)
    if isinstance(val, (int, np.integer)):
        return np.asarray(int(val), dtype=np.int64)
    if isinstance(val, (float, np.floating)):
        return np.asarray(float(val), dtype=np.float64)
    raise TypeError(f"Unsupported type for strict allow_pickle=False npz export: {type(val)!r}")


def write_npz_strict_addb(out_path: Path, payload: Mapping[str, Any]) -> Path:
    """Write compressed `.npz` with zero `dtype=object` arrays and verify `allow_pickle=False` reload."""
    p = Path(out_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    strict_map: Dict[str, np.ndarray] = {k: _to_strict_npz_array(v) for k, v in payload.items()}
    for k, arr in strict_map.items():
        if arr.dtype == object:
            raise TypeError(f"Refusing object dtype array for key {k!r} in {p}")
    np.savez_compressed(p, **strict_map)
    with np.load(p, allow_pickle=False) as chk:
        _ = list(chk.files)
    return p


def _json_safe(obj: Any) -> Any:
    """Recursively convert NumPy scalars/arrays and handle non-finite floats cleanly for JSON serialization."""
    if isinstance(obj, Mapping):
        return {str(k): _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(x) for x in obj]
    if isinstance(obj, np.ndarray):
        return _json_safe(obj.tolist())
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if isinstance(obj, (np.floating, float)):
        v = float(obj)
        if not math.isfinite(v):
            return None
        return v
    if isinstance(obj, Path):
        return str(obj)
    return obj


def _write_json(path: Path, payload: Any) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(_json_safe(payload), indent=2, sort_keys=True), encoding="utf-8")
    return p


def build_occam_delta_payload_addb(
    delta_pred: np.ndarray,
    target_names: Sequence[str],
    gene_names: Sequence[str],
    resolved_pert_keys: Sequence[str],
    resolved_pert_ids: Sequence[int],
    input_axis_mode: str,
    baseline_mode: str,
    control_source: str,
    control_barcodes: Sequence[str],
    pair_counterpart: str,
    n_active_rhj_genes: int,
    axis_mask_sha256: str,
) -> Dict[str, Any]:
    """Build exact 15-key `preds_*.npz` payload matching `ADDB_ARRAY_SCHEMA.md §4.1`."""
    dp = np.ascontiguousarray(np.asarray(delta_pred, dtype=np.float64))
    p_len, g_len = len(target_names), len(gene_names)
    if dp.shape != (p_len, g_len):
        raise ValueError(f"delta_pred shape {dp.shape} != expected ({p_len}, {g_len})")
    if len(resolved_pert_keys) != p_len or len(resolved_pert_ids) != p_len:
        raise ValueError("resolved_pert_keys/ids length mismatch vs target_names")
    return {
        "delta_pred": dp,
        "target_names": [str(x) for x in target_names],
        "gene_names": [str(x) for x in gene_names],
        "resolved_pert_keys": [str(x) for x in resolved_pert_keys],
        "resolved_pert_ids": np.asarray(resolved_pert_ids, dtype=np.int64),
        "model_key": OCCAM_MODEL_KEY,
        "checkpoint_id": OCCAM_CHECKPOINT_ID,
        "input_axis_mode": str(input_axis_mode),
        "baseline_mode": str(baseline_mode),
        "control_source": str(control_source),
        "control_barcodes": [str(x) for x in control_barcodes],
        "pair_counterpart": str(pair_counterpart),
        "n_active_rhj_genes": int(n_active_rhj_genes),
        "axis_mask_sha256": str(axis_mask_sha256),
        "delta_sha256": sha256_array_bytes(dp),
    }


def export_observed_npz_addb(
    out_dir: Path,
    ref: CommittedR16ReferenceBundle,
) -> Dict[str, str]:
    """Export exact 18-key `observed.npz` (plus Group-3 `preds_SENS_FCT_*` if `gse295214_h1`)."""
    files_written: Dict[str, str] = {}
    mu_t = ref.mu_target if ref.mu_target is not None else np.zeros(len(ref.eval_genes), dtype=np.float64)
    mu_s = (
        ref.mu_sources
        if ref.mu_sources is not None
        else np.zeros((5 if ref.target_line == "gse295214_h1" else 6, len(ref.eval_genes)), dtype=np.float64)
    )
    e1_d = ref.e1_constmed if ref.e1_constmed is not None else np.zeros_like(ref.delta_obs, dtype=np.float64)

    obs_payload: Dict[str, Any] = {
        "delta_obs": np.asarray(ref.delta_obs, dtype=np.float64),
        "delta_half_a": np.asarray(ref.delta_half_a, dtype=np.float64),
        "delta_half_b": np.asarray(ref.delta_half_b, dtype=np.float64),
        "half_valid": np.asarray(ref.half_valid, dtype=bool),
        "delta_transfer": np.asarray(ref.delta_transfer, dtype=np.float64),
        "source_grand_mean": np.asarray(ref.source_grand_mean, dtype=np.float64),
        "fct_delta_pred": np.asarray(ref.delta_fct, dtype=np.float64),
        "e1_constmed_delta_pred": np.asarray(e1_d, dtype=np.float64),
        "target_indices": np.asarray(ref.target_gene_idx, dtype=np.int64),
        "target_names": list(ref.eval_perts),
        "gene_names": list(ref.eval_genes),
        "mu_target": np.asarray(mu_t, dtype=np.float64),
        "mu_sources": np.asarray(mu_s, dtype=np.float64),
        "matched_ctrl_barcodes": list(ref.matched_barcodes),
        "resample_ctrl_barcodes": list(ref.resample_barcodes),
        "n_ctrl_total": int(ref.n_ctrl_total),
        "n_ctrl_sampled": int(len(ref.matched_barcodes)),
        "ctrl_sampling_fraction": float(ref.ctrl_sampling_fraction),
    }
    obs_fp = write_npz_strict_addb(Path(out_dir) / "observed.npz", obs_payload)
    files_written[obs_fp.name] = sha256_file(obs_fp)

    if ref.target_line == "gse295214_h1" and ref.extra_fct_full_r9_6 is not None:
        ef = ref.extra_fct_full_r9_6
        sfp = write_npz_strict_addb(
            Path(out_dir) / "preds_SENS_FCT_full_r9_6_incl_vcc2025_h1.npz",
            {
                "delta_pred": np.asarray(ef["delta_pred"], dtype=np.float64),
                "delta_transfer": np.asarray(ef["delta_transfer"], dtype=np.float64),
                "source_grand_mean": np.asarray(ef["source_grand_mean"], dtype=np.float64),
                "target_names": list(ref.eval_perts),
                "gene_names": list(ref.eval_genes),
                "model_key": "FCT_full_r9_6_incl_vcc2025_h1",
            },
        )
        files_written[sfp.name] = sha256_file(sfp)

    return files_written


def export_z_ctrl_latent_npz(
    z_ctrl_dir: Path,
    pool_key: str,
    z_ctrl: np.ndarray,
    control_barcodes: Sequence[str],
    seed: int,
    pool_name: str,
    input_axis_mode: str,
    n_active_rhj_genes: int,
    axis_mask_sha256: str,
    audit_meta: Mapping[str, Any],
) -> Tuple[str, str]:
    """Export one 9-key `z_ctrl_<pool_key>.npz` file (`ADDB_ARRAY_SCHEMA.md §5.1`)."""
    zc = np.ascontiguousarray(np.asarray(z_ctrl, dtype=np.float32))
    z_sha = sha256_array_bytes(zc)
    full_audit = dict(audit_meta)
    full_audit["z_ctrl_sha256"] = z_sha
    fname = f"z_ctrl_{pool_key}.npz"
    fp = write_npz_strict_addb(
        Path(z_ctrl_dir) / fname,
        {
            "Z_ctrl": zc,
            "control_barcodes": [str(x) for x in control_barcodes],
            "seed": int(seed),
            "pool_name": str(pool_name),
            "input_axis_mode": str(input_axis_mode),
            "n_active_rhj_genes": int(n_active_rhj_genes),
            "axis_mask_sha256": str(axis_mask_sha256),
            "z_ctrl_sha256": z_sha,
            "audit_json": json.dumps(_json_safe(full_audit), sort_keys=True),
        },
    )
    return fname, sha256_file(fp)


def _predict_deltas_both_baselines_from_z(
    adapter: Any,
    x_44728: np.ndarray,
    z_ctrl: np.ndarray,
    pert_ids: Sequence[int],
    out_cols: np.ndarray,
    batch_size: int = 256,
    compute_real_baseline: bool = False,
) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """
    Run ONE forward decoding pass per perturbation over cached `z_ctrl` (`(512, 512)`),
    producing `baseline="reconstructed"` AND optionally `baseline="real"` with zero extra passes.
    """
    n_cells = int(z_ctrl.shape[0])
    cols = np.asarray(out_cols, dtype=np.int64)
    ctrl_sub = np.asarray(x_44728[:, cols], dtype=np.float32)

    base_recon = adapter.predict_reconstructed_baseline(z_ctrl, batch_size=batch_size, out_cols=cols)
    base_recon_mean = base_recon.mean(axis=0, dtype=np.float64)
    base_real_mean = ctrl_sub.mean(axis=0, dtype=np.float64) if compute_real_baseline else None

    deltas_recon = np.empty((len(pert_ids), cols.size), dtype=np.float64)
    deltas_real = np.empty((len(pert_ids), cols.size), dtype=np.float64) if compute_real_baseline else None

    for row_i, pid in enumerate(pert_ids):
        pred_cells = adapter.predict_per_cell(
            z_ctrl,
            [[int(pid)]] * n_cells,
            batch_size=batch_size,
            out_cols=cols,
        )
        bad = ~np.isfinite(pred_cells)
        if bad.any():
            pred_cells[bad] = ctrl_sub[bad]
        mean_pred = pred_cells.mean(axis=0, dtype=np.float64)
        deltas_recon[row_i] = mean_pred - base_recon_mean
        if deltas_real is not None and base_real_mean is not None:
            deltas_real[row_i] = mean_pred - base_real_mean

    return deltas_recon, deltas_real


def _mask_sha256_uint8(mask_bool: Sequence[bool]) -> str:
    arr = np.ascontiguousarray(np.asarray(mask_bool, dtype=np.uint8))
    return sha256_array_bytes(arr)


def _format_latent_audit(norm_audit: Mapping[str, Any], raw_pool: ExtractedRawControlCells) -> Dict[str, Any]:
    return {
        "n_present_rhj": int(norm_audit["n_rhj_present_in_raw"]),
        "n_absent_rhj": int(norm_audit["n_rhj_absent_in_raw"]),
        "n_present_but_zeroed_by_mask": int(norm_audit["n_rhj_masked_out"]),
        "n_active_injected": int(norm_audit["n_active_injected"]),
        "n_zero_libsize_cells": int(np.sum(np.asarray(raw_pool.row_sum_umis, dtype=np.float64) <= 0.0)),
    }


def _run_export_pipeline_core(
    adapter: Any,
    rhj_axis_ensembl: Sequence[str],
    rhj_axis_symbols: Sequence[str],
    hgnc_bridge: Any,
    pert_resolver: Any,
    ref_bundles: Mapping[str, CommittedR16ReferenceBundle],
    locked_barcodes: Mapping[str, Any],
    primary_pools: Mapping[str, ExtractedRawControlCells],
    resample_pools: Mapping[str, ExtractedRawControlCells],
    output_dir: Path,
    gates_verified: Mapping[str, str],
    device: str = LOCKED_DEVICE,
    num_threads: int = LOCKED_NUM_THREADS,
    batch_size: int = 256,
    n_boot: int = 4000,
    seed: int = DEFAULT_R16_SEED,
    perturbations_csv: Path = DEFAULT_PERTURBATIONS_CSV,
) -> Dict[str, Any]:
    """
    Shared clean-room export + evaluation pipeline invoked identically by both
    `--phase real_addb` and `--phase synthetic_e2e` (`ADDB_ARRAY_SCHEMA.md`).
    """
    out_dir = Path(output_dir)
    z_ctrl_dir = out_dir / "z_ctrl"
    arrays_root = out_dir / "arrays"
    out_dir.mkdir(parents=True, exist_ok=True)
    z_ctrl_dir.mkdir(parents=True, exist_ok=True)
    arrays_root.mkdir(parents=True, exist_ok=True)

    # 1. Control barcodes JSON
    cb_json = format_control_barcodes_addb_json(locked_barcodes)
    _write_json(out_dir / "control_barcodes_addb.json", cb_json)
    _log_progress("Exported control_barcodes_addb.json")

    enc_genes = list(adapter.gene_names)
    target_sum = float(adapter.target_sum)
    n_rhj_total = len(rhj_axis_ensembl)
    all_ones_mask = np.ones(n_rhj_total, dtype=bool)
    all_ones_sha = _mask_sha256_uint8(all_ones_mask)

    z_manifest: Dict[str, str] = {}
    # Cached tensors & latents for reuse whenever mask zeroes 0 extra genes (`ADDB_ARRAY_SCHEMA §5.2`)
    primary_rhj_cache: Dict[str, Tuple[np.ndarray, np.ndarray, Dict[str, Any]]] = {}
    resample_rhj_cache: Dict[str, Tuple[np.ndarray, np.ndarray, Dict[str, Any]]] = {}
    primary_full_cache: Dict[str, Tuple[np.ndarray, np.ndarray, int]] = {}
    resample_full_cache: Dict[str, Tuple[np.ndarray, np.ndarray, int]] = {}

    # 2A. Encode & save primary + resample pools on rhj_6983_axis
    for pname in sorted(primary_pools.keys()):
        raw_p = primary_pools[pname]
        x44, _, norm_aud = normalize_counts_to_rhj_44728(
            raw_p.raw_counts, raw_p.raw_ensembl_ids, rhj_axis_ensembl, enc_genes,
            target_sum=target_sum, active_rhj_mask=all_ones_mask,
        )
        z = adapter.encode_controls(x44, batch_size=batch_size)
        aud = _format_latent_audit(norm_aud, raw_p)
        fname, fsha = export_z_ctrl_latent_npz(
            z_ctrl_dir, f"rhj_6983_axis__{pname}__primary", z, raw_p.barcodes,
            seed=DEFAULT_R16_SEED, pool_name=pname, input_axis_mode="rhj_6983_axis",
            n_active_rhj_genes=aud["n_active_injected"], axis_mask_sha256=all_ones_sha, audit_meta=aud,
        )
        z_manifest[fname] = fsha
        primary_rhj_cache[pname] = (x44, z, aud)
        _log_progress(f"Encoded Z_ctrl rhj_6983_axis__{pname}__primary (active={aud['n_active_injected']}/{n_rhj_total})")

    for tname in sorted(resample_pools.keys()):
        raw_r = resample_pools[tname]
        x44, _, norm_aud = normalize_counts_to_rhj_44728(
            raw_r.raw_counts, raw_r.raw_ensembl_ids, rhj_axis_ensembl, enc_genes,
            target_sum=target_sum, active_rhj_mask=all_ones_mask,
        )
        z = adapter.encode_controls(x44, batch_size=batch_size)
        aud = _format_latent_audit(norm_aud, raw_r)
        fname, fsha = export_z_ctrl_latent_npz(
            z_ctrl_dir, f"rhj_6983_axis__{tname}__resample", z, raw_r.barcodes,
            seed=RESAMPLE_SEED, pool_name=tname, input_axis_mode="rhj_6983_axis",
            n_active_rhj_genes=aud["n_active_injected"], axis_mask_sha256=all_ones_sha, audit_meta=aud,
        )
        z_manifest[fname] = fsha
        resample_rhj_cache[tname] = (x44, z, aud)

    # 2B. Encode & save unique pairwise_common_axis primary latents (`n_present_but_zeroed_by_mask > 0` only!)
    p4_latent_cache: Dict[Tuple[str, str, str], Tuple[np.ndarray, np.ndarray, int, str, bool]] = {}
    for tname, ref in ref_bundles.items():
        raw_t = primary_pools[tname]
        for sname in ref.seen_sources + ref.unseen_sources:
            raw_s = primary_pools[sname]
            la, lb = (tname, sname) if tname < sname else (sname, tname)
            cmask = compute_pairwise_common_rhj_mask(raw_t.raw_ensembl_ids, raw_s.raw_ensembl_ids, rhj_axis_ensembl)
            n_common = int(np.sum(cmask))
            n_union_abs = int(n_rhj_total - n_common)
            if n_union_abs == 0:
                continue
            cmask_sha = _mask_sha256_uint8(cmask)
            for pname in (tname, sname):
                key = (pname, la, lb)
                if key in p4_latent_cache:
                    continue
                raw_p = primary_pools[pname]
                x44_p, _, norm_aud = normalize_counts_to_rhj_44728(
                    raw_p.raw_counts, raw_p.raw_ensembl_ids, rhj_axis_ensembl, enc_genes,
                    target_sum=target_sum, active_rhj_mask=cmask,
                )
                aud = _format_latent_audit(norm_aud, raw_p)
                if aud["n_present_but_zeroed_by_mask"] == 0:
                    # Bit-for-bit identical input tensor & Z_ctrl to primary rhj_6983_axis; do not duplicate z_ctrl file
                    x44_base, z_base, aud_base = primary_rhj_cache[pname]
                    p4_latent_cache[key] = (x44_base, z_base, aud_base["n_active_injected"], cmask_sha, True)
                else:
                    z_p = adapter.encode_controls(x44_p, batch_size=batch_size)
                    fname, fsha = export_z_ctrl_latent_npz(
                        z_ctrl_dir, f"pairwise_common_axis__{pname}__pair_{la}__{lb}__primary",
                        z_p, raw_p.barcodes, seed=DEFAULT_R16_SEED, pool_name=pname,
                        input_axis_mode="pairwise_common_axis",
                        n_active_rhj_genes=aud["n_active_injected"],
                        axis_mask_sha256=cmask_sha, audit_meta=aud,
                    )
                    z_manifest[fname] = fsha
                    p4_latent_cache[key] = (x44_p, z_p, aud["n_active_injected"], cmask_sha, False)
                    _log_progress(
                        f"Encoded Z_ctrl pairwise_common_axis__{pname}__pair_{la}__{lb}__primary "
                        f"(zeroed={aud['n_present_but_zeroed_by_mask']}, active={aud['n_active_injected']})"
                    )

    # 2C. Encode & save full_44728_api primary + resample pools
    for pname in sorted(primary_pools.keys()):
        raw_p = primary_pools[pname]
        x44, aud_f = normalize_counts_full_44728_api(
            raw_p.raw_counts, raw_p.raw_ensembl_ids, enc_genes, target_sum=target_sum,
        )
        z = adapter.encode_controls(x44, batch_size=batch_size)
        n_match = int(aud_f["n_matched_44728"])
        n_zero = int(np.sum(np.asarray(raw_p.row_sum_umis, dtype=np.float64) <= 0.0))
        fname, fsha = export_z_ctrl_latent_npz(
            z_ctrl_dir, f"full_44728_api__{pname}__primary", z, raw_p.barcodes,
            seed=DEFAULT_R16_SEED, pool_name=pname, input_axis_mode="full_44728_api",
            n_active_rhj_genes=n_match, axis_mask_sha256="",
            audit_meta={"n_matched_44728": n_match, "n_zero_libsize_cells": n_zero},
        )
        z_manifest[fname] = fsha
        primary_full_cache[pname] = (x44, z, n_match)

    for tname in sorted(resample_pools.keys()):
        raw_r = resample_pools[tname]
        x44, aud_f = normalize_counts_full_44728_api(
            raw_r.raw_counts, raw_r.raw_ensembl_ids, enc_genes, target_sum=target_sum,
        )
        z = adapter.encode_controls(x44, batch_size=batch_size)
        n_match = int(aud_f["n_matched_44728"])
        n_zero = int(np.sum(np.asarray(raw_r.row_sum_umis, dtype=np.float64) <= 0.0))
        fname, fsha = export_z_ctrl_latent_npz(
            z_ctrl_dir, f"full_44728_api__{tname}__resample", z, raw_r.barcodes,
            seed=RESAMPLE_SEED, pool_name=tname, input_axis_mode="full_44728_api",
            n_active_rhj_genes=n_match, axis_mask_sha256="",
            audit_meta={"n_matched_44728": n_match, "n_zero_libsize_cells": n_zero},
        )
        z_manifest[fname] = fsha
        resample_full_cache[tname] = (x44, z, n_match)

    _write_json(z_ctrl_dir / "z_ctrl_manifest_sha256.json", z_manifest)
    _log_progress(f"Exported z_ctrl_manifest_sha256.json ({len(z_manifest)} latent files)")

    # 3. Decode & export all arrays per target line, then evaluate
    line_entries: Dict[str, Any] = {}
    r16_se_shas: Dict[str, Dict[str, str]] = {}

    for tname, ref in ref_bundles.items():
        ldir = arrays_root / tname
        ldir.mkdir(parents=True, exist_ok=True)
        arr_manifest = export_observed_npz_addb(ldir, ref)
        if ref.r16_reference_sha256:
            r16_se_shas[tname] = dict(ref.r16_reference_sha256)

        # Resolve Two-Tier pert2id & RHJ_6983 output projection
        pert_ids, pert_rows = pert_resolver.resolve_many(ref.eval_perts)
        pert_keys = [r["resolved_key"] for r in pert_rows]
        t2_targets = [r["eval_target"] for r in pert_rows if r["tier"] != "tier1_rhj2381"]
        bridged = [{"eval": r["eval_target"], "resolved": r["resolved_key"]} for r in pert_rows if r["bridged"]]
        pert_audit = {
            "n_tier1": int(sum(1 for r in pert_rows if r["tier"] == "tier1_rhj2381")),
            "n_tier2": int(len(t2_targets)),
            "tier2_targets": t2_targets,
            "bridged_pairs": bridged,
            "unresolved_targets": [],
        }
        out_proj = resolve_output_projection(
            ref.eval_genes, rhj_axis_ensembl, rhj_axis_symbols, enc_genes, hgnc_bridge=hgnc_bridge,
        )
        if out_proj["n_mapped"] != len(ref.eval_genes) or out_proj["missing_genes"]:
            raise RuntimeError(f"Output projection incomplete on {tname}: missing={out_proj['missing_genes'][:5]}")
        out_cols = out_proj["enc_cols"]
        gene_audit = {
            "n_target_genes": int(len(ref.eval_genes)),
            "n_mapped": int(out_proj["n_mapped"]),
            "n_missing": int(len(out_proj["missing_genes"])),
            "rhj_hgnc_collisions": 0,
        }

        def _write_pred(fname: str, delta_mat: np.ndarray, axis_mode: str, base_mode: str,
                        c_src: str, c_bcs: Sequence[str], counterpart: str, n_act: int, m_sha: str) -> None:
            payload = build_occam_delta_payload_addb(
                delta_pred=delta_mat, target_names=ref.eval_perts, gene_names=ref.eval_genes,
                resolved_pert_keys=pert_keys, resolved_pert_ids=pert_ids,
                input_axis_mode=axis_mode, baseline_mode=base_mode,
                control_source=c_src, control_barcodes=c_bcs,
                pair_counterpart=counterpart, n_active_rhj_genes=n_act, axis_mask_sha256=m_sha,
            )
            fp = write_npz_strict_addb(ldir / fname, payload)
            arr_manifest[fname] = sha256_file(fp)

        # 3A. Axis 1: rhj_6983_axis (`baseline="reconstructed"` + `"real"` zero extra forward passes!)
        x_m, z_m, aud_m = primary_rhj_cache[tname]
        d_m_rec, d_m_real = _predict_deltas_both_baselines_from_z(
            adapter, x_m, z_m, pert_ids, out_cols, batch_size=batch_size, compute_real_baseline=True,
        )
        assert d_m_real is not None
        _write_pred("preds_Occam_rhj_6983_axis_reconstructed_matched.npz", d_m_rec, "rhj_6983_axis",
                    "reconstructed", tname, ref.matched_barcodes, "", aud_m["n_active_injected"], all_ones_sha)
        _write_pred("preds_Occam_rhj_6983_axis_real_matched.npz", d_m_real, "rhj_6983_axis",
                    "real", tname, ref.matched_barcodes, "", aud_m["n_active_injected"], all_ones_sha)

        x_r, z_r, aud_r = resample_rhj_cache[tname]
        d_r_rec, d_r_real = _predict_deltas_both_baselines_from_z(
            adapter, x_r, z_r, pert_ids, out_cols, batch_size=batch_size, compute_real_baseline=True,
        )
        assert d_r_real is not None
        res_src_label = f"{tname}_within_line_resample_seed{RESAMPLE_SEED}"
        _write_pred("preds_Occam_rhj_6983_axis_reconstructed_resample.npz", d_r_rec, "rhj_6983_axis",
                    "reconstructed", res_src_label, ref.resample_barcodes, "", aud_r["n_active_injected"], all_ones_sha)
        _write_pred("preds_Occam_rhj_6983_axis_real_resample.npz", d_r_real, "rhj_6983_axis",
                    "real", res_src_label, ref.resample_barcodes, "", aud_r["n_active_injected"], all_ones_sha)

        d_mis_rec: Dict[str, np.ndarray] = {}
        d_mis_real: Dict[str, np.ndarray] = {}
        for sname in ref.seen_sources + ref.unseen_sources:
            x_s, z_s, aud_s = primary_rhj_cache[sname]
            ds_rec, ds_real = _predict_deltas_both_baselines_from_z(
                adapter, x_s, z_s, pert_ids, out_cols, batch_size=batch_size, compute_real_baseline=True,
            )
            assert ds_real is not None
            d_mis_rec[sname] = ds_rec
            d_mis_real[sname] = ds_real
            _write_pred(f"preds_Occam_rhj_6983_axis_reconstructed_mis_{sname}.npz", ds_rec, "rhj_6983_axis",
                        "reconstructed", sname, ref.mismatched_barcodes[sname], "", aud_s["n_active_injected"], all_ones_sha)
            _write_pred(f"preds_Occam_rhj_6983_axis_real_mis_{sname}.npz", ds_real, "rhj_6983_axis",
                        "real", sname, ref.mismatched_barcodes[sname], "", aud_s["n_active_injected"], all_ones_sha)
        _log_progress(f"Target {tname}: completed Axis 1 rhj_6983_axis (matched + resample + {len(d_mis_rec)} mis)")

        # 3B. Axis 2: pairwise_common_axis (`P4`)
        p4_pairs_eval: Dict[str, Any] = {}
        raw_t = primary_pools[tname]
        for sname in ref.seen_sources + ref.unseen_sources:
            raw_s = primary_pools[sname]
            la, lb = (tname, sname) if tname < sname else (sname, tname)
            cmask = compute_pairwise_common_rhj_mask(raw_t.raw_ensembl_ids, raw_s.raw_ensembl_ids, rhj_axis_ensembl)
            n_common = int(np.sum(cmask))
            n_union_abs = int(n_rhj_total - n_common)
            cmask_sha = _mask_sha256_uint8(cmask)
            if n_union_abs == 0:
                # Full-coverage pair: strictly reuse rhj_6983_axis arrays without duplicate files (`ADDB_ARRAY_SCHEMA §4.2`)
                p4_pairs_eval[sname] = {
                    "matched": d_m_rec,
                    "mismatched": d_mis_rec[sname],
                    "n_common": n_common,
                    "n_union_absent": 0,
                    "reused_primary_rhj_6983": True,
                    "mask_sha256": cmask_sha,
                }
                continue

            xt_p, zt_p, nt_act, _, reused_t = p4_latent_cache[(tname, la, lb)]
            xs_p, zs_p, ns_act, _, reused_s = p4_latent_cache[(sname, la, lb)]
            if reused_t:
                pm_common = d_m_rec
            else:
                pm_common, _ = _predict_deltas_both_baselines_from_z(
                    adapter, xt_p, zt_p, pert_ids, out_cols, batch_size=batch_size, compute_real_baseline=False,
                )
            if reused_s:
                ps_common = d_mis_rec[sname]
            else:
                ps_common, _ = _predict_deltas_both_baselines_from_z(
                    adapter, xs_p, zs_p, pert_ids, out_cols, batch_size=batch_size, compute_real_baseline=False,
                )

            _write_pred(
                f"preds_Occam_pairwise_common_axis_pair_{sname}_reconstructed_matched.npz",
                pm_common, "pairwise_common_axis", "reconstructed", tname,
                ref.matched_barcodes, sname, nt_act, cmask_sha,
            )
            _write_pred(
                f"preds_Occam_pairwise_common_axis_pair_{sname}_reconstructed_mis_{sname}.npz",
                ps_common, "pairwise_common_axis", "reconstructed", sname,
                ref.mismatched_barcodes[sname], tname, ns_act, cmask_sha,
            )
            p4_pairs_eval[sname] = {
                "matched": pm_common,
                "mismatched": ps_common,
                "n_common": n_common,
                "n_union_absent": n_union_abs,
                "reused_primary_rhj_6983": False,
                "mask_sha256": cmask_sha,
            }
        _log_progress(f"Target {tname}: completed Axis 2 pairwise_common_axis")

        # 3C. Axis 3: full_44728_api
        xf_m, zf_m, nf_m = primary_full_cache[tname]
        df_m, _ = _predict_deltas_both_baselines_from_z(
            adapter, xf_m, zf_m, pert_ids, out_cols, batch_size=batch_size, compute_real_baseline=False,
        )
        _write_pred("preds_Occam_full_44728_api_reconstructed_matched.npz", df_m, "full_44728_api",
                    "reconstructed", tname, ref.matched_barcodes, "", nf_m, "")

        xf_r, zf_r, nf_r = resample_full_cache[tname]
        df_r, _ = _predict_deltas_both_baselines_from_z(
            adapter, xf_r, zf_r, pert_ids, out_cols, batch_size=batch_size, compute_real_baseline=False,
        )
        _write_pred("preds_Occam_full_44728_api_reconstructed_resample.npz", df_r, "full_44728_api",
                    "reconstructed", res_src_label, ref.resample_barcodes, "", nf_r, "")

        df_mis: Dict[str, np.ndarray] = {}
        for sname in ref.seen_sources + ref.unseen_sources:
            xf_s, zf_s, nf_s = primary_full_cache[sname]
            dfs, _ = _predict_deltas_both_baselines_from_z(
                adapter, xf_s, zf_s, pert_ids, out_cols, batch_size=batch_size, compute_real_baseline=False,
            )
            df_mis[sname] = dfs
            _write_pred(f"preds_Occam_full_44728_api_reconstructed_mis_{sname}.npz", dfs, "full_44728_api",
                        "reconstructed", sname, ref.mismatched_barcodes[sname], "", nf_s, "")
        _log_progress(f"Target {tname}: completed Axis 3 full_44728_api ({len(arr_manifest)} files total)")

        _write_json(ldir / "array_manifest_sha256.json", arr_manifest)

        # 4. Compute line metrics + adjudications
        line_entries[tname] = evaluate_addb_line_complete(
            ref=ref,
            occam_matched=d_m_rec,
            occam_resample=d_r_rec,
            occam_mismatched=d_mis_rec,
            p4_common_pairs=p4_pairs_eval,
            occam_baseline_real_matched=d_m_real,
            occam_baseline_real_resample=d_r_real,
            occam_baseline_real_mismatched=d_mis_real,
            occam_full44728_matched=df_m,
            occam_full44728_resample=df_r,
            occam_full44728_mismatched=df_mis,
            pert_resolution_audit=pert_audit,
            output_gene_projection_audit=gene_audit,
            array_manifest_sha256=arr_manifest,
            n_boot=n_boot,
            seed=seed,
        )

    summary = {
        "protocol": "PREREG_I_ADDENDUM_B_v2",
        "model_key": OCCAM_MODEL_KEY,
        "checkpoint_id": OCCAM_CHECKPOINT_ID,
        "ckpt_sha256": str(adapter.contract_summary.get("ckpt_sha256", EXPECTED_BEST_MODEL_SHA256)),
        "device": str(device),
        "torch_num_threads": int(num_threads),
        "seed": int(seed),
        "resample_seed": RESAMPLE_SEED,
        "n_boot": int(n_boot),
        "delta_tost": 0.04,
        "gates_verified": dict(gates_verified),
        "control_barcodes_file": "control_barcodes_addb.json",
        "z_ctrl_dir": "z_ctrl",
        "z_ctrl_manifest_sha256": z_manifest,
        "r16_state_se_reference_sha256": r16_se_shas,
        "group_adjudications": compute_group_adjudications_addb(line_entries),
        "lines": line_entries,
    }
    _write_json(out_dir / "summary_addb.json", summary)
    _log_progress(f"Wrote summary_addb.json cleanly to {out_dir / 'summary_addb.json'}")
    return summary


# ==============================================================================
# Stage-4 Independent From-Arrays Recomputation & Parity Checker (<= 1e-10 / <= 1e-12)
# ==============================================================================


def _compare_json_trees(
    expected: Any,
    got: Any,
    path: str = "",
    tol: float = 1e-10,
) -> Dict[str, Any]:
    """Compare nested JSON dicts/lists recursively; return max_abs_diff & first mismatches."""
    max_diff = 0.0
    worst_path = ""
    n_numeric = 0
    mismatches: List[str] = []

    def _walk(a: Any, b: Any, p: str) -> None:
        nonlocal max_diff, worst_path, n_numeric
        if isinstance(a, Mapping) and isinstance(b, Mapping):
            ak = set(a.keys())
            bk = set(b.keys())
            if ak != bk:
                mismatches.append(f"{p}: key diff missing={sorted(ak - bk)}, extra={sorted(bk - ak)}")
            for k in sorted(ak & bk):
                _walk(a[k], b[k], f"{p}.{k}" if p else str(k))
            return
        if isinstance(a, list) and isinstance(b, list):
            if len(a) != len(b):
                mismatches.append(f"{p}: len {len(a)} != {len(b)}")
                return
            for i, (av, bv) in enumerate(zip(a, b)):
                _walk(av, bv, f"{p}[{i}]")
            return
        if isinstance(a, (bool, np.bool_)) or isinstance(b, (bool, np.bool_)):
            if bool(a) != bool(b):
                mismatches.append(f"{p}: bool {a!r} != {b!r}")
            return
        if isinstance(a, (int, float, np.number)) and isinstance(b, (int, float, np.number)):
            fa, fb = float(a), float(b)
            n_numeric += 1
            if not (math.isfinite(fa) and math.isfinite(fb)):
                if math.isfinite(fa) != math.isfinite(fb):
                    mismatches.append(f"{p}: finiteness {fa} vs {fb}")
                return
            d = abs(fa - fb)
            if d > max_diff:
                max_diff = d
                worst_path = p
            if d > tol:
                mismatches.append(f"{p}: |{fa} - {fb}| = {d:.3e} > {tol:.1e}")
            return
        if a is None and b is None:
            return
        if a != b:
            mismatches.append(f"{p}: {a!r} != {b!r}")

    _walk(expected, got, path)
    return {
        "n_numeric_leaves": int(n_numeric),
        "max_abs_diff": float(max_diff),
        "worst_path": worst_path,
        "mismatches": mismatches,
        "passed": bool(len(mismatches) == 0 and max_diff <= tol),
    }


def _infer_pert_audit_from_npz(
    eval_perts: Sequence[str],
    resolved_keys: Sequence[str],
    eligible_rhj2381: Optional[Sequence[str]],
) -> Dict[str, Any]:
    elig_set = set(eligible_rhj2381) if eligible_rhj2381 is not None else None
    t2_targets: List[str] = []
    bridged: List[Dict[str, str]] = []
    n_t1 = 0
    for ep, rk in zip(eval_perts, resolved_keys):
        if rk != ep:
            bridged.append({"eval": str(ep), "resolved": str(rk)})
        if elig_set is not None:
            is_t1 = str(rk) in elig_set
        else:
            is_t1 = str(ep) != "RBM14"
        if is_t1:
            n_t1 += 1
        else:
            t2_targets.append(str(ep))
    return {
        "n_tier1": int(n_t1),
        "n_tier2": int(len(t2_targets)),
        "tier2_targets": t2_targets,
        "bridged_pairs": bridged,
        "unresolved_targets": [],
    }


def _reconstruct_ref_from_exported_dir(
    target_line: str,
    ldir: Path,
    se_audit_dir: Path,
    main_r16_dir: Path,
    main_sens_dir: Path,
    repo_root: Path,
) -> CommittedR16ReferenceBundle:
    """Load observed matrices from `<results_dir>/arrays/<line>/observed.npz` and paired STATE_SE arrays."""
    with np.load(ldir / "observed.npz", allow_pickle=False) as oz:
        required_obs = {
            "delta_obs", "delta_half_a", "delta_half_b", "half_valid", "delta_transfer",
            "source_grand_mean", "fct_delta_pred", "e1_constmed_delta_pred", "target_indices",
            "target_names", "gene_names", "mu_target", "mu_sources",
            "matched_ctrl_barcodes", "resample_ctrl_barcodes",
            "n_ctrl_total", "n_ctrl_sampled", "ctrl_sampling_fraction",
        }
        missing_obs = required_obs - set(oz.files)
        if missing_obs:
            raise KeyError(f"{ldir / 'observed.npz'} missing required keys: {sorted(missing_obs)}")
        eval_genes = [str(x) for x in oz["gene_names"].tolist()]
        eval_perts = [str(x) for x in oz["target_names"].tolist()]
        tg = np.asarray(oz["target_indices"], dtype=np.int64)
        d_obs = np.asarray(oz["delta_obs"], dtype=np.float64)
        ha = np.asarray(oz["delta_half_a"], dtype=np.float64)
        hb = np.asarray(oz["delta_half_b"], dtype=np.float64)
        hv = np.asarray(oz["half_valid"], dtype=bool)
        tm = np.asarray(oz["delta_transfer"], dtype=np.float64)
        gm = np.asarray(oz["source_grand_mean"], dtype=np.float64)
        fct_d = np.asarray(oz["fct_delta_pred"], dtype=np.float64)
        e1_d = np.asarray(oz["e1_constmed_delta_pred"], dtype=np.float64)
        mu_t = np.asarray(oz["mu_target"], dtype=np.float64)
        mu_s = np.asarray(oz["mu_sources"], dtype=np.float64)
        m_bcs = [str(x) for x in oz["matched_ctrl_barcodes"].tolist()]
        r_bcs = [str(x) for x in oz["resample_ctrl_barcodes"].tolist()]
        n_ctrl_tot = int(oz["n_ctrl_total"])
        samp_frac = float(oz["ctrl_sampling_fraction"])

    # Load matching committed STATE_SE deltas without recomputing R16 source lines unnecessarily
    se_ref = load_committed_r16_reference_and_se(
        target_line,
        audit_results_dir=se_audit_dir,
        main_r16_dir=main_r16_dir,
        main_sens_dir=main_sens_dir,
        repo_root=repo_root,
        load_observed_extras=False,
    )
    if se_ref.eval_perts != eval_perts or se_ref.eval_genes != eval_genes:
        raise ValueError(f"Recompute target/gene mismatch vs STATE_SE reference on {target_line}")
    if float(np.max(np.abs(se_ref.delta_obs - d_obs))) > 1e-12:
        raise ValueError(f"observed.npz delta_obs diverges from STATE_SE reference on {target_line}")

    extra_fct: Optional[Dict[str, Any]] = None
    sfp = ldir / "preds_SENS_FCT_full_r9_6_incl_vcc2025_h1.npz"
    if target_line == "gse295214_h1" and sfp.is_file():
        with np.load(sfp, allow_pickle=False) as sz:
            extra_fct = {
                "source_lines_fct": list(se_ref.extra_fct_full_r9_6.get("source_lines_fct", [])) if se_ref.extra_fct_full_r9_6 else [
                    "nadig_hepg2", "nadig_jurkat", "nourreddine_kolf21j", "replogle_k562_essential", "replogle_rpe1", "vcc2025_h1"
                ],
                "delta_pred": np.asarray(sz["delta_pred"], dtype=np.float64),
                "delta_transfer": np.asarray(sz["delta_transfer"], dtype=np.float64),
                "source_grand_mean": np.asarray(sz["source_grand_mean"], dtype=np.float64),
            }

    return CommittedR16ReferenceBundle(
        target_line=target_line,
        role=ADDB_LINE_ROLES[target_line],
        eval_perts=eval_perts,
        eval_genes=eval_genes,
        target_gene_idx=tg,
        delta_obs=d_obs,
        delta_half_a=ha,
        delta_half_b=hb,
        half_valid=hv,
        delta_transfer=tm,
        source_grand_mean=gm,
        delta_fct=fct_d,
        seen_sources=list(ADDB_SEEN_SOURCES),
        unseen_sources=list(ADDB_UNSEEN_SOURCES[target_line]),
        matched_barcodes=m_bcs,
        resample_barcodes=r_bcs,
        mismatched_barcodes=se_ref.mismatched_barcodes,
        se_delta_matched=se_ref.se_delta_matched,
        se_delta_resample=se_ref.se_delta_resample,
        se_delta_mismatched=se_ref.se_delta_mismatched,
        mu_target=mu_t,
        mu_sources=mu_s,
        e1_constmed=e1_d,
        n_ctrl_total=n_ctrl_tot,
        ctrl_sampling_fraction=samp_frac,
        r16_reference_sha256=se_ref.r16_reference_sha256,
        extra_fct_full_r9_6=extra_fct,
    )


def recompute_addb_from_exported_arrays(
    results_dir: Path,
    se_audit_dir: Path = DEFAULT_LOCKED_SE_AUDIT_DIR,
    main_r16_dir: Path = DEFAULT_LOCKED_MAIN_R16_DIR,
    main_sens_dir: Path = DEFAULT_LOCKED_MAIN_SENS_DIR,
    repo_root: Path = REPO_ROOT,
    tol: float = 1e-10,
    perturbations_csv: Path = DEFAULT_PERTURBATIONS_CSV,
) -> Dict[str, Any]:
    """
    Stage-4 contract (`ADDB_ARRAY_SCHEMA.md §8`): recompute EVERY metric, bootstrap CI,
    `H-I1_Occam`, `H-I2_Occam`, `H-I3a`, `H-I3b`, `pairwise_common_axis`, `baseline_real`,
    `full_44728_api`, `tier1_trained_only_sensitivity`, and `extra_fct_sensitivity` directly
    from `<results_dir>/arrays/` `.npz` files (`allow_pickle=False`), verifying manifests and
    asserting agreement with `<results_dir>/summary_addb.json` within `tol` (`<= 1e-10`).
    """
    rdir = Path(results_dir)
    arrays_root = rdir / "arrays"
    summary_path = rdir / "summary_addb.json"
    cb_path = rdir / "control_barcodes_addb.json"
    z_manifest_path = rdir / "z_ctrl" / "z_ctrl_manifest_sha256.json"
    orig_summary = json.loads(summary_path.read_text(encoding="utf-8"))
    cb_json = json.loads(cb_path.read_text(encoding="utf-8"))
    eligible_2381 = load_rhj_eligible_perturbations(perturbations_csv) if Path(perturbations_csv).is_file() else None

    if z_manifest_path.is_file():
        z_man = json.loads(z_manifest_path.read_text(encoding="utf-8"))
        for fname, want_sha in z_man.items():
            zfp = rdir / "z_ctrl" / fname
            if sha256_file(zfp) != want_sha:
                raise ValueError(f"z_ctrl SHA-256 mismatch: {zfp}")
            with np.load(zfp, allow_pickle=False) as zz:
                if sha256_array_bytes(np.ascontiguousarray(zz["Z_ctrl"], dtype=np.float32)) != str(zz["z_ctrl_sha256"]):
                    raise ValueError(f"Z_ctrl internal z_ctrl_sha256 mismatch: {zfp}")

    n_boot = int(orig_summary.get("n_boot", 4000))
    seed = int(orig_summary.get("seed", DEFAULT_R16_SEED))

    def _load_pred(fp: Path, ref: CommittedR16ReferenceBundle) -> Tuple[np.ndarray, int, str, List[str]]:
        req_keys = {
            "delta_pred", "target_names", "gene_names", "resolved_pert_keys", "resolved_pert_ids",
            "model_key", "checkpoint_id", "input_axis_mode", "baseline_mode", "control_source",
            "control_barcodes", "pair_counterpart", "n_active_rhj_genes", "axis_mask_sha256", "delta_sha256",
        }
        with np.load(fp, allow_pickle=False) as z:
            miss = req_keys - set(z.files)
            if miss:
                raise KeyError(f"{fp} missing keys: {sorted(miss)}")
            if [str(x) for x in z["target_names"].tolist()] != ref.eval_perts:
                raise ValueError(f"target_names mismatch in {fp}")
            if [str(x) for x in z["gene_names"].tolist()] != ref.eval_genes:
                raise ValueError(f"gene_names mismatch in {fp}")
            dp = np.ascontiguousarray(np.asarray(z["delta_pred"], dtype=np.float64))
            if sha256_array_bytes(dp) != str(z["delta_sha256"]):
                raise ValueError(f"delta_sha256 mismatch in {fp}")
            rkeys = [str(x) for x in z["resolved_pert_keys"].tolist()]
            return dp, int(z["n_active_rhj_genes"]), str(z["axis_mask_sha256"]), rkeys

    recomputed_lines: Dict[str, Any] = {}
    for tname in sorted(orig_summary["lines"].keys()):
        ldir = arrays_root / tname
        man = json.loads((ldir / "array_manifest_sha256.json").read_text(encoding="utf-8"))
        for fname, want_sha in man.items():
            afp = ldir / fname
            if sha256_file(afp) != want_sha:
                raise ValueError(f"Array manifest SHA-256 mismatch: {afp}")

        ref = _reconstruct_ref_from_exported_dir(
            tname, ldir, se_audit_dir, main_r16_dir, main_sens_dir, repo_root,
        )
        # Verify locked barcode consistency against control_barcodes_addb.json
        if cb_json["pools"][tname]["barcodes"] != ref.matched_barcodes:
            raise ValueError(f"Matched barcodes mismatch vs control_barcodes_addb.json on {tname}")

        dm_rec, _, _, rkeys = _load_pred(ldir / "preds_Occam_rhj_6983_axis_reconstructed_matched.npz", ref)
        dr_rec, _, _, _ = _load_pred(ldir / "preds_Occam_rhj_6983_axis_reconstructed_resample.npz", ref)
        dmis_rec = {
            s: _load_pred(ldir / f"preds_Occam_rhj_6983_axis_reconstructed_mis_{s}.npz", ref)[0]
            for s in ref.seen_sources + ref.unseen_sources
        }

        dm_real, _, _, _ = _load_pred(ldir / "preds_Occam_rhj_6983_axis_real_matched.npz", ref)
        dr_real, _, _, _ = _load_pred(ldir / "preds_Occam_rhj_6983_axis_real_resample.npz", ref)
        dmis_real = {
            s: _load_pred(ldir / f"preds_Occam_rhj_6983_axis_real_mis_{s}.npz", ref)[0]
            for s in ref.seen_sources + ref.unseen_sources
        }

        df_m, _, _, _ = _load_pred(ldir / "preds_Occam_full_44728_api_reconstructed_matched.npz", ref)
        df_r, _, _, _ = _load_pred(ldir / "preds_Occam_full_44728_api_reconstructed_resample.npz", ref)
        df_mis = {
            s: _load_pred(ldir / f"preds_Occam_full_44728_api_reconstructed_mis_{s}.npz", ref)[0]
            for s in ref.seen_sources + ref.unseen_sources
        }

        # Infer total RHJ axis width & all-ones SHA directly from exported full-coverage z_ctrl latent
        ref_z_fp = rdir / "z_ctrl" / "z_ctrl_rhj_6983_axis__nadig_hepg2__primary.npz"
        with np.load(ref_z_fp, allow_pickle=False) as rzz:
            n_rhj_axis_len = int(rzz["n_active_rhj_genes"])
            all_ones_sha = str(rzz["axis_mask_sha256"])

        p4_pairs: Dict[str, Any] = {}
        for s in ref.seen_sources + ref.unseen_sources:
            p_m_fp = ldir / f"preds_Occam_pairwise_common_axis_pair_{s}_reconstructed_matched.npz"
            p_s_fp = ldir / f"preds_Occam_pairwise_common_axis_pair_{s}_reconstructed_mis_{s}.npz"
            if p_m_fp.is_file() and p_s_fp.is_file():
                pm, n_c_m, msha, _ = _load_pred(p_m_fp, ref)
                ps, n_c_s, _, _ = _load_pred(p_s_fp, ref)
                n_common = min(n_c_m, n_c_s)
                p4_pairs[s] = {
                    "matched": pm,
                    "mismatched": ps,
                    "n_common": int(n_common),
                    "n_union_absent": int(n_rhj_axis_len - n_common),
                    "reused_primary_rhj_6983": False,
                    "mask_sha256": msha,
                }
            else:
                p4_pairs[s] = {
                    "matched": dm_rec,
                    "mismatched": dmis_rec[s],
                    "n_common": int(n_rhj_axis_len),
                    "n_union_absent": 0,
                    "reused_primary_rhj_6983": True,
                    "mask_sha256": all_ones_sha,
                }

        pert_aud = _infer_pert_audit_from_npz(ref.eval_perts, rkeys, eligible_2381)
        gene_aud = {
            "n_target_genes": int(len(ref.eval_genes)),
            "n_mapped": int(len(ref.eval_genes)),
            "n_missing": 0,
            "rhj_hgnc_collisions": 0,
        }
        recomputed_lines[tname] = evaluate_addb_line_complete(
            ref=ref,
            occam_matched=dm_rec,
            occam_resample=dr_rec,
            occam_mismatched=dmis_rec,
            p4_common_pairs=p4_pairs,
            occam_baseline_real_matched=dm_real,
            occam_baseline_real_resample=dr_real,
            occam_baseline_real_mismatched=dmis_real,
            occam_full44728_matched=df_m,
            occam_full44728_resample=df_r,
            occam_full44728_mismatched=df_mis,
            pert_resolution_audit=pert_aud,
            output_gene_projection_audit=gene_aud,
            array_manifest_sha256=man,
            n_boot=n_boot,
            seed=seed,
        )

    recomputed_summary = {
        "group_adjudications": compute_group_adjudications_addb(recomputed_lines),
        "lines": recomputed_lines,
    }
    orig_comparable = {
        "group_adjudications": orig_summary["group_adjudications"],
        "lines": orig_summary["lines"],
    }
    cmp_res = _compare_json_trees(_json_safe(orig_comparable), _json_safe(recomputed_summary), tol=tol)
    if not cmp_res["passed"]:
        raise AssertionError(
            f"Stage-4 recompute parity failed (max_abs={cmp_res['max_abs_diff']}, mismatches={cmp_res['mismatches'][:5]})"
        )
    return cmp_res


# ==============================================================================
# Synthetic Hermetic E2E Fixture & Driver (`--phase synthetic_e2e`)
# ==============================================================================


class _MockOccamCleanRoomAdapter:
    """Deterministic lightweight CPU mock reproducing Occam's per-cell encoder/FiLM structure."""

    def __init__(self, encoder_gene_names: Sequence[str], ckpt_pert2id: Mapping[str, int], target_sum: float = 2000.0) -> None:
        self.gene_names = list(encoder_gene_names)
        self.pert2id = {str(k): int(v) for k, v in ckpt_pert2id.items()}
        self.target_sum = float(target_sum)
        self.contract_summary = {"ckpt_sha256": EXPECTED_BEST_MODEL_SHA256, "mock": True}

    def encode_controls(self, X_ctrl_enc: np.ndarray, batch_size: int = 256) -> np.ndarray:
        # Linear projection (512, G_enc) -> (512, 512) float32 sensitive to gene missingness & pool identity
        x = np.asarray(X_ctrl_enc, dtype=np.float32)
        n = x.shape[0]
        rng = np.random.default_rng(777)
        proj = rng.normal(0.0, 0.05, size=(x.shape[1], 512)).astype(np.float32)
        return (x @ proj).astype(np.float32)

    def predict_reconstructed_baseline(self, z_ctrl: np.ndarray, batch_size: int = 256, out_cols: Optional[np.ndarray] = None) -> np.ndarray:
        cols = np.arange(len(self.gene_names)) if out_cols is None else np.asarray(out_cols, dtype=np.int64)
        z = np.asarray(z_ctrl, dtype=np.float32)
        dec_w = np.sin(np.arange(512 * cols.size, dtype=np.float32).reshape(512, cols.size) * 0.01) * 0.02
        return np.log1p(np.exp(np.clip(z @ dec_w, -4.0, 4.0))).astype(np.float32)

    def predict_per_cell(self, z_ctrl: np.ndarray, pert_indices: Sequence[Sequence[int]], batch_size: int = 256, out_cols: Optional[np.ndarray] = None) -> np.ndarray:
        base = self.predict_reconstructed_baseline(z_ctrl, batch_size=batch_size, out_cols=out_cols)
        cols = np.arange(len(self.gene_names)) if out_cols is None else np.asarray(out_cols, dtype=np.int64)
        pid = int(pert_indices[0][0])
        z_mean = np.mean(np.asarray(z_ctrl, dtype=np.float32), axis=1, keepdims=True)
        pert_vec = np.cos((np.arange(cols.size, dtype=np.float32) + 1.0) * (pid + 1.0) * 0.07)[None, :]
        # Combine control-sensitive FiLM shift + perturbation signal
        shift = 0.08 * pert_vec * (1.0 + 0.35 * np.tanh(z_mean))
        return np.maximum(base + shift, 1e-6).astype(np.float32)


class _MockHGNCBridge:
    def canonicalize(self, sym: str) -> str:
        s = str(sym).strip()
        if s == "AARS":
            return "AARS1"
        return s


class _MockTwoTierPertResolver:
    def __init__(self, ckpt_pert2id: Mapping[str, int]) -> None:
        self.ckpt_pert2id = dict(ckpt_pert2id)

    def resolve_many(self, target_symbols: Sequence[str]) -> Tuple[List[int], List[Dict[str, Any]]]:
        ids: List[int] = []
        rows: List[Dict[str, Any]] = []
        for sym in target_symbols:
            s = str(sym)
            if s == "AARS1" and "AARS" in self.ckpt_pert2id:
                rk, tier = "AARS", "tier1_rhj2381"
            elif s == "RBM14":
                rk, tier = "RBM14", "tier2_direct"
            else:
                rk, tier = s, "tier1_rhj2381"
            ids.append(self.ckpt_pert2id[rk])
            rows.append({
                "eval_target": s,
                "resolved_key": rk,
                "pert_id": int(self.ckpt_pert2id[rk]),
                "tier": tier,
                "bridged": bool(rk != s),
            })
        return ids, rows


def run_synthetic_e2e_addb(output_dir: Path, n_boot: int = 200, seed: int = DEFAULT_R16_SEED) -> Dict[str, Any]:
    """
    Execute complete pure-CPU synthetic E2E Addendum B pipeline without touching real data,
    generating exact `ADDB_ARRAY_SCHEMA.md` directory tree + standalone mock STATE_SE dirs
    so `--recompute-from-arrays` succeeds out-of-the-box in a separate process.
    """
    import audit_addb_controls_and_metrics as mod_cm
    import audit_addb_occam_adapter as mod_adp

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    se_mock_dir = out_dir / "_mock_se_audit"
    main_r16_mock = out_dir / "_mock_main_r16"
    main_sens_mock = out_dir / "_mock_main_sens"
    (se_mock_dir / "group3").mkdir(parents=True, exist_ok=True)
    main_r16_mock.mkdir(parents=True, exist_ok=True)
    (main_sens_mock / "arrays" / "feng_ipsc").mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(int(seed))
    n_enc = 160
    n_rhj = 80
    n_eval_g = 40
    n_eval_p = 54  # >= 50 so valid_frac == 1.0 cleanly in unit/E2E TOST
    n_ctrl = 64

    enc_genes = [f"ENSG{100000 + i:06d}" for i in range(n_enc)]
    rhj_ens = enc_genes[:n_rhj]
    rhj_syms = [f"GENE{i:03d}" for i in range(n_rhj)]
    eval_genes = sorted(rhj_syms[:n_eval_g])

    if DEFAULT_PERTURBATIONS_CSV.is_file():
        elig = [k for k in load_rhj_eligible_perturbations(DEFAULT_PERTURBATIONS_CSV) if k not in ("AARS", "AARS1", "RBM14")]
        filler_perts = elig[: max(0, n_eval_p - 2)]
        filler_extra = elig[max(0, n_eval_p - 2)]
    else:
        filler_perts = [f"PERT{i:03d}" for i in range(n_eval_p - 2)]
        filler_extra = "PERT999"

    base_perts = sorted(["AARS1", "RBM14"] + filler_perts)
    ckpt_p2i = {p: idx for idx, p in enumerate(["AARS"] + base_perts)}

    absent_idx = {
        "a549_ko": {50, 51, 52},
        "cd4t_rest": {50, 53, 54, 55, 56, 57, 58, 59},
        "gse295214_h1": {60},
        "dld1": set(),
        "feng_ipsc": set(),
        "nadig_hepg2": set(),
        "nadig_jurkat": set(),
    }

    primary_pools: Dict[str, ExtractedRawControlCells] = {}
    resample_pools: Dict[str, ExtractedRawControlCells] = {}
    locked_pools: Dict[str, Dict[str, Any]] = {}

    for p_idx, pname in enumerate(ADDB_ALL_CONTROL_POOLS):
        keep_ens = [e for idx, e in enumerate(enc_genes) if idx not in absent_idx[pname]]
        p_bcs = [f"{pname}_pri_{k:04d}" for k in range(n_ctrl)]
        r_bcs = [f"{pname}_res_{k:04d}" for k in range(n_ctrl)] if pname in ADDB_ALL_TARGET_LINES else None
        raw_p = rng.poisson(lam=2.0 + 0.3 * p_idx, size=(n_ctrl, len(keep_ens))).astype(np.float32)
        primary_pools[pname] = ExtractedRawControlCells(
            pool_name=pname, raw_counts=raw_p, raw_ensembl_ids=keep_ens,
            barcodes=p_bcs, row_sum_umis=raw_p.sum(axis=1, dtype=np.float64), n_ctrl_total=5000,
        )
        if r_bcs is not None:
            raw_r = rng.poisson(lam=2.0 + 0.3 * p_idx, size=(n_ctrl, len(keep_ens))).astype(np.float32)
            resample_pools[pname] = ExtractedRawControlCells(
                pool_name=pname, raw_counts=raw_r, raw_ensembl_ids=keep_ens,
                barcodes=r_bcs, row_sum_umis=raw_r.sum(axis=1, dtype=np.float64), n_ctrl_total=5000,
            )
        locked_pools[pname] = {"primary_barcodes": p_bcs, "resample_barcodes": r_bcs, "n_ctrl_total": 5000}

    locked_barcodes = {
        "pools": locked_pools,
        "provenance": {"primary_seed": int(seed), "resample_seed": RESAMPLE_SEED, "r16_sources_sha256": {"synthetic": "0" * 64}},
    }

    orig_proj = sys.modules[__name__].resolve_output_projection
    orig_occ_num = mod_adp.OCCAM_NUM_ENCODER_GENES

    ref_bundles: Dict[str, CommittedR16ReferenceBundle] = {}
    for tname in ADDB_ALL_TARGET_LINES:
        ep = list(base_perts) if tname == "a549_ko" else [p for p in base_perts if p != "RBM14"] + [filler_extra]
        ep = sorted(ep)
        ckpt_p2i.setdefault(filler_extra, len(ckpt_p2i))
        p_len = len(ep)
        tg = np.full(p_len, -1, dtype=np.int64)
        tg[:min(5, p_len)] = np.arange(min(5, p_len), dtype=np.int64)
        signal = rng.normal(0.0, 0.35, size=(p_len, n_eval_g))
        d_obs = signal + rng.normal(0.0, 0.12, size=(p_len, n_eval_g))
        ha = signal + rng.normal(0.0, 0.15, size=(p_len, n_eval_g))
        hb = signal + rng.normal(0.0, 0.15, size=(p_len, n_eval_g))
        hv = np.ones(p_len, dtype=bool)
        tm = rng.normal(0.0, 0.15, size=(p_len, n_eval_g))
        gm = rng.normal(0.0, 0.08, size=n_eval_g)
        fct_d = tm + 0.35 * signal
        e1_d = tm + 0.10 * signal
        mu_t = rng.uniform(0.2, 3.0, size=n_eval_g)
        mu_s = rng.uniform(0.2, 3.0, size=(5 if tname == "gse295214_h1" else 6, n_eval_g))
        se_m = tm + 0.20 * signal
        se_r = se_m + rng.normal(0.0, 0.01, size=(p_len, n_eval_g))
        seen = list(ADDB_SEEN_SOURCES)
        unseen = list(ADDB_UNSEEN_SOURCES[tname])
        se_mis = {s: se_m + rng.normal(0.0, 0.05, size=(p_len, n_eval_g)) for s in seen + unseen}
        mis_bcs = {s: locked_pools[s]["primary_barcodes"] for s in seen + unseen}

        extra_fct = None
        if tname == "gse295214_h1":
            extra_fct = {
                "source_lines_fct": ["nadig_hepg2", "nadig_jurkat", "nourreddine_kolf21j", "replogle_k562_essential", "replogle_rpe1", "vcc2025_h1"],
                "delta_pred": fct_d * 1.05,
                "delta_transfer": tm * 1.02,
                "source_grand_mean": gm * 1.01,
            }

        # Export standalone mock R16 STATE_SE bundles for independent offline Stage-4 CLI verification
        ref_shas: Dict[str, str] = {}
        if tname == "feng_ipsc":
            fdir = main_sens_mock / "arrays" / "feng_ipsc"
            obs_fp = write_npz_strict_addb(
                fdir / "observed.npz",
                {
                    "delta_obs": d_obs, "delta_half_a": ha, "delta_half_b": hb, "half_valid": hv,
                    "delta_transfer": tm, "source_grand_mean": gm, "fct_delta_pred": fct_d,
                    "e1_constmed_delta_pred": e1_d, "target_indices": tg,
                    "target_names": ep, "gene_names": eval_genes,
                    "mu_target": mu_t, "mu_sources": mu_s,
                    "matched_ctrl_barcodes": locked_pools[tname]["primary_barcodes"],
                    "resample_ctrl_barcodes": locked_pools[tname]["resample_barcodes"],
                    "n_ctrl_total": 5000, "n_ctrl_sampled": n_ctrl,
                    "ctrl_sampling_fraction": float(n_ctrl) / 5000.0,
                },
            )
            ref_shas["observed.npz"] = sha256_file(obs_fp)
            mfp = write_npz_strict_addb(fdir / "preds_STATE_SE_fixed_index_zero_matched.npz", {"delta_pred": se_m})
            rfp = write_npz_strict_addb(fdir / "preds_STATE_SE_fixed_index_zero_resample.npz", {"delta_pred": se_r})
            ref_shas["matched"] = sha256_file(mfp)
            ref_shas["resample"] = sha256_file(rfp)
            for s in seen + unseen:
                sfp = write_npz_strict_addb(
                    fdir / f"preds_STATE_SE_fixed_index_zero_mis_{s}.npz",
                    {"delta_pred": se_mis[s], "control_barcodes": mis_bcs[s]},
                )
                ref_shas[f"mis_{s}"] = sha256_file(sfp)
        else:
            se_base = (se_mock_dir / "group3") if tname == "gse295214_h1" else se_mock_dir
            se_payload: Dict[str, Any] = {
                "eval_genes": eval_genes, "eval_perts": ep, "target_gene_idx": tg,
                "delta_obs": d_obs, "delta_half_a": ha, "delta_half_b": hb, "half_valid": hv,
                "delta_transfer": tm, "source_grand_mean": gm, "delta_fct": fct_d,
                "ctrl_mean_libnorm": mu_t,
                "delta_pred_matched": se_m, "delta_pred_resample": se_r,
                "matched_barcodes": locked_pools[tname]["primary_barcodes"],
                "resample_barcodes": locked_pools[tname]["resample_barcodes"],
            }
            for s in seen + unseen:
                se_payload[f"delta_pred_mismatched__{s}"] = se_mis[s]
                se_payload[f"mismatched_barcodes__{s}"] = mis_bcs[s]
            se_fp = write_npz_strict_addb(se_base / f"STATE_SE_ensemble4__{tname}.npz", se_payload)
            ref_shas[se_fp.name] = sha256_file(se_fp)

        ref_bundles[tname] = CommittedR16ReferenceBundle(
            target_line=tname, role=ADDB_LINE_ROLES[tname],
            eval_perts=ep, eval_genes=list(eval_genes), target_gene_idx=tg,
            delta_obs=d_obs, delta_half_a=ha, delta_half_b=hb, half_valid=hv,
            delta_transfer=tm, source_grand_mean=gm, delta_fct=fct_d,
            seen_sources=seen, unseen_sources=unseen,
            matched_barcodes=locked_pools[tname]["primary_barcodes"],
            resample_barcodes=locked_pools[tname]["resample_barcodes"],
            mismatched_barcodes=mis_bcs,
            se_delta_matched=se_m, se_delta_resample=se_r, se_delta_mismatched=se_mis,
            mu_target=mu_t, mu_sources=mu_s, e1_constmed=e1_d,
            n_ctrl_total=5000, ctrl_sampling_fraction=float(n_ctrl) / 5000.0,
            r16_reference_sha256=ref_shas, extra_fct_full_r9_6=extra_fct,
        )

    def _mock_resolve_out(eg: Sequence[str], rhj_e: Sequence[str], rhj_s: Sequence[str], enc_g: Sequence[str], hgnc_bridge: Any = None) -> Dict[str, Any]:
        sym_to_pos = {s: i for i, s in enumerate(rhj_s)}
        rhj_cols = [sym_to_pos[g] for g in eg]
        return {
            "eval_genes": list(eg), "n_eval_genes": len(eg), "n_mapped": len(eg),
            "missing_genes": [],
            "rhj_cols": np.asarray(rhj_cols, dtype=np.int64),
            "enc_cols": np.asarray(rhj_cols, dtype=np.int64),
        }

    try:
        globals()["resolve_output_projection"] = _mock_resolve_out
        mod_adp.OCCAM_NUM_ENCODER_GENES = n_enc

        mock_adapter = _MockOccamCleanRoomAdapter(enc_genes, ckpt_p2i)
        mock_hb = _MockHGNCBridge()
        mock_resolver = _MockTwoTierPertResolver(ckpt_p2i)

        summary = _run_export_pipeline_core(
            adapter=mock_adapter,
            rhj_axis_ensembl=rhj_ens,
            rhj_axis_symbols=rhj_syms,
            hgnc_bridge=mock_hb,
            pert_resolver=mock_resolver,
            ref_bundles=ref_bundles,
            locked_barcodes=locked_barcodes,
            primary_pools=primary_pools,
            resample_pools=resample_pools,
            output_dir=out_dir,
            gates_verified={"mode": "synthetic_e2e_mock"},
            device=LOCKED_DEVICE,
            num_threads=LOCKED_NUM_THREADS,
            batch_size=64,
            n_boot=n_boot,
            seed=seed,
        )

        cmp_res = recompute_addb_from_exported_arrays(
            results_dir=out_dir,
            se_audit_dir=se_mock_dir,
            main_r16_dir=main_r16_mock,
            main_sens_dir=main_sens_mock,
            tol=1e-10,
            perturbations_csv=out_dir / "__nonexistent_pert__.csv",
        )
        return {"summary": summary, "recompute_parity": cmp_res}
    finally:
        globals()["resolve_output_projection"] = orig_proj
        mod_adp.OCCAM_NUM_ENCODER_GENES = orig_occ_num


# ==============================================================================
# Real Addendum B Execution Driver (`--phase real_addb`, Gate-Protected)
# ==============================================================================


def run_real_addb(
    output_dir: Path,
    addb_root: Optional[Path] = None,
    device: str = LOCKED_DEVICE,
    num_threads: int = LOCKED_NUM_THREADS,
    batch_size: int = 256,
    n_boot: int = 4000,
    seed: int = DEFAULT_R16_SEED,
    nadig_raw_dir: Path = DEFAULT_NADIG_RAW_DIR,
    gse295214_raw_dir: Path = DEFAULT_GSE295214_RAW_DIR,
    se_audit_dir: Path = DEFAULT_LOCKED_SE_AUDIT_DIR,
    main_r16_dir: Path = DEFAULT_LOCKED_MAIN_R16_DIR,
    main_sens_dir: Path = DEFAULT_LOCKED_MAIN_SENS_DIR,
    gate5b_summary_path: Optional[Path] = None,
) -> Dict[str, Any]:
    """Run full gate-protected Addendum B evaluation on locked CPU float32 (4 threads)."""
    gates = assert_freeze_gates_addb(addb_root=addb_root)
    if str(device).strip().lower() != LOCKED_DEVICE:
        raise ValueError(f"Gate 5A.3 locked device is {LOCKED_DEVICE!r}; refusing --device {device!r}")
    torch.set_num_threads(int(num_threads))

    _log_progress("Step 1/5: Verifying Stage 5A asset checksums and Stage 5B HepG2 reproduction gate")
    verify_addb_checksums()
    adapter = OccamCleanRoomAdapter(device=LOCKED_DEVICE, num_threads=int(num_threads), verify_sha256=True)

    if gate5b_summary_path is not None and Path(gate5b_summary_path).is_file():
        g5b = json.loads(Path(gate5b_summary_path).read_text(encoding="utf-8"))
        if not g5b.get("overall_cpu_gates_pass", False):
            raise RuntimeError(f"Gate 5B summary did not pass: {gate5b_summary_path}")
        _log_progress(f"Verified existing Gate 5B CPU summary pass: {gate5b_summary_path}")
    else:
        _log_progress("Running Stage 5B HepG2 benchmark gate on CPU prior to any target scoring...")
        g5b_dir = Path(output_dir) / "_gate5b_precheck"
        g5b = run_stage_5b_cpu(adapter, output_dir=g5b_dir, batch_size=batch_size)
        if not g5b.get("passed", False):
            raise RuntimeError("Stage 5B benchmark gate FAILED — halting prior to any target scoring!")
        _log_progress("Stage 5B benchmark gate PASSED cleanly")

    _log_progress("Step 2/5: Loading locked Round 16 control barcodes and committed STATE_SE references")
    locked_barcodes = load_locked_control_barcodes_addb(
        ADDB_ALL_TARGET_LINES,
        audit_results_dir=se_audit_dir,
        main_r16_dir=main_r16_dir,
        main_sens_dir=main_sens_dir,
        repo_root=REPO_ROOT,
    )
    primary_pools, resample_pools = extract_all_required_raw_pools_addb(
        locked_barcodes,
        nadig_raw_dir=nadig_raw_dir,
        gse295214_raw_dir=gse295214_raw_dir,
    )
    ref_bundles = {
        ln: load_committed_r16_reference_and_se(
            ln,
            audit_results_dir=se_audit_dir,
            main_r16_dir=main_r16_dir,
            main_sens_dir=main_sens_dir,
            repo_root=REPO_ROOT,
        )
        for ln in ADDB_ALL_TARGET_LINES
    }

    rhj_ens, rhj_syms, _ = load_rhj_gene_axis(DEFAULT_GENE_AXIS_CSV, verify_sha256=True)
    hgnc_bridge = load_hgnc_bridge()
    pert_resolver = TwoTierPertResolver(adapter.pert2id, DEFAULT_PERTURBATIONS_CSV, hgnc_bridge=hgnc_bridge)

    _log_progress("Step 3/5: Running clean-room CPU float32 inference across all axes and arms")
    summary = _run_export_pipeline_core(
        adapter=adapter,
        rhj_axis_ensembl=rhj_ens,
        rhj_axis_symbols=rhj_syms,
        hgnc_bridge=hgnc_bridge,
        pert_resolver=pert_resolver,
        ref_bundles=ref_bundles,
        locked_barcodes=locked_barcodes,
        primary_pools=primary_pools,
        resample_pools=resample_pools,
        output_dir=Path(output_dir),
        gates_verified=gates,
        device=LOCKED_DEVICE,
        num_threads=int(num_threads),
        batch_size=batch_size,
        n_boot=n_boot,
        seed=seed,
    )

    _log_progress("Step 5/5: Running Stage-4 independent from-arrays recomputation check")
    recompute_res = recompute_addb_from_exported_arrays(
        results_dir=Path(output_dir),
        se_audit_dir=se_audit_dir,
        main_r16_dir=main_r16_dir,
        main_sens_dir=main_sens_dir,
        repo_root=REPO_ROOT,
        tol=1e-10,
    )
    _log_progress(
        f"Stage-4 recompute parity PASSED: n_numeric={recompute_res['n_numeric_leaves']}, "
        f"max_abs_diff={recompute_res['max_abs_diff']:.3e}"
    )
    return {"summary": summary, "recompute_parity": recompute_res}


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Clean-room Round 16 Addendum B (Occam) audit runner")
    p.add_argument(
        "--phase",
        choices=("real_addb", "synthetic_e2e", "recompute_only"),
        default="recompute_only",
        help="Execution mode: real_addb (gate-locked), synthetic_e2e, or recompute_only.",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=ADDB_ROOT / "audit" / "results" / "addb_run",
        help="Output results directory for real_addb or synthetic_e2e.",
    )
    p.add_argument(
        "--recompute-from-arrays",
        type=Path,
        default=None,
        help="Recompute summary_addb.json directly from exported .npz arrays in this results dir.",
    )
    p.add_argument("--device", type=str, default=LOCKED_DEVICE)
    p.add_argument("--num-threads", type=int, default=LOCKED_NUM_THREADS)
    p.add_argument("--batch-size", type=int, default=256)
    p.add_argument("--n-boot", type=int, default=4000)
    p.add_argument("--seed", type=int, default=DEFAULT_R16_SEED)
    p.add_argument("--gate5b-summary", type=Path, default=None)
    p.add_argument("--se-audit-dir", type=Path, default=DEFAULT_LOCKED_SE_AUDIT_DIR)
    p.add_argument("--main-r16-dir", type=Path, default=DEFAULT_LOCKED_MAIN_R16_DIR)
    p.add_argument("--main-sens-dir", type=Path, default=DEFAULT_LOCKED_MAIN_SENS_DIR)
    p.add_argument("--nadig-raw-dir", type=Path, default=DEFAULT_NADIG_RAW_DIR)
    p.add_argument("--gse295214-raw-dir", type=Path, default=DEFAULT_GSE295214_RAW_DIR)
    p.add_argument("--tol", type=float, default=1e-10)
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    if args.recompute_from_arrays is not None:
        res = recompute_addb_from_exported_arrays(
            results_dir=args.recompute_from_arrays,
            se_audit_dir=args.se_audit_dir,
            main_r16_dir=args.main_r16_dir,
            main_sens_dir=args.main_sens_dir,
            repo_root=REPO_ROOT,
            tol=float(args.tol),
        )
        print(json.dumps({"mode": "recompute_from_arrays", **res}, indent=2, sort_keys=True))
        return 0

    if args.phase == "synthetic_e2e":
        res = run_synthetic_e2e_addb(output_dir=args.output_dir, n_boot=args.n_boot, seed=args.seed)
        print(json.dumps({
            "mode": "synthetic_e2e",
            "output_dir": str(args.output_dir),
            "recompute_parity": res["recompute_parity"],
            "lines": sorted(res["summary"]["lines"].keys()),
        }, indent=2, sort_keys=True))
        return 0

    if args.phase == "real_addb":
        res = run_real_addb(
            output_dir=args.output_dir,
            device=args.device,
            num_threads=args.num_threads,
            batch_size=args.batch_size,
            n_boot=args.n_boot,
            seed=args.seed,
            nadig_raw_dir=args.nadig_raw_dir,
            gse295214_raw_dir=args.gse295214_raw_dir,
            se_audit_dir=args.se_audit_dir,
            main_r16_dir=args.main_r16_dir,
            main_sens_dir=args.main_sens_dir,
            gate5b_summary_path=args.gate5b_summary,
        )
        print(json.dumps({
            "mode": "real_addb",
            "output_dir": str(args.output_dir),
            "recompute_parity": res["recompute_parity"],
        }, indent=2, sort_keys=True))
        return 0

    raise ValueError("Pass --phase real_addb, --phase synthetic_e2e, or --recompute-from-arrays <dir>")


if __name__ == "__main__":
    raise SystemExit(main())
