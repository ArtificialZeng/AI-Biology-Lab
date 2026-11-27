#!/usr/bin/env python3
"""
Round 16 Addendum B (Occam) Gate 5A Item 3 (MPS-vs-CPU float32 parity) and
conditional Gate 5B MPS rerun (`PREREG_I_ADDENDUM_B_DRAFT_v2.md §5.1.3 & §5.2`).

Locked rules:
  1. Compare identical synthetic inputs on CPU float32 (`torch.set_num_threads(4)`)
     vs. MPS float32 (`allow_mps=True`).
  2. Lock `mps` iff:
       max|z_ctrl_mps - z_ctrl_cpu| <= 1e-5
       AND max|delta_mps - delta_cpu| <= 1e-4
     (checked across both `baseline="reconstructed"` primary and `baseline="real"` sensitivity);
     otherwise lock `cpu` for both parties.
  3. If `mps` is locked, rerun Gate 5B (`nadig_hepg2`, exact Case A criteria) on `mps`.
     If `cpu` is locked, the already-passed CPU Gate 5B (`gates_cpu_summary.json`) stands.
  4. ZERO Group 2, Group 3, or `feng_ipsc` Occam inference.
"""

from __future__ import annotations

import argparse
import gc
import json
import os
import subprocess
import sys
import tempfile
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict

import numpy as np
import torch

ADDB_SRC_DIR = Path(__file__).resolve().parent
ADDB_ROOT = ADDB_SRC_DIR.parent
if str(ADDB_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(ADDB_SRC_DIR))

from occam_adapter_addb import (  # noqa: E402
    DEFAULT_OCCAM_CHECKSUMS_FILE,
    DEFAULT_OCCAM_CKPT_PATH,
    DEFAULT_OCCAM_GENE_AXIS_CSV,
    DEFAULT_OCCAM_PERTURBATIONS_CSV,
    DEFAULT_OCCAM_SPLITS_MANIFEST_JSON,
    DEFAULT_OCCAM_SRC_DIR,
    OccamInferenceAdapter,
    configure_occam_runtime,
)
from raw_controls_addb import DEFAULT_LOCAL_NADIG_DIR  # noqa: E402
from run_gates_cpu_addb import run_gate5b_hepg2_cpu  # noqa: E402

Z_CTRL_MAX_ABS_TOL = 1e-5
DELTA_MAX_ABS_TOL = 1e-4


def _now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _run_synthetic_forward_in_subprocess(
    *,
    device: str,
    num_threads: int,
    n_cells: int,
    n_perts: int,
    seed: int,
    ckpt_path: Path,
    occam_src_dir: Path,
    gene_axis_csv: Path,
    perturbations_csv: Path,
    out_npz: Path,
    pythonhashseed: str = "1",
) -> Dict[str, Any]:
    """
    Execute identical synthetic forward pass in a fresh isolated subprocess so CPU
    and MPS have zero shared global PyTorch deterministic/allocator state.
    """
    code = f"""
import json, sys, time, numpy as np, torch
from pathlib import Path
sys.path.insert(0, {str(ADDB_SRC_DIR)!r})
from occam_adapter_addb import OccamInferenceAdapter
from run_gates_cpu_addb import _make_synthetic_controls

t0 = time.time()
dev = {device!r}
if dev == "mps":
    torch.use_deterministic_algorithms(False)
    if not torch.backends.mps.is_available():
        raise RuntimeError("torch.backends.mps.is_available() is False")
    if hasattr(torch, "mps") and hasattr(torch.mps, "manual_seed"):
        torch.mps.manual_seed({int(seed)})

adp, _ = OccamInferenceAdapter.load_verified(
    ckpt_path=Path({str(ckpt_path)!r}),
    occam_src_dir=Path({str(occam_src_dir)!r}),
    gene_axis_csv=Path({str(gene_axis_csv)!r}),
    perturbations_csv=Path({str(perturbations_csv)!r}),
    device=dev,
    num_threads={int(num_threads)},
    allow_mps=(dev == "mps"),
    batch_size=32,
    verify_file_hashes=False,
)
x = _make_synthetic_controls({int(n_cells)}, len(adp.encoder_gene_names), adp.rhj_enc_cols, seed={int(seed)}, scale=1.0)
sample_perts = sorted(adp.rhj_eligible_keys)[:{int(n_perts)}]
sample_pids = [adp.pert2id[p] for p in sample_perts]

res1 = adp.predict_deltas_from_controls(x, sample_pids, adp.rhj_enc_cols, batch_size=32)
res2 = adp.predict_deltas_from_controls(x, sample_pids, adp.rhj_enc_cols, batch_size=32)
if dev == "mps":
    torch.mps.synchronize()

repeat_z = float(np.max(np.abs(res1["z_ctrl"] - res2["z_ctrl"])))
repeat_d_recon = float(np.max(np.abs(res1["delta_reconstructed"] - res2["delta_reconstructed"])))
repeat_d_real = float(np.max(np.abs(res1["delta_real"] - res2["delta_real"])))

np.savez(
    {str(out_npz)!r},
    z_ctrl=res1["z_ctrl"].astype(np.float32),
    delta_reconstructed=res1["delta_reconstructed"].astype(np.float64),
    delta_real=res1["delta_real"].astype(np.float64),
    pred_means=res1["pred_means"].astype(np.float64),
    recon_ctrl_mean=res1["recon_ctrl_mean"].astype(np.float64),
)
print(json.dumps({{
    "device": dev,
    "torch_num_threads": int(torch.get_num_threads()),
    "sample_perts": sample_perts,
    "elapsed_s": round(time.time() - t0, 3),
    "repeat_max_abs_z_ctrl": repeat_z,
    "repeat_max_abs_delta_reconstructed": repeat_d_recon,
    "repeat_max_abs_delta_real": repeat_d_real,
    "z_all_finite": bool(np.all(np.isfinite(res1["z_ctrl"]))),
    "delta_all_finite": bool(np.all(np.isfinite(res1["delta_reconstructed"])) and np.all(np.isfinite(res1["delta_real"]))),
}}))
"""
    env = os.environ.copy()
    env["PYTHONHASHSEED"] = str(pythonhashseed)
    env["OMP_NUM_THREADS"] = str(int(num_threads))
    env["MKL_NUM_THREADS"] = str(int(num_threads))
    env["OPENBLAS_NUM_THREADS"] = str(int(num_threads))
    proc = subprocess.run(
        [sys.executable, "-c", code],
        env=env,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"Subprocess failed (device={device}, rc={proc.returncode}):\n"
            f"STDOUT:\n{proc.stdout[-1500:]}\nSTDERR:\n{proc.stderr[-2000:]}"
        )
    last_line = [ln for ln in proc.stdout.strip().splitlines() if ln.strip()][-1]
    return json.loads(last_line)


def run_gate5a_item3_and_conditional_5b(
    *,
    num_threads: int = 4,
    n_cells: int = 32,
    n_perts: int = 8,
    seed: int = 20261601,
    batch_size: int = 512,
    ckpt_path: Path = DEFAULT_OCCAM_CKPT_PATH,
    occam_src_dir: Path = DEFAULT_OCCAM_SRC_DIR,
    gene_axis_csv: Path = DEFAULT_OCCAM_GENE_AXIS_CSV,
    perturbations_csv: Path = DEFAULT_OCCAM_PERTURBATIONS_CSV,
    splits_manifest_json: Path = DEFAULT_OCCAM_SPLITS_MANIFEST_JSON,
    checksums_file: Path = DEFAULT_OCCAM_CHECKSUMS_FILE,
    nadig_hepg2_h5ad: Path = DEFAULT_LOCAL_NADIG_DIR / "GSE264667_hepg2_raw_singlecell_01.h5ad",
    out_json: Path = ADDB_ROOT / "results" / "gate5a3_mps_summary.json",
) -> Dict[str, Any]:
    t_start = time.time()
    print(
        f"[{_now_iso()}] Starting Addendum B Gate 5A Item 3 (CPU threads={num_threads} vs MPS float32) | "
        f"python={sys.version.split()[0]} torch={torch.__version__} mps_avail={torch.backends.mps.is_available()}",
        flush=True,
    )

    with tempfile.TemporaryDirectory(prefix="occam_g5a3_") as td:
        td_path = Path(td)
        cpu_npz = td_path / "cpu_synth.npz"
        mps_npz = td_path / "mps_synth.npz"
        mps_h2_npz = td_path / "mps_synth_h99991.npz"

        cpu_meta = _run_synthetic_forward_in_subprocess(
            device="cpu",
            num_threads=num_threads,
            n_cells=n_cells,
            n_perts=n_perts,
            seed=seed,
            ckpt_path=ckpt_path,
            occam_src_dir=occam_src_dir,
            gene_axis_csv=gene_axis_csv,
            perturbations_csv=perturbations_csv,
            out_npz=cpu_npz,
            pythonhashseed="1",
        )
        print(f"[{_now_iso()}] CPU synthetic forward done in {cpu_meta['elapsed_s']}s.", flush=True)

        mps_error = None
        mps_meta = None
        mps_cross_proc_max_abs = None
        max_abs_z_ctrl = float("inf")
        max_abs_delta_recon = float("inf")
        max_abs_delta_real = float("inf")
        max_abs_delta_overall = float("inf")
        max_abs_pred_means = float("inf")
        max_abs_recon_ctrl = float("inf")

        if not torch.backends.mps.is_available():
            mps_error = "torch.backends.mps.is_available() returned False"
        else:
            try:
                mps_meta = _run_synthetic_forward_in_subprocess(
                    device="mps",
                    num_threads=num_threads,
                    n_cells=n_cells,
                    n_perts=n_perts,
                    seed=seed,
                    ckpt_path=ckpt_path,
                    occam_src_dir=occam_src_dir,
                    gene_axis_csv=gene_axis_csv,
                    perturbations_csv=perturbations_csv,
                    out_npz=mps_npz,
                    pythonhashseed="1",
                )
                print(f"[{_now_iso()}] MPS synthetic forward done in {mps_meta['elapsed_s']}s.", flush=True)

                # Cross-process PYTHONHASHSEED check on MPS (1 vs 99991)
                _run_synthetic_forward_in_subprocess(
                    device="mps",
                    num_threads=num_threads,
                    n_cells=n_cells,
                    n_perts=n_perts,
                    seed=seed,
                    ckpt_path=ckpt_path,
                    occam_src_dir=occam_src_dir,
                    gene_axis_csv=gene_axis_csv,
                    perturbations_csv=perturbations_csv,
                    out_npz=mps_h2_npz,
                    pythonhashseed="99991",
                )

                with np.load(cpu_npz) as c_arr, np.load(mps_npz) as m_arr, np.load(mps_h2_npz) as m2_arr:
                    max_abs_z_ctrl = float(np.max(np.abs(m_arr["z_ctrl"].astype(np.float64) - c_arr["z_ctrl"].astype(np.float64))))
                    max_abs_delta_recon = float(np.max(np.abs(m_arr["delta_reconstructed"] - c_arr["delta_reconstructed"])))
                    max_abs_delta_real = float(np.max(np.abs(m_arr["delta_real"] - c_arr["delta_real"])))
                    max_abs_delta_overall = float(max(max_abs_delta_recon, max_abs_delta_real))
                    max_abs_pred_means = float(np.max(np.abs(m_arr["pred_means"] - c_arr["pred_means"])))
                    max_abs_recon_ctrl = float(np.max(np.abs(m_arr["recon_ctrl_mean"] - c_arr["recon_ctrl_mean"])))
                    mps_cross_proc_max_abs = float(
                        max(
                            np.max(np.abs(m_arr["z_ctrl"] - m2_arr["z_ctrl"])),
                            np.max(np.abs(m_arr["delta_reconstructed"] - m2_arr["delta_reconstructed"])),
                            np.max(np.abs(m_arr["delta_real"] - m2_arr["delta_real"])),
                        )
                    )
            except Exception as exc:
                mps_error = f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}"
                print(f"[{_now_iso()}] MPS synthetic forward raised exception; locking CPU.\n{mps_error}", flush=True)

    pass_z = bool(np.isfinite(max_abs_z_ctrl) and max_abs_z_ctrl <= Z_CTRL_MAX_ABS_TOL)
    pass_delta_recon = bool(np.isfinite(max_abs_delta_recon) and max_abs_delta_recon <= DELTA_MAX_ABS_TOL)
    pass_delta_real = bool(np.isfinite(max_abs_delta_real) and max_abs_delta_real <= DELTA_MAX_ABS_TOL)
    pass_mps_det = bool(
        mps_meta is not None
        and mps_meta["z_all_finite"]
        and mps_meta["delta_all_finite"]
        and mps_meta["repeat_max_abs_z_ctrl"] <= 1e-6
        and mps_meta["repeat_max_abs_delta_reconstructed"] <= 1e-6
        and mps_meta["repeat_max_abs_delta_real"] <= 1e-6
        and (mps_cross_proc_max_abs is not None and mps_cross_proc_max_abs <= 1e-6)
    )

    mps_parity_passed = bool(mps_error is None and pass_z and pass_delta_recon and pass_delta_real and pass_mps_det)
    locked_device = "mps" if mps_parity_passed else "cpu"

    print(
        f"[{_now_iso()}] Gate 5A Item 3 Result: max|z_ctrl_mps - z_ctrl_cpu|={max_abs_z_ctrl:.3e} (<=1e-5: {pass_z}), "
        f"max|delta_recon_mps - delta_recon_cpu|={max_abs_delta_recon:.3e} (<=1e-4: {pass_delta_recon}), "
        f"max|delta_real_mps - delta_real_cpu|={max_abs_delta_real:.3e} (<=1e-4: {pass_delta_real}), "
        f"mps_det_ok={pass_mps_det} -> LOCKED DEVICE = {locked_device.upper()}",
        flush=True,
    )

    gate5b_mps_result: Any = None
    gate5b_mps_vs_cpu_diff: Any = None
    if locked_device == "mps":
        print(f"[{_now_iso()}] MPS locked -> rerunning Gate 5B (nadig_hepg2) on MPS...", flush=True)
        torch.use_deterministic_algorithms(False)
        configure_occam_runtime(device="mps", num_threads=num_threads, seed=seed, allow_mps=True)
        if hasattr(torch, "mps") and hasattr(torch.mps, "manual_seed"):
            torch.mps.manual_seed(int(seed))
        adapter_mps, _ = OccamInferenceAdapter.load_verified(
            ckpt_path=ckpt_path,
            occam_src_dir=occam_src_dir,
            gene_axis_csv=gene_axis_csv,
            perturbations_csv=perturbations_csv,
            checksums_file=checksums_file,
            device="mps",
            num_threads=num_threads,
            allow_mps=True,
            batch_size=batch_size,
            verify_file_hashes=True,
        )
        gate5b_mps_result = run_gate5b_hepg2_cpu(
            adapter_mps,
            h5ad_path=nadig_hepg2_h5ad,
            gene_axis_csv=gene_axis_csv,
            perturbations_csv=perturbations_csv,
            splits_manifest_json=splits_manifest_json,
            batch_size=batch_size,
        )
        torch.mps.synchronize()
        del adapter_mps
        gc.collect()
        torch.mps.empty_cache()

        # Cross-check against existing CPU Gate 5B summary if present
        cpu_sum_path = ADDB_ROOT / "results" / "gates_cpu_summary.json"
        if cpu_sum_path.is_file():
            cpu_g5b = json.loads(cpu_sum_path.read_text(encoding="utf-8")).get("gate5b", {}).get("metrics", {})
            mps_m = gate5b_mps_result["metrics"]
            if cpu_g5b:
                gate5b_mps_vs_cpu_diff = {
                    "abs_diff_pearson_delta_mps_vs_cpu": abs(
                        float(mps_m["adapter_pearson_delta_mean"]) - float(cpu_g5b["adapter_pearson_delta_mean"])
                    ),
                    "abs_diff_mae_mps_vs_cpu": abs(
                        float(mps_m["adapter_mae_mean"]) - float(cpu_g5b["adapter_mae_mean"])
                    ),
                    "abs_diff_mae_delta_mps_vs_cpu": abs(
                        float(mps_m["adapter_mae_delta_mean"]) - float(cpu_g5b["adapter_mae_delta_mean"])
                    ),
                }
    else:
        gate5b_mps_result = "NOT_RUN_CPU_LOCKED_EXISTING_CPU_GATE5B_STANDS"

    overall_passed = bool(
        (locked_device == "cpu")
        or (locked_device == "mps" and isinstance(gate5b_mps_result, dict) and gate5b_mps_result.get("passed") is True)
    )

    summary = {
        "timestamp_utc": _now_iso(),
        "elapsed_total_s": round(time.time() - t_start, 2),
        "python_version": sys.version.split()[0],
        "torch_version": str(torch.__version__),
        "cpu_num_threads": int(num_threads),
        "mps_available": bool(torch.backends.mps.is_available()),
        "synthetic_spec": {
            "n_cells": int(n_cells),
            "n_perts": int(n_perts),
            "seed": int(seed),
            "sample_perts": cpu_meta["sample_perts"],
        },
        "gate5a_item3": {
            "tol_z_ctrl": Z_CTRL_MAX_ABS_TOL,
            "tol_delta": DELTA_MAX_ABS_TOL,
            "max_abs_z_ctrl_mps_minus_cpu": max_abs_z_ctrl,
            "max_abs_delta_reconstructed_mps_minus_cpu": max_abs_delta_recon,
            "max_abs_delta_real_mps_minus_cpu": max_abs_delta_real,
            "max_abs_delta_overall_mps_minus_cpu": max_abs_delta_overall,
            "max_abs_pred_means_mps_minus_cpu": max_abs_pred_means,
            "max_abs_recon_ctrl_mps_minus_cpu": max_abs_recon_ctrl,
            "pass_z_ctrl_le_1e-5": pass_z,
            "pass_delta_reconstructed_le_1e-4": pass_delta_recon,
            "pass_delta_real_le_1e-4": pass_delta_real,
            "pass_mps_determinism_le_1e-6": pass_mps_det,
            "mps_cross_process_pythonhashseed_max_abs": mps_cross_proc_max_abs,
            "cpu_run_meta": cpu_meta,
            "mps_run_meta": mps_meta,
            "mps_error": mps_error,
            "mps_parity_passed": mps_parity_passed,
            "locked_device": locked_device,
        },
        "gate5b_mps_rerun": gate5b_mps_result,
        "gate5b_mps_vs_cpu_diff": gate5b_mps_vs_cpu_diff,
        "overall_passed": overall_passed,
    }

    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)
    return summary


def main() -> int:
    ap = argparse.ArgumentParser(description="Run Round 16 Addendum B Gate 5A Item 3 & conditional Gate 5B MPS rerun.")
    ap.add_argument("--num-threads", type=int, default=4)
    ap.add_argument("--n-cells", type=int, default=32)
    ap.add_argument("--n-perts", type=int, default=8)
    ap.add_argument("--seed", type=int, default=20261601)
    ap.add_argument("--batch-size", type=int, default=512)
    ap.add_argument(
        "--out-json",
        type=Path,
        default=ADDB_ROOT / "results" / "gate5a3_mps_summary.json",
    )
    args = ap.parse_args()
    res = run_gate5a_item3_and_conditional_5b(
        num_threads=args.num_threads,
        n_cells=args.n_cells,
        n_perts=args.n_perts,
        seed=args.seed,
        batch_size=args.batch_size,
        out_json=args.out_json,
    )
    return 0 if res["overall_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
