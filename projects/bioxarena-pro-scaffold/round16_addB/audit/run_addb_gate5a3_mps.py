#!/usr/bin/env python3
"""
Clean-Room Addendum B Stage 5A item 3 (`MPS float32` vs `CPU float32` Parity & Lock)
and conditional Stage 5B MPS Rerun (`round16_addB/audit/run_addb_gate5a3_mps.py`, `A-T3`).

Strictly enforces `PREREG_I_ADDENDUM_B_DRAFT_v2.md §5.1.3`:
  - Sets `torch.set_num_threads(4)`.
  - On identical synthetic inputs (`n = 32` cells, `seed = 20261601`), compares `CPU float32` vs `MPS float32`.
  - Locks `mps` iff `max|z_ctrl_mps - z_ctrl_cpu| <= 1e-5` AND `max|delta_mps - delta_cpu| <= 1e-4`;
    otherwise locks `cpu`.
  - If `mps` is locked, reruns Stage 5B (`run_stage_5b_cpu`) on `adapter_mps`.
  - NEVER runs inference on Group 2, Group 3, or `feng_ipsc` (`R16_ADDB_GO.txt` required).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import traceback
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import torch

AUDIT_DIR = Path(__file__).resolve().parent
if str(AUDIT_DIR) not in sys.path:
    sys.path.insert(0, str(AUDIT_DIR))

from audit_addb_occam_adapter import (  # noqa: E402
    DEFAULT_CKPT_PATH,
    DEFAULT_GENE_AXIS_CSV,
    DEFAULT_OCCAM_SRC_ROOT,
    OCCAM_NUM_ENCODER_GENES,
    OccamCleanRoomAdapter,
    load_rhj_gene_axis,
    sha256_file,
)
from run_addb_gates_cpu import run_stage_5b_cpu  # noqa: E402

Z_CTRL_PARITY_ATOL = 1e-5
DELTA_PARITY_ATOL = 1e-4


def _ts() -> str:
    t = os.times()
    return f"[{time.strftime('%H:%M:%S')} cpu={t.user + t.system:.1f}s]"


def decide_addb_locked_device(
    max_abs_z_ctrl: float,
    max_abs_delta: float,
    mps_available: bool = True,
    mps_error: Optional[str] = None,
) -> Tuple[str, bool]:
    """
    Apply locked §5.1.3 rule:
      lock `'mps'` iff `mps_available`, no MPS runtime error,
      `max|z_ctrl_mps - z_ctrl_cpu| <= 1e-5`, and `max|delta_mps - delta_cpu| <= 1e-4`;
      else lock `'cpu'`.
    """
    ok = bool(
        mps_available
        and mps_error is None
        and np.isfinite(max_abs_z_ctrl)
        and np.isfinite(max_abs_delta)
        and float(max_abs_z_ctrl) <= Z_CTRL_PARITY_ATOL
        and float(max_abs_delta) <= DELTA_PARITY_ATOL
    )
    return ("mps" if ok else "cpu"), ok


def build_synthetic_parity_inputs(
    adapter_cpu: OccamCleanRoomAdapter,
    seed: int = 20261601,
    n_cells: int = 32,
    n_perts: int = 8,
) -> Tuple[np.ndarray, List[int], np.ndarray]:
    """Construct identical pure-synthetic `float32` control counts and perturbation IDs."""
    axis_ens, _ = load_rhj_gene_axis(DEFAULT_GENE_AXIS_CSV)
    enc_cols = np.asarray([adapter_cpu.gene_to_index[g] for g in axis_ens], dtype=np.int64)

    rng = np.random.default_rng(int(seed))
    x_ctrl = np.zeros((int(n_cells), OCCAM_NUM_ENCODER_GENES), dtype=np.float32)
    x_ctrl[:, enc_cols] = rng.uniform(0.0, 3.0, size=(int(n_cells), enc_cols.size)).astype(np.float32)

    all_pids = sorted(int(v) for v in adapter_cpu.pert2id.values() if int(v) >= 0)
    pick_idx = np.linspace(0, len(all_pids) - 1, num=int(n_perts), dtype=np.int64)
    pert_ids = [all_pids[int(i)] for i in pick_idx]
    return x_ctrl, pert_ids, enc_cols


def run_gate_5a3_and_conditional_5b(
    num_threads: int = 4,
    batch_size: int = 256,
    out_json: Path = AUDIT_DIR / "results" / "gate5a3_mps_summary.json",
    cpu_5b_reference_json: Optional[Path] = AUDIT_DIR / "results" / "gates_cpu_summary.json",
) -> Dict[str, Any]:
    """Execute Stage 5A item 3 CPU vs MPS parity and rerun Gate 5B on MPS iff MPS locks."""
    torch.set_num_threads(int(num_threads))
    print(f"{_ts()} [5A.3] Initializing Occam CPU adapter (threads={num_threads})...", flush=True)
    adapter_cpu = OccamCleanRoomAdapter(
        ckpt_path=DEFAULT_CKPT_PATH,
        device="cpu",
        num_threads=int(num_threads),
        occam_src_root=DEFAULT_OCCAM_SRC_ROOT,
        gene_axis_csv=DEFAULT_GENE_AXIS_CSV,
        verify_sha256=True,
    )

    x_ctrl, pert_ids, enc_cols = build_synthetic_parity_inputs(adapter_cpu)
    print(f"{_ts()} [5A.3] Running CPU float32 synthetic parity pass ({x_ctrl.shape[0]} cells x {len(pert_ids)} perts)...", flush=True)
    z_cpu = adapter_cpu.encode_controls(x_ctrl, batch_size=batch_size)
    delta_cpu, _, _ = adapter_cpu.predict_deltas(
        x_ctrl, pert_ids, output_col_idx=enc_cols, baseline="reconstructed", batch_size=batch_size, z_ctrl=z_cpu
    )

    mps_avail = bool(torch.backends.mps.is_available())
    mps_error: Optional[str] = None
    adapter_mps: Optional[OccamCleanRoomAdapter] = None
    max_abs_z = float("nan")
    max_abs_delta = float("nan")
    mps_repeat_max_abs = float("nan")

    if mps_avail:
        try:
            print(f"{_ts()} [5A.3] Initializing Occam MPS float32 adapter...", flush=True)
            adapter_mps = OccamCleanRoomAdapter(
                ckpt_path=DEFAULT_CKPT_PATH,
                device="mps",
                num_threads=int(num_threads),
                occam_src_root=DEFAULT_OCCAM_SRC_ROOT,
                gene_axis_csv=DEFAULT_GENE_AXIS_CSV,
                verify_sha256=False,
            )
            print(f"{_ts()} [5A.3] Running MPS float32 synthetic parity pass...", flush=True)
            z_mps = adapter_mps.encode_controls(x_ctrl, batch_size=batch_size)
            z_mps_2 = adapter_mps.encode_controls(x_ctrl, batch_size=batch_size)
            mps_repeat_max_abs = float(np.max(np.abs(z_mps - z_mps_2)))
            delta_mps, _, _ = adapter_mps.predict_deltas(
                x_ctrl, pert_ids, output_col_idx=enc_cols, baseline="reconstructed", batch_size=batch_size, z_ctrl=z_mps
            )
            max_abs_z = float(np.max(np.abs(z_mps - z_cpu)))
            max_abs_delta = float(np.max(np.abs(delta_mps - delta_cpu)))
        except Exception:
            mps_error = traceback.format_exc()
    else:
        mps_error = "torch.backends.mps.is_available() == False"

    locked_device, mps_locked = decide_addb_locked_device(
        max_abs_z_ctrl=max_abs_z,
        max_abs_delta=max_abs_delta,
        mps_available=mps_avail,
        mps_error=mps_error,
    )
    print(
        f"{_ts()} [5A.3] Parity: max|z_ctrl_mps - z_ctrl_cpu|={max_abs_z:.3e} (tol<={Z_CTRL_PARITY_ATOL:.1e}), "
        f"max|delta_mps - delta_cpu|={max_abs_delta:.3e} (tol<={DELTA_PARITY_ATOL:.1e}) -> locked_device={locked_device}",
        flush=True,
    )

    stage_5b_mps: Optional[Dict[str, Any]] = None
    vs_cpu_5b: Optional[Dict[str, Any]] = None
    if mps_locked and adapter_mps is not None:
        print(f"{_ts()} [5B-MPS] Rerunning Stage 5B HepG2 reproduction gate on MPS...", flush=True)
        stage_5b_mps = run_stage_5b_cpu(adapter=adapter_mps, batch_size=batch_size)
        if cpu_5b_reference_json is not None and Path(cpu_5b_reference_json).is_file():
            cpu_ref = json.loads(Path(cpu_5b_reference_json).read_text(encoding="utf-8"))["stage_5b"]
            vs_cpu_5b = {
                "pearson_delta_abs_diff": abs(
                    float(stage_5b_mps["adapter_mean_pearson_delta"]) - float(cpu_ref["adapter_mean_pearson_delta"])
                ),
                "mae_abs_diff": abs(
                    float(stage_5b_mps["adapter_mean_mae"]) - float(cpu_ref["adapter_mean_mae"])
                ),
                "experimental_pearson_delta_abs_diff": abs(
                    float(stage_5b_mps["local_experimental_pearson_delta"])
                    - float(cpu_ref["local_experimental_pearson_delta"])
                ),
            }

    payload: Dict[str, Any] = {
        "torch_version": torch.__version__,
        "torch_num_threads": int(num_threads),
        "mps_available": mps_avail,
        "mps_error": mps_error,
        "synthetic_shape": {
            "n_cells": int(x_ctrl.shape[0]),
            "n_perts": len(pert_ids),
            "n_output_genes": int(enc_cols.size),
        },
        "parity": {
            "max_abs_z_ctrl": max_abs_z,
            "z_ctrl_tol": Z_CTRL_PARITY_ATOL,
            "z_ctrl_pass": bool(np.isfinite(max_abs_z) and max_abs_z <= Z_CTRL_PARITY_ATOL),
            "max_abs_delta": max_abs_delta,
            "delta_tol": DELTA_PARITY_ATOL,
            "delta_pass": bool(np.isfinite(max_abs_delta) and max_abs_delta <= DELTA_PARITY_ATOL),
            "mps_repeat_z_max_abs": mps_repeat_max_abs,
        },
        "locked_device": locked_device,
        "mps_locked": mps_locked,
        "stage_5b_mps": stage_5b_mps,
        "stage_5b_mps_vs_cpu_diff": vs_cpu_5b,
        "gate_5b_mps_pass": None if stage_5b_mps is None else bool(stage_5b_mps["gate_5b_pass"]),
    }

    out_path = Path(out_json)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    out_sha = sha256_file(out_path)
    print(f"{_ts()} [5A.3] Saved {out_path} sha256={out_sha}", flush=True)
    return payload


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Run Addendum B Stage 5A.3 (CPU vs MPS parity) & conditional 5B MPS.")
    ap.add_argument("--num-threads", type=int, default=4)
    ap.add_argument("--batch-size", type=int, default=256)
    ap.add_argument("--output-json", type=Path, default=AUDIT_DIR / "results" / "gate5a3_mps_summary.json")
    args = ap.parse_args(argv)

    res = run_gate_5a3_and_conditional_5b(
        num_threads=args.num_threads,
        batch_size=args.batch_size,
        out_json=args.output_json,
    )
    if res["mps_locked"] and not res["gate_5b_mps_pass"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
