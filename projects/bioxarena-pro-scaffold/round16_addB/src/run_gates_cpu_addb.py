#!/usr/bin/env python3
"""
CPU-only Gate 5A (Items 1–2) & Gate 5B (`nadig_hepg2` 5-split author reproduction)
runner for Round 16 Addendum B (`Occam`, `PREREG_I_ADDENDUM_B_DRAFT_v2.md §5.1–§5.2`).

Strict guardrails:
  - CPU ONLY (`device="cpu"`, `torch.set_num_threads(6)` so A-T3 keeps its cores).
  - ZERO MPS/GPU calls (Gate 5A Item 3 waits until Controller confirms A-T3 finished R16).
  - ZERO Group 2, Group 3, or `feng_ipsc` inference.
  - Halts immediately (non-zero exit) on ANY Gate 5A or Gate 5B failure; zero threshold changes.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import torch

ADDB_SRC_DIR = Path(__file__).resolve().parent
ADDB_ROOT = ADDB_SRC_DIR.parent
if str(ADDB_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(ADDB_SRC_DIR))

from gate5b_panel_addb import (  # noqa: E402
    GATE5B_AUTHOR_EXP_PEARSON_DELTA,
    SPLITHALF_CONTROL_LABEL,
    build_gate5b_hepg2_panel,
)
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
from preprocess_addb import load_rhj_gene_axis, load_rhj_perturbations_table  # noqa: E402
from raw_controls_addb import DEFAULT_LOCAL_NADIG_DIR  # noqa: E402

GATE5B_AUTHOR_MODEL_PEARSON_DELTA = 0.394178
GATE5B_AUTHOR_MODEL_MAE = 0.030157
GATE5B_AUTHOR_MODEL_MAE_DELTA = 0.030284


def _now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _cpu_time_s() -> float:
    t = os.times()
    return float(t.user + t.system)


def _make_synthetic_controls(
    n_cells: int,
    n_enc: int,
    rhj_enc_cols: np.ndarray,
    seed: int,
    scale: float = 1.0,
) -> np.ndarray:
    rng = np.random.default_rng(int(seed))
    x = np.zeros((int(n_cells), int(n_enc)), dtype=np.float32)
    raw = rng.poisson(lam=float(scale) * 1.5, size=(int(n_cells), len(rhj_enc_cols))).astype(np.float32)
    rs = np.maximum(raw.sum(axis=1, keepdims=True, dtype=np.float64), 1.0)
    x[:, rhj_enc_cols] = np.log1p((raw / rs) * 2000.0).astype(np.float32)
    return np.ascontiguousarray(x, dtype=np.float32)


def run_gate5a_item2_synthetic(
    adapter: OccamInferenceAdapter,
    *,
    num_threads: int = 6,
    seed: int = 20261601,
    ckpt_path: Path = DEFAULT_OCCAM_CKPT_PATH,
    occam_src_dir: Path = DEFAULT_OCCAM_SRC_DIR,
) -> Dict[str, Any]:
    """
    Gate 5A Item 2 (`PREREG_I_ADDENDUM_B_DRAFT_v2.md §5.1`):
      - n = 32 synthetic cells, seed = 20261601;
      - Strictly positive finite softplus outputs & dead_gene_frac < 0.01;
      - Exact CPU forward determinism (`max_abs == 0.0`);
      - Encoder-cache parity (`model(z_ctrl)` vs `model(X_ctrl)` == `0.0`);
      - Cross-process `PYTHONHASHSEED` (`1` vs `99991`) determinism (`max_abs == 0.0`);
      - Control-swap sensitivity (`X_A != X_B` -> `max_abs_delta_diff > 1e-4`);
      - Single-perturbation identity (`decode_seq_mean == sum` within `<= 1e-6`).
    """
    n_cells = 32
    n_enc = len(adapter.encoder_gene_names)
    rhj_cols = adapter.rhj_enc_cols

    x_a = _make_synthetic_controls(n_cells, n_enc, rhj_cols, seed=seed, scale=1.0)
    x_b = _make_synthetic_controls(n_cells, n_enc, rhj_cols, seed=seed + 777, scale=4.0)

    # Pick 4 deterministic trained RHJ perturbations
    sample_perts = sorted(adapter.rhj_eligible_keys)[:4]
    sample_pids = [adapter.pert2id[p] for p in sample_perts]
    pidx_one = [[sample_pids[0]]] * n_cells

    adapter.model.eval()
    with torch.no_grad():
        xt_a = torch.from_numpy(x_a).to("cpu")
        z_a = adapter.model.encoder(xt_a)
        out_sum_1 = adapter.model.forward_composed(xt_a, pidx_one, "sum")["prediction"].float().numpy()
        out_sum_2 = adapter.model.forward_composed(xt_a, pidx_one, "sum")["prediction"].float().numpy()
        out_dsm = adapter.model.forward_composed(xt_a, pidx_one, "decode_seq_mean")["prediction"].float().numpy()
        out_from_z = adapter.model.forward_composed(z_a, pidx_one, "decode_seq_mean")["prediction"].float().numpy()

    rhj_out = out_sum_1[:, rhj_cols]
    all_finite = bool(np.all(np.isfinite(out_sum_1)))
    min_val = float(out_sum_1.min())
    strictly_positive = bool(min_val > 0.0)
    dead_gene_frac_all = float(np.mean(out_sum_1.std(axis=0) == 0.0))
    dead_gene_frac_rhj = float(np.mean(rhj_out.std(axis=0) == 0.0))

    repeat_max_abs = float(np.max(np.abs(out_sum_1 - out_sum_2)))
    compose_max_abs = float(np.max(np.abs(out_dsm - out_sum_1)))
    z_cache_max_abs = float(np.max(np.abs(out_from_z - out_dsm)))

    # Control-swap sensitivity across the 4 perturbations (both reconstructed & real baselines)
    res_a = adapter.predict_deltas_from_controls(x_a, sample_pids, rhj_cols, batch_size=32)
    res_b = adapter.predict_deltas_from_controls(x_b, sample_pids, rhj_cols, batch_size=32)
    swap_recon_max_abs = float(np.max(np.abs(res_a["delta_reconstructed"] - res_b["delta_reconstructed"])))
    swap_real_max_abs = float(np.max(np.abs(res_a["delta_real"] - res_b["delta_real"])))

    # Cross-process PYTHONHASHSEED=1 vs 99991 check
    sub_outputs: List[np.ndarray] = []
    sub_code = (
        "import sys, numpy as np, torch; "
        f"sys.path.insert(0, {str(ADDB_SRC_DIR)!r}); "
        "from occam_adapter_addb import OccamInferenceAdapter; "
        "from run_gates_cpu_addb import _make_synthetic_controls; "
        f"adp, _ = OccamInferenceAdapter.load_verified(ckpt_path={str(ckpt_path)!r}, "
        f"occam_src_dir={str(occam_src_dir)!r}, device='cpu', num_threads={int(num_threads)}, "
        "verify_file_hashes=False); "
        f"x = _make_synthetic_controls({n_cells}, len(adp.encoder_gene_names), adp.rhj_enc_cols, seed={seed}); "
        f"res = adp.predict_deltas_from_controls(x, {sample_pids[:2]!r}, adp.rhj_enc_cols, batch_size=32); "
        "np.save(sys.argv[1], res['delta_reconstructed'])"
    )
    import tempfile

    with tempfile.TemporaryDirectory() as td:
        for hseed in ("1", "99991"):
            out_npy = Path(td) / f"delta_h{hseed}.npy"
            env = os.environ.copy()
            env["PYTHONHASHSEED"] = hseed
            env["OMP_NUM_THREADS"] = str(int(num_threads))
            env["MKL_NUM_THREADS"] = str(int(num_threads))
            proc = subprocess.run(
                [sys.executable, "-c", sub_code, str(out_npy)],
                env=env,
                capture_output=True,
                text=True,
            )
            if proc.returncode != 0:
                raise RuntimeError(
                    f"Gate 5A cross-process PYTHONHASHSEED={hseed} subprocess failed (rc={proc.returncode}): "
                    f"{proc.stderr[-1000:]}"
                )
            sub_outputs.append(np.load(out_npy))

    cross_proc_max_abs = float(np.max(np.abs(sub_outputs[0] - sub_outputs[1])))

    checks = {
        "all_finite": all_finite,
        "strictly_positive_softplus": strictly_positive,
        "dead_gene_frac_all_lt_0.01": dead_gene_frac_all < 0.01,
        "dead_gene_frac_rhj_lt_0.01": dead_gene_frac_rhj < 0.01,
        "cpu_repeat_determinism_exact_0.0": repeat_max_abs == 0.0,
        "z_cache_parity_exact_0.0": z_cache_max_abs == 0.0,
        "cross_process_pythonhashseed_exact_0.0": cross_proc_max_abs == 0.0,
        "control_swap_recon_gt_1e-4": swap_recon_max_abs > 1e-4,
        "control_swap_real_gt_1e-4": swap_real_max_abs > 1e-4,
        "decode_seq_mean_eq_sum_le_1e-6": compose_max_abs <= 1e-6,
    }
    passed = all(checks.values())
    return {
        "passed": passed,
        "n_cells": n_cells,
        "seed": int(seed),
        "sample_perts": sample_perts,
        "min_output_value": min_val,
        "dead_gene_frac_all": dead_gene_frac_all,
        "dead_gene_frac_rhj": dead_gene_frac_rhj,
        "repeat_max_abs": repeat_max_abs,
        "z_cache_max_abs": z_cache_max_abs,
        "cross_process_pythonhashseed_max_abs": cross_proc_max_abs,
        "control_swap_recon_max_abs": swap_recon_max_abs,
        "control_swap_real_max_abs": swap_real_max_abs,
        "compose_decode_seq_mean_vs_sum_max_abs": compose_max_abs,
        "checks": checks,
    }


def run_gate5b_hepg2_cpu(
    adapter: OccamInferenceAdapter,
    *,
    h5ad_path: Path = DEFAULT_LOCAL_NADIG_DIR / "GSE264667_hepg2_raw_singlecell_01.h5ad",
    gene_axis_csv: Path = DEFAULT_OCCAM_GENE_AXIS_CSV,
    perturbations_csv: Path = DEFAULT_OCCAM_PERTURBATIONS_CSV,
    splits_manifest_json: Path = DEFAULT_OCCAM_SPLITS_MANIFEST_JSON,
    batch_size: int = 512,
    base_seed: int = 42,
) -> Dict[str, Any]:
    """
    Execute Gate 5B (`PREREG_I_ADDENDUM_B_DRAFT_v2.md §5.2`, `Lock L1`) on CPU:
      - Exact `611` cohort perts / `4,976` controls / `42,402` panel rows;
      - Split-hash parity (`/20`) selects Case A (`0.005 / 0.002 / 0.002`) vs Case B (`0.015 / 0.005 / 0.015`);
      - Reproduces Occam's `splithalf` `ABL_N_ALL` `seed_42` `nocal` `injctrl` (`baseline="real"`) protocol
        across all 5 splits (`split_seed = 42 + 1000 * s`, `decode_seq_mean`) and checks all 3 locked metrics.
    """
    t0 = time.time()
    print(f"[{_now_iso()}] [Gate 5B] Building local nadig_hepg2 panel (cpu_t={_cpu_time_s():.1f}s)...", flush=True)
    axis = load_rhj_gene_axis(gene_axis_csv, verify_sha256=True)
    ptab = load_rhj_perturbations_table(perturbations_csv)
    manifest = json.loads(Path(splits_manifest_json).read_text(encoding="utf-8"))

    panel = build_gate5b_hepg2_panel(
        h5ad_path=h5ad_path,
        axis_ensembl=axis.ensembl_ids,
        eligible_test_perts=ptab["test_keys"],
        splits_manifest=manifest,
        load_expression=True,
        compute_experimental=True,
        enforce_lock_l1_counts=True,
    )
    assert panel.X_panel is not None and panel.experimental_metrics is not None

    # Precompute z_ctrl ONCE for all 4,976 unique control cells in the panel
    ctrl_panel_rows = np.flatnonzero(panel.labels == SPLITHALF_CONTROL_LABEL)
    row_to_ctrl_pos = np.full(panel.n_panel_rows, -1, dtype=np.int64)
    row_to_ctrl_pos[ctrl_panel_rows] = np.arange(len(ctrl_panel_rows), dtype=np.int64)

    t_enc0 = time.time()
    print(
        f"[{_now_iso()}] [Gate 5B] Encoding all {len(ctrl_panel_rows)} unique HepG2 control cells once "
        f"(cpu_t={_cpu_time_s():.1f}s)...",
        flush=True,
    )
    z_ctrl_all = adapter.encode_rhj_axis_matrix(panel.X_panel[ctrl_panel_rows], batch_size=batch_size)
    print(
        f"[{_now_iso()}] [Gate 5B] Encoded {z_ctrl_all.shape} in {time.time() - t_enc0:.1f}s "
        f"(cpu_t={_cpu_time_s():.1f}s).",
        flush=True,
    )

    # Verify all 611 cohort perturbations resolve cleanly in ckpt pert2id
    missing_p = [p for p in panel.cohort_perts if p not in adapter.pert2id]
    if missing_p:
        raise KeyError(f"Gate 5B cohort pert(s) absent from ckpt pert2id: {missing_p[:10]}")

    per_split_records: List[Dict[str, Any]] = []
    for sp_obj in panel.splits:
        ts = time.time()
        s_idx = int(sp_obj.split_index)
        split_seed = int(base_seed) + 1000 * s_idx

        # Exact splithalf.panel_for_split stable ordering of half-A perturbed rows
        a_rows_raw = np.asarray(sp_obj.pert_A, dtype=np.int64)
        a_order = np.argsort(panel.labels[a_rows_raw], kind="stable")
        a_rows = a_rows_raw[a_order]
        arm_labels = panel.labels[a_rows]
        arm_pids = [adapter.pert2id[str(p)] for p in arm_labels.tolist()]

        ctrl_b_rows = np.asarray(sp_obj.ctrl_B, dtype=np.int64)
        ctrl_b_z_pos = row_to_ctrl_pos[ctrl_b_rows]
        if np.any(ctrl_b_z_pos < 0):
            raise RuntimeError(f"Split {s_idx}: ctrl_B contains non-control panel rows.")

        # Exact predict_panel control draw (`seed = 42 + 1000 * s`)
        rng = np.random.default_rng(split_seed)
        cin_idx = rng.integers(0, len(ctrl_b_rows), size=len(a_rows))
        z_sel_pos = ctrl_b_z_pos[cin_idx]

        pred_rhj = adapter.predict_cell_pairs_from_z(
            z_ctrl_pool=z_ctrl_all,
            ctrl_indices=z_sel_pos,
            pert_ids=arm_pids,
            output_enc_cols=adapter.rhj_enc_cols,
            compose_mode="decode_seq_mean",
            batch_size=batch_size,
            progress_label=f"gate5b_split_{s_idx}",
        )

        # Score exact injctrl / baseline="real" nocal pearson_delta, mae, and mae_delta
        ctrl_a_mean = panel.X_panel[sp_obj.ctrl_A].mean(axis=0, dtype=np.float64)
        ctrl_b_mean = panel.X_panel[sp_obj.ctrl_B].mean(axis=0, dtype=np.float64)
        real_a_mat = panel.X_panel[a_rows]

        # Group indices once in O(N)
        loc_by_p: Dict[str, List[int]] = {}
        for i_loc, p_name in enumerate(arm_labels.tolist()):
            loc_by_p.setdefault(p_name, []).append(i_loc)

        pds: List[float] = []
        maes: List[float] = []
        mae_deltas: List[float] = []
        for p in panel.cohort_perts:
            locs = loc_by_p[p]
            pb = pred_rhj[locs].mean(axis=0, dtype=np.float64)
            rb = real_a_mat[locs].mean(axis=0, dtype=np.float64)
            dp = pb - ctrl_b_mean
            dr = rb - ctrl_a_mean
            dp_c = dp - dp.mean()
            dr_c = dr - dr.mean()
            denom = float(np.linalg.norm(dp_c) * np.linalg.norm(dr_c))
            if denom > 1e-12:
                pds.append(float(np.dot(dp_c, dr_c) / denom))
            maes.append(float(np.mean(np.abs(pb - rb))))
            mae_deltas.append(float(np.mean(np.abs(dp - dr))))

        sp_pd = float(np.mean(pds))
        sp_mae = float(np.mean(maes))
        sp_mae_d = float(np.mean(mae_deltas))
        dt_s = time.time() - ts
        print(
            f"[{_now_iso()}] [Gate 5B] split {s_idx}/4 done in {dt_s:.1f}s (cpu_t={_cpu_time_s():.1f}s): "
            f"pearson_delta={sp_pd:.6f}, mae={sp_mae:.6f}, mae_delta={sp_mae_d:.6f}",
            flush=True,
        )
        per_split_records.append(
            {
                "split": s_idx,
                "split_seed": split_seed,
                "n_pert_rows": int(len(a_rows)),
                "n_perts_scored": len(pds),
                "pearson_delta": sp_pd,
                "mae": sp_mae,
                "mae_delta": sp_mae_d,
                "elapsed_s": round(dt_s, 2),
            }
        )

    mean_pd = float(np.mean([r["pearson_delta"] for r in per_split_records]))
    mean_mae = float(np.mean([r["mae"] for r in per_split_records]))
    mean_mae_delta = float(np.mean([r["mae_delta"] for r in per_split_records]))
    exp_pd = float(panel.experimental_metrics["experimental_pearson_delta_mean"])

    mcmp = panel.manifest_comparison
    tols = mcmp["gate5b_tolerances"]
    diff_pd = abs(mean_pd - GATE5B_AUTHOR_MODEL_PEARSON_DELTA)
    diff_mae = abs(mean_mae - GATE5B_AUTHOR_MODEL_MAE)
    diff_exp = abs(exp_pd - GATE5B_AUTHOR_EXP_PEARSON_DELTA)

    pass_pd = diff_pd <= float(tols["pearson_delta"])
    pass_mae = diff_mae <= float(tols["mae"])
    pass_exp = diff_exp <= float(tols["experimental_pearson_delta"])
    passed = bool(pass_pd and pass_mae and pass_exp and mcmp.get("max_ctrl_batch_abs_diff", 99) <= 1)

    return {
        "passed": passed,
        "elapsed_total_s": round(time.time() - t0, 2),
        "h5ad_path": str(h5ad_path),
        "n_cohort_perts": len(panel.cohort_perts),
        "n_controls": int(panel.n_controls),
        "n_panel_rows": int(panel.n_panel_rows),
        "cohort_sha256": panel.cohort_sha256,
        "tolerance_case": str(mcmp["tolerance_case"]),
        "split_hash_matches": int(mcmp["n_hash_matches"]),
        "split_hash_total": int(mcmp["n_hash_total"]),
        "max_ctrl_batch_abs_diff": int(mcmp["max_ctrl_batch_abs_diff"]),
        "tolerances": tols,
        "metrics": {
            "adapter_pearson_delta_mean": mean_pd,
            "target_pearson_delta_mean": GATE5B_AUTHOR_MODEL_PEARSON_DELTA,
            "abs_diff_pearson_delta": diff_pd,
            "pass_pearson_delta": pass_pd,
            "adapter_mae_mean": mean_mae,
            "target_mae_mean": GATE5B_AUTHOR_MODEL_MAE,
            "abs_diff_mae": diff_mae,
            "pass_mae": pass_mae,
            "adapter_mae_delta_mean": mean_mae_delta,
            "target_mae_delta_mean": GATE5B_AUTHOR_MODEL_MAE_DELTA,
            "experimental_pearson_delta_mean": exp_pd,
            "target_experimental_pearson_delta_mean": GATE5B_AUTHOR_EXP_PEARSON_DELTA,
            "abs_diff_experimental_pearson_delta": diff_exp,
            "pass_experimental_pearson_delta": pass_exp,
        },
        "per_split": per_split_records,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description="Run Round 16 Addendum B CPU Gates 5A (items 1-2) and 5B.")
    ap.add_argument("--num-threads", type=int, default=6)
    ap.add_argument("--batch-size", type=int, default=512)
    ap.add_argument("--ckpt-path", type=Path, default=DEFAULT_OCCAM_CKPT_PATH)
    ap.add_argument("--occam-src-dir", type=Path, default=DEFAULT_OCCAM_SRC_DIR)
    ap.add_argument("--gene-axis-csv", type=Path, default=DEFAULT_OCCAM_GENE_AXIS_CSV)
    ap.add_argument("--perturbations-csv", type=Path, default=DEFAULT_OCCAM_PERTURBATIONS_CSV)
    ap.add_argument("--splits-manifest-json", type=Path, default=DEFAULT_OCCAM_SPLITS_MANIFEST_JSON)
    ap.add_argument("--checksums-file", type=Path, default=DEFAULT_OCCAM_CHECKSUMS_FILE)
    ap.add_argument(
        "--nadig-hepg2-h5ad",
        type=Path,
        default=DEFAULT_LOCAL_NADIG_DIR / "GSE264667_hepg2_raw_singlecell_01.h5ad",
    )
    ap.add_argument(
        "--out-json",
        type=Path,
        default=ADDB_ROOT / "results" / "gates_cpu_summary.json",
    )
    args = ap.parse_args()

    configure_occam_runtime(device="cpu", num_threads=args.num_threads, allow_mps=False)
    print(
        f"[{_now_iso()}] Starting Addendum B CPU Gates 5A (items 1-2) + 5B | "
        f"python={sys.version.split()[0]} torch={torch.__version__} threads={torch.get_num_threads()}",
        flush=True,
    )

    adapter, g5a_1 = OccamInferenceAdapter.load_verified(
        ckpt_path=args.ckpt_path,
        occam_src_dir=args.occam_src_dir,
        gene_axis_csv=args.gene_axis_csv,
        perturbations_csv=args.perturbations_csv,
        checksums_file=args.checksums_file,
        device="cpu",
        num_threads=args.num_threads,
        allow_mps=False,
        batch_size=args.batch_size,
        verify_file_hashes=True,
    )
    print(f"[{_now_iso()}] Gate 5A Item 1 PASSED.", flush=True)

    g5a_2 = run_gate5a_item2_synthetic(
        adapter,
        num_threads=args.num_threads,
        ckpt_path=args.ckpt_path,
        occam_src_dir=args.occam_src_dir,
    )
    if not g5a_2["passed"]:
        summary_fail = {
            "timestamp_utc": _now_iso(),
            "device": "cpu",
            "num_threads": int(args.num_threads),
            "python_version": sys.version.split()[0],
            "torch_version": str(torch.__version__),
            "gate5a_item1": dataclasses.asdict(g5a_1),
            "gate5a_item2": g5a_2,
            "gate5b": None,
            "overall_passed": False,
        }
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(json.dumps(summary_fail, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print("GATE 5A ITEM 2 FAILED — HALTING.", file=sys.stderr)
        return 1
    print(f"[{_now_iso()}] Gate 5A Item 2 PASSED.", flush=True)

    g5b = run_gate5b_hepg2_cpu(
        adapter,
        h5ad_path=args.nadig_hepg2_h5ad,
        gene_axis_csv=args.gene_axis_csv,
        perturbations_csv=args.perturbations_csv,
        splits_manifest_json=args.splits_manifest_json,
        batch_size=args.batch_size,
    )

    overall_ok = bool(g5a_1.passed and g5a_2["passed"] and g5b["passed"])
    summary = {
        "timestamp_utc": _now_iso(),
        "device": "cpu",
        "num_threads": int(args.num_threads),
        "python_version": sys.version.split()[0],
        "torch_version": str(torch.__version__),
        "gate5a_item1": dataclasses.asdict(g5a_1),
        "gate5a_item2": g5a_2,
        "gate5a_item3_mps_parity": "DEFERRED_UNTIL_AT3_FINISHES_R16",
        "gate5b": g5b,
        "overall_passed": overall_ok,
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if overall_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
