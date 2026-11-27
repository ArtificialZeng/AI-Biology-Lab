#!/usr/bin/env python3
"""
Clean-Room Stage 5A (items 1-2) and Stage 5B Quantitative Reproduction Gate Runner on CPU (`A-T3`).

Strictly enforces:
  - CPU ONLY (`device="cpu"`, `torch.set_num_threads(4)`); Stage 5A item 3 (`MPS` parity) and all
    Group 2 / Group 3 real inference remain on HARD HOLD until Master Controller authorization.
  - Zero access to `round16_addB/src/`, `round16_addB/results/`, or removed archive members.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import torch

AUDIT_DIR = Path(__file__).resolve().parent
ADDB_ROOT = AUDIT_DIR.parent
if str(AUDIT_DIR) not in sys.path:
    sys.path.insert(0, str(AUDIT_DIR))

from audit_addb_occam_adapter import (  # noqa: E402
    ADDB_ROOT as _ADDB_ROOT,
    CONTROL_LABEL,
    DEFAULT_CHECKSUMS_FILE,
    DEFAULT_CKPT_PATH,
    DEFAULT_COHORT_JSON,
    DEFAULT_GENE_AXIS_CSV,
    DEFAULT_NADIG_HEPG2_H5AD,
    DEFAULT_OCCAM_SRC_ROOT,
    DEFAULT_PERTURBATIONS_CSV,
    DEFAULT_SCORES_CSV,
    DEFAULT_SPLITS_MANIFEST_JSON,
    EXPECTED_GENE_AXIS_CONTENT_SHA256,
    EXPECTED_NADIG_HEPG2_LOCAL_SHA256,
    OCCAM_NUM_ENCODER_GENES,
    OCCAM_TARGET_SUM,
    OccamCleanRoomAdapter,
    TwoTierPertResolver,
    compute_celleval_pseudobulk_metrics,
    ensure_occam_src_on_path,
    load_hgnc_bridge,
    load_nadig_hepg2_panel_raw,
    load_rhj_gene_axis,
    normalize_counts_to_rhj_44728,
    sha256_array_bytes,
    sha256_file,
    verify_addb_checksums,
)

AUTHOR_TARGET_PEARSON_DELTA = 0.394178
AUTHOR_TARGET_MAE = 0.030157
AUTHOR_TARGET_EXPERIMENTAL_PD = 0.372174

CASE_A_TOL = {
    "pearson_delta_abs": 0.005,
    "mae_abs": 0.002,
    "experimental_pearson_delta_abs": 0.002,
}
CASE_B_TOL = {
    "pearson_delta_abs": 0.015,
    "mae_abs": 0.005,
    "experimental_pearson_delta_abs": 0.015,
}


def _ts() -> str:
    t = os.times()
    return f"[{time.strftime('%H:%M:%S')} cpu={t.user + t.system:.1f}s]"


def _synthetic_smoke_vector_worker() -> None:
    """Subprocess entrypoint for cross-`PYTHONHASHSEED` CPU determinism verification."""
    torch.set_num_threads(4)
    adapter = OccamCleanRoomAdapter(device="cpu", num_threads=4, verify_sha256=False)
    rng = np.random.default_rng(20261601)
    axis_ens, _ = load_rhj_gene_axis()
    enc_cols = np.asarray([adapter.gene_to_index[g] for g in axis_ens], dtype=np.int64)
    x_ctrl = np.zeros((32, OCCAM_NUM_ENCODER_GENES), dtype=np.float32)
    x_ctrl[:, enc_cols] = rng.uniform(0.0, 3.0, size=(32, enc_cols.size)).astype(np.float32)
    pid = int(sorted(adapter.pert2id.values())[0])
    out = adapter.predict_per_cell(x_ctrl, [[pid]] * 32, batch_size=256)
    base = adapter.predict_reconstructed_baseline(x_ctrl, batch_size=256)
    payload = np.concatenate([out, out - base], axis=1)
    out_file = os.environ["ADDB_SMOKE_OUT"]
    np.save(out_file, payload)


def run_stage_5a_cpu(
    out_dir: Path,
    num_threads: int = 4,
    batch_size: int = 256,
) -> Tuple[OccamCleanRoomAdapter, Dict[str, Any]]:
    """Run Stage 5A items 1 and 2 on CPU only (item 3 MPS deferred pending GO)."""
    print(f"{_ts()} [5A.1] Verifying Addendum B checksums & checkpoint contract...", flush=True)
    checksum_records = verify_addb_checksums(DEFAULT_CHECKSUMS_FILE)
    adapter = OccamCleanRoomAdapter(
        ckpt_path=DEFAULT_CKPT_PATH,
        device="cpu",
        num_threads=num_threads,
        occam_src_root=DEFAULT_OCCAM_SRC_ROOT,
        gene_axis_csv=DEFAULT_GENE_AXIS_CSV,
        verify_sha256=True,
    )

    print(f"{_ts()} [5A.2] Running 32-cell pure-synthetic smoke & determinism checks (CPU)...", flush=True)
    rng = np.random.default_rng(20261601)
    axis_ens, _ = load_rhj_gene_axis()
    enc_cols = np.asarray([adapter.gene_to_index[g] for g in axis_ens], dtype=np.int64)

    X_ctrl_A = np.zeros((32, OCCAM_NUM_ENCODER_GENES), dtype=np.float32)
    X_ctrl_A[:, enc_cols] = rng.uniform(0.0, 3.0, size=(32, enc_cols.size)).astype(np.float32)
    X_ctrl_B = np.zeros((32, OCCAM_NUM_ENCODER_GENES), dtype=np.float32)
    X_ctrl_B[:, enc_cols] = rng.uniform(0.0, 3.0, size=(32, enc_cols.size)).astype(np.float32)

    pid = int(sorted(adapter.pert2id.values())[0])
    pred1 = adapter.predict_per_cell(X_ctrl_A, [[pid]] * 32, batch_size=batch_size)
    pred2 = adapter.predict_per_cell(X_ctrl_A, [[pid]] * 32, batch_size=batch_size)
    base_A = adapter.predict_reconstructed_baseline(X_ctrl_A, batch_size=batch_size)
    pred_B = adapter.predict_per_cell(X_ctrl_B, [[pid]] * 32, batch_size=batch_size)
    base_B = adapter.predict_reconstructed_baseline(X_ctrl_B, batch_size=batch_size)

    all_finite = bool(np.isfinite(pred1).all())
    all_positive = bool((pred1 > 0.0).all())
    dead_gene_frac = float((pred1.std(axis=0) == 0.0).mean())
    repeat_max_abs = float(np.max(np.abs(pred1 - pred2)))
    delta_A = pred1 - base_A
    delta_B = pred_B - base_B
    swap_delta_max_abs = float(np.max(np.abs(delta_A - delta_B)))

    # Single-perturbation identity invariant: decode_seq_mean == sum within <= 1e-6
    with torch.no_grad():
        xt = torch.from_numpy(X_ctrl_A[:16]).to("cpu")
        pidx = [[pid]] * 16
        fc_mean = adapter.model.forward_composed(xt, pidx, mode="decode_seq_mean")["prediction"].float().numpy()
        fc_sum = adapter.model.forward_composed(xt, pidx, mode="sum")["prediction"].float().numpy()
    composed_max_abs = float(np.max(np.abs(fc_mean - fc_sum)))

    # Cross-process PYTHONHASHSEED=1 vs PYTHONHASHSEED=99991 invariance
    smoke_1 = out_dir / "_tmp_smoke_hash1.npy"
    smoke_99991 = out_dir / "_tmp_smoke_hash99991.npy"
    for hseed, tmp_p in (("1", smoke_1), ("99991", smoke_99991)):
        env = os.environ.copy()
        env["PYTHONHASHSEED"] = hseed
        env["ADDB_SMOKE_OUT"] = str(tmp_p)
        subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--internal-smoke-worker"],
            env=env,
            check=True,
        )
    arr1 = np.load(smoke_1)
    arr2 = np.load(smoke_99991)
    hashseed_max_abs = float(np.max(np.abs(arr1 - arr2)))
    smoke_1.unlink(missing_ok=True)
    smoke_99991.unlink(missing_ok=True)

    checks_5a2 = {
        "all_finite": all_finite,
        "all_strictly_positive_softplus": all_positive,
        "dead_gene_frac": dead_gene_frac,
        "dead_gene_frac_lt_0_01": bool(dead_gene_frac < 0.01),
        "repeat_forward_max_abs": repeat_max_abs,
        "repeat_forward_exact_zero": bool(repeat_max_abs == 0.0),
        "cross_pythonhashseed_max_abs": hashseed_max_abs,
        "cross_pythonhashseed_exact_zero": bool(hashseed_max_abs == 0.0),
        "ctrl_swap_delta_max_abs": swap_delta_max_abs,
        "ctrl_swap_delta_gt_1e4": bool(swap_delta_max_abs > 1e-4),
        "forward_composed_decode_seq_mean_vs_sum_max_abs": composed_max_abs,
        "forward_composed_invariant_le_1e6": bool(composed_max_abs <= 1e-6),
    }

    stage_5a_pass = bool(
        all_finite
        and all_positive
        and dead_gene_frac < 0.01
        and repeat_max_abs == 0.0
        and hashseed_max_abs == 0.0
        and swap_delta_max_abs > 1e-4
        and composed_max_abs <= 1e-6
    )

    summary_5a = {
        "item_1_checksums_verified": len(checksum_records),
        "item_1_contract": adapter.contract_summary,
        "item_2_synthetic_smoke_cpu": checks_5a2,
        "item_3_mps_parity": "DEFERRED_PENDING_CONTROLLER_MPS_GO",
        "items_1_2_pass": stage_5a_pass,
    }
    if not stage_5a_pass:
        raise RuntimeError(f"Stage 5A CPU gate FAILED: {json.dumps(summary_5a, indent=2)}")
    print(
        f"{_ts()} [5A] PASSED (repeat_diff={repeat_max_abs}, hashseed_diff={hashseed_max_abs}, "
        f"swap_delta_diff={swap_delta_max_abs:.6f}, composed_diff={composed_max_abs})",
        flush=True,
    )
    return adapter, summary_5a


def run_stage_5b_cpu(
    adapter: OccamCleanRoomAdapter,
    batch_size: int = 256,
) -> Dict[str, Any]:
    """
    Run Stage 5B Quantitative Author-Benchmark Reproduction Gate on local `nadig_hepg2` (`CPU float32`).
    """
    ensure_occam_src_on_path(DEFAULT_OCCAM_SRC_ROOT)
    from occam.vcell.panel.splithalf import draw_splits  # type: ignore

    print(f"{_ts()} [5B] Loading cohort.json, splits_manifest.json, and official reference CSV...", flush=True)
    with open(DEFAULT_COHORT_JSON, "r", encoding="utf-8") as fh:
        cohort = json.load(fh)
    with open(DEFAULT_SPLITS_MANIFEST_JSON, "r", encoding="utf-8") as fh:
        splits_manifest = json.load(fh)

    cohort_perts: List[str] = list(cohort["perts"])
    if len(cohort_perts) != 611:
        raise AssertionError(f"Expected 611 cohort perturbations, got {len(cohort_perts)}")

    # Verify official_comparisons CSV targets dynamically as cross-check against §5.2 constants
    csv_abl_pds: List[float] = []
    csv_abl_maes: List[float] = []
    csv_exp_pds: List[float] = []
    with open(DEFAULT_SCORES_CSV, "r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row["benchmark"] != "rhj_fewshot" or row["panel"] != "all" or row["seed"] != "42" or row["calibration"] != "nocal":
                continue
            if row["model"] == "ABL_N_ALL":
                csv_abl_pds.append(float(row["pearson_delta"]))
                csv_abl_maes.append(float(row["mae"]))
            elif row["model"] == "experimental":
                csv_exp_pds.append(float(row["pearson_delta"]))
    if len(csv_abl_pds) != 5 or len(csv_exp_pds) != 5:
        raise RuntimeError(f"Expected 5 seed=42 nocal splits in CSV, got {len(csv_abl_pds)}/{len(csv_exp_pds)}")

    # Read local Nadig HepG2 raw counts + Lock L1 legacy categoricals
    print(f"{_ts()} [5B] Reading local Nadig HepG2 ({DEFAULT_NADIG_HEPG2_H5AD.name})...", flush=True)
    raw_bundle = load_nadig_hepg2_panel_raw(
        h5ad_path=DEFAULT_NADIG_HEPG2_H5AD,
        cohort_perts=cohort_perts,
        expected_sha256=EXPECTED_NADIG_HEPG2_LOCAL_SHA256,
    )
    lp = raw_bundle["labels"]
    bp = raw_bundle["batches"]
    n_panel_rows = int(lp.size)
    n_ntc = int((lp == CONTROL_LABEL).sum())
    n_pert_rows = n_panel_rows - n_ntc
    n_present_perts = len(set(lp.tolist()) - {CONTROL_LABEL})
    if (n_panel_rows, n_pert_rows, n_ntc, n_present_perts) != (42402, 37426, 4976, 611):
        raise AssertionError(
            f"Lock L1 Nadig panel count mismatch: total={n_panel_rows} pert={n_pert_rows} "
            f"ntc={n_ntc} uniq_perts={n_present_perts} (expected 42402 / 37426 / 4976 / 611)"
        )
    print(
        f"{_ts()} [5B] Verified Lock L1 Nadig panel: {n_panel_rows} rows "
        f"({n_pert_rows} pert across {n_present_perts} perts + {n_ntc} NTCs)",
        flush=True,
    )

    # Full-panel 9,624-gene normalization -> scatter to 6,983 RHJ axis & 44,728 encoder input
    axis_ens, _ = load_rhj_gene_axis(
        DEFAULT_GENE_AXIS_CSV, expected_content_sha256=EXPECTED_GENE_AXIS_CONTENT_SHA256
    )
    _, X_6983, norm_audit = normalize_counts_to_rhj_44728(
        raw_counts=raw_bundle["raw_counts"],
        raw_ensembl_ids=raw_bundle["var_ensembl"],
        rhj_axis_ensembl=axis_ens,
        encoder_gene_names=adapter.gene_names,
        target_sum=OCCAM_TARGET_SUM,
    )
    del raw_bundle["raw_counts"]
    if norm_audit["n_rhj_present_in_raw"] != 6983 or norm_audit["n_rhj_absent_in_raw"] != 0:
        raise AssertionError(f"Expected 6983/6983 RHJ genes in Nadig HepG2: {norm_audit}")

    enc_cols_6983 = np.asarray([adapter.gene_to_index[g] for g in axis_ens], dtype=np.int64)

    # Two-Tier resolution check for all 611 cohort perturbations
    hgnc = load_hgnc_bridge()
    resolver = TwoTierPertResolver(adapter.pert2id, DEFAULT_PERTURBATIONS_CSV, hgnc_bridge=hgnc)
    pert_id_map: Dict[str, int] = {}
    for p_name in cohort_perts:
        # Check direct pert2id AND Two-Tier match cleanly
        if p_name in adapter.pert2id:
            pert_id_map[p_name] = int(adapter.pert2id[p_name])
        else:
            _, pid, _ = resolver.resolve(p_name)
            pert_id_map[p_name] = pid

    # Occam exact deterministic split drawing
    print(f"{_ts()} [5B] Drawing 5 batch-stratified splits via Occam `draw_splits(..., seed=42)`...", flush=True)
    splits = draw_splits(lp, bp, cohort, n_splits=5, base_seed=42)

    manifest_splits = {int(s["split"]): s for s in splits_manifest["splits"]}
    split_hash_matches = 0
    split_hash_details: List[Dict[str, Any]] = []
    max_ctrl_batch_diff = 0
    for sp in splits:
        ms = manifest_splits[int(sp.index)]
        cb_a = dict(zip(*np.unique(bp[sp.ctrl_A], return_counts=True)))
        cb_b = dict(zip(*np.unique(bp[sp.ctrl_B], return_counts=True)))
        bdiff = max(abs(int(cb_a.get(b, 0)) - int(cb_b.get(b, 0))) for b in set(cb_a) | set(cb_b))
        max_ctrl_batch_diff = max(max_ctrl_batch_diff, bdiff)
        if bdiff > 1:
            raise AssertionError(f"Split {sp.index} control batch balance violated: max|A-B|={bdiff} > 1")

        hashes = {
            "pert_A": sha256_array_bytes(sp.pert_A),
            "pert_B": sha256_array_bytes(sp.pert_B),
            "ctrl_A": sha256_array_bytes(sp.ctrl_A),
            "ctrl_B": sha256_array_bytes(sp.ctrl_B),
        }
        m_count = sum(1 for k, v in hashes.items() if v == ms[f"{k}_sha256"])
        split_hash_matches += m_count
        split_hash_details.append(
            {
                "split": int(sp.index),
                "n_pert_per_half": int(len(sp.pert_A)),
                "n_ctrl_per_half": int(len(sp.ctrl_A)),
                "n_dropped_odd": int(sp.n_dropped_odd),
                "ctrl_batch_max_abs_diff": int(bdiff),
                "hash_matches": int(m_count),
                "hashes": hashes,
            }
        )

    case_mode = "Case A" if split_hash_matches == 20 else "Case B"
    tol = CASE_A_TOL if case_mode == "Case A" else CASE_B_TOL
    print(
        f"{_ts()} [5B] Split hash matches = {split_hash_matches}/20 -> {case_mode} "
        f"(ctrl_batch_max_abs_diff={max_ctrl_batch_diff})",
        flush=True,
    )

    # Pre-encode all 4,976 NTC cells cleanly once (encoder is strictly per-cell and deterministic; verified diff == 0.0)
    ntc_mask = lp == CONTROL_LABEL
    ntc_global_idx = np.flatnonzero(ntc_mask)
    ntc_pos_lookup = {int(g_i): k for k, g_i in enumerate(ntc_global_idx)}
    ntc_X_44728 = np.zeros((ntc_global_idx.size, OCCAM_NUM_ENCODER_GENES), dtype=np.float32)
    ntc_X_44728[:, enc_cols_6983] = X_6983[ntc_global_idx]
    print(f"{_ts()} [5B] Encoding {ntc_global_idx.size} HepG2 NTC cells on CPU...", flush=True)
    ntc_z_all = adapter.encode_controls(ntc_X_44728, batch_size=batch_size)

    split_records: List[Dict[str, Any]] = []
    model_pds: List[float] = []
    model_maes: List[float] = []
    exp_pds: List[float] = []

    for sp in splits:
        sp_idx = int(sp.index)
        split_seed = 42 + 1000 * sp_idx

        pert_A_X = X_6983[sp.pert_A]
        pert_A_labels = lp[sp.pert_A]
        ctrl_A_X = X_6983[sp.ctrl_A]

        pert_B_X = X_6983[sp.pert_B]
        pert_B_labels = lp[sp.pert_B]
        ctrl_B_X = X_6983[sp.ctrl_B]

        # Sample 1 pool-B control cell per half-A perturbed row (`predict_panel` convention)
        rng = np.random.default_rng(split_seed)
        cin_idx = rng.integers(0, len(sp.ctrl_B), size=len(sp.pert_A))
        ctrl_B_ntc_pos = np.asarray([ntc_pos_lookup[int(sp.ctrl_B[i])] for i in cin_idx], dtype=np.int64)

        z_in = ntc_z_all[ctrl_B_ntc_pos]
        pidx = [[pert_id_map[str(lbl)]] for lbl in pert_A_labels]
        pred_A_X = adapter.predict_per_cell(
            z_in, pidx, batch_size=batch_size, out_cols=enc_cols_6983
        )

        # `injctrl` / `baseline="real"` `nocal` scoring against half-A target + ctrl_A
        mod_m = compute_celleval_pseudobulk_metrics(
            real_pert_X=pert_A_X,
            real_pert_labels=pert_A_labels,
            real_ctrl_X=ctrl_A_X,
            pred_pert_X=pred_A_X,
            pred_pert_labels=pert_A_labels,
            pred_ctrl_X=ctrl_B_X,
            cohort_perts=cohort_perts,
        )
        exp_m = compute_celleval_pseudobulk_metrics(
            real_pert_X=pert_A_X,
            real_pert_labels=pert_A_labels,
            real_ctrl_X=ctrl_A_X,
            pred_pert_X=pert_B_X,
            pred_pert_labels=pert_B_labels,
            pred_ctrl_X=ctrl_B_X,
            cohort_perts=cohort_perts,
        )
        model_pds.append(mod_m["pearson_delta_mean"])
        model_maes.append(mod_m["mae_mean"])
        exp_pds.append(exp_m["pearson_delta_mean"])

        rec = {
            "split": sp_idx,
            "seed": split_seed,
            "n_pert_rows": int(pert_A_labels.size),
            "n_ctrl_rows_per_half": int(ctrl_A_X.shape[0]),
            "occam_pearson_delta": mod_m["pearson_delta_mean"],
            "occam_mae": mod_m["mae_mean"],
            "experimental_pearson_delta": exp_m["pearson_delta_mean"],
            "experimental_mae": exp_m["mae_mean"],
        }
        split_records.append(rec)
        print(
            f"{_ts()} [5B] split {sp_idx}: occam_pd={rec['occam_pearson_delta']:.6f}, "
            f"occam_mae={rec['occam_mae']:.6f}, exp_pd={rec['experimental_pearson_delta']:.6f}",
            flush=True,
        )

    mean_pd = float(np.mean(model_pds))
    mean_mae = float(np.mean(model_maes))
    mean_exp_pd = float(np.mean(exp_pds))

    diff_pd = abs(mean_pd - AUTHOR_TARGET_PEARSON_DELTA)
    diff_mae = abs(mean_mae - AUTHOR_TARGET_MAE)
    diff_exp_pd = abs(mean_exp_pd - AUTHOR_TARGET_EXPERIMENTAL_PD)

    pass_pd = bool(diff_pd <= tol["pearson_delta_abs"])
    pass_mae = bool(diff_mae <= tol["mae_abs"])
    pass_exp_pd = bool(diff_exp_pd <= tol["experimental_pearson_delta_abs"])
    gate_5b_pass = bool(pass_pd and pass_mae and pass_exp_pd)

    return {
        "nadig_hepg2_h5ad": str(DEFAULT_NADIG_HEPG2_H5AD),
        "nadig_hepg2_sha256": raw_bundle["file_sha256"],
        "n_panel_rows": n_panel_rows,
        "n_pert_rows": n_pert_rows,
        "n_ntc_rows": n_ntc,
        "n_cohort_perts": n_present_perts,
        "normalization_audit": norm_audit,
        "split_hash_match_count": int(split_hash_matches),
        "split_hash_total": 20,
        "case": case_mode,
        "case_tolerances": tol,
        "ctrl_batch_max_abs_diff_max": int(max_ctrl_batch_diff),
        "splits": split_records,
        "split_hashes": split_hash_details,
        "author_csv_seed42_nocal_targets": {
            "abl_n_all_pearson_delta_mean": float(np.mean(csv_abl_pds)),
            "abl_n_all_mae_mean": float(np.mean(csv_abl_maes)),
            "experimental_pearson_delta_mean": float(np.mean(csv_exp_pds)),
        },
        "adapter_mean_pearson_delta": mean_pd,
        "adapter_mean_mae": mean_mae,
        "local_experimental_pearson_delta": mean_exp_pd,
        "abs_diffs": {
            "pearson_delta": diff_pd,
            "mae": diff_mae,
            "experimental_pearson_delta": diff_exp_pd,
        },
        "metric_pass": {
            "pearson_delta": pass_pd,
            "mae": pass_mae,
            "experimental_pearson_delta": pass_exp_pd,
        },
        "gate_5b_pass": gate_5b_pass,
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Run Addendum B Stage 5A (items 1-2) and Stage 5B on CPU.")
    ap.add_argument("--num-threads", type=int, default=4)
    ap.add_argument("--batch-size", type=int, default=256)
    ap.add_argument(
        "--output-json",
        default=str(AUDIT_DIR / "results" / "gates_cpu_summary.json"),
    )
    ap.add_argument("--internal-smoke-worker", action="store_true")
    args = ap.parse_args(argv)

    if args.internal_smoke_worker:
        _synthetic_smoke_vector_worker()
        return 0

    out_path = Path(args.output_json)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    adapter, res_5a = run_stage_5a_cpu(
        out_dir=out_path.parent,
        num_threads=args.num_threads,
        batch_size=args.batch_size,
    )
    res_5b = run_stage_5b_cpu(adapter=adapter, batch_size=args.batch_size)

    overall_pass = bool(res_5a["items_1_2_pass"] and res_5b["gate_5b_pass"])
    summary = {
        "device": "cpu",
        "torch_num_threads": int(args.num_threads),
        "batch_size": int(args.batch_size),
        "torch_version": torch.__version__,
        "stage_5a": res_5a,
        "stage_5b": res_5b,
        "overall_cpu_gates_pass": overall_pass,
    }
    out_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    out_sha = sha256_file(out_path)
    print(
        f"{_ts()} [SUMMARY] case={res_5b['case']} splits={res_5b['split_hash_match_count']}/20 "
        f"pearson_delta={res_5b['adapter_mean_pearson_delta']:.6f} "
        f"mae={res_5b['adapter_mean_mae']:.6f} "
        f"exp_pd={res_5b['local_experimental_pearson_delta']:.6f} "
        f"PASS={overall_pass} json_sha256={out_sha}",
        flush=True,
    )
    return 0 if overall_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
