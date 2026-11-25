#!/usr/bin/env python3
"""
Clean-Room Round 16 Real-Data & Synthetic-Parity Evaluation Runner (`A-T3`, `PREREG_I`).

Enforces strict runtime blinding & schema hard gates BEFORE accessing real data:
  1. `PREREG_I_COMMIT.txt` exists and is non-empty.
  2. `CODE_FREEZE_R16.txt` exists and is non-empty.
  3. `R16_GO.txt` exists, is non-empty, and postdates `R16_GATE2_SCHEMA.md` when present.

Implements all Gate-2 (`R16_GATE2_SCHEMA.md`) structural facts:
  - `SE-600M` direct `safetensors` load + strict `float32` on Apple `MPS` (or CPU fallback).
  - `ST-SE`: `clip(gene_decoder(preds), 0.0, 14.0)`.
  - `ST-HVG`: `clip(preds, 0.0, 14.0)` directly (`gene_decoder` never trained).
  - `16 × 32` cell sets with `padded=True` and fixed `batch = 0` identical across matched & mismatched arms.
  - `delta_pred(p) = pred(p) - pred("non-targeting")` on identical sampled input controls.
  - Frozen R10 `HGNCBridge` symbol alignment across genes and targets.
  - Exports per-line, per-arm `.npz` bundles (predicted delta matrices, `delta_obs`, split-half
    `delta_half_a`/`delta_half_b`, `Transfer_Mean`, `FCT`, aligned gene/pert arrays) for
    three-way bitwise comparison alongside full metric & adjudicator JSON reports.

Once `PREREG_I v5` locks, only `PREREG_I_V5_CONFIG` constants (`N_eval`, `G_eval`, `hvg_preproc_mode`)
need updating. Strictly NEVER reads or imports `round16/src/`.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")

AUDIT_DIR = Path(__file__).resolve().parent
REPO_ROOT = AUDIT_DIR.parents[1]
if str(AUDIT_DIR) not in sys.path:
    sys.path.insert(0, str(AUDIT_DIR))

from audit_r16_adjudicate import (  # noqa: E402
    DEFAULT_TOST_DELTA,
    adjudicate_claim_i1_line,
    adjudicate_claim_i2_line,
    adjudicate_group_r16,
)
from audit_r16_metrics import (  # noqa: E402
    DEFAULT_R16_SEED,
    compute_orth_suite_r16,
    compute_paired_contrast_r16,
    compute_paired_model_vs_all_mismatched_r16,
    compute_winner_curse_weight_shift_r16,
)
from audit_r16_raw_controls import (  # noqa: E402
    extract_raw_control_pool_r16,
    load_line_eval_reference_r16,
)
from audit_r16_sampling import (  # noqa: E402
    LOCAL_RAW_CONTROL_SOURCES,
    build_mismatched_strata_r16,
)
from audit_r16_state_adapter import (  # noqa: E402
    CleanRoomArcStateAdapter,
    CleanRoomSE600MEncoder,
    MODELS_ROOT_DEFAULT,
    assert_freeze_gates_r16,
    get_hgnc_bridge,
)

STATUS_DIR_DEFAULT = REPO_ROOT / "round16" / "status"
OUTPUT_DIR_DEFAULT = AUDIT_DIR / "results"


# Pending final lock once `PREREG_I v5` is committed by Master Controller / Main Line.
PREREG_I_V5_CONFIG: Dict[str, Any] = {
    "prereg_version": "v4-GATE2-PENDING-V5",
    "n_ctrl": 512,
    "cell_set_len": 32,
    "seed": DEFAULT_R16_SEED,
    "n_boot": 4000,
    "tost_delta": DEFAULT_TOST_DELTA,
    "neutral_batch_mode_primary": "index_zero",  # fixed batch = 0
    "delta_mode_primary": "minus_pred_ntc",      # pred(p) - pred("non-targeting")
    # Update after PREREG_I v5 lock ("official_median_log1p" vs "libnorm1e4_log1p" per Controller ruling)
    "hvg_preproc_primary": "official_median_log1p",
    "group1_lines": ("nadig_hepg2", "nadig_jurkat"),
    "group2_crispri_lines": ("dld1", "cd4t_rest"),
    "group2_crisprko_lines": ("a549_ko",),
    "control_pool_lines": LOCAL_RAW_CONTROL_SOURCES,
    # Expected canonical counts from R16_GATE2_SCHEMA.md §2.2 (|G_eval^R16|, N_eval P5 baseline pending Gate 5 verification)
    "expected_g_eval_canon": {
        "nadig_hepg2": 1843,
        "nadig_jurkat": 1843,
        "a549_ko": 1837,
        "dld1": 1832,
        "cd4t_rest": 1830,
        "feng_ipsc": 1575,
    },
    "expected_n_eval_r16_min": {
        "nadig_hepg2": 50,
        "nadig_jurkat": 50,
        "a549_ko": 50,
        "dld1": 50,
        "cd4t_rest": 50,
        "feng_ipsc": 50,
    },
    # Exact v5 lock slots (populated when PREREG_I v5 commits exact N_eval counts)
    "locked_n_eval_exact": None,
}


def checkpoint_for_line(target_line: str) -> str:
    if target_line == "nadig_hepg2":
        return "zeroshot/hepg2"
    if target_line == "nadig_jurkat":
        return "zeroshot/jurkat"
    return "zeroshot/ensemble4"


def project_model_deltas_to_eval_genes(
    delta_model: np.ndarray,
    model_genes: Sequence[str],
    eval_genes: Sequence[str],
    hgnc_bridge: Any = None,
) -> np.ndarray:
    """
    Reindex `(N_eval, G_model)` delta predictions cleanly onto exact `eval_genes` (`G_eval^R16(t)`)
    order using frozen R10 `HGNCBridge`.
    """
    canon_pos: Dict[str, int] = {}
    for j, mg in enumerate(model_genes):
        ms = str(mg).strip()
        cm = hgnc_bridge.canonicalize(ms) if hgnc_bridge is not None else ms
        if cm and (cm not in canon_pos or ms == cm):
            canon_pos[cm] = j

    col_idx: List[int] = []
    missing: List[str] = []
    for eg in eval_genes:
        if eg in canon_pos:
            col_idx.append(canon_pos[eg])
        else:
            missing.append(eg)
    if missing:
        raise KeyError(f"eval_genes missing from model output genes: {missing[:5]}")
    return np.asarray(delta_model[:, np.asarray(col_idx, dtype=np.int64)], dtype=np.float64)


def export_line_arm_npz(
    out_path: Path,
    model_key: str,
    checkpoint_id: str,
    target_line: str,
    ref_bundle: Dict[str, Any],
    delta_pred_matched_eval: np.ndarray,
    delta_pred_mis_eval: Dict[str, np.ndarray],
    seen_sources: Sequence[str],
    unseen_sources: Sequence[str],
    neutral_batch_mode: str,
    delta_mode: str,
    hvg_preproc_mode: str,
    arm_metadata: Dict[str, Any],
) -> Path:
    """
    Save per-line, per-arm per-perturbation prediction deltas + observed split-half deltas
    into an uncompressed `.npz` bundle for exact three-way bitwise comparison.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    ordered_mis = sorted(delta_pred_mis_eval.keys())
    arm_names = ["matched"] + [f"mismatched__{s}" for s in ordered_mis]
    seen_set = set(seen_sources)
    arm_strata = ["matched"] + [
        "seen_by_checkpoint" if s in seen_set else "unseen_by_checkpoint" for s in ordered_mis
    ]
    stacked_arms = np.stack(
        [np.asarray(delta_pred_matched_eval, dtype=np.float64)]
        + [np.asarray(delta_pred_mis_eval[s], dtype=np.float64) for s in ordered_mis],
        axis=0,
    )

    payload: Dict[str, Any] = {
        "model_key": np.asarray(str(model_key)),
        "checkpoint_id": np.asarray(str(checkpoint_id)),
        "target_line": np.asarray(str(target_line)),
        "neutral_batch_mode": np.asarray(str(neutral_batch_mode)),
        "delta_mode": np.asarray(str(delta_mode)),
        "hvg_preproc_mode": np.asarray(str(hvg_preproc_mode)),
        "eval_genes": np.asarray(ref_bundle["eval_genes"], dtype=str),
        "eval_perts": np.asarray(ref_bundle["eval_perts"], dtype=str),
        "target_gene_idx": np.asarray(ref_bundle["target_gene_idx"], dtype=np.int64),
        "delta_obs": np.asarray(ref_bundle["delta_obs"], dtype=np.float64),
        "delta_half_a": np.asarray(ref_bundle["delta_half_a"], dtype=np.float64),
        "delta_half_b": np.asarray(ref_bundle["delta_half_b"], dtype=np.float64),
        "half_valid": np.asarray(ref_bundle["half_valid"], dtype=bool),
        "delta_transfer": np.asarray(ref_bundle["delta_transfer"], dtype=np.float64),
        "source_grand_mean": np.asarray(ref_bundle["source_grand_mean"], dtype=np.float64),
        "delta_fct": np.asarray(ref_bundle["delta_fct"], dtype=np.float64),
        "ctrl_mean_libnorm": np.asarray(ref_bundle["ctrl_mean_libnorm"], dtype=np.float64),
        "arm_names": np.asarray(arm_names, dtype=str),
        "arm_strata": np.asarray(arm_strata, dtype=str),
        "seen_sources": np.asarray(list(seen_sources), dtype=str),
        "unseen_sources": np.asarray(list(unseen_sources), dtype=str),
        "delta_pred_arms": stacked_arms,
        "delta_pred_matched": np.asarray(delta_pred_matched_eval, dtype=np.float64),
        "arm_metadata_json": np.asarray(json.dumps(arm_metadata, sort_keys=True)),
    }
    for s in ordered_mis:
        payload[f"delta_pred_mismatched__{s}"] = np.asarray(delta_pred_mis_eval[s], dtype=np.float64)

    np.savez(out_path, **payload)
    return out_path


def evaluate_and_adjudicate_line(
    ref_bundle: Dict[str, Any],
    delta_pred_matched_eval: np.ndarray,
    delta_pred_mis_eval: Dict[str, np.ndarray],
    seen_sources: Sequence[str],
    unseen_sources: Sequence[str],
    n_boot: int = 4000,
    seed: int = DEFAULT_R16_SEED,
    delta: float = DEFAULT_TOST_DELTA,
) -> Dict[str, Any]:
    """Run all `PREREG_I` metrics (`r_orth_cal`, `r_dorth_cal`, `H-I1`, `H-I2`, `M2`) on aligned eval arrays."""
    obs = ref_bundle["delta_obs"]
    ha = ref_bundle["delta_half_a"]
    hb = ref_bundle["delta_half_b"]
    hv = ref_bundle["half_valid"]
    t_mean = ref_bundle["delta_transfer"]
    grand_mean = ref_bundle["source_grand_mean"]
    t_gidx = ref_bundle["target_gene_idx"]
    d_fct = ref_bundle["delta_fct"]

    orth_matched = compute_orth_suite_r16(
        obs, delta_pred_matched_eval, t_mean, t_gidx,
        source_grand_mean=None, delta_half_a=ha, delta_half_b=hb, half_valid=hv,
        mode="all_nontarget", n_boot=n_boot, seed=seed,
    )
    orth_paired = compute_paired_model_vs_all_mismatched_r16(
        obs, delta_pred_matched_eval, delta_pred_mis_eval, t_mean, t_gidx,
        source_grand_mean=None, delta_half_a=ha, delta_half_b=hb, half_valid=hv,
        mode="all_nontarget", n_boot=n_boot, seed=seed,
    )
    h_i1 = adjudicate_claim_i1_line(
        orth_suite=orth_matched,
        paired_mismatched=orth_paired,
        seen_sources=seen_sources,
        unseen_sources=unseen_sources,
        delta=delta,
    )

    dorth_matched = compute_orth_suite_r16(
        obs, delta_pred_matched_eval, t_mean, t_gidx,
        source_grand_mean=grand_mean, delta_half_a=ha, delta_half_b=hb, half_valid=hv,
        mode="all_nontarget", n_boot=n_boot, seed=seed,
    )
    dorth_paired = compute_paired_model_vs_all_mismatched_r16(
        obs, delta_pred_matched_eval, delta_pred_mis_eval, t_mean, t_gidx,
        source_grand_mean=grand_mean, delta_half_a=ha, delta_half_b=hb, half_valid=hv,
        mode="all_nontarget", n_boot=n_boot, seed=seed,
    )
    h_i1_dorth = adjudicate_claim_i1_line(
        orth_suite=dorth_matched,
        paired_mismatched=dorth_paired,
        seen_sources=seen_sources,
        unseen_sources=unseen_sources,
        delta=delta,
    )

    vs_fct_orth = compute_paired_contrast_r16(
        obs, delta_pred_matched_eval, d_fct, t_mean, t_gidx,
        source_grand_mean=None, delta_half_a=ha, delta_half_b=hb, half_valid=hv,
        mode="all_nontarget", n_boot=n_boot, seed=seed,
    )
    h_i2 = adjudicate_claim_i2_line(vs_fct_orth)

    vs_fct_dorth = compute_paired_contrast_r16(
        obs, delta_pred_matched_eval, d_fct, t_mean, t_gidx,
        source_grand_mean=grand_mean, delta_half_a=ha, delta_half_b=hb, half_valid=hv,
        mode="all_nontarget", n_boot=n_boot, seed=seed,
    )

    resp_top50_matched = compute_orth_suite_r16(
        obs, delta_pred_matched_eval, t_mean, t_gidx,
        source_grand_mean=None, delta_half_a=ha, delta_half_b=hb, half_valid=hv,
        mode="resp_top50", top_k=50, n_boot=n_boot, seed=seed,
    )
    m2_weight_shift = compute_winner_curse_weight_shift_r16(
        ctrl_mean_libnorm=ref_bundle["ctrl_mean_libnorm"],
        delta_obs=obs,
        target_gene_idx=t_gidx,
        top_k=50,
        n_boot=n_boot,
        seed=seed,
    )

    return {
        "n_eval": len(ref_bundle["eval_perts"]),
        "g_eval": len(ref_bundle["eval_genes"]),
        "orth_matched": orth_matched,
        "orth_paired": orth_paired,
        "h_i1": h_i1,
        "dorth_matched": dorth_matched,
        "dorth_paired": dorth_paired,
        "h_i1_dorth": h_i1_dorth,
        "vs_fct_orth": vs_fct_orth,
        "vs_fct_dorth": vs_fct_dorth,
        "h_i2": h_i2,
        "resp_top50_matched": resp_top50_matched,
        "m2_weight_shift": m2_weight_shift,
    }


def _check_v5_line_constants(target_line: str, ref_bundle: Dict[str, Any]) -> None:
    g_n = len(ref_bundle["eval_genes"])
    p_n = len(ref_bundle["eval_perts"])
    exp_g = PREREG_I_V5_CONFIG["expected_g_eval_canon"].get(target_line)
    if exp_g is not None and g_n != exp_g:
        raise ValueError(f"{target_line}: |G_eval^R16| = {g_n} != Gate-2 constant {exp_g}")
    locked_exact = PREREG_I_V5_CONFIG["locked_n_eval_exact"]
    if locked_exact is not None and target_line in locked_exact:
        if p_n != int(locked_exact[target_line]):
            raise ValueError(f"{target_line}: N_eval^R16 = {p_n} != PREREG_I v5 locked {locked_exact[target_line]}")
    min_n = PREREG_I_V5_CONFIG["expected_n_eval_r16_min"].get(target_line, 50)
    if p_n < min_n:
        raise ValueError(f"{target_line}: N_eval^R16 = {p_n} < {min_n}")


def run_real_r16(
    models: Sequence[str],
    target_lines: Sequence[str],
    status_dir: Path,
    output_dir: Path,
    models_root: Path = MODELS_ROOT_DEFAULT,
    raw_nadig_dir: Optional[Path] = None,
    device: str = "auto",
    neutral_batch_mode: str = "index_zero",
    hvg_preproc_mode: str = "official_median_log1p",
    n_boot: int = 4000,
    seed: int = DEFAULT_R16_SEED,
    allow_npz_fallback_for_hvg: bool = False,
) -> Dict[str, Any]:
    """
    Execute Gate 5 real-cell evaluation strictly behind `assert_freeze_gates_r16(status_dir)`.
    Never called in synthetic tests.
    """
    t0 = time.time()
    gate_receipt = assert_freeze_gates_r16(status_dir)
    hgnc = get_hgnc_bridge()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Determine all required control-source pools upfront
    needed_ctrl_lines: set[str] = set(target_lines)
    for tl in target_lines:
        ckpt = checkpoint_for_line(tl)
        strata = build_mismatched_strata_r16(tl, ckpt)
        needed_ctrl_lines.update(strata["seen_by_checkpoint"])
        needed_ctrl_lines.update(strata["unseen_by_checkpoint"])

    raw_pools: Dict[str, Any] = {}
    for cl in sorted(needed_ctrl_lines):
        raw_pools[cl] = extract_raw_control_pool_r16(
            cl,
            custom_nadig_dir=raw_nadig_dir,
            n_ctrl=int(PREREG_I_V5_CONFIG["n_ctrl"]),
            cell_set_len=int(PREREG_I_V5_CONFIG["cell_set_len"]),
            seed=seed,
            allow_npz_fallback_for_hvg=allow_npz_fallback_for_hvg,
        )

    se_encoder: Optional[CleanRoomSE600MEncoder] = None
    cached_se_embeddings: Dict[str, Tuple[np.ndarray, Dict[str, Any], Dict[str, Any]]] = {}
    if "STATE_SE" in models:
        se_encoder = CleanRoomSE600MEncoder(
            se600m_dir=models_root / "SE-600M", device=device, batch_size=16, hgnc_bridge=hgnc
        )
        for cl in sorted(needed_ctrl_lines):
            rp = raw_pools[cl]
            if rp.is_log1p:
                raise ValueError(f"STATE_SE requires raw counts for {cl}, got log1p fallback.")
            emb, map_stats = se_encoder.encode_cells(
                ctrl_counts=rp.raw_counts,
                source_genes=rp.gene_names,
                source_ensembl_ids=rp.gene_ids,
                seed=seed,
                hgnc_bridge=hgnc,
            )
            cached_se_embeddings[cl] = (emb, rp.sampling_audit, map_stats)

    summary_manifest: Dict[str, Any] = {
        "gate_receipt": gate_receipt,
        "device": device,
        "neutral_batch_mode": neutral_batch_mode,
        "hvg_preproc_mode": hvg_preproc_mode,
        "seed": int(seed),
        "n_boot": int(n_boot),
        "models": {},
    }

    for model_key in models:
        per_line_hi1: Dict[str, Any] = {}
        per_line_hi2: Dict[str, Any] = {}
        model_entries: Dict[str, Any] = {}

        for tl in target_lines:
            ckpt = checkpoint_for_line(tl)
            ckpt_short = ckpt.split("/")[-1]
            adapter = CleanRoomArcStateAdapter(
                model_key=model_key,
                checkpoint_id=ckpt,
                models_root=models_root,
                status_dir=status_dir,
                device=device,
                neutral_batch_mode=neutral_batch_mode,
                delta_mode="minus_pred_ntc",
                hvg_preproc_mode=hvg_preproc_mode,
                allow_synthetic_smoke=False,
                se600m_encoder=se_encoder,
                hgnc_bridge=hgnc,
            )
            adapter.load_weights()
            ref_bundle_mgenes = next(iter(adapter._bundles.values()))["gene_names"]
            ref_bundle_perts = sorted(k for k in next(iter(adapter._bundles.values()))["pert_lookup"].keys() if k != "non-targeting")

            ref_data = load_line_eval_reference_r16(
                target_line=tl,
                model_genes=ref_bundle_mgenes,
                model_pert_names=ref_bundle_perts,
                hgnc_bridge=hgnc,
            )
            _check_v5_line_constants(tl, ref_data)

            strata = build_mismatched_strata_r16(tl, ckpt)
            mis_sources = list(strata["seen_by_checkpoint"]) + list(strata["unseen_by_checkpoint"])

            matched_rp = raw_pools[tl]
            mis_ctrls = {
                s: (None if model_key == "STATE_SE" else raw_pools[s].raw_counts, raw_pools[s].gene_names, raw_pools[s].batch_labels)
                for s in mis_sources
            }
            umi_map = {
                cl: (raw_pools[cl].cell_total_umis, raw_pools[cl].median_ntc_umi)
                for cl in ([tl] + mis_sources)
            }
            pre_emb = (
                {cl: cached_se_embeddings[cl] for cl in ([tl] + mis_sources)}
                if model_key == "STATE_SE"
                else None
            )

            arm_out = adapter.run_matched_and_mismatched_arms(
                target_line=tl,
                matched_ctrl_matrix=None if model_key == "STATE_SE" else matched_rp.raw_counts,
                matched_genes=matched_rp.gene_names,
                mismatched_ctrls=mis_ctrls,
                eval_perts=ref_data["eval_perts"],
                matched_batch_labels=matched_rp.batch_labels,
                input_is_log1p=bool(matched_rp.is_log1p),
                seed=seed,
                n_ctrl=int(PREREG_I_V5_CONFIG["n_ctrl"]),
                cell_set_len=int(PREREG_I_V5_CONFIG["cell_set_len"]),
                is_synthetic_data=False,
                source_total_umis=umi_map,
                precomputed_embeddings=pre_emb,
                hgnc_bridge=hgnc,
                pre_sampled=True,
                pre_sampled_audits={cl: raw_pools[cl].sampling_audit for cl in ([tl] + mis_sources)},
            )

            delta_m_eval = project_model_deltas_to_eval_genes(
                arm_out["matched"]["delta_pred"],
                arm_out["matched"]["model_genes"],
                ref_data["eval_genes"],
                hgnc_bridge=hgnc,
            )
            delta_mis_eval = {
                s: project_model_deltas_to_eval_genes(
                    arm_out["mismatched"][s]["delta_pred"],
                    arm_out["mismatched"][s]["model_genes"],
                    ref_data["eval_genes"],
                    hgnc_bridge=hgnc,
                )
                for s in mis_sources
            }

            arm_meta = {
                "matched_sampling": {
                    "n_ntc_valid": raw_pools[tl].sampling_audit.get("n_ntc_valid"),
                    "n_sampled": raw_pools[tl].sampling_audit.get("n_sampled"),
                    "sampling_mode": raw_pools[tl].sampling_audit.get("sampling_mode"),
                },
                "gene_mapping": {
                    "matched": arm_out["matched"]["gene_mapping"],
                    **{s: arm_out["mismatched"][s]["gene_mapping"] for s in mis_sources},
                },
                "control_centroid_diagnostics": arm_out["control_centroid_diagnostics"],
            }

            npz_path = output_dir / f"{model_key}_{ckpt_short}__{tl}.npz"
            export_line_arm_npz(
                out_path=npz_path,
                model_key=model_key,
                checkpoint_id=ckpt,
                target_line=tl,
                ref_bundle=ref_data,
                delta_pred_matched_eval=delta_m_eval,
                delta_pred_mis_eval=delta_mis_eval,
                seen_sources=strata["seen_by_checkpoint"],
                unseen_sources=strata["unseen_by_checkpoint"],
                neutral_batch_mode=neutral_batch_mode,
                delta_mode="minus_pred_ntc",
                hvg_preproc_mode=hvg_preproc_mode if model_key == "STATE_HVG" else "se600m_raw_log1p",
                arm_metadata=arm_meta,
            )

            eval_res = evaluate_and_adjudicate_line(
                ref_bundle=ref_data,
                delta_pred_matched_eval=delta_m_eval,
                delta_pred_mis_eval=delta_mis_eval,
                seen_sources=strata["seen_by_checkpoint"],
                unseen_sources=strata["unseen_by_checkpoint"],
                n_boot=n_boot,
                seed=seed,
                delta=float(PREREG_I_V5_CONFIG["tost_delta"]),
            )
            per_line_hi1[tl] = eval_res["h_i1"]
            per_line_hi2[tl] = eval_res["h_i2"]

            metrics_json = output_dir / f"{model_key}_{ckpt_short}__{tl}_metrics.json"
            metrics_json.write_text(
                json.dumps(
                    {
                        "model_key": model_key,
                        "checkpoint_id": ckpt,
                        "target_line": tl,
                        "npz_path": str(npz_path),
                        "arm_metadata": arm_meta,
                        "results": eval_res,
                    },
                    indent=2,
                ),
                encoding="utf-8",
            )
            model_entries[tl] = {"npz_path": str(npz_path), "metrics_json": str(metrics_json), "h_i1_summary": eval_res["h_i1"]["line_summary"], "h_i2_class": eval_res["h_i2"]["classification"]}

        grp = adjudicate_group_r16(per_line_hi1, per_line_hi2)
        summary_manifest["models"][model_key] = {"lines": model_entries, "group_rollup": grp}

    summary_manifest["elapsed_sec"] = round(float(time.time() - t0), 2)
    (output_dir / "manifest_r16_audit.json").write_text(json.dumps(summary_manifest, indent=2), encoding="utf-8")
    return summary_manifest


def run_synthetic_runner_smoke(
    output_dir: Optional[Path] = None,
    device: str = "auto",
    n_sets: int = 4,
    n_perts: int = 10,
    seed: int = DEFAULT_R16_SEED,
    include_real_se600m: bool = False,
) -> Dict[str, Any]:
    """
    End-to-end pure-synthetic runner verification (ZERO real cells opened):
      - Verifies Gate-2 `ST-HVG` (direct `preds`, `16 × 32` `padded=True`, fixed `batch = 0`,
        `HGNCBridge` canonical remaps `AARS1->AARS`, `DARS1->DARS`, `H2AZ1->H2AFZ`, `TAFAZZIN->TAZ`),
        `official_median_log1p` & `libnorm1e4_log1p`, `.npz` export schema, and adjudicator outputs.
      - When `include_real_se600m=True`, also runs real `SE-600M` (`model.safetensors` on `MPS float32`)
        + `ST-SE` (`clip(gene_decoder(preds), 0, 14)`) strictly on synthetic Poisson count matrices.
    """
    hgnc = get_hgnc_bridge() if Path(REPO_ROOT / "round10/src/data_loader_r10.py").exists() else None
    n_ctrl = int(n_sets) * 32
    adapter_hvg = CleanRoomArcStateAdapter(
        model_key="STATE_HVG",
        checkpoint_id="zeroshot/hepg2",
        device=device,
        neutral_batch_mode="index_zero",
        delta_mode="minus_pred_ntc",
        hvg_preproc_mode="official_median_log1p",
        allow_synthetic_smoke=True,
        hgnc_bridge=hgnc,
    )
    adapter_hvg.load_weights()
    g_model = adapter_hvg._bundles["hepg2"]["gene_names"]
    # Build source gene panel replacing legacy symbols with updated canonical symbols to prove HGNCBridge aligns them
    remap_pairs = {"AARS": "AARS1", "DARS": "DARS1", "H2AFZ": "H2AZ1", "TAZ": "TAFAZZIN"}
    src_genes = [remap_pairs.get(g, g) for g in g_model]
    # Include both legacy-renamed canonical targets and standard targets
    canonical_targets = [hgnc.canonicalize(t) if hgnc else t for t in ("AARS", "DARS", "H2AFZ", "TAZ", "AAAS", "AAMP", "AAR2", "AARS2", "AATF", "ABCB7")][:n_perts]

    rng = np.random.default_rng(int(seed))
    G = len(src_genes)
    matched_raw = rng.poisson(2.0, size=(n_ctrl + 64, G)).astype(np.float32)
    matched_b = ["B1"] * ((n_ctrl + 64) // 2) + ["B2"] * (n_ctrl + 64 - (n_ctrl + 64) // 2)
    mis_a_raw = rng.poisson(7.0, size=(n_ctrl + 32, G)).astype(np.float32)
    mis_b_raw = rng.poisson(0.5, size=(n_ctrl + 48, G)).astype(np.float32)
    mis_a_b = ["S1"] * (n_ctrl + 32)
    mis_b_b = ["U1"] * (n_ctrl + 48)

    arm_hvg = adapter_hvg.run_matched_and_mismatched_arms(
        target_line="nadig_hepg2",
        matched_ctrl_matrix=matched_raw,
        matched_genes=src_genes,
        mismatched_ctrls={
            "nadig_jurkat": (mis_a_raw, src_genes, mis_a_b),
            "cd4t_rest": (mis_b_raw, src_genes, mis_b_b),
        },
        eval_perts=canonical_targets,
        matched_batch_labels=matched_b,
        input_is_log1p=False,
        seed=seed,
        n_ctrl=n_ctrl,
        cell_set_len=32,
        is_synthetic_data=True,
        hgnc_bridge=hgnc,
    )

    eval_genes = [hgnc.canonicalize(g) if hgnc else g for g in g_model[:240]]
    d_m_eval = project_model_deltas_to_eval_genes(arm_hvg["matched"]["delta_pred"], g_model, eval_genes, hgnc_bridge=hgnc)
    d_mis_eval = {
        k: project_model_deltas_to_eval_genes(v["delta_pred"], g_model, eval_genes, hgnc_bridge=hgnc)
        for k, v in arm_hvg["mismatched"].items()
    }
    P, Ge = d_m_eval.shape
    t_mean = rng.normal(0.0, 0.15, size=(P, Ge))
    grand_mean = rng.normal(0.0, 0.25, size=Ge)
    obs = d_m_eval + rng.normal(0.0, 0.08, size=(P, Ge))
    ha = obs + rng.normal(0.0, 0.04, size=(P, Ge))
    hb = obs + rng.normal(0.0, 0.04, size=(P, Ge))
    ref_syn = {
        "eval_genes": eval_genes,
        "eval_perts": canonical_targets,
        "target_gene_idx": np.asarray([eval_genes.index(p) if p in eval_genes else -1 for p in canonical_targets], dtype=np.int64),
        "delta_obs": obs,
        "delta_half_a": ha,
        "delta_half_b": hb,
        "half_valid": np.ones(P, dtype=bool),
        "delta_transfer": t_mean,
        "source_grand_mean": grand_mean,
        "delta_fct": t_mean + 0.25 * (obs - t_mean),
        "ctrl_mean_libnorm": np.linspace(0.1, 5.0, Ge),
    }

    npz_written: Optional[Path] = None
    if output_dir is not None:
        npz_written = export_line_arm_npz(
            out_path=Path(output_dir) / "SYNTHETIC_STATE_HVG_hepg2__nadig_hepg2.npz",
            model_key="STATE_HVG",
            checkpoint_id="zeroshot/hepg2",
            target_line="nadig_hepg2",
            ref_bundle=ref_syn,
            delta_pred_matched_eval=d_m_eval,
            delta_pred_mis_eval=d_mis_eval,
            seen_sources=("nadig_jurkat",),
            unseen_sources=("cd4t_rest",),
            neutral_batch_mode="index_zero",
            delta_mode="minus_pred_ntc",
            hvg_preproc_mode="official_median_log1p",
            arm_metadata={"synthetic": True},
        )
    ev_hvg = evaluate_and_adjudicate_line(
        ref_syn, d_m_eval, d_mis_eval, seen_sources=("nadig_jurkat",), unseen_sources=("cd4t_rest",),
        n_boot=150, seed=seed,
    )

    out: Dict[str, Any] = {
        "st_hvg_smoke": {
            "delta_matched_shape": d_m_eval.shape,
            "gene_mapping": arm_hvg["matched"]["gene_mapping"],
            "h_i1_summary": ev_hvg["h_i1"]["line_summary"],
            "npz_path": str(npz_written) if npz_written else None,
        }
    }

    if include_real_se600m:
        se_enc = CleanRoomSE600MEncoder(device=device, batch_size=16, hgnc_bridge=hgnc)
        se_enc.load()
        pe_sample_genes = se_enc._pe_keys[:2500]
        raw_se_m = rng.poisson(3.0, size=(n_ctrl, len(pe_sample_genes))).astype(np.float32)
        raw_se_s = rng.poisson(9.0, size=(n_ctrl, len(pe_sample_genes))).astype(np.float32)
        adapter_se = CleanRoomArcStateAdapter(
            model_key="STATE_SE",
            checkpoint_id="zeroshot/hepg2",
            device=device,
            neutral_batch_mode="index_zero",
            delta_mode="minus_pred_ntc",
            allow_synthetic_smoke=True,
            se600m_encoder=se_enc,
            hgnc_bridge=hgnc,
        )
        arm_se = adapter_se.run_matched_and_mismatched_arms(
            target_line="nadig_hepg2",
            matched_ctrl_matrix=raw_se_m,
            matched_genes=pe_sample_genes,
            mismatched_ctrls={"cd4t_rest": (raw_se_s, pe_sample_genes, ["U1"] * n_ctrl)},
            eval_perts=canonical_targets,
            matched_batch_labels=["B1"] * n_ctrl,
            input_is_log1p=False,
            seed=seed,
            n_ctrl=n_ctrl,
            cell_set_len=32,
            is_synthetic_data=True,
            hgnc_bridge=hgnc,
        )
        d_se_m = arm_se["matched"]["delta_pred"]
        d_se_s = arm_se["mismatched"]["cd4t_rest"]["delta_pred"]
        out["se600m_st_se_smoke"] = {
            "device": se_enc.device_str,
            "delta_matched_shape": d_se_m.shape,
            "l2_ctrl_swap": float(np.linalg.norm(d_se_m - d_se_s)),
            "mapping": arm_se["matched"]["gene_mapping"],
        }

    return out


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Clean-Room Round 16 Real-Data & Synthetic Runner (A-T3).")
    p.add_argument("--phase", choices=["real_r16", "synthetic_smoke"], default="real_r16")
    p.add_argument("--models", nargs="+", default=["STATE_SE", "STATE_HVG"], choices=["STATE_SE", "STATE_HVG"])
    p.add_argument(
        "--lines",
        nargs="+",
        default=["nadig_hepg2", "nadig_jurkat", "a549_ko", "dld1", "cd4t_rest"],
    )
    p.add_argument("--status-dir", type=Path, default=STATUS_DIR_DEFAULT)
    p.add_argument("--output-dir", type=Path, default=OUTPUT_DIR_DEFAULT)
    p.add_argument("--models-root", type=Path, default=MODELS_ROOT_DEFAULT)
    p.add_argument("--raw-nadig-dir", type=Path, default=None)
    p.add_argument("--device", choices=["auto", "mps", "cpu"], default="auto")
    p.add_argument("--neutral-batch-mode", default="index_zero", choices=["index_zero", "bypass_zero", "mean_batch_emb", "multi_batch_ens8"])
    p.add_argument("--hvg-preproc", default="official_median_log1p", choices=["official_median_log1p", "libnorm1e4_log1p"])
    p.add_argument("--n-boot", type=int, default=4000)
    p.add_argument("--seed", type=int, default=DEFAULT_R16_SEED)
    p.add_argument("--include-real-se600m-synthetic", action="store_true")
    p.add_argument("--allow-npz-fallback-for-hvg", action="store_true")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = parse_args(argv)
    if args.phase == "real_r16":
        run_real_r16(
            models=args.models,
            target_lines=args.lines,
            status_dir=args.status_dir,
            output_dir=args.output_dir,
            models_root=args.models_root,
            raw_nadig_dir=args.raw_nadig_dir,
            device=args.device,
            neutral_batch_mode=args.neutral_batch_mode,
            hvg_preproc_mode=args.hvg_preproc,
            n_boot=args.n_boot,
            seed=args.seed,
            allow_npz_fallback_for_hvg=args.allow_npz_fallback_for_hvg,
        )
    else:
        res = run_synthetic_runner_smoke(
            output_dir=args.output_dir,
            device=args.device,
            include_real_se600m=args.include_real_se600m,
            seed=args.seed,
        )
        print(json.dumps(res, indent=2, default=str))


if __name__ == "__main__":
    main()
