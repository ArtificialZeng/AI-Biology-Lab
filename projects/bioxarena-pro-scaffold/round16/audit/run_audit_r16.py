#!/usr/bin/env python3
"""
Clean-Room Round 16 Real-Data & Synthetic-Parity Evaluation Runner
(`A-T3`, strictly aligned to `PREREG_I.md v5-FINAL-PRE-LOCK` & `Addendum A`).

Enforces strict runtime blinding & schema hard gates BEFORE accessing real data:
  1. `PREREG_I_COMMIT.txt` exists and is non-empty.
  2. `CODE_FREEZE_R16.txt` exists and is non-empty.
  3. `R16_GO.txt` exists, is non-empty, and postdates `R16_GATE2_SCHEMA.md` when present.

Aligned Structural & Addendum A Rules:
  - `SE-600M`: exact gene-symbol string lookup (`use_hgnc_rescue=False`, `v5 §2.1.4 L81`), `MPS float32`.
  - `ST-HVG`: `official_median_log1p` primary (`sc.pp.normalize_total(target_sum=None)` computed strictly
    on each sampled control pool's full measured transcriptome -> `log1p`), `libnorm1e4_log1p` sensitivity,
    and `Perm100` null (`PREREG_E §4.1` & `Addendum A §4.2`: `np.random.default_rng(k).permutation(G)`
    on matched-arm predicted delta gene columns in `G_eval^R16`).
  - `CellularAI`: explicitly excluded per Addendum A (`cellularai_included=False`).
  - `EXACT_INVARIANCE` (`Addendum A Item 10`): `max_abs_diff <= 1e-12` treated as `EQUIVALENT` at line level.
  - Same-line non-overlapping `resample` arm (`seed=20261602`, `Addendum A Items 8–9`), with sampled
    `cell_barcodes` exported per control source for exact 1-to-1 cross-party verification.
  - Group 3 (`gse295214_h1`): separate reporting, primary `FCT` = `r9_6 \\ {vcc2025_h1}`, secondary
    `full r9_6` strictly evaluated on identical locked primary targets & genes (`Addendum A`).
  - `G4` 2×2 decomposition (`Group 1`: positional alignment `B_M/B_S`; `Group 2/3`: fixed `B_M=0, B_S=17`),
    neutral-batch robustness panel, `nadig_hepg2` `adata_pred.h5ad` check (`median>=0.90, pooled>=0.95`),
    and cross-party array/barcode verifier (`barcodes exact`, `pred<=1e-4`, `recompute<=1e-10`,
    `cross_party<=1e-5`, `verdict 100%`).

Strictly NEVER reads or imports `round16/src/`.
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
    compute_gene_label_perm_null_r16,
    compute_orth_suite_r16,
    compute_paired_contrast_r16,
    compute_paired_diff_vs_resample_r16,
    compute_paired_model_vs_all_mismatched_r16,
    compute_practical_metrics_r16,
    compute_winner_curse_weight_shift_r16,
    summarize_perm100_r16,
)
from audit_r16_raw_controls import (  # noqa: E402
    extract_raw_control_pool_r16,
    load_line_eval_reference_r16,
)
from audit_r16_sampling import (  # noqa: E402
    LOCAL_RAW_CONTROL_SOURCES,
    RESAMPLE_R16_SEED,
    build_mismatched_strata_r16,
)
from audit_r16_state_adapter import (  # noqa: E402
    CleanRoomArcStateAdapter,
    CleanRoomSE600MEncoder,
    MODELS_ROOT_DEFAULT,
    assert_freeze_gates_r16,
    get_hgnc_bridge,
    map_batch_labels_to_indices,
)

STATUS_DIR_DEFAULT = REPO_ROOT / "round16" / "status"
OUTPUT_DIR_DEFAULT = AUDIT_DIR / "results"
MAIN_ARRAYS_DIR_DEFAULT = REPO_ROOT / "round16" / "results" / "arrays"


PREREG_I_V5_CONFIG: Dict[str, Any] = {
    "prereg_version": "v5-FINAL-PRE-LOCK+AddendumA",
    "prereg_commit": "c1fc7c6d24dd128bc5e2d969d8c78224a2cbc118",
    "prereg_sha256": "43d242e460fd64db2e4af1e628e886cf1d94fb4ccd1e95925928a99cb96ee4a8",
    "n_ctrl": 512,
    "cell_set_len": 32,
    "seed": DEFAULT_R16_SEED,
    "resample_seed": RESAMPLE_R16_SEED,
    "n_boot": 4000,
    "tost_delta": DEFAULT_TOST_DELTA,
    "exact_invariance_max_abs_tol": 1e-12,
    "neutral_batch_mode_primary": "index_zero",  # fixed batch = 0
    "delta_mode_primary": "minus_pred_ntc",      # pred(p) - pred("non-targeting")
    "hvg_preproc_primary": "official_median_log1p",
    "se600m_gene_lookup": "exact_string_lookup",
    "cellularai_included": False,                # Addendum A: CellularAI excluded
    "adata_pred_parity_line": "nadig_hepg2",     # Addendum A: nadig_hepg2 only
    "adata_pred_parity_thresholds": {
        "median_per_pert_pearson": 0.90,
        "pooled_pearson": 0.95,
    },
    "g4_pairing": {
        "group1": "positional_0_to_511_native_gem_group",
        "group2_and_group3": {"b_m": 0, "b_s": 17, "spread_indices": (0, 1, 17, 35)},
    },
    "cross_party_tolerances": {
        "barcodes_exact_match": True,
        "se600m_x_state_abs_diff": 1e-5,
        "pred_abs_diff": 1e-4,
        "same_array_recompute_metric_diff": 1e-10,
        "cross_party_metric_diff": 1e-5,
        "verdict_agreement": 1.0,
    },
    "se600m_collator": "official_FilteredGenesCounts+VCIDatasetSentenceCollator(training=False)",
    "se600m_batch_size": 16,
    "group1_lines": ("nadig_hepg2", "nadig_jurkat"),
    "group2_crispri_lines": ("dld1", "cd4t_rest"),
    "group2_crisprko_lines": ("a549_ko",),
    "group3_prospective_lines": ("gse295214_h1",),
    "control_pool_lines": LOCAL_RAW_CONTROL_SOURCES,
    "expected_g_eval_canon": {
        "nadig_hepg2": 1843,
        "nadig_jurkat": 1843,
        "a549_ko": 1837,
        "dld1": 1832,
        "cd4t_rest": 1830,
        "feng_ipsc": 1575,
        "cd4t_stim48": 1836,
        "a549_ko_cas9pos": 1837,
        "a549_myllia": 1836,
        "hek293_xucao": 1838,
    },
    "expected_n_eval_r16_min": {
        "nadig_hepg2": 50,
        "nadig_jurkat": 50,
        "a549_ko": 50,
        "dld1": 50,
        "cd4t_rest": 50,
        "feng_ipsc": 50,
        "gse295214_h1": 50,
    },
    "locked_n_eval_exact": {
        "nadig_hepg2": 740,
        "nadig_jurkat": 873,
        "a549_ko": 403,
        "dld1": 75,
        "cd4t_rest": 245,
        "feng_ipsc": 88,
        "cd4t_stim48": 401,
        "a549_ko_cas9pos": 403,
        "a549_myllia": 20,
        "hek293_xucao": 16,
    },
    "locked_half_valid_exact": {
        "nadig_hepg2": 740,
        "nadig_jurkat": 873,
        "a549_ko": 403,
        "dld1": 75,
        "cd4t_rest": 230,
        "feng_ipsc": 88,
        "cd4t_stim48": 333,
        "a549_ko_cas9pos": 403,
        "a549_myllia": 20,
        "hek293_xucao": 16,
    },
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


def _json_default_r16(obj: Any) -> Any:
    """Serialize NumPy arrays/scalars, Paths, and sets/tuples deterministically without altering values."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, (set, tuple)):
        if all(isinstance(x, str) for x in obj):
            return sorted(obj)
        return list(obj)
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


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
    matched_barcodes: Optional[Sequence[str]] = None,
    mismatched_barcodes: Optional[Dict[str, Sequence[str]]] = None,
    delta_pred_resample_eval: Optional[np.ndarray] = None,
    resample_barcodes: Optional[Sequence[str]] = None,
) -> Path:
    """
    Save per-line, per-arm per-perturbation prediction deltas, observed split-half deltas,
    and exact sampled control cell barcodes into an uncompressed `.npz` bundle.
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
        "arm_metadata_json": np.asarray(json.dumps(arm_metadata, sort_keys=True, default=_json_default_r16)),
    }
    if matched_barcodes is not None:
        arr_m_bc = np.asarray(list(matched_barcodes), dtype=str)
        payload["matched_barcodes"] = arr_m_bc
        payload["matched_ctrl_barcodes"] = arr_m_bc
    if delta_pred_resample_eval is not None:
        payload["delta_pred_resample"] = np.asarray(delta_pred_resample_eval, dtype=np.float64)
    if resample_barcodes is not None:
        arr_r_bc = np.asarray(list(resample_barcodes), dtype=str)
        payload["resample_barcodes"] = arr_r_bc
        payload["resample_ctrl_barcodes"] = arr_r_bc
    mis_bcs = mismatched_barcodes or {}
    for s in ordered_mis:
        payload[f"delta_pred_mismatched__{s}"] = np.asarray(delta_pred_mis_eval[s], dtype=np.float64)
        if s in mis_bcs and mis_bcs[s] is not None:
            payload[f"mismatched_barcodes__{s}"] = np.asarray(list(mis_bcs[s]), dtype=str)

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
    delta_pred_resample_eval: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """Run all `PREREG_I v5` & `Addendum A` metrics + mechanical adjudication on aligned eval arrays."""
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
    orth_pos_ci975 = bool(
        np.isfinite(orth_matched.get("valid_frac", 0.0))
        and orth_matched["valid_frac"] >= 0.95
        and np.isfinite(orth_matched.get("rho_ab_orth", float("nan")))
        and orth_matched["rho_ab_orth"] > 0.0
        and np.isfinite(orth_matched.get("ci975_lo", float("nan")))
        and orth_matched["ci975_lo"] > 0.0
    )
    resample_orth = (
        compute_paired_diff_vs_resample_r16(
            obs, delta_pred_matched_eval, delta_pred_resample_eval, delta_pred_mis_eval,
            t_mean, t_gidx, source_grand_mean=None, delta_half_a=ha, delta_half_b=hb,
            half_valid=hv, mode="all_nontarget", n_boot=n_boot, seed=seed,
            matched_positive_ci975=orth_pos_ci975,
        )
        if delta_pred_resample_eval is not None
        else None
    )
    h_i1 = adjudicate_claim_i1_line(
        orth_suite=orth_matched,
        paired_mismatched=orth_paired,
        seen_sources=seen_sources,
        unseen_sources=unseen_sources,
        delta=delta,
        resample_comparison=resample_orth,
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
    dorth_pos_ci975 = bool(
        np.isfinite(dorth_matched.get("valid_frac", 0.0))
        and dorth_matched["valid_frac"] >= 0.95
        and np.isfinite(dorth_matched.get("rho_ab_orth", float("nan")))
        and dorth_matched["rho_ab_orth"] > 0.0
        and np.isfinite(dorth_matched.get("ci975_lo", float("nan")))
        and dorth_matched["ci975_lo"] > 0.0
    )
    resample_dorth = (
        compute_paired_diff_vs_resample_r16(
            obs, delta_pred_matched_eval, delta_pred_resample_eval, delta_pred_mis_eval,
            t_mean, t_gidx, source_grand_mean=grand_mean, delta_half_a=ha, delta_half_b=hb,
            half_valid=hv, mode="all_nontarget", n_boot=n_boot, seed=seed,
            matched_positive_ci975=dorth_pos_ci975,
        )
        if delta_pred_resample_eval is not None
        else None
    )
    h_i1_dorth = adjudicate_claim_i1_line(
        orth_suite=dorth_matched,
        paired_mismatched=dorth_paired,
        seen_sources=seen_sources,
        unseen_sources=unseen_sources,
        delta=delta,
        resample_comparison=resample_dorth,
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
    practical_matched = compute_practical_metrics_r16(obs, delta_pred_matched_eval, t_gidx, top_k=50)
    practical_fct = compute_practical_metrics_r16(obs, d_fct, t_gidx, top_k=50)

    out: Dict[str, Any] = {
        "n_eval": len(ref_bundle["eval_perts"]),
        "g_eval": len(ref_bundle["eval_genes"]),
        "orth_matched": orth_matched,
        "orth_paired": orth_paired,
        "resample_orth": resample_orth,
        "h_i1": h_i1,
        "dorth_matched": dorth_matched,
        "dorth_paired": dorth_paired,
        "resample_dorth": resample_dorth,
        "h_i1_dorth": h_i1_dorth,
        "vs_fct_orth": vs_fct_orth,
        "vs_fct_dorth": vs_fct_dorth,
        "h_i2": h_i2,
        "resp_top50_matched": resp_top50_matched,
        "m2_weight_shift": m2_weight_shift,
        "practical_metrics_matched": practical_matched,
        "practical_metrics_fct": practical_fct,
    }

    # Addendum A §1: Group 3 secondary `full r9_6` strictly evaluated on identical primary targets & genes
    if "secondary_full_r9_6_same_targets" in ref_bundle:
        sens = ref_bundle["secondary_full_r9_6_same_targets"]
        st_orth = compute_orth_suite_r16(
            obs, delta_pred_matched_eval, sens["delta_transfer"], t_gidx,
            source_grand_mean=None, delta_half_a=ha, delta_half_b=hb, half_valid=hv,
            mode="all_nontarget", n_boot=n_boot, seed=seed,
        )
        fct_orth = compute_orth_suite_r16(
            obs, sens["delta_fct"], sens["delta_transfer"], t_gidx,
            source_grand_mean=None, delta_half_a=ha, delta_half_b=hb, half_valid=hv,
            mode="all_nontarget", n_boot=n_boot, seed=seed,
        )
        vs_fct_full = compute_paired_contrast_r16(
            obs, delta_pred_matched_eval, sens["delta_fct"], sens["delta_transfer"], t_gidx,
            source_grand_mean=None, delta_half_a=ha, delta_half_b=hb, half_valid=hv,
            mode="all_nontarget", n_boot=n_boot, seed=seed,
        )
        st_dorth = compute_orth_suite_r16(
            obs, delta_pred_matched_eval, sens["delta_transfer"], t_gidx,
            source_grand_mean=sens["source_grand_mean"], delta_half_a=ha, delta_half_b=hb, half_valid=hv,
            mode="all_nontarget", n_boot=n_boot, seed=seed,
        )
        fct_dorth = compute_orth_suite_r16(
            obs, sens["delta_fct"], sens["delta_transfer"], t_gidx,
            source_grand_mean=sens["source_grand_mean"], delta_half_a=ha, delta_half_b=hb, half_valid=hv,
            mode="all_nontarget", n_boot=n_boot, seed=seed,
        )
        vs_fct_full_dorth = compute_paired_contrast_r16(
            obs, delta_pred_matched_eval, sens["delta_fct"], sens["delta_transfer"], t_gidx,
            source_grand_mean=sens["source_grand_mean"], delta_half_a=ha, delta_half_b=hb, half_valid=hv,
            mode="all_nontarget", n_boot=n_boot, seed=seed,
        )
        out["secondary_full_r9_6_same_targets"] = {
            "source_lines_fct": sens["source_lines_fct"],
            "n_eval": len(sens["eval_perts"]),
            "state_orth": st_orth,
            "fct_orth": fct_orth,
            "vs_fct_orth": vs_fct_full,
            "state_dorth": st_dorth,
            "fct_dorth": fct_dorth,
            "vs_fct_dorth": vs_fct_full_dorth,
            "h_i2": adjudicate_claim_i2_line(vs_fct_full),
            "h_i2_dorth": adjudicate_claim_i2_line(vs_fct_full_dorth),
        }
    return out


def verify_nadig_hepg2_adata_pred_parity(
    adapter_delta_mbm: np.ndarray,
    adata_pred_delta: np.ndarray,
    target_gene_idx: Sequence[int],
    min_median_per_pert_pearson: float = 0.90,
    min_pooled_pearson: float = 0.95,
) -> Dict[str, Any]:
    """
    Addendum A (`适配器核对`): compare `nadig_hepg2` native-batch `(M, B_M)` adapter deltas
    against archived `adata_pred.h5ad` deltas, requiring:
      - median per-perturbation Pearson >= 0.90
      - pooled Pearson >= 0.95
    """
    pm = compute_practical_metrics_r16(
        delta_obs=adata_pred_delta,
        delta_pred=adapter_delta_mbm,
        target_gene_idx=target_gene_idx,
        top_k=50,
    )
    med_r = float(pm["per_pert_pearson_median"])
    pool_r = float(pm["raw_pearson_delta"])
    passed = bool(
        np.isfinite(med_r)
        and np.isfinite(pool_r)
        and med_r >= float(min_median_per_pert_pearson)
        and pool_r >= float(min_pooled_pearson)
    )
    return {
        "line": "nadig_hepg2",
        "median_per_pert_pearson": med_r,
        "pooled_pearson": pool_r,
        "cosine_similarity": float(pm["cosine_similarity"]),
        "threshold_median": float(min_median_per_pert_pearson),
        "threshold_pooled": float(min_pooled_pearson),
        "passed": passed,
    }


def _collect_nested_numeric_leaves(d: Any, prefix: str = "") -> Dict[str, float]:
    out: Dict[str, float] = {}
    if isinstance(d, dict):
        for k, v in d.items():
            p = f"{prefix}.{k}" if prefix else str(k)
            out.update(_collect_nested_numeric_leaves(v, p))
    elif isinstance(d, (int, float, np.integer, np.floating)) and not isinstance(d, bool):
        fv = float(d)
        if np.isfinite(fv):
            out[prefix] = fv
    return out


def verify_cross_party_arrays_and_barcodes(
    audit_npz_path: Path,
    main_pred_arrays: Dict[str, np.ndarray],
    main_barcodes: Optional[Dict[str, Sequence[str]]] = None,
    main_metrics_summary: Optional[Dict[str, Any]] = None,
    main_verdict: Optional[Dict[str, Any]] = None,
    n_boot: int = 4000,
    seed: int = DEFAULT_R16_SEED,
    delta: float = DEFAULT_TOST_DELTA,
    pred_tol: float = 1e-4,
    same_array_tol: float = 1e-10,
    cross_party_tol: float = 1e-5,
    audit_se_npz_path: Optional[Path] = None,
    main_se_embeddings: Optional[Dict[str, np.ndarray]] = None,
    se_emb_tol: float = 1e-5,
) -> Dict[str, Any]:
    """
    Addendum A strict cross-party parity & tolerance verifier:
      - Control barcodes per source match 100% identically
      - Pre-ST-SE `2058`-dim `X_state` embeddings match `<= 1e-5` (`Addendum A §11`)
      - Max absolute prediction delta difference `<= 1e-4`
      - Same-array recomputed metric difference `<= 1e-10`
      - Cross-party metric difference `<= 1e-5`
      - `100%` verdict parity (`H-I1`, `H-I2`, and pair categories)
    """
    z = np.load(audit_npz_path, allow_pickle=False)
    ref_bundle = {
        "eval_genes": z["eval_genes"].tolist(),
        "eval_perts": z["eval_perts"].tolist(),
        "target_gene_idx": z["target_gene_idx"],
        "delta_obs": z["delta_obs"],
        "delta_half_a": z["delta_half_a"],
        "delta_half_b": z["delta_half_b"],
        "half_valid": z["half_valid"],
        "delta_transfer": z["delta_transfer"],
        "source_grand_mean": z["source_grand_mean"],
        "delta_fct": z["delta_fct"],
        "ctrl_mean_libnorm": z["ctrl_mean_libnorm"],
    }
    seen_s = z["seen_sources"].tolist()
    unseen_s = z["unseen_sources"].tolist()
    audit_matched = np.asarray(z["delta_pred_matched"], dtype=np.float64)
    audit_mis = {str(s): np.asarray(z[f"delta_pred_mismatched__{s}"], dtype=np.float64) for s in (seen_s + unseen_s)}
    audit_resamp = np.asarray(z["delta_pred_resample"], dtype=np.float64) if "delta_pred_resample" in z else None

    # 1. Exact barcode verification
    bc_mismatches: List[str] = []
    if main_barcodes is not None:
        if "matched" in main_barcodes and "matched_barcodes" in z:
            if list(main_barcodes["matched"]) != z["matched_barcodes"].tolist():
                bc_mismatches.append("matched")
        if "resample" in main_barcodes and "resample_barcodes" in z:
            if list(main_barcodes["resample"]) != z["resample_barcodes"].tolist():
                bc_mismatches.append("resample")
        for s in (seen_s + unseen_s):
            k_z = f"mismatched_barcodes__{s}"
            if s in main_barcodes and k_z in z:
                if list(main_barcodes[s]) != z[k_z].tolist():
                    bc_mismatches.append(s)
    barcodes_ok = len(bc_mismatches) == 0

    # 2. Prediction array tolerance <= 1e-4
    pred_diffs: Dict[str, float] = {}
    pred_diffs["matched"] = float(np.max(np.abs(audit_matched - np.asarray(main_pred_arrays["matched"], dtype=np.float64))))
    for s in (seen_s + unseen_s):
        key = s if s in main_pred_arrays else f"mismatched__{s}"
        pred_diffs[s] = float(np.max(np.abs(audit_mis[s] - np.asarray(main_pred_arrays[key], dtype=np.float64))))
    if audit_resamp is not None and "resample" in main_pred_arrays:
        pred_diffs["resample"] = float(np.max(np.abs(audit_resamp - np.asarray(main_pred_arrays["resample"], dtype=np.float64))))
    max_pred_diff = max(pred_diffs.values()) if pred_diffs else 0.0
    pred_ok = max_pred_diff <= float(pred_tol)

    # 3. Same-array metric recomputation consistency <= 1e-10 (Audit vs Recomputed on Main Arrays)
    m_main_matched = np.asarray(main_pred_arrays["matched"], dtype=np.float64)
    m_main_mis = {
        s: np.asarray(main_pred_arrays[s if s in main_pred_arrays else f"mismatched__{s}"], dtype=np.float64)
        for s in (seen_s + unseen_s)
    }
    m_main_resamp = np.asarray(main_pred_arrays["resample"], dtype=np.float64) if "resample" in main_pred_arrays else None

    eval_on_main_pass1 = evaluate_and_adjudicate_line(
        ref_bundle, m_main_matched, m_main_mis, seen_s, unseen_s,
        n_boot=n_boot, seed=seed, delta=delta, delta_pred_resample_eval=m_main_resamp,
    )
    eval_on_main_pass2 = evaluate_and_adjudicate_line(
        ref_bundle, m_main_matched, m_main_mis, seen_s, unseen_s,
        n_boot=n_boot, seed=seed, delta=delta, delta_pred_resample_eval=m_main_resamp,
    )
    leaves_1 = _collect_nested_numeric_leaves(eval_on_main_pass1)
    leaves_2 = _collect_nested_numeric_leaves(eval_on_main_pass2)
    same_arr_diff = max((abs(leaves_1[k] - leaves_2[k]) for k in leaves_1 if k in leaves_2), default=0.0)
    if main_metrics_summary is not None:
        ext_leaves = _collect_nested_numeric_leaves(main_metrics_summary)
        same_arr_diff = max(same_arr_diff, max((abs(leaves_1[k] - ext_leaves[k]) for k in ext_leaves if k in leaves_1), default=0.0))
    same_array_ok = same_arr_diff <= float(same_array_tol)

    # 4. Cross-party metric difference <= 1e-5 & 100% verdict parity (Audit Arrays vs Main Arrays)
    eval_on_audit = evaluate_and_adjudicate_line(
        ref_bundle, audit_matched, audit_mis, seen_s, unseen_s,
        n_boot=n_boot, seed=seed, delta=delta, delta_pred_resample_eval=audit_resamp,
    )
    leaves_a = _collect_nested_numeric_leaves(eval_on_audit)
    cross_metric_diff = max(
        (abs(leaves_a[k] - leaves_1[k]) for k in ("orth_matched.estimate", "orth_matched.ci975_lo", "orth_matched.ci975_hi", "vs_fct_orth.estimate", "h_i1.matched_rho_ab_orth") if k in leaves_a and k in leaves_1),
        default=0.0,
    )
    for k in leaves_a:
        if k.startswith("orth_paired.per_source.") or k.startswith("h_i1.per_source."):
            if k in leaves_1:
                cross_metric_diff = max(cross_metric_diff, abs(leaves_a[k] - leaves_1[k]))
    cross_ok = cross_metric_diff <= float(cross_party_tol)

    verdict_ok = (
        eval_on_audit["h_i1"]["line_summary"] == eval_on_main_pass1["h_i1"]["line_summary"]
        and eval_on_audit["h_i2"]["classification"] == eval_on_main_pass1["h_i2"]["classification"]
        and all(
            eval_on_audit["h_i1"]["per_source"][s]["category"] == eval_on_main_pass1["h_i1"]["per_source"][s]["category"]
            for s in (seen_s + unseen_s)
        )
    )
    if main_verdict is not None:
        if "h_i1_summary" in main_verdict:
            verdict_ok &= bool(eval_on_audit["h_i1"]["line_summary"] == main_verdict["h_i1_summary"])
        if "h_i2_class" in main_verdict:
            verdict_ok &= bool(eval_on_audit["h_i2"]["classification"] == main_verdict["h_i2_class"])

    # 5. Pre-ST-SE 2058-dim X_state embedding tolerance <= 1e-5 (Addendum A §11)
    se_emb_max_diff: Optional[float] = None
    se_emb_ok: Optional[bool] = None
    if audit_se_npz_path is not None and main_se_embeddings is not None:
        with np.load(audit_se_npz_path, allow_pickle=False) as sz:
            se_diffs: List[float] = []
            for pool_k, m_arr in sorted(main_se_embeddings.items()):
                z_key = f"{pool_k}__X_state" if f"{pool_k}__X_state" in sz else pool_k
                if z_key in sz:
                    diff_k = float(
                        np.max(
                            np.abs(
                                np.asarray(sz[z_key], dtype=np.float32)
                                - np.asarray(m_arr, dtype=np.float32)
                            )
                        )
                    )
                    se_diffs.append(diff_k)
            se_emb_max_diff = max(se_diffs) if se_diffs else 0.0
            se_emb_ok = bool(se_emb_max_diff <= float(se_emb_tol))

    all_ok = bool(barcodes_ok and pred_ok and same_array_ok and cross_ok and verdict_ok)
    if se_emb_ok is not None:
        all_ok &= se_emb_ok

    return {
        "all_passed": all_ok,
        "barcodes_exact_match": barcodes_ok,
        "barcode_mismatch_arms": bc_mismatches,
        "se_emb_max_abs_diff": se_emb_max_diff,
        "se_emb_within_tol": se_emb_ok,
        "max_pred_abs_diff": max_pred_diff,
        "pred_within_tol": pred_ok,
        "same_array_metric_max_diff": same_arr_diff,
        "same_array_within_tol": same_array_ok,
        "cross_party_metric_max_diff": cross_metric_diff,
        "cross_party_within_tol": cross_ok,
        "verdict_100pct_match": verdict_ok,
    }


def _check_v5_line_constants(target_line: str, ref_bundle: Dict[str, Any]) -> None:
    g_n = len(ref_bundle["eval_genes"])
    p_n = len(ref_bundle["eval_perts"])
    hv_n = int(np.sum(ref_bundle["half_valid"]))
    exp_g = PREREG_I_V5_CONFIG["expected_g_eval_canon"].get(target_line)
    if exp_g is not None and g_n != int(exp_g):
        raise ValueError(f"{target_line}: |G_eval^R16| = {g_n} != PREREG_I v5 locked {exp_g}")
    locked_exact = PREREG_I_V5_CONFIG["locked_n_eval_exact"]
    if target_line in locked_exact and p_n != int(locked_exact[target_line]):
        raise ValueError(f"{target_line}: N_eval^R16 = {p_n} != PREREG_I v5 locked {locked_exact[target_line]}")
    locked_hv = PREREG_I_V5_CONFIG["locked_half_valid_exact"]
    if target_line in locked_hv and hv_n != int(locked_hv[target_line]):
        raise ValueError(f"{target_line}: half_valid = {hv_n} != PREREG_I v5 locked {locked_hv[target_line]}")
    min_n = PREREG_I_V5_CONFIG["expected_n_eval_r16_min"].get(target_line, 50)
    if target_line in PREREG_I_V5_CONFIG["expected_n_eval_r16_min"] and p_n < min_n:
        raise ValueError(f"{target_line}: N_eval^R16 = {p_n} < {min_n}")


def run_g4_factorial_for_line(
    target_line: str,
    ref_bundle: Dict[str, Any],
    raw_pools: Dict[str, Any],
    se_embeddings: Dict[str, Tuple[np.ndarray, Dict[str, Any], Dict[str, Any]]],
    models_root: Path,
    status_dir: Path,
    device: str,
    n_boot: int,
    seed: int,
    hgnc: Any,
    is_synthetic_data: bool = False,
    allow_synthetic_smoke: bool = False,
) -> Dict[str, Any]:
    """
    Secondary `G4` 2×2 `(X_ctrl, batch)` decomposition (`PREREG_I v5 §4.4` + `Addendum A G4 配对`):
      - Group 1 (`nadig_hepg2`, `nadig_jurkat`): positional `0..511` native `gem_group` alignment
        between `target_line` and its Nadig counterpart.
      - Group 2 & Group 3: fixed `(B_M = 0, B_S = 17)` across matched vs mismatched `X_ctrl`,
        plus spread across `{0, 1, 17, 35}`.
    """
    ckpt = checkpoint_for_line(target_line)
    strata = build_mismatched_strata_r16(target_line, ckpt)
    mis_all = list(strata["seen_by_checkpoint"]) + list(strata["unseen_by_checkpoint"])
    eval_perts = ref_bundle["eval_perts"]
    eval_genes = ref_bundle["eval_genes"]

    adapter = CleanRoomArcStateAdapter(
        model_key="STATE_SE",
        checkpoint_id=ckpt,
        models_root=models_root,
        status_dir=status_dir,
        device=device,
        neutral_batch_mode="index_zero",
        allow_synthetic_smoke=allow_synthetic_smoke,
        hgnc_bridge=hgnc,
    )
    adapter.load_weights()
    first_bundle = next(iter(adapter._bundles.values()))
    bmap = first_bundle.get("batch_onehot_idx_map", {})
    m_emb = se_embeddings[target_line][0]

    def _run_cell(ctrl_line: str, b_vec: Optional[Sequence[int]], b_scalar: int = 0) -> np.ndarray:
        c_emb = se_embeddings[ctrl_line][0]
        batch_indices = list(b_vec) if b_vec is not None else [int(b_scalar)] * int(c_emb.shape[0])
        res = adapter.predict_arm_deltas(
            target_line=target_line,
            control_source_line=ctrl_line,
            ctrl_matrix=None,
            source_genes=raw_pools[ctrl_line].gene_names,
            eval_perts=eval_perts,
            is_synthetic_data=is_synthetic_data,
            precomputed_ctrl_embed=c_emb,
            precomputed_sampling_audit=raw_pools[ctrl_line].sampling_audit,
            precomputed_gene_mapping=se_embeddings[ctrl_line][2],
            hgnc_bridge=hgnc,
            per_cell_batch_indices=batch_indices,
        )
        return project_model_deltas_to_eval_genes(res["delta_pred"], res["model_genes"], eval_genes, hgnc_bridge=hgnc)

    if target_line in ("nadig_hepg2", "nadig_jurkat"):
        counterpart = "nadig_jurkat" if target_line == "nadig_hepg2" else "nadig_hepg2"
        b_m = map_batch_labels_to_indices(raw_pools[target_line].batch_labels, bmap, fallback_idx=0).tolist()
        b_s_raw = map_batch_labels_to_indices(raw_pools[counterpart].batch_labels, bmap, fallback_idx=0).tolist()
        n_pos = m_emb.shape[0]
        if len(b_s_raw) < n_pos:
            raise ValueError(f"G4 Group 1 positional alignment length mismatch: {len(b_s_raw)} < {n_pos}")
        # Addendum A: Group 1 positional alignment 0..511
        b_s = b_s_raw[:n_pos]
        d_mm = _run_cell(target_line, b_m)
        d_ms = _run_cell(target_line, b_s)
        d_sm = _run_cell(counterpart, b_m)
        d_ss = _run_cell(counterpart, b_s)
        cells = {"MM": d_mm, "MS": d_ms, "SM": d_sm, "SS": d_ss}
        contrasts_in = {counterpart: cells}
        spread = None
    else:
        # Addendum A: Group 2 & Group 3 fixed B_M = 0, B_S = 17 + spread {0, 1, 17, 35}
        d_mm = _run_cell(target_line, None, b_scalar=0)
        d_ms = _run_cell(target_line, None, b_scalar=17)
        contrasts_in = {}
        for s in mis_all:
            contrasts_in[s] = {
                "MM": d_mm,
                "MS": d_ms,
                "SM": _run_cell(s, None, b_scalar=0),
                "SS": _run_cell(s, None, b_scalar=17),
            }
        spread = {}
        for bi in (0, 1, 17, 35):
            d_bi = _run_cell(target_line, None, b_scalar=bi)
            spread[str(bi)] = compute_orth_suite_r16(
                ref_bundle["delta_obs"], d_bi, ref_bundle["delta_transfer"], ref_bundle["target_gene_idx"],
                source_grand_mean=None, delta_half_a=ref_bundle["delta_half_a"], delta_half_b=ref_bundle["delta_half_b"],
                half_valid=ref_bundle["half_valid"], mode="all_nontarget", n_boot=n_boot, seed=seed,
            )

    def _decompose_pair(mm: np.ndarray, ms: np.ndarray, sm: np.ndarray, ss: np.ndarray) -> Dict[str, Any]:
        d_batch = compute_paired_contrast_r16(obs, mm, ms, tm, tg, None, ha, hb, hv, n_boot=n_boot, seed=seed)
        d_ctrl = compute_paired_contrast_r16(obs, mm, sm, tm, tg, None, ha, hb, hv, n_boot=n_boot, seed=seed)
        d_both = compute_paired_contrast_r16(obs, mm, ss, tm, tg, None, ha, hb, hv, n_boot=n_boot, seed=seed)
        b_est, c_est, both_est = float(d_batch["estimate"]), float(d_ctrl["estimate"]), float(d_both["estimate"])
        denom = abs(c_est) + abs(b_est)
        d_batch_dorth = compute_paired_contrast_r16(obs, mm, ms, tm, tg, gm, ha, hb, hv, n_boot=n_boot, seed=seed)
        d_ctrl_dorth = compute_paired_contrast_r16(obs, mm, sm, tm, tg, gm, ha, hb, hv, n_boot=n_boot, seed=seed)
        return {
            "delta_batch_orth": d_batch,
            "delta_ctrl_orth": d_ctrl,
            "delta_both_orth": d_both,
            "delta_interaction_orth": float(both_est - c_est - b_est),
            "batch_share_orth": float(abs(b_est) / denom) if denom > 1e-12 else float("nan"),
            "delta_batch_dorth": d_batch_dorth,
            "delta_ctrl_dorth": d_ctrl_dorth,
        }

    decomps: Dict[str, Any] = {}
    symmetry_decomps: Dict[str, Any] = {}
    obs = ref_bundle["delta_obs"]
    ha, hb, hv = ref_bundle["delta_half_a"], ref_bundle["delta_half_b"], ref_bundle["half_valid"]
    tm, gm, tg = ref_bundle["delta_transfer"], ref_bundle["source_grand_mean"], ref_bundle["target_gene_idx"]
    for s, mat in contrasts_in.items():
        decomps[s] = _decompose_pair(mat["MM"], mat["MS"], mat["SM"], mat["SS"])
        if target_line not in ("nadig_hepg2", "nadig_jurkat"):
            # Addendum A §6.2: mandatory (17, 0) symmetry check reusing exact same four cell matrices
            symmetry_decomps[s] = _decompose_pair(mat["MS"], mat["MM"], mat["SS"], mat["SM"])

    def _stratum_g4_summary(src_list: Sequence[str], src_map: Dict[str, Any]) -> Dict[str, Any]:
        vals = [src_map[s] for s in src_list if s in src_map]
        if not vals:
            return {"n_sources": 0}
        b_ests = [float(v["delta_batch_orth"]["estimate"]) for v in vals]
        c_ests = [float(v["delta_ctrl_orth"]["estimate"]) for v in vals]
        shares = [float(v["batch_share_orth"]) for v in vals if np.isfinite(v["batch_share_orth"])]
        return {
            "n_sources": len(vals),
            "median_delta_batch_orth": float(np.median(b_ests)),
            "median_delta_ctrl_orth": float(np.median(c_ests)),
            "median_batch_share_orth": float(np.median(shares)) if shares else float("nan"),
        }

    return {
        "target_line": target_line,
        "pairing": PREREG_I_V5_CONFIG["g4_pairing"],
        "decompositions": decomps,
        "symmetry_check_b17_vs_b0": symmetry_decomps if symmetry_decomps else None,
        "stratum_summaries": {
            "seen_by_checkpoint": _stratum_g4_summary(strata["seen_by_checkpoint"], decomps),
            "unseen_by_checkpoint": _stratum_g4_summary(strata["unseen_by_checkpoint"], decomps),
        },
        "batch_spread": spread,
    }


def run_real_r16(
    models: Sequence[str],
    target_lines: Sequence[str],
    status_dir: Path,
    output_dir: Path,
    models_root: Path = MODELS_ROOT_DEFAULT,
    raw_nadig_dir: Optional[Path] = None,
    raw_gse295214_dir: Optional[Path] = None,
    device: str = "auto",
    neutral_batch_mode: str = "index_zero",
    hvg_preproc_mode: str = "official_median_log1p",
    n_boot: int = 4000,
    seed: int = DEFAULT_R16_SEED,
    allow_npz_fallback_for_hvg: bool = False,
    include_resample_arm: bool = True,
    include_per_ckpt_group2: bool = True,
    include_hvg_perm100: bool = False,
    include_g4: bool = False,
) -> Dict[str, Any]:
    """
    Execute Gate 5 real-cell evaluation strictly behind `assert_freeze_gates_r16(status_dir)`.
    """
    t0 = time.time()
    gate_receipt = assert_freeze_gates_r16(status_dir)
    hgnc = get_hgnc_bridge()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    needed_ctrl_lines: set[str] = set(target_lines)
    for tl in target_lines:
        ckpt = checkpoint_for_line(tl)
        strata = build_mismatched_strata_r16(tl, ckpt)
        needed_ctrl_lines.update(strata["seen_by_checkpoint"])
        needed_ctrl_lines.update(strata["unseen_by_checkpoint"])
    if include_g4 and any(tl in ("nadig_hepg2", "nadig_jurkat") for tl in target_lines):
        needed_ctrl_lines.update(["nadig_hepg2", "nadig_jurkat"])

    raw_pools: Dict[str, Any] = {}
    resample_pools: Dict[str, Optional[Any]] = {}
    for cl in sorted(needed_ctrl_lines):
        want_res = bool(include_resample_arm and cl in target_lines)
        extracted = extract_raw_control_pool_r16(
            cl,
            custom_nadig_dir=raw_nadig_dir,
            custom_gse295214_dir=raw_gse295214_dir,
            n_ctrl=int(PREREG_I_V5_CONFIG["n_ctrl"]),
            cell_set_len=int(PREREG_I_V5_CONFIG["cell_set_len"]),
            seed=seed,
            allow_npz_fallback_for_hvg=allow_npz_fallback_for_hvg,
            return_resample=want_res,
            resample_seed=RESAMPLE_R16_SEED,
        )
        if want_res:
            pri_p, res_p = extracted  # type: ignore[misc]
            raw_pools[cl] = pri_p
            resample_pools[cl] = res_p
        else:
            raw_pools[cl] = extracted

    # Point 3: Export exact sampled cell barcodes for every control pool (+ non-overlapping resample)
    barcode_registry = {
        "seed_primary": int(seed),
        "seed_resample": int(RESAMPLE_R16_SEED),
        "primary_control_pools": {
            cl: {
                "n_sampled": int(len(rp.cell_barcodes)),
                "barcodes": rp.cell_barcodes.tolist(),
                "batches": rp.batch_labels.tolist(),
                "median_ntc_umi_sampled": float(rp.median_ntc_umi),
            }
            for cl, rp in sorted(raw_pools.items())
        },
        "resample_control_pools": {
            cl: (
                None
                if rp2 is None
                else {
                    "n_sampled": int(len(rp2.cell_barcodes)),
                    "barcodes": rp2.cell_barcodes.tolist(),
                    "batches": rp2.batch_labels.tolist(),
                    "median_ntc_umi_sampled": float(rp2.median_ntc_umi),
                    "n_overlap_with_primary": int(rp2.sampling_audit.get("n_overlap_with_primary", 0)),
                }
            )
            for cl, rp2 in sorted(resample_pools.items())
        },
    }
    bc_text = json.dumps(barcode_registry, indent=2, default=_json_default_r16)
    (output_dir / "sampled_control_barcodes_r16.json").write_text(bc_text, encoding="utf-8")
    (output_dir / "control_barcodes_r16.json").write_text(bc_text, encoding="utf-8")

    se_encoder: Optional[CleanRoomSE600MEncoder] = None
    cached_se_embeddings: Dict[str, Tuple[np.ndarray, Dict[str, Any], Dict[str, Any]]] = {}
    cached_se_resamp_embeddings: Dict[str, Tuple[np.ndarray, Dict[str, Any], Dict[str, Any]]] = {}
    se_emb_sha256_map: Dict[str, str] = {}
    se_emb_npz_path: Optional[Path] = None
    if "STATE_SE" in models:
        se_encoder = CleanRoomSE600MEncoder(
            se600m_dir=models_root / "SE-600M", device=device, batch_size=16, hgnc_bridge=hgnc
        )
        se_payload: Dict[str, Any] = {}
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
                use_hgnc_rescue=False,
                dataset_name=str(cl),
            )
            cached_se_embeddings[cl] = (emb, rp.sampling_audit, map_stats)
            sha_pri = str(map_stats["emb_sha256"])
            se_emb_sha256_map[cl] = sha_pri
            barcode_registry["primary_control_pools"][cl]["se600m_x_state_sha256"] = sha_pri
            se_payload[f"{cl}__X_state"] = np.asarray(emb, dtype=np.float32)
            se_payload[f"{cl}__barcodes"] = np.asarray(rp.cell_barcodes, dtype=str)
            se_payload[f"{cl}__sha256"] = np.asarray(sha_pri)
            rp2 = resample_pools.get(cl)
            if rp2 is not None and not rp2.is_log1p:
                emb2, ms2 = se_encoder.encode_cells(
                    ctrl_counts=rp2.raw_counts,
                    source_genes=rp2.gene_names,
                    source_ensembl_ids=rp2.gene_ids,
                    seed=RESAMPLE_R16_SEED,
                    hgnc_bridge=hgnc,
                    use_hgnc_rescue=False,
                    dataset_name=f"{cl}__resample",
                )
                cached_se_resamp_embeddings[cl] = (emb2, rp2.sampling_audit, ms2)
                sha_res = str(ms2["emb_sha256"])
                res_key = f"{cl}__resample"
                se_emb_sha256_map[res_key] = sha_res
                if barcode_registry["resample_control_pools"].get(cl) is not None:
                    barcode_registry["resample_control_pools"][cl]["se600m_x_state_sha256"] = sha_res
                se_payload[f"{res_key}__X_state"] = np.asarray(emb2, dtype=np.float32)
                se_payload[f"{res_key}__barcodes"] = np.asarray(rp2.cell_barcodes, dtype=str)
                se_payload[f"{res_key}__sha256"] = np.asarray(sha_res)

        se_emb_npz_path = output_dir / "se600m_embeddings_r16.npz"
        np.savez(se_emb_npz_path, **se_payload)
        bc_text = json.dumps(barcode_registry, indent=2, default=_json_default_r16)
        (output_dir / "sampled_control_barcodes_r16.json").write_text(bc_text, encoding="utf-8")
        (output_dir / "control_barcodes_r16.json").write_text(bc_text, encoding="utf-8")

    summary_manifest: Dict[str, Any] = {
        "prereg_version": PREREG_I_V5_CONFIG["prereg_version"],
        "gate_receipt": gate_receipt,
        "device": device,
        "neutral_batch_mode": neutral_batch_mode,
        "hvg_preproc_mode": hvg_preproc_mode,
        "seed": int(seed),
        "resample_seed": int(RESAMPLE_R16_SEED),
        "n_boot": int(n_boot),
        "barcodes_file": str(output_dir / "sampled_control_barcodes_r16.json"),
        "se600m_embeddings_file": str(se_emb_npz_path) if se_emb_npz_path is not None else None,
        "se600m_embeddings_sha256": se_emb_sha256_map,
        "sampled_control_barcodes": barcode_registry,
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
                se600m_use_hgnc_rescue=False,
            )
            adapter.load_weights()
            ref_bundle_mgenes = next(iter(adapter._bundles.values()))["gene_names"]
            ref_bundle_perts = sorted(
                k for k in next(iter(adapter._bundles.values()))["pert_lookup"].keys() if k != "non-targeting"
            )

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
                s: (
                    None if model_key == "STATE_SE" else raw_pools[s].raw_counts,
                    raw_pools[s].gene_names,
                    raw_pools[s].batch_labels,
                )
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
                matched_barcodes=matched_rp.cell_barcodes,
                mismatched_barcodes={s: raw_pools[s].cell_barcodes for s in mis_sources},
            )

            # Optional same-line non-overlapping resample arm (Addendum A Item 8)
            rp_res = resample_pools.get(tl)
            res_out = None
            delta_res_eval = None
            if rp_res is not None:
                res_pre = cached_se_resamp_embeddings.get(tl) if model_key == "STATE_SE" else None
                res_out = adapter.predict_arm_deltas(
                    target_line=tl,
                    control_source_line=f"{tl}__resample",
                    ctrl_matrix=None if model_key == "STATE_SE" else rp_res.raw_counts,
                    source_genes=rp_res.gene_names,
                    eval_perts=ref_data["eval_perts"],
                    batch_labels=rp_res.batch_labels,
                    input_is_log1p=bool(rp_res.is_log1p),
                    seed=RESAMPLE_R16_SEED,
                    n_ctrl=int(len(rp_res.cell_barcodes)),
                    cell_set_len=int(PREREG_I_V5_CONFIG["cell_set_len"]),
                    is_synthetic_data=False,
                    cell_total_umis=rp_res.cell_total_umis,
                    target_median_umi=rp_res.median_ntc_umi,
                    precomputed_ctrl_embed=res_pre[0] if res_pre else None,
                    precomputed_sampling_audit=res_pre[1] if res_pre else rp_res.sampling_audit,
                    precomputed_gene_mapping=res_pre[2] if res_pre else None,
                    hgnc_bridge=hgnc,
                    pre_sampled=True,
                    cell_barcodes=rp_res.cell_barcodes,
                )
                delta_res_eval = project_model_deltas_to_eval_genes(
                    res_out["delta_pred"], res_out["model_genes"], ref_data["eval_genes"], hgnc_bridge=hgnc
                )

            delta_m_eval = project_model_deltas_to_eval_genes(
                arm_out["matched"]["delta_pred"], arm_out["matched"]["model_genes"], ref_data["eval_genes"], hgnc_bridge=hgnc
            )
            delta_mis_eval = {
                s: project_model_deltas_to_eval_genes(
                    arm_out["mismatched"][s]["delta_pred"], arm_out["mismatched"][s]["model_genes"], ref_data["eval_genes"], hgnc_bridge=hgnc
                )
                for s in mis_sources
            }

            arm_meta = {
                "matched_sampling": raw_pools[tl].sampling_audit,
                "resample_sampling": rp_res.sampling_audit if rp_res is not None else None,
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
                matched_barcodes=matched_rp.cell_barcodes,
                mismatched_barcodes={s: raw_pools[s].cell_barcodes for s in mis_sources},
                delta_pred_resample_eval=delta_res_eval,
                resample_barcodes=rp_res.cell_barcodes if rp_res is not None else None,
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
                delta_pred_resample_eval=delta_res_eval,
            )

            # Group 2/3 per-checkpoint breakdown (§4.3.2)
            if include_per_ckpt_group2 and len(arm_out["matched"]["per_ckpt_deltas"]) > 1:
                per_ckpt_eval: Dict[str, Any] = {}
                for ck_name, dm_ck in sorted(arm_out["matched"]["per_ckpt_deltas"].items()):
                    dm_ck_e = project_model_deltas_to_eval_genes(dm_ck, ref_bundle_mgenes, ref_data["eval_genes"], hgnc_bridge=hgnc)
                    dmis_ck_e = {
                        s: project_model_deltas_to_eval_genes(
                            arm_out["mismatched"][s]["per_ckpt_deltas"][ck_name],
                            ref_bundle_mgenes,
                            ref_data["eval_genes"],
                            hgnc_bridge=hgnc,
                        )
                        for s in mis_sources
                    }
                    per_ckpt_eval[ck_name] = evaluate_and_adjudicate_line(
                        ref_bundle=ref_data,
                        delta_pred_matched_eval=dm_ck_e,
                        delta_pred_mis_eval=dmis_ck_e,
                        seen_sources=strata["seen_by_checkpoint"],
                        unseen_sources=strata["unseen_by_checkpoint"],
                        n_boot=n_boot,
                        seed=seed,
                        delta=float(PREREG_I_V5_CONFIG["tost_delta"]),
                    )
                eval_res["per_checkpoint_breakdown"] = per_ckpt_eval

            if include_hvg_perm100 and model_key == "STATE_HVG":
                eval_res["hvg_perm100"] = compute_gene_label_perm_null_r16(
                    delta_obs=ref_data["delta_obs"],
                    pred_matched=delta_m_eval,
                    delta_transfer=ref_data["delta_transfer"],
                    target_gene_idx=ref_data["target_gene_idx"],
                    source_grand_mean=ref_data["source_grand_mean"],
                    delta_half_a=ref_data["delta_half_a"],
                    delta_half_b=ref_data["delta_half_b"],
                    half_valid=ref_data["half_valid"],
                    n_perm=100,
                    matched_orth_suite=eval_res["orth_matched"],
                    seed=seed,
                )

            if include_g4 and model_key == "STATE_SE":
                eval_res["g4_factorial"] = run_g4_factorial_for_line(
                    target_line=tl,
                    ref_bundle=ref_data,
                    raw_pools=raw_pools,
                    se_embeddings=cached_se_embeddings,
                    models_root=models_root,
                    status_dir=status_dir,
                    device=device,
                    n_boot=n_boot,
                    seed=seed,
                    hgnc=hgnc,
                    is_synthetic_data=False,
                    allow_synthetic_smoke=False,
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
                    default=_json_default_r16,
                ),
                encoding="utf-8",
            )
            model_entries[tl] = {
                "npz_path": str(npz_path),
                "metrics_json": str(metrics_json),
                "h_i1_summary": eval_res["h_i1"]["line_summary"],
                "h_i2_class": eval_res["h_i2"]["classification"],
            }

        grp = adjudicate_group_r16(per_line_hi1, per_line_hi2)
        summary_manifest["models"][model_key] = {"lines": model_entries, "group_rollup": grp}

    summary_manifest["elapsed_sec"] = round(float(time.time() - t0), 2)
    (output_dir / "manifest_r16_audit.json").write_text(json.dumps(summary_manifest, indent=2, default=_json_default_r16), encoding="utf-8")
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
      - Verifies Gate-2/v5 `ST-HVG` (direct `preds`, `16 × 32` `padded=True`, `official_median_log1p`
        default, fixed `batch = 0`, `Perm100` null), `EXACT_INVARIANCE`, `resample` arm, barcode export,
        `.npz` schema, cross-party tolerance verifier, and adjudicator outputs.
      - When `include_real_se600m=True`, also runs real `SE-600M` (`MPS float32`, exact symbol lookup)
        + `ST-SE` strictly on synthetic Poisson count matrices.
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
    remap_pairs = {"AARS": "AARS1", "DARS": "DARS1", "H2AFZ": "H2AZ1", "TAZ": "TAFAZZIN"}
    src_genes = [remap_pairs.get(g, g) for g in g_model]
    canonical_targets = [
        hgnc.canonicalize(t) if hgnc else t
        for t in ("AARS", "DARS", "H2AFZ", "TAZ", "AAAS", "AAMP", "AAR2", "AARS2", "AATF", "ABCB7")
    ][:n_perts]

    rng = np.random.default_rng(int(seed))
    g_len = len(src_genes)
    matched_raw = rng.poisson(2.0, size=(n_ctrl + 64, g_len)).astype(np.float32)
    matched_b = ["B1"] * ((n_ctrl + 64) // 2) + ["B2"] * (n_ctrl + 64 - (n_ctrl + 64) // 2)
    matched_bcs = [f"SYN_M_{i:04d}" for i in range(n_ctrl + 64)]
    mis_a_raw = rng.poisson(7.0, size=(n_ctrl + 32, g_len)).astype(np.float32)
    mis_b_raw = rng.poisson(1.0, size=(n_ctrl + 32, g_len)).astype(np.float32)
    mis_a_bcs = [f"SYN_SA_{i:04d}" for i in range(n_ctrl + 32)]
    mis_b_bcs = [f"SYN_SB_{i:04d}" for i in range(n_ctrl + 32)]

    arm_res = adapter_hvg.run_matched_and_mismatched_arms(
        target_line="a549_ko",
        matched_ctrl_matrix=matched_raw,
        matched_genes=src_genes,
        mismatched_ctrls={
            "nadig_jurkat": (mis_a_raw, src_genes, ["L1"] * (n_ctrl + 32)),
            "dld1": (mis_b_raw, src_genes, ["U1"] * (n_ctrl + 32)),
        },
        eval_perts=canonical_targets,
        matched_batch_labels=matched_b,
        input_is_log1p=False,
        seed=seed,
        n_ctrl=n_ctrl,
        cell_set_len=32,
        is_synthetic_data=True,
        hgnc_bridge=hgnc,
        matched_barcodes=matched_bcs,
        mismatched_barcodes={"nadig_jurkat": mis_a_bcs, "dld1": mis_b_bcs},
    )

    perm5 = adapter_hvg.predict_hvg_perm100_deltas(
        target_line="a549_ko",
        ctrl_sets_embed=arm_res["matched"]["ctrl_sets_embed"],
        eval_perts=canonical_targets,
        n_perm=5,
        is_synthetic_data=True,
        hgnc_bridge=hgnc,
    )

    d_m_raw = arm_res["matched"]["delta_pred"]
    eval_g_canon: List[str] = []
    seen_eg: set[str] = set()
    for g in g_model:
        cg = hgnc.canonicalize(g) if hgnc else g
        if cg and cg not in seen_eg:
            seen_eg.add(cg)
            eval_g_canon.append(cg)
        if len(eval_g_canon) >= 240:
            break

    d_m = project_model_deltas_to_eval_genes(d_m_raw, g_model, eval_g_canon, hgnc_bridge=hgnc)
    d_mis = {
        k: project_model_deltas_to_eval_genes(arm_res["mismatched"][k]["delta_pred"], g_model, eval_g_canon, hgnc_bridge=hgnc)
        for k in ("nadig_jurkat", "dld1")
    }
    p_len = len(canonical_targets)
    fake_ref = {
        "eval_genes": eval_g_canon,
        "eval_perts": canonical_targets,
        "target_gene_idx": np.full(p_len, -1, dtype=np.int64),
        "delta_obs": d_m + rng.normal(0.0, 0.02, size=d_m.shape),
        "delta_half_a": d_m + rng.normal(0.0, 0.03, size=d_m.shape),
        "delta_half_b": d_m + rng.normal(0.0, 0.03, size=d_m.shape),
        "half_valid": np.ones(p_len, dtype=bool),
        "delta_transfer": rng.normal(0.0, 0.01, size=d_m.shape),
        "source_grand_mean": np.zeros(len(eval_g_canon), dtype=np.float64),
        "delta_fct": d_m * 0.5,
        "ctrl_mean_libnorm": rng.uniform(0.1, 3.0, size=len(eval_g_canon)),
    }
    eval_summary = evaluate_and_adjudicate_line(
        fake_ref, d_m, d_mis, seen_sources=["nadig_jurkat"], unseen_sources=["dld1"],
        n_boot=200, seed=seed, delta_pred_resample_eval=d_m.copy(),
    )
    out: Dict[str, Any] = {
        "device": adapter_hvg.device,
        "delta_shape": d_m.shape,
        "perm5_shape": perm5.shape,
        "l2_mis_jurkat": float(np.linalg.norm(d_m - d_mis["nadig_jurkat"])),
        "l2_mis_dld1": float(np.linalg.norm(d_m - d_mis["dld1"])),
        "matched_barcodes_len": len(arm_res["matched"]["selected_barcodes"] or []),
        "hgnc_rescued_hvg": arm_res["matched"]["gene_mapping"]["n_hgnc_canon_rescued"],
        "h_i1_line_summary": eval_summary["h_i1"]["line_summary"],
    }

    if output_dir is not None:
        npz_out = export_line_arm_npz(
            Path(output_dir) / "synthetic_smoke_hvg__a549_ko.npz",
            model_key="STATE_HVG",
            checkpoint_id="zeroshot/hepg2",
            target_line="a549_ko",
            ref_bundle=fake_ref,
            delta_pred_matched_eval=d_m,
            delta_pred_mis_eval=d_mis,
            seen_sources=["nadig_jurkat"],
            unseen_sources=["dld1"],
            neutral_batch_mode="index_zero",
            delta_mode="minus_pred_ntc",
            hvg_preproc_mode="official_median_log1p",
            arm_metadata={"synthetic": True},
            matched_barcodes=arm_res["matched"]["selected_barcodes"],
            mismatched_barcodes={
                "nadig_jurkat": arm_res["mismatched"]["nadig_jurkat"]["selected_barcodes"],
                "dld1": arm_res["mismatched"]["dld1"]["selected_barcodes"],
            },
            delta_pred_resample_eval=d_m,
            resample_barcodes=arm_res["matched"]["selected_barcodes"],
        )
        out["npz_path"] = str(npz_out)

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
            se600m_use_hgnc_rescue=False,
        )
        arm_se = adapter_se.run_matched_and_mismatched_arms(
            target_line="a549_ko",
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
            matched_barcodes=[f"SE_M_{i}" for i in range(n_ctrl)],
            mismatched_barcodes={"cd4t_rest": [f"SE_S_{i}" for i in range(n_ctrl)]},
        )
        d_se_m = arm_se["matched"]["delta_pred"]
        d_se_s = arm_se["mismatched"]["cd4t_rest"]["delta_pred"]
        out["se600m_st_se_smoke"] = {
            "device": se_enc.device_str,
            "delta_matched_shape": d_se_m.shape,
            "l2_ctrl_swap": float(np.linalg.norm(d_se_m - d_se_s)),
            "mapping": arm_se["matched"]["gene_mapping"],
            "centroid_diag": arm_se["control_centroid_diagnostics"]["cd4t_rest"],
        }

    return out


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Clean-Room Round 16 Real-Data & Synthetic Runner (A-T3, PREREG_I v5 + Addendum A).")
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
    p.add_argument("--raw-gse295214-dir", type=Path, default=None)
    p.add_argument("--device", choices=["auto", "mps", "cpu"], default="auto")
    p.add_argument(
        "--neutral-batch-mode",
        default="index_zero",
        choices=["index_zero", "bypass_zero", "mean_batch_emb", "multi_batch_ens8", "per_cell_native"],
    )
    p.add_argument("--hvg-preproc", default="official_median_log1p", choices=["official_median_log1p", "libnorm1e4_log1p"])
    p.add_argument("--n-boot", type=int, default=4000)
    p.add_argument("--seed", type=int, default=DEFAULT_R16_SEED)
    p.add_argument("--include-real-se600m-synthetic", action="store_true")
    p.add_argument("--allow-npz-fallback-for-hvg", action="store_true")
    p.add_argument("--no-resample-arm", action="store_true")
    p.add_argument("--include-hvg-perm100", action="store_true")
    p.add_argument("--include-g4", action="store_true")
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
            raw_gse295214_dir=args.raw_gse295214_dir,
            device=args.device,
            neutral_batch_mode=args.neutral_batch_mode,
            hvg_preproc_mode=args.hvg_preproc,
            n_boot=args.n_boot,
            seed=args.seed,
            allow_npz_fallback_for_hvg=args.allow_npz_fallback_for_hvg,
            include_resample_arm=not args.no_resample_arm,
            include_hvg_perm100=args.include_hvg_perm100,
            include_g4=args.include_g4,
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
