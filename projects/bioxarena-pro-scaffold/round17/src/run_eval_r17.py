#!/usr/bin/env python3
"""
Round 17 Pre-Registered Confirmatory Evaluation Runner (`PREREG_J` v2).

Implements the exact `PREREG_J_DRAFT_v2.md` protocol (hash-locked at commit `9789b0b`):
  1. Hard Execution Gate (`assert_all_execution_gates_r17`): Refuses ANY real-cell data
     access or model inference unless BOTH `round17/status/R17_CODE_FREEZE.txt` AND
     `round17/status/R17_GO.txt` (plus `round17/prereg/PREREG_J_LOCK.md`) exist and are
     non-empty.
  2. Lines (`PREREG_J §2.1, §2.2`, Controller Dispatch Item 2):
       - Confirmatory (`K_conf = 3`): `cd4t_stim48`, `gse295214_h1`, `pfizer_telohaec_basal`.
       - Within-family replication: `pfizer_telohaec_cytokine` (never counted in `K_conf` or `k`).
       - Secondary fresh descriptive: `a549_myllia`, `hek293_xucao`.
       - Descriptive micro-N: `telohaec`, `calu3`, `mac_lps`, `mac_ctrl`.
       - R16 subset sensitivity: `a549_ko_cas9pos`.
       - Explicitly omitted: `wtc11_i3neuron` (`N = 1`, Controller Ruling 2).
  3. Auditable `.npz` Array & Embedding Exports (`round17/src/R17_ARRAY_SCHEMA.md`):
     Exports per-line `observed.npz`, per-model `preds_{STATE_SE,STATE_HVG}_sweep56.npz` &
     `preds_{STATE_SE,STATE_HVG}_g4_mismatched.npz`, `array_manifest_sha256.json`,
     `embeddings/x_state_<pool>.npz`, `control_barcodes_r17.json`, and `summary_r17.json`
     so Master Controller and `A-T3` can recompute EVERY number from arrays alone
     (`recompute_cohort_from_exported_arrays_r17`).
  4. `--synthetic-smoke`: Runs the complete end-to-end evaluation, adjudication, Holm-Bonferroni,
     pooled-Holm, `J5a` cross-fitted `G4`, `J5b` Spearman + permutation null, `.npz`/`.json`
     export, and Stage 4 from-disk recomputation parity check strictly on CPU with 100%
     synthetic tensors (ZERO real cells touched).
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence, Set, Tuple

import numpy as np

ROUND17_SRC_DIR = Path(__file__).resolve().parent
ROUND17_ROOT = ROUND17_SRC_DIR.parent
REPO_ROOT = ROUND17_ROOT.parent
R16_SRC_DIR = REPO_ROOT / "round16" / "src"
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
for _p in (str(ROUND17_SRC_DIR), str(R16_SRC_DIR), str(R10_SRC_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from data_loader_r10 import (  # noqa: E402
    DEFAULT_HGNC_FILE,
    HGNCBridge,
    LineDatasetR10,
    R9_COMMON_GENES_FILE,
    R9_LINES_DIR,
    load_r9_common_genes,
    load_standardized_line_r10,
)
from fct_pipeline import (  # noqa: E402
    FctSourceState,
    build_source_input_r10,
    compute_source_transfer_and_fc_bar,
    fct_from_fc_bar,
)
from inference_r17 import (  # noqa: E402
    run_state_g4_mismatched_r17,
    run_state_sweep56_r17,
)
from model_adapters_r16 import (  # noqa: E402
    DEFAULT_MODELS_DIR,
    MULTI_BATCH_ENS8_INDICES,
    STATE_CHECKPOINT_FOLDS,
    compute_embedding_separation_diagnostics,
    encode_controls_se600m_r16,
    preprocess_hvg_counts_r16,
)
from ntc_extractors_r17 import (  # noqa: E402
    ALL_EVALUATED_LINES_R17,
    CONFIRMATORY_LINES_R17,
    DEFAULT_GSE295214_RAW_DIR,
    DEFAULT_PFIZER_RAW_PATH,
    DEFAULT_R17_PREREG_LOCK_PATH,
    DEFAULT_R17_RESULTS_DIR,
    DEFAULT_R17_STATUS_DIR,
    LINE_SPECS_R17,
    MICRO_N_LINES_R17,
    SECONDARY_FRESH_LINES_R17,
    SUBSET_SENSITIVITY_LINES_R17,
    WITHIN_FAMILY_LINES_R17,
    assert_all_execution_gates_r17,
    build_mismatched_strata_r17,
    canonicalize_gene_r17,
    extract_control_ntc_bundle_r17,
    validate_target_lines_r17,
)
from run_eval_r16 import (  # noqa: E402
    load_state_vocab_r16,
    save_control_embeddings_r16,
)
from stats_r17 import (  # noqa: E402
    DEFAULT_PRACTICAL_DELTA,
    DEFAULT_R17_N_BOOT,
    DEFAULT_R17_N_BOOT_BONF56,
    DEFAULT_R17_N_PERM,
    DEFAULT_R17_N_SPLITS_DIAG,
    DEFAULT_R17_RESAMPLE_SEED,
    DEFAULT_R17_SEED,
    N_BATCH_INDICES_STATE,
    adjudicate_cohort_r17,
    compute_cross_fitted_spread_r17,
    compute_sweep56_moments_r17,
    evaluate_model_on_line_r17,
    make_stratified_half_split_r17,
)

MODELS_R17: Tuple[str, ...] = ("STATE_SE", "STATE_HVG")


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ==============================================================================
# Line Alignment (`PREREG_I` / `PREREG_J` Parity, Allowing Micro-N Gene Counts)
# ==============================================================================


def align_line_and_sources_r17(
    line_name: str,
    fct_source_names: Sequence[str],
    line_dir: Path,
    g_model_2000: Sequence[str],
    p_model_canon_to_key: Mapping[str, str],
    hgnc_bridge: HGNCBridge,
    common_genes_r9: Sequence[str],
    fixed_target_names: Optional[Sequence[str]] = None,
    fixed_gene_names: Optional[Sequence[str]] = None,
) -> Dict[str, Any]:
    """
    Align target line and locked `FCT` sources onto `(T_eval^R17, G_eval^STATE)`, canonically
    sorted by target symbol. Uses `min_gene_overlap=1` so micro-N `telohaec` (`|G_eval^R9|=4869`)
    loads cleanly while preserving bit-exact `FCT`, `Transfer_Mean`, and `Source_Grand_Mean`
    parity with frozen `round16/src/run_eval_r16.py::_align_line_and_sources_r16`.
    """
    ds_t = load_standardized_line_r10(
        line_dir,
        representation="log1p",
        common_genes=common_genes_r9,
        hgnc_bridge=hgnc_bridge,
        min_gene_overlap=1,
    )
    src_ds_map: Dict[str, LineDatasetR10] = {
        sn: load_standardized_line_r10(
            R9_LINES_DIR / sn,
            representation="log1p",
            common_genes=common_genes_r9,
            hgnc_bridge=hgnc_bridge,
            min_gene_overlap=1,
        )
        for sn in fct_source_names
    }

    shared_r9_genes_set = set(ds_t.genes).intersection(*(set(src_ds_map[sn].genes) for sn in fct_source_names))
    if fixed_gene_names is not None:
        missing_g = [g for g in fixed_gene_names if g not in shared_r9_genes_set]
        if missing_g:
            raise ValueError(
                f"[{line_name}] fixed_gene_names has {len(missing_g)} genes absent from source pool "
                f"{list(fct_source_names)} (e.g. {missing_g[:5]})."
            )
        g_eval = list(fixed_gene_names)
        g_r9_set = set(shared_r9_genes_set)
    else:
        g_r9_aligned = [g for g in ds_t.genes if g in shared_r9_genes_set]
        g_r9_set = set(g_r9_aligned)
        g_model_canon_set = {canonicalize_gene_r17(hgnc_bridge, g) for g in g_model_2000}
        g_eval = [g for g in g_r9_aligned if g in g_model_canon_set]

    g_eval_map = {g: i for i, g in enumerate(g_eval)}
    t_pos = {g: i for i, g in enumerate(ds_t.genes)}
    t_cols = np.asarray([t_pos[g] for g in g_eval], dtype=np.int64)

    if fixed_target_names is not None:
        t_row_map = {t: i for i, t in enumerate(ds_t.target_names)}
        missing_t = [t for t in fixed_target_names if t not in t_row_map]
        if missing_t:
            raise ValueError(f"[{line_name}] fixed_target_names missing {missing_t[:5]}.")
        row_arr = np.asarray([t_row_map[t] for t in fixed_target_names], dtype=np.int64)
        q_names = list(fixed_target_names)
    else:
        src_repr_targets: Set[str] = set()
        for sn in fct_source_names:
            s_ds = src_ds_map[sn]
            if s_ds.repr_mask is not None:
                src_repr_targets.update(t for t, ok in zip(s_ds.target_names, s_ds.repr_mask) if ok)
            else:
                src_repr_targets.update(s_ds.target_names)
        keep_pairs = [
            (ds_t.target_names[i], i)
            for i, t in enumerate(ds_t.target_names)
            if t in g_r9_set
            and (ds_t.repr_mask is None or bool(ds_t.repr_mask[i]))
            and t in src_repr_targets
            and t in p_model_canon_to_key
        ]
        # Canonical ascending target-symbol sort (`PREREG_J §3.2`)
        keep_pairs.sort(key=lambda x: (str(x[0]), int(x[1])))
        q_names = [p[0] for p in keep_pairs]
        row_arr = np.asarray([p[1] for p in keep_pairs], dtype=np.int64)

    q_gidx = np.asarray([g_eval_map.get(t, -1) for t in q_names], dtype=np.int64)

    src_inputs = []
    for sn in fct_source_names:
        s_ds = src_ds_map[sn]
        s_pos = {g: i for i, g in enumerate(s_ds.genes)}
        s_cols = np.asarray([s_pos[g] for g in g_eval], dtype=np.int64)
        s_tidx = np.asarray([g_eval_map.get(t, -1) for t in s_ds.target_names], dtype=np.int64)
        s_sub = LineDatasetR10(
            name=s_ds.name,
            genes=list(g_eval),
            target_names=list(s_ds.target_names),
            target_indices=s_tidx,
            deltas=s_ds.deltas[:, s_cols],
            ctrl_cells=s_ds.ctrl_cells[:, s_cols],
            ctrl_mean_log1p=s_ds.ctrl_mean_log1p[s_cols] if s_ds.ctrl_mean_log1p is not None else None,
            ctrl_mean_libnorm=s_ds.ctrl_mean_libnorm[s_cols] if s_ds.ctrl_mean_libnorm is not None else None,
            ctrl_mean_libnorm_full=(
                s_ds.ctrl_mean_libnorm_full[s_cols] if s_ds.ctrl_mean_libnorm_full is not None else None
            ),
            repr_mask=s_ds.repr_mask,
            half_a_deltas=s_ds.half_a_deltas[:, s_cols] if s_ds.half_a_deltas is not None else None,
            half_b_deltas=s_ds.half_b_deltas[:, s_cols] if s_ds.half_b_deltas is not None else None,
            half_valid=s_ds.half_valid,
            n_cells=s_ds.n_cells,
            knockdown_depth_log1p=s_ds.knockdown_depth_log1p,
            knockdown_fraction=s_ds.knockdown_fraction,
            target_ensembl_ids=s_ds.target_ensembl_ids,
            representation=s_ds.representation,
            metadata=s_ds.metadata,
        )
        src_inputs.append(build_source_input_r10(s_sub, only_reproducible=True, compute_cov=False))

    fct_state = FctSourceState.fit(src_inputs, mask_self=True)
    delta_transfer, fc_bar, has_hits = compute_source_transfer_and_fc_bar(fct_state, q_names)
    mu_target = np.maximum(np.asarray(ds_t.ctrl_mean_libnorm[t_cols], dtype=np.float64), 0.0)
    mu_sources = np.vstack([np.maximum(s.ctrl_mean_libnorm, 0.0) for s in src_inputs])
    fct_pred = fct_from_fc_bar(mu_target, fc_bar, delta_transfer, has_hits)
    c_med = float(np.median(mu_target))
    e1_pred = fct_from_fc_bar(np.full_like(mu_target, c_med), fc_bar, delta_transfer, has_hits)

    return {
        "gene_names": g_eval,
        "target_names": q_names,
        "target_indices": q_gidx,
        "delta_obs": ds_t.deltas[row_arr][:, t_cols],
        "delta_half_a": ds_t.half_a_deltas[row_arr][:, t_cols],
        "delta_half_b": ds_t.half_b_deltas[row_arr][:, t_cols],
        "half_valid": ds_t.half_valid[row_arr],
        "delta_transfer": delta_transfer,
        "source_grand_mean": fct_state.grand_mean,
        "mu_target": mu_target,
        "mu_sources": mu_sources,
        "fc_bar": fc_bar,
        "has_hits": has_hits,
        "fct_delta_pred": fct_pred,
        "e1_constmed_delta_pred": e1_pred,
    }


# ==============================================================================
# Auditable `.npz` Export & Independent From-Disk Stage 4 Recompute
# ==============================================================================


def save_line_arrays_r17(
    arrays_dir: Path,
    line_name: str,
    observed_payload: Mapping[str, Any],
    model_exports: Mapping[str, Mapping[str, Any]],
    extra_fct_sensitivity: Optional[Mapping[str, Mapping[str, np.ndarray]]] = None,
) -> Dict[str, str]:
    """
    Save auditable `.npz` files under `arrays_dir / line_name` per `round17/src/R17_ARRAY_SCHEMA.md`:
      - `observed.npz`
      - `preds_<MODEL>_sweep56.npz`
      - `preds_<MODEL>_g4_mismatched.npz` (when present)
      - `SENS_<key>.npz` (optional, e.g., `gse295214_h1` full `r9_6` sensitivity)
      - `array_manifest_sha256.json`
    """
    ldir = Path(arrays_dir) / line_name
    ldir.mkdir(parents=True, exist_ok=True)
    manifest: Dict[str, str] = {}

    obs_path = ldir / "observed.npz"
    np.savez_compressed(obs_path, **dict(observed_payload))
    manifest[obs_path.name] = _sha256_file(obs_path)

    for mkey, payload in sorted(model_exports.items()):
        sw_path = ldir / f"preds_{mkey}_sweep56.npz"
        np.savez_compressed(sw_path, **dict(payload["sweep56"]))
        manifest[sw_path.name] = _sha256_file(sw_path)

        if payload.get("g4_mismatched") is not None:
            g4_path = ldir / f"preds_{mkey}_g4_mismatched.npz"
            np.savez_compressed(g4_path, **dict(payload["g4_mismatched"]))
            manifest[g4_path.name] = _sha256_file(g4_path)

    if extra_fct_sensitivity:
        for skey, spay in sorted(extra_fct_sensitivity.items()):
            spath = ldir / f"SENS_{skey}.npz"
            np.savez_compressed(
                spath,
                fct_delta_pred=np.asarray(spay["fct_delta_pred"], dtype=np.float64),
                delta_transfer=np.asarray(spay["delta_transfer"], dtype=np.float64),
                source_grand_mean=np.asarray(spay["source_grand_mean"], dtype=np.float64),
                target_names=np.asarray(observed_payload["target_names"], dtype=str),
                gene_names=np.asarray(observed_payload["gene_names"], dtype=str),
            )
            manifest[spath.name] = _sha256_file(spath)

    (ldir / "array_manifest_sha256.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def evaluate_cohort_from_arrays_r17(
    line_bundles: Mapping[str, Dict[str, Any]],
    arrays_dir: Optional[Path] = None,
    n_boot: int = DEFAULT_R17_N_BOOT,
    n_boot_bonf56: int = DEFAULT_R17_N_BOOT_BONF56,
    n_perm: int = DEFAULT_R17_N_PERM,
    seed: int = DEFAULT_R17_SEED,
    n_splits_diag: int = DEFAULT_R17_N_SPLITS_DIAG,
) -> Dict[str, Any]:
    """
    Pure array-only Round 17 cohort evaluator & auditable `.npz` exporter.
    Used identically by real-data execution, `--synthetic-smoke`, and Stage 4 from-disk recomputation.
    """
    per_line_results: Dict[str, Dict[str, Any]] = {}
    line_roles: Dict[str, str] = {}
    line_n_eval: Dict[str, int] = {}
    line_manifests: Dict[str, Dict[str, str]] = {}
    line_meta_out: Dict[str, Dict[str, Any]] = {}

    for ln in sorted(line_bundles.keys()):
        b = line_bundles[ln]
        role = str(b["role"])
        line_roles[ln] = role
        target_names = [str(t) for t in b["target_names"]]
        gene_names = [str(g) for g in b["gene_names"]]
        line_n_eval[ln] = len(target_names)

        delta_obs = np.asarray(b["delta_obs"], dtype=np.float64)
        delta_half_a = np.asarray(b["delta_half_a"], dtype=np.float64)
        delta_half_b = np.asarray(b["delta_half_b"], dtype=np.float64)
        half_valid = np.asarray(b["half_valid"], dtype=bool)
        delta_transfer = np.asarray(b["delta_transfer"], dtype=np.float64)
        source_grand_mean = np.asarray(b["source_grand_mean"], dtype=np.float64)
        target_indices = np.asarray(b["target_indices"], dtype=np.int64)
        fct_delta_pred = np.asarray(b["fct_delta_pred"], dtype=np.float64)
        e1_pred = np.asarray(b.get("e1_constmed_delta_pred", np.zeros_like(fct_delta_pred)), dtype=np.float64)
        mu_target = np.asarray(b.get("mu_target", np.zeros(len(gene_names))), dtype=np.float64)
        mu_sources = np.asarray(b.get("mu_sources", np.zeros((1, len(gene_names)))), dtype=np.float64)
        obs_ntc_mean = np.asarray(b["obs_ntc_mean_official_median_log1p"], dtype=np.float64)
        seen_s = tuple(b.get("seen_sources", ()))
        unseen_s = tuple(b.get("unseen_sources", ()))
        same_s = tuple(b.get("same_lineage_sources", ()))
        extra_fct_sens = b.get("extra_fct_sensitivity")

        p1_idx, p2_idx = make_stratified_half_split_r17(target_names, half_valid, seed=int(seed))

        per_line_results[ln] = {}
        model_exports: Dict[str, Dict[str, Any]] = {}

        for mkey, mdata in sorted(b["models"].items()):
            delta_pred_56 = np.asarray(mdata["delta_pred_56"], dtype=np.float64)
            ntc_pb_ens = np.asarray(mdata["ntc_pb_ensemble_56"], dtype=np.float64)
            ntc_pb_pf = (
                np.asarray(mdata["ntc_pb_per_fold_56"], dtype=np.float64)
                if mdata.get("ntc_pb_per_fold_56") is not None
                else None
            )
            fold_names = list(mdata.get("fold_names", [f"zeroshot/{f}" for f in STATE_CHECKPOINT_FOLDS]))
            per_fold_r56 = (
                np.asarray(mdata["per_fold_r_orth_cal_56"], dtype=np.float64)
                if mdata.get("per_fold_r_orth_cal_56") is not None
                else None
            )
            medoid_idx = mdata.get("medoid_batch_idx")
            sec_comps: Dict[str, np.ndarray] = {}
            for cname in ("bypass_zero", "mean_batch_emb", "multi_batch_ens8", "medoid_idx"):
                arr_c = mdata.get(f"delta_pred_{cname}")
                if arr_c is not None:
                    sec_comps[cname] = np.asarray(arr_c, dtype=np.float64)

            mis_map = mdata.get("mismatched_by_source_and_idx")
            mres = evaluate_model_on_line_r17(
                line_name=ln,
                model_key=mkey,
                delta_obs=delta_obs,
                delta_half_a=delta_half_a,
                delta_half_b=delta_half_b,
                half_valid=half_valid,
                delta_transfer=delta_transfer,
                source_grand_mean=source_grand_mean,
                target_indices=target_indices,
                target_names=target_names,
                fct_delta_pred=fct_delta_pred,
                delta_pred_56=delta_pred_56,
                ntc_pb_ensemble_56=ntc_pb_ens,
                obs_ntc_mean_eval=obs_ntc_mean,
                ntc_pb_per_fold_56=ntc_pb_pf,
                fold_names=fold_names,
                per_fold_r_orth_cal_56=per_fold_r56,
                secondary_comparators=sec_comps,
                medoid_batch_idx=medoid_idx,
                mismatched_by_source_and_idx=mis_map,
                seen_sources=seen_s,
                unseen_sources=unseen_s,
                same_lineage_sources=same_s,
                extra_fct_sensitivity=extra_fct_sens,
                n_boot=n_boot,
                n_boot_bonf56=n_boot_bonf56,
                seed=seed,
                n_splits_diag=n_splits_diag,
            )
            per_line_results[ln][mkey] = mres

            if arrays_dir is not None:
                b_ntc = int(mres["l2_ntc_selection"]["selected_batch_idx"])
                sw_payload: Dict[str, Any] = {
                    "delta_pred_56": delta_pred_56,
                    "batch_indices": np.arange(N_BATCH_INDICES_STATE, dtype=np.int64),
                    "ntc_pb_ensemble_56": ntc_pb_ens,
                    "fold_names": np.asarray(fold_names, dtype=str),
                    "delta_pred_L1_ens56": np.mean(delta_pred_56, axis=0),
                    "delta_pred_L2_ntc_recon_match": delta_pred_56[b_ntc],
                    "l2_selected_batch_idx": np.int64(b_ntc),
                    "delta_pred_multi_batch_ens8": sec_comps.get(
                        "multi_batch_ens8", np.mean(delta_pred_56[list(MULTI_BATCH_ENS8_INDICES)], axis=0)
                    ),
                    "target_names": np.asarray(target_names, dtype=str),
                    "gene_names": np.asarray(gene_names, dtype=str),
                    "model_key": np.asarray(mkey),
                    "checkpoint_id": np.asarray(str(b.get("checkpoint_id", "zeroshot/ensemble4"))),
                    "hvg_norm_mode": np.asarray("official_median_log1p" if mkey == "STATE_HVG" else ""),
                    "control_source": np.asarray(ln),
                    "control_barcodes": np.asarray(list(b.get("matched_ctrl_barcodes", [])), dtype=str),
                }
                if ntc_pb_pf is not None:
                    sw_payload["ntc_pb_per_fold_56"] = ntc_pb_pf
                if per_fold_r56 is not None:
                    sw_payload["per_fold_r_orth_cal_56"] = per_fold_r56
                if "bypass_zero" in sec_comps:
                    sw_payload["delta_pred_bypass_zero"] = sec_comps["bypass_zero"]
                if "mean_batch_emb" in sec_comps:
                    sw_payload["delta_pred_mean_batch_emb"] = sec_comps["mean_batch_emb"]
                if medoid_idx is not None:
                    sw_payload["medoid_batch_idx"] = np.int64(int(medoid_idx))
                    sw_payload["delta_pred_medoid_idx"] = sec_comps.get("medoid_idx", delta_pred_56[int(medoid_idx)])

                g4_payload: Optional[Dict[str, Any]] = None
                if mis_map:
                    src_sorted = sorted(mis_map.keys())
                    b_sorted = sorted({int(k) for s_d in mis_map.values() for k in s_d.keys()})
                    stacked = np.empty(
                        (len(src_sorted), len(b_sorted), len(target_names), len(gene_names)), dtype=np.float64
                    )
                    for si, sk in enumerate(src_sorted):
                        for bi, bv in enumerate(b_sorted):
                            stacked[si, bi] = np.asarray(mis_map[sk][bv], dtype=np.float64)
                    cf_s0 = mres["cross_fitting"]["split_s0"]
                    g4_payload = {
                        "mismatched_sources": np.asarray(src_sorted, dtype=str),
                        "seen_sources": np.asarray(list(seen_s), dtype=str),
                        "unseen_sources": np.asarray(list(unseen_s), dtype=str),
                        "same_lineage_sources": np.asarray(list(same_s), dtype=str),
                        "g4_batch_indices": np.asarray(b_sorted, dtype=np.int64),
                        "delta_pred_mismatched": stacked,
                        "b_star_1": np.int64(int(cf_s0["b_star_1"])),
                        "b_dag_1": np.int64(int(cf_s0["b_dag_1"])),
                        "b_star_2": np.int64(int(cf_s0["b_star_2"])),
                        "b_dag_2": np.int64(int(cf_s0["b_dag_2"])),
                        "target_names": np.asarray(target_names, dtype=str),
                        "gene_names": np.asarray(gene_names, dtype=str),
                        "model_key": np.asarray(mkey),
                        "checkpoint_id": np.asarray(str(b.get("checkpoint_id", "zeroshot/ensemble4"))),
                    }
                model_exports[mkey] = {"sweep56": sw_payload, "g4_mismatched": g4_payload}

        if arrays_dir is not None:
            ctrl_a = b.get("ctrl_audit") or {}
            obs_payload = {
                "line_name": np.asarray(ln),
                "role": np.asarray(role),
                "modality": np.asarray(str(b.get("modality", ""))),
                "checkpoint_id": np.asarray(str(b.get("checkpoint_id", "zeroshot/ensemble4"))),
                "delta_obs": delta_obs,
                "delta_half_a": delta_half_a,
                "delta_half_b": delta_half_b,
                "half_valid": half_valid,
                "delta_transfer": delta_transfer,
                "source_grand_mean": source_grand_mean,
                "fct_delta_pred": fct_delta_pred,
                "e1_constmed_delta_pred": e1_pred,
                "target_indices": target_indices,
                "target_names": np.asarray(target_names, dtype=str),
                "gene_names": np.asarray(gene_names, dtype=str),
                "mu_target": mu_target,
                "mu_sources": mu_sources,
                "obs_ntc_mean_official_median_log1p": obs_ntc_mean,
                "split_s0_p1_indices": np.asarray(p1_idx, dtype=np.int64),
                "split_s0_p2_indices": np.asarray(p2_idx, dtype=np.int64),
                "seen_sources": np.asarray(list(seen_s), dtype=str),
                "unseen_sources": np.asarray(list(unseen_s), dtype=str),
                "same_lineage_sources": np.asarray(list(same_s), dtype=str),
                "matched_ctrl_barcodes": np.asarray(list(b.get("matched_ctrl_barcodes", [])), dtype=str),
                "matched_ctrl_batch_labels": np.asarray(list(b.get("matched_ctrl_batch_labels", [])), dtype=str),
                "n_ctrl_total": np.int64(ctrl_a.get("n_ctrl_total", -1)),
                "n_ctrl_sampled": np.int64(ctrl_a.get("n_ctrl_sampled", -1)),
                "ctrl_sampling_fraction": np.float64(ctrl_a.get("ctrl_sampling_fraction", float("nan"))),
            }
            line_manifests[ln] = save_line_arrays_r17(
                arrays_dir=arrays_dir,
                line_name=ln,
                observed_payload=obs_payload,
                model_exports=model_exports,
                extra_fct_sensitivity=extra_fct_sens,
            )

        line_meta_out[ln] = {
            "role": role,
            "modality": str(b.get("modality", "")),
            "checkpoint_id": str(b.get("checkpoint_id", "zeroshot/ensemble4")),
            "n_eval_targets": len(target_names),
            "n_eval_genes": len(gene_names),
            "n_half_valid": int(np.sum(half_valid)),
            "n_self_masked": int(np.sum(target_indices >= 0)),
            "strata_sources": {
                "seen_by_checkpoint": list(seen_s),
                "unseen_by_checkpoint": list(unseen_s),
                "same_lineage_comparators": list(same_s),
            },
            "disclosures": list(b.get("disclosures", [])),
            "ctrl_audit": b.get("ctrl_audit", {}),
            "embedding_diagnostics": b.get("embedding_diagnostics", {}),
            "array_manifest_sha256": line_manifests.get(ln, {}),
        }

    cohort = adjudicate_cohort_r17(
        per_line_results=per_line_results,
        line_roles=line_roles,
        line_n_eval=line_n_eval,
        n_perm=n_perm,
        seed=seed,
    )

    lines_combined: Dict[str, Any] = {}
    for ln in sorted(line_bundles.keys()):
        lines_combined[ln] = {
            **line_meta_out[ln],
            "models": cohort["lines_clean"][ln],
        }

    return {
        "protocol": "PREREG_J_v2",
        "seed": int(seed),
        "n_boot": int(n_boot),
        "n_boot_bonf56": int(n_boot_bonf56),
        "n_perm": int(n_perm),
        "n_splits_diag": int(n_splits_diag),
        "delta_practical": float(DEFAULT_PRACTICAL_DELTA),
        "K_conf": cohort["K_conf"],
        "confirmatory_units": cohort["confirmatory_units"],
        "separate_model_holm_fwer05": cohort["separate_model_holm_fwer05"],
        "pooled_model_holm_sensitivity_fwer05": cohort["pooled_model_holm_sensitivity_fwer05"],
        "j5b_cross_line_rank_correlations": cohort["j5b_cross_line_rank_correlations"],
        "parametric_power_fill_in_s43": cohort["parametric_power_fill_in_s43"],
        "lines": lines_combined,
    }


def recompute_cohort_from_exported_arrays_r17(
    arrays_dir: Path,
    n_boot: int = DEFAULT_R17_N_BOOT,
    n_boot_bonf56: int = DEFAULT_R17_N_BOOT_BONF56,
    n_perm: int = DEFAULT_R17_N_PERM,
    seed: int = DEFAULT_R17_SEED,
    n_splits_diag: int = DEFAULT_R17_N_SPLITS_DIAG,
) -> Dict[str, Any]:
    """
    Stage 4 Controller Verification Helper (`PREREG_J §5.2`):
    Reconstruct and recompute EVERY statistic, CI, Holm step-down call, `J5a` `G4` decomposition,
    and `J5b` Spearman + permutation p-value STRICTLY from saved `.npz` files under `arrays_dir`.
    """
    adir = Path(arrays_dir)
    line_dirs = sorted(p for p in adir.iterdir() if p.is_dir() and (p / "observed.npz").is_file())
    if not line_dirs:
        raise FileNotFoundError(f"No exported line directories with observed.npz found under {adir}")

    bundles: Dict[str, Dict[str, Any]] = {}
    for ldir in line_dirs:
        ln = ldir.name
        with np.load(ldir / "observed.npz", allow_pickle=False) as obs:
            entry: Dict[str, Any] = {
                "role": str(obs["role"].item()),
                "modality": str(obs["modality"].item()),
                "checkpoint_id": str(obs["checkpoint_id"].item()),
                "delta_obs": obs["delta_obs"].copy(),
                "delta_half_a": obs["delta_half_a"].copy(),
                "delta_half_b": obs["delta_half_b"].copy(),
                "half_valid": obs["half_valid"].copy(),
                "delta_transfer": obs["delta_transfer"].copy(),
                "source_grand_mean": obs["source_grand_mean"].copy(),
                "fct_delta_pred": obs["fct_delta_pred"].copy(),
                "e1_constmed_delta_pred": obs["e1_constmed_delta_pred"].copy(),
                "target_indices": obs["target_indices"].copy(),
                "target_names": [str(x) for x in obs["target_names"].tolist()],
                "gene_names": [str(x) for x in obs["gene_names"].tolist()],
                "mu_target": obs["mu_target"].copy(),
                "mu_sources": obs["mu_sources"].copy(),
                "obs_ntc_mean_official_median_log1p": obs["obs_ntc_mean_official_median_log1p"].copy(),
                "seen_sources": [str(x) for x in obs["seen_sources"].tolist()],
                "unseen_sources": [str(x) for x in obs["unseen_sources"].tolist()],
                "same_lineage_sources": [str(x) for x in obs["same_lineage_sources"].tolist()],
                "matched_ctrl_barcodes": [str(x) for x in obs["matched_ctrl_barcodes"].tolist()],
                "matched_ctrl_batch_labels": [str(x) for x in obs["matched_ctrl_batch_labels"].tolist()],
                "models": {},
            }

        # Optional sensitivity FCT files (`SENS_*.npz`)
        sens_files = sorted(ldir.glob("SENS_*.npz"))
        if sens_files:
            sens_dict: Dict[str, Dict[str, np.ndarray]] = {}
            for sf in sens_files:
                skey = sf.stem[len("SENS_") :]
                with np.load(sf, allow_pickle=False) as sd:
                    sens_dict[skey] = {
                        "fct_delta_pred": sd["fct_delta_pred"].copy(),
                        "delta_transfer": sd["delta_transfer"].copy(),
                        "source_grand_mean": sd["source_grand_mean"].copy(),
                    }
            entry["extra_fct_sensitivity"] = sens_dict

        for sw_file in sorted(ldir.glob("preds_*_sweep56.npz")):
            mkey = sw_file.stem[len("preds_") : -len("_sweep56")]
            with np.load(sw_file, allow_pickle=False) as sw:
                m_entry: Dict[str, Any] = {
                    "delta_pred_56": sw["delta_pred_56"].copy(),
                    "ntc_pb_ensemble_56": sw["ntc_pb_ensemble_56"].copy(),
                    "fold_names": [str(x) for x in sw["fold_names"].tolist()],
                }
                if "ntc_pb_per_fold_56" in sw:
                    m_entry["ntc_pb_per_fold_56"] = sw["ntc_pb_per_fold_56"].copy()
                if "per_fold_r_orth_cal_56" in sw:
                    m_entry["per_fold_r_orth_cal_56"] = sw["per_fold_r_orth_cal_56"].copy()
                for cname in ("bypass_zero", "mean_batch_emb", "multi_batch_ens8", "medoid_idx"):
                    k = f"delta_pred_{cname}"
                    if k in sw:
                        m_entry[k] = sw[k].copy()
                if "medoid_batch_idx" in sw:
                    m_entry["medoid_batch_idx"] = int(sw["medoid_batch_idx"].item())

            g4_file = ldir / f"preds_{mkey}_g4_mismatched.npz"
            if g4_file.is_file():
                with np.load(g4_file, allow_pickle=False) as g4:
                    srcs = [str(x) for x in g4["mismatched_sources"].tolist()]
                    b_idxs = [int(x) for x in g4["g4_batch_indices"].tolist()]
                    d_mis = g4["delta_pred_mismatched"]
                    mis_map: Dict[str, Dict[int, np.ndarray]] = {}
                    for si, sk in enumerate(srcs):
                        mis_map[sk] = {bv: d_mis[si, bi].copy() for bi, bv in enumerate(b_idxs)}
                    m_entry["mismatched_by_source_and_idx"] = mis_map

            entry["models"][mkey] = m_entry

        bundles[ln] = entry

    return evaluate_cohort_from_arrays_r17(
        line_bundles=bundles,
        arrays_dir=None,
        n_boot=n_boot,
        n_boot_bonf56=n_boot_bonf56,
        n_perm=n_perm,
        seed=seed,
        n_splits_diag=n_splits_diag,
    )


# ==============================================================================
# Real-Data Orchestrator (Strictly Gated on `R17_CODE_FREEZE.txt` + `R17_GO.txt`)
# ==============================================================================


def run_real_evaluation_r17(
    lines_to_run: Sequence[str],
    status_dir: Path = DEFAULT_R17_STATUS_DIR,
    prereg_lock_path: Optional[Path] = DEFAULT_R17_PREREG_LOCK_PATH,
    models_dir: Path = DEFAULT_MODELS_DIR,
    results_dir: Path = DEFAULT_R17_RESULTS_DIR,
    device: str = "mps",
    n_boot: int = DEFAULT_R17_N_BOOT,
    n_boot_bonf56: int = DEFAULT_R17_N_BOOT_BONF56,
    n_perm: int = DEFAULT_R17_N_PERM,
    seed: int = DEFAULT_R17_SEED,
    ctrl_seed: int = DEFAULT_R17_SEED,
    ctrl_resample_seed: int = DEFAULT_R17_RESAMPLE_SEED,
    n_splits_diag: int = DEFAULT_R17_N_SPLITS_DIAG,
    gse295214_raw_dir: Optional[Path] = None,
    pfizer_raw_path: Optional[Path] = None,
    raw_nadig_dir: Optional[Path] = None,
    include_p2_to_p1_g4: bool = True,
) -> Dict[str, Any]:
    """
    Real-cell Round 17 evaluation orchestrator (`PREREG_J` v2).
    STRICTLY refuses execution unless `R17_CODE_FREEZE.txt` and `R17_GO.txt` (plus
    `PREREG_J_LOCK.md`) exist and are non-empty.
    """
    gate_paths = assert_all_execution_gates_r17(status_dir=status_dir, prereg_lock_path=prereg_lock_path)
    target_lines = validate_target_lines_r17(lines_to_run)

    os.environ.setdefault("WANDB_MODE", "disabled")
    os.environ.setdefault("HF_HUB_OFFLINE", "1")
    os.environ.setdefault("TRANSFORMERS_OFFLINE", "1")

    t0 = time.time()
    bridge = HGNCBridge.from_tsv(DEFAULT_HGNC_FILE)
    common_genes_r9 = load_r9_common_genes(hgnc_bridge=bridge, common_genes_file=R9_COMMON_GENES_FILE)
    g_model_2000, p_model_canon_to_key = load_state_vocab_r16(models_dir, bridge)
    gene_canon_map = {g: canonicalize_gene_r17(bridge, g) for g in g_model_2000}

    # Determine all unique control pools needed across requested lines + their J5a mismatched strata
    strata_by_line: Dict[str, Dict[str, Any]] = {}
    ctrl_pool_names: Set[str] = set()
    for ln in target_lines:
        st = build_mismatched_strata_r17(ln)
        strata_by_line[ln] = st
        ctrl_pool_names.add(ln)
        ctrl_pool_names.update(st["seen_by_checkpoint"])
        ctrl_pool_names.update(st["unseen_by_checkpoint"])
        ctrl_pool_names.update(st["same_lineage_comparators"])

    # Step 1: Extract raw NTC counts + barcodes ONCE per pool & cache SE-600M (`X_state`) + ST-HVG (`X_hvg`)
    ctrl_cache: Dict[str, Dict[str, Any]] = {}
    control_source_barcodes: Dict[str, Dict[str, Any]] = {}
    control_embeddings_to_save: Dict[str, Dict[str, Any]] = {}

    for pool_name in sorted(ctrl_pool_names):
        bundle = extract_control_ntc_bundle_r17(
            pool_name=pool_name,
            status_dir=status_dir,
            prereg_lock_path=prereg_lock_path,
            seed=int(ctrl_seed),
            resample_seed=int(ctrl_resample_seed),
            include_resample=False,
            gse295214_raw_dir=gse295214_raw_dir,
            pfizer_raw_path=pfizer_raw_path,
            raw_nadig_dir=raw_nadig_dir,
        )
        raw_mat = bundle["raw_counts"]
        raw_genes = bundle["gene_names"]
        batch_labels = bundle["batch_labels"]
        barcodes = bundle["barcodes"]
        c_audit = bundle["ctrl_audit"]

        x_state, se_audit = encode_controls_se600m_r16(
            raw_counts=raw_mat,
            raw_gene_names=raw_genes,
            se_dir=Path(models_dir) / "SE-600M",
            device=device,
            batch_size=16,
            dataset_name=pool_name,
            seed=int(ctrl_seed),
        )
        control_embeddings_to_save[pool_name] = {
            "X_state": x_state,
            "control_barcodes": list(barcodes),
            "seed": int(ctrl_seed),
            "dataset_name": pool_name,
            "se_audit": se_audit,
        }

        raw_canon = {g: canonicalize_gene_r17(bridge, g) for g in raw_genes}
        full_canon = {**raw_canon, **gene_canon_map}
        x_hvg_official, hvg_off_audit = preprocess_hvg_counts_r16(
            raw_mat, raw_genes, g_model_2000, norm_mode="official_median_log1p", canon_map=full_canon
        )
        # Observed unlabeled NTC mean expression in official_median_log1p space on ordered G_model (2000)
        obs_ntc_mean_g_model = np.mean(x_hvg_official.astype(np.float64), axis=0)

        ctrl_cache[pool_name] = {
            "x_state": x_state,
            "x_hvg_official": x_hvg_official,
            "obs_ntc_mean_g_model": obs_ntc_mean_g_model,
            "batch_labels": batch_labels,
            "barcodes": barcodes,
            "ctrl_audit": c_audit,
            "se_audit": se_audit,
            "hvg_official_audit": hvg_off_audit,
        }
        control_source_barcodes[pool_name] = {
            "ctrl_seed": int(ctrl_seed),
            "n_ctrl_total": int(c_audit["n_ctrl_total"]),
            "n_ctrl_sampled": int(c_audit["n_ctrl_sampled"]),
            "barcodes": list(barcodes),
            "batch_labels": list(batch_labels),
        }
        del raw_mat, bundle
        gc.collect()

    rdir = Path(results_dir)
    rdir.mkdir(parents=True, exist_ok=True)
    barcodes_json_path = rdir / "control_barcodes_r17.json"
    barcodes_json_path.write_text(
        json.dumps(control_source_barcodes, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    embeddings_dir = rdir / "embeddings"
    embeddings_manifest = save_control_embeddings_r16(embeddings_dir, control_embeddings_to_save)

    g_model_canon_to_idx: Dict[str, int] = {}
    for idx, gm in enumerate(g_model_2000):
        g_model_canon_to_idx.setdefault(gene_canon_map[gm], idx)

    # Step 2: Phase 1 (56-index sweep + neutral comparators) & Phase 2 (J5a G4 mismatched) per line & model
    line_bundles: Dict[str, Dict[str, Any]] = {}
    for ln in target_lines:
        spec = LINE_SPECS_R17[ln]
        st = strata_by_line[ln]
        seen_s = st["seen_by_checkpoint"]
        unseen_s = st["unseen_by_checkpoint"]
        same_s = st["same_lineage_comparators"]
        all_mis = sorted(set(seen_s) | set(unseen_s) | set(same_s))

        aligned = align_line_and_sources_r17(
            line_name=ln,
            fct_source_names=spec["fct_sources"],
            line_dir=spec["line_dir"],
            g_model_2000=g_model_2000,
            p_model_canon_to_key=p_model_canon_to_key,
            hgnc_bridge=bridge,
            common_genes_r9=common_genes_r9,
        )
        q_names = aligned["target_names"]
        g_eval = aligned["gene_names"]
        g_eval_in_model_idx = np.asarray([g_model_canon_to_idx[g] for g in g_eval], dtype=np.int64)
        obs_ntc_mean_eval = ctrl_cache[ln]["obs_ntc_mean_g_model"][g_eval_in_model_idx]

        extra_fct_sens: Optional[Dict[str, Dict[str, np.ndarray]]] = None
        if "fct_sources_sensitivity" in spec:
            sens_al = align_line_and_sources_r17(
                line_name=ln,
                fct_source_names=spec["fct_sources_sensitivity"],
                line_dir=spec["line_dir"],
                g_model_2000=g_model_2000,
                p_model_canon_to_key=p_model_canon_to_key,
                hgnc_bridge=bridge,
                common_genes_r9=common_genes_r9,
                fixed_target_names=q_names,
                fixed_gene_names=g_eval,
            )
            extra_fct_sens = {
                "FCT_full_r9_6_incl_vcc2025_h1": {
                    "fct_delta_pred": sens_al["fct_delta_pred"],
                    "delta_transfer": sens_al["delta_transfer"],
                    "source_grand_mean": sens_al["source_grand_mean"],
                }
            }

        emb_diag = {
            "matched_se_audit": ctrl_cache[ln]["se_audit"],
            "vs_mismatched_centroid_separation": {
                sk: compute_embedding_separation_diagnostics(ctrl_cache[ln]["x_state"], ctrl_cache[sk]["x_state"])
                for sk in all_mis
            },
        }

        models_out: Dict[str, Dict[str, Any]] = {}
        for mkey in MODELS_R17:
            ctrl_key = "x_state" if mkey == "STATE_SE" else "x_hvg_official"
            sw_res = run_state_sweep56_r17(
                variant=mkey,
                target_line=ln,
                ctrl_tensor_np=ctrl_cache[ln][ctrl_key],
                pert_names=q_names,
                gene_names=g_eval,
                pert_canon_to_model=p_model_canon_to_key,
                gene_canon_map=gene_canon_map,
                models_dir=models_dir,
                status_dir=status_dir,
                prereg_lock_path=prereg_lock_path,
                device=device,
                seed=int(seed),
                obs_for_per_fold_diag=aligned,
            )

            # Determine cross-fitted (b*_1, b_dag_1) [and optional (b*_2, b_dag_2)] on locked split s=0
            moms_tmp = compute_sweep56_moments_r17(
                delta_obs=aligned["delta_obs"],
                delta_pred_56=sw_res["delta_pred_56"],
                fct_delta_pred=aligned["fct_delta_pred"],
                delta_transfer=aligned["delta_transfer"],
                target_gene_idx=aligned["target_indices"],
                source_grand_mean=None,
                delta_half_a=aligned["delta_half_a"],
                delta_half_b=aligned["delta_half_b"],
                half_valid=aligned["half_valid"],
            )
            cf_tmp = compute_cross_fitted_spread_r17(
                moments_bundle=moms_tmp,
                target_names=q_names,
                half_valid=aligned["half_valid"],
                n_boot=2,  # Only need exact deterministic split s=0 selected indices here
                seed=int(seed),
                n_splits_diag=1,
            )
            s0_tmp = cf_tmp["split_s0"]
            g4_needed = {0, 17, 35, int(s0_tmp["b_star_1"]), int(s0_tmp["b_dag_1"])}
            if include_p2_to_p1_g4:
                g4_needed.update({int(s0_tmp["b_star_2"]), int(s0_tmp["b_dag_2"])})

            g4_res = run_state_g4_mismatched_r17(
                variant=mkey,
                target_line=ln,
                mismatched_ctrl_tensors={sk: ctrl_cache[sk][ctrl_key] for sk in all_mis},
                g4_batch_indices=sorted(g4_needed),
                pert_names=q_names,
                gene_names=g_eval,
                pert_canon_to_model=p_model_canon_to_key,
                gene_canon_map=gene_canon_map,
                models_dir=models_dir,
                status_dir=status_dir,
                prereg_lock_path=prereg_lock_path,
                device=device,
                seed=int(seed),
            )

            models_out[mkey] = {
                "delta_pred_56": sw_res["delta_pred_56"],
                "ntc_pb_ensemble_56": sw_res["ntc_pb_ensemble_56"],
                "ntc_pb_per_fold_56": sw_res["ntc_pb_per_fold_56"],
                "fold_names": sw_res["fold_names"],
                "per_fold_r_orth_cal_56": sw_res["per_fold_r_orth_cal_56"],
                "delta_pred_bypass_zero": sw_res["delta_pred_bypass_zero"],
                "delta_pred_mean_batch_emb": sw_res["delta_pred_mean_batch_emb"],
                "delta_pred_multi_batch_ens8": sw_res["delta_pred_multi_batch_ens8"],
                "delta_pred_medoid_idx": sw_res["delta_pred_medoid_idx"],
                "medoid_batch_idx": int(sw_res["medoid_info"]["medoid_idx"]),
                "mismatched_by_source_and_idx": g4_res["by_source_and_idx"],
            }

        line_bundles[ln] = {
            "role": spec["role"],
            "modality": spec["modality"],
            "checkpoint_id": spec["checkpoint_id"],
            "disclosures": spec["disclosures"],
            "seen_sources": seen_s,
            "unseen_sources": unseen_s,
            "same_lineage_sources": same_s,
            "obs_ntc_mean_official_median_log1p": obs_ntc_mean_eval,
            "matched_ctrl_barcodes": ctrl_cache[ln]["barcodes"],
            "matched_ctrl_batch_labels": ctrl_cache[ln]["batch_labels"],
            "ctrl_audit": ctrl_cache[ln]["ctrl_audit"],
            "embedding_diagnostics": emb_diag,
            "extra_fct_sensitivity": extra_fct_sens,
            "models": models_out,
            **aligned,
        }

    arrays_dir = rdir / "arrays"
    summary = evaluate_cohort_from_arrays_r17(
        line_bundles=line_bundles,
        arrays_dir=arrays_dir,
        n_boot=n_boot,
        n_boot_bonf56=n_boot_bonf56,
        n_perm=n_perm,
        seed=seed,
        n_splits_diag=n_splits_diag,
    )
    summary.update(
        {
            "ctrl_seed": int(ctrl_seed),
            "ctrl_resample_seed": int(ctrl_resample_seed),
            "gates_verified": {k: str(v) for k, v in gate_paths.items()},
            "elapsed_sec": round(float(time.time() - t0), 2),
            "control_barcodes_file": str(barcodes_json_path),
            "embeddings_dir": str(embeddings_dir),
            "embeddings_manifest_sha256": embeddings_manifest,
        }
    )
    (rdir / "summary_r17.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


# ==============================================================================
# Pure-Synthetic End-to-End CPU Smoke Runner (ZERO Real Cells Touched)
# ==============================================================================


def _make_synthetic_line_bundle_r17(
    line_name: str,
    role: str,
    modality: str,
    n_perts: int,
    n_genes: int,
    regime_se: str,
    regime_hvg: str,
    seed: int,
    frac_half_valid: float = 1.0,
    include_fct_sensitivity: bool = False,
) -> Dict[str, Any]:
    """
    Generate one synthetic line across regimes:
      - `strong_batch`: strong true 56-index variation + sign reversal vs FCT
      - `invariant_nonzero`: positive detectable signal with near-zero between-index spread
      - `low_signal`: near-zero signal floor triggering `J2` `INCONCLUSIVE_LOW_SIGNAL`
    """
    rng = np.random.default_rng(int(seed))
    P, G = int(n_perts), int(n_genes)
    gene_names = [f"SYN_G{j:04d}" for j in range(G)]
    target_names = [f"SYN_G{i:04d}" for i in range(P)]
    target_indices = np.arange(P, dtype=np.int64)
    half_valid = np.ones(P, dtype=bool)
    if frac_half_valid < 1.0 and P >= 4:
        n_false = max(1, int(round(P * (1.0 - frac_half_valid))))
        half_valid[-n_false:] = False

    b_s = rng.normal(0.0, 0.12, size=G)
    tm = b_s[None, :] + rng.normal(0.0, 0.16, size=(P, G))
    true_spec = rng.normal(0.0, 0.30, size=(P, G))
    ha = tm + true_spec + rng.normal(0.0, 0.05, size=(P, G))
    hb = tm + true_spec + rng.normal(0.0, 0.05, size=(P, G))
    obs = 0.5 * (ha + hb)

    mu_sources = np.exp(rng.normal(0.0, 0.7, size=(5, G)))
    mu_target = np.exp(rng.normal(0.2, 0.7, size=G))
    fc_bar = np.exp(rng.normal(0.0, 0.25, size=(P, G)))
    has_hits = np.ones(P, dtype=bool)
    # Give FCT a moderate positive residual correlation (~0.25) so H-J2 sign reversals can occur cleanly
    fct_pred = fct_from_fc_bar(mu_target, fc_bar, tm, has_hits) + 0.22 * true_spec
    e1_pred = fct_from_fc_bar(np.full(G, float(np.median(mu_target))), fc_bar, tm, has_hits)

    obs_ntc_mean = np.abs(rng.normal(1.2, 0.6, size=G))
    st = build_mismatched_strata_r17(line_name) if line_name in LINE_SPECS_R17 else {
        "seen_by_checkpoint": ("nadig_hepg2", "nadig_jurkat"),
        "unseen_by_checkpoint": ("a549_ko", "dld1", "feng_ipsc"),
        "same_lineage_comparators": (),
    }
    all_mis = sorted(set(st["seen_by_checkpoint"]) | set(st["unseen_by_checkpoint"]) | set(st["same_lineage_comparators"]))

    def _gen_model(regime: str, m_offset: int) -> Dict[str, Any]:
        mrng = np.random.default_rng(int(seed) + m_offset)
        p56 = np.empty((N_BATCH_INDICES_STATE, P, G), dtype=np.float64)
        ntc_pf = np.empty((4, N_BATCH_INDICES_STATE, G), dtype=np.float64)

        if regime == "strong_batch":
            # Smooth index quality gradient from -0.02 (b=0..9) up to +0.55 (b=45..55)
            gains = np.linspace(-0.02, 0.55, N_BATCH_INDICES_STATE)
            # Keep b=0 below FCT so b=0 -> FCT_SUPERIOR and high b -> MODEL_SUPERIOR (full sign reversal)
            gains[0] = 0.04
            gains[35] = 0.48
        elif regime == "invariant_nonzero":
            gains = np.full(N_BATCH_INDICES_STATE, 0.45, dtype=np.float64)
        elif regime == "low_signal":
            gains = np.full(N_BATCH_INDICES_STATE, 0.0, dtype=np.float64)
        else:
            raise ValueError(f"Unknown synthetic regime: {regime!r}")

        shared_noise = mrng.normal(0.0, 0.20, size=(P, G))
        for b in range(N_BATCH_INDICES_STATE):
            if regime == "invariant_nonzero":
                p56[b] = tm + gains[b] * true_spec + 0.35 * shared_noise + mrng.normal(0.0, 0.002, size=(P, G))
            elif regime == "low_signal":
                p56[b] = tm + mrng.normal(0.0, 0.25, size=(P, G))
            else:
                p56[b] = tm + gains[b] * true_spec + (0.35 - 0.25 * max(gains[b], 0.0)) * shared_noise + mrng.normal(0.0, 0.04, size=(P, G))

            # Construct synthetic NTC reconstructions where high-gain index 42 aligns best with obs_ntc_mean
            align_w = 0.96 if (regime == "strong_batch" and b == 42) else (0.5 + 0.3 * (b / 55.0))
            for fi in range(4):
                ntc_pf[fi, b] = np.clip(
                    align_w * obs_ntc_mean + (1.0 - align_w) * mrng.normal(1.2, 0.6, size=G) + mrng.normal(0.0, 0.03, size=G),
                    0.0,
                    14.0,
                )

        ntc_ens = np.mean(ntc_pf, axis=0)
        # Find exact locked s=0 split indices so we only synthesize mismatched arrays for needed G4 indices
        moms_t = compute_sweep56_moments_r17(
            obs, p56, fct_pred, tm, target_indices, None, ha, hb, half_valid
        )
        cf_t = compute_cross_fitted_spread_r17(moms_t, target_names, half_valid, n_boot=2, seed=int(DEFAULT_R17_SEED), n_splits_diag=1)
        s0_t = cf_t["split_s0"]
        g4_idxs = sorted({0, 17, 35, int(s0_t["b_star_1"]), int(s0_t["b_dag_1"]), int(s0_t["b_star_2"]), int(s0_t["b_dag_2"])})
        mis_map: Dict[str, Dict[int, np.ndarray]] = {}
        for sk in all_mis:
            mis_map[sk] = {}
            for bv in g4_idxs:
                mis_map[sk][bv] = p56[bv] + mrng.normal(0.0, 0.03, size=(P, G))

        # Per-fold 56-index synthetic diagnostic
        pf_r = np.tile(np.clip(gains, -1.0, 1.0)[None, :], (4, 1)) + mrng.normal(0.0, 0.005, size=(4, N_BATCH_INDICES_STATE))

        return {
            "delta_pred_56": p56,
            "ntc_pb_ensemble_56": ntc_ens,
            "ntc_pb_per_fold_56": ntc_pf,
            "fold_names": [f"zeroshot/{f}" for f in STATE_CHECKPOINT_FOLDS],
            "per_fold_r_orth_cal_56": pf_r,
            "delta_pred_bypass_zero": np.mean(p56, axis=0) + mrng.normal(0.0, 0.01, size=(P, G)),
            "delta_pred_mean_batch_emb": np.mean(p56, axis=0) + mrng.normal(0.0, 0.015, size=(P, G)),
            "delta_pred_multi_batch_ens8": np.mean(p56[list(MULTI_BATCH_ENS8_INDICES)], axis=0),
            "delta_pred_medoid_idx": p56[12].copy(),
            "medoid_batch_idx": 12,
            "mismatched_by_source_and_idx": mis_map,
        }

    extra_sens = None
    if include_fct_sensitivity:
        sens_tm = tm + rng.normal(0.0, 0.01, size=(P, G))
        sens_fct = fct_pred + rng.normal(0.0, 0.01, size=(P, G))
        extra_sens = {
            "FCT_full_r9_6_incl_vcc2025_h1": {
                "fct_delta_pred": sens_fct,
                "delta_transfer": sens_tm,
                "source_grand_mean": b_s + rng.normal(0.0, 0.005, size=G),
            }
        }

    bcs = [f"SYN_{line_name.upper()}_{i:04d}" for i in range(512)]
    return {
        "role": role,
        "modality": modality,
        "checkpoint_id": "zeroshot/ensemble4",
        "disclosures": [f"Synthetic smoke fixture for {line_name} (ZERO real cells touched)."],
        "seen_sources": st["seen_by_checkpoint"],
        "unseen_sources": st["unseen_by_checkpoint"],
        "same_lineage_sources": st["same_lineage_comparators"],
        "gene_names": gene_names,
        "target_names": target_names,
        "target_indices": target_indices,
        "delta_obs": obs,
        "delta_half_a": ha,
        "delta_half_b": hb,
        "half_valid": half_valid,
        "delta_transfer": tm,
        "source_grand_mean": b_s,
        "mu_target": mu_target,
        "mu_sources": mu_sources,
        "fct_delta_pred": fct_pred,
        "e1_constmed_delta_pred": e1_pred,
        "obs_ntc_mean_official_median_log1p": obs_ntc_mean,
        "matched_ctrl_barcodes": bcs,
        "matched_ctrl_batch_labels": [f"batch_{i % 4}" for i in range(512)],
        "ctrl_audit": {"n_ctrl_total": 2048, "n_ctrl_sampled": 512, "ctrl_sampling_fraction": 0.25},
        "extra_fct_sensitivity": extra_sens,
        "models": {
            "STATE_SE": _gen_model(regime_se, 1001),
            "STATE_HVG": _gen_model(regime_hvg, 2002),
        },
    }


def run_synthetic_smoke_r17(
    out_dir: Path,
    n_boot: int = 800,
    n_boot_bonf56: int = 1000,
    n_perm: int = 1000,
    seed: int = DEFAULT_R17_SEED,
) -> Dict[str, Any]:
    """
    Pure-CPU, zero-real-cell end-to-end smoke run exercising EVERY Round 17 component:
      - 3 confirmatory units (`cd4t_stim48` [mixed `half_valid` + `same_lineage_resting_comparator`],
        `gse295214_h1` [with full `r9_6` `FCT` sensitivity], `pfizer_telohaec_basal` [low-signal floor
        on `STATE_HVG` triggering `J2` `INCONCLUSIVE_LOW_SIGNAL`!])
      - 1 within-family replication (`pfizer_telohaec_cytokine`)
      - 1 secondary descriptive (`a549_myllia`, `N=20`)
      - 1 micro-N descriptive (`mac_ctrl`, `N=3`)
      - Full `.npz` array export, `embeddings/` export, `control_barcodes_r17.json`, `summary_r17.json`,
        AND Stage 4 from-disk recomputation parity verification (`max_abs_diff <= 1e-12`).
    """
    t0 = time.time()
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    bundles: Dict[str, Dict[str, Any]] = {
        "cd4t_stim48": _make_synthetic_line_bundle_r17(
            "cd4t_stim48", "confirmatory", "CRISPRi", 60, 120,
            regime_se="strong_batch", regime_hvg="strong_batch", seed=seed + 1, frac_half_valid=0.85,
        ),
        "gse295214_h1": _make_synthetic_line_bundle_r17(
            "gse295214_h1", "confirmatory", "CRISPRi", 64, 120,
            regime_se="strong_batch", regime_hvg="invariant_nonzero", seed=seed + 2, include_fct_sensitivity=True,
        ),
        "pfizer_telohaec_basal": _make_synthetic_line_bundle_r17(
            "pfizer_telohaec_basal", "confirmatory", "CRISPRi", 54, 120,
            regime_se="invariant_nonzero", regime_hvg="low_signal", seed=seed + 3,
        ),
        "pfizer_telohaec_cytokine": _make_synthetic_line_bundle_r17(
            "pfizer_telohaec_cytokine", "within_family_replication", "CRISPRi", 54, 120,
            regime_se="invariant_nonzero", regime_hvg="low_signal", seed=seed + 4,
        ),
        "a549_myllia": _make_synthetic_line_bundle_r17(
            "a549_myllia", "secondary_fresh", "CRISPRko", 20, 100,
            regime_se="strong_batch", regime_hvg="strong_batch", seed=seed + 5,
        ),
        "mac_ctrl": _make_synthetic_line_bundle_r17(
            "mac_ctrl", "micro_n_descriptive", "CRISPRko", 3, 100,
            regime_se="strong_batch", regime_hvg="low_signal", seed=seed + 6,
        ),
    }

    arrays_dir = out_path / "arrays"
    summary = evaluate_cohort_from_arrays_r17(
        line_bundles=bundles,
        arrays_dir=arrays_dir,
        n_boot=n_boot,
        n_boot_bonf56=n_boot_bonf56,
        n_perm=n_perm,
        seed=seed,
        n_splits_diag=DEFAULT_R17_N_SPLITS_DIAG,
    )

    # Export synthetic `2,058`-dim `X_state` embeddings & barcodes
    rng_emb = np.random.default_rng(int(seed))
    syn_embeddings: Dict[str, Dict[str, Any]] = {}
    syn_barcodes: Dict[str, Dict[str, Any]] = {}
    for ln in sorted(bundles.keys()):
        x_arr = rng_emb.normal(0.0, 0.5, size=(512, 2058)).astype(np.float32)
        bcs = bundles[ln]["matched_ctrl_barcodes"]
        syn_embeddings[ln] = {
            "X_state": x_arr,
            "control_barcodes": bcs,
            "seed": int(seed),
            "dataset_name": ln,
            "se_audit": {"device": "cpu_synthetic", "dtype": "float32", "n_cells": 512, "embed_dim": 2058, "all_finite": True},
        }
        syn_barcodes[ln] = {"ctrl_seed": int(seed), "n_ctrl_sampled": 512, "barcodes": bcs}

    emb_manifest = save_control_embeddings_r16(out_path / "embeddings", syn_embeddings)
    (out_path / "control_barcodes_r17.json").write_text(
        json.dumps(syn_barcodes, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    # Stage 4 independent recompute directly from exported `.npz` files!
    recomp = recompute_cohort_from_exported_arrays_r17(
        arrays_dir=arrays_dir,
        n_boot=n_boot,
        n_boot_bonf56=n_boot_bonf56,
        n_perm=n_perm,
        seed=seed,
        n_splits_diag=DEFAULT_R17_N_SPLITS_DIAG,
    )
    max_recomp_diff = _max_numeric_diff(summary, recomp)

    summary["mode"] = "synthetic_smoke"
    summary["real_cells_touched"] = 0
    summary["device"] = "cpu"
    summary["embeddings_manifest_sha256"] = emb_manifest
    summary["stage4_from_disk_recompute_max_abs_diff"] = float(max_recomp_diff)
    summary["stage4_from_disk_recompute_passed"] = bool(np.isfinite(max_recomp_diff) and max_recomp_diff <= 1e-12)
    summary["elapsed_sec"] = round(float(time.time() - t0), 2)
    (out_path / "summary_r17.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def _max_numeric_diff(a: Any, b: Any) -> float:
    """Recursively compute max absolute numeric diff between two matching summary dicts (skipping metadata keys)."""
    skip_keys = {
        "disclosures", "ctrl_audit", "embedding_diagnostics", "array_manifest_sha256",
        "mode", "real_cells_touched", "device", "embeddings_manifest_sha256",
        "stage4_from_disk_recompute_max_abs_diff", "stage4_from_disk_recompute_passed", "elapsed_sec",
    }
    if isinstance(a, dict) and isinstance(b, dict):
        mx = 0.0
        for k in set(a.keys()) & set(b.keys()):
            if k in skip_keys:
                continue
            d = _max_numeric_diff(a[k], b[k])
            if np.isnan(d):
                return float("nan")
            mx = max(mx, d)
        return mx
    if isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        mx = 0.0
        for x, y in zip(a, b):
            d = _max_numeric_diff(x, y)
            if np.isnan(d):
                return float("nan")
            mx = max(mx, d)
        return mx
    if isinstance(a, (int, float)) and not isinstance(a, bool) and isinstance(b, (int, float)) and not isinstance(b, bool):
        fa, fb = float(a), float(b)
        if np.isnan(fa) and np.isnan(fb):
            return 0.0
        if np.isnan(fa) or np.isnan(fb):
            return float("nan")
        return abs(fa - fb)
    if a != b:
        return float("inf")
    return 0.0


# ==============================================================================
# CLI Entrypoint
# ==============================================================================


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Round 17 Evaluation Runner (PREREG_J v2).")
    p.add_argument(
        "--synthetic-smoke",
        action="store_true",
        help="Run end-to-end CPU synthetic smoke + array export + Stage 4 recompute check (zero real cells).",
    )
    p.add_argument(
        "--cohort",
        choices=["confirmatory", "within_family", "descriptive", "subset_sensitivity", "all"],
        default="all",
        help="Which line cohort to run in real-data mode.",
    )
    p.add_argument("--lines", nargs="*", default=None, help="Explicit subset of Round 17 line keys.")
    p.add_argument("--status-dir", type=Path, default=DEFAULT_R17_STATUS_DIR)
    p.add_argument("--prereg-lock-path", type=Path, default=DEFAULT_R17_PREREG_LOCK_PATH)
    p.add_argument("--models-dir", type=Path, default=DEFAULT_MODELS_DIR)
    p.add_argument("--results-dir", type=Path, default=DEFAULT_R17_RESULTS_DIR)
    p.add_argument("--gse295214-raw-dir", type=Path, default=DEFAULT_GSE295214_RAW_DIR)
    p.add_argument("--pfizer-raw-path", type=Path, default=DEFAULT_PFIZER_RAW_PATH)
    p.add_argument("--raw-nadig-dir", type=Path, default=None)
    p.add_argument("--device", choices=["mps", "cpu"], default="mps")
    p.add_argument("--n-boot", type=int, default=DEFAULT_R17_N_BOOT)
    p.add_argument("--n-boot-bonf56", type=int, default=DEFAULT_R17_N_BOOT_BONF56)
    p.add_argument("--n-perm", type=int, default=DEFAULT_R17_N_PERM)
    p.add_argument("--n-splits-diag", type=int, default=DEFAULT_R17_N_SPLITS_DIAG)
    p.add_argument("--seed", type=int, default=DEFAULT_R17_SEED)
    p.add_argument("--ctrl-seed", type=int, default=DEFAULT_R17_SEED)
    p.add_argument("--ctrl-resample-seed", type=int, default=DEFAULT_R17_RESAMPLE_SEED)
    p.add_argument("--skip-p2-to-p1-g4", action="store_true")
    p.add_argument("--smoke-out-dir", type=Path, default=None)
    p.add_argument("--recompute-from-arrays", type=Path, default=None)
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)

    if args.recompute_from_arrays is not None:
        res = recompute_cohort_from_exported_arrays_r17(
            arrays_dir=args.recompute_from_arrays,
            n_boot=args.n_boot,
            n_boot_bonf56=args.n_boot_bonf56,
            n_perm=args.n_perm,
            seed=args.seed,
            n_splits_diag=args.n_splits_diag,
        )
        print(json.dumps({"recomputed_lines": sorted(res["lines"].keys()), "K_conf": res["K_conf"]}, indent=2))
        return 0

    if args.synthetic_smoke:
        if args.smoke_out_dir is not None:
            res = run_synthetic_smoke_r17(
                args.smoke_out_dir,
                n_boot=min(args.n_boot, 800),
                n_boot_bonf56=min(args.n_boot_bonf56, 1000),
                n_perm=min(args.n_perm, 1000),
                seed=args.seed,
            )
        else:
            with tempfile.TemporaryDirectory(prefix="r17_synth_smoke_") as tmpd:
                res = run_synthetic_smoke_r17(
                    Path(tmpd),
                    n_boot=min(args.n_boot, 800),
                    n_boot_bonf56=min(args.n_boot_bonf56, 1000),
                    n_perm=min(args.n_perm, 1000),
                    seed=args.seed,
                )
        se_coh = res["separate_model_holm_fwer05"]["cohort_adjudications"]["STATE_SE"]
        hvg_coh = res["separate_model_holm_fwer05"]["cohort_adjudications"]["STATE_HVG"]
        compact = {
            "mode": res["mode"],
            "protocol": res["protocol"],
            "real_cells_touched": res["real_cells_touched"],
            "device": res["device"],
            "elapsed_sec": res["elapsed_sec"],
            "K_conf": res["K_conf"],
            "STATE_SE_k_over_K": se_coh["k_over_K"],
            "STATE_SE_general_claim": se_coh["general_cross_context_claim_permitted"],
            "STATE_SE_conf_verdicts": {
                ln: v["H_J1_holm"]["verdict"] for ln, v in se_coh["per_confirmatory_line"].items()
            },
            "STATE_HVG_k_over_K": hvg_coh["k_over_K"],
            "STATE_HVG_general_claim": hvg_coh["general_cross_context_claim_permitted"],
            "STATE_HVG_conf_verdicts": {
                ln: v["H_J1_holm"]["verdict"] for ln, v in hvg_coh["per_confirmatory_line"].items()
            },
            "stage4_from_disk_recompute_max_abs_diff": res["stage4_from_disk_recompute_max_abs_diff"],
            "stage4_from_disk_recompute_passed": res["stage4_from_disk_recompute_passed"],
            "n_exported_lines": len(res["lines"]),
            "n_exported_embeddings": len(res["embeddings_manifest_sha256"]),
        }
        print(json.dumps(compact, indent=2))
        return 0

    # STRICT EXECUTION GATE BEFORE ANY REAL-CELL ACCESS
    assert_all_execution_gates_r17(status_dir=args.status_dir, prereg_lock_path=args.prereg_lock_path)

    if args.lines:
        target_lines = list(args.lines)
    elif args.cohort == "confirmatory":
        target_lines = list(CONFIRMATORY_LINES_R17)
    elif args.cohort == "within_family":
        target_lines = list(WITHIN_FAMILY_LINES_R17)
    elif args.cohort == "descriptive":
        target_lines = list(SECONDARY_FRESH_LINES_R17) + list(MICRO_N_LINES_R17)
    elif args.cohort == "subset_sensitivity":
        target_lines = list(SUBSET_SENSITIVITY_LINES_R17)
    else:
        target_lines = list(ALL_EVALUATED_LINES_R17)

    run_real_evaluation_r17(
        lines_to_run=target_lines,
        status_dir=args.status_dir,
        prereg_lock_path=args.prereg_lock_path,
        models_dir=args.models_dir,
        results_dir=args.results_dir,
        device=args.device,
        n_boot=args.n_boot,
        n_boot_bonf56=args.n_boot_bonf56,
        n_perm=args.n_perm,
        seed=args.seed,
        ctrl_seed=args.ctrl_seed,
        ctrl_resample_seed=args.ctrl_resample_seed,
        n_splits_diag=args.n_splits_diag,
        gse295214_raw_dir=args.gse295214_raw_dir,
        pfizer_raw_path=args.pfizer_raw_path,
        raw_nadig_dir=args.raw_nadig_dir,
        include_p2_to_p1_g4=not args.skip_p2_to_p1_g4,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
