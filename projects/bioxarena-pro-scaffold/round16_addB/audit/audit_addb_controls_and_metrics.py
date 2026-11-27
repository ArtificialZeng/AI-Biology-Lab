#!/usr/bin/env python3
"""
Clean-Room Addendum B Control Barcode Extraction, Gate Guard, `H-I1_Occam`, `H-I2_Occam`,
Mandatory `P4` (`pairwise_common_axis`), Paired `STATE_SE` vs `Occam` Contrasts (`H-I3a`/`H-I3b`),
and §4.3 Advance Matrix Adjudicator (`round16_addB/audit/audit_addb_controls_and_metrics.py`, `A-T3`).

Strictly Python 3.8.19 compatible (`round16_addB/.venv`). Zero imports of `statsmodels`,
zero reads of `round16_addB/src/*.py`, `round16_addB/results/`, `round16/results/group3/`,
or `MAIN.txt` (beyond schema marker when released).
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Set, Tuple, Union

import h5py
import numpy as np
import scipy.sparse as sp

AUDIT_ADDB_DIR = Path(__file__).resolve().parent
ADDB_ROOT = AUDIT_ADDB_DIR.parent
REPO_ROOT = ADDB_ROOT.parent
R16_AUDIT_DIR = REPO_ROOT / "round16" / "audit"
R10_SRC_DIR = REPO_ROOT / "round10" / "src"

for _p in (str(AUDIT_ADDB_DIR), str(R16_AUDIT_DIR), str(R10_SRC_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from audit_addb_occam_adapter import (  # noqa: E402
    CONTROL_LABEL,
    DEFAULT_NADIG_HEPG2_H5AD,
    DEFAULT_NADIG_RAW_DIR,
    read_obs_or_var_col_legacy,
    sha256_file,
    sha256_str_lines,
    strip_ensembl_version,
)
from audit_r16_adjudicate import (  # noqa: E402
    DEFAULT_TOST_DELTA,
    MIN_BOOT_VALID_FRAC,
    _summarize_stratum_r16,
    adjudicate_claim_i1_line,
    adjudicate_claim_i2_line,
    classify_pair_tost_r16,
)
from audit_r16_metrics import (  # noqa: E402
    DEFAULT_R16_SEED,
    _extract_pert_gene_indices,
    _moments,
    _moments_to_corr_and_slope,
    _project_streams,
    compute_orth_suite_r16,
    compute_paired_contrast_r16,
    compute_paired_diff_vs_resample_r16,
    compute_paired_model_vs_all_mismatched_r16,
    compute_practical_metrics_r16,
    compute_winner_curse_weight_shift_r16,
    disattenuate,
    spearman_brown_ceiling,
    summarize_boot_r16,
)

EXPECTED_ADDB_V2_SHA256 = "537467db51b00b8e178ad5d7638e25b5d809f350c8feae229d86cdb96acb281e"
DEFAULT_LOCKED_SE_AUDIT_DIR = REPO_ROOT / "round16" / "audit" / "results"
DEFAULT_LOCKED_MAIN_SENS_DIR = REPO_ROOT / "round16" / "results" / "sensitivity"
DEFAULT_LOCKED_MAIN_R16_DIR = REPO_ROOT / "round16" / "results"
DEFAULT_GSE295214_RAW_DIR = Path("/Users/mac/.cache/r16_gse295214/raw")

ADDB_GROUP2_LINES: Tuple[str, ...] = ("a549_ko", "dld1", "cd4t_rest")
ADDB_GROUP3_LINES: Tuple[str, ...] = ("gse295214_h1",)
ADDB_SENSITIVITY_LINES: Tuple[str, ...] = ("feng_ipsc",)
ADDB_ALL_TARGET_LINES: Tuple[str, ...] = ADDB_GROUP2_LINES + ADDB_GROUP3_LINES + ADDB_SENSITIVITY_LINES
ADDB_ALL_CONTROL_POOLS: Tuple[str, ...] = (
    "a549_ko",
    "dld1",
    "cd4t_rest",
    "nadig_hepg2",
    "nadig_jurkat",
    "feng_ipsc",
    "gse295214_h1",
)

ADDB_SEEN_SOURCES: Tuple[str, ...] = ("nadig_hepg2", "nadig_jurkat")

ADDB_LINE_GROUP_META: Dict[str, Dict[str, str]] = {
    "a549_ko": {"group_key": "group2_crisprko", "role": "primary", "modality": "CRISPRko"},
    "dld1": {"group_key": "group2_crispri", "role": "primary", "modality": "CRISPRi"},
    "cd4t_rest": {"group_key": "group2_crispri", "role": "primary", "modality": "CRISPRi"},
    "gse295214_h1": {"group_key": "group3_crispri", "role": "prospective_group3", "modality": "CRISPRi"},
    "feng_ipsc": {"group_key": "sensitivity_only", "role": "sensitivity_only", "modality": "CRISPRi"},
}

ADDB_LINE_ROLES: Dict[str, str] = {
    "a549_ko": "group2_primary_crisprko",
    "dld1": "group2_primary_crispri",
    "cd4t_rest": "group2_primary_crispri",
    "gse295214_h1": "group3_prospective_holdout",
    "feng_ipsc": "sensitivity_only",
}

ADDB_UNSEEN_SOURCES: Dict[str, Tuple[str, ...]] = {
    "a549_ko": ("dld1", "cd4t_rest", "feng_ipsc"),
    "dld1": ("a549_ko", "cd4t_rest", "feng_ipsc"),
    "cd4t_rest": ("a549_ko", "dld1", "feng_ipsc"),
    "gse295214_h1": ("a549_ko", "cd4t_rest", "dld1", "feng_ipsc"),
    "feng_ipsc": ("a549_ko", "cd4t_rest", "dld1"),
}

R9_6_ALL_LINES: Tuple[str, ...] = (
    "nadig_hepg2",
    "nadig_jurkat",
    "nourreddine_kolf21j",
    "replogle_k562_essential",
    "replogle_rpe1",
    "vcc2025_h1",
)
GROUP3_PRIMARY_FCT_SOURCES: Tuple[str, ...] = (
    "nadig_hepg2",
    "nadig_jurkat",
    "nourreddine_kolf21j",
    "replogle_k562_essential",
    "replogle_rpe1",
)

# Pre-locked §3.3.4 P4 lines lacking some RHJ_6983 genes (`cd4t_rest` 413 absent, `a549_ko` 33 absent)
ADDB_P4_PARTIAL_AXIS_POOLS: Tuple[str, ...] = ("a549_ko", "cd4t_rest", "gse295214_h1")


class AddBBlindingGateError(RuntimeError):
    """Raised when Addendum B real-target freeze or GO authorization gates are absent."""


@dataclass
class ExtractedRawControlCells:
    pool_name: str
    raw_counts: np.ndarray        # (n_cells, n_raw_genes) float32
    raw_ensembl_ids: List[str]    # stripped version-free Ensembl IDs (GRCh38 RNA order)
    barcodes: List[str]           # exact ordered 512 barcodes
    row_sum_umis: np.ndarray      # (n_cells,) float64 full raw feature sums (`DEVIATION_R16_D8`)
    n_ctrl_total: Optional[int] = None


@dataclass
class CommittedR16ReferenceBundle:
    target_line: str
    role: str
    eval_perts: List[str]
    eval_genes: List[str]
    target_gene_idx: np.ndarray        # (P,) int64
    delta_obs: np.ndarray              # (P, G) float64
    delta_half_a: np.ndarray           # (P, G) float64
    delta_half_b: np.ndarray           # (P, G) float64
    half_valid: np.ndarray             # (P,) bool
    delta_transfer: np.ndarray         # (P, G) float64
    source_grand_mean: np.ndarray      # (G,) float64
    delta_fct: np.ndarray              # (P, G) float64
    seen_sources: List[str]
    unseen_sources: List[str]
    matched_barcodes: List[str]
    resample_barcodes: List[str]
    mismatched_barcodes: Dict[str, List[str]]
    se_delta_matched: np.ndarray       # (P, G) float64
    se_delta_resample: np.ndarray      # (P, G) float64
    se_delta_mismatched: Dict[str, np.ndarray]
    mu_target: Optional[np.ndarray] = None              # (G,) float64
    mu_sources: Optional[np.ndarray] = None             # (S_src, G) float64
    e1_constmed: Optional[np.ndarray] = None            # (P, G) float64
    n_ctrl_total: int = 512
    ctrl_sampling_fraction: float = 1.0
    r16_reference_sha256: Optional[Dict[str, str]] = None
    extra_fct_full_r9_6: Optional[Dict[str, Any]] = None


def assert_freeze_gates_addb(addb_root: Optional[Path] = None) -> Dict[str, str]:
    """
    Hard runtime freeze & authorization gate for real Addendum B inference:
      Refuses unless BOTH `R16_ADDB_CODE_FREEZE_AT3.txt` and `R16_ADDB_GO_AT3.txt`
      exist and are non-empty under `addb_root` (default `round16_addB/`), AND
      `PREREG_I_ADDENDUM_B_DRAFT_v2.md` matches locked SHA-256.
    """
    root = Path(addb_root) if addb_root is not None else ADDB_ROOT
    freeze_fp = root / "R16_ADDB_CODE_FREEZE_AT3.txt"
    go_fp = root / "R16_ADDB_GO_AT3.txt"
    lock_fp = root / "PREREG_I_ADDENDUM_B_LOCK.md"
    draft_fp = root / "PREREG_I_ADDENDUM_B_DRAFT_v2.md"
    approval_fp = root / "R16_ADDB_HUMAN_APPROVAL.txt"

    for label, fp in (
        ("R16_ADDB_CODE_FREEZE_AT3.txt", freeze_fp),
        ("R16_ADDB_GO_AT3.txt", go_fp),
        ("PREREG_I_ADDENDUM_B_LOCK.md", lock_fp),
    ):
        if not fp.is_file():
            raise AddBBlindingGateError(f"Missing required Addendum B gate file ({label}): {fp}")
        if not fp.read_text(encoding="utf-8").strip():
            raise AddBBlindingGateError(f"Empty required Addendum B gate file ({label}): {fp}")

    if draft_fp.is_file():
        got_sha = sha256_file(draft_fp)
        if got_sha != EXPECTED_ADDB_V2_SHA256:
            raise AddBBlindingGateError(
                f"Addendum B draft v2 SHA-256 mismatch: {got_sha} != {EXPECTED_ADDB_V2_SHA256}"
            )

    return {
        "code_freeze_path": str(freeze_fp),
        "go_path": str(go_fp),
        "lock_path": str(lock_fp),
        "human_approval_path": str(approval_fp) if approval_fp.is_file() else "",
        "r16_addb_code_freeze_at3": freeze_fp.read_text(encoding="utf-8").strip().splitlines()[0],
        "r16_addb_go_at3": go_fp.read_text(encoding="utf-8").strip().splitlines()[0],
        "addb_v2_sha256": EXPECTED_ADDB_V2_SHA256,
    }


def _resolve_r16_source_line_dir(line_name: str, repo_root: Path = REPO_ROOT) -> Path:
    r9_set = {
        "nadig_hepg2",
        "nadig_jurkat",
        "replogle_k562_essential",
        "replogle_rpe1",
        "vcc2025_h1",
        "nourreddine_kolf21j",
    }
    if line_name in r9_set:
        return Path(repo_root) / "round9" / "p2" / "shared" / "lines" / line_name
    if line_name in ("dld1", "cd4t_rest"):
        return Path(repo_root) / "round10" / "shared" / "lines" / line_name
    if line_name in ("feng_ipsc", "cd4t_stim48"):
        return Path(repo_root) / "round11" / "shared" / "lines" / line_name
    if line_name in ("a549_ko", "hek293_xucao"):
        return Path(repo_root) / "round12" / "shared" / "lines" / line_name
    if line_name == "gse295214_h1":
        return Path(repo_root) / "round16" / "shared" / "lines" / "gse295214_h1"
    raise ValueError(f"Unknown standardized line: {line_name!r}")


def _build_gse295214_observed_extras(
    eval_genes: Sequence[str],
    eval_perts: Sequence[str],
    target_gene_idx: np.ndarray,
    mu_target: np.ndarray,
    hgnc_path: Path = Path("/Users/mac/.cache/r9_p2_raw/hgnc_complete_set.txt"),
    repo_root: Path = REPO_ROOT,
) -> Dict[str, Any]:
    """
    Reconstruct `mu_sources` (`(5, G)`), `e1_constmed_delta_pred` (`(P, G)`), `n_ctrl_total`,
    and `preds_SENS_FCT_full_r9_6_incl_vcc2025_h1` (`(P, G)`) for `gse295214_h1` cleanly from
    standardized R9/R10/R16 inputs without ever reading `round16/results/group3/`.
    """
    from data_loader_r10 import (
        HGNCBridge,
        LineDatasetR10,
        load_r9_common_genes,
        load_standardized_line_r10,
    )
    from fct_pipeline import (
        FctSourceState,
        build_source_input_r10,
        predict_group_c_suite,
    )

    hb = HGNCBridge.from_tsv(Path(hgnc_path))
    r9_common = load_r9_common_genes(hgnc_bridge=hb)
    eval_g_list = list(eval_genes)
    eval_g_to_idx = {g: i for i, g in enumerate(eval_g_list)}
    tg_arr = np.asarray(target_gene_idx, dtype=np.int64)
    mu_t = np.maximum(np.asarray(mu_target, dtype=np.float64), 0.0)

    src_ds_cache: Dict[str, Any] = {}
    for sn in R9_6_ALL_LINES:
        ld = _resolve_r16_source_line_dir(sn, repo_root=repo_root)
        sd = load_standardized_line_r10(
            ld,
            representation="log1p",
            common_genes=r9_common,
            hgnc_bridge=hb,
            min_gene_overlap=1000,
        )
        pos = {g: i for i, g in enumerate(sd.genes)}
        cols = np.asarray([pos[g] for g in eval_g_list], dtype=np.int64)
        sub_sd = LineDatasetR10(
            name=sd.name,
            genes=eval_g_list,
            target_names=sd.target_names,
            target_indices=np.asarray([eval_g_to_idx.get(t, -1) for t in sd.target_names], dtype=np.int64),
            deltas=sd.deltas[:, cols],
            ctrl_cells=sd.ctrl_cells[:, cols],
            ctrl_mean_log1p=sd.ctrl_mean_log1p[cols] if sd.ctrl_mean_log1p is not None else None,
            ctrl_mean_libnorm=sd.ctrl_mean_libnorm[cols] if sd.ctrl_mean_libnorm is not None else None,
            ctrl_mean_libnorm_full=sd.ctrl_mean_libnorm_full[cols] if sd.ctrl_mean_libnorm_full is not None else None,
            repr_mask=sd.repr_mask,
            half_a_deltas=sd.half_a_deltas[:, cols] if sd.half_a_deltas is not None else None,
            half_b_deltas=sd.half_b_deltas[:, cols] if sd.half_b_deltas is not None else None,
            half_valid=sd.half_valid,
            n_cells=sd.n_cells,
            knockdown_depth_log1p=sd.knockdown_depth_log1p,
            knockdown_fraction=sd.knockdown_fraction,
            target_ensembl_ids=sd.target_ensembl_ids,
            representation=sd.representation,
            metadata=sd.metadata,
        )
        src_ds_cache[sn] = (
            build_source_input_r10(sub_sd, only_reproducible=True, compute_cov=False),
            np.maximum(np.asarray(sd.ctrl_mean_libnorm[cols], dtype=np.float64), 0.0),
        )

    pri_inputs = [src_ds_cache[sn][0] for sn in GROUP3_PRIMARY_FCT_SOURCES]
    mu_sources = np.stack([src_ds_cache[sn][1] for sn in GROUP3_PRIMARY_FCT_SOURCES], axis=0)
    fst_pri = FctSourceState.fit(pri_inputs, mask_self=True)
    c_med = float(np.median(mu_t))
    preds_e1, _ = predict_group_c_suite(
        state=fst_pri,
        target_ctrl_mean_libnorm=np.full_like(mu_t, c_med),
        query_target_names=list(eval_perts),
        query_target_gene_indices=tg_arr,
        include_r9_anchors=False,
    )

    full_inputs = [src_ds_cache[sn][0] for sn in R9_6_ALL_LINES]
    fst_full = FctSourceState.fit(full_inputs, mask_self=True)
    preds_full, _ = predict_group_c_suite(
        state=fst_full,
        target_ctrl_mean_libnorm=mu_t,
        query_target_names=list(eval_perts),
        query_target_gene_indices=tg_arr,
        include_r9_anchors=False,
    )

    meta_fp = Path(repo_root) / "round16" / "shared" / "lines" / "gse295214_h1" / "line_metadata.json"
    n_ctrl_total = 512
    if meta_fp.is_file():
        meta_j = json.loads(meta_fp.read_text(encoding="utf-8"))
        n_ctrl_total = int(meta_j.get("n_nt_total") or meta_j.get("adapter_audit", {}).get("n_kept_ntc", 512))

    return {
        "mu_sources": np.asarray(mu_sources, dtype=np.float64),
        "e1_constmed": np.asarray(preds_e1["FCT"], dtype=np.float64),
        "n_ctrl_total": int(n_ctrl_total),
        "extra_fct_full_r9_6": {
            "source_lines_fct": list(R9_6_ALL_LINES),
            "delta_pred": np.asarray(preds_full["FCT"], dtype=np.float64),
            "delta_transfer": np.asarray(preds_full["Transfer_Mean"], dtype=np.float64),
            "source_grand_mean": np.asarray(fst_full.grand_mean, dtype=np.float64),
        },
    }


def load_locked_control_barcodes_addb(
    target_lines: Sequence[str],
    audit_results_dir: Path = DEFAULT_LOCKED_SE_AUDIT_DIR,
    main_r16_dir: Path = DEFAULT_LOCKED_MAIN_R16_DIR,
    main_sens_dir: Path = DEFAULT_LOCKED_MAIN_SENS_DIR,
    repo_root: Path = REPO_ROOT,
) -> Dict[str, Any]:
    """
    Load exact ordered `512`-cell primary (`seed=20261601`) and within-line resample
    (`seed=20261602`, `0` overlap) control barcode lists verbatim from R16 locked JSONs
    (`PREREG_I_ADDENDUM_B_DRAFT_v2.md §3.2`) + provenance SHA-256s.
    Never reads `round16/results/group3/`.
    """
    needed_primary: Set[str] = set()
    needed_resample: Set[str] = set()
    for ln in target_lines:
        needed_primary.add(ln)
        needed_resample.add(ln)
        for src in ADDB_SEEN_SOURCES + ADDB_UNSEEN_SOURCES[ln]:
            needed_primary.add(src)

    audit_cb_path = Path(audit_results_dir) / "control_barcodes_r16.json"
    g3_cb_path = Path(audit_results_dir) / "group3" / "control_barcodes_r16.json"
    main_sens_cb_path = Path(main_sens_dir) / "control_barcodes_r16.json"
    main_cb_path = Path(main_r16_dir) / "control_barcodes_r16.json"

    audit_cb = json.loads(audit_cb_path.read_text(encoding="utf-8"))
    g3_cb = json.loads(g3_cb_path.read_text(encoding="utf-8")) if g3_cb_path.is_file() else {}
    main_sens_cb = json.loads(main_sens_cb_path.read_text(encoding="utf-8")) if main_sens_cb_path.is_file() else {}
    main_cb = json.loads(main_cb_path.read_text(encoding="utf-8")) if main_cb_path.is_file() else {}

    sources_sha256: Dict[str, str] = {}
    for fp in (audit_cb_path, g3_cb_path, main_cb_path, main_sens_cb_path):
        if fp.is_file():
            try:
                rel = str(fp.relative_to(repo_root))
            except ValueError:
                rel = str(fp)
            sources_sha256[rel] = sha256_file(fp)

    if g3_cb:
        for k, gval in g3_cb.get("primary_control_pools", {}).items():
            if k in audit_cb.get("primary_control_pools", {}):
                if [str(x) for x in gval["barcodes"]] != [str(x) for x in audit_cb["primary_control_pools"][k]["barcodes"]]:
                    raise ValueError(f"Primary barcode mismatch between audit & g3 audit for {k!r}")

    pools_out: Dict[str, Dict[str, Any]] = {}

    for pname in sorted(needed_primary):
        entry_p = None
        if pname == "gse295214_h1":
            entry_p = g3_cb.get("primary_control_pools", {}).get(pname)
        else:
            entry_p = audit_cb.get("primary_control_pools", {}).get(pname)
            if entry_p is None and pname in main_cb:
                entry_p = main_cb[pname]
        if entry_p is None or not entry_p.get("barcodes"):
            raise FileNotFoundError(f"Primary locked barcodes missing for pool {pname!r}")
        p_bcs = [str(b) for b in entry_p["barcodes"]]
        if main_cb.get(pname, {}).get("barcodes") and pname != "gse295214_h1":
            main_bcs = [str(b) for b in main_cb[pname]["barcodes"]]
            if p_bcs != main_bcs:
                raise ValueError(f"Primary barcode mismatch between audit & main R16 for {pname!r}")

        n_ctrl_total = entry_p.get("n_ctrl_total", entry_p.get("n_valid_pool"))
        if n_ctrl_total is None and pname in main_cb:
            n_ctrl_total = main_cb[pname].get("n_ctrl_total")
        if n_ctrl_total is None and pname == "feng_ipsc" and pname in main_sens_cb:
            n_ctrl_total = main_sens_cb[pname].get("n_ctrl_total")
        if n_ctrl_total is None and pname == "gse295214_h1":
            meta_fp = Path(repo_root) / "round16" / "shared" / "lines" / "gse295214_h1" / "line_metadata.json"
            if meta_fp.is_file():
                meta_j = json.loads(meta_fp.read_text(encoding="utf-8"))
                n_ctrl_total = meta_j.get("n_nt_total") or meta_j.get("adapter_audit", {}).get("n_kept_ntc")
        if n_ctrl_total is None:
            n_ctrl_total = len(p_bcs)

        res_bcs: Optional[List[str]] = None
        if pname in needed_resample:
            if pname == "gse295214_h1":
                rentry = g3_cb.get("resample_control_pools", {}).get(pname)
            elif pname == "feng_ipsc":
                rentry = main_sens_cb.get("feng_ipsc", {}).get("within_line_resample")
                if rentry is None:
                    obs_p = Path(main_sens_dir) / "arrays" / "feng_ipsc" / "observed.npz"
                    if obs_p.is_file():
                        with np.load(obs_p, allow_pickle=False) as oz:
                            if "resample_ctrl_barcodes" in oz.files:
                                rentry = {"barcodes": [str(b) for b in oz["resample_ctrl_barcodes"].tolist()]}
            else:
                rentry = audit_cb.get("resample_control_pools", {}).get(pname)
                if rentry is None:
                    rentry = main_cb.get(pname, {}).get("within_line_resample")
            if rentry is None or not rentry.get("barcodes"):
                raise FileNotFoundError(f"Resample locked barcodes missing for target {pname!r}")
            res_bcs = [str(b) for b in rentry["barcodes"]]
            overlap = set(p_bcs) & set(res_bcs)
            if overlap:
                raise ValueError(
                    f"Disjoint resample violated on {pname!r}: {len(overlap)} overlapping barcodes"
                )

        pools_out[pname] = {
            "primary_barcodes": p_bcs,
            "resample_barcodes": res_bcs,
            "n_ctrl_total": int(n_ctrl_total),
        }

    return {
        "pools": pools_out,
        "provenance": {
            "primary_seed": 20261601,
            "resample_seed": 20261602,
            "r16_sources_sha256": sources_sha256,
        },
    }


def format_control_barcodes_addb_json(locked_data: Mapping[str, Any]) -> Dict[str, Any]:
    """Serialize loaded primary + resample pools into exact `control_barcodes_addb.json` schema."""
    pools_in = locked_data["pools"] if "pools" in locked_data else locked_data
    prov = locked_data.get(
        "provenance",
        {"primary_seed": 20261601, "resample_seed": 20261602, "r16_sources_sha256": {}},
    )
    formatted_pools: Dict[str, Any] = {}
    for pname in sorted(pools_in.keys()):
        info = pools_in[pname]
        p_bcs = [str(b) for b in info["primary_barcodes"]]
        r_bcs = info.get("resample_barcodes")
        if r_bcs is not None:
            r_list = [str(b) for b in r_bcs]
            w_block: Dict[str, Any] = {
                "eligible": True,
                "seed": 20261602,
                "n_resample": len(r_list),
                "overlap_with_primary": len(set(p_bcs) & set(r_list)),
                "regime": "full_512_disjoint" if len(r_list) == 512 else f"full_{len(r_list)}_disjoint",
                "barcodes": r_list,
                "barcodes_sha256": sha256_str_lines(r_list),
            }
        else:
            w_block = {"eligible": False}
        formatted_pools[pname] = {
            "primary_seed": 20261601,
            "n_ctrl_total": int(info.get("n_ctrl_total", len(p_bcs))),
            "n_ctrl_sampled": len(p_bcs),
            "barcodes": p_bcs,
            "barcodes_sha256": sha256_str_lines(p_bcs),
            "within_line_resample": w_block,
        }
    return {
        "provenance": prov,
        "pools": formatted_pools,
    }


# ==============================================================================
# Targeted Barcode-Exact Raw Single-Cell Extractors (Py3.8, Zero statsmodels)
# ==============================================================================


def _read_obs_index_barcodes_h5ad(obs_grp: h5py.Group) -> List[str]:
    if "cell_barcode" in obs_grp:
        bcs = read_obs_or_var_col_legacy(obs_grp, "cell_barcode")
        if bcs and any(b for b in bcs):
            return bcs
    idx_attr = obs_grp.attrs.get("_index", "_index")
    if isinstance(idx_attr, (bytes, np.bytes_)):
        idx_attr = idx_attr.decode("utf-8", errors="replace")
    col = str(idx_attr) if str(idx_attr) in obs_grp else ("_index" if "_index" in obs_grp else "index")
    return read_obs_or_var_col_legacy(obs_grp, col)


def _extract_nadig_by_barcodes(
    h5ad_path: Path,
    pool_name: str,
    barcodes: Sequence[str],
    expected_sha256: Optional[str] = None,
) -> ExtractedRawControlCells:
    if not h5ad_path.is_file():
        raise FileNotFoundError(f"Missing Nadig raw file: {h5ad_path}")
    if expected_sha256:
        got = sha256_file(h5ad_path)
        if got != expected_sha256:
            raise ValueError(f"SHA-256 mismatch on {h5ad_path}: {got} != {expected_sha256}")

    with h5py.File(str(h5ad_path), "r") as f:
        obs_bcs = _read_obs_index_barcodes_h5ad(f["obs"])
        var_ens = [strip_ensembl_version(g) for g in read_obs_or_var_col_legacy(f["var"], "gene_id")]
        bc_to_row = {b: i for i, b in enumerate(obs_bcs)}
        req_rows: List[int] = []
        missing: List[str] = []
        for bc in barcodes:
            ri = bc_to_row.get(str(bc))
            if ri is None:
                missing.append(str(bc))
            else:
                req_rows.append(ri)
        if missing:
            raise KeyError(f"{pool_name}: {len(missing)} locked barcodes not found in {h5ad_path} (e.g. {missing[:3]})")

        req_arr = np.asarray(req_rows, dtype=np.int64)
        order = np.argsort(req_arr, kind="mergesort")
        inv_order = np.empty_like(order)
        inv_order[order] = np.arange(order.size)
        sorted_rows = req_arr[order]

        x_node = f["X"]
        n_vars = len(var_ens)
        sorted_mat = np.empty((sorted_rows.size, n_vars), dtype=np.float32)
        if isinstance(x_node, h5py.Dataset):
            sorted_mat[:] = np.asarray(x_node[sorted_rows, :], dtype=np.float32)
        else:
            indptr = np.asarray(x_node["indptr"][:], dtype=np.int64)
            data = x_node["data"]
            indices = x_node["indices"]
            sorted_mat.fill(0.0)
            for out_i, r in enumerate(sorted_rows.tolist()):
                p0, p1 = int(indptr[r]), int(indptr[r + 1])
                if p1 > p0:
                    sorted_mat[out_i, indices[p0:p1]] = np.asarray(data[p0:p1], dtype=np.float32)
        raw_mat = sorted_mat[inv_order]

    return ExtractedRawControlCells(
        pool_name=pool_name,
        raw_counts=raw_mat,
        raw_ensembl_ids=var_ens,
        barcodes=[str(b) for b in barcodes],
        row_sum_umis=raw_mat.sum(axis=1, dtype=np.float64),
        n_ctrl_total=len(obs_bcs),
    )


def _extract_10x_h5_subset(h5_path: Path, need_bcs: Sequence[str]) -> Tuple[np.ndarray, List[str], List[str]]:
    """Read ONLY requested `need_bcs` columns and `Gene Expression` features from one 10x `.h5`."""
    with h5py.File(str(h5_path), "r") as f:
        m = f["matrix"]
        bcs_raw = [
            x.decode("utf-8", errors="replace") if isinstance(x, (bytes, np.bytes_)) else str(x)
            for x in m["barcodes"][:]
        ]
        f_ids = [
            strip_ensembl_version(x.decode("utf-8", errors="replace") if isinstance(x, (bytes, np.bytes_)) else str(x))
            for x in m["features"]["id"][:]
        ]
        f_types = [
            x.decode("utf-8", errors="replace") if isinstance(x, (bytes, np.bytes_)) else str(x)
            for x in m["features"]["feature_type"][:]
        ]
        rna_idx = np.asarray([i for i, ft in enumerate(f_types) if ft == "Gene Expression"], dtype=np.int64)
        if rna_idx.size == 0:
            rna_idx = np.arange(len(f_ids), dtype=np.int64)
        rna_ens = [f_ids[i] for i in rna_idx.tolist()]
        rna_pos_map = np.full(len(f_ids), -1, dtype=np.int64)
        rna_pos_map[rna_idx] = np.arange(rna_idx.size, dtype=np.int64)

        bc_map = {b: i for i, b in enumerate(bcs_raw)}
        data = m["data"]
        indices = m["indices"]
        indptr = np.asarray(m["indptr"][:], dtype=np.int64)

        mat = np.zeros((len(need_bcs), rna_idx.size), dtype=np.float32)
        found_bcs: List[str] = []
        for row_i, bc in enumerate(need_bcs):
            ci = bc_map.get(bc)
            if ci is None:
                raise KeyError(f"10x barcode {bc!r} missing from {h5_path}")
            p0, p1 = int(indptr[ci]), int(indptr[ci + 1])
            if p1 > p0:
                feat_idx = np.asarray(indices[p0:p1], dtype=np.int64)
                vals = np.asarray(data[p0:p1], dtype=np.float32)
                mapped = rna_pos_map[feat_idx]
                keep = mapped >= 0
                if np.any(keep):
                    np.add.at(mat[row_i], mapped[keep], vals[keep])
            found_bcs.append(bc)
    return mat, rna_ens, found_bcs


def _resolve_10x_channel_file(raw_dir: Path, channel_key: str, pool_name: str) -> Path:
    cands = sorted(p for p in Path(raw_dir).glob("*.h5") if ".hashing." not in p.name and ".crispr." not in p.name)
    if pool_name == "cd4t_rest":
        parts = channel_key.split("_")
        donor = next((t for t in parts if len(t) == 2 and t.startswith("D") and t[1:].isdigit()), "")
        lane_idx = channel_key.find("CD4")
        lane = channel_key[lane_idx:] if lane_idx >= 0 else channel_key
        hits = [p for p in cands if f"_{lane}_" in p.name and (not donor or f"_{donor}_" in p.name)]
    else:
        hits = [p for p in cands if p.name.startswith(f"{channel_key}_")]
    if len(hits) != 1:
        raise FileNotFoundError(
            f"Expected unique 10x .h5 file for {pool_name}:{channel_key} in {raw_dir}, found {[p.name for p in hits]}"
        )
    return hits[0]


def _extract_10x_by_barcodes(
    raw_dir: Path,
    pool_name: str,
    barcodes: Sequence[str],
) -> ExtractedRawControlCells:
    """Group prefixed `{prefix}:{bc}` barcodes by channel `.h5` file and slice exact ordered rows."""
    by_ch: Dict[str, List[Tuple[int, str]]] = {}
    for pos, full_bc in enumerate(barcodes):
        if ":" not in str(full_bc):
            raise ValueError(f"Expected prefixed barcode 'prefix:bc' for {pool_name}, got {full_bc!r}")
        ch, raw_bc = str(full_bc).split(":", 1)
        by_ch.setdefault(ch, []).append((pos, raw_bc))

    ref_ens: Optional[List[str]] = None
    out_mat: Optional[np.ndarray] = None
    for ch, items in sorted(by_ch.items()):
        hp = _resolve_10x_channel_file(raw_dir, ch, pool_name)
        ch_mat, ch_ens, _ = _extract_10x_h5_subset(hp, [bc for _, bc in items])
        if ref_ens is None:
            ref_ens = ch_ens
            out_mat = np.empty((len(barcodes), len(ref_ens)), dtype=np.float32)
        elif ch_ens != ref_ens:
            raise ValueError(f"Feature ID mismatch inside {hp} relative to first channel of {pool_name}")
        assert out_mat is not None
        for local_i, (pos, _) in enumerate(items):
            out_mat[pos] = ch_mat[local_i]

    assert out_mat is not None and ref_ens is not None
    return ExtractedRawControlCells(
        pool_name=pool_name,
        raw_counts=out_mat,
        raw_ensembl_ids=ref_ens,
        barcodes=[str(b) for b in barcodes],
        row_sum_umis=out_mat.sum(axis=1, dtype=np.float64),
    )


def _extract_a549_multi_pools_by_barcodes(
    raw_dir: Path,
    barcode_lists: Mapping[str, Sequence[str]],
) -> Dict[str, ExtractedRawControlCells]:
    """
    Extract exact requested `cell_id` sets (e.g. primary + resample `a549_ko`) in ONE
    single streaming pass over `GSE345058_raw_counts.mtx.gz` keeping `genome != 'CAS9'` GRCh38 rows.
    """
    cell_ann = Path(raw_dir) / "GSE345058_cell_annotations.csv.gz"
    gene_ann = Path(raw_dir) / "GSE345058_gene_annotations.csv.gz"
    mtx_gz = Path(raw_dir) / "GSE345058_raw_counts.mtx.gz"

    wanted: Set[str] = set()
    for bcs in barcode_lists.values():
        wanted.update(str(b) for b in bcs)

    cell_to_row1: Dict[str, int] = {}
    with gzip.open(cell_ann, "rt", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row_idx_1, row in enumerate(reader, start=1):
            cid = str(row["cell_id"]).strip()
            if cid in wanted:
                cell_to_row1[cid] = row_idx_1

    missing_c = wanted - set(cell_to_row1.keys())
    if missing_c:
        raise KeyError(f"a549_ko: {len(missing_c)} barcodes absent from cell_annotations (e.g. {sorted(missing_c)[:3]})")

    col1_to_gidx: Dict[int, int] = {}
    grch38_ens: List[str] = []
    with gzip.open(gene_ann, "rt", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for col1, row in enumerate(reader, start=1):
            if str(row.get("genome", "")).strip().upper() == "CAS9":
                continue
            eid = strip_ensembl_version(str(row["ensembl_id"]))
            col1_to_gidx[col1] = len(grch38_ens)
            grch38_ens.append(eid)

    unique_rows1 = sorted(set(cell_to_row1.values()))
    row1_to_buf = {r1: idx for idx, r1 in enumerate(unique_rows1)}
    buf = np.zeros((len(unique_rows1), len(grch38_ens)), dtype=np.float32)

    seen_dims = False
    row1_set = set(row1_to_buf.keys())
    with gzip.open(mtx_gz, "rt", encoding="ascii", errors="replace") as fh:
        for line in fh:
            if not line or line[0] == "%":
                continue
            if not seen_dims:
                seen_dims = True
                continue
            parts = line.split()
            r1 = int(parts[0])
            if r1 not in row1_set:
                continue
            c1 = int(parts[1])
            gi = col1_to_gidx.get(c1)
            if gi is not None:
                buf[row1_to_buf[r1], gi] += float(parts[2])

    out: Dict[str, ExtractedRawControlCells] = {}
    for pname, bcs in barcode_lists.items():
        idx_arr = np.asarray([row1_to_buf[cell_to_row1[str(b)]] for b in bcs], dtype=np.int64)
        mat = np.ascontiguousarray(buf[idx_arr], dtype=np.float32)
        out[pname] = ExtractedRawControlCells(
            pool_name=pname,
            raw_counts=mat,
            raw_ensembl_ids=list(grch38_ens),
            barcodes=[str(b) for b in bcs],
            row_sum_umis=mat.sum(axis=1, dtype=np.float64),
        )
    return out


def _extract_feng_multi_pools_by_barcodes(
    raw_dir: Path,
    barcode_lists: Mapping[str, Sequence[str]],
) -> Dict[str, ExtractedRawControlCells]:
    """
    Extract exact requested barcodes from `TargetedScreen_RNA-UMI-Counts.csv.gz`
    in a single streaming pass (`PREREG_I_ADDENDUM_B §3.3.1`).
    """
    csv_gz = Path(raw_dir) / "TargetedScreen_RNA-UMI-Counts.csv.gz"
    if not csv_gz.is_file():
        raise FileNotFoundError(f"Missing feng_ipsc RNA matrix: {csv_gz}")

    wanted: Set[str] = set()
    for bcs in barcode_lists.values():
        wanted.update(str(b) for b in bcs)

    rows_ens: List[str] = []
    cols_data: List[List[float]] = []

    with gzip.open(csv_gz, "rt", encoding="utf-8", errors="replace") as fh:
        header_line = fh.readline().rstrip("\r\n")
        headers = header_line.split(",")
        bc_to_col: Dict[str, int] = {}
        for ci, tok in enumerate(headers[1:], start=1):
            t = tok.strip()
            if t in wanted:
                bc_to_col[t] = ci
        missing = wanted - set(bc_to_col.keys())
        if missing:
            raise KeyError(f"feng_ipsc: {len(missing)} barcodes absent from header (e.g. {sorted(missing)[:3]})")

        unique_bcs = sorted(bc_to_col.keys())
        sel_cols = [bc_to_col[b] for b in unique_bcs]
        bc_to_bufcol = {b: j for j, b in enumerate(unique_bcs)}

        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\r\n").split(",")
            feat = parts[0].strip()
            eid = strip_ensembl_version(feat.split(":", 1)[0])
            rows_ens.append(eid)
            cols_data.append([float(parts[ci]) for ci in sel_cols])

    buf = np.asarray(cols_data, dtype=np.float32).T
    out: Dict[str, ExtractedRawControlCells] = {}
    for pname, bcs in barcode_lists.items():
        idx_arr = np.asarray([bc_to_bufcol[str(b)] for b in bcs], dtype=np.int64)
        mat = np.ascontiguousarray(buf[idx_arr], dtype=np.float32)
        out[pname] = ExtractedRawControlCells(
            pool_name=pname,
            raw_counts=mat,
            raw_ensembl_ids=list(rows_ens),
            barcodes=[str(b) for b in bcs],
            row_sum_umis=mat.sum(axis=1, dtype=np.float64),
        )
    return out


def extract_all_required_raw_pools_addb(
    locked_barcodes: Mapping[str, Any],
    nadig_raw_dir: Path = DEFAULT_NADIG_RAW_DIR,
    gse295214_raw_dir: Path = DEFAULT_GSE295214_RAW_DIR,
    r10_raw_dir: Path = REPO_ROOT / "round10" / "shared" / "raw",
    r11_raw_dir: Path = REPO_ROOT / "round11" / "shared" / "raw",
    r12_raw_dir: Path = REPO_ROOT / "round12" / "shared" / "raw",
) -> Tuple[Dict[str, ExtractedRawControlCells], Dict[str, ExtractedRawControlCells]]:
    """
    Extract exact raw count pools for all locked primary and resample barcode sets.
    Returns `(primary_pools, resample_pools)` keyed by pool name.
    """
    pools_map = locked_barcodes["pools"] if "pools" in locked_barcodes else locked_barcodes
    primary: Dict[str, ExtractedRawControlCells] = {}
    resample: Dict[str, ExtractedRawControlCells] = {}

    if "nadig_hepg2" in pools_map:
        primary["nadig_hepg2"] = _extract_nadig_by_barcodes(
            Path(nadig_raw_dir) / "GSE264667_hepg2_raw_singlecell_01.h5ad",
            "nadig_hepg2",
            pools_map["nadig_hepg2"]["primary_barcodes"],
        )
    if "nadig_jurkat" in pools_map:
        primary["nadig_jurkat"] = _extract_nadig_by_barcodes(
            Path(nadig_raw_dir) / "GSE264667_jurkat_raw_singlecell_01.h5ad",
            "nadig_jurkat",
            pools_map["nadig_jurkat"]["primary_barcodes"],
        )
    if "dld1" in pools_map:
        primary["dld1"] = _extract_10x_by_barcodes(
            Path(r10_raw_dir) / "dld1", "dld1", pools_map["dld1"]["primary_barcodes"]
        )
        if pools_map["dld1"].get("resample_barcodes"):
            resample["dld1"] = _extract_10x_by_barcodes(
                Path(r10_raw_dir) / "dld1", "dld1", pools_map["dld1"]["resample_barcodes"]
            )
    if "cd4t_rest" in pools_map:
        primary["cd4t_rest"] = _extract_10x_by_barcodes(
            Path(r10_raw_dir) / "cd4t_rest", "cd4t_rest", pools_map["cd4t_rest"]["primary_barcodes"]
        )
        if pools_map["cd4t_rest"].get("resample_barcodes"):
            resample["cd4t_rest"] = _extract_10x_by_barcodes(
                Path(r10_raw_dir) / "cd4t_rest", "cd4t_rest", pools_map["cd4t_rest"]["resample_barcodes"]
            )
    if "gse295214_h1" in pools_map:
        primary["gse295214_h1"] = _extract_10x_by_barcodes(
            Path(gse295214_raw_dir), "gse295214_h1", pools_map["gse295214_h1"]["primary_barcodes"]
        )
        if pools_map["gse295214_h1"].get("resample_barcodes"):
            resample["gse295214_h1"] = _extract_10x_by_barcodes(
                Path(gse295214_raw_dir), "gse295214_h1", pools_map["gse295214_h1"]["resample_barcodes"]
            )
    if "a549_ko" in pools_map:
        req_a549: Dict[str, Sequence[str]] = {"primary": pools_map["a549_ko"]["primary_barcodes"]}
        if pools_map["a549_ko"].get("resample_barcodes"):
            req_a549["resample"] = pools_map["a549_ko"]["resample_barcodes"]
        ext_a549 = _extract_a549_multi_pools_by_barcodes(Path(r12_raw_dir) / "a549_ko", req_a549)
        primary["a549_ko"] = ext_a549["primary"]
        primary["a549_ko"].pool_name = "a549_ko"
        if "resample" in ext_a549:
            resample["a549_ko"] = ext_a549["resample"]
            resample["a549_ko"].pool_name = "a549_ko"
    if "feng_ipsc" in pools_map:
        req_feng: Dict[str, Sequence[str]] = {"primary": pools_map["feng_ipsc"]["primary_barcodes"]}
        if pools_map["feng_ipsc"].get("resample_barcodes"):
            req_feng["resample"] = pools_map["feng_ipsc"]["resample_barcodes"]
        ext_feng = _extract_feng_multi_pools_by_barcodes(Path(r11_raw_dir) / "feng_ipsc", req_feng)
        primary["feng_ipsc"] = ext_feng["primary"]
        primary["feng_ipsc"].pool_name = "feng_ipsc"
        if "resample" in ext_feng:
            resample["feng_ipsc"] = ext_feng["resample"]
            resample["feng_ipsc"].pool_name = "feng_ipsc"

    return primary, resample


# ==============================================================================
# Committed R16 Reference + STATE_SE Delta Bundle Loader (H-I3a / H-I3b Side)
# ==============================================================================


def _load_r16_observed_extras_group2(
    main_obs_path: Path,
    eval_perts: Sequence[str],
    eval_genes: Sequence[str],
    d_obs: np.ndarray,
    d_fct: np.ndarray,
) -> Dict[str, Any]:
    """Load `mu_target`, `mu_sources`, `e1_constmed`, and `n_ctrl_total` from Group 2 R16 `observed.npz` and verify parity."""
    with np.load(main_obs_path, allow_pickle=False) as oz:
        obs_p = [str(x) for x in oz["target_names"].tolist()]
        obs_g = [str(x) for x in oz["gene_names"].tolist()]
        if obs_p != list(eval_perts) or obs_g != list(eval_genes):
            raise ValueError(f"Observed target/gene order mismatch vs audit on {main_obs_path}")
        if float(np.max(np.abs(np.asarray(oz["delta_obs"], dtype=np.float64) - d_obs))) > 1e-12:
            raise ValueError(f"delta_obs mismatch vs audit on {main_obs_path}")
        if float(np.max(np.abs(np.asarray(oz["fct_delta_pred"], dtype=np.float64) - d_fct))) > 1e-12:
            raise ValueError(f"fct_delta_pred mismatch vs audit on {main_obs_path}")
        return {
            "mu_target": np.asarray(oz["mu_target"], dtype=np.float64),
            "mu_sources": np.asarray(oz["mu_sources"], dtype=np.float64),
            "e1_constmed": np.asarray(oz["e1_constmed_delta_pred"], dtype=np.float64),
            "n_ctrl_total": int(oz["n_ctrl_total"]),
            "ctrl_sampling_fraction": float(oz["ctrl_sampling_fraction"]),
        }


def load_committed_r16_reference_and_se(
    target_line: str,
    audit_results_dir: Path = DEFAULT_LOCKED_SE_AUDIT_DIR,
    main_r16_dir: Path = DEFAULT_LOCKED_MAIN_R16_DIR,
    main_sens_dir: Path = DEFAULT_LOCKED_MAIN_SENS_DIR,
    repo_root: Path = REPO_ROOT,
    load_observed_extras: bool = True,
) -> CommittedR16ReferenceBundle:
    """
    Load committed R16 ground-truth matrices AND `STATE_SE` (`fixed_index_zero`) deltas
    for `target_line`:
      - Group 2 (`a549_ko`, `dld1`, `cd4t_rest`) from `round16/audit/results/STATE_SE_ensemble4__<line>.npz`
        plus `round16/results/arrays/<line>/observed.npz` (`mu_target`, `mu_sources`, `e1_constmed`)
      - Group 3 (`gse295214_h1`) strictly from `round16/audit/results/group3/STATE_SE_ensemble4__gse295214_h1.npz`
        + R9/R10 source reconstruction (NEVER reads `round16/results/group3/`)
      - `feng_ipsc` from `round16/results/sensitivity/arrays/feng_ipsc/`
    """
    ref_shas: Dict[str, str] = {}
    if target_line == "feng_ipsc":
        fdir = Path(main_sens_dir) / "arrays" / "feng_ipsc"
        obs_p = fdir / "observed.npz"
        ref_shas["observed.npz"] = sha256_file(obs_p)
        with np.load(obs_p, allow_pickle=False) as oz:
            eval_genes = [str(x) for x in oz["gene_names"].tolist()]
            eval_perts = [str(x) for x in oz["target_names"].tolist()]
            t_gidx = np.asarray(oz["target_indices"], dtype=np.int64)
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
            n_ctrl_tot = int(oz["n_ctrl_total"])
            samp_frac = float(oz["ctrl_sampling_fraction"])
            m_bcs = [str(x) for x in oz["matched_ctrl_barcodes"].tolist()]
            r_bcs = [str(x) for x in oz["resample_ctrl_barcodes"].tolist()]

        m_fp = fdir / "preds_STATE_SE_fixed_index_zero_matched.npz"
        r_fp = fdir / "preds_STATE_SE_fixed_index_zero_resample.npz"
        ref_shas["matched"] = sha256_file(m_fp)
        ref_shas["resample"] = sha256_file(r_fp)
        with np.load(m_fp, allow_pickle=False) as mz:
            se_matched = np.asarray(mz["delta_pred"], dtype=np.float64)
        with np.load(r_fp, allow_pickle=False) as rz:
            se_resample = np.asarray(rz["delta_pred"], dtype=np.float64)

        se_mis: Dict[str, np.ndarray] = {}
        mis_bcs: Dict[str, List[str]] = {}
        for src in ADDB_SEEN_SOURCES + ADDB_UNSEEN_SOURCES["feng_ipsc"]:
            sfp = fdir / f"preds_STATE_SE_fixed_index_zero_mis_{src}.npz"
            ref_shas[f"mis_{src}"] = sha256_file(sfp)
            with np.load(sfp, allow_pickle=False) as sz:
                se_mis[src] = np.asarray(sz["delta_pred"], dtype=np.float64)
                mis_bcs[src] = [str(x) for x in sz["control_barcodes"].tolist()]

        return CommittedR16ReferenceBundle(
            target_line=target_line,
            role=ADDB_LINE_ROLES[target_line],
            eval_perts=eval_perts,
            eval_genes=eval_genes,
            target_gene_idx=t_gidx,
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
            mismatched_barcodes=mis_bcs,
            se_delta_matched=se_matched,
            se_delta_resample=se_resample,
            se_delta_mismatched=se_mis,
            mu_target=mu_t,
            mu_sources=mu_s,
            e1_constmed=e1_d,
            n_ctrl_total=n_ctrl_tot,
            ctrl_sampling_fraction=samp_frac,
            r16_reference_sha256=ref_shas,
            extra_fct_full_r9_6=None,
        )

    base = Path(audit_results_dir) / "group3" if target_line == "gse295214_h1" else Path(audit_results_dir)
    npz_path = base / f"STATE_SE_ensemble4__{target_line}.npz"
    if not npz_path.is_file():
        raise FileNotFoundError(f"Committed audit STATE_SE file missing: {npz_path}")
    ref_shas[f"STATE_SE_ensemble4__{target_line}.npz"] = sha256_file(npz_path)

    with np.load(npz_path, allow_pickle=False) as z:
        eval_genes = [str(x) for x in z["eval_genes"].tolist()]
        eval_perts = [str(x) for x in z["eval_perts"].tolist()]
        t_gidx = np.asarray(z["target_gene_idx"], dtype=np.int64)
        d_obs = np.asarray(z["delta_obs"], dtype=np.float64)
        ha = np.asarray(z["delta_half_a"], dtype=np.float64)
        hb = np.asarray(z["delta_half_b"], dtype=np.float64)
        hv = np.asarray(z["half_valid"], dtype=bool)
        tm = np.asarray(z["delta_transfer"], dtype=np.float64)
        gm = np.asarray(z["source_grand_mean"], dtype=np.float64)
        fct_d = np.asarray(z["delta_fct"], dtype=np.float64)
        mu_t = np.asarray(z["ctrl_mean_libnorm"], dtype=np.float64) if "ctrl_mean_libnorm" in z.files else None
        se_matched = np.asarray(z["delta_pred_matched"], dtype=np.float64)
        se_resample = np.asarray(z["delta_pred_resample"], dtype=np.float64)
        m_bcs = [str(x) for x in z["matched_barcodes"].tolist()]
        r_bcs = [str(x) for x in z["resample_barcodes"].tolist()]
        se_mis = {}
        mis_bcs = {}
        for src in ADDB_SEEN_SOURCES + ADDB_UNSEEN_SOURCES[target_line]:
            se_mis[src] = np.asarray(z[f"delta_pred_mismatched__{src}"], dtype=np.float64)
            mis_bcs[src] = [str(x) for x in z[f"mismatched_barcodes__{src}"].tolist()]

    mu_s: Optional[np.ndarray] = None
    e1_d: Optional[np.ndarray] = None
    n_ctrl_tot = len(m_bcs)
    samp_frac = 1.0
    extra_fct: Optional[Dict[str, Any]] = None

    if load_observed_extras and target_line == "gse295214_h1":
        if mu_t is None:
            raise KeyError("ctrl_mean_libnorm missing in audit Group 3 STATE_SE .npz")
        extras = _build_gse295214_observed_extras(
            eval_genes=eval_genes,
            eval_perts=eval_perts,
            target_gene_idx=t_gidx,
            mu_target=mu_t,
            repo_root=repo_root,
        )
        mu_s = extras["mu_sources"]
        e1_d = extras["e1_constmed"]
        n_ctrl_tot = extras["n_ctrl_total"]
        samp_frac = float(len(m_bcs)) / float(max(n_ctrl_tot, 1))
        extra_fct = extras["extra_fct_full_r9_6"]
    elif load_observed_extras:
        main_obs = Path(main_r16_dir) / "arrays" / target_line / "observed.npz"
        if main_obs.is_file():
            ref_shas["r16_observed.npz"] = sha256_file(main_obs)
            ex = _load_r16_observed_extras_group2(main_obs, eval_perts, eval_genes, d_obs, fct_d)
            if mu_t is not None and float(np.max(np.abs(ex["mu_target"] - mu_t))) > 1e-12:
                raise ValueError(f"mu_target mismatch vs audit ctrl_mean_libnorm on {target_line}")
            mu_t = ex["mu_target"]
            mu_s = ex["mu_sources"]
            e1_d = ex["e1_constmed"]
            n_ctrl_tot = ex["n_ctrl_total"]
            samp_frac = ex["ctrl_sampling_fraction"]

    return CommittedR16ReferenceBundle(
        target_line=target_line,
        role=ADDB_LINE_ROLES[target_line],
        eval_perts=eval_perts,
        eval_genes=eval_genes,
        target_gene_idx=t_gidx,
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
        mismatched_barcodes=mis_bcs,
        se_delta_matched=se_matched,
        se_delta_resample=se_resample,
        se_delta_mismatched=se_mis,
        mu_target=mu_t,
        mu_sources=mu_s,
        e1_constmed=e1_d,
        n_ctrl_total=int(n_ctrl_tot),
        ctrl_sampling_fraction=float(samp_frac),
        r16_reference_sha256=ref_shas,
        extra_fct_full_r9_6=extra_fct,
    )


# ==============================================================================
# Statistical Endpoints: H-I1_Occam, H-I2_Occam, P4, H-I3a, H-I3b, §4.3 Matrix
# ==============================================================================


def _compute_boot_vector_for_pred(
    delta_obs: np.ndarray,
    delta_pred: np.ndarray,
    delta_transfer: np.ndarray,
    target_gene_idx: Sequence[int],
    delta_half_a: Optional[np.ndarray],
    delta_half_b: Optional[np.ndarray],
    half_valid: Optional[np.ndarray],
    kept: List[Tuple[int, np.ndarray]],
    draws: np.ndarray,
    b_ceil: np.ndarray,
    pt_ceil: float,
    source_grand_mean: Optional[np.ndarray] = None,
) -> Tuple[float, np.ndarray]:
    """Compute exact disattenuated point estimate and `(n_boot,)` bootstrap vector sharing `draws`."""
    st = _project_streams(
        kept=kept,
        delta_obs=delta_obs,
        pred=delta_pred,
        t_mean=delta_transfer,
        source_grand_mean=source_grand_mean,
        half_a_obs=delta_half_a,
        half_b_obs=delta_half_b,
        half_valid=half_valid,
    )
    is_dorth = source_grand_mean is not None
    mom = _moments(st["xo_proj"], st["yo_proj"])
    pt_r = float(_moments_to_corr_and_slope(*mom, zero_degenerate_y_as_zero=is_dorth)[0])
    b_r = _moments_to_corr_and_slope(*mom, draws=draws, zero_degenerate_y_as_zero=is_dorth)[0]

    def _dis(rv: float, cv: float) -> float:
        if is_dorth and np.isfinite(rv) and abs(rv) <= 1e-12 and np.isfinite(cv) and cv > 1e-8:
            return 0.0
        return disattenuate(rv, cv)

    pt_cal = _dis(pt_r, pt_ceil)
    b_cal = np.asarray([_dis(float(rv), float(cv)) for rv, cv in zip(b_r, b_ceil)], dtype=np.float64)
    return pt_cal, b_cal


def _classify_hi3a_summary(contrast: Mapping[str, Any]) -> str:
    vf = float(contrast.get("valid_frac", 0.0))
    lo = float(contrast.get("ci975_lo", float("nan")))
    hi = float(contrast.get("ci975_hi", float("nan")))
    if not (np.isfinite(vf) and vf >= 0.95 and np.isfinite(lo) and np.isfinite(hi)):
        return "INDISTINGUISHABLE"
    if lo == 0.0 or hi == 0.0:
        return "INDISTINGUISHABLE"
    if lo > 0.0:
        return "OCCAM_SUPERIOR_TO_STATE_SE"
    if hi < 0.0:
        return "STATE_SE_SUPERIOR_TO_OCCAM"
    return "INDISTINGUISHABLE"


def _classify_hi3b_summary(contrast: Mapping[str, Any]) -> str:
    vf = float(contrast.get("valid_frac", 0.0))
    lo = float(contrast.get("ci975_lo", float("nan")))
    hi = float(contrast.get("ci975_hi", float("nan")))
    if not (np.isfinite(vf) and vf >= 0.95 and np.isfinite(lo) and np.isfinite(hi)):
        return "INDISTINGUISHABLE"
    if lo == 0.0 or hi == 0.0:
        return "INDISTINGUISHABLE"
    if lo > 0.0:
        return "OCCAM_MORE_CONTROL_SPECIFIC"
    if hi < 0.0:
        return "STATE_SE_MORE_CONTROL_SPECIFIC"
    return "INDISTINGUISHABLE"


def compute_hi3a_and_hi3b_line(
    ref: CommittedR16ReferenceBundle,
    occam_matched: np.ndarray,
    occam_resample: np.ndarray,
    occam_mismatched: Mapping[str, np.ndarray],
    n_boot: int = 4000,
    seed: int = DEFAULT_R16_SEED,
    source_grand_mean: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """
    Compute pre-registered paired `STATE_SE` vs `Occam` contrasts (`PREREG_I_ADDENDUM_B §4.2`):
      - `H-I3a`: matched accuracy `r(Occam_matched) - r(STATE_SE_matched)`
      - `H-I3b`: per-source and per-stratum (`seen_by_checkpoint`, `unseen_by_checkpoint`)
        specificity contrast `Delta_spec = Delta_pair(Occam) - Delta_pair(STATE_SE)`, plus
        resample-adjusted `Delta_spec_excess`.
    """
    hi3a = compute_paired_contrast_r16(
        delta_obs=ref.delta_obs,
        pred_a=occam_matched,
        pred_b=ref.se_delta_matched,
        delta_transfer=ref.delta_transfer,
        target_gene_idx=ref.target_gene_idx,
        source_grand_mean=source_grand_mean,
        delta_half_a=ref.delta_half_a,
        delta_half_b=ref.delta_half_b,
        half_valid=ref.half_valid,
        mode="all_nontarget",
        n_boot=n_boot,
        seed=seed,
    )
    hi3a["classification"] = _classify_hi3a_summary(hi3a)
    hi3a["occam_summary"] = compute_orth_suite_r16(
        ref.delta_obs, occam_matched, ref.delta_transfer, ref.target_gene_idx,
        source_grand_mean=source_grand_mean, delta_half_a=ref.delta_half_a,
        delta_half_b=ref.delta_half_b, half_valid=ref.half_valid,
        mode="all_nontarget", n_boot=n_boot, seed=seed,
    )
    hi3a["state_se_summary"] = compute_orth_suite_r16(
        ref.delta_obs, ref.se_delta_matched, ref.delta_transfer, ref.target_gene_idx,
        source_grand_mean=source_grand_mean, delta_half_a=ref.delta_half_a,
        delta_half_b=ref.delta_half_b, half_valid=ref.half_valid,
        mode="all_nontarget", n_boot=n_boot, seed=seed,
    )

    kept = _extract_pert_gene_indices(
        delta_obs=ref.delta_obs,
        t_mean=ref.delta_transfer,
        pred=occam_matched,
        target_cols=ref.target_gene_idx,
        mode="all_nontarget",
        extra_pred=ref.se_delta_matched,
        source_grand_mean=source_grand_mean,
    )
    m = len(kept)
    if m == 0:
        return {"h_i3a": hi3a, "h_i3b": {"per_source": {}, "per_stratum_mean": {}, "n_eval_perts": 0}}

    st_ref = _project_streams(
        kept=kept,
        delta_obs=ref.delta_obs,
        pred=occam_matched,
        t_mean=ref.delta_transfer,
        source_grand_mean=source_grand_mean,
        half_a_obs=ref.delta_half_a,
        half_b_obs=ref.delta_half_b,
        half_valid=ref.half_valid,
    )
    rng = np.random.default_rng(int(seed))
    draws = rng.integers(0, m, size=(int(n_boot), m))
    if st_ref["any_hv"]:
        mh = _moments(st_ref["ha_proj"], st_ref["hb_proj"])
        pt_rho = float(_moments_to_corr_and_slope(*mh)[0])
        pt_ceil = spearman_brown_ceiling(pt_rho)
        b_rho = _moments_to_corr_and_slope(*mh, draws=draws)[0]
        b_ceil = np.asarray([spearman_brown_ceiling(float(v)) for v in b_rho], dtype=np.float64)
    else:
        pt_ceil = float("nan")
        b_ceil = np.full(int(n_boot), float("nan"), dtype=np.float64)

    def _vec(p: np.ndarray) -> Tuple[float, np.ndarray]:
        return _compute_boot_vector_for_pred(
            ref.delta_obs, p, ref.delta_transfer, ref.target_gene_idx,
            ref.delta_half_a, ref.delta_half_b, ref.half_valid,
            kept=kept, draws=draws, b_ceil=b_ceil, pt_ceil=pt_ceil,
            source_grand_mean=source_grand_mean,
        )

    pt_om, b_om = _vec(occam_matched)
    pt_or, b_or = _vec(occam_resample)
    pt_sm, b_sm = _vec(ref.se_delta_matched)
    pt_sr, b_sr = _vec(ref.se_delta_resample)

    all_sources = sorted(set(ref.seen_sources) | set(ref.unseen_sources))
    per_source_spec: Dict[str, Any] = {}
    spec_pt_map: Dict[str, float] = {}
    spec_b_map: Dict[str, np.ndarray] = {}
    ex_pt_map: Dict[str, float] = {}
    ex_b_map: Dict[str, np.ndarray] = {}
    occ_dp_map: Dict[str, float] = {}
    se_dp_map: Dict[str, float] = {}

    for src in all_sources:
        pt_os, b_os = _vec(occam_mismatched[src])
        pt_ss, b_ss = _vec(ref.se_delta_mismatched[src])
        dp_o_pt = pt_om - pt_os
        dp_s_pt = pt_sm - pt_ss
        dp_o_b = b_om - b_os
        dp_s_b = b_sm - b_ss

        spec_pt = dp_o_pt - dp_s_pt
        spec_b = dp_o_b - dp_s_b
        spec_sum = summarize_boot_r16(spec_pt, spec_b, m, int(n_boot))
        spec_sum["classification"] = _classify_hi3b_summary(spec_sum)
        spec_sum["occam_delta_pair"] = float(dp_o_pt)
        spec_sum["state_se_delta_pair"] = float(dp_s_pt)

        ex_o_pt = dp_o_pt - (pt_om - pt_or)
        ex_s_pt = dp_s_pt - (pt_sm - pt_sr)
        ex_pt = ex_o_pt - ex_s_pt
        ex_b = (dp_o_b - (b_om - b_or)) - (dp_s_b - (b_sm - b_sr))
        ex_sum = summarize_boot_r16(ex_pt, ex_b, m, int(n_boot))
        ex_sum["classification"] = _classify_hi3b_summary(ex_sum)
        ex_sum["occam_delta_excess"] = float(ex_o_pt)
        ex_sum["state_se_delta_excess"] = float(ex_s_pt)

        stratum = "seen_by_checkpoint" if src in ref.seen_sources else "unseen_by_checkpoint"
        per_source_spec[src] = {
            "stratum": stratum,
            "delta_spec": spec_sum,
            "delta_spec_excess": ex_sum,
        }
        spec_pt_map[src] = spec_pt
        spec_b_map[src] = spec_b
        ex_pt_map[src] = ex_pt
        ex_b_map[src] = ex_b
        occ_dp_map[src] = dp_o_pt
        se_dp_map[src] = dp_s_pt

    per_stratum_mean: Dict[str, Any] = {}
    for sname, ssources in (("seen_by_checkpoint", ref.seen_sources), ("unseen_by_checkpoint", ref.unseen_sources)):
        avail = [s for s in ssources if s in spec_b_map]
        if not avail:
            per_stratum_mean[sname] = {"n_sources": 0, "sources": []}
            continue
        mean_pt = float(np.mean([spec_pt_map[s] for s in avail]))
        mean_b = np.mean(np.stack([spec_b_map[s] for s in avail], axis=0), axis=0)
        msum = summarize_boot_r16(mean_pt, mean_b, m, int(n_boot))
        msum["classification"] = _classify_hi3b_summary(msum)

        ex_mpt = float(np.mean([ex_pt_map[s] for s in avail]))
        ex_mb = np.mean(np.stack([ex_b_map[s] for s in avail], axis=0), axis=0)
        mex_sum = summarize_boot_r16(ex_mpt, ex_mb, m, int(n_boot))
        mex_sum["classification"] = _classify_hi3b_summary(mex_sum)
        per_stratum_mean[sname] = {
            "n_sources": len(avail),
            "sources": avail,
            "occam_stratum_mean": float(np.mean([occ_dp_map[s] for s in avail])),
            "state_se_stratum_mean": float(np.mean([se_dp_map[s] for s in avail])),
            "delta_spec": msum,
            "delta_spec_excess": mex_sum,
        }

    return {
        "h_i3a": hi3a,
        "h_i3b": {
            "n_eval_perts": int(m),
            "per_source": per_source_spec,
            "per_stratum_mean": per_stratum_mean,
        },
    }


def classify_advance_outcome_matrix(
    se_h_i1: Mapping[str, Any],
    occam_h_i1: Mapping[str, Any],
    occam_h_i2: Mapping[str, Any],
    hi3a: Mapping[str, Any],
    hi3b: Mapping[str, Any],
    p4_any_verdict_change: bool = False,
) -> Dict[str, Any]:
    """
    Classify the joint `(STATE_SE, Occam, FCT)` result against the Pre-Locked Advance
    Outcome Interpretation Matrix (`PREREG_I_ADDENDUM_B_DRAFT_v2.md §4.3 Cases 1–6`).
    """
    def _stratum_cat(h1: Mapping[str, Any]) -> str:
        lsum = str(h1.get("line_summary", ""))
        counts = h1.get("counts", {})
        if lsum in ("INCONCLUSIVE", "MISSING_STRATUM", "MIXED_OR_INCONCLUSIVE") or counts.get("INCONCLUSIVE", 0) > 0:
            return "INCONCLUSIVE"
        if lsum == "LINE_MISMATCHED_BETTER":
            return "MISMATCHED_BETTER"
        rc = h1.get("resample_reference") or h1.get("resample_control") or {}
        if lsum == "LINE_MATCHED_SPECIFIC":
            if rc and not rc.get("any_source_exceeds_resample", True):
                return "SUB_RESAMPLE"
            return "MATCHED_SPECIFIC"
        return "EQUIVALENT_OR_SUB_DELTA"

    se_cat = _stratum_cat(se_h_i1)
    oc_cat = _stratum_cat(occam_h_i1)
    h3a_cls = str(hi3a.get("classification", "INDISTINGUISHABLE"))

    if se_cat in ("INCONCLUSIVE", "MISMATCHED_BETTER") or oc_cat in ("INCONCLUSIVE", "MISMATCHED_BETTER"):
        case_id = "Case_6_Power_Limited_Or_OOD"
    elif se_cat == "MATCHED_SPECIFIC" and oc_cat == "MATCHED_SPECIFIC":
        case_id = "Case_5_Cross_Architecture_Control_Utilization"
    elif oc_cat == "MATCHED_SPECIFIC" and se_cat != "MATCHED_SPECIFIC":
        case_id = "Case_2_Occam_Control_Sensitive_SE_Not"
    elif se_cat == "MATCHED_SPECIFIC" and oc_cat != "MATCHED_SPECIFIC":
        case_id = "Case_3_STATE_SE_Control_Sensitive_Occam_Not"
    elif h3a_cls == "OCCAM_SUPERIOR_TO_STATE_SE":
        case_id = "Case_4_Occam_Prior_Accuracy_Gain_Without_Control_Specificity"
    else:
        case_id = "Case_1_Cross_Architecture_Control_Underutilization"

    return {
        "matrix_case": case_id,
        "state_se_h_i1_bucket": se_cat,
        "occam_h_i1_bucket": oc_cat,
        "occam_h_i2_call": occam_h_i2.get("classification"),
        "h_i3a_call": h3a_cls,
        "p4_verdict_changed": bool(p4_any_verdict_change),
    }


def _evaluate_arm_suite_addb(
    ref: CommittedR16ReferenceBundle,
    pred_matched: np.ndarray,
    pred_resample: Optional[np.ndarray],
    pred_mismatched: Optional[Mapping[str, np.ndarray]],
    source_grand_mean: Optional[np.ndarray],
    n_boot: int,
    seed: int,
    h_i1_key: str,
    h_i2_key: str,
    matched_key: str,
) -> Dict[str, Any]:
    """Evaluate matched / mismatched / FCT / H-I1 / H-I2 / resample suite for one arm convention."""
    d_obs = ref.delta_obs
    ha = ref.delta_half_a
    hb = ref.delta_half_b
    hv = ref.half_valid
    tm = ref.delta_transfer
    tg = ref.target_gene_idx
    fct = ref.delta_fct

    orth_m = compute_orth_suite_r16(
        d_obs, pred_matched, tm, tg, source_grand_mean=source_grand_mean,
        delta_half_a=ha, delta_half_b=hb, half_valid=hv,
        mode="all_nontarget", n_boot=n_boot, seed=seed,
    )
    vs_fct = compute_paired_contrast_r16(
        d_obs, pred_matched, fct, tm, tg, source_grand_mean=source_grand_mean,
        delta_half_a=ha, delta_half_b=hb, half_valid=hv,
        mode="all_nontarget", n_boot=n_boot, seed=seed,
    )
    h_i2 = adjudicate_claim_i2_line(vs_fct)

    out: Dict[str, Any] = {
        matched_key: orth_m,
        "Occam_vs_FCT": vs_fct,
    }

    if pred_mismatched is not None:
        paired_mis = compute_paired_model_vs_all_mismatched_r16(
            d_obs, pred_matched, dict(pred_mismatched), tm, tg,
            source_grand_mean=source_grand_mean,
            delta_half_a=ha, delta_half_b=hb, half_valid=hv,
            mode="all_nontarget", n_boot=n_boot, seed=seed,
        )
        pos_ci975 = bool(
            np.isfinite(orth_m.get("valid_frac", 0.0))
            and float(orth_m["valid_frac"]) >= MIN_BOOT_VALID_FRAC
            and np.isfinite(orth_m.get("rho_ab_orth", float("nan")))
            and float(orth_m["rho_ab_orth"]) > 0.0
            and np.isfinite(orth_m.get("ci975_lo", float("nan")))
            and float(orth_m["ci975_lo"]) > 0.0
        )
        res_ref = (
            compute_paired_diff_vs_resample_r16(
                d_obs, pred_matched, pred_resample, dict(pred_mismatched), tm, tg,
                source_grand_mean=source_grand_mean,
                delta_half_a=ha, delta_half_b=hb, half_valid=hv,
                mode="all_nontarget", n_boot=n_boot, seed=seed,
                matched_positive_ci975=pos_ci975,
            )
            if pred_resample is not None
            else None
        )
        h_i1 = adjudicate_claim_i1_line(
            orth_m, paired_mis, ref.seen_sources, ref.unseen_sources,
            delta=DEFAULT_TOST_DELTA, resample_comparison=res_ref,
        )
        # Mechanical M6.2 anti-masking indicator: winning/tying FCT on weak-FCT or specificity-failing line
        fct_pt = float(compute_orth_suite_r16(
            d_obs, fct, tm, tg, source_grand_mean=source_grand_mean,
            delta_half_a=ha, delta_half_b=hb, half_valid=hv,
            mode="all_nontarget", n_boot=0, seed=seed,
        ).get("estimate", float("nan")))
        m62_triggered = bool(
            h_i2["classification"] != "FCT_SUPERIOR_TO_MODEL"
            and (
                h_i1["line_summary"] != "LINE_MATCHED_SPECIFIC"
                or ref.target_line == "cd4t_rest"
                or (np.isfinite(fct_pt) and fct_pt < 0.10)
            )
        )
        h_i2_full = {**h_i2, "anti_masking_m62_triggered": m62_triggered}
        res_block = dict(res_ref) if res_ref is not None else None
        if res_block is not None:
            res_block["within_line_resample_status"] = res_block.get("overall_status")
            res_block["mandatory_narrative_rule"] = res_block.get("overall_interpretation")
        out["Occam_paired_vs_mismatched"] = paired_mis
        out[h_i1_key] = h_i1
        out[h_i2_key] = h_i2_full
        out["within_line_resample_reference"] = res_block
    else:
        h_i2_full = {**h_i2, "anti_masking_m62_triggered": False}
        out[h_i2_key] = h_i2_full
        out["practical_accuracy"] = compute_practical_metrics_r16(d_obs, pred_matched, tg)

    return out


def _rollup_common_axis_h_i1(
    primary_orth_matched: Mapping[str, Any],
    pair_calls_by_src: Mapping[str, Mapping[str, Any]],
    seen_sources: Sequence[str],
    unseen_sources: Sequence[str],
    delta: float = DEFAULT_TOST_DELTA,
) -> Dict[str, Any]:
    """
    CB1.2 + Amendment CB1.2a roll-up for `H_I1_Occam_common_axis`:
      - Per-pair categories come from `pair_calls_by_src` (gated by validity-only of each
        pair's own masked matched arm, or reused primary category when `n_union_absent == 0`);
      - Line-level `matched_arm_valid` / `matched_positive_ci975` gate uses the primary
        `rhj_6983_axis` matched arm (`primary_orth_matched`);
      - Exact `PREREG_I` 8-state cascade applied across strata.
    """
    seen_set = set(seen_sources)
    unseen_set = set(unseen_sources)
    if set(pair_calls_by_src.keys()) != (seen_set | unseen_set):
        raise ValueError("pair_calls_by_src keys must match seen | unseen exactly")

    vf = float(primary_orth_matched.get("valid_frac", 0.0))
    rho_orth = float(primary_orth_matched.get("rho_ab_orth", float("nan")))
    lo975 = float(primary_orth_matched.get("ci975_lo", float("nan")))
    hi975 = float(primary_orth_matched.get("ci975_hi", float("nan")))
    matched_arm_valid = bool(
        np.isfinite(vf)
        and vf >= MIN_BOOT_VALID_FRAC
        and np.isfinite(rho_orth)
        and rho_orth > 0.0
        and np.isfinite(lo975)
        and np.isfinite(hi975)
    )
    matched_positive_ci975 = bool(matched_arm_valid and lo975 > 0.0)

    def _stratum_from_calls(sources: Sequence[str]) -> Dict[str, Any]:
        counts = {
            "MATCHED_BETTER": 0,
            "EQUIVALENT": 0,
            "EXACT_INVARIANCE": 0,
            "INCONCLUSIVE": 0,
            "MISMATCHED_BETTER": 0,
        }
        p_calls: Dict[str, Dict[str, Any]] = {}
        for sk in sorted(sources):
            c = dict(pair_calls_by_src[sk])
            p_calls[sk] = c
            counts[str(c["category"])] += 1
        n_pairs = len(p_calls)
        all_w = bool(n_pairs > 0 and all(bool(c["within_delta"]) for c in p_calls.values()))
        any_w = bool(n_pairs > 0 and any(bool(c["within_delta"]) for c in p_calls.values()))
        eq_like = counts["EQUIVALENT"] + counts["EXACT_INVARIANCE"]
        if n_pairs == 0:
            st = "STRATUM_EMPTY"
        elif not matched_arm_valid:
            st = "STRATUM_INCONCLUSIVE"
        elif matched_positive_ci975 and counts["MATCHED_BETTER"] == n_pairs:
            st = "STRATUM_MATCHED_SPECIFIC_SUB_DELTA" if all_w else "STRATUM_MATCHED_SPECIFIC"
        elif eq_like == n_pairs:
            st = "STRATUM_EQUIVALENT_NONZERO_SIGNAL" if matched_positive_ci975 else "STRATUM_EQUIVALENT_NULL_SIGNAL"
        elif counts["MISMATCHED_BETTER"] == n_pairs:
            st = "STRATUM_MISMATCHED_BETTER"
        else:
            st = "STRATUM_MIXED_OR_INCONCLUSIVE"
        return {
            "n_pairs": n_pairs,
            "sources": sorted(sources),
            "counts": counts,
            "equiv_like_count": int(eq_like),
            "all_within_delta": all_w,
            "any_within_delta": any_w,
            "per_source": p_calls,
            "stratum_summary": st,
        }

    seen_res = _stratum_from_calls(seen_sources)
    unseen_res = _stratum_from_calls(unseen_sources)
    combined_counts = {
        k: seen_res["counts"][k] + unseen_res["counts"][k]
        for k in ("MATCHED_BETTER", "EQUIVALENT", "EXACT_INVARIANCE", "INCONCLUSIVE", "MISMATCHED_BETTER")
    }
    combined_per_source: Dict[str, Dict[str, Any]] = {}
    for sk, val in seen_res["per_source"].items():
        combined_per_source[sk] = {**val, "stratum": "seen_by_checkpoint"}
    for sk, val in unseen_res["per_source"].items():
        combined_per_source[sk] = {**val, "stratum": "unseen_by_checkpoint"}

    n_seen = seen_res["n_pairs"]
    n_unseen = unseen_res["n_pairs"]
    n_total = n_seen + n_unseen
    all_within_delta = bool(n_total > 0 and all(bool(v["within_delta"]) for v in combined_per_source.values()))
    any_within_delta = bool(n_total > 0 and any(bool(v["within_delta"]) for v in combined_per_source.values()))

    if n_seen == 0 or n_unseen == 0:
        line_summary = "MISSING_STRATUM"
    elif not matched_arm_valid:
        line_summary = "INCONCLUSIVE"
    elif matched_positive_ci975 and seen_res["counts"]["MATCHED_BETTER"] == n_seen and unseen_res["counts"]["MATCHED_BETTER"] == n_unseen:
        line_summary = "LINE_MATCHED_SPECIFIC_SUB_DELTA" if all_within_delta else "LINE_MATCHED_SPECIFIC"
    elif seen_res["equiv_like_count"] == n_seen and unseen_res["equiv_like_count"] == n_unseen:
        line_summary = "LINE_EQUIVALENT_NONZERO_SIGNAL" if matched_positive_ci975 else "LINE_EQUIVALENT_NULL_SIGNAL"
    elif seen_res["counts"]["MISMATCHED_BETTER"] == n_seen and unseen_res["counts"]["MISMATCHED_BETTER"] == n_unseen:
        line_summary = "LINE_MISMATCHED_BETTER"
    else:
        line_summary = "MIXED_OR_INCONCLUSIVE"

    interpretation_note = (
        "能检测到，但小于事先规定的最小有意义效应 δ=0.04（须与同线重抽参照臂结果合并报告，不可单独写成‘模型用到了对照’）"
        if line_summary == "LINE_MATCHED_SPECIFIC_SUB_DELTA"
        else None
    )
    return {
        "matched_arm_valid": matched_arm_valid,
        "matched_positive_ci975": matched_positive_ci975,
        "matched_rho_ab_orth": rho_orth,
        "matched_valid_frac": vf,
        "n_pairs": n_total,
        "delta": float(delta),
        "all_within_delta": all_within_delta,
        "any_within_delta": any_within_delta,
        "strata": {"seen_by_checkpoint": seen_res, "unseen_by_checkpoint": unseen_res},
        "counts": combined_counts,
        "per_source": combined_per_source,
        "line_summary": line_summary,
        "interpretation_note": interpretation_note,
    }


def _evaluate_p4_common_axis_block(
    ref: CommittedR16ReferenceBundle,
    p4_common_pairs: Mapping[str, Any],
    primary_h_i1: Mapping[str, Any],
    n_boot: int,
    seed: int,
    primary_orth_matched: Optional[Mapping[str, Any]] = None,
) -> Tuple[Dict[str, Any], bool]:
    """Evaluate mandatory P4 `pairwise_common_axis` contrasts, CB1.2 line H-I1 rollup, and H-I3b crosscheck."""
    d_obs = ref.delta_obs
    ha = ref.delta_half_a
    hb = ref.delta_half_b
    hv = ref.half_valid
    tm = ref.delta_transfer
    tg = ref.target_gene_idx

    per_src_out: Dict[str, Any] = {}
    pair_calls_for_h1: Dict[str, Dict[str, Any]] = {}
    any_changed = False

    # H-I3b common-axis setup (shared draws on all_nontarget)
    kept_common = _extract_pert_gene_indices(
        delta_obs=d_obs, t_mean=tm, pred=ref.se_delta_matched, target_cols=tg, mode="all_nontarget",
    )
    m_k = len(kept_common)
    hi3b_per_src: Dict[str, Any] = {}
    hi3b_per_stratum: Dict[str, Any] = {}
    spec_pt_map: Dict[str, float] = {}
    spec_b_map: Dict[str, np.ndarray] = {}
    occ_dp_map: Dict[str, float] = {}
    se_dp_map: Dict[str, float] = {}

    if m_k > 0:
        st_ref = _project_streams(
            kept=kept_common, delta_obs=d_obs, pred=ref.se_delta_matched, t_mean=tm,
            source_grand_mean=None, half_a_obs=ha, half_b_obs=hb, half_valid=hv,
        )
        draws = np.random.default_rng(int(seed)).integers(0, m_k, size=(int(n_boot), m_k))
        if st_ref["any_hv"]:
            mh = _moments(st_ref["ha_proj"], st_ref["hb_proj"])
            pt_ceil = spearman_brown_ceiling(float(_moments_to_corr_and_slope(*mh)[0]))
            b_rho = _moments_to_corr_and_slope(*mh, draws=draws)[0]
            b_ceil = np.asarray([spearman_brown_ceiling(float(v)) for v in b_rho], dtype=np.float64)
        else:
            pt_ceil = float("nan")
            b_ceil = np.full(int(n_boot), float("nan"), dtype=np.float64)
        pt_sm, b_sm = _compute_boot_vector_for_pred(
            d_obs, ref.se_delta_matched, tm, tg, ha, hb, hv,
            kept=kept_common, draws=draws, b_ceil=b_ceil, pt_ceil=pt_ceil, source_grand_mean=None,
        )
    else:
        draws = np.zeros((0, 0), dtype=np.int64)
        pt_ceil = float("nan")
        b_ceil = np.zeros(0, dtype=np.float64)
        pt_sm, b_sm = float("nan"), np.zeros(0, dtype=np.float64)

    for src in sorted(p4_common_pairs.keys()):
        entry = p4_common_pairs[src]
        if isinstance(entry, tuple):
            p_m, p_s, n_c = entry
            n_union_abs = max(0, 6983 - int(n_c))
            reused = bool(n_union_abs == 0)
            mask_sha = ""
        else:
            p_m = entry["matched"]
            p_s = entry["mismatched"]
            n_c = int(entry["n_common"])
            n_union_abs = int(entry.get("n_union_absent", max(0, 6983 - n_c)))
            reused = bool(entry.get("reused_primary_rhj_6983", n_union_abs == 0))
            mask_sha = str(entry.get("mask_sha256", ""))

        m_orth = compute_orth_suite_r16(
            d_obs, p_m, tm, tg, source_grand_mean=None,
            delta_half_a=ha, delta_half_b=hb, half_valid=hv,
            mode="all_nontarget", n_boot=n_boot, seed=seed,
        )
        s_orth = compute_orth_suite_r16(
            d_obs, p_s, tm, tg, source_grand_mean=None,
            delta_half_a=ha, delta_half_b=hb, half_valid=hv,
            mode="all_nontarget", n_boot=n_boot, seed=seed,
        )
        pc = compute_paired_contrast_r16(
            d_obs, p_m, p_s, tm, tg, source_grand_mean=None,
            delta_half_a=ha, delta_half_b=hb, half_valid=hv,
            mode="all_nontarget", n_boot=n_boot, seed=seed,
        )
        prim_pair = primary_h_i1.get("per_source", {}).get(src) or {}
        prim_cat = prim_pair.get("category")
        if reused and prim_pair:
            # CB1.2(c): pairs reusing primary axis take their primary category
            tost = {k: v for k, v in prim_pair.items() if k != "stratum"}
            tost["masked_matched_arm_gate_passed"] = bool(primary_h_i1.get("matched_arm_valid", False))
        else:
            # Amendment CB1.2(a): INCONCLUSIVE only if this pair's own masked matched arm fails
            # validity (valid_frac >= 0.95, rho_ab_orth > 0, finite 97.5% bounds); positivity is
            # checked only at the line level via CB1.2(b) on primary rhj_6983_axis.
            tost = classify_pair_tost_r16(pc, delta=DEFAULT_TOST_DELTA)
            vf_m = float(m_orth.get("valid_frac", 0.0))
            rho_m = float(m_orth.get("rho_ab_orth", float("nan")))
            lo975_m = float(m_orth.get("ci975_lo", float("nan")))
            hi975_m = float(m_orth.get("ci975_hi", float("nan")))
            pair_m_ok = bool(
                np.isfinite(vf_m)
                and vf_m >= MIN_BOOT_VALID_FRAC
                and np.isfinite(rho_m)
                and rho_m > 0.0
                and np.isfinite(lo975_m)
                and np.isfinite(hi975_m)
            )
            if not pair_m_ok:
                tost["category"] = "INCONCLUSIVE"
                tost["within_delta"] = False
            tost["masked_matched_arm_gate_passed"] = pair_m_ok

        changed = bool(prim_cat is not None and tost["category"] != prim_cat)
        if changed:
            any_changed = True
        pair_calls_for_h1[src] = tost
        per_src_out[src] = {
            "n_common": int(n_c),
            "n_union_absent": int(n_union_abs),
            "reused_primary_rhj_6983": bool(reused),
            "mask_sha256": mask_sha,
            "matched_orth": m_orth,
            "mismatched_orth": s_orth,
            "paired_delta": pc,
            "pair_classification": tost,
            "primary_pair_classification": prim_cat,
            "verdict_changed_vs_primary": changed,
        }

        if m_k > 0 and src in ref.se_delta_mismatched:
            pt_pm, b_pm = _compute_boot_vector_for_pred(
                d_obs, p_m, tm, tg, ha, hb, hv,
                kept=kept_common, draws=draws, b_ceil=b_ceil, pt_ceil=pt_ceil, source_grand_mean=None,
            )
            pt_ps, b_ps = _compute_boot_vector_for_pred(
                d_obs, p_s, tm, tg, ha, hb, hv,
                kept=kept_common, draws=draws, b_ceil=b_ceil, pt_ceil=pt_ceil, source_grand_mean=None,
            )
            pt_ss, b_ss = _compute_boot_vector_for_pred(
                d_obs, ref.se_delta_mismatched[src], tm, tg, ha, hb, hv,
                kept=kept_common, draws=draws, b_ceil=b_ceil, pt_ceil=pt_ceil, source_grand_mean=None,
            )
            dp_o_pt = pt_pm - pt_ps
            dp_s_pt = pt_sm - pt_ss
            spec_pt = dp_o_pt - dp_s_pt
            spec_b = (b_pm - b_ps) - (b_sm - b_ss)
            ssum = summarize_boot_r16(spec_pt, spec_b, m_k, int(n_boot))
            ssum["classification"] = _classify_hi3b_summary(ssum)
            ssum["occam_delta_pair_common"] = float(dp_o_pt)
            ssum["state_se_delta_pair"] = float(dp_s_pt)
            stratum = "seen_by_checkpoint" if src in ref.seen_sources else "unseen_by_checkpoint"
            hi3b_per_src[src] = {"stratum": stratum, "delta_spec": ssum}
            spec_pt_map[src] = spec_pt
            spec_b_map[src] = spec_b
            occ_dp_map[src] = dp_o_pt
            se_dp_map[src] = dp_s_pt

    for sname, ssources in (("seen_by_checkpoint", ref.seen_sources), ("unseen_by_checkpoint", ref.unseen_sources)):
        avail = [s for s in ssources if s in spec_b_map]
        if not avail:
            hi3b_per_stratum[sname] = {"n_sources": 0, "sources": []}
            continue
        mean_pt = float(np.mean([spec_pt_map[s] for s in avail]))
        mean_b = np.mean(np.stack([spec_b_map[s] for s in avail], axis=0), axis=0)
        msum = summarize_boot_r16(mean_pt, mean_b, m_k, int(n_boot))
        msum["classification"] = _classify_hi3b_summary(msum)
        hi3b_per_stratum[sname] = {
            "n_sources": len(avail),
            "sources": avail,
            "occam_stratum_mean_common": float(np.mean([occ_dp_map[s] for s in avail])),
            "state_se_stratum_mean": float(np.mean([se_dp_map[s] for s in avail])),
            "delta_spec": msum,
        }

    # CB1.2(b)-(c): Line-level H_I1_Occam_common_axis roll-up using primary rhj_6983_axis matched arm
    expected_keys = set(ref.seen_sources) | set(ref.unseen_sources)
    h_i1_common: Optional[Dict[str, Any]] = None
    line_changed = False
    if set(pair_calls_for_h1.keys()) == expected_keys:
        prim_m_gate = (
            primary_orth_matched
            if primary_orth_matched is not None
            else {
                "valid_frac": primary_h_i1.get("matched_valid_frac", 0.0),
                "rho_ab_orth": primary_h_i1.get("matched_rho_ab_orth", float("nan")),
                "ci975_lo": 1.0 if primary_h_i1.get("matched_positive_ci975") else 0.0,
                "ci975_hi": 1.0 if primary_h_i1.get("matched_positive_ci975") else 0.0,
            }
        )
        h_i1_common = _rollup_common_axis_h_i1(
            primary_orth_matched=prim_m_gate,
            pair_calls_by_src=pair_calls_for_h1,
            seen_sources=ref.seen_sources,
            unseen_sources=ref.unseen_sources,
            delta=DEFAULT_TOST_DELTA,
        )
        line_changed = bool(h_i1_common["line_summary"] != primary_h_i1.get("line_summary"))

    any_flagged = bool(any_changed or line_changed)
    return {
        "per_source": per_src_out,
        "H_I1_Occam_common_axis": h_i1_common,
        "H_I3b_common_axis_crosscheck": {
            "n_eval_perts": int(m_k),
            "per_source": hi3b_per_src,
            "per_stratum_mean": hi3b_per_stratum,
        },
        "pair_verdict_changed_vs_primary": bool(any_changed),
        "line_verdict_changed_vs_primary": bool(line_changed),
        "any_verdict_changed_vs_primary": any_flagged,
    }, any_flagged


def _build_mandatory_disclosures_addb(
    target_line: str,
    pert_resolution_audit: Optional[Mapping[str, Any]],
) -> List[str]:
    disclosures = [
        "M6.1 TOST delta=0.04 hindsight caveat: report delta as fraction of line FCT r_orth_cal; INCONCLUSIVE means insufficient power, never equivalence ruled out.",
        "M6.2 Anti-masking rule: report H-I2 line-by-line; winning/tying FCT on weak-FCT or specificity-failing lines never masks strong-FCT inferiority.",
    ]
    if target_line == "a549_ko":
        disclosures.append(
            "CRISPRko modality shift vs Occam rhj_fewshot CRISPRi training; Tier-2 RBM14 untrained node (see 402-target sensitivity); 33/6983 RHJ genes absent (pairwise_common_axis mandatory)."
        )
    elif target_line == "cd4t_rest":
        disclosures.append(
            "413/6983 RHJ input genes absent (pairwise_common_axis mandatory); post-2026-02 release overlaps DiverseRNA-1.4M unsupervised pretraining cutoff window (7/9 blocks unattributed)."
        )
    elif target_line == "dld1":
        disclosures.append(
            "Low/borderline TOST power at N=75 (M4); pairwise_common_axis required for a549_ko/cd4t_rest/gse295214_h1 pairs."
        )
    elif target_line == "gse295214_h1":
        disclosures.extend([
            "shared H1 hESC lineage with Round 9 vcc2025_h1",
            "shared Replogle/Weissman laboratory provenance with Replogle-4",
            "5/6983 RHJ input genes absent (pairwise_common_axis mandatory); 2026-05-05 release unexcluded from DiverseRNA-1.4M unsupervised pretraining; Group 3 never pooled with Group 1/2.",
        ])
    elif target_line == "feng_ipsc":
        disclosures.append("SENSITIVITY_ONLY Q4 targeted 1,575-gene screen; never merged with Group 2 primary adjudication.")
    if pert_resolution_audit and int(pert_resolution_audit.get("n_tier2", 0)) > 0:
        disclosures.append(f"Tier-2 pert2id fallback targets: {pert_resolution_audit.get('tier2_targets')}")
    return disclosures


def evaluate_addb_line_complete(
    ref: CommittedR16ReferenceBundle,
    occam_matched: np.ndarray,
    occam_resample: np.ndarray,
    occam_mismatched: Mapping[str, np.ndarray],
    p4_common_pairs: Optional[Mapping[str, Any]] = None,
    occam_baseline_real_matched: Optional[np.ndarray] = None,
    occam_baseline_real_resample: Optional[np.ndarray] = None,
    occam_baseline_real_mismatched: Optional[Mapping[str, np.ndarray]] = None,
    occam_full44728_matched: Optional[np.ndarray] = None,
    occam_full44728_resample: Optional[np.ndarray] = None,
    occam_full44728_mismatched: Optional[Mapping[str, np.ndarray]] = None,
    pert_resolution_audit: Optional[Mapping[str, Any]] = None,
    output_gene_projection_audit: Optional[Mapping[str, Any]] = None,
    array_manifest_sha256: Optional[Mapping[str, str]] = None,
    n_boot: int = 4000,
    seed: int = DEFAULT_R16_SEED,
    subset_row_mask: Optional[np.ndarray] = None,
    include_tier1_sensitivity: bool = True,
) -> Dict[str, Any]:
    """
    Evaluate all primary & secondary Addendum B metrics (`§4.1`–`§4.3`) for one target line
    strictly following `ADDB_ARRAY_SCHEMA.md §7.2` on `(N_eval, G_eval)` matrices.
    """
    if subset_row_mask is not None:
        rm = np.asarray(subset_row_mask, dtype=bool)
        sub_idx = np.flatnonzero(rm)
        ref = CommittedR16ReferenceBundle(
            target_line=ref.target_line,
            role=ref.role,
            eval_perts=[ref.eval_perts[i] for i in sub_idx.tolist()],
            eval_genes=ref.eval_genes,
            target_gene_idx=ref.target_gene_idx[rm],
            delta_obs=ref.delta_obs[rm],
            delta_half_a=ref.delta_half_a[rm],
            delta_half_b=ref.delta_half_b[rm],
            half_valid=ref.half_valid[rm],
            delta_transfer=ref.delta_transfer[rm],
            source_grand_mean=ref.source_grand_mean,
            delta_fct=ref.delta_fct[rm],
            seen_sources=ref.seen_sources,
            unseen_sources=ref.unseen_sources,
            matched_barcodes=ref.matched_barcodes,
            resample_barcodes=ref.resample_barcodes,
            mismatched_barcodes=ref.mismatched_barcodes,
            se_delta_matched=ref.se_delta_matched[rm],
            se_delta_resample=ref.se_delta_resample[rm],
            se_delta_mismatched={k: v[rm] for k, v in ref.se_delta_mismatched.items()},
            mu_target=ref.mu_target,
            mu_sources=ref.mu_sources,
            e1_constmed=ref.e1_constmed[rm] if ref.e1_constmed is not None else None,
            n_ctrl_total=ref.n_ctrl_total,
            ctrl_sampling_fraction=ref.ctrl_sampling_fraction,
            r16_reference_sha256=ref.r16_reference_sha256,
            extra_fct_full_r9_6=None,
        )
        occam_matched = occam_matched[rm]
        occam_resample = occam_resample[rm]
        occam_mismatched = {k: v[rm] for k, v in occam_mismatched.items()}
        if occam_baseline_real_matched is not None:
            occam_baseline_real_matched = occam_baseline_real_matched[rm]
        if occam_baseline_real_resample is not None:
            occam_baseline_real_resample = occam_baseline_real_resample[rm]
        if occam_baseline_real_mismatched is not None:
            occam_baseline_real_mismatched = {k: v[rm] for k, v in occam_baseline_real_mismatched.items()}
        if occam_full44728_matched is not None:
            occam_full44728_matched = occam_full44728_matched[rm]
        if occam_full44728_resample is not None:
            occam_full44728_resample = occam_full44728_resample[rm]
        if occam_full44728_mismatched is not None:
            occam_full44728_mismatched = {k: v[rm] for k, v in occam_full44728_mismatched.items()}
        if p4_common_pairs is not None:
            new_p4: Dict[str, Any] = {}
            for s, val in p4_common_pairs.items():
                if isinstance(val, tuple):
                    pm, ps, nc = val
                    new_p4[s] = (pm[rm], ps[rm], nc)
                else:
                    new_val = dict(val)
                    new_val["matched"] = val["matched"][rm]
                    new_val["mismatched"] = val["mismatched"][rm]
                    new_p4[s] = new_val
            p4_common_pairs = new_p4

    d_obs = ref.delta_obs
    ha = ref.delta_half_a
    hb = ref.delta_half_b
    hv = ref.half_valid
    tm = ref.delta_transfer
    gm = ref.source_grand_mean
    tg = ref.target_gene_idx
    fct = ref.delta_fct

    # Primary rhj_6983_axis, baseline="reconstructed", r_orth_cal
    primary_block = _evaluate_arm_suite_addb(
        ref=ref, pred_matched=occam_matched, pred_resample=occam_resample,
        pred_mismatched=occam_mismatched, source_grand_mean=None,
        n_boot=n_boot, seed=seed,
        h_i1_key="H_I1_Occam", h_i2_key="H_I2_Occam", matched_key="Occam_matched_orth",
    )
    fct_orth = compute_orth_suite_r16(
        d_obs, fct, tm, tg, source_grand_mean=None,
        delta_half_a=ha, delta_half_b=hb, half_valid=hv,
        mode="all_nontarget", n_boot=n_boot, seed=seed,
    )
    primary_block["FCT_matched_orth"] = fct_orth

    # Secondary r_dorth_cal
    dorth_block = _evaluate_arm_suite_addb(
        ref=ref, pred_matched=occam_matched, pred_resample=occam_resample,
        pred_mismatched=occam_mismatched, source_grand_mean=gm,
        n_boot=n_boot, seed=seed,
        h_i1_key="H_I1_Occam_dorth", h_i2_key="H_I2_Occam_dorth", matched_key="Occam_matched_dorth",
    )

    # Committed STATE_SE reference calls
    se_orth_m = compute_orth_suite_r16(
        d_obs, ref.se_delta_matched, tm, tg, source_grand_mean=None,
        delta_half_a=ha, delta_half_b=hb, half_valid=hv,
        mode="all_nontarget", n_boot=n_boot, seed=seed,
    )
    se_paired_mis = compute_paired_model_vs_all_mismatched_r16(
        d_obs, ref.se_delta_matched, ref.se_delta_mismatched, tm, tg,
        source_grand_mean=None, delta_half_a=ha, delta_half_b=hb, half_valid=hv,
        mode="all_nontarget", n_boot=n_boot, seed=seed,
    )
    se_pos = bool(
        np.isfinite(se_orth_m.get("valid_frac", 0.0))
        and float(se_orth_m["valid_frac"]) >= MIN_BOOT_VALID_FRAC
        and np.isfinite(se_orth_m.get("rho_ab_orth", float("nan")))
        and float(se_orth_m["rho_ab_orth"]) > 0.0
        and np.isfinite(se_orth_m.get("ci975_lo", float("nan")))
        and float(se_orth_m["ci975_lo"]) > 0.0
    )
    se_res = compute_paired_diff_vs_resample_r16(
        d_obs, ref.se_delta_matched, ref.se_delta_resample, ref.se_delta_mismatched, tm, tg,
        source_grand_mean=None, delta_half_a=ha, delta_half_b=hb, half_valid=hv,
        mode="all_nontarget", n_boot=n_boot, seed=seed, matched_positive_ci975=se_pos,
    )
    se_h_i1 = adjudicate_claim_i1_line(
        se_orth_m, se_paired_mis, ref.seen_sources, ref.unseen_sources,
        delta=DEFAULT_TOST_DELTA, resample_comparison=se_res,
    )
    se_vs_fct = compute_paired_contrast_r16(
        d_obs, ref.se_delta_matched, fct, tm, tg, source_grand_mean=None,
        delta_half_a=ha, delta_half_b=hb, half_valid=hv,
        mode="all_nontarget", n_boot=n_boot, seed=seed,
    )
    se_h_i2 = adjudicate_claim_i2_line(se_vs_fct)

    # Mandatory P4 pairwise_common_axis
    p4_block: Optional[Dict[str, Any]] = None
    p4_any_changed = False
    if p4_common_pairs:
        p4_block, p4_any_changed = _evaluate_p4_common_axis_block(
            ref=ref, p4_common_pairs=p4_common_pairs,
            primary_h_i1=primary_block["H_I1_Occam"], n_boot=n_boot, seed=seed,
            primary_orth_matched=primary_block["Occam_matched_orth"],
        )

    # H-I3a & H-I3b
    hi3_orth = compute_hi3a_and_hi3b_line(
        ref, occam_matched, occam_resample, occam_mismatched,
        n_boot=n_boot, seed=seed, source_grand_mean=None,
    )
    hi3_dorth = compute_hi3a_and_hi3b_line(
        ref, occam_matched, occam_resample, occam_mismatched,
        n_boot=n_boot, seed=seed, source_grand_mean=gm,
    )
    adv_matrix = classify_advance_outcome_matrix(
        se_h_i1=se_h_i1,
        occam_h_i1=primary_block["H_I1_Occam"],
        occam_h_i2=primary_block["H_I2_Occam"],
        hi3a=hi3_orth["h_i3a"],
        hi3b=hi3_orth["h_i3b"],
        p4_any_verdict_change=p4_any_changed,
    )

    # Secondary baseline="real" and full_44728_api
    baseline_real_block: Optional[Dict[str, Any]] = None
    if occam_baseline_real_matched is not None:
        baseline_real_block = _evaluate_arm_suite_addb(
            ref=ref,
            pred_matched=occam_baseline_real_matched,
            pred_resample=occam_baseline_real_resample,
            pred_mismatched=occam_baseline_real_mismatched,
            source_grand_mean=None, n_boot=n_boot, seed=seed,
            h_i1_key="H_I1_Occam_real", h_i2_key="H_I2_Occam_real", matched_key="Occam_matched_orth",
        )

    full44728_block: Optional[Dict[str, Any]] = None
    if occam_full44728_matched is not None:
        full44728_block = _evaluate_arm_suite_addb(
            ref=ref,
            pred_matched=occam_full44728_matched,
            pred_resample=occam_full44728_resample,
            pred_mismatched=occam_full44728_mismatched,
            source_grand_mean=None, n_boot=n_boot, seed=seed,
            h_i1_key="H_I1_Occam_full44728", h_i2_key="H_I2_Occam_full44728", matched_key="Occam_matched_orth",
        )

    # Practical accuracy + M2 cross-mask shift + resp_top50
    practical_occam = compute_practical_metrics_r16(d_obs, occam_matched, tg, top_k=50)
    practical_fct = compute_practical_metrics_r16(d_obs, fct, tg, top_k=50)
    m2_shift = (
        compute_winner_curse_weight_shift_r16(ref.mu_target, d_obs, tg, top_k=50, n_boot=n_boot, seed=seed)
        if ref.mu_target is not None
        else None
    )
    resp50_orth = compute_orth_suite_r16(
        d_obs, occam_matched, tm, tg, source_grand_mean=None,
        delta_half_a=ha, delta_half_b=hb, half_valid=hv,
        mode="resp_top50", top_k=50, n_boot=n_boot, seed=seed,
    )
    practical_block: Dict[str, Any] = {
        **practical_occam,
        "m2_cross_mask_shift": m2_shift,
        "resp_top50_matched": resp50_orth,
        "practical_accuracy_FCT": practical_fct,
    }

    # Tier-1 trained-only sensitivity (`a549_ko` N=402 excluding Tier-2 RBM14)
    tier1_sens: Optional[Dict[str, Any]] = None
    if include_tier1_sensitivity:
        t2_set: Set[str] = set()
        if pert_resolution_audit and pert_resolution_audit.get("tier2_targets"):
            t2_set.update(str(x) for x in pert_resolution_audit["tier2_targets"])
        elif ref.target_line == "a549_ko" and "RBM14" in ref.eval_perts:
            t2_set.add("RBM14")
        if t2_set:
            t1_mask = np.asarray([p not in t2_set for p in ref.eval_perts], dtype=bool)
            if 0 < int(np.sum(t1_mask)) < len(ref.eval_perts):
                sub_res = evaluate_addb_line_complete(
                    ref=ref,
                    occam_matched=occam_matched,
                    occam_resample=occam_resample,
                    occam_mismatched=occam_mismatched,
                    p4_common_pairs=p4_common_pairs,
                    occam_baseline_real_matched=occam_baseline_real_matched,
                    occam_baseline_real_resample=occam_baseline_real_resample,
                    occam_baseline_real_mismatched=occam_baseline_real_mismatched,
                    occam_full44728_matched=occam_full44728_matched,
                    occam_full44728_resample=occam_full44728_resample,
                    occam_full44728_mismatched=occam_full44728_mismatched,
                    pert_resolution_audit={"n_tier1": int(np.sum(t1_mask)), "n_tier2": 0, "tier2_targets": []},
                    output_gene_projection_audit=output_gene_projection_audit,
                    array_manifest_sha256=array_manifest_sha256,
                    n_boot=n_boot,
                    seed=seed,
                    subset_row_mask=t1_mask,
                    include_tier1_sensitivity=False,
                )
                tier1_sens = {
                    "excluded_targets": sorted(t2_set),
                    "n_eval_targets": int(np.sum(t1_mask)),
                    "primary": sub_res["primary"],
                    "paired_vs_STATE_SE": sub_res["paired_vs_STATE_SE"],
                }

    # Group 3 extra FCT sensitivity (`preds_SENS_FCT_full_r9_6_incl_vcc2025_h1.npz` fixed to N=406)
    extra_fct_sens: Optional[Dict[str, Any]] = None
    if ref.extra_fct_full_r9_6 is not None:
        ef = ref.extra_fct_full_r9_6
        ef_fct = np.asarray(ef["delta_pred"], dtype=np.float64)
        ef_tm = np.asarray(ef["delta_transfer"], dtype=np.float64)
        ef_gm = np.asarray(ef["source_grand_mean"], dtype=np.float64)
        o_full_orth = compute_orth_suite_r16(
            d_obs, occam_matched, ef_tm, tg, source_grand_mean=None,
            delta_half_a=ha, delta_half_b=hb, half_valid=hv, mode="all_nontarget", n_boot=n_boot, seed=seed,
        )
        f_full_orth = compute_orth_suite_r16(
            d_obs, ef_fct, ef_tm, tg, source_grand_mean=None,
            delta_half_a=ha, delta_half_b=hb, half_valid=hv, mode="all_nontarget", n_boot=n_boot, seed=seed,
        )
        vs_full_orth = compute_paired_contrast_r16(
            d_obs, occam_matched, ef_fct, ef_tm, tg, source_grand_mean=None,
            delta_half_a=ha, delta_half_b=hb, half_valid=hv, mode="all_nontarget", n_boot=n_boot, seed=seed,
        )
        o_full_dorth = compute_orth_suite_r16(
            d_obs, occam_matched, ef_tm, tg, source_grand_mean=ef_gm,
            delta_half_a=ha, delta_half_b=hb, half_valid=hv, mode="all_nontarget", n_boot=n_boot, seed=seed,
        )
        f_full_dorth = compute_orth_suite_r16(
            d_obs, ef_fct, ef_tm, tg, source_grand_mean=ef_gm,
            delta_half_a=ha, delta_half_b=hb, half_valid=hv, mode="all_nontarget", n_boot=n_boot, seed=seed,
        )
        vs_full_dorth = compute_paired_contrast_r16(
            d_obs, occam_matched, ef_fct, ef_tm, tg, source_grand_mean=ef_gm,
            delta_half_a=ha, delta_half_b=hb, half_valid=hv, mode="all_nontarget", n_boot=n_boot, seed=seed,
        )
        extra_fct_sens = {
            "source_lines_fct": list(ef.get("source_lines_fct", R9_6_ALL_LINES)),
            "n_eval": len(ref.eval_perts),
            "Occam_orth": o_full_orth,
            "FCT_orth": f_full_orth,
            "Occam_vs_FCT_orth": vs_full_orth,
            "H_I2_Occam": adjudicate_claim_i2_line(vs_full_orth),
            "Occam_dorth": o_full_dorth,
            "FCT_dorth": f_full_dorth,
            "Occam_vs_FCT_dorth": vs_full_dorth,
            "H_I2_Occam_dorth": adjudicate_claim_i2_line(vs_full_dorth),
        }

    gmeta = ADDB_LINE_GROUP_META.get(
        ref.target_line, {"group_key": ref.role, "role": ref.role, "modality": "CRISPRi"}
    )
    n_self_masked = int(np.sum(np.asarray(tg) >= 0))
    n_half_valid = int(np.sum(np.asarray(hv, dtype=bool)))
    fct_est = float(fct_orth.get("estimate", float("nan")))
    m61_rel = (
        float(DEFAULT_TOST_DELTA / fct_est)
        if (np.isfinite(fct_est) and abs(fct_est) > 1e-12)
        else float("nan")
    )

    return {
        "line_name": ref.target_line,
        "group_key": gmeta["group_key"],
        "role": gmeta["role"],
        "modality": gmeta["modality"],
        "checkpoint_id": "rhj_fewshot/all_sources/seed_42",
        "n_eval_targets": len(ref.eval_perts),
        "n_eval_genes": len(ref.eval_genes),
        "n_half_valid": n_half_valid,
        "n_self_masked": n_self_masked,
        "m61_rel_delta_share": m61_rel,
        "strata_sources": {
            "seen_by_checkpoint": list(ref.seen_sources),
            "unseen_by_checkpoint": list(ref.unseen_sources),
        },
        "pert_resolution_audit": dict(pert_resolution_audit) if pert_resolution_audit is not None else {},
        "output_gene_projection_audit": dict(output_gene_projection_audit) if output_gene_projection_audit is not None else {},
        "mandatory_disclosures": _build_mandatory_disclosures_addb(ref.target_line, pert_resolution_audit),
        "array_manifest_sha256": dict(array_manifest_sha256) if array_manifest_sha256 is not None else {},
        "primary": primary_block,
        "paired_vs_STATE_SE": {
            "H_I3a": {"r_orth_cal": hi3_orth["h_i3a"], "r_dorth_cal": hi3_dorth["h_i3a"]},
            "H_I3b": {"r_orth_cal": hi3_orth["h_i3b"], "r_dorth_cal": hi3_dorth["h_i3b"]},
            "STATE_SE_reference_calls": {
                "H_I1_line_summary": se_h_i1["line_summary"],
                "H_I2_classification": se_h_i2["classification"],
                "matched_r_orth_cal": se_orth_m,
            },
            "advance_interpretation_case": adv_matrix,
        },
        "secondary": {
            "r_dorth_cal": dorth_block,
            "pairwise_common_axis": p4_block,
            "baseline_real": baseline_real_block,
            "full_44728_api": full44728_block,
            "practical_accuracy_Occam": practical_block,
            "tier1_trained_only_sensitivity": tier1_sens,
            "extra_fct_sensitivity": extra_fct_sens,
        },
    }


def compute_group_adjudications_addb(
    line_entries: Mapping[str, Mapping[str, Any]],
) -> Dict[str, Any]:
    """
    Aggregate strictly separated Addendum B group verdicts (`ADDB_ARRAY_SCHEMA.md §7.1`).
    Never pools CRISPRko vs CRISPRi, Group 2 vs Group 3, or Sensitivity lines.
    """
    specs: Tuple[Tuple[str, Tuple[str, ...], str, List[str]], ...] = (
        ("group2_crisprko", ("a549_ko",), "RETROSPECTIVE_ADJUDICATED", ["CRISPRko modality shift; Tier-2 RBM14 sensitivity."]),
        ("group2_crispri", ("dld1", "cd4t_rest"), "RETROSPECTIVE_ADJUDICATED", ["M4 low-power caveat on dld1 (N=75); P4 413-gene asymmetry on cd4t_rest."]),
        ("group3_crispri", ("gse295214_h1",), "PROSPECTIVE_HOLDOUT_SEPARATE", [
            "shared H1 hESC lineage with Round 9 vcc2025_h1",
            "shared Replogle/Weissman laboratory provenance with Replogle-4",
            "strictly separated from Group 1/2",
        ]),
        ("sensitivity_only", ("feng_ipsc",), "SENSITIVITY_ONLY", ["Targeted 1,575-gene screen; sensitivity only."]),
    )
    out: Dict[str, Any] = {}
    for gkey, glines, gstatus, gdisc in specs:
        present = [ln for ln in glines if ln in line_entries]
        h1_counts: Dict[str, int] = {}
        h2_counts: Dict[str, int] = {}
        h1_per: Dict[str, str] = {}
        h2_per: Dict[str, str] = {}
        k_spec = 0
        for ln in present:
            le = line_entries[ln]
            v1 = str(le["primary"]["H_I1_Occam"]["line_summary"])
            v2 = str(le["primary"]["H_I2_Occam"]["classification"])
            h1_per[ln] = v1
            h2_per[ln] = v2
            h1_counts[v1] = h1_counts.get(v1, 0) + 1
            h2_counts[v2] = h2_counts.get(v2, 0) + 1
            if v1 in ("LINE_MATCHED_SPECIFIC", "LINE_MATCHED_SPECIFIC_SUB_DELTA"):
                k_spec += 1
        k_len = len(present)
        out[gkey] = {
            "K": int(k_len),
            "lines": present,
            "adjudication_status": gstatus if k_len > 0 else "NOT_EVALUATED",
            "H_I1_Occam_verdict_counts": h1_counts,
            "H_I1_Occam_k_over_K": f"{k_spec}/{k_len}",
            "H_I1_Occam_per_line": h1_per,
            "H_I2_Occam_verdict_counts": h2_counts,
            "H_I2_Occam_per_line": h2_per,
            "mandatory_disclosures": gdisc,
        }
    return out
