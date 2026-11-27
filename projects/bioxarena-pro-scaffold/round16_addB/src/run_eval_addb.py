#!/usr/bin/env python3
"""
Round 16 Addendum B (`Occam`) Pre-Registered Evaluation Engine (`PREREG_I_ADDENDUM_B_DRAFT_v2.md`).

Implements the complete locked Addendum B v2 protocol (Lock L1/L2, DEVIATION_R16_ADDB_D1,
Gate 5A Item 3 `CPU float32` device lock, and `ADDB_ARRAY_SCHEMA.md`):
  1. Hard Execution Gate (`assert_addb_execution_gates`): Refuses ANY real-target execution
     unless BOTH `round16_addB/R16_ADDB_CODE_FREEZE.txt` AND `round16_addB/R16_ADDB_GO.txt`
     (plus `PREREG_I_ADDENDUM_B_LOCK.md` and `R16_ADDB_HUMAN_APPROVAL.txt`) exist and are
     non-empty. Strictly locks `device="cpu"` (`CPU float32`, `torch.set_num_threads(4)`).
  2. Target Lines & `§2.3` Strata (Group 1 `nadig_hepg2`/`nadig_jurkat` strictly refused as
     targets and retained as `seen_by_checkpoint` mismatched sources):
       - Group 2 Primary CRISPRko: `a549_ko` (`N=403`, `G=1837`, `RBM14` Tier-2 `402` sensitivity)
       - Group 2 Primary CRISPRi:  `dld1` (`N=75`, `G=1832`), `cd4t_rest` (`N=245`, `G=1830`)
       - Group 3 Prospective:      `gse295214_h1` (`N=406`, `G=1837`, `FCT_full_r9_6` sensitivity)
       - Sensitivity Only (`Q4`):  `feng_ipsc` (`N=88`, `G=1575`)
  3. Scoring Axes & Baselines (`§3.3–§3.4`, `ADDB_ARRAY_SCHEMA.md`):
       - Primary `rhj_6983_axis` (`baseline="reconstructed"` primary + `baseline="real"` sensitivity)
       - Mandatory `P4` `pairwise_common_axis` (`baseline="reconstructed"`, every pair with `n_union_absent > 0`)
       - Secondary `full_44728_api` (`baseline="reconstructed"`)
  4. Metrics, Cascades, & Contrasts (`§4.1–§4.3`):
       - Exact Round 16 `r_orth_cal` / `r_dorth_cal`, `H-I1_Occam`, `H-I2_Occam`, within-line
         disjoint resample reference, practical accuracy, and `M2` cross-mask shift imported
         read-only from frozen `round16/src/{metrics_r16,adjudicate_r16}.py`.
       - Pre-registered paired `STATE_SE` vs. `Occam` contrasts (`H-I3a` matched accuracy,
         `H-I3b` per-source & per-stratum-mean control specificity + resample-excess).
  5. `--synthetic-smoke`: Pure-CPU, ZERO-real-cell end-to-end test + `z_ctrl`/`arrays` export +
     Stage 4 from-disk recomputation parity check (`max_abs_diff <= 1e-12`).
"""

from __future__ import annotations

import argparse
import gc
import json
import os
import shutil
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Set, Tuple

import numpy as np

ADDB_SRC_DIR = Path(__file__).resolve().parent
ADDB_ROOT = ADDB_SRC_DIR.parent
REPO_ROOT = ADDB_ROOT.parent
R16_SRC_DIR = REPO_ROOT / "round16" / "src"
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
for _p in (str(ADDB_SRC_DIR), str(R16_SRC_DIR), str(R10_SRC_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from adjudicate_r16 import (  # noqa: E402
    MIN_BOOT_VALID_FRAC,
    SUB_DELTA_INTERPRETATION_NOTE,
    adjudicate_claim_i1_line,
    adjudicate_claim_i2_line,
    adjudicate_group_r16,
    classify_pair_tost_r16,
    summarize_within_line_resample_r16,
)
from data_loader_r10 import DEFAULT_HGNC_FILE, HGNCBridge  # noqa: E402
from export_addb import (  # noqa: E402
    OCCAM_CHECKPOINT_ID,
    OCCAM_LATENT_DIM,
    OCCAM_MODEL_KEY,
    build_occam_delta_payload,
    save_occam_line_deltas,
    save_occam_z_ctrl_latents,
    sha256_file,
)
from metrics_r16 import (  # noqa: E402
    DEFAULT_R16_SEED,
    DEFAULT_RESAMPLE_SEED,
    DEFAULT_TOST_DELTA,
    _array_max_abs_diff_and_exact_inv,
    _collect_mode_moments,
    _eval_boot_cal,
    _extract_masked_idx_per_pert,
    compute_joint_cross_mask_r16,
    compute_orth_suite_r16,
    compute_paired_contrast_r16,
    compute_paired_model_vs_all_mismatched_r16,
    compute_paired_with_within_line_resample_r16,
    summarize_boot_r16,
)
from occam_adapter_addb import (  # noqa: E402
    DEFAULT_OCCAM_CHECKSUMS_FILE,
    DEFAULT_OCCAM_CKPT_PATH,
    DEFAULT_OCCAM_GENE_AXIS_CSV,
    DEFAULT_OCCAM_PERTURBATIONS_CSV,
    DEFAULT_OCCAM_SRC_DIR,
    LOCKED_BEST_MODEL_SHA256,
    OccamInferenceAdapter,
    configure_occam_runtime,
)
from preprocess_addb import (  # noqa: E402
    build_output_gene_projection,
    build_pairwise_common_axis_mask,
    compute_string_list_sha256,
    load_rhj_gene_axis,
    load_rhj_perturbations_table,
    resolve_perturbation_keys_two_tier,
    scatter_full_44728_api,
    scatter_rhj_6983_to_encoder,
)
from raw_controls_addb import (  # noqa: E402
    DEFAULT_LOCAL_NADIG_DIR,
    DEFAULT_PRIMARY_BARCODES_JSON,
    DEFAULT_SENSITIVITY_BARCODES_JSON,
    RawControlPool,
    load_10x_multi_h5_controls_by_barcodes,
    load_raw_controls_by_locked_barcodes,
)

DEFAULT_GROUP3_BARCODES_JSON = REPO_ROOT / "round16" / "results" / "group3" / "control_barcodes_r16.json"
DEFAULT_GSE295214_RAW_DIR = Path("/Users/mac/.cache/r16_gse295214/raw")
DEFAULT_ADDB_RESULTS_DIR = ADDB_ROOT / "results"

DEFAULT_CODE_FREEZE_PATH = ADDB_ROOT / "R16_ADDB_CODE_FREEZE.txt"
DEFAULT_GO_PATH = ADDB_ROOT / "R16_ADDB_GO.txt"
DEFAULT_LOCK_PATH = ADDB_ROOT / "PREREG_I_ADDENDUM_B_LOCK.md"
DEFAULT_HUMAN_APPROVAL_PATH = ADDB_ROOT / "R16_ADDB_HUMAN_APPROVAL.txt"

LOCKED_DEVICE = "cpu"
LOCKED_CPU_NUM_THREADS = 4
LOCKED_BATCH_SIZE = 256
LOCKED_N_BOOT = 4000

EXCLUDED_GROUP1_TRAINING_LINES: Tuple[str, ...] = ("nadig_hepg2", "nadig_jurkat")

ALL_EVALUATED_LINES_ADDB: Tuple[str, ...] = (
    "a549_ko",
    "dld1",
    "cd4t_rest",
    "gse295214_h1",
    "feng_ipsc",
)
GROUP2_LINES_ADDB: Tuple[str, ...] = ("a549_ko", "dld1", "cd4t_rest")
GROUP3_LINES_ADDB: Tuple[str, ...] = ("gse295214_h1",)
SENSITIVITY_LINES_ADDB: Tuple[str, ...] = ("feng_ipsc",)

OBSERVED_NPZ_REQUIRED_KEYS: Tuple[str, ...] = (
    "delta_obs",
    "delta_half_a",
    "delta_half_b",
    "half_valid",
    "delta_transfer",
    "source_grand_mean",
    "fct_delta_pred",
    "e1_constmed_delta_pred",
    "target_indices",
    "target_names",
    "gene_names",
    "mu_target",
    "mu_sources",
    "matched_ctrl_barcodes",
    "resample_ctrl_barcodes",
    "n_ctrl_total",
    "n_ctrl_sampled",
    "ctrl_sampling_fraction",
)

ADDB_LINE_SPECS: Dict[str, Dict[str, Any]] = {
    "a549_ko": {
        "group_key": "group2_crisprko",
        "r16_group_key": "group2",
        "role": "primary",
        "modality": "CRISPRko",
        "seen_by_checkpoint": ("nadig_hepg2", "nadig_jurkat"),
        "unseen_by_checkpoint": ("cd4t_rest", "dld1", "feng_ipsc"),
        "r16_arrays_rel": "round16/results/arrays/a549_ko",
        "disclosures": (
            "CRISPRko-vs-CRISPRi modality shift vs rhj_fewshot training (I11d).",
            "33/6,983 RHJ input genes absent in raw A549 (P4 pairwise_common_axis mandatory).",
            "1/403 target (RBM14) resolved via Tier-2 graph extrapolation (402-target Tier-1 sensitivity reported).",
            "Strictly post-dates DiverseRNA (arXiv:2605.31562).",
        ),
    },
    "dld1": {
        "group_key": "group2_crispri",
        "r16_group_key": "group2",
        "role": "primary",
        "modality": "CRISPRi",
        "seen_by_checkpoint": ("nadig_hepg2", "nadig_jurkat"),
        "unseen_by_checkpoint": ("a549_ko", "cd4t_rest", "feng_ipsc"),
        "r16_arrays_rel": "round16/results/arrays/dld1",
        "disclosures": (
            "Low/borderline TOST power at N=75 (M4; INCONCLUSIVE reported strictly as insufficient power).",
            "P4 pairwise_common_axis reported for pairs involving a549_ko or cd4t_rest.",
            "Strictly post-dates DiverseRNA (arXiv:2605.31562).",
        ),
    },
    "cd4t_rest": {
        "group_key": "group2_crispri",
        "r16_group_key": "group2",
        "role": "primary",
        "modality": "CRISPRi",
        "seen_by_checkpoint": ("nadig_hepg2", "nadig_jurkat"),
        "unseen_by_checkpoint": ("a549_ko", "dld1", "feng_ipsc"),
        "r16_arrays_rel": "round16/results/arrays/cd4t_rest",
        "disclosures": (
            "413/6,983 RHJ input genes absent in raw cd4t_rest (P4 pairwise_common_axis mandatory for all pairs).",
            "GEO 2026-02-06 precedes DiverseRNA May 2026 (7/9 encoder pretraining blocks unattributed; supervised perturbation unseen).",
            "Fragile relative TOST delta=0.04 vs weak FCT (M6.1/M6.2 anti-masking rule applies).",
        ),
    },
    "gse295214_h1": {
        "group_key": "group3_crispri",
        "r16_group_key": "group3",
        "role": "prospective_group3",
        "modality": "CRISPRi",
        "seen_by_checkpoint": ("nadig_hepg2", "nadig_jurkat"),
        "unseen_by_checkpoint": ("a549_ko", "cd4t_rest", "dld1", "feng_ipsc"),
        "r16_arrays_rel": "round16/results/group3/arrays/gse295214_h1",
        "disclosures": (
            "Reported strictly separately from Group 2 (never pooled).",
            "Shared H1 hESC cell lineage with Round 9 vcc2025_h1 (excluded from primary FCT pool; full r9_6 FCT sensitivity reported).",
            "Shared Replogle/Weissman laboratory provenance with Replogle-4 training lines.",
            "GEO 2026-05-05 overlaps DiverseRNA May 2026 (encoder pretraining absence unverifiable; supervised perturbation unseen).",
            "5/6,983 RHJ input genes absent in raw GSE295214 (P4 pairwise_common_axis mandatory for all pairs).",
        ),
    },
    "feng_ipsc": {
        "group_key": "sensitivity_only",
        "r16_group_key": "sensitivity_only",
        "role": "sensitivity_only",
        "modality": "CRISPRi",
        "seen_by_checkpoint": ("nadig_hepg2", "nadig_jurkat"),
        "unseen_by_checkpoint": ("a549_ko", "cd4t_rest", "dld1"),
        "r16_arrays_rel": "round16/results/sensitivity/arrays/feng_ipsc",
        "disclosures": (
            "SENSITIVITY_ONLY target (Q4; ENA 2024-11-08 predates both SE-600M and DiverseRNA); never pooled with Group 2.",
            "Retained as unseen_by_checkpoint mismatched source for Group 2 and Group 3.",
            "P4 pairwise_common_axis reported for pairs involving a549_ko or cd4t_rest.",
        ),
    },
}


class AddBExecutionGateError(RuntimeError):
    """Raised when Addendum B real-target execution is attempted without required freeze/GO gates."""


def _now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _cpu_time_s() -> float:
    t = os.times()
    return float(t.user + t.system)


def assert_addb_execution_gates(
    *,
    code_freeze_path: Path = DEFAULT_CODE_FREEZE_PATH,
    go_path: Path = DEFAULT_GO_PATH,
    lock_path: Path = DEFAULT_LOCK_PATH,
    human_approval_path: Path = DEFAULT_HUMAN_APPROVAL_PATH,
) -> Dict[str, str]:
    """
    Hard pre-execution gate refusing ANY real-target Occam inference unless ALL required
    governance files exist and are non-empty (`R16_ADDB_CODE_FREEZE.txt`, `R16_ADDB_GO.txt`,
    `PREREG_I_ADDENDUM_B_LOCK.md`, `R16_ADDB_HUMAN_APPROVAL.txt`).
    """
    required = {
        "code_freeze_path": Path(code_freeze_path),
        "go_path": Path(go_path),
        "lock_path": Path(lock_path),
        "human_approval_path": Path(human_approval_path),
    }
    missing: List[str] = []
    out: Dict[str, str] = {}
    for label, fpath in required.items():
        if not fpath.is_file() or fpath.stat().st_size <= 0:
            missing.append(f"{label}={fpath}")
        else:
            out[label] = str(fpath)
    if missing:
        raise AddBExecutionGateError(
            "Refusing Addendum B real-target execution — required gate file(s) missing or empty: "
            + ", ".join(missing)
        )
    return out


def validate_target_lines_addb(lines: Sequence[str]) -> List[str]:
    """Validate requested target lines and strictly refuse Group 1 training lines."""
    if not lines:
        raise ValueError("At least one Addendum B target line must be specified.")
    out: List[str] = []
    seen: Set[str] = set()
    for raw_ln in lines:
        ln = str(raw_ln).strip()
        if ln in EXCLUDED_GROUP1_TRAINING_LINES:
            raise ValueError(
                f"Group 1 training line {ln!r} is strictly EXCLUDED as an Occam evaluation target "
                f"(PREREG_I_ADDENDUM_B_DRAFT_v2.md §2.1; used only as seen_by_checkpoint mismatched pool)."
            )
        if ln not in ADDB_LINE_SPECS:
            raise ValueError(
                f"Unknown Addendum B target line {ln!r}; allowed: {list(ALL_EVALUATED_LINES_ADDB)}."
            )
        if ln not in seen:
            seen.add(ln)
            out.append(ln)
    return out


def build_gse295214_prefix_map(raw_dir: Path = DEFAULT_GSE295214_RAW_DIR) -> Tuple[Dict[str, Path], str]:
    """
    Build `{GSM_id: path_to_h5}` prefix map for `gse295214_h1` (`85` 10x `.h5` files)
    matching Round 16's `GSM89437XX:<barcode>` locked control barcode format.
    """
    base = Path(raw_dir)
    if not base.is_dir():
        raise FileNotFoundError(f"GSE295214 raw directory not found: {base}")
    h5_files = sorted(base.glob("*.h5"))
    if not h5_files:
        raise FileNotFoundError(f"No .h5 files found under {base}")
    pmap: Dict[str, Path] = {}
    for fp in h5_files:
        pref = fp.name.split("_", 1)[0]
        if pref in pmap:
            raise ValueError(f"Duplicate GSE295214 GSM prefix {pref!r}: {pmap[pref]} vs {fp}")
        pmap[pref] = fp
    ref_pref = sorted(pmap.keys())[0]
    return pmap, ref_pref


def load_raw_control_pool_addb(
    pool_name: str,
    barcodes: Sequence[str],
    *,
    pool_kind: str = "primary",
    seed: int = DEFAULT_R16_SEED,
    repo_root: Path = REPO_ROOT,
    raw_nadig_dir: Optional[Path] = None,
    gse295214_raw_dir: Optional[Path] = None,
) -> RawControlPool:
    """Unified read-only raw control loader covering all 6 local pools + `gse295214_h1`."""
    if pool_name == "gse295214_h1":
        raw_d = Path(gse295214_raw_dir) if gse295214_raw_dir is not None else DEFAULT_GSE295214_RAW_DIR
        pmap, ref_pref = build_gse295214_prefix_map(raw_d)
        return load_10x_multi_h5_controls_by_barcodes(
            pmap,
            barcodes,
            line_name="gse295214_h1",
            canonical_ref_prefix=ref_pref,
            pool_kind=pool_kind,
            seed=seed,
        )
    return load_raw_controls_by_locked_barcodes(
        pool_name,
        barcodes,
        pool_kind=pool_kind,
        seed=seed,
        repo_root=repo_root,
        raw_nadig_dir=raw_nadig_dir,
        verify_nadig_ntc=True,
    )


def load_consolidated_r16_barcodes_addb(
    target_lines: Sequence[str],
    *,
    primary_json: Path = DEFAULT_PRIMARY_BARCODES_JSON,
    sensitivity_json: Path = DEFAULT_SENSITIVITY_BARCODES_JSON,
    group3_json: Path = DEFAULT_GROUP3_BARCODES_JSON,
) -> Dict[str, Any]:
    """
    Load and cross-verify exact locked Round 16 control barcodes (`seed=20261601` primary &
    `seed=20261602` disjoint resample) across all needed pools with ZERO new control sampling.
    """
    needed_primary: Set[str] = set()
    needed_resample: Set[str] = set()
    for ln in target_lines:
        spec = ADDB_LINE_SPECS[ln]
        needed_primary.add(ln)
        needed_resample.add(ln)
        needed_primary.update(spec["seen_by_checkpoint"])
        needed_primary.update(spec["unseen_by_checkpoint"])

    src_paths: List[Path] = [Path(primary_json)]
    if "feng_ipsc" in needed_resample:
        src_paths.append(Path(sensitivity_json))
    if "gse295214_h1" in needed_primary or "gse295214_h1" in needed_resample:
        src_paths.append(Path(group3_json))

    loaded_sources: Dict[str, Dict[str, Any]] = {}
    sources_sha: Dict[str, str] = {}
    for sp in src_paths:
        if not sp.is_file():
            raise FileNotFoundError(f"Required locked Round 16 barcode file not found: {sp}")
        rel = str(sp.relative_to(REPO_ROOT)) if str(sp.resolve()).startswith(str(REPO_ROOT.resolve())) else str(sp)
        sources_sha[rel] = sha256_file(sp)
        loaded_sources[rel] = json.loads(sp.read_text(encoding="utf-8"))

    pools_out: Dict[str, Dict[str, Any]] = {}
    for pool_name in sorted(needed_primary):
        primary_bcs: Optional[List[str]] = None
        n_ctrl_total = 0
        resamp_entry: Optional[Dict[str, Any]] = None

        for rel, blob in loaded_sources.items():
            if pool_name not in blob:
                continue
            entry = blob[pool_name]
            bcs = [str(b) for b in entry["barcodes"]]
            if primary_bcs is None:
                primary_bcs = bcs
                n_ctrl_total = int(entry.get("n_ctrl_total", len(bcs)))
            elif bcs != primary_bcs:
                raise ValueError(
                    f"Locked Round 16 primary barcode mismatch for {pool_name!r} across R16 JSON files ({rel})."
                )

            wr = entry.get("within_line_resample") or {}
            if wr.get("eligible") and len(wr.get("barcodes") or []) > 0:
                rbcs = [str(b) for b in wr["barcodes"]]
                if resamp_entry is None:
                    resamp_entry = {
                        "eligible": True,
                        "seed": int(wr.get("seed", DEFAULT_RESAMPLE_SEED)),
                        "n_resample": len(rbcs),
                        "overlap_with_primary": int(len(set(bcs) & set(rbcs))),
                        "regime": str(wr.get("regime", "full_512_disjoint")),
                        "barcodes": rbcs,
                        "barcodes_sha256": compute_string_list_sha256(rbcs),
                    }
                    if resamp_entry["overlap_with_primary"] != 0:
                        raise ValueError(
                            f"Within-line resample for {pool_name!r} has non-zero overlap "
                            f"{resamp_entry['overlap_with_primary']} with primary barcodes."
                        )
                elif rbcs != resamp_entry["barcodes"]:
                    raise ValueError(
                        f"Locked Round 16 resample barcode mismatch for {pool_name!r} across R16 JSON files ({rel})."
                    )

        if primary_bcs is None or len(primary_bcs) == 0:
            raise KeyError(f"Primary locked barcodes for pool {pool_name!r} missing from all R16 barcode JSONs.")
        if pool_name in needed_resample and resamp_entry is None:
            raise KeyError(f"Within-line resample locked barcodes for target {pool_name!r} missing from R16 JSONs.")

        pools_out[pool_name] = {
            "primary_seed": int(DEFAULT_R16_SEED),
            "n_ctrl_total": int(n_ctrl_total),
            "n_ctrl_sampled": int(len(primary_bcs)),
            "barcodes": primary_bcs,
            "barcodes_sha256": compute_string_list_sha256(primary_bcs),
            "within_line_resample": resamp_entry
            if resamp_entry is not None
            else {
                "eligible": False,
                "seed": int(DEFAULT_RESAMPLE_SEED),
                "n_resample": 0,
                "overlap_with_primary": 0,
                "regime": "",
                "barcodes": [],
                "barcodes_sha256": "",
            },
        }

    return {
        "provenance": {
            "primary_seed": int(DEFAULT_R16_SEED),
            "resample_seed": int(DEFAULT_RESAMPLE_SEED),
            "r16_sources_sha256": sources_sha,
        },
        "pools": pools_out,
    }


# ==============================================================================
# Statistical Helpers: Practical Accuracy, H-I3a, H-I3b, Pairwise Common Axis
# ==============================================================================


def compute_practical_accuracy_metrics_addb(
    delta_obs: np.ndarray,
    delta_pred: np.ndarray,
    target_indices: Sequence[int],
    top_k: int = 50,
) -> Dict[str, Any]:
    """Descriptive practical accuracy metrics (`raw_pearson`, `centered_de50_pearson`, `cosine`) matching R16."""
    obs = np.asarray(delta_obs, dtype=np.float64)
    pred = np.asarray(delta_pred, dtype=np.float64)
    n_p, n_g = obs.shape
    raw_rs: List[float] = []
    de50_rs: List[float] = []
    cos_vals: List[float] = []
    for i in range(n_p):
        mask = np.ones(n_g, dtype=bool)
        ti = int(target_indices[i])
        if 0 <= ti < n_g:
            mask[ti] = False
        x = obs[i, mask]
        y = pred[i, mask]
        sx, sy = float(np.std(x)), float(np.std(y))
        if sx > 1e-12 and sy > 1e-12:
            raw_rs.append(float(np.corrcoef(x, y)[0, 1]))
        nx, ny = float(np.linalg.norm(x)), float(np.linalg.norm(y))
        if nx > 1e-12 and ny > 1e-12:
            cos_vals.append(float(np.dot(x, y) / (nx * ny)))
        k = min(int(top_k), x.size)
        if k >= 3:
            top_idx = np.argpartition(np.abs(x), -k)[-k:]
            xt, yt = x[top_idx], y[top_idx]
            sxt, syt = float(np.std(xt)), float(np.std(yt))
            if sxt > 1e-12 and syt > 1e-12:
                de50_rs.append(float(np.corrcoef(xt, yt)[0, 1]))
    return {
        "mean_raw_pearson": float(np.mean(raw_rs)) if raw_rs else float("nan"),
        "mean_centered_de50_pearson": float(np.mean(de50_rs)) if de50_rs else float("nan"),
        "mean_cosine_similarity": float(np.mean(cos_vals)) if cos_vals else float("nan"),
        "n_perts": int(n_p),
    }


def _classify_two_sided_975(
    summary_dict: Mapping[str, Any],
    *,
    pos_label: str,
    neg_label: str,
    indist_label: str = "INDISTINGUISHABLE",
) -> str:
    """Symmetric 97.5% CI classification with strict zero-tie & valid_frac >= 0.95 guards (`§4.2`)."""
    vf = float(summary_dict.get("valid_frac", 0.0))
    lo = float(summary_dict.get("ci975_lo", float("nan")))
    hi = float(summary_dict.get("ci975_hi", float("nan")))
    if vf < MIN_BOOT_VALID_FRAC or not (np.isfinite(lo) and np.isfinite(hi)):
        return indist_label
    if lo == 0.0 or hi == 0.0:
        return indist_label
    if lo > 0.0:
        return pos_label
    if hi < 0.0:
        return neg_label
    return indist_label


def compute_h_i3a_contrast_addb(
    delta_obs: np.ndarray,
    occam_matched: np.ndarray,
    state_se_matched: np.ndarray,
    delta_transfer: np.ndarray,
    target_indices: Sequence[int],
    source_grand_mean: Optional[np.ndarray],
    delta_half_a: np.ndarray,
    delta_half_b: np.ndarray,
    half_valid: np.ndarray,
    *,
    n_boot: int = LOCKED_N_BOOT,
    seed: int = DEFAULT_R16_SEED,
) -> Dict[str, Any]:
    """
    `PREREG_I_ADDENDUM_B_DRAFT_v2.md §4.2` (`H-I3a`):
    Matched-control accuracy contrast `Delta_acc = r(Occam_matched) - r(STATE_SE_matched)`
    on both primary `r_orth_cal` and descriptive `r_dorth_cal`.
    """
    out: Dict[str, Any] = {}
    for metric_key, b_s in (("r_orth_cal", None), ("r_dorth_cal", source_grand_mean)):
        c = compute_paired_contrast_r16(
            delta_obs=delta_obs,
            delta_pred_a=occam_matched,
            delta_pred_b=state_se_matched,
            delta_transfer=delta_transfer,
            target_gene_idx=target_indices,
            source_grand_mean=b_s,
            delta_half_a=delta_half_a,
            delta_half_b=delta_half_b,
            half_valid=half_valid,
            mode="all_nontarget",
            n_boot=n_boot,
            seed=seed,
        )
        cls = _classify_two_sided_975(
            c,
            pos_label="OCCAM_SUPERIOR_TO_STATE_SE",
            neg_label="STATE_SE_SUPERIOR_TO_OCCAM",
        )
        out[metric_key] = {
            "classification": cls,
            "estimate": float(c["estimate"]),
            "ci975_lo": float(c["ci975_lo"]),
            "ci975_hi": float(c["ci975_hi"]),
            "ci95_lo": float(c["ci95_lo"]),
            "ci95_hi": float(c["ci95_hi"]),
            "boot_se": float(c["boot_se"]),
            "valid_frac": float(c["valid_frac"]),
            "n_items": int(c["n_items"]),
            "n_boot": int(c["n_boot"]),
            "rho_ab_orth": float(c["rho_ab_orth"]),
            "ceiling": float(c["ceiling"]),
            "double_orth": bool(c["double_orth"]),
            "occam_summary": c["a_summary"],
            "state_se_summary": c["b_summary"],
        }
    return out


def compute_h_i3b_contrast_addb(
    delta_obs: np.ndarray,
    occam_matched: np.ndarray,
    occam_resample: np.ndarray,
    occam_mismatched: Mapping[str, np.ndarray],
    state_se_matched: np.ndarray,
    state_se_resample: np.ndarray,
    state_se_mismatched: Mapping[str, np.ndarray],
    seen_sources: Sequence[str],
    unseen_sources: Sequence[str],
    delta_transfer: np.ndarray,
    target_indices: Sequence[int],
    source_grand_mean: Optional[np.ndarray],
    delta_half_a: np.ndarray,
    delta_half_b: np.ndarray,
    half_valid: np.ndarray,
    *,
    n_boot: int = LOCKED_N_BOOT,
    seed: int = DEFAULT_R16_SEED,
    occam_matched_by_source: Optional[Mapping[str, np.ndarray]] = None,
) -> Dict[str, Any]:
    """
    `PREREG_I_ADDENDUM_B_DRAFT_v2.md §4.2` (`H-I3b`):
    Paired control-specificity contrast sharing the exact same joint perturbation bootstrap draws:
      - Per mismatched source `s'`:
          `Delta_spec(s') = Delta_pair_Occam(s') - Delta_pair_STATE_SE(s')`
          `Delta_spec_excess(s') = Delta_excess_Occam(s') - Delta_excess_STATE_SE(s')`
      - Per stratum mean (`seen_by_checkpoint`, `unseen_by_checkpoint` separately):
          Mean over stratum sources computed INSIDE each bootstrap draw for exact pairing.

    When `occam_matched_by_source` is supplied (for `pairwise_common_axis` cross-check),
    `Delta_pair_Occam(s') = r(occam_matched_by_source[s']) - r(occam_mismatched[s'])`.
    """
    all_sources = sorted(set(seen_sources) | set(unseen_sources))
    d_obs = np.asarray(delta_obs, dtype=np.float64)
    tm = np.asarray(delta_transfer, dtype=np.float64)
    r_obs = d_obs - tm
    ha_res = np.asarray(delta_half_a, dtype=np.float64) - tm
    hb_res = np.asarray(delta_half_b, dtype=np.float64) - tm

    rows, idx_list = _extract_masked_idx_per_pert(
        d_obs, target_indices, gene_keep_mask=None, mode="all_nontarget", top_k=50
    )
    m = len(rows)
    if m == 0:
        raise ValueError("Zero valid perturbations for H-I3b contrast.")

    draws = np.random.default_rng(int(seed)).integers(0, m, size=(int(n_boot), m))

    # Build ordered list of residual prediction matrices so every arm shares exact moments & draws
    pred_mats: List[np.ndarray] = [
        np.asarray(occam_matched, dtype=np.float64) - tm,
        np.asarray(occam_resample, dtype=np.float64) - tm,
        np.asarray(state_se_matched, dtype=np.float64) - tm,
        np.asarray(state_se_resample, dtype=np.float64) - tm,
    ]
    occ_mis_idx: Dict[str, int] = {}
    se_mis_idx: Dict[str, int] = {}
    occ_m_src_idx: Dict[str, int] = {}
    for sk in all_sources:
        occ_mis_idx[sk] = len(pred_mats)
        pred_mats.append(np.asarray(occam_mismatched[sk], dtype=np.float64) - tm)
        se_mis_idx[sk] = len(pred_mats)
        pred_mats.append(np.asarray(state_se_mismatched[sk], dtype=np.float64) - tm)
        if occam_matched_by_source is not None and sk in occam_matched_by_source:
            occ_m_src_idx[sk] = len(pred_mats)
            pred_mats.append(np.asarray(occam_matched_by_source[sk], dtype=np.float64) - tm)

    out: Dict[str, Any] = {}
    for metric_key, b_s in (("r_orth_cal", None), ("r_dorth_cal", source_grand_mean)):
        moms, hm = _collect_mode_moments(
            r_obs, pred_mats, tm, rows, idx_list, b_s, ha_res, hb_res, half_valid
        )
        zero_null = b_s is not None
        pt_om, b_om, _, pt_rho, pt_ceil = _eval_boot_cal(moms[0], hm, draws, single_half=False, zero_null_to_zero=zero_null)
        pt_or, b_or, _, _, _ = _eval_boot_cal(moms[1], hm, draws, single_half=False, zero_null_to_zero=zero_null)
        pt_sm, b_sm, _, _, _ = _eval_boot_cal(moms[2], hm, draws, single_half=False, zero_null_to_zero=zero_null)
        pt_sr, b_sr, _, _, _ = _eval_boot_cal(moms[3], hm, draws, single_half=False, zero_null_to_zero=zero_null)

        pt_occ_resamp = float(pt_om - pt_or)
        b_occ_resamp = b_om - b_or
        pt_se_resamp = float(pt_sm - pt_sr)
        b_se_resamp = b_sm - b_sr

        per_source_out: Dict[str, Dict[str, Any]] = {}
        dp_occ_pt: Dict[str, float] = {}
        dp_occ_b: Dict[str, np.ndarray] = {}
        dp_se_pt: Dict[str, float] = {}
        dp_se_b: Dict[str, np.ndarray] = {}
        ex_occ_pt: Dict[str, float] = {}
        ex_occ_b: Dict[str, np.ndarray] = {}
        ex_se_pt: Dict[str, float] = {}
        ex_se_b: Dict[str, np.ndarray] = {}

        for sk in all_sources:
            pt_os, b_os, _, _, _ = _eval_boot_cal(
                moms[occ_mis_idx[sk]], hm, draws, single_half=False, zero_null_to_zero=zero_null
            )
            pt_ss, b_ss, _, _, _ = _eval_boot_cal(
                moms[se_mis_idx[sk]], hm, draws, single_half=False, zero_null_to_zero=zero_null
            )
            if sk in occ_m_src_idx:
                pt_om_s, b_om_s, _, _, _ = _eval_boot_cal(
                    moms[occ_m_src_idx[sk]], hm, draws, single_half=False, zero_null_to_zero=zero_null
                )
            else:
                pt_om_s, b_om_s = pt_om, b_om

            dp_o_p = float(pt_om_s - pt_os)
            dp_o_b = b_om_s - b_os
            dp_s_p = float(pt_sm - pt_ss)
            dp_s_b = b_sm - b_ss

            ex_o_p = float(dp_o_p - pt_occ_resamp)
            ex_o_b = dp_o_b - b_occ_resamp
            ex_s_p = float(dp_s_p - pt_se_resamp)
            ex_s_b = dp_s_b - b_se_resamp

            dp_occ_pt[sk], dp_occ_b[sk] = dp_o_p, dp_o_b
            dp_se_pt[sk], dp_se_b[sk] = dp_s_p, dp_s_b
            ex_occ_pt[sk], ex_occ_b[sk] = ex_o_p, ex_o_b
            ex_se_pt[sk], ex_se_b[sk] = ex_s_p, ex_s_b

            spec_sum = summarize_boot_r16(float(dp_o_p - dp_s_p), dp_o_b - dp_s_b, m, n_boot)
            spec_sum["classification"] = _classify_two_sided_975(
                spec_sum,
                pos_label="OCCAM_MORE_CONTROL_SPECIFIC",
                neg_label="STATE_SE_MORE_CONTROL_SPECIFIC",
            )
            spec_sum["occam_delta_pair"] = dp_o_p
            spec_sum["state_se_delta_pair"] = dp_s_p

            spec_ex_sum = summarize_boot_r16(float(ex_o_p - ex_s_p), ex_o_b - ex_s_b, m, n_boot)
            spec_ex_sum["classification"] = _classify_two_sided_975(
                spec_ex_sum,
                pos_label="OCCAM_MORE_CONTROL_SPECIFIC",
                neg_label="STATE_SE_MORE_CONTROL_SPECIFIC",
            )
            spec_ex_sum["occam_delta_excess"] = ex_o_p
            spec_ex_sum["state_se_delta_excess"] = ex_s_p

            stratum = "seen_by_checkpoint" if sk in seen_sources else "unseen_by_checkpoint"
            per_source_out[sk] = {
                "stratum": stratum,
                "delta_spec": spec_sum,
                "delta_spec_excess": spec_ex_sum,
            }

        per_stratum_out: Dict[str, Dict[str, Any]] = {}
        for st_name, st_srcs in (
            ("seen_by_checkpoint", sorted(seen_sources)),
            ("unseen_by_checkpoint", sorted(unseen_sources)),
        ):
            if not st_srcs:
                continue
            occ_st_pt = float(np.mean([dp_occ_pt[s] for s in st_srcs]))
            se_st_pt = float(np.mean([dp_se_pt[s] for s in st_srcs]))
            occ_st_b = np.mean(np.stack([dp_occ_b[s] for s in st_srcs], axis=0), axis=0)
            se_st_b = np.mean(np.stack([dp_se_b[s] for s in st_srcs], axis=0), axis=0)

            st_spec = summarize_boot_r16(float(occ_st_pt - se_st_pt), occ_st_b - se_st_b, m, n_boot)
            st_spec["classification"] = _classify_two_sided_975(
                st_spec,
                pos_label="OCCAM_MORE_CONTROL_SPECIFIC",
                neg_label="STATE_SE_MORE_CONTROL_SPECIFIC",
            )
            st_spec["occam_stratum_mean"] = occ_st_pt
            st_spec["state_se_stratum_mean"] = se_st_pt

            occ_ex_pt = float(np.mean([ex_occ_pt[s] for s in st_srcs]))
            se_ex_pt = float(np.mean([ex_se_pt[s] for s in st_srcs]))
            occ_ex_b = np.mean(np.stack([ex_occ_b[s] for s in st_srcs], axis=0), axis=0)
            se_ex_b = np.mean(np.stack([ex_se_b[s] for s in st_srcs], axis=0), axis=0)

            st_ex = summarize_boot_r16(float(occ_ex_pt - se_ex_pt), occ_ex_b - se_ex_b, m, n_boot)
            st_ex["classification"] = _classify_two_sided_975(
                st_ex,
                pos_label="OCCAM_MORE_CONTROL_SPECIFIC",
                neg_label="STATE_SE_MORE_CONTROL_SPECIFIC",
            )
            st_ex["occam_stratum_excess_mean"] = occ_ex_pt
            st_ex["state_se_stratum_excess_mean"] = se_ex_pt

            per_stratum_out[st_name] = {
                "sources": st_srcs,
                "delta_spec": st_spec,
                "delta_spec_excess": st_ex,
            }

        out[metric_key] = {
            "rho_ab_orth": float(pt_rho),
            "ceiling": float(pt_ceil),
            "double_orth": bool(zero_null),
            "per_source": per_source_out,
            "per_stratum_mean": per_stratum_out,
        }
    return out


def classify_advance_interpretation_case_addb(
    *,
    h_i1_state_se: str,
    h_i1_occam: str,
    h_i3a_call: str,
    occam_exceeds_resample: bool,
    state_se_exceeds_resample: bool,
) -> str:
    """
    Read-only descriptive mapper into the pre-locked `§4.3` Advance Outcome Interpretation Matrix
    (Cases 1–6). Strictly descriptive; never modifies `H-I1_Occam`, `H-I2_Occam`, `H-I3a`, or `H-I3b`.
    """
    occ_specific = (h_i1_occam == "LINE_MATCHED_SPECIFIC") and occam_exceeds_resample
    se_specific = (h_i1_state_se == "LINE_MATCHED_SPECIFIC") and state_se_exceeds_resample
    fail_states = (
        "LINE_EQUIVALENT_NONZERO_SIGNAL",
        "LINE_EQUIVALENT_NULL_SIGNAL",
        "LINE_MATCHED_SPECIFIC_SUB_DELTA",
    )
    occ_failed = (h_i1_occam in fail_states) and not occ_specific
    se_failed = (h_i1_state_se in fail_states) and not se_specific

    if occ_specific and se_specific:
        return "Case_5_Both_Matched_Specific"
    if occ_specific and not se_specific:
        return "Case_2_Occam_Specific_Not_STATE_SE"
    if se_specific and not occ_specific:
        return "Case_3_STATE_SE_Specific_Not_Occam"
    if h_i3a_call == "OCCAM_SUPERIOR_TO_STATE_SE" and occ_failed:
        return "Case_4_Occam_Higher_Accuracy_Without_Control_Specificity"
    if occ_failed and se_failed:
        return "Case_1_Cross_Architecture_Control_Underutilization"
    return "Case_6_Mixed_Inconclusive_Or_Mismatched_Better"


def _evaluate_one_axis_suite_addb(
    *,
    delta_obs: np.ndarray,
    delta_half_a: np.ndarray,
    delta_half_b: np.ndarray,
    half_valid: np.ndarray,
    delta_transfer: np.ndarray,
    source_grand_mean: Optional[np.ndarray],
    fct_delta_pred: np.ndarray,
    target_indices: Sequence[int],
    pred_matched: np.ndarray,
    pred_mismatched: Mapping[str, np.ndarray],
    pred_resample: Optional[np.ndarray],
    seen_sources: Sequence[str],
    unseen_sources: Sequence[str],
    n_boot: int,
    seed: int,
    delta_tost: float,
    h_i1_key: str,
    h_i2_key: str,
) -> Dict[str, Any]:
    """Evaluate matched, paired vs mismatched, H-I1, FCT, H-I2, and optional within-line resample for one arm."""
    occ_m = compute_orth_suite_r16(
        delta_obs, pred_matched, delta_transfer, target_indices,
        source_grand_mean=source_grand_mean, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
        half_valid=half_valid, mode="all_nontarget", n_boot=n_boot, seed=seed,
    )
    paired_mis = compute_paired_model_vs_all_mismatched_r16(
        delta_obs, pred_matched, dict(pred_mismatched), delta_transfer, target_indices,
        source_grand_mean=source_grand_mean, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
        half_valid=half_valid, n_boot=n_boot, seed=seed,
    )
    h_i1 = adjudicate_claim_i1_line(
        occ_m, paired_mis, seen_sources=seen_sources, unseen_sources=unseen_sources, delta=delta_tost,
    )
    fct_m = compute_orth_suite_r16(
        delta_obs, fct_delta_pred, delta_transfer, target_indices,
        source_grand_mean=source_grand_mean, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
        half_valid=half_valid, mode="all_nontarget", n_boot=n_boot, seed=seed,
    )
    vs_fct = compute_paired_contrast_r16(
        delta_obs, pred_matched, fct_delta_pred, delta_transfer, target_indices,
        source_grand_mean=source_grand_mean, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
        half_valid=half_valid, mode="all_nontarget", n_boot=n_boot, seed=seed,
    )
    h_i2 = adjudicate_claim_i2_line(vs_fct)
    h_i2["fct_matched_r_orth_cal"] = float(fct_m["estimate"])
    h_i2["fct_matched_ci975_lo"] = float(fct_m["ci975_lo"])
    h_i2["fct_matched_ci975_hi"] = float(fct_m["ci975_hi"])
    h_i2["anti_masking_m62_note"] = (
        "PREREG_I §3.2 M6.2: H-I2 MUST be reported line-by-line alongside FCT; never let weak-FCT "
        "results mask strong-FCT inferiority or H-I1 control-specificity failure."
    )

    resamp_block: Optional[Dict[str, Any]] = None
    if pred_resample is not None:
        raw_res = compute_paired_with_within_line_resample_r16(
            delta_obs, pred_matched, pred_resample, dict(pred_mismatched), delta_transfer, target_indices,
            source_grand_mean=source_grand_mean, delta_half_a=delta_half_a, delta_half_b=delta_half_b,
            half_valid=half_valid, n_boot=n_boot, seed=seed,
        )
        resamp_block = summarize_within_line_resample_r16(
            raw_res,
            seen_sources=seen_sources,
            unseen_sources=unseen_sources,
            matched_positive_ci975=bool(h_i1["matched_positive_ci975"]),
        )

    return {
        "Occam_matched_orth": occ_m,
        "Occam_paired_vs_mismatched": paired_mis,
        h_i1_key: h_i1,
        "FCT_matched_orth": fct_m,
        "Occam_vs_FCT": vs_fct,
        h_i2_key: h_i2,
        "within_line_resample_reference": resamp_block,
    }


def _summarize_common_axis_stratum_cb1_2(
    stratum_sources: Sequence[str],
    per_pair_calls: Mapping[str, Dict[str, Any]],
    matched_arm_valid: bool,
    matched_positive: bool,
) -> Dict[str, Any]:
    """Summarize one mismatched-control stratum for `H_I1_Occam_common_axis` (`CB1.2(c)`)."""
    counts = {
        "MATCHED_BETTER": 0,
        "EQUIVALENT": 0,
        "EXACT_INVARIANCE": 0,
        "INCONCLUSIVE": 0,
        "MISMATCHED_BETTER": 0,
    }
    pair_calls: Dict[str, Dict[str, Any]] = {}
    for sk in sorted(stratum_sources):
        c = dict(per_pair_calls[sk])
        pair_calls[sk] = c
        counts[str(c["category"])] += 1

    n_pairs = len(pair_calls)
    all_within = bool(n_pairs > 0 and all(bool(c["within_delta"]) for c in pair_calls.values()))
    any_within = bool(any(bool(c["within_delta"]) for c in pair_calls.values()))
    n_equiv_effective = counts["EQUIVALENT"] + counts["EXACT_INVARIANCE"]

    if n_pairs == 0:
        status = "STRATUM_EMPTY"
    elif not matched_arm_valid:
        status = "STRATUM_INCONCLUSIVE"
    elif matched_positive and counts["MATCHED_BETTER"] == n_pairs:
        status = "STRATUM_MATCHED_SPECIFIC_SUB_DELTA" if all_within else "STRATUM_MATCHED_SPECIFIC"
    elif n_equiv_effective == n_pairs:
        status = "STRATUM_EQUIVALENT_NONZERO_SIGNAL" if matched_positive else "STRATUM_EQUIVALENT_NULL_SIGNAL"
    elif counts["MISMATCHED_BETTER"] == n_pairs:
        status = "STRATUM_MISMATCHED_BETTER"
    else:
        status = "STRATUM_MIXED_OR_INCONCLUSIVE"

    return {
        "n_pairs": n_pairs,
        "sources": sorted(stratum_sources),
        "counts": counts,
        "all_within_delta": all_within,
        "any_within_delta": any_within,
        "per_source": pair_calls,
        "stratum_summary": status,
    }


def _adjudicate_common_axis_h_i1_cb1_2(
    primary_matched_orth: Dict[str, Any],
    per_pair_calls: Mapping[str, Dict[str, Any]],
    seen_sources: Sequence[str],
    unseen_sources: Sequence[str],
    delta: float = DEFAULT_TOST_DELTA,
) -> Dict[str, Any]:
    """
    Adjudicate line-level `H_I1_Occam_common_axis` per binding ruling `CB1.2` + `Amendment CB1.2a`:
      (a) Per-pair categories in `per_pair_calls` have already been gated by each pair's own
          masked matched arm validity (`valid_frac >= 0.95`, `rho_ab_orth > 0`, finite `97.5%` CI;
          positivity is NOT checked per pair) for `n_union_absent > 0`, while pairs reusing the
          primary axis (`n_union_absent == 0`) retain their primary category.
      (b) The line-level matched-arm gate (validity + positivity) uses the primary `rhj_6983_axis`
          matched arm.
      (c) The exact `PREREG_I` + Addendum A stratum/line cascade is applied to `per_pair_calls`.
    """
    seen_set = set(seen_sources)
    unseen_set = set(unseen_sources)
    overlap = seen_set & unseen_set
    if overlap:
        raise ValueError(f"Mismatched strata must be disjoint; overlap={sorted(overlap)}")

    observed_keys = set(per_pair_calls.keys())
    expected_keys = seen_set | unseen_set
    if observed_keys != expected_keys:
        missing = sorted(expected_keys - observed_keys)
        extra = sorted(observed_keys - expected_keys)
        raise ValueError(
            "per_pair_calls keys must match seen_sources | unseen_sources exactly "
            f"(missing={missing}, extra={extra})."
        )

    vf = float(primary_matched_orth.get("valid_frac", 0.0))
    rho = float(primary_matched_orth.get("rho_ab_orth", float("nan")))
    lo975 = float(primary_matched_orth.get("ci975_lo", float("nan")))
    hi975 = float(primary_matched_orth.get("ci975_hi", float("nan")))
    matched_arm_valid = bool(
        vf >= MIN_BOOT_VALID_FRAC
        and np.isfinite(rho)
        and rho > 0.0
        and np.isfinite(lo975)
        and np.isfinite(hi975)
    )
    matched_positive = bool(matched_arm_valid and lo975 > 0.0)

    seen_res = _summarize_common_axis_stratum_cb1_2(
        seen_sources, per_pair_calls, matched_arm_valid, matched_positive
    )
    unseen_res = _summarize_common_axis_stratum_cb1_2(
        unseen_sources, per_pair_calls, matched_arm_valid, matched_positive
    )

    combined_counts = {
        k: seen_res["counts"][k] + unseen_res["counts"][k]
        for k in ("MATCHED_BETTER", "EQUIVALENT", "EXACT_INVARIANCE", "INCONCLUSIVE", "MISMATCHED_BETTER")
    }
    combined_per_source: Dict[str, Dict[str, Any]] = {}
    for sk, val in seen_res["per_source"].items():
        combined_per_source[sk] = {**val, "stratum": "seen_by_checkpoint"}
    for sk, val in unseen_res["per_source"].items():
        combined_per_source[sk] = {**val, "stratum": "unseen_by_checkpoint"}

    n_total = seen_res["n_pairs"] + unseen_res["n_pairs"]
    all_within_delta = bool(
        n_total > 0 and all(bool(v["within_delta"]) for v in combined_per_source.values())
    )
    any_within_delta = bool(any(bool(v["within_delta"]) for v in combined_per_source.values()))

    s_stat = seen_res["stratum_summary"]
    u_stat = unseen_res["stratum_summary"]
    mb_states = ("STRATUM_MATCHED_SPECIFIC", "STRATUM_MATCHED_SPECIFIC_SUB_DELTA")

    if s_stat == "STRATUM_EMPTY" or u_stat == "STRATUM_EMPTY":
        line_summary = "MISSING_STRATUM"
    elif not matched_arm_valid:
        line_summary = "INCONCLUSIVE"
    elif s_stat in mb_states and u_stat in mb_states:
        line_summary = (
            "LINE_MATCHED_SPECIFIC_SUB_DELTA" if all_within_delta else "LINE_MATCHED_SPECIFIC"
        )
    elif s_stat == "STRATUM_EQUIVALENT_NONZERO_SIGNAL" and u_stat == "STRATUM_EQUIVALENT_NONZERO_SIGNAL":
        line_summary = "LINE_EQUIVALENT_NONZERO_SIGNAL"
    elif s_stat == "STRATUM_EQUIVALENT_NULL_SIGNAL" and u_stat == "STRATUM_EQUIVALENT_NULL_SIGNAL":
        line_summary = "LINE_EQUIVALENT_NULL_SIGNAL"
    elif s_stat == "STRATUM_MISMATCHED_BETTER" and u_stat == "STRATUM_MISMATCHED_BETTER":
        line_summary = "LINE_MISMATCHED_BETTER"
    else:
        line_summary = "MIXED_OR_INCONCLUSIVE"

    if line_summary == "LINE_MATCHED_SPECIFIC_SUB_DELTA":
        interp_note = SUB_DELTA_INTERPRETATION_NOTE
    elif line_summary == "LINE_EQUIVALENT_NONZERO_SIGNAL":
        interp_note = (
            "Non-zero matched signal with matched-vs-mismatched equivalence within ±delta "
            "(EXACT_INVARIANCE pairs counted as EQUIVALENT per Addendum A Item 10)."
        )
    elif line_summary == "LINE_EQUIVALENT_NULL_SIGNAL":
        interp_note = (
            "Trivial null-signal equivalence (matched arm not strictly positive; NEVER report as "
            "non-zero signal ignoring target controls)."
        )
    elif line_summary in ("INCONCLUSIVE", "MIXED_OR_INCONCLUSIVE"):
        interp_note = (
            "Any INCONCLUSIVE pair reflects insufficient statistical power, NEVER equivalence ruled out "
            "(PREREG_I §3.3)."
        )
    else:
        interp_note = ""

    return {
        "matched_arm_valid": matched_arm_valid,
        "matched_positive_ci975": matched_positive,
        "rho_ab_orth": rho,
        "n_pairs": n_total,
        "delta": float(delta),
        "all_within_delta": all_within_delta,
        "any_within_delta": any_within_delta,
        "strata": {
            "seen_by_checkpoint": seen_res,
            "unseen_by_checkpoint": unseen_res,
        },
        "counts": combined_counts,
        "per_source": combined_per_source,
        "line_summary": line_summary,
        "interpretation_note": interp_note,
    }


def evaluate_addb_line_from_arrays(
    *,
    line_name: str,
    observed: Mapping[str, Any],
    occam_arrays: Mapping[str, Any],
    state_se_arrays: Mapping[str, Any],
    pert_resolution_audit: Optional[Mapping[str, Any]] = None,
    output_gene_projection_audit: Optional[Mapping[str, Any]] = None,
    extra_fct_sens: Optional[Mapping[str, Any]] = None,
    n_boot: int = LOCKED_N_BOOT,
    seed: int = DEFAULT_R16_SEED,
    delta_tost: float = DEFAULT_TOST_DELTA,
) -> Dict[str, Any]:
    """
    Evaluate every primary, paired `STATE_SE`, and secondary Addendum B statistic for a single
    target line strictly from arrays (zero model or raw-count I/O).
    """
    spec = ADDB_LINE_SPECS[line_name]
    seen_sources = list(spec["seen_by_checkpoint"])
    unseen_sources = list(spec["unseen_by_checkpoint"])
    all_sources = sorted(set(seen_sources) | set(unseen_sources))

    delta_obs = np.asarray(observed["delta_obs"], dtype=np.float64)
    delta_ha = np.asarray(observed["delta_half_a"], dtype=np.float64)
    delta_hb = np.asarray(observed["delta_half_b"], dtype=np.float64)
    half_valid = np.asarray(observed["half_valid"], dtype=bool)
    delta_tm = np.asarray(observed["delta_transfer"], dtype=np.float64)
    b_s = np.asarray(observed["source_grand_mean"], dtype=np.float64)
    fct_pred = np.asarray(observed["fct_delta_pred"], dtype=np.float64)
    target_indices = np.asarray(observed["target_indices"], dtype=np.int64)
    target_names = [str(x) for x in observed["target_names"]]
    gene_names = [str(x) for x in observed["gene_names"]]

    # 1. Primary: rhj_6983_axis, baseline="reconstructed", r_orth_cal
    prim_occ = occam_arrays["rhj_6983_axis_reconstructed"]
    occ_m = np.asarray(prim_occ["matched"], dtype=np.float64)
    occ_r = np.asarray(prim_occ["resample"], dtype=np.float64)
    occ_mis = {sk: np.asarray(prim_occ["mismatched"][sk], dtype=np.float64) for sk in all_sources}

    primary_suite = _evaluate_one_axis_suite_addb(
        delta_obs=delta_obs,
        delta_half_a=delta_ha,
        delta_half_b=delta_hb,
        half_valid=half_valid,
        delta_transfer=delta_tm,
        source_grand_mean=None,
        fct_delta_pred=fct_pred,
        target_indices=target_indices,
        pred_matched=occ_m,
        pred_mismatched=occ_mis,
        pred_resample=occ_r,
        seen_sources=seen_sources,
        unseen_sources=unseen_sources,
        n_boot=n_boot,
        seed=seed,
        delta_tost=delta_tost,
        h_i1_key="H_I1_Occam",
        h_i2_key="H_I2_Occam",
    )

    # 2. Pre-registered paired STATE_SE vs Occam contrasts (§4.2: H-I3a & H-I3b)
    se_m = np.asarray(state_se_arrays["matched"], dtype=np.float64)
    se_r = np.asarray(state_se_arrays["resample"], dtype=np.float64)
    se_mis = {sk: np.asarray(state_se_arrays["mismatched"][sk], dtype=np.float64) for sk in all_sources}

    se_ref_suite = _evaluate_one_axis_suite_addb(
        delta_obs=delta_obs,
        delta_half_a=delta_ha,
        delta_half_b=delta_hb,
        half_valid=half_valid,
        delta_transfer=delta_tm,
        source_grand_mean=None,
        fct_delta_pred=fct_pred,
        target_indices=target_indices,
        pred_matched=se_m,
        pred_mismatched=se_mis,
        pred_resample=se_r,
        seen_sources=seen_sources,
        unseen_sources=unseen_sources,
        n_boot=n_boot,
        seed=seed,
        delta_tost=delta_tost,
        h_i1_key="H_I1_STATE_SE",
        h_i2_key="H_I2_STATE_SE",
    )

    h_i3a = compute_h_i3a_contrast_addb(
        delta_obs, occ_m, se_m, delta_tm, target_indices, b_s, delta_ha, delta_hb, half_valid,
        n_boot=n_boot, seed=seed,
    )
    h_i3b = compute_h_i3b_contrast_addb(
        delta_obs, occ_m, occ_r, occ_mis, se_m, se_r, se_mis,
        seen_sources, unseen_sources, delta_tm, target_indices, b_s, delta_ha, delta_hb, half_valid,
        n_boot=n_boot, seed=seed,
    )
    adv_case = classify_advance_interpretation_case_addb(
        h_i1_state_se=str(se_ref_suite["H_I1_STATE_SE"]["line_summary"]),
        h_i1_occam=str(primary_suite["H_I1_Occam"]["line_summary"]),
        h_i3a_call=str(h_i3a["r_orth_cal"]["classification"]),
        occam_exceeds_resample=bool(
            primary_suite["within_line_resample_reference"]["resample_comparison_status"]
            == "EXCEEDS_WITHIN_LINE_RESAMPLE"
        ),
        state_se_exceeds_resample=bool(
            se_ref_suite["within_line_resample_reference"]["resample_comparison_status"]
            == "EXCEEDS_WITHIN_LINE_RESAMPLE"
        ),
    )

    # 3. Secondary: r_dorth_cal mirror
    dorth_suite = _evaluate_one_axis_suite_addb(
        delta_obs=delta_obs,
        delta_half_a=delta_ha,
        delta_half_b=delta_hb,
        half_valid=half_valid,
        delta_transfer=delta_tm,
        source_grand_mean=b_s,
        fct_delta_pred=fct_pred,
        target_indices=target_indices,
        pred_matched=occ_m,
        pred_mismatched=occ_mis,
        pred_resample=occ_r,
        seen_sources=seen_sources,
        unseen_sources=unseen_sources,
        n_boot=n_boot,
        seed=seed,
        delta_tost=delta_tost,
        h_i1_key="H_I1_Occam_dorth",
        h_i2_key="H_I2_Occam_dorth",
    )

    # 4. Secondary: P4 mandatory pairwise_common_axis (§3.3.4 + binding ruling CB1.2)
    p4_in = occam_arrays.get("pairwise_common_axis_reconstructed", {})
    p4_per_source: Dict[str, Dict[str, Any]] = {}
    p4_pair_calls_for_adj: Dict[str, Dict[str, Any]] = {}
    p4_occ_m_by_src: Dict[str, np.ndarray] = {}
    p4_occ_mis_by_src: Dict[str, np.ndarray] = {}
    any_verdict_changed = False

    for sk in all_sources:
        prim_pair_dict = primary_suite["H_I1_Occam"]["per_source"][sk]
        prim_pair_cat = str(prim_pair_dict["category"])
        if sk in p4_in:
            pair_rec = p4_in[sk]
            pm_arr = np.asarray(pair_rec["matched"], dtype=np.float64)
            ps_arr = np.asarray(pair_rec["mismatched"], dtype=np.float64)
            n_common = int(pair_rec["n_active_rhj_genes"])
            n_union_abs = int(6983 - n_common) if n_common <= 6983 else 0
            mask_sha = str(pair_rec.get("axis_mask_sha256", ""))
            reused_prim = False
        else:
            pm_arr = occ_m
            ps_arr = occ_mis[sk]
            n_common = 6983
            n_union_abs = 0
            mask_sha = ""
            reused_prim = True

        p4_occ_m_by_src[sk] = pm_arr
        p4_occ_mis_by_src[sk] = ps_arr

        c_pair = compute_paired_contrast_r16(
            delta_obs, pm_arr, ps_arr, delta_tm, target_indices,
            source_grand_mean=None, delta_half_a=delta_ha, delta_half_b=delta_hb,
            half_valid=half_valid, mode="all_nontarget", n_boot=n_boot, seed=seed,
        )
        mad_p4, exact_p4 = _array_max_abs_diff_and_exact_inv(pm_arr, ps_arr)
        pair_sum_entry = {
            "estimate": float(c_pair["estimate"]),
            "ci95_lo": float(c_pair["ci95_lo"]),
            "ci95_hi": float(c_pair["ci95_hi"]),
            "ci975_lo": float(c_pair["ci975_lo"]),
            "ci975_hi": float(c_pair["ci975_hi"]),
            "boot_se": float(c_pair["boot_se"]),
            "valid_frac": float(c_pair["valid_frac"]),
            "n_items": int(c_pair["n_items"]),
            "n_boot": int(c_pair["n_boot"]),
            "mismatched_estimate": float(c_pair["b_summary"]["estimate"]),
            "max_abs_diff": mad_p4,
            "exact_invariance": exact_p4,
        }

        m_arm = c_pair["a_summary"]
        m_vf = float(m_arm.get("valid_frac", 0.0))
        m_rho = float(c_pair.get("rho_ab_orth", float("nan")))
        m_lo975 = float(m_arm.get("ci975_lo", float("nan")))
        m_hi975 = float(m_arm.get("ci975_hi", float("nan")))
        masked_matched_arm_valid = bool(
            m_vf >= MIN_BOOT_VALID_FRAC
            and np.isfinite(m_rho)
            and m_rho > 0.0
            and np.isfinite(m_lo975)
            and np.isfinite(m_hi975)
        )
        masked_matched_arm_positive = bool(masked_matched_arm_valid and m_lo975 > 0.0)

        if reused_prim:
            # CB1.2(c): Pairs that reuse the primary axis (n_union_absent == 0) take their primary category.
            p4_cls = {k: v for k, v in prim_pair_dict.items() if k != "stratum"}
            p4_cls["cb1_2_masked_matched_gate_failed"] = False
            raw_cat = prim_pair_cat
            cb1_2_forced_inconclusive = False
        else:
            # Corrected CB1.2(a) (Amendment CB1.2a): Classify from pairwise_common_axis contrast;
            # force INCONCLUSIVE iff own masked matched arm fails validity (`valid_frac >= 0.95`,
            # `rho_ab_orth > 0`, finite `97.5%` CI). Positivity is NOT checked per pair (it enters
            # only through the line-level primary matched-arm gate of CB1.2(b)).
            raw_p4_cls = classify_pair_tost_r16(pair_sum_entry, delta=delta_tost)
            raw_cat = str(raw_p4_cls["category"])
            p4_cls = dict(raw_p4_cls)
            if not masked_matched_arm_valid:
                cb1_2_forced_inconclusive = bool(raw_cat != "INCONCLUSIVE")
                p4_cls["category"] = "INCONCLUSIVE"
                p4_cls["within_delta"] = False
                p4_cls["cb1_2_masked_matched_gate_failed"] = True
            else:
                cb1_2_forced_inconclusive = False
                p4_cls["cb1_2_masked_matched_gate_failed"] = False

        p4_pair_calls_for_adj[sk] = p4_cls
        changed = bool(str(p4_cls["category"]) != prim_pair_cat)
        if changed:
            any_verdict_changed = True
        p4_per_source[sk] = {
            "stratum": "seen_by_checkpoint" if sk in seen_sources else "unseen_by_checkpoint",
            "n_common": n_common,
            "n_union_absent": n_union_abs,
            "reused_primary_rhj_6983": reused_prim,
            "mask_sha256": mask_sha,
            "matched_orth": {**c_pair["a_summary"], "rho_ab_orth": m_rho},
            "mismatched_orth": c_pair["b_summary"],
            "masked_matched_arm_valid": masked_matched_arm_valid,
            "masked_matched_arm_positive": masked_matched_arm_positive,
            "cb1_2_raw_pair_category": raw_cat,
            "cb1_2_forced_inconclusive": cb1_2_forced_inconclusive,
            "paired_delta": pair_sum_entry,
            "pair_classification": p4_cls,
            "primary_pair_classification": prim_pair_cat,
            "verdict_changed_vs_primary": changed,
        }

    h_i1_common = _adjudicate_common_axis_h_i1_cb1_2(
        primary_suite["Occam_matched_orth"],
        p4_pair_calls_for_adj,
        seen_sources=seen_sources,
        unseen_sources=unseen_sources,
        delta=delta_tost,
    )
    if str(h_i1_common["line_summary"]) != str(primary_suite["H_I1_Occam"]["line_summary"]):
        any_verdict_changed = True

    h_i3b_common = compute_h_i3b_contrast_addb(
        delta_obs, occ_m, occ_r, p4_occ_mis_by_src, se_m, se_r, se_mis,
        seen_sources, unseen_sources, delta_tm, target_indices, b_s, delta_ha, delta_hb, half_valid,
        n_boot=n_boot, seed=seed, occam_matched_by_source=p4_occ_m_by_src,
    )

    p4_block = {
        "cb1_2_rule": "per_pair_masked_matched_validity_gate_with_primary_line_gate_cb1_2a",
        "per_source": p4_per_source,
        "H_I1_Occam_common_axis": h_i1_common,
        "H_I3b_common_axis_crosscheck": h_i3b_common,
        "any_verdict_changed_vs_primary": any_verdict_changed,
    }

    # 5. Secondary: baseline="real" on rhj_6983_axis (§3.4)
    real_occ = occam_arrays["rhj_6983_axis_real"]
    baseline_real_suite = _evaluate_one_axis_suite_addb(
        delta_obs=delta_obs,
        delta_half_a=delta_ha,
        delta_half_b=delta_hb,
        half_valid=half_valid,
        delta_transfer=delta_tm,
        source_grand_mean=None,
        fct_delta_pred=fct_pred,
        target_indices=target_indices,
        pred_matched=np.asarray(real_occ["matched"], dtype=np.float64),
        pred_mismatched={sk: np.asarray(real_occ["mismatched"][sk], dtype=np.float64) for sk in all_sources},
        pred_resample=np.asarray(real_occ["resample"], dtype=np.float64),
        seen_sources=seen_sources,
        unseen_sources=unseen_sources,
        n_boot=n_boot,
        seed=seed,
        delta_tost=delta_tost,
        h_i1_key="H_I1_Occam_real",
        h_i2_key="H_I2_Occam_real",
    )

    # 6. Secondary: full_44728_api (§3.3.5)
    full_occ = occam_arrays["full_44728_api_reconstructed"]
    full_api_suite = _evaluate_one_axis_suite_addb(
        delta_obs=delta_obs,
        delta_half_a=delta_ha,
        delta_half_b=delta_hb,
        half_valid=half_valid,
        delta_transfer=delta_tm,
        source_grand_mean=None,
        fct_delta_pred=fct_pred,
        target_indices=target_indices,
        pred_matched=np.asarray(full_occ["matched"], dtype=np.float64),
        pred_mismatched={sk: np.asarray(full_occ["mismatched"][sk], dtype=np.float64) for sk in all_sources},
        pred_resample=np.asarray(full_occ["resample"], dtype=np.float64),
        seen_sources=seen_sources,
        unseen_sources=unseen_sources,
        n_boot=n_boot,
        seed=seed,
        delta_tost=delta_tost,
        h_i1_key="H_I1_Occam_full44728",
        h_i2_key="H_I2_Occam_full44728",
    )

    # 7. Secondary: practical accuracy + M2 cross-mask shift
    prac_acc = compute_practical_accuracy_metrics_addb(delta_obs, occ_m, target_indices)
    prac_acc["m2_joint_cross_mask"] = compute_joint_cross_mask_r16(
        delta_obs, occ_m, delta_tm, target_indices,
        source_grand_mean=None, delta_half_a=delta_ha, delta_half_b=delta_hb,
        half_valid=half_valid, top_k=50, n_boot=n_boot, seed=seed,
    )

    # 8. Secondary: Tier-1 trained-only sensitivity (e.g., a549_ko N=402 excluding Tier-2 RBM14)
    tier1_sens: Optional[Dict[str, Any]] = None
    tier2_targets = list((pert_resolution_audit or {}).get("tier2_targets", []))
    if tier2_targets:
        t2_set = set(tier2_targets)
        keep_idx = np.asarray([i for i, t in enumerate(target_names) if t not in t2_set], dtype=np.int64)
        if keep_idx.size > 0 and keep_idx.size < len(target_names):
            t1_suite = _evaluate_one_axis_suite_addb(
                delta_obs=delta_obs[keep_idx],
                delta_half_a=delta_ha[keep_idx],
                delta_half_b=delta_hb[keep_idx],
                half_valid=half_valid[keep_idx],
                delta_transfer=delta_tm[keep_idx],
                source_grand_mean=None,
                fct_delta_pred=fct_pred[keep_idx],
                target_indices=target_indices[keep_idx],
                pred_matched=occ_m[keep_idx],
                pred_mismatched={sk: occ_mis[sk][keep_idx] for sk in all_sources},
                pred_resample=occ_r[keep_idx],
                seen_sources=seen_sources,
                unseen_sources=unseen_sources,
                n_boot=n_boot,
                seed=seed,
                delta_tost=delta_tost,
                h_i1_key="H_I1_Occam_tier1",
                h_i2_key="H_I2_Occam_tier1",
            )
            t1_i3a = compute_h_i3a_contrast_addb(
                delta_obs[keep_idx], occ_m[keep_idx], se_m[keep_idx], delta_tm[keep_idx],
                target_indices[keep_idx], b_s, delta_ha[keep_idx], delta_hb[keep_idx], half_valid[keep_idx],
                n_boot=n_boot, seed=seed,
            )
            t1_i3b = compute_h_i3b_contrast_addb(
                delta_obs[keep_idx], occ_m[keep_idx], occ_r[keep_idx],
                {sk: occ_mis[sk][keep_idx] for sk in all_sources},
                se_m[keep_idx], se_r[keep_idx],
                {sk: se_mis[sk][keep_idx] for sk in all_sources},
                seen_sources, unseen_sources, delta_tm[keep_idx], target_indices[keep_idx],
                b_s, delta_ha[keep_idx], delta_hb[keep_idx], half_valid[keep_idx],
                n_boot=n_boot, seed=seed,
            )
            tier1_sens = {
                "n_targets_tier1": int(keep_idx.size),
                "excluded_tier2_targets": sorted(t2_set),
                **t1_suite,
                "H_I3a_tier1": t1_i3a,
                "H_I3b_tier1": t1_i3b,
            }

    # 9. Secondary: extra FCT sensitivity (gse295214_h1 FCT_full_r9_6_incl_vcc2025_h1)
    extra_fct_out: Optional[Dict[str, Any]] = None
    if extra_fct_sens:
        extra_fct_out = {}
        for sens_key, sens_data in sorted(extra_fct_sens.items()):
            s_fct = np.asarray(sens_data["delta_pred"], dtype=np.float64)
            s_tm = np.asarray(sens_data["delta_transfer"], dtype=np.float64)
            s_bs = np.asarray(sens_data["source_grand_mean"], dtype=np.float64)
            fct_o = compute_orth_suite_r16(
                delta_obs, s_fct, s_tm, target_indices,
                source_grand_mean=None, delta_half_a=delta_ha, delta_half_b=delta_hb,
                half_valid=half_valid, mode="all_nontarget", n_boot=n_boot, seed=seed,
            )
            occ_o = compute_orth_suite_r16(
                delta_obs, occ_m, s_tm, target_indices,
                source_grand_mean=None, delta_half_a=delta_ha, delta_half_b=delta_hb,
                half_valid=half_valid, mode="all_nontarget", n_boot=n_boot, seed=seed,
            )
            occ_vs_sfct = compute_paired_contrast_r16(
                delta_obs, occ_m, s_fct, s_tm, target_indices,
                source_grand_mean=None, delta_half_a=delta_ha, delta_half_b=delta_hb,
                half_valid=half_valid, mode="all_nontarget", n_boot=n_boot, seed=seed,
            )
            occ_vs_sfct_d = compute_paired_contrast_r16(
                delta_obs, occ_m, s_fct, s_tm, target_indices,
                source_grand_mean=s_bs, delta_half_a=delta_ha, delta_half_b=delta_hb,
                half_valid=half_valid, mode="all_nontarget", n_boot=n_boot, seed=seed,
            )
            extra_fct_out[sens_key] = {
                "fixed_to_primary_n_eval": int(delta_obs.shape[0]),
                "fct_orth": fct_o,
                "occam_orth": occ_o,
                "occam_vs_sens_fct": occ_vs_sfct,
                "H_I2_Occam_sens_fct": adjudicate_claim_i2_line(occ_vs_sfct),
                "occam_vs_sens_fct_dorth": occ_vs_sfct_d,
            }

    return {
        "line_name": line_name,
        "group_key": spec["group_key"],
        "role": spec["role"],
        "modality": spec["modality"],
        "checkpoint_id": OCCAM_CHECKPOINT_ID,
        "n_eval_targets": int(len(target_names)),
        "n_eval_genes": int(len(gene_names)),
        "n_half_valid": int(np.sum(half_valid)),
        "n_self_masked": int(np.sum(target_indices >= 0)),
        "strata_sources": {
            "seen_by_checkpoint": seen_sources,
            "unseen_by_checkpoint": unseen_sources,
        },
        "pert_resolution_audit": dict(pert_resolution_audit) if pert_resolution_audit else {},
        "output_gene_projection_audit": dict(output_gene_projection_audit) if output_gene_projection_audit else {},
        "mandatory_disclosures": list(spec["disclosures"]),
        "primary": primary_suite,
        "paired_vs_STATE_SE": {
            "H_I3a": h_i3a,
            "H_I3b": h_i3b,
            "STATE_SE_reference_calls": {
                "H_I1_line_summary": str(se_ref_suite["H_I1_STATE_SE"]["line_summary"]),
                "H_I2_classification": str(se_ref_suite["H_I2_STATE_SE"]["classification"]),
                "matched_r_orth_cal": float(se_ref_suite["Occam_matched_orth"]["estimate"]),
                "resample_comparison_status": str(
                    se_ref_suite["within_line_resample_reference"]["resample_comparison_status"]
                ),
            },
            "advance_interpretation_case": adv_case,
        },
        "secondary": {
            "r_dorth_cal": dorth_suite,
            "pairwise_common_axis": p4_block,
            "baseline_real": baseline_real_suite,
            "full_44728_api": full_api_suite,
            "practical_accuracy_Occam": prac_acc,
            "tier1_trained_only_sensitivity": tier1_sens,
            "extra_fct_sensitivity": extra_fct_out,
        },
    }


def build_group_adjudications_addb(lines_results: Mapping[str, Mapping[str, Any]]) -> Dict[str, Any]:
    """
    Build strictly separated group adjudications (`§4.1`: Group 2 CRISPRko, Group 2 CRISPRi,
    Group 3 CRISPRi, and Sensitivity-Only `feng_ipsc` — NEVER pooled).
    """
    groups_map: Dict[str, List[str]] = {
        "group2_crisprko": [],
        "group2_crispri": [],
        "group3_crispri": [],
        "sensitivity_only": [],
    }
    for ln in sorted(lines_results.keys()):
        gk = str(lines_results[ln]["group_key"])
        groups_map.setdefault(gk, []).append(ln)

    out: Dict[str, Any] = {}
    for gk, lns in groups_map.items():
        if not lns:
            continue
        verdicts_i1 = {ln: str(lines_results[ln]["primary"]["H_I1_Occam"]["line_summary"]) for ln in lns}
        verdicts_i2 = {ln: str(lines_results[ln]["primary"]["H_I2_Occam"]["classification"]) for ln in lns}
        mods = {ln: str(lines_results[ln]["modality"]) for ln in lns}

        # Use frozen adjudicate_group_r16 per separated group (group2 safeguard applies to total group2 count)
        n_g2_total = len(groups_map.get("group2_crisprko", [])) + len(groups_map.get("group2_crispri", []))
        r16_gkey = "group2" if gk.startswith("group2") else gk
        g_adj = adjudicate_group_r16(
            r16_gkey,
            verdicts_i1,
            mods,
            min_group2_lines=1 if n_g2_total >= 2 else 2,
        )
        i2_counts: Dict[str, int] = {}
        for v in verdicts_i2.values():
            i2_counts[v] = i2_counts.get(v, 0) + 1
        out[gk] = {
            "group_key": gk,
            "K": len(lns),
            "lines": sorted(lns),
            "H_I1_Occam_group_adjudication": g_adj,
            "H_I1_Occam_per_line": verdicts_i1,
            "H_I2_Occam_per_line": verdicts_i2,
            "H_I2_Occam_verdict_counts": i2_counts,
            "H_I2_Occam_k_over_K": {k: f"{v}/{len(lns)}" for k, v in sorted(i2_counts.items())},
            "mandatory_disclosures": sorted(
                {d for ln in lns for d in lines_results[ln].get("mandatory_disclosures", [])}
            ),
        }
    return out


# ==============================================================================
# Stage 4 Independent From-Disk Recomputation (`ADDB_ARRAY_SCHEMA.md §8`)
# ==============================================================================


def _load_r16_state_se_arrays_for_line(
    r16_line_dir: Path,
    all_sources: Sequence[str],
) -> Tuple[Dict[str, Any], Dict[str, str]]:
    """Read frozen Round 16 `STATE_SE` (`fixed_index_zero`) arrays + SHA-256 provenance."""
    ldir = Path(r16_line_dir)
    m_path = ldir / "preds_STATE_SE_fixed_index_zero_matched.npz"
    r_path = ldir / "preds_STATE_SE_fixed_index_zero_resample.npz"
    if not m_path.is_file() or not r_path.is_file():
        raise FileNotFoundError(f"Missing frozen R16 STATE_SE matched/resample files under {ldir}")

    shas: Dict[str, str] = {
        m_path.name: sha256_file(m_path),
        r_path.name: sha256_file(r_path),
    }
    with np.load(m_path, allow_pickle=False) as dm:
        matched_arr = dm["delta_pred"].copy()
    with np.load(r_path, allow_pickle=False) as dr:
        resample_arr = dr["delta_pred"].copy()

    mis_map: Dict[str, np.ndarray] = {}
    for sk in sorted(all_sources):
        sp = ldir / f"preds_STATE_SE_fixed_index_zero_mis_{sk}.npz"
        if not sp.is_file():
            raise FileNotFoundError(f"Missing frozen R16 STATE_SE mismatched file: {sp}")
        shas[sp.name] = sha256_file(sp)
        with np.load(sp, allow_pickle=False) as ds:
            mis_map[sk] = ds["delta_pred"].copy()

    return {"matched": matched_arr, "resample": resample_arr, "mismatched": mis_map}, shas


def recompute_addb_from_exported_arrays(
    arrays_dir: Path,
    r16_state_se_dirs: Mapping[str, Path],
    *,
    tier2_targets_by_line: Optional[Mapping[str, Sequence[str]]] = None,
    n_boot: int = LOCKED_N_BOOT,
    seed: int = DEFAULT_R16_SEED,
    delta_tost: float = DEFAULT_TOST_DELTA,
) -> Dict[str, Any]:
    """
    Stage 4 Controller Verification Helper (`ADDB_ARRAY_SCHEMA.md §8`):
    Recompute EVERY Addendum B metric, CI, `H-I1_Occam`, `H-I2_Occam`, `H-I3a`, `H-I3b`,
    `pairwise_common_axis`, `baseline_real`, `full_44728_api`, `tier1_trained_only_sensitivity`,
    and `extra_fct_sensitivity` STRICTLY from saved `.npz` arrays (zero Occam or raw-count reads).
    """
    adir = Path(arrays_dir)
    line_dirs = sorted(p for p in adir.iterdir() if p.is_dir() and (p / "observed.npz").is_file())
    if not line_dirs:
        raise FileNotFoundError(f"No exported Addendum B line directories found under {adir}")

    lines_out: Dict[str, Any] = {}
    r16_shas_all: Dict[str, Dict[str, str]] = {}

    for ldir in line_dirs:
        ln = ldir.name
        if ln not in ADDB_LINE_SPECS:
            continue
        spec = ADDB_LINE_SPECS[ln]
        all_sources = sorted(set(spec["seen_by_checkpoint"]) | set(spec["unseen_by_checkpoint"]))

        with np.load(ldir / "observed.npz", allow_pickle=False) as zobs:
            obs_dict = {k: zobs[k].copy() for k in zobs.files}

        def _load_arm_set(prefix: str) -> Dict[str, Any]:
            with np.load(ldir / f"preds_{prefix}_matched.npz", allow_pickle=False) as zm:
                m_arr = zm["delta_pred"].copy()
            with np.load(ldir / f"preds_{prefix}_resample.npz", allow_pickle=False) as zr:
                r_arr = zr["delta_pred"].copy()
            mis_d: Dict[str, np.ndarray] = {}
            for sk in all_sources:
                with np.load(ldir / f"preds_{prefix}_mis_{sk}.npz", allow_pickle=False) as zs:
                    mis_d[sk] = zs["delta_pred"].copy()
            return {"matched": m_arr, "resample": r_arr, "mismatched": mis_d}

        p4_map: Dict[str, Dict[str, Any]] = {}
        for sk in all_sources:
            pm_f = ldir / f"preds_Occam_pairwise_common_axis_pair_{sk}_reconstructed_matched.npz"
            ps_f = ldir / f"preds_Occam_pairwise_common_axis_pair_{sk}_reconstructed_mis_{sk}.npz"
            if pm_f.is_file() and ps_f.is_file():
                with np.load(pm_f, allow_pickle=False) as zpm, np.load(ps_f, allow_pickle=False) as zps:
                    p4_map[sk] = {
                        "matched": zpm["delta_pred"].copy(),
                        "mismatched": zps["delta_pred"].copy(),
                        "n_active_rhj_genes": int(zpm["n_active_rhj_genes"].item()),
                        "axis_mask_sha256": str(zpm["axis_mask_sha256"].item()),
                    }

        occam_arrays = {
            "rhj_6983_axis_reconstructed": _load_arm_set("Occam_rhj_6983_axis_reconstructed"),
            "rhj_6983_axis_real": _load_arm_set("Occam_rhj_6983_axis_real"),
            "full_44728_api_reconstructed": _load_arm_set("Occam_full_44728_api_reconstructed"),
            "pairwise_common_axis_reconstructed": p4_map,
        }

        if ln not in r16_state_se_dirs:
            raise KeyError(f"Missing frozen R16 STATE_SE directory for line {ln!r}.")
        se_arrays, se_shas = _load_r16_state_se_arrays_for_line(Path(r16_state_se_dirs[ln]), all_sources)
        r16_shas_all[ln] = se_shas

        extra_fct: Dict[str, Dict[str, np.ndarray]] = {}
        for sf in sorted(ldir.glob("preds_SENS_FCT_*.npz")):
            skey = sf.stem[len("preds_SENS_") :]
            with np.load(sf, allow_pickle=False) as zs:
                extra_fct[skey] = {
                    "delta_pred": zs["delta_pred"].copy(),
                    "delta_transfer": zs["delta_transfer"].copy(),
                    "source_grand_mean": zs["source_grand_mean"].copy(),
                }

        t2_list = list((tier2_targets_by_line or {}).get(ln, []))
        pert_audit = {"tier2_targets": t2_list} if t2_list else None

        lines_out[ln] = evaluate_addb_line_from_arrays(
            line_name=ln,
            observed=obs_dict,
            occam_arrays=occam_arrays,
            state_se_arrays=se_arrays,
            pert_resolution_audit=pert_audit,
            extra_fct_sens=extra_fct if extra_fct else None,
            n_boot=n_boot,
            seed=seed,
            delta_tost=delta_tost,
        )

    return {
        "protocol": "PREREG_I_ADDENDUM_B_v2",
        "model_key": OCCAM_MODEL_KEY,
        "checkpoint_id": OCCAM_CHECKPOINT_ID,
        "seed": int(seed),
        "resample_seed": int(DEFAULT_RESAMPLE_SEED),
        "n_boot": int(n_boot),
        "delta_tost": float(delta_tost),
        "r16_state_se_reference_sha256": r16_shas_all,
        "group_adjudications": build_group_adjudications_addb(lines_out),
        "lines": lines_out,
    }


# ==============================================================================
# Real-Target Evaluation Orchestrator (Strictly Gated + CPU float32 Locked)
# ==============================================================================


def run_real_evaluation_addb(
    target_lines: Sequence[str],
    *,
    results_dir: Path = DEFAULT_ADDB_RESULTS_DIR,
    code_freeze_path: Path = DEFAULT_CODE_FREEZE_PATH,
    go_path: Path = DEFAULT_GO_PATH,
    lock_path: Path = DEFAULT_LOCK_PATH,
    human_approval_path: Path = DEFAULT_HUMAN_APPROVAL_PATH,
    ckpt_path: Path = DEFAULT_OCCAM_CKPT_PATH,
    occam_src_dir: Path = DEFAULT_OCCAM_SRC_DIR,
    gene_axis_csv: Path = DEFAULT_OCCAM_GENE_AXIS_CSV,
    perturbations_csv: Path = DEFAULT_OCCAM_PERTURBATIONS_CSV,
    checksums_file: Path = DEFAULT_OCCAM_CHECKSUMS_FILE,
    raw_nadig_dir: Optional[Path] = DEFAULT_LOCAL_NADIG_DIR,
    gse295214_raw_dir: Optional[Path] = DEFAULT_GSE295214_RAW_DIR,
    device: str = LOCKED_DEVICE,
    num_threads: int = LOCKED_CPU_NUM_THREADS,
    batch_size: int = LOCKED_BATCH_SIZE,
    n_boot: int = LOCKED_N_BOOT,
    seed: int = DEFAULT_R16_SEED,
    delta_tost: float = DEFAULT_TOST_DELTA,
) -> Dict[str, Any]:
    """
    Execute the locked Round 16 Addendum B (`Occam`) evaluation on real targets.
    Strictly enforces `assert_addb_execution_gates` and `device == "cpu"` (`CPU float32`).
    """
    gates_verified = assert_addb_execution_gates(
        code_freeze_path=code_freeze_path,
        go_path=go_path,
        lock_path=lock_path,
        human_approval_path=human_approval_path,
    )
    if str(device).strip().lower() != LOCKED_DEVICE:
        raise ValueError(
            f"Addendum B device is locked to {LOCKED_DEVICE!r} (CPU float32) per Gate 5A Item 3; "
            f"refusing device={device!r}."
        )
    lines_run = validate_target_lines_addb(target_lines)

    t0 = time.time()
    out_dir = Path(results_dir)
    arrays_dir = out_dir / "arrays"
    z_ctrl_dir = out_dir / "z_ctrl"
    arrays_dir.mkdir(parents=True, exist_ok=True)
    z_ctrl_dir.mkdir(parents=True, exist_ok=True)

    configure_occam_runtime(device="cpu", num_threads=int(num_threads), seed=int(seed), allow_mps=False)
    print(
        f"[{_now_iso()}] Starting Addendum B real evaluation on {lines_run} | "
        f"device=cpu threads={num_threads} batch_size={batch_size} n_boot={n_boot} (cpu_t={_cpu_time_s():.1f}s)",
        flush=True,
    )

    # Step 1: Load frozen R10 HGNC bridge, Occam model (Gate 5A Item 1 verified), and locked R16 barcodes
    bridge = HGNCBridge.from_tsv(DEFAULT_HGNC_FILE)
    rhj_axis = load_rhj_gene_axis(gene_axis_csv, verify_sha256=True)
    ptab = load_rhj_perturbations_table(perturbations_csv)

    barcodes_doc = load_consolidated_r16_barcodes_addb(lines_run)
    (out_dir / "control_barcodes_addb.json").write_text(
        json.dumps(barcodes_doc, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    adapter, _ = OccamInferenceAdapter.load_verified(
        ckpt_path=ckpt_path,
        occam_src_dir=occam_src_dir,
        gene_axis_csv=gene_axis_csv,
        perturbations_csv=perturbations_csv,
        checksums_file=checksums_file,
        device="cpu",
        num_threads=int(num_threads),
        allow_mps=False,
        batch_size=int(batch_size),
        verify_file_hashes=True,
    )
    print(f"[{_now_iso()}] Occam checkpoint loaded & Gate 5A.1 verified (cpu_t={_cpu_time_s():.1f}s).", flush=True)

    # Step 2: Load all unique raw control pools (primary + resample) by exact locked R16 barcodes
    needed_primary_pools: Set[str] = set()
    needed_resample_pools: Set[str] = set()
    for ln in lines_run:
        spec = ADDB_LINE_SPECS[ln]
        needed_primary_pools.add(ln)
        needed_resample_pools.add(ln)
        needed_primary_pools.update(spec["seen_by_checkpoint"])
        needed_primary_pools.update(spec["unseen_by_checkpoint"])

    raw_primary: Dict[str, RawControlPool] = {}
    raw_resample: Dict[str, RawControlPool] = {}
    for pn in sorted(needed_primary_pools):
        bcs = barcodes_doc["pools"][pn]["barcodes"]
        raw_primary[pn] = load_raw_control_pool_addb(
            pn, bcs, pool_kind="primary", seed=int(DEFAULT_R16_SEED),
            raw_nadig_dir=raw_nadig_dir, gse295214_raw_dir=gse295214_raw_dir,
        )
        print(
            f"[{_now_iso()}] Loaded raw primary pool {pn}: {raw_primary[pn].raw_counts.shape} "
            f"(cpu_t={_cpu_time_s():.1f}s)",
            flush=True,
        )
    for pn in sorted(needed_resample_pools):
        rbcs = barcodes_doc["pools"][pn]["within_line_resample"]["barcodes"]
        raw_resample[pn] = load_raw_control_pool_addb(
            pn, rbcs, pool_kind="resample", seed=int(DEFAULT_RESAMPLE_SEED),
            raw_nadig_dir=raw_nadig_dir, gse295214_raw_dir=gse295214_raw_dir,
        )
        print(
            f"[{_now_iso()}] Loaded raw resample pool {pn}: {raw_resample[pn].raw_counts.shape} "
            f"(cpu_t={_cpu_time_s():.1f}s)",
            flush=True,
        )

    # Step 3: Encode & cache ALL unique z_ctrl latents across all 3 scoring axes
    z_ctrl_export_entries: Dict[str, Dict[str, Any]] = {}
    encoded_cache: Dict[str, Dict[str, Any]] = {}  # key -> {z_ctrl, x_enc_mean, audit, barcodes, ...}

    def _encode_and_cache(
        cache_key: str,
        pool: RawControlPool,
        x_enc: np.ndarray,
        audit: Dict[str, Any],
        *,
        input_axis_mode: str,
        n_active_genes: int,
        mask_sha: str,
    ) -> Dict[str, Any]:
        z = adapter.encode_controls(x_enc, batch_size=int(batch_size))
        x_mean = x_enc.mean(axis=0, dtype=np.float64)
        rec = {
            "z_ctrl": z,
            "x_enc_mean": x_mean,
            "audit": audit,
            "barcodes": list(pool.barcodes),
            "n_active_genes": int(n_active_genes),
            "mask_sha": str(mask_sha),
        }
        encoded_cache[cache_key] = rec
        z_ctrl_export_entries[cache_key] = {
            "Z_ctrl": z,
            "control_barcodes": pool.barcodes,
            "seed": int(pool.seed),
            "pool_name": pool.line_name,
            "input_axis_mode": input_axis_mode,
            "n_active_rhj_genes": int(n_active_genes),
            "axis_mask_sha256": str(mask_sha),
            "audit": audit,
        }
        return rec

    # 3a. Primary rhj_6983_axis (primary + resample pools)
    for pn in sorted(needed_primary_pools):
        pool = raw_primary[pn]
        x_enc, aud = scatter_rhj_6983_to_encoder(
            pool.raw_counts, pool.raw_ensembl_ids, rhj_axis.ensembl_ids, adapter.encoder_gene_names
        )
        _encode_and_cache(
            f"rhj_6983_axis__{pn}__primary", pool, x_enc, aud,
            input_axis_mode="rhj_6983_axis",
            n_active_genes=int(aud["n_active_rhj_populated"]),
            mask_sha=str(aud["active_rhj_mask_sha256"]),
        )
        del x_enc
    for pn in sorted(needed_resample_pools):
        pool = raw_resample[pn]
        x_enc, aud = scatter_rhj_6983_to_encoder(
            pool.raw_counts, pool.raw_ensembl_ids, rhj_axis.ensembl_ids, adapter.encoder_gene_names
        )
        _encode_and_cache(
            f"rhj_6983_axis__{pn}__resample", pool, x_enc, aud,
            input_axis_mode="rhj_6983_axis",
            n_active_genes=int(aud["n_active_rhj_populated"]),
            mask_sha=str(aud["active_rhj_mask_sha256"]),
        )
        del x_enc

    # 3b. P4 pairwise_common_axis (only genuinely masked primary pools with n_present_but_zeroed_by_mask > 0)
    pair_masks: Dict[Tuple[str, str], Any] = {}
    for ln in lines_run:
        spec = ADDB_LINE_SPECS[ln]
        for sk in sorted(set(spec["seen_by_checkpoint"]) | set(spec["unseen_by_checkpoint"])):
            a, b = (ln, sk) if ln < sk else (sk, ln)
            if (a, b) in pair_masks:
                continue
            pmask = build_pairwise_common_axis_mask(
                rhj_axis.ensembl_ids,
                raw_primary[a].raw_ensembl_ids,
                raw_primary[b].raw_ensembl_ids,
                target_line=a,
                source_line=b,
            )
            pair_masks[(a, b)] = pmask
            if not pmask.requires_common_axis_arm:
                continue
            for side_name in (a, b):
                pool = raw_primary[side_name]
                x_enc, aud = scatter_rhj_6983_to_encoder(
                    pool.raw_counts,
                    pool.raw_ensembl_ids,
                    rhj_axis.ensembl_ids,
                    adapter.encoder_gene_names,
                    active_rhj_mask=pmask.mask,
                )
                if int(aud["n_present_but_zeroed_by_mask"]) > 0:
                    ckey = f"pairwise_common_axis__{side_name}__pair_{a}__{b}__primary"
                    _encode_and_cache(
                        ckey, pool, x_enc, aud,
                        input_axis_mode="pairwise_common_axis",
                        n_active_genes=int(pmask.n_common),
                        mask_sha=str(pmask.mask_sha256),
                    )
                del x_enc

    # 3c. Secondary full_44728_api (primary + resample pools)
    for pn in sorted(needed_primary_pools):
        pool = raw_primary[pn]
        x_enc, aud = scatter_full_44728_api(
            pool.raw_counts, pool.raw_ensembl_ids, adapter.encoder_gene_names
        )
        _encode_and_cache(
            f"full_44728_api__{pn}__primary", pool, x_enc, aud,
            input_axis_mode="full_44728_api",
            n_active_genes=int(aud["n_matched_encoder_genes"]),
            mask_sha="",
        )
        del x_enc
    for pn in sorted(needed_resample_pools):
        pool = raw_resample[pn]
        x_enc, aud = scatter_full_44728_api(
            pool.raw_counts, pool.raw_ensembl_ids, adapter.encoder_gene_names
        )
        _encode_and_cache(
            f"full_44728_api__{pn}__resample", pool, x_enc, aud,
            input_axis_mode="full_44728_api",
            n_active_genes=int(aud["n_matched_encoder_genes"]),
            mask_sha="",
        )
        del x_enc

    z_ctrl_manifest = save_occam_z_ctrl_latents(z_ctrl_dir, z_ctrl_export_entries)
    print(
        f"[{_now_iso()}] Exported {len(z_ctrl_manifest)} unique z_ctrl latent pools "
        f"(cpu_t={_cpu_time_s():.1f}s).",
        flush=True,
    )
    # Raw count matrices are no longer needed once all z_ctrl latents and x_enc_means are cached!
    del raw_primary, raw_resample
    gc.collect()

    # Step 4: Decode per-target arms, save observed.npz + preds_*.npz, and evaluate
    lines_results: Dict[str, Any] = {}
    r16_se_shas_all: Dict[str, Dict[str, str]] = {}

    def _decode_from_cached_latent(
        cache_key: str,
        pert_ids: Sequence[int],
        out_cols: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        rec = encoded_cache[cache_key]
        z = rec["z_ctrl"]
        recon_mean = adapter.predict_zero_pert_from_z(z, out_cols, batch_size=int(batch_size)).mean(
            axis=0, dtype=np.float64
        )
        real_mean = rec["x_enc_mean"][out_cols]
        pred_means = np.empty((len(pert_ids), len(out_cols)), dtype=np.float64)
        for pi, pid in enumerate(pert_ids):
            pred_means[pi] = adapter.predict_single_pert_mean_from_z(
                z, int(pid), out_cols, compose_mode="decode_seq_mean", batch_size=int(batch_size)
            )
        return pred_means - recon_mean[None, :], pred_means - real_mean[None, :]

    for ln in lines_run:
        t_ln = time.time()
        spec = ADDB_LINE_SPECS[ln]
        all_sources = sorted(set(spec["seen_by_checkpoint"]) | set(spec["unseen_by_checkpoint"]))
        r16_ldir = REPO_ROOT / spec["r16_arrays_rel"]
        r16_obs_path = r16_ldir / "observed.npz"
        if not r16_obs_path.is_file():
            raise FileNotFoundError(f"Frozen Round 16 observed.npz not found: {r16_obs_path}")

        with np.load(r16_obs_path, allow_pickle=False) as zobs:
            obs_dict = {k: zobs[k].copy() for k in zobs.files}
        for rk in OBSERVED_NPZ_REQUIRED_KEYS:
            if rk not in obs_dict:
                raise KeyError(f"Frozen Round 16 {r16_obs_path} missing required key {rk!r}.")

        # Verify exact barcode agreement with consolidated locked R16 barcodes
        obs_m_bcs = [str(x) for x in obs_dict["matched_ctrl_barcodes"].tolist()]
        obs_r_bcs = [str(x) for x in obs_dict["resample_ctrl_barcodes"].tolist()]
        if obs_m_bcs != barcodes_doc["pools"][ln]["barcodes"]:
            raise ValueError(f"[{ln}] observed.npz matched_ctrl_barcodes != control_barcodes_r16.json!")
        if obs_r_bcs != barcodes_doc["pools"][ln]["within_line_resample"]["barcodes"]:
            raise ValueError(f"[{ln}] observed.npz resample_ctrl_barcodes != control_barcodes_r16.json!")

        target_names = [str(x) for x in obs_dict["target_names"].tolist()]
        gene_names = [str(x) for x in obs_dict["gene_names"].tolist()]

        # Resolve targets (Two-Tier §3.5) and output genes (RHJ -> G_eval^R16)
        pert_res = resolve_perturbation_keys_two_tier(
            target_names, ptab["eligible_keys"], adapter.pert2id, bridge.canonicalize, require_all_resolved=True
        )
        out_cols, out_gene_audit = build_output_gene_projection(
            rhj_axis, adapter.encoder_gene_names, gene_names, bridge.canonicalize, require_full_coverage=True
        )

        print(
            f"[{_now_iso()}] [{ln}] Decoding N={len(target_names)} (T1={len(pert_res.tier1_targets)}, "
            f"T2={len(pert_res.tier2_targets)}) x G={len(gene_names)} across 3 axes...",
            flush=True,
        )

        pred_arms_to_save: Dict[str, Dict[str, Any]] = {}
        occam_eval_arrays: Dict[str, Any] = {
            "rhj_6983_axis_reconstructed": {"mismatched": {}},
            "rhj_6983_axis_real": {"mismatched": {}},
            "pairwise_common_axis_reconstructed": {},
            "full_44728_api_reconstructed": {"mismatched": {}},
        }

        def _pack_delta(
            d_arr: np.ndarray,
            *,
            ctrl_source: str,
            ctrl_bcs: Sequence[str],
            axis_mode: str,
            base_mode: str,
            pair_cp: str = "",
            n_act: int = 6983,
            mask_sha: str = "",
        ) -> Dict[str, Any]:
            return build_occam_delta_payload(
                delta_pred=d_arr,
                target_names=target_names,
                gene_names=gene_names,
                control_source=ctrl_source,
                control_barcodes=ctrl_bcs,
                input_axis_mode=axis_mode,
                baseline_mode=base_mode,
                resolved_pert_keys=pert_res.resolved_keys,
                resolved_pert_ids=pert_res.resolved_ids,
                pair_counterpart=pair_cp,
                n_active_rhj_genes=n_act,
                axis_mask_sha256=mask_sha,
            )

        # Axis 1: rhj_6983_axis (primary reconstructed + zero-extra-cost real)
        ck_m = f"rhj_6983_axis__{ln}__primary"
        d_rec_m, d_real_m = _decode_from_cached_latent(ck_m, pert_res.resolved_ids, out_cols)
        occam_eval_arrays["rhj_6983_axis_reconstructed"]["matched"] = d_rec_m
        occam_eval_arrays["rhj_6983_axis_real"]["matched"] = d_real_m
        pred_arms_to_save["Occam_rhj_6983_axis_reconstructed_matched"] = _pack_delta(
            d_rec_m, ctrl_source=ln, ctrl_bcs=encoded_cache[ck_m]["barcodes"],
            axis_mode="rhj_6983_axis", base_mode="reconstructed",
            n_act=encoded_cache[ck_m]["n_active_genes"], mask_sha=encoded_cache[ck_m]["mask_sha"],
        )
        pred_arms_to_save["Occam_rhj_6983_axis_real_matched"] = _pack_delta(
            d_real_m, ctrl_source=ln, ctrl_bcs=encoded_cache[ck_m]["barcodes"],
            axis_mode="rhj_6983_axis", base_mode="real",
            n_act=encoded_cache[ck_m]["n_active_genes"], mask_sha=encoded_cache[ck_m]["mask_sha"],
        )

        ck_r = f"rhj_6983_axis__{ln}__resample"
        d_rec_r, d_real_r = _decode_from_cached_latent(ck_r, pert_res.resolved_ids, out_cols)
        occam_eval_arrays["rhj_6983_axis_reconstructed"]["resample"] = d_rec_r
        occam_eval_arrays["rhj_6983_axis_real"]["resample"] = d_real_r
        resamp_src_lbl = f"{ln}_within_line_resample_seed{DEFAULT_RESAMPLE_SEED}"
        pred_arms_to_save["Occam_rhj_6983_axis_reconstructed_resample"] = _pack_delta(
            d_rec_r, ctrl_source=resamp_src_lbl, ctrl_bcs=encoded_cache[ck_r]["barcodes"],
            axis_mode="rhj_6983_axis", base_mode="reconstructed",
            n_act=encoded_cache[ck_r]["n_active_genes"], mask_sha=encoded_cache[ck_r]["mask_sha"],
        )
        pred_arms_to_save["Occam_rhj_6983_axis_real_resample"] = _pack_delta(
            d_real_r, ctrl_source=resamp_src_lbl, ctrl_bcs=encoded_cache[ck_r]["barcodes"],
            axis_mode="rhj_6983_axis", base_mode="real",
            n_act=encoded_cache[ck_r]["n_active_genes"], mask_sha=encoded_cache[ck_r]["mask_sha"],
        )

        for sk in all_sources:
            ck_s = f"rhj_6983_axis__{sk}__primary"
            d_rec_s, d_real_s = _decode_from_cached_latent(ck_s, pert_res.resolved_ids, out_cols)
            occam_eval_arrays["rhj_6983_axis_reconstructed"]["mismatched"][sk] = d_rec_s
            occam_eval_arrays["rhj_6983_axis_real"]["mismatched"][sk] = d_real_s
            pred_arms_to_save[f"Occam_rhj_6983_axis_reconstructed_mis_{sk}"] = _pack_delta(
                d_rec_s, ctrl_source=sk, ctrl_bcs=encoded_cache[ck_s]["barcodes"],
                axis_mode="rhj_6983_axis", base_mode="reconstructed",
                n_act=encoded_cache[ck_s]["n_active_genes"], mask_sha=encoded_cache[ck_s]["mask_sha"],
            )
            pred_arms_to_save[f"Occam_rhj_6983_axis_real_mis_{sk}"] = _pack_delta(
                d_real_s, ctrl_source=sk, ctrl_bcs=encoded_cache[ck_s]["barcodes"],
                axis_mode="rhj_6983_axis", base_mode="real",
                n_act=encoded_cache[ck_s]["n_active_genes"], mask_sha=encoded_cache[ck_s]["mask_sha"],
            )
        print(f"[{_now_iso()}] [{ln}] Axis 1 (rhj_6983_axis) done (cpu_t={_cpu_time_s():.1f}s).", flush=True)

        # Axis 2: P4 pairwise_common_axis (reusing bit-identical primary predictions whenever zero genes zeroed!)
        for sk in all_sources:
            a, b = (ln, sk) if ln < sk else (sk, ln)
            pmask = pair_masks[(a, b)]
            if not pmask.requires_common_axis_arm:
                continue
            ck_p4_m = f"pairwise_common_axis__{ln}__pair_{a}__{b}__primary"
            ck_p4_s = f"pairwise_common_axis__{sk}__pair_{a}__{b}__primary"

            if ck_p4_m in encoded_cache:
                d_p4_m, _ = _decode_from_cached_latent(ck_p4_m, pert_res.resolved_ids, out_cols)
                bcs_p4_m = encoded_cache[ck_p4_m]["barcodes"]
            else:
                d_p4_m = d_rec_m.copy()
                bcs_p4_m = encoded_cache[ck_m]["barcodes"]

            if ck_p4_s in encoded_cache:
                d_p4_s, _ = _decode_from_cached_latent(ck_p4_s, pert_res.resolved_ids, out_cols)
                bcs_p4_s = encoded_cache[ck_p4_s]["barcodes"]
            else:
                d_p4_s = occam_eval_arrays["rhj_6983_axis_reconstructed"]["mismatched"][sk].copy()
                bcs_p4_s = encoded_cache[f"rhj_6983_axis__{sk}__primary"]["barcodes"]

            occam_eval_arrays["pairwise_common_axis_reconstructed"][sk] = {
                "matched": d_p4_m,
                "mismatched": d_p4_s,
                "n_active_rhj_genes": int(pmask.n_common),
                "axis_mask_sha256": str(pmask.mask_sha256),
            }
            pred_arms_to_save[f"Occam_pairwise_common_axis_pair_{sk}_reconstructed_matched"] = _pack_delta(
                d_p4_m, ctrl_source=ln, ctrl_bcs=bcs_p4_m,
                axis_mode="pairwise_common_axis", base_mode="reconstructed",
                pair_cp=sk, n_act=int(pmask.n_common), mask_sha=str(pmask.mask_sha256),
            )
            pred_arms_to_save[f"Occam_pairwise_common_axis_pair_{sk}_reconstructed_mis_{sk}"] = _pack_delta(
                d_p4_s, ctrl_source=sk, ctrl_bcs=bcs_p4_s,
                axis_mode="pairwise_common_axis", base_mode="reconstructed",
                pair_cp=ln, n_act=int(pmask.n_common), mask_sha=str(pmask.mask_sha256),
            )
        print(f"[{_now_iso()}] [{ln}] Axis 2 (pairwise_common_axis) done (cpu_t={_cpu_time_s():.1f}s).", flush=True)

        # Axis 3: full_44728_api (secondary descriptive)
        ck_fm = f"full_44728_api__{ln}__primary"
        d_full_m, _ = _decode_from_cached_latent(ck_fm, pert_res.resolved_ids, out_cols)
        occam_eval_arrays["full_44728_api_reconstructed"]["matched"] = d_full_m
        pred_arms_to_save["Occam_full_44728_api_reconstructed_matched"] = _pack_delta(
            d_full_m, ctrl_source=ln, ctrl_bcs=encoded_cache[ck_fm]["barcodes"],
            axis_mode="full_44728_api", base_mode="reconstructed",
            n_act=encoded_cache[ck_fm]["n_active_genes"], mask_sha="",
        )

        ck_fr = f"full_44728_api__{ln}__resample"
        d_full_r, _ = _decode_from_cached_latent(ck_fr, pert_res.resolved_ids, out_cols)
        occam_eval_arrays["full_44728_api_reconstructed"]["resample"] = d_full_r
        pred_arms_to_save["Occam_full_44728_api_reconstructed_resample"] = _pack_delta(
            d_full_r, ctrl_source=resamp_src_lbl, ctrl_bcs=encoded_cache[ck_fr]["barcodes"],
            axis_mode="full_44728_api", base_mode="reconstructed",
            n_act=encoded_cache[ck_fr]["n_active_genes"], mask_sha="",
        )

        for sk in all_sources:
            ck_fs = f"full_44728_api__{sk}__primary"
            d_full_s, _ = _decode_from_cached_latent(ck_fs, pert_res.resolved_ids, out_cols)
            occam_eval_arrays["full_44728_api_reconstructed"]["mismatched"][sk] = d_full_s
            pred_arms_to_save[f"Occam_full_44728_api_reconstructed_mis_{sk}"] = _pack_delta(
                d_full_s, ctrl_source=sk, ctrl_bcs=encoded_cache[ck_fs]["barcodes"],
                axis_mode="full_44728_api", base_mode="reconstructed",
                n_act=encoded_cache[ck_fs]["n_active_genes"], mask_sha="",
            )
        print(f"[{_now_iso()}] [{ln}] Axis 3 (full_44728_api) done (cpu_t={_cpu_time_s():.1f}s).", flush=True)

        # Save arrays/<ln>/observed.npz, optional SENS_FCT copy, preds_*.npz, & array_manifest_sha256.json
        out_ldir = arrays_dir / ln
        out_ldir.mkdir(parents=True, exist_ok=True)
        obs_out_path = out_ldir / "observed.npz"
        shutil.copy2(r16_obs_path, obs_out_path)

        extra_fct_sens: Dict[str, Dict[str, np.ndarray]] = {}
        extra_manifest: Dict[str, str] = {"observed.npz": sha256_file(obs_out_path)}
        for sf in sorted(r16_ldir.glob("preds_SENS_FCT_*.npz")):
            dst_sf = out_ldir / sf.name
            shutil.copy2(sf, dst_sf)
            extra_manifest[sf.name] = sha256_file(dst_sf)
            skey = sf.stem[len("preds_SENS_") :]
            with np.load(dst_sf, allow_pickle=False) as zs:
                extra_fct_sens[skey] = {
                    "delta_pred": zs["delta_pred"].copy(),
                    "delta_transfer": zs["delta_transfer"].copy(),
                    "source_grand_mean": zs["source_grand_mean"].copy(),
                }

        pred_manifest = save_occam_line_deltas(arrays_dir, ln, pred_arms_to_save)
        full_manifest = {**extra_manifest, **pred_manifest}
        (out_ldir / "array_manifest_sha256.json").write_text(
            json.dumps(full_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )

        # Evaluate line strictly from arrays + frozen R16 STATE_SE arrays
        se_arrays, se_shas = _load_r16_state_se_arrays_for_line(r16_ldir, all_sources)
        r16_se_shas_all[ln] = se_shas

        pert_audit_dict = {
            "n_tier1": len(pert_res.tier1_targets),
            "n_tier2": len(pert_res.tier2_targets),
            "tier2_targets": list(pert_res.tier2_targets),
            "bridged_pairs": dict(pert_res.bridged_pairs),
            "unresolved_targets": list(pert_res.unresolved_targets),
        }
        line_res = evaluate_addb_line_from_arrays(
            line_name=ln,
            observed=obs_dict,
            occam_arrays=occam_eval_arrays,
            state_se_arrays=se_arrays,
            pert_resolution_audit=pert_audit_dict,
            output_gene_projection_audit=out_gene_audit,
            extra_fct_sens=extra_fct_sens if extra_fct_sens else None,
            n_boot=int(n_boot),
            seed=int(seed),
            delta_tost=float(delta_tost),
        )
        line_res["array_manifest_sha256"] = full_manifest
        lines_results[ln] = line_res
        print(
            f"[{_now_iso()}] [{ln}] DONE in {time.time() - t_ln:.1f}s | "
            f"H-I1_Occam={line_res['primary']['H_I1_Occam']['line_summary']} | "
            f"H-I2_Occam={line_res['primary']['H_I2_Occam']['classification']} | "
            f"H-I3a={line_res['paired_vs_STATE_SE']['H_I3a']['r_orth_cal']['classification']}",
            flush=True,
        )

    summary = {
        "protocol": "PREREG_I_ADDENDUM_B_v2",
        "model_key": OCCAM_MODEL_KEY,
        "checkpoint_id": OCCAM_CHECKPOINT_ID,
        "ckpt_sha256": LOCKED_BEST_MODEL_SHA256,
        "device": LOCKED_DEVICE,
        "torch_num_threads": int(num_threads),
        "batch_size": int(batch_size),
        "seed": int(seed),
        "resample_seed": int(DEFAULT_RESAMPLE_SEED),
        "n_boot": int(n_boot),
        "delta_tost": float(delta_tost),
        "elapsed_sec": round(time.time() - t0, 2),
        "gates_verified": gates_verified,
        "control_barcodes_file": "control_barcodes_addb.json",
        "z_ctrl_dir": "z_ctrl",
        "z_ctrl_manifest_sha256": z_ctrl_manifest,
        "r16_state_se_reference_sha256": r16_se_shas_all,
        "group_adjudications": build_group_adjudications_addb(lines_results),
        "lines": lines_results,
    }
    (out_dir / "summary_addb.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"[{_now_iso()}] Addendum B ALL DONE in {summary['elapsed_sec']}s.", flush=True)
    return summary


# ==============================================================================
# Pure-CPU Zero-Real-Cell Synthetic Smoke + Stage 4 From-Disk Recomputation
# ==============================================================================


def _max_numeric_diff_addb(a: Any, b: Any) -> float:
    skip_keys = {
        "pert_resolution_audit",
        "output_gene_projection_audit",
        "array_manifest_sha256",
        "z_ctrl_manifest_sha256",
        "gates_verified",
        "elapsed_sec",
        "mode",
        "real_cells_touched",
        "stage4_from_disk_recompute_max_abs_diff",
        "stage4_from_disk_recompute_passed",
    }
    if isinstance(a, dict) and isinstance(b, dict):
        mx = 0.0
        for k in set(a.keys()) & set(b.keys()):
            if k in skip_keys:
                continue
            d = _max_numeric_diff_addb(a[k], b[k])
            if np.isnan(d):
                return float("nan")
            mx = max(mx, d)
        return mx
    if isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        mx = 0.0
        for x, y in zip(a, b):
            d = _max_numeric_diff_addb(x, y)
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


def run_synthetic_smoke_addb(
    out_dir: Path,
    *,
    n_boot: int = 400,
    seed: int = DEFAULT_R16_SEED,
    delta_tost: float = DEFAULT_TOST_DELTA,
) -> Dict[str, Any]:
    """
    Pure-CPU, ZERO-real-cell end-to-end synthetic smoke test exercising EVERY Addendum B
    component across all 5 line roles (`a549_ko` w/ Tier-2 `RBM14`, `dld1`, `cd4t_rest`,
    `gse295214_h1` w/ `SENS_FCT`, `feng_ipsc` `sensitivity_only`), all 3 axes (`rhj_6983_axis`,
    `pairwise_common_axis`, `full_44728_api`), both baselines, `z_ctrl` exports, `control_barcodes_addb.json`,
    `H-I1_Occam`, `H-I2_Occam`, `H-I3a`, `H-I3b`, AND Stage 4 from-disk recomputation parity (`<= 1e-12`).
    """
    t0 = time.time()
    out_path = Path(out_dir)
    arrays_dir = out_path / "arrays"
    z_ctrl_dir = out_path / "z_ctrl"
    syn_r16_dir = out_path / "synthetic_r16_state_se"
    arrays_dir.mkdir(parents=True, exist_ok=True)
    z_ctrl_dir.mkdir(parents=True, exist_ok=True)
    syn_r16_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(int(seed))
    line_configs = [
        ("a549_ko", 56, 120, "occam_specific", ["RBM14"], False),
        ("dld1", 52, 120, "occam_equiv_nonzero", [], False),
        ("cd4t_rest", 54, 120, "occam_sub_delta", [], False),
        ("gse295214_h1", 56, 120, "occam_mixed", [], True),
        ("feng_ipsc", 50, 110, "occam_equiv_null", [], False),
    ]

    # Export synthetic z_ctrl latents for all 7 pools (primary) + 5 targets (resample)
    z_entries: Dict[str, Dict[str, Any]] = {}
    bc_pools: Dict[str, Dict[str, Any]] = {}
    all_pool_names = ("a549_ko", "cd4t_rest", "dld1", "feng_ipsc", "gse295214_h1", "nadig_hepg2", "nadig_jurkat")
    for pn in all_pool_names:
        bcs_p = [f"SYN_{pn.upper()}_P_{i:04d}" for i in range(512)]
        bcs_r = [f"SYN_{pn.upper()}_R_{i:04d}" for i in range(512)]
        z_entries[f"rhj_6983_axis__{pn}__primary"] = {
            "Z_ctrl": rng.normal(0.0, 0.4, size=(512, OCCAM_LATENT_DIM)).astype(np.float32),
            "control_barcodes": bcs_p,
            "seed": int(seed),
            "pool_name": pn,
            "input_axis_mode": "rhj_6983_axis",
            "n_active_rhj_genes": 6570 if pn == "cd4t_rest" else (6950 if pn == "a549_ko" else 6983),
            "axis_mask_sha256": "0" * 64,
            "audit": {"synthetic": True},
        }
        has_r = pn in ALL_EVALUATED_LINES_ADDB
        if has_r:
            z_entries[f"rhj_6983_axis__{pn}__resample"] = {
                "Z_ctrl": rng.normal(0.0, 0.4, size=(512, OCCAM_LATENT_DIM)).astype(np.float32),
                "control_barcodes": bcs_r,
                "seed": int(DEFAULT_RESAMPLE_SEED),
                "pool_name": pn,
                "input_axis_mode": "rhj_6983_axis",
                "n_active_rhj_genes": 6570 if pn == "cd4t_rest" else (6950 if pn == "a549_ko" else 6983),
                "axis_mask_sha256": "0" * 64,
                "audit": {"synthetic": True},
            }
        bc_pools[pn] = {
            "primary_seed": int(seed),
            "n_ctrl_total": 2048,
            "n_ctrl_sampled": 512,
            "barcodes": bcs_p,
            "barcodes_sha256": compute_string_list_sha256(bcs_p),
            "within_line_resample": {
                "eligible": has_r,
                "seed": int(DEFAULT_RESAMPLE_SEED),
                "n_resample": 512 if has_r else 0,
                "overlap_with_primary": 0,
                "regime": "full_512_disjoint" if has_r else "",
                "barcodes": bcs_r if has_r else [],
                "barcodes_sha256": compute_string_list_sha256(bcs_r) if has_r else "",
            },
        }

    z_manifest = save_occam_z_ctrl_latents(z_ctrl_dir, z_entries)
    (out_path / "control_barcodes_addb.json").write_text(
        json.dumps({"provenance": {"primary_seed": int(seed), "resample_seed": int(DEFAULT_RESAMPLE_SEED)}, "pools": bc_pools}, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    lines_results: Dict[str, Any] = {}
    r16_se_dirs: Dict[str, Path] = {}
    t2_by_line: Dict[str, List[str]] = {}

    for li, (ln, n_p, n_g, regime, t2_list, has_sens_fct) in enumerate(line_configs):
        lrng = np.random.default_rng(int(seed) + 100 * (li + 1))
        spec = ADDB_LINE_SPECS[ln]
        all_sources = sorted(set(spec["seen_by_checkpoint"]) | set(spec["unseen_by_checkpoint"]))

        t_names = [f"SYN_T{i:03d}" for i in range(n_p)]
        if t2_list:
            t_names[-1] = t2_list[0]
            t2_by_line[ln] = list(t2_list)
        t_names = sorted(t_names)
        g_names = [f"SYN_G{j:04d}" for j in range(n_g)]
        t_idx = np.full(n_p, -1, dtype=np.int64)
        t_idx[: min(8, n_p)] = np.arange(min(8, n_p), dtype=np.int64)

        tm = lrng.normal(0.0, 0.35, size=(n_p, n_g))
        b_s = lrng.normal(0.0, 0.25, size=n_g)
        true_spec = lrng.normal(0.0, 0.55, size=(n_p, n_g))
        noise_a = lrng.normal(0.0, 0.22, size=(n_p, n_g))
        noise_b = lrng.normal(0.0, 0.22, size=(n_p, n_g))
        delta_ha = tm + true_spec + noise_a
        delta_hb = tm + true_spec + noise_b
        delta_obs = 0.5 * (delta_ha + delta_hb)
        half_valid = np.ones(n_p, dtype=bool)

        fct_pred = tm + 0.35 * true_spec + lrng.normal(0.0, 0.25, size=(n_p, n_g))
        e1_pred = tm + 0.10 * true_spec + lrng.normal(0.0, 0.30, size=(n_p, n_g))

        # Construct synthetic STATE_SE arrays (moderate signal, mostly control-invariant)
        se_shared_noise = lrng.normal(0.0, 0.28, size=(n_p, n_g))
        se_m = tm + 0.20 * true_spec + se_shared_noise
        se_r = se_m + lrng.normal(0.0, 0.005, size=(n_p, n_g))
        se_mis = {sk: se_m + lrng.normal(0.0, 0.015, size=(n_p, n_g)) for sk in all_sources}

        # Construct synthetic Occam arrays per regime
        occ_noise = lrng.normal(0.0, 0.18, size=(n_p, n_g))
        if regime == "occam_specific":
            occ_m = tm + 0.55 * true_spec + occ_noise
            occ_r = occ_m + lrng.normal(0.0, 0.008, size=(n_p, n_g))
            occ_mis = {sk: tm + 0.15 * true_spec + occ_noise + lrng.normal(0.0, 0.08, size=(n_p, n_g)) for sk in all_sources}
        elif regime == "occam_equiv_nonzero":
            occ_m = tm + 0.45 * true_spec + occ_noise
            occ_r = occ_m + lrng.normal(0.0, 0.002, size=(n_p, n_g))
            occ_mis = {sk: occ_m + lrng.normal(0.0, 0.003, size=(n_p, n_g)) for sk in all_sources}
        elif regime == "occam_sub_delta":
            occ_m = tm + 0.45 * true_spec + occ_noise
            occ_r = occ_m + lrng.normal(0.0, 0.002, size=(n_p, n_g))
            occ_mis = {sk: tm + 0.42 * true_spec + occ_noise + lrng.normal(0.0, 0.01, size=(n_p, n_g)) for sk in all_sources}
        elif regime == "occam_mixed":
            occ_m = tm + 0.45 * true_spec + occ_noise
            occ_r = occ_m + lrng.normal(0.0, 0.005, size=(n_p, n_g))
            occ_mis = {}
            for si, sk in enumerate(all_sources):
                gain = 0.15 if si % 2 == 0 else 0.45
                occ_mis[sk] = tm + gain * true_spec + occ_noise + lrng.normal(0.0, 0.01, size=(n_p, n_g))
        else:  # occam_equiv_null
            occ_m = tm + lrng.normal(0.0, 0.30, size=(n_p, n_g))
            occ_r = occ_m.copy()
            occ_mis = {sk: occ_m.copy() for sk in all_sources}  # tests EXACT_INVARIANCE -> LINE_EQUIVALENT_NULL_SIGNAL!

        occ_real_m = occ_m + lrng.normal(0.0, 0.02, size=(n_p, n_g))
        occ_real_r = occ_r + lrng.normal(0.0, 0.02, size=(n_p, n_g))
        occ_real_mis = {sk: occ_mis[sk] + lrng.normal(0.0, 0.02, size=(n_p, n_g)) for sk in all_sources}

        occ_full_m = occ_m + lrng.normal(0.0, 0.015, size=(n_p, n_g))
        occ_full_r = occ_r + lrng.normal(0.0, 0.015, size=(n_p, n_g))
        occ_full_mis = {sk: occ_mis[sk] + lrng.normal(0.0, 0.015, size=(n_p, n_g)) for sk in all_sources}

        obs_dict = {
            "delta_obs": delta_obs,
            "delta_half_a": delta_ha,
            "delta_half_b": delta_hb,
            "half_valid": half_valid,
            "delta_transfer": tm,
            "source_grand_mean": b_s,
            "fct_delta_pred": fct_pred,
            "e1_constmed_delta_pred": e1_pred,
            "target_indices": t_idx,
            "target_names": np.asarray(t_names, dtype=str),
            "gene_names": np.asarray(g_names, dtype=str),
            "mu_target": lrng.uniform(0.0, 3.0, size=n_g),
            "mu_sources": lrng.uniform(0.0, 3.0, size=(6, n_g)),
            "matched_ctrl_barcodes": np.asarray(bc_pools[ln]["barcodes"], dtype=str),
            "resample_ctrl_barcodes": np.asarray(bc_pools[ln]["within_line_resample"]["barcodes"], dtype=str),
            "n_ctrl_total": np.int64(2048),
            "n_ctrl_sampled": np.int64(512),
            "ctrl_sampling_fraction": np.float64(0.25),
        }

        out_ldir = arrays_dir / ln
        out_ldir.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(out_ldir / "observed.npz", **obs_dict)

        # Save synthetic R16 STATE_SE arrays for H-I3a/b & Stage 4 recompute
        se_ldir = syn_r16_dir / ln
        se_ldir.mkdir(parents=True, exist_ok=True)
        r16_se_dirs[ln] = se_ldir
        np.savez_compressed(se_ldir / "preds_STATE_SE_fixed_index_zero_matched.npz", delta_pred=se_m)
        np.savez_compressed(se_ldir / "preds_STATE_SE_fixed_index_zero_resample.npz", delta_pred=se_r)
        for sk in all_sources:
            np.savez_compressed(se_ldir / f"preds_STATE_SE_fixed_index_zero_mis_{sk}.npz", delta_pred=se_mis[sk])

        extra_fct_sens = None
        if has_sens_fct:
            sf_pred = fct_pred + lrng.normal(0.0, 0.01, size=(n_p, n_g))
            sf_tm = tm + lrng.normal(0.0, 0.005, size=(n_p, n_g))
            sf_bs = b_s + lrng.normal(0.0, 0.005, size=n_g)
            np.savez_compressed(
                out_ldir / "preds_SENS_FCT_full_r9_6_incl_vcc2025_h1.npz",
                delta_pred=sf_pred,
                delta_transfer=sf_tm,
                source_grand_mean=sf_bs,
                target_names=np.asarray(t_names, dtype=str),
                gene_names=np.asarray(g_names, dtype=str),
                model_key=np.asarray("FCT_full_r9_6_incl_vcc2025_h1"),
            )
            extra_fct_sens = {
                "FCT_full_r9_6_incl_vcc2025_h1": {
                    "delta_pred": sf_pred,
                    "delta_transfer": sf_tm,
                    "source_grand_mean": sf_bs,
                }
            }

        # Pack & save Occam arms matching ADDB_ARRAY_SCHEMA.md
        arms_save: Dict[str, Dict[str, Any]] = {}
        r_ids = list(range(n_p))
        p4_eval: Dict[str, Dict[str, Any]] = {}

        def _pk(arr: np.ndarray, src: str, bcs: Sequence[str], ax: str, bm: str, cp: str = "", na: int = 6983) -> Dict[str, Any]:
            return build_occam_delta_payload(
                delta_pred=arr, target_names=t_names, gene_names=g_names,
                control_source=src, control_barcodes=bcs, input_axis_mode=ax, baseline_mode=bm,
                resolved_pert_keys=t_names, resolved_pert_ids=r_ids, pair_counterpart=cp,
                n_active_rhj_genes=na, axis_mask_sha256="a" * 64 if ax != "full_44728_api" else "",
            )

        res_lbl = f"{ln}_within_line_resample_seed{DEFAULT_RESAMPLE_SEED}"
        arms_save["Occam_rhj_6983_axis_reconstructed_matched"] = _pk(occ_m, ln, bc_pools[ln]["barcodes"], "rhj_6983_axis", "reconstructed")
        arms_save["Occam_rhj_6983_axis_reconstructed_resample"] = _pk(occ_r, res_lbl, bc_pools[ln]["within_line_resample"]["barcodes"], "rhj_6983_axis", "reconstructed")
        arms_save["Occam_rhj_6983_axis_real_matched"] = _pk(occ_real_m, ln, bc_pools[ln]["barcodes"], "rhj_6983_axis", "real")
        arms_save["Occam_rhj_6983_axis_real_resample"] = _pk(occ_real_r, res_lbl, bc_pools[ln]["within_line_resample"]["barcodes"], "rhj_6983_axis", "real")
        arms_save["Occam_full_44728_api_reconstructed_matched"] = _pk(occ_full_m, ln, bc_pools[ln]["barcodes"], "full_44728_api", "reconstructed", na=15000)
        arms_save["Occam_full_44728_api_reconstructed_resample"] = _pk(occ_full_r, res_lbl, bc_pools[ln]["within_line_resample"]["barcodes"], "full_44728_api", "reconstructed", na=15000)

        for sk in all_sources:
            arms_save[f"Occam_rhj_6983_axis_reconstructed_mis_{sk}"] = _pk(occ_mis[sk], sk, bc_pools[sk]["barcodes"], "rhj_6983_axis", "reconstructed")
            arms_save[f"Occam_rhj_6983_axis_real_mis_{sk}"] = _pk(occ_real_mis[sk], sk, bc_pools[sk]["barcodes"], "rhj_6983_axis", "real")
            arms_save[f"Occam_full_44728_api_reconstructed_mis_{sk}"] = _pk(occ_full_mis[sk], sk, bc_pools[sk]["barcodes"], "full_44728_api", "reconstructed", na=15000)

            # Exercise P4 pairs on pairs involving a549_ko, cd4t_rest, or gse295214_h1
            if any(x in ("a549_ko", "cd4t_rest", "gse295214_h1") for x in (ln, sk)):
                n_c = 6557 if {ln, sk} == {"a549_ko", "cd4t_rest"} else (6570 if "cd4t_rest" in (ln, sk) else 6950)
                p4_m = occ_m + lrng.normal(0.0, 0.005, size=(n_p, n_g))
                p4_s = occ_mis[sk] + lrng.normal(0.0, 0.005, size=(n_p, n_g))
                p4_eval[sk] = {"matched": p4_m, "mismatched": p4_s, "n_active_rhj_genes": n_c, "axis_mask_sha256": "a" * 64}
                arms_save[f"Occam_pairwise_common_axis_pair_{sk}_reconstructed_matched"] = _pk(p4_m, ln, bc_pools[ln]["barcodes"], "pairwise_common_axis", "reconstructed", cp=sk, na=n_c)
                arms_save[f"Occam_pairwise_common_axis_pair_{sk}_reconstructed_mis_{sk}"] = _pk(p4_s, sk, bc_pools[sk]["barcodes"], "pairwise_common_axis", "reconstructed", cp=ln, na=n_c)

        pred_manifest = save_occam_line_deltas(arrays_dir, ln, arms_save)
        obs_sha = sha256_file(out_ldir / "observed.npz")
        (out_ldir / "array_manifest_sha256.json").write_text(
            json.dumps({"observed.npz": obs_sha, **pred_manifest}, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

        lines_results[ln] = evaluate_addb_line_from_arrays(
            line_name=ln,
            observed=obs_dict,
            occam_arrays={
                "rhj_6983_axis_reconstructed": {"matched": occ_m, "resample": occ_r, "mismatched": occ_mis},
                "rhj_6983_axis_real": {"matched": occ_real_m, "resample": occ_real_r, "mismatched": occ_real_mis},
                "pairwise_common_axis_reconstructed": p4_eval,
                "full_44728_api_reconstructed": {"matched": occ_full_m, "resample": occ_full_r, "mismatched": occ_full_mis},
            },
            state_se_arrays={"matched": se_m, "resample": se_r, "mismatched": se_mis},
            pert_resolution_audit={"n_tier1": n_p - len(t2_list), "n_tier2": len(t2_list), "tier2_targets": t2_list},
            extra_fct_sens=extra_fct_sens,
            n_boot=int(n_boot),
            seed=int(seed),
            delta_tost=float(delta_tost),
        )

    summary: Dict[str, Any] = {
        "mode": "synthetic_smoke",
        "real_cells_touched": 0,
        "protocol": "PREREG_I_ADDENDUM_B_v2",
        "model_key": OCCAM_MODEL_KEY,
        "checkpoint_id": OCCAM_CHECKPOINT_ID,
        "device": LOCKED_DEVICE,
        "torch_num_threads": LOCKED_CPU_NUM_THREADS,
        "seed": int(seed),
        "resample_seed": int(DEFAULT_RESAMPLE_SEED),
        "n_boot": int(n_boot),
        "delta_tost": float(delta_tost),
        "z_ctrl_manifest_sha256": z_manifest,
        "group_adjudications": build_group_adjudications_addb(lines_results),
        "lines": lines_results,
    }

    # Stage 4 independent from-disk recomputation check!
    recomp = recompute_addb_from_exported_arrays(
        arrays_dir=arrays_dir,
        r16_state_se_dirs=r16_se_dirs,
        tier2_targets_by_line=t2_by_line,
        n_boot=int(n_boot),
        seed=int(seed),
        delta_tost=float(delta_tost),
    )
    max_diff = _max_numeric_diff_addb(summary, recomp)
    summary["stage4_from_disk_recompute_max_abs_diff"] = float(max_diff)
    summary["stage4_from_disk_recompute_passed"] = bool(np.isfinite(max_diff) and max_diff <= 1e-12)
    summary["elapsed_sec"] = round(time.time() - t0, 2)
    (out_path / "summary_addb.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


# ==============================================================================
# CLI Entrypoint
# ==============================================================================


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Round 16 Addendum B (Occam) Evaluation Runner (CPU float32 locked).")
    ap.add_argument("--synthetic-smoke", action="store_true", help="Run zero-real-cell CPU synthetic smoke + Stage 4 recompute.")
    ap.add_argument("--cohort", choices=["all", "group2", "group3", "sensitivity"], default="all")
    ap.add_argument("--lines", nargs="*", default=None, help="Explicit subset of Addendum B target lines.")
    ap.add_argument("--results-dir", type=Path, default=DEFAULT_ADDB_RESULTS_DIR)
    ap.add_argument("--code-freeze-path", type=Path, default=DEFAULT_CODE_FREEZE_PATH)
    ap.add_argument("--go-path", type=Path, default=DEFAULT_GO_PATH)
    ap.add_argument("--lock-path", type=Path, default=DEFAULT_LOCK_PATH)
    ap.add_argument("--human-approval-path", type=Path, default=DEFAULT_HUMAN_APPROVAL_PATH)
    ap.add_argument("--raw-nadig-dir", type=Path, default=DEFAULT_LOCAL_NADIG_DIR)
    ap.add_argument("--gse295214-raw-dir", type=Path, default=DEFAULT_GSE295214_RAW_DIR)
    ap.add_argument("--device", choices=["cpu"], default=LOCKED_DEVICE)
    ap.add_argument("--num-threads", type=int, default=LOCKED_CPU_NUM_THREADS)
    ap.add_argument("--batch-size", type=int, default=LOCKED_BATCH_SIZE)
    ap.add_argument("--n-boot", type=int, default=LOCKED_N_BOOT)
    ap.add_argument("--seed", type=int, default=DEFAULT_R16_SEED)
    ap.add_argument("--delta-tost", type=float, default=DEFAULT_TOST_DELTA)
    ap.add_argument("--smoke-out-dir", type=Path, default=None)
    return ap.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.synthetic_smoke:
        if args.smoke_out_dir is not None:
            res = run_synthetic_smoke_addb(
                args.smoke_out_dir, n_boot=min(args.n_boot, 400), seed=args.seed, delta_tost=args.delta_tost
            )
        else:
            with tempfile.TemporaryDirectory(prefix="r16_addb_smoke_") as tmpd:
                res = run_synthetic_smoke_addb(
                    Path(tmpd), n_boot=min(args.n_boot, 400), seed=args.seed, delta_tost=args.delta_tost
                )
        compact = {
            "mode": res["mode"],
            "protocol": res["protocol"],
            "real_cells_touched": res["real_cells_touched"],
            "device": res["device"],
            "elapsed_sec": res["elapsed_sec"],
            "lines": {
                ln: {
                    "H_I1_Occam": v["primary"]["H_I1_Occam"]["line_summary"],
                    "H_I2_Occam": v["primary"]["H_I2_Occam"]["classification"],
                    "H_I3a": v["paired_vs_STATE_SE"]["H_I3a"]["r_orth_cal"]["classification"],
                    "advance_case": v["paired_vs_STATE_SE"]["advance_interpretation_case"],
                    "p4_any_changed": v["secondary"]["pairwise_common_axis"]["any_verdict_changed_vs_primary"],
                    "has_tier1_sens": v["secondary"]["tier1_trained_only_sensitivity"] is not None,
                    "has_extra_fct_sens": v["secondary"]["extra_fct_sensitivity"] is not None,
                }
                for ln, v in res["lines"].items()
            },
            "n_z_ctrl_exported": len(res["z_ctrl_manifest_sha256"]),
            "stage4_from_disk_recompute_max_abs_diff": res["stage4_from_disk_recompute_max_abs_diff"],
            "stage4_from_disk_recompute_passed": res["stage4_from_disk_recompute_passed"],
        }
        print(json.dumps(compact, indent=2))
        return 0 if res["stage4_from_disk_recompute_passed"] else 1

    # STRICT EXECUTION GATE BEFORE ANY REAL-TARGET ACCESS
    assert_addb_execution_gates(
        code_freeze_path=args.code_freeze_path,
        go_path=args.go_path,
        lock_path=args.lock_path,
        human_approval_path=args.human_approval_path,
    )

    if args.lines:
        target_lines = list(args.lines)
    elif args.cohort == "group2":
        target_lines = list(GROUP2_LINES_ADDB)
    elif args.cohort == "group3":
        target_lines = list(GROUP3_LINES_ADDB)
    elif args.cohort == "sensitivity":
        target_lines = list(SENSITIVITY_LINES_ADDB)
    else:
        target_lines = list(ALL_EVALUATED_LINES_ADDB)

    run_real_evaluation_addb(
        target_lines=target_lines,
        results_dir=args.results_dir,
        code_freeze_path=args.code_freeze_path,
        go_path=args.go_path,
        lock_path=args.lock_path,
        human_approval_path=args.human_approval_path,
        raw_nadig_dir=args.raw_nadig_dir,
        gse295214_raw_dir=args.gse295214_raw_dir,
        device=args.device,
        num_threads=args.num_threads,
        batch_size=args.batch_size,
        n_boot=args.n_boot,
        seed=args.seed,
        delta_tost=args.delta_tost,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
