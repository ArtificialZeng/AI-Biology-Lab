#!/usr/bin/env python3
"""
Round 17 Line Registry, Execution Gates, Mismatched Control Strata, &
Read-Only Raw-Count NTC + Barcode Extractors (`PREREG_J` v2).

Hard Boundaries & Reproducibility Guarantees:
  1. Execution Gate (`assert_all_execution_gates_r17`): Refuses ANY real-cell data read or
     inference unless BOTH `round17/status/R17_CODE_FREEZE.txt` AND `round17/status/R17_GO.txt`
     exist and are non-empty (plus `round17/prereg/PREREG_J_LOCK.md` when using default repo
     paths), raising frozen `BlindingGateError`.
  2. Strict Zero-Leakage & Omission Enforcement:
       - `wtc11_i3neuron` (`N = 1`) is explicitly OMITTED per Controller Ruling 2.
       - All 6 Round 16 discovery lines (`nadig_hepg2`, `nadig_jurkat`, `a549_ko`, `dld1`,
         `cd4t_rest`, `feng_ipsc`) are strictly quarantined as mismatched control pools ONLY
         and are NEVER allowed as Round 17 target lines.
  3. Read-Only Raw-Count NTC Extractors:
       - Imports frozen `round16/src/run_eval_r16.py` & `round16/src/model_adapters_r16.py`
         for stratified control sampling (`sample_stratified_control_sets_r16`,
         `sample_disjoint_resample_control_sets_r16`), `gse295214_h1`, and all 6 R16
         mismatched control pools (`nadig_hepg2`, `nadig_jurkat`, `a549_ko`, `dld1`,
         `cd4t_rest`, `feng_ipsc`).
       - Mirrors frozen Round 10–14 readers strictly in read-only mode (ZERO writes under
         `round9/`–`round16/` or locked `round17/` files) for `cd4t_stim48`, `a549_myllia`,
         `hek293_xucao`, `telohaec`, `calu3`, `mac_ctrl`, `mac_lps`, `a549_ko_cas9pos`, and
         reuses locked `pfizer_telohaec_reader_r17.py` helpers for `pfizer_telohaec_basal`
         and `pfizer_telohaec_cytokine`.
"""

from __future__ import annotations

import csv
import re
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROUND17_SRC_DIR = Path(__file__).resolve().parent
ROUND17_ROOT = ROUND17_SRC_DIR.parent
REPO_ROOT = ROUND17_ROOT.parent
R16_SRC_DIR = REPO_ROOT / "round16" / "src"
R14_SRC_DIR = REPO_ROOT / "round14" / "src"
R13_SRC_DIR = REPO_ROOT / "round13" / "src"
R12_SRC_DIR = REPO_ROOT / "round12" / "src"
R11_SRC_DIR = REPO_ROOT / "round11" / "src"
R11_SCRIPTS_DIR = REPO_ROOT / "round11" / "shared" / "scripts"
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
R10_SCRIPTS_DIR = REPO_ROOT / "round10" / "shared" / "scripts"
for _p in (
    str(ROUND17_SRC_DIR),
    str(R16_SRC_DIR),
    str(R14_SRC_DIR),
    str(R13_SRC_DIR),
    str(R12_SRC_DIR),
    str(R11_SRC_DIR),
    str(R11_SCRIPTS_DIR),
    str(R10_SRC_DIR),
    str(R10_SCRIPTS_DIR),
):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from data_loader_r10 import (  # noqa: E402
    DEFAULT_HGNC_FILE,
    HGNCBridge,
    R9_SOURCE_LINE_ORDER,
    R10_LINES_DIR,
    decode_str_array,
)
from model_adapters_r16 import (  # noqa: E402
    DEFAULT_MODELS_DIR,
    GROUP3_PRIMARY_FCT_SOURCES,
    BlindingGateError,
    assert_zero_training_leakage,
    sample_disjoint_resample_control_sets_r16,
    sample_stratified_control_sets_r16,
)
from pfizer_telohaec_reader_r17 import (  # noqa: E402
    PFIZER_TELOHAEC_FILENAME,
    _assign_hardware_split_halves_r17,
    _compute_rna_umi_and_ngenes_r17,
    _extract_var_rna_features_r17,
    _make_h5ad_subset_readers_r17,
    _pick_obs_col_r17,
    _resolve_h5ad_path_r17,
    classify_pfizer_condition_r17,
    classify_pfizer_guide_r17,
    read_h5ad_table_col_r17,
)
from run_eval_r16 import (  # noqa: E402
    _A549_NTC_TOKENS,
    _A549_UNASSIGNED_TOKENS,
    _CD4_EXCEL_GUIDE_REMAP,
    _MULTI_GUIDE_DELIM_RE,
    LightweightPartSpec,
    _apply_min_ntc_batch_fallback_r16,
    _build_cd4t_rest_ntc_lightweight_parts_r16,
    _build_dld1_ntc_lightweight_parts_r16,
    _extract_a549_ntc_bundle_r16,
    _extract_feng_ipsc_ntc_bundle_r16,
    _extract_top2_csc_r16,
    _find_file_r16,
    _fit_poisson_gaussian_thresholds_r16,
    _make_10x_csc_subset_readers_r16,
    _pick_col_r16,
    _read_mtx_header_r16,
    _resolve_nadig_h5ad_path_r16,
    build_gse295214_parts,
    extract_nadig_h5ad_ntc_bundle_r16,
    extract_ntc_bundle_from_parts_r16,
)
from stats_r17 import DEFAULT_R17_RESAMPLE_SEED, DEFAULT_R17_SEED  # noqa: E402

DEFAULT_R17_STATUS_DIR = ROUND17_ROOT / "status"
DEFAULT_R17_PREREG_LOCK_PATH = ROUND17_ROOT / "prereg" / "PREREG_J_LOCK.md"
DEFAULT_R17_SHARED_DIR = ROUND17_ROOT / "shared"
DEFAULT_R17_RESULTS_DIR = ROUND17_ROOT / "results"
DEFAULT_PFIZER_RAW_PATH = Path("/Users/mac/.cache/r17_telohaec_pfizer") / PFIZER_TELOHAEC_FILENAME
DEFAULT_GSE295214_RAW_DIR = Path("/Users/mac/.cache/r16_gse295214/raw")

R16_DISCOVERY_QUARANTINED_LINES: Tuple[str, ...] = (
    "nadig_hepg2",
    "nadig_jurkat",
    "a549_ko",
    "dld1",
    "cd4t_rest",
    "feng_ipsc",
)
OMITTED_LINES_R17: Tuple[str, ...] = ("wtc11_i3neuron",)

CONFIRMATORY_LINES_R17: Tuple[str, ...] = (
    "cd4t_stim48",
    "gse295214_h1",
    "pfizer_telohaec_basal",
)
WITHIN_FAMILY_LINES_R17: Tuple[str, ...] = ("pfizer_telohaec_cytokine",)
SECONDARY_FRESH_LINES_R17: Tuple[str, ...] = ("a549_myllia", "hek293_xucao")
MICRO_N_LINES_R17: Tuple[str, ...] = ("telohaec", "calu3", "mac_lps", "mac_ctrl")
SUBSET_SENSITIVITY_LINES_R17: Tuple[str, ...] = ("a549_ko_cas9pos",)

ALL_EVALUATED_LINES_R17: Tuple[str, ...] = (
    CONFIRMATORY_LINES_R17
    + WITHIN_FAMILY_LINES_R17
    + SECONDARY_FRESH_LINES_R17
    + MICRO_N_LINES_R17
    + SUBSET_SENSITIVITY_LINES_R17
)

SEEN_MISMATCHED_POOL_R17: Tuple[str, ...] = ("nadig_hepg2", "nadig_jurkat")
UNSEEN_MISMATCHED_CANDIDATES_R17: Tuple[str, ...] = ("a549_ko", "cd4t_rest", "dld1", "feng_ipsc")

LINE_SPECS_R17: Dict[str, Dict[str, Any]] = {
    "cd4t_stim48": {
        "role": "confirmatory",
        "modality": "CRISPRi",
        "checkpoint_id": "zeroshot/ensemble4",
        "fct_sources": R9_SOURCE_LINE_ORDER,
        "line_dir": REPO_ROOT / "round11" / "shared" / "lines" / "cd4t_stim48",
        "lineage_group": "primary_cd4_t",
        "disclosures": [
            "Disjoint donors (D2/D3 vs D1/D4), disjoint lanes (R2L25-48), and 48h stimulated vs resting state relative to R16 discovery cd4t_rest, but shares primary CD4+ T lineage, sgRNA library, and GSE314342 study.",
        ],
    },
    "gse295214_h1": {
        "role": "confirmatory",
        "modality": "CRISPRi",
        "checkpoint_id": "zeroshot/ensemble4",
        "fct_sources": GROUP3_PRIMARY_FCT_SOURCES,
        "fct_sources_sensitivity": R9_SOURCE_LINE_ORDER,
        "line_dir": REPO_ROOT / "round16" / "shared" / "lines" / "gse295214_h1",
        "lineage_group": "h1_hesc",
        "disclosures": [
            "Strictly prospective relative to PREREG_J lock (zero prior STATE inference before lock).",
            "Shares H1 hESC lineage with Round 9 vcc2025_h1 (excluded from primary FCT pool; full r9_6 reported as fixed-N sensitivity) and Replogle/Weissman laboratory provenance.",
        ],
    },
    "pfizer_telohaec_basal": {
        "role": "confirmatory",
        "modality": "CRISPRi",
        "checkpoint_id": "zeroshot/ensemble4",
        "fct_sources": R9_SOURCE_LINE_ORDER,
        "line_dir": DEFAULT_R17_SHARED_DIR / "lines" / "pfizer_telohaec_basal",
        "lineage_group": "telohaec",
        "disclosures": [
            "Zenodo 18792213 (imru_full.h5ad, CC BY 4.0, 2026-03-04, Wong et al.); Untreated basal condition is sole primary confirmatory unit for this study.",
            "Shares TeloHAEC lineage with local micro-N line telohaec (GSE210681, N=9, never run through STATE prior to R17).",
        ],
    },
    "pfizer_telohaec_cytokine": {
        "role": "within_family_replication",
        "modality": "CRISPRi",
        "checkpoint_id": "zeroshot/ensemble4",
        "fct_sources": R9_SOURCE_LINE_ORDER,
        "line_dir": DEFAULT_R17_SHARED_DIR / "lines" / "pfizer_telohaec_cytokine",
        "lineage_group": "telohaec",
        "disclosures": [
            "Treated (IL-1b + TNF-a) condition from Zenodo 18792213 (alias pfizer_telohaec_treated in PREREG_J prose); strictly secondary within-family replication, NEVER counted in K_conf or k.",
        ],
    },
    "a549_myllia": {
        "role": "secondary_fresh",
        "modality": "CRISPRko",
        "checkpoint_id": "zeroshot/ensemble4",
        "fct_sources": R9_SOURCE_LINE_ORDER,
        "line_dir": REPO_ROOT / "round14" / "shared" / "lines" / "a549_myllia",
        "lineage_group": "a549",
        "disclosures": [
            "Secondary descriptive fresh line (N=20 < 50); shares A549 cell line with R16 a549_ko (excluded from G4 unseen mismatched pool per PREREG_J §3.5).",
        ],
    },
    "hek293_xucao": {
        "role": "secondary_fresh",
        "modality": "CRISPRi",
        "checkpoint_id": "zeroshot/ensemble4",
        "fct_sources": R9_SOURCE_LINE_ORDER,
        "line_dir": REPO_ROOT / "round12" / "shared" / "lines" / "hek293_xucao",
        "lineage_group": "hek293",
        "disclosures": [
            "Secondary descriptive fresh line (N=16 < 50); predates SE-600M (2023).",
        ],
    },
    "telohaec": {
        "role": "micro_n_descriptive",
        "modality": "CRISPRi",
        "checkpoint_id": "zeroshot/ensemble4",
        "fct_sources": R9_SOURCE_LINE_ORDER,
        "line_dir": REPO_ROOT / "round11" / "shared" / "lines" / "telohaec",
        "lineage_group": "telohaec",
        "disclosures": [
            "Descriptive micro-N line (N=9); shares TeloHAEC lineage with Pfizer Unit 3; predates SE-600M (2023-24).",
        ],
    },
    "calu3": {
        "role": "micro_n_descriptive",
        "modality": "CRISPRi",
        "checkpoint_id": "zeroshot/ensemble4",
        "fct_sources": R9_SOURCE_LINE_ORDER,
        "line_dir": R10_LINES_DIR / "calu3",
        "lineage_group": "calu3",
        "disclosures": ["Descriptive micro-N line (N=8); predates SE-600M (2023)."],
    },
    "mac_lps": {
        "role": "micro_n_descriptive",
        "modality": "CRISPRko",
        "checkpoint_id": "zeroshot/ensemble4",
        "fct_sources": R9_SOURCE_LINE_ORDER,
        "line_dir": REPO_ROOT / "round13" / "shared" / "lines" / "mac_lps",
        "lineage_group": "primary_macrophage_gse327644",
        "disclosures": ["Descriptive micro-N line (N=4); paired condition with mac_ctrl (GSE327644)."],
    },
    "mac_ctrl": {
        "role": "micro_n_descriptive",
        "modality": "CRISPRko",
        "checkpoint_id": "zeroshot/ensemble4",
        "fct_sources": R9_SOURCE_LINE_ORDER,
        "line_dir": REPO_ROOT / "round13" / "shared" / "lines" / "mac_ctrl",
        "lineage_group": "primary_macrophage_gse327644",
        "disclosures": ["Descriptive micro-N line (N=3); paired condition with mac_lps (GSE327644)."],
    },
    "a549_ko_cas9pos": {
        "role": "r16_subset_sensitivity",
        "modality": "CRISPRko",
        "checkpoint_id": "zeroshot/ensemble4",
        "fct_sources": R9_SOURCE_LINE_ORDER,
        "line_dir": REPO_ROOT / "round12" / "shared" / "lines" / "a549_ko_cas9pos",
        "lineage_group": "a549",
        "disclosures": [
            "R16 subset sensitivity ONLY (strict 73.3% CAS9>0 cell subset of R16 discovery a549_ko); strictly excluded from confirmatory family and K_conf.",
        ],
    },
}

_POOL_LINEAGE_GROUP: Dict[str, str] = {
    "nadig_hepg2": "hepg2",
    "nadig_jurkat": "jurkat",
    "a549_ko": "a549",
    "cd4t_rest": "primary_cd4_t",
    "dld1": "dld1",
    "feng_ipsc": "feng_ipsc",
}


# ==============================================================================
# Execution Gates & Leakage/Lineage Guards
# ==============================================================================


def assert_all_execution_gates_r17(
    status_dir: Path = DEFAULT_R17_STATUS_DIR,
    prereg_lock_path: Optional[Path] = DEFAULT_R17_PREREG_LOCK_PATH,
) -> Dict[str, Path]:
    """
    Hard execution gate refusing ANY real-cell data access or model inference unless:
      - `<status_dir>/R17_CODE_FREEZE.txt` exists and is non-empty;
      - `<status_dir>/R17_GO.txt` exists and is non-empty;
      - `prereg_lock_path` (`round17/prereg/PREREG_J_LOCK.md` by default, if provided) exists
        and is non-empty.
    Raises frozen `BlindingGateError` on any missing/empty gate.
    """
    sdir = Path(status_dir)
    freeze_path = sdir / "R17_CODE_FREEZE.txt"
    go_path = sdir / "R17_GO.txt"

    if prereg_lock_path is not None:
        lp = Path(prereg_lock_path)
        if not lp.is_file() or lp.stat().st_size == 0:
            raise BlindingGateError(
                f"Refusing Round 17 real-cell execution: missing or empty PREREG_J lock ({lp})."
            )
    if not freeze_path.is_file() or freeze_path.stat().st_size == 0:
        raise BlindingGateError(
            f"Refusing Round 17 real-cell execution: missing or empty freeze gate ({freeze_path}). "
            "Master Controller must write R17_CODE_FREEZE.txt first."
        )
    if not go_path.is_file() or go_path.stat().st_size == 0:
        raise BlindingGateError(
            f"Refusing Round 17 real-cell execution: missing or empty GO gate ({go_path}). "
            "Master Controller must authorize execution by writing R17_GO.txt."
        )
    out = {"code_freeze": freeze_path, "r17_go": go_path}
    if prereg_lock_path is not None:
        out["prereg_lock"] = Path(prereg_lock_path)
    return out


def validate_target_lines_r17(lines: Sequence[str]) -> List[str]:
    """Validate target line list, refusing omitted (`wtc11_i3neuron`) and quarantined R16 lines."""
    out: List[str] = []
    for ln in lines:
        if ln in OMITTED_LINES_R17:
            raise ValueError(
                f"Line {ln!r} is explicitly OMITTED from Round 17 inference per Controller Ruling 2 (N=1)."
            )
        if ln in R16_DISCOVERY_QUARANTINED_LINES:
            raise ValueError(
                f"Line {ln!r} is a quarantined Round 16 discovery line and cannot run as a Round 17 target."
            )
        if ln not in LINE_SPECS_R17:
            raise ValueError(
                f"Unknown Round 17 line {ln!r}; supported lines={list(ALL_EVALUATED_LINES_R17)}."
            )
        spec = LINE_SPECS_R17[ln]
        assert_zero_training_leakage(
            "STATE_SE",
            ln,
            checkpoint_id=spec["checkpoint_id"],
            source_pool_lines=spec["fct_sources"],
        )
        out.append(ln)
    return out


def build_mismatched_strata_r17(target_line: str) -> Dict[str, Tuple[str, ...]]:
    """
    Build exact `PREREG_J §3.5` `J5a` mismatched-control pools for `target_line`:
      - `seen_by_checkpoint`: `('nadig_hepg2', 'nadig_jurkat')`
      - `unseen_by_checkpoint`: `('a549_ko', 'cd4t_rest', 'dld1', 'feng_ipsc')` excluding any
        pool sharing canonical lineage with `target_line` (`a549_ko` excluded on `a549_myllia`
        and `a549_ko_cas9pos`; `cd4t_rest` excluded from foreign unseen pool on `cd4t_stim48`).
      - `same_lineage_comparators`: `('cd4t_rest',)` on `cd4t_stim48` (`same_lineage_resting_comparator`),
        empty tuple elsewhere.
    """
    if target_line not in LINE_SPECS_R17:
        raise ValueError(f"Unknown Round 17 target line: {target_line!r}")
    t_lin = LINE_SPECS_R17[target_line]["lineage_group"]

    seen = tuple(sorted(s for s in SEEN_MISMATCHED_POOL_R17 if _POOL_LINEAGE_GROUP[s] != t_lin))
    unseen = tuple(sorted(s for s in UNSEEN_MISMATCHED_CANDIDATES_R17 if _POOL_LINEAGE_GROUP[s] != t_lin))
    same_lin: Tuple[str, ...] = ("cd4t_rest",) if target_line == "cd4t_stim48" else ()
    return {
        "seen_by_checkpoint": seen,
        "unseen_by_checkpoint": unseen,
        "same_lineage_comparators": same_lin,
    }


def canonicalize_gene_r17(bridge: HGNCBridge, raw_gene: str) -> str:
    """
    Canonicalize a raw gene identifier via frozen `HGNCBridge`, routing `ENSG...` IDs through
    `ensembl_id` so Ensembl-indexed datasets (such as Pfizer TeloHAEC) map cleanly to approved symbols.
    """
    s = str(raw_gene).strip()
    if not s:
        return s
    ens_arg = s.split(".")[0] if s.upper().startswith("ENSG") else None
    return bridge.canonicalize(s, ensembl_id=ens_arg) or s


# ==============================================================================
# Read-Only Raw-Count NTC + Barcode Extractors for Fresh R17 Lines
# ==============================================================================


def _build_cd4t_stim48_ntc_lightweight_parts_r17(raw_dir: Path) -> List[LightweightPartSpec]:
    """
    Read-only exact mirror of `dataset_readers_r11.build_cd4t_stim48_parts` strictly on
    QC-passing NTC cells (`MT <= 5.0%`, `ngenes >= 500`, per-donor Poisson-Gaussian EM,
    strict `NTC-` guides, Half A = `CE0008678_D3_Stim48hr`, Half B = `CE0010866_D2_Stim48hr`),
    preserving exact `f"{donor_id}_{lane_tag}:{raw_bc}"` barcodes with ZERO writes.
    """
    cdir = Path(raw_dir) / "cd4t_stim48" if (Path(raw_dir) / "cd4t_stim48").is_dir() else Path(raw_dir)
    csv_name = "sgrna_library_metadata.suppl_table.csv"
    sg_csv = (cdir / csv_name) if (cdir / csv_name).exists() else (REPO_ROOT / "round10" / "shared" / "raw" / "cd4t_rest" / csv_name)
    if not sg_csv.exists():
        raise FileNotFoundError(f"Missing {csv_name} for cd4t_stim48.")

    ntc_guide_set: Set[str] = set()
    with open(sg_csv, "r", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh):
            g = str(row.get("sgRNA", "")).strip()
            if g.startswith("NTC-"):
                ntc_guide_set.add(g)

    d3_files = sorted(cdir.glob("GSM93938*_D3_Stim48hr_*.h5"))
    d2_files = sorted(cdir.glob("GSM93938*_D2_Stim48hr_*.h5"))
    if not d3_files or not d2_files:
        raise FileNotFoundError(f"Missing cd4t_stim48 D3/D2 .h5 files in {cdir}")

    parts: List[LightweightPartSpec] = []
    for donor_id, half_val, files in [
        ("CE0008678_D3_Stim48hr", 1, d3_files),
        ("CE0010866_D2_Stim48hr", 2, d2_files),
    ]:
        lane_cache = []
        hist: Optional[np.ndarray] = None
        crp_names_ref: Optional[List[str]] = None
        for lp in files:
            parts_name = lp.name.split("_")
            lane_tag = f"{parts_name[1]}_{parts_name[2]}" if len(parts_name) >= 3 else lp.stem
            with h5py.File(lp, "r") as f:
                m = f["matrix"]
                n_feat, n_col = int(m["shape"][0]), int(m["shape"][1])
                raw_bcs = decode_str_array(m["barcodes"][:])
                ftypes = np.array(decode_str_array(m["features/feature_type"][:]))
                fids = np.array(decode_str_array(m["features/id"][:]))
                fnames = np.array(decode_str_array(m["features/name"][:]))
                rna_rows = np.where(ftypes == "Gene Expression")[0]
                crp_rows = np.where(ftypes != "Gene Expression")[0]
                if crp_names_ref is None:
                    crp_names_ref = [_CD4_EXCEL_GUIDE_REMAP.get(x, x) for x in fnames[crp_rows]]
                    hist = np.zeros((len(crp_rows), 66), dtype=np.int64)
                mt_raw = np.zeros(n_feat, dtype=bool)
                mt_raw[rna_rows] = [str(x).startswith("MT-") for x in fnames[rna_rows]]

                is_crp = np.zeros(n_feat, dtype=bool)
                is_crp[crp_rows] = True
                crp_remap_idx = np.full(n_feat, -1, dtype=np.int32)
                crp_remap_idx[crp_rows] = np.arange(len(crp_rows), dtype=np.int32)
                is_rna_feat = ~is_crp

                umi_l = np.zeros(n_col, dtype=np.float64)
                ng_l = np.zeros(n_col, dtype=np.int32)
                mt_umi_l = np.zeros(n_col, dtype=np.float64)
                crp_idx_parts, crp_dat_parts, crp_col_parts = [], [], []

                indptr_all = m["indptr"][:]
                data_ds = m["data"]
                indices_ds = m["indices"]
                for cs in range(0, n_col, 25000):
                    ce = min(cs + 25000, n_col)
                    p0, p1 = int(indptr_all[cs]), int(indptr_all[ce])
                    if p1 == p0:
                        continue
                    dat = data_ds[p0:p1]
                    idx = indices_ds[p0:p1]
                    col_lens = np.diff(indptr_all[cs : ce + 1])
                    col_ids = np.repeat(np.arange(cs, ce, dtype=np.int32), col_lens)
                    m_r = is_rna_feat[idx]
                    if np.any(m_r):
                        r_dat = dat[m_r].astype(np.float64, copy=False)
                        r_idx = idx[m_r]
                        r_col = col_ids[m_r] - cs
                        umi_l[cs:ce] = np.bincount(r_col, weights=r_dat, minlength=ce - cs)
                        ng_l[cs:ce] = np.bincount(r_col, minlength=ce - cs)
                        mt_m = mt_raw[r_idx]
                        if np.any(mt_m):
                            mt_umi_l[cs:ce] = np.bincount(r_col[mt_m], weights=r_dat[mt_m], minlength=ce - cs)
                    m_c = ~m_r
                    if np.any(m_c):
                        crp_idx_parts.append(crp_remap_idx[idx[m_c]])
                        crp_dat_parts.append(dat[m_c].astype(np.int64, copy=False))
                        crp_col_parts.append(col_ids[m_c])

            mt_l = 100.0 * mt_umi_l / np.maximum(umi_l, 1.0)
            ok_init = (mt_l < 20.0) & (ng_l >= 200)
            ok_hq = (mt_l <= 5.0) & (ng_l >= 500)
            c_idx = np.concatenate(crp_idx_parts) if crp_idx_parts else np.empty(0, dtype=np.int32)
            c_dat = np.concatenate(crp_dat_parts) if crp_dat_parts else np.empty(0, dtype=np.int64)
            c_col = np.concatenate(crp_col_parts) if crp_col_parts else np.empty(0, dtype=np.int32)
            if c_col.size > 0 and hist is not None:
                ok_init_mask = ok_init[c_col]
                if np.any(ok_init_mask):
                    np.add.at(hist, (c_idx[ok_init_mask], np.clip(c_dat[ok_init_mask], 1, 65)), 1)

            hq_cols = np.where(ok_hq)[0]
            col_map_hq = np.full(n_col, -1, dtype=np.int32)
            col_map_hq[hq_cols] = np.arange(len(hq_cols), dtype=np.int32)
            if c_col.size > 0 and len(hq_cols) > 0:
                hq_m = col_map_hq[c_col] >= 0
                C_hq = sp.csr_matrix(
                    (c_dat[hq_m].astype(np.float64), (c_idx[hq_m], col_map_hq[c_col[hq_m]])),
                    shape=(len(crp_rows), len(hq_cols)),
                ).tocsc()
            else:
                C_hq = sp.csc_matrix((len(crp_rows), len(hq_cols)), dtype=np.float64)
            lane_cache.append(
                (lp, lane_tag, n_feat, rna_rows, [str(fids[i]) for i in rna_rows], [str(fnames[i]) for i in rna_rows], raw_bcs, umi_l, hq_cols, C_hq)
            )

        assert hist is not None and crp_names_ref is not None
        thrs = _fit_poisson_gaussian_thresholds_r16(hist, min_thr=5, max_thr=20, default_thr=7)
        for lp, lane_tag, n_feat, rna_rows, rids, rnames, raw_bcs, umi_l, hq_cols, C_hq in lane_cache:
            t1, t2, a1 = _extract_top2_csc_r16(C_hq)
            pos_mask = C_hq.data >= thrs[C_hq.indices]
            pos_csc = sp.csc_matrix((pos_mask.astype(np.int8), C_hq.indices, C_hq.indptr), shape=C_hq.shape)
            n_pos = np.asarray(pos_csc.sum(axis=0)).ravel()
            s1_mask = (n_pos == 1) & (t1 >= 5.0) & ((t2 < 3.0) | (t1 / np.maximum(t2, 0.5) >= 3.0))

            keep_cols = [
                int(hq_cols[loc_j])
                for loc_j in np.where(s1_mask)[0]
                if crp_names_ref[int(a1[loc_j])] in ntc_guide_set
            ]
            keep_arr = np.asarray(keep_cols, dtype=np.int64)
            bcs_kept = [f"{donor_id}_{lane_tag}:{raw_bcs[int(c)]}" for c in keep_arr]
            read_chunk, read_rows, read_chunk_cols = _make_10x_csc_subset_readers_r16(lp, keep_arr, n_feat, rna_rows)
            parts.append(
                LightweightPartSpec(
                    tag=lane_tag,
                    path=lp,
                    n_cells=len(keep_arr),
                    n_genes_raw=len(rna_rows),
                    nt_mask=np.ones(len(keep_arr), dtype=bool),
                    batches=np.full(len(keep_arr), f"{donor_id}_{lane_tag}"),
                    half_flags=np.full(len(keep_arr), half_val, dtype=np.int8),
                    umi=umi_l[keep_arr],
                    pert_keys=[""] * len(keep_arr),
                    pert_meta={},
                    var_gene_ids=rids,
                    var_gene_names=rnames,
                    read_chunk_fn=read_chunk,
                    read_rows_fn=read_rows,
                    read_chunk_cols_fn=read_chunk_cols,
                    cell_barcodes=bcs_kept,
                )
            )
    return parts


def extract_pfizer_telohaec_ntc_bundle_r17(
    raw_path: Path,
    condition: str = "untreated",
    seed: int = DEFAULT_R17_SEED,
    resample_seed: int = DEFAULT_R17_RESAMPLE_SEED,
    target_n: int = 512,
    include_resample: bool = False,
    min_umi: float = 300.0,
    min_ngenes: int = 200,
    min_batch_nt: int = 15,
) -> Dict[str, Any]:
    """
    Bounded-RAM fast NTC + barcode extractor for Pfizer TeloHAEC (`imru_full.h5ad`)
    producing bit-identical QC-passing NTC set, batch labels, barcodes, and sampled raw counts
    as locked `pfizer_telohaec_reader_r17.py` without reading non-NTC rows from `/X`.
    """
    from dataset_readers_r12 import _resolve_h5ad_index_col  # noqa: E402
    from dataset_readers_r14 import (  # noqa: E402
        _get_h5ad_slot_node,
        _get_var_group_for_slot,
        _inspect_slot_shape_and_layout,
        _resolve_h5ad_integer_slot,
    )
    from prepare_r10_lines import build_hgnc_lookup  # noqa: E402

    cond_req = str(condition).strip().lower()
    if cond_req not in ("untreated", "treated"):
        raise ValueError(f"Unsupported Pfizer condition: {condition!r}")

    hgnc = build_hgnc_lookup()
    h5_path = _resolve_h5ad_path_r17(Path(raw_path))
    with h5py.File(h5_path, "r") as f:
        chosen_slot, _, _ = _resolve_h5ad_integer_slot(f, x_layer="auto", file_label=h5_path.name)
        obs_grp = f["obs"]
        obs_keys = sorted(k for k in obs_grp.keys() if k != "__categories")
        obs_idx_col = _resolve_h5ad_index_col(obs_grp, group_name="obs")
        obs_index = [x.strip() for x in read_h5ad_table_col_r17(obs_grp, obs_idx_col)]
        n_obs = len(obs_index)

        var_grp = _get_var_group_for_slot(f, chosen_slot)
        node = _get_h5ad_slot_node(f, chosen_slot)
        n_var = len(read_h5ad_table_col_r17(var_grp, _resolve_h5ad_index_col(var_grp, group_name="var")))
        (_, mat_c), layout = _inspect_slot_shape_and_layout(node, n_obs, n_var)
        var_meta = _extract_var_rna_features_r17(var_grp, n_var_expected=mat_c)

        c_cond = _pick_obs_col_r17(obs_keys, ["condition", "cell_treat", "treatment", "stim", "stimulation"])
        c_tgt = _pick_obs_col_r17(obs_keys, ["gene_target", "target_gene", "perturbation", "gene", "target"])
        c_guide = _pick_obs_col_r17(obs_keys, ["guide", "sgRNA", "sgrna", "guide_id", "sgrna_id"])
        c_cellid = _pick_obs_col_r17(obs_keys, ["cellID", "cell_id", "barcode"])
        c_flask = _pick_obs_col_r17(obs_keys, ["flask", "replicate", "batch", "gem_group", "channel", "lane"])
        c_well = _pick_obs_col_r17(obs_keys, ["welltag", "well", "plate_well"])
        c_nguide = _pick_obs_col_r17(obs_keys, ["num_features", "n_guides", "guide_count", "n_sgrna"])
        if c_cond is None or c_tgt is None:
            raise KeyError(f"Pfizer TeloHAEC .obs missing condition/gene_target: {obs_keys}")

        cond_vals = read_h5ad_table_col_r17(obs_grp, c_cond)
        tgt_vals = read_h5ad_table_col_r17(obs_grp, c_tgt)
        guide_vals = read_h5ad_table_col_r17(obs_grp, c_guide) if c_guide else [""] * n_obs
        cellid_vals = read_h5ad_table_col_r17(obs_grp, c_cellid) if c_cellid else obs_index
        flask_vals = read_h5ad_table_col_r17(obs_grp, c_flask) if c_flask else None
        well_vals = read_h5ad_table_col_r17(obs_grp, c_well) if c_well else None
        nguide_vals = read_h5ad_table_col_r17(obs_grp, c_nguide) if c_nguide else None

    # Filter strictly to candidate NTC cells (skipping 416k targeting cells for instant RNA QC!)
    ntc_cand_rows: List[int] = []
    ntc_cand_bcs: List[str] = []
    ntc_cand_flasks: List[str] = []
    ntc_cand_wells: List[str] = []
    for i in range(n_obs):
        if classify_pfizer_condition_r17(cond_vals[i]) != cond_req:
            continue
        if nguide_vals is not None:
            try:
                ng_i = float(nguide_vals[i])
                if np.isfinite(ng_i) and int(round(ng_i)) != 1:
                    continue
            except ValueError:
                pass
        cls, _, _ = classify_pfizer_guide_r17(tgt_vals[i], guide_vals[i], hgnc)
        if cls != "ntc":
            continue
        ntc_cand_rows.append(i)
        bc_base = str(cellid_vals[i]).strip() or obs_index[i]
        ntc_cand_bcs.append(f"{cond_req}:{bc_base}")
        ntc_cand_flasks.append(flask_vals[i] if flask_vals is not None else "")
        ntc_cand_wells.append(well_vals[i] if well_vals is not None else "")

    cand_arr = np.asarray(ntc_cand_rows, dtype=np.int64)
    umi_c, ng_c = _compute_rna_umi_and_ngenes_r17(h5_path, chosen_slot, cand_arr, var_meta["rna_cols"])
    ok = (umi_c >= float(min_umi)) & (ng_c >= int(min_ngenes))
    keep_local = np.flatnonzero(ok)
    keep_rows = cand_arr[keep_local]
    ntc_bcs = [ntc_cand_bcs[int(k)] for k in keep_local]
    flasks_k = [ntc_cand_flasks[int(k)] for k in keep_local] if flask_vals is not None else None
    wells_k = [ntc_cand_wells[int(k)] for k in keep_local] if well_vals is not None else None

    # Note: hardware split-half assignment depends solely on NTC counts per flask/welltag!
    batches_arr, half_flags, half_rule, _, _, _ = _assign_hardware_split_halves_r17(
        nt_mask=np.ones(len(keep_rows), dtype=bool),
        flasks=flasks_k,
        welltags=wells_k,
        seed=seed,
        min_batch_nt=min_batch_nt,
    )
    read_chunk_fn, read_rows_fn, read_chunk_cols_fn = _make_h5ad_subset_readers_r17(
        h5_path=h5_path,
        slot=chosen_slot,
        keep_rows=keep_rows,
        rna_cols=var_meta["rna_cols"],
        n_obs_total=n_obs,
        n_var_total=mat_c,
        layout=layout,
    )
    part = LightweightPartSpec(
        tag=f"Pfizer_TeloHAEC_{cond_req}",
        path=h5_path,
        n_cells=len(keep_rows),
        n_genes_raw=len(var_meta["rna_syms"]),
        nt_mask=np.ones(len(keep_rows), dtype=bool),
        batches=batches_arr,
        half_flags=half_flags,
        umi=umi_c[keep_local],
        pert_keys=[""] * len(keep_rows),
        pert_meta={},
        var_gene_ids=var_meta["rna_ids"],
        var_gene_names=var_meta["rna_syms"],
        read_chunk_fn=read_chunk_fn,
        read_rows_fn=read_rows_fn,
        read_chunk_cols_fn=read_chunk_cols_fn,
        cell_barcodes=ntc_bcs,
    )
    bundle = extract_ntc_bundle_from_parts_r16(
        [part],
        target_n=target_n,
        cell_set_len=32,
        min_ctrl=128,
        seed=seed,
        resample_seed=resample_seed,
        include_resample=include_resample,
    )
    bundle["ctrl_audit"]["half_rule_used"] = half_rule
    return bundle


def _extract_a549_cas9pos_ntc_bundle_r17(
    raw_dir: Path,
    seed: int = DEFAULT_R17_SEED,
    resample_seed: int = DEFAULT_R17_RESAMPLE_SEED,
    target_n: int = 512,
    include_resample: bool = False,
    bucket_size: int = 40000,
    chunk_lines: int = 1_500_000,
    min_batch_nt: int = 15,
    max_pct_mt: float = 20.0,
    min_ngenes: int = 500,
) -> Dict[str, Any]:
    """
    Read-only exact mirror of `dataset_readers_r12.build_a549_parts(..., exclude_zero_cas9=True)`
    strictly on NTC cells (`total_counts_CAS9 > 0` or MTX CAS9 UMI > 0; ZERO writes to `round12/`).
    """
    base = Path(raw_dir)
    adir = base / "a549_ko" if (base / "a549_ko").is_dir() else (base / "a549" if (base / "a549").is_dir() else base)
    cell_path = _find_file_r16(
        adir, ["GSE345058_cell_annotations.csv.gz", "GSE345058_cell_annotations.csv"], "*cell_annotations*.csv*"
    )
    gene_path = _find_file_r16(
        adir, ["GSE345058_gene_annotations.csv.gz", "GSE345058_gene_annotations.csv"], "*gene_annotations*.csv*"
    )
    mtx_path = _find_file_r16(
        adir, ["GSE345058_raw_counts.mtx.gz", "GSE345058_raw_counts.mtx"], "*raw_counts*.mtx*"
    )

    cell_df = pd.read_csv(cell_path, dtype=str, keep_default_na=False)
    gene_df = pd.read_csv(gene_path, dtype=str, keep_default_na=False)
    n_annot_cells = len(cell_df)
    n_annot_genes = len(gene_df)

    c_cellid_col = _pick_col_r16(cell_df.columns, ["cell_id", "Cell_ID", "barcode", "Barcode", "cell"])
    c_batch_col = _pick_col_r16(cell_df.columns, ["batch", "Batch", "batch_id"])
    c_unit_col = _pick_col_r16(cell_df.columns, ["Lane", "lane", "channel", "sample", "batch", "Batch"])
    c_sgrna_col = _pick_col_r16(cell_df.columns, ["sgRNA", "sgrna", "guide_id", "guide"])
    c_pert_col = _pick_col_r16(cell_df.columns, ["Perturbation", "perturbation", "target_gene", "gene_symbol", "gene"])
    c_cond_col = _pick_col_r16(cell_df.columns, ["Condition", "condition"])
    c_cas9_col = _pick_col_r16(cell_df.columns, ["total_counts_CAS9", "cas9_counts", "CAS9_counts", "CAS9"])
    if c_batch_col is None or c_pert_col is None or c_sgrna_col is None:
        raise KeyError(f"A549 cell annotations missing required columns in {cell_path}")

    g_sym_col = _pick_col_r16(gene_df.columns, ["gene_symbol", "gene_name", "symbol", "Symbol", "name", "gene"]) or gene_df.columns[0]
    g_ens_col = _pick_col_r16(gene_df.columns, ["ensembl_id", "gene_id", "gene_ids", "Ensembl", "ensembl"])
    g_type_col = _pick_col_r16(gene_df.columns, ["feature_types", "feature_type", "type"])
    g_gen_col = _pick_col_r16(gene_df.columns, ["genome", "species", "reference"])
    raw_g_syms = [str(x).strip() for x in gene_df[g_sym_col].tolist()]
    raw_g_ens = [str(x).strip().split(".")[0] for x in gene_df[g_ens_col].tolist()] if g_ens_col else [""] * n_annot_genes
    raw_g_gen = [str(x).strip().upper() for x in gene_df[g_gen_col].tolist()] if g_gen_col else [""] * n_annot_genes

    def _is_cas9(i_g: int) -> bool:
        return raw_g_syms[i_g].upper() == "CAS9" or raw_g_ens[i_g].upper() == "CAS9" or raw_g_gen[i_g] == "CAS9"

    if g_type_col is not None:
        g_types = [str(x).strip() for x in gene_df[g_type_col].tolist()]
        rna_gene_raw_idx = [i for i, t in enumerate(g_types) if t in ("Gene Expression", "RNA", "") and not _is_cas9(i)]
        if not rna_gene_raw_idx:
            rna_gene_raw_idx = [i for i in range(n_annot_genes) if not _is_cas9(i)]
    else:
        rna_gene_raw_idx = [i for i in range(n_annot_genes) if not _is_cas9(i)]

    raw_to_rna = np.full(n_annot_genes, -1, dtype=np.int32)
    rna_names: List[str] = []
    for k, ri in enumerate(rna_gene_raw_idx):
        raw_to_rna[ri] = k
        rna_names.append(raw_g_syms[ri])
    n_genes_rna = len(rna_names)
    mt_mask_rna = np.array([gn.startswith("MT-") for gn in rna_names], dtype=bool)
    cas9_raw_idx = next((i for i in range(n_annot_genes) if _is_cas9(i)), -1)
    cas9_vals = (
        pd.to_numeric(cell_df[c_cas9_col], errors="coerce").to_numpy(dtype=np.float64)
        if c_cas9_col is not None
        else None
    )

    cids_raw = (
        cell_df[c_cellid_col].astype(str).str.strip().tolist()
        if c_cellid_col is not None
        else [f"a549_row_{i}" for i in range(n_annot_cells)]
    )
    batches_raw = cell_df[c_batch_col].astype(str).str.strip().tolist()
    units_raw = cell_df[c_unit_col].astype(str).str.strip().tolist() if c_unit_col else batches_raw
    sgrnas_raw = cell_df[c_sgrna_col].astype(str).str.strip().tolist()
    perts_raw = cell_df[c_pert_col].astype(str).str.strip().tolist()
    conds_raw = cell_df[c_cond_col].astype(str).str.strip().tolist() if c_cond_col is not None else None

    raw_to_ntc = np.full(n_annot_cells, -1, dtype=np.int32)
    ntc_pre_cids: List[str] = []
    ntc_pre_units: List[str] = []
    ntc_pre_halves: List[int] = []
    for ci in range(n_annot_cells):
        sg = sgrnas_raw[ci]
        pr = perts_raw[ci]
        if not sg or sg.lower() in _A549_UNASSIGNED_TOKENS or _MULTI_GUIDE_DELIM_RE.search(sg) or " " in sg:
            continue
        if _MULTI_GUIDE_DELIM_RE.search(pr):
            continue
        m_b = re.search(r"(\d+)$", batches_raw[ci])
        b_num = int(m_b.group(1)) if m_b else None
        if b_num not in (1, 2, 3, 4):
            continue
        if cas9_vals is not None and not (np.isfinite(cas9_vals[ci]) and cas9_vals[ci] > 0.0):
            continue
        half_flag = 1 if b_num in (1, 3) else 2
        pert_is_ntc = pr.lower() in _A549_NTC_TOKENS
        if conds_raw is not None:
            cond_low = conds_raw[ci].lower()
            if cond_low in _A549_UNASSIGNED_TOKENS:
                continue
            cond_is_ntc = cond_low in _A549_NTC_TOKENS
            if pert_is_ntc != cond_is_ntc:
                continue
            is_ntc = pert_is_ntc and cond_is_ntc
        else:
            is_ntc = pert_is_ntc
        if is_ntc:
            raw_to_ntc[ci] = len(ntc_pre_cids)
            ntc_pre_cids.append(cids_raw[ci])
            ntc_pre_units.append(units_raw[ci])
            ntc_pre_halves.append(half_flag)

    n_pre_ntc = len(ntc_pre_cids)
    if n_pre_ntc == 0:
        raise RuntimeError(f"A549 CAS9>0 metadata yielded zero NTC cells in {adir}")

    n_mtx_rows, n_mtx_cols, _, n_skip_lines, sep = _read_mtx_header_r16(mtx_path)
    cell_first = (n_mtx_rows == n_annot_cells and n_mtx_cols == n_annot_genes)
    if not cell_first and not (n_mtx_rows == n_annot_genes and n_mtx_cols == n_annot_cells):
        raise RuntimeError(f"A549 MTX shape ({n_mtx_rows}, {n_mtx_cols}) incompatible with annotations.")

    n_buckets = max(1, (n_pre_ntc + bucket_size - 1) // bucket_size)
    cas9_from_mtx = np.zeros(n_pre_ntc, dtype=np.float32) if cas9_vals is None else None
    with tempfile.TemporaryDirectory(prefix="r17_a549_c9_ntc_") as tmp_dir:
        tmp_path = Path(tmp_dir)
        tmp_h5 = tmp_path / "a549_c9_ntc_csc.h5"
        f_g = [open(tmp_path / f"b_{bi}.g.bin", "wb") for bi in range(n_buckets)]
        f_c = [open(tmp_path / f"b_{bi}.c.bin", "wb") for bi in range(n_buckets)]
        f_v = [open(tmp_path / f"b_{bi}.v.bin", "wb") for bi in range(n_buckets)]
        try:
            reader = pd.read_csv(
                mtx_path, sep=sep, header=None, skiprows=n_skip_lines,
                names=["r", "c", "v"], dtype={"r": np.int32, "c": np.int32, "v": np.float32},
                chunksize=chunk_lines, engine="c",
            )
            for chunk in reader:
                rr = chunk["r"].to_numpy(dtype=np.int64) - 1
                cc = chunk["c"].to_numpy(dtype=np.int64) - 1
                vv = chunk["v"].to_numpy(dtype=np.float32)
                c_raw = rr if cell_first else cc
                g_raw = cc if cell_first else rr
                ntc_ci = raw_to_ntc[c_raw]
                if cas9_from_mtx is not None and cas9_raw_idx >= 0:
                    m_c9 = (ntc_ci >= 0) & (g_raw == cas9_raw_idx) & (vv > 0)
                    if np.any(m_c9):
                        np.add.at(cas9_from_mtx, ntc_ci[m_c9], vv[m_c9])
                g_rna = raw_to_rna[g_raw]
                m_ok = (ntc_ci >= 0) & (g_rna >= 0) & (vv != 0.0)
                if not np.any(m_ok):
                    continue
                k_arr = ntc_ci[m_ok].astype(np.int32)
                gr_arr = g_rna[m_ok].astype(np.int32)
                vr_arr = vv[m_ok].astype(np.float32)
                b_ids = k_arr // bucket_size
                for bi in np.unique(b_ids):
                    mb = b_ids == bi
                    gr_arr[mb].tofile(f_g[int(bi)])
                    (k_arr[mb] - int(bi) * bucket_size).astype(np.int32).tofile(f_c[int(bi)])
                    vr_arr[mb].tofile(f_v[int(bi)])
        finally:
            for fh in f_g + f_c + f_v:
                fh.close()

        umi_pre = np.zeros(n_pre_ntc, dtype=np.float64)
        ng_pre = np.zeros(n_pre_ntc, dtype=np.int32)
        mt_umi_pre = np.zeros(n_pre_ntc, dtype=np.float64)
        with h5py.File(tmp_h5, "w") as hf:
            mg = hf.create_group("matrix")
            mg.create_dataset("shape", data=np.array([n_genes_rna, n_pre_ntc], dtype=np.int64))
            data_ds = mg.create_dataset("data", shape=(0,), maxshape=(None,), dtype=np.float32, chunks=(65536,))
            idx_ds = mg.create_dataset("indices", shape=(0,), maxshape=(None,), dtype=np.int32, chunks=(65536,))
            indptr_ds = mg.create_dataset("indptr", shape=(n_pre_ntc + 1,), dtype=np.int64)
            indptr_ds[0] = 0
            cur_nnz = 0
            for bi in range(n_buckets):
                cs = bi * bucket_size
                ce = min(cs + bucket_size, n_pre_ntc)
                nc_b = ce - cs
                ga = np.fromfile(tmp_path / f"b_{bi}.g.bin", dtype=np.int32)
                ca = np.fromfile(tmp_path / f"b_{bi}.c.bin", dtype=np.int32)
                va = np.fromfile(tmp_path / f"b_{bi}.v.bin", dtype=np.float32)
                if ga.size > 0:
                    csc_b = sp.coo_matrix((va, (ga, ca)), shape=(n_genes_rna, nc_b)).tocsc()
                    umi_pre[cs:ce] = np.asarray(csc_b.sum(axis=0)).ravel()
                    ng_pre[cs:ce] = np.diff(csc_b.indptr).astype(np.int32)
                    if np.any(mt_mask_rna):
                        mt_umi_pre[cs:ce] = np.asarray(csc_b[mt_mask_rna, :].sum(axis=0)).ravel()
                    add_nnz = int(csc_b.nnz)
                    if add_nnz > 0:
                        data_ds.resize((cur_nnz + add_nnz,))
                        idx_ds.resize((cur_nnz + add_nnz,))
                        data_ds[cur_nnz : cur_nnz + add_nnz] = csc_b.data
                        idx_ds[cur_nnz : cur_nnz + add_nnz] = csc_b.indices
                    indptr_ds[cs + 1 : ce + 1] = cur_nnz + csc_b.indptr[1:].astype(np.int64)
                    cur_nnz += add_nnz
                else:
                    indptr_ds[cs + 1 : ce + 1] = cur_nnz

        pct_mt_pre = 100.0 * mt_umi_pre / np.maximum(umi_pre, 1.0)
        qc_ok = (umi_pre > 0.0) & (ng_pre >= int(min_ngenes)) & (pct_mt_pre <= float(max_pct_mt))
        if cas9_from_mtx is not None:
            qc_ok &= cas9_from_mtx > 0.0
        keep_idx = np.where(qc_ok)[0].astype(np.int64)
        ntc_cids = [ntc_pre_cids[int(i)] for i in keep_idx]
        ntc_units = np.asarray([ntc_pre_units[int(i)] for i in keep_idx], dtype=str)
        ntc_halves = np.asarray([ntc_pre_halves[int(i)] for i in keep_idx], dtype=np.int8)
        ntc_batches_arr, _ = _apply_min_ntc_batch_fallback_r16(
            ntc_units, np.ones(len(keep_idx), dtype=bool), ntc_halves, min_batch_nt=min_batch_nt
        )
        ntc_batches = [str(b) for b in ntc_batches_arr]

        sets = sample_stratified_control_sets_r16(
            n_total=len(keep_idx), batch_labels=ntc_batches, target_n=target_n, cell_set_len=32, min_ctrl=128, seed=seed
        )
        flat_loc = sets.ravel()
        _, read_rows, _ = _make_10x_csc_subset_readers_r16(tmp_h5, keep_idx, n_genes_rna, None)
        raw_mat = np.asarray(read_rows(flat_loc), dtype=np.float32)
        sel_batches = [ntc_batches[int(i)] for i in flat_loc]
        sel_barcodes = [ntc_cids[int(i)] for i in flat_loc]

        resamp_raw = resamp_batches = resamp_barcodes = None
        resamp_audit: Dict[str, Any] = {"eligible": False, "reason": "not_requested"}
        if include_resample:
            resamp_sets, resamp_audit = sample_disjoint_resample_control_sets_r16(
                n_total=len(keep_idx), primary_sets=sets, batch_labels=ntc_batches,
                target_n=target_n, cell_set_len=32, min_ctrl=128, seed=resample_seed,
            )
            if resamp_sets is not None:
                r_loc = resamp_sets.ravel()
                resamp_raw = np.asarray(read_rows(r_loc), dtype=np.float32)
                resamp_batches = [ntc_batches[int(i)] for i in r_loc]
                resamp_barcodes = [ntc_cids[int(i)] for i in r_loc]
                resamp_audit["sampled_barcodes"] = resamp_barcodes
                resamp_audit["sampled_batch_labels"] = resamp_batches

    return {
        "raw_counts": raw_mat,
        "gene_names": rna_names,
        "batch_labels": sel_batches,
        "barcodes": sel_barcodes,
        "ctrl_audit": {
            "source_file": str(mtx_path),
            "exclude_zero_cas9": True,
            "n_ctrl_total": int(len(keep_idx)),
            "n_ctrl_sampled": int(len(sel_barcodes)),
            "n_sets": int(sets.shape[0]),
            "cell_set_len": int(sets.shape[1]),
            "ctrl_sampling_fraction": float(len(sel_barcodes) / max(len(keep_idx), 1)),
            "n_raw_genes": int(n_genes_rna),
            "seed": int(seed),
            "sampled_barcodes": sel_barcodes,
            "sampled_batch_labels": sel_batches,
            "within_line_resample": resamp_audit,
        },
        "resample_raw_counts": resamp_raw,
        "resample_batch_labels": resamp_batches,
        "resample_barcodes": resamp_barcodes,
        "resample_audit": resamp_audit,
    }


class _CommonAxisPartProxyR17:
    """Read-only proxy restricting a PartSpec's gene axis onto an exact common-key column map (DEVIATION_R17_D1)."""

    def __init__(
        self,
        part: Any,
        col_idx: np.ndarray,
        common_gene_names: List[str],
        common_gene_ids: Optional[List[str]],
    ) -> None:
        self._part = part
        self._col_idx = np.asarray(col_idx, dtype=np.int64)
        self.var_gene_names = list(common_gene_names)
        if common_gene_ids is not None:
            self.var_gene_ids = list(common_gene_ids)
        self.n_genes_raw = int(len(common_gene_names))

        orig_rows_fn = getattr(part, "read_rows_fn", None)
        if orig_rows_fn is not None:
            def _read_rows(
                idx_arr: np.ndarray,
                _fn: Any = orig_rows_fn,
                _c: np.ndarray = self._col_idx,
            ) -> np.ndarray:
                return np.asarray(_fn(idx_arr), dtype=np.float32)[:, _c]

            self.read_rows_fn = _read_rows
        else:
            self.read_rows_fn = None

        orig_chunk_fn = getattr(part, "read_chunk_fn", None)
        if orig_chunk_fn is not None:
            def _read_chunk(
                s: int,
                e: int,
                _fn: Any = orig_chunk_fn,
                _c: np.ndarray = self._col_idx,
            ) -> np.ndarray:
                return np.asarray(_fn(s, e), dtype=np.float32)[:, _c]

            self.read_chunk_fn = _read_chunk
        else:
            self.read_chunk_fn = None

    def __getattr__(self, name: str) -> Any:
        return getattr(self._part, name)


def _restrict_parts_to_common_axis_r17(parts: Sequence[Any]) -> Tuple[List[Any], Dict[str, Any]]:
    """
    `DEVIATION_R17_D1`: Restrict multi-part lines whose files/channels have unequal gene sets
    (specifically `a549_myllia`: `ch_1` 18,711 vs `ch_2` 18,808 -> 18,480 common keys) onto the
    exact Round 10 `process_r10_line` common-key intersection (`key = stripped gene name or gene id`,
    in part-0 order) before calling frozen `extract_ntc_bundle_from_parts_r16`.
    Parts already matching `parts[0]` or `DEVIATION_R16_D1` are returned untouched.
    """
    parts_list = list(parts)
    if len(parts_list) <= 1:
        return parts_list, {"applied": False, "n_genes_multifile_dropped": 0}

    ref_genes = list(parts_list[0].var_gene_names)
    ref_set = set(ref_genes)
    all_compatible_r16_d1 = True
    for part in parts_list[1:]:
        pg = list(part.var_gene_names)
        if pg == ref_genes:
            continue
        if len(set(pg)) == len(pg) == len(ref_genes) == len(ref_set) and set(pg) == ref_set:
            continue
        all_compatible_r16_d1 = False
        break

    if all_compatible_r16_d1:
        return parts_list, {"applied": False, "n_genes_multifile_dropped": 0}

    sym_maps: List[Dict[str, int]] = []
    for part in parts_list:
        raw_names = list(part.var_gene_names)
        raw_ids = list(getattr(part, "var_gene_ids", None) or ([""] * len(raw_names)))
        smap: Dict[str, int] = {}
        for j, (gid, gnm) in enumerate(zip(raw_ids, raw_names)):
            key = str(gnm).strip() or str(gid).strip()
            if key and key not in smap:
                smap[key] = j
        sym_maps.append(smap)

    common_keys = [k for k in sym_maps[0] if all(k in sm for sm in sym_maps[1:])]
    if not common_keys:
        raise ValueError("Multi-part common-key gene intersection is empty (DEVIATION_R17_D1).")
    col_maps = [np.asarray([sm[k] for k in common_keys], dtype=np.int64) for sm in sym_maps]
    p0_ids = getattr(parts_list[0], "var_gene_ids", None)
    common_ids = [p0_ids[int(j)] for j in col_maps[0]] if p0_ids is not None else None

    wrapped = [
        _CommonAxisPartProxyR17(part, cmap, common_keys, common_ids)
        for part, cmap in zip(parts_list, col_maps)
    ]
    raw_total_p0 = len(set(str(x).strip() for x in parts_list[0].var_gene_names))
    return wrapped, {
        "applied": True,
        "rule": "DEVIATION_R17_D1_R10_common_key_intersection_part0_order",
        "n_genes_part0": int(len(parts_list[0].var_gene_names)),
        "n_genes_common": int(len(common_keys)),
        "n_genes_multifile_dropped": int(raw_total_p0 - len(common_keys)),
    }


def _extract_bundle_from_frozen_reader_r17(
    line_name: str,
    seed: int = DEFAULT_R17_SEED,
    resample_seed: int = DEFAULT_R17_RESAMPLE_SEED,
    target_n: int = 512,
    include_resample: bool = False,
) -> Dict[str, Any]:
    """
    Read-only NTC + barcode extractor for smaller fresh local lines (`a549_myllia`, `hek293_xucao`,
    `telohaec`, `calu3`, `mac_ctrl`, `mac_lps`), calling frozen Round 10–14 `build_*_parts` in
    read-only mode (ZERO files written to `round10/`–`round14/`).
    """
    from prepare_r10_lines import build_hgnc_lookup  # noqa: E402

    hgnc_l = build_hgnc_lookup()
    if line_name == "a549_myllia":
        from dataset_readers_r14 import build_gse337804_parts  # noqa: E402

        parts, audit = build_gse337804_parts(REPO_ROOT / "round14" / "shared" / "raw", hgnc_l, line="a549_myllia")
    elif line_name == "hek293_xucao":
        from dataset_readers_r12 import build_hek293_parts  # noqa: E402

        parts, audit = build_hek293_parts(REPO_ROOT / "round12" / "shared" / "raw", hgnc_l)
    elif line_name == "telohaec":
        # DEVIATION_R17_D2: use the exact Round 11 observed-data builder (`prepare_r11_lines.py`).
        from prepare_r11_lines import build_telohaec_parts  # noqa: E402

        parts, audit = build_telohaec_parts(REPO_ROOT / "round11" / "shared" / "raw", hgnc_l)
    elif line_name == "calu3":
        from prepare_r10_lines import build_calu3_parts  # noqa: E402

        parts, audit = build_calu3_parts(REPO_ROOT / "round10" / "shared" / "raw", hgnc_l)
    elif line_name in ("mac_ctrl", "mac_lps"):
        from dataset_readers_r13 import build_gse327644_parts  # noqa: E402

        parts, audit = build_gse327644_parts(REPO_ROOT / "round13" / "shared" / "raw", hgnc_l, line=line_name)
    else:
        raise ValueError(f"Unsupported small fresh line: {line_name!r}")

    parts, d1_axis_audit = _restrict_parts_to_common_axis_r17(parts)
    bundle = extract_ntc_bundle_from_parts_r16(
        parts=parts,
        target_n=target_n,
        cell_set_len=32,
        min_ctrl=128,
        seed=seed,
        resample_seed=resample_seed,
        include_resample=include_resample,
    )
    bundle["ctrl_audit"]["reader_audit"] = audit
    if d1_axis_audit.get("applied"):
        bundle["ctrl_audit"]["deviation_r17_d1_axis_audit"] = d1_axis_audit
    return bundle


def extract_control_ntc_bundle_r17(
    pool_name: str,
    status_dir: Path = DEFAULT_R17_STATUS_DIR,
    prereg_lock_path: Optional[Path] = DEFAULT_R17_PREREG_LOCK_PATH,
    seed: int = DEFAULT_R17_SEED,
    resample_seed: int = DEFAULT_R17_RESAMPLE_SEED,
    include_resample: bool = False,
    gse295214_raw_dir: Optional[Path] = None,
    pfizer_raw_path: Optional[Path] = None,
    raw_nadig_dir: Optional[Path] = None,
) -> Dict[str, Any]:
    """
    Strictly gated, read-only raw-count NTC + barcode extractor covering ALL 11 evaluated Round 17
    target lines AND all 6 Round 16 mismatched control pools.
    """
    assert_all_execution_gates_r17(status_dir=status_dir, prereg_lock_path=prereg_lock_path)

    # 1. Round 16 mismatched control pools (reuse frozen R16 extractors directly, zero writes)
    if pool_name in ("nadig_hepg2", "nadig_jurkat"):
        p = _resolve_nadig_h5ad_path_r16(pool_name, raw_nadig_dir=raw_nadig_dir)
        return extract_nadig_h5ad_ntc_bundle_r16(
            p, seed=seed, resample_seed=resample_seed, include_resample=include_resample
        )
    if pool_name == "a549_ko":
        return _extract_a549_ntc_bundle_r16(
            REPO_ROOT / "round12" / "shared" / "raw",
            seed=seed, resample_seed=resample_seed, include_resample=include_resample,
        )
    if pool_name == "dld1":
        parts = _build_dld1_ntc_lightweight_parts_r16(REPO_ROOT / "round10" / "shared" / "raw")
        return extract_ntc_bundle_from_parts_r16(
            parts, seed=seed, resample_seed=resample_seed, include_resample=include_resample
        )
    if pool_name == "cd4t_rest":
        parts = _build_cd4t_rest_ntc_lightweight_parts_r16(REPO_ROOT / "round10" / "shared" / "raw")
        return extract_ntc_bundle_from_parts_r16(
            parts, seed=seed, resample_seed=resample_seed, include_resample=include_resample
        )
    if pool_name == "feng_ipsc":
        return _extract_feng_ipsc_ntc_bundle_r16(
            REPO_ROOT / "round11" / "shared" / "raw",
            seed=seed, resample_seed=resample_seed, include_resample=include_resample,
        )

    # 2. Round 17 fresh & confirmatory lines
    if pool_name == "cd4t_stim48":
        parts = _build_cd4t_stim48_ntc_lightweight_parts_r17(REPO_ROOT / "round11" / "shared" / "raw")
        return extract_ntc_bundle_from_parts_r16(
            parts, seed=seed, resample_seed=resample_seed, include_resample=include_resample
        )
    if pool_name == "gse295214_h1":
        raw_d = gse295214_raw_dir or DEFAULT_GSE295214_RAW_DIR
        bridge = HGNCBridge.from_tsv(DEFAULT_HGNC_FILE)
        parts, reader_audit = build_gse295214_parts(raw_dir=raw_d, hgnc=bridge, use_r10_partspec=False)
        b = extract_ntc_bundle_from_parts_r16(
            parts, seed=seed, resample_seed=resample_seed, include_resample=include_resample
        )
        b["ctrl_audit"]["reader_audit"] = reader_audit
        return b
    if pool_name == "pfizer_telohaec_basal":
        p_raw = pfizer_raw_path or DEFAULT_PFIZER_RAW_PATH
        return extract_pfizer_telohaec_ntc_bundle_r17(
            p_raw, condition="untreated", seed=seed, resample_seed=resample_seed, include_resample=include_resample
        )
    if pool_name == "pfizer_telohaec_cytokine":
        p_raw = pfizer_raw_path or DEFAULT_PFIZER_RAW_PATH
        return extract_pfizer_telohaec_ntc_bundle_r17(
            p_raw, condition="treated", seed=seed, resample_seed=resample_seed, include_resample=include_resample
        )
    if pool_name == "a549_ko_cas9pos":
        return _extract_a549_cas9pos_ntc_bundle_r17(
            REPO_ROOT / "round12" / "shared" / "raw",
            seed=seed, resample_seed=resample_seed, include_resample=include_resample,
        )
    if pool_name in ("a549_myllia", "hek293_xucao", "telohaec", "calu3", "mac_ctrl", "mac_lps"):
        return _extract_bundle_from_frozen_reader_r17(
            pool_name, seed=seed, resample_seed=resample_seed, include_resample=include_resample
        )
    raise ValueError(f"Unsupported control pool: {pool_name!r}")


# Keep unused-import-free exports explicit
__all__ = [
    "ALL_EVALUATED_LINES_R17",
    "CONFIRMATORY_LINES_R17",
    "DEFAULT_GSE295214_RAW_DIR",
    "DEFAULT_MODELS_DIR",
    "DEFAULT_PFIZER_RAW_PATH",
    "DEFAULT_R17_PREREG_LOCK_PATH",
    "DEFAULT_R17_RESULTS_DIR",
    "DEFAULT_R17_SHARED_DIR",
    "DEFAULT_R17_STATUS_DIR",
    "LINE_SPECS_R17",
    "MICRO_N_LINES_R17",
    "OMITTED_LINES_R17",
    "R16_DISCOVERY_QUARANTINED_LINES",
    "SECONDARY_FRESH_LINES_R17",
    "SUBSET_SENSITIVITY_LINES_R17",
    "WITHIN_FAMILY_LINES_R17",
    "BlindingGateError",
    "assert_all_execution_gates_r17",
    "build_mismatched_strata_r17",
    "canonicalize_gene_r17",
    "extract_control_ntc_bundle_r17",
    "extract_pfizer_telohaec_ntc_bundle_r17",
    "validate_target_lines_r17",
]
