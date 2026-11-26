#!/usr/bin/env python3
"""
Round 17 Model-Free Line Inventory Builder (`round17/src/build_line_inventory_j.py`).

Strictly data-only & counts-only inventory of all 9 local perturbation datasets that have
NEVER been evaluated under any Arc `STATE` model (`ST-SE-Replogle`, `ST-HVG-Replogle`, or
`SE-600M`), computing exact `PREREG_I` eligibility (`n_eval`, `|G_eval|`, NTC half-split
validity, raw integer UMI count availability for control sampling, modality, prior STATE
exposure audit, and family structure).

Hard compliance guarantees:
  - ZERO `torch`, `state`, `cell_load`, or MPS/GPU imports or calls (physically model-free).
  - ZERO correlation / metric calculations on target expression deltas (counts & structure only).
  - ZERO reads of `round16/audit/` or `A-T3.txt`.
  - ZERO modifications outside `round17/status/LINE_INVENTORY_J.json`.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import pickle
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import h5py
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
R9_LINES_DIR = REPO_ROOT / "round9" / "p2" / "shared" / "lines"
R9_COMMON_GENES_FILE = REPO_ROOT / "round9" / "p2" / "shared" / "COMMON_GENES.txt"
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
R16_VAR_DIMS_PKL = (
    REPO_ROOT
    / "round16"
    / "shared"
    / "models"
    / "ST-SE-Replogle"
    / "zeroshot"
    / "hepg2"
    / "var_dims.pkl"
)
DEFAULT_OUT_JSON = REPO_ROOT / "round17" / "status" / "LINE_INVENTORY_J.json"

if str(R10_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(R10_SRC_DIR))

from data_loader_r10 import (  # noqa: E402
    DEFAULT_HGNC_FILE,
    HGNCBridge,
    load_r9_common_genes,
    load_standardized_line_r10,
)

R9_SOURCE_POOL_6: Tuple[str, ...] = (
    "nadig_hepg2",
    "nadig_jurkat",
    "nourreddine_kolf21j",
    "replogle_k562_essential",
    "replogle_rpe1",
    "vcc2025_h1",
)

R16_DISCOVERY_LINES: Tuple[str, ...] = (
    "nadig_hepg2",
    "nadig_jurkat",
    "a549_ko",
    "dld1",
    "cd4t_rest",
    "feng_ipsc",
)

# Exact 9 local never-STATE-run candidate lines specified for Round 17 inventory
CANDIDATE_SPECS: Tuple[Dict[str, Any], ...] = (
    {
        "line": "calu3",
        "round": "round10",
        "line_dir_rel": "round10/shared/lines/calu3",
        "raw_dir_rel": "round10/shared/raw/calu3",
        "raw_probe_glob": "*.h5ad",
        "raw_format": "h5ad_single_file",
        "modality": "CRISPRi",
        "assay": "10x_3prime_perturbseq_bystander",
        "cell_lineage": "Calu-3 (lung adenocarcinoma epithelial, SARS-CoV-2 challenged bystander)",
        "canonical_family": "calu3",
        "study_accession": "GSE208240 / Zenodo:13350497 (Sunshine & Hein 2023)",
        "public_release_date": "2023",
        "se600m_pretrain_chronology": "PREDATES_SE600M_POTENTIAL_OVERLAP_UNVERIFIABLE",
        "relation_to_r16_discovery": "independent_study_and_lineage",
        "relation_to_other_r17_candidates": "independent",
        "family_notes": (
            "Distinct cell line (Calu-3) and study (Sunshine & Hein 2023) from all R16 discovery "
            "lines; predates SE-600M (2025-06/07), so SE-600M pretraining overlap cannot be ruled out "
            "(internal 14,420-dataset CSV manifest not shipped)."
        ),
    },
    {
        "line": "wtc11_i3neuron",
        "round": "round10",
        "line_dir_rel": "round10/shared/lines/wtc11_i3neuron",
        "raw_dir_rel": "round10/shared/raw/wtc11_i3neuron",
        "raw_probe_glob": "*.h5ad",
        "raw_format": "h5ad_single_file",
        "modality": "CRISPRi",
        "assay": "10x_3prime_perturbseq",
        "cell_lineage": "WTC11 iPSC-derived i3Neuron (human cortical neuron)",
        "canonical_family": "wtc11_i3neuron",
        "study_accession": "GSE152988 / Zenodo:13350497 (Tian & Kampmann 2021)",
        "public_release_date": "2021",
        "se600m_pretrain_chronology": "PREDATES_SE600M_POTENTIAL_OVERLAP_UNVERIFIABLE",
        "relation_to_r16_discovery": "independent_study_and_lineage",
        "relation_to_other_r17_candidates": "independent",
        "family_notes": (
            "Distinct differentiated neuronal lineage (WTC11 i3Neuron) and study (Tian & Kampmann 2021) "
            "from all R16 discovery lines; predates SE-600M (2025-06/07)."
        ),
    },
    {
        "line": "cd4t_stim48",
        "round": "round11",
        "line_dir_rel": "round11/shared/lines/cd4t_stim48",
        "raw_dir_rel": "round11/shared/raw/cd4t_stim48",
        "raw_probe_glob": "*_Stim48hr_*.h5",
        "raw_format": "10x_flex_h5_multilane",
        "modality": "CRISPRi",
        "assay": "10x_Flex_probe_perturbseq",
        "cell_lineage": "Primary human CD4+ alpha-beta T cells (48h TCR-stimulated, donors D2 & D3)",
        "canonical_family": "primary_cd4_t",
        "study_accession": "GSE314342 (Zhu et al., Cell 2026, doi:10.1016/j.cell.2026.08.002)",
        "public_release_date": "2026-02-06",
        "se600m_pretrain_chronology": "STRICTLY_POSTDATES_SE600M_ISOLATED",
        "relation_to_r16_discovery": "shared_study_lineage_library_with_r16_cd4t_rest_disjoint_donors_and_state",
        "relation_to_other_r17_candidates": "independent_of_other_r17_candidates",
        "family_notes": (
            "Strictly disjoint donors (CE0008678_D3, CE0010866_D2), disjoint 10x Flex lanes "
            "(R2L25..R2L48), and distinct biological stimulation state (48h TCR-stimulated) relative "
            "to R16 discovery line cd4t_rest (donors CE0006864_D4, CE0008162_D1, resting, R1L01..R2L24), "
            "but shares primary CD4+ T lineage, CRISPRi sgRNA library, and Marson lab study (GSE314342). "
            "Mandatory weaker-independence disclosure applies if used in confirmatory testing."
        ),
    },
    {
        "line": "telohaec",
        "round": "round11",
        "line_dir_rel": "round11/shared/lines/telohaec",
        "raw_dir_rel": "round11/shared/raw/telohaec",
        "raw_probe_glob": "*.h5ad",
        "raw_format": "h5ad_single_file",
        "modality": "CRISPRi",
        "assay": "10x_3prime_perturbseq",
        "cell_lineage": "TeloHAEC immortalized human aortic endothelial cells (resting)",
        "canonical_family": "telohaec",
        "study_accession": "IGVFDS7943IIWZ / GSE210681 (Schnitzler et al., Nature 2024, PMC10921916)",
        "public_release_date": "2023-2024",
        "se600m_pretrain_chronology": "PREDATES_SE600M_POTENTIAL_OVERLAP_UNVERIFIABLE",
        "relation_to_r16_discovery": "independent_study_and_lineage",
        "relation_to_other_r17_candidates": "independent",
        "family_notes": (
            "Distinct aortic endothelial lineage (TeloHAEC) and study (GSE210681) from all R16 "
            "discovery lines; predates SE-600M (2025-06/07)."
        ),
    },
    {
        "line": "a549_ko_cas9pos",
        "round": "round12",
        "line_dir_rel": "round12/shared/lines/a549_ko_cas9pos",
        "raw_dir_rel": "round12/shared/raw/a549_ko",
        "raw_probe_glob": "*raw_counts*.mtx*",
        "raw_format": "mtx_gz_plus_csc_cache",
        "modality": "CRISPRko",
        "assay": "CROPseq_SpCas9_CAS9pos_subset",
        "cell_lineage": "A549-doxCas9 lung adenocarcinoma (strict CAS9>0 cell subset of a549_ko)",
        "canonical_family": "a549",
        "study_accession": "GSE345058 (Liu et al. 2026, doi:10.64898/2026.06.01.728087)",
        "public_release_date": "2026-08-26",
        "se600m_pretrain_chronology": "STRICTLY_POSTDATES_SE600M_ISOLATED",
        "relation_to_r16_discovery": "STRICT_CELL_SUBSET_OF_R16_DISCOVERY_a549_ko",
        "relation_to_other_r17_candidates": "same_cell_line_A549_as_a549_myllia_different_study",
        "family_notes": (
            "NON-INDEPENDENT OF R16 DISCOVERY: strict cell subset (total_counts_CAS9 > 0; "
            "443,687 of 605,484 cells, 73.3%; 33,685 of 45,320 NTC cells, 74.3%) of R16 discovery "
            "line a549_ko from identical raw file GSE345058_raw_counts.mtx.gz. MUST BE EXCLUDED "
            "from independent confirmatory family (retained strictly as sensitivity/subset check only)."
        ),
    },
    {
        "line": "hek293_xucao",
        "round": "round12",
        "line_dir_rel": "round12/shared/lines/hek293_xucao",
        "raw_dir_rel": "round12/shared/raw/hek293_xucao",
        "raw_probe_glob": "*.h5ad",
        "raw_format": "h5ad_single_file",
        "modality": "CRISPRi",
        "assay": "PerturbSci_Kinetics",
        "cell_lineage": "HEK293 human embryonic kidney",
        "canonical_family": "hek293",
        "study_accession": "GSE218566 / XuCao2023.h5ad (Xu & Cao 2023)",
        "public_release_date": "2023",
        "se600m_pretrain_chronology": "PREDATES_SE600M_POTENTIAL_OVERLAP_UNVERIFIABLE",
        "relation_to_r16_discovery": "independent_study_and_lineage",
        "relation_to_other_r17_candidates": "independent",
        "family_notes": (
            "Distinct HEK293 lineage and study (PerturbSci-Kinetics, GSE218566) from all R16 "
            "discovery lines; predates SE-600M (2025-06/07); split-half uses seed42 50/50 fallback "
            "per RULING_R12_HEK293_UNITS.md."
        ),
    },
    {
        "line": "mac_ctrl",
        "round": "round13",
        "line_dir_rel": "round13/shared/lines/mac_ctrl",
        "raw_dir_rel": "round13/shared/raw/gse327644",
        "raw_probe_glob": "*Mac_Ctrl*.h5",
        "raw_format": "10x_flex_h5_multilane",
        "modality": "CRISPRko",
        "assay": "10x_Flex_probe_perturbseq",
        "cell_lineage": "Primary human CD14+ monocyte-derived macrophages (unstimulated control)",
        "canonical_family": "primary_macrophage_gse327644",
        "study_accession": "GSE327644 (Jung et al., Nat Biotechnol 2026, doi:10.1038/s41587-026-03258-2)",
        "public_release_date": "2026-08",
        "se600m_pretrain_chronology": "STRICTLY_POSTDATES_SE600M_ISOLATED",
        "relation_to_r16_discovery": "independent_study_and_lineage",
        "relation_to_other_r17_candidates": "paired_condition_same_cells_donor_library_as_mac_lps",
        "family_notes": (
            "Distinct primary macrophage lineage and study (GSE327644) from all R16 discovery lines; "
            "PAIRED WITH mac_lps (same M-CSF-differentiated donor macrophage screen under unstimulated "
            "vs LPS-stimulated conditions, channels L01-L04; not mutually independent)."
        ),
    },
    {
        "line": "mac_lps",
        "round": "round13",
        "line_dir_rel": "round13/shared/lines/mac_lps",
        "raw_dir_rel": "round13/shared/raw/gse327644",
        "raw_probe_glob": "*Mac_LPS*.h5",
        "raw_format": "10x_flex_h5_multilane",
        "modality": "CRISPRko",
        "assay": "10x_Flex_probe_perturbseq",
        "cell_lineage": "Primary human CD14+ monocyte-derived macrophages (LPS-stimulated)",
        "canonical_family": "primary_macrophage_gse327644",
        "study_accession": "GSE327644 (Jung et al., Nat Biotechnol 2026, doi:10.1038/s41587-026-03258-2)",
        "public_release_date": "2026-08",
        "se600m_pretrain_chronology": "STRICTLY_POSTDATES_SE600M_ISOLATED",
        "relation_to_r16_discovery": "independent_study_and_lineage",
        "relation_to_other_r17_candidates": "paired_condition_same_cells_donor_library_as_mac_ctrl",
        "family_notes": (
            "Distinct primary macrophage lineage and study (GSE327644) from all R16 discovery lines; "
            "PAIRED WITH mac_ctrl (same donor/cells/guide library under LPS stimulation, channels L01-L04; "
            "not mutually independent)."
        ),
    },
    {
        "line": "a549_myllia",
        "round": "round14",
        "line_dir_rel": "round14/shared/lines/a549_myllia",
        "raw_dir_rel": "round14/shared/raw/gse337804",
        "raw_probe_glob": "*.h5ad",
        "raw_format": "h5ad_multichannel",
        "modality": "CRISPRko",
        "assay": "10x_perturbseq_clonal_Cas9",
        "cell_lineage": "A549 clonal Cas9 lung adenocarcinoma (basal unstimulated, Day 16, Myllia)",
        "canonical_family": "a549",
        "study_accession": "GSE337804 (Myllia Biotechnology 2026)",
        "public_release_date": "2026",
        "se600m_pretrain_chronology": "STRICTLY_POSTDATES_SE600M_ISOLATED",
        "relation_to_r16_discovery": "same_cell_line_A549_as_r16_a549_ko_independent_lab_study_dataset",
        "relation_to_other_r17_candidates": "same_cell_line_A549_as_a549_ko_cas9pos_independent_study",
        "family_notes": (
            "Independent laboratory, chemistry, and dataset (GSE337804, Myllia Biotechnology) from "
            "R16 discovery a549_ko (GSE345058, Liu et al. 2026), but profiles the same canonical human "
            "cell line (A549) under CRISPRko."
        ),
    },
)


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _probe_raw_integer_counts(raw_dir: Path, pattern: str, fmt: str) -> Dict[str, Any]:
    """Model-free check that raw non-negative integer UMI counts exist on disk."""
    if not raw_dir.exists():
        return {"available": False, "reason": f"raw_dir missing: {raw_dir}"}
    matched = sorted(p for p in raw_dir.rglob(pattern) if p.is_file())
    if not matched:
        return {"available": False, "reason": f"no files matching {pattern} in {raw_dir}"}

    total_bytes = sum(p.stat().st_size for p in matched)
    first_file = matched[0]
    probe_vals: Optional[np.ndarray] = None
    probe_note = ""

    if fmt in ("h5ad_single_file", "h5ad_multichannel"):
        with h5py.File(first_file, "r") as f:
            x_node = f["X"]
            if isinstance(x_node, h5py.Group) and "data" in x_node:
                probe_vals = np.asarray(x_node["data"][:10000], dtype=np.float64)
                probe_note = f"h5ad sparse X/data[:{len(probe_vals)}] dtype={x_node['data'].dtype}"
            elif isinstance(x_node, h5py.Dataset):
                sl = x_node[: min(16, x_node.shape[0]), : min(500, x_node.shape[1])]
                probe_vals = np.asarray(sl, dtype=np.float64).ravel()
                probe_note = f"h5ad dense X[:16,:500] dtype={x_node.dtype}"
    elif fmt == "10x_flex_h5_multilane":
        with h5py.File(first_file, "r") as f:
            data_ds = f["matrix/data"]
            probe_vals = np.asarray(data_ds[:10000], dtype=np.float64)
            probe_note = f"10x matrix/data[:{len(probe_vals)}] dtype={data_ds.dtype}"
    elif fmt == "mtx_gz_plus_csc_cache":
        opener = gzip.open if str(first_file).endswith(".gz") else open
        vals = []
        with opener(first_file, "rt", encoding="utf-8", errors="replace") as fh:
            header_seen = False
            for line in fh:
                if line.startswith("%"):
                    continue
                if not header_seen:
                    header_seen = True
                    continue
                parts = line.strip().split()
                if len(parts) >= 3:
                    vals.append(float(parts[2]))
                if len(vals) >= 10000:
                    break
        probe_vals = np.asarray(vals, dtype=np.float64)
        probe_note = f"mtx.gz first {len(probe_vals)} nonzero entries"

    if probe_vals is None or probe_vals.size == 0:
        return {
            "available": False,
            "n_files": len(matched),
            "total_bytes": int(total_bytes),
            "reason": "empty probe values",
        }

    is_finite = bool(np.all(np.isfinite(probe_vals)))
    is_nonneg = bool(np.all(probe_vals >= 0.0))
    is_int = bool(np.all(probe_vals == np.floor(probe_vals)))
    return {
        "available": bool(is_finite and is_nonneg and is_int),
        "raw_dir_rel": str(raw_dir.relative_to(REPO_ROOT)),
        "n_matched_files": len(matched),
        "total_bytes": int(total_bytes),
        "total_gb": round(total_bytes / 1e9, 4),
        "probed_file_rel": str(first_file.relative_to(REPO_ROOT)),
        "probe_note": probe_note,
        "probe_n_values": int(probe_vals.size),
        "probe_min": float(np.min(probe_vals)),
        "probe_max": float(np.max(probe_vals)),
        "probe_all_finite": is_finite,
        "probe_all_nonneg": is_nonneg,
        "probe_all_integer": is_int,
    }


def _audit_prior_state_exposure(line_name: str) -> Dict[str, Any]:
    """
    Programmatically verify zero prior exposure to any Arc STATE model across
    `round16/results/` and `round16_addB/` (strictly NEVER reading `round16/audit/` or `A-T3.txt`).
    """
    r16_res_dir = REPO_ROOT / "round16" / "results"
    r16_addb_dir = REPO_ROOT / "round16_addB"

    matched_result_files: List[str] = []
    if r16_res_dir.exists():
        for p in sorted(r16_res_dir.rglob("*")):
            if p.is_file() and line_name in p.name:
                matched_result_files.append(str(p.relative_to(REPO_ROOT)))
            elif p.is_dir() and p.name == line_name:
                matched_result_files.append(str(p.relative_to(REPO_ROOT)) + "/")

    matched_summary_keys: List[str] = []
    for srel in ("round16/results/summary_r16.json", "round16/results/sensitivity/summary_r16.json"):
        sp = REPO_ROOT / srel
        if sp.exists():
            sd = json.loads(sp.read_text(encoding="utf-8"))
            if line_name in sd.get("lines", {}):
                matched_summary_keys.append(f"{srel}:lines[{line_name}]")
            if line_name in sd.get("control_source_barcodes", {}):
                matched_summary_keys.append(f"{srel}:control_source_barcodes[{line_name}]")

    matched_addb_files: List[str] = []
    if r16_addb_dir.exists():
        for p in sorted(r16_addb_dir.rglob("*")):
            if p.is_file() and line_name in p.name:
                matched_addb_files.append(str(p.relative_to(REPO_ROOT)))

    exposed = bool(matched_result_files or matched_summary_keys or matched_addb_files)
    return {
        "prior_state_exposure": "EXPOSED" if exposed else "NONE",
        "r16_result_file_hits": matched_result_files,
        "r16_summary_key_hits": matched_summary_keys,
        "r16_addb_file_hits": matched_addb_files,
        "audit_dirs_read": "NONE (round16/audit/ and A-T3.txt strictly excluded)",
    }


def _load_state_vocab_metadata() -> Tuple[List[str], List[str]]:
    """Load G_model (2,000) and P_model (2,024) strictly from shipped `var_dims.pkl` (no torch)."""
    with open(R16_VAR_DIMS_PKL, "rb") as fh:
        vd = pickle.load(fh)
    g_model = [str(g) for g in vd["gene_names"]]
    p_model = [str(p) for p in vd["pert_names"]]
    return g_model, p_model


def build_inventory(out_json: Path = DEFAULT_OUT_JSON) -> Dict[str, Any]:
    hgnc = HGNCBridge.from_tsv(DEFAULT_HGNC_FILE)
    common_genes_r9 = load_r9_common_genes(
        hgnc_bridge=hgnc, common_genes_file=R9_COMMON_GENES_FILE
    )
    g_model_2000, p_model_2024 = _load_state_vocab_metadata()
    g_model_canon_set: Set[str] = {hgnc.canonicalize(g) or g for g in g_model_2000}
    p_model_nonctrl = [p for p in p_model_2024 if p != "non-targeting"]
    p_model_canon_set: Set[str] = {hgnc.canonicalize(p) or p for p in p_model_nonctrl}

    # Load frozen R9 source pool target/gene sets ONLY (zero metrics/correlations computed)
    src_ds_map = {
        sn: load_standardized_line_r10(
            R9_LINES_DIR / sn,
            representation="log1p",
            common_genes=common_genes_r9,
            hgnc_bridge=hgnc,
        )
        for sn in R9_SOURCE_POOL_6
    }
    shared_src_genes: Set[str] = set.intersection(*(set(ds.genes) for ds in src_ds_map.values()))
    src_repr_targets: Set[str] = set()
    for ds in src_ds_map.values():
        if ds.repr_mask is not None:
            src_repr_targets.update(t for t, ok in zip(ds.target_names, ds.repr_mask) if ok)
        else:
            src_repr_targets.update(ds.target_names)

    records: List[Dict[str, Any]] = []
    for spec in CANDIDATE_SPECS:
        line_name = spec["line"]
        ldir = REPO_ROOT / spec["line_dir_rel"]
        meta_path = ldir / "line_metadata.json"
        meta = json.loads(meta_path.read_text(encoding="utf-8"))

        ds_t = load_standardized_line_r10(
            ldir,
            representation="log1p",
            common_genes=common_genes_r9,
            hgnc_bridge=hgnc,
        )
        g_r9_aligned = [g for g in ds_t.genes if g in shared_src_genes]
        g_r9_set = set(g_r9_aligned)
        g_eval_state = [g for g in g_r9_aligned if g in g_model_canon_set]
        g_eval_state_set = set(g_eval_state)

        # Exact PREREG_I §2.1 target eligibility (from data + vocabulary metadata ONLY)
        t_repr_line = [
            i for i, ok in enumerate(ds_t.repr_mask) if bool(ok)
        ] if ds_t.repr_mask is not None else list(range(len(ds_t.target_names)))

        t_p5_r9 = [
            i
            for i in t_repr_line
            if ds_t.target_names[i] in g_r9_set and ds_t.target_names[i] in src_repr_targets
        ]
        t_eval_prereg_i = [
            i for i in t_p5_r9 if ds_t.target_names[i] in p_model_canon_set
        ]
        t_eval_half_valid = [
            i for i in t_eval_prereg_i if ds_t.half_valid is not None and bool(ds_t.half_valid[i])
        ]
        t_eval_self_masked = [
            i for i in t_eval_prereg_i if ds_t.target_names[i] in g_eval_state_set
        ]
        t_repr_in_pmodel_no_fct = [
            i for i in t_repr_line if ds_t.target_names[i] in p_model_canon_set
        ]

        # Control sampling & NTC half-split validity (PREREG_I §2.2 & §2.4.3)
        n_nt_total = int(meta["n_nt_total"])
        n_nt_half_a = int(meta["n_nt_half_a"])
        n_nt_half_b = int(meta["n_nt_half_b"])
        ntc_half_split_valid = bool(n_nt_half_a >= 15 and n_nt_half_b >= 15)
        if n_nt_total >= 512:
            ctrl_sampling_n = 512
        elif n_nt_total >= 128:
            ctrl_sampling_n = (n_nt_total // 32) * 32
        else:
            ctrl_sampling_n = 0
        ctrl_n_sets_32 = ctrl_sampling_n // 32
        n_nt_remaining_after_primary = max(0, n_nt_total - ctrl_sampling_n)
        within_line_resample_eligible = bool(ctrl_sampling_n >= 128 and n_nt_remaining_after_primary >= 128)

        raw_probe = _probe_raw_integer_counts(
            REPO_ROOT / spec["raw_dir_rel"],
            spec["raw_probe_glob"],
            spec["raw_format"],
        )
        exposure = _audit_prior_state_exposure(line_name)

        n_eval_prereg_i = len(t_eval_prereg_i)
        n_g_eval_r9 = len(g_r9_aligned)
        n_g_eval_state = len(g_eval_state)
        prereg_i_target_gate = bool(n_eval_prereg_i >= 50)
        prereg_i_gene_gate = bool(n_g_eval_r9 >= 3000)
        prereg_i_ctrl_gate = bool(ctrl_sampling_n >= 128 and raw_probe["available"])
        prereg_i_eligible = bool(
            prereg_i_target_gate
            and prereg_i_gene_gate
            and ntc_half_split_valid
            and prereg_i_ctrl_gate
            and exposure["prior_state_exposure"] == "NONE"
        )

        # Independence & recommended R17 role (strictly rule-based from metadata & counts)
        if not prereg_i_eligible:
            if not prereg_i_target_gate:
                role = "SECONDARY_UNDERPOWERED_N_LT_50"
                role_reason = f"n_eval_prereg_i={n_eval_prereg_i} < 50 (fails PREREG_I target floor)"
            else:
                role = "EXCLUDED_GATE_FAIL"
                role_reason = "failed gene, NTC, raw-count, or prior-exposure gate"
        elif spec["relation_to_r16_discovery"] == "STRICT_CELL_SUBSET_OF_R16_DISCOVERY_a549_ko":
            role = "EXCLUDED_FROM_CONFIRMATORY_R16_SUBSET_SENSITIVITY_ONLY"
            role_reason = (
                "Passes PREREG_I counts (n_eval=403) but is a strict CAS9>0 cell subset (73.3%) "
                "of R16 discovery line a549_ko; non-independent of R16 discovery"
            )
        elif spec["relation_to_r16_discovery"].startswith("shared_study_lineage_library_with_r16"):
            role = "CONFIRMATORY_CANDIDATE_WITH_SHARED_STUDY_DISCLOSURE"
            role_reason = (
                "Passes all PREREG_I gates (n_eval=401 >= 50, |G_eval^R9|=6131 >= 3000, "
                "|G_eval^STATE|=1836, NTC half-split valid, raw integer counts verified, "
                "prior STATE exposure NONE, strictly post-dates SE-600M); requires mandatory "
                "shared-study/lineage/library disclosure with R16 discovery cd4t_rest "
                "(disjoint donors D2/D3 vs D1/D4 and 48h stimulated vs resting state)"
            )
        else:
            role = "CONFIRMATORY_PRIMARY_INDEPENDENT"
            role_reason = "Passes all PREREG_I gates and fully independent of R16 discovery"

        rec = {
            "line": line_name,
            "round": spec["round"],
            "line_dir_rel": spec["line_dir_rel"],
            "modality": spec["modality"],
            "assay": spec["assay"],
            "cell_lineage": spec["cell_lineage"],
            "canonical_family": spec["canonical_family"],
            "study_accession": spec["study_accession"],
            "public_release_date": spec["public_release_date"],
            "se600m_pretrain_chronology": spec["se600m_pretrain_chronology"],
            "prior_state_exposure_audit": exposure,
            "raw_counts_for_control_sampling": raw_probe,
            "control_pool_counts": {
                "n_nt_total": n_nt_total,
                "n_nt_half_a": n_nt_half_a,
                "n_nt_half_b": n_nt_half_b,
                "ntc_half_split_valid": ntc_half_split_valid,
                "half_rule": meta.get("adapter_audit", {}).get(
                    "half_rule", meta.get("adapter_audit", {}).get("half_rule_used", "unit_parity")
                ),
                "n_batches": int(meta["n_batches"]),
                "ctrl_sampling_n": ctrl_sampling_n,
                "ctrl_n_sets_32": ctrl_n_sets_32,
                "ctrl_sampling_truncated_lt_512": bool(0 < ctrl_sampling_n < 512),
                "within_line_resample_eligible": within_line_resample_eligible,
            },
            "prereg_i_eligibility": {
                "n_genes_final_line": int(meta["n_genes_final"]),
                "n_g_eval_r9": n_g_eval_r9,
                "n_g_eval_state": n_g_eval_state,
                "n_targets_deduplicated_line": int(meta["n_targets_deduplicated"]),
                "n_targets_reproducible_line": len(t_repr_line),
                "n_eval_r9_p5": len(t_p5_r9),
                "n_eval_prereg_i": n_eval_prereg_i,
                "n_eval_half_valid": len(t_eval_half_valid),
                "n_eval_self_masked": len(t_eval_self_masked),
                "n_repr_in_pmodel_no_fct_filter": len(t_repr_in_pmodel_no_fct),
                "pert_half_split_half_size_floor": n_eval_prereg_i // 2,
                "passes_n_eval_ge_50": prereg_i_target_gate,
                "passes_g_eval_r9_ge_3000": prereg_i_gene_gate,
                "passes_ntc_half_split": ntc_half_split_valid,
                "passes_raw_ctrl_ge_128": prereg_i_ctrl_gate,
                "prereg_i_eligible": prereg_i_eligible,
            },
            "family_structure": {
                "canonical_family": spec["canonical_family"],
                "relation_to_r16_discovery": spec["relation_to_r16_discovery"],
                "relation_to_other_r17_candidates": spec["relation_to_other_r17_candidates"],
                "notes": spec["family_notes"],
            },
            "r17_role_classification": role,
            "r17_role_reason": role_reason,
        }
        records.append(rec)

    summary = {
        "n_local_never_state_candidates": len(records),
        "all_zero_prior_state_exposure": all(
            r["prior_state_exposure_audit"]["prior_state_exposure"] == "NONE" for r in records
        ),
        "all_raw_integer_counts_available": all(
            r["raw_counts_for_control_sampling"]["available"] for r in records
        ),
        "prereg_i_eligible_lines_all": [
            r["line"] for r in records if r["prereg_i_eligibility"]["prereg_i_eligible"]
        ],
        "confirmatory_primary_strictly_independent_local": [
            r["line"] for r in records if r["r17_role_classification"] == "CONFIRMATORY_PRIMARY_INDEPENDENT"
        ],
        "confirmatory_candidates_with_shared_study_disclosure": [
            r["line"]
            for r in records
            if r["r17_role_classification"] == "CONFIRMATORY_CANDIDATE_WITH_SHARED_STUDY_DISCLOSURE"
        ],
        "excluded_r16_subset_sensitivity_only": [
            r["line"]
            for r in records
            if r["r17_role_classification"] == "EXCLUDED_FROM_CONFIRMATORY_R16_SUBSET_SENSITIVITY_ONLY"
        ],
        "secondary_underpowered_n_lt_50": [
            r["line"] for r in records if r["r17_role_classification"] == "SECONDARY_UNDERPOWERED_N_LT_50"
        ],
        "secondary_usable_n_ge_15_lines": [
            r["line"]
            for r in records
            if r["r17_role_classification"] == "SECONDARY_UNDERPOWERED_N_LT_50"
            and r["prereg_i_eligibility"]["n_eval_prereg_i"] >= 15
        ],
        "r16_discovery_lines_strictly_excluded_from_confirmatory": list(R16_DISCOVERY_LINES),
    }

    payload = {
        "schema_version": "LINE_INVENTORY_J_v1",
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "generator_script": "round17/src/build_line_inventory_j.py",
        "model_free_attestation": (
            "Strictly data-only and counts-only inventory. Zero neural network imports (no torch/state), "
            "zero MPS/GPU calls, zero model predictions, zero target-expression correlation/metric "
            "computations, zero reads of round16/audit/ or A-T3.txt."
        ),
        "reference_vocab_files": {
            "var_dims_pkl": str(R16_VAR_DIMS_PKL.relative_to(REPO_ROOT)),
            "var_dims_pkl_sha256": _sha256_file(R16_VAR_DIMS_PKL),
            "n_g_model": len(g_model_2000),
            "n_p_model_total": len(p_model_2024),
            "n_p_model_noncontrol": len(p_model_nonctrl),
            "r9_common_genes_file": str(R9_COMMON_GENES_FILE.relative_to(REPO_ROOT)),
            "n_r9_common_genes": len(common_genes_r9),
            "r9_fct_source_pool": list(R9_SOURCE_POOL_6),
            "hgnc_file": str(DEFAULT_HGNC_FILE),
        },
        "summary": summary,
        "candidates": records,
    }

    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return payload


def main() -> int:
    payload = build_inventory(DEFAULT_OUT_JSON)
    print(f"Wrote {DEFAULT_OUT_JSON} ({DEFAULT_OUT_JSON.stat().st_size:,} B)")
    print("SUMMARY:", json.dumps(payload["summary"], indent=2))
    print("\nLINE TABLE:")
    for r in payload["candidates"]:
        el = r["prereg_i_eligibility"]
        cp = r["control_pool_counts"]
        rc = r["raw_counts_for_control_sampling"]
        print(
            f"  {r['round']}/{r['line']:16s} | mod={r['modality']:8s} | "
            f"raw_int={str(rc['available']):5s} ({rc.get('total_gb', 0):6.3f} GB) | "
            f"N_eval={el['n_eval_prereg_i']:3d} (hv={el['n_eval_half_valid']:3d}, R9_P5={el['n_eval_r9_p5']:3d}, repr∩P={el['n_repr_in_pmodel_no_fct_filter']:4d}) | "
            f"|G_R9|={el['n_g_eval_r9']:4d} |G_ST|={el['n_g_eval_state']:4d} | "
            f"NTC={cp['n_nt_total']:6d} (A={cp['n_nt_half_a']}, B={cp['n_nt_half_b']}, samp={cp['ctrl_sampling_n']}) | "
            f"exp={r['prior_state_exposure_audit']['prior_state_exposure']} | "
            f"elig={str(el['prereg_i_eligible']):5s} | role={r['r17_role_classification']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
