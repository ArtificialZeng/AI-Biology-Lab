#!/usr/bin/env python3
"""
Unit & Integration Tests for Round 13 (`PREREG_F`, `GSE327644`) Modules (`round13/src/`).

Runs cleanly exclusively on synthetic Cell Ranger Flex `.h5` / companion fixtures without
touching any real Round 13 data or modifying any file under `round10/`, `round11/`, or `round12/`.
"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

R13_SRC_DIR = Path(__file__).resolve().parent
REPO_ROOT = R13_SRC_DIR.parents[1]
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
R10_SCRIPTS_DIR = REPO_ROOT / "round10" / "shared" / "scripts"
R11_SRC_DIR = REPO_ROOT / "round11" / "src"
R12_SRC_DIR = REPO_ROOT / "round12" / "src"
for p in (str(R10_SRC_DIR), str(R10_SCRIPTS_DIR), str(R11_SRC_DIR), str(R12_SRC_DIR), str(R13_SRC_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from data_loader_r10 import load_standardized_line_r10  # noqa: E402
from prepare_r10_lines import HGNCLookup, build_hgnc_lookup  # noqa: E402
from test_unit_r11 import _make_planted_synthetic_cohort  # noqa: E402
from run_eval_r12 import evaluate_target_against_sources_r12  # noqa: E402

from adjudicate_r13 import adjudicate_r13_results, render_r13_markdown_report  # noqa: E402
from data_loader_r13 import DEFAULT_R13_SEED, check_prereg_f_gate  # noqa: E402
from dataset_readers_r13 import (  # noqa: E402
    build_gse327644_parts,
    check_gse327644_feasibility,
    load_gse327644_guide_metadata,
    prepare_r13_line,
)

_CACHED_HGNC: Optional[HGNCLookup] = None


def _get_hgnc() -> HGNCLookup:
    global _CACHED_HGNC
    if _CACHED_HGNC is None:
        _CACHED_HGNC = build_hgnc_lookup()
    return _CACHED_HGNC


def _pick_hgnc_genes(n_genes: int) -> List[str]:
    hgnc = _get_hgnc()
    picked = [g for g in sorted(hgnc.approved) if g and not g.startswith("MT-")][:n_genes]
    if len(picked) < n_genes:
        raise RuntimeError(f"Requested {n_genes} HGNC genes, found {len(picked)}")
    return picked


def _write_flex_10x_h5(
    path: Path,
    mat_csc: sp.csc_matrix,
    barcodes: Sequence[str],
    feat_ids: Sequence[str],
    feat_names: Sequence[str],
    feat_types: Sequence[str],
    float_data: bool = False,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        mg = f.create_group("matrix")
        mg.create_dataset("shape", data=np.array(mat_csc.shape, dtype=np.int64))
        mg.create_dataset("barcodes", data=np.array([str(x).encode("utf-8") for x in barcodes]))
        dt = np.float32 if float_data else np.int32
        mg.create_dataset("data", data=mat_csc.data.astype(dt))
        mg.create_dataset("indices", data=mat_csc.indices.astype(np.int32))
        mg.create_dataset("indptr", data=mat_csc.indptr.astype(np.int64))
        fg = mg.create_group("features")
        fg.create_dataset("id", data=np.array([str(x).encode("utf-8") for x in feat_ids]))
        fg.create_dataset("name", data=np.array([str(x).encode("utf-8") for x in feat_names]))
        fg.create_dataset("feature_type", data=np.array([str(x).encode("utf-8") for x in feat_types]))
        fg.create_dataset("genome", data=np.array([b"GRCh38"] * len(feat_ids)))


def _make_synthetic_flex_matrix(
    rng: np.random.Generator,
    genes: Sequence[str],
    guide_specs: Sequence[tuple],  # (guide_id, target_row_or_neg1)
    guide_cell_counts: Dict[str, int],
    prefix: str,
    low_umi_cells: int = 2,
    doublet_cells: int = 2,
    mt_gene: bool = False,
    high_mt_fraction: bool = False,
) -> tuple:
    """Create synthetic Flex (features x cells) integer count matrix with known QC / doublets."""
    gene_list = list(genes) + (["MT-CO1"] if mt_gene else [])
    n_rna = len(gene_list)
    guide_ids = [g[0] for g in guide_specs]
    feat_names = gene_list + guide_ids
    feat_ids = feat_names[:]
    feat_types = ["Gene Expression"] * n_rna + ["CRISPR Guide Capture"] * len(guide_ids)
    n_feat = len(feat_names)

    g_row_map = {gid: n_rna + k for k, gid in enumerate(guide_ids)}
    g_tgt_map = {gid: tgt for gid, tgt in guide_specs}

    bcs: List[str] = []
    cols_meta: List[tuple] = []  # (primary_row_or_-1, second_row_or_-1, tgt_idx, low_umi)
    cid = 0
    for gid, n_c in guide_cell_counts.items():
        grow = g_row_map[gid]
        trow = g_tgt_map[gid]
        for _ in range(n_c):
            bcs.append(f"{prefix}_C{cid:04d}-1")
            cols_meta.append((grow, -1, trow, False))
            cid += 1
    for _ in range(low_umi_cells):
        bcs.append(f"{prefix}_LOW{cid:04d}-1")
        cols_meta.append((g_row_map[guide_ids[0]], -1, -1, True))
        cid += 1
    for _ in range(doublet_cells):
        bcs.append(f"{prefix}_DBL{cid:04d}-1")
        cols_meta.append((g_row_map[guide_ids[0]], g_row_map[guide_ids[-1]], g_tgt_map[guide_ids[-1]], False))
        cid += 1

    mat = np.zeros((n_feat, len(bcs)), dtype=np.int32)
    for j, (gr1, gr2, trow, is_low) in enumerate(cols_meta):
        if is_low:
            mat[:40, j] = 1  # < 300 UMI & < 200 genes
        else:
            mat[: len(genes), j] = rng.integers(2, 5, size=len(genes), dtype=np.int32)  # >= 500 UMI & >= 250 genes
            if mt_gene and high_mt_fraction:
                mat[len(genes), j] = 600  # ~40-50% MT (must NOT be excluded when max_pct_mt=None)
        if j % 3 == 0:
            mat[n_rna:n_feat, j] = 1
        mat[gr1, j] = 24
        if gr2 >= 0:
            mat[gr2, j] = 22  # doublet: n_pos == 2 -> excluded
        if trow >= 0 and not is_low:
            mat[:15, j] += 6
            mat[trow, j] = 1
    return sp.csc_matrix(mat), bcs, feat_ids, feat_names, feat_types


class Round13UnitTests(unittest.TestCase):
    def test_1_blinding_gate_and_adjudicator_r13(self) -> None:
        self.assertEqual(DEFAULT_R13_SEED, 20261223)
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)
            with self.assertRaises(RuntimeError):
                check_prereg_f_gate(p)
            (p / "PREREG_F_COMMIT.txt").write_text("prereg_f_111\n", encoding="utf-8")
            with self.assertRaises(RuntimeError):
                check_prereg_f_gate(p)
            (p / "CODE_FREEZE_R13.txt").write_text("freeze_f_222\n", encoding="utf-8")
            pr, fr = check_prereg_f_gate(p)
            self.assertEqual(pr, "prereg_f_111")
            self.assertEqual(fr, "freeze_f_222")

        src6, tgt_syn = _make_planted_synthetic_cohort(seed=131, n_genes=150, n_targets=55, n_sources=6)
        tgt_ctrl = {**tgt_syn.__dict__, "name": "mac_ctrl"}
        from data_loader_r10 import LineDatasetR10  # noqa: E402

        ds_ctrl = LineDatasetR10(**tgt_ctrl)
        ds_lps = LineDatasetR10(**{**tgt_syn.__dict__, "name": "mac_lps"})
        l1, b1, ex1 = evaluate_target_against_sources_r12(
            ds_ctrl, src6, include_r9_anchors=False, min_gene_overlap=100, n_boot=120, seed=DEFAULT_R13_SEED, n_perm=10
        )
        l2, b2, ex2 = evaluate_target_against_sources_r12(
            ds_lps, src6, include_r9_anchors=False, min_gene_overlap=100, n_boot=120, seed=DEFAULT_R13_SEED, n_perm=10
        )
        df_l = pd.DataFrame(l1 + l2)
        df_b = pd.DataFrame(b1 + b2)
        summary = {
            "phase": "holdout_r13",
            "source_mode": "r9_6",
            "source_lines": sorted(src6.keys()),
            "prereg_commit_hash": "prereg_f",
            "code_freeze_hash": "freeze_r13",
            "bootstrap_seed": DEFAULT_R13_SEED,
            "n_boot": 120,
            "fct_perm100": {"mac_ctrl": ex1["perm100"], "mac_lps": ex2["perm100"]},
            "n5_cone_metrics": {"mac_ctrl": ex1["n5_cone_metrics"], "mac_lps": ex2["n5_cone_metrics"]},
            "n5_cone_split_metrics": {"mac_ctrl": ex1["n5_cone_split_metrics"], "mac_lps": ex2["n5_cone_split_metrics"]},
            "n5_perm100": {"mac_ctrl": ex1["n5_perm100"], "mac_lps": ex2["n5_perm100"]},
        }

        # Guard: mac_lps CANNOT be confirmatory under PREREG_F §2
        with self.assertRaises(ValueError):
            adjudicate_r13_results(df_l, df_b, summary, confirmatory_lines=["mac_lps"])
        with self.assertRaises(ValueError):
            adjudicate_r13_results(df_l, df_b, summary, confirmatory_lines=["mac_ctrl", "mac_lps"])

        adj0 = adjudicate_r13_results(df_l, df_b, summary, confirmatory_lines=[])
        self.assertEqual(adj0["confirmatory_status"], "NO_CONFIRMATORY_LINES_DESCRIPTIVE_ONLY")
        adj1 = adjudicate_r13_results(df_l, df_b, summary, confirmatory_lines=["mac_ctrl"])
        self.assertEqual(adj1["confirmatory_status"], "CONFIRMATORY_SINGLE_LINE_K1")
        self.assertTrue(adj1["no_cross_round_pooling_notice"])
        self.assertEqual(adj1["secondary_lines"], ["mac_lps"])
        md = render_r13_markdown_report(df_l, df_b, summary, adj1, "unit_r13")
        self.assertIn("PREREG_F", md)
        self.assertIn("20261223", md)

    def test_2_gse327644_flex_h5_channel_parity_safetarget_and_process_r10(self) -> None:
        rng = np.random.default_rng(1302)
        hgnc = _get_hgnc()
        genes = _pick_hgnc_genes(260)  # Flex panel without MT genes (>= 200 genes, >= 300 UMI)
        g0, g1 = genes[0], genes[1]
        guide_specs = [
            ("NTC_g1", -1),
            ("safe-targeting_g1", -1),
            (f"{g0}_g1", 0),
            (f"{g1}_g1", 1),
        ]

        with tempfile.TemporaryDirectory() as td:
            raw_dir = Path(td) / "raw"
            gdir = raw_dir / "gse327644"
            gdir.mkdir(parents=True)

            # Guide template distinguishing non-targeting vs safe-targeting -> ONLY non-targeting = NTC
            pd.DataFrame(
                [
                    {"guide_id": "NTC_g1", "name": "NTC_g1", "guide_type": "non-targeting", "target_gene_name": ""},
                    {"guide_id": "safe-targeting_g1", "name": "safe-targeting_g1", "guide_type": "safe-targeting", "target_gene_name": ""},
                    {"guide_id": f"{g0}_g1", "name": f"{g0}_g1", "guide_type": "targeting", "target_gene_name": g0},
                    {"guide_id": f"{g1}_g1", "name": f"{g1}_g1", "guide_type": "targeting", "target_gene_name": g1},
                ]
            ).to_csv(gdir / "GSE327644_guide_template.txt.gz", sep="\t", index=False, compression="gzip")

            # L01(20 NTC), L02(18 NTC), L03(6 NTC < 15 -> falls back to half_1_pooled!), L04(18 NTC)
            lane_ntc = {"L01": 20, "L02": 18, "L03": 6, "L04": 18}
            for ln, n_nt in lane_ntc.items():
                counts = {"NTC_g1": n_nt, "safe-targeting_g1": 4, f"{g0}_g1": 10, f"{g1}_g1": 10}
                csc, bcs, fids, fnames, ftypes = _make_synthetic_flex_matrix(
                    rng, genes, guide_specs, counts, prefix=f"CTRL_{ln}"
                )
                _write_flex_10x_h5(
                    gdir / f"GSM900001_Mac_{ln}_Ctrl_sample_filtered_feature_bc_matrix.h5",
                    csc, bcs, fids, fnames, ftypes,
                )
                # Also write a raw_feature_bc_matrix (must be ignored) & an LPS file (must be excluded for mac_ctrl)
                _write_flex_10x_h5(
                    gdir / f"GSM900001_Mac_{ln}_Ctrl_raw_feature_bc_matrix.h5",
                    csc, bcs, fids, fnames, ftypes, float_data=True,
                )
                csc_lps, bcs_lps, *_ = _make_synthetic_flex_matrix(
                    rng, genes, guide_specs, counts, prefix=f"LPS_{ln}"
                )
                _write_flex_10x_h5(
                    gdir / f"GSM900002_Mac_{ln}_LPS_sample_filtered_feature_bc_matrix.h5",
                    csc_lps, bcs_lps, fids, fnames, ftypes,
                )

            feas = check_gse327644_feasibility(raw_dir, line="mac_ctrl", hgnc=hgnc)
            self.assertTrue(feas["feasible"], feas.get("reason"))
            self.assertEqual(feas["n_h5_files"], 4)
            self.assertEqual(feas["guide_call_source"], "h5_crispr_umi_em")
            self.assertTrue(feas["has_safe_targeting"])
            self.assertEqual(feas["ntc_rule_used"], "only_nontargeting_excluding_safetargeting")
            self.assertEqual(feas["lanes_detected"], ["L01", "L02", "L03", "L04"])

            parts, audit = build_gse327644_parts(raw_dir, hgnc, line="mac_ctrl")
            self.assertFalse(audit["used_seed42_fallback"])
            self.assertEqual(audit["n_safe_targeting_excluded"], 4 * 4)
            self.assertEqual(audit["n_rna_qc_excluded"], 4 * 2)  # 2 low-UMI cells per channel excluded in RNA QC
            self.assertEqual(audit["n_multi_or_zero_guide_excluded"], 4 * 2)  # 2 doublets per channel excluded in guide calling
            self.assertEqual(audit["n_kept_half_a_ntc"], 20 + 6)
            self.assertEqual(audit["n_kept_half_b_ntc"], 18 + 18)
            self.assertEqual(set(audit["channels_used"]), {"L01", "L02", "half_1_pooled", "L04"})
            self.assertEqual(parts[0].half_flags[0], 1)  # L01 -> half 1
            self.assertEqual(parts[1].half_flags[0], 2)  # L02 -> half 2
            self.assertEqual(parts[2].batches[0], "half_1_pooled")  # L03 (<15 NTC) -> half_1_pooled
            self.assertEqual(parts[3].half_flags[0], 2)  # L04 -> half 2

            # Integrate end-to-end with frozen process_r10_line -> identical standardized schema
            out_root = Path(td) / "lines"
            prepare_r13_line("mac_ctrl", raw_dir, out_root, hgnc=hgnc, workers=1)
            ds = load_standardized_line_r10(out_root / "mac_ctrl", min_gene_overlap=50)
            self.assertEqual(len(ds.target_names), 2)
            self.assertEqual(len(ds.genes), len(genes))

    def test_3_gse327644_companion_protospacer_csv_and_all_controls_ntc(self) -> None:
        rng = np.random.default_rng(1303)
        hgnc = _get_hgnc()
        genes = _pick_hgnc_genes(260)
        g0, g1 = genes[0], genes[1]
        # No safe-targeting distinction in template -> ALL control guides count as NTC
        guide_specs = [("Ctrl_1", -1), ("Control_2", -1), (f"{g0}_sg", 0), (f"{g1}_sg", 1)]

        with tempfile.TemporaryDirectory() as td:
            raw_dir = Path(td) / "raw"
            gdir = raw_dir / "gse327644"
            gdir.mkdir(parents=True)
            pd.DataFrame(
                [
                    {"guide_id": "Ctrl_1", "guide_type": "control", "target_gene_name": "Control"},
                    {"guide_id": "Control_2", "guide_type": "control", "target_gene_name": "Control"},
                    {"guide_id": f"{g0}_sg", "guide_type": "targeting", "target_gene_name": g0},
                    {"guide_id": f"{g1}_sg", "guide_type": "targeting", "target_gene_name": g1},
                ]
            ).to_csv(gdir / "GSE327644_guide_template.txt.gz", sep="\t", index=False, compression="gzip")

            for ln, lps_name in (("L01", "Mac_L01_LPS"), ("L02", "Mac_L02_LPS")):
                counts = {"Ctrl_1": 10, "Control_2": 8, f"{g0}_sg": 8, f"{g1}_sg": 8}
                csc, bcs, fids, fnames, ftypes = _make_synthetic_flex_matrix(
                    rng, genes, guide_specs, counts, prefix=ln, low_umi_cells=0, doublet_cells=0
                )
                sub = gdir / lps_name
                h5_p = sub / f"{lps_name}_sample_filtered_feature_bc_matrix.h5"
                _write_flex_10x_h5(h5_p, csc, bcs, fids, fnames, ftypes)

                # Companion protospacer_calls_per_cell.csv assigning 1 multiplet & 1 unassigned to exclude
                rows_csv = []
                idx_cursor = 0
                for gid, n_c in counts.items():
                    for k in range(n_c):
                        if gid == f"{g1}_sg" and k == 0:
                            rows_csv.append({"cell_barcode": bcs[idx_cursor], "num_features": 2, "feature_call": f"{g0}_sg|{g1}_sg", "num_umis": "15|12"})
                        elif gid == f"{g0}_sg" and k == 0:
                            rows_csv.append({"cell_barcode": bcs[idx_cursor], "num_features": 0, "feature_call": "None", "num_umis": "0"})
                        else:
                            rows_csv.append({"cell_barcode": bcs[idx_cursor], "num_features": 1, "feature_call": gid, "num_umis": "25"})
                        idx_cursor += 1
                pd.DataFrame(rows_csv).to_csv(sub / "protospacer_calls_per_cell.csv", index=False)

            feas = check_gse327644_feasibility(raw_dir, line="mac_lps", hgnc=hgnc)
            self.assertTrue(feas["feasible"])
            self.assertFalse(feas["has_safe_targeting"])
            self.assertEqual(feas["ntc_rule_used"], "all_controls_no_safetargeting_distinction")
            self.assertEqual(feas["guide_call_source"], "protospacer_calls_per_cell_csv")

            parts, audit = build_gse327644_parts(raw_dir, hgnc, line="mac_lps")
            self.assertEqual(audit["guide_call_mode"], "protospacer_calls_per_cell_csv")
            self.assertEqual(audit["n_kept_ntc"], 2 * 18)  # BOTH Ctrl_1 & Control_2 counted as NTC
            self.assertEqual(audit["n_multi_or_zero_guide_excluded"], 2 * 2)  # multiplet + None per lane
            self.assertEqual(audit["n_kept_half_a_ntc"], 18)
            self.assertEqual(audit["n_kept_half_b_ntc"], 18)

    def test_4_gse327644_stage2_seed42_high_mt_retention_and_feasibility_failures(self) -> None:
        rng = np.random.default_rng(1304)
        hgnc = _get_hgnc()
        genes = _pick_hgnc_genes(260)
        g0 = genes[0]
        guide_specs = [("NTC_1", -1), (f"{g0}_1", 0)]

        with tempfile.TemporaryDirectory() as td:
            raw_dir = Path(td) / "raw"
            gdir = raw_dir / "gse327644"
            gdir.mkdir(parents=True)
            pd.DataFrame(
                [
                    {"guide_id": "NTC_1", "guide_type": "non-targeting", "target_gene_name": ""},
                    {"guide_id": f"{g0}_1", "guide_type": "targeting", "target_gene_name": g0},
                ]
            ).to_csv(gdir / "GSE327644_guide_template.txt.gz", sep="\t", index=False, compression="gzip")

            # Single file with NO channel/lane tag + ~45% mitochondrial fraction -> MUST use seed-42 & retain high-MT cells
            counts = {"NTC_1": 24, f"{g0}_1": 16}
            csc, bcs, fids, fnames, ftypes = _make_synthetic_flex_matrix(
                rng, genes, guide_specs, counts, prefix="NOLANE", low_umi_cells=0, doublet_cells=0,
                mt_gene=True, high_mt_fraction=True,
            )
            h5_p = gdir / "Mac_Ctrl_filtered_feature_bc_matrix.h5"
            _write_flex_10x_h5(h5_p, csc, bcs, fids, fnames, ftypes)

            parts, audit = build_gse327644_parts(raw_dir, hgnc, line="mac_ctrl")
            self.assertTrue(audit["used_seed42_fallback"])
            self.assertEqual(audit["half_rule_used"], "seed42_stratified_50_50")
            self.assertIsNone(audit["max_pct_mt"])
            self.assertEqual(audit["n_qc_kept_cells"], 40)  # all 40 high-MT cells kept!
            self.assertTrue(np.array_equal(parts[0].batches, np.where(parts[0].half_flags == 1, "half_1", "half_2")))

            # Float normalized counts in /matrix/data MUST fail feasibility -> EXCLUDED_FEASIBILITY
            _write_flex_10x_h5(h5_p, csc, bcs, fids, fnames, ftypes, float_data=False)
            with h5py.File(h5_p, "r+") as f:
                d = f["matrix/data"][:].astype(np.float32) + 0.25
                del f["matrix/data"]
                f["matrix"].create_dataset("data", data=d)
            feas_bad = check_gse327644_feasibility(raw_dir, line="mac_ctrl", hgnc=hgnc)
            self.assertFalse(feas_bad["feasible"])
            with self.assertRaisesRegex(RuntimeError, "EXCLUDED_FEASIBILITY"):
                build_gse327644_parts(raw_dir, hgnc, line="mac_ctrl")

    def test_5_prereg_f_target_symbols_hla_ppan_p2ry11_jun_and_ngenes200_threshold_fit(self) -> None:
        rng = np.random.default_rng(1305)
        hgnc = _get_hgnc()
        base_genes = [g for g in _pick_hgnc_genes(260) if g not in {"HLA-A", "PPAN-P2RY11", "PPAN", "JUN"}]
        genes = base_genes[:256] + ["HLA-A", "PPAN-P2RY11", "PPAN", "JUN"]
        idx_hla, idx_pp, idx_ppan, idx_jun = 256, 257, 258, 259
        guide_specs = [
            ("NTC-001", -1),
            ("HLA-A-1", idx_hla),
            ("PPAN-P2RY11-1", idx_pp),
            ("PPAN-1", idx_ppan),
            ("1-Jun", idx_jun),
        ]

        with tempfile.TemporaryDirectory() as td:
            raw_dir = Path(td) / "raw"
            gdir = raw_dir / "gse327644"
            gdir.mkdir(parents=True)

            # Exact PREREG_F §3.3/§3.4 columns + truncated artifacts & Excel 1-Jun artifact
            pd.DataFrame(
                [
                    {
                        "guide_id": "NTC-001",
                        "guide_protospacer": "AAAA",
                        "guide_role": "Control",
                        "guide_target_gene_id": "n/a",
                        "guide_PAM": "NGG",
                        "guide_target_gene_name": "n/a",
                    },
                    {
                        "guide_id": "HLA-A-1",
                        "guide_protospacer": "CCCC",
                        "guide_role": "Targeting a Gene",
                        "guide_target_gene_id": "ENSG00000206503",
                        "guide_PAM": "NGG",
                        "guide_target_gene_name": "HLA",
                    },
                    {
                        "guide_id": "PPAN-P2RY11-1",
                        "guide_protospacer": "GGGG",
                        "guide_role": "Targeting a Gene",
                        "guide_target_gene_id": "ENSG00000272031",
                        "guide_PAM": "NGG",
                        "guide_target_gene_name": "PPAN",
                    },
                    {
                        "guide_id": "PPAN-1",
                        "guide_protospacer": "TTTT",
                        "guide_role": "Targeting a Gene",
                        "guide_target_gene_id": "ENSG00000130810",
                        "guide_PAM": "NGG",
                        "guide_target_gene_name": "PPAN",
                    },
                    {
                        "guide_id": "1-Jun",
                        "guide_protospacer": "ACGT",
                        "guide_role": "Targeting a Gene",
                        "guide_target_gene_id": "ENSG00000177606",
                        "guide_PAM": "NGG",
                        "guide_target_gene_name": "JUN",
                    },
                ]
            ).to_csv(gdir / "GSE327644_guide_template.txt.gz", sep="\t", index=False, compression="gzip")

            # Build 2 channels with <= 12 ok_qc cells each per guide (< 50 total in ok_qc), PLUS
            # 60 cells per channel having ngenes = 220 >= 200 and umi = 220 < 300 so bimodal EM
            # fits on `ngenes >= 200` (PREREG_F §3.2) without changing final RNA QC cell retention.
            for ln in ("L01", "L02"):
                qc_counts = {"NTC-001": 6, "HLA-A-1": 4, "PPAN-P2RY11-1": 3, "PPAN-1": 2, "1-Jun": 4}
                csc_base, bcs_base, fids, fnames, ftypes = _make_synthetic_flex_matrix(
                    rng, genes, guide_specs, qc_counts, prefix=ln, low_umi_cells=0, doublet_cells=0
                )
                n_rna = len(genes)
                n_fit_only = 60
                fit_mat = np.zeros((len(fids), n_fit_only), dtype=np.int32)
                bcs_fit = []
                for j in range(n_fit_only):
                    bcs_fit.append(f"{ln}_FIT200_{j:03d}-1")
                    fit_mat[:220, j] = 1  # ngenes = 220 >= 200, total RNA UMI = 220 < 300
                    for g_row in range(n_rna, len(fids)):
                        fit_mat[g_row, j] = 1 if j < 42 else 28  # bimodal UMI across guides (>= 50 entries)
                full_mat = np.hstack([csc_base.toarray(), fit_mat])
                h5_p = gdir / f"GSM900099_Mac_{ln}_Ctrl_sample_filtered_feature_bc_matrix.h5"
                _write_flex_10x_h5(
                    h5_p,
                    sp.csc_matrix(full_mat),
                    list(bcs_base) + bcs_fit,
                    fids,
                    fnames,
                    ftypes,
                )

            # Verify .h5-only fallback also resolves 1-Jun / HLA-A-1 / PPAN-P2RY11-1 -> JUN / HLA-A / PPAN-P2RY11
            h5_sample = gdir / "GSM900099_Mac_L01_Ctrl_sample_filtered_feature_bc_matrix.h5"
            meta_h5 = load_gse327644_guide_metadata(None, h5_sample, hgnc)
            self.assertEqual(meta_h5["guide_to_target"]["HLA-A-1"][1], "HLA-A")
            self.assertEqual(meta_h5["guide_to_target"]["PPAN-P2RY11-1"][1], "PPAN-P2RY11")
            self.assertEqual(meta_h5["guide_to_target"]["PPAN-1"][1], "PPAN")
            self.assertEqual(meta_h5["guide_to_target"]["1-Jun"][1], "JUN")
            self.assertEqual(meta_h5["guide_to_target"]["JUN-1"][1], "JUN")

            parts, audit = build_gse327644_parts(raw_dir, hgnc, line="mac_ctrl")
            # Because fit_hist pools ngenes >= 200 (120+ entries/guide vs < 50 in ok_qc),
            # EM threshold settles at 6.0 rather than falling back to default_thr = 7.0.
            self.assertNotEqual(audit["median_guide_threshold"], 7.0)
            self.assertEqual(audit["median_guide_threshold"], 6.0)
            self.assertEqual(audit["n_qc_kept_cells"], 2 * 19)
            self.assertEqual(audit["n_rna_qc_excluded"], 2 * 60)

            kept_pkeys = [pk for p in parts for pk, is_nt in zip(p.pert_keys, p.nt_mask) if not is_nt]
            from collections import Counter

            counts_by_target = Counter(kept_pkeys)
            self.assertEqual(
                dict(counts_by_target),
                {"HLA-A": 8, "PPAN-P2RY11": 6, "PPAN": 4, "JUN": 8},
            )


if __name__ == "__main__":
    unittest.main(verbosity=2)
