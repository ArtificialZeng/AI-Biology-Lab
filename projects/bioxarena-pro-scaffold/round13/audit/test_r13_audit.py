#!/usr/bin/env python3
"""
Clean-Room Unit & Synthetic Cell Ranger Flex H5 Test Suite for `round13/audit/` (A-T3).

Strictly NEVER reads `round13/src/`, `round13/results/`, or any real `GSE327644` matrices.
Run strictly via:
    `python3 -I round13/audit/test_r13_audit.py -v`
"""

from __future__ import annotations

import gzip
import sys
import tempfile
import unittest
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

AUDIT_R13_DIR = Path(__file__).resolve().parent
REPO_ROOT = AUDIT_R13_DIR.parents[1]
R10_SCRIPTS_DIR = REPO_ROOT / "round10" / "shared" / "scripts"
R11_AUDIT_DIR = REPO_ROOT / "round11" / "audit"
R12_AUDIT_DIR = REPO_ROOT / "round12" / "audit"
for _p in (str(AUDIT_R13_DIR), str(R10_SCRIPTS_DIR), str(R11_AUDIT_DIR), str(R12_AUDIT_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from independent_r11_loader import DEFAULT_HGNC_TSV  # noqa: E402
from independent_r13_loader import (  # noqa: E402
    R13_RAW_DIR,
    R13_SHARED_DIR,
    build_gse327644_audit_parts,
    check_gse327644_feasibility,
    enforce_r13_gate,
    load_gse327644_guide_template,
    parse_confirmatory_lines_r13,
    standardize_gse327644_to_audit_line,
)
from prepare_r10_lines import HGNCLookup, build_hgnc_lookup, process_r10_line  # noqa: E402
from run_r13_audit import DEFAULT_R13_SEED, adjudicate_round13, render_r13_audit_markdown  # noqa: E402


def _make_mini_hgnc(genes: List[str]) -> HGNCLookup:
    return HGNCLookup(
        approved=set(genes),
        ens_to_sym={f"ENSG{i:011d}": g for i, g in enumerate(genes, start=1)},
        sym_to_ens={g: f"ENSG{i:011d}" for i, g in enumerate(genes, start=1)},
        prev_to_sym={},
        alias_to_sym={},
    )


def _write_synthetic_guide_template(
    path: Path,
    include_safe_targeting: bool = False,
) -> None:
    rows = [
        ("NTC-001", "AAA001", "Control", "n/a", "3' NGG", "n/a"),
        ("NTC-002", "AAA002", "Control", "n/a", "3' NGG", "n/a"),
        ("HLA-A-1", "TTT001", "Targeting a Gene", "ENSG00000235657", "3' NGG", "HLA"),
        ("HLA-B-1", "TTT002", "Targeting a Gene", "ENSG00000234745", "3' NGG", "HLA"),
        ("ARPC4-TTLL3-1", "TTT003", "Targeting a Gene", "ENSG00000241553", "3' NGG", "ARPC4"),
        ("PPAN-P2RY11-1", "TTT004", "Targeting a Gene", "ENSG00000130810", "3' NGG", "PPAN"),
        ("1-Jun", "TTT005", "Targeting a Gene", "ENSG00000177606", "3' NGG", "JUN"),
        ("STAT1-1", "TTT006", "Targeting a Gene", "ENSG00000115415", "3' NGG", "STAT1"),
    ]
    if include_safe_targeting:
        rows.append(("SAFE-001", "SSS001", "Control", "safe-targeting", "3' NGG", "safe-targeting"))
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        fh.write("guide_id\tguide_protospacer\tguide_role\tguide_target_gene_id\tguide_PAM\tguide_target_gene_name\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")


def _write_synthetic_cellranger_h5(
    path: Path,
    rna_counts: np.ndarray,
    crp_counts: np.ndarray,
    rna_genes: List[str],
    crp_guides: List[str],
    barcodes: List[str],
    float_rna: bool = False,
    embedded_calls: Optional[List[Tuple[int, str]]] = None,
) -> None:
    """
    Write 10x Cell Ranger CSC `/matrix` HDF5 `(n_feat, n_cells)` matching GEO `GSE327644`.
    """
    n_cells = len(barcodes)
    assert rna_counts.shape == (n_cells, len(rna_genes))
    assert crp_counts.shape == (n_cells, len(crp_guides))

    full_dense = np.vstack([rna_counts.T, crp_counts.T])  # (n_feat, n_cells)
    csc = sp.csc_matrix(full_dense if float_rna else full_dense.astype(np.int32))
    ftypes = (["Gene Expression"] * len(rna_genes)) + (["CRISPR Guide Capture"] * len(crp_guides))
    fids = [f"ENSG_G_{i:04d}" for i in range(len(rna_genes))] + list(crp_guides)
    fnames = list(rna_genes) + list(crp_guides)

    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as hf:
        mg = hf.create_group("matrix")
        mg.create_dataset("shape", data=np.asarray([len(ftypes), n_cells], dtype=np.int64))
        mg.create_dataset("data", data=csc.data.astype(np.float32 if float_rna else np.int32))
        mg.create_dataset("indices", data=csc.indices.astype(np.int32))
        mg.create_dataset("indptr", data=csc.indptr.astype(np.int32))
        mg.create_dataset("barcodes", data=np.asarray([b.encode("utf-8") for b in barcodes], dtype="S"))

        fg = mg.create_group("features")
        fg.create_dataset("feature_type", data=np.asarray([x.encode("utf-8") for x in ftypes], dtype="S"))
        fg.create_dataset("id", data=np.asarray([x.encode("utf-8") for x in fids], dtype="S"))
        fg.create_dataset("name", data=np.asarray([x.encode("utf-8") for x in fnames], dtype="S"))

        if embedded_calls is not None:
            pg = mg.create_group("protospacer_calls")
            pg.create_dataset("barcode", data=np.asarray([b.encode("utf-8") for b in barcodes], dtype="S"))
            pg.create_dataset("num_features", data=np.asarray([int(x[0]) for x in embedded_calls], dtype=np.int32))
            pg.create_dataset("feature_call", data=np.asarray([x[1].encode("utf-8") for x in embedded_calls], dtype="S"))


def _build_synthetic_rna_matrix(
    rng: np.random.Generator,
    n_cells: int,
    n_genes: int,
    umi_per_cell: int = 1200,
    high_mt_mask: Optional[np.ndarray] = None,
    effect_shift_genes: Optional[np.ndarray] = None,
    low_umi_rows: Sequence[int] = (),
    low_gene_rows: Sequence[int] = (),
) -> np.ndarray:
    mat = rng.integers(1, 5, size=(n_cells, n_genes), dtype=np.int32)
    extra = max(umi_per_cell - int(mat.sum(axis=1).mean()), 10)
    mat[:, 1:6] += extra // 5
    if high_mt_mask is not None:
        # Put large MT- gene count in index 0 (named MT-CO1) to confirm max_pct_mt = None retention
        mat[high_mt_mask, 0] += 8000
    if effect_shift_genes is not None:
        mat[effect_shift_genes, :60] += 8
    for r in low_umi_rows:
        mat[r, :] = 0
        mat[r, :205] = 1  # 205 genes, 205 UMIs (< 300 UMI -> drop)
    for r in low_gene_rows:
        mat[r, :] = 0
        mat[r, :50] = 20  # 1,000 UMIs, 50 genes (< 200 genes -> drop)
    return mat


class TestRound13CleanRoomAudit(unittest.TestCase):
    def test_1_blinding_gates_and_confirmatory_parser(self) -> None:
        # Cannot access round13/shared during synthetic_test or dev_r12
        with self.assertRaises(RuntimeError):
            enforce_r13_gate("synthetic_test", extra_dirs=[R13_RAW_DIR])
        with self.assertRaises(RuntimeError):
            enforce_r13_gate("dev_r12", target_lines=["mac_ctrl"])
        with self.assertRaises(RuntimeError):
            enforce_r13_gate("dev_r12", extra_dirs=[R13_SHARED_DIR])

        with tempfile.TemporaryDirectory() as td:
            sd = Path(td) / "status"
            sd.mkdir(parents=True)
            with self.assertRaises(RuntimeError):
                enforce_r13_gate("reader_xcheck_r13", status_dir=sd)
            (sd / "PREREG_F_COMMIT.txt").write_text("5b51e5b\n", encoding="utf-8")
            # Still blocked without R13_CODE_PATCHED
            with self.assertRaises(RuntimeError):
                enforce_r13_gate("reader_xcheck_r13", status_dir=sd)
            (sd / "R13_CODE_PATCHED").write_text("ok\n", encoding="utf-8")
            meta_x = enforce_r13_gate("reader_xcheck_r13", status_dir=sd)
            self.assertEqual(meta_x["prereg_commit"], "5b51e5b")

            # holdout_r13 still blocked without CODE_FREEZE_R13.txt & CONFIRMATORY_LINES_R13.txt
            with self.assertRaises(RuntimeError):
                enforce_r13_gate("holdout_r13", status_dir=sd)
            (sd / "CODE_FREEZE_R13.txt").write_text("deadbeef\n", encoding="utf-8")
            (sd / "CONFIRMATORY_LINES_R13.txt").write_text("mac_ctrl=CONFIRMATORY\nmac_lps=SECONDARY_ONLY\n", encoding="utf-8")
            meta_h = enforce_r13_gate("holdout_r13", status_dir=sd)
            self.assertEqual(meta_h["code_freeze"], "deadbeef")

            # mac_lps can NEVER be CONFIRMATORY
            bad_conf = sd / "BAD_CONF.txt"
            bad_conf.write_text("mac_lps=CONFIRMATORY\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                parse_confirmatory_lines_r13(bad_conf)

    def test_2_guide_template_ntc_and_target_canonicalization(self) -> None:
        base_genes = ["HLA-A", "HLA-B", "ARPC4", "ARPC4-TTLL3", "PPAN", "PPAN-P2RY11", "JUN", "STAT1"]
        hgnc = build_hgnc_lookup(DEFAULT_HGNC_TSV) if DEFAULT_HGNC_TSV.exists() else _make_mini_hgnc(base_genes)

        with tempfile.TemporaryDirectory() as td:
            t_p1 = Path(td) / "tpl_indist.txt.gz"
            _write_synthetic_guide_template(t_p1, include_safe_targeting=False)
            res1 = load_gse327644_guide_template(t_p1, hgnc)
            self.assertFalse(res1["distinguishes_safe_vs_ntc"])
            self.assertEqual(res1["ntc_guide_ids"], {"NTC-001", "NTC-002"})
            self.assertEqual(res1["guide_to_target"]["HLA-A-1"][1], "HLA-A")
            self.assertEqual(res1["guide_to_target"]["HLA-B-1"][1], "HLA-B")
            self.assertEqual(res1["guide_to_target"]["ARPC4-TTLL3-1"][1], "ARPC4-TTLL3")
            self.assertEqual(res1["guide_to_target"]["PPAN-P2RY11-1"][1], "PPAN-P2RY11")
            self.assertEqual(res1["guide_to_target"]["1-Jun"][1], "JUN")

            # Safe-targeting separation test
            t_p2 = Path(td) / "tpl_safe.txt.gz"
            _write_synthetic_guide_template(t_p2, include_safe_targeting=True)
            res2 = load_gse327644_guide_template(t_p2, hgnc)
            self.assertTrue(res2["distinguishes_safe_vs_ntc"])
            self.assertEqual(res2["ntc_guide_ids"], {"NTC-001", "NTC-002"})
            self.assertEqual(res2["excluded_safe_ids"], {"SAFE-001"})

    def test_3_fallback_poisson_gaussian_6ch_parity_no_mt_limit_and_lt15_fallback(self) -> None:
        enforce_r13_gate("synthetic_test", extra_dirs=[])
        rng = np.random.default_rng(7)
        rna_genes = ["MT-CO1"] + [f"GENE{i:04d}" for i in range(1, 240)] + ["HLA-A", "JUN", "STAT1"]
        hgnc = _make_mini_hgnc(rna_genes + ["HLA-B", "ARPC4-TTLL3", "PPAN-P2RY11"])
        crp_guides = ["NTC-001", "NTC-002", "HLA-A-1", "1-Jun", "STAT1-1"]

        with tempfile.TemporaryDirectory() as td:
            raw_dir = Path(td) / "raw" / "gse327644"
            tpl_path = raw_dir / "GSE327644_guide_template.txt.gz"
            _write_synthetic_guide_template(tpl_path, include_safe_targeting=False)

            # Reject unfiltered raw_feature_bc_matrix even if matching condition/channel
            bad_raw_h5 = raw_dir / "GSM000_Mac_L01_Mac_Ctrl_raw_feature_bc_matrix.h5"
            _write_synthetic_cellranger_h5(
                bad_raw_h5,
                np.ones((4, len(rna_genes)), dtype=np.int32),
                np.ones((4, len(crp_guides)), dtype=np.int32),
                rna_genes,
                crp_guides,
                [f"BAD_{i}" for i in range(4)],
            )

            # Create 6 channels for Mac_Ctrl (tests PREREG_F §9.1 >4 channels parity)
            # L01, L03 (20 NTC each), L05 (only 5 NTC -> <15 merges to half_1_pooled) -> Half A
            # L02, L04, L06 (20 NTC each) -> Half B
            ntc_per_ch = {1: 20, 2: 20, 3: 20, 4: 20, 5: 5, 6: 20}
            for ch in range(1, 7):
                n_nt = ntc_per_ch[ch]
                n_stat1 = 8
                n_jun = 8
                n_doublet = 3
                n_low_umi = 2
                n_low_gene = 2
                n_cells = n_nt + n_stat1 + n_jun + n_doublet + n_low_umi + n_low_gene

                mt_mask = np.zeros(n_cells, dtype=bool)
                mt_mask[0] = True  # cell 0 is NTC with >80% MT (must NOT be dropped since max_pct_mt=None)
                eff_mask = np.zeros(n_cells, dtype=bool)
                eff_mask[n_nt : n_nt + n_stat1] = True

                low_umi_idx = [n_cells - 4, n_cells - 3]
                low_gene_idx = [n_cells - 2, n_cells - 1]
                rna_mat = _build_synthetic_rna_matrix(
                    rng,
                    n_cells,
                    len(rna_genes),
                    high_mt_mask=mt_mask,
                    effect_shift_genes=eff_mask,
                    low_umi_rows=low_umi_idx,
                    low_gene_rows=low_gene_idx,
                )

                crp_mat = np.zeros((n_cells, len(crp_guides)), dtype=np.int32)
                # Background ambient noise on all guides = 1..2
                crp_mat[:, :] = rng.integers(0, 3, size=crp_mat.shape, dtype=np.int32)
                # Clear signal UMIs >= 15 for true assigned guide
                crp_mat[:n_nt, 0] = 20  # NTC-001
                crp_mat[n_nt : n_nt + n_stat1, 4] = 22  # STAT1-1
                crp_mat[n_nt + n_stat1 : n_nt + n_stat1 + n_jun, 3] = 24  # 1-Jun -> JUN
                # Multi-guide doublets (top1=20, top2=18 -> rejected!)
                d_start = n_nt + n_stat1 + n_jun
                crp_mat[d_start : d_start + n_doublet, 0] = 20
                crp_mat[d_start : d_start + n_doublet, 4] = 18
                # Low-UMI/Low-gene cells even with strong guide signal must be filtered by >=300 UMI & >=200 genes
                crp_mat[low_umi_idx[0] : low_gene_idx[-1] + 1, 0] = 25

                bcs = [f"L{ch:02d}_CELL_{i:04d}-1" for i in range(n_cells)]
                _write_synthetic_cellranger_h5(
                    raw_dir / f"GSM96625{50+ch}_Mac_L{ch:02d}_Mac_Ctrl_sample_filtered_feature_bc_matrix.h5",
                    rna_mat,
                    crp_mat,
                    rna_genes,
                    crp_guides,
                    bcs,
                )

            # Also write 2 channels of Mac_LPS with fewer cells to ensure line isolation
            for ch in (1, 2):
                n_c = 25
                rna_m = _build_synthetic_rna_matrix(rng, n_c, len(rna_genes))
                crp_m = np.zeros((n_c, len(crp_guides)), dtype=np.int32)
                crp_m[:18, 1] = 25  # NTC-002
                crp_m[18:, 2] = 25  # HLA-A-1
                _write_synthetic_cellranger_h5(
                    raw_dir / f"GSM96625{50+ch}_Mac_L{ch:02d}_Mac_LPS_sample_filtered_feature_bc_matrix.h5",
                    rna_m,
                    crp_m,
                    rna_genes,
                    crp_guides,
                    [f"LPS_L{ch:02d}_{i:04d}" for i in range(n_c)],
                )

            parts, meta = build_gse327644_audit_parts("mac_ctrl", raw_dir, hgnc, guide_template_path=tpl_path)
            self.assertEqual(meta["guide_call_mode"], "H5_CRISPR_POISSON_GAUSSIAN_LINE_FIT")
            self.assertEqual(meta["split_branch"], "STAGE1_CHANNEL_PARITY")
            self.assertEqual(meta["n_channels_found"], 6)
            self.assertEqual(meta["max_pct_mt"], None)
            self.assertEqual(meta["n_ntc_half_a"], 45)  # 20 (L01) + 20 (L03) + 5 (L05)
            self.assertEqual(meta["n_ntc_half_b"], 60)  # 20 (L02) + 20 (L04) + 20 (L06)

            # Check L05 (< 15 NTC) merged to half_1_pooled while L01/L03 retained L01/L03
            l05_part = [p for p in parts if p.tag.endswith("L05")][0]
            l01_part = [p for p in parts if p.tag.endswith("L01")][0]
            l06_part = [p for p in parts if p.tag.endswith("L06")][0]
            self.assertTrue(np.all(l05_part.batches == "half_1_pooled"))
            self.assertTrue(np.all(l05_part.half_flags == 1))
            self.assertTrue(np.all(l01_part.batches == "L01"))
            self.assertTrue(np.all(l06_part.batches == "L06"))
            self.assertTrue(np.all(l06_part.half_flags == 2))

            # Verify high-%MT NTC cell 0 in L01 WAS RETAINED (20 NTC + 16 targeting = 36 cells kept)
            self.assertEqual(l01_part.n_cells, 36)
            # Verify JUN normalization mapped 1-Jun -> JUN cleanly in pert_meta
            self.assertEqual(l01_part.pert_meta["1-Jun"][1], "JUN")
            self.assertEqual(l01_part.pert_meta["STAT1-1"][1], "STAT1")

            # Verify mac_lps isolated (2 channels only)
            lps_parts, lps_meta = build_gse327644_audit_parts("mac_lps", raw_dir, hgnc, guide_template_path=tpl_path)
            self.assertEqual(lps_meta["n_channels_found"], 2)
            self.assertEqual(sum(p.n_cells for p in lps_parts), 50)
            self.assertEqual(lps_parts[0].pert_meta["HLA-A-1"][1], "HLA-A")

    def test_4_protospacer_call_precedence_and_stage2_seed42(self) -> None:
        rng = np.random.default_rng(11)
        rna_genes = [f"GENE{i:04d}" for i in range(250)]
        hgnc = _make_mini_hgnc(rna_genes + ["STAT1"])
        crp_guides = ["NTC-001", "STAT1-1"]

        with tempfile.TemporaryDirectory() as td:
            raw_dir = Path(td) / "raw" / "gse327644"
            tpl_path = raw_dir / "guide_template.txt.gz"
            _write_synthetic_guide_template(tpl_path)

            # Single file with NO channel label -> forces STAGE2_STRATIFIED_SEED42
            # Plus companion protospacer_calls_per_cell.csv overrides CRISPR UMIs
            n_c = 50
            rna_mat = _build_synthetic_rna_matrix(rng, n_c, len(rna_genes))
            crp_mat = np.zeros((n_c, len(crp_guides)), dtype=np.int32)
            crp_mat[:40, 0] = 30
            crp_mat[40:, 1] = 30
            bcs = [f"CELL_{i:04d}-1" for i in range(n_c)]
            h5_p = raw_dir / "Mac_Ctrl_pooled_sample_filtered_feature_bc_matrix.h5"
            _write_synthetic_cellranger_h5(h5_p, rna_mat, crp_mat, rna_genes, crp_guides, bcs)

            csv_p = raw_dir / "Mac_Ctrl_pooled_protospacer_calls_per_cell.csv"
            rows_csv = []
            for i, bc in enumerate(bcs):
                if i < 40:
                    rows_csv.append((bc, "1", "NTC-001"))
                elif i < 48:
                    rows_csv.append((bc, "1", "STAT1-1"))
                else:
                    # Marked num_features=2 even though CRISPR UMI looked clean -> must be rejected!
                    rows_csv.append((bc, "2", "NTC-001|STAT1-1"))
            pd.DataFrame(rows_csv, columns=["cell_barcode", "num_features", "feature_call"]).to_csv(csv_p, index=False)

            parts, meta = build_gse327644_audit_parts("mac_ctrl", raw_dir, hgnc, guide_template_path=tpl_path)
            self.assertEqual(meta["guide_call_mode"], "CELLRANGER_PROTOSPACER_CALLS")
            self.assertEqual(meta["split_branch"], "STAGE2_STRATIFIED_SEED42")
            self.assertEqual(meta["n_kept_cells"], 48)  # 2 num_features==2 excluded!
            self.assertEqual(meta["n_ntc_half_a"], 20)
            self.assertEqual(meta["n_ntc_half_b"], 20)
            self.assertEqual(set(parts[0].batches.tolist()), {"half_1", "half_2"})

    def test_5_feasibility_rejection_process_r10_and_adjudication_rules(self) -> None:
        rng = np.random.default_rng(19)
        rna_genes = [f"GENE{i:04d}" for i in range(250)]
        hgnc = _make_mini_hgnc(rna_genes + ["STAT1"])
        crp_guides = ["NTC-001", "STAT1-1"]

        with tempfile.TemporaryDirectory() as td:
            raw_dir = Path(td) / "raw" / "gse327644"
            tpl_path = raw_dir / "guide_template.txt.gz"
            _write_synthetic_guide_template(tpl_path)

            # Float non-integer RNA counts must fail pre-metric feasibility -> EXCLUDED_FEASIBILITY
            bad_h5 = raw_dir / "GSM1_Mac_L01_Mac_Ctrl_sample_filtered_feature_bc_matrix.h5"
            rna_f = np.full((10, len(rna_genes)), 1.5, dtype=np.float32)
            crp_i = np.ones((10, len(crp_guides)), dtype=np.int32) * 10
            _write_synthetic_cellranger_h5(bad_h5, rna_f, crp_i, rna_genes, crp_guides, [f"C{i}" for i in range(10)], float_rna=True)
            feas = check_gse327644_feasibility("mac_ctrl", raw_dir, hgnc, guide_template_path=tpl_path)
            self.assertFalse(feas["feasible"])
            self.assertEqual(feas["status"], "EXCLUDED_FEASIBILITY")

            # Clean integer embedded-protospacer 2-channel fixture through full process_r10_line
            bad_h5.unlink()
            out_root = Path(td) / "lines"
            for ch in (1, 2):
                n_c = 36
                rna_m = _build_synthetic_rna_matrix(rng, n_c, len(rna_genes), effect_shift_genes=(np.arange(n_c) >= 20))
                crp_m = np.zeros((n_c, len(crp_guides)), dtype=np.int32)
                calls = [(1, "NTC-001") if i < 20 else (1, "STAT1-1") for i in range(n_c)]
                _write_synthetic_cellranger_h5(
                    raw_dir / f"GSM{ch}_Mac_L{ch:02d}_Mac_Ctrl_sample_filtered_feature_bc_matrix.h5",
                    rna_m,
                    crp_m,
                    rna_genes,
                    crp_guides,
                    [f"BC_L{ch}_{i:04d}" for i in range(n_c)],
                    embedded_calls=calls,
                )
            parts, meta = build_gse327644_audit_parts("mac_ctrl", raw_dir, hgnc, guide_template_path=tpl_path)
            self.assertEqual(meta["guide_call_mode"], "CELLRANGER_PROTOSPACER_CALLS")
            line_dir = process_r10_line("mac_ctrl", parts, hgnc=hgnc, out_root=out_root, n_workers=1)
            self.assertTrue((line_dir / "pseudobulk_delta.npz").exists())
            self.assertTrue((line_dir / "ctrl_cells.npz").exists())
            with np.load(line_dir / "pseudobulk_delta.npz", allow_pickle=False) as pb:
                self.assertEqual(list(pb["targets"]), ["STAT1"])
                self.assertEqual(int(pb["n_cells"][0]), 32)  # 16 per channel across 2 channels

            # Exercise full standardize_gse327644_to_audit_line end-to-end (with real HGNC table if cached)
            if DEFAULT_HGNC_TSV.exists():
                adata, std_meta = standardize_gse327644_to_audit_line(
                    line_name="mac_ctrl",
                    raw_dir=raw_dir,
                    cache_out_root=Path(td) / "std_lines",
                    aligned_genes=rna_genes[:200],
                    hgnc_tsv=DEFAULT_HGNC_TSV,
                    guide_template_path=tpl_path,
                    min_gene_overlap=1,
                    n_workers=1,
                )
                self.assertEqual(std_meta["guide_call_mode"], "CELLRANGER_PROTOSPACER_CALLS")
                self.assertEqual(list(adata.target_names), ["STAT1"])
                self.assertEqual(int(adata.ctrl_cells.shape[0]), 40)

        # Adjudication rules verification (seed 20261223, K_conf=0 vs K_conf=1 pass vs fail)
        self.assertEqual(DEFAULT_R13_SEED, 20261223)
        adj_k0 = adjudicate_round13({}, {"roles": {"mac_ctrl": "SECONDARY_ONLY", "mac_lps": "SECONDARY_ONLY"}, "confirmatory_lines": []})
        self.assertEqual(adj_k0["k_conf"], 0)
        self.assertEqual(adj_k0["claim1_he1_verdict"], "NO_CONFIRMATORY_ADJUDICATION_GATE_FAILED")

        adj_pass = adjudicate_round13(
            {"mac_ctrl": {"he1_pass": True, "he2_pass": True}},
            {"roles": {"mac_ctrl": "CONFIRMATORY", "mac_lps": "SECONDARY_ONLY"}, "confirmatory_lines": ["mac_ctrl"]},
        )
        self.assertEqual(adj_pass["claim1_he1_verdict"], "CONFIRMED")
        self.assertIn("two separately preregistered prospective tests", adj_pass["no_cross_round_pooling_notice"])

        adj_fail = adjudicate_round13(
            {"mac_ctrl": {"he1_pass": True, "he2_pass": False}},
            {"roles": {"mac_ctrl": "CONFIRMATORY", "mac_lps": "SECONDARY_ONLY"}, "confirmatory_lines": ["mac_ctrl"]},
        )
        self.assertEqual(adj_fail["claim1_he1_verdict"], "CONFIRMED")
        self.assertEqual(adj_fail["claim2_he2_verdict"], "FAILED_TO_REPLICATE")

        # Markdown report smoke check
        md = render_r13_audit_markdown({
            "phase": "holdout_r13",
            "source_mode": "r9_6",
            "seed": 20261223,
            "n_boot": 4000,
            "n_perm": 100,
            "gate_meta": {"prereg_commit": "5b51e5b"},
            "adjudication": adj_k0,
            "per_line": {},
        })
        self.assertIn("target-baseline-free transfer predictor", md)


if __name__ == "__main__":
    unittest.main()
