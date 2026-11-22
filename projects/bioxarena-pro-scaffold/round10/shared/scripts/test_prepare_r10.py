#!/usr/bin/env python3
"""Comprehensive unit and integration tests for R10 preprocessing, sharded AD/Energy, adapters, and verifier."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from prepare_r10_lines import (
    DLD1_TARGET_RE,
    DEFAULT_HGNC_PATH,
    HGNCLookup,
    PartSpec,
    build_hgnc_lookup,
    compute_10x_csc_rna_qc,
    extract_top2_csc,
    fit_poisson_gaussian_thresholds_from_hist,
    make_10x_csc_subset_readers,
    process_r10_line,
    select_canonical_genes,
)
from verify_r10_lines import evaluate_confirmatory_lines, verify_and_count_line


class TestHGNCLookup(unittest.TestCase):
    def test_synthetic_resolution_priority_and_ambiguity(self):
        with tempfile.TemporaryDirectory() as td:
            tsv = Path(td) / "hgnc.tsv"
            tsv.write_text(
                "\t".join(["symbol", "status", "ensembl_gene_id", "prev_symbol", "alias_symbol"])
                + "\n"
                + "\t".join(["TP53", "Approved", "ENSG00000141510", "P53", "TRP53|SHARED"])
                + "\n"
                + "\t".join(["BRCA1", "Approved", "ENSG00000012048", "IRIS", "RNF53|SHARED"])
                + "\n"
                + "\t".join(["OLDW", "Withdrawn", "ENSG00000000001", "", ""])
                + "\n",
                encoding="utf-8",
            )
            lk = build_hgnc_lookup(tsv)
            self.assertEqual(lk.resolve("TP53", ""), ("TP53", "approved"))
            self.assertEqual(lk.resolve("NOT_A_SYM", "ENSG00000141510.12"), ("TP53", "ensg"))
            self.assertEqual(lk.resolve("P53", ""), ("TP53", "prev_symbol"))
            self.assertEqual(lk.resolve("RNF53", ""), ("BRCA1", "alias_symbol"))
            self.assertEqual(lk.resolve("SHARED", ""), ("SHARED", "unmapped"))
            self.assertEqual(lk.resolve("OLDW", ""), ("OLDW", "unmapped"))

    @unittest.skipUnless(DEFAULT_HGNC_PATH.exists(), "Frozen HGNC table not present")
    def test_frozen_hgnc_r9_common_legacy_symbols(self):
        lk = build_hgnc_lookup(DEFAULT_HGNC_PATH)
        self.assertEqual(lk.resolve("TP53", "")[0], "TP53")
        self.assertEqual(lk.resolve("H3F3A", ""), ("H3-3A", "prev_symbol"))
        self.assertEqual(lk.resolve("", "ENSG00000141510"), ("TP53", "ensg"))


class TestSelectCanonicalGenes(unittest.TestCase):
    def test_collision_and_force_keep_cap(self):
        lk = HGNCLookup(
            approved={"G1", "G2", "G3", "G4"},
            ens_to_sym={"ENSG1": "G1", "ENSG2": "G2"},
            sym_to_ens={"G1": "ENSG1", "G2": "ENSG2", "G3": "G3_E", "G4": "G4_E"},
            prev_to_sym={"OLD_G1": "G1"},
            alias_to_sym={},
        )
        var_ids = ["ENSG1", "ENSG_DUP", "ENSG2", "ENSG3", "ENSG4", ""]
        var_names = ["OLD_G1", "G1", "G2", "G3", "G4", ""]
        finite = np.array([True, True, True, False, True, True])
        cmean = np.array([50.0, 10.0, 20.0, 999.0, 5.0, 100.0])

        sel_idx, sel_names, sel_eids, audit = select_canonical_genes(
            var_ids,
            var_names,
            finite,
            cmean,
            lk,
            max_genes=2,
            extra_keep_symbols={"G4"},
        )
        self.assertEqual(sel_idx.tolist(), [0, 4])
        self.assertEqual(sel_names, ["G1", "G4"])
        self.assertEqual(audit["n_genes_duplicate_symbol_dropped"], 1)
        self.assertEqual(audit["n_genes_empty_symbol_dropped"], 1)
        self.assertEqual(audit["n_genes_low_expr_cap_dropped"], 1)


def _build_synthetic_fixture(td_p: Path):
    rng = np.random.RandomState(7)
    G_raw = 65
    var_names = [f"GENE{i:02d}" for i in range(G_raw)]
    var_names[0] = "TP53"
    var_names[1] = "BRCA1"
    var_ids = [f"ENSG{i:011d}" for i in range(G_raw)]

    lk = HGNCLookup(
        approved=set(var_names),
        ens_to_sym={eid: g for eid, g in zip(var_ids, var_names)},
        sym_to_ens={g: eid for eid, g in zip(var_ids, var_names)},
        prev_to_sym={"P53_OLD": "TP53"},
        alias_to_sym={},
    )

    def make_cells(n_cells: int, base_row: np.ndarray) -> np.ndarray:
        noise = rng.normal(0.0, 0.25, size=(n_cells, G_raw))
        X = np.maximum(base_row[None, :] + noise, 0.05)
        X = X * (10000.0 / X.sum(axis=1, keepdims=True))
        return X

    base_a = np.full(G_raw, 100.0)
    base_a[0] = 400.0
    base_a[1] = 300.0
    base_b = base_a.copy()
    base_b[:10] += 80.0

    Xa_nt0 = make_cells(20, base_a)
    Xa_nt1 = make_cells(5, base_a)
    Xa_tp53 = make_cells(15, base_a)
    Xa_tp53[:, 0] *= 0.25
    Xa_tp53[:, 2:60] += 120.0
    Xa_tp53 *= 10000.0 / Xa_tp53.sum(axis=1, keepdims=True)
    Xa_brca = make_cells(24, base_a)
    Xa_brca[:, 2:60] += 90.0
    Xa_brca *= 10000.0 / Xa_brca.sum(axis=1, keepdims=True)
    Xa_u10 = make_cells(5, base_a)

    X_part_a = np.vstack([Xa_nt0, Xa_nt1, Xa_tp53, Xa_brca, Xa_u10])
    nt_a = np.array([True] * 25 + [False] * (15 + 24 + 5))
    batches_a = np.array(["A0"] * 20 + ["A1"] * 5 + ["A0"] * 15 + ["A0"] * 24 + ["A0"] * 5)
    hf_a = np.full(len(X_part_a), 1, dtype=np.int8)
    pk_a = [""] * 25 + ["TP53_g1"] * 15 + ["BRCA1_g1"] * 24 + ["U10_g1"] * 5

    Xb_nt = make_cells(24, base_b)
    Xb_tp53 = make_cells(12, base_b)
    Xb_tp53[:, 0] *= 0.25
    Xb_tp53[:, 2:60] += 120.0
    Xb_tp53 *= 10000.0 / Xb_tp53.sum(axis=1, keepdims=True)
    Xb_brca = make_cells(1, base_b)
    Xb_b_only = make_cells(12, base_b)

    X_part_b = np.vstack([Xb_nt, Xb_tp53, Xb_brca, Xb_b_only])
    nt_b = np.array([True] * 24 + [False] * (12 + 1 + 12))
    batches_b = np.array(["B0"] * len(X_part_b))
    hf_b = np.full(len(X_part_b), 2, dtype=np.int8)
    pk_b = [""] * 24 + ["TP53_g1"] * 12 + ["BRCA1_g1"] * 1 + ["BONLY_g1"] * 12

    pmeta = {
        "TP53_g1": ("ENSG00000000000", "TP53", "P53_OLD", "prev_symbol"),
        "BRCA1_g1": ("ENSG00000000001", "BRCA1", "BRCA1", "approved"),
        "U10_g1": ("", "GENE10", "GENE10", "approved"),
        "BONLY_g1": ("", "GENE11", "GENE11", "approved"),
    }

    pa_path = td_p / "part_a.h5"
    pb_path = td_p / "part_b.h5"
    pa_path.write_bytes(b"a")
    pb_path.write_bytes(b"b")

    part_a = PartSpec(
        tag="part_a", path=pa_path, n_cells=len(X_part_a), n_genes_raw=G_raw,
        nt_mask=nt_a, batches=batches_a, half_flags=hf_a, umi=X_part_a.sum(axis=1),
        pert_keys=pk_a, pert_meta=pmeta, var_gene_ids=var_ids, var_gene_names=var_names,
        read_chunk_fn=lambda s, e, arr=X_part_a: arr[s:e],
        read_rows_fn=lambda rows, arr=X_part_a: arr[rows],
    )
    part_b = PartSpec(
        tag="part_b", path=pb_path, n_cells=len(X_part_b), n_genes_raw=G_raw,
        nt_mask=nt_b, batches=batches_b, half_flags=hf_b, umi=X_part_b.sum(axis=1),
        pert_keys=pk_b, pert_meta=pmeta, var_gene_ids=var_ids, var_gene_names=var_names,
        read_chunk_fn=lambda s, e, arr=X_part_b: arr[s:e],
        read_rows_fn=lambda rows, arr=X_part_b: arr[rows],
    )
    return lk, [part_a, part_b], (Xa_nt0, Xa_nt1, Xa_tp53, Xb_nt, Xb_tp53, G_raw)


class TestProcessR10LineSynthetic(unittest.TestCase):
    def test_disjoint_split_halves_s7_and_sharded_mmap_exact_parity(self):
        with tempfile.TemporaryDirectory() as td:
            td_p = Path(td)
            lk, parts, arrays = _build_synthetic_fixture(td_p)
            Xa_nt0, Xa_nt1, Xa_tp53, Xb_nt, Xb_tp53, G_raw = arrays

            # 1. Run in-RAM single-worker baseline
            out_ram = process_r10_line(
                "synth_ram", parts, lk, max_ctrl_cells=10, chunk_size=17,
                seed=42, out_root=td_p / "out1", n_workers=1,
            )
            # 2. Force sharded SSD binary stream path (`max_ram_z_bytes=1`) + spawn workers and compare bit-for-bit!
            out_shd = process_r10_line(
                "synth_shd", parts, lk, max_ctrl_cells=10, chunk_size=17,
                seed=42, out_root=td_p / "out2", max_ram_z_bytes=1, shard_max_bytes=5000, n_workers=2,
            )

            # Assert zero leftover _tmp_shard_* temp files!
            self.assertEqual(list(out_shd.glob("_tmp_shard_*")), [])

            ctrl = np.load(out_ram / "ctrl_cells.npz")
            pb = np.load(out_ram / "pseudobulk_delta.npz")
            ctrl_s = np.load(out_shd / "ctrl_cells.npz")
            pb_s = np.load(out_shd / "pseudobulk_delta.npz")
            meta = json.loads((out_ram / "line_metadata.json").read_text(encoding="utf-8"))

            for k in ("delta_log1p", "delta_libnorm", "ad_counts", "energy_test_p_value", "reproducible", "half_valid"):
                self.assertTrue(np.array_equal(pb[k], pb_s[k]), f"Shard parity mismatch on {k}")
            self.assertTrue(np.array_equal(ctrl["X_log1p"], ctrl_s["X_log1p"]))

            # Disjoint halves & S7 math
            self.assertEqual(ctrl["X_libnorm"].shape, (10, G_raw))
            self.assertEqual(int(ctrl["n_ctrl_total"]), 49)
            self.assertAlmostEqual(float(ctrl["half_a_ctrl_mean_libnorm"][0]), float(np.vstack([Xa_nt0, Xa_nt1])[:, 0].mean()), places=2)
            self.assertAlmostEqual(float(ctrl["half_b_ctrl_mean_libnorm"][0]), float(Xb_nt[:, 0].mean()), places=2)
            self.assertFalse(np.allclose(ctrl["ctrl_mean_libnorm_full"], ctrl["ctrl_mean_libnorm_sub20k"]))

            self.assertEqual(pb["targets"].tolist(), ["BRCA1", "GENE11", "TP53"])
            self.assertFalse(bool(pb["half_valid"][0]))
            self.assertFalse(bool(pb["half_valid"][1]))
            self.assertTrue(bool(pb["half_valid"][2]))
            self.assertTrue(bool(pb["reproducible"][0]))
            self.assertFalse(bool(pb["reproducible"][1]))
            self.assertTrue(bool(pb["reproducible"][2]))
            self.assertAlmostEqual(float(pb["half_a_delta_libnorm"][2, 0]), float(Xa_tp53[:, 0].mean() - Xa_nt0[:, 0].mean()), places=3)
            self.assertAlmostEqual(float(pb["half_b_delta_libnorm"][2, 0]), float(Xb_tp53[:, 0].mean() - Xb_nt[:, 0].mean()), places=3)
            self.assertAlmostEqual(float(pb["knockdown_depth_log1p"][2]), float(pb["delta_log1p"][2, 0]), places=5)
            self.assertAlmostEqual(float(pb["knockdown_fraction"][2]), float(-pb["delta_libnorm"][2, 0] / ctrl["ctrl_mean_libnorm_full"][0]), places=5)
            self.assertEqual(meta["n_targets_half_valid"], 1)
            self.assertEqual(meta["n_targets_reproducible"], 2)

            # Verifier contract check on output directory
            r9_ref = {
                "r9_common_raw_set": set(pb["targets"].tolist()) | {"GENE02"},
                "r9_common_canon_set": {"TP53", "BRCA1", "GENE11", "GENE02"},
                "r9_rep_raw_union": {"TP53"},
                "r9_rep_canon_union": {"TP53"},
                "per_r9_canon": {},
            }
            v_res = verify_and_count_line(out_shd, lk, r9_ref)
            self.assertEqual(v_res["n_shared_reproducible_eval_canon"], 1)
            self.assertEqual(v_res["n_targets_reproducible_half_valid"], 1)


class Test10xCSCHelpersAndAdapters(unittest.TestCase):
    def test_10x_csc_qc_subset_reader_and_guide_em(self):
        # Synthetic 10x CSC matrix containing 4 GEX rows (1 MT-) + 2 CRISPR rows across 6 cells
        X = np.array([
            [10, 50, 20, 40, 30, 60],  # G1 (GEX)
            [0,   2,  1,  0,  3,  4],  # G2 (GEX)
            [5,   1, 90,  2,  1,  1],  # MT-CO1 (GEX) -> cell 2 is high MT!
            [20, 30, 10, 80, 15, 25],  # G4 (GEX)
            [15,  1, 12,  0, 20, 18],  # GuideA (CRISPR)
            [0,  22,  0, 19,  1, 25],  # GuideB (CRISPR)
        ], dtype=np.float32)
        C = sp.csc_matrix(X)
        with tempfile.TemporaryDirectory() as td:
            hp = Path(td) / "test10x.h5"
            with h5py.File(hp, "w") as f:
                m = f.create_group("matrix")
                m.create_dataset("data", data=C.data)
                m.create_dataset("indices", data=C.indices)
                m.create_dataset("indptr", data=C.indptr)
                m.create_dataset("shape", data=np.array(C.shape, dtype=np.int32))
                ft = m.create_group("features")
                ft.create_dataset("name", data=np.array([b"G1", b"G2", b"MT-CO1", b"G4", b"GuideA", b"GuideB"]))
                ft.create_dataset("id", data=np.array([b"E1", b"E2", b"E3", b"E4", b"GuideA", b"GuideB"]))
                ft.create_dataset("feature_type", data=np.array([
                    b"Gene Expression", b"Gene Expression", b"Gene Expression", b"Gene Expression",
                    b"CRISPR Guide Capture", b"CRISPR Guide Capture",
                ]))

            rna_rows = np.array([0, 1, 2, 3], dtype=np.int64)
            mt_raw = np.array([False, False, True, False, False, False], dtype=bool)
            umi, ng, pm = compute_10x_csc_rna_qc(hp, rna_rows, mt_raw, batch_cols=2)
            # Cell 2 has MT-CO1=90 out of 121 -> ~74.4% MT! Cell 1 has GEX sum=83, ngenes=4
            self.assertAlmostEqual(float(umi[1]), 83.0)
            self.assertEqual(int(ng[1]), 4)
            self.assertGreater(float(pm[2]), 70.0)
            self.assertLess(float(pm[0]), 20.0)

            keep_cols = np.array([0, 3, 5], dtype=np.int64)
            read_chunk, read_rows = make_10x_csc_subset_readers(hp, keep_cols, 6, rna_rows)
            dense_all = read_chunk(0, 3)
            self.assertEqual(dense_all.shape, (3, 4))
            np.testing.assert_allclose(dense_all, X[:4, [0, 3, 5]].T)
            dense_sub = read_rows(np.array([0, 2], dtype=np.int64))  # sparse mode (<15% density across span)
            np.testing.assert_allclose(dense_sub, X[:4, [0, 5]].T)

            # Test extract_top2_csc & bimodal Poisson-Gaussian EM threshold finder
            t1, t2, a1 = extract_top2_csc(C[4:, :])
            self.assertEqual(t1[5], 25.0)
            self.assertEqual(t2[5], 18.0)
            self.assertEqual(a1[5], 1)

            hist = np.zeros((3, 66), dtype=np.int64)
            # Guide 0: ambient mode at UMI 1-2 (400 cells) + signal mode at UMI 25-35 (200 cells)
            hist[0, 1] = 280; hist[0, 2] = 120; hist[0, 25] = 100; hist[0, 35] = 100
            # Guide 1: too few entries -> default 7
            hist[1, 1] = 10
            thrs = fit_poisson_gaussian_thresholds_from_hist(hist, min_thr=5, max_thr=20, default_thr=7)
            self.assertEqual(int(thrs[1]), 7)
            self.assertGreaterEqual(int(thrs[0]), 5)
            self.assertLessEqual(int(thrs[0]), 20)

    def test_dld1_regex_and_confirmatory_gate(self):
        self.assertEqual(DLD1_TARGET_RE.match("BRD4_ENST00000263377(W)_1").groups(), ("BRD4", None))
        self.assertEqual(DLD1_TARGET_RE.match("DIABLO_ENSG00000184047_P1_2").groups(), ("DIABLO", "ENSG00000184047"))
        counts = {
            "dld1": {"g_eval_canon": 5500, "n_shared_reproducible_eval_canon": 60},
            "cd4t_rest": {"g_eval_canon": 5200, "n_shared_reproducible_eval_canon": 50},
            "wtc11_i3neuron": {"g_eval_canon": 5900, "n_shared_reproducible_eval_canon": 80},
            "calu3": {"g_eval_canon": 5800, "n_shared_reproducible_eval_canon": 55},
        }
        passing, txt = evaluate_confirmatory_lines(counts)
        self.assertEqual(passing, ["dld1", "cd4t_rest"])
        self.assertIn("CONFIRMATORY_LINES=dld1,cd4t_rest", txt)
        self.assertIn("K_CONF=2", txt)


if __name__ == "__main__":
    unittest.main(verbosity=2)
