#!/usr/bin/env python3
"""Unit tests for prepare_lines.py (synthetic HDF5 fixtures only; zero held-out access)."""

from __future__ import annotations

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

import prepare_lines as pl


def _write_str_ds(grp: h5py.Group, name: str, vals):
    dt = h5py.string_dtype(encoding="utf-8")
    grp.create_dataset(name, data=np.array(list(vals), dtype=object), dtype=dt)


def _write_csr_h5ad(path: Path, mat: sp.csr_matrix, obs: dict, var_syms, gene_ids=None):
    with h5py.File(path, "w") as f:
        gx = f.create_group("X")
        gx.attrs["encoding-type"] = "csr_matrix"
        gx.attrs["shape"] = np.array(mat.shape, dtype=np.int64)
        gx.create_dataset("data", data=mat.data.astype(np.float32))
        gx.create_dataset("indices", data=mat.indices.astype(np.int32))
        gx.create_dataset("indptr", data=mat.indptr.astype(np.int64))
        go = f.create_group("obs")
        for k, v in obs.items():
            _write_str_ds(go, k, v)
        gv = f.create_group("var")
        _write_str_ds(gv, "_index", var_syms)
        if gene_ids is not None:
            _write_str_ds(gv, "gene_id", gene_ids)


class TestPrepareLines(unittest.TestCase):
    def test_csc_open_and_umi_with_zero_rows(self):
        rng = np.random.RandomState(0)
        dense = rng.poisson(2.0, size=(6, 5)).astype(np.float32)
        dense[2, :] = 0.0  # empty row edge case
        csc = sp.csc_matrix(dense)
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "csc.h5ad"
            with h5py.File(p, "w") as f:
                gl = f.create_group("layers")
                gc = gl.create_group("counts")
                gc.attrs["encoding-type"] = "csc_matrix"
                gc.attrs["shape"] = np.array(csc.shape, dtype=np.int64)
                gc.create_dataset("data", data=csc.data)
                gc.create_dataset("indices", data=csc.indices)
                gc.create_dataset("indptr", data=csc.indptr)
            with h5py.File(p, "r") as f:
                xn, shp = pl.open_x_node(f, "counts")
                self.assertEqual(shp, (6, 5))
                out = pl.read_h5ad_chunk_dense(xn, 0, 6, 5)
                np.testing.assert_allclose(out, dense)
                umis = pl.compute_row_umi_sums(xn, 6, 5, chunk_size=3)
                np.testing.assert_allclose(umis, dense.sum(axis=1))

    def test_csr_umi_zero_row_reduceat_fallback(self):
        dense = np.array([[1, 2, 0], [0, 0, 0], [4, 0, 5], [0, 3, 0]], dtype=np.float32)
        csr = sp.csr_matrix(dense)
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "csr.h5ad"
            _write_csr_h5ad(p, csr, {"_index": ["c0", "c1", "c2", "c3"]}, ["G1", "G2", "G3"])
            with h5py.File(p, "r") as f:
                xn, shp = pl.open_x_node(f, "X")
                umis = pl.compute_row_umi_sums(xn, shp[0], shp[1], chunk_size=2)
                np.testing.assert_allclose(umis, [3.0, 0.0, 9.0, 3.0])

    def test_stratified_split_halves_disjoint(self):
        idx = np.arange(50)
        batches = np.repeat(np.arange(5), 10)
        ha, hb = pl.stratified_split_halves(idx, batches, seed=42)
        self.assertEqual(len(set(ha) & set(hb)), 0)
        self.assertEqual(set(ha) | set(hb), set(idx))
        for b in range(5):
            self.assertEqual(np.sum(batches[ha] == b), 5)
            self.assertEqual(np.sum(batches[hb] == b), 5)

    def test_vcc_adapter_nt_dedup_and_gene_id(self):
        mat1 = sp.csr_matrix(np.ones((4, 3), dtype=np.float32))
        mat2 = sp.csr_matrix(np.ones((3, 3), dtype=np.float32))
        with tempfile.TemporaryDirectory() as tmp:
            p1, p2 = Path(tmp) / "t.h5ad", Path(tmp) / "v.h5ad"
            _write_csr_h5ad(
                p1,
                mat1,
                {
                    "_index": ["nt1", "nt2", "pA1", "pA2"],
                    "target_gene": ["non-targeting", "non-targeting", "TP53", "TP53"],
                    "guide_id": ["nt", "nt", "TP53_1", "TP53_1"],
                    "batch": ["B1", "B2", "B1", "B2"],
                },
                ["TP53", "MYC", "GAPDH"],
                ["ENSG1", "ENSG2", "ENSG3"],
            )
            _write_csr_h5ad(
                p2,
                mat2,
                {
                    "_index": ["nt1", "nt2", "pB1"],
                    "target_gene": ["non-targeting", "non-targeting", "MYC"],
                    "guide_id": ["nt", "nt", "MYC_1"],
                    "batch": ["B1", "B2", "B1"],
                },
                ["TP53", "MYC", "GAPDH"],
            )
            seen = set()
            with h5py.File(p1, "r") as f1, h5py.File(p2, "r") as f2:
                m1 = pl.vcc_adapter(f1, "train", sym_to_ens={"MYC": "ENS_HGNC"}, seen_nt_barcodes=seen)
                m2 = pl.vcc_adapter(f2, "val", sym_to_ens={"MYC": "ENS_HGNC"}, seen_nt_barcodes=seen)
            self.assertEqual(int(m1["nt_mask"].sum()), 2)
            self.assertEqual(m1["n_duplicate_nt_skipped"], 0)
            self.assertEqual(m1["var_gene_ids"], ["ENSG1", "ENSG2", "ENSG3"])
            self.assertEqual(int(m2["nt_mask"].sum()), 0)
            self.assertEqual(m2["n_duplicate_nt_skipped"], 2)
            self.assertEqual(m2["pert_keys"], ["", "", "val_MYC_1"])
            self.assertEqual(m2["batches"], ["B1", "B2", "B1"])

    def test_build_gene_keep_mask_force_and_dedup(self):
        vnames = ["G0", "G1", "G1", "TARGET_LOW", "G3"]
        finite = np.array([True, True, True, True, True])
        cmean = np.array([10.0, 2.0, 8.0, 0.1, 5.0])
        sel, audit = pl.build_gene_keep_mask(
            vnames, finite, cmean, extra_keep_symbols={"TARGET_LOW"}, max_genes=3
        )
        self.assertEqual(audit["n_genes_duplicate_symbol_dropped"], 1)
        self.assertEqual(audit["n_genes_low_expr_cap_dropped"], 1)
        # Expect G0 (idx 0), G1 best (idx 2), TARGET_LOW (idx 3)
        self.assertEqual(sel.tolist(), [0, 2, 3])

    def test_end_to_end_process_singlecell_line_exact_math(self):
        rng = np.random.RandomState(123)
        n_g = 6
        nt0 = rng.poisson(20, size=(16, n_g)).astype(np.float32) + 1.0
        nt1 = rng.poisson(35, size=(16, n_g)).astype(np.float32) + 2.0
        pt0 = rng.poisson(20, size=(16, n_g)).astype(np.float32) + 1.0
        pt0[:, 0] += 15.0  # strong upward effect on G0
        pt1 = rng.poisson(35, size=(16, n_g)).astype(np.float32) + 2.0
        pt1[:, 0] += 25.0
        full = np.vstack([nt0, nt1, pt0, pt1])  # 64 unique cells total: 32 NT + 32 TP53
        batches = ["B0"] * 16 + ["B1"] * 16 + ["B0"] * 16 + ["B1"] * 16
        tg = ["non-targeting"] * 32 + ["TP53"] * 32
        gd = ["nt"] * 32 + ["TP53_g1"] * 32
        bcs = [f"c{i}" for i in range(len(tg))]
        genes = [f"G{i}" for i in range(n_g)]
        genes[0] = "TP53"

        with tempfile.TemporaryDirectory() as tmp:
            tmp_p = Path(tmp)
            p_h5 = tmp_p / "syn.h5ad"
            # File 1 gets all 32 NT cells + 4 perturbed cells
            _write_csr_h5ad(
                p_h5,
                sp.csr_matrix(full[:36]),
                {"_index": bcs[:36], "target_gene": tg[:36], "guide_id": gd[:36], "batch": batches[:36]},
                genes,
            )
            p_h5_2 = tmp_p / "syn2.h5ad"
            # File 2 includes duplicate NT cell "c0" + remaining 28 perturbed cells + extra unshared gene "EXTRA"
            f2_mat = np.hstack([np.vstack([full[:1], full[36:]]), np.zeros((29, 1), dtype=np.float32)])
            _write_csr_h5ad(
                p_h5_2,
                sp.csr_matrix(f2_mat),
                {
                    "_index": [bcs[0]] + bcs[36:],
                    "target_gene": [tg[0]] + tg[36:],
                    "guide_id": [gd[0]] + gd[36:],
                    "batch": [batches[0]] + batches[36:],
                },
                genes + ["EXTRA"],
            )
            old_lines_dir = pl.LINES_DIR
            try:
                pl.LINES_DIR = tmp_p / "lines"
                seen = set()
                out_dir = pl.process_singlecell_line(
                    "syn_line",
                    [(p_h5, "f1"), (p_h5_2, "f2")],
                    lambda f, tag: pl.vcc_adapter(f, tag=tag, sym_to_ens={}, seen_nt_barcodes=seen),
                    precomputed_stats_fn=lambda: {"f1_TP53_g1": (0, 10.0, 0.5), "f2_TP53_g1": (1, 55.0, 0.0001)},
                )
            finally:
                pl.LINES_DIR = old_lines_dir

            meta = __import__("json").loads((out_dir / "line_metadata.json").read_text())
            self.assertEqual(meta["n_duplicate_nt_skipped"], 1)
            self.assertEqual(meta["n_genes_multifile_dropped"], 1)
            self.assertEqual(meta["raw_shape"], [64, 6])
            self.assertEqual(meta["n_nt_total"], 32)
            with np.load(out_dir / "ctrl_cells.npz") as zc:
                self.assertEqual(zc["X_log1p"].shape, (32, 6))
                np.testing.assert_allclose(zc["X_zscore"].mean(axis=0), 0.0, atol=1e-5)
            with np.load(out_dir / "pseudobulk_delta.npz") as zd:
                self.assertTrue(bool(zd["half_valid"][0]))
                self.assertTrue(bool(zd["reproducible"][0]))
                self.assertEqual(int(zd["n_cells"][0]), 28)
                self.assertGreater(float(zd["delta_log1p"][0, 0]), 0.0)
                # Hand-verify f2_TP53_g1 (rows 36:64, batch B0 for 36:48, B1 for 48:64)
                lib_all = full * (1e4 / full.sum(axis=1, keepdims=True))
                mu_b0 = lib_all[:16].mean(axis=0)
                mu_b1 = lib_all[16:32].mean(axis=0)
                deltas_f2 = np.vstack([lib_all[36:48] - mu_b0, lib_all[48:64] - mu_b1])
                np.testing.assert_allclose(zd["delta_libnorm"][0], deltas_f2.mean(axis=0), rtol=1e-5, atol=1e-4)


if __name__ == "__main__":
    unittest.main(verbosity=2)
