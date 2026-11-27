#!/usr/bin/env python3
"""
Unit tests for model-free Round 16 Addendum B (`Occam`) adapter components.

Uses strictly:
  - Synthetic data fixtures written inside `tempfile.TemporaryDirectory()`;
  - Small public text sidecars already read (`rhj_gene_axis.csv`, `rhj_perturbations.csv`,
    `occam_master_gene_set.csv`, `cohort.json`, `splits_manifest.json`) plus read-only locked
    Round 16 barcode/target/gene artifacts.

Zero model weights, zero Occam code imports, zero GPU/MPS, zero downloads, zero writes outside temp dirs.
"""

from __future__ import annotations

import ast
import csv
import gzip
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import List, Optional, Sequence

import h5py
import numpy as np
import scipy.sparse as sp

ADDB_SRC_DIR = Path(__file__).resolve().parent
ADDB_ROOT = ADDB_SRC_DIR.parent
REPO_ROOT = ADDB_ROOT.parent
R10_SRC_DIR = REPO_ROOT / "round10" / "src"

for _p in (str(ADDB_SRC_DIR), str(R10_SRC_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from export_addb import (  # noqa: E402
    DELTA_CROSS_PARTY_ATOL,
    Z_CTRL_CROSS_PARTY_ATOL,
    build_occam_delta_payload,
    compare_exported_npz_arrays,
    save_occam_line_deltas,
    save_occam_z_ctrl_latents,
    verify_exported_manifest,
)
from gate5b_panel_addb import (  # noqa: E402
    GATE5B_EXPECTED_COHORT_SHA256,
    GATE5B_EXPECTED_TOO_SHALLOW,
    Gate5BPanelCountError,
    build_gate5b_hepg2_panel,
    compare_splits_against_manifest,
    draw_splithalf_splits,
)
from preprocess_addb import (  # noqa: E402
    OCCAM_TARGET_SUM,
    RHJ_GENE_AXIS_SHA256,
    RHJ_N_AXIS_GENES,
    RhjGeneAxis,
    build_output_gene_projection,
    build_pairwise_common_axis_mask,
    dedup_leading_ids,
    load_master_pert2id_csv,
    load_rhj_gene_axis,
    load_rhj_perturbations_table,
    normalize_counts_target_sum_log1p,
    resolve_perturbation_keys_two_tier,
    scatter_full_44728_api,
    scatter_rhj_6983_to_encoder,
)
from raw_controls_addb import (  # noqa: E402
    DEFAULT_PRIMARY_BARCODES_JSON,
    DEFAULT_SENSITIVITY_BARCODES_JSON,
    load_10x_multi_h5_controls_by_barcodes,
    load_a549_mtx_controls_by_barcodes,
    load_feng_csv_controls_by_barcodes,
    load_locked_r16_barcodes,
    load_nadig_controls_by_barcodes,
    read_h5ad_table_col,
)


def _find_public_sidecar(candidates: Sequence[str]) -> Optional[Path]:
    """Locate an already-read small public text file without any network access."""
    roots: List[Path] = []
    if os.environ.get("ADDB_SIDECAR_DIR"):
        roots.append(Path(os.environ["ADDB_SIDECAR_DIR"]))
    roots.append(ADDB_ROOT / "assets" / "sidecars")
    roots.append(
        Path(
            "/Users/mac/.gemini/antigravity/brain/"
            "19f7514f-ce6a-4297-92cb-9803e648afd1/scratch/addb_inspect"
        )
    )
    for root in roots:
        if not root.is_dir():
            continue
        for name in candidates:
            p = root / name
            if p.is_file():
                return p
    return None


def _write_synthetic_h5ad(
    path: Path,
    X: np.ndarray,
    barcodes: Sequence[str],
    obs_gene: Sequence[str],
    obs_batch: Sequence[str],
    var_gene_id: Sequence[str],
    var_gene_name: Sequence[str],
    *,
    cat_mode: str = "legacy",  # "legacy", "modern", "plain", or "corrupt_raw_int"
    x_format: str = "dense",  # "dense", "csr", "csc"
) -> None:
    str_dt = h5py.string_dtype(encoding="utf-8")
    with h5py.File(path, "w") as f:
        if x_format == "dense":
            f.create_dataset("X", data=np.asarray(X, dtype=np.float32))
        elif x_format == "csr":
            mat = sp.csr_matrix(np.asarray(X, dtype=np.float32))
            gx = f.create_group("X")
            gx.attrs["encoding-type"] = "csr_matrix"
            gx.attrs["shape"] = np.asarray(mat.shape, dtype=np.int64)
            gx.create_dataset("data", data=mat.data.astype(np.float32))
            gx.create_dataset("indices", data=mat.indices.astype(np.int32))
            gx.create_dataset("indptr", data=mat.indptr.astype(np.int64))
        elif x_format == "csc":
            mat = sp.csc_matrix(np.asarray(X, dtype=np.float32))
            gx = f.create_group("X")
            gx.attrs["encoding-type"] = "csc_matrix"
            gx.attrs["shape"] = np.asarray(mat.shape, dtype=np.int64)
            gx.create_dataset("data", data=mat.data.astype(np.float32))
            gx.create_dataset("indices", data=mat.indices.astype(np.int32))
            gx.create_dataset("indptr", data=mat.indptr.astype(np.int64))
        else:
            raise ValueError(x_format)

        obs = f.create_group("obs")
        obs.attrs["_index"] = "cell_barcode"
        obs.create_dataset("cell_barcode", data=np.asarray(barcodes, dtype=object), dtype=str_dt)
        obs.create_dataset("gem_group", data=np.asarray(obs_batch, dtype=object), dtype=str_dt)

        var = f.create_group("var")
        var.attrs["_index"] = "gene_id"
        var.create_dataset("gene_id", data=np.asarray(var_gene_id, dtype=object), dtype=str_dt)

        def _write_col(grp: h5py.Group, name: str, vals: Sequence[str]) -> None:
            cats = sorted(set(vals))
            c2i = {c: i for i, c in enumerate(cats)}
            codes = np.asarray([c2i[v] for v in vals], dtype=np.int32)
            if cat_mode == "legacy":
                grp.create_dataset(name, data=codes)
                cgrp = grp.require_group("__categories")
                cgrp.create_dataset(name, data=np.asarray(cats, dtype=object), dtype=str_dt)
            elif cat_mode == "modern":
                cg = grp.create_group(name)
                cg.create_dataset("categories", data=np.asarray(cats, dtype=object), dtype=str_dt)
                cg.create_dataset("codes", data=codes)
            elif cat_mode == "plain":
                grp.create_dataset(name, data=np.asarray(vals, dtype=object), dtype=str_dt)
            elif cat_mode == "corrupt_raw_int":
                grp.create_dataset(name, data=codes)
            else:
                raise ValueError(cat_mode)

        _write_col(obs, "gene", obs_gene)
        _write_col(var, "gene_name", var_gene_name)


def _write_synthetic_10x_h5(
    path: Path,
    X_cells_by_feats: np.ndarray,
    barcodes: Sequence[str],
    feat_ids: Sequence[str],
    feat_names: Sequence[str],
    feat_types: Sequence[str],
) -> None:
    str_dt = h5py.string_dtype(encoding="utf-8")
    # 10x stores features x cells in CSC
    csc = sp.csc_matrix(np.asarray(X_cells_by_feats, dtype=np.float32).T)
    with h5py.File(path, "w") as f:
        m = f.create_group("matrix")
        m.create_dataset("shape", data=np.asarray(csc.shape, dtype=np.int64))
        m.create_dataset("data", data=csc.data.astype(np.float32))
        m.create_dataset("indices", data=csc.indices.astype(np.int32))
        m.create_dataset("indptr", data=csc.indptr.astype(np.int64))
        m.create_dataset("barcodes", data=np.asarray(barcodes, dtype=object), dtype=str_dt)
        fg = m.create_group("features")
        fg.create_dataset("id", data=np.asarray(feat_ids, dtype=object), dtype=str_dt)
        fg.create_dataset("name", data=np.asarray(feat_names, dtype=object), dtype=str_dt)
        fg.create_dataset("feature_type", data=np.asarray(feat_types, dtype=object), dtype=str_dt)


class TestAddendumBModelFreeAdapter(unittest.TestCase):
    """Comprehensive model-free unit tests for Round 16 Addendum B (`Occam`)."""

    def test_01_d3_l1_legacy_categorical_and_nadig_x_formats(self) -> None:
        rng = np.random.default_rng(101)
        X = rng.poisson(3.0, size=(20, 6)).astype(np.float32)
        bcs = [f"CELL_{i:03d}-1" for i in range(20)]
        genes = ["non-targeting"] * 12 + ["TP53"] * 8
        batches = ["1", "2"] * 10
        v_ids = [f"ENSG0000000000{j}.2" for j in range(6)]
        v_syms = [f"SYM{j}" for j in range(6)]
        req_bcs = ["CELL_009-1", "CELL_001-1", "CELL_011-1", "CELL_000-1"]
        expected_rows = [9, 1, 11, 0]

        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            ref_counts = None
            for cat_mode in ("legacy", "modern", "plain"):
                for x_fmt in ("dense", "csr", "csc"):
                    fp = tdp / f"nadig_{cat_mode}_{x_fmt}.h5ad"
                    _write_synthetic_h5ad(
                        fp, X, bcs, genes, batches, v_ids, v_syms, cat_mode=cat_mode, x_format=x_fmt
                    )
                    pool = load_nadig_controls_by_barcodes(
                        fp, req_bcs, line_name="nadig_hepg2", verify_ntc=True
                    )
                    self.assertEqual(pool.barcodes, req_bcs)
                    self.assertEqual(pool.raw_ensembl_ids, [f"ENSG0000000000{j}" for j in range(6)])
                    self.assertEqual(pool.raw_gene_symbols, v_syms)
                    np.testing.assert_array_equal(pool.raw_counts, X[expected_rows])
                    if ref_counts is None:
                        ref_counts = pool.raw_counts
                    else:
                        np.testing.assert_array_equal(pool.raw_counts, ref_counts)

            # Non-NTC barcode must fail verify_ntc=True
            fp_leg = tdp / "nadig_legacy_dense.h5ad"
            with self.assertRaises(ValueError):
                load_nadig_controls_by_barcodes(
                    fp_leg, ["CELL_000-1", "CELL_015-1"], line_name="nadig_hepg2", verify_ntc=True
                )

            # Missing barcode must fail
            with self.assertRaises(KeyError):
                load_nadig_controls_by_barcodes(
                    fp_leg, ["CELL_000-1", "DOES_NOT_EXIST"], line_name="nadig_hepg2", verify_ntc=True
                )

            # Lock L1: corrupt undecoded integer dataset must raise ValueError
            fp_bad = tdp / "nadig_corrupt.h5ad"
            _write_synthetic_h5ad(
                fp_bad, X, bcs, genes, batches, v_ids, v_syms, cat_mode="corrupt_raw_int", x_format="dense"
            )
            with h5py.File(fp_bad, "r") as f:
                with self.assertRaises(ValueError):
                    read_h5ad_table_col(f["obs"], "gene", forbid_raw_int_codes=True)
            with self.assertRaises(ValueError):
                load_nadig_controls_by_barcodes(
                    fp_bad, req_bcs, line_name="nadig_hepg2", verify_ntc=True
                )

    def test_02_raw_control_loaders_10x_d1_a549_feng(self) -> None:
        rng = np.random.default_rng(202)
        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)

            # (A) Multi-lane 10x .h5 with D1 gene permutation + CRISPR Guide Capture exclusion
            gex_ids_1 = ["ENSG00001.1", "ENSG00002.3", "ENSG00003.1", "ENSG00004.2"]
            gex_syms_1 = ["G1", "G2", "G3", "G4"]
            x1_gex = rng.poisson(4.0, size=(5, 4)).astype(np.float32)
            x1_full = np.hstack([x1_gex, np.full((5, 1), 99.0, dtype=np.float32)])
            p1 = tdp / "lane1.h5"
            _write_synthetic_10x_h5(
                p1,
                x1_full,
                [f"BC1_{i}" for i in range(5)],
                gex_ids_1 + ["GUIDE_1"],
                gex_syms_1 + ["GUIDE_1"],
                ["Gene Expression"] * 4 + ["CRISPR Guide Capture"],
            )

            # Lane 2 has reversed Gene Expression column order (D1)
            perm = [2, 0, 3, 1]
            x2_gex_canonical = rng.poisson(5.0, size=(4, 4)).astype(np.float32)
            x2_full_perm = np.hstack(
                [np.full((4, 1), 77.0, dtype=np.float32), x2_gex_canonical[:, perm]]
            )
            p2 = tdp / "lane2.h5"
            _write_synthetic_10x_h5(
                p2,
                x2_full_perm,
                [f"BC2_{i}" for i in range(4)],
                ["GUIDE_1"] + [gex_ids_1[k] for k in perm],
                ["GUIDE_1"] + [gex_syms_1[k] for k in perm],
                ["CRISPR Guide Capture"] + ["Gene Expression"] * 4,
            )

            req_10x = ["L2:BC2_3", "L1:BC1_1", "L2:BC2_0", "L1:BC1_4"]
            pool_10x = load_10x_multi_h5_controls_by_barcodes(
                {"L1": p1, "L2": p2},
                req_10x,
                line_name="cd4t_rest",
                canonical_ref_prefix="L1",
            )
            self.assertEqual(pool_10x.barcodes, req_10x)
            self.assertEqual(pool_10x.raw_ensembl_ids, ["ENSG00001", "ENSG00002", "ENSG00003", "ENSG00004"])
            self.assertEqual(pool_10x.audit["reordered_lane_prefixes"], ["L2"])
            expected_10x = np.vstack([x2_gex_canonical[3], x1_gex[1], x2_gex_canonical[0], x1_gex[4]])
            np.testing.assert_array_equal(pool_10x.raw_counts, expected_10x)

            # Lane 3 with different gene set MUST raise ValueError
            p3 = tdp / "lane3_bad.h5"
            _write_synthetic_10x_h5(
                p3,
                x1_gex,
                [f"BC3_{i}" for i in range(5)],
                ["ENSG00001", "ENSG00002", "ENSG00003", "ENSG99999"],
                ["G1", "G2", "G3", "BAD"],
                ["Gene Expression"] * 4,
            )
            with self.assertRaises(ValueError):
                load_10x_multi_h5_controls_by_barcodes(
                    {"L1": p1, "L3": p3}, ["L3:BC3_0"], line_name="cd4t_rest", canonical_ref_prefix="L1"
                )

            # (B) A549 MTX loader (CAS9 excluded, both orientations)
            for cell_first in (True, False):
                adir = tdp / f"a549_{cell_first}"
                adir.mkdir()
                n_c, n_g_tot = 6, 5  # col 2 is Cas9
                mat_full = rng.poisson(2.0, size=(n_c, n_g_tot)).astype(np.float32)
                mat_full[:, 2] = 50.0  # Cas9 transgene counts must be dropped
                with gzip.open(adir / "GSE345058_cell_annotations.csv.gz", "wt") as fh:
                    w = csv.writer(fh)
                    w.writerow(["cell_id", "Batch", "sgRNA", "Perturbation"])
                    for i in range(n_c):
                        w.writerow([f"A549_C{i}", "batch1", "sgNTC", "non-targeting"])
                with gzip.open(adir / "GSE345058_gene_annotations.csv.gz", "wt") as fh:
                    w = csv.writer(fh)
                    w.writerow(["index", "ensembl_id", "gene_symbol", "genome"])
                    w.writerow(["G0", "ENSG0010.1", "G0", "GRCh38"])
                    w.writerow(["G1", "ENSG0020.2", "G1", "GRCh38"])
                    w.writerow(["Cas9", "Cas9", "Cas9", "CAS9"])
                    w.writerow(["G2", "ENSG0030.3", "G2", "GRCh38"])
                    w.writerow(["G3", "ENSG0040.4", "G3", "GRCh38"])
                mtx_mat = mat_full if cell_first else mat_full.T
                coo = sp.coo_matrix(mtx_mat)
                with gzip.open(adir / "GSE345058_raw_counts.mtx.gz", "wt") as fh:
                    fh.write("%%MatrixMarket matrix coordinate real general\n")
                    fh.write(f"{mtx_mat.shape[0]} {mtx_mat.shape[1]} {coo.nnz}\n")
                    for r, c, v in zip(coo.row, coo.col, coo.data):
                        fh.write(f"{r + 1} {c + 1} {int(v)}\n")
                req_a = ["A549_C4", "A549_C1", "A549_C0"]
                pool_a = load_a549_mtx_controls_by_barcodes(adir, req_a, chunk_lines=3)
                self.assertEqual(pool_a.raw_ensembl_ids, ["ENSG0010", "ENSG0020", "ENSG0030", "ENSG0040"])
                self.assertEqual(pool_a.audit["n_excluded_features"], 1)
                np.testing.assert_array_equal(pool_a.raw_counts, mat_full[[4, 1, 0]][:, [0, 1, 3, 4]])

            # (C) Feng CSV loader (CRISPR-Guide-Capture excluded, col_offset 0 and 1)
            for col_offset in (0, 1):
                fdir = tdp / f"feng_off{col_offset}"
                fdir.mkdir()
                c_bcs = ["P1_I1_A", "P2_I2_B", "P3_I3_C"]
                f_mat = np.array(
                    [
                        [1, 2, 3],  # ENSG01
                        [9, 9, 9],  # guide row -> must be skipped
                        [4, 5, 6],  # ENSG02
                    ],
                    dtype=np.float32,
                )
                with gzip.open(fdir / "TargetedScreen_RNA-UMI-Counts.csv.gz", "wt", newline="") as fh:
                    w = csv.writer(fh)
                    w.writerow(([""] + c_bcs) if col_offset == 0 else c_bcs)
                    w.writerow(["ENSG0001.1:SYM1:Gene-Expression"] + [int(v) for v in f_mat[0]])
                    w.writerow(["GUIDE1:GUIDE1:CRISPR-Guide-Capture"] + [int(v) for v in f_mat[1]])
                    w.writerow(["ENSG0002.9:SYM2:Gene-Expression"] + [int(v) for v in f_mat[2]])
                pool_f = load_feng_csv_controls_by_barcodes(fdir, ["P3_I3_C", "P1_I1_A"])
                self.assertEqual(pool_f.raw_ensembl_ids, ["ENSG0001", "ENSG0002"])
                self.assertEqual(pool_f.raw_gene_symbols, ["SYM1", "SYM2"])
                expected_f = np.array([[3, 6], [1, 4]], dtype=np.float32)
                np.testing.assert_array_equal(pool_f.raw_counts, expected_f)

    def test_03_normalization_target_sum_2000_log1p(self) -> None:
        raw = np.array(
            [
                [10.0, 30.0, 60.0],  # sum = 100 -> scaled = [200, 600, 1200]
                [0.0, 0.0, 0.0],  # zero-UMI cell -> must remain 0.0 (no NaN/Inf)
            ],
            dtype=np.float32,
        )
        normed, n_zero = normalize_counts_target_sum_log1p(raw, target_sum=OCCAM_TARGET_SUM)
        self.assertEqual(n_zero, 1)
        expected0 = np.log1p(np.array([200.0, 600.0, 1200.0], dtype=np.float32))
        np.testing.assert_allclose(normed[0], expected0, rtol=1e-6, atol=1e-6)
        np.testing.assert_array_equal(normed[1], np.zeros(3, dtype=np.float32))

        with self.assertRaises(ValueError):
            normalize_counts_target_sum_log1p(np.array([[-1.0, 2.0]], dtype=np.float32))

    def test_04_rhj_axis_scatter_dedup_d7_pairwise_common_axis_and_full_44728(self) -> None:
        enc_genes = [f"ENSG{i:04d}" for i in range(10)]
        rhj_axis = ["ENSG0001", "ENSG0003", "ENSG0005", "ENSG0007"]

        # Duplicate ENSG0001 (with .1 and .2 versions): D7 must keep ONLY leading column (value 10, NOT 999)
        raw_t_ens = ["ENSG0001.1", "ENSG0001.2", "ENSG0003.4", "ENSG0009.1", "ENSG9999.1"]
        raw_t = np.array([[10.0, 999.0, 20.0, 30.0, 40.0]], dtype=np.float32)  # full L = 1099.0

        x_enc_t, audit_t = scatter_rhj_6983_to_encoder(raw_t, raw_t_ens, rhj_axis, enc_genes)
        self.assertEqual(audit_t["n_present_rhj"], 2)  # ENSG0001, ENSG0003
        self.assertEqual(audit_t["n_absent_rhj"], 2)  # ENSG0005, ENSG0007
        self.assertEqual(audit_t["n_dropped_duplicate_raw_ensembl"], 1)
        # Verify leading copy (10.0) used with full-panel row sum 1099.0 (NOT axis-only sum!)
        exp_g1 = np.log1p(np.float32(10.0) / np.float32(1099.0) * np.float32(2000.0))
        exp_g3 = np.log1p(np.float32(20.0) / np.float32(1099.0) * np.float32(2000.0))
        self.assertAlmostEqual(float(x_enc_t[0, 1]), float(exp_g1), places=6)
        self.assertAlmostEqual(float(x_enc_t[0, 3]), float(exp_g3), places=6)
        # Non-axis column ENSG0009 (idx 9) and absent RHJ columns (5, 7) MUST be strictly 0.0
        for zero_col in (0, 2, 4, 5, 6, 7, 8, 9):
            self.assertEqual(float(x_enc_t[0, zero_col]), 0.0)

        # Pairwise common axis mask (`P4`): source has ENSG0001, ENSG0005, ENSG0007 (missing ENSG0003)
        raw_s_ens = ["ENSG0001", "ENSG0005", "ENSG0007", "ENSG0002"]
        raw_s = np.array([[5.0, 15.0, 25.0, 55.0]], dtype=np.float32)  # full L = 100.0

        pmask_ts = build_pairwise_common_axis_mask(
            rhj_axis, raw_t_ens, raw_s_ens, target_line="t", source_line="s"
        )
        pmask_st = build_pairwise_common_axis_mask(
            rhj_axis, raw_s_ens, raw_t_ens, target_line="s", source_line="t"
        )
        np.testing.assert_array_equal(pmask_ts.mask, pmask_st.mask)
        self.assertEqual(pmask_ts.mask_sha256, pmask_st.mask_sha256)
        self.assertEqual(pmask_ts.n_common, 1)  # only ENSG0001
        self.assertEqual(pmask_ts.n_union_absent, 3)
        self.assertTrue(pmask_ts.requires_common_axis_arm)

        x_t_cmn, aud_t_cmn = scatter_rhj_6983_to_encoder(
            raw_t, raw_t_ens, rhj_axis, enc_genes, active_rhj_mask=pmask_ts.mask
        )
        x_s_cmn, aud_s_cmn = scatter_rhj_6983_to_encoder(
            raw_s, raw_s_ens, rhj_axis, enc_genes, active_rhj_mask=pmask_ts.mask
        )
        self.assertEqual(aud_t_cmn["input_axis_mode"], "pairwise_common_axis")
        self.assertEqual(aud_t_cmn["n_active_rhj_populated"], 1)
        self.assertEqual(aud_t_cmn["n_present_but_zeroed_by_mask"], 1)
        self.assertEqual(aud_s_cmn["n_active_rhj_populated"], 1)
        self.assertEqual(aud_s_cmn["n_present_but_zeroed_by_mask"], 2)
        # Non-common RHJ genes MUST be 0.0 in BOTH target and source
        self.assertEqual(float(x_t_cmn[0, 3]), 0.0)
        self.assertEqual(float(x_s_cmn[0, 5]), 0.0)
        self.assertEqual(float(x_s_cmn[0, 7]), 0.0)
        # Common ENSG0001 retains exact full-panel-normalized value
        self.assertAlmostEqual(float(x_t_cmn[0, 1]), float(exp_g1), places=6)
        self.assertAlmostEqual(
            float(x_s_cmn[0, 1]),
            float(np.log1p(np.float32(5.0) / np.float32(100.0) * np.float32(2000.0))),
            places=6,
        )

        # Secondary `full_44728_api` normalizes strictly over encoder-matched genes (drops ENSG9999 = 40.0 and dup = 999.0)
        x_full_t, aud_full_t = scatter_full_44728_api(raw_t, raw_t_ens, enc_genes)
        self.assertEqual(aud_full_t["n_matched_encoder_genes"], 3)  # ENSG0001(10), ENSG0003(20), ENSG0009(30) -> sum=60
        exp_full_g9 = np.log1p(np.float32(30.0) / np.float32(60.0) * np.float32(2000.0))
        self.assertAlmostEqual(float(x_full_t[0, 9]), float(exp_full_g9), places=6)

        # DEDUP direct check
        keep_idx, uniq_ids, dups = dedup_leading_ids(["ENSG1.1", "", "ENSG2.1", "ENSG1.3"])
        self.assertEqual(keep_idx.tolist(), [0, 2])
        self.assertEqual(uniq_ids, ["ENSG1", "ENSG2"])
        self.assertEqual(dups, ["ENSG1"])

    def test_05_two_tier_perturbation_resolution_and_output_projection_synthetic(self) -> None:
        alias_map = {"AARS": "AARS1", "TAZ": "TAFAZZIN", "H2AFZ": "H2AZ1"}
        canon_fn = lambda s: alias_map.get(s, s)  # noqa: E731

        # Both legacy AND approved nodes exist in ckpt_pert2id with DIFFERENT IDs
        ckpt_pert2id = {
            "AARS": 10,
            "AARS1": 9010,  # untrained duplicate node!
            "TAZ": 20,
            "TAFAZZIN": 9020,  # untrained duplicate node!
            "H2AFZ": 30,
            "H2AZ1": 9030,  # untrained duplicate node!
            "TP53": 40,
            "RBM14": 50,
        }
        rhj_elig = ["AARS", "TAZ", "H2AFZ", "TP53"]  # trained legacy keys only

        res = resolve_perturbation_keys_two_tier(
            ["AARS1", "TAFAZZIN", "H2AZ1", "TP53", "RBM14"],
            rhj_elig,
            ckpt_pert2id,
            canon_fn,
        )
        self.assertEqual(res.resolved_keys, ["AARS", "TAZ", "H2AFZ", "TP53", "RBM14"])
        self.assertEqual(res.resolved_ids, [10, 20, 30, 40, 50])  # NEVER 9010/9020/9030!
        self.assertEqual(res.tiers, [1, 1, 1, 1, 2])
        self.assertEqual(res.tier2_targets, ["RBM14"])
        self.assertEqual(res.bridged_pairs, {"AARS1": "AARS", "TAFAZZIN": "TAZ", "H2AZ1": "H2AFZ"})

        # Unresolved gene raises KeyError
        with self.assertRaises(KeyError):
            resolve_perturbation_keys_two_tier(["BOGUS_GENE"], rhj_elig, ckpt_pert2id, canon_fn)

        # Collision in Tier-1 eligible keys raises ValueError
        with self.assertRaises(ValueError):
            resolve_perturbation_keys_two_tier(["AARS1"], ["AARS", "AARS1"], ckpt_pert2id, canon_fn)

        # Output gene projection
        axis = RhjGeneAxis(
            ensembl_ids=["ENSG1", "ENSG2", "ENSG3"],
            gene_symbols=["AARS", "TP53", "TAZ"],
            axis_sha256="dummy",
        )
        enc_genes = ["ENSG0", "ENSG3", "ENSG1", "ENSG2"]
        cols, p_aud = build_output_gene_projection(axis, enc_genes, ["TAFAZZIN", "AARS1"], canon_fn)
        self.assertEqual(cols.tolist(), [1, 2])
        self.assertEqual(p_aud["n_mapped"], 2)
        with self.assertRaises(KeyError):
            build_output_gene_projection(axis, enc_genes, ["MISSING_GENE"], canon_fn)

    def test_06_gate5b_panel_builder_synthetic_l1_and_determinism(self) -> None:
        n_genes_total = 12
        axis_ens = [f"ENSG{g:05d}" for g in range(8)]
        all_ens = [f"ENSG{g:05d}" for g in range(n_genes_total)]
        all_syms = [f"G{g}" for g in range(n_genes_total)]

        # 10 kept test perts (12-15 cells each), 2 shallow test perts (8 cells each), 40 controls across 4 batches
        perts_kept = [f"P{i:02d}" for i in range(10)]
        perts_shallow = ["SHALLOW_A", "SHALLOW_B"]
        eligible_test = perts_kept + perts_shallow

        obs_gene: List[str] = []
        obs_batch: List[str] = []
        rows_X: List[np.ndarray] = []
        base_ctrl = np.linspace(5.0, 20.0, n_genes_total, dtype=np.float32)

        # 40 controls
        for i in range(40):
            obs_gene.append("non-targeting")
            obs_batch.append(str((i % 4) + 1))
            rows_X.append(base_ctrl.copy())

        # 10 kept test perts with deterministic noiseless shifts -> experimental pearson_delta == 1.0!
        n_kept_cells = 0
        for idx, p in enumerate(perts_kept):
            n_c = 12 + (idx % 3)  # mix of even & odd counts
            n_kept_cells += n_c
            shift = np.zeros(n_genes_total, dtype=np.float32)
            shift[idx % 8] = 10.0 + idx
            shift[(idx + 3) % 8] = 5.0
            for j in range(n_c):
                obs_gene.append(p)
                obs_batch.append(str((j % 4) + 1))
                rows_X.append(base_ctrl + shift)

        # 2 shallow perts (< 10 cells) + 15 training-only cells (must be excluded from panel)
        for p in perts_shallow:
            for j in range(8):
                obs_gene.append(p)
                obs_batch.append(str((j % 4) + 1))
                rows_X.append(base_ctrl + 2.0)
        for j in range(15):
            obs_gene.append("TRAIN_ONLY_PERT")
            obs_batch.append("1")
            rows_X.append(base_ctrl + 1.0)

        X_all = np.vstack(rows_X).astype(np.float32)
        bcs = [f"BC_{i:04d}" for i in range(len(obs_gene))]
        expected_panel_rows = 40 + n_kept_cells

        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            fp_dense = tdp / "hepg2_synth_dense.h5ad"
            fp_csr = tdp / "hepg2_synth_csr.h5ad"
            _write_synthetic_h5ad(
                fp_dense, X_all, bcs, obs_gene, obs_batch, all_ens, all_syms, cat_mode="legacy", x_format="dense"
            )
            _write_synthetic_h5ad(
                fp_csr, X_all, bcs, obs_gene, obs_batch, all_ens, all_syms, cat_mode="legacy", x_format="csr"
            )

            panel_d = build_gate5b_hepg2_panel(
                fp_dense,
                axis_ens,
                eligible_test,
                enforce_lock_l1_counts=True,
                expected_n_perts=10,
                expected_n_controls=40,
                expected_n_panel_rows=expected_panel_rows,
                expected_cohort_sha256=None,
            )
            panel_s = build_gate5b_hepg2_panel(
                fp_csr,
                axis_ens,
                eligible_test,
                enforce_lock_l1_counts=True,
                expected_n_perts=10,
                expected_n_controls=40,
                expected_n_panel_rows=expected_panel_rows,
                expected_cohort_sha256=None,
            )
            self.assertEqual(panel_d.too_shallow, perts_shallow)
            self.assertEqual(panel_d.n_controls, 40)
            self.assertEqual(panel_d.n_panel_rows, expected_panel_rows)
            np.testing.assert_allclose(panel_d.X_panel, panel_s.X_panel, rtol=1e-6, atol=1e-6)
            self.assertEqual(
                [sp_o.pert_A_sha256 for sp_o in panel_d.splits],
                [sp_o.pert_A_sha256 for sp_o in panel_s.splits],
            )
            # Noiseless identical A/B cells -> experimental pearson_delta == 1.0 and mae == 0.0
            self.assertIsNotNone(panel_d.experimental_metrics)
            self.assertAlmostEqual(panel_d.experimental_metrics["experimental_pearson_delta_mean"], 1.0, places=6)
            self.assertAlmostEqual(panel_d.experimental_metrics["experimental_mae_mean"], 0.0, places=6)

            # Lock L1 hard halt when expected counts mismatch
            with self.assertRaises(Gate5BPanelCountError):
                build_gate5b_hepg2_panel(
                    fp_dense,
                    axis_ens,
                    eligible_test,
                    enforce_lock_l1_counts=True,
                    expected_n_perts=611,
                    expected_n_controls=4976,
                    expected_n_panel_rows=42402,
                )

            # Per-group blake2b independence: removing P09 must NOT alter P00..P08 halves
            sp_full = draw_splithalf_splits(panel_d.labels, panel_d.batches, perts_kept, n_splits=2, base_seed=42)
            sp_sub = draw_splithalf_splits(panel_d.labels, panel_d.batches, perts_kept[:-1], n_splits=2, base_seed=42)
            for sf, ss in zip(sp_full, sp_sub):
                p0_f = sf.pert_A[panel_d.labels[sf.pert_A] == "P00"]
                p0_s = ss.pert_A[panel_d.labels[ss.pert_A] == "P00"]
                np.testing.assert_array_equal(p0_f, p0_s)
                np.testing.assert_array_equal(sf.ctrl_A, ss.ctrl_A)

            # Manifest comparison exact Case A vs Case B
            synth_manifest = {
                "splits": [
                    {
                        "split": s.split_index,
                        "n_pert_per_half": len(s.pert_A),
                        "n_ctrl_per_half": len(s.ctrl_A),
                        "n_dropped_odd": s.n_dropped_odd,
                        "pert_A_sha256": s.pert_A_sha256,
                        "pert_B_sha256": s.pert_B_sha256,
                        "ctrl_A_sha256": s.ctrl_A_sha256,
                        "ctrl_B_sha256": s.ctrl_B_sha256,
                    }
                    for s in panel_d.splits
                ]
            }
            cmp_a = compare_splits_against_manifest(panel_d.splits, synth_manifest)
            self.assertEqual(cmp_a["tolerance_case"], "CASE_A_EXACT_SPLIT_MATCH")
            self.assertEqual(cmp_a["n_hash_matches"], 20)

            synth_manifest["splits"][0]["pert_A_sha256"] = "0" * 64
            cmp_b = compare_splits_against_manifest(panel_d.splits, synth_manifest)
            self.assertEqual(cmp_b["tolerance_case"], "CASE_B_FRESH_DETERMINISTIC_SPLITS")
            self.assertEqual(cmp_b["n_hash_matches"], 19)

            # Cross-process PYTHONHASHSEED determinism check
            hashes_by_seed = []
            for hseed in ("1", "99991"):
                env = dict(os.environ)
                env["PYTHONHASHSEED"] = hseed
                script = (
                    f"import sys, json; sys.path.insert(0, {str(ADDB_SRC_DIR)!r}); "
                    f"from pathlib import Path; "
                    f"from gate5b_panel_addb import build_gate5b_hepg2_panel; "
                    f"p = build_gate5b_hepg2_panel(Path({str(fp_dense)!r}), {axis_ens!r}, {eligible_test!r}, "
                    f"load_expression=False, enforce_lock_l1_counts=False); "
                    f"print(json.dumps([s.pert_A_sha256 + s.ctrl_B_sha256 for s in p.splits]))"
                )
                r = subprocess.run([sys.executable, "-I", "-c", script], capture_output=True, text=True, env=env)
                self.assertEqual(r.returncode, 0, msg=r.stderr)
                hashes_by_seed.append(json.loads(r.stdout.strip()))
            self.assertEqual(hashes_by_seed[0], hashes_by_seed[1])

    def test_07_export_z_ctrl_and_deltas_manifests_and_tolerances(self) -> None:
        rng = np.random.default_rng(707)
        bcs = [f"BC_{i:03d}" for i in range(16)]
        z_a = rng.standard_normal((16, 512)).astype(np.float32)
        z_b_ok = (z_a + 1e-6).astype(np.float32)
        z_b_bad = (z_a + 1e-3).astype(np.float32)

        d_a = rng.standard_normal((5, 8)).astype(np.float64)
        d_b_ok = d_a + 1e-5
        d_b_bad = d_a + 1e-2

        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            man_z_a = save_occam_z_ctrl_latents(
                tdp / "lat_a",
                {
                    "a549_ko": {
                        "Z_ctrl": z_a,
                        "control_barcodes": bcs,
                        "seed": 20261601,
                        "input_axis_mode": "rhj_6983_axis",
                        "n_active_rhj_genes": 6950,
                    }
                },
            )
            self.assertIn("z_ctrl_a549_ko.npz", man_z_a)
            ver_z = verify_exported_manifest(tdp / "lat_a", "z_ctrl_manifest_sha256.json")
            self.assertTrue(ver_z["ok"])

            save_occam_z_ctrl_latents(
                tdp / "lat_b_ok",
                {"a549_ko": {"Z_ctrl": z_b_ok, "control_barcodes": bcs, "seed": 20261601}},
            )
            save_occam_z_ctrl_latents(
                tdp / "lat_b_bad",
                {"a549_ko": {"Z_ctrl": z_b_bad, "control_barcodes": bcs, "seed": 20261601}},
            )
            save_occam_z_ctrl_latents(
                tdp / "lat_b_perm",
                {"a549_ko": {"Z_ctrl": z_a, "control_barcodes": list(reversed(bcs)), "seed": 20261601}},
            )

            ok_z, _ = compare_exported_npz_arrays(
                tdp / "lat_a" / "z_ctrl_a549_ko.npz",
                tdp / "lat_b_ok" / "z_ctrl_a549_ko.npz",
                array_key="Z_ctrl",
                atol=Z_CTRL_CROSS_PARTY_ATOL,
            )
            self.assertTrue(ok_z)
            bad_z, _ = compare_exported_npz_arrays(
                tdp / "lat_a" / "z_ctrl_a549_ko.npz",
                tdp / "lat_b_bad" / "z_ctrl_a549_ko.npz",
                array_key="Z_ctrl",
                atol=Z_CTRL_CROSS_PARTY_ATOL,
            )
            self.assertFalse(bad_z)
            perm_z, _ = compare_exported_npz_arrays(
                tdp / "lat_a" / "z_ctrl_a549_ko.npz",
                tdp / "lat_b_perm" / "z_ctrl_a549_ko.npz",
                array_key="Z_ctrl",
                atol=Z_CTRL_CROSS_PARTY_ATOL,
            )
            self.assertFalse(perm_z)

            # Predicted deltas export across primary + pairwise_common_axis + real baseline
            t_names = [f"T{i}" for i in range(5)]
            g_names = [f"G{j}" for j in range(8)]
            p_main = build_occam_delta_payload(
                delta_pred=d_a,
                target_names=t_names,
                gene_names=g_names,
                control_source="a549_ko",
                control_barcodes=bcs,
                input_axis_mode="rhj_6983_axis",
                baseline_mode="reconstructed",
                resolved_pert_keys=t_names,
                resolved_pert_ids=[1, 2, 3, 4, 5],
            )
            p_cmn = build_occam_delta_payload(
                delta_pred=d_a,
                target_names=t_names,
                gene_names=g_names,
                control_source="a549_ko",
                control_barcodes=bcs,
                input_axis_mode="pairwise_common_axis",
                baseline_mode="reconstructed",
                pair_counterpart="cd4t_rest",
                n_active_rhj_genes=6557,
            )
            p_real = build_occam_delta_payload(
                delta_pred=d_b_ok,
                target_names=t_names,
                gene_names=g_names,
                control_source="a549_ko",
                control_barcodes=bcs,
                input_axis_mode="full_44728_api",
                baseline_mode="real",
            )
            man_d = save_occam_line_deltas(
                tdp / "arrays_a",
                "a549_ko",
                {
                    "Occam_rhj_6983_axis_reconstructed_matched": p_main,
                    "Occam_pairwise_common_vs_cd4t_rest_matched": p_cmn,
                    "Occam_full_44728_api_real_matched": p_real,
                },
            )
            self.assertEqual(len(man_d), 3)
            self.assertTrue(verify_exported_manifest(tdp / "arrays_a" / "a549_ko", "array_manifest_sha256.json")["ok"])

            save_occam_line_deltas(
                tdp / "arrays_b",
                "a549_ko",
                {
                    "Occam_rhj_6983_axis_reconstructed_matched": build_occam_delta_payload(
                        delta_pred=d_b_ok,
                        target_names=t_names,
                        gene_names=g_names,
                        control_source="a549_ko",
                        control_barcodes=bcs,
                    ),
                    "Occam_bad": build_occam_delta_payload(
                        delta_pred=d_b_bad,
                        target_names=t_names,
                        gene_names=g_names,
                        control_source="a549_ko",
                        control_barcodes=bcs,
                    ),
                },
            )
            ok_d, _ = compare_exported_npz_arrays(
                tdp / "arrays_a" / "a549_ko" / "preds_Occam_rhj_6983_axis_reconstructed_matched.npz",
                tdp / "arrays_b" / "a549_ko" / "preds_Occam_rhj_6983_axis_reconstructed_matched.npz",
                array_key="delta_pred",
                atol=DELTA_CROSS_PARTY_ATOL,
                check_ordered_keys=("control_barcodes", "target_names", "gene_names"),
            )
            self.assertTrue(ok_d)
            bad_d, _ = compare_exported_npz_arrays(
                tdp / "arrays_a" / "a549_ko" / "preds_Occam_rhj_6983_axis_reconstructed_matched.npz",
                tdp / "arrays_b" / "a549_ko" / "preds_Occam_bad.npz",
                array_key="delta_pred",
                atol=DELTA_CROSS_PARTY_ATOL,
                check_ordered_keys=("control_barcodes", "target_names", "gene_names"),
            )
            self.assertFalse(bad_d)

            # Negative validation guards
            with self.assertRaises(ValueError):
                build_occam_delta_payload(
                    delta_pred=d_a,
                    target_names=t_names,
                    gene_names=g_names,
                    control_source="a549_ko",
                    control_barcodes=bcs,
                    input_axis_mode="pairwise_common_axis",
                    pair_counterpart="",
                )
            d_nan = d_a.copy()
            d_nan[0, 0] = np.nan
            with self.assertRaises(ValueError):
                build_occam_delta_payload(
                    delta_pred=d_nan,
                    target_names=t_names,
                    gene_names=g_names,
                    control_source="a549_ko",
                    control_barcodes=bcs,
                )
            with self.assertRaises(ValueError):
                save_occam_z_ctrl_latents(
                    tdp / "lat_bad_dim",
                    {"bad": {"Z_ctrl": np.zeros((4, 511), dtype=np.float32), "control_barcodes": bcs[:4], "seed": 1}},
                )

    def test_08_public_text_sidecars_and_locked_r16_barcodes(self) -> None:
        axis_csv = _find_public_sidecar(["gene_axis.csv", "rhj_gene_axis.csv"])
        perts_csv = _find_public_sidecar(["perturbations.csv", "rhj_perturbations.csv"])
        master_csv = _find_public_sidecar(["master_gene_set_sorted.csv", "occam_master_gene_set.csv"])
        cohort_json = _find_public_sidecar(
            ["cohort.json", "hf_benchmarks__rhj_fewshot__splithalf__cohort.json"]
        )
        if not (axis_csv and perts_csv and master_csv and cohort_json):
            self.skipTest("Small public text sidecars not present locally.")

        # 1. RHJ gene axis hash & count
        axis = load_rhj_gene_axis(axis_csv, verify_sha256=True)
        self.assertEqual(len(axis.ensembl_ids), RHJ_N_AXIS_GENES)
        self.assertEqual(axis.axis_sha256, RHJ_GENE_AXIS_SHA256)

        # 2. RHJ perturbations & cohort.json 621 vs 611 exact check
        ptab = load_rhj_perturbations_table(perts_csv)
        self.assertEqual(len(ptab["eligible_keys"]), 2381)
        self.assertEqual(len(ptab["test_keys"]), 621)
        self.assertEqual(len(ptab["train_keys"]), 1657)
        self.assertEqual(len(ptab["val_keys"]), 103)

        cohort = json.loads(cohort_json.read_text(encoding="utf-8"))
        self.assertEqual(cohort["n_perts"], 611)
        self.assertEqual(cohort["n_controls"], 4976)
        self.assertEqual(cohort["cohort_sha256"], GATE5B_EXPECTED_COHORT_SHA256)
        shallow = sorted(set(ptab["test_keys"]) - set(cohort["perts"]))
        self.assertEqual(tuple(shallow), GATE5B_EXPECTED_TOO_SHALLOW)

        # 3. Master pert2id (18,216 nodes) + frozen R10 HGNCBridge on locked R16 targets/genes
        master_pert2id = load_master_pert2id_csv(master_csv)
        self.assertEqual(len(master_pert2id), 18216)

        from data_loader_r10 import DEFAULT_HGNC_FILE, HGNCBridge  # noqa: E402

        bridge = HGNCBridge.from_tsv(DEFAULT_HGNC_FILE)
        line_expectations = {
            "a549_ko": (
                REPO_ROOT / "round16" / "results" / "arrays" / "a549_ko" / "observed.npz",
                402,
                1,
                14,
                1837,
                ["RBM14"],
            ),
            "dld1": (
                REPO_ROOT / "round16" / "results" / "arrays" / "dld1" / "observed.npz",
                75,
                0,
                3,
                1832,
                [],
            ),
            "cd4t_rest": (
                REPO_ROOT / "round16" / "results" / "arrays" / "cd4t_rest" / "observed.npz",
                245,
                0,
                1,
                1830,
                [],
            ),
            "feng_ipsc": (
                REPO_ROOT / "round16" / "results" / "sensitivity" / "arrays" / "feng_ipsc" / "observed.npz",
                88,
                0,
                0,
                1575,
                [],
            ),
        }
        # Build synthetic encoder gene list containing all 6,983 RHJ genes cleanly
        enc_genes_synth = list(axis.ensembl_ids)
        for line_name, (obs_path, exp_t1, exp_t2, exp_bridged, exp_g, exp_t2_list) in line_expectations.items():
            if not obs_path.is_file():
                continue
            z = np.load(obs_path, allow_pickle=False)
            t_names = [str(x) for x in z["target_names"].tolist()]
            g_names = [str(x) for x in z["gene_names"].tolist()]
            res = resolve_perturbation_keys_two_tier(
                t_names, ptab["eligible_keys"], master_pert2id, bridge.canonicalize, require_all_resolved=True
            )
            self.assertEqual(len(res.tier1_targets), exp_t1, msg=line_name)
            self.assertEqual(len(res.tier2_targets), exp_t2, msg=line_name)
            self.assertEqual(res.tier2_targets, exp_t2_list, msg=line_name)
            self.assertEqual(len(res.bridged_pairs), exp_bridged, msg=line_name)
            # Confirm bridged keys map to the trained rhj_fewshot pert2id, never a distinct ID
            for t_canon, occam_k in res.bridged_pairs.items():
                self.assertEqual(master_pert2id[occam_k], ptab["pert2id"][occam_k], msg=f"{line_name}:{t_canon}")

            proj_cols, p_aud = build_output_gene_projection(
                axis, enc_genes_synth, g_names, bridge.canonicalize, require_full_coverage=True
            )
            self.assertEqual(len(proj_cols), exp_g, msg=line_name)
            self.assertEqual(p_aud["n_missing"], 0, msg=line_name)

        # 4. Locked R16 control barcodes JSON check (read-only)
        if DEFAULT_PRIMARY_BARCODES_JSON.is_file() and DEFAULT_SENSITIVITY_BARCODES_JSON.is_file():
            for ln in ("a549_ko", "dld1", "cd4t_rest", "feng_ipsc", "nadig_hepg2", "nadig_jurkat"):
                bcs, seed, _ = load_locked_r16_barcodes(ln, pool="primary")
                self.assertEqual(len(bcs), 512, msg=ln)
                self.assertEqual(seed, 20261601, msg=ln)
            for ln in ("a549_ko", "dld1", "cd4t_rest", "feng_ipsc"):
                rbcs, rseed, raud = load_locked_r16_barcodes(ln, pool="resample")
                self.assertEqual(len(rbcs), 512, msg=ln)
                self.assertEqual(rseed, 20261602, msg=ln)
                self.assertEqual(raud["overlap_with_primary"], 0, msg=ln)

    def test_09_model_free_modules_zero_occam_imports(self) -> None:
        """Guardrail: verify zero `import occam` / `from occam` in all model-free modules."""
        model_free_names = (
            "__init__.py",
            "raw_controls_addb.py",
            "preprocess_addb.py",
            "gate5b_panel_addb.py",
            "export_addb.py",
            "test_unit_addb.py",
        )
        for name in model_free_names:
            fp = ADDB_SRC_DIR / name
            self.assertTrue(fp.is_file(), msg=f"Missing model-free module {name}")
            tree = ast.parse(fp.read_text(encoding="utf-8"), filename=str(fp))
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    for alias in node.names:
                        self.assertFalse(
                            alias.name == "occam" or alias.name.startswith("occam."),
                            msg=f"Forbidden Occam import {alias.name!r} in {fp.name}",
                        )
                elif isinstance(node, ast.ImportFrom):
                    mod = node.module or ""
                    self.assertFalse(
                        mod == "occam" or mod.startswith("occam."),
                        msg=f"Forbidden Occam import from {mod!r} in {fp.name}",
                    )


if __name__ == "__main__":
    unittest.main(verbosity=2)
