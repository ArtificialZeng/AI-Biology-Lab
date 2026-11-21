#!/usr/bin/env python3
"""
Unit & end-to-end integration tests for the VCC 2026 prospective prediction pipeline:
  - `vcc2026_format.py` (format validation, dispersed integer count emission, `.vcc` packaging,
    SHA-256 timestamp bundle, upload -> cell-eval2 scoring control-attachment bridge)
  - `prospective_predict.py` (standalone `.h5ad` source mode + `p2/shared/lines/` Replogle dev mode)
"""

from __future__ import annotations

import json
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

PROSPECTIVE_DIR = Path(__file__).resolve().parent
if str(PROSPECTIVE_DIR) not in sys.path:
    sys.path.insert(0, str(PROSPECTIVE_DIR))

from prospective_predict import build_parser, run_prospective
from vcc2026_format import (
    DEFAULT_CONTEXT_COL,
    DEFAULT_NTC_NAME,
    DEFAULT_PERT_COL,
    attach_real_controls_for_celleval2,
    emit_dispersed_counts_for_context,
    load_gene_names,
    load_pert_counts,
    package_vcc_archive,
    sha256_file,
    validate_vcc2026_anndata,
)


def _make_synthetic_raw_h5ad(
    path: Path,
    genes: list[str],
    perts: list[str],
    n_ctrl: int = 80,
    n_per_pert: int = 20,
    seed: int = 0,
    include_perturbed: bool = True,
) -> None:
    rng = np.random.default_rng(seed)
    n_genes = len(genes)
    base_rates = rng.gamma(4.0, 5.0, size=n_genes)

    blocks = []
    labels = []
    X_ctrl = rng.poisson(base_rates[None, :], size=(n_ctrl, n_genes)).astype(np.float32)
    blocks.append(X_ctrl)
    labels.extend([DEFAULT_NTC_NAME] * n_ctrl)

    if include_perturbed:
        for j, p in enumerate(perts):
            mult = np.ones(n_genes, dtype=np.float64)
            # Perturb 16 downstream genes strongly enough to clear cell-eval2's >=10 DE gene gate
            for k in range(16):
                mult[(j * 7 + k) % n_genes] = 3.5
            if p in genes:
                mult[genes.index(p)] = 0.05
            X_p = rng.poisson(base_rates[None, :] * mult[None, :], size=(n_per_pert, n_genes)).astype(np.float32)
            blocks.append(X_p)
            labels.extend([p] * n_per_pert)

    X_all = sp.csr_matrix(np.vstack(blocks), dtype=np.float32)
    X_all.eliminate_zeros()
    obs = pd.DataFrame({DEFAULT_PERT_COL: labels}, index=[f"c_{i}" for i in range(X_all.shape[0])])
    var = pd.DataFrame(index=genes)
    ad.AnnData(X=X_all, obs=obs, var=var).write_h5ad(path)


class TestVCC2026Format(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.genes = [f"G{i:03d}" for i in range(60)]
        self.perts = ["G000", "G001", "G002"]
        self.genes_csv = self.root / "gene_names.csv"
        self.genes_csv.write_text("gene_name\n" + "\n".join(self.genes) + "\n", encoding="utf-8")
        self.perts_csv = self.root / "pert_counts.csv"
        self.perts_csv.write_text("target_gene\n" + "\n".join(self.perts) + "\n", encoding="utf-8")

    def tearDown(self) -> None:
        self.tmp.cleanup()

    def test_load_gene_names_and_pert_counts(self) -> None:
        g = load_gene_names(self.genes_csv)
        self.assertEqual(g, self.genes)
        pc = load_pert_counts(self.perts_csv, default_n_cells=400)
        self.assertEqual(pc[None], {"G000": 400, "G001": 400, "G002": 400})

    def test_emit_dispersed_counts_enforces_cap_and_integrality(self) -> None:
        rng = np.random.default_rng(7)
        # High counts to test the 1,000,000 per-cell cap enforcement
        ctrl = sp.csr_matrix(rng.integers(500, 1500, size=(10, 60)).astype(np.float32))
        ratios = np.full((2, 60), 25.0, dtype=np.float64)
        X_out = emit_dispersed_counts_for_context(
            ctrl, ["G000", "G001"], ratios, cells_per_pert=12, max_counts_per_cell=100_000, seed=1
        )
        self.assertTrue(sp.issparse(X_out))
        self.assertEqual(X_out.shape, (24, 60))
        self.assertTrue(np.all(X_out.data > 0))  # no explicit zeros
        self.assertTrue(np.all(np.abs(X_out.data - np.trunc(X_out.data)) < 1e-6))
        totals = np.asarray(X_out.sum(axis=1)).ravel()
        self.assertTrue(np.all(totals <= 100_000))
        self.assertTrue(np.all(totals > 0))

    def _build_valid_adata(self, n_per_pert: int = 8) -> ad.AnnData:
        rng = np.random.default_rng(123)
        ctrl = sp.csr_matrix(rng.poisson(5.0, size=(20, len(self.genes))).astype(np.float32))
        ratios = np.ones((len(self.perts), len(self.genes)), dtype=np.float64)
        blocks, obs_p, obs_c = [], [], []
        for c_idx, ctx in enumerate(["D", "E", "F"]):
            blocks.append(
                emit_dispersed_counts_for_context(
                    ctrl, self.perts, ratios, cells_per_pert=n_per_pert, seed=10 + c_idx
                )
            )
            for p in self.perts:
                obs_p.extend([p] * n_per_pert)
                obs_c.extend([ctx] * n_per_pert)
        X = sp.vstack(blocks, format="csr")
        obs = pd.DataFrame({DEFAULT_PERT_COL: obs_p, DEFAULT_CONTEXT_COL: obs_c},
                           index=[str(i) for i in range(X.shape[0])])
        return ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=self.genes))

    def test_validation_catches_all_error_modes(self) -> None:
        pc = load_pert_counts(self.perts_csv, default_n_cells=8)
        base_kw = dict(
            gene_names=self.genes, pert_counts=pc, required_contexts=("D", "E", "F"),
            cells_per_pert=8, expected_gene_dim=None, max_cell_dim=None,
        )
        ok_ad = self._build_valid_adata(n_per_pert=8)
        rep_ok = validate_vcc2026_anndata(ok_ad, **base_kw)
        self.assertTrue(rep_ok.passed, msg=rep_ok.errors)

        # 1. Fractional values
        bad = ok_ad.copy()
        bad.X.data[0] = 1.5
        self.assertFalse(validate_vcc2026_anndata(bad, **base_kw).passed)

        # 2. Negative values
        bad = ok_ad.copy()
        bad.X.data[0] = -1.0
        self.assertFalse(validate_vcc2026_anndata(bad, **base_kw).passed)

        # 3. non-targeting control rows included
        bad = ok_ad.copy()
        bad.obs.iloc[0, bad.obs.columns.get_loc(DEFAULT_PERT_COL)] = DEFAULT_NTC_NAME
        self.assertFalse(validate_vcc2026_anndata(bad, **base_kw).passed)

        # 4. Guide suffix instead of gene symbol
        bad = ok_ad.copy()
        col_i = bad.obs.columns.get_loc(DEFAULT_PERT_COL)
        mask_g0 = bad.obs[DEFAULT_PERT_COL] == "G000"
        bad.obs.iloc[np.where(mask_g0)[0][0], col_i] = "G000-1"
        rep_guide = validate_vcc2026_anndata(bad, **base_kw)
        self.assertFalse(rep_guide.passed)
        self.assertTrue(any("guide" in e.lower() for e in rep_guide.errors))

        # 5. Missing context
        bad = ok_ad[ok_ad.obs[DEFAULT_CONTEXT_COL] != "F"].copy()
        self.assertFalse(validate_vcc2026_anndata(bad, **base_kw).passed)

        # 6. Dense matrix rejected
        bad = ok_ad.copy()
        bad.X = bad.X.toarray()
        self.assertFalse(validate_vcc2026_anndata(bad, **base_kw).passed)

        # 7. Per-cell sum > cap rejected
        bad = ok_ad.copy()
        bad.X.data[0] = 2_000_000.0
        self.assertFalse(validate_vcc2026_anndata(bad, **base_kw).passed)

        # 8. Out-of-order genes rejected
        bad = ok_ad[:, list(reversed(self.genes))].copy()
        self.assertFalse(validate_vcc2026_anndata(bad, **base_kw).passed)

    def test_vcc_archive_packaging(self) -> None:
        ok_ad = self._build_valid_adata(n_per_pert=8)
        vcc_path = self.root / "out.vcc"
        nnz = package_vcc_archive(
            ok_ad, vcc_path, contexts=["D", "E", "F"], perts_path=self.perts_csv, panel_id="test-panel"
        )
        self.assertEqual(nnz, int(ok_ad.X.nnz))
        with tarfile.open(vcc_path, "r:*") as tar:
            names = tar.getnames()
            self.assertIn("meta.json", names)
            self.assertIn("pred.h5ad.zst", names)
            meta = json.loads(tar.extractfile("meta.json").read().decode("utf-8"))
            self.assertEqual(meta["contexts"], ["D", "E", "F"])
            self.assertEqual(meta["nnz"], nnz)
            self.assertEqual(meta["panel_id"], "test-panel")


class TestProspectiveEndToEnd(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self) -> None:
        self.tmp.cleanup()

    def test_end_to_end_bundle_and_celleval2_bridge(self) -> None:
        genes = [f"GENE{i:03d}" for i in range(80)]
        perts = ["GENE000", "GENE001", "GENE002", "GENE003", "GENE004"]

        # Create 2 source lines (Replogle/Nadig-style raw h5ads)
        src1 = self.root / "src_K562.h5ad"
        src2 = self.root / "src_RPE1.h5ad"
        _make_synthetic_raw_h5ad(src1, genes, perts, n_ctrl=60, n_per_pert=20, seed=1, include_perturbed=True)
        _make_synthetic_raw_h5ad(src2, genes, perts, n_ctrl=60, n_per_pert=20, seed=2, include_perturbed=True)

        # Create target bundle (D, E, F) with unperturbed controls only + manifest + gene_names + pert_counts
        bdir = self.root / "bundle"
        bdir.mkdir()
        (bdir / "manifest.json").write_text(
            json.dumps({"contexts": ["D", "E", "F"], "partition": "test", "panel_id": "sim-test"}),
            encoding="utf-8",
        )
        (bdir / "gene_names.csv").write_text("gene_name\n" + "\n".join(genes) + "\n", encoding="utf-8")
        (bdir / "pert_counts.csv").write_text("target_gene\n" + "\n".join(perts) + "\n", encoding="utf-8")

        for idx, ctx in enumerate(["D", "E", "F"]):
            _make_synthetic_raw_h5ad(
                bdir / f"context_{ctx}.h5ad", genes, perts, n_ctrl=60, n_per_pert=10,
                seed=10 + idx, include_perturbed=False
            )

        outdir = self.root / "prospective_out"
        args = build_parser().parse_args([
            "--bundle-dir", str(bdir),
            "--source-h5ad", f"K562={src1}",
            "--source-h5ad", f"RPE1={src2}",
            "--method", "Hybrid_ProjShared",
            "--cells-per-pert", "30",
            "--simulation",
            "--package-vcc",
            "--outdir", str(outdir),
            "--seed", "42",
        ])
        manifest = run_prospective(args)
        self.assertTrue(manifest["validation"]["passed"])
        self.assertEqual(manifest["validation"]["n_cells"], 3 * len(perts) * 30)
        self.assertEqual(manifest["validation"]["n_genes"], len(genes))
        self.assertTrue((outdir / "prediction.h5ad").exists())
        self.assertTrue((outdir / "prediction.prep.vcc").exists())
        self.assertTrue((outdir / "SHA256SUMS.txt").exists())
        self.assertTrue((outdir / "TIMESTAMP_DECLARATION.md").exists())

        # Verify SHA256SUMS.txt independently
        sums_lines = [ln.strip() for ln in (outdir / "SHA256SUMS.txt").read_text().splitlines() if ln.strip()]
        for ln in sums_lines:
            digest, fname = ln.split()
            self.assertEqual(digest, sha256_file(outdir / fname))

        # Verify bridging context D to cell-eval2 scoring against a synthetic real context D
        real_d = self.root / "real_D.h5ad"
        _make_synthetic_raw_h5ad(real_d, genes, perts, n_ctrl=40, n_per_pert=30, seed=99, include_perturbed=True)
        bridged_d = self.root / "pred_D_for_ce2.h5ad"
        attach_real_controls_for_celleval2(
            outdir / "prediction.h5ad",
            real_d,
            bridged_d,
            context_col=DEFAULT_CONTEXT_COL,
            context_value="D",
        )
        from dataclasses import replace
        from cell_eval2 import EvalConfig, aggregate_metrics, compute_metrics

        cfg = replace(
            EvalConfig.from_preset("vcc2026"),
            pert_col="target_gene",
            device="cpu",
        )
        cfg = replace(cfg, de=replace(cfg.de, backend="pdex"))
        df = compute_metrics(str(bridged_d), str(real_d), config=cfg)
        agg = aggregate_metrics(df)
        self.assertEqual(len(agg), 10)  # 6 scored + 4 expr diagnostics

    def test_shared_lines_loader_and_prereg_gate(self) -> None:
        shared_dir = self.root / "shared_lines"
        genes = [f"G{i:03d}" for i in range(60)]
        perts = ["G000", "G001", "G002"]
        rng = np.random.default_rng(5)

        for ln in ["Replogle_K562", "Replogle_RPE1"]:
            ldir = shared_dir / ln
            ldir.mkdir(parents=True)
            (ldir / "genes.txt").write_text("\n".join(genes) + "\n", encoding="utf-8")
            (ldir / "targets.txt").write_text("\n".join(perts) + "\n", encoding="utf-8")
            np.savez_compressed(ldir / "ctrl_cells.npz", X_log1p=rng.normal(2.0, 0.4, size=(40, 60)).astype(np.float32))
            np.savez_compressed(
                ldir / "pseudobulk_delta.npz",
                delta_log1p=rng.normal(0.0, 0.3, size=(3, 60)).astype(np.float32),
                reproducible=np.array([True, True, True]),
            )

        ctrl_h5ad = self.root / "target_ctrl.h5ad"
        _make_synthetic_raw_h5ad(ctrl_h5ad, genes, perts, n_ctrl=40, include_perturbed=False, seed=77)
        targets_txt = self.root / "targets.txt"
        targets_txt.write_text("\n".join(perts) + "\n", encoding="utf-8")

        outdir = self.root / "out_shared"
        args = build_parser().parse_args([
            "--controls", f"SIM={ctrl_h5ad}",
            "--perts", str(targets_txt),
            "--source-lines-dir", str(shared_dir),
            "--method", "Hybrid_ProjShared",
            "--cells-per-pert", "15",
            "--simulation",
            "--outdir", str(outdir),
        ])
        manifest = run_prospective(args)
        self.assertTrue(manifest["validation"]["passed"])
        self.assertEqual(manifest["method"]["source_mode"], "shared_lines")
        self.assertEqual(manifest["method"]["source_lines"], ["Replogle_K562", "Replogle_RPE1"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
