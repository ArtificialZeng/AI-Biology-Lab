#!/usr/bin/env python3
"""Unit tests on synthetic fixtures for clean-room Addendum B adapter."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

AUDIT_DIR = Path(__file__).resolve().parent
if str(AUDIT_DIR) not in sys.path:
    sys.path.insert(0, str(AUDIT_DIR))

from audit_addb_occam_adapter import (  # noqa: E402
    DEFAULT_GENE_AXIS_CSV,
    DEFAULT_PERTURBATIONS_CSV,
    OCCAM_NUM_ENCODER_GENES,
    TwoTierPertResolver,
    compute_celleval_pseudobulk_metrics,
    compute_pairwise_common_rhj_mask,
    dedup_leading_ensembl_positions,
    load_hgnc_bridge,
    load_rhj_gene_axis,
    normalize_counts_full_44728_api,
    normalize_counts_to_rhj_44728,
    read_obs_or_var_col_legacy,
    save_delta_bundle_npz,
    strip_ensembl_version,
)


class TestAddBAdapterUnit(unittest.TestCase):
    def test_legacy_hdf5_decoding_and_dedup(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "legacy.h5ad"
            with h5py.File(str(p), "w") as f:
                obs = f.create_group("obs")
                cat_grp = obs.create_group("__categories")
                cat_grp.create_dataset("gene", data=np.asarray([b"non-targeting", b"TP53"]))
                obs.create_dataset("gene", data=np.asarray([0, 1, -1, 1], dtype=np.int32))
                obs.create_dataset("gem_group", data=np.asarray([b"b1", b"b2", b"b1", b"b3"]))

                var = f.create_group("var")
                mod_grp = var.create_group("gene_id")
                mod_grp.attrs["encoding-type"] = "categorical"
                mod_grp.create_dataset("categories", data=np.asarray([b"ENSG001.2", b"ENSG001.5"]))
                mod_grp.create_dataset("codes", data=np.asarray([0, 1], dtype=np.int32))

            with h5py.File(str(p), "r") as f:
                genes = read_obs_or_var_col_legacy(f["obs"], "gene")
                batches = read_obs_or_var_col_legacy(f["obs"], "gem_group")
                v_ids = read_obs_or_var_col_legacy(f["var"], "gene_id")

            self.assertEqual(genes, ["non-targeting", "TP53", "", "TP53"])
            self.assertEqual(batches, ["b1", "b2", "b1", "b3"])
            self.assertEqual([strip_ensembl_version(x) for x in v_ids], ["ENSG001", "ENSG001"])
            dedup = dedup_leading_ensembl_positions(v_ids)
            self.assertEqual(dedup, {"ENSG001": 0})

    def test_rhj_normalization_and_pairwise_mask(self) -> None:
        axis_ens, _ = load_rhj_gene_axis(DEFAULT_GENE_AXIS_CSV)
        enc_names = list(axis_ens) + [f"EXTRA_{i}" for i in range(OCCAM_NUM_ENCODER_GENES - len(axis_ens))]
        raw_ids = [axis_ens[0] + ".1", axis_ens[1], axis_ens[0] + ".9", "UNMAPPED_EXTRA"]
        raw = np.asarray(
            [
                [10.0, 20.0, 999.0, 70.0],  # full raw sum = 1099.0 (NOT 30.0)
                [0.0, 0.0, 0.0, 0.0],
            ],
            dtype=np.float32,
        )
        X_44728, X_6983, audit = normalize_counts_to_rhj_44728(
            raw, raw_ids, axis_ens, enc_names, target_sum=2000.0
        )
        self.assertEqual(X_44728.shape, (2, OCCAM_NUM_ENCODER_GENES))
        self.assertEqual(X_6983.shape, (2, 6983))
        self.assertEqual(audit["n_rhj_present_in_raw"], 2)
        expected_0 = np.float32(np.log1p(np.float32(10.0) / np.float32(1099.0) * np.float32(2000.0)))
        self.assertAlmostEqual(float(X_6983[0, 0]), float(expected_0), places=6)
        self.assertEqual(float(X_6983[1, 0]), 0.0)

        # Pairwise common mask drops axis_ens[1] when partner lacks it
        mask = compute_pairwise_common_rhj_mask(raw_ids, [axis_ens[0]], axis_ens)
        self.assertEqual(int(mask.sum()), 1)
        X_44728_m, X_6983_m, audit_m = normalize_counts_to_rhj_44728(
            raw, raw_ids, axis_ens, enc_names, target_sum=2000.0, active_rhj_mask=mask
        )
        self.assertEqual(float(X_6983_m[0, 1]), 0.0)
        self.assertEqual(audit_m["n_rhj_masked_out"], 1)
        self.assertEqual(float(X_44728_m[0, 0]), float(expected_0))

        # Full 44728 API excludes unmapped raw gene from library sum
        X_api, api_audit = normalize_counts_full_44728_api(raw, raw_ids, enc_names, target_sum=2000.0)
        self.assertEqual(api_audit["n_matched_44728"], 2)
        expected_api_0 = np.float32(np.log1p(np.float32(10.0) / np.float32(30.0) * np.float32(2000.0)))
        self.assertAlmostEqual(float(X_api[0, 0]), float(expected_api_0), places=6)

    def test_two_tier_pert_resolver(self) -> None:
        hb = load_hgnc_bridge()
        # Construct minimal mock pert2id containing rhj_2381 eligible keys + RBM14
        import csv

        with open(DEFAULT_PERTURBATIONS_CSV, "r", encoding="utf-8") as fh:
            elig = [r["gene_name"] for r in csv.DictReader(fh) if int(r["eligible"]) == 1]
        mock_pert2id = {k: i for i, k in enumerate(elig)}
        mock_pert2id["RBM14"] = 99999
        resolver = TwoTierPertResolver(mock_pert2id, DEFAULT_PERTURBATIONS_CSV, hgnc_bridge=hb)

        for sym, want_key in (("AARS1", "AARS"), ("DARS1", "DARS"), ("TAFAZZIN", "TAZ"), ("H2AZ1", "H2AFZ")):
            k, pid, tier = resolver.resolve(sym)
            self.assertEqual(k, want_key)
            self.assertEqual(tier, "tier1_rhj2381")
            self.assertEqual(pid, mock_pert2id[want_key])

        k2, pid2, tier2 = resolver.resolve("RBM14")
        self.assertEqual((k2, pid2, tier2), ("RBM14", 99999, "tier2_direct"))
        with self.assertRaises(KeyError):
            resolver.resolve("DOES_NOT_EXIST_XYZ")

    def test_pseudobulk_metrics_and_npz_roundtrip(self) -> None:
        real_ctrl = np.asarray([[1.0, 2.0, 3.0], [1.0, 2.0, 3.0]], dtype=np.float32)
        pred_ctrl = np.asarray([[0.5, 1.5, 2.5]], dtype=np.float32)
        real_pert = np.asarray([[2.0, 4.0, 3.0], [2.0, 4.0, 3.0]], dtype=np.float32)
        pred_pert = np.asarray([[1.5, 3.5, 2.5]], dtype=np.float32)
        m = compute_celleval_pseudobulk_metrics(
            real_pert, np.asarray(["P1", "P1"]), real_ctrl,
            pred_pert, np.asarray(["P1"]), pred_ctrl, ["P1"],
        )
        self.assertAlmostEqual(m["pearson_delta_mean"], 1.0, places=6)
        self.assertAlmostEqual(m["mae_mean"], 0.5, places=6)

        with tempfile.TemporaryDirectory() as tmp:
            npz_p = save_delta_bundle_npz(
                Path(tmp) / "bundle.npz",
                delta=np.zeros((1, 3), dtype=np.float32),
                z_ctrl=np.ones((2, 4), dtype=np.float32),
                target_genes=["P1"],
                eval_genes=["G1", "G2", "G3"],
                control_barcodes=["BC1", "BC2"],
                meta={"line": "test", "mode": "rhj_6983_axis"},
            )
            z = np.load(npz_p, allow_pickle=True)
            self.assertEqual(z["delta"].shape, (1, 3))
            self.assertEqual(z["z_ctrl"].shape, (2, 4))
            self.assertEqual(json.loads(str(z["meta_json"][0]))["mode"], "rhj_6983_axis")


if __name__ == "__main__":
    unittest.main(verbosity=2)
