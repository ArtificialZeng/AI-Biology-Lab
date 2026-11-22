#!/usr/bin/env python3
"""
Unit tests for Round 11 counts-only verification, Q3 feasibility gate,
NTC-only vs safe-targeting/unassigned exclusions, TeloHAEC plain-target filtering,
Feng disk-sharded CSV->CSR streaming, and stratified seed=42 split fallback (A-T1).
Strictly tests synthetic fixtures and rules only; touches ZERO R11 expression data.
"""

from __future__ import annotations

import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

SCRIPTS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = Path(__file__).resolve().parents[3]
R10_SCRIPTS = PROJECT_ROOT / "round10" / "shared" / "scripts"
for p in (SCRIPTS_DIR, R10_SCRIPTS):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from inspect_hcasmc_feasibility import run_feasibility  # noqa: E402
from prepare_r10_lines import HGNCLookup  # noqa: E402
from prepare_r11_lines import (  # noqa: E402
    _stream_feng_csv_to_csr_h5,
    apply_stratified_seed42_split,
    is_eligible_telohaec_plain_target,
    make_csr_h5_subset_readers,
)
from verify_r11_lines import evaluate_confirmatory_lines_r11  # noqa: E402


class TestRound11Rules(unittest.TestCase):
    def test_confirmatory_gate_excludes_feng_and_handles_q3_and_p1(self) -> None:
        counts = {
            "telohaec": {"g_eval_canon": 5500, "n_shared_reproducible_eval_canon": 75},
            "hcasmc": {"g_eval_canon": 5800, "n_shared_reproducible_eval_canon": 60},
            "cd4t_stim48": {"g_eval_canon": 6100, "n_shared_reproducible_eval_canon": 280},
            "feng_ipsc": {"g_eval_canon": 6000, "n_shared_reproducible_eval_canon": 250},
        }
        passing, txt = evaluate_confirmatory_lines_r11(counts, {"feasible": True})
        self.assertEqual(passing, ["telohaec", "hcasmc", "cd4t_stim48"])
        self.assertIn("feng_ipsc: SECONDARY_ONLY", txt)
        self.assertNotIn("\nfeng_ipsc\n", txt)

        passing_q3, txt_q3 = evaluate_confirmatory_lines_r11(counts, {"feasible": False, "reason": "No guide rows"})
        self.assertEqual(passing_q3, ["telohaec", "cd4t_stim48"])
        self.assertIn("EXCLUDED_FEASIBILITY_Q3", txt_q3)
        self.assertIn("K_CONF=2", txt_q3)

        counts_low = {
            "telohaec": {"g_eval_canon": 5500, "n_shared_reproducible_eval_canon": 75},
            "cd4t_stim48": {"g_eval_canon": 6100, "n_shared_reproducible_eval_canon": 40},
            "feng_ipsc": {"g_eval_canon": 6000, "n_shared_reproducible_eval_canon": 250},
        }
        passing_p1, txt_p1 = evaluate_confirmatory_lines_r11(counts_low, {"feasible": False, "reason": "OOM"})
        self.assertEqual(passing_p1, ["telohaec"])
        self.assertIn("P1_GUARD", txt_p1)

    def test_hcasmc_q3_feasibility_and_safe_exclusion(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tdir = Path(td)
            gdf = pd.DataFrame([
                {"guide_id": "non-targeting_1", "spacer": "AAAA", "targeting": False, "type": "non-targeting", "intended_target_name": "non-targeting"},
                {"guide_id": "0Safe_safe_1", "spacer": "CCCC", "targeting": False, "type": "safe-targeting", "intended_target_name": "safe-targeting"},
                {"guide_id": "TP53_1", "spacer": "GGGG", "targeting": True, "type": "targeting", "intended_target_name": "ENSG00000141510"},
            ])
            gpath = tdir / "guides.tsv.gz"
            gdf.to_csv(gpath, sep="\t", index=False, compression="gzip")

            n_rna, n_cells = 3005, 60
            rng = np.random.default_rng(7)
            rna_mat = sp.random(n_rna, n_cells, density=0.3, format="csc", random_state=rng, data_rvs=lambda n: np.full(n, 2, dtype=np.int32))
            crp = np.zeros((3, n_cells), dtype=np.int32)
            crp[0, :15] = 15
            crp[1, 15:25] = 15
            crp[2, 25:55] = 15
            crp[0, 55:] = 10
            crp[2, 55:] = 10
            full_csc = sp.vstack([rna_mat, sp.csc_matrix(crp)], format="csc")

            h5p = tdir / "synth_hcasmc.h5"
            with h5py.File(h5p, "w") as f:
                m = f.create_group("matrix")
                m.create_dataset("shape", data=np.array(full_csc.shape, dtype=np.int64))
                m.create_dataset("data", data=full_csc.data.astype(np.uint32))
                m.create_dataset("indices", data=full_csc.indices.astype(np.uint32))
                m.create_dataset("indptr", data=full_csc.indptr.astype(np.uint32))
                bcs = [f"AAAC-BBBB-CCCC__{'3A' if i % 2 == 0 else '3B'}".encode() for i in range(n_cells)]
                m.create_dataset("barcodes", data=np.array(bcs))
                fg = m.create_group("features")
                fids = [f"ENSG{i:011d}".encode() for i in range(n_rna)] + [b"non-targeting_1", b"0Safe_safe_1", b"TP53_1"]
                fnms = [f"GENE{i}".encode() for i in range(n_rna)] + [b"non-targeting_1", b"0Safe_safe_1", b"TP53_1"]
                ftps = [b"Gene Expression"] * n_rna + [b"CRISPR Guide Capture"] * 3
                fg.create_dataset("id", data=np.array(fids))
                fg.create_dataset("name", data=np.array(fnms))
                fg.create_dataset("feature_type", data=np.array(ftps))

            res = run_feasibility(h5p, gpath, tdir / "feas.json", batch_cols=20, mt_max_pct=10.0, min_genes=100)
            self.assertTrue(res["feasible"])
            self.assertEqual(res["n_single_sgrna_ntc"], 15)
            self.assertEqual(res["n_single_sgrna_safe_excluded"], 10)
            self.assertEqual(res["n_single_sgrna_targeting"], 30)

            h5_bad = tdir / "bad_no_guides.h5"
            rna_only = sp.csc_matrix(rna_mat)
            with h5py.File(h5_bad, "w") as f:
                m = f.create_group("matrix")
                m.create_dataset("shape", data=np.array(rna_only.shape, dtype=np.int64))
                m.create_dataset("data", data=rna_only.data.astype(np.uint32))
                m.create_dataset("indices", data=rna_only.indices.astype(np.uint32))
                m.create_dataset("indptr", data=rna_only.indptr.astype(np.uint32))
                m.create_dataset("barcodes", data=np.array([f"BC{i}".encode() for i in range(n_cells)]))
                fg = m.create_group("features")
                fg.create_dataset("id", data=np.array([f"ENSG{i:011d}".encode() for i in range(n_rna)]))
                fg.create_dataset("name", data=np.array([f"GENE{i}".encode() for i in range(n_rna)]))
                fg.create_dataset("feature_type", data=np.array([b"Gene Expression"] * n_rna))

            res_bad = run_feasibility(h5_bad, gpath, tdir / "feas_bad.json", batch_cols=20, min_genes=100)
            self.assertFalse(res_bad["feasible"])

    def test_telohaec_plain_target_and_seed42_split_and_feng_stream(self) -> None:
        hgnc = HGNCLookup(
            approved={"TP53", "BABAM2"},
            ens_to_sym={"ENSG00000141510": "TP53", "ENSG00000158019": "BABAM2"},
            sym_to_ens={"TP53": "ENSG00000141510", "BABAM2": "ENSG00000158019"},
            prev_to_sym={"BRE": "BABAM2"},
            alias_to_sym={},
        )
        ok_tp53, ens, sym, st = is_eligible_telohaec_plain_target("TP53", hgnc)
        self.assertTrue(ok_tp53 and sym == "TP53" and ens == "ENSG00000141510")
        ok_bre, ens_b, sym_b, st_b = is_eligible_telohaec_plain_target("BRE", hgnc)
        self.assertTrue(ok_bre and sym_b == "BABAM2" and st_b == "prev_symbol")
        for bad in ("TP53-TSS2", "Enhancer-at-CAD-SNP-rs590121", "9p21-add-AL354709-1", "FAM118B-and-RPUSD4", "negative-control", "safe-targeting", "Sox17"):
            self.assertFalse(is_eligible_telohaec_plain_target(bad, hgnc)[0], f"Expected reject for {bad}")

        # Stratified seed=42 fallback
        nt = np.array([True] * 10 + [False] * 12)
        pks = [""] * 10 + ["G1"] * 6 + ["G2"] * 6
        h1 = apply_stratified_seed42_split(pks, nt, seed=42)
        h2 = apply_stratified_seed42_split(pks, nt, seed=42)
        self.assertTrue(np.array_equal(h1, h2))
        self.assertEqual(int((h1[nt] == 1).sum()), 5)
        self.assertEqual(int((h1[nt] == 2).sum()), 5)
        self.assertEqual(int((h1[10:16] == 1).sum()), 3)
        self.assertEqual(int((h1[16:] == 2).sum()), 3)

        # Feng transposed CSV stream with shard size = 2
        with tempfile.TemporaryDirectory() as td:
            tdir = Path(td)
            csv_p = tdir / "rna.csv.gz"
            with gzip.open(csv_p, "wt", encoding="utf-8") as f:
                f.write('"","C_KOLF","C1","C_UN","C2","C3"\n')
                f.write('"ENSG00000141510:TP53:Gene-Expression",9,3,4,0,5\n')
                f.write('"ENSG00000198888:MT-ND1:Gene-Expression",1,2,0,4,0\n')
                f.write('"GUIDE_X:CRISPR-Guide-Capture",50,10,0,0,0\n')
            h5_out = tdir / "feng_csr.h5"
            gids, gnames, umi, ng, pct_mt = _stream_feng_csv_to_csr_h5(csv_p, ["C1", "C2", "C3"], h5_out, cell_shard_size=2)
            self.assertEqual(gnames, ["TP53", "MT-ND1"])
            self.assertEqual(gids, ["ENSG00000141510", "ENSG00000198888"])
            np.testing.assert_allclose(umi, [5.0, 4.0, 5.0])
            np.testing.assert_array_equal(ng, [2, 1, 1])
            np.testing.assert_allclose(pct_mt, [40.0, 100.0, 0.0])
            rc, rr, rcc = make_csr_h5_subset_readers(h5_out, "X", np.array([0, 2]), 2)
            mat = rc(0, 2)
            np.testing.assert_allclose(mat, [[3.0, 2.0], [5.0, 0.0]])
            sub_col = rcc(0, 2, np.array([0]))
            np.testing.assert_allclose(sub_col, [[3.0], [5.0]])


if __name__ == "__main__":
    unittest.main()
