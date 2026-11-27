#!/usr/bin/env python3
"""
Unit & End-to-End Test Suite for Clean-Room Round 16 Addendum B Audit
(`round16_addB/audit/test_audit_addb.py`, `A-T3`).

Covers:
  1. Freeze & GO runtime gate enforcement (`AddBBlindingGateError` + CLI `--phase real_addb` refusal)
  2. Targeted barcode-exact raw count extraction across Nadig legacy h5ad, 10x `.h5` (`dld1`/`cd4t_rest`/
     `gse295214_h1`), `a549_ko` `.mtx.gz` excluding `CAS9`, and `feng_ipsc` wide CSV
  3. Leading duplicate-Ensembl deduplication, full-transcriptome normalization, and `pairwise_common_axis`
  4. `H-I1_Occam`, `H-I2_Occam`, `H-I3a`, `H-I3b` per-source + `per_stratum_mean`, `P4` verdict-change
     detection, and `§4.3` Advance Outcome Matrix Cases 1–6
  5. Full standalone CLI `--phase synthetic_e2e` schema verification (`allow_pickle=False`, all 15
     keys, `148` `preds_Occam_*.npz` files, `reused_primary_rhj_6983`, Tier-2 `RBM14`, `gse295214_h1`
     `extra_fct_sensitivity`) AND independent standalone `--recompute-from-arrays` (`<= 1e-12`).
"""

from __future__ import annotations

import csv
import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

AUDIT_ADDB_DIR = Path(__file__).resolve().parent
ADDB_ROOT = AUDIT_ADDB_DIR.parent
REPO_ROOT = ADDB_ROOT.parent

if str(AUDIT_ADDB_DIR) not in sys.path:
    sys.path.insert(0, str(AUDIT_ADDB_DIR))

from audit_addb_controls_and_metrics import (  # noqa: E402
    ADDB_ALL_CONTROL_POOLS,
    ADDB_ALL_TARGET_LINES,
    ADDB_SEEN_SOURCES,
    ADDB_UNSEEN_SOURCES,
    EXPECTED_ADDB_V2_SHA256,
    AddBBlindingGateError,
    CommittedR16ReferenceBundle,
    _extract_10x_by_barcodes,
    _extract_a549_multi_pools_by_barcodes,
    _extract_feng_multi_pools_by_barcodes,
    _extract_nadig_by_barcodes,
    assert_freeze_gates_addb,
    classify_advance_outcome_matrix,
    compute_hi3a_and_hi3b_line,
    evaluate_addb_line_complete,
    sha256_file,
)
from audit_addb_occam_adapter import (  # noqa: E402
    compute_pairwise_common_rhj_mask,
    normalize_counts_to_rhj_44728,
    sha256_array_bytes,
    sha256_str_lines,
)
from run_audit_addb import write_npz_strict_addb  # noqa: E402


class TestAddBFreezeGatesAndExtractors(unittest.TestCase):
    def test_freeze_gate_requires_nonempty_freeze_and_go_files(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            (root / "PREREG_I_ADDENDUM_B_LOCK.md").write_text("LOCK\n", encoding="utf-8")
            with self.assertRaises(AddBBlindingGateError):
                assert_freeze_gates_addb(root)

            (root / "R16_ADDB_CODE_FREEZE_AT3.txt").write_text("FREEZE_OK\n", encoding="utf-8")
            with self.assertRaises(AddBBlindingGateError):
                assert_freeze_gates_addb(root)

            (root / "R16_ADDB_GO_AT3.txt").write_text("   \n", encoding="utf-8")
            with self.assertRaises(AddBBlindingGateError):
                assert_freeze_gates_addb(root)

            (root / "R16_ADDB_GO_AT3.txt").write_text("GO_OK\n", encoding="utf-8")
            res = assert_freeze_gates_addb(root)
            self.assertEqual(res["r16_addb_code_freeze_at3"], "FREEZE_OK")
            self.assertEqual(res["r16_addb_go_at3"], "GO_OK")

    def test_freeze_gate_exercises_real_addb_draft_v2_hash_and_d2_fixes(self) -> None:
        from audit_addb_occam_adapter import DEFAULT_GENE_AXIS_CSV, load_rhj_gene_axis
        from audit_addb_controls_and_metrics import load_locked_control_barcodes_addb

        # Exercise hash verification against real round16_addB/ files
        res = assert_freeze_gates_addb(ADDB_ROOT)
        actual_v2_sha = sha256_file(ADDB_ROOT / "PREREG_I_ADDENDUM_B_DRAFT_v2.md")
        self.assertEqual(res["addb_v2_sha256"], EXPECTED_ADDB_V2_SHA256)
        self.assertEqual(res["addb_v2_sha256"], actual_v2_sha)

        # Negative check: mismatched draft v2 file inside tempdir must raise AddBBlindingGateError
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            (root / "PREREG_I_ADDENDUM_B_LOCK.md").write_text("LOCK\n", encoding="utf-8")
            (root / "R16_ADDB_CODE_FREEZE_AT3.txt").write_text("FREEZE_OK\n", encoding="utf-8")
            (root / "R16_ADDB_GO_AT3.txt").write_text("GO_OK\n", encoding="utf-8")
            (root / "PREREG_I_ADDENDUM_B_DRAFT_v2.md").write_text("CORRUPT_V2\n", encoding="utf-8")
            with self.assertRaises(AddBBlindingGateError):
                assert_freeze_gates_addb(root)

        # DEVIATION_R16_ADDB_D2 item 2: load_rhj_gene_axis 2-tuple signature & gse295214_h1 n_ctrl_total=58153
        rhj_ens, rhj_syms = load_rhj_gene_axis(DEFAULT_GENE_AXIS_CSV)
        self.assertEqual(len(rhj_ens), 6983)
        self.assertEqual(len(rhj_syms), 6983)
        cb = load_locked_control_barcodes_addb(["gse295214_h1"])
        self.assertEqual(int(cb["pools"]["gse295214_h1"]["n_ctrl_total"]), 58153)

    def test_cli_real_addb_refuses_without_freeze_and_go_gates(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            empty_root = Path(td) / "empty_addb"
            empty_root.mkdir()
            out_dir = Path(td) / "out"
            py_snippet = (
                f"import sys; sys.path.insert(0, {str(AUDIT_ADDB_DIR)!r}); "
                "from pathlib import Path; from run_audit_addb import run_real_addb; "
                f"run_real_addb(output_dir=Path({str(out_dir)!r}), addb_root=Path({str(empty_root)!r}))"
            )
            proc = subprocess.run([sys.executable, "-c", py_snippet], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("AddBBlindingGateError", proc.stderr)

    def test_nadig_legacy_csr_and_10x_extractors_exact_row_sums(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            # 1. Legacy Nadig .h5ad with categorical byte array & CSR X
            h5p = tmp / "nadig_test.h5ad"
            with h5py.File(str(h5p), "w") as f:
                obs = f.create_group("obs")
                cb = obs.create_group("cell_barcode")
                cb.attrs["encoding-type"] = "categorical"
                cb.create_dataset("categories", data=np.asarray([b"BC_C", b"BC_A", b"BC_B"]))
                cb.create_dataset("codes", data=np.asarray([1, 2, 0], dtype=np.int8))  # rows: BC_A, BC_B, BC_C
                var = f.create_group("var")
                var.create_dataset("gene_id", data=np.asarray([b"ENSG1.1", b"ENSG2.2", b"ENSG3.3"]))
                xg = f.create_group("X")
                # Row 0 (BC_A): [1, 0, 4], Row 1 (BC_B): [0, 5, 0], Row 2 (BC_C): [2, 3, 1]
                xg.create_dataset("indptr", data=np.asarray([0, 2, 3, 6], dtype=np.int64))
                xg.create_dataset("indices", data=np.asarray([0, 2, 1, 0, 1, 2], dtype=np.int64))
                xg.create_dataset("data", data=np.asarray([1, 4, 5, 2, 3, 1], dtype=np.float32))

            ext_n = _extract_nadig_by_barcodes(h5p, "nadig_hepg2", ["BC_C", "BC_A"])
            self.assertEqual(ext_n.barcodes, ["BC_C", "BC_A"])
            self.assertEqual(ext_n.raw_ensembl_ids, ["ENSG1", "ENSG2", "ENSG3"])
            np.testing.assert_allclose(ext_n.raw_counts, np.asarray([[2, 3, 1], [1, 0, 4]], dtype=np.float32))
            np.testing.assert_allclose(ext_n.row_sum_umis, np.asarray([6.0, 5.0]))

            # 2. 10x CD4T h5 excluding CRISPR Guide Capture from RNA sum & mapping CE0006864_D4_CD4i_R2L15
            cdir = tmp / "cd4t"
            cdir.mkdir()
            h10 = cdir / "GSM999_CD4i_R2L15_D4_filtered.h5"
            with h5py.File(str(h10), "w") as f:
                m = f.create_group("matrix")
                m.create_dataset("barcodes", data=np.asarray([b"AA-1", b"BB-1"]))
                fg = m.create_group("features")
                fg.create_dataset("id", data=np.asarray([b"ENSG10.1", b"ENSG20.2", b"GUIDE_1"]))
                fg.create_dataset(
                    "feature_type",
                    data=np.asarray([b"Gene Expression", b"Gene Expression", b"CRISPR Guide Capture"]),
                )
                # col 0 (AA-1): [7, 1, 999]; col 1 (BB-1): [0, 3, 500]
                m.create_dataset("indptr", data=np.asarray([0, 3, 5], dtype=np.int64))
                m.create_dataset("indices", data=np.asarray([0, 1, 2, 1, 2], dtype=np.int64))
                m.create_dataset("data", data=np.asarray([7, 1, 999, 3, 500], dtype=np.float32))

            ext_10x = _extract_10x_by_barcodes(cdir, "cd4t_rest", ["CE0006864_D4_CD4i_R2L15:BB-1", "CE0006864_D4_CD4i_R2L15:AA-1"])
            self.assertEqual(ext_10x.raw_ensembl_ids, ["ENSG10", "ENSG20"])
            np.testing.assert_allclose(ext_10x.raw_counts, np.asarray([[0, 3], [7, 1]], dtype=np.float32))
            np.testing.assert_allclose(ext_10x.row_sum_umis, np.asarray([3.0, 8.0]))

    def test_10x_permuted_channel_order_aligns_and_mismatch_raises(self) -> None:
        def _write_10x(path: Path, bcs: list[bytes], fids: list[bytes], ftypes: list[bytes], cols_by_cell: np.ndarray) -> None:
            indptr = [0]
            indices = []
            data = []
            for ci in range(cols_by_cell.shape[0]):
                nz = np.flatnonzero(cols_by_cell[ci])
                indices.extend(nz.tolist())
                data.extend(cols_by_cell[ci, nz].tolist())
                indptr.append(len(indices))
            with h5py.File(str(path), "w") as f:
                m = f.create_group("matrix")
                m.create_dataset("barcodes", data=np.asarray(bcs))
                fg = m.create_group("features")
                fg.create_dataset("id", data=np.asarray(fids))
                fg.create_dataset("feature_type", data=np.asarray(ftypes))
                m.create_dataset("indptr", data=np.asarray(indptr, dtype=np.int64))
                m.create_dataset("indices", data=np.asarray(indices, dtype=np.int64))
                m.create_dataset("data", data=np.asarray(data, dtype=np.float32))

        with tempfile.TemporaryDirectory() as td:
            cdir = Path(td)
            # Channel 1 (R2L01): canonical order ENSG10, ENSG20, ENSG30 (+ CRISPR guide)
            _write_10x(
                cdir / "GSM998_CD4i_R2L01_D4_filtered.h5",
                [b"AA-1", b"BB-1"],
                [b"ENSG10.1", b"ENSG20.2", b"ENSG30.1", b"GUIDE_1"],
                [b"Gene Expression", b"Gene Expression", b"Gene Expression", b"CRISPR Guide Capture"],
                np.asarray([[7, 1, 2, 999], [0, 3, 5, 500]], dtype=np.float32),
            )
            # Channel 2 (R2L02): same unique RNA set in permuted order [ENSG30, ENSG10, ENSG20]
            _write_10x(
                cdir / "GSM999_CD4i_R2L02_D4_filtered.h5",
                [b"CC-1"],
                [b"ENSG30.9", b"ENSG10.4", b"ENSG20.8", b"GUIDE_2"],
                [b"Gene Expression", b"Gene Expression", b"Gene Expression", b"CRISPR Guide Capture"],
                np.asarray([[11, 4, 6, 888]], dtype=np.float32),
            )
            # Channel 3 (R2L03): genuine gene-set mismatch (ENSG99 != ENSG30)
            _write_10x(
                cdir / "GSM997_CD4i_R2L03_D4_filtered.h5",
                [b"DD-1"],
                [b"ENSG10.1", b"ENSG20.2", b"ENSG99.1"],
                [b"Gene Expression", b"Gene Expression", b"Gene Expression"],
                np.asarray([[1, 2, 3]], dtype=np.float32),
            )
            # Channel 4 (R2L04): same length as Channel 1 but duplicate Ensembl ID
            _write_10x(
                cdir / "GSM996_CD4i_R2L04_D4_filtered.h5",
                [b"EE-1"],
                [b"ENSG10.1", b"ENSG10.2", b"ENSG20.2"],
                [b"Gene Expression", b"Gene Expression", b"Gene Expression"],
                np.asarray([[1, 2, 3]], dtype=np.float32),
            )

            ext = _extract_10x_by_barcodes(
                cdir,
                "cd4t_rest",
                [
                    "CE0006864_D4_CD4i_R2L02:CC-1",
                    "CE0006864_D4_CD4i_R2L01:AA-1",
                    "CE0006864_D4_CD4i_R2L01:BB-1",
                ],
            )
            self.assertEqual(ext.raw_ensembl_ids, ["ENSG10", "ENSG20", "ENSG30"])
            np.testing.assert_allclose(
                ext.raw_counts,
                np.asarray([[4, 6, 11], [7, 1, 2], [0, 3, 5]], dtype=np.float32),
            )
            np.testing.assert_allclose(ext.row_sum_umis, np.asarray([21.0, 10.0, 8.0]))

            with self.assertRaises(ValueError):
                _extract_10x_by_barcodes(
                    cdir,
                    "cd4t_rest",
                    ["CE0006864_D4_CD4i_R2L01:AA-1", "CE0006864_D4_CD4i_R2L03:DD-1"],
                )
            with self.assertRaises(ValueError):
                _extract_10x_by_barcodes(
                    cdir,
                    "cd4t_rest",
                    ["CE0006864_D4_CD4i_R2L01:AA-1", "CE0006864_D4_CD4i_R2L04:EE-1"],
                )

    def test_a549_excludes_cas9_and_feng_wide_csv(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            a549_dir = tmp / "a549"
            a549_dir.mkdir()
            with gzip.open(a549_dir / "GSE345058_cell_annotations.csv.gz", "wt", encoding="utf-8", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["cell_id", "batch"])
                w.writerow(["C1", "B1"])
                w.writerow(["C2", "B1"])
                w.writerow(["C3", "B2"])
            with gzip.open(a549_dir / "GSE345058_gene_annotations.csv.gz", "wt", encoding="utf-8", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["ensembl_id", "genome"])
                w.writerow(["ENSG1.1", "GRCh38"])
                w.writerow(["CAS9_FEAT", "CAS9"])
                w.writerow(["ENSG2.4", "GRCh38"])
            with gzip.open(a549_dir / "GSE345058_raw_counts.mtx.gz", "wt", encoding="ascii") as fh:
                fh.write("%%MatrixMarket matrix coordinate real general\n")
                fh.write("3 3 5\n")
                fh.write("1 1 5\n1 2 99\n2 3 4\n3 1 2\n3 3 6\n")

            ext = _extract_a549_multi_pools_by_barcodes(a549_dir, {"pri": ["C3"], "res": ["C1"]})
            self.assertEqual(ext["pri"].raw_ensembl_ids, ["ENSG1", "ENSG2"])
            np.testing.assert_allclose(ext["pri"].raw_counts, np.asarray([[2, 6]], dtype=np.float32))
            np.testing.assert_allclose(ext["pri"].row_sum_umis, np.asarray([8.0]))  # CAS9 col 2 excluded
            np.testing.assert_allclose(ext["res"].raw_counts, np.asarray([[5, 0]], dtype=np.float32))

            feng_dir = tmp / "feng"
            feng_dir.mkdir()
            with gzip.open(feng_dir / "TargetedScreen_RNA-UMI-Counts.csv.gz", "wt", encoding="utf-8") as fh:
                fh.write("feature,FB_2,FB_1\n")
                fh.write("ENSG100.2:GENE1,4,1\n")
                fh.write("ENSG200.1:GENE2,2,9\n")
            ext_f = _extract_feng_multi_pools_by_barcodes(feng_dir, {"pri": ["FB_1", "FB_2"]})
            self.assertEqual(ext_f["pri"].raw_ensembl_ids, ["ENSG100", "ENSG200"])
            np.testing.assert_allclose(ext_f["pri"].raw_counts, np.asarray([[1, 9], [4, 2]], dtype=np.float32))
            np.testing.assert_allclose(ext_f["pri"].row_sum_umis, np.asarray([10.0, 6.0]))

    def test_dedup_leading_ensembl_and_pairwise_common_mask(self) -> None:
        raw_ids = ["ENSG1", "ENSG2", "ENSG1", "ENSG3"]
        raw = np.asarray([[2.0, 4.0, 99.0, 4.0]], dtype=np.float32)  # full row sum = 109.0
        rhj = ["ENSG1", "ENSG2", "ENSG4"]
        enc = ["ENSG1", "ENSG2", "ENSG3", "ENSG4"]
        mask = compute_pairwise_common_rhj_mask(raw_ids, ["ENSG1", "ENSG3", "ENSG4"], rhj)
        np.testing.assert_array_equal(mask, np.asarray([True, False, False]))

        x_enc, x_rhj, aud = normalize_counts_to_rhj_44728(raw, raw_ids, rhj, enc, target_sum=2000.0, active_rhj_mask=mask)
        self.assertEqual(aud["n_rhj_present_in_raw"], 2)
        self.assertEqual(aud["n_rhj_absent_in_raw"], 1)
        self.assertEqual(aud["n_rhj_masked_out"], 1)
        self.assertEqual(aud["n_active_injected"], 1)
        self.assertAlmostEqual(float(x_rhj[0, 0]), float(np.log1p(2.0 / 109.0 * 2000.0)), places=5)
        self.assertEqual(float(x_enc[0, 1]), 0.0)


class TestAddBMetricsAndAdvanceMatrix(unittest.TestCase):
    def _make_ref(self, seed: int = 123, n_p: int = 60, n_g: int = 35) -> CommittedR16ReferenceBundle:
        rng = np.random.default_rng(seed)
        sig = rng.normal(0.0, 0.45, size=(n_p, n_g))
        d_obs = sig + rng.normal(0.0, 0.1, size=(n_p, n_g))
        ha = sig + rng.normal(0.0, 0.12, size=(n_p, n_g))
        hb = sig + rng.normal(0.0, 0.12, size=(n_p, n_g))
        hv = np.ones(n_p, dtype=bool)
        tm = rng.normal(0.0, 0.12, size=(n_p, n_g))
        gm = rng.normal(0.0, 0.05, size=n_g)
        fct = tm + 0.50 * sig
        se_m = tm + 0.30 * sig
        se_r = se_m + rng.normal(0.0, 0.01, size=(n_p, n_g))
        seen = list(ADDB_SEEN_SOURCES)
        unseen = list(ADDB_UNSEEN_SOURCES["cd4t_rest"])
        se_mis = {s: se_m.copy() for s in seen + unseen}
        return CommittedR16ReferenceBundle(
            target_line="cd4t_rest",
            role="group2_primary_crispri",
            eval_perts=[f"P{i:02d}" for i in range(n_p)],
            eval_genes=[f"G{i:02d}" for i in range(n_g)],
            target_gene_idx=np.full(n_p, -1, dtype=np.int64),
            delta_obs=d_obs,
            delta_half_a=ha,
            delta_half_b=hb,
            half_valid=hv,
            delta_transfer=tm,
            source_grand_mean=gm,
            delta_fct=fct,
            seen_sources=seen,
            unseen_sources=unseen,
            matched_barcodes=[f"m{i}" for i in range(16)],
            resample_barcodes=[f"r{i}" for i in range(16)],
            mismatched_barcodes={s: [f"{s}_{i}" for i in range(16)] for s in seen + unseen},
            se_delta_matched=se_m,
            se_delta_resample=se_r,
            se_delta_mismatched=se_mis,
            mu_target=rng.uniform(0.3, 2.0, size=n_g),
            mu_sources=rng.uniform(0.3, 2.0, size=(6, n_g)),
            e1_constmed=tm + 0.15 * sig,
            n_ctrl_total=1000,
            ctrl_sampling_fraction=0.512,
        )

    def test_hi3a_hi3b_p4_verdict_change_and_advance_cases(self) -> None:
        ref = self._make_ref()
        rng = np.random.default_rng(456)
        sig = ref.delta_obs - ref.delta_transfer
        # Case where Occam beats STATE_SE accurately and specifically under primary,
        # yet P4 common axis drops matched-vs-mis difference on one source -> verdict change!
        occam_m = ref.delta_transfer + 0.70 * sig
        occam_r = occam_m + rng.normal(0.0, 0.01, size=occam_m.shape)
        occam_mis = {s: ref.delta_transfer + rng.normal(0.0, 0.04, size=occam_m.shape) for s in ref.seen_sources + ref.unseen_sources}

        hi3 = compute_hi3a_and_hi3b_line(ref, occam_m, occam_r, occam_mis, n_boot=200, seed=20261601)
        self.assertEqual(hi3["h_i3a"]["classification"], "OCCAM_SUPERIOR_TO_STATE_SE")
        self.assertIn("occam_summary", hi3["h_i3a"])
        self.assertIn("state_se_summary", hi3["h_i3a"])
        for sk in ("seen_by_checkpoint", "unseen_by_checkpoint"):
            sm = hi3["h_i3b"]["per_stratum_mean"][sk]
            self.assertEqual(sm["delta_spec"]["classification"], "OCCAM_MORE_CONTROL_SPECIFIC")

        p4_pairs = {}
        for idx, s in enumerate(ref.seen_sources + ref.unseen_sources):
            if idx == 0:
                # Identical matched vs mis on common axis -> EXACT_INVARIANCE -> flips MATCHED_BETTER
                p4_pairs[s] = {"matched": occam_m, "mismatched": occam_m.copy(), "n_common": 6570, "n_union_absent": 413, "reused_primary_rhj_6983": False, "mask_sha256": "a" * 64}
            else:
                p4_pairs[s] = {"matched": occam_m, "mismatched": occam_mis[s], "n_common": 6570, "n_union_absent": 413, "reused_primary_rhj_6983": False, "mask_sha256": "b" * 64}

        res = evaluate_addb_line_complete(
            ref=ref, occam_matched=occam_m, occam_resample=occam_r, occam_mismatched=occam_mis,
            p4_common_pairs=p4_pairs, occam_baseline_real_matched=occam_m * 0.9,
            occam_baseline_real_resample=occam_r * 0.9,
            occam_baseline_real_mismatched={k: v * 0.9 for k, v in occam_mis.items()},
            occam_full44728_matched=occam_m, occam_full44728_resample=occam_r, occam_full44728_mismatched=occam_mis,
            n_boot=200, seed=20261601,
        )
        self.assertTrue(res["secondary"]["pairwise_common_axis"]["any_verdict_changed_vs_primary"])
        first_src = (ref.seen_sources + ref.unseen_sources)[0]
        self.assertTrue(res["secondary"]["pairwise_common_axis"]["per_source"][first_src]["verdict_changed_vs_primary"])

        # Advance Outcome Matrix Cases 1-6 checks
        spec = {"line_summary": "LINE_MATCHED_SPECIFIC", "counts": {"INCONCLUSIVE": 0}, "resample_reference": {"any_source_exceeds_resample": True}}
        eq = {"line_summary": "LINE_EQUIVALENT_NONZERO_SIGNAL", "counts": {"INCONCLUSIVE": 0}}
        mixed = {"line_summary": "MIXED_OR_INCONCLUSIVE", "counts": {"INCONCLUSIVE": 0}}
        mis_b = {"line_summary": "LINE_MISMATCHED_BETTER", "counts": {"INCONCLUSIVE": 0}}
        self.assertEqual(classify_advance_outcome_matrix(eq, eq, {}, {"classification": "INDISTINGUISHABLE"}, {})["matrix_case"], "Case_1_Cross_Architecture_Control_Underutilization")
        self.assertEqual(classify_advance_outcome_matrix(eq, spec, {}, {"classification": "INDISTINGUISHABLE"}, {})["matrix_case"], "Case_2_Occam_Control_Sensitive_SE_Not")
        self.assertEqual(classify_advance_outcome_matrix(spec, eq, {}, {"classification": "INDISTINGUISHABLE"}, {})["matrix_case"], "Case_3_STATE_SE_Control_Sensitive_Occam_Not")
        self.assertEqual(classify_advance_outcome_matrix(eq, eq, {}, {"classification": "OCCAM_SUPERIOR_TO_STATE_SE"}, {})["matrix_case"], "Case_4_Occam_Prior_Accuracy_Gain_Without_Control_Specificity")
        self.assertEqual(classify_advance_outcome_matrix(spec, spec, {}, {"classification": "INDISTINGUISHABLE"}, {})["matrix_case"], "Case_5_Cross_Architecture_Control_Utilization")
        self.assertEqual(classify_advance_outcome_matrix(mixed, spec, {}, {"classification": "INDISTINGUISHABLE"}, {})["matrix_case"], "Case_6_Power_Limited_Or_OOD")
        self.assertEqual(classify_advance_outcome_matrix(eq, mis_b, {}, {"classification": "INDISTINGUISHABLE"}, {})["matrix_case"], "Case_6_Power_Limited_Or_OOD")

    def test_cb12_common_axis_rollup_rules(self) -> None:
        from audit_addb_controls_and_metrics import _evaluate_p4_common_axis_block

        ref = self._make_ref(seed=123)
        sig = ref.delta_obs - ref.delta_transfer
        rng = np.random.default_rng(456)
        occam_m = ref.delta_transfer + 0.70 * sig
        occam_mis = {s: ref.delta_transfer + rng.normal(0.0, 0.4, size=occam_m.shape) for s in ref.seen_sources + ref.unseen_sources}

        prim_eval = evaluate_addb_line_complete(
            ref=ref, occam_matched=occam_m, occam_resample=occam_m.copy(), occam_mismatched=occam_mis,
            n_boot=150, seed=20261601,
        )
        prim_h1 = prim_eval["primary"]["H_I1_Occam"]
        prim_m_orth = prim_eval["primary"]["Occam_matched_orth"]
        self.assertEqual(prim_h1["line_summary"], "LINE_MATCHED_SPECIFIC")

        # Pair 0 ("nadig_hepg2"): masked matched arm fails validity (pred == t_mean -> zero residual
        # -> NaN Pearson / valid_frac=0.0), while primary rhj_6983_axis matched arm stays positive/valid;
        # Pair 1 reuses primary axis.
        srcs = list(ref.seen_sources) + list(ref.unseen_sources)
        p4_pairs = {}
        for idx, s in enumerate(srcs):
            if idx == 0:
                p4_pairs[s] = {"matched": ref.delta_transfer.copy(), "mismatched": occam_mis[s], "n_common": 6570, "n_union_absent": 413, "reused_primary_rhj_6983": False, "mask_sha256": "a" * 64}
            elif idx == 1:
                p4_pairs[s] = {"matched": occam_m, "mismatched": occam_mis[s], "n_common": 6983, "n_union_absent": 0, "reused_primary_rhj_6983": True, "mask_sha256": "f" * 64}
            else:
                p4_pairs[s] = {"matched": occam_m, "mismatched": occam_mis[s], "n_common": 6570, "n_union_absent": 413, "reused_primary_rhj_6983": False, "mask_sha256": "b" * 64}

        p4_res, flagged = _evaluate_p4_common_axis_block(
            ref=ref, p4_common_pairs=p4_pairs, primary_h_i1=prim_h1, n_boot=150, seed=20261601,
            primary_orth_matched=prim_m_orth,
        )
        self.assertTrue(flagged)
        self.assertTrue(p4_res["line_verdict_changed_vs_primary"])
        # CB1.2(a): failing masked matched arm forces INCONCLUSIVE for pair 0
        self.assertEqual(p4_res["per_source"][srcs[0]]["pair_classification"]["category"], "INCONCLUSIVE")
        self.assertFalse(p4_res["per_source"][srcs[0]]["pair_classification"]["masked_matched_arm_gate_passed"])
        # CB1.2(c): reused pair 1 takes primary category
        self.assertEqual(p4_res["per_source"][srcs[1]]["pair_classification"]["category"], prim_h1["per_source"][srcs[1]]["category"])
        # CB1.2(b): line-level matched_arm_valid/positive come from primary rhj_6983_axis, cascade yields MIXED_OR_INCONCLUSIVE
        h1c = p4_res["H_I1_Occam_common_axis"]
        self.assertTrue(h1c["matched_arm_valid"])
        self.assertTrue(h1c["matched_positive_ci975"])
        self.assertAlmostEqual(h1c["matched_rho_ab_orth"], prim_m_orth["rho_ab_orth"], places=12)
        self.assertEqual(h1c["line_summary"], "MIXED_OR_INCONCLUSIVE")

        # Also verify end-to-end wiring inside evaluate_addb_line_complete
        full_with_p4 = evaluate_addb_line_complete(
            ref=ref, occam_matched=occam_m, occam_resample=occam_m.copy(), occam_mismatched=occam_mis,
            p4_common_pairs=p4_pairs, n_boot=150, seed=20261601,
        )
        full_h1c = full_with_p4["secondary"]["pairwise_common_axis"]["H_I1_Occam_common_axis"]
        self.assertEqual(full_h1c["line_summary"], "MIXED_OR_INCONCLUSIVE")
        self.assertAlmostEqual(full_h1c["matched_rho_ab_orth"], prim_m_orth["rho_ab_orth"], places=12)

        # Case B (Amendment CB1.2a): null-signal line where masked matched arm passes validity
        # (valid_frac >= 0.95, rho_ab_orth > 0, finite bounds) but NOT positivity (ci975_lo <= 0).
        # Per-pair gate checks validity only, so masked pairs keep their TOST category and the
        # cascade yields LINE_EQUIVALENT_NULL_SIGNAL rather than INCONCLUSIVE.
        null_pred = ref.delta_transfer + rng.normal(0.0, 0.04, size=occam_m.shape)
        null_mis = {s: null_pred.copy() for s in srcs}
        null_p4 = {
            s: (
                {
                    "matched": null_pred.copy(),
                    "mismatched": null_mis[s],
                    "n_common": 6983,
                    "n_union_absent": 0,
                    "reused_primary_rhj_6983": True,
                    "mask_sha256": "f" * 64,
                }
                if idx == 1
                else {
                    "matched": null_pred.copy(),
                    "mismatched": null_mis[s],
                    "n_common": 6570,
                    "n_union_absent": 413,
                    "reused_primary_rhj_6983": False,
                    "mask_sha256": "c" * 64,
                }
            )
            for idx, s in enumerate(srcs)
        }
        null_full = evaluate_addb_line_complete(
            ref=ref, occam_matched=null_pred, occam_resample=null_pred.copy(), occam_mismatched=null_mis,
            p4_common_pairs=null_p4, n_boot=150, seed=20261601,
        )
        null_h1 = null_full["primary"]["H_I1_Occam"]
        self.assertTrue(null_h1["matched_arm_valid"])
        self.assertFalse(null_h1["matched_positive_ci975"])
        self.assertEqual(null_h1["line_summary"], "LINE_EQUIVALENT_NULL_SIGNAL")
        null_p4_out = null_full["secondary"]["pairwise_common_axis"]
        for idx, s in enumerate(srcs):
            pcls = null_p4_out["per_source"][s]["pair_classification"]
            self.assertEqual(pcls["category"], "EXACT_INVARIANCE")
            self.assertTrue(pcls["masked_matched_arm_gate_passed"])
        null_h1c = null_p4_out["H_I1_Occam_common_axis"]
        self.assertTrue(null_h1c["matched_arm_valid"])
        self.assertFalse(null_h1c["matched_positive_ci975"])
        self.assertEqual(null_h1c["line_summary"], "LINE_EQUIVALENT_NULL_SIGNAL")
        self.assertFalse(null_p4_out["any_verdict_changed_vs_primary"])


class TestAddBSyntheticE2ECLISchemaAndRecompute(unittest.TestCase):
    def test_synthetic_e2e_cli_schema_allow_pickle_false_and_standalone_recompute(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            syn_out = Path(td) / "syn_e2e"
            cmd_e2e = [
                sys.executable,
                str(AUDIT_ADDB_DIR / "run_audit_addb.py"),
                "--phase",
                "synthetic_e2e",
                "--output-dir",
                str(syn_out),
                "--n-boot",
                "120",
            ]
            p1 = subprocess.run(cmd_e2e, capture_output=True, text=True)
            self.assertEqual(p1.returncode, 0, msg=f"synthetic_e2e failed:\n{p1.stdout}\n{p1.stderr}")

            # 1. Verify top-level files & control_barcodes_addb.json
            cb = json.loads((syn_out / "control_barcodes_addb.json").read_text(encoding="utf-8"))
            self.assertEqual(set(cb["pools"].keys()), set(ADDB_ALL_CONTROL_POOLS))
            for pname, pentry in cb["pools"].items():
                self.assertEqual(pentry["barcodes_sha256"], sha256_str_lines(pentry["barcodes"]))
                if pname in ADDB_ALL_TARGET_LINES:
                    self.assertTrue(pentry["within_line_resample"]["eligible"])
                    self.assertEqual(pentry["within_line_resample"]["overlap_with_primary"], 0)
                else:
                    self.assertFalse(pentry["within_line_resample"]["eligible"])

            # 2. Verify z_ctrl manifest, keys, allow_pickle=False, & no redundant zero-masked latents
            z_man = json.loads((syn_out / "z_ctrl" / "z_ctrl_manifest_sha256.json").read_text(encoding="utf-8"))
            self.assertEqual(len(list((syn_out / "z_ctrl").glob("z_ctrl_*.npz"))), len(z_man))
            for fn in z_man:
                with np.load(syn_out / "z_ctrl" / fn, allow_pickle=False) as zz:
                    self.assertEqual(zz["Z_ctrl"].dtype, np.float32)
                    self.assertEqual(zz["Z_ctrl"].shape[1], 512)
                    self.assertEqual(
                        sha256_array_bytes(np.ascontiguousarray(zz["Z_ctrl"], dtype=np.float32)),
                        str(zz["z_ctrl_sha256"]),
                    )
                    if "pairwise_common_axis__" in fn:
                        aud = json.loads(str(zz["audit_json"]))
                        self.assertGreater(aud["n_present_but_zeroed_by_mask"], 0)

            # 3. Verify all 5 array dirs, file counts (148 preds_Occam_*.npz total + 1 SENS_FCT), & 15 keys
            total_preds = 0
            for tline in ADDB_ALL_TARGET_LINES:
                tdir = syn_out / "arrays" / tline
                aman = json.loads((tdir / "array_manifest_sha256.json").read_text(encoding="utf-8"))
                self.assertIn("observed.npz", aman)
                if tline == "gse295214_h1":
                    self.assertIn("preds_SENS_FCT_full_r9_6_incl_vcc2025_h1.npz", aman)
                else:
                    self.assertNotIn("preds_SENS_FCT_full_r9_6_incl_vcc2025_h1.npz", aman)
                occam_npzs = sorted(tdir.glob("preds_Occam_*.npz"))
                total_preds += len(occam_npzs)
                for pfp in occam_npzs:
                    with np.load(pfp, allow_pickle=False) as pz:
                        self.assertEqual(len(pz.files), 15)
                        self.assertEqual(pz["delta_pred"].dtype, np.float64)
                        self.assertEqual(
                            sha256_array_bytes(np.ascontiguousarray(pz["delta_pred"], dtype=np.float64)),
                            str(pz["delta_sha256"]),
                        )
            self.assertEqual(total_preds, 148)

            # 4. Verify summary_addb.json strict schema structure & Tier-2 / Group-3 sensitivities
            summ = json.loads((syn_out / "summary_addb.json").read_text(encoding="utf-8"))
            self.assertEqual(
                set(summ["group_adjudications"].keys()),
                {"group2_crisprko", "group2_crispri", "group3_crispri", "sensitivity_only"},
            )
            self.assertIsNotNone(summ["lines"]["a549_ko"]["secondary"]["tier1_trained_only_sensitivity"])
            self.assertIsNone(summ["lines"]["dld1"]["secondary"]["tier1_trained_only_sensitivity"])
            self.assertIsNotNone(summ["lines"]["gse295214_h1"]["secondary"]["extra_fct_sensitivity"])
            self.assertIsNone(summ["lines"]["a549_ko"]["secondary"]["extra_fct_sensitivity"])
            self.assertTrue(
                summ["lines"]["dld1"]["secondary"]["pairwise_common_axis"]["per_source"]["nadig_hepg2"]["reused_primary_rhj_6983"]
            )
            self.assertFalse(
                summ["lines"]["dld1"]["secondary"]["pairwise_common_axis"]["per_source"]["a549_ko"]["reused_primary_rhj_6983"]
            )
            self.assertIsNotNone(
                summ["lines"]["dld1"]["secondary"]["pairwise_common_axis"]["H_I1_Occam_common_axis"]["line_summary"]
            )

            # 5. Standalone fresh-process CLI --recompute-from-arrays verification
            cmd_rec = [
                sys.executable,
                str(AUDIT_ADDB_DIR / "run_audit_addb.py"),
                "--recompute-from-arrays",
                str(syn_out),
                "--se-audit-dir",
                str(syn_out / "_mock_se_audit"),
                "--main-r16-dir",
                str(syn_out / "_mock_main_r16"),
                "--main-sens-dir",
                str(syn_out / "_mock_main_sens"),
                "--tol",
                "1e-12",
            ]
            p2 = subprocess.run(cmd_rec, capture_output=True, text=True)
            self.assertEqual(p2.returncode, 0, msg=f"recompute-from-arrays failed:\n{p2.stdout}\n{p2.stderr}")
            rec_out = json.loads(p2.stdout)
            self.assertTrue(rec_out["passed"])
            self.assertLessEqual(rec_out["max_abs_diff"], 1e-12)


if __name__ == "__main__":
    unittest.main(verbosity=2)
