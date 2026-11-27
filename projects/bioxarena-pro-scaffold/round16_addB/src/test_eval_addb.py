#!/usr/bin/env python3
"""
CPU-only, zero-real-cell unit test suite for Round 16 Addendum B (`Occam`)
evaluation engine (`round16_addB/src/run_eval_addb.py`).

Tests:
  1. Hard execution gate refusal (`R16_ADDB_CODE_FREEZE.txt` + `R16_ADDB_GO.txt` + lock + approval)
     and strict refusal of `device != "cpu"`.
  2. Strict refusal of Group 1 training lines (`nadig_hepg2`, `nadig_jurkat`) as targets and
     exact verification of `§2.3` strata for all 5 evaluated lines.
  3. `H-I3a` matched-control accuracy contrast (`OCCAM_SUPERIOR_TO_STATE_SE`,
     `STATE_SE_SUPERIOR_TO_OCCAM`, `INDISTINGUISHABLE`, and exact zero-tie guard).
  4. `H-I3b` per-source & per-stratum-mean control-specificity contrast + resample-excess.
  5. Pre-locked `§4.3` Advance Outcome Interpretation Matrix (`Case_1`–`Case_6`).
  6. `P4` `pairwise_common_axis` verdict-change flagging & full-6983 primary reuse.
  7. Strictly separated group adjudications (`group2_crisprko`, `group2_crispri`,
     `group3_crispri`, `sensitivity_only`) with zero pooling.
  8. Synthetic `gse295214_h1` multi-channel 10x `.h5` prefix map + permuted-gene alignment.
  9. Full end-to-end CPU synthetic smoke (`run_synthetic_smoke_addb`, zero real cells),
     `ADDB_ARRAY_SCHEMA.md` file completeness, manifest SHA-256 verification, and
     Stage 4 from-disk recomputation parity (`max_abs_diff <= 1e-12`).
"""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

ADDB_SRC_DIR = Path(__file__).resolve().parent
if str(ADDB_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(ADDB_SRC_DIR))

from export_addb import verify_exported_manifest  # noqa: E402
from raw_controls_addb import load_10x_multi_h5_controls_by_barcodes  # noqa: E402
from run_eval_addb import (  # noqa: E402
    ADDB_LINE_SPECS,
    ALL_EVALUATED_LINES_ADDB,
    EXCLUDED_GROUP1_TRAINING_LINES,
    AddBExecutionGateError,
    _adjudicate_common_axis_h_i1_cb1_2,
    assert_addb_execution_gates,
    build_group_adjudications_addb,
    build_gse295214_prefix_map,
    classify_advance_interpretation_case_addb,
    compute_h_i3a_contrast_addb,
    compute_h_i3b_contrast_addb,
    evaluate_addb_line_from_arrays,
    run_real_evaluation_addb,
    run_synthetic_smoke_addb,
    validate_target_lines_addb,
)


def _make_synthetic_observed_block(n_p: int = 44, n_g: int = 90, seed: int = 777):
    rng = np.random.default_rng(seed)
    tm = rng.normal(0.0, 0.3, size=(n_p, n_g))
    b_s = rng.normal(0.0, 0.2, size=n_g)
    true_spec = rng.normal(0.0, 0.55, size=(n_p, n_g))
    ha = tm + true_spec + rng.normal(0.0, 0.20, size=(n_p, n_g))
    hb = tm + true_spec + rng.normal(0.0, 0.20, size=(n_p, n_g))
    obs = 0.5 * (ha + hb)
    hv = np.ones(n_p, dtype=bool)
    t_idx = np.full(n_p, -1, dtype=np.int64)
    t_idx[:5] = np.arange(5, dtype=np.int64)
    fct = tm + 0.35 * true_spec + rng.normal(0.0, 0.25, size=(n_p, n_g))
    return {
        "delta_obs": obs,
        "delta_half_a": ha,
        "delta_half_b": hb,
        "half_valid": hv,
        "delta_transfer": tm,
        "source_grand_mean": b_s,
        "fct_delta_pred": fct,
        "target_indices": t_idx,
        "target_names": [f"T{i:03d}" for i in range(n_p)],
        "gene_names": [f"G{j:04d}" for j in range(n_g)],
        "true_spec": true_spec,
        "rng": rng,
    }


class TestExecutionGatesAndLineSpecs(unittest.TestCase):
    def test_01_execution_gates_refuse_missing_or_empty(self) -> None:
        with tempfile.TemporaryDirectory(prefix="addb_gate_") as td:
            p = Path(td)
            cf = p / "R16_ADDB_CODE_FREEZE.txt"
            go = p / "R16_ADDB_GO.txt"
            lk = p / "PREREG_I_ADDENDUM_B_LOCK.md"
            ha = p / "R16_ADDB_HUMAN_APPROVAL.txt"
            lk.write_text("locked\n", encoding="utf-8")
            ha.write_text("approved\n", encoding="utf-8")

            # Missing both freeze and GO -> refuse
            with self.assertRaises(AddBExecutionGateError):
                assert_addb_execution_gates(code_freeze_path=cf, go_path=go, lock_path=lk, human_approval_path=ha)

            # Empty freeze -> refuse
            cf.write_text("", encoding="utf-8")
            go.write_text("GO\n", encoding="utf-8")
            with self.assertRaises(AddBExecutionGateError):
                assert_addb_execution_gates(code_freeze_path=cf, go_path=go, lock_path=lk, human_approval_path=ha)

            # Missing GO -> refuse
            cf.write_text("sha256 file\n", encoding="utf-8")
            go.unlink()
            with self.assertRaises(AddBExecutionGateError):
                assert_addb_execution_gates(code_freeze_path=cf, go_path=go, lock_path=lk, human_approval_path=ha)

            # All present and non-empty -> gate passes, but device="mps" STILL refused!
            go.write_text("GO\n", encoding="utf-8")
            gates = assert_addb_execution_gates(code_freeze_path=cf, go_path=go, lock_path=lk, human_approval_path=ha)
            self.assertEqual(len(gates), 4)
            with self.assertRaises(ValueError):
                run_real_evaluation_addb(
                    ["a549_ko"],
                    code_freeze_path=cf,
                    go_path=go,
                    lock_path=lk,
                    human_approval_path=ha,
                    device="mps",
                )

    def test_02_group1_training_lines_strictly_refused(self) -> None:
        for g1 in EXCLUDED_GROUP1_TRAINING_LINES:
            with self.assertRaises(ValueError):
                validate_target_lines_addb([g1])
        with self.assertRaises(ValueError):
            validate_target_lines_addb(["unknown_line"])

        validated = validate_target_lines_addb(list(ALL_EVALUATED_LINES_ADDB))
        self.assertEqual(tuple(validated), ALL_EVALUATED_LINES_ADDB)
        for ln in ALL_EVALUATED_LINES_ADDB:
            spec = ADDB_LINE_SPECS[ln]
            seen = set(spec["seen_by_checkpoint"])
            unseen = set(spec["unseen_by_checkpoint"])
            self.assertEqual(seen, {"nadig_hepg2", "nadig_jurkat"})
            self.assertEqual(len(seen & unseen), 0)
            self.assertNotIn(ln, seen | unseen)


class TestHI3ContrastsAndP4(unittest.TestCase):
    def test_03_h_i3a_matched_accuracy_contrast_branches(self) -> None:
        blk = _make_synthetic_observed_block(n_p=48, n_g=90, seed=101)
        rng = blk["rng"]
        tm, ts = blk["delta_transfer"], blk["true_spec"]
        noise = rng.normal(0.0, 0.15, size=tm.shape)

        strong = tm + 0.65 * ts + noise
        weak = tm + 0.15 * ts + noise

        # Occam stronger than STATE_SE
        r1 = compute_h_i3a_contrast_addb(
            blk["delta_obs"], strong, weak, tm, blk["target_indices"],
            blk["source_grand_mean"], blk["delta_half_a"], blk["delta_half_b"], blk["half_valid"],
            n_boot=300, seed=20261601,
        )
        self.assertEqual(r1["r_orth_cal"]["classification"], "OCCAM_SUPERIOR_TO_STATE_SE")
        self.assertGreater(r1["r_orth_cal"]["ci975_lo"], 0.0)
        self.assertEqual(r1["r_dorth_cal"]["classification"], "OCCAM_SUPERIOR_TO_STATE_SE")

        # STATE_SE stronger than Occam
        r2 = compute_h_i3a_contrast_addb(
            blk["delta_obs"], weak, strong, tm, blk["target_indices"],
            blk["source_grand_mean"], blk["delta_half_a"], blk["delta_half_b"], blk["half_valid"],
            n_boot=300, seed=20261601,
        )
        self.assertEqual(r2["r_orth_cal"]["classification"], "STATE_SE_SUPERIOR_TO_OCCAM")
        self.assertLess(r2["r_orth_cal"]["ci975_hi"], 0.0)

        # Identical arrays -> exact zero tie -> INDISTINGUISHABLE
        r3 = compute_h_i3a_contrast_addb(
            blk["delta_obs"], strong, strong.copy(), tm, blk["target_indices"],
            blk["source_grand_mean"], blk["delta_half_a"], blk["delta_half_b"], blk["half_valid"],
            n_boot=300, seed=20261601,
        )
        self.assertEqual(r3["r_orth_cal"]["classification"], "INDISTINGUISHABLE")
        self.assertEqual(r3["r_orth_cal"]["estimate"], 0.0)

    def test_04_h_i3b_control_specificity_and_stratum_means(self) -> None:
        blk = _make_synthetic_observed_block(n_p=48, n_g=90, seed=202)
        rng = blk["rng"]
        tm, ts = blk["delta_transfer"], blk["true_spec"]
        noise = rng.normal(0.0, 0.15, size=tm.shape)

        seen = ["nadig_hepg2", "nadig_jurkat"]
        unseen = ["a549_ko", "cd4t_rest", "feng_ipsc"]
        all_s = seen + unseen

        # Occam has big matched-vs-mismatched drop; STATE_SE is invariant across controls
        occ_m = tm + 0.65 * ts + noise
        occ_r = occ_m + rng.normal(0.0, 0.005, size=tm.shape)
        occ_mis = {sk: tm + 0.10 * ts + noise + rng.normal(0.0, 0.05, size=tm.shape) for sk in all_s}

        se_m = tm + 0.35 * ts + noise
        se_r = se_m + rng.normal(0.0, 0.005, size=tm.shape)
        se_mis = {sk: se_m + rng.normal(0.0, 0.005, size=tm.shape) for sk in all_s}

        res = compute_h_i3b_contrast_addb(
            blk["delta_obs"], occ_m, occ_r, occ_mis, se_m, se_r, se_mis,
            seen, unseen, tm, blk["target_indices"], blk["source_grand_mean"],
            blk["delta_half_a"], blk["delta_half_b"], blk["half_valid"],
            n_boot=300, seed=20261601,
        )
        ro = res["r_orth_cal"]
        for sk in all_s:
            self.assertEqual(ro["per_source"][sk]["delta_spec"]["classification"], "OCCAM_MORE_CONTROL_SPECIFIC")
            self.assertEqual(ro["per_source"][sk]["delta_spec_excess"]["classification"], "OCCAM_MORE_CONTROL_SPECIFIC")
        self.assertEqual(ro["per_stratum_mean"]["seen_by_checkpoint"]["delta_spec"]["classification"], "OCCAM_MORE_CONTROL_SPECIFIC")
        self.assertEqual(ro["per_stratum_mean"]["unseen_by_checkpoint"]["delta_spec"]["classification"], "OCCAM_MORE_CONTROL_SPECIFIC")
        self.assertEqual(ro["per_stratum_mean"]["seen_by_checkpoint"]["delta_spec_excess"]["classification"], "OCCAM_MORE_CONTROL_SPECIFIC")

        # Swap Occam & STATE_SE -> STATE_SE_MORE_CONTROL_SPECIFIC everywhere
        res_sw = compute_h_i3b_contrast_addb(
            blk["delta_obs"], se_m, se_r, se_mis, occ_m, occ_r, occ_mis,
            seen, unseen, tm, blk["target_indices"], blk["source_grand_mean"],
            blk["delta_half_a"], blk["delta_half_b"], blk["half_valid"],
            n_boot=300, seed=20261601,
        )
        self.assertEqual(
            res_sw["r_orth_cal"]["per_stratum_mean"]["unseen_by_checkpoint"]["delta_spec"]["classification"],
            "STATE_SE_MORE_CONTROL_SPECIFIC",
        )

    def test_05_advance_interpretation_cases_1_through_6(self) -> None:
        self.assertTrue(
            classify_advance_interpretation_case_addb(
                h_i1_state_se="LINE_EQUIVALENT_NONZERO_SIGNAL",
                h_i1_occam="LINE_EQUIVALENT_NONZERO_SIGNAL",
                h_i3a_call="INDISTINGUISHABLE",
                occam_exceeds_resample=False,
                state_se_exceeds_resample=False,
            ).startswith("Case_1")
        )
        self.assertTrue(
            classify_advance_interpretation_case_addb(
                h_i1_state_se="LINE_EQUIVALENT_NONZERO_SIGNAL",
                h_i1_occam="LINE_MATCHED_SPECIFIC",
                h_i3a_call="OCCAM_SUPERIOR_TO_STATE_SE",
                occam_exceeds_resample=True,
                state_se_exceeds_resample=False,
            ).startswith("Case_2")
        )
        self.assertTrue(
            classify_advance_interpretation_case_addb(
                h_i1_state_se="LINE_MATCHED_SPECIFIC",
                h_i1_occam="LINE_EQUIVALENT_NONZERO_SIGNAL",
                h_i3a_call="STATE_SE_SUPERIOR_TO_OCCAM",
                occam_exceeds_resample=False,
                state_se_exceeds_resample=True,
            ).startswith("Case_3")
        )
        self.assertTrue(
            classify_advance_interpretation_case_addb(
                h_i1_state_se="LINE_EQUIVALENT_NONZERO_SIGNAL",
                h_i1_occam="LINE_EQUIVALENT_NONZERO_SIGNAL",
                h_i3a_call="OCCAM_SUPERIOR_TO_STATE_SE",
                occam_exceeds_resample=False,
                state_se_exceeds_resample=False,
            ).startswith("Case_4")
        )
        self.assertTrue(
            classify_advance_interpretation_case_addb(
                h_i1_state_se="LINE_MATCHED_SPECIFIC",
                h_i1_occam="LINE_MATCHED_SPECIFIC",
                h_i3a_call="INDISTINGUISHABLE",
                occam_exceeds_resample=True,
                state_se_exceeds_resample=True,
            ).startswith("Case_5")
        )
        self.assertTrue(
            classify_advance_interpretation_case_addb(
                h_i1_state_se="MIXED_OR_INCONCLUSIVE",
                h_i1_occam="LINE_MISMATCHED_BETTER",
                h_i3a_call="INDISTINGUISHABLE",
                occam_exceeds_resample=False,
                state_se_exceeds_resample=False,
            ).startswith("Case_6")
        )

    def test_06_p4_pairwise_common_axis_verdict_change_and_reuse(self) -> None:
        blk = _make_synthetic_observed_block(n_p=48, n_g=90, seed=303)
        rng = blk["rng"]
        tm, ts = blk["delta_transfer"], blk["true_spec"]
        noise = rng.normal(0.0, 0.15, size=tm.shape)
        spec = ADDB_LINE_SPECS["dld1"]
        all_s = sorted(set(spec["seen_by_checkpoint"]) | set(spec["unseen_by_checkpoint"]))

        # On primary rhj_6983_axis, cd4t_rest is much worse (MATCHED_BETTER), others equivalent
        occ_m = tm + 0.55 * ts + noise
        occ_r = occ_m.copy()
        occ_mis = {sk: occ_m.copy() for sk in all_s}
        occ_mis["cd4t_rest"] = tm + 0.10 * ts + noise + rng.normal(0.0, 0.05, size=tm.shape)

        # Under pairwise_common_axis, cd4t_rest becomes identical to matched -> flips MATCHED_BETTER -> EXACT_INVARIANCE!
        p4_in = {
            "a549_ko": {"matched": occ_m.copy(), "mismatched": occ_m.copy(), "n_active_rhj_genes": 6950, "axis_mask_sha256": "a" * 64},
            "cd4t_rest": {"matched": occ_m.copy(), "mismatched": occ_m.copy(), "n_active_rhj_genes": 6570, "axis_mask_sha256": "b" * 64},
        }
        se_arr = {"matched": occ_m.copy(), "resample": occ_m.copy(), "mismatched": {sk: occ_m.copy() for sk in all_s}}

        res = evaluate_addb_line_from_arrays(
            line_name="dld1",
            observed=blk,
            occam_arrays={
                "rhj_6983_axis_reconstructed": {"matched": occ_m, "resample": occ_r, "mismatched": occ_mis},
                "rhj_6983_axis_real": {"matched": occ_m, "resample": occ_r, "mismatched": occ_mis},
                "pairwise_common_axis_reconstructed": p4_in,
                "full_44728_api_reconstructed": {"matched": occ_m, "resample": occ_r, "mismatched": occ_mis},
            },
            state_se_arrays=se_arr,
            n_boot=250,
            seed=20261601,
        )
        p4 = res["secondary"]["pairwise_common_axis"]
        self.assertTrue(p4["any_verdict_changed_vs_primary"])
        self.assertTrue(p4["per_source"]["cd4t_rest"]["verdict_changed_vs_primary"])
        self.assertEqual(p4["per_source"]["cd4t_rest"]["n_common"], 6570)
        self.assertTrue(p4["per_source"]["nadig_hepg2"]["reused_primary_rhj_6983"])
        self.assertEqual(p4["per_source"]["nadig_hepg2"]["n_common"], 6983)

    def test_07_group_adjudications_strictly_separated_no_pooling(self) -> None:
        mock_lines = {
            "a549_ko": {"group_key": "group2_crisprko", "modality": "CRISPRko", "mandatory_disclosures": ["d1"], "primary": {"H_I1_Occam": {"line_summary": "LINE_MATCHED_SPECIFIC"}, "H_I2_Occam": {"classification": "FCT_SUPERIOR_TO_MODEL"}}},
            "dld1": {"group_key": "group2_crispri", "modality": "CRISPRi", "mandatory_disclosures": ["d2"], "primary": {"H_I1_Occam": {"line_summary": "LINE_EQUIVALENT_NONZERO_SIGNAL"}, "H_I2_Occam": {"classification": "INDISTINGUISHABLE"}}},
            "cd4t_rest": {"group_key": "group2_crispri", "modality": "CRISPRi", "mandatory_disclosures": ["d3"], "primary": {"H_I1_Occam": {"line_summary": "MIXED_OR_INCONCLUSIVE"}, "H_I2_Occam": {"classification": "MODEL_SUPERIOR_TO_FCT"}}},
            "gse295214_h1": {"group_key": "group3_crispri", "modality": "CRISPRi", "mandatory_disclosures": ["d4"], "primary": {"H_I1_Occam": {"line_summary": "MIXED_OR_INCONCLUSIVE"}, "H_I2_Occam": {"classification": "FCT_SUPERIOR_TO_MODEL"}}},
            "feng_ipsc": {"group_key": "sensitivity_only", "modality": "CRISPRi", "mandatory_disclosures": ["d5"], "primary": {"H_I1_Occam": {"line_summary": "LINE_EQUIVALENT_NULL_SIGNAL"}, "H_I2_Occam": {"classification": "FCT_SUPERIOR_TO_MODEL"}}},
        }
        gadj = build_group_adjudications_addb(mock_lines)
        self.assertEqual(set(gadj.keys()), {"group2_crisprko", "group2_crispri", "group3_crispri", "sensitivity_only"})
        self.assertEqual(gadj["group2_crisprko"]["lines"], ["a549_ko"])
        self.assertEqual(gadj["group2_crispri"]["lines"], ["cd4t_rest", "dld1"])
        self.assertNotIn("feng_ipsc", gadj["group2_crispri"]["lines"])
        self.assertEqual(gadj["sensitivity_only"]["lines"], ["feng_ipsc"])


class TestGse295214LoaderAndEndToEndSmoke(unittest.TestCase):
    def test_08_gse295214_prefix_map_and_h5_alignment(self) -> None:
        with tempfile.TemporaryDirectory(prefix="gse295214_syn_") as td:
            d = Path(td)
            genes_a = ["ENSG000001.1", "ENSG000002.2", "ENSG000003.3"]
            syms_a = ["G1", "G2", "G3"]
            # Second channel has reversed gene order (Lock D1 alignment check!)
            genes_b = list(reversed(genes_a))
            syms_b = list(reversed(syms_a))

            def _write_10x_h5(fp: Path, mat_cg: np.ndarray, bcs, gids, gsyms):
                # 10x HDF5 stores CSC (n_features, n_cells)
                csc = sp.csc_matrix(mat_cg.T.astype(np.int32))
                with h5py.File(fp, "w") as fh:
                    mg = fh.create_group("matrix")
                    mg.create_dataset("data", data=csc.data)
                    mg.create_dataset("indices", data=csc.indices)
                    mg.create_dataset("indptr", data=csc.indptr)
                    mg.create_dataset("shape", data=np.asarray(csc.shape, dtype=np.int64))
                    mg.create_dataset("barcodes", data=np.asarray([b.encode("utf-8") for b in bcs]))
                    fg = mg.create_group("features")
                    fg.create_dataset("id", data=np.asarray([g.encode("utf-8") for g in gids]))
                    fg.create_dataset("name", data=np.asarray([s.encode("utf-8") for s in gsyms]))
                    fg.create_dataset("feature_type", data=np.asarray([b"Gene Expression"] * len(gids)))

            m_a = np.array([[10, 20, 30], [1, 2, 3]], dtype=np.float32)
            m_b = np.array([[300, 200, 100], [6, 5, 4]], dtype=np.float32)  # reversed cols -> [100, 200, 300], [4, 5, 6]
            _write_10x_h5(d / "GSM8943713_lane1_filtered.h5", m_a, ["AA-1", "BB-1"], genes_a, syms_a)
            _write_10x_h5(d / "GSM8943714_lane2_filtered.h5", m_b, ["CC-1", "DD-1"], genes_b, syms_b)

            pmap, ref_pref = build_gse295214_prefix_map(d)
            self.assertEqual(ref_pref, "GSM8943713")
            req = ["GSM8943714:CC-1", "GSM8943713:AA-1", "GSM8943714:DD-1"]
            pool = load_10x_multi_h5_controls_by_barcodes(
                pmap, req, line_name="gse295214_h1", canonical_ref_prefix=ref_pref
            )
            self.assertEqual(pool.barcodes, req)
            self.assertEqual(pool.raw_ensembl_ids, ["ENSG000001", "ENSG000002", "ENSG000003"])
            np.testing.assert_allclose(pool.raw_counts[0], [100.0, 200.0, 300.0])
            np.testing.assert_allclose(pool.raw_counts[1], [10.0, 20.0, 30.0])
            np.testing.assert_allclose(pool.raw_counts[2], [4.0, 5.0, 6.0])

    def test_09_synthetic_end_to_end_and_stage4_from_disk_recompute(self) -> None:
        with tempfile.TemporaryDirectory(prefix="addb_e2e_smoke_") as td:
            out_dir = Path(td)
            res = run_synthetic_smoke_addb(out_dir, n_boot=300, seed=20261601)
            self.assertEqual(res["real_cells_touched"], 0)
            self.assertEqual(res["device"], "cpu")
            self.assertTrue(res["stage4_from_disk_recompute_passed"])
            self.assertLessEqual(res["stage4_from_disk_recompute_max_abs_diff"], 1e-12)

            # Verify all 5 lines evaluated and separated properly
            self.assertEqual(set(res["lines"].keys()), set(ALL_EVALUATED_LINES_ADDB))
            self.assertEqual(res["lines"]["a549_ko"]["primary"]["H_I1_Occam"]["line_summary"], "LINE_MATCHED_SPECIFIC")
            self.assertEqual(res["lines"]["dld1"]["primary"]["H_I1_Occam"]["line_summary"], "LINE_EQUIVALENT_NONZERO_SIGNAL")
            self.assertEqual(res["lines"]["feng_ipsc"]["primary"]["H_I1_Occam"]["line_summary"], "LINE_EQUIVALENT_NULL_SIGNAL")
            self.assertIsNotNone(res["lines"]["a549_ko"]["secondary"]["tier1_trained_only_sensitivity"])
            self.assertEqual(res["lines"]["a549_ko"]["secondary"]["tier1_trained_only_sensitivity"]["excluded_tier2_targets"], ["RBM14"])
            self.assertIsNotNone(res["lines"]["gse295214_h1"]["secondary"]["extra_fct_sensitivity"])

            # Verify manifests on disk (z_ctrl + every line's array_manifest_sha256.json)
            z_chk = verify_exported_manifest(out_dir / "z_ctrl", "z_ctrl_manifest_sha256.json")
            self.assertTrue(z_chk["ok"])
            self.assertEqual(z_chk["n_files"], 12)
            for ln in ALL_EVALUATED_LINES_ADDB:
                l_chk = verify_exported_manifest(out_dir / "arrays" / ln, "array_manifest_sha256.json")
                self.assertTrue(l_chk["ok"], f"Manifest check failed for {ln}: {l_chk}")
            self.assertTrue((out_dir / "control_barcodes_addb.json").is_file())
            self.assertTrue((out_dir / "summary_addb.json").is_file())
            doc = json.loads((out_path := out_dir / "summary_addb.json").read_text(encoding="utf-8"))
            self.assertEqual(doc["protocol"], "PREREG_I_ADDENDUM_B_v2")
            self.assertTrue(out_path.is_file())

    def test_10_cb1_2_common_axis_masked_matched_gate(self) -> None:
        blk = _make_synthetic_observed_block(n_p=48, n_g=90, seed=404)
        rng = blk["rng"]
        tm, ts = blk["delta_transfer"], blk["true_spec"]
        noise = rng.normal(0.0, 0.15, size=tm.shape)
        spec = ADDB_LINE_SPECS["dld1"]
        seen = list(spec["seen_by_checkpoint"])
        unseen = list(spec["unseen_by_checkpoint"])
        all_s = sorted(set(seen) | set(unseen))

        # Primary rhj_6983_axis: strong positive matched signal, MATCHED_BETTER across all 5 sources -> LINE_MATCHED_SPECIFIC
        occ_m = tm + 0.60 * ts + noise
        occ_r = occ_m + rng.normal(0.0, 0.005, size=tm.shape)
        occ_mis = {sk: tm + 0.10 * ts + noise + rng.normal(0.0, 0.04, size=tm.shape) for sk in all_s}
        se_arr = {"matched": occ_m.copy(), "resample": occ_r.copy(), "mismatched": {sk: occ_m.copy() for sk in all_s}}

        # Case A (Corrected CB1.2(a) — validity-only per-pair gate):
        #  - cd4t_rest: masked matched arm has constant zero residual (pm == tm -> valid_frac=0.0, NaN CI -> FAILS validity)
        #    -> forced INCONCLUSIVE (within_delta=False, cb1_2_masked_matched_gate_failed=True) -> flips line to MIXED_OR_INCONCLUSIVE.
        #  - a549_ko: masked matched arm is VALID (valid_frac >= 0.95, rho_ab_orth > 0, finite CI) but NON-POSITIVE
        #    (ci975_lo <= 0) while raw matched > mis -> keeps its TOST category (MATCHED_BETTER, NOT forced INCONCLUSIVE).
        p4_in = {
            "a549_ko": {
                "matched": tm + rng.normal(0.0, 0.25, size=tm.shape),  # valid, non-positive (ci975_lo <= 0)
                "mismatched": tm - 0.50 * ts + rng.normal(0.0, 0.15, size=tm.shape),  # negative -> raw matched > mis!
                "n_active_rhj_genes": 6950,
                "axis_mask_sha256": "a" * 64,
            },
            "cd4t_rest": {
                "matched": tm.copy(),  # constant zero residual -> undefined Pearson r -> valid_frac == 0.0!
                "mismatched": occ_mis["cd4t_rest"].copy(),
                "n_active_rhj_genes": 6570,
                "axis_mask_sha256": "b" * 64,
            },
        }
        res = evaluate_addb_line_from_arrays(
            line_name="dld1",
            observed=blk,
            occam_arrays={
                "rhj_6983_axis_reconstructed": {"matched": occ_m, "resample": occ_r, "mismatched": occ_mis},
                "rhj_6983_axis_real": {"matched": occ_m, "resample": occ_r, "mismatched": occ_mis},
                "pairwise_common_axis_reconstructed": p4_in,
                "full_44728_api_reconstructed": {"matched": occ_m, "resample": occ_r, "mismatched": occ_mis},
            },
            state_se_arrays=se_arr,
            n_boot=250,
            seed=20261601,
        )
        self.assertEqual(res["primary"]["H_I1_Occam"]["line_summary"], "LINE_MATCHED_SPECIFIC")
        p4 = res["secondary"]["pairwise_common_axis"]
        cd_entry = p4["per_source"]["cd4t_rest"]
        self.assertFalse(cd_entry["masked_matched_arm_valid"])
        self.assertTrue(cd_entry["pair_classification"]["cb1_2_masked_matched_gate_failed"])
        self.assertEqual(cd_entry["pair_classification"]["category"], "INCONCLUSIVE")
        self.assertFalse(cd_entry["pair_classification"]["within_delta"])
        self.assertTrue(cd_entry["verdict_changed_vs_primary"])

        # a549_ko passes masked validity despite non-positive masked matched arm -> keeps MATCHED_BETTER
        a549_entry = p4["per_source"]["a549_ko"]
        self.assertTrue(a549_entry["masked_matched_arm_valid"])
        self.assertFalse(a549_entry["masked_matched_arm_positive"])
        self.assertFalse(a549_entry["cb1_2_forced_inconclusive"])
        self.assertFalse(a549_entry["pair_classification"]["cb1_2_masked_matched_gate_failed"])
        self.assertEqual(a549_entry["pair_classification"]["category"], "MATCHED_BETTER")

        # Reused pairs retain primary MATCHED_BETTER; line flips due to invalid cd4t_rest masked matched arm
        self.assertEqual(p4["per_source"]["nadig_hepg2"]["pair_classification"]["category"], "MATCHED_BETTER")
        self.assertEqual(p4["H_I1_Occam_common_axis"]["per_source"]["cd4t_rest"]["category"], "INCONCLUSIVE")
        self.assertEqual(p4["H_I1_Occam_common_axis"]["line_summary"], "MIXED_OR_INCONCLUSIVE")
        self.assertTrue(p4["any_verdict_changed_vs_primary"])

        # Case B (Amendment CB1.2a null-signal line): Both reused pairs (n_union_absent == 0) AND
        # valid masked pairs (n_union_absent > 0) keep their TOST category (EXACT_INVARIANCE / EQUIVALENT),
        # so H_I1_Occam_common_axis preserves LINE_EQUIVALENT_NULL_SIGNAL via primary line-level gate CB1.2(b).
        null_m = tm + rng.normal(0.0, 0.30, size=tm.shape)
        null_mis = {sk: null_m.copy() for sk in all_s}
        p4_null = {
            "cd4t_rest": {
                "matched": null_m.copy(),
                "mismatched": null_m.copy(),
                "n_active_rhj_genes": 6570,
                "axis_mask_sha256": "c" * 64,
            },
        }
        res_null = evaluate_addb_line_from_arrays(
            line_name="dld1",
            observed=blk,
            occam_arrays={
                "rhj_6983_axis_reconstructed": {"matched": null_m, "resample": null_m.copy(), "mismatched": null_mis},
                "rhj_6983_axis_real": {"matched": null_m, "resample": null_m.copy(), "mismatched": null_mis},
                "pairwise_common_axis_reconstructed": p4_null,
                "full_44728_api_reconstructed": {"matched": null_m, "resample": null_m.copy(), "mismatched": null_mis},
            },
            state_se_arrays=se_arr,
            n_boot=250,
            seed=20261601,
        )
        self.assertEqual(res_null["primary"]["H_I1_Occam"]["line_summary"], "LINE_EQUIVALENT_NULL_SIGNAL")
        p4_n = res_null["secondary"]["pairwise_common_axis"]
        self.assertEqual(p4_n["per_source"]["nadig_hepg2"]["pair_classification"]["category"], "EXACT_INVARIANCE")
        self.assertFalse(p4_n["per_source"]["nadig_hepg2"]["cb1_2_forced_inconclusive"])
        self.assertTrue(p4_n["per_source"]["cd4t_rest"]["masked_matched_arm_valid"])
        self.assertFalse(p4_n["per_source"]["cd4t_rest"]["masked_matched_arm_positive"])
        self.assertEqual(p4_n["per_source"]["cd4t_rest"]["cb1_2_raw_pair_category"], "EXACT_INVARIANCE")
        self.assertEqual(p4_n["per_source"]["cd4t_rest"]["pair_classification"]["category"], "EXACT_INVARIANCE")
        self.assertFalse(p4_n["per_source"]["cd4t_rest"]["cb1_2_forced_inconclusive"])
        self.assertFalse(p4_n["per_source"]["cd4t_rest"]["pair_classification"]["cb1_2_masked_matched_gate_failed"])
        self.assertEqual(p4_n["H_I1_Occam_common_axis"]["line_summary"], "LINE_EQUIVALENT_NULL_SIGNAL")
        self.assertFalse(p4_n["any_verdict_changed_vs_primary"])

        # Case C (CB1.2(b)): Line-level matched-arm gate uses primary rhj_6983_axis matched arm.
        eq_calls = {sk: {"category": "EQUIVALENT", "within_delta": True} for sk in all_s}
        inv_prim = {"valid_frac": 0.80, "rho_ab_orth": 0.5, "ci975_lo": 0.2, "ci975_hi": 0.4}
        self.assertEqual(
            _adjudicate_common_axis_h_i1_cb1_2(inv_prim, eq_calls, seen, unseen)["line_summary"],
            "INCONCLUSIVE",
        )
        null_prim = {"valid_frac": 1.0, "rho_ab_orth": 0.5, "ci975_lo": -0.05, "ci975_hi": 0.05}
        self.assertEqual(
            _adjudicate_common_axis_h_i1_cb1_2(null_prim, eq_calls, seen, unseen)["line_summary"],
            "LINE_EQUIVALENT_NULL_SIGNAL",
        )
        pos_prim = {"valid_frac": 1.0, "rho_ab_orth": 0.5, "ci975_lo": 0.15, "ci975_hi": 0.35}
        mb_calls = {sk: {"category": "MATCHED_BETTER", "within_delta": False} for sk in all_s}
        self.assertEqual(
            _adjudicate_common_axis_h_i1_cb1_2(pos_prim, mb_calls, seen, unseen)["line_summary"],
            "LINE_MATCHED_SPECIFIC",
        )


if __name__ == "__main__":
    unittest.main()
