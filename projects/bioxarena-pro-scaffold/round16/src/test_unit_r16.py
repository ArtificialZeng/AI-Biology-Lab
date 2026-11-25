#!/usr/bin/env python3
"""
Unit & Invariance Test Suite for Round 16 (`PREREG_I v5`, incorporating `I14`, `M1–M6`,
`Q4`, Gate-2 locked constants, triple execution gate, and pre-download `GSE295214` reader).

Runs strictly on 100% synthetic fixtures (touches ZERO real perturbation datasets,
downloads ZERO weights, opens ZERO archived `adata_pred.h5ad` files).
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

SRC_DIR = Path(__file__).resolve().parent
REPO_ROOT = SRC_DIR.parents[1]
R12_SRC_DIR = REPO_ROOT / "round12" / "src"
for _p in (str(SRC_DIR), str(R12_SRC_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from adjudicate_r16 import (  # noqa: E402
    DEFAULT_TOST_DELTA,
    adjudicate_claim_i1_line,
    adjudicate_claim_i2_line,
    adjudicate_group_r16,
    classify_pair_tost_r16,
)
from metrics_r12 import compute_joint_sa8_split_minus_all_r12  # noqa: E402
from metrics_r16 import (  # noqa: E402
    DEFAULT_R16_SEED,
    compute_joint_cross_mask_r16,
    compute_orth_suite_r16,
    compute_paired_contrast_r16,
    compute_paired_model_vs_all_mismatched_r16,
)
from model_adapters_r16 import (  # noqa: E402
    DEFAULT_MODELS_DIR,
    GROUP2_UNSEEN_PRIMARY_LINES,
    GROUP3_PRIMARY_FCT_SOURCES,
    ArcStateAdapter,
    BlindingGateError,
    CellularAIOpenCodePublicRefitAdapter,
    InsufficientControlCellsError,
    TrainingLeakageError,
    assert_all_execution_gates_r16,
    assert_code_freeze_r16,
    assert_prereg_i_locked,
    assert_r16_go,
    assert_zero_training_leakage,
    build_mismatched_strata_r16,
    compute_embedding_separation_diagnostics,
    compute_pseudobulk_log1p_shift,
    preprocess_hvg_counts_r16,
    sample_control_indices,
    sample_stratified_control_sets_r16,
)
from run_eval_r16 import (  # noqa: E402
    build_gse295214_parts,
    extract_gse295214_ntc_raw_counts,
    parse_gse295214_guide_r16,
    run_real_evaluation_r16,
    run_synthetic_smoke_r16,
)
from theory_sim_r16 import (  # noqa: E402
    illustrate_mechanism_m2_synthetic,
    verify_property_p1_synthetic,
)

_HAS_ARC_STATE = (
    importlib.util.find_spec("torch") is not None
    and importlib.util.find_spec("state") is not None
)


class _StubHGNC:
    """Minimal synthetic HGNC lookup for zero-disk-dependency reader tests."""

    def __init__(self, approved: set[str]):
        self.approved = set(approved)
        self.sym_to_ens = {g: f"ENSG000000{i:05d}" for i, g in enumerate(sorted(approved), start=1)}

    def resolve(self, raw_sym: str, ens_id: str = "") -> tuple[str, str]:
        s = str(raw_sym).strip()
        if s in self.approved:
            return s, "approved"
        return s, "unmapped"


def _write_synthetic_10x_h5(
    path: Path,
    rna_counts: np.ndarray,
    guide_counts: np.ndarray,
    gene_names: list[str],
    guide_names: list[str],
    barcodes: list[str],
) -> None:
    """Write a standard 10x CellRanger `/matrix` CSC HDF5 fixture."""
    full = np.vstack([rna_counts, guide_counts])  # (n_features, n_cells)
    csc = sp.csc_matrix(full)
    n_g = len(gene_names)
    n_crp = len(guide_names)
    f_ids = [f"ENSG_SYN_{i:04d}" for i in range(n_g)] + list(guide_names)
    f_names = list(gene_names) + list(guide_names)
    f_types = ["Gene Expression"] * n_g + ["CRISPR Guide Capture"] * n_crp

    with h5py.File(path, "w") as f:
        m = f.create_group("matrix")
        m.create_dataset("shape", data=np.array(csc.shape, dtype=np.int64))
        m.create_dataset("data", data=csc.data.astype(np.float32))
        m.create_dataset("indices", data=csc.indices.astype(np.int32))
        m.create_dataset("indptr", data=csc.indptr.astype(np.int64))
        m.create_dataset("barcodes", data=np.array(barcodes, dtype=h5py.string_dtype("utf-8")))
        fg = m.create_group("features")
        fg.create_dataset("id", data=np.array(f_ids, dtype=h5py.string_dtype("utf-8")))
        fg.create_dataset("name", data=np.array(f_names, dtype=h5py.string_dtype("utf-8")))
        fg.create_dataset("feature_type", data=np.array(f_types, dtype=h5py.string_dtype("utf-8")))


class TestRound16ScaffoldV5(unittest.TestCase):
    def test_1_primary_r_orth_exact_parity_vs_round12(self) -> None:
        """Primary `r_orth_cal` (`source_grand_mean=None`) must match frozen Round 12 <= 1e-12."""
        rng = np.random.default_rng(101)
        P, G = 45, 300
        tm = rng.normal(0.0, 0.2, size=(P, G))
        true_d = tm + rng.normal(0.0, 0.25, size=(P, G))
        ha = true_d + rng.normal(0.0, 0.05, size=(P, G))
        hb = true_d + rng.normal(0.0, 0.05, size=(P, G))
        obs = 0.5 * (ha + hb)
        pred = tm + 0.4 * (true_d - tm) + rng.normal(0.0, 0.1, size=(P, G))
        t_idx = np.arange(P, dtype=np.int64)
        hv = np.ones(P, dtype=bool)

        r12 = compute_joint_sa8_split_minus_all_r12(
            obs, pred, tm, t_idx, ha, hb, hv, top_k=30, n_boot=300, seed=DEFAULT_R16_SEED
        )
        r16 = compute_joint_cross_mask_r16(
            obs, pred, tm, t_idx, source_grand_mean=None,
            delta_half_a=ha, delta_half_b=hb, half_valid=hv,
            top_k=30, n_boot=300, seed=DEFAULT_R16_SEED,
        )
        self.assertLessEqual(abs(r12["estimate"] - r16["estimate"]), 1e-12)
        self.assertLessEqual(abs(r12["ci975_lo"] - r16["ci975_lo"]), 1e-12)
        self.assertLessEqual(abs(r12["ci975_hi"] - r16["ci975_hi"]), 1e-12)
        self.assertLessEqual(
            abs(r12["split_summary"]["estimate"] - r16["split_summary"]["estimate"]), 1e-12
        )
        self.assertLessEqual(
            abs(r12["all_common_summary"]["estimate"] - r16["all_common_summary"]["estimate"]), 1e-12
        )

    def test_2_property_p1_dorth_null_invariance_by_construction(self) -> None:
        """Audit I3: Null-space invariance of secondary `r_dorth_cal` verified here as a unit test (<= 1e-14)."""
        res = verify_property_p1_synthetic(seed=DEFAULT_R16_SEED, n_random_affine=50, n_boot=150)
        self.assertTrue(res["property_p1_pass"], f"Property P1 failed: {res}")
        self.assertGreater(abs(res["named"]["Source_Grand_Mean"]["r_orth_cal"]), 0.15)
        self.assertLessEqual(abs(res["named"]["Source_Grand_Mean"]["r_dorth_cal"]), 1e-14)
        self.assertLessEqual(abs(res["named"]["Transfer_Mean_GlobalScale"]["r_dorth_cal"]), 1e-14)
        self.assertLessEqual(abs(res["named"]["Transfer_Mean"]["r_dorth_cal"]), 1e-14)
        self.assertLessEqual(res["max_abs_random_affine_r_dorth_cal"], 1e-14)

    def test_3_mechanism_m2_synthetic_illustration(self) -> None:
        """Audit I4: Supporting synthetic illustration for secondary descriptive Mechanism Hypothesis M2."""
        res = illustrate_mechanism_m2_synthetic(seed=DEFAULT_R16_SEED, n_boot=250)
        self.assertTrue(res["illustration_consistent"], f"M2 illustration failed: {res}")
        self.assertGreater(res["weight_shift_ci975_lo"], 0.0)
        self.assertGreater(res["e1_orth_jump_ci975_lo"], 0.0)

    def test_4_zero_leakage_and_triple_hard_execution_gates(self) -> None:
        """Leakage guard blocks training/H1 overlaps; triple hard gate blocks pre-gate access (`I9`, `I10`, Ruling 7)."""
        assert_zero_training_leakage("STATE_SE", "nadig_hepg2", checkpoint_id="zeroshot/hepg2")
        assert_zero_training_leakage("STATE_HVG", "a549_ko", checkpoint_id="zeroshot/ensemble4")
        assert_zero_training_leakage("STATE_SE", "dld1", checkpoint_id="zeroshot/ensemble4")
        assert_zero_training_leakage(
            "STATE_SE", "gse295214_h1", checkpoint_id="zeroshot/ensemble4",
            source_pool_lines=GROUP3_PRIMARY_FCT_SOURCES,
        )

        # Group 3 H1 lineage leakage guard: `vcc2025_h1` MUST be blocked from primary FCT pool for `gse295214_h1`
        with self.assertRaises(TrainingLeakageError):
            assert_zero_training_leakage(
                "STATE_SE", "gse295214_h1", checkpoint_id="zeroshot/ensemble4",
                source_pool_lines=GROUP3_PRIMARY_FCT_SOURCES + ("vcc2025_h1",),
            )
        with self.assertRaises(TrainingLeakageError):
            assert_zero_training_leakage("STATE_SE", "nadig_jurkat", checkpoint_id="zeroshot/hepg2")
        with self.assertRaises(TrainingLeakageError):
            assert_zero_training_leakage("STATE_HVG", "jurkat", checkpoint_id="zeroshot/ensemble4")
        with self.assertRaises(TrainingLeakageError):
            assert_zero_training_leakage("STATE_SE", "replogle_rpe1", checkpoint_id="zeroshot/k562")
        with self.assertRaises(TrainingLeakageError):
            assert_zero_training_leakage(
                "CellularAI_OpenCode_PublicRefit", "vcc2025_h1",
                source_pool_lines=["nadig_hepg2", "vcc2025_h1"],
            )

        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            fake_models = tdir / "no_models_here"
            with self.assertRaises(BlindingGateError):
                assert_prereg_i_locked(tdir)
            with self.assertRaises(BlindingGateError):
                assert_code_freeze_r16(tdir)
            with self.assertRaises(BlindingGateError):
                assert_r16_go(tdir)
            with self.assertRaises(BlindingGateError):
                assert_all_execution_gates_r16(tdir)

            adapter = ArcStateAdapter(
                "STATE_SE", checkpoint_id="zeroshot/ensemble4",
                models_dir=fake_models, status_dir=tdir,
            )
            # Missing all 3 -> blocked in BOTH adapter & run_real_evaluation_r16
            with self.assertRaises(BlindingGateError):
                adapter.predict_deltas("a549_ko", "a549_ko", np.ones((32, 5)), ["G1"], ["A", "B", "C", "D", "E"])
            with self.assertRaises(BlindingGateError):
                run_real_evaluation_r16(["a549_ko"], status_dir=tdir, models_dir=fake_models)

            # 1 of 3 (`PREREG_I_COMMIT.txt`) -> STILL blocked
            (tdir / "PREREG_I_COMMIT.txt").write_text("LOCKED_TEST\n", encoding="utf-8")
            self.assertEqual(assert_prereg_i_locked(tdir), tdir / "PREREG_I_COMMIT.txt")
            with self.assertRaises(BlindingGateError):
                adapter.predict_deltas("a549_ko", "a549_ko", np.ones((32, 5)), ["G1"], ["A", "B", "C", "D", "E"])

            # 2 of 3 (`+ CODE_FREEZE_R16.txt`, missing `R16_GO.txt`) -> STILL blocked
            (tdir / "CODE_FREEZE_R16.txt").write_text("FROZEN_TEST\n", encoding="utf-8")
            self.assertEqual(assert_code_freeze_r16(tdir), tdir / "CODE_FREEZE_R16.txt")
            with self.assertRaises(BlindingGateError):
                adapter.predict_deltas("a549_ko", "a549_ko", np.ones((32, 5)), ["G1"], ["A", "B", "C", "D", "E"])
            with self.assertRaises(BlindingGateError):
                run_real_evaluation_r16(["a549_ko"], status_dir=tdir, models_dir=fake_models)

            # Empty `R16_GO.txt` -> STILL blocked
            (tdir / "R16_GO.txt").write_text("", encoding="utf-8")
            with self.assertRaises(BlindingGateError):
                assert_all_execution_gates_r16(tdir)

            # All 3 non-empty -> passes triple gate and proceeds to missing checkpoint FileNotFoundError
            (tdir / "R16_GO.txt").write_text("GO_TEST\n", encoding="utf-8")
            gates = assert_all_execution_gates_r16(tdir)
            self.assertEqual(set(gates.keys()), {"prereg_lock", "code_freeze", "r16_go"})
            with self.assertRaises(FileNotFoundError):
                adapter.predict_deltas("a549_ko", "a549_ko", np.ones((32, 5)), ["G1"], ["A", "B", "C", "D", "E"])

            cell_ai = CellularAIOpenCodePublicRefitAdapter(status_dir=tdir)
            with self.assertRaises(TrainingLeakageError):
                cell_ai.predict_deltas(
                    "nadig_hepg2", "nadig_hepg2", np.ones((10, 5)), ["G1"], ["A", "B", "C", "D", "E"],
                    source_pool_lines=["nadig_hepg2", "nadig_jurkat"],
                )

    def test_5_h_i1_m1_m2_m3_tost_and_two_stratum_adjudication(self) -> None:
        """Audit I6, I14, M1, M2, M3 & 2026-11-25 tie ruling: verify `rho_ab_orth`, `within_delta`, SUB_DELTA, EQUIV split."""
        mb = {"estimate": 0.08, "ci975_lo": 0.05, "ci975_hi": 0.11, "ci95_lo": 0.055, "ci95_hi": 0.105, "valid_frac": 1.0}
        mb_small = {"estimate": 0.015, "ci975_lo": 0.004, "ci975_hi": 0.026, "ci95_lo": 0.006, "ci95_hi": 0.024, "valid_frac": 1.0}
        eq = {"estimate": 0.005, "ci975_lo": -0.020, "ci975_hi": 0.030, "ci95_lo": -0.015, "ci95_hi": 0.025, "valid_frac": 1.0}
        inc = {"estimate": 0.020, "ci975_lo": -0.025, "ci975_hi": 0.065, "ci95_lo": -0.018, "ci95_hi": 0.058, "valid_frac": 1.0}
        mis = {"estimate": -0.07, "ci975_lo": -0.10, "ci975_hi": -0.04, "ci95_lo": -0.095, "ci95_hi": -0.045, "valid_frac": 1.0}

        low_vf = {"estimate": 0.001, "ci975_lo": -0.01, "ci975_hi": 0.01, "ci95_lo": -0.008, "ci95_hi": 0.008, "valid_frac": 0.50}
        c_low_vf = classify_pair_tost_r16(low_vf, DEFAULT_TOST_DELTA)
        self.assertEqual(c_low_vf["category"], "INCONCLUSIVE")
        self.assertFalse(c_low_vf["within_delta"])

        nan_ci = {"estimate": 0.0, "ci975_lo": float("nan"), "ci975_hi": 0.01, "ci95_lo": -0.01, "ci95_hi": 0.01, "valid_frac": 1.0}
        c_nan = classify_pair_tost_r16(nan_ci, DEFAULT_TOST_DELTA)
        self.assertEqual(c_nan["category"], "INCONCLUSIVE")
        self.assertFalse(c_nan["within_delta"])

        for tie in (
            {"estimate": 0.02, "ci975_lo": 0.0, "ci975_hi": 0.05, "ci95_lo": 0.005, "ci95_hi": 0.045, "valid_frac": 1.0},
            {"estimate": 0.0, "ci975_lo": -0.05, "ci975_hi": 0.05, "ci95_lo": -0.04, "ci95_hi": 0.02, "valid_frac": 1.0},
            {"estimate": 0.0, "ci975_lo": -0.05, "ci975_hi": 0.05, "ci95_lo": -0.02, "ci95_hi": 0.04, "valid_frac": 1.0},
        ):
            ct = classify_pair_tost_r16(tie, DEFAULT_TOST_DELTA)
            self.assertEqual(ct["category"], "INCONCLUSIVE")
            self.assertFalse(ct["within_delta"])

        tie_lo_inside = {
            "estimate": 0.015, "ci975_lo": 0.0, "ci975_hi": 0.030,
            "ci95_lo": 0.003, "ci95_hi": 0.027, "valid_frac": 1.0,
        }
        tie_hi_inside = {
            "estimate": -0.015, "ci975_lo": -0.030, "ci975_hi": 0.0,
            "ci95_lo": -0.027, "ci95_hi": -0.003, "valid_frac": 1.0,
        }
        tie_both_zero = {
            "estimate": 0.0, "ci975_lo": 0.0, "ci975_hi": 0.0,
            "ci95_lo": 0.0, "ci95_hi": 0.0, "valid_frac": 1.0,
        }
        for tie_in in (tie_lo_inside, tie_hi_inside, tie_both_zero):
            ct_in = classify_pair_tost_r16(tie_in, DEFAULT_TOST_DELTA)
            self.assertEqual(ct_in["category"], "INCONCLUSIVE")
            self.assertTrue(ct_in["within_delta"])

        orth_pos = {"ci975_lo": 0.05, "ci975_hi": 0.25, "valid_frac": 1.0, "rho_ab_orth": 0.6}
        orth_null = {"ci975_lo": -0.02, "ci975_hi": 0.04, "valid_frac": 1.0, "rho_ab_orth": 0.6}
        seen = ("nadig_hepg2", "nadig_jurkat")
        unseen = ("cd4t_rest", "dld1", "feng_ipsc")

        r_tie_line = adjudicate_claim_i1_line(
            orth_pos,
            {"per_source": {**{k: eq for k in seen}, "cd4t_rest": eq, "dld1": eq, "feng_ipsc": tie_lo_inside}},
            seen_sources=seen, unseen_sources=unseen, delta=DEFAULT_TOST_DELTA,
        )
        self.assertEqual(r_tie_line["line_summary"], "MIXED_OR_INCONCLUSIVE")
        self.assertEqual(r_tie_line["counts"]["INCONCLUSIVE"], 1)

        for bad_orth in (
            {"ci975_lo": 0.05, "ci975_hi": 0.25, "valid_frac": 1.0, "rho_ab_orth": 0.0},
            {"ci975_lo": 0.05, "ci975_hi": 0.25, "valid_frac": 1.0, "rho_ab_orth": -0.15},
            {"ci975_lo": 0.05, "ci975_hi": 0.25, "valid_frac": 1.0, "rho_ab": 0.90},
            {"ci975_lo": 0.05, "ci975_hi": 0.25, "valid_frac": 0.90, "rho_ab_orth": 0.6},
        ):
            r_bad_mb = adjudicate_claim_i1_line(
                bad_orth, {"per_source": {k: mb for k in seen + unseen}},
                seen_sources=seen, unseen_sources=unseen, delta=DEFAULT_TOST_DELTA,
            )
            self.assertEqual(r_bad_mb["line_summary"], "INCONCLUSIVE")
            self.assertFalse(r_bad_mb["matched_arm_valid"])

            r_bad_eq = adjudicate_claim_i1_line(
                bad_orth, {"per_source": {k: eq for k in seen + unseen}},
                seen_sources=seen, unseen_sources=unseen, delta=DEFAULT_TOST_DELTA,
            )
            self.assertEqual(r_bad_eq["line_summary"], "INCONCLUSIVE")

        r_sub = adjudicate_claim_i1_line(
            orth_pos, {"per_source": {k: mb_small for k in seen + unseen}},
            seen_sources=seen, unseen_sources=unseen, delta=DEFAULT_TOST_DELTA,
        )
        self.assertEqual(r_sub["line_summary"], "LINE_MATCHED_SPECIFIC_SUB_DELTA")
        self.assertTrue(r_sub["all_within_delta"])
        self.assertTrue(r_sub["any_within_delta"])

        r_mix_mb = adjudicate_claim_i1_line(
            orth_pos, {"per_source": {**{k: mb_small for k in seen}, **{k: mb for k in unseen}}},
            seen_sources=seen, unseen_sources=unseen, delta=DEFAULT_TOST_DELTA,
        )
        self.assertEqual(r_mix_mb["line_summary"], "LINE_MATCHED_SPECIFIC")
        self.assertFalse(r_mix_mb["all_within_delta"])
        self.assertTrue(r_mix_mb["any_within_delta"])

        r_eq_nz = adjudicate_claim_i1_line(
            orth_pos, {"per_source": {k: eq for k in seen + unseen}},
            seen_sources=seen, unseen_sources=unseen, delta=DEFAULT_TOST_DELTA,
        )
        self.assertEqual(r_eq_nz["line_summary"], "LINE_EQUIVALENT_NONZERO_SIGNAL")

        r_eq_null = adjudicate_claim_i1_line(
            orth_null, {"per_source": {k: eq for k in seen + unseen}},
            seen_sources=seen, unseen_sources=unseen, delta=DEFAULT_TOST_DELTA,
        )
        self.assertEqual(r_eq_null["line_summary"], "LINE_EQUIVALENT_NULL_SIGNAL")

        r_mis = adjudicate_claim_i1_line(
            orth_pos, {"per_source": {k: mis for k in seen + unseen}},
            seen_sources=seen, unseen_sources=unseen, delta=DEFAULT_TOST_DELTA,
        )
        self.assertEqual(r_mis["line_summary"], "LINE_MISMATCHED_BETTER")

        r_split = adjudicate_claim_i1_line(
            orth_pos, {"per_source": {**{k: mb for k in seen}, **{k: eq for k in unseen}}},
            seen_sources=seen, unseen_sources=unseen, delta=DEFAULT_TOST_DELTA,
        )
        self.assertEqual(r_split["line_summary"], "MIXED_OR_INCONCLUSIVE")

        r_inc = adjudicate_claim_i1_line(
            orth_pos, {"per_source": {**{k: eq for k in seen}, "cd4t_rest": eq, "dld1": eq, "feng_ipsc": inc}},
            seen_sources=seen, unseen_sources=unseen, delta=DEFAULT_TOST_DELTA,
        )
        self.assertEqual(r_inc["line_summary"], "MIXED_OR_INCONCLUSIVE")

        self.assertEqual(
            adjudicate_claim_i1_line(
                orth_pos, {"per_source": {k: mb for k in seen}},
                seen_sources=seen, unseen_sources=(), delta=DEFAULT_TOST_DELTA,
            )["line_summary"],
            "MISSING_STRATUM",
        )
        self.assertEqual(
            adjudicate_claim_i1_line(
                orth_pos, {"per_source": {k: mb for k in unseen}},
                seen_sources=(), unseen_sources=unseen, delta=DEFAULT_TOST_DELTA,
            )["line_summary"],
            "MISSING_STRATUM",
        )

        with self.assertRaises(ValueError):
            adjudicate_claim_i1_line(
                orth_pos, {"per_source": {"nadig_hepg2": mb}},
                seen_sources=("nadig_hepg2",), unseen_sources=("nadig_hepg2",),
            )
        with self.assertRaises(ValueError):
            adjudicate_claim_i1_line(
                orth_pos, {"per_source": {"nadig_hepg2": mb, "unexpected": mb}},
                seen_sources=("nadig_hepg2",), unseen_sources=("dld1",),
            )

        g2_ok = adjudicate_group_r16(
            "group2",
            {"a549_ko": "LINE_MATCHED_SPECIFIC", "dld1": "MIXED_OR_INCONCLUSIVE", "cd4t_rest": "LINE_EQUIVALENT_NONZERO_SIGNAL"},
            {"a549_ko": "CRISPRko", "dld1": "CRISPRi", "cd4t_rest": "CRISPRi"},
        )
        self.assertTrue(g2_ok["line_level_verdicts_awarded"])
        self.assertEqual(g2_ok["by_modality"]["CRISPRko"]["k_over_K"]["LINE_MATCHED_SPECIFIC"], "1/1")
        self.assertEqual(g2_ok["by_modality"]["CRISPRi"]["K"], 2)

        g2_deg = adjudicate_group_r16(
            "group2",
            {"a549_ko": "LINE_MATCHED_SPECIFIC"},
            {"a549_ko": "CRISPRko"},
        )
        self.assertFalse(g2_deg["line_level_verdicts_awarded"])
        self.assertEqual(g2_deg["group_status"], "SINGLE_STRATUM_DESCRIPTIVE_ONLY")

    def test_6_i14_q4_group3_local_mismatched_strata_builder(self) -> None:
        """Audit I14, Q4, & Group 3: verify Group 1, Group 2, and Group 3 `gse295214_h1` strata."""
        self.assertEqual(GROUP2_UNSEEN_PRIMARY_LINES, ("a549_ko", "dld1", "cd4t_rest"))

        g2_a549 = build_mismatched_strata_r16("a549_ko", "zeroshot/ensemble4")
        self.assertEqual(g2_a549["seen_by_checkpoint"], ("nadig_hepg2", "nadig_jurkat"))
        self.assertEqual(g2_a549["unseen_by_checkpoint"], ("cd4t_rest", "dld1", "feng_ipsc"))

        g1_hepg2 = build_mismatched_strata_r16("nadig_hepg2", "zeroshot/hepg2")
        self.assertEqual(g1_hepg2["seen_by_checkpoint"], ("nadig_jurkat",))
        self.assertEqual(g1_hepg2["unseen_by_checkpoint"], ("a549_ko", "cd4t_rest", "dld1", "feng_ipsc"))

        g1_jurkat = build_mismatched_strata_r16("nadig_jurkat", "zeroshot/jurkat")
        self.assertEqual(g1_jurkat["seen_by_checkpoint"], ("nadig_hepg2",))
        self.assertEqual(g1_jurkat["unseen_by_checkpoint"], ("a549_ko", "cd4t_rest", "dld1", "feng_ipsc"))

        # Group 3 prospective holdout (`PREREG_I v5 §2.4`)
        g3_h1 = build_mismatched_strata_r16("gse295214_h1", "zeroshot/ensemble4")
        self.assertEqual(g3_h1["seen_by_checkpoint"], ("nadig_hepg2", "nadig_jurkat"))
        self.assertEqual(g3_h1["unseen_by_checkpoint"], ("a549_ko", "cd4t_rest", "dld1", "feng_ipsc"))

        for bad in ("nourreddine_kolf21j", "replogle_k562_essential", "replogle_rpe1", "vcc2025_h1"):
            with self.assertRaises(ValueError):
                build_mismatched_strata_r16(
                    "a549_ko", "zeroshot/ensemble4",
                    group2_qualifying_lines=("a549_ko", "dld1", bad),
                )

    def test_7_paired_contrasts_and_h_i1_h_i2_end_to_end(self) -> None:
        """End-to-end check: `rho_ab_orth` propagation, stratified `H-I1`, and all 3 `H-I2` branches."""
        rng = np.random.default_rng(777)
        P, G = 80, 400
        tm = rng.normal(0.0, 0.2, size=(P, G))
        true_spec = rng.normal(0.0, 0.3, size=(P, G))
        ha = tm + true_spec + rng.normal(0.0, 0.04, size=(P, G))
        hb = tm + true_spec + rng.normal(0.0, 0.04, size=(P, G))
        obs = 0.5 * (ha + hb)
        t_idx = np.arange(P, dtype=np.int64)
        hv = np.ones(P, dtype=bool)

        pred_good = tm + 0.8 * true_spec + rng.normal(0.0, 0.05, size=(P, G))
        pred_fct = tm + 0.4 * true_spec + rng.normal(0.0, 0.08, size=(P, G))
        strata = build_mismatched_strata_r16("a549_ko", "zeroshot/ensemble4")
        all_mis_names = strata["seen_by_checkpoint"] + strata["unseen_by_checkpoint"]
        mis_preds = {
            name: tm + 0.05 * true_spec + rng.normal(0.0, 0.25, size=(P, G))
            for name in all_mis_names
        }

        orth_good = compute_orth_suite_r16(
            obs, pred_good, tm, t_idx, source_grand_mean=None,
            delta_half_a=ha, delta_half_b=hb, half_valid=hv, n_boot=300, seed=777
        )
        self.assertIn("rho_ab_orth", orth_good)
        self.assertNotIn("rho_ab", orth_good)
        self.assertGreater(orth_good["rho_ab_orth"], 0.5)

        paired_mis = compute_paired_model_vs_all_mismatched_r16(
            obs, pred_good, mis_preds, tm, t_idx, source_grand_mean=None,
            delta_half_a=ha, delta_half_b=hb, half_valid=hv, n_boot=300, seed=777
        )
        self.assertIn("rho_ab_orth", paired_mis)
        self.assertNotIn("rho_ab", paired_mis)

        c1 = adjudicate_claim_i1_line(
            orth_good, paired_mis,
            seen_sources=strata["seen_by_checkpoint"],
            unseen_sources=strata["unseen_by_checkpoint"],
            delta=DEFAULT_TOST_DELTA,
        )
        self.assertEqual(c1["line_summary"], "LINE_MATCHED_SPECIFIC")
        self.assertEqual(c1["strata"]["seen_by_checkpoint"]["counts"]["MATCHED_BETTER"], 2)
        self.assertEqual(c1["strata"]["unseen_by_checkpoint"]["counts"]["MATCHED_BETTER"], 3)

        vs_fct = compute_paired_contrast_r16(
            obs, pred_good, pred_fct, tm, t_idx, source_grand_mean=None,
            delta_half_a=ha, delta_half_b=hb, half_valid=hv, n_boot=300, seed=777
        )
        self.assertIn("rho_ab_orth", vs_fct)
        self.assertEqual(adjudicate_claim_i2_line(vs_fct)["classification"], "MODEL_SUPERIOR_TO_FCT")

        fct_better = compute_paired_contrast_r16(
            obs, pred_fct, pred_good, tm, t_idx, source_grand_mean=None,
            delta_half_a=ha, delta_half_b=hb, half_valid=hv, n_boot=300, seed=777
        )
        self.assertEqual(adjudicate_claim_i2_line(fct_better)["classification"], "FCT_SUPERIOR_TO_MODEL")

        vs_self = compute_paired_contrast_r16(
            obs, pred_fct, pred_fct, tm, t_idx, source_grand_mean=None,
            delta_half_a=ha, delta_half_b=hb, half_valid=hv, n_boot=300, seed=777
        )
        self.assertEqual(adjudicate_claim_i2_line(vs_self)["classification"], "INDISTINGUISHABLE")
        self.assertEqual(
            adjudicate_claim_i2_line({"estimate": 0.2, "ci975_lo": 0.1, "ci975_hi": 0.3, "valid_frac": 0.90})["classification"],
            "INDISTINGUISHABLE",
        )
        self.assertEqual(
            adjudicate_claim_i2_line({"estimate": 0.2, "ci975_lo": float("nan"), "ci975_hi": 0.3, "valid_frac": 1.0})["classification"],
            "INDISTINGUISHABLE",
        )

    def test_8_m5_stratified_control_sampler_and_pseudobulk_shift(self) -> None:
        """Audit M5.1: 512 batch-stratified shared 16x32 sets, <512 truncation to mult-of-32, <128 exclusion."""
        idx1 = sample_control_indices(1000, 512, seed=DEFAULT_R16_SEED)
        idx2 = sample_control_indices(1000, 512, seed=DEFAULT_R16_SEED)
        self.assertTrue(np.array_equal(idx1, idx2))
        self.assertEqual(len(np.unique(idx1)), 512)

        batches = ["B0"] * 500 + ["B1"] * 250 + ["B2"] * 150 + ["B3"] * 100
        sets_a = sample_stratified_control_sets_r16(1000, batch_labels=batches, seed=DEFAULT_R16_SEED)
        sets_b = sample_stratified_control_sets_r16(1000, batch_labels=batches, seed=DEFAULT_R16_SEED)
        self.assertTrue(np.array_equal(sets_a, sets_b))
        self.assertEqual(sets_a.shape, (16, 32))
        flat = sets_a.ravel()
        self.assertEqual(len(np.unique(flat)), 512)
        sampled_batches = [batches[i] for i in flat]
        self.assertEqual(sampled_batches.count("B0"), 256)
        self.assertEqual(sampled_batches.count("B1"), 128)
        self.assertEqual(sampled_batches.count("B2"), 77)
        self.assertEqual(sampled_batches.count("B3"), 51)

        sets_small = sample_stratified_control_sets_r16(200, seed=DEFAULT_R16_SEED)
        self.assertEqual(sets_small.shape, (6, 32))
        self.assertEqual(len(np.unique(sets_small.ravel())), 192)

        with self.assertRaises(InsufficientControlCellsError):
            sample_stratified_control_sets_r16(127, seed=DEFAULT_R16_SEED)

        ctrl = np.array([[10.0, 20.0, 70.0], [30.0, 10.0, 60.0]])
        shift = compute_pseudobulk_log1p_shift(3.0 * ctrl, ctrl, input_is_log1p=False)
        self.assertLessEqual(float(np.max(np.abs(shift))), 1e-13)

    def test_9_hvg_preprocessing_and_embedding_diagnostics(self) -> None:
        """Verify official `normalize_total(median)->log1p` primary vs `libnorm(1e4)->log1p` & centroid diagnostics."""
        raw = np.array(
            [
                [10.0, 20.0, 70.0, 0.0],   # sum = 100
                [40.0, 80.0, 280.0, 0.0],  # sum = 400
                [20.0, 40.0, 140.0, 0.0],  # sum = 200 -> median = 200
            ],
            dtype=np.float32,
        )
        raw_genes = ["G1", "G2", "G3", "G4"]
        g_model = ["G2", "MISSING_GENE", "G1"]

        x_off, a_off = preprocess_hvg_counts_r16(raw, raw_genes, g_model, norm_mode="official_median_log1p")
        self.assertAlmostEqual(a_off["target_sum"], 200.0, places=6)
        self.assertEqual(a_off["n_g_model_mapped"], 2)
        self.assertTrue(np.allclose(x_off[:, 1], 0.0))
        # All 3 proportional rows map identically to log1p(200 * [0.2, 0.1]) = log1p([40, 20])
        self.assertTrue(np.allclose(x_off[:, 0], np.log1p(40.0), atol=1e-5))
        self.assertTrue(np.allclose(x_off[:, 2], np.log1p(20.0), atol=1e-5))

        x_1e4, a_1e4 = preprocess_hvg_counts_r16(raw, raw_genes, g_model, norm_mode="libnorm1e4_log1p")
        self.assertAlmostEqual(a_1e4["target_sum"], 1e4, places=6)
        self.assertTrue(np.allclose(x_1e4[:, 0], np.log1p(2000.0), atol=1e-4))
        self.assertTrue(np.allclose(x_1e4[:, 2], np.log1p(1000.0), atol=1e-4))

        e_a = np.hstack([np.ones((8, 2048)), np.zeros((8, 10))])
        e_b = np.hstack([-np.ones((8, 2048)), np.ones((8, 10))])
        diag = compute_embedding_separation_diagnostics(e_a, e_b)
        self.assertAlmostEqual(diag["centroid_cos_2048"], -1.0, places=6)
        self.assertGreater(diag["centroid_l2_2048"], 1.0)
        self.assertGreater(diag["centroid_l2_full"], diag["centroid_l2_2048"])

    def test_10_gse295214_pre_download_synthetic_10x_h5_reader(self) -> None:
        """Verify pre-download Group 3 `GSE295214` 10x `.h5` reader strictly on synthetic fixtures (`§2.4`)."""
        hgnc = _StubHGNC({"TP53", "MYC", "HLA-A", "BRCA1"})
        self.assertEqual(parse_gse295214_guide_r16("non-targeting_0001", "non-targeting_0001", hgnc)[0], "ntc")
        self.assertEqual(parse_gse295214_guide_r16("safe_targeting_1", "safe_targeting_1", hgnc)[0], "exclude")
        self.assertEqual(parse_gse295214_guide_r16("TP53_MYC_1", "TP53_MYC_1", hgnc)[0], "exclude")
        t_kind, t_info = parse_gse295214_guide_r16("HLA-A_P1P2_1", "HLA-A_P1P2_1", hgnc)
        self.assertEqual(t_kind, "target")
        self.assertEqual(t_info[1], "HLA-A")
        d_kind, d_info = parse_gse295214_guide_r16("TP53_1|TP53_2", "TP53_1|TP53_2", hgnc)
        self.assertEqual(d_kind, "target")
        self.assertEqual(d_info[1], "TP53")
        self.assertEqual(parse_gse295214_guide_r16("TP53_1|MYC_1", "TP53_1|MYC_1", hgnc)[0], "exclude")

        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            n_genes = 260  # >= 200 genes so QC passes; include high-MT genes to verify max_pct_mt=None!
            gene_names = ["MT-CO1", "MT-ND1"] + [f"G{i:04d}" for i in range(n_genes - 2)]
            guide_names = [
                "non-targeting_0001",
                "non-targeting_0002",
                "TP53_P1P2_1",
                "MYC_1",
                "safe_targeting_01",
            ]
            # 4 GSM channels: odd (1, 3) -> Half A; even (2, 4) -> Half B
            # Note: GSM8943716 has only 5 NTCs (< 15) to verify within-half fallback `half_2_pooled`!
            chan_specs = [
                ("GSM8943713", 45, 10),  # Half A (>=15 NTC)
                ("GSM8943714", 45, 10),  # Half B (>=15 NTC)
                ("GSM8943715", 45, 10),  # Half A (>=15 NTC)
                ("GSM8943716", 5, 10),   # Half B (<15 NTC -> fallback to half_2_pooled!)
            ]
            for gsm, n_ntc_c, n_tgt_c in chan_specs:
                # Extra 3 bad cells per channel: (a) low UMI/genes, (b) doublet 2 guides, (c) safe-targeting
                n_cells = n_ntc_c + n_tgt_c + 3
                rna = np.zeros((n_genes, n_cells), dtype=np.float32)
                gd = np.zeros((len(guide_names), n_cells), dtype=np.float32)

                # Valid cells: 220 expressed genes, total UMI >= 400, WITH 60% MT UMIs (must PASS since max_pct_mt=None!)
                rna[2:222, : n_ntc_c + n_tgt_c + 2] = 2.0  # 440 UMIs, 220 genes
                rna[0, : n_ntc_c + n_tgt_c + 2] = 660.0    # 60% MT-CO1 UMI!
                # Cell -1: low UMI (100) and low genes (50) -> MUST be filtered out
                rna[2:52, -1] = 2.0

                # NTC singlets: top1 = 15, ambient background = 1
                gd[0, :n_ntc_c] = 15.0
                gd[1, :n_ntc_c] = 1.0
                # Targeting singlets: TP53 & MYC
                half_t = n_tgt_c // 2
                gd[2, n_ntc_c : n_ntc_c + half_t] = 18.0
                gd[3, n_ntc_c + half_t : n_ntc_c + n_tgt_c] = 16.0
                # Doublet cell: two guides above threshold -> MUST be rejected (n_pos == 2)
                gd[2, n_ntc_c + n_tgt_c] = 15.0
                gd[3, n_ntc_c + n_tgt_c] = 14.0
                # Safe-targeting cell -> MUST be excluded
                gd[4, n_ntc_c + n_tgt_c + 1] = 20.0
                # Low-QC cell with valid guide -> MUST be excluded by RNA QC
                gd[0, -1] = 20.0

                bcs = [f"{gsm}_BC{i:04d}-1" for i in range(n_cells)]
                _write_synthetic_10x_h5(tdir / f"{gsm}_sample.h5", rna, gd, gene_names, guide_names, bcs)

            parts, audit = build_gse295214_parts(tdir, hgnc=hgnc, use_r10_partspec=False)
            self.assertEqual(audit["n_channels"], 4)
            self.assertEqual(audit["qc_spec"]["max_pct_mt"], None)
            self.assertEqual(audit["n_kept_ntc"], 45 + 45 + 45 + 5)  # 140 NTCs (all high-MT passed!)
            self.assertEqual(audit["n_kept_targeting"], 40)
            self.assertEqual(audit["n_kept_half_a_ntc"], 90)
            self.assertEqual(audit["n_kept_half_b_ntc"], 50)
            self.assertTrue(audit["per_half_ntc_ge15"])
            # GSM8943716 (5 NTC + 10 tgt = 15 cells) must fall back to `half_2_pooled` (never Half A!)
            self.assertEqual(audit["n_low_ntc_channel_fallback_cells"], 15)
            p_by_tag = {p.tag: p for p in parts}
            self.assertTrue(np.all(p_by_tag["GSM8943713"].half_flags == 1))
            self.assertTrue(np.all(p_by_tag["GSM8943714"].half_flags == 2))
            self.assertTrue(np.all(p_by_tag["GSM8943715"].half_flags == 1))
            self.assertTrue(np.all(p_by_tag["GSM8943716"].half_flags == 2))
            self.assertTrue(np.all(p_by_tag["GSM8943716"].batches == "half_2_pooled"))

            # NTC raw-count sampler on 140 NTCs -> truncates to largest multiple of 32 = 128 (`4 x 32`)
            raw_ntc, r_genes, _, c_aud = extract_gse295214_ntc_raw_counts(
                tdir, hgnc=hgnc, seed=DEFAULT_R16_SEED, parts_override=parts
            )
            self.assertEqual(raw_ntc.shape, (128, n_genes))
            self.assertEqual(r_genes, gene_names)
            self.assertEqual(c_aud["n_ctrl_total"], 140)
            self.assertEqual(c_aud["n_ctrl_sampled"], 128)

            # Non-integer RNA count rejection check
            bad_dir = tdir / "non_int"
            bad_dir.mkdir()
            rna_bad = np.full((n_genes, 5), 1.5, dtype=np.float32)
            gd_bad = np.full((len(guide_names), 5), 10.0, dtype=np.float32)
            _write_synthetic_10x_h5(bad_dir / "GSM9999999.h5", rna_bad, gd_bad, gene_names, guide_names, [f"B{i}" for i in range(5)])
            with self.assertRaises(ValueError):
                build_gse295214_parts(bad_dir, hgnc=hgnc, use_r10_partspec=False)

    def test_11_run_eval_r16_synthetic_smoke_and_independent_npz_parity(self) -> None:
        """Verify end-to-end `run_eval_r16` synthetic smoke + independent `.npz` recomputation parity <= 1e-12."""
        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp)
            res = run_synthetic_smoke_r16(out_dir, n_boot=250, seed=DEFAULT_R16_SEED)
            line_r = res["lines"]["synthetic_a549_ko"]
            self.assertEqual(line_r["primary"]["H_I1_adjudication"]["line_summary"], "LINE_MATCHED_SPECIFIC")
            self.assertIn("b0_vs_b17_mis_nadig_hepg2", line_r["secondary"]["g4_factorial_decomposition"])
            self.assertIn("bypass_zero", line_r["secondary"]["neutral_batch_robustness_panel"])

            ldir = out_dir / "arrays" / "synthetic_a549_ko"
            manifest = json.loads((ldir / "array_manifest_sha256.json").read_text(encoding="utf-8"))
            self.assertIn("observed.npz", manifest)
            self.assertIn("preds_STATE_SE_fixed_index_zero_matched.npz", manifest)

            # Independently reload exported `.npz` arrays and recompute `H-I1` & `H-I2` to <= 1e-12
            obs_npz = np.load(ldir / "observed.npz", allow_pickle=False)
            se_m_npz = np.load(ldir / "preds_STATE_SE_fixed_index_zero_matched.npz", allow_pickle=False)
            strata = build_mismatched_strata_r16("a549_ko", "zeroshot/ensemble4")
            all_mis = strata["seen_by_checkpoint"] + strata["unseen_by_checkpoint"]
            mis_reload = {
                sk: np.load(ldir / f"preds_STATE_SE_fixed_index_zero_mis_{sk}.npz", allow_pickle=False)["delta_pred"]
                for sk in all_mis
            }

            re_orth = compute_orth_suite_r16(
                obs_npz["delta_obs"], se_m_npz["delta_pred"], obs_npz["delta_transfer"],
                obs_npz["target_indices"], source_grand_mean=None,
                delta_half_a=obs_npz["delta_half_a"], delta_half_b=obs_npz["delta_half_b"],
                half_valid=obs_npz["half_valid"], n_boot=250, seed=DEFAULT_R16_SEED,
            )
            re_mis = compute_paired_model_vs_all_mismatched_r16(
                obs_npz["delta_obs"], se_m_npz["delta_pred"], mis_reload, obs_npz["delta_transfer"],
                obs_npz["target_indices"], source_grand_mean=None,
                delta_half_a=obs_npz["delta_half_a"], delta_half_b=obs_npz["delta_half_b"],
                half_valid=obs_npz["half_valid"], n_boot=250, seed=DEFAULT_R16_SEED,
            )
            re_adj = adjudicate_claim_i1_line(
                re_orth, re_mis,
                seen_sources=strata["seen_by_checkpoint"],
                unseen_sources=strata["unseen_by_checkpoint"],
                delta=DEFAULT_TOST_DELTA,
            )
            self.assertLessEqual(
                abs(re_orth["estimate"] - line_r["primary"]["STATE_SE_matched_orth"]["estimate"]), 1e-12
            )
            self.assertEqual(re_adj["line_summary"], line_r["primary"]["H_I1_adjudication"]["line_summary"])

    @unittest.skipUnless(_HAS_ARC_STATE, "Requires round16/.venv with arc-state + torch installed")
    def test_12_arc_state_adapter_real_hvg_weights_on_pure_synthetic_cells(self) -> None:
        """
        Verify `ArcStateAdapter('STATE_HVG')` against real downloaded `ST-HVG-Replogle/zeroshot/hepg2`
        weights using 100% synthetic Poisson count tensors (ZERO real cells!), checking `M5.3`
        non-zero `X_ctrl` sensitivity, all 4 neutral batch modes, and both HVG normalizations.
        """
        os.environ.setdefault("WANDB_MODE", "disabled")
        os.environ.setdefault("HF_HUB_OFFLINE", "1")
        os.environ.setdefault("TRANSFORMERS_OFFLINE", "1")
        ckpt = DEFAULT_MODELS_DIR / "ST-HVG-Replogle" / "zeroshot" / "hepg2" / "checkpoints" / "final.ckpt"
        if not ckpt.is_file():
            self.skipTest(f"ST-HVG checkpoint not present at {ckpt}")

        import pickle

        with open(ckpt.parent.parent / "var_dims.pkl", "rb") as fh:
            g_model = list(pickle.load(fh)["gene_names"])

        rng = np.random.default_rng(DEFAULT_R16_SEED)
        # Purely synthetic Poisson counts (64 cells = 2 x 32 sets, 2,000 G_model genes)
        rates_a = rng.uniform(0.1, 3.0, size=(1, len(g_model)))
        rates_b = rng.uniform(0.1, 3.0, size=(1, len(g_model)))
        raw_a = rng.poisson(rates_a, size=(64, len(g_model))).astype(np.float32)
        raw_b = rng.poisson(rates_b, size=(64, len(g_model))).astype(np.float32)
        perts = ["AAAS", "AAMP", "AARS"]
        eval_genes = g_model[:20]

        preds_by_mode = {}
        for bmode in ("fixed_index_zero", "bypass_zero", "mean_batch_emb", "multi_batch_ens8"):
            ad_a = ArcStateAdapter(
                "STATE_HVG", checkpoint_id="zeroshot/hepg2",
                models_dir=DEFAULT_MODELS_DIR, device="cpu", batch_mode=bmode,
                hvg_norm_mode="official_median_log1p",
            )
            out_a = ad_a.predict_deltas(
                target_line="nadig_hepg2", control_line_label="synth_A",
                ctrl_cells=raw_a, pert_names=perts, gene_names=eval_genes,
                ctrl_raw_gene_names=g_model, allow_synthetic_without_freeze=True,
            )
            self.assertEqual(out_a.shape, (3, 20))
            self.assertTrue(np.isfinite(out_a).all())
            preds_by_mode[bmode] = out_a

        # M5.3 synthetic X_ctrl sensitivity check under primary `fixed_index_zero`
        ad_b = ArcStateAdapter(
            "STATE_HVG", checkpoint_id="zeroshot/hepg2",
            models_dir=DEFAULT_MODELS_DIR, device="cpu", batch_mode="fixed_index_zero",
            hvg_norm_mode="official_median_log1p",
        )
        out_b = ad_b.predict_deltas(
            target_line="nadig_hepg2", control_line_label="synth_B",
            ctrl_cells=raw_b, pert_names=perts, gene_names=eval_genes,
            ctrl_raw_gene_names=g_model, allow_synthetic_without_freeze=True,
        )
        self.assertGreater(float(np.linalg.norm(preds_by_mode["fixed_index_zero"] - out_b)), 1e-4)

        # Secondary `libnorm1e4_log1p` normalization runs cleanly and yields finite outputs
        ad_1e4 = ArcStateAdapter(
            "STATE_HVG", checkpoint_id="zeroshot/hepg2",
            models_dir=DEFAULT_MODELS_DIR, device="cpu", batch_mode="fixed_index_zero",
            hvg_norm_mode="libnorm1e4_log1p",
        )
        out_1e4 = ad_1e4.predict_deltas(
            target_line="nadig_hepg2", control_line_label="synth_A",
            ctrl_cells=raw_a, pert_names=perts, gene_names=eval_genes,
            ctrl_raw_gene_names=g_model, allow_synthetic_without_freeze=True,
        )
        self.assertEqual(out_1e4.shape, (3, 20))
        self.assertTrue(np.isfinite(out_1e4).all())


if __name__ == "__main__":
    unittest.main()
