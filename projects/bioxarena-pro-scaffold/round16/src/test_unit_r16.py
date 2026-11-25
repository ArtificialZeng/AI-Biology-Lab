#!/usr/bin/env python3
"""
Unit & Invariance Test Suite for Round 16 (`PREREG_I v5 + Addendum A`, incorporating `I14`,
`M1–M6`, `Q4`, Gate-2 locked constants, triple execution gate, pre-download `GSE295214` reader,
and Addendum A Items 1–11).

Runs strictly on 100% synthetic fixtures (touches ZERO real perturbation datasets,
downloads ZERO weights, opens ZERO archived `adata_pred.h5ad` files).
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import inspect
import json
import os
import subprocess
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
    EXACT_INVARIANCE_TOL,
    RESAMPLE_EXCEEDS_WORDING,
    RESAMPLE_NOT_EXCEEDS_WORDING,
    SUB_DELTA_INTERPRETATION_NOTE,
    adjudicate_claim_i1_line,
    adjudicate_claim_i2_line,
    adjudicate_group_r16,
    classify_pair_tost_r16,
    summarize_within_line_resample_r16,
)
from metrics_r12 import compute_joint_sa8_split_minus_all_r12  # noqa: E402
from metrics_r16 import (  # noqa: E402
    DEFAULT_R16_SEED,
    DEFAULT_RESAMPLE_SEED,
    compute_gene_label_perm_null_r16,
    compute_joint_cross_mask_r16,
    compute_orth_suite_r16,
    compute_paired_contrast_r16,
    compute_paired_model_vs_all_mismatched_r16,
    compute_paired_with_within_line_resample_r16,
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
    encode_controls_se600m_r16,
    preprocess_hvg_counts_r16,
    sample_control_indices,
    sample_disjoint_resample_control_sets_r16,
    sample_stratified_control_sets_r16,
)
from run_eval_r16 import (  # noqa: E402
    ARC_ADAPTER_CHECK_LINE,
    _extract_a549_ntc_bundle_r16,
    _extract_feng_ipsc_ntc_bundle_r16,
    _resolve_nadig_h5ad_path_r16,
    build_gse295214_parts,
    compare_predicted_deltas_vs_arc_archive_r16,
    evaluate_line_from_arrays_r16,
    extract_gse295214_ntc_raw_counts,
    extract_nadig_h5ad_ntc_bundle_r16,
    parse_gse295214_guide_r16,
    run_real_evaluation_r16,
    run_synthetic_smoke_r16,
    save_control_embeddings_r16,
    verify_arc_archived_adata_pred_r16,
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

        # Exact integer Hamilton largest-remainder tie-breaking (3 equal 2/3 remainders -> lexicographic B0, B1)
        tie_batches = ["B0"] * 300 + ["B1"] * 300 + ["B2"] * 300
        tie_sets = sample_stratified_control_sets_r16(900, batch_labels=tie_batches, seed=DEFAULT_R16_SEED)
        tie_sampled = [tie_batches[i] for i in tie_sets.ravel()]
        self.assertEqual(tie_sampled.count("B0"), 171)
        self.assertEqual(tie_sampled.count("B1"), 171)
        self.assertEqual(tie_sampled.count("B2"), 170)

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
        hgnc = _StubHGNC({"TP53", "MYC", "HLA-A", "BRCA1", "P2RX7"})
        self.assertEqual(parse_gse295214_guide_r16("non-targeting_0001", "non-targeting_0001", hgnc)[0], "ntc")
        self.assertEqual(parse_gse295214_guide_r16("safe_targeting_1", "safe_targeting_1", hgnc)[0], "exclude")
        self.assertEqual(parse_gse295214_guide_r16("intergenic_001", "intergenic_001", hgnc)[0], "exclude")
        self.assertEqual(parse_gse295214_guide_r16("TP53-and-MYC", "TP53-and-MYC", hgnc)[0], "exclude")
        self.assertEqual(parse_gse295214_guide_r16("TP53_MYC_1", "TP53_MYC_1", hgnc)[0], "exclude")
        self.assertEqual(parse_gse295214_guide_r16("TP53+MYC", "TP53+MYC", hgnc)[0], "exclude")
        self.assertEqual(parse_gse295214_guide_r16("TP53_1+MYC_1", "TP53_1+MYC_1", hgnc)[0], "exclude")
        self.assertEqual(parse_gse295214_guide_r16("TP53_+_MYC", "TP53_+_MYC", hgnc)[0], "exclude")
        self.assertEqual(parse_gse295214_guide_r16("TP53_P2RX7_1", "TP53_P2RX7_1", hgnc)[0], "exclude")
        self.assertEqual(parse_gse295214_guide_r16("TP53-P2RX7", "TP53-P2RX7", hgnc)[0], "exclude")
        wp_kind, wp_info = parse_gse295214_guide_r16("TP53_+_7687490.23-P1P2", "TP53_+_7687490.23-P1P2", hgnc)
        self.assertEqual(wp_kind, "target")
        self.assertEqual(wp_info[1], "TP53")
        wm_kind, wm_info = parse_gse295214_guide_r16("TP53_-_7687490.23-P1P2", "TP53_-_7687490.23-P1P2", hgnc)
        self.assertEqual(wm_kind, "target")
        self.assertEqual(wm_info[1], "TP53")
        sc_kind, sc_info = parse_gse295214_guide_r16("TP53_+12345_P1", "TP53_+12345_P1", hgnc)
        self.assertEqual(sc_kind, "target")
        self.assertEqual(sc_info[1], "TP53")
        p2_kind, p2_info = parse_gse295214_guide_r16("P2RX7_P1P2_1", "P2RX7_P1P2_1", hgnc)
        self.assertEqual(p2_kind, "target")
        self.assertEqual(p2_info[1], "P2RX7")
        tss_h_kind, tss_h_info = parse_gse295214_guide_r16("TP53-TSS1", "TP53-TSS1", hgnc)
        self.assertEqual(tss_h_kind, "target")
        self.assertEqual(tss_h_info[1], "TP53")
        tss_u_kind, tss_u_info = parse_gse295214_guide_r16("MYC_TSS1_1", "MYC_TSS1_1", hgnc)
        self.assertEqual(tss_u_kind, "target")
        self.assertEqual(tss_u_info[1], "MYC")
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
            self.assertEqual(len(c_aud["sampled_barcodes"]), 128)
            self.assertEqual(len(set(c_aud["sampled_barcodes"])), 128)
            self.assertTrue(all(bc.startswith("GSM89437") and ":" in bc for bc in c_aud["sampled_barcodes"]))

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
            self.assertIn("b17_vs_b0_sym_mis_nadig_hepg2", line_r["secondary"]["g4_factorial_decomposition"])
            self.assertEqual(len(line_r["secondary"]["g4_factorial_decomposition"]), 10)
            g4_st = line_r["secondary"]["g4_stratum_summary"]
            self.assertIn("b0_vs_b17", g4_st)
            self.assertIn("b17_vs_b0_sym", g4_st)
            self.assertEqual(g4_st["b0_vs_b17"]["seen_by_checkpoint"]["n_sources"], 2)
            self.assertEqual(g4_st["b0_vs_b17"]["unseen_by_checkpoint"]["n_sources"], 3)
            self.assertEqual(
                set(line_r["secondary"]["g4_single_index_spread"]["per_index_r_orth_cal"].keys()),
                {"b0", "b1", "b17", "b35"},
            )
            self.assertEqual(line_r["secondary"]["STATE_HVG"]["perm100_r_orth_cal"]["n_finite_perms"], 100)
            self.assertIsNotNone(line_r["secondary"]["within_line_resample_reference"])
            self.assertIn("bypass_zero", line_r["secondary"]["neutral_batch_robustness_panel"])
            self.assertTrue((out_dir / "control_barcodes_synthetic_smoke.json").is_file())

            ldir = out_dir / "arrays" / "synthetic_a549_ko"
            manifest = json.loads((ldir / "array_manifest_sha256.json").read_text(encoding="utf-8"))
            self.assertIn("observed.npz", manifest)
            self.assertIn("preds_STATE_SE_fixed_index_zero_matched.npz", manifest)
            self.assertIn("preds_STATE_SE_fixed_index_zero_resample.npz", manifest)
            g4_s_npz = np.load(ldir / "preds_STATE_SE_G4_b0_vs_b17_mis_nadig_hepg2_S_BM.npz", allow_pickle=False)
            self.assertEqual(str(g4_s_npz["control_source"]), "nadig_hepg2")
            self.assertEqual(len(g4_s_npz["control_barcodes"]), 512)
            self.assertTrue(all(str(x).startswith("SYN_NADIG_HEPG2_") for x in g4_s_npz["control_barcodes"]))

            # Independently reload exported `.npz` arrays and recompute `H-I1` & `H-I2` to <= 1e-12
            obs_npz = np.load(ldir / "observed.npz", allow_pickle=False)
            self.assertEqual(len(obs_npz["matched_ctrl_barcodes"]), 512)
            self.assertEqual(len(obs_npz["resample_ctrl_barcodes"]), 512)
            self.assertEqual(
                len(set(obs_npz["matched_ctrl_barcodes"].tolist()) & set(obs_npz["resample_ctrl_barcodes"].tolist())),
                0,
            )
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

            # Addendum A Item 11: verify exported `X_state` embeddings + SHA-256 manifest
            emb_dir = out_dir / "embeddings"
            emb_manifest_path = emb_dir / "embeddings_manifest_sha256.json"
            self.assertTrue(emb_manifest_path.is_file())
            emb_manifest = json.loads(emb_manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(emb_manifest, res["embeddings_manifest_sha256"])
            expected_emb_files = {
                "x_state_synthetic_a549_ko.npz",
                "x_state_synthetic_a549_ko_resample.npz",
                *(f"x_state_{sk}.npz" for sk in all_mis),
            }
            self.assertEqual(set(emb_manifest.keys()), expected_emb_files)
            for efname, expected_sha in sorted(emb_manifest.items()):
                epath = emb_dir / efname
                self.assertTrue(epath.is_file())
                self.assertEqual(hashlib.sha256(epath.read_bytes()).hexdigest(), expected_sha)
                ez = np.load(epath, allow_pickle=False)
                self.assertEqual(ez["X_state"].shape, (512, 2058))
                self.assertEqual(ez["X_state"].dtype, np.float32)
                self.assertTrue(np.all(np.isfinite(ez["X_state"])))
                self.assertEqual(len(ez["control_barcodes"]), 512)
                expected_pool_key = efname[len("x_state_") : -len(".npz")]
                self.assertEqual(str(ez["dataset_name"]), expected_pool_key)
                if efname == "x_state_synthetic_a549_ko.npz":
                    self.assertEqual(
                        ez["control_barcodes"].tolist(), obs_npz["matched_ctrl_barcodes"].tolist()
                    )
                elif efname == "x_state_synthetic_a549_ko_resample.npz":
                    self.assertEqual(
                        ez["control_barcodes"].tolist(), obs_npz["resample_ctrl_barcodes"].tolist()
                    )
                expected_seed = (
                    DEFAULT_RESAMPLE_SEED
                    if efname == "x_state_synthetic_a549_ko_resample.npz"
                    else DEFAULT_R16_SEED
                )
                self.assertEqual(int(ez["seed"]), expected_seed)
                aud_obj = json.loads(str(ez["se_audit_json"]))
                self.assertEqual(aud_obj["seed"], expected_seed)
                self.assertEqual(aud_obj["batch_size"], 16)
                self.assertFalse(aud_obj["shuffle"])
                self.assertEqual(aud_obj["num_workers"], 0)
                self.assertFalse(aud_obj["collator_training"])
                self.assertEqual(
                    aud_obj["x_state_sha256"],
                    hashlib.sha256(np.ascontiguousarray(ez["X_state"], dtype=np.float32).tobytes()).hexdigest(),
                )

            # Edge cases on `save_control_embeddings_r16`, `encode_controls_se600m_r16`, and `ArcStateAdapter`
            none_opt_dir = out_dir / "emb_none_opt"
            none_opt_man = save_control_embeddings_r16(
                none_opt_dir,
                {
                    "opt_pool": {
                        "X_state": np.zeros((4, 2058), dtype=np.float32),
                        "control_barcodes": None,
                        "seed": DEFAULT_R16_SEED,
                        "se_audit": None,
                    }
                },
            )
            self.assertEqual(set(none_opt_man.keys()), {"x_state_opt_pool.npz"})
            ez_opt = np.load(none_opt_dir / "x_state_opt_pool.npz", allow_pickle=False)
            self.assertEqual(len(ez_opt["control_barcodes"]), 0)
            self.assertEqual(json.loads(str(ez_opt["se_audit_json"])), {})

            with self.assertRaises(ValueError):
                save_control_embeddings_r16(
                    out_dir / "bad_emb_dim",
                    {"bad": {"X_state": np.zeros((4, 128), dtype=np.float32), "seed": DEFAULT_R16_SEED}},
                )
            with self.assertRaises(ValueError):
                save_control_embeddings_r16(
                    out_dir / "bad_emb_empty",
                    {"bad": {"X_state": np.zeros((0, 2058), dtype=np.float32), "seed": DEFAULT_R16_SEED}},
                )
            bad_nan = np.zeros((4, 2058), dtype=np.float32)
            bad_nan[0, 0] = np.nan
            with self.assertRaises(ValueError):
                save_control_embeddings_r16(
                    out_dir / "bad_emb_nan",
                    {"bad": {"X_state": bad_nan, "seed": DEFAULT_R16_SEED}},
                )
            with self.assertRaises(ValueError):
                save_control_embeddings_r16(
                    out_dir / "bad_emb_bc",
                    {
                        "bad": {
                            "X_state": np.zeros((4, 2058), dtype=np.float32),
                            "control_barcodes": ["B0", "B1", "B2"],
                            "seed": DEFAULT_R16_SEED,
                        }
                    },
                )
            with self.assertRaises(ValueError):
                save_control_embeddings_r16(
                    out_dir / "bad_emb_dup_fname",
                    {
                        "pool/a": {"X_state": np.zeros((4, 2058), dtype=np.float32), "seed": DEFAULT_R16_SEED},
                        "pool a": {"X_state": np.zeros((4, 2058), dtype=np.float32), "seed": DEFAULT_R16_SEED},
                    },
                )
            for bad_key in ("", "   ", "..."):
                with self.assertRaises(ValueError):
                    save_control_embeddings_r16(
                        out_dir / "bad_emb_key",
                        {bad_key: {"X_state": np.zeros((4, 2058), dtype=np.float32), "seed": DEFAULT_R16_SEED}},
                    )
            for bad_seed in (None, True, False):
                with self.assertRaises(ValueError):
                    save_control_embeddings_r16(
                        out_dir / "bad_emb_seed",
                        {"p": {"X_state": np.zeros((4, 2058), dtype=np.float32), "seed": bad_seed}},
                    )
            with self.assertRaises(ValueError):
                encode_controls_se600m_r16(np.zeros((0, 5), dtype=np.float32), ["A", "B", "C", "D", "E"])
            with self.assertRaises(ValueError):
                encode_controls_se600m_r16(np.zeros((2, 5), dtype=np.float32), ["A", "B"])
            with self.assertRaises(ValueError):
                encode_controls_se600m_r16(np.zeros((2, 2), dtype=np.float32), ["A", "B"], batch_size=0)
            with self.assertRaises(ValueError):
                encode_controls_se600m_r16(np.zeros((2, 2), dtype=np.float32), ["A", "B"], batch_size=True)
            with self.assertRaises(ValueError):
                encode_controls_se600m_r16(np.zeros((2, 2), dtype=np.float32), ["A", "B"], seed=True)
            with self.assertRaises(ValueError):
                encode_controls_se600m_r16(np.array([[1.0, np.nan]], dtype=np.float32), ["A", "B"])
            sig_pred = inspect.signature(ArcStateAdapter._predict_impl)
            self.assertIn("se_encode_seed", sig_pred.parameters)
            self.assertEqual(sig_pred.parameters["se_encode_seed"].default, DEFAULT_R16_SEED)

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

    def test_13_addendum_a_item10_exact_invariance_vs_zero_tie(self) -> None:
        """Addendum A Item 10: `max_abs_diff <= 1e-12` -> `EXACT_INVARIANCE` (counted as EQUIVALENT); non-identical [0,0] stays INCONCLUSIVE."""
        rng = np.random.default_rng(1301)
        P, G = 70, 200
        tm = rng.normal(0.0, 0.2, size=(P, G))
        true_spec = rng.normal(0.0, 0.3, size=(P, G))
        ha = tm + true_spec + rng.normal(0.0, 0.04, size=(P, G))
        hb = tm + true_spec + rng.normal(0.0, 0.04, size=(P, G))
        obs = 0.5 * (ha + hb)
        t_idx = np.arange(P, dtype=np.int64)
        hv = np.ones(P, dtype=bool)
        seen = ("nadig_hepg2", "nadig_jurkat")
        unseen = ("cd4t_rest", "dld1", "feng_ipsc")
        all_mis = seen + unseen

        # Case 1: Positive matched arm + identical mismatched arrays (exact 0.0 and 5e-13 <= 1e-12)
        pred_pos = tm + 0.7 * true_spec + rng.normal(0.0, 0.04, size=(P, G))
        mis_exact = {sk: pred_pos.copy() for sk in all_mis}
        mis_exact["dld1"] = pred_pos + 5e-13

        orth_pos = compute_orth_suite_r16(
            obs, pred_pos, tm, t_idx, source_grand_mean=None,
            delta_half_a=ha, delta_half_b=hb, half_valid=hv, n_boot=250, seed=DEFAULT_R16_SEED,
        )
        paired_exact = compute_paired_model_vs_all_mismatched_r16(
            obs, pred_pos, mis_exact, tm, t_idx, source_grand_mean=None,
            delta_half_a=ha, delta_half_b=hb, half_valid=hv, n_boot=250, seed=DEFAULT_R16_SEED,
        )
        for sk in all_mis:
            self.assertTrue(paired_exact["per_source"][sk]["exact_invariance"])
            self.assertLessEqual(paired_exact["per_source"][sk]["max_abs_diff"], EXACT_INVARIANCE_TOL)
        adj_pos = adjudicate_claim_i1_line(
            orth_pos, paired_exact, seen_sources=seen, unseen_sources=unseen, delta=DEFAULT_TOST_DELTA
        )
        self.assertEqual(adj_pos["line_summary"], "LINE_EQUIVALENT_NONZERO_SIGNAL")
        self.assertEqual(adj_pos["counts"]["EXACT_INVARIANCE"], 5)
        self.assertEqual(adj_pos["counts"]["INCONCLUSIVE"], 0)

        # Case 2: Null-signal matched arm + identical mismatched arrays -> LINE_EQUIVALENT_NULL_SIGNAL
        pred_null = tm - 0.15 * true_spec + rng.normal(0.0, 0.25, size=(P, G))
        orth_null = compute_orth_suite_r16(
            obs, pred_null, tm, t_idx, source_grand_mean=None,
            delta_half_a=ha, delta_half_b=hb, half_valid=hv, n_boot=250, seed=DEFAULT_R16_SEED,
        )
        self.assertLessEqual(orth_null["ci975_lo"], 0.0)
        paired_null = compute_paired_model_vs_all_mismatched_r16(
            obs, pred_null, {sk: pred_null.copy() for sk in all_mis}, tm, t_idx,
            source_grand_mean=None, delta_half_a=ha, delta_half_b=hb, half_valid=hv,
            n_boot=250, seed=DEFAULT_R16_SEED,
        )
        adj_null = adjudicate_claim_i1_line(
            orth_null, paired_null, seen_sources=seen, unseen_sources=unseen, delta=DEFAULT_TOST_DELTA
        )
        self.assertEqual(adj_null["line_summary"], "LINE_EQUIVALENT_NULL_SIGNAL")
        self.assertEqual(adj_null["counts"]["EXACT_INVARIANCE"], 5)

        # Case 3: Mix of EXACT_INVARIANCE and EQUIVALENT counts cleanly as EQUIVALENT at line level
        eq_pair = {
            "estimate": 0.002, "ci975_lo": -0.015, "ci975_hi": 0.019,
            "ci95_lo": -0.012, "ci95_hi": 0.016, "valid_frac": 1.0,
            "max_abs_diff": 0.10, "exact_invariance": False,
        }
        exact_pair = {
            "estimate": 0.0, "ci975_lo": 0.0, "ci975_hi": 0.0,
            "ci95_lo": 0.0, "ci95_hi": 0.0, "valid_frac": 1.0,
            "max_abs_diff": 0.0, "exact_invariance": True,
        }
        adj_mix = adjudicate_claim_i1_line(
            {"ci975_lo": 0.10, "ci975_hi": 0.30, "valid_frac": 1.0, "rho_ab_orth": 0.7},
            {"per_source": {**{k: exact_pair for k in seen}, **{k: eq_pair for k in unseen}}},
            seen_sources=seen, unseen_sources=unseen, delta=DEFAULT_TOST_DELTA,
        )
        self.assertEqual(adj_mix["line_summary"], "LINE_EQUIVALENT_NONZERO_SIGNAL")
        self.assertEqual(adj_mix["counts"]["EXACT_INVARIANCE"], 2)
        self.assertEqual(adj_mix["counts"]["EQUIVALENT"], 3)

        # Case 4: Non-identical arrays (`max_abs_diff = 1e-6 > 1e-12`) with [0, 0] CI MUST stay INCONCLUSIVE!
        non_id_zero = {
            "estimate": 0.0, "ci975_lo": 0.0, "ci975_hi": 0.0,
            "ci95_lo": 0.0, "ci95_hi": 0.0, "valid_frac": 1.0,
            "max_abs_diff": 1e-6, "exact_invariance": False,
        }
        c_non_id = classify_pair_tost_r16(non_id_zero, DEFAULT_TOST_DELTA)
        self.assertEqual(c_non_id["category"], "INCONCLUSIVE")
        self.assertFalse(c_non_id["exact_invariance"])

        # Case 5: Low valid_frac or rho_ab_orth <= 0 still forces INCONCLUSIVE even when max_abs_diff == 0
        c_low_vf = classify_pair_tost_r16({**exact_pair, "valid_frac": 0.50}, DEFAULT_TOST_DELTA)
        self.assertEqual(c_low_vf["category"], "INCONCLUSIVE")
        adj_bad_rho = adjudicate_claim_i1_line(
            {"ci975_lo": 0.10, "ci975_hi": 0.30, "valid_frac": 1.0, "rho_ab_orth": 0.0},
            {"per_source": {k: exact_pair for k in all_mis}},
            seen_sources=seen, unseen_sources=unseen, delta=DEFAULT_TOST_DELTA,
        )
        self.assertEqual(adj_bad_rho["line_summary"], "INCONCLUSIVE")

    def test_14_addendum_a_item8_item9_within_line_resample_and_sub_delta(self) -> None:
        """Addendum A Items 8 & 9: disjoint `seed=20261602` resample sampler, paired excess CI, & locked wording."""
        # 1. Disjoint resample sampler boundaries (`>=1024`, `1024`, `640..1023`, `<640`)
        b1500 = ["B0"] * 750 + ["B1"] * 450 + ["B2"] * 300
        p1500 = sample_stratified_control_sets_r16(1500, batch_labels=b1500, seed=DEFAULT_R16_SEED)
        r1500_a, aud1500 = sample_disjoint_resample_control_sets_r16(
            1500, p1500, batch_labels=b1500, seed=DEFAULT_RESAMPLE_SEED
        )
        r1500_b, _ = sample_disjoint_resample_control_sets_r16(
            1500, p1500, batch_labels=b1500, seed=DEFAULT_RESAMPLE_SEED
        )
        self.assertIsNotNone(r1500_a)
        self.assertTrue(np.array_equal(r1500_a, r1500_b))
        self.assertEqual(r1500_a.shape, (16, 32))
        self.assertEqual(aud1500["overlap_with_primary"], 0)
        self.assertEqual(len(set(p1500.ravel().tolist()) & set(r1500_a.ravel().tolist())), 0)
        self.assertEqual(aud1500["regime"], "full_512_disjoint")

        p1024 = sample_stratified_control_sets_r16(1024, seed=DEFAULT_R16_SEED)
        r1024, aud1024 = sample_disjoint_resample_control_sets_r16(1024, p1024, seed=DEFAULT_RESAMPLE_SEED)
        self.assertEqual(r1024.shape, (16, 32))
        self.assertEqual(len(set(p1024.ravel().tolist()) & set(r1024.ravel().tolist())), 0)
        self.assertEqual(aud1024["regime"], "full_512_disjoint")

        p700 = sample_stratified_control_sets_r16(700, seed=DEFAULT_R16_SEED)
        r700, aud700 = sample_disjoint_resample_control_sets_r16(700, p700, seed=DEFAULT_RESAMPLE_SEED)
        self.assertEqual(r700.shape, (5, 32))  # (700 - 512) // 32 * 32 = 160
        self.assertEqual(len(set(p700.ravel().tolist()) & set(r700.ravel().tolist())), 0)
        self.assertEqual(aud700["regime"], "truncated_disjoint_complement")

        p640 = sample_stratified_control_sets_r16(640, seed=DEFAULT_R16_SEED)
        r640, aud640 = sample_disjoint_resample_control_sets_r16(640, p640, seed=DEFAULT_RESAMPLE_SEED)
        self.assertEqual(r640.shape, (4, 32))
        self.assertTrue(aud640["eligible"])

        p639 = sample_stratified_control_sets_r16(639, seed=DEFAULT_R16_SEED)
        r639, aud639 = sample_disjoint_resample_control_sets_r16(639, p639, seed=DEFAULT_RESAMPLE_SEED)
        self.assertIsNone(r639)
        self.assertFalse(aud639["eligible"])

        # 2. Controller synthetic SUB_DELTA noise scenario: isotropic noise on mismatched arms == resample noise
        rng = np.random.default_rng(1402)
        P, G = 200, 200
        tm = rng.normal(0.0, 0.20, size=(P, G))
        true_spec = rng.normal(0.0, 0.30, size=(P, G))
        ha = tm + true_spec + rng.normal(0.0, 0.03, size=(P, G))
        hb = tm + true_spec + rng.normal(0.0, 0.03, size=(P, G))
        obs = 0.5 * (ha + hb)
        t_idx = np.arange(P, dtype=np.int64)
        hv = np.ones(P, dtype=bool)
        seen = ("nadig_hepg2", "nadig_jurkat")
        unseen = ("cd4t_rest", "dld1", "feng_ipsc")
        all_mis = seen + unseen

        pred_m = tm + 0.70 * true_spec + rng.normal(0.0, 0.06, size=(P, G))
        pred_resamp_noisy = pred_m + rng.normal(0.0, 0.045, size=(P, G))
        mis_noisy = {sk: pred_m + rng.normal(0.0, 0.035, size=(P, G)) for sk in all_mis}

        orth_m = compute_orth_suite_r16(
            obs, pred_m, tm, t_idx, delta_half_a=ha, delta_half_b=hb, half_valid=hv,
            n_boot=300, seed=DEFAULT_R16_SEED,
        )
        paired_noisy = compute_paired_model_vs_all_mismatched_r16(
            obs, pred_m, mis_noisy, tm, t_idx, delta_half_a=ha, delta_half_b=hb, half_valid=hv,
            n_boot=300, seed=DEFAULT_R16_SEED,
        )
        adj_sub = adjudicate_claim_i1_line(
            orth_m, paired_noisy, seen_sources=seen, unseen_sources=unseen, delta=DEFAULT_TOST_DELTA
        )
        self.assertEqual(adj_sub["line_summary"], "LINE_MATCHED_SPECIFIC_SUB_DELTA")
        self.assertEqual(adj_sub["interpretation_note"], SUB_DELTA_INTERPRETATION_NOTE)

        resamp_cmp_noisy = compute_paired_with_within_line_resample_r16(
            obs, pred_m, pred_resamp_noisy, mis_noisy, tm, t_idx,
            delta_half_a=ha, delta_half_b=hb, half_valid=hv, n_boot=300, seed=DEFAULT_R16_SEED,
        )
        interp_noisy = summarize_within_line_resample_r16(
            resamp_cmp_noisy, seen_sources=seen, unseen_sources=unseen,
            matched_positive_ci975=adj_sub["matched_positive_ci975"],
        )
        self.assertTrue(interp_noisy["none_exceeds_resample_loss"])
        self.assertEqual(interp_noisy["interpretation_zh"], RESAMPLE_NOT_EXCEEDS_WORDING)
        for sk in all_mis:
            self.assertIn("excess_ci95_lo", interp_noisy["per_source"][sk])
            self.assertIn("excess_ci95_hi", interp_noisy["per_source"][sk])

        # 3. Genuine target-specific control usage: mismatched loss >> within-line resample loss
        mis_true_loss = {sk: tm + 0.10 * true_spec + rng.normal(0.0, 0.20, size=(P, G)) for sk in all_mis}
        resamp_cmp_real = compute_paired_with_within_line_resample_r16(
            obs, pred_m, pred_resamp_noisy, mis_true_loss, tm, t_idx,
            delta_half_a=ha, delta_half_b=hb, half_valid=hv, n_boot=300, seed=DEFAULT_R16_SEED,
        )
        interp_real = summarize_within_line_resample_r16(
            resamp_cmp_real, seen_sources=seen, unseen_sources=unseen, matched_positive_ci975=True
        )
        self.assertTrue(interp_real["all_exceeds_resample_loss"])
        self.assertEqual(interp_real["interpretation_zh"], RESAMPLE_EXCEEDS_WORDING)

        # 4. Edge-case guards: all_exceed with matched_positive=False & invalid/low-valid_frac CI
        interp_null_m = summarize_within_line_resample_r16(
            resamp_cmp_real, seen_sources=seen, unseen_sources=unseen, matched_positive_ci975=False
        )
        self.assertEqual(interp_null_m["resample_comparison_status"], "EXCEEDS_RESAMPLE_BUT_MATCHED_NOT_POSITIVE")
        bad_resamp_dict = {
            "delta_pair_resample": {},
            "per_source": {
                sk: {
                    "delta_pair_excess_over_resample": {
                        "estimate": float("nan"), "ci975_lo": float("nan"), "ci975_hi": float("nan"), "valid_frac": 0.50
                    }
                }
                for sk in all_mis
            },
        }
        interp_inv = summarize_within_line_resample_r16(
            bad_resamp_dict, seen_sources=seen, unseen_sources=unseen, matched_positive_ci975=True
        )
        self.assertEqual(interp_inv["resample_comparison_status"], "INCONCLUSIVE_VS_WITHIN_LINE_RESAMPLE")
        self.assertFalse(interp_inv["none_exceeds_resample_loss"])

    def test_15_addendum_a_item4_sthvg_perm100_and_item5_cellularai_dropped(self) -> None:
        """Addendum A Item 4 (`ST-HVG` `Perm100`, `k=0..99`, `default_rng(k)`) & Item 5 (`CellularAI` dropped)."""
        rng = np.random.default_rng(1501)
        P, G = 45, 160
        tm = rng.normal(0.0, 0.2, size=(P, G))
        true_spec = rng.normal(0.0, 0.3, size=(P, G))
        ha = tm + true_spec + rng.normal(0.0, 0.04, size=(P, G))
        hb = tm + true_spec + rng.normal(0.0, 0.04, size=(P, G))
        obs = 0.5 * (ha + hb)
        pred = tm + 0.65 * true_spec + rng.normal(0.0, 0.06, size=(P, G))
        t_idx = np.arange(P, dtype=np.int64)
        hv = np.ones(P, dtype=bool)

        p100 = compute_gene_label_perm_null_r16(
            obs, pred, tm, t_idx, source_grand_mean=None,
            delta_half_a=ha, delta_half_b=hb, half_valid=hv, n_perm=100, perm_seed_offset=0,
        )
        self.assertEqual(p100["n_perm"], 100)
        self.assertEqual(p100["n_finite_perms"], 100)
        self.assertGreater(p100["observed_estimate"], p100["q975"])
        self.assertLess(abs(p100["median"]), 0.15)

        # Bit-for-bit parity vs direct `compute_orth_suite_r16` on k=0, 17, 99
        for k in (0, 17, 99):
            perm_k = np.random.default_rng(k).permutation(G)
            direct = compute_orth_suite_r16(
                obs, pred[:, perm_k], tm, t_idx, source_grand_mean=None,
                delta_half_a=ha, delta_half_b=hb, half_valid=hv, n_boot=20, seed=DEFAULT_R16_SEED,
            )
            self.assertLessEqual(abs(direct["estimate"] - p100["perm_estimates"][k]), 1e-14)

        # Item 5: CellularAI_OpenCode_PublicRefit explicitly dropped pre-data
        cell_ai = CellularAIOpenCodePublicRefitAdapter()
        self.assertEqual(cell_ai.spec.role, "dropped_per_addendum_a")
        with self.assertRaises(NotImplementedError):
            cell_ai.predict_deltas(
                "a549_ko", "a549_ko", np.ones((8, 5)), ["G1"], ["A", "B", "C", "D", "E"],
                allow_synthetic_without_freeze=True,
            )

    def test_16_addendum_a_item3_arc_adata_pred_check_and_item1_fixed_targets(self) -> None:
        """Addendum A Item 3 (`nadig_hepg2` `adata_pred.h5ad` check >=0.90/>=0.95) & Item 1 fixed-target FCT sensitivity."""
        rng = np.random.default_rng(1601)
        P, G = 35, 120
        perts = [f"GENE_{i:03d}" for i in range(P)]
        genes = [f"OUT_{j:03d}" for j in range(G)]
        base = rng.normal(0.0, 0.4, size=(P, G))
        arc_close = base + rng.normal(0.0, 0.04, size=(P, G))
        arc_bad = 0.4 * base + rng.normal(0.0, 0.50, size=(P, G))

        res_pass = compare_predicted_deltas_vs_arc_archive_r16(base, perts, genes, arc_close, perts, genes)
        self.assertTrue(res_pass["passed"])
        self.assertFalse(res_pass["stop_for_controller_diagnosis"])
        self.assertGreaterEqual(res_pass["median_per_pert_pearson"], 0.90)
        self.assertGreaterEqual(res_pass["pooled_pearson"], 0.95)

        res_fail = compare_predicted_deltas_vs_arc_archive_r16(base, perts, genes, arc_bad, perts, genes)
        self.assertFalse(res_fail["passed"])
        self.assertTrue(res_fail["stop_for_controller_diagnosis"])

        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            h5p = tdir / "adata_pred.h5ad"
            # Write synthetic archived adata_pred.h5ad (20 NTC + 4 cells per pert)
            ntc_base = rng.uniform(1.0, 3.0, size=(1, G))
            rows_x = [ntc_base + rng.normal(0.0, 0.01, size=(20, G))]
            obs_perts = ["non-targeting"] * 20
            for i, p_name in enumerate(perts):
                rows_x.append(ntc_base + arc_close[i : i + 1, :] + rng.normal(0.0, 0.01, size=(4, G)))
                obs_perts.extend([p_name] * 4)
            X_full = np.vstack(rows_x).astype(np.float32)
            with h5py.File(h5p, "w") as f:
                f.create_dataset("X", data=X_full)
                og = f.create_group("obs")
                og.create_dataset("gene", data=np.array(obs_perts, dtype=h5py.string_dtype("utf-8")))
                vg = f.create_group("var")
                vg.create_dataset("gene_name", data=np.array(genes, dtype=h5py.string_dtype("utf-8")))

            # Pre-gate MUST fail with BlindingGateError
            with self.assertRaises(BlindingGateError):
                verify_arc_archived_adata_pred_r16(h5p, base, perts, genes, status_dir=tdir)

            # Unlock mock gates; non-nadig_hepg2 MUST fail with ValueError
            for gfn in ("PREREG_I_COMMIT.txt", "CODE_FREEZE_R16.txt", "R16_GO.txt"):
                (tdir / gfn).write_text("OK\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                verify_arc_archived_adata_pred_r16(
                    h5p, base, perts, genes, line_name="nadig_jurkat", status_dir=tdir
                )
            chk = verify_arc_archived_adata_pred_r16(
                h5p, base, perts, genes, line_name=ARC_ADAPTER_CHECK_LINE, status_dir=tdir
            )
            self.assertTrue(chk["passed"])
            self.assertEqual(chk["n_common_perts"], P)
            self.assertEqual(chk["n_common_genes"], G)

            # Verify CSC-encoded archived adata_pred.h5ad gives exact same result (hoisted CSC->CSR path)
            h5p_csc = tdir / "adata_pred_csc.h5ad"
            X_csc = sp.csc_matrix(X_full)
            with h5py.File(h5p_csc, "w") as f:
                xg = f.create_group("X")
                xg.attrs["encoding-type"] = "csc_matrix"
                xg.attrs["shape"] = np.array(X_csc.shape, dtype=np.int64)
                xg.create_dataset("data", data=X_csc.data)
                xg.create_dataset("indices", data=X_csc.indices)
                xg.create_dataset("indptr", data=X_csc.indptr)
                og = f.create_group("obs")
                og.create_dataset("gene", data=np.array(obs_perts, dtype=h5py.string_dtype("utf-8")))
                vg = f.create_group("var")
                vg.create_dataset("gene_name", data=np.array(genes, dtype=h5py.string_dtype("utf-8")))
            chk_csc = verify_arc_archived_adata_pred_r16(
                h5p_csc, base, perts, genes, line_name=ARC_ADAPTER_CHECK_LINE, status_dir=tdir
            )
            self.assertTrue(chk_csc["passed"])
            self.assertAlmostEqual(chk_csc["pooled_pearson"], chk["pooled_pearson"], places=10)
            self.assertAlmostEqual(chk_csc["median_per_pert_pearson"], chk["median_per_pert_pearson"], places=10)

            # Verify CSR-encoded archived adata_pred.h5ad gives exact same result (single-pass CSR chunk path)
            h5p_csr = tdir / "adata_pred_csr.h5ad"
            X_csr = sp.csr_matrix(X_full)
            with h5py.File(h5p_csr, "w") as f:
                xg = f.create_group("X")
                xg.attrs["encoding-type"] = "csr_matrix"
                xg.attrs["shape"] = np.array(X_csr.shape, dtype=np.int64)
                xg.create_dataset("data", data=X_csr.data)
                xg.create_dataset("indices", data=X_csr.indices)
                xg.create_dataset("indptr", data=X_csr.indptr)
                og = f.create_group("obs")
                og.create_dataset("gene", data=np.array(obs_perts, dtype=h5py.string_dtype("utf-8")))
                vg = f.create_group("var")
                vg.create_dataset("gene_name", data=np.array(genes, dtype=h5py.string_dtype("utf-8")))
            chk_csr = verify_arc_archived_adata_pred_r16(
                h5p_csr, base, perts, genes, line_name=ARC_ADAPTER_CHECK_LINE, status_dir=tdir
            )
            self.assertTrue(chk_csr["passed"])
            self.assertAlmostEqual(chk_csr["pooled_pearson"], chk["pooled_pearson"], places=10)
            self.assertAlmostEqual(chk_csr["median_per_pert_pearson"], chk["median_per_pert_pearson"], places=10)

            # run_real_evaluation_r16 refuses --verify-arc-adata-pred if nadig_hepg2 is omitted from lines_to_run
            with self.assertRaises(ValueError):
                run_real_evaluation_r16(
                    ["nadig_jurkat"], status_dir=tdir, models_dir=tdir / "no_models",
                    verify_arc_adata_pred_path=h5p,
                )

            # Item 1: extra_fct_sensitivity uses its own delta_transfer & source_grand_mean on fixed primary targets
            tm_prim = rng.normal(0.0, 0.2, size=(P, G))
            tm_sens = tm_prim + rng.normal(0.0, 0.15, size=(P, G))
            bs_prim = tm_prim.mean(axis=0)
            bs_sens = tm_sens.mean(axis=0)
            true_spec = rng.normal(0.0, 0.3, size=(P, G))
            ha = tm_prim + true_spec + rng.normal(0.0, 0.04, size=(P, G))
            hb = tm_prim + true_spec + rng.normal(0.0, 0.04, size=(P, G))
            obs = 0.5 * (ha + hb)
            se_m = tm_prim + 0.7 * true_spec + rng.normal(0.0, 0.15, size=(P, G))
            strata = build_mismatched_strata_r16("gse295214_h1", "zeroshot/ensemble4")
            se_mis = {sk: tm_prim + 0.05 * true_spec for sk in strata["seen_by_checkpoint"] + strata["unseen_by_checkpoint"]}
            ev = evaluate_line_from_arrays_r16(
                line_name="syn_gse295214_h1",
                group_key="group3",
                modality="CRISPRi",
                checkpoint_id="zeroshot/ensemble4",
                seen_sources=strata["seen_by_checkpoint"],
                unseen_sources=strata["unseen_by_checkpoint"],
                delta_obs=obs,
                delta_half_a=ha,
                delta_half_b=hb,
                half_valid=np.ones(P, dtype=bool),
                delta_transfer=tm_prim,
                source_grand_mean=bs_prim,
                target_indices=np.arange(P, dtype=np.int64),
                target_names=perts,
                gene_names=genes,
                mu_target=np.ones(G),
                mu_sources=np.ones((5, G)),
                fct_delta_pred=tm_prim + 0.4 * true_spec,
                e1_constmed_delta_pred=tm_prim,
                state_se_matched=se_m,
                state_se_mismatched=se_mis,
                extra_fct_sensitivity={
                    "FCT_full_r9_6_incl_vcc2025_h1": {
                        "fct_delta_pred": tm_sens + 0.45 * true_spec,
                        "delta_transfer": tm_sens,
                        "source_grand_mean": bs_sens,
                    }
                },
                arrays_dir=tdir / "arr",
                n_boot=100,
                seed=DEFAULT_R16_SEED,
            )
            s_out = ev["secondary"]["extra_fct_sensitivity"]["FCT_full_r9_6_incl_vcc2025_h1"]
            self.assertEqual(s_out["fixed_to_primary_n_eval"], P)
            direct_se_sens = compute_orth_suite_r16(
                obs, se_m, tm_sens, np.arange(P, dtype=np.int64), source_grand_mean=None,
                delta_half_a=ha, delta_half_b=hb, half_valid=np.ones(P, dtype=bool),
                n_boot=100, seed=DEFAULT_R16_SEED,
            )
            self.assertLessEqual(abs(s_out["state_se_orth"]["estimate"] - direct_se_sens["estimate"]), 1e-12)
            self.assertIn("state_se_dorth", s_out)
            self.assertNotAlmostEqual(
                s_out["state_se_orth"]["estimate"],
                ev["primary"]["STATE_SE_matched_orth"]["estimate"],
                places=4,
            )
            self.assertIn("preds_SENS_FCT_full_r9_6_incl_vcc2025_h1.npz", ev["array_manifest_sha256"])

    def test_17_addendum_a_item7_barcodes_nadig_a549_feng_and_nadig_dir_override(self) -> None:
        """Addendum A Item 7 & 8: barcode + disjoint resample extractors (`nadig`, `a549`, `feng`) & `NADIG_DATA_DIR` override."""
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)

            # 1. `--raw-nadig-dir` and `NADIG_DATA_DIR` resolution + synthetic Nadig `.h5ad`
            ndir = tdir / "nadig_local"
            ndir.mkdir()
            hepg2_p = ndir / "GSE264667_hepg2_raw_singlecell_01.h5ad"
            n_ntc, n_tgt, n_g = 1100, 40, 20
            X_n = sp.csr_matrix(np.ones((n_ntc + n_tgt, n_g), dtype=np.float32))
            with h5py.File(hepg2_p, "w") as f:
                xg = f.create_group("X")
                xg.attrs["encoding-type"] = "csr_matrix"
                xg.attrs["shape"] = np.array(X_n.shape, dtype=np.int64)
                xg.create_dataset("data", data=X_n.data)
                xg.create_dataset("indices", data=X_n.indices)
                xg.create_dataset("indptr", data=X_n.indptr)
                og = f.create_group("obs")
                og.attrs["_index"] = "_index"
                og.create_dataset(
                    "_index",
                    data=np.array([f"NADIG_BC_{i:05d}-1" for i in range(n_ntc + n_tgt)], dtype=h5py.string_dtype("utf-8")),
                )
                og.create_dataset(
                    "gene",
                    data=np.array(["non-targeting"] * n_ntc + ["TP53"] * n_tgt, dtype=h5py.string_dtype("utf-8")),
                )
                og.create_dataset(
                    "gem_group",
                    data=np.array([str((i % 4) + 1) for i in range(n_ntc + n_tgt)], dtype=h5py.string_dtype("utf-8")),
                )
                vg = f.create_group("var")
                vg.create_dataset(
                    "gene_name",
                    data=np.array([f"G{j:03d}" for j in range(n_g)], dtype=h5py.string_dtype("utf-8")),
                )

            self.assertEqual(_resolve_nadig_h5ad_path_r16("nadig_hepg2", raw_nadig_dir=ndir), hepg2_p)
            old_env = os.environ.get("NADIG_DATA_DIR")
            try:
                os.environ["NADIG_DATA_DIR"] = str(ndir)
                self.assertEqual(_resolve_nadig_h5ad_path_r16("nadig_hepg2"), hepg2_p)
            finally:
                if old_env is None:
                    os.environ.pop("NADIG_DATA_DIR", None)
                else:
                    os.environ["NADIG_DATA_DIR"] = old_env

            nb = extract_nadig_h5ad_ntc_bundle_r16(hepg2_p, seed=DEFAULT_R16_SEED, include_resample=True)
            self.assertEqual(nb["raw_counts"].shape, (512, n_g))
            self.assertEqual(nb["resample_raw_counts"].shape, (512, n_g))
            self.assertEqual(len(nb["barcodes"]), 512)
            self.assertEqual(len(nb["resample_barcodes"]), 512)
            self.assertEqual(len(set(nb["barcodes"]) & set(nb["resample_barcodes"])), 0)
            self.assertTrue(all(bc.startswith("NADIG_BC_") for bc in nb["barcodes"]))

            # 2. Synthetic A549 `.mtx.gz` read-only NTC + barcode extractor (ZERO files written to raw_dir)
            adir = tdir / "a549_raw" / "a549_ko"
            adir.mkdir(parents=True)
            n_a_ntc, n_a_tgt, n_a_genes = 280, 10, 520  # 520 >= 500 min_ngenes + 1 CAS9 row excluded
            with gzip.open(adir / "GSE345058_gene_annotations.csv.gz", "wt", encoding="utf-8") as fg:
                fg.write("gene_symbol,ensembl_id,feature_types,genome\n")
                for j in range(n_a_genes):
                    fg.write(f"AG{j:04d},ENSG{j:08d},Gene Expression,GRCh38\n")
                fg.write("Cas9,Cas9,Gene Expression,CAS9\n")
            with gzip.open(adir / "GSE345058_cell_annotations.csv.gz", "wt", encoding="utf-8") as fc:
                fc.write("cell_id,batch,Lane,sgRNA,Perturbation,Condition\n")
                for i in range(n_a_ntc):
                    b_num = (i % 4) + 1
                    fc.write(f"A549_CID_{i:05d},batch{b_num},Lane1-{b_num},NTC_sg1,Nontargeting,Nontargeting\n")
                for i in range(n_a_tgt):
                    fc.write(f"A549_TGT_{i:05d},batch1,Lane1-1,TP53_sg1,TP53,Targeting\n")
            n_a_cells = n_a_ntc + n_a_tgt
            nnz_a = n_a_cells * (n_a_genes + 1)
            with gzip.open(adir / "GSE345058_raw_counts.mtx.gz", "wt", encoding="utf-8") as fm:
                fm.write("%%MatrixMarket matrix coordinate integer general\n")
                fm.write(f"{n_a_cells} {n_a_genes + 1} {nnz_a}\n")
                for i in range(1, n_a_cells + 1):
                    for j in range(1, n_a_genes + 2):
                        fm.write(f"{i} {j} 2\n")

            files_before_a = sorted(p.name for p in adir.iterdir())
            ab = _extract_a549_ntc_bundle_r16(
                tdir / "a549_raw", seed=DEFAULT_R16_SEED, resample_seed=DEFAULT_RESAMPLE_SEED,
                target_n=128, include_resample=True,
            )
            self.assertEqual(sorted(p.name for p in adir.iterdir()), files_before_a)
            self.assertEqual(ab["raw_counts"].shape, (128, n_a_genes))  # CAS9 excluded!
            self.assertEqual(ab["resample_raw_counts"].shape, (128, n_a_genes))
            self.assertEqual(len(set(ab["barcodes"]) & set(ab["resample_barcodes"])), 0)
            self.assertTrue(all(bc.startswith("A549_CID_") for bc in ab["barcodes"]))

            # 3. Synthetic Feng iPSC read-only NTC + barcode extractor (`kolf_2` & `unassigned` excluded)
            fdir = tdir / "feng_raw" / "feng_ipsc"
            fdir.mkdir(parents=True)
            f_cids = [f"FENG_NTC_{i:04d}" for i in range(280)] + ["FENG_KOLF_1", "FENG_UNASS_1", "FENG_TGT_1"]
            with gzip.open(fdir / "TargetedScreen_Cell-Metadata.tsv.gz", "wt", encoding="utf-8") as fm:
                fm.write("Cell_ID\tCell_Line\tGuide_Call\tTarget_Gene\tInlet\n")
                for i in range(280):
                    ln_s = "line_A" if i % 2 == 0 else "line_B"
                    fm.write(f"FENG_NTC_{i:04d}\t{ln_s}\tNonTarget_01\t\t1\n")
                fm.write("FENG_KOLF_1\tkolf_2\tNonTarget_01\t\t1\n")
                fm.write("FENG_UNASS_1\tline_A\tunassigned\t\t1\n")
                fm.write("FENG_TGT_1\tline_A\tTP53_1\tTP53\t1\n")
            with gzip.open(fdir / "TargetedScreen_RNA-UMI-Counts.csv.gz", "wt", encoding="utf-8") as fc:
                fc.write("," + ",".join(f_cids) + "\n")
                for j in range(12):
                    fc.write(f"ENSG{j:08d}:FG{j:02d}:Gene-Expression," + ",".join(["3"] * len(f_cids)) + "\n")
                fc.write("GUIDE_1:NonTarget_01:CRISPR-Guide-Capture," + ",".join(["10"] * len(f_cids)) + "\n")

            files_before_f = sorted(p.name for p in fdir.iterdir())
            fb = _extract_feng_ipsc_ntc_bundle_r16(
                tdir / "feng_raw", seed=DEFAULT_R16_SEED, resample_seed=DEFAULT_RESAMPLE_SEED,
                target_n=128, include_resample=True,
            )
            self.assertEqual(sorted(p.name for p in fdir.iterdir()), files_before_f)
            self.assertEqual(fb["raw_counts"].shape, (128, 12))
            self.assertEqual(fb["resample_raw_counts"].shape, (128, 12))
            self.assertEqual(len(set(fb["barcodes"]) & set(fb["resample_barcodes"])), 0)
            self.assertTrue(all(bc.startswith("FENG_NTC_") for bc in fb["barcodes"]))

            # Also verify R-style write.csv header (`col_offset == 1`, no leading empty header field)
            fdir_r = tdir / "feng_raw_rstyle" / "feng_ipsc"
            fdir_r.mkdir(parents=True)
            (fdir_r / "TargetedScreen_Cell-Metadata.tsv.gz").write_bytes(
                (fdir / "TargetedScreen_Cell-Metadata.tsv.gz").read_bytes()
            )
            with gzip.open(fdir_r / "TargetedScreen_RNA-UMI-Counts.csv.gz", "wt", encoding="utf-8") as fc:
                fc.write(",".join(f_cids) + "\n")
                for j in range(12):
                    fc.write(f"ENSG{j:08d}:FG{j:02d}:Gene-Expression," + ",".join(["3"] * len(f_cids)) + "\n")
            fb_r = _extract_feng_ipsc_ntc_bundle_r16(
                tdir / "feng_raw_rstyle", seed=DEFAULT_R16_SEED, resample_seed=DEFAULT_RESAMPLE_SEED,
                target_n=128, include_resample=True,
            )
            self.assertEqual(fb_r["barcodes"], fb["barcodes"])
            self.assertEqual(fb_r["resample_barcodes"], fb["resample_barcodes"])
            self.assertTrue(np.array_equal(fb_r["raw_counts"], fb["raw_counts"]))

    def test_18_cross_process_pythonhashseed_determinism(self) -> None:
        """
        Pre-freeze reproducibility gate (`2026-11-26 01:01` & Addendum A Item 11 `01:31`):
        Verify across two fresh independent subprocesses with distinct `PYTHONHASHSEED`
        (`1` vs `99991`) strictly inside `tempfile.TemporaryDirectory()`:
          1. Python hash randomization genuinely differs (`hash('nadig_hepg2')` differs);
          2. Control sampling (`20261601` primary + `20261602` resample) and full synthetic
             `run_eval_r16` pipeline export bit-for-bit identical `.npz` arrays, `embeddings/x_state_*.npz`,
             barcodes, and metric/adjudication summaries (`max_abs == 0.0`);
          3. When `arc-state` + `SE-600M/model.safetensors` weights are present, official
             `encode_controls_se600m_r16` (`FilteredGenesCounts` + `VCIDatasetSentenceCollator`,
             `training=False`, `batch_size=16`, `shuffle=False`, `num_workers=0`) on synthetic
             integer UMI counts with zero-count genes (actively exercising `loader.py:529-533`
             `torch.randperm` tie-breaking) produces bit-for-bit identical (`max_abs == 0.0`)
             `2,058`-dim `X_state` embeddings across processes for BOTH `seed=20261601` and
             `seed=20261602` on CPU (and `max_abs <= 1e-6` on MPS when available), while
             `seed=20261601` and `seed=20261602` genuinely differ on zero-count tie-breaks;
          4. When `arc-state` + `ST-{HVG,SE}-Replogle/zeroshot/hepg2` weights are present, full
             `ArcStateAdapter('STATE_HVG')` and primary `ArcStateAdapter('STATE_SE')` inference
             on identical synthetic inputs is bit-for-bit identical on CPU (`0.0`) and `<= 1e-6`
             on MPS (when available).
        """
        with tempfile.TemporaryDirectory(prefix="r16_det_test_") as tmp:
            tdir = Path(tmp)
            ckpt_hvg = DEFAULT_MODELS_DIR / "ST-HVG-Replogle" / "zeroshot" / "hepg2" / "checkpoints" / "final.ckpt"
            ckpt_se = DEFAULT_MODELS_DIR / "ST-SE-Replogle" / "zeroshot" / "hepg2" / "checkpoints" / "final.ckpt"
            se_model_file = DEFAULT_MODELS_DIR / "SE-600M" / "model.safetensors"
            se_pe_file = DEFAULT_MODELS_DIR / "SE-600M" / "protein_embeddings.pt"

            counts_npz_arg = ""
            if _HAS_ARC_STATE and ckpt_hvg.is_file():
                import pickle

                with open(ckpt_hvg.parent.parent / "var_dims.pkl", "rb") as fh:
                    g_model = list(pickle.load(fh)["gene_names"])
                rng = np.random.default_rng(DEFAULT_R16_SEED)
                rates = rng.uniform(0.1, 3.0, size=(1, len(g_model)))
                raw_counts = rng.poisson(rates, size=(64, len(g_model))).astype(np.float32)
                # Synthetic 2,058-dim X_state embeddings ([L2-norm 2048 || 10]) for primary ST-SE check
                emb_2048 = rng.normal(0.0, 1.0, size=(64, 2048)).astype(np.float32)
                emb_2048 /= np.maximum(np.linalg.norm(emb_2048, axis=1, keepdims=True), 1e-12)
                ds_10 = rng.normal(0.0, 0.15, size=(64, 10)).astype(np.float32)
                x_state_syn = np.concatenate([emb_2048, ds_10], axis=1).astype(np.float32)
                counts_npz_path = tdir / "synth_counts.npz"
                np.savez_compressed(
                    counts_npz_path,
                    raw_counts=raw_counts,
                    x_state_syn=x_state_syn,
                    has_se_ckpt=np.asarray(bool(ckpt_se.is_file())),
                    g_model=np.asarray(g_model, dtype=str),
                )
                counts_npz_arg = str(counts_npz_path)

            se_npz_arg = ""
            if _HAS_ARC_STATE and se_model_file.is_file() and se_pe_file.is_file():
                import torch

                pe = torch.load(se_pe_file, map_location="cpu", weights_only=False)
                mapped_genes = list(pe.keys())[:2500]
                dup_gene = next(g for g in mapped_genes if f"{g}-1" not in pe)
                del pe
                # Include 1 duplicate mapped symbol (exercising `var_names_make_unique()`) + 10 unmapped genes
                se_gene_names = mapped_genes + [dup_gene] + [f"UNMAPPED_SYN_{k:02d}" for k in range(10)]
                rng_se = np.random.default_rng(DEFAULT_R16_SEED)
                # 4 synthetic cells (fast CPU forward pass, batch_size=16) with 800 nonzero integer UMI
                # genes and 1,700 zero-count mapped genes so `loader.py:529-533` `torch.randperm`
                # tie-breaking across zero-count genes is actively exercised!
                se_raw = np.zeros((4, len(se_gene_names)), dtype=np.float32)
                for ci in range(4):
                    nz = rng_se.choice(2500, size=800, replace=False)
                    se_raw[ci, nz] = rng_se.integers(1, 50, size=800).astype(np.float32)
                    se_raw[ci, nz[0]] = 50.0  # guarantees is_raw_integer_counts() == True (>35)
                    se_raw[ci, 2500:] = 3.0   # exercises duplicate var_names_make_unique + valid_gene_mask
                se_npz_path = tdir / "synth_se600m.npz"
                np.savez_compressed(
                    se_npz_path,
                    raw_counts=se_raw,
                    gene_names=np.asarray(se_gene_names, dtype=str),
                )
                se_npz_arg = str(se_npz_path)

            worker_py = tdir / "determinism_worker.py"
            worker_py.write_text(
                "\n".join(
                    [
                        "import json, os, random, sys",
                        "from pathlib import Path",
                        "import numpy as np",
                        "src_dir, out_dir_s, counts_npz_s, models_dir_s, se_npz_s = sys.argv[1:6]",
                        "sys.path.insert(0, src_dir)",
                        "from metrics_r16 import DEFAULT_R16_SEED, DEFAULT_RESAMPLE_SEED",
                        "from model_adapters_r16 import (",
                        "    ArcStateAdapter,",
                        "    encode_controls_se600m_r16,",
                        "    sample_disjoint_resample_control_sets_r16,",
                        "    sample_stratified_control_sets_r16,",
                        ")",
                        "from run_eval_r16 import run_synthetic_smoke_r16, save_control_embeddings_r16",
                        "out_dir = Path(out_dir_s)",
                        "out_dir.mkdir(parents=True, exist_ok=True)",
                        "batches = ['B0'] * 400 + ['B1'] * 400 + ['B2'] * 400",
                        "p_sets = sample_stratified_control_sets_r16(1200, batch_labels=batches, seed=DEFAULT_R16_SEED)",
                        "r_sets, _ = sample_disjoint_resample_control_sets_r16(1200, p_sets, batch_labels=batches, seed=DEFAULT_RESAMPLE_SEED)",
                        "np.savez_compressed(out_dir / 'sampling.npz', primary_sets=p_sets, resample_sets=r_sets)",
                        "run_synthetic_smoke_r16(out_dir / 'smoke', n_boot=120, seed=DEFAULT_R16_SEED)",
                        "ran_hvg_cpu = False",
                        "ran_hvg_mps = False",
                        "ran_se_cpu = False",
                        "ran_se_mps = False",
                        "ran_se600m_cpu = False",
                        "ran_se600m_mps = False",
                        "se600m_primary_audit = {}",
                        "se600m_resample_audit = {}",
                        "se600m_mps_audit = {}",
                        "if counts_npz_s:",
                        "    import torch",
                        "    zc = np.load(counts_npz_s, allow_pickle=False)",
                        "    raw_c = zc['raw_counts']",
                        "    x_st = zc['x_state_syn']",
                        "    has_se = bool(zc['has_se_ckpt'])",
                        "    g_mod = [str(x) for x in zc['g_model'].tolist()]",
                        "    perts = ['AAAS', 'AAMP', 'AAR2', 'AARS', 'AARS2']",
                        "    ev_g = g_mod[:20]",
                        "    ad_cpu = ArcStateAdapter(",
                        "        'STATE_HVG', checkpoint_id='zeroshot/hepg2',",
                        "        models_dir=Path(models_dir_s), device='cpu',",
                        "        batch_mode='fixed_index_zero', hvg_norm_mode='official_median_log1p',",
                        "    )",
                        "    d_cpu = ad_cpu.predict_deltas(",
                        "        target_line='nadig_hepg2', control_line_label='synth_det',",
                        "        ctrl_cells=raw_c, pert_names=perts, gene_names=ev_g,",
                        "        ctrl_raw_gene_names=g_mod, allow_synthetic_without_freeze=True,",
                        "    )",
                        "    np.savez_compressed(out_dir / 'hvg_delta_cpu.npz', delta_pred=d_cpu)",
                        "    ran_hvg_cpu = True",
                        "    if torch.backends.mps.is_available():",
                        "        ad_mps = ArcStateAdapter(",
                        "            'STATE_HVG', checkpoint_id='zeroshot/hepg2',",
                        "            models_dir=Path(models_dir_s), device='mps',",
                        "            batch_mode='fixed_index_zero', hvg_norm_mode='official_median_log1p',",
                        "        )",
                        "        d_mps = ad_mps.predict_deltas(",
                        "            target_line='nadig_hepg2', control_line_label='synth_det',",
                        "            ctrl_cells=raw_c, pert_names=perts, gene_names=ev_g,",
                        "            ctrl_raw_gene_names=g_mod, allow_synthetic_without_freeze=True,",
                        "        )",
                        "        np.savez_compressed(out_dir / 'hvg_delta_mps.npz', delta_pred=d_mps)",
                        "        ran_hvg_mps = True",
                        "    if has_se:",
                        "        se_cpu = ArcStateAdapter(",
                        "            'STATE_SE', checkpoint_id='zeroshot/hepg2',",
                        "            models_dir=Path(models_dir_s), device='cpu',",
                        "            batch_mode='fixed_index_zero',",
                        "        )",
                        "        ds_cpu = se_cpu.predict_deltas(",
                        "            target_line='nadig_hepg2', control_line_label='synth_det',",
                        "            ctrl_cells=x_st, pert_names=perts, gene_names=ev_g,",
                        "            ctrl_is_embedded_x_state=True, allow_synthetic_without_freeze=True,",
                        "        )",
                        "        np.savez_compressed(out_dir / 'se_delta_cpu.npz', delta_pred=ds_cpu)",
                        "        ran_se_cpu = True",
                        "        if torch.backends.mps.is_available():",
                        "            se_mps = ArcStateAdapter(",
                        "                'STATE_SE', checkpoint_id='zeroshot/hepg2',",
                        "                models_dir=Path(models_dir_s), device='mps',",
                        "                batch_mode='fixed_index_zero',",
                        "            )",
                        "            ds_mps = se_mps.predict_deltas(",
                        "                target_line='nadig_hepg2', control_line_label='synth_det',",
                        "                ctrl_cells=x_st, pert_names=perts, gene_names=ev_g,",
                        "                ctrl_is_embedded_x_state=True, allow_synthetic_without_freeze=True,",
                        "            )",
                        "            np.savez_compressed(out_dir / 'se_delta_mps.npz', delta_pred=ds_mps)",
                        "            ran_se_mps = True",
                        "if se_npz_s:",
                        "    import torch",
                        "    zs = np.load(se_npz_s, allow_pickle=False)",
                        "    se_raw_c = zs['raw_counts']",
                        "    se_gnames = [str(x) for x in zs['gene_names'].tolist()]",
                        "    se_dir = Path(models_dir_s) / 'SE-600M'",
                        "    h_int = int(os.environ.get('PYTHONHASHSEED', '1'))",
                        "    def _scramble_rng(tag_n):",
                        "        s_val = (h_int * 31 + tag_n * 17) & 0x7FFFFFFF",
                        "        random.seed(s_val)",
                        "        np.random.seed(s_val)",
                        "        torch.manual_seed(s_val)",
                        "        _ = torch.rand(17 if h_int == 1 else 251)",
                        "    _scramble_rng(1)",
                        "    x_se_p, se600m_primary_audit = encode_controls_se600m_r16(",
                        "        raw_counts=se_raw_c, raw_gene_names=se_gnames, se_dir=se_dir,",
                        "        device='cpu', batch_size=16, dataset_name='synth_se_primary',",
                        "        seed=DEFAULT_R16_SEED,",
                        "    )",
                        "    _scramble_rng(2)",
                        "    x_se_r, se600m_resample_audit = encode_controls_se600m_r16(",
                        "        raw_counts=se_raw_c, raw_gene_names=se_gnames, se_dir=se_dir,",
                        "        device='cpu', batch_size=16, dataset_name='synth_se_resample',",
                        "        seed=DEFAULT_RESAMPLE_SEED,",
                        "    )",
                        "    se_bcs = [f'SYN_SE_{i:02d}' for i in range(se_raw_c.shape[0])]",
                        "    save_control_embeddings_r16(",
                        "        out_dir / 'se600m_emb',",
                        "        {",
                        "            'synth_se_primary': {",
                        "                'X_state': x_se_p, 'control_barcodes': se_bcs,",
                        "                'seed': DEFAULT_R16_SEED, 'dataset_name': 'synth_se_primary',",
                        "                'se_audit': se600m_primary_audit,",
                        "            },",
                        "            'synth_se_resample': {",
                        "                'X_state': x_se_r, 'control_barcodes': se_bcs,",
                        "                'seed': DEFAULT_RESAMPLE_SEED, 'dataset_name': 'synth_se_resample',",
                        "                'se_audit': se600m_resample_audit,",
                        "            },",
                        "        },",
                        "    )",
                        "    ran_se600m_cpu = True",
                        "    if torch.backends.mps.is_available():",
                        "        _scramble_rng(3)",
                        "        x_se_m, se600m_mps_audit = encode_controls_se600m_r16(",
                        "            raw_counts=se_raw_c, raw_gene_names=se_gnames, se_dir=se_dir,",
                        "            device='mps', batch_size=16, dataset_name='synth_se_mps',",
                        "            seed=DEFAULT_R16_SEED,",
                        "        )",
                        "        np.savez_compressed(out_dir / 'se600m_mps_primary.npz', X_state=x_se_m)",
                        "        ran_se600m_mps = True",
                        "meta = {",
                        "    'pythonhashseed_env': os.environ.get('PYTHONHASHSEED'),",
                        "    'probe_hash': hash('nadig_hepg2'),",
                        "    'ran_hvg_cpu': ran_hvg_cpu,",
                        "    'ran_hvg_mps': ran_hvg_mps,",
                        "    'ran_se_cpu': ran_se_cpu,",
                        "    'ran_se_mps': ran_se_mps,",
                        "    'ran_se600m_cpu': ran_se600m_cpu,",
                        "    'ran_se600m_mps': ran_se600m_mps,",
                        "    'se600m_primary_audit': se600m_primary_audit,",
                        "    'se600m_resample_audit': se600m_resample_audit,",
                        "    'se600m_mps_audit': se600m_mps_audit,",
                        "}",
                        "(out_dir / 'meta.json').write_text(json.dumps(meta, indent=2), encoding='utf-8')",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            run_dirs: dict[str, Path] = {}
            metas: dict[str, dict] = {}
            for hseed in ("1", "99991"):
                rdir = tdir / f"run_{hseed}"
                env = os.environ.copy()
                env.pop("PYTHONPATH", None)
                env["PYTHONHASHSEED"] = hseed
                env["WANDB_MODE"] = "disabled"
                env["HF_HUB_OFFLINE"] = "1"
                env["TRANSFORMERS_OFFLINE"] = "1"
                proc = subprocess.run(
                    [
                        sys.executable,
                        "-s",
                        str(worker_py),
                        str(SRC_DIR),
                        str(rdir),
                        counts_npz_arg,
                        str(DEFAULT_MODELS_DIR),
                        se_npz_arg,
                    ],
                    env=env,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertEqual(
                    proc.returncode,
                    0,
                    f"Determinism worker failed for PYTHONHASHSEED={hseed}:\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}",
                )
                run_dirs[hseed] = rdir
                metas[hseed] = json.loads((rdir / "meta.json").read_text(encoding="utf-8"))

            self.assertEqual(metas["1"]["pythonhashseed_env"], "1")
            self.assertEqual(metas["99991"]["pythonhashseed_env"], "99991")
            self.assertNotEqual(metas["1"]["probe_hash"], metas["99991"]["probe_hash"])

            # 1. Exact sampling parity across PYTHONHASHSEED=1 vs 99991
            s1 = np.load(run_dirs["1"] / "sampling.npz", allow_pickle=False)
            s2 = np.load(run_dirs["99991"] / "sampling.npz", allow_pickle=False)
            self.assertTrue(np.array_equal(s1["primary_sets"], s2["primary_sets"]))
            self.assertTrue(np.array_equal(s1["resample_sets"], s2["resample_sets"]))

            # 2. Exact full synthetic pipeline parity (barcodes, summaries, embeddings, and every .npz array bit-for-bit)
            bc1 = (run_dirs["1"] / "smoke" / "control_barcodes_synthetic_smoke.json").read_text(encoding="utf-8")
            bc2 = (run_dirs["99991"] / "smoke" / "control_barcodes_synthetic_smoke.json").read_text(encoding="utf-8")
            self.assertEqual(bc1, bc2)

            sum1 = json.loads((run_dirs["1"] / "smoke" / "summary_synthetic_smoke.json").read_text(encoding="utf-8"))
            sum2 = json.loads((run_dirs["99991"] / "smoke" / "summary_synthetic_smoke.json").read_text(encoding="utf-8"))
            # Strip ZIP-mtime-dependent .npz file SHA-256 digests after verifying identical exported file sets
            m1_keys = sorted(sum1["lines"]["synthetic_a549_ko"].pop("array_manifest_sha256").keys())
            m2_keys = sorted(sum2["lines"]["synthetic_a549_ko"].pop("array_manifest_sha256").keys())
            self.assertEqual(m1_keys, m2_keys)
            e1_keys = sorted(sum1.pop("embeddings_manifest_sha256").keys())
            e2_keys = sorted(sum2.pop("embeddings_manifest_sha256").keys())
            self.assertEqual(e1_keys, e2_keys)
            self.assertGreater(len(e1_keys), 0)
            self.assertEqual(json.dumps(sum1, sort_keys=True), json.dumps(sum2, sort_keys=True))

            for subdir_rel in (Path("arrays") / "synthetic_a549_ko", Path("embeddings")):
                arr_dir_1 = run_dirs["1"] / "smoke" / subdir_rel
                arr_dir_2 = run_dirs["99991"] / "smoke" / subdir_rel
                npz_files_1 = sorted(p.name for p in arr_dir_1.glob("*.npz"))
                npz_files_2 = sorted(p.name for p in arr_dir_2.glob("*.npz"))
                self.assertEqual(npz_files_1, npz_files_2)
                self.assertGreater(len(npz_files_1), 0)
                for fname in npz_files_1:
                    z1 = np.load(arr_dir_1 / fname, allow_pickle=False)
                    z2 = np.load(arr_dir_2 / fname, allow_pickle=False)
                    self.assertEqual(sorted(z1.files), sorted(z2.files))
                    for k in z1.files:
                        a1, a2 = z1[k], z2[k]
                        if np.issubdtype(a1.dtype, np.floating):
                            self.assertTrue(
                                np.array_equal(a1, a2, equal_nan=True), f"Float mismatch in {subdir_rel}/{fname}:{k}"
                            )
                        else:
                            self.assertTrue(
                                np.array_equal(a1, a2), f"Array mismatch in {subdir_rel}/{fname}:{k}"
                            )

            # 3. Full real-weight ArcStateAdapter (`STATE_HVG` & primary `STATE_SE`) determinism (CPU bit-exact; MPS <= 1e-6)
            if counts_npz_arg:
                self.assertTrue(metas["1"]["ran_hvg_cpu"] and metas["99991"]["ran_hvg_cpu"])
                dc1 = np.load(run_dirs["1"] / "hvg_delta_cpu.npz", allow_pickle=False)["delta_pred"]
                dc2 = np.load(run_dirs["99991"] / "hvg_delta_cpu.npz", allow_pickle=False)["delta_pred"]
                self.assertEqual(dc1.shape, (5, 20))
                self.assertTrue(np.all(np.isfinite(dc1)))
                self.assertTrue(np.array_equal(dc1, dc2), f"HVG CPU max_abs={float(np.max(np.abs(dc1 - dc2)))}")
                if metas["1"]["ran_hvg_mps"] and metas["99991"]["ran_hvg_mps"]:
                    dm1 = np.load(run_dirs["1"] / "hvg_delta_mps.npz", allow_pickle=False)["delta_pred"]
                    dm2 = np.load(run_dirs["99991"] / "hvg_delta_mps.npz", allow_pickle=False)["delta_pred"]
                    self.assertTrue(np.all(np.isfinite(dm1)) and np.all(np.isfinite(dm2)))
                    self.assertLessEqual(float(np.max(np.abs(dm1 - dm2))), 1e-6)
                if metas["1"]["ran_se_cpu"] and metas["99991"]["ran_se_cpu"]:
                    sc1 = np.load(run_dirs["1"] / "se_delta_cpu.npz", allow_pickle=False)["delta_pred"]
                    sc2 = np.load(run_dirs["99991"] / "se_delta_cpu.npz", allow_pickle=False)["delta_pred"]
                    self.assertEqual(sc1.shape, (5, 20))
                    self.assertTrue(np.all(np.isfinite(sc1)))
                    self.assertTrue(np.array_equal(sc1, sc2), f"ST-SE CPU max_abs={float(np.max(np.abs(sc1 - sc2)))}")
                if metas["1"]["ran_se_mps"] and metas["99991"]["ran_se_mps"]:
                    sm1 = np.load(run_dirs["1"] / "se_delta_mps.npz", allow_pickle=False)["delta_pred"]
                    sm2 = np.load(run_dirs["99991"] / "se_delta_mps.npz", allow_pickle=False)["delta_pred"]
                    self.assertTrue(np.all(np.isfinite(sm1)) and np.all(np.isfinite(sm2)))
                    self.assertLessEqual(float(np.max(np.abs(sm1 - sm2))), 1e-6)

            # 4. Addendum A Item 11: Real-weight `SE-600M` official-collator cross-process determinism
            if se_npz_arg:
                self.assertTrue(metas["1"]["ran_se600m_cpu"] and metas["99991"]["ran_se600m_cpu"])
                for hseed in ("1", "99991"):
                    pa = metas[hseed]["se600m_primary_audit"]
                    ra = metas[hseed]["se600m_resample_audit"]
                    self.assertEqual(pa["seed"], DEFAULT_R16_SEED)
                    self.assertEqual(ra["seed"], DEFAULT_RESAMPLE_SEED)
                    self.assertEqual(pa["batch_size"], 16)
                    self.assertFalse(pa["shuffle"])
                    self.assertEqual(pa["num_workers"], 0)
                    self.assertFalse(pa["collator_training"])
                    self.assertEqual(pa["n_mapped_to_pe19790"], 2500)
                    self.assertEqual(pa["n_mapped_to_pe19790_raw_rows"], 2501)
                    self.assertEqual(ra["n_mapped_to_pe19790"], 2500)
                    self.assertEqual(ra["n_mapped_to_pe19790_raw_rows"], 2501)

                    emb_man_p = run_dirs[hseed] / "se600m_emb" / "embeddings_manifest_sha256.json"
                    self.assertTrue(emb_man_p.is_file())
                    emb_man = json.loads(emb_man_p.read_text(encoding="utf-8"))
                    self.assertEqual(
                        set(emb_man.keys()),
                        {"x_state_synth_se_primary.npz", "x_state_synth_se_resample.npz"},
                    )
                    for efname, esha in emb_man.items():
                        self.assertEqual(
                            hashlib.sha256((run_dirs[hseed] / "se600m_emb" / efname).read_bytes()).hexdigest(),
                            esha,
                        )

                ep1 = np.load(run_dirs["1"] / "se600m_emb" / "x_state_synth_se_primary.npz", allow_pickle=False)["X_state"]
                ep2 = np.load(run_dirs["99991"] / "se600m_emb" / "x_state_synth_se_primary.npz", allow_pickle=False)["X_state"]
                er1 = np.load(run_dirs["1"] / "se600m_emb" / "x_state_synth_se_resample.npz", allow_pickle=False)["X_state"]
                er2 = np.load(run_dirs["99991"] / "se600m_emb" / "x_state_synth_se_resample.npz", allow_pickle=False)["X_state"]
                self.assertEqual(ep1.shape, (4, 2058))
                self.assertEqual(er1.shape, (4, 2058))
                self.assertEqual(ep1.dtype, np.float32)
                self.assertEqual(er1.dtype, np.float32)
                self.assertTrue(np.all(np.isfinite(ep1)) and np.all(np.isfinite(er1)))
                self.assertTrue(
                    np.array_equal(ep1, ep2),
                    f"SE-600M CPU primary (20261601) max_abs={float(np.max(np.abs(ep1 - ep2)))}",
                )
                self.assertTrue(
                    np.array_equal(er1, er2),
                    f"SE-600M CPU resample (20261602) max_abs={float(np.max(np.abs(er1 - er2)))}",
                )
                self.assertEqual(
                    metas["1"]["se600m_primary_audit"]["x_state_sha256"],
                    metas["99991"]["se600m_primary_audit"]["x_state_sha256"],
                )
                self.assertEqual(
                    metas["1"]["se600m_primary_audit"]["x_state_sha256"],
                    hashlib.sha256(np.ascontiguousarray(ep1, dtype=np.float32).tobytes()).hexdigest(),
                )
                self.assertEqual(
                    metas["1"]["se600m_resample_audit"]["x_state_sha256"],
                    metas["99991"]["se600m_resample_audit"]["x_state_sha256"],
                )
                self.assertEqual(
                    metas["1"]["se600m_resample_audit"]["x_state_sha256"],
                    hashlib.sha256(np.ascontiguousarray(er1, dtype=np.float32).tobytes()).hexdigest(),
                )
                # Confirm zero-count `torch.randperm` tie-breaking is genuinely exercised (different seeds differ)
                self.assertFalse(
                    np.array_equal(ep1, er1),
                    "Expected seed=20261601 vs 20261602 to differ under zero-count torch.randperm tie-breaking.",
                )
                self.assertNotEqual(
                    metas["1"]["se600m_primary_audit"]["x_state_sha256"],
                    metas["1"]["se600m_resample_audit"]["x_state_sha256"],
                )
                if metas["1"]["ran_se600m_mps"] and metas["99991"]["ran_se600m_mps"]:
                    for hseed in ("1", "99991"):
                        ma = metas[hseed]["se600m_mps_audit"]
                        self.assertEqual(ma["device"], "mps")
                        self.assertEqual(ma["seed"], DEFAULT_R16_SEED)
                        self.assertEqual(ma["batch_size"], 16)
                        self.assertFalse(ma["shuffle"])
                        self.assertEqual(ma["num_workers"], 0)
                        self.assertFalse(ma["collator_training"])
                        self.assertEqual(ma["n_mapped_to_pe19790"], 2500)
                        self.assertEqual(ma["n_mapped_to_pe19790_raw_rows"], 2501)
                    em1 = np.load(run_dirs["1"] / "se600m_mps_primary.npz", allow_pickle=False)["X_state"]
                    em2 = np.load(run_dirs["99991"] / "se600m_mps_primary.npz", allow_pickle=False)["X_state"]
                    self.assertEqual(em1.shape, (4, 2058))
                    self.assertTrue(np.all(np.isfinite(em1)) and np.all(np.isfinite(em2)))
                    mps_max_abs = float(np.max(np.abs(em1 - em2)))
                    self.assertLessEqual(
                        mps_max_abs,
                        1e-6,
                        f"SE-600M MPS primary (20261601) max_abs={mps_max_abs}",
                    )


if __name__ == "__main__":
    unittest.main()
