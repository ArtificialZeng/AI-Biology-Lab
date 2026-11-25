#!/usr/bin/env python3
"""
Clean-Room Unit, Parity, & Real-Weight Synthetic Smoke Suite (`round16/audit/`, `A-T3`,
strictly aligned to `PREREG_I.md v5-FINAL-PRE-LOCK` & `Addendum A`).

Strictly NEVER reads, imports, or references `round16/src/` or any real biological cell data.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import h5py
import numpy as np
import scipy.sparse as sp

AUDIT_DIR = Path(__file__).resolve().parent
REPO_ROOT = AUDIT_DIR.parents[1]
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
R10_SCRIPTS_DIR = REPO_ROOT / "round10" / "shared" / "scripts"
for _p in (str(AUDIT_DIR), str(R10_SRC_DIR), str(R10_SCRIPTS_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from audit_r16_adjudicate import (  # noqa: E402
    DEFAULT_TOST_DELTA,
    adjudicate_claim_i1_line,
    adjudicate_claim_i2_line,
    adjudicate_group_r16,
    classify_pair_tost_r16,
)
from audit_r16_metrics import (  # noqa: E402
    DEFAULT_R16_SEED,
    compute_gene_label_perm_null_r16,
    compute_orth_suite_r16,
    compute_paired_contrast_r16,
    compute_paired_diff_vs_resample_r16,
    compute_practical_metrics_r16,
    compute_winner_curse_weight_shift_r16,
    summarize_perm100_r16,
)
from audit_r16_raw_controls import (  # noqa: E402
    GROUP3_PRIMARY_FCT_SOURCES,
    build_fct_sensitivity_on_locked_targets_r16,
    build_gse295214_parts,
    evaluate_gse295214_counts_gate_r16,
    extract_gse295214_ntc_raw_counts,
)
from audit_r16_sampling import (  # noqa: E402
    RESAMPLE_R16_SEED,
    build_mismatched_strata_r16,
    libnorm_log1p,
    sample_nonoverlapping_resample_sets_r16,
    sample_stratified_control_sets_r16,
)
from audit_r16_state_adapter import (  # noqa: E402
    BlindingGateError,
    CleanRoomArcStateAdapter,
    CleanRoomSE600MEncoder,
    TrainingLeakageError,
    align_hvg_expression_to_model,
    assert_freeze_gates_r16,
    assert_zero_training_leakage,
    encode_with_loaded_se_model,
    get_hgnc_bridge,
    map_batch_labels_to_indices,
    map_genes_to_se600m_pe,
    preprocess_hvg_expression,
)
from run_audit_r16 import (  # noqa: E402
    PREREG_I_V5_CONFIG,
    evaluate_and_adjudicate_line,
    export_line_arm_npz,
    main as runner_main,
    run_g4_factorial_for_line,
    run_synthetic_runner_smoke,
    verify_cross_party_arrays_and_barcodes,
    verify_nadig_hepg2_adata_pred_parity,
)
from metrics_r10 import (  # noqa: E402
    compute_joint_paired_r_orth_diff_r10,
    compute_residual_suite_r10,
)
from prepare_r10_lines import (  # noqa: E402
    build_hgnc_lookup,
    process_r10_line,
)
from synthetic_fixtures_r16 import (  # noqa: E402
    evaluate_synthetic_fixture_r16,
    export_canonical_parity_pack,
    make_batch_imbalanced_control_fixture,
    make_synthetic_pert_cohort,
)


class TestRound16CleanRoomAudit(unittest.TestCase):
    def test_1_exact_parity_vs_frozen_round10_r_orth_cal(self) -> None:
        """Primary `r_orth_cal` and paired contrast must match frozen `round10/src/metrics_r10.py` <= 1e-12."""
        fx = make_synthetic_pert_cohort(n_perts=75, n_genes=160, regime="matched_better_strong", seed=9101)
        mis_first = fx["mis_preds"][fx["seen_sources"][0]]

        r10_suite = compute_residual_suite_r10(
            delta_obs=fx["delta_obs"],
            delta_pred=fx["pred_matched"],
            delta_transfer=fx["t_mean"],
            target_gene_idx=fx["target_cols"],
            delta_half_a=fx["half_a_obs"],
            delta_half_b=fx["half_b_obs"],
            half_valid=fx["half_valid"],
            n_boot=250,
            seed=DEFAULT_R16_SEED,
        )["all_nontarget"]

        r16_suite = compute_orth_suite_r16(
            delta_obs=fx["delta_obs"],
            delta_pred=fx["pred_matched"],
            delta_transfer=fx["t_mean"],
            target_gene_idx=fx["target_cols"],
            source_grand_mean=None,
            delta_half_a=fx["half_a_obs"],
            delta_half_b=fx["half_b_obs"],
            half_valid=fx["half_valid"],
            mode="all_nontarget",
            n_boot=250,
            seed=DEFAULT_R16_SEED,
        )
        self.assertLessEqual(abs(r10_suite["r_orth_cal"]["estimate"] - r16_suite["estimate"]), 1e-12)
        self.assertLessEqual(abs(r10_suite["r_orth_cal"]["ci975_lo"] - r16_suite["ci975_lo"]), 1e-12)
        self.assertLessEqual(abs(r10_suite["r_orth_cal"]["ci975_hi"] - r16_suite["ci975_hi"]), 1e-12)
        self.assertLessEqual(abs(r10_suite["r_orth_cal"]["ci95_lo"] - r16_suite["ci95_lo"]), 1e-12)
        self.assertLessEqual(abs(r10_suite["r_orth_cal"]["ci95_hi"] - r16_suite["ci95_hi"]), 1e-12)
        self.assertLessEqual(abs(r10_suite["rho_half_ab_orth"] - r16_suite["rho_ab_orth"]), 1e-12)

        r10_pair = compute_joint_paired_r_orth_diff_r10(
            delta_obs=fx["delta_obs"],
            delta_pred_a=fx["pred_matched"],
            delta_pred_b=mis_first,
            delta_transfer=fx["t_mean"],
            target_gene_idx=fx["target_cols"],
            delta_half_a=fx["half_a_obs"],
            delta_half_b=fx["half_b_obs"],
            half_valid=fx["half_valid"],
            mode="all_nontarget",
            n_boot=250,
            seed=DEFAULT_R16_SEED,
        )
        r16_pair = compute_paired_contrast_r16(
            delta_obs=fx["delta_obs"],
            pred_a=fx["pred_matched"],
            pred_b=mis_first,
            delta_transfer=fx["t_mean"],
            target_gene_idx=fx["target_cols"],
            source_grand_mean=None,
            delta_half_a=fx["half_a_obs"],
            delta_half_b=fx["half_b_obs"],
            half_valid=fx["half_valid"],
            mode="all_nontarget",
            n_boot=250,
            seed=DEFAULT_R16_SEED,
        )
        self.assertLessEqual(abs(r10_pair["estimate"] - r16_pair["estimate"]), 1e-12)
        self.assertLessEqual(abs(r10_pair["ci975_lo"] - r16_pair["ci975_lo"]), 1e-12)
        self.assertLessEqual(abs(r10_pair["ci95_hi"] - r16_pair["ci95_hi"]), 1e-12)

    def test_2_property_p1_r_dorth_null_invariance(self) -> None:
        """Secondary `r_dorth_cal` annihilates all affine combinations in `span{1, t_mean, b_S}` <= 1e-14."""
        fx = make_synthetic_pert_cohort(n_perts=60, n_genes=140, regime="matched_better_strong", seed=9102)
        P, G = fx["n_perts"], fx["n_genes"]
        rng = np.random.default_rng(4242)

        sgm_only = np.broadcast_to(fx["grand_mean"][None, :], (P, G)).copy()
        tm_global_scale = 0.65 * fx["t_mean"]
        # Random per-perturbation intercept + slope on t_mean + slope on b_S
        a0 = rng.normal(0.0, 0.5, size=(P, 1))
        a1 = rng.normal(1.0, 0.4, size=(P, 1))
        a2 = rng.normal(0.5, 0.3, size=(P, 1))
        rand_affine = a0 + a1 * fx["t_mean"] + a2 * sgm_only

        # First confirm Source_Grand_Mean genuinely has non-zero r_orth_cal on this fixture
        orth_sgm = compute_orth_suite_r16(
            fx["delta_obs"], sgm_only, fx["t_mean"], fx["target_cols"],
            source_grand_mean=None, delta_half_a=fx["half_a_obs"], delta_half_b=fx["half_b_obs"],
            half_valid=fx["half_valid"], n_boot=120, seed=DEFAULT_R16_SEED,
        )
        self.assertGreater(abs(orth_sgm["estimate"]), 0.15)

        for name, pred in (
            ("Transfer_Mean", fx["t_mean"]),
            ("Transfer_Mean_GlobalScale", tm_global_scale),
            ("Source_Grand_Mean", sgm_only),
            ("Random_Affine", rand_affine),
        ):
            res = compute_orth_suite_r16(
                fx["delta_obs"], pred, fx["t_mean"], fx["target_cols"],
                source_grand_mean=fx["grand_mean"], delta_half_a=fx["half_a_obs"],
                delta_half_b=fx["half_b_obs"], half_valid=fx["half_valid"],
                n_boot=120, seed=DEFAULT_R16_SEED,
            )
            self.assertLessEqual(abs(res["estimate"]), 1e-14, msg=f"{name} leaked under r_dorth_cal: {res}")

        # Meanwhile real matched signal survives r_dorth_cal cleanly
        dorth_good = compute_orth_suite_r16(
            fx["delta_obs"], fx["pred_matched"], fx["t_mean"], fx["target_cols"],
            source_grand_mean=fx["grand_mean"], delta_half_a=fx["half_a_obs"],
            delta_half_b=fx["half_b_obs"], half_valid=fx["half_valid"],
            n_boot=120, seed=DEFAULT_R16_SEED,
        )
        self.assertGreater(dorth_good["ci975_lo"], 0.4)

    def test_3_batch_imbalanced_control_sampling_and_strata(self) -> None:
        """Verify batch-stratified 16x32 sampling, tie-breaking, truncation, exclusion floor, and strata."""
        bf = make_batch_imbalanced_control_fixture(seed=DEFAULT_R16_SEED)
        s1 = sample_stratified_control_sets_r16(
            qc_mask=bf["large_imbalanced"]["qc_mask"],
            batch_labels=bf["large_imbalanced"]["batch_labels"],
            seed=DEFAULT_R16_SEED,
        )
        s1_repeat = sample_stratified_control_sets_r16(
            qc_mask=bf["large_imbalanced"]["qc_mask"],
            batch_labels=bf["large_imbalanced"]["batch_labels"],
            seed=DEFAULT_R16_SEED,
        )
        self.assertFalse(s1["excluded"])
        self.assertEqual(s1["sampling_mode"], "stratified_by_batch")
        self.assertEqual(s1["n_sampled"], 512)
        self.assertEqual(s1["n_sets"], 16)
        self.assertEqual(s1["sets"].shape, (16, 32))
        self.assertEqual(len(np.unique(s1["selected_indices"])), 512)
        self.assertTrue(np.array_equal(s1["sets"], s1_repeat["sets"]))
        self.assertTrue(np.all(bf["large_imbalanced"]["qc_mask"][s1["selected_indices"]]))
        # Check exact largest-remainder allocation: totals [701, 317, 103, 61, 18] out of 1200
        self.assertEqual(sum(s1["batch_counts_sampled"].values()), 512)
        for b, nb in bf["large_imbalanced"]["expected_valid_counts"].items():
            exact = 512.0 * nb / 1200.0
            self.assertLessEqual(abs(s1["batch_counts_sampled"][b] - exact), 1.0)

        # Lexicographic remainder tie-breaker (`B_1` beats `B_2` on exact tie `170.496`)
        s_tie = sample_stratified_control_sets_r16(
            n_cells=len(bf["tie_remainder"]["batch_labels"]),
            batch_labels=bf["tie_remainder"]["batch_labels"],
            seed=DEFAULT_R16_SEED,
        )
        self.assertEqual(s_tie["batch_counts_sampled"], bf["tie_remainder"]["expected_alloc"])

        # Small NTC truncation (200 -> 192 = 6 x 32)
        s_small = sample_stratified_control_sets_r16(
            n_cells=bf["small_truncated"]["n_cells"],
            batch_labels=bf["small_truncated"]["batch_labels"],
            seed=DEFAULT_R16_SEED,
        )
        self.assertEqual(s_small["sampling_mode"], "small_ntc_truncated")
        self.assertEqual(s_small["n_sampled"], 192)
        self.assertEqual(s_small["sets"].shape, (6, 32))

        # Below exclusion floor (< 128)
        s_floor = sample_stratified_control_sets_r16(n_cells=bf["below_floor"]["n_cells"], seed=DEFAULT_R16_SEED)
        self.assertTrue(s_floor["excluded"])
        self.assertEqual(s_floor["n_sets"], 0)

        # Two-stratum local builder (`I14`, `Q4`: feng_ipsc retained in unseen!)
        st_a549 = build_mismatched_strata_r16("a549_ko", "zeroshot/ensemble4")
        self.assertEqual(st_a549["seen_by_checkpoint"], ("nadig_hepg2", "nadig_jurkat"))
        self.assertEqual(st_a549["unseen_by_checkpoint"], ("cd4t_rest", "dld1", "feng_ipsc"))
        st_hepg2 = build_mismatched_strata_r16("nadig_hepg2", "zeroshot/hepg2")
        self.assertEqual(st_hepg2["seen_by_checkpoint"], ("nadig_jurkat",))
        self.assertEqual(st_hepg2["unseen_by_checkpoint"], ("a549_ko", "cd4t_rest", "dld1", "feng_ipsc"))
        with self.assertRaises(ValueError):
            build_mismatched_strata_r16("a549_ko", "zeroshot/ensemble4", group2_qualifying_lines=("nourreddine_kolf21j",))

    def test_4_adjudicator_edge_cases_and_canonical_fixtures(self) -> None:
        """Verify all M1-M4 edge cases and canonical multi-regime synthetic fixtures."""
        # Low valid_frac forces within_delta=False & INCONCLUSIVE (M2)
        low_vf = {"estimate": 0.01, "ci95_lo": -0.01, "ci95_hi": 0.03, "ci975_lo": -0.015, "ci975_hi": 0.035, "valid_frac": 0.80}
        res_vf = classify_pair_tost_r16(low_vf, DEFAULT_TOST_DELTA)
        self.assertEqual(res_vf["category"], "INCONCLUSIVE")
        self.assertFalse(res_vf["within_delta"])

        # Open boundary equality (== +/- delta) forces INCONCLUSIVE
        at_bound = {"estimate": 0.02, "ci95_lo": -0.01, "ci95_hi": 0.04, "ci975_lo": -0.02, "ci975_hi": 0.05, "valid_frac": 1.0}
        self.assertEqual(classify_pair_tost_r16(at_bound, DEFAULT_TOST_DELTA)["category"], "INCONCLUSIVE")
        self.assertFalse(classify_pair_tost_r16(at_bound, DEFAULT_TOST_DELTA)["within_delta"])

        # Conservative tie rule (PREREG_I §3.1): ci975_lo == 0.0 or ci975_hi == 0.0 MUST return INCONCLUSIVE
        # BEFORE EQUIVALENT, even when the 95% CI falls strictly inside (-delta, +delta).
        tie_lo_in_delta = {
            "estimate": 0.014, "ci975_lo": 0.0, "ci975_hi": 0.030,
            "ci95_lo": 0.002, "ci95_hi": 0.025, "valid_frac": 1.0,
        }
        tie_hi_in_delta = {
            "estimate": -0.014, "ci975_lo": -0.030, "ci975_hi": 0.0,
            "ci95_lo": -0.025, "ci95_hi": -0.002, "valid_frac": 1.0,
        }
        tie_neg_zero = {
            "estimate": 0.010, "ci975_lo": -0.0, "ci975_hi": 0.025,
            "ci95_lo": 0.001, "ci95_hi": 0.020, "valid_frac": 1.0,
        }
        for tie_case in (tie_lo_in_delta, tie_hi_in_delta, tie_neg_zero):
            res_tie = classify_pair_tost_r16(tie_case, DEFAULT_TOST_DELTA)
            self.assertEqual(res_tie["category"], "INCONCLUSIVE")
            self.assertTrue(res_tie["within_delta"])

        # At line level, one 97.5% CI 0.0 tie among EQUIVALENT pairs prevents LINE_EQUIVALENT_* -> MIXED_OR_INCONCLUSIVE
        eq_clean = {"estimate": 0.005, "ci975_lo": -0.015, "ci975_hi": 0.025, "ci95_lo": -0.010, "ci95_hi": 0.020, "valid_frac": 1.0}
        line_with_tie = adjudicate_claim_i1_line(
            {"ci975_lo": 0.30, "ci975_hi": 0.45, "valid_frac": 1.0, "rho_ab_orth": 0.7},
            {"per_source": {"nadig_hepg2": eq_clean, "cd4t_rest": tie_lo_in_delta}},
            seen_sources=("nadig_hepg2",), unseen_sources=("cd4t_rest",),
        )
        self.assertEqual(line_with_tie["line_summary"], "MIXED_OR_INCONCLUSIVE")
        self.assertEqual(line_with_tie["counts"]["INCONCLUSIVE"], 1)

        # Symmetric check on H-I2 boundary tie
        self.assertEqual(
            adjudicate_claim_i2_line({"estimate": 0.01, "ci975_lo": 0.0, "ci975_hi": 0.02, "valid_frac": 1.0})["classification"],
            "INDISTINGUISHABLE",
        )
        self.assertEqual(
            adjudicate_claim_i2_line({"estimate": -0.01, "ci975_lo": -0.02, "ci975_hi": 0.0, "valid_frac": 1.0})["classification"],
            "INDISTINGUISHABLE",
        )

        # Key contract (M1): passing legacy `rho_ab` WITHOUT `rho_ab_orth` fails matched_arm_valid -> INCONCLUSIVE
        mb = {"estimate": 0.09, "ci975_lo": 0.06, "ci975_hi": 0.12, "ci95_lo": 0.065, "ci95_hi": 0.115, "valid_frac": 1.0}
        missing_orth_key = adjudicate_claim_i1_line(
            {"ci975_lo": 0.30, "ci975_hi": 0.45, "valid_frac": 1.0, "rho_ab": 0.7},
            {"per_source": {"nadig_hepg2": mb, "cd4t_rest": mb}},
            seen_sources=("nadig_hepg2",), unseen_sources=("cd4t_rest",),
        )
        self.assertFalse(missing_orth_key["matched_arm_valid"])
        self.assertEqual(missing_orth_key["line_summary"], "INCONCLUSIVE")

        # Degenerate Gate-2 Group-2 safeguard (`|G_2_qual| < 2`)
        grp_ok = adjudicate_group_r16(
            per_line_h_i1={
                "nadig_hepg2": {"line_summary": "LINE_EQUIVALENT_NONZERO_SIGNAL"},
                "a549_ko": {"line_summary": "LINE_MATCHED_SPECIFIC"},
                "cd4t_rest": {"line_summary": "LINE_EQUIVALENT_NONZERO_SIGNAL"},
            },
            per_line_h_i2={
                "nadig_hepg2": {"classification": "FCT_SUPERIOR_TO_MODEL"},
                "a549_ko": {"classification": "FCT_SUPERIOR_TO_MODEL"},
                "cd4t_rest": {"classification": "INDISTINGUISHABLE"},
            },
        )
        self.assertFalse(grp_ok["group2_downgraded_single_stratum"])
        self.assertEqual(grp_ok["group2_crispri"]["h_i1_counts"]["LINE_EQUIVALENT_NONZERO_SIGNAL"], 1)
        self.assertEqual(grp_ok["group2_crisprko"]["h_i1_counts"]["LINE_MATCHED_SPECIFIC"], 1)

        pack = export_canonical_parity_pack(output_json=None, n_boot=300, seed=DEFAULT_R16_SEED)
        sc = pack["scenarios"]
        self.assertEqual(sc["N300_matched_better_strong"]["h_i1_line_summary"], "LINE_MATCHED_SPECIFIC")
        self.assertFalse(sc["N300_matched_better_strong"]["h_i1_all_within_delta"])
        self.assertEqual(sc["N400_matched_better_sub_delta"]["h_i1_line_summary"], "LINE_MATCHED_SPECIFIC_SUB_DELTA")
        self.assertTrue(sc["N400_matched_better_sub_delta"]["h_i1_all_within_delta"])
        self.assertEqual(sc["N200_nonzero_signal_equivalent"]["h_i1_line_summary"], "LINE_EQUIVALENT_NONZERO_SIGNAL")
        self.assertEqual(sc["N200_null_signal_equivalent"]["h_i1_line_summary"], "LINE_EQUIVALENT_NULL_SIGNAL")
        self.assertEqual(sc["N120_mismatched_better"]["h_i1_line_summary"], "LINE_MISMATCHED_BETTER")
        self.assertEqual(sc["N120_stratum_split"]["h_i1_line_summary"], "MIXED_OR_INCONCLUSIVE")
        self.assertTrue(pack["degenerate_group2_check"]["group2_downgraded_single_stratum"])

        # Varying N TOST precision check: same noisy DGP passes EQUIVALENT at N=300 yet loses TOST power at N=16
        self.assertEqual(sc["N300_noisy_wellpowered_equiv"]["h_i1_line_summary"], "LINE_EQUIVALENT_NONZERO_SIGNAL")
        self.assertEqual(sc["N016_underpowered_equiv"]["h_i1_line_summary"], "MIXED_OR_INCONCLUSIVE")
        self.assertGreater(sc["N016_underpowered_equiv"]["h_i1_counts"]["INCONCLUSIVE"], 0)
        ci_w_300 = max(
            abs(v["ci95"][1] - v["ci95"][0])
            for v in sc["N300_noisy_wellpowered_equiv"]["h_i1_per_source_est"].values()
        )
        ci_w_16 = max(
            abs(v["ci95"][1] - v["ci95"][0])
            for v in sc["N016_underpowered_equiv"]["h_i1_per_source_est"].values()
        )
        self.assertGreater(ci_w_16, ci_w_300)

        # H-I2 negative branch (`FCT_SUPERIOR_TO_MODEL`) + M2 weight shift check
        self.assertEqual(
            adjudicate_claim_i2_line({"estimate": -0.1, "ci975_lo": -0.15, "ci975_hi": -0.05, "valid_frac": 1.0})["classification"],
            "FCT_SUPERIOR_TO_MODEL",
        )
        wc = compute_winner_curse_weight_shift_r16(
            ctrl_mean_libnorm=np.linspace(0.05, 8.0, 120),
            delta_obs=np.random.default_rng(7).normal(0.0, 0.3, size=(30, 120)),
            target_gene_idx=np.arange(30),
            n_boot=100,
        )
        self.assertTrue(np.isfinite(wc["estimate"]))

    def test_5_blinding_and_training_leakage_guards(self) -> None:
        """Freeze & leakage gates must block pre-freeze real execution and training-line leakage."""
        assert_zero_training_leakage("STATE_SE", "nadig_hepg2", "zeroshot/hepg2")
        assert_zero_training_leakage("STATE_HVG", "a549_ko", "zeroshot/ensemble4")
        with self.assertRaises(TrainingLeakageError):
            assert_zero_training_leakage("STATE_SE", "nadig_jurkat", "zeroshot/hepg2")
        with self.assertRaises(TrainingLeakageError):
            assert_zero_training_leakage("STATE_HVG", "replogle_rpe1", "zeroshot/k562")

        with tempfile.TemporaryDirectory() as td:
            sdir = Path(td)
            with self.assertRaises(BlindingGateError):
                assert_freeze_gates_r16(sdir)
            with self.assertRaises(BlindingGateError):
                runner_main(["--phase", "real_r16", "--status-dir", str(sdir), "--output-dir", str(sdir / "out")])

            (sdir / "PREREG_I_COMMIT.txt").write_text("LOCK\n", encoding="utf-8")
            with self.assertRaises(BlindingGateError):
                assert_freeze_gates_r16(sdir)

            (sdir / "CODE_FREEZE_R16.txt").write_text("FREEZE\n", encoding="utf-8")
            # Missing R16_GO.txt MUST still fail cleanly
            with self.assertRaises(BlindingGateError):
                assert_freeze_gates_r16(sdir)

            # Empty R16_GO.txt MUST fail cleanly
            (sdir / "R16_GO.txt").write_text("   \n", encoding="utf-8")
            with self.assertRaises(BlindingGateError):
                assert_freeze_gates_r16(sdir)

            # R16_GO.txt predating R16_GATE2_SCHEMA.md MUST fail cleanly
            gate2_file = sdir / "R16_GATE2_SCHEMA.md"
            go_file = sdir / "R16_GO.txt"
            go_file.write_text("GO\n", encoding="utf-8")
            gate2_file.write_text("# schema\n", encoding="utf-8")
            os.utime(go_file, (1_000_000.0, 1_000_000.0))
            os.utime(gate2_file, (2_000_000.0, 2_000_000.0))
            with self.assertRaises(BlindingGateError):
                assert_freeze_gates_r16(sdir)

            # Fresh R16_GO.txt after R16_GATE2_SCHEMA.md passes all 3 gates
            os.utime(go_file, (3_000_000.0, 3_000_000.0))
            receipt = assert_freeze_gates_r16(sdir)
            self.assertEqual(receipt["r16_go"], "GO")

            # Adapter with allow_synthetic_smoke=False still blocks even synthetic call pre-freeze
            adapter_strict = CleanRoomArcStateAdapter(
                model_key="STATE_HVG", checkpoint_id="zeroshot/hepg2", status_dir=sdir
            )
            go_file.unlink()
            with self.assertRaises(BlindingGateError):
                adapter_strict.predict_arm_deltas(
                    "nadig_hepg2", "nadig_hepg2", np.ones((200, 10)), [f"G{i}" for i in range(10)], ["AAAS"],
                    is_synthetic_data=True,
                )

    def test_6_st_hvg_real_weights_end_to_end_synthetic_smoke(self) -> None:
        """
        Execute real `ST-HVG-Replogle` (`zeroshot/hepg2`) weights strictly on randomly generated
        synthetic counts across matched and batch-imbalanced mismatched synthetic control pools,
        then evaluate end-to-end through `r_orth_cal`, `r_dorth_cal`, paired bootstrap, TOST & `H-I1`/`H-I2`.
        Verifies control sensitivity under BOTH `bypass_zero` and `index_zero` neutral covariates.
        """
        adapter_bypass = CleanRoomArcStateAdapter(
            model_key="STATE_HVG",
            checkpoint_id="zeroshot/hepg2",
            device="cpu",
            neutral_batch_mode="bypass_zero",
            delta_mode="minus_pred_ntc",
            allow_synthetic_smoke=True,
            batch_size_sets=16,
        )
        adapter_bypass.load_weights()
        model_genes = adapter_bypass._bundles["hepg2"]["gene_names"]
        pert_lookup = adapter_bypass._bundles["hepg2"]["pert_lookup"]
        eval_perts = [p for p in sorted(pert_lookup.keys()) if p != "non-targeting"][:25]
        G = len(model_genes)
        P = len(eval_perts)

        rng = np.random.default_rng(DEFAULT_R16_SEED)
        # Synthesize raw integer counts for matched and 2 imbalanced mismatched synthetic cohorts (zero real cells!)
        matched_batches = (["B1"] * 420) + (["B2"] * 180) + (["B3"] * 100)
        matched_ctrl_raw = rng.poisson(1.5, size=(len(matched_batches), G)).astype(np.float32)

        mis_seen_raw = rng.poisson(6.0, size=(650, G)).astype(np.float32)
        mis_seen_batches = (["S1"] * 500) + (["S2"] * 150)
        mis_unseen_raw = rng.poisson(0.4, size=(580, G)).astype(np.float32)
        mis_unseen_batches = (["U1"] * 450) + (["U2"] * 130)

        out_bypass = adapter_bypass.run_matched_and_mismatched_arms(
            target_line="nadig_hepg2",
            matched_ctrl_matrix=matched_ctrl_raw,
            matched_genes=model_genes,
            mismatched_ctrls={
                "nadig_jurkat": (mis_seen_raw, model_genes, mis_seen_batches),
                "cd4t_rest": (mis_unseen_raw, model_genes, mis_unseen_batches),
            },
            eval_perts=eval_perts,
            matched_batch_labels=matched_batches,
            input_is_log1p=False,
            seed=DEFAULT_R16_SEED,
            is_synthetic_data=True,
        )
        d_m = out_bypass["matched"]["delta_pred"]
        d_s1 = out_bypass["mismatched"]["nadig_jurkat"]["delta_pred"]
        d_s2 = out_bypass["mismatched"]["cd4t_rest"]["delta_pred"]
        self.assertEqual(d_m.shape, (P, G))
        self.assertTrue(np.all(np.isfinite(d_m)))
        self.assertEqual(out_bypass["matched"]["sampling_audit"]["sets"].shape, (16, 32))
        self.assertEqual(out_bypass["matched"]["gene_mapping"]["n_mapped_genes"], 2000)
        # Sensitivity proof (`PREREG_I v4 §2.1 M5.3`): changing X_ctrl changes predictions cleanly!
        self.assertGreater(float(np.linalg.norm(d_m - d_s1)), 1e-2)
        self.assertGreater(float(np.linalg.norm(d_m - d_s2)), 1e-2)

        # Verify `index_zero` neutral covariate mode also works cleanly and shows X_ctrl sensitivity
        adapter_idx0 = CleanRoomArcStateAdapter(
            model_key="STATE_HVG",
            checkpoint_id="zeroshot/hepg2",
            device="cpu",
            neutral_batch_mode="index_zero",
            delta_mode="minus_pred_ntc",
            allow_synthetic_smoke=True,
            batch_size_sets=16,
        )
        res_idx0_m = adapter_idx0.predict_arm_deltas(
            "nadig_hepg2", "nadig_hepg2", matched_ctrl_raw, model_genes, eval_perts,
            batch_labels=matched_batches, input_is_log1p=False, is_synthetic_data=True,
        )
        res_idx0_s = adapter_idx0.predict_arm_deltas(
            "nadig_hepg2", "cd4t_rest", mis_unseen_raw, model_genes, eval_perts,
            batch_labels=mis_unseen_batches, input_is_log1p=False, is_synthetic_data=True,
        )
        self.assertGreater(float(np.linalg.norm(res_idx0_m["delta_pred"] - res_idx0_s["delta_pred"])), 1e-2)

        # Feed real ST-HVG outputs into clean-room evaluation & two-stratum adjudicator
        t_mean = rng.normal(0.0, 0.15, size=(P, G))
        grand_mean = rng.normal(0.0, 0.25, size=G)
        obs = d_m + rng.normal(0.0, 0.10, size=(P, G))
        ha = obs + rng.normal(0.0, 0.04, size=(P, G))
        hb = obs + rng.normal(0.0, 0.04, size=(P, G))
        target_cols = np.arange(P, dtype=np.int64)
        hv = np.ones(P, dtype=bool)
        pred_fct = t_mean + 0.3 * (obs - t_mean)

        fx_pipe = {
            "regime": "st_hvg_synthetic_smoke",
            "n_perts": P,
            "n_genes": G,
            "seen_sources": ("nadig_jurkat",),
            "unseen_sources": ("cd4t_rest",),
            "delta_obs": obs,
            "half_a_obs": ha,
            "half_b_obs": hb,
            "half_valid": hv,
            "t_mean": t_mean,
            "grand_mean": grand_mean,
            "target_cols": target_cols,
            "pred_matched": d_m,
            "pred_fct": pred_fct,
            "mis_preds": {"nadig_jurkat": d_s1, "cd4t_rest": d_s2},
        }
        ev = evaluate_synthetic_fixture_r16(fx_pipe, n_boot=200, seed=DEFAULT_R16_SEED)
        self.assertIn(
            ev["h_i1"]["line_summary"],
            {
                "LINE_MATCHED_SPECIFIC",
                "LINE_MATCHED_SPECIFIC_SUB_DELTA",
                "LINE_EQUIVALENT_NONZERO_SIGNAL",
                "LINE_EQUIVALENT_NULL_SIGNAL",
                "LINE_MISMATCHED_BETTER",
                "MIXED_OR_INCONCLUSIVE",
                "INCONCLUSIVE",
            },
        )
        self.assertIn(
            ev["h_i2"]["classification"],
            {"MODEL_SUPERIOR_TO_FCT", "FCT_SUPERIOR_TO_MODEL", "INDISTINGUISHABLE"},
        )
        self.assertTrue(np.isfinite(ev["dorth_matched"]["estimate"]))

    def test_7_gate2_hgnc_hvg_se_forward_and_npz_runner(self) -> None:
        """
        Verify Gate-2 structural contracts on purely synthetic inputs:
          - Frozen R10 `HGNCBridge` canonical alignment across ST-HVG, P_model, & SE-600M PE keys
          - `official_median_log1p` vs `libnorm1e4_log1p` preprocessing
          - `ST-HVG` reads direct `clip(preds, 0, 14)` and differs from untrained `gene_decoder(preds)`
          - `padded=True` matches manual 32-cell `padded=False` windows, and `pre_sampled=True` is exact
          - Lightweight synthetic `StateEmbeddingModel` forward pipeline
          - End-to-end synthetic runner `.npz` export schema & finite values
        """
        import torch
        from omegaconf import OmegaConf
        from state.emb.nn.model import StateEmbeddingModel

        hb = get_hgnc_bridge()
        # 1. Verify HVG preprocessing modes on synthetic counts
        raw_syn = np.array([[10.0, 20.0, 0.0], [5.0, 5.0, 10.0], [100.0, 0.0, 0.0]], dtype=np.float32)
        p_1e4 = preprocess_hvg_expression(raw_syn, mode="libnorm1e4_log1p")
        p_med = preprocess_hvg_expression(raw_syn, mode="official_median_log1p")
        # Row totals [30, 20, 100] -> median 30 -> row 0 unchanged before log1p
        np.testing.assert_allclose(p_med[0], np.log1p(raw_syn[0]), atol=1e-5)
        self.assertGreater(float(np.max(np.abs(p_1e4 - p_med))), 0.1)

        # 2. Verify SE-600M gene mapping: PREREG_I v5 §2.1.4 L81 default exact lookup (use_hgnc_rescue=False) vs optional rescue
        pe_fake = ["A1BG", "AARS", "DDI2", "RSC1A1", "MYO18A", "TIAF1", "TAZ"]
        src_g = ["AARS1", "DDI2", "TIAF1", "MISSING_GENE"]
        src_exact, pe_exact, m_exact = map_genes_to_se600m_pe(src_g, pe_fake, hgnc_bridge=hb)
        self.assertEqual(m_exact["mapping_mode"], "exact_string_lookup")
        self.assertEqual(m_exact["n_hgnc_canon_rescued"], 0)
        self.assertEqual(m_exact["n_raw_exact"], 2)
        self.assertEqual(pe_exact.tolist(), [2, 5])
        self.assertEqual(src_exact.tolist(), [1, 2])

        src_c, pe_g, mstats = map_genes_to_se600m_pe(src_g, pe_fake, hgnc_bridge=hb, use_hgnc_rescue=True)
        self.assertEqual(mstats["n_hgnc_canon_rescued"], 1)  # AARS1 -> AARS when explicit rescue enabled
        self.assertEqual(mstats["n_raw_exact"], 2)
        self.assertEqual(pe_g.tolist(), [1, 2, 5])
        self.assertEqual(src_c.tolist(), [0, 1, 2])

        # Verify Adapter & preprocess defaults follow PREREG_I v5 primary `official_median_log1p`
        self.assertEqual(
            CleanRoomArcStateAdapter(model_key="STATE_HVG", checkpoint_id="zeroshot/hepg2", allow_synthetic_smoke=True).hvg_preproc_mode,
            "official_median_log1p",
        )
        np.testing.assert_allclose(preprocess_hvg_expression(raw_syn)[0], p_med[0], atol=1e-6)

        # 3. Verify synthetic StateEmbeddingModel forward pipeline matches Gate 2 [emb || ds_emb]
        cfg = OmegaConf.create({
            "dataset": {
                "current": "syn",
                "syn": {"num_datasets": 3},
                "pad_length": 16,
                "cls_token_idx": 3,
                "P": 4,
                "N": 4,
                "S": 4,
            },
            "model": {"counts": True, "rda": True, "dataset_correction": True, "d_hid": 16, "emsize": 8, "nhead": 2, "nlayers": 1, "output_dim": 8, "dropout": 0.0, "batch_size": 4},
            "embeddings": {"current": "e1", "e1": {"num": 16, "size": 16, "valid_genes_masks": None}},
            "tokenizer": {"token_dim": 16},
            "loss": {"name": "tabular"},
        })
        se_tiny = StateEmbeddingModel(
            token_dim=16, d_model=8, nhead=2, d_hid=16, nlayers=1, output_dim=8, dropout=0.0, cfg=cfg
        )
        pe_tensor = torch.randn(16, 16, dtype=torch.float32)
        se_tiny.pe_embedding = torch.nn.Embedding.from_pretrained(pe_tensor)
        syn_counts = np.random.default_rng(7).poisson(4.0, size=(8, 7)).astype(np.float32)
        syn_counts[:, -1] = 0.0
        emb_tiny, st_tiny = encode_with_loaded_se_model(
            se_tiny, pe_keys=pe_fake + [f"G{i}" for i in range(9)],
            ctrl_counts=syn_counts, source_genes=pe_fake, hgnc_bridge=hb, use_hgnc_rescue=False,
            pad_length=16, cls_token_idx=3, batch_size=4, seed=DEFAULT_R16_SEED, device=torch.device("cpu"), cfg=cfg,
        )
        self.assertEqual(emb_tiny.shape, (8, 18))  # 8 emb + 10 ds_emb
        self.assertEqual(st_tiny["n_mapped_genes"], 7)
        self.assertEqual(st_tiny["collator"], "official_VCIDatasetSentenceCollator_eval")
        self.assertEqual(len(st_tiny["emb_sha256"]), 64)
        np.testing.assert_allclose(np.linalg.norm(emb_tiny[:, :8], axis=1), np.ones(8), atol=1e-5)

        # 4. ST-HVG direct preds vs untrained gene_decoder, padded=True parity, and pre_sampled=True parity
        adapter = CleanRoomArcStateAdapter(
            model_key="STATE_HVG", checkpoint_id="zeroshot/hepg2", device="cpu",
            neutral_batch_mode="index_zero", delta_mode="minus_pred_ntc",
            hvg_preproc_mode="libnorm1e4_log1p", allow_synthetic_smoke=True, hgnc_bridge=hb,
        )
        adapter.load_weights()
        b = adapter._bundles["hepg2"]
        m_genes = b["gene_names"]
        rng = np.random.default_rng(123)
        raw_x = rng.poisson(2.0, size=(128, len(m_genes))).astype(np.float32)
        samp = sample_stratified_control_sets_r16(
            n_cells=128, n_ctrl=128, cell_set_len=32, min_ctrl_floor=128, seed=DEFAULT_R16_SEED,
            cell_barcodes=[f"BC_{i}" for i in range(128)],
        )
        res_full = adapter.predict_arm_deltas(
            "nadig_hepg2", "nadig_hepg2", raw_x, m_genes, ["AAAS", "AARS1"],
            n_ctrl=128, cell_set_len=32, is_synthetic_data=True, hgnc_bridge=hb,
            cell_barcodes=[f"BC_{i}" for i in range(128)],
        )
        res_presamp = adapter.predict_arm_deltas(
            "nadig_hepg2", "nadig_hepg2", raw_x[samp["selected_indices"]], m_genes, ["AAAS", "AARS1"],
            n_ctrl=128, cell_set_len=32, is_synthetic_data=True, pre_sampled=True,
            precomputed_sampling_audit=samp, hgnc_bridge=hb,
        )
        np.testing.assert_allclose(res_full["delta_pred"], res_presamp["delta_pred"], atol=1e-7)
        self.assertEqual(res_full["selected_barcodes"], samp["selected_barcodes"])

        x_in = preprocess_hvg_expression(raw_x[samp["selected_indices"]], mode="libnorm1e4_log1p")
        x_t = torch.tensor(x_in, dtype=torch.float32)
        pvec = b["pert_lookup"]["AAAS"].unsqueeze(0).expand(128, -1)
        b0 = torch.zeros(128, dtype=torch.long)
        with torch.no_grad():
            step_pad = b["model"].predict_step({"ctrl_cell_emb": x_t, "pert_emb": pvec, "batch": b0}, batch_idx=0, padded=True)
            win_outs = [
                b["model"].predict_step(
                    {"ctrl_cell_emb": x_t[k * 32 : (k + 1) * 32], "pert_emb": pvec[k * 32 : (k + 1) * 32], "batch": b0[k * 32 : (k + 1) * 32]},
                    batch_idx=0, padded=False,
                )["preds"]
                for k in range(4)
            ]
        cat_preds = torch.cat(win_outs, dim=0).numpy()
        np.testing.assert_allclose(step_pad["preds"].numpy(), cat_preds, atol=1e-6)
        self.assertGreater(float(np.max(np.abs(step_pad["preds"].numpy() - step_pad["pert_cell_counts_preds"].numpy()))), 0.1)

        # 5. End-to-end synthetic runner + .npz schema + barcode + Perm5 check
        with tempfile.TemporaryDirectory() as td:
            res_run = run_synthetic_runner_smoke(output_dir=Path(td), device="cpu", n_sets=4, n_perts=6, include_real_se600m=False)
            npz_f = Path(res_run["npz_path"])
            self.assertTrue(npz_f.exists())
            self.assertEqual(res_run["perm5_shape"], (5, 6, 2000))
            self.assertEqual(res_run["matched_barcodes_len"], 128)
            self.assertGreater(res_run["hgnc_rescued_hvg"], 0)
            with np.load(npz_f, allow_pickle=False) as z:
                for req in (
                    "model_key", "checkpoint_id", "target_line", "neutral_batch_mode",
                    "eval_genes", "eval_perts", "target_gene_idx",
                    "delta_obs", "delta_half_a", "delta_half_b", "half_valid",
                    "delta_transfer", "delta_fct", "arm_names", "arm_strata",
                    "delta_pred_arms", "delta_pred_matched", "delta_pred_resample",
                    "matched_barcodes", "resample_barcodes",
                    "delta_pred_mismatched__nadig_jurkat", "mismatched_barcodes__nadig_jurkat",
                    "delta_pred_mismatched__dld1", "mismatched_barcodes__dld1",
                ):
                    self.assertIn(req, z.files)
                self.assertEqual(z["delta_pred_matched"].shape, (6, 240))
                self.assertEqual(len(z["matched_barcodes"]), 128)
                self.assertEqual(str(z["neutral_batch_mode"]), "index_zero")
                self.assertEqual(str(z["hvg_preproc_mode"]), "official_median_log1p")

    def test_8_real_se600m_safetensors_mps_fp32_to_st_se_synthetic_smoke(self) -> None:
        """
        End-to-end v5 verification using REAL `SE-600M` (`model.safetensors`, Apple `MPS fp32`
        or CPU fallback) + REAL `ST-SE-Replogle` (`zeroshot/hepg2`, decoded via `clip(gene_decoder(preds), 0, 14)`)
        strictly on 32 (`1 x 32`) synthetic Poisson raw integer counts (ZERO real cells used).
        Verifies exact string lookup (`use_hgnc_rescue=False`), both `2048` & `2058` centroid metrics, & MPS fp32 guard.
        """
        import torch

        hb = get_hgnc_bridge()
        dev_str = "mps" if torch.backends.mps.is_available() else "cpu"
        se_enc = CleanRoomSE600MEncoder(device=dev_str, batch_size=16, hgnc_bridge=hb)
        se_enc.load()
        self.assertEqual(len(se_enc._pe_keys), 19790)

        rng = np.random.default_rng(DEFAULT_R16_SEED)
        pe_sub = se_enc._pe_keys[:2400]
        raw_a = rng.poisson(2.5, size=(32, len(pe_sub))).astype(np.float32)
        raw_b = rng.poisson(4.5, size=(32, len(pe_sub))).astype(np.float32)

        emb_a, st_a = se_enc.encode_cells(raw_a, pe_sub, seed=DEFAULT_R16_SEED, use_hgnc_rescue=False)
        self.assertEqual(emb_a.shape, (32, 2058))
        self.assertEqual(st_a["mapping_mode"], "exact_string_lookup")
        self.assertEqual(st_a["collator"], "official_VCIDatasetSentenceCollator_eval")
        self.assertEqual(st_a["batch_size"], 16)
        self.assertEqual(len(st_a["emb_sha256"]), 64)
        self.assertEqual(st_a["n_hgnc_canon_rescued"], 0)
        self.assertEqual(st_a["n_mapped_genes"], 2400)
        self.assertTrue(np.all(np.isfinite(emb_a)))
        np.testing.assert_allclose(np.linalg.norm(emb_a[:, :2048], axis=1), np.ones(32), atol=1e-4)

        if dev_str == "mps":
            se_enc._model.half()
            try:
                with self.assertRaises(ValueError):
                    se_enc.encode_cells(raw_a[:2], pe_sub, seed=DEFAULT_R16_SEED)
            finally:
                se_enc._model.float()

        adapter_se = CleanRoomArcStateAdapter(
            model_key="STATE_SE",
            checkpoint_id="zeroshot/hepg2",
            device=dev_str,
            neutral_batch_mode="index_zero",
            delta_mode="minus_pred_ntc",
            allow_synthetic_smoke=True,
            se600m_encoder=se_enc,
            hgnc_bridge=hb,
            se600m_use_hgnc_rescue=False,
        )
        perts = ["AAAS", "DARS1", "H2AZ1", "TAFAZZIN"]
        arm_res = adapter_se.run_matched_and_mismatched_arms(
            target_line="nadig_hepg2",
            matched_ctrl_matrix=None,
            matched_genes=pe_sub,
            mismatched_ctrls={"cd4t_rest": (raw_b, pe_sub, ["U1"] * 32)},
            eval_perts=perts,
            matched_batch_labels=["B1"] * 32,
            input_is_log1p=False,
            seed=DEFAULT_R16_SEED,
            n_ctrl=32,
            cell_set_len=32,
            is_synthetic_data=True,
            precomputed_embeddings={
                "nadig_hepg2": (emb_a, {"n_sampled": 32, "n_sets": 1, "sampling_mode": "pre_sampled", "selected_barcodes": [f"M_{i}" for i in range(32)]}, st_a)
            },
            hgnc_bridge=hb,
            pre_sampled=True,
            mismatched_barcodes={"cd4t_rest": [f"S_{i}" for i in range(32)]},
        )
        d_m = arm_res["matched"]["delta_pred"]
        d_s = arm_res["mismatched"]["cd4t_rest"]["delta_pred"]
        self.assertEqual(d_m.shape, (4, 2000))
        self.assertTrue(np.all(np.isfinite(d_m)))
        self.assertGreater(float(np.linalg.norm(d_m - d_s)), 1e-2)
        cdiag = arm_res["control_centroid_diagnostics"]["cd4t_rest"]
        self.assertIn("dim_2048", cdiag)
        self.assertIn("dim_2058", cdiag)
        self.assertIn("cosine", cdiag["dim_2048"])
        self.assertIn("centered_cosine", cdiag["dim_2048"])

    def test_9_addendum_a_exact_invariance_resample_and_sub_delta_rules(self) -> None:
        """
        Verify Addendum A Items 8–10 strictly:
          - Item 10: `max_pred_abs_diff <= 1e-12` (`ci975=[0,0]`) -> `EXACT_INVARIANCE` (pair)
            and counts as `EQUIVALENT` in line-level adjudication (`LINE_EQUIVALENT_NONZERO_SIGNAL`),
            whereas `max_abs_diff > 1e-12` with `ci975_lo==0` remains `INCONCLUSIVE`.
          - Items 8–9: non-overlapping resample sampling (`>=1024`, `128..511`, `<128`), exact 0 overlap,
            `compute_paired_diff_vs_resample_r16` interpretation strings, & `SUB_DELTA` interpretation note.
        """
        fx = make_synthetic_pert_cohort(n_perts=120, n_genes=90, regime="matched_better_strong", seed=20261610)
        obs = fx["delta_obs"]
        tm = fx["t_mean"]
        ha, hb, hv = fx["half_a_obs"], fx["half_b_obs"], fx["half_valid"]
        tc = fx["target_cols"]

        pred_ident = fx["pred_matched"].copy()
        ev_ident = evaluate_and_adjudicate_line(
            {
                "eval_genes": [f"G{i}" for i in range(90)],
                "eval_perts": [f"P{i}" for i in range(120)],
                "target_gene_idx": tc,
                "delta_obs": obs,
                "delta_half_a": ha,
                "delta_half_b": hb,
                "half_valid": hv,
                "delta_transfer": tm,
                "source_grand_mean": fx["grand_mean"],
                "delta_fct": fx["pred_fct"],
                "ctrl_mean_libnorm": np.ones(90),
            },
            delta_pred_matched_eval=pred_ident,
            delta_pred_mis_eval={
                "nadig_jurkat": pred_ident + 5e-13,
                "dld1": pred_ident.copy(),
            },
            seen_sources=("nadig_jurkat",),
            unseen_sources=("dld1",),
            n_boot=200,
            seed=DEFAULT_R16_SEED,
        )
        self.assertEqual(ev_ident["h_i1"]["per_source"]["nadig_jurkat"]["category"], "EXACT_INVARIANCE")
        self.assertEqual(ev_ident["h_i1"]["per_source"]["dld1"]["category"], "EXACT_INVARIANCE")
        self.assertEqual(ev_ident["h_i1"]["line_summary"], "LINE_EQUIVALENT_NONZERO_SIGNAL")
        self.assertEqual(ev_ident["h_i1"]["counts"]["EXACT_INVARIANCE"], 2)

        non_inv_tie = {
            "estimate": 0.01, "ci975_lo": 0.0, "ci975_hi": 0.02,
            "ci95_lo": 0.002, "ci95_hi": 0.018, "valid_frac": 1.0,
            "max_pred_abs_diff": 1e-4, "exact_invariance": False,
        }
        self.assertEqual(classify_pair_tost_r16(non_inv_tie, DEFAULT_TOST_DELTA)["category"], "INCONCLUSIVE")

        for total_n, exp_mode, exp_sec_n in ((1200, "stratified_by_batch", 512), (720, "small_ntc_truncated", 192), (600, "excluded_below_floor", 0)):
            bcs = [f"CELL_{i:04d}" for i in range(total_n)]
            bl = [f"B_{i % 3}" for i in range(total_n)]
            p_s = sample_stratified_control_sets_r16(n_cells=total_n, batch_labels=bl, seed=DEFAULT_R16_SEED, cell_barcodes=bcs)
            r_s = sample_nonoverlapping_resample_sets_r16(p_s, n_cells=total_n, batch_labels=bl, resample_seed=RESAMPLE_R16_SEED, cell_barcodes=bcs)
            self.assertEqual(r_s["sampling_mode"], exp_mode)
            self.assertEqual(r_s["n_sampled"], exp_sec_n)
            self.assertEqual(r_s["n_overlap_with_primary"], 0)
            if not r_s["excluded"]:
                self.assertTrue(set(p_s["selected_barcodes"]).isdisjoint(set(r_s["selected_barcodes"])))

        # Sub-delta scenario where resample loss exceeds/matches mis loss -> locked non-specificity text + Item 9 text
        fx_sub = make_synthetic_pert_cohort(n_perts=350, n_genes=110, regime="matched_better_sub_delta", seed=20261611)
        rng_j = np.random.default_rng(777)
        resamp_degraded = fx_sub["pred_matched"] + rng_j.normal(0.0, 0.25, size=fx_sub["pred_matched"].shape)
        ev_jit = evaluate_and_adjudicate_line(
            {
                "eval_genes": [f"G{i}" for i in range(110)],
                "eval_perts": [f"P{i}" for i in range(350)],
                "target_gene_idx": fx_sub["target_cols"],
                "delta_obs": fx_sub["delta_obs"],
                "delta_half_a": fx_sub["half_a_obs"],
                "delta_half_b": fx_sub["half_b_obs"],
                "half_valid": fx_sub["half_valid"],
                "delta_transfer": fx_sub["t_mean"],
                "source_grand_mean": fx_sub["grand_mean"],
                "delta_fct": fx_sub["pred_fct"],
                "ctrl_mean_libnorm": np.ones(110),
            },
            delta_pred_matched_eval=fx_sub["pred_matched"],
            delta_pred_mis_eval=fx_sub["mis_preds"],
            seen_sources=fx_sub["seen_sources"],
            unseen_sources=fx_sub["unseen_sources"],
            n_boot=250,
            seed=DEFAULT_R16_SEED,
            delta_pred_resample_eval=resamp_degraded,
        )
        self.assertEqual(ev_jit["h_i1"]["line_summary"], "LINE_MATCHED_SPECIFIC_SUB_DELTA")
        self.assertIn("δ=0.04", ev_jit["h_i1"]["interpretation_note"])
        self.assertFalse(ev_jit["resample_orth"]["all_sources_exceed_resample"])
        self.assertEqual(ev_jit["resample_orth"]["overall_interpretation"], "换成其他细胞系的对照，影响不超过重抽本细胞系的对照")

        # True line specificity: mis is much worse than near-identical resample -> exceeds_resample_loss = True!
        rng = np.random.default_rng(778)
        mis_bad = {s: fx_sub["t_mean"] + rng.normal(0.0, 0.30, size=fx_sub["pred_matched"].shape) for s in fx_sub["seen_sources"] + fx_sub["unseen_sources"]}
        res_vs_bad = compute_paired_diff_vs_resample_r16(
            fx_sub["delta_obs"], fx_sub["pred_matched"], fx_sub["pred_matched"] + rng.normal(0.0, 0.01, size=fx_sub["pred_matched"].shape),
            mis_bad, fx_sub["t_mean"], fx_sub["target_cols"], delta_half_a=fx_sub["half_a_obs"],
            delta_half_b=fx_sub["half_b_obs"], half_valid=fx_sub["half_valid"], n_boot=200, seed=DEFAULT_R16_SEED,
        )
        self.assertTrue(res_vs_bad["all_sources_exceed_resample"])
        self.assertEqual(res_vs_bad["overall_interpretation"], "模型用到了目标线特有的对照信息")

        # Mixed resample outcome: one mismatched source worse than resample, one equal -> mixed per-source narrative!
        mis_mixed = {
            "nadig_jurkat": fx_sub["pred_matched"] + rng.normal(0.0, 0.002, size=fx_sub["pred_matched"].shape),
            "dld1": fx_sub["t_mean"] + rng.normal(0.0, 0.30, size=fx_sub["pred_matched"].shape),
        }
        res_mixed = compute_paired_diff_vs_resample_r16(
            fx_sub["delta_obs"], fx_sub["pred_matched"], fx_sub["pred_matched"] + rng.normal(0.0, 0.01, size=fx_sub["pred_matched"].shape),
            mis_mixed, fx_sub["t_mean"], fx_sub["target_cols"], delta_half_a=fx_sub["half_a_obs"],
            delta_half_b=fx_sub["half_b_obs"], half_valid=fx_sub["half_valid"], n_boot=200, seed=DEFAULT_R16_SEED,
        )
        self.assertEqual(res_mixed["overall_status"], "MIXED_SOURCE_BY_SOURCE")
        self.assertIn("不得作全细胞系特异性结论", res_mixed["overall_interpretation"])

    def test_10_gse295214_synthetic_h5_reader_qc_counts_and_fct_same_targets(self) -> None:
        """
        Test Group 3 (`PREREG_I v5 §2.4` + `Addendum A`) strictly on synthetic CellRanger 10x `.h5`
        fixtures: integer validation, pooled bimodal threshold, `non-targeting_*`, `max_pct_mt=None`,
        lexicographical odd/even split-half, non-overlapping resample + barcodes, `R16_GSE295214_QC_COUNTS.md`,
        `r9_6 \\ {vcc2025_h1}` primary FCT, and fixed-target `full r9_6` secondary sensitivity.
        """
        from data_loader_r10 import load_r9_common_genes

        hb = get_hgnc_bridge()
        hgnc_l = build_hgnc_lookup()
        r9_all = load_r9_common_genes(hgnc_bridge=hb)
        r9_genes = r9_all[:235]
        for req_g in ("TP53", "MYC", "AAAS", "AAMP", "AAR2"):
            if req_g not in r9_genes:
                r9_genes.append(req_g)
        rna_names = ["MT-CO1"] + r9_genes
        rna_ids = [f"ENSG{i:011d}" for i in range(len(rna_names))]
        guide_names = [
            "non-targeting_0001", "non-targeting_0002",
            "TP53_g1", "MYC_g1", "AAAS_g1", "safe-targeting_01", "intergenic_01",
        ]
        # Verify Addendum A §2.1 exact guide parsing rules directly before building mock 10x HDF5s
        from audit_r16_raw_controls import parse_gse295214_guide_r16

        self.assertEqual(parse_gse295214_guide_r16("non-targeting_1|non-targeting_2", "", hgnc_l)[0], "ntc")
        self.assertEqual(parse_gse295214_guide_r16("TP53_1|TP53_P1P2_2", "", hgnc_l)[1][1], "TP53")
        self.assertEqual(parse_gse295214_guide_r16("TP53_1|MYC_1", "", hgnc_l)[0], "exclude")
        self.assertEqual(parse_gse295214_guide_r16("TP53_MYC_1", "", hgnc_l)[0], "exclude")
        self.assertEqual(parse_gse295214_guide_r16("HLA-A_P1P2_1", "", hgnc_l)[1][1], "HLA-A")
        self.assertEqual(parse_gse295214_guide_r16("MYC-sg1", "", hgnc_l)[1][1], "MYC")
        self.assertEqual(parse_gse295214_guide_r16("gene_desert_1", "", hgnc_l)[0], "exclude")
        self.assertEqual(parse_gse295214_guide_r16("TP53-and-MYC", "", hgnc_l)[0], "exclude")
        guide_ids = list(guide_names)
        all_feat_names = rna_names + guide_names
        all_feat_ids = rna_ids + guide_ids
        all_feat_types = (["Gene Expression"] * len(rna_names)) + (["CRISPR Guide Capture"] * len(guide_names))
        n_rna = len(rna_names)
        n_feat = len(all_feat_names)

        rng = np.random.default_rng(20261699)

        def _build_h5(path: Path, gsm_tag: str, non_int: bool = False, n_ntc: int = 55) -> None:
            n_targ = 24
            n_bad = 6
            n_col = n_ntc + n_targ + n_bad
            mat = np.zeros((n_feat, n_col), dtype=np.float64 if non_int else np.int32)
            for c in range(n_ntc + n_targ):
                g_choices = rng.choice(np.arange(1, n_rna), size=215, replace=False)
                mat[g_choices, c] = 2
                if c == 0:
                    mat[0, c] = 300  # >40% MT UMI; max_pct_mt=None MUST keep!
                if c < n_ntc:
                    gi = n_rna + (c % 2)
                    mat[gi, c] = 18
                    mat[n_rna + (1 - (c % 2)), c] = 1
                else:
                    ti = (c - n_ntc) % 3
                    mat[n_rna + 2 + ti, c] = 20

            c_low_umi = n_ntc + n_targ
            mat[1:210, c_low_umi] = 1
            mat[n_rna, c_low_umi] = 20
            c_low_ng = c_low_umi + 1
            mat[1:190, c_low_ng] = 3
            mat[n_rna, c_low_ng] = 20
            c_multi = c_low_umi + 2
            mat[1:216, c_multi] = 2
            mat[n_rna + 2, c_multi] = 15
            mat[n_rna + 3, c_multi] = 15
            c_safe = c_low_umi + 3
            mat[1:216, c_safe] = 2
            mat[n_rna + 5, c_safe] = 25
            c_inter = c_low_umi + 4
            mat[1:216, c_inter] = 2
            mat[n_rna + 6, c_inter] = 25
            c_amb = c_low_umi + 5
            mat[1:216, c_amb] = 2
            mat[n_rna, c_amb] = 2
            if non_int:
                mat[5, 0] = 1.5

            csc = sp.csc_matrix(mat)
            with h5py.File(path, "w") as f:
                m = f.create_group("matrix")
                m.create_dataset("data", data=csc.data)
                m.create_dataset("indices", data=csc.indices.astype(np.int32))
                m.create_dataset("indptr", data=csc.indptr.astype(np.int32))
                m.create_dataset("shape", data=np.array([n_feat, n_col], dtype=np.int32))
                m.create_dataset("barcodes", data=np.array([f"BC_{gsm_tag}_{k:03d}-1".encode() for k in range(n_col)]))
                fg = m.create_group("features")
                fg.create_dataset("id", data=np.array([x.encode() for x in all_feat_ids]))
                fg.create_dataset("name", data=np.array([x.encode() for x in all_feat_names]))
                fg.create_dataset("feature_type", data=np.array([x.encode() for x in all_feat_types]))

        with tempfile.TemporaryDirectory() as td:
            raw_g3 = Path(td) / "raw_g3"
            raw_g3.mkdir()
            for gsm, ntc_n in (("GSM8943713", 55), ("GSM8943714", 55), ("GSM8943715", 55), ("GSM8943716", 10)):
                _build_h5(raw_g3 / f"{gsm}_H1_filtered_feature_bc_matrix.h5", gsm, non_int=False, n_ntc=ntc_n)

            bad_dir = Path(td) / "bad_float"
            bad_dir.mkdir()
            _build_h5(bad_dir / "GSM8943713_bad.h5", "GSM8943713", non_int=True, n_ntc=20)
            with self.assertRaises(ValueError):
                build_gse295214_parts(bad_dir, hgnc=hgnc_l, allow_partial_channels=True)

            parts, p_audit = build_gse295214_parts(raw_g3, hgnc=hgnc_l, allow_partial_channels=True)
            self.assertEqual(len(parts), 4)
            self.assertTrue(p_audit["integer_rna_verified"])
            self.assertGreater(p_audit["n_low_ntc_fallback_cells"], 0)
            self.assertEqual([int(p.half_flags[0]) for p in parts], [1, 2, 1, 2])
            self.assertEqual(p_audit["n_ntc_half_a"], 110)
            self.assertEqual(p_audit["n_ntc_half_b"], 65)
            self.assertTrue(bool(parts[0].nt_mask[0]))  # >40% MT cell retained!

            pri_p, res_p = extract_gse295214_ntc_raw_counts(
                raw_g3, n_ctrl=64, cell_set_len=16, min_ctrl_floor=32, seed=DEFAULT_R16_SEED,
                return_resample=True, resample_seed=RESAMPLE_R16_SEED, allow_partial_channels=True,
            )
            self.assertEqual(pri_p.raw_counts.shape[0], 64)
            self.assertIsNotNone(res_p)
            self.assertEqual(res_p.raw_counts.shape[0], 64)
            self.assertTrue(set(pri_p.cell_barcodes.tolist()).isdisjoint(set(res_p.cell_barcodes.tolist())))
            self.assertTrue(all(":BC_" in bc and bc.startswith("GSM") for bc in pri_p.cell_barcodes))
            self.assertAlmostEqual(pri_p.median_ntc_umi, float(np.median(pri_p.cell_total_umis)), places=6)

            std_root = Path(td) / "std_lines"
            std_g3_dir = process_r10_line("gse295214_h1", parts, hgnc=hgnc_l, chunk_size=128, out_root=std_root, n_workers=1)
            qc_md = Path(td) / "R16_GSE295214_QC_COUNTS.md"
            gate_res = evaluate_gse295214_counts_gate_r16(
                std_line_dir=std_g3_dir,
                model_genes=r9_genes[:40],
                model_pert_names=["AAAS", "TP53", "MYC"],
                ntc_half_a=p_audit["n_ntc_half_a"],
                ntc_half_b=p_audit["n_ntc_half_b"],
                hgnc_bridge=hb,
                out_md_path=qc_md,
            )
            self.assertTrue(qc_md.exists())
            self.assertEqual(gate_res["status"], "SENSITIVITY_ONLY")
            self.assertTrue(gate_res["pass_per_half_ntc_ge_15"])
            self.assertFalse(gate_res["pass_g_eval_r9_ge_3000"])
            self.assertEqual(gate_res["fct_primary_pool"], list(GROUP3_PRIMARY_FCT_SOURCES))

            mock_ref = {
                "eval_genes": r9_genes[:80],
                "eval_perts": ["AAAS", "AAMP"],
                "target_gene_idx": np.array([r9_genes[:80].index("AAAS") if "AAAS" in r9_genes[:80] else -1, -1], dtype=np.int64),
                "ctrl_mean_libnorm": np.ones(80, dtype=np.float64),
            }
            sens_locked = build_fct_sensitivity_on_locked_targets_r16(mock_ref, hgnc_bridge=hb)
            self.assertEqual(sens_locked["eval_perts"], ["AAAS", "AAMP"])
            self.assertEqual(sens_locked["eval_genes"], r9_genes[:80])
            self.assertIn("vcc2025_h1", sens_locked["source_lines_fct"])

            grp_g3 = adjudicate_group_r16(
                per_line_h_i1={"nadig_hepg2": {"line_summary": "LINE_EQUIVALENT_NONZERO_SIGNAL"}, "gse295214_h1": {"line_summary": "LINE_MATCHED_SPECIFIC"}},
                per_line_h_i2={"nadig_hepg2": {"classification": "FCT_SUPERIOR_TO_MODEL"}, "gse295214_h1": {"classification": "INDISTINGUISHABLE"}},
            )
            self.assertFalse(grp_g3["group3_prospective"]["merged_with_group1_or_group2"])
            self.assertEqual(grp_g3["group3_prospective"]["K"], 1)
            self.assertEqual(len(grp_g3["group3_prospective"]["mandatory_disclosures"]), 2)

            # Confirm degenerate Group 2 uses PREREG_I v5 literal SINGLE_STRATUM_DESCRIPTIVE_ONLY
            grp_deg = adjudicate_group_r16(
                per_line_h_i1={"a549_ko": {"line_summary": "MISSING_STRATUM"}},
                per_line_h_i2={"a549_ko": {"classification": "INDISTINGUISHABLE"}},
                group1_lines=(), group2_crispri_lines=(), group2_crisprko_lines=("a549_ko",),
            )
            self.assertEqual(grp_deg["group2_crisprko"]["adjudication_status"], "SINGLE_STRATUM_DESCRIPTIVE_ONLY")

    def test_11_addendum_a_per_cell_batch_nadig_hepg2_parity_g4_and_perm100(self) -> None:
        """
        Verify:
          - `batch_onehot_map.torch` per-cell native `gem_group` vectors run cleanly in `CleanRoomArcStateAdapter`
            and differ from `index_zero`
          - `run_g4_factorial_for_line` executes Group 1 positional `B_M/B_S` and Group 2 `B_M=0, B_S=17` cleanly
          - `verify_nadig_hepg2_adata_pred_parity` enforces `median >= 0.90` & `pooled >= 0.95`
          - `summarize_perm100_r16` & `compute_practical_metrics_r16` match expectations
        """
        hb = get_hgnc_bridge()
        adapter = CleanRoomArcStateAdapter(
            model_key="STATE_HVG", checkpoint_id="zeroshot/hepg2", device="cpu",
            neutral_batch_mode="per_cell_native", delta_mode="minus_pred_ntc",
            allow_synthetic_smoke=True, hgnc_bridge=hb,
        )
        adapter.load_weights()
        b = adapter._bundles["hepg2"]
        bmap = b["batch_onehot_idx_map"]
        self.assertGreater(len(bmap), 10)
        b_labels = ["1", "10", "17", "35"] * 16
        mapped_idx = map_batch_labels_to_indices(b_labels, bmap)
        self.assertEqual(len(mapped_idx), 64)

        rng = np.random.default_rng(555)
        raw_m = rng.poisson(3.0, size=(64, len(b["gene_names"]))).astype(np.float32)
        res_native = adapter.predict_arm_deltas(
            "nadig_hepg2", "nadig_hepg2", raw_m, b["gene_names"], ["AAAS", "AAMP", "AAR2"],
            batch_labels=b_labels, n_ctrl=64, cell_set_len=32, is_synthetic_data=True, pre_sampled=True,
            per_cell_batch_labels=b_labels, hgnc_bridge=hb,
        )
        res_explicit = CleanRoomArcStateAdapter(
            model_key="STATE_HVG", checkpoint_id="zeroshot/hepg2", device="cpu",
            neutral_batch_mode="index_zero", delta_mode="minus_pred_ntc",
            allow_synthetic_smoke=True, hgnc_bridge=hb,
        ).predict_arm_deltas(
            "nadig_hepg2", "nadig_hepg2", raw_m, b["gene_names"], ["AAAS", "AAMP", "AAR2"],
            batch_labels=b_labels, n_ctrl=64, cell_set_len=32, is_synthetic_data=True, pre_sampled=True,
            per_cell_batch_indices=mapped_idx.tolist(), hgnc_bridge=hb,
        )
        np.testing.assert_allclose(res_native["delta_pred"], res_explicit["delta_pred"], atol=1e-6)

        res_idx0 = CleanRoomArcStateAdapter(
            model_key="STATE_HVG", checkpoint_id="zeroshot/hepg2", device="cpu",
            neutral_batch_mode="index_zero", delta_mode="minus_pred_ntc",
            allow_synthetic_smoke=True, hgnc_bridge=hb,
        ).predict_arm_deltas(
            "nadig_hepg2", "nadig_hepg2", raw_m, b["gene_names"], ["AAAS", "AAMP", "AAR2"],
            batch_labels=b_labels, n_ctrl=64, cell_set_len=32, is_synthetic_data=True, pre_sampled=True,
            hgnc_bridge=hb,
        )
        self.assertGreater(float(np.max(np.abs(res_native["delta_pred"] - res_idx0["delta_pred"]))), 1e-4)

        # Synthetic Group 1 G4 2x2 factorial (`nadig_hepg2` <-> `nadig_jurkat` positional B_M/B_S)
        class _DummyPool:
            def __init__(self, blabels, gnames):
                self.batch_labels = np.asarray(blabels, dtype=str)
                self.gene_names = gnames
                self.sampling_audit = {"n_sampled": 64, "n_sets": 2, "sampling_mode": "synthetic"}

        syn_eval_g = [hb.canonicalize(g) for g in b["gene_names"][:60]]
        syn_perts = [k for k in sorted(b["pert_lookup"].keys()) if k != "non-targeting"][:24]
        syn_sig = rng.normal(0.0, 0.4, size=(len(syn_perts), len(syn_eval_g)))
        syn_ref = {
            "eval_genes": syn_eval_g,
            "eval_perts": syn_perts,
            "target_gene_idx": np.full(len(syn_perts), -1, dtype=np.int64),
            "delta_obs": syn_sig,
            "delta_half_a": syn_sig + rng.normal(0.0, 0.04, size=syn_sig.shape),
            "delta_half_b": syn_sig + rng.normal(0.0, 0.04, size=syn_sig.shape),
            "half_valid": np.ones(len(syn_perts), dtype=bool),
            "delta_transfer": rng.normal(0.0, 0.05, size=syn_sig.shape),
            "source_grand_mean": np.zeros(len(syn_eval_g), dtype=np.float64),
        }
        g4_pools = {
            "nadig_hepg2": _DummyPool(["1", "10"] * 32, b["gene_names"]),
            "nadig_jurkat": _DummyPool(["17", "35"] * 32, b["gene_names"]),
        }
        g4_embs = {
            "nadig_hepg2": (rng.normal(0.0, 0.1, size=(64, 2058)).astype(np.float32), g4_pools["nadig_hepg2"].sampling_audit, {"source": "syn"}),
            "nadig_jurkat": (rng.normal(0.1, 0.1, size=(64, 2058)).astype(np.float32), g4_pools["nadig_jurkat"].sampling_audit, {"source": "syn"}),
        }
        with tempfile.TemporaryDirectory() as td_g4:
            g4_res = run_g4_factorial_for_line(
                target_line="nadig_hepg2",
                ref_bundle=syn_ref,
                raw_pools=g4_pools,
                se_embeddings=g4_embs,
                models_root=REPO_ROOT / "round16" / "shared" / "models",
                status_dir=Path(td_g4),
                device="cpu",
                n_boot=200,
                seed=DEFAULT_R16_SEED,
                hgnc=hb,
                is_synthetic_data=True,
                allow_synthetic_smoke=True,
            )
        self.assertIn("nadig_jurkat", g4_res["decompositions"])
        g4_g1_batch = g4_res["decompositions"]["nadig_jurkat"]["delta_batch_orth"]
        self.assertGreater(float(g4_g1_batch["rho_ab_orth"]), 0.5)
        self.assertGreaterEqual(float(g4_g1_batch["valid_frac"]), 0.95)
        self.assertTrue(np.isfinite(float(g4_g1_batch["estimate"])))
        self.assertTrue(np.isfinite(float(g4_res["decompositions"]["nadig_jurkat"]["batch_share_orth"])))

        # verify_nadig_hepg2_adata_pred_parity check (Addendum A)
        ref_mat = rng.normal(0.0, 0.5, size=(20, 200))
        ok_p = verify_nadig_hepg2_adata_pred_parity(ref_mat + rng.normal(0.0, 0.02, size=ref_mat.shape), ref_mat, np.arange(20))
        self.assertTrue(ok_p["passed"])
        self.assertGreaterEqual(ok_p["median_per_pert_pearson"], 0.90)
        self.assertGreaterEqual(ok_p["pooled_pearson"], 0.95)
        bad_p = verify_nadig_hepg2_adata_pred_parity(rng.normal(0.0, 0.5, size=ref_mat.shape), ref_mat, np.arange(20))
        self.assertFalse(bad_p["passed"])

        # Perm100 summary & compute_gene_label_perm_null_r16 check (Addendum A §4.2)
        p_sum = summarize_perm100_r16({"estimate": 0.40, "ci975_lo": 0.25}, [0.02 * np.sin(k) for k in range(100)])
        self.assertEqual(p_sum["n_perm"], 100)
        self.assertTrue(p_sum["matched_gt_q975"])
        self.assertAlmostEqual(p_sum["p_emp"], 1.0 / 101.0, places=9)
        self.assertEqual(len(p_sum["perm_scores"]), 100)
        fx_p = make_synthetic_pert_cohort(n_perts=40, n_genes=60, regime="matched_better_strong", seed=20261655)
        perm_full = compute_gene_label_perm_null_r16(
            fx_p["delta_obs"], fx_p["pred_matched"], fx_p["t_mean"], fx_p["target_cols"],
            source_grand_mean=fx_p["grand_mean"], delta_half_a=fx_p["half_a_obs"],
            delta_half_b=fx_p["half_b_obs"], half_valid=fx_p["half_valid"], n_perm=10, seed=DEFAULT_R16_SEED,
        )
        self.assertEqual(perm_full["n_perm"], 10)
        self.assertIn("dorth_summary", perm_full)
        self.assertIn("mean", perm_full)

        # Addendum A §2.2 raw-UMI collision tie-break check in ST-HVG column mapping
        aligned_tie, _ = align_hvg_expression_to_model(
            np.array([[1.0, 9.0, 3.0, 7.0]], dtype=np.float32),
            source_genes=["AARS", "AARS1", "DARS", "DARS1"],
            model_genes=["AARS1", "DARS1"],
            hgnc_bridge=hb,
            col_mean_scores=np.array([10.0, 50.0, 25.0, 25.0], dtype=np.float64),
        )
        self.assertAlmostEqual(float(aligned_tie[0, 0]), 9.0, places=5)
        self.assertAlmostEqual(float(aligned_tie[0, 1]), 3.0, places=5)

        # Group 2 G4 mandatory (17, 0) symmetry check + spread {0, 1, 17, 35}
        g4_g2_pools = {
            "a549_ko": _DummyPool(["0"] * 64, b["gene_names"]),
            "nadig_hepg2": _DummyPool(["0"] * 64, b["gene_names"]),
            "nadig_jurkat": _DummyPool(["0"] * 64, b["gene_names"]),
            "dld1": _DummyPool(["0"] * 64, b["gene_names"]),
            "cd4t_rest": _DummyPool(["0"] * 64, b["gene_names"]),
            "feng_ipsc": _DummyPool(["0"] * 64, b["gene_names"]),
        }
        g4_g2_embs = {
            k: (rng.normal(0.05 if k == "a549_ko" else 0.0, 0.1, size=(64, 2058)).astype(np.float32), g4_g2_pools[k].sampling_audit, {"source": "syn"})
            for k in g4_g2_pools
        }
        with tempfile.TemporaryDirectory() as td_g2:
            g4_g2 = run_g4_factorial_for_line(
                target_line="a549_ko", ref_bundle=syn_ref, raw_pools=g4_g2_pools, se_embeddings=g4_g2_embs,
                models_root=REPO_ROOT / "round16" / "shared" / "models", status_dir=Path(td_g2),
                device="cpu", n_boot=200, seed=DEFAULT_R16_SEED, hgnc=hb, is_synthetic_data=True, allow_synthetic_smoke=True,
            )
        self.assertIsNotNone(g4_g2["symmetry_check_b17_vs_b0"])
        self.assertEqual(sorted(g4_g2["batch_spread"].keys()), ["0", "1", "17", "35"])
        est_0_17 = float(g4_g2["decompositions"]["dld1"]["delta_batch_orth"]["estimate"])
        est_17_0 = float(g4_g2["symmetry_check_b17_vs_b0"]["dld1"]["delta_batch_orth"]["estimate"])
        self.assertGreater(float(g4_g2["decompositions"]["dld1"]["delta_batch_orth"]["rho_ab_orth"]), 0.5)
        self.assertGreaterEqual(float(g4_g2["decompositions"]["dld1"]["delta_batch_orth"]["valid_frac"]), 0.95)
        self.assertTrue(np.isfinite(est_0_17) and np.isfinite(est_17_0))
        self.assertAlmostEqual(est_0_17, -est_17_0, places=6)

    def test_12_addendum_a_cross_party_tolerances_and_barcodes(self) -> None:
        """
        Verify Addendum A tolerance checker:
          - identical barcodes, pred diff <= 1e-4, same-array recompute <= 1e-10,
            cross-party metrics <= 1e-5, 100% verdict parity
          - rejects any barcode mismatch or pred diff > 1e-4 cleanly.
        """
        fx = make_synthetic_pert_cohort(n_perts=80, n_genes=100, regime="matched_better_strong", seed=20261620)
        ref = {
            "eval_genes": [f"G{i}" for i in range(100)],
            "eval_perts": [f"P{i}" for i in range(80)],
            "target_gene_idx": fx["target_cols"],
            "delta_obs": fx["delta_obs"],
            "delta_half_a": fx["half_a_obs"],
            "delta_half_b": fx["half_b_obs"],
            "half_valid": fx["half_valid"],
            "delta_transfer": fx["t_mean"],
            "source_grand_mean": fx["grand_mean"],
            "delta_fct": fx["pred_fct"],
            "ctrl_mean_libnorm": np.ones(100),
        }
        m_bcs = [f"MB_{i:03d}" for i in range(64)]
        mis_bcs = {s: [f"{s}_{i:03d}" for i in range(64)] for s in fx["mis_preds"]}
        with tempfile.TemporaryDirectory() as td:
            npz_p = export_line_arm_npz(
                Path(td) / "p.npz", "STATE_SE", "zeroshot/hepg2", "a549_ko", ref,
                fx["pred_matched"], fx["mis_preds"], fx["seen_sources"], fx["unseen_sources"],
                "index_zero", "minus_pred_ntc", "se600m_raw_log1p", {"t": 1},
                matched_barcodes=m_bcs, mismatched_barcodes=mis_bcs,
                delta_pred_resample_eval=fx["pred_matched"], resample_barcodes=m_bcs,
            )
            tiny_eps = 1e-8
            main_preds = {
                "matched": fx["pred_matched"] + tiny_eps,
                "resample": fx["pred_matched"] + tiny_eps,
                **{s: v + tiny_eps for s, v in fx["mis_preds"].items()},
            }
            pass_chk = verify_cross_party_arrays_and_barcodes(
                npz_p, main_preds, main_barcodes={"matched": m_bcs, "resample": m_bcs, **mis_bcs}, n_boot=150, seed=DEFAULT_R16_SEED
            )
            self.assertTrue(pass_chk["all_passed"])
            self.assertLessEqual(pass_chk["same_array_metric_max_diff"], 1e-10)
            self.assertLessEqual(pass_chk["cross_party_metric_max_diff"], 1e-5)

            bad_bc_chk = verify_cross_party_arrays_and_barcodes(
                npz_p, main_preds, main_barcodes={"matched": ["WRONG"] + m_bcs[1:], "resample": m_bcs, **mis_bcs}, n_boot=150
            )
            self.assertFalse(bad_bc_chk["all_passed"])
            self.assertFalse(bad_bc_chk["barcodes_exact_match"])

            se_npz = Path(td) / "se.npz"
            fake_emb = np.random.default_rng(1).normal(0.0, 0.05, size=(64, 2058)).astype(np.float32)
            np.savez(se_npz, a549_ko__X_state=fake_emb)
            se_ok = verify_cross_party_arrays_and_barcodes(
                npz_p, main_preds, main_barcodes={"matched": m_bcs, "resample": m_bcs, **mis_bcs},
                n_boot=150, seed=DEFAULT_R16_SEED,
                audit_se_npz_path=se_npz, main_se_embeddings={"a549_ko": fake_emb + 5e-6},
            )
            self.assertTrue(se_ok["all_passed"])
            self.assertTrue(se_ok["se_emb_within_tol"])
            self.assertLessEqual(se_ok["se_emb_max_abs_diff"], 1e-5)

            se_bad = verify_cross_party_arrays_and_barcodes(
                npz_p, main_preds, main_barcodes={"matched": m_bcs, "resample": m_bcs, **mis_bcs},
                n_boot=150, seed=DEFAULT_R16_SEED,
                audit_se_npz_path=se_npz, main_se_embeddings={"a549_ko": fake_emb + 2e-5},
            )
            self.assertFalse(se_bad["all_passed"])
            self.assertFalse(se_bad["se_emb_within_tol"])

            bad_pred_chk = verify_cross_party_arrays_and_barcodes(
                npz_p, {**main_preds, "matched": fx["pred_matched"] + 2e-4}, n_boot=150
            )
            self.assertFalse(bad_pred_chk["all_passed"])
            self.assertFalse(bad_pred_chk["pred_within_tol"])

    def test_13_cross_process_pythonhashseed_reproducibility(self) -> None:
        """
        Launch two independent fresh Python subprocesses with distinct `PYTHONHASHSEED` values (`17` vs `999983`),
        running full end-to-end inference + `.npz` export strictly on identical synthetic Poisson raw counts:
          - CPU (`ST-HVG-Replogle` + `SE-600M -> ST-SE-Replogle`): exported predicted delta arrays MUST be
            bitwise identical (`max_abs_diff == 0.0`) and sampled barcodes identical.
          - MPS (`float32`, when available): exported predicted delta arrays MUST match within `max_abs_diff <= 1e-6`.
        Tests write strictly inside independent temporary directories.
        """
        import textwrap
        import torch

        worker_code = textwrap.dedent(
            f"""
            import json
            import sys
            from pathlib import Path
            import numpy as np
            import torch

            audit_dir = {repr(str(AUDIT_DIR))}
            sys.path.insert(0, audit_dir)
            from audit_r16_sampling import DEFAULT_R16_SEED, RESAMPLE_R16_SEED, sample_nonoverlapping_resample_sets_r16, sample_stratified_control_sets_r16
            from audit_r16_state_adapter import CleanRoomArcStateAdapter, CleanRoomSE600MEncoder, get_hgnc_bridge
            from run_audit_r16 import export_line_arm_npz, project_model_deltas_to_eval_genes

            out_dir = Path(sys.argv[1])
            hb = get_hgnc_bridge()
            rng = np.random.default_rng(DEFAULT_R16_SEED)

            # Build fixed synthetic raw integer count matrices (`160` cells for HVG, `64` for SE-600M)
            ad_ref = CleanRoomArcStateAdapter("STATE_HVG", "zeroshot/hepg2", device="cpu", allow_synthetic_smoke=True, hgnc_bridge=hb)
            ad_ref.load_weights()
            g_model = ad_ref._bundles["hepg2"]["gene_names"]
            eval_perts = ["AAAS", "AAMP", "AAR2", "AATF"]
            eval_genes = [hb.canonicalize(g) for g in g_model[:120]]
            se_tmp = CleanRoomSE600MEncoder(device="cpu", batch_size=16, hgnc_bridge=hb)
            se_tmp.load()
            pe_sub = se_tmp._pe_keys[:2400]

            hvg_m_raw = rng.poisson(2.5, size=(256, len(g_model))).astype(np.float32)
            hvg_m_b = [f"B{{i % 3}}" for i in range(256)]
            hvg_m_bc = [f"HVG_M_{{i:04d}}" for i in range(256)]
            hvg_s_raw = rng.poisson(6.0, size=(160, len(g_model))).astype(np.float32)
            hvg_s_b = [f"L{{i % 2}}" for i in range(160)]
            hvg_s_bc = [f"HVG_S_{{i:04d}}" for i in range(160)]

            se_m_raw = rng.poisson(2.5, size=(160, len(pe_sub))).astype(np.float32)
            se_m_b = [f"B{{i % 2}}" for i in range(160)]
            se_m_bc = [f"SE_M_{{i:04d}}" for i in range(160)]
            se_s_raw = rng.poisson(4.5, size=(160, len(pe_sub))).astype(np.float32)
            se_s_b = [f"U{{i % 2}}" for i in range(160)]
            se_s_bc = [f"SE_S_{{i:04d}}" for i in range(160)]

            dummy_sig = rng.normal(0.0, 0.3, size=(len(eval_perts), len(eval_genes)))
            ref_bundle = {{
                "eval_genes": eval_genes,
                "eval_perts": eval_perts,
                "target_gene_idx": np.full(len(eval_perts), -1, dtype=np.int64),
                "delta_obs": dummy_sig,
                "delta_half_a": dummy_sig,
                "delta_half_b": dummy_sig,
                "half_valid": np.ones(len(eval_perts), dtype=bool),
                "delta_transfer": np.zeros_like(dummy_sig),
                "source_grand_mean": np.zeros(len(eval_genes), dtype=np.float64),
                "delta_fct": np.zeros_like(dummy_sig),
                "ctrl_mean_libnorm": np.ones(len(eval_genes), dtype=np.float64),
            }}

            devices = ["cpu"] + (["mps"] if torch.backends.mps.is_available() else [])
            for dev in devices:
                # 1. STATE_HVG end-to-end
                ad_hvg = CleanRoomArcStateAdapter(
                    model_key="STATE_HVG", checkpoint_id="zeroshot/hepg2", device=dev,
                    neutral_batch_mode="index_zero", hvg_preproc_mode="official_median_log1p",
                    allow_synthetic_smoke=True, hgnc_bridge=hb,
                )
                arm_hvg = ad_hvg.run_matched_and_mismatched_arms(
                    target_line="nadig_hepg2",
                    matched_ctrl_matrix=hvg_m_raw,
                    matched_genes=g_model,
                    mismatched_ctrls={{"nadig_jurkat": (hvg_s_raw, g_model, hvg_s_b)}},
                    eval_perts=eval_perts,
                    matched_batch_labels=hvg_m_b,
                    seed=DEFAULT_R16_SEED,
                    n_ctrl=64,
                    cell_set_len=32,
                    is_synthetic_data=True,
                    hgnc_bridge=hb,
                    matched_barcodes=hvg_m_bc,
                    mismatched_barcodes={{"nadig_jurkat": hvg_s_bc}},
                )
                res_samp = sample_nonoverlapping_resample_sets_r16(
                    arm_hvg["matched"]["sampling_audit"], n_cells=len(hvg_m_bc), batch_labels=hvg_m_b,
                    n_ctrl=64, cell_set_len=32, min_ctrl_floor=32, resample_seed=RESAMPLE_R16_SEED, cell_barcodes=hvg_m_bc,
                )
                res_idx = res_samp["selected_indices"]
                res_hvg = ad_hvg.predict_arm_deltas(
                    "nadig_hepg2", "nadig_hepg2__resample", hvg_m_raw[res_idx], g_model, eval_perts,
                    n_ctrl=len(res_idx), cell_set_len=32, is_synthetic_data=True, hgnc_bridge=hb,
                    pre_sampled=True, precomputed_sampling_audit=res_samp, cell_barcodes=res_samp["selected_barcodes"],
                )
                dm_hvg = project_model_deltas_to_eval_genes(arm_hvg["matched"]["delta_pred"], g_model, eval_genes, hgnc_bridge=hb)
                ds_hvg = project_model_deltas_to_eval_genes(arm_hvg["mismatched"]["nadig_jurkat"]["delta_pred"], g_model, eval_genes, hgnc_bridge=hb)
                dr_hvg = project_model_deltas_to_eval_genes(res_hvg["delta_pred"], g_model, eval_genes, hgnc_bridge=hb)
                export_line_arm_npz(
                    out_dir / f"{{dev}}_hvg.npz", "STATE_HVG", "zeroshot/hepg2", "nadig_hepg2", ref_bundle,
                    dm_hvg, {{"nadig_jurkat": ds_hvg}}, ("nadig_jurkat",), (), "index_zero", "minus_pred_ntc", "official_median_log1p",
                    {{"dev": dev}}, matched_barcodes=arm_hvg["matched"]["selected_barcodes"],
                    mismatched_barcodes={{"nadig_jurkat": arm_hvg["mismatched"]["nadig_jurkat"]["selected_barcodes"]}},
                    delta_pred_resample_eval=dr_hvg, resample_barcodes=res_samp["selected_barcodes"],
                )

                # 2. SE-600M -> STATE_SE end-to-end
                se_enc = CleanRoomSE600MEncoder(device=dev, batch_size=16, hgnc_bridge=hb)
                ad_se = CleanRoomArcStateAdapter(
                    model_key="STATE_SE", checkpoint_id="zeroshot/hepg2", device=dev,
                    neutral_batch_mode="index_zero", allow_synthetic_smoke=True,
                    se600m_encoder=se_enc, hgnc_bridge=hb, se600m_use_hgnc_rescue=False,
                )
                if dev == "cpu":
                    # Encode 2 matched + 2 mismatched cells through full SE-600M on CPU (~13s)
                    # then tile to 1x32 and run full ST-SE inference on CPU (bitwise exact across processes)
                    samp_m = sample_stratified_control_sets_r16(
                        n_cells=len(se_m_bc), batch_labels=se_m_b, n_ctrl=32, cell_set_len=32,
                        seed=DEFAULT_R16_SEED, cell_barcodes=se_m_bc,
                    )
                    samp_s = sample_stratified_control_sets_r16(
                        n_cells=len(se_s_bc), batch_labels=se_s_b, n_ctrl=32, cell_set_len=32,
                        seed=DEFAULT_R16_SEED, cell_barcodes=se_s_bc,
                    )
                    emb_m2, st_m = se_enc.encode_cells(
                        se_m_raw[samp_m["selected_indices"][:2]], pe_sub, seed=DEFAULT_R16_SEED, use_hgnc_rescue=False, dataset_name="nadig_hepg2"
                    )
                    emb_s2, st_s = se_enc.encode_cells(
                        se_s_raw[samp_s["selected_indices"][:2]], pe_sub, seed=DEFAULT_R16_SEED, use_hgnc_rescue=False, dataset_name="cd4t_rest"
                    )
                    emb_r2, st_r = se_enc.encode_cells(
                        se_m_raw[samp_m["selected_indices"][2:4]], pe_sub, seed=RESAMPLE_R16_SEED, use_hgnc_rescue=False, dataset_name="nadig_hepg2__resample"
                    )
                    emb_m32 = np.tile(emb_m2, (16, 1))
                    emb_s32 = np.tile(emb_s2, (16, 1))
                    arm_se = ad_se.run_matched_and_mismatched_arms(
                        target_line="nadig_hepg2",
                        matched_ctrl_matrix=None,
                        matched_genes=pe_sub,
                        mismatched_ctrls={{"cd4t_rest": (None, pe_sub, ["U0"] * 32)}},
                        eval_perts=eval_perts,
                        matched_batch_labels=["B0"] * 32,
                        seed=DEFAULT_R16_SEED,
                        n_ctrl=32,
                        cell_set_len=32,
                        is_synthetic_data=True,
                        precomputed_embeddings={{
                            "nadig_hepg2": (emb_m32, samp_m, st_m),
                            "cd4t_rest": (emb_s32, samp_s, st_s),
                        }},
                        hgnc_bridge=hb,
                        pre_sampled=True,
                        matched_barcodes=samp_m["selected_barcodes"],
                        mismatched_barcodes={{"cd4t_rest": samp_s["selected_barcodes"]}},
                    )
                    np.savez(
                        out_dir / "cpu_se_emb.npz",
                        emb_m=emb_m2, emb_s=emb_s2, emb_r=emb_r2,
                        sha_m=np.asarray(st_m["emb_sha256"]),
                        sha_s=np.asarray(st_s["emb_sha256"]),
                        sha_r=np.asarray(st_r["emb_sha256"]),
                    )
                else:
                    arm_se = ad_se.run_matched_and_mismatched_arms(
                        target_line="nadig_hepg2",
                        matched_ctrl_matrix=se_m_raw,
                        matched_genes=pe_sub,
                        mismatched_ctrls={{"cd4t_rest": (se_s_raw, pe_sub, se_s_b)}},
                        eval_perts=eval_perts,
                        matched_batch_labels=se_m_b,
                        seed=DEFAULT_R16_SEED,
                        n_ctrl=32,
                        cell_set_len=32,
                        is_synthetic_data=True,
                        hgnc_bridge=hb,
                        matched_barcodes=se_m_bc,
                        mismatched_barcodes={{"cd4t_rest": se_s_bc}},
                    )
                    np.savez(
                        out_dir / "mps_se_emb.npz",
                        emb_m=arm_se["matched"]["ctrl_sets_embed"],
                        emb_s=arm_se["mismatched"]["cd4t_rest"]["ctrl_sets_embed"],
                        sha_m=np.asarray(arm_se["matched"]["gene_mapping"]["emb_sha256"]),
                        sha_s=np.asarray(arm_se["mismatched"]["cd4t_rest"]["gene_mapping"]["emb_sha256"]),
                    )
                dm_se = project_model_deltas_to_eval_genes(arm_se["matched"]["delta_pred"], g_model, eval_genes, hgnc_bridge=hb)
                ds_se = project_model_deltas_to_eval_genes(arm_se["mismatched"]["cd4t_rest"]["delta_pred"], g_model, eval_genes, hgnc_bridge=hb)
                export_line_arm_npz(
                    out_dir / f"{{dev}}_se.npz", "STATE_SE", "zeroshot/hepg2", "nadig_hepg2", ref_bundle,
                    dm_se, {{"cd4t_rest": ds_se}}, (), ("cd4t_rest",), "index_zero", "minus_pred_ntc", "se600m_raw_log1p",
                    {{"dev": dev}}, matched_barcodes=arm_se["matched"]["selected_barcodes"],
                    mismatched_barcodes={{"cd4t_rest": arm_se["mismatched"]["cd4t_rest"]["selected_barcodes"]}},
                )
            """
        )

        with tempfile.TemporaryDirectory() as td1, tempfile.TemporaryDirectory() as td2:
            for hash_seed, out_td in (("17", td1), ("999983", td2)):
                env = dict(os.environ)
                env["PYTHONHASHSEED"] = hash_seed
                subprocess.run([sys.executable, "-c", worker_code, out_td], env=env, check=True)

            # 1. CPU: exact SHA-256 & bitwise equality (`max_abs_diff == 0.0`) on X_state embeddings + downstream arrays
            with np.load(Path(td1) / "cpu_se_emb.npz", allow_pickle=False) as e1, np.load(
                Path(td2) / "cpu_se_emb.npz", allow_pickle=False
            ) as e2:
                self.assertEqual(str(e1["sha_m"]), str(e2["sha_m"]))
                self.assertEqual(str(e1["sha_s"]), str(e2["sha_s"]))
                self.assertEqual(str(e1["sha_r"]), str(e2["sha_r"]))
                self.assertNotEqual(str(e1["sha_m"]), str(e1["sha_r"]))
                np.testing.assert_array_equal(e1["emb_m"], e2["emb_m"])
                np.testing.assert_array_equal(e1["emb_s"], e2["emb_s"])
                np.testing.assert_array_equal(e1["emb_r"], e2["emb_r"])

            for tag, mis_k in (("hvg", "nadig_jurkat"), ("se", "cd4t_rest")):
                with np.load(Path(td1) / f"cpu_{tag}.npz", allow_pickle=False) as z1, np.load(
                    Path(td2) / f"cpu_{tag}.npz", allow_pickle=False
                ) as z2:
                    self.assertEqual(z1["matched_barcodes"].tolist(), z2["matched_barcodes"].tolist())
                    self.assertEqual(
                        z1[f"mismatched_barcodes__{mis_k}"].tolist(),
                        z2[f"mismatched_barcodes__{mis_k}"].tolist(),
                    )
                    np.testing.assert_array_equal(z1["delta_pred_matched"], z2["delta_pred_matched"])
                    np.testing.assert_array_equal(z1["delta_pred_arms"], z2["delta_pred_arms"])
                    np.testing.assert_array_equal(
                        z1[f"delta_pred_mismatched__{mis_k}"], z2[f"delta_pred_mismatched__{mis_k}"]
                    )
                    if "delta_pred_resample" in z1:
                        self.assertEqual(z1["resample_barcodes"].tolist(), z2["resample_barcodes"].tolist())
                        np.testing.assert_array_equal(z1["delta_pred_resample"], z2["delta_pred_resample"])

            # 2. MPS (when available): X_state <= 1e-5 (Addendum A §11) & predicted delta max_abs_diff <= 1e-6
            if torch.backends.mps.is_available():
                with np.load(Path(td1) / "mps_se_emb.npz", allow_pickle=False) as e1, np.load(
                    Path(td2) / "mps_se_emb.npz", allow_pickle=False
                ) as e2:
                    self.assertLessEqual(float(np.max(np.abs(e1["emb_m"] - e2["emb_m"]))), 1e-5)
                    self.assertLessEqual(float(np.max(np.abs(e1["emb_s"] - e2["emb_s"]))), 1e-5)
                for tag, mis_k in (("hvg", "nadig_jurkat"), ("se", "cd4t_rest")):
                    with np.load(Path(td1) / f"mps_{tag}.npz", allow_pickle=False) as z1, np.load(
                        Path(td2) / f"mps_{tag}.npz", allow_pickle=False
                    ) as z2:
                        self.assertEqual(z1["matched_barcodes"].tolist(), z2["matched_barcodes"].tolist())
                        max_d = float(np.max(np.abs(z1["delta_pred_arms"] - z2["delta_pred_arms"])))
                        self.assertLessEqual(max_d, 1e-6)


if __name__ == "__main__":
    unittest.main()
