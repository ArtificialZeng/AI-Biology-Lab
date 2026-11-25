#!/usr/bin/env python3
"""
Clean-Room Unit, Parity, & Real-Weight ST-HVG Synthetic Smoke Suite (`round16/audit/`, `A-T3`).

Strictly NEVER reads, imports, or references `round16/src/` or any real biological cell data.
Verifies:
  1. Exact numerical parity (`<= 1e-12`) against frozen `round10/src/metrics_r10.py`
  2. `r_dorth_cal` Property P1 null invariance (`<= 1e-14`) for affine predictors in `span{1, t_mean, b_S}`
  3. Batch-imbalanced `16 × 32 = 512` stratified control sampling (`seed = 20261601`), tie-breaking, truncation, exclusion floor
  4. All TOST 4-way + `within_delta` (`valid_frac` gated) + `rho_ab_orth` + two-stratum + degenerate Group-2 rules (`PREREG_I v4`)
  5. Canonical multi-regime synthetic fixtures (varying N, null-signal equiv, nonzero-signal equiv, matched_better strong & sub_delta, mismatched_better, stratum split)
  6. Leakage & freeze gates (`BlindingGateError`, `TrainingLeakageError`)
  7. End-to-end real-weight `ST-HVG-Replogle` synthetic count smoke test (matched vs. batch-imbalanced mismatched synthetic controls, proving control sensitivity under both neutral covariate modes `bypass_zero` & `index_zero`)
"""

from __future__ import annotations

import os
import sys
import tempfile
import unittest
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

AUDIT_DIR = Path(__file__).resolve().parent
REPO_ROOT = AUDIT_DIR.parents[1]
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
for _p in (str(AUDIT_DIR), str(R10_SRC_DIR)):
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
    compute_orth_suite_r16,
    compute_paired_contrast_r16,
    compute_winner_curse_weight_shift_r16,
)
from audit_r16_sampling import (  # noqa: E402
    build_mismatched_strata_r16,
    libnorm_log1p,
    sample_stratified_control_sets_r16,
)
from audit_r16_state_adapter import (  # noqa: E402
    BlindingGateError,
    CleanRoomArcStateAdapter,
    CleanRoomSE600MEncoder,
    TrainingLeakageError,
    assert_freeze_gates_r16,
    assert_zero_training_leakage,
    encode_with_loaded_se_model,
    get_hgnc_bridge,
    map_genes_to_se600m_pe,
    preprocess_hvg_expression,
)
from run_audit_r16 import (  # noqa: E402
    main as runner_main,
    run_synthetic_runner_smoke,
)
from metrics_r10 import (  # noqa: E402
    compute_joint_paired_r_orth_diff_r10,
    compute_residual_suite_r10,
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

        # 2. Verify SE-600M gene mapping with HGNCBridge canonical rescues
        pe_fake = ["A1BG", "AARS", "DDI2", "RSC1A1", "MYO18A", "TIAF1", "TAZ"]
        src_g = ["AARS1", "DDI2", "TIAF1", "MISSING_GENE"]
        src_c, pe_g, mstats = map_genes_to_se600m_pe(src_g, pe_fake, hgnc_bridge=hb)
        self.assertEqual(mstats["n_hgnc_canon_rescued"], 1)  # AARS1 -> AARS
        self.assertEqual(mstats["n_raw_exact"], 2)           # DDI2, TIAF1 exact over collisions
        self.assertEqual(pe_g.tolist(), [1, 2, 5])
        self.assertEqual(src_c.tolist(), [0, 1, 2])

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
            "embeddings": {"current": "e1", "e1": {"num": 16, "size": 16}},
            "tokenizer": {"token_dim": 16},
            "loss": {"name": "tabular"},
        })
        se_tiny = StateEmbeddingModel(
            token_dim=16, d_model=8, nhead=2, d_hid=16, nlayers=1, output_dim=8, dropout=0.0, cfg=cfg
        )
        pe_tensor = torch.randn(16, 16, dtype=torch.float32)
        se_tiny.pe_embedding = torch.nn.Embedding.from_pretrained(pe_tensor)
        syn_counts = np.random.default_rng(7).poisson(4.0, size=(8, 7)).astype(np.float32)
        emb_tiny, st_tiny = encode_with_loaded_se_model(
            se_tiny, pe_keys=pe_fake + [f"G{i}" for i in range(9)],
            ctrl_counts=syn_counts, source_genes=pe_fake, hgnc_bridge=hb,
            pad_length=16, cls_token_idx=3, batch_size=4, seed=DEFAULT_R16_SEED, device=torch.device("cpu"),
        )
        self.assertEqual(emb_tiny.shape, (8, 18))  # 8 emb + 10 ds_emb
        self.assertEqual(st_tiny["n_mapped_genes"], 7)
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
        samp = sample_stratified_control_sets_r16(n_cells=128, n_ctrl=128, cell_set_len=32, min_ctrl_floor=128, seed=DEFAULT_R16_SEED)
        res_full = adapter.predict_arm_deltas(
            "nadig_hepg2", "nadig_hepg2", raw_x, m_genes, ["AAAS", "AARS1"],
            n_ctrl=128, cell_set_len=32, is_synthetic_data=True, hgnc_bridge=hb,
        )
        res_presamp = adapter.predict_arm_deltas(
            "nadig_hepg2", "nadig_hepg2", raw_x[samp["selected_indices"]], m_genes, ["AAAS", "AARS1"],
            n_ctrl=128, cell_set_len=32, is_synthetic_data=True, pre_sampled=True,
            precomputed_sampling_audit=samp, hgnc_bridge=hb,
        )
        np.testing.assert_allclose(res_full["delta_pred"], res_presamp["delta_pred"], atol=1e-7)

        # Confirm padded=True matches 4 disjoint padded=False 32-cell forwards & ignores ST-HVG untrained gene_decoder
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

        # 5. End-to-end synthetic runner + .npz schema check
        with tempfile.TemporaryDirectory() as td:
            res_run = run_synthetic_runner_smoke(output_dir=Path(td), device="cpu", n_sets=4, n_perts=6, include_real_se600m=False)
            npz_f = Path(res_run["st_hvg_smoke"]["npz_path"])
            self.assertTrue(npz_f.exists())
            with np.load(npz_f, allow_pickle=False) as z:
                for req in (
                    "model_key", "checkpoint_id", "target_line", "neutral_batch_mode",
                    "eval_genes", "eval_perts", "target_gene_idx",
                    "delta_obs", "delta_half_a", "delta_half_b", "half_valid",
                    "delta_transfer", "delta_fct", "arm_names", "arm_strata",
                    "delta_pred_arms", "delta_pred_matched",
                    "delta_pred_mismatched__nadig_jurkat", "delta_pred_mismatched__cd4t_rest",
                ):
                    self.assertIn(req, z.files)
                self.assertEqual(z["delta_pred_matched"].shape, (6, 240))
                self.assertEqual(str(z["neutral_batch_mode"]), "index_zero")
                self.assertGreater(res_run["st_hvg_smoke"]["gene_mapping"]["n_hgnc_canon_rescued"], 0)

    def test_8_real_se600m_safetensors_mps_fp32_to_st_se_synthetic_smoke(self) -> None:
        """
        End-to-end Gate-2 verification using REAL `SE-600M` (`model.safetensors`, Apple `MPS fp32`
        or CPU fallback) + REAL `ST-SE-Replogle` (`zeroshot/hepg2`, decoded via `clip(gene_decoder(preds), 0, 14)`)
        strictly on 32 (`1 x 32`) synthetic Poisson raw integer counts (ZERO real cells used).
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
        raw_b = rng.poisson(11.0, size=(32, len(pe_sub))).astype(np.float32)

        emb_a, st_a = se_enc.encode_cells(raw_a, pe_sub, seed=DEFAULT_R16_SEED)
        self.assertEqual(emb_a.shape, (32, 2058))
        self.assertTrue(np.all(np.isfinite(emb_a)))
        np.testing.assert_allclose(np.linalg.norm(emb_a[:, :2048], axis=1), np.ones(32), atol=1e-4)
        self.assertEqual(st_a["n_mapped_genes"], 2400)

        # Reject non-float32 tensor dtype guard on MPS if MPS active
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
                "nadig_hepg2": (emb_a, {"n_sampled": 32, "n_sets": 1, "sampling_mode": "pre_sampled"}, st_a)
            },
            hgnc_bridge=hb,
            pre_sampled=True,
        )
        d_m = arm_res["matched"]["delta_pred"]
        d_s = arm_res["mismatched"]["cd4t_rest"]["delta_pred"]
        self.assertEqual(d_m.shape, (4, 2000))
        self.assertTrue(np.all(np.isfinite(d_m)))
        self.assertGreater(float(np.linalg.norm(d_m - d_s)), 1e-2)


if __name__ == "__main__":
    unittest.main()
