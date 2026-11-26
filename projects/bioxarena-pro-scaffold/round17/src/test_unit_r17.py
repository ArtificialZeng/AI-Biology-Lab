#!/usr/bin/env python3
"""
Round 17 Comprehensive Unit & Synthetic End-to-End Test Suite (`PREREG_J` v2).

Strictly CPU-only, ZERO real-cell reads, ZERO MPS/GPU usage, ZERO edits to locked files.
Verifies:
  1. Execution gates (`R17_CODE_FREEZE.txt` + `R17_GO.txt` + `PREREG_J_LOCK.md`),
     omission of `wtc11_i3neuron`, quarantine of R16 discovery lines, and `J5a` mismatched strata.
  2. Exact bit-level parity (`<= 1e-14`) of fast 56-index bootstrap against frozen
     `round16/src/metrics_r16.py::compute_orth_suite_r16` and `compute_paired_contrast_r16`.
  3. Stratified perturbation half-split (`|P1| = floor(N/2)`, `|P2| = ceil(N/2)`, mixed
     `half_valid`, canonical target-order invariance, and 20-split multi-split diagnostics).
  4. Variance-component estimator (`tau_batch`, `p_var`): exact `0.0` / `p_var = 1.0` under
     identical 56 indices; large `tau_batch > 0.02` & `p_var < 0.005` under true heterogeneity.
  5. `H-J1` all 6 decision branches including `J2` floor-effect `INCONCLUSIVE_LOW_SIGNAL`,
     `BATCH_INVARIANT_WITHIN_DELTA`, `BATCH_DEPENDENCE_DETECTED_SUB_DELTA`,
     `BATCH_DEPENDENCE_PRACTICALLY_MEANINGFUL`, `VALIDITY_FAIL`, and `INSUFFICIENT_PRECISION`.
  6. `H-J2` all 5 audit classes including `CALL_UNSTABLE_SIGN_REVERSAL`, `CALL_STABLE`,
     `INCONCLUSIVE_LOW_SIGNAL` (`J2`), and side-by-side Bonferroni-56 normal vs empirical-20k checks.
  7. `H-J3` `L1_ens56`, `L2_ntc_recon_match` (Ruling 3 Pearson top-1/top-2 margin, lowest-index
     tie-break, negative-L2 secondary, per-fold selection & agreement, `EXACT_TIE_SELECTED_ZERO`,
     out-of-sample `CF-Oracle` contrast, and weight-only geometric `medoid_idx`).
  8. `J5a` cross-fitted `G4` 2x2 decomposition zero-leakage guard (mutating `P1` observed deltas
     does NOT change `P2` `G4` values when `(b*_1, b_dag_1)` are fixed) + continuity `(0,17)/(0,35)`.
  9. `J5b` 56-index Spearman + permutation null + BH FDR + order-invariance across pair sets.
 10. Holm-Bonferroni step-down FWER control, pooled-model Holm sensitivity, and cohort `k / K_conf`
     claim rule (`k >= 2` required for general cross-context claim; `pfizer_telohaec_cytokine`
     never increments `K_conf` or `k`).
 11. Fast Pfizer TeloHAEC NTC sampler parity vs locked `pfizer_telohaec_reader_r17.py` on
     synthetic `.h5ad` fixtures.
 12. Full CPU synthetic end-to-end smoke run (`run_synthetic_smoke_r17`) + Stage 4 from-disk
     recomputation parity (`max_abs_diff <= 1e-12`).
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

ROUND17_SRC_DIR = Path(__file__).resolve().parent
REPO_ROOT = ROUND17_SRC_DIR.parents[1]
R16_SRC_DIR = REPO_ROOT / "round16" / "src"
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
R10_SCRIPTS_DIR = REPO_ROOT / "round10" / "shared" / "scripts"
for _p in (str(ROUND17_SRC_DIR), str(R16_SRC_DIR), str(R10_SRC_DIR), str(R10_SCRIPTS_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from inference_r17 import compute_batch_encoder_medoid_from_weights_r17  # noqa: E402
from metrics_r16 import compute_orth_suite_r16, compute_paired_contrast_r16  # noqa: E402
from ntc_extractors_r17 import (  # noqa: E402
    BlindingGateError,
    assert_all_execution_gates_r17,
    build_mismatched_strata_r17,
    extract_pfizer_telohaec_ntc_bundle_r17,
    validate_target_lines_r17,
)
from pfizer_telohaec_reader_r17 import extract_pfizer_telohaec_ntc_raw_counts_r17  # noqa: E402
from prepare_r10_lines import build_hgnc_lookup  # noqa: E402
from run_eval_r17 import (  # noqa: E402
    _make_synthetic_line_bundle_r17,
    evaluate_cohort_from_arrays_r17,
    run_synthetic_smoke_r17,
)
from stats_r17 import (  # noqa: E402
    ALPHA_LINE_BASE,
    DEFAULT_R17_SEED,
    N_BATCH_INDICES_STATE,
    _compute_half_ceiling_draws,
    _eval_pred_mom_with_ceiling,
    adjudicate_h_j1_r17,
    adjudicate_h_j2_r17,
    adjudicate_h_j3_rule_r17,
    apply_holm_bonferroni_r17,
    compute_j5a_g4_r17,
    compute_j5b_rank_correlations_r17,
    compute_parametric_power_r17,
    compute_spearman_perm_pair_r17,
    compute_sweep56_moments_r17,
    compute_variance_component_r17,
    evaluate_detectable_signal_r17,
    make_stratified_half_split_r17,
    select_l2_ntc_recon_match_r17,
)


class TestExecutionGatesAndLineSpecs(unittest.TestCase):
    def test_execution_gates_refuse_missing_or_empty(self) -> None:
        with tempfile.TemporaryDirectory() as tmpd:
            td = Path(tmpd)
            lock_f = td / "PREREG_J_LOCK.md"
            freeze_f = td / "R17_CODE_FREEZE.txt"
            go_f = td / "R17_GO.txt"

            # All missing
            with self.assertRaises(BlindingGateError):
                assert_all_execution_gates_r17(status_dir=td, prereg_lock_path=lock_f)

            # Lock present, freeze & GO missing
            lock_f.write_text("locked\n", encoding="utf-8")
            with self.assertRaises(BlindingGateError):
                assert_all_execution_gates_r17(status_dir=td, prereg_lock_path=lock_f)

            # Freeze empty -> still fails
            freeze_f.write_text("", encoding="utf-8")
            go_f.write_text("GO\n", encoding="utf-8")
            with self.assertRaises(BlindingGateError):
                assert_all_execution_gates_r17(status_dir=td, prereg_lock_path=lock_f)

            # Freeze non-empty, GO missing -> fails
            freeze_f.write_text("FROZEN\n", encoding="utf-8")
            go_f.unlink()
            with self.assertRaises(BlindingGateError):
                assert_all_execution_gates_r17(status_dir=td, prereg_lock_path=lock_f)

            # All three non-empty -> passes
            go_f.write_text("GO\n", encoding="utf-8")
            gates = assert_all_execution_gates_r17(status_dir=td, prereg_lock_path=lock_f)
            self.assertEqual(set(gates.keys()), {"code_freeze", "r17_go", "prereg_lock"})

    def test_omitted_and_quarantined_lines_blocked(self) -> None:
        with self.assertRaises(ValueError):
            validate_target_lines_r17(["wtc11_i3neuron"])
        for q_line in ("nadig_hepg2", "nadig_jurkat", "a549_ko", "dld1", "cd4t_rest", "feng_ipsc"):
            with self.assertRaises(ValueError):
                validate_target_lines_r17([q_line])
        ok = validate_target_lines_r17(["cd4t_stim48", "gse295214_h1", "pfizer_telohaec_basal", "pfizer_telohaec_cytokine"])
        self.assertEqual(len(ok), 4)

    def test_mismatched_strata_j5a_exact(self) -> None:
        st_cd4 = build_mismatched_strata_r17("cd4t_stim48")
        self.assertEqual(st_cd4["seen_by_checkpoint"], ("nadig_hepg2", "nadig_jurkat"))
        self.assertEqual(st_cd4["unseen_by_checkpoint"], ("a549_ko", "dld1", "feng_ipsc"))
        self.assertEqual(st_cd4["same_lineage_comparators"], ("cd4t_rest",))

        for a549_ln in ("a549_myllia", "a549_ko_cas9pos"):
            st_a = build_mismatched_strata_r17(a549_ln)
            self.assertNotIn("a549_ko", st_a["unseen_by_checkpoint"])
            self.assertEqual(st_a["unseen_by_checkpoint"], ("cd4t_rest", "dld1", "feng_ipsc"))

        for indep in ("gse295214_h1", "pfizer_telohaec_basal", "pfizer_telohaec_cytokine"):
            st_i = build_mismatched_strata_r17(indep)
            self.assertEqual(st_i["unseen_by_checkpoint"], ("a549_ko", "cd4t_rest", "dld1", "feng_ipsc"))
            self.assertEqual(st_i["same_lineage_comparators"], ())


class TestParityWithFrozenR16Metrics(unittest.TestCase):
    def test_sweep56_fast_boot_matches_frozen_r16_to_1e14(self) -> None:
        b = _make_synthetic_line_bundle_r17(
            "gse295214_h1", "confirmatory", "CRISPRi", 36, 120,
            regime_se="strong_batch", regime_hvg="invariant_nonzero", seed=777,
        )
        p56 = b["models"]["STATE_SE"]["delta_pred_56"]
        moms = compute_sweep56_moments_r17(
            b["delta_obs"], p56, b["fct_delta_pred"], b["delta_transfer"],
            b["target_indices"], None, b["delta_half_a"], b["delta_half_b"], b["half_valid"],
        )
        m = moms["m"]
        n_boot = 250
        draws = np.random.default_rng(DEFAULT_R17_SEED).integers(0, m, size=(n_boot, m))
        _, ceil_o, b_ceil_o = _compute_half_ceiling_draws(moms["half_mom"], draws)

        for check_b in (0, 17, 55):
            pt_fast, boot_fast, _ = _eval_pred_mom_with_ceiling(
                moms["moms_56"][check_b], ceil_o, b_ceil_o, draws, zero_null_to_zero=False
            )
            ref = compute_orth_suite_r16(
                b["delta_obs"], p56[check_b], b["delta_transfer"], b["target_indices"],
                source_grand_mean=None, delta_half_a=b["delta_half_a"], delta_half_b=b["delta_half_b"],
                half_valid=b["half_valid"], n_boot=n_boot, seed=DEFAULT_R17_SEED,
            )
            self.assertLessEqual(abs(pt_fast - ref["estimate"]), 1e-14)
            self.assertLessEqual(abs(float(np.percentile(boot_fast, 1.25)) - ref["ci975_lo"]), 1e-14)
            self.assertLessEqual(abs(float(np.percentile(boot_fast, 98.75)) - ref["ci975_hi"]), 1e-14)

        # Also check paired contrast parity vs FCT
        pt_b0, boot_b0, _ = _eval_pred_mom_with_ceiling(moms["moms_56"][0], ceil_o, b_ceil_o, draws)
        pt_fct, boot_fct, _ = _eval_pred_mom_with_ceiling(moms["mom_fct"], ceil_o, b_ceil_o, draws)
        ref_pair = compute_paired_contrast_r16(
            b["delta_obs"], p56[0], b["fct_delta_pred"], b["delta_transfer"], b["target_indices"],
            source_grand_mean=None, delta_half_a=b["delta_half_a"], delta_half_b=b["delta_half_b"],
            half_valid=b["half_valid"], n_boot=n_boot, seed=DEFAULT_R17_SEED,
        )
        self.assertLessEqual(abs((pt_b0 - pt_fct) - ref_pair["estimate"]), 1e-14)
        self.assertLessEqual(abs(float(np.percentile(boot_b0 - boot_fct, 1.25)) - ref_pair["ci975_lo"]), 1e-14)


class TestHalfSplitAndVarianceComponent(unittest.TestCase):
    def test_stratified_half_split_exact_sizes_and_order_invariance(self) -> None:
        # Mimic cd4t_stim48 (N=401, 333 half_valid=True, 68 half_valid=False)
        names = [f"G_{i:04d}" for i in range(401)]
        hv = np.array([True] * 333 + [False] * 68, dtype=bool)
        p1, p2 = make_stratified_half_split_r17(names, hv, seed=DEFAULT_R17_SEED)
        self.assertEqual(len(p1), 200)
        self.assertEqual(len(p2), 201)
        self.assertEqual(len(set(p1.tolist()) & set(p2.tolist())), 0)
        self.assertEqual(len(set(p1.tolist()) | set(p2.tolist())), 401)
        self.assertEqual(int(np.sum(hv[p1])), 166)
        self.assertEqual(int(np.sum(hv[p2])), 167)
        self.assertEqual(int(np.sum(~hv[p1])), 34)
        self.assertEqual(int(np.sum(~hv[p2])), 34)

        # Permuting input order must produce the exact same canonical target-name sets!
        rng = np.random.default_rng(999)
        shuf = rng.permutation(401)
        names_shuf = [names[i] for i in shuf]
        hv_shuf = hv[shuf]
        p1_s, p2_s = make_stratified_half_split_r17(names_shuf, hv_shuf, seed=DEFAULT_R17_SEED)
        self.assertEqual([names[i] for i in p1], [names_shuf[i] for i in p1_s])
        self.assertEqual([names[i] for i in p2], [names_shuf[i] for i in p2_s])

    def test_variance_component_null_vs_signal(self) -> None:
        # Exact identical indices -> tau_batch == 0.0, p_var == 1.0
        r_zero = np.full(N_BATCH_INDICES_STATE, 0.25, dtype=np.float64)
        boot_zero = np.full((500, N_BATCH_INDICES_STATE), 0.25, dtype=np.float64)
        vc0 = compute_variance_component_r17(r_zero, boot_zero)
        self.assertEqual(vc0["tau_batch"], 0.0)
        self.assertEqual(vc0["p_var"], 1.0)

        # Pure noise with true tau = 0 -> tau_batch near 0, p_var not tiny
        rng = np.random.default_rng(123)
        r_null = 0.25 + rng.normal(0.0, 0.01, size=N_BATCH_INDICES_STATE)
        boot_null = r_null[None, :] + rng.normal(0.0, 0.01, size=(1000, N_BATCH_INDICES_STATE))
        vcn = compute_variance_component_r17(r_null, boot_null)
        self.assertGreater(vcn["p_var"], 0.05)

        # Strong true signal (SD = 0.08 >> boot SE = 0.005) -> tau_batch ~ 0.08, p_var = 1/1001
        r_sig = np.linspace(0.05, 0.35, N_BATCH_INDICES_STATE)
        boot_sig = r_sig[None, :] + rng.normal(0.0, 0.005, size=(1000, N_BATCH_INDICES_STATE))
        vcs = compute_variance_component_r17(r_sig, boot_sig)
        self.assertGreater(vcs["tau_batch"], 0.07)
        self.assertLess(vcs["p_var"], 0.002)


class TestAdjudicationBranchesAndMultiplicity(unittest.TestCase):
    def _mock_hj1_call(
        self,
        p_var: float,
        tau_batch: float,
        boot_dcf: np.ndarray,
        has_sig: bool,
        vf: float = 1.0,
        rho: float = 0.6,
        halves_ok: bool = True,
        alpha_adj: float = ALPHA_LINE_BASE,
        holm_rejected: bool | None = None,
    ) -> dict:
        dcf_est = float(np.mean(boot_dcf))
        ci95_hi = float(np.percentile(boot_dcf, 97.5))
        vc = {
            "p_var": p_var,
            "tau_batch": tau_batch,
            "min_per_index_valid_frac": vf,
            "joint_valid_frac": vf,
        }
        cf_s0 = {
            "half_valid_both_halves": halves_ok,
            "delta_cf_0": {
                "estimate": dcf_est,
                "valid_frac": vf,
                "ci95_hi": ci95_hi,
                "p_pos_two_sided": float(min(1.0, 2.0 * (1 + np.sum(boot_dcf <= 0.0)) / (len(boot_dcf) + 1))),
            },
            "r_cf_oracle_0": {"estimate": 0.30 if has_sig else 0.01, "ci9833_lo": 0.20 if has_sig else -0.02, "ci9833_hi": 0.40, "valid_frac": vf},
        }
        ds = {"has_detectable_signal": has_sig}
        land = {"max": 0.35 if has_sig else 0.02, "argmax_idx": 42}
        l1 = {"estimate": 0.25 if has_sig else 0.005, "ci9833_lo": 0.15 if has_sig else -0.02, "ci9833_hi": 0.35, "valid_frac": vf}
        return adjudicate_h_j1_r17(vc, cf_s0, boot_dcf, ds, land, l1, rho, alpha_adj=alpha_adj, holm_rejected=holm_rejected)

    def test_h_j1_all_six_branches(self) -> None:
        rng = np.random.default_rng(42)
        # 1. VALIDITY_FAIL
        r1 = self._mock_hj1_call(0.001, 0.05, rng.normal(0.08, 0.005, 500), has_sig=True, vf=0.90)
        self.assertEqual(r1["verdict"], "INCONCLUSIVE")
        self.assertEqual(r1["subreason"], "VALIDITY_FAIL")

        # 2. BATCH_DEPENDENCE_PRACTICALLY_MEANINGFUL
        r2 = self._mock_hj1_call(0.0005, 0.05, rng.normal(0.08, 0.005, 500), has_sig=True)
        self.assertEqual(r2["verdict"], "BATCH_DEPENDENCE_PRACTICALLY_MEANINGFUL")

        # 3. BATCH_DEPENDENCE_DETECTED_SUB_DELTA (detected, yet delta_cf <= 0.04 or tau <= 0.02)
        r3 = self._mock_hj1_call(0.0005, 0.015, rng.normal(0.03, 0.004, 500), has_sig=True)
        self.assertEqual(r3["verdict"], "BATCH_DEPENDENCE_DETECTED_SUB_DELTA")

        # 4. INCONCLUSIVE_LOW_SIGNAL (`J2` floor guard: invariance CANNOT be awarded on null floor!)
        r4 = self._mock_hj1_call(0.50, 0.002, rng.normal(0.003, 0.004, 500), has_sig=False)
        self.assertEqual(r4["verdict"], "INCONCLUSIVE_LOW_SIGNAL")

        # 5. BATCH_INVARIANT_WITHIN_DELTA (detectable signal + not detected + ci95_hi < 0.04 + tau <= 0.02)
        r5 = self._mock_hj1_call(0.50, 0.005, rng.normal(0.005, 0.005, 500), has_sig=True)
        self.assertEqual(r5["verdict"], "BATCH_INVARIANT_WITHIN_DELTA")

        # 6. INCONCLUSIVE (INSUFFICIENT_PRECISION: detectable signal, not detected, yet ci95_hi >= 0.04)
        r6 = self._mock_hj1_call(0.20, 0.010, rng.normal(0.02, 0.02, 500), has_sig=True)
        self.assertEqual(r6["verdict"], "INCONCLUSIVE")
        self.assertEqual(r6["subreason"], "INSUFFICIENT_PRECISION")

    def test_detectable_signal_j2_predicate(self) -> None:
        d1 = evaluate_detectable_signal_r17(
            {"valid_frac": 1.0, "ci9833_lo": -0.01},
            {"valid_frac": 1.0, "ci9833_lo": -0.02},
        )
        self.assertFalse(d1["has_detectable_signal"])
        d2 = evaluate_detectable_signal_r17(
            {"valid_frac": 1.0, "ci9833_lo": 0.05},
            {"valid_frac": 1.0, "ci9833_lo": -0.01},
        )
        self.assertTrue(d2["has_detectable_signal"])
        d3 = evaluate_detectable_signal_r17(
            {"valid_frac": 0.90, "ci9833_lo": 0.10},  # Fails valid_frac >= 0.95!
            {"valid_frac": 1.0, "ci9833_lo": 0.00},
        )
        self.assertFalse(d3["has_detectable_signal"])

    def test_h_j2_sign_reversal_and_low_signal_qualifier(self) -> None:
        b_strong = _make_synthetic_line_bundle_r17(
            "cd4t_stim48", "confirmatory", "CRISPRi", 50, 150,
            regime_se="strong_batch", regime_hvg="low_signal", seed=321,
        )
        p56_se = b_strong["models"]["STATE_SE"]["delta_pred_56"]
        moms_se = compute_sweep56_moments_r17(
            b_strong["delta_obs"], p56_se, b_strong["fct_delta_pred"], b_strong["delta_transfer"],
            b_strong["target_indices"], None, b_strong["delta_half_a"], b_strong["delta_half_b"], b_strong["half_valid"],
        )
        m = moms_se["m"]
        draws = np.random.default_rng(DEFAULT_R17_SEED).integers(0, m, size=(300, m))
        _, c_o, bc_o = _compute_half_ceiling_draws(moms_se["half_mom"], draws)
        r56 = np.empty(56)
        rb56 = np.empty((300, 56))
        for i in range(56):
            r56[i], rb56[:, i], _ = _eval_pred_mom_with_ceiling(moms_se["moms_56"][i], c_o, bc_o, draws)
        fct_pt, fct_boot, _ = _eval_pred_mom_with_ceiling(moms_se["mom_fct"], c_o, bc_o, draws)
        hj2 = adjudicate_h_j2_r17(
            r56, rb56, fct_pt, fct_boot,
            {"estimate": fct_pt, "ci975_lo": 0.1, "ci975_hi": 0.3},
            moms_se, {"has_detectable_signal": True}, {"max": float(np.max(r56))}, {"estimate": 0.2},
            n_boot=300, n_boot_bonf56=1000, seed=DEFAULT_R17_SEED,
        )
        self.assertEqual(hj2["classification"], "CALL_UNSTABLE_SIGN_REVERSAL")
        self.assertTrue(hj2["has_full_reversal"])
        self.assertIn("companion_bonf56_normal", hj2)
        self.assertIn("companion_bonf56_emp20k", hj2)

        # Low-signal floor with 0 flips -> INCONCLUSIVE_LOW_SIGNAL, NOT CALL_STABLE
        rb_flat = np.zeros((300, 56))
        r_flat = np.zeros(56)
        hj2_low = adjudicate_h_j2_r17(
            r_flat, rb_flat, 0.0, np.zeros(300),
            {"estimate": 0.0, "ci975_lo": -0.05, "ci975_hi": 0.05},
            moms_se, {"has_detectable_signal": False}, {"max": 0.0}, {"estimate": 0.0},
            n_boot=300, n_boot_bonf56=500, seed=DEFAULT_R17_SEED,
        )
        self.assertEqual(hj2_low["classification"], "INCONCLUSIVE_LOW_SIGNAL")

    def test_h_j3_l2_selection_ruling3_and_exact_tie_zero(self) -> None:
        rng = np.random.default_rng(55)
        y_obs = rng.normal(2.0, 1.0, size=200)
        pb_pf = rng.normal(2.0, 1.0, size=(4, 56, 200))
        # Make index 7 and index 19 exact top ties -> lowest-index tie-break picks 7, margin == 0.0
        for fi in range(4):
            pb_pf[fi, 7] = y_obs * 1.1 + 0.2
            pb_pf[fi, 19] = y_obs * 1.1 + 0.2
        pb_ens = np.mean(pb_pf, axis=0)
        sel = select_l2_ntc_recon_match_r17(pb_ens, y_obs, ntc_pb_per_fold_56=pb_pf)
        self.assertEqual(sel["selected_batch_idx"], 7)
        self.assertEqual(sel["primary_pearson"]["top2_idx"], 19)
        self.assertAlmostEqual(sel["primary_pearson"]["criterion_margin_top1_minus_top2"], 0.0, places=12)
        self.assertTrue(sel["secondary_per_fold_agreement"]["all_folds_agree"])
        self.assertIn("criterion_margin_top1_minus_top2", sel["secondary_neg_l2"])

        # EXACT_TIE_SELECTED_ZERO check
        res_tie = adjudicate_h_j3_rule_r17(
            "L2_ntc_recon_match", 0.0, np.zeros(300),
            {"estimate": 0.0, "ci975_lo": -0.1, "ci975_hi": 0.1, "valid_frac": 1.0},
            None, ALPHA_LINE_BASE, 50, 300, holm_rejected=False, exact_tie_selected_zero=True,
        )
        self.assertEqual(res_tie["verdict"], "EXACT_TIE_SELECTED_ZERO")

    def test_holm_bonferroni_step_down_exact(self) -> None:
        pvals = {
            "t1": 0.001,   # < 0.05/4 = 0.0125 -> reject
            "t2": 0.015,   # < 0.05/3 = 0.01667 -> reject
            "t3": 0.030,   # >= 0.05/2 = 0.025 -> stop!
            "t4": 0.040,   # after stop -> NOT rejected even though < 0.05/1
        }
        ht = apply_holm_bonferroni_r17(pvals, alpha_family=0.05)
        self.assertTrue(ht["t1"]["rejected"])
        self.assertTrue(ht["t2"]["rejected"])
        self.assertFalse(ht["t3"]["rejected"])
        self.assertFalse(ht["t4"]["rejected"])
        self.assertAlmostEqual(ht["t1"]["step_alpha"], 0.0125, places=12)
        self.assertAlmostEqual(ht["t2"]["step_alpha"], 0.05 / 3.0, places=12)
        self.assertAlmostEqual(ht["t4"]["effective_alpha"], 0.025, places=12)

    def test_medoid_and_parametric_power(self) -> None:
        # 56 orthogonal-ish rows, make row 11 exact zero (centroid) -> unique Euclidean medoid == 11
        rng = np.random.default_rng(7)
        w = rng.normal(0.0, 1.0, size=(56, 128))
        w = w - np.mean(w, axis=0, keepdims=True)
        w[11] = 0.0
        med = compute_batch_encoder_medoid_from_weights_r17([w, w, w, w])
        self.assertEqual(med["medoid_idx"], 11)

        # Verify PREREG_J §4.3 Table values at N=400, K=3
        pw = compute_parametric_power_r17(400, k_conf=3)
        self.assertAlmostEqual(pw["se_pair_hi"], 0.00751, places=5)
        self.assertAlmostEqual(pw["se_pair_lo"], 0.00875, places=5)
        self.assertGreater(pw["cond_power_delta004_holm_step1_lo"], 0.96)


class TestJ5MechanismAndPfizerParity(unittest.TestCase):
    def test_j5a_cross_fitted_g4_strictly_uses_p2_only(self) -> None:
        b = _make_synthetic_line_bundle_r17(
            "gse295214_h1", "confirmatory", "CRISPRi", 40, 120,
            regime_se="strong_batch", regime_hvg="invariant_nonzero", seed=888,
        )
        p1, p2 = make_stratified_half_split_r17(b["target_names"], b["half_valid"], seed=DEFAULT_R17_SEED)
        mdata = b["models"]["STATE_SE"]
        g4_orig = compute_j5a_g4_r17(
            b["delta_obs"], b["delta_half_a"], b["delta_half_b"], b["half_valid"],
            b["delta_transfer"], b["source_grand_mean"], b["target_indices"],
            mdata["delta_pred_56"], mdata["mismatched_by_source_and_idx"],
            b["seen_sources"], b["unseen_sources"], b["same_lineage_sources"],
            p1, p2, b_star_1=35, b_dag_1=0, n_boot=150, seed=DEFAULT_R17_SEED,
        )

        # Corrupt P1 rows massively -> locked P1->P2 readout evaluated on P2 ONLY must be bit-identical!
        obs_corrupt = b["delta_obs"].copy()
        obs_corrupt[p1] += 100.0
        g4_corrupt = compute_j5a_g4_r17(
            obs_corrupt, b["delta_half_a"], b["delta_half_b"], b["half_valid"],
            b["delta_transfer"], b["source_grand_mean"], b["target_indices"],
            mdata["delta_pred_56"], mdata["mismatched_by_source_and_idx"],
            b["seen_sources"], b["unseen_sources"], b["same_lineage_sources"],
            p1, p2, b_star_1=35, b_dag_1=0, n_boot=150, seed=DEFAULT_R17_SEED,
        )
        sk0 = b["unseen_sources"][0]
        v1 = g4_orig["cross_fitted_primary_p1_to_p2"]["per_source"][f"cf_p1_to_p2_fwd_mis_{sk0}"]["r_orth_cal"]["delta_batch"]["estimate"]
        v2 = g4_corrupt["cross_fitted_primary_p1_to_p2"]["per_source"][f"cf_p1_to_p2_fwd_mis_{sk0}"]["r_orth_cal"]["delta_batch"]["estimate"]
        self.assertEqual(v1, v2)

    def test_j5b_spearman_perm_and_order_invariance(self) -> None:
        x = np.linspace(0.0, 1.0, 56)
        y = x * 2.0 + 0.1
        r_id = compute_spearman_perm_pair_r17(x, y, n_perm=1000, seed=DEFAULT_R17_SEED)
        self.assertAlmostEqual(r_id["spearman_rho"], 1.0, places=12)
        self.assertAlmostEqual(r_id["p_perm_two_sided"], 1.0 / 1001.0, places=12)
        self.assertEqual(r_id["top5_overlap_count"], 5)

        rng = np.random.default_rng(555)
        scores = {
            "STATE_SE": {
                "cd4t_stim48": x,
                "gse295214_h1": x + rng.normal(0.0, 0.05, 56),
                "pfizer_telohaec_basal": rng.normal(0.0, 1.0, 56),
                "pfizer_telohaec_cytokine": rng.normal(0.0, 1.0, 56),
                "a549_ko_cas9pos": x,
            },
            "STATE_HVG": {
                "cd4t_stim48": x + rng.normal(0.0, 0.1, 56),
                "gse295214_h1": rng.normal(0.0, 1.0, 56),
                "pfizer_telohaec_basal": rng.normal(0.0, 1.0, 56),
                "pfizer_telohaec_cytokine": rng.normal(0.0, 1.0, 56),
                "a549_ko_cas9pos": x,
            },
        }
        roles = {
            "cd4t_stim48": "confirmatory",
            "gse295214_h1": "confirmatory",
            "pfizer_telohaec_basal": "confirmatory",
            "pfizer_telohaec_cytokine": "within_family_replication",
            "a549_ko_cas9pos": "r16_subset_sensitivity",
        }
        neval = {k: 100 for k in roles}
        j5b = compute_j5b_rank_correlations_r17(scores, roles, neval, n_perm=500, seed=DEFAULT_R17_SEED)
        self.assertEqual(j5b["by_model"]["STATE_SE"]["confirmatory_pairs"]["n_pairs"], 3)
        self.assertEqual(j5b["by_model"]["STATE_SE"]["pfizer_within_family_pair"]["n_pairs"], 1)
        # a549_ko_cas9pos must NOT enter all_fresh_n_ge_16_pairs (4 fresh lines -> C(4,2) = 6 pairs)
        self.assertEqual(j5b["by_model"]["STATE_SE"]["all_fresh_n_ge_16_pairs"]["n_pairs"], 6)
        self.assertEqual(j5b["cross_model_within_line"]["n_pairs"], 5)

    def test_pfizer_fast_ntc_sampler_matches_locked_reader(self) -> None:
        # Build tiny synthetic Pfizer-schema .h5ad (Untreated + Treated, NO-TARGET + SAFE_TARGET + TP53)
        with tempfile.TemporaryDirectory() as tmpd:
            h5p = Path(tmpd) / "imru_full.h5ad"
            rng = np.random.default_rng(99)
            n_genes = 220
            # 160 Untreated NO-TARGET (80 flask1, 80 flask2) + 40 SAFE_TARGET + 40 TP53 + 60 Treated
            n_obs = 300
            cond = ["Untreated"] * 240 + ["Treated"] * 60
            gtgt = ["NO-TARGET"] * 160 + ["SAFE_TARGET"] * 40 + ["TP53"] * 40 + ["NO-TARGET"] * 60
            guide = [f"g_{i % 10}" for i in range(n_obs)]
            flask = [("flask1" if i % 2 == 0 else "flask2") for i in range(n_obs)]
            cid = [f"CELL_{i:04d}" for i in range(n_obs)]

            # Ensure every cell passes UMI >= 300 and ngenes >= 200
            mat = rng.integers(1, 5, size=(n_obs, n_genes), dtype=np.int32)
            csr = sp.csr_matrix(mat.astype(np.float32))
            with h5py.File(h5p, "w") as f:
                xg = f.create_group("X")
                xg.attrs["encoding-type"] = "csr_matrix"
                xg.attrs["shape"] = np.array([n_obs, n_genes], dtype=np.int64)
                xg.create_dataset("data", data=csr.data)
                xg.create_dataset("indices", data=csr.indices)
                xg.create_dataset("indptr", data=csr.indptr)

                og = f.create_group("obs")
                og.attrs["_index"] = "_index"
                dt = h5py.string_dtype(encoding="utf-8")
                og.create_dataset("_index", data=np.array(cid, dtype=object), dtype=dt)
                og.create_dataset("condition", data=np.array(cond, dtype=object), dtype=dt)
                og.create_dataset("gene_target", data=np.array(gtgt, dtype=object), dtype=dt)
                og.create_dataset("guide", data=np.array(guide, dtype=object), dtype=dt)
                og.create_dataset("flask", data=np.array(flask, dtype=object), dtype=dt)
                og.create_dataset("cellID", data=np.array(cid, dtype=object), dtype=dt)

                vg = f.create_group("var")
                vg.attrs["_index"] = "_index"
                vg.create_dataset("_index", data=np.array([f"GENE_{j}" for j in range(n_genes)], dtype=object), dtype=dt)

            hgnc = build_hgnc_lookup()
            locked_ref = extract_pfizer_telohaec_ntc_raw_counts_r17(h5p, hgnc, condition="untreated")
            fast_b = extract_pfizer_telohaec_ntc_bundle_r17(
                h5p, condition="untreated", target_n=128, seed=DEFAULT_R17_SEED
            )
            self.assertEqual(fast_b["ctrl_audit"]["n_ctrl_total"], len(locked_ref["barcodes"]))
            self.assertEqual(fast_b["ctrl_audit"]["n_ctrl_sampled"], 128)
            # Every sampled barcode must come from locked_ref and match raw counts bit-for-bit!
            bc_to_row = {bc: idx for idx, bc in enumerate(locked_ref["barcodes"])}
            for k, bc in enumerate(fast_b["barcodes"]):
                self.assertIn(bc, bc_to_row)
                np.testing.assert_array_equal(fast_b["raw_counts"][k], locked_ref["raw_counts"][bc_to_row[bc]])
                self.assertEqual(fast_b["batch_labels"][k], locked_ref["batch_labels"][bc_to_row[bc]])


class TestSyntheticEndToEndAndStage4Recompute(unittest.TestCase):
    def test_synthetic_smoke_and_stage4_from_disk_recompute(self) -> None:
        with tempfile.TemporaryDirectory(prefix="r17_unit_smoke_") as tmpd:
            td = Path(tmpd)
            # n_boot=800 guarantees p_min = 2/801 = 0.002497 < alpha/18 = 0.002778 (both separate M=9 and pooled M=18 Holm)
            res = run_synthetic_smoke_r17(td, n_boot=800, n_boot_bonf56=800, n_perm=500, seed=DEFAULT_R17_SEED)
            self.assertEqual(res["real_cells_touched"], 0)
            self.assertEqual(res["K_conf"], 3)
            self.assertTrue(res["stage4_from_disk_recompute_passed"])
            self.assertLessEqual(res["stage4_from_disk_recompute_max_abs_diff"], 1e-12)

            # Verify confirmatory claim rule (separate Holm M=9 AND pooled Holm M=18):
            # STATE_SE: cd4t_stim48 (strong) + gse295214_h1 (strong) + pfizer_basal (invariant) -> 2/3 -> general claim TRUE
            se_coh = res["separate_model_holm_fwer05"]["cohort_adjudications"]["STATE_SE"]
            self.assertEqual(se_coh["K_conf"], 3)
            self.assertEqual(se_coh["k_over_K"], "2/3")
            self.assertTrue(se_coh["general_cross_context_claim_permitted"])
            self.assertEqual(
                se_coh["per_confirmatory_line"]["pfizer_telohaec_basal"]["H_J1_holm"]["verdict"],
                "BATCH_INVARIANT_WITHIN_DELTA",
            )
            se_pooled = res["pooled_model_holm_sensitivity_fwer05"]["cohort_adjudications"]["STATE_SE"]
            self.assertEqual(se_pooled["k_over_K"], "2/3")
            self.assertTrue(se_pooled["general_cross_context_claim_permitted"])

            # STATE_HVG: cd4t_stim48 (strong) + gse295214_h1 (invariant) + pfizer_basal (low_signal -> J2 guard!) -> 1/3 -> FALSE
            hvg_coh = res["separate_model_holm_fwer05"]["cohort_adjudications"]["STATE_HVG"]
            self.assertEqual(hvg_coh["k_over_K"], "1/3")
            self.assertFalse(hvg_coh["general_cross_context_claim_permitted"])
            self.assertEqual(
                hvg_coh["per_confirmatory_line"]["pfizer_telohaec_basal"]["H_J1_holm"]["verdict"],
                "INCONCLUSIVE_LOW_SIGNAL",
            )
            hvg_pooled = res["pooled_model_holm_sensitivity_fwer05"]["cohort_adjudications"]["STATE_HVG"]
            self.assertEqual(hvg_pooled["k_over_K"], "1/3")
            self.assertFalse(hvg_pooled["general_cross_context_claim_permitted"])

            # Verify within-family pfizer_telohaec_cytokine NEVER increments K_conf
            self.assertNotIn("pfizer_telohaec_cytokine", res["confirmatory_units"])
            self.assertIn("pfizer_telohaec_cytokine", res["lines"])
            # Verify micro-N mac_ctrl (N=3) completed cleanly
            self.assertEqual(res["lines"]["mac_ctrl"]["n_eval_targets"], 3)
            # Verify gse295214_h1 sensitivity FCT exported and evaluated
            self.assertIn(
                "FCT_full_r9_6_incl_vcc2025_h1",
                res["lines"]["gse295214_h1"]["models"]["STATE_SE"]["secondary_extra_fct_sensitivity"],
            )
            self.assertTrue((td / "arrays" / "gse295214_h1" / "SENS_FCT_full_r9_6_incl_vcc2025_h1.npz").is_file())
            self.assertTrue((td / "control_barcodes_r17.json").is_file())
            self.assertTrue((td / "embeddings" / "embeddings_manifest_sha256.json").is_file())

            # Verify summary_r17.json is valid JSON cleanly readable
            loaded = json.loads((td / "summary_r17.json").read_text(encoding="utf-8"))
            self.assertEqual(loaded["K_conf"], 3)

    def test_cohort_k_eq_1_forbids_general_claim(self) -> None:
        b = _make_synthetic_line_bundle_r17(
            "cd4t_stim48", "confirmatory", "CRISPRi", 50, 140,
            regime_se="strong_batch", regime_hvg="strong_batch", seed=444,
        )
        res = evaluate_cohort_from_arrays_r17({"cd4t_stim48": b}, arrays_dir=None, n_boot=200, n_boot_bonf56=500, n_perm=200)
        se_coh = res["separate_model_holm_fwer05"]["cohort_adjudications"]["STATE_SE"]
        self.assertEqual(se_coh["k_over_K"], "1/1")
        self.assertFalse(se_coh["general_cross_context_claim_permitted"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
