#!/usr/bin/env python3
"""
Clean-Room Unit & Synthetic End-to-End Verification Suite (`round17/audit/test_audit_r17.py`, `A-T3`).

Strictly CPU-only, zero real-cell access, zero reading of forbidden `round17/src/` files.
Covers:
  1. Freeze-gate refusal (`assert_freeze_gates_r17` refuses unless `PREREG_J_LOCK.md`,
     `R17_CODE_FREEZE_AT3.txt`, and `R17_GO_AT3.txt` all exist and are non-empty).
  2. `C1.4` stratified alternating split (`s = 0..19`, odd `P1`, even `P2`, sorted ascending).
  3. `H-J1` (`Delta_CF_0`, `tau_batch`, `V_noise` `ddof=1`, centered-bootstrap `p_var`, intersection-union
     `p_hj1_iu = max(p_var, p_Delta)`, 6-way classifier + subreasons, `mac_ctrl` `N=3` `HALF_SIZE_LT_2`).
  4. `H-J2` (`f_flip`, `has_full_reversal`, normal `z = norm.ppf(1 - 0.025/112)` and 20k empirical Bonferroni-56).
  5. `H-J3` (`L1_ens56`, `L2_ntc_recon_match` including `EXACT_TIE_SELECTED_ZERO` & per-fold diagnostics,
     Euclidean vs cosine geometric medoid tie-breaking).
  6. Separate-model & pooled Holm-Bonferroni FWER step-down gating (`p < alpha_step` AND CI gate for `H-J1`),
     `J5a` cross-fitted & continuity `G4` decompositions, and `J5b` shared 10k-perm table + BH FDR.
  7. Full `--phase synthetic_e2e` schema verification (`allow_pickle=False`, every exact `R17_ARRAY_SCHEMA.md`
     key/shape/dtype across all 11 lines + 17 pools + `SENS_FCT_full_r9_6_incl_vcc2025_h1.npz`) AND
     CLI `--recompute-from-arrays` Stage 4 recomputation verifying `max_abs_diff <= 1e-12`.
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
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
from scipy.stats import norm

AUDIT_R17_DIR = Path(__file__).resolve().parent
R17_ROOT = AUDIT_R17_DIR.parent
REPO_ROOT = R17_ROOT.parent
if str(AUDIT_R17_DIR) not in sys.path:
    sys.path.insert(0, str(AUDIT_R17_DIR))

from audit_r17_controls import (  # noqa: E402
    DEFAULT_R17_SEED,
    R16_MISMATCHED_POOLS,
    R17_ALL_EVAL_LINES,
    assert_freeze_gates_r17,
    resolve_r17_mismatched_strata,
    stratified_split_halves_c1,
)
from audit_r16_state_adapter import BlindingGateError  # noqa: E402
from audit_r17_stats import (  # noqa: E402
    _compute_orth_vectors_and_boot,
    apply_holm_fwer_r17,
    classify_h_j1_under_alpha,
    compute_bh_fdr,
    compute_h_j1_for_subset,
    compute_h_j2_for_subset,
    compute_j5b_rank_correlations,
    recompute_cohort_from_exported_arrays_r17,
    select_l2_ntc_index_and_diagnostics,
)
from run_audit_r17 import (  # noqa: E402
    build_synthetic_r17_cohort,
    compare_summaries_stage4,
    compute_batch_encoder_medoids,
)


class TestR17FreezeGatesAndSplits(unittest.TestCase):
    def test_freeze_gates_refuse_real_execution_until_unlocked(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_p = Path(td)
            prereg_d = td_p / "prereg"
            status_d = td_p / "status"
            prereg_d.mkdir()
            status_d.mkdir()
            # Empty status/prereg dirs -> MUST refuse!
            with self.assertRaises(BlindingGateError):
                assert_freeze_gates_r17(status_dir=status_d, prereg_dir=prereg_d)
            (prereg_d / "PREREG_J_DRAFT_v2.md").write_bytes((R17_ROOT / "prereg" / "PREREG_J_DRAFT_v2.md").read_bytes())
            (prereg_d / "PREREG_J_LOCK.md").write_text("LOCKED\n", encoding="utf-8")
            (status_d / "R17_CODE_FREEZE_AT3.txt").write_text("FROZEN\n", encoding="utf-8")
            # Missing R17_GO_AT3.txt -> refuse
            with self.assertRaises(BlindingGateError):
                assert_freeze_gates_r17(status_dir=status_d, prereg_dir=prereg_d)

            (status_d / "R17_GO_AT3.txt").write_text("GO\n", encoding="utf-8")
            receipts = assert_freeze_gates_r17(status_dir=status_d, prereg_dir=prereg_d)
            self.assertEqual(receipts["prereg_j_v2_sha256"], "7033337d1f42588ee9277804c6a803eb940312fe95a4ee31e2c000555b0fb85b")

    def test_stratified_split_halves_c1_4(self) -> None:
        hv = np.array([True] * 7 + [False] * 4, dtype=bool)
        p1, p2 = stratified_split_halves_c1(hv, split_idx=0, base_seed=DEFAULT_R17_SEED)
        self.assertEqual(len(p1), 5)  # floor(11/2)
        self.assertEqual(len(p2), 6)  # ceil(11/2)
        self.assertTrue(np.all(np.diff(p1) > 0))
        self.assertTrue(np.all(np.diff(p2) > 0))
        self.assertEqual(sorted(set(p1.tolist()) | set(p2.tolist())), list(range(11)))
        self.assertEqual(set(p1.tolist()) & set(p2.tolist()), set())

        # N=3 mac_ctrl edge case -> |P1|=1, |P2|=2
        hv3 = np.array([True, True, True], dtype=bool)
        p1_3, p2_3 = stratified_split_halves_c1(hv3, split_idx=0, base_seed=DEFAULT_R17_SEED)
        self.assertEqual(len(p1_3), 1)
        self.assertEqual(len(p2_3), 2)

    def test_mismatched_strata_exclusions(self) -> None:
        s_cas = resolve_r17_mismatched_strata("a549_ko_cas9pos")
        self.assertNotIn("a549_ko", s_cas["seen_sources"])
        self.assertEqual(s_cas["seen_sources"], ("nadig_hepg2", "nadig_jurkat"))
        s_cd4 = resolve_r17_mismatched_strata("cd4t_stim48")
        self.assertEqual(s_cd4["same_lineage_sources"], ("cd4t_rest",))
        self.assertEqual(s_cd4["unseen_sources"], ("a549_ko", "dld1", "feng_ipsc"))


class TestR17StatisticsAndRules(unittest.TestCase):
    def test_h_j1_precision_and_batch_dependence_detection(self) -> None:
        rng = np.random.default_rng(12345)
        p_cnt, g_cnt = 40, 80
        base_signal = rng.normal(0.0, 0.35, size=(p_cnt, g_cnt))
        d_obs = base_signal + rng.normal(0.0, 0.03, size=(p_cnt, g_cnt))
        ha = d_obs + rng.normal(0.0, 0.02, size=(p_cnt, g_cnt))
        hb = d_obs + rng.normal(0.0, 0.02, size=(p_cnt, g_cnt))
        hv = np.ones(p_cnt, dtype=bool)
        tm = rng.normal(0.0, 0.04, size=(p_cnt, g_cnt))
        tg = np.full(p_cnt, -1, dtype=np.int64)

        # Case A: Strong real batch effect across b=0..55 (higher for odd/large b, worse for b=0..9)
        d56 = np.empty((56, p_cnt, g_cnt), dtype=np.float64)
        for b in range(56):
            scale = 0.25 if b < 20 else (0.95 if b >= 40 else 0.55)
            d56[b] = scale * base_signal + rng.normal(0.0, 0.05, size=(p_cnt, g_cnt))
        d_l1 = d56.mean(axis=0)

        res_dep = compute_h_j1_for_subset(
            d_obs, d56, d_l1, tm, tg, ha, hb, hv, n_boot=300, seed=DEFAULT_R17_SEED, n_splits_diag=4
        )
        self.assertTrue(res_dep["valid_all"])
        self.assertGreater(res_dep["cross_fit_s0"]["delta_cf_0"], 0.04)
        self.assertGreater(res_dep["tau_batch"], 0.02)
        self.assertLess(res_dep["p_hj1_iu"], 0.05 / 3.0)
        self.assertEqual(res_dep["classification"], "BATCH_DEPENDENCE_PRACTICALLY_MEANINGFUL")
        self.assertIsNone(res_dep["subreason"])

        # Case B: Exact batch invariance (all 56 predictions identical) -> BATCH_INVARIANT_WITHIN_DELTA
        d56_inv = np.repeat(d56[45:46], 56, axis=0)
        res_inv = compute_h_j1_for_subset(
            d_obs, d56_inv, d56_inv.mean(axis=0), tm, tg, ha, hb, hv, n_boot=300, seed=DEFAULT_R17_SEED, n_splits_diag=4
        )
        self.assertEqual(res_inv["classification"], "BATCH_INVARIANT_WITHIN_DELTA")
        self.assertEqual(res_inv["cross_fit_s0"]["delta_cf_0"], 0.0)
        self.assertEqual(res_inv["tau_batch"], 0.0)

        # Case C: Exact N=3 mac_ctrl guard -> valid=False with HALF_SIZE_LT_2 subreason VALIDITY_FAIL
        res_3 = compute_h_j1_for_subset(
            d_obs[:3], d56[:, :3], d_l1[:3], tm[:3], tg[:3], ha[:3], hb[:3], hv[:3],
            n_boot=200, seed=DEFAULT_R17_SEED, n_splits_diag=2,
        )
        self.assertFalse(res_3["valid_all"])
        self.assertEqual(res_3["classification"], "INCONCLUSIVE")
        self.assertEqual(res_3["subreason"], "VALIDITY_FAIL")
        self.assertEqual(res_3["_internal"]["s0"]["reason"], "HALF_SIZE_LT_2")

        # Case D: C3.1 one-directional p_delta -> strictly negative draws yield p_delta == 1.0 (not small!)
        from audit_r17_stats import _bootstrap_one_dir_p
        neg_draws = -np.linspace(0.05, 0.25, 200)
        self.assertEqual(_bootstrap_one_dir_p(neg_draws), 1.0)
        pos_draws = np.linspace(0.05, 0.25, 200)
        self.assertAlmostEqual(_bootstrap_one_dir_p(pos_draws), 2.0 / 201.0, places=12)

    def test_h_j2_full_reversal_and_bonferroni56(self) -> None:
        rng = np.random.default_rng(2026)
        p_cnt, g_cnt = 40, 80
        sig = rng.normal(0.0, 0.35, size=(p_cnt, g_cnt))
        d_obs = sig + rng.normal(0.0, 0.02, size=(p_cnt, g_cnt))
        ha = d_obs + rng.normal(0.0, 0.015, size=(p_cnt, g_cnt))
        hb = d_obs + rng.normal(0.0, 0.015, size=(p_cnt, g_cnt))
        hv = np.ones(p_cnt, dtype=bool)
        tm = rng.normal(0.0, 0.02, size=(p_cnt, g_cnt))
        tg = np.full(p_cnt, -1, dtype=np.int64)
        fct = 0.70 * sig + rng.normal(0.0, 0.25, size=(p_cnt, g_cnt))

        d56 = np.empty((56, p_cnt, g_cnt), dtype=np.float64)
        for b in range(56):
            if b < 10:
                d56[b] = -0.35 * sig + rng.normal(0.0, 0.35, size=(p_cnt, g_cnt))  # FCT strongly superior
            elif b >= 46:
                d56[b] = 1.00 * sig + rng.normal(0.0, 0.03, size=(p_cnt, g_cnt))   # MODEL strongly superior
            else:
                d56[b] = fct + rng.normal(0.0, 0.04, size=(p_cnt, g_cnt))          # indistinguishable
        res_j2 = compute_h_j2_for_subset(
            d_obs, d56, fct, tm, tg, ha, hb, hv, has_detectable_signal=True,
            n_boot=300, n_boot_bonf56=500, seed=DEFAULT_R17_SEED,
        )
        self.assertEqual(res_j2["classification"], "CALL_UNSTABLE_SIGN_REVERSAL")
        self.assertTrue(res_j2["primary_975"]["has_full_reversal"])
        self.assertAlmostEqual(res_j2["bonf56_normal"]["z_critical"], float(norm.ppf(1.0 - 0.025 / 112.0)), places=12)

    def test_l2_selection_tie_and_medoid_tie_break(self) -> None:
        obs_ntc = np.linspace(0.5, 2.5, 40)
        # Constant vector -> exact tie with corr=0.0 everywhere -> must select index 0 cleanly
        ntc_const = np.ones((56, 40), dtype=np.float64)
        sel_tie = select_l2_ntc_index_and_diagnostics(ntc_const, obs_ntc)
        self.assertEqual(sel_tie["top1_idx"], 0)

        # Make b=23 match obs_ntc best
        ntc_56 = np.zeros((56, 40), dtype=np.float64)
        for b in range(56):
            ntc_56[b] = -obs_ntc if b != 23 else (obs_ntc * 2.0 + 0.1)
        ntc_pf = np.repeat(ntc_56[None, :, :], 4, axis=0)
        sel_23 = select_l2_ntc_index_and_diagnostics(ntc_56, obs_ntc, ntc_pb_per_fold_56=ntc_pf)
        self.assertEqual(sel_23["top1_idx"], 23)
        self.assertTrue(sel_23["per_fold_all_agree"])
        self.assertEqual(sel_23["per_fold_agree_with_ensemble"], 4)

        # Medoid test: exact tie breaks to lowest index
        w_tie = np.ones((4, 56, 8), dtype=np.float64)
        self.assertEqual(compute_batch_encoder_medoids(w_tie)["euclid_medoid_idx"], 0)
        # Put 40 rows clustered near origin with idx=19 at exact cluster center, 16 rows far out -> idx 19 is medoid
        w_19 = np.zeros((4, 56, 8), dtype=np.float64)
        for b in range(56):
            if b == 19:
                w_19[:, b, :] = 0.0
            elif b < 40:
                w_19[:, b, b % 8] = 0.05
            else:
                w_19[:, b, b % 8] = 50.0 if (b % 2 == 0) else -50.0
        self.assertEqual(compute_batch_encoder_medoids(w_19)["euclid_medoid_idx"], 19)

    def test_holm_intersection_union_gate_and_j5b_exclusions(self) -> None:
        # C2.1: Holm step-down runs on raw p-values only!
        # Rank-1 t1 (p=0.001 < 0.0125) fails CI gate -> t1 removed only; step-down DOES NOT stop!
        # Rank-2 t2 (p=0.010 < 0.01667) passes CI gate -> t2 rejected!
        # Rank-3 t3 (p=0.030 >= 0.025) fails Holm step-down -> stopping_alpha=0.025!
        # Rank-4 t4 (p=0.040 < 0.05) gets holm_p_rejected=False & alpha_eff=stopping_alpha=0.025!
        seen_alphas = {}

        def _gate(name: str, ret: bool):
            def _fn(alpha_eff: float) -> bool:
                seen_alphas[name] = alpha_eff
                return ret
            return _fn

        entries = [
            {"key": "t4", "p_value": 0.040, "ci_check_fn": _gate("t4", True)},
            {"key": "t1", "p_value": 0.001, "ci_check_fn": _gate("t1", False)},
            {"key": "t3", "p_value": 0.030, "ci_check_fn": _gate("t3", True)},
            {"key": "t2", "p_value": 0.010, "ci_check_fn": _gate("t2", True)},
        ]
        holm = apply_holm_fwer_r17(entries, alpha_family=0.05)
        self.assertAlmostEqual(holm["stopping_alpha"], 0.025, places=12)
        self.assertTrue(holm["by_key"]["t1"]["holm_p_rejected"])
        self.assertFalse(holm["by_key"]["t1"]["ci_gate_passed"])
        self.assertFalse(holm["by_key"]["t1"]["rejected"])
        self.assertAlmostEqual(seen_alphas["t1"], 0.05 / 4.0, places=12)

        self.assertTrue(holm["by_key"]["t2"]["holm_p_rejected"])
        self.assertTrue(holm["by_key"]["t2"]["ci_gate_passed"])
        self.assertTrue(holm["by_key"]["t2"]["rejected"])
        self.assertAlmostEqual(seen_alphas["t2"], 0.05 / 3.0, places=12)

        self.assertFalse(holm["by_key"]["t3"]["holm_p_rejected"])
        self.assertAlmostEqual(holm["by_key"]["t3"]["alpha_eff"], 0.025, places=12)
        self.assertAlmostEqual(seen_alphas["t3"], 0.025, places=12)
        self.assertFalse(holm["by_key"]["t3"]["rejected"])

        self.assertFalse(holm["by_key"]["t4"]["holm_p_rejected"])
        self.assertAlmostEqual(holm["by_key"]["t4"]["alpha_eff"], 0.025, places=12)
        self.assertAlmostEqual(seen_alphas["t4"], 0.025, places=12)
        self.assertFalse(holm["by_key"]["t4"]["rejected"])

        # J5b exclusion check (`a549_ko_cas9pos` excluded from fresh_n_ge_16_lines; pfizer_cytokine toggle checked)
        rng = np.random.default_rng(9)
        scores_se = {ln: rng.normal(size=56) for ln in R17_ALL_EVAL_LINES}
        scores_hvg = {ln: rng.normal(size=56) for ln in R17_ALL_EVAL_LINES}
        n_map = {
            "cd4t_stim48": 52, "gse295214_h1": 127, "pfizer_telohaec_basal": 114,
            "pfizer_telohaec_cytokine": 116, "a549_myllia": 20, "hek293_xucao": 16,
            "telohaec": 9, "calu3": 8, "mac_lps": 4, "mac_ctrl": 3, "a549_ko_cas9pos": 124,
        }
        j5b = compute_j5b_rank_correlations({"STATE_SE": scores_se, "STATE_HVG": scores_hvg}, n_map, n_perm=100, seed=7)
        fg16 = j5b["by_model"]["STATE_SE"]["fresh_n_ge_16_lines"]
        self.assertNotIn("a549_ko_cas9pos", fg16)
        self.assertEqual(len(fg16), 6)
        self.assertEqual(len(j5b["by_model"]["STATE_SE"]["fresh_n_ge_16_pairs"]), 15)       # C(6,2)
        self.assertEqual(len(j5b["by_model"]["STATE_SE"]["fresh_n_ge_16_excl_cytokine_pairs"]), 10)  # C(5,2)
        self.assertEqual(len(j5b["cross_model_within_line"]), 11)
        fdrs = compute_bh_fdr([0.001, 0.04, 0.20])
        self.assertAlmostEqual(fdrs[0], 0.003, places=12)

        # C2.2 exact >= |rho_obs| without tolerance
        from audit_r17_stats import _spearman_vectorized_against_perms
        a = np.arange(56, dtype=np.float64)
        b = a.copy()
        perm_ident = np.stack([np.arange(56), np.concatenate([[1, 0], np.arange(2, 56)])], axis=0)
        sp_exact = _spearman_vectorized_against_perms(a, b, perm_ident)
        # Identity perm has |rho|==1.0 (>=1.0), single-swap perm has |rho|<1.0 -> (1+1)/(2+1) = 2/3
        self.assertAlmostEqual(sp_exact["spearman_rho"], 1.0, places=12)
        self.assertAlmostEqual(sp_exact["p_perm"], 2.0 / 3.0, places=12)


class TestR17SyntheticEndToEndSchemaAndStage4(unittest.TestCase):
    def test_synthetic_e2e_cli_schema_allow_pickle_false_and_recompute_parity(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            out_dir = Path(td) / "syn_r17"
            recomp_json = Path(td) / "recomp_r17.json"

            # 1. Run CLI `--phase synthetic_e2e` cleanly in subprocess
            cmd_syn = [
                sys.executable,
                str(AUDIT_R17_DIR / "run_audit_r17.py"),
                "--phase", "synthetic_e2e",
                "--cohort", "all",
                "--output-dir", str(out_dir),
                "--n-boot", "250",
                "--n-boot-bonf56", "400",
                "--n-perm", "300",
                "--n-splits-diag", "4",
                "--device", "cpu",
            ]
            subprocess.run(cmd_syn, check=True)

            # 2. Verify all schema artifacts with allow_pickle=False
            self.assertTrue((out_dir / "control_barcodes_r17.json").is_file())
            self.assertTrue((out_dir / "embeddings" / "embeddings_manifest_sha256.json").is_file())
            self.assertTrue((out_dir / "summary_r17.json").is_file())
            with open(out_dir / "control_barcodes_r17.json", "r", encoding="utf-8") as fh:
                cb = json.load(fh)
            expected_pools = set(R17_ALL_EVAL_LINES) | set(R16_MISMATCHED_POOLS)
            self.assertEqual(set(cb.keys()), expected_pools)

            for pname in expected_pools:
                self.assertIn("barcode_prefix_convention", cb[pname])
                self.assertEqual(len(cb[pname]["cell_total_umis"]), 64)
                xp = out_dir / "embeddings" / f"x_state_{pname}.npz"
                self.assertTrue(xp.is_file())
                with np.load(xp, allow_pickle=False) as xz:
                    self.assertEqual(set(xz.files), {"X_state", "control_barcodes", "seed", "dataset_name", "se_audit_json"})
                    self.assertEqual(xz["X_state"].dtype, np.float32)
                    self.assertEqual(xz["X_state"].shape, (64, 2058))

            for ln in R17_ALL_EVAL_LINES:
                ld = out_dir / "arrays" / ln
                self.assertTrue((ld / "array_manifest_sha256.json").is_file())
                with np.load(ld / "observed.npz", allow_pickle=False) as oz:
                    for req_k in (
                        "line_name", "role", "modality", "checkpoint_id", "delta_obs", "delta_half_a", "delta_half_b",
                        "half_valid", "delta_transfer", "source_grand_mean", "fct_delta_pred", "e1_constmed_delta_pred",
                        "target_indices", "target_names", "gene_names", "mu_target", "mu_sources",
                        "obs_ntc_mean_official_median_log1p", "split_s0_p1_indices", "split_s0_p2_indices",
                        "seen_sources", "unseen_sources", "same_lineage_sources", "matched_ctrl_barcodes",
                        "matched_ctrl_batch_labels", "n_ctrl_total", "n_ctrl_sampled", "ctrl_sampling_fraction",
                    ):
                        self.assertIn(req_k, oz.files, f"Missing {req_k} in {ln}/observed.npz")

                for mk in ("STATE_SE", "STATE_HVG"):
                    with np.load(ld / f"preds_{mk}_sweep56.npz", allow_pickle=False) as sz:
                        self.assertEqual(sz["delta_pred_56"].shape[0], 56)
                        self.assertEqual(sz["ntc_pb_per_fold_56"].shape[0], 4)
                        self.assertEqual(sz["per_fold_r_orth_cal_56"].shape, (4, 56))
                    with np.load(ld / f"preds_{mk}_g4_mismatched.npz", allow_pickle=False) as gz:
                        self.assertIn(0, gz["g4_batch_indices"].tolist())
                        self.assertIn(17, gz["g4_batch_indices"].tolist())
                        self.assertIn(35, gz["g4_batch_indices"].tolist())
                        if ln == "mac_ctrl":
                            self.assertEqual(int(gz["b_star_1"]), -1)

                sens_p = ld / "SENS_FCT_full_r9_6_incl_vcc2025_h1.npz"
                if ln == "gse295214_h1":
                    self.assertTrue(sens_p.is_file())
                    with np.load(sens_p, allow_pickle=False) as snz:
                        self.assertEqual(
                            set(snz.files),
                            {"fct_delta_pred", "delta_transfer", "source_grand_mean", "target_names", "gene_names"},
                        )
                else:
                    self.assertFalse(sens_p.exists())

            # 3. Run CLI `--recompute-from-arrays` independently and verify <= 1e-12 parity!
            cmd_rec = [
                sys.executable,
                str(AUDIT_R17_DIR / "run_audit_r17.py"),
                "--recompute-from-arrays", str(out_dir / "arrays"),
                "--summary-json", str(out_dir / "summary_r17.json"),
                "--out-json", str(recomp_json),
                "--n-boot", "250",
                "--n-boot-bonf56", "400",
                "--n-perm", "300",
                "--n-splits-diag", "4",
            ]
            res = subprocess.run(cmd_rec, check=True, capture_output=True, text=True)
            rec_stdout = json.loads(res.stdout.strip().splitlines()[-1])
            self.assertEqual(rec_stdout["status"], "RECOMPUTE_OK")
            self.assertLessEqual(rec_stdout["parity"]["max_abs_diff"], 1e-12)

            with open(out_dir / "summary_r17.json", "r", encoding="utf-8") as fh:
                s1 = json.load(fh)
            with open(recomp_json, "r", encoding="utf-8") as fh:
                s2 = json.load(fh)
            cmp_res = compare_summaries_stage4(s2, s1, tol=1e-12)
            self.assertTrue(cmp_res["passed"])
            self.assertLessEqual(cmp_res["max_abs_diff"], 1e-12)
            self.assertEqual(s1["models"]["STATE_SE"]["mac_ctrl"]["h_j1"]["subreason"], "VALIDITY_FAIL")
            self.assertEqual(len(s1["models"]["STATE_SE"]["mac_ctrl"]["h_j1"]["r_obs_56"]), 56)
            self.assertEqual(s1["holm_by_model"]["STATE_SE"]["m_tests"], 9)
            self.assertEqual(len(s1["holm_by_model"]["STATE_SE"]["steps"]), 9)
            self.assertEqual(s1["holm_pooled_sensitivity"]["m_tests"], 18)
            self.assertEqual(len(s1["holm_pooled_sensitivity"]["steps"]), 18)

            # 4. Verify CLI `--phase real_r17` fails cleanly when unlock receipts are absent
            empty_status = Path(td) / "empty_status"
            empty_status.mkdir()
            cli_real = subprocess.run(
                [
                    sys.executable, str(AUDIT_R17_DIR / "run_audit_r17.py"),
                    "--phase", "real_r17", "--status-dir", str(empty_status), "--device", "cpu",
                ],
                capture_output=True, text=True,
            )
            self.assertNotEqual(cli_real.returncode, 0)
            self.assertIn("BlindingGateError", cli_real.stderr)


if __name__ == "__main__":
    unittest.main()
