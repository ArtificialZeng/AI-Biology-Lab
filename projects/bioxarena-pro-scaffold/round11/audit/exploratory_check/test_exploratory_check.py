#!/usr/bin/env python3
"""
Unit & adversarial tests for the Round 11 Exploratory Review Audit (`exploratory_check/`).
Run strictly via `python3 -I round11/audit/exploratory_check/test_exploratory_check.py -v`.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

EXPLORATORY_CHECK_DIR = Path(__file__).resolve().parent
AUDIT_DIR = EXPLORATORY_CHECK_DIR.parent
for p in (str(EXPLORATORY_CHECK_DIR), str(AUDIT_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from exploratory_core import evaluate_cross_mask_split_minus_all  # noqa: E402
from independent_r11_predictors import fct_from_fc_bar  # noqa: E402
from run_n5_audit import (  # noqa: E402
    _eval_n5_from_moments,
    _pred_n5_moments,
    build_n5_basis,
    prepare_n5_per_pert_data,
)
from run_noiseless_sim import _eval_noiseless_line  # noqa: E402
from run_q5a_blocked import (  # noqa: E402
    _directed_block_tests,
    _dyadic_cluster_test,
    _exact_line_sign_flip_p,
    max_independent_matching_size,
)


class TestExploratoryAuditSuite(unittest.TestCase):
    def test_1_n5_zero_by_construction_vs_orthogonal_signal(self) -> None:
        rng = np.random.default_rng(7)
        n_p, n_g = 25, 400
        t_mean = rng.normal(0.0, 0.25, size=(n_p, n_g))
        grand_mean = rng.normal(0.05, 0.15, size=n_g)
        fc_bar = np.exp(rng.normal(0.0, 0.35, size=(n_p, n_g)))

        # Construct genuine biological residual strictly orthogonal to N_5(p)
        target_cols = np.arange(n_p, dtype=np.int64)
        true_orth_sig = np.zeros((n_p, n_g), dtype=np.float64)
        for p in range(n_p):
            q = build_n5_basis(t_mean[p], grand_mean, fc_bar[p])
            self.assertLessEqual(q.shape[1], 6)
            raw_v = rng.normal(0.0, 0.2, size=n_g)
            orth_v = raw_v - q @ (q.T @ raw_v)
            true_orth_sig[p] = orth_v

        # Pure artifact + genuine orthogonal signal in noisy target observations
        noise_a = rng.normal(0.0, 0.05, size=(n_p, n_g))
        noise_b = rng.normal(0.0, 0.05, size=(n_p, n_g))
        obs_sig = t_mean + 0.5 * grand_mean[None, :] + 0.4 * (fct_from_fc_bar(fc_bar, np.full(n_g, 0.7)) - t_mean) + true_orth_sig
        ha = obs_sig + noise_a
        hb = obs_sig + noise_b
        delta_obs = 0.5 * (ha + hb)

        per_p, xx_o, haa_o, hbb_o, hab_o, meta = prepare_n5_per_pert_data(
            delta_obs=delta_obs,
            t_mean=t_mean,
            grand_mean=grand_mean,
            fc_bar=fc_bar,
            target_cols=target_cols,
            half_a_obs=ha,
            half_b_obs=hb,
            half_valid=np.ones(n_p, dtype=bool),
        )
        self.assertEqual(meta["n_eval_perts"], n_p)
        draws = rng.integers(0, n_p, size=(200, n_p))

        # GlobalScale & SGM must have eta == 0 and r_cone_sp == 0 by exact construction
        for name, pred in [
            ("GS", 1.35 * t_mean),
            ("SGM", np.broadcast_to(grand_mean[None, :], (n_p, n_g)).copy()),
        ]:
            xy, yy_o, yy_c = _pred_n5_moments(per_p, delta_obs, pred, t_mean)
            res = _eval_n5_from_moments(xy, yy_o, yy_c, xx_o, haa_o, hbb_o, hab_o, draws, 200)
            self.assertLess(res["eta_est"], 1e-12, msg=name)
            self.assertAlmostEqual(res["r_cone_sp_est"], 0.0, places=12, msg=name)
            self.assertAlmostEqual(res["r2_in_n5"], 1.0, places=12, msg=name)

        # Constant-baseline E1(c=0.7) is suppressed inside N_5 cone despite huge E1 drift in obs_sig!
        e1_pred = fct_from_fc_bar(fc_bar, np.full(n_g, 0.7))
        xy, yy_o, yy_c = _pred_n5_moments(per_p, delta_obs, e1_pred, t_mean)
        res_e1 = _eval_n5_from_moments(xy, yy_o, yy_c, xx_o, haa_o, hbb_o, hab_o, draws, 200)
        self.assertLess(res_e1["eta_est"], 0.02)
        self.assertGreater(res_e1["r2_in_n5"], 0.999)
        self.assertLess(abs(res_e1["r_cone_sp_est"]), 0.005)

        # Genuine orthogonal predictor preserves eta > 0.98 and recovers r_cone_sp ~ 1.0
        true_pred = t_mean + true_orth_sig
        xy, yy_o, yy_c = _pred_n5_moments(per_p, delta_obs, true_pred, t_mean)
        res_true = _eval_n5_from_moments(xy, yy_o, yy_c, xx_o, haa_o, hbb_o, hab_o, draws, 200)
        self.assertGreater(res_true["eta_est"], 0.98)
        self.assertGreater(res_true["r_cone_sp_est"], 0.90)

    def test_2_cross_mask_split_minus_all(self) -> None:
        rng = np.random.default_rng(19)
        n_p, n_g = 30, 300
        t_mean = rng.normal(0.0, 0.1, size=(n_p, n_g))
        fc_bar = np.ones((n_p, n_g), dtype=np.float64)
        # First 40 genes have coherent strong signal, remaining 260 are background noise
        sig = np.zeros((n_p, n_g), dtype=np.float64)
        sig[:, :40] = rng.normal(0.0, 0.8, size=(n_p, 40))
        fc_bar[:, :40] = np.exp(sig[:, :40])
        ha = t_mean + fct_from_fc_bar(fc_bar, np.full(n_g, 0.5)) - t_mean + rng.normal(0.0, 0.05, size=(n_p, n_g))
        hb = t_mean + fct_from_fc_bar(fc_bar, np.full(n_g, 0.5)) - t_mean + rng.normal(0.0, 0.05, size=(n_p, n_g))
        delta_obs = 0.5 * (ha + hb)
        pred = fct_from_fc_bar(fc_bar, np.full(n_g, 0.5))
        pred[:, 40:] = t_mean[:, 40:] - 0.5 * sig[:, :1].repeat(260, axis=1) * 0.0  # pure zero-signal mismatch on rest
        # Add opposite tail artifact on background so split top50 > all
        pred[:, :40] = t_mean[:, :40] + 0.8 * (delta_obs[:, :40] - t_mean[:, :40])
        pred[:, 40:] = t_mean[:, 40:] - 0.3 * (delta_obs[:, 40:] - t_mean[:, 40:])

        out = evaluate_cross_mask_split_minus_all(
            delta_obs=delta_obs,
            pred=pred,
            t_mean=t_mean,
            target_cols=np.arange(n_p),
            half_a_obs=ha,
            half_b_obs=hb,
            half_valid=np.ones(n_p, dtype=bool),
            top_k=50,
            n_boot=300,
            seed=20261124,
        )
        self.assertEqual(out["n_eval_perts"], n_p)
        self.assertGreater(out["split_est"], out["all_est"])
        self.assertGreater(out["diff_ci975_lo"], 0.0)

    def test_3_dyadic_cluster_and_blocked_perm(self) -> None:
        # Exact 2^4 = 16 sign-flip verification on 4 strictly positive blocks -> p_one = 1/16, p_two = 2/16
        p2, p1 = _exact_line_sign_flip_p(np.array([0.2, 0.4, 0.1, 0.3]))
        self.assertAlmostEqual(p1, 1.0 / 16.0, places=12)
        self.assertAlmostEqual(p2, 2.0 / 16.0, places=12)

        # Complete graph K_5 has max matching size 2
        nodes = ["A", "B", "C", "D", "E"]
        k5_pairs = [(nodes[i], nodes[j]) for i in range(5) for j in range(i + 1, 5)]
        self.assertEqual(max_independent_matching_size(k5_pairs, nodes), 2)

        # Simulated hub dependency: hub line H shifts all 4 incident edges by +0.5 while non-hub edges scatter around 0
        # Naive edge Wilcoxon can look too tiny, while line-level tests reflect N_lines=5!
        hub_diffs = np.array([0.5, 0.5, 0.5, 0.5, -0.05, 0.02, -0.01, 0.03, -0.02, 0.01])
        dy = _dyadic_cluster_test(k5_pairs, hub_diffs, nodes)
        self.assertLess(dy["line_exact_signflip_p_two"], 1.0)
        self.assertGreater(dy["dyadic_cgm_p_two"], dy["naive_wilcoxon_p_two"])

        dir_pairs = [(u, v) for (u, v) in k5_pairs] + [(v, u) for (u, v) in k5_pairs]
        dir_diffs = np.concatenate([hub_diffs, hub_diffs])
        db = _directed_block_tests(dir_pairs, dir_diffs, nodes)
        self.assertEqual(db["n_directed_pairs"], 20)
        self.assertEqual(db["n_blocks"], 5)

    def test_4_noiseless_theorem1_sign_behavior(self) -> None:
        rng = np.random.default_rng(99)
        n_p, n_g = 20, 500
        # Test both regimes:
        # (A) Moderately correlated downshifted mu_t: a_all < 0 (E1(all) < 0) while a_top50 > 0 (E1(top50) > 0)
        mu_s = np.exp(rng.normal(0.5, 1.2, size=n_g))
        mu_t_mid = 0.45 * mu_s * np.exp(rng.normal(0.0, 1.0, size=n_g))
        delta = rng.normal(0.0, 0.4, size=(n_p, n_g))
        fc_mat = np.exp(delta)
        mask_fin = np.ones((n_p, n_g), dtype=bool)
        eval_targets = [f"G{p}" for p in range(n_p)]
        sim = _eval_noiseless_line(
            fc_mat=fc_mat,
            mu_t=mu_t_mid,
            source_names=["S1"],
            source_ctrl_mu={"S1": mu_s},
            source_deltas={"S1": {t: np.zeros(n_g) for t in eval_targets}},
            eval_targets=eval_targets,
            target_cols=np.arange(n_p, dtype=np.int64),
            mask_finite=mask_fin,
        )
        # First-order Theorem 1 identity: simulated E1 sign matches analytical a_thm1_delta2 sign,
        # and tail conditioning strictly raises E1(top50) > E1(all).
        self.assertEqual(np.sign(sim["all_nontarget"]["e1_r_orth"]), np.sign(sim["all_nontarget"]["a_thm1_delta2"]))
        self.assertLess(sim["all_nontarget"]["e1_r_orth"], 0.0)
        self.assertGreater(sim["resp_top50_noiseless_split"]["e1_r_orth"], sim["all_nontarget"]["e1_r_orth"])

        # (B) Realistic uncorrelated baseline mismatch (toy-like): both all and top50 become positive
        mu_t_ind = 0.5 * np.exp(rng.normal(0.5, 1.2, size=n_g))
        sim_ind = _eval_noiseless_line(
            fc_mat=fc_mat,
            mu_t=mu_t_ind,
            source_names=["S1"],
            source_ctrl_mu={"S1": mu_s},
            source_deltas={"S1": {t: np.zeros(n_g) for t in eval_targets}},
            eval_targets=eval_targets,
            target_cols=np.arange(n_p, dtype=np.int64),
            mask_finite=mask_fin,
        )
        self.assertGreater(sim_ind["all_nontarget"]["e1_r_orth"], 0.0)
        self.assertGreater(sim_ind["resp_top50_noiseless_split"]["e1_r_orth"], sim_ind["all_nontarget"]["e1_r_orth"])


if __name__ == "__main__":
    unittest.main()
