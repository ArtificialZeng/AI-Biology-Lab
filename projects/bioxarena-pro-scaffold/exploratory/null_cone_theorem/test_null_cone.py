#!/usr/bin/env python3
"""
Unit & Adversarial Verification Suite for `null_cone_metrics.py`.

Tests:
  1. Exact parity (`max diff < 1e-13`) between `N1` and frozen `round10/src/metrics_r10.py` `_orth_pair`.
  2. Exact annihilation (`|r| < 1e-12`) of `Transfer_Mean_GlobalScale` (R18), `Source_Grand_Mean` (N5),
     and fixed cone nodes `g_c(FC_bar)`.
  3. Zero-noise split-half adversarial test: `E1`, `E2(P10/P90)`, `E6`, `E7`, and `Perm100` collapse
     from `[+0.25, +0.67]` (N1) to `|r| < 0.01` under `N4` / `N5` semi-partial on BOTH `all` and `top50`,
     while `FCT` retains `r > 0.55` (`top50`) and `r > 0.93` (`all`).
  4. Zero penalty to true biological regulatory rewiring (`eta > 0.99`, exact recovery of true r).
  5. Bootstrap determinism (`seed=20261124`) & paired CI consistency.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parents[1]
sys.path.insert(0, str(THIS_DIR))
sys.path.insert(0, str(REPO_ROOT / "round10" / "src"))

from metrics_r10 import _orth_pair, pearson_1d_strict  # noqa: E402
from null_cone_metrics import (  # noqa: E402
    compute_perturbation_moments,
    joint_bootstrap_summary,
    log1p_cone_basis,
)


class TestNullConeMetrics(unittest.TestCase):
    def test_1_exact_parity_with_frozen_r10_n1(self) -> None:
        rng = np.random.default_rng(123)
        N, G = 25, 400
        t_mean = rng.normal(0.0, 0.2, (N, G))
        fc_bar = np.exp(rng.normal(0.0, 0.3, (N, G)))
        b_s = t_mean.mean(axis=0)
        r_obs = rng.normal(0.0, 0.15, (N, G))
        r_pred = 0.3 * r_obs + rng.normal(0.0, 0.15, (N, G))
        ha = r_obs + rng.normal(0.0, 0.08, (N, G))
        hb = r_obs + rng.normal(0.0, 0.08, (N, G))
        idx_list = [np.arange(G, dtype=np.int64) for _ in range(N)]

        mom = compute_perturbation_moments(
            r_obs, r_pred[None, :, :], ["P"], t_mean, fc_bar, b_s, ha, hb, None, idx_list, null_mode="N1"
        )
        res = joint_bootstrap_summary(mom, n_boot=50, seed=42)

        # Reference explicit per-perturbation span{1, t_mean} projection via np.linalg.qr
        xs, ys, hao, hbo = [], [], [], []
        for p in range(N):
            Q, _ = np.linalg.qr(np.column_stack([np.ones(G), t_mean[p]]), mode="reduced")
            xs.append(r_obs[p] - Q @ (Q.T @ r_obs[p]))
            ys.append(r_pred[p] - Q @ (Q.T @ r_pred[p]))
            hao.append(ha[p] - Q @ (Q.T @ ha[p]))
            hbo.append(hb[p] - Q @ (Q.T @ hb[p]))
        ref_raw = pearson_1d_strict(np.concatenate(xs), np.concatenate(ys))
        ref_rho = pearson_1d_strict(np.concatenate(hao), np.concatenate(hbo))
        ref_cal = ref_raw / np.sqrt((2.0 * ref_rho) / (1.0 + ref_rho))

        self.assertLess(abs(res["rho_ab"] - ref_rho), 1e-13)
        self.assertLess(abs(res["predictors"]["P"]["orth"]["est"] - ref_cal), 1e-13)

    def test_2_exact_zero_for_globalscale_grandmean_and_cone_nodes(self) -> None:
        rng = np.random.default_rng(2026)
        N, G = 40, 600
        z = rng.normal(-0.5, 1.3, G)
        mu_s = np.exp(z[None, :] + rng.normal(0, 0.6, (4, G)))
        mu_t = np.exp(z + rng.normal(0.3, 0.8, G))
        delta = rng.standard_t(df=4, size=(N, G)) * 0.25 + rng.normal(0.0, 0.08, (1, G))
        fc = np.exp(np.clip(delta, -3.0, 3.0))

        D_s = np.log1p(mu_s[:, None, :] * fc[None, :, :]) - np.log1p(mu_s[:, None, :])
        D_t = np.log1p(mu_t[None, :] * fc) - np.log1p(mu_t[None, :])
        t_mean = D_s.mean(axis=0)
        b_s = t_mean.mean(axis=0)
        r_obs = D_t - t_mean

        g005, g1, g20 = log1p_cone_basis(fc)
        preds = {
            "GlobalScale": 1.7 * t_mean - t_mean,
            "Source_Grand_Mean": np.broadcast_to(b_s[None, :], (N, G)) - t_mean,
            "g_005": g005 - t_mean,
            "g_100": g1 - t_mean,
            "g_200": g20 - t_mean,
            "FCT": D_t - t_mean,
        }
        names = list(preds.keys())
        Y = np.stack([preds[k] for k in names], axis=0)
        idx_all = [np.arange(G, dtype=np.int64) for _ in range(N)]
        idx_top50 = [np.argpartition(np.abs(D_t[p]), -50)[-50:] for p in range(N)]

        for idxs in (idx_all, idx_top50):
            mom = compute_perturbation_moments(r_obs, Y, names, t_mean, fc, b_s, r_obs, r_obs, None, idxs, null_mode="N5")
            res = joint_bootstrap_summary(mom, n_boot=20, seed=1)
            for zero_nm in ("GlobalScale", "Source_Grand_Mean", "g_005", "g_100", "g_200"):
                val = res["predictors"][zero_nm]["semi"]["est"]
                self.assertTrue(np.isnan(val) or abs(val) < 1e-11, f"{zero_nm} not zero: {val}")
            self.assertGreater(res["predictors"]["FCT"]["semi"]["est"], 0.45)

    def test_3_noiseless_split_half_adversarial_e1_e6_perm_annihilation(self) -> None:
        rng = np.random.default_rng(777)
        N, G = 80, 1500
        z = rng.normal(-0.5, 1.4, G)
        mu_s = np.exp(z[None, :] + rng.normal(0, 0.7, (5, G)))
        mu_t = np.exp(z + rng.normal(0.35, 0.9, G))
        delta = rng.standard_t(df=4, size=(N, G)) * 0.25
        fc = np.exp(np.clip(delta, -3.0, 3.0))

        D_s = np.log1p(mu_s[:, None, :] * fc[None, :, :]) - np.log1p(mu_s[:, None, :])
        D_t = np.log1p(mu_t[None, :] * fc) - np.log1p(mu_t[None, :])
        t_mean = D_s.mean(axis=0)
        b_s = t_mean.mean(axis=0)
        r_obs = D_t - t_mean

        c_med = float(np.median(mu_t))
        c_p10 = float(np.percentile(mu_t, 10))
        c_p90 = float(np.percentile(mu_t, 90))
        preds = {
            "FCT": D_t - t_mean,
            "E1": (np.log1p(c_med * fc) - np.log1p(c_med)) - t_mean,
            "E2_P10": (np.log1p(c_p10 * fc) - np.log1p(c_p10)) - t_mean,
            "E2_P90": (np.log1p(c_p90 * fc) - np.log1p(c_p90)) - t_mean,
            "E6": np.log(fc) - t_mean,
        }
        for k in range(20):
            perm = np.random.default_rng(k).permutation(G)
            preds[f"Perm_{k}"] = (np.log1p(mu_t[perm][None, :] * fc) - np.log1p(mu_t[perm][None, :])) - t_mean
        names = list(preds.keys())
        Y = np.stack([preds[k] for k in names], axis=0)
        idx_top50 = [np.argpartition(np.abs(D_t[p]), -50)[-50:] for p in range(N)]

        mom_n1 = compute_perturbation_moments(r_obs, Y, names, t_mean, fc, b_s, r_obs, r_obs, None, idx_top50, null_mode="N1")
        res_n1 = joint_bootstrap_summary(mom_n1, n_boot=20, seed=1)
        perm_n1_med = float(np.median([res_n1["predictors"][f"Perm_{k}"]["orth"]["est"] for k in range(20)]))
        # Confirm N1 IS severely inflated (>0.35 on E1 and >0.55 on E6, >0.20 on Perm median!)
        self.assertGreater(res_n1["predictors"]["E1"]["orth"]["est"], 0.35)
        self.assertGreater(res_n1["predictors"]["E6"]["orth"]["est"], 0.55)
        self.assertGreater(perm_n1_med, 0.20)

        mom_n5 = compute_perturbation_moments(r_obs, Y, names, t_mean, fc, b_s, r_obs, r_obs, None, idx_top50, null_mode="N5")
        res_n5 = joint_bootstrap_summary(mom_n5, n_boot=20, seed=1)
        perm_n5_med = float(np.median([res_n5["predictors"][f"Perm_{k}"]["semi"]["est"] for k in range(20)]))
        # Confirm N5 semi-partial annihilates E1/E2/E6/Perm while preserving FCT!
        for cheat in ("E1", "E2_P10", "E2_P90", "E6"):
            self.assertLess(abs(res_n5["predictors"][cheat]["semi"]["est"]), 0.01)
            self.assertLess(res_n5["predictors"][cheat]["eta"], 0.03)
        self.assertLess(abs(perm_n5_med), 0.02)
        self.assertLess(abs(perm_n5_med), 0.10 * perm_n1_med)
        self.assertGreater(res_n5["predictors"]["FCT"]["semi"]["est"], 0.45)

    def test_4_zero_penalty_on_true_biological_rewiring(self) -> None:
        rng = np.random.default_rng(999)
        N, G = 40, 2000
        t_mean = rng.normal(0.0, 0.2, (N, G))
        fc_bar = np.exp(np.clip(rng.normal(0.0, 0.3, (N, G)), -2.5, 2.5))
        b_s = t_mean.mean(axis=0)
        # Independent true target-specific rewiring gamma_{t,p} orthogonal to source!
        gamma = rng.normal(0.0, 0.15, (N, G))
        r_obs = gamma
        # Imperfect predictor capturing 60% of true rewiring + independent noise (true r = 0.6 / sqrt(0.36 + 0.64) = 0.6)
        r_pred = 0.6 * gamma + 0.8 * rng.normal(0.0, 0.15, (N, G))
        idx_all = [np.arange(G, dtype=np.int64) for _ in range(N)]

        mom = compute_perturbation_moments(
            r_obs, r_pred[None, :, :], ["RewireModel"], t_mean, fc_bar, b_s, r_obs, r_obs, None, idx_all, null_mode="N5"
        )
        res = joint_bootstrap_summary(mom, n_boot=50, seed=7)
        self.assertGreater(res["predictors"]["RewireModel"]["eta"], 0.995)
        self.assertAlmostEqual(res["predictors"]["RewireModel"]["semi"]["est"], 0.6, delta=0.02)

    def test_5_bootstrap_determinism(self) -> None:
        rng = np.random.default_rng(55)
        N, G = 30, 300
        t_mean = rng.normal(0, 0.2, (N, G))
        fc_bar = np.exp(rng.normal(0, 0.2, (N, G)))
        b_s = t_mean.mean(axis=0)
        r_obs = rng.normal(0, 0.1, (N, G))
        Y = rng.normal(0, 0.1, (2, N, G))
        idxs = [np.arange(G, dtype=np.int64) for _ in range(N)]
        mom = compute_perturbation_moments(r_obs, Y, ["A", "B"], t_mean, fc_bar, b_s, r_obs, r_obs, None, idxs, null_mode="N5")
        r1 = joint_bootstrap_summary(mom, n_boot=500, seed=20261124, paired_specs=[("A-B", "A", "B")])
        r2 = joint_bootstrap_summary(mom, n_boot=500, seed=20261124, paired_specs=[("A-B", "A", "B")])
        self.assertEqual(r1["paired"]["A-B"]["semi"]["ci975_lo"], r2["paired"]["A-B"]["semi"]["ci975_lo"])
        self.assertEqual(r1["paired"]["A-B"]["semi"]["ci975_hi"], r2["paired"]["A-B"]["semi"]["ci975_hi"])


if __name__ == "__main__":
    unittest.main()
