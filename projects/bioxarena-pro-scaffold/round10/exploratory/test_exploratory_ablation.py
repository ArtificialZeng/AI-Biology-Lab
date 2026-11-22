#!/usr/bin/env python3
"""Unit tests for R10-E exploratory baseline ablation script (`round10/exploratory/`)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

EXP_DIR = Path(__file__).resolve().parent
R10_SRC_DIR = EXP_DIR.parent / "src"
for p in (EXP_DIR, R10_SRC_DIR):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from fct_pipeline import FctSourceState, SourceLineInputR10, fct_from_fc_bar
from metrics_r10 import compute_residual_suite_r10
from run_baseline_ablation import build_ablation_predictors, compute_fct_perm100_for_mode


def test_ablation_identities_and_limits() -> None:
    rng = np.random.default_rng(77)
    N, G = 35, 80
    q_names = [f"G{i}" for i in range(N)]
    q_gidx = np.arange(N, dtype=np.int64)
    log_fc = rng.normal(0.0, 0.3, size=(N, G))

    srcs = []
    for s in range(6):
        mu_s = rng.uniform(0.2, 5.0, size=G)
        deltas = np.log1p(mu_s[None, :] * np.exp(log_fc)) - np.log1p(mu_s[None, :])
        srcs.append(SourceLineInputR10(f"S{s}", q_names, q_gidx, deltas, None, mu_s))
    state = FctSourceState.fit(srcs, mask_self=True)

    mu_t = rng.uniform(0.3, 6.0, size=G)
    preds, fc_bundle, _ = build_ablation_predictors(state, mu_t, q_names)
    t_mean, fc_bar, has_hits = fc_bundle

    # E5 keeps exact sorted marginal distribution of mu_t and matches rank order of mu_e3
    mu_e3 = np.mean([s.ctrl_mean_libnorm for s in srcs], axis=0)
    mu_e5 = np.sort(mu_t)[np.argsort(np.argsort(mu_e3, kind="mergesort"), kind="mergesort")]
    assert np.allclose(np.sort(mu_e5), np.sort(mu_t))
    assert np.allclose(preds["E5_RankMatchedMu"], fct_from_fc_bar(mu_e5, fc_bar, t_mean, has_hits))

    # Asymptotic limit check under r_orth_cal: large c -> E6 log(fc_bar), small c -> E7 (fc_bar - 1)
    obs = preds["E6_LogFCBar"] + rng.normal(0.0, 0.05, size=(N, G))
    ha = obs + rng.normal(0.0, 0.01, size=(N, G))
    hb = obs + rng.normal(0.0, 0.01, size=(N, G))
    c_huge = fct_from_fc_bar(np.full(G, 1e9), fc_bar, t_mean, has_hits)
    c_tiny = fct_from_fc_bar(np.full(G, 1e-8), fc_bar, t_mean, has_hits)

    s_huge = compute_residual_suite_r10(obs, c_huge, t_mean, q_gidx, ha, hb, n_boot=50, seed=1)
    s_e6 = compute_residual_suite_r10(obs, preds["E6_LogFCBar"], t_mean, q_gidx, ha, hb, n_boot=50, seed=1)
    s_tiny = compute_residual_suite_r10(obs, c_tiny, t_mean, q_gidx, ha, hb, n_boot=50, seed=1)
    s_e7 = compute_residual_suite_r10(obs, preds["E7_FCBarMinus1"], t_mean, q_gidx, ha, hb, n_boot=50, seed=1)

    assert abs(s_huge["resp_top50"]["r_orth_cal"]["estimate"] - s_e6["resp_top50"]["r_orth_cal"]["estimate"]) < 1e-4
    assert abs(s_tiny["resp_top50"]["r_orth_cal"]["estimate"] - s_e7["resp_top50"]["r_orth_cal"]["estimate"]) < 1e-4

    # Check perm100 median helper under all_nontarget vs resp_top50
    p_top = compute_fct_perm100_for_mode(
        state, mu_t, q_names, obs, t_mean, q_gidx, 0.9, 0.9, fc_bundle, mode="resp_top50", n_perm=10
    )
    p_all = compute_fct_perm100_for_mode(
        state, mu_t, q_names, obs, t_mean, q_gidx, 0.9, 0.9, fc_bundle, mode="all_nontarget", n_perm=10
    )
    assert np.isfinite(p_top["r_orth_cal"]["median"]) and np.isfinite(p_all["r_orth_cal"]["median"])


if __name__ == "__main__":
    test_ablation_identities_and_limits()
    print("R10E_EXPLORATORY_UNIT_TESTS_PASSED")
