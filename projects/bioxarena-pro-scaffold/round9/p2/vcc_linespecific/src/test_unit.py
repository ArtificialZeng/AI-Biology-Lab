#!/usr/bin/env python3
"""Unit tests for predictors and metrics using synthetic planted-truth systems."""

import numpy as np

from adjudicate import adjudicate_results, render_markdown_report
from data_loader import LineDataset, PREREG_COMMIT_FILE, check_prereg_gate, is_replogle_dev_line
from metrics import (
    cis_non_overlapping,
    compute_line_specific_metrics_and_boots,
    compute_row_metrics,
    drf_score,
    paired_row_bootstrap_ci,
    single_array_bootstrap_ci,
    spearman_brown_ceiling,
    stratified_row_bootstrap_ci,
    tie_aware_disc_rows,
)
from predictors import (
    SourceLineInput,
    SourceModelState,
    compute_control_covariance,
    predict_zero_shot_suite,
)
from run_eval import control_splithalf_covariance_reliability, evaluate_datasets


def _make_synthetic_system(seed: int = 0):
    rng = np.random.RandomState(seed)
    G = 80
    N_ctrl = 400
    N_pert = 30
    target_indices = np.arange(N_pert, dtype=np.int64)
    target_names = tuple(f"G{i}" for i in target_indices)

    # Shared perturbation signature across cell lines (low-rank stress + per-gene shared component)
    shared_gene_prog = rng.randn(N_pert, G) * 0.25
    shared_base = rng.randn(G) * 0.20

    # Source lines S1, S2 and Target line T with distinct control factor loadings
    def gen_ctrl_and_deltas(line_seed: int, true_lam: float):
        lrng = np.random.RandomState(line_seed)
        # Rank-10 latent gene regulatory covariance
        B = lrng.randn(G, 10) * 0.35
        Z = lrng.randn(N_ctrl, 10)
        noise = lrng.randn(N_ctrl, G) * 0.15
        X_ctrl = Z @ B.T + noise
        C_true = compute_control_covariance(X_ctrl)
        # True line-specific effect follows linear response on C_true
        deltas = np.empty((N_pert, G), dtype=np.float64)
        for i, g in enumerate(target_indices):
            deltas[i] = shared_base + shared_gene_prog[i] + true_lam * (C_true[:, g] / C_true[g, g])
        return X_ctrl, C_true, deltas

    X_s1, C_s1, D_s1 = gen_ctrl_and_deltas(101, true_lam=-1.5)
    X_s2, C_s2, D_s2 = gen_ctrl_and_deltas(202, true_lam=-1.5)
    X_t, C_t, D_t = gen_ctrl_and_deltas(303, true_lam=-1.5)
    return (
        [
            SourceLineInput("S1", C_s1, D_s1, target_indices, target_names),
            SourceLineInput("S2", C_s2, D_s2, target_indices, target_names),
        ],
        C_t,
        D_t,
        target_indices,
        target_names,
    )


def test_cipher_and_hybrid_planted_recovery():
    src_lines, C_t, D_t, target_indices, target_names = _make_synthetic_system(seed=7)
    # Also construct a pure linear-response system (zero shared offset) to test exact lambda fit
    G = C_t.shape[0]
    true_lam = -1.75
    pure_srcs = []
    for s in src_lines:
        D_pure = np.array([true_lam * s.ctrl_cov[:, g] / s.ctrl_cov[g, g] for g in s.target_indices])
        pure_srcs.append(SourceLineInput(s.name, s.ctrl_cov, D_pure, s.target_indices, s.target_names))
    state_pure = SourceModelState.fit(pure_srcs, mask_self=True)
    assert abs(state_pure.lam - true_lam) < 1e-10, f"Expected lambda={true_lam}, got {state_pure.lam}"

    # Full system with shared + line-specific linear response
    state = SourceModelState.fit(src_lines, mask_self=True)
    preds = predict_zero_shot_suite(
        state=state,
        target_ctrl_cov=C_t,
        query_target_names=target_names,
        query_target_gene_indices=target_indices,
        mismatched_ctrl_cov=src_lines[0].ctrl_cov,
        shuffle_seed=42,
    )

    m_transfer = compute_row_metrics(preds["Transfer_Mean"], D_t, target_indices, centroid=state.grand_mean)
    m_hybrid = compute_row_metrics(preds["Hybrid_Raw"], D_t, target_indices, centroid=state.grand_mean)
    m_mis = compute_row_metrics(
        preds["NegCtrl_Hybrid_MismatchedCov"], D_t, target_indices, centroid=state.grand_mean
    )
    m_shuf = compute_row_metrics(
        preds["NegCtrl_Hybrid_ShuffledGenes"], D_t, target_indices, centroid=state.grand_mean
    )

    # Hybrid must beat Transfer_Mean and negative controls on planted linear-response line-specific signal
    assert m_hybrid["mse"].mean() < m_transfer["mse"].mean()
    assert m_hybrid["delta_r"].mean() > m_transfer["delta_r"].mean()
    assert m_hybrid["delta_r"].mean() > m_mis["delta_r"].mean()
    assert m_hybrid["delta_r"].mean() > m_shuf["delta_r"].mean()
    assert m_hybrid["disc_raw"].mean() > 0.9


def test_metrics_and_bootstrap():
    # Constant predictions across rows must average to 0.5; zero-effect predictions score 0.5 per row
    rng = np.random.RandomState(11)
    T = rng.randn(20, 50)
    P_const = np.tile(rng.randn(1, 50), (20, 1))
    disc = tie_aware_disc_rows(P_const, T)
    assert abs(float(disc.mean()) - 0.5) < 1e-12, f"Constant pred mean should be 0.5, got {disc.mean()}"
    disc_zero = tie_aware_disc_rows(np.zeros_like(T), T)
    assert np.allclose(disc_zero, 0.5), f"Zero pred should give 0.5 per row, got {disc_zero}"

    # DRF check
    assert abs(drf_score(0.8, 0.5) - 0.6) < 1e-12
    assert abs(drf_score(0.1, 0.5, lower_is_better=True) - 0.8) < 1e-12

    # Paired bootstrap check
    a = np.linspace(0.5, 1.0, 40)
    b = np.linspace(0.2, 0.7, 40)
    res = paired_row_bootstrap_ci(a, b, n_boot=1000, seed=12345)
    assert abs(res["delta"] - 0.3) < 1e-10
    assert res["ci_gt_zero"] is True

    s_res = single_array_bootstrap_ci(a, n_boot=1000, seed=12345)
    assert s_res["ci_gt_zero"] is True

    strat = stratified_row_bootstrap_ci(a, b, ["L1"] * 20 + ["L2"] * 20, n_boot=1000, seed=12345)
    assert strat["ci_gt_zero"] is True

    assert cis_non_overlapping((0.0, 0.1), (0.2, 0.3)) is True
    assert cis_non_overlapping((0.0, 0.25), (0.2, 0.3)) is False


def test_prereg_gate_and_splithalf():
    assert is_replogle_dev_line("replogle_k562") is True
    assert is_replogle_dev_line("replogle_rpe1") is True
    assert is_replogle_dev_line("nadig_hepg2") is False

    # Replogle dev lines must pass gate even without PREREG_COMMIT.txt
    check_prereg_gate(["replogle_k562", "replogle_rpe1"])
    if not PREREG_COMMIT_FILE.exists():
        blocked = False
        try:
            check_prereg_gate(["nadig_hepg2"])
        except RuntimeError:
            blocked = True
        assert blocked, "Non-Replogle line must be blocked before PREREG_COMMIT.txt exists!"

    rng = np.random.RandomState(5)
    Z = rng.randn(500, 10)
    B = rng.randn(40, 10)
    X_ctrl = Z @ B.T + rng.randn(500, 40) * 0.1
    rel = control_splithalf_covariance_reliability(X_ctrl, np.arange(20), seed=42)
    assert rel["ctrl_cov_splithalf_r_mean"] > 0.8


def test_end_to_end_synthetic_loo():
    rng = np.random.RandomState(99)
    G = 70
    N_ctrl = 450
    N_pert = 25
    genes = [f"G{i}" for i in range(G)]
    g2i = {g: i for i, g in enumerate(genes)}
    t_names = tuple(genes[:N_pert])
    t_idx = np.arange(N_pert, dtype=np.int64)
    shared_prog = rng.randn(N_pert, G) * 0.25

    datasets = {}
    for s_idx, lname in enumerate(["Replogle_K562", "Replogle_RPE1", "Synthetic_L3"]):
        lrng = np.random.RandomState(500 + s_idx)
        B = lrng.randn(G, 8) * 0.4
        X_ctrl = lrng.randn(N_ctrl, 8) @ B.T + lrng.randn(N_ctrl, G) * 0.1
        C = compute_control_covariance(X_ctrl)
        deltas = np.array([shared_prog[i] - 1.4 * C[:, g] / C[g, g] for i, g in enumerate(t_idx)])
        ha = deltas + lrng.randn(*deltas.shape) * 0.05
        hb = deltas + lrng.randn(*deltas.shape) * 0.05
        datasets[lname] = LineDataset(
            name=lname,
            genes=genes,
            gene_to_idx=g2i,
            ctrl_cells=X_ctrl,
            deltas=deltas,
            target_names=t_names,
            target_indices=t_idx,
            repr_mask=np.ones(N_pert, dtype=bool),
            representation="log1p",
            half_a_deltas=ha,
            half_b_deltas=hb,
            half_valid=np.ones(N_pert, dtype=bool),
        )

    df_lines, df_boot, summary = evaluate_datasets(
        datasets=datasets,
        n_boot=200,
        seed=12345,
        require_official_spec=False,
    )
    assert len(df_lines) > 0 and len(df_boot) > 0
    macro = summary["macro_means_by_method"]
    assert abs(macro["Transfer_Mean"]["ls_slope_ols_top50"]) < 1e-12
    assert macro["Hybrid_Raw"]["ls_slope_ols_top50"] > 0.7
    assert macro["Hybrid_Raw"]["ls_resid_corr_cal_top50"] > 0.8
    assert len(summary["pert_noise_ceilings"]) == 3
    assert len(summary["drf_calibration"]) == 3
    adj = adjudicate_results(df_lines, df_boot, summary, primary_method="Hybrid_Raw", min_lines_pass=2)
    assert adj["H-A1"]["verdict"] == "CONFIRMED"
    assert adj["H-A2"]["verdict"] == "CONFIRMED"
    md = render_markdown_report(adj, df_lines, summary, "unit_test")
    assert "CONFIRMED" in md
    assert macro["Hybrid_Raw"]["ls_slope_ols_top50"] > macro["NegCtrl_Hybrid_MismatchedCov"]["ls_slope_ols_top50"]
    assert macro["Hybrid_Raw"]["ls_slope_ols_top50"] > macro["NegCtrl_Hybrid_ShuffledGenes"]["ls_slope_ols_top50"]

    # Parity check between vectorized sufficient-stats bootstrap and literal concatenation
    rng_p = np.random.RandomState(33)
    d_obs = rng_p.randn(15, 70)
    d_src = rng_p.randn(15, 70)
    d_pred = d_src + 0.5 * (d_obs - d_src) + rng_p.randn(15, 70) * 0.1
    t_i = np.arange(15, dtype=np.int64)
    ls = compute_line_specific_metrics_and_boots(d_obs, d_pred, d_src, t_i, top_k=20, n_boot=50, seed=77)
    rs_ref = np.random.RandomState(77)
    draws_ref = rs_ref.choice(15, size=(50, 15), replace=True)
    r_true = d_obs - d_src
    r_pred = d_pred - d_src
    xs_lit, ys_lit = [], []
    for p in draws_ref[0]:
        valid = np.ones(70, dtype=bool)
        valid[p] = False
        idx = np.flatnonzero(valid)
        idx = idx[np.argsort(np.abs(d_obs[p, idx]))[-20:]]
        xs_lit.append(r_true[p, idx])
        ys_lit.append(r_pred[p, idx])
    xc = np.concatenate(xs_lit)
    yc = np.concatenate(ys_lit)
    slope_lit = float(np.dot(xc - xc.mean(), yc - yc.mean()) / np.dot(xc - xc.mean(), xc - xc.mean()))
    assert abs(ls["resp_top50"]["boot_slope_ols"][0] - slope_lit) < 1e-12
    assert abs(spearman_brown_ceiling(0.5)[0] - (2.0 / 3.0)) < 1e-12
    assert np.isnan(spearman_brown_ceiling(-0.1)[0])


if __name__ == "__main__":
    test_cipher_and_hybrid_planted_recovery()
    test_metrics_and_bootstrap()
    test_prereg_gate_and_splithalf()
    test_end_to_end_synthetic_loo()
    print("ALL_UNIT_TESTS_PASSED")
