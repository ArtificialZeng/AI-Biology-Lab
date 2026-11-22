#!/usr/bin/env python3
"""Unit tests for predictors and metrics using synthetic planted-truth systems."""

import numpy as np

from adjudicate import adjudicate_results, render_markdown_report
from data_loader import CODE_FREEZE_FILE, LineDataset, PREREG_COMMIT_FILE, check_prereg_gate, is_replogle_dev_line
from metrics import (
    cis_non_overlapping,
    compute_line_specific_metrics_and_boots,
    compute_residual_variance_decomposition,
    compute_row_metrics,
    compute_splithalf_noise_ceiling,
    compute_tmvs_shuffled_var_perm100,
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
    fit_shared_subspace,
    predict_zero_shot_suite,
    project_out_basis,
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

    # Replogle dev lines must pass gate even without PREREG_COMMIT.txt / CODE_FREEZE.txt
    check_prereg_gate(["replogle_k562", "replogle_rpe1"])
    if not (
        PREREG_COMMIT_FILE.exists()
        and PREREG_COMMIT_FILE.read_text(encoding="utf-8").strip()
        and CODE_FREEZE_FILE.exists()
        and CODE_FREEZE_FILE.read_text(encoding="utf-8").strip()
    ):
        blocked = False
        try:
            check_prereg_gate(["nadig_hepg2"])
        except RuntimeError:
            blocked = True
        assert blocked, "Non-Replogle line must be blocked before PREREG_COMMIT.txt & CODE_FREEZE.txt exist!"

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
    assert "NegCtrl_HybridProj_SrcCorrTgtVar" in macro
    assert "Transfer_Mean_VarScaled" in macro

    adj = adjudicate_results(
        df_lines, df_boot, summary, primary_method="Hybrid_ProjShared", min_lines_pass=2
    )
    assert adj["H-A1"]["verdict"] == "CONFIRMED"
    assert adj["H-A2"]["verdict"] == "CONFIRMED"
    assert adj["H-A3"]["verdict"] == "CONFIRMED"
    assert adj["H-A3"]["per_line"]["Synthetic_L3"]["n_sources_checked"] == 2
    assert "R4_Wording_Checks" in adj
    md = render_markdown_report(adj, df_lines, summary, "unit_test")
    assert "CONFIRMED" in md
    assert macro["Hybrid_ProjShared"]["ls_slope_ols_top50"] > macro["NegCtrl_HybridProj_ShuffledGenes"]["ls_slope_ols_top50"]

    # R6 & R11 identities check on Hybrid_ProjShared rows
    hp = df_lines[df_lines["method"] == "Hybrid_ProjShared"]
    for _, row in hp.iterrows():
        assert abs(row["ls_slope_ols_top50"] - row["ls_resid_corr_top50"] * row["ls_sd_ratio_top50"]) < 1e-12
        assert abs(row["ls_resid_corr_cal_top50_sq"] - row["ls_resid_corr_cal_top50"] ** 2) < 1e-12
        assert row["ls_slope_deming_top50"] > 0.0
        assert row["ls_slope_obs_on_pred_top50"] > 0.0
        assert np.isfinite(row["ctrl_disc_raw_panel"]) and np.isfinite(row["cent_disc_raw_panel"])

    # Parity check between vectorized sufficient-stats bootstrap and literal concatenation
    rng_p = np.random.RandomState(33)
    d_obs = rng_p.randn(15, 70)
    d_src = rng_p.randn(15, 70)
    d_pred = d_src + 0.5 * (d_obs - d_src) + rng_p.randn(15, 70) * 0.1
    t_i = np.arange(15, dtype=np.int64)
    ls = compute_line_specific_metrics_and_boots(d_obs, d_pred, d_src, t_i, top_k=20, n_boot=50, seed=77)
    rs_ref = np.random.default_rng(77)
    draws_ref = rs_ref.integers(0, 15, size=(50, 15))
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
    xc_c = xc - xc.mean()
    yc_c = yc - yc.mean()
    sxx = float(np.dot(xc_c, xc_c))
    syy = float(np.dot(yc_c, yc_c))
    sxy = float(np.dot(xc_c, yc_c))
    slope_lit = sxy / sxx
    deming_lit = float(((syy - sxx) + np.sqrt((syy - sxx) ** 2 + 4 * sxy ** 2)) / (2 * sxy))
    assert abs(ls["resp_top50"]["boot_slope_ols"][0] - slope_lit) < 1e-12
    assert abs(ls["resp_top50"]["boot_slope_deming"][0] - deming_lit) < 1e-12
    assert abs(ls["resp_top50"]["boot_sd_ratio"][0] - np.sqrt(syy / sxx)) < 1e-12
    assert abs(ls["resp_top50"]["boot_slope_obs_on_pred"][0] - (sxy / syy)) < 1e-12
    assert abs(spearman_brown_ceiling(0.5)[0] - (2.0 / 3.0)) < 1e-12
    assert np.isnan(spearman_brown_ceiling(-0.1)[0])


def test_r4_controls_and_confirmatory_threshold():
    src_lines, C_t, D_t, target_indices, target_names = _make_synthetic_system(seed=19)
    state = SourceModelState.fit(src_lines, mask_self=True)
    # Double marginal variance on target line covariance -> Transfer_Mean_VarScaled should reflect ratio sqrt(diag_t / diag_s)
    preds = predict_zero_shot_suite(
        state=state,
        target_ctrl_cov=4.0 * src_lines[0].ctrl_cov,
        query_target_names=target_names,
        query_target_gene_indices=target_indices,
        mismatched_ctrl_cov=None,
    )
    assert "NegCtrl_HybridProj_SrcCorrTgtVar" in preds
    assert "Transfer_Mean_VarScaled" in preds
    assert "NegCtrl_HybridProj_MismatchedCov(S1)" in preds
    assert "NegCtrl_HybridProj_MismatchedCov(S2)" in preds
    # Check when S1 and S2 have identical covariance, sqrt(4 C_s / C_s) == 2.0
    state_same = SourceModelState.fit(
        [
            src_lines[0],
            SourceLineInput("S2_clone", src_lines[0].ctrl_cov, src_lines[1].deltas, target_indices, target_names),
        ],
        mask_self=True,
    )
    preds_same = predict_zero_shot_suite(
        state=state_same,
        target_ctrl_cov=4.0 * src_lines[0].ctrl_cov,
        query_target_names=target_names,
        query_target_gene_indices=target_indices,
    )
    assert np.allclose(preds_same["Transfer_Mean_VarScaled"], 2.0 * preds_same["Transfer_Mean"], atol=1e-7)

    # Strict holdout confirmatory auto-selection test: replogle_* & <50 target lines become secondary
    import pandas as pd
    rows = []
    boots = []
    for ln, n_t in [("replogle_k562", 500), ("nadig_hepg2", 500), ("nadig_jurkat", 500), ("nourreddine_kolf21j", 18), ("vcc2025_h1", 80)]:
        for m, sl, cal, vfrac in [
            ("Hybrid_ProjShared", 0.20, 0.50, 1.0),
            ("CIPHER_Linear", 0.30, 0.40, 1.0),
            ("NegCtrl_HybridProj_MismatchedCov(src)", 0.01, 0.02, 1.0),
            ("NegCtrl_HybridProj_ShuffledGenes", 0.00, 0.00, 1.0),
            ("NegCtrl_HybridProj_SrcCorrTgtVar", 0.05, 0.10, 1.0),
            ("Transfer_Mean_VarScaled", 0.02, 0.05, 1.0),
        ]:
            rows.append({
                "target_line": ln, "n_eval_targets": n_t, "method": m,
                "ls_slope_ols_top50": sl, "ls_slope_ols_top50_ci95_lo": sl - 0.05, "ls_slope_ols_top50_ci95_hi": sl + 0.05,
                "ls_slope_noint_top50": sl, "ls_slope_per_pert_top50": sl,
                "ls_resid_corr_cal_top50": cal, "ls_resid_corr_cal_top50_sq": cal ** 2,
                "ls_resid_corr_cal_top50_ci95_lo": cal - 0.1, "ls_resid_corr_cal_top50_resamp_ci95_lo": cal - 0.1,
                "ls_resid_corr_cal_top50_resamp_valid_frac": vfrac,
                "cent_disc_raw": 0.8, "cent_de50_r": 0.4, "ctrl_disc_raw": 0.8, "ctrl_de50_r": 0.4, "mse": 0.01, "ls_slope_ols_all": sl, "ls_resid_corr_top50": cal * 0.5,
            })
        boots.append({"target_line": ln, "metric": "ls_slope_ols_top50", "method_a": "Hybrid_ProjShared", "method_b": "NegCtrl_HybridProj_MismatchedCov(src)", "disjoint_ci95": True, "a_ci95_lo": 0.15, "a_ci95_hi": 0.25, "b_ci95_lo": -0.04, "b_ci95_hi": 0.06, "mean_a": 0.20, "mean_b": 0.01, "delta": 0.19, "ci95_lo": 0.10, "ci95_hi": 0.28, "ci_gt_zero": True, "ci_lt_zero": False})
        boots.append({"target_line": ln, "metric": "ls_slope_ols_top50", "method_a": "Hybrid_ProjShared", "method_b": "NegCtrl_HybridProj_ShuffledGenes", "disjoint_ci95": True, "a_ci95_lo": 0.15, "a_ci95_hi": 0.25, "b_ci95_lo": -0.04, "b_ci95_hi": 0.04, "mean_a": 0.20, "mean_b": 0.00, "delta": 0.20, "ci95_lo": 0.10, "ci95_hi": 0.28, "ci_gt_zero": True, "ci_lt_zero": False})
        boots.append({"target_line": ln, "metric": "ls_slope_ols_top50", "method_a": "Hybrid_ProjShared", "method_b": "NegCtrl_HybridProj_SrcCorrTgtVar", "disjoint_ci95": True, "a_ci95_lo": 0.15, "a_ci95_hi": 0.25, "b_ci95_lo": 0.01, "b_ci95_hi": 0.09, "mean_a": 0.20, "mean_b": 0.05, "delta": 0.15, "ci95_lo": 0.06, "ci95_hi": 0.24, "ci_gt_zero": True, "ci_lt_zero": False})
        boots.append({"target_line": ln, "metric": "ls_resid_corr_cal_top50", "method_a": "Hybrid_ProjShared", "method_b": "NegCtrl_HybridProj_SrcCorrTgtVar", "disjoint_ci95": True, "a_ci95_lo": 0.40, "a_ci95_hi": 0.60, "b_ci95_lo": 0.05, "b_ci95_hi": 0.15, "mean_a": 0.50, "mean_b": 0.10, "delta": 0.40, "ci95_lo": 0.25, "ci95_hi": 0.55, "ci_gt_zero": True, "ci_lt_zero": False})
        boots.append({"target_line": ln, "metric": "ls_resid_corr_cal_top50", "method_a": "Hybrid_ProjShared", "method_b": "Transfer_Mean_VarScaled", "disjoint_ci95": True, "a_ci95_lo": 0.40, "a_ci95_hi": 0.60, "b_ci95_lo": 0.00, "b_ci95_hi": 0.10, "mean_a": 0.50, "mean_b": 0.05, "delta": 0.45, "ci95_lo": 0.30, "ci95_hi": 0.60, "ci_gt_zero": True, "ci_lt_zero": False})
        for met in ["ctrl_disc_raw", "ctrl_de50_r", "ctrl_disc_raw_panel", "cent_disc_raw", "cent_de50_r", "cent_disc_raw_panel"]:
            # Make one secondary line inferior to ensure secondary exclusions work
            lt_zero = (ln == "nourreddine_kolf21j" and met == "ctrl_disc_raw")
            boots.append({"target_line": ln, "metric": met, "method_a": "Hybrid_ProjShared", "method_b": "Transfer_Mean", "disjoint_ci95": False, "a_ci95_lo": 0.7, "a_ci95_hi": 0.9, "b_ci95_lo": 0.7, "b_ci95_hi": 0.9, "mean_a": 0.8, "mean_b": 0.8, "delta": -0.05 if lt_zero else 0.01, "ci95_lo": -0.08 if lt_zero else -0.01, "ci95_hi": -0.02 if lt_zero else 0.03, "ci_gt_zero": False, "ci_lt_zero": lt_zero})
    adj_h = adjudicate_results(pd.DataFrame(rows), pd.DataFrame(boots), {"phase": "holdout", "lines": ["nadig_hepg2", "nadig_jurkat", "nourreddine_kolf21j", "replogle_k562", "vcc2025_h1"]})
    assert adj_h["confirmatory_lines"] == ["nadig_hepg2", "nadig_jurkat", "vcc2025_h1"]
    assert adj_h["H-A1"]["required"] == 3 and adj_h["H-A4"]["max_inferior_allowed"] == 0
    assert adj_h["OVERALL_VERDICT"] == "CONFIRMED"
    assert adj_h["R4_Wording_Checks"]["allowed_cofluctuation_phrase"] == "细胞系特异的共波动结构"


def test_v3_group_b_fct_hybrid_on_tmvs_s5_and_ci975():
    rng = np.random.RandomState(2026)
    G = 80
    N_ctrl = 300
    N_pert = 24
    genes = [f"G{i}" for i in range(G)]
    g2i = {g: i for i, g in enumerate(genes)}
    t_names = tuple(genes[:N_pert])
    t_idx = np.arange(N_pert, dtype=np.int64)

    # Plant distinct source vs target control libnorm means & variances -> test FCT geometric mean + TMVS
    mu_s1 = rng.uniform(0.5, 4.0, size=G)
    mu_s2 = rng.uniform(0.5, 4.0, size=G)
    mu_t = rng.uniform(0.5, 4.0, size=G)
    fc_true = rng.uniform(0.4, 1.8, size=(N_pert, G))
    d_s1 = np.log1p(mu_s1[None, :] * fc_true) - np.log1p(mu_s1[None, :])
    d_s2 = np.log1p(mu_s2[None, :] * fc_true) - np.log1p(mu_s2[None, :])
    d_t_fct_expected = np.log1p(mu_t[None, :] * fc_true) - np.log1p(mu_t[None, :])

    cov_s1 = np.diag(rng.uniform(0.1, 1.5, size=G))
    cov_s2 = np.diag(rng.uniform(0.1, 1.5, size=G))
    cov_t = np.diag(rng.uniform(0.1, 1.5, size=G))

    state = SourceModelState.fit(
        [
            SourceLineInput("S1", cov_s1, d_s1, t_idx, t_names, ctrl_mean_libnorm=mu_s1),
            SourceLineInput("S2", cov_s2, d_s2, t_idx, t_names, ctrl_mean_libnorm=mu_s2),
        ],
        mask_self=True,
    )
    preds = predict_zero_shot_suite(
        state=state,
        target_ctrl_cov=cov_t,
        query_target_names=t_names,
        query_target_gene_indices=t_idx,
        target_ctrl_mean_libnorm=mu_t,
    )
    assert "FCT" in preds
    assert "Hybrid_on_TMVS" in preds
    assert "NegCtrl_TMVS_ShuffledVar" in preds
    assert "NegCtrl_TMVS_MismatchedVar(S1)" in preds
    assert "NegCtrl_TMVS_MismatchedVar(S2)" in preds
    assert np.allclose(preds["FCT"], d_t_fct_expected, atol=1e-10)

    # Check ci975 wider than or equal to ci95 + S5 & Group B end-to-end synthetic LOO
    datasets = {}
    shared_prog = rng.randn(N_pert, G) * 0.25
    for idx_l, (lname, mu_l, scale_sd) in enumerate(
        [("Replogle_K562", mu_s1, 0.7), ("Replogle_RPE1", mu_s2, 1.0), ("Conf_L3", mu_t, 1.8)]
    ):
        lrng = np.random.RandomState(700 + idx_l)
        gene_sds = lrng.uniform(0.15, 0.8, size=G) * scale_sd
        B = (lrng.randn(G, 6) * 0.25) * gene_sds[:, None]
        X_ctrl = lrng.randn(N_ctrl, 6) @ B.T + lrng.randn(N_ctrl, G) * (0.1 * gene_sds[None, :])
        # True perturbations scaled by line gene std (exact TMVS mechanism)
        deltas = shared_prog * gene_sds[None, :]
        ha = deltas + lrng.randn(*deltas.shape) * 0.02
        hb = deltas + lrng.randn(*deltas.shape) * 0.02
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
            ctrl_mean_libnorm=mu_l,
        )

    df_lines, df_boot, summary = evaluate_datasets(
        datasets=datasets, n_boot=200, seed=12345, require_official_spec=False
    )
    assert "s5_expressed_top75_var" in summary
    assert summary["s5_expressed_top75_var"]["Conf_L3"]["n_genes_kept"] == G - (G // 4)
    # 97.5% CI must enclose 95% CI
    tmvs_row = df_lines[(df_lines["target_line"] == "Conf_L3") & (df_lines["method"] == "Transfer_Mean_VarScaled")].iloc[0]
    assert tmvs_row["ls_slope_ols_top50_ci975_lo"] <= tmvs_row["ls_slope_ols_top50_ci95_lo"] + 1e-12
    assert tmvs_row["ls_slope_ols_top50_ci975_hi"] >= tmvs_row["ls_slope_ols_top50_ci95_hi"] - 1e-12
    assert np.isfinite(tmvs_row["ls_resid_corr_cal_top50_resamp_ci975_lo"])

    adj_b = adjudicate_results(
        df_lines,
        df_boot,
        summary,
        confirmatory_lines=["Conf_L3"],
        min_lines_pass=1,
        max_lines_fail_ha4=0,
        ci_level=0.975,
    )
    assert adj_b["ci_level_used"] == 0.975
    assert adj_b["H-B1"]["verdict"] == "CONFIRMED"
    assert adj_b["H-B2"]["verdict"] == "CONFIRMED"
    assert adj_b["H-B3"]["verdict"] == "CONFIRMED"
    assert adj_b["H-B3"]["per_line"]["Conf_L3"]["n_mismatch_controls_checked"] == 2
    assert "NegCtrl_TMVS_ShuffledVar_perm100" in adj_b["H-B3"]["per_line"]["Conf_L3"]["by_control"]
    assert adj_b["Group_B_Verdict"] == "CONFIRMED"
    assert "Group A" in adj_b["R13_Dual_Summary"] and "Group B" in adj_b["R13_Dual_Summary"]
    md = render_markdown_report(adj_b, df_lines, summary, "v4_unit")
    assert "S5:" in md and "S6:" in md and "97.5%" in md and "R18:" in md and "R19:" in md


def test_v4_r17_r19_orth_globalscale_perm100_vardecomp():
    """Verify R17 (r_orth_cal, rho_AB_orth), R18 (GlobalScale, perm100, H-A3/B3), and R19 (3-layer decomp)."""
    rng = np.random.RandomState(404)
    P, G = 40, 80
    t_idx = np.arange(P, dtype=np.int64)
    t_names = tuple(f"G{i}" for i in t_idx)
    genes = tuple(f"G{i}" for i in range(G))
    g2i = {g: i for i, g in enumerate(genes)}

    X = rng.randn(P, G) * 0.4
    # Planted: pure global amplitude gain beta = 1.8 (zero per-gene variance pattern)
    Y_global = 1.8 * X
    # Check R19 pure global gain recovery
    vd_glob = compute_residual_variance_decomposition(
        delta_obs=Y_global,
        delta_src_mean=X,
        target_gene_idx=t_idx,
        top_k=50,
        n_boot=100,
        seed=7,
    )
    for mode in ("resp_top50", "all_nontarget"):
        d_m = vd_glob[mode]
        f1 = d_m["frac_line_gain"]["estimate"]
        f2 = d_m["frac_gene_gain"]["estimate"]
        f3 = d_m["frac_remaining"]["estimate"]
        assert abs((f1 + f2 + f3) - 1.0) < 1e-10
        assert f1 > 0.99 and d_m["beta_mult"]["estimate"] > 1.79

    # Planted: per-gene multiplicative gain m_g = 1.6 * var_ratio_g + remaining
    gene_m = rng.uniform(0.6, 2.4, size=G)
    Y_gene = X * gene_m[None, :]
    vd_gene = compute_residual_variance_decomposition(
        delta_obs=Y_gene,
        delta_src_mean=X,
        target_gene_idx=t_idx,
        top_k=50,
        n_boot=100,
        seed=8,
    )
    assert abs(
        vd_gene["all_nontarget"]["frac_line_gain"]["estimate"]
        + vd_gene["all_nontarget"]["frac_gene_gain"]["estimate"]
        + vd_gene["all_nontarget"]["frac_remaining"]["estimate"]
        - 1.0
    ) < 1e-10
    assert (
        vd_gene["all_nontarget"]["frac_line_gain"]["estimate"]
        + vd_gene["all_nontarget"]["frac_gene_gain"]["estimate"]
    ) > 0.99

    # Verify Transfer_Mean_GlobalScale gets positive raw resid_corr yet orthogonalizes to zero/nan!
    cov_s = np.diag(np.ones(G))
    cov_t = np.diag(np.square(gene_m))
    s1 = SourceLineInput(
        name="S1",
        ctrl_cov=cov_s,
        deltas=X,
        target_indices=t_idx,
        target_names=t_names,
    )
    state = SourceModelState.fit([s1])
    preds = predict_zero_shot_suite(state, cov_t, t_names, t_idx)
    expected_c_glob = float(np.sqrt(np.mean(np.diag(cov_t)) / np.mean(np.diag(cov_s))))
    assert np.allclose(preds["Transfer_Mean_GlobalScale"], expected_c_glob * X, atol=1e-10)

    ls_gs = compute_line_specific_metrics_and_boots(
        delta_obs=Y_gene,
        delta_pred=preds["Transfer_Mean_GlobalScale"],
        delta_src_mean=X,
        target_gene_idx=t_idx,
        top_k=50,
        n_boot=100,
        seed=9,
    )["resp_top50"]
    assert ls_gs["resid_corr"] > 0.25  # Raw r_resid sees global amplitude mismatch!
    assert (not np.isfinite(ls_gs["resid_corr_orth"])) or abs(ls_gs["resid_corr_orth"]) < 1e-8

    # Check perm100 distribution monotonicity
    ceil = compute_splithalf_noise_ceiling(
        delta_half_a=Y_gene + rng.randn(P, G) * 0.02,
        delta_half_b=Y_gene + rng.randn(P, G) * 0.02,
        delta_obs_full=Y_gene,
        delta_src_mean=X,
        target_gene_idx=t_idx,
        top_k=50,
        n_boot=100,
        seed=9,
    )
    assert ceil["resp_top50"]["rho_half_ab_orth"] > 0.9
    p100 = compute_tmvs_shuffled_var_perm100(
        delta_obs=Y_gene,
        delta_src_mean=X,
        target_gene_idx=t_idx,
        target_diag_var=np.diag(cov_t),
        mean_src_var=np.diag(cov_s),
        ceil_r_top50=ceil["resp_top50"]["ceiling_r"],
        ceil_r_orth_top50=ceil["resp_top50"]["ceiling_r_orth"],
        n_perm=100,
        top_k=50,
        seed=42,
    )
    ro = p100["r_orth_cal"]
    assert p100["n_perm"] == 100
    assert ro["p2_5"] <= ro["median"] <= ro["p97_5"] <= ro["p98_75"] <= ro["max"] + 1e-12

    # Check adjudicate R17/R18 negative gates (rho_ab_orth<=0 fails H-B2; perm100 p97.5 fail breaks H-B3; H-A3 AND)
    import pandas as pd
    df_lines = pd.DataFrame([
        {
            "target_line": "C1", "n_eval_targets": 60, "method": "Hybrid_ProjShared",
            "ls_slope_ols_top50": 0.20, "ls_slope_ols_top50_ci95_lo": 0.16, "ls_slope_ols_top50_ci95_hi": 0.24,
            "ls_slope_ols_top50_ci975_lo": 0.15, "ls_slope_ols_top50_ci975_hi": 0.25,
            "ls_resid_corr_cal_top50": 0.20, "ls_resid_corr_cal_top50_ci95_lo": 0.10, "ls_resid_corr_cal_top50_resamp_ci95_lo": 0.10, "ls_resid_corr_cal_top50_resamp_ci975_lo": 0.08, "ls_resid_corr_cal_top50_resamp_valid_frac": 1.0,
            "ls_resid_corr_orth_cal_top50": 0.15, "ls_resid_corr_orth_cal_top50_ci95_lo": 0.08, "ls_resid_corr_orth_cal_top50_resamp_ci95_lo": 0.08, "ls_resid_corr_orth_cal_top50_resamp_ci975_lo": 0.06, "ls_resid_corr_orth_cal_top50_resamp_valid_frac": 1.0,
            "ls_resid_corr_top50": 0.18, "ls_resid_corr_orth_top50": 0.14, "ctrl_disc_raw": 0.7, "ctrl_de50_r": 0.4, "mse": 0.1, "ls_slope_ols_all": 0.1,
        },
        {
            "target_line": "C1", "n_eval_targets": 60, "method": "CIPHER_Linear",
            "ls_slope_ols_top50": 0.18, "ls_slope_ols_top50_ci95_lo": 0.14, "ls_slope_ols_top50_ci95_hi": 0.22,
            "ls_slope_ols_top50_ci975_lo": 0.13, "ls_slope_ols_top50_ci975_hi": 0.23,
            "ls_resid_corr_cal_top50": 0.15, "ls_resid_corr_cal_top50_ci95_lo": 0.05, "ls_resid_corr_cal_top50_resamp_ci95_lo": 0.05, "ls_resid_corr_cal_top50_resamp_ci975_lo": 0.03, "ls_resid_corr_cal_top50_resamp_valid_frac": 1.0,
            "ls_resid_corr_orth_cal_top50": 0.10, "ls_resid_corr_orth_cal_top50_ci95_lo": 0.02, "ls_resid_corr_orth_cal_top50_resamp_ci95_lo": 0.02, "ls_resid_corr_orth_cal_top50_resamp_ci975_lo": 0.01, "ls_resid_corr_orth_cal_top50_resamp_valid_frac": 1.0,
            "ls_resid_corr_top50": 0.14, "ls_resid_corr_orth_top50": 0.09, "ctrl_disc_raw": 0.7, "ctrl_de50_r": 0.4, "mse": 0.1, "ls_slope_ols_all": 0.1,
        },
        {
            "target_line": "C1", "n_eval_targets": 60, "method": "Transfer_Mean_VarScaled",
            "ls_slope_ols_top50": 0.20, "ls_slope_ols_top50_ci95_lo": 0.16, "ls_slope_ols_top50_ci95_hi": 0.24,
            "ls_slope_ols_top50_ci975_lo": 0.15, "ls_slope_ols_top50_ci975_hi": 0.25,
            "ls_resid_corr_cal_top50": 0.25, "ls_resid_corr_cal_top50_ci95_lo": 0.15, "ls_resid_corr_cal_top50_resamp_ci95_lo": 0.15, "ls_resid_corr_cal_top50_resamp_ci975_lo": 0.12, "ls_resid_corr_cal_top50_resamp_valid_frac": 1.0,
            "ls_resid_corr_orth_cal_top50": 0.15, "ls_resid_corr_orth_cal_top50_ci95_lo": 0.08, "ls_resid_corr_orth_cal_top50_resamp_ci95_lo": 0.08, "ls_resid_corr_orth_cal_top50_resamp_ci975_lo": 0.06, "ls_resid_corr_orth_cal_top50_resamp_valid_frac": 1.0,
            "ls_resid_corr_top50": 0.23, "ls_resid_corr_orth_top50": 0.14, "ctrl_disc_raw": 0.7, "ctrl_de50_r": 0.4, "mse": 0.1, "ls_slope_ols_all": 0.1,
        },
        {
            "target_line": "C1", "n_eval_targets": 60, "method": "Transfer_Mean_GlobalScale",
            "ls_slope_ols_top50": 0.30, "ls_slope_ols_top50_ci95_lo": 0.25, "ls_slope_ols_top50_ci95_hi": 0.35,
            "ls_slope_ols_top50_ci975_lo": 0.24, "ls_slope_ols_top50_ci975_hi": 0.36,
            "ls_resid_corr_cal_top50": 0.50, "ls_resid_corr_cal_top50_ci95_lo": 0.40, "ls_resid_corr_cal_top50_resamp_ci95_lo": 0.40, "ls_resid_corr_cal_top50_resamp_ci975_lo": 0.38, "ls_resid_corr_cal_top50_resamp_valid_frac": 1.0,
            "ls_resid_corr_orth_cal_top50": 0.0, "ls_resid_corr_orth_cal_top50_ci95_lo": 0.0, "ls_resid_corr_orth_cal_top50_resamp_ci95_lo": 0.0, "ls_resid_corr_orth_cal_top50_resamp_ci975_lo": 0.0, "ls_resid_corr_orth_cal_top50_resamp_valid_frac": 1.0,
            "ls_resid_corr_top50": 0.48, "ls_resid_corr_orth_top50": 0.0, "ctrl_disc_raw": 0.7, "ctrl_de50_r": 0.4, "mse": 0.1, "ls_slope_ols_all": 0.1,
        },
    ])
    # H-A3 beats slope (0.15 > 0.04) but fails r_orth_cal (0.06 <= 0.07) -> MUST fail H-A3!
    # H-B3 perm100 p97_5 = 0.07 >= 0.06 -> MUST fail H-B3!
    # H-B2 fails when rho_half_ab_orth <= 0!
    df_boot = pd.DataFrame([
        {"target_line": "C1", "metric": "ls_slope_ols_top50", "method_a": "Hybrid_ProjShared", "method_b": "NegCtrl_HybridProj_MismatchedCov(S1)", "disjoint_ci95": True, "disjoint_ci975": True, "a_ci95_lo": 0.16, "a_ci95_hi": 0.24, "b_ci95_lo": 0.0, "b_ci95_hi": 0.03, "a_ci975_lo": 0.15, "a_ci975_hi": 0.25, "b_ci975_lo": -0.01, "b_ci975_hi": 0.04, "mean_a": 0.20, "mean_b": 0.01, "delta": 0.19, "ci95_lo": 0.12, "ci95_hi": 0.26, "ci975_lo": 0.11, "ci975_hi": 0.27, "ci_gt_zero": True, "ci_lt_zero": False, "ci975_gt_zero": True, "ci975_lt_zero": False},
        {"target_line": "C1", "metric": "ls_resid_corr_orth_cal_top50", "method_a": "Hybrid_ProjShared", "method_b": "NegCtrl_HybridProj_MismatchedCov(S1)", "disjoint_ci95": True, "disjoint_ci975": False, "a_ci95_lo": 0.08, "a_ci95_hi": 0.22, "b_ci95_lo": -0.02, "b_ci95_hi": 0.06, "a_ci975_lo": 0.06, "a_ci975_hi": 0.24, "b_ci975_lo": -0.03, "b_ci975_hi": 0.07, "mean_a": 0.15, "mean_b": 0.02, "delta": 0.13, "ci95_lo": 0.03, "ci95_hi": 0.23, "ci975_lo": 0.02, "ci975_hi": 0.24, "ci_gt_zero": True, "ci_lt_zero": False, "ci975_gt_zero": True, "ci975_lt_zero": False},
        {"target_line": "C1", "metric": "ls_slope_ols_top50", "method_a": "Transfer_Mean_VarScaled", "method_b": "NegCtrl_TMVS_MismatchedVar(S1)", "disjoint_ci95": True, "disjoint_ci975": True, "a_ci95_lo": 0.16, "a_ci95_hi": 0.24, "b_ci95_lo": 0.0, "b_ci95_hi": 0.03, "a_ci975_lo": 0.15, "a_ci975_hi": 0.25, "b_ci975_lo": -0.01, "b_ci975_hi": 0.04, "mean_a": 0.20, "mean_b": 0.01, "delta": 0.19, "ci95_lo": 0.12, "ci95_hi": 0.26, "ci975_lo": 0.11, "ci975_hi": 0.27, "ci_gt_zero": True, "ci_lt_zero": False, "ci975_gt_zero": True, "ci975_lt_zero": False},
        {"target_line": "C1", "metric": "ls_resid_corr_orth_cal_top50", "method_a": "Transfer_Mean_VarScaled", "method_b": "NegCtrl_TMVS_MismatchedVar(S1)", "disjoint_ci95": True, "disjoint_ci975": True, "a_ci95_lo": 0.08, "a_ci95_hi": 0.22, "b_ci95_lo": -0.03, "b_ci95_hi": 0.03, "a_ci975_lo": 0.06, "a_ci975_hi": 0.24, "b_ci975_lo": -0.04, "b_ci975_hi": 0.04, "mean_a": 0.15, "mean_b": 0.00, "delta": 0.15, "ci95_lo": 0.05, "ci95_hi": 0.25, "ci975_lo": 0.04, "ci975_hi": 0.26, "ci_gt_zero": True, "ci_lt_zero": False, "ci975_gt_zero": True, "ci975_lt_zero": False},
        {"target_line": "C1", "metric": "ctrl_disc_raw", "method_a": "Hybrid_ProjShared", "method_b": "Transfer_Mean", "disjoint_ci95": False, "disjoint_ci975": False, "a_ci95_lo": 0.65, "a_ci95_hi": 0.75, "b_ci95_lo": 0.65, "b_ci95_hi": 0.75, "a_ci975_lo": 0.64, "a_ci975_hi": 0.76, "b_ci975_lo": 0.64, "b_ci975_hi": 0.76, "mean_a": 0.7, "mean_b": 0.7, "delta": 0.0, "ci95_lo": -0.01, "ci95_hi": 0.01, "ci975_lo": -0.02, "ci975_hi": 0.02, "ci_gt_zero": False, "ci_lt_zero": False, "ci975_gt_zero": False, "ci975_lt_zero": False},
        {"target_line": "C1", "metric": "ctrl_disc_raw", "method_a": "Transfer_Mean_VarScaled", "method_b": "Transfer_Mean", "disjoint_ci95": False, "disjoint_ci975": False, "a_ci95_lo": 0.65, "a_ci95_hi": 0.75, "b_ci95_lo": 0.65, "b_ci95_hi": 0.75, "a_ci975_lo": 0.64, "a_ci975_hi": 0.76, "b_ci975_lo": 0.64, "b_ci975_hi": 0.76, "mean_a": 0.7, "mean_b": 0.7, "delta": 0.0, "ci95_lo": -0.01, "ci95_hi": 0.01, "ci975_lo": -0.02, "ci975_hi": 0.02, "ci_gt_zero": False, "ci_lt_zero": False, "ci975_gt_zero": False, "ci975_lt_zero": False},
    ])
    summary_neg = {
        "pert_noise_ceilings": {"C1": {"resp_top50": {"rho_half_ab_orth": -0.05}}},
        "tmvs_shuffled_var_perm100": {"C1": {"n_perm": 100, "r_orth_cal": {"p2_5": -0.04, "median": 0.01, "p97_5": 0.07, "p98_75": 0.08, "max": 0.09}}},
    }
    adj_neg = adjudicate_results(df_lines, df_boot, summary_neg, confirmatory_lines=["C1"], min_lines_pass=1, max_lines_fail_ha4=0, ci_level=0.975)
    assert adj_neg["H-A2"]["verdict"] == "未确认"
    assert adj_neg["H-B2"]["verdict"] == "未确认"
    assert adj_neg["H-A3"]["verdict"] == "未确认"
    assert adj_neg["H-B3"]["verdict"] == "未确认"
    assert "Transfer_Mean_GlobalScale" in adj_neg["R18_GlobalScale_Summary"]["summary_statement"]


def test_s2b_rank1_q1_projection_mode() -> None:
    rng = np.random.default_rng(123)
    N, G = 40, 24
    Y = rng.normal(0.0, 0.2, size=(N, G)) + np.linspace(0.5, -0.5, G)[None, :]
    b_def, a_def, Q_def = fit_shared_subspace(Y)
    b_sp3, a_sp3, Q_sp3 = fit_shared_subspace(Y, proj_mode="span3")
    assert np.array_equal(b_def, b_sp3)
    assert np.array_equal(a_def, a_sp3)
    assert np.array_equal(Q_def, Q_sp3)
    assert Q_def.shape == (G, 3)

    b_r1, q1, Q_r1 = fit_shared_subspace(Y, proj_mode="rank1_q1")
    assert np.array_equal(b_r1, b_def)
    assert Q_r1.shape == (G, 1)
    assert np.allclose(np.linalg.norm(q1), 1.0, atol=1e-12)
    assert np.allclose(Q_r1[:, 0], q1, atol=1e-12)
    _, _, Vt = np.linalg.svd(Y - Y.mean(axis=0, keepdims=True), full_matrices=False)
    expected_q1 = Vt[0].copy()
    if expected_q1[np.argmax(np.abs(expected_q1))] < 0:
        expected_q1 = -expected_q1
    assert np.allclose(q1, expected_q1, atol=1e-12)
    assert float(np.linalg.norm(project_out_basis(q1[None, :], Q_r1))) < 1e-12

    cov_s = np.eye(G) + 0.2 * np.outer(q1, q1)
    cov_t = np.diag(rng.uniform(0.5, 1.5, size=G)) + 0.3 * np.outer(q1, q1)
    src = SourceLineInput(
        name="S1",
        target_names=[f"T{i}" for i in range(N)],
        target_indices=np.arange(N, dtype=np.int64) % G,
        deltas=Y,
        ctrl_cov=cov_s,
    )
    st_sp3 = SourceModelState.fit([src], mask_self=True, proj_mode="span3")
    st_r1 = SourceModelState.fit([src], mask_self=True, proj_mode="rank1_q1")
    preds_sp3 = predict_zero_shot_suite(st_sp3, cov_t, [f"T{i}" for i in range(10)], np.arange(10, dtype=np.int64))
    preds_r1 = predict_zero_shot_suite(st_r1, cov_t, [f"T{i}" for i in range(10)], np.arange(10, dtype=np.int64))
    for inv_k in ("Transfer_Mean", "CIPHER_Linear", "Hybrid_Raw", "Transfer_Mean_VarScaled", "Transfer_Mean_GlobalScale"):
        assert np.allclose(preds_sp3[inv_k], preds_r1[inv_k], atol=1e-12)
    assert not np.allclose(preds_sp3["Hybrid_ProjShared"], preds_r1["Hybrid_ProjShared"])
    assert not np.allclose(preds_sp3["Hybrid_on_TMVS"], preds_r1["Hybrid_on_TMVS"])


if __name__ == "__main__":
    test_cipher_and_hybrid_planted_recovery()
    test_metrics_and_bootstrap()
    test_prereg_gate_and_splithalf()
    test_end_to_end_synthetic_loo()
    test_r4_controls_and_confirmatory_threshold()
    test_v3_group_b_fct_hybrid_on_tmvs_s5_and_ci975()
    test_v4_r17_r19_orth_globalscale_perm100_vardecomp()
    test_s2b_rank1_q1_projection_mode()
    print("ALL_UNIT_TESTS_PASSED")
