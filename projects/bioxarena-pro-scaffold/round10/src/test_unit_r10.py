#!/usr/bin/env python3
"""Unit tests for Round 10 Group C (`FCT`), HDF5 stream loaders, HGNC bridge, R19, S7 (`P6`), and `PREREG_C` gates."""

from __future__ import annotations

import gzip
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

from adjudicate_r10 import SIX_CONVENTIONS, adjudicate_r10_results
from data_loader_r10 import (
    DEFAULT_HGNC_FILE,
    HGNCBridge,
    LineDatasetR10,
    check_prereg_c_gate,
    count_igvf_single_sgrna_from_metadata,
    iter_h5ad_csr_chunks_r10,
    load_dld1_gsm_channel,
    parse_dld1_guide_gene,
    read_10x_feature_bc_matrix_h5,
    read_cd4t_flex_h5,
    split_units_even_odd,
)
from fct_pipeline import (
    FctSourceState,
    SourceLineInputR10,
    compute_source_transfer_and_fc_bar,
    iter_fct_perm_predictions,
    predict_group_c_suite,
)
from metrics_r10 import (
    compute_fct_perm100_null,
    compute_r19_variance_decomposition_r10,
    compute_residual_suite_r10,
    compute_row_metrics_r10,
    compute_s7_cellcount_and_kd_analysis,
    paired_row_bootstrap_ci_r10,
)

R9_SRC_DIR = Path(__file__).resolve().parents[2] / "round9" / "p2" / "vcc_linespecific" / "src"


def _write_10x_h5(path: Path, X_csr: sp.csr_matrix, barcodes: list[str], feat_names: list[str], feat_ids: list[str], feat_types: list[str]) -> None:
    """Create a CellRanger-layout CSC HDF5 matrix (n_features x n_cells) for loader unit tests."""
    csc = X_csr.T.tocsc()  # (n_features, n_cells)
    with h5py.File(path, "w") as f:
        m = f.create_group("matrix")
        m.create_dataset("barcodes", data=np.asarray(barcodes, dtype="S"))
        m.create_dataset("data", data=csc.data.astype(np.int32))
        m.create_dataset("indices", data=csc.indices.astype(np.int32))
        m.create_dataset("indptr", data=csc.indptr.astype(np.int64))
        m.create_dataset("shape", data=np.asarray(csc.shape, dtype=np.int64))
        feat = m.create_group("features")
        feat.create_dataset("name", data=np.asarray(feat_names, dtype="S"))
        feat.create_dataset("id", data=np.asarray(feat_ids, dtype="S"))
        feat.create_dataset("feature_type", data=np.asarray(feat_types, dtype="S"))


def test_hgnc_bridge_and_gate() -> None:
    if DEFAULT_HGNC_FILE.exists():
        hgnc = HGNCBridge.from_tsv(DEFAULT_HGNC_FILE)
        assert hgnc.canonicalize("AARS") == "AARS1"
        assert hgnc.canonicalize("ADAL") == "MAPDA"
        assert hgnc.canonicalize("OLD_SYM", ensembl_id="ENSG00000141510.18") == "TP53"
        assert hgnc.canonicalize("TP53") == "TP53"
        syms, idx, _ = hgnc.dedup_map_indices(
            ["AARS", "AARS1", "ADAL", "UNKNOWN_X"],
            priority_scores=[1.0, 5.0, 2.0, 0.5],
        )
        assert syms == ["AARS1", "MAPDA", "UNKNOWN_X"]
        assert idx.tolist() == [1, 2, 3]

    with tempfile.TemporaryDirectory() as tmp:
        st_dir = Path(tmp)
        try:
            check_prereg_c_gate(st_dir)
            raise AssertionError("Gate must fail when lock files missing")
        except RuntimeError:
            pass
        (st_dir / "PREREG_C_COMMIT.txt").write_text("dc4b4a5\n", encoding="utf-8")
        try:
            check_prereg_c_gate(st_dir)
            raise AssertionError("Gate must fail without CODE_FREEZE_R10.txt")
        except RuntimeError:
            pass
        (st_dir / "CODE_FREEZE_R10.txt").write_text("deadbeef\n", encoding="utf-8")
        p_h, f_h = check_prereg_c_gate(st_dir)
        assert p_h == "dc4b4a5" and f_h == "deadbeef"


def test_hdf5_stream_loaders_and_splits() -> None:
    # Unit split-half even vs odd
    ha, hb = split_units_even_odd([f"SAM_{i:02d}" for i in range(1, 7)])
    assert ha == {"SAM_01", "SAM_03", "SAM_05"}
    assert hb == {"SAM_02", "SAM_04", "SAM_06"}

    # DLD-1 regex parser
    assert parse_dld1_guide_gene("NTC_ctrl_12") == (True, None)
    assert parse_dld1_guide_gene("BRD4_ENST00000263377(W)_1") == (False, "BRD4")
    assert parse_dld1_guide_gene("MYC_P1P2_3") == (False, "MYC")

    with tempfile.TemporaryDirectory() as tmp:
        td = Path(tmp)
        # 1. DLD-1 GSM trio: 4 cells; cell0=doublet hash -> drop, cell1=double guide -> drop, cell2=TP53 singlet, cell3=NTC singlet
        bcs = ["BC0", "BC1", "BC2", "BC3"]
        gex_csr = sp.csr_matrix(np.arange(12, dtype=np.int32).reshape(4, 3) + 1)
        _write_10x_h5(
            td / "gex.h5", gex_csr, bcs,
            ["TP53", "MYC", "GAPDH"], ["E1", "E2", "E3"], ["Gene Expression"] * 3,
        )
        crispr_mat = np.array([[5, 0, 0], [4, 4, 0], [8, 0, 0], [0, 0, 7]], dtype=np.int32)
        _write_10x_h5(
            td / "crispr.h5", sp.csr_matrix(crispr_mat), bcs,
            ["TP53_P1_1", "MYC_P1_1", "NTC_01"], ["g1", "g2", "g3"], ["CRISPR Guide Capture"] * 3,
        )
        hash_mat = np.array([[6, 5], [7, 0], [9, 0], [0, 8]], dtype=np.int32)
        _write_10x_h5(
            td / "hash.h5", sp.csr_matrix(hash_mat), bcs,
            ["CMO1", "CMO2"], ["h1", "h2"], ["Multiplexing Capture"] * 2,
        )
        bundle, is_ntc, raw_tgts, chans = load_dld1_gsm_channel(
            td / "gex.h5", td / "crispr.h5", td / "hash.h5", "SAM24482910_1"
        )
        assert bundle.barcodes == ["SAM24482910_1:BC2", "SAM24482910_1:BC3"]
        assert is_ntc == [False, True]
        assert raw_tgts == ["TP53", None]
        assert chans == ["SAM24482910_1", "SAM24482910_1"]
        assert bundle.X_csr.shape == (2, 3)

        # 2. CD4+ T Flex combined `.h5` reader
        flex_X = sp.hstack([
            gex_csr,
            sp.csr_matrix(np.array([[10, 0], [5, 5], [0, 12], [2, 0]], dtype=np.int32)),
        ]).tocsr()
        _write_10x_h5(
            td / "flex.h5", flex_X, bcs,
            ["TP53", "MYC", "GAPDH", "SG_CD4_1", "SG_NTC_1"],
            ["E1", "E2", "E3", "SG_CD4_1", "SG_NTC_1"],
            ["Gene Expression"] * 3 + ["CRISPR Guide Capture"] * 2,
        )
        sgrna_df = pd.DataFrame([
            {"sgrna_id": "SG_CD4_1", "target_gene_name": "TP53", "target_gene_id": "ENSG00000141510", "guide_type": "targeting"},
            {"sgrna_id": "SG_NTC_1", "target_gene_name": "NTC", "target_gene_id": "NTC", "guide_type": "non-targeting"},
        ])
        f_bun, f_ntc, f_tgts, _ = read_cd4t_flex_h5(td / "flex.h5", sgrna_meta_df=sgrna_df, min_guide_umi=3)
        # BC1 is multi-guide -> dropped; BC3 has only UMI=2 < 3 -> dropped; keeps BC0 (TP53) and BC2 (NTC)
        assert f_bun.barcodes == ["BC0", "BC2"]
        assert f_ntc.tolist() == [False, True]
        assert f_tgts == ["TP53", None]

        # 3. AnnData CSR `.h5ad` streaming chunks
        csr_test = sp.csr_matrix(np.arange(30, dtype=np.float32).reshape(6, 5))
        with h5py.File(td / "test.h5ad", "w") as hf:
            xg = hf.create_group("X")
            xg.attrs["encoding-type"] = "csr_matrix"
            xg.attrs["shape"] = np.asarray([6, 5], dtype=np.int64)
            xg.create_dataset("data", data=csr_test.data)
            xg.create_dataset("indices", data=csr_test.indices)
            xg.create_dataset("indptr", data=csr_test.indptr)
        chunks = [blk for _, _, blk in iter_h5ad_csr_chunks_r10(td / "test.h5ad", row_chunk=2)]
        assert len(chunks) == 3
        assert np.allclose(np.vstack(chunks), csr_test.toarray())

        # 4. IGVF P2 metadata pre-check without expression data
        igvf_csv = td / "igvf_meta.csv.gz"
        with gzip.open(igvf_csv, "wt", encoding="utf-8") as gz:
            gz.write("guide_id,intended_target_name,type\n")
            for i in range(52):
                gz.write(f"g_{i},T_{i},targeting\n")
            gz.write("g_ntc,non-targeting,non-targeting\n")
        if DEFAULT_HGNC_FILE.exists():
            hgnc = HGNCBridge.from_tsv(DEFAULT_HGNC_FILE)
            cells = [f"g_{i}" for i in range(51) for _ in range(30)] + ["g_51"] * 10 + ["g_ntc"] * 40
            r9_set = {f"T_{i}" for i in range(52)}
            p2 = count_igvf_single_sgrna_from_metadata(igvf_csv, cells, r9_set, hgnc, min_cells=25)
            assert p2["n_shared_ge_min_cells"] == 51 and p2["passes_p2_50_target_gate"]


def test_fct_parity_with_round9_and_adjudication() -> None:
    rng = np.random.default_rng(20261122)
    N, G = 45, 90
    q_names = [f"G{i}" for i in range(N)]
    q_gidx = np.arange(N, dtype=np.int64)

    # Plant 6 source lines + 1 confirmatory target line where target response follows exact FCT fold-change
    src_inputs = []
    mu_sources = [rng.uniform(0.5, 4.0, size=G) for _ in range(6)]
    log_fc_true = rng.normal(0.0, 0.25, size=(N, G))

    for s_i in range(6):
        mu_s = mu_sources[s_i]
        delta_s = np.log1p(mu_s[None, :] * np.exp(log_fc_true)) - np.log1p(mu_s[None, :])
        delta_s += rng.normal(0.0, 0.015, size=(N, G))
        cov_s = np.diag(rng.uniform(0.2, 1.0, size=G))
        src_inputs.append(
            SourceLineInputR10(
                name=f"SRC_{s_i}",
                target_names=q_names,
                target_indices=q_gidx,
                deltas=delta_s,
                ctrl_cov=cov_s,
                ctrl_mean_libnorm=mu_s,
            )
        )

    mu_t = rng.uniform(0.2, 8.0, size=G)
    signal_t = np.log1p(mu_t[None, :] * np.exp(log_fc_true)) - np.log1p(mu_t[None, :])
    ha_t = signal_t + rng.normal(0.0, 0.02, size=(N, G))
    hb_t = signal_t + rng.normal(0.0, 0.02, size=(N, G))
    obs_t = 0.5 * (ha_t + hb_t)
    cov_t = np.diag(rng.uniform(0.2, 1.0, size=G))

    state = FctSourceState.fit(src_inputs, mask_self=True)
    preds, fc_bundle = predict_group_c_suite(
        state=state,
        target_ctrl_mean_libnorm=mu_t,
        query_target_names=q_names,
        query_target_gene_indices=q_gidx,
        target_ctrl_cov=cov_t,
        include_r9_anchors=True,
    )

    # Verify exact <1e-12 numerical parity with Round 9 FCT implementation
    if R9_SRC_DIR.exists():
        sys.path.insert(0, str(R9_SRC_DIR))
        from predictors import SourceLineInput as R9SrcIn, SourceModelState as R9State, predict_zero_shot_suite as r9_pred
        r9_srcs = [
            R9SrcIn(
                name=s.name,
                ctrl_cov=s.ctrl_cov,
                deltas=s.deltas,
                target_indices=s.target_indices,
                target_names=tuple(s.target_names),
                ctrl_mean_libnorm=s.ctrl_mean_libnorm,
            )
            for s in src_inputs
        ]
        r9_p = r9_pred(R9State.fit(r9_srcs), cov_t, q_names, q_gidx, target_ctrl_mean_libnorm=mu_t)
        assert float(np.max(np.abs(preds["FCT"] - r9_p["FCT"]))) < 1e-12
        assert float(np.max(np.abs(preds["Hybrid_ProjShared"] - r9_p["Hybrid_ProjShared"]))) < 1e-12

    # Check GlobalScale orthogonalizes cleanly to NaN (r_orth_raw is ~0 or degenerate NaN)
    gs_suite = compute_residual_suite_r10(
        obs_t, preds["Transfer_Mean_GlobalScale"], preds["Transfer_Mean"], q_gidx,
        delta_half_a=ha_t, delta_half_b=hb_t, n_boot=50, seed=20261122
    )
    assert np.isnan(gs_suite["resp_top50"]["r_orth_cal"]["estimate"])

    # Compute FCT + 6 mismatched means + 100 perm + 6 conventions -> must CONFIRM H-C1/H-C2/H-C3
    rows_list = []
    boot_list = []
    m_row_cache = {}
    for m_name, p_mat in preds.items():
        suite = compute_residual_suite_r10(
            obs_t, p_mat, preds["Transfer_Mean"], q_gidx,
            delta_half_a=ha_t, delta_half_b=hb_t, n_boot=200, seed=20261122
        )
        t50 = suite["resp_top50"]
        ro = t50["r_orth_cal"]
        for c_line in ("CONF_1", "CONF_2"):
            rows_list.append({
                "target_line": c_line,
                "method": m_name,
                "n_eval_targets": N,
                "r_orth_cal": ro["estimate"],
                "r_orth_cal_ci975_lo": ro["ci975_lo"],
                "r_orth_cal_ci975_hi": ro["ci975_hi"],
                "r_orth_cal_valid_frac": ro["valid_frac"],
                "rho_half_ab_orth": t50["rho_half_ab_orth"],
            })
        rm_c = compute_row_metrics_r10(p_mat, obs_t, q_gidx, centroid=None)
        rm_m = compute_row_metrics_r10(p_mat, obs_t, q_gidx, centroid=state.grand_mean)
        m_row_cache[m_name] = {
            "ctrl_disc_raw": rm_c["disc_raw"],
            "ctrl_disc_raw_panel": rm_c["disc_raw_panel"],
            "ctrl_de50_r": rm_c["de50_r"],
            "cent_disc_raw": rm_m["disc_raw"],
            "cent_disc_raw_panel": rm_m["disc_raw_panel"],
            "cent_de50_r": rm_m["de50_r"],
        }
    for conv in SIX_CONVENTIONS:
        ci = paired_row_bootstrap_ci_r10(m_row_cache["FCT"][conv], m_row_cache["Transfer_Mean"][conv], n_boot=200, seed=20261122)
        for c_line in ("CONF_1", "CONF_2"):
            boot_list.append({"target_line": c_line, "metric": conv, "method_a": "FCT", "method_b": "Transfer_Mean", **ci})

    fct_r = next(r for r in rows_list if r["method"] == "FCT")
    p100 = compute_fct_perm100_null(
        state, mu_t, q_names, obs_t, preds["Transfer_Mean"], q_gidx,
        ceil_r_top50=0.9, ceil_r_orth_top50=0.9, fc_bundle=fc_bundle, n_perm=100
    )
    assert p100["n_perm"] == 100
    summary = {"phase": "test", "fct_perm100": {"CONF_1": p100, "CONF_2": p100}}
    adj_pass = adjudicate_r10_results(pd.DataFrame(rows_list), pd.DataFrame(boot_list), summary, ["CONF_1", "CONF_2"])
    assert adj_pass["Group_C_Verdict"] == "CONFIRMED"
    assert adj_pass["H-C1"]["verdict"] == "CONFIRMED"
    assert adj_pass["H-C2"]["verdict"] == "CONFIRMED"
    assert adj_pass["H-C3"]["verdict"] == "CONFIRMED"
    assert adj_pass["k_conf_disclosure"] == "K_conf = 2 confirmatory cell lines"

    # P1 K_conf < 2 guard check
    adj_k1 = adjudicate_r10_results(pd.DataFrame(rows_list), pd.DataFrame(boot_list), summary, ["CONF_1"])
    assert adj_k1["confirmatory_status"] == "INSUFFICIENT_LINES"
    assert "确认性线不足" in adj_k1["k_conf_disclosure"]

    # Null signal must NOT confirm H-C1/H-C2
    null_rows = [dict(r) for r in rows_list]
    for r in null_rows:
        if r["method"] == "FCT":
            r["r_orth_cal"] = 0.01
            r["r_orth_cal_ci975_lo"] = -0.02
    adj_fail = adjudicate_r10_results(pd.DataFrame(null_rows), pd.DataFrame(boot_list), summary, ["CONF_1", "CONF_2"])
    assert adj_fail["Group_C_Verdict"] == "未确认"

    # 100-perm generator parity
    perm_stream = list(iter_fct_perm_predictions(state, mu_t, q_names, n_perm=3, precomputed=fc_bundle))
    assert np.allclose(perm_stream[0][1], preds["NegCtrl_FCT_Perm0"], atol=1e-12)


def test_s7_cellcount_kd_and_r19_unbiased() -> None:
    rng = np.random.default_rng(999)
    N, G = 36, 70
    q_names = [f"G{i}" for i in range(N)]
    q_gidx = np.arange(N, dtype=np.int64)
    mu_s = rng.uniform(1.0, 4.0, size=G)
    mu_t = rng.uniform(0.5, 6.0, size=G)
    t_mean = rng.normal(0.0, 0.3, size=(N, G))

    # Plant per-perturbation residual signal that scales strongly with n_cells & |kd| + true off-transfer pattern
    nc_vec = np.linspace(30.0, 3000.0, N)
    kd_vec = np.linspace(0.2, 2.0, N)
    base_dir = rng.normal(0.0, 1.0, size=(N, G))
    pred_dir = base_dir.copy()
    obs = np.empty((N, G), dtype=np.float64)
    pred = t_mean + 0.1 * pred_dir
    for p in range(N):
        sig_strength = 0.05 + 0.45 * (p / float(N - 1))
        noise_strength = 0.25 * (1.0 - p / float(N - 1)) + 0.02
        obs[p] = t_mean[p] + sig_strength * base_dir[p] + rng.normal(0.0, noise_strength, size=G)
        obs[p, q_gidx[p]] = -kd_vec[p]

    # Build dummy state for tertile perm100
    src = SourceLineInputR10("S0", q_names, q_gidx, t_mean, np.eye(G), mu_s)
    state = FctSourceState.fit([src])
    fc_bundle = compute_source_transfer_and_fc_bar(state, q_names)
    s7 = compute_s7_cellcount_and_kd_analysis(
        state=state,
        query_target_names=q_names,
        target_ctrl_mean_libnorm=mu_t,
        fc_bundle=fc_bundle,
        delta_obs=obs,
        delta_pred_fct=pred,
        delta_transfer=t_mean,
        target_gene_idx=q_gidx,
        n_cells=nc_vec,
        knockdown_depth=kd_vec,
        delta_half_a=obs,
        delta_half_b=obs,
        mismatch_preds={"NegCtrl_FCT_MismatchedMean(S0)": t_mean},
        n_boot=200,
        seed=20261122,
    )
    assert s7["spearman_r_orth_p_vs_ln_ncells"]["estimate"] > 0.5
    ols_b = s7["ols_adjusted_both"]
    assert abs(ols_b["adjusted_intercept"]["estimate"] - ols_b["unadjusted_mean_r_orth_p"]["estimate"]) < 1e-10
    assert list(s7["n_cells_tertiles"].keys()) == ["T1_low", "T2_mid", "T3_high"]
    assert s7["n_cells_tertiles"]["T3_high"]["n_cells_min"] > s7["n_cells_tertiles"]["T1_low"]["n_cells_max"]

    # R19 cluster-df W^2 bias test: pure Tier-1 signal + i.i.d. noise -> L2 stays small near 0
    X_tm = rng.normal(0.0, 0.2, size=(60, 80))
    Y_t1 = 0.05 + 1.3 * X_tm + rng.normal(0.0, 0.15, size=(60, 80))
    r19 = compute_r19_variance_decomposition_r10(
        Y_t1, X_tm, np.arange(60, dtype=np.int64), n_boot=100, seed=20261122
    )
    l1 = r19["all_nontarget"]["L1_whole_line_frac"]["estimate"]
    l2 = r19["all_nontarget"]["L2_gene_gain_frac"]["estimate"]
    l3 = r19["all_nontarget"]["L3_remaining_resid_frac"]["estimate"]
    assert abs((l1 + l2 + l3) - 1.0) < 1e-9
    assert l2 < 0.05


if __name__ == "__main__":
    test_hgnc_bridge_and_gate()
    test_hdf5_stream_loaders_and_splits()
    test_fct_parity_with_round9_and_adjudication()
    test_s7_cellcount_kd_and_r19_unbiased()
    print("R10_ALL_UNIT_TESTS_PASSED")
