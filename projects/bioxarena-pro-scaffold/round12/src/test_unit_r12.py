#!/usr/bin/env python3
"""
Unit & Integration Tests for Round 12 (`PREREG_E`) Modules (`round12/src/`).

Runs cleanly on synthetic fixtures and cross-checks exact numerical parity against
`round11/audit/exploratory_check/` and frozen Round 10/11 primitives without accessing
any unseen Round 12 held-out data.
"""

from __future__ import annotations

import gzip
import sys
import tempfile
import unittest
from pathlib import Path
from typing import List

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

R12_SRC_DIR = Path(__file__).resolve().parent
REPO_ROOT = R12_SRC_DIR.parents[1]
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
R10_SCRIPTS_DIR = REPO_ROOT / "round10" / "shared" / "scripts"
R11_SRC_DIR = REPO_ROOT / "round11" / "src"
R11_AUDIT_DIR = REPO_ROOT / "round11" / "audit"
R11_EXPL_DIR = R11_AUDIT_DIR / "exploratory_check"
for p in (str(R10_SRC_DIR), str(R10_SCRIPTS_DIR), str(R11_SRC_DIR), str(R11_AUDIT_DIR), str(R11_EXPL_DIR), str(R12_SRC_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from exploratory_core import evaluate_cross_mask_split_minus_all  # noqa: E402
from run_n5_audit import (  # noqa: E402
    _eval_n5_from_moments as audit_eval_n5_from_moments,
    _pred_n5_moments as audit_pred_n5_moments,
    prepare_n5_per_pert_data as audit_prepare_n5,
)

from data_loader_r10 import LineDatasetR10, load_standardized_line_r10  # noqa: E402
from fct_pipeline import FctSourceState, SourceLineInputR10, fct_from_fc_bar  # noqa: E402
from prepare_r10_lines import HGNCLookup, build_hgnc_lookup, process_r10_line  # noqa: E402
from predictors_r11 import predict_group_d_suite  # noqa: E402
from test_unit_r11 import _make_planted_synthetic_cohort  # noqa: E402

from adjudicate_r12 import adjudicate_r12_results, render_r12_markdown_report  # noqa: E402
from data_loader_r12 import check_prereg_e_gate, parse_confirmatory_lines_r12  # noqa: E402
from dataset_readers_r12 import (  # noqa: E402
    build_a549_parts,
    build_hek293_parts,
    check_hek293_feasibility,
)
from metrics_r12 import (  # noqa: E402
    compute_joint_sa8_split_minus_all_r12,
    compute_n5_perm100_r12,
    compute_n5_suite_r12,
    compute_paired_fct_minus_all_mismatched_r12,
)
from run_eval_r12 import evaluate_target_against_sources_r12  # noqa: E402

_CACHED_HGNC: HGNCLookup | None = None


def _get_hgnc() -> HGNCLookup:
    global _CACHED_HGNC
    if _CACHED_HGNC is None:
        _CACHED_HGNC = build_hgnc_lookup()
    return _CACHED_HGNC


def _pick_hgnc_genes(n_genes: int) -> List[str]:
    hgnc = _get_hgnc()
    picked = [g for g in sorted(hgnc.approved) if g and not g.startswith("MT-")][:n_genes]
    if len(picked) < n_genes:
        raise RuntimeError(f"Requested {n_genes} HGNC genes, found {len(picked)}")
    return picked


class Round12UnitTests(unittest.TestCase):
    def test_1_blinding_gate_and_confirmatory_file_r12(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)
            with self.assertRaises(RuntimeError):
                check_prereg_e_gate(p)
            (p / "PREREG_E_COMMIT.txt").write_text("commit_e_111\n", encoding="utf-8")
            with self.assertRaises(RuntimeError):
                check_prereg_e_gate(p)
            (p / "CODE_FREEZE_R12.txt").write_text("freeze_e_222\n", encoding="utf-8")
            pr, fr = check_prereg_e_gate(p)
            self.assertEqual(pr, "commit_e_111")
            self.assertEqual(fr, "freeze_e_222")

            cf = p / "CONFIRMATORY_LINES_R12.txt"
            cf.write_text("# comment\n\na549\n# hek293: SECONDARY_ONLY\n", encoding="utf-8")
            self.assertEqual(parse_confirmatory_lines_r12(cf), ["a549"])

    def test_2_sa8_joint_split_minus_all_parity_vs_audit_and_paired_mis(self) -> None:
        src_map, tgt_ds = _make_planted_synthetic_cohort(seed=121, n_genes=150, n_targets=55, n_sources=6)
        # Drop 5 targets from half_valid to test exact common-perturbation rows_common handling
        hv = tgt_ds.half_valid.copy()
        hv[:5] = False
        tgt_part = LineDatasetR10(**{**tgt_ds.__dict__, "half_valid": hv})

        src_inputs = [
            SourceLineInputR10(
                name=k,
                target_names=v.target_names,
                target_indices=v.target_indices,
                deltas=v.deltas,
                ctrl_cov=None,
                ctrl_mean_libnorm=v.ctrl_mean_libnorm,
            )
            for k, v in sorted(src_map.items())
        ]
        state = FctSourceState.fit(src_inputs, mask_self=True)
        preds, _, _ = predict_group_d_suite(
            state=state,
            target_ctrl_mean_libnorm=tgt_part.ctrl_mean_libnorm,
            query_target_names=tgt_part.target_names,
            include_r9_anchors=False,
        )

        # 1) Joint SA-8 split - all vs independent audit implementation in round11/audit/exploratory_check
        seed_test = 20261201
        r12_diff = compute_joint_sa8_split_minus_all_r12(
            delta_obs=tgt_part.deltas,
            delta_pred=preds["E1_ConstMed"],
            delta_transfer=preds["Transfer_Mean"],
            target_gene_idx=tgt_part.target_indices,
            delta_half_a=tgt_part.half_a_deltas,
            delta_half_b=tgt_part.half_b_deltas,
            half_valid=tgt_part.half_valid,
            top_k=50,
            n_boot=250,
            seed=seed_test,
        )
        audit_diff = evaluate_cross_mask_split_minus_all(
            delta_obs=tgt_part.deltas,
            pred=preds["E1_ConstMed"],
            t_mean=preds["Transfer_Mean"],
            target_cols=tgt_part.target_indices,
            half_a_obs=tgt_part.half_a_deltas,
            half_b_obs=tgt_part.half_b_deltas,
            half_valid=tgt_part.half_valid,
            top_k=50,
            n_boot=250,
            seed=seed_test,
        )
        self.assertEqual(r12_diff["n_common_perts"], 50)
        self.assertEqual(r12_diff["n_common_perts"], audit_diff["n_eval_perts"])
        self.assertAlmostEqual(r12_diff["estimate"], audit_diff["diff_est"], places=12)
        self.assertAlmostEqual(r12_diff["ci95_lo"], audit_diff["diff_ci95_lo"], places=12)
        self.assertAlmostEqual(r12_diff["ci975_lo"], audit_diff["diff_ci975_lo"], places=12)
        self.assertAlmostEqual(r12_diff["ci975_hi"], audit_diff["diff_ci975_hi"], places=12)

        # 2) Paired FCT - Mis_s across ALL 6 mismatched sources
        paired = compute_paired_fct_minus_all_mismatched_r12(
            delta_obs=tgt_part.deltas,
            preds=preds,
            delta_transfer=preds["Transfer_Mean"],
            target_gene_idx=tgt_part.target_indices,
            delta_half_a=tgt_part.half_a_deltas,
            delta_half_b=tgt_part.half_b_deltas,
            half_valid=tgt_part.half_valid,
            n_boot=120,
            seed=seed_test,
        )
        self.assertEqual(sorted(paired.keys()), sorted(src_map.keys()))
        for _, res in paired.items():
            self.assertGreater(res["ci975_lo"], 0.30)

    def test_3_n5_null_cone_parity_vs_audit_and_construction_roles(self) -> None:
        src_map, tgt_ds = _make_planted_synthetic_cohort(seed=131, n_genes=140, n_targets=50, n_sources=4)
        eye_cov = np.eye(len(tgt_ds.genes), dtype=np.float64) * 0.05
        src_inputs = [
            SourceLineInputR10(
                name=k,
                target_names=v.target_names,
                target_indices=v.target_indices,
                deltas=v.deltas,
                ctrl_cov=eye_cov,
                ctrl_mean_libnorm=v.ctrl_mean_libnorm,
            )
            for k, v in sorted(src_map.items())
        ]
        state = FctSourceState.fit(src_inputs, mask_self=True)
        preds, fc_bundle, _ = predict_group_d_suite(
            state=state,
            target_ctrl_mean_libnorm=tgt_ds.ctrl_mean_libnorm,
            query_target_names=tgt_ds.target_names,
            query_target_gene_indices=tgt_ds.target_indices,
            include_r9_anchors=True,
            include_exploratory_e=True,
            target_ctrl_cov=eye_cov,
        )
        _, fc_bar, _ = fc_bundle
        seed_test = 20261201
        n5_r12 = compute_n5_suite_r12(
            delta_obs=tgt_ds.deltas,
            preds=preds,
            delta_transfer=preds["Transfer_Mean"],
            grand_mean=state.grand_mean,
            fc_bar=fc_bar,
            target_gene_idx=tgt_ds.target_indices,
            delta_half_a=tgt_ds.half_a_deltas,
            delta_half_b=tgt_ds.half_b_deltas,
            half_valid=tgt_ds.half_valid,
            n_boot=150,
            seed=seed_test,
        )
        aud_perp, aud_xx, aud_ha, aud_hb, aud_hab, _ = audit_prepare_n5(
            delta_obs=tgt_ds.deltas,
            t_mean=preds["Transfer_Mean"],
            grand_mean=state.grand_mean,
            fc_bar=fc_bar,
            target_cols=tgt_ds.target_indices,
            half_a_obs=tgt_ds.half_a_deltas,
            half_b_obs=tgt_ds.half_b_deltas,
            half_valid=tgt_ds.half_valid,
            mode="all_nontarget",
        )
        m_aud = len(aud_perp)
        aud_draws = np.random.default_rng(seed_test).integers(0, m_aud, size=(150, m_aud))
        for mn in ("FCT", "Source_Grand_Mean", "Transfer_Mean_GlobalScale", "E1_ConstMed", "E6_LogFCBar", "NegCtrl_FCT_MismatchedMean(src_0)"):
            r12_m = n5_r12["methods"][mn]
            xy, yy_o, yy_c = audit_pred_n5_moments(aud_perp, tgt_ds.deltas, preds[mn], preds["Transfer_Mean"])
            aud_m = audit_eval_n5_from_moments(xy, yy_o, yy_c, aud_xx, aud_ha, aud_hb, aud_hab, aud_draws, 150)
            self.assertAlmostEqual(r12_m["r_cone_sp"]["estimate"], aud_m["r_cone_sp_est"], places=12)
            self.assertAlmostEqual(r12_m["r_cone_sp"]["ci975_lo"], aud_m["r_cone_sp_ci975_lo"], places=12)
            self.assertAlmostEqual(r12_m["eta"], aud_m["eta_est"], places=12)
            self.assertAlmostEqual(r12_m["r2_n5"], aud_m["r2_in_n5"], places=12)

        # Verify construction-zero & cone-suppressed roles explicitly
        for z_m in ("Transfer_Mean", "Source_Grand_Mean", "Transfer_Mean_GlobalScale"):
            zm = n5_r12["methods"][z_m]
            self.assertEqual(zm["construction_status"], "ZERO_BY_EXACT_CONSTRUCTION")
            self.assertLess(zm["eta"], 1e-12)
            self.assertEqual(zm["r_cone_sp"]["estimate"], 0.0)
        for s_m in ("E1_ConstMed", "E2_ConstP10", "E2_ConstP90", "E6_LogFCBar", "E7_FCBarMinus1"):
            sm = n5_r12["methods"][s_m]
            self.assertEqual(sm["construction_status"], "SUPPRESSED_BY_CONE_CONSTRUCTION")
            self.assertGreaterEqual(sm["r2_n5"], 0.99)
        self.assertEqual(n5_r12["methods"]["FCT"]["construction_status"], "INFORMATIVE")
        self.assertGreater(n5_r12["methods"]["FCT"]["r_cone_sp"]["ci975_lo"], 0.40)

        # Verify N_5 perm100 matches audit perm100 single-draw logic
        p_r12 = compute_n5_perm100_r12(
            state=state,
            target_ctrl_mean_libnorm=tgt_ds.ctrl_mean_libnorm,
            query_target_names=tgt_ds.target_names,
            delta_obs=tgt_ds.deltas,
            delta_transfer=preds["Transfer_Mean"],
            grand_mean=state.grand_mean,
            target_gene_idx=tgt_ds.target_indices,
            delta_half_a=tgt_ds.half_a_deltas,
            delta_half_b=tgt_ds.half_b_deltas,
            half_valid=tgt_ds.half_valid,
            fc_bundle=fc_bundle,
            n_perm=10,
        )
        n_g = tgt_ds.deltas.shape[1]
        p_aud_cal, p_aud_eta = [], []
        for k in range(10):
            perm = np.random.default_rng(k).permutation(n_g)
            pred_pk = fct_from_fc_bar(tgt_ds.ctrl_mean_libnorm[perm], fc_bundle[1], fc_bundle[0], fc_bundle[2])
            xy, yy_o, yy_c = audit_pred_n5_moments(aud_perp, tgt_ds.deltas, pred_pk, preds["Transfer_Mean"])
            r_k = audit_eval_n5_from_moments(xy, yy_o, yy_c, aud_xx, aud_ha, aud_hb, aud_hab, draws=None, n_boot=1)
            p_aud_cal.append(r_k["r_cone_sp_est"])
            p_aud_eta.append(r_k["eta_est"])
        self.assertAlmostEqual(p_r12["r_cone_sp"]["p97_5"], float(np.percentile(p_aud_cal, 97.5)), places=12)
        self.assertAlmostEqual(p_r12["med_eta"], float(np.median(p_aud_eta)), places=12)
        self.assertGreater(p_r12["med_eta"], 0.30)

    def test_4_adjudicator_r12_k0_k1_k2_and_iut_guards(self) -> None:
        src6, tgt1 = _make_planted_synthetic_cohort(seed=41, n_genes=150, n_targets=60, n_sources=6)
        tgt1_named = LineDatasetR10(**{**tgt1.__dict__, "name": "a549"})

        rng2 = np.random.default_rng(42)
        n_t, n_g = tgt1.deltas.shape
        mu_t2 = np.exp(rng2.uniform(-1.5, 3.0, size=n_g))
        mu_s0 = src6["src_0"].ctrl_mean_libnorm
        fc0 = np.maximum(np.expm1(np.log1p(mu_s0)[None, :] + src6["src_0"].deltas), 1e-8) / mu_s0[None, :]
        d_clean2 = np.log1p(np.maximum(mu_t2[None, :] * fc0, 1e-8)) - np.log1p(mu_t2[None, :])
        ha2 = d_clean2 + rng2.normal(0.0, 0.05, size=d_clean2.shape)
        hb2 = d_clean2 + rng2.normal(0.0, 0.05, size=d_clean2.shape)
        tgt2 = LineDatasetR10(
            **{
                **tgt1.__dict__,
                "name": "hek293",
                "deltas": 0.5 * (ha2 + hb2),
                "ctrl_mean_log1p": np.log1p(mu_t2),
                "ctrl_mean_libnorm": mu_t2,
                "ctrl_mean_libnorm_full": mu_t2,
                "half_a_deltas": ha2,
                "half_b_deltas": hb2,
            }
        )

        l1, b1, ex1 = evaluate_target_against_sources_r12(
            tgt1_named, src6, include_r9_anchors=False, min_gene_overlap=100, n_boot=200, seed=20261201, n_perm=15
        )
        l2, b2, ex2 = evaluate_target_against_sources_r12(
            tgt2, src6, include_r9_anchors=False, min_gene_overlap=100, n_boot=200, seed=20261201, n_perm=15
        )
        df_l = pd.DataFrame(l1 + l2)
        df_b = pd.DataFrame(b1 + b2)
        summary = {
            "phase": "holdout_r12",
            "source_mode": "r9_6",
            "source_lines": sorted(src6.keys()),
            "prereg_commit_hash": "prereg_e",
            "code_freeze_hash": "freeze_r12",
            "bootstrap_seed": 20261201,
            "n_boot": 200,
            "fct_perm100": {"a549": ex1["perm100"], "hek293": ex2["perm100"]},
            "n5_cone_metrics": {"a549": ex1["n5_cone_metrics"], "hek293": ex2["n5_cone_metrics"]},
            "n5_perm100": {"a549": ex1["n5_perm100"], "hek293": ex2["n5_perm100"]},
        }

        # Case 0: K_conf = 0 aborts confirmatory verdicts
        adj0 = adjudicate_r12_results(df_l, df_b, summary, confirmatory_lines=[])
        self.assertEqual(adj0["confirmatory_status"], "NO_CONFIRMATORY_LINES")
        self.assertEqual(adj0["Group_E_Claim1_Verdict"], "NOT_EVALUATED_NO_CONFIRMATORY_LINES")
        self.assertEqual(adj0["Group_E_Claim2_Verdict"], "NOT_EVALUATED_NO_CONFIRMATORY_LINES")

        # Case 1: K_conf = 1 (a549 only, hek293 secondary) -> IS confirmatory + limited generalizability flag
        adj1 = adjudicate_r12_results(df_l, df_b, summary, confirmatory_lines=["a549"])
        self.assertEqual(adj1["confirmatory_status"], "CONFIRMATORY_SINGLE_LINE_K1")
        self.assertTrue(adj1["limited_single_line_generalizability_notice"])
        self.assertEqual(adj1["secondary_lines"], ["hek293"])
        self.assertEqual(adj1["Group_E_Claim1_Verdict"], "CONFIRMED")
        self.assertEqual(adj1["Group_E_Claim2_Verdict"], "CONFIRMED")
        self.assertEqual(adj1["H-E1"]["per_line"]["a549"]["H-E1b_n_sources_checked"], 6)
        md1 = render_r12_markdown_report(df_l, df_b, summary, adj1, "unit_k1")
        self.assertIn("single-line generalizability is limited", md1)
        self.assertIn("ZERO_BY_EXACT_CONSTRUCTION", md1)

        # Case 2: K_conf = 2 passes both Claim 1 & Claim 2 unanimously
        adj2 = adjudicate_r12_results(df_l, df_b, summary, confirmatory_lines=["a549", "hek293"])
        self.assertEqual(adj2["confirmatory_status"], "CONFIRMATORY_K2")
        self.assertFalse(adj2["limited_single_line_generalizability_notice"])
        self.assertEqual(adj2["Group_E_Claim1_Verdict"], "CONFIRMED")
        self.assertEqual(adj2["Group_E_Claim2_Verdict"], "CONFIRMED")

        # Case 3: Fail JUST 1/6 paired FCT - Mis_s on hek293 -> Claim 1 MUST be 未确认
        df_b_fail_1src = df_b.copy()
        m_1src = (df_b_fail_1src["target_line"] == "hek293") & (
            df_b_fail_1src["comparison"] == "H-E1b_FCT_minus_Mis_src_0_all_nontarget"
        )
        df_b_fail_1src.loc[m_1src, "ci975_lo"] = -0.001
        adj3 = adjudicate_r12_results(df_l, df_b_fail_1src, summary, confirmatory_lines=["a549", "hek293"])
        self.assertEqual(adj3["H-E1"]["sub_hypotheses"]["H-E1b"]["verdict"], "未确认")
        self.assertEqual(adj3["Group_E_Claim1_Verdict"], "未确认")
        self.assertEqual(adj3["Group_E_Claim2_Verdict"], "CONFIRMED")

        # Case 4: Fail Addendum A1 valid_frac < 0.95 on a549 -> Claim 1 MUST be 未确认
        df_l_a1 = df_l.copy()
        df_l_a1.loc[(df_l_a1["target_line"] == "a549") & (df_l_a1["method"] == "FCT"), "r_orth_cal_valid_frac"] = 0.91
        adj4 = adjudicate_r12_results(df_l_a1, df_b, summary, confirmatory_lines=["a549"])
        self.assertTrue(adj4["H-E1"]["per_line"]["a549"]["addendum_a1_ci_gt0_but_low_valid_frac"])
        self.assertEqual(adj4["Group_E_Claim1_Verdict"], "未确认")

        # Case 5: Fail H-E2b on a549 -> Claim 2 MUST be 未确认 while Claim 1 stays CONFIRMED
        df_b_fail_h2b = df_b.copy()
        m_h2b = (df_b_fail_h2b["target_line"] == "a549") & (df_b_fail_h2b["comparison"] == "H-E2b_E1_split_minus_all")
        df_b_fail_h2b.loc[m_h2b, "ci975_lo"] = -0.002
        adj5 = adjudicate_r12_results(df_l, df_b_fail_h2b, summary, confirmatory_lines=["a549"])
        self.assertEqual(adj5["Group_E_Claim1_Verdict"], "CONFIRMED")
        self.assertEqual(adj5["Group_E_Claim2_Verdict"], "未确认")

    def test_5_dataset_reader_a549_synthetic(self) -> None:
        rng = np.random.default_rng(1201)
        hgnc = _get_hgnc()
        genes = _pick_hgnc_genes(520)
        gene_syms = genes + ["MT-CO1", "CAS9"]
        n_genes_all = len(gene_syms)

        # Build 4 batches (1,3 -> Half A; 2,4 -> Half B) + ~27% CAS9==0 cells + multi-guide exclusions
        cell_rows = []
        cell_kinds = []  # (kind, tgt_idx_or_-1, cas9_cnt)
        batches = ["batch1", "batch2", "batch3", "batch4"]
        cid = 0
        for b in batches:
            for i in range(18):
                c9 = 0 if (i % 4 == 0) else 3
                cell_rows.append({"cell_id": f"C{cid}", "batch": b, "lane": f"{b}_L1", "sgRNA": f"NTC_sg{i % 2}", "Perturbation": "Nontargeting", "total_counts_CAS9": c9})
                cell_kinds.append(("ntc", -1, c9))
                cid += 1
            for ti, gsym in enumerate((genes[0], genes[1])):
                for i in range(16):
                    c9 = 0 if (i % 4 == 0) else 2
                    cell_rows.append({"cell_id": f"C{cid}", "batch": b, "lane": f"{b}_L1", "sgRNA": f"{gsym}_sg{i % 2}", "Perturbation": gsym, "total_counts_CAS9": c9})
                    cell_kinds.append(("tgt", ti, c9))
                    cid += 1
        # Excluded rows (multi-guide & unassigned)
        cell_rows.append({"cell_id": f"C{cid}", "batch": "batch1", "lane": "batch1_L1", "sgRNA": "g1,g2", "Perturbation": f"{genes[0]},{genes[1]}", "total_counts_CAS9": 3})
        cell_kinds.append(("excl", -1, 3))
        cell_rows.append({"cell_id": f"C{cid + 1}", "batch": "batch2", "lane": "batch2_L1", "sgRNA": "unassigned", "Perturbation": "unassigned", "total_counts_CAS9": 3})
        cell_kinds.append(("excl", -1, 3))

        n_cells = len(cell_rows)
        X = np.zeros((n_cells, n_genes_all), dtype=np.int32)
        X[:, : len(genes)] = rng.integers(1, 5, size=(n_cells, len(genes)), dtype=np.int32)
        X[:, len(genes)] = 5  # low MT-CO1
        for ci, (kind, ti, c9) in enumerate(cell_kinds):
            X[ci, len(genes) + 1] = c9
            if kind == "tgt":
                X[ci, :25] += 8
                X[ci, ti] = 1
        X_coo = sp.coo_matrix(X)

        def _write_mtx(path: Path, coo: sp.coo_matrix, sep: str) -> None:
            with gzip.open(path, "wt", encoding="utf-8") as f:
                f.write("%%MatrixMarket matrix coordinate integer general\n%\n")
                f.write(f"{coo.shape[0]} {coo.shape[1]} {coo.nnz}\n")
                for r, c, v in zip(coo.row + 1, coo.col + 1, coo.data):
                    f.write(f"{r}{sep}{c}{sep}{int(v)}\n")

        with tempfile.TemporaryDirectory() as td:
            raw_dir = Path(td) / "raw"
            adir = raw_dir / "a549_ko"
            adir.mkdir(parents=True)
            pd.DataFrame(cell_rows).to_csv(adir / "GSE345058_cell_annotations.csv.gz", index=False, compression="gzip")
            pd.DataFrame(
                [
                    {
                        "gene_symbol": ("Cas9" if gs == "CAS9" else gs),
                        "ensembl_id": ("Cas9" if gs == "CAS9" else hgnc.sym_to_ens.get(gs, gs)),
                        "genome": ("CAS9" if gs == "CAS9" else "GRCh38"),
                        "feature_types": "Gene Expression",
                    }
                    for gs in gene_syms
                ]
            ).to_csv(adir / "GSE345058_gene_annotations.csv.gz", index=False, compression="gzip")
            pd.DataFrame(
                [
                    {"sgRNA": "NTC_sg0", "Perturbation": "Nontargeting"},
                    {"sgRNA": f"{genes[0]}_sg0", "Perturbation": genes[0]},
                ]
            ).to_csv(adir / "GSE345058_sgRNA.csv.gz", index=False, compression="gzip")
            _write_mtx(adir / "GSE345058_raw_counts.mtx.gz", X_coo, sep=" ")

            # 1) Full A549 (cells x genes MTX) -> Half A=batch1+3, Half B=batch2+4 -> standardize
            parts, audit = build_a549_parts(raw_dir, hgnc, exclude_zero_cas9=False, bucket_size=90, chunk_lines=25000)
            self.assertEqual(audit["n_multiguide_or_unassigned_excluded"], 2)
            self.assertEqual(audit["n_qc_kept_cells"], (18 + 32) * 4)
            self.assertEqual(set(np.unique(parts[0].half_flags).tolist()), {1, 2})

            out_root = Path(td) / "lines"
            process_r10_line("a549_ko", parts, hgnc, out_root=out_root, extra_meta={"adapter_audit": audit}, n_workers=1)
            ds = load_standardized_line_r10(out_root / "a549_ko", min_gene_overlap=50)
            self.assertEqual(len(ds.target_names), 2)
            self.assertTrue(np.any(ds.repr_mask))

            # 2) Sensitivity exclude_zero_cas9=True drops CAS9==0 subset
            parts_c9, audit_c9 = build_a549_parts(raw_dir, hgnc, exclude_zero_cas9=True, bucket_size=90)
            self.assertGreater(audit_c9["n_cas9_zero_excluded"], 0)
            self.assertLess(audit_c9["n_qc_kept_cells"], audit["n_qc_kept_cells"])

            # 3) Transposed (genes x cells) tab-delimited MTX + feature CAS9 fallback without metadata column
            cell_df_noc9 = pd.DataFrame(cell_rows).drop(columns=["total_counts_CAS9"])
            cell_df_noc9.to_csv(adir / "GSE345058_cell_annotations.csv.gz", index=False, compression="gzip")
            _write_mtx(adir / "GSE345058_raw_counts.mtx.gz", X_coo.transpose().tocoo(), sep="\t")
            parts_tr, audit_tr = build_a549_parts(raw_dir, hgnc, exclude_zero_cas9=True, bucket_size=90)
            self.assertEqual(audit_tr["n_qc_kept_cells"], audit_c9["n_qc_kept_cells"])
            self.assertEqual(int(parts_tr[0].nt_mask.sum()), int(parts_c9[0].nt_mask.sum()))

    def test_6_dataset_reader_hek293_synthetic(self) -> None:
        rng = np.random.default_rng(1301)
        hgnc = _get_hgnc()
        genes = _pick_hgnc_genes(520)
        var_names = genes + ["MT-CO1", "dCas9-KRAB-MeCP2"]
        n_var = len(var_names)

        # Build 4 replicates ('repA'..'repD', 1-based odd->Half A, even->Half B) + low-KD targets + multi-sgRNA exclusions
        obs_records = []
        for rep_lbl in ("repA", "repB", "repC", "repD"):
            for i in range(18):
                obs_records.append({"pert": "control", "nperts": 0, "rep": rep_lbl, "kd_pct": 0.0, "tgt_i": -1})
            for ti, gsym in enumerate((genes[0], genes[1])):
                for i in range(16):
                    kd = 0.25 if (i % 2 == 0) else 0.80
                    obs_records.append({"pert": gsym, "nperts": 1, "rep": rep_lbl, "kd_pct": kd, "tgt_i": ti})
        # Excluded doublets / unassigned
        obs_records.append({"pert": f"{genes[0]},{genes[1]}", "nperts": 2, "rep": "repA", "kd_pct": 0.9, "tgt_i": -1})
        obs_records.append({"pert": "unassigned", "nperts": 1, "rep": "repB", "kd_pct": 0.9, "tgt_i": -1})

        n_cells = len(obs_records)
        X = rng.integers(1, 5, size=(n_cells, n_var), dtype=np.int32)
        X[:, -2] = 5  # low MT-CO1
        X[:, -1] = 2  # synthetic dCas9 vector row (MUST be excluded)
        for ci, rec in enumerate(obs_records):
            if rec["tgt_i"] >= 0:
                X[ci, :25] += 8
                X[ci, rec["tgt_i"]] = 1
        X_csr = sp.csr_matrix(X, dtype=np.float32)

        def _write_hek_h5ad(path: Path, mat: sp.spmatrix, enc: str, rep_col_vals: List[str]) -> None:
            with h5py.File(path, "w") as f:
                xg = f.create_group("X")
                xg.attrs["encoding-type"] = enc
                xg.attrs["shape"] = np.array(mat.shape, dtype=np.int64)
                xg.create_dataset("data", data=mat.data.astype(np.float32))
                xg.create_dataset("indices", data=mat.indices.astype(np.int32))
                xg.create_dataset("indptr", data=mat.indptr.astype(np.int64))
                og = f.create_group("obs")
                og.create_dataset("_index", data=np.array([f"cell_{i}" for i in range(n_cells)], dtype="S"))
                og.create_dataset("perturbation", data=np.array([r["pert"] for r in obs_records], dtype="S"))
                og.create_dataset("nperts", data=np.array([r["nperts"] for r in obs_records], dtype=np.int32))
                og.create_dataset("replicate", data=np.array(rep_col_vals, dtype="S"))
                og.create_dataset("kd_pct", data=np.array([r["kd_pct"] for r in obs_records], dtype=np.float32))
                vg = f.create_group("var")
                vg.create_dataset("_index", data=np.array(var_names, dtype="S"))
                vg.create_dataset("gene_symbol", data=np.array(var_names, dtype="S"))
                vg.create_dataset("feature_types", data=np.array(["Gene Expression"] * n_var, dtype="S"))

        with tempfile.TemporaryDirectory() as td:
            raw_dir = Path(td) / "raw"
            hdir = raw_dir / "hek293_xucao"
            hdir.mkdir(parents=True)
            _write_hek_h5ad(hdir / "XuCao2023.h5ad", X_csr, "csr_matrix", [r["rep"] for r in obs_records])

            feas = check_hek293_feasibility(raw_dir, hgnc=hgnc)
            self.assertTrue(feas["feasible"])
            self.assertEqual(feas["int_count_sources"], ["X"])
            self.assertEqual(feas["n_unique_canon_targets_meta"], 2)

            # Default: ALL single-sgRNA cells (apply_author_kd40_filter=False), obs replicate parity split
            parts, audit = build_hek293_parts(raw_dir, hgnc, apply_author_kd40_filter=False)
            self.assertFalse(audit["used_seed42_fallback"])
            self.assertTrue(audit["half_rule_used"].startswith("obs_replicate_parity"))
            self.assertEqual(parts[0].n_genes_raw, len(genes) + 1)  # dCas9 excluded
            self.assertEqual(audit["n_multiguide_or_unassigned_excluded"], 2)
            self.assertEqual(audit["n_qc_kept_cells"], (18 + 32) * 4)

            out_root = Path(td) / "lines"
            process_r10_line("hek293_xucao", parts, hgnc, out_root=out_root, extra_meta={"adapter_audit": audit}, n_workers=1)
            ds = load_standardized_line_r10(out_root / "hek293_xucao", min_gene_overlap=50)
            self.assertEqual(len(ds.target_names), 2)

            # Sensitivity: apply_author_kd40_filter=True drops <=40% KD target cells while preserving NTCs
            parts_kd, audit_kd = build_hek293_parts(raw_dir, hgnc, apply_author_kd40_filter=True)
            self.assertEqual(audit_kd["n_kd40_excluded"], 4 * 2 * 8)
            self.assertEqual(audit_kd["n_kept_ntc"], audit["n_kept_ntc"])

            # CSC `.h5ad` + single replicate triggers deterministic seed=42 stratified 50/50 fallback
            _write_hek_h5ad(hdir / "XuCao2023.h5ad", X_csr.tocsc(), "csc_matrix", ["batch0"] * n_cells)
            parts_bkp, audit_bkp = build_hek293_parts(raw_dir, hgnc)
            self.assertTrue(audit_bkp["used_seed42_fallback"])
            self.assertEqual(audit_bkp["half_rule_used"], "seed42_stratified_50_50")
            self.assertEqual(audit_bkp["n_kept_half_a_ntc"], audit_bkp["n_kept_half_b_ntc"])

            # E3 float-only X fails feasibility; layers/counts passes and loads identically
            X_log_float = sp.csr_matrix(np.log1p(X).astype(np.float32))
            _write_hek_h5ad(hdir / "XuCao2023.h5ad", X_log_float, "csr_matrix", [r["rep"] for r in obs_records])
            self.assertFalse(check_hek293_feasibility(raw_dir, hgnc=hgnc)["feasible"])
            with h5py.File(hdir / "XuCao2023.h5ad", "a") as fh:
                lg = fh.create_group("layers").create_group("counts")
                lg.attrs["encoding-type"] = "csr_matrix"
                lg.attrs["shape"] = np.array(X_csr.shape, dtype=np.int64)
                lg.create_dataset("data", data=X_csr.data.astype(np.float32))
                lg.create_dataset("indices", data=X_csr.indices.astype(np.int32))
                lg.create_dataset("indptr", data=X_csr.indptr.astype(np.int64))
            feas_lay = check_hek293_feasibility(raw_dir, hgnc=hgnc)
            self.assertTrue(feas_lay["feasible"])
            self.assertEqual(feas_lay["int_count_sources"], ["layers/counts"])
            parts_lay, _ = build_hek293_parts(raw_dir, hgnc, x_layer="counts")
            self.assertTrue(np.allclose(parts_lay[0].umi, parts[0].umi))


if __name__ == "__main__":
    unittest.main(verbosity=2)
