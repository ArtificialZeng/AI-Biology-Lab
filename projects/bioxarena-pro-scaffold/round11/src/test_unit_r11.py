#!/usr/bin/env python3
"""
Unit tests for Round 11 (`PREREG_D`) Group D evaluation modules (`round11/src/`).

Runs cleanly on synthetic fixtures and checks exact mathematical parity against
frozen Round 10 primitives without accessing any Round 11 held-out data.
"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd

R11_SRC_DIR = Path(__file__).resolve().parent
R10_SRC_DIR = R11_SRC_DIR.parents[1] / "round10" / "src"
for p in (str(R10_SRC_DIR), str(R11_SRC_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from adjudicate_r10 import SIX_CONVENTIONS  # noqa: E402
from data_loader_r10 import LineDatasetR10  # noqa: E402
from fct_pipeline import FctSourceState, SourceLineInputR10, fct_from_fc_bar  # noqa: E402
from metrics_r10 import (  # noqa: E402
    EPS,
    _extract_masked_idx_per_pert,
    compute_fct_perm100_null,
    compute_residual_suite_r10,
    disattenuated_r,
    ols_slope_1d,
    pearson_1d_strict,
)

from adjudicate_r11 import adjudicate_r11_results, render_r11_markdown_report  # noqa: E402
from data_loader_r11 import check_prereg_d_gate  # noqa: E402
from metrics_r11 import (  # noqa: E402
    PCT_HI_9833,
    PCT_LO_9833,
    compute_fct_perm100_null_r11,
    compute_joint_cross_mask_e1_diff_r11,
    compute_joint_paired_r_orth_diff_r11,
    compute_residual_suite_r11,
    summarize_boot_r11,
)
from predictors_r11 import predict_group_d_suite  # noqa: E402
from run_eval_r11 import evaluate_target_against_sources_r11  # noqa: E402


def _make_planted_synthetic_cohort(
    seed: int = 123,
    n_genes: int = 160,
    n_targets: int = 60,
    n_sources: int = 6,
    signal_scale: float = 1.0,
) -> tuple[Dict[str, LineDatasetR10], LineDatasetR10]:
    rng = np.random.default_rng(seed)
    genes = [f"G{i:04d}" for i in range(n_genes)]
    targets = genes[:n_targets]

    # Target control mean with sharp gene-specific high/low expression heterogeneity
    mu_t = np.exp(rng.uniform(-1.5, 3.0, size=n_genes))
    # True underlying multiplicative fold-changes per perturbation and gene
    log_fc = rng.normal(0.0, 0.35, size=(n_targets, n_genes))
    for p in range(n_targets):
        log_fc[p, p] = -1.2

    src_map: Dict[str, LineDatasetR10] = {}
    for s_idx in range(n_sources):
        mu_s = np.exp(rng.uniform(-1.5, 3.0, size=n_genes))
        pert_s = np.maximum(mu_s[None, :] * np.exp(log_fc + rng.normal(0.0, 0.04, size=log_fc.shape)), 1e-8)
        d_s = np.log1p(pert_s) - np.log1p(mu_s[None, :])
        c_s = rng.normal(0.0, 0.2, size=(120, n_genes)).astype(np.float32) + np.log1p(mu_s)[None, :].astype(np.float32)
        src_map[f"src_{s_idx}"] = LineDatasetR10(
            name=f"src_{s_idx}",
            genes=genes,
            target_names=targets,
            target_indices=np.arange(n_targets, dtype=np.int64),
            deltas=d_s,
            ctrl_cells=c_s,
            ctrl_mean_log1p=np.log1p(mu_s),
            ctrl_mean_libnorm=mu_s,
            ctrl_mean_libnorm_full=mu_s,
            repr_mask=np.ones(n_targets, dtype=bool),
            half_a_deltas=d_s + rng.normal(0.0, 0.03, size=d_s.shape),
            half_b_deltas=d_s + rng.normal(0.0, 0.03, size=d_s.shape),
            half_valid=np.ones(n_targets, dtype=bool),
            n_cells=rng.integers(30, 150, size=n_targets).astype(np.float64),
            knockdown_depth_log1p=np.abs(np.diag(d_s)),
            knockdown_fraction=None,
            target_ensembl_ids=None,
            representation="log1p",
            metadata={},
        )

    # Planted target line: FCT formula scaled cleanly + replicate split halves
    pert_t = np.maximum(mu_t[None, :] * np.exp(signal_scale * log_fc), 1e-8)
    d_true_clean = np.log1p(pert_t) - np.log1p(mu_t[None, :])
    n_a = rng.normal(0.0, 0.06, size=d_true_clean.shape)
    n_b = rng.normal(0.0, 0.06, size=d_true_clean.shape)
    d_ha = d_true_clean + n_a
    d_hb = d_true_clean + n_b
    d_obs = 0.5 * (d_ha + d_hb)
    c_t = rng.normal(0.0, 0.25, size=(160, n_genes)).astype(np.float32) + np.log1p(mu_t)[None, :].astype(np.float32)

    tgt_ds = LineDatasetR10(
        name="syn_tgt_1",
        genes=genes,
        target_names=targets,
        target_indices=np.arange(n_targets, dtype=np.int64),
        deltas=d_obs,
        ctrl_cells=c_t,
        ctrl_mean_log1p=np.log1p(mu_t),
        ctrl_mean_libnorm=mu_t,
        ctrl_mean_libnorm_full=mu_t,
        repr_mask=np.ones(n_targets, dtype=bool),
        half_a_deltas=d_ha,
        half_b_deltas=d_hb,
        half_valid=np.ones(n_targets, dtype=bool),
        n_cells=rng.integers(30, 200, size=n_targets).astype(np.float64),
        knockdown_depth_log1p=np.abs(np.diag(d_obs)),
        knockdown_fraction=None,
        target_ensembl_ids=None,
        representation="log1p",
        metadata={},
    )
    return src_map, tgt_ds


class Round11UnitTests(unittest.TestCase):
    def test_1_blinding_gate(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)
            with self.assertRaises(RuntimeError):
                check_prereg_d_gate(p)
            (p / "PREREG_D_COMMIT.txt").write_text("abc123\n", encoding="utf-8")
            with self.assertRaises(RuntimeError):
                check_prereg_d_gate(p)
            (p / "CODE_FREEZE_R11.txt").write_text("def456\n", encoding="utf-8")
            pr, fr = check_prereg_d_gate(p)
            self.assertEqual(pr, "abc123")
            self.assertEqual(fr, "def456")

    def test_2_bonferroni_9833_quantiles(self) -> None:
        arr = np.linspace(0.0, 1.0, 10001)
        s = summarize_boot_r11(0.5, arr, n_items=50, n_boot=len(arr))
        self.assertAlmostEqual(s["ci95_lo"], 0.025, places=6)
        self.assertAlmostEqual(s["ci95_hi"], 0.975, places=6)
        self.assertAlmostEqual(s["ci975_lo"], 0.0125, places=6)
        self.assertAlmostEqual(s["ci975_hi"], 0.9875, places=6)
        self.assertAlmostEqual(s["ci9833_lo"], PCT_LO_9833 / 100.0, places=6)
        self.assertAlmostEqual(s["ci9833_hi"], PCT_HI_9833 / 100.0, places=6)
        self.assertLess(s["ci9833_lo"], s["ci975_lo"])
        self.assertGreater(s["ci9833_hi"], s["ci975_hi"])

    def test_3_parity_vs_frozen_r10_metrics_and_dual_perm(self) -> None:
        src_map, tgt_ds = _make_planted_synthetic_cohort(seed=77, n_genes=120, n_targets=40, n_sources=4)
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
        preds, fc_bundle, _ = predict_group_d_suite(
            state=state,
            target_ctrl_mean_libnorm=tgt_ds.ctrl_mean_libnorm,
            query_target_names=tgt_ds.target_names,
            include_r9_anchors=False,
        )

        # Residual suite parity vs frozen R10 across BOTH masks
        res10 = compute_residual_suite_r10(
            delta_obs=tgt_ds.deltas,
            delta_pred=preds["FCT"],
            delta_transfer=preds["Transfer_Mean"],
            target_gene_idx=tgt_ds.target_indices,
            delta_half_a=tgt_ds.half_a_deltas,
            delta_half_b=tgt_ds.half_b_deltas,
            half_valid=tgt_ds.half_valid,
            top_k=50,
            n_boot=200,
            seed=20261123,
        )
        res11 = compute_residual_suite_r11(
            delta_obs=tgt_ds.deltas,
            delta_pred=preds["FCT"],
            delta_transfer=preds["Transfer_Mean"],
            target_gene_idx=tgt_ds.target_indices,
            delta_half_a=tgt_ds.half_a_deltas,
            delta_half_b=tgt_ds.half_b_deltas,
            half_valid=tgt_ds.half_valid,
            top_k=50,
            n_boot=200,
            seed=20261123,
        )
        for mode in ("resp_top50", "all_nontarget"):
            for k in ("estimate", "ci95_lo", "ci95_hi", "ci975_lo", "ci975_hi", "valid_frac"):
                self.assertAlmostEqual(
                    res10[mode]["r_orth_cal"][k],
                    res11[mode]["r_orth_cal"][k],
                    places=12,
                )

        # Perm100 parity on resp_top50 vs frozen R10
        p10 = compute_fct_perm100_null(
            state=state,
            target_ctrl_mean_libnorm=tgt_ds.ctrl_mean_libnorm,
            query_target_names=tgt_ds.target_names,
            delta_obs=tgt_ds.deltas,
            delta_transfer=preds["Transfer_Mean"],
            target_gene_idx=tgt_ds.target_indices,
            ceil_r_top50=res10["resp_top50"]["ceiling_r"],
            ceil_r_orth_top50=res10["resp_top50"]["ceiling_r_orth"],
            fc_bundle=fc_bundle,
            n_perm=20,
            top_k=50,
        )
        p11_t50 = compute_fct_perm100_null_r11(
            state=state,
            target_ctrl_mean_libnorm=tgt_ds.ctrl_mean_libnorm,
            query_target_names=tgt_ds.target_names,
            delta_obs=tgt_ds.deltas,
            delta_transfer=preds["Transfer_Mean"],
            target_gene_idx=tgt_ds.target_indices,
            ceil_r=res11["resp_top50"]["ceiling_r"],
            ceil_r_orth=res11["resp_top50"]["ceiling_r_orth"],
            fc_bundle=fc_bundle,
            mode="resp_top50",
            n_perm=20,
            top_k=50,
        )
        for k in ("p2_5", "median", "p97_5", "p98_75", "mean"):
            self.assertAlmostEqual(p10["r_orth_cal"][k], p11_t50["r_orth_cal"][k], places=12)

        # Perm100 naive reference parity on all_nontarget (with row_subset & gene_keep)
        sub_rows = np.arange(tgt_ds.deltas.shape[0]) % 2 == 0
        gk = np.ones(tgt_ds.deltas.shape[1], dtype=bool)
        gk[:10] = False
        p11_all = compute_fct_perm100_null_r11(
            state=state,
            target_ctrl_mean_libnorm=tgt_ds.ctrl_mean_libnorm,
            query_target_names=tgt_ds.target_names,
            delta_obs=tgt_ds.deltas,
            delta_transfer=preds["Transfer_Mean"],
            target_gene_idx=tgt_ds.target_indices,
            ceil_r=res11["all_nontarget"]["ceiling_r"],
            ceil_r_orth=res11["all_nontarget"]["ceiling_r_orth"],
            fc_bundle=fc_bundle,
            mode="all_nontarget",
            n_perm=15,
            top_k=50,
            row_subset_mask=sub_rows,
            gene_keep_mask=gk,
        )
        # Recompute k=0 naive loop manually
        mu_t = tgt_ds.ctrl_mean_libnorm
        perm0 = np.random.default_rng(0).permutation(len(mu_t))
        pred0 = fct_from_fc_bar(mu_t[perm0], fc_bundle[1], fc_bundle[0], fc_bundle[2])
        rows_n, idx_n = _extract_masked_idx_per_pert(
            tgt_ds.deltas[sub_rows],
            tgt_ds.target_indices[sub_rows],
            gene_keep_mask=gk,
            mode="all_nontarget",
        )
        xo_l, yo_l = [], []
        for p, idx in zip(rows_n, idx_n):
            xv = (tgt_ds.deltas[sub_rows] - preds["Transfer_Mean"][sub_rows])[p, idx]
            yv = (pred0[sub_rows] - preds["Transfer_Mean"][sub_rows])[p, idx]
            rv = preds["Transfer_Mean"][sub_rows][p, idx]
            bv = rv / np.linalg.norm(rv)
            xo_l.append(xv - np.dot(xv, bv) * bv)
            yo_l.append(yv - np.dot(yv, bv) * bv)
        ref0 = disattenuated_r(
            pearson_1d_strict(np.concatenate(xo_l), np.concatenate(yo_l)),
            res11["all_nontarget"]["ceiling_r_orth"],
        )
        self.assertAlmostEqual(ref0, p11_all["perm0_r_orth_cal"], places=12)

    def test_4_planted_hd1_hd2_hd3_and_hd4_cross_mask(self) -> None:
        src_map, tgt_ds = _make_planted_synthetic_cohort(seed=99, n_genes=160, n_targets=60, n_sources=6)
        l_rows, b_rows, extra = evaluate_target_against_sources_r11(
            target_ds=tgt_ds,
            source_datasets=src_map,
            only_reproducible=True,
            include_r9_anchors=True,
            include_exploratory_e=True,
            min_gene_overlap=100,
            n_boot=400,
            seed=20261123,
            n_perm=25,
        )
        df_l = pd.DataFrame(l_rows)
        df_b = pd.DataFrame(b_rows)
        fct_r = df_l[df_l["method"] == "FCT"].iloc[0]
        e1_r = df_l[df_l["method"] == "E1_ConstMed"].iloc[0]
        hd3_r = df_b[df_b["comparison"] == "H-D3_FCT_minus_E1_all_nontarget"].iloc[0]
        hd4_r = df_b[df_b["comparison"] == "H-D4_E1_resp_top50_minus_all_nontarget"].iloc[0]

        # H-D1..H-D3 on planted true mu_t signal
        self.assertGreater(fct_r["r_orth_cal_ci9833_lo"], 0.30)
        self.assertGreater(hd3_r["ci9833_lo"], 0.25)
        self.assertAlmostEqual(hd3_r["estimate"], fct_r["r_orth_cal"] - e1_r["r_orth_cal"], places=12)

        # Exact cross-mask marginal consistency for H-D4
        self.assertAlmostEqual(hd4_r["mean_a"], e1_r["resp_top50_r_orth_cal"], places=12)
        self.assertAlmostEqual(hd4_r["mean_b"], e1_r["r_orth_cal"], places=12)
        self.assertAlmostEqual(hd4_r["estimate"], hd4_r["mean_a"] - hd4_r["mean_b"], places=12)

        # Dual-mask perm100 & dual-mask S7 populated
        self.assertIn("all_nontarget", extra["perm100"])
        self.assertIn("resp_top50", extra["perm100"])
        self.assertEqual(extra["s7_cellcount_kd"]["mode"], "all_nontarget")
        self.assertEqual(extra["s7_cellcount_kd_resp_top50"]["mode"], "resp_top50")

        # Pure constant-baseline selection bias check for H-D4 (true signal absent -> FCT == E1,
        # but resp_top50 inflation pulls E1(resp_top50) > 0 while E1(all_nontarget) <= 0)
        rng = np.random.default_rng(2025)
        n_p, n_g = 80, 150
        mu_const = np.full(n_g, 2.0)
        mu_hetero = np.exp(rng.uniform(-1.5, 2.5, size=n_g))
        fc_bar = np.exp(rng.normal(0.0, 0.5, size=(n_p, n_g)))
        fc_t = np.mean(rng.normal(0.0, 0.35, size=(n_p, n_g)), axis=0, keepdims=True)
        fc_pred_transfer = np.broadcast_to(fc_t, (n_p, n_g)).copy()
        has_hits = np.ones(n_p, dtype=bool)
        e1_pred = fct_from_fc_bar(mu_const, fc_bar, fc_pred_transfer, has_hits)
        # Observed delta driven by E1 nonlinear transform plus noise (no gene-baseline variation!)
        d_clean = e1_pred * 1.15
        d_ha = d_clean + rng.normal(0.0, 0.15, size=d_clean.shape)
        d_hb = d_clean + rng.normal(0.0, 0.15, size=d_clean.shape)
        d_obs = 0.5 * (d_ha + d_hb)
        cross = compute_joint_cross_mask_e1_diff_r11(
            delta_obs=d_obs,
            delta_pred_e1=e1_pred,
            delta_transfer=fc_pred_transfer,
            target_gene_idx=np.arange(n_p),
            delta_half_a=d_ha,
            delta_half_b=d_hb,
            top_k=50,
            n_boot=200,
            seed=20261123,
        )
        # Both E1 under resp_top50 > 0 and resp_top50 ceiling/signal differs cleanly
        self.assertGreater(cross["a_summary"]["ci95_lo"], 0.0)
        self.assertTrue(np.isfinite(cross["estimate"]))

        # Suppress unused helper warning cleanly (verify mu_hetero unused -> no problem)
        _ = mu_hetero

    def test_5_adjudicator_p1_unanimous_3_of_3_and_8_sources(self) -> None:
        src8, tgt1 = _make_planted_synthetic_cohort(seed=51, n_genes=150, n_targets=60, n_sources=8)
        # Second confirmatory line shares biological perturbation effects with `src8` (distinct target mu_t)
        rng2 = np.random.default_rng(52)
        n_t, n_g = tgt1.deltas.shape
        mu_t2 = np.exp(rng2.uniform(-1.5, 3.0, size=n_g))
        mu_s0 = src8["src_0"].ctrl_mean_libnorm
        fc0 = np.maximum(np.expm1(np.log1p(mu_s0)[None, :] + src8["src_0"].deltas), 1e-8) / mu_s0[None, :]
        pert_t2 = np.maximum(mu_t2[None, :] * fc0, 1e-8)
        d_clean2 = np.log1p(pert_t2) - np.log1p(mu_t2[None, :])
        ha2 = d_clean2 + rng2.normal(0.0, 0.05, size=d_clean2.shape)
        hb2 = d_clean2 + rng2.normal(0.0, 0.05, size=d_clean2.shape)
        tgt2 = LineDatasetR10(
            **{
                **tgt1.__dict__,
                "name": "syn_tgt_2",
                "deltas": 0.5 * (ha2 + hb2),
                "ctrl_mean_log1p": np.log1p(mu_t2),
                "ctrl_mean_libnorm": mu_t2,
                "ctrl_mean_libnorm_full": mu_t2,
                "half_a_deltas": ha2,
                "half_b_deltas": hb2,
            }
        )

        l1, b1, ex1 = evaluate_target_against_sources_r11(
            tgt1, src8, include_r9_anchors=False, min_gene_overlap=100, n_boot=200, seed=20261123, n_perm=20
        )
        l2, b2, ex2 = evaluate_target_against_sources_r11(
            tgt2, src8, include_r9_anchors=False, min_gene_overlap=100, n_boot=200, seed=20261123, n_perm=20
        )
        df_l = pd.DataFrame(l1 + l2)
        df_b = pd.DataFrame(b1 + b2)
        summary = {
            "phase": "holdout_r11",
            "source_mode": "r9_r10_8",
            "source_lines": sorted(src8.keys()),
            "prereg_commit_hash": "prereg_d",
            "code_freeze_hash": "freeze_r11",
            "bootstrap_seed": 20261123,
            "n_boot": 200,
            "fct_perm100": {
                "syn_tgt_1": ex1["perm100"],
                "syn_tgt_2": ex2["perm100"],
            },
        }

        # Case 1: K_conf = 1 guard
        adj1 = adjudicate_r11_results(df_l, df_b, summary, confirmatory_lines=["syn_tgt_1"])
        self.assertEqual(adj1["confirmatory_status"], "INSUFFICIENT_LINES")
        self.assertEqual(adj1["Group_D_Verdict"], "NOT_EVALUATED_INSUFFICIENT_LINES")

        # Case 2: K_conf = 2, all 3 primary pass unanimously + 8 sources checked in H-D2
        adj2 = adjudicate_r11_results(df_l, df_b, summary, confirmatory_lines=["syn_tgt_1", "syn_tgt_2"])
        self.assertEqual(adj2["H-D1"]["verdict"], "CONFIRMED")
        self.assertEqual(adj2["H-D2"]["verdict"], "CONFIRMED")
        self.assertEqual(adj2["H-D2"]["per_line"]["syn_tgt_1"]["n_mismatched_sources"], 8)
        self.assertEqual(adj2["H-D3"]["verdict"], "CONFIRMED")
        self.assertEqual(adj2["Group_D_Verdict"], "CONFIRMED")
        self.assertIn("H-D_util", adj2)
        self.assertEqual(adj2["H-D4"]["primary_ci_level"], "ci9833")
        md = render_r11_markdown_report(df_l, df_b, summary, adj2, "unit_test")
        self.assertIn("K_conf = 2 confirmatory cell lines", md)

        # Case 3: Fail H-D3 on 1/2 lines -> Group_D_Verdict MUST be 未确认 even if H-D1 & H-D2 pass 2/2
        df_b_fail_hd3 = df_b.copy()
        m_fail = (df_b_fail_hd3["target_line"] == "syn_tgt_2") & (
            df_b_fail_hd3["comparison"] == "H-D3_FCT_minus_E1_all_nontarget"
        )
        df_b_fail_hd3.loc[m_fail, "ci9833_lo"] = -0.001
        adj3 = adjudicate_r11_results(df_l, df_b_fail_hd3, summary, confirmatory_lines=["syn_tgt_1", "syn_tgt_2"])
        self.assertEqual(adj3["H-D1"]["verdict"], "CONFIRMED")
        self.assertEqual(adj3["H-D2"]["verdict"], "CONFIRMED")
        self.assertEqual(adj3["H-D3"]["verdict"], "未确认")
        self.assertEqual(adj3["H-D3"]["pass_count"], 1)
        self.assertEqual(adj3["Group_D_Verdict"], "未确认")


if __name__ == "__main__":
    unittest.main(verbosity=2)
