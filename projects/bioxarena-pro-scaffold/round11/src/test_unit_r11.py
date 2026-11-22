#!/usr/bin/env python3
"""
Unit tests for Round 11 (`PREREG_D`) Group D evaluation modules (`round11/src/`).

Runs cleanly on synthetic fixtures and checks exact mathematical parity against
frozen Round 10 primitives without accessing any Round 11 held-out data.
"""

from __future__ import annotations

import csv
import gzip
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Dict, List

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

R11_SRC_DIR = Path(__file__).resolve().parent
REPO_ROOT = R11_SRC_DIR.parents[1]
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
R10_SCRIPTS_DIR = REPO_ROOT / "round10" / "shared" / "scripts"
for p in (str(R10_SRC_DIR), str(R10_SCRIPTS_DIR), str(R11_SRC_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from data_loader_r10 import LineDatasetR10, load_standardized_line_r10  # noqa: E402
from fct_pipeline import FctSourceState, SourceLineInputR10, fct_from_fc_bar  # noqa: E402
from metrics_r10 import (  # noqa: E402
    _extract_masked_idx_per_pert,
    compute_fct_perm100_null,
    compute_residual_suite_r10,
    disattenuated_r,
    pearson_1d_strict,
)
from prepare_r10_lines import HGNCLookup, build_hgnc_lookup, process_r10_line  # noqa: E402

from adjudicate_r11 import adjudicate_r11_results, render_r11_markdown_report  # noqa: E402
from data_loader_r11 import check_prereg_d_gate  # noqa: E402
from dataset_readers_r11 import (  # noqa: E402
    build_cd4t_stim48_parts,
    build_feng_ipsc_parts,
    build_hcasmc_parts,
    build_telohaec_parts,
    check_hcasmc_feasibility,
)
from metrics_r11 import (  # noqa: E402
    PCT_HI_9833,
    PCT_LO_9833,
    compute_fct_perm100_null_r11,
    compute_fct_perm100_split_r11,
    compute_joint_cross_mask_e1_diff_r11,
    compute_residual_suite_r11,
    compute_resp_top50_split_suite_r11,
    single_half_ceiling,
    summarize_boot_r11,
)
from predictors_r11 import predict_group_d_suite  # noqa: E402
from run_eval_r11 import evaluate_target_against_sources_r11  # noqa: E402


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


def _write_10x_csc_h5(
    out_path: Path,
    csc: sp.csc_matrix,
    barcodes: List[str],
    feat_ids: List[str],
    feat_names: List[str],
    feat_types: List[str],
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    csc_u = csc.tocsc()
    with h5py.File(out_path, "w") as f:
        m = f.create_group("matrix")
        m.create_dataset("shape", data=np.array(csc_u.shape, dtype=np.int64))
        m.create_dataset("barcodes", data=np.array(barcodes, dtype="S"))
        m.create_dataset("data", data=csc_u.data.astype(np.float32))
        m.create_dataset("indices", data=csc_u.indices.astype(np.int32))
        m.create_dataset("indptr", data=csc_u.indptr.astype(np.int64))
        fg = m.create_group("features")
        fg.create_dataset("id", data=np.array(feat_ids, dtype="S"))
        fg.create_dataset("name", data=np.array(feat_names, dtype="S"))
        fg.create_dataset("feature_type", data=np.array(feat_types, dtype="S"))


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

        # H-D1..H-D3 on planted true mu_t signal under PREREG_D v2 97.5% CIs
        self.assertGreater(fct_r["r_orth_cal_ci975_lo"], 0.30)
        self.assertGreater(hd3_r["ci975_lo"], 0.25)
        self.assertAlmostEqual(hd3_r["estimate"], fct_r["r_orth_cal"] - e1_r["r_orth_cal"], places=12)

        # Exact cross-mask marginal consistency for H-D4
        self.assertAlmostEqual(hd4_r["mean_a"], e1_r["resp_top50_r_orth_cal"], places=12)
        self.assertAlmostEqual(hd4_r["mean_b"], e1_r["r_orth_cal"], places=12)
        self.assertAlmostEqual(hd4_r["estimate"], hd4_r["mean_a"] - hd4_r["mean_b"], places=12)

        # SA-8 resp_top50_split columns + tri-mask perm100 & dual-mask S7 populated
        self.assertTrue(np.isfinite(fct_r["resp_top50_split_r_orth_cal"]))
        self.assertIn("all_nontarget", extra["perm100"])
        self.assertIn("resp_top50", extra["perm100"])
        self.assertIn("resp_top50_split", extra["perm100"])
        self.assertEqual(extra["s7_cellcount_kd"]["mode"], "all_nontarget")
        self.assertEqual(extra["s7_cellcount_kd_resp_top50"]["mode"], "resp_top50")

        # Targeted panel min_gene_overlap guard: >= 3000 fails on 160-gene panel, >= 50 passes
        with self.assertRaises(RuntimeError):
            evaluate_target_against_sources_r11(
                target_ds=tgt_ds,
                source_datasets=src_map,
                min_gene_overlap=3000,
                n_boot=50,
                n_perm=5,
            )

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
        self.assertGreater(cross["a_summary"]["ci975_lo"], 0.0)
        self.assertTrue(np.isfinite(cross["estimate"]))
        _ = mu_hetero

    def test_5_adjudicator_v2_and_addendum_a1_and_hd5(self) -> None:
        src8, tgt1 = _make_planted_synthetic_cohort(seed=51, n_genes=150, n_targets=60, n_sources=8)
        # Rename first source to nourreddine_kolf21j so H-D5 can also be exercised
        src_named = {"nourreddine_kolf21j": src8["src_0"]}
        for k, v in src8.items():
            if k != "src_0":
                src_named[k] = v
        src_named["nourreddine_kolf21j"].name = "nourreddine_kolf21j"

        rng2 = np.random.default_rng(52)
        n_t, n_g = tgt1.deltas.shape
        mu_t2 = np.exp(rng2.uniform(-1.5, 3.0, size=n_g))
        mu_s0 = src8["src_1"].ctrl_mean_libnorm
        fc0 = np.maximum(np.expm1(np.log1p(mu_s0)[None, :] + src8["src_1"].deltas), 1e-8) / mu_s0[None, :]
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
        # Secondary feng_ipsc synthetic line on same synthetic cohort
        mu_f = np.exp(rng2.uniform(-1.5, 3.0, size=n_g))
        pert_f = np.maximum(mu_f[None, :] * fc0, 1e-8)
        d_clean_f = np.log1p(pert_f) - np.log1p(mu_f[None, :])
        ha_f = d_clean_f + rng2.normal(0.0, 0.05, size=d_clean_f.shape)
        hb_f = d_clean_f + rng2.normal(0.0, 0.05, size=d_clean_f.shape)
        tgt_feng = LineDatasetR10(
            **{
                **tgt1.__dict__,
                "name": "feng_ipsc",
                "deltas": 0.5 * (ha_f + hb_f),
                "ctrl_mean_log1p": np.log1p(mu_f),
                "ctrl_mean_libnorm": mu_f,
                "ctrl_mean_libnorm_full": mu_f,
                "half_a_deltas": ha_f,
                "half_b_deltas": hb_f,
            }
        )

        l1, b1, ex1 = evaluate_target_against_sources_r11(
            tgt1, src_named, include_r9_anchors=False, min_gene_overlap=100, n_boot=200, seed=20261123, n_perm=20
        )
        l2, b2, ex2 = evaluate_target_against_sources_r11(
            tgt2, src_named, include_r9_anchors=False, min_gene_overlap=100, n_boot=200, seed=20261123, n_perm=20
        )
        lf, bf, exf = evaluate_target_against_sources_r11(
            tgt_feng, src_named, include_r9_anchors=False, min_gene_overlap=50, n_boot=200, seed=20261123, n_perm=20
        )
        df_l = pd.DataFrame(l1 + l2 + lf)
        df_b = pd.DataFrame(b1 + b2 + bf)
        summary = {
            "phase": "holdout_r11",
            "source_mode": "r9_r10_8",
            "source_lines": sorted(src_named.keys()),
            "prereg_commit_hash": "prereg_d_v2",
            "code_freeze_hash": "freeze_r11",
            "bootstrap_seed": 20261123,
            "n_boot": 200,
            "fct_perm100": {
                "syn_tgt_1": ex1["perm100"],
                "syn_tgt_2": ex2["perm100"],
                "feng_ipsc": exf["perm100"],
            },
        }

        # Case 1: K_conf = 1 guard
        adj1 = adjudicate_r11_results(df_l, df_b, summary, confirmatory_lines=["syn_tgt_1"])
        self.assertEqual(adj1["confirmatory_status"], "INSUFFICIENT_LINES")
        self.assertEqual(adj1["Group_D_Verdict"], "NOT_EVALUATED_INSUFFICIENT_LINES")

        # Case 2: K_conf = 2, all 3 primary pass unanimously + 8 sources in H-D2 + H-D5 present
        adj2 = adjudicate_r11_results(df_l, df_b, summary, confirmatory_lines=["syn_tgt_1", "syn_tgt_2"])
        self.assertEqual(adj2["H-D1"]["verdict"], "CONFIRMED")
        self.assertEqual(adj2["H-D2"]["verdict"], "CONFIRMED")
        self.assertEqual(adj2["H-D2"]["per_line"]["syn_tgt_1"]["n_mismatched_sources"], 8)
        self.assertEqual(adj2["H-D3"]["verdict"], "CONFIRMED")
        self.assertEqual(adj2["Group_D_Verdict"], "CONFIRMED")
        self.assertIn("H-D_util", adj2)
        self.assertEqual(adj2["H-D4"]["primary_ci_level"], "ci975")
        self.assertTrue(adj2["descriptive_H_D5_feng_ipsc"]["present"])
        self.assertTrue(adj2["descriptive_H_D5_feng_ipsc"]["fct_strictly_exceeds_kolf21j_paired_975"])
        md = render_r11_markdown_report(df_l, df_b, summary, adj2, "unit_test")
        self.assertIn("K_conf = 2 confirmatory cell lines", md)
        self.assertIn("`SA-8` Split-Half Out-of-Sample Top-50", md)
        self.assertIn("H-D5", md)

        # Case 3: Fail H-D3 on 1/2 lines using ci975_lo -> Group_D_Verdict MUST be 未确认
        df_b_fail_hd3 = df_b.copy()
        m_fail = (df_b_fail_hd3["target_line"] == "syn_tgt_2") & (
            df_b_fail_hd3["comparison"] == "H-D3_FCT_minus_E1_all_nontarget"
        )
        df_b_fail_hd3.loc[m_fail, "ci975_lo"] = -0.001
        adj3 = adjudicate_r11_results(df_l, df_b_fail_hd3, summary, confirmatory_lines=["syn_tgt_1", "syn_tgt_2"])
        self.assertEqual(adj3["H-D1"]["verdict"], "CONFIRMED")
        self.assertEqual(adj3["H-D2"]["verdict"], "CONFIRMED")
        self.assertEqual(adj3["H-D3"]["verdict"], "未确认")
        self.assertEqual(adj3["H-D3"]["pass_count"], 1)
        self.assertEqual(adj3["Group_D_Verdict"], "未确认")

        # Case 4 (Addendum A1): ci975_lo > 0 but valid_frac < 0.95 -> H-D1 MUST fail & set Addendum A1 flag
        df_l_a1_vf = df_l.copy()
        m_fct2 = (df_l_a1_vf["target_line"] == "syn_tgt_2") & (df_l_a1_vf["method"] == "FCT")
        df_l_a1_vf.loc[m_fct2, "r_orth_cal_valid_frac"] = 0.90
        adj4 = adjudicate_r11_results(df_l_a1_vf, df_b, summary, confirmatory_lines=["syn_tgt_1", "syn_tgt_2"])
        self.assertEqual(adj4["H-D1"]["verdict"], "未确认")
        self.assertTrue(adj4["H-D1"]["per_line"]["syn_tgt_2"]["addendum_a1_ci_gt0_but_low_valid_frac"])
        self.assertFalse(adj4["H-D1"]["per_line"]["syn_tgt_2"]["valid_frac_ge_095"])
        self.assertEqual(adj4["Group_D_Verdict"], "未确认")

        # Case 5 (Addendum A1): ci975_lo > 0 and valid_frac == 1.0, but rho_half_ab_orth <= 0 -> H-D1 MUST fail
        df_l_a1_rho = df_l.copy()
        df_l_a1_rho.loc[m_fct2, "rho_half_ab_orth"] = 0.0
        adj5 = adjudicate_r11_results(df_l_a1_rho, df_b, summary, confirmatory_lines=["syn_tgt_1", "syn_tgt_2"])
        self.assertEqual(adj5["H-D1"]["verdict"], "未确认")
        self.assertFalse(adj5["H-D1"]["per_line"]["syn_tgt_2"]["rho_ab_orth_gt_0"])
        self.assertEqual(adj5["Group_D_Verdict"], "未确认")

    def test_6_sa8_resp_top50_split_single_half_ceiling_and_perm(self) -> None:
        src_map, tgt_ds = _make_planted_synthetic_cohort(seed=66, n_genes=120, n_targets=40, n_sources=4)
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

        top_k = 30
        sa8 = compute_resp_top50_split_suite_r11(
            delta_pred=preds["FCT"],
            delta_transfer=preds["Transfer_Mean"],
            target_gene_idx=tgt_ds.target_indices,
            delta_half_a=tgt_ds.half_a_deltas,
            delta_half_b=tgt_ds.half_b_deltas,
            half_valid=tgt_ds.half_valid,
            top_k=top_k,
            n_boot=120,
            seed=20261123,
        )

        # Manual reference loop verifying single-half ceiling sqrt(rho_AB_orth) != Spearman-Brown
        P, G = tgt_ds.deltas.shape
        t_idx = tgt_ds.target_indices
        tm = preds["Transfer_Mean"]
        pred_r = preds["FCT"] - tm
        ha_r = tgt_ds.half_a_deltas - tm
        hb_r = tgt_ds.half_b_deltas - tm

        def _manual_dir(sel_mat: np.ndarray, eval_mat: np.ndarray, other_mat: np.ndarray) -> tuple[float, float, float]:
            xo_l, yo_l, hao_l, hbo_l = [], [], [], []
            for p in range(P):
                tg = int(t_idx[p])
                cand = [j for j in range(G) if j != tg]
                order = sorted(cand, key=lambda j: -abs(float(sel_mat[p, j])))[:top_k]
                idx = np.asarray(order, dtype=np.int64)
                rv = tm[p, idx]
                bv = rv / np.linalg.norm(rv)
                for src_m, dst_l in ((eval_mat, xo_l), (pred_r, yo_l), (eval_mat, hao_l), (other_mat, hbo_l)):
                    vv = src_m[p, idx]
                    dst_l.append(vv - float(np.dot(vv, bv)) * bv)
            r_raw = pearson_1d_strict(np.concatenate(xo_l), np.concatenate(yo_l))
            rho_ab = pearson_1d_strict(np.concatenate(hao_l), np.concatenate(hbo_l))
            ceil_single = float(np.sqrt(rho_ab)) if rho_ab > 0.0 else float("nan")
            self.assertAlmostEqual(ceil_single, single_half_ceiling(rho_ab), places=12)
            # Must NOT equal Spearman-Brown full-average ceiling
            sb = float(np.sqrt((2.0 * rho_ab) / (1.0 + rho_ab)))
            self.assertNotAlmostEqual(ceil_single, sb, places=3)
            return disattenuated_r(r_raw, ceil_single), rho_ab, ceil_single

        cal_a, rho_a, ceil_a = _manual_dir(tgt_ds.half_a_deltas, hb_r, ha_r)
        cal_b, rho_b, ceil_b = _manual_dir(tgt_ds.half_b_deltas, ha_r, hb_r)
        self.assertAlmostEqual(sa8["r_orth_cal_a_to_b"], cal_a, places=12)
        self.assertAlmostEqual(sa8["r_orth_cal_b_to_a"], cal_b, places=12)
        self.assertAlmostEqual(sa8["ceiling_single_a"], ceil_a, places=12)
        self.assertAlmostEqual(sa8["ceiling_single_b"], ceil_b, places=12)
        self.assertAlmostEqual(sa8["r_orth_cal"]["estimate"], 0.5 * (cal_a + cal_b), places=12)

        # Perm100 under SA-8 perm0 must match SA-8 FCT estimate exactly
        p100_split = compute_fct_perm100_split_r11(
            state=state,
            target_ctrl_mean_libnorm=tgt_ds.ctrl_mean_libnorm,
            query_target_names=tgt_ds.target_names,
            delta_transfer=preds["Transfer_Mean"],
            target_gene_idx=tgt_ds.target_indices,
            delta_half_a=tgt_ds.half_a_deltas,
            delta_half_b=tgt_ds.half_b_deltas,
            half_valid=tgt_ds.half_valid,
            fc_bundle=fc_bundle,
            n_perm=10,
            top_k=top_k,
        )
        # Recompute unpermuted point check against perm0 manually on mu_t[rng(0).permutation(G)]
        rng0 = np.random.default_rng(0).permutation(G)
        pred_tr, fc_bar_mat, has_hits_arr = fc_bundle
        pred_perm0 = fct_from_fc_bar(tgt_ds.ctrl_mean_libnorm[rng0], fc_bar_mat, pred_tr, has_hits_arr)
        sa8_p0 = compute_resp_top50_split_suite_r11(
            delta_pred=pred_perm0,
            delta_transfer=preds["Transfer_Mean"],
            target_gene_idx=tgt_ds.target_indices,
            delta_half_a=tgt_ds.half_a_deltas,
            delta_half_b=tgt_ds.half_b_deltas,
            half_valid=tgt_ds.half_valid,
            top_k=top_k,
            n_boot=20,
            seed=1,
        )
        self.assertAlmostEqual(p100_split["perm0_r_orth_raw"], sa8_p0["r_orth_raw"]["estimate"], places=12)

        # Pure selection-noise null fixture: E[delta_true] = 0 -> naive resp_top50 inflates over SA-8 split
        rng_n = np.random.default_rng(999)
        n_p_n, n_g_n = 80, 200
        tm_zero = rng_n.normal(0.0, 0.15, size=(n_p_n, n_g_n))
        pred_bias = tm_zero + rng_n.normal(0.0, 0.25, size=(n_p_n, n_g_n))
        # Halves share zero biological signal, only modest shared selection noise artifact vs pure half noise
        d_ha_n = rng_n.normal(0.0, 0.35, size=(n_p_n, n_g_n)) + 0.12 * (pred_bias - tm_zero)
        d_hb_n = rng_n.normal(0.0, 0.35, size=(n_p_n, n_g_n)) + 0.12 * (pred_bias - tm_zero)
        # Force E1-like constant-baseline direction alignment on extreme random rows of avg(ha, hb)
        d_avg_n = 0.5 * (d_ha_n + d_hb_n)
        pred_aligned = tm_zero + 0.5 * d_avg_n
        suite_std = compute_residual_suite_r11(
            delta_obs=d_avg_n,
            delta_pred=pred_aligned,
            delta_transfer=tm_zero,
            target_gene_idx=np.arange(n_p_n),
            delta_half_a=d_ha_n,
            delta_half_b=d_hb_n,
            top_k=20,
        )
        suite_split = compute_resp_top50_split_suite_r11(
            delta_pred=pred_aligned,
            delta_transfer=tm_zero,
            target_gene_idx=np.arange(n_p_n),
            delta_half_a=d_ha_n,
            delta_half_b=d_hb_n,
            top_k=20,
            n_boot=40,
            seed=20261123,
        )
        self.assertGreater(
            suite_std["resp_top50"]["r_orth_raw"]["estimate"],
            suite_split["r_orth_raw"]["estimate"] + 0.05,
        )

    def test_7_dataset_reader_telohaec_synthetic(self) -> None:
        rng = np.random.default_rng(701)
        hgnc = _get_hgnc()
        genes = _pick_hgnc_genes(520)
        var_names = genes + ["MT-CO1", "GuideCaptureFeat"]
        var_types = ["Gene Expression"] * (len(genes) + 1) + ["CRISPR Guide Capture"]
        n_var = len(var_names)

        # Build 4 lanes (1,3 odd -> Half A; 2,4 even -> Half B) with NTC & 2 targets
        # plus excluded safe-targeting, -TSS2, Enhancer-1, low-P4, and high-MT cells
        cell_specs = []
        for lane in (1, 2, 3, 4):
            for i in range(18):
                cell_specs.append((f"NTC_{lane}_{i:02d}-{lane}", lane, "ntc", f"g_ntc_{i % 2}", "", 40, 0, False))
            for tgt_idx, gsym in enumerate((genes[0], genes[1])):
                for i in range(16):
                    cell_specs.append((f"T{tgt_idx}_{lane}_{i:02d}-{lane}", lane, "tgt", f"g_{gsym}", gsym, 35, 1, False))
        # Exclusions that MUST NOT survive filtering
        cell_specs.append(("SAFE_1-1", 1, "safe", "g_safe_1", "safe_target", 40, 0, False))
        cell_specs.append(("TSS_1-1", 1, "tss", "g_tss_1", f"{genes[2]}-TSS2", 40, 0, False))
        cell_specs.append(("ENH_1-2", 2, "enh", "g_enh_1", "Enhancer-1", 40, 0, False))
        cell_specs.append(("LOWP4_1-1", 1, "ntc", "g_ntc_0", "", 10, 6, False))
        cell_specs.append(("HIGHMT_1-1", 1, "ntc", "g_ntc_0", "", 40, 0, True))

        n_cells = len(cell_specs)
        # Positive integer RNA counts across all 520 HGNC genes so ngenes>=500 & sd_lib>0 across every lane
        X = rng.integers(1, 5, size=(n_cells, n_var), dtype=np.int32)
        for ci, (_, _, kind, _, gsym, _, _, high_mt) in enumerate(cell_specs):
            X[ci, len(genes)] = 2000 if high_mt else 5  # MT-CO1
            X[ci, -1] = 10  # non-RNA feature
            if kind == "tgt":
                # Consistent planted perturbation response across halves
                X[ci, :25] += 8
                g_pos = genes.index(gsym)
                X[ci, g_pos] = 1
        X_csr = sp.csr_matrix(X, dtype=np.float32)

        def _write_h5ad(path: Path, mat: sp.spmatrix, enc: str) -> None:
            with h5py.File(path, "w") as f:
                xg = f.create_group("X")
                xg.attrs["encoding-type"] = enc
                xg.attrs["shape"] = np.array(mat.shape, dtype=np.int64)
                xg.create_dataset("data", data=mat.data.astype(np.float32))
                xg.create_dataset("indices", data=mat.indices.astype(np.int32))
                xg.create_dataset("indptr", data=mat.indptr.astype(np.int64))
                og = f.create_group("obs")
                og.create_dataset("_index", data=np.array([c[0] for c in cell_specs], dtype="S"))
                vg = f.create_group("var")
                vg.create_dataset("_index", data=np.array(var_names, dtype="S"))
                vg.create_dataset("gene_symbol", data=np.array(var_names, dtype="S"))
                vg.create_dataset("feature_types", data=np.array(var_types, dtype="S"))

        with tempfile.TemporaryDirectory() as td:
            raw_dir = Path(td) / "raw"
            tdir = raw_dir / "telohaec"
            tdir.mkdir(parents=True)
            _write_h5ad(tdir / "IGVFFI8645YMGX.h5ad", X_csr, "csr_matrix")

            # Guide TSV with safe-targeting intentionally tagged targeting=True
            guide_rows = [
                {"guide_id": "g_ntc_0", "type": "non-targeting", "targeting": "False", "intended_target_name": ""},
                {"guide_id": "g_ntc_1", "type": "non-targeting", "targeting": "False", "intended_target_name": ""},
                {"guide_id": "g_safe_1", "type": "safe-targeting", "targeting": "True", "intended_target_name": "safe"},
                {"guide_id": f"g_{genes[0]}", "type": "targeting", "targeting": "True", "intended_target_name": genes[0]},
                {"guide_id": f"g_{genes[1]}", "type": "targeting", "targeting": "True", "intended_target_name": genes[1]},
                {"guide_id": "g_tss_1", "type": "targeting", "targeting": "True", "intended_target_name": f"{genes[2]}-TSS2"},
                {"guide_id": "g_enh_1", "type": "targeting", "targeting": "True", "intended_target_name": "Enhancer-1"},
            ]
            pd.DataFrame(guide_rows).to_csv(tdir / "IGVFFI9754AGFB.tsv.gz", sep="\t", index=False, compression="gzip")

            # Dialout table with digit-prefixed `2nd_place_umi_ct` column name
            dial_rows = [
                {
                    "CBC": cbc,
                    "guide_id": gid,
                    "gene_symbol": gsym,
                    "max_umi_ct": m1,
                    "2nd_place_umi_ct": m2,
                }
                for (cbc, _, _, gid, gsym, m1, m2, _) in cell_specs
            ]
            pd.DataFrame(dial_rows).to_csv(
                tdir / "GSE210681_Dialout_guides_per_CBC.txt.gz", sep="\t", index=False, compression="gzip"
            )

            parts, audit = build_telohaec_parts(raw_dir, hgnc)
            self.assertEqual(audit["n_safe_targeting_excluded"], 1)
            self.assertEqual(audit["n_tss_excluded"], 1)
            self.assertEqual(audit["n_nonplain_or_unapproved_excluded"], 1)
            self.assertEqual(audit["n_qc_kept_cells"], (18 + 32) * 4)
            p = parts[0]
            self.assertEqual(p.n_genes_raw, len(genes) + 1)
            self.assertEqual(set(np.unique(p.half_flags).tolist()), {1, 2})

            # End-to-end standardization check
            out_root = Path(td) / "lines"
            process_r10_line("telohaec", parts, hgnc, out_root=out_root, extra_meta={"adapter_audit": audit}, n_workers=1)
            ds = load_standardized_line_r10(out_root / "telohaec", min_gene_overlap=50)
            self.assertEqual(len(ds.target_names), 2)
            self.assertTrue(np.any(ds.repr_mask))

            # Also verify CSC `.h5ad` sparse QC + read path cleanly
            _write_h5ad(tdir / "IGVFFI8645YMGX.h5ad", X_csr.tocsc(), "csc_matrix")
            parts_csc, audit_csc = build_telohaec_parts(raw_dir, hgnc)
            self.assertEqual(audit_csc["n_qc_kept_cells"], audit["n_qc_kept_cells"])
            self.assertTrue(np.allclose(parts_csc[0].umi, p.umi))

    def test_8_dataset_reader_hcasmc_synthetic(self) -> None:
        rng = np.random.default_rng(801)
        hgnc = _get_hgnc()
        genes = _pick_hgnc_genes(520)
        feat_names = genes + ["MT-CO1", "non-targeting_1", "0Safe_1", f"g_{genes[0]}", f"g_{genes[1]}"]
        feat_ids = feat_names[:]
        feat_types = ["Gene Expression"] * (len(genes) + 1) + ["CRISPR Guide Capture"] * 4
        n_feat = len(feat_names)
        idx_ntc = len(genes) + 1
        idx_safe = len(genes) + 2
        idx_g0 = len(genes) + 3
        idx_g1 = len(genes) + 4

        sublibs = ["3A", "3B", "4A", "4B"]
        barcodes: List[str] = []
        col_assign: List[tuple[int, int]] = []  # (guide_row_idx, tgt_gene_idx_or_-1)
        for sl in sublibs:
            for i in range(18):
                barcodes.append(f"AA{i:04d}__{sl}")
                col_assign.append((idx_ntc, -1))
            for i in range(16):
                barcodes.append(f"BB{i:04d}__{sl}")
                col_assign.append((idx_g0, 0))
            for i in range(16):
                barcodes.append(f"CC{i:04d}__{sl}")
                col_assign.append((idx_g1, 1))
        # Extra safe-targeting cells that MUST be excluded even though guide has >=20 UMI & targeting=True in TSV
        for i in range(6):
            barcodes.append(f"SF{i:04d}__3A")
            col_assign.append((idx_safe, -1))

        n_cells = len(barcodes)
        mat = np.zeros((n_feat, n_cells), dtype=np.float32)
        mat[: len(genes), :] = rng.integers(1, 5, size=(len(genes), n_cells), dtype=np.int32)
        mat[len(genes), :] = 5  # MT-CO1 low (<10%)
        # Low background noise (UMI=1) across half the guide matrix so bimodal EM finds clean threshold
        for ci, (g_row, tgt_i) in enumerate(col_assign):
            if ci % 2 == 0:
                mat[idx_ntc:idx_g1 + 1, ci] = 1
            mat[g_row, ci] = 25
            if tgt_i >= 0:
                mat[:25, ci] += 8
                mat[tgt_i, ci] = 1

        with tempfile.TemporaryDirectory() as td:
            raw_dir = Path(td) / "raw"
            hdir = raw_dir / "hcasmc"
            _write_10x_csc_h5(hdir / "IGVFFI5999NXAF.h5", sp.csc_matrix(mat), barcodes, feat_ids, feat_names, feat_types)
            pd.DataFrame(
                [
                    {"guide_id": "non-targeting_1", "type": "non-targeting", "targeting": "False", "intended_target_name": ""},
                    {"guide_id": "0Safe_1", "type": "safe-targeting", "targeting": "True", "intended_target_name": "safe"},
                    {"guide_id": f"g_{genes[0]}", "type": "targeting", "targeting": "True", "intended_target_name": genes[0]},
                    {"guide_id": f"g_{genes[1]}", "type": "targeting", "targeting": "True", "intended_target_name": genes[1]},
                ]
            ).to_csv(hdir / "IGVFFI_guides.tsv.gz", sep="\t", index=False, compression="gzip")

            feas = check_hcasmc_feasibility(raw_dir, hgnc=hgnc, batch_cols=60)
            self.assertTrue(feas["feasible"])
            self.assertGreaterEqual(feas["n_safe_guides_excluded_meta"], 1)
            self.assertEqual(feas["n_matched_ntc_guides"], 1)

            parts, audit = build_hcasmc_parts(raw_dir, hgnc, batch_cols=60)
            self.assertFalse(audit["used_stratified_seed42_backup"])
            self.assertEqual(audit["n_safe_targeting_excluded"], 6)
            self.assertEqual(audit["n_qc_kept_cells"], 4 * (18 + 16 + 16))
            self.assertEqual(set(np.unique(parts[0].half_flags).tolist()), {1, 2})

            out_root = Path(td) / "lines"
            process_r10_line("hcasmc", parts, hgnc, out_root=out_root, extra_meta={"adapter_audit": audit}, n_workers=1)
            ds = load_standardized_line_r10(out_root / "hcasmc", min_gene_overlap=50)
            self.assertEqual(len(ds.target_names), 2)

            # Single-sublib backup check -> used_stratified_seed42_backup == True
            bcs_one = [b.rsplit("__", 1)[0] + "__3A" for b in barcodes]
            _write_10x_csc_h5(hdir / "IGVFFI5999NXAF.h5", sp.csc_matrix(mat), bcs_one, feat_ids, feat_names, feat_types)
            parts_bkp, audit_bkp = build_hcasmc_parts(raw_dir, hgnc, batch_cols=60)
            self.assertTrue(audit_bkp["used_stratified_seed42_backup"])
            self.assertGreater(audit_bkp["n_kept_half_a_ntc"], 0)
            self.assertGreater(audit_bkp["n_kept_half_b_ntc"], 0)
            self.assertEqual(set(np.unique(parts_bkp[0].half_flags).tolist()), {1, 2})

            # Q3 negative feasibility check on RNA-only HDF5 (missing guide capture rows)
            _write_10x_csc_h5(
                hdir / "IGVFFI5999NXAF.h5",
                sp.csc_matrix(mat[: len(genes) + 1, :]),
                barcodes,
                feat_ids[: len(genes) + 1],
                feat_names[: len(genes) + 1],
                feat_types[: len(genes) + 1],
            )
            feas_bad = check_hcasmc_feasibility(raw_dir, hgnc=hgnc, batch_cols=60)
            self.assertFalse(feas_bad["feasible"])

    def test_9_dataset_reader_cd4t_stim48_synthetic(self) -> None:
        rng = np.random.default_rng(901)
        hgnc = _get_hgnc()
        genes = _pick_hgnc_genes(520)
        feat_names = genes + ["MT-CO1", "NTC-1", f"g_{genes[0]}", f"g_{genes[1]}"]
        feat_ids = feat_names[:]
        feat_types = ["Gene Expression"] * (len(genes) + 1) + ["CRISPR Guide Capture"] * 3
        n_feat = len(feat_names)
        idx_ntc = len(genes) + 1
        idx_g0 = len(genes) + 2
        idx_g1 = len(genes) + 3

        def _make_donor_h5(path: Path, tag: str) -> None:
            bcs, assign = [], []
            for i in range(22):
                bcs.append(f"{tag}_N{i:03d}")
                assign.append((idx_ntc, -1))
            for i in range(20):
                bcs.append(f"{tag}_A{i:03d}")
                assign.append((idx_g0, 0))
            for i in range(20):
                bcs.append(f"{tag}_B{i:03d}")
                assign.append((idx_g1, 1))
            mat = np.zeros((n_feat, len(bcs)), dtype=np.float32)
            mat[: len(genes), :] = rng.integers(1, 5, size=(len(genes), len(bcs)), dtype=np.int32)
            mat[len(genes), :] = 5  # low MT (< 5%)
            for ci, (g_row, tgt_i) in enumerate(assign):
                if ci % 2 == 0:
                    mat[idx_ntc : idx_g1 + 1, ci] = 1
                mat[g_row, ci] = 25
                if tgt_i >= 0:
                    mat[:25, ci] += 8
                    mat[tgt_i, ci] = 1
            _write_10x_csc_h5(path, sp.csc_matrix(mat), bcs, feat_ids, feat_names, feat_types)

        with tempfile.TemporaryDirectory() as td:
            raw_dir = Path(td) / "raw"
            cdir = raw_dir / "cd4t_stim48"
            cdir.mkdir(parents=True)
            _make_donor_h5(cdir / "GSM9393801_D3_Stim48hr_CD4i_R2L25.h5", "D3")
            _make_donor_h5(cdir / "GSM9393802_D2_Stim48hr_CD4i_R2L26.h5", "D2")
            pd.DataFrame(
                [
                    {"sgRNA": "NTC-1", "target_gene_name": "", "target_gene_id": ""},
                    {"sgRNA": f"g_{genes[0]}", "target_gene_name": genes[0], "target_gene_id": hgnc.sym_to_ens.get(genes[0], "")},
                    {"sgRNA": f"g_{genes[1]}", "target_gene_name": genes[1], "target_gene_id": hgnc.sym_to_ens.get(genes[1], "")},
                ]
            ).to_csv(cdir / "sgrna_library_metadata.suppl_table.csv", index=False)

            # Full-lane guard (requires 24+24 unless allow_partial_lanes=True)
            with self.assertRaises(RuntimeError):
                build_cd4t_stim48_parts(raw_dir, hgnc, allow_partial_lanes=False)

            parts, audit = build_cd4t_stim48_parts(raw_dir, hgnc, allow_partial_lanes=True)
            self.assertEqual(len(parts), 2)
            self.assertEqual(set(parts[0].half_flags.tolist()), {1})
            self.assertEqual(set(parts[1].half_flags.tolist()), {2})

            out_root = Path(td) / "lines"
            process_r10_line("cd4t_stim48", parts, hgnc, out_root=out_root, extra_meta={"adapter_audit": audit}, n_workers=1)
            ds = load_standardized_line_r10(out_root / "cd4t_stim48", min_gene_overlap=50)
            self.assertEqual(len(ds.target_names), 2)

            # Q6 check: presence of any Rest / D1 / CE0006864 file MUST raise RuntimeError
            _make_donor_h5(cdir / "GSM9393700_D1_Rest_CD4i_CE0006864.h5", "BAD")
            with self.assertRaises(RuntimeError):
                build_cd4t_stim48_parts(raw_dir, hgnc, allow_partial_lanes=True)

    def test_10_dataset_reader_feng_ipsc_synthetic(self) -> None:
        rng = np.random.default_rng(1001)
        hgnc = _get_hgnc()
        genes = _pick_hgnc_genes(60)  # targeted gene panel (< 3,000 genes)
        tgts = genes[:3]

        meta_rows = []
        # 4 non-KOLF cell lines: line_a, line_c (index 0,2 -> Half A); line_b, line_d (index 1,3 -> Half B)
        for ln in ("line_a", "line_b", "line_c", "line_d"):
            # Main inlet with >= 15 NTCs + all 3 perturbations
            for i in range(18):
                meta_rows.append({"Cell_ID": f"{ln}_main_nt_{i}", "Cell_Line": ln, "Guide_Call": "NonTarget_1", "Target_Gene": "NonTarget", "Inlet": "1"})
            for g in tgts:
                for i in range(8):
                    meta_rows.append({"Cell_ID": f"{ln}_main_{g}_{i}", "Cell_Line": ln, "Guide_Call": f"{g}_g1", "Target_Gene": g, "Inlet": "1"})
            # Small inlet with ONLY 3 NTCs (< 15 -> triggers Cell_Line / half pooled NTC fallback!)
            for i in range(3):
                meta_rows.append({"Cell_ID": f"{ln}_small_nt_{i}", "Cell_Line": ln, "Guide_Call": "NonTarget_1", "Target_Gene": "NonTarget", "Inlet": "2"})
            for i in range(5):
                meta_rows.append({"Cell_ID": f"{ln}_small_{tgts[0]}_{i}", "Cell_Line": ln, "Guide_Call": f"{tgts[0]}_g1", "Target_Gene": tgts[0], "Inlet": "2"})

        # Excluded Q6 & P4 rows
        meta_rows.extend([
            {"Cell_ID": "kolf2_cell", "Cell_Line": "kolf_2", "Guide_Call": "NonTarget_1", "Target_Gene": "NonTarget", "Inlet": "1"},
            {"Cell_ID": "kolf3_cell", "Cell_Line": "kolf_3", "Guide_Call": f"{tgts[0]}_g1", "Target_Gene": tgts[0], "Inlet": "1"},
            {"Cell_ID": "unassigned_cell", "Cell_Line": "line_a", "Guide_Call": "unassigned", "Target_Gene": "", "Inlet": "1"},
            {"Cell_ID": "doublet_cell", "Cell_Line": "line_b", "Guide_Call": f"{tgts[0]}_g1,{tgts[1]}_g1", "Target_Gene": tgts[0], "Inlet": "1"},
        ])

        cell_ids = [r["Cell_ID"] for r in meta_rows]
        with tempfile.TemporaryDirectory() as td:
            raw_dir = Path(td) / "raw"
            fdir = raw_dir / "feng_ipsc"
            fdir.mkdir(parents=True)
            pd.DataFrame(meta_rows).to_csv(
                fdir / "TargetedScreen_Cell-Metadata.tsv.gz", sep="\t", index=False, compression="gzip"
            )

            # Transposed R `write.csv`-style header (leading blank `""` cell -> col_offset=1 auto-detection)
            # + skipped `:CRISPR-Guide-Capture` rows + planted consistent perturbation deltas
            with gzip.open(fdir / "TargetedScreen_RNA-UMI-Counts.csv.gz", "wt", newline="", encoding="utf-8") as gz:
                w = csv.writer(gz)
                w.writerow([""] + cell_ids)
                for gi, gsym in enumerate(genes):
                    ens = hgnc.sym_to_ens.get(gsym, f"ENSG{gi:011d}")
                    row_label = f"{ens}:{gsym}:Gene-Expression"
                    vals = rng.integers(2, 10, size=len(meta_rows), dtype=np.int32)
                    for ci, r in enumerate(meta_rows):
                        if r["Guide_Call"].endswith("_g1") and "," not in r["Guide_Call"]:
                            if gi < 18:
                                vals[ci] += 15
                            if gsym == r["Target_Gene"]:
                                vals[ci] = 1
                    w.writerow([row_label] + vals.tolist())
                # Guide capture rows that MUST be skipped cleanly during streaming
                w.writerow(["GUIDE_1:NonTarget_1:CRISPR-Guide-Capture"] + [5] * len(meta_rows))

            parts, audit = build_feng_ipsc_parts(raw_dir, hgnc, bucket_size=45)
            self.assertEqual(audit["n_kolf_cells_excluded"], 2)
            self.assertEqual(audit["n_unassigned_excluded"], 1)
            self.assertEqual(audit["n_multiguide_excluded"], 1)
            self.assertEqual(audit["n_guide_rows_skipped_in_csv"], 1)
            self.assertEqual(audit["n_non_kolf_lines"], 4)
            self.assertEqual(audit["half_a_lines"], ["line_a", "line_c"])
            self.assertEqual(audit["half_b_lines"], ["line_b", "line_d"])
            self.assertGreater(audit["n_low_ntc_unit_fallback_cells"], 0)

            out_root = Path(td) / "lines"
            process_r10_line("feng_ipsc", parts, hgnc, out_root=out_root, extra_meta={"adapter_audit": audit}, n_workers=1)
            ds_feng = load_standardized_line_r10(out_root / "feng_ipsc", min_gene_overlap=50)
            self.assertEqual(ds_feng.name, "feng_ipsc")
            self.assertEqual(len(ds_feng.genes), 60)
            self.assertEqual(len(ds_feng.target_names), 3)

            # Evaluate feng_ipsc on 60-gene panel against a matching source named `nourreddine_kolf21j`
            # (min_gene_overlap=50) and verify H-D5 output row is created cleanly
            src_mat = np.zeros_like(ds_feng.deltas)
            src_mat[:, :18] = 1.2
            for ti, tg_idx in enumerate(ds_feng.target_indices):
                src_mat[ti, tg_idx] = -1.5
            mu_kolf = np.exp(rng.uniform(-1.0, 2.5, size=len(ds_feng.genes)))
            c_kolf = rng.normal(0.0, 0.2, size=(24, len(ds_feng.genes))).astype(np.float32) + np.log1p(mu_kolf)[None, :].astype(np.float32)
            src_kolf = LineDatasetR10(
                name="nourreddine_kolf21j",
                genes=list(ds_feng.genes),
                target_names=list(ds_feng.target_names),
                target_indices=ds_feng.target_indices.copy(),
                deltas=src_mat,
                ctrl_cells=c_kolf,
                ctrl_mean_log1p=np.log1p(mu_kolf),
                ctrl_mean_libnorm=mu_kolf,
                ctrl_mean_libnorm_full=mu_kolf,
                repr_mask=np.ones(len(ds_feng.target_names), dtype=bool),
                half_a_deltas=src_mat.copy(),
                half_b_deltas=src_mat.copy(),
                half_valid=np.ones(len(ds_feng.target_names), dtype=bool),
                n_cells=np.full(len(ds_feng.target_names), 40.0),
                knockdown_depth_log1p=np.full(len(ds_feng.target_names), 1.5),
                knockdown_fraction=None,
                target_ensembl_ids=None,
                representation="log1p",
                metadata={},
            )
            _, b_rows, _ = evaluate_target_against_sources_r11(
                target_ds=ds_feng,
                source_datasets={"nourreddine_kolf21j": src_kolf},
                only_reproducible=False,
                include_r9_anchors=False,
                include_exploratory_e=False,
                min_gene_overlap=50,
                n_boot=60,
                seed=20261123,
                n_perm=5,
            )
            comps = {r["comparison"] for r in b_rows}
            self.assertIn("H-D5_FCT_minus_Mis_kolf21j_all_nontarget", comps)


if __name__ == "__main__":
    unittest.main(verbosity=2)
