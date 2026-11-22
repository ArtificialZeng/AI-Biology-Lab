#!/usr/bin/env python3
"""
Unit & Adversarial Verification Suite for Independent Round 11 Audit (`round11/audit/`).

Strictly tests ONLY synthetic fixtures and never accesses `round11/src/`, `round11/results/`,
or any unseen Round 11 dataset.
"""

from __future__ import annotations

import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

AUDIT_DIR = Path(__file__).resolve().parent
if str(AUDIT_DIR) not in sys.path:
    sys.path.insert(0, str(AUDIT_DIR))

from independent_r11_loader import (  # noqa: E402
    AuditLineData,
    enforce_r11_gate,
    load_line_dataset,
    parse_confirmatory_lines_r11,
)
from independent_r11_metrics import (  # noqa: E402
    DEFAULT_R11_SEED,
    adjudicate_round11_group_d,
    ceiling_single_half,
    disattenuate,
    evaluate_cross_mask_cal_diff,
    evaluate_fct_perm_null,
    evaluate_joint_residuals,
    evaluate_paired_cal_diff,
    evaluate_sa8_perm_null,
    evaluate_sa8_resp_top50_split,
    orthogonalize_to_ref,
    pearson_vec,
)
from run_r11_audit import evaluate_target_line_r11  # noqa: E402


def _make_synthetic_cohort(
    n_genes: int = 120,
    n_targets: int = 28,
    noise_sd: float = 0.04,
    seed: int = 77,
    noisy_split_frac: float = 0.0,
):
    """Build synthetic R9-like sources + unseen target exhibiting nonlinear log1p baseline coupling."""
    rng = np.random.default_rng(seed)
    genes = [f"G{i}" for i in range(n_genes)]
    targets = genes[:n_targets]

    # Realistic skewed libnorm control expression across target and 6 sources
    tgt_mu = rng.lognormal(mean=-0.3, sigma=1.1, size=n_genes)
    tgt_ctrl_cells = np.log1p(tgt_mu)[None, :] + rng.normal(0.0, 0.15, size=(60, n_genes))

    src_names = (
        "nadig_hepg2",
        "nadig_jurkat",
        "nourreddine_kolf21j",
        "replogle_k562_essential",
        "replogle_rpe1",
        "vcc2025_h1",
    )
    true_fc = np.exp(rng.normal(0.0, 0.45, size=(n_targets, n_genes)))

    sources = []
    for s_idx, s_name in enumerate(src_names):
        s_mu = rng.lognormal(mean=-0.35 + 0.1 * s_idx, sigma=1.0, size=n_genes)
        s_mu = np.maximum(s_mu, 1e-4)
        s_delta = np.log1p(s_mu[None, :] * true_fc) - np.log1p(s_mu)[None, :]
        s_delta += rng.normal(0.0, 0.02, size=s_delta.shape)
        s_ctrl = np.log1p(s_mu)[None, :] + rng.normal(0.0, 0.15, size=(50, n_genes))
        sources.append(
            AuditLineData(
                name=s_name,
                genes=genes,
                target_names=targets,
                target_indices=np.arange(n_targets, dtype=np.int64),
                deltas=s_delta,
                ctrl_cells=s_ctrl,
                ctrl_mean_libnorm=s_mu,
                ctrl_mean_libnorm_full=s_mu,
                repr_mask=np.ones(n_targets, dtype=bool),
                half_a_deltas=None,
                half_b_deltas=None,
                half_valid=None,
                n_cells=np.full(n_targets, 100.0),
                kd_depth_log1p=np.full(n_targets, 0.4),
                metadata={},
            )
        )

    signal = np.log1p(tgt_mu[None, :] * true_fc) - np.log1p(tgt_mu)[None, :]
    noise_a = rng.normal(0.0, noise_sd, size=signal.shape)
    noise_b = rng.normal(0.0, noise_sd, size=signal.shape)
    if noisy_split_frac > 0.0:
        # Inject anti-correlated split-half noise into a subset of rows so
        # point rho_ab_orth > 0 and CI975_lo > 0, while ~20% of bootstrap draws cross rho <= 0.
        n_bad = max(1, int(round(noisy_split_frac * n_targets)))
        flip = rng.normal(0.0, 0.13, size=(n_bad, n_genes))
        noise_a[:n_bad] += flip
        noise_b[:n_bad] -= flip
    half_a = signal + noise_a
    half_b = signal + noise_b
    delta_obs = 0.5 * (half_a + half_b)

    tgt_data = AuditLineData(
        name="synth_conf_line",
        genes=genes,
        target_names=targets,
        target_indices=np.arange(n_targets, dtype=np.int64),
        deltas=delta_obs,
        ctrl_cells=tgt_ctrl_cells,
        ctrl_mean_libnorm=tgt_mu,
        ctrl_mean_libnorm_full=tgt_mu,
        repr_mask=np.ones(n_targets, dtype=bool),
        half_a_deltas=half_a,
        half_b_deltas=half_b,
        half_valid=np.ones(n_targets, dtype=bool),
        n_cells=np.linspace(40.0, 200.0, n_targets),
        kd_depth_log1p=np.linspace(0.2, 0.8, n_targets),
        metadata={},
    )
    return tgt_data, sources


class TestRound11Audit(unittest.TestCase):
    def test_blinding_gate_and_confirmatory_parser(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sdir = Path(tmp) / "status"
            sdir.mkdir()
            with self.assertRaises(RuntimeError):
                enforce_r11_gate(["telohaec"], dev_only=False, status_dir=sdir)
            (sdir / "PREREG_D_COMMIT.txt").write_text("abc123\n", encoding="utf-8")
            with self.assertRaises(RuntimeError):
                enforce_r11_gate(["telohaec"], dev_only=False, status_dir=sdir)
            (sdir / "CODE_FREEZE_R11.txt").write_text("def456\n", encoding="utf-8")
            h1, h2 = enforce_r11_gate(["telohaec"], dev_only=False, status_dir=sdir)
            self.assertEqual((h1, h2), ("abc123", "def456"))

            # dev_only refuses Round 11 lines and round11/shared/ directories
            with self.assertRaises(RuntimeError):
                enforce_r11_gate(["telohaec"], dev_only=True, status_dir=sdir)
            with self.assertRaises(RuntimeError):
                enforce_r11_gate(["feng_ipsc"], dev_only=True, status_dir=sdir)
            with self.assertRaises(RuntimeError):
                enforce_r11_gate(
                    ["dld1"],
                    dev_only=True,
                    status_dir=sdir,
                    extra_lines_dirs=[Path(__file__).resolve().parents[1] / "shared" / "lines" / "telohaec"],
                )
            # parse_confirmatory_lines_r11 ignores feng_ipsc even if present
            cf = sdir / "CONFIRMATORY_LINES_R11.txt"
            cf.write_text("# comment\ntelohaec\nfeng_ipsc\ncd4t_stim48\n", encoding="utf-8")
            self.assertEqual(parse_confirmatory_lines_r11(cf), ["telohaec", "cd4t_stim48"])

    def test_planted_primary_and_naive_elementwise_parity(self) -> None:
        tgt_data, sources = _make_synthetic_cohort(n_genes=120, n_targets=30, noise_sd=0.03, seed=101)
        res = evaluate_target_line_r11(
            target_data=tgt_data,
            source_list=sources,
            n_boot=200,
            n_perm=20,
            seed=DEFAULT_R11_SEED,
            include_sa3_anchors=True,
            include_r19=True,
            include_s7=True,
        )
        fa = res["method_eval"]["fct"]["all_nontarget"]
        self.assertTrue(res["hd1_pass"])
        self.assertTrue(res["hd2_pass"])
        self.assertTrue(res["hd3_pass"])
        self.assertTrue(res["hd_util_pass"])
        self.assertGreater(fa["r_orth_cal_valid_frac"], 0.99)
        self.assertGreater(fa["rho_ab_orth"], 0.0)

        # Dual-mask permutations computed cleanly
        self.assertEqual(res["perm100_all_nontarget"]["mode"], "all_nontarget")
        self.assertEqual(res["perm100_resp_top50"]["mode"], "resp_top50")
        self.assertLess(res["perm100_all_nontarget"]["r_orth_cal"]["p975"], fa["r_orth_cal_ci975_lo"])

        # Naive elementwise check for all_nontarget point estimate (zero reliance on _moments)
        t_mean = res["method_eval"]["transferred_mean"]  # triggers existence check
        self.assertIsNotNone(t_mean)
        from independent_r11_predictors import align_source_and_target, compute_source_fc_bar, fct_from_fc_bar

        cohort = align_source_and_target(sources, tgt_data, only_reproducible=True)
        fc_bar, tm, _ = compute_source_fc_bar(
            cohort.eval_targets, cohort.source_names, cohort.source_deltas_by_line, cohort.source_ctrl_mean_libnorm
        )
        pred_fct = fct_from_fc_bar(fc_bar, cohort.tgt_ctrl_mean_libnorm)
        x_list, y_list, ha_list, hb_list = [], [], [], []
        for p in range(len(cohort.eval_targets)):
            idx = [g for g in range(len(cohort.eval_genes)) if g != int(cohort.target_cols[p])]
            ref = tm[p, idx]
            x_list.append(orthogonalize_to_ref(cohort.delta_obs[p, idx] - ref, ref))
            y_list.append(orthogonalize_to_ref(pred_fct[p, idx] - ref, ref))
            ha_list.append(orthogonalize_to_ref(cohort.half_a_obs[p, idx] - ref, ref))
            hb_list.append(orthogonalize_to_ref(cohort.half_b_obs[p, idx] - ref, ref))
        r_raw_naive = pearson_vec(np.concatenate(x_list), np.concatenate(y_list))
        rho_ab_naive = pearson_vec(np.concatenate(ha_list), np.concatenate(hb_list))
        ceil_naive = math.sqrt((2.0 * rho_ab_naive) / (1.0 + rho_ab_naive))
        r_cal_naive = disattenuate(r_raw_naive, ceil_naive)
        self.assertAlmostEqual(fa["r_orth_raw_est"], r_raw_naive, places=12)
        self.assertAlmostEqual(fa["rho_ab_orth"], rho_ab_naive, places=12)
        self.assertAlmostEqual(fa["r_orth_cal_est"], r_cal_naive, places=12)

    def test_addendum_a1_validity_and_adjudication(self) -> None:
        # Plant a line with CI_lo > 0 but boot_valid_frac < 0.95 via anti-correlated halves
        tgt_data, sources = _make_synthetic_cohort(
            n_genes=100, n_targets=18, noise_sd=0.015, seed=202, noisy_split_frac=0.45
        )
        from independent_r11_predictors import align_source_and_target, compute_source_fc_bar, fct_from_fc_bar

        cohort = align_source_and_target(sources, tgt_data, only_reproducible=True)
        fc_bar, tm, _ = compute_source_fc_bar(
            cohort.eval_targets, cohort.source_names, cohort.source_deltas_by_line, cohort.source_ctrl_mean_libnorm
        )
        pred_fct = fct_from_fc_bar(fc_bar, cohort.tgt_ctrl_mean_libnorm)
        ev = evaluate_joint_residuals(
            cohort.delta_obs,
            pred_fct,
            tm,
            cohort.target_cols,
            mode="all_nontarget",
            n_boot=200,
            seed=DEFAULT_R11_SEED,
            half_a_obs=cohort.half_a_obs,
            half_b_obs=cohort.half_b_obs,
            half_valid=cohort.half_valid,
        )
        self.assertLess(ev["r_orth_cal_valid_frac"], 0.95)
        self.assertGreater(ev["r_orth_cal_ci975_lo"], 0.0)

        res = evaluate_target_line_r11(
            tgt_data, sources, n_boot=200, n_perm=10, seed=DEFAULT_R11_SEED, include_sa3_anchors=False, include_r19=False, include_s7=False
        )
        self.assertTrue(res["hd1_ci975_gt_0_only"])
        self.assertFalse(res["hd1_pass"])
        self.assertTrue(res["hd1_a1_deviation_flag"])

        # Adjudicate P1 (< 2), deviation flagging, & independent Claim 1 vs Claim 2
        p1 = adjudicate_round11_group_d({"telohaec": res}, ["telohaec"])
        self.assertEqual(p1["status"], "NO_CONFIRMATORY_ADJUDICATION")
        self.assertIn("确认性线不足", p1["abstract_notice"])

        good = {
            "hd1_pass": True, "hd1_ci975_gt_0_only": True, "hd1_a1_deviation_flag": False,
            "hd2_pass": True, "hd3_pass": True, "hd4a_pass": False, "hd4b_pass": True, "hd_util_pass": True,
        }
        adj = adjudicate_round11_group_d({"telohaec": good, "hcasmc": res}, ["telohaec", "hcasmc"])
        self.assertEqual(adj["hd1_pass_count"], "1/2")
        self.assertEqual(adj["hd1_ci975_gt_0_only_count"], "2/2")
        self.assertEqual(adj["hd1_a1_deviation_lines"], ["hcasmc"])
        self.assertEqual(adj["primary_verdict"], "未确认 (NOT CONFIRMED)")
        self.assertEqual(adj["hd4_verdict"], "未确认 (NOT CONFIRMED)")

        # Independent Claim 1 PASS while Claim 2 FAIL, and vice versa
        adj_c1_only = adjudicate_round11_group_d(
            {"telohaec": good, "cd4t_stim48": dict(good)}, ["telohaec", "cd4t_stim48"]
        )
        self.assertEqual(adj_c1_only["primary_verdict"], "已确认 (CONFIRMED)")
        self.assertEqual(adj_c1_only["hd4_verdict"], "未确认 (NOT CONFIRMED)")
        self.assertTrue(adj_c1_only["cd4t_stim48_independence_required"])

    def test_sa8_split_and_nicol_artifact_removal(self) -> None:
        # Pure-noise split halves around zero true signal + nonlinear FC_bar shows how
        # standard resp_top50 inflates E1_ConstMed vs SA-8 collapsing back near zero,
        # AND checks exact elementwise parity of SA-8 A->B and B->A.
        rng = np.random.default_rng(555)
        n_t, n_g = 24, 150
        mu_t = np.exp(rng.normal(0.0, 1.0, size=n_g))
        fc_bar = np.exp(rng.normal(0.0, 0.35, size=(n_t, n_g)))
        t_mean = np.log1p(mu_t[None, :] * 1.0) * 0.0 + rng.normal(0.0, 0.1, size=(n_t, n_g))
        e1_pred = np.log1p(float(np.median(mu_t)) * fc_bar) - math.log1p(float(np.median(mu_t)))
        t_cols = np.arange(n_t, dtype=np.int64)

        # Shared finite signal in one half + correlated amplitude shift to yield positive rho_ab
        shared_sig = 0.25 * (np.log1p(mu_t[None, :] * fc_bar) - np.log1p(mu_t)[None, :])
        half_a = shared_sig + rng.normal(0.0, 0.08, size=(n_t, n_g))
        half_b = shared_sig + rng.normal(0.0, 0.08, size=(n_t, n_g))
        delta_obs = 0.5 * (half_a + half_b)

        sa8 = evaluate_sa8_resp_top50_split(
            delta_obs=delta_obs,
            pred=e1_pred,
            t_mean=t_mean,
            target_cols=t_cols,
            half_a_obs=half_a,
            half_b_obs=half_b,
            top_k=30,
            n_boot=200,
            seed=DEFAULT_R11_SEED,
        )
        # Elementwise manual verification of A->B, B->A single-half ceiling & average
        xa_l, ya_l, ha1_l, hb1_l = [], [], [], []
        xb_l, yb_l, ha2_l, hb2_l = [], [], [], []
        for p in range(n_t):
            avail = np.asarray([g for g in range(n_g) if g != p], dtype=np.int64)
            ia = avail[np.argsort(np.abs(half_a[p, avail]), kind="mergesort")[-30:]]
            ib = avail[np.argsort(np.abs(half_b[p, avail]), kind="mergesort")[-30:]]
            ra = t_mean[p, ia]
            xa_l.append(orthogonalize_to_ref(half_b[p, ia] - ra, ra))
            ya_l.append(orthogonalize_to_ref(e1_pred[p, ia] - ra, ra))
            ha1_l.append(orthogonalize_to_ref(half_a[p, ia] - ra, ra))
            hb1_l.append(orthogonalize_to_ref(half_b[p, ia] - ra, ra))

            rb = t_mean[p, ib]
            xb_l.append(orthogonalize_to_ref(half_a[p, ib] - rb, rb))
            yb_l.append(orthogonalize_to_ref(e1_pred[p, ib] - rb, rb))
            ha2_l.append(orthogonalize_to_ref(half_a[p, ib] - rb, rb))
            hb2_l.append(orthogonalize_to_ref(half_b[p, ib] - rb, rb))

        r_ab = pearson_vec(np.concatenate(xa_l), np.concatenate(ya_l))
        r_ba = pearson_vec(np.concatenate(xb_l), np.concatenate(yb_l))
        c_ab = ceiling_single_half(pearson_vec(np.concatenate(ha1_l), np.concatenate(hb1_l)))
        c_ba = ceiling_single_half(pearson_vec(np.concatenate(ha2_l), np.concatenate(hb2_l)))
        self.assertAlmostEqual(sa8["r_orth_raw_split_est"], 0.5 * (r_ab + r_ba), places=12)
        self.assertAlmostEqual(sa8["r_orth_cal_split_est"], 0.5 * (disattenuate(r_ab, c_ab) + disattenuate(r_ba, c_ba)), places=12)

        sa8_perm = evaluate_sa8_perm_null(
            delta_obs=delta_obs,
            t_mean=t_mean,
            fc_bar=fc_bar,
            tgt_ctrl_mean_libnorm=mu_t,
            target_cols=t_cols,
            half_a_obs=half_a,
            half_b_obs=half_b,
            ceil_single_ab=sa8["ceil_single_half_AtoB"],
            ceil_single_ba=sa8["ceil_single_half_BtoA"],
            top_k=30,
            n_perm=15,
        )
        self.assertTrue(np.isfinite(sa8_perm["r_orth_cal_split"]["median"]))
        self.assertTrue(np.isfinite(sa8_perm["r_orth_cal_split"]["p975"]))

        # Cross-mask E1 diff runs cleanly
        cm = evaluate_cross_mask_cal_diff(
            delta_obs, e1_pred, t_mean, t_cols,
            mode_a="resp_top50", mode_b="all_nontarget",
            n_boot=100, seed=DEFAULT_R11_SEED,
            half_a_obs=half_a, half_b_obs=half_b,
        )
        self.assertTrue(np.isfinite(cm["diff_est"]))
        self.assertGreater(cm["diff_valid_frac"], 0.9)

    def test_hd5_targeted_panel_overlap_and_kolf21j_comparison(self) -> None:
        # Verify targeted panel (< 3000 genes) fails at min_gene_overlap=3000 & passes at 100
        # and verifies H-D5 output block when `nourreddine_kolf21j` is present.
        with tempfile.TemporaryDirectory() as tmp:
            ldir = Path(tmp) / "feng_ipsc"
            ldir.mkdir()
            n_g, n_t = 320, 20
            genes = [f"G{i}" for i in range(n_g)]
            targets = genes[:n_t]
            (ldir / "genes.txt").write_text("\n".join(genes) + "\n", encoding="utf-8")
            (ldir / "targets.txt").write_text("\n".join(targets) + "\n", encoding="utf-8")
            (ldir / "line_metadata.json").write_text(json.dumps({"role": "secondary"}), encoding="utf-8")

            tgt_data, sources = _make_synthetic_cohort(n_genes=n_g, n_targets=n_t, noise_sd=0.03, seed=303)
            np.savez(
                ldir / "pseudobulk_delta.npz",
                targets=np.array(targets),
                delta_log1p=tgt_data.deltas,
                reproducible=np.ones(n_t, dtype=bool),
                half_a_delta_log1p=tgt_data.half_a_deltas,
                half_b_delta_log1p=tgt_data.half_b_deltas,
                half_valid=np.ones(n_t, dtype=bool),
                n_cells=tgt_data.n_cells,
                knockdown_depth_log1p=tgt_data.kd_depth_log1p,
            )
            np.savez(
                ldir / "ctrl_cells.npz",
                X_log1p=tgt_data.ctrl_cells,
                ctrl_mean_libnorm=tgt_data.ctrl_mean_libnorm,
            )

            with self.assertRaises(RuntimeError):
                load_line_dataset(ldir, aligned_genes=genes, hgnc=None, min_gene_overlap=3000)
            loaded = load_line_dataset(ldir, aligned_genes=genes, hgnc=None, min_gene_overlap=100)
            res = evaluate_target_line_r11(
                loaded, sources, n_boot=150, n_perm=10, seed=DEFAULT_R11_SEED, include_sa3_anchors=False, include_r19=False, include_s7=False
            )
            hd5 = res["hd5_kolf21j"]
            self.assertIsNotNone(hd5)
            self.assertTrue(hd5["is_feng_ipsc"])
            self.assertTrue(hd5["fct_ci975_lo_gt_kolf21j_ci975_hi"])
            self.assertTrue(hd5["paired_diff_ci975_gt_0"])


if __name__ == "__main__":
    unittest.main()
