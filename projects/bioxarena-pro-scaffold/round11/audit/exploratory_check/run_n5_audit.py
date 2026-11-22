#!/usr/bin/env python3
"""
Exploratory Check 2: Independent Implementation & Basis Audit of Null Cone `N_5(p)`
and Noise-Calibrated Semi-Partial Correlation `r_cone_sp`.

Strictly follows `README.md` and `THEOREM_AND_PROOFS.md` definitions:
  - N_5(p) = span{1, t_mean[p], b_S, g_{0.05}(FC_bar[p]), g_{1.0}(FC_bar[p]), g_{20.0}(FC_bar[p])}
  - g_c(FC) = log(1 + c * FC) - log(1 + c)
  - r_cone_sp(Delta_obs, pred) =
        (1 / c_{N5}) * (<P_{N5}^perp r_obs, P_{N5}^perp r_pred>) / (||P_{N5}^perp r_obs|| * ||(r_pred)_cent||)
    = eta(r_pred) * r_cone_orth_cal(r_obs, r_pred)

Explicitly audits:
  1. Construction-level nulling (`eta` and `R^2_{N5}`): distinguishes exact algebraic
     zeroing (`Transfer_Mean_GlobalScale`, `Source_Grand_Mean`: eta == 0 by construction)
     and convex-cone spanning (`E1`, `E2_P10`, `E2_P90`, `E6_LogFCBar`, `E7_FCBarMinus1`:
     R^2_{N5} > 0.997) from genuine non-cone validation tests.
  2. Single-draw gene permutation distribution (`Perm0..99`, NOT expected mean E_pi).
  3. All 6 mismatched baselines `Mis_s` across all 12 lines.
  4. Spotlight paired comparison: `FCT` vs `Mis(vcc2025_h1)` on `cd4t_rest` and `cd4t_stim48`
     (both unpaired CI975_lo vs CI975_hi AND joint perturbation bootstrap difference `97.5% CI`).
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

EXPLORATORY_CHECK_DIR = Path(__file__).resolve().parent
AUDIT_DIR = EXPLORATORY_CHECK_DIR.parent
for p in (str(EXPLORATORY_CHECK_DIR), str(AUDIT_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from exploratory_core import (  # noqa: E402
    DEFAULT_AUDIT_SEED,
    LINES_12,
    build_r9_lolo_cohort,
    load_12_lines,
)
from independent_r11_metrics import (  # noqa: E402
    EPS,
    ceiling_from_halves,
    disattenuate,
    summarize_boot,
)
from independent_r11_predictors import (  # noqa: E402
    build_predictions,
    compute_source_fc_bar,
    fct_from_fc_bar,
)

CONE_NODES: Tuple[float, ...] = (0.05, 1.0, 20.0)


def g_const(fc: np.ndarray, c: float) -> np.ndarray:
    """Normalized log1p response curve g_c(FC) = log(1 + c * FC) - log(1 + c)."""
    return np.log1p(c * np.asarray(fc, dtype=np.float64)) - np.log1p(c)


def build_n5_basis(
    t_mean_p: np.ndarray,
    grand_mean_p: np.ndarray,
    fc_bar_p: np.ndarray,
    rcond: float = 1e-12,
) -> np.ndarray:
    """
    Orthonormal SVD basis Q_p (m_p x r_p) for:
        N_5(p) = span{1, t_mean_p, b_S, g_{0.05}(fc_bar_p), g_{1.0}(fc_bar_p), g_{20.0}(fc_bar_p)}.
    Intercept `1` is ALWAYS included as first unit basis column so any vector
    projected by `v - Q_p @ (Q_p.T @ v)` is guaranteed zero-mean per perturbation.
    """
    m = t_mean_p.shape[0]
    ones = np.full(m, 1.0 / np.sqrt(float(m)), dtype=np.float64)
    cols_raw = [
        t_mean_p.astype(np.float64),
        grand_mean_p.astype(np.float64),
        g_const(fc_bar_p, CONE_NODES[0]),
        g_const(fc_bar_p, CONE_NODES[1]),
        g_const(fc_bar_p, CONE_NODES[2]),
    ]
    # Remove intercept from non-constant columns first, then orthogonalize via SVD
    orth_cols = []
    for v in cols_raw:
        vp = v - float(np.dot(ones, v)) * ones
        nrm = float(np.linalg.norm(vp))
        if nrm > EPS:
            orth_cols.append(vp / nrm)
    if not orth_cols:
        return ones[:, None]
    mat = np.column_stack(orth_cols)
    u, sv, _ = np.linalg.svd(mat, full_matrices=False)
    keep = sv > max(EPS, rcond * float(sv[0]))
    if not np.any(keep):
        return ones[:, None]
    return np.column_stack([ones, u[:, keep]])


def prepare_n5_per_pert_data(
    delta_obs: np.ndarray,
    t_mean: np.ndarray,
    grand_mean: np.ndarray,
    fc_bar: np.ndarray,
    target_cols: Sequence[int],
    half_a_obs: Optional[np.ndarray],
    half_b_obs: Optional[np.ndarray],
    half_valid: Optional[np.ndarray],
    mode: str = "all_nontarget",
    top_k: int = 50,
) -> Tuple[List[Tuple[int, np.ndarray, np.ndarray]], np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Precompute per-perturbation N_5 bases `Q_p`, observed projected residuals
    `P_{N5}^perp r_obs`, and split-half sufficient moments for fast joint bootstrap.
    """
    n_t, n_g = delta_obs.shape
    per_p: List[Tuple[int, np.ndarray, np.ndarray]] = []
    ro_ss: List[float] = []
    haa_ss: List[float] = []
    hbb_ss: List[float] = []
    hab_ss: List[float] = []

    has_halves = half_a_obs is not None and half_b_obs is not None
    ranks: List[int] = []

    for p in range(n_t):
        valid = np.isfinite(delta_obs[p]) & np.isfinite(t_mean[p]) & np.isfinite(fc_bar[p]) & np.isfinite(grand_mean)
        tc = int(target_cols[p])
        if 0 <= tc < n_g:
            valid[tc] = False
        idx = np.flatnonzero(valid)
        if idx.size == 0:
            continue
        if mode == "resp_top50" and idx.size > top_k:
            idx = idx[np.argsort(np.abs(delta_obs[p, idx]))[-top_k:]]

        tm_p = t_mean[p, idx].astype(np.float64)
        gm_p = grand_mean[idx].astype(np.float64)
        fcb_p = fc_bar[p, idx].astype(np.float64)
        Q_p = build_n5_basis(tm_p, gm_p, fcb_p)
        ranks.append(Q_p.shape[1])

        r_obs = delta_obs[p, idx].astype(np.float64) - tm_p
        r_obs_orth = r_obs - Q_p @ (Q_p.T @ r_obs)
        per_p.append((p, idx, Q_p))
        ro_ss.append(float(np.dot(r_obs_orth, r_obs_orth)))

        if has_halves and (half_valid is None or bool(half_valid[p])):
            ha = half_a_obs[p, idx].astype(np.float64)
            hb = half_b_obs[p, idx].astype(np.float64)
            if np.all(np.isfinite(ha) & np.isfinite(hb)):
                ra = ha - tm_p
                rb = hb - tm_p
                ra_o = ra - Q_p @ (Q_p.T @ ra)
                rb_o = rb - Q_p @ (Q_p.T @ rb)
                haa_ss.append(float(np.dot(ra_o, ra_o)))
                hbb_ss.append(float(np.dot(rb_o, rb_o)))
                hab_ss.append(float(np.dot(ra_o, rb_o)))
                continue
        haa_ss.append(0.0)
        hbb_ss.append(0.0)
        hab_ss.append(0.0)

    meta = {
        "n_eval_perts": len(per_p),
        "rank_min": min(ranks) if ranks else 0,
        "rank_max": max(ranks) if ranks else 0,
        "rank_mean": float(np.mean(ranks)) if ranks else 0.0,
    }
    return (
        per_p,
        np.asarray(ro_ss, dtype=np.float64),
        np.asarray(haa_ss, dtype=np.float64),
        np.asarray(hbb_ss, dtype=np.float64),
        np.asarray(hab_ss, dtype=np.float64),
        meta,
    )


def _pred_n5_moments(
    per_p: List[Tuple[int, np.ndarray, np.ndarray]],
    delta_obs: np.ndarray,
    pred: np.ndarray,
    t_mean: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-perturbation <r_obs_orth, r_pred_orth>, ||r_pred_orth||^2, ||r_pred_cent||^2."""
    m = len(per_p)
    xy = np.empty(m, dtype=np.float64)
    yy_o = np.empty(m, dtype=np.float64)
    yy_c = np.empty(m, dtype=np.float64)
    for k, (p, idx, Q_p) in enumerate(per_p):
        tm_p = t_mean[p, idx].astype(np.float64)
        r_obs = delta_obs[p, idx].astype(np.float64) - tm_p
        r_obs_o = r_obs - Q_p @ (Q_p.T @ r_obs)
        r_pr = pred[p, idx].astype(np.float64) - tm_p
        r_pr_c = r_pr - float(np.mean(r_pr))
        r_pr_o = r_pr - Q_p @ (Q_p.T @ r_pr)
        xy[k] = float(np.dot(r_obs_o, r_pr_o))
        yy_o[k] = float(np.dot(r_pr_o, r_pr_o))
        yy_c[k] = float(np.dot(r_pr_c, r_pr_c))
    return xy, yy_o, yy_c


def _eval_n5_from_moments(
    xy: np.ndarray,
    yy_o: np.ndarray,
    yy_c: np.ndarray,
    xx_o: np.ndarray,
    haa_o: np.ndarray,
    hbb_o: np.ndarray,
    hab_o: np.ndarray,
    draws: Optional[np.ndarray],
    n_boot: int,
) -> Dict[str, Any]:
    """Compute r_cone_sp, r_cone_sp_raw, r_cone_orth_cal, eta, and R2_N5 + bootstrap CIs."""
    m = xy.size

    def _calc(
        s_xy: Any, s_xx: Any, s_yyc: Any, s_yyo: Any,
        s_ha: Any, s_hb: Any, s_hab: Any,
    ) -> Tuple[Any, Any, Any, Any, Any, Any]:
        eta = np.where(s_yyc > EPS, np.sqrt(np.clip(s_yyo / np.maximum(s_yyc, EPS), 0.0, 1.0)), 0.0)
        denom_sp = np.sqrt(np.maximum(s_xx, 0.0) * np.maximum(s_yyc, 0.0))
        r_sp_raw = np.where((denom_sp > EPS) & (s_yyo > 1e-24), np.clip(s_xy / np.maximum(denom_sp, EPS), -1.0, 1.0), 0.0)
        denom_orth = np.sqrt(np.maximum(s_xx, 0.0) * np.maximum(s_yyo, 0.0))
        r_orth_raw = np.where((denom_orth > EPS) & (s_yyo > 1e-24), np.clip(s_xy / np.maximum(denom_orth, EPS), -1.0, 1.0), 0.0)

        hden = np.sqrt(np.maximum(s_ha, 0.0) * np.maximum(s_hb, 0.0))
        rho = np.where(hden > EPS, np.clip(s_hab / np.maximum(hden, EPS), -1.0, 1.0), 0.0)
        ceil = np.where(rho > 0.0, np.sqrt((2.0 * rho) / (1.0 + rho)), 0.0)
        r_sp_cal = np.where(ceil > EPS, r_sp_raw / np.maximum(ceil, EPS), float("nan"))
        return r_sp_cal, r_sp_raw, r_orth_raw, eta, rho, ceil

    pt_cal, pt_raw, pt_orth_raw, pt_eta, pt_rho, pt_ceil = _calc(
        float(np.sum(xy)), float(np.sum(xx_o)), float(np.sum(yy_c)), float(np.sum(yy_o)),
        float(np.sum(haa_o)), float(np.sum(hbb_o)), float(np.sum(hab_o)),
    )
    pt_r2_n5 = float(max(0.0, 1.0 - float(pt_eta) ** 2)) if np.sum(yy_c) > EPS else 1.0

    if draws is not None and m > 0:
        b_cal, b_raw, _, b_eta, _, _ = _calc(
            np.sum(xy[draws], axis=1),
            np.sum(xx_o[draws], axis=1),
            np.sum(yy_c[draws], axis=1),
            np.sum(yy_o[draws], axis=1),
            np.sum(haa_o[draws], axis=1),
            np.sum(hbb_o[draws], axis=1),
            np.sum(hab_o[draws], axis=1),
        )
    else:
        b_cal = np.array([float(pt_cal)])
        b_raw = np.array([float(pt_raw)])
        b_eta = np.array([float(pt_eta)])

    out: Dict[str, Any] = {
        "n_eval_perts": m,
        "rho_ab_n5": float(pt_rho),
        "ceiling_n5": float(pt_ceil),
        "r_cone_orth_raw": float(pt_orth_raw),
        "r2_in_n5": pt_r2_n5,
        "_b_cal": b_cal,  # dropped before json serialization
    }
    out.update(summarize_boot(float(pt_cal), b_cal, m, n_boot).to_dict("r_cone_sp"))
    out.update(summarize_boot(float(pt_raw), b_raw, m, n_boot).to_dict("r_cone_sp_raw"))
    out.update(summarize_boot(float(pt_eta), b_eta, m, n_boot).to_dict("eta"))
    return out


def run_check2(
    seed: int = DEFAULT_AUDIT_SEED,
    n_boot: int = 4000,
    n_perm: int = 100,
    out_dir: Path = EXPLORATORY_CHECK_DIR / "outputs",
) -> Dict[str, Any]:
    t0 = time.time()
    out_dir.mkdir(parents=True, exist_ok=True)
    lines = load_12_lines()

    method_rows: List[Dict[str, Any]] = []
    spotlight_rows: List[Dict[str, Any]] = []
    perm_rows: List[Dict[str, Any]] = []
    per_line: Dict[str, Any] = {}

    construction_classes = {
        "transfer_mean_global_scale": "ZERO_BY_EXACT_CONSTRUCTION (r_pred in span{t_mean} subset N_5)",
        "source_grand_mean": "ZERO_BY_EXACT_CONSTRUCTION (r_pred = b_S - t_mean in N_5)",
        "e1_const_med": "SUPPRESSED_BY_CONE_CONSTRUCTION (g_{c_med} in conv{g_0.05, g_1.0, g_20.0})",
        "e2_const_p10": "SUPPRESSED_BY_CONE_CONSTRUCTION (g_{p10} in conv{g_0.05, g_1.0, g_20.0})",
        "e2_const_p90": "SUPPRESSED_BY_CONE_CONSTRUCTION (g_{p90} in conv{g_0.05, g_1.0, g_20.0})",
        "e6_log_fc_bar": "SUPPRESSED_BY_CONE_CONSTRUCTION (c -> inf boundary ray of N_5 cone)",
        "e7_fc_minus_1": "SUPPRESSED_BY_CONE_CONSTRUCTION (c -> 0+ boundary ray of N_5 cone)",
        "e3_src_mean_ctrl": "SOURCE_BASELINE_REFERENCE (mu_S weighted)",
        "e5_rank_matched_mu": "TARGET_RANK_BASELINE (quantile-matched mu_t)",
        "fct": "INFORMATIVE_TARGET_BASELINE_TEST",
        "transfer_mean_var_scaled": "INFORMATIVE_NON_CONE_METHOD",
        "hybrid": "INFORMATIVE_NON_CONE_METHOD",
    }

    for ln, grp in LINES_12:
        cohort = build_r9_lolo_cohort(lines, ln)
        fc_bar, t_mean, grand_mean = compute_source_fc_bar(
            eval_targets=cohort.eval_targets,
            source_names=cohort.source_names,
            source_deltas_by_line=cohort.source_deltas_by_line,
            source_ctrl_mean_libnorm=cohort.source_ctrl_mean_libnorm,
        )
        preds, _, _, _ = build_predictions(
            cohort,
            include_sa3_anchors=True,
            include_exploratory_e=True,
        )
        per_p, xx_o, haa_o, hbb_o, hab_o, n5_meta = prepare_n5_per_pert_data(
            delta_obs=cohort.delta_obs,
            t_mean=t_mean,
            grand_mean=grand_mean,
            fc_bar=fc_bar,
            target_cols=cohort.target_cols,
            half_a_obs=cohort.half_a_obs,
            half_b_obs=cohort.half_b_obs,
            half_valid=cohort.half_valid,
            mode="all_nontarget",
        )
        m = len(per_p)
        rng = np.random.default_rng(seed)
        draws = rng.integers(0, m, size=(n_boot, m)) if m > 0 else None

        eval_names = [
            "fct",
            "transfer_mean_global_scale",
            "source_grand_mean",
            "e1_const_med",
            "e2_const_p10",
            "e2_const_p90",
            "e6_log_fc_bar",
            "e7_fc_minus_1",
            "e3_src_mean_ctrl",
            "e5_rank_matched_mu",
            "transfer_mean_var_scaled",
            "hybrid",
        ] + [f"fct_mismatch::{s}" for s in cohort.source_names]

        line_methods: Dict[str, Any] = {}
        boot_arrays: Dict[str, np.ndarray] = {}
        for mn in eval_names:
            xy, yy_o, yy_c = _pred_n5_moments(per_p, cohort.delta_obs, preds[mn], t_mean)
            res = _eval_n5_from_moments(xy, yy_o, yy_c, xx_o, haa_o, hbb_o, hab_o, draws, n_boot)
            b_arr = res.pop("_b_cal")
            boot_arrays[mn] = b_arr
            cls_tag = construction_classes.get(mn, "INFORMATIVE_MISMATCH_BASELINE")
            res["audit_classification"] = cls_tag
            line_methods[mn] = res
            method_rows.append({
                "line": ln,
                "group": grp,
                "n_eval": m,
                "method": mn,
                "audit_classification": cls_tag,
                "r_cone_sp_est": res["r_cone_sp_est"],
                "r_cone_sp_ci975_lo": res["r_cone_sp_ci975_lo"],
                "r_cone_sp_ci975_hi": res["r_cone_sp_ci975_hi"],
                "r_cone_sp_raw_est": res["r_cone_sp_raw_est"],
                "r_cone_orth_raw": res["r_cone_orth_raw"],
                "eta_est": res["eta_est"],
                "r2_in_n5": res["r2_in_n5"],
                "rho_ab_n5": res["rho_ab_n5"],
                "valid_frac": res["r_cone_sp_valid_frac"],
            })

        # Single-draw gene permutations Perm0..99 (non-expected)
        mu_safe = np.maximum(np.asarray(cohort.tgt_ctrl_mean_libnorm, dtype=np.float64), 0.0)
        n_g = cohort.delta_obs.shape[1]
        p_cal = np.empty(n_perm, dtype=np.float64)
        p_raw = np.empty(n_perm, dtype=np.float64)
        p_eta = np.empty(n_perm, dtype=np.float64)
        for k in range(n_perm):
            perm = np.random.default_rng(k).permutation(n_g)
            pred_pk = fct_from_fc_bar(fc_bar, mu_safe[perm])
            xy, yy_o, yy_c = _pred_n5_moments(per_p, cohort.delta_obs, pred_pk, t_mean)
            r_k = _eval_n5_from_moments(xy, yy_o, yy_c, xx_o, haa_o, hbb_o, hab_o, draws=None, n_boot=1)
            r_k.pop("_b_cal", None)
            p_cal[k] = r_k["r_cone_sp_est"]
            p_raw[k] = r_k["r_cone_sp_raw_est"]
            p_eta[k] = r_k["eta_est"]

        fin_c = p_cal[np.isfinite(p_cal)]
        perm_summary = {
            "n_perm": n_perm,
            "min": float(np.min(fin_c)) if fin_c.size else float("nan"),
            "p025": float(np.percentile(fin_c, 2.5)) if fin_c.size else float("nan"),
            "median": float(np.median(fin_c)) if fin_c.size else float("nan"),
            "mean": float(np.mean(fin_c)) if fin_c.size else float("nan"),
            "p975": float(np.percentile(fin_c, 97.5)) if fin_c.size else float("nan"),
            "max": float(np.max(fin_c)) if fin_c.size else float("nan"),
            "median_abs_r_cone_sp": float(np.median(np.abs(fin_c))) if fin_c.size else float("nan"),
            "median_eta": float(np.median(p_eta)),
            "min_eta": float(np.min(p_eta)),
            "max_eta": float(np.max(p_eta)),
            "single_draws": [float(v) for v in p_cal],
        }
        perm_rows.append({
            "line": ln,
            "group": grp,
            "n_eval": m,
            "perm_min": perm_summary["min"],
            "perm_p025": perm_summary["p025"],
            "perm_median": perm_summary["median"],
            "perm_mean": perm_summary["mean"],
            "perm_p975": perm_summary["p975"],
            "perm_max": perm_summary["max"],
            "perm_median_abs": perm_summary["median_abs_r_cone_sp"],
            "perm_median_eta": perm_summary["median_eta"],
        })

        # Pair FCT vs all 6 Mis_s (paired joint bootstrap diff)
        mis_details: Dict[str, Any] = {}
        for s in cohort.source_names:
            mk = f"fct_mismatch::{s}"
            pt_d = line_methods["fct"]["r_cone_sp_est"] - line_methods[mk]["r_cone_sp_est"]
            b_d = boot_arrays["fct"] - boot_arrays[mk]
            ds = summarize_boot(pt_d, b_d, m, n_boot).to_dict("diff")
            mis_details[s] = {
                "mis_est": line_methods[mk]["r_cone_sp_est"],
                "mis_ci975_lo": line_methods[mk]["r_cone_sp_ci975_lo"],
                "mis_ci975_hi": line_methods[mk]["r_cone_sp_ci975_hi"],
                "mis_eta": line_methods[mk]["eta_est"],
                "fct_minus_mis": ds,
                "unpaired_fct_lo_gt_mis_hi": bool(line_methods["fct"]["r_cone_sp_ci975_lo"] > line_methods[mk]["r_cone_sp_ci975_hi"]),
                "paired_ci975_gt_0": bool(ds["diff_ci975_lo"] > 0.0),
            }

        max_mis_name = max(cohort.source_names, key=lambda s: line_methods[f"fct_mismatch::{s}"]["r_cone_sp_ci975_hi"])
        vcc_cmp = mis_details.get("vcc2025_h1")
        if vcc_cmp is not None:
            spotlight_rows.append({
                "line": ln,
                "group": grp,
                "n_eval": m,
                "n5_fct_est": line_methods["fct"]["r_cone_sp_est"],
                "n5_fct_ci975_lo": line_methods["fct"]["r_cone_sp_ci975_lo"],
                "n5_fct_ci975_hi": line_methods["fct"]["r_cone_sp_ci975_hi"],
                "n5_vcc_mis_est": vcc_cmp["mis_est"],
                "n5_vcc_mis_ci975_lo": vcc_cmp["mis_ci975_lo"],
                "n5_vcc_mis_ci975_hi": vcc_cmp["mis_ci975_hi"],
                "unpaired_fct_lo_minus_vcc_hi": line_methods["fct"]["r_cone_sp_ci975_lo"] - vcc_cmp["mis_ci975_hi"],
                "paired_fct_minus_vcc_est": vcc_cmp["fct_minus_mis"]["diff_est"],
                "paired_fct_minus_vcc_ci975_lo": vcc_cmp["fct_minus_mis"]["diff_ci975_lo"],
                "paired_fct_minus_vcc_ci975_hi": vcc_cmp["fct_minus_mis"]["diff_ci975_hi"],
                "unpaired_pass": vcc_cmp["unpaired_fct_lo_gt_mis_hi"],
                "paired_pass": vcc_cmp["paired_ci975_gt_0"],
                "perm100_p975": perm_summary["p975"],
                "max_mis_source": max_mis_name,
                "max_mis_ci975_hi": mis_details[max_mis_name]["mis_ci975_hi"],
            })

        per_line[ln] = {
            "group": grp,
            "n5_meta": n5_meta,
            "methods": line_methods,
            "perm100_single_draws": perm_summary,
            "mismatched_vs_fct": mis_details,
            "max_mis_source": max_mis_name,
            "max_mis_ci975_hi": mis_details[max_mis_name]["mis_ci975_hi"],
        }

        print(
            f"[Check2 N5] {ln:24s} N={m:4d} rank={n5_meta['rank_min']}-{n5_meta['rank_max']} | "
            f"FCT={line_methods['fct']['r_cone_sp_est']:+.4f} [{line_methods['fct']['r_cone_sp_ci975_lo']:+.4f},{line_methods['fct']['r_cone_sp_ci975_hi']:+.4f}] "
            f"(eta={line_methods['fct']['eta_est']:.3f}) | "
            f"GS={line_methods['transfer_mean_global_scale']['r_cone_sp_est']:+.1e}(eta={line_methods['transfer_mean_global_scale']['eta_est']:.1e}) "
            f"SGM={line_methods['source_grand_mean']['r_cone_sp_est']:+.1e}(eta={line_methods['source_grand_mean']['eta_est']:.1e}) "
            f"E1={line_methods['e1_const_med']['r_cone_sp_est']:+.4f}(eta={line_methods['e1_const_med']['eta_est']:.4f}) | "
            f"Perm97.5={perm_summary['p975']:+.4f} MaxMis={max_mis_name}({mis_details[max_mis_name]['mis_ci975_hi']:+.4f})"
        )

    payload = {
        "check": "N_5_independent_audit",
        "seed": seed,
        "n_boot": n_boot,
        "n_perm": n_perm,
        "elapsed_sec": round(time.time() - t0, 2),
        "per_line": per_line,
        "spotlight_vs_vcc2025_h1": spotlight_rows,
    }
    with open(out_dir / "n5_audit.json", "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    with open(out_dir / "n5_methods_12lines.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(method_rows[0].keys()))
        w.writeheader()
        w.writerows(method_rows)
    with open(out_dir / "n5_perm100_12lines.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(perm_rows[0].keys()))
        w.writeheader()
        w.writerows(perm_rows)
    with open(out_dir / "n5_spotlight_vs_vcc2025_h1.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(spotlight_rows[0].keys()))
        w.writeheader()
        w.writerows(spotlight_rows)
    return payload


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=DEFAULT_AUDIT_SEED)
    parser.add_argument("--n-boot", type=int, default=4000)
    parser.add_argument("--n-perm", type=int, default=100)
    parser.add_argument("--out-dir", type=Path, default=EXPLORATORY_CHECK_DIR / "outputs")
    args = parser.parse_args()
    run_check2(seed=args.seed, n_boot=args.n_boot, n_perm=args.n_perm, out_dir=args.out_dir)
