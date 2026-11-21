#!/usr/bin/env python3
"""
Round 8 (R8-2 & R8-3 Exploratory):
  - R8-2: PC1 Axis vs. Off-Axis Mechanism Decomposition across all 3 datasets
    (`xatlas`, `Replogle 2022`, `Nadig 2025`, all 10 splits `SEEDS_10`).
    Decomposes `delta_r` into:
      * `axis_r`     : across-perturbation Pearson between true and predicted projections on training PC1 axis `a_c`.
      * `off_axis_r` : per-perturbation Pearson after removing training PC1 axis `a_c`.
    Adjudicates the pre-declared rule (ROUND8_PREREGISTRATION.md Section 2.2):
      mean(off_axis_r_AE - off_axis_r_S5) <= (1/3) * mean(delta_r_AE - delta_r_S5).

  - R8-3 Exploratory: Parameter-Free Orthogonal Decomposition Fusion (`OrthoFusion_scGPT_S5`)
    on already-seen datasets (`Replogle 2022`, `Nadig 2025`, all 10 splits `SEEDS_10`,
    strictly labeled EXPLORATORY):
      y_hat_fusion = b_c + P_S(y_hat_scGPT - b_c) + P_perp(y_hat_S5 - b_c),
    where S = span{1, b_c, a_c} (QR orthonormal basis fit strictly on training set, ZERO tunable weights).
    Evaluates H4 and H5 with 4,000-resample cell-line-stratified row-level paired bootstrap 95% CIs.
"""

import argparse
import gc
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.linear_model import Ridge

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from compute_row_bootstrap_ci import row_scores_v5, stratified_row_bootstrap_ci
from replogle_external_validation import (
    load_replogle_bundle,
    make_strict_union_split,
    run_ahlmann_eltze_bilinear,
)
from round6_t2_discrim_bootstrap_ci import compute_row_metrics_for_split
from round6_t2_multi_metric_bench import evaluate_5metrics, load_xatlas_bundle, make_xatlas_split
from xatlas_panhuman_solver import (
    REPORTS_DIR,
    SEEDS_10,
    PanPerturbBridgeModel,
    fit_baseline_and_axis,
)


def resolve_scgpt_dir() -> Path:
    return Path(os.environ.get("SCGPT_WEIGHTS_DIR", str(Path.home() / ".cache/r7_scgpt_weights")))


def pearson_1d(u: np.ndarray, v: np.ndarray) -> float:
    u0 = u - u.mean()
    v0 = v - v.mean()
    d = float(np.linalg.norm(u0) * np.linalg.norm(v0))
    return float(np.clip(np.dot(u0, v0) / d, -1.0, 1.0)) if d > 1e-12 else 0.0


def orthonormal_subspace_basis(b: np.ndarray, a: np.ndarray) -> np.ndarray:
    """Orthonormal QR basis Q for S = span{1, b, a} (strictly training-set derived)."""
    ones = np.ones_like(b, dtype=np.float64)
    M = np.column_stack([ones, b.astype(np.float64), a.astype(np.float64)])
    Q, R = np.linalg.qr(M, mode="reduced")
    keep = np.abs(np.diag(R)) > 1e-10
    return Q[:, keep]


def ortho_fusion_predict(p_scgpt: dict, p_s5: dict, val_df: pd.DataFrame, bases: dict, axes: dict) -> dict:
    """
    Exact parameter-free orthogonal decomposition fusion (ROUND8_PREREGISTRATION.md Section 3.1):
      y_hat = b + P_S(y_hat_scGPT - b) + P_perp(y_hat_S5 - b),  S = span{1, b, a}.
    """
    out = {}
    Q_map = {src: orthonormal_subspace_basis(bases[src], axes[src]) for src in ["HCT116", "HEK293T"]}
    for _, row in val_df.iterrows():
        rid = row["id"]
        src = row["source"]
        b = bases[src].astype(np.float64)
        Q = Q_map[src]
        d_sc = p_scgpt[rid].astype(np.float64) - b
        d_s5 = p_s5[rid].astype(np.float64) - b
        proj_sc = Q @ (Q.T @ d_sc)
        perp_s5 = d_s5 - Q @ (Q.T @ d_s5)
        out[rid] = b + proj_sc + perp_s5
    return out


def bootstrap_axis_r_diff(
    z_true_by_src: dict,
    z_a_by_src: dict,
    z_b_by_src: dict,
    n_boot: int = 4000,
    seed: int = 12345,
) -> dict:
    """4,000-resample cell-line-stratified row-level paired bootstrap for axis_r(A) - axis_r(B)."""
    rng = np.random.RandomState(seed)
    srcs = ["HCT116", "HEK293T"]
    r_a_obs = 0.5 * sum(pearson_1d(z_true_by_src[s], z_a_by_src[s]) for s in srcs)
    r_b_obs = 0.5 * sum(pearson_1d(z_true_by_src[s], z_b_by_src[s]) for s in srcs)
    delta_obs = r_a_obs - r_b_obs

    boots = np.empty(n_boot, dtype=np.float64)
    idx_cache = {s: rng.randint(0, len(z_true_by_src[s]), size=(n_boot, len(z_true_by_src[s]))) for s in srcs}
    for b_i in range(n_boot):
        da = 0.0
        for s in srcs:
            ii = idx_cache[s][b_i]
            zt = z_true_by_src[s][ii]
            da += 0.5 * (pearson_1d(zt, z_a_by_src[s][ii]) - pearson_1d(zt, z_b_by_src[s][ii]))
        boots[b_i] = da
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return {
        "score_a": float(r_a_obs),
        "score_b": float(r_b_obs),
        "delta": float(delta_obs),
        "ci95_lo": float(lo),
        "ci95_hi": float(hi),
        "boot_se": float(boots.std(ddof=1)),
        "ci_gt_zero": bool(lo > 0.0),
        "ci_lt_zero": bool(hi < 0.0),
    }


def compute_r82_decomposition(
    preds_by_method: dict,
    tr_df: pd.DataFrame,
    val_df: pd.DataFrame,
    all_y: dict,
    score_indices: dict,
    gene2idx: dict,
) -> tuple[dict, dict, dict]:
    """
    Compute per-perturbation `delta_r`, per-perturbation `off_axis_r`, and across-perturbation `axis_r`
    using strictly leak-free training `(b_c, a_c)` fit on `tr_df[reproducible == 1]`.
    """
    methods = list(preds_by_method.keys())
    n = len(val_df)
    row_delta = {m: np.empty(n, dtype=np.float64) for m in methods}
    row_off = {m: np.empty(n, dtype=np.float64) for m in methods}
    z_true_by_src = {}
    z_pred_by_method = {m: {} for m in methods}

    for src in ["HCT116", "HEK293T"]:
        pos = np.where(val_df["source"].values == src)[0]
        sub = val_df.iloc[pos].reset_index(drop=True)
        s_idx = score_indices[f"{src}_indices"]

        # Fit strictly leak-free training centroid b_c and training PC1 unit vector a_c
        tr_rep = tr_df[(tr_df["source"] == src) & (tr_df["reproducible"] == 1)]
        b_full, a_full = fit_baseline_and_axis(tr_rep, all_y, s_idx)
        # Unit-norm PC1 restricted to scored genes s_idx
        a_s = a_full[s_idx].astype(np.float64)
        a_norm = float(np.linalg.norm(a_s))
        a_s = a_s / a_norm if a_norm > 1e-12 else a_s
        b_s = b_full[s_idx].astype(np.float64)

        tids = [gene2idx[src].get(p) for p in sub["perturbation"]]
        Y_true_s = np.array([all_y[rid].astype(np.float64)[s_idx] for rid in sub["id"]], dtype=np.float64)
        D_true_s = Y_true_s - b_s[None, :]
        z_true = D_true_s @ a_s
        z_true_by_src[src] = z_true

        for m in methods:
            preds = preds_by_method[m]
            Y_pred_s = np.array([preds[rid].astype(np.float64)[s_idx] for rid in sub["id"]], dtype=np.float64)
            D_pred_s = Y_pred_s - b_s[None, :]
            z_pred = D_pred_s @ a_s
            z_pred_by_method[m][src] = z_pred

            D_true_perp = D_true_s - np.outer(z_true, a_s)
            D_pred_perp = D_pred_s - np.outer(z_pred, a_s)

            r_d = np.empty(len(sub), dtype=np.float64)
            r_o = np.empty(len(sub), dtype=np.float64)
            # Map s_idx positions for target-gene masking
            idx_to_pos = {int(g): k for k, g in enumerate(s_idx)}
            for i in range(len(sub)):
                t = tids[i]
                t_pos = idx_to_pos.get(t) if t is not None else None
                if t_pos is not None:
                    mask = np.ones(len(s_idx), dtype=bool)
                    mask[t_pos] = False
                    r_d[i] = pearson_1d(D_true_s[i, mask], D_pred_s[i, mask])
                    r_o[i] = pearson_1d(D_true_perp[i, mask], D_pred_perp[i, mask])
                else:
                    r_d[i] = pearson_1d(D_true_s[i], D_pred_s[i])
                    r_o[i] = pearson_1d(D_true_perp[i], D_pred_perp[i])
            row_delta[m][pos] = r_d
            row_off[m][pos] = r_o

    return row_delta, row_off, {"z_true": z_true_by_src, "z_pred": z_pred_by_method}


def run_dataset(dataset: str, emb_mat=None, vocab=None):
    t0 = time.time()
    print(f"\n================ Round 8 R8-2 (+R8-3 Exploratory) on {dataset} ================", flush=True)
    if dataset == "xatlas":
        train_df, all_y, bases_pub, axes_pub, sidx, sh0, sh1, g2i_map = load_xatlas_bundle()
    elif dataset == "replogle":
        train_df, all_y, _, sh0, sh1, g2i_map, sidx = load_replogle_bundle()
    elif dataset == "nadig":
        from prepare_nadig2025_cache import load_nadig_bundle
        train_df, all_y, _, sh0, sh1, g2i_map, sidx = load_nadig_bundle()
    else:
        raise ValueError(f"Unknown dataset: {dataset}")

    has_scgpt = (dataset in ("replogle", "nadig")) and (emb_mat is not None)
    pad_id = vocab["<pad>"] if has_scgpt else None
    model = PanPerturbBridgeModel(g2i_map, sh0, sh1, sidx)

    r82_rows = []
    r82_ci_rows = []
    r83_rows = []
    r83_ci_rows = []

    for seed in SEEDS_10:
        ts = time.time()
        if dataset == "xatlas":
            val_df, vr_set, vt_line = make_xatlas_split(train_df, g2i_map, seed)
            bases, axes = bases_pub, axes_pub
            df_hct = train_df[(train_df["source"] == "HCT116") & (~train_df["id"].isin(vr_set))]
            df_hek = train_df[(train_df["source"] == "HEK293T") & (~train_df["id"].isin(vr_set))]
            df_oth_for_hct = train_df[(train_df["source"] == "HEK293T") & (~train_df["perturbation"].isin(vt_line["HCT116"]))]
            df_oth_for_hek = train_df[(train_df["source"] == "HCT116") & (~train_df["perturbation"].isin(vt_line["HEK293T"]))]
            tr_df = pd.concat([df_hct, df_hek], ignore_index=True)
        else:
            val_df, tr_df, _, bases, axes = make_strict_union_split(train_df, g2i_map, all_y, sidx, seed)
            df_hct = tr_df[tr_df["source"] == "HCT116"]
            df_hek = tr_df[tr_df["source"] == "HEK293T"]
            df_oth_for_hct = df_hek
            df_oth_for_hek = df_hct

        q_hct = val_df[val_df["source"] == "HCT116"]
        q_hek = val_df[val_df["source"] == "HEK293T"]

        # AE + optional scGPT_FrozenEmb_Ridge
        p_ae, p_scgpt = {}, {}
        for src, df_s, q_s in [("HCT116", df_hct, q_hct), ("HEK293T", df_hek, q_hek)]:
            g2i = g2i_map[src]
            own_rep = df_s[df_s["reproducible"] == 1]
            tr_sub = own_rep[own_rep["perturbation"].isin(g2i.keys())]
            tr_idx = [g2i[p] for p in tr_sub["perturbation"].tolist()]
            q_idx = [g2i[p] for p in q_s["perturbation"].tolist()]
            tr_Y_raw = np.array([all_y[rid] for rid in tr_sub["id"]], dtype=np.float32)
            Y_all_raw = np.array([all_y[rid] for rid in df_s["id"]], dtype=np.float32)
            arr_ae = run_ahlmann_eltze_bilinear(Y_all_raw, tr_idx, q_idx, tr_Y_raw, rank=100, lam=10.0)
            for i, rid in enumerate(q_s["id"].tolist()):
                p_ae[rid] = arr_ae[i]

            if has_scgpt:
                tr_syms = own_rep["gene_symbol"].astype(str).tolist()
                q_syms = q_s["gene_symbol"].astype(str).tolist()
                X_tr = np.array([emb_mat[vocab.get(s, pad_id)] for s in tr_syms], dtype=np.float32)
                X_q = np.array([emb_mat[vocab.get(s, pad_id)] for s in q_syms], dtype=np.float32)
                Y_tr_all = np.array([all_y[rid] for rid in own_rep["id"]], dtype=np.float32)
                p_arr = Ridge(alpha=10.0).fit(X_tr, Y_tr_all).predict(X_q)
                for i, rid in enumerate(q_s["id"].tolist()):
                    p_scgpt[rid] = p_arr[i].astype(np.float64)

        r0 = model.predict_for_source_all_stages("HCT116", df_hct, df_oth_for_hct, q_hct, all_y, bases, axes, ("S0", "S5", "S7"))
        r1 = model.predict_for_source_all_stages("HEK293T", df_hek, df_oth_for_hek, q_hek, all_y, bases, axes, ("S0", "S5", "S7"))

        preds_map = {
            "AE": p_ae,
            "S0": {**r0["S0"], **r1["S0"]},
            "S5": {**r0["S5"], **r1["S5"]},
            "S7": {**r0["S7"], **r1["S7"]},
        }
        if has_scgpt:
            preds_map["scGPT_FrozenEmb_Ridge"] = p_scgpt
            preds_map["OrthoFusion_scGPT_S5"] = ortho_fusion_predict(p_scgpt, preds_map["S5"], val_df, bases, axes)

        # ---- R8-2: PC1 Axis vs Off-Axis Decomposition ----
        row_delta, row_off, z_pack = compute_r82_decomposition(preds_map, tr_df, val_df, all_y, sidx, g2i_map)
        src_arr = val_df["source"].values
        i0, i1 = np.where(src_arr == "HCT116")[0], np.where(src_arr == "HEK293T")[0]

        for m in preds_map:
            dr = float(0.5 * (row_delta[m][i0].mean() + row_delta[m][i1].mean()))
            offr = float(0.5 * (row_off[m][i0].mean() + row_off[m][i1].mean()))
            axr_l0 = pearson_1d(z_pack["z_true"]["HCT116"], z_pack["z_pred"][m]["HCT116"])
            axr_l1 = pearson_1d(z_pack["z_true"]["HEK293T"], z_pack["z_pred"][m]["HEK293T"])
            axr = float(0.5 * (axr_l0 + axr_l1))
            r82_rows.append({
                "dataset": dataset, "seed": seed, "method": m,
                "delta_r": round(dr, 6),
                "axis_r": round(axr, 6), "axis_r_L0": round(axr_l0, 6), "axis_r_L1": round(axr_l1, 6),
                "off_axis_r": round(offr, 6),
                "off_axis_r_L0": round(float(row_off[m][i0].mean()), 6),
                "off_axis_r_L1": round(float(row_off[m][i1].mean()), 6),
            })

        comps_r82 = [("AE_vs_S5", "AE", "S5"), ("S0_vs_S5", "S0", "S5")]
        if has_scgpt:
            comps_r82.append(("scGPT_Probe_vs_S5", "scGPT_FrozenEmb_Ridge", "S5"))
            comps_r82.append(("scGPT_Probe_vs_AE", "scGPT_FrozenEmb_Ridge", "AE"))

        for cname, ma, mb in comps_r82:
            for mname, rdict in [("delta_r", row_delta), ("off_axis_r", row_off)]:
                ci = stratified_row_bootstrap_ci(rdict[ma], rdict[mb], src_arr, n_boot=4000, seed=12345)
                r82_ci_rows.append({
                    "dataset": dataset, "seed": seed, "metric": mname, "comparison": cname,
                    **{k: round(v, 6) if isinstance(v, float) else v for k, v in ci.items()},
                })
            ci_ax = bootstrap_axis_r_diff(z_pack["z_true"], z_pack["z_pred"][ma], z_pack["z_pred"][mb], n_boot=4000, seed=12345)
            r82_ci_rows.append({
                "dataset": dataset, "seed": seed, "metric": "axis_r", "comparison": cname,
                **{k: round(v, 6) if isinstance(v, float) else v for k, v in ci_ax.items()},
            })

        # ---- R8-3 Exploratory: Parameter-Free Orthogonal Fusion (Replogle & Nadig) ----
        if has_scgpt:
            eval_models = ["OrthoFusion_scGPT_S5", "scGPT_FrozenEmb_Ridge", "S5", "AE"]
            sub_preds = {m: preds_map[m] for m in eval_models}
            rm = compute_row_metrics_for_split(sub_preds, val_df, all_y, bases, axes, sidx, g2i_map)
            v5_map = {m: row_scores_v5(sub_preds[m], val_df, all_y, bases, axes, sidx, g2i_map) for m in eval_models}
            # Verify exact algebraic identity v5(OrthoFusion) == v5(S5)
            max_v5_diff = float(np.max(np.abs(v5_map["OrthoFusion_scGPT_S5"] - v5_map["S5"])))
            assert max_v5_diff < 1e-10, f"Orthogonal fusion v5 parity violated: {max_v5_diff}"

            for m in eval_models:
                ev = evaluate_5metrics(sub_preds[m], val_df, all_y, bases, axes, sidx, g2i_map)
                r83_rows.append({
                    "dataset": dataset, "status": "EXPLORATORY_SEEN_DATASET", "seed": seed, "method": m,
                    **{k: round(v, 6) for k, v in ev.items()},
                    "disc_raw": round(float(0.5 * (rm[m]["disc_raw"][i0].mean() + rm[m]["disc_raw"][i1].mean())), 6),
                    "disc_v5": round(float(0.5 * (rm[m]["disc_v5"][i0].mean() + rm[m]["disc_v5"][i1].mean())), 6),
                    "v5_parity_vs_S5_max_abs": round(max_v5_diff, 12) if m == "OrthoFusion_scGPT_S5" else 0.0,
                })

            for cname, ma, mb in [
                ("OrthoFusion_vs_scGPT_Probe", "OrthoFusion_scGPT_S5", "scGPT_FrozenEmb_Ridge"),
                ("OrthoFusion_vs_S5", "OrthoFusion_scGPT_S5", "S5"),
                ("OrthoFusion_vs_AE", "OrthoFusion_scGPT_S5", "AE"),
            ]:
                for mname, a_vec, b_vec in [
                    ("v5", v5_map[ma], v5_map[mb]),
                    ("delta_r", rm[ma]["delta_r"], rm[mb]["delta_r"]),
                    ("disc_raw", rm[ma]["disc_raw"], rm[mb]["disc_raw"]),
                    ("disc_v5", rm[ma]["disc_v5"], rm[mb]["disc_v5"]),
                ]:
                    ci = stratified_row_bootstrap_ci(a_vec, b_vec, src_arr, n_boot=4000, seed=12345)
                    r83_ci_rows.append({
                        "dataset": dataset, "status": "EXPLORATORY_SEEN_DATASET",
                        "seed": seed, "metric": mname, "comparison": cname,
                        **{k: round(v, 6) if isinstance(v, float) else v for k, v in ci.items()},
                    })

        print(f"  [{dataset}] seed={seed} done in {time.time()-ts:.1f}s", flush=True)

    # Save R8-2 outputs & summary
    df_r82 = pd.DataFrame(r82_rows)
    df_r82_ci = pd.DataFrame(r82_ci_rows)
    df_r82.to_csv(REPORTS_DIR / f"r82_axis_decomposition_{dataset}_metrics.csv", index=False)
    df_r82_ci.to_csv(REPORTS_DIR / f"r82_axis_decomposition_{dataset}_row_ci.csv", index=False)

    r82_summary = {"dataset": dataset, "n_splits": len(SEEDS_10), "methods": {}, "comparisons": {}}
    for m, g in df_r82.groupby("method"):
        r82_summary["methods"][m] = {
            col: {"mean": round(float(g[col].mean()), 6), "std": round(float(g[col].std(ddof=1)), 6)}
            for col in ["delta_r", "axis_r", "off_axis_r"]
        }
    for (met, comp), g in df_r82_ci.groupby(["metric", "comparison"]):
        r82_summary["comparisons"].setdefault(met, {})[comp] = {
            "mean_delta": round(float(g["delta"].mean()), 6),
            "std_delta": round(float(g["delta"].std(ddof=1)), 6),
            "mean_ci95_lo": round(float(g["ci95_lo"].mean()), 6),
            "mean_ci95_hi": round(float(g["ci95_hi"].mean()), 6),
            "splits_ci_gt_zero": f"{int(g['ci_gt_zero'].sum())}/{len(g)}",
            "splits_ci_lt_zero": f"{int(g['ci_lt_zero'].sum())}/{len(g)}",
        }

    # Pre-declared R8-2 adjudication
    d_raw_ae = r82_summary["comparisons"]["delta_r"]["AE_vs_S5"]["mean_delta"]
    d_off_ae = r82_summary["comparisons"]["off_axis_r"]["AE_vs_S5"]["mean_delta"]
    ratio_ae = d_off_ae / d_raw_ae if abs(d_raw_ae) > 1e-12 else float("nan")
    passed_ae = bool(d_off_ae <= (1.0 / 3.0) * d_raw_ae)
    r82_summary["adjudication"] = {
        "criterion": "mean(off_axis_r_AE - off_axis_r_S5) <= (1/3) * mean(delta_r_AE - delta_r_S5)",
        "AE_minus_S5_delta_r": round(d_raw_ae, 6),
        "AE_minus_S5_off_axis_r": round(d_off_ae, 6),
        "off_to_raw_ratio": round(ratio_ae, 4),
        "passed": passed_ae,
        "verdict": "CONFIRMED" if passed_ae else "NOT_CONFIRMED",
    }
    if has_scgpt:
        d_raw_sc = r82_summary["comparisons"]["delta_r"]["scGPT_Probe_vs_S5"]["mean_delta"]
        d_off_sc = r82_summary["comparisons"]["off_axis_r"]["scGPT_Probe_vs_S5"]["mean_delta"]
        r82_summary["adjudication"]["scGPT_minus_S5_delta_r"] = round(d_raw_sc, 6)
        r82_summary["adjudication"]["scGPT_minus_S5_off_axis_r"] = round(d_off_sc, 6)
        r82_summary["adjudication"]["scGPT_off_to_raw_ratio"] = round(d_off_sc / d_raw_sc, 4) if abs(d_raw_sc) > 1e-12 else None

    with open(REPORTS_DIR / f"r82_axis_decomposition_{dataset}_summary.json", "w") as f:
        json.dump(r82_summary, f, indent=2)

    print(f"  [R8-2 {dataset}] AE-S5: delta_r={d_raw_ae:+.4f} -> off_axis_r={d_off_ae:+.4f} (ratio={ratio_ae:.3f}) => {r82_summary['adjudication']['verdict']}", flush=True)
    print(f"    axis_r: AE={r82_summary['methods']['AE']['axis_r']['mean']:.4f}, S5={r82_summary['methods']['S5']['axis_r']['mean']:.4f}", flush=True)
    if has_scgpt:
        print(f"    scGPT-S5: delta_r={d_raw_sc:+.4f} -> off_axis_r={d_off_sc:+.4f} (ratio={d_off_sc/d_raw_sc:.3f}), axis_r_scGPT={r82_summary['methods']['scGPT_FrozenEmb_Ridge']['axis_r']['mean']:.4f}", flush=True)

    # Save R8-3 Exploratory outputs & summary
    if has_scgpt:
        df_r83 = pd.DataFrame(r83_rows)
        df_r83_ci = pd.DataFrame(r83_ci_rows)
        df_r83.to_csv(REPORTS_DIR / f"r83_exploratory_orthofusion_{dataset}_metrics.csv", index=False)
        df_r83_ci.to_csv(REPORTS_DIR / f"r83_exploratory_orthofusion_{dataset}_row_ci.csv", index=False)

        r83_summary = {
            "dataset": dataset,
            "status": "EXPLORATORY_SEEN_DATASET",
            "formula": "y_hat = b + P_S(y_hat_scGPT - b) + P_perp(y_hat_S5 - b), S = span{1, b, a}, ZERO tunable weights",
            "n_splits": len(SEEDS_10),
            "methods": {},
            "comparisons": {},
        }
        for m, g in df_r83.groupby("method"):
            r83_summary["methods"][m] = {
                col: {"mean": round(float(g[col].mean()), 6), "std": round(float(g[col].std(ddof=1)), 6)}
                for col in ["v5", "delta_r", "de50_r", "de50_sign", "disc_raw", "disc_v5"]
            }
        for (met, comp), g in df_r83_ci.groupby(["metric", "comparison"]):
            r83_summary["comparisons"].setdefault(met, {})[comp] = {
                "mean_delta": round(float(g["delta"].mean()), 6),
                "std_delta": round(float(g["delta"].std(ddof=1)), 6),
                "mean_ci95_lo": round(float(g["ci95_lo"].mean()), 6),
                "mean_ci95_hi": round(float(g["ci95_hi"].mean()), 6),
                "splits_ci_gt_zero": f"{int(g['ci_gt_zero'].sum())}/{len(g)}",
                "splits_ci_lt_zero": f"{int(g['ci_lt_zero'].sum())}/{len(g)}",
            }

        c_dr = r83_summary["comparisons"]["delta_r"]["OrthoFusion_vs_scGPT_Probe"]
        c_v5 = r83_summary["comparisons"]["v5"]["OrthoFusion_vs_S5"]
        c_d1 = r83_summary["comparisons"]["disc_raw"]["OrthoFusion_vs_scGPT_Probe"]
        c_d2 = r83_summary["comparisons"]["disc_raw"]["OrthoFusion_vs_S5"]

        h4_dr_neg = int(c_dr["splits_ci_lt_zero"].split("/")[0])
        h4_v5_neg = int(c_v5["splits_ci_lt_zero"].split("/")[0])
        h5_sc_pos = int(c_d1["splits_ci_gt_zero"].split("/")[0])
        h5_s5_pos = int(c_d2["splits_ci_gt_zero"].split("/")[0])

        h4_pass = (h4_dr_neg <= 2) and (h4_v5_neg <= 2)
        h5_pass = (h5_sc_pos >= 7) and (c_d1["mean_delta"] > 0) and (h5_s5_pos >= 7) and (c_d2["mean_delta"] > 0)

        r83_summary["exploratory_adjudication"] = {
            "note": "EXPLORATORY ONLY on already-seen dataset; confirmatory test requires unseen 4th dataset.",
            "H4_non_inferiority": {
                "delta_r_vs_scGPT": c_dr,
                "v5_vs_S5": c_v5,
                "passed": h4_pass,
            },
            "H5_disc_raw_synergy": {
                "disc_raw_vs_scGPT": c_d1,
                "disc_raw_vs_S5": c_d2,
                "passed": h5_pass,
            },
        }
        with open(REPORTS_DIR / f"r83_exploratory_orthofusion_{dataset}_summary.json", "w") as f:
            json.dump(r83_summary, f, indent=2)

        of = r83_summary["methods"]["OrthoFusion_scGPT_S5"]
        print(
            f"  [R8-3 EXPLORATORY {dataset}] OrthoFusion: delta_r={of['delta_r']['mean']:.4f} (vs scGPT {c_dr['mean_delta']:+.4f}, neg={c_dr['splits_ci_lt_zero']}), "
            f"v5={of['v5']['mean']:.4f} (vs S5 {c_v5['mean_delta']:+.4f}), "
            f"disc_raw={of['disc_raw']['mean']:.4f} (vs scGPT {c_d1['mean_delta']:+.4f} pos={c_d1['splits_ci_gt_zero']}; vs S5 {c_d2['mean_delta']:+.4f} pos={c_d2['splits_ci_gt_zero']}) "
            f"=> H4={h4_pass}, H5={h5_pass}",
            flush=True,
        )

    print(f"Done {dataset} in {time.time()-t0:.1f}s", flush=True)
    gc.collect()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="+", default=["xatlas", "replogle", "nadig"])
    args = ap.parse_args()

    emb_mat, vocab = None, None
    if any(d in ("replogle", "nadig") for d in args.datasets):
        wdir = resolve_scgpt_dir()
        vocab = json.load(open(wdir / "scgpt_vocab.json"))
        sd = torch.load(wdir / "scgpt_gh_repo_original_model.bin", map_location="cpu", weights_only=False)
        emb_mat = sd["encoder.embedding.weight"].float().numpy()
        del sd

    for ds in args.datasets:
        run_dataset(ds, emb_mat=emb_mat, vocab=vocab)


if __name__ == "__main__":
    main()
