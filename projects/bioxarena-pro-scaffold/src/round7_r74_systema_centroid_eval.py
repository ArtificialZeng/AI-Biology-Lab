#!/usr/bin/env python3
"""
Round 7 (R7-4): Systema-Style Perturbation-Centroid-Centered Evaluation Across All 3 Datasets
(xatlas-orion-perturb, Replogle 2022, Nadig 2025).

Metric Definition (Systema / Viñas Torné et al., Nature Biotechnology 2025):
  Instead of correlating raw delta expression (y_i - baseline), subtract the
  mean perturbation response across perturbations (the perturbation centroid)
  from BOTH true and predicted profiles before computing Pearson correlation:
    - Primary (strictly leak-free): `centroid_train_r`
        centroid = mean_{j in train_reproducible(cell_line)} y_j
        r_i = Pearson(y_i - centroid, y_hat_i - centroid) on s_idx \ {t_idx}
    - Companion (evaluated-set centroid): `centroid_val_r`
        centroid_val = mean_{j in val(cell_line)} y_j
        r_i = Pearson(y_i - centroid_val, y_hat_i - centroid_val) on s_idx \ {t_idx}
  Also computes 4,000-resample cell-line-stratified row-level paired bootstrap 95% CIs
  for S5 vs AE, S7 vs AE, S5 vs S0, S7 vs S0, and (exploratory) S7_no_S6 vs AE.
"""

import argparse
import gc
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from compute_row_bootstrap_ci import stratified_row_bootstrap_ci
from replogle_external_validation import (
    load_replogle_bundle,
    make_strict_union_split,
    run_ahlmann_eltze_bilinear,
    run_knn30,
)
from round6_t2_multi_metric_bench import load_xatlas_bundle, make_xatlas_split
from round6_t4_leave_one_out_ablation import FULL_CFG, run_variant_for_source
from xatlas_panhuman_solver import (
    REPORTS_DIR,
    SEEDS_10,
    PanPerturbBridgeModel,
)


def compute_centroid_rows(preds_by_method: dict, tr_df: pd.DataFrame, val_df: pd.DataFrame,
                          all_y: dict, score_indices: dict, gene2idx: dict) -> dict:
    """Return {method: {'centroid_train_r': arr, 'centroid_val_r': arr}} aligned with val_df."""
    methods = list(preds_by_method.keys())
    n = len(val_df)
    out = {
        m: {"centroid_train_r": np.empty(n, dtype=np.float64), "centroid_val_r": np.empty(n, dtype=np.float64)}
        for m in methods
    }

    for src in ["HCT116", "HEK293T"]:
        mask = (val_df["source"].values == src)
        pos_idx = np.where(mask)[0]
        sub = val_df.iloc[pos_idx].reset_index(drop=True)
        s_idx = score_indices[f"{src}_indices"]

        # Training perturbation centroid (strictly leak-free: mean over training reproducible perturbations)
        tr_rep = tr_df[(tr_df["source"] == src) & (tr_df["reproducible"] == 1)]
        tr_mat = np.array([all_y[rid] for rid in tr_rep["id"]], dtype=np.float64)
        c_train = tr_mat.mean(axis=0)

        # Held-out perturbation centroid (companion check matching in-set Systema centering)
        val_mat = np.array([all_y[rid] for rid in sub["id"]], dtype=np.float64)
        c_val = val_mat.mean(axis=0)

        tids = [gene2idx[src].get(p) for p in sub["perturbation"]]
        for m in methods:
            preds = preds_by_method[m]
            r_tr = np.empty(len(sub), dtype=np.float64)
            r_va = np.empty(len(sub), dtype=np.float64)
            for i, rid in enumerate(sub["id"]):
                t = tids[i]
                idx = s_idx[s_idx != t] if t is not None else s_idx
                yt = all_y[rid].astype(np.float64)[idx]
                yp = preds[rid].astype(np.float64)[idx]

                # Train-centroid centered Pearson
                a = yt - c_train[idx]
                b = yp - c_train[idx]
                a = a - a.mean()
                b = b - b.mean()
                da = np.linalg.norm(a) * np.linalg.norm(b)
                r_tr[i] = float(np.clip(np.dot(a, b) / da, -1.0, 1.0)) if da > 1e-12 else 0.0

                # Val-centroid centered Pearson
                a2 = yt - c_val[idx]
                b2 = yp - c_val[idx]
                a2 = a2 - a2.mean()
                b2 = b2 - b2.mean()
                da2 = np.linalg.norm(a2) * np.linalg.norm(b2)
                r_va[i] = float(np.clip(np.dot(a2, b2) / da2, -1.0, 1.0)) if da2 > 1e-12 else 0.0

            out[m]["centroid_train_r"][pos_idx] = r_tr
            out[m]["centroid_val_r"][pos_idx] = r_va
    return out


def run_dataset(dataset: str):
    t0 = time.time()
    print(f"=== R7-4 Systema-Style Centroid-Centered Eval on {dataset} (10 splits) ===", flush=True)
    if dataset == "xatlas":
        train_df, all_y, bases_pub, axes_pub, sidx, sh0, sh1, g2i_map = load_xatlas_bundle()
    elif dataset == "replogle":
        train_df, all_y, _, sh0, sh1, g2i_map, sidx = load_replogle_bundle()
    elif dataset == "nadig":
        from prepare_nadig2025_cache import load_nadig_bundle
        train_df, all_y, _, sh0, sh1, g2i_map, sidx = load_nadig_bundle()
    else:
        raise ValueError(f"Unknown dataset: {dataset}")

    model = PanPerturbBridgeModel(g2i_map, sh0, sh1, sidx)
    methods_order = ["Constant_Mean", "kNN_30", "AE", "S0", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S7_no_S6_exploratory"]

    score_rows = []
    ci_rows = []

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

        # Constant_Mean, kNN_30, AE
        p_const, p_knn, p_ae = {}, {}, {}
        for src, df_s, q_s in [("HCT116", df_hct, q_hct), ("HEK293T", df_hek, q_hek)]:
            g2i = g2i_map[src]
            own_rep = df_s[df_s["reproducible"] == 1]
            mu = np.mean([all_y[rid] for rid in own_rep["id"].tolist()], axis=0).astype(np.float64)
            for rid in q_s["id"]:
                p_const[rid] = mu

            tr_sub = own_rep[own_rep["perturbation"].isin(g2i.keys())]
            tr_idx = [g2i[p] for p in tr_sub["perturbation"].tolist()]
            q_idx = [g2i[p] for p in q_s["perturbation"].tolist()]
            tr_Y_raw = np.array([all_y[rid] for rid in tr_sub["id"]], dtype=np.float32)
            Y_all_raw = np.array([all_y[rid] for rid in df_s["id"]], dtype=np.float32)
            arr_knn = run_knn30(Y_all_raw, tr_idx, q_idx, tr_Y_raw, k=30)
            arr_ae = run_ahlmann_eltze_bilinear(Y_all_raw, tr_idx, q_idx, tr_Y_raw, rank=100, lam=10.0)
            for i, rid in enumerate(q_s["id"].tolist()):
                p_knn[rid] = arr_knn[i]
                p_ae[rid] = arr_ae[i]

        p0 = model.predict_for_source_all_stages("HCT116", df_hct, df_oth_for_hct, q_hct, all_y, bases, axes)
        p1 = model.predict_for_source_all_stages("HEK293T", df_hek, df_oth_for_hek, q_hek, all_y, bases, axes)

        # Exploratory S7_no_S6
        cfg_no_s6 = {**FULL_CFG, "use_reciprocal": False}
        l0 = run_variant_for_source("HCT116", df_hct, df_oth_for_hct, q_hct, all_y, bases, axes, g2i_map, sh0, sh1, sidx, cfg_no_s6)
        l1 = run_variant_for_source("HEK293T", df_hek, df_oth_for_hek, q_hek, all_y, bases, axes, g2i_map, sh0, sh1, sidx, cfg_no_s6)

        preds_by_method = {
            "Constant_Mean": p_const,
            "kNN_30": p_knn,
            "AE": p_ae,
            **{st: {**p0[st], **p1[st]} for st in ["S0", "S1", "S2", "S3", "S4", "S5", "S6", "S7"]},
            "S7_no_S6_exploratory": {**l0, **l1},
        }

        row_m = compute_centroid_rows(preds_by_method, tr_df, val_df, all_y, sidx, g2i_map)
        sources = val_df["source"].values
        i0 = np.where(sources == "HCT116")[0]
        i1 = np.where(sources == "HEK293T")[0]

        for m in methods_order:
            for met in ["centroid_train_r", "centroid_val_r"]:
                arr = row_m[m][met]
                m0 = float(arr[i0].mean())
                m1 = float(arr[i1].mean())
                score_rows.append({
                    "dataset": dataset,
                    "seed": seed,
                    "method": m,
                    "metric": met,
                    "macro_mean": round(0.5 * (m0 + m1), 6),
                    "line0_mean": round(m0, 6),
                    "line1_mean": round(m1, 6),
                })

        for met in ["centroid_train_r", "centroid_val_r"]:
            for comp, ma, mb in [
                ("S5_vs_AE", "S5", "AE"),
                ("S7_vs_AE", "S7", "AE"),
                ("S5_vs_S0", "S5", "S0"),
                ("S7_vs_S0", "S7", "S0"),
                ("S7_no_S6_vs_AE", "S7_no_S6_exploratory", "AE"),
                ("S6_vs_S5", "S6", "S5"),
            ]:
                res = stratified_row_bootstrap_ci(row_m[ma][met], row_m[mb][met], sources)
                ci_rows.append({
                    "dataset": dataset,
                    "seed": seed,
                    "metric": met,
                    "comparison": comp,
                    **res,
                })

        gc.collect()
        print(
            f"  [seed={seed}] centroid_train_r: AE={0.5*(row_m['AE']['centroid_train_r'][i0].mean()+row_m['AE']['centroid_train_r'][i1].mean()):.4f}, "
            f"S0={0.5*(row_m['S0']['centroid_train_r'][i0].mean()+row_m['S0']['centroid_train_r'][i1].mean()):.4f}, "
            f"S5={0.5*(row_m['S5']['centroid_train_r'][i0].mean()+row_m['S5']['centroid_train_r'][i1].mean()):.4f}, "
            f"S7={0.5*(row_m['S7']['centroid_train_r'][i0].mean()+row_m['S7']['centroid_train_r'][i1].mean()):.4f} ({time.time()-ts:.1f}s)",
            flush=True,
        )

    df_scores = pd.DataFrame(score_rows)
    df_ci = pd.DataFrame(ci_rows)
    df_scores.to_csv(REPORTS_DIR / f"r74_systema_{dataset}_scores.csv", index=False)
    df_ci.to_csv(REPORTS_DIR / f"r74_systema_{dataset}_row_ci.csv", index=False)

    summary = {"dataset": dataset, "n_splits": len(SEEDS_10), "methods": {}, "comparisons": {}}
    for met in ["centroid_train_r", "centroid_val_r"]:
        summary["methods"][met] = {}
        for m in methods_order:
            sub = df_scores[(df_scores["method"] == m) & (df_scores["metric"] == met)]
            summary["methods"][met][m] = {
                "mean": round(float(sub["macro_mean"].mean()), 6),
                "std": round(float(sub["macro_mean"].std(ddof=1)), 6),
            }
        summary["comparisons"][met] = {}
        for comp in ["S5_vs_AE", "S7_vs_AE", "S5_vs_S0", "S7_vs_S0", "S7_no_S6_vs_AE", "S6_vs_S5"]:
            sub = df_ci[(df_ci["comparison"] == comp) & (df_ci["metric"] == met)]
            summary["comparisons"][met][comp] = {
                "mean_delta": round(float(sub["delta"].mean()), 6),
                "std_delta": round(float(sub["delta"].std(ddof=1)), 6),
                "mean_ci95_lo": round(float(sub["ci95_lo"].mean()), 6),
                "mean_ci95_hi": round(float(sub["ci95_hi"].mean()), 6),
                "splits_ci_gt_zero": f"{int(sub['ci_gt_zero'].sum())}/{len(sub)}",
                "splits_ci_lt_zero": f"{int((sub['ci95_hi'] < 0).sum())}/{len(sub)}",
            }

    with open(REPORTS_DIR / f"r74_systema_{dataset}_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n=== R7-4 Summary ({dataset}, {time.time()-t0:.1f}s) ===", flush=True)
    for met in ["centroid_train_r", "centroid_val_r"]:
        print(f"  [{met}]", flush=True)
        for m in ["Constant_Mean", "kNN_30", "AE", "S0", "S5", "S6", "S7", "S7_no_S6_exploratory"]:
            v = summary["methods"][met][m]
            print(f"    {m:22s}: {v['mean']:.4f} +- {v['std']:.4f}", flush=True)
        for comp in ["S5_vs_AE", "S7_vs_AE", "S5_vs_S0", "S7_vs_S0", "S7_no_S6_vs_AE"]:
            c = summary["comparisons"][met][comp]
            print(f"    {comp:16s}: delta={c['mean_delta']:+.4f} (CI>0: {c['splits_ci_gt_zero']}, CI<0: {c['splits_ci_lt_zero']})", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", choices=["xatlas", "replogle", "nadig"], required=True)
    args = ap.parse_args()
    run_dataset(args.dataset)
