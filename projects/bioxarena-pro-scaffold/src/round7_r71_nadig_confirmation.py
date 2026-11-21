#!/usr/bin/env python3
"""
Round 7 (R7-1): Pre-Registered Unseen External Confirmation on Nadig et al. (Nature Genetics 2025, GSE264667)
HepG2 (slot HCT116) + Jurkat (slot HEK293T) Essential Perturb-seq.

Strictly adheres to ROUND7_PREREGISTRATION.md (commits c5294c0, 87f4bc2, b9fe2bb):
  - Frozen solver (`src/xatlas_panhuman_solver.py` @ 3609811, zero constant changes)
  - Primary method: S5; Secondary: S7; Exploratory only: S7_no_S6
  - Locked SEEDS_10 = [42, 0, 1, 2, 3, 7, 11, 99, 2024, 31337]
  - Strict simultaneous union holdout (`make_strict_union_split`)
  - 5 metrics + disc_raw + disc_v5 + 4,000-resample row-level stratified paired bootstrap 95% CIs
  - Automatic adjudication of pre-registered H1, H2, Main Thesis (H1 AND H2), and H3.

Modes (process-isolated to avoid macOS libomp/Accelerate conflicts):
  --mode dl  : Simplified_ResMLP & Simplified_GNN (torch-first, CPU), saves per-row metric arrays
  --mode pan : Constant_Mean, kNN_30, AE, S0..S7, S7_no_S6_exploratory + merges DL + row CIs + H1-H3 summary
"""

import argparse
import gc
import json
import sys
import time
from pathlib import Path

if "--mode" in sys.argv and "dl" in sys.argv:
    import torch
    torch.set_num_threads(4)

import numpy as np
import pandas as pd
from sklearn.decomposition import TruncatedSVD

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from compute_row_bootstrap_ci import row_scores_v5, stratified_row_bootstrap_ci
from prepare_nadig2025_cache import load_nadig_bundle
from replogle_external_validation import (
    make_strict_union_split,
    run_ahlmann_eltze_bilinear,
    run_knn30,
)
from round6_t2_discrim_bootstrap_ci import compute_row_metrics_for_split
from round6_t2_multi_metric_bench import evaluate_5metrics
from round6_t4_leave_one_out_ablation import FULL_CFG, run_variant_for_source
from xatlas_panhuman_solver import (
    REPORTS_DIR,
    SEEDS_10,
    PanPerturbBridgeModel,
    build_spectral_pyramid_Gs,
    project_rows_5000,
)

TMP_DL_DIR = REPORTS_DIR / "_tmp_nadig_dl_rows"


def run_dl():
    t0 = time.time()
    TMP_DL_DIR.mkdir(parents=True, exist_ok=True)
    import xatlas_sota_dl_benchmark as _dl_mod
    _dl_mod.DEVICE = torch.device("cpu")
    from xatlas_sota_dl_benchmark import train_predict_gnn, train_predict_mlp

    print("=== R7-1 (Nadig 2025 DL): Simplified ResMLP & GNN on 10 Strict Union Splits ===", flush=True)
    train_df, all_y, _, _, _, g2i_map, sidx = load_nadig_bundle()
    rows = []

    for seed in SEEDS_10:
        ts = time.time()
        val_df, df_tr, _, bases, axes = make_strict_union_split(train_df, g2i_map, all_y, sidx, seed)
        preds = {"Simplified_ResMLP": {}, "Simplified_GNN": {}}
        for src in ["HCT116", "HEK293T"]:
            g2i = g2i_map[src]
            s_i = sidx[f"{src}_indices"]
            df_own = df_tr[df_tr["source"] == src]
            q_df = val_df[val_df["source"] == src]
            q_idx = [g2i[p] for p in q_df["perturbation"].tolist()]
            own_rep = df_own[df_own["reproducible"] == 1]
            tr_sub = own_rep[own_rep["perturbation"].isin(g2i.keys())]
            tr_idx = [g2i[p] for p in tr_sub["perturbation"].tolist()]
            tr_Y_raw = np.array([all_y[rid] for rid in tr_sub["id"]], dtype=np.float32)
            tr_Y_p = project_rows_5000(tr_Y_raw, bases[src], axes[src])
            Y_all_raw = np.array([all_y[rid] for rid in df_own["id"]], dtype=np.float32)
            Y_clean = Y_all_raw.copy()
            for i, p in enumerate(df_own["perturbation"].tolist()):
                if p in g2i:
                    Y_clean[i, g2i[p]] = 0.0
            Y_clean = Y_clean - bases[src].astype(np.float32)[None, :]
            X_200 = TruncatedSVD(n_components=200, random_state=42).fit_transform(Y_clean.T).astype(np.float32)
            X_200 = X_200 / (np.linalg.norm(X_200, axis=1, keepdims=True) + 1e-8)

            p_mlp = train_predict_mlp(X_200, tr_idx, q_idx, tr_Y_p, bases[src], axes[src], s_i, seed=seed)
            Gs_own = build_spectral_pyramid_Gs(df_own, all_y, g2i, bases[src], pc_levels=(0, 1))
            G_single = 0.5 * (Gs_own[0] + Gs_own[1])
            p_gnn = train_predict_gnn(X_200, G_single, tr_idx, q_idx, tr_Y_p, bases[src], axes[src], s_i, seed=seed)
            for i, rid in enumerate(q_df["id"].tolist()):
                preds["Simplified_ResMLP"][rid] = p_mlp[i]
                preds["Simplified_GNN"][rid] = p_gnn[i]

        rm = compute_row_metrics_for_split(preds, val_df, all_y, bases, axes, sidx, g2i_map)
        src_arr = val_df["source"].values
        i0, i1 = np.where(src_arr == "HCT116")[0], np.where(src_arr == "HEK293T")[0]
        save_dict = {}
        for m in ["Simplified_ResMLP", "Simplified_GNN"]:
            ev = evaluate_5metrics(preds[m], val_df, all_y, bases, axes, sidx, g2i_map)
            v5_arr = row_scores_v5(preds[m], val_df, all_y, bases, axes, sidx, g2i_map)
            dr_raw = rm[m]["disc_raw"]
            dr_v5 = rm[m]["disc_v5"]
            save_dict[f"{m}__v5"] = v5_arr
            save_dict[f"{m}__delta_r"] = rm[m]["delta_r"]
            save_dict[f"{m}__disc_raw"] = dr_raw
            save_dict[f"{m}__disc_v5"] = dr_v5
            rows.append({
                "dataset": "nadig2025",
                "seed": seed,
                "method": m,
                **{k: round(v, 6) for k, v in ev.items()},
                "disc_raw": round(float(0.5 * (dr_raw[i0].mean() + dr_raw[i1].mean())), 6),
                "disc_v5": round(float(0.5 * (dr_v5[i0].mean() + dr_v5[i1].mean())), 6),
            })
        np.savez(TMP_DL_DIR / f"seed_{seed}.npz", **save_dict)
        g_r = [r for r in rows if r["seed"] == seed and r["method"] == "Simplified_GNN"][0]
        m_r = [r for r in rows if r["seed"] == seed and r["method"] == "Simplified_ResMLP"][0]
        print(
            f"  [nadig_dl] Seed {seed:5d} | GNN v5={g_r['v5']:.4f} disc_raw={g_r['disc_raw']:.4f} | "
            f"MLP v5={m_r['v5']:.4f} disc_raw={m_r['disc_raw']:.4f} [{time.time()-ts:.1f}s]",
            flush=True,
        )
        gc.collect()

    out_csv = REPORTS_DIR / "r71_nadig_5metric_dl.csv"
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    print(f"Saved {out_csv} in {time.time()-t0:.1f}s", flush=True)


def run_pan():
    t0 = time.time()
    print("=== R7-1 (Nadig 2025 PAN + Pre-Registered Confirmation H1-H3): 10 Strict Union Splits ===", flush=True)
    train_df, all_y, shared, sh0, sh1, g2i_map, sidx = load_nadig_bundle()
    model = PanPerturbBridgeModel(g2i_map, sh0, sh1, sidx)
    stages = ["S0", "S1", "S2", "S3", "S4", "S5", "S6", "S7"]
    methods = ["Constant_Mean", "kNN_30", "AE"] + stages + ["S7_no_S6_exploratory"]
    cfg_no_s6 = {**FULL_CFG, "use_reciprocal": False}

    metric_rows = []
    ci_rows = []
    split_stats = []

    for seed in SEEDS_10:
        ts = time.time()
        val_df, df_tr, union_val, bases, axes = make_strict_union_split(train_df, g2i_map, all_y, sidx, seed)
        df_hct = df_tr[df_tr["source"] == "HCT116"]
        df_hek = df_tr[df_tr["source"] == "HEK293T"]
        q_hct = val_df[val_df["source"] == "HCT116"]
        q_hek = val_df[val_df["source"] == "HEK293T"]
        split_stats.append({
            "seed": seed,
            "n_val_hepg2": int(len(q_hct)),
            "n_val_jurkat": int(len(q_hek)),
            "n_union_val_targets": int(len(union_val)),
            "n_tr_hepg2_rep": int((df_hct["reproducible"] == 1).sum()),
            "n_tr_jurkat_rep": int((df_hek["reproducible"] == 1).sum()),
        })

        preds = {m: {} for m in methods}
        for src, df_s, df_o, q_s in [
            ("HCT116", df_hct, df_hek, q_hct),
            ("HEK293T", df_hek, df_hct, q_hek),
        ]:
            g2i = g2i_map[src]
            res = model.predict_for_source_all_stages(src, df_s, df_o, q_s, all_y, bases, axes, tuple(stages))
            for st in stages:
                preds[st].update(res[st])

            # Exploratory S7_no_S6
            preds["S7_no_S6_exploratory"].update(
                run_variant_for_source(src, df_s, df_o, q_s, all_y, bases, axes, g2i_map, sh0, sh1, sidx, cfg_no_s6)
            )

            own_rep = df_s[df_s["reproducible"] == 1]
            mu = np.mean([all_y[rid] for rid in own_rep["id"].tolist()], axis=0).astype(np.float64)
            tr_sub = own_rep[own_rep["perturbation"].isin(g2i.keys())]
            tr_idx = [g2i[p] for p in tr_sub["perturbation"].tolist()]
            q_idx = [g2i[p] for p in q_s["perturbation"].tolist()]
            tr_Y_raw = np.array([all_y[rid] for rid in tr_sub["id"]], dtype=np.float32)
            Y_all_raw = np.array([all_y[rid] for rid in df_s["id"]], dtype=np.float32)
            p_knn = run_knn30(Y_all_raw, tr_idx, q_idx, tr_Y_raw, k=30)
            p_ae = run_ahlmann_eltze_bilinear(Y_all_raw, tr_idx, q_idx, tr_Y_raw, rank=100, lam=10.0)
            for i, rid in enumerate(q_s["id"].tolist()):
                preds["Constant_Mean"][rid] = mu
                preds["kNN_30"][rid] = p_knn[i]
                preds["AE"][rid] = p_ae[i]

        rm = compute_row_metrics_for_split(preds, val_df, all_y, bases, axes, sidx, g2i_map)
        src_arr = val_df["source"].values
        i0, i1 = np.where(src_arr == "HCT116")[0], np.where(src_arr == "HEK293T")[0]

        v5_rows = {}
        for m in methods:
            ev = evaluate_5metrics(preds[m], val_df, all_y, bases, axes, sidx, g2i_map)
            v5_rows[m] = row_scores_v5(preds[m], val_df, all_y, bases, axes, sidx, g2i_map)
            dr_raw = rm[m]["disc_raw"]
            dr_v5 = rm[m]["disc_v5"]
            metric_rows.append({
                "dataset": "nadig2025",
                "seed": seed,
                "method": m,
                **{k: round(v, 6) for k, v in ev.items()},
                "disc_raw": round(float(0.5 * (dr_raw[i0].mean() + dr_raw[i1].mean())), 6),
                "disc_v5": round(float(0.5 * (dr_v5[i0].mean() + dr_v5[i1].mean())), 6),
            })

        # Load DL row scores if available
        dl_npz = TMP_DL_DIR / f"seed_{seed}.npz"
        has_dl = dl_npz.exists()
        if has_dl:
            dl_d = np.load(dl_npz)
            for dm in ["Simplified_ResMLP", "Simplified_GNN"]:
                v5_rows[dm] = dl_d[f"{dm}__v5"]
                rm[dm] = {
                    "delta_r": dl_d[f"{dm}__delta_r"],
                    "disc_raw": dl_d[f"{dm}__disc_raw"],
                    "disc_v5": dl_d[f"{dm}__disc_v5"],
                }

        comps = [
            ("S5_vs_AE", "S5", "AE"),
            ("S7_vs_AE", "S7", "AE"),
            ("S7_no_S6_vs_AE", "S7_no_S6_exploratory", "AE"),
            ("S5_vs_S0", "S5", "S0"),
            ("S7_vs_S0", "S7", "S0"),
            ("S6_vs_S5", "S6", "S5"),
        ]
        if has_dl:
            comps += [
                ("S5_vs_Simplified_GNN", "S5", "Simplified_GNN"),
                ("S7_vs_Simplified_GNN", "S7", "Simplified_GNN"),
                ("S5_vs_Simplified_ResMLP", "S5", "Simplified_ResMLP"),
            ]

        for comp, ka, kb in comps:
            # v5
            ci_rows.append({
                "dataset": "nadig2025", "seed": seed, "metric": "v5", "comparison": comp,
                **stratified_row_bootstrap_ci(v5_rows[ka], v5_rows[kb], src_arr),
            })
            # delta_r, disc_raw, disc_v5
            for met in ["delta_r", "disc_raw", "disc_v5"]:
                ci_rows.append({
                    "dataset": "nadig2025", "seed": seed, "metric": met, "comparison": comp,
                    **stratified_row_bootstrap_ci(rm[ka][met], rm[kb][met], src_arr),
                })

        r_disc = [r for r in ci_rows if r["seed"] == seed and r["metric"] == "disc_raw" and r["comparison"] == "S5_vs_AE"][0]
        r_dr = [r for r in ci_rows if r["seed"] == seed and r["metric"] == "delta_r" and r["comparison"] == "S5_vs_AE"][0]
        r_v5 = [r for r in ci_rows if r["seed"] == seed and r["metric"] == "v5" and r["comparison"] == "S5_vs_AE"][0]
        print(
            f"  [Seed {seed:5d}] S5-AE | disc_raw={r_disc['delta']:+.4f} [{r_disc['ci95_lo']:+.4f},{r_disc['ci95_hi']:+.4f}] | "
            f"delta_r={r_dr['delta']:+.4f} | v5={r_v5['delta']:+.4f} [{r_v5['ci95_lo']:+.4f},{r_v5['ci95_hi']:+.4f}] ({time.time()-ts:.1f}s)",
            flush=True,
        )
        gc.collect()

    df_pan = pd.DataFrame(metric_rows)
    df_pan.to_csv(REPORTS_DIR / "r71_nadig_5metric_pan.csv", index=False)
    dl_csv = REPORTS_DIR / "r71_nadig_5metric_dl.csv"
    if dl_csv.exists():
        df_all = pd.concat([df_pan, pd.read_csv(dl_csv)], ignore_index=True)
    else:
        df_all = df_pan.copy()
    df_all.to_csv(REPORTS_DIR / "r71_nadig_5metric_combined.csv", index=False)

    df_ci = pd.DataFrame(ci_rows)
    df_ci.to_csv(REPORTS_DIR / "r71_nadig_row_ci.csv", index=False)

    # Aggregate summary & strict pre-registered H1/H2/H3 adjudication
    summary = {
        "dataset": "Nadig et al. (Nature Genetics 2025, GSE264667) HepG2 + Jurkat",
        "shared_readout_genes": int(len(shared)),
        "n_splits": len(SEEDS_10),
        "split_stats": split_stats,
        "methods": {},
        "comparisons": {},
    }
    for m in df_all["method"].unique():
        sub = df_all[df_all["method"] == m]
        summary["methods"][m] = {}
        for col in ["v5", "delta_r", "de50_r", "de50_sign", "disc_raw", "disc_v5", "discrim_delta_pct", "discrim_v5_pct"]:
            summary["methods"][m][col] = {
                "mean": round(float(sub[col].mean()), 6),
                "std": round(float(sub[col].std(ddof=1)), 6),
            }

    for met in ["v5", "delta_r", "disc_raw", "disc_v5"]:
        summary["comparisons"][met] = {}
        for comp in df_ci["comparison"].unique():
            sub = df_ci[(df_ci["metric"] == met) & (df_ci["comparison"] == comp)]
            summary["comparisons"][met][comp] = {
                "mean_delta": round(float(sub["delta"].mean()), 6),
                "std_delta": round(float(sub["delta"].std(ddof=1)), 6),
                "mean_ci95_lo": round(float(sub["ci95_lo"].mean()), 6),
                "mean_ci95_hi": round(float(sub["ci95_hi"].mean()), 6),
                "splits_ci_gt_zero": int(sub["ci_gt_zero"].sum()),
                "splits_ci_lt_zero": int((sub["ci95_hi"] < 0).sum()),
                "n_splits": int(len(sub)),
            }

    h1_c = summary["comparisons"]["disc_raw"]["S5_vs_AE"]
    h2_c = summary["comparisons"]["delta_r"]["S5_vs_AE"]
    h3_c = summary["comparisons"]["v5"]["S5_vs_AE"]
    h1_pass = bool(h1_c["mean_delta"] > 0.0 and h1_c["splits_ci_gt_zero"] >= 7)
    h2_pass = bool(h2_c["mean_delta"] < 0.0)
    h3_pass = bool(h3_c["splits_ci_gt_zero"] >= 7)
    main_confirmed = bool(h1_pass and h2_pass)

    summary["preregistered_adjudication"] = {
        "H1_disc_raw_S5_vs_AE": {
            "criterion": "mean_delta > 0 AND >= 7/10 splits row-level 95% CI > 0",
            "mean_delta": h1_c["mean_delta"],
            "splits_ci_gt_zero": f"{h1_c['splits_ci_gt_zero']}/10",
            "verdict": "确认 (CONFIRMED)" if h1_pass else "未确认 (NOT CONFIRMED)",
            "passed": h1_pass,
        },
        "H2_delta_r_metric_inversion_S5_vs_AE": {
            "criterion": "mean_delta < 0 (metric inversion on unprojected delta_r)",
            "mean_delta": h2_c["mean_delta"],
            "splits_ci_lt_zero": f"{h2_c['splits_ci_lt_zero']}/10",
            "verdict": "确认 (CONFIRMED)" if h2_pass else "未确认 (NOT CONFIRMED)",
            "passed": h2_pass,
        },
        "MAIN_THESIS_JOINT_H1_AND_H2": {
            "criterion": "BOTH H1 AND H2 must pass",
            "verdict": "确认 (CONFIRMED)" if main_confirmed else "未确认 (NOT CONFIRMED)",
            "passed": main_confirmed,
        },
        "H3_secondary_v5_S5_vs_AE": {
            "criterion": ">= 7/10 splits row-level 95% CI > 0",
            "mean_delta": h3_c["mean_delta"],
            "splits_ci_gt_zero": f"{h3_c['splits_ci_gt_zero']}/10",
            "verdict": "确认 (CONFIRMED)" if h3_pass else "未确认 (NOT CONFIRMED)",
            "passed": h3_pass,
        },
    }

    with open(REPORTS_DIR / "r71_nadig_confirmation_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n=== R7-1 NADIG 2025 PRE-REGISTERED CONFIRMATION SUMMARY ({time.time()-t0:.1f}s) ===", flush=True)
    for k, v in summary["preregistered_adjudication"].items():
        print(f"  {k}: {v['verdict']} ({v})", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["dl", "pan"], required=True)
    args = ap.parse_args()
    if args.mode == "dl":
        run_dl()
    else:
        run_pan()
