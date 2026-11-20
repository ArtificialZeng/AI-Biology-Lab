#!/usr/bin/env python3
"""
T3: Row-Level Stratified Paired Bootstrap 95% Confidence Intervals (4,000 Resamples)
for BioXArena-Pro `xatlas-orion-perturb` and Replogle et al. (Cell 2022).

Replaces invalid split-level paired t-test p-values (whose splits share overlapping training sets)
with exact per-split row-level paired bootstrap 95% CIs stratified by cell line, matching the
independent audit protocol (`rowboot.py`, 4,000 resamples, RandomState(12345)).

Modes (process-isolated to avoid macOS libomp/Accelerate conflicts):
  1. `--mode xatlas`      : xatlas 10 splits (Benchmark A & Strict Benchmark B):
                            S7 vs S0, S7 vs AE, S1 vs S0, S6 vs S5
  2. `--mode replogle_gnn`: Replogle 10 strict union splits (torch-first): saves per-row GNN v5 scores
  3. `--mode replogle`    : Replogle 10 strict union splits (pure NumPy/Sklearn, loads saved GNN rows):
                            S7 vs AE, S5 vs AE, S7 vs GNN, S5 vs GNN, S1 vs S0, S6 vs S5, S7 vs S0
"""

import argparse
import gc
import json
import os
import sys
import time
from pathlib import Path

if "replogle_gnn" in sys.argv:
    import torch
    torch.set_num_threads(4)

import numpy as np
import pandas as pd
from sklearn.decomposition import TruncatedSVD

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from xatlas_panhuman_solver import (
    DATA_DIR,
    REPORTS_DIR,
    SEEDS_10,
    PanPerturbBridgeModel,
    build_spectral_pyramid_Gs,
    fit_baseline_and_axis,
    project_rows_5000,
    score_prediction,
)
from replogle_external_validation import (
    load_replogle_bundle,
    make_strict_union_split,
    run_ahlmann_eltze_bilinear,
)

TMP_GNN_DIR = REPORTS_DIR / "_tmp_replogle_gnn_rows"


def stratified_row_bootstrap_ci(scores_a: np.ndarray, scores_b: np.ndarray, sources: np.ndarray,
                                line0: str = "HCT116", line1: str = "HEK293T",
                                n_boot: int = 4000, seed: int = 12345) -> dict:
    rs = np.random.RandomState(seed)
    d = scores_a - scores_b
    i0 = np.where(sources == line0)[0]
    i1 = np.where(sources == line1)[0]
    point_a = 0.5 * (scores_a[i0].mean() + scores_a[i1].mean())
    point_b = 0.5 * (scores_b[i0].mean() + scores_b[i1].mean())
    point_d = 0.5 * (d[i0].mean() + d[i1].mean())

    # Vectorized 4,000 stratified resamples for speed and exactness
    b0 = d[rs.choice(i0, size=(n_boot, len(i0)), replace=True)].mean(axis=1)
    b1 = d[rs.choice(i1, size=(n_boot, len(i1)), replace=True)].mean(axis=1)
    boots = 0.5 * (b0 + b1)
    ci_lo, ci_hi = np.percentile(boots, [2.5, 97.5])
    return {
        "score_a": round(float(point_a), 6),
        "score_b": round(float(point_b), 6),
        "delta": round(float(point_d), 6),
        "ci95_lo": round(float(ci_lo), 6),
        "ci95_hi": round(float(ci_hi), 6),
        "boot_se": round(float(boots.std()), 6),
        "ci_gt_zero": bool(ci_lo > 0.0),
        "ci_lt_zero": bool(ci_hi < 0.0),
    }


def row_scores_v5(preds: dict, val_df: pd.DataFrame, all_y: dict, bases: dict, axes: dict, score_indices: dict, gene2idx: dict) -> np.ndarray:
    arr = np.empty(len(val_df), dtype=np.float64)
    for i, (_, r) in enumerate(val_df.iterrows()):
        src = r["source"]
        rid = r["id"]
        t_idx = gene2idx[src].get(r["perturbation"])
        arr[i] = score_prediction(
            all_y[rid].astype(np.float64),
            preds[rid].astype(np.float64),
            bases[src], axes[src],
            score_indices[f"{src}_indices"], t_idx,
        )
    return arr


def run_xatlas_bootstrap():
    t0 = time.time()
    print("=== T3: Row-Level Stratified Paired Bootstrap 95% CI on xatlas (10 splits, Bench A & Strict Bench B) ===", flush=True)
    train_df = pd.read_csv(DATA_DIR / "train.csv")
    genes_npz = np.load(DATA_DIR / "genes.npz")
    baselines_pub = {k: np.asarray(v, dtype=np.float64) for k, v in np.load(DATA_DIR / "context_baseline.npz").items()}
    axes_pub = {k: np.asarray(v, dtype=np.float64) for k, v in np.load(DATA_DIR / "residual_axis.npz").items()}
    score_indices = {k: np.asarray(v, dtype=np.int64) for k, v in np.load(DATA_DIR / "score_gene_indices.npz").items()}
    all_train_y = {k: np.asarray(v, dtype=np.float32) for k, v in np.load(DATA_DIR / "train_y.npz").items()}

    hct_genes = list(genes_npz["HCT116_genes"])
    hek_genes = list(genes_npz["HEK293T_genes"])
    shared = sorted(list(set(hct_genes) & set(hek_genes)))
    gene2idx = {
        "HCT116": {g: i for i, g in enumerate(hct_genes)},
        "HEK293T": {g: i for i, g in enumerate(hek_genes)},
    }
    model = PanPerturbBridgeModel(
        gene2idx,
        np.array([hct_genes.index(g) for g in shared]),
        np.array([hek_genes.index(g) for g in shared]),
        score_indices,
    )

    records = []
    for bench in ["BenchA_Standard", "BenchB_StrictUnionRefit"]:
        print(f"\n--- {bench} ---", flush=True)
        for seed in SEEDS_10:
            ts = time.time()
            rng = np.random.RandomState(seed)
            val_rids, vt_line, union_vt = [], {}, set()
            for src, n_val in [("HCT116", 129), ("HEK293T", 247)]:
                sub = train_df[
                    (train_df["source"] == src)
                    & (train_df["reproducible"] == 1)
                    & (train_df["perturbation"].isin(gene2idx[src].keys()))
                ]
                vt = set(rng.permutation(sub["perturbation"].unique())[:n_val])
                vt_line[src] = vt
                union_vt |= vt
                val_rids.extend(sub[sub["perturbation"].isin(vt)]["id"].tolist())
            val_df = train_df[train_df["id"].isin(val_rids)].copy().reset_index(drop=True)
            vr_set = set(val_rids)

            if bench == "BenchA_Standard":
                bases, axes = baselines_pub, axes_pub
            else:
                df_tr_u = train_df[~train_df["perturbation"].isin(union_vt)]
                bases, axes = {}, {}
                for src in ["HCT116", "HEK293T"]:
                    rep_s = df_tr_u[(df_tr_u["source"] == src) & (df_tr_u["reproducible"] == 1)]
                    bases[src], axes[src] = fit_baseline_and_axis(rep_s, all_train_y, score_indices[f"{src}_indices"])

            stages = ("S0", "S1", "S5", "S6", "S7")
            preds = {k: {} for k in list(stages) + ["AE"]}
            for src in ["HCT116", "HEK293T"]:
                oth = "HEK293T" if src == "HCT116" else "HCT116"
                g2i = gene2idx[src]
                if bench == "BenchA_Standard":
                    df_own = train_df[(train_df["source"] == src) & (~train_df["id"].isin(vr_set))]
                    df_oth = train_df[(train_df["source"] == oth) & (~train_df["perturbation"].isin(vt_line[src]))]
                else:
                    df_own = df_tr_u[df_tr_u["source"] == src]
                    df_oth = df_tr_u[df_tr_u["source"] == oth]
                q_df = val_df[val_df["source"] == src]

                res = model.predict_for_source_all_stages(src, df_own, df_oth, q_df, all_train_y, bases, axes, stages)
                for st in stages:
                    preds[st].update(res[st])

                own_rep = df_own[df_own["reproducible"] == 1]
                tr_sub = own_rep[own_rep["perturbation"].isin(g2i.keys())]
                tr_idx = [g2i[p] for p in tr_sub["perturbation"].tolist()]
                q_idx = [g2i[p] for p in q_df["perturbation"].tolist()]
                tr_Y_raw = np.array([all_train_y[rid] for rid in tr_sub["id"]], dtype=np.float32)
                Y_all_raw = np.array([all_train_y[rid] for rid in df_own["id"]], dtype=np.float32)
                p_ae = run_ahlmann_eltze_bilinear(Y_all_raw, tr_idx, q_idx, tr_Y_raw, rank=100, lam=10.0)
                for i, rid in enumerate(q_df["id"].tolist()):
                    preds["AE"][rid] = p_ae[i]

            src_arr = val_df["source"].values
            sc = {k: row_scores_v5(preds[k], val_df, all_train_y, bases, axes, score_indices, gene2idx) for k in preds}
            for comp, ka, kb in [("S7_vs_S0", "S7", "S0"), ("S7_vs_AE", "S7", "AE"), ("S1_vs_S0", "S1", "S0"), ("S5_vs_AE", "S5", "AE"), ("S6_vs_S5", "S6", "S5")]:
                st_ci = stratified_row_bootstrap_ci(sc[ka], sc[kb], src_arr)
                records.append({"dataset": "xatlas", "protocol": bench, "seed": seed, "comparison": comp, **st_ci})
            r_s7s0 = [r for r in records if r["protocol"] == bench and r["seed"] == seed and r["comparison"] == "S7_vs_S0"][0]
            r_s7ae = [r for r in records if r["protocol"] == bench and r["seed"] == seed and r["comparison"] == "S7_vs_AE"][0]
            print(
                f"  Seed {seed:5d} | S7-S0={r_s7s0['delta']:+.4f} [{r_s7s0['ci95_lo']:+.4f}, {r_s7s0['ci95_hi']:+.4f}] "
                f"| S7-AE={r_s7ae['delta']:+.4f} [{r_s7ae['ci95_lo']:+.4f}, {r_s7ae['ci95_hi']:+.4f}] [{time.time()-ts:.1f}s]",
                flush=True,
            )
            gc.collect()

    df_out = pd.DataFrame(records)
    df_out.to_csv(REPORTS_DIR / "row_bootstrap_ci_xatlas.csv", index=False)
    print(f"\nSaved reports/row_bootstrap_ci_xatlas.csv in {time.time()-t0:.1f}s!", flush=True)


def run_replogle_gnn_rows():
    t0 = time.time()
    import xatlas_sota_dl_benchmark as _dl_mod
    _dl_mod.DEVICE = torch.device("cpu")
    from xatlas_sota_dl_benchmark import train_predict_gnn

    print("=== T3 (Replogle GNN Rows): Saving per-row GNN v5 scores across 10 strict splits ===", flush=True)
    os.makedirs(TMP_GNN_DIR, exist_ok=True)
    train_df, all_y, _, _, _, gene2idx, score_indices = load_replogle_bundle()
    for seed in SEEDS_10:
        ts = time.time()
        val_df, df_tr, _, bases, axes = make_strict_union_split(train_df, gene2idx, all_y, score_indices, seed)
        preds_gnn = {}
        for src in ["HCT116", "HEK293T"]:
            g2i = gene2idx[src]
            s_idx = score_indices[f"{src}_indices"]
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
            Gs_own = build_spectral_pyramid_Gs(df_own, all_y, g2i, bases[src], pc_levels=(0, 1))
            G_single = 0.5 * (Gs_own[0] + Gs_own[1])
            p_gnn = train_predict_gnn(X_200, G_single, tr_idx, q_idx, tr_Y_p, bases[src], axes[src], s_idx, seed=seed)
            for i, rid in enumerate(q_df["id"].tolist()):
                preds_gnn[rid] = p_gnn[i]
        sc_gnn = row_scores_v5(preds_gnn, val_df, all_y, bases, axes, score_indices, gene2idx)
        np.save(TMP_GNN_DIR / f"gnn_rows_seed_{seed}.npy", sc_gnn)
        print(f"  [GNN] Seed {seed:5d} saved (macro v5={0.5*(sc_gnn[val_df.source=='HCT116'].mean()+sc_gnn[val_df.source=='HEK293T'].mean()):.4f}) [{time.time()-ts:.1f}s]", flush=True)
        gc.collect()
    print(f"Saved all 10 Replogle GNN per-row scores in {time.time()-t0:.1f}s!", flush=True)


def run_replogle_bootstrap():
    t0 = time.time()
    print("=== T3: Row-Level Stratified Paired Bootstrap 95% CI on Replogle 2022 (10 strict union splits) ===", flush=True)
    train_df, all_y, _, rpe_sh, k562_sh, gene2idx, score_indices = load_replogle_bundle()
    model = PanPerturbBridgeModel(gene2idx, rpe_sh, k562_sh, score_indices)

    records = []
    for seed in SEEDS_10:
        ts = time.time()
        val_df, df_tr, _, bases, axes = make_strict_union_split(train_df, gene2idx, all_y, score_indices, seed)
        stages = ("S0", "S1", "S5", "S6", "S7")
        preds = {k: {} for k in list(stages) + ["AE"]}
        for src in ["HCT116", "HEK293T"]:
            oth = "HEK293T" if src == "HCT116" else "HCT116"
            g2i = gene2idx[src]
            df_own = df_tr[df_tr["source"] == src]
            df_oth = df_tr[df_tr["source"] == oth]
            q_df = val_df[val_df["source"] == src]
            res = model.predict_for_source_all_stages(src, df_own, df_oth, q_df, all_y, bases, axes, stages)
            for st in stages:
                preds[st].update(res[st])

            own_rep = df_own[df_own["reproducible"] == 1]
            tr_sub = own_rep[own_rep["perturbation"].isin(g2i.keys())]
            tr_idx = [g2i[p] for p in tr_sub["perturbation"].tolist()]
            q_idx = [g2i[p] for p in q_df["perturbation"].tolist()]
            tr_Y_raw = np.array([all_y[rid] for rid in tr_sub["id"]], dtype=np.float32)
            Y_all_raw = np.array([all_y[rid] for rid in df_own["id"]], dtype=np.float32)
            p_ae = run_ahlmann_eltze_bilinear(Y_all_raw, tr_idx, q_idx, tr_Y_raw, rank=100, lam=10.0)
            for i, rid in enumerate(q_df["id"].tolist()):
                preds["AE"][rid] = p_ae[i]

        src_arr = val_df["source"].values
        sc = {k: row_scores_v5(preds[k], val_df, all_y, bases, axes, score_indices, gene2idx) for k in preds}
        gnn_path = TMP_GNN_DIR / f"gnn_rows_seed_{seed}.npy"
        if gnn_path.exists():
            sc["GNN"] = np.load(gnn_path)

        pairs = [("S7_vs_AE", "S7", "AE"), ("S5_vs_AE", "S5", "AE"), ("S1_vs_S0", "S1", "S0"), ("S6_vs_S5", "S6", "S5"), ("S7_vs_S0", "S7", "S0")]
        if "GNN" in sc:
            pairs += [("S7_vs_GNN", "S7", "GNN"), ("S5_vs_GNN", "S5", "GNN")]
        for comp, ka, kb in pairs:
            st_ci = stratified_row_bootstrap_ci(sc[ka], sc[kb], src_arr)
            records.append({"dataset": "replogle2022", "protocol": "StrictUnionRefit", "seed": seed, "comparison": comp, **st_ci})
        r7ae = [r for r in records if r["seed"] == seed and r["comparison"] == "S7_vs_AE"][0]
        r65 = [r for r in records if r["seed"] == seed and r["comparison"] == "S6_vs_S5"][0]
        print(
            f"  Seed {seed:5d} | S7-AE={r7ae['delta']:+.4f} [{r7ae['ci95_lo']:+.4f}, {r7ae['ci95_hi']:+.4f}] "
            f"| S6-S5={r65['delta']:+.4f} [{r65['ci95_lo']:+.4f}, {r65['ci95_hi']:+.4f}] [{time.time()-ts:.1f}s]",
            flush=True,
        )
        gc.collect()

    df_out = pd.DataFrame(records)
    df_out.to_csv(REPORTS_DIR / "row_bootstrap_ci_replogle.csv", index=False)

    # Aggregate summary JSON across both xatlas and replogle
    summary = {}
    for fn, key in [(REPORTS_DIR / "row_bootstrap_ci_xatlas.csv", "xatlas"), (REPORTS_DIR / "row_bootstrap_ci_replogle.csv", "replogle2022")]:
        if not fn.exists():
            continue
        df = pd.read_csv(fn)
        summary[key] = {}
        for (proto, comp), grp in df.groupby(["protocol", "comparison"]):
            k = f"{proto}__{comp}"
            summary[key][k] = {
                "mean_delta": round(float(grp["delta"].mean()), 4),
                "median_delta": round(float(grp["delta"].median()), 4),
                "median_ci95_lo": round(float(grp["ci95_lo"].median()), 4),
                "median_ci95_hi": round(float(grp["ci95_hi"].median()), 4),
                "mean_ci95_lo": round(float(grp["ci95_lo"].mean()), 4),
                "mean_ci95_hi": round(float(grp["ci95_hi"].mean()), 4),
                "min_ci95_lo": round(float(grp["ci95_lo"].min()), 4),
                "max_ci95_hi": round(float(grp["ci95_hi"].max()), 4),
                "mean_boot_se": round(float(grp["boot_se"].mean()), 4),
                "splits_ci_strictly_gt_zero": f"{int(grp['ci_gt_zero'].sum())}/{len(grp)}",
                "splits_ci_strictly_lt_zero": f"{int(grp['ci_lt_zero'].sum())}/{len(grp)}",
            }
    with open(REPORTS_DIR / "row_bootstrap_ci_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("\n=== T3 ROW-LEVEL BOOTSTRAP 95% CI SUMMARY ===")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["xatlas", "replogle_gnn", "replogle"], required=True)
    args = ap.parse_args()
    if args.mode == "xatlas":
        run_xatlas_bootstrap()
    elif args.mode == "replogle_gnn":
        run_replogle_gnn_rows()
    else:
        run_replogle_bootstrap()
