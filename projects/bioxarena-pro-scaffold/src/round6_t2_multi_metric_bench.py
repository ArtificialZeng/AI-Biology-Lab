#!/usr/bin/env python3
"""
T2: Comprehensive 5-Metric Benchmark on BOTH xatlas (10 dev + 10 additional resplits)
and Replogle et al. (Cell 2022, 10 strict union splits).

Evaluates ALL methods (Constant_Mean, kNN_30, AE_Bilinear, S0..S7, simplified ResMLP, simplified GNN)
on FIVE metrics simultaneously (including unfavorable metrics in main tables):
  1. v5_orth_r         : Official v5 orthogonal-residual Pearson (span{1, baseline, axis} projected out)
  2. delta_r           : Standard perturbation Delta-expression Pearson (Ahlmann-Eltze et al. 2025)
  3. top50_de_r        : Top-50 differentially expressed gene Pearson (by |y - baseline|)
  4. top50_sign        : Top-50 DE gene sign / direction accuracy
  5. discrim_delta_pct : Perturbation Discrimination Rank Percentile (independent of v5!):
                         For each held-out perturbation i, rank standard Pearson(y_hat_i - b, y_j - b)
                         against all M held-out true perturbations j in the same cell line;
                         percentile = (# of j != i with corr < corr(i,i)) / (M - 1) in [0, 1]
                         (0.50 = random, 1.00 = perfect #1 retrieval; also reports discrim_v5_pct).

Modes (process-isolated to avoid macOS libomp/Accelerate conflicts):
  --mode xatlas_pan    : Constant_Mean, kNN_30, AE_Bilinear, S0..S7 on 10 dev + 10 resplits (20 splits)
  --mode xatlas_dl     : Simplified ResMLP & GNN on 10 dev splits (torch-first)
  --mode replogle_pan  : Constant_Mean, kNN_30, AE_Bilinear, S0..S7 on 10 strict union splits
  --mode replogle_dl   : Simplified ResMLP & GNN on 10 strict union splits (torch-first)
"""

import argparse
import gc
import json
import os
import sys
import time
from pathlib import Path

if "dl" in " ".join(sys.argv):
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
    _project_out,
    build_spectral_pyramid_Gs,
    project_rows_5000,
    score_prediction,
)
from replogle_external_validation import (
    load_replogle_bundle,
    make_strict_union_split,
    run_ahlmann_eltze_bilinear,
    run_knn30,
)

NEW_SEEDS_10 = [101, 102, 103, 104, 105, 106, 107, 108, 109, 110]


def evaluate_5metrics(preds: dict, val_df: pd.DataFrame, all_y: dict, bases: dict, axes: dict,
                      score_indices: dict, gene2idx: dict) -> dict:
    """Compute all 5 metrics macro-averaged across the two cell lines."""
    by_line = {}
    for src in ["HCT116", "HEK293T"]:
        sub = val_df[val_df["source"] == src].reset_index(drop=True)
        m = len(sub)
        base = bases[src]
        ax = axes[src]
        s_idx = score_indices[f"{src}_indices"]

        # Precompute line-level Q for v5 discrimination matrix
        Q_full, Rm = np.linalg.qr(np.stack([np.ones(len(s_idx)), base[s_idx], ax[s_idx]], axis=1))
        Q_full = Q_full[:, np.abs(np.diag(Rm)) > 1e-10 * max(1.0, np.abs(Rm).max())]

        v5_list, dr_list, de50_r_list, de50_s_list = [], [], [], []
        T_delta = np.empty((m, len(s_idx)), dtype=np.float64)
        P_delta = np.empty((m, len(s_idx)), dtype=np.float64)
        T_v5 = np.empty((m, len(s_idx)), dtype=np.float64)
        P_v5 = np.empty((m, len(s_idx)), dtype=np.float64)

        for i, (_, r) in enumerate(sub.iterrows()):
            rid = r["id"]
            t_idx = gene2idx[src].get(r["perturbation"])
            true = all_y[rid].astype(np.float64)
            pred = preds[rid].astype(np.float64)

            v5_list.append(score_prediction(true, pred, base, ax, s_idx, t_idx))

            idx = s_idx[s_idx != t_idx] if (t_idx is not None) else s_idx
            dt = true[idx] - base[idx]
            dp = pred[idx] - base[idx]

            dt_c = dt - dt.mean()
            dp_c = dp - dp.mean()
            denom = np.linalg.norm(dt_c) * np.linalg.norm(dp_c)
            dr_list.append(float(np.clip(np.dot(dt_c, dp_c) / denom, -1.0, 1.0)) if denom > 1e-12 else 0.0)

            top50 = np.argsort(np.abs(dt))[-50:]
            dt50, dp50 = dt[top50], dp[top50]
            dt50_c, dp50_c = dt50 - dt50.mean(), dp50 - dp50.mean()
            d50 = np.linalg.norm(dt50_c) * np.linalg.norm(dp50_c)
            de50_r_list.append(float(np.clip(np.dot(dt50_c, dp50_c) / d50, -1.0, 1.0)) if d50 > 1e-12 else 0.0)
            de50_s_list.append(float(np.mean(np.sign(dt50) == np.sign(dp50))))

            # Store full s_idx vectors (masking target gene i to 0 for fair cross-perturbation retrieval)
            t_vec = true[s_idx] - base[s_idx]
            p_vec = pred[s_idx] - base[s_idx]
            if t_idx is not None:
                pos = np.where(s_idx == t_idx)[0]
                if len(pos) > 0:
                    t_vec[pos[0]] = 0.0
                    p_vec[pos[0]] = 0.0
            T_delta[i] = t_vec - t_vec.mean()
            P_delta[i] = p_vec - p_vec.mean()
            T_v5[i] = _project_out(t_vec, Q_full)
            P_v5[i] = _project_out(p_vec, Q_full)

        # Perturbation Discrimination Rank Percentile (standard unprojected Delta & v5-projected)
        def calc_discrim(P_mat, T_mat):
            Pn = P_mat / np.maximum(np.linalg.norm(P_mat, axis=1, keepdims=True), 1e-12)
            Tn = T_mat / np.maximum(np.linalg.norm(T_mat, axis=1, keepdims=True), 1e-12)
            # Zero-norm predictions (e.g. Constant_Mean) get random chance 0.50
            p_zero = np.linalg.norm(P_mat, axis=1) < 1e-10
            S = Pn @ Tn.T  # (m, m): row i = pred i vs all true j
            pcts = np.empty(m, dtype=np.float64)
            top1 = np.empty(m, dtype=np.float64)
            for i in range(m):
                if p_zero[i]:
                    pcts[i] = 0.5
                    top1[i] = 1.0 / m
                else:
                    diag = S[i, i]
                    others = np.delete(S[i], i)
                    # Tie-aware percentile: strictly less + 0.5 * equal
                    pcts[i] = (np.sum(others < diag) + 0.5 * np.sum(others == diag)) / (m - 1)
                    top1[i] = 1.0 if np.all(diag > others) else 0.0
            return float(pcts.mean()), float(top1.mean())

        dd_pct, dd_top1 = calc_discrim(P_delta, T_delta)
        dv_pct, dv_top1 = calc_discrim(P_v5, T_v5)

        by_line[src] = {
            "v5": float(np.mean(v5_list)),
            "delta_r": float(np.mean(dr_list)),
            "de50_r": float(np.mean(de50_r_list)),
            "de50_sign": float(np.mean(de50_s_list)),
            "discrim_delta_pct": dd_pct,
            "discrim_delta_top1": dd_top1,
            "discrim_v5_pct": dv_pct,
            "discrim_v5_top1": dv_top1,
        }

    out = {}
    for k in by_line["HCT116"].keys():
        out[k] = 0.5 * (by_line["HCT116"][k] + by_line["HEK293T"][k])
        out[f"{k}_L0"] = by_line["HCT116"][k]
        out[f"{k}_L1"] = by_line["HEK293T"][k]
    return out


def load_xatlas_bundle():
    train_df = pd.read_csv(DATA_DIR / "train.csv")
    genes_npz = np.load(DATA_DIR / "genes.npz")
    baselines_pub = {k: np.asarray(v, dtype=np.float64) for k, v in np.load(DATA_DIR / "context_baseline.npz").items()}
    axes_pub = {k: np.asarray(v, dtype=np.float64) for k, v in np.load(DATA_DIR / "residual_axis.npz").items()}
    score_indices = {k: np.asarray(v, dtype=np.int64) for k, v in np.load(DATA_DIR / "score_gene_indices.npz").items()}
    all_train_y = {k: np.asarray(v, dtype=np.float32) for k, v in np.load(DATA_DIR / "train_y.npz").items()}
    hct_genes = list(genes_npz["HCT116_genes"])
    hek_genes = list(genes_npz["HEK293T_genes"])
    shared = sorted(list(set(hct_genes) & set(hek_genes)))
    hct_sh = np.array([hct_genes.index(g) for g in shared])
    hek_sh = np.array([hek_genes.index(g) for g in shared])
    gene2idx = {
        "HCT116": {g: i for i, g in enumerate(hct_genes)},
        "HEK293T": {g: i for i, g in enumerate(hek_genes)},
    }
    return train_df, all_train_y, baselines_pub, axes_pub, score_indices, hct_sh, hek_sh, gene2idx


def make_xatlas_split(train_df, gene2idx, seed):
    rng = np.random.RandomState(seed)
    val_rids, vt_line = [], {}
    for src, n_val in [("HCT116", 129), ("HEK293T", 247)]:
        sub = train_df[
            (train_df["source"] == src)
            & (train_df["reproducible"] == 1)
            & (train_df["perturbation"].isin(gene2idx[src].keys()))
        ]
        vt = set(rng.permutation(sub["perturbation"].unique())[:n_val])
        vt_line[src] = vt
        val_rids.extend(sub[sub["perturbation"].isin(vt)]["id"].tolist())
    val_df = train_df[train_df["id"].isin(val_rids)].copy().reset_index(drop=True)
    return val_df, set(val_rids), vt_line


def run_xatlas_pan():
    t0 = time.time()
    print("=== T2 (xatlas PAN): 5-Metric Benchmark on 10 Dev Splits + 10 Additional Resplits ===", flush=True)
    train_df, all_y, bases, axes, sidx, hct_sh, hek_sh, g2i_map = load_xatlas_bundle()
    model = PanPerturbBridgeModel(g2i_map, hct_sh, hek_sh, sidx)
    stages = ["S0", "S1", "S2", "S3", "S4", "S5", "S6", "S7"]
    methods = ["Constant_Mean", "kNN_30", "AE_Bilinear"] + stages

    rows = []
    for group, seeds in [("dev10", SEEDS_10), ("resplit10", NEW_SEEDS_10)]:
        for seed in seeds:
            ts = time.time()
            val_df, vr_set, vt_line = make_xatlas_split(train_df, g2i_map, seed)
            preds = {m: {} for m in methods}
            for src in ["HCT116", "HEK293T"]:
                oth = "HEK293T" if src == "HCT116" else "HCT116"
                g2i = g2i_map[src]
                df_own = train_df[(train_df["source"] == src) & (~train_df["id"].isin(vr_set))]
                df_oth = train_df[(train_df["source"] == oth) & (~train_df["perturbation"].isin(vt_line[src]))]
                q_df = val_df[val_df["source"] == src]
                q_idx = [g2i[p] for p in q_df["perturbation"].tolist()]
                q_rids = q_df["id"].tolist()

                res = model.predict_for_source_all_stages(src, df_own, df_oth, q_df, all_y, bases, axes, stages)
                for st in stages:
                    preds[st].update(res[st])

                own_rep = df_own[df_own["reproducible"] == 1]
                tr_sub = own_rep[own_rep["perturbation"].isin(g2i.keys())]
                tr_idx = [g2i[p] for p in tr_sub["perturbation"].tolist()]
                tr_Y_raw = np.array([all_y[rid] for rid in tr_sub["id"]], dtype=np.float32)
                Y_all_raw = np.array([all_y[rid] for rid in df_own["id"]], dtype=np.float32)

                mean_v = tr_Y_raw.mean(axis=0).astype(np.float64)
                p_knn = run_knn30(Y_all_raw, tr_idx, q_idx, tr_Y_raw, k=30)
                p_ae = run_ahlmann_eltze_bilinear(Y_all_raw, tr_idx, q_idx, tr_Y_raw, rank=100, lam=10.0)
                for i, rid in enumerate(q_rids):
                    preds["Constant_Mean"][rid] = mean_v
                    preds["kNN_30"][rid] = p_knn[i]
                    preds["AE_Bilinear"][rid] = p_ae[i]

            for m in methods:
                ev = evaluate_5metrics(preds[m], val_df, all_y, bases, axes, sidx, g2i_map)
                rows.append({"dataset": "xatlas", "split_group": group, "seed": seed, "method": m, **{k: round(v, 6) for k, v in ev.items()}})
            s7_r = [r for r in rows if r["seed"] == seed and r["method"] == "S7"][0]
            ae_r = [r for r in rows if r["seed"] == seed and r["method"] == "AE_Bilinear"][0]
            print(
                f"  [{group:9s}] Seed {seed:5d} | v5: AE={ae_r['v5']:.4f} S7={s7_r['v5']:.4f} | "
                f"DeltaR: AE={ae_r['delta_r']:.4f} S7={s7_r['delta_r']:.4f} | "
                f"DE50: AE={ae_r['de50_r']:.4f} S7={s7_r['de50_r']:.4f} | "
                f"Discrim(Delta): AE={ae_r['discrim_delta_pct']:.4f} S7={s7_r['discrim_delta_pct']:.4f} [{time.time()-ts:.1f}s]",
                flush=True,
            )
            gc.collect()

    pd.DataFrame(rows).to_csv(REPORTS_DIR / "t2_xatlas_5metric_pan.csv", index=False)
    print(f"Saved t2_xatlas_5metric_pan.csv in {time.time()-t0:.1f}s!", flush=True)


def run_xatlas_dl():
    t0 = time.time()
    import xatlas_sota_dl_benchmark as _dl_mod
    _dl_mod.DEVICE = torch.device("cpu")
    from xatlas_sota_dl_benchmark import train_predict_gnn, train_predict_mlp

    print("=== T2 (xatlas DL): Simplified ResMLP & GNN on 10 Dev Splits (5 Metrics) ===", flush=True)
    train_df, all_y, bases, axes, sidx, _, _, g2i_map = load_xatlas_bundle()
    rows = []
    for seed in SEEDS_10:
        ts = time.time()
        val_df, vr_set, _ = make_xatlas_split(train_df, g2i_map, seed)
        preds = {"Simplified_ResMLP": {}, "Simplified_GNN": {}}
        for src in ["HCT116", "HEK293T"]:
            g2i = g2i_map[src]
            s_i = sidx[f"{src}_indices"]
            df_own = train_df[(train_df["source"] == src) & (~train_df["id"].isin(vr_set))]
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
        for m in preds:
            ev = evaluate_5metrics(preds[m], val_df, all_y, bases, axes, sidx, g2i_map)
            rows.append({"dataset": "xatlas", "split_group": "dev10", "seed": seed, "method": m, **{k: round(v, 6) for k, v in ev.items()}})
        print(f"  [xatlas_dl] Seed {seed:5d} done [{time.time()-ts:.1f}s]", flush=True)
        gc.collect()
    pd.DataFrame(rows).to_csv(REPORTS_DIR / "t2_xatlas_5metric_dl.csv", index=False)


def run_replogle_pan():
    t0 = time.time()
    print("=== T2 (Replogle PAN): 5-Metric Benchmark on 10 Strict Union Splits ===", flush=True)
    train_df, all_y, _, rpe_sh, k562_sh, g2i_map, sidx = load_replogle_bundle()
    model = PanPerturbBridgeModel(g2i_map, rpe_sh, k562_sh, sidx)
    stages = ["S0", "S1", "S2", "S3", "S4", "S5", "S6", "S7"]
    methods = ["Constant_Mean", "kNN_30", "AE_Bilinear"] + stages
    rows = []
    for seed in SEEDS_10:
        ts = time.time()
        val_df, df_tr, _, bases, axes = make_strict_union_split(train_df, g2i_map, all_y, sidx, seed)
        preds = {m: {} for m in methods}
        for src in ["HCT116", "HEK293T"]:
            oth = "HEK293T" if src == "HCT116" else "HCT116"
            g2i = g2i_map[src]
            df_own = df_tr[df_tr["source"] == src]
            df_oth = df_tr[df_tr["source"] == oth]
            q_df = val_df[val_df["source"] == src]
            q_idx = [g2i[p] for p in q_df["perturbation"].tolist()]
            q_rids = q_df["id"].tolist()

            res = model.predict_for_source_all_stages(src, df_own, df_oth, q_df, all_y, bases, axes, stages)
            for st in stages:
                preds[st].update(res[st])

            own_rep = df_own[df_own["reproducible"] == 1]
            tr_sub = own_rep[own_rep["perturbation"].isin(g2i.keys())]
            tr_idx = [g2i[p] for p in tr_sub["perturbation"].tolist()]
            tr_Y_raw = np.array([all_y[rid] for rid in tr_sub["id"]], dtype=np.float32)
            Y_all_raw = np.array([all_y[rid] for rid in df_own["id"]], dtype=np.float32)
            mean_v = tr_Y_raw.mean(axis=0).astype(np.float64)
            p_knn = run_knn30(Y_all_raw, tr_idx, q_idx, tr_Y_raw, k=30)
            p_ae = run_ahlmann_eltze_bilinear(Y_all_raw, tr_idx, q_idx, tr_Y_raw, rank=100, lam=10.0)
            for i, rid in enumerate(q_rids):
                preds["Constant_Mean"][rid] = mean_v
                preds["kNN_30"][rid] = p_knn[i]
                preds["AE_Bilinear"][rid] = p_ae[i]
        for m in methods:
            ev = evaluate_5metrics(preds[m], val_df, all_y, bases, axes, sidx, g2i_map)
            rows.append({"dataset": "replogle2022", "split_group": "strict10", "seed": seed, "method": m, **{k: round(v, 6) for k, v in ev.items()}})
        s5_r = [r for r in rows if r["seed"] == seed and r["method"] == "S5"][0]
        s7_r = [r for r in rows if r["seed"] == seed and r["method"] == "S7"][0]
        ae_r = [r for r in rows if r["seed"] == seed and r["method"] == "AE_Bilinear"][0]
        print(
            f"  [replogle] Seed {seed:5d} | v5: AE={ae_r['v5']:.4f} S5={s5_r['v5']:.4f} S7={s7_r['v5']:.4f} | "
            f"Discrim(Delta): AE={ae_r['discrim_delta_pct']:.4f} S5={s5_r['discrim_delta_pct']:.4f} S7={s7_r['discrim_delta_pct']:.4f} | "
            f"Discrim(v5): AE={ae_r['discrim_v5_pct']:.4f} S7={s7_r['discrim_v5_pct']:.4f} [{time.time()-ts:.1f}s]",
            flush=True,
        )
        gc.collect()
    pd.DataFrame(rows).to_csv(REPORTS_DIR / "t2_replogle_5metric_pan.csv", index=False)


def run_replogle_dl():
    t0 = time.time()
    import xatlas_sota_dl_benchmark as _dl_mod
    _dl_mod.DEVICE = torch.device("cpu")
    from xatlas_sota_dl_benchmark import train_predict_gnn, train_predict_mlp

    print("=== T2 (Replogle DL): Simplified ResMLP & GNN on 10 Strict Union Splits (5 Metrics) ===", flush=True)
    train_df, all_y, _, _, _, g2i_map, sidx = load_replogle_bundle()
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
        for m in preds:
            ev = evaluate_5metrics(preds[m], val_df, all_y, bases, axes, sidx, g2i_map)
            rows.append({"dataset": "replogle2022", "split_group": "strict10", "seed": seed, "method": m, **{k: round(v, 6) for k, v in ev.items()}})
        print(f"  [replogle_dl] Seed {seed:5d} done [{time.time()-ts:.1f}s]", flush=True)
        gc.collect()
    pd.DataFrame(rows).to_csv(REPORTS_DIR / "t2_replogle_5metric_dl.csv", index=False)


def build_summary():
    dfs = []
    for fn in ["t2_xatlas_5metric_pan.csv", "t2_xatlas_5metric_dl.csv", "t2_replogle_5metric_pan.csv", "t2_replogle_5metric_dl.csv"]:
        p = REPORTS_DIR / fn
        if p.exists():
            dfs.append(pd.read_csv(p))
    df_all = pd.concat(dfs, ignore_index=True)
    df_all.to_csv(REPORTS_DIR / "t2_all_5metric_combined.csv", index=False)

    summary = {}
    metrics = ["v5", "delta_r", "de50_r", "de50_sign", "discrim_delta_pct", "discrim_v5_pct"]
    for (ds, grp), sub in df_all.groupby(["dataset", "split_group"]):
        key = f"{ds}__{grp}"
        summary[key] = {}
        for m, mdf in sub.groupby("method"):
            summary[key][m] = {f"{k}_mean": round(float(mdf[k].mean()), 4) for k in metrics}
            for k in metrics:
                summary[key][m][f"{k}_sd"] = round(float(mdf[k].std()), 4)
    with open(REPORTS_DIR / "t2_all_5metric_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("\n=== T2 5-METRIC SUMMARY ===")
    for k, md in summary.items():
        print(f"\n[{k}]")
        for m in ["Constant_Mean", "kNN_30", "AE_Bilinear", "Simplified_ResMLP", "Simplified_GNN", "S0", "S1", "S2", "S3", "S4", "S5", "S6", "S7"]:
            if m in md:
                r = md[m]
                print(
                    f"  {m:18s} | v5={r['v5_mean']:.4f} | DeltaR={r['delta_r_mean']:.4f} | DE50={r['de50_r_mean']:.4f} "
                    f"| Sign={r['de50_sign_mean']:.4f} | Discrim(D)={r['discrim_delta_pct_mean']:.4f} | Discrim(v5)={r['discrim_v5_pct_mean']:.4f}"
                )


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["xatlas_pan", "xatlas_dl", "replogle_pan", "replogle_dl", "summarize"], required=True)
    args = ap.parse_args()
    if args.mode == "xatlas_pan":
        run_xatlas_pan()
    elif args.mode == "xatlas_dl":
        run_xatlas_dl()
    elif args.mode == "replogle_pan":
        run_replogle_pan()
    elif args.mode == "replogle_dl":
        run_replogle_dl()
    else:
        build_summary()
