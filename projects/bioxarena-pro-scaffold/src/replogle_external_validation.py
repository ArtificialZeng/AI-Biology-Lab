#!/usr/bin/env python3
"""
Zero-Shot External Validation of Frozen PanPerturb-Bridge on
Replogle et al. (Cell 2022) Genome-Scale Perturb-seq (K562 & RPE1).

CRITICAL SCIENTIFIC RIGOR:
  - ZERO hyperparameter retuning: imports the EXACT frozen `PanPerturbBridgeModel`
    class from `xatlas_panhuman_solver.py` (hyperparameters frozen on Orion HCT116/HEK293T).
  - Strict simultaneous union holdout: 20% of unique strong in-panel targets per line
    are removed simultaneously from BOTH K562 and RPE1 (zero target leakage across lines).
  - Strict training-only fitting of `context_baseline` and `residual_axis` per split.
  - Evaluated across 10 independent random splits on THREE standard field metrics:
      (1) BioXArena v5 Orthogonal-Residual Pearson (span{1, baseline, axis} projected out)
      (2) Standard Perturbation Delta-Expression Pearson (Ahlmann-Eltze Nature Methods 2025)
      (3) Top-50 Differentially Expressed (DE) Gene Pearson & Sign Accuracy (GEARS/scGPT/STATE)
  - Supports `--mode pan` (pure NumPy/SciPy/Sklearn: Constant_Mean, kNN_30, AhlmannEltze_Bilinear,
    PanPerturb-Bridge S0..S7) and `--mode dl` (dedicated PyTorch MPS process: ResMLP, GEARS_GNN,
    merging into the final unified CSV & JSON to avoid macOS OpenMP/Accelerate threadpool conflicts).
"""

import argparse
import gc
import json
import os
import sys
import time
from pathlib import Path

# Initialize PyTorch before NumPy/SciPy/Sklearn in DL mode to prevent macOS libomp/Accelerate init conflict
if "dl" in sys.argv:
    import torch
    torch.set_num_threads(4)

import numpy as np
import pandas as pd
from scipy.stats import ttest_rel
from sklearn.decomposition import TruncatedSVD

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from xatlas_panhuman_solver import (
    REPORTS_DIR,
    SEEDS_10,
    PanPerturbBridgeModel,
    build_spectral_pyramid_Gs,
    fit_baseline_and_axis,
    project_rows_5000,
    score_prediction,
)

def _resolve_replogle_dir() -> Path:
    cands = [
        *( [Path(os.environ["REPLOGLE_DATA_DIR"])] if os.environ.get("REPLOGLE_DATA_DIR") else [] ),
        Path.cwd() / "external_data/replogle2022",
        Path(__file__).resolve().parents[4] / "external_data/replogle2022",
        Path(__file__).resolve().parents[3] / "external_data/replogle2022",
    ]
    return next((c for c in cands if c.exists()), cands[0])


REPLOGLE_DIR = _resolve_replogle_dir()


def load_replogle_cell_line(h5_path: Path, source_name: str, n_genes: int = 5000):
    """Load pre-extracted Replogle 2022 top-5000 expressed genes cache."""
    line_tag = "RPE1" if "rpe1" in h5_path.name.lower() else "K562"
    cache_dir = h5_path.parent / "cache"
    meta_df = pd.read_csv(cache_dir / f"{line_tag}_meta.csv")
    genes_df = pd.read_csv(cache_dir / f"{line_tag}_genes.csv")
    X_mat = np.load(cache_dir / f"{line_tag}_X.npy")
    meta_df["source"] = source_name
    meta_df["id"] = [f"{source_name}_{ens}" for ens in meta_df["perturbation"]]
    y_dict = {rid: X_mat[i] for i, rid in enumerate(meta_df["id"])}
    gids = genes_df["gene_id"].tolist()
    gnames = genes_df["gene_name"].tolist()
    return meta_df, y_dict, gids, dict(zip(gids, gnames))


def run_ahlmann_eltze_bilinear(Y_all_raw, tr_idx, query_idx, tr_Y_raw, rank=100, lam=10.0):
    y_mean = tr_Y_raw.mean(axis=0, keepdims=True)
    Y_c = tr_Y_raw - y_mean
    G = TruncatedSVD(n_components=rank, random_state=42).fit_transform(Y_all_raw.T)
    r_p = min(rank, len(tr_Y_raw) - 1)
    svd_p = TruncatedSVD(n_components=r_p, random_state=42).fit(Y_c)
    P = svd_p.components_.T
    G_tr = G[tr_idx]
    A = G_tr.T @ G_tr + lam * np.eye(G_tr.shape[1])
    W = np.linalg.solve(A, G_tr.T @ (Y_c @ P))
    return y_mean + (G[query_idx] @ W) @ P.T


def run_knn30(Y_all_raw, tr_idx, query_idx, tr_Y_raw, k=30):
    G = Y_all_raw.T.astype(np.float64)
    G_u = G / (np.linalg.norm(G, axis=1, keepdims=True) + 1e-12)
    sim = G_u[query_idx] @ G_u[tr_idx].T
    out = np.zeros((len(query_idx), tr_Y_raw.shape[1]), dtype=np.float64)
    for i in range(len(query_idx)):
        top = np.argsort(sim[i])[-k:]
        w = np.maximum(0.0, sim[i, top])
        out[i] = tr_Y_raw[top].mean(axis=0) if w.sum() < 1e-12 else (w[:, None] * tr_Y_raw[top]).sum(axis=0) / w.sum()
    return out


def evaluate_replogle_multi_metric(preds, val_df, all_y, baselines, axes, score_indices, gene2idx):
    metrics_by_line = {s: {"v5": [], "delta": [], "de50_r": [], "de50_sign": []} for s in ["HCT116", "HEK293T"]}
    for _, r in val_df.iterrows():
        rid = r["id"]
        src = r["source"]
        t_idx = gene2idx[src].get(r["perturbation"])
        true = all_y[rid].astype(np.float64)
        pred = preds[rid].astype(np.float64)
        base = baselines[src]
        ax = axes[src]
        s_idx = score_indices[f"{src}_indices"]

        v5 = score_prediction(true, pred, base, ax, s_idx, t_idx)
        metrics_by_line[src]["v5"].append(v5)

        idx = s_idx[s_idx != t_idx] if (t_idx is not None) else s_idx
        dt = true[idx] - base[idx]
        dp = pred[idx] - base[idx]

        dt_c = dt - dt.mean()
        dp_c = dp - dp.mean()
        denom = np.linalg.norm(dt_c) * np.linalg.norm(dp_c)
        dr = float(np.clip(np.dot(dt_c, dp_c) / denom, -1.0, 1.0)) if denom > 1e-12 else 0.0
        metrics_by_line[src]["delta"].append(dr)

        top50 = np.argsort(np.abs(dt))[-50:]
        dt50, dp50 = dt[top50], dp[top50]
        dt50_c, dp50_c = dt50 - dt50.mean(), dp50 - dp50.mean()
        denom50 = np.linalg.norm(dt50_c) * np.linalg.norm(dp50_c)
        r50 = float(np.clip(np.dot(dt50_c, dp50_c) / denom50, -1.0, 1.0)) if denom50 > 1e-12 else 0.0
        sign50 = float(np.mean(np.sign(dt50) == np.sign(dp50)))
        metrics_by_line[src]["de50_r"].append(r50)
        metrics_by_line[src]["de50_sign"].append(sign50)

    out = {}
    for k in ["v5", "delta", "de50_r", "de50_sign"]:
        rpe_m = float(np.mean(metrics_by_line["HCT116"][k]))
        k562_m = float(np.mean(metrics_by_line["HEK293T"][k]))
        out[k] = 0.5 * (rpe_m + k562_m)
        out[f"{k}_RPE1"] = rpe_m
        out[f"{k}_K562"] = k562_m
    return out


def make_strict_union_split(train_df, gene2idx, all_y, score_indices, seed):
    rng = np.random.RandomState(seed)
    val_rids = []
    union_val_targets = set()
    for src, frac in [("HCT116", 0.20), ("HEK293T", 0.20)]:
        sub = train_df[
            (train_df["source"] == src)
            & (train_df["reproducible"] == 1)
            & (train_df["perturbation"].isin(gene2idx[src].keys()))
        ]
        uniq = sub["perturbation"].unique()
        n_val = int(round(len(uniq) * frac))
        vt = set(rng.permutation(uniq)[:n_val])
        union_val_targets |= vt
        val_rids.extend(sub[sub["perturbation"].isin(vt)]["id"].tolist())
    val_df = train_df[train_df["id"].isin(val_rids)].copy().reset_index(drop=True)
    df_tr = train_df[~train_df["perturbation"].isin(union_val_targets)].copy()
    assert len(set(df_tr["perturbation"]) & union_val_targets) == 0
    bases, axes = {}, {}
    for src in ["HCT116", "HEK293T"]:
        rep_s = df_tr[(df_tr["source"] == src) & (df_tr["reproducible"] == 1)]
        bases[src], axes[src] = fit_baseline_and_axis(rep_s, all_y, score_indices[f"{src}_indices"])
    return val_df, df_tr, union_val_targets, bases, axes


def load_replogle_bundle():
    df_rpe, y_rpe, rpe_genes, _ = load_replogle_cell_line(REPLOGLE_DIR / "rpe1_normalized_bulk_01.h5ad", "HCT116", 5000)
    df_k562, y_k562, k562_genes, _ = load_replogle_cell_line(REPLOGLE_DIR / "K562_gwps_normalized_bulk_01.h5ad", "HEK293T", 5000)
    train_df = pd.concat([df_rpe, df_k562], ignore_index=True)
    all_y = {**y_rpe, **y_k562}
    shared_genes = sorted(list(set(rpe_genes) & set(k562_genes)))
    rpe_sh_idx = np.array([rpe_genes.index(g) for g in shared_genes])
    k562_sh_idx = np.array([k562_genes.index(g) for g in shared_genes])
    gene2idx = {
        "HCT116": {g: i for i, g in enumerate(rpe_genes)},
        "HEK293T": {g: i for i, g in enumerate(k562_genes)},
    }
    score_indices = {
        "HCT116_indices": np.arange(5000, dtype=np.int64),
        "HEK293T_indices": np.arange(5000, dtype=np.int64),
    }
    return train_df, all_y, shared_genes, rpe_sh_idx, k562_sh_idx, gene2idx, score_indices


def run_pan_mode():
    t_start = time.time()
    print("=" * 90)
    print("REPLOGLE 2022 EXTERNAL VALIDATION — PHASE 1 (Frozen PanPerturb-Bridge S0..S7 + Linear/kNN)")
    print("=" * 90, flush=True)
    train_df, all_y, shared_genes, rpe_sh_idx, k562_sh_idx, gene2idx, score_indices = load_replogle_bundle()
    print(f"Loaded Replogle 2022: total={len(train_df)} rows | Shared genes={len(shared_genes)}/5000", flush=True)

    frozen_model = PanPerturbBridgeModel(gene2idx, rpe_sh_idx, k562_sh_idx, score_indices)
    stages = ["S0", "S1", "S2", "S3", "S4", "S5", "S6", "S7"]
    methods = ["Constant_Mean", "kNN_30", "AhlmannEltze_Bilinear"] + stages

    records = []
    for seed in SEEDS_10:
        t0 = time.time()
        val_df, df_tr, union_val, bases, axes = make_strict_union_split(train_df, gene2idx, all_y, score_indices, seed)
        preds_all = {m: {} for m in methods}
        for src in ["HCT116", "HEK293T"]:
            oth = "HEK293T" if src == "HCT116" else "HCT116"
            g2i = gene2idx[src]
            df_own = df_tr[df_tr["source"] == src]
            df_oth = df_tr[df_tr["source"] == oth]
            q_df = val_df[val_df["source"] == src]
            query_idx = [g2i[p] for p in q_df["perturbation"].tolist()]
            query_rids = q_df["id"].tolist()

            res_pan = frozen_model.predict_for_source_all_stages(src, df_own, df_oth, q_df, all_y, bases, axes, stages)
            for st in stages:
                preds_all[st].update(res_pan[st])

            own_rep = df_own[df_own["reproducible"] == 1]
            tr_sub = own_rep[own_rep["perturbation"].isin(g2i.keys())]
            tr_idx = [g2i[p] for p in tr_sub["perturbation"].tolist()]
            tr_Y_raw = np.array([all_y[rid] for rid in tr_sub["id"]], dtype=np.float32)
            Y_all_raw = np.array([all_y[rid] for rid in df_own["id"]], dtype=np.float32)

            mean_vec = tr_Y_raw.mean(axis=0).astype(np.float64)
            p_knn = run_knn30(Y_all_raw, tr_idx, query_idx, tr_Y_raw, k=30)
            p_ae = run_ahlmann_eltze_bilinear(Y_all_raw, tr_idx, query_idx, tr_Y_raw, rank=100, lam=10.0)
            for i, rid in enumerate(query_rids):
                preds_all["Constant_Mean"][rid] = mean_vec
                preds_all["kNN_30"][rid] = p_knn[i]
                preds_all["AhlmannEltze_Bilinear"][rid] = p_ae[i]
            del Y_all_raw

        row = {"seed": seed, "n_union_heldout": len(union_val)}
        for m in methods:
            ev = evaluate_replogle_multi_metric(preds_all[m], val_df, all_y, bases, axes, score_indices, gene2idx)
            for k, v in ev.items():
                row[f"{m}_{k}"] = round(v, 6)
        records.append(row)
        print(
            f"  [PAN] Seed {seed:5d} (held={len(union_val)}) | v5: SVD={row['S0_v5']:.4f} AE={row['AhlmannEltze_Bilinear_v5']:.4f} "
            f"-> S2={row['S2_v5']:.4f} -> S5={row['S5_v5']:.4f} -> S6={row['S6_v5']:.4f} -> S7={row['S7_v5']:.4f} "
            f"(RPE1:{row['S7_v5_RPE1']:.4f}, K562:{row['S7_v5_K562']:.4f}) | "
            f"DeltaR:{row['S0_delta']:.4f}->{row['S7_delta']:.4f} | DE50-R:{row['S0_de50_r']:.4f}->{row['S7_de50_r']:.4f} [{time.time()-t0:.1f}s]",
            flush=True,
        )
        gc.collect()

    df_pan = pd.DataFrame(records)
    df_pan.to_csv(REPORTS_DIR / "replogle2022_external_validation_10split.csv", index=False)
    print(f"Phase 1 saved in {time.time()-t_start:.1f}s!", flush=True)


def run_dl_mode():
    t_start = time.time()
    # Import PyTorch ONLY in dedicated DL process using thread-bounded CPU to avoid macOS Metal/Accelerate collision
    import torch
    torch.set_num_threads(4)
    import xatlas_sota_dl_benchmark as _dl_mod
    _dl_mod.DEVICE = torch.device("cpu")
    DEVICE = _dl_mod.DEVICE
    from xatlas_sota_dl_benchmark import train_predict_gnn, train_predict_mlp

    print("=" * 90)
    print(f"REPLOGLE 2022 EXTERNAL VALIDATION — PHASE 2 (PyTorch ResMLP & GEARS-GNN | Device: {DEVICE})")
    print("=" * 90, flush=True)
    train_df, all_y, _, _, _, gene2idx, score_indices = load_replogle_bundle()
    df_base = pd.read_csv(REPORTS_DIR / "replogle2022_external_validation_10split.csv").set_index("seed")

    dl_methods = ["PyTorch_ResMLP", "PyTorch_GEARS_GNN"]
    for seed in SEEDS_10:
        t0 = time.time()
        val_df, df_tr, union_val, bases, axes = make_strict_union_split(train_df, gene2idx, all_y, score_indices, seed)
        preds_dl = {m: {} for m in dl_methods}
        for src in ["HCT116", "HEK293T"]:
            g2i = gene2idx[src]
            s_idx = score_indices[f"{src}_indices"]
            df_own = df_tr[df_tr["source"] == src]
            q_df = val_df[val_df["source"] == src]
            query_idx = [g2i[p] for p in q_df["perturbation"].tolist()]
            query_rids = q_df["id"].tolist()

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

            p_mlp = train_predict_mlp(X_200, tr_idx, query_idx, tr_Y_p, bases[src], axes[src], s_idx, seed=seed)
            Gs_own = build_spectral_pyramid_Gs(df_own, all_y, g2i, bases[src], pc_levels=(0, 1))
            G_single = 0.5 * (Gs_own[0] + Gs_own[1])
            p_gnn = train_predict_gnn(X_200, G_single, tr_idx, query_idx, tr_Y_p, bases[src], axes[src], s_idx, seed=seed)

            for i, rid in enumerate(query_rids):
                preds_dl["PyTorch_ResMLP"][rid] = p_mlp[i]
                preds_dl["PyTorch_GEARS_GNN"][rid] = p_gnn[i]
            del Y_all_raw, Y_clean, Gs_own, G_single

        for m in dl_methods:
            ev = evaluate_replogle_multi_metric(preds_dl[m], val_df, all_y, bases, axes, score_indices, gene2idx)
            for k, v in ev.items():
                df_base.loc[seed, f"{m}_{k}"] = round(v, 6)
        print(
            f"  [DL ] Seed {seed:5d} | v5: ResMLP={df_base.loc[seed, 'PyTorch_ResMLP_v5']:.4f}, "
            f"GEARS-GNN={df_base.loc[seed, 'PyTorch_GEARS_GNN_v5']:.4f} | "
            f"DeltaR: MLP={df_base.loc[seed, 'PyTorch_ResMLP_delta']:.4f}, GNN={df_base.loc[seed, 'PyTorch_GEARS_GNN_delta']:.4f} [{time.time()-t0:.1f}s]",
            flush=True,
        )
        gc.collect()

    df_out = df_base.reset_index()
    df_out.to_csv(REPORTS_DIR / "replogle2022_external_validation_10split.csv", index=False)

    all_methods = ["Constant_Mean", "kNN_30", "AhlmannEltze_Bilinear", "PyTorch_ResMLP", "PyTorch_GEARS_GNN",
                   "S0", "S1", "S2", "S3", "S4", "S5", "S6", "S7"]
    print("\n" + "=" * 108)
    print("REPLOGLE ET AL. (CELL 2022) ZERO-SHOT EXTERNAL VALIDATION SUMMARY (10 Strict Union Splits)")
    print("=" * 108)
    summary = {}
    s0_v5 = df_out["S0_v5"].values
    s7_v5 = df_out["S7_v5"].values
    for m in all_methods:
        v5 = df_out[f"{m}_v5"].values
        dr = df_out[f"{m}_delta"].values
        de = df_out[f"{m}_de50_r"].values
        sg = df_out[f"{m}_de50_sign"].values
        rpe = df_out[f"{m}_v5_RPE1"].values
        k562 = df_out[f"{m}_v5_K562"].values
        _, p_v = (0.0, 1.0) if m == "S7" else ttest_rel(s7_v5, v5)
        d_vs_s0 = v5 - s0_v5
        print(
            f"  {m:22s} | v5: {v5.mean():.4f}±{v5.std():.4f} (RPE1:{rpe.mean():.4f}, K562:{k562.mean():.4f}) "
            f"| DeltaR: {dr.mean():.4f}±{dr.std():.4f} | DE50-R: {de.mean():.4f} | Sign: {sg.mean():.4f} "
            f"| S7 vs m p={p_v:.2e}",
            flush=True,
        )
        summary[m] = {
            "v5_mean": round(float(v5.mean()), 4), "v5_sd": round(float(v5.std()), 4),
            "v5_RPE1_mean": round(float(rpe.mean()), 4), "v5_K562_mean": round(float(k562.mean()), 4),
            "delta_r_mean": round(float(dr.mean()), 4), "delta_r_sd": round(float(dr.std()), 4),
            "de50_r_mean": round(float(de.mean()), 4), "de50_r_sd": round(float(de.std()), 4),
            "de50_sign_mean": round(float(sg.mean()), 4),
            "delta_v5_vs_S0": round(float(d_vs_s0.mean()), 4),
            "paired_p_S7_vs_m": float(p_v),
        }

    with open(REPORTS_DIR / "replogle2022_external_validation_summary.json", "w") as f:
        json.dump({
            "dataset": "Replogle et al. (Cell 2022) K562 GWPS + RPE1 Essential (Figshare+ 20029387)",
            "hyperparameter_transfer": "ZERO-SHOT FROZEN (identical PanPerturbBridgeModel from Orion HCT116/HEK293T)",
            "split_protocol": "10 independent strict simultaneous union-held-out splits (20% per line) + training-only refitted baseline & residual_axis",
            "methods": summary,
        }, f, indent=2)
    print(f"\nMerged and saved replogle2022_external_validation_10split.csv & summary JSON in {time.time()-t_start:.1f}s!", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["pan", "dl"], required=True)
    args = ap.parse_args()
    if args.mode == "pan":
        run_pan_mode()
    else:
        run_dl_mode()
