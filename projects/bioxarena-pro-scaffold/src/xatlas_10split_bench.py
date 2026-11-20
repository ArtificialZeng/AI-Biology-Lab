#!/usr/bin/env python3
"""
Rigorous 10-Split Benchmarking Harness for xatlas-orion-perturb.
Directly uses official strict_group_mean_from_row_scores.
Computes multi-seed Mean +/- SD across 10 independent random splits.
"""

import os
import sys
import numpy as np
import pandas as pd
from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import Ridge

# Import official scoring logic
sys.path.insert(0, "BioXArena-Pro")
from tasks.utils import strict_group_mean_from_row_scores

DATA_DIR = "BioXArena-Pro/data/perturbation/xatlas-orion-perturb/public"

def _project_out(vec, basis):
    return vec - basis @ (basis.T @ vec)

def score_prediction(true, pred, baseline, axis, scoring_indices, target_idx=None):
    indices = scoring_indices
    if target_idx is not None and target_idx in indices:
        indices = indices[indices != target_idx]
    Q, Rm = np.linalg.qr(np.stack([np.ones(len(indices)), baseline[indices], axis[indices]], axis=1))
    Q = Q[:, np.abs(np.diag(Rm)) > 1e-10 * max(1.0, np.abs(Rm).max())]
    a = _project_out(true[indices] - baseline[indices], Q)
    b = _project_out(pred[indices] - baseline[indices], Q)
    scale = np.linalg.norm(pred[indices]) + np.linalg.norm(baseline[indices])
    if np.linalg.norm(b) <= 1e-4 * scale + 1e-12 or np.linalg.norm(a) == 0:
        return 0.0
    with np.errstate(all="ignore"):
        v = float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))
    return float(np.clip(v, -1.0, 1.0)) if np.isfinite(v) else 0.0

def evaluate_predictions(preds, val_df, train_y, baselines, axes, score_indices, gene2idx):
    row_scores = []
    for _, r in val_df.iterrows():
        rid = r["id"]
        source = r["source"]
        target_gene = r["perturbation"]
        target_idx = gene2idx[source].get(target_gene)
        s = score_prediction(
            train_y[rid], preds[rid], baselines[source], axes[source],
            score_indices[f"{source}_indices"], target_idx
        )
        row_scores.append(s)
    return strict_group_mean_from_row_scores(val_df, row_scores, "source")

def main():
    print("Loading data...")
    train_df = pd.read_csv(f"{DATA_DIR}/train.csv")
    genes_npz = np.load(f"{DATA_DIR}/genes.npz")
    baselines_npz = np.load(f"{DATA_DIR}/context_baseline.npz")
    axes_npz = np.load(f"{DATA_DIR}/residual_axis.npz")
    score_idx_npz = np.load(f"{DATA_DIR}/score_gene_indices.npz")
    train_y_npz = np.load(f"{DATA_DIR}/train_y.npz")
    
    # Pre-cache all train_y vectors in a dict
    print("Caching train_y vectors in memory...")
    all_train_y = {k: np.asarray(train_y_npz[k], dtype=np.float32) for k in train_y_npz.files}

    gene2idx = {
        "HCT116": {g: i for i, g in enumerate(genes_npz["HCT116_genes"])},
        "HEK293T": {g: i for i, g in enumerate(genes_npz["HEK293T_genes"])}
    }
    baselines = {k: np.asarray(baselines_npz[k], dtype=np.float64) for k in baselines_npz.files}
    axes = {k: np.asarray(axes_npz[k], dtype=np.float64) for k in axes_npz.files}
    score_indices = {k: np.asarray(score_idx_npz[k], dtype=np.int64) for k in score_idx_npz.files}

    seeds = [42, 0, 1, 2, 3, 7, 11, 99, 2024, 31337]
    print(f"Benchmarking on {len(seeds)} independent random splits: {seeds}")

    results = {
        "seed": [],
        "constant_mean": [],
        "knn_30": [],
        "svd100_ridge": [],
        "cosine_ridge": []
    }

    for seed in seeds:
        rng = np.random.RandomState(seed)
        val_rids = []
        train_meas_rids = {"HCT116": [], "HEK293T": []}

        for source, n_val in [("HCT116", 129), ("HEK293T", 247)]:
            sub_df = train_df[
                (train_df["source"] == source) &
                (train_df["reproducible"] == 1) &
                (train_df["perturbation"].isin(gene2idx[source].keys()))
            ]
            unique_targets = sub_df["perturbation"].unique()
            val_targets = set(rng.permutation(unique_targets)[:n_val])
            for _, r in sub_df.iterrows():
                if r["perturbation"] in val_targets:
                    val_rids.append(r["id"])
                else:
                    train_meas_rids[source].append(r["id"])

        val_df = train_df[train_df["id"].isin(val_rids)].copy().reset_index(drop=True)

        # 1. Constant Mean
        mean_effects = {s: np.mean([all_train_y[rid] for rid in train_meas_rids[s]], axis=0) for s in ["HCT116", "HEK293T"]}
        preds_const = {r["id"]: mean_effects[r["source"]] for _, r in val_df.iterrows()}
        s_const = evaluate_predictions(preds_const, val_df, all_train_y, baselines, axes, score_indices, gene2idx)

        # Prepare per-line matrices
        preds_knn = {}
        preds_svd = {}
        preds_cos = {}

        for source in ["HCT116", "HEK293T"]:
            line_train_df = train_df[(train_df["source"] == source) & (~train_df["id"].isin(val_rids))]
            tr_rids = line_train_df["id"].tolist()
            Y_all = np.array([all_train_y[rid] for rid in tr_rids], dtype=np.float32) # (N_tr, 5000)

            # Training instances for regression / kNN
            tr_meas_sub = train_df[train_df["id"].isin(train_meas_rids[source])]
            tr_targets = tr_meas_sub["perturbation"].tolist()
            tr_target_idx = [gene2idx[source][t] for t in tr_targets]
            tr_Y = np.array([all_train_y[rid] for rid in tr_meas_sub["id"]], dtype=np.float32)

            val_sub = val_df[val_df["source"] == source]
            val_targets = val_sub["perturbation"].tolist()
            val_target_idx = [gene2idx[source][t] for t in val_targets]

            # 2. 30-kNN co-expression
            C = Y_all.T # (5000, N_tr)
            C_norm = C - C.mean(axis=1, keepdims=True)
            C_norm /= (np.linalg.norm(C_norm, axis=1, keepdims=True) + 1e-12)
            C_tr = C_norm[tr_target_idx]
            C_val = C_norm[val_target_idx]
            corr = C_val @ C_tr.T

            k = 30
            for i, (_, r) in enumerate(val_sub.iterrows()):
                c_i = corr[i]
                top_k = np.argsort(c_i)[-k:]
                w = np.maximum(0.0, c_i[top_k])
                if w.sum() > 0:
                    w /= w.sum()
                    pred = np.sum(tr_Y[top_k] * w[:, None], axis=0)
                else:
                    pred = np.mean(tr_Y[top_k], axis=0)
                preds_knn[r["id"]] = pred.astype(np.float64)

            # 3. 100-SVD + Ridge(alpha=10.0)
            svd = TruncatedSVD(n_components=100, random_state=42)
            Z_genes = svd.fit_transform(Y_all.T)
            X_tr_svd = Z_genes[tr_target_idx]
            ridge_svd = Ridge(alpha=10.0)
            ridge_svd.fit(X_tr_svd, tr_Y)
            X_val_svd = Z_genes[val_target_idx]
            pred_val_svd = ridge_svd.predict(X_val_svd)
            for i, (_, r) in enumerate(val_sub.iterrows()):
                preds_svd[r["id"]] = pred_val_svd[i].astype(np.float64)

            # 4. Cosine Similarity Ridge(alpha=5.0)
            # Feature is cosine similarity to all 5000 genes
            C_unit = Y_all.T / (np.linalg.norm(Y_all.T, axis=1, keepdims=True) + 1e-12)
            # Gene-gene kernel
            K_tr = C_unit[tr_target_idx] @ C_unit.T # (N_meas, 5000)
            K_val = C_unit[val_target_idx] @ C_unit.T # (N_val, 5000)
            ridge_cos = Ridge(alpha=5.0)
            ridge_cos.fit(K_tr, tr_Y)
            pred_val_cos = ridge_cos.predict(K_val)
            for i, (_, r) in enumerate(val_sub.iterrows()):
                preds_cos[r["id"]] = pred_val_cos[i].astype(np.float64)

        s_knn = evaluate_predictions(preds_knn, val_df, all_train_y, baselines, axes, score_indices, gene2idx)
        s_svd = evaluate_predictions(preds_svd, val_df, all_train_y, baselines, axes, score_indices, gene2idx)
        s_cos = evaluate_predictions(preds_cos, val_df, all_train_y, baselines, axes, score_indices, gene2idx)

        results["seed"].append(seed)
        results["constant_mean"].append(s_const)
        results["knn_30"].append(s_knn)
        results["svd100_ridge"].append(s_svd)
        results["cosine_ridge"].append(s_cos)

        print(f"Seed {seed:5d} | Const: {s_const:.4f} | kNN: {s_knn:.4f} | SVD+Ridge: {s_svd:.4f} | Cosine Ridge: {s_cos:.4f} | Delta: {s_cos - s_svd:+.4f}")

    df_res = pd.DataFrame(results)
    print("\n" + "="*70)
    print("10-SPLIT BENCHMARK SUMMARY:")
    for col in ["constant_mean", "knn_30", "svd100_ridge", "cosine_ridge"]:
        vals = df_res[col].values
        print(f"  {col:15s}: Mean = {vals.mean():.4f} +/- {vals.std():.4f} (Min: {vals.min():.4f}, Max: {vals.max():.4f})")
    
    delta = df_res["cosine_ridge"] - df_res["svd100_ridge"]
    print(f"\n  Cosine Ridge Delta over SVD+Ridge: Mean = {delta.mean():+.4f} +/- {delta.std():.4f} (Positive: {(delta > 0).sum()}/10)")
    print("="*70)

if __name__ == "__main__":
    main()
