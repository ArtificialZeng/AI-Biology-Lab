#!/usr/bin/env python3
"""
TASK-007: Official Baseline Reproduction for xatlas-orion-perturb.
Strictly adheres to official scoring logic and validation split definition.
Imports official strict_group_mean_from_row_scores from tasks.utils.
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

def evaluate_model_on_val(val_preds, val_df, train_y, baselines, axes, score_indices, gene2idx):
    row_scores = []
    for _, r in val_df.iterrows():
        rid = r["id"]
        source = r["source"]
        target_gene = r["perturbation"]
        target_idx = gene2idx[source].get(target_gene)
        s = score_prediction(
            train_y[rid], val_preds[rid], baselines[source], axes[source],
            score_indices[f"{source}_indices"], target_idx
        )
        row_scores.append(s)
    return strict_group_mean_from_row_scores(val_df, row_scores, "source")

def main():
    print("Loading xatlas-orion-perturb data...")
    train_df = pd.read_csv(f"{DATA_DIR}/train.csv")
    genes_npz = np.load(f"{DATA_DIR}/genes.npz")
    baselines_npz = np.load(f"{DATA_DIR}/context_baseline.npz")
    axes_npz = np.load(f"{DATA_DIR}/residual_axis.npz")
    score_idx_npz = np.load(f"{DATA_DIR}/score_gene_indices.npz")
    train_y_npz = np.load(f"{DATA_DIR}/train_y.npz")
    
    gene2idx = {
        "HCT116": {g: i for i, g in enumerate(genes_npz["HCT116_genes"])},
        "HEK293T": {g: i for i, g in enumerate(genes_npz["HEK293T_genes"])}
    }
    
    baselines = {k: np.asarray(baselines_npz[k], dtype=np.float64) for k in baselines_npz.files}
    axes = {k: np.asarray(axes_npz[k], dtype=np.float64) for k in axes_npz.files}
    score_indices = {k: np.asarray(score_idx_npz[k], dtype=np.int64) for k in score_idx_npz.files}
    all_train_y = {k: np.asarray(train_y_npz[k], dtype=np.float32) for k in train_y_npz.files}
    
    # Validation split construction (seed 42 reference split: 129 HCT116, 247 HEK293T)
    val_rids = []
    train_meas_rids = {"HCT116": [], "HEK293T": []}
    
    rng = np.random.RandomState(42)
    for source, n_val in [("HCT116", 129), ("HEK293T", 247)]:
        source_genes = set(gene2idx[source].keys())
        sub_df = train_df[
            (train_df["source"] == source) &
            (train_df["reproducible"] == 1) &
            (train_df["perturbation"].isin(source_genes))
        ].copy()
        
        unique_targets = sub_df["perturbation"].unique()
        perm = rng.permutation(unique_targets)
        val_targets = set(perm[:n_val])
        
        for _, r in sub_df.iterrows():
            if r["perturbation"] in val_targets:
                val_rids.append(r["id"])
            else:
                train_meas_rids[source].append(r["id"])
                
    val_df = train_df[train_df["id"].isin(val_rids)].copy().reset_index(drop=True)
    print(f"Validation set: {len(val_df)} rows (HCT116: {sum(val_df['source'] == 'HCT116')}, HEK293T: {sum(val_df['source'] == 'HEK293T')})")
    
    # 1. Baseline: Constant Baseline Vector
    pred_base = {r["id"]: baselines[r["source"]] for _, r in val_df.iterrows()}
    s_base = evaluate_model_on_val(pred_base, val_df, all_train_y, baselines, axes, score_indices, gene2idx)
    print(f"Baseline (Zeros/Baseline): Overall = {s_base:.4f} (Ref: 0.000)")
    
    # 2. Baseline: Constant Training Mean of Reproducible Measured Targets
    mean_effects = {s: np.mean([all_train_y[rid] for rid in train_meas_rids[s]], axis=0) for s in ["HCT116", "HEK293T"]}
    pred_const = {r["id"]: mean_effects[r["source"]] for _, r in val_df.iterrows()}
    s_const = evaluate_model_on_val(pred_const, val_df, all_train_y, baselines, axes, score_indices, gene2idx)
    print(f"Constant Mean: Overall = {s_const:.4f} (Ref: 0.009)")

    # 3. Baseline: 30-kNN co-expression
    pred_knn = {}
    for source in ["HCT116", "HEK293T"]:
        line_train_df = train_df[(train_df["source"] == source) & (~train_df["id"].isin(val_rids))]
        tr_rids = line_train_df["id"].tolist()
        Y_all = np.array([all_train_y[rid] for rid in tr_rids], dtype=np.float32)
        
        C = Y_all.T
        C_norm = C - C.mean(axis=1, keepdims=True)
        C_norm /= (np.linalg.norm(C_norm, axis=1, keepdims=True) + 1e-12)
        
        tr_meas_sub = train_df[train_df["id"].isin(train_meas_rids[source])]
        tr_targets = tr_meas_sub["perturbation"].tolist()
        tr_target_idx = [gene2idx[source][t] for t in tr_targets]
        tr_Y = np.array([all_train_y[rid] for rid in tr_meas_sub["id"]], dtype=np.float32)
        
        val_sub = val_df[val_df["source"] == source]
        val_targets = val_sub["perturbation"].tolist()
        val_target_idx = [gene2idx[source][t] for t in val_targets]
        
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
            pred_knn[r["id"]] = pred.astype(np.float64)
            
    s_knn = evaluate_model_on_val(pred_knn, val_df, all_train_y, baselines, axes, score_indices, gene2idx)
    print(f"30-kNN Co-expression: Overall = {s_knn:.4f} (Ref: 0.047)")

    # 4. Baseline: 100-component SVD + Ridge(alpha=10)
    pred_ridge = {}
    for source in ["HCT116", "HEK293T"]:
        line_train_df = train_df[(train_df["source"] == source) & (~train_df["id"].isin(val_rids))]
        tr_rids = line_train_df["id"].tolist()
        Y_all = np.array([all_train_y[rid] for rid in tr_rids], dtype=np.float32)
        
        svd = TruncatedSVD(n_components=100, random_state=42)
        Z_genes = svd.fit_transform(Y_all.T)
        
        tr_meas_sub = train_df[train_df["id"].isin(train_meas_rids[source])]
        tr_targets = tr_meas_sub["perturbation"].tolist()
        tr_target_idx = [gene2idx[source][t] for t in tr_targets]
        tr_Y = np.array([all_train_y[rid] for rid in tr_meas_sub["id"]], dtype=np.float32)
        
        X_tr = Z_genes[tr_target_idx]
        val_sub = val_df[val_df["source"] == source]
        val_targets = val_sub["perturbation"].tolist()
        val_target_idx = [gene2idx[source][t] for t in val_targets]
        X_val = Z_genes[val_target_idx]
        
        ridge = Ridge(alpha=10.0)
        ridge.fit(X_tr, tr_Y)
        p_val = ridge.predict(X_val)
        for i, (_, r) in enumerate(val_sub.iterrows()):
            pred_ridge[r["id"]] = p_val[i].astype(np.float64)
                
    s_ridge = evaluate_model_on_val(pred_ridge, val_df, all_train_y, baselines, axes, score_indices, gene2idx)
    print(f"100-SVD + Ridge (alpha=10): Overall = {s_ridge:.4f} (Ref: 0.051; 10-split mean is 0.0613)")

if __name__ == "__main__":
    main()
