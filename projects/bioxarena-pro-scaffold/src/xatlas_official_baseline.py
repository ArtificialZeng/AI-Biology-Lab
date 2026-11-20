#!/usr/bin/env python3
"""
TASK-007: Official Baseline Reproduction for xatlas-orion-perturb.
Strictly adheres to official scoring logic and validation split definition.
"""

import os
import sys
import numpy as np
import pandas as pd
from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import Ridge

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

def evaluate_model_on_val(val_preds, val_rows, train_y, baselines, axes, score_indices, gene2idx):
    scores_by_line = {"HCT116": [], "HEK293T": []}
    for r_id, pred in val_preds.items():
        row = val_rows[r_id]
        source = row["source"]
        target_gene = row["perturbation"]
        target_idx = gene2idx[source].get(target_gene)
        
        baseline = baselines[source]
        axis = axes[source]
        s_idx = score_indices[f"{source}_indices"]
        true = np.asarray(train_y[r_id], dtype=np.float64)
        
        s = score_prediction(true, pred, baseline, axis, s_idx, target_idx)
        scores_by_line[source].append(s)
        
    res = {}
    line_means = []
    for line, s_list in scores_by_line.items():
        m = float(np.mean(s_list))
        se = float(np.std(s_list) / np.sqrt(len(s_list)))
        res[line] = (m, se)
        line_means.append(max(0.0, m))
    res["overall"] = float(np.mean(line_means))
    return res

def main():
    print("Loading xatlas-orion-perturb data...")
    train_df = pd.read_csv(f"{DATA_DIR}/train.csv")
    genes_npz = np.load(f"{DATA_DIR}/genes.npz")
    baselines_npz = np.load(f"{DATA_DIR}/context_baseline.npz")
    axes_npz = np.load(f"{DATA_DIR}/residual_axis.npz")
    score_idx_npz = np.load(f"{DATA_DIR}/score_gene_indices.npz")
    train_y = np.load(f"{DATA_DIR}/train_y.npz")
    
    gene2idx = {
        "HCT116": {g: i for i, g in enumerate(genes_npz["HCT116_genes"])},
        "HEK293T": {g: i for i, g in enumerate(genes_npz["HEK293T_genes"])}
    }
    
    baselines = {k: np.asarray(baselines_npz[k], dtype=np.float64) for k in baselines_npz.files}
    axes = {k: np.asarray(axes_npz[k], dtype=np.float64) for k in axes_npz.files}
    score_indices = {k: np.asarray(score_idx_npz[k], dtype=np.int64) for k in score_idx_npz.files}
    
    # Validation split construction:
    # 20% of reproducible measured-target rows (129 HCT116, 247 HEK293T)
    val_rows = {}
    train_meas_rows = {"HCT116": [], "HEK293T": []}
    
    rng = np.random.RandomState(42)
    for source, n_val in [("HCT116", 129), ("HEK293T", 247)]:
        source_genes = set(gene2idx[source].keys())
        sub_df = train_df[(train_df["source"] == source) & (train_df["reproducible"] == 1) & (train_df["perturbation"].isin(source_genes))].copy()
        
        # Unique target genes (1 row per target gene)
        unique_targets = sub_df["perturbation"].unique()
        perm = rng.permutation(unique_targets)
        val_targets = set(perm[:n_val])
        
        for _, r in sub_df.iterrows():
            rid = r["id"]
            if r["perturbation"] in val_targets:
                val_rows[rid] = r.to_dict()
            else:
                train_meas_rows[source].append(rid)
                
    print(f"Validation set: {len(val_rows)} rows (HCT116: {sum(1 for r in val_rows.values() if r['source'] == 'HCT116')}, HEK293T: {sum(1 for r in val_rows.values() if r['source'] == 'HEK293T')})")
    
    # 1. Baseline: Constant Baseline Vector
    pred_base = {rid: baselines[row["source"]] for rid, row in val_rows.items()}
    res_base = evaluate_model_on_val(pred_base, val_rows, train_y, baselines, axes, score_indices, gene2idx)
    print(f"Baseline (Zeros/Baseline): Overall = {res_base['overall']:.4f} (Ref: 0.000)")
    
    # 2. Baseline: Constant Training Mean of Reproducible Measured Targets
    mean_effects = {}
    for source in ["HCT116", "HEK293T"]:
        rids = train_meas_rows[source]
        vecs = [train_y[r] for r in rids]
        mean_effects[source] = np.mean(vecs, axis=0).astype(np.float64)
        
    pred_const = {rid: mean_effects[row["source"]] for rid, row in val_rows.items()}
    res_const = evaluate_model_on_val(pred_const, val_rows, train_y, baselines, axes, score_indices, gene2idx)
    print(f"Constant Mean: Overall = {res_const['overall']:.4f} | HCT116 = {res_const['HCT116'][0]:.4f}±{res_const['HCT116'][1]:.4f} | HEK293T = {res_const['HEK293T'][0]:.4f}±{res_const['HEK293T'][1]:.4f} (Ref: 0.009)")

    # 3. Baseline: 100-component SVD + Ridge(alpha=10)
    print("\nTraining 100-SVD + Ridge baseline...")
    pred_ridge = {}
    for source in ["HCT116", "HEK293T"]:
        # Build training effect matrix for the line (all training rows of the line that are NOT in validation)
        line_train_df = train_df[(train_df["source"] == source) & (~train_df["id"].isin(val_rows.keys()))]
        tr_rids = line_train_df["id"].tolist()
        print(f"  [{source}] Loading {len(tr_rids)} training knock-downs for gene embedding...")
        
        # Load effect matrix in batches
        Y_all = np.array([train_y[rid] for rid in tr_rids], dtype=np.float32)  # (N_tr, 5000)
        
        # SVD of the column matrix (5000 genes x N_tr)
        print(f"  [{source}] Computing 100-component SVD of gene columns...")
        svd = TruncatedSVD(n_components=100, random_state=42)
        Z_genes = svd.fit_transform(Y_all.T)  # (5000, 100)
        
        # Training instances: reproducible measured targets
        tr_meas_rids = [rid for rid in train_meas_rows[source]]
        tr_target_indices = [gene2idx[source][train_df.loc[train_df["id"] == rid, "perturbation"].values[0]] for rid in tr_meas_rids]
        
        X_tr = Z_genes[tr_target_indices]  # (N_meas, 100)
        Y_tr = np.array([train_y[rid] for rid in tr_meas_rids], dtype=np.float32)  # (N_meas, 5000)
        
        ridge = Ridge(alpha=10.0)
        ridge.fit(X_tr, Y_tr)
        
        # Predict on validation rows
        for rid, row in val_rows.items():
            if row["source"] == source:
                t_idx = gene2idx[source][row["perturbation"]]
                x_val = Z_genes[t_idx:t_idx+1]
                pred_ridge[rid] = ridge.predict(x_val)[0].astype(np.float64)
                
    res_ridge = evaluate_model_on_val(pred_ridge, val_rows, train_y, baselines, axes, score_indices, gene2idx)
    print(f"\n100-SVD + Ridge (alpha=10): Overall = {res_ridge['overall']:.4f} | HCT116 = {res_ridge['HCT116'][0]:.4f}±{res_ridge['HCT116'][1]:.4f} | HEK293T = {res_ridge['HEK293T'][0]:.4f}±{res_ridge['HEK293T'][1]:.4f} (Ref: 0.051 [HCT: 0.060±0.013, HEK: 0.042±0.008])")

if __name__ == "__main__":
    main()
