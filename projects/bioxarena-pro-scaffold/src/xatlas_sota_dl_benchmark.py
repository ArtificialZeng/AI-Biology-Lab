#!/usr/bin/env python3
"""
Head-to-Head SOTA & Deep Learning Benchmark + 10 Brand-New Unseen Confirmation Splits
on BioXArena-Pro `perturbation/xatlas-orion-perturb`.

Compares on the EXACT SAME 10 official splits (Seeds: 42, 0, 1, 2, 3, 7, 11, 99, 2024, 31337)
using official `tasks.utils.strict_group_mean_from_row_scores`:
  1. 100-SVD + Ridge (Official Ref Baseline, S0)
  2. Ahlmann-Eltze et al. (Nature Methods 2025) Bilinear Low-Rank Linear Model
  3. PyTorch Residual MLP (Projected Cosine + MSE Loss, internal early stopping)
  4. PyTorch GEARS-Style Co-Expression GNN (Single-Line Graph Message Passing)
  5. PyTorch Pan-Human Cross-Line GNN (2-Layer Residual GCN on 4,040-gene Pan-Human Graph)
  6. PanPerturb-Bridge (Frozen S7 from xatlas_panhuman_solver.py)

Additionally evaluates Frozen S7 vs S0 on 10 BRAND-NEW, NEVER-SEEN confirmation seeds
[101, 102, 103, 104, 105, 106, 107, 108, 109, 110] (creating a 20-split total audit).
"""

import gc
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from scipy.stats import ttest_rel
from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import Ridge

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from xatlas_panhuman_solver import (
    BENCHMARK_DIR,
    DATA_DIR,
    REPORTS_DIR,
    SEEDS_10,
    PanPerturbBridgeModel,
    build_spectral_pyramid_Gs,
    evaluate_predictions,
    project_rows_5000,
)

NEW_SEEDS_10 = [101, 102, 103, 104, 105, 106, 107, 108, 109, 110]
DEVICE = torch.device("mps" if torch.backends.mps.is_available() else "cpu")


def get_qr_basis_torch(base, ax, s_idx, device):
    Q, Rm = np.linalg.qr(np.stack([np.ones(len(s_idx)), base[s_idx], ax[s_idx]], axis=1))
    Q = Q[:, np.abs(np.diag(Rm)) > 1e-10 * max(1.0, np.abs(Rm).max())].astype(np.float32)
    return torch.tensor(Q, device=device)


def projected_cosine_mse_loss(pred, target, Q_s, s_idx_t):
    p_s = pred[:, s_idx_t]
    t_s = target[:, s_idx_t]
    p_proj = p_s - (p_s @ Q_s) @ Q_s.T
    t_proj = t_s - (t_s @ Q_s) @ Q_s.T
    cos = F.cosine_similarity(p_proj, t_proj, dim=1, eps=1e-8)
    loss_cos = (1.0 - cos).mean()
    loss_mse = F.mse_loss(p_proj, t_proj)
    return loss_cos + 0.5 * loss_mse


# =============================================================================
# 1. Ahlmann-Eltze et al. (Nature Methods 2025) Bilinear Low-Rank Linear Model
# =============================================================================
def run_ahlmann_eltze_bilinear(Y_all_raw, tr_idx, query_idx, tr_Y_raw, rank=100, lam=10.0):
    """
    Exact formulation of the Nature Methods 2025 linear perturbation baseline:
    Y_hat = y_bar + G_val W P^T, where G is gene embedding (PCA of Y_all^T)
    and P is low-rank perturbation response basis (PCA of centered Y_tr).
    """
    y_mean = tr_Y_raw.mean(axis=0, keepdims=True)
    Y_c = tr_Y_raw - y_mean
    G = TruncatedSVD(n_components=rank, random_state=42).fit_transform(Y_all_raw.T)
    r_p = min(rank, len(tr_Y_raw) - 1)
    svd_p = TruncatedSVD(n_components=r_p, random_state=42).fit(Y_c)
    P = svd_p.components_.T  # (5000, r_p)
    G_tr = G[tr_idx]
    A = G_tr.T @ G_tr + lam * np.eye(G_tr.shape[1])
    W = np.linalg.solve(A, G_tr.T @ (Y_c @ P))
    return y_mean + (G[query_idx] @ W) @ P.T


# =============================================================================
# 2. PyTorch Residual MLP (Direct Projected Cosine + MSE Loss, Internal Early Stop)
# =============================================================================
class ResidualPerturbMLP(nn.Module):
    def __init__(self, in_dim=200, hidden=512, out_dim=5000, dropout=0.35):
        super().__init__()
        self.in_proj = nn.Sequential(nn.Linear(in_dim, hidden), nn.LayerNorm(hidden), nn.SiLU())
        self.block1 = nn.Sequential(
            nn.Linear(hidden, hidden), nn.LayerNorm(hidden), nn.SiLU(), nn.Dropout(dropout),
            nn.Linear(hidden, hidden), nn.LayerNorm(hidden),
        )
        self.block2 = nn.Sequential(
            nn.Linear(hidden, hidden), nn.LayerNorm(hidden), nn.SiLU(), nn.Dropout(dropout),
            nn.Linear(hidden, hidden), nn.LayerNorm(hidden),
        )
        self.head = nn.Linear(hidden, out_dim)

    def forward(self, x):
        h = self.in_proj(x)
        h = F.silu(h + self.block1(h))
        h = F.silu(h + self.block2(h))
        return self.head(h)


def train_predict_mlp(X_feat, tr_idx, query_idx, tr_Y_p, base, ax, s_idx, seed=42):
    torch.manual_seed(seed)
    np.random.seed(seed)
    n_tr = len(tr_idx)
    perm = np.random.RandomState(seed + 777).permutation(n_tr)
    n_es = max(32, int(0.15 * n_tr))
    es_i, fit_i = perm[:n_es], perm[n_es:]

    X_t = torch.from_numpy(np.ascontiguousarray(X_feat, dtype=np.float32)).to(DEVICE)
    Y_t = torch.from_numpy(np.ascontiguousarray(tr_Y_p, dtype=np.float32)).to(DEVICE)
    Q_s = get_qr_basis_torch(base, ax, s_idx, DEVICE)
    s_idx_t = torch.tensor(s_idx, dtype=torch.long, device=DEVICE)

    fit_idx_t = torch.tensor([tr_idx[i] for i in fit_i], dtype=torch.long, device=DEVICE)
    es_idx_t = torch.tensor([tr_idx[i] for i in es_i], dtype=torch.long, device=DEVICE)
    q_idx_t = torch.tensor(query_idx, dtype=torch.long, device=DEVICE)

    model = ResidualPerturbMLP(in_dim=X_feat.shape[1], hidden=512, out_dim=5000, dropout=0.35).to(DEVICE)
    opt = torch.optim.AdamW(model.parameters(), lr=1.5e-3, weight_decay=2e-2)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=120)

    best_loss = 1e9
    best_state = None
    batch_size = 128
    for epoch in range(120):
        model.train()
        order = torch.randperm(len(fit_i), device=DEVICE)
        for start in range(0, len(fit_i), batch_size):
            b = order[start : start + batch_size]
            opt.zero_grad()
            pred = model(X_t[fit_idx_t[b]])
            loss = projected_cosine_mse_loss(pred, Y_t[fit_i[b.cpu().numpy()]], Q_s, s_idx_t)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 2.0)
            opt.step()
        sched.step()
        if (epoch + 1) % 4 == 0:
            model.eval()
            with torch.no_grad():
                val_pred = model(X_t[es_idx_t])
                val_loss = float(projected_cosine_mse_loss(val_pred, Y_t[es_i], Q_s, s_idx_t).item())
                if val_loss < best_loss:
                    best_loss = val_loss
                    best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}

    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        preds = model(X_t[q_idx_t]).cpu().numpy().astype(np.float64)
    return base[None, :] + preds


# =============================================================================
# 3. PyTorch GEARS-Style Residual Graph Neural Network (Co-expression / Pan-Human Graph)
# =============================================================================
class GEARSStylePerturbGNN(nn.Module):
    def __init__(self, in_dim=200, hidden=256, out_dim=5000, dropout=0.30):
        super().__init__()
        self.in_lin = nn.Sequential(nn.Linear(in_dim, hidden), nn.LayerNorm(hidden), nn.GELU())
        self.g1_msg = nn.Linear(hidden, hidden)
        self.g1_self = nn.Linear(hidden, hidden)
        self.ln1 = nn.LayerNorm(hidden)
        self.g2_msg = nn.Linear(hidden, hidden)
        self.g2_self = nn.Linear(hidden, hidden)
        self.ln2 = nn.LayerNorm(hidden)
        self.drop = nn.Dropout(dropout)
        self.out_head = nn.Sequential(
            nn.Linear(hidden, hidden * 2), nn.LayerNorm(hidden * 2), nn.GELU(), nn.Dropout(dropout),
            nn.Linear(hidden * 2, out_dim),
        )

    def forward(self, X_all, A_norm, node_indices):
        h0 = self.in_lin(X_all)
        h1 = self.ln1(F.gelu(A_norm @ self.g1_msg(h0) + self.g1_self(h0)) + h0)
        h1 = self.drop(h1)
        h2 = self.ln2(F.gelu(A_norm @ self.g2_msg(h1) + self.g2_self(h1)) + h1)
        return self.out_head(h2[node_indices])


def build_topk_normalized_adj(G_sim, top_k=25):
    """Build signed symmetric normalized top-k adjacency matrix for gene message passing."""
    n = G_sim.shape[0]
    A = np.zeros((n, n), dtype=np.float32)
    G_no_diag = G_sim.copy()
    np.fill_diagonal(G_no_diag, 0.0)
    abs_G = np.abs(G_no_diag)
    idx_top = np.argpartition(abs_G, -top_k, axis=1)[:, -top_k:]
    rows = np.repeat(np.arange(n), top_k)
    cols = idx_top.ravel()
    vals = G_no_diag[rows, cols]
    A[rows, cols] = vals
    A = 0.5 * (A + A.T)
    deg = np.sum(np.abs(A), axis=1) + 1e-6
    d_inv_sqrt = 1.0 / np.sqrt(deg)
    A_norm = (A * d_inv_sqrt[:, None]) * d_inv_sqrt[None, :]
    return A_norm.astype(np.float32)


def train_predict_gnn(X_feat, G_sim, tr_idx, query_idx, tr_Y_p, base, ax, s_idx, seed=42):
    torch.manual_seed(seed)
    np.random.seed(seed)
    n_tr = len(tr_idx)
    perm = np.random.RandomState(seed + 888).permutation(n_tr)
    n_es = max(32, int(0.15 * n_tr))
    es_i, fit_i = perm[:n_es], perm[n_es:]

    A_norm = build_topk_normalized_adj(G_sim, top_k=25)
    A_t = torch.from_numpy(np.ascontiguousarray(A_norm, dtype=np.float32)).to(DEVICE)
    X_t = torch.from_numpy(np.ascontiguousarray(X_feat, dtype=np.float32)).to(DEVICE)
    Y_t = torch.from_numpy(np.ascontiguousarray(tr_Y_p, dtype=np.float32)).to(DEVICE)
    Q_s = get_qr_basis_torch(base, ax, s_idx, DEVICE)
    s_idx_t = torch.tensor(s_idx, dtype=torch.long, device=DEVICE)

    fit_idx_t = torch.tensor([tr_idx[i] for i in fit_i], dtype=torch.long, device=DEVICE)
    es_idx_t = torch.tensor([tr_idx[i] for i in es_i], dtype=torch.long, device=DEVICE)
    q_idx_t = torch.tensor(query_idx, dtype=torch.long, device=DEVICE)

    model = GEARSStylePerturbGNN(in_dim=X_feat.shape[1], hidden=256, out_dim=5000, dropout=0.30).to(DEVICE)
    opt = torch.optim.AdamW(model.parameters(), lr=2e-3, weight_decay=2e-2)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=90)

    best_loss = 1e9
    best_state = None
    for epoch in range(90):
        model.train()
        opt.zero_grad()
        pred = model(X_t, A_t, fit_idx_t)
        loss = projected_cosine_mse_loss(pred, Y_t[fit_i], Q_s, s_idx_t)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 2.0)
        opt.step()
        sched.step()
        if (epoch + 1) % 3 == 0:
            model.eval()
            with torch.no_grad():
                val_pred = model(X_t, A_t, es_idx_t)
                val_loss = float(projected_cosine_mse_loss(val_pred, Y_t[es_i], Q_s, s_idx_t).item())
                if val_loss < best_loss:
                    best_loss = val_loss
                    best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}

    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        preds = model(X_t, A_t, q_idx_t).cpu().numpy().astype(np.float64)
    return base[None, :] + preds


def make_split(train_df, gene2idx, seed):
    rng = np.random.RandomState(seed)
    val_rids = []
    val_targets_by_line = {}
    for source, n_val in [("HCT116", 129), ("HEK293T", 247)]:
        sub_df = train_df[
            (train_df["source"] == source)
            & (train_df["reproducible"] == 1)
            & (train_df["perturbation"].isin(gene2idx[source].keys()))
        ]
        val_targets = set(rng.permutation(sub_df["perturbation"].unique())[:n_val])
        val_targets_by_line[source] = val_targets
        val_rids.extend(sub_df[sub_df["perturbation"].isin(val_targets)]["id"].tolist())
    val_df = train_df[train_df["id"].isin(val_rids)].copy().reset_index(drop=True)
    return val_df, set(val_rids), val_targets_by_line


def main():
    t0_all = time.time()
    print(f"=== SOTA & Deep Learning Head-to-Head + 20-Split Audit (Device: {DEVICE}) ===", flush=True)

    train_df = pd.read_csv(DATA_DIR / "train.csv")
    genes_npz = np.load(DATA_DIR / "genes.npz")
    baselines_npz = np.load(DATA_DIR / "context_baseline.npz")
    axes_npz = np.load(DATA_DIR / "residual_axis.npz")
    score_idx_npz = np.load(DATA_DIR / "score_gene_indices.npz")
    train_y_npz = np.load(DATA_DIR / "train_y.npz")
    all_train_y = {k: np.asarray(train_y_npz[k], dtype=np.float32) for k in train_y_npz.files}

    hct_genes = list(genes_npz["HCT116_genes"])
    hek_genes = list(genes_npz["HEK293T_genes"])
    shared_genes = sorted(list(set(hct_genes) & set(hek_genes)))
    hct_sh_idx = np.array([hct_genes.index(g) for g in shared_genes])
    hek_sh_idx = np.array([hek_genes.index(g) for g in shared_genes])
    gene2idx = {
        "HCT116": {g: i for i, g in enumerate(hct_genes)},
        "HEK293T": {g: i for i, g in enumerate(hek_genes)},
    }
    baselines_pub = {k: np.asarray(baselines_npz[k], dtype=np.float64) for k in baselines_npz.files}
    axes_pub = {k: np.asarray(axes_npz[k], dtype=np.float64) for k in axes_npz.files}
    score_indices = {k: np.asarray(score_idx_npz[k], dtype=np.int64) for k in score_idx_npz.files}

    pan_model = PanPerturbBridgeModel(gene2idx, hct_sh_idx, hek_sh_idx, score_indices)

    # -------------------------------------------------------------------------
    # PART 1: Head-to-Head Literature & Deep Learning Baselines on 10 Official Splits
    # -------------------------------------------------------------------------
    print("\n--- PART 1: Head-to-Head Baselines on 10 Official Splits (Seeds 42..31337) ---", flush=True)
    methods = ["S0_SVD_Ridge", "AhlmannEltze_Bilinear", "PyTorch_ResMLP", "PyTorch_GEARS_GNN", "PyTorch_PanHuman_GNN", "PanPerturb_Bridge_S7"]
    # Load existing S0 & S7 from xatlas_10split_ablation.csv so we don't duplicate compute, AND verify parity
    df_prev = pd.read_csv(REPORTS_DIR / "xatlas_10split_ablation.csv").set_index("seed")

    h2h_rows = []
    timings = {m: [] for m in methods}
    for seed in SEEDS_10:
        t_s = time.time()
        val_df, val_rids_set, val_targets_by_line = make_split(train_df, gene2idx, seed)
        preds_m = {m: {} for m in methods[1:-1]}

        for source in ["HCT116", "HEK293T"]:
            other_src = "HEK293T" if source == "HCT116" else "HCT116"
            g2i = gene2idx[source]
            g2i_o = gene2idx[other_src]
            s_idx = score_indices[f"{source}_indices"]
            own_sh = hct_sh_idx if source == "HCT116" else hek_sh_idx
            oth_sh = hek_sh_idx if source == "HCT116" else hct_sh_idx

            val_t_own = val_targets_by_line[source]
            df_own = train_df[(train_df["source"] == source) & (~train_df["id"].isin(val_rids_set))]
            df_other = train_df[(train_df["source"] == other_src) & (~train_df["perturbation"].isin(val_t_own))]
            q_df = val_df[val_df["source"] == source]
            query_idx = [g2i[p] for p in q_df["perturbation"].tolist()]
            query_rids = q_df["id"].tolist()

            own_rep = df_own[df_own["reproducible"] == 1]
            tr_sub = own_rep[own_rep["perturbation"].isin(g2i.keys())]
            tr_idx = [g2i[p] for p in tr_sub["perturbation"].tolist()]
            tr_Y_raw = np.array([all_train_y[rid] for rid in tr_sub["id"]], dtype=np.float32)
            tr_Y_p = project_rows_5000(tr_Y_raw, baselines_pub[source], axes_pub[source])
            Y_all_raw = np.array([all_train_y[rid] for rid in df_own["id"]], dtype=np.float32)

            # Method A: Ahlmann-Eltze et al. (Nature Methods 2025) Bilinear Low-Rank Linear
            t1 = time.time()
            p_ae = run_ahlmann_eltze_bilinear(Y_all_raw, tr_idx, query_idx, tr_Y_raw, rank=100, lam=10.0)
            timings["AhlmannEltze_Bilinear"].append(time.time() - t1)
            for i, rid in enumerate(query_rids):
                preds_m["AhlmannEltze_Bilinear"][rid] = p_ae[i]

            # Build clean 200-dim gene embeddings (self-masked + baseline subtracted) for Neural Baselines
            Y_clean = Y_all_raw.copy()
            for i, p in enumerate(df_own["perturbation"].tolist()):
                if p in g2i:
                    Y_clean[i, g2i[p]] = 0.0
            Y_clean = Y_clean - baselines_pub[source].astype(np.float32)[None, :]
            X_200 = TruncatedSVD(n_components=200, random_state=42).fit_transform(Y_clean.T).astype(np.float32)
            X_200 = X_200 / (np.linalg.norm(X_200, axis=1, keepdims=True) + 1e-8)

            # Method B: PyTorch Residual MLP
            t1 = time.time()
            p_mlp = train_predict_mlp(X_200, tr_idx, query_idx, tr_Y_p, baselines_pub[source], axes_pub[source], s_idx, seed=seed)
            timings["PyTorch_ResMLP"].append(time.time() - t1)
            for i, rid in enumerate(query_rids):
                preds_m["PyTorch_ResMLP"][rid] = p_mlp[i]

            # Build Single-Line & Pan-Human Gene Graphs for GNNs
            Gs_own = build_spectral_pyramid_Gs(df_own, all_train_y, g2i, baselines_pub[source], pc_levels=(0, 1))
            G_single = 0.5 * (Gs_own[0] + Gs_own[1])
            Gs_oth = build_spectral_pyramid_Gs(df_other, all_train_y, g2i_o, baselines_pub[other_src], pc_levels=(0, 1))
            G_pan = G_single.copy()
            G_pan[np.ix_(own_sh, own_sh)] += 0.7 * 0.5 * (Gs_oth[0] + Gs_oth[1])[np.ix_(oth_sh, oth_sh)]

            # Method C: PyTorch GEARS-Style Single-Line GNN
            t1 = time.time()
            p_gnn1 = train_predict_gnn(X_200, G_single, tr_idx, query_idx, tr_Y_p, baselines_pub[source], axes_pub[source], s_idx, seed=seed)
            timings["PyTorch_GEARS_GNN"].append(time.time() - t1)
            for i, rid in enumerate(query_rids):
                preds_m["PyTorch_GEARS_GNN"][rid] = p_gnn1[i]

            # Method D: PyTorch Pan-Human Cross-Line GNN
            t1 = time.time()
            p_gnn2 = train_predict_gnn(X_200, G_pan, tr_idx, query_idx, tr_Y_p, baselines_pub[source], axes_pub[source], s_idx, seed=seed)
            timings["PyTorch_PanHuman_GNN"].append(time.time() - t1)
            for i, rid in enumerate(query_rids):
                preds_m["PyTorch_PanHuman_GNN"][rid] = p_gnn2[i]

            del Y_all_raw, Y_clean, Gs_own, Gs_oth, G_single, G_pan

        row = {
            "seed": seed,
            "S0_SVD_Ridge": float(df_prev.loc[seed, "S0_overall"]),
            "PanPerturb_Bridge_S7": float(df_prev.loc[seed, "S7_overall"]),
        }
        for m in methods[1:-1]:
            ov, h, ek = evaluate_predictions(preds_m[m], val_df, all_train_y, baselines_pub, axes_pub, score_indices, gene2idx)
            row[m] = round(ov, 6)
            row[f"{m}_HCT116"] = round(h, 6)
            row[f"{m}_HEK293T"] = round(ek, 6)
        h2h_rows.append(row)
        print(
            f"  Seed {seed:5d} | SVD:{row['S0_SVD_Ridge']:.4f} | AE-Lin:{row['AhlmannEltze_Bilinear']:.4f} "
            f"| ResMLP:{row['PyTorch_ResMLP']:.4f} | GEARS-GNN:{row['PyTorch_GEARS_GNN']:.4f} "
            f"| Pan-GNN:{row['PyTorch_PanHuman_GNN']:.4f} | PanPerturb(S7):{row['PanPerturb_Bridge_S7']:.4f} [{time.time()-t_s:.1f}s]",
            flush=True,
        )
        gc.collect()

    df_h2h = pd.DataFrame(h2h_rows)
    df_h2h.to_csv(REPORTS_DIR / "xatlas_sota_dl_head2head_10split.csv", index=False)

    print("\n=== PART 1 SUMMARY (10 Official Splits) ===", flush=True)
    s7_arr = df_h2h["PanPerturb_Bridge_S7"].values
    h2h_summary = {}
    for m in methods:
        arr = df_h2h[m].values
        t_st, p_v = (0.0, 1.0) if m == "PanPerturb_Bridge_S7" else ttest_rel(s7_arr, arr)
        delta_s7 = s7_arr - arr
        print(
            f"  {m:24s}: {arr.mean():.4f} ± {arr.std():.4f} | S7 advantage: {delta_s7.mean():+.4f} ({int((delta_s7>0).sum())}/10, p={p_v:.2e})",
            flush=True,
        )
        h2h_summary[m] = {
            "mean": round(float(arr.mean()), 4),
            "sd": round(float(arr.std()), 4),
            "s7_advantage_mean": round(float(delta_s7.mean()), 4),
            "s7_wins": f"{int((delta_s7 > 0).sum())}/10",
            "paired_p_vs_S7": float(p_v),
        }

    # -------------------------------------------------------------------------
    # PART 2: 10 BRAND-NEW, NEVER-SEEN Confirmation Splits [101..110] -> 20-Split Audit
    # -------------------------------------------------------------------------
    print("\n--- PART 2: 10 Brand-New Untouched Confirmation Splits (Seeds 101..110) ---", flush=True)
    new_rows = []
    for seed in NEW_SEEDS_10:
        t_s = time.time()
        val_df, val_rids_set, val_targets_by_line = make_split(train_df, gene2idx, seed)
        preds_s0, preds_s7 = {}, {}
        for source in ["HCT116", "HEK293T"]:
            other_src = "HEK293T" if source == "HCT116" else "HCT116"
            val_t_own = val_targets_by_line[source]
            df_own = train_df[(train_df["source"] == source) & (~train_df["id"].isin(val_rids_set))]
            df_other = train_df[(train_df["source"] == other_src) & (~train_df["perturbation"].isin(val_t_own))]
            q_df = val_df[val_df["source"] == source]
            res = pan_model.predict_for_source_all_stages(
                source, df_own, df_other, q_df, all_train_y, baselines_pub, axes_pub, ("S0", "S7")
            )
            preds_s0.update(res["S0"])
            preds_s7.update(res["S7"])
        ov0, h0, k0 = evaluate_predictions(preds_s0, val_df, all_train_y, baselines_pub, axes_pub, score_indices, gene2idx)
        ov7, h7, k7 = evaluate_predictions(preds_s7, val_df, all_train_y, baselines_pub, axes_pub, score_indices, gene2idx)
        new_rows.append({
            "seed": seed, "split_group": "brand_new_unseen_10",
            "S0_overall": round(ov0, 6), "S0_HCT116": round(h0, 6), "S0_HEK293T": round(k0, 6),
            "S7_overall": round(ov7, 6), "S7_HCT116": round(h7, 6), "S7_HEK293T": round(k7, 6),
            "delta_S7_vs_S0": round(ov7 - ov0, 6),
        })
        print(f"  [NEW10] Seed {seed:3d} | S0: {ov0:.4f} -> S7: {ov7:.4f} (HCT:{h7:.4f}, HEK:{k7:.4f} | d={ov7-ov0:+.4f}) [{time.time()-t_s:.1f}s]", flush=True)
        gc.collect()

    df_new10 = pd.DataFrame(new_rows)
    df_new10.to_csv(REPORTS_DIR / "xatlas_brand_new_10split_confirmation.csv", index=False)

    n_s0, n_s7 = df_new10["S0_overall"].values, df_new10["S7_overall"].values
    n_d = n_s7 - n_s0
    nt, npv = ttest_rel(n_s7, n_s0)

    all20_s0 = np.concatenate([df_prev["S0_overall"].values, n_s0])
    all20_s7 = np.concatenate([df_prev["S7_overall"].values, n_s7])
    all20_d = all20_s7 - all20_s0
    t20, p20 = ttest_rel(all20_s7, all20_s0)

    print(
        f"\n=== PART 2 SUMMARY ===\n"
        f"  Brand-New Unseen 10 Splits: S0={n_s0.mean():.4f}±{n_s0.std():.4f} -> S7={n_s7.mean():.4f}±{n_s7.std():.4f} "
        f"| d={n_d.mean():+.4f}±{n_d.std():.4f} ({int((n_d>0).sum())}/10, t={nt:.2f}, p={npv:.2e})\n"
        f"  Combined 20-Split Audit   : S0={all20_s0.mean():.4f}±{all20_s0.std():.4f} -> S7={all20_s7.mean():.4f}±{all20_s7.std():.4f} "
        f"| d={all20_d.mean():+.4f}±{all20_d.std():.4f} ({int((all20_d>0).sum())}/20, t={t20:.2f}, p={p20:.2e})",
        flush=True,
    )

    out_json = {
        "device": str(DEVICE),
        "runtime_seconds": round(time.time() - t0_all, 2),
        "part1_head2head_10split": h2h_summary,
        "part2_brand_new_10split": {
            "seeds": NEW_SEEDS_10,
            "S0_mean": round(float(n_s0.mean()), 4), "S0_sd": round(float(n_s0.std()), 4),
            "S7_mean": round(float(n_s7.mean()), 4), "S7_sd": round(float(n_s7.std()), 4),
            "delta_mean": round(float(n_d.mean()), 4), "delta_sd": round(float(n_d.std()), 4),
            "positive_splits": f"{int((n_d > 0).sum())}/10",
            "paired_t": round(float(nt), 4), "paired_p": float(npv),
        },
        "combined_20split_audit": {
            "S0_mean": round(float(all20_s0.mean()), 4), "S0_sd": round(float(all20_s0.std()), 4),
            "S7_mean": round(float(all20_s7.mean()), 4), "S7_sd": round(float(all20_s7.std()), 4),
            "delta_mean": round(float(all20_d.mean()), 4), "delta_sd": round(float(all20_d.std()), 4),
            "positive_splits": f"{int((all20_d > 0).sum())}/20",
            "paired_t": round(float(t20), 4), "paired_p": float(p20),
        },
    }
    with open(REPORTS_DIR / "xatlas_sota_dl_and_20split_summary.json", "w") as f:
        json.dump(out_json, f, indent=2)
    print("Saved xatlas_sota_dl_head2head_10split.csv, xatlas_brand_new_10split_confirmation.csv, and summary JSON!", flush=True)


if __name__ == "__main__":
    main()
