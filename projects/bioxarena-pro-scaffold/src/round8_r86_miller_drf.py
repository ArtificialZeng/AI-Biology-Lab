#!/usr/bin/env python3
"""
Round 8 (R8-6): Miller et al. (Nature Biotechnology 2026; bioRxiv 10.1101/2025.10.20.683304)
Dynamic Range Fraction (DRF) Metric Calibration Framework on Replogle 2022 and Nadig 2025.

Computes exact Miller DRF for every evaluation metric (`MSE`, `delta_r`, `de50_r`, `de50_sign`,
`v5`, `disc_raw`, `disc_v5`):
  - Higher-is-better (perfect = 1.0): DRF = (m(PosCtrl) - m(NegCtrl)) / (1.0 - m(NegCtrl))
  - Lower-is-better  (MSE, perfect = 0.0): DRF = (m(NegCtrl) - m(PosCtrl)) / m(NegCtrl)
where NegCtrl = training reproducible perturbation mean (`Constant_Mean`), and PosCtrl is:
  1. **Replogle 2022 (`K562`)**: True independent experimental replicate screen
     (`K562_essential_normalized_bulk_01.h5ad`, SHA-256 `c1ca6456...`, Day 6 essential screen)
     predicting `K562_gwps_normalized_bulk_01.h5ad` (Day 8 genome-wide screen) on all held-out
     shared K562 targets across `SEEDS_10` (both Raw Duplicate and Miller Interpolated Duplicate).
     Note: RPE1 has no second bulk replicate screen in the published normalized bulk release;
     transparently reported as K562 replicate DRF.
  2. **Nadig 2025 (`HepG2` + `Jurkat`)**: Exact single-cell **split-half** positive controls
     constructed from raw `GSE264667_{hepg2,jurkat}_raw_singlecell_01.h5ad` using locked Addendum B
     per-`gem_group` library-size normalization and control z-scoring (50/50 split of every
     perturbation's single cells into Half A & Half B; Half A predicts Half B, evaluated both as
     Raw Split-Half Duplicate and Miller Interpolated Duplicate blending Half A toward training mean
     by per-gene BH-FDR differential expression weight `w_g = clip(1 - q_g, 0, 1)`).

Also evaluates and ranks `Constant_Mean`, `AE`, `S0`, `S5`, `S7`, `scGPT_FrozenEmb_Ridge`, and
exploratory `OrthoFusion_scGPT_S5` across ALL 7 metrics (including `MSE`) on the exact same
`SEEDS_10` splits, highlighting method rankings on well-calibrated (high-DRF) metrics without
pre-judging the outcome.
"""

import gc
import json
import os
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import torch
from scipy.stats import ttest_ind
from sklearn.linear_model import Ridge
from statsmodels.stats.multitest import multipletests

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from compute_row_bootstrap_ci import row_scores_v5, stratified_row_bootstrap_ci
from prepare_nadig2025_cache import load_nadig_bundle, resolve_nadig_dir
from replogle_external_validation import (
    load_replogle_bundle,
    make_strict_union_split,
    run_ahlmann_eltze_bilinear,
)
from round6_t2_discrim_bootstrap_ci import compute_row_metrics_for_split
from round6_t2_multi_metric_bench import evaluate_5metrics
from round8_r82_and_r83_exploratory import ortho_fusion_predict, pearson_1d, resolve_scgpt_dir
from xatlas_panhuman_solver import REPORTS_DIR, SEEDS_10, PanPerturbBridgeModel


def drf_score(m_pos: float, m_neg: float, lower_is_better: bool = False) -> float:
    if lower_is_better:
        return float((m_neg - m_pos) / max(abs(m_neg), 1e-12))
    return float((m_pos - m_neg) / max(1.0 - m_neg, 1e-12))


def row_mse(preds: dict, val_df: pd.DataFrame, all_y: dict, sidx: dict, g2i_map: dict) -> np.ndarray:
    out = np.empty(len(val_df), dtype=np.float64)
    for i, r in enumerate(val_df.itertuples(index=False)):
        src = r.source
        s = sidx[f"{src}_indices"]
        t = g2i_map[src].get(r.perturbation)
        idx = s[s != t] if t is not None else s
        diff = preds[r.id].astype(np.float64)[idx] - all_y[r.id].astype(np.float64)[idx]
        out[i] = float(np.mean(diff * diff))
    return out


def eval_subset_metrics(
    preds: dict,
    sub_df: pd.DataFrame,
    y_true_map: dict,
    base: np.ndarray,
    axis: np.ndarray,
    s_idx: np.ndarray,
    g2i: dict,
) -> dict:
    """Evaluate all 7 metrics on a single cell line's subset DataFrame."""
    src = sub_df["source"].iloc[0]
    n = len(sub_df)
    rids = sub_df["id"].tolist()
    perts = sub_df["perturbation"].tolist()

    # Build orthonormal Q for v5
    ones = np.ones(len(s_idx), dtype=np.float64)
    M = np.column_stack([ones, base[s_idx], axis[s_idx]])
    Q, R = np.linalg.qr(M, mode="reduced")
    Q = Q[:, np.abs(np.diag(R)) > 1e-10]

    mses, delta_rs, de50_rs, de50_signs, v5s = [], [], [], [], []
    Y_t_raw = np.empty((n, len(s_idx)), dtype=np.float64)
    Y_p_raw = np.empty((n, len(s_idx)), dtype=np.float64)
    Y_t_v5 = np.empty((n, len(s_idx)), dtype=np.float64)
    Y_p_v5 = np.empty((n, len(s_idx)), dtype=np.float64)

    for i, (rid, pert) in enumerate(zip(rids, perts)):
        t = g2i.get(pert)
        yt_full = y_true_map[rid].astype(np.float64)[s_idx]
        yp_full = preds[rid].astype(np.float64)[s_idx]
        dt_full = yt_full - base[s_idx]
        dp_full = yp_full - base[s_idx]
        Y_t_raw[i] = dt_full
        Y_p_raw[i] = dp_full

        vt = yt_full - Q @ (Q.T @ yt_full)
        vp = yp_full - Q @ (Q.T @ yp_full)
        Y_t_v5[i] = vt
        Y_p_v5[i] = vp

        mask = np.ones(len(s_idx), dtype=bool)
        if t is not None and t < len(s_idx):
            mask[t] = False
        dt, dp = dt_full[mask], dp_full[mask]
        mses.append(float(np.mean((yt_full[mask] - yp_full[mask]) ** 2)))
        delta_rs.append(pearson_1d(dt, dp))
        top50 = np.argsort(np.abs(dt))[-50:]
        de50_rs.append(pearson_1d(dt[top50], dp[top50]))
        s_dt, s_dp = np.sign(dt[top50]), np.sign(dp[top50])
        de50_signs.append(float(np.mean(np.where(s_dp == 0, 0.5, (s_dt == s_dp).astype(np.float64)))))
        v5s.append(pearson_1d(vt[mask], vp[mask]))

    # Unprojected and v5 perturbation discrimination (tie-aware mid-rank percentile)
    def disc_pct(A: np.ndarray, B: np.ndarray) -> float:
        A0 = A - A.mean(axis=1, keepdims=True)
        B0 = B - B.mean(axis=1, keepdims=True)
        An = A0 / np.maximum(np.linalg.norm(A0, axis=1, keepdims=True), 1e-12)
        Bn = B0 / np.maximum(np.linalg.norm(B0, axis=1, keepdims=True), 1e-12)
        C = Bn @ An.T  # C[i, j] = corr(pred_i, true_j)
        diag = np.diag(C)[:, None]
        strictly_better = np.sum(C > diag + 1e-12, axis=1)
        tied_others = np.sum(np.abs(C - diag) <= 1e-12, axis=1) - 1
        eff_better = strictly_better + 0.5 * np.maximum(tied_others, 0)
        return float(np.mean(1.0 - eff_better / max(n - 1, 1)))

    return {
        "MSE": float(np.mean(mses)),
        "delta_r": float(np.mean(delta_rs)),
        "de50_r": float(np.mean(de50_rs)),
        "de50_sign": float(np.mean(de50_signs)),
        "v5": float(np.mean(v5s)),
        "disc_raw": disc_pct(Y_t_raw, Y_p_raw),
        "disc_v5": disc_pct(Y_t_v5, Y_p_v5),
    }


def build_or_load_nadig_splithalf(line: str, h5ad_name: str, seed: int = 42) -> tuple[dict, dict, dict]:
    """
    Build exact Addendum-B normalized 50/50 split-half pseudobulk profiles (Half A, Half B) and
    gene-wise BH-FDR DE weights `w_g = clip(1 - q_g, 0, 1)` from raw single-cell Nadig h5ad.
    Returns ({ensg: halfA_5000}, {ensg: halfB_5000}, {ensg: w_interp_5000}).
    """
    nadig_dir = resolve_nadig_dir()
    cache_dir = nadig_dir / "cache"
    sh_dir = nadig_dir / "split_half_cache"
    sh_dir.mkdir(parents=True, exist_ok=True)
    out_npz = sh_dir / f"{line}_splithalf_s{seed}.npz"
    if out_npz.exists():
        print(f"  [Split-Half Cache Hit] {line}", flush=True)
        z = np.load(out_npz, allow_pickle=True)
        ensgs = z["ensgs"].tolist()
        return (
            {g: z["halfA"][i] for i, g in enumerate(ensgs)},
            {g: z["halfB"][i] for i, g in enumerate(ensgs)},
            {g: z["w_interp"][i] for i, g in enumerate(ensgs)},
        )

    t0 = time.time()
    print(f"  [Split-Half Extract] Building Addendum-B 50/50 split-half for Nadig {line}...", flush=True)
    meta_df = pd.read_csv(cache_dir / f"{line}_meta.csv")
    rep_ensgs = set(meta_df[meta_df["reproducible"] == 1]["perturbation"].astype(str).tolist())
    genes_df = pd.read_csv(cache_dir / f"{line}_genes.csv")
    readout_ensg = genes_df["gene_id"].astype(str).tolist()

    fp = nadig_dir / h5ad_name
    with h5py.File(fp, "r") as f:
        n_cells_total, _ = f["X"].shape
        umi = f["obs/UMI_count"][:].astype(np.float64)
        gems = f["obs/gem_group"][:]
        g_cats = [x.decode() if isinstance(x, bytes) else str(x) for x in f["obs/__categories/gene"][:]]
        gid_cats = [x.decode() if isinstance(x, bytes) else str(x) for x in f["obs/__categories/gene_id"][:]]
        g_codes = f["obs/gene"][:]
        gid_codes = f["obs/gene_id"][:]
        var_gids = [x.decode() if isinstance(x, bytes) else str(x) for x in f["var/gene_id"][:]]
        var_pos = {g: i for i, g in enumerate(var_gids)}
        col_idx = np.array([var_pos[g] for g in readout_ensg], dtype=np.int64)

        uniq_gems = np.unique(gems)
        med_umi_by_gem = {int(g): float(np.median(umi[gems == g])) for g in uniq_gems}
        gem_to_idx = {int(g): gi for gi, g in enumerate(uniq_gems)}
        gem_idx_per_cell = np.array([gem_to_idx[int(g)] for g in gems], dtype=np.int32)
        n_gems = len(uniq_gems)
        nt_code = g_cats.index("non-targeting")
        nt_mask = (g_codes == nt_code)

        # Assign every cell of each reproducible perturbation deterministically to Half A (1) or Half B (2)
        rng = np.random.RandomState(seed)
        gid_name_to_code = {name: code for code, name in enumerate(gid_cats)}
        rep_codes = {gid_name_to_code[e]: e for e in sorted(rep_ensgs) if e in gid_name_to_code}
        cells_by_code = {c: [] for c in rep_codes}
        for ci, c in enumerate(gid_codes):
            ci_c = int(c)
            if ci_c in cells_by_code:
                cells_by_code[ci_c].append(ci)

        ensg_list = sorted(rep_codes.values())
        ensg_to_ti = {e: i for i, e in enumerate(ensg_list)}
        n_t = len(ensg_list)
        cell_ti = np.full(n_cells_total, -1, dtype=np.int32)
        cell_half = np.zeros(n_cells_total, dtype=np.int8)  # 1 = Half A, 2 = Half B
        for code, ensg in rep_codes.items():
            arr = np.asarray(cells_by_code[code], dtype=np.int64)
            if len(arr) < 10:
                continue
            perm = rng.permutation(arr)
            mid = len(perm) // 2
            ti = ensg_to_ti[ensg]
            cell_ti[perm[:mid]] = ti
            cell_half[perm[:mid]] = 1
            cell_ti[perm[mid:]] = ti
            cell_half[perm[mid:]] = 2

        # Pass 1: per-gem_group NT mean & std on 5,000 readout genes
        chunk_size = 8192
        X_ds = f["X"]
        n_g = len(col_idx)
        nt_sum = np.zeros((n_gems, n_g), dtype=np.float64)
        nt_sq = np.zeros((n_gems, n_g), dtype=np.float64)
        nt_cnt = np.zeros(n_gems, dtype=np.int64)
        for s in range(0, n_cells_total, chunk_size):
            e = min(s + chunk_size, n_cells_total)
            m = nt_mask[s:e]
            if not np.any(m):
                continue
            Xc = X_ds[s:e, :].astype(np.float64)
            rows = np.where(m)[0]
            X_nt = Xc[rows][:, col_idx]
            gi = gem_idx_per_cell[s:e][rows]
            u = umi[s:e][rows]
            med = np.array([med_umi_by_gem[int(uniq_gems[k])] for k in gi], dtype=np.float64)
            X_nt *= (med / np.maximum(u, 1.0))[:, None]
            np.add.at(nt_sum, gi, X_nt)
            np.add.at(nt_sq, gi, X_nt * X_nt)
            np.add.at(nt_cnt, gi, 1)

        mu_g = nt_sum / np.maximum(nt_cnt[:, None], 1)
        sd_g = np.sqrt(np.maximum((nt_sq - nt_cnt[:, None] * (mu_g ** 2)) / np.maximum(nt_cnt[:, None] - 1, 1), 1e-12))

        # Pass 2: accumulate Half A mean/sq (for Welch t-test vs NT) and Half B mean
        sumA = np.zeros((n_t, n_g), dtype=np.float64)
        sqA = np.zeros((n_t, n_g), dtype=np.float64)
        cntA = np.zeros(n_t, dtype=np.int64)
        sumB = np.zeros((n_t, n_g), dtype=np.float64)
        cntB = np.zeros(n_t, dtype=np.int64)
        nt_z_n = int(nt_cnt.sum())  # By construction z-scored NT has mean=0, var=1

        for s in range(0, n_cells_total, chunk_size):
            e = min(s + chunk_size, n_cells_total)
            h_loc = cell_half[s:e]
            if not np.any(h_loc > 0):
                continue
            Xc = X_ds[s:e, :].astype(np.float64)
            rows = np.where(h_loc > 0)[0]
            X_p = Xc[rows][:, col_idx]
            gi = gem_idx_per_cell[s:e][rows]
            u = umi[s:e][rows]
            med = np.array([med_umi_by_gem[int(uniq_gems[k])] for k in gi], dtype=np.float64)
            X_p *= (med / np.maximum(u, 1.0))[:, None]
            Z_p = (X_p - mu_g[gi]) / sd_g[gi]

            tis = cell_ti[s:e][rows]
            hs = h_loc[rows]
            mA = (hs == 1)
            mB = (hs == 2)
            if np.any(mA):
                np.add.at(sumA, tis[mA], Z_p[mA])
                np.add.at(sqA, tis[mA], Z_p[mA] * Z_p[mA])
                np.add.at(cntA, tis[mA], 1)
            if np.any(mB):
                np.add.at(sumB, tis[mB], Z_p[mB])
                np.add.at(cntB, tis[mB], 1)

    halfA = (sumA / np.maximum(cntA[:, None], 1)).astype(np.float32)
    halfB = (sumB / np.maximum(cntB[:, None], 1)).astype(np.float32)
    varA = np.maximum((sqA - cntA[:, None] * (halfA.astype(np.float64) ** 2)) / np.maximum(cntA[:, None] - 1, 1), 1e-12)

    # Vectorized Welch two-sample t-test (Half A vs NT z-score mean=0, var=1, n=nt_z_n) -> BH FDR -> w_interp
    from scipy.stats import norm
    se = np.sqrt(varA / np.maximum(cntA[:, None], 1) + 1.0 / max(nt_z_n, 1))
    t_stat = np.abs(halfA.astype(np.float64)) / np.maximum(se, 1e-12)
    pvals = 2.0 * norm.sf(t_stat)
    w_interp = np.empty((n_t, n_g), dtype=np.float32)
    for ti in range(n_t):
        _, q, _, _ = multipletests(pvals[ti], method="fdr_bh")
        w_interp[ti] = np.clip(1.0 - q, 0.0, 1.0).astype(np.float32)

    np.savez_compressed(out_npz, ensgs=np.array(ensg_list, dtype=object), halfA=halfA, halfB=halfB, w_interp=w_interp)
    print(f"  [Split-Half Done] {line}: {n_t} reproducible targets in {time.time()-t0:.1f}s", flush=True)
    return (
        {g: halfA[i] for i, g in enumerate(ensg_list)},
        {g: halfB[i] for i, g in enumerate(ensg_list)},
        {g: w_interp[i] for i, g in enumerate(ensg_list)},
    )


def load_replogle_k562_essential_replicate() -> dict:
    """Load K562_essential_normalized_bulk_01.h5ad aligned to cached K562 5,000 readout genes."""
    rdir = Path(os.environ["REPLOGLE_DATA_DIR"])
    cache_dir = rdir / "cache"
    k562_genes = pd.read_csv(cache_dir / "K562_genes.csv")["gene_id"].astype(str).tolist()
    fp = rdir / "K562_essential_normalized_bulk_01.h5ad"
    if not fp.exists():
        raise FileNotFoundError(f"Missing {fp}")
    with h5py.File(fp, "r") as f:
        X = f["X"][:].astype(np.float32)
        var_gids = [x.decode() if isinstance(x, bytes) else str(x) for x in f["var/gene_id"][:]]
        gt = [x.decode() if isinstance(x, bytes) else str(x) for x in f["obs/gene_transcript"][:]]
        nc = f["obs/num_cells_filtered"][:] if "num_cells_filtered" in f["obs"] else np.ones(len(gt))
        ad_c = f["obs/anderson_darling_counts"][:] if "anderson_darling_counts" in f["obs"] else np.zeros(len(gt))
    var_pos = {g: i for i, g in enumerate(var_gids)}
    col_idx = np.array([var_pos[g] for g in k562_genes], dtype=np.int64)
    X_sel = X[:, col_idx]
    # Deduplicate by Ensembl gene_id embedded in gene_transcript (`<num>_<sym>_<pos>_<ENSG>`)
    best = {}
    for i, s in enumerate(gt):
        parts = s.split("_")
        ensg = next((p for p in parts if p.startswith("ENSG")), None)
        if ensg is None:
            continue
        score = (float(ad_c[i]), float(nc[i]))
        if ensg not in best or score > best[ensg][0]:
            best[ensg] = (score, X_sel[i])
    return {e: v[1] for e, v in best.items()}


def run_calibration_and_ranking(dataset: str, emb_mat: np.ndarray, vocab: dict):
    t0 = time.time()
    print(f"\n================ R8-6 Miller DRF Calibration & Model Ranking on {dataset} ================", flush=True)
    pad_id = vocab["<pad>"]
    if dataset == "replogle":
        train_df, all_y, _, sh0, sh1, g2i_map, sidx = load_replogle_bundle()
        k562_rep_map = load_replogle_k562_essential_replicate()
        sh_maps = None
    elif dataset == "nadig":
        train_df, all_y, _, sh0, sh1, g2i_map, sidx = load_nadig_bundle()
        hA_hep, hB_hep, w_hep = build_or_load_nadig_splithalf("HepG2", "GSE264667_hepg2_raw_singlecell_01.h5ad")
        hA_jur, hB_jur, w_jur = build_or_load_nadig_splithalf("Jurkat", "GSE264667_jurkat_raw_singlecell_01.h5ad")
        sh_maps = {
            "HCT116": (hA_hep, hB_hep, w_hep),
            "HEK293T": (hA_jur, hB_jur, w_jur),
        }
        k562_rep_map = None
    else:
        raise ValueError(dataset)

    model = PanPerturbBridgeModel(g2i_map, sh0, sh1, sidx)
    metrics_list = ["MSE", "delta_r", "de50_r", "de50_sign", "v5", "disc_raw", "disc_v5"]
    lower_set = {"MSE"}

    drf_rows = []
    model_rows = []
    ci_rows = []

    for seed in SEEDS_10:
        val_df, tr_df, _, bases, axes = make_strict_union_split(train_df, g2i_map, all_y, sidx, seed)
        df_hct, df_hek = tr_df[tr_df["source"] == "HCT116"], tr_df[tr_df["source"] == "HEK293T"]
        q_hct, q_hek = val_df[val_df["source"] == "HCT116"], val_df[val_df["source"] == "HEK293T"]

        # ---- Part A: DRF Calibration ----
        if dataset == "nadig":
            # Evaluate split-half positive control (Half A -> Half B) and negative control (b_c -> Half B) on both lines
            line_drf_raw, line_drf_interp = [], []
            for src, q_s in [("HCT116", q_hct), ("HEK293T", q_hek)]:
                hA, hB, wI = sh_maps[src]
                sub = q_s[q_s["perturbation"].isin(hA.keys())].copy().reset_index(drop=True)
                y_true_B = {r.id: hB[r.perturbation] for r in sub.itertuples(index=False)}
                p_neg = {r.id: bases[src] for r in sub.itertuples(index=False)}
                p_pos_raw = {r.id: hA[r.perturbation] for r in sub.itertuples(index=False)}
                p_pos_interp = {
                    r.id: bases[src] + wI[r.perturbation] * (hA[r.perturbation] - bases[src])
                    for r in sub.itertuples(index=False)
                }
                m_neg = eval_subset_metrics(p_neg, sub, y_true_B, bases[src], axes[src], sidx[f"{src}_indices"], g2i_map[src])
                m_pr = eval_subset_metrics(p_pos_raw, sub, y_true_B, bases[src], axes[src], sidx[f"{src}_indices"], g2i_map[src])
                m_pi = eval_subset_metrics(p_pos_interp, sub, y_true_B, bases[src], axes[src], sidx[f"{src}_indices"], g2i_map[src])
                line_drf_raw.append({k: drf_score(m_pr[k], m_neg[k], k in lower_set) for k in metrics_list})
                line_drf_interp.append({k: drf_score(m_pi[k], m_neg[k], k in lower_set) for k in metrics_list})
                for k in metrics_list:
                    drf_rows.append({
                        "dataset": dataset, "seed": seed, "line": "HepG2" if src == "HCT116" else "Jurkat",
                        "control_type": "split_half_single_cell", "n_targets": len(sub), "metric": k,
                        "m_neg": round(m_neg[k], 6), "m_pos_raw": round(m_pr[k], 6), "m_pos_interp": round(m_pi[k], 6),
                        "DRF_raw": round(line_drf_raw[-1][k], 6), "DRF_interp": round(line_drf_interp[-1][k], 6),
                    })
        else:
            # Replogle K562 true independent replicate screen (`K562_essential` -> `K562_gwps`)
            sub = q_hek[q_hek["perturbation"].isin(k562_rep_map.keys())].copy().reset_index(drop=True)
            y_true_gwps = {r.id: all_y[r.id] for r in sub.itertuples(index=False)}
            p_neg = {r.id: bases["HEK293T"] for r in sub.itertuples(index=False)}
            p_pos_raw = {r.id: k562_rep_map[r.perturbation] for r in sub.itertuples(index=False)}
            # Interpolated duplicate proxy on bulk z-scores: w_g = clip(|z_rep| / 3.0, 0, 1) (strong |z|>=3 uses replicate)
            p_pos_interp = {
                r.id: bases["HEK293T"] + np.clip(np.abs(k562_rep_map[r.perturbation]) / 3.0, 0.0, 1.0) * (k562_rep_map[r.perturbation] - bases["HEK293T"])
                for r in sub.itertuples(index=False)
            }
            m_neg = eval_subset_metrics(p_neg, sub, y_true_gwps, bases["HEK293T"], axes["HEK293T"], sidx["HEK293T_indices"], g2i_map["HEK293T"])
            m_pr = eval_subset_metrics(p_pos_raw, sub, y_true_gwps, bases["HEK293T"], axes["HEK293T"], sidx["HEK293T_indices"], g2i_map["HEK293T"])
            m_pi = eval_subset_metrics(p_pos_interp, sub, y_true_gwps, bases["HEK293T"], axes["HEK293T"], sidx["HEK293T_indices"], g2i_map["HEK293T"])
            for k in metrics_list:
                drf_rows.append({
                    "dataset": dataset, "seed": seed, "line": "K562_essential_vs_gwps",
                    "control_type": "independent_replicate_screen", "n_targets": len(sub), "metric": k,
                    "m_neg": round(m_neg[k], 6), "m_pos_raw": round(m_pr[k], 6), "m_pos_interp": round(m_pi[k], 6),
                    "DRF_raw": round(drf_score(m_pr[k], m_neg[k], k in lower_set), 6),
                    "DRF_interp": round(drf_score(m_pi[k], m_neg[k], k in lower_set), 6),
                })

        # ---- Part B: Full 7-metric Model Ranking (including MSE) ----
        p_const, p_ae, p_scgpt = {}, {}, {}
        for src, df_s, q_s in [("HCT116", df_hct, q_hct), ("HEK293T", df_hek, q_hek)]:
            g2i = g2i_map[src]
            own_rep = df_s[df_s["reproducible"] == 1]
            mu = np.mean([all_y[rid] for rid in own_rep["id"].tolist()], axis=0).astype(np.float64)
            for rid in q_s["id"].tolist():
                p_const[rid] = mu

            tr_sub = own_rep[own_rep["perturbation"].isin(g2i.keys())]
            tr_idx = [g2i[p] for p in tr_sub["perturbation"].tolist()]
            q_idx = [g2i[p] for p in q_s["perturbation"].tolist()]
            tr_Y_raw = np.array([all_y[rid] for rid in tr_sub["id"]], dtype=np.float32)
            Y_all_raw = np.array([all_y[rid] for rid in df_s["id"]], dtype=np.float32)
            arr_ae = run_ahlmann_eltze_bilinear(Y_all_raw, tr_idx, q_idx, tr_Y_raw, rank=100, lam=10.0)
            for i, rid in enumerate(q_s["id"].tolist()):
                p_ae[rid] = arr_ae[i]

            tr_syms = own_rep["gene_symbol"].astype(str).tolist()
            q_syms = q_s["gene_symbol"].astype(str).tolist()
            X_tr = np.array([emb_mat[vocab.get(s, pad_id)] for s in tr_syms], dtype=np.float32)
            X_q = np.array([emb_mat[vocab.get(s, pad_id)] for s in q_syms], dtype=np.float32)
            Y_tr_all = np.array([all_y[rid] for rid in own_rep["id"]], dtype=np.float32)
            p_arr = Ridge(alpha=10.0).fit(X_tr, Y_tr_all).predict(X_q)
            for i, rid in enumerate(q_s["id"].tolist()):
                p_scgpt[rid] = p_arr[i].astype(np.float64)

        r0 = model.predict_for_source_all_stages("HCT116", df_hct, df_hek, q_hct, all_y, bases, axes, ("S0", "S5", "S7"))
        r1 = model.predict_for_source_all_stages("HEK293T", df_hek, df_hct, q_hek, all_y, bases, axes, ("S0", "S5", "S7"))
        p_s5 = {**r0["S5"], **r1["S5"]}
        preds_map = {
            "Constant_Mean": p_const,
            "AE": p_ae,
            "S0": {**r0["S0"], **r1["S0"]},
            "S5": p_s5,
            "S7": {**r0["S7"], **r1["S7"]},
            "scGPT_FrozenEmb_Ridge": p_scgpt,
            "OrthoFusion_scGPT_S5": ortho_fusion_predict(p_scgpt, p_s5, val_df, bases, axes),
        }
        rm = compute_row_metrics_for_split(preds_map, val_df, all_y, bases, axes, sidx, g2i_map)
        v5_map = {m: row_scores_v5(preds_map[m], val_df, all_y, bases, axes, sidx, g2i_map) for m in preds_map}
        mse_map = {m: row_mse(preds_map[m], val_df, all_y, sidx, g2i_map) for m in preds_map}
        src_arr = val_df["source"].values
        i0, i1 = np.where(src_arr == "HCT116")[0], np.where(src_arr == "HEK293T")[0]

        for m in preds_map:
            ev = evaluate_5metrics(preds_map[m], val_df, all_y, bases, axes, sidx, g2i_map)
            model_rows.append({
                "dataset": dataset, "seed": seed, "method": m,
                "MSE": round(float(0.5 * (mse_map[m][i0].mean() + mse_map[m][i1].mean())), 6),
                **{k: round(v, 6) for k, v in ev.items()},
                "disc_raw": round(float(0.5 * (rm[m]["disc_raw"][i0].mean() + rm[m]["disc_raw"][i1].mean())), 6),
                "disc_v5": round(float(0.5 * (rm[m]["disc_v5"][i0].mean() + rm[m]["disc_v5"][i1].mean())), 6),
            })

        for cname, ma, mb in [
            ("S5_vs_AE", "S5", "AE"),
            ("S5_vs_scGPT_Probe", "S5", "scGPT_FrozenEmb_Ridge"),
            ("scGPT_Probe_vs_AE", "scGPT_FrozenEmb_Ridge", "AE"),
            ("OrthoFusion_vs_scGPT_Probe", "OrthoFusion_scGPT_S5", "scGPT_FrozenEmb_Ridge"),
            ("OrthoFusion_vs_S5", "OrthoFusion_scGPT_S5", "S5"),
        ]:
            ci = stratified_row_bootstrap_ci(mse_map[ma], mse_map[mb], src_arr, n_boot=4000, seed=12345)
            ci_rows.append({
                "dataset": dataset, "seed": seed, "metric": "MSE", "comparison": cname,
                **{k: round(v, 6) if isinstance(v, float) else v for k, v in ci.items()},
            })

    df_drf = pd.DataFrame(drf_rows)
    df_mod = pd.DataFrame(model_rows)
    df_ci = pd.DataFrame(ci_rows)
    df_drf.to_csv(REPORTS_DIR / f"r86_miller_drf_{dataset}.csv", index=False)
    df_mod.to_csv(REPORTS_DIR / f"r86_model_7metrics_{dataset}.csv", index=False)
    df_ci.to_csv(REPORTS_DIR / f"r86_mse_row_ci_{dataset}.csv", index=False)

    # Aggregate DRF summary across seeds (macro-averaged over lines if Nadig)
    summary = {
        "dataset": dataset,
        "reference": "Miller et al., Nature Biotechnology 2026 (bioRxiv 10.1101/2025.10.20.683304)",
        "control_construction": (
            "Exact Addendum-B 50/50 single-cell split-half on both HepG2 and Jurkat (Half A -> Half B)"
            if dataset == "nadig"
            else "True independent experimental replicate screen on K562 (K562_essential Day 6 -> K562_gwps Day 8; RPE1 has no second normalized bulk replicate screen locally)"
        ),
        "drf_by_metric": {},
        "model_means_7metrics": {},
        "mse_comparisons": {},
    }
    for k in metrics_list:
        sub = df_drf[df_drf["metric"] == k]
        # Average across lines per seed first, then mean/std across 10 seeds
        per_seed = sub.groupby("seed")[["m_neg", "m_pos_raw", "m_pos_interp", "DRF_raw", "DRF_interp"]].mean()
        summary["drf_by_metric"][k] = {
            col: {"mean": round(float(per_seed[col].mean()), 4), "std": round(float(per_seed[col].std(ddof=1)), 4)}
            for col in ["m_neg", "m_pos_raw", "m_pos_interp", "DRF_raw", "DRF_interp"]
        }
    for m, g in df_mod.groupby("method"):
        summary["model_means_7metrics"][m] = {
            k: {"mean": round(float(g[k].mean()), 4), "std": round(float(g[k].std(ddof=1)), 4)}
            for k in metrics_list
        }
    for comp, g in df_ci.groupby("comparison"):
        summary["mse_comparisons"][comp] = {
            "mean_delta": round(float(g["delta"].mean()), 6),
            "splits_ci_gt_zero": f"{int(g['ci_gt_zero'].sum())}/{len(g)}",
            "splits_ci_lt_zero": f"{int(g['ci_lt_zero'].sum())}/{len(g)}",
        }

    with open(REPORTS_DIR / f"r86_miller_drf_{dataset}_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"=== Miller DRF Summary ({dataset}) ===", flush=True)
    for k in metrics_list:
        s = summary["drf_by_metric"][k]
        print(f"  {k:10s}: neg={s['m_neg']['mean']:.4f}, pos_raw={s['m_pos_raw']['mean']:.4f}, DRF_raw={s['DRF_raw']['mean']:.4f}, DRF_interp={s['DRF_interp']['mean']:.4f}", flush=True)
    print(f"Done {dataset} in {time.time()-t0:.1f}s", flush=True)


def main():
    wdir = resolve_scgpt_dir()
    vocab = json.load(open(wdir / "scgpt_vocab.json"))
    sd = torch.load(wdir / "scgpt_gh_repo_original_model.bin", map_location="cpu", weights_only=False)
    emb_mat = sd["encoder.embedding.weight"].float().numpy()
    del sd
    run_calibration_and_ranking("replogle", emb_mat, vocab)
    run_calibration_and_ranking("nadig", emb_mat, vocab)


if __name__ == "__main__":
    main()
