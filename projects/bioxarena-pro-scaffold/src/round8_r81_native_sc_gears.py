#!/usr/bin/env python3
"""
Round 8 (R8-1): Native Single-Cell GEARS (`cell-gears==0.1.2`, Roohani et al., Nature Biotechnology 2023)
on Nadig et al. (2025, GSE264667 HepG2 + Jurkat) with Official Default Pipeline (`skip_calc_de=False`).

Pre-declared protocol (reports/ROUND8_PREREGISTRATION.md Section 1):
  - Exact same `make_strict_union_split` on seeds `[42, 0, 1]`.
  - Per-`gem_group` library-size normalization to `gem_group` median UMI (matching Addendum B) + `log1p`.
  - Up to `max_cells_per_pert` (default 64; may be lowered to >=32 only if single seed >3h on MPS)
    real single cells per perturbation + up to 5,000 `non-targeting` control cells.
  - Native `PertData.new_data_process(..., skip_calc_de=False)` using TRAIN + CTRL cells ONLY
    (ZERO test cells ever enter AnnData, DE ranking, or co-expression graph).
  - Official `GEARS` default hyperparameters (`hidden_size=64, epochs=20, lr=1e-3, weight_decay=5e-4,
    batch_size=32`), 90/10 train/internal-val split.
  - Converts predicted post-perturbation `log1p` expression to parameter-free control-standardized
    pseudobulk z-score `(pred_log - mu_ctrl_log) / sd_ctrl_log` (plus unstandardized `pred_log - mu_ctrl_log`
    sensitivity check), evaluated on identical metrics and 4,000-resample row-level paired bootstrap 95% CIs
    vs S5 and AE.
"""

import argparse
import gc
import json
import os
import pickle
import shutil
import sys
import tempfile
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from compute_row_bootstrap_ci import row_scores_v5, stratified_row_bootstrap_ci
from gears import GEARS, PertData
from prepare_nadig2025_cache import load_nadig_bundle, resolve_nadig_dir, stratified_nt_sample
from replogle_external_validation import make_strict_union_split, run_ahlmann_eltze_bilinear
from round6_t2_discrim_bootstrap_ci import compute_row_metrics_for_split
from round6_t2_multi_metric_bench import evaluate_5metrics
from xatlas_panhuman_solver import REPORTS_DIR, PanPerturbBridgeModel


def make_unique_gene_names(names: list[str]) -> list[str]:
    seen = {}
    out = []
    for g in names:
        if g not in seen:
            seen[g] = 0
            out.append(g)
        else:
            seen[g] += 1
            out.append(f"{g}__dup{seen[g]}")
    return out


def extract_or_load_singlecell_cache(
    line: str,
    h5ad_name: str,
    max_cells_per_pert: int = 64,
    max_ctrl_cells: int = 5000,
    sample_seed: int = 42,
) -> tuple[sp.csr_matrix, np.ndarray, np.ndarray, list[str], np.ndarray, np.ndarray]:
    """
    Sequentially scan raw Nadig h5ad once and cache subsampled per-gem_group libnorm+log1p single cells
    for the exact deduplicated `gene_transcript` perturbations in `{line}_meta.csv` and exact 5,000 readout genes.
    Returns (X_csr_log1p, cell_gt_labels, cell_is_ctrl, gene_names_5000, mu_ctrl_log, sd_ctrl_log).
    """
    nadig_dir = resolve_nadig_dir()
    cache_dir = nadig_dir / "cache"
    sc_dir = nadig_dir / "sc_gears_cache"
    sc_dir.mkdir(parents=True, exist_ok=True)
    tag = f"{line}_cap{max_cells_per_pert}_ctrl{max_ctrl_cells}_s{sample_seed}"
    npz_path = sc_dir / f"{tag}_X.npz"
    meta_path = sc_dir / f"{tag}_meta.json"
    labels_path = sc_dir / f"{tag}_labels.npy"
    stats_path = sc_dir / f"{tag}_ctrl_stats.npz"

    genes_df = pd.read_csv(cache_dir / f"{line}_genes.csv")
    readout_ensg = genes_df["gene_id"].astype(str).tolist()
    readout_names = genes_df["gene_name"].astype(str).tolist()

    if npz_path.exists() and meta_path.exists() and labels_path.exists() and stats_path.exists():
        print(f"  [SC Cache Hit] Loading {tag}...", flush=True)
        X_csr = sp.load_npz(npz_path)
        labels = np.load(labels_path, allow_pickle=True)
        st = np.load(stats_path)
        is_ctrl = (labels == "ctrl")
        return X_csr, labels, is_ctrl, readout_names, st["mu_ctrl_log"], st["sd_ctrl_log"]

    t0 = time.time()
    fp = nadig_dir / h5ad_name
    if not fp.exists():
        raise FileNotFoundError(f"Missing Nadig raw h5ad: {fp}")
    meta_df = pd.read_csv(cache_dir / f"{line}_meta.csv")
    valid_ensgs = set(meta_df["perturbation"].astype(str).tolist())

    print(f"  [SC Extract] Scanning {h5ad_name} (cap={max_cells_per_pert}/pert, ctrl={max_ctrl_cells})...", flush=True)
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
        med_umi_per_cell = np.array([med_umi_by_gem[int(g)] for g in gems], dtype=np.float64)

        # Select control cells (gem_group-stratified)
        nt_code = g_cats.index("non-targeting")
        nt_all = np.where(g_codes == nt_code)[0]
        ctrl_idx = stratified_nt_sample(nt_all, gems, max_ctrl_cells, seed=sample_seed)

        # Select up to max_cells_per_pert per Ensembl gene_id in meta_df
        rng = np.random.RandomState(sample_seed)
        keep_rows = [ctrl_idx]
        keep_labels = [np.full(len(ctrl_idx), "ctrl", dtype=object)]

        gid_name_to_code = {name: code for code, name in enumerate(gid_cats)}
        target_codes = {gid_name_to_code[ensg]: ensg for ensg in valid_ensgs if ensg in gid_name_to_code}
        cells_by_code = {c: [] for c in target_codes}
        for ci, c in enumerate(gid_codes):
            ci_code = int(c)
            if ci_code in cells_by_code:
                cells_by_code[ci_code].append(ci)

        for code in sorted(cells_by_code.keys()):
            arr = np.asarray(cells_by_code[code], dtype=np.int64)
            if len(arr) == 0:
                continue
            if len(arr) > max_cells_per_pert:
                arr = np.sort(rng.choice(arr, size=max_cells_per_pert, replace=False))
            keep_rows.append(arr)
            keep_labels.append(np.full(len(arr), target_codes[code], dtype=object))

        all_keep_idx = np.concatenate(keep_rows)
        all_keep_lbl = np.concatenate(keep_labels)
        order = np.argsort(all_keep_idx)
        all_keep_idx = all_keep_idx[order]
        all_keep_lbl = all_keep_lbl[order]

        # Build fast lookup boolean mask for sequential chunk reading
        keep_mask = np.zeros(n_cells_total, dtype=bool)
        keep_mask[all_keep_idx] = True

        chunk_size = 8192
        X_ds = f["X"]
        csr_chunks = []
        ctrl_sum = np.zeros(len(col_idx), dtype=np.float64)
        ctrl_sq = np.zeros(len(col_idx), dtype=np.float64)
        ctrl_n = 0

        # Cursor into sorted all_keep_lbl
        cursor = 0
        for s in range(0, n_cells_total, chunk_size):
            e = min(s + chunk_size, n_cells_total)
            m_loc = keep_mask[s:e]
            n_loc = int(np.sum(m_loc))
            if n_loc == 0:
                continue
            Xc = X_ds[s:e, :].astype(np.float32)
            rows_loc = np.where(m_loc)[0]
            X_sub = Xc[rows_loc][:, col_idx]
            u_sub = umi[s:e][rows_loc].astype(np.float32)
            med_sub = med_umi_per_cell[s:e][rows_loc].astype(np.float32)
            scale = med_sub / np.maximum(u_sub, 1.0)
            X_sub *= scale[:, None]
            np.log1p(X_sub, out=X_sub)

            lbl_chunk = all_keep_lbl[cursor : cursor + n_loc]
            cursor += n_loc
            c_mask = (lbl_chunk == "ctrl")
            if np.any(c_mask):
                Xc_ctrl = X_sub[c_mask].astype(np.float64)
                ctrl_sum += Xc_ctrl.sum(axis=0)
                ctrl_sq += (Xc_ctrl * Xc_ctrl).sum(axis=0)
                ctrl_n += int(np.sum(c_mask))

            csr_chunks.append(sp.csr_matrix(X_sub))

    X_csr = sp.vstack(csr_chunks, format="csr")
    mu_ctrl_log = ctrl_sum / max(ctrl_n, 1)
    var_ctrl_log = np.maximum((ctrl_sq - ctrl_n * (mu_ctrl_log ** 2)) / max(ctrl_n - 1, 1), 1e-12)
    sd_ctrl_log = np.sqrt(var_ctrl_log)

    sp.save_npz(npz_path, X_csr)
    np.save(labels_path, all_keep_lbl)
    np.savez(stats_path, mu_ctrl_log=mu_ctrl_log, sd_ctrl_log=sd_ctrl_log, ctrl_n=ctrl_n)
    with open(meta_path, "w") as fm:
        json.dump({
            "line": line, "h5ad": h5ad_name, "n_cells": int(X_csr.shape[0]),
            "n_genes": int(X_csr.shape[1]), "ctrl_n": int(ctrl_n),
            "max_cells_per_pert": max_cells_per_pert, "extract_sec": round(time.time() - t0, 1),
        }, fm, indent=2)

    print(f"  [SC Extract Done] {line}: shape={X_csr.shape}, nnz={X_csr.nnz:,}, ctrl={ctrl_n}, elapsed={time.time()-t0:.1f}s", flush=True)
    return X_csr, all_keep_lbl, (all_keep_lbl == "ctrl"), readout_names, mu_ctrl_log, sd_ctrl_log


def run_native_gears_one_line(
    src_name: str,
    line_name: str,
    df_own: pd.DataFrame,
    q_df: pd.DataFrame,
    all_y: dict,
    X_sc: sp.csr_matrix,
    sc_labels: np.ndarray,
    gene_names: list[str],
    mu_ctrl_log: np.ndarray,
    sd_ctrl_log: np.ndarray,
    seed: int,
    gears_cache: Path,
    device: str = "mps",
    epochs: int = 20,
    max_cells_per_pert: int = 64,
) -> tuple[dict, dict, dict]:
    """
    Train native single-cell official `cell-gears==0.1.2` (`skip_calc_de=False`) on one cell line & split.
    Returns (preds_zscore, preds_rawlog_delta, meta_stats).
    """
    t0 = time.time()
    torch.manual_seed(seed)
    np.random.seed(seed)
    rng = np.random.RandomState(seed)

    own_rep = df_own[df_own["reproducible"] == 1].copy().reset_index(drop=True)
    mu_train_z = np.mean([all_y[rid] for rid in own_rep["id"].tolist()], axis=0).astype(np.float64)

    # Load official GO graph gene list
    go_pkl = gears_cache / "gene2go_all.pkl"
    if not go_pkl.exists():
        from gears.utils import dataverse_download
        dataverse_download("https://dataverse.harvard.edu/api/access/datafile/6153417", str(go_pkl))
    with open(go_pkl, "rb") as f:
        gene2go = pickle.load(f)
    go_genes = set(gene2go.keys())

    # Index single-cell rows by label (Ensembl ID or "ctrl")
    rows_by_lbl = {}
    for idx_i, lbl in enumerate(sc_labels):
        rows_by_lbl.setdefault(lbl, []).append(idx_i)

    # Filter reproducible training rows present in GO graph AND having >= 10 cells in sc cache
    tr_seen_sym = set()
    tr_items = []
    for _, r in own_rep.iterrows():
        sym = str(r["gene_symbol"])
        ensg = str(r["perturbation"])
        if sym in go_genes and sym not in tr_seen_sym and ensg in rows_by_lbl and len(rows_by_lbl[ensg]) >= 10:
            tr_seen_sym.add(sym)
            tr_items.append((r["id"], sym, ensg))

    q_in_go, q_out_go = [], []
    for _, r in q_df.iterrows():
        sym = str(r["gene_symbol"])
        if sym in go_genes and sym not in tr_seen_sym:
            q_in_go.append((r["id"], sym))
        else:
            q_out_go.append((r["id"], sym))

    # Build AnnData with STRICTLY TRAIN + CTRL real single cells (ZERO test cells!)
    ctrl_indices = np.asarray(rows_by_lbl["ctrl"], dtype=np.int64)
    sel_indices = [ctrl_indices]
    conds = ["ctrl"] * len(ctrl_indices)
    tr_conds = []

    for rid, sym, gt in tr_items:
        arr = np.asarray(rows_by_lbl[gt], dtype=np.int64)
        if len(arr) > max_cells_per_pert:
            arr = np.sort(rng.choice(arr, size=max_cells_per_pert, replace=False))
        c = f"{sym}+ctrl"
        tr_conds.append(c)
        sel_indices.append(arr)
        conds.extend([c] * len(arr))

    # Add 1 minimal dummy test condition sampled from ctrl cells ONLY if PertData requires non-empty test split
    # (uses ctrl cells ONLY — NEVER true test cells; test_loader is popped before training anyway!)
    dummy_test_conds = []
    if len(q_in_go) > 0:
        d_sym = q_in_go[0][1]
        d_c = f"{d_sym}+ctrl"
        dummy_test_conds.append(d_c)
        dummy_idx = rng.choice(ctrl_indices, size=min(10, len(ctrl_indices)), replace=False)
        sel_indices.append(dummy_idx)
        conds.extend([d_c] * len(dummy_idx))

    all_sel = np.concatenate(sel_indices)
    X_sub = X_sc[all_sel].tocsr().astype(np.float32)

    obs_df = pd.DataFrame({
        "condition": conds,
        "cell_type": [line_name] * len(conds),
    })
    obs_df.index = [f"cell_{i}" for i in range(len(conds))]
    uniq_gnames = make_unique_gene_names(gene_names)
    var_df = pd.DataFrame({"gene_name": uniq_gnames}, index=uniq_gnames)
    adata = ad.AnnData(X=X_sub, obs=obs_df, var=var_df)

    # 90% train / 10% internal val split of training conditions
    perm_tr = rng.permutation(len(tr_conds))
    n_ival = max(10, int(round(0.10 * len(tr_conds))))
    ival_conds = [tr_conds[i] for i in perm_tr[:n_ival]]
    itr_conds = [tr_conds[i] for i in perm_tr[n_ival:]] + ["ctrl"]
    split_dict = {"train": itr_conds, "val": ival_conds, "test": dummy_test_conds}

    print(
        f"  [Native GEARS {line_name} seed={seed}] cells={len(conds):,} (ctrl={len(ctrl_indices):,}, "
        f"tr_perts={len(tr_items)}, cap={max_cells_per_pert}), running native calc_de...",
        flush=True,
    )

    tmp_root = Path(tempfile.mkdtemp(prefix=f"native_gears_{line_name}_s{seed}_"))
    try:
        for fname in ["gene2go_all.pkl", "essential_all_data_pert_genes.pkl"]:
            src_f = gears_cache / fname
            if src_f.exists():
                os.symlink(src_f, tmp_root / fname)

        pd_obj = PertData(str(tmp_root), default_pert_graph=True)
        # NATIVE GEARS preprocessing WITH skip_calc_de=False (official scanpy rank_genes_groups!)
        t_prep = time.time()
        pd_obj.new_data_process(dataset_name="split_ds", adata=adata, skip_calc_de=False)
        print(f"    [Native calc_de + PertData done in {time.time()-t_prep:.1f}s]", flush=True)

        for fname in ["gene2go_all.pkl", "essential_all_data_pert_genes.pkl"]:
            f_tmp = tmp_root / fname
            f_dst = gears_cache / fname
            if f_tmp.exists() and not f_tmp.is_symlink() and not f_dst.exists():
                shutil.copy2(f_tmp, f_dst)

        split_pkl = tmp_root / "custom_split.pkl"
        with open(split_pkl, "wb") as f:
            pickle.dump(split_dict, f)
        pd_obj.prepare_split(split="custom", seed=seed, split_dict_path=str(split_pkl))
        pd_obj.get_dataloader(batch_size=32, test_batch_size=64)
        pd_obj.dataloader.pop("test_loader", None)

        gears_model = GEARS(pd_obj, device=device, weight_bias_track=False)
        gears_model.model_initialize(hidden_size=64)
        t_tr = time.time()
        gears_model.train(epochs=epochs, lr=1e-3, weight_decay=5e-4)
        train_sec = time.time() - t_tr
        assert gears_model.best_model is not None, "GEARS best_model was not set"
        print(f"    [Native GEARS train done in {train_sec:.1f}s ({train_sec/epochs:.1f}s/ep)]", flush=True)

        preds_z = {}
        preds_logd = {}
        uniq_q_syms = sorted(set(sym for _, sym in q_in_go))
        if len(uniq_q_syms) > 0:
            raw_pred = gears_model.predict([[s] for s in uniq_q_syms])
            for rid, sym in q_in_go:
                x_log = np.asarray(raw_pred[sym], dtype=np.float64)
                # Primary pre-registered control-standardized pseudobulk z-score
                preds_z[rid] = (x_log - mu_ctrl_log) / np.maximum(sd_ctrl_log, 1e-6)
                # Sensitivity check: unstandardized log1p delta
                preds_logd[rid] = x_log - mu_ctrl_log
        for rid, _ in q_out_go:
            preds_z[rid] = mu_train_z.copy()
            preds_logd[rid] = np.zeros_like(mu_train_z)

        meta = {
            "line": line_name,
            "n_cells_train_adata": int(len(conds)),
            "n_ctrl": int(len(ctrl_indices)),
            "max_cells_per_pert": max_cells_per_pert,
            "n_train_in_go": len(tr_items),
            "n_train_total": len(own_rep),
            "n_test_in_go": len(q_in_go),
            "n_test_total": len(q_df),
            "test_go_coverage": round(len(q_in_go) / max(len(q_df), 1), 4),
            "in_go_rids": [rid for rid, _ in q_in_go],
            "train_sec": round(train_sec, 1),
            "total_sec": round(time.time() - t0, 1),
        }
    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)
        if device == "mps":
            torch.mps.empty_cache()
        gc.collect()

    return preds_z, preds_logd, meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, nargs="+", default=[42, 0, 1])
    ap.add_argument("--max-cells", type=int, default=64)
    ap.add_argument("--max-ctrl", type=int, default=5000)
    ap.add_argument("--epochs", type=int, default=20)
    ap.add_argument("--device", default="mps" if torch.backends.mps.is_available() else "cpu")
    ap.add_argument("--resume", action="store_true", help="Resume from existing metrics/row_ci CSVs by skipping completed seeds")
    args = ap.parse_args()

    t_all = time.time()
    print(
        f"=== R8-1 Native Single-Cell GEARS on Nadig 2025 (seeds={args.seeds}, "
        f"cap={args.max_cells}, ctrl={args.max_ctrl}, device={args.device}, resume={args.resume}) ===",
        flush=True,
    )
    train_df, all_y, _, sh0, sh1, g2i_map, sidx = load_nadig_bundle()
    model = PanPerturbBridgeModel(g2i_map, sh0, sh1, sidx)
    gears_cache = Path(os.environ.get("GEARS_CACHE_DIR", str(Path.home() / ".cache/r7_gears_data")))
    gears_cache.mkdir(parents=True, exist_ok=True)

    # Extract/load native single-cell matrices once
    X_hep, lbl_hep, _, gnames_hep, mu_hep, sd_hep = extract_or_load_singlecell_cache(
        "HepG2", "GSE264667_hepg2_raw_singlecell_01.h5ad",
        max_cells_per_pert=args.max_cells, max_ctrl_cells=args.max_ctrl, sample_seed=42,
    )
    X_jur, lbl_jur, _, gnames_jur, mu_jur, sd_jur = extract_or_load_singlecell_cache(
        "Jurkat", "GSE264667_jurkat_raw_singlecell_01.h5ad",
        max_cells_per_pert=args.max_cells, max_ctrl_cells=args.max_ctrl, sample_seed=42,
    )

    rows = []
    ci_rows = []
    cov_list = []
    done_seeds = set()
    metrics_csv = REPORTS_DIR / "r81_native_sc_gears_nadig_metrics.csv"
    ci_csv = REPORTS_DIR / "r81_native_sc_gears_nadig_row_ci.csv"
    if args.resume and metrics_csv.exists() and ci_csv.exists():
        prev_m = pd.read_csv(metrics_csv)
        prev_c = pd.read_csv(ci_csv)
        if not prev_m.empty and not prev_c.empty:
            prev_m = prev_m[
                (prev_m["max_cells_per_pert"] == args.max_cells)
                & (prev_m["max_ctrl_cells"] == args.max_ctrl)
            ]
            done_seeds = set(prev_m["seed"].astype(int).tolist()) & set(prev_c["seed"].astype(int).tolist())
            prev_m = prev_m[prev_m["seed"].astype(int).isin(done_seeds)]
            prev_c = prev_c[prev_c["seed"].astype(int).isin(done_seeds)]
            rows = prev_m.to_dict(orient="records")
            ci_rows = prev_c.to_dict(orient="records")
            print(f"  [Resume] Loaded {len(rows)} metric rows & {len(ci_rows)} CI rows for completed seeds={sorted(done_seeds)}", flush=True)

    for seed in args.seeds:
        if seed in done_seeds:
            print(f"  [Resume] Skipping already completed seed={seed}", flush=True)
            continue
        t_s = time.time()
        val_df, df_tr, _, bases, axes = make_strict_union_split(train_df, g2i_map, all_y, sidx, seed)
        df_hct = df_tr[df_tr["source"] == "HCT116"]
        df_hek = df_tr[df_tr["source"] == "HEK293T"]
        q_hct = val_df[val_df["source"] == "HCT116"]
        q_hek = val_df[val_df["source"] == "HEK293T"]

        # Native single-cell GEARS per cell line
        pg0_z, pg0_ld, m0 = run_native_gears_one_line(
            "HCT116", "HepG2", df_hct, q_hct, all_y, X_hep, lbl_hep, gnames_hep,
            mu_hep, sd_hep, seed, gears_cache, device=args.device, epochs=args.epochs,
            max_cells_per_pert=args.max_cells,
        )
        pg1_z, pg1_ld, m1 = run_native_gears_one_line(
            "HEK293T", "Jurkat", df_hek, q_hek, all_y, X_jur, lbl_jur, gnames_jur,
            mu_jur, sd_jur, seed, gears_cache, device=args.device, epochs=args.epochs,
            max_cells_per_pert=args.max_cells,
        )
        p_gears_z = {**pg0_z, **pg1_z}
        p_gears_ld = {**pg0_ld, **pg1_ld}

        # Reference S5 and AE on exact same split
        p_ae = {}
        for src, df_s, q_s in [("HCT116", df_hct, q_hct), ("HEK293T", df_hek, q_hek)]:
            g2i = g2i_map[src]
            own_rep = df_s[df_s["reproducible"] == 1]
            tr_sub = own_rep[own_rep["perturbation"].isin(g2i.keys())]
            tr_idx = [g2i[p] for p in tr_sub["perturbation"].tolist()]
            q_idx = [g2i[p] for p in q_s["perturbation"].tolist()]
            tr_Y_raw = np.array([all_y[rid] for rid in tr_sub["id"]], dtype=np.float32)
            Y_all_raw = np.array([all_y[rid] for rid in df_s["id"]], dtype=np.float32)
            arr_ae = run_ahlmann_eltze_bilinear(Y_all_raw, tr_idx, q_idx, tr_Y_raw, rank=100, lam=10.0)
            for i, rid in enumerate(q_s["id"].tolist()):
                p_ae[rid] = arr_ae[i]

        r0 = model.predict_for_source_all_stages("HCT116", df_hct, df_hek, q_hct, all_y, bases, axes, ("S5",))
        r1 = model.predict_for_source_all_stages("HEK293T", df_hek, df_hct, q_hek, all_y, bases, axes, ("S5",))
        p_s5 = {**r0["S5"], **r1["S5"]}

        preds_map = {
            "Native_SC_GEARS": p_gears_z,
            "Native_SC_GEARS_rawlog": p_gears_ld,
            "S5": p_s5,
            "AE": p_ae,
        }
        rm = compute_row_metrics_for_split(preds_map, val_df, all_y, bases, axes, sidx, g2i_map)
        v5_map = {m: row_scores_v5(preds_map[m], val_df, all_y, bases, axes, sidx, g2i_map) for m in preds_map}
        src_arr = val_df["source"].values
        i0, i1 = np.where(src_arr == "HCT116")[0], np.where(src_arr == "HEK293T")[0]

        in_go_set = set(m0["in_go_rids"]) | set(m1["in_go_rids"])
        go_mask = val_df["id"].isin(in_go_set).values
        i0_go = np.where((src_arr == "HCT116") & go_mask)[0]
        i1_go = np.where((src_arr == "HEK293T") & go_mask)[0]

        for mname in ["Native_SC_GEARS", "Native_SC_GEARS_rawlog"]:
            ev = evaluate_5metrics(preds_map[mname], val_df, all_y, bases, axes, sidx, g2i_map)
            rows.append({
                "dataset": "nadig", "seed": seed, "method": mname,
                "max_cells_per_pert": args.max_cells, "max_ctrl_cells": args.max_ctrl,
                **{k: round(v, 6) for k, v in ev.items()},
                "disc_raw": round(float(0.5 * (rm[mname]["disc_raw"][i0].mean() + rm[mname]["disc_raw"][i1].mean())), 6),
                "disc_v5": round(float(0.5 * (rm[mname]["disc_v5"][i0].mean() + rm[mname]["disc_v5"][i1].mean())), 6),
                "v5_go_only": round(float(0.5 * (v5_map[mname][i0_go].mean() + v5_map[mname][i1_go].mean())), 6),
                "disc_raw_go_only": round(float(0.5 * (rm[mname]["disc_raw"][i0_go].mean() + rm[mname]["disc_raw"][i1_go].mean())), 6),
                "go_cov_hepg2": m0["test_go_coverage"],
                "go_cov_jurkat": m1["test_go_coverage"],
                "seed_elapsed_sec": round(time.time() - t_s, 1),
            })

        for cname, ma, mb in [
            ("S5_vs_Native_SC_GEARS", "S5", "Native_SC_GEARS"),
            ("Native_SC_GEARS_vs_AE", "Native_SC_GEARS", "AE"),
            ("S5_vs_Native_SC_GEARS_rawlog", "S5", "Native_SC_GEARS_rawlog"),
            ("Native_SC_GEARS_rawlog_vs_AE", "Native_SC_GEARS_rawlog", "AE"),
        ]:
            for met, a_v, b_v in [
                ("v5", v5_map[ma], v5_map[mb]),
                ("delta_r", rm[ma]["delta_r"], rm[mb]["delta_r"]),
                ("disc_raw", rm[ma]["disc_raw"], rm[mb]["disc_raw"]),
                ("disc_v5", rm[ma]["disc_v5"], rm[mb]["disc_v5"]),
            ]:
                ci = stratified_row_bootstrap_ci(a_v, b_v, src_arr, n_boot=4000, seed=12345)
                ci_rows.append({
                    "dataset": "nadig", "seed": seed, "metric": met, "comparison": cname,
                    **{k: round(v, 6) if isinstance(v, float) else v for k, v in ci.items()},
                })

        cov_list.append({
            "seed": seed,
            "HepG2": {k: v for k, v in m0.items() if k != "in_go_rids"},
            "Jurkat": {k: v for k, v in m1.items() if k != "in_go_rids"},
        })

        # Incremental checkpoint save after every seed
        pd.DataFrame(rows).to_csv(REPORTS_DIR / "r81_native_sc_gears_nadig_metrics.csv", index=False)
        pd.DataFrame(ci_rows).to_csv(REPORTS_DIR / "r81_native_sc_gears_nadig_row_ci.csv", index=False)
        r_primary = [r for r in rows if r["method"] == "Native_SC_GEARS"][-1]
        print(
            f"  [R8-1 seed={seed} DONE in {(time.time()-t_s)/60:.1f}m] "
            f"Native_SC_GEARS: v5={r_primary['v5']:.4f}, delta_r={r_primary['delta_r']:.4f}, "
            f"disc_raw={r_primary['disc_raw']:.4f}, de50_sign={r_primary['de50_sign']:.4f}",
            flush=True,
        )

    # Final summary JSON + pre-declared adjudication
    df_m = pd.DataFrame(rows)
    df_c = pd.DataFrame(ci_rows)
    summary = {
        "dataset": "nadig",
        "package": "cell-gears==0.1.2 (native single-cell, skip_calc_de=False)",
        "max_cells_per_pert": args.max_cells,
        "max_ctrl_cells": args.max_ctrl,
        "epochs": args.epochs,
        "seeds": args.seeds,
        "device": args.device,
        "total_elapsed_sec": round(time.time() - t_all, 1),
        "coverage_and_timing": cov_list,
        "metrics": {},
        "comparisons": {},
    }
    for mname, g in df_m.groupby("method"):
        summary["metrics"][mname] = {
            col: {"mean": round(float(g[col].mean()), 6), "std": round(float(g[col].std(ddof=1)) if len(g) > 1 else 0.0, 6)}
            for col in ["v5", "delta_r", "de50_r", "de50_sign", "disc_raw", "disc_v5", "v5_go_only", "disc_raw_go_only"]
        }
    for (met, comp), g in df_c.groupby(["metric", "comparison"]):
        summary["comparisons"].setdefault(met, {})[comp] = {
            "mean_delta": round(float(g["delta"].mean()), 6),
            "mean_ci95_lo": round(float(g["ci95_lo"].mean()), 6),
            "mean_ci95_hi": round(float(g["ci95_hi"].mean()), 6),
            "splits_ci_gt_zero": f"{int(g['ci_gt_zero'].sum())}/{len(g)}",
            "splits_ci_lt_zero": f"{int(g['ci_lt_zero'].sum())}/{len(g)}",
        }
    with open(REPORTS_DIR / "r81_native_sc_gears_nadig_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"=== R8-1 Native Single-Cell GEARS Complete in {(time.time()-t_all)/60:.1f}m ===", flush=True)


if __name__ == "__main__":
    main()
