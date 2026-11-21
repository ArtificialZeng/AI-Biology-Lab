#!/usr/bin/env python3
"""
Round 7 (R7-2): Official GEARS (`cell-gears==0.1.2`, Roohani et al., Nature Biotechnology 2023)
Benchmark on Replogle 2022 and Nadig 2025 using Official Default Hyperparameters
(`hidden_size=64, num_go_gnn_layers=1, num_gene_gnn_layers=1, decoder_hidden_size=16,
  epochs=20, lr=1e-3, weight_decay=5e-4, batch_size=32`).

Strictly evaluated on the EXACT SAME 10 strict simultaneous union-held-out splits (`SEEDS_10`),
exact same 5 metrics + `disc_raw` + `disc_v5` + Systema `centroid_train_r`, and 4,000-resample
cell-line-stratified row-level paired bootstrap 95% CIs against S5, S7, and AE.

Notes on Official GEARS Integration:
  - Runs in an isolated venv with `cell-gears==0.1.2` (e.g., `$HOME/.cache/r7_gears_venv`).
  - Official GEARS is single-cell-line by design; trained independently per cell line on
    each split's training set (90% train / 10% internal val for GEARS early stopping;
    benchmark held-out perturbations strictly in `test`, zero test-label leakage).
  - Official GEARS requires target gene symbols to exist in its GO similarity graph
    (`gene2go_all.pkl`); any held-out target absent from the GO graph receives the
    training reproducible mean fallback (coverage % and GO-covered-only scores both reported).
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
import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch

torch.set_num_threads(6)

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from gears import GEARS, PertData
from gears.utils import get_similarity_network

from compute_row_bootstrap_ci import row_scores_v5, stratified_row_bootstrap_ci
from replogle_external_validation import (
    load_replogle_bundle,
    make_strict_union_split,
    run_ahlmann_eltze_bilinear,
)
from round6_t2_discrim_bootstrap_ci import compute_row_metrics_for_split
from round6_t2_multi_metric_bench import evaluate_5metrics
from round7_r74_systema_centroid_eval import compute_centroid_rows
from xatlas_panhuman_solver import (
    REPORTS_DIR,
    SEEDS_10,
    PanPerturbBridgeModel,
)


def resolve_gears_cache_dir() -> Path:
    d = Path(os.environ.get("GEARS_CACHE_DIR", str(Path.home() / ".cache/r7_gears_data")))
    d.mkdir(parents=True, exist_ok=True)
    return d


def load_bundle(dataset: str):
    if dataset == "replogle":
        train_df, all_y, _, sh0, sh1, g2i_map, sidx = load_replogle_bundle()
        cache_dir = Path(os.environ["REPLOGLE_DATA_DIR"]) / "cache"
        g_line0 = pd.read_csv(cache_dir / "RPE1_genes.csv")["gene_name"].astype(str).tolist()
        g_line1 = pd.read_csv(cache_dir / "K562_genes.csv")["gene_name"].astype(str).tolist()
        return train_df, all_y, sh0, sh1, g2i_map, sidx, {"HCT116": g_line0, "HEK293T": g_line1}
    elif dataset == "nadig":
        from prepare_nadig2025_cache import load_nadig_bundle, resolve_nadig_dir
        train_df, all_y, _, sh0, sh1, g2i_map, sidx = load_nadig_bundle()
        cache_dir = resolve_nadig_dir() / "cache"
        g_line0 = pd.read_csv(cache_dir / "HepG2_genes.csv")["gene_name"].astype(str).tolist()
        g_line1 = pd.read_csv(cache_dir / "Jurkat_genes.csv")["gene_name"].astype(str).tolist()
        return train_df, all_y, sh0, sh1, g2i_map, sidx, {"HCT116": g_line0, "HEK293T": g_line1}
    else:
        raise ValueError(f"Unsupported dataset: {dataset}")


def make_unique_gene_names(names: list[str]) -> list[str]:
    seen = {}
    out = []
    for nm in names:
        if nm not in seen:
            seen[nm] = 0
            out.append(nm)
        else:
            seen[nm] += 1
            out.append(f"{nm}_dup{seen[nm]}")
    return out


def run_official_gears_one_line(
    src_name: str,
    df_own: pd.DataFrame,
    q_df: pd.DataFrame,
    all_y: dict,
    base: np.ndarray,
    gene_names: list[str],
    seed: int,
    gears_cache: Path,
    device: str = "cpu",
    epochs: int = 20,
) -> tuple[dict, dict]:
    """
    Train official `cell-gears` (default hyperparams) on one cell line & split,
    returning ({rid: pred_5000}, coverage_stats).
    """
    torch.manual_seed(seed)
    np.random.seed(seed)

    own_rep = df_own[df_own["reproducible"] == 1].copy().reset_index(drop=True)
    mu_train = np.mean([all_y[rid] for rid in own_rep["id"].tolist()], axis=0).astype(np.float64)

    # Ensure PertData GO graph metadata is downloaded once into gears_cache
    pert_data_probe = PertData(str(gears_cache), default_pert_graph=True)
    go_pkl = gears_cache / "gene2go_all.pkl"
    if not go_pkl.exists():
        # Trigger official download via helper
        from gears.utils import dataverse_download
        dataverse_download("https://dataverse.harvard.edu/api/access/datafile/6153417", str(go_pkl))
    with open(go_pkl, "rb") as f:
        gene2go = pickle.load(f)
    go_genes = set(gene2go.keys())

    # Filter training reproducible rows to symbols in official GO graph (deduplicate symbols)
    tr_seen_sym = set()
    tr_rows = []
    for _, r in own_rep.iterrows():
        sym = str(r["gene_symbol"])
        if sym in go_genes and sym not in tr_seen_sym:
            tr_seen_sym.add(sym)
            tr_rows.append((r["id"], sym))

    q_in_go, q_out_go = [], []
    q_seen_sym = set()
    for _, r in q_df.iterrows():
        sym = str(r["gene_symbol"])
        if sym in go_genes and sym not in tr_seen_sym and sym not in q_seen_sym:
            q_seen_sym.add(sym)
            q_in_go.append((r["id"], sym))
        elif sym in go_genes and sym not in tr_seen_sym:
            q_in_go.append((r["id"], sym))
        else:
            q_out_go.append((r["id"], sym))

    # Build temporary AnnData for this split (reproducible train + ctrl replicates + test placeholders)
    # Note: test profiles are ONLY used in test_loader evaluation inside GEARS; GEARS never trains on 'test'.
    # For zero leakage even in DE/coexpression graph construction, we compute DE & coexpression on train+ctrl only!
    n_ctrl = 10
    X_rows = []
    conds = []
    # Control cells = split training baseline `base` (small tiny jitter 1e-6 to avoid zero std in rank_genes_groups)
    rng = np.random.RandomState(seed)
    for _ in range(n_ctrl):
        X_rows.append((base + rng.normal(0.0, 1e-5, size=len(base))).astype(np.float32))
        conds.append("ctrl")

    # 1 pseudobulk profile per training perturbation (skip_calc_de=True avoids rank_genes_groups)
    tr_conds = []
    for rid, sym in tr_rows:
        y = all_y[rid].astype(np.float32)
        c = f"{sym}+ctrl"
        tr_conds.append(c)
        X_rows.append(y)
        conds.append(c)

    # Add dummy/placeholder test rows with baseline expression (NEVER true y_test!) so zero test leakage can occur
    # anywhere in PertData preprocessing or co-expression graph construction!
    test_conds = []
    for rid, sym in q_in_go:
        c = f"{sym}+ctrl"
        if c not in test_conds:
            test_conds.append(c)
            X_rows.append(base.astype(np.float32))
            conds.append(c)

    X_mat = sp.csr_matrix(np.vstack(X_rows).astype(np.float32))
    obs_df = pd.DataFrame({
        "condition": conds,
        "cell_type": [src_name] * len(conds),
        "condition_name": [f"{src_name}_{c}_1+1" for c in conds],
    })
    obs_df.index = [f"cell_{i}" for i in range(len(conds))]
    uniq_gnames = make_unique_gene_names(gene_names)
    var_df = pd.DataFrame({"gene_name": uniq_gnames}, index=uniq_gnames)
    adata = ad.AnnData(X=X_mat, obs=obs_df, var=var_df)

    # Pre-populate GEARS uns fields cleanly from TRAIN+CTRL ONLY (zero test leakage, valid for real-valued z-scores)
    all_idx = np.arange(len(uniq_gnames), dtype=np.int64)
    non_zeros_dict = {}
    top20_de_dict = {}
    rgg_dict = {}
    for rid, sym in tr_rows:
        cname = f"{src_name}_{sym}+ctrl_1+1"
        non_zeros_dict[cname] = all_idx
        dt = np.abs(all_y[rid].astype(np.float64) - base)
        top20 = np.argsort(dt)[-20:][::-1]
        top20_de_dict[cname] = top20
        rgg_dict[cname] = np.array([uniq_gnames[j] for j in top20])
    for rid, sym in q_in_go:
        cname = f"{src_name}_{sym}+ctrl_1+1"
        non_zeros_dict[cname] = all_idx
        top20 = np.arange(20, dtype=np.int64)
        top20_de_dict[cname] = top20
        rgg_dict[cname] = np.array([uniq_gnames[j] for j in top20])
    adata.uns["non_zeros_gene_idx"] = non_zeros_dict
    adata.uns["top_non_dropout_de_20"] = top20_de_dict
    adata.uns["rank_genes_groups_cov_all"] = rgg_dict

    # Split training conditions 90/10 into GEARS internal train / internal val (zero test leakage)
    perm_tr = rng.permutation(len(tr_conds))
    n_ival = max(10, int(round(0.10 * len(tr_conds))))
    ival_conds = [tr_conds[i] for i in perm_tr[:n_ival]]
    itr_conds = [tr_conds[i] for i in perm_tr[n_ival:]] + ["ctrl"]
    split_dict = {"train": itr_conds, "val": ival_conds, "test": test_conds}

    tmp_root = Path(tempfile.mkdtemp(prefix=f"gears_{src_name}_s{seed}_"))
    try:
        # Symlink cached GO files into tmp_root so PertData reuses them without re-downloading
        for fname in ["gene2go_all.pkl", "essential_all_data_pert_genes.pkl"]:
            src_f = gears_cache / fname
            if src_f.exists():
                os.symlink(src_f, tmp_root / fname)

        pd_obj = PertData(str(tmp_root), default_pert_graph=True)
        # Skip scanpy rank_genes_groups DE recomputation (we already set exact train-only DE dicts)
        pd_obj.new_data_process(dataset_name="split_ds", adata=adata, skip_calc_de=True)
        # Copy back any newly downloaded GO pickle to gears_cache for reuse
        for fname in ["gene2go_all.pkl", "essential_all_data_pert_genes.pkl"]:
            f_tmp = tmp_root / fname
            f_dst = gears_cache / fname
            if f_tmp.exists() and not f_tmp.is_symlink() and not f_dst.exists():
                shutil.copy2(f_tmp, f_dst)

        # Ensure non_zeros_gene_idx / top_non_dropout_de_20 are intact after new_data_process
        pd_obj.adata.uns["non_zeros_gene_idx"] = non_zeros_dict
        pd_obj.adata.uns["top_non_dropout_de_20"] = top20_de_dict

        split_pkl = tmp_root / "custom_split.pkl"
        with open(split_pkl, "wb") as f:
            pickle.dump(split_dict, f)
        pd_obj.prepare_split(split="custom", seed=seed, split_dict_path=str(split_pkl))
        pd_obj.get_dataloader(batch_size=32, test_batch_size=64)
        # Pop dummy test_loader so GEARS.train() cleanly sets self.best_model and returns
        # before calling deeper_analysis on placeholder test rows (zero effect on training/val)
        pd_obj.dataloader.pop("test_loader", None)

        gears_model = GEARS(pd_obj, device=device, weight_bias_track=False)
        # Official default hyperparameters from Roohani et al. (2023) / cell-gears==0.1.2
        gears_model.model_initialize(hidden_size=64)
        gears_model.train(epochs=epochs, lr=1e-3, weight_decay=5e-4)
        assert gears_model.best_model is not None, "GEARS best_model was not set"

        # Predict held-out GO-covered targets via official GEARS.predict()
        preds_out = {}
        uniq_q_syms = sorted(set(sym for _, sym in q_in_go))
        if len(uniq_q_syms) > 0:
            raw_pred = gears_model.predict([[s] for s in uniq_q_syms])
            for rid, sym in q_in_go:
                preds_out[rid] = np.asarray(raw_pred[sym], dtype=np.float64)
        for rid, _ in q_out_go:
            preds_out[rid] = mu_train.copy()

        cov = {
            "n_train_in_go": len(tr_rows),
            "n_train_total": len(own_rep),
            "n_test_in_go": len(q_in_go),
            "n_test_total": len(q_df),
            "test_go_coverage": round(len(q_in_go) / max(len(q_df), 1), 4),
            "in_go_rids": [rid for rid, _ in q_in_go],
        }
        return preds_out, cov
    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)
        gc.collect()


def run_benchmark(dataset: str, seeds: list[int], device: str = "cpu", epochs: int = 20):
    t0 = time.time()
    print(f"=== R7-2 Official GEARS (cell-gears==0.1.2, default hyperparams) on {dataset} | seeds={seeds} | device={device} ===", flush=True)
    train_df, all_y, sh0, sh1, g2i_map, sidx, gnames_map = load_bundle(dataset)
    model = PanPerturbBridgeModel(g2i_map, sh0, sh1, sidx)
    gears_cache = resolve_gears_cache_dir()

    rows = []
    ci_rows = []
    cov_records = []

    for seed in seeds:
        ts = time.time()
        val_df, df_tr, _, bases, axes = make_strict_union_split(train_df, g2i_map, all_y, sidx, seed)
        df_hct = df_tr[df_tr["source"] == "HCT116"]
        df_hek = df_tr[df_tr["source"] == "HEK293T"]
        q_hct = val_df[val_df["source"] == "HCT116"]
        q_hek = val_df[val_df["source"] == "HEK293T"]

        # Run official GEARS on both cell lines
        p_gears_0, cov0 = run_official_gears_one_line(
            "HCT116", df_hct, q_hct, all_y, bases["HCT116"], gnames_map["HCT116"], seed, gears_cache, device=device, epochs=epochs
        )
        p_gears_1, cov1 = run_official_gears_one_line(
            "HEK293T", df_hek, q_hek, all_y, bases["HEK293T"], gnames_map["HEK293T"], seed, gears_cache, device=device, epochs=epochs
        )
        p_gears = {**p_gears_0, **p_gears_1}
        cov_records.append({"seed": seed, "line0": {k: v for k, v in cov0.items() if k != "in_go_rids"}, "line1": {k: v for k, v in cov1.items() if k != "in_go_rids"}})

        # Compute S5, S7, AE on exact same split for exact paired row-level CIs
        r0 = model.predict_for_source_all_stages("HCT116", df_hct, df_hek, q_hct, all_y, bases, axes, ("S5", "S7"))
        r1 = model.predict_for_source_all_stages("HEK293T", df_hek, df_hct, q_hek, all_y, bases, axes, ("S5", "S7"))
        p_s5 = {**r0["S5"], **r1["S5"]}
        p_s7 = {**r0["S7"], **r1["S7"]}

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

        preds_map = {"Official_GEARS": p_gears, "S5": p_s5, "S7": p_s7, "AE": p_ae}
        rm = compute_row_metrics_for_split(preds_map, val_df, all_y, bases, axes, sidx, g2i_map)
        cm = compute_centroid_rows(preds_map, df_tr, val_df, all_y, sidx, g2i_map)
        src_arr = val_df["source"].values
        i0, i1 = np.where(src_arr == "HCT116")[0], np.where(src_arr == "HEK293T")[0]
        v5_map = {m: row_scores_v5(preds_map[m], val_df, all_y, bases, axes, sidx, g2i_map) for m in preds_map}

        ev_g = evaluate_5metrics(p_gears, val_df, all_y, bases, axes, sidx, g2i_map)
        dr_raw_g = rm["Official_GEARS"]["disc_raw"]
        dr_v5_g = rm["Official_GEARS"]["disc_v5"]
        ct_g = cm["Official_GEARS"]["centroid_train_r"]

        # Also evaluate GO-covered subset only
        in_go_set = set(cov0["in_go_rids"] + cov1["in_go_rids"])
        go_mask = val_df["id"].isin(in_go_set).values
        i0_go = np.where((src_arr == "HCT116") & go_mask)[0]
        i1_go = np.where((src_arr == "HEK293T") & go_mask)[0]
        v5_go_only = float(0.5 * (v5_map["Official_GEARS"][i0_go].mean() + v5_map["Official_GEARS"][i1_go].mean()))
        disc_go_only = float(0.5 * (dr_raw_g[i0_go].mean() + dr_raw_g[i1_go].mean()))

        rows.append({
            "dataset": dataset,
            "seed": seed,
            "method": "Official_GEARS",
            **{k: round(v, 6) for k, v in ev_g.items()},
            "disc_raw": round(float(0.5 * (dr_raw_g[i0].mean() + dr_raw_g[i1].mean())), 6),
            "disc_v5": round(float(0.5 * (dr_v5_g[i0].mean() + dr_v5_g[i1].mean())), 6),
            "centroid_train_r": round(float(0.5 * (ct_g[i0].mean() + ct_g[i1].mean())), 6),
            "go_coverage_line0": cov0["test_go_coverage"],
            "go_coverage_line1": cov1["test_go_coverage"],
            "v5_go_covered_only": round(v5_go_only, 6),
            "disc_raw_go_covered_only": round(disc_go_only, 6),
        })

        for comp, ka, kb in [
            ("S5_vs_Official_GEARS", "S5", "Official_GEARS"),
            ("S7_vs_Official_GEARS", "S7", "Official_GEARS"),
            ("Official_GEARS_vs_AE", "Official_GEARS", "AE"),
        ]:
            ci_rows.append({
                "dataset": dataset, "seed": seed, "metric": "v5", "comparison": comp,
                **stratified_row_bootstrap_ci(v5_map[ka], v5_map[kb], src_arr),
            })
            for met in ["delta_r", "disc_raw", "disc_v5"]:
                ci_rows.append({
                    "dataset": dataset, "seed": seed, "metric": met, "comparison": comp,
                    **stratified_row_bootstrap_ci(rm[ka][met], rm[kb][met], src_arr),
                })
            ci_rows.append({
                "dataset": dataset, "seed": seed, "metric": "centroid_train_r", "comparison": comp,
                **stratified_row_bootstrap_ci(cm[ka]["centroid_train_r"], cm[kb]["centroid_train_r"], src_arr),
            })

        print(
            f"  [GEARS {dataset} seed={seed}] v5={ev_g['v5']:.4f} (GO-only {v5_go_only:.4f}) | "
            f"delta_r={ev_g['delta_r']:.4f} | disc_raw={rows[-1]['disc_raw']:.4f} | "
            f"centroid_r={rows[-1]['centroid_train_r']:.4f} | cov=({cov0['test_go_coverage']:.2f},{cov1['test_go_coverage']:.2f}) [{time.time()-ts:.1f}s]",
            flush=True,
        )
        # Incremental persist after every completed split
        pd.DataFrame(rows).to_csv(REPORTS_DIR / f"r72_official_gears_{dataset}_metrics.csv", index=False)
        pd.DataFrame(ci_rows).to_csv(REPORTS_DIR / f"r72_official_gears_{dataset}_row_ci.csv", index=False)
        gc.collect()

    df_m = pd.DataFrame(rows)
    df_ci = pd.DataFrame(ci_rows)
    tag = f"r72_official_gears_{dataset}"
    df_m.to_csv(REPORTS_DIR / f"{tag}_metrics.csv", index=False)
    df_ci.to_csv(REPORTS_DIR / f"{tag}_row_ci.csv", index=False)

    summary = {
        "dataset": dataset,
        "package": "cell-gears==0.1.2",
        "hyperparameters": {"hidden_size": 64, "epochs": epochs, "lr": 1e-3, "weight_decay": 5e-4, "batch_size": 32},
        "seeds": seeds,
        "coverage": cov_records,
        "metrics": {
            col: {
                "mean": round(float(df_m[col].mean()), 6),
                "std": round(float(df_m[col].std(ddof=1)), 6) if len(df_m) > 1 else 0.0,
            }
            for col in ["v5", "delta_r", "de50_r", "de50_sign", "disc_raw", "disc_v5", "discrim_delta_pct", "centroid_train_r", "v5_go_covered_only", "disc_raw_go_covered_only"]
        },
        "comparisons": {},
    }
    for met in ["v5", "delta_r", "disc_raw", "disc_v5", "centroid_train_r"]:
        summary["comparisons"][met] = {}
        for comp in ["S5_vs_Official_GEARS", "S7_vs_Official_GEARS", "Official_GEARS_vs_AE"]:
            sub = df_ci[(df_ci["metric"] == met) & (df_ci["comparison"] == comp)]
            summary["comparisons"][met][comp] = {
                "mean_delta": round(float(sub["delta"].mean()), 6),
                "mean_ci95_lo": round(float(sub["ci95_lo"].mean()), 6),
                "mean_ci95_hi": round(float(sub["ci95_hi"].mean()), 6),
                "splits_ci_gt_zero": f"{int(sub['ci_gt_zero'].sum())}/{len(sub)}",
                "splits_ci_lt_zero": f"{int((sub['ci95_hi'] < 0).sum())}/{len(sub)}",
            }

    with open(REPORTS_DIR / f"{tag}_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nSaved {tag}_metrics.csv / _row_ci.csv / _summary.json in {time.time()-t0:.1f}s!", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", choices=["replogle", "nadig"], required=True)
    ap.add_argument("--seeds", type=int, nargs="+", default=SEEDS_10)
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--epochs", type=int, default=20)
    args = ap.parse_args()
    run_benchmark(args.dataset, args.seeds, device=args.device, epochs=args.epochs)
