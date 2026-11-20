#!/usr/bin/env python3
"""
Biological Mechanism Discovery & Live Reactome Pathway Validation for PanPerturb-Bridge.

Runs a 5-fold leak-free out-of-fold evaluation covering ALL reproducible in-panel
target genes in HCT116 (n=646) and HEK293T (n=1,236) to quantify per-gene gains from:
  - S0 (Official 100-SVD + Ridge)
  - S2 (+ Multi-Scale Spectral Deflation Pyramid, single-line)
  - S5 (+ Pan-Human Cross-Line Bridge)
  - S6 (+ Thresholded Reciprocal Causal Law, recording exact unmeasured donor genes k with |z_{k,g}| > 4.0)
  - S7 (Full PanPerturb-Bridge)

Then queries the live Reactome Analysis Service (Homo sapiens, projection to human pathways,
FDR < 0.05) on:
  1. Top Pan-Human Cross-Line Bridge winners (S5 - S2)
  2. Top Reciprocal Causal Law winners (S6 - S5)
  3. Top Full PanPerturb-Bridge high-accuracy targets (r_S7 >= 0.30)
  4. Lineage-specific / low-predictability boundary targets (r_S7 < 0.02)
And saves concrete single-gene case-study projected transcriptome vectors for publication figures.
"""

import gc
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from xatlas_panhuman_solver import (
    DATA_DIR,
    REPORTS_DIR,
    PanPerturbBridgeModel,
    _project_out,
    score_prediction,
)

REACTOME_CLI = Path(os.environ["REACTOME_CLI_PATH"]) if os.environ.get("REACTOME_CLI_PATH") else None


def get_projected_vectors(true, pred, baseline, axis, scoring_indices, target_idx=None):
    indices = scoring_indices
    if target_idx is not None and target_idx in indices:
        indices = indices[indices != target_idx]
    Q, Rm = np.linalg.qr(np.stack([np.ones(len(indices)), baseline[indices], axis[indices]], axis=1))
    Q = Q[:, np.abs(np.diag(Rm)) > 1e-10 * max(1.0, np.abs(Rm).max())]
    a = _project_out(true[indices] - baseline[indices], Q)
    b = _project_out(pred[indices] - baseline[indices], Q)
    return indices, a, b


def run_reactome_enrichment(genes: list[str], out_json: Path, label: str):
    if not genes or REACTOME_CLI is None or not REACTOME_CLI.exists():
        return []
    gene_str = ",".join(sorted(set(genes)))
    cmd = [
        "uv", "run", str(REACTOME_CLI),
        "analyze-projection",
        "--data", gene_str,
        "--fdr", "0.05",
        "--limit", "25",
        "--output", str(out_json),
    ]
    env = os.environ.copy()
    env.pop("http_proxy", None)
    env.pop("https_proxy", None)
    env.pop("HTTP_PROXY", None)
    env.pop("HTTPS_PROXY", None)
    res = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=60)
    if res.returncode != 0 or not out_json.exists():
        print(f"  [Reactome WARN] {label}: rc={res.returncode}, stderr={res.stderr[:200]}", flush=True)
        return []
    with open(out_json) as f:
        data = json.load(f)
    pathways = data.get("pathways", [])
    top = []
    for p in pathways[:15]:
        ent = p.get("entities", {})
        top.append({
            "stId": p.get("stId"),
            "name": p.get("name"),
            "found": ent.get("found"),
            "total": ent.get("total"),
            "pValue": ent.get("pValue"),
            "fdr": ent.get("fdr"),
        })
    print(f"  [Reactome] {label}: {len(genes)} query genes -> {len(pathways)} enriched pathways (FDR<0.05); top={top[0]['name'] if top else 'None'} (FDR={top[0]['fdr']:.2e})" if top else f"  [Reactome] {label}: 0 pathways at FDR<0.05", flush=True)
    return top


def main():
    t0 = time.time()
    print("=== Biological Mechanism Discovery & Live Reactome Validation (5-Fold Full OOF) ===", flush=True)

    train_df = pd.read_csv(DATA_DIR / "train.csv")
    genes_npz = np.load(DATA_DIR / "genes.npz")
    baselines_npz = np.load(DATA_DIR / "context_baseline.npz")
    axes_npz = np.load(DATA_DIR / "residual_axis.npz")
    score_idx_npz = np.load(DATA_DIR / "score_gene_indices.npz")
    train_y_npz = np.load(DATA_DIR / "train_y.npz")
    all_train_y = {k: np.asarray(train_y_npz[k], dtype=np.float32) for k in train_y_npz.files}

    hct_genes = list(genes_npz["HCT116_genes"])
    hek_genes = list(genes_npz["HEK293T_genes"])
    genes_by_line = {"HCT116": hct_genes, "HEK293T": hek_genes}
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

    model = PanPerturbBridgeModel(gene2idx, hct_sh_idx, hek_sh_idx, score_indices)
    stages = ("S0", "S2", "S5", "S6", "S7")

    # Build 5-fold partition over reproducible in-panel genes per line
    rng = np.random.RandomState(42)
    folds_by_line = {}
    for source in ["HCT116", "HEK293T"]:
        sub = train_df[
            (train_df["source"] == source)
            & (train_df["reproducible"] == 1)
            & (train_df["perturbation"].isin(gene2idx[source].keys()))
        ]
        uniq = rng.permutation(sub["perturbation"].unique())
        folds_by_line[source] = np.array_split(uniq, 5)

    oof_rows = []
    oof_preds_store = {}  # rid -> (S0, S7) for case studies

    for fold in range(5):
        tf = time.time()
        for source in ["HCT116", "HEK293T"]:
            other_src = "HEK293T" if source == "HCT116" else "HCT116"
            val_t = set(folds_by_line[source][fold])
            q_df = train_df[
                (train_df["source"] == source)
                & (train_df["reproducible"] == 1)
                & (train_df["perturbation"].isin(val_t))
            ].copy()
            val_rids = set(q_df["id"])
            df_own = train_df[(train_df["source"] == source) & (~train_df["id"].isin(val_rids))]
            df_other = train_df[(train_df["source"] == other_src) & (~train_df["perturbation"].isin(val_t))]
            assert len(set(df_own["perturbation"]) & val_t) == 0
            assert len(set(df_other["perturbation"]) & val_t) == 0

            res = model.predict_for_source_all_stages(
                source, df_own, df_other, q_df, all_train_y, baselines_pub, axes_pub, stages
            )

            # Also extract top reciprocal causal donor genes (|z| > 4.0) from df_own reproducible rows
            g2i = gene2idx[source]
            own_rep = df_own[df_own["reproducible"] == 1].reset_index(drop=True)
            Y_raw = np.array([all_train_y[rid] for rid in own_rep["id"]], dtype=np.float32)
            Z = Y_raw.copy()
            rep_perts = own_rep["perturbation"].tolist()
            for i, p in enumerate(rep_perts):
                if p in g2i:
                    Z[i, g2i[p]] = 0.0
            Q5, Rm = np.linalg.qr(np.stack([np.ones(5000), baselines_pub[source], axes_pub[source]], axis=1))
            Q5 = Q5[:, np.abs(np.diag(Rm)) > 1e-10 * max(1.0, np.abs(Rm).max())].astype(np.float32)
            Z = Z - baselines_pub[source].astype(np.float32)[None, :]
            Z = Z - (Z @ Q5) @ Q5.T
            Z_z = Z / (np.std(Z, axis=0, keepdims=True) + 1e-6)

            for _, r in q_df.iterrows():
                rid = r["id"]
                pert = r["perturbation"]
                t_idx = g2i[pert]
                scores = {}
                for st in stages:
                    scores[st] = score_prediction(
                        all_train_y[rid], res[st][rid],
                        baselines_pub[source], axes_pub[source],
                        score_indices[f"{source}_indices"], t_idx,
                    )
                # Find strongest reciprocal causal donors (|z_{k,g}| > 4.0)
                col_z = Z_z[:, t_idx]
                abs_z = np.abs(col_z)
                hits = np.where(abs_z > 4.0)[0]
                hits_sorted = hits[np.argsort(abs_z[hits])[::-1]]
                n_hits = len(hits_sorted)
                n_unmeasured_hits = sum(1 for h in hits_sorted if rep_perts[h] not in g2i)
                top_donors = [
                    f"{rep_perts[h]}({col_z[h]:+.1f}z{'*' if rep_perts[h] not in g2i else ''})"
                    for h in hits_sorted[:5]
                ]
                oof_rows.append({
                    "fold": fold,
                    "source": source,
                    "id": rid,
                    "perturbation": pert,
                    "in_shared_4040": pert in set(shared_genes),
                    "n_cells": int(r["n_cells"]),
                    "r_S0": round(scores["S0"], 5),
                    "r_S2": round(scores["S2"], 5),
                    "r_S5": round(scores["S5"], 5),
                    "r_S6": round(scores["S6"], 5),
                    "r_S7": round(scores["S7"], 5),
                    "gain_total_S7_S0": round(scores["S7"] - scores["S0"], 5),
                    "gain_bridge_S5_S2": round(scores["S5"] - scores["S2"], 5),
                    "gain_recip_S6_S5": round(scores["S6"] - scores["S5"], 5),
                    "n_recip_hits_gt4z": int(n_hits),
                    "n_unmeasured_recip_hits": int(n_unmeasured_hits),
                    "top_recip_donors": "; ".join(top_donors),
                })
                oof_preds_store[rid] = (source, pert, res["S0"][rid].astype(np.float32), res["S7"][rid].astype(np.float32))
        print(f"  Completed OOF fold {fold+1}/5 in {time.time()-tf:.1f}s", flush=True)
        gc.collect()

    df_oof = pd.DataFrame(oof_rows)
    df_oof.to_csv(REPORTS_DIR / "xatlas_5fold_oof_per_gene_mechanisms.csv", index=False)

    # Quantitative mechanism breakdown
    has_unm = df_oof["n_unmeasured_recip_hits"] > 0
    has_any = df_oof["n_recip_hits_gt4z"] > 0
    print("\n=== OOF Per-Gene Mechanism Summary (n=1,882 reproducible perturbations) ===", flush=True)
    for src in ["HCT116", "HEK293T"]:
        sub = df_oof[df_oof["source"] == src]
        print(f"  {src:7s} (n={len(sub):4d}): S0={sub['r_S0'].mean():.4f} -> S2={sub['r_S2'].mean():.4f} -> S5={sub['r_S5'].mean():.4f} -> S6={sub['r_S6'].mean():.4f} -> S7={sub['r_S7'].mean():.4f}")
    print(f"  With >=1 Reciprocal hit (|z|>4, n={has_any.sum()}): mean S6-S5 gain = {df_oof.loc[has_any, 'gain_recip_S6_S5'].mean():+.4f} vs no hit (n={(~has_any).sum()}): {df_oof.loc[~has_any, 'gain_recip_S6_S5'].mean():+.4f}")
    print(f"  With >=1 UNMEASURED Reciprocal donor (n={has_unm.sum()}): mean S6-S5 gain = {df_oof.loc[has_unm, 'gain_recip_S6_S5'].mean():+.4f}, mean S7={df_oof.loc[has_unm, 'r_S7'].mean():.4f} (vs S0={df_oof.loc[has_unm, 'r_S0'].mean():.4f})")

    # Select top gene cohorts for live Reactome pathway enrichment
    top_bridge_genes = df_oof.sort_values("gain_bridge_S5_S2", ascending=False)["perturbation"].drop_duplicates().head(120).tolist()
    top_recip_genes = df_oof[df_oof["n_recip_hits_gt4z"] > 0].sort_values("gain_recip_S6_S5", ascending=False)["perturbation"].drop_duplicates().head(120).tolist()
    top_s7_genes = df_oof[df_oof["r_S7"] >= 0.30]["perturbation"].drop_duplicates().tolist()
    hard_genes = df_oof[df_oof["r_S7"] < 0.02].sort_values("r_S7", ascending=True)["perturbation"].drop_duplicates().head(120).tolist()

    print("\n--- Querying Live Reactome Analysis Service (Homo sapiens) ---", flush=True)
    reactome_results = {
        "bridge_winners_S5_vs_S2": run_reactome_enrichment(top_bridge_genes, REPORTS_DIR / "reactome_bridge_winners.json", "Bridge Winners (Top 120)"),
        "reciprocal_winners_S6_vs_S5": run_reactome_enrichment(top_recip_genes, REPORTS_DIR / "reactome_reciprocal_winners.json", "Reciprocal Winners (Top 120)"),
        "high_accuracy_S7_gt030": run_reactome_enrichment(top_s7_genes, REPORTS_DIR / "reactome_high_accuracy_s7.json", f"High-Accuracy S7>=0.30 (n={len(top_s7_genes)})"),
        "hard_boundary_S7_lt002": run_reactome_enrichment(hard_genes, REPORTS_DIR / "reactome_hard_boundary.json", "Hard Boundary S7<0.02 (Top 120)"),
    }

    # Select 4 showcase single-gene knockdown case studies (2 HCT116 + 2 HEK293T with strong gain & clear unmeasured donor)
    case_studies = []
    for src in ["HCT116", "HEK293T"]:
        cand = df_oof[(df_oof["source"] == src) & (df_oof["n_unmeasured_recip_hits"] > 0) & (df_oof["r_S7"] > 0.45)].sort_values("gain_total_S7_S0", ascending=False)
        if len(cand) < 2:
            cand = df_oof[df_oof["source"] == src].sort_values("gain_total_S7_S0", ascending=False)
        for _, row in cand.head(2).iterrows():
            rid = row["id"]
            _, pert, p_s0, p_s7 = oof_preds_store[rid]
            t_idx = gene2idx[src][pert]
            idx, a_true, b_s0 = get_projected_vectors(all_train_y[rid], p_s0, baselines_pub[src], axes_pub[src], score_indices[f"{src}_indices"], t_idx)
            _, _, b_s7 = get_projected_vectors(all_train_y[rid], p_s7, baselines_pub[src], axes_pub[src], score_indices[f"{src}_indices"], t_idx)
            top_de = np.argsort(np.abs(a_true))[-8:][::-1]
            top_de_names = [genes_by_line[src][idx[i]] for i in top_de]
            case_studies.append({
                "source": src,
                "id": rid,
                "perturbation": pert,
                "r_S0": float(row["r_S0"]),
                "r_S5": float(row["r_S5"]),
                "r_S6": float(row["r_S6"]),
                "r_S7": float(row["r_S7"]),
                "top_recip_donors": row["top_recip_donors"],
                "top_de_genes": top_de_names,
                "a_true": [round(float(x), 5) for x in a_true],
                "b_s0": [round(float(x), 5) for x in b_s0],
                "b_s7": [round(float(x), 5) for x in b_s7],
            })
            print(f"  Case Study [{src}] {pert}: S0={row['r_S0']:.3f} -> S5={row['r_S5']:.3f} -> S6={row['r_S6']:.3f} -> S7={row['r_S7']:.3f} | donors={row['top_recip_donors']}", flush=True)

    summary = {
        "n_total_oof_perturbations": int(len(df_oof)),
        "hct116_n": int((df_oof["source"] == "HCT116").sum()),
        "hek293t_n": int((df_oof["source"] == "HEK293T").sum()),
        "oof_means": {
            "S0": round(float(df_oof.groupby("source")["r_S0"].mean().mean()), 4),
            "S2": round(float(df_oof.groupby("source")["r_S2"].mean().mean()), 4),
            "S5": round(float(df_oof.groupby("source")["r_S5"].mean().mean()), 4),
            "S6": round(float(df_oof.groupby("source")["r_S6"].mean().mean()), 4),
            "S7": round(float(df_oof.groupby("source")["r_S7"].mean().mean()), 4),
        },
        "reciprocal_stratification": {
            "n_with_any_hit_gt4z": int(has_any.sum()),
            "mean_S6_minus_S5_with_hit": round(float(df_oof.loc[has_any, "gain_recip_S6_S5"].mean()), 4),
            "mean_S6_minus_S5_no_hit": round(float(df_oof.loc[~has_any, "gain_recip_S6_S5"].mean()), 4),
            "n_with_unmeasured_donor_gt4z": int(has_unm.sum()),
            "mean_S6_minus_S5_with_unmeasured_donor": round(float(df_oof.loc[has_unm, "gain_recip_S6_S5"].mean()), 4),
            "mean_S0_with_unmeasured_donor": round(float(df_oof.loc[has_unm, "r_S0"].mean()), 4),
            "mean_S7_with_unmeasured_donor": round(float(df_oof.loc[has_unm, "r_S7"].mean()), 4),
        },
        "reactome_enrichments": reactome_results,
        "case_studies_meta": [{k: v for k, v in c.items() if k not in ("a_true", "b_s0", "b_s7")} for c in case_studies],
        "runtime_seconds": round(time.time() - t0, 2),
    }
    with open(REPORTS_DIR / "xatlas_biological_mechanism_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    with open(REPORTS_DIR / "xatlas_case_study_vectors.json", "w") as f:
        json.dump(case_studies, f)
    print(f"\nSaved xatlas_5fold_oof_per_gene_mechanisms.csv, Reactome enrichments, and case studies in {time.time()-t0:.1f}s!", flush=True)


if __name__ == "__main__":
    main()
