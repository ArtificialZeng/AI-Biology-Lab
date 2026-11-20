#!/usr/bin/env python3
"""
Real Named-Gene Biological Mechanism Discovery & Live Reactome Pathway Enrichment
on Replogle et al. (Cell 2022) K562 & RPE1 Genome-Scale Perturb-seq.

Because BioXArena-Pro `xatlas-orion-perturb` anonymizes gene names into cryptographic
hashes (`g_<hex>`), external biological pathway enrichment cannot be run on Orion hashes
(which simultaneously guarantees zero external database leakage on Orion).

Here, we run a leak-free 5-fold strict union out-of-fold audit of frozen PanPerturb-Bridge
on Replogle et al. (Cell 2022) with REAL human HGNC gene symbols & Ensembl IDs, then query
the live Reactome Analysis Service (`Homo sapiens`, FDR < 0.05) to test which curated human
pathways and macromolecular complexes drive:
  1. Cross-Line Bridge gains (S5 - S2) across myeloid leukemia (K562) and non-cancer RPE1;
  2. Thresholded Reciprocal Causal Law gains (S6 - S5), including exact unmeasured
     panel-outside subunit knockdowns that unlock held-out complex partners;
  3. High-accuracy vs. hard lineage-specific biological boundaries.
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

from replogle_external_validation import REPLOGLE_DIR, load_replogle_cell_line
from xatlas_panhuman_solver import (
    REPORTS_DIR,
    PanPerturbBridgeModel,
    _project_out,
    fit_baseline_and_axis,
    score_prediction,
)

REACTOME_CLI = Path(os.environ["REACTOME_CLI_PATH"]) if os.environ.get("REACTOME_CLI_PATH") else None


def run_reactome_enrichment(genes: list[str], out_json: Path, label: str):
    if not genes or REACTOME_CLI is None or not REACTOME_CLI.exists():
        return []
    clean_genes = sorted(set(g for g in genes if g and not g.startswith("g_")))
    gene_str = ",".join(clean_genes)
    cmd = [
        "uv", "run", str(REACTOME_CLI),
        "analyze-projection",
        "--data", gene_str,
        "--fdr", "0.05",
        "--limit", "25",
        "--output", str(out_json),
    ]
    env = os.environ.copy()
    for k in ("http_proxy", "https_proxy", "HTTP_PROXY", "HTTPS_PROXY"):
        env.pop(k, None)
    res = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=90)
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
    if top:
        print(f"  [Reactome] {label}: {len(clean_genes)} genes -> {len(pathways)} enriched pathways (FDR<0.05); #1: {top[0]['name']} (found={top[0]['found']}/{top[0]['total']}, FDR={top[0]['fdr']:.2e})", flush=True)
    else:
        print(f"  [Reactome] {label}: {len(clean_genes)} genes -> 0 pathways at FDR<0.05", flush=True)
    return top


def get_projected_vectors(true, pred, baseline, axis, s_idx, t_idx):
    indices = s_idx[s_idx != t_idx] if (t_idx is not None) else s_idx
    Q, Rm = np.linalg.qr(np.stack([np.ones(len(indices)), baseline[indices], axis[indices]], axis=1))
    Q = Q[:, np.abs(np.diag(Rm)) > 1e-10 * max(1.0, np.abs(Rm).max())]
    a = _project_out(true[indices] - baseline[indices], Q)
    b = _project_out(pred[indices] - baseline[indices], Q)
    return indices, a, b


def main():
    t0 = time.time()
    print("=== Replogle 2022 Real Named-Gene Mechanism Discovery & Live Reactome Audit ===", flush=True)

    df_rpe, y_rpe, rpe_genes, rpe_gmap = load_replogle_cell_line(REPLOGLE_DIR / "rpe1_normalized_bulk_01.h5ad", "HCT116", 5000)
    df_k562, y_k562, k562_genes, k562_gmap = load_replogle_cell_line(REPLOGLE_DIR / "K562_gwps_normalized_bulk_01.h5ad", "HEK293T", 5000)

    train_df = pd.concat([df_rpe, df_k562], ignore_index=True)
    all_y = {**y_rpe, **y_k562}
    ens2sym = dict(zip(train_df["perturbation"], train_df["gene_symbol"]))
    ens2sym.update(rpe_gmap)
    ens2sym.update(k562_gmap)

    shared_genes = sorted(list(set(rpe_genes) & set(k562_genes)))
    rpe_sh_idx = np.array([rpe_genes.index(g) for g in shared_genes])
    k562_sh_idx = np.array([k562_genes.index(g) for g in shared_genes])
    gene2idx = {
        "HCT116": {g: i for i, g in enumerate(rpe_genes)},
        "HEK293T": {g: i for i, g in enumerate(k562_genes)},
    }
    genes_by_slot = {"HCT116": rpe_genes, "HEK293T": k562_genes}
    slot2line = {"HCT116": "RPE1", "HEK293T": "K562"}
    score_indices = {
        "HCT116_indices": np.arange(5000, dtype=np.int64),
        "HEK293T_indices": np.arange(5000, dtype=np.int64),
    }

    model = PanPerturbBridgeModel(gene2idx, rpe_sh_idx, k562_sh_idx, score_indices)
    stages = ("S0", "S2", "S5", "S6", "S7")

    # 5-fold strict union partition over ALL unique strong in-panel genes across RPE1 & K562
    rng = np.random.RandomState(42)
    all_strong_in = sorted(list(
        set(df_rpe[(df_rpe["reproducible"] == 1) & (df_rpe["perturbation"].isin(gene2idx["HCT116"].keys()))]["perturbation"])
        | set(df_k562[(df_k562["reproducible"] == 1) & (df_k562["perturbation"].isin(gene2idx["HEK293T"].keys()))]["perturbation"])
    ))
    folds = np.array_split(rng.permutation(all_strong_in), 5)

    oof_rows = []
    case_store = {}
    for fold in range(5):
        tf = time.time()
        val_union = set(folds[fold])
        df_tr = train_df[~train_df["perturbation"].isin(val_union)].copy()
        assert len(set(df_tr["perturbation"]) & val_union) == 0

        bases, axes = {}, {}
        for slot in ["HCT116", "HEK293T"]:
            rep_s = df_tr[(df_tr["source"] == slot) & (df_tr["reproducible"] == 1)]
            bases[slot], axes[slot] = fit_baseline_and_axis(rep_s, all_y, score_indices[f"{slot}_indices"])

        for slot in ["HCT116", "HEK293T"]:
            oth = "HEK293T" if slot == "HCT116" else "HCT116"
            g2i = gene2idx[slot]
            q_df = train_df[
                (train_df["source"] == slot)
                & (train_df["reproducible"] == 1)
                & (train_df["perturbation"].isin(val_union))
                & (train_df["perturbation"].isin(g2i.keys()))
            ].copy()
            if len(q_df) == 0:
                continue
            df_own = df_tr[df_tr["source"] == slot]
            df_oth = df_tr[df_tr["source"] == oth]
            res = model.predict_for_source_all_stages(slot, df_own, df_oth, q_df, all_y, bases, axes, stages)

            # Extract named reciprocal donors (|z| > 4.0)
            own_rep = df_own[df_own["reproducible"] == 1].reset_index(drop=True)
            Y_raw = np.array([all_y[rid] for rid in own_rep["id"]], dtype=np.float32)
            Z = Y_raw.copy()
            rep_perts = own_rep["perturbation"].tolist()
            rep_syms = own_rep["gene_symbol"].tolist()
            for i, p in enumerate(rep_perts):
                if p in g2i:
                    Z[i, g2i[p]] = 0.0
            Q5, Rm = np.linalg.qr(np.stack([np.ones(5000), bases[slot], axes[slot]], axis=1))
            Q5 = Q5[:, np.abs(np.diag(Rm)) > 1e-10 * max(1.0, np.abs(Rm).max())].astype(np.float32)
            Z = Z - bases[slot].astype(np.float32)[None, :]
            Z = Z - (Z @ Q5) @ Q5.T
            Z_z = Z / (np.std(Z, axis=0, keepdims=True) + 1e-6)

            for _, r in q_df.iterrows():
                rid = r["id"]
                ens = r["perturbation"]
                sym = r["gene_symbol"]
                t_idx = g2i[ens]
                sc = {}
                for st in stages:
                    sc[st] = score_prediction(all_y[rid], res[st][rid], bases[slot], axes[slot], score_indices[f"{slot}_indices"], t_idx)
                col_z = Z_z[:, t_idx]
                abs_z = np.abs(col_z)
                hits = np.where(abs_z > 4.0)[0]
                hits_s = hits[np.argsort(abs_z[hits])[::-1]]
                n_unm = sum(1 for h in hits_s if rep_perts[h] not in g2i)
                top_donors = [
                    f"{rep_syms[h]}({col_z[h]:+.1f}z{'*' if rep_perts[h] not in g2i else ''})"
                    for h in hits_s[:5]
                ]
                oof_rows.append({
                    "fold": fold,
                    "cell_line": slot2line[slot],
                    "id": rid,
                    "ensembl_id": ens,
                    "gene_symbol": sym,
                    "r_S0": round(sc["S0"], 5),
                    "r_S2": round(sc["S2"], 5),
                    "r_S5": round(sc["S5"], 5),
                    "r_S6": round(sc["S6"], 5),
                    "r_S7": round(sc["S7"], 5),
                    "gain_total_S7_S0": round(sc["S7"] - sc["S0"], 5),
                    "gain_bridge_S5_S2": round(sc["S5"] - sc["S2"], 5),
                    "gain_recip_S6_S5": round(sc["S6"] - sc["S5"], 5),
                    "n_recip_hits_gt4z": int(len(hits_s)),
                    "n_unmeasured_recip_hits": int(n_unm),
                    "top_recip_donors": "; ".join(top_donors),
                })
                if sc["S7"] > 0.50 and (sc["S7"] - sc["S0"]) > 0.25:
                    case_store[rid] = (slot, sym, ens, bases[slot].copy(), axes[slot].copy(), res["S0"][rid].astype(np.float32), res["S7"][rid].astype(np.float32))
        print(f"  Completed Replogle OOF fold {fold+1}/5 in {time.time()-tf:.1f}s", flush=True)
        gc.collect()

    df_oof = pd.DataFrame(oof_rows)
    df_oof.to_csv(REPORTS_DIR / "replogle2022_5fold_oof_per_gene_mechanisms.csv", index=False)

    print("\n=== Replogle 2022 5-Fold Strict-Union OOF Summary ===", flush=True)
    for cl in ["RPE1", "K562"]:
        sub = df_oof[df_oof["cell_line"] == cl]
        print(f"  {cl:4s} (n={len(sub):4d}): S0={sub['r_S0'].mean():.4f} -> S2={sub['r_S2'].mean():.4f} -> S5={sub['r_S5'].mean():.4f} -> S6={sub['r_S6'].mean():.4f} -> S7={sub['r_S7'].mean():.4f}", flush=True)

    # Live Reactome pathway enrichment on REAL HUMAN HGNC GENE SYMBOLS!
    top_bridge = df_oof.sort_values("gain_bridge_S5_S2", ascending=False)["gene_symbol"].drop_duplicates().head(150).tolist()
    top_recip = df_oof[df_oof["n_recip_hits_gt4z"] > 0].sort_values("gain_recip_S6_S5", ascending=False)["gene_symbol"].drop_duplicates().head(150).tolist()
    top_s7 = df_oof[df_oof["r_S7"] >= 0.45]["gene_symbol"].drop_duplicates().tolist()
    hard_genes = df_oof[df_oof["r_S7"] < 0.05].sort_values("r_S7", ascending=True)["gene_symbol"].drop_duplicates().head(150).tolist()

    print("\n--- Querying Live Reactome Analysis Service on Real Human HGNC Symbols ---", flush=True)
    reactome_res = {
        "bridge_winners_S5_vs_S2": run_reactome_enrichment(top_bridge, REPORTS_DIR / "reactome_replogle_bridge_winners.json", "Replogle Bridge Winners (Top 150)"),
        "reciprocal_winners_S6_vs_S5": run_reactome_enrichment(top_recip, REPORTS_DIR / "reactome_replogle_reciprocal_winners.json", "Replogle Reciprocal Winners (Top 150)"),
        "high_accuracy_S7_gt045": run_reactome_enrichment(top_s7, REPORTS_DIR / "reactome_replogle_high_accuracy.json", f"Replogle High-Accuracy S7>=0.45 (n={len(top_s7)})"),
        "hard_boundary_S7_lt005": run_reactome_enrichment(hard_genes, REPORTS_DIR / "reactome_replogle_hard_boundary.json", "Replogle Hard Boundary S7<0.05 (Top 150)"),
    }

    # Named real-gene showcase case studies
    case_studies = []
    for cl in ["RPE1", "K562"]:
        cand = df_oof[(df_oof["cell_line"] == cl) & (df_oof["id"].isin(case_store.keys())) & (df_oof["n_unmeasured_recip_hits"] > 0)].sort_values("gain_total_S7_S0", ascending=False)
        if len(cand) < 2:
            cand = df_oof[(df_oof["cell_line"] == cl) & (df_oof["id"].isin(case_store.keys()))].sort_values("gain_total_S7_S0", ascending=False)
        for _, row in cand.head(2).iterrows():
            rid = row["id"]
            slot, sym, ens, base_s, ax_s, p_s0, p_s7 = case_store[rid]
            t_idx = gene2idx[slot][ens]
            idx, a_true, b_s0 = get_projected_vectors(all_y[rid], p_s0, base_s, ax_s, score_indices[f"{slot}_indices"], t_idx)
            _, _, b_s7 = get_projected_vectors(all_y[rid], p_s7, base_s, ax_s, score_indices[f"{slot}_indices"], t_idx)
            top_de = np.argsort(np.abs(a_true))[-8:][::-1]
            top_de_syms = [ens2sym.get(genes_by_slot[slot][idx[i]], genes_by_slot[slot][idx[i]]) for i in top_de]
            case_studies.append({
                "cell_line": cl, "id": rid, "gene_symbol": sym, "ensembl_id": ens,
                "r_S0": float(row["r_S0"]), "r_S2": float(row["r_S2"]),
                "r_S5": float(row["r_S5"]), "r_S6": float(row["r_S6"]), "r_S7": float(row["r_S7"]),
                "top_recip_donors": row["top_recip_donors"],
                "top_de_genes": top_de_syms,
                "a_true": [round(float(x), 4) for x in a_true],
                "b_s0": [round(float(x), 4) for x in b_s0],
                "b_s7": [round(float(x), 4) for x in b_s7],
            })
            print(f"  Named Case Study [{cl}] {sym}: S0={row['r_S0']:.3f} -> S5={row['r_S5']:.3f} -> S6={row['r_S6']:.3f} -> S7={row['r_S7']:.3f} | donors={row['top_recip_donors']}", flush=True)

    summary = {
        "dataset": "Replogle et al. (Cell 2022) K562 + RPE1",
        "n_total_oof": int(len(df_oof)),
        "oof_means_by_line": {
            cl: {st: round(float(df_oof[df_oof["cell_line"] == cl][f"r_{st}"].mean()), 4) for st in stages}
            for cl in ["RPE1", "K562"]
        },
        "reactome_enrichments": reactome_res,
        "case_studies_meta": [{k: v for k, v in c.items() if k not in ("a_true", "b_s0", "b_s7")} for c in case_studies],
        "runtime_seconds": round(time.time() - t0, 2),
    }
    with open(REPORTS_DIR / "replogle2022_biological_mechanism_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    with open(REPORTS_DIR / "replogle2022_case_study_vectors.json", "w") as f:
        json.dump(case_studies, f)
    print(f"Saved Replogle real-gene Reactome mechanism audit in {time.time()-t0:.1f}s!", flush=True)


if __name__ == "__main__":
    main()
