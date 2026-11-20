#!/usr/bin/env python3
"""
T5: Local ReactomePathways.gmt Hypergeometric Enrichment + Benjamini-Hochberg FDR
using the EXACT Evaluated Gene Universe as Background (Replogle 2022 ONLY).

Fixes the Round 5 audit flaw where online Reactome used the ~20,000 whole-genome
background while the evaluated Replogle essential/reproducible gene set is already
strongly pre-enriched for housekeeping complexes.

Mandatory Comparisons (all against the evaluated gene universe background):
  1. Top-N genes by S7 accuracy (r_S7 >= 0.45) vs Top-N genes by S0 accuracy (identical N)
  2. Top-150 genes by (S7 - S0) gain
  3. Top-150 genes by (S6 - S5) gain
  4. (Auxiliary) Top-150 genes by (S5 - S2) Cross-Line Bridge gain
"""

import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

PROJECT_DIR = Path(__file__).resolve().parents[1]
REPORTS_DIR = PROJECT_DIR / "reports"


def resolve_gmt_path() -> Path:
    cands = [
        *( [Path(os.environ["REACTOME_GMT_PATH"])] if os.environ.get("REACTOME_GMT_PATH") else [] ),
        Path.cwd() / "external_data/reference/ReactomePathways.gmt",
        Path(__file__).resolve().parents[4] / "external_data/reference/ReactomePathways.gmt",
        Path(__file__).resolve().parents[3] / "external_data/reference/ReactomePathways.gmt",
    ]
    return next((c for c in cands if c.exists()), cands[0])


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_gmt(path: Path) -> dict[str, tuple[str, set[str]]]:
    """Return {stId: (pathway_name, set_of_hgnc_symbols)}."""
    pathways = {}
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name, stid = parts[0], parts[1]
            genes = set(g.strip() for g in parts[2:] if g.strip())
            if genes:
                pathways[stid] = (name, genes)
    return pathways


def hypergeom_bh(query_genes: set[str], background_genes: set[str], pathways: dict, min_bg_overlap: int = 5) -> pd.DataFrame:
    """
    One-sided hypergeometric overrepresentation test strictly within `background_genes`.
    M = |background_genes|
    n = |pathway ∩ background_genes|
    N = |query_genes ∩ background_genes|
    k = |query_genes ∩ pathway ∩ background_genes|
    """
    bg = set(background_genes)
    q = set(query_genes) & bg
    M = len(bg)
    N = len(q)
    rows = []
    for stid, (name, pgenes) in pathways.items():
        p_bg = pgenes & bg
        n = len(p_bg)
        if n < min_bg_overlap:
            continue
        k = len(q & p_bg)
        if k == 0:
            pval = 1.0
        else:
            pval = float(hypergeom.sf(k - 1, M, n, N))
        expected = (N * n) / max(M, 1)
        fold = k / expected if expected > 0 else 0.0
        rows.append({
            "stId": stid,
            "name": name,
            "k_in_query": int(k),
            "N_query": int(N),
            "n_in_bg": int(n),
            "M_bg": int(M),
            "fold_enrichment": round(fold, 4),
            "p_value": pval,
        })
    df = pd.DataFrame(rows)
    if len(df) > 0:
        _, fdrs, _, _ = multipletests(df["p_value"].values, method="fdr_bh")
        df["fdr_bh"] = fdrs
        df = df.sort_values(["fdr_bh", "p_value", "fold_enrichment"], ascending=[True, True, False]).reset_index(drop=True)
    return df


def main():
    gmt_path = resolve_gmt_path()
    if not gmt_path.exists():
        raise FileNotFoundError(f"ReactomePathways.gmt not found at {gmt_path}")
    gmt_sha = sha256_file(gmt_path)
    pathways = load_gmt(gmt_path)

    df_oof = pd.read_csv(REPORTS_DIR / "replogle2022_5fold_oof_per_gene_mechanisms.csv")
    # Aggregate per unique HGNC gene symbol across RPE1 & K562 (mean score per unique evaluated gene)
    df_gene = df_oof.groupby("gene_symbol", as_index=False).agg({
        "r_S0": "mean",
        "r_S2": "mean",
        "r_S5": "mean",
        "r_S6": "mean",
        "r_S7": "mean",
        "gain_total_S7_S0": "mean",
        "gain_bridge_S5_S2": "mean",
        "gain_recip_S6_S5": "mean",
    })
    bg_all = set(df_gene["gene_symbol"].tolist())
    # Also compute GMT-annotated subset of evaluated genes as secondary sensitivity check
    gmt_universe = set().union(*[g for _, g in pathways.values()])
    bg_gmt = bg_all & gmt_universe

    # Comparison 1: Same N for S7 top vs S0 top (N = count of unique genes with mean r_S7 >= 0.45, or matching raw row count)
    n_top = int((df_gene["r_S7"] >= 0.45).sum())
    if n_top < 50:
        n_top = 200
    s7_top = set(df_gene.sort_values("r_S7", ascending=False).head(n_top)["gene_symbol"])
    s0_top = set(df_gene.sort_values("r_S0", ascending=False).head(n_top)["gene_symbol"])
    overlap_s7_s0 = len(s7_top & s0_top)

    # Comparison 2: Top-150 by S7 - S0 gain
    s7_minus_s0_150 = set(df_gene.sort_values("gain_total_S7_S0", ascending=False).head(150)["gene_symbol"])
    # Comparison 3: Top-150 by S6 - S5 gain
    s6_minus_s5_150 = set(df_gene.sort_values("gain_recip_S6_S5", ascending=False).head(150)["gene_symbol"])
    # Auxiliary: Top-150 by S5 - S2 gain
    s5_minus_s2_150 = set(df_gene.sort_values("gain_bridge_S5_S2", ascending=False).head(150)["gene_symbol"])

    print(f"=== T5: Local Reactome GMT Enrichment (Evaluated Universe Background) ===")
    print(f"GMT: {gmt_path.name} (SHA-256: {gmt_sha[:16]}..., {len(pathways)} pathways)")
    print(f"Evaluated unique genes (M_bg): {len(bg_all)} (in Reactome GMT: {len(bg_gmt)})")
    print(f"Top-N comparison N={n_top}: |S7_top ∩ S0_top| = {overlap_s7_s0}/{n_top} ({100*overlap_s7_s0/n_top:.1f}%)\n")

    cohorts = {
        f"S7_top_{n_top}": s7_top,
        f"S0_top_{n_top}": s0_top,
        "S7_minus_S0_top150": s7_minus_s0_150,
        "S6_minus_S5_top150": s6_minus_s5_150,
        "S5_minus_S2_top150": s5_minus_s2_150,
    }

    summary = {
        "gmt_file": gmt_path.name,
        "gmt_license": "CC BY 4.0 (https://reactome.org/license)",
        "gmt_sha256": gmt_sha,
        "n_gmt_pathways": len(pathways),
        "n_evaluated_unique_genes": len(bg_all),
        "n_evaluated_in_gmt": len(bg_gmt),
        "top_N": n_top,
        "overlap_S7_top_vs_S0_top": overlap_s7_s0,
        "cohorts_full_evaluated_bg": {},
        "cohorts_gmt_annotated_bg": {},
    }

    all_tables = []
    for cname, qset in cohorts.items():
        df_res = hypergeom_bh(qset, bg_all, pathways, min_bg_overlap=5)
        n_sig = int((df_res["fdr_bh"] < 0.05).sum()) if len(df_res) > 0 else 0
        top10 = df_res.head(10).to_dict(orient="records")
        summary["cohorts_full_evaluated_bg"][cname] = {
            "n_query": len(qset),
            "n_pathways_tested": len(df_res),
            "n_fdr_lt_0_05": n_sig,
            "top10": top10,
        }
        df_res_g = hypergeom_bh(qset, bg_gmt, pathways, min_bg_overlap=5)
        n_sig_g = int((df_res_g["fdr_bh"] < 0.05).sum()) if len(df_res_g) > 0 else 0
        summary["cohorts_gmt_annotated_bg"][cname] = {
            "n_query_in_gmt": len(qset & bg_gmt),
            "n_pathways_tested": len(df_res_g),
            "n_fdr_lt_0_05": n_sig_g,
            "top5": df_res_g.head(5).to_dict(orient="records"),
        }
        df_save = df_res.head(20).copy()
        df_save.insert(0, "cohort", cname)
        all_tables.append(df_save)
        top_str = f"{top10[0]['name']} (k={top10[0]['k_in_query']}/{top10[0]['n_in_bg']}, fold={top10[0]['fold_enrichment']}x, FDR={top10[0]['fdr_bh']:.2e})" if top10 else "None"
        print(f"  {cname:22s} | FDR<0.05 (full bg): {n_sig:3d} (gmt bg: {n_sig_g:3d}) | #1: {top_str}")

    # Direct verdict per T5 pre-registered rule
    sig_s7s0 = summary["cohorts_full_evaluated_bg"]["S7_minus_S0_top150"]["n_fdr_lt_0_05"]
    sig_s6s5 = summary["cohorts_full_evaluated_bg"]["S6_minus_S5_top150"]["n_fdr_lt_0_05"]
    summary["t5_verdict"] = (
        "GAIN_ENRICHMENT_FOUND" if (sig_s7s0 > 0 or sig_s6s5 > 0)
        else "NO_SIGNIFICANT_GAIN_ENRICHMENT_DELETE_MECHANISM_VALIDATION_CLAIM"
    )
    print(f"\nT5 VERDICT: {summary['t5_verdict']} (S7-S0 FDR<0.05: {sig_s7s0}, S6-S5 FDR<0.05: {sig_s6s5})")

    pd.concat(all_tables, ignore_index=True).to_csv(REPORTS_DIR / "t5_local_reactome_enrichment.csv", index=False)
    with open(REPORTS_DIR / "t5_local_reactome_enrichment_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("Saved reports/t5_local_reactome_enrichment.csv & summary JSON!")


if __name__ == "__main__":
    main()
