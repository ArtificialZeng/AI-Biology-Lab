#!/usr/bin/env python3
"""
T6: CRISPRi Genomic Proximity / Bidirectional Promoter Hypothesis Test
(STRICTLY on Replogle et al., Cell 2022 — NEVER touches xatlas!).

Uses official GENCODE v46 GRCh38 GTF (`gencode.v46.basic.annotation.gtf.gz`) to extract
TSS coordinates for Replogle Ensembl genes, then tests on 5-fold strict-union OOF splits:
  1. TSS distance distribution of all active S6 donor-query pairs (|z| > 4.0) vs 100,000 random pairs;
  2. Proportion of S6 donor weight & positive S6-S5 gains arising from:
     - < 1 kb (bidirectional promoters / overlapping TSS)
     - < 10 kb (CRISPRi local cis-repression window)
     - Same chromosome < 100 kb
     - Same chromosome >= 100 kb vs different chromosome
  3. Breakdown by donor z-score sign (z < -4 repression vs z > +4 activation);
  4. Direct counterfactual ablation: mask out all donors within < 10 kb (and < 100 kb)
     of the query gene's TSS and recompute S6 across all 5 OOF folds to quantify how much
     of S6's gain remains (overall and on top-40 positive S6-S5 genes).
"""

import gc
import gzip
import hashlib
import json
import os
import re
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
    compute_relu_reciprocal,
    fit_baseline_and_axis,
    project_rows_5000,
    score_prediction,
)


def resolve_gtf_path() -> Path:
    cands = [
        *( [Path(os.environ["GRCH38_GTF_PATH"])] if os.environ.get("GRCH38_GTF_PATH") else [] ),
        Path.cwd() / "external_data/reference/gencode.v46.basic.annotation.gtf.gz",
        Path(__file__).resolve().parents[4] / "external_data/reference/gencode.v46.basic.annotation.gtf.gz",
        Path(__file__).resolve().parents[3] / "external_data/reference/gencode.v46.basic.annotation.gtf.gz",
    ]
    return next((c for c in cands if c.exists()), cands[0])


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_gtf_tss(gtf_path: Path) -> dict[str, tuple[str, int, str]]:
    """Return {ensembl_id_no_version: (chrom, tss, strand)} from gene features in GTF."""
    re_gid = re.compile(r'gene_id "(ENSG\d+)')
    tss_map = {}
    with gzip.open(gtf_path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split("\t", 9)
            if len(parts) < 9 or parts[2] != "gene":
                continue
            m = re_gid.search(parts[8])
            if not m:
                continue
            gid = m.group(1)
            chrom = parts[0]
            start, end, strand = int(parts[3]), int(parts[4]), parts[6]
            tss = start if strand == "+" else end
            if gid not in tss_map:
                tss_map[gid] = (chrom, tss, strand)
    return tss_map


def tss_distance(ens_a: str, ens_b: str, tss_map: dict):
    if ens_a not in tss_map or ens_b not in tss_map:
        return None
    ca, ta, _ = tss_map[ens_a]
    cb, tb, _ = tss_map[ens_b]
    if ca != cb:
        return np.inf
    return abs(ta - tb)


def main():
    t0 = time.time()
    gtf_path = resolve_gtf_path()
    if not gtf_path.exists():
        raise FileNotFoundError(f"GRCh38 GTF not found at {gtf_path}")
    gtf_sha = sha256_file(gtf_path)
    print(f"=== T6: CRISPRi Genomic Proximity Test on Replogle 2022 (GTF: {gtf_path.name}) ===", flush=True)
    tss_map = parse_gtf_tss(gtf_path)
    print(f"Parsed {len(tss_map)} GRCh38 gene TSS entries (SHA-256: {gtf_sha[:16]}...)", flush=True)

    df_rpe, y_rpe, rpe_genes, _ = load_replogle_cell_line(REPLOGLE_DIR / "rpe1_normalized_bulk_01.h5ad", "HCT116", 5000)
    df_k562, y_k562, k562_genes, _ = load_replogle_cell_line(REPLOGLE_DIR / "K562_gwps_normalized_bulk_01.h5ad", "HEK293T", 5000)
    train_df = pd.concat([df_rpe, df_k562], ignore_index=True)
    all_y = {**y_rpe, **y_k562}

    shared_genes = sorted(list(set(rpe_genes) & set(k562_genes)))
    rpe_sh_idx = np.array([rpe_genes.index(g) for g in shared_genes])
    k562_sh_idx = np.array([k562_genes.index(g) for g in shared_genes])
    gene2idx = {
        "HCT116": {g: i for i, g in enumerate(rpe_genes)},
        "HEK293T": {g: i for i, g in enumerate(k562_genes)},
    }
    slot2line = {"HCT116": "RPE1", "HEK293T": "K562"}
    score_indices = {
        "HCT116_indices": np.arange(5000, dtype=np.int64),
        "HEK293T_indices": np.arange(5000, dtype=np.int64),
    }

    all_ens = sorted(set(train_df["perturbation"]) | set(rpe_genes) | set(k562_genes))
    mapped_ens = [g for g in all_ens if g in tss_map]
    print(f"Mapped {len(mapped_ens)}/{len(all_ens)} Replogle Ensembl genes to GRCh38 TSS ({100*len(mapped_ens)/len(all_ens):.1f}%)", flush=True)

    # Null distribution: 100,000 random distinct gene pairs from mapped Replogle universe
    rng_null = np.random.RandomState(42)
    idx_a = rng_null.randint(0, len(mapped_ens), size=100000)
    idx_b = rng_null.randint(0, len(mapped_ens), size=100000)
    mask_ne = idx_a != idx_b
    idx_a, idx_b = idx_a[mask_ne], idx_b[mask_ne]
    null_dists = np.array([tss_distance(mapped_ens[a], mapped_ens[b], tss_map) for a, b in zip(idx_a, idx_b)], dtype=np.float64)
    null_stats = {
        "n_pairs": int(len(null_dists)),
        "frac_lt_1kb": round(float(np.mean(null_dists < 1_000)), 6),
        "frac_lt_10kb": round(float(np.mean(null_dists < 10_000)), 6),
        "frac_lt_100kb": round(float(np.mean(null_dists < 100_000)), 6),
        "frac_same_chrom": round(float(np.mean(np.isfinite(null_dists))), 6),
    }

    # Run 5-fold strict-union OOF matching replogle_mechanism_reactome.py (seed 42)
    # and compare S5, original S6, S6_no_lt10kb, S6_no_lt100kb!
    model = PanPerturbBridgeModel(gene2idx, rpe_sh_idx, k562_sh_idx, score_indices)
    rng = np.random.RandomState(42)
    all_strong_in = sorted(list(
        set(df_rpe[(df_rpe["reproducible"] == 1) & (df_rpe["perturbation"].isin(gene2idx["HCT116"].keys()))]["perturbation"])
        | set(df_k562[(df_k562["reproducible"] == 1) & (df_k562["perturbation"].isin(gene2idx["HEK293T"].keys()))]["perturbation"])
    ))
    folds = np.array_split(rng.permutation(all_strong_in), 5)

    pair_records = []
    gene_records = []

    for fold in range(5):
        val_union = set(folds[fold])
        df_tr = train_df[~train_df["perturbation"].isin(val_union)].copy()
        bases, axes = {}, {}
        for slot in ["HCT116", "HEK293T"]:
            rep_s = df_tr[(df_tr["source"] == slot) & (df_tr["reproducible"] == 1)]
            bases[slot], axes[slot] = fit_baseline_and_axis(rep_s, all_y, score_indices[f"{slot}_indices"])

        for slot in ["HCT116", "HEK293T"]:
            oth = "HEK293T" if slot == "HCT116" else "HCT116"
            g2i = gene2idx[slot]
            s_idx = score_indices[f"{slot}_indices"]
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
            res = model.predict_for_source_all_stages(slot, df_own, df_oth, q_df, all_y, bases, axes, ("S5", "S6"))

            # Build own-line reciprocal matrices identically to PanPerturbBridgeModel
            own_rep = df_own[df_own["reproducible"] == 1].reset_index(drop=True)
            Y_raw = np.array([all_y[rid] for rid in own_rep["id"]], dtype=np.float32)
            Y_rep_p = project_rows_5000(Y_raw, bases[slot], axes[slot])
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

            # Evaluate own-line S6 counterfactuals: S5 + own_rec (all vs excluding <10kb vs excluding <100kb)
            rec_w = 0.55
            for _, r in q_df.iterrows():
                rid = r["id"]
                q_ens = r["perturbation"]
                q_sym = r["gene_symbol"]
                t_idx = g2i[q_ens]
                r_s5 = score_prediction(all_y[rid], res["S5"][rid], bases[slot], axes[slot], s_idx, t_idx)
                r_s6 = score_prediction(all_y[rid], res["S6"][rid], bases[slot], axes[slot], s_idx, t_idx)
                base_s5_proj = res["S5"][rid] - bases[slot]

                col_z = Z_z[:, t_idx]
                abs_z = np.abs(col_z)
                hits = np.where(abs_z > 4.0)[0]

                # Compute exact donor distances and weights
                dists = []
                weights = []
                for h in hits:
                    d_tss = tss_distance(q_ens, rep_perts[h], tss_map)
                    w_val = float((abs_z[h] - 4.0) ** 2)
                    dists.append(d_tss)
                    weights.append(w_val)
                    pair_records.append({
                        "cell_line": slot2line[slot],
                        "query_ens": q_ens,
                        "query_sym": q_sym,
                        "donor_ens": rep_perts[h],
                        "donor_sym": rep_syms[h],
                        "donor_z": round(float(col_z[h]), 4),
                        "weight": round(w_val, 4),
                        "tss_dist_bp": float(d_tss) if (d_tss is not None and np.isfinite(d_tss)) else (-1.0 if d_tss is None else 1e12),
                        "same_chrom": bool(d_tss is not None and np.isfinite(d_tss)),
                        "lt_1kb": bool(d_tss is not None and d_tss < 1_000),
                        "lt_10kb": bool(d_tss is not None and d_tss < 10_000),
                        "lt_100kb": bool(d_tss is not None and d_tss < 100_000),
                        "is_top_donor": bool(abs_z[h] == abs_z[hits].max()),
                        "gene_gain_S6_S5": round(r_s6 - r_s5, 5),
                    })

                def apply_own_rec(filter_min_bp=0.0):
                    out_v = base_s5_proj.copy()
                    if len(hits) == 0:
                        return score_prediction(all_y[rid], bases[slot] + out_v, bases[slot], axes[slot], s_idx, t_idx)
                    keep = []
                    for h, d_tss in zip(hits, dists):
                        if filter_min_bp > 0 and (d_tss is not None) and (d_tss < filter_min_bp):
                            continue
                        keep.append(h)
                    if not keep:
                        return r_s5
                    keep = np.array(keep, dtype=int)
                    wk = (abs_z[keep] - 4.0) ** 2
                    conf = float(np.sqrt(np.sum(wk)))
                    rec_vec = (wk @ Y_rep_p[keep]) / np.sum(wk)
                    no = np.linalg.norm(out_v[s_idx]) + 1e-8
                    nr = np.linalg.norm(rec_vec[s_idx]) + 1e-8
                    gate = np.clip(conf / 3.0, 0.2, 2.0)
                    out_v[s_idx] += (rec_w * gate) * (no / nr) * rec_vec[s_idx]
                    return score_prediction(all_y[rid], bases[slot] + out_v, bases[slot], axes[slot], s_idx, t_idx)

                r_own_all = apply_own_rec(0.0)
                r_own_no10k = apply_own_rec(10_000.0)
                r_own_no100k = apply_own_rec(100_000.0)

                w_tot = sum(weights) if weights else 0.0
                w_1k = sum(w for w, d in zip(weights, dists) if d is not None and d < 1_000)
                w_10k = sum(w for w, d in zip(weights, dists) if d is not None and d < 10_000)
                w_100k = sum(w for w, d in zip(weights, dists) if d is not None and d < 100_000)
                top_d = dists[int(np.argmax(weights))] if weights else None
                top_z = float(col_z[hits[int(np.argmax(weights))]]) if weights else 0.0

                gene_records.append({
                    "cell_line": slot2line[slot],
                    "query_ens": q_ens,
                    "query_sym": q_sym,
                    "r_S5": round(r_s5, 5),
                    "r_S6_full": round(r_s6, 5),
                    "r_S5_plus_ownRec_all": round(r_own_all, 5),
                    "r_S5_plus_ownRec_no10kb": round(r_own_no10k, 5),
                    "r_S5_plus_ownRec_no100kb": round(r_own_no100k, 5),
                    "gain_S6_S5": round(r_s6 - r_s5, 5),
                    "gain_ownRec_all": round(r_own_all - r_s5, 5),
                    "gain_ownRec_no10kb": round(r_own_no10k - r_s5, 5),
                    "gain_ownRec_no100kb": round(r_own_no100k - r_s5, 5),
                    "n_hits": int(len(hits)),
                    "frac_weight_lt1kb": round(w_1k / w_tot, 4) if w_tot > 0 else 0.0,
                    "frac_weight_lt10kb": round(w_10k / w_tot, 4) if w_tot > 0 else 0.0,
                    "frac_weight_lt100kb": round(w_100k / w_tot, 4) if w_tot > 0 else 0.0,
                    "top_donor_z": round(top_z, 3),
                    "top_donor_lt10kb": bool(top_d is not None and top_d < 10_000),
                    "top_donor_lt100kb": bool(top_d is not None and top_d < 100_000),
                })
        gc.collect()

    df_pairs = pd.DataFrame(pair_records)
    df_genes = pd.DataFrame(gene_records)
    df_genes.to_csv(REPORTS_DIR / "t6_replogle_crispri_proximity_genes.csv", index=False)

    # Compute comprehensive statistics
    p_valid = df_pairs[df_pairs["tss_dist_bp"] >= 0]
    p_neg = p_valid[p_valid["donor_z"] < -4.0]
    p_pos = p_valid[p_valid["donor_z"] > 4.0]
    p_top = p_valid[p_valid["is_top_donor"]]
    top40_pos = df_genes.sort_values("gain_ownRec_all", ascending=False).head(40)
    top40_s6 = df_genes.sort_values("gain_S6_S5", ascending=False).head(40)

    def pair_breakdown(sub):
        if len(sub) == 0:
            return {}
        w = sub["weight"].values
        return {
            "n_pairs": int(len(sub)),
            "frac_lt_1kb": round(float(sub["lt_1kb"].mean()), 5),
            "frac_lt_10kb": round(float(sub["lt_10kb"].mean()), 5),
            "frac_lt_100kb": round(float(sub["lt_100kb"].mean()), 5),
            "frac_same_chrom": round(float(sub["same_chrom"].mean()), 5),
            "weight_share_lt_1kb": round(float(w[sub["lt_1kb"].values].sum() / max(w.sum(), 1e-12)), 5),
            "weight_share_lt_10kb": round(float(w[sub["lt_10kb"].values].sum() / max(w.sum(), 1e-12)), 5),
            "weight_share_lt_100kb": round(float(w[sub["lt_100kb"].values].sum() / max(w.sum(), 1e-12)), 5),
        }

    summary = {
        "gtf_file": gtf_path.name,
        "gtf_sha256": gtf_sha,
        "mapped_genes": f"{len(mapped_ens)}/{len(all_ens)}",
        "null_random_pairs": null_stats,
        "all_active_pairs_abs_z_gt4": pair_breakdown(p_valid),
        "negative_z_pairs_z_lt_minus4": pair_breakdown(p_neg),
        "positive_z_pairs_z_gt_plus4": pair_breakdown(p_pos),
        "top_donor_per_query": pair_breakdown(p_top),
        "counterfactual_overall_all_2614": {
            "mean_S5": round(float(df_genes["r_S5"].mean()), 5),
            "mean_S5_plus_ownRec_all": round(float(df_genes["r_S5_plus_ownRec_all"].mean()), 5),
            "mean_S5_plus_ownRec_no10kb": round(float(df_genes["r_S5_plus_ownRec_no10kb"].mean()), 5),
            "mean_S5_plus_ownRec_no100kb": round(float(df_genes["r_S5_plus_ownRec_no100kb"].mean()), 5),
            "delta_all": round(float(df_genes["gain_ownRec_all"].mean()), 5),
            "delta_no10kb": round(float(df_genes["gain_ownRec_no10kb"].mean()), 5),
            "delta_no100kb": round(float(df_genes["gain_ownRec_no100kb"].mean()), 5),
        },
        "counterfactual_top40_ownRec_winners": {
            "n": 40,
            "frac_top_donor_lt10kb": round(float(top40_pos["top_donor_lt10kb"].mean()), 4),
            "frac_top_donor_lt100kb": round(float(top40_pos["top_donor_lt100kb"].mean()), 4),
            "mean_weight_share_lt10kb": round(float(top40_pos["frac_weight_lt10kb"].mean()), 4),
            "mean_gain_all": round(float(top40_pos["gain_ownRec_all"].mean()), 4),
            "mean_gain_no10kb": round(float(top40_pos["gain_ownRec_no10kb"].mean()), 4),
            "mean_gain_no100kb": round(float(top40_pos["gain_ownRec_no100kb"].mean()), 4),
            "pct_gain_retained_after_removing_10kb": round(100.0 * float(top40_pos["gain_ownRec_no10kb"].mean()) / max(float(top40_pos["gain_ownRec_all"].mean()), 1e-8), 2),
        },
        "counterfactual_top40_S6_S5_winners": {
            "n": 40,
            "frac_top_donor_lt10kb": round(float(top40_s6["top_donor_lt10kb"].mean()), 4),
            "frac_top_donor_lt100kb": round(float(top40_s6["top_donor_lt100kb"].mean()), 4),
            "frac_top_donor_neg_z": round(float((top40_s6["top_donor_z"] < 0).mean()), 4),
            "mean_gain_ownRec_all": round(float(top40_s6["gain_ownRec_all"].mean()), 4),
            "mean_gain_ownRec_no10kb": round(float(top40_s6["gain_ownRec_no10kb"].mean()), 4),
        },
        "runtime_seconds": round(time.time() - t0, 2),
    }

    with open(REPORTS_DIR / "t6_replogle_crispri_proximity_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
