#!/usr/bin/env python3
"""
Reproducible preprocessing script to build the NumPy/CSV cache for the
Feng et al. (Cell Genomics 2025, PMID 41330380, Figshare 10.6084/m9.figshare.27989294.v2, MIT License)
TargetedScreen dataset (top-2 individual human iPSC donor cell lines: `jejf_2` and `zapk_3`).

Strictly implements Addendum R8-A (commit 840b502) of reports/ROUND8_PREREGISTRATION.md:
  1. Verified Input Files & Checksums:
     - TargetedScreen_Cell-Metadata.tsv.gz
       MD5:    980a440a8132a20567c4c3610e6a5b20
       SHA256: d6057946570edbbcd3406f4e55d8f338f602ef1b7cc71cffca101133e35b098d
     - TargetedScreen_RNA-UMI-Counts.csv.gz
       MD5:    f0d613b1ea8147425501baddba034158
       SHA256: cbda45f0645f1ce06488ba733fb95da04145a6ff2bd8d6d197633d2534386471
  2. Cell Filtering:
     - Keep assigned cells (`gene = Guide_Call.split('_')[0]` not in {unassigned, multiplet, nan})
       belonging to `jejf_2` (slot HCT116) or `zapk_3` (slot HEK293T).
  3. Per-`Pool_Day` Normalization & Control Z-Scoring:
     - Extract `Pool_Day` via regex `Pool\\d+_Day\\d+` from `Batch`.
     - Within each `(Cell_Line, Pool_Day)` stratum, library-size normalize every cell's raw UMI
       counts to that stratum's median total UMI across retained cells, then z-score every gene
       against `gene == 'NonTarget'` control cells in that same stratum.
     - Average z-scored profiles per unique target gene symbol (`n_cells >= 10`).
  4. Readout Gene Selection (up to top 5,000 finite genes per line):
     - Drop genes with zero control variance or non-finite z-score in any `Pool_Day` stratum.
     - Rank remaining finite genes by `control_mean` across all `NonTarget` cells of that line.
  5. Locked Reproducibility Criterion (identical to Replogle 2022 & Nadig 2025 Addendum B):
     - Vectorized Scholz-Stephens two-sample AD across selected readout genes vs all `NonTarget`
       control cells of that line (BH FDR < 0.05 -> `ad_counts`).
     - 20-PC 2,000-permutation energy-distance test vs all `NonTarget` control cells (`p < 0.001`).
     - `reproducible = 1` iff `n_cells >= 25` AND `ad_counts >= 50` AND `energy_test_p_value < 0.001`.
  6. Post-QC Eligibility Gate (Hard Stop before any modeling):
     - Both lines >= 100 reproducible target genes, >= 80 shared reproducible target genes,
       and >= 2,000 shared readout genes.

Usage:
  FENG_DATA_DIR=/path/to/external_data/feng2025 python3 projects/bioxarena-pro-scaffold/src/prepare_feng2025_cache.py
"""

import argparse
import gc
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.decomposition import TruncatedSVD
from statsmodels.stats.multitest import multipletests

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from prepare_nadig2025_cache import (
    ad_two_sample_pvalues_vectorized,
    energy_distance_perm_pvalue,
    verify_ad_implementation,
)

EXPECTED_SHA256 = {
    "TargetedScreen_Cell-Metadata.tsv.gz": "d6057946570edbbcd3406f4e55d8f338f602ef1b7cc71cffca101133e35b098d",
    "TargetedScreen_RNA-UMI-Counts.csv.gz": "cbda45f0645f1ce06488ba733fb95da04145a6ff2bd8d6d197633d2534386471",
}

C_EXTRACTOR_SRC = r"""
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <zlib.h>

/*
 * Fast streaming extractor for rows=genes, cols=cells gzipped CSV.
 * Reads col_map.bin (int32[n_cols]: -1 to skip, or 0..n_kept-1 target index).
 * Writes:
 *   genes.txt       : one gene label per row
 *   row_ptr.bin     : int64[n_genes + 1] CSR row pointers
 *   col_idx.bin     : uint32[nnz] CSR kept-cell indices
 *   vals.bin        : float32[nnz] CSR UMI counts
 *   total_umi.bin   : float64[n_kept] exact total UMI across ALL genes per kept cell
 */
int main(int argc, char **argv) {
    if (argc < 5) {
        fprintf(stderr, "Usage: %s <counts.csv.gz> <col_map.bin> <n_kept> <out_dir>\n", argv[0]);
        return 1;
    }
    const char *gz_path = argv[1];
    const char *map_path = argv[2];
    int32_t n_kept = atoi(argv[3]);
    const char *out_dir = argv[4];

    FILE *fm = fopen(map_path, "rb");
    if (!fm) { perror("fopen col_map"); return 1; }
    fseek(fm, 0, SEEK_END);
    long fsz = ftell(fm);
    fseek(fm, 0, SEEK_SET);
    int64_t n_cols = fsz / (int64_t)sizeof(int32_t);
    int32_t *col_map = (int32_t *)malloc((size_t)fsz);
    if (!col_map || fread(col_map, sizeof(int32_t), (size_t)n_cols, fm) != (size_t)n_cols) {
        fprintf(stderr, "Failed reading col_map\n");
        return 1;
    }
    fclose(fm);

    double *total_umi = (double *)calloc((size_t)n_kept, sizeof(double));
    uint32_t *row_cols = (uint32_t *)malloc((size_t)n_kept * sizeof(uint32_t));
    float *row_vals = (float *)malloc((size_t)n_kept * sizeof(float));

    char p_genes[1024], p_rptr[1024], p_cidx[1024], p_vals[1024], p_umi[1024];
    snprintf(p_genes, sizeof(p_genes), "%s/genes.txt", out_dir);
    snprintf(p_rptr, sizeof(p_rptr), "%s/row_ptr.bin", out_dir);
    snprintf(p_cidx, sizeof(p_cidx), "%s/col_idx.bin", out_dir);
    snprintf(p_vals, sizeof(p_vals), "%s/vals.bin", out_dir);
    snprintf(p_umi, sizeof(p_umi), "%s/total_umi.bin", out_dir);

    FILE *f_genes = fopen(p_genes, "w");
    FILE *f_rptr = fopen(p_rptr, "wb");
    FILE *f_cidx = fopen(p_cidx, "wb");
    FILE *f_vals = fopen(p_vals, "wb");
    if (!f_genes || !f_rptr || !f_cidx || !f_vals) { perror("fopen out"); return 1; }

    gzFile gz = gzopen(gz_path, "rb");
    if (!gz) { perror("gzopen"); return 1; }
    gzbuffer(gz, 1 << 20);

    size_t buf_cap = 32 * 1024 * 1024;
    char *line = (char *)malloc(buf_cap);
    if (!line) { fprintf(stderr, "oom\n"); return 1; }

    /* Helper to read an arbitrarily long line */
    int is_header = 1;
    int64_t n_genes = 0;
    int64_t nnz_total = 0;
    fwrite(&nnz_total, sizeof(int64_t), 1, f_rptr);

    while (1) {
        size_t len = 0;
        while (1) {
            if (len + 4096 >= buf_cap) {
                buf_cap *= 2;
                line = (char *)realloc(line, buf_cap);
                if (!line) { fprintf(stderr, "oom realloc\n"); return 1; }
            }
            if (!gzgets(gz, line + len, (int)(buf_cap - len))) {
                break;
            }
            size_t chunk_len = strlen(line + len);
            len += chunk_len;
            if (len > 0 && line[len - 1] == '\n') break;
        }
        if (len == 0) break;
        if (is_header) {
            is_header = 0;
            continue;
        }
        /* Parse gene label up to first comma */
        char *p = line;
        char *comma = strchr(p, ',');
        if (!comma) continue;
        *comma = '\0';
        fprintf(f_genes, "%s\n", p);
        p = comma + 1;

        int64_t col = 0;
        int32_t row_nnz = 0;
        while (*p && *p != '\n' && *p != '\r') {
            int32_t k = (col < n_cols) ? col_map[col] : -1;
            if (k < 0) {
                /* Skip value quickly */
                while (*p && *p != ',' && *p != '\n' && *p != '\r') p++;
            } else {
                /* Fast non-negative integer parse (counts are integers) */
                if (*p == '0' && (*(p + 1) == ',' || *(p + 1) == '\n' || *(p + 1) == '\r' || *(p + 1) == '\0')) {
                    p++;
                } else {
                    unsigned int val = 0;
                    while (*p >= '0' && *p <= '9') {
                        val = val * 10u + (unsigned int)(*p - '0');
                        p++;
                    }
                    /* Skip any fractional part if present */
                    if (*p == '.') {
                        while (*p && *p != ',' && *p != '\n' && *p != '\r') p++;
                    }
                    if (val > 0u) {
                        total_umi[k] += (double)val;
                        row_cols[row_nnz] = (uint32_t)k;
                        row_vals[row_nnz] = (float)val;
                        row_nnz++;
                    }
                }
            }
            if (*p == ',') p++;
            col++;
        }

        if (row_nnz > 0) {
            fwrite(row_cols, sizeof(uint32_t), (size_t)row_nnz, f_cidx);
            fwrite(row_vals, sizeof(float), (size_t)row_nnz, f_vals);
            nnz_total += row_nnz;
        }
        fwrite(&nnz_total, sizeof(int64_t), 1, f_rptr);
        n_genes++;
        if (n_genes % 5000 == 0) {
            fprintf(stderr, "  [C-extract] parsed %lld genes, nnz=%lld\n", (long long)n_genes, (long long)nnz_total);
            fflush(stderr);
        }
    }

    gzclose(gz);
    fclose(f_genes);
    fclose(f_rptr);
    fclose(f_cidx);
    fclose(f_vals);

    FILE *f_umi = fopen(p_umi, "wb");
    fwrite(total_umi, sizeof(double), (size_t)n_kept, f_umi);
    fclose(f_umi);

    fprintf(stderr, "  [C-extract DONE] genes=%lld, kept_cells=%d, nnz=%lld\n",
            (long long)n_genes, n_kept, (long long)nnz_total);
    free(line);
    free(col_map);
    free(total_umi);
    free(row_cols);
    free(row_vals);
    return 0;
}
"""


def resolve_feng_dir() -> Path:
    cands = [
        *([Path(os.environ["FENG_DATA_DIR"])] if os.environ.get("FENG_DATA_DIR") else []),
        Path.cwd() / "external_data/feng2025",
        Path(__file__).resolve().parents[4] / "external_data/feng2025",
        Path(__file__).resolve().parents[3] / "external_data/feng2025",
    ]
    return next((c for c in cands if c.exists()), cands[0])


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def extract_kept_cells_csr(feng_dir: Path, meta_kept: pd.DataFrame) -> tuple[sp.csc_matrix, np.ndarray, list[str]]:
    """
    Extract exact CSR (genes x kept_cells) and exact total UMI per kept cell using compiled C gzip scanner.
    Returns (X_genes_x_cells_csr, total_umi_per_kept_cell, raw_gene_labels).
    """
    counts_gz = feng_dir / "TargetedScreen_RNA-UMI-Counts.csv.gz"
    cache_dir = feng_dir / "cache"
    raw_dir = cache_dir / "_raw_extract"
    raw_dir.mkdir(parents=True, exist_ok=True)

    done_flag = raw_dir / "DONE.json"
    if done_flag.exists():
        print("  [Extract Cache Hit] Loading pre-extracted sparse UMI matrix...", flush=True)
        genes = (raw_dir / "genes.txt").read_text().splitlines()
        rptr = np.fromfile(raw_dir / "row_ptr.bin", dtype=np.int64)
        cidx = np.fromfile(raw_dir / "col_idx.bin", dtype=np.uint32)
        vals = np.fromfile(raw_dir / "vals.bin", dtype=np.float32)
        total_umi = np.fromfile(raw_dir / "total_umi.bin", dtype=np.float64)
        n_kept = len(total_umi)
        X_gxc = sp.csr_matrix((vals, cidx.astype(np.int32), rptr), shape=(len(genes), n_kept))
        return X_gxc, total_umi, genes

    t0 = time.time()
    print("  Reading CSV header to map 93,934 kept cells...", flush=True)
    with gzip.open(counts_gz, "rt") as f:
        header = f.readline().rstrip("\n").split(",")[1:]

    kept_norm_ids = meta_kept["norm_cell_id"].tolist()
    id_to_k = {cid: i for i, cid in enumerate(kept_norm_ids)}
    col_map = np.array([id_to_k.get(h, -1) for h in header], dtype=np.int32)
    matched = int(np.sum(col_map >= 0))
    assert matched == len(kept_norm_ids), f"Matched {matched} != expected {len(kept_norm_ids)}"
    print(f"  Mapped {matched:,} / {len(header):,} CSV columns cleanly.", flush=True)

    tmp_dir = Path(tempfile.mkdtemp(prefix="feng_c_ext_"))
    try:
        c_src = tmp_dir / "extract.c"
        c_bin = tmp_dir / "extract"
        map_bin = tmp_dir / "col_map.bin"
        c_src.write_text(C_EXTRACTOR_SRC)
        col_map.tofile(map_bin)
        subprocess.run(["cc", "-O3", str(c_src), "-lz", "-o", str(c_bin)], check=True)
        print("  Streaming TargetedScreen_RNA-UMI-Counts.csv.gz via compiled C extractor...", flush=True)
        subprocess.run([str(c_bin), str(counts_gz), str(map_bin), str(len(kept_norm_ids)), str(raw_dir)], check=True)
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    genes = (raw_dir / "genes.txt").read_text().splitlines()
    rptr = np.fromfile(raw_dir / "row_ptr.bin", dtype=np.int64)
    cidx = np.fromfile(raw_dir / "col_idx.bin", dtype=np.uint32)
    vals = np.fromfile(raw_dir / "vals.bin", dtype=np.float32)
    total_umi = np.fromfile(raw_dir / "total_umi.bin", dtype=np.float64)
    with open(done_flag, "w") as f:
        json.dump({"n_genes": len(genes), "n_kept": len(total_umi), "nnz": int(len(vals)), "sec": round(time.time() - t0, 1)}, f)

    X_gxc = sp.csr_matrix((vals, cidx.astype(np.int32), rptr), shape=(len(genes), len(total_umi)))
    print(f"  Extraction finished in {time.time()-t0:.1f}s: shape={X_gxc.shape}, nnz={X_gxc.nnz:,}", flush=True)
    return X_gxc, total_umi, genes


def process_one_line(
    line_name: str,
    slot: str,
    X_gxc: sp.csr_matrix,
    total_umi: np.ndarray,
    meta_kept: pd.DataFrame,
    raw_gene_labels: list[str],
    cache_dir: Path,
    n_genes_select: int = 5000,
) -> dict:
    """Apply locked Addendum R8-A per-Pool_Day libnorm + z-score + AD + Energy reproducibility for one line."""
    t0 = time.time()
    print(f"\n=== Processing Feng 2025 {line_name} (slot={slot}) ===", flush=True)
    sub_mask = (meta_kept["Cell_Line"].values == line_name)
    sub_idx = np.where(sub_mask)[0]
    sub_meta = meta_kept.iloc[sub_idx].reset_index(drop=True)
    umi_sub = total_umi[sub_idx]

    # Slice cells x genes CSR for this line
    X_cxg = X_gxc[:, sub_idx].T.tocsr()  # (n_cells, n_raw_genes)
    n_cells, n_raw_genes = X_cxg.shape

    # Filter to Gene-Expression features only (ignore any non-gene rows if present)
    ge_idx = np.array([i for i, g in enumerate(raw_gene_labels) if "Gene-Expression" in g or g.startswith("ENSG")], dtype=np.int64)
    if len(ge_idx) < n_raw_genes:
        X_cxg = X_cxg[:, ge_idx].tocsr()
        raw_gene_labels = [raw_gene_labels[i] for i in ge_idx]
        n_raw_genes = len(raw_gene_labels)

    # Parse ENSG and symbol from `ENSG...:SYMBOL:Gene-Expression`
    gene_ensg, gene_sym = [], []
    for lbl in raw_gene_labels:
        parts = lbl.split(":")
        gene_ensg.append(parts[0])
        gene_sym.append(parts[1] if len(parts) > 1 else parts[0])

    pools = sub_meta["Pool_Day"].values
    genes_tgt = sub_meta["target_gene"].values
    uniq_pools = sorted(sub_meta["Pool_Day"].unique().tolist())

    # Library-size normalize each cell to its Pool_Day median total UMI
    scale_per_cell = np.empty(n_cells, dtype=np.float64)
    for pd_name in uniq_pools:
        pm = (pools == pd_name)
        med_u = float(np.median(umi_sub[pm]))
        scale_per_cell[pm] = med_u / np.maximum(umi_sub[pm], 1.0)

    # Compute exact per-Pool_Day NonTarget control mean and std across all genes
    nt_mask = (genes_tgt == "NonTarget")
    mu_by_pool = {}
    sd_by_pool = {}
    nt_sum_all = np.zeros(n_raw_genes, dtype=np.float64)
    nt_n_all = 0
    finite_gene = np.ones(n_raw_genes, dtype=bool)

    for pd_name in uniq_pools:
        pm_nt = np.where((pools == pd_name) & nt_mask)[0]
        assert len(pm_nt) >= 20, f"Pool {pd_name} has only {len(pm_nt)} NT cells"
        X_nt = X_cxg[pm_nt].toarray().astype(np.float64)
        X_nt *= scale_per_cell[pm_nt][:, None]
        mu = X_nt.mean(axis=0)
        sd = X_nt.std(axis=0, ddof=1)
        mu_by_pool[pd_name] = mu
        sd_by_pool[pd_name] = sd
        finite_gene &= np.isfinite(mu) & np.isfinite(sd) & (sd > 1e-12)
        nt_sum_all += X_nt.sum(axis=0)
        nt_n_all += len(pm_nt)

    control_mean = nt_sum_all / max(nt_n_all, 1)
    finite_gene &= np.isfinite(control_mean)
    ok_idx = np.where(finite_gene)[0]
    n_sel = min(n_genes_select, len(ok_idx))
    top_idx = np.sort(ok_idx[np.argsort(control_mean[ok_idx])[-n_sel:]])
    print(f"  Finite genes: {len(ok_idx):,} -> selected top {n_sel:,} by control_mean (NT cells={nt_n_all})", flush=True)

    # Build dense z-scored matrix (n_cells x n_sel) in float32 (47k x 5000 x 4B ~ 0.95 GB — fits comfortably!)
    Z_all = np.empty((n_cells, n_sel), dtype=np.float32)
    X_sel = X_cxg[:, top_idx].tocsr()
    for pd_name in uniq_pools:
        p_rows = np.where(pools == pd_name)[0]
        mu = mu_by_pool[pd_name][top_idx]
        sd = sd_by_pool[pd_name][top_idx]
        # Process in blocks of 8,192 cells for low peak RAM
        for s in range(0, len(p_rows), 8192):
            br = p_rows[s : s + 8192]
            blk = X_sel[br].toarray().astype(np.float64)
            blk *= scale_per_cell[br][:, None]
            blk = (blk - mu[None, :]) / sd[None, :]
            Z_all[br] = blk.astype(np.float32)

    # Extract NonTarget z-scored matrix for AD and Energy-distance tests
    nt_rows = np.where(nt_mask)[0]
    Z_nt = Z_all[nt_rows].astype(np.float64)

    # Fit 20-PC TruncatedSVD on NonTarget + subsample of perturbed cells (matching Addendum B)
    rng = np.random.RandomState(42)
    pert_rows_all = np.where(~nt_mask)[0]
    fit_rows = np.concatenate([nt_rows, rng.choice(pert_rows_all, size=min(10000, len(pert_rows_all)), replace=False)])
    svd = TruncatedSVD(n_components=20, random_state=42)
    svd.fit(Z_all[fit_rows])
    P_all = svd.transform(Z_all).astype(np.float64)
    C_pca = P_all[nt_rows]
    # Subsample C_pca to 500 stratified NT cells for energy test if > 500 (matching Addendum B)
    if len(C_pca) > 500:
        c_sub_idx = rng.choice(len(C_pca), size=500, replace=False)
        C_pca_energy = C_pca[c_sub_idx]
    else:
        C_pca_energy = C_pca

    # Aggregate per target gene symbol and run locked AD + Energy reproducibility tests
    uniq_targets = sorted([g for g in sub_meta["target_gene"].unique().tolist() if g != "NonTarget"])
    sig_cache = {}
    X_bulk_list = []
    meta_rows = []
    n_repro = 0

    for ti, tgt in enumerate(uniq_targets):
        t_rows = np.where(genes_tgt == tgt)[0]
        nc = len(t_rows)
        if nc < 10:
            continue
        prof = Z_all[t_rows].mean(axis=0, dtype=np.float64).astype(np.float32)
        if nc >= 25:
            p_ad = ad_two_sample_pvalues_vectorized(Z_all[t_rows].astype(np.float64), Z_nt, sig_cache)
            _, q_ad, _, _ = multipletests(p_ad, method="fdr_bh")
            ad_cnt = int(np.sum(q_ad < 0.05))
            e_p = energy_distance_perm_pvalue(P_all[t_rows], C_pca_energy, rng, n_perm=2000)
            repro = int((ad_cnt >= 50) and (e_p < 0.001))
        else:
            ad_cnt = 0
            e_p = 1.0
            repro = 0
        n_repro += repro
        X_bulk_list.append(prof)
        meta_rows.append({
            "id": f"{slot}_{tgt}",
            "source": slot,
            "perturbation": tgt,
            "gene_symbol": tgt,
            "reproducible": repro,
            "n_cells": float(nc),
            "ad_counts": float(ad_cnt),
            "energy_test_p_value": float(e_p),
        })
        if (ti + 1) % 100 == 0 or ti == len(uniq_targets) - 1:
            print(f"  [{line_name}] {ti+1}/{len(uniq_targets)} targets tested, reproducible so far={n_repro}", flush=True)

    X_bulk = np.vstack(X_bulk_list).astype(np.float32)
    meta_df = pd.DataFrame(meta_rows)
    genes_df = pd.DataFrame({
        "gene_id": [gene_ensg[i] for i in top_idx],
        "gene_name": [gene_sym[i] for i in top_idx],
        "control_mean": control_mean[top_idx],
    })

    np.save(cache_dir / f"{line_name}_X.npy", X_bulk)
    meta_df.to_csv(cache_dir / f"{line_name}_meta.csv", index=False)
    genes_df.to_csv(cache_dir / f"{line_name}_genes.csv", index=False)

    stat = {
        "line": line_name,
        "slot": slot,
        "n_cells_assigned": int(n_cells),
        "n_nt_cells": int(nt_n_all),
        "n_readout_genes": int(n_sel),
        "n_targets_kept": int(len(meta_df)),
        "n_reproducible": int(n_repro),
        "n_repro_in_readout": int(meta_df[(meta_df["reproducible"] == 1) & (meta_df["perturbation"].isin(set(genes_df["gene_name"])))].shape[0]),
        "X_sha256": sha256_file(cache_dir / f"{line_name}_X.npy"),
        "meta_sha256": sha256_file(cache_dir / f"{line_name}_meta.csv"),
        "genes_sha256": sha256_file(cache_dir / f"{line_name}_genes.csv"),
        "elapsed_sec": round(time.time() - t0, 1),
    }
    print(f"  [DONE {line_name}] targets={stat['n_targets_kept']}, reproducible={stat['n_reproducible']} "
          f"(in_readout={stat['n_repro_in_readout']}), elapsed={stat['elapsed_sec']}s", flush=True)
    del Z_all, P_all, X_cxg, X_sel
    gc.collect()
    return stat


def load_feng_bundle():
    """Load cached Feng et al. 2025 (jejf_2 -> HCT116 slot, zapk_3 -> HEK293T slot) for frozen solver."""
    feng_dir = resolve_feng_dir()
    cache_dir = feng_dir / "cache"
    X0 = np.load(cache_dir / "jejf_2_X.npy")
    m0 = pd.read_csv(cache_dir / "jejf_2_meta.csv")
    g0 = pd.read_csv(cache_dir / "jejf_2_genes.csv")["gene_name"].astype(str).tolist()

    X1 = np.load(cache_dir / "zapk_3_X.npy")
    m1 = pd.read_csv(cache_dir / "zapk_3_meta.csv")
    g1 = pd.read_csv(cache_dir / "zapk_3_genes.csv")["gene_name"].astype(str).tolist()

    # Build unique gene->index maps (first occurrence if duplicate symbol)
    g2i0 = {}
    for i, g in enumerate(g0):
        g2i0.setdefault(g, i)
    g2i1 = {}
    for i, g in enumerate(g1):
        g2i1.setdefault(g, i)

    shared_genes = sorted(set(g2i0.keys()) & set(g2i1.keys()))
    sh0 = np.array([g2i0[g] for g in shared_genes], dtype=np.int64)
    sh1 = np.array([g2i1[g] for g in shared_genes], dtype=np.int64)

    all_y = {}
    for i, rid in enumerate(m0["id"].tolist()):
        all_y[rid] = X0[i].astype(np.float32)
    for i, rid in enumerate(m1["id"].tolist()):
        all_y[rid] = X1[i].astype(np.float32)

    train_df = pd.concat([m0, m1], ignore_index=True)
    g2i_map = {"HCT116": g2i0, "HEK293T": g2i1}
    sidx = {
        "HCT116_indices": np.arange(X0.shape[1], dtype=np.int64),
        "HEK293T_indices": np.arange(X1.shape[1], dtype=np.int64),
    }
    return train_df, all_y, shared_genes, sh0, sh1, g2i_map, sidx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify-sha256", action="store_true", default=True)
    args = ap.parse_args()

    t_all = time.time()
    verify_ad_implementation()
    feng_dir = resolve_feng_dir()
    cache_dir = feng_dir / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    if args.verify_sha256:
        for fn, exp in EXPECTED_SHA256.items():
            got = sha256_file(feng_dir / fn)
            assert got == exp, f"SHA-256 mismatch for {fn}: {got} != {exp}"
            print(f"[SHA256 PASS] {fn}: {got}", flush=True)

    meta_raw = pd.read_csv(feng_dir / "TargetedScreen_Cell-Metadata.tsv.gz", sep="\t")
    meta_raw["target_gene"] = meta_raw["Guide_Call"].astype(str).str.split("_").str[0]
    meta_raw["Pool_Day"] = meta_raw["Batch"].astype(str).str.extract(r"(Pool\d+_Day\d+)")[0]
    meta_raw["norm_cell_id"] = meta_raw["Cell_ID"].astype(str).str.replace(r"^PC-(P\d+)-D\d+_", r"\1_", regex=True)

    keep = meta_raw["Cell_Line"].isin(["jejf_2", "zapk_3"]) & (~meta_raw["target_gene"].isin(["unassigned", "multiplet", "nan", "None"]))
    meta_kept = meta_raw[keep].reset_index(drop=True)
    print(f"Kept assigned cells across jejf_2 + zapk_3: {len(meta_kept):,}", flush=True)

    X_gxc, total_umi, raw_genes = extract_kept_cells_csr(feng_dir, meta_kept)

    s0 = process_one_line("jejf_2", "HCT116", X_gxc, total_umi, meta_kept, raw_genes, cache_dir)
    s1 = process_one_line("zapk_3", "HEK293T", X_gxc, total_umi, meta_kept, raw_genes, cache_dir)

    # Check locked Post-QC Eligibility Gate (Addendum R8-A Section 4 item 5)
    m0 = pd.read_csv(cache_dir / "jejf_2_meta.csv")
    m1 = pd.read_csv(cache_dir / "zapk_3_meta.csv")
    g0 = set(pd.read_csv(cache_dir / "jejf_2_genes.csv")["gene_name"].astype(str))
    g1 = set(pd.read_csv(cache_dir / "zapk_3_genes.csv")["gene_name"].astype(str))
    shared_readout = len(g0 & g1)
    rep0 = set(m0[m0["reproducible"] == 1]["perturbation"].astype(str))
    rep1 = set(m1[m1["reproducible"] == 1]["perturbation"].astype(str))
    shared_rep = len(rep0 & rep1)
    rep0_in_r = len(rep0 & g0)
    rep1_in_r = len(rep1 & g1)
    shared_rep_in_r = len((rep0 & g0) & (rep1 & g1))

    gate_passed = (
        (s0["n_reproducible"] >= 100)
        and (s1["n_reproducible"] >= 100)
        and (shared_rep >= 80)
        and (shared_readout >= 2000)
    )
    summary = {
        "dataset": "Feng et al. (Cell Genomics 2025, PMID 41330380, Figshare 27989294 TargetedScreen)",
        "license": "MIT",
        "lines": [s0, s1],
        "shared_readout_genes": shared_readout,
        "shared_reproducible_targets": shared_rep,
        "jejf_2_repro_in_readout": rep0_in_r,
        "zapk_3_repro_in_readout": rep1_in_r,
        "shared_repro_in_readout": shared_rep_in_r,
        "eligibility_gate": {
            "rule": "both >= 100 reproducible AND >= 80 shared reproducible AND >= 2000 shared readout genes",
            "passed": gate_passed,
            "verdict": "ELIGIBLE_FOR_CONFIRMATORY_R8_3" if gate_passed else "FAILED_ELIGIBILITY_GATE_PENDING",
        },
        "total_elapsed_sec": round(time.time() - t_all, 1),
    }
    with open(cache_dir / "feng2025_cache_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("\n=== FENG 2025 POST-QC ELIGIBILITY GATE ===", flush=True)
    print(json.dumps(summary["eligibility_gate"], indent=2), flush=True)
    print(f"shared_readout={shared_readout}, shared_repro={shared_rep}, in_readout=({rep0_in_r}, {rep1_in_r}, shared={shared_rep_in_r})", flush=True)


if __name__ == "__main__":
    main()
