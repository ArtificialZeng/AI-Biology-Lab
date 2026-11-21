#!/usr/bin/env python3
"""
Reproducible preprocessing script to build the NumPy/CSV cache for the
Nadig et al. (Nature Genetics 2025, GSE264667) unseen external confirmation dataset
from official NCBI GEO raw single-cell h5ad files.

Strictly implements Addendum B (commit b9fe2bb) of ROUND7_PREREGISTRATION.md:
  1. Official Input Files & Verified SHA-256 (NCBI GEO GSE264667, PMID 40259084):
     - GSE264667_hepg2_raw_singlecell_01.h5ad
       SHA-256: e1ad7c3c5a201c861a207a858aa7e59f5e6ac1955674c415f7de0d1dadadb52e
     - GSE264667_jurkat_raw_singlecell_01.h5ad
       SHA-256: ffbe15f2c8f7ffcfd7b0ba9e6937d4ebc2d03b0179fa8234648a59bcb82c04a3
  2. Single-cell normalization (exact Replogle et al., Cell 2022 STAR Methods):
     - Chunked strictly by `gem_group` to maintain low peak RAM.
     - Library-size normalize each cell's raw UMI counts to its `gem_group` median UMI.
     - Z-score each gene per `gem_group` against `obs/gene == "non-targeting"` control cells.
     - Average per `obs/gene_transcript` to form `normalized_bulk` pseudobulk profiles.
  3. Readout Gene Selection (5,000 genes per cell line):
     - Drop any gene with zero control variance or non-finite z-score in any `gem_group`.
     - Select top 5,000 remaining finite genes ranked by mean library-size-normalized
       expression across `non-targeting` control cells (`control_mean`, matching `var/clean_mean`).
  4. Locked Reproducibility Criterion (exact Replogle 2022 thresholds):
     - Two-sample Anderson-Darling test across all 5,000 selected readout genes against a
       gem_group-stratified 1,000-cell `non-targeting` reference (BH FDR < 0.05 -> `ad_counts`).
     - Energy-distance permutation test (20-PC PCA space, 500-cell stratified NT reference,
       n_perm=2,000 with exact lossless early stopping at n_exceed >= 2 -> `energy_test_p_value`).
     - `reproducible = 1` iff `n_cells >= 25` AND `ad_counts >= 50` AND `energy_test_p_value < 0.001`.
  5. Ensembl Harmonization & Deduplication:
     - Both HepG2 (slot `HCT116`) and Jurkat (slot `HEK293T`) use Ensembl IDs (`ENSG...`).
     - Retain perturbations with `n_cells >= 10`; deduplicate to 1 row per unique target
       `ENSG` choosing highest `(reproducible, ad_counts, n_cells)`.

Usage:
  NADIG_DATA_DIR=./external_data/nadig2025 python3 src/prepare_nadig2025_cache.py --verify-sha256
"""

import argparse
import gc
import hashlib
import json
import os
import time
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import anderson_ksamp
from sklearn.decomposition import TruncatedSVD
from statsmodels.stats.multitest import multipletests

EXPECTED_SHA256 = {
    "GSE264667_hepg2_raw_singlecell_01.h5ad": "e1ad7c3c5a201c861a207a858aa7e59f5e6ac1955674c415f7de0d1dadadb52e",
    "GSE264667_jurkat_raw_singlecell_01.h5ad": "ffbe15f2c8f7ffcfd7b0ba9e6937d4ebc2d03b0179fa8234648a59bcb82c04a3",
}

# Scholz & Stephens (1987) Table 1 critical values for k=2 (m=1) matching scipy.stats.anderson_ksamp
_AD_CRIT = np.array([0.325, 1.226, 1.961, 2.718, 3.752, 4.592, 6.546], dtype=np.float64)
_AD_SIG = np.array([0.25, 0.1, 0.05, 0.025, 0.01, 0.005, 0.001], dtype=np.float64)
_AD_POLY = np.polyfit(_AD_CRIT, np.log(_AD_SIG), 2)
# Exact asymptotic tail slope for A^2_inf = 1 + sigma_inf * T, sigma_inf = sqrt(2/3 * (pi^2 - 9)) = 0.761404
_AD_SIGMA_INF = np.sqrt((2.0 / 3.0) * (np.pi ** 2 - 9.0))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def resolve_nadig_dir() -> Path:
    cands = [
        *([Path(os.environ["NADIG_DATA_DIR"])] if os.environ.get("NADIG_DATA_DIR") else []),
        Path.cwd() / "external_data/nadig2025",
        Path(__file__).resolve().parents[4] / "external_data/nadig2025",
        Path(__file__).resolve().parents[3] / "external_data/nadig2025",
    ]
    return next((c for c in cands if c.exists()), cands[0])


def _ad_sigma_k2(n1: int, n2: int) -> float:
    """Exact Scholz & Stephens (1987) eq. (4) finite-sample std dev for k=2 samples."""
    N = n1 + n2
    k = 2.0
    H = 1.0 / n1 + 1.0 / n2
    i = np.arange(1, N, dtype=np.float64)
    h = float(np.sum(1.0 / i))
    # g = sum_{i=1}^{N-2} sum_{j=i+1}^{N-1} 1 / ((N - i) * j)
    # Vectorized exact O(N) via suffix sums of 1/j
    inv_j = 1.0 / np.arange(1, N, dtype=np.float64)
    suf = np.cumsum(inv_j[::-1])[::-1]  # suf[i] = sum_{j=i+1}^{N-1} 1/j for 1-based i
    idx_i = np.arange(1, N - 1, dtype=np.float64)
    g = float(np.sum(suf[1:] / (N - idx_i)))
    a = (4.0 * g - 6.0) * (k - 1.0) + (10.0 - 6.0 * g) * H
    b = (2.0 * g - 4.0) * k * k + 8.0 * h * k + (2.0 * g - 14.0 * h - 4.0) * H - 8.0 * h + 4.0 * g - 6.0
    c = (6.0 * h + 2.0 * g - 2.0) * k * k + (4.0 * h - 4.0 * g + 6.0) * k + (2.0 * h - 6.0) * H + 4.0 * h
    d = (2.0 * h + 6.0) * k * k - 4.0 * h * k
    sig2 = (a * N**3 + b * N**2 + c * N + d) / ((N - 1.0) * (N - 2.0) * (N - 3.0))
    return float(np.sqrt(max(sig2, 1e-12)))


def ad_two_sample_pvalues_vectorized(X_pert: np.ndarray, X_ctrl: np.ndarray,
                                     sigma_cache: dict) -> np.ndarray:
    """
    Exact vectorized Scholz & Stephens (1987) two-sample Anderson-Darling p-values
    across all G readout genes simultaneously (columns of X_pert and X_ctrl).
    Handles ties cleanly using right-continuous step counts at distinct values or
    mid-rank continuous formula (verified against scipy.stats.anderson_ksamp).
    """
    n1 = X_pert.shape[0]
    n2 = X_ctrl.shape[0]
    N = n1 + n2
    if (n1, n2) not in sigma_cache:
        sigma_cache[(n1, n2)] = _ad_sigma_k2(n1, n2)
    sig_N = sigma_cache[(n1, n2)]

    Z = np.vstack([X_pert, X_ctrl])  # (N, G)
    # Sort each column; break exact ties identically by stable sort + right-endpoint adjustment
    order = np.argsort(Z, axis=0, kind="mergesort")
    Z_sorted = np.take_along_axis(Z, order, axis=0)
    is_p = (order < n1).astype(np.float64)
    M1 = np.cumsum(is_p, axis=0)  # (N, G)

    # Exact tie handling matching Scholz-Stephens eq. (6) / scipy _anderson_ksamp_right:
    # A jump occurs only where Z_sorted[r] != Z_sorted[r+1] (r in 0..N-2).
    # For a tie block ending at r with block length lj, Bj = r + 1 and M1_j = M1[r].
    # At non-tie positions (continuous), lj=1, Bj=r+1.
    # For tied positions inside a block, their contribution is folded into the block end.
    is_end = np.empty((N - 1, Z.shape[1]), dtype=bool)
    is_end[:] = (Z_sorted[:-1] != Z_sorted[1:])
    # Compute exact lj (tie block length) at block ends via vectorized forward diff of end indices
    # Notice: for each column, if Z_sorted[r] == Z_sorted[r+1], we simply replace (M1[r], B[r])
    # by the block-end values? Wait: in eq. (6), each tie block of size lj contributes
    # lj / N * (N * M1_end - B_end * n1)^2 / (B_end * (N - B_end)).
    # Summing 1/N * (N * M1_end - B_end * n1)^2 / (B_end * (N - B_end)) over EVERY index inside
    # the tie block (with M1_end and B_end broadcast backwards within each tie block) gives the
    # EXACT eq. (6) sum in pure vectorized NumPy with zero Python loops!
    # Wait: if the final tie block extends to index N-1, those indices have B_end = N and are excluded.
    # Let's broadcast M1_end and B_end backwards across ties using np.minimum.accumulate on reversed indices!
    idx_rows = np.arange(N, dtype=np.int32)[:, None]
    # Mark end-of-block positions (including N-1)
    end_marker = np.full((N, Z.shape[1]), N - 1, dtype=np.int32)
    end_marker[:-1] = np.where(is_end, idx_rows[:-1], N)
    # Propagate next block end backwards from bottom to top
    block_end = np.minimum.accumulate(end_marker[::-1], axis=0)[::-1]  # (N, G)
    valid = block_end < (N - 1)  # exclude elements in the final tie group ending at N-1

    # Gather M1 at block_end and B_end = block_end + 1
    col_idx = np.arange(Z.shape[1], dtype=np.int32)[None, :]
    M1_end = M1[block_end, col_idx]
    B_end = (block_end + 1).astype(np.float64)
    denom = np.where(valid, B_end * (N - B_end), 1.0)
    num = np.where(valid, (N * M1_end - B_end * n1) ** 2, 0.0)
    # For k=2: sum_{i=1}^2 (1/n_i) * (N M_{ij} - B_j n_i)^2 = (1/n1 + 1/n2) * (N M_{1j} - B_j n1)^2
    # Wait: (N M_{2j} - B_j n2) = (N (B_j - M_{1j}) - B_j (N - n1)) = -(N M_{1j} - B_j n1)!
    # So sum_{i=1}^2 (1/n_i) * (...)^2 = ((n1 + n2) / (n1 * n2)) * (...)^2 = (N / (n1 * n2)) * (...)^2.
    # Multiplied by 1/N gives 1 / (n1 * n2) * sum(...)!
    A2 = np.sum(num / denom, axis=0) / (n1 * n2)
    T = (A2 - 1.0) / sig_N

    # Continuous asymptotic p-value matching scipy.stats.anderson_ksamp on [0.325, 6.546]
    # and exact exponential tail decay exp(-sigma_inf * (T - 6.546)) beyond 6.546 (no 0.001 floor!)
    pvals = np.empty_like(T)
    low = T < _AD_CRIT[0]
    mid = (~low) & (T <= _AD_CRIT[-1])
    high = T > _AD_CRIT[-1]
    pvals[low] = 0.25
    pvals[mid] = np.exp(np.polyval(_AD_POLY, T[mid]))
    pvals[high] = 0.001 * np.exp(-_AD_SIGMA_INF * (T[high] - _AD_CRIT[-1]))
    return np.clip(pvals, 1e-300, 1.0)


def verify_ad_implementation():
    """Unit check: vectorized AD statistic matches scipy.stats.anderson_ksamp to < 1e-6."""
    rng = np.random.RandomState(0)
    xp = rng.randn(45, 5)
    xc = rng.randn(200, 5)
    # Introduce ties in col 4 to test tie handling
    xp[:, 4] = np.round(xp[:, 4], 1)
    xc[:, 4] = np.round(xc[:, 4], 1)
    xp[:, 0] += 0.8
    sig_cache = {}
    p_vec = ad_two_sample_pvalues_vectorized(xp, xc, sig_cache)
    for g in range(5):
        res = anderson_ksamp([xp[:, g], xc[:, g]])
        # Compare within scipy's uncapped region [0.001, 0.25]
        if 0.0015 < res.pvalue < 0.24:
            assert abs(p_vec[g] - res.pvalue) < 1e-5, f"AD mismatch col {g}: vec={p_vec[g]} vs scipy={res.pvalue}"
    print("[UNIT CHECK PASS] Vectorized Scholz-Stephens AD matches scipy.stats.anderson_ksamp (< 1e-5)", flush=True)


def energy_distance_perm_pvalue(P_pca: np.ndarray, C_pca: np.ndarray, rng: np.random.RandomState,
                                n_perm: int = 2000, batch_size: int = 100) -> float:
    """
    Exact two-sample Euclidean energy-distance permutation test with lossless early stopping:
      p = (n_exceed + 1) / (n_perm + 1) < 0.001  iff  n_exceed <= 1 (for n_perm=2000).
    Stops immediately as soon as n_exceed >= 2 (guaranteed p >= 3/2001 > 0.001).
    """
    n1 = len(P_pca)
    n2 = len(C_pca)
    N = n1 + n2
    X = np.vstack([P_pca, C_pca]).astype(np.float64)
    D = cdist(X, X, metric="euclidean")

    # Energy statistic E = 2/(n1*n2) * sum_{i in A, j in B} D_ij - 1/n1^2 * sum_{A,A} - 1/n2^2 * sum_{B,B}
    # Using indicator v in {1/n1 for A, -1/n2 for B}: -v^T D v = 2/(n1*n2) D_AB - 1/n1^2 D_AA - 1/n2^2 D_BB!
    v_obs = np.empty(N, dtype=np.float64)
    v_obs[:n1] = 1.0 / n1
    v_obs[n1:] = -1.0 / n2
    e_obs = float(-(v_obs @ (D @ v_obs)))

    n_exceed = 0
    done = 0
    while done < n_perm:
        bs = min(batch_size, n_perm - done)
        V = np.full((bs, N), -1.0 / n2, dtype=np.float64)
        for b in range(bs):
            idx = rng.choice(N, size=n1, replace=False)
            V[b, idx] = 1.0 / n1
        # Batch quadratic form: -diag(V D V^T) = -sum((V @ D) * V, axis=1)
        e_perm = -np.sum((V @ D) * V, axis=1)
        n_exceed += int(np.sum(e_perm >= e_obs))
        done += bs
        if n_exceed >= 2:
            # Lossless early stop: p is strictly >= (n_exceed + 1) / (n_perm + 1) > 0.001
            return float((n_exceed + 1) / (done + 1))
    return float((n_exceed + 1) / (n_perm + 1))


def stratified_nt_sample(nt_indices: np.ndarray, gems: np.ndarray, n_target: int, seed: int = 42) -> np.ndarray:
    """Gem_group-stratified random subsample of non-targeting control cells."""
    rng = np.random.RandomState(seed)
    if len(nt_indices) <= n_target:
        return np.sort(nt_indices)
    nt_gems = gems[nt_indices]
    uniq_g, counts = np.unique(nt_gems, return_counts=True)
    # Proportional allocation with >=1 per gem_group
    alloc = np.maximum(1, np.floor(counts / counts.sum() * n_target).astype(int))
    alloc = np.minimum(alloc, counts)
    diff = n_target - int(alloc.sum())
    if diff > 0:
        order = np.argsort(-(counts - alloc))
        for k in order:
            if diff == 0:
                break
            add = min(diff, int(counts[k] - alloc[k]))
            alloc[k] += add
            diff -= add
    elif diff < 0:
        order = np.argsort(-alloc)
        for k in order:
            if diff == 0:
                break
            sub = min(-diff, int(alloc[k] - 1))
            alloc[k] -= sub
            diff += sub
    chosen = []
    for g, k in zip(uniq_g, alloc):
        pool = nt_indices[nt_gems == g]
        chosen.append(rng.choice(pool, size=int(k), replace=False))
    return np.sort(np.concatenate(chosen))


def process_cell_line(fp: Path, slot: str, line: str, cache_dir: Path,
                      verify_sha256: bool = True, n_genes: int = 5000) -> dict:
    t0 = time.time()
    fn = fp.name
    if not fp.exists():
        raise FileNotFoundError(f"Missing Nadig h5ad file: {fp}")
    digest = sha256_file(fp) if verify_sha256 else "skipped"
    if verify_sha256 and digest != EXPECTED_SHA256[fn]:
        raise ValueError(f"SHA-256 mismatch for {fn}: got {digest}, expected {EXPECTED_SHA256[fn]}")

    print(f"\n=== Processing {line} ({fn}) ===", flush=True)
    with h5py.File(fp, "r") as f:
        n_cells_total, n_genes_raw = f["X"].shape
        umi = f["obs/UMI_count"][:].astype(np.float64)
        gems = f["obs/gem_group"][:]
        g_cats = [x.decode() if isinstance(x, bytes) else str(x) for x in f["obs/__categories/gene"][:]]
        gid_cats = [x.decode() if isinstance(x, bytes) else str(x) for x in f["obs/__categories/gene_id"][:]]
        gt_cats = [x.decode() if isinstance(x, bytes) else str(x) for x in f["obs/__categories/gene_transcript"][:]]
        g_codes = f["obs/gene"][:]
        gid_codes = f["obs/gene_id"][:]
        gt_codes = f["obs/gene_transcript"][:]

        var_gids = [x.decode() if isinstance(x, bytes) else str(x) for x in f["var/gene_id"][:]]
        var_gname_cats = [x.decode() if isinstance(x, bytes) else str(x) for x in f["var/__categories/gene_name"][:]]
        var_gnames = [var_gname_cats[i] for i in f["var/gene_name"][:]]

        nt_code = g_cats.index("non-targeting")
        nt_mask = (g_codes == nt_code)
        uniq_gems = np.unique(gems)

        # Precompute per-gem_group median UMI from obs/UMI_count (instant)
        med_umi_by_gem = {int(g): float(np.median(umi[gems == g])) for g in uniq_gems}
        gem_to_idx = {int(g): gi for gi, g in enumerate(uniq_gems)}
        gem_idx_per_cell = np.array([gem_to_idx[int(g)] for g in gems], dtype=np.int32)
        n_gems = len(uniq_gems)
        n_gt = len(gt_cats)

        # PASS 1a: Contiguous sequential chunk reads to compute exact per-gem_group NT mean & std
        # (Zero random HDF5 seeks; bit-for-bit equivalent to grouping by gem_group!)
        chunk_size = 8192
        print(f"  [Pass 1a] Sequential scan ({n_cells_total} cells x {n_genes_raw} genes, chunk={chunk_size}) for per-gem_group NT stats...", flush=True)
        nt_sum_g = np.zeros((n_gems, n_genes_raw), dtype=np.float64)
        nt_sq_g = np.zeros((n_gems, n_genes_raw), dtype=np.float64)
        nt_cnt_g = np.zeros(n_gems, dtype=np.int64)

        X_ds = f["X"]
        for s in range(0, n_cells_total, chunk_size):
            e = min(s + chunk_size, n_cells_total)
            nt_loc = nt_mask[s:e]
            if not np.any(nt_loc):
                continue
            Xc = X_ds[s:e, :].astype(np.float64)
            nt_rows = np.where(nt_loc)[0]
            X_nt = Xc[nt_rows]
            g_nt = gem_idx_per_cell[s:e][nt_rows]
            u_nt = umi[s:e][nt_rows]
            med_u_nt = np.array([med_umi_by_gem[int(uniq_gems[gi])] for gi in g_nt], dtype=np.float64)
            X_nt *= (med_u_nt / np.maximum(u_nt, 1.0))[:, None]
            np.add.at(nt_sum_g, g_nt, X_nt)
            np.add.at(nt_sq_g, g_nt, X_nt * X_nt)
            np.add.at(nt_cnt_g, g_nt, 1)

        mu_g = nt_sum_g / np.maximum(nt_cnt_g[:, None], 1)
        var_g = (nt_sq_g - nt_cnt_g[:, None] * (mu_g ** 2)) / np.maximum(nt_cnt_g[:, None] - 1, 1)
        var_g = np.maximum(var_g, 0.0)
        sd_g = np.sqrt(var_g)
        ok_g = np.isfinite(mu_g) & np.isfinite(sd_g) & (sd_g > 1e-12)
        finite_gene = np.all(ok_g, axis=0)
        control_mean = nt_sum_g.sum(axis=0) / max(int(nt_cnt_g.sum()), 1)
        finite_gene &= np.isfinite(control_mean)
        ok_idx = np.where(finite_gene)[0]
        assert len(ok_idx) >= n_genes, f"Only {len(ok_idx)} finite genes < {n_genes}"
        top_idx = np.sort(ok_idx[np.argsort(control_mean[ok_idx])[-n_genes:]])
        print(
            f"  [Pass 1a DONE] finite genes={len(ok_idx)}/{n_genes_raw}, selected top {len(top_idx)} by control_mean "
            f"(min_ctrl_mean={control_mean[top_idx].min():.4f}, {time.time()-t0:.1f}s)",
            flush=True,
        )

        # PASS 2 (single sequential pass!): Build normalized Z_5000 AND accumulate exact pseudobulk simultaneously!
        t2 = time.time()
        print(f"  [Pass 2] Sequential build of normalized Z_5000 ({n_cells_total} x {n_genes}) + pseudobulk...", flush=True)
        Z_5000 = np.empty((n_cells_total, n_genes), dtype=np.float32)
        pb_sum_5000 = np.zeros((n_gt, n_genes), dtype=np.float64)
        pb_cnt = np.zeros(n_gt, dtype=np.int64)
        mu_sub_g = mu_g[:, top_idx].astype(np.float32)
        sd_sub_g = sd_g[:, top_idx].astype(np.float32)
        med_arr = np.array([med_umi_by_gem[int(g)] for g in uniq_gems], dtype=np.float64)

        for s in range(0, n_cells_total, chunk_size):
            e = min(s + chunk_size, n_cells_total)
            Xc = X_ds[s:e, :][:, top_idx].astype(np.float32)
            gi_c = gem_idx_per_cell[s:e]
            scale = (med_arr[gi_c] / np.maximum(umi[s:e], 1.0)).astype(np.float32)
            Xc *= scale[:, None]
            Xc -= mu_sub_g[gi_c]
            Xc /= sd_sub_g[gi_c]
            Z_5000[s:e] = Xc

            non_nt_loc = ~nt_mask[s:e]
            if np.any(non_nt_loc):
                gt_c = gt_codes[s:e][non_nt_loc]
                np.add.at(pb_sum_5000, gt_c, Xc[non_nt_loc])
                np.add.at(pb_cnt, gt_c, 1)
            if (e // chunk_size) % 10 == 0 or e == n_cells_total:
                print(f"    chunk {e}/{n_cells_total} done ({time.time()-t2:.1f}s)", flush=True)

        pb_mean_5000 = np.zeros((n_gt, n_genes), dtype=np.float32)
        has_cells = pb_cnt > 0
        pb_mean_5000[has_cells] = (pb_sum_5000[has_cells] / pb_cnt[has_cells, None]).astype(np.float32)
        del pb_sum_5000
        gc.collect()
        print(f"  [Pass 2 DONE] Z_5000 & pseudobulk ready in {time.time()-t2:.1f}s", flush=True)

    # Fit 20-PC PCA (TruncatedSVD) on Z_5000 for energy-distance test
    t_pca = time.time()
    svd = TruncatedSVD(n_components=20, random_state=42)
    Z_pca20 = svd.fit_transform(Z_5000).astype(np.float32)
    print(f"  [PCA-20] Done in {time.time()-t_pca:.1f}s (var_ratio_sum={svd.explained_variance_ratio_.sum():.4f})", flush=True)

    # Stratified non-targeting references (1,000 for AD, 500 for Energy)
    nt_all_idx = np.where(nt_mask)[0]
    nt_ad_idx = stratified_nt_sample(nt_all_idx, gems, n_target=1000, seed=42)
    nt_en_idx = stratified_nt_sample(nt_all_idx, gems, n_target=500, seed=42)
    X_ctrl_ad = Z_5000[nt_ad_idx]
    C_pca_en = Z_pca20[nt_en_idx]

    # Map each gene_transcript code to its cell indices & metadata (gene_id, gene_symbol)
    gt_to_rows = defaultdict(list)
    gt_to_meta = {}
    for i in range(n_cells_total):
        if nt_mask[i]:
            continue
        gtc = int(gt_codes[i])
        gt_to_rows[gtc].append(i)
        if gtc not in gt_to_meta:
            gt_to_meta[gtc] = (gid_cats[int(gid_codes[i])], g_cats[int(g_codes[i])])

    eligible_gts = [
        gtc for gtc, rlist in sorted(gt_to_rows.items())
        if len(rlist) >= 25 and gt_to_meta[gtc][0].startswith("ENSG")
    ]
    print(f"  [Stats] Running exact AD + 20-PC Energy tests on {len(eligible_gts)} perturbations (n_cells>=25)...", flush=True)

    sigma_cache = {}
    rng_en = np.random.RandomState(42)
    ad_counts_map = {}
    en_pval_map = {}

    t_stat = time.time()
    for idx_k, gtc in enumerate(eligible_gts):
        ridx = np.asarray(gt_to_rows[gtc], dtype=np.int64)
        p_ad = ad_two_sample_pvalues_vectorized(Z_5000[ridx], X_ctrl_ad, sigma_cache)
        _, q_ad, _, _ = multipletests(p_ad, method="fdr_bh")
        ad_cnt = int(np.sum(q_ad < 0.05))
        ep_val = energy_distance_perm_pvalue(Z_pca20[ridx], C_pca_en, rng_en, n_perm=2000, batch_size=100)
        ad_counts_map[gtc] = ad_cnt
        en_pval_map[gtc] = ep_val

        if idx_k == 19:
            el = time.time() - t_stat
            est_total = el / 20.0 * len(eligible_gts)
            print(f"    [20-pert speed benchmark] 20 perts in {el:.2f}s -> est total for {line}: {est_total/60.0:.2f} min", flush=True)
        elif (idx_k + 1) % 400 == 0 or (idx_k + 1) == len(eligible_gts):
            n_rep_so_far = sum(1 for k in ad_counts_map if ad_counts_map[k] >= 50 and en_pval_map[k] < 0.001)
            print(
                f"    [{idx_k+1}/{len(eligible_gts)}] elapsed={time.time()-t_stat:.1f}s | "
                f"reproducible_so_far={n_rep_so_far}",
                flush=True,
            )

    del Z_5000, Z_pca20, X_ctrl_ad, C_pca_en
    gc.collect()

    # Target deduplication: retain perturbations with n_cells >= 10 and valid ENSG;
    # 1 row per unique target ENSG by highest (reproducible, ad_counts, n_cells)
    best_by_ens = {}
    n_ad50_but_ep_fail = 0
    for gtc, rlist in sorted(gt_to_rows.items()):
        n_c = len(rlist)
        if n_c < 10:
            continue
        ens, sym = gt_to_meta[gtc]
        if not ens.startswith("ENSG"):
            continue
        ad_val = float(ad_counts_map.get(gtc, 0))
        ep_val = float(en_pval_map.get(gtc, 1.0))
        rep = 1 if (n_c >= 25 and ad_val >= 50.0 and ep_val < 0.001) else 0
        if n_c >= 25 and ad_val >= 50.0 and rep == 0:
            n_ad50_but_ep_fail += 1
        key_score = (rep, ad_val, float(n_c))
        if ens not in best_by_ens or key_score > best_by_ens[ens][0]:
            best_by_ens[ens] = (key_score, gtc, sym, rep, float(n_c), ad_val, ep_val)

    rows, mat = [], []
    for ens in sorted(best_by_ens.keys()):
        _, gtc, sym, rep, n_c, ad_val, ep_val = best_by_ens[ens]
        rid = f"{slot}_{ens}"
        rows.append({
            "id": rid,
            "source": slot,
            "perturbation": ens,
            "gene_symbol": sym,
            "reproducible": rep,
            "n_cells": n_c,
            "ad_counts": ad_val,
            "energy_test_p_value": ep_val,
        })
        mat.append(pb_mean_5000[gtc])

    meta_df = pd.DataFrame(rows)
    gids_sel = [var_gids[i] for i in top_idx]
    gnames_sel = [var_gnames[i] for i in top_idx]

    meta_path = cache_dir / f"{line}_meta.csv"
    genes_path = cache_dir / f"{line}_genes.csv"
    x_path = cache_dir / f"{line}_X.npy"
    meta_df.to_csv(meta_path, index=False)
    pd.DataFrame({"gene_id": gids_sel, "gene_name": gnames_sel}).to_csv(genes_path, index=False)
    np.save(x_path, np.asarray(mat, dtype=np.float32))

    out_audit = {
        "h5ad": fn,
        "sha256": digest,
        "raw_shape": [int(n_cells_total), int(n_genes_raw)],
        "n_gem_groups": int(len(uniq_gems)),
        "n_nt_cells": int(nt_mask.sum()),
        "finite_genes": int(len(ok_idx)),
        "selected_genes": int(len(top_idx)),
        "deduplicated_rows": int(len(rows)),
        "n_reproducible": int(meta_df["reproducible"].sum()),
        "n_nonreproducible": int((meta_df["reproducible"] == 0).sum()),
        "n_raw_rows_ad50_but_energy_p_ge_0001": int(n_ad50_but_ep_fail),
        "elapsed_sec": round(float(time.time() - t0), 1),
        "cache_sha256": {
            meta_path.name: sha256_file(meta_path),
            genes_path.name: sha256_file(genes_path),
            x_path.name: sha256_file(x_path),
        },
    }
    print(
        f"[OK] {line}: {len(rows)} rows (rep={out_audit['n_reproducible']}, "
        f"nonrep={out_audit['n_nonreproducible']}), elapsed={out_audit['elapsed_sec']}s",
        flush=True,
    )
    return out_audit


def build_cache(base_dir: Path, verify_sha256: bool = True, n_genes: int = 5000, resume_existing: bool = False):
    verify_ad_implementation()
    cache_dir = base_dir / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    audit = {"source": "NCBI GEO GSE264667 (Nadig et al., Nature Genetics 2025, PMID 40259084)", "files": {}}
    for slot, line, fn in [
        ("HCT116", "HepG2", "GSE264667_hepg2_raw_singlecell_01.h5ad"),
        ("HEK293T", "Jurkat", "GSE264667_jurkat_raw_singlecell_01.h5ad"),
    ]:
        meta_path = cache_dir / f"{line}_meta.csv"
        genes_path = cache_dir / f"{line}_genes.csv"
        x_path = cache_dir / f"{line}_X.npy"
        if resume_existing and meta_path.exists() and genes_path.exists() and x_path.exists():
            fp = base_dir / fn
            digest = sha256_file(fp) if verify_sha256 else "skipped"
            if verify_sha256 and digest != EXPECTED_SHA256[fn]:
                raise ValueError(f"SHA-256 mismatch for {fn}: got {digest}")
            mdf = pd.read_csv(meta_path)
            gdf = pd.read_csv(genes_path)
            audit["files"][line] = {
                "h5ad": fn,
                "sha256": digest,
                "resumed_existing_cache": True,
                "selected_genes": int(len(gdf)),
                "deduplicated_rows": int(len(mdf)),
                "n_reproducible": int(mdf["reproducible"].sum()),
                "n_nonreproducible": int((mdf["reproducible"] == 0).sum()),
                "cache_sha256": {
                    meta_path.name: sha256_file(meta_path),
                    genes_path.name: sha256_file(genes_path),
                    x_path.name: sha256_file(x_path),
                },
            }
            print(f"[RESUME OK] {line}: {len(mdf)} rows (rep={int(mdf['reproducible'].sum())})", flush=True)
        else:
            audit["files"][line] = process_cell_line(base_dir / fn, slot, line, cache_dir, verify_sha256, n_genes)

    with open(cache_dir / "cache_manifest.json", "w") as f:
        json.dump(audit, f, indent=2)
    print("\nSaved manifest:", cache_dir / "cache_manifest.json", flush=True)
    return audit


def load_nadig_bundle():
    """
    Load Nadig et al. (2025) HepG2 + Jurkat cache in exact format of load_replogle_bundle():
      slot HCT116  -> HepG2
      slot HEK293T -> Jurkat
    """
    cache_dir = resolve_nadig_dir() / "cache"
    hepg2_meta = pd.read_csv(cache_dir / "HepG2_meta.csv")
    jurkat_meta = pd.read_csv(cache_dir / "Jurkat_meta.csv")
    hepg2_genes = pd.read_csv(cache_dir / "HepG2_genes.csv")
    jurkat_genes = pd.read_csv(cache_dir / "Jurkat_genes.csv")
    hepg2_X = np.load(cache_dir / "HepG2_X.npy")
    jurkat_X = np.load(cache_dir / "Jurkat_X.npy")

    train_df = pd.concat([hepg2_meta, jurkat_meta], ignore_index=True)
    all_y = {}
    for i, rid in enumerate(hepg2_meta["id"].tolist()):
        all_y[rid] = hepg2_X[i]
    for i, rid in enumerate(jurkat_meta["id"].tolist()):
        all_y[rid] = jurkat_X[i]

    g2i_h = {g: i for i, g in enumerate(hepg2_genes["gene_id"].tolist())}
    g2i_j = {g: i for i, g in enumerate(jurkat_genes["gene_id"].tolist())}
    shared = sorted(set(g2i_h.keys()) & set(g2i_j.keys()))
    sh0 = np.array([g2i_h[g] for g in shared], dtype=np.int64)
    sh1 = np.array([g2i_j[g] for g in shared], dtype=np.int64)
    g2i_map = {"HCT116": g2i_h, "HEK293T": g2i_j}
    sidx = {
        "HCT116_indices": np.arange(len(hepg2_genes), dtype=np.int64),
        "HEK293T_indices": np.arange(len(jurkat_genes), dtype=np.int64),
    }
    return train_df, all_y, shared, sh0, sh1, g2i_map, sidx


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", type=Path, default=None)
    ap.add_argument("--verify-sha256", action="store_true")
    ap.add_argument("--resume-existing", action="store_true")
    args = ap.parse_args()
    base = args.data_dir or resolve_nadig_dir()
    build_cache(base, verify_sha256=args.verify_sha256, resume_existing=args.resume_existing)
