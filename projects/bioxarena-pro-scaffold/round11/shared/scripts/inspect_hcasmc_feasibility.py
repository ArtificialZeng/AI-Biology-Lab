#!/usr/bin/env python3
"""
PREREG_D §2.2.2 Q3 Mandatory Pre-Evaluation Feasibility Gate for HCASMC (`IGVFFI5999NXAF.h5`).

Counts & structure ONLY (zero target-gene correlations, residuals, slopes, or FCT values):
  1. Verifies Cell Ranger HDF5 `/matrix` layout and chunked column-sliced reading on 36 GiB RAM.
  2. Verifies guide-capture features (`feature_type != 'Gene Expression'`) & alignment with
     `IGVFFI4575UMXX.tsv.gz` (`non-targeting_`, `0Safe_`, targeting guides).
  3. Checks Parse sublibrary suffix recoverability in `/matrix/barcodes` and runs chunked RNA QC
     + single-sgRNA P4 counts across `batch_cols = 20,000` slices without loading full matrix into RAM.
"""

from __future__ import annotations

import argparse
import json
import resource
import sys
import time
from collections import Counter
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

PROJECT_ROOT = Path(__file__).resolve().parents[3]
R10_SCRIPTS = PROJECT_ROOT / "round10" / "shared" / "scripts"
if str(R10_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(R10_SCRIPTS))

from prepare_r10_lines import (  # noqa: E402
    extract_top2_csc,
    fit_poisson_gaussian_thresholds_from_hist,
)


def max_rss_gb() -> float:
    # On macOS, ru_maxrss is in bytes
    return float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1e9


def run_feasibility(
    h5_path: Path,
    guide_tsv_path: Path,
    out_json: Path,
    batch_cols: int = 20000,
    mt_max_pct: float = 10.0,
    min_genes: int = 500,
) -> dict:
    t0 = time.time()
    result = {
        "h5_path": str(h5_path),
        "feasible": False,
        "reason": "",
    }
    try:
        if not h5_path.exists():
            raise FileNotFoundError(f"Missing {h5_path}")
        gdf = pd.read_csv(guide_tsv_path, sep="\t", low_memory=False)
        ref_ntc = set(gdf.loc[gdf["type"] == "non-targeting", "guide_id"].astype(str))
        ref_safe = set(gdf.loc[gdf["type"] == "safe-targeting", "guide_id"].astype(str))
        ref_tgt = set(gdf.loc[gdf["type"] == "targeting", "guide_id"].astype(str))

        with h5py.File(h5_path, "r") as f:
            if "matrix" not in f:
                flat_req = ("cells", "genes", "indptr", "indices", "data", "shape")
                if all(k in f for k in flat_req):
                    sh = [int(x) for x in f["shape"][:]]
                    n_c, n_g = (sh[0], sh[1]) if len(sh) == 2 else (int(f["cells"].shape[0]), int(f["genes"].shape[0]))
                    genes_dec = [x.decode("utf-8", errors="replace") for x in f["genes"][:]]
                    cells_dec = [x.decode("utf-8", errors="replace") for x in f["cells"][:]]
                    us_cnt = Counter(b.split("__")[-1] if "__" in b else "" for b in cells_dec)
                    us_cnt.pop("", None)
                    n_ensg = sum(1 for g in genes_dec if g.startswith("ENSG"))
                    non_ensg = [g for g in genes_dec if not g.startswith("ENSG")]
                    all_guides = ref_ntc | ref_safe | ref_tgt
                    g_hits = sum(1 for g in genes_dec if g in all_guides)
                    result.update({
                        "h5_layout": "parse_flat_root_csr",
                        "shape": sh,
                        "nnz": int(f["data"].shape[0]),
                        "indptr_len": int(f["indptr"].shape[0]),
                        "indptr_dtype": str(f["indptr"].dtype),
                        "data_dtype": str(f["data"].dtype),
                        "data_chunks": str(f["data"].chunks),
                        "data_compression": str(f["data"].compression),
                        "n_cells": len(cells_dec),
                        "n_genes_entries": len(genes_dec),
                        "n_rna_ensg": int(n_ensg),
                        "non_ensg_sample": non_ensg[:20],
                        "n_guide_features": int(g_hits),
                        "guide_tsv_n_rows": int(len(gdf)),
                        "guide_tsv_types": {str(k): int(v) for k, v in gdf["type"].value_counts().items()},
                        "barcode_samples": cells_dec[:5] + cells_dec[-5:],
                        "barcode_us_suffix_counts": dict(us_cnt.most_common(20)),
                    })
                    raise RuntimeError(
                        f"Parse flat-root CSR matrix ({n_c:,} cells x {n_g:,} ENSG RNA features, nnz={int(f['data'].shape[0]):,}) "
                        f"contains ZERO CRISPR guide-capture features (guide_matches={g_hits}, non_ensg={non_ensg[:5]}); "
                        f"manifest guide TSVs ({len(gdf)} rows) are guide library tables lacking per-cell barcodes/assignments."
                    )
                raise RuntimeError(f"HDF5 missing /matrix group (root keys={list(f.keys())})")
            m = f["matrix"]
            for req in ("shape", "barcodes", "indptr", "indices", "data", "features"):
                if req not in m:
                    raise RuntimeError(f"HDF5 missing /matrix/{req}")

            n_feat, n_col = int(m["shape"][0]), int(m["shape"][1])
            indptr = m["indptr"][:]
            if len(indptr) != n_col + 1 or int(indptr[0]) != 0 or int(indptr[-1]) != int(m["data"].shape[0]):
                raise RuntimeError(
                    f"Malformed CSC indptr: len={len(indptr)} vs n_col+1={n_col+1}, "
                    f"end={indptr[-1]} vs data={m['data'].shape[0]}"
                )

            ftypes = np.array([x.decode("utf-8", errors="replace") for x in m["features/feature_type"][:]])
            fids = [x.decode("utf-8", errors="replace") for x in m["features/id"][:]]
            fnames = [x.decode("utf-8", errors="replace") for x in m["features/name"][:]]
            ftype_counts = {str(k): int(v) for k, v in zip(*np.unique(ftypes, return_counts=True))}

            rna_rows = np.where(ftypes == "Gene Expression")[0]
            crp_rows = np.where(ftypes != "Gene Expression")[0]
            result["shape"] = [n_feat, n_col]
            result["nnz"] = int(indptr[-1])
            result["indptr_dtype"] = str(indptr.dtype)
            result["data_dtype"] = str(m["data"].dtype)
            result["data_chunks"] = str(m["data"].chunks)
            result["data_compression"] = str(m["data"].compression)
            result["feature_type_counts"] = ftype_counts
            result["n_rna_features"] = int(len(rna_rows))
            result["n_non_gex_features"] = int(len(crp_rows))

            if len(rna_rows) < 3000:
                raise RuntimeError(f"Too few Gene Expression rows: {len(rna_rows)}")
            if len(crp_rows) == 0:
                raise RuntimeError("No guide-capture features (`feature_type != 'Gene Expression'`) found inside HDF5")

            # Match guide rows against IGVFFI4575UMXX by id or name
            crp_matched_key = []
            match_mode = None
            id_hits = sum((fids[i] in ref_ntc or fids[i] in ref_safe or fids[i] in ref_tgt) for i in crp_rows)
            nm_hits = sum((fnames[i] in ref_ntc or fnames[i] in ref_safe or fnames[i] in ref_tgt) for i in crp_rows)
            if nm_hits >= id_hits and nm_hits > 0:
                match_mode = "features/name"
                crp_matched_key = [fnames[i] for i in crp_rows]
            elif id_hits > 0:
                match_mode = "features/id"
                crp_matched_key = [fids[i] for i in crp_rows]
            else:
                raise RuntimeError(
                    f"Guide rows ({len(crp_rows)}) do not match IGVFFI4575UMXX guide_id "
                    f"(sample ids={[fids[i] for i in crp_rows[:5]]}, names={[fnames[i] for i in crp_rows[:5]]})"
                )

            n_ntc_ref = sum(k in ref_ntc for k in crp_matched_key)
            n_safe_ref = sum(k in ref_safe for k in crp_matched_key)
            n_tgt_ref = sum(k in ref_tgt for k in crp_matched_key)
            result["guide_match_mode"] = match_mode
            result["guide_feature_matches"] = {
                "ntc": int(n_ntc_ref),
                "safe_targeting": int(n_safe_ref),
                "targeting": int(n_tgt_ref),
            }
            if n_ntc_ref == 0 or n_tgt_ref == 0:
                raise RuntimeError(f"Missing NTC ({n_ntc_ref}) or targeting ({n_tgt_ref}) guide features")

            # Inspect barcode format & sublibrary suffixes cleanly
            bcs = [x.decode("utf-8", errors="replace") for x in m["barcodes"][:]]
            sample_bcs = bcs[:5] + bcs[-5:]
            # Parse sublibrary suffix candidates: `__<sublib>` or `-<sublib>` or `_s<N>`
            suffix_us = [b.split("__")[-1] if "__" in b else "" for b in bcs]
            us_counts = Counter(x for x in suffix_us if x)
            suffix_dash = [b.split("-")[-1] if "-" in b else "" for b in bcs]
            dash_counts = Counter(x for x in suffix_dash if x)
            result["barcode_samples"] = sample_bcs
            result["barcode_us_suffix_counts"] = dict(us_counts.most_common(20))
            result["barcode_dash_suffix_counts"] = dict(dash_counts.most_common(20))

            # Chunked streaming RNA QC + CRISPR sparse capture extraction
            mt_raw = np.zeros(n_feat, dtype=bool)
            mt_raw[rna_rows] = [fnames[i].startswith("MT-") for i in rna_rows]
            is_crp = np.zeros(n_feat, dtype=bool)
            is_crp[crp_rows] = True
            crp_remap_idx = np.full(n_feat, -1, dtype=np.int32)
            crp_remap_idx[crp_rows] = np.arange(len(crp_rows), dtype=np.int32)
            is_rna_feat = ~is_crp

            umi = np.zeros(n_col, dtype=np.float64)
            ngenes = np.zeros(n_col, dtype=np.int32)
            mt_umi = np.zeros(n_col, dtype=np.float64)
            hist = np.zeros((len(crp_rows), 66), dtype=np.int64)
            c_idx_list, c_dat_list, c_col_list = [], [], []

            data_ds = m["data"]
            indices_ds = m["indices"]
            for cs in range(0, n_col, batch_cols):
                ce = min(cs + batch_cols, n_col)
                p0, p1 = int(indptr[cs]), int(indptr[ce])
                if p1 == p0:
                    continue
                dat = data_ds[p0:p1]
                idx = indices_ds[p0:p1]
                col_lens = np.diff(indptr[cs : ce + 1])
                col_ids = np.repeat(np.arange(cs, ce, dtype=np.int32), col_lens)

                m_r = is_rna_feat[idx]
                if np.any(m_r):
                    r_dat = dat[m_r].astype(np.float64, copy=False)
                    r_idx = idx[m_r]
                    r_col = col_ids[m_r] - cs
                    umi[cs:ce] = np.bincount(r_col, weights=r_dat, minlength=ce - cs)
                    ngenes[cs:ce] = np.bincount(r_col, minlength=ce - cs)
                    mt_m = mt_raw[r_idx]
                    if np.any(mt_m):
                        mt_umi[cs:ce] = np.bincount(r_col[mt_m], weights=r_dat[mt_m], minlength=ce - cs)

                m_c = ~m_r
                if np.any(m_c):
                    c_idx_list.append(crp_remap_idx[idx[m_c]])
                    c_dat_list.append(dat[m_c].astype(np.int64, copy=False))
                    c_col_list.append(col_ids[m_c])

        pct_mt = 100.0 * mt_umi / np.maximum(umi, 1.0)
        ok_init = (pct_mt < 20.0) & (ngenes >= 200)
        ok_hq = (pct_mt <= mt_max_pct) & (ngenes >= min_genes)

        c_idx = np.concatenate(c_idx_list) if c_idx_list else np.empty(0, dtype=np.int32)
        c_dat = np.concatenate(c_dat_list) if c_dat_list else np.empty(0, dtype=np.int64)
        c_col = np.concatenate(c_col_list) if c_col_list else np.empty(0, dtype=np.int32)
        del c_idx_list, c_dat_list, c_col_list

        ok_init_mask = ok_init[c_col]
        bins = np.clip(c_dat[ok_init_mask], 1, 65)
        np.add.at(hist, (c_idx[ok_init_mask], bins), 1)
        thrs = fit_poisson_gaussian_thresholds_from_hist(hist, min_thr=5, max_thr=20, default_thr=7)

        hq_cols = np.where(ok_hq)[0]
        col_map_hq = np.full(n_col, -1, dtype=np.int32)
        col_map_hq[hq_cols] = np.arange(len(hq_cols), dtype=np.int32)
        hq_m = col_map_hq[c_col] >= 0
        C_hq = sp.csr_matrix(
            (c_dat[hq_m].astype(np.float64), (c_idx[hq_m], col_map_hq[c_col[hq_m]])),
            shape=(len(crp_rows), len(hq_cols)),
        ).tocsc()
        t1, t2, a1 = extract_top2_csc(C_hq)
        pos_mask = C_hq.data >= thrs[C_hq.indices]
        pos_csc = sp.csc_matrix((pos_mask.astype(np.int8), C_hq.indices, C_hq.indptr), shape=C_hq.shape)
        n_pos = np.asarray(pos_csc.sum(axis=0)).ravel()
        s1_mask = (n_pos == 1) & (t1 >= 5.0) & ((t2 < 3.0) | (t1 / np.maximum(t2, 0.5) >= 3.0))

        # Check also simple floor>=3 / 4 vs Poisson-Gaussian
        n_ntc = n_safe = n_tgt = 0
        for loc_j in np.where(s1_mask)[0]:
            gk = crp_matched_key[int(a1[loc_j])]
            if gk in ref_ntc:
                n_ntc += 1
            elif gk in ref_safe:
                n_safe += 1
            elif gk in ref_tgt:
                n_tgt += 1

        result.update({
            "n_hq_cells": int(len(hq_cols)),
            "median_guide_threshold": float(np.median(thrs)),
            "n_single_sgrna_ntc": int(n_ntc),
            "n_single_sgrna_safe_excluded": int(n_safe),
            "n_single_sgrna_targeting": int(n_tgt),
            "max_rss_gb": round(max_rss_gb(), 3),
            "elapsed_sec": round(time.time() - t0, 2),
        })
        if n_ntc == 0 or n_tgt == 0:
            raise RuntimeError(f"Insufficient P4 single-sgRNA counts: ntc={n_ntc}, tgt={n_tgt}")

        result["feasible"] = True
        result["reason"] = "Chunked HDF5 read + guide-capture extraction succeeded within 36 GiB RAM"
    except Exception as e:
        result["feasible"] = False
        result["reason"] = f"{type(e).__name__}: {e}"
        result["max_rss_gb"] = round(max_rss_gb(), 3)
        result["elapsed_sec"] = round(time.time() - t0, 2)

    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return result


def main() -> None:
    p = argparse.ArgumentParser(description="HCASMC Q3 pre-evaluation HDF5 feasibility checker (counts & structure only).")
    p.add_argument("--h5", type=Path, default=PROJECT_ROOT / "round11" / "shared" / "raw" / "hcasmc" / "IGVFFI5999NXAF.h5")
    p.add_argument("--guide-tsv", type=Path, default=PROJECT_ROOT / "round11" / "shared" / "raw" / "hcasmc" / "IGVFFI4575UMXX.tsv.gz")
    p.add_argument("--out-json", type=Path, default=PROJECT_ROOT / "round11" / "shared" / "hcasmc_feasibility.json")
    p.add_argument("--mt-max-pct", type=float, default=10.0)
    p.add_argument("--min-genes", type=int, default=500)
    args = p.parse_args()
    res = run_feasibility(args.h5, args.guide_tsv, args.out_json, mt_max_pct=args.mt_max_pct, min_genes=args.min_genes)
    if not res["feasible"]:
        print("[Q3 FEASIBILITY EXCLUSION TRIGGERED]", res["reason"], file=sys.stderr)


if __name__ == "__main__":
    main()
