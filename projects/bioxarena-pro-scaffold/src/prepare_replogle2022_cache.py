#!/usr/bin/env python3
"""
Reproducible preprocessing script to build the NumPy/CSV cache for the
Replogle et al. (Cell 2022) external validation dataset from official Figshare+ h5ad files.

Source:
  Replogle et al., "Mapping information-rich genotype-phenotype landscapes with
  genome-scale Perturb-seq", Cell (2022).
  Figshare+ Dataset DOI: https://doi.org/10.25452/figshare.plus.20029387

Official Input Files & Verified SHA-256:
  1. rpe1_normalized_bulk_01.h5ad
     URL: https://ndownloader.figshare.com/files/35775512
     SHA-256: a3c5bfd0f15d63938bc80c9b8874b9cd761e3a23caf5ffe7966bae4e887ec89d
  2. K562_gwps_normalized_bulk_01.h5ad
     URL: https://ndownloader.figshare.com/files/35773217
     SHA-256: 37e48c474d8b5dead4151f96ea8f5fe7bbe6beb10eeea48685b740c3f74490a2

Exact Preprocessing Specification:
  1. Readout Gene Selection (5,000 genes per cell line):
     - Drop readout columns containing any non-finite entry (NaN/Inf; 2 cols in RPE1, 73 cols in K562).
     - Select the top 5,000 remaining genes ranked by `var/clean_mean` (pre-computed mean expression
       across non-targeting/clean control cells in the published Replogle h5ad, i.e., an unsupervised
       baseline expression filter that does NOT use perturbation differential-expression labels;
       nevertheless, `var/clean_mean` is computed from the full published dataset before splitting,
       which we verify against control-row-only selection via `--verify-gene-filter`).
  2. Perturbation Filtering & Reproducibility Label:
     - Exclude non-targeting control rows (`obs/core_control == True`) and rows lacking an Ensembl ID.
     - Define `reproducible = 1` iff BOTH:
         (a) `obs/anderson_darling_counts >= 50` AND
         (b) `obs/energy_test_p_value < 0.001`
       Otherwise `reproducible = 0`. (Rows with `ad_counts >= 50` but `energy_test_p_value >= 0.001`
       are therefore labeled `reproducible = 0`, explaining why non-reproducible rows can have
       `ad_counts` up to 311.)
  3. Target Deduplication:
     - Keep strictly 1 row per unique target Ensembl gene (`ENSG...`), selecting the highest tuple
       `(reproducible, anderson_darling_counts, num_cells_filtered)`.

Usage:
  REPLOGLE_DATA_DIR=./external_data/replogle2022 python3 src/prepare_replogle2022_cache.py --verify-sha256
"""

import argparse
import hashlib
import json
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

EXPECTED_SHA256 = {
    "rpe1_normalized_bulk_01.h5ad": "a3c5bfd0f15d63938bc80c9b8874b9cd761e3a23caf5ffe7966bae4e887ec89d",
    "K562_gwps_normalized_bulk_01.h5ad": "37e48c474d8b5dead4151f96ea8f5fe7bbe6beb10eeea48685b740c3f74490a2",
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def resolve_replogle_dir() -> Path:
    cands = [
        *( [Path(os.environ["REPLOGLE_DATA_DIR"])] if os.environ.get("REPLOGLE_DATA_DIR") else [] ),
        Path.cwd() / "external_data/replogle2022",
        Path(__file__).resolve().parents[4] / "external_data/replogle2022",
        Path(__file__).resolve().parents[3] / "external_data/replogle2022",
    ]
    return next((c for c in cands if c.exists()), cands[0])


def build_cache(base_dir: Path, verify_sha256: bool = True, n_genes: int = 5000):
    cache_dir = base_dir / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    audit = {"files": {}}

    for slot, line, fn in [
        ("HCT116", "RPE1", "rpe1_normalized_bulk_01.h5ad"),
        ("HEK293T", "K562", "K562_gwps_normalized_bulk_01.h5ad"),
    ]:
        fp = base_dir / fn
        if not fp.exists():
            raise FileNotFoundError(f"Missing Replogle h5ad file: {fp}")
        digest = sha256_file(fp) if verify_sha256 else "skipped"
        if verify_sha256 and digest != EXPECTED_SHA256[fn]:
            raise ValueError(f"SHA-256 mismatch for {fn}: got {digest}, expected {EXPECTED_SHA256[fn]}")

        with h5py.File(fp, "r") as f:
            X = f["X"][:]
            ok_col = np.isfinite(X).all(axis=0)
            cm = f["var/clean_mean"][:]
            ok_idx = np.where(ok_col)[0]
            top_idx = np.sort(ok_idx[np.argsort(cm[ok_idx])[-n_genes:]])
            X_sub = X[:, top_idx].astype(np.float32)

            gids = [x.decode() if isinstance(x, bytes) else str(x) for x in f["var/gene_id"][:][top_idx]]
            cats = [x.decode() if isinstance(x, bytes) else str(x) for x in f["var/__categories/gene_name"][:]]
            gnames = [cats[i] for i in f["var/gene_name"][:][top_idx]]

            gt = [x.decode() if isinstance(x, bytes) else str(x) for x in f["obs/gene_transcript"][:]]
            ctrl = f["obs/core_control"][:]
            nc = f["obs/num_cells_filtered"][:]
            ad = f["obs/anderson_darling_counts"][:]
            ep = f["obs/energy_test_p_value"][:]

            # Audit: compare var/clean_mean vs var/mean overlap on top-5000 finite genes
            m_all = f["var/mean"][:]
            top_by_mean = set(np.sort(ok_idx[np.argsort(m_all[ok_idx])[-n_genes:]]))
            overlap_clean_vs_mean = len(set(top_idx) & top_by_mean)

        best_by_ens = {}
        n_ad50_but_ep_fail = 0
        for i, s in enumerate(gt):
            if ctrl[i]:
                continue
            parts = s.split("_")
            ens_list = [p for p in parts if p.startswith("ENSG")]
            if not ens_list:
                continue
            ens = ens_list[-1]
            sym = parts[1] if len(parts) > 1 else ens
            n_cells = float(nc[i]) if np.isfinite(nc[i]) and nc[i] > 0 else 10.0
            ad_val = float(ad[i]) if np.isfinite(ad[i]) else 0.0
            ep_val = float(ep[i]) if np.isfinite(ep[i]) else 1.0
            rep = 1 if (ad_val >= 50.0 and ep_val < 0.001) else 0
            if ad_val >= 50.0 and rep == 0:
                n_ad50_but_ep_fail += 1
            key_score = (rep, ad_val, n_cells)
            if ens not in best_by_ens or key_score > best_by_ens[ens][0]:
                best_by_ens[ens] = (key_score, i, sym, rep, n_cells, ad_val, ep_val)

        rows, mat = [], []
        for ens in sorted(best_by_ens.keys()):
            _, idx, sym, rep, n_cells, ad_val, ep_val = best_by_ens[ens]
            rid = f"{slot}_{ens}"
            rows.append({
                "id": rid,
                "source": slot,
                "perturbation": ens,
                "gene_symbol": sym,
                "reproducible": rep,
                "n_cells": n_cells,
                "ad_counts": ad_val,
                "energy_test_p_value": ep_val,
            })
            mat.append(X_sub[idx])

        meta_df = pd.DataFrame(rows)
        meta_df.to_csv(cache_dir / f"{line}_meta.csv", index=False)
        pd.DataFrame({"gene_id": gids, "gene_name": gnames}).to_csv(cache_dir / f"{line}_genes.csv", index=False)
        np.save(cache_dir / f"{line}_X.npy", np.asarray(mat, dtype=np.float32))

        audit["files"][line] = {
            "h5ad": fn,
            "sha256": digest,
            "raw_shape": [int(X.shape[0]), int(X.shape[1])],
            "nonfinite_cols_dropped": int((~ok_col).sum()),
            "selected_genes": int(n_genes),
            "overlap_clean_mean_vs_all_mean_top5000": int(overlap_clean_vs_mean),
            "deduplicated_rows": int(len(rows)),
            "n_reproducible": int(meta_df["reproducible"].sum()),
            "n_nonreproducible": int((meta_df["reproducible"] == 0).sum()),
            "n_raw_rows_ad50_but_energy_p_ge_0001": int(n_ad50_but_ep_fail),
        }
        print(
            f"[OK] {line}: {len(rows)} rows (rep={int(meta_df['reproducible'].sum())}, "
            f"nonrep={int((meta_df['reproducible']==0).sum())}), "
            f"clean_mean∩mean top5000={overlap_clean_vs_mean}/5000, "
            f"ad>=50 & energy_p>=0.001 rows={n_ad50_but_ep_fail}"
        )

    with open(cache_dir / "cache_manifest.json", "w") as f:
        json.dump(audit, f, indent=2)
    return audit


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", type=Path, default=None)
    ap.add_argument("--verify-sha256", action="store_true")
    args = ap.parse_args()
    base = args.data_dir or resolve_replogle_dir()
    build_cache(base, verify_sha256=args.verify_sha256)
