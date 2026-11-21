#!/usr/bin/env python3
"""
Verify `p2/shared/lines/<line>/` integrity, compute `COMMON_GENES.txt`,
and generate `p2/shared/README.md` with per-line stats and pairwise target overlaps.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Set

import numpy as np
import pandas as pd

SHARED_DIR = Path(__file__).resolve().parents[1]
LINES_DIR = SHARED_DIR / "lines"


def verify_and_summarize() -> None:
    line_dirs = sorted(p for p in LINES_DIR.iterdir() if p.is_dir())
    if not line_dirs:
        raise RuntimeError(f"No line directories found in {LINES_DIR}")

    lines: List[str] = []
    gene_sets: Dict[str, Set[str]] = {}
    gene_lists: Dict[str, List[str]] = {}
    target_all: Dict[str, Set[str]] = {}
    target_rep: Dict[str, Set[str]] = {}
    metas: Dict[str, dict] = {}

    for ld in line_dirs:
        if not (ld / "line_metadata.json").exists():
            print(f"[SKIP incomplete/in-progress] {ld.name}")
            continue
        name = ld.name
        genes = [x.strip() for x in (ld / "genes.txt").read_text(encoding="utf-8").splitlines() if x.strip()]
        targets = [x.strip() for x in (ld / "targets.txt").read_text(encoding="utf-8").splitlines() if x.strip()]
        meta = json.loads((ld / "line_metadata.json").read_text(encoding="utf-8"))

        if len(genes) != len(set(genes)):
            raise ValueError(f"{name}: duplicate symbols in genes.txt!")
        if len(targets) != len(set(targets)):
            raise ValueError(f"{name}: duplicate symbols in targets.txt!")

        with np.load(ld / "ctrl_cells.npz") as zc:
            for k in ("X_log1p", "X_zscore", "X_libnorm"):
                arr = zc[k]
                if arr.shape[1] != len(genes):
                    raise ValueError(f"{name} {k} col mismatch: {arr.shape[1]} != {len(genes)}")
                if not np.isfinite(arr).all():
                    raise ValueError(f"{name} {k} has non-finite entries!")
            n_ctrl = int(zc["X_log1p"].shape[0])

        with np.load(ld / "pseudobulk_delta.npz") as zd:
            for k in (
                "delta_log1p", "delta_zscore", "delta_libnorm",
                "half_a_delta_log1p", "half_b_delta_log1p",
                "half_a_delta_zscore", "half_b_delta_zscore",
                "half_a_delta_libnorm", "half_b_delta_libnorm",
            ):
                arr = zd[k]
                if arr.shape != (len(targets), len(genes)):
                    raise ValueError(f"{name} {k} shape {arr.shape} != ({len(targets)}, {len(genes)})")
                if not np.isfinite(arr).all():
                    raise ValueError(f"{name} {k} has non-finite entries!")
            rep_mask = np.asarray(zd["reproducible"], dtype=bool)
            t_arr = [str(x) for x in zd["targets"].tolist()]
            if t_arr != targets:
                raise ValueError(f"{name}: targets.txt != pseudobulk_delta.npz['targets']")

        lines.append(name)
        gene_lists[name] = genes
        gene_sets[name] = set(genes)
        target_all[name] = set(targets)
        target_rep[name] = {t for t, r in zip(targets, rep_mask) if r}
        metas[name] = meta
        meta["n_ctrl_verified"] = n_ctrl

    common_genes = sorted(set.intersection(*(gene_sets[n] for n in lines)))
    common_set = set(common_genes)
    (SHARED_DIR / "COMMON_GENES.txt").write_text("\n".join(common_genes) + "\n", encoding="utf-8")

    # Build pairwise overlap matrices
    ov_all = pd.DataFrame(index=lines, columns=lines, dtype=int)
    ov_rep = pd.DataFrame(index=lines, columns=lines, dtype=int)
    for a in lines:
        for b in lines:
            ov_all.loc[a, b] = len(target_all[a] & target_all[b])
            ov_rep.loc[a, b] = len(target_rep[a] & target_rep[b])

    # Build README.md
    md = [
        "# Shared CRISPRi Cell-Line Inputs (`round9/p2/shared/`)",
        "",
        "Unified single-cell control matrices and pseudobulk perturbation deltas prepared by **A-T1** using Addendum-B preprocessing.",
        "See `DATA_SOURCES.md` for dataset URLs, licenses, byte sizes, and checksums.",
        "",
        f"**Common Aligned Readout Genes (`COMMON_GENES.txt`)**: **{len(common_genes)}** HGNC gene symbols shared across all {len(lines)} lines.",
        "",
        "## 1. Per-Line Summary & Gene-Symbol Alignment Loss",
        "",
        "| Line | Raw Cells × Raw Genes | NT Ctrl (Sampled / Total) | Genes (Finite → Dedup Final) | Alignment Loss (MultiFile / Empty / Dup / Non-Finite / Cap) | Targets (Dedup / Reproducible) | Targets in `COMMON_GENES` (All / Rep) |",
        "|---|---|---|---|---|---|---|",
    ]
    for name in lines:
        m = metas[name]
        raw_s = f"{m['raw_shape'][0]:,} × {m['raw_shape'][1]:,}"
        nt_s = f"{m['n_nt_sampled']:,} / {m['n_nt_total']:,}"
        g_s = f"{m['n_genes_finite_all_batches']:,} → {m['n_genes_final']:,}"
        non_fin = m["raw_shape"][1] - m["n_genes_finite_all_batches"]
        mf_drop = m.get("n_genes_multifile_dropped", 0)
        cap_drop = m.get("n_genes_low_expr_cap_dropped", max(0, m["n_genes_finite_all_batches"] - m["n_genes_empty_symbol_dropped"] - m["n_genes_duplicate_symbol_dropped"] - m["n_genes_final"]))
        loss_s = f"{mf_drop:,} / {m['n_genes_empty_symbol_dropped']} / {m['n_genes_duplicate_symbol_dropped']} / {non_fin:,} / {cap_drop:,}"
        t_s = f"{m['n_targets_deduplicated']:,} / {m['n_targets_reproducible']:,}"
        t_in_cg = len(target_all[name] & common_set)
        tr_in_cg = len(target_rep[name] & common_set)
        md.append(f"| `{name}` | {raw_s} | {nt_s} | {g_s} | {loss_s} | {t_s} | {t_in_cg:,} / {tr_in_cg:,} |")

    md.extend([
        "",
        "## 2. Pairwise Shared Target Counts (All Deduplicated Targets)",
        "",
        ov_all.to_markdown(),
        "",
        "## 3. Pairwise Shared Target Counts (Reproducible Targets: `n_cells >= 25`, `ad_counts >= 50`, `energy_p < 0.001`)",
        "",
        ov_rep.to_markdown(),
        "",
        "## 4. Leave-One-Line-Out Shared Target Counts (Target Line vs. Union of Other Lines)",
        "",
        "| Held-Out Line | Targets Shared w/ ≥1 Source Line (All) | Targets Shared w/ ≥1 Source Line (Reproducible) | Reproducible & in `COMMON_GENES` & ≥1 Source |",
        "|---|---|---|---|",
    ])
    for name in lines:
        others_all = set().union(*(target_all[o] for o in lines if o != name)) if len(lines) > 1 else set()
        others_rep = set().union(*(target_rep[o] for o in lines if o != name)) if len(lines) > 1 else set()
        n_sh_all = len(target_all[name] & others_all)
        n_sh_rep = len(target_rep[name] & others_rep)
        n_sh_rep_cg = len(target_rep[name] & others_rep & common_set)
        md.append(f"| `{name}` | {n_sh_all:,} | {n_sh_rep:,} | {n_sh_rep_cg:,} |")

    md.extend([
        "",
        "## 5. Notes, HGNC Alignment & Budget Disclosures",
        "- All datasets share 10x/CellRanger Ensembl-v98-compatible symbols (`497/503` legacy symbols in `vcc2025_h1` match `nourreddine_kolf21j`, `nadig_*`, and `replogle_*` identically without renaming); where target aliases broke self-gene lookup or cross-line matching (`nourreddine_kolf21j`: `KAT3A -> CREBBP` via `ENSG00000005339`, `MSK1 -> RPS6KA5` via `ENSG00000100784`), official HGNC symbols are restored via `hgnc_complete_set.txt`.",
        "- `nourreddine_kolf21j` combines the two CC BY 4.0 Strong Perturbation sub-atlases (`KOLF_Chromatin_Modifiers_Strong_Perturbations.h5ad` + `KOLF_Metabolic_Enzymes_Strong_Perturbations.h5ad`, `2.855 GB`), yielding `95` targets (`67` reproducible, `14` reproducible shared in LOLO), because `KOLF_Strong_Perturbations.h5ad` (`46.719 GB`) and `KOLF_Pan_Genome_QC_Filtered.h5ad` (`189.393 GB`) exceed the `60.00 GB` team download cap.",
        "- Disjoint 50/50 split-halves (`half_a_delta_*`, `half_b_delta_*`) split both non-targeting controls and perturbed cells independently per batch (strictly zero shared cells between Half A and Half B, as required for noise-ceiling estimation).",
        "",
        "## 6. Saved File Schema (`p2/shared/lines/<line>/`)",
        "- `genes.txt`: Aligned HGNC gene symbols matching column axis `0..G-1`.",
        "- `targets.txt`: Deduplicated target HGNC gene symbols matching row axis `0..P-1`.",
        "- `ctrl_cells.npz`: `X_log1p` (`log1p(1e4 * x / UMI)`), `X_zscore` (Addendum-B per-batch NT z-score), `X_libnorm` (`1e4 * x / UMI`), `ctrl_mean_log1p`, `ctrl_mean_libnorm`, `n_ctrl_total`, `batch_labels`.",
        "- `pseudobulk_delta.npz`: `delta_log1p`, `delta_zscore`, `delta_libnorm` (all relative to per-batch NT control means), disjoint 50/50 split-halves (`half_a_delta_*`, `half_b_delta_*`, `half_valid`), `targets`, `target_ensembl_ids`, `reproducible`, `n_cells`, `ad_counts`, `energy_test_p_value`.",
        "- `line_metadata.json`: Full provenance, SHA-256 of raw input `.h5ad`, and alignment/filtering loss counts.",
    ])

    (SHARED_DIR / "README.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"[VERIFIED] {len(lines)} lines, COMMON_GENES={len(common_genes)}")
    print(ov_all)


if __name__ == "__main__":
    verify_and_summarize()
