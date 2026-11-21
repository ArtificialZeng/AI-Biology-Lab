# Shared CRISPRi Cell-Line Inputs (`round9/p2/shared/`)

Unified single-cell control matrices and pseudobulk perturbation deltas prepared by **A-T1** using Addendum-B preprocessing.
See `DATA_SOURCES.md` for dataset URLs, licenses, byte sizes, and checksums.

**Common Aligned Readout Genes (`COMMON_GENES.txt`)**: **7501** HGNC gene symbols shared across all 2 lines.

## 1. Per-Line Summary & Gene-Symbol Alignment Loss

| Line | Raw Cells × Raw Genes | NT Ctrl (Sampled / Total) | Genes (Finite → Dedup Final) | Alignment Loss (Empty / Dup / Non-Finite) | Targets (Dedup / Reproducible) | Targets in `COMMON_GENES` (All / Rep) |
|---|---|---|---|---|---|---|
| `nadig_hepg2` | 145,473 × 9,624 | 4,976 / 4,976 | 9,433 → 9,432 | 0 / 1 / 191 | 2,343 / 1,000 | 2,052 / 951 |
| `nadig_jurkat` | 262,956 × 8,882 | 12,013 / 12,013 | 8,827 → 8,826 | 0 / 1 / 55 | 2,350 / 1,262 | 2,055 / 1,177 |

## 2. Pairwise Shared Target Counts (All Deduplicated Targets)

|              |   nadig_hepg2 |   nadig_jurkat |
|:-------------|--------------:|---------------:|
| nadig_hepg2  |          2343 |           2318 |
| nadig_jurkat |          2318 |           2350 |

## 3. Pairwise Shared Target Counts (Reproducible Targets: `n_cells >= 25`, `ad_counts >= 50`, `energy_p < 0.001`)

|              |   nadig_hepg2 |   nadig_jurkat |
|:-------------|--------------:|---------------:|
| nadig_hepg2  |          1000 |            762 |
| nadig_jurkat |           762 |           1262 |

## 4. Leave-One-Line-Out Shared Target Counts (Target Line vs. Union of Other Lines)

| Held-Out Line | Targets Shared w/ ≥1 Source Line (All) | Targets Shared w/ ≥1 Source Line (Reproducible) | Reproducible & in `COMMON_GENES` & ≥1 Source |
|---|---|---|---|
| `nadig_hepg2` | 2,318 | 762 | 721 |
| `nadig_jurkat` | 2,318 | 762 | 721 |

## 5. Saved File Schema (`p2/shared/lines/<line>/`)
- `genes.txt`: Aligned HGNC gene symbols matching column axis `0..G-1`.
- `targets.txt`: Deduplicated target HGNC gene symbols matching row axis `0..P-1`.
- `ctrl_cells.npz`: `X_log1p` (`log1p(1e4 * x / UMI)`), `X_zscore` (Addendum-B per-batch NT z-score), `X_libnorm` (`1e4 * x / UMI`), `ctrl_mean_log1p`, `ctrl_mean_libnorm`, `n_ctrl_total`, `batch_labels`.
- `pseudobulk_delta.npz`: `delta_log1p`, `delta_zscore`, `delta_libnorm` (all relative to per-batch NT control means), disjoint 50/50 split-halves (`half_a_delta_*`, `half_b_delta_*`, `half_valid`), `targets`, `target_ensembl_ids`, `reproducible`, `n_cells`, `ad_counts`, `energy_test_p_value`.
- `line_metadata.json`: Full provenance, SHA-256 of raw input `.h5ad`, and alignment/filtering loss counts.
