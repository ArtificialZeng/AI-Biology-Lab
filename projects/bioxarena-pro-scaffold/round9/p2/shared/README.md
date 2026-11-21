# Shared CRISPRi Cell-Line Inputs (`round9/p2/shared/`)

Unified single-cell control matrices and pseudobulk perturbation deltas prepared by **A-T1** using Addendum-B preprocessing.
See `DATA_SOURCES.md` for dataset URLs, licenses, byte sizes, and checksums.

**Common Aligned Readout Genes (`COMMON_GENES.txt`)**: **6584** HGNC gene symbols shared across all 5 lines.

## 1. Per-Line Summary & Gene-Symbol Alignment Loss

| Line | Raw Cells × Raw Genes | NT Ctrl (Sampled / Total) | Genes (Finite → Dedup Final) | Alignment Loss (MultiFile / Empty / Dup / Non-Finite / Cap) | Targets (Dedup / Reproducible) | Targets in `COMMON_GENES` (All / Rep) |
|---|---|---|---|---|---|---|
| `nadig_hepg2` | 145,473 × 9,624 | 4,976 / 4,976 | 9,433 → 9,432 | 0 / 0 / 1 / 191 / 0 | 2,343 / 1,000 | 1,989 / 932 |
| `nadig_jurkat` | 262,956 × 8,882 | 12,013 / 12,013 | 8,827 → 8,826 | 0 / 0 / 1 / 55 / 0 | 2,350 / 1,262 | 1,992 / 1,147 |
| `nourreddine_kolf21j` | 38,322 × 34,583 | 17,895 / 17,895 | 24,116 → 15,000 | 484 / 0 / 0 / 10,467 / 9,116 | 95 / 67 | 71 / 51 |
| `replogle_k562_essential` | 310,385 × 8,563 | 10,691 / 10,691 | 8,563 → 8,561 | 0 / 0 / 2 / 0 / 0 | 2,049 / 1,113 | 1,762 / 1,054 |
| `replogle_rpe1` | 247,914 × 8,749 | 11,485 / 11,485 | 8,747 → 8,746 | 0 / 0 / 1 / 2 / 0 | 2,351 / 1,575 | 1,994 / 1,414 |

## 2. Pairwise Shared Target Counts (All Deduplicated Targets)

|                         |   nadig_hepg2 |   nadig_jurkat |   nourreddine_kolf21j |   replogle_k562_essential |   replogle_rpe1 |
|:------------------------|--------------:|---------------:|----------------------:|--------------------------:|----------------:|
| nadig_hepg2             |          2343 |           2318 |                    25 |                      2005 |            2316 |
| nadig_jurkat            |          2318 |           2350 |                    25 |                      2011 |            2322 |
| nourreddine_kolf21j     |            25 |             25 |                    95 |                        17 |              24 |
| replogle_k562_essential |          2005 |           2011 |                    17 |                      2049 |            2017 |
| replogle_rpe1           |          2316 |           2322 |                    24 |                      2017 |            2351 |

## 3. Pairwise Shared Target Counts (Reproducible Targets: `n_cells >= 25`, `ad_counts >= 50`, `energy_p < 0.001`)

|                         |   nadig_hepg2 |   nadig_jurkat |   nourreddine_kolf21j |   replogle_k562_essential |   replogle_rpe1 |
|:------------------------|--------------:|---------------:|----------------------:|--------------------------:|----------------:|
| nadig_hepg2             |          1000 |            762 |                     8 |                       750 |             845 |
| nadig_jurkat            |           762 |           1262 |                    10 |                       791 |             995 |
| nourreddine_kolf21j     |             8 |             10 |                    67 |                         9 |              13 |
| replogle_k562_essential |           750 |            791 |                     9 |                      1113 |             872 |
| replogle_rpe1           |           845 |            995 |                    13 |                       872 |            1575 |

## 4. Leave-One-Line-Out Shared Target Counts (Target Line vs. Union of Other Lines)

| Held-Out Line | Targets Shared w/ ≥1 Source Line (All) | Targets Shared w/ ≥1 Source Line (Reproducible) | Reproducible & in `COMMON_GENES` & ≥1 Source |
|---|---|---|---|
| `nadig_hepg2` | 2,343 | 984 | 918 |
| `nadig_jurkat` | 2,350 | 1,178 | 1,073 |
| `nourreddine_kolf21j` | 25 | 14 | 13 |
| `replogle_k562_essential` | 2,039 | 1,041 | 989 |
| `replogle_rpe1` | 2,350 | 1,287 | 1,172 |

## 5. Notes, HGNC Alignment & Budget Disclosures
- All datasets share 10x/CellRanger Ensembl-v98-compatible symbols (`497/503` legacy symbols in `vcc2025_h1` match `nourreddine_kolf21j`, `nadig_*`, and `replogle_*` identically without renaming); where target aliases broke self-gene lookup or cross-line matching (`nourreddine_kolf21j`: `KAT3A -> CREBBP` via `ENSG00000005339`, `MSK1 -> RPS6KA5` via `ENSG00000100784`), official HGNC symbols are restored via `hgnc_complete_set.txt`.
- `nourreddine_kolf21j` combines the two CC BY 4.0 Strong Perturbation sub-atlases (`KOLF_Chromatin_Modifiers_Strong_Perturbations.h5ad` + `KOLF_Metabolic_Enzymes_Strong_Perturbations.h5ad`, `2.855 GB`), yielding `95` targets (`67` reproducible, `14` reproducible shared in LOLO), because `KOLF_Strong_Perturbations.h5ad` (`46.719 GB`) and `KOLF_Pan_Genome_QC_Filtered.h5ad` (`189.393 GB`) exceed the `60.00 GB` team download cap.
- Disjoint 50/50 split-halves (`half_a_delta_*`, `half_b_delta_*`) split both non-targeting controls and perturbed cells independently per batch (strictly zero shared cells between Half A and Half B, as required for noise-ceiling estimation).

## 6. Saved File Schema (`p2/shared/lines/<line>/`)
- `genes.txt`: Aligned HGNC gene symbols matching column axis `0..G-1`.
- `targets.txt`: Deduplicated target HGNC gene symbols matching row axis `0..P-1`.
- `ctrl_cells.npz`: `X_log1p` (`log1p(1e4 * x / UMI)`), `X_zscore` (Addendum-B per-batch NT z-score), `X_libnorm` (`1e4 * x / UMI`), `ctrl_mean_log1p`, `ctrl_mean_libnorm`, `n_ctrl_total`, `batch_labels`.
- `pseudobulk_delta.npz`: `delta_log1p`, `delta_zscore`, `delta_libnorm` (all relative to per-batch NT control means), disjoint 50/50 split-halves (`half_a_delta_*`, `half_b_delta_*`, `half_valid`), `targets`, `target_ensembl_ids`, `reproducible`, `n_cells`, `ad_counts`, `energy_test_p_value`.
- `line_metadata.json`: Full provenance, SHA-256 of raw input `.h5ad`, and alignment/filtering loss counts.
