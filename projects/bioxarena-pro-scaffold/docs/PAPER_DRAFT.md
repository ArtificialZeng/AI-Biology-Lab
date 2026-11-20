# Residual-Space Cross-Line Perturbation Modeling for Unseen CRISPRi Knockdowns: Metric Trade-Offs and External Generalization

**Working Manuscript Draft (Round 6 Phase A — T1/T3/T7 Complete; Phase B T2/T4 Pending)**  
**Primary Benchmark**: BioXArena-Pro `perturbation/xatlas-orion-perturb` (HCT116, HEK293T)  
**External Independent Dataset**: Replogle et al. (*Cell* 2022) Perturb-seq (K562 GWPS, RPE1 Essential; Figshare+ 20029387)  
**Official Codabench Private Test (Single Frozen Submission `1062221`)**: `0.0998` (leading within the `0.02` noise band; human-AI multi-round iteration, NOT directly comparable to single-shot 2h autonomous agents)

---

## Abstract

Evaluating transcriptional response predictors for unseen CRISPRi gene knockdowns depends critically on whether dominant systemic stress axes are retained or projected out. Under standard unprojected $\Delta$-expression Pearson correlation, simple low-rank linear baselines outperform many deep-learning models (Ahlmann-Eltze et al., *Nature Methods* 2025). Conversely, the BioXArena-Pro `xatlas-orion-perturb` v5 evaluation metric projects both predictions and ground truth onto the orthogonal complement of $\mathrm{span}\{\mathbf{1}, b_c, a_c\}$, isolating perturbation-specific residual responses.

We evaluate a frozen, cell-line-agnostic kernel and cross-line bridging pipeline (`Cross-Line Perturb-Bridge`, implemented in `PanPerturbBridgeModel`) across two human Perturb-seq atlases (HCT116/HEK293T in `xatlas-orion-perturb`; K562/RPE1 in Replogle et al., *Cell* 2022). Because random splits share overlapping training sets, we quantify statistical uncertainty strictly via **per-split row-level stratified paired bootstrap 95% confidence intervals (4,000 resamples)**.

On BioXArena-Pro `xatlas-orion-perturb`, across 10 development splits, S7 increases official v5 orthogonal-residual Pearson correlation from `0.0613 ± 0.0040` (`100-SVD + Ridge`, S0) and `0.0672 ± 0.0052` (Ahlmann-Eltze bilinear linear baseline, AE) to **`0.1099 ± 0.0087`** (S7 vs S0 mean Δ `+0.0485`, median 95% CI `[+0.0318, +0.0677]`, 10/10 splits CI > 0; S7 vs AE mean Δ `+0.0427`, median 95% CI `[+0.0256, +0.0606]`, 10/10 splits CI > 0). Under strict simultaneous union-held-out + refitted-axis evaluation, S7 reaches `0.1115 ± 0.0089` vs S0 `0.0658 ± 0.0033` (mean Δ `+0.0457`, median 95% CI `[+0.0276, +0.0635]`, 10/10 splits CI > 0). On a single official Codabench private test submission (`Submission ID 1062221`, frozen at commit `2a8b2e0`), S7 scored **`0.0998`** (measured validation-to-test shift `−0.0101`), `+0.0173` above the best 2-hour autonomous agent run (`0.0825`), which lies within the official `0.02` noise threshold and represents a different development category (human-AI multi-round iteration with ~25 tuned constants).

Evaluated **zero-shot with 100% frozen hyperparameters** on the independent **Replogle et al. (*Cell* 2022) K562 & RPE1** atlas (10 strict simultaneous union-held-out splits with training-only refitted axes), we report **S5 and S7 as co-primary external variants** alongside all four computed metrics:
1. **Residual vs. traditional metric trade-off**: On the v5 residual metric, frozen **S5** (`0.1401 ± 0.0057`; S5 vs AE mean Δ `+0.0434`, median 95% CI `[+0.0323, +0.0549]`, 10/10 CI > 0; S5 vs simplified GNN `+0.0168`, median 95% CI `[+0.0038, +0.0302]`, 7/10 CI > 0) and **S7** (`0.1359 ± 0.0076`; S7 vs AE `+0.0392`, median 95% CI `[+0.0240, +0.0537]`, 10/10 CI > 0; S7 vs simplified GNN `+0.0126`, median 95% CI `[-0.0041, +0.0286]`, only 3/10 CI > 0) substantially outperform linear baselines (`S0 0.0883`, `AE 0.0967`, `30-kNN 0.0971`). However, on **traditional unprojected metrics**, S7 **underperforms** linear baselines: standard $\Delta$-expression Pearson (`0.1162` vs `30-kNN 0.1479`, `AE 0.1470`, `S0 0.1408`), Top-50 DE Pearson (`0.2013` vs `S0 0.2119`, `AE 0.2026`; S5 reaches `0.2170`), and Top-50 sign accuracy (`0.5962` vs `S0 0.6256`, `AE 0.6255`).
2. **Module transferability**: Externally on Replogle 2022, v5 gains transfer primarily from **S1** (`+0.0356`, median 95% CI `[+0.0256, +0.0450]`, 10/10 CI > 0) and **S5 cross-line bridging** (`+0.0134`), whereas **S2 spectral deflation (`−0.0024`)** and **S6 Reciprocal Response Borrowing (`−0.0109`, median 95% CI `[-0.0169, -0.0050]`, 9/10 CI < 0)** are **negative contributors**.

---

## 1. Main Quantitative Results (All Metrics Including Unfavorable Outcomes)

### 1.1 BioXArena-Pro `xatlas-orion-perturb` (HCT116 & HEK293T)

#### Table 1. Progressive Ablation on `xatlas-orion-perturb` (10 Development Splits, Seeds `42..31337`)

| Stage | Component Added | 10-Split v5 Mean ± SD | HCT116 | HEK293T | Δ vs S0 | Row-Level Median 95% CI (Splits > 0) |
|:---:|:---|:---:|:---:|:---:|:---:|:---:|
| **S0** | `100-SVD + Ridge (α=10)` (Official Ref) | 0.0613 ± 0.0040 | 0.0767 | 0.0460 | — | — |
| **S1** | `+ SelfMask + BaseSub + √n + HybridTOM(1.5)` | 0.0677 ± 0.0068 | 0.0866 | 0.0487 | +0.0063 | `[-0.0015, +0.0174]` (4/10) |
| **S2** | `+ Multi-Scale Spectral Deflation (d0..d10)` | 0.0722 ± 0.0061 | 0.0935 | 0.0509 | +0.0109 | — |
| **S3** | `+ Cross-Line Gene Graph (4,040 shared genes)` | 0.0748 ± 0.0053 | 0.0941 | 0.0555 | +0.0135 | — |
| **S4** | `+ Empirical-SNR Non-Reproducible Anchors` | 0.0811 ± 0.0038 | 0.1033 | 0.0589 | +0.0197 | — |
| **S5** | `+ GeneSlope-CosKRR Cross-Line Bridge` | 0.0887 ± 0.0046 | 0.1138 | 0.0635 | +0.0273 | vs AE: `[+0.0091, +0.0317]` (10/10) |
| **S6** | `+ Reciprocal Response Borrowing (6,018 rep)` | 0.1027 ± 0.0065 | 0.1371 | 0.0682 | +0.0413 | vs S5: `[+0.0050, +0.0244]` (9/10) |
| **S7** | `+ SVD Wiener Denoise + Manifold Snapping` | **0.1099 ± 0.0087** | **0.1462** | **0.0736** | **+0.0485** | vs S0: **`[+0.0318, +0.0677]`** (10/10) |

#### Table 2. Baseline Comparisons, Additional Random Resplits, Strict Refit Audit & Official Private Test (`xatlas-orion-perturb`)

*(Neural baselines below are self-implemented, untuned simplified baselines, NOT official GEARS / scGPT / STATE implementations. The 10 development splits were also used during hyperparameter selection.)*

| Method / Evaluation Protocol | v5 Score | Δ (`S7` vs Method) | Row-Level Median 95% CI (Splits > 0) |
|:---|:---:|:---:|:---:|
| Self-implemented untuned simplified ResMLP (10 dev splits) | 0.0482 ± 0.0052 | +0.0617 | — |
| Self-implemented untuned simplified single-line GNN (10 dev splits) | 0.0528 ± 0.0088 | +0.0571 | — |
| Self-implemented untuned simplified cross-line GNN (10 dev splits) | 0.0538 ± 0.0086 | +0.0561 | — |
| Official `100-SVD + Ridge` (10 dev splits) | 0.0613 ± 0.0040 | +0.0485 | `[+0.0318, +0.0677]` (10/10) |
| Ahlmann-Eltze Bilinear Linear Baseline (10 dev splits) | 0.0672 ± 0.0052 | +0.0427 | `[+0.0256, +0.0606]` (10/10) |
| **S7 (10 dev splits)** | **0.1099 ± 0.0087** | — | — |
| **额外随机重切分（与开发分割共享基因池，Seeds 101–110）：S0 → S7** | `0.0596` → `0.1073 ± 0.0072` | +0.0477 | — |
| **Strict Simultaneous Union (~360 genes) + Refit Axis (10 splits): S0 → S7** | `0.0658` → `0.1115 ± 0.0089` | +0.0457 | `[+0.0276, +0.0635]` (10/10) |
| **Official Codabench Private Test (Single Submission `1062221`, human-AI multi-round)** | **`0.0998`** | `+0.0173` vs best 2h agent `0.0825` | Within `0.02` noise margin (different dev category) |

---

### 1.2 External Independent Evaluation on Replogle et al. (*Cell* 2022) K562 & RPE1 (Zero Hyperparameter Retuning)

**Preprocessing Disclosure (`src/prepare_replogle2022_cache.py`)**: Top 5,000 finite readout genes per line were selected via `var/clean_mean` from the published `h5ad` (unsupervised expression-level filter; `4,991/5,000` RPE1 and `4,996/5,000` K562 overlap with `var/mean`). Perturbations are labeled `reproducible = 1` iff `anderson_darling_counts >= 50` AND `energy_test_p_value < 0.001`. When multiple guides targeted the same Ensembl gene, a single row was retained maximizing `(reproducible, ad_counts, num_cells_filtered)`—note that this deduplication rule uses perturbation-level summary statistics across the screen and is disclosed here for full methodological transparency.

#### Table 3. Full Four-Metric External Validation on Replogle et al. (*Cell* 2022, 10 Strict Union Splits; S5 & S7 Co-Primary)

| Method (Frozen / Untuned) | v5 Residual $r$ | Δ-Expr Pearson $r$ | Top-50 DE Pearson $r$ | Top-50 Sign Acc | Row-Level Median 95% CI (v5) |
|:---|:---:|:---:|:---:|:---:|:---:|
| Constant Mean | 0.0009 ± 0.0074 | -0.0322 ± 0.0166 | -0.0531 | 0.4992 | — |
| `100-SVD + Ridge (S0)` | 0.0883 ± 0.0048 | 0.1408 ± 0.0055 | 0.2119 | **0.6256** | — |
| Ahlmann-Eltze Bilinear Linear (AE) | 0.0967 ± 0.0053 | 0.1470 ± 0.0068 | 0.2026 | 0.6255 | — |
| `30-kNN` Co-expression | 0.0971 ± 0.0059 | **0.1479 ± 0.0080** | 0.1822 | 0.6095 | — |
| Self-implemented untuned simplified ResMLP | 0.1091 ± 0.0067 | 0.0876 ± 0.0062 | 0.1896 | 0.5901 | — |
| Self-implemented untuned simplified GNN | 0.1233 ± 0.0061 | 0.0971 ± 0.0055 | 0.1892 | 0.5999 | — |
| **S1** (`+ SelfMask + BaseSub + √n + HybridTOM`) | 0.1239 ± 0.0066 | 0.1054 ± 0.0062 | 0.1917 | 0.5968 | S1 vs S0: `[+0.0256, +0.0450]` (10/10) |
| **S2** (`+ Spectral Deflation`, **negative −0.0024**) | 0.1215 ± 0.0052 | 0.1033 ± 0.0050 | 0.1918 | 0.5922 | — |
| **S3** (`+ Cross-Line Graph`) | 0.1267 ± 0.0051 | 0.1079 ± 0.0049 | 0.2067 | 0.5996 | — |
| **S4** (`+ Non-Rep Anchors`, **flat 0.0000**) | 0.1267 ± 0.0057 | 0.1077 ± 0.0053 | 0.2088 | 0.5995 | — |
| **S5** (`+ Cross-Line Bridge`, **co-primary / external peak**) | **0.1401 ± 0.0057** | 0.1193 ± 0.0053 | **0.2170** | 0.6064 | S5 vs AE: `[+0.0323, +0.0549]` (10/10); S5 vs GNN: `[+0.0038, +0.0302]` (7/10) |
| **S6** (`+ Reciprocal Borrowing`, **negative −0.0109**) | 0.1292 ± 0.0064 | 0.1106 ± 0.0057 | 0.2155 | 0.6020 | S6 vs S5: `[-0.0169, -0.0050]` (9/10 CI < 0) |
| **S7** (`+ Wiener + Manifold Snap`, **co-primary / full frozen**) | 0.1359 ± 0.0076 | 0.1162 ± 0.0066 | 0.2013 | 0.5962 | S7 vs AE: `[+0.0240, +0.0537]` (10/10); S7 vs GNN: `[-0.0041, +0.0286]` (3/10) |

---

## 2. Pathway Enrichment Re-Evaluation (T5) & CRISPRi Proximity Audit (T6)

- **T5 Local `ReactomePathways.gmt` Enrichment (`M_bg = 1,784` evaluated unique genes)**: Confirms that `S7_top_164` and `S0_top_164` overlap heavily (`113/164`, `68.9%`) and both enrich the same amino-acid/ribosome/proteasome pathways (prior whole-genome enrichment claims withdrawn). Strictly against the 1,784-gene evaluated background, `S7 − S0` Top-150 yields 6 pathways at BH-FDR < 0.05 (`mRNA Splicing - Major Pathway`, `38/145`, `3.12×`, `FDR = 1.18e-08`), `S6 − S5` Top-150 yields 21 pathways (`mRNA Splicing`, `34/150`, `2.70×`, `FDR = 9.64e-06`), and `S5 − S2` Top-150 yields 5 pathways (`Mitochondrial translation`, `20/82`, `2.90×`, `FDR = 2.76e-03`).
- **T6 CRISPRi Genomic Proximity Audit (Replogle 2022 + GENCODE v46 GRCh38 GTF, 10,530/10,539 genes mapped)**: Active `z < −4` donor-query pairs (`n = 3,939`) are enriched **474×** for `< 10 kb` TSS proximity (`1.422%` vs `0.003%` random null, `3.32%` of negative-z weight), confirming real CRISPRi local cis-repression/bidirectional-promoter confounding in negative-z pairs. On Replogle itself, removing all `< 10 kb` (and `< 100 kb`) donors retains **`97.41%`** of top-40 ownRec gains (`0.2536 → 0.2471`), while overall `S6` remains negative (`−0.0105 → −0.0108`). Because xatlas gene identities cannot be examined, local cis-repression remains an unexcluded confounder for xatlas's `S6` gain (where 78% of top-40 winners have negative-z top donors).
