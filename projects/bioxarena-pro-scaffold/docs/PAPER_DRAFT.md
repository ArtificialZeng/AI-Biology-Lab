# Residual-Space Cross-Line Perturbation Modeling for Unseen CRISPRi Knockdowns: Metric Trade-Offs and External Generalization

**Working Manuscript Draft (Round 6 Revision — T1 Remediated; `[T3-CI 待填]` / Phase B Results Pending Backfill)**  
**Primary Benchmark**: BioXArena-Pro `perturbation/xatlas-orion-perturb` (HCT116, HEK293T)  
**External Independent Dataset**: Replogle et al. (*Cell* 2022) Perturb-seq (K562 GWPS, RPE1 Essential; Figshare+ 20029387)  
**Frozen Core Solver**: `src/xatlas_panhuman_solver.py` (hyperparameters frozen; path resolution parameterized via `BIOXARENA_DIR`)

---

## Abstract

Evaluating transcriptional response predictors for unseen CRISPRi gene knockdowns depends critically on whether dominant systemic stress axes are retained or projected out. Under standard unprojected $\Delta$-expression Pearson correlation, simple low-rank linear baselines outperform many deep-learning models (Ahlmann-Eltze et al., *Nature Methods* 2025), in part because systemic variation shared across perturbations can yield high correlation without recovering pathway-specific structure. Conversely, the BioXArena-Pro `xatlas-orion-perturb` v5 evaluation metric projects both predictions and ground truth onto the orthogonal complement of $\mathrm{span}\{\mathbf{1}, b_c, a_c\}$, isolating perturbation-specific residual responses.

We study a cell-line-agnostic kernel and cross-line bridging pipeline (`Cross-Line Perturb-Bridge`, implemented in `PanPerturbBridgeModel`) combining: (S1) self-knockdown spike masking, baseline subtraction, $\sqrt{n_{\text{cells}}}$ row weighting, and soft-thresholded Hybrid 1-hop/2-hop Topological Overlap (`HybridTOM`, power 1.5); (S2) a multi-scale spectral deflation pyramid (`d0..d10`); (S3) cross-line gene-graph fusion across shared genes; (S4) empirical-SNR weighted non-reproducible anchors; (S5) a `GeneSlope-CosKRR` cross-line response translator and bidirectional bridge; (S6) thresholded Reciprocal Response Borrowing over reproducible knockdowns; and (S7) output SVD Wiener shrinkage with empirical manifold snapping.

On BioXArena-Pro `xatlas-orion-perturb`, across 10 development splits, S7 increases official v5 orthogonal-residual Pearson correlation from `0.0613 ± 0.0040` (`100-SVD + Ridge`, S0) and `0.0672 ± 0.0052` (Ahlmann-Eltze bilinear linear baseline) to **`0.1099 ± 0.0087`** (row-level 95% CI: `[T3-CI 待填]`). On 10 additional random resplits sharing the same reproducible gene pool (`Seeds 101–110`), S7 yields `0.1073 ± 0.0072` vs S0 `0.0596 ± 0.0056`; under a strict simultaneous union-held-out + refitted-axis protocol, S7 reaches `0.1115 ± 0.0089` vs S0 `0.0658 ± 0.0033` (row-level 95% CI: `[T3-CI 待填]`).

Evaluating the **exact frozen model (zero hyperparameter retuning)** on an independent external atlas—**Replogle et al. (*Cell* 2022) K562 & RPE1** (10 strict simultaneous union-held-out splits with training-only refitted axes)—reveals a sharp, instructive **metric and module transferability split**:
1. **Residual vs. traditional metrics**: On the v5 orthogonal-residual metric, frozen S5 (`0.1401 ± 0.0057`) and S7 (`0.1359 ± 0.0076`) substantially outperform `100-SVD + Ridge` (`0.0883 ± 0.0048`), the Ahlmann-Eltze bilinear baseline (`0.0967 ± 0.0053`, row-level 95% CI `[T3-CI 待填]`), `30-kNN` (`0.0971 ± 0.0059`), and self-implemented untuned simplified neural baselines (`ResMLP 0.1091`, `GNN 0.1233`, row-level 95% CI `[T3-CI 待填]`). However, on **traditional unprojected metrics**, S7 **underperforms** linear baselines: standard $\Delta$-expression Pearson (`0.1162` vs `30-kNN 0.1479`, `AE-Bilinear 0.1470`, `S0 0.1408`), Top-50 DE Pearson (`0.2013` vs `S0 0.2119`, `AE 0.2026`; S5 reaches `0.2170`), and Top-50 sign accuracy (`0.5962` vs `S0 0.6256`, `AE 0.6255`).
2. **Module transferability**: Externally on Replogle 2022, most of the v5 gain transfers from **S1** (`+0.0356`) and **S5 cross-line bridging** (`+0.0134`), whereas **S2 spectral deflation (`-0.0024`)** and **S6 Reciprocal Response Borrowing (`-0.0109`)** are **negative contributors** under frozen hyperparameters.

---

## 1. Main Quantitative Results (All Metrics Including Unfavorable Outcomes)

### 1.1 BioXArena-Pro `xatlas-orion-perturb` (HCT116 & HEK293T)

All scores below use the official BioXArena-Pro `tasks.utils.strict_group_mean_from_row_scores` grader and exact `grade.py` QR orthogonal projection. Because random splits share overlapping training sets, split-level paired $t$-tests are invalid; statistical uncertainty is quantified strictly by **per-split row-level stratified paired bootstrap 95% confidence intervals (4,000 resamples)**.

#### Table 1. Progressive Ablation on `xatlas-orion-perturb` (10 Development Splits, Seeds `42..31337`)

| Stage | Component Added | 10-Split v5 Mean ± SD | HCT116 | HEK293T | Δ vs S0 | Row-Level 95% CI |
|:---:|:---|:---:|:---:|:---:|:---:|:---:|
| **S0** | `100-SVD + Ridge (α=10)` (Official Ref) | 0.0613 ± 0.0040 | 0.0767 | 0.0460 | — | — |
| **S1** | `+ SelfMask + BaseSub + √n + HybridTOM(1.5)` | 0.0677 ± 0.0068 | 0.0866 | 0.0487 | +0.0063 | `[T3-CI 待填]` |
| **S2** | `+ Multi-Scale Spectral Deflation (d0..d10)` | 0.0722 ± 0.0061 | 0.0935 | 0.0509 | +0.0109 | — |
| **S3** | `+ Cross-Line Gene Graph (4,040 shared genes)` | 0.0748 ± 0.0053 | 0.0941 | 0.0555 | +0.0135 | — |
| **S4** | `+ Empirical-SNR Non-Reproducible Anchors` | 0.0811 ± 0.0038 | 0.1033 | 0.0589 | +0.0197 | — |
| **S5** | `+ GeneSlope-CosKRR Cross-Line Bridge` | 0.0887 ± 0.0046 | 0.1138 | 0.0635 | +0.0273 | `[T3-CI 待填]` |
| **S6** | `+ Reciprocal Response Borrowing (6,018 rep)` | 0.1027 ± 0.0065 | 0.1371 | 0.0682 | +0.0413 | `[T3-CI 待填]` |
| **S7** | `+ SVD Wiener Denoise + Manifold Snapping` | **0.1099 ± 0.0087** | **0.1462** | **0.0736** | **+0.0485** | `[T3-CI 待填]` |

#### Table 2. Baseline Comparisons, Additional Random Resplits, and Strict Refit Audit (`xatlas-orion-perturb`)

*(Note: Neural networks below are self-implemented, untuned simplified baselines trained on small in-line sample sizes, NOT official GEARS / scGPT / STATE implementations. The 10 development splits were also used during hyperparameter design; comparison against untuned baselines on those splits favors S7.)*

| Method / Evaluation Protocol | v5 Mean ± SD | Δ (`S7` vs Method) | Row-Level 95% CI |
|:---|:---:|:---:|:---:|
| Self-implemented untuned simplified ResMLP (10 dev splits) | 0.0482 ± 0.0052 | +0.0617 | — |
| Self-implemented untuned simplified single-line GNN (10 dev splits) | 0.0528 ± 0.0088 | +0.0571 | — |
| Self-implemented untuned simplified cross-line GNN (10 dev splits) | 0.0538 ± 0.0086 | +0.0561 | — |
| Official `100-SVD + Ridge` (10 dev splits) | 0.0613 ± 0.0040 | +0.0485 | `[T3-CI 待填]` |
| Ahlmann-Eltze Bilinear Linear Baseline (10 dev splits) | 0.0672 ± 0.0052 | +0.0427 | `[T3-CI 待填]` |
| **S7 (10 dev splits)** | **0.1099 ± 0.0087** | — | — |
| **额外随机重切分（与开发分割共享基因池，Seeds 101–110）：S0 → S7** | `0.0596` → **`0.1073 ± 0.0072`** | `+0.0477` | — |
| **合并 20 次随机切分（共享基因池）：S0 → S7** | `0.0605` → **`0.1086 ± 0.0080`** | `+0.0481` | — |
| **Strict Simultaneous Union (~360 genes) + Refit Axis (10 splits): S0 → S7** | `0.0658` → **`0.1115 ± 0.0089`** | `+0.0457` | `[T3-CI 待填]` |

---

### 1.2 External Independent Evaluation on Replogle et al. (*Cell* 2022) K562 & RPE1 (Zero Hyperparameter Retuning)

Frozen `PanPerturbBridgeModel` was evaluated without modifying a single constant across 10 strict simultaneous union-held-out splits (20% per cell line, ~488 unique genes removed from both RPE1 and K562) with training-only refitted `context_baseline` and `residual_axis`. **All four computed metrics are reported in Table 3**, including the three traditional metrics where S7 underperforms linear baselines (a 5th v5-independent perturbation discrimination metric is scheduled in Phase B T2).

#### Table 3. Full Four-Metric External Validation on Replogle et al. (*Cell* 2022, 10 Strict Union Splits)

| Method (Frozen / Untuned) | v5 Residual $r$ (Mean ± SD) | Δ-Expr Pearson $r$ (Ahlmann-Eltze) | Top-50 DE Pearson $r$ | Top-50 Sign Accuracy | Row-Level 95% CI (v5) |
|:---|:---:|:---:|:---:|:---:|:---:|
| Constant Mean | 0.0009 ± 0.0074 | -0.0322 ± 0.0166 | -0.0531 | 0.4992 | — |
| `100-SVD + Ridge (S0)` | 0.0883 ± 0.0048 | 0.1408 ± 0.0055 | 0.2119 | **0.6256** | — |
| Ahlmann-Eltze Bilinear Linear (AE) | 0.0967 ± 0.0053 | 0.1470 ± 0.0068 | 0.2026 | 0.6255 | — |
| `30-kNN` Co-expression | 0.0971 ± 0.0059 | **0.1479 ± 0.0080** | 0.1822 | 0.6095 | — |
| Self-implemented untuned simplified ResMLP | 0.1091 ± 0.0067 | 0.0876 ± 0.0062 | 0.1896 | 0.5901 | — |
| Self-implemented untuned simplified GNN | 0.1233 ± 0.0061 | 0.0971 ± 0.0055 | 0.1892 | 0.5999 | — |
| **S1** (`+ SelfMask + BaseSub + √n + HybridTOM`) | 0.1239 ± 0.0066 | 0.1054 ± 0.0062 | 0.1917 | 0.5968 | S1 vs S0: `[T3-CI 待填]` |
| **S2** (`+ Spectral Deflation`, **negative**) | 0.1215 ± 0.0052 | 0.1033 ± 0.0050 | 0.1918 | 0.5922 | — |
| **S3** (`+ Cross-Line Graph`) | 0.1267 ± 0.0051 | 0.1079 ± 0.0049 | 0.2067 | 0.5996 | — |
| **S4** (`+ Non-Rep Anchors`, **flat**) | 0.1267 ± 0.0057 | 0.1077 ± 0.0053 | 0.2088 | 0.5995 | — |
| **S5** (`+ Cross-Line Bridge`, **external peak**) | **0.1401 ± 0.0057** | 0.1193 ± 0.0053 | **0.2170** | 0.6064 | S5 vs AE: `[T3-CI 待填]` |
| **S6** (`+ Reciprocal Borrowing`, **negative**) | 0.1292 ± 0.0064 | 0.1106 ± 0.0057 | 0.2155 | 0.6020 | S6 vs S5: `[T3-CI 待填]` |
| **S7** (`+ Wiener + Manifold Snap`, full frozen) | 0.1359 ± 0.0076 | 0.1162 ± 0.0066 | 0.2013 | 0.5962 | S7 vs AE: `[T3-CI 待填]`; S7 vs GNN: `[T3-CI 待填]` |

**Key Empirical Takeaways (Honest Scope)**:
1. **Metric Trade-Off**: Frozen S5/S7 clearly beats linear and simplified neural baselines on the residualized **v5** metric (note that adding identical S1 self-mask/baseline preprocessing to the bilinear baseline does not close the gap), yet **loses** to linear baselines on unprojected $\Delta$-expression Pearson (`0.1162` vs `0.1470–0.1479`), Top-50 DE Pearson (`0.2013` vs `0.2119`), and Top-50 sign accuracy (`0.5962` vs `0.6256`). Because S1–S7 explicitly project training targets orthogonal to `residual_axis`, they discard the shared stress component that contributes substantially to unprojected correlation and sign direction on top DE genes. Full 5-metric analysis including perturbation discrimination is performed in Phase B (T2).
2. **Module Transferability**: On Replogle 2022, external transfer is driven primarily by **S1 (`+0.0356`)** and **S5 cross-line bridging (`+0.0134`)**; **S2 (`-0.0024`)** and **S6 (`-0.0109`)** reduce v5 performance externally. Leave-one-module-out ablation across both datasets is scheduled in Phase B (T4).

---

## 2. Biological Enrichment & Case Studies (Audit Corrections)

- **Correction to Round 5 Reactome Analysis**: Online Reactome queries in commit `a9ab830` used the ~20,000-gene whole-genome background, whereas the evaluated Replogle gene pool (~1,784 genes) is already pre-enriched 8–10× for essential macromolecular complexes. Moreover, the top-236 genes predicted by baseline S0 (`100-SVD + Ridge`) exhibit virtually identical proteasome (`12.3%` vs S7 `11.9%`) and ribosomal (`23.7%` vs S7 `21.2%`) proportions (172/236 overlap). Therefore, prior whole-genome enrichment claims are **withdrawn**; proper evaluated-universe hypergeometric + BH enrichment (T5) and CRISPRi genomic proximity testing (T6) are executed in Phase B.
- **Qualitative Single-Gene Examples (Selected from All Evaluated Genes)**: While S6 is negative overall on Replogle (`-0.0109`), individual high-gain cases correspond to known biological complexes (e.g., `DICER1` borrowing `DGCR8`/`DROSHA` Microprocessor profiles; `NUDT21` borrowing `CPSF6` CFIm complex; `MYC` borrowing `TELO2`/`TTI1`/`TRRAP`). These are strictly qualitative case illustrations, not global mechanism proof, and genomic cis-proximity confounding is directly tested in Phase B (T6).

---

## 3. Limitations & Strict Comparability Statement

1. **No Private Codabench Test Score Yet**: All xatlas results are on internal public-train splits. The frozen Codabench submission pack (`commit 2a8b2e0`, SHA-256 verified, preflight passed) awaits a single upload by the human PI; no private test score is predicted or claimed.
2. **Comparison Baselines**: Official implementations of GEARS, scGPT, and STATE have not yet been benchmarked; neural baselines in this draft are self-implemented, untuned simplified networks.
