# Round 7 Pre-Registration: Unseen-Dataset Confirmatory Validation (R7-0)

- **Registration Date**: 2026-11-21
- **Status**: **LOCKED PRIOR TO ANY DATA ACQUISITION OR INSPECTION OF CONFIRMATION EXPRESSION MATRICES**
- **Frozen Codebase Reference Commit**: `3609811` (`main`, verified by Master Controller in `审核报告_第六轮.md`)
- **Master Controller**: Claude | **Executor**: Antigravity (Gemini) | **External Release Authority**: Human PI

---

## 1. Frozen Methods (Zero Hyperparameter or Constant Modification Permitted)

All solver constants and algorithmic implementations remain strictly frozen at commit `3609811` (`src/xatlas_panhuman_solver.py`, `PanPerturbBridgeModel`, SHA-256 unchanged from Round 6 T0/T7 freeze). External datasets are strictly evaluation-only.

| Role | Method Name | Exact Frozen Implementation (`src/xatlas_panhuman_solver.py` @ `3609811`) | Pre-Registered Interpretation Rule |
|---|---|---|---|
| **Primary Method** | **`S5`** (`Cross-Line Bridge`) | `PanPerturbBridgeModel.fit_predict(..., stage=5)` (includes S1 HybridTOM $p=1.5$, S2 Multi-Scale Deflation, S3 Precision Weighting, S4 Non-Reproducible Anchors, S5 Cross-Line Bridge; **excludes** S6 Reciprocal Borrowing and S7 Wiener/Manifold Snapping) | Sole method used to adjudicate Primary Hypotheses **H1, H2, H3** |
| **Secondary Method** | **`S7`** (`Full Pipeline`) | `PanPerturbBridgeModel.fit_predict(..., stage=7)` (full Codabench-frozen pipeline including S6 and S7) | Reported in all main tables alongside S5 for full comparability with xatlas and Replogle |
| **Exploratory Only** | **`S7_no_S6`** (`LOO_S6_Reciprocal`) | `round6_t4_leave_one_out_ablation.py` `LOO_S6_Reciprocal` (`S7` with `use_s6=False`) | **Strictly exploratory**; because its v5 advantage (`0.1512`) was discovered after observing Replogle 2022 in Round 6 T4, it MUST NOT be used to claim confirmatory validation |
| **Reference Baselines** | `Constant_Mean`, `kNN_30`, `AE` (`AE_Bilinear`), `S0` (`SVD100_Ridge`), `Simplified_ResMLP`, `Simplified_GNN` | Identical to `src/round6_t2_multi_metric_bench.py` and `src/replogle_external_validation.py` @ `3609811` | Evaluated on identical splits and metrics; `AE` (Ahlmann-Eltze et al. *Nature Methods* 2025 bilinear baseline) is the primary comparator for H1–H3 |

> Note: Official external deep learning / foundation model comparisons (`R7-2`, official `cell-gears` package and official STATE/scGPT if feasible) will be run separately in `R7-2` on the exact same splits and metrics, and never conflated with `Simplified_ResMLP` / `Simplified_GNN`.

---

## 2. Unseen Confirmation Dataset & Pre-Registered Fallback Rule

### 2.1 Primary Target Dataset
- **Nadig et al. (2025, *Nature Genetics*)**: Transcriptome-wide essential-gene CRISPRi Perturb-seq across two human cell lines — **HepG2** and **Jurkat** (paired into the two-line framework slot `(line0=HepG2, line1=Jurkat)`, analogous to `(HCT116, HEK293T)` on xatlas and `(RPE1, K562)` on Replogle 2022).
- **Zero-Look Guarantee**: At the time of this commit, **no Nadig et al. (2025) expression matrix has been downloaded, opened, inspected, or evaluated**. Upon download from the official repository (Figshare / GEO / Zenodo / scPerturb official release), the exact download URL, dataset version, license, and SHA-256 checksum will be recorded prior to running preprocessing.

### 2.2 Pre-Registered Fallback Rule (Triggered Only If Primary Target Is Inaccessible)
If official HepG2 + Jurkat pseudo-bulk or single-cell Perturb-seq expression matrices from Nadig et al. (2025) cannot be downloaded, are license-restricted/broken, or lack compatible gene/perturbation metadata required to form a clean two-cell-line paired perturbation dataset:
1. Document the exact technical/access failure transparently in `reports/ROUND7_CHECKPOINT_1.md` and `reports/ROUND7_EXECUTION_REPORT.md`.
2. Fall back to the **Replogle et al. (2022) `K562_essential` screen** (independent essential-gene screen from Figshare+ DOI `10.25452/figshare.plus.20029387`, distinct from `K562_gwps_normalized_bulk_01.h5ad`), paired with `RPE1` or evaluated under the exact same strict union holdout protocol, explicitly noting the partial cell-line overlap limitation if fallback is invoked.

### 2.3 Preprocessing Specification (Locked to `prepare_replogle2022_cache.py` @ `3609811`)
`src/prepare_nadig2025_cache.py` will mirror `src/prepare_replogle2022_cache.py` (`3609811`) with zero dataset-specific retuning:
1. **Readout Gene Selection (`N = 5,000` per cell line)**: Drop columns with any non-finite entry (`NaN`/`Inf`), then select the top 5,000 finite genes ranked by baseline/control mean expression (`var/clean_mean` if present, otherwise control-cell mean expression; if pseudo-bulk z-scores lack `var/clean_mean`, use `var/mean` or expression-mean metadata as available in the official file, documenting the exact column used and verifying zero perturbation-label usage).
2. **Perturbation Filtering & Reproducibility Label**: Exclude non-targeting controls and rows without valid gene identifiers (Ensembl ID `ENSG...` or official HGNC symbol harmonized consistently across both lines). Define `reproducible = 1` using the identical statistical criterion (`anderson_darling_counts >= 50` AND `energy_test_p_value < 0.001`) when those official columns are provided; if Nadig 2025 provides different official QC/DE-count columns (e.g., `mitopercent` / `n_cells` / `num_de_genes` / `energy_distance`), map the exact equivalent official strong-perturbation threshold (e.g., `>= 50` differentially expressed genes and `p < 0.001` where available, or top strong perturbations strictly by official cell-count/DE-test metadata without looking at model scores) and record the exact rule BEFORE running any model.
3. **Multi-Guide / Duplicate Target Deduplication**: Strictly 1 row per unique target gene per cell line, retaining the single highest-ranked guide/entry by `(reproducible, de_or_ad_count, n_cells)`.

---

## 3. Pre-Registered Hypotheses & Decision Criteria (Strictly Locked)

All three hypotheses compare **Primary Method `S5`** against **`AE` (Ahlmann-Eltze et al. 2025 bilinear baseline)** across the **10 strict simultaneous union-held-out splits** (`SEEDS_10 = [0, 1, 2, 3, 7, 13, 21, 42, 99, 2024]`).

| Hypothesis | Priority | Metric | Exact Quantitative Criterion | Outcome Rule |
|---|---|---|---|---|
| **H1** | **Primary** | Raw-space Perturbation Discrimination (`disc_raw`) | Mean `Δ(S5 − AE) > 0` across 10 splits **AND** `≥ 7 / 10` splits have row-level stratified paired bootstrap 95% CI strictly `> 0` (`ci95_lo > 0`) | Pass / Fail |
| **H2** | **Primary** | Standard Unprojected Δ-Expression Pearson (`delta_r`) | Mean `Δ(S5 − AE) < 0` across 10 splits (confirming **metric inversion**: linear baseline outperforms S5 on unprojected Δ-expression Pearson while losing on perturbation-specific discrimination) | Pass / Fail |
| **Main Thesis Confirmation** | **Joint Primary Decision** | **H1 AND H2 simultaneously** | **BOTH H1 AND H2 must pass** to declare the core thesis **「确认 (Confirmed)」** on the unseen dataset; if either H1 or H2 fails, declare **「未确认 (Not Confirmed)」** | Confirmed / Not Confirmed |
| **H3** | **Secondary** | Orthogonal-Residual Pearson (`v5_orth_r`) | `≥ 7 / 10` splits have row-level stratified paired bootstrap 95% CI strictly `> 0` (`ci95_lo > 0`) for `Δ(S5 − AE)` | Pass / Fail |

### Non-Negotiable Reporting Rules
1. **No post-hoc switching**: Regardless of whether H1, H2, or H3 pass or fail, ALL methods (`Constant_Mean`, `kNN_30`, `AE`, `S0`, `S5`, `S7`, `S7_no_S6 [exploratory]`, `Simplified_ResMLP`, `Simplified_GNN`) and ALL metrics (`v5_orth_r`, `delta_r`, `top50_de_r`, `top50_sign`, `disc_raw`, `disc_v5`, plus target-gene-masked `discrim_delta_pct`) MUST appear in the main tables of `ROUND7_CHECKPOINT_1.md`, `ROUND7_EXECUTION_REPORT.md`, and `docs/PREPRINT_DRAFT.md`.
2. **Failure honesty**: If H1 or H2 fails, state **「未确认」** directly in `ROUND7_CHECKPOINT_1.md` and the preprint Abstract; never swap primary method (`S5`), comparator (`AE`), or metric post-hoc.

---

## 4. Locked Split Protocol, Metrics, and Statistical Inference (Exact Commit `3609811` References)

1. **Split Protocol (`make_strict_union_split` in `src/replogle_external_validation.py` @ `3609811`)**:
   - Seeds: `SEEDS_10 = [0, 1, 2, 3, 7, 13, 21, 42, 99, 2024]` (`src/xatlas_panhuman_solver.py` @ `3609811`).
   - Strictly simultaneous union holdout: 20% of unique target genes across the union of both cell lines are held out simultaneously from **both** cell lines (zero cross-line leakage of held-out targets); evaluation is computed on held-out reproducible perturbations, macro-averaged equally (`0.5 * line0 + 0.5 * line1`) across the two cell lines.
   - Training baseline (`base`) and generic stress axis (`axis`) are fit strictly on training rows of each split (`fit_baseline_and_axis`).
2. **Metric Definitions (Exact Match to Round 6 Verified Scripts @ `3609811`)**:
   - **`v5_orth_r`**: `score_prediction` in `src/xatlas_panhuman_solver.py` @ `3609811` (Pearson correlation after QR projection orthogonal to `span{1, baseline, axis}`, excluding perturbed target gene `t_idx`).
   - **`delta_r`**: Standard perturbation Δ-expression Pearson on `s_idx \ {t_idx}` (`src/round6_t2_multi_metric_bench.py` & `src/round6_t2_discrim_bootstrap_ci.py` @ `3609811`).
   - **`top50_de_r` & `top50_sign`**: Pearson correlation and sign accuracy on the top-50 differentially expressed genes ranked by `|y_true - baseline|` (`src/round6_t2_multi_metric_bench.py` @ `3609811`).
   - **`disc_raw` (Primary H1 Definition) & `disc_v5`**: Exact row-level perturbation discrimination percentile `1.0 - (# of j with corr(p_i, t_j) > corr(p_i, t_i)) / (M - 1)` computed in `compute_row_metrics_for_split` (`src/round6_t2_discrim_bootstrap_ci.py` @ `3609811`, verified to 4 decimal places against Master Controller's `replogle_rowlevel.py`). The target-gene-masked variant `discrim_delta_pct` (`src/round6_t2_multi_metric_bench.py` @ `3609811`) will also be reported alongside `disc_raw` for full transparency.
3. **Statistical Inference (`stratified_row_bootstrap_ci` in `src/compute_row_bootstrap_ci.py` @ `3609811`)**:
   - **4,000-resample cell-line-stratified row-level paired bootstrap 95% CI** (`n_boot=4000`, `np.random.RandomState(12345)`, macro-averaged `0.5 * mean(line0) + 0.5 * mean(line1)`).
   - **Split-level paired $t$-test $p$-values are permanently prohibited** because the 10 random splits share overlapping training sets.

---

**Pre-Registration Signed & Locked**: `R7-0` (standalone commit prior to any Nadig et al. 2025 data download or inspection).


---

## Addendum A: Pre-Registration Errata & Metadata-Only Lock Rule (Prior to Data Download)

- **Date**: 2026-11-21 (following Master Controller audit of commit `c5294c0`; original text above preserved unmodified for full auditability)
- **Data-State Audit Note**: Following `c5294c0`, two background `curl` downloads of `GSE264667` were briefly initiated (~10 seconds), immediately cancelled upon receiving Master Controller's Addendum A instruction, and completely deleted (`rm -rf external_data/nadig2025`). Zero Nadig et al. (2025) files exist on disk and zero expression data has been opened or inspected prior to this Addendum A commit.

### Errata (Code-Authoritative Corrections)
1. **Exact Split Seeds (`src/xatlas_panhuman_solver.py:58`)**:
   - Original text erroneously listed `[0, 1, 2, 3, 7, 13, 21, 42, 99, 2024]`.
   - **Authoritative code value**: `SEEDS_10 = [42, 0, 1, 2, 3, 7, 11, 99, 2024, 31337]`. All 10 confirmatory splits will use `SEEDS_10` imported directly from `src/xatlas_panhuman_solver.py`.
2. **Exact Strict Union Holdout Split Definition (`src/replogle_external_validation.py:150-165`)**:
   - Original text loosely described "20% of unique target genes across the union".
   - **Authoritative code behavior (`make_strict_union_split`)**: For each cell line independently (`line0`, `line1`), sample `round(0.20 * len(unique_reproducible_in_readout_targets))` unique target genes as that line's validation targets; then take the **union** of the two cell lines' validation target sets (`union_val_targets`) and simultaneously exclude `union_val_targets` from **both** cell lines' training sets (zero cross-line target leakage).
3. **Exact Stage Naming & Method Signature (`src/xatlas_panhuman_solver.py:282, 543`)**:
   - Original text erroneously referred to `S3` as "Precision Weighting" and the method as `fit_predict(..., stage=5)`.
   - **Authoritative code definitions**:
     - `S3` is **Pan-Human Cross-Line Gene Graph** (`+ Pan-Human Cross-Line Gene Graph` over shared readout genes).
     - The frozen solver method is `PanPerturbBridgeModel.predict_for_source_all_stages(..., stages_to_run=("S0", ..., "S7"))`.
     - Primary method `S5` therefore consists of: `S0` (baseline ridge) + `S1` (HybridTOM $p=1.5$) + `S2` (Multi-Scale Spectral Deflation) + `S3` (Pan-Human Cross-Line Gene Graph) + `S4` (Non-Reproducible Anchors) + `S5` (Cross-Line Perturb-Bridge); excluding `S6` (Reciprocal Borrowing) and `S7` (Output SVD Wiener + Manifold Snapping).

### Additional Lock Rule: Reproducibility Label Definition on Nadig et al. (2025)
4. **Metadata-Only Lock Before Modeling**:
   - If the official Nadig et al. (2025) files do **not** include the exact `anderson_darling_counts` and `energy_test_p_value` columns present in Replogle et al. (2022), the equivalent `reproducible` threshold (and readout-gene selection column) **MUST** be decided strictly by inspecting file schema / `obs` & `var` metadata only — **without inspecting any expression matrix values or running any model** — and **MUST be locked in a dedicated standalone git commit (`Addendum B`) BEFORE building predictions or running any evaluation script**.
