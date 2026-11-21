# Round 8 Pre-Registration & Pre-Declared Adjudication Rules (`R8-1`, `R8-2`, `R8-3`)

**Date**: 2026-11-21  
**Governance**: Locked and committed locally strictly **BEFORE** executing `R8-1`, `R8-2`, `R8-3` exploratory/confirmatory runs, and strictly **BEFORE** downloading any file of any candidate 4th dataset.  
**Frozen Solver Invariant**: `src/xatlas_panhuman_solver.py` remains 100% untouched from commit `3609811`.  
**Remote Push Policy**: Local commits only (zero `git push` to `origin`).

---

## 1. Pre-Declared Rules for `R8-1` — Native Single-Cell GEARS on Nadig et al. (2025)

### 1.1 Protocol
- **Dataset**: Local raw single-cell Nadig et al. (2025) `.h5ad` files (`GSE264667_hepg2_raw_singlecell_01.h5ad`, `GSE264667_jurkat_raw_singlecell_01.h5ad`, SHA-256 already locked in Round 7 Addendum B).
- **Splits**: Exact same `make_strict_union_split` (`SEEDS_10` subset: first 3 seeds **`[42, 0, 1]`**), independent single-cell GEARS model per cell line (`HepG2`, `Jurkat`) per split (zero test leakage; 90%/10% train-internal validation split for GEARS early stopping).
- **Single-Cell Sampling & Normalization**:
  - Readout genes: exact locked 5,000 genes per cell line (`HepG2_genes.csv`, `Jurkat_genes.csv`).
  - Perturbed cells: deterministically sample up to **`64` cells per perturbation** (RNG seeded by split seed). If a single split exceeds **3 hours** on MPS, reduce cap to no lower than **`32` cells per perturbation** and record the exact cap in `ROUND8_EXECUTION_REPORT.md`.
  - Control cells: deterministically sample up to **`5,000` `non-targeting` control cells** per cell line.
  - Normalization consistent with Addendum B: per-`gem_group` library-size normalization to `gem_group` median UMI count, followed by `log1p` (standard non-negative single-cell input required by native `cell-gears==0.1.2` for `non_zeros_gene_idx` and `calc_de`).
  - Official GEARS package (`cell-gears==0.1.2`) default hyperparameters (`hidden_size=64, num_go_gnn_layers=1, num_gene_gnn_layers=1, decoder_hidden_size=16, epochs=20, lr=1e-3, weight_decay=5e-4, batch_size=32`) with **native `calc_de` enabled** (`skip_calc_de=False`).
- **Conversion to Pseudobulk $\Delta$ & Evaluation**:
  - For each held-out perturbation $p$, obtain native GEARS predicted post-perturbation expression $\hat{\mathbf{x}}_{p,\log}$ and convert to parameter-free control-standardized pseudobulk z-score $\hat{\mathbf{y}}_{p,z} = (\hat{\mathbf{x}}_{p,\log} - \boldsymbol{\mu}_{\text{ctrl},\log}) / \boldsymbol{\sigma}_{\text{ctrl},\log}$ (using strictly the 5,000 sampled `non-targeting` control cells; also report unstandardized $\hat{\mathbf{x}}_{p,\log} - \boldsymbol{\mu}_{\text{ctrl},\log}$ as a sensitivity check).
  - Evaluate against the exact locked Nadig 2025 ground-truth cache (`HepG2_X.npy`, `Jurkat_X.npy`) on the identical metric suite (`v5`, `delta_r`, `de50_r`, `de50_sign`, `disc_raw`, `disc_v5`) and 4,000-resample cell-line-stratified row-level paired bootstrap 95% CIs vs `S5` and `AE` (reporting both full and GO-covered-only metrics).

### 1.2 Pre-Declared Decision Rule (`R8-1`)
- **GEARS Underperformance Confirmed**: Native single-cell GEARS still trails both `S5` and `AE` on `v5` and `disc_raw` (`3/3` splits CI `< 0` vs `S5` on `v5` and `disc_raw`, and mean `v5` and `disc_raw` below `AE`) → manuscript statement that official GEARS trails `S5` and `AE` holds even under native single-cell training.
- **Native Single-Cell GEARS Improves Substantially**: If native single-cell GEARS matches or beats `AE` or `S5` on `v5` or `disc_raw` (or closes > 50% of the pseudobulk gap on `v5` / `disc_raw`), report native single-cell GEARS prominently in the main table and revise all manuscript comparisons accordingly.

---

## 2. Pre-Declared Rules for `R8-2` — PC1 Axis vs. Off-Axis Mechanism Decomposition

### 2.1 Definitions (Using Existing Frozen Predictions Only)
Across all three datasets (`xatlas`, `Replogle 2022`, `Nadig 2025`), all 10 splits (`SEEDS_10`), and both cell lines $c \in \{0, 1\}$, let $\mathbf{b}_c$ be the reproducible training perturbation centroid and $\mathbf{a}_c$ ($\|\mathbf{a}_c\|_2 = 1$) be the training PC1 axis (`fit_baseline_and_axis`, training set only). For held-out perturbations $i = 1, \dots, M_c$ with $\mathbf{d}_i = \mathbf{y}_i - \mathbf{b}_c$ and $\hat{\mathbf{d}}_i = \hat{\mathbf{y}}_i - \mathbf{b}_c$:
1. **Axis Component (`axis_r`)**: Project true and predicted $\Delta$ onto $\mathbf{a}_c$: $z_i = \mathbf{d}_i^\top \mathbf{a}_c$, $\hat{z}_i = \hat{\mathbf{d}}_i^\top \mathbf{a}_c$. Compute Pearson correlation $\mathrm{corr}(\{z_i\}, \{\hat{z}_i\})$ across the $M_c$ held-out perturbations (macro-averaged across both lines; 4,000-resample cell-line-stratified row-level paired bootstrap 95% CIs).
2. **Off-Axis Component (`off_axis_r`)**: Project out $\mathbf{a}_c$: $\mathbf{d}_i^\perp = \mathbf{d}_i - z_i \mathbf{a}_c$, $\hat{\mathbf{d}}_i^\perp = \hat{\mathbf{d}}_i - \hat{z}_i \mathbf{a}_c$. Compute per-perturbation Pearson correlation $\mathrm{corr}(\mathbf{d}_i^\perp, \hat{\mathbf{d}}_i^\perp)$ (macro-averaged across both lines; 4,000-resample row-level paired bootstrap 95% CIs).
3. Compare against raw `delta_r` (`= centroid_train_r`) for `AE`, `S0`, `S5`, `S7` on all 3 datasets, plus `scGPT_FrozenEmb_Ridge` on Replogle 2022 and Nadig 2025.

### 2.2 Pre-Declared Decision Rule (`R8-2`)
- **Hypothesis**: The `delta_r` advantage of `AE` (and `scGPT_FrozenEmb_Ridge`) over `S5` arises primarily from the training PC1 axis component.
- **Pass Criterion (per dataset)**: Let $\Delta_{\text{raw}} = \mathrm{mean}(\mathrm{delta\_r}_{\text{AE}} - \mathrm{delta\_r}_{\text{S5}}) > 0$ and $\Delta_{\text{off}} = \mathrm{mean}(\mathrm{off\_axis\_r}_{\text{AE}} - \mathrm{off\_axis\_r}_{\text{S5}})$. The mechanism hypothesis is **CONFIRMED** on a dataset iff:
  $$\Delta_{\text{off}} \le \frac{1}{3}\,\Delta_{\text{raw}} \quad (\text{including sign flip } \Delta_{\text{off}} \le 0).$$
- **Manuscript Rule**: State the PC1-amplitude mechanism as an empirically supported finding in the Abstract **only if** confirmed on all 3 datasets; otherwise keep it qualified as a partial/unsupported mechanism and report the exact decomposition numbers.

---

## 3. Pre-Registration for `R8-3` — Parameter-Free Orthogonal Fusion (`scGPT_Probe ⊕ S5`) & Fourth-Dataset Confirmation

### 3.1 Parameter-Free Orthogonal Decomposition Fusion Rule (ZERO Tunable Weights)
For each split $s \in \mathrm{SEEDS}_{10}$ and each cell line $c \in \{0, 1\}$, let $Q_c \in \mathbb{R}^{G \times k}$ ($k \le 3$) be the orthonormal QR basis of the training shared subspace $S_c = \mathrm{span}\{\mathbf{1}, \mathbf{b}_c, \mathbf{a}_c\}$ (computed strictly from the split's training set, identical to `fit_baseline_and_axis` and `project_rows_5000`), with orthogonal projector $P_{S_c}(\mathbf{u}) = Q_c(Q_c^\top \mathbf{u})$ and orthogonal complement $P_{\perp,c}(\mathbf{u}) = \mathbf{u} - P_{S_c}(\mathbf{u})$.

Given the frozen predictions of **`scGPT_FrozenEmb_Ridge`** ($\hat{\mathbf{y}}_{\text{scGPT}}$) and frozen **`S5`** ($\hat{\mathbf{y}}_{\text{S5}}$), the **Parameter-Free Orthogonal Fusion (`OrthoFusion_scGPT_S5`)** is defined strictly without weights or tuning as:
$$\hat{\mathbf{y}}_{\text{fusion}} = \mathbf{b}_c + P_{S_c}\!\big(\hat{\mathbf{y}}_{\text{scGPT}} - \mathbf{b}_c\big) + P_{\perp,c}\!\big(\hat{\mathbf{y}}_{\text{S5}} - \mathbf{b}_c\big).$$

*No blending weights, norm rescaling, ridge re-fitting, or post-hoc adjustments are permitted.*

### 3.2 Pre-Registered Hypotheses (`H4`, `H5`) — Evaluated with 4,000-Resample Row-Level Paired Bootstrap 95% CIs across `SEEDS_10`
- **H4 (Non-Inferiority on Both Orthogonal Subspaces)**:
  1. `delta_r` (`OrthoFusion_scGPT_S5 − scGPT_FrozenEmb_Ridge`): no significant degradation — **$\le 2/10$ splits with row-level 95% CI $< 0$** (and mean delta $\ge 0$ or within non-inferiority noise); **AND**
  2. `v5` (`OrthoFusion_scGPT_S5 − S5`): no significant degradation — **$\le 2/10$ splits with row-level 95% CI $< 0$**. *(Algebraically guaranteed `0/10` CI `< 0` since $P_{\perp,c}(\hat{\mathbf{y}}_{\text{fusion}}) \equiv P_{\perp,c}(\hat{\mathbf{y}}_{\text{S5}})$.)*
- **H5 (Strict Synergy on Unprojected Perturbation Identification `disc_raw`)**:
  1. `disc_raw` (`OrthoFusion_scGPT_S5 − scGPT_FrozenEmb_Ridge`): **$\ge 7/10$ splits with row-level 95% CI $> 0$** (and mean delta $> 0$); **AND**
  2. `disc_raw` (`OrthoFusion_scGPT_S5 − S5`): **$\ge 7/10$ splits with row-level 95% CI $> 0$** (and mean delta $> 0$).
- **Adjudication**: If either clause of H4 or H5 fails, report **`未确认 (NOT CONFIRMED)`** prominently with exact split counts and CIs.

### 3.3 Exploratory vs. Confirmatory Separation
1. **Exploratory Evaluation (Already-Seen Datasets)**: Run `OrthoFusion_scGPT_S5` on **Replogle 2022** and **Nadig 2025** (`SEEDS_10`) and label these results explicitly as **`EXPLORATORY (seen datasets)`** — never as confirmation.
2. **Confirmatory Evaluation (Unseen 4th Dataset)**: Evaluate on a 4th independent dataset **ONLY** if a candidate passes the strict pre-registered eligibility gate below.

### 3.4 Fourth-Dataset Candidate List, Exclusion Rules, and Locked Eligibility Gate

**Hard Requirements (Non-Negotiable)**:
- Organism: Human (*Homo sapiens*).
- Design: **Two individual human cell lines** from the **same experimental study and matched condition** (no pooling of unrelated donors into synthetic cohorts, no mixing distinct cytokine stimulations as cell-line replicates, no combining two different laboratories' single-line screens).
- Perturbation modality: Single-gene **CRISPRi or CRISPR KO** Perturb-seq.
- Provenance: Publicly downloadable with open license; **distinct experiment** from Replogle et al. (2022) and Nadig et al. (2025).
- **Absolute Prohibition**: **ZERO data from X-Atlas/Orion or any derivative thereof** (benchmark rule).

**Locked Minimum Statistical Power Eligibility Gate** (checked strictly on downloaded metadata / QC counts prior to running any model; **never lowered** to force a dataset):
1. Both selected individual cell lines have **`>= 150` unique target genes with `n_cells >= 25`** (in the same matched experimental condition) and **`>= 100` reproducible target genes** under the locked Replogle/Nadig criterion (`n_cells >= 25`, `ad_counts >= 50` at BH FDR `< 0.05`, 20-PC 2,000-perm energy distance `p < 0.001`).
2. **`>= 80` shared reproducible target genes** observed in both cell lines (required for `S5` cross-line translator fitting and simultaneous union holdout).
3. **`>= 2,000` shared finite readout genes** (selecting up to top `5,000` finite genes per line by `control_mean`, identical to Addendum B).

**Ranked Candidate Queue (Metadata-Only Audit Conducted Prior to Any Download)**:

| Rank | Candidate Study | Source & License | Reported Scope | Pre-Registered Selection Rule & Risk |
|---|---|---|---|---|
| **1 (Primary)** | **Feng et al., *Cell Genomics* 2025** (Wellcome Sanger / EMBL / DKFZ; PMID `41330380`; *"A genome-scale single cell CRISPRi map of trans gene regulation across many human pluripotent stem cell lines"*) | Figshare `10.6084/m9.figshare.27989294.v2` (**MIT License**). Files: `TargetedScreen_Cell-Metadata.tsv.gz` (md5 `980a440a...`), `TargetedScreen_Guide-UMI-Counts.csv.gz` (md5 `96342f75...`), `TargetedScreen_RNA-UMI-Counts.csv.gz` (md5 `f0d613b1...`), `GenomeWideScreen_FitnessGenes_*` (md5 `15423c32...`, `1d49fc95...`, `ba724359...`). | 34 human iPSC lines (26 distinct donors), single-gene CRISPRi (`7,226` genes genome-wide + targeted screen, ~2M cells). | Inspect cell-metadata TSV first (no expression read): pick the **top-2 individual human donor iPSC lines** by total assigned cell count in `TargetedScreen` (or `GenomeWideScreen_FitnessGenes` if `TargetedScreen` fails gate). **Risk**: individual donor lines may have too few cells/perturbation (`n_cells < 25`) across 34 pooled lines; also same cell type (iPSC), different genetic backgrounds. If individual donor lines fail the gate, **DO NOT pool donors into synthetic meta-lines** — reject Candidate 1. |
| **2 (Secondary)** | **Jiang, Dalgarno et al., *Nature Cell Biology* 2025** (Satija Lab; GEO `GSE281048`; *"Systematic reconstruction of molecular pathway signatures using scalable single-cell perturbation screens"*) | Zenodo `10.5281/zenodo.14518762` (**CC-BY-4.0**). Files: `Seurat_object_{IFNB,IFNG,INS,TGFB,TNFA}_Perturb_seq.rds` (md5 recorded in Zenodo API). | 6 distinct human cancer cell lines (`A549`, `BxPC-3`, `HAP1`, `HT-29`, `K562`, `MCF-7`), single-gene CRISPRi (`218` unique signaling genes across 5 cytokine contexts, `1,656` pert-condition pairs). | Inspect metadata in the single largest cytokine context (`INS`, `TNFA`, or `IFNB` by perturbation count): pick the **top-2 distinct cell lines** by perturbation count in that single matched context. **Risk**: each single cytokine file likely targets only ~40–60 genes per line (`1,656 / 30 ≈ 55`), failing the `>= 150` tested / `>= 100` reproducible target gate. If so, **DO NOT pool different cytokine stimulations** — reject Candidate 2. |

**Explicitly Excluded Studies (Not Eligible)**:
- `X-Atlas / Orion` (and any derivative): Strictly prohibited.
- `Replogle et al. (2022)` & `Nadig et al. (2025)`: Already seen (Atlases 1 and 2).
- `Drepanos et al. (2026)` (`Zenodo 20722064`, `K562 + A549`, CC-BY-4.0): Only `25` target genes — fails minimum scale gate.
- `Tian et al. (2019)` (`GSE124703`, iPSC + Day-7 neurons): Small targeted panel (~27 genes) — fails minimum scale gate.
- `Doctor et al. (2026)` (`KOLF2.1J`, Figshare `27261219`), `GSE295214` (hESC), `VIPerturb-seq` (`Zenodo 18460279`): Single cell line per genome-wide screen.
- `GSE320250` (`Figshare 31276981`, HSCs): Mouse (*Mus musculus*), not human.
- `Tahoe-100M` / `sci-Plex`: Chemical small-molecule perturbations, not CRISPRi/KO.

**Fail-Safe Clause**: If **neither** Candidate 1 nor Candidate 2 passes the locked Eligibility Gate above, we **do not** lower standards or force an underpowered/confounded dataset; we document the exact quantitative failure counts for every screened candidate and formally record the 4th-dataset confirmatory verdict as **`待定 (PENDING — no publicly available 4th human 2-cell-line single-gene CRISPRi/KO atlas met pre-registered power criteria)`**. If a candidate passes, we record its exact file SHA-256 checksums and locked preprocessing parameters in `Addendum R8-A` **prior to reading any expression counts**, then execute the exact pre-registered `H4`/`H5` evaluation.

---

## 4. Addendum R8-A — Locked Candidate 1 (`Feng et al. 2025` `TargetedScreen`) Selection & Preprocessing (Committed Prior to Downloading Expression Counts)

**Metadata-Only Screening Audit (Zero Expression Counts Downloaded or Read)**:
- Inspected `GenomeWideScreen_FitnessGenes_Cell-Metadata.tsv.gz`, `GenomeWideScreen_NonFitnessGenes_Cell-Metadata.tsv.gz`, and `TargetedScreen_Cell-Metadata.tsv.gz` from Figshare `10.6084/m9.figshare.27989294.v2` (**MIT License**):
  - `GenomeWideScreen_FitnessGenes`: top-2 lines `wahn_1` (`6` genes $\ge 25$ cells) & `siqu_4` (`2` genes $\ge 25$ cells) → **FAILS** minimum cell-depth gate.
  - `GenomeWideScreen_NonFitnessGenes`: top-2 lines `kolf_2` (`0` genes $\ge 25$ cells) & `miaj_4` (`0` genes $\ge 25$ cells) → **FAILS** minimum cell-depth gate.
  - **`TargetedScreen`** (`1,161,864` cells, `19` human iPSC donor lines): top-2 individual human donor cell lines by assigned cell count are:
    - **Line 0 (`jejf_2`)**: `47,558` assigned cells (`605` `NonTarget` control cells across `Pool4_Day3` [`203`], `Pool5_Day3` [`187`], `Pool7_Day3` [`148`], `Pool6_Day3` [`67`]); `444` unique target genes, **`442` genes with `n_cells >= 25`**.
    - **Line 1 (`zapk_3`)**: `46,376` assigned cells (`622` `NonTarget` control cells across `Pool2_Day3` [`270`], `Pool3_Day3` [`177`], `Pool1_Day3` [`175`]); `444` unique target genes, **`442` genes with `n_cells >= 25`**.
    - **Shared target genes with `n_cells >= 25`**: **`442`** (completely disjoint experimental pools: `Pool4–7` vs `Pool1–3`, zero batch overlap).
    - *Transparent scope note*: Both `jejf_2` and `zapk_3` are human iPSC cell lines from two unrelated healthy HipSci donors (distinct genetic backgrounds and distinct experimental pools, though same pluripotent cell state rather than two distinct tissue lineages; this limitation will be explicitly stated in the manuscript).

**Locked Official Files & Published MD5 Checksums (Before Expression Download)**:
- `TargetedScreen_Cell-Metadata.tsv.gz` (`13.5 MB`, Figshare MD5 `980a440a8132a20567c4c3610e6a5b20`)
- `TargetedScreen_RNA-UMI-Counts.csv.gz` (`4,481.0 MB`, Figshare MD5 `f0d613b1ea8147425501baddba034158`)
- Exact SHA-256 of both downloaded files will be verified and logged in `prepare_feng2025_cache.py` and `ROUND8_EXECUTION_REPORT.md`.

**Locked Preprocessing Protocol (Zero Deviations Permitted)**:
1. **Cell Filtering**: Retain cells belonging to `Cell_Line in {"jejf_2", "zapk_3"}` whose `Guide_Call` is a single assigned target (`gene = Guide_Call.split("_")[0]`, excluding `unassigned`, `multiplet`, `nan`).
2. **Per-`Pool_Day` Normalization & Control Z-Scoring** (exact analog of Addendum B `gem_group`, using `Pool_Day` because individual `Inlet`s contain as few as 2 `NonTarget` cells whereas every `Pool_Day` contains `67–270` `NonTarget` cells):
   - Extract `Pool_Day` via regex `Pool\d+_Day\d+` from `Batch`.
   - Within each `(Cell_Line, Pool_Day)` stratum, library-size normalize every cell's raw UMI counts to that stratum's median total UMI count across retained cells.
   - Z-score every gene within each `(Cell_Line, Pool_Day)` stratum against cells with `gene == "NonTarget"` in that same stratum.
   - Average z-scored profiles across all retained cells targeting each unique gene symbol (merging guides targeting the same gene) to form the pseudobulk profile per `(Cell_Line, target_gene)`.
3. **Readout Gene Selection (Up to Top `5,000` Finite Genes per Line)**:
   - Drop any gene with zero control variance or non-finite z-score in any `Pool_Day` stratum of that line.
   - Rank remaining finite genes by mean library-size-normalized expression across all `NonTarget` control cells (`control_mean`) and select the top `5,000` genes per line (or all finite genes if `< 5,000`).
4. **Locked Reproducibility Criterion (Identical to Replogle 2022 & Nadig 2025 Addendum B)**:
   - Vectorized Scholz-Stephens two-sample Anderson-Darling test across all selected readout genes against all `NonTarget` control cells of that line (`605` in `jejf_2`, `622` in `zapk_3`; BH FDR `< 0.05` → `ad_counts`).
   - 20-PC 2,000-permutation energy-distance test against all `NonTarget` control cells of that line (`energy_test_p_value < 0.001`).
   - `reproducible = 1` iff `n_cells >= 25` AND `ad_counts >= 50` AND `energy_test_p_value < 0.001`.
5. **Post-QC Eligibility Gate (Hard Stop Before Running Any Model)**:
   - Proceed to `H4`/`H5` confirmatory testing **ONLY** if both `jejf_2` and `zapk_3` have `>= 100` reproducible target genes, `>= 80` shared reproducible target genes, and `>= 2,000` shared readout genes. If this gate fails, **do not** relax thresholds; mark 4th-dataset confirmation **`待定 (PENDING)`**.

