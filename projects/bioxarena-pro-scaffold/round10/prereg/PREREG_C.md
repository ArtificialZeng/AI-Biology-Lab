# Round 10 Confirmatory Pre-Registration Protocol (`PREREG_C.md`, v1 R1-Revised)

- **Project:** AI-Biology-Lab Round 10 — Prospective Zero-Shot Verification of Amplitude-Orthogonalized Fold-Change Transfer (`FCT`) on Completely Unseen Human CRISPRi Cell Lines
- **Protocol Author:** Agent `A-T3` (Team A · Pre-registration & Statistical Audit), per Master Controller Directive [`R10_第十轮行动指令.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/controller/directives/R10_第十轮行动指令.md) (`R10-1`, `2026-11-22`)
- **Upstream Context & Blinding Basis:**
  - Round 9 Confirmatory Verdict (`PREREG_A v4`, commit `4cc89ebae67f0ed0748b80f994d26dcc9c2b95ef`; `CODE_FREEZE: 1ca9f60f49adf1ac7cce8266e40533273411db30`): Both Group A (`Hybrid_ProjShared`) and Group B (`Transfer_Mean_VarScaled` / `TMVS`) were **NOT CONFIRMED (`未确认`)**.
  - Post-hoc Exploratory Clue ([`EXPLORATORY_FCT.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2/vcc_linespecific/results/EXPLORATORY_FCT.md)): On Round 9's secondary `FCT` baseline, post-hoc specificity testing observed `r_orth_cal = +0.2611` (`[+0.2420, +0.2807]`) on `nadig_hepg2` (`True`) and `+0.3505` (`[+0.2784, +0.4245]`) on `vcc2025_h1` (`True`), while failing permutation nulls on `nadig_jurkat` (`+0.1043` `[+0.0849, +0.1238] <= +0.1201`). Because that check followed unblinding of Round 9 held-out expression matrices, Round 9 cannot treat `FCT` as confirmatory.
  - Literature / Competition Audit ([`PRIOR_ART_C.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2/vcc_linespecific/audit/LIT/PRIOR_ART_C.md)): Base `FCT` (`mu_t * FC`) already exists in VCC 2026 #1 `Illumina/BioInsight CellularAI` and VCC 2026 #82 `AtlasShift` (**not a novel method**). However, Q1 confirmed **`没有找到`** any study testing whether an `FCT`-class method recovers **de-amplituded (`beta_t`-removed) line-specific gene patterns (`r_orth_cal`) above mismatched target control means and gene-label permutations**.
  - Zero-Data-Access Guarantee: As of drafting `PREREG_C.md`, **zero Round 10 expression matrices have been downloaded or inspected** (only public API/supplementary metadata & guide-target symbol lists were screened in `R10-0`, [`R10_DATA_SCREEN.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2/shared/R10_DATA_SCREEN.md)).

---

## 1. Strict Six-Step Execution Gate (`R10-1` Order Lock)

No downstream step may commence until the preceding gate file exists in `round10/status/`:

1. **Step 1 (Pre-Registration Lock):** `A-T3` drafts `round10/prereg/PREREG_C.md` and submits `round10/status/COMMIT_REQUEST_A-T3.txt`.
2. **Step 2 (`R10_GO.txt` Issuance):** Master Controller (`Claude`) audits `PREREG_C.md`, locks its git commit hash into `round10/status/PREREG_C_COMMIT.txt` (and `PREREG_C_LOCK_PUSHED.txt`), and writes `round10/status/R10_GO.txt` together with `round10/shared/R10_DOWNLOAD_PLAN.md` sign-off. **Strict Prohibition:** Zero new Round 10 data files may be downloaded before `round10/status/R10_GO.txt` exists.
3. **Step 3 (`A-T1` Download & Counts-Only QC):** `A-T1` downloads candidate datasets strictly within the `60.000 GB` Round 10 quota, applies the pre-locked normalization, Addendum-B QC, and split-half partition, and emits **counts-only metadata** (`n_cells`, `n_nt_cells`, `n_genes_finite`, `n_targets_dedup`, `n_targets_reproducible`, `n_shared_reproducible_with_r9`, `half_valid` counts; **zero** prediction metrics or correlations computed). Meanwhile `A-T2` and `A-T3` implement and unit-test `round10/src/` and `round10/audit/` strictly on Round 9's six already-unblinded lines (`round9/p2/shared/lines/`).
4. **Step 4 (Mechanical Confirmatory Line Determination):** Apply §2.2's pre-locked threshold (`>= 50` reproducible shared targets with Round 9 sources) mechanically from `A-T1`'s counts-only table and record eligible lines in `round10/status/CONFIRMATORY_LINES.txt` without human post-hoc overrides.
5. **Step 5 (Code Freeze):** Master Controller verifies code, `PREREG_C` parity, and unit tests, locking commit hash in `round10/status/CODE_FREEZE_R10.txt`. Runtime gate (`check_prereg_gate`) refuses any Round 10 evaluation unless both `PREREG_C_COMMIT.txt` and `CODE_FREEZE_R10.txt` exist and are non-empty.
6. **Step 6 (Unblinding & Independent Dual-Team Evaluation):** `A-T2` executes primary Round 10 evaluation into `round10/results/` and `A-T3` executes completely independent recomputation + `C3` cross-check into `round10/audit/`. Unfavorable results are reported verbatim.

---

## 2. Candidate Datasets, Eligibility Rules, & Fixed Split-Half Partitions

### 2.1 Source Cohort (Fixed, Zero Leakage)

All Round 10 predictors take their training source pool strictly from all six harmonized Round 9 human CRISPRi lines (`round9/p2/shared/lines/`, `6,148` aligned common genes `COMMON_GENES.txt`, `2,655` `union6` deduplicated target genes):
1. `nadig_hepg2` (`2,343` targets)
2. `nadig_jurkat` (`2,350` targets)
3. `vcc2025_h1` (`300` targets)
4. `replogle_k562_essential` (`2,049` targets)
5. `replogle_rpe1` (`2,351` targets)
6. `nourreddine_kolf21j` (`95` targets)

Because Round 10 evaluates strictly **brand-new, distinct human cell lines**, zero target-line perturbation data ever appear in the source pool $\mathcal{S} = \{1,\dots,6\}$.

### 2.2 Confirmatory Candidate Pool & Eligibility Rule (Locked Prior to Download)

Per `R10_第十轮行动指令.md` (`R10-0 主控复查结论`, `R10-1`, and `PREREG_C_AUDIT_R1.md` `P1–P5`) and [`R10_DATA_SCREEN.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2/shared/R10_DATA_SCREEN.md), candidate confirmatory lines, single-sgRNA cell rules, pre-download exclusions, and target eligibility are fixed **before** downloading or accessing any Round 10 expression values:

- **Single-sgRNA / Low-MOI Cell Filtering Lock (`P4`):**
  - **`dld1`:** Strictly **Low-MOI** library ONLY (`GSE337988_sublib2_processed_objects_Low_assays.h5`, excluding Med/High MOI), retaining **only perturbed cells uniquely assigned to exactly `1` sgRNA** (plus QC-valid non-targeting control cells).
  - **`cd4t_rest`, `wtc11_i3neuron`, and `calu3`:** Analogously retain **only cells uniquely assigned to `1` targeting sgRNA** (single-gene elements only on `calu3`, excluding dual-gene pairs) and valid NTC control cells.
- **Candidate Confirmatory Pool (Pre-Specified Candidates):**
  1. **`dld1`** (Human colorectal adenocarcinoma epithelial, `Zim3-dCas9` CRISPRi, GEO `GSE337988`, Low-MOI `Low_assays.h5` + guide reference, `4.621–5.291 GB`, `504` NTC guides, `4,954` target genes, `668` raw shared with R9 `union6`).
  2. **`cd4t_rest` (`P3` Pre-Locked 2-Donor Subset):** Primary normal human CD4+ $\alpha\beta$ T cells resting state, `ZIM3KRAB-dCas9` CRISPRi, GEO `GSE314342` `Rest` single-cell `.h5` lanes strictly from the **two lexicographically earliest human donors in GEO (`CE0006864` and `CE0008162`)**, subject to `A-T1` confirming exact GEO accession/lane IDs in `R10_DOWNLOAD_PLAN.md` before `R10_GO.txt` (~`27.3 GB` combined; exact 2 donors locked pre-download so split halves partition strictly by donor).
  3. **`h9_neuron` (`P2` Pre-Download Metadata-Gated Candidate):** `WA09` H9 female hESC-derived Day-6 immature cortical neurons (`IGVFDS9825FOBP`, split `.h5`+`.pkl`+`.csv` = `50.999 GB`, `MOI ~ 4` with `~5%` single-sgRNA cells). **Mandatory pre-download decision rule (`zero expression data touched`):**
     - Before downloading the `51 GB` H9 expression/guide matrix, `A-T1` inspects only the lightweight IGVF guide/cell metadata table (`IGVFFI7085MUTI.csv.gz`, `0.24 MB` or file metadata) to count single-sgRNA cells per target (strictly a QC count, NOT an expression outcome).
     - **Condition (a):** If `< 50` targets shared with Round 9 sources have **`>= 25` single-sgRNA cells**, `h9_neuron` is **excluded pre-download and NOT downloaded**.
     - **Condition (b):** Even if `>= 50` shared targets reach `>= 25` single-sgRNA cells, the 2-donor `cd4t_rest` + `dld1` + `i3Neurons` + `Calu-3` core consumes `~33–34 GB`; adding `h9_neuron` (`51.00 GB`) would total `~85 GB`, exceeding the **`60.000 GB` Round 10 download cap**. Under that budget overrun, `h9_neuron` is likewise **excluded pre-download and NOT downloaded**, and the Round 10 manuscript explicitly states: *"H9 cortical neurons were excluded pre-download due to the 60 GB download budget (因预算排除)."*
- **Mandatory Secondary-Only Cohort (Never Enter Confirmatory `K_conf`):**
  - **`wtc11_i3neuron`** (`TianKampmann2021_CRISPRi.h5ad`, Zenodo `13350497` / `GSE152988`, `0.289 GB`, `70` raw shared targets, CRISPRi single-sgRNA subset only) -> **Secondary report only** (`R10-1`).
  - **`calu3`** (`SunshineHein2023.h5ad`, Zenodo `13350497` / `GSE208240`, `0.744 GB`, `63` raw shared single-gene targets, single-gene elements only, SARS-CoV-2 challenged bystander paracrine state) -> **Secondary report only** (`R10-1`).
- **Complete Evaluated Target Eligibility (`P5`, `only_reproducible = True`):**
  Matches Round 9's primary pipeline verbatim: for any target line $t$, an evaluated target gene $p$ MUST satisfy **all four criteria**:
  1. Reproducible on target line $t$ under Addendum B (`n_cells >= 25` single-sgRNA cells, `>= 50` top-5,000 genes with Anderson–Darling `BH FDR < 0.05`, and `20-PC` energy-distance permutation `p < 0.001`),
  2. Available AND reproducible under Addendum B on **at least 1** Round 9 source line ($|\mathcal{S}_p^{\mathrm{rep}}| \ge 1$),
  3. Target gene symbol $g_p$ belongs to aligned evaluation gene set $G_{\mathrm{eval}}(t)$ (while target self-gene $g_p$ itself is masked out of correlation/residuals during scoring), and
  4. Non-empty overlap in gene space (`n_eval_targets` counted after criteria 1–3).
- **Mechanical `K_conf` Assignment & `< 2` Confirmatory Abortion Rule (`P1`):**
  - Qualifying lines from `{dld1, cd4t_rest[, h9_neuron if ever downloaded]}` achieving **`n_eval_targets >= 50`** enter Confirmatory Cohort `K_conf`.
  - **If `K_conf >= 2`:** All confirmatory claims require **unanimous passage (`K_conf / K_conf`, 100%)** (`3/3` if `K_conf = 3`; `2/2` if `K_conf = 2`). Whenever `K_conf = 2`, abstracts and summary tables **must explicitly state `"K_conf = 2 confirmatory cell lines"`**.
  - **If `K_conf < 2` (`P1`):** Round 10 **performs NO confirmatory hypothesis adjudication (`不做确认性判定`)**, reports strictly **descriptive results only**, and states verbatim in the Abstract: *"Insufficient confirmatory cell lines (`K_conf < 2`); confirmatory hypothesis testing was not performed (`确认性线不足，未做确认性检验`)."*

### 2.3 Locked Split-Half Partitions (`R10-1` & `P3` Exact Batch Definitions)

To eliminate shared-channel/shared-donor noise leakage between Half A and Half B noise-ceiling estimates, split halves partition **disjoint biological/experimental units** rather than mixing cells within the same unit:

| Target Line | Experimental Unit Partitioned | Exact Deterministic Assignment (`A` vs `B`) |
| :--- | :--- | :--- |
| **`dld1`** | **10x Chromium Channels** (`24` Low-MOI channels) | Sort distinct channel IDs lexicographically; even indices `[0, 2, 4, ...]` $\to$ Half A; odd indices `[1, 3, 5, ...]` $\to$ Half B. Within each half, subtract channel-matched NT control mean $\mu_{\mathrm{ctrl}, \mathrm{channel}}$ before pooling single-sgRNA perturbed cells inside that half (valid if perturbation has `>= 2` cells in Half A AND `>= 2` cells in Half B, and each contributing channel has `>= 15` NT cells; fallback to half-pooled NT mean/sd if a channel has `< 15` NT cells). |
| **`cd4t_rest`** | **Human Donors (`CE0006864`, `CE0008162`)** | Pre-locked 2-donor split: **`Half A = CE0006864`**, **`Half B = CE0008162`** (within each half, subtract 10x processing-lane NT control mean using strictly that donor's NT cells). **Because halves span distinct human donors, inter-donor biological differences are incorporated directly into the split-half noise ceiling (`供体间的生物学差异会算进噪声天花板`).** |
| **`h9_neuron`** *(if ever downloaded)* | **HTO Superloaded Libraries** (`24` libraries `LW558..LW581`) | Sort distinct HTO library IDs lexicographically; even indices $\to$ Half A (`12` libs), odd indices $\to$ Half B (`12` libs), using strictly each half's own library NT controls. |
| `wtc11_i3neuron` *(Sec.)* | **10x Chromium Lanes** (`lib1..lib4`) | Sort lane IDs lexicographically; `{even}` $\to$ Half A, `{odd}` $\to$ Half B. |
| `calu3` *(Sec.)* | **10x Channel / Replicate Batch** | Sort batch/lane `.obs` IDs lexicographically `{even}` vs `{odd}` if multi-batch metadata exists; otherwise deterministic stratified 50/50 cell split (`seed = 42`, explicitly flagged). |

### 2.4 Standardized Expression Normalization, Gene Alignment, & Addendum-B QC

- **Expression Representations:** Identical to `round9/p2/shared/scripts/prepare_lines.py`:
  1. Library-size normalize raw UMI counts per cell to sum $10^4$ (`X_libnorm`).
  2. Transform via $\operatorname{log1p}(X_{\mathrm{libnorm}})$; subtract unit-matched (channel/lane/batch) NT control `log1p` mean $\mu_{t,\mathrm{batch}}$ per perturbed cell prior to averaging per-target pseudobulk delta $\Delta^{\mathrm{obs}}_{t,p}$ (`delta_log1p`).
  3. Compute unperturbed target control libnorm mean $\mu_{t,g}$ across batch-stratified sampled NT controls (subsampled up to `20,000` cells matching Round 9, storing full NT batch means across all QC-passed NT cells for $\Delta^{\mathrm{obs}}$ centering and `FCT` $\mu_t$; sensitivity check compares full-NT vs `20k`-NT $\mu_t$ parity) and store $C_t$ (`X_log1p` covariance) for Group A/B comparator baselines.
- **Gene Alignment (No Leakage):** Target gene symbols are mapped to canonical HGNC (`Ensembl` cross-checked; drop empty/duplicate/non-finite-across-batches symbols). For each target line $t$, intersect target finite HGNC genes $G_t^{\mathrm{fin}}$ with Round 9's fixed `6,148` `COMMON_GENES.txt`, forming evaluation gene set $G_{\mathrm{eval}}(t) = G_{\mathrm{R9\_6148}} \cap G_t^{\mathrm{fin}}$ (`require >= 3,000 genes`; target self-gene excluded from evaluation masks). (Sensitivity/baseline comparators use only genes present in $G_{\mathrm{eval}}(t)$, where all 6 Round 9 source lines already cover every gene in `COMMON_GENES.txt`.)
- **Addendum-B Reproducibility Filter:** Applied identically to Round 9 (`n_cells >= 25`, `>= 50` top-5,000 finite genes by target control mean passing vectorized Scholz–Stephens two-sample Anderson–Darling `BH FDR < 0.05` vs `1,000` stratified NT cells, `20-PC` energy-distance permutation `p < 0.001` over `2,000` permutations vs `500` NT cells; deduplicate 1 row per target gene symbol prioritizing reproducible, AD count, `n_cells`).

---

## 3. Predictors & Pre-Locked Model Specification (`Group C = FCT`)

### 3.1 Primary Confirmatory Predictor (`Group C`: Fold-Change Transfer `FCT`)

Strictly **unchanged from `PREREG_A v4` §4.6** (`round9/p2/vcc_linespecific/src/predictors.py:382-409`), zero tuning on any Round 10 or held-out data:

For each evaluated shared target $p$ (observed in source subset $\mathcal{S}_p \subseteq \{1,\dots,6\}$):
$$\mu_{s,g} = \max(\mu^{\mathrm{libnorm}}_{s,g},\, 10^{-8}), \qquad \hat{p}_{s,g} = \max\!\big(\operatorname{expm1}(\operatorname{log1p}(\mu_{s,g}) + \Delta^{\mathrm{obs}}_{s,p,g}),\, 10^{-8}\big),$$
$$\overline{\mathrm{FC}}_{p,g} = \exp\!\left(\frac{1}{|\mathcal{S}_p|} \sum_{s \in \mathcal{S}_p} \ln\!\left(\frac{\hat{p}_{s,g}}{\mu_{s,g}}\right)\right), \qquad \hat{\Delta}^{\mathrm{FCT}}_{t,p,g} = \operatorname{log1p}\!\big(\mu_{t,g}\, \overline{\mathrm{FC}}_{p,g}\big) - \operatorname{log1p}(\mu_{t,g}).$$

**Methodological Attribution Lock (`PRIOR_ART_B` / `C`):** `FCT` implements target baseline mean-expression multiplicative fold-change scaling, matching the formulation class deployed in VCC 2026 #1 (`CellularAI`) and VCC 2026 #82 (`AtlasShift`). No algorithmic novelty is claimed for the `FCT` predictor; Round 10 tests solely the pre-registered empirical question of whether `FCT` achieves **post-amplitude (`beta_t`-removed) line-specific gene-pattern specificity** on independent human cell lines.

### 3.2 Mandatory Negative & Comparator Baselines

1. **Primary `Transfer_Mean`:** Unscaled cross-source average $\hat{\Delta}^{\mathrm{TMean}}_p = \frac{1}{|\mathcal{S}_p|}\sum_{s\in \mathcal{S}_p} \Delta^{\mathrm{obs}}_{s,p}$ (anchor for defining residuals $r^{\mathrm{obs}}_p = \Delta^{\mathrm{obs}}_{t,p} - \hat{\Delta}^{\mathrm{TMean}}_p$, $r^{\mathrm{pred}}_p = \hat{\Delta}_{t,p} - \hat{\Delta}^{\mathrm{TMean}}_p$, and `H-C3` non-inferiority).
2. **`NegCtrl_FCT_MismatchedMean(s')` (`6` controls):** For each of the `6` Round 9 source lines $s' \in \mathcal{S}$, replace $\mu_{t}$ with source line $s'$'s libnorm control mean $\mu_{s'}$ inside `FCT`:
   $$\hat{\Delta}^{\mathrm{FCT\_Mis}(s')}_{t,p,g} = \operatorname{log1p}(\mu_{s',g}\overline{\mathrm{FC}}_{p,g}) - \operatorname{log1p}(\mu_{s',g}).$$
3. **`NegCtrl_FCT_Perm(k)` (`100` permutations, `k = 0..99`):** Permute target control mean $\mu_t$ across genes $g \in \{1,\dots,G_{\mathrm{eval}}(t)\}$ using `np.random.default_rng(k).permutation(G)` while keeping $\overline{\mathrm{FC}}_p$ unpermuted, recomputing `FCT` residual predictions for each `k = 0..99` (`C2` parity).
4. **R18 Global-Amplitude Control (`Transfer_Mean_GlobalScale`):** Single per-line scalar variance-ratio rescaling $c_{\mathrm{global}} \hat{\Delta}^{\mathrm{TMean}}_p$ to verify whether raw un-orthogonalized `r_cal` remains susceptible to whole-line scalar artifacts.
5. **Round 9 Anchors (Secondary Benchmark Context, `SA-3`):** `Hybrid_ProjShared` (Group A), `Transfer_Mean_VarScaled` (`TMVS`, Group B), `Hybrid_on_TMVS` (`R14`), `CIPHER_Linear`, and `Source_Grand_Mean`.

---

## 4. Endpoints, Statistical Inference, & Adjudication Rules

### 4.1 De-Amplituded Primary Endpoint (`r_orth_cal`)

Identical to `PREREG_A v4` (`R17`):
- For each perturbation $p$, excluding self-target gene $g_p$, form evaluation gene set (`resp_top50`: top 50 finite genes in $G_{\mathrm{eval}}(t) \setminus \{g_p\}$ ranked by $|\Delta^{\mathrm{obs}}_{t,p,g}|$).
- Orthogonally project both true residual $r^{\mathrm{obs}}_p$ and predicted residual $r^{\mathrm{pred}}_p$ perpendicular to `Transfer_Mean` ($\hat{\Delta}^{\mathrm{TMean}}_p$) within that perturbation's 50-gene masked vector:
  $$r^{\mathrm{obs}}_{\perp, p} = r^{\mathrm{obs}}_p - \frac{\langle r^{\mathrm{obs}}_p, \hat{\Delta}^{\mathrm{TMean}}_p\rangle}{\|\hat{\Delta}^{\mathrm{TMean}}_p\|_2^2}\hat{\Delta}^{\mathrm{TMean}}_p, \qquad r^{\mathrm{pred}}_{\perp, p} = r^{\mathrm{pred}}_p - \frac{\langle r^{\mathrm{pred}}_p, \hat{\Delta}^{\mathrm{TMean}}_p\rangle}{\|\hat{\Delta}^{\mathrm{TMean}}_p\|_2^2}\hat{\Delta}^{\mathrm{TMean}}_p.$$
- Pool across all evaluated perturbations to compute Pearson correlation $r^{\mathrm{raw}}_{\perp}$ and disattenuate by the orthogonalized split-half Spearman–Brown noise ceiling $c^{\perp}_{\mathrm{ceil}} = \sqrt{2\rho^{\perp}_{AB}/(1+\rho^{\perp}_{AB})}$ computed on the pre-locked channel/donor/HTO split halves ($\hat{r}^{\mathrm{orth\_cal}} = \operatorname{clip}(r^{\mathrm{raw}}_{\perp} / c^{\perp}_{\mathrm{ceil}}, -1, 1)$ whenever $c^{\perp}_{\mathrm{ceil}} > 0$).
- **Joint Perturbation-Level Bootstrap (`C2`):** `n_boot = 4,000` resamples over perturbations using `np.random.default_rng(20261122)`. Crucially, **both $r^{\mathrm{raw}}_{\perp}$ AND $c^{\perp}_{\mathrm{ceil}}$ are recomputed jointly inside each bootstrap draw** to incorporate finite-sample uncertainty of the noise ceiling directly into `r_orth_cal` intervals.
- **Bonferroni Family Alpha:** Two co-primary claims (`H-C1`, `H-C2`) + parallel usability (`H-C3`) use **two-sided `97.5%` percentile CIs** (`[1.25%, 98.75%]`, matching Round 9 `R13`).

### 4.2 Confirmatory Decision Rules (Unanimous Across All `K_conf` Lines)

#### Confirmatory Cohort Threshold (`P1`)
Confirmatory hypothesis testing is conducted **only if `K_conf >= 2`**. If `K_conf < 2`, Round 10 **performs NO confirmatory adjudication (`不做确认性判定`)**, reports strictly **descriptive metrics only**, and explicitly states in the Abstract: *"Insufficient confirmatory cell lines (`K_conf < 2`); confirmatory hypothesis testing was not performed (`确认性线不足，未做确认性检验`)."*

#### Primary Scientific Claim: *Target-Line-Specific Gene-Level Perturbation Signal (`FCT`)*
When `K_conf >= 2`, established **if and only if BOTH `H-C1` AND `H-C2` pass on 100% (`K_conf / K_conf`) of confirmatory lines**; otherwise declared **`未确认 (NOT CONFIRMED)`**:
- **`H-C1` (Positive De-Amplituded Signal):** `FCT` `r_orth_cal` joint bootstrap `97.5% CI_lo > 0`.
- **`H-C2` (Specificity Over Mismatched Controls & Permutations):** `FCT` `r_orth_cal` `97.5% CI_lo` strictly exceeds BOTH:
  1. The maximum upper `97.5%` CI bound across all `6` Round 9 source mismatched control means (`max_{s'=1..6} NegCtrl_FCT_MismatchedMean(s') CI97.5_hi`), AND
  2. The 97.5th percentile across `100` target-control gene-label permutations (`Q_{0.975}(NegCtrl_FCT_Perm(0..99))`).

#### Parallel Practical Utility Verdict: *Non-Inferiority Relative to `Transfer_Mean` (`H-C3`)*
Reported **independently and co-featured in the Abstract** alongside (`H-C1`, `H-C2`):
- **`H-C3` (`AND` across all 6 aggregate metrics):** For paired differences (`FCT` minus `Transfer_Mean`) with `4,000`-resample `97.5%` CIs across `(1) ctrl_disc_panel`, `(2) ctrl_disc_ownmasked`, `(3) ctrl_de50_r`, `(4) cent_disc_panel` (relative to `Source_Grand_Mean`), `(5) cent_disc_ownmasked`, and `(6) cent_de50_r`:
  - **Non-inferior (`H-C3 PASS`):** Across **all 6** conventions on **all `K_conf`** confirmatory lines, no `97.5%` CI lies strictly below zero (`diff_ci975_hi >= 0` everywhere; i.e., `0` inferior conventions on `0` confirmatory lines).
  - **Strict superiority (additional descriptive flag):** Count of conventions/lines with `diff_ci975_lo > 0`. (`pds_cosine_panel` also reported descriptively.)

---

## 5. Mandatory Secondary & Sensitivity Analyses (`SA-1`–`SA-6` Including `S7` Cell-Count Controls)

All secondary analyses below are pre-specified (`non-confirmatory` diagnostics):
- **`SA-1` (Slope vs. `0.13` Descriptive Check):** Report `FCT` `slope_ols` (`[97.5% CI]`), `sd_pred / sd_obs`, and `slope_deming` against Arc VCC 2026 Wrap-Up `[0.02, 0.13]`, retaining the `R6` limitation note that Arc's slope estimator is unpublished in `cell-eval2 v0.22.0`.
- **`SA-2` (Raw `r_cal` & `R18 GlobalScale` Callout):** Report raw un-orthogonalized `r_cal` and `Transfer_Mean_GlobalScale` alongside `r_orth_cal` (`0.0000` collapse check).
- **`SA-3` (Round 9 Benchmark Anchors — `P7` Code-Freeze Definition):** Evaluate `Hybrid_ProjShared`, `Transfer_Mean_VarScaled` (`TMVS`), `Hybrid_on_TMVS` (`R14`), `CIPHER_Linear`, and `Source_Grand_Mean` cleanly on the new lines for direct cross-round comparability. Specifically, `Hybrid_ProjShared` follows Round 9's frozen code definition (`commit 1ca9f60f`, `src/predictors.py::fit_shared_subspace`), projecting out the **3-dimensional source subspace $\operatorname{span}\{\mathbf{1}, \bar{b}_S, a_S\}$** (noting explicitly, as documented in Round 9 `PAPER_REVIEW_R1.md` B1 / `DEVIATIONS.md`, that `PREREG_A v4` prose described rank-1 `q1` rather than the frozen 3D implementation).
- **`SA-4` (`R19` 3-Tier Multiplicative Gain Decomposition):** Partition $\mathrm{SST}(r^{\mathrm{obs}})$ into whole-line (`L1` $\beta_t$), per-gene incremental (`L2` $m_{t,g}$ for genes observed in $\ge 5$ distinct perturbations with cluster-df bias correction), and unexplained residual (`L3`) on both `resp_top50` and `all_nontarget`.
- **`SA-5` (`S6` All-Non-Target, `S5` Top-75% Genes, & `S1` All-Targets):** Quantify `resp_top50` ceiling upward selection bias (`all_nontarget` $\rho^{\perp}_{AB}$ and `r_orth_cal`), robustness on top-75% expressed/variable genes (`S5`), and inclusion of non-reproducible perturbations (`S1`, `--all-targets`). Also check sensitivity of `FCT` if $\mu_t$ uses all QC-passed NT cells vs `20k`-stratified NT subsample (`SA-5b`).
- **`SA-6` / `S7` (Cell-Count & Knockdown-Depth Confounding Controls — `arXiv:2608.00152v3`, `P6`):** *Strictly secondary analyses (`全部是次要分析`)* following Shoeibi & Yousefi (2026):
  1. **Per-Perturbation Orthogonalized Agreement (`r_orth_p`):** For each evaluated perturbation $p$ individually, compute `r_orth_p` as the Pearson correlation between $r^{\mathrm{obs}}_{\perp, p}$ and $r^{\mathrm{pred}}_{\perp, p}$ across perturbation $p$'s masked 50-gene vector (`resp_top50`), with knockdown depth $|k_p| = |\Delta^{\mathrm{obs}}_{t,p,g_p}|$ (or self-gene availability flagged if $g_p$ absent) and single-sgRNA cell count `n_cells_p`.
  2. **Rank Correlation with Cell Count and Knockdown Depth:** Report **Spearman correlation** (with `97.5%` perturbation-level bootstrap CIs) across evaluated perturbations between `r_orth_p` and $\ln(n_{\mathrm{cells},p})$, and between `r_orth_p` and $|k_p|$.
  3. **Covariate-Adjusted OLS Intercept (`[1, ln(n_cells), |k_p|]`):** Fit OLS of `r_orth_p` on mean-centered covariates `[1, ln(n_cells_p) - mean(ln(n_cells)), |k_p| - mean(|k|)]` and report the **intercept** (the average orthogonalized `FCT` signal after controlling for perturbation cell count and knockdown depth) with its `97.5%` bootstrap CI (alongside the regression slopes on `ln(n_cells)` and `|k_p|`, noting invariance of slopes and exact centering interpretation relative to raw uncentered `[1, ln(n_cells), |k_p|]`).
  4. **`n_cells` Tertile Stratification:** Within each confirmatory line, split evaluated targets into three balanced tertiles by per-perturbation `n_cells` (`T1` lowest, `T2` middle, `T3` highest), recomputing pooled `FCT` `r_orth_cal [97.5% CI]` and mismatched/permutation comparisons within each tertile to examine whether signal persists in low-depth perturbations (`T1`).

---

## 6. Software Governance, Code Freeze (`C1`), Cross-Check (`C3`), & Reporting Rules

- **Strict Blinding (`C1`):** `round10/src/` (implemented by `A-T2`) and `round10/audit/` (implemented independently by `A-T3`) may ONLY be developed/tested using Round 9's six harmonized lines (`round9/p2/shared/lines/`) until `round10/status/CODE_FREEZE_R10.txt` is locked. Any post-freeze bug fix requires a public `round10/prereg/DEVIATIONS_R10.md` amendment approved by Master Controller before re-running.
- **Deterministic RNGs (`C2`):** `A-T2` and `A-T3` must both use `np.random.default_rng(20261122)` for `n_boot = 4,000` bootstrap resamples and `np.random.default_rng(k)` (`k = 0..99`) for gene permutations, plus deterministic `seed = 42` for subsampling.
- **Dual-Team `C3` Tolerance:** `max_point_abs_diff <= 1e-4` on deterministic metrics, `max_ci_abs_diff <= 0.01`, and `100%` consensus on `H-C1`, `H-C2`, `H-C3` across all `K_conf` lines (`1e-10` target parity on `FCT` matched point values).
- **Unvarnished Disclosure:** Negative or mixed findings (`未确认` on `H-C1/C2` or failure on `H-C3`) must be reported accurately across Abstract, Main Tables, and Conclusions without post-hoc subgroup promotion. Cite `arXiv:2608.00152v3` (Shoeibi & Yousefi, 2026) and Round 9 prior-art boundaries (`CellularAI`, `AtlasShift`, `CGC`, `AMM-SimWMag`, `CASTER`, `CIPHER v2`, `Molina & Zhang`) transparently.
