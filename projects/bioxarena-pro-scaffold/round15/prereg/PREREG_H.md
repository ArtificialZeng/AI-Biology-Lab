# Round 15 Pre-Release Prospective Multi-Line Validation Addendum (`PREREG_H.md`, `v2`)

> **v2 change log (Master Controller Claude, 2026-11-23 16:54 +08:00, per `round15/status/PREREG_H_AUDIT_R1.md`; no change to claims, tests, thresholds, endpoint, masks, sources, multiplicity or seed):** H1 §10.6 Fan et al. 2026 now read and verified (partial overlap with Claim 1); H2 new §3.1 line-identity rule for multi-line files; H3 new §8.4 release-wave rule; H4 §1.2 ground-truth wording limited to what is publicly verifiable; H5 §8.2 single-laboratory/platform scope phrase; H6 §11.2 cap addendum may change bytes only.

- **Project:** AI-Biology-Lab Round 15 — Pre-Release Preregistration for Arc Virtual Cell Challenge 2026 Open-License Public Release (`Cal-33`, `HeLa`, `AC16`, `MOLM13`, `U-87MG`; `Jurkat` strictly excluded)
- **Protocol Author:** Agent `A-T3` (Team A · Pre-Registration & Statistical Audit)
- **Parent Protocol (Locked, Unmodified):** [`PREREG_E.md` (`v2`)](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/prereg/PREREG_E.md) (`commit 4c606f4732ea080768e27757c689d40d7892cd3d`, `sha256 d10f92ced16a4519c43ea2b92b8b7ee881ab3d1764e39ffd10580473e4108924`; with Round 13 addendum [`PREREG_F.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round13/prereg/PREREG_F.md) `5b51e5b` and Round 14 addendum [`PREREG_G.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round14/prereg/PREREG_G.md) `d34560d`).
- **Authorizing Directive & Pre-Lock Inputs:** Master Controller Directive [`R15_第十五轮阶段2_PREREG_H指令.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/controller/directives/R15_第十五轮阶段2_PREREG_H指令.md), Feasibility Assessment [`R15_FEASIBILITY_VCC2026.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round15/shared/R15_FEASIBILITY_VCC2026.md), Prior-Art Supplement [`PRIOR_ART_H.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round15/shared/PRIOR_ART_H.md), [`A-T1.txt`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round15/status/A-T1.txt), and [`A-T5.txt`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round15/status/A-T5.txt).

---

## §0. Verbatim Inheritance Statement & New Seed (`沿用声明`)

`PREREG_H.md` is strictly an **addendum** to [`PREREG_E.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/prereg/PREREG_E.md). **Not a single word (`一字不改`)** is modified in the two confirmatory claims, sub-hypotheses, predictor equations, primary endpoint, gene masks, primary source pool, statistical tests, count/reproducibility thresholds, or multiplicity rules of `PREREG_E.md` (§0.1, §0.3, §2.1, §3.3, §4, §5, §7, §8), EXCEPT assigning the new Round 15 RNG seed **`20261231`** and specifying the pre-release governance, access gate, sunset clause, candidate lines, dual-schema ingestion, two-tier split-half, multi-line adjudication, secondary analyses, prior disclosures, and download rules below:

1. **Two Confirmatory Claims (Identical to `PREREG_E §4.2–§4.3`):**
   - **Claim 1 (`H-E1 = H-E1a ∩ H-E1b ∩ H-E1c`, Paired Target-Line Baseline Specificity):** Evaluated under `all_nontarget` ($G_{\mathrm{eval}}(t) \setminus \{g_p\}$) on every qualifying line $t \in \mathcal{K}_{\mathrm{conf}}$:
     - **`H-E1a`:** `FCT` achieves `r_orth_cal > 0` (**`97.5% CI_lo > 0`**) under `all_nontarget`, with `boot valid_frac >= 0.95` and `rho_ab_orth > 0`;
     - **`H-E1b`:** Across **each** of the `6` `r9_6` source lines $s'$ jointly (IUT over `6` pairs), the joint perturbation bootstrap paired contrast $\Delta^{\mathrm{pair}}(s') = r^{\mathrm{orth\_cal}}(\mathrm{FCT}) - r^{\mathrm{orth\_cal}}(\mathrm{NegCtrl\_FCT\_MismatchedMean}(s'))$ sharing identical perturbation draws and denominator $c^\perp_{\mathrm{ceil}}$ satisfies **`97.5% CI_lo > 0`**;
     - **`H-E1c`:** `FCT` `r_orth_cal` **`97.5% CI_lo > Q_{0.975}(\mathrm{NegCtrl\_FCT\_Perm}(0..99))`** under `all_nontarget`.
   - **Claim 2 (`H-E2 = H-E2a ∩ H-E2b`, Empirical Split-Half Top-50 Constant-Baseline Selection Inflation `SA-8`):** Tested independently of Claim 1 on every qualifying line $t \in \mathcal{K}_{\mathrm{conf}}$:
     - **`H-E2a`:** Under `resp_top50_split`, **`E1_ConstMed`** ($\hat{\Delta}^{\mathrm{E1}}_{t,p,g} = \operatorname{log1p}(c_{\mathrm{med}}\overline{\mathrm{FC}}_{p,g}) - \operatorname{log1p}(c_{\mathrm{med}})$, $c_{\mathrm{med}} = \operatorname{median}_{g \in G_{\mathrm{eval}}(t)} \max(\mu_{t,g}, 0)$) achieves `r_orth_cal_split` with **`97.5% CI_lo > 0`**;
     - **`H-E2b`:** Across the common perturbation set `split_valid ∩ all_valid`, the joint cross-mask perturbation bootstrap difference $\Delta^{\mathrm{cross}}_{\mathrm{SA8}} = r^{\mathrm{orth\_cal}}_{\mathrm{split}}(\mathrm{E1\_ConstMed}) - r^{\mathrm{orth\_cal}}_{\mathrm{all}}(\mathrm{E1\_ConstMed})$ (recomputing each mask's noise ceiling within identical draws) satisfies **`97.5% CI_lo > 0`**.
2. **Primary Endpoint, Masks, & Sources (Identical to `PREREG_E §2.1, §4.1`):**
   - Primary endpoint: **`r_orth_cal`** (`Transfer_Mean` reference, `orthogonalize_to_ref`, Spearman–Brown disattenuation clipped to `[-1, 1]`).
   - Masks: **`all_nontarget`** ($G_{\mathrm{eval}}(t) \setminus \{g_p\}$) and **`resp_top50_split`** (symmetric `A->B` / `B->A` out-of-sample top-50 responding mask); gene universe $G_{\mathrm{eval}}(t) = G_{\mathrm{R9\_6148}} \cap G_t^{\mathrm{fin}}$.
   - Primary training source cohort (**`r9_6`**): strictly the `6` harmonized Round 9 human `CRISPRi` cell lines (`nadig_hepg2`, `nadig_jurkat`, `vcc2025_h1`, `replogle_k562_essential`, `replogle_rpe1`, `nourreddine_kolf21j`; `6,148` `R9_COMMON` genes; `2,117` line-reproducible `u6_rep` symbols). Sensitivity pool **`SA-Src8`** (`--source-mode r9_r10_8`, adding `dld1` & `cd4t_rest`) remains descriptive only in `§9`.
3. **Multiplicity, Resampling (`seed = 20261231`), & `P5` Standardization:**
   - **Intersection–Union Test (IUT) + 2-Claim Bonferroni**: $\alpha_{\mathrm{claim}} = 0.025$ per claim, two-sided **`97.5%` percentile bootstrap confidence intervals (`[1.25%, 98.75%]`)**.
   - **`n_boot = 4,000`** joint perturbation bootstraps using **new Round 15 RNG seed `20261231`** (`np.random.default_rng(20261231)`), and **`n_perm = 100`** (`k = 0..99`, `np.random.default_rng(k)`) gene-label permutations on $\mu_t$.
   - Standardization (`10^4` libnorm -> `log1p` -> matched-unit NTC subtraction), `CanonicalHGNC` mapping, and `P5` eligibility (`n_cells >= 25`, `>= 50` top-5,000 finite genes with Anderson–Darling `BH FDR < 0.05` vs `1,000` stratified NT, `20-PC` energy-distance `p < 0.001` over `2,000` perms vs `500` NT, reproducible in `>= 1` `r9_6` source, self-gene in $G_{\mathrm{eval}}(t)$, `half_valid` `>= 2` cells/half) remain **identical word-for-word to `PREREG_E §3.3`**, alongside all `PREREG_E §5` descriptive analyses (`H-D2`, `FCT - E1_ConstMed`, `N_5`, `H-D_util`, `SA-1`, `SA-4/R19`, `SA-5`, `SA-6/S7`).
   - Clean-room dual-implementation tolerances (identical to `PREREG_E §7`): `max_point_abs_diff <= 1e-10`, `max_ci_abs_diff <= 1e-10` under aligned seed, and `100%` verdict parity.

---

## §1. Pre-Release Timing & Zero-Data-Access Blinding Guarantee (`时间和盲态`)

1. **Public GCS Bucket & Platform State Prior to Lock Commit (`2026-11-23`):**
   - Verified directly by Master Controller (`16:30 +08:00`) and independently by `A-T1` (`2026-11-23T08:38:50Z` / `16:38:50 +08:00`, final multi-platform re-check **`2026-11-23T08:41:41Z`** / `16:41:41 +08:00`, logged in [`R15_FEASIBILITY_VCC2026.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round15/shared/R15_FEASIBILITY_VCC2026.md)):
     - **`storage.googleapis.com/storage/v1/b/arc-institute-virtual-cell-atlas/o?delimiter=/&prefix=virtual-cell-challenge%2F2026%2F`**: returned **`prefixes = None, items = []` (`0` objects / `0` subdirectories)**.
     - **`virtual-cell-challenge/`**: contains ONLY `virtual-cell-challenge/2025/` (`adata_Training.h5ad`, `adata_Validation.h5ad`, `adata_Test.h5ad`, `gene_names.csv`, `pert_counts_*.csv`, uploaded `2025-12-16T17:08–17:36Z`, plus `FASTQ/`).
     - HuggingFace `arcinstitute` (`22` datasets, latest revision `2026-10-28T21:52:33Z`), NCBI GEO/BioProject, Zenodo, and Figshare hold **`0`** Arc VCC 2026 expression datasets or target tables. `arcinstitute.org/news` stops at `2026-11-19` (*Virtual Cell Challenge 2026: Winners and Wrap-Up*) with no open release date published.
2. **Strict Zero-Account & Zero-File Guarantee:**
   - This project holds **zero `virtualcellchallenge.org` or `vcc-cli` accounts**, has **never clicked or accepted any competition Terms or Persona selfie verification**, and has **downloaded zero 2026 files — including neither `pert_counts.csv` nor any `context_*.h5ad` control packages**.
   - This project has no access to any VCC 2026 perturbation data or ground truth. Per Arc's public pages as of the lock date, the 2026 ground truth has not been publicly released (Arc staff necessarily hold it; we make no claim about anyone else's access). Pushing this locked commit to the public repository establishes a verifiable timestamp for the analysis plan prior to open dataset release.

---

## §2. Open-License Access Gate & Sunset Clause (`取数条件与失效日期`)

1. **Strict Open-Release Trigger (`取数条件`):**
   Download and ingestion are permitted **IF AND ONLY IF ALL** of the following hold:
   - Arc Institute officially publishes the 2026 perturbation expression dataset in the open public GCS bucket (`gs://arc-institute-virtual-cell-atlas/virtual-cell-challenge/2026/`) or the Virtual Cell Atlas (`https://arcinstitute.org/tools/virtualcellatlas`) **without login or identity verification**; **AND**
   - An official Arc webpage explicitly states an open license — **`CC0`** (`CC0 1.0`) OR **`CC BY`** (`CC BY 4.0`) — authorizing unrestricted scientific download, derivative computation, and public dissemination.
2. **Mandatory Exclusion of Gated Channels (`EXCLUDED_ACCESS`):**
   - Any release channel requiring account registration, login credentials, Persona/anti-bot/human verification, or acceptance of restrictive Terms (including `https://virtualcellchallenge.org/terms`, effective `2026-08-04`, which restricts use solely to Challenge participation and prohibits derivative works or publication — even via human account holders) is **strictly classified `EXCLUDED_ACCESS`**: perform **ZERO downloads and ZERO analyses**.
3. **Hard Expiration / Sunset Clause (`失效日期: 2027-06-30`):**
   - If the VCC 2026 datasets have **NOT** been published openly under an explicit `CC0` or `CC BY` license on or before **`2027-06-30`**, **`PREREG_H.md` automatically expires (`自动失效`)**.
   - Upon expiration, ZERO post-hoc extensions or gated data acquisition may be attempted, and the manuscript MUST record truthfully: **「已于真值公开前完成预注册（附公开仓库锁定提交哈希），截至 2027-06-30 数据未以开放许可公开发布，未作确认性检验」**.

---

## §3. Candidate Cell Lines & Hard `Jurkat` Exclusion (`细胞系`)

| Line Key | Biological Context | VCC 2026 Role | Pre-Locked `PREREG_H` Status & Rule |
|---|---|---|---|
| **`cal33`** | Head and neck squamous carcinoma (`Cal-33`) | Validation (`Context A/B/C`) | **Candidate (`1/5`)** — independent evaluation |
| **`hela`** | Cervical adenocarcinoma (`HeLa`) | Validation (`Context A/B/C`) | **Candidate (`2/5`)** — independent evaluation |
| **`ac16`** | Cardiomyocyte (`AC16`) | Test (`Context D/E/F`) | **Candidate (`3/5`)** — independent evaluation |
| **`molm13`** | Acute myeloid leukemia (`MOLM13`) | Test (`Context D/E/F`) | **Candidate (`4/5`)** — independent evaluation |
| **`u87mg`** | Glioblastoma (`U-87MG`) | Test (`Context D/E/F`) | **Candidate (`5/5`)** — independent evaluation |
| **`jurkat`** | T-cell leukemia (`Jurkat`) | Validation (`Context A/B/C`) | **STRICTLY EXCLUDED** — homologous to `nadig_jurkat` inside `r9_6` training sources |

- **No Cross-Line Pooling:** Each surviving candidate line among `{cal33, hela, ac16, molm13, u87mg}` is ingested independently, passes or fails the `§7` count gate independently, and is evaluated and adjudicated per line (`max K_conf = 5`). Data MUST NEVER be pooled across cell lines.

### §3.1 Line-Identity Rule for Files Containing Several Lines (v2)
- Released files may hold several cell lines (the 2025 release was split by train/validation/test, not by line). Each cell's line is taken, in this order, from: (a) an explicit `.obs` column naming the cell line (e.g. `cell_line`, `cell_type`, `line`), matched case-insensitively after removing `-`, `_` and spaces to `cal33`, `hela`, `jurkat`, `ac16`, `molm13`, `u87mg`; else (b) a context column (e.g. `context`) mapped to lines only through a context-to-line table published by Arc with the release or in its official documentation (URL and quote recorded in `R15_QC_COUNTS.md` before any metric).
- If neither yields an unambiguous line for every cell, the affected lines are `EXCLUDED_FEASIBILITY`. Line labels are never inferred from expression.
- All `jurkat` cells are dropped at ingestion, before any QC count, centering or metric.

---

## §4. Primary vs Full-Screen Data Form Priority (`数据形式，事先写死`)

1. **Primary Analysis Rule (`First Open Release` / `Challenge Precedence`):**
   - For each candidate line, the **confirmatory primary analysis (`主分析`) MUST use the earliest open-licensed data form publicly released for that line**.
   - **Pre-Registered Co-Release Rule:** If both (a) the curated `300`-target competition dataset and (b) the complete unscreened (`>10,000`-target) screen are published simultaneously under an open `CC0`/`CC BY` license:
     - **Confirmatory Primary (`主分析`):** Strictly uses the **`300`-target challenge dataset** (smaller volume matching the `2025/` release format and directly commensurate with `vcc2025_h1`).
     - **Full-Screen Secondary (`§9`):** The **`>10,000`-target full screen** is evaluated strictly as a pre-registered descriptive secondary analysis (`只描述、不判定`).
2. **Pre-Locked Accounting of Arc Competition Curation:** All evaluations on challenge-curated subsets MUST explicitly disclose upstream Arc curation (removing targets failing median knockdown `>= 80%`, selecting the better-performing of 2 sgRNAs per gene, and uniform downsampling to `400` cells/perturbation).

---

## §5. Integer Raw Counts, NTC, Single-Construct, Gene Symbols, & Cell QC (`细胞和计数`)

1. **Mechanical Pre-Metric Integer Raw-UMI Verification (`EXCLUDED_FEASIBILITY`):**
   - Immediately upon authorized download (and **before computing any model prediction or evaluation metric**), mechanically inspect `/raw/X`, `/X`, and all non-excluded `/layers/*` (excluding subfraction or normalized layers such as `nascent`, `unspliced`, `spliced`, `ambiguous`, `velocity`, `intron`, `exon`, `lognorm`, `log1p`, `scale`, `norm`, `cpm`, `tpm`, `vst`, `sct`):
     - Verify across **100% of stored non-zero entries** (chunked streaming, zero sampling) that all values are finite, non-negative **exact integers** (`>= 0`, `x == floor(x)`).
     - If `>= 2` non-excluded slots pass integer verification (e.g., `/X` and `/raw/X`), verify exact shape, `/var` gene axis alignment, and bit-exact matrix identity across all rows; if ANY non-excluded integer slots disagree in shape, gene axis, or matrix values -> **STOP and mark `EXCLUDED_FEASIBILITY`** (never pick one arbitrarily).
     - If ZERO integer raw-UMI slots exist -> **`EXCLUDED_FEASIBILITY`**, record exact diagnostic reason in `round15/shared/R15_QC_COUNTS.md`, and skip all downstream evaluation on that line.
2. **Schema-Agnostic `.obs` Ingestion (Supports Both 2025 & 2026 Published Schemas):**
   - Reader supports both `2025`-style `.obs` (`target_gene`, `guide_id`, `batch`) and `2026`-controls-style `.obs` (`target_gene`, `context`, `ntc_id`, optional `batch`/`channel`) and resolves target symbols strictly through frozen `hgnc_complete_set.txt` (`HGNCLookup.resolve`, `pert_key = canon_sym`), stripping any non-human or synthetic Cas9 vector features prior to mRNA depth QC.
3. **Locked NTC Rule:**
   - **NTC = strictly `target_gene == "non-targeting"`** (matching normalized `"non-targeting"` strings while excluding `"unassigned"`, `"multiplet"`, `"doublet"`, `"safe-targeting"`, `"gene_desert"`, or empty/NaN labels).
4. **Locked Single-Construct / Single-sgRNA Rule:**
   - Challenge data assigns a single construct per gene (`target_gene`, optionally accompanied by `guide_id` or `ntc_id`).
   - Any cell with heterogeneous multi-target or multi-construct calls (`target_gene` containing multi-gene delimiters `[,;|+\s]`, explicit `n_guides > 1`/`nperts > 1`, `multiplet`, `doublet`, or `unassigned`) is strictly excluded; single-construct dual-promoter cassette IDs in `guide_id` (e.g., `2025`-style `non-targeting_A|non-targeting_B` or `GENE_1|GENE_2` representing a single targeted construct) are permitted **iff** `target_gene` unambiguously designates a single gene or `"non-targeting"` (and any remaining multi-construct `guide_id` with conflicting targets is excluded).
   - If a full-screen or un-called release lacks per-cell construct assignments and exposes a raw guide count matrix, apply **`PREREG_F §3.2`** Poisson-Gaussian histogram thresholding across cells with `>= 200` genes (`fit_poisson_gaussian_thresholds_from_hist`, `min_thr=5, max_thr=20, default_thr=7`, retaining `n_pos == 1`, `top1 >= 5`, `top2 < 3 or top1 / max(top2, 0.5) >= 3.0`).
5. **Cell QC (`max_pct_mt = None`):**
   - Retain valid single-construct cells satisfying **`n_umis >= 300`** AND **`n_genes >= 200`**, with **NO mitochondrial upper limit (`max_pct_mt = None`)**.

---

## §6. Two-Tier Disjoint Split-Half & NTC Subtraction (`分半和扣除`)

Deterministically select split tier **per candidate line** strictly prior to metric evaluation and record exact metadata into **`round15/shared/R15_QC_COUNTS.md`**:

1. **Tier 1 — Hardware Channel / Batch Parity (Preferred when available):**
   - Inspect `.obs` for channel/batch columns (`batch`, `flex_channel`, `channel`, `lane`, `replicate`, `sample`, ordered priority; NOT global line/context label when identical across all cells of a line).
   - Trigger Tier 1 **iff** a channel/batch column exists with **`>= 2` distinct units** yielding `>= 15` QC-valid NTC cells per half:
     - Sort unit labels in standard Python **lexicographical order** (`1`-based position `i = 1, 2, ...`);
     - **Odd positions (`i % 2 == 1`) -> `Half A` (`half_1`)**; **Even positions (`i % 2 == 0`) -> `Half B` (`half_2`)**.
     - Subtract **unit-matched** NTC `log1p` mean per channel/batch unit (for any individual unit holding `< 15` NTCs while its assigned half holds `>= 15`, fall back cleanly to within-half pooled NTC mean/SD via `_apply_min_ntc_batch_fallback` — zero cross-half control sharing).
2. **Tier 2 — Stratified Barcode-Sorted Seed-`20261231` Fallback (When Channel/Batch Metadata Absent or <2 Units):**
   - Triggered automatically whenever `.obs` lacks `>= 2` channel/batch units (e.g., if open 2026 `.h5ad` omits `batch` and includes only `target_gene, context, ntc_id`) or fails Tier 1 `>= 15` NTC/half under parity:
     - Partition cells into strata by **`canon_sym`** (for targeting perturbations) and **`"__NTC__"`** (for controls), iterating strata keys in **Python ascending string order**;
     - Within each stratum, sort cell indices ascending by **`(obs_barcode, raw_row_idx)`**;
     - Using one generator **`rng = np.random.default_rng(20261231)`**, permute each stratum (`perm = rng.permutation(n)`): **`perm[:n // 2] -> Half A (1)`**, **`perm[n // 2:] -> Half B (2)`**;
     - Center strictly per half (`batches = ["half_1", "half_2"]`, zero cross-half control sharing).

---

## §7. Per-Line Confirmatory Count Threshold (`门槛`, Identical to `PREREG_E §2.2`)

Candidate line $t \in \{\text{cal33, hela, ac16, molm13, u87mg}\}$ enters confirmatory adjudication (`CONFIRMATORY_LINES_R15.txt`) **if and only if ALL four conditions hold simultaneously**:
1. Technical open-access & integer raw-count feasibility (`§2`, `§5`) passes (`feasible == PASS`);
2. **`n_eval_targets >= 50`** (`P5`-reproducible on line $t$, reproducible in `>= 1` `r9_6` source, self-gene in $G_{\mathrm{eval}}(t)$, and `half_valid`);
3. **`|G_eval(t)| >= 3,000`** (`|G_R9_6148 ∩ G_fin(t)| >= 3,000`);
4. **`per-half NTC >= 15`** (`>= 15` QC-valid `"non-targeting"` singlets in `Half A` AND `>= 15` in `Half B`).

- Any line failing `n_eval_targets >= 50`, `|G_eval(t)| >= 3,000`, or `per-half NTC >= 15` while feasible becomes **`SECONDARY_ONLY`** (descriptive reporting only, excluded from confirmatory `K_conf`). Any infeasible line becomes **`EXCLUDED_FEASIBILITY`** (zero downstream metrics).

---

## §8. Multi-Line IUT + Bonferroni Adjudication (`多线判定`)

Let **`K_conf ∈ {0, 1, 2, 3, 4, 5}`** denote the exact count of candidate lines passing `§7` (`CONFIRMATORY_LINES_R15.txt`):

1. **Line-Level & Cross-Line IUT Rule:**
   - Every qualifying line $t \in \mathcal{K}_{\mathrm{conf}}$ is evaluated and reported individually (`est [97.5% CI]`, `valid_frac`, `rho_ab_orth`, `6` `H-E1b` paired contrasts, `Perm100`, `H-E2a/b`).
   - Within each claim (`H-E1` or `H-E2`), **confirmation requires unanimous `K_conf / K_conf` (`100%`) passage across ALL `K_conf` qualifying lines** (Intersection–Union Test across confirmatory lines; 2-Claim Bonferroni `[1.25%, 98.75%]` CI unchanged).
2. **Pre-Locked Scope Phrases by `K_conf` (Non-Negotiable):**
   - **`K_conf == 0`:** Execute **ZERO confirmatory adjudication (`CONFIRMATORY_ABORTED_K0`)**; emit descriptive `SECONDARY_ONLY` metrics only (if feasible) alongside prior non-confirmatory disclosures.
   - **`K_conf == 1`:** Confirmatory adjudication executes on that single line, AND every summary/abstract **MUST explicitly disclose `K_conf = 1` single-line confirmation with limited cross-line generalizability (`单线，推广性有限`)** — NEVER claim multi-line validation.
   - **`K_conf >= 2` (and unanimous `K_conf / K_conf` passage):** Claims passing `K_conf / K_conf` MAY be stated as **"confirmed across `K_conf` unseen human cell lines from a single laboratory and platform (Arc Institute VCC 2026, 10x Flex CRISPRi)"**. The single-laboratory/platform qualifier is mandatory, and lines sharing a challenge panel share the same target set, so they are not independent target samples.
3. **Strict Anti-Flexibility Lock:** Post-hoc lowering of thresholds (`50` / `3,000` / `15`), post-hoc switching of primary data regime (`§4`), or selective omission of a failing qualifying line from `K_conf` is strictly prohibited.
4. **Release-Wave Rule (v2):** If lines are released openly at different times (for example validation lines before test lines), each release wave is processed as soon as it meets `§2`, with its own counts-only gate, freeze check and per-wave verdict, reported immediately and in full. The **overall** statement for each claim covers every line that qualified in any wave up to the end of the study: a claim is confirmed overall only if it passed on all qualifying lines of all waves (IUT across waves), so a later wave can only keep or remove an earlier confirmation, never add one that an earlier wave failed. The study ends at the first of (a) all five candidate lines adjudicated or excluded, or (b) the `§2.3` sunset date; the overall statement is written then, and any per-wave statement made earlier is labelled "interim, more lines pending".

---

## §9. Pre-Registered Secondary Analyses (`事先登记的次要分析`, Strictly Descriptive / `只描述、不判定`)

All items below are computed identically under `seed = 20261231` and reported in full alongside confirmatory results **without altering Claim 1 or Claim 2 adjudication**:

1. **`SA-H1` (Leave-`vcc2025_h1`-Out Source Sensitivity):** Exclude `vcc2025_h1` (same Arc Institute + 10x Flex + UG100 + `cyto` + CRISPRi platform as VCC 2026) from `r9_6`, retaining strictly the `5` remaining `r9_6` sources (`nadig_hepg2`, `nadig_jurkat`, `replogle_k562_essential`, `replogle_rpe1`, `nourreddine_kolf21j`); recompute `FCT r_orth_cal`, the `5` paired `H-E1b` contrasts `FCT - Mis_s`, `Perm100`, and `H-E2` (`SA-8`) on every evaluated line to quantify platform-matching inflation vs true cross-platform biology.
2. **`SA-Panel` (`129` Shared Val–Test Targets vs Line-Specific Targets):** Partition `eval_targets` of each line into (a) targets belonging to the `129` shared perturbations between the VCC 2026 validation (`Cal-33`, `HeLa`) and test (`AC16`, `MOLM13`, `U-87MG`) panels vs (b) the remaining non-shared panel targets (excluding the `129` overlap), reporting subgroup counts, `FCT r_orth_cal [97.5% CI]`, and `H-E1b` paired contrasts per subset (`只描述、不判定`).
3. **Full Unscreened Screen (`>10,000` targets):** Evaluated descriptively under identical metrics whenever co-released with a primary `300`-target challenge dataset (`§4`).
4. **`PREREG_E §5` Full Descriptive Suite:** `N_5` semi-partial analysis (including explicit non-evidentiary algebraic-zero tags for `Transfer_Mean_GlobalScale` / `Source_Grand_Mean` and cone-suppression tags for `E1/E2/E6/E7`), `SA-Src8` (`--source-mode r9_r10_8`), legacy `H-D2`, `FCT - E1_ConstMed`, `H-D_util`, `SA-1` slopes, `SA-4/R19`, `SA-5`, and `SA-6/S7` (`N_cells` & `|KD|` confounders).

---

## §10. Mandatory Prior Disclosures (`事先披露`)

Before observing any VCC 2026 counts, targets, or matrices:
1. **Post-R12–R14 Design Provenance:** Round 15 was preregistered after observing Round 12 (`a549_ko` confirmed; `hek293` count-gate failed), Round 13 (`mac_ctrl` `3`-target gate fail; `mac_lps` excluded), and Round 14 (`a549_myllia` `32`-target gate fail), in addition to the post-hoc R9–R11 provenance of `H-E1` and `H-E2`.
2. **Unpublished Target List & Scenario Yield Projections ([`R15_FEASIBILITY_VCC2026.md §3`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round15/shared/R15_FEASIBILITY_VCC2026.md#L58-L71)):**
   - The `300`-target identities (`pert_counts.csv`) and `>10,000` full target list are **NOT publicly accessible without gated login**; only public aggregate overlaps (`277/300` test targets in Replogle 2022 K562 genome-wide screen, `298/300` in X-Atlas/Orion) are known. **`u6_rep ∩ R9_COMMON` projections remain unreliable until open release.**
   - Using `vcc2025_h1` leave-one-out (`300` targets vs `u5_rep = 2,004`, `u5_rep ∩ R9_COMMON = 1,610`) as reference, `95/300` overlapped `u5_rep` (`84` in `u5 ∩ R9_COMMON`) with `91/95 = 95.8%` (`80/84 = 95.2%`) passing Addendum-B. For VCC 2026 `300`-target panels (assuming `~95–150` in `u6_rep` / `~80–125` in `u6_rep ∩ R9_COMMON`): a `20–40%` rate yields `19–60` (`16–50` in `R9_COMMON`, borderline `< 50`); conservative `HEK293` `9.76%` (`16/164`) yields `9.3–14.6` (`< 50`); optimistic `A549` `83.96%` (`445/530`) yields `79.8–125.9` (`67–105` in `R9_COMMON`, likely `>= 50`). If full `>10,000` screens open, yield could reach `117–156+` (`9.76%`) to `240–640+` (`20–40%`).
3. **Arc Challenge Curation Bias:** Upstream `>= 80%` median knockdown selection, choosing the better of `2` sgRNAs per gene, and `400`-cell downsampling truncate weak/moderate perturbation tails and may inflate effect size/reproducibility relative to genome-wide screens.
4. **Shared Platform Confounding with `vcc2025_h1`:** Shared Arc Institute + 10x Flex (`18,533` probes) + Ultima UG100 + `cyto` + CRISPRi pipeline requires mandatory `SA-H1` (`leave-vcc2025_h1-out`) disclosure.
5. **Cross-Panel Target Overlap (`129`) & `Jurkat` Exclusion:** Validation and test `300`-target panels share **`129`** targets (`43%`; also `21` targets in `U-87MG` reportedly overlap McFaline-Figueroa 2024); **`Jurkat` is hard-excluded** due to `nadig_jurkat` in `r9_6`.
6. **Prior-Art & Unread-Paper Boundary ([`PRIOR_ART_H.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round15/shared/PRIOR_ART_H.md)):**
   - Cite **Kawsar et al. (2026)** (`CanPert/CASTER`, `doi:10.64898/2026.11.15.766815`) for per-target paired wrong-line control comparison on `K562/RPE1/HepG2/Jurkat` (partially overlapping Claim 1), and note VCC 2026 competitor denoising-PC checks (`jam-sudo`, `1minute99`).
   - **Fan et al. 2026 (v2, read and verified):** Xingyu Fan, Jinghao Wang, Chang-Yu Hsieh, Chunbin Gu, Pheng-Ann Heng, *Zero-Shot Cross-Cell-Line Perturbation Prediction Needs Controls*, NeurIPS 2026 Evaluations and Datasets Track (poster; OpenReview `lnpzb1HvwX`; CC BY 4.0; code `github.com/xyfan22/crossline-controls`). The controller read the abstract, meta-review, all reviews and author responses on 2026-11-23 (PDF full text not read); notes in `round15/status/PRIOR_ART_H_CONTROLLER_FAN2026.md`. **Partial overlap with Claim 1:** they already show that target-line controls drive cross-line scores and introduce target-control-copy and shifted-control (target controls plus a source-derived mean shift, the additive analogue of `Transfer_Mean` / `FCT`) as required baselines, on Replogle/Nadig essential genes with Jurkat held out. Not done there (per abstract and reviews): replacing the target line's controls with mismatched lines' controls in a paired test, amplitude orthogonalization (`r_orth_cal`), noise-ceiling calibration, or a preregistered prospective test on unseen lines. Claim 1 MUST be described as building on Fan et al. and Kawsar et al., never as first.
   - **Unread Dobin Preprint (`doi:10.64898/2026.11.17.767214`):** Returns HTTP 404/403 across APIs (**`未能读取`**); readable `cell-eval2` docs & Wrap-Up (`2026-11-19`) do not define the `0.02–0.13 vs 1.02` line-specific slope formula and do not evaluate `FCT`-class unanchored predictors.

---

## §11. Locked Download & Storage Plan (`下载方案`)

Zero files may be downloaded until: (i) lock commit of `PREREG_H.md`, (ii) verified open `CC0`/`CC BY` public release per `§2`, and (iii) Master Controller `R15_GO.txt` authorization.

1. **Authorized Download Channels:** Strictly official unauthenticated public GCS (`gs://arc-institute-virtual-cell-atlas/virtual-cell-challenge/2026/` / HTTPS JSON API) or Virtual Cell Atlas official open URLs (`https://arcinstitute.org/tools/virtualcellatlas`). Reject any `virtualcellchallenge.org/terms`-gated or credentialed URL (`EXCLUDED_ACCESS`).
2. **Volume Cap & Storage Path:**
   - Hard cumulative download cap: **`60.0 GB` (`60,000,000,000` bytes)** across primary challenge files (for reference, the entire `2025/` 3-file release totaled `34.36 GB`; if co-released `>10,000` full screens exceed `60 GB`, prioritize primary `300`-target files per `§4` and seek explicit controller cap addendum before touching larger full-screen archives). A cap addendum, if any, is written and committed before the files concerned are downloaded and may change only the byte cap and file list, never which data form is primary or any analysis rule.
   - Store files strictly under **`round15/shared/raw/vcc2026/`**.
   - Record exact byte counts and SHA-256 digests of every downloaded file in **`round15/shared/R15_CHECKSUMS.sha256`** and **`round15/shared/download_manifest_r15.json`**.
   - Inspect and ingest exclusively via **`python3 -I`**, running counts-only integer/QC checks (`round15/shared/R15_QC_COUNTS.md`) first before any evaluation or metric computation. Never execute `git`, never create accounts, and never write outward communications.
