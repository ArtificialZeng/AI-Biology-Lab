# Round 14 Independent-Dataset Same-Line Replication Pre-Registration Addendum (`PREREG_G.md`, `v1`)

- **Project:** AI-Biology-Lab Round 14 — Independent-Lab / Independent-Chemistry Replication Test in `A549` CRISPRko (`GSE337804`, Myllia Biotechnology / Adduri et al., *Cell* 2026)
- **Protocol Author:** Agent `A-T3` (Team A · Pre-Registration & Statistical Audit)
- **Parent Protocol (Locked, Unmodified):** [`PREREG_E.md` (`v2`)](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/prereg/PREREG_E.md) (`commit 4c606f4732ea080768e27757c689d40d7892cd3d`, `sha256 d10f92ced16a4519c43ea2b92b8b7ee881ab3d1764e39ffd10580473e4108924`; Round 12 freeze commit `11d91c72025264b838af091abcf0eac71588eb91`), alongside locked Round 13 addendum [`PREREG_F.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round13/prereg/PREREG_F.md) (`commit 5b51e5b079c2c01b44408b9ed9f8d6027b4ae4f5`).
- **Authorizing Directive & Pre-Lock Inputs:** Master Controller Directive [`R14_第十四轮阶段2_PREREG_G指令.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/controller/directives/R14_第十四轮阶段2_PREREG_G指令.md), Feasibility Assessment [`R14_FEASIBILITY_GSE337804.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round14/shared/R14_FEASIBILITY_GSE337804.md), [`A-T1.txt`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round14/status/A-T1.txt), and public `GSE337804_sgRNA_lib.csv.gz` (`12,261` bytes, `sha256 96f6eb62449ba276d36591af2f7e0fca47ffedb9f19d711bcc27671d50882df6`).

> **Strict Zero-Data-Access Blinding Statement:** As of submitting this addendum, **zero `.h5ad` expression matrices (`GSE337804_A549_ch_1.RDS.h5ad`, `563,808,587` B; `GSE337804_A549_ch_2.RDS.h5ad`, `608,829,330` B) have been downloaded or opened**, and `round14/shared/raw/` does not exist. Only public GEO/SRA/PubMed metadata and small non-expression files (`GSE337804_2cl_description.pdf`, `GSE337804_sgRNA_lib.csv.gz`) have been inspected. All Round 10–13 frozen/status files remain unmodified.

---

## §1. Verbatim Inheritance Statement & Round 14 Seed (`沿用声明`)

`PREREG_G.md` is strictly an **addendum** to [`PREREG_E.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/prereg/PREREG_E.md). **Not a single word (`一字不改`)** is altered in any claim, hypothesis, predictor formula, residual endpoint, gene mask, source cohort, inference rule, multiplicity correction, eligibility criterion, or descriptive sensitivity analysis from `PREREG_E.md` (§0.1, §0.3, §2.1, §3.3, §4, §5, §7, §8), EXCEPT setting the new Round 14 bootstrap RNG seed and specifying the `GSE337804` ingestion, interpretive, subgroup-descriptive, disclosure, and download rules below:

1. **Two Confirmatory Claims (Identical to `PREREG_E §4.2–§4.3`):**
   - **Claim 1 (`H-E1 = H-E1a ∩ H-E1b ∩ H-E1c`, Paired Target-Line Baseline Specificity):** Confirmed iff ALL three hold simultaneously on `a549_myllia` (assuming it clears §6 gate):
     - **`H-E1a`:** `FCT` achieves `r_orth_cal > 0` (`97.5% CI_lo > 0`) under `all_nontarget`, with `boot valid_frac >= 0.95` and `rho_ab_orth > 0`;
     - **`H-E1b`:** Across **each** of the `6` source lines `s' in r9_6` jointly (IUT over all `6` pairs), the joint perturbation bootstrap paired contrast $\Delta^{\mathrm{pair}}(s') = r^{\mathrm{orth\_cal}}(\mathrm{FCT}) - r^{\mathrm{orth\_cal}}(\mathrm{NegCtrl\_FCT\_MismatchedMean}(s'))$ under `all_nontarget` (sharing identical perturbation draws and denominator $c^\perp_{\mathrm{ceil}}$ within each draw) satisfies **`97.5% CI_lo > 0`**;
     - **`H-E1c`:** `FCT` `r_orth_cal` **`97.5% CI_lo > Q_{0.975}(\mathrm{NegCtrl\_FCT\_Perm}(0..99))`** under `all_nontarget`.
   - **Claim 2 (`H-E2 = H-E2a ∩ H-E2b`, Empirical Split-Half Top-50 Constant-Baseline Selection Inflation `SA-8`):** Confirmed independently of Claim 1 iff BOTH hold simultaneously on `a549_myllia` (assuming it clears §6 gate):
     - **`H-E2a`:** Under `resp_top50_split`, **`E1_ConstMed`** achieves `r_orth_cal_split` with **`97.5% CI_lo > 0`**;
     - **`H-E2b`:** Across the common perturbation set `split_valid ∩ all_valid`, the joint cross-mask perturbation bootstrap difference $\Delta^{\mathrm{cross}}_{\mathrm{SA8}} = r^{\mathrm{orth\_cal}}_{\mathrm{split}}(\mathrm{E1\_ConstMed}) - r^{\mathrm{orth\_cal}}_{\mathrm{all}}(\mathrm{E1\_ConstMed})$ (recomputing each mask's noise ceiling inside identical draws) satisfies **`97.5% CI_lo > 0`**.
2. **Predictors, Primary Endpoint, Masks, & Sources (Identical to `PREREG_E §2.1, §4.1`):**
   - Primary endpoint: **`r_orth_cal`** (`Transfer_Mean` reference, `orthogonalize_to_ref`, Spearman–Brown split-half ceiling disattenuation clipped to `[-1, 1]`). **`E1_ConstMed`** ($\hat{\Delta}^{\mathrm{E1}}_{t,p,g} = \operatorname{log1p}(c_{\mathrm{med}}\overline{\mathrm{FC}}_{p,g}) - \operatorname{log1p}(c_{\mathrm{med}})$ with scalar $c_{\mathrm{med}} = \operatorname{median}_{g \in G_{\mathrm{eval}}(t)} \max(\mu_{t,g}, 0)$) is strictly a **target-baseline-free transfer predictor** (`不使用目标线基线的迁移预测器` — NOT independent of source-line fold changes $\overline{\mathrm{FC}}_{p,g}$).
   - Masks: **`all_nontarget`** ($G_{\mathrm{eval}}(t) \setminus \{g_p\}$) and **`resp_top50_split`** (symmetric `A->B` / `B->A` out-of-sample top-50 responding mask); gene universe $G_{\mathrm{eval}}(t) = G_{\mathrm{R9\_6148}} \cap G_t^{\mathrm{fin}}$.
   - Training source cohort: Strictly the **`6`** Round 9 harmonized human **`CRISPRi`** cell lines (**`r9_6`**: `nadig_hepg2`, `nadig_jurkat`, `vcc2025_h1`, `replogle_k562_essential`, `replogle_rpe1`, `nourreddine_kolf21j`; `6,148` `R9_COMMON` genes; `2,117` line-reproducible `u6_rep` symbols). Sensitivity pool **`SA-Src8`** (`--source-mode r9_r10_8`, adding Round 10 `dld1` & `cd4t_rest`) remains strictly descriptive in `§5` (`只报告、不判定`). No `A549` dataset exists in `r9_6` or `r9_r10_8`.
3. **Inference, Resampling, Addendum-B `P5`, & New Seed:**
   - Multiplicity: **Intersection–Union Test (IUT) + 2-Claim Bonferroni** ($\alpha_{\mathrm{claim}} = 0.025$ per claim, two-sided **`97.5%` percentile bootstrap confidence intervals `[1.25%, 98.75%]`**).
   - Resampling: **`n_boot = 4,000`** joint perturbation bootstraps using **new Round 14 RNG seed `20261224`** (`np.random.default_rng(20261224)`), and **`n_perm = 100`** (`k = 0..99`, `np.random.default_rng(k)`) gene-label permutations on $\mu_t$.
   - Standardization & `P5` target eligibility: Identical to `PREREG_E §3.3` (`10^4` library-size normalization, channel-matched NTC `log1p` subtraction, Addendum-B `n_cells >= 25`, `>= 50` top-5,000 finite genes with Anderson–Darling `BH FDR < 0.05` vs `1,000` stratified NT, `20-PC` energy-distance `p < 0.001` over `2,000` perms vs `500` NT, reproducible in `>= 1` `r9_6` source, self-gene in $G_{\mathrm{eval}}(t)$, `half_valid` `>= 2` cells in Half A AND `>= 2` in Half B).
   - All `PREREG_E §5` descriptive analyses (`H-D2`, `FCT - E1_ConstMed`, `N_5`, `SA-Src8`, in-sample `resp_top50`, `H-D_util`, `SA-1`, `SA-4/R19`, `SA-5`, `SA-6/S7`), clean-room cross-check tolerances (`max_point_abs_diff <= 1e-10`, `max_ci_abs_diff <= 1e-10`, `100%` verdict parity), and literature/resource provenance boundaries:
     - **Kawsar et al. (2026)** (`CanPert/CASTER`, `doi:10.64898/2026.11.15.766815`) for prior per-target paired correct-vs-mismatched control comparison;
     - **Nicol, Shivakumar & Irizarry (2026)** (`doi:10.64898/2026.05.07.723486`) for finite-sample DEG selection inflation (noting Claim 2 `SA-8` is empirical on `r_orth_cal`, **unknown mechanism** — noiseless simulation predicted wrong sign on `8/12` lines under `all_nontarget` — and does NOT extrapolate to DEG metrics);
     - **Liu et al. (2026)** (`GSE345058`, `doi:10.64898/2026.06.01.728087`) strictly as Round 12 `a549_ko` data source;
     - **Adduri et al. (2026)** (*Cell* `2026 Sep 17`, Epub `2026 Aug 31`, PMID `42673963`, `doi:10.1016/j.cell.2026.07.052`, BioProject `PRJNA1491103`) strictly as the `GSE337804` data source (`Myllia Biotechnology × UCSF/Arc`). Confirmed **NOT** `X-Atlas/Orion`, **NOT** Arc Virtual Cell Challenge (`VCC 2026` Flex CRISPRi), and **NOT** Myllia's "Echoes of Silenced Genes" CRISPRi competition.

---

## §2. Line Definition (`GSE337804`, Independent-Lab `A549` `CRISPRko` / `CRISPRn`)

| Line Key | Cohort & Samples | Platform & Perturbation Design | Pre-Locked Role Rule |
|---|---|---|---|
| **`a549_myllia`** | **Basal unstimulated clonal Cas9 `A549` (`GSE337804`, Day 16)** across `2` independent channels: `ch_1` (`GSM9860781` WTS / `GSM9860782` sgRNA) & `ch_2` (`GSM9860783` WTS / `GSM9860784` sgRNA) | Low-MOI (`~6% EGFP+`, puro-selected) `SpCas9` **CRISPRko (`CRISPRn`)** (`218` genes × `3` guides + `20` `ctrl_NT` + `20` `ctrl_T`), profiled on **10x Next GEM `5' HT v2` + Jumpcode `DepleteX`** (GENCODE v39 GRCh38) | **SOLE confirmatory candidate line (`唯一的确认性候选线`, max `K_conf = 1`)** |

---

## §3. Pre-Locked Interpretive Scope (`同一细胞系的独立数据集复现`)

Written **strictly prior to downloading or viewing any `.h5ad` expression matrix**:

1. **Same-Line, Independent-Dataset Scope ONLY:** Round 14 tests replication on an **independent dataset of the SAME human cell line (`A549`) generated by an independent laboratory (`Myllia Biotechnology` vs Round 12 `Liu et al.`), using a different CRISPRko library (`218` genes × `3` distinct sgRNAs), different 10x chemistry (`5' HT v2 + DepleteX` vs `3' CROP-seq`), and mostly unseen target genes (`46/61 = 75.4%` not in R12 `a549_ko`)**. Passing Round 14 does **NOT** establish cross-cell-line generalization beyond `A549`.
2. **Permitted Claim Wordings (`No Cross-Round Pooling`):**
   - **If `a549_myllia` passes the §6 count gate AND both Claim 1 (`H-E1a/b/c`) and Claim 2 (`H-E2a/b`) pass:** Results MUST ONLY be stated as **"confirmed independently on two separate A549 datasets in two separately preregistered prospective tests"** (`「在两个独立的 A549 数据集上、两轮分别预注册的前瞻检验中各自确认」`: Round 12 `GSE345058` and Round 14 `GSE337804`). **Round 12 and Round 14 data MUST NEVER be pooled into a joint hypothesis test or joint confidence interval.**
   - **If `a549_myllia` passes the §6 count gate AND either Claim 1 or Claim 2 fails:** Report honestly and prominently that the unconfirmed claim **failed to replicate on the independent `A549` dataset (`在独立数据集上没有复现`)**, noting Round 12 `GSE345058` conclusions remain unchanged and examining the pre-registered §7 subgroup comparison without altering confirmatory verdicts.
   - **If `a549_myllia` fails the §6 count gate (`K_conf = 0`):** Report honestly as a gate failure / descriptive-only result (`不做确认性判定`), leaving Round 12 unchanged.

---

## §4. Integer Raw-Count Gate, Single-sgRNA, Locked NTC, Target Aggregation, & Cell QC (`细胞和 guide`)

1. **Mandatory Pre-Metric Integer Raw-UMI Feasibility Gate (`EXCLUDED_FEASIBILITY`):**
   - Because `GSE337804_A549_ch_1.RDS.h5ad` and `GSE337804_A549_ch_2.RDS.h5ad` are `.h5ad` exports from Seurat v3, immediately upon download and **strictly before computing any predictor or metric**, inspect `/X`, `/raw/X`, and `/layers/*` across both files mechanically (`A-T1` + independent `A-T3` verification):
     - Select the raw integer count matrix unambiguously via pre-locked priority: inspect `/raw/X`, `/X`, and non-subfraction `/layers/*` (excluding subfraction layers named `nascent`, `unspliced`, `spliced`, `ambiguous`, `velocity`, `intron`, `exon`, `lognorm`, `scale`, matching `PREREG_E §3.2` / `READER_XCHECK_R12.md`) and verify finite, non-negative, **exact integer raw UMI values** across all non-zero entries (preferring raw `/raw/X` or `/X` over derived layers when both are integer raw UMIs and non-negative integer counts match mRNA raw depth; record the chosen HDF5 dataset path and integer proof in `round14/shared/R14_QC_COUNTS.md` before any downstream evaluation).
     - **If NO matrix slot contains verified non-negative finite integer raw UMI counts**, OR if per-cell single-sgRNA & explicit `ctrl_NT` vs `ctrl_T` labels cannot be verified unambiguously across both channels -> **`a549_myllia` is immediately classified `EXCLUDED_FEASIBILITY`** with zero downstream metrics computed.
2. **Single-sgRNA Cell Calling:**
   - **Primary rule:** Assign single-sgRNA identity from the authors' per-cell guide columns in `.obs` (`sgrna_id` and `sgrna_symbol`, verified against `GSE337804_sgRNA_lib.csv.gz`), keeping **ONLY** cells assigned to **strictly `1` guide** (`!has_multi_delim`, not empty/`nan`/`unassigned`/`multiplet`, exact match to one valid `id` in `GSE337804_sgRNA_lib.csv.gz`).
   - **Pre-registered fallback (applied ONLY if per-cell guide calls are absent and raw per-guide count features exist):** Apply `PREREG_F §3.2` per-line Poisson-Gaussian histogram thresholds across cells with `>= 200` detected genes (`fit_poisson_gaussian_thresholds_from_hist` with `min_thr=5, max_thr=20, default_thr=7`, retaining `n_pos == 1`, `top1_UMI >= 5`, `(top2_UMI < 3 or top1_UMI / max(top2_UMI, 0.5) >= 3.0)`). Record branch deterministically in `R14_QC_COUNTS.md`.
3. **Locked NTC Rule (Strict Exclusion of `ctrl_T` Safe-Targeting Controls):**
   - Verified inspection of `GSE337804_sgRNA_lib.csv.gz` (`12,261` B, `694` rows across `id, symbol, rank, class, entrez, guide, count`) confirms `2` cleanly separated control classes:
     - **`ctrl_NT` (`20` guides):** `id == "non-targeting_1".."non-targeting_20"`, `symbol == "non-targeting"`;
     - **`ctrl_T` (`20` guides):** `id == "gene_desert_1".."gene_desert_20"`, `symbol == "gene_desert"` (genomic desert safe-targeting Cas9 cutting controls).
   - **Pre-Locked Rule:** **ONLY the `20` `ctrl_NT` guides (`non-targeting_1..non-targeting_20`) count as `NTC`. ALL `ctrl_T` (`gene_desert_*`) cells are strictly EXCLUDED from both NTC and targeting perturbations.**
4. **Target Gene Symbol Resolution & Pooling (`PREREG_E §3.1` Parity):**
   - Map all `class == "gene"` targets (`218` unique symbols × `3` guides = `654` targeting guides, all `218/218` verified HGNC-approved in `sgRNA_lib.csv.gz`) through frozen `hgnc_complete_set.txt` (`HGNCLookup.resolve`), and **pool all QC-valid single-sgRNA cells targeting the same canonical HGNC symbol (`pert_key = canon_sym`, identical to Round 12 `a549_ko`)** before Addendum-B reproducibility testing (`process_r10_line`).
5. **Cell QC Thresholds (`max_pct_mt = None`):**
   - Retain single-sgRNA cells passing BOTH **`n_umis >= 300`** AND **`n_genes >= 200`** (computed over selected raw integer mRNA UMI counts);
   - **NO mitochondrial percentage upper bound (`max_pct_mt = None`):** `GSE337804` employs Jumpcode `DepleteX` (which depletes high-abundance `MT-*` & `RPL/RPS` transcripts), so no `%MT` ceiling is imposed.

---

## §5. Disjoint Split-Half & NTC Subtraction Units (`分半和扣除单位`)

- **`Half A` (`half_1`):** **`ch_1`** (`GSE337804_A549_ch_1.RDS.h5ad`, `GSM9860781/82`);
- **`Half B` (`half_2`):** **`ch_2`** (`GSE337804_A549_ch_2.RDS.h5ad`, `GSM9860783/84`);
- **Centering / Subtraction / Batch Unit:** Individual **channel** (`batches = ["ch_1", "ch_2"]`). Each half has its own dedicated channel and its own disjoint `ctrl_NT` NTC pool (`zero` cross-half control sharing; no `<15` half-pooled fallback needed if both channels hold `>= 15` NTCs, and failure of `>= 15` NTC in either channel fails the §6 gate directly).

---

## §6. Confirmatory Count Gate (`PREREG_E §2.2` Identical)

`a549_myllia` enters confirmatory adjudication (`CONFIRMATORY_LINES_R14.txt`: `K_conf = 1`) **if and only if ALL four conditions hold simultaneously**:
1. Technical feasibility gate (§4.1: non-negative integer raw UMI counts + unambiguous single-sgRNA & `ctrl_NT`/`ctrl_T` labels across both `.h5ad` channels) passes;
2. **`n_eval_targets >= 50`** (`P5`-reproducible targets on `a549_myllia`, reproducible in `>= 1` `r9_6` source, self-gene in $G_{\mathrm{eval}}$, and `half_valid`);
3. **`|G_eval(a549_myllia)| >= 3,000`** (`|G_R9_6148 ∩ G_fin(a549_myllia)| >= 3,000`);
4. **`per-half NTC >= 15`** (`>= 15` QC-valid `ctrl_NT` singlets in `ch_1` AND `>= 15` in `ch_2`).

- **Gate Failure Consequence:** If `a549_myllia` fails ANY condition above (`EXCLUDED_FEASIBILITY`, `< 50` targets, `< 3,000` genes, or `< 15` NTC in `ch_1` or `ch_2`), **Round 14 executes ZERO confirmatory adjudication (`K_conf = 0`, 不做确认性判定)** and reports `a549_myllia` (if feasible, `SECONDARY_ONLY`) strictly descriptively under `§7` and `PREREG_E §5` (leaving Round 12 unchanged).

---

## §7. Pre-Registered Descriptive Subgroup Analysis: `R12-Overlap` vs `Novel` Targets (`只描述、不判定`)

Pre-locked **before inspecting `.h5ad` expression counts** to verify whether Round 14 `A549` transfer specificity comes solely from genes already evaluated in Round 12 `a549_ko` (`GSE345058`) or generalizes equally to unseen `A549` targets:

1. **Exact Pre-Registered Partition of the `61` `u6_rep ∩ R9_COMMON` Candidate Library Targets (Cross-Checked Against `round12/shared/lines/a549_ko/targets.txt`):**
   - **`overlap_r12_a549` (`15` candidates present in R12 `a549_ko`, all `15/15` reproducible in R12 `a549_ko`):**
     `BMS1`, `CMTR1`, `CNOT10`, `CREBBP`, `DHDDS`, `EIF4G1`, `GMPPB`, `HDAC3`, `KDM6A`, `PTEN`, `RIOK2`, `SEC61A1`, `SLC35B1`, `SMC3`, `SPCS3`.
   - **`novel_not_in_r12_a549` (`46` candidates absent from R12 `a549_ko`, `75.4%` unseen in R12):**
     `BRPF1`, `CALR`, `CBX3`, `CTCF`, `DDX1`, `DHX29`, `DHX36`, `DOT1L`, `EHMT1`, `EP300`, `EWSR1`, `FASTKD5`, `GSK3B`, `HDAC1`, `HSD17B10`, `INTS14`, `IRF2BP2`, `KAT6A`, `KDM1A`, `KDM2B`, `KDM3B`, `KEAP1`, `MSL1`, `MTPAP`, `NEDD4`, `PIK3R4`, `PRORP`, `RAF1`, `RPP14`, `RPUSD3`, `RPUSD4`, `RRP9`, `SHOC2`, `SIN3B`, `SMARCA5`, `SRF`, `SRFBP1`, `SUPT20H`, `SUPV3L1`, `TFAM`, `TFPT`, `TRMT10C`, `UBA3`, `UBA5`, `UTP6`, `ZNF598`.
2. **Pre-Specified Metrics Per Subgroup (`Strictly Descriptive / 只描述、不判定`):**
   - Using identical primary settings (`--source-mode r9_6`, `all_nontarget`, `r_orth_cal`, `n_boot = 4,000` with `seed = 20261224`, two-sided `97.5%` CIs), evaluate and report separately on `eval_targets ∩ overlap_r12_a549` and `eval_targets ∩ novel_not_in_r12_a549`:
     1. Exact target count `n_sub` and target symbols;
     2. **`FCT` `r_orth_cal` `[97.5% CI]`** (`valid_frac`, `rho_ab_orth`);
     3. **All `6` individual paired contrasts $\Delta^{\mathrm{pair}}(s') = r^{\mathrm{orth\_cal}}(\mathrm{FCT}) - r^{\mathrm{orth\_cal}}(\mathrm{NegCtrl\_FCT\_MismatchedMean}(s'))$** (`[97.5% CI]`) across `s' in r9_6` (plus descriptive `Perm Q97.5` and `E1_ConstMed`).
   - *(Note: Subgroup outputs NEVER alter Claim 1 / Claim 2 confirmatory adjudication.)*

---

## §8. Mandatory Prior Disclosures (`事先披露`)

Before downloading or inspecting `.h5ad` expression matrices:
1. **Post-R12 & Post-R13 Timing:** Round 14 (`GSE337804`) was designed after knowing Round 12 (`a549_ko` confirmed, `hek293` gate fail) AND Round 13 (`mac_ctrl` gate fail `n_eval = 3`) outcomes.
2. **`A-T1` Target Yield Estimates & Narrow Margin (`Ceiling = 61` vs Floor `= 50`):**
   - R12 `A549` `R9C` survival calibration (`445 / 452 = 98.45%`): **`60.1`** (`445 / 451 -> 60.2`);
   - R12 `A549` `u6_rep` survival calibration (`445 / 530 = 83.96%`): **`54.6`** (`445 / 529 -> 54.7`);
   - General `20%–40%` rule: **`12.2–24.4`** (fails `< 50`);
   - R12 `HEK293` calibration (`16 / 143 = 11.19%`): **`6.8`** (fails `< 50`);
   - R13 `mac_ctrl` calibration (`3 / 91 = 3.30%`): **`2.0`** (fails `< 50`).
   - With theoretical ceiling `61`, passing `>= 50` allows at most `11` target dropouts, while expected cell depth (`~318` targeting cells/gene from `~72.6k` singlets) is roughly `57%` of R12 `GSE345058`'s `565` median.
3. **Same Cell Line (`A549`):** Round 14 tests independent-dataset replication on `A549`, NOT a second lineage; cross-lineage generality remains bounded even if Round 14 confirms.
4. **Seurat `.RDS.h5ad` Raw-Count Storage Uncertainty:** Exact HDF5 slot holding integer raw UMIs (`/X`, `/raw/X`, or `/layers/*`), actual QC singlet counts, and `DepleteX` impact on `|G_eval|` vs `3,000` await post-lock verification.

---

## §9. Locked Download Plan (`下载方案`)

**Strict Execution Gate:** Zero downloads of `.h5ad` files may begin until Master Controller locks `PREREG_G.md` (`round14/status/PREREG_G_COMMIT.txt`) AND issues `round14/status/R14_GO.txt`.

1. **Exact Allowed Files (`3` files ONLY, `1,172,650,178 B` = `1.173 GB` total):**
   - **File 1 — `GSE337804_A549_ch_1.RDS.h5ad`:**
     - **URL:** `https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE337804&format=file&file=GSE337804%5FA549%5Fch%5F1%2ERDS%2Eh5ad`
     - **Verified HTTP `HEAD` `Content-Length`:** **`563,808,587` bytes** (`563.8 MB`)
   - **File 2 — `GSE337804_A549_ch_2.RDS.h5ad`:**
     - **URL:** `https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE337804&format=file&file=GSE337804%5FA549%5Fch%5F2%2ERDS%2Eh5ad`
     - **Verified HTTP `HEAD` `Content-Length`:** **`608,829,330` bytes** (`608.8 MB`)
   - **File 3 — `GSE337804_sgRNA_lib.csv.gz`:**
     - **URL:** `https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE337804&format=file&file=GSE337804%5FsgRNA%5Flib%2Ecsv%2Egz`
     - **Verified Byte Size:** **`12,261` bytes** (`12.3 KB`), pre-verified **`sha256: 96f6eb62449ba276d36591af2f7e0fca47ffedb9f19d711bcc27671d50882df6`**
2. **Storage & Verification Constraints:**
   - Save strictly under `round14/shared/raw/gse337804/`; assert exact byte lengths (`563,808,587`, `608,829,330`, `12,261`) and record all SHA-256 digests into `round14/shared/R14_CHECKSUMS.sha256` and `round14/shared/download_manifest_r14.json`;
   - Inspect and ingest exclusively via **`python3 -I`** (isolated Python runtime);
   - Never modify frozen `round10/`–`round13/` files (`round13/` freeze untouched), never register accounts or click terms, and never run `git`.
