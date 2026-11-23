# Round 13 Prospective Replication Pre-Registration Addendum (`PREREG_F.md`, `v1`)

- **Project:** AI-Biology-Lab Round 13 — Prospective Single-Line Replication Test on Unseen Primary Human Macrophage CRISPRko (`GSE327644`, Jung et al. 2026)
- **Protocol Author:** Agent `A-T3` (Team A · Pre-Registration & Statistical Audit)
- **Parent Protocol (Locked, Unmodified):** [`PREREG_E.md` (`v2`)](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/prereg/PREREG_E.md) (`commit 4c606f4732ea080768e27757c689d40d7892cd3d`, `sha256 d10f92ced16a4519c43ea2b92b8b7ee881ab3d1764e39ffd10580473e4108924`; Round 12 freeze commit `11d91c72025264b838af091abcf0eac71588eb91`)
- **Authorizing Directive & Pre-Lock Inputs:** Master Controller Directive [`R13_第十三轮阶段2_GSE327644预注册指令.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/controller/directives/R13_第十三轮阶段2_GSE327644预注册指令.md), Target Screening Report [`R13_DATA_SCREEN.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round13/shared/R13_DATA_SCREEN.md), Prior-Art Audit [`PRIOR_ART_G.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round13/shared/PRIOR_ART_G.md) + Controller Audit [`PRIOR_ART_G_REVIEW.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round13/status/PRIOR_ART_G_REVIEW.md), and Stage-2 Fallback Precedent [`RULING_R12_HEK293_UNITS.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/status/RULING_R12_HEK293_UNITS.md).

> **Strict Zero-Data-Access Blinding Statement:** As of submitting this addendum, **zero expression matrices (`GSE327644_RAW.tar`, `6,579,169,280` bytes) have been downloaded or opened**, and `round13/shared/raw/` does not exist. Only public GEO/PubMed metadata and the pre-approved `38,166`-byte sgRNA feature template (`GSE327644_PRJNA1450051_marson_mac_perturbseq_guide_template.txt.gz`, `sha256 66bc28c8f6338709a72964399861855f494aca2676f9ccd57ae27c2c7311c1fb`) have been inspected. All Round 10, Round 11, and Round 12 frozen files remain strictly unmodified.

---

## §1. Verbatim Inheritance Statement & Round 13 Seed (`沿用声明`)

`PREREG_F.md` is strictly an **addendum** to [`PREREG_E.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/prereg/PREREG_E.md). **Not a single word (`一个字都不改`)** is altered in any claim, hypothesis, predictor formula, residual endpoint, gene mask, source cohort, inference rule, multiplicity correction, eligibility criterion, or descriptive sensitivity analysis from `PREREG_E.md` (§0.1, §0.3, §2.1, §3.3, §4, §5, §7, §8), EXCEPT updating the bootstrap RNG seed and adding the dataset-specific `GSE327644` specifications below:

1. **Two Confirmatory Claims (Identical to `PREREG_E §4.2–§4.3`):**
   - **Claim 1 (`H-E1 = H-E1a ∩ H-E1b ∩ H-E1c`, Paired Target-Line Baseline Specificity):** Confirmed iff ALL three hold simultaneously on `mac_ctrl` (if it clears the §5 gate):
     - **`H-E1a`:** `FCT` achieves `r_orth_cal > 0` (`97.5% CI_lo > 0`) under `all_nontarget`, with `boot valid_frac >= 0.95` and `rho_ab_orth > 0`;
     - **`H-E1b`:** Across **each** of the `6` source lines `s' in r9_6` jointly (IUT over all `6` pairs), the joint perturbation bootstrap paired contrast $\Delta^{\mathrm{pair}}(s') = r^{\mathrm{orth\_cal}}(\mathrm{FCT}) - r^{\mathrm{orth\_cal}}(\mathrm{NegCtrl\_FCT\_MismatchedMean}(s'))$ under `all_nontarget` (sharing identical perturbation draws and denominator $c^\perp_{\mathrm{ceil}}$ within each draw) satisfies **`97.5% CI_lo > 0`**;
     - **`H-E1c`:** `FCT` `r_orth_cal` **`97.5% CI_lo > Q_{0.975}(\mathrm{NegCtrl\_FCT\_Perm}(0..99))`** under `all_nontarget`.
   - **Claim 2 (`H-E2 = H-E2a ∩ H-E2b`, Empirical Split-Half Top-50 Constant-Baseline Selection Inflation `SA-8`):** Confirmed independently of Claim 1 iff BOTH hold simultaneously on `mac_ctrl` (if it clears the §5 gate):
     - **`H-E2a`:** Under `resp_top50_split`, **`E1_ConstMed`** achieves `r_orth_cal_split` with **`97.5% CI_lo > 0`**;
     - **`H-E2b`:** Across the common perturbation set `split_valid ∩ all_valid`, the joint cross-mask perturbation bootstrap difference $\Delta^{\mathrm{cross}}_{\mathrm{SA8}} = r^{\mathrm{orth\_cal}}_{\mathrm{split}}(\mathrm{E1\_ConstMed}) - r^{\mathrm{orth\_cal}}_{\mathrm{all}}(\mathrm{E1\_ConstMed})$ (recomputing each mask's noise ceiling inside identical draws) satisfies **`97.5% CI_lo > 0`**.
2. **Predictors, Primary Endpoint, Masks, & Sources (Identical to `PREREG_E §2.1, §4.1`):**
   - Primary endpoint: **`r_orth_cal`** (using `Transfer_Mean` reference, `orthogonalize_to_ref`, and Spearman–Brown split-half ceiling disattenuation clipped to `[-1, 1]`). Note carefully that **`E1_ConstMed`** ($\hat{\Delta}^{\mathrm{E1}}_{t,p,g} = \operatorname{log1p}(c_{\mathrm{med}}\overline{\mathrm{FC}}_{p,g}) - \operatorname{log1p}(c_{\mathrm{med}})$, where $c_{\mathrm{med}} = \operatorname{median}_{g \in G_{\mathrm{eval}}(t)} \max(\mu_{t,g}, 0)$) is strictly a **target-baseline-free transfer predictor** (`不使用目标线基线的迁移预测器` — it still depends on source-line fold changes $\overline{\mathrm{FC}}_{p,g}$ and must NEVER be described as "target-independent").
   - Masks: **`all_nontarget`** ($G_{\mathrm{eval}}(t) \setminus \{g_p\}$) for Claim 1 and the `H-E2b` whole-transcriptome comparator; **`resp_top50_split`** (disjoint symmetric `A->B` / `B->A` top-50 responding mask) for Claim 2; gene universe $G_{\mathrm{eval}}(t) = G_{\mathrm{R9\_6148}} \cap G_t^{\mathrm{fin}}$.
   - Training source cohort: Strictly the **`6`** Round 9 harmonized human **`CRISPRi`** cell lines (**`r9_6`**: `nadig_hepg2`, `nadig_jurkat`, `vcc2025_h1`, `replogle_k562_essential`, `replogle_rpe1`, `nourreddine_kolf21j`; `6,148` `R9_COMMON` genes; `2,117` line-reproducible `u6_rep` symbols). Sensitivity pool **`SA-Src8`** (`--source-mode r9_r10_8`, adding Round 10 `dld1` & `cd4t_rest`) remains strictly descriptive in `§5` (`只报告、不判定`). Neither macrophage cohort exists in `r9_6` or `r9_r10_8`.
3. **Inference, Resampling, Addendum-B `P5`, & New Seed:**
   - Multiplicity: **Intersection–Union Test (IUT) + 2-Claim Bonferroni** ($\alpha_{\mathrm{claim}} = 0.025$ per claim, two-sided **`97.5%` percentile bootstrap confidence intervals `[1.25%, 98.75%]`**).
   - Resampling: **`n_boot = 4,000`** joint perturbation bootstraps using **new Round 13 RNG seed `20261223`** (`np.random.default_rng(20261223)`), and **`n_perm = 100`** (`k = 0..99`, `np.random.default_rng(k)`) gene-label permutations on $\mu_t$.
   - Standardization & `P5` target eligibility: Identical to `PREREG_E §3.3` (library-size $10^4$ normalization, unit-matched NTC `log1p` subtraction for $\Delta^{\mathrm{obs}}_{t,p}$, Addendum-B reproducibility `n_cells >= 25`, `>= 50` top-5,000 genes with two-sample Anderson–Darling `BH FDR < 0.05` vs `1,000` stratified NT, `20-PC` energy-distance `p < 0.001` over `2,000` perms vs `500` NT, reproducible in `>= 1` `r9_6` source, self-gene in $G_{\mathrm{eval}}(t)$, `half_valid` `>= 2` cells in Half A AND `>= 2` in Half B).
   - All `PREREG_E §5` descriptive analyses (`H-D2` unpaired, `FCT - E1_ConstMed`, `N_5` semi-partial `r_cone_sp` including construction-null notes, `SA-Src8`, in-sample `resp_top50`, `H-D_util`, `SA-1`, `SA-4/R19`, `SA-5` full-NT sensitivity, `SA-6/S7`), clean-room tolerances (`max_point_abs_diff <= 1e-10`, `max_ci_abs_diff <= 1e-10`, `100%` verdict parity), and `PREREG_E §0.3`/`§8` prior-art locks remain completely unchanged:
     - **Kawsar et al. (2026)** (`CanPert/CASTER`, `doi:10.64898/2026.11.15.766815`) for prior per-target paired correct-vs-mismatched control comparison;
     - **Nicol, Shivakumar & Irizarry (2026)** (`doi:10.64898/2026.05.07.723486`) for finite-sample DEG selection inflation (explicitly noting Claim 2 here applies strictly to `r_orth_cal` `resp_top50_split`, with **unknown mechanism** — Round 11 noiseless simulation failed sign on `8/12` lines under `all_nontarget` — and does NOT extrapolate to DEG metrics);
     - **Liu et al. (2026)** (`GSE345058`, `doi:10.64898/2026.06.01.728087`) as the `A549` data source in Round 12 (never misattributed);
     - **Jung et al. (2026)** (*Nat Biotechnol*, Epub `2026-08-17`, PMID `42608566`, `doi:10.1038/s41587-026-03258-2`, BioProject `PRJNA1450051`) strictly as the `GSE327644` dataset source (zero prior cross-line transcriptomic prediction work found on `GSE327644` per [`PRIOR_ART_G.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round13/shared/PRIOR_ART_G.md)). Zero use of `X-Atlas/Orion` or `Arc Virtual Cell Challenge (VCC)`.

---

## §2. Line Definitions (`GSE327644`, Primary Human `CD14+` Macrophages)

Upstream GEO metadata (`GSM9662558–GSM9662561`) profile M-CSF-differentiated primary human `CD14+` monocyte-derived macrophages from donor `RG5581`, transduced with a VPX-lentiviral sgRNA library (`583` myeloid-relevant genes, `4` guides/gene + `593` `Control` guides) and Cas9-eVLP ribonucleoprotein delivery (`SLICeVLP` **CRISPR knockout [`CRISPRko`]**, noting GEO's structured submission field `genetic perturbation_strategy: CRISPR interference` is a documented metadata artifact contradicted by `!Sample_title`, `!Series_summary`, and `!Sample_extract_protocol_ch1`), profiled via 10x **GEM-X Flex v2** across `4` channels (`Mac_L01–Mac_L04`), processed by Cell Ranger v9.0.1 (`Chromium Human Transcriptome Probe Set v1.1.0: 2024-A`, GRCh38), and demultiplexed per channel by Flex probe barcodes into two conditions (`8` `.h5` files):

| Line Key | Condition & Channels (GEO GSMs) | Flex Probe Barcodes | Pre-Locked Role Rule |
|---|---|---|---|
| **`mac_ctrl`** | **Unstimulated macrophages (`Mac_Ctrl`)**, channels **`L01–L04`** (`GSM9662558`–`GSM9662561`, `*_Mac_Ctrl_sample_filtered_feature_bc_matrix.h5`) | `BC001+CR001` .. `BC008+CR008` | **SOLE confirmatory candidate line (`唯一的确认性候选线`)** — enters confirmatory adjudication iff §5 technical feasibility + count gates pass (`K_conf = 1`). |
| **`mac_lps`** | **LPS-stimulated macrophages (`Mac_LPS`, `50 ng/ml` ultrapure LPS for `5 h`)**, channels **`L01–L04`** (`GSM9662558`–`GSM9662561`, `*_Mac_LPS_sample_filtered_feature_bc_matrix.h5`) | `BC009+CR009` .. `BC016+CR016` | **Strictly descriptive secondary line (`只作描述性的次要线`)** — **NEVER enters confirmatory adjudication** regardless of cell/target counts or scores; processed identically if technically feasible and reported strictly under `§5`. |

---

## §3. Cell, Single-sgRNA, Locked NTC, Target Symbol, & QC Rules (`细胞和 guide`)

1. **Cell Ranger Filtered Input Matrix:**
   - Ingest strictly the Cell Ranger v9.0.1 **filtered** feature-barcode matrices (`*sample_filtered_feature_bc_matrix.h5`; `raw_matrix_background_barcodes_included: False`) per channel (`L01–L04`) and condition (`Mac_Ctrl`, `Mac_LPS`). Never switch to unfiltered/raw background barcodes.
   - Retain strictly `Gene Expression` features (`feature_type == "Gene Expression"`) in the expression matrix (casting HDF5 `indptr`/`data` to `int64` to prevent `int32` overflow per `DEVIATIONS_R12.md D2`), and map expression gene features through frozen `hgnc_complete_set.txt` identically to `process_r10_line`.
2. **Single-sgRNA Cell Calling (`Exact 1 Guide`):**
   - **Primary rule (per directive):** Call single-sgRNA identity from Cell Ranger's protospacer call (`protospacer_calls_per_cell.csv` if bundled in `GSE327644_RAW.tar`, or an embedded per-cell protospacer call field in `.h5`), keeping **ONLY** barcodes whose Cell Ranger protospacer call has **strictly `1` guide** (`num_features == 1` and its single `feature_call` matches a valid `guide_id` in `guide_template.txt.gz`).
   - **Pre-registered deterministic fallback (applied ONLY if `GSE327644_RAW.tar` ships exclusively the `8` `_sample_filtered_feature_bc_matrix.h5` count matrices without `protospacer_calls_per_cell.csv` or embedded call fields — matching the upstream Lattice/Marson GEM-X Flex `cd4t_rest` schema in `round10/shared/scripts/prepare_r10_lines.py`):** Call per-guide positivity within each line (`mac_ctrl` and `mac_lps` independently, across cells with `>= 200` detected genes) via the frozen Round 10 Poisson-Gaussian histogram threshold estimator (`fit_poisson_gaussian_thresholds_from_hist` with `min_thr=5, max_thr=20, default_thr=7` on non-`Gene Expression` CRISPR guide UMI counts), and retain **ONLY** cells with **strictly `1` positive guide (`n_pos == 1`)**, **`top1_UMI >= 5`**, and **`(top2_UMI < 3 or top1_UMI / max(top2_UMI, 0.5) >= 3.0)`**. Which branch executed must be logged in `round13/shared/R13_QC_COUNTS.md` prior to computing any downstream evaluation correlation.
3. **Locked NTC Definition (Set Pre-Lock From `guide_template.txt.gz` Inspection):**
   - Direct inspection of `GSE327644_PRJNA1450051_marson_mac_perturbseq_guide_template.txt.gz` (`38,166` bytes, `sha256 66bc28c8f6338709a72964399861855f494aca2676f9ccd57ae27c2c7311c1fb`, `2,925` data rows across `6` columns `guide_id`, `guide_protospacer`, `guide_role`, `guide_target_gene_id`, `guide_PAM`, `guide_target_gene_name`) establishes:
     - Exactly `2` `guide_role` categories exist: **`'Targeting a Gene'`** (`2,332` rows) and **`'Control'`** (`593` rows);
     - **ALL `593` `'Control'` rows are uniformly named `NTC-001` through `NTC-593`** (`guide_target_gene_id == 'n/a'`, `guide_target_gene_name == 'n/a'`), with **zero distinction between `non-targeting` and `safe-targeting` controls** (and `0` `safe`/`AAVS1`/`OR`/intergenic control guides anywhere in the template).
   - **Pre-Locked NTC Rule:** Because the guide template does not separate `non-targeting` from `safe-targeting`, **ALL `593` `Control` guides (`NTC-001`..`NTC-593`) count as `NTC`**. Cells assigned to targeting guides, multiple guides, or zero guides are strictly excluded from `NTC`.
4. **Locked Target Gene Symbol Canonicalization (Guide Template Audit):**
   - Upstream inspection of `guide_template.txt.gz` revealed two export artifacts if naively trusting `guide_target_gene_name`: (a) `guide_target_gene_name` truncated all `19` multi-hyphen `guide_id`s at the first hyphen (collapsing `14` distinct `HLA-*` genes `HLA-A..HLA-G` into invalid symbol `"HLA"`, and collapsing readthrough targets `ARPC4-TTLL3`, `KLRC4-KLRK1`, `PPAN-P2RY11`, `TMED7-TICAM2`, `TMEM56-RWDD3` into `ARPC4`, `KLRC4`, `PPAN`, `TMED7`, `TMEM56`); (b) Excel auto-formatted `JUN-1..JUN-4` into `guide_id`s `'1-Jun'`..`'4-Jun'` (`guide_target_gene_name == 'JUN'`, `ENSG00000177606`).
   - **Pre-Locked Mapping Rule:** For each `'Targeting a Gene'` row, form raw target stem `stem = guide_id.rsplit("-", 1)[0]` (remapping `{"1", "2", "3", "4"}` where `guide_target_gene_name == "JUN"` to `"JUN"`; likewise matching `'1-Jun'..'4-Jun'` if present in `.h5` `features/name`), then resolve `(stem, guide_target_gene_id)` through frozen `HGNCLookup.resolve` (`approved -> ensg -> prev_symbol -> alias_symbol`). This restores all **`583`** target genes described by Jung et al. (`582` HGNC-approved + `CRIPAK` unmapped, preserving `HLA-A..G` and readthrough genes distinctly without false collision of `PPAN-P2RY11` into `PPAN`).
5. **Cell QC Thresholds (`No %MT Ceiling`):**
   - Retain single-sgRNA cells passing BOTH **`n_umis >= 300`** AND **`n_genes >= 200`** (computed over `Gene Expression` UMI counts).
   - **NO mitochondrial percentage upper bound (`max_pct_mt = None`):** GEM-X Flex probe panels (`2024-A`) may omit `MT-` probes; consistent with `PREREG_E` HEK293 (`DEVIATIONS_R12.md D1`), no `%MT` ceiling is imposed.

---

## §4. Disjoint Split-Half & NTC Subtraction Units (`分半和扣除单位`)

1. **Stage 1 (Default Primary — Channel Parity & Channel-Matched Subtraction):**
   - Both `mac_ctrl` and `mac_lps` span `4` independent 10x Flex channels (`L01`, `L02`, `L03`, `L04` = `GSM9662558..GSM9662561`):
     - **`Half A` (`half_1`):** **`L01` + `L03`** (`GSM9662558`, `GSM9662560`);
     - **`Half B` (`half_2`):** **`L02` + `L04`** (`GSM9662559`, `GSM9662561`).
   - **Centering / Subtraction / Batch Unit:** Each individual **channel** (`L01`, `L02`, `L03`, `L04`), matching `a549_ko` `Lane`-matched NTC subtraction in `PREREG_E §3.1/§3.3`. Because each channel lies strictly within one half, zero NTC cells are shared across halves.
   - **Channel `< 15` NTC Fallback Rule:** Using the frozen Round 11/12 `_apply_min_ntc_batch_fallback(..., min_batch_nt=15)` semantics prior to `process_r10_line`, if any channel (`L01..L04`) contains **`< 15`** QC-valid single-sgRNA NTC cells, its cells (NTC + perturbed) merge for NTC mean/SD subtraction into that channel's **own half-pooled NTC bucket (`half_1_pooled` for `L01/L03`; `half_2_pooled` for `L02/L04`)** — NEVER across halves.
2. **Stage 2 (Locked Fallback ONLY When Channel Metadata Is Missing):**
   - Invoked **strictly and exclusively** if per-cell channel labels (`L01–L04`) cannot be identified from filenames or metadata (not expected given GEO GSM `Mac_L01..L04` filenames). Follows **`RULING_R12_HEK293_UNITS.md §2`** verbatim:
     1. Instantiate **one** RNG: `rng = np.random.default_rng(42)`;
     2. Define stratum key as the canonical HGNC target symbol for perturbed cells and literal `"__NTC__"` for NTC cells; iterate over strata sorted in Python default lexicographical string order;
     3. Within each stratum, sort cells in ascending order of `(true cell_barcode, kept_row_idx)`;
     4. Draw `perm = rng.permutation(n)` once per stratum and assign `perm[:n//2]` -> **`Half A` (`half_1`)** and `perm[n//2:]` -> **`Half B` (`half_2`)**, setting centering/batch unit `batches = ["half_1", "half_2"]` so NTC subtraction occurs strictly within each half.
   - Which stage executed (`Stage 1` vs `Stage 2`, plus any `< 15` NTC channel fallbacks) is recorded deterministically into `round13/shared/R13_QC_COUNTS.md` before unblinding.

---

## §5. Pre-Metric Technical Feasibility & Confirmatory Count Gate (`PREREG_E §2.2`)

1. **Pre-Metric Technical Feasibility Check (`EXCLUDED_FEASIBILITY`):**
   - Immediately upon ingestion of `GSE327644_RAW.tar` (before computing any predictor correlation or metric), `A-T1` verifies on raw counts and barcodes:
     1. Non-negative integer raw UMI counts exist in the `Gene Expression` HDF5 matrix (`np.all(data >= 0)` and integer-valued);
     2. Per-cell guide assignment succeeds (either Cell Ranger protospacer calls or CRISPR Guide Capture UMI rows matching `guide_template.txt.gz`).
   - If EITHER condition fails on `mac_ctrl` (or `mac_lps`), that line is classified **`EXCLUDED_FEASIBILITY`** with exact technical failure logs recorded in `round13/shared/R13_QC_COUNTS.md` and `CONFIRMATORY_LINES_R13.txt`, entering neither confirmatory nor secondary metrics.
2. **Strict Confirmatory Count Gate (`mac_ctrl` Only):**
   - `mac_ctrl` clears confirmatory eligibility (`CONFIRMATORY_LINES_R13.txt`: `K_conf = 1`) **if and only if ALL THREE conditions hold simultaneously**:
     1. **`n_eval_targets >= 50`** (`P5`-reproducible targets on `mac_ctrl`, reproducible in `>= 1` `r9_6` source, self-gene in `G_eval(mac_ctrl)`, and `half_valid`);
     2. **`|G_eval(mac_ctrl)| >= 3,000`** (`|G_R9_6148 ∩ G_fin(mac_ctrl)| >= 3,000`);
     3. **`per-half NTC >= 15`** QC-valid single-sgRNA NTC cells in BOTH `Half A` (`L01+L03`) and `Half B` (`L02+L04`).
3. **Consequence of Gate Failure (`Early Stop at Gate`):**
   - If `mac_ctrl` fails **ANY** of the three gates above (`n_eval_targets < 50`, `|G_eval| < 3,000`, or `< 15` NTC in either half), **Round 13 performs ZERO confirmatory adjudication (`K_conf = 0`, 不做确认性判定)**, terminates at descriptive counts/metrics, and reports `mac_ctrl` (if technically feasible, `SECONDARY_ONLY`) alongside `mac_lps` strictly descriptively (`只作描述性报告`, Round 12 `A549` conclusions remain unchanged).

---

## §6. Pre-Locked Interpretive Rules (`写在看数据之前`)

1. **Strict `max K_conf = 1` & Ban on Cross-Round Pooling:**
   - Within Round 13 itself, maximum `K_conf = 1` (`mac_ctrl` only).
   - **If `mac_ctrl` passes the §5 count gate AND BOTH Claim 1 (`H-E1a/b/c`) and Claim 2 (`H-E2a/b`) pass (`97.5% CI_lo > 0` across all components):** Results MUST ONLY be stated as **"confirmed independently on two separate cell lines in two separately preregistered prospective tests"** (`「在两轮分别预注册的前瞻检验中、两条独立细胞系上各自确认」`: Round 12 `A549` `PREREG_E` and Round 13 `mac_ctrl` `PREREG_F`). **Round 12 and Round 13 data MUST NEVER be pooled into any joint hypothesis test or joint confidence interval.**
   - **If `mac_ctrl` passes the §5 count gate AND EITHER Claim 1 or Claim 2 fails:** Report honestly and prominently that the unconfirmed claim **failed to replicate on the second cell line (`在第二条细胞系上未能复现`)**, noting Round 12 `A549` results remain unchanged and disclosing confounding between macrophage biology, Flex probe panel, single donor (`RG5581`), CRISPRko modality, and statistical power.
   - **If `mac_ctrl` fails the §5 gate (`K_conf = 0`):** Report honestly as a prospective gate failure / descriptive-only result (`不做确认性判定`), leaving Round 12 `A549` (`K_conf = 1`, limited cross-line generalizability) unchanged.
2. **Mandatory Invariant Reporting (`PREREG_E §7`):** Every Round 13 report must emit full line-by-line numbers (`estimate [97.5% CI]`, `valid_frac`, `rho_ab_orth`, all `6` individual `H-E1b` pairs, `H-E1c`, `H-E2a`, `H-E2b`), modest whole-transcriptome $R^2$, complete `§5` descriptive metrics (`H-D2`, `FCT - E1`, `N_5`, `SA-Src8`, `mac_lps`, `S7`), and all non-mechanistic/prior-art boundaries (`E1_ConstMed` target-baseline-free wording; Claim 2 `SA-8` unknown mechanism & non-DEG scope; Kawsar et al. 2026; Nicol et al. 2026; Jung et al. 2026; Liu et al. 2026).

---

## §7. Mandatory Prior Disclosures (`事先披露`)

Before downloading or opening any expression matrix from `GSE327644`, the following five risks and provenance facts are permanently on record:

1. **Post-Round-12 Timing:** Round 13 (`PREREG_F.md`) was conceived and designed **after Round 12 (`PREREG_E`) results were already known**. Systematic candidate screening across public portals ([`R13_DATA_SCREEN.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round13/shared/R13_DATA_SCREEN.md)) yielded no safe large-scale second line and identified `GSE327644` as the sole borderline candidate worth testing prospectively under locked count gates.
2. **Metadata-Only Prospective Target Yield Estimates & Pre-Recorded High Failure Probability:**
   - Using raw `guide_target_gene_name` from `guide_template.txt.gz`, `u6_rep` contains `108` symbols (`92` inside `R9_COMMON`); after fixing the `-`-truncation artifact (`§3.4`), `u6_rep` contains `107` symbols (`91` inside `R9_COMMON`, dropping only the false collision `PPAN-P2RY11 -> PPAN`).
   - Applying `A-T1`'s pre-download `R13_DATA_SCREEN.md` empirical yield rates to `108` `u6_rep` (`92` `R9_COMMON`):
     - `HEK293` calibration (`9.76%`): **`10.5`** (`10.3` on `R9_COMMON`) — fails `>= 50`;
     - General `20%–40%` rule of thumb: **`21.6–43.2`** (`18.4–36.8` on `R9_COMMON`) — fails `>= 50`;
     - `A549` calibration (`83.96%`): **`90.7`** (`90.6` on `R9_COMMON`) — clears `>= 50`.
   - Structurally, `GSE327644` (`4` channels enabling channel-parity splitting and channel-level NTC subtraction) avoids the `HEK293` single-pool failure mode and matches `A549`'s lane structure, but its immune-focused `583`-gene library provides narrower `u6_rep` overlap (`107`/`108` vs `A549`'s `530`). **Master Controller and Team A record beforehand that failure probability at the `>= 50` count gate (or in statistical power) is high.**
3. **Single Biological Donor (`RG5581`):** All macrophages derive from one healthy 27-year-old European American male donor (`RG5581`), restricting biological donor diversity.
4. **10x GEM-X Flex Targeted Probe Panel (`2024-A`):** Fixed probe hybridization (`Chromium Human Transcriptome Probe Set v1.1.0`) differs from whole-transcriptome poly-A capture and may alter `G_fin ∩ R9_COMMON` (`|G_eval|` vs `3,000` floor), target self-gene coverage, or `MT-` gene inclusion.
5. **Undisclosed MOI & Single-sgRNA Cell Counts:** Library transduction MOI, single-sgRNA fraction vs multi-guide/unassigned droplets, per-channel sequenced singlets (`L01–L04`), per-channel NTC counts, and biological perturbation effect sizes (`ad_counts >= 50`) in unstimulated macrophages remain unknown prior to download.

---

## §8. Locked Download Plan (`下载方案`)

**Strict Execution Gate:** Zero downloads of `GSE327644_RAW.tar` may begin until Master Controller locks `PREREG_F.md` (recording commit hash in `round13/status/PREREG_F_COMMIT.txt`) AND issues `round13/status/R13_GO.txt`.

1. **Exact Allowed Files (`2` files ONLY, `6,579,207,446 B` = `6.579 GB` total; zero extra series/raw fastq downloads):**
   - **File 1 — `GSE327644_RAW.tar`:**
     - **URL:** `https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE327644&format=file&file=GSE327644%5FRAW%2Etar` (HTTP HEAD redirect target: `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE327nnn/GSE327644/suppl/GSE327644_RAW.tar`)
     - **Verified Upstream Byte Size (`Content-Length`):** **`6,579,169,280` bytes** (`~6.579 GB` / `6.127 GiB`)
     - **Verification Rule:** Assert downloaded file size equals `6,579,169,280` bytes cleanly; compute and log full SHA-256 of `GSE327644_RAW.tar` AND every extracted member file (`GSM9662558..GSM9662561` `.h5` matrices and any bundled `protospacer_calls_per_cell.csv`) into `round13/shared/R13_CHECKSUMS.sha256` and `round13/shared/download_manifest_r13.json`.
   - **File 2 — `GSE327644_PRJNA1450051_marson_mac_perturbseq_guide_template.txt.gz`:**
     - **URL:** `https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE327644&format=file&file=GSE327644%5FPRJNA1450051%5Fmarson%5Fmac%5Fperturbseq%5Fguide%5Ftemplate%2Etxt%2Egz`
     - **Verified Upstream Byte Size:** **`38,166` bytes** (`38.2 KB`)
     - **Pre-Verified SHA-256:** **`66bc28c8f6338709a72964399861855f494aca2676f9ccd57ae27c2c7311c1fb`** (verify exact byte + SHA-256 match upon archiving into `round13/shared/raw/` post-`R13_GO.txt`).
2. **Storage & Ingestion Constraints:**
   - Archive both raw downloads strictly under `round13/shared/raw/gse327644/`;
   - Read and process all raw `.tar`/`.h5`/`.gz` files strictly with **`python3 -I`** (isolated Python runtime);
   - Never register external accounts, click clickwrap agreements, or contact external authors/registries; never modify frozen `round10/`, `round11/`, or `round12/` files.

---

## §9. Controller Lock Note (added by Claude, master controller, 2026-11-23 08:30 +08:00, before lock and before any expression data; no claim, test, threshold, metric, mask or source changed)

1. **Channel parity is defined over all channels present, not a fixed list.** The R12 screen recorded `Mac_L01..L08` / 16 `.h5` files, while the R13 screen and GEO GSM list give `Mac_L01..L04` / 8 files. To remove this ambiguity before data: sort the channel labels found in `GSE327644_RAW.tar` by channel number; odd-numbered channels (L01, L03, L05, …) → Half A (`half_1`), even-numbered (L02, L04, L06, …) → Half B (`half_2`); every channel present is used; the centering unit is the individual channel with the `< 15` NTC half-pooled fallback of §4.1. The number of channels found is logged in `R13_QC_COUNTS.md`. With exactly 4 channels this is identical to §4.1.
2. **Guide counts** are read from the non-`Gene Expression` feature rows of each `.h5` (e.g. `CRISPR Guide Capture`), matched to `guide_template.txt.gz` `guide_id` (including the `1-Jun`..`4-Jun` artifacts of §3.4). Which source of guide calls was used (Cell Ranger protospacer calls or the §3.2 fallback) is logged before any metric.
3. **Gate failure still allows the descriptive report** (§5.3): `mac_ctrl` as `SECONDARY_ONLY` and `mac_lps` are evaluated with the same frozen code and reported descriptively; no confirmatory language.
4. **Disk and data handling:** the raw tar (6.58 GB) and extracted members go under `round13/shared/raw/gse327644/` and are never committed to git; their sha256 go into the round-13 freeze manifest.
