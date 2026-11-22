# Round 11 Confirmatory Pre-Registration Protocol (`PREREG_D.md`, v2 R1-Revised)

- **Project:** AI-Biology-Lab Round 11 — Prospective Zero-Shot Verification of Whole-Transcriptome De-Amplituded Fold-Change Transfer (`FCT`) and Confirmatory Audit of the `resp_top50` Selection Artifact (`PREREG_D`)
- **Protocol Author:** Agent `A-T3` (Team A · Pre-Registration & Statistical Audit), revised per Master Controller Directive [`R11_第十一轮行动指令.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/controller/directives/R11_第十一轮行动指令.md), Pre-Lock Audit [`PREREG_D_AUDIT_R1.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round11/status/PREREG_D_AUDIT_R1.md) (`Q1–Q10`), Verified Download Plan [`R11_DOWNLOAD_PLAN.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round11/shared/R11_DOWNLOAD_PLAN.md), and Literature Audit [`PRIOR_ART_D.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2/vcc_linespecific/audit/LIT/PRIOR_ART_D.md).
- **Revision History:**
  - `v1` (`sha256: 87e428b48ee7d87d0ecd90887d0a0e77a664cda764bb5ac4d1fbdc418b0c9f1f`): Initial prospective protocol.
  - **`v2 (R1-Revised, this document)`:** Incorporates all `Q1–Q10` audit revisions (`97.5%` Bonferroni CIs via IUT + 2-claim familywise control, `80.000 GB` cap + `hcasmc` pre-eval feasibility gate, removal of all fallback phrasing, verified `3 + 1` dataset manifest including `feng_ipsc` `17`-line/`9`-donor non-KOLF scope, strict NTC-only vs `safe-targeting`/`unassigned` exclusions, `cd4t_stim48` disjointness & independence disclosures, *Nicol, Shivakumar & Irizarry (2026)* priority lock, `SA-8` `resp_top50_split`, descriptive `H-D5` for `feng_ipsc`, recomputed `§6` power numbers under $z_{97.5} = 2.2414$, and verbatim biological claim wording & abstract disclosures).

---

## §0. Scientific Scope, Exact Primary Claim Wording (`Q10`), Prior-Art Lock (`Q9`), & Blinding Guarantee

### 0.1 Exact Primary Biological Claim & Magnitude Limitation (`Q10`)

Round 11 tests prospectively on completely unseen Perturb-seq / CRISPRi datasets whether whole-transcriptome `FCT` (`all_nontarget` evaluation mask, free of target-magnitude gene selection) carries genuine target-baseline-driven, line-specific signal beyond amplitude rescaling.

The exact pre-registered Primary Biological Claim is stated verbatim as:

> **「在全基因口径下，用目标线自身基线把跨线平均倍数变化映射回 `log1p` 空间（`FCT`），能预测一部分去幅度后的线特异响应，并且胜过错配基线、置换和常数基线。」**
>
> *(Under whole-transcriptome evaluation, mapping cross-line average fold changes back to `log1p` space using the target line's own baseline (`FCT`) predicts a modest portion of de-amplituded line-specific response, outperforming mismatched baselines, permutations, and constant baselines.)*

**Explicit Variance-Explained Limitation (`Q10`):** In Round 10-E exploratory whole-transcriptome (`all_nontarget`) ablation analyses ([`AUDIT_R10E_REPORT.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round10/audit/exploratory/AUDIT_R10E_REPORT.md)), `FCT` achieved de-amplituded noise-corrected correlation $r_{\mathrm{orth\_cal}} \in [+0.0556, +0.1848]$ across unblinded R9/R10 lines and **`[+0.0556, +0.1355]` (`[0.056, 0.136]`)** on Round 10's two new confirmatory lines (`cd4t_rest`, `dld1`), corresponding to only **$R^2 \sim 0.3\%\text{–}1.8\%$** of de-amplituded line-specific variance. Round 11 **does not claim large variance explained**; rather, it tests rigorously whether this modest whole-transcriptome baseline-dependent signal replicates out-of-sample above negative controls and target-agnostic nonlinear scaling.

### 0.2 Literature Attribution & Priority Lock (`Q9`)

- **Base Method (`FCT`):** Unchanged from `PREREG_A v4 §4.6` / `PREREG_C §3.1`. Multiplicative target-baseline fold-change scaling (`mu_t * FC_bar`) already exists in VCC 2026 #1 (`Illumina/BioInsight CellularAI`) and #82 (`AtlasShift`); **zero novelty is claimed for the `FCT` predictor formula.**
- **Selection-Induced Shared-Control Inflation (`Nicol, Shivakumar & Irizarry 2026`):** Round 11 explicitly cites **Nicol, Shivakumar & Irizarry (2026)** (`doi:10.64898/2026.05.07.723486`), which established that selecting top-$k$ differentially expressed / largest-magnitude genes from shared control samples inflates downstream effect and correlation metrics through finite-sample shared-control noise. **Round 11 does NOT claim priority on DEG/top-$k$ selection inflation.**
- **Scoped Methodological Question of `H-D4` (`§4.2`) and `SA-8` (`§5`):** Building directly on that mechanism, Round 11 tests specifically (`H-D4` in `§4.2`) whether the Perturb-seq `resp_top50` mask creates a large artificial boost (`E1_ConstMed(resp_top50) > 0` and `E1_ConstMed(resp_top50) - E1_ConstMed(all_nontarget) > 0`) **even for a deterministic, gene-invariant constant-baseline nonlinear `log1p(c * FC)` scaling predictor (`E1_ConstMed`)**—alongside `Perm100` null elevation—and (`SA-8`, `resp_top50_split`) quantifies how `resp_top50` correlation behaves when top-50 gene selection is performed independently across disjoint split halves (`Half A` selects top 50, `Half B` evaluates, averaged with `B -> A`).
- **Cell-Count & Knockdown-Depth Confounding:** Preserves attribution to **Shoeibi & Yousefi (2026)** (`arXiv:2608.00152v3`) via pre-specified secondary analysis `SA-6 / S7`.

### 0.3 Strict Zero-Data-Access Blinding Guarantee

As of locking `PREREG_D.md v2`:
1. **Zero Round 11 expression matrices** have been downloaded, opened, or inspected (`round11/shared/raw/` and `round11/shared/lines/` do not exist).
2. All dataset eligibility, `P4` single-sgRNA definitions, NTC-only rules, split-half partitions, and `kolf_2`/`kolf_3` and `cd4t_rest` disjointness checks were verified solely from lightweight public API/repository metadata JSON/TSVs and existing Round 10 cached CSVs recorded in `R11_DOWNLOAD_PLAN.md`.
3. Every `.py` file listed in `round10/status/FREEZE_MANIFEST_R10.txt` remains verified intact (`22/22 OK`, zero modifications).

---

## §1. Strict Six-Step Execution Gate (`R11` Order Lock)

No downstream step may begin until its preceding gate file exists in `round11/status/`:

1. **Step 1 (Pre-Registration Revision & Code Parity Alignment):** `A-T3` submits `PREREG_D.md v2 (R1-Revised)` and `COMMIT_REQUEST_A-T3.txt`. Prior to `CODE_FREEZE_R11`, `A-T2` aligns `round11/src/` with `v2` (`97.5%` CIs, `SA-8` `resp_top50_split`, `H-D5` `nourreddine_kolf21j` on `feng_ipsc`, etc.) and `A-T3` independently implements `round11/audit/`, testing strictly against Round 9/10 unblinded lines and synthetic fixtures.
2. **Step 2 (`R11_GO.txt` Issuance):** Master Controller (`Claude`) verifies `PREREG_D.md v2` and `R11_DOWNLOAD_PLAN.md`, locks the git commit hash into `round11/status/PREREG_D_COMMIT.txt`, and issues `round11/status/R11_GO.txt`. **Zero Round 11 data files may be downloaded before `R11_GO.txt` exists.**
3. **Step 3 (`A-T1` Download, Pre-Eval `hcasmc` Feasibility Gate, & Counts-Only QC):**
   - `A-T1` downloads only the locked dataset files in `§2.2` strictly within the **`80.000 GB` (`80,000,000,000 B`)** Round 11 download cap (`76.254 GB` core / `76.338 GB` including optional Feng Guide-UMI & cached CD4i CSVs).
   - **Pre-Evaluation Technical Feasibility Gate (`Q3`, `hcasmc` only):** Immediately upon completing download and **strictly BEFORE computing any expression correlation, `r_orth_cal`, or `CONFIRMATORY_LINES_R11.txt`**, `A-T1` executes a bounded-memory HDF5 read + guide-capture extraction test on the `36 GiB` host (see `§2.2.2`). Any technical failure triggers immediate **pre-evaluation feasibility exclusion (`可行性排除`)** with root cause recorded in `R11_QC_COUNTS.md`.
   - **`Q6` Pre-QC Disjointness Enforcement:** Prior to computing Addendum-B QC, verify/enforce (i) `cd4t_stim48` (`CE0008678`, `CE0010866`, `Stim48hr`) has zero donor/lane/timepoint overlap with Round 10 `cd4t_rest` (`CE0006864`, `CE0008162`, `Rest`), and (ii) strict exclusion of `kolf_2` and `kolf_3` from `feng_ipsc`.
   - `A-T1` emits **counts-only** metadata in `round11/status/R11_QC_COUNTS.md` (cell counts, NTC counts, finite aligned genes, reproducible targets, shared targets, split-half valid counts; **strictly zero** predictor correlations or evaluation outcomes).
4. **Step 4 (Mechanical Confirmatory Cohort Assignment):** Apply `§2.3`'s deterministic thresholds (`candidate in {telohaec, hcasmc, cd4t_stim48}`, `passed_pre_eval_feasibility == True`, `|G_eval(t)| >= 3,000`, `n_eval_targets >= 50`) mechanically from `R11_QC_COUNTS.md` to emit `round11/status/CONFIRMATORY_LINES_R11.txt` without post-hoc overrides.
5. **Step 5 (`CODE_FREEZE_R11.txt`):** Master Controller audits `round11/src/` and `round11/audit/` parity against `PREREG_D.md v2` (`97.5%` CIs, `SA-8`, `H-D5`) and locks `round11/status/CODE_FREEZE_R11.txt`. Runtime gate `check_prereg_d_gate()` blocks any `--phase holdout_r11` run unless both `PREREG_D_COMMIT.txt` and `CODE_FREEZE_R11.txt` exist and are non-empty.
6. **Step 6 (Unblinding, Dual-Team Independent Evaluation, & `C3` Cross-Check):** `A-T2` executes primary evaluation into `round11/results/` and `A-T3` executes fully independent recomputation + `C3` cross-check into `round11/audit/`. Negative, mixed, or null results are reported verbatim alongside all mandatory `Q7` abstract disclosures.

---

## §2. Source Cohort, Locked `80.000 GB` Dataset Pool (`3 + 1`), Exact `P4` / NTC Rules, & Split-Half Partitions

### 2.1 Training Source Cohort (Fixed, Zero Leakage)

- **Primary Confirmatory Source Pool (`--source-mode r9_6`):** Strictly the **`6`** harmonized Round 9 human CRISPRi cell lines (`round9/p2/shared/lines/`, `6,148` aligned genes in `COMMON_GENES.txt`, `2,655` `union6` deduplicated target genes):
  1. `nadig_hepg2` (`2,343` targets)
  2. `nadig_jurkat` (`2,350` targets)
  3. `vcc2025_h1` (`300` targets)
  4. `replogle_k562_essential` (`2,049` targets)
  5. `replogle_rpe1` (`2,351` targets)
  6. `nourreddine_kolf21j` (`95` targets)
  Because all Round 11 confirmatory candidates and `feng_ipsc` (after `kolf_2`/`kolf_3` exclusion) are distinct cell lines/donors absent from $\mathcal{S} = \{1,\dots,6\}$, zero target-line perturbation data enter the primary predictor.
- **Secondary Source-Expansion Sensitivity (`SA-Src8`, `--source-mode r9_r10_8`, Strictly Non-Confirmatory):** Adds Round 10 `dld1` and `cd4t_rest` (`8` sources total) to quantify source-pool scaling. Results on `cd4t_stim48` under `SA-Src8` are explicitly flagged as **non-independent (`非独立灵敏度分析`)** because `cd4t_rest` shares cell type, sgRNA library, and lab with `cd4t_stim48`.

### 2.2 Locked Round 11 Dataset Pool (`80.000 GB` Cap, No Fallbacks — `Q3`, `Q4`)

All fallback options and open-ended placeholders are permanently removed (`Q4`). `X-Atlas/Orion` and `Cas9 KO` datasets are excluded. Total Round 11 download budget is capped at **`80.000 GB` (`80,000,000,000 B`)** (`Q3`). The locked cohort comprises **`3` confirmatory candidates + `1` secondary cohort** totaling **`76.254 GB` (`76,253,859,765 B`, `56` core files)** (`76.338 GB` including optional Feng Guide-UMI and cached Round 10 CD4i CSVs):

| Cohort Role | Dataset Key | Accession / Study | Cell Type / State / Lab | Core Download Size | Pre-Eval / Independence Gate |
| :--- | :--- | :--- | :--- | :---: | :--- |
| Confirmatory Candidate 1 | **`telohaec`** | `IGVFDS7943IIWZ` / `GSE210681` (`TeloHAEC_2285`, `PMC10921916` *Nat* 2024) | Resting unstimulated female TeloHAEC aortic endothelial (new lineage & lab) | **`1.487 GB`** (`3` files) | `-TSS*` & cis-element exclusion; `safe-targeting` strictly excluded (`P4`) |
| Confirmatory Candidate 2 | **`hcasmc`** | `IGVFDS6248NDUW` (Quertermous Lab, MOI `0.1–0.2`, ZIM3-KRAB) | `HCASMC-hTERT` coronary artery smooth muscle (new lineage & lab) | **`42.885 GB`** (`3` files) | **`Q3` Pre-Eval Feasibility Gate** on `36 GiB` host (`Parse v1.2.1` CSC HDF5 + guides); `safe-targeting` strictly excluded |
| Confirmatory Candidate 3 | **`cd4t_stim48`** | `GSE314342 Stim48hr` (Marson Lab, `CE0008678` & `CE0010866`) | Stimulated (`48h`) primary human `CD4+ αβ T` (`2` donors) | **`27.387 GB`** (`48` `.h5` + cached R10 CSVs) | **`Q6`** donor/timepoint disjoint from R10 `cd4t_rest`; **`Q2`/`Q7` weaker independence** (shared cell type, library, lab) |
| **Mandatory Secondary Only** | **`feng_ipsc`** | Figshare `27989294 v2` (`Feng et al. 2026 TargetedScreen`) | **`17` non-KOLF** HipSci human iPSC lines (`9` donors, `444` targets) | **`4.494 GB`** (`2` core files; `+73.78 MB` opt.) | **Never enters `K_conf`**; **`Q6`** exclude `kolf_2` & `kolf_3`; exclude `unassigned`; targeted gene panel; tests **`H-D5`** |
| **Locked Total** | **4 datasets** | — | — | **`76.254 GB`** (`56` core files; `76.338 GB` w/ opt.) | Within **`80.000 GB`** cap (`3.662 GB` headroom) |

Exact pre-locked files, single-sgRNA (`P4`), NTC-only, and disjoint split-half specifications for each dataset:

#### 2.2.1 Confirmatory Candidate 1 — `telohaec` (`TeloHAEC_2285`, `IGVFDS7943IIWZ` / `GSE210681`, `1.487 GB` Bundle / `1.473 GB` `.h5ad`)
- **Locked Files:** `IGVFFI8645YMGX.h5ad` (`1,473,190,952 B`), `IGVFFI9754AGFB.tsv.gz` (`729,687 B`, IGVF guide metadata), `GSE210681_Dialout_guides_per_CBC.txt.gz` (`12,857,614 B`, cell-barcode guide UMI table).
- **`P4` Single-sgRNA & NTC-Only Lock:** Assign single-sgRNA cells strictly where **`max_umi_ct >= 4` AND `X2nd_place_umi_ct < 4`**.
  - **Primary NTC ONLY:** Restrict NTC controls strictly to **`type == "non-targeting"`** (`400` guides, prefix `negative_control__`, Dialout `negative-control`, `2,656` pre-QC cells). **Strictly exclude `type == "safe-targeting"`** (`600` guides, prefix `safe_targeting__`, of which `538/600` carry `targeting == True` and `3,912` pre-QC cells — never pool `safe-targeting` into NTC!).
  - **Target Symbol Purification:** Retain **HGNC-resolvable single-gene plain symbols ONLY** (`2,251` plain symbols among `2,618` Dialout `P4` targeting symbols -> `2,133` directly HGNC-approved, `1,064` in `R9_COMMON`). Strictly exclude `199` alternative-promoter guides matching **`-TSS*`** (`192` collide with plain gene symbols and would otherwise create ambiguous promoter mixing) and exclude `168` non-coding/dual elements (`Enhancer-*`, `9p21-*`, `*-and-*`). *(An exploratory `-TSS*` pooled sensitivity check may be run post-eval under `SA-5`, never altering primary `telohaec` targets.)*
- **Split-Half Partition:** Sort the `20` 10x channel/lane IDs (`1..20` spanning `11AMDox_1..6`, `2PMDox1_1..7`, `2PMDox2_1..7`) lexicographically by numeric/canonical lane index and partition by parity: **`Half A` = `10` odd lanes (`1,364` NTC)** vs **`Half B` = `10` even lanes (`1,292` NTC)** (`20/20` lanes have `>= 31` pre-QC NTCs, exceeding the `>= 15` lane threshold).

#### 2.2.2 Confirmatory Candidate 2 — `hcasmc` (`HCASMC-hTERT`, `IGVFDS6248NDUW`, `42.885 GB`)
- **Locked Files:** `IGVFFI5999NXAF.h5` (`filtered_feature_barcode_matrix.hdf5`, `42,885,130,097 B`, Parse Biosciences `v1.2.1` CSC HDF5 layout; excludes the `87.29 GB` `IGVFFI9007DBZO.h5ad`) + `IGVFFI4575UMXX.tsv.gz` (`80,945 B`) + `IGVFFI9546UECF.tsv.gz` (`76,842 B`).
- **`P4` Single-sgRNA & NTC-Only Lock:** RNA QC (`MT <= 5–10%`, `ngenes >= 500`) plus Perturb-seq Poisson-Gaussian guide thresholding + dominance on guide-capture rows (`feature_type != "Gene Expression"`): retain cells with `n_pos == 1`, `t1 >= 5.0`, and (`t2 < 3.0` OR `t1 / max(t2, 0.5) >= 3.0`).
  - **CRITICAL NTC-Only Filter:** Metadata inspection shows BOTH `non-targeting` (`125` guides, `non-targeting_`) AND `safe-targeting` (`125` guides, `0Safe_`) carry `targeting == False`. Therefore NTC filtering MUST match strictly on **`type == "non-targeting"`** ONLY, **strictly excluding `safe-targeting`**.
- **Split-Half Partition (With Pre-Locked Missing-Unit Rule):** Extract sublibrary IDs (`3A..3H`, `4A..4F`, `14` Parse sublibraries) from `/matrix/barcodes` suffix `__<sublib>`: **`Half A = {3A, 3C, 3E, 3G, 4A, 4C, 4E}` (`7`)** vs **`Half B = {3B, 3D, 3F, 3H, 4B, 4D, 4F}` (`7`)** (using half-pooled NTC centering if any sublibrary has `< 15` NTC).
  - **Pre-Locked Missing-Unit Split Rule (`Q2`):** If `/matrix/barcodes` inside `IGVFFI5999NXAF.h5` (or any dataset lacking usable batch-unit metadata upon ingestion) contains `< 2` parseable sublibrary suffixes after prefix/suffix stripping, apply a deterministic stratified 50/50 split (`np.random.default_rng(42)`) across QC-valid cells per target/NTC and explicitly flag that condition in `R11_QC_COUNTS.md`.
- **`Q3` Mandatory Pre-Evaluation Technical Feasibility Gate:**
  Immediately after downloading `IGVFFI5999NXAF.h5` and **strictly BEFORE computing any evaluation metrics (`r_orth_cal`, `FCT`, nulls, or `CONFIRMATORY_LINES_R11.txt`)**, `A-T1` MUST verify two hard technical conditions on the `36 GiB` execution host:
  1. Streaming/chunked column reading of `IGVFFI5999NXAF.h5` (`Parse v1.2.1` CSC `/matrix/{data,indices,indptr,features,barcodes}`) completes end-to-end within host memory without OOM or HDF5 corruption; AND
  2. `feature_type != "Gene Expression"` rows exist inside `.h5` and yield non-empty extractable single-sgRNA guide assignments satisfying `P4`.
  If **either** condition fails, `hcasmc` is excluded **prior to evaluation** as a **feasibility exclusion (`可行性排除`)**, removed permanently from `K_conf`, and its exact technical failure trace/root cause is recorded in `R11_QC_COUNTS.md`. `hcasmc` may NEVER be excluded after computing or inspecting any evaluation outcome.

#### 2.2.3 Confirmatory Candidate 3 — `cd4t_stim48` (`GSE314342 Stim48hr`, `27.387 GB`)
- **Locked Files:** `48` GEO `.h5` lane matrices `GSM9393875..GSM9393898` (`CD4i_R2L25..R2L48`, `24` lanes per donor: `D3 = 13.591 GB`, `D2 = 13.797 GB`) + reuse cached `round10/shared/raw/cd4t_rest/*.csv` (`9.999 MB`).
- **`P4` Single-sgRNA & NTC-Only Lock:** Identical to Round 10 `cd4t_rest` (`MT <= 5.0%`, `n_genes >= 500`, per-donor Poisson-Gaussian `min_thr` over `[5..20]`, `top1 / top2 >= 3.0`, `CD4_EXCEL_GUIDE_REMAP`). Primary NTC ONLY: `sgRNA.startswith("NTC-")` (`992` NTC guides, `0` safe-targeting; strictly exclude `315` blank/mismatch guides).
- **Split-Half Partition:** Exact lexicographic disjoint human donor split: **`Half A = CE0008678`** (`D3_Stim48hr`, `24` lanes, `77,334` pre-QC NTCs, `1,808,899` targeting) vs **`Half B = CE0010866`** (`D2_Stim48hr`, `24` lanes, `87,747` pre-QC NTCs, `1,965,311` targeting); `48/48` lanes have `>= 3,027` NTCs (`>= 15`). Inter-donor biological variation enters the split-half noise ceiling identically to Round 10 `cd4t_rest`.
- **`Q2` + `Q6` Pre-QC Disjointness & Weakened Independence Disclosure:**
  - Verified **pre-download (`Q6`)** that `CE0008678` (`D3_Stim48hr`) and `CE0010866` (`D2_Stim48hr`) are **100% disjoint** in human donor ID, 10x processing lanes (`R2L25..R2L48`), and stimulation state (`Stim48hr`) from Round 10 `cd4t_rest` (`CE0006864` / `CE0008162`, `Rest`).
  - **Explicit Independence Limitation (`Q2`, `Q7`):** Unlike `telohaec` and `hcasmc` (completely unused lineages & labs), `cd4t_stim48` shares primary `CD4+ αβ T` cell type, sgRNA library, and Marson lab with Round 10 `cd4t_rest` (**weaker biological/laboratory independence**). Although `--source-mode r9_6` does not train on `cd4t_rest`, if `cd4t_stim48` enters `K_conf`, this shared cell-type/library/lab limitation MUST be disclosed verbatim in the Abstract (`Q7`).

#### 2.2.4 Mandatory Secondary Cohort (Never Enters `K_conf`) — `feng_ipsc` (`4.494 GB` Core / `4.568 GB` w/ Opt.)
- **Locked Files:** Figshare `27989294 v2`: `TargetedScreen_RNA-UMI-Counts.csv.gz` (`51072902`, `4,480,969,391 B` — must be streamed row-by-row on `36 GiB` host) + `TargetedScreen_Cell-Metadata.tsv.gz` (`51051809`, `13,491,635 B`) (+ optional `TargetedScreen_Guide-UMI-Counts.csv.gz` `51051815`, `73.78 MB`). (`Fitness` & `NonFitness` screens excluded due to insufficient NTC counts `32`/`10`.)
- **Exact `kolf_2` & `kolf_3` Exclusion (`Q6`) & Line Scope:** Whereas earlier 3-screen screening (`R11_DATA_SCREEN.md`) counted `43` non-KOLF lines across all screens combined, metadata inspection (`R11_DOWNLOAD_PLAN.md`) confirms `TargetedScreen` (`Pica`) contains **`19` HipSci iPSC lines (`10` donors)** total (`1,161,864` cells); **`kolf_2` (`64,009` cells) and `kolf_3` (`69,622` cells) MUST be strictly excluded pre-QC** (`KOLF2.1J` in Round 9 source pool originates from `kolf_2`), leaving **`17` non-KOLF HipSci lines (`9` donors)** covering `444` target genes (`333` in `R9_COMMON`).
- **`P4` Single-sgRNA & NTC-Only Lock:** Retain single-sgRNA cells called in `Guide_Call` (`freq_max_guide > th_freq` & `max_UMI >= 3`). **CRITICAL:** Strictly exclude `470,561` `"unassigned"` non-KOLF cells (**never pool `"unassigned"` with NTC controls!**), retaining strictly **`7,266` `NonTarget_*`** (`20` guides, `0` safe-targeting) + **`550,406` single-guide targeting cells**.
- **Split-Half Partition:** Sort the `17` non-KOLF `Cell_Line` IDs lexicographically and split by `1`-based index parity:
  - **`Half A` (`9` lines, `9` donors):** `{eipl_1, fiaj_1, iudw_1, jejf_2, oikd_2, paab_3, pipw_4, tolg_4, zapk_3}` (`4,091` NTC).
  - **`Half B` (`8` lines, `8` donors):** `{eipl_3, fiaj_3, iudw_4, jejf_3, oikd_5, paab_4, pipw_5, tolg_6}` (`3,175` NTC).
  `17/17` cell lines have `>= 86` NTC (`>= 15`), and `444/444` targets have `>= 25` total cells and `>= 10` per half. Within each half, center perturbed cells by `Cell_Line x Inlet` when that unit has `>= 15` NTC cells, otherwise by line/half-pooled NTC mean (`§2.3`).
- **Role Restriction & Descriptive `H-D5`:** `feng_ipsc` is permanently **secondary-only** (never enters `K_conf`) due both to sharing the iPSC lineage class with Round 9 `nourreddine_kolf21j`/`vcc2025_h1` and to `TargetedScreen` using a targeted gene panel below the `>= 3,000` whole-transcriptome confirmatory floor (secondary evaluation runs on the full finite overlap $G_{\mathrm{eval}}(\mathrm{feng\_ipsc}) = G_{\mathrm{R9\_6148}} \cap G_{\mathrm{feng\_ipsc}}^{\mathrm{fin}}$ without requiring `>= 3,000`). On `feng_ipsc`, evaluate descriptive question **`H-D5`** (`§5`): compare `FCT` (using `feng_ipsc` non-KOLF $\mu_t$) against `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)`.

### 2.3 Standardized Normalization, Gene Alignment, & Addendum-B Target Eligibility

- **Normalization & Expression:** Identical to Round 9/10 (`round9/p2/shared/scripts/prepare_lines.py` & `PREREG_C §2.4`):
  1. Library-size normalize raw UMI counts to $10^4$ per cell (`X_libnorm`).
  2. Transform via $\operatorname{log1p}(X_{\mathrm{libnorm}})$; subtract unit/batch-matched (`channel`/`lane`/`sublibrary`/`Cell_Line`) NTC `log1p` mean $\mu_{t,\mathrm{batch}}$ per perturbed cell prior to averaging per-target pseudobulk $\Delta^{\mathrm{obs}}_{t,p}$.
  3. Store batch-stratified unperturbed control libnorm mean $\mu_{t,g}$ (`20,000`-cell stratified NT subsample + full-NT check `SA-5b`, using full batch-matched NT means for $\Delta^{\mathrm{obs}}$ centering and `FCT`) and `C_t` covariance for Round 9 anchor comparators.
- **Gene Alignment & Target Eligibility (`P5`, `only_reproducible = True`):** Identical to `PREREG_C §2.2 & §2.4`:
  - Map gene symbols to canonical HGNC (`Ensembl` cross-checked), drop non-finite/duplicate/non-constant-panel symbols, and form $G_{\mathrm{eval}}(t) = G_{\mathrm{R9\_6148}} \cap G_t^{\mathrm{fin}}$. Confirmatory lines require `|G_eval(t)| >= 3,000`.
  - Evaluated targets $p$ must satisfy all four conditions: (1) reproducible on target line $t$ under Addendum B (`n_cells >= 25` single-sgRNA cells, `>= 50` top-5,000 genes passing two-sample Anderson–Darling `BH FDR < 0.05` vs `1,000` stratified NT, `20-PC` energy-distance permutation `p < 0.001` over `2,000` perms vs `500` NT; 1 deduplicated row per symbol); (2) reproducible on `>= 1` source line in the Round 9 source pool ($|\mathcal{S}_p^{\mathrm{rep}}| \ge 1$); (3) target self-gene $g_p \in G_{\mathrm{eval}}(t)$ (self-gene $g_p$ excluded from evaluation masks); (4) split-half valid (`>= 2` cells in Half A AND `>= 2` in Half B, with `>= 15` NTC per contributing batch unit or half-pooled NT backup rule when unit `< 15`).
- **Mechanical `K_conf` Assignment & `P1` Abortion Guard:**
  - Confirmatory candidates `{telohaec, hcasmc [if feasible], cd4t_stim48}` passing `|G_eval(t)| >= 3,000` and **`n_eval_targets >= 50`** enter `K_conf` mechanically via `CONFIRMATORY_LINES_R11.txt`.
  - **If `K_conf >= 2`:** All confirmatory claims require **unanimous passage (`K_conf / K_conf`, `100%`)** (`3/3` or `2/2`). If `K_conf = 2`, Abstracts MUST state `"K_conf = 2 confirmatory cell lines"`.
  - **If `K_conf < 2` (`P1`):** Round 11 **performs NO confirmatory hypothesis adjudication (`不做确认性判定`)**, reports strictly **descriptive results only**, and states verbatim in the Abstract: *"Insufficient confirmatory cell lines (`K_conf < 2`); confirmatory hypothesis testing was not performed (`确认性线不足，未做确认性检验`)."*

---

## §3. Predictors & Null Baselines (`FCT`, `E1_ConstMed`, `NegCtrl`, `Perm100`)

1. **Primary Confirmatory Predictor (`FCT`, Unchanged from `PREREG_A v4` / `PREREG_C §3.1`):**
   $$\mu_{s,g} = \max(\mu^{\mathrm{libnorm}}_{s,g},\, 10^{-8}), \qquad \hat{p}_{s,g} = \max\!\big(\operatorname{expm1}(\operatorname{log1p}(\mu_{s,g}) + \Delta^{\mathrm{obs}}_{s,p,g}),\, 10^{-8}\big),$$
   $$\overline{\mathrm{FC}}_{p,g} = \exp\!\left(\frac{1}{|\mathcal{S}_p|} \sum_{s \in \mathcal{S}_p} \ln\!\left(\frac{\hat{p}_{s,g}}{\mu_{s,g}}\right)\right), \qquad \hat{\Delta}^{\mathrm{FCT}}_{t,p,g} = \operatorname{log1p}\!\big(\mu_{t,g}\, \overline{\mathrm{FC}}_{p,g}\big) - \operatorname{log1p}(\mu_{t,g}).$$
2. **Pre-Registered Target-Baseline-Free Constant Null (`E1_ConstMed`, for `H-D3` & `H-D4`):**
   Compute the single gene-invariant median control expression across $G_{\mathrm{eval}}(t)$, $c_{\mathrm{med}} = \operatorname{median}_{g \in G_{\mathrm{eval}}(t)} \max(\mu^{\mathrm{libnorm}}_{t,g}, 0)$, and pass constant vector $\mu^{\mathrm{E1}}_{t,g} \equiv c_{\mathrm{med}}$ through identical `FCT` nonlinear `log1p` compression:
   $$\hat{\Delta}^{\mathrm{E1}}_{t,p,g} = \operatorname{log1p}\!\big(c_{\mathrm{med}}\, \overline{\mathrm{FC}}_{p,g}\big) - \operatorname{log1p}(c_{\mathrm{med}}).$$
3. **`NegCtrl_FCT_MismatchedMean(s')` (`6` Primary Mismatched Controls, `s' = 1..6`):**
   Replace $\mu_t$ inside `FCT` with source line $s'$'s libnorm control mean $\mu_{s'}$ for each of the `6` Round 9 source lines (`nadig_hepg2`, `nadig_jurkat`, `vcc2025_h1`, `replogle_k562_essential`, `replogle_rpe1`, `nourreddine_kolf21j`):
   $$\hat{\Delta}^{\mathrm{FCT\_Mis}(s')}_{t,p,g} = \operatorname{log1p}(\mu_{s',g}\overline{\mathrm{FC}}_{p,g}) - \operatorname{log1p}(\mu_{s',g}).$$
4. **`NegCtrl_FCT_Perm(k)` (`100` Gene-Label Permutations, `k = 0..99`):**
   Permute target control mean $\mu_t$ across $G_{\mathrm{eval}}(t)$ using `np.random.default_rng(k).permutation(G)` while holding $\overline{\mathrm{FC}}_p$ fixed, evaluated identically under BOTH **`all_nontarget`** (confirmatory `H-D2`) and **`resp_top50`** (`SA-7` / `H-D4` context).
5. **Primary Anchor (`Transfer_Mean`):** Unscaled cross-source average $\hat{\Delta}^{\mathrm{TMean}}_p = \frac{1}{|\mathcal{S}_p|}\sum_{s \in \mathcal{S}_p} \Delta^{\mathrm{obs}}_{s,p}$, defining residuals $r^{\mathrm{obs}}_p = \Delta^{\mathrm{obs}}_{t,p} - \hat{\Delta}^{\mathrm{TMean}}_p$ and $r^{\mathrm{pred}}_p = \hat{\Delta}_{t,p} - \hat{\Delta}^{\mathrm{TMean}}_p$, plus `H-D_util` non-inferiority.
6. **Comparator & Exploratory Diagnostics (`SA-2`, `SA-3`, `SA-5`, Strictly Non-Confirmatory):** `Transfer_Mean_GlobalScale`, `Hybrid_ProjShared` (frozen 3D subspace $\operatorname{span}\{\mathbf{1}, \bar{b}_S, a_S\}$, `commit 1ca9f60f`), `Transfer_Mean_VarScaled` (`TMVS`), `Hybrid_on_TMVS`, `CIPHER_Linear`, `Source_Grand_Mean`, plus constant/scaling ablations `E2_ConstP10`, `E2_ConstP90`, `E3_SrcMeanCtrl`, `E5_RankMatchedMu`, `E6_LogFCBar`, and `E7_FCMinusOne`.

---

## §4. Endpoints, Multiplicity Control (`97.5%` CIs via IUT + 2-Claim Bonferroni), & Adjudication

### 4.1 Exact Multiplicity Justification (`Q1`, `Q5`, `Q8`) & De-Amplituded Primary Endpoint (`r_orth_cal`)

**Exact Familywise Error Control & `97.5%` Confidence Interval Rationale (`Q1`, `Q5`, `Q8`):**
- Round 11 tests **`2` separate confirmatory claims**:
  1. **Primary Biological Claim (`H-D1 ∩ H-D2 ∩ H-D3`):** Because all three conditions (`H-D1`, `H-D2`, and `H-D3`) are required simultaneously (`conjunctive / AND`), Claim 1 is an **Intersection–Union Test (IUT)**. By standard IUT theory, testing each component hypothesis individually at level $\alpha_{\mathrm{claim}}$ guarantees overall Type-I error $\le \alpha_{\mathrm{claim}}$ without dividing $\alpha$ across the three conjunctive conditions.
  2. **Independent Methodological Claim (`H-D4 = H-D4a ∩ H-D4b`, `Q8`):** Also tested via IUT (`H-D4a` AND `H-D4b`), adjudicated independently of Claim 1 across all `K_conf` lines.
- Applying Bonferroni correction across the **`2` confirmatory claims** yields per-claim $\alpha_{\mathrm{claim}} = 0.05 / 2 = 0.025$, i.e., **two-sided `97.5%` percentile bootstrap confidence intervals (`[1.25%, 98.75%]`)** matching Round 9 and `PREREG_C`.
- Per `Q5`, all `6` mismatched control upper bounds $\max_{s'=1..6} \mathrm{CI}_{97.5,\mathrm{hi}}(\mathrm{NegCtrl\_FCT\_MismatchedMean}(s'))$ inside `H-D2` likewise use **`97.5%` (`[1.25%, 98.75%]`)** (while `Perm100` uses empirical $Q_{0.975}$ across `k = 0..99`). Parallel utility non-inferiority (`H-D_util`, `H-C3` continuation) likewise uses two-sided **`97.5%` CIs (`[1.25%, 98.75%]`)**.

**De-Amplituded Endpoint (`r_orth_cal`) & Primary vs. Secondary Evaluation Masks:**

For each evaluated perturbation $p$ and evaluation gene mask $M_p \subseteq G_{\mathrm{eval}}(t) \setminus \{g_p\}$:
- **Primary Evaluation Mask (`all_nontarget`, Confirmatory `H-D1`, `H-D2`, `H-D3`, & `H-D5`):**
  $$M_p^{\mathrm{all}} = G_{\mathrm{eval}}(t) \setminus \{g_p\}$$
  All finite non-self genes in the common aligned panel are retained without thresholding on $|\Delta^{\mathrm{obs}}_{t,p,g}|$, completely eliminating observation-driven gene selection bias.
- **Secondary Evaluation Mask (`resp_top50`, Confirmatory Trap Audit `H-D4` & `SA-7`):**
  $$M_p^{\mathrm{top50}} = \operatorname{Top50}_{g \in G_{\mathrm{eval}}(t) \setminus \{g_p\}} |\Delta^{\mathrm{obs}}_{t,p,g}|$$
  Retains Round 9/10's top 50 strongest observed responses per perturbation strictly to test `H-D4` and secondary comparability.

Within mask $M_p$, both observed and predicted residuals relative to `Transfer_Mean` are orthogonally projected perpendicular to $\hat{\Delta}^{\mathrm{TMean}}_p$:
$$r^{\mathrm{obs}}_{\perp, p} = r^{\mathrm{obs}}_p - \frac{\langle r^{\mathrm{obs}}_p, \hat{\Delta}^{\mathrm{TMean}}_p\rangle}{\|\hat{\Delta}^{\mathrm{TMean}}_p\|_2^2}\hat{\Delta}^{\mathrm{TMean}}_p, \qquad r^{\mathrm{pred}}_{\perp, p} = r^{\mathrm{pred}}_p - \frac{\langle r^{\mathrm{pred}}_p, \hat{\Delta}^{\mathrm{TMean}}_p\rangle}{\|\hat{\Delta}^{\mathrm{TMean}}_p\|_2^2}\hat{\Delta}^{\mathrm{TMean}}_p.$$
Pooling orthogonalized entries across all $N_{\mathrm{eval}}$ perturbations yields raw Pearson correlation $r^{\mathrm{raw}}_{\perp}$ and Spearman–Brown split-half ceiling $c^{\perp}_{\mathrm{ceil}} = \sqrt{2\rho^{\perp}_{AB}/(1+\rho^{\perp}_{AB})}$ computed on the `§2.2` disjoint split halves under the identical mask $M_p$, giving:
$$r^{\mathrm{orth\_cal}} = \operatorname{clip}\!\left(\frac{r^{\mathrm{raw}}_{\perp}}{c^{\perp}_{\mathrm{ceil}}},\, -1,\, 1\right) \quad (\text{whenever } c^{\perp}_{\mathrm{ceil}} > 0).$$
All intervals use **joint perturbation-level bootstrap (`n_boot = 4,000`, `np.random.default_rng(20261123)`)**, recomputing both $r^{\mathrm{raw}}_{\perp}$ and $c^{\perp}_{\mathrm{ceil}}$ inside every bootstrap draw. Paired contrasts (`FCT - E1_ConstMed` in `H-D3`, and cross-mask `E1(resp_top50) - E1(all_nontarget)` in `H-D4b`) are evaluated across **identical perturbation bootstrap draws**, incorporating each mask's split-half noise ceiling uncertainty jointly.

### 4.2 Confirmatory Decision Rules (`K_conf >= 2`, Unanimous `100%` Passage)

- **`P1` (`K_conf < 2`) Guard:** No confirmatory claims tested (`不做确认性判定`); descriptive reporting only.
- **Confirmatory Claim 1 — Primary Biological Claim (`FCT` under `all_nontarget`):**
  Declared **CONFIRMED (`确认`)** if and only if **all three conditions** hold on **`K_conf / K_conf` (`100%`)** confirmatory cell lines; otherwise declared **`NOT CONFIRMED (未确认)`**:
  - **`H-D1` (Whole-Transcriptome De-Amplituded Signal):** `FCT` `r_orth_cal` (`all_nontarget`) **`97.5% CI_lo > 0`**.
  - **`H-D2` (Target-Baseline Specificity Over Mismatches & Permutations):** `FCT` `r_orth_cal` (`all_nontarget`) **`97.5% CI_lo`** strictly exceeds BOTH:
    1. $\max_{s'=1..6} \mathrm{CI}_{97.5,\mathrm{hi}}(\mathrm{NegCtrl\_FCT\_MismatchedMean}(s'))$ under `all_nontarget`, AND
    2. $Q_{0.975}(\mathrm{NegCtrl\_FCT\_Perm}(0..99))$ under `all_nontarget`.
  - **`H-D3` (Superiority Over Target-Agnostic Constant Scaling `E1_ConstMed`):** Joint paired difference `(FCT - E1_ConstMed)` on `r_orth_cal` (`all_nontarget`) has **`97.5% CI_lo > 0`**.
- **Confirmatory Claim 2 — Independent Methodological Claim (`H-D4`: `resp_top50` Selection Trap, `Q8`, `Q9`):**
  Declared **CONFIRMED (`确认`)** independently if and only if **both conditions** hold on **`K_conf / K_conf` (`100%`)** confirmatory cell lines; otherwise declared **`NOT CONFIRMED (未确认)`**:
  - **`H-D4a`:** Under `resp_top50`, the target-baseline-free constant predictor `E1_ConstMed` exhibits positive `r_orth_cal` with **`97.5% CI_lo > 0`**.
  - **`H-D4b`:** Joint cross-mask paired difference `E1_ConstMed(resp_top50) - E1_ConstMed(all_nontarget)` on `r_orth_cal` has **`97.5% CI_lo > 0`**.
  *(Consistent with Nicol, Shivakumar & Irizarry 2026 [`10.64898/2026.05.07.723486`], `H-D4` isolates `resp_top50`'s artificial elevation of deterministic target-baseline-free `E1_ConstMed` scaling relative to `all_nontarget`.)*
- **Parallel Practical Utility Verdict (`H-D_util`, Non-Inferiority Relative to `Transfer_Mean`):**
  Co-featured alongside Claims 1 and 2 across all **`6` aggregate conventions** (`ctrl_disc_panel`, `ctrl_disc_ownmasked`, `ctrl_de50_r`, `cent_disc_panel`, `cent_disc_ownmasked`, `cent_de50_r`):
  - **Non-inferior (`H-D_util PASS`):** `diff_ci975_hi >= 0` across all `6` metrics on `K_conf / K_conf` lines (`0` inferior metrics). Also record strict superiority counts (`diff_ci975_lo > 0`) and descriptive `pds_cosine_panel`.

---

## §5. Pre-Specified Secondary & Sensitivity Analyses (`SA-1`–`SA-8`, `H-D5`)

All analyses below are strictly **non-confirmatory** (`次要/敏感性分析`):

- **`SA-1` (Amplitude Slope Check):** Report `FCT` `slope_ols` (`[97.5% CI]`), `sd_pred / sd_obs`, and `slope_deming` vs `[0.02, 0.13]` under both `all_nontarget` and `resp_top50`.
- **`SA-2` (Raw `r_cal` & `GlobalScale` Check):** Report un-orthogonalized `r_cal` and `Transfer_Mean_GlobalScale` alongside `r_orth_cal`.
- **`SA-3` (Round 9 Anchors & `SA-Src8` 8-Source Sensitivity):** Report `Hybrid_ProjShared` (3D frozen source subspace `commit 1ca9f60f`), `TMVS`, `Hybrid_on_TMVS`, `CIPHER_Linear`, and `Source_Grand_Mean` under `--source-mode r9_6`, plus **`SA-Src8`** (`--source-mode r9_r10_8`, adding R10 `dld1` & `cd4t_rest`, flagged non-independent on `cd4t_stim48`).
- **`SA-4` (`R19` 3-Tier Variance Decomposition):** Partition $\mathrm{SST}(r^{\mathrm{obs}})$ into whole-line (`L1`), per-gene (`L2`, genes observed in $\ge 5$ perturbations), and residual (`L3`) on both `all_nontarget` and `resp_top50`.
- **`SA-5` (Universe & Control-Subsampling Sensitivity):** `S5` top-75% expressed/variable genes, `S1` all-targets (`--all-targets`), `SA-5b` full-NT vs `20k`-stratified NT $\mu_t$ parity, exploratory `E2_ConstP10`/`E2_ConstP90`/`E3_SrcMeanCtrl`/`E5_RankMatchedMu`/`E6_LogFCBar`/`E7_FCMinusOne` across both masks, plus `-TSS*` exclusion sensitivity on `telohaec`.
- **`SA-6 / S7` (`arXiv:2608.00152v3` Cell-Count & Knockdown-Depth Controls):** Run identically under BOTH `all_nontarget` and `resp_top50`: per-perturbation `r_orth_p`, Spearman rank correlation with $\ln(n_{\mathrm{cells},p})$ and $|k_p|$ (`97.5%` CIs), mean-centered OLS intercept on `[1, ln(n_cells) - mean, |k| - mean]`, and `n_cells` tertiles (`T1`, `T2`, `T3`).
- **`SA-7` (`resp_top50` Parallel Audit & `Perm100` Inflation):** Report `FCT`, `E1_ConstMed`, `6` mismatched controls, and `Perm100` (`median`, `Q_{0.975}`) under `resp_top50` side-by-side with `all_nontarget`.
- **`SA-8` (Split-Half Selected Top-50 Validation `resp_top50_split`, `Q9`):**
  To test directly whether `resp_top50` score inflation disappears when top-50 gene selection is out-of-sample relative to the evaluated observed response:
  - **Direction 1 (`A -> B`):** For each perturbation $p$, rank $G_{\mathrm{eval}}(t) \setminus \{g_p\}$ using **`half_a_obs` ONLY** (`|\Delta^{\mathrm{obs}}_{t,\mathrm{HalfA},p,g}|`) to select mask $M_p^{(A)}$ (`top 50`). Orthogonally project `pred` residual ($\hat{\Delta}^{\mathrm{pred}}_{t,p} - \hat{\Delta}^{\mathrm{TMean}}_p$) and **`half_b_obs` residual** ($\Delta^{\mathrm{obs}}_{t,\mathrm{HalfB},p} - \hat{\Delta}^{\mathrm{TMean}}_p$) perpendicular to $\hat{\Delta}^{\mathrm{TMean}}_p$ on $M_p^{(A)}$, yielding pooled Pearson correlation $r^{\mathrm{raw}}_{\perp, A\to B}$ and single-half split-ceiling disattenuated correlation:
    $$r^{\mathrm{orth\_cal}}_{A\to B} = \operatorname{clip}\!\left(\frac{r^{\mathrm{raw}}_{\perp, A\to B}}{\sqrt{\rho^{\perp}_{AB}(M^{(A)})}},\, -1,\, 1\right).$$
  - **Direction 2 (`B -> A`):** Symmetric reverse: rank top 50 on **`half_b_obs`** ($M_p^{(B)}$), evaluate `pred` residual against **`half_a_obs`** residual, disattenuated by $\sqrt{\rho^{\perp}_{AB}(M^{(B)})}$.
  - **Symmetric Endpoint:** Report averaged `r_orth_raw_split = (r_raw_{A->B} + r_raw_{B->A}) / 2` and `r_orth_cal_split = (r_cal_{A->B} + r_cal_{B->A}) / 2` with `97.5%` joint perturbation bootstrap CIs (`n_boot = 4,000`, `seed = 20261123`) for **`FCT`**, **`E1_ConstMed`**, and **`Perm100` (`median` & `Q_{0.975}`)**.
- **`H-D5` (`feng_ipsc` Specific Non-Confirmatory Secondary Test):**
  Under `all_nontarget` on `feng_ipsc` (`17` non-KOLF lines, `9` donors), compare **`FCT`** (using `feng_ipsc` own non-KOLF $\mu_t$) vs **`NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)`** (`KOLF2.1J` control mean, same human iPSC cell type, `[97.5% CI]`) to distinguish **specific cell-line / genetic-donor baseline fidelity (`FCT > Mis(kolf21j)`)** from **broad cell-type baseline fidelity (`Mis(kolf21j) ≈ FCT`)**.

---

## §6. Prospective Statistical Power Analysis (`97.5%` CIs, $z_{97.5} = 2.2414027$, $z_{0.80} = 0.8416212$)

Using observed Round 10-E `all_nontarget` asymptotic SE scaling ($\widehat{\mathrm{SE}}(N) \approx c/\sqrt{N}$, matching joint perturbation bootstrap `97.5%` CI half-widths at $z_{97.5} = 2.2414027$ and 80% one-sided detection shift $z_{0.80} = 0.8416212$):
$$N_{50} = \left(\frac{z_{97.5}\, c}{\Delta_{\mathrm{effect}}}\right)^2, \qquad N_{80} = \left(\frac{(z_{97.5} + z_{0.80})\, c}{\Delta_{\mathrm{effect}}}\right)^2,$$
with mismatched/difference SEs calibrated under conservative independent variance ($\rho = 0$, $c_{\mathrm{diff}} = \sqrt{c_F^2 + c_M^2}$) and moderately correlated bootstrap variance ($\rho = 0.5$, $c_{\mathrm{diff}} = \sqrt{c_F^2 + c_M^2 - c_F c_M}$):

| Benchmark Reference (`R10-E` `all_nontarget`) | Hypothesis | Effect $\Delta_{\mathrm{effect}}$ | SE Scale $c$ | $N_{50}$ (Targets) | $N_{80}$ (Targets) |
| :--- | :--- | :---: | :---: | :---: | :---: |
| **`dld1`-like Clear Margin** (`r_FCT = +0.1355`) | **`H-D1`** (`FCT > 0`) | `+0.1355` | `c_F = 0.1175` | **`3.8`** | **`7.1`** |
| (`FCT - E1 = +0.1969`) | **`H-D3`** (`FCT > E1`) | `+0.196871` | `c_Delta = 0.1201` | **`1.9`** | **`3.5`** |
| (`FCT - max Mis/Perm = +0.080210`) | **`H-D2`** (`FCT > max Null`) | `+0.080210` | `c_F = 0.1175`, `c_M = 0.0983` | **`36.4`** | **`51.5`** ($\rho=0.5$) / **`58.3`** ($\rho=0$) |
| **`cd4t_rest`-like Narrow Margin** (`r_FCT = +0.0556`) | **`H-D1`** (`FCT > 0`) | `+0.0556` | `c_F = 0.1413` | **`32.5`** | **`61.5`** |
| (`FCT - E1 = +0.0426`) | **`H-D3`** (`FCT > E1`) | `+0.042610` | `c_Delta = 0.1385` | **`53.1`** | **`100.4`** |
| (`FCT - max Mis/Perm = +0.034728`) | **`H-D2`** (`FCT > max Null`) | `+0.034728` | `c_F = 0.1413`, `c_M = 0.1482` | **`349.1`** | **`492.7`** ($\rho=0.5$) / **`559.2`** ($\rho=0$) |

**Practical Interpretation:**
- On lines with a clear baseline margin (`dld1`-like), the confirmatory gate `n_eval_targets >= 50` provides $>99\%$ power for `H-D1` and `H-D3` and $\approx 80\%$ power for `H-D2` even against two separate $97.5\%$ bounds.
- On primary-donor lines with inter-donor noise ceilings (`cd4t`-like), `>= 50` targets is well-powered for `H-D1` at `n >= 62` (`302` shared targets on R10 `cd4t_rest`), borderline for `H-D3` (`N_80 = 100.4`), and **`H-D2` requires `~350–560` targets** to exceed mismatched controls under narrow margins. Unfavorable `H-D2` power outcomes resulting from genuinely small margins (`0.035`) are expected biological signal limits and must be reported transparently rather than post-hoc adjusted.

---

## §7. Software Governance, Code Freeze (`C1–C3`), & Mandatory Disclosures (`Q7`)

1. **Strict Blinding & Frozen Dependencies (`C1`):** `round11/src/` (`A-T2`) and `round11/audit/` (`A-T3`) import Round 10's frozen primitives without touching any file in `round10/` (`FREEZE_MANIFEST_R10.txt` `22/22 OK`). Only Round 9/10 unblinded datasets and synthetic fixtures are permitted before `CODE_FREEZE_R11.txt` and `R11_GO.txt`.
2. **Deterministic RNGs (`C2`):** `np.random.default_rng(20261123)` (`n_boot = 4,000`), `np.random.default_rng(k)` (`k = 0..99`) for gene permutations, and `seed = 42` for deterministic subsampling / missing-unit backup splits.
3. **Dual-Team `C3` Cross-Check Tolerance:** `max_point_abs_diff <= 1e-4` on deterministic metrics (`<= 1e-10` target on primary `FCT`), `max_ci_abs_diff <= 0.01`, and `100%` boolean agreement across all `H-D1`, `H-D2`, `H-D3`, `H-D4`, `H-D_util` verdicts.
4. **Mandatory Abstract Disclosures (`Q7` — Non-Negotiable):** Every Round 11 summary and manuscript Abstract MUST explicitly state ALL four items:
   1. Individual line-by-line verdicts for **`H-D1`**, **`H-D2`**, and **`H-D3`** across all confirmatory lines (and the overall `K_conf / K_conf` Primary Claim verdict);
   2. Separate unanimous verdict for **`H-D4`** (`H-D4a` & `H-D4b`) alongside `SA-8` (`resp_top50_split`) interpretation and **Nicol, Shivakumar & Irizarry (2026)** attribution (`Q9`);
   3. Exact final count **`K_conf`** (together with any pre-eval `hcasmc` feasibility exclusion per `Q3`, or the verbatim `P1` notice if `K_conf < 2`, or `"K_conf = 2 confirmatory cell lines"` if `K_conf = 2`);
   4. Explicit independence limitation for **`cd4t_stim48`** (shared primary `CD4+ αβ T` cell type, sgRNA library, and Marson lab with Round 10 `cd4t_rest`) whenever `cd4t_stim48` enters `K_conf`, alongside exact whole-transcriptome variance-explained range (`Q10`).
