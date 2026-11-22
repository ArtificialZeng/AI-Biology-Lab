# Round 12 Prospective Validation Pre-Registration Protocol (`PREREG_E.md`, `v2`)

- **Project:** AI-Biology-Lab Round 12 — Prospective Confirmation of Paired Cell-Line Baseline Specificity (`FCT - Mis_s`) and Split-Half Selection Inflation (`SA-8`) Across Perturbation Modalities (`CRISPRko` & `CRISPRi`)
- **Protocol Author:** Agent `A-T3` (Team A · Pre-Registration & Statistical Audit)
- **Authorizing Directives & Inputs:** Master Controller Directive [`R12_第十二轮行动指令.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/controller/directives/R12_第十二轮行动指令.md), Master Controller Review [`PREREG_E_AUDIT_R1.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/status/PREREG_E_AUDIT_R1.md) (addressing `E1–E6`), Locked Download & Ingestion Manifest [`R12_DOWNLOAD_PLAN.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/shared/R12_DOWNLOAD_PLAN.md), Round 12 Metadata Screening [`R12_DATA_SCREEN.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round11/shared/R12_DATA_SCREEN.md), Round 11 Independent Exploratory Audit [`REPORT.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round11/audit/exploratory_check/REPORT.md) + [`paired_rorth_12lines.csv`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round11/audit/exploratory_check/outputs/paired_rorth_12lines.csv), Prior-Art Audit [`PRIOR_ART_E.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round11/shared/PRIOR_ART_E.md), and Protocol Structural Template [`PREREG_D.md` (`v2`)](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round11/prereg/PREREG_D.md).

---

## §0. Scientific Scope, Key Round 12 Design Changes, & Blinding Guarantee

### 0.1 Why Round 12 Is Conducted (Post-Hoc Discovery Provenance & Exact Prior-Round Metrics)

**Post-Hoc Hypothesis Provenance (`E6`):** Both Claim 1 (`H-E1`, paired `FCT - Mis_s`) and Claim 2 (`H-E2`, `SA-8` `resp_top50_split`) hypotheses were formulated **after inspecting Rounds 9–11 data** (`round11/status/EXPLORATORY_REVIEW_R11.md` §D; [`REPORT.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round11/audit/exploratory_check/REPORT.md)) — neither can be treated as prospectively confirmed on Rounds 9–11 data, which is why blinded prospective validation on completely unseen Round 12 data is required:

1. **Paired Target-Baseline Specificity (`FCT - NegCtrl_FCT_MismatchedMean(s')`, Separating `r_orth_cal` From Post-Hoc `N_5`):**
   - Under `PREREG_D`'s *unpaired* interval-separation rule (`H-D2`: `FCT CI97.5_lo > max(Mis_s CI97.5_hi, Perm Q97.5)`), only `7/12` R9–R11 lines passed under BOTH `r_orth_cal` and `N_5` `r_cone_sp` (`cd4t_rest`, `cd4t_stim48`, `calu3`, `telohaec`, and `nourreddine_kolf21j` failed).
   - **Clarification of Prior Metrics (`E1`):** In [`REPORT.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round11/audit/exploratory_check/REPORT.md) Check 2B, `11/12` passing paired `FCT - Mis_s` was computed under the post-hoc **`N_5` semi-partial metric `r_cone_sp`** (`nourreddine_kolf21j` failed `1/5` contrast `k562`). Because `N_5` was constructed after inspecting R9–R11 data, **`H-E1b` strictly remains on `r_orth_cal` (`all_nontarget`)** — matching Rounds 9–11 — and `N_5` paired contrasts are placed strictly in `§5` descriptive analyses (`只报告、不判定`).
   - **Exact 12-Line Baseline Under Primary `r_orth_cal` ([`paired_rorth_12lines.csv`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round11/audit/exploratory_check/outputs/paired_rorth_12lines.csv)):** Recomputing paired `FCT - Mis_s` using frozen `round11/audit/` `r_orth_cal` (`evaluate_paired_cal_diff`, `all_nontarget`, locked `P5` cohorts, R9 leave-one-out `5` sources / R10–R11 `6` sources, `seed = 20261124`, `4,000` joint bootstraps, `97.5%` CIs) gives **`11/12` lines (`6/6` unseen R10+R11 lines passing all `6/6` pairs; `5/6` R9 LOO lines passing all `5/5` pairs; only `nourreddine_kolf21j`, `N=18`, fails `2/5` pairs: `nadig_jurkat` `[-0.0008, +0.1162]` and `replogle_k562_essential` `[-0.0105, +0.0917]`)**:

| Line (`Grp`, `N`) | `FCT r_orth_cal` `[97.5% CI]` | Unpaired `H-D2` | Paired `> 0` | Tightest Paired Contrast (`min CI_lo` or `min pt`, Source) | Source Audit |
|---|---|---:|:---:|---|---|
| `nadig_hepg2` (R9 LOO, `797`) | `+0.1848` `[+0.1772, +0.1929]` | `PASS` | **`5/5`** | `+0.1179` `[+0.1108, +0.1255]` (`k562`; `rpe1 +0.1212`) | `A-T3` |
| `nadig_jurkat` (R9 LOO, `931`) | `+0.1081` `[+0.0975, +0.1185]` | `PASS` | **`5/5`** | `+0.0285` `[+0.0231, +0.0342]` (`k562`) | `A-T3` |
| `nourreddine_kolf21j` (R9 LOO, `18`) | `+0.1142` `[+0.0399, +0.1809]` | `FAIL` | **`3/5`** | `+0.0436` `[-0.0105, +0.0917]` (`k562`; `jurkat -0.0008`) | `A-T3` |
| `replogle_k562_essential` (R9 LOO, `858`) | `+0.1614` `[+0.1506, +0.1720]` | `PASS` | **`5/5`** | `+0.1191` `[+0.1103, +0.1277]` (`hepg2`) | `A-T3` |
| `replogle_rpe1` (R9 LOO, `1068`) | `+0.1510` `[+0.1429, +0.1590]` | `PASS` | **`5/5`** | `+0.0668` `[+0.0614, +0.0723]` (`jurkat`) | `A-T3` |
| `vcc2025_h1` (R9 LOO, `80`) | `+0.1427` `[+0.1168, +0.1707]` | `PASS` | **`5/5`** | `+0.0890` `[+0.0717, +0.1077]` (`kolf21j`) | `A-T3` |
| `dld1` (R10, `86`) | `+0.1355` `[+0.1083, +0.1639]` | `PASS` | **`6/6`** | `+0.0802` `[+0.0626, +0.0994]` (`rpe1`) | `A-T3` + Ctrl |
| `cd4t_rest` (R10, `302`) | `+0.0556` `[+0.0373, +0.0736]` | `FAIL` | **`6/6`** | `+0.0347` `[+0.0262, +0.0431]` (`vcc`; `rpe1 lo +0.0188`) | `A-T3` + Ctrl |
| `calu3` (R10, `8`) | `+0.0820` `[+0.0381, +0.1591]` | `FAIL` | **`6/6`** | `+0.0716` `[+0.0045, +0.1468]` (`vcc`; `rpe1 +0.0451`) | `A-T3` |
| `cd4t_stim48` (R11, `503`) | `+0.0532` `[+0.0346, +0.0725]` | `FAIL` | **`6/6`** | `+0.0165` `[+0.0086, +0.0247]` (`vcc`) | `A-T3` + Ctrl |
| `telohaec` (R11, `9`) | `+0.2173` `[+0.1047, +0.3559]` | `FAIL` | **`6/6`** | `+0.1335` `[+0.0378, +0.2270]` (`hepg2`) | `A-T3` |
| `feng_ipsc` (R11, `93`) | `+0.1212` `[+0.1037, +0.1377]` | `PASS` | **`6/6`** | `+0.0543` `[+0.0407, +0.0674]` (`vcc`) | `A-T3` + Ctrl |

   - **Independent Master Controller Replication (`E1`):** Master Controller (`controller/checks/r12_paired_check.py`, `seed = 777`, `1,000` bootstraps) independently replicated all `4` unseen cohorts above (`dld1`: `FCT +0.1355 [+0.1065, +0.1632]`, min `+0.0802 [+0.0603, +0.0996]` vs `rpe1`; `cd4t_rest`: `FCT +0.0556 [+0.0384, +0.0734]`, min `+0.0347 [+0.0267, +0.0438]` vs `vcc`; `cd4t_stim48`: `FCT +0.0532 [+0.0359, +0.0739]`, min `+0.0165 [+0.0090, +0.0250]` vs `vcc`; `feng_ipsc`: `FCT +0.1212 [+0.1040, +0.1363]`, min `+0.0543 [+0.0422, +0.0671]` vs `vcc`), confirming **`6/6` pairs `> 0`** on each line and noting `vcc2025_h1` (`H1` ESC) is consistently the tightest contrast on `cd4t_rest`, `cd4t_stim48`, and `feng_ipsc` (reason unknown).
   - **Post-Hoc Reading Requiring Prospective Confirmation:** Both T-cell lines (`cd4t_rest`, `cd4t_stim48`) failed unpaired `H-D2` in Rounds 10–11 yet clear `6/6` paired comparisons under `r_orth_cal` (because `FCT` and `Mis_s` share $\overline{\mathrm{FC}}_p$ and co-vary strongly), indicating prior `H-D2` failures reflected low unpaired statistical power rather than lack of specificity. Because this interpretation was formed *after* inspecting R9–R11 results, Round 12 tests paired `H-E1b` (`r_orth_cal`) prospectively.

2. **Split-Half Selection Inflation of Target-Baseline-Free Constant Scaling (`SA-8`, Strictly Empirical Scope [`E6`]):**
   - Even when top-50 responding genes are selected out-of-sample on disjoint split Half A and evaluated on Half B (averaged with `B -> A`, `resp_top50_split`), the deterministic target-baseline-free constant-median transform `E1_ConstMed` retains positive `r_orth_cal_split` on `12/12` lines (`9/12` strict `97.5% CI > 0` for `E1(split) - E1(all)` under locked `P5` cohorts). While **Nicol, Shivakumar & Irizarry (2026)** (`doi:10.64898/2026.05.07.723486`) showed shared-control random noise inflates `top-k` DEG metrics and `Split DEG + control` removes zero-mean noise bias, systematic elevation of deterministic constant-baseline `log1p` scaling under out-of-sample split-half selection had not been documented ([`PRIOR_ART_E.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round11/shared/PRIOR_ART_E.md)).
   - **Strict Non-Mechanistic & Metric-Scope Boundaries (`E6`):** Claim 2 (`SA-8`) is strictly an **empirical observation**, **NOT mechanistic proof**: in Round 11 noiseless real-$\mu$ simulation (`REPORT.md` Check 4), Theorem 1 predicted the wrong sign (`+` instead of empirical `-`) on **`8/12` lines** under `all_nontarget`, leaving the underlying mechanism (finite-sample source noise + biological rewiring + baseline scaling) unexplained; Claim 2 MUST NOT be presented as mechanistic validation. Moreover, Claim 2's applicability is **strictly confined to this study's residual metric `r_orth_cal`** and does **NOT** extrapolate to DEG-based differential-expression metrics analyzed by Nicol et al. (2026).

**Critical Non-Claims (Scoped Out of Round 12):** Consistent with `EXPLORATORY_REVIEW_R11.md` §D, neither the zero-noise Theorem 1 mechanistic attribution (falsified as a complete sign predictor on `all_nontarget`, where real finite-sample source noise + biological rewiring flips `8/12` lines negative while noiseless simulation predicts `+` on `12/12`), nor the `N_5` semi-partial cone projection as a "cure," nor `Q5a` cross-line conservation is claimed as a Round 12 confirmatory hypothesis. All `N_5` metrics are reported strictly descriptively (`§5`).

### 0.2 Key Structural Departures From `PREREG_D` (Mandatory Design & Interpretive Rules)

1. **Perturbation Modality Crossing (`CRISPRi` Sources -> `CRISPRko` Primary Target `A549`):**
   `PREREG_D` excluded Cas9 knockout (`CRISPRko`) datasets to preserve modality homogeneity with Round 9 `CRISPRi` sources. Comprehensive portal screening across IGVF, GEO (through `2026-11-22`), `scPerturb`, and Arc VCA ([`R12_DATA_SCREEN.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round11/shared/R12_DATA_SCREEN.md)) showed that only one unseen public dataset clears target-count power with high certainty (`GSE345058`, `A549` CROP-seq dox-Cas9 **`CRISPRko`**), alongside one borderline **`CRISPRi`** candidate (`XuCao2023.h5ad`, `HEK293`). Round 12 **intentionally tests cross-modality generalization onto `A549` CRISPRko** under pre-locked asymmetric interpretive rules:
   - **If `A549` fails confirmatory testing (`H-E1` or `H-E2`):** Cell-line biology (`A549`) and perturbation modality (`CRISPRko` vs `CRISPRi` sources — e.g., bi-allelic knockout dosage vs. partial transcriptional knockdown, DNA damage response, incomplete editing) are **strictly confounded (`无法区分是细胞系差异还是扰动方式差异`)** and MUST be reported as an unresolvable confound.
   - **If `A549` passes confirmatory testing:** Passing on `A549` provides stronger out-of-domain generalization evidence across BOTH an unseen human lineage and an unseen perturbation modality.
2. **Asymmetric Confirmatory Eligibility (`A549` Mandatory; `HEK293` Count-Gated; `K_conf >= 1` Permitted):**
   - Unlike `PREREG_D`'s uniform `K_conf >= 2` guard, Round 12 locks **before inspecting any expression data**:
     - **`A549` (`GSE345058`) is MANDATORY (`必需线`):** `A549` must clear the Round 11 target/gene count gate (`n_eval_targets >= 50` AND `|G_eval| >= 3,000`) and ingest cleanly. **If `A549` fails the count gate or technical ingestion, Round 12 stops confirmatory testing (`不做确认性判定`)** and reports descriptive results only, even if `HEK293` passes counts.
     - **`HEK293` (`XuCao2023.h5ad`) is FEASIBILITY- & COUNT-GATED (`技术可行性 + 过门槛才算确认性线`):** Upon ingestion (before computing any metric), `HEK293` must first clear the raw-UMI & label feasibility check (`§3.2`, non-negative integer raw UMI counts in `/X` or `/layers` + per-cell `sgRNA` & NTC labels -> else **`EXCLUDED_FEASIBILITY_Q3`**, logged in `round12/shared/R12_QC_COUNTS.md` & `CONFIRMATORY_LINES_R12.txt` without entering confirmatory or secondary metrics). If technically feasible, `HEK293` joins `K_conf` as the second confirmatory line **if and only if** `A549` passes AND `HEK293` independently clears `n_eval_targets >= 50` AND `|G_eval(HEK293)| >= 3,000` (with sufficient NTCs per split half). Otherwise `HEK293` is deterministically classified **`SECONDARY_ONLY`** (`非确认性线`) and reported descriptively.
   - **Resulting `K_conf` Adjudication:**
     - If `K_conf = 2` (`{A549, HEK293}`): Unanimous **`2/2` (`100%`)** passage across both lines is required per claim (IUT across confirmatory lines).
     - If `K_conf = 1` (`{A549}` only): Confirmatory adjudication still executes on `A549`, **AND every summary and manuscript Abstract MUST explicitly state that confirmation relies on a single line (`K_conf = 1`, `A549` CRISPRko) with limited cross-line generalizability (`推广性有限`)**, alongside complete `HEK293` secondary metrics (or exact `EXCLUDED_FEASIBILITY_Q3` reason if excluded pre-metric).

### 0.3 Literature Attribution & Exact Biological Scope

- **Zero-Shot Transfer (`FCT`):** Formula unchanged from `PREREG_A v4 §4.6` / `PREREG_C §3.1` / `PREREG_D §3`; zero novelty claimed for multiplicative fold-change baseline transfer. Expected whole-transcriptome variance explained remains modest ($r_{\mathrm{orth\_cal}} \sim 0.05\text{–}0.19$, $R^2 \sim 0.3\%\text{–}3.6\%$ on R9–R11 lines; `PREREG_D §0.1`).
- **Prior Literature Locks:** **Nicol, Shivakumar & Irizarry (2026)** (`doi:10.64898/2026.05.07.723486`) for finite-sample shared-control DEG selection inflation; **Roohani et al. (2025) / `VCC 2026`** (`doi:10.1016/j.cell.2026.08.004`) for contrast-gating artifacts and split-half anchors; **Shoeibi & Yousefi (2026)** (`arXiv:2608.00152v3`) for cell-count and knockdown confounding (`S7`); **Liu et al. (2026)** (`GSE345058`, `doi:10.64898/2026.06.01.728087`) and **Xu et al. (2023)** (`GSE218566`, `doi:10.1038/s41587-023-01948-9`) for source cohorts.
- **Forbidden Resources Compliance (`E5`):** Strictly zero use of `X-Atlas/Orion` or `Arc Virtual Cell Challenge (VCC)` data, zero account registrations or terms click-throughs, and zero edits to `round10/` or frozen `round11/` pipelines (`FREEZE_MANIFEST_R10.txt` `22/22 OK`, `FREEZE_MANIFEST_R11.txt` **`44/44` manifest lines verified (`41/41` file SHA-256 entries + `3` header lines, `0` failures)**).

### 0.4 Strict Zero-Data-Access Blinding Guarantee

As of submitting `PREREG_E.md`:
1. **`round12/shared/` contains strictly `R12_DOWNLOAD_PLAN.md`** — zero `GSE345058` count matrices (`.mtx.gz`) or `XuCao2023.h5ad` matrices have been downloaded or opened (`round12/shared/raw/` and `round12/shared/lines/` do not exist).
2. Dataset facts, verified byte counts, and metadata columns in `§3` derive strictly from `A-T1`'s metadata-only audit in [`R12_DATA_SCREEN.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round11/shared/R12_DATA_SCREEN.md) and [`R12_DOWNLOAD_PLAN.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/shared/R12_DOWNLOAD_PLAN.md).

---

## §1. Strict Execution Gate (`R12` Order Lock)

No downstream step may start until its prerequisite gate exists in `round12/status/`:

1. **Step 1 (Pre-Registration & Independent Code Prep):** `A-T3` submits `PREREG_E.md` and `COMMIT_REQUEST_A-T3.txt`; `A-T1` prepares `round12/shared/R12_DOWNLOAD_PLAN.md` (exact URLs, byte sizes, SHA-256 verification plan, split-half & NTC locks — zero matrix downloads); `A-T5` delivers `round12/shared/PRIOR_ART_F.md` verifying *Liu et al. 2026* & paired baseline-specificity precedent; `A-T2` prepares `round12/src/` (A549 MTX/CSV reader, HEK293 `.h5ad` reader, paired `FCT - Mis_s`, `H-E2b` common-target `SA-8` joint bootstrap, `N_5`, `SA-Cas9`), testing strictly on Round 9–11 data and synthetic fixtures (alongside `A-T3` independent verification in `round12/audit/`).
2. **Step 2 (`PREREG_E_COMMIT.txt` & `R12_GO.txt`):** Master Controller (`Claude`) audits `PREREG_E.md` + `R12_DOWNLOAD_PLAN.md`, commits and writes the locked commit hash into `round12/status/PREREG_E_COMMIT.txt`, and authorizes download (`R12_GO.txt`). **Zero matrix downloads prior to lock.**
3. **Step 3 (`A-T1` Download Within `15.000 GB` Cap & Counts-Only QC):** `A-T1` downloads ONLY the locked files in `§3.2`, verifies technical ingestion and executes the **`HEK293` pre-metric feasibility check (`§3.2`: non-negative integer raw UMI counts + per-cell `sgRNA`/NTC labels)**, determines which `HEK293` split branch applies (`§3.2` Stage 1 vs Stage 2), enforces `P4`/`P5` eligibility without reading or computing any predictor correlation, and emits counts-only `round12/shared/R12_QC_COUNTS.md`.
4. **Step 4 (Mechanical `CONFIRMATORY_LINES_R12.txt`):** Applying `§2.2` / `§3.2` / `§3.3` rules (`A549` mandatory; `HEK293` `EXCLUDED_FEASIBILITY_Q3` if failing raw-UMI/label gate, else `CONFIRMATORY` iff `n_eval_targets >= 50` and `|G_eval| >= 3,000`, else `SECONDARY_ONLY`) mechanically from `round12/shared/R12_QC_COUNTS.md`, write `round12/status/CONFIRMATORY_LINES_R12.txt` (`K_conf in {0, 1, 2}`).
5. **Step 5 (`CODE_FREEZE_R12.txt`):** Master Controller audits `round12/src/` & `round12/audit/` against `PREREG_E.md` and locks `round12/status/CODE_FREEZE_R12.txt` + `FREEZE_MANIFEST_R12.txt`. Runtime blinding gate blocks any Round 12 evaluation without both `PREREG_E_COMMIT.txt` and `CODE_FREEZE_R12.txt`.
6. **Step 6 (Evaluation & Triple Cross-Check):** `A-T2` runs primary evaluation into `round12/results/`, `A-T3` runs blind independent recomputation into `round12/audit/`, and results (including negative/ambiguous outcomes and `K_conf` / modality disclosures) are reported verbatim.

---

## §2. Source Cohort & Confirmatory Role Logic

### 2.1 Training Source Cohort (`r9_6` Locked Primary)

- **Primary Confirmatory Source Pool (`--source-mode r9_6`):** Strictly the **`6`** harmonized Round 9 human **`CRISPRi`** cell lines (`round9/p2/shared/lines/`, `6,148` genes in `R9_COMMON`, `2,117` line-reproducible `union6_rep` symbols / `1,698` in `R9_COMMON`):
  `nadig_hepg2`, `nadig_jurkat`, `vcc2025_h1`, `replogle_k562_essential`, `replogle_rpe1`, and `nourreddine_kolf21j`.
  Neither `A549` nor `HEK293` exists in `r9_6` (or `r9_r10_8`), ensuring zero target-line leakage.
- **Secondary Source-Expansion Sensitivity (`SA-Src8`, `--source-mode r9_r10_8`):** Adds Round 10 `dld1` (`colon`) and `cd4t_rest` (`CD4+ T`) strictly as a non-confirmatory sensitivity analysis (`§5`).

### 2.2 Target Portfolio Role Assignment (`A549` Mandatory, `HEK293` Count-Gated)

| Dataset Key | Source / Accession | Modality & Lineage | Confirmatory Gate (`n_eval_targets >= 50` & `|G_eval| >= 3,000`) | Pre-Locked Role Rule |
|---|---|---|---|---|
| **`a549_ko`** | **`GSE345058`** (Liu et al. 2026, `doi:10.64898/2026.06.01.728087`) | **`CRISPRko`** (CROP-seq dox-Cas9, lung adenocarcinoma `A549`, completely unseen line & new modality) | **Mandatory Gate:** Failure (`<50` or `<3,000` or ingestion fail) aborts Round 12 confirmatory adjudication (`K_conf = 0` -> descriptive only) | **Primary Confirmatory Line (`必需线`)** |
| **`hek293_xucao`** | **`XuCao2023.h5ad`** (`Zenodo 13350497` / `GSE218566`, Xu et al. 2023) | **`CRISPRi`** (`dCas9-KRAB-MeCP2`, `4sU` PerturbSci-Kinetics, `HEK293`, unseen line) | **Pre-Metric Feasibility + Conditional Count Gate:** Excluded pre-metric (`EXCLUDED_FEASIBILITY_Q3`) if lacking non-negative integer raw UMI counts or `sgRNA`/NTC labels; else enters `K_conf` **only if** `a549_ko` passed AND `hek293_xucao` passes `>=50` targets & `>=3,000` genes (`+` >=15 NTC/half unit) | **Conditional Confirmatory** if gate passes; **`SECONDARY_ONLY`** if feasible & `<50`/`<3000`; **`EXCLUDED_FEASIBILITY_Q3`** if infeasible |

---

## §3. Locked Downloads (`<= 15.000 GB`), Single-sgRNA, NTC, Split-Half, & Normalization Rules

All Round 12 downloads strictly remain under the **`15.000 GB` (`15,000,000,000 B`)** cap — `a549_ko` (`10,997,306,975 B` = `10.997 GB`, of which `GSE345058_raw_counts.mtx.gz` is **`10,956,707,618 B` = `10.957 GB`**) plus `hek293_xucao` (`336,637,209 B` = `0.337 GB`), for a grand total of **`11,333,944,184 B` (`~11.33 GB` / `11.334 GB` / `10.556 GiB`, `5` files, `3.666 GB` headroom)** per [`R12_DOWNLOAD_PLAN.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/shared/R12_DOWNLOAD_PLAN.md) — and occur **only after `PREREG_E_COMMIT.txt` and `R12_GO.txt`**.

### 3.1 `a549_ko` (`A549`, `GSE345058`, `10,997,306,975 B` = `10.997 GB`)

- **Exact Download Files (`.h5ad` `25.05 GB` & `normalized_counts.mtx.gz` `12 GB` strictly excluded):**
  1. `GSE345058_raw_counts.mtx.gz` (**`10,956,707,618 B` = `10.957 GB`**)
  2. `GSE345058_cell_annotations.csv.gz` (`40,152,708 B`)
  3. `GSE345058_gene_annotations.csv.gz` (`378,455 B`)
  4. `GSE345058_sgRNA.csv.gz` (`68,194 B`)
- **Single-sgRNA Rule (`P4`), Join Key, & Feature Exclusion:** Join strictly by unique **`cell_id`** (`606,075` unique rows; **NEVER** join by `cell_barcode`, which collides across the `32` lanes -> `559,858` unique). Verified `MOI = 0.1` + multiplet/Scrublet filtering leaves `606,075 / 606,075` rows with strictly `1` clean `sgRNA` and `Perturbation` entry (`1,000` genes × `4` sgRNAs/gene; `min=73`, `median=565` cells/gene). Strip Ensembl `.version` IDs, map via canonical HGNC, and **strictly exclude synthetic vector row `gene_symbol == "Cas9"` (`genome == "CAS9"`, `ensembl_id == "Cas9"`)**, leaving `29,891` `GRCh38` features.
  - **Primary `CAS9` Dropout Rule (`No Outcome/Transcript Subsetting`):** All QC-annotated singlets are retained regardless of `total_counts_CAS9` (`161,953 / 606,075 = 26.7%` have `total_counts_CAS9 == 0`). Excluding `CAS9 == 0` cells (`a549_ko_cas9pos`, `444,122` cells) is strictly deferred to pre-specified descriptive sensitivity **`SA-Cas9`** (`§5`).
- **NTC Rule:** Restrict control cells strictly to **`Condition == "Nontargeting"`** / **`Perturbation == "Nontargeting"`** (`45,320` cells across `211` NTC sgRNAs: `batch1=11,839`, `batch2=11,111`, `batch3=11,919`, `batch4=10,451`; zero safe-targeting or sequence-overlapping guides).
- **Disjoint Split-Half Partition & `Lane`-Matched NTC Subtraction (`E2`):**
  - **`Half A` = `batch1 + batch3`** (`16` lanes: `Lane1-1..Lane8-1` & `Lane1-3..Lane8-3`; `318,447` cells, `23,758` NTCs).
  - **`Half B` = `batch2 + batch4`** (`16` lanes: `Lane1-2..Lane8-2` & `Lane1-4..Lane8-4`; `287,628` cells, `21,562` NTCs).
  - **Unified `Lane`-Level Subtraction (`E2`):** Subtract **`Lane`-matched** NTC `log1p` control means (NOT `batch`-level means) prior to averaging pseudobulks. All `32` `Lane`s hold `1,244–1,577` NTCs (`>> 15` minimum) and nest strictly inside the `4` batches and inside `Half A` vs `Half B`, guaranteeing zero cross-half control sharing and exact methodological parity with Round 10/11 channel/lane subtraction.

### 3.2 `hek293_xucao` (`HEK293`, `XuCao2023.h5ad`, `Zenodo 13350497` / `GSE218566`, `336,637,209 B` = `0.337 GB`)

- **Exact Download File:** `XuCao2023.h5ad` (`336,637,209 B` = `0.337 GB`, `CC-BY-4.0`, upstream verified **`md5:3efb25a2927c5c228f04f783f18004d8`**; read directly via `h5py` since `anndata` is absent, excluding any non-human spike-in or synthetic lentiviral `dCas9-KRAB-MeCP2` rows).
- **Mandatory Post-Download, Pre-Metric Feasibility Exclusion Rule (`E3`):**
  - Immediately after downloading `XuCao2023.h5ad` and **before computing any model prediction or evaluation metric**, inspect the HDF5 schema (`A-T1`):
    1. Verify whether **integer raw UMI counts** exist in `/X` or any `/layers/*` dataset (required prior to library-size normalization to compute unlogged $\mu_t$ and `FCT`);
    2. Verify whether per-cell **`sgRNA` assignment and explicit NTC labels** exist in `.obs`.
  - **If EITHER condition fails** (no non-negative integer raw UMI matrix in both `/X` and `/layers`, OR missing per-cell `sgRNA` and NTC annotations), **`hek293_xucao` is immediately excluded on technical feasibility grounds (`EXCLUDED_FEASIBILITY_Q3` — enters NEITHER confirmatory adjudication NOR §5 secondary evaluation, with zero downstream metrics computed)**; the exact schema failure reason is logged cleanly in `round12/shared/R12_QC_COUNTS.md` & `round12/status/CONFIRMATORY_LINES_R12.txt`, and Round 12 proceeds with `K_conf = 1` (`A549` only).
- **Author Single-sgRNA Rule & Locked NO Post-Hoc Knockdown Filter:**
  - If feasible, apply the **authors' single-sgRNA assignment rule** as annotated in `.obs` (`nperts == 1` / top sgRNA `>= 60%` of guide UMIs, `>= 3x` 2nd-highest sgRNA, `>= 10` sgRNA UMIs, `>= 300` RNA UMIs, `>= 200` genes; exclude unassigned/multiplet cells).
  - **Pre-Registered Exclusion of Post-Hoc `> 40%` Knockdown Filtering (`Directive §3`):** Primary ingestion **retains ALL QC-valid single-sgRNA cells and sgRNAs** and **MUST NOT** filter cells or guides by `> 40%` on-target knockdown (`>40% KD` filtering uses perturbation response magnitude to select samples — *按结果选样本*). Upon ingestion (`A-T1`), record strictly in `R12_QC_COUNTS.md` whether upstream `XuCao2023.h5ad` holds the full `~126,271` singlet pool or was exported after authors' `~98,315` filter, without applying any further KD subsetting.
- **NTC Rule:** Match strictly on the explicit non-targeting control perturbation label per `scPerturb` `.obs` schema (`15` NTC sgRNAs in library; exclude `"unassigned"` / doublet / empty labels; never pool unassigned cells into NTC).
- **Disjoint Split-Half Partition (`Two-Stage Deterministic Rule Locked Pre-Download`):**
  1. **Stage 1 (`.obs` Replicate Parity):** Inspect `.obs` replicate columns (`replicate`, `rep`, `sample`, `batch`). If `>= 2` distinct replicate IDs exist (e.g., `4` replicates spanning `MOI = 0.1` and `MOI = 0.2`) with `>= 15` NTCs per half, sort replicate labels lexicographically and partition by `1`-based parity: **`Half A` = odd replicates (`1st, 3rd, ...`)** vs **`Half B` = even replicates (`2nd, 4th, ...`)**, subtracting replicate-matched NTC `log1p` means (with half-pooled NTC backup if any single replicate has `< 15` NTCs).
  2. **Stage 2 (Locked Fallback When Replicate Metadata Is Absent):** Because GEO `GSE218566` pooled replicates under `GSM6752591`, if `XuCao2023.h5ad` `.obs` contains `< 2` distinct replicate units, partition all QC-valid single-sgRNA cells using a deterministic **stratified 50/50 split with fixed `np.random.default_rng(42)`** independently within each canonical perturbation label and within the NTC stratum (on barcode-sorted cell indices), centering each half strictly against its own half's NTCs. (`A-T1` records which stage triggered into `R12_QC_COUNTS.md` before unblinding.)

### 3.3 Standardized Processing, Gene Alignment, & `P5` Eligibility (Identical to `PREREG_D §2.3`)

1. **Normalization:** Raw integer UMI counts library-size normalized to $10^4$ per cell (`X_libnorm`) -> $\operatorname{log1p}(X_{\mathrm{libnorm}})$ -> subtract unit-matched (**`Lane` for `A549`**; `replicate`/half for `HEK293`) NTC `log1p` mean prior to averaging per-target pseudobulk $\Delta^{\mathrm{obs}}_{t,p}$. Store unperturbed libnorm control mean $\mu_{t,g}$ (`20,000`-cell stratified NT subsample for fast primary matching + full NT sensitivity `SA-5b`; `15`-NTC minimum unit rule; using full matched NT pool for $\Delta^{\mathrm{obs}}$ and `FCT`) and control covariance $C_t$.
2. **Gene Alignment (`CanonicalHGNC`) & Panel Floor:** Intersect onto `R9_COMMON` (`6,148` symbols) -> $G_{\mathrm{eval}}(t) = G_{\mathrm{R9\_6148}} \cap G_t^{\mathrm{fin}}$. Confirmatory status requires `|G_eval(t)| >= 3,000` (if `HEK293` falls below `3,000` due to shallow `~2,155` UMI/cell depth, it automatically becomes **`SECONDARY_ONLY`** and evaluates on all available finite overlap).
3. **Evaluated Target Eligibility (`P5`, `only_reproducible = True`):** Target $p$ enters evaluation iff: (i) Addendum-B reproducible on $t$ (`n_cells >= 25` single-sgRNA cells, `>= 50` of top-5,000 genes with two-sample Anderson–Darling `BH FDR < 0.05` vs `1,000` stratified NT, `20-PC` energy-distance permutation `p < 0.001` over `2,000` perms vs `500` NT; deduplicated to 1 symbol); (ii) reproducible on `>= 1` source line in `r9_6`; (iii) target self-gene $g_p \in G_{\mathrm{eval}}(t)$ (self-gene masked from masks); (iv) `half_valid` (`>= 2` cells in Half A AND `>= 2` cells in Half B, with `>= 15` NTCs in contributing units / half-pooled backup).

---

## §4. Predictors, Two Confirmatory Claims, & Multiplicity (`IUT + 2-Claim Bonferroni`)

### 4.1 Predictors & Primary Residual Endpoint (`r_orth_cal`)

All predictor formulas follow `PREREG_D §3` identically on `--source-mode r9_6`:
- **`FCT`**: $\hat{\Delta}^{\mathrm{FCT}}_{t,p,g} = \operatorname{log1p}(\mu_{t,g}\overline{\mathrm{FC}}_{p,g}) - \operatorname{log1p}(\mu_{t,g})$.
- **`E1_ConstMed`**: $c_{\mathrm{med}} = \operatorname{median}_{g \in G_{\mathrm{eval}}(t)} \max(\mu_{t,g}, 0)$, $\hat{\Delta}^{\mathrm{E1}}_{t,p,g} = \operatorname{log1p}(c_{\mathrm{med}}\overline{\mathrm{FC}}_{p,g}) - \operatorname{log1p}(c_{\mathrm{med}})$.
- **`NegCtrl_FCT_MismatchedMean(s')`** (`6` mismatched source controls, `s' in r9_6`): replace $\mu_t$ by source control mean $\mu_{s'}$.
- **`NegCtrl_FCT_Perm(k)`** (`100` gene-label permutations, `k = 0..99`, `np.random.default_rng(k).permutation(G)` on $\mu_t$).
- **`Transfer_Mean`** ($\hat{\Delta}^{\mathrm{TMean}}_p = \frac{1}{|\mathcal{S}_p|} \sum_{s \in \mathcal{S}_p} \Delta^{\mathrm{obs}}_{s,p}$): primary residual reference defining $r^{\mathrm{obs}}_p = \Delta^{\mathrm{obs}}_{t,p} - \hat{\Delta}^{\mathrm{TMean}}_p$, $r^{\mathrm{pred}}_p = \hat{\Delta}_{t,p} - \hat{\Delta}^{\mathrm{TMean}}_p$.

Within mask $M_p \subseteq G_{\mathrm{eval}}(t) \setminus \{g_p\}$ (primary: **`all_nontarget`**, $M_p^{\mathrm{all}} = G_{\mathrm{eval}}(t) \setminus \{g_p\}$), residuals are orthogonally projected perpendicular to $\mathrm{span}\{\mathbf{1}, \hat{\Delta}^{\mathrm{TMean}}_p\}$ (`orthogonalize_to_ref`), pooled across perturbations, and disattenuated by the identically masked Spearman–Brown noise ceiling $c^\perp_{\mathrm{ceil}} = \sqrt{2\rho^\perp_{AB}/(1+\rho^\perp_{AB})}$:
$$r^{\mathrm{orth\_cal}} = \operatorname{clip}\!\left(\frac{r^{\mathrm{raw}}_\perp}{c^\perp_{\mathrm{ceil}}},\, -1,\, 1\right).$$
For **`SA-8`** (`resp_top50_split`), top-50 genes per perturbation are ranked out-of-sample on `|half_a_obs[p]|` and evaluated on `half_b_obs[p]` disattenuated by single-half ceiling $\sqrt{\rho^\perp_{AB}(M^{(A)})}$, averaged symmetrically with `B -> A` (`r_orth_cal_split = 0.5 * (r_cal_{A->B} + r_cal_{B->A})`).

### 4.2 Exact Familywise Error Control (`IUT + 2-Claim Bonferroni`, `97.5%` CIs, `seed = 20261201`)

Round 12 tests **`2` separate confirmatory claims**, each structured as an **Intersection–Union Test (IUT)** across conjunctive (`AND`) sub-hypotheses and confirmatory cell lines (`K_conf in {1, 2}`):
1. **Claim 1 (`H-E1 = H-E1a ∩ H-E1b ∩ H-E1c`, across all `K_conf` lines):** Every sub-hypothesis (including all `6` individual paired contrasts `s' = 1..6` inside `H-E1b`) must hold simultaneously (`AND`), so by IUT theory no alpha-splitting is needed within Claim 1.
2. **Claim 2 (`H-E2 = H-E2a ∩ H-E2b`, across all `K_conf` lines):** Both split-half components (`H-E2a` AND `H-E2b`) must hold simultaneously (`AND`), likewise tested via IUT.
3. **2-Claim Bonferroni:** Dividing global $\alpha = 0.05$ evenly across the `2` claims gives $\alpha_{\mathrm{claim}} = 0.025$ per claim, i.e., **two-sided `97.5%` percentile bootstrap confidence intervals (`[1.25%, 98.75%]`)** computed with **`n_boot = 4,000` joint perturbation resamples under new seed `np.random.default_rng(20261201)`**.

### 4.3 Confirmatory Decision Rules (Applied Only When `a549_ko` Passes Count Gate, `K_conf >= 1`)

- **Confirmatory Claim 1 — Paired Target-Line Baseline Specificity (`all_nontarget`):**
  Declared **CONFIRMED (`确认`)** iff ALL three conditions hold simultaneously on **`K_conf / K_conf` (`100%`)** confirmatory lines (`A549`, plus `HEK293` if it enters `K_conf`):
  - **`H-E1a` (De-Amplituded Positive Signal & Addendum A1 Validity):** `FCT` `r_orth_cal` (`all_nontarget`) has **`97.5% CI_lo > 0`** **AND** `boot valid_frac >= 0.95` **AND** `rho_ab_orth > 0`.
  - **`H-E1b` (Paired Specificity Over ALL `6` Mismatched Baselines):** For **EACH** source line $s' \in \text{r9\_6}$ (`6` tests jointly via IUT), the **joint perturbation bootstrap paired difference**
    $$\Delta^{\mathrm{pair}}(s') = r^{\mathrm{orth\_cal}}(\mathrm{FCT}) - r^{\mathrm{orth\_cal}}(\mathrm{NegCtrl\_FCT\_MismatchedMean}(s'))$$
    under `all_nontarget` (sharing identical bootstrap perturbation draws and denominator $c^\perp_{\mathrm{ceil}}$ within each draw) has **`97.5% CI_lo > 0`**.
  - **`H-E1c` (Specificity Over Random Gene Permutations):** Under `all_nontarget`, `FCT` `r_orth_cal` **`97.5% CI_lo > Q_{0.975}(\mathrm{NegCtrl\_FCT\_Perm}(0..99))`**.

- **Confirmatory Claim 2 — Split-Half Top-50 Constant-Baseline Selection Inflation (`SA-8`):**
  Declared **CONFIRMED (`确认`)** independently of Claim 1 iff BOTH conditions hold simultaneously on **`K_conf / K_conf` (`100%`)** confirmatory lines:
  - **`H-E2a` (`E1` Positive Under Split-Half Top-50):** Under `resp_top50_split`, `E1_ConstMed` achieves `r_orth_cal_split` with **`97.5% CI_lo > 0`**.
  - **`H-E2b` (Cross-Mask Elevation Over Whole-Transcriptome `all_nontarget`):** Across the **common perturbation set valid in both `resp_top50_split` and `all_nontarget`** (`split_valid ∩ all_valid`), the joint cross-mask perturbation bootstrap difference
    $$\Delta^{\mathrm{cross}}_{\mathrm{SA8}} = r^{\mathrm{orth\_cal}}_{\mathrm{split}}(\mathrm{E1\_ConstMed}) - r^{\mathrm{orth\_cal}}_{\mathrm{all}}(\mathrm{E1\_ConstMed})$$
    recomputing each mask's noise ceiling ($\sqrt{\rho_{AB}}$ for split; $\sqrt{2\rho_{AB}/(1+\rho_{AB})}$ for all) inside identical perturbation draws, has **`97.5% CI_lo > 0`**.

---

## §5. Strictly Pre-Specified Non-Confirmatory Descriptive & Sensitivity Analyses (`只报告、不判定`)

Every item below is **reported completely alongside confirmatory results without affecting adjudication (`只报告、不判定`)**:

1. **Legacy Unpaired `H-D2` Criterion:** Unpaired check `FCT CI97.5_lo > max_{s'=1..6} Mis_s CI97.5_hi` under `all_nontarget` (`PREREG_D` continuity check).
2. **`FCT - E1_ConstMed` (`all_nontarget`, `H-D3` continuity):** Joint paired bootstrap difference `[97.5% CI]`.
3. **Null-Cone $\mathcal{N}_5(p)$ Semi-Partial Analysis (`r_cone_sp`, `η`, $R^2_{\mathcal{N}_5}$; Strictly Non-Confirmatory [`E1`]):** Because `N_5` was constructed after observing R9–R11 data, **ALL `N_5` outputs — including `N_5` paired `FCT - Mis_s` contrasts (`11/12` in Check 2B) — belong strictly here in `§5` (`只报告、不判定`) and NEVER replace `r_orth_cal` in confirmatory `H-E1b`**. Run `N_5` independently on `all_nontarget` (and `resp_top50_split`) across `FCT`, `6` mismatched `Mis_s` (unpaired & paired vs `FCT`), **single-draw** `Perm0..99` distribution (median, `med|r|`, `Q2.5/Q97.5`, `med_eta` — never expected $\mathbb{E}_{\mathrm{perm}}$), `TMVS`, `Hybrid_ProjShared`, `E3`, `E5`, and construction-nulled controls. Reports MUST **explicitly state that `Transfer_Mean_GlobalScale` ($\eta=0$) and `Source_Grand_Mean` ($\eta=0$) are zero by exact algebraic construction, and `E1_ConstMed`, `E2_ConstP10/P90`, `E6_LogFCBar`, `E7_FCBarMinus1` are suppressed by cone projection ($g_c \in \mathcal{N}_5$), and therefore MUST NOT be counted as empirical validation evidence**.
4. **8-Source Expansion Sensitivity (`SA-Src8`, `--source-mode r9_r10_8`):** Add Round 10 `dld1` & `cd4t_rest` (`8` CRISPRi sources) across `H-E1a/b/c` and `H-E2a/b`.
5. **A549 `CAS9`-Expressing Subset Sensitivity (`SA-Cas9`):** Exclude `total_counts_CAS9 == 0` cells (`~26.7%` of A549 singlets in metadata) from `GSE345058` prior to pseudobulking and Addendum-B QC, re-running `FCT`, `H-E1a/b/c`, and `H-E2a/b` strictly as a descriptive sensitivity check on Cas9 expression dropout vs effect dilution.
6. **In-Sample `resp_top50` Comparison:** Report `E1_ConstMed`, `FCT`, `Perm100` (`median`, `Q0.975`), and `E1(resp_top50) - E1(all_nontarget)` under in-sample `resp_top50` alongside `SA-8`.
7. **Complete Reporting for Non-Confirmatory Lines & Standard R11 Diagnostics:**
   - Whether `HEK293` (`XuCao2023.h5ad`) is classified `CONFIRMATORY` or `SECONDARY_ONLY` (or if `A549` aborts confirmatory adjudication), emit full metrics (`FCT`, paired `FCT - Mis_s`, `Perm100`, `SA-8`, `N_5`, `H-D_util` 6+1 conventions vs `Transfer_Mean`, `SA-1` OLS/Deming slopes, `SA-4/R19` 3-tier variance, `SA-5` S5/full-NT/ablation controls, and `SA-6/S7` `N_cells` & `|KD|` confounder regressions).

---

## §6. Prospective Statistical Power Estimation (Extrapolated via $\sqrt{N}$ to `N ≈ 130` & `N = 50`)

Using observed Round 9–11 unblinded development data (`round11/audit/`, `seed = 20261124`, `4,000` joint bootstraps, $z_{97.5} = 2.2414027$, $z_{0.80} = 0.8416212$), we calibrate the empirical `97.5%` CI half-width $w_{97.5}(N) = \frac{1}{2}(\mathrm{CI}_{97.5,\mathrm{hi}} - \mathrm{CI}_{97.5,\mathrm{lo}}) \approx z_{97.5}\, c / \sqrt{N}$ and extrapolate by $\sqrt{N_{\mathrm{dev}}/N_{\mathrm{new}}}$ to **`A549`'s calibrated expected yield (`N = 130`, matching `126–134` in `R12_DATA_SCREEN.md`)** and **`HEK293`'s confirmatory floor (`N = 50`)**:

| Dev Benchmark Line (`N_dev`) | Hypothesis & Contrast | Observed Effect (`[97.5% CI]`) | Dev Half-Width $w(N)$ (`c` scale) | Projected Half-Width @ `N=130` (Proj. `CI_lo`) | Projected Half-Width @ `N=50` | Required $N_{80}$ |
|---|---|---|---|---|---|---|
| **`dld1` (R10 Cancer, `86`)** | `H-E1b` Bottleneck vs `rpe1` (`min` diff) | **`+0.0802`** `[+0.0626, +0.0994]` | `0.0184` (`c = 0.0761`) | **`0.0150`** (`+0.0652 > 0`) | `0.0241` (`+0.0561`) | **`8.6`** |
| `dld1` (`86`) | `H-E1b` Max CI Width (`vs hepg2`) | `+0.0966` `[+0.0719, +0.1238]` | `0.0259` (`c = 0.1074`) | **`0.0211`** (`+0.0755 > 0`) | `0.0340` (`+0.0626`) | **`11.7`** |
| `dld1` (`86`) | `H-E2a` `E1(split)` | `+0.1294` `[+0.0903, +0.1706]` | `0.0401` (`c = 0.1659`) | **`0.0326`** (`+0.0968 > 0`) | `0.0526` (`+0.0768`) | `15.6` |
| `dld1` (`86`) | `H-E2b` `E1(split) − E1(all)` | `+0.1907` `[+0.1500, +0.2356]` | `0.0428` (`c = 0.1772`) | **`0.0348`** (`+0.1559 > 0`) | `0.0562` (`+0.1345`) | **`8.2`** |
| **`nadig_hepg2` (R9 Cancer, `797`)** | `H-E1b` Widest CI among 5 src (`vs rpe1`; min diff `vs k562` `+0.1179`) | `+0.1212` `[+0.1126, +0.1301]` | `0.0087` (`c = 0.1100`) | **`0.0216`** (`+0.0996 > 0`) | `0.0349` (`+0.0863`) | **`7.8`** |
| `nadig_hepg2` (`797`) | `H-E2b` `E1(split) − E1(all)` | `+0.2447` `[+0.2241, +0.2660]` | `0.0210` (`c = 0.2643`) | **`0.0519`** (`+0.1928 > 0`) | `0.0837` (`+0.1610`) | **`11.1`** |
| **`cd4t_rest` (R10 Primary T, `302`)** | `H-E1b` Tightest Margin (`vs vcc`) | `+0.0347` `[+0.0262, +0.0431]` | `0.0084` (`c = 0.0653`) | **`0.0128`** (`+0.0219 > 0`) | `0.0207` (`+0.0140`) | **`33.6`** |
| `cd4t_rest` (`302`) | `H-E1b` Max $N_{80}$ (`vs rpe1` / `hepg2`) | `+0.0398` / `+0.0486` | `0.0201` / `0.0231` (`c = 0.1559` / `0.1792`) | **`0.0307` / `0.0352`** (`+0.0091` / `+0.0134`) | `0.0494` / `0.0568` *(borderline)* | **`146.1` / `129.1`** |
| `cd4t_rest` (`284` split) | `H-E2b` `E1(split) − E1(all)` | `+0.0383` `[-0.0061, +0.0864]` | `0.0463` (`c = 0.3484`) | `0.0684` *(< 0 lower bound)* | `0.1104` | `782.7` *(underpowered)* |
| **`cd4t_stim48` (R11 Flex Stim, `503`)** | `H-E1b` Tightest Margin (`vs vcc`) | `+0.0165` `[+0.0086, +0.0247]` | `0.0081` (`c = 0.0808`) | **`0.0159`** (`+0.0006 > 0`) | `0.0256` *(-0.0091)* | **`227.0`** |
| `cd4t_stim48` (`503`) | `H-E1b` Widest CI (`vs hepg2` / `k562`) | `+0.0456` / `+0.0356` | `0.0206` / `0.0151` (`c = 0.2058` / `0.1510`) | **`0.0405` / `0.0297`** (`+0.0052` / `+0.0059`) | `0.0652` / `0.0479` | **`193.5` / `171.0`** |
| `cd4t_stim48` (`422` split) | `H-E2b` `E1(split) − E1(all)` | `+0.0730` `[+0.0384, +0.1098]` | `0.0357` (`c = 0.3274`) | **`0.0643`** (`+0.0087 > 0`) | `0.1038` *(-0.0308)* | **`190.5`** |

**Practical Power Conclusions (`A549` vs `HEK293`):**
1. **`A549` (`N ≈ 130`, `~106–212` calibrated range):** Under a cancer-line effect regime (`dld1` / `nadig_hepg2` — expected half-width `w(130) = 0.015–0.022` for paired `FCT - Mis_s`, `~0.033–0.052` for `SA-8`), **`N = 130` achieves `> 99%` power** across every sub-hypothesis of Claims 1 and 2 (`N_80 <= 16`). Under an attenuated primary/cross-modality regime (`cd4t_rest` / `cd4t_stim48`, where `CRISPRko` vs `CRISPRi` difference could narrow margins to `+0.017–+0.046`), `N = 130` sits near **`60%–80%` power** (`N_80 = 129–227`) on the tightest paired contrasts (`H-E1b`) and **`~65%`** (`N_80 = 191`) on `H-E2b` (`cd4t_stim48`-like), or underpowered if `dE1 ~ +0.038` (`cd4t_rest`).
2. **`HEK293` (`N ≈ 29–57` near the `50` floor, shallow `2,155` median UMIs/cell):** If `HEK293` clears `n_eval_targets >= 50`, a `dld1`-magnitude margin clears `> 95%` power (`w(50) = 0.024–0.034`, `N_80 <= 16`), whereas shallow UMI depth can decrease split-half ceiling `c_ceil` (`SE(cal) ~ c_raw / c_ceil`) and raise `N_80`; failed `HEK293` power or `K_conf = 1` fallback will be reported transparently.

---

## §7. Software Governance, Cross-Check Tolerances, & Mandatory Abstract Disclosures

1. **Dual-Team Clean-Room Execution (`C1–C3`):** `round12/src/` (`A-T2`) and `round12/audit/` (`A-T3`) implement the Round 12 readers (`a549_ko`, `hek293_xucao`), `H-E1a/b/c`, and `H-E2a/b` independently, verified on R9–R11 data and unit tests **prior to `CODE_FREEZE_R12.txt`**. Cross-check acceptance criteria (`round12/audit/` vs `round12/results/`): `max_point_abs_diff <= 1e-10` on shared primary metrics (`FCT`, `E1`, `H-E1b` paired diffs, `H-E2a/b`, `Perm100`), `max_ci_abs_diff <= 1e-10` under identical RNG stream (`seed = 20261201`), and `100%` verdict parity.
2. **Mandatory Disclosures in Any Report & Abstract (Non-Negotiable):**
   - **Post-Hoc Hypothesis Provenance (`E6`):** Explicitly disclose that BOTH confirmatory hypotheses (`H-E1` paired specificity and `H-E2` `SA-8`) — as well as the post-hoc reading that prior T-cell unpaired `H-D2` failures reflected low unpaired power — were formulated after inspecting Rounds 9–11 data, which is why Round 12 prospective confirmation is strictly required;
   - **Empirical Non-Mechanistic Scope of Claim 2 (`E6`):** Explicitly disclose that Claim 2 (`SA-8`) is purely an **empirical observation**, **NOT mechanistic validation** (Round 11 noiseless simulation predicted the wrong sign on `8/12` lines under `all_nontarget`, leaving mechanism unexplained), and that Claim 2 scope is **strictly limited to `r_orth_cal` (`resp_top50_split`) and does NOT extrapolate to Nicol et al. (2026) DEG metrics**;
   - Exact line-by-line results (`estimate [97.5% CI]`, `valid_frac`, `rho_ab_orth`) for **`H-E1a`**, all **`6`** individual **`H-E1b`** pairs (`r_orth_cal`), **`H-E1c`**, **`H-E2a`**, and **`H-E2b`**, plus overall Claim 1 & Claim 2 verdicts;
   - Exact **`K_conf`** (`0`, `1`, or `2`), whether `HEK293` passed pre-metric feasibility (`§3.2`: non-negative integer raw UMI counts + `sgRNA`/NTC labels -> else `EXCLUDED_FEASIBILITY_Q3` logged in `round12/shared/R12_QC_COUNTS.md` & `CONFIRMATORY_LINES_R12.txt`), cleared the `>=50` / `>=3,000` confirmatory gate, or became `SECONDARY_ONLY` (plus which split branch fired: replicate parity vs `seed=42`), AND explicit **`K_conf = 1` limited generalizability (`推广性有限`)** disclosure if only `A549` enters `K_conf`;
   - Explicit **`CRISPRi` sources -> `A549` `CRISPRko` modality boundary disclosure**: confounding of cell-line vs modality if `A549` fails, or broader cross-modality scope if `A549` passes;
   - Explicit modest whole-transcriptome variance explained ($R^2$) and full non-confirmatory results (`H-D2` unpaired, `FCT - E1`, `N_5` `§5` including `N_5` paired contrasts, `SA-Src8`, `SA-Cas9`, `HEK293` secondary if feasible & `<50` / `<3000`).

---

## §8. Controller Lock Note (added by Claude, master controller, before lock; no hypothesis, threshold, dataset or mask changed)

- **Prior-art attribution for Claim 1 (`round12/shared/PRIOR_ART_F.md`, item 2: "partially done").** Kawsar et al. 2026 (`CanPert/CASTER`, doi:10.64898/2026.11.15.766815) already compared, per target and paired, a neural model fed the correct target-line control against the same model fed one mismatched control, and found the mismatched input lowers scores. Any report of Claim 1 must cite this. It must not claim to be the first correct-vs-mismatched-control comparison. Claim 1 may be described only as: a preregistered, prospective test that a baseline-only transform (`FCT`) beats each of six mismatched-line baselines on a de-amplituded residual metric (`r_orth_cal`) in a never-used cell line and modality.
- **Prior-art attribution for the dataset (item 1: "not found").** No prior cross-line transcriptomic prediction work using GSE345058 was found. Liu et al. 2026 is cited as the data source only.
- **Timestamps.** Controller documents dated "23:00" or "23:05" on 2026-11-22 (`PREREG_E_AUDIT_R1.md` and the matching agent messages) were actually written about 22:40–22:50 local time; the order of events is unchanged.
