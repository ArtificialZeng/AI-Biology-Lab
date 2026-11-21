# Pre-Registration A (`PREREG_A`, v2): Zero-Shot Cross-Cell-Line CRISPRi Perturbation Prediction via Target-Line Unperturbed Control Co-Fluctuations

**Study Team**: BioXArena Team A (`A-T1` Data, `A-T2` Solver/Eval, `A-T3` Pre-Registration & Statistical Audit, `A-T4` Prospective, `A-T5` Literature & Metrics)  
**Human lead**: Zijian Zeng (Round 9 Phase 2, AI-Biology-Lab)  
**Preregistration Version & Status**: Version 2 (Revised per Master Controller Audit `PREREG_AUDIT_R1.md`; awaiting remote commit verification in `round9/p2/status/PREREG_COMMIT.txt`)  
**Upstream Reference Rules**: `round9/p2/P2_RULES.md`, `round9/p2/cards/A-T3.md`  
**Prior-Art & Official Metric Specs**: `round9/p2/vcc_linespecific/audit/LIT/{LINE_SPECIFIC_METRIC.md,CIPHER_METHOD.md,PRIOR_ART_A.md,COMPETITION_MONITOR.md}`

---

## §0 Revision History (v1 -> v2)

- **Version 1 (`git commit 222cf75`, `2026-11-21 20:56:39 +08:00`)**: Submitted and pushed before Master Controller Audit R1 arrived (`PREREG_COMMIT.txt` had not yet been written, and **zero held-out cell lines had been evaluated**—only `Replogle K562 <-> RPE1` development checks ran under Rule 3). Version 1 remains preserved in git history without rewriting commits.
- **Version 2 (This Document)**: Revised per all 11 mandatory Master Controller requirements (`round9/p2/status/PREREG_AUDIT_R1.md`, **R1–R11**):
  - **R1**: Standardized role attribution to `"Human lead: Zijian Zeng"` throughout the document.
  - **R2**: Excluded dev selection lines (`replogle_k562_essential`, `replogle_rpe1`) from confirmatory hypothesis counting (still reported in full 6-line tables).
  - **R3**: Expanded `NegCtrl_HybridProj_MismatchedCov(s')` to evaluate **all $K-1 = 5$ source cell lines**, requiring `matched_ci95_lo > max_{s'} mismatch_ci95_hi(s')`.
  - **R4**: Added `NegCtrl_HybridProj_SrcCorrTgtVar` ($D_t^{1/2}\bar{R}_S D_t^{1/2}$) and `Transfer_Mean_VarScaled` ($\bar{\Delta}^{\mathrm{src}} \sqrt{C_{t,gg}/\bar{C}_{S,gg}}$) plus strict wording constraints distinguishing co-fluctuation structure from marginal gene variance.
  - **R5**: Required H-A4 non-inferiority (`0/3` inferior) jointly (`AND`) across **all 6** metric conventions (`ctrl_*` & `cent_*`; both A-T2 own-target-zeroed and A-T3 all-panel-masked `disc_raw`, plus `de50_r`), using A-T2 own-target `ctrl_*` in main tables and others in appendix tables.
  - **R6**: Explicitly disclosed that the Arc VCC 2026 Wrap-Up `0.02–0.13` slope estimator formula is not published in `cell-eval2` or accessible preprints, explained OLS scale dependence ($\beta = r \cdot \mathrm{SD}_{\mathrm{pred}}/\mathrm{SD}_{\mathrm{obs}} \approx r^2$ under MSE calibration), and added secondary `r`, `sd_pred/sd_obs`, `slope_deming`, and `slope_obs_on_pred`.
  - **R7**: Clarified scientific boundary with CIPHER v2 (which already established co-fluctuations > variance alone and within-cell-type transferability) and noted Luan et al. (2026) (`10.64898/2026.10.20.760514`, `KOLF2.1J -> H1`).
  - **R8**: Disclosed Rounds 7–8 prior use of `nadig_hepg2` / `nadig_jurkat` and verified zero prior codebase use of control covariance $C_t$.
  - **R9**: Recorded `2026-11-21 19:28 +08:00` email from Human lead (`Zijian Zeng`) to `help@virtualcellchallenge.org` requesting post-deadline `controls-final` research permission and pre-GT hash publication.
  - **R10**: Applied `>= 50` target floor excluding `nourreddine_kolf21j` (`18` targets) to secondary reporting; locked **$K_{\mathrm{conf}} = 3$** confirmatory lines (`nadig_hepg2`, `nadig_jurkat`, `vcc2025_h1`) with **`3/3`** (`ceil(0.8 * 3)`) for H-A1..H-A3 and **`0/3`** (`floor(0.2 * 3)`) for H-A4.
  - **R11**: Added descriptive metric **`r_cal_sq` ($r_{\mathrm{cal}}^2$) + 95% CI** (proportion of recoverable line-specific variance explained).
- **Binding Adjudication Rule**: **All confirmatory hypothesis verdicts (H-A1..H-A4), mechanism wording constraints, and independent audits follow Version 2 exclusively.**

---

## 1. Falsifiable Research Question & Scientific Boundary

In the Arc Virtual Cell Challenge 2026 Wrap-Up (2026-11-19), top-10 challenge models achieved high bulk/global agreement yet recovered a reported **cell-line-specific perturbation slope of only `0.02–0.13`** (versus **`1.02`** for experimental replicates). Existing cross-context predictors predominantly recover shared cross-line responses while missing cell-line-specific propagation.

**Boundary with Prior Art (R7)**:
- **CIPHER v2** (`bioRxiv 10.1101/2025.06.27.661814 v2`, Kuznets-Speck et al. 2025/2026) already demonstrated in its abstract and text that removing off-diagonal gene–gene co-fluctuations (retaining only diagonal marginal variances) drops performance several-fold, and that co-fluctuation correlation patterns transfer across studies within the *same* cell type. Accordingly, **"unperturbed control covariance carries perturbation information" is NOT a novel finding claimed by Team A**.
- **Luan et al. (2026)** (`bioRxiv 10.64898/2026.10.20.760514`, indexed 2026-10-25) evaluated centered cross-context transfer from `KOLF2.1J` to `H1` (overlapping two cell lines used here; audited alongside VCC 2026 methods in `PRIOR_ART_A.md`), using target-line training perturbations rather than unperturbed control covariance $C_t$.

**Core Falsifiable Question**: Under **strict zero-shot cross-cell-line generalization**—using **zero perturbation labels or perturbed cells from the target cell line**—can adding target unperturbed single-cell control co-fluctuations $C_t$ onto shared cross-line transfer recover the **cell-line-specific perturbation residual** $\Delta_{p,t} - \bar{\Delta}^{\mathrm{src}}_{p,(-t)}$ cleanly above the external `0.13` benchmark and all marginal-variance / mismatched-covariance controls without degrading overall perturbation discrimination (`disc_raw`) or DE correlation (`de50_r`)?

---

## 2. Rule 3 Zero-Look, Prior-Round Disclosure (R8), & Rule 5 Clean Data

1. **Zero-Look Prior to Remote Lock**: Prior to Master Controller verification in `round9/p2/status/PREREG_COMMIT.txt`, **no team member or script evaluated, tuned, inspected, or loaded expression matrices (`.npz`) from any held-out cell line** (`nadig_hepg2`, `nadig_jurkat`, `nourreddine_kolf21j`, `vcc2025_h1`). Configuration selection (`log1p`, `shrink=0.0`, `mask_self=True`, `Hybrid_ProjShared`, `only_reproducible=True`) occurred exclusively on `Replogle K562 <-> RPE1` (`dev_sweep.py` / `--dev-only`).
2. **Prior-Round & Codebase Disclosure (R8)**:
   - Public datasets **`nadig_hepg2`** and **`nadig_jurkat`** (`GSE264667`) were previously analyzed in Rounds 7 and 8 of this project using other cross-cell-line baselines (e.g., S5 cross-line transfer methods).
   - Comprehensive search of the Round 1–8 repository confirmed that **no prior code in this repository used unperturbed target single-cell control covariance $C_t$**.
3. **Clean Data & Zero X-Atlas (Rule 5)**: Every dataset is open-access (`CC BY 4.0`, `CC0 1.0`, or public NCBI GEO), checksummed in `round9/p2/shared/DATA_SOURCES.md` (total new downloads **`56.597333` GB** = `56,597,332,700` B `<= 60.000` GB cap). Zero X-Atlas / Orion data or panels are used.

---

## 3. Dataset Cohorts, Target-Count Thresholds (R2/R10), & Locked Config

### 3.1 Full 6-Line Shared Dataset Inventory (`COMMON_GENES = 6,148`)

| Cell Line (`p2/shared/lines/`) | Adjudication Role (**R2/R10**) | Raw Input SHA-256 *(abbreviated if multi-file; exact in `DATA_SOURCES.md`)* | Raw Cells × Genes | NT Ctrl | Final Genes | Dedup Targets (All / Rep) | In `COMMON_GENES` (All / Rep) | LOLO >=1 Src (All / Rep-in-CG) |
|---|---|---|---|---:|---:|---:|---:|---:|
| **`nadig_hepg2`** | **Confirmatory (`>=50`)** | `e1ad7c3c5a201c861a207a858aa7e59f5e6ac1955674c415f7de0d1dadadb52e` | 145,473 × 9,624 | 4,976 | 9,432 | 2,343 / 1,000 | 1,777 / 811 | 2,343 / **799** |
| **`nadig_jurkat`** | **Confirmatory (`>=50`)** | `ffbe15f2c8f7ffcfd7b0ba9e6937d4ebc2d03b0179fa8234648a59bcb82c04a3` | 262,956 × 8,882 | 12,013 | 8,826 | 2,350 / 1,262 | 1,779 / 995 | 2,350 / **934** |
| **`vcc2025_h1`** | **Confirmatory (`>=50`)** | `a0997710...`, `376f0bab...`, `ba5ce667...` | 414,694 × 18,080 *(76,352 dup NT removed)* | 20,000 / 38,176 | 14,656 | 300 / 204 | 228 / 168 | 106 / **80** |
| `replogle_k562_essential` | **Dev Selection *(excluded from H-A verdict)*** | `3e5a63a9e892b21029bb55fca4e12517a49aad7af6c14133ca63d12cf68c6cee` | 310,385 × 8,563 | 10,691 | 8,561 | 2,049 / 1,113 | 1,558 / 910 | 2,039 / 859 |
| `replogle_rpe1` | **Dev Selection *(excluded from H-A verdict)*** | `9b05ef1f81526216fa008d677e9e0d03dce9a2f7a95499a4fb81e505e9d88ef1` | 247,914 × 8,749 | 11,485 | 8,746 | 2,351 / 1,575 | 1,780 / 1,290 | 2,350 / 1,070 |
| `nourreddine_kolf21j` | **Secondary Low-$N$ (`18 < 50`, *excluded*)** | `741674b8...` (`CM`), `d86fa326...` (`ME`) | 38,322 × 34,583 | 17,895 | 15,000 | 95 / 67 | 69 / 51 | 32 / **18** |

- **Partition Rule (R2 & R10)**:
  - **Confirmatory Held-Out Cohort ($K_{\mathrm{conf}} = 3$)**: Only cell lines untouched before `PREREG_COMMIT.txt` **and** having $\ge 50$ evaluable targets in metadata (`nadig_hepg2`, `nadig_jurkat`, `vcc2025_h1`). All 6 lines are run cleanly as 6-line LOO-CL (with `replogle_*` serving only as source training/dev-reference lines and `nourreddine_kolf21j` reported separately as an underpowered $<50$-target sub-atlas).
- **Preprocessing Summary**: Batch-stratified 50/50 disjoint split-halves (`half_a`, `half_b`; 100% `half_valid = True`); `reproducible` requires `n_cells >= 25`, `ad_counts >= 50`, `energy_test_p_value < 0.001`; HGNC unification via `hgnc_complete_set.txt` (`6,148` common genes).

### 3.2 Locked Primary Configuration & Dev Verification (`Replogle K562 <-> RPE1` Only)

- **Primary Representation (`--repr` / `--space`)**: **`log1p`** (`log1p(1e4 * UMI / libsize)` minus per-batch NT mean)
- **Primary Covariance Shrinkage (`--shrink`)**: **`0.0`** (diagonal floor `max(diag(C), 1e-8)`)
- **Control Sample Scope (`cov_max_cells`)**: All NT control cells in `ctrl_cells.npz` (`4,976–20,000` per line)
- **Self-Target Masking (`mask_self`)**: **`True`** ($v_\ell(p)_{g_p} = 0$)
- **Primary Evaluation Target Filter (`only_reproducible`)**: **`True`** (targets $p$ in `COMMON_GENES` with `reproducible == True` in target line $t$ and $\ge 1$ `reproducible == True` source line $s \ne t$)
- **Primary Model**: **`Hybrid_ProjShared`**
- **Pre-Registered Sensitivity Checks *(secondary, non-confirmatory)***: **S1** `only_reproducible = False` (`all` targets); **S2** `Hybrid_Raw` vs `Hybrid_ProjShared`; **S3** `shrink = 1e-3`; **S4** `cov_max_cells = 10,000` and `zscore`.

**Replogle `K562 <-> RPE1` Dev Verification Summary (including R3/R4/R6/R11)**:  
On `Replogle K562 <-> RPE1` (`7,225` common genes, `832` reproducible shared targets; `1,761` all):
- `log1p` outperformed `zscore` (`Hybrid_ProjShared` slope `0.0011–0.0016`) and `libnorm` (MSE `1.46–1.72`, `120x` $\lambda$ swing).
- On `log1p` (`rep1`): `Hybrid_ProjShared` yielded `ls_slope_ols_top50 = 0.0050` (`K562`) / `0.0085` (`RPE1`), `r_cal = +0.0523 [0.0324, 0.0710]` (`r_cal_sq = 0.0027`) and `+0.1054 [0.0599, 0.1546]` (`r_cal_sq = 0.0111`), `sd_pred/sd_obs = 0.0964` / `0.0823`, `slope_obs_on_pred = +0.5356` / `+1.2622`, `slope_deming = +0.0050` / `+0.0086`.
- However, under **R4**, **`Transfer_Mean_VarScaled` (`tmvs`)** achieved **`r_cal = +0.2198` (`K562`) and `+0.2145` (`RPE1`)** (slope `0.0737` / `0.0267`)—more than double `Hybrid_ProjShared`—and `SrcCorrTgtVar` (`sctv`) reached `r_cal = +0.1144` on `RPE1`! And H-A1 stayed `未确认` (`0.005–0.009 < 0.13`, faithfully preserved per Rule 8).

---

## 4. Exact LOO-CL Zero-Shot Predictors, All Mismatched Controls (R3), & Variance Controls (R4)

For each line $t \in \{1, \dots, 6\}$ with source lines $\mathcal{S} = \{1, \dots, 6\} \setminus \{t\}$:

1. **Unperturbed Covariance Linear Response**: $v_\ell(p) = C_{\ell,:,g_p} / \max(C_{\ell,g_p,g_p}, 10^{-8})$, with $v_\ell(p)_{g_p} \leftarrow 0$.
2. **Signed Source-Only $\hat{\lambda}_{(-t)}$**:
   $$\hat{\lambda}_{(-t)} = \mathrm{clip}\!\left(\frac{\sum_{s \in \mathcal{S}} \sum_{p \in \mathcal{P}_s} \langle \Delta_{p,s},\, v_s(p) \rangle}{\max\!\left(\sum_{s \in \mathcal{S}} \sum_{p \in \mathcal{P}_s} \|v_s(p)\|_2^2,\, 10^{-12}\right)},\; -10.0,\; 10.0\right).$$
3. **Source-Only Shared-Axis Projection**: Orthonormal rank-1 basis (`Q_shared`, $q_1$ from source deltas) defines $P_\perp = I_G - q_1 q_1^T$.
4. **Core Predictors** *(with $\mathcal{S}_p = \{s \in \mathcal{S} : p \text{ available & reproducible}\}$, $\bar{\Delta}^{\mathrm{src}}_{p,(-t)} = \frac{1}{|\mathcal{S}_p|}\sum_{s \in \mathcal{S}_p} \Delta_{p,s}$, $\bar{v}_{\mathcal{S}_p}(p) = \frac{1}{|\mathcal{S}_p|}\sum_{s \in \mathcal{S}_p} v_s(p)$)*:
   - **`Transfer_Mean`**: $\bar{\Delta}^{\mathrm{src}}_{p,(-t)}$
   - **`CIPHER_Linear`**: $\hat{\lambda}_{(-t)}\, v_t(p)$
   - **`Hybrid_Raw`**: $\bar{\Delta}^{\mathrm{src}}_{p,(-t)} + \hat{\lambda}_{(-t)}\big(v_t(p) - \bar{v}_{\mathcal{S}_p}(p)\big)$
   - **`Hybrid_ProjShared` *(Primary Locked)***: $\bar{\Delta}^{\mathrm{src}}_{p,(-t)} + \hat{\lambda}_{(-t)}\, P_\perp\big(v_t(p) - \bar{v}_{\mathcal{S}_p}(p)\big)$
5. **All $K-1$ Mismatched-Covariance Specificity Controls (R3)**:
   - For **every source cell line** $s' \in \mathcal{S}$ (`5` lines in 6-line LOO):
     $$\text{NegCtrl\_HybridProj\_MismatchedCov}(s'):\quad \bar{\Delta}^{\mathrm{src}}_{p,(-t)} + \hat{\lambda}_{(-t)}\, P_\perp\big(v_{s'}(p) - \bar{v}_{\mathcal{S}_p}(p)\big).$$
   - All 5 mismatched slopes and 95% CIs are computed and reported; H-A3 requires matched CI lower bound strictly above **all 5** mismatched CI upper bounds (`matched_ci95_lo > max_{s'} mismatch_ci95_hi(s')`).
6. **Target Marginal Variance vs Co-Fluctuation Controls (R4, *Secondary / Wording-Governing*)**:
   - **(a) `NegCtrl_HybridProj_SrcCorrTgtVar`**: Form target correlation-replaced covariance $C^{\mathrm{SCTV}}_t = D_t^{1/2} \bar{R}_S D_t^{1/2}$, where $D_t = \mathrm{diag}(\max(\mathrm{diag}(C_t), 10^{-8}))$ and $\bar{R}_S = \frac{1}{|\mathcal{S}|}\sum_{s \in \mathcal{S}} R_s$ ($R_s$ = Pearson correlation matrix of $C_s$), then apply `Hybrid_ProjShared` using $C^{\mathrm{SCTV}}_t$ in place of $C_t$.
   - **(b) `Transfer_Mean_VarScaled`**: Pure gene-level marginal variance scaling without gene–gene covariance:
     $$\hat{\Delta}^{\mathrm{TMVS}}_{p,t,g} = \bar{\Delta}^{\mathrm{src}}_{p,(-t),g} \cdot \sqrt{\frac{\max(C_{t,gg},\, 10^{-8})}{\max\!\left(\frac{1}{|\mathcal{S}|}\sum_{s \in \mathcal{S}} C_{s,gg},\, 10^{-8}\right)}}.$$
7. **Additional Null Controls**: `NegCtrl_HybridProj_ShuffledGenes` (permute $C_t$ gene labels), `NegCtrl_HybridProj_ShuffledAllCov` (permute both $C_t$ and $C_s$ in A-T3 audit), and `NegCtrl_HybridProj_DiagCov` ($C_t \leftarrow \mathrm{diag}(C_t)$).

---

## 5. Metrics, Slope Limitations (R6), Variance Explained (R11), & Noise Ceiling

1. **Line-Specific Residuals & Primary Slope (`ls_slope_ols_top50`)**: $r^{\mathrm{obs}}_{p,t} = \Delta^{\mathrm{obs}}_{p,t} - \bar{\Delta}^{\mathrm{src}}_{p,(-t)}$, $r^{\mathrm{pred}}_{p,t} = \hat{\Delta}_{p,t} - \bar{\Delta}^{\mathrm{src}}_{p,(-t)}$, evaluated across pooled pairs in `resp_top50` (top 50 genes by $|\Delta^{\mathrm{obs}}_{p,g}|$ excluding $g_p$) and secondary `all_nontarget`.
2. **Scope & Limitations of Slope Metric (R6)**:
   - **Verified External Gap**: Neither the Arc VCC 2026 Wrap-Up nor `ArcInstitute/cell-eval2 v0.22.0` publishes the estimator equation for the reported `0.02–0.13` top-10 line-specific slope (`Dobin 10.64898/2026.11.17.767214` unindexed across Crossref/bioRxiv/OpenAlex/EuropePMC/PubMed APIs as of 2026-11-21). H-A1 therefore benchmarks against an external figure whose exact code estimator is unverified.
   - **Scale Dependence**: OLS `slope(r_pred ~ r_obs)` equals $r \cdot \frac{\mathrm{SD}(r^{\mathrm{pred}})}{\mathrm{SD}(r^{\mathrm{obs}})}$, so an MSE-optimal scalar prediction achieves slope $\approx r^2$. Consequently, every line **must also report**: pooled Pearson **`r`** (`resid_corr`), **`sd_pred / sd_obs`**, orthogonal equal-error-variance **Deming slope** (`slope_deming`, $\delta = 1$), **reverse OLS slope** (`slope_obs_on_pred`), zero-intercept slope, per-perturbation mean slope, and `cos_agreement_top50`.
3. **Noise Ceiling & `r_cal_sq` (R11)**:
   - From disjoint 50/50 split-halves (`half_a`, `half_b`), $\rho_{AB} = \mathrm{corr}(r^{(A)}, r^{(B)})$; if $\rho_{AB} > 0$, $\rho_{\mathrm{SB}} = \frac{2\rho_{AB}}{1+\rho_{AB}}$, $c_{\mathrm{ceil}} = \sqrt{\rho_{\mathrm{SB}}}$, and $r_{\mathrm{cal}} = r_{\mathrm{resid}} / c_{\mathrm{ceil}}$. If $\rho_{AB} \le 0$ on full sample, H-A2 fails automatically.
   - **Joint Perturbation Bootstrap (`n_boot = 4000`)**: Resample perturbation indices $I^{(b)}$ identically for both $r_{\mathrm{resid}}^{(b)}$ and $\rho_{AB}^{(b)} \to c_{\mathrm{ceil}}^{(b)}$; require $\ge 95\%$ finite draws.
   - **Explainable Line-Specific Variance Fraction (`r_cal_sq`, R11)**: Report $r_{\mathrm{cal}}^2$ (`r_cal_sq` + 95% joint-bootstrap CI, plus signed $\mathrm{sign}(r_{\mathrm{cal}}) r_{\mathrm{cal}}^2$; interpret variance explained only when `r_cal` CI > 0) on every line for transparent reporting.
4. **Global Accuracy Metrics (H-A4, R5)**:
   - Tie-aware `disc_raw` evaluated in **both** A-T2's own-target-zeroed convention (`main table`) and A-T3's all-panel-targets-masked convention (`appendix`, plus `pds_cosine_panel`), alongside `de50_r` (own target masked), evaluated under **both** unperturbed control (`ctrl_*`, main table) and `Source_Grand_Mean` centroid (`cent_*`, appendix).

---

## 6. Locked Confirmatory Adjudication ($K_{\mathrm{conf}} = 3$) & Wording Governance (R4)

Confirmatory Cohort ($K_{\mathrm{conf}} = 3$): **`nadig_hepg2`**, **`nadig_jurkat`**, **`vcc2025_h1`**.  
Applying `A-T3.md` proportions with ceiling/floor rounding (`only stricter`):  
$$N_{\mathrm{pass}} = \lceil 0.8 \times 3 \rceil = \mathbf{3/3}, \qquad N_{\mathrm{inf}} = \lfloor 0.2 \times 3 \rfloor = \mathbf{0/3}.$$

1. **H-A1 (Line-Specific Slope > 0.13)**: On **3/3 confirmatory lines**, locked primary **`Hybrid_ProjShared`** achieves `ls_slope_ols_top50` **95% CI lower bound $> 0.13$**. *(Pure `CIPHER_Linear` and `Hybrid_Raw` cannot substitute for `Hybrid_ProjShared` due to $-\bar{\Delta}^{\mathrm{src}}$ / $-\bar{v}_{\mathcal{S}_p}$ artifacts verified in Replogle dev.)*
2. **H-A2 (Noise-Calibrated Correlation > 0)**: On **3/3 confirmatory lines**, $\rho_{AB} > 0$, $\ge 95\%$ finite joint-bootstrap draws, and $r_{\mathrm{cal}}$ **joint-bootstrap 95% CI lower bound $> 0$**.
3. **H-A3 (Specificity Against ALL 5 Source Covariances, R3)**: On **3/3 confirmatory lines**, `Hybrid_ProjShared` beats **every** mismatched source covariance ($s' \in \mathcal{S}$, all 5 lines): `matched_ci95_lo > max_{s'} mismatch_ci95_hi(s')`.
4. **H-A4 (Strict Multi-Convention Non-Inferiority, R5)**: On **0/3 confirmatory lines**, `Hybrid_ProjShared - Transfer_Mean` has paired bootstrap 95% CI strictly below zero (`diff_ci95_hi < 0`) across **ALL SIX** conventions simultaneously (`AND`):
   `ctrl_disc_ownmasked`, `ctrl_disc_panel`, `ctrl_de50_r`, `cent_disc_ownmasked`, `cent_disc_panel`, and `cent_de50_r`.

### Hard-Coded Decision & Wording Rules
- **Overall Verdict**: **`CONFIRMED`** iff **H-A1 AND H-A2 AND H-A3 AND H-A4** all pass on 3/3 confirmatory lines (`0/3` inferior on H-A4). If **ANY** criterion fails, all conclusions **MUST state `未确认` (NOT CONFIRMED)** verbatim. Post-hoc changes to parameters, hypotheses, cohorts, or thresholds are forbidden.
- **Locked R4 Wording Rules (Binding Regardless of Pass/Fail)**:
  1. Only if **`Hybrid_ProjShared`** beats `(a) NegCtrl_HybridProj_SrcCorrTgtVar` on **BOTH** `ls_slope_ols_top50` and `noise_cal_corr` with non-overlapping 95% CIs on **3/3 confirmatory lines** may papers/summaries refer to **`"cell-line-specific co-fluctuation structure"`**; otherwise conclusions **MUST** use **`"target-line control statistics"`** (or equivalent non-structural phrasing).
  2. If `Hybrid_ProjShared` does not outperform `(b) Transfer_Mean_VarScaled` on `noise_cal_corr` (`matched_cal_ci95_lo > tmvs_cal_ci95_hi`) across **3/3 confirmatory lines**, conclusions **MUST** state explicitly that **observed improvements may be explained by target-line marginal gene variance**.

---

## 7. Prospective VCC 2026 Generalization Protocol (`AC16`, `MOLM13`, `U-87MG`, R9)

Aligned with `round9/p2/vcc_linespecific/src/prospective/{HUMAN_STEPS.md,CELL_EVAL2_REPRO.md,vcc2026_format.py,prospective_predict.py}` and `P2_RULES.md` Rule 6:
1. **Compliance & Outreach Update (R9)**: On **`2026-11-21 19:28 (+08:00)`**, Human lead (`Zijian Zeng`) emailed **`help@virtualcellchallenge.org`** requesting written clarification/permission to run offline research predictions on post-deadline `controls-final` and publish prediction hashes prior to GT release; currently awaiting response. Agents strictly **do not** download gated Challenge files, complete KYC verification, or push externally.
2. **Pre-GT Cryptographic Lock**: Once compliant `controls-final` (via explicit written permission or CC0 1.0 GCS release prior to GT release) are obtained by Human lead, `prospective_predict.py` runs the exact locked configuration (`log1p`, `shrink=0.0`, `mask_self=True`, `Hybrid_ProjShared`) fitted strictly on our 6 public lines + target `non-targeting` controls (`18,533` panel genes, `400` cells/pert/context). Prior to GT release, commit only `SHA256SUMS.txt`, `prediction_manifest.json`, `TIMESTAMP_DECLARATION.md`, and code (never restricted input matrices) to public GitHub.
3. **Post-GT Evaluation**: Score unmodified predictions upon CC0 GT release with **`ArcInstitute/cell-eval2 v0.22.0`** (`--preset vcc2026`) and the identical Version-2 **H-A1..H-A4 + R3/R4/R5/R6/R11** protocol, reporting negative results verbatim as `未确认`.

---

## 8. Independent Statistical Audit Plan (`A-T3`)

After `round9/p2/status/PREREG_COMMIT.txt` appears and `A-T2` completes evaluation in `round9/p2/vcc_linespecific/results/`:
1. Execute standalone `audit/audit_runner.py` (**zero imports from `src/`**) with `--space log1p --shrink 0.0 --cov-max-cells 0 --project-shared-axis --confirmatory-lines nadig_hepg2 nadig_jurkat vcc2025_h1 --status-dir round9/p2/status`, evaluating all 6 lines, all $5$ source mismatches per line (R3), both R4 controls, all 6 H-A4 sub-metrics (R5), R6 scale diagnostics, and R11 $r_{\mathrm{cal}}^2$, plus S1–S4 sensitivity checks.
2. Verify zero-leakage invariance and numerical parity (`tol <= 1e-4`) against `A-T2` outputs via `audit/cross_check.py`, writing complete unvarnished findings and wording-governance flags to `round9/p2/vcc_linespecific/audit/AUDIT_REPORT.md`.

---

## 9. API-Verified Bibliography

All citations verified directly via Crossref, bioRxiv, or arXiv APIs on 2026-11-21:
1. **Kuznets-Speck B et al. (19 authors)** (2025/2026). *Co-fluctuations of genes in single cells predict transcriptome-wide outcomes to perturbations* (`CIPHER`). bioRxiv `doi:10.1101/2025.06.27.661814` (v1: 2025-07-01; v2: 2026-09-03), PMID `40631127`.
2. **Replogle JM et al. (17 authors)** (2022). *Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq*. *Cell*, `doi:10.1016/j.cell.2022.05.013`, PMID `35688146`.
3. **Nadig A et al. (8 authors)** (2025). *Transcriptome-wide analysis of differential expression in perturbation atlases*. *Nature Genetics*, `doi:10.1038/s41588-025-02169-3`, PMID `40259084`.
4. **Nourreddine S et al. (22 authors)** (2026). *A genome-scale CRISPRi perturbation atlas of human induced pluripotent stem cells*. *Nature Biotechnology*, `doi:10.1038/s41587-026-03199-w`, PMID `42386990`.
5. **Miller HE et al. (6 authors)** (2026). *Deep learning perturbation models can outperform baselines on calibrated metrics*. *Nature Biotechnology*, `doi:10.1038/s41587-026-03307-w`.
6. **Adduri A et al. (81 authors)** (2026). *Virtual Cell Challenge 2026: Benchmarking zero-shot generalization across cellular contexts*. *Cell*, `doi:10.1016/j.cell.2026.08.004`.
7. **Adduri AK et al. (32 authors)** (2026). *Predicting cellular responses to perturbation across diverse contexts with State*. *Cell*, `doi:10.1016/j.cell.2026.07.052`.
8. **Viñas Torné R et al. (8 authors)** (2025). *Systema: a framework for evaluating genetic perturbation response prediction beyond systematic variation*. *Nature Biotechnology*, `doi:10.1038/s41587-025-02777-8`.
9. **Kawsar M et al. (6 authors)** (2026). *Transferring CRISPRi knockdown responses to new cell lines requires screening part of each line* (`CASTER`). bioRxiv `doi:10.64898/2026.11.15.766815` (2026-11-17).
10. **Huang Y et al. (4 authors)** (2026). *The recoverable resolution of cellular perturbation-response prediction* (`CGC`). bioRxiv `doi:10.64898/2026.09.30.755756` (2026-10-01).
11. **Dhaka B et al. (3 authors)** (2026). *HyperMap: An Efficient Framework for Transferring Perturbation Responses Across Diverse Biological Contexts*. bioRxiv `doi:10.64898/2026.04.23.720505` (2026-04-27).
12. **Luan et al. (4 authors)** (2026). *Centered cross-context perturbation transfer on the public 2025 Virtual Cell Challenge benchmark*. bioRxiv `doi:10.64898/2026.10.20.760514` (2026-10-25).
13. **Molina et al. (2 authors)** (2026). *Perturbation response decomposition enables biologically aligned generalization to unseen perturbations and cellular contexts*. bioRxiv `doi:10.64898/2026.07.24.740459` (2026-07-27).
14. **Gao J** (2026). *The Delta Gap: Cellular World Models for Perturbation Response Prediction Are Far Harder Than Global Correlation Implies*. `arXiv:2611.10658v1` (2026-11-07).
15. **Zhang Z, Huang Y, Feng F et al. (4 authors)** (2026). *A low-rank decoder bottleneck bounds reliability in foundation-model perturbation prediction*. `arXiv:2611.05849v1` (2026-11-04).
16. **Arc Institute (`cell-eval2` & VCC 2026 Wrap-Up)**: `ArcInstitute/cell-eval2` `v0.22.0` (`--preset vcc2026`) and Arc VCC 2026 Wrap-Up (`2026-11-19`). *Unindexed DOI Disclosure*: Dobin metrics preprint (`10.64898/2026.11.17.767214`) returned HTTP 404 / 0 hits via Crossref, bioRxiv, OpenAlex, EuropePMC, and PubMed APIs as of 2026-11-21 (`www.biorxiv.org` HTML Cloudflare-gated, not bypassed).
