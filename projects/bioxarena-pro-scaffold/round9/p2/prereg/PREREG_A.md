# Pre-Registration A (`PREREG_A`, v3): Zero-Shot Cross-Cell-Line CRISPRi Perturbation Prediction via Target-Line Unperturbed Control Statistics

**Study Team**: BioXArena Team A (`A-T1` Data, `A-T2` Solver/Eval, `A-T3` Pre-Registration & Statistical Audit, `A-T4` Prospective, `A-T5` Literature & Metrics)  
**Human lead**: Zijian Zeng (Round 9 Phase 2, AI-Biology-Lab)  
**Preregistration Version & Status**: Version 3 (Revised per Master Controller Audit `round9/p2/status/PREREG_AUDIT_R2.md`; awaiting remote commit verification in `round9/p2/status/PREREG_COMMIT.txt`)  
**Upstream Reference Rules**: `round9/p2/P2_RULES.md`, `round9/p2/cards/A-T3.md`  
**Prior-Art & Official Metric Specs**: `round9/p2/vcc_linespecific/audit/LIT/{LINE_SPECIFIC_METRIC.md,CIPHER_METHOD.md,PRIOR_ART_A.md,PRIOR_ART_B.md,COMPETITION_MONITOR.md}`  
**Dev Audit Evidence (`Replogle K562 <-> RPE1` Only, Rule 7)**: `round9/p2/vcc_linespecific/audit/dev/{audit_per_line.csv,audit_summary.json,AUDIT_REPORT.md}`

---

## §0 Revision History (v1 -> v2 -> v3)

- **Version 1 (`git commit 222cf75`, `2026-11-21 20:56:39 +08:00`)**: Submitted and pushed before Master Controller Audit R1 arrived (`PREREG_COMMIT.txt` had not yet been written, and **zero held-out cell lines had been evaluated**—only `Replogle K562 <-> RPE1` development checks ran under Rule 3). Preserved unmodified in git history.
- **Version 2 (`sha256 d9273cf6104a08fd498cb3ad7923a664b93696d50922d0a1ace39bb3ae869a48`)**: Addressed all 11 Master Controller requirements in `round9/p2/status/PREREG_AUDIT_R1.md` (**R1–R11**): `"Human lead: Zijian Zeng"` attribution (R1); exclusion of dev lines (`replogle_*`) and $<50$-target line (`nourreddine_kolf21j`, `18` targets) leaving $K_{\mathrm{conf}} = 3$ confirmatory lines (`3/3` pass, `0/3` inferior) (R2, R10); all $K-1 = 5$ mismatched covariance controls (R3); `SrcCorrTgtVar` & `Transfer_Mean_VarScaled` controls and wording locks (R4); strict `AND` across all 6 H-A4 conventions (R5); R6 VCC slope estimator disclosure and scale diagnostics (`r`, `sd_pred/sd_obs`, `slope_deming`, `slope_obs_on_pred`) (R6); CIPHER v2 / Luan et al. (2026) boundary (R7); Rounds 7–8 disclosure (R8); `19:28 +08:00` VCC permission email disclosure (R9); and $r_{\mathrm{cal}}^2$ variance-explained reporting (R11). Verified on `PREREG_AUDIT_R2.md` §1 (`0` held-out evaluations prior to v3).
- **Version 3 (This Document, per `round9/p2/status/PREREG_AUDIT_R2.md`, **R12–R16** + Rule 7)**:
  - **Explicit Post-Dev Disclosure**: **Confirmatory Group B (`Transfer_Mean_VarScaled`, R12) is added after observing development-pair results on `Replogle K562 <-> RPE1`** (where `Transfer_Mean_VarScaled` achieved `r_cal = +0.2198` on `K562` and `+0.2145` on `RPE1`, slope `0.0737` and `0.0267`, i.e., `2–4x` stronger correlation than primary `Hybrid_ProjShared`, `0.0523` / `0.1054`), **strictly before inspecting or evaluating any held-out cell line** (`nadig_hepg2`, `nadig_jurkat`, `vcc2025_h1`, `nourreddine_kolf21j`) per `PREREG_AUDIT_R2.md`.
  - **Rule 7 Dev Artifacts**: Saved complete `--dev-only` audit outputs to `round9/p2/vcc_linespecific/audit/dev/{audit_per_line.csv,audit_summary.json,AUDIT_REPORT.md}` so every dev number in §3.2 traces 1-to-1 to disk artifacts.
  - **R12**: Added **Confirmatory Group B (`H-B1..H-B4`)** testing zero-parameter `Transfer_Mean_VarScaled` in parallel with and separately adjudicated from Group A (`H-A1..H-A4`), including 6 specificity controls (all 5 mismatched source variances + shuffled-gene target variance) and descriptive comparison to `[0.02, 0.13]`.
  - **R13**: Tightened **ALL confirmatory CIs (both Group A and Group B)** from 95% to **97.5% two-sided (`[1.25%, 98.75%]`, Bonferroni correction for 2 hypothesis groups)** and mandated reporting both Group A and Group B verdicts in every abstract and summary.
  - **R14**: Added secondary mechanism test **`Hybrid_on_TMVS`** checking whether gene–gene co-fluctuations add positive `r_cal` increment over marginal-variance scaling alone, locking mechanism phrasing accordingly.
  - **R15**: Added secondary **`FCT`** (Fold-Change Transfer) baseline, **S5** expressed-genes sensitivity (dropping lowest 25% of target control variance), and **S6** unselected `all_nontarget` split-half $\rho_{AB}$ reporting to quantify `resp_top50` ceiling selection bias.
  - **R16**: Locked prior-art compliance on marginal-variance / expression-dependent scaling via `A-T5`'s `audit/LIT/PRIOR_ART_B.md` (and `ALERT_B` flag in `status/A-T5.txt`), forbidding methodological novelty claims if already described in prior literature.
- **Binding Adjudication Rule**: **All confirmatory hypothesis verdicts (Group A `H-A1..H-A4` and Group B `H-B1..H-B4` at 97.5% CI), R4/R14 wording constraints, R16 novelty limits, and independent audits follow Version 3 exclusively.**

---

## 1. Falsifiable Research Questions & Scientific Boundaries

In the Arc Virtual Cell Challenge 2026 Wrap-Up (2026-11-19), top-10 challenge models achieved high bulk/global agreement yet recovered a reported **cell-line-specific perturbation slope of only `0.02–0.13`** (versus **`1.02`** for experimental replicates). Existing cross-context predictors predominantly recover shared cross-line responses while missing cell-line-specific propagation.

**Boundaries with Prior Art (R7 & R16)**:
- **CIPHER v2 (R7)** (`bioRxiv 10.1101/2025.06.27.661814 v2`, Kuznets-Speck et al. 2025/2026) already demonstrated in its abstract and text that removing off-diagonal gene–gene co-fluctuations drops performance several-fold, and that co-fluctuation correlation patterns transfer across studies within the *same* cell type. Accordingly, **"unperturbed control covariance carries perturbation information" is NOT a novel finding claimed by Team A**.
- **Luan et al. (2026)** (`bioRxiv 10.64898/2026.10.20.760514`, 2026-10-25) evaluated centered transfer from `KOLF2.1J` to `H1` using target training perturbations rather than unperturbed target controls $C_t$.
- **Marginal-Variance / Expression-Dependent Scaling Boundary (R16)**: Prior-art screening across bioRxiv, arXiv, PubMed, Scholar, GitHub, VCC 2025/2026 solutions, CASTER, CGC, Systema, and Molina & Zhang (2026) (`10.64898/2026.07.24.740459`) is documented by `A-T5` in `audit/LIT/PRIOR_ART_B.md` (and `status/A-T5.txt`). Regardless of whether `A-T5` raises `ALERT_B`, Group B is tested identically on held-out lines; if `ALERT_B` or prior scaling art is documented, Team A claims **solely** that *"we systematically tested target marginal-variance scaling on zero-shot line-specific residuals"*, with zero claim of inventing a new method.

**Two Parallel Falsifiable Questions (Tested Strictly Zero-Shot with Zero Target Perturbations)**:
1. **Group A (`Hybrid_ProjShared`)**: Can adding target unperturbed single-cell control co-fluctuations $C_t$ onto shared cross-line transfer recover cell-line-specific residuals $\Delta_{p,t} - \bar{\Delta}^{\mathrm{src}}_{p,(-t)}$ strictly above `0.13`, `rho_AB > 0` calibrated signal, all 5 mismatched covariances, and `Transfer_Mean` non-inferiority?
2. **Group B (`Transfer_Mean_VarScaled`, Added Post-Dev per R12)**: Can zero-parameter target control gene-variance scaling $\bar{\Delta}^{\mathrm{src}}_{p,g} \sqrt{C_{t,gg}/\bar{C}_{S,gg}}$ (using zero gene–gene covariance and zero training) recover cell-line-specific residuals above `0.13`, `r_cal > 0`, all 6 variance specificity controls, and `Transfer_Mean` non-inferiority?

---

## 2. Rule 3 Zero-Look, Prior-Round Disclosure (R8), & Rule 5 Clean Data

1. **Zero-Look Prior to Remote Lock**: Prior to Master Controller verification in `round9/p2/status/PREREG_COMMIT.txt`, **no team member or script evaluated, tuned, inspected, or loaded expression matrices (`.npz`) from any held-out cell line** (`nadig_hepg2`, `nadig_jurkat`, `nourreddine_kolf21j`, `vcc2025_h1`). All dev verification (`PREREG_AUDIT_R2.md` §2 and §3.2) used exclusively `Replogle K562 <-> RPE1`.
2. **Prior-Round & Codebase Disclosure (R8)**: Public datasets `nadig_hepg2` and `nadig_jurkat` (`GSE264667`) were analyzed in Rounds 7–8 with non-covariance baselines; repository search confirmed zero prior codebase use of target unperturbed control covariance $C_t$.
3. **Clean Data & Zero X-Atlas (Rule 5)**: Every dataset is open-access (`CC BY 4.0`, `CC0 1.0`, or public NCBI GEO), checksummed in `round9/p2/shared/DATA_SOURCES.md` (`56.597333` GB `<= 60.000` GB cap). Zero X-Atlas / Orion data are used.

---

## 3. Dataset Cohorts, Target-Count Thresholds (R2/R10), & Dev Verification

### 3.1 Full 6-Line Shared Dataset Inventory (`COMMON_GENES = 6,148`)

| Cell Line (`p2/shared/lines/`) | Adjudication Role (**R2/R10**) | Raw Input SHA-256 *(abbreviated if multi-file; exact in `DATA_SOURCES.md`)* | Raw Cells × Genes | NT Ctrl | Final Genes | Dedup Targets (All / Rep) | In `COMMON_GENES` (All / Rep) | LOLO >=1 Src (All / Rep-in-CG) |
|---|---|---|---|---:|---:|---:|---:|---:|
| **`nadig_hepg2`** | **Confirmatory (`>=50`)** | `e1ad7c3c5a201c861a207a858aa7e59f5e6ac1955674c415f7de0d1dadadb52e` | 145,473 × 9,624 | 4,976 | 9,432 | 2,343 / 1,000 | 1,777 / 811 | 2,343 / **799** |
| **`nadig_jurkat`** | **Confirmatory (`>=50`)** | `ffbe15f2c8f7ffcfd7b0ba9e6937d4ebc2d03b0179fa8234648a59bcb82c04a3` | 262,956 × 8,882 | 12,013 | 8,826 | 2,350 / 1,262 | 1,779 / 995 | 2,350 / **934** |
| **`vcc2025_h1`** | **Confirmatory (`>=50`)** | `a0997710...`, `376f0bab...`, `ba5ce667...` | 414,694 × 18,080 *(76,352 dup NT removed)* | 20,000 / 38,176 | 14,656 | 300 / 204 | 228 / 168 | 106 / **80** |
| `replogle_k562_essential` | **Dev Selection *(excluded from Confirmatory)*** | `3e5a63a9e892b21029bb55fca4e12517a49aad7af6c14133ca63d12cf68c6cee` | 310,385 × 8,563 | 10,691 | 8,561 | 2,049 / 1,113 | 1,558 / 910 | 2,039 / 859 |
| `replogle_rpe1` | **Dev Selection *(excluded from Confirmatory)*** | `9b05ef1f81526216fa008d677e9e0d03dce9a2f7a95499a4fb81e505e9d88ef1` | 247,914 × 8,749 | 11,485 | 8,746 | 2,351 / 1,575 | 1,780 / 1,290 | 2,350 / 1,070 |
| `nourreddine_kolf21j` | **Secondary Low-$N$ (`18 < 50`, *excluded*)** | `741674b8...` (`CM`), `d86fa326...` (`ME`) | 38,322 × 34,583 | 17,895 | 15,000 | 95 / 67 | 69 / 51 | 32 / **18** |

- **Confirmatory Held-Out Cohort ($K_{\mathrm{conf}} = 3$)**: Only `nadig_hepg2` (`799`), `nadig_jurkat` (`934`), and `vcc2025_h1` (`80`). Full 6-line LOO-CL is evaluated cleanly, while dev lines (`replogle_*`) and low-$N$ (`nourreddine_kolf21j`, `18 < 50`) are reported separately and excluded from confirmatory verdicts.
- **Preprocessing Summary**: Batch-stratified 50/50 disjoint split-halves (`half_a`, `half_b`, 100% `half_valid`); `reproducible` requires `n_cells >= 25`, `ad_counts >= 50`, `energy_test_p_value < 0.001`; `ctrl_mean_libnorm` stored alongside CP10K `log1p`/`zscore`/`libnorm` deltas.

### 3.2 Locked Primary Configuration & Reproducible Dev Audit (`audit/dev/`, Rule 7)

- **Locked Configuration**: `--space log1p`, `--shrink 0.0`, `cov_max_cells = all` (`4,976–20,000` per line), `mask_self = True`, `reproducible_only = True`; **Group A Primary**: `Hybrid_ProjShared`; **Group B Primary**: `Transfer_Mean_VarScaled`.
- **Pre-Registered Sensitivities**: **S1** `all_targets` (`only_reproducible = False`); **S2** `Hybrid_Raw` vs `Hybrid_ProjShared`; **S3** `shrink = 1e-3`; **S4** `cov_max_cells = 10,000` & `zscore`; **S5** expressed genes only (top 75% target `diag(C_t)`, R15); **S6** unselected `all_nontarget` split-half $\rho_{AB}$ vs `resp_top50` (R15).

**Exact Traceable `Replogle K562 <-> RPE1` Dev Numbers (`7,225` common genes, `832` reproducible shared targets)**  
*(Recorded directly in `round9/p2/vcc_linespecific/audit/dev/{audit_per_line.csv,audit_summary.json,AUDIT_REPORT.md}`)*:

- **`replogle_k562_essential <- replogle_rpe1`** (`rho_AB(top50) = 0.9503`, `c_ceil(top50) = 0.9872`; S6 `rho_AB(all) = 0.8421`, `c_ceil(all) = 0.9562`):
  - **Group A `Hybrid_ProjShared`**: `slope = +0.0050` (`97.5% CI [0.0028, 0.0072]`), `r_cal = +0.0523` (`97.5% CI [0.0302, 0.0740]`), `r_cal_sq = 0.0027`, `sd_pred/sd_obs = 0.0964`, `slope_deming = +0.0050`, `slope_obs_on_pred = +0.5356`. `SrcCorrTgtVar`: `r_cal = +0.0292`. S5 (`5,418` expressed genes): `r_cal = +0.0549`.
  - **Group B `Transfer_Mean_VarScaled`**: `slope = +0.0737` (`97.5% CI [0.0626, 0.0848]`), `r_cal = +0.2198` (`97.5% CI [0.1844, 0.2544]`), `r_cal_sq = 0.0483`, `sd_pred/sd_obs = 0.3396`, `slope_deming = +0.0827`, S5 (`5,418` genes): `r_cal = +0.2306`, `slope = +0.0778`; S6 `r_cal(all) = +0.4659`. **H-B3 Fail Artifact Found**: `tmvs_shuffled_var` achieved `slope = +0.2083` (`97.5% CI [0.1832, 0.2323]`) and `r_cal = +0.3626` (`97.5% CI [0.3143, 0.4102]`) because shuffling $D_t$ against unshuffled large $\bar{D}_S$ amplifies negative regression-to-the-mean shrinkage ($\mathrm{corr}(K562 - RPE1, RPE1) = -0.8325$)!
  - **R14 & FCT**: `Hybrid_on_TMVS` `r_cal = +0.2254` (`delta = +0.0056`, `95% CI [+0.0003, +0.0106]`); `FCT` `r_cal = +0.1394`, `slope = +0.0462`.
- **`replogle_rpe1 <- replogle_k562_essential`** (`rho_AB(top50) = 0.9441`, `c_ceil(top50) = 0.9855`; S6 `rho_AB(all) = 0.7160`, `c_ceil(all) = 0.9135`):
  - **Group A `Hybrid_ProjShared`**: `slope = +0.0085` (`97.5% CI [0.0044, 0.0125]`), `r_cal = +0.1054` (`97.5% CI [0.0526, 0.1602]`), `r_cal_sq = 0.0111`, `sd_pred/sd_obs = 0.0823`, `slope_deming = +0.0086`, `slope_obs_on_pred = +1.2622`. `SrcCorrTgtVar`: `r_cal = +0.1144`. S5 (`5,418` expressed genes): `r_cal = +0.1063`.
  - **Group B `Transfer_Mean_VarScaled`**: `slope = +0.0267` (`97.5% CI [0.0243, 0.0289]`), `r_cal = +0.2145` (`97.5% CI [0.1962, 0.2327]`), `r_cal_sq = 0.0460`, `sd_pred/sd_obs = 0.1262`, `slope_deming = +0.0271`, S5: `r_cal = +0.2327`, `slope = +0.0285`; S6 `r_cal(all) = +0.1214`. `tmvs_shuffled_var` fell to `slope = -0.0144`, `r_cal = -0.0893` (passes H-B3 on RPE1).
  - **R14 & FCT**: `Hybrid_on_TMVS` `r_cal = +0.2351` (`delta = +0.0205`, `95% CI [-0.0041, +0.0443]` crosses zero); `FCT` `r_cal = +0.0903`, `slope = +0.0111`.
- Both Group A and Group B adjudicate **`未确认 (NOT CONFIRMED)`** on dev.

---

## 4. Exact LOO-CL Zero-Shot Predictors & Controls (R3, R4, R12, R14, R15)

For each target cell line $t \in \{1,\dots,6\}$ with source lines $\mathcal{S} = \{1,\dots,6\} \setminus \{t\}$:

1. **Unperturbed Covariance Linear Response**: $v_\ell(p) = C_{\ell,:,g_p} / \max(C_{\ell,g_p,g_p}, 10^{-8})$, with $v_\ell(p)_{g_p} \leftarrow 0$.
2. **Signed Source-Only $\hat{\lambda}_{(-t)}$** and **Shared-Axis Projector $P_\perp = I_G - q_1 q_1^T$** fit strictly on $\mathcal{S}$:
   $$\hat{\lambda}_{(-t)} = \mathrm{clip}\!\left(\frac{\sum_{s \in \mathcal{S}} \sum_{p \in \mathcal{P}_s} \langle \Delta_{p,s},\, v_s(p) \rangle}{\max\!\left(\sum_{s \in \mathcal{S}} \sum_{p \in \mathcal{P}_s} \|v_s(p)\|_2^2,\, 10^{-12}\right)},\; -10.0,\; 10.0\right).$$
3. **Baselines & Group A Models** *(for $\mathcal{S}_p = \{s \in \mathcal{S} : p \text{ available & reproducible}\}$)*:
   - **`Transfer_Mean`**: $\bar{\Delta}^{\mathrm{src}}_{p,(-t)} = \frac{1}{|\mathcal{S}_p|}\sum_{s \in \mathcal{S}_p} \Delta_{p,s}$
   - **`CIPHER_Linear`**: $\hat{\lambda}_{(-t)} v_t(p)$
   - **`Hybrid_Raw`**: $\bar{\Delta}^{\mathrm{src}}_{p,(-t)} + \hat{\lambda}_{(-t)}(v_t(p) - \bar{v}_{\mathcal{S}_p}(p))$
   - **`Hybrid_ProjShared` *(Group A Primary)***: $\bar{\Delta}^{\mathrm{src}}_{p,(-t)} + \hat{\lambda}_{(-t)} P_\perp(v_t(p) - \bar{v}_{\mathcal{S}_p}(p))$
   - **`NegCtrl_HybridProj_MismatchedCov(s')` *(R3, all 5 source lines)***: Replace $v_t(p)$ with $v_{s'}(p)$ for each $s' \in \mathcal{S}$.
   - **`NegCtrl_HybridProj_SrcCorrTgtVar` *(R4a)***: Replace $C_t$ with $D_t^{1/2} \bar{R}_S D_t^{1/2}$ ($\bar{R}_S = \frac{1}{|\mathcal{S}|}\sum_s R_s$).
   - **Additional Covariance Nulls**: `NegCtrl_HybridProj_ShuffledGenes`, `NegCtrl_HybridProj_ShuffledAllCov`, `NegCtrl_HybridProj_DiagCov`.
4. **Group B Predictor & Specificity Controls (R12)**:
   - **`Transfer_Mean_VarScaled` *(Group B Primary, `tmvs`)***:
     $$\hat{\Delta}^{\mathrm{TMVS}}_{p,t,g} = \bar{\Delta}^{\mathrm{src}}_{p,(-t),g} \cdot \sqrt{\frac{\max(C_{t,gg},\, 10^{-8})}{\max\!\left(\frac{1}{|\mathcal{S}|}\sum_{s \in \mathcal{S}} C_{s,gg},\, 10^{-8}\right)}}.$$
   - **`NegCtrl_TMVS_MismatchedVar(s')` *(all 5 source lines)***: Replace $C_{t,gg}$ in numerator with $C_{s',gg}$ for each $s' \in \mathcal{S}$.
   - **`NegCtrl_TMVS_ShuffledVar`**: Permute gene indices of target variance diagonal vector $\max(\mathrm{diag}(C_t), 10^{-8})$ across the $G$ common genes while preserving source mean variances $\bar{C}_{S,gg}$.
5. **Secondary Mechanism Predictor (R14)**:
   - **`Hybrid_on_TMVS`**:
     $$\hat{\Delta}^{\mathrm{Hybrid\_on\_TMVS}}_{p,t} = \hat{\Delta}^{\mathrm{TMVS}}_{p,t} + \hat{\lambda}_{(-t)}\, P_\perp\big(v_t(p) - \bar{v}_{\mathcal{S}_p}(p)\big),$$
     using the exact same $\hat{\lambda}_{(-t)}$ and $P_\perp$ fit on $\mathcal{S}$.
6. **Secondary Fold-Change Transfer Baseline (R15)**:
   - **`FCT`**: Using `ctrl_mean_libnorm` ($\mu_s, \mu_t$):
     $$\mathrm{FC}_{p,s} = \frac{\max\!\big(\mathrm{expm1}(\mathrm{log1p}(\max(\mu_s, 10^{-8})) + \Delta^{\mathrm{log1p}}_{p,s}),\, 10^{-8}\big)}{\max(\mu_s,\, 10^{-8})}, \qquad \overline{\mathrm{FC}}^{\mathrm{src}}_p = \exp\!\left(\frac{1}{|\mathcal{S}_p|}\sum_{s \in \mathcal{S}_p}\log\mathrm{FC}_{p,s}\right),$$
     $$\hat{\Delta}^{\mathrm{FCT}}_{p,t} = \mathrm{log1p}\!\big(\max(\mu_t, 0.0)\cdot \overline{\mathrm{FC}}^{\mathrm{src}}_p\big) - \mathrm{log1p}\!\big(\max(\mu_t, 0.0)\big).$$

---

## 5. Metrics, Slope Scope & Ceiling Selection-Bias Limitations (R5, R6, R11, S5, S6)

1. **Line-Specific Residuals & Primary Metrics**: $r^{\mathrm{obs}}_{p,t} = \Delta^{\mathrm{obs}}_{p,t} - \bar{\Delta}^{\mathrm{src}}_{p,(-t)}$, $r^{\mathrm{pred}}_{p,t} = \hat{\Delta}_{p,t} - \bar{\Delta}^{\mathrm{src}}_{p,(-t)}$ on `resp_top50` (top 50 genes by $|\Delta^{\mathrm{obs}}_{p,g}|$, $g \ne g_p$) and unselected `all_nontarget`.
2. **Scope & Limitations of Slope Metric (R6)**: Neither Arc VCC 2026 Wrap-Up nor `cell-eval2 v0.22.0` publishes the `0.02–0.13` top-10 slope formula (`Dobin 10.64898/2026.11.17.767214` unindexed across Crossref/bioRxiv/OpenAlex/EuropePMC/PubMed APIs as of 2026-11-21). Moreover OLS `slope(r_pred ~ r_obs) = r * SD_pred / SD_obs`. Every line reports secondary `r`, `sd_pred/sd_obs`, orthogonal **Deming slope** ($\delta=1$), reverse OLS `slope_obs_on_pred`, zero-intercept slope, per-perturbation mean slope, and `cos_agreement_top50`.
3. **Noise Ceiling, `r_cal_sq` (R11), & S6 Selection-Bias Disclosure**:
   - From disjoint 50/50 split-halves, $\rho_{AB} = \mathrm{corr}(r^{(A)}, r^{(B)})$, $c_{\mathrm{ceil}} = \sqrt{\frac{2\rho_{AB}}{1+\rho_{AB}}}$, $r_{\mathrm{cal}} = r_{\mathrm{resid}} / c_{\mathrm{ceil}}$, $r_{\mathrm{cal}}^2$ (`r_cal_sq` + joint-bootstrap CI), using **4,000** joint perturbation resamples ($\ge 95\%$ finite draws).
   - **S6 Selection-Bias Limitation**: Selecting `resp_top50` by full-sample $|\Delta^{\mathrm{obs}}|$ conditions on extreme shared sample fluctuations and **upward-biases $\rho_{AB}$** (verified in Replogle dev: `0.944 -> 0.716` on `RPE1`, `0.950 -> 0.842` on `K562`), whereas `r_cal` also depends on out-of-sample vs top-50 residual SNR. Both `resp_top50` and unselected `all_nontarget` ($\rho_{AB}$, $c_{\mathrm{ceil}}$, $r_{\mathrm{cal}}$) MUST be reported side by side in all reports.
4. **S5 Expressed-Genes Sensitivity (R15)**: Drop the lowest 25% of common genes ranked strictly by unperturbed target control variance `diag(C_t) <= Q_{0.25}` (no perturbation info), re-select `resp_top50` within the remaining 75% expressed/variable genes, and recompute `tmvs` and `Hybrid_ProjShared` $r_{\mathrm{cal}}$ and slope to separate abundance/dropout on–off detection effects from continuous magnitude prediction.
5. **Global Accuracy Conventions (R5)**: Simultaneous non-inferiority across all **6** conventions: `ctrl_disc_ownmasked`, `ctrl_disc_panel`, `ctrl_de50_r`, `cent_disc_ownmasked`, `cent_disc_panel`, and `cent_de50_r` (with `ctrl_*` own-target in main tables and others in appendix).

---

## 6. Locked Confirmatory Groups A & B ($K_{\mathrm{conf}} = 3$, 97.5% CI) & Wording Governance

- **Confirmatory Cohort ($K_{\mathrm{conf}} = 3$)**: **`nadig_hepg2`**, **`nadig_jurkat`**, **`vcc2025_h1`**.  
- **Thresholds**: $N_{\mathrm{pass}} = \lceil 0.8 \times 3 \rceil = \mathbf{3/3}$ (`H-A1..A3`, `H-B1..B3`), $N_{\mathrm{inf}} = \lfloor 0.2 \times 3 \rfloor = \mathbf{0/3}$ (`H-A4`, `H-B4`).
- **R13 Multiple-Comparison Level**: **ALL confirmatory criteria below use two-sided 97.5% percentile confidence intervals (`[1.25%, 98.75%]`, Bonferroni over 2 groups)**. **Both Group A and Group B outcomes MUST be reported jointly in every abstract and summary.**

### Confirmatory Group A (`Hybrid_ProjShared`)
1. **H-A1**: On 3/3 confirmatory lines, `Hybrid_ProjShared` `ls_slope_ols_top50` **97.5% CI lower bound $> 0.13$**. *(Pure `CIPHER_Linear` / `Hybrid_Raw` cannot substitute.)*
2. **H-A2**: On 3/3 confirmatory lines, $\rho_{AB} > 0$, $\ge 95\%$ finite joint-bootstrap draws, and `Hybrid_ProjShared` $r_{\mathrm{cal}}$ **97.5% CI lower bound $> 0$**.
3. **H-A3**: On 3/3 confirmatory lines, `Hybrid_ProjShared` beats **all 5** source mismatched covariance controls (`matched_slope_ci975_lo > max_{s'} mismatch_slope_ci975_hi(s')`).
4. **H-A4**: On 0/3 confirmatory lines, `Hybrid_ProjShared - Transfer_Mean` **97.5% CI upper bound $< 0$** across **ALL 6** conventions (`AND`).

### Confirmatory Group B (`Transfer_Mean_VarScaled`, Added Post-Dev per R12)
1. **H-B1**: On 3/3 confirmatory lines, `Transfer_Mean_VarScaled` `ls_slope_ols_top50` **97.5% CI lower bound $> 0.13$**.
2. **H-B2**: On 3/3 confirmatory lines, $\rho_{AB} > 0$, $\ge 95\%$ finite draws, and `Transfer_Mean_VarScaled` $r_{\mathrm{cal}}$ **97.5% CI lower bound $> 0$**.
3. **H-B3 (Specificity Against 6 Variance Controls)**: On 3/3 confirmatory lines, `Transfer_Mean_VarScaled` beats **all 5** `NegCtrl_TMVS_MismatchedVar(s')` AND `NegCtrl_TMVS_ShuffledVar` on **BOTH** `ls_slope_ols_top50` **AND** $r_{\mathrm{cal}}$ with non-overlapping **97.5% CIs** (`tmvs_slope_ci975_lo > max_ctrl_slope_ci975_hi` **AND** `tmvs_cal_r_ci975_lo > max_ctrl_cal_r_ci975_hi`).
4. **H-B4**: On 0/3 confirmatory lines, `Transfer_Mean_VarScaled - Transfer_Mean` **97.5% CI upper bound $< 0$** across **ALL 6** conventions (`AND`).
5. **Descriptive Group B Check**: Report whether `tmvs` slope point estimate lies in `[0.02, 0.13]` with R6 estimator/scale caveat.

### Hard-Coded Decision, Mechanism Wording, & Novelty Governance (R4, R14, R16)
- **Group Verdicts**: Each group (`Group A`, `Group B`) independently yields **`CONFIRMED`** iff all 4 of its criteria (`1 AND 2 AND 3 AND 4`) pass at 97.5% CI; if any criterion in a group fails, that group **MUST state `未确认` (NOT CONFIRMED)** verbatim.
- **R4 Locks (Group A)**:
  1. Use **`"cell-line-specific co-fluctuation structure"`** only if `Hybrid_ProjShared` beats `SrcCorrTgtVar` on BOTH slope and $r_{\mathrm{cal}}$ (disjoint CIs, 3/3 confirmatory lines); otherwise use **`"target-line control statistics"`**.
  2. Unless `Hybrid_ProjShared` beats `Transfer_Mean_VarScaled` on $r_{\mathrm{cal}}$ (3/3 confirmatory lines), state explicitly that **gains over `Transfer_Mean` may be explained by target-line marginal gene variance**.
- **R14 Mechanism Lock**: Only if `Hybrid_on_TMVS - Transfer_Mean_VarScaled` $r_{\mathrm{cal}}$ paired bootstrap 95% CI lower bound is **`> 0` on 3/3 confirmatory lines** may reports state *"co-fluctuation structure provides incremental information beyond marginal variance"*; otherwise they **MUST** state *"In this dataset, co-fluctuation structure exhibits no detectable increment beyond marginal variance"*.
- **R16 Novelty Lock**: Governed by `audit/LIT/PRIOR_ART_B.md` and `status/A-T5.txt` (`ALERT_B`): if prior literature includes expression/variance-dependent scaling, describe Group B strictly as *"a systematic zero-shot test on line-specific residuals"* and claim zero methodological novelty.

---

## 7. Prospective VCC 2026 Protocol (`AC16`, `MOLM13`, `U-87MG`, R9)

On **`2026-11-21 19:28 (+08:00)`**, Human lead (`Zijian Zeng`) emailed **`help@virtualcellchallenge.org`** requesting permission for offline `controls-final` research predictions (`prospective_predict.py`) and pre-GT hash commits (`CELL_EVAL2_REPRO.md` / `HUMAN_STEPS.md`). Agents never download gated data, complete KYC, or push externally. Upon compliant `controls-final` access and post-CC0 GT release, evaluate unmodified predictions via `cell-eval2 v0.22.0` (`--preset vcc2026`) under the exact Version-3 **Group A + Group B (`97.5%` CI) + R4/R14/R15/R16** rules, reporting any failed group verbatim as `未确认`.

---

## 8. Independent Statistical Audit Plan (`A-T3`)

Strictly AFTER `round9/p2/status/PREREG_COMMIT.txt` exists and `A-T2` completes 6-line held-out + sensitivity runs in `round9/p2/vcc_linespecific/results/`:
1. Run standalone `audit/audit_runner.py` (**zero `src/` imports**) with locked primary parameters (`--space log1p --shrink 0.0 --cov-max-cells 0 --project-shared-axis --confirmatory-lines nadig_hepg2 nadig_jurkat vcc2025_h1`), auditing all 6 lines, Group A + Group B at 97.5% CI, R3/R12 all-source & shuffled controls, R4/R14 mechanism/wording checks, R15 `FCT` + `S5` + `S6`, and S1–S4 sensitivity runs.
2. Verify zero-leakage invariance and numerical parity (`tol <= 1e-4`) against `A-T2` via `audit/cross_check.py`, publishing the unvarnished final report in `round9/p2/vcc_linespecific/audit/AUDIT_REPORT.md`.

---

## 9. API-Verified Bibliography

All citations verified directly via Crossref, bioRxiv, or arXiv APIs on 2026-11-21 *(R16 Group B literature expansion recorded separately by `A-T5` in `audit/LIT/PRIOR_ART_B.md`)*:
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
