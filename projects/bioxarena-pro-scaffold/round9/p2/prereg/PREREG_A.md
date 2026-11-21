# Pre-Registration A (`PREREG_A`, v4): Zero-Shot Cross-Cell-Line CRISPRi Perturbation Prediction via Target-Line Unperturbed Control Statistics

**Study Team**: BioXArena Team A (`A-T1` Data, `A-T2` Solver/Eval, `A-T3` Pre-Registration & Statistical Audit, `A-T4` Prospective, `A-T5` Literature & Metrics)  
**Human lead**: Zijian Zeng (Round 9 Phase 2, AI-Biology-Lab)  
**Preregistration Version & Status**: Version 4 (Final Lock Revision per Master Controller Audit `round9/p2/status/PREREG_AUDIT_R3.md`; awaiting remote commit verification in `round9/p2/status/PREREG_COMMIT.txt`)  
**Upstream Reference Rules**: `round9/p2/P2_RULES.md`, `round9/p2/cards/A-T3.md`  
**Prior-Art & Official Metric Specs**: `round9/p2/vcc_linespecific/audit/LIT/{LINE_SPECIFIC_METRIC.md,CIPHER_METHOD.md,PRIOR_ART_A.md,PRIOR_ART_B.md,COMPETITION_MONITOR.md}`  
**Dev Audit Evidence (`Replogle K562 <-> RPE1` Only, Rule 7)**: `round9/p2/vcc_linespecific/audit/dev/{audit_per_line.csv,audit_summary.json,AUDIT_REPORT.md}`

---

## §0 Revision History (v1 -> v2 -> v3 -> v4)

- **Version 1 (`git commit 222cf75`, `2026-11-21 20:56:39 +08:00`)**: Submitted and pushed before Master Controller Audit R1 arrived (`PREREG_COMMIT.txt` had not yet been written, and **zero held-out cell lines had been evaluated**—only `Replogle K562 <-> RPE1` development checks ran under Rule 3). Preserved unmodified in git history.
- **Version 2 (`git commit 11aca10`, `2026-11-21 21:11:48 +08:00`, `sha256 d9273cf6104a08fd498cb3ad7923a664b93696d50922d0a1ace39bb3ae869a48`)**: Addressed **R1–R11** (`PREREG_AUDIT_R1.md`): `"Human lead: Zijian Zeng"` attribution (R1); exclusion of dev lines (`replogle_*`) and $<50$-target line (`nourreddine_kolf21j`, `18` targets) leaving $K_{\mathrm{conf}} = 3$ confirmatory lines (`3/3` pass, `0/3` inferior) (R2, R10); all $K-1 = 5$ mismatched covariance controls (R3); `SrcCorrTgtVar` & `Transfer_Mean_VarScaled` controls and wording locks (R4); strict `AND` across all 6 H-A4 conventions (R5); R6 VCC slope estimator disclosure and scale diagnostics (R6); CIPHER v2 / Luan et al. (2026) boundary (R7); Rounds 7–8 disclosure (R8); `19:28 +08:00` VCC permission email disclosure (R9); and $r_{\mathrm{cal}}^2$ variance-explained reporting (R11). Preserved unmodified in git history.
- **Version 3 (`sha256 cce7e2855e87bea47861542d2b7fe91366e54b41427f191e568ba8f7c27cc528`)**: Addressed **R12–R16** + Rule 7 (`PREREG_AUDIT_R2.md`): post-dev disclosure adding Confirmatory Group B (`Transfer_Mean_VarScaled`, `R12`), Bonferroni 97.5% two-sided CIs across both groups (`R13`), secondary `Hybrid_on_TMVS` increment test (`R14`), secondary `FCT` baseline + S5/S6 (`R15`), and prior-art boundary lock via `PRIOR_ART_B.md` (`R16`). Verified clean in `PREREG_AUDIT_R3.md` §1 (`0` held-out evaluations prior to v4).
- **Version 4 (This Document, per `round9/p2/status/PREREG_AUDIT_R3.md`, **R17–R20** + `PRIOR_ART_B.md` §5.2)**:
  - **Explicit Post-Dev Whole-Line Amplitude Artifact Disclosure (`v3 -> v4`)**: After Master Controller verification on `Replogle K562 <-> RPE1` development runs (`PREREG_AUDIT_R3.md`)—**strictly before touching or evaluating any held-out cell line**—we observed that scaling the shared source response by **a single whole-line amplitude scalar per target line (`Transfer_Mean_GlobalScale`)** alone achieved raw residual correlation **`r_resid = +0.5546` (`K562`)** and **`+0.3034` (`RPE1`)** (`r_cal = +0.5618` / `+0.3078`), explaining why raw un-orthogonalized $r_{\mathrm{cal}}$ mixes a 1-scalar whole-line shrinkage artifact with cross-gene residual structure. Notice that even single-seed shuffled-variance `tmvs` preserved denominator scaling $\propto 1/\sqrt{\bar{C}_{S,gg}}$ (`K562 perm100 median r_orth_cal = +0.0902`). To isolate genuine gene-pattern specificity cleanly before held-out unblinding:
  - **R17**: Defined per-perturbation amplitude orthogonalization (`r_orth_cal`, $\rho_{AB}^{\perp}$); promoted **`r_orth_cal` (`97.5% CI > 0`)** as primary correlation endpoint for **H-A2, H-B2, H-A3, H-B3, R4, and R14**, keeping `ls_slope_ols_top50 > 0.13` (`97.5% CI`) for H-A1/H-B1 and downgrading raw `r_cal` to secondary.
  - **R18**: Locked **`Transfer_Mean_GlobalScale`** control, mandatory abstract/summary global-amplitude artifact callout, and **100-permutation `NegCtrl_TMVS_ShuffledVar`** null (`median`, `2.5th`, `97.5th` pct; **H-B3** requires `tmvs_r_orth_cal_ci975_lo > max(all 5 mismatched var ci975_hi, perm100_p975)`).
  - **R19**: Added **descriptive 3-tier multiplicative gain variance decomposition** ($\beta_t$ whole-line gain $\to m_{t,g} - \beta_t$ per-gene gain increment $\to$ remaining off-diagonal residual + 95% CIs on both `resp_top50` and `all_nontarget`).
  - **R20 & `PRIOR_ART_B.md` §5.2**: Incorporated `A-T5`'s `ALERT_B` prior-art findings (`CGC`, `AMM-SimWMag`, VCC 2026 `CellularAI`/`AtlasShift`, `CASTER`, `Molina & Zhang 2026`), framing Group B `tmvs` strictly as *"a parameter-free variance-ratio variant of expression-aware transfer"* and framing Team A's contribution strictly around the 3-layer zero-shot residual decomposition + 3-tier multiplicative partition.
- **Binding Adjudication Rule**: **All confirmatory hypothesis verdicts (Group A `H-A1..H-A4` and Group B `H-B1..H-B4` at 97.5% CI), R4/R14 wording constraints, R16/R20 novelty limits, R18 callout, and independent audits follow Version 4 exclusively.**

---

## 1. Falsifiable Research Questions & Scientific Boundaries

In the Arc Virtual Cell Challenge 2026 Wrap-Up (2026-11-19), top-10 challenge models achieved high bulk/global agreement yet recovered a reported **cell-line-specific perturbation slope of only `0.02–0.13`** (versus **`1.02`** for experimental replicates). Existing cross-context predictors predominantly recover shared cross-line responses while missing cell-line-specific propagation.

**Boundaries with Prior Art (R7, R16, R20, `PRIOR_ART_B.md` §5.2)**:
- **CIPHER v2 (R7)** (`bioRxiv 10.1101/2025.06.27.661814 v2`, Kuznets-Speck et al. 2025/2026) already showed that removing off-diagonal gene–gene co-fluctuations drops performance several-fold and transfers across studies within the *same* cell type. **"Unperturbed control covariance carries perturbation information" is NOT claimed as a novel finding by Team A.**
- **Luan et al. (2026)** (`bioRxiv 10.64898/2026.10.20.760514`, 2026-10-25) evaluated centered transfer from `KOLF2.1J` to `H1` using target training perturbations rather than unperturbed target controls $C_t$.
- **`CGC` & `AMM-SimWMag` (R20 Tier-1 Prior Art)**: Tier-1 global amplitude orthogonalization $T - \beta_t S$ is already established in **`CGC`** (`Huang et al. 2026`, `bioRxiv 10.64898/2026.09.30.755756`) and **`AMM-SimWMag`** (`Sun et al. 2026`, `Genes 10.3390/genes17070816`).
- **Expression- / Variance-Aware Transfer (`ALERT_B`, R16/R20)**: Target control mean-expression scaling (`FCT` class) is already used by VCC 2026 solutions **`CellularAI`** (#1) and **`AtlasShift`** (#82), as well as supervised **`CASTER`** (`Kawsar et al. 2026`, `10.64898/2026.11.15.766815`, which found scaling transfer fails zero-shot). Consequently, Team A claims **zero methodological invention** for marginal-variance scaling and refers to **`Transfer_Mean_VarScaled` (`tmvs`) strictly as *"a parameter-free variance-ratio variant of expression-aware transfer"*** whose empirical zero-shot residual behavior is audited here.
- **Intended Methodological Contribution Scope (`PRIOR_ART_B.md` §5.2)**: If any signal survives confirmatory testing, Team A's contribution is framed strictly around **(a)** a 3-layer zero-shot control-statistics residual decomposition (`whole-line scale -> marginal gene variance -> off-diagonal co-fluctuation`) with `r_orth_cal` and `perm100` null controls, and **(b)** a descriptive 3-tier multiplicative gain variance accounting framework ($\beta_t \to m_{t,g} \to \text{remaining}$), explicitly contrasted against additive 2-way ANOVA (`Molina & Zhang 2026`, `10.64898/2026.07.24.740459`), which absorbs multiplicative whole-line and per-gene gains into additive interaction terms $\gamma_{c,p}$.

**Two Parallel Falsifiable Questions (Tested Strictly Zero-Shot with Zero Target Perturbations)**:
1. **Group A (`Hybrid_ProjShared`)**: Does adding target unperturbed single-cell control co-fluctuations $C_t$ onto shared cross-line transfer recover cell-line-specific residuals $\Delta_{p,t} - \bar{\Delta}^{\mathrm{src}}_{p,(-t)}$ with slope `> 0.13`, positive de-amplituded co-fluctuation `r_orth_cal > 0`, specificity over all 5 mismatched covariances on both slope and `r_orth_cal`, and `Transfer_Mean` non-inferiority?
2. **Group B (`Transfer_Mean_VarScaled`, Parameter-Free Variance-Ratio Transfer)**: Does zero-parameter target control gene-variance scaling $\bar{\Delta}^{\mathrm{src}}_{p,g} \sqrt{C_{t,gg}/\bar{C}_{S,gg}}$ (zero gene–gene covariance, zero training) recover cell-line-specific residuals with slope `> 0.13`, de-amplituded `r_orth_cal > 0`, specificity over all 5 mismatched source variances and 100 shuffled-gene permutations (`r_orth_cal`), and `Transfer_Mean` non-inferiority?

---

## 2. Rule 3 Zero-Look, Prior-Round Disclosure (R8), & Rule 5 Clean Data

1. **Zero-Look Prior to Remote Lock**: Prior to Master Controller verification in `round9/p2/status/PREREG_COMMIT.txt`, **no team member or script evaluated, tuned, inspected, or loaded expression matrices (`.npz`) from any held-out cell line** (`nadig_hepg2`, `nadig_jurkat`, `nourreddine_kolf21j`, `vcc2025_h1`). All dev verification (`PREREG_AUDIT_R3.md` and §3.2) used exclusively `Replogle K562 <-> RPE1`.
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

### 3.2 Locked Primary Configuration & Reproducible Dev Audit (`audit/dev/`, Rule 7 + R17–R19)

- **Locked Configuration**: `--space log1p`, `--shrink 0.0`, `cov_max_cells = all` (`4,976–20,000` per line), `mask_self = True`, `reproducible_only = True`; **Group A Primary**: `Hybrid_ProjShared`; **Group B Primary**: `Transfer_Mean_VarScaled`.
- **Pre-Registered Sensitivities**: **S1** `all_targets` (`only_reproducible = False`); **S2** `Hybrid_Raw` vs `Hybrid_ProjShared`; **S3** `shrink = 1e-3`; **S4** `cov_max_cells = 10,000` & `zscore`; **S5** expressed genes only (top 75% target `diag(C_t)`, R15); **S6** unselected `all_nontarget` split-half $\rho_{AB}, \rho_{AB}^{\perp}$ vs `resp_top50` (R15/R17).

**Exact Traceable `Replogle K562 <-> RPE1` Dev Numbers (`7,225` common genes, `832` reproducible shared targets)**  
*(Recorded directly in `round9/p2/vcc_linespecific/audit/dev/{audit_per_line.csv,audit_summary.json,AUDIT_REPORT.md}`)*:

- **`replogle_k562_essential <- replogle_rpe1`** ($\rho_{AB}^{\perp}\mathrm{(top50)} = 0.9273$, $c_{\mathrm{ceil}}^{\perp}\mathrm{(top50)} = 0.9810$; raw $\rho_{AB}\mathrm{(top50)} = 0.9503$, S6 $\rho_{AB}^{\perp}\mathrm{(all)} = 0.5418$, raw $\rho_{AB}\mathrm{(all)} = 0.8421$):
  - **R18 Global-Amplitude Artifact (`Transfer_Mean_GlobalScale`)**: Single scalar $\hat{c}_t$ gets raw `r_resid = +0.5546` (`r_cal = +0.5618`, `slope = +0.0308`), collapsing to **`r_orth_cal = 0.0000`** (`0.0008` on `all_nontarget`) once per-perturbation amplitude is removed.
  - **Group A `Hybrid_ProjShared`**: `slope = +0.0050` (`97.5% CI [0.0028, 0.0072]`), **`r_orth_cal = +0.0838` (`97.5% CI [0.0593, 0.1073]`)**, secondary raw `r_cal = +0.0523`, `r_cal_sq = 0.0027`, `sd_pred/sd_obs = 0.0964`, `slope_deming = +0.0050`. `SrcCorrTgtVar`: `r_orth_cal = +0.0235`, `slope = +0.0011`. S5 (`5,418` genes): `r_orth_cal = +0.0861`.
  - **Group B `Transfer_Mean_VarScaled`**: `slope = +0.0737` (`97.5% CI [0.0626, 0.0848]`), **`r_orth_cal = +0.1325` (`97.5% CI [0.1058, 0.1578]`)**, secondary raw `r_cal = +0.2198`, `r_cal_sq = 0.0483`, `sd_pred/sd_obs = 0.3396`, `slope_deming = +0.0827`, S5: `r_orth_cal = +0.1339`, S6 `r_orth_cal(all) = +0.1193`. **R18 100-Perm Shuffled Null (`H-B3` Fails on K562)**: `perm100 r_orth_cal median = +0.0902` `[p025 = +0.0605, p975 = +0.1207]` (`tmvs ci975_lo = +0.1058 < +0.1207`) because leaving unshuffled $\bar{C}_{S,gg}$ in denominator retains $\propto 1/\sqrt{\bar{C}_{S,gg}}$ scaling!
  - **R14 & FCT**: `Hybrid_on_TMVS` `r_orth_cal = +0.1500` (**`R14 delta = +0.0175`, `95% CI [+0.0118, +0.0233]`**); `FCT` `r_orth_cal = +0.2357`, `r_cal = +0.1394`, `slope = +0.0462`.
  - **R19 Multiplicative Gain Variance Split**: `resp_top50`: $\beta_t = 0.5611$ `[0.540, 0.581]`, `frac_whole = 0.3076` `[0.279, 0.338]`, `frac_per_gene_inc = 0.1425` `[0.130, 0.152]`, `frac_remaining = 0.5500` `[0.524, 0.579]`; `all_nontarget`: $\beta_t = 0.1989$, `0.6964 / 0.0259 / 0.2777`.
- **`replogle_rpe1 <- replogle_k562_essential`** ($\rho_{AB}^{\perp}\mathrm{(top50)} = 0.9447$, $c_{\mathrm{ceil}}^{\perp}\mathrm{(top50)} = 0.9857$; raw $\rho_{AB}\mathrm{(top50)} = 0.9441$, S6 $\rho_{AB}^{\perp}\mathrm{(all)} = 0.6965$, raw $\rho_{AB}\mathrm{(all)} = 0.7160$):
  - **R18 Global-Amplitude Artifact (`Transfer_Mean_GlobalScale`)**: Raw `r_resid = +0.3034` (`r_cal = +0.3078`, `slope = +0.0058`), collapsing to **`r_orth_cal = 0.0000`** (`0.0020` on `all_nontarget`).
  - **Group A `Hybrid_ProjShared`**: `slope = +0.0085` (`97.5% CI [0.0044, 0.0125]`), **`r_orth_cal = +0.1527` (`97.5% CI [0.1013, 0.2057]`)**, secondary raw `r_cal = +0.1054`, `r_cal_sq = 0.0111`, `sd_pred/sd_obs = 0.0823`, `slope_deming = +0.0086`. `SrcCorrTgtVar`: `r_orth_cal = +0.1211`, `slope = +0.0018`. S5 (`5,418` genes): `r_orth_cal = +0.1569`.
  - **Group B `Transfer_Mean_VarScaled`**: `slope = +0.0267` (`97.5% CI [0.0243, 0.0289]`), **`r_orth_cal = +0.0979` (`97.5% CI [0.0831, 0.1124]`)**, secondary raw `r_cal = +0.2145`, `r_cal_sq = 0.0460`, `sd_pred/sd_obs = 0.1262`, `slope_deming = +0.0271`, S5: `r_orth_cal = +0.0968`, S6 `r_orth_cal(all) = +0.1502`. **R18 100-Perm Shuffled Null (`H-B3` Passes on RPE1)**: `perm100 r_orth_cal median = +0.0320` `[p025 = -0.0053, p975 = +0.0734]` (`tmvs ci975_lo = +0.0831 > +0.0734`).
  - **R14 & FCT**: `Hybrid_on_TMVS` `r_orth_cal = +0.1664` (**`R14 delta = +0.0684`, `95% CI [+0.0437, +0.0924]`**); `FCT` `r_orth_cal = +0.0532`, `r_cal = +0.0903`, `slope = +0.0111`.
  - **R19 Multiplicative Gain Variance Split**: `resp_top50`: $\beta_t = 1.7340$ `[1.662, 1.812]`, `frac_whole = 0.0920` `[0.078, 0.107]`, `frac_per_gene_inc = 0.1195` `[0.106, 0.129]`, `frac_remaining = 0.7885` `[0.773, 0.808]`; `all_nontarget`: $\beta_t = 0.6230$, `0.0493 / 0.1103 / 0.8404`.
- Both Group A and Group B adjudicate **`未确认 (NOT CONFIRMED)`** on dev (`H-A1`/`H-B1` fail slope threshold; `H-A4` fails parity on dev; `H-B3` fails `K562`). Note that once global amplitude is orthogonalized (`r_orth_cal`), co-fluctuation structure shows consistent `R14` increment (`+0.0175` & `+0.0684`) on dev; whether this generalizes to the 3 held-out lines remains an open empirical test.

---

## 4. Exact LOO-CL Zero-Shot Predictors & Controls (R3, R4, R12, R14, R15, R18)

For each target cell line $t \in \{1,\dots,6\}$ with source lines $\mathcal{S} = \{1,\dots,6\} \setminus \{t\}$:

1. **Unperturbed Covariance Linear Response**: $v_\ell(p) = C_{\ell,:,g_p} / \max(C_{\ell,g_p,g_p}, 10^{-8})$, with $v_\ell(p)_{g_p} \leftarrow 0$.
2. **Signed Source-Only $\hat{\lambda}_{(-t)}$** and **Shared-Axis Projector $P_\perp = I_G - q_1 q_1^T$** fit strictly on $\mathcal{S}$:
   $$\hat{\lambda}_{(-t)} = \mathrm{clip}\!\left(\frac{\sum_{s \in \mathcal{S}} \sum_{p \in \mathcal{P}_s} \langle \Delta_{p,s},\, v_s(p) \rangle}{\max\!\left(\sum_{s \in \mathcal{S}} \sum_{p \in \mathcal{P}_s} \|v_s(p)\|_2^2,\, 10^{-12}\right)},\; -10.0,\; 10.0\right).$$
3. **Baselines & Group A Models** *(for $\mathcal{S}_p = \{s \in \mathcal{S} : p \text{ available & reproducible}\}$)*:
   - **`Transfer_Mean`**: $\bar{\Delta}^{\mathrm{src}}_{p,(-t)} = \frac{1}{|\mathcal{S}_p|}\sum_{s \in \mathcal{S}_p} \Delta_{p,s}$
   - **`Transfer_Mean_GlobalScale` *(R18 Single-Scalar Whole-Line Artifact Control)***:
     $$\hat{\Delta}^{\mathrm{GlobalScale}}_{p,t} = \hat{c}_t\, \bar{\Delta}^{\mathrm{src}}_{p,(-t)}, \qquad \hat{c}_t = \sqrt{\frac{\frac{1}{G}\sum_{g=1}^G \max(C_{t,gg},\, 10^{-8})}{\frac{1}{G}\sum_{g=1}^G \max(\bar{C}_{S,gg},\, 10^{-8})}},$$
     using strictly **1 scalar** per target line computed from unperturbed control variances.
   - **`CIPHER_Linear`**: $\hat{\lambda}_{(-t)} v_t(p)$
   - **`Hybrid_Raw`**: $\bar{\Delta}^{\mathrm{src}}_{p,(-t)} + \hat{\lambda}_{(-t)}(v_t(p) - \bar{v}_{\mathcal{S}_p}(p))$
   - **`Hybrid_ProjShared` *(Group A Primary)***: $\bar{\Delta}^{\mathrm{src}}_{p,(-t)} + \hat{\lambda}_{(-t)} P_\perp(v_t(p) - \bar{v}_{\mathcal{S}_p}(p))$
   - **`NegCtrl_HybridProj_MismatchedCov(s')` *(R3, all 5 source lines)***: Replace $v_t(p)$ with $v_{s'}(p)$ for each $s' \in \mathcal{S}$.
   - **`NegCtrl_HybridProj_SrcCorrTgtVar` *(R4a)***: Replace $C_t$ with $D_t^{1/2} \bar{R}_S D_t^{1/2}$ ($\bar{R}_S = \frac{1}{|\mathcal{S}|}\sum_s R_s$).
   - **Additional Covariance Nulls**: `NegCtrl_HybridProj_ShuffledGenes`, `NegCtrl_HybridProj_ShuffledAllCov`, `NegCtrl_HybridProj_DiagCov`.
4. **Group B Predictor & Specificity Controls (R12 & R18)**:
   - **`Transfer_Mean_VarScaled` *(Group B Primary, `tmvs` — parameter-free variance-ratio transfer)***:
     $$\hat{\Delta}^{\mathrm{TMVS}}_{p,t,g} = \bar{\Delta}^{\mathrm{src}}_{p,(-t),g} \cdot \sqrt{\frac{\max(C_{t,gg},\, 10^{-8})}{\max\!\left(\frac{1}{|\mathcal{S}|}\sum_{s \in \mathcal{S}} C_{s,gg},\, 10^{-8}\right)}}.$$
   - **`NegCtrl_TMVS_MismatchedVar(s')` *(all 5 source lines)***: Replace $C_{t,gg}$ in numerator with $C_{s',gg}$ for each $s' \in \mathcal{S}$.
   - **`NegCtrl_TMVS_ShuffledVar` (100-Permutation Null Distribution, R18)**: Generate **100 independent gene-label permutations** of the target variance diagonal vector $\max(\mathrm{diag}(C_t), 10^{-8})$ across the $G$ common genes while keeping source mean variances $\bar{C}_{S,gg}$ unshuffled, recording exact empirical **median, 2.5th, and 97.5th percentiles** of `r_orth_cal`, `r_cal`, and `slope_ols` (plus single-seed `tmvs_shuffled_var` diagnostic).
5. **Secondary Mechanism Predictor (R14)**:
   - **`Hybrid_on_TMVS`**:
     $$\hat{\Delta}^{\mathrm{Hybrid\_on\_TMVS}}_{p,t} = \hat{\Delta}^{\mathrm{TMVS}}_{p,t} + \hat{\lambda}_{(-t)}\, P_\perp\big(v_t(p) - \bar{v}_{\mathcal{S}_p}(p)\big),$$
     using the exact same $\hat{\lambda}_{(-t)}$ and $P_\perp$ fit on $\mathcal{S}$.
6. **Secondary Fold-Change Transfer Baseline (R15)**:
   - **`FCT`**: Using `ctrl_mean_libnorm` ($\mu_s, \mu_t$):
     $$\mathrm{FC}_{p,s} = \frac{\max\!\big(\mathrm{expm1}(\mathrm{log1p}(\max(\mu_s, 10^{-8})) + \Delta^{\mathrm{log1p}}_{p,s}),\, 10^{-8}\big)}{\max(\mu_s,\, 10^{-8})}, \qquad \overline{\mathrm{FC}}^{\mathrm{src}}_p = \exp\!\left(\frac{1}{|\mathcal{S}_p|}\sum_{s \in \mathcal{S}_p}\log\mathrm{FC}_{p,s}\right),$$
     $$\hat{\Delta}^{\mathrm{FCT}}_{p,t} = \mathrm{log1p}\!\big(\max(\mu_t, 0.0)\cdot \overline{\mathrm{FC}}^{\mathrm{src}}_p\big) - \mathrm{log1p}\!\big(\max(\mu_t, 0.0)\big).$$

---

## 5. Metrics, De-Amplituded `r_orth_cal`, Gain Decomposition, & Ceiling Limits (R5, R6, R11, R15, R17, R19)

1. **Line-Specific Residuals**: $r^{\mathrm{obs}}_{p,t} = \Delta^{\mathrm{obs}}_{p,t} - \bar{\Delta}^{\mathrm{src}}_{p,(-t)}$, $r^{\mathrm{pred}}_{p,t} = \hat{\Delta}_{p,t} - \bar{\Delta}^{\mathrm{src}}_{p,(-t)}$ evaluated on masked `resp_top50` (top 50 genes by $|\Delta^{\mathrm{obs}}_{p,g}|$, $g \ne g_p$) and unselected `all_nontarget`.
2. **Primary De-Amplituded Noise-Calibrated Correlation (`r_orth_cal`, R17)**:
   To eliminate whole-line/per-perturbation amplitude shrinkage along the shared transfer axis (on which 1-scalar `GlobalScale` alone gets `+0.555 / +0.303` raw correlation), let $b_p = \bar{\Delta}^{\mathrm{src}}_{p, M_p} / \max(\|\bar{\Delta}^{\mathrm{src}}_{p, M_p}\|_2, 10^{-12})$ on each perturbation's masked gene set $M_p$. Define amplitude-orthogonalized vectors:
   $$r^{\mathrm{obs},\perp}_p = r^{\mathrm{obs}}_{p, M_p} - \langle r^{\mathrm{obs}}_{p, M_p}, b_p \rangle\, b_p, \qquad r^{\mathrm{pred},\perp}_p = r^{\mathrm{pred}}_{p, M_p} - \langle r^{\mathrm{pred}}_{p, M_p}, b_p \rangle\, b_p,$$
   and analogously $(r^{(A),\perp}_p, r^{(B),\perp}_p)$ for disjoint 50/50 split-halves. Compute orthogonalized split-half ceiling $\rho_{AB}^{\perp} = \mathrm{corr}(r^{(A),\perp}, r^{(B),\perp})$, $c_{\mathrm{ceil}}^{\perp} = \sqrt{\frac{2\rho_{AB}^{\perp}}{1+\rho_{AB}^{\perp}}}$, and **$r_{\mathrm{orth\_cal}} = \mathrm{corr}(r^{\mathrm{pred},\perp}, r^{\mathrm{obs},\perp}) / c_{\mathrm{ceil}}^{\perp}$** (plus $r_{\mathrm{orth\_cal}}^2$), using **4,000** joint perturbation resamples ($\ge 95\%$ finite draws). Raw un-orthogonalized $r_{\mathrm{cal}}$ and $r_{\mathrm{cal}}^2$ (`R11`) remain reported as secondary diagnostics.
3. **Scope & Limitations of Slope Metric (R6)**: Neither Arc VCC 2026 Wrap-Up nor `cell-eval2 v0.22.0` publishes the `0.02–0.13` top-10 slope formula (`Dobin 10.64898/2026.11.17.767214` unindexed across Crossref/bioRxiv/OpenAlex/EuropePMC/PubMed APIs as of 2026-11-21). Moreover OLS `slope(r_pred ~ r_obs) = r * SD_pred / SD_obs`. Every line reports secondary `r`, `sd_pred/sd_obs`, orthogonal **Deming slope** ($\delta=1$), reverse OLS `slope_obs_on_pred`, zero-intercept slope, per-perturbation mean slope, and `cos_agreement_top50`.
4. **Descriptive 3-Tier Multiplicative Gain Variance Decomposition (R19)**:
   On both `resp_top50` and `all_nontarget`, decompose total cross-line residual variance $\mathrm{SST}(r^{\mathrm{obs}})$ ($\sum_p \sum_{g \in M_p} (r^{\mathrm{obs}}_{p,g} - \bar{r})^2$) into three exhaustive non-negative fractions (`sum = 1.0`) with perturbation-level 95% CIs (`C_distinct >= 5` perturbations per gene plus exact cluster-$W^2$ noise-overfit subtraction so finite-$df$ bootstrap replicates remain unbiased):
   - **Tier 1 (`frac_whole_line`)**: Variance captured by pooled whole-line multiplicative gain $\hat{\alpha} + \hat{\beta}_t \bar{\Delta}^{\mathrm{src}}$.
   - **Tier 2 (`frac_per_gene_inc`)**: Incremental variance captured by per-gene multiplicative slopes $\hat{m}_{t,g} - \hat{\beta}_t$ ($\hat{y}_{p,g} = \hat{\alpha} + \hat{m}_{t,g} \bar{\Delta}^{\mathrm{src}}_{p,g}$ for eligible genes).
   - **Tier 3 (`frac_remaining`)**: Unexplained cross-gene / perturbation-specific residual fraction.
5. **S5 Expressed-Genes Sensitivity (R15) & S6 Selection-Bias Disclosure**:
   - **S5**: Drop the lowest 25% of common genes ranked strictly by unperturbed target control variance `diag(C_t) <= Q_{0.25}`, re-select `resp_top50` within the remaining 75% expressed genes, and recompute `tmvs` and `Hybrid_ProjShared` `r_orth_cal`, `r_cal`, and slope.
   - **S6**: Selecting `resp_top50` by full-sample $|\Delta^{\mathrm{obs}}|$ upward-biases $\rho_{AB}$ and $\rho_{AB}^{\perp}$ (verified on `Replogle` dev: $\rho_{AB}^{\perp}$ drops `0.9273 -> 0.5418` on `K562`, `0.9447 -> 0.6965` on `RPE1`). Unselected `all_nontarget` ($\rho_{AB}^{\perp}, c_{\mathrm{ceil}}^{\perp}, r_{\mathrm{orth\_cal}}, r_{\mathrm{cal}}$) MUST always be reported alongside `resp_top50`.
6. **Global Accuracy Conventions (R5)**: Simultaneous non-inferiority across all **6** conventions: `ctrl_disc_ownmasked`, `ctrl_disc_panel`, `ctrl_de50_r`, `cent_disc_ownmasked`, `cent_disc_panel`, and `cent_de50_r` (with `ctrl_*` own-target in main tables and others in appendix).

---

## 6. Locked Confirmatory Groups A & B ($K_{\mathrm{conf}} = 3$, 97.5% CI) & Wording Governance

- **Confirmatory Cohort ($K_{\mathrm{conf}} = 3$)**: **`nadig_hepg2`**, **`nadig_jurkat`**, **`vcc2025_h1`**.  
- **Thresholds**: $N_{\mathrm{pass}} = \lceil 0.8 \times 3 \rceil = \mathbf{3/3}$ (`H-A1..A3`, `H-B1..B3`), $N_{\mathrm{inf}} = \lfloor 0.2 \times 3 \rfloor = \mathbf{0/3}$ (`H-A4`, `H-B4`).
- **R13 Multiple-Comparison Level**: **ALL confirmatory criteria below use two-sided 97.5% percentile confidence intervals (`[1.25%, 98.75%]`, Bonferroni over 2 groups)**. **Both Group A and Group B outcomes MUST be reported jointly in every abstract and summary, accompanied by the mandatory R18 `Transfer_Mean_GlobalScale` callout.**

### Confirmatory Group A (`Hybrid_ProjShared`)
1. **H-A1**: On 3/3 confirmatory lines, `Hybrid_ProjShared` `ls_slope_ols_top50` **97.5% CI lower bound $> 0.13$**. *(Pure `CIPHER_Linear` / `Hybrid_Raw` cannot substitute.)*
2. **H-A2 (R17)**: On 3/3 confirmatory lines, $\rho_{AB}^{\perp} > 0$, $\ge 95\%$ finite joint-bootstrap draws, and `Hybrid_ProjShared` de-amplituded **`r_orth_cal` 97.5% CI lower bound $> 0$**.
3. **H-A3 (R3 + R17)**: On 3/3 confirmatory lines, `Hybrid_ProjShared` beats **all 5** source mismatched covariance controls on **BOTH** **`r_orth_cal`** **AND** `ls_slope_ols_top50` (`matched_r_orth_cal_ci975_lo > max_{s'} mismatch_r_orth_cal_ci975_hi(s')` **AND** `matched_slope_ci975_lo > max_{s'} mismatch_slope_ci975_hi(s')`).
4. **H-A4 (R5)**: On 0/3 confirmatory lines, `Hybrid_ProjShared - Transfer_Mean` **97.5% CI upper bound $< 0$** across **ALL 6** conventions (`AND`).

### Confirmatory Group B (`Transfer_Mean_VarScaled`, Parameter-Free Variance-Ratio Transfer)
1. **H-B1**: On 3/3 confirmatory lines, `Transfer_Mean_VarScaled` `ls_slope_ols_top50` **97.5% CI lower bound $> 0.13$**.
2. **H-B2 (R17)**: On 3/3 confirmatory lines, $\rho_{AB}^{\perp} > 0$, $\ge 95\%$ finite draws, and `Transfer_Mean_VarScaled` de-amplituded **`r_orth_cal` 97.5% CI lower bound $> 0$**.
3. **H-B3 (R18 100-Perm & Mismatched Specificity on `r_orth_cal`)**: On 3/3 confirmatory lines, `Transfer_Mean_VarScaled` **`r_orth_cal` 97.5% CI lower bound** strictly exceeds `max` over **all 5** `NegCtrl_TMVS_MismatchedVar(s')` `r_orth_cal` 97.5% CI upper bounds **AND** the **100-permutation `NegCtrl_TMVS_ShuffledVar` `r_orth_cal` 97.5th percentile** (`tmvs_r_orth_cal_ci975_lo > max(max_mismatch_r_orth_cal_ci975_hi, perm100_r_orth_cal_p975)`).
4. **H-B4 (R5)**: On 0/3 confirmatory lines, `Transfer_Mean_VarScaled - Transfer_Mean` **97.5% CI upper bound $< 0$** across **ALL 6** conventions (`AND`).
5. **Descriptive Group B Check**: Report whether `tmvs` slope point estimate lies in `[0.02, 0.13]` with R6 estimator/scale caveat.

### Hard-Coded Decision, Mechanism Wording, GlobalScale Callout, & Novelty Governance (R4, R14, R16, R18, R20)
- **Group Verdicts**: Each group (`Group A`, `Group B`) independently yields **`CONFIRMED`** iff all 4 of its criteria (`1 AND 2 AND 3 AND 4`) pass at 97.5% CI; if any criterion in a group fails, that group **MUST state `未确认` (NOT CONFIRMED)** verbatim.
- **R4 Locks (Group A, Updated to `r_orth_cal` per R17)**:
  1. Use **`"cell-line-specific co-fluctuation structure"`** only if `Hybrid_ProjShared` beats `SrcCorrTgtVar` on BOTH slope and `r_orth_cal` (disjoint CIs, 3/3 confirmatory lines); otherwise use **`"target-line control statistics"`**.
  2. Unless `Hybrid_ProjShared` beats `Transfer_Mean_VarScaled` on `r_orth_cal` (3/3 confirmatory lines), state explicitly that **gains over `Transfer_Mean` may be explained by target-line marginal gene variance**.
- **R14 Mechanism Lock (Updated to `r_orth_cal` per R17)**: Only if `Hybrid_on_TMVS - Transfer_Mean_VarScaled` **`r_orth_cal`** paired bootstrap 95% CI lower bound is **`> 0` on 3/3 confirmatory lines** may reports state *"co-fluctuation structure provides incremental information beyond marginal variance"*; otherwise they **MUST** state *"In this dataset, co-fluctuation structure exhibits no detectable increment beyond marginal variance"*.
- **R18 Mandatory Abstract / Summary Callout**: Every report MUST state explicitly: *"Using a single scalar per cell line (`Transfer_Mean_GlobalScale`) achieves raw residual correlation X (and collapses near 0 under de-amplituded `r_orth_cal`)."*
- **R16/R20 Novelty Lock**: Per `PRIOR_ART_B.md` §5.2 (`ALERT_B`), never claim methodological novelty for `tmvs`, `FCT`, or Tier-1 global amplitude orthogonalization (`CGC`, `AMM-SimWMag`, `CellularAI`, `AtlasShift`, `CASTER`); describe `tmvs` strictly as *"a parameter-free variance-ratio variant of expression-aware transfer"*.

---

## 7. Prospective VCC 2026 Protocol (`AC16`, `MOLM13`, `U-87MG`, R9)

On **`2026-11-21 19:28 (+08:00)`**, Human lead (`Zijian Zeng`) emailed **`help@virtualcellchallenge.org`** requesting permission for offline `controls-final` research predictions (`prospective_predict.py`) and pre-GT hash commits (`CELL_EVAL2_REPRO.md` / `HUMAN_STEPS.md`). Agents never download gated data, complete KYC, or push externally. Upon compliant `controls-final` access and post-CC0 GT release, evaluate unmodified predictions via `cell-eval2 v0.22.0` (`--preset vcc2026`) under the exact Version-4 **Group A + Group B (`97.5%` CI, `r_orth_cal`, `perm100`) + R4/R14/R16/R18–R20** rules, reporting any failed group verbatim as `未确认`.

---

## 8. Independent Statistical Audit Plan (`A-T3`)

Strictly AFTER `round9/p2/status/PREREG_COMMIT.txt` exists and `A-T2` completes 6-line held-out + sensitivity runs in `round9/p2/vcc_linespecific/results/`:
1. Run standalone `audit/audit_runner.py` (**zero `src/` imports**) with locked primary parameters (`--space log1p --shrink 0.0 --cov-max-cells 0 --project-shared-axis --confirmatory-lines nadig_hepg2 nadig_jurkat vcc2025_h1`), auditing all 6 lines, Group A + Group B at 97.5% CI on `r_orth_cal` + slope, R18 `GlobalScale` + 100-perm `NegCtrl_TMVS_ShuffledVar`, R19 3-tier variance decomposition, R4/R14 wording checks, R15 `FCT` + `S5` + `S6`, and S1–S4 sensitivity runs.
2. Verify zero-leakage invariance and numerical parity (`tol <= 1e-4`) against `A-T2` via `audit/cross_check.py`, publishing the unvarnished final report in `round9/p2/vcc_linespecific/audit/AUDIT_REPORT.md`.

---

## 9. API-Verified Bibliography

All citations verified directly via Crossref, bioRxiv, or arXiv APIs on 2026-11-21 *(R16/R20 Group B screening recorded by `A-T5` in `audit/LIT/PRIOR_ART_B.md`)*:
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
11. **Sun B, Zheng H, Chen J, Jin J** (2026). *A Neural-Network-Free Calibration Matches or Beats Deep Single-Cell Perturbation Response Models Across Four Datasets* (`AMM-SimWMag`). *Genes*, `doi:10.3390/genes17070816` (2026-07-17).
12. **Dhaka B et al. (3 authors)** (2026). *HyperMap: An Efficient Framework for Transferring Perturbation Responses Across Diverse Biological Contexts*. bioRxiv `doi:10.64898/2026.04.23.720505` (2026-04-27).
13. **Luan et al. (4 authors)** (2026). *Centered cross-context perturbation transfer on the public 2025 Virtual Cell Challenge benchmark*. bioRxiv `doi:10.64898/2026.10.20.760514` (2026-10-25).
14. **Molina et al. (2 authors)** (2026). *Perturbation response decomposition enables biologically aligned generalization to unseen perturbations and cellular contexts*. bioRxiv `doi:10.64898/2026.07.24.740459` (2026-07-27).
15. **Gao J** (2026). *The Delta Gap: Cellular World Models for Perturbation Response Prediction Are Far Harder Than Global Correlation Implies*. `arXiv:2611.10658v1` (2026-11-07).
16. **Zhang Z, Huang Y, Feng F et al. (4 authors)** (2026). *A low-rank decoder bottleneck bounds reliability in foundation-model perturbation prediction*. `arXiv:2611.05849v1` (2026-11-04).
17. **Arc Institute (`cell-eval2` & VCC 2026 Wrap-Up)**: `ArcInstitute/cell-eval2` `v0.22.0` (`--preset vcc2026`) and Arc VCC 2026 Wrap-Up (`2026-11-19`, including VCC 2026 `#1 CellularAI` and `#82 AtlasShift` solutions). *Unindexed DOI Disclosure*: Dobin metrics preprint (`10.64898/2026.11.17.767214`) returned HTTP 404 / 0 hits via Crossref, bioRxiv, OpenAlex, EuropePMC, and PubMed APIs as of 2026-11-21 (`www.biorxiv.org` HTML Cloudflare-gated, not bypassed).
