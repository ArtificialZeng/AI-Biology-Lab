# Round 11 Confirmatory Group D (`FCT`, Primary Mask `all_nontarget`) Report (`dev_r9_r9_6_log1psmoke_v2_r9`)

- **Phase**: `dev_r9` | **Source Mode**: `r9_6` (`6` sources) | **PREREG_D Commit**: `` | **Code Freeze**: ``
- **Bootstrap**: `np.random.default_rng(20261123)` x `400` (Bonferroni 2-claim + IUT two-sided `97.5%` CIs `[1.25%, 98.75%]`)
- **Cohort Disclosure**: **Insufficient confirmatory cell lines (K_conf < 2); confirmatory hypothesis testing was not performed (确认性线不足，未做确认性检验).** (`[]`)
- **Secondary Lines**: `['vcc2025_h1']`
- **Confirmatory Claim 1 — Primary Group D (`FCT` 全基因线特异信号) Verdict**: **NOT_EVALUATED_INSUFFICIENT_LINES**
- **Confirmatory Claim 2 — Methodological Control (`H-D4` `resp_top50` Selection Artifact)**: **NOT_EVALUATED_INSUFFICIENT_LINES** (all lines descriptive: `CONFIRMED`)
- **Parallel Practical Utility (`H-D_util` / `H-C3`-style Non-Inferiority)**: **NOT_EVALUATED_INSUFFICIENT_LINES**

## Per-Line Primary Mask (`all_nontarget`, Joint-Bootstrap 97.5% CIs) & Secondary (`resp_top50`)

| Target Line | Role | Method | N | Primary `r_orth_cal` | Primary 97.5% CI | `resp_top50` `r_orth_cal` | `resp_top50` 97.5% CI | `rho_AB_orth` |
|---|---|---|---:|---:|---|---:|---|---:|
| vcc2025_h1 | Secondary | `Source_Grand_Mean` | 80 | +0.0446 | [+0.0165, +0.0777] | +0.1086 | [+0.0486, +0.1694] | +0.8886 |
| vcc2025_h1 | Secondary | `Transfer_Mean` | 80 | +nan | [+nan, +nan] | +nan | [+nan, +nan] | +0.8886 |
| vcc2025_h1 | Secondary | `FCT` | 80 | +0.1427 | [+0.1189, +0.1726] | +0.3505 | [+0.2876, +0.4248] | +0.8886 |
| vcc2025_h1 | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 80 | +0.0451 | [+0.0343, +0.0557] | +0.1240 | [+0.0811, +0.1676] | +0.8886 |
| vcc2025_h1 | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 80 | +0.0180 | [+0.0085, +0.0284] | +0.0550 | [+0.0077, +0.0978] | +0.8886 |
| vcc2025_h1 | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 80 | +0.0537 | [+0.0429, +0.0654] | +0.1453 | [+0.0960, +0.1977] | +0.8886 |
| vcc2025_h1 | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 80 | +0.0281 | [+0.0183, +0.0378] | +0.0536 | [+0.0157, +0.0934] | +0.8886 |
| vcc2025_h1 | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 80 | -0.0134 | [-0.0290, -0.0006] | -0.0541 | [-0.1064, -0.0042] | +0.8886 |
| vcc2025_h1 | Secondary | `NegCtrl_FCT_Perm0` | 80 | +0.0318 | [+0.0204, +0.0448] | +0.1802 | [+0.1319, +0.2216] | +0.8886 |
| vcc2025_h1 | Secondary | `Transfer_Mean_GlobalScale` | 80 | +nan | [+nan, +nan] | +nan | [+nan, +nan] | +0.8886 |
| vcc2025_h1 | Secondary | `Transfer_Mean_VarScaled` | 80 | +0.0323 | [-0.0044, +0.0710] | -0.0772 | [-0.1863, +0.0317] | +0.8886 |
| vcc2025_h1 | Secondary | `CIPHER_Linear` | 80 | +0.0022 | [-0.0243, +0.0293] | -0.0468 | [-0.1501, +0.0379] | +0.8886 |
| vcc2025_h1 | Secondary | `Hybrid_ProjShared` | 80 | +0.0067 | [-0.0155, +0.0269] | -0.0403 | [-0.1329, +0.0344] | +0.8886 |
| vcc2025_h1 | Secondary | `Hybrid_on_TMVS` | 80 | +0.0208 | [-0.0055, +0.0428] | -0.0819 | [-0.1870, +0.0054] | +0.8886 |
| vcc2025_h1 | Secondary | `E1_ConstMed` | 80 | +0.0437 | [+0.0288, +0.0619] | +0.3076 | [+0.2473, +0.3802] | +0.8886 |
| vcc2025_h1 | Secondary | `E2_ConstP10` | 80 | +0.0435 | [+0.0282, +0.0636] | +0.2796 | [+0.2154, +0.3544] | +0.8886 |
| vcc2025_h1 | Secondary | `E2_ConstP90` | 80 | +0.0437 | [+0.0291, +0.0606] | +0.3189 | [+0.2660, +0.3887] | +0.8886 |
| vcc2025_h1 | Secondary | `E3_SrcMeanCtrl` | 80 | +0.0429 | [+0.0327, +0.0526] | +0.1182 | [+0.0773, +0.1644] | +0.8886 |
| vcc2025_h1 | Secondary | `E5_RankMatchedMu` | 80 | +0.0592 | [+0.0450, +0.0744] | +0.2108 | [+0.1588, +0.2682] | +0.8886 |
| vcc2025_h1 | Secondary | `E6_LogFCBar` | 80 | +0.0419 | [+0.0283, +0.0601] | +0.3117 | [+0.2582, +0.3753] | +0.8886 |
| vcc2025_h1 | Secondary | `E7_FCBarMinus1` | 80 | +0.0430 | [+0.0282, +0.0637] | +0.2559 | [+0.1881, +0.3304] | +0.8886 |

## `H-D3` (`FCT − E1`) & `H-D4` (`E1(resp_top50) − E1(all_nontarget)`) Paired Comparisons

| Target Line | Comparison | Estimate | 95% CI | 97.5% CI |
|---|---|---:|---|---|
| vcc2025_h1 | `H-D3_FCT_minus_E1_all_nontarget` | +0.0989 | [+0.0852, +0.1143] | [+0.0836, +0.1174] |
| vcc2025_h1 | `diag_FCT_minus_E1_resp_top50` | +0.0429 | [+0.0157, +0.0728] | [+0.0094, +0.0760] |
| vcc2025_h1 | `H-D4_E1_resp_top50_minus_all_nontarget` | +0.2639 | [+0.2164, +0.3147] | [+0.2120, +0.3217] |
| vcc2025_h1 | `H-D5_FCT_minus_Mis_kolf21j_all_nontarget` | +0.0890 | [+0.0735, +0.1064] | [+0.0731, +0.1106] |

## `SA-8` Split-Half Out-of-Sample Top-50 (`resp_top50_split`, Single-Half Ceiling `sqrt(rho_AB)`)

| Target Line | Method | N | `r_orth_raw_split` | `r_orth_cal_split` | 97.5% CI | `Perm100 split` Med / p97.5 |
|---|---|---:|---:|---:|---|---|
| vcc2025_h1 | `FCT` | 80 | +0.3510 | +0.3533 | [+0.2889, +0.4263] | +0.2120 / +0.2423 |
| vcc2025_h1 | `E1_ConstMed` | 80 | +0.3079 | +0.3100 | [+0.2509, +0.3852] | — |

## 100-Perm Null (`all_nontarget` Primary, `resp_top50`, & `SA-8 resp_top50_split`)

- **vcc2025_h1**: `all_nontarget` median=+0.0319, p97.5=+0.0363; `resp_top50` median=+0.2088, p97.5=+0.2381; `resp_top50_split` median=+0.2120, p97.5=+0.2423
