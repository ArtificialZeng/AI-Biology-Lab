# Round 11 Confirmatory Group D (`FCT`, Primary Mask `all_nontarget`) Report (`dev_r10_r9_6_log1p_smoke_v2_r10`)

- **Phase**: `dev_r10` | **Source Mode**: `r9_6` (`6` sources) | **PREREG_D Commit**: `` | **Code Freeze**: ``
- **Bootstrap**: `np.random.default_rng(20261123)` x `4000` (Bonferroni 2-claim + IUT two-sided `97.5%` CIs `[1.25%, 98.75%]`)
- **Cohort Disclosure**: **Insufficient confirmatory cell lines (K_conf < 2); confirmatory hypothesis testing was not performed (确认性线不足，未做确认性检验).** (`[]`)
- **Secondary Lines**: `['dld1']`
- **Confirmatory Claim 1 — Primary Group D (`FCT` 全基因线特异信号) Verdict**: **NOT_EVALUATED_INSUFFICIENT_LINES**
- **Confirmatory Claim 2 — Methodological Control (`H-D4` `resp_top50` Selection Artifact)**: **NOT_EVALUATED_INSUFFICIENT_LINES** (all lines descriptive: `CONFIRMED`)
- **Parallel Practical Utility (`H-D_util` / `H-C3`-style Non-Inferiority)**: **NOT_EVALUATED_INSUFFICIENT_LINES**

## Per-Line Primary Mask (`all_nontarget`, Joint-Bootstrap 97.5% CIs) & Secondary (`resp_top50`)

| Target Line | Role | Method | N | Primary `r_orth_cal` | Primary 97.5% CI | `resp_top50` `r_orth_cal` | `resp_top50` 97.5% CI | `rho_AB_orth` |
|---|---|---|---:|---:|---|---:|---|---:|
| dld1 | Secondary | `Source_Grand_Mean` | 86 | +0.1788 | [+0.1088, +0.2473] | +0.2222 | [+0.1134, +0.3355] | +0.2533 |
| dld1 | Secondary | `Transfer_Mean` | 86 | +nan | [+nan, +nan] | +nan | [+nan, +nan] | +0.2533 |
| dld1 | Secondary | `FCT` | 86 | +0.1355 | [+0.1079, +0.1620] | +0.1401 | [+0.0850, +0.1902] | +0.2533 |
| dld1 | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 86 | +0.0389 | [+0.0242, +0.0534] | +0.0294 | [-0.0163, +0.0707] | +0.2533 |
| dld1 | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 86 | +0.0063 | [-0.0130, +0.0256] | -0.0225 | [-0.0718, +0.0276] | +0.2533 |
| dld1 | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 86 | +0.0215 | [+0.0076, +0.0360] | +0.0300 | [-0.0135, +0.0729] | +0.2533 |
| dld1 | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 86 | +0.0197 | [+0.0039, +0.0363] | -0.0058 | [-0.0570, +0.0462] | +0.2533 |
| dld1 | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 86 | +0.0553 | [+0.0323, +0.0781] | +0.0463 | [-0.0175, +0.1027] | +0.2533 |
| dld1 | Secondary | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 86 | -0.0379 | [-0.0492, -0.0273] | +0.0707 | [+0.0330, +0.1065] | +0.2533 |
| dld1 | Secondary | `NegCtrl_FCT_Perm0` | 86 | -0.0416 | [-0.0534, -0.0310] | +0.1159 | [+0.0777, +0.1499] | +0.2533 |
| dld1 | Secondary | `FCT_FullNTMean` | 86 | +0.1355 | [+0.1079, +0.1620] | +0.1401 | [+0.0850, +0.1902] | +0.2533 |
| dld1 | Secondary | `Transfer_Mean_GlobalScale` | 86 | +nan | [+nan, +nan] | +nan | [+nan, +nan] | +0.2533 |
| dld1 | Secondary | `Transfer_Mean_VarScaled` | 86 | +0.1389 | [+0.1198, +0.1579] | +0.0888 | [+0.0236, +0.1549] | +0.2533 |
| dld1 | Secondary | `CIPHER_Linear` | 86 | +0.0235 | [-0.0326, +0.0785] | -0.0314 | [-0.1054, +0.0501] | +0.2533 |
| dld1 | Secondary | `Hybrid_ProjShared` | 86 | +0.0138 | [-0.0192, +0.0459] | -0.0429 | [-0.1517, +0.0669] | +0.2533 |
| dld1 | Secondary | `Hybrid_on_TMVS` | 86 | +0.0955 | [+0.0665, +0.1249] | +0.0373 | [-0.0657, +0.1367] | +0.2533 |
| dld1 | Secondary | `E1_ConstMed` | 86 | -0.0614 | [-0.0774, -0.0467] | +0.1567 | [+0.1076, +0.2033] | +0.2533 |
| dld1 | Secondary | `E2_ConstP10` | 86 | -0.0573 | [-0.0733, -0.0424] | +0.1364 | [+0.0827, +0.1875] | +0.2533 |
| dld1 | Secondary | `E2_ConstP90` | 86 | -0.0673 | [-0.0837, -0.0521] | +0.1719 | [+0.1265, +0.2165] | +0.2533 |
| dld1 | Secondary | `E3_SrcMeanCtrl` | 86 | +0.0069 | [-0.0102, +0.0239] | +0.0755 | [+0.0174, +0.1290] | +0.2533 |
| dld1 | Secondary | `E5_RankMatchedMu` | 86 | +0.0372 | [+0.0188, +0.0554] | +0.0418 | [-0.0182, +0.0959] | +0.2533 |
| dld1 | Secondary | `E6_LogFCBar` | 86 | -0.0713 | [-0.0881, -0.0554] | +0.1492 | [+0.1078, +0.1914] | +0.2533 |
| dld1 | Secondary | `E7_FCBarMinus1` | 86 | -0.0541 | [-0.0702, -0.0392] | +0.1165 | [+0.0578, +0.1737] | +0.2533 |

## `H-D3` (`FCT − E1`) & `H-D4` (`E1(resp_top50) − E1(all_nontarget)`) Paired Comparisons

| Target Line | Comparison | Estimate | 95% CI | 97.5% CI |
|---|---|---:|---|---|
| dld1 | `H-D3_FCT_minus_E1_all_nontarget` | +0.1969 | [+0.1709, +0.2219] | [+0.1672, +0.2251] |
| dld1 | `diag_FCT_minus_E1_resp_top50` | -0.0166 | [-0.0660, +0.0313] | [-0.0718, +0.0385] |
| dld1 | `H-D4_E1_resp_top50_minus_all_nontarget` | +0.2180 | [+0.1705, +0.2634] | [+0.1639, +0.2697] |
| dld1 | `H-D5_FCT_minus_Mis_kolf21j_all_nontarget` | +0.1140 | [+0.0962, +0.1311] | [+0.0942, +0.1333] |

## `SA-8` Split-Half Out-of-Sample Top-50 (`resp_top50_split`, Single-Half Ceiling `sqrt(rho_AB)`)

| Target Line | Method | N | `r_orth_raw_split` | `r_orth_cal_split` | 97.5% CI | `Perm100 split` Med / p97.5 |
|---|---|---:|---:|---:|---|---|
| dld1 | `FCT` | 86 | +0.1145 | +0.1407 | [+0.0722, +0.2014] | +0.0689 / +0.1008 |
| dld1 | `E1_ConstMed` | 86 | +0.1053 | +0.1294 | [+0.0887, +0.1692] | — |

## 100-Perm Null (`all_nontarget` Primary, `resp_top50`, & `SA-8 resp_top50_split`)

- **dld1**: `all_nontarget` median=-0.0428, p97.5=-0.0352; `resp_top50` median=+0.0892, p97.5=+0.1179; `resp_top50_split` median=+0.0689, p97.5=+0.1008
