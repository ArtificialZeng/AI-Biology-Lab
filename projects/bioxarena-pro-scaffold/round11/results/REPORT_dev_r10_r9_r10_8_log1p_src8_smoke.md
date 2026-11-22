# Round 11 Confirmatory Group D (`FCT`, Primary Mask `all_nontarget`) Report (`dev_r10_r9_r10_8_log1p_src8_smoke`)

- **Phase**: `dev_r10` | **Source Mode**: `r9_r10_8` (`8` sources) | **PREREG_D Commit**: `` | **Code Freeze**: ``
- **Bootstrap**: `np.random.default_rng(20261123)` x `400` (Bonferroni co-primary `98.33%` CIs `[0.8333%, 99.1667%]`)
- **Cohort Disclosure**: **Insufficient confirmatory cell lines (K_conf < 2); confirmatory hypothesis testing was not performed (确认性线不足，未做确认性检验).** (`[]`)
- **Secondary Lines**: `['calu3']`
- **Primary Group D (`FCT` 全基因线特异信号) Verdict**: **NOT_EVALUATED_INSUFFICIENT_LINES**
- **Secondary Methodological Hypothesis (`H-D4`)**: **NOT_EVALUATED_INSUFFICIENT_LINES** (all lines descriptive: `未确认`)
- **Parallel Practical Utility (`H-D_util` / `H-C3`-style Non-Inferiority)**: **NOT_EVALUATED_INSUFFICIENT_LINES**

## Per-Line Primary Mask (`all_nontarget`, Joint-Bootstrap 98.33% CIs) & Secondary (`resp_top50`)

| Target Line | Role | Method | N | Primary `r_orth_cal` | Primary 98.33% CI | `resp_top50` `r_orth_cal` | `resp_top50` 98.33% CI | `rho_AB_orth` |
|---|---|---|---:|---:|---|---:|---|---:|
| calu3 | Secondary | `Source_Grand_Mean` | 9 | -0.0004 | [-0.3203, +0.2978] | -0.0988 | [-0.3934, +0.2751] | +0.3621 |
| calu3 | Secondary | `Transfer_Mean` | 9 | +nan | [+nan, +nan] | +nan | [+nan, +nan] | +0.3621 |
| calu3 | Secondary | `FCT` | 9 | +0.0597 | [-0.0085, +0.1302] | +0.0169 | [-0.1070, +0.1776] | +0.3621 |
| calu3 | Secondary | `NegCtrl_FCT_MismatchedMean(cd4t_rest)` | 9 | +0.0286 | [-0.0238, +0.0834] | +0.0504 | [-0.0641, +0.1320] | +0.3621 |
| calu3 | Secondary | `NegCtrl_FCT_MismatchedMean(dld1)` | 9 | -0.0082 | [-0.0517, +0.0502] | -0.0995 | [-0.1854, +0.0109] | +0.3621 |
| calu3 | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 9 | -0.0029 | [-0.0438, +0.0282] | -0.0106 | [-0.1032, +0.0628] | +0.3621 |
| calu3 | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 9 | -0.0441 | [-0.0800, -0.0003] | -0.1046 | [-0.1929, -0.0142] | +0.3621 |
| calu3 | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 9 | -0.0238 | [-0.0549, +0.0159] | -0.1140 | [-0.1786, -0.0250] | +0.3621 |
| calu3 | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 9 | -0.0424 | [-0.0668, -0.0044] | -0.0802 | [-0.1455, +0.0165] | +0.3621 |
| calu3 | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 9 | +0.0012 | [-0.0607, +0.0648] | -0.0922 | [-0.1972, +0.0800] | +0.3621 |
| calu3 | Secondary | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 9 | -0.0040 | [-0.0402, +0.0439] | +0.0170 | [-0.1319, +0.1220] | +0.3621 |
| calu3 | Secondary | `NegCtrl_FCT_Perm0` | 9 | +0.0122 | [-0.0321, +0.0587] | +0.0577 | [-0.0344, +0.1769] | +0.3621 |
| calu3 | Secondary | `FCT_FullNTMean` | 9 | +0.0597 | [-0.0085, +0.1302] | +0.0169 | [-0.1070, +0.1776] | +0.3621 |
| calu3 | Secondary | `Transfer_Mean_GlobalScale` | 9 | +nan | [+nan, +nan] | +nan | [+nan, +nan] | +0.3621 |
| calu3 | Secondary | `Transfer_Mean_VarScaled` | 9 | +0.0675 | [-0.0264, +0.1328] | -0.0104 | [-0.1074, +0.0895] | +0.3621 |
| calu3 | Secondary | `CIPHER_Linear` | 9 | +0.3371 | [-0.2738, +0.6839] | +0.3603 | [-0.2890, +0.7976] | +0.3621 |
| calu3 | Secondary | `Hybrid_ProjShared` | 9 | +0.2501 | [-0.2714, +0.5159] | +0.3130 | [-0.3097, +0.7166] | +0.3621 |
| calu3 | Secondary | `Hybrid_on_TMVS` | 9 | +0.2447 | [-0.2497, +0.5041] | +0.2436 | [-0.2928, +0.6340] | +0.3621 |
| calu3 | Secondary | `E1_ConstMed` | 9 | +0.0210 | [-0.0271, +0.0812] | +0.0881 | [-0.1092, +0.2731] | +0.3621 |
| calu3 | Secondary | `E2_ConstP10` | 9 | +0.0245 | [-0.0203, +0.0810] | +0.0842 | [-0.1167, +0.2628] | +0.3621 |
| calu3 | Secondary | `E2_ConstP90` | 9 | +0.0147 | [-0.0394, +0.0822] | +0.0857 | [-0.0975, +0.2730] | +0.3621 |
| calu3 | Secondary | `E3_SrcMeanCtrl` | 9 | -0.0072 | [-0.0406, +0.0480] | -0.0737 | [-0.1840, +0.0730] | +0.3621 |
| calu3 | Secondary | `E5_RankMatchedMu` | 9 | -0.0221 | [-0.0651, +0.0163] | -0.1259 | [-0.2228, +0.0121] | +0.3621 |
| calu3 | Secondary | `E6_LogFCBar` | 9 | +0.0100 | [-0.0483, +0.0825] | +0.0728 | [-0.0831, +0.2544] | +0.3621 |
| calu3 | Secondary | `E7_FCBarMinus1` | 9 | +0.0264 | [-0.0169, +0.0814] | +0.0803 | [-0.1226, +0.2534] | +0.3621 |

## `H-D3` (`FCT − E1`) & `H-D4` (`E1(resp_top50) − E1(all_nontarget)`) Paired Comparisons

| Target Line | Comparison | Estimate | 95% CI | 98.33% CI |
|---|---|---:|---|---|
| calu3 | `H-D3_FCT_minus_E1_all_nontarget` | +0.0387 | [-0.0391, +0.0958] | [-0.0618, +0.1029] |
| calu3 | `diag_FCT_minus_E1_resp_top50` | -0.0713 | [-0.1366, +0.0397] | [-0.1524, +0.0779] |
| calu3 | `H-D4_E1_resp_top50_minus_all_nontarget` | +0.0672 | [-0.1079, +0.1969] | [-0.1198, +0.2562] |

## 100-Perm Null (`all_nontarget` Primary & `resp_top50` Secondary)

- **calu3**: `all_nontarget` median=+0.0104, p97.5=+0.0190; `resp_top50` median=+0.0476, p97.5=+0.0939
