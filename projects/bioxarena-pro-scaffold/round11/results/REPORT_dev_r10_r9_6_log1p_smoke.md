# Round 11 Confirmatory Group D (`FCT`, Primary Mask `all_nontarget`) Report (`dev_r10_r9_6_log1p_smoke`)

- **Phase**: `dev_r10` | **Source Mode**: `r9_6` (`6` sources) | **PREREG_D Commit**: `` | **Code Freeze**: ``
- **Bootstrap**: `np.random.default_rng(20261123)` x `4000` (Bonferroni co-primary `98.33%` CIs `[0.8333%, 99.1667%]`)
- **Cohort Disclosure**: **Insufficient confirmatory cell lines (K_conf < 2); confirmatory hypothesis testing was not performed (确认性线不足，未做确认性检验).** (`[]`)
- **Secondary Lines**: `['cd4t_rest', 'dld1']`
- **Primary Group D (`FCT` 全基因线特异信号) Verdict**: **NOT_EVALUATED_INSUFFICIENT_LINES**
- **Secondary Methodological Hypothesis (`H-D4`)**: **NOT_EVALUATED_INSUFFICIENT_LINES** (all lines descriptive: `CONFIRMED`)
- **Parallel Practical Utility (`H-C3`-style Non-Inferiority)**: **NOT_EVALUATED_INSUFFICIENT_LINES**

## Per-Line Primary Mask (`all_nontarget`, Joint-Bootstrap 98.33% CIs) & Secondary (`resp_top50`)

| Target Line | Role | Method | N | Primary `r_orth_cal` | Primary 98.33% CI | `resp_top50` `r_orth_cal` | `resp_top50` 98.33% CI | `rho_AB_orth` |
|---|---|---|---:|---:|---|---:|---|---:|
| dld1 | Secondary | `Source_Grand_Mean` | 86 | +0.1788 | [+0.1053, +0.2535] | +0.2222 | [+0.1092, +0.3413] | +0.2533 |
| dld1 | Secondary | `Transfer_Mean` | 86 | +nan | [+nan, +nan] | +nan | [+nan, +nan] | +0.2533 |
| dld1 | Secondary | `FCT` | 86 | +0.1355 | [+0.1059, +0.1645] | +0.1401 | [+0.0817, +0.1933] | +0.2533 |
| dld1 | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 86 | +0.0389 | [+0.0235, +0.0546] | +0.0294 | [-0.0188, +0.0735] | +0.2533 |
| dld1 | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 86 | +0.0063 | [-0.0143, +0.0270] | -0.0225 | [-0.0740, +0.0307] | +0.2533 |
| dld1 | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 86 | +0.0215 | [+0.0067, +0.0371] | +0.0300 | [-0.0169, +0.0748] | +0.2533 |
| dld1 | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 86 | +0.0197 | [+0.0022, +0.0371] | -0.0058 | [-0.0616, +0.0471] | +0.2533 |
| dld1 | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 86 | +0.0553 | [+0.0309, +0.0796] | +0.0463 | [-0.0218, +0.1050] | +0.2533 |
| dld1 | Secondary | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 86 | -0.0379 | [-0.0501, -0.0267] | +0.0707 | [+0.0291, +0.1085] | +0.2533 |
| dld1 | Secondary | `NegCtrl_FCT_Perm0` | 86 | -0.0416 | [-0.0541, -0.0305] | +0.1159 | [+0.0758, +0.1525] | +0.2533 |
| dld1 | Secondary | `FCT_FullNTMean` | 86 | +0.1355 | [+0.1059, +0.1645] | +0.1401 | [+0.0817, +0.1933] | +0.2533 |
| dld1 | Secondary | `Transfer_Mean_GlobalScale` | 86 | +nan | [+nan, +nan] | +nan | [+nan, +nan] | +0.2533 |
| dld1 | Secondary | `Transfer_Mean_VarScaled` | 86 | +0.1389 | [+0.1184, +0.1590] | +0.0888 | [+0.0207, +0.1583] | +0.2533 |
| dld1 | Secondary | `CIPHER_Linear` | 86 | +0.0235 | [-0.0357, +0.0822] | -0.0314 | [-0.1101, +0.0576] | +0.2533 |
| dld1 | Secondary | `Hybrid_ProjShared` | 86 | +0.0138 | [-0.0217, +0.0480] | -0.0429 | [-0.1574, +0.0728] | +0.2533 |
| dld1 | Secondary | `Hybrid_on_TMVS` | 86 | +0.0955 | [+0.0645, +0.1263] | +0.0373 | [-0.0727, +0.1431] | +0.2533 |
| dld1 | Secondary | `E1_ConstMed` | 86 | -0.0614 | [-0.0780, -0.0459] | +0.1567 | [+0.1051, +0.2059] | +0.2533 |
| cd4t_rest | Secondary | `Source_Grand_Mean` | 302 | +0.0277 | [-0.0509, +0.1025] | +0.0018 | [-0.0914, +0.0922] | +0.1468 |
| cd4t_rest | Secondary | `Transfer_Mean` | 302 | +nan | [+nan, +nan] | +nan | [+nan, +nan] | +0.1468 |
| cd4t_rest | Secondary | `FCT` | 302 | +0.0556 | [+0.0358, +0.0763] | +0.0560 | [+0.0134, +0.0988] | +0.1468 |
| cd4t_rest | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 302 | +0.0070 | [-0.0074, +0.0225] | +0.0238 | [-0.0161, +0.0629] | +0.1468 |
| cd4t_rest | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 302 | +0.0163 | [-0.0053, +0.0391] | +0.0398 | [-0.0014, +0.0786] | +0.1468 |
| cd4t_rest | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 302 | +0.0017 | [-0.0107, +0.0147] | +0.0213 | [-0.0115, +0.0564] | +0.1468 |
| cd4t_rest | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 302 | +0.0042 | [-0.0123, +0.0217] | +0.0348 | [-0.0037, +0.0728] | +0.1468 |
| cd4t_rest | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 302 | +0.0158 | [-0.0006, +0.0350] | +0.0353 | [-0.0131, +0.0761] | +0.1468 |
| cd4t_rest | Secondary | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 302 | +0.0208 | [-0.0003, +0.0414] | +0.0294 | [-0.0147, +0.0729] | +0.1468 |
| cd4t_rest | Secondary | `NegCtrl_FCT_Perm0` | 302 | +0.0115 | [-0.0000, +0.0222] | +0.0405 | [+0.0131, +0.0683] | +0.1468 |
| cd4t_rest | Secondary | `FCT_FullNTMean` | 302 | +0.0556 | [+0.0358, +0.0763] | +0.0560 | [+0.0134, +0.0988] | +0.1468 |
| cd4t_rest | Secondary | `Transfer_Mean_GlobalScale` | 302 | +nan | [-0.0052, +0.0030] | +nan | [+nan, +nan] | +0.1468 |
| cd4t_rest | Secondary | `Transfer_Mean_VarScaled` | 302 | +0.0658 | [+0.0367, +0.0958] | +0.0198 | [-0.0192, +0.0551] | +0.1468 |
| cd4t_rest | Secondary | `CIPHER_Linear` | 302 | +0.0219 | [-0.0872, +0.1216] | -0.0324 | [-0.1442, +0.0771] | +0.1468 |
| cd4t_rest | Secondary | `Hybrid_ProjShared` | 302 | +0.0137 | [-0.0421, +0.0680] | -0.0244 | [-0.1060, +0.0552] | +0.1468 |
| cd4t_rest | Secondary | `Hybrid_on_TMVS` | 302 | +0.0521 | [+0.0033, +0.1013] | -0.0040 | [-0.0675, +0.0565] | +0.1468 |
| cd4t_rest | Secondary | `E1_ConstMed` | 302 | +0.0130 | [-0.0007, +0.0259] | +0.0495 | [+0.0108, +0.0884] | +0.1468 |

## `H-D3` (`FCT − E1`) & `H-D4` (`E1(resp_top50) − E1(all_nontarget)`) Paired Comparisons

| Target Line | Comparison | Estimate | 95% CI | 98.33% CI |
|---|---|---:|---|---|
| dld1 | `H-D3_FCT_minus_E1_all_nontarget` | +0.1969 | [+0.1709, +0.2219] | [+0.1657, +0.2270] |
| dld1 | `diag_FCT_minus_E1_resp_top50` | -0.0166 | [-0.0660, +0.0313] | [-0.0758, +0.0418] |
| dld1 | `H-D4_E1_resp_top50_minus_all_nontarget` | +0.2180 | [+0.1705, +0.2634] | [+0.1606, +0.2724] |
| cd4t_rest | `H-D3_FCT_minus_E1_all_nontarget` | +0.0426 | [+0.0268, +0.0593] | [+0.0234, +0.0627] |
| cd4t_rest | `diag_FCT_minus_E1_resp_top50` | +0.0065 | [-0.0114, +0.0242] | [-0.0154, +0.0289] |
| cd4t_rest | `H-D4_E1_resp_top50_minus_all_nontarget` | +0.0365 | [+0.0038, +0.0704] | [-0.0038, +0.0778] |

## 100-Perm Null (`all_nontarget` Primary & `resp_top50` Secondary)

- **dld1**: `all_nontarget` median=-0.0428, p97.5=-0.0352; `resp_top50` median=+0.0892, p97.5=+0.1179
- **cd4t_rest**: `all_nontarget` median=+0.0096, p97.5=+0.0131; `resp_top50` median=+0.0347, p97.5=+0.0554
