# Round 11 Confirmatory Group D (`FCT`, Primary Mask `all_nontarget`) Report (`holdout_r11_r9_r10_8_log1p`)

- **Phase**: `holdout_r11` | **Source Mode**: `r9_r10_8` (`8` sources) | **PREREG_D Commit**: `72fc0c937358b089e78664b464f9db5aafb790b1` | **Code Freeze**: `cc834d863105c66874acb4c27cc214fd79704955`
- **Bootstrap**: `np.random.default_rng(20261123)` x `4000` (Bonferroni 2-claim + IUT two-sided `97.5%` CIs `[1.25%, 98.75%]`)
- **Cohort Disclosure**: **Insufficient confirmatory cell lines (K_conf < 2); confirmatory hypothesis testing was not performed (确认性线不足，未做确认性检验).** (`['cd4t_stim48']`)
- **Secondary Lines**: `['feng_ipsc', 'telohaec']`
- **Confirmatory Claim 1 — Primary Group D (`FCT` 全基因线特异信号) Verdict**: **NOT_EVALUATED_INSUFFICIENT_LINES**
- **Confirmatory Claim 2 — Methodological Control (`H-D4` `resp_top50` Selection Artifact)**: **NOT_EVALUATED_INSUFFICIENT_LINES** (all lines descriptive: `CONFIRMED`)
- **Parallel Practical Utility (`H-D_util` / `H-C3`-style Non-Inferiority)**: **NOT_EVALUATED_INSUFFICIENT_LINES**

- **Mandatory `cd4t_stim48` Disclosures (`Q7` & `Addendum A2`)**: `cd4t_stim48` shares primary human CD4+ alpha-beta T lineage, sgRNA library, and Marson lab with Round 10 `cd4t_rest` (weaker biological/laboratory independence). Furthermore, both donors (`CE0008678`, `CE0010866`) share 24 processing channels (`CD4i_R2L25..R2L48`), so residual channel noise slightly elevates split-half noise ceilings (`r_orth_cal` biased conservatively downward).

## Per-Line Primary Mask (`all_nontarget`, Joint-Bootstrap 97.5% CIs) & Secondary (`resp_top50`)

| Target Line | Role | Method | N | Primary `r_orth_cal` | Primary 97.5% CI | `resp_top50` `r_orth_cal` | `resp_top50` 97.5% CI | `rho_AB_orth` |
|---|---|---|---:|---:|---|---:|---|---:|
| cd4t_stim48 | Confirmatory | `Source_Grand_Mean` | 1406 | +0.2204 | [+0.1831, +0.2554] | +0.1779 | [+0.1473, +0.2069] | +0.0614 |
| cd4t_stim48 | Confirmatory | `Transfer_Mean` | 1406 | +nan | [+nan, +nan] | +nan | [+nan, +nan] | +0.0614 |
| cd4t_stim48 | Confirmatory | `FCT` | 1406 | +0.0547 | [+0.0403, +0.0700] | +0.0388 | [+0.0219, +0.0559] | +0.0614 |
| cd4t_stim48 | Confirmatory | `NegCtrl_FCT_MismatchedMean(cd4t_rest)` | 1406 | +0.0214 | [+0.0098, +0.0333] | +0.0260 | [+0.0086, +0.0432] | +0.0614 |
| cd4t_stim48 | Confirmatory | `NegCtrl_FCT_MismatchedMean(dld1)` | 1406 | +0.0159 | [+0.0080, +0.0237] | +0.0134 | [-0.0010, +0.0270] | +0.0614 |
| cd4t_stim48 | Confirmatory | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 1406 | +0.0181 | [+0.0096, +0.0260] | +0.0195 | [+0.0056, +0.0336] | +0.0614 |
| cd4t_stim48 | Confirmatory | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 1406 | +0.0222 | [+0.0115, +0.0329] | +0.0116 | [-0.0030, +0.0263] | +0.0614 |
| cd4t_stim48 | Confirmatory | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 1406 | +0.0100 | [+0.0028, +0.0171] | +0.0113 | [-0.0040, +0.0264] | +0.0614 |
| cd4t_stim48 | Confirmatory | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 1406 | +0.0240 | [+0.0127, +0.0354] | +0.0148 | [-0.0001, +0.0292] | +0.0614 |
| cd4t_stim48 | Confirmatory | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 1406 | +0.0149 | [+0.0076, +0.0221] | +0.0108 | [-0.0018, +0.0244] | +0.0614 |
| cd4t_stim48 | Confirmatory | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 1406 | +0.0317 | [+0.0193, +0.0443] | +0.0262 | [+0.0089, +0.0426] | +0.0614 |
| cd4t_stim48 | Confirmatory | `NegCtrl_FCT_Perm0` | 1406 | -0.0166 | [-0.0217, -0.0117] | +0.0369 | [+0.0235, +0.0491] | +0.0614 |
| cd4t_stim48 | Confirmatory | `FCT_FullNTMean` | 1406 | +0.0547 | [+0.0403, +0.0700] | +0.0388 | [+0.0219, +0.0559] | +0.0614 |
| cd4t_stim48 | Confirmatory | `Transfer_Mean_GlobalScale` | 1406 | +0.0058 | [+0.0016, +0.0102] | +nan | [+nan, +nan] | +0.0614 |
| cd4t_stim48 | Confirmatory | `Transfer_Mean_VarScaled` | 1406 | +0.0564 | [+0.0457, +0.0673] | +0.0054 | [-0.0093, +0.0206] | +0.0614 |
| cd4t_stim48 | Confirmatory | `CIPHER_Linear` | 1406 | +0.2385 | [+0.1687, +0.3047] | +0.1704 | [+0.1232, +0.2163] | +0.0614 |
| cd4t_stim48 | Confirmatory | `Hybrid_ProjShared` | 1406 | +0.1622 | [+0.1181, +0.2061] | +0.1276 | [+0.0901, +0.1649] | +0.0614 |
| cd4t_stim48 | Confirmatory | `Hybrid_on_TMVS` | 1406 | +0.1510 | [+0.1187, +0.1842] | +0.0927 | [+0.0640, +0.1214] | +0.0614 |
| cd4t_stim48 | Confirmatory | `E1_ConstMed` | 1406 | -0.0212 | [-0.0275, -0.0152] | +0.0414 | [+0.0252, +0.0567] | +0.0614 |
| cd4t_stim48 | Confirmatory | `E2_ConstP10` | 1406 | -0.0174 | [-0.0238, -0.0113] | +0.0369 | [+0.0199, +0.0528] | +0.0614 |
| cd4t_stim48 | Confirmatory | `E2_ConstP90` | 1406 | -0.0237 | [-0.0301, -0.0174] | +0.0415 | [+0.0255, +0.0569] | +0.0614 |
| cd4t_stim48 | Confirmatory | `E3_SrcMeanCtrl` | 1406 | +0.0214 | [+0.0128, +0.0300] | +0.0203 | [+0.0058, +0.0346] | +0.0614 |
| cd4t_stim48 | Confirmatory | `E5_RankMatchedMu` | 1406 | +0.0209 | [+0.0122, +0.0298] | +0.0255 | [+0.0104, +0.0401] | +0.0614 |
| cd4t_stim48 | Confirmatory | `E6_LogFCBar` | 1406 | -0.0245 | [-0.0305, -0.0184] | +0.0379 | [+0.0220, +0.0532] | +0.0614 |
| cd4t_stim48 | Confirmatory | `E7_FCBarMinus1` | 1406 | -0.0125 | [-0.0187, -0.0065] | +0.0263 | [+0.0079, +0.0464] | +0.0614 |
| feng_ipsc | Secondary | `Source_Grand_Mean` | 95 | -0.1698 | [-0.2171, -0.1223] | -0.2597 | [-0.3309, -0.1885] | +0.5536 |
| feng_ipsc | Secondary | `Transfer_Mean` | 95 | +nan | [+nan, +nan] | +nan | [+nan, +nan] | +0.5536 |
| feng_ipsc | Secondary | `FCT` | 95 | +0.1210 | [+0.1043, +0.1378] | +0.2327 | [+0.1771, +0.2879] | +0.5536 |
| feng_ipsc | Secondary | `NegCtrl_FCT_MismatchedMean(cd4t_rest)` | 95 | +0.0208 | [+0.0068, +0.0341] | +0.0847 | [+0.0514, +0.1179] | +0.5536 |
| feng_ipsc | Secondary | `NegCtrl_FCT_MismatchedMean(dld1)` | 95 | -0.0033 | [-0.0114, +0.0052] | +0.0058 | [-0.0213, +0.0338] | +0.5536 |
| feng_ipsc | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 95 | +0.0477 | [+0.0352, +0.0613] | +0.0916 | [+0.0487, +0.1355] | +0.5536 |
| feng_ipsc | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 95 | +0.0027 | [-0.0091, +0.0152] | +0.0221 | [-0.0215, +0.0617] | +0.5536 |
| feng_ipsc | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 95 | +0.0568 | [+0.0384, +0.0748] | +0.1767 | [+0.1302, +0.2192] | +0.5536 |
| feng_ipsc | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 95 | +0.0113 | [+0.0000, +0.0232] | +0.0010 | [-0.0419, +0.0435] | +0.5536 |
| feng_ipsc | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 95 | -0.0170 | [-0.0306, -0.0032] | -0.0489 | [-0.0948, -0.0136] | +0.5536 |
| feng_ipsc | Secondary | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 95 | +0.0663 | [+0.0504, +0.0823] | +0.1622 | [+0.1223, +0.2030] | +0.5536 |
| feng_ipsc | Secondary | `NegCtrl_FCT_Perm0` | 95 | +0.0403 | [+0.0228, +0.0587] | +0.1706 | [+0.1324, +0.2103] | +0.5536 |
| feng_ipsc | Secondary | `FCT_FullNTMean` | 95 | +0.1210 | [+0.1043, +0.1378] | +0.2327 | [+0.1771, +0.2879] | +0.5536 |
| feng_ipsc | Secondary | `Transfer_Mean_GlobalScale` | 95 | +nan | [+nan, +nan] | +nan | [+nan, +nan] | +0.5536 |
| feng_ipsc | Secondary | `Transfer_Mean_VarScaled` | 95 | +0.1290 | [+0.1100, +0.1492] | +0.2054 | [+0.1484, +0.2650] | +0.5536 |
| feng_ipsc | Secondary | `CIPHER_Linear` | 95 | -0.0161 | [-0.0489, +0.0155] | -0.0168 | [-0.0606, +0.0270] | +0.5536 |
| feng_ipsc | Secondary | `Hybrid_ProjShared` | 95 | +0.0141 | [-0.0043, +0.0322] | +0.0259 | [-0.0212, +0.0717] | +0.5536 |
| feng_ipsc | Secondary | `Hybrid_on_TMVS` | 95 | +0.1014 | [+0.0865, +0.1165] | +0.1926 | [+0.1370, +0.2499] | +0.5536 |
| feng_ipsc | Secondary | `E1_ConstMed` | 95 | +0.0523 | [+0.0281, +0.0767] | +0.2415 | [+0.1891, +0.2945] | +0.5536 |
| feng_ipsc | Secondary | `E2_ConstP10` | 95 | +0.0517 | [+0.0280, +0.0753] | +0.2227 | [+0.1714, +0.2767] | +0.5536 |
| feng_ipsc | Secondary | `E2_ConstP90` | 95 | +0.0544 | [+0.0292, +0.0806] | +0.2582 | [+0.2047, +0.3098] | +0.5536 |
| feng_ipsc | Secondary | `E3_SrcMeanCtrl` | 95 | +0.0485 | [+0.0347, +0.0627] | +0.1656 | [+0.1213, +0.2088] | +0.5536 |
| feng_ipsc | Secondary | `E5_RankMatchedMu` | 95 | +0.0351 | [+0.0219, +0.0477] | +0.1006 | [+0.0603, +0.1404] | +0.5536 |
| feng_ipsc | Secondary | `E6_LogFCBar` | 95 | +0.0572 | [+0.0304, +0.0849] | +0.2459 | [+0.1952, +0.2951] | +0.5536 |
| feng_ipsc | Secondary | `E7_FCBarMinus1` | 95 | +0.0514 | [+0.0286, +0.0743] | +0.1925 | [+0.1407, +0.2484] | +0.5536 |
| telohaec | Secondary | `Source_Grand_Mean` | 9 | -0.0680 | [-0.1626, +0.0261] | -0.0232 | [-0.1890, +0.1474] | +0.1285 |
| telohaec | Secondary | `Transfer_Mean` | 9 | +nan | [+nan, +nan] | +nan | [+nan, +nan] | +0.1285 |
| telohaec | Secondary | `FCT` | 9 | +0.2256 | [+0.1198, +0.3587] | +0.1919 | [+0.0408, +0.3699] | +0.1285 |
| telohaec | Secondary | `NegCtrl_FCT_MismatchedMean(cd4t_rest)` | 9 | -0.0344 | [-0.1276, +0.0156] | -0.0098 | [-0.1339, +0.1230] | +0.1285 |
| telohaec | Secondary | `NegCtrl_FCT_MismatchedMean(dld1)` | 9 | +0.0430 | [+0.0169, +0.0864] | +0.0559 | [-0.0658, +0.1924] | +0.1285 |
| telohaec | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 9 | +0.0985 | [+0.0389, +0.1881] | +0.0665 | [-0.0089, +0.1509] | +0.1285 |
| telohaec | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 9 | -0.0210 | [-0.0780, +0.0539] | -0.0745 | [-0.1233, +0.0039] | +0.1285 |
| telohaec | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 9 | +0.0128 | [-0.0186, +0.0534] | +0.0677 | [-0.0391, +0.2086] | +0.1285 |
| telohaec | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 9 | +0.0099 | [-0.0506, +0.1104] | -0.0330 | [-0.1343, +0.0941] | +0.1285 |
| telohaec | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 9 | +0.0976 | [+0.0239, +0.1920] | -0.0121 | [-0.1858, +0.1837] | +0.1285 |
| telohaec | Secondary | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 9 | -0.0661 | [-0.1774, -0.0127] | +0.0472 | [-0.1181, +0.2114] | +0.1285 |
| telohaec | Secondary | `NegCtrl_FCT_Perm0` | 9 | -0.0442 | [-0.1109, -0.0122] | +0.0829 | [+0.0002, +0.1977] | +0.1285 |
| telohaec | Secondary | `FCT_FullNTMean` | 9 | +0.2256 | [+0.1198, +0.3587] | +0.1919 | [+0.0408, +0.3699] | +0.1285 |
| telohaec | Secondary | `Transfer_Mean_GlobalScale` | 9 | +nan | [+nan, +nan] | +nan | [+nan, +nan] | +0.1285 |
| telohaec | Secondary | `Transfer_Mean_VarScaled` | 9 | +0.1203 | [+0.0525, +0.2119] | +0.1784 | [+0.0754, +0.2762] | +0.1285 |
| telohaec | Secondary | `CIPHER_Linear` | 9 | -0.0005 | [-0.0785, +0.1222] | +0.0363 | [-0.2103, +0.2525] | +0.1285 |
| telohaec | Secondary | `Hybrid_ProjShared` | 9 | +0.0339 | [-0.0535, +0.1346] | +0.0601 | [-0.2192, +0.2971] | +0.1285 |
| telohaec | Secondary | `Hybrid_on_TMVS` | 9 | +0.1169 | [+0.0359, +0.1909] | +0.1801 | [-0.0138, +0.3118] | +0.1285 |
| telohaec | Secondary | `E1_ConstMed` | 9 | -0.0423 | [-0.1615, +0.0211] | +0.1808 | [+0.0428, +0.3637] | +0.1285 |
| telohaec | Secondary | `E2_ConstP10` | 9 | -0.0379 | [-0.1622, +0.0289] | +0.1591 | [+0.0175, +0.3535] | +0.1285 |
| telohaec | Secondary | `E2_ConstP90` | 9 | -0.0483 | [-0.1523, +0.0088] | +0.2078 | [+0.0634, +0.3758] | +0.1285 |
| telohaec | Secondary | `E3_SrcMeanCtrl` | 9 | +0.0262 | [-0.0178, +0.0785] | +0.1009 | [-0.0153, +0.2270] | +0.1285 |
| telohaec | Secondary | `E5_RankMatchedMu` | 9 | +0.0777 | [+0.0085, +0.1940] | +0.0297 | [-0.0699, +0.1483] | +0.1285 |
| telohaec | Secondary | `E6_LogFCBar` | 9 | -0.0457 | [-0.1216, -0.0029] | +0.2032 | [+0.0528, +0.3604] | +0.1285 |
| telohaec | Secondary | `E7_FCBarMinus1` | 9 | -0.0340 | [-0.1604, +0.0342] | +0.1396 | [-0.0070, +0.3443] | +0.1285 |

## `H-D3` (`FCT − E1`) & `H-D4` (`E1(resp_top50) − E1(all_nontarget)`) Paired Comparisons

| Target Line | Comparison | Estimate | 95% CI | 97.5% CI |
|---|---|---:|---|---|
| cd4t_stim48 | `H-D3_FCT_minus_E1_all_nontarget` | +0.0759 | [+0.0633, +0.0887] | [+0.0617, +0.0911] |
| cd4t_stim48 | `diag_FCT_minus_E1_resp_top50` | -0.0027 | [-0.0094, +0.0041] | [-0.0102, +0.0052] |
| cd4t_stim48 | `H-D4_E1_resp_top50_minus_all_nontarget` | +0.0626 | [+0.0483, +0.0762] | [+0.0465, +0.0776] |
| cd4t_stim48 | `H-D5_FCT_minus_Mis_kolf21j_all_nontarget` | +0.0447 | [+0.0343, +0.0554] | [+0.0332, +0.0570] |
| feng_ipsc | `H-D3_FCT_minus_E1_all_nontarget` | +0.0687 | [+0.0480, +0.0901] | [+0.0450, +0.0942] |
| feng_ipsc | `diag_FCT_minus_E1_resp_top50` | -0.0089 | [-0.0362, +0.0185] | [-0.0406, +0.0228] |
| feng_ipsc | `H-D4_E1_resp_top50_minus_all_nontarget` | +0.1892 | [+0.1449, +0.2355] | [+0.1392, +0.2412] |
| feng_ipsc | `H-D5_FCT_minus_Mis_kolf21j_all_nontarget` | +0.0642 | [+0.0535, +0.0757] | [+0.0518, +0.0774] |
| telohaec | `H-D3_FCT_minus_E1_all_nontarget` | +0.2680 | [+0.1670, +0.4240] | [+0.1551, +0.4585] |
| telohaec | `diag_FCT_minus_E1_resp_top50` | +0.0112 | [-0.1089, +0.1463] | [-0.1195, +0.1627] |
| telohaec | `H-D4_E1_resp_top50_minus_all_nontarget` | +0.2231 | [+0.0805, +0.4275] | [+0.0651, +0.4590] |
| telohaec | `H-D5_FCT_minus_Mis_kolf21j_all_nontarget` | +0.2128 | [+0.1220, +0.3154] | [+0.1067, +0.3311] |

## `SA-8` Split-Half Out-of-Sample Top-50 (`resp_top50_split`, Single-Half Ceiling `sqrt(rho_AB)`)

| Target Line | Method | N | `r_orth_raw_split` | `r_orth_cal_split` | 97.5% CI | `Perm100 split` Med / p97.5 |
|---|---|---:|---:|---:|---|---|
| cd4t_stim48 | `FCT` | 1406 | +0.0232 | +0.0549 | [+0.0280, +0.0828] | +0.0352 / +0.0481 |
| cd4t_stim48 | `E1_ConstMed` | 1406 | +0.0212 | +0.0501 | [+0.0286, +0.0720] | — |
| feng_ipsc | `FCT` | 95 | +0.2278 | +0.2406 | [+0.1857, +0.2943] | +0.1442 / +0.1806 |
| feng_ipsc | `E1_ConstMed` | 95 | +0.2297 | +0.2427 | [+0.1921, +0.2935] | — |
| telohaec | `FCT` | 9 | +0.1931 | +0.2866 | [+0.0271, +0.5790] | +0.1168 / +0.2060 |
| telohaec | `E1_ConstMed` | 9 | +0.1578 | +0.2351 | [+0.0062, +0.4512] | — |

## Descriptive `H-D5` (`feng_ipsc` Non-KOLF vs `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)`, Non-Confirmatory)

- **`feng_ipsc`** (N=95, G_eval=5297): `FCT` = +0.1210 [+0.1043, +0.1378] vs `Mis(kolf21j)` = +0.0568 [+0.0384, +0.0748]; paired `FCT − Mis(kolf21j)` = +0.0642 [+0.0518, +0.0774].
- *Addendum A3 Notice:* Split halves share genetic donors across sibling iPSC clones; noise ceiling reflects clone/channel variability rather than inter-donor variation.

## 100-Perm Null (`all_nontarget` Primary, `resp_top50`, & `SA-8 resp_top50_split`)

- **cd4t_stim48**: `all_nontarget` median=-0.0176, p97.5=-0.0148; `resp_top50` median=+0.0272, p97.5=+0.0361; `resp_top50_split` median=+0.0352, p97.5=+0.0481
- **feng_ipsc**: `all_nontarget` median=+0.0361, p97.5=+0.0442; `resp_top50` median=+0.1447, p97.5=+0.1779; `resp_top50_split` median=+0.1442, p97.5=+0.1806
- **telohaec**: `all_nontarget` median=-0.0286, p97.5=-0.0112; `resp_top50` median=+0.1001, p97.5=+0.1689; `resp_top50_split` median=+0.1168, p97.5=+0.2060
