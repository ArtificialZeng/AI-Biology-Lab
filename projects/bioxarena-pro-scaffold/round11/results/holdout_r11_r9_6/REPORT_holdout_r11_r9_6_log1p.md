# Round 11 Confirmatory Group D (`FCT`, Primary Mask `all_nontarget`) Report (`holdout_r11_r9_6_log1p`)

- **Phase**: `holdout_r11` | **Source Mode**: `r9_6` (`6` sources) | **PREREG_D Commit**: `72fc0c937358b089e78664b464f9db5aafb790b1` | **Code Freeze**: `cc834d863105c66874acb4c27cc214fd79704955`
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
| cd4t_stim48 | Confirmatory | `Source_Grand_Mean` | 503 | +0.2173 | [+0.1684, +0.2667] | +0.1862 | [+0.1411, +0.2314] | +0.0784 |
| cd4t_stim48 | Confirmatory | `Transfer_Mean` | 503 | +nan | [+nan, +nan] | +nan | [+nan, +nan] | +0.0784 |
| cd4t_stim48 | Confirmatory | `FCT` | 503 | +0.0532 | [+0.0337, +0.0731] | +0.0553 | [+0.0292, +0.0808] | +0.0784 |
| cd4t_stim48 | Confirmatory | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 503 | +0.0076 | [-0.0024, +0.0175] | +0.0296 | [+0.0093, +0.0491] | +0.0784 |
| cd4t_stim48 | Confirmatory | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 503 | +0.0158 | [-0.0009, +0.0329] | +0.0169 | [-0.0063, +0.0407] | +0.0784 |
| cd4t_stim48 | Confirmatory | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 503 | +0.0016 | [-0.0089, +0.0121] | +0.0081 | [-0.0140, +0.0297] | +0.0784 |
| cd4t_stim48 | Confirmatory | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 503 | +0.0177 | [+0.0013, +0.0346] | +0.0138 | [-0.0077, +0.0347] | +0.0784 |
| cd4t_stim48 | Confirmatory | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 503 | +0.0015 | [-0.0094, +0.0127] | +0.0150 | [-0.0059, +0.0357] | +0.0784 |
| cd4t_stim48 | Confirmatory | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 503 | +0.0367 | [+0.0161, +0.0572] | +0.0362 | [+0.0082, +0.0631] | +0.0784 |
| cd4t_stim48 | Confirmatory | `NegCtrl_FCT_Perm0` | 503 | -0.0175 | [-0.0272, -0.0080] | +0.0347 | [+0.0133, +0.0553] | +0.0784 |
| cd4t_stim48 | Confirmatory | `FCT_FullNTMean` | 503 | +0.0532 | [+0.0337, +0.0731] | +0.0553 | [+0.0292, +0.0808] | +0.0784 |
| cd4t_stim48 | Confirmatory | `Transfer_Mean_GlobalScale` | 503 | +0.0108 | [+0.0045, +0.0175] | +nan | [+nan, +nan] | +0.0784 |
| cd4t_stim48 | Confirmatory | `Transfer_Mean_VarScaled` | 503 | +0.0851 | [+0.0660, +0.1053] | +0.0280 | [+0.0064, +0.0519] | +0.0784 |
| cd4t_stim48 | Confirmatory | `CIPHER_Linear` | 503 | +0.2543 | [+0.1631, +0.3526] | +0.1965 | [+0.1281, +0.2677] | +0.0784 |
| cd4t_stim48 | Confirmatory | `Hybrid_ProjShared` | 503 | +0.1870 | [+0.1279, +0.2478] | +0.1622 | [+0.1079, +0.2184] | +0.0784 |
| cd4t_stim48 | Confirmatory | `Hybrid_on_TMVS` | 503 | +0.1932 | [+0.1417, +0.2448] | +0.1481 | [+0.1020, +0.1955] | +0.0784 |
| cd4t_stim48 | Confirmatory | `E1_ConstMed` | 503 | -0.0203 | [-0.0326, -0.0087] | +0.0522 | [+0.0257, +0.0792] | +0.0784 |
| cd4t_stim48 | Confirmatory | `E2_ConstP10` | 503 | -0.0178 | [-0.0302, -0.0062] | +0.0507 | [+0.0239, +0.0780] | +0.0784 |
| cd4t_stim48 | Confirmatory | `E2_ConstP90` | 503 | -0.0226 | [-0.0350, -0.0104] | +0.0501 | [+0.0235, +0.0757] | +0.0784 |
| cd4t_stim48 | Confirmatory | `E3_SrcMeanCtrl` | 503 | +0.0196 | [+0.0029, +0.0366] | +0.0296 | [+0.0045, +0.0539] | +0.0784 |
| cd4t_stim48 | Confirmatory | `E5_RankMatchedMu` | 503 | +0.0141 | [-0.0028, +0.0309] | +0.0356 | [+0.0097, +0.0614] | +0.0784 |
| cd4t_stim48 | Confirmatory | `E6_LogFCBar` | 503 | -0.0238 | [-0.0362, -0.0119] | +0.0436 | [+0.0179, +0.0683] | +0.0784 |
| cd4t_stim48 | Confirmatory | `E7_FCBarMinus1` | 503 | -0.0159 | [-0.0285, -0.0038] | +0.0453 | [+0.0191, +0.0724] | +0.0784 |
| feng_ipsc | Secondary | `Source_Grand_Mean` | 93 | -0.1852 | [-0.2333, -0.1374] | -0.2827 | [-0.3512, -0.2129] | +0.5570 |
| feng_ipsc | Secondary | `Transfer_Mean` | 93 | +nan | [+nan, +nan] | +nan | [+nan, +nan] | +0.5570 |
| feng_ipsc | Secondary | `FCT` | 93 | +0.1212 | [+0.1043, +0.1371] | +0.2375 | [+0.1846, +0.2883] | +0.5570 |
| feng_ipsc | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 93 | +0.0437 | [+0.0320, +0.0563] | +0.0928 | [+0.0515, +0.1349] | +0.5570 |
| feng_ipsc | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 93 | -0.0024 | [-0.0131, +0.0097] | +0.0167 | [-0.0277, +0.0561] | +0.5570 |
| feng_ipsc | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 93 | +0.0562 | [+0.0378, +0.0744] | +0.1789 | [+0.1330, +0.2200] | +0.5570 |
| feng_ipsc | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 93 | +0.0069 | [-0.0041, +0.0179] | -0.0019 | [-0.0443, +0.0391] | +0.5570 |
| feng_ipsc | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 93 | -0.0234 | [-0.0367, -0.0112] | -0.0548 | [-0.1035, -0.0177] | +0.5570 |
| feng_ipsc | Secondary | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 93 | +0.0669 | [+0.0507, +0.0829] | +0.1617 | [+0.1200, +0.2030] | +0.5570 |
| feng_ipsc | Secondary | `NegCtrl_FCT_Perm0` | 93 | +0.0325 | [+0.0161, +0.0502] | +0.1216 | [+0.0855, +0.1591] | +0.5570 |
| feng_ipsc | Secondary | `FCT_FullNTMean` | 93 | +0.1212 | [+0.1043, +0.1371] | +0.2375 | [+0.1846, +0.2883] | +0.5570 |
| feng_ipsc | Secondary | `Transfer_Mean_GlobalScale` | 93 | +nan | [+nan, +nan] | +nan | [+nan, +nan] | +0.5570 |
| feng_ipsc | Secondary | `Transfer_Mean_VarScaled` | 93 | +0.1330 | [+0.1146, +0.1540] | +0.2143 | [+0.1590, +0.2738] | +0.5570 |
| feng_ipsc | Secondary | `CIPHER_Linear` | 93 | -0.0204 | [-0.0544, +0.0106] | -0.0255 | [-0.0735, +0.0198] | +0.5570 |
| feng_ipsc | Secondary | `Hybrid_ProjShared` | 93 | +0.0152 | [-0.0033, +0.0321] | +0.0307 | [-0.0155, +0.0743] | +0.5570 |
| feng_ipsc | Secondary | `Hybrid_on_TMVS` | 93 | +0.0805 | [+0.0651, +0.0950] | +0.1732 | [+0.1188, +0.2257] | +0.5570 |
| feng_ipsc | Secondary | `E1_ConstMed` | 93 | +0.0530 | [+0.0287, +0.0785] | +0.2479 | [+0.2003, +0.2975] | +0.5570 |
| feng_ipsc | Secondary | `E2_ConstP10` | 93 | +0.0523 | [+0.0284, +0.0775] | +0.2280 | [+0.1803, +0.2791] | +0.5570 |
| feng_ipsc | Secondary | `E2_ConstP90` | 93 | +0.0552 | [+0.0301, +0.0825] | +0.2654 | [+0.2174, +0.3144] | +0.5570 |
| feng_ipsc | Secondary | `E3_SrcMeanCtrl` | 93 | +0.0579 | [+0.0425, +0.0737] | +0.1721 | [+0.1249, +0.2161] | +0.5570 |
| feng_ipsc | Secondary | `E5_RankMatchedMu` | 93 | +0.0473 | [+0.0328, +0.0614] | +0.1232 | [+0.0761, +0.1650] | +0.5570 |
| feng_ipsc | Secondary | `E6_LogFCBar` | 93 | +0.0582 | [+0.0320, +0.0875] | +0.2519 | [+0.2047, +0.2983] | +0.5570 |
| feng_ipsc | Secondary | `E7_FCBarMinus1` | 93 | +0.0519 | [+0.0285, +0.0766] | +0.1961 | [+0.1469, +0.2517] | +0.5570 |
| telohaec | Secondary | `Source_Grand_Mean` | 9 | -0.0524 | [-0.1485, +0.0338] | -0.0306 | [-0.1913, +0.1344] | +0.1248 |
| telohaec | Secondary | `Transfer_Mean` | 9 | +nan | [+nan, +nan] | +nan | [+nan, +nan] | +0.1248 |
| telohaec | Secondary | `FCT` | 9 | +0.2173 | [+0.0958, +0.3602] | +0.1726 | [-0.0320, +0.3783] | +0.1248 |
| telohaec | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 9 | +0.0837 | [+0.0289, +0.1712] | +0.0310 | [-0.0543, +0.1348] | +0.1248 |
| telohaec | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 9 | -0.0435 | [-0.0943, +0.0140] | -0.1171 | [-0.1785, -0.0462] | +0.1248 |
| telohaec | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 9 | +0.0042 | [-0.0253, +0.0441] | +0.0676 | [-0.0350, +0.2151] | +0.1248 |
| telohaec | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 9 | -0.0103 | [-0.0566, +0.0703] | -0.0667 | [-0.1539, +0.0540] | +0.1248 |
| telohaec | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 9 | +0.0770 | [+0.0116, +0.1633] | -0.0575 | [-0.2212, +0.1611] | +0.1248 |
| telohaec | Secondary | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 9 | -0.0679 | [-0.1800, -0.0144] | +0.0579 | [-0.0841, +0.2152] | +0.1248 |
| telohaec | Secondary | `NegCtrl_FCT_Perm0` | 9 | -0.0456 | [-0.1210, -0.0066] | +0.0389 | [-0.0315, +0.1318] | +0.1248 |
| telohaec | Secondary | `FCT_FullNTMean` | 9 | +0.2173 | [+0.0958, +0.3602] | +0.1726 | [-0.0320, +0.3783] | +0.1248 |
| telohaec | Secondary | `Transfer_Mean_GlobalScale` | 9 | +nan | [+nan, +nan] | +nan | [+nan, +nan] | +0.1248 |
| telohaec | Secondary | `Transfer_Mean_VarScaled` | 9 | +0.1309 | [+0.0710, +0.2212] | +0.1165 | [+0.0576, +0.1790] | +0.1248 |
| telohaec | Secondary | `CIPHER_Linear` | 9 | -0.0062 | [-0.0796, +0.1118] | +0.0306 | [-0.2120, +0.2422] | +0.1248 |
| telohaec | Secondary | `Hybrid_ProjShared` | 9 | +0.0180 | [-0.0556, +0.1103] | +0.0413 | [-0.2196, +0.2868] | +0.1248 |
| telohaec | Secondary | `Hybrid_on_TMVS` | 9 | +0.1049 | [+0.0215, +0.1880] | +0.1199 | [-0.0443, +0.2456] | +0.1248 |
| telohaec | Secondary | `E1_ConstMed` | 9 | -0.0475 | [-0.1700, +0.0187] | +0.1798 | [+0.0334, +0.3756] | +0.1248 |
| telohaec | Secondary | `E2_ConstP10` | 9 | -0.0432 | [-0.1712, +0.0250] | +0.1598 | [+0.0177, +0.3670] | +0.1248 |
| telohaec | Secondary | `E2_ConstP90` | 9 | -0.0531 | [-0.1637, +0.0060] | +0.2047 | [+0.0462, +0.3899] | +0.1248 |
| telohaec | Secondary | `E3_SrcMeanCtrl` | 9 | +0.0073 | [-0.0352, +0.0614] | +0.0402 | [-0.1020, +0.1931] | +0.1248 |
| telohaec | Secondary | `E5_RankMatchedMu` | 9 | +0.0506 | [-0.0117, +0.1606] | -0.0342 | [-0.1454, +0.0680] | +0.1248 |
| telohaec | Secondary | `E6_LogFCBar` | 9 | -0.0490 | [-0.1354, -0.0039] | +0.1970 | [+0.0286, +0.3723] | +0.1248 |
| telohaec | Secondary | `E7_FCBarMinus1` | 9 | -0.0396 | [-0.1685, +0.0295] | +0.1420 | [-0.0009, +0.3560] | +0.1248 |

## `H-D3` (`FCT − E1`) & `H-D4` (`E1(resp_top50) − E1(all_nontarget)`) Paired Comparisons

| Target Line | Comparison | Estimate | 95% CI | 97.5% CI |
|---|---|---:|---|---|
| cd4t_stim48 | `H-D3_FCT_minus_E1_all_nontarget` | +0.0735 | [+0.0585, +0.0892] | [+0.0561, +0.0917] |
| cd4t_stim48 | `diag_FCT_minus_E1_resp_top50` | +0.0030 | [-0.0056, +0.0119] | [-0.0069, +0.0133] |
| cd4t_stim48 | `H-D4_E1_resp_top50_minus_all_nontarget` | +0.0725 | [+0.0506, +0.0943] | [+0.0477, +0.0982] |
| cd4t_stim48 | `H-D5_FCT_minus_Mis_kolf21j_all_nontarget` | +0.0516 | [+0.0380, +0.0655] | [+0.0360, +0.0672] |
| feng_ipsc | `H-D3_FCT_minus_E1_all_nontarget` | +0.0682 | [+0.0468, +0.0896] | [+0.0427, +0.0923] |
| feng_ipsc | `diag_FCT_minus_E1_resp_top50` | -0.0104 | [-0.0402, +0.0182] | [-0.0444, +0.0210] |
| feng_ipsc | `H-D4_E1_resp_top50_minus_all_nontarget` | +0.1949 | [+0.1510, +0.2372] | [+0.1457, +0.2429] |
| feng_ipsc | `H-D5_FCT_minus_Mis_kolf21j_all_nontarget` | +0.0650 | [+0.0543, +0.0766] | [+0.0528, +0.0784] |
| telohaec | `H-D3_FCT_minus_E1_all_nontarget` | +0.2648 | [+0.1588, +0.4286] | [+0.1443, +0.4636] |
| telohaec | `diag_FCT_minus_E1_resp_top50` | -0.0073 | [-0.1523, +0.1475] | [-0.1684, +0.1704] |
| telohaec | `H-D4_E1_resp_top50_minus_all_nontarget` | +0.2273 | [+0.0805, +0.4393] | [+0.0659, +0.4838] |
| telohaec | `H-D5_FCT_minus_Mis_kolf21j_all_nontarget` | +0.2130 | [+0.1111, +0.3271] | [+0.0953, +0.3476] |

## `SA-8` Split-Half Out-of-Sample Top-50 (`resp_top50_split`, Single-Half Ceiling `sqrt(rho_AB)`)

| Target Line | Method | N | `r_orth_raw_split` | `r_orth_cal_split` | 97.5% CI | `Perm100 split` Med / p97.5 |
|---|---|---:|---:|---:|---|---|
| cd4t_stim48 | `FCT` | 503 | +0.0308 | +0.0621 | [+0.0231, +0.1031] | +0.0381 / +0.0571 |
| cd4t_stim48 | `E1_ConstMed` | 503 | +0.0268 | +0.0541 | [+0.0193, +0.0896] | — |
| feng_ipsc | `FCT` | 93 | +0.2309 | +0.2440 | [+0.1909, +0.2936] | +0.1482 / +0.1874 |
| feng_ipsc | `E1_ConstMed` | 93 | +0.2335 | +0.2468 | [+0.1982, +0.2981] | — |
| telohaec | `FCT` | 9 | +0.1719 | +0.2579 | [-0.0858, +0.6027] | +0.1070 / +0.2343 |
| telohaec | `E1_ConstMed` | 9 | +0.1459 | +0.2188 | [-0.0366, +0.4558] | — |

## Descriptive `H-D5` (`feng_ipsc` Non-KOLF vs `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)`, Non-Confirmatory)

- **`feng_ipsc`** (N=93, G_eval=5326): `FCT` = +0.1212 [+0.1043, +0.1371] vs `Mis(kolf21j)` = +0.0562 [+0.0378, +0.0744]; paired `FCT − Mis(kolf21j)` = +0.0650 [+0.0528, +0.0784].
- *Addendum A3 Notice:* Split halves share genetic donors across sibling iPSC clones; noise ceiling reflects clone/channel variability rather than inter-donor variation.

## 100-Perm Null (`all_nontarget` Primary, `resp_top50`, & `SA-8 resp_top50_split`)

- **cd4t_stim48**: `all_nontarget` median=-0.0164, p97.5=-0.0129; `resp_top50` median=+0.0373, p97.5=+0.0481; `resp_top50_split` median=+0.0381, p97.5=+0.0571
- **feng_ipsc**: `all_nontarget` median=+0.0380, p97.5=+0.0452; `resp_top50` median=+0.1477, p97.5=+0.1876; `resp_top50_split` median=+0.1482, p97.5=+0.1874
- **telohaec**: `all_nontarget` median=-0.0306, p97.5=-0.0132; `resp_top50` median=+0.0908, p97.5=+0.1723; `resp_top50_split` median=+0.1070, p97.5=+0.2343
