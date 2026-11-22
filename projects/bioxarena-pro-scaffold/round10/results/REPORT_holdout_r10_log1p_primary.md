# Round 10 Confirmatory Group C (`FCT`) Report (`holdout_r10_log1p_primary`)

- **Phase**: `holdout_r10` | **PREREG_C Commit**: `PREREG_C_COMMIT（主控核验通过）` | **Code Freeze**: `f8d2ebfc238d83a0c3a6c86696950ed292bf90e9`
- **Bootstrap**: `np.random.default_rng(20261122)` x `4000` (Bonferroni two-sided `97.5%` CIs)
- **Cohort Disclosure**: **K_conf = 2 confirmatory cell lines** (`['dld1', 'cd4t_rest']`)
- **Secondary Lines**: `['calu3', 'wtc11_i3neuron']`
- **Primary Group C (`FCT`) Verdict**: **未确认**
- **Parallel Practical Utility (`H-C3` Non-Inferiority)**: **CONFIRMED**

## Confirmatory Hypothesis Adjudication (97.5% Bonferroni CIs)

| Hypothesis | Pass / Total | Criterion | Verdict |
|---|---:|---|---|
| **H-C1** (`FCT` `r_orth_cal` 97.5% CI lo > 0) | 2/2 | Unanimous K_conf | **CONFIRMED** |
| **H-C2** (`FCT` `r_orth_cal` lo > max 6 mismatched hi & perm100 p97.5) | 0/2 | Unanimous K_conf | **未确认** |
| **H-C3** (Non-inferiority vs `Transfer_Mean` across 6 conventions) | 12/12 non-inf (`1` strict sup) | 0 inferior pairs | **CONFIRMED** |

## Per-Line Results (`resp_top50`, Joint-Bootstrap 97.5% CIs)

| Target Line | Role | Method | N | r_orth_cal | 97.5% CI | Slope OLS | 97.5% CI | r_cal (raw) | rho_AB_orth | pds_cos |
|---|---|---|---:|---:|---|---:|---|---:|---:|---:|
| calu3 | Secondary | `Source_Grand_Mean` | 8 | -0.1231 | [-0.4372, +0.2562] | +0.0045 | [-0.0594, +0.0823] | +0.0104 | +0.8776 | 0.5000 |
| calu3 | Secondary | `Transfer_Mean` | 8 | +nan | [+nan, +nan] | +0.0000 | [+0.0000, +0.0000] | +nan | +0.8776 | 0.6429 |
| calu3 | Secondary | `FCT` | 8 | +0.0879 | [-0.0689, +0.2743] | +0.0068 | [-0.0181, +0.0400] | +0.0381 | +0.8776 | 0.6964 |
| calu3 | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 8 | +0.0794 | [+0.0212, +0.1486] | +0.0077 | [+0.0011, +0.0160] | +0.0846 | +0.8776 | 0.6429 |
| calu3 | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 8 | -0.0657 | [-0.1408, +0.0030] | -0.0049 | [-0.0188, +0.0088] | -0.0439 | +0.8776 | 0.6607 |
| calu3 | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 8 | +0.0099 | [-0.1332, +0.1394] | +0.0072 | [-0.0081, +0.0270] | +0.0613 | +0.8776 | 0.6429 |
| calu3 | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 8 | +0.0364 | [-0.1106, +0.1956] | +0.0078 | [-0.0073, +0.0249] | +0.0776 | +0.8776 | 0.6250 |
| calu3 | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 8 | +0.0365 | [-0.1087, +0.1573] | +0.0023 | [-0.0052, +0.0126] | +0.0343 | +0.8776 | 0.6250 |
| calu3 | Secondary | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 8 | +0.1134 | [-0.0999, +0.2702] | +0.0217 | [-0.0057, +0.0494] | +0.1189 | +0.8776 | 0.7143 |
| calu3 | Secondary | `NegCtrl_FCT_Perm0` | 8 | +0.1372 | [-0.0315, +0.2845] | +0.0423 | [+0.0079, +0.0852] | +0.1287 | +0.8776 | 0.7321 |
| calu3 | Secondary | `FCT_FullNTMean` | 8 | +0.0879 | [-0.0689, +0.2743] | +0.0068 | [-0.0181, +0.0400] | +0.0381 | +0.8776 | 0.6964 |
| calu3 | Secondary | `Transfer_Mean_GlobalScale` | 8 | +nan | [+nan, +nan] | -0.0063 | [-0.0180, +0.0028] | -0.0822 | +0.8776 | 0.6429 |
| calu3 | Secondary | `Transfer_Mean_VarScaled` | 8 | -0.0352 | [-0.1762, +0.1212] | -0.0513 | [-0.1377, +0.0294] | -0.0999 | +0.8776 | 0.6250 |
| calu3 | Secondary | `CIPHER_Linear` | 8 | +0.2973 | [-0.3593, +0.7573] | +0.2265 | [-0.1487, +0.5966] | +0.3055 | +0.8776 | 0.5893 |
| calu3 | Secondary | `Hybrid_ProjShared` | 8 | +0.2584 | [-0.3067, +0.6849] | +0.1287 | [-0.1513, +0.3706] | +0.2508 | +0.8776 | 0.6607 |
| calu3 | Secondary | `Hybrid_on_TMVS` | 8 | +0.2027 | [-0.3292, +0.6104] | +0.0775 | [-0.2447, +0.3633] | +0.1028 | +0.8776 | 0.6786 |
| cd4t_rest | Confirmatory | `Source_Grand_Mean` | 302 | +0.0018 | [-0.0837, +0.0908] | +0.0731 | [+0.0382, +0.1126] | +0.1928 | +0.6479 | 0.5000 |
| cd4t_rest | Confirmatory | `Transfer_Mean` | 302 | +nan | [+nan, +nan] | +0.0000 | [+0.0000, +0.0000] | +nan | +0.6479 | 0.6378 |
| cd4t_rest | Confirmatory | `FCT` | 302 | +0.0560 | [+0.0171, +0.0938] | +0.0156 | [+0.0026, +0.0313] | +0.0802 | +0.6479 | 0.6378 |
| cd4t_rest | Confirmatory | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 302 | +0.0238 | [-0.0132, +0.0599] | +0.0091 | [+0.0038, +0.0153] | +0.0834 | +0.6479 | 0.6372 |
| cd4t_rest | Confirmatory | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 302 | +0.0398 | [+0.0002, +0.0761] | -0.0006 | [-0.0063, +0.0059] | -0.0050 | +0.6479 | 0.6284 |
| cd4t_rest | Confirmatory | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 302 | +0.0213 | [-0.0095, +0.0548] | +0.0154 | [+0.0093, +0.0219] | +0.1148 | +0.6479 | 0.6355 |
| cd4t_rest | Confirmatory | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 302 | +0.0348 | [-0.0015, +0.0691] | +0.0082 | [+0.0036, +0.0138] | +0.0750 | +0.6479 | 0.6290 |
| cd4t_rest | Confirmatory | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 302 | +0.0353 | [-0.0106, +0.0714] | +0.0030 | [-0.0010, +0.0082] | +0.0352 | +0.6479 | 0.6353 |
| cd4t_rest | Confirmatory | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 302 | +0.0294 | [-0.0094, +0.0680] | +0.0247 | [+0.0044, +0.0517] | +0.1130 | +0.6479 | 0.6314 |
| cd4t_rest | Confirmatory | `NegCtrl_FCT_Perm0` | 302 | +0.0405 | [+0.0156, +0.0656] | +0.0641 | [+0.0397, +0.0942] | +0.2407 | +0.6479 | 0.6418 |
| cd4t_rest | Confirmatory | `FCT_FullNTMean` | 302 | +0.0560 | [+0.0171, +0.0938] | +0.0156 | [+0.0026, +0.0313] | +0.0802 | +0.6479 | 0.6378 |
| cd4t_rest | Confirmatory | `Transfer_Mean_GlobalScale` | 302 | +nan | [+nan, +nan] | -0.0452 | [-0.0637, -0.0281] | -0.2616 | +0.6479 | 0.6378 |
| cd4t_rest | Confirmatory | `Transfer_Mean_VarScaled` | 302 | +0.0198 | [-0.0162, +0.0549] | -0.0896 | [-0.1350, -0.0497] | -0.2445 | +0.6479 | 0.6392 |
| cd4t_rest | Confirmatory | `CIPHER_Linear` | 302 | -0.0324 | [-0.1365, +0.0678] | +0.0954 | [+0.0493, +0.1446] | +0.2182 | +0.6479 | 0.4869 |
| cd4t_rest | Confirmatory | `Hybrid_ProjShared` | 302 | -0.0244 | [-0.1019, +0.0504] | -0.0047 | [-0.0160, +0.0061] | -0.0286 | +0.6479 | 0.6197 |
| cd4t_rest | Confirmatory | `Hybrid_on_TMVS` | 302 | -0.0040 | [-0.0640, +0.0531] | -0.0943 | [-0.1402, -0.0542] | -0.2267 | +0.6479 | 0.6263 |
| dld1 | Confirmatory | `Source_Grand_Mean` | 86 | +0.2222 | [+0.1120, +0.3297] | -0.0576 | [-0.1105, -0.0078] | -0.1093 | +0.8277 | 0.5000 |
| dld1 | Confirmatory | `Transfer_Mean` | 86 | +nan | [+nan, +nan] | +0.0000 | [+0.0000, +0.0000] | +nan | +0.8277 | 0.9023 |
| dld1 | Confirmatory | `FCT` | 86 | +0.1401 | [+0.0870, +0.1917] | +0.0290 | [+0.0194, +0.0396] | +0.2067 | +0.8277 | 0.8969 |
| dld1 | Confirmatory | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 86 | +0.0294 | [-0.0174, +0.0726] | +0.0031 | [-0.0029, +0.0091] | +0.0240 | +0.8277 | 0.9014 |
| dld1 | Confirmatory | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 86 | -0.0225 | [-0.0739, +0.0297] | -0.0008 | [-0.0068, +0.0058] | -0.0067 | +0.8277 | 0.8907 |
| dld1 | Confirmatory | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 86 | +0.0300 | [-0.0138, +0.0759] | -0.0050 | [-0.0124, +0.0031] | -0.0321 | +0.8277 | 0.8967 |
| dld1 | Confirmatory | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 86 | -0.0058 | [-0.0557, +0.0447] | -0.0048 | [-0.0111, +0.0014] | -0.0421 | +0.8277 | 0.8911 |
| dld1 | Confirmatory | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 86 | +0.0463 | [-0.0196, +0.1032] | +0.0089 | [+0.0014, +0.0171] | +0.0816 | +0.8277 | 0.8999 |
| dld1 | Confirmatory | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 86 | +0.0707 | [+0.0331, +0.1092] | -0.0046 | [-0.0206, +0.0135] | -0.0161 | +0.8277 | 0.8967 |
| dld1 | Confirmatory | `NegCtrl_FCT_Perm0` | 86 | +0.1159 | [+0.0765, +0.1509] | -0.0542 | [-0.0823, -0.0279] | -0.1346 | +0.8277 | 0.9059 |
| dld1 | Confirmatory | `FCT_FullNTMean` | 86 | +0.1401 | [+0.0870, +0.1917] | +0.0290 | [+0.0194, +0.0396] | +0.2067 | +0.8277 | 0.8969 |
| dld1 | Confirmatory | `Transfer_Mean_GlobalScale` | 86 | +nan | [+nan, +nan] | +0.0369 | [+0.0248, +0.0503] | +0.2404 | +0.8277 | 0.9023 |
| dld1 | Confirmatory | `Transfer_Mean_VarScaled` | 86 | +0.0888 | [+0.0233, +0.1529] | +0.0536 | [+0.0311, +0.0773] | +0.2149 | +0.8277 | 0.9047 |
| dld1 | Confirmatory | `CIPHER_Linear` | 86 | -0.0314 | [-0.1049, +0.0539] | -0.1590 | [-0.2145, -0.1085] | -0.2441 | +0.8277 | 0.4707 |
| dld1 | Confirmatory | `Hybrid_ProjShared` | 86 | -0.0429 | [-0.1479, +0.0675] | -0.0053 | [-0.0203, +0.0053] | -0.0461 | +0.8277 | 0.8886 |
| dld1 | Confirmatory | `Hybrid_on_TMVS` | 86 | +0.0373 | [-0.0657, +0.1361] | +0.0483 | [+0.0198, +0.0756] | +0.1726 | +0.8277 | 0.8943 |
| wtc11_i3neuron | Secondary | `Source_Grand_Mean` | 2 | -0.1816 | [-0.3404, +0.1441] | +0.2667 | [+0.2667, +0.4511] | +0.4154 | +0.8562 | 0.5000 |
| wtc11_i3neuron | Secondary | `Transfer_Mean` | 2 | +nan | [+nan, +nan] | +0.0000 | [+0.0000, +0.0000] | +nan | +0.8562 | 1.0000 |
| wtc11_i3neuron | Secondary | `FCT` | 2 | -0.0313 | [-0.0334, +0.2649] | -0.0329 | [-0.0684, +0.0235] | -0.1547 | +0.8562 | 0.5000 |
| wtc11_i3neuron | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 2 | +0.2347 | [+0.2347, +0.3224] | +0.0288 | [+0.0288, +0.0544] | +0.3683 | +0.8562 | 1.0000 |
| wtc11_i3neuron | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 2 | -0.0619 | [-0.3131, -0.0531] | -0.0083 | [-0.0679, +0.0218] | -0.0638 | +0.8562 | 1.0000 |
| wtc11_i3neuron | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 2 | -0.2284 | [-0.3453, -0.0071] | -0.0104 | [-0.0276, +0.0198] | -0.0717 | +0.8562 | 1.0000 |
| wtc11_i3neuron | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 2 | -0.1127 | [-0.2492, +0.0938] | -0.0142 | [-0.0665, +0.0447] | -0.1147 | +0.8562 | 1.0000 |
| wtc11_i3neuron | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 2 | +0.1313 | [+0.0687, +0.3130] | +0.0168 | [+0.0031, +0.0737] | +0.0996 | +0.8562 | 0.5000 |
| wtc11_i3neuron | Secondary | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 2 | -0.0915 | [-0.2308, +0.0369] | +0.0222 | [+0.0105, +0.0222] | +0.0734 | +0.8562 | 1.0000 |
| wtc11_i3neuron | Secondary | `NegCtrl_FCT_Perm0` | 2 | -0.0647 | [-0.0660, +0.0429] | +0.2285 | [+0.2285, +0.5362] | +0.4206 | +0.8562 | 0.5000 |
| wtc11_i3neuron | Secondary | `FCT_FullNTMean` | 2 | -0.0313 | [-0.0334, +0.2649] | -0.0329 | [-0.0684, +0.0235] | -0.1547 | +0.8562 | 0.5000 |
| wtc11_i3neuron | Secondary | `Transfer_Mean_GlobalScale` | 2 | +nan | [+nan, +nan] | -0.0110 | [-0.0221, -0.0110] | -0.5286 | +0.8562 | 1.0000 |
| wtc11_i3neuron | Secondary | `Transfer_Mean_VarScaled` | 2 | -0.1069 | [-0.1069, -0.0068] | -0.0584 | [-0.0834, -0.0584] | -0.3205 | +0.8562 | 0.5000 |
| wtc11_i3neuron | Secondary | `CIPHER_Linear` | 2 | -0.5754 | [-0.6267, +0.2638] | +0.2670 | [+0.2670, +0.8812] | +0.3089 | +0.8562 | 0.5000 |
| wtc11_i3neuron | Secondary | `Hybrid_ProjShared` | 2 | -0.5066 | [-0.5201, +0.1724] | -0.1167 | [-0.1181, +0.0587] | -0.3388 | +0.8562 | 0.5000 |
| wtc11_i3neuron | Secondary | `Hybrid_on_TMVS` | 2 | -0.4892 | [-0.4892, +0.1240] | -0.1750 | [-0.1877, -0.0247] | -0.4516 | +0.8562 | 0.5000 |

## 100-Perm `NegCtrl_FCT_Perm` (`k=0..99`) & Mismatched Mean Summary

- **calu3**: perm100 `r_orth_cal` median=+0.0776, p97.5=+0.1519, p98.75=+0.1596; max mismatched hi=+0.2702
- **cd4t_rest**: perm100 `r_orth_cal` median=+0.0347, p97.5=+0.0554, p98.75=+0.0570; max mismatched hi=+0.0761
- **dld1**: perm100 `r_orth_cal` median=+0.0892, p97.5=+0.1179, p98.75=+0.1187; max mismatched hi=+0.1092
- **wtc11_i3neuron**: perm100 `r_orth_cal` median=-0.0219, p97.5=+0.1681, p98.75=+0.1796; max mismatched hi=+0.3224

## `SA-3` Round 9 Anchors Notice (`P7`) & `SA-6 / S7` Cell-Count Controls (`P6`)

- *Benchmark Anchor Note (`P7`):* `Hybrid_ProjShared` follows frozen Round 9 commit `1ca9f60f` 3D `span{1, b_S, a_S}` projection (whereas Round 9 `PREREG_A v4` prose described rank-1 `q_1`).
- **calu3 S7**: Spearman(`r_orth_p`, ln n_cells)=-0.0952 [-0.9747, +0.9259], Spearman(`r_orth_p`, |k_p|)=-0.0952 [-0.8519, +0.9753], OLS-adjusted intercept=+0.1003 [-0.1496, +0.3616]
- **cd4t_rest S7**: Spearman(`r_orth_p`, ln n_cells)=-0.1428 [-0.2637, -0.0167], Spearman(`r_orth_p`, |k_p|)=+0.0828 [-0.0454, +0.2072], OLS-adjusted intercept=+0.0607 [+0.0288, +0.0920]
  * `T1_low` (N=101, n_cells=[25..86]): FCT `r_orth_cal`=+0.0988 [+0.0443, +0.1529], max_mis_hi=+0.1203, perm100_p975=+0.0745, HC1=True, HC2=False
  * `T2_mid` (N=101, n_cells=[87..181]): FCT `r_orth_cal`=-0.0030 [-0.0967, +0.0797], max_mis_hi=+0.1174, perm100_p975=+0.0316, HC1=False, HC2=False
  * `T3_high` (N=100, n_cells=[181..1786]): FCT `r_orth_cal`=+0.0369 [-0.0310, +0.1082], max_mis_hi=+0.1257, perm100_p975=+0.0630, HC1=False, HC2=False
- **dld1 S7**: Spearman(`r_orth_p`, ln n_cells)=-0.0005 [-0.2382, +0.2399], Spearman(`r_orth_p`, |k_p|)=+0.1322 [-0.1066, +0.3649], OLS-adjusted intercept=+0.1400 [+0.0978, +0.1829]
  * `T1_low` (N=29, n_cells=[25..33]): FCT `r_orth_cal`=+0.1403 [+0.0654, +0.2137], max_mis_hi=+0.1538, perm100_p975=+0.1339, HC1=True, HC2=False
  * `T2_mid` (N=29, n_cells=[34..47]): FCT `r_orth_cal`=+0.1456 [+0.0434, +0.2384], max_mis_hi=+0.1567, perm100_p975=+0.1277, HC1=True, HC2=False
  * `T3_high` (N=28, n_cells=[48..89]): FCT `r_orth_cal`=+0.1290 [+0.0171, +0.2282], max_mis_hi=+0.1763, perm100_p975=+0.1284, HC1=True, HC2=False
- **wtc11_i3neuron S7**: Spearman(`r_orth_p`, ln n_cells)=+nan [+nan, +nan], Spearman(`r_orth_p`, |k_p|)=+nan [+nan, +nan], OLS-adjusted intercept=+nan [+nan, +nan]
