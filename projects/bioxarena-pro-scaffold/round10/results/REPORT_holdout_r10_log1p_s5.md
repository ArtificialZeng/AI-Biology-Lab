# Round 10 Confirmatory Group C (`FCT`) Report (`holdout_r10_log1p_s5`)

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
| calu3 | Secondary | `Source_Grand_Mean` | 8 | -0.1235 | [-0.4381, +0.2563] | +0.0040 | [-0.0598, +0.0814] | +0.0092 | +0.8777 | 0.5000 |
| calu3 | Secondary | `Transfer_Mean` | 8 | +nan | [+nan, +nan] | +0.0000 | [+0.0000, +0.0000] | +nan | +0.8777 | 0.6429 |
| calu3 | Secondary | `FCT` | 8 | +0.0897 | [-0.0721, +0.2768] | +0.0066 | [-0.0185, +0.0395] | +0.0374 | +0.8777 | 0.6964 |
| calu3 | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 8 | +0.0799 | [+0.0226, +0.1486] | +0.0076 | [+0.0007, +0.0164] | +0.0843 | +0.8777 | 0.6429 |
| calu3 | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 8 | -0.0659 | [-0.1408, +0.0037] | -0.0049 | [-0.0190, +0.0090] | -0.0437 | +0.8777 | 0.6607 |
| calu3 | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 8 | +0.0101 | [-0.1351, +0.1420] | +0.0070 | [-0.0077, +0.0272] | +0.0610 | +0.8777 | 0.6429 |
| calu3 | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 8 | +0.0361 | [-0.1133, +0.1956] | +0.0078 | [-0.0074, +0.0248] | +0.0778 | +0.8777 | 0.6250 |
| calu3 | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 8 | +0.0365 | [-0.1087, +0.1570] | +0.0023 | [-0.0051, +0.0126] | +0.0340 | +0.8777 | 0.6250 |
| calu3 | Secondary | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 8 | +0.1134 | [-0.0992, +0.2704] | +0.0216 | [-0.0057, +0.0495] | +0.1186 | +0.8777 | 0.7143 |
| calu3 | Secondary | `NegCtrl_FCT_Perm0` | 8 | +0.1376 | [-0.0315, +0.2845] | +0.0420 | [+0.0086, +0.0847] | +0.1286 | +0.8777 | 0.7321 |
| calu3 | Secondary | `FCT_FullNTMean` | 8 | +0.0897 | [-0.0721, +0.2768] | +0.0066 | [-0.0185, +0.0395] | +0.0374 | +0.8777 | 0.6964 |
| calu3 | Secondary | `Transfer_Mean_GlobalScale` | 8 | +nan | [+nan, +nan] | -0.0062 | [-0.0171, +0.0022] | -0.0817 | +0.8777 | 0.6429 |
| calu3 | Secondary | `Transfer_Mean_VarScaled` | 8 | -0.0363 | [-0.1799, +0.1212] | -0.0513 | [-0.1381, +0.0301] | -0.1000 | +0.8777 | 0.6250 |
| calu3 | Secondary | `CIPHER_Linear` | 8 | +0.2971 | [-0.3634, +0.7573] | +0.2260 | [-0.1512, +0.5963] | +0.3051 | +0.8777 | 0.5893 |
| calu3 | Secondary | `Hybrid_ProjShared` | 8 | +0.2583 | [-0.3074, +0.6849] | +0.1288 | [-0.1513, +0.3706] | +0.2510 | +0.8777 | 0.6607 |
| calu3 | Secondary | `Hybrid_on_TMVS` | 8 | +0.2039 | [-0.3306, +0.6113] | +0.0774 | [-0.2447, +0.3632] | +0.1028 | +0.8777 | 0.6786 |
| cd4t_rest | Confirmatory | `Source_Grand_Mean` | 302 | -0.0002 | [-0.0865, +0.0885] | +0.0723 | [+0.0378, +0.1118] | +0.1913 | +0.6483 | 0.5000 |
| cd4t_rest | Confirmatory | `Transfer_Mean` | 302 | +nan | [+nan, +nan] | +0.0000 | [+0.0000, +0.0000] | +nan | +0.6483 | 0.6378 |
| cd4t_rest | Confirmatory | `FCT` | 302 | +0.0543 | [+0.0151, +0.0918] | +0.0151 | [+0.0022, +0.0304] | +0.0778 | +0.6483 | 0.6378 |
| cd4t_rest | Confirmatory | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 302 | +0.0163 | [-0.0162, +0.0508] | +0.0085 | [+0.0036, +0.0139] | +0.0829 | +0.6483 | 0.6372 |
| cd4t_rest | Confirmatory | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 302 | +0.0341 | [-0.0024, +0.0697] | -0.0012 | [-0.0065, +0.0046] | -0.0111 | +0.6483 | 0.6284 |
| cd4t_rest | Confirmatory | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 302 | +0.0152 | [-0.0136, +0.0460] | +0.0145 | [+0.0088, +0.0204] | +0.1136 | +0.6483 | 0.6355 |
| cd4t_rest | Confirmatory | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 302 | +0.0297 | [-0.0037, +0.0635] | +0.0077 | [+0.0034, +0.0124] | +0.0743 | +0.6483 | 0.6290 |
| cd4t_rest | Confirmatory | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 302 | +0.0318 | [-0.0102, +0.0688] | +0.0025 | [-0.0011, +0.0068] | +0.0324 | +0.6483 | 0.6353 |
| cd4t_rest | Confirmatory | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 302 | +0.0270 | [-0.0126, +0.0659] | +0.0237 | [+0.0033, +0.0505] | +0.1099 | +0.6483 | 0.6314 |
| cd4t_rest | Confirmatory | `NegCtrl_FCT_Perm0` | 302 | +0.0339 | [+0.0096, +0.0583] | +0.0634 | [+0.0388, +0.0936] | +0.2406 | +0.6483 | 0.6418 |
| cd4t_rest | Confirmatory | `FCT_FullNTMean` | 302 | +0.0543 | [+0.0151, +0.0918] | +0.0151 | [+0.0022, +0.0304] | +0.0778 | +0.6483 | 0.6378 |
| cd4t_rest | Confirmatory | `Transfer_Mean_GlobalScale` | 302 | +nan | [+nan, +nan] | -0.0451 | [-0.0636, -0.0281] | -0.2621 | +0.6483 | 0.6378 |
| cd4t_rest | Confirmatory | `Transfer_Mean_VarScaled` | 302 | +0.0202 | [-0.0126, +0.0544] | -0.0906 | [-0.1360, -0.0505] | -0.2464 | +0.6483 | 0.6392 |
| cd4t_rest | Confirmatory | `CIPHER_Linear` | 302 | -0.0355 | [-0.1406, +0.0652] | +0.0946 | [+0.0484, +0.1441] | +0.2168 | +0.6483 | 0.4869 |
| cd4t_rest | Confirmatory | `Hybrid_ProjShared` | 302 | -0.0274 | [-0.1050, +0.0474] | -0.0052 | [-0.0166, +0.0057] | -0.0315 | +0.6483 | 0.6197 |
| cd4t_rest | Confirmatory | `Hybrid_on_TMVS` | 302 | -0.0063 | [-0.0656, +0.0518] | -0.0958 | [-0.1418, -0.0554] | -0.2295 | +0.6483 | 0.6263 |
| dld1 | Confirmatory | `Source_Grand_Mean` | 86 | +0.2216 | [+0.1122, +0.3287] | -0.0572 | [-0.1102, -0.0075] | -0.1087 | +0.8275 | 0.5000 |
| dld1 | Confirmatory | `Transfer_Mean` | 86 | +nan | [+nan, +nan] | +0.0000 | [+0.0000, +0.0000] | +nan | +0.8275 | 0.9023 |
| dld1 | Confirmatory | `FCT` | 86 | +0.1433 | [+0.0888, +0.1971] | +0.0294 | [+0.0198, +0.0400] | +0.2111 | +0.8275 | 0.8969 |
| dld1 | Confirmatory | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 86 | +0.0254 | [-0.0216, +0.0694] | +0.0023 | [-0.0036, +0.0080] | +0.0179 | +0.8275 | 0.9014 |
| dld1 | Confirmatory | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 86 | -0.0221 | [-0.0737, +0.0303] | -0.0007 | [-0.0066, +0.0060] | -0.0055 | +0.8275 | 0.8907 |
| dld1 | Confirmatory | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 86 | +0.0307 | [-0.0132, +0.0765] | -0.0047 | [-0.0122, +0.0033] | -0.0305 | +0.8275 | 0.8967 |
| dld1 | Confirmatory | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 86 | -0.0079 | [-0.0574, +0.0427] | -0.0052 | [-0.0114, +0.0009] | -0.0453 | +0.8275 | 0.8911 |
| dld1 | Confirmatory | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 86 | +0.0474 | [-0.0186, +0.1040] | +0.0091 | [+0.0016, +0.0173] | +0.0840 | +0.8275 | 0.8999 |
| dld1 | Confirmatory | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 86 | +0.0696 | [+0.0318, +0.1081] | -0.0048 | [-0.0207, +0.0134] | -0.0167 | +0.8275 | 0.8967 |
| dld1 | Confirmatory | `NegCtrl_FCT_Perm0` | 86 | +0.1154 | [+0.0761, +0.1509] | -0.0542 | [-0.0825, -0.0275] | -0.1345 | +0.8275 | 0.9059 |
| dld1 | Confirmatory | `FCT_FullNTMean` | 86 | +0.1433 | [+0.0888, +0.1971] | +0.0294 | [+0.0198, +0.0400] | +0.2111 | +0.8275 | 0.8969 |
| dld1 | Confirmatory | `Transfer_Mean_GlobalScale` | 86 | +nan | [+nan, +nan] | +0.0368 | [+0.0247, +0.0501] | +0.2397 | +0.8275 | 0.9023 |
| dld1 | Confirmatory | `Transfer_Mean_VarScaled` | 86 | +0.0933 | [+0.0279, +0.1571] | +0.0540 | [+0.0314, +0.0775] | +0.2164 | +0.8275 | 0.9047 |
| dld1 | Confirmatory | `CIPHER_Linear` | 86 | -0.0316 | [-0.1051, +0.0542] | -0.1585 | [-0.2140, -0.1078] | -0.2434 | +0.8275 | 0.4707 |
| dld1 | Confirmatory | `Hybrid_ProjShared` | 86 | -0.0424 | [-0.1470, +0.0695] | -0.0053 | [-0.0202, +0.0053] | -0.0456 | +0.8275 | 0.8886 |
| dld1 | Confirmatory | `Hybrid_on_TMVS` | 86 | +0.0402 | [-0.0640, +0.1407] | +0.0488 | [+0.0202, +0.0762] | +0.1741 | +0.8275 | 0.8943 |
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

- **calu3**: perm100 `r_orth_cal` median=+0.0779, p97.5=+0.1524, p98.75=+0.1599; max mismatched hi=+0.2704
- **cd4t_rest**: perm100 `r_orth_cal` median=+0.0297, p97.5=+0.0514, p98.75=+0.0519; max mismatched hi=+0.0697
- **dld1**: perm100 `r_orth_cal` median=+0.0874, p97.5=+0.1187, p98.75=+0.1193; max mismatched hi=+0.1081
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
