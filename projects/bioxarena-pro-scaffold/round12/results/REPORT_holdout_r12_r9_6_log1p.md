# Round 12 Prospective Group E Report (`holdout_r12_r9_6_log1p`)

- **Phase**: `holdout_r12` | **Source Mode**: `r9_6` (`6` sources) | **PREREG_E Commit**: `4c606f4732ea080768e27757c689d40d7892cd3d` | **Code Freeze**: `11d91c72025264b838af091abcf0eac71588eb91`
- **Bootstrap**: `np.random.default_rng(20261201)` x `4000` (Bonferroni 2-claim + IUT two-sided `97.5%` CIs `[1.25%, 98.75%]`)
- **Cohort Disclosure**: **K_conf = 1 confirmatory line** (`['a549_ko']`) — confirmatory inference permitted under `PREREG_E §1`; **mandatory disclosure: single-line generalizability is limited.**
- **Secondary / Non-Confirmatory Lines (Reporting-Only)**: `['a549_ko_cas9pos', 'hek293_xucao']`
- **Confirmatory Claim 1 (`H-E1` Target-Line Specificity)**: **CONFIRMED**
- **Confirmatory Claim 2 (`H-E2` Persistent `SA-8` Constant-Baseline Inflation)**: **CONFIRMED**

## Confirmatory Hypothesis Adjudication (`97.5%` Two-Sided CIs via IUT + 2-Claim Bonferroni)

| Claim / Hypothesis | Pass / Total | Criterion | Verdict |
|---|---:|---|---|
| **Claim 1 (`H-E1`)** | 1/1 | `H-E1a ∩ H-E1b ∩ H-E1c` (Unanimous) | **CONFIRMED** |
| ├─ `H-E1a` (`FCT` all_nt 97.5% lo > 0, valid >= 0.95, rho_orth > 0) | 1/1 | Addendum A1 | **CONFIRMED** |
| ├─ `H-E1b` (Paired `FCT − Mis_s` 97.5% lo > 0, ALL sources) | 1/1 | Source-wise IUT | **CONFIRMED** |
| └─ `H-E1c` (`FCT` 97.5% lo > Perm100 `all_nontarget` Q0.975) | 1/1 | Perm null | **CONFIRMED** |
| **Claim 2 (`H-E2`)** | 1/1 | `H-E2a ∩ H-E2b` (Unanimous) | **CONFIRMED** |
| ├─ `H-E2a` (`E1_ConstMed` under `resp_top50_split` 97.5% lo > 0) | 1/1 | SA-8 inflation | **CONFIRMED** |
| └─ `H-E2b` (Joint `E1(resp_top50_split) − E1(all_nontarget)` 97.5% lo > 0) | 1/1 | Common-pert IUT | **CONFIRMED** |

## `H-E1b`: Paired Specificity Contrasts `FCT − NegCtrl_FCT_MismatchedMean(s)` (`all_nontarget`, 97.5% CI)

| Target Line | Role | Mismatched Source | Estimate | 97.5% CI | `lo > 0` |
|---|---|---|---:|---|---|
| a549_ko | Confirmatory | `nadig_hepg2` | +0.1869 | [+0.1693, +0.2049] | **PASS** |
| a549_ko | Confirmatory | `nadig_jurkat` | +0.0796 | [+0.0693, +0.0900] | **PASS** |
| a549_ko | Confirmatory | `nourreddine_kolf21j` | +0.1069 | [+0.0985, +0.1156] | **PASS** |
| a549_ko | Confirmatory | `replogle_k562_essential` | +0.1279 | [+0.1150, +0.1410] | **PASS** |
| a549_ko | Confirmatory | `replogle_rpe1` | +0.0804 | [+0.0682, +0.0917] | **PASS** |
| a549_ko | Confirmatory | `vcc2025_h1` | +0.1505 | [+0.1402, +0.1605] | **PASS** |
| hek293_xucao | Secondary | `nadig_hepg2` | +0.1914 | [+0.1013, +0.2945] | **PASS** |
| hek293_xucao | Secondary | `nadig_jurkat` | +0.1914 | [+0.1223, +0.2607] | **PASS** |
| hek293_xucao | Secondary | `nourreddine_kolf21j` | +0.1110 | [+0.0666, +0.1588] | **PASS** |
| hek293_xucao | Secondary | `replogle_k562_essential` | +0.1611 | [+0.1047, +0.2234] | **PASS** |
| hek293_xucao | Secondary | `replogle_rpe1` | +0.1997 | [+0.1139, +0.2923] | **PASS** |
| hek293_xucao | Secondary | `vcc2025_h1` | +0.1798 | [+0.1043, +0.2601] | **PASS** |
| a549_ko_cas9pos | Secondary | `nadig_hepg2` | +0.1907 | [+0.1737, +0.2096] | **PASS** |
| a549_ko_cas9pos | Secondary | `nadig_jurkat` | +0.0849 | [+0.0739, +0.0957] | **PASS** |
| a549_ko_cas9pos | Secondary | `nourreddine_kolf21j` | +0.1112 | [+0.1028, +0.1200] | **PASS** |
| a549_ko_cas9pos | Secondary | `replogle_k562_essential` | +0.1330 | [+0.1193, +0.1464] | **PASS** |
| a549_ko_cas9pos | Secondary | `replogle_rpe1` | +0.0835 | [+0.0719, +0.0959] | **PASS** |
| a549_ko_cas9pos | Secondary | `vcc2025_h1` | +0.1563 | [+0.1465, +0.1665] | **PASS** |

## `H-E2` / `SA-8`: Split-Half Top-50 (`resp_top50_split`) & Joint Cross-Mask `E1(split) − E1(all)`

| Target Line | Role | N_common | `E1(all_common)` | `E1(split)` [97.5% CI] | Joint `E1(split) − E1(all)` [97.5% CI] | `FCT(split)` [97.5% CI] | `Perm100 split` Q97.5 |
|---|---|---:|---:|---|---|---|---:|
| a549_ko | Confirmatory | 445 | -0.0201 | +0.2565 [+0.2195, +0.2981] | +0.2766 [+0.2412, +0.3164] | +0.2965 [+0.2633, +0.3310] | +0.2108 |
| a549_ko_cas9pos | Secondary | 445 | -0.0213 | +0.2523 [+0.2172, +0.2918] | +0.2735 [+0.2398, +0.3112] | +0.3013 [+0.2701, +0.3376] | +0.2137 |
| hek293_xucao | Secondary | 16 | -0.0182 | +0.0157 [-0.0988, +0.1350] | +0.0339 [-0.0868, +0.1797] | +0.2074 [-0.0221, +0.3985] | +0.1301 |

## Reporting-Only Legacy Contrasts (Unpaired `H-D2` & Paired `FCT − E1_ConstMed`)

| Target Line | Role | Unpaired `FCT_lo − max(Mis_hi)` | `FCT − E1 (all_nt)` [97.5% CI] | `FCT − E1 (resp_top50)` [97.5% CI] |
|---|---|---:|---|---|
| a549_ko | Confirmatory | +0.0504 (PASS) | +0.2099 [+0.1975, +0.2223] | +0.0207 [-0.0078, +0.0516] |
| a549_ko_cas9pos | Secondary | +0.0547 (PASS) | +0.2163 [+0.2045, +0.2283] | +0.0356 [+0.0071, +0.0649] |
| hek293_xucao | Secondary | -0.0487 (FAIL) | +0.1743 [+0.0853, +0.2653] | +0.2124 [+0.0917, +0.3095] |

## Reporting-Only Full Line Metrics (`all_nontarget` & Naive `resp_top50`, 97.5% CIs)

| Target Line | Role | Method | N | `all_nontarget r_orth_cal` | 97.5% CI | Naive `resp_top50 r_orth_cal` | `rho_AB_orth` | `valid_frac` |
|---|---|---|---:|---:|---|---:|---:|---:|
| a549_ko | Confirmatory | `Source_Grand_Mean` | 445 | +0.4051 | [+0.3814, +0.4287] | +0.4520 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `Transfer_Mean` | 445 | +nan | [+nan, +nan] | +nan | +0.7887 | 0.0000 |
| a549_ko | Confirmatory | `FCT` | 445 | +0.1898 | [+0.1742, +0.2055] | +0.3073 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 445 | +0.0029 | [-0.0053, +0.0114] | +0.0608 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 445 | +0.1102 | [+0.0977, +0.1237] | +0.0957 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 445 | +0.0829 | [+0.0734, +0.0927] | +0.1139 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 445 | +0.0619 | [+0.0519, +0.0723] | -0.0184 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 445 | +0.1094 | [+0.0966, +0.1231] | +0.1284 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 445 | +0.0393 | [+0.0303, +0.0483] | +0.2244 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `NegCtrl_FCT_Perm0` | 445 | -0.0012 | [-0.0080, +0.0057] | +0.2068 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `FCT_FullNTMean` | 445 | +0.1898 | [+0.1742, +0.2055] | +0.3073 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `Transfer_Mean_GlobalScale` | 445 | +nan | [+nan, +nan] | +nan | +0.7887 | 0.0000 |
| a549_ko | Confirmatory | `Transfer_Mean_VarScaled` | 445 | +0.2217 | [+0.2071, +0.2359] | +0.1871 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `CIPHER_Linear` | 445 | +0.2817 | [+0.2408, +0.3238] | +0.3128 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `Hybrid_ProjShared` | 445 | +0.1554 | [+0.1303, +0.1820] | +0.2307 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `Hybrid_on_TMVS` | 445 | +0.2373 | [+0.2127, +0.2625] | +0.2767 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `E1_ConstMed` | 445 | -0.0201 | [-0.0293, -0.0108] | +0.2866 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `E2_ConstP10` | 445 | -0.0229 | [-0.0318, -0.0137] | +0.2412 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `E2_ConstP90` | 445 | -0.0141 | [-0.0239, -0.0042] | +0.3257 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `E3_SrcMeanCtrl` | 445 | +0.0794 | [+0.0674, +0.0916] | +0.1845 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `E5_RankMatchedMu` | 445 | +0.1022 | [+0.0912, +0.1142] | +0.1224 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `E6_LogFCBar` | 445 | -0.0056 | [-0.0163, +0.0048] | +0.3167 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `E7_FCBarMinus1` | 445 | -0.0241 | [-0.0329, -0.0150] | +0.2009 | +0.7887 | 1.0000 |
| hek293_xucao | Secondary | `Source_Grand_Mean` | 16 | -0.1415 | [-0.2391, -0.0289] | -0.1147 | +0.3005 | 1.0000 |
| hek293_xucao | Secondary | `Transfer_Mean` | 16 | +nan | [+nan, +nan] | +nan | +0.3005 | 0.0000 |
| hek293_xucao | Secondary | `FCT` | 16 | +0.1561 | [+0.0595, +0.2593] | +0.1722 | +0.3005 | 1.0000 |
| hek293_xucao | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 16 | -0.0354 | [-0.0630, -0.0036] | -0.1246 | +0.3005 | 1.0000 |
| hek293_xucao | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 16 | -0.0354 | [-0.0796, +0.0117] | -0.0439 | +0.3005 | 1.0000 |
| hek293_xucao | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 16 | +0.0450 | [-0.0139, +0.1083] | -0.0876 | +0.3005 | 1.0000 |
| hek293_xucao | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 16 | -0.0050 | [-0.0532, +0.0457] | -0.0731 | +0.3005 | 1.0000 |
| hek293_xucao | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 16 | -0.0437 | [-0.0657, -0.0099] | +0.0095 | +0.3005 | 1.0000 |
| hek293_xucao | Secondary | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 16 | -0.0237 | [-0.0769, +0.0199] | -0.1066 | +0.3005 | 1.0000 |
| hek293_xucao | Secondary | `NegCtrl_FCT_Perm0` | 16 | -0.0119 | [-0.0429, +0.0130] | -0.0686 | +0.3005 | 1.0000 |
| hek293_xucao | Secondary | `FCT_FullNTMean` | 16 | +0.1561 | [+0.0595, +0.2593] | +0.1722 | +0.3005 | 1.0000 |
| hek293_xucao | Secondary | `Transfer_Mean_GlobalScale` | 16 | +nan | [+nan, +nan] | +nan | +0.3005 | 0.0000 |
| hek293_xucao | Secondary | `Transfer_Mean_VarScaled` | 16 | +0.1972 | [+0.1163, +0.2841] | +0.0216 | +0.3005 | 1.0000 |
| hek293_xucao | Secondary | `CIPHER_Linear` | 16 | -0.1292 | [-0.2175, -0.0162] | -0.1667 | +0.3005 | 1.0000 |
| hek293_xucao | Secondary | `Hybrid_ProjShared` | 16 | -0.0395 | [-0.0849, -0.0009] | -0.0861 | +0.3005 | 1.0000 |
| hek293_xucao | Secondary | `Hybrid_on_TMVS` | 16 | +0.1633 | [+0.0807, +0.2503] | +0.0010 | +0.3005 | 1.0000 |
| hek293_xucao | Secondary | `E1_ConstMed` | 16 | -0.0182 | [-0.0582, +0.0148] | -0.0402 | +0.3005 | 1.0000 |
| hek293_xucao | Secondary | `E2_ConstP10` | 16 | -0.0192 | [-0.0568, +0.0134] | -0.0734 | +0.3005 | 1.0000 |
| hek293_xucao | Secondary | `E2_ConstP90` | 16 | -0.0182 | [-0.0599, +0.0150] | +0.0041 | +0.3005 | 1.0000 |
| hek293_xucao | Secondary | `E3_SrcMeanCtrl` | 16 | -0.0197 | [-0.0670, +0.0323] | -0.1419 | +0.3005 | 1.0000 |
| hek293_xucao | Secondary | `E5_RankMatchedMu` | 16 | -0.0139 | [-0.0619, +0.0387] | -0.1440 | +0.3005 | 1.0000 |
| hek293_xucao | Secondary | `E6_LogFCBar` | 16 | -0.0172 | [-0.0575, +0.0142] | +0.0474 | +0.3005 | 1.0000 |
| hek293_xucao | Secondary | `E7_FCBarMinus1` | 16 | -0.0209 | [-0.0560, +0.0103] | -0.0909 | +0.3005 | 1.0000 |
| a549_ko_cas9pos | Secondary | `Source_Grand_Mean` | 445 | +0.4101 | [+0.3873, +0.4327] | +0.4615 | +0.7388 | 1.0000 |
| a549_ko_cas9pos | Secondary | `Transfer_Mean` | 445 | +nan | [+nan, +nan] | +nan | +0.7388 | 0.0000 |
| a549_ko_cas9pos | Secondary | `FCT` | 445 | +0.1950 | [+0.1797, +0.2108] | +0.3192 | +0.7388 | 1.0000 |
| a549_ko_cas9pos | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 445 | +0.0043 | [-0.0039, +0.0129] | +0.0579 | +0.7388 | 1.0000 |
| a549_ko_cas9pos | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 445 | +0.1102 | [+0.0976, +0.1225] | +0.0983 | +0.7388 | 1.0000 |
| a549_ko_cas9pos | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 445 | +0.0838 | [+0.0745, +0.0933] | +0.1155 | +0.7388 | 1.0000 |
| a549_ko_cas9pos | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 445 | +0.0620 | [+0.0523, +0.0722] | -0.0181 | +0.7388 | 1.0000 |
| a549_ko_cas9pos | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 445 | +0.1115 | [+0.0987, +0.1249] | +0.1337 | +0.7388 | 1.0000 |
| a549_ko_cas9pos | Secondary | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 445 | +0.0387 | [+0.0298, +0.0476] | +0.2280 | +0.7388 | 1.0000 |
| a549_ko_cas9pos | Secondary | `NegCtrl_FCT_Perm0` | 445 | -0.0017 | [-0.0086, +0.0051] | +0.2061 | +0.7388 | 1.0000 |
| a549_ko_cas9pos | Secondary | `FCT_FullNTMean` | 445 | +0.1950 | [+0.1797, +0.2108] | +0.3192 | +0.7388 | 1.0000 |
| a549_ko_cas9pos | Secondary | `Transfer_Mean_GlobalScale` | 445 | +nan | [+nan, +nan] | +nan | +0.7388 | 0.0000 |
| a549_ko_cas9pos | Secondary | `Transfer_Mean_VarScaled` | 445 | +0.2297 | [+0.2152, +0.2443] | +0.1829 | +0.7388 | 1.0000 |
| a549_ko_cas9pos | Secondary | `CIPHER_Linear` | 445 | +0.2979 | [+0.2556, +0.3394] | +0.3273 | +0.7388 | 1.0000 |
| a549_ko_cas9pos | Secondary | `Hybrid_ProjShared` | 445 | +0.1627 | [+0.1373, +0.1874] | +0.2416 | +0.7388 | 1.0000 |
| a549_ko_cas9pos | Secondary | `Hybrid_on_TMVS` | 445 | +0.2470 | [+0.2236, +0.2705] | +0.2812 | +0.7388 | 1.0000 |
| a549_ko_cas9pos | Secondary | `E1_ConstMed` | 445 | -0.0213 | [-0.0308, -0.0119] | +0.2836 | +0.7388 | 1.0000 |
| a549_ko_cas9pos | Secondary | `E2_ConstP10` | 445 | -0.0239 | [-0.0330, -0.0148] | +0.2372 | +0.7388 | 1.0000 |
| a549_ko_cas9pos | Secondary | `E2_ConstP90` | 445 | -0.0155 | [-0.0254, -0.0057] | +0.3240 | +0.7388 | 1.0000 |
| a549_ko_cas9pos | Secondary | `E3_SrcMeanCtrl` | 445 | +0.0794 | [+0.0677, +0.0912] | +0.1843 | +0.7388 | 1.0000 |
| a549_ko_cas9pos | Secondary | `E5_RankMatchedMu` | 445 | +0.1026 | [+0.0918, +0.1138] | +0.1242 | +0.7388 | 1.0000 |
| a549_ko_cas9pos | Secondary | `E6_LogFCBar` | 445 | -0.0072 | [-0.0178, +0.0033] | +0.3164 | +0.7388 | 1.0000 |
| a549_ko_cas9pos | Secondary | `E7_FCBarMinus1` | 445 | -0.0250 | [-0.0341, -0.0161] | +0.1965 | +0.7388 | 1.0000 |

## Reporting-Only `N_5` Null-Cone Analysis — `all_nontarget` (`r_cone_sp`, 97.5% CIs)

*Note: `Transfer_Mean`, `Transfer_Mean_GlobalScale`, and `Source_Grand_Mean` are **zero by exact mathematical construction** under `N_5`, and `E1_ConstMed`, `E2_ConstP10`, `E2_ConstP90`, `E6_LogFCBar`, `E7_FCBarMinus1` are **suppressed by cone construction** (`R2_N5 >= 0.99`); only `FCT`, `Mis_s`, single-draw `Perm0..99`, `E3`/`E5`, and non-cone target models are empirically informative.*

| Target Line | Method | `N_5 r_cone_sp` | 97.5% CI | `eta` | `R2_N5` | Structural Role |
|---|---|---:|---|---:|---:|---|
| a549_ko | `Source_Grand_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko | `Transfer_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko | `FCT` | +0.1992 | [+0.1850, +0.2133] | 0.9350 | 0.125788 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | -0.0114 | [-0.0167, -0.0057] | 0.9336 | 0.128458 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | +0.0528 | [+0.0474, +0.0586] | 0.8569 | 0.265742 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | +0.0638 | [+0.0581, +0.0698] | 0.8980 | 0.193679 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | +0.0175 | [+0.0121, +0.0228] | 0.9276 | 0.139574 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | +0.0716 | [+0.0637, +0.0799] | 0.9355 | 0.124827 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | +0.0236 | [+0.0214, +0.0260] | 0.7829 | 0.387123 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_Perm0` | +0.0091 | [+0.0078, +0.0103] | 0.7029 | 0.505865 | `INFORMATIVE` |
| a549_ko | `FCT_FullNTMean` | +0.1992 | [+0.1850, +0.2133] | 0.9350 | 0.125788 | `INFORMATIVE` |
| a549_ko | `Transfer_Mean_GlobalScale` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko | `Transfer_Mean_VarScaled` | +0.1273 | [+0.1190, +0.1354] | 0.5880 | 0.654241 | `INFORMATIVE` |
| a549_ko | `CIPHER_Linear` | +0.0706 | [+0.0581, +0.0840] | 0.3199 | 0.897687 | `INFORMATIVE` |
| a549_ko | `Hybrid_ProjShared` | +0.1709 | [+0.1427, +0.2018] | 0.9926 | 0.014715 | `INFORMATIVE` |
| a549_ko | `Hybrid_on_TMVS` | +0.2049 | [+0.1821, +0.2282] | 0.8251 | 0.319223 | `INFORMATIVE` |
| a549_ko | `E1_ConstMed` | +0.0001 | [+0.0000, +0.0001] | 0.0028 | 0.999992 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko | `E2_ConstP10` | +0.0000 | [+0.0000, +0.0000] | 0.0007 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko | `E2_ConstP90` | -0.0000 | [-0.0000, -0.0000] | 0.0020 | 0.999996 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko | `E3_SrcMeanCtrl` | +0.0440 | [+0.0396, +0.0487] | 0.7678 | 0.410512 | `INFORMATIVE` |
| a549_ko | `E5_RankMatchedMu` | +0.0555 | [+0.0500, +0.0613] | 0.9400 | 0.116366 | `INFORMATIVE` |
| a549_ko | `E6_LogFCBar` | +0.0000 | [-0.0000, +0.0000] | 0.0123 | 0.999849 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko | `E7_FCBarMinus1` | -0.0001 | [-0.0002, -0.0001] | 0.0052 | 0.999973 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko_cas9pos | `Source_Grand_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko_cas9pos | `Transfer_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko_cas9pos | `FCT` | +0.2049 | [+0.1912, +0.2190] | 0.9341 | 0.127517 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | -0.0109 | [-0.0166, -0.0052] | 0.9334 | 0.128821 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | +0.0527 | [+0.0475, +0.0581] | 0.8567 | 0.265983 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | +0.0649 | [+0.0595, +0.0708] | 0.8979 | 0.193829 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | +0.0170 | [+0.0119, +0.0224] | 0.9274 | 0.139837 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | +0.0732 | [+0.0657, +0.0813] | 0.9354 | 0.125006 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | +0.0237 | [+0.0214, +0.0259] | 0.7828 | 0.387231 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_Perm0` | +0.0093 | [+0.0082, +0.0105] | 0.7024 | 0.506572 | `INFORMATIVE` |
| a549_ko_cas9pos | `FCT_FullNTMean` | +0.2049 | [+0.1912, +0.2190] | 0.9341 | 0.127517 | `INFORMATIVE` |
| a549_ko_cas9pos | `Transfer_Mean_GlobalScale` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko_cas9pos | `Transfer_Mean_VarScaled` | +0.1358 | [+0.1275, +0.1443] | 0.6064 | 0.632331 | `INFORMATIVE` |
| a549_ko_cas9pos | `CIPHER_Linear` | +0.0749 | [+0.0624, +0.0875] | 0.3196 | 0.897868 | `INFORMATIVE` |
| a549_ko_cas9pos | `Hybrid_ProjShared` | +0.1805 | [+0.1519, +0.2082] | 0.9927 | 0.014459 | `INFORMATIVE` |
| a549_ko_cas9pos | `Hybrid_on_TMVS` | +0.2171 | [+0.1949, +0.2389] | 0.8362 | 0.300795 | `INFORMATIVE` |
| a549_ko_cas9pos | `E1_ConstMed` | +0.0001 | [+0.0000, +0.0001] | 0.0028 | 0.999992 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko_cas9pos | `E2_ConstP10` | +0.0000 | [+0.0000, +0.0000] | 0.0007 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko_cas9pos | `E2_ConstP90` | -0.0000 | [-0.0001, -0.0000] | 0.0020 | 0.999996 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko_cas9pos | `E3_SrcMeanCtrl` | +0.0440 | [+0.0398, +0.0483] | 0.7677 | 0.410585 | `INFORMATIVE` |
| a549_ko_cas9pos | `E5_RankMatchedMu` | +0.0557 | [+0.0504, +0.0610] | 0.9383 | 0.119528 | `INFORMATIVE` |
| a549_ko_cas9pos | `E6_LogFCBar` | +0.0000 | [-0.0000, +0.0000] | 0.0123 | 0.999849 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko_cas9pos | `E7_FCBarMinus1` | -0.0001 | [-0.0002, -0.0001] | 0.0052 | 0.999973 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| hek293_xucao | `Source_Grand_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| hek293_xucao | `Transfer_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| hek293_xucao | `FCT` | +0.1966 | [+0.1047, +0.2837] | 0.8523 | 0.273609 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | -0.0217 | [-0.0416, -0.0033] | 0.9657 | 0.067414 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | +0.0091 | [-0.0171, +0.0406] | 0.8840 | 0.218492 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | +0.0826 | [+0.0304, +0.1337] | 0.9140 | 0.164526 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | +0.0293 | [-0.0066, +0.0648] | 0.9464 | 0.104330 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | -0.0074 | [-0.0232, +0.0087] | 0.9521 | 0.093476 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | +0.0068 | [-0.0132, +0.0259] | 0.7791 | 0.393017 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_Perm0` | -0.0027 | [-0.0071, +0.0022] | 0.6759 | 0.543189 | `INFORMATIVE` |
| hek293_xucao | `FCT_FullNTMean` | +0.1966 | [+0.1047, +0.2837] | 0.8523 | 0.273609 | `INFORMATIVE` |
| hek293_xucao | `Transfer_Mean_GlobalScale` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| hek293_xucao | `Transfer_Mean_VarScaled` | +0.1468 | [+0.0688, +0.2281] | 0.6531 | 0.573450 | `INFORMATIVE` |
| hek293_xucao | `CIPHER_Linear` | -0.0299 | [-0.0490, -0.0088] | 0.2130 | 0.954639 | `INFORMATIVE` |
| hek293_xucao | `Hybrid_ProjShared` | -0.0458 | [-0.0874, -0.0089] | 0.9929 | 0.014238 | `INFORMATIVE` |
| hek293_xucao | `Hybrid_on_TMVS` | +0.1267 | [+0.0522, +0.2074] | 0.6888 | 0.525497 | `INFORMATIVE` |
| hek293_xucao | `E1_ConstMed` | +0.0000 | [-0.0000, +0.0000] | 0.0016 | 0.999997 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| hek293_xucao | `E2_ConstP10` | +0.0000 | [-0.0000, +0.0000] | 0.0005 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| hek293_xucao | `E2_ConstP90` | -0.0000 | [-0.0000, +0.0000] | 0.0010 | 0.999999 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| hek293_xucao | `E3_SrcMeanCtrl` | +0.0214 | [-0.0050, +0.0494] | 0.7879 | 0.379155 | `INFORMATIVE` |
| hek293_xucao | `E5_RankMatchedMu` | +0.0256 | [-0.0052, +0.0570] | 0.9303 | 0.134507 | `INFORMATIVE` |
| hek293_xucao | `E6_LogFCBar` | -0.0000 | [-0.0001, +0.0000] | 0.0013 | 0.999998 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| hek293_xucao | `E7_FCBarMinus1` | -0.0000 | [-0.0000, +0.0000] | 0.0029 | 0.999991 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |

## Reporting-Only `N_5` Null-Cone Analysis — `resp_top50_split` (`r_cone_sp`, 97.5% CIs)

*Note: `Transfer_Mean`, `Transfer_Mean_GlobalScale`, and `Source_Grand_Mean` are **zero by exact mathematical construction** under `N_5`, and `E1_ConstMed`, `E2_ConstP10`, `E2_ConstP90`, `E6_LogFCBar`, `E7_FCBarMinus1` are **suppressed by cone construction** (`R2_N5 >= 0.99`); only `FCT`, `Mis_s`, single-draw `Perm0..99`, `E3`/`E5`, and non-cone target models are empirically informative.*

| Target Line | Method | `N_5 r_cone_sp` | 97.5% CI | `eta` | `R2_N5` | Structural Role |
|---|---|---:|---|---:|---:|---|
| a549_ko | `Source_Grand_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko | `Transfer_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko | `FCT` | +0.1301 | [+0.1099, +0.1518] | 0.6330 | 0.599272 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | -0.0374 | [-0.0509, -0.0237] | 0.7452 | 0.444621 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | +0.0170 | [-0.0006, +0.0330] | 0.8191 | 0.328820 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | +0.0495 | [+0.0359, +0.0634] | 0.8178 | 0.331200 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | +0.0023 | [-0.0127, +0.0169] | 0.8053 | 0.351417 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | +0.0041 | [-0.0133, +0.0216] | 0.7021 | 0.507010 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | +0.0438 | [+0.0324, +0.0550] | 0.7157 | 0.487785 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_Perm0` | +0.0432 | [+0.0354, +0.0515] | 0.5730 | 0.671689 | `INFORMATIVE` |
| a549_ko | `FCT_FullNTMean` | +0.1301 | [+0.1099, +0.1518] | 0.6330 | 0.599272 | `INFORMATIVE` |
| a549_ko | `Transfer_Mean_GlobalScale` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko | `Transfer_Mean_VarScaled` | +0.1458 | [+0.1228, +0.1719] | 0.7164 | 0.486631 | `INFORMATIVE` |
| a549_ko | `CIPHER_Linear` | +0.0599 | [+0.0462, +0.0753] | 0.2838 | 0.919480 | `INFORMATIVE` |
| a549_ko | `Hybrid_ProjShared` | +0.1649 | [+0.1275, +0.2024] | 0.8701 | 0.242934 | `INFORMATIVE` |
| a549_ko | `Hybrid_on_TMVS` | +0.2070 | [+0.1726, +0.2424] | 0.8067 | 0.349245 | `INFORMATIVE` |
| a549_ko | `E1_ConstMed` | -0.0001 | [-0.0002, -0.0001] | 0.0015 | 0.999998 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko | `E2_ConstP10` | -0.0000 | [-0.0001, -0.0000] | 0.0004 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko | `E2_ConstP90` | +0.0002 | [+0.0001, +0.0003] | 0.0023 | 0.999995 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko | `E3_SrcMeanCtrl` | +0.0078 | [-0.0065, +0.0220] | 0.6963 | 0.515032 | `INFORMATIVE` |
| a549_ko | `E5_RankMatchedMu` | +0.0122 | [-0.0041, +0.0282] | 0.7856 | 0.382769 | `INFORMATIVE` |
| a549_ko | `E6_LogFCBar` | -0.0001 | [-0.0002, +0.0000] | 0.0026 | 0.999993 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko | `E7_FCBarMinus1` | +0.0004 | [+0.0002, +0.0007] | 0.0073 | 0.999946 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko_cas9pos | `Source_Grand_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko_cas9pos | `Transfer_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko_cas9pos | `FCT` | +0.1352 | [+0.1146, +0.1569] | 0.6298 | 0.603329 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | -0.0446 | [-0.0587, -0.0313] | 0.7420 | 0.449409 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | +0.0190 | [+0.0026, +0.0362] | 0.8136 | 0.337877 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | +0.0499 | [+0.0363, +0.0641] | 0.8152 | 0.335358 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | -0.0032 | [-0.0186, +0.0122] | 0.8073 | 0.348173 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | +0.0057 | [-0.0116, +0.0235] | 0.6950 | 0.516886 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | +0.0465 | [+0.0349, +0.0584] | 0.7121 | 0.492927 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_Perm0` | +0.0456 | [+0.0374, +0.0543] | 0.5702 | 0.674838 | `INFORMATIVE` |
| a549_ko_cas9pos | `FCT_FullNTMean` | +0.1352 | [+0.1146, +0.1569] | 0.6298 | 0.603329 | `INFORMATIVE` |
| a549_ko_cas9pos | `Transfer_Mean_GlobalScale` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko_cas9pos | `Transfer_Mean_VarScaled` | +0.1571 | [+0.1324, +0.1852] | 0.7278 | 0.470293 | `INFORMATIVE` |
| a549_ko_cas9pos | `CIPHER_Linear` | +0.0626 | [+0.0494, +0.0764] | 0.2809 | 0.921109 | `INFORMATIVE` |
| a549_ko_cas9pos | `Hybrid_ProjShared` | +0.1775 | [+0.1414, +0.2132] | 0.8671 | 0.248140 | `INFORMATIVE` |
| a549_ko_cas9pos | `Hybrid_on_TMVS` | +0.2219 | [+0.1889, +0.2555] | 0.8082 | 0.346867 | `INFORMATIVE` |
| a549_ko_cas9pos | `E1_ConstMed` | -0.0001 | [-0.0002, -0.0001] | 0.0015 | 0.999998 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko_cas9pos | `E2_ConstP10` | -0.0000 | [-0.0001, -0.0000] | 0.0005 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko_cas9pos | `E2_ConstP90` | +0.0002 | [+0.0001, +0.0003] | 0.0024 | 0.999994 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko_cas9pos | `E3_SrcMeanCtrl` | +0.0061 | [-0.0087, +0.0200] | 0.6907 | 0.522911 | `INFORMATIVE` |
| a549_ko_cas9pos | `E5_RankMatchedMu` | +0.0107 | [-0.0061, +0.0263] | 0.7772 | 0.395992 | `INFORMATIVE` |
| a549_ko_cas9pos | `E6_LogFCBar` | -0.0001 | [-0.0002, -0.0000] | 0.0026 | 0.999993 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko_cas9pos | `E7_FCBarMinus1` | +0.0004 | [+0.0002, +0.0007] | 0.0072 | 0.999948 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| hek293_xucao | `Source_Grand_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| hek293_xucao | `Transfer_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| hek293_xucao | `FCT` | +0.1382 | [+0.0229, +0.2581] | 0.5267 | 0.722613 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | -0.0056 | [-0.0804, +0.0955] | 0.7187 | 0.483109 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | -0.0678 | [-0.1229, -0.0020] | 0.7752 | 0.399047 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | -0.0218 | [-0.1065, +0.0513] | 0.7014 | 0.507646 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | -0.0353 | [-0.1023, +0.0273] | 0.6822 | 0.533508 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | +0.0383 | [-0.0235, +0.1193] | 0.6658 | 0.556543 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | -0.0166 | [-0.1304, +0.0305] | 0.4280 | 0.816555 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_Perm0` | -0.0122 | [-0.0410, +0.0409] | 0.3126 | 0.902100 | `INFORMATIVE` |
| hek293_xucao | `FCT_FullNTMean` | +0.1382 | [+0.0229, +0.2581] | 0.5267 | 0.722613 | `INFORMATIVE` |
| hek293_xucao | `Transfer_Mean_GlobalScale` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| hek293_xucao | `Transfer_Mean_VarScaled` | -0.0543 | [-0.1422, +0.1116] | 0.3648 | 0.866918 | `INFORMATIVE` |
| hek293_xucao | `CIPHER_Linear` | -0.0111 | [-0.0384, +0.0025] | 0.0984 | 0.990317 | `INFORMATIVE` |
| hek293_xucao | `Hybrid_ProjShared` | -0.1124 | [-0.2139, +0.0048] | 0.9311 | 0.132880 | `INFORMATIVE` |
| hek293_xucao | `Hybrid_on_TMVS` | -0.0650 | [-0.1456, +0.0911] | 0.3725 | 0.861213 | `INFORMATIVE` |
| hek293_xucao | `E1_ConstMed` | +0.0000 | [-0.0001, +0.0001] | 0.0005 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| hek293_xucao | `E2_ConstP10` | -0.0000 | [-0.0000, +0.0000] | 0.0002 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| hek293_xucao | `E2_ConstP90` | -0.0001 | [-0.0003, +0.0001] | 0.0008 | 0.999999 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| hek293_xucao | `E3_SrcMeanCtrl` | -0.0594 | [-0.1247, +0.0143] | 0.5705 | 0.673575 | `INFORMATIVE` |
| hek293_xucao | `E5_RankMatchedMu` | -0.0592 | [-0.1342, +0.0164] | 0.6734 | 0.544460 | `INFORMATIVE` |
| hek293_xucao | `E6_LogFCBar` | +0.0015 | [-0.0001, +0.0049] | 0.0065 | 0.999957 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| hek293_xucao | `E7_FCBarMinus1` | +0.0001 | [-0.0001, +0.0003] | 0.0018 | 0.999997 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |

## 100-Perm Empirical Null (`all_nontarget`, `SA-8`, & Single-Draw `N_5`)

- **a549_ko**: `all_nontarget` med=-0.0106, Q97.5=+0.0018; `SA-8 split` med=+0.1347, Q97.5=+0.2108; `N_5` single-draw med=+0.0004, med|r|=+0.0036, Q97.5=+0.0112 (med_eta=0.6975)
- **a549_ko_cas9pos**: `all_nontarget` med=-0.0114, Q97.5=+0.0004; `SA-8 split` med=+0.1337, Q97.5=+0.2137; `N_5` single-draw med=+0.0005, med|r|=+0.0038, Q97.5=+0.0112 (med_eta=0.6971)
- **hek293_xucao**: `all_nontarget` med=-0.0126, Q97.5=+0.0015; `SA-8 split` med=+0.0091, Q97.5=+0.1301; `N_5` single-draw med=-0.0002, med|r|=+0.0039, Q97.5=+0.0157 (med_eta=0.6807)

## Mandatory Governance, Provenance & Prior-Art Disclosures (`PREREG_E §0.2, §5, §7.2, §8`)

- **Post-Hoc Provenance**: `H-E1` (paired `FCT - Mis_s` on `all_nontarget`) and `H-E2` (`SA-8 resp_top50_split` constant-baseline inflation) were specified in `PREREG_E` after unblinding Rounds 9-11 (`N=12` lines) and evaluated zero-shot on unseen Round 12 lines.
- **Claim 2 Non-Mechanistic Scope**: Confirmation of `H-E2` establishes an **empirical evaluation artifact strictly under `r_orth_cal`** (`resp_top50_split` vs `all_nontarget`) — namely that single-line constant baseline `E1_ConstMed` attains inflated disattenuated residual correlation under magnitude-selected subsets — and **must NOT be extrapolated** as a mechanistic biological claim or applied directly to raw DEG discovery benchmarks without independent validation.
- **Modality Boundary & Single-Line Generalizability**: Source cohorts (`r9_6`, `r9_r10_8`) comprise **CRISPRi** perturbations whereas `a549_ko` (`GSE345058`) is **CRISPRko**. Failure on `a549_ko` is confounded with CRISPRi->CRISPRko modality shift rather than purely same-modality failure; passage constitutes cross-modality generalization evidence. Under `K_conf = 1`, population generalizability across additional lineages remains explicitly limited.
- **Prior-Art & Data Attribution**: Baseline-aware perturbation transfer context acknowledges **Kawsar et al. (2026), `CanPert/CASTER` (`doi:10.64898/2026.11.15.766815`)**; `a549_ko` CRISPRko dataset derives from **Liu et al. (2026), `GSE345058`**; `hek293_xucao` CRISPRi derives from **Xu & Cao (2023), PerturbSci-Kinetics**.
- **Reporting-Only §5 Scope**: Legacy unpaired `H-D2`, paired `FCT - E1_ConstMed`, `N_5` null-cone suites (`all_nontarget` & `resp_top50_split`), `SA-Src8` (`--source-mode r9_r10_8`), `SA-Cas9` (`a549_ko_cas9pos`), naive `resp_top50`, and `H-D_util` are strictly descriptive and carry zero multiplicity weight.
