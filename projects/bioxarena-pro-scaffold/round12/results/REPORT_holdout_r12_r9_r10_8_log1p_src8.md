# Round 12 Prospective Group E Report (`holdout_r12_r9_r10_8_log1p_src8`)

- **Phase**: `holdout_r12` | **Source Mode**: `r9_r10_8` (`8` sources) | **PREREG_E Commit**: `4c606f4732ea080768e27757c689d40d7892cd3d` | **Code Freeze**: `11d91c72025264b838af091abcf0eac71588eb91`
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
| a549_ko | Confirmatory | `cd4t_rest` | +0.1451 | [+0.1332, +0.1575] | **PASS** |
| a549_ko | Confirmatory | `dld1` | +0.0423 | [+0.0372, +0.0480] | **PASS** |
| a549_ko | Confirmatory | `nadig_hepg2` | +0.1445 | [+0.1283, +0.1611] | **PASS** |
| a549_ko | Confirmatory | `nadig_jurkat` | +0.0598 | [+0.0506, +0.0695] | **PASS** |
| a549_ko | Confirmatory | `nourreddine_kolf21j` | +0.0816 | [+0.0734, +0.0901] | **PASS** |
| a549_ko | Confirmatory | `replogle_k562_essential` | +0.0996 | [+0.0873, +0.1122] | **PASS** |
| a549_ko | Confirmatory | `replogle_rpe1` | +0.0745 | [+0.0657, +0.0839] | **PASS** |
| a549_ko | Confirmatory | `vcc2025_h1` | +0.1176 | [+0.1068, +0.1290] | **PASS** |
| hek293_xucao | Secondary | `cd4t_rest` | +0.1784 | [+0.0879, +0.2706] | **PASS** |
| hek293_xucao | Secondary | `dld1` | +0.1430 | [+0.0821, +0.2091] | **PASS** |
| hek293_xucao | Secondary | `nadig_hepg2` | +0.1783 | [+0.0914, +0.2733] | **PASS** |
| hek293_xucao | Secondary | `nadig_jurkat` | +0.1794 | [+0.1122, +0.2471] | **PASS** |
| hek293_xucao | Secondary | `nourreddine_kolf21j` | +0.1060 | [+0.0629, +0.1524] | **PASS** |
| hek293_xucao | Secondary | `replogle_k562_essential` | +0.1501 | [+0.0980, +0.2075] | **PASS** |
| hek293_xucao | Secondary | `replogle_rpe1` | +0.1898 | [+0.1052, +0.2814] | **PASS** |
| hek293_xucao | Secondary | `vcc2025_h1` | +0.1767 | [+0.1016, +0.2549] | **PASS** |
| a549_ko_cas9pos | Secondary | `cd4t_rest` | +0.1503 | [+0.1374, +0.1629] | **PASS** |
| a549_ko_cas9pos | Secondary | `dld1` | +0.0457 | [+0.0399, +0.0515] | **PASS** |
| a549_ko_cas9pos | Secondary | `nadig_hepg2` | +0.1489 | [+0.1323, +0.1658] | **PASS** |
| a549_ko_cas9pos | Secondary | `nadig_jurkat` | +0.0645 | [+0.0550, +0.0741] | **PASS** |
| a549_ko_cas9pos | Secondary | `nourreddine_kolf21j` | +0.0856 | [+0.0766, +0.0946] | **PASS** |
| a549_ko_cas9pos | Secondary | `replogle_k562_essential` | +0.1045 | [+0.0915, +0.1168] | **PASS** |
| a549_ko_cas9pos | Secondary | `replogle_rpe1` | +0.0777 | [+0.0680, +0.0874] | **PASS** |
| a549_ko_cas9pos | Secondary | `vcc2025_h1` | +0.1234 | [+0.1119, +0.1349] | **PASS** |

## `H-E2` / `SA-8`: Split-Half Top-50 (`resp_top50_split`) & Joint Cross-Mask `E1(split) − E1(all)`

| Target Line | Role | N_common | `E1(all_common)` | `E1(split)` [97.5% CI] | Joint `E1(split) − E1(all)` [97.5% CI] | `FCT(split)` [97.5% CI] | `Perm100 split` Q97.5 |
|---|---|---:|---:|---|---|---|---:|
| a549_ko | Confirmatory | 528 | -0.0160 | +0.2302 [+0.1955, +0.2677] | +0.2462 [+0.2132, +0.2812] | +0.2341 [+0.2021, +0.2677] | +0.1961 |
| a549_ko_cas9pos | Secondary | 526 | -0.0170 | +0.2291 [+0.1961, +0.2643] | +0.2461 [+0.2149, +0.2782] | +0.2416 [+0.2103, +0.2744] | +0.1975 |
| hek293_xucao | Secondary | 16 | -0.0189 | +0.0115 [-0.0960, +0.1468] | +0.0305 [-0.0798, +0.1882] | +0.2136 [-0.0014, +0.3910] | +0.1312 |

## Reporting-Only Legacy Contrasts (Unpaired `H-D2` & Paired `FCT − E1_ConstMed`)

| Target Line | Role | Unpaired `FCT_lo − max(Mis_hi)` | `FCT − E1 (all_nt)` [97.5% CI] | `FCT − E1 (resp_top50)` [97.5% CI] |
|---|---|---:|---|---|
| a549_ko | Confirmatory | +0.0149 (PASS) | +0.1723 [+0.1591, +0.1858] | -0.0151 [-0.0416, +0.0128] |
| a549_ko_cas9pos | Secondary | +0.0182 (PASS) | +0.1787 [+0.1647, +0.1924] | -0.0062 [-0.0307, +0.0210] |
| hek293_xucao | Secondary | -0.0478 (FAIL) | +0.1734 [+0.0826, +0.2641] | +0.2177 [+0.1029, +0.3030] |

## Reporting-Only Full Line Metrics (`all_nontarget` & Naive `resp_top50`, 97.5% CIs)

| Target Line | Role | Method | N | `all_nontarget r_orth_cal` | 97.5% CI | Naive `resp_top50 r_orth_cal` | `rho_AB_orth` | `valid_frac` |
|---|---|---|---:|---:|---|---:|---:|---:|
| a549_ko | Confirmatory | `Source_Grand_Mean` | 528 | +0.4209 | [+0.3989, +0.4427] | +0.4778 | +0.7776 | 1.0000 |
| a549_ko | Confirmatory | `Transfer_Mean` | 528 | +nan | [+nan, +nan] | +nan | +0.7776 | 0.0000 |
| a549_ko | Confirmatory | `FCT` | 528 | +0.1563 | [+0.1414, +0.1721] | +0.2392 | +0.7776 | 1.0000 |
| a549_ko | Confirmatory | `NegCtrl_FCT_MismatchedMean(cd4t_rest)` | 528 | +0.0112 | [+0.0026, +0.0201] | +0.2044 | +0.7776 | 1.0000 |
| a549_ko | Confirmatory | `NegCtrl_FCT_MismatchedMean(dld1)` | 528 | +0.1140 | [+0.1013, +0.1265] | +0.0930 | +0.7776 | 1.0000 |
| a549_ko | Confirmatory | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 528 | +0.0119 | [+0.0054, +0.0183] | +0.0516 | +0.7776 | 1.0000 |
| a549_ko | Confirmatory | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 528 | +0.0965 | [+0.0840, +0.1088] | +0.0827 | +0.7776 | 1.0000 |
| a549_ko | Confirmatory | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 528 | +0.0748 | [+0.0657, +0.0840] | +0.0994 | +0.7776 | 1.0000 |
| a549_ko | Confirmatory | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 528 | +0.0567 | [+0.0471, +0.0662] | +0.0034 | +0.7776 | 1.0000 |
| a549_ko | Confirmatory | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 528 | +0.0819 | [+0.0715, +0.0926] | +0.0943 | +0.7776 | 1.0000 |
| a549_ko | Confirmatory | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 528 | +0.0388 | [+0.0303, +0.0476] | +0.2056 | +0.7776 | 1.0000 |
| a549_ko | Confirmatory | `NegCtrl_FCT_Perm0` | 528 | -0.0143 | [-0.0197, -0.0088] | +0.1620 | +0.7776 | 1.0000 |
| a549_ko | Confirmatory | `FCT_FullNTMean` | 528 | +0.1563 | [+0.1414, +0.1721] | +0.2392 | +0.7776 | 1.0000 |
| a549_ko | Confirmatory | `Transfer_Mean_GlobalScale` | 528 | +nan | [+nan, +nan] | +nan | +0.7776 | 0.0000 |
| a549_ko | Confirmatory | `Transfer_Mean_VarScaled` | 528 | +0.2031 | [+0.1884, +0.2174] | +0.1887 | +0.7776 | 1.0000 |
| a549_ko | Confirmatory | `CIPHER_Linear` | 528 | +0.2710 | [+0.2333, +0.3094] | +0.3023 | +0.7776 | 1.0000 |
| a549_ko | Confirmatory | `Hybrid_ProjShared` | 528 | +0.1569 | [+0.1338, +0.1808] | +0.2543 | +0.7776 | 1.0000 |
| a549_ko | Confirmatory | `Hybrid_on_TMVS` | 528 | +0.2423 | [+0.2219, +0.2626] | +0.2860 | +0.7776 | 1.0000 |
| a549_ko | Confirmatory | `E1_ConstMed` | 528 | -0.0160 | [-0.0240, -0.0077] | +0.2543 | +0.7776 | 1.0000 |
| a549_ko | Confirmatory | `E2_ConstP10` | 528 | -0.0181 | [-0.0260, -0.0101] | +0.2148 | +0.7776 | 1.0000 |
| a549_ko | Confirmatory | `E2_ConstP90` | 528 | -0.0107 | [-0.0194, -0.0021] | +0.2860 | +0.7776 | 1.0000 |
| a549_ko | Confirmatory | `E3_SrcMeanCtrl` | 528 | +0.0644 | [+0.0541, +0.0750] | +0.1500 | +0.7776 | 1.0000 |
| a549_ko | Confirmatory | `E5_RankMatchedMu` | 528 | +0.0868 | [+0.0764, +0.0974] | +0.1030 | +0.7776 | 1.0000 |
| a549_ko | Confirmatory | `E6_LogFCBar` | 528 | -0.0034 | [-0.0130, +0.0060] | +0.2740 | +0.7776 | 1.0000 |
| a549_ko | Confirmatory | `E7_FCBarMinus1` | 528 | -0.0187 | [-0.0264, -0.0109] | +0.1785 | +0.7776 | 1.0000 |
| hek293_xucao | Secondary | `Source_Grand_Mean` | 16 | -0.1479 | [-0.2340, -0.0514] | -0.1189 | +0.2990 | 1.0000 |
| hek293_xucao | Secondary | `Transfer_Mean` | 16 | +nan | [+nan, +nan] | +nan | +0.2990 | 0.0000 |
| hek293_xucao | Secondary | `FCT` | 16 | +0.1544 | [+0.0600, +0.2509] | +0.1841 | +0.2990 | 1.0000 |
| hek293_xucao | Secondary | `NegCtrl_FCT_MismatchedMean(cd4t_rest)` | 16 | -0.0239 | [-0.0492, +0.0008] | -0.0797 | +0.2990 | 1.0000 |
| hek293_xucao | Secondary | `NegCtrl_FCT_MismatchedMean(dld1)` | 16 | +0.0114 | [-0.0296, +0.0558] | -0.0111 | +0.2990 | 1.0000 |
| hek293_xucao | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 16 | -0.0239 | [-0.0533, +0.0157] | -0.1075 | +0.2990 | 1.0000 |
| hek293_xucao | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 16 | -0.0250 | [-0.0633, +0.0185] | -0.0370 | +0.2990 | 1.0000 |
| hek293_xucao | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 16 | +0.0484 | [-0.0090, +0.1079] | -0.0775 | +0.2990 | 1.0000 |
| hek293_xucao | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 16 | +0.0043 | [-0.0431, +0.0540] | -0.0521 | +0.2990 | 1.0000 |
| hek293_xucao | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 16 | -0.0353 | [-0.0624, +0.0081] | +0.0084 | +0.2990 | 1.0000 |
| hek293_xucao | Secondary | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 16 | -0.0222 | [-0.0682, +0.0173] | -0.1006 | +0.2990 | 1.0000 |
| hek293_xucao | Secondary | `NegCtrl_FCT_Perm0` | 16 | -0.0093 | [-0.0280, +0.0071] | +0.0098 | +0.2990 | 1.0000 |
| hek293_xucao | Secondary | `FCT_FullNTMean` | 16 | +0.1544 | [+0.0600, +0.2509] | +0.1841 | +0.2990 | 1.0000 |
| hek293_xucao | Secondary | `Transfer_Mean_GlobalScale` | 16 | +nan | [+nan, +nan] | +nan | +0.2990 | 0.0000 |
| hek293_xucao | Secondary | `Transfer_Mean_VarScaled` | 16 | +0.2042 | [+0.1149, +0.2958] | +0.0861 | +0.2990 | 1.0000 |
| hek293_xucao | Secondary | `CIPHER_Linear` | 16 | -0.1309 | [-0.2178, -0.0198] | -0.1648 | +0.2990 | 1.0000 |
| hek293_xucao | Secondary | `Hybrid_ProjShared` | 16 | -0.0341 | [-0.0831, +0.0100] | -0.0672 | +0.2990 | 1.0000 |
| hek293_xucao | Secondary | `Hybrid_on_TMVS` | 16 | +0.1854 | [+0.0949, +0.2805] | +0.0725 | +0.2990 | 1.0000 |
| hek293_xucao | Secondary | `E1_ConstMed` | 16 | -0.0189 | [-0.0538, +0.0094] | -0.0336 | +0.2990 | 1.0000 |
| hek293_xucao | Secondary | `E2_ConstP10` | 16 | -0.0194 | [-0.0520, +0.0088] | -0.0684 | +0.2990 | 1.0000 |
| hek293_xucao | Secondary | `E2_ConstP90` | 16 | -0.0194 | [-0.0559, +0.0084] | +0.0124 | +0.2990 | 1.0000 |
| hek293_xucao | Secondary | `E3_SrcMeanCtrl` | 16 | -0.0155 | [-0.0522, +0.0278] | -0.1033 | +0.2990 | 1.0000 |
| hek293_xucao | Secondary | `E5_RankMatchedMu` | 16 | -0.0099 | [-0.0465, +0.0337] | -0.1146 | +0.2990 | 1.0000 |
| hek293_xucao | Secondary | `E6_LogFCBar` | 16 | -0.0196 | [-0.0543, +0.0080] | +0.0560 | +0.2990 | 1.0000 |
| hek293_xucao | Secondary | `E7_FCBarMinus1` | 16 | -0.0208 | [-0.0520, +0.0073] | -0.0877 | +0.2990 | 1.0000 |
| a549_ko_cas9pos | Secondary | `Source_Grand_Mean` | 526 | +0.4246 | [+0.4029, +0.4463] | +0.4853 | +0.7267 | 1.0000 |
| a549_ko_cas9pos | Secondary | `Transfer_Mean` | 526 | +nan | [+nan, +nan] | +nan | +0.7267 | 0.0000 |
| a549_ko_cas9pos | Secondary | `FCT` | 526 | +0.1616 | [+0.1464, +0.1773] | +0.2445 | +0.7267 | 1.0000 |
| a549_ko_cas9pos | Secondary | `NegCtrl_FCT_MismatchedMean(cd4t_rest)` | 526 | +0.0113 | [+0.0027, +0.0199] | +0.2058 | +0.7267 | 1.0000 |
| a549_ko_cas9pos | Secondary | `NegCtrl_FCT_MismatchedMean(dld1)` | 526 | +0.1160 | [+0.1044, +0.1282] | +0.0923 | +0.7267 | 1.0000 |
| a549_ko_cas9pos | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 526 | +0.0128 | [+0.0061, +0.0196] | +0.0475 | +0.7267 | 1.0000 |
| a549_ko_cas9pos | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 526 | +0.0971 | [+0.0855, +0.1090] | +0.0823 | +0.7267 | 1.0000 |
| a549_ko_cas9pos | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 526 | +0.0760 | [+0.0673, +0.0850] | +0.0997 | +0.7267 | 1.0000 |
| a549_ko_cas9pos | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 526 | +0.0572 | [+0.0483, +0.0663] | +0.0022 | +0.7267 | 1.0000 |
| a549_ko_cas9pos | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 526 | +0.0840 | [+0.0740, +0.0949] | +0.0962 | +0.7267 | 1.0000 |
| a549_ko_cas9pos | Secondary | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 526 | +0.0382 | [+0.0301, +0.0468] | +0.2072 | +0.7267 | 1.0000 |
| a549_ko_cas9pos | Secondary | `NegCtrl_FCT_Perm0` | 526 | -0.0148 | [-0.0202, -0.0095] | +0.1642 | +0.7267 | 1.0000 |
| a549_ko_cas9pos | Secondary | `FCT_FullNTMean` | 526 | +0.1616 | [+0.1464, +0.1773] | +0.2445 | +0.7267 | 1.0000 |
| a549_ko_cas9pos | Secondary | `Transfer_Mean_GlobalScale` | 526 | +nan | [+nan, +nan] | +nan | +0.7267 | 0.0000 |
| a549_ko_cas9pos | Secondary | `Transfer_Mean_VarScaled` | 526 | +0.2109 | [+0.1962, +0.2254] | +0.1880 | +0.7267 | 1.0000 |
| a549_ko_cas9pos | Secondary | `CIPHER_Linear` | 526 | +0.2882 | [+0.2507, +0.3267] | +0.3183 | +0.7267 | 1.0000 |
| a549_ko_cas9pos | Secondary | `Hybrid_ProjShared` | 526 | +0.1651 | [+0.1426, +0.1883] | +0.2664 | +0.7267 | 1.0000 |
| a549_ko_cas9pos | Secondary | `Hybrid_on_TMVS` | 526 | +0.2523 | [+0.2318, +0.2723] | +0.2914 | +0.7267 | 1.0000 |
| a549_ko_cas9pos | Secondary | `E1_ConstMed` | 526 | -0.0170 | [-0.0251, -0.0089] | +0.2507 | +0.7267 | 1.0000 |
| a549_ko_cas9pos | Secondary | `E2_ConstP10` | 526 | -0.0191 | [-0.0267, -0.0111] | +0.2106 | +0.7267 | 1.0000 |
| a549_ko_cas9pos | Secondary | `E2_ConstP90` | 526 | -0.0120 | [-0.0205, -0.0034] | +0.2826 | +0.7267 | 1.0000 |
| a549_ko_cas9pos | Secondary | `E3_SrcMeanCtrl` | 526 | +0.0651 | [+0.0553, +0.0753] | +0.1474 | +0.7267 | 1.0000 |
| a549_ko_cas9pos | Secondary | `E5_RankMatchedMu` | 526 | +0.0880 | [+0.0783, +0.0987] | +0.1016 | +0.7267 | 1.0000 |
| a549_ko_cas9pos | Secondary | `E6_LogFCBar` | 526 | -0.0048 | [-0.0141, +0.0045] | +0.2700 | +0.7267 | 1.0000 |
| a549_ko_cas9pos | Secondary | `E7_FCBarMinus1` | 526 | -0.0195 | [-0.0270, -0.0118] | +0.1741 | +0.7267 | 1.0000 |

## Reporting-Only `N_5` Null-Cone Analysis — `all_nontarget` (`r_cone_sp`, 97.5% CIs)

*Note: `Transfer_Mean`, `Transfer_Mean_GlobalScale`, and `Source_Grand_Mean` are **zero by exact mathematical construction** under `N_5`, and `E1_ConstMed`, `E2_ConstP10`, `E2_ConstP90`, `E6_LogFCBar`, `E7_FCBarMinus1` are **suppressed by cone construction** (`R2_N5 >= 0.99`); only `FCT`, `Mis_s`, single-draw `Perm0..99`, `E3`/`E5`, and non-cone target models are empirically informative.*

| Target Line | Method | `N_5 r_cone_sp` | 97.5% CI | `eta` | `R2_N5` | Structural Role |
|---|---|---:|---|---:|---:|---|
| a549_ko | `Source_Grand_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko | `Transfer_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko | `FCT` | +0.1556 | [+0.1408, +0.1710] | 0.8541 | 0.270532 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(cd4t_rest)` | +0.0090 | [+0.0063, +0.0117] | 0.7628 | 0.418186 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(dld1)` | +0.0870 | [+0.0781, +0.0964] | 0.8512 | 0.275387 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | -0.0054 | [-0.0097, -0.0009] | 0.8276 | 0.315007 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | +0.0427 | [+0.0375, +0.0482] | 0.8154 | 0.335156 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | +0.0535 | [+0.0479, +0.0596] | 0.8383 | 0.297240 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | +0.0143 | [+0.0096, +0.0188] | 0.8394 | 0.295438 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | +0.0478 | [+0.0411, +0.0548] | 0.8239 | 0.321243 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | +0.0200 | [+0.0179, +0.0224] | 0.7848 | 0.384067 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_Perm0` | -0.0066 | [-0.0078, -0.0054] | 0.6914 | 0.522032 | `INFORMATIVE` |
| a549_ko | `FCT_FullNTMean` | +0.1556 | [+0.1408, +0.1710] | 0.8541 | 0.270532 | `INFORMATIVE` |
| a549_ko | `Transfer_Mean_GlobalScale` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko | `Transfer_Mean_VarScaled` | +0.0916 | [+0.0853, +0.0980] | 0.4605 | 0.787906 | `INFORMATIVE` |
| a549_ko | `CIPHER_Linear` | +0.0404 | [+0.0339, +0.0476] | 0.1910 | 0.963536 | `INFORMATIVE` |
| a549_ko | `Hybrid_ProjShared` | +0.1762 | [+0.1501, +0.2048] | 0.9931 | 0.013747 | `INFORMATIVE` |
| a549_ko | `Hybrid_on_TMVS` | +0.1597 | [+0.1451, +0.1748] | 0.6284 | 0.605066 | `INFORMATIVE` |
| a549_ko | `E1_ConstMed` | +0.0000 | [+0.0000, +0.0001] | 0.0037 | 0.999986 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko | `E2_ConstP10` | +0.0000 | [+0.0000, +0.0000] | 0.0010 | 0.999999 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko | `E2_ConstP90` | -0.0000 | [-0.0000, -0.0000] | 0.0027 | 0.999993 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko | `E3_SrcMeanCtrl` | +0.0349 | [+0.0305, +0.0393] | 0.7053 | 0.502556 | `INFORMATIVE` |
| a549_ko | `E5_RankMatchedMu` | +0.0448 | [+0.0391, +0.0506] | 0.8796 | 0.226389 | `INFORMATIVE` |
| a549_ko | `E6_LogFCBar` | +0.0000 | [-0.0000, +0.0000] | 0.0098 | 0.999904 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko | `E7_FCBarMinus1` | -0.0001 | [-0.0001, -0.0000] | 0.0088 | 0.999923 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko_cas9pos | `Source_Grand_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko_cas9pos | `Transfer_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko_cas9pos | `FCT` | +0.1614 | [+0.1466, +0.1764] | 0.8546 | 0.269586 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_MismatchedMean(cd4t_rest)` | +0.0096 | [+0.0069, +0.0122] | 0.7626 | 0.418373 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_MismatchedMean(dld1)` | +0.0889 | [+0.0801, +0.0979] | 0.8525 | 0.273235 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | -0.0051 | [-0.0096, -0.0005] | 0.8292 | 0.312478 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | +0.0429 | [+0.0379, +0.0480] | 0.8160 | 0.334108 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | +0.0548 | [+0.0493, +0.0605] | 0.8396 | 0.295016 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | +0.0142 | [+0.0100, +0.0187] | 0.8407 | 0.293289 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | +0.0494 | [+0.0430, +0.0563] | 0.8250 | 0.319360 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | +0.0200 | [+0.0179, +0.0224] | 0.7846 | 0.384341 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_Perm0` | -0.0064 | [-0.0076, -0.0052] | 0.6910 | 0.522506 | `INFORMATIVE` |
| a549_ko_cas9pos | `FCT_FullNTMean` | +0.1614 | [+0.1466, +0.1764] | 0.8546 | 0.269586 | `INFORMATIVE` |
| a549_ko_cas9pos | `Transfer_Mean_GlobalScale` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko_cas9pos | `Transfer_Mean_VarScaled` | +0.0975 | [+0.0908, +0.1041] | 0.4730 | 0.776301 | `INFORMATIVE` |
| a549_ko_cas9pos | `CIPHER_Linear` | +0.0430 | [+0.0360, +0.0502] | 0.1906 | 0.963688 | `INFORMATIVE` |
| a549_ko_cas9pos | `Hybrid_ProjShared` | +0.1872 | [+0.1607, +0.2151] | 0.9932 | 0.013585 | `INFORMATIVE` |
| a549_ko_cas9pos | `Hybrid_on_TMVS` | +0.1704 | [+0.1555, +0.1857] | 0.6411 | 0.589041 | `INFORMATIVE` |
| a549_ko_cas9pos | `E1_ConstMed` | +0.0001 | [+0.0000, +0.0001] | 0.0037 | 0.999986 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko_cas9pos | `E2_ConstP10` | +0.0000 | [+0.0000, +0.0000] | 0.0010 | 0.999999 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko_cas9pos | `E2_ConstP90` | -0.0000 | [-0.0000, -0.0000] | 0.0027 | 0.999993 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko_cas9pos | `E3_SrcMeanCtrl` | +0.0354 | [+0.0315, +0.0398] | 0.7052 | 0.502704 | `INFORMATIVE` |
| a549_ko_cas9pos | `E5_RankMatchedMu` | +0.0456 | [+0.0404, +0.0514] | 0.8795 | 0.226476 | `INFORMATIVE` |
| a549_ko_cas9pos | `E6_LogFCBar` | +0.0000 | [-0.0000, +0.0000] | 0.0098 | 0.999903 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko_cas9pos | `E7_FCBarMinus1` | -0.0001 | [-0.0001, -0.0000] | 0.0088 | 0.999923 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| hek293_xucao | `Source_Grand_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| hek293_xucao | `Transfer_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| hek293_xucao | `FCT` | +0.1962 | [+0.1039, +0.2827] | 0.8501 | 0.277287 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_MismatchedMean(cd4t_rest)` | -0.0069 | [-0.0248, +0.0086] | 0.7754 | 0.398702 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_MismatchedMean(dld1)` | +0.0499 | [+0.0150, +0.0867] | 0.9449 | 0.107106 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | -0.0071 | [-0.0293, +0.0185] | 0.9572 | 0.083797 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | +0.0209 | [-0.0090, +0.0535] | 0.8852 | 0.216438 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | +0.0860 | [+0.0325, +0.1378] | 0.9126 | 0.167098 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | +0.0388 | [+0.0002, +0.0761] | 0.9447 | 0.107523 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | +0.0020 | [-0.0186, +0.0259] | 0.9432 | 0.110466 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | +0.0067 | [-0.0130, +0.0266] | 0.7887 | 0.377971 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_Perm0` | +0.0035 | [-0.0041, +0.0118] | 0.6960 | 0.515585 | `INFORMATIVE` |
| hek293_xucao | `FCT_FullNTMean` | +0.1962 | [+0.1039, +0.2827] | 0.8501 | 0.277287 | `INFORMATIVE` |
| hek293_xucao | `Transfer_Mean_GlobalScale` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| hek293_xucao | `Transfer_Mean_VarScaled` | +0.1758 | [+0.0960, +0.2542] | 0.7313 | 0.465217 | `INFORMATIVE` |
| hek293_xucao | `CIPHER_Linear` | -0.0179 | [-0.0290, -0.0053] | 0.1274 | 0.983766 | `INFORMATIVE` |
| hek293_xucao | `Hybrid_ProjShared` | -0.0402 | [-0.0858, +0.0036] | 0.9924 | 0.015185 | `INFORMATIVE` |
| hek293_xucao | `Hybrid_on_TMVS` | +0.1620 | [+0.0820, +0.2421] | 0.7459 | 0.443587 | `INFORMATIVE` |
| hek293_xucao | `E1_ConstMed` | +0.0001 | [-0.0000, +0.0002] | 0.0028 | 0.999992 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| hek293_xucao | `E2_ConstP10` | +0.0000 | [-0.0000, +0.0001] | 0.0008 | 0.999999 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| hek293_xucao | `E2_ConstP90` | -0.0000 | [-0.0001, +0.0000] | 0.0019 | 0.999996 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| hek293_xucao | `E3_SrcMeanCtrl` | +0.0247 | [+0.0019, +0.0499] | 0.7440 | 0.446419 | `INFORMATIVE` |
| hek293_xucao | `E5_RankMatchedMu` | +0.0297 | [+0.0024, +0.0592] | 0.9317 | 0.131973 | `INFORMATIVE` |
| hek293_xucao | `E6_LogFCBar` | +0.0000 | [-0.0001, +0.0002] | 0.0026 | 0.999993 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| hek293_xucao | `E7_FCBarMinus1` | -0.0001 | [-0.0003, +0.0000] | 0.0043 | 0.999982 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |

## Reporting-Only `N_5` Null-Cone Analysis — `resp_top50_split` (`r_cone_sp`, 97.5% CIs)

*Note: `Transfer_Mean`, `Transfer_Mean_GlobalScale`, and `Source_Grand_Mean` are **zero by exact mathematical construction** under `N_5`, and `E1_ConstMed`, `E2_ConstP10`, `E2_ConstP90`, `E6_LogFCBar`, `E7_FCBarMinus1` are **suppressed by cone construction** (`R2_N5 >= 0.99`); only `FCT`, `Mis_s`, single-draw `Perm0..99`, `E3`/`E5`, and non-cone target models are empirically informative.*

| Target Line | Method | `N_5 r_cone_sp` | 97.5% CI | `eta` | `R2_N5` | Structural Role |
|---|---|---:|---|---:|---:|---|
| a549_ko | `Source_Grand_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko | `Transfer_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko | `FCT` | +0.0987 | [+0.0813, +0.1175] | 0.5464 | 0.701274 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(cd4t_rest)` | +0.0220 | [+0.0107, +0.0339] | 0.7266 | 0.471989 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(dld1)` | -0.0178 | [-0.0279, -0.0072] | 0.6406 | 0.589421 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | -0.0319 | [-0.0441, -0.0204] | 0.6253 | 0.608947 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | +0.0114 | [-0.0021, +0.0251] | 0.6721 | 0.548307 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | +0.0412 | [+0.0292, +0.0532] | 0.7164 | 0.486709 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | +0.0012 | [-0.0128, +0.0146] | 0.7091 | 0.496706 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | -0.0025 | [-0.0159, +0.0107] | 0.5724 | 0.672233 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | +0.0396 | [+0.0303, +0.0488] | 0.7091 | 0.497217 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_Perm0` | -0.0221 | [-0.0293, -0.0148] | 0.4829 | 0.766758 | `INFORMATIVE` |
| a549_ko | `FCT_FullNTMean` | +0.0987 | [+0.0813, +0.1175] | 0.5464 | 0.701274 | `INFORMATIVE` |
| a549_ko | `Transfer_Mean_GlobalScale` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko | `Transfer_Mean_VarScaled` | +0.1229 | [+0.1046, +0.1430] | 0.6060 | 0.632717 | `INFORMATIVE` |
| a549_ko | `CIPHER_Linear` | +0.0371 | [+0.0300, +0.0458] | 0.1805 | 0.967415 | `INFORMATIVE` |
| a549_ko | `Hybrid_ProjShared` | +0.1671 | [+0.1344, +0.2011] | 0.8616 | 0.257680 | `INFORMATIVE` |
| a549_ko | `Hybrid_on_TMVS` | +0.1845 | [+0.1592, +0.2116] | 0.7042 | 0.504069 | `INFORMATIVE` |
| a549_ko | `E1_ConstMed` | -0.0001 | [-0.0002, -0.0001] | 0.0015 | 0.999998 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko | `E2_ConstP10` | -0.0000 | [-0.0001, -0.0000] | 0.0004 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko | `E2_ConstP90` | +0.0002 | [+0.0001, +0.0003] | 0.0022 | 0.999995 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko | `E3_SrcMeanCtrl` | -0.0055 | [-0.0166, +0.0055] | 0.5106 | 0.739210 | `INFORMATIVE` |
| a549_ko | `E5_RankMatchedMu` | +0.0016 | [-0.0113, +0.0146] | 0.6045 | 0.634506 | `INFORMATIVE` |
| a549_ko | `E6_LogFCBar` | -0.0001 | [-0.0002, +0.0000] | 0.0025 | 0.999994 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko | `E7_FCBarMinus1` | +0.0004 | [+0.0001, +0.0006] | 0.0070 | 0.999951 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko_cas9pos | `Source_Grand_Mean` | +0.0000 | [-0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko_cas9pos | `Transfer_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko_cas9pos | `FCT` | +0.1047 | [+0.0870, +0.1227] | 0.5469 | 0.700545 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_MismatchedMean(cd4t_rest)` | +0.0227 | [+0.0108, +0.0349] | 0.7258 | 0.473101 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_MismatchedMean(dld1)` | -0.0207 | [-0.0325, -0.0090] | 0.6513 | 0.575180 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | -0.0364 | [-0.0483, -0.0238] | 0.6316 | 0.600917 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | +0.0147 | [-0.0002, +0.0287] | 0.6866 | 0.528434 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | +0.0423 | [+0.0297, +0.0549] | 0.7238 | 0.475975 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | -0.0030 | [-0.0167, +0.0107] | 0.7192 | 0.482000 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | +0.0006 | [-0.0127, +0.0142] | 0.5766 | 0.667459 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | +0.0424 | [+0.0319, +0.0530] | 0.7064 | 0.500942 | `INFORMATIVE` |
| a549_ko_cas9pos | `NegCtrl_FCT_Perm0` | -0.0208 | [-0.0282, -0.0133] | 0.4792 | 0.770328 | `INFORMATIVE` |
| a549_ko_cas9pos | `FCT_FullNTMean` | +0.1047 | [+0.0870, +0.1227] | 0.5469 | 0.700545 | `INFORMATIVE` |
| a549_ko_cas9pos | `Transfer_Mean_GlobalScale` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko_cas9pos | `Transfer_Mean_VarScaled` | +0.1333 | [+0.1141, +0.1548] | 0.6128 | 0.624413 | `INFORMATIVE` |
| a549_ko_cas9pos | `CIPHER_Linear` | +0.0388 | [+0.0308, +0.0470] | 0.1783 | 0.968223 | `INFORMATIVE` |
| a549_ko_cas9pos | `Hybrid_ProjShared` | +0.1787 | [+0.1465, +0.2120] | 0.8562 | 0.266854 | `INFORMATIVE` |
| a549_ko_cas9pos | `Hybrid_on_TMVS` | +0.1984 | [+0.1712, +0.2265] | 0.7076 | 0.499285 | `INFORMATIVE` |
| a549_ko_cas9pos | `E1_ConstMed` | -0.0001 | [-0.0002, -0.0001] | 0.0015 | 0.999998 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko_cas9pos | `E2_ConstP10` | -0.0000 | [-0.0001, -0.0000] | 0.0004 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko_cas9pos | `E2_ConstP90` | +0.0002 | [+0.0001, +0.0003] | 0.0023 | 0.999995 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko_cas9pos | `E3_SrcMeanCtrl` | -0.0060 | [-0.0176, +0.0053] | 0.5194 | 0.730100 | `INFORMATIVE` |
| a549_ko_cas9pos | `E5_RankMatchedMu` | +0.0003 | [-0.0136, +0.0135] | 0.6130 | 0.624014 | `INFORMATIVE` |
| a549_ko_cas9pos | `E6_LogFCBar` | -0.0001 | [-0.0002, +0.0000] | 0.0026 | 0.999993 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko_cas9pos | `E7_FCBarMinus1` | +0.0004 | [+0.0002, +0.0007] | 0.0070 | 0.999952 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| hek293_xucao | `Source_Grand_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| hek293_xucao | `Transfer_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| hek293_xucao | `FCT` | +0.1158 | [+0.0072, +0.2253] | 0.5370 | 0.711584 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_MismatchedMean(cd4t_rest)` | -0.0194 | [-0.1691, +0.0675] | 0.6003 | 0.639595 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_MismatchedMean(dld1)` | -0.0369 | [-0.2236, +0.1176] | 0.7581 | 0.425332 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | -0.0409 | [-0.0949, +0.0449] | 0.6626 | 0.560564 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | -0.0938 | [-0.1512, -0.0287] | 0.7094 | 0.496761 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | -0.0567 | [-0.1402, +0.0256] | 0.7202 | 0.480585 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | -0.0632 | [-0.1287, +0.0027] | 0.6437 | 0.585258 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | -0.0097 | [-0.0632, +0.0668] | 0.6249 | 0.609060 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | -0.0163 | [-0.1192, +0.0307] | 0.4544 | 0.793089 | `INFORMATIVE` |
| hek293_xucao | `NegCtrl_FCT_Perm0` | -0.0085 | [-0.0518, +0.0365] | 0.4389 | 0.807034 | `INFORMATIVE` |
| hek293_xucao | `FCT_FullNTMean` | +0.1158 | [+0.0072, +0.2253] | 0.5370 | 0.711584 | `INFORMATIVE` |
| hek293_xucao | `Transfer_Mean_GlobalScale` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| hek293_xucao | `Transfer_Mean_VarScaled` | -0.0446 | [-0.1608, +0.1195] | 0.3963 | 0.842960 | `INFORMATIVE` |
| hek293_xucao | `CIPHER_Linear` | -0.0065 | [-0.0231, +0.0026] | 0.0616 | 0.996199 | `INFORMATIVE` |
| hek293_xucao | `Hybrid_ProjShared` | -0.1017 | [-0.1978, +0.0062] | 0.9334 | 0.128613 | `INFORMATIVE` |
| hek293_xucao | `Hybrid_on_TMVS` | -0.0532 | [-0.1640, +0.1077] | 0.4028 | 0.837731 | `INFORMATIVE` |
| hek293_xucao | `E1_ConstMed` | +0.0000 | [-0.0001, +0.0001] | 0.0005 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| hek293_xucao | `E2_ConstP10` | -0.0000 | [-0.0000, +0.0000] | 0.0002 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| hek293_xucao | `E2_ConstP90` | -0.0001 | [-0.0004, +0.0001] | 0.0009 | 0.999999 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| hek293_xucao | `E3_SrcMeanCtrl` | -0.0807 | [-0.1676, -0.0063] | 0.5512 | 0.695483 | `INFORMATIVE` |
| hek293_xucao | `E5_RankMatchedMu` | -0.0977 | [-0.2098, +0.0001] | 0.7063 | 0.499519 | `INFORMATIVE` |
| hek293_xucao | `E6_LogFCBar` | +0.0012 | [-0.0001, +0.0038] | 0.0047 | 0.999978 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| hek293_xucao | `E7_FCBarMinus1` | +0.0001 | [-0.0001, +0.0002] | 0.0017 | 0.999997 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |

## 100-Perm Empirical Null (`all_nontarget`, `SA-8`, & Single-Draw `N_5`)

- **a549_ko**: `all_nontarget` med=-0.0093, Q97.5=+0.0023; `SA-8 split` med=+0.1318, Q97.5=+0.1961; `N_5` single-draw med=+0.0007, med|r|=+0.0042, Q97.5=+0.0112 (med_eta=0.6890)
- **a549_ko_cas9pos**: `all_nontarget` med=-0.0100, Q97.5=+0.0016; `SA-8 split` med=+0.1300, Q97.5=+0.1975; `N_5` single-draw med=+0.0006, med|r|=+0.0044, Q97.5=+0.0113 (med_eta=0.6886)
- **hek293_xucao**: `all_nontarget` med=-0.0135, Q97.5=+0.0002; `SA-8 split` med=+0.0181, Q97.5=+0.1312; `N_5` single-draw med=-0.0010, med|r|=+0.0042, Q97.5=+0.0131 (med_eta=0.6846)

## Mandatory Governance, Provenance & Prior-Art Disclosures (`PREREG_E §0.2, §5, §7.2, §8`)

- **Post-Hoc Provenance**: `H-E1` (paired `FCT - Mis_s` on `all_nontarget`) and `H-E2` (`SA-8 resp_top50_split` constant-baseline inflation) were specified in `PREREG_E` after unblinding Rounds 9-11 (`N=12` lines) and evaluated zero-shot on unseen Round 12 lines.
- **Claim 2 Non-Mechanistic Scope**: Confirmation of `H-E2` establishes an **empirical evaluation artifact strictly under `r_orth_cal`** (`resp_top50_split` vs `all_nontarget`) — namely that single-line constant baseline `E1_ConstMed` attains inflated disattenuated residual correlation under magnitude-selected subsets — and **must NOT be extrapolated** as a mechanistic biological claim or applied directly to raw DEG discovery benchmarks without independent validation.
- **Modality Boundary & Single-Line Generalizability**: Source cohorts (`r9_6`, `r9_r10_8`) comprise **CRISPRi** perturbations whereas `a549_ko` (`GSE345058`) is **CRISPRko**. Failure on `a549_ko` is confounded with CRISPRi->CRISPRko modality shift rather than purely same-modality failure; passage constitutes cross-modality generalization evidence. Under `K_conf = 1`, population generalizability across additional lineages remains explicitly limited.
- **Prior-Art & Data Attribution**: Baseline-aware perturbation transfer context acknowledges **Kawsar et al. (2026), `CanPert/CASTER` (`doi:10.64898/2026.11.15.766815`)**; `a549_ko` CRISPRko dataset derives from **Liu et al. (2026), `GSE345058`**; `hek293_xucao` CRISPRi derives from **Xu & Cao (2023), PerturbSci-Kinetics**.
- **Reporting-Only §5 Scope**: Legacy unpaired `H-D2`, paired `FCT - E1_ConstMed`, `N_5` null-cone suites (`all_nontarget` & `resp_top50_split`), `SA-Src8` (`--source-mode r9_r10_8`), `SA-Cas9` (`a549_ko_cas9pos`), naive `resp_top50`, and `H-D_util` are strictly descriptive and carry zero multiplicity weight.
