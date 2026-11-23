# Round 13 Prospective Replication Report (`PREREG_F` — `holdout_r13_r9_6_log1p`)

- **Phase**: `holdout_r13` | **Source Mode**: `r9_6` (`6` sources) | **PREREG_E Commit**: `5b51e5b079c2c01b44408b9ed9f8d6027b4ae4f5` | **Code Freeze**: `540ba1896fab3bd8f240190175befe4a8ef5fe7d`
- **Bootstrap**: `np.random.default_rng(20261223)` x `4000` (Bonferroni 2-claim + IUT two-sided `97.5%` CIs `[1.25%, 98.75%]`)
- **Cohort Disclosure**: **K_conf = 0 (count gate failed); descriptive reporting only.**
- **Secondary / Non-Confirmatory Lines (Reporting-Only)**: `['mac_ctrl', 'mac_lps']`
- **Confirmatory Claim 1 (`H-E1` Target-Line Specificity)**: **NOT_EVALUATED_NO_CONFIRMATORY_LINES**
- **Confirmatory Claim 2 (`H-E2` Persistent `SA-8` Constant-Baseline Inflation)**: **NOT_EVALUATED_NO_CONFIRMATORY_LINES**

## `H-E1b`: Paired Specificity Contrasts `FCT − NegCtrl_FCT_MismatchedMean(s)` (`all_nontarget`, 97.5% CI)

| Target Line | Role | Mismatched Source | Estimate | 97.5% CI | `lo > 0` |
|---|---|---|---:|---|---|
| mac_ctrl | Secondary | `nadig_hepg2` | +0.1634 | [+0.0297, +0.3799] | **PASS** |
| mac_ctrl | Secondary | `nadig_jurkat` | +0.1843 | [+0.1249, +0.3364] | **PASS** |
| mac_ctrl | Secondary | `nourreddine_kolf21j` | +0.1966 | [+0.0988, +0.3867] | **PASS** |
| mac_ctrl | Secondary | `replogle_k562_essential` | +0.1954 | [+0.1088, +0.3993] | **PASS** |
| mac_ctrl | Secondary | `replogle_rpe1` | +0.1743 | [+0.0441, +0.4002] | **PASS** |
| mac_ctrl | Secondary | `vcc2025_h1` | +0.1126 | [+0.0578, +0.1849] | **PASS** |
| mac_lps | Secondary | `nadig_hepg2` | +0.1949 | [+0.0932, +0.3373] | **PASS** |
| mac_lps | Secondary | `nadig_jurkat` | +0.1952 | [+0.1042, +0.3359] | **PASS** |
| mac_lps | Secondary | `nourreddine_kolf21j` | +0.1888 | [+0.0966, +0.3307] | **PASS** |
| mac_lps | Secondary | `replogle_k562_essential` | +0.1845 | [+0.0823, +0.3375] | **PASS** |
| mac_lps | Secondary | `replogle_rpe1` | +0.1992 | [+0.0869, +0.3527] | **PASS** |
| mac_lps | Secondary | `vcc2025_h1` | +0.0876 | [-0.0067, +0.1452] | **FAIL** |

## `H-E2` / `SA-8`: Split-Half Top-50 (`resp_top50_split`) & Joint Cross-Mask `E1(split) − E1(all)`

| Target Line | Role | N_common | `E1(all_common)` | `E1(split)` [97.5% CI] | Joint `E1(split) − E1(all)` [97.5% CI] | `FCT(split)` [97.5% CI] | `Perm100 split` Q97.5 |
|---|---|---:|---:|---|---|---|---:|
| mac_ctrl | Secondary | 3 | +0.0323 | +0.1102 [-0.0380, +0.3878] | +0.0780 [-0.0201, +0.2123] | +0.1310 [-0.0799, +0.4055] | +0.1734 |
| mac_lps | Secondary | 5 | +0.0440 | +0.1510 [-0.0264, +0.4773] | +0.1070 [+0.0000, +0.3543] | +0.1548 [-0.0047, +0.4997] | +0.1707 |

## Reporting-Only Legacy Contrasts (Unpaired `H-D2` & Paired `FCT − E1_ConstMed`)

| Target Line | Role | Unpaired `FCT_lo − max(Mis_hi)` | `FCT − E1 (all_nt)` [97.5% CI] | `FCT − E1 (resp_top50)` [97.5% CI] |
|---|---|---:|---|---|
| mac_ctrl | Secondary | -0.1456 (FAIL) | +0.1488 [+0.0977, +0.2349] | +0.0121 [-0.0453, +0.0555] |
| mac_lps | Secondary | -0.1597 (FAIL) | +0.1445 [+0.0798, +0.2024] | -0.0266 [-0.0431, +0.0352] |

## Reporting-Only Full Line Metrics (`all_nontarget` & Naive `resp_top50`, 97.5% CIs)

| Target Line | Role | Method | N | `all_nontarget r_orth_cal` | 97.5% CI | Naive `resp_top50 r_orth_cal` | `rho_AB_orth` | `valid_frac` |
|---|---|---|---:|---:|---|---:|---:|---:|
| mac_ctrl | Secondary | `Source_Grand_Mean` | 3 | -0.0818 | [-0.2327, +0.0015] | +0.0304 | +0.3748 | 1.0000 |
| mac_ctrl | Secondary | `Transfer_Mean` | 3 | +nan | [+nan, +nan] | +nan | +0.3748 | 0.0000 |
| mac_ctrl | Secondary | `FCT` | 3 | +0.1810 | [+0.0799, +0.4104] | +0.1089 | +0.3748 | 1.0000 |
| mac_ctrl | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 3 | +0.0176 | [+0.0008, +0.0502] | +0.0409 | +0.3748 | 1.0000 |
| mac_ctrl | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 3 | -0.0032 | [-0.0449, +0.0740] | +0.0262 | +0.3748 | 1.0000 |
| mac_ctrl | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 3 | -0.0156 | [-0.0370, +0.0237] | -0.0979 | +0.3748 | 1.0000 |
| mac_ctrl | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 3 | -0.0144 | [-0.0289, +0.0111] | -0.0077 | +0.3748 | 1.0000 |
| mac_ctrl | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 3 | +0.0067 | [-0.0081, +0.0358] | -0.2135 | +0.3748 | 1.0000 |
| mac_ctrl | Secondary | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 3 | +0.0684 | [+0.0103, +0.2255] | -0.0484 | +0.3748 | 1.0000 |
| mac_ctrl | Secondary | `NegCtrl_FCT_Perm0` | 3 | +0.0243 | [-0.0064, +0.0817] | +0.0223 | +0.3748 | 1.0000 |
| mac_ctrl | Secondary | `FCT_FullNTMean` | 3 | +0.1810 | [+0.0799, +0.4104] | +0.1089 | +0.3748 | 1.0000 |
| mac_ctrl | Secondary | `Transfer_Mean_GlobalScale` | 3 | +nan | [+nan, +nan] | +nan | +0.3748 | 0.0000 |
| mac_ctrl | Secondary | `Transfer_Mean_VarScaled` | 3 | +0.1942 | [+0.1047, +0.3637] | +0.1552 | +0.3748 | 1.0000 |
| mac_ctrl | Secondary | `CIPHER_Linear` | 3 | +0.0768 | [-0.0752, +0.2631] | +0.0739 | +0.3748 | 1.0000 |
| mac_ctrl | Secondary | `Hybrid_ProjShared` | 3 | +0.0537 | [-0.0626, +0.2011] | +0.0590 | +0.3748 | 1.0000 |
| mac_ctrl | Secondary | `Hybrid_on_TMVS` | 3 | +0.1880 | [+0.1015, +0.3901] | +0.1402 | +0.3748 | 1.0000 |
| mac_ctrl | Secondary | `E1_ConstMed` | 3 | +0.0323 | [-0.0246, +0.1755] | +0.0968 | +0.3748 | 1.0000 |
| mac_ctrl | Secondary | `E2_ConstP10` | 3 | +0.0370 | [-0.0215, +0.1819] | +0.1040 | +0.3748 | 1.0000 |
| mac_ctrl | Secondary | `E2_ConstP90` | 3 | +0.0260 | [-0.0281, +0.1638] | +0.0844 | +0.3748 | 1.0000 |
| mac_ctrl | Secondary | `E3_SrcMeanCtrl` | 3 | +0.0221 | [-0.0198, +0.1349] | -0.0827 | +0.3748 | 1.0000 |
| mac_ctrl | Secondary | `E5_RankMatchedMu` | 3 | +0.0116 | [-0.0165, +0.0749] | -0.1415 | +0.3748 | 1.0000 |
| mac_ctrl | Secondary | `E6_LogFCBar` | 3 | +0.0182 | [-0.0318, +0.1457] | +0.0652 | +0.3748 | 1.0000 |
| mac_ctrl | Secondary | `E7_FCBarMinus1` | 3 | +0.0387 | [-0.0204, +0.1834] | +0.1057 | +0.3748 | 1.0000 |
| mac_lps | Secondary | `Source_Grand_Mean` | 5 | -0.0834 | [-0.1834, -0.0323] | -0.0121 | +0.4009 | 1.0000 |
| mac_lps | Secondary | `Transfer_Mean` | 5 | +nan | [+nan, +nan] | +nan | +0.4009 | 0.0000 |
| mac_lps | Secondary | `FCT` | 5 | +0.1885 | [+0.0897, +0.3522] | +0.1418 | +0.4009 | 1.0000 |
| mac_lps | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 5 | -0.0064 | [-0.0292, +0.0261] | -0.0333 | +0.4009 | 1.0000 |
| mac_lps | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 5 | -0.0066 | [-0.0297, +0.0269] | -0.0521 | +0.4009 | 1.0000 |
| mac_lps | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 5 | -0.0003 | [-0.0260, +0.0558] | -0.0204 | +0.4009 | 1.0000 |
| mac_lps | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 5 | +0.0041 | [-0.0232, +0.0785] | -0.0356 | +0.4009 | 1.0000 |
| mac_lps | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 5 | -0.0107 | [-0.0390, +0.0125] | -0.1508 | +0.4009 | 1.0000 |
| mac_lps | Secondary | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 5 | +0.1009 | [+0.0158, +0.2494] | +0.0690 | +0.4009 | 1.0000 |
| mac_lps | Secondary | `NegCtrl_FCT_Perm0` | 5 | +0.0191 | [-0.0273, +0.0877] | +0.0798 | +0.4009 | 1.0000 |
| mac_lps | Secondary | `FCT_FullNTMean` | 5 | +0.1885 | [+0.0897, +0.3522] | +0.1418 | +0.4009 | 1.0000 |
| mac_lps | Secondary | `Transfer_Mean_GlobalScale` | 5 | +nan | [+nan, +nan] | +nan | +0.4009 | 0.0000 |
| mac_lps | Secondary | `Transfer_Mean_VarScaled` | 5 | +0.2101 | [+0.1184, +0.3398] | +0.1745 | +0.4009 | 1.0000 |
| mac_lps | Secondary | `CIPHER_Linear` | 5 | +0.0292 | [-0.1404, +0.2069] | +0.0851 | +0.4009 | 1.0000 |
| mac_lps | Secondary | `Hybrid_ProjShared` | 5 | -0.0147 | [-0.1362, +0.1344] | +0.0681 | +0.4009 | 1.0000 |
| mac_lps | Secondary | `Hybrid_on_TMVS` | 5 | +0.1628 | [+0.0624, +0.3017] | +0.1758 | +0.4009 | 1.0000 |
| mac_lps | Secondary | `E1_ConstMed` | 5 | +0.0440 | [-0.0275, +0.1525] | +0.1684 | +0.4009 | 1.0000 |
| mac_lps | Secondary | `E2_ConstP10` | 5 | +0.0476 | [-0.0250, +0.1565] | +0.1630 | +0.4009 | 1.0000 |
| mac_lps | Secondary | `E2_ConstP90` | 5 | +0.0387 | [-0.0309, +0.1429] | +0.1687 | +0.4009 | 1.0000 |
| mac_lps | Secondary | `E3_SrcMeanCtrl` | 5 | +0.0371 | [-0.0221, +0.1550] | +0.0252 | +0.4009 | 1.0000 |
| mac_lps | Secondary | `E5_RankMatchedMu` | 5 | +0.0195 | [-0.0189, +0.0992] | -0.0598 | +0.4009 | 1.0000 |
| mac_lps | Secondary | `E6_LogFCBar` | 5 | +0.0321 | [-0.0348, +0.1295] | +0.1616 | +0.4009 | 1.0000 |
| mac_lps | Secondary | `E7_FCBarMinus1` | 5 | +0.0486 | [-0.0242, +0.1562] | +0.1603 | +0.4009 | 1.0000 |

## Reporting-Only `N_5` Null-Cone Analysis — `all_nontarget` (`r_cone_sp`, 97.5% CIs)

*Note: `Transfer_Mean`, `Transfer_Mean_GlobalScale`, and `Source_Grand_Mean` are **zero by exact mathematical construction** under `N_5`, and `E1_ConstMed`, `E2_ConstP10`, `E2_ConstP90`, `E6_LogFCBar`, `E7_FCBarMinus1` are **suppressed by cone construction** (`R2_N5 >= 0.99`); only `FCT`, `Mis_s`, single-draw `Perm0..99`, `E3`/`E5`, and non-cone target models are empirically informative.*

| Target Line | Method | `N_5 r_cone_sp` | 97.5% CI | `eta` | `R2_N5` | Structural Role |
|---|---|---:|---|---:|---:|---|
| mac_ctrl | `Source_Grand_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| mac_ctrl | `Transfer_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| mac_ctrl | `FCT` | +0.1401 | [+0.0759, +0.2690] | 0.7709 | 0.405651 | `INFORMATIVE` |
| mac_ctrl | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | +0.0224 | [+0.0114, +0.0532] | 0.9373 | 0.121504 | `INFORMATIVE` |
| mac_ctrl | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | +0.0027 | [-0.0337, +0.0597] | 0.9100 | 0.171984 | `INFORMATIVE` |
| mac_ctrl | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | -0.0102 | [-0.0203, +0.0071] | 0.8900 | 0.207895 | `INFORMATIVE` |
| mac_ctrl | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | -0.0048 | [-0.0168, +0.0099] | 0.9483 | 0.100770 | `INFORMATIVE` |
| mac_ctrl | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | +0.0094 | [-0.0055, +0.0418] | 0.9746 | 0.050201 | `INFORMATIVE` |
| mac_ctrl | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | +0.0636 | [+0.0346, +0.1617] | 0.7914 | 0.373759 | `INFORMATIVE` |
| mac_ctrl | `NegCtrl_FCT_Perm0` | +0.0051 | [-0.0281, +0.0218] | 0.7139 | 0.490287 | `INFORMATIVE` |
| mac_ctrl | `FCT_FullNTMean` | +0.1401 | [+0.0759, +0.2690] | 0.7709 | 0.405651 | `INFORMATIVE` |
| mac_ctrl | `Transfer_Mean_GlobalScale` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| mac_ctrl | `Transfer_Mean_VarScaled` | +0.1552 | [+0.0926, +0.2600] | 0.8933 | 0.201932 | `INFORMATIVE` |
| mac_ctrl | `CIPHER_Linear` | +0.0096 | [-0.0180, +0.0685] | 0.1546 | 0.976103 | `INFORMATIVE` |
| mac_ctrl | `Hybrid_ProjShared` | +0.0571 | [-0.0623, +0.1990] | 0.9959 | 0.008220 | `INFORMATIVE` |
| mac_ctrl | `Hybrid_on_TMVS` | +0.1570 | [+0.0902, +0.2985] | 0.9234 | 0.147257 | `INFORMATIVE` |
| mac_ctrl | `E1_ConstMed` | +0.0000 | [-0.0000, +0.0001] | 0.0010 | 0.999999 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| mac_ctrl | `E2_ConstP10` | +0.0000 | [-0.0000, +0.0000] | 0.0001 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| mac_ctrl | `E2_ConstP90` | -0.0000 | [-0.0000, +0.0000] | 0.0008 | 0.999999 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| mac_ctrl | `E3_SrcMeanCtrl` | +0.0219 | [+0.0022, +0.0840] | 0.7974 | 0.364137 | `INFORMATIVE` |
| mac_ctrl | `E5_RankMatchedMu` | +0.0241 | [-0.0029, +0.0993] | 0.9662 | 0.066517 | `INFORMATIVE` |
| mac_ctrl | `E6_LogFCBar` | +0.0000 | [-0.0000, +0.0000] | 0.0007 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| mac_ctrl | `E7_FCBarMinus1` | -0.0000 | [-0.0001, +0.0000] | 0.0008 | 0.999999 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| mac_lps | `Source_Grand_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| mac_lps | `Transfer_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| mac_lps | `FCT` | +0.1457 | [+0.0806, +0.2403] | 0.7779 | 0.394816 | `INFORMATIVE` |
| mac_lps | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | -0.0023 | [-0.0247, +0.0227] | 0.9458 | 0.105413 | `INFORMATIVE` |
| mac_lps | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | -0.0034 | [-0.0153, +0.0051] | 0.8984 | 0.192859 | `INFORMATIVE` |
| mac_lps | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | -0.0025 | [-0.0115, +0.0224] | 0.8925 | 0.203517 | `INFORMATIVE` |
| mac_lps | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | +0.0095 | [-0.0050, +0.0520] | 0.9421 | 0.112374 | `INFORMATIVE` |
| mac_lps | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | -0.0011 | [-0.0263, +0.0237] | 0.9661 | 0.066581 | `INFORMATIVE` |
| mac_lps | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | +0.0860 | [+0.0293, +0.1765] | 0.7872 | 0.380311 | `INFORMATIVE` |
| mac_lps | `NegCtrl_FCT_Perm0` | -0.0083 | [-0.0138, +0.0026] | 0.7191 | 0.482955 | `INFORMATIVE` |
| mac_lps | `FCT_FullNTMean` | +0.1457 | [+0.0806, +0.2403] | 0.7779 | 0.394816 | `INFORMATIVE` |
| mac_lps | `Transfer_Mean_GlobalScale` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| mac_lps | `Transfer_Mean_VarScaled` | +0.1713 | [+0.1020, +0.2626] | 0.8788 | 0.227736 | `INFORMATIVE` |
| mac_lps | `CIPHER_Linear` | -0.0004 | [-0.0358, +0.0524] | 0.1873 | 0.964927 | `INFORMATIVE` |
| mac_lps | `Hybrid_ProjShared` | -0.0137 | [-0.1389, +0.1347] | 0.9941 | 0.011742 | `INFORMATIVE` |
| mac_lps | `Hybrid_on_TMVS` | +0.1337 | [+0.0515, +0.2372] | 0.9127 | 0.166913 | `INFORMATIVE` |
| mac_lps | `E1_ConstMed` | +0.0000 | [-0.0000, +0.0000] | 0.0009 | 0.999999 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| mac_lps | `E2_ConstP10` | +0.0000 | [-0.0000, +0.0000] | 0.0000 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| mac_lps | `E2_ConstP90` | -0.0000 | [-0.0000, +0.0000] | 0.0007 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| mac_lps | `E3_SrcMeanCtrl` | +0.0298 | [+0.0007, +0.0883] | 0.7930 | 0.371104 | `INFORMATIVE` |
| mac_lps | `E5_RankMatchedMu` | +0.0348 | [-0.0018, +0.1048] | 0.9687 | 0.061622 | `INFORMATIVE` |
| mac_lps | `E6_LogFCBar` | +0.0000 | [-0.0000, +0.0000] | 0.0006 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| mac_lps | `E7_FCBarMinus1` | -0.0000 | [-0.0000, +0.0000] | 0.0007 | 0.999999 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |

## Reporting-Only `N_5` Null-Cone Analysis — `resp_top50_split` (`r_cone_sp`, 97.5% CIs)

*Note: `Transfer_Mean`, `Transfer_Mean_GlobalScale`, and `Source_Grand_Mean` are **zero by exact mathematical construction** under `N_5`, and `E1_ConstMed`, `E2_ConstP10`, `E2_ConstP90`, `E6_LogFCBar`, `E7_FCBarMinus1` are **suppressed by cone construction** (`R2_N5 >= 0.99`); only `FCT`, `Mis_s`, single-draw `Perm0..99`, `E3`/`E5`, and non-cone target models are empirically informative.*

| Target Line | Method | `N_5 r_cone_sp` | 97.5% CI | `eta` | `R2_N5` | Structural Role |
|---|---|---:|---|---:|---:|---|
| mac_ctrl | `Source_Grand_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| mac_ctrl | `Transfer_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| mac_ctrl | `FCT` | +0.0165 | [-0.0583, +0.0431] | 0.2931 | 0.914018 | `INFORMATIVE` |
| mac_ctrl | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | -0.0688 | [-0.1207, +0.0851] | 0.7674 | 0.410749 | `INFORMATIVE` |
| mac_ctrl | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | +0.0397 | [+0.0113, +0.1268] | 0.7343 | 0.452756 | `INFORMATIVE` |
| mac_ctrl | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | -0.0835 | [-0.1534, +0.0572] | 0.8076 | 0.347238 | `INFORMATIVE` |
| mac_ctrl | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | -0.0384 | [-0.1093, +0.1293] | 0.8251 | 0.318463 | `INFORMATIVE` |
| mac_ctrl | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | -0.1412 | [-0.2308, +0.1080] | 0.8170 | 0.331825 | `INFORMATIVE` |
| mac_ctrl | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | -0.0247 | [-0.0705, +0.1041] | 0.6136 | 0.622809 | `INFORMATIVE` |
| mac_ctrl | `NegCtrl_FCT_Perm0` | -0.0320 | [-0.0624, -0.0287] | 0.4825 | 0.766219 | `INFORMATIVE` |
| mac_ctrl | `FCT_FullNTMean` | +0.0165 | [-0.0583, +0.0431] | 0.2931 | 0.914018 | `INFORMATIVE` |
| mac_ctrl | `Transfer_Mean_GlobalScale` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| mac_ctrl | `Transfer_Mean_VarScaled` | -0.0026 | [-0.1661, +0.0577] | 0.5256 | 0.723239 | `INFORMATIVE` |
| mac_ctrl | `CIPHER_Linear` | -0.0144 | [-0.0884, +0.0529] | 0.2113 | 0.955363 | `INFORMATIVE` |
| mac_ctrl | `Hybrid_ProjShared` | -0.0283 | [-0.1992, +0.3084] | 0.7971 | 0.361531 | `INFORMATIVE` |
| mac_ctrl | `Hybrid_on_TMVS` | -0.0170 | [-0.0717, +0.0940] | 0.5697 | 0.675377 | `INFORMATIVE` |
| mac_ctrl | `E1_ConstMed` | +0.0000 | [-0.0000, +0.0001] | 0.0004 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| mac_ctrl | `E2_ConstP10` | +0.0000 | [-0.0000, +0.0000] | 0.0000 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| mac_ctrl | `E2_ConstP90` | -0.0000 | [-0.0000, +0.0000] | 0.0003 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| mac_ctrl | `E3_SrcMeanCtrl` | -0.0686 | [-0.1293, +0.0785] | 0.6259 | 0.605390 | `INFORMATIVE` |
| mac_ctrl | `E5_RankMatchedMu` | -0.0715 | [-0.1561, +0.0903] | 0.7798 | 0.390780 | `INFORMATIVE` |
| mac_ctrl | `E6_LogFCBar` | +0.0000 | [+0.0000, +0.0000] | 0.0002 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| mac_ctrl | `E7_FCBarMinus1` | -0.0000 | [-0.0000, +0.0000] | 0.0003 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| mac_lps | `Source_Grand_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| mac_lps | `Transfer_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| mac_lps | `FCT` | +0.0363 | [-0.0093, +0.1105] | 0.3701 | 0.862762 | `INFORMATIVE` |
| mac_lps | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | -0.0713 | [-0.1734, +0.0614] | 0.8170 | 0.332495 | `INFORMATIVE` |
| mac_lps | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | -0.0641 | [-0.1924, +0.0149] | 0.7800 | 0.388857 | `INFORMATIVE` |
| mac_lps | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | -0.0483 | [-0.1574, +0.0788] | 0.8459 | 0.284400 | `INFORMATIVE` |
| mac_lps | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | -0.0766 | [-0.1970, +0.0590] | 0.8213 | 0.324951 | `INFORMATIVE` |
| mac_lps | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | -0.1140 | [-0.3429, +0.1507] | 0.8493 | 0.278556 | `INFORMATIVE` |
| mac_lps | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | -0.1006 | [-0.2156, +0.0319] | 0.6068 | 0.631609 | `INFORMATIVE` |
| mac_lps | `NegCtrl_FCT_Perm0` | +0.0034 | [-0.0655, +0.0513] | 0.6413 | 0.586746 | `INFORMATIVE` |
| mac_lps | `FCT_FullNTMean` | +0.0363 | [-0.0093, +0.1105] | 0.3701 | 0.862762 | `INFORMATIVE` |
| mac_lps | `Transfer_Mean_GlobalScale` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| mac_lps | `Transfer_Mean_VarScaled` | +0.1273 | [+0.0721, +0.1745] | 0.5910 | 0.647512 | `INFORMATIVE` |
| mac_lps | `CIPHER_Linear` | -0.0131 | [-0.1024, +0.1110] | 0.2323 | 0.946011 | `INFORMATIVE` |
| mac_lps | `Hybrid_ProjShared` | -0.0473 | [-0.3524, +0.3676] | 0.8932 | 0.202094 | `INFORMATIVE` |
| mac_lps | `Hybrid_on_TMVS` | +0.0729 | [-0.1336, +0.2880] | 0.6854 | 0.526798 | `INFORMATIVE` |
| mac_lps | `E1_ConstMed` | +0.0000 | [-0.0000, +0.0001] | 0.0005 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| mac_lps | `E2_ConstP10` | +0.0000 | [-0.0000, +0.0000] | 0.0000 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| mac_lps | `E2_ConstP90` | -0.0000 | [-0.0001, +0.0000] | 0.0004 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| mac_lps | `E3_SrcMeanCtrl` | -0.1155 | [-0.2506, +0.0507] | 0.6661 | 0.555889 | `INFORMATIVE` |
| mac_lps | `E5_RankMatchedMu` | -0.1498 | [-0.3246, +0.0470] | 0.8057 | 0.350198 | `INFORMATIVE` |
| mac_lps | `E6_LogFCBar` | +0.0000 | [-0.0000, +0.0000] | 0.0003 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| mac_lps | `E7_FCBarMinus1` | -0.0000 | [-0.0000, +0.0000] | 0.0005 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |

## 100-Perm Empirical Null (`all_nontarget`, `SA-8`, & Single-Draw `N_5`)

- **mac_ctrl**: `all_nontarget` med=+0.0191, Q97.5=+0.0361; `SA-8 split` med=+0.0776, Q97.5=+0.1734; `N_5` single-draw med=-0.0003, med|r|=+0.0056, Q97.5=+0.0161 (med_eta=0.7149)
- **mac_lps**: `all_nontarget` med=+0.0270, Q97.5=+0.0384; `SA-8 split` med=+0.0875, Q97.5=+0.1707; `N_5` single-draw med=+0.0004, med|r|=+0.0039, Q97.5=+0.0103 (med_eta=0.7211)

## Mandatory Governance, Provenance & Prior-Art Disclosures (`PREREG_E §0.2, §5, §7.2, §8`)

- **Post-Hoc Provenance**: `H-E1` (paired `FCT - Mis_s` on `all_nontarget`) and `H-E2` (`SA-8 resp_top50_split` constant-baseline inflation) were specified in `PREREG_E` after unblinding Rounds 9-11 (`N=12` lines) and evaluated zero-shot on unseen Round 12 lines.
- **Claim 2 Non-Mechanistic Scope**: Confirmation of `H-E2` establishes an **empirical evaluation artifact strictly under `r_orth_cal`** (`resp_top50_split` vs `all_nontarget`) — namely that single-line constant baseline `E1_ConstMed` attains inflated disattenuated residual correlation under magnitude-selected subsets — and **must NOT be extrapolated** as a mechanistic biological claim or applied directly to raw DEG discovery benchmarks without independent validation.
- **Modality Boundary & Single-Line Generalizability**: Source cohorts (`r9_6`, `r9_r10_8`) comprise **CRISPRi** perturbations whereas `a549_ko` (`GSE345058`) is **CRISPRko**. Failure on `a549_ko` is confounded with CRISPRi->CRISPRko modality shift rather than purely same-modality failure; passage constitutes cross-modality generalization evidence. Under `K_conf = 1`, population generalizability across additional lineages remains explicitly limited.
- **Prior-Art & Data Attribution**: Baseline-aware perturbation transfer context acknowledges **Kawsar et al. (2026), `CanPert/CASTER` (`doi:10.64898/2026.11.15.766815`)**; `a549_ko` CRISPRko dataset derives from **Liu et al. (2026), `GSE345058`**; `hek293_xucao` CRISPRi derives from **Xu & Cao (2023), PerturbSci-Kinetics**.
- **Reporting-Only §5 Scope**: Legacy unpaired `H-D2`, paired `FCT - E1_ConstMed`, `N_5` null-cone suites (`all_nontarget` & `resp_top50_split`), `SA-Src8` (`--source-mode r9_r10_8`), `SA-Cas9` (`a549_ko_cas9pos`), naive `resp_top50`, and `H-D_util` are strictly descriptive and carry zero multiplicity weight.

## Round 13 (`PREREG_F`) Governance & Interpretation Rules

- **Protocol Inheritance**: Claims 1–2 (`H-E1`, `H-E2`), `r_orth_cal`, `all_nontarget` / `resp_top50_split`, `r9_6`, 2-claim Bonferroni + IUT (`97.5%` CIs), `4,000` bootstraps, and `100` permutations are inherited identically from locked `PREREG_E`; seed=`20261223`.
- **Separate-Round Independent Scope**: Round 13 has `K_conf <= 1` (`mac_ctrl` ONLY; `mac_lps` descriptive only). Confirmation in Round 13 MUST be reported strictly as independent prospective replication across two separately pre-registered rounds/lines without pooling Round 12 and Round 13 data; failure on any sub-criterion in `mac_ctrl` MUST be reported verbatim as non-replication on the second line.
- **Technical Boundaries**: Single donor (`RG5581`), 10x GEM-X Flex probe panel (`~18k` targeted probes, no mitochondrial cutoff), channel parity (`L01/L03` vs `L02/L04`).
