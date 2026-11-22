# Line-Specific Zero-Shot Perturbation Evaluation Report (`holdout_log1p_sens_s1_rep0`)

- **Code Freeze Commit (C1)**: `1ca9f60f49adf1ac7cce8266e40533273411db30` | **Bootstrap Seed (C2)**: `20261121`
- **Confirmatory Lines (3)**: nadig_hepg2, nadig_jurkat, vcc2025_h1
- **Secondary Lines (3)**: nourreddine_kolf21j, replogle_k562_essential, replogle_rpe1
- **Dual Confirmatory Summary (R13, 97.5% Bonferroni CIs)**: **Group A (Hybrid_ProjShared): 未确认 | Group B (Transfer_Mean_VarScaled): 未确认 (Bonferroni 97.5% CIs)**
- **Group A Primary (`Hybrid_ProjShared`)**: **未确认**
- **Group B Primary (`Transfer_Mean_VarScaled`)**: **未确认**
- **R17 Primary Correlation Metric**: `ls_resid_corr_orth_cal_top50` (amplitude-orthogonalized; raw `r_cal` reported as secondary)
- **R18 Global-Scale Scalar Baseline**: 单个标量 Transfer_Mean_GlobalScale 即可得到原始残差相关 [-0.205, +0.109]
- **R4 Co-fluctuation Structure Claim Allowed**: `目标线对照统计量` (`all beat SrcCorrTgtVar disjoint=False`)
- **R4 Marginal Variance Note**: `增益可由边际方差解释`
- **R14 Mechanism Statement (`Hybrid_on_TMVS` vs `TMVS` on `r_orth_cal`, 95% CI)**: `在本数据上，共波动结构没有超出边际方差的可检测增量` (`95%: 1/3`, `97.5%: 1/3`)
- **R20 Prior-Art Novelty Scope**: FCT is equivalent to VCC 2026 CellularAI (mu_t * 2^log2FC); CASTER performed supervised control-variance scaling. TMVS is framed strictly as a parameter-free variance-ratio variant of expression-aware transfer (not a new method); contribution is the pre-registered 3-layer decomposition of recoverable zero-shot line-specific residuals into (1) line-wide amplitude, (2) gene marginal variance, and (3) off-diagonal co-fluctuations.

## Hypothesis Adjudication Summary (Confirmatory Cohort, 97.5% CIs)

| Hypothesis | Group | Pass / Confirmatory | Threshold | Verdict |
|---|---|---:|---|---|
| **H-A1** (Slope 97.5% CI lo > 0.13) | A (`Hybrid_ProjShared`) | 0/3 | >= 3 | **未确认** |
| **H-A2** (Joint-boot `r_orth_cal` 97.5% CI lo > 0, valid>=95%, rho_AB_orth>0) | A (`Hybrid_ProjShared`) | 2/3 (fixed95: 2/3) | >= 3 | **未确认** |
| **H-A3** (Matched > ALL source mismatched cov on slope AND `r_orth_cal`) | A (`Hybrid_ProjShared`) | 0/3 | >= 3 | **未确认** |
| **H-A4** (Strict AND non-inferiority across 6 conventions) | A (`Hybrid_ProjShared`) | 3/3 inf (ctrl_disc_raw=3/3, ctrl_de50_r=2/3, ctrl_disc_raw_panel=3/3, cent_disc_raw=2/3, cent_de50_r=2/3, cent_disc_raw_panel=2/3; union=3) | <= 0 inf | **未确认** |
| **H-B1** (Slope 97.5% CI lo > 0.13) | B (`Transfer_Mean_VarScaled`) | 0/3 | >= 3 | **未确认** |
| **H-B2** (Joint-boot `r_orth_cal` 97.5% CI lo > 0, valid>=95%, rho_AB_orth>0) | B (`Transfer_Mean_VarScaled`) | 2/3 | >= 3 | **未确认** |
| **H-B3** (Matched `r_orth_cal` lo > max(ALL mismatched vars CI hi, 100-perm p97.5)) | B (`Transfer_Mean_VarScaled`) | 0/3 | >= 3 | **未确认** |
| **H-B4** (Strict AND non-inferiority across 6 conventions) | B (`Transfer_Mean_VarScaled`) | 3/3 inf (ctrl_disc_raw=1/3, ctrl_de50_r=1/3, ctrl_disc_raw_panel=2/3, cent_disc_raw=2/3, cent_de50_r=2/3, cent_disc_raw_panel=3/3; union=3) | <= 0 inf | **未确认** |

## Main Table: Primary View (A-T2 Own-Target-Zeroed Control Reference, 97.5% CIs, R5/R6/R11/R13/R17/R18)

| Cohort | Target Line | Method | N | Slope OLS (Top50) | 97.5% CI | r_orth_cal | r_orth_cal^2 | r_orth_cal^2 Joint 97.5% CI | Resid r_orth | r_cal (raw) | Resid r (raw) | sd_pred/sd_obs | Deming (delta=1) | Obs~Pred | Ctrl disc_raw | Ctrl de50_r | MSE |
|---|---|---|---:|---:|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Confirmatory | nadig_hepg2 | `Source_Grand_Mean` | 1777 | -0.0549 | [-0.0731, -0.0360] | +0.2072 | 0.0429 | [0.0340, 0.0529] | +0.1955 | -0.0935 | -0.0887 | 0.6195 | -0.0885 | -0.1432 | 0.5001 | 0.4287 | 0.0064 |
| Confirmatory | nadig_hepg2 | `Transfer_Mean` | 1777 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.7706 | 0.5535 | 0.0066 |
| Confirmatory | nadig_hepg2 | `CIPHER_Linear` | 1777 | -0.1179 | [-0.1380, -0.0977] | +0.0655 | 0.0043 | [0.0018, 0.0078] | +0.0618 | -0.1699 | -0.1611 | 0.7317 | -0.2392 | -0.2202 | 0.5001 | 0.1545 | 0.0071 |
| Confirmatory | nadig_hepg2 | `Hybrid_Raw` | 1777 | +0.0222 | [+0.0179, +0.0267] | +0.0580 | 0.0034 | [0.0012, 0.0068] | +0.0547 | +0.1430 | +0.1356 | 0.1640 | +0.0228 | +0.8273 | 0.7449 | 0.5441 | 0.0069 |
| Confirmatory | nadig_hepg2 | `Hybrid_ProjShared` | 1777 | +0.0062 | [+0.0041, +0.0085] | +0.0196 | 0.0004 | [0.0000, 0.0015] | +0.0185 | +0.0576 | +0.0547 | 0.1142 | +0.0063 | +0.4786 | 0.7577 | 0.5421 | 0.0068 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_MismatchedCov(nadig_jurkat)` | 1777 | -0.0002 | [-0.0007, +0.0003] | -0.0060 | 0.0000 | [0.0000, 0.0003] | -0.0057 | -0.0044 | -0.0041 | 0.0485 | -0.0002 | -0.0855 | 0.7695 | 0.5498 | 0.0067 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_MismatchedCov(nourreddine_kolf21j)` | 1777 | +0.0045 | [+0.0036, +0.0055] | +0.0352 | 0.0012 | [0.0006, 0.0021] | +0.0332 | +0.0544 | +0.0516 | 0.0882 | +0.0046 | +0.5846 | 0.7670 | 0.5480 | 0.0067 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 1777 | +0.0022 | [+0.0017, +0.0028] | +0.0223 | 0.0005 | [0.0002, 0.0010] | +0.0211 | +0.0432 | +0.0410 | 0.0547 | +0.0022 | +0.7499 | 0.7677 | 0.5524 | 0.0067 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 1777 | -0.0025 | [-0.0030, -0.0019] | -0.0212 | 0.0004 | [0.0001, 0.0010] | -0.0200 | -0.0433 | -0.0411 | 0.0602 | -0.0025 | -0.6824 | 0.7708 | 0.5463 | 0.0067 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_MismatchedCov(vcc2025_h1)` | 1777 | +0.0014 | [-0.0003, +0.0030] | +0.0115 | 0.0001 | [0.0000, 0.0005] | +0.0108 | +0.0104 | +0.0099 | 0.1437 | +0.0015 | +0.0690 | 0.7641 | 0.5331 | 0.0071 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_ShuffledGenes` | 1777 | +0.0043 | [+0.0031, +0.0056] | +0.0169 | 0.0003 | [0.0000, 0.0008] | +0.0160 | +0.0377 | +0.0358 | 0.1200 | +0.0044 | +0.2983 | 0.7681 | 0.5403 | 0.0071 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_DiagCov` | 1777 | +0.0047 | [+0.0036, +0.0059] | +0.0248 | 0.0006 | [0.0002, 0.0014] | +0.0234 | +0.0577 | +0.0547 | 0.0861 | +0.0047 | +0.6353 | 0.7685 | 0.5492 | 0.0067 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 1777 | +0.0030 | [+0.0024, +0.0035] | +0.0326 | 0.0011 | [0.0005, 0.0018] | +0.0307 | +0.0632 | +0.0600 | 0.0493 | +0.0030 | +1.2164 | 0.7677 | 0.5522 | 0.0066 |
| Confirmatory | nadig_hepg2 | `Transfer_Mean_VarScaled` | 1777 | +0.0424 | [+0.0377, +0.0472] | +0.0324 | 0.0011 | [0.0003, 0.0022] | +0.0306 | +0.1835 | +0.1740 | 0.2437 | +0.0450 | +0.7143 | 0.7645 | 0.5404 | 0.0064 |
| Confirmatory | nadig_hepg2 | `Transfer_Mean_GlobalScale` | 1777 | -0.0118 | [-0.0133, -0.0102] | +nan | nan | [nan, nan] | +nan | -0.2156 | -0.2046 | 0.0576 | -0.0118 | -3.5512 | 0.7706 | 0.5535 | 0.0064 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_MismatchedVar(nadig_jurkat)` | 1777 | +0.0144 | [+0.0124, +0.0165] | -0.0603 | 0.0036 | [0.0026, 0.0049] | -0.0569 | +0.1140 | +0.1081 | 0.1336 | +0.0147 | +0.8092 | 0.7698 | 0.5452 | 0.0071 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_MismatchedVar(nourreddine_kolf21j)` | 1777 | +0.0081 | [+0.0057, +0.0104] | +0.0484 | 0.0023 | [0.0015, 0.0034] | +0.0457 | +0.0549 | +0.0521 | 0.1564 | +0.0083 | +0.3331 | 0.7754 | 0.5503 | 0.0072 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 1777 | +0.0009 | [-0.0020, +0.0042] | +0.0866 | 0.0075 | [0.0055, 0.0100] | +0.0817 | +0.0077 | +0.0073 | 0.1245 | +0.0009 | +0.0587 | 0.7717 | 0.5534 | 0.0066 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 1777 | +0.0285 | [+0.0243, +0.0328] | -0.0366 | 0.0013 | [0.0004, 0.0030] | -0.0345 | +0.1628 | +0.1544 | 0.1846 | +0.0295 | +0.8365 | 0.7627 | 0.5405 | 0.0069 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_MismatchedVar(vcc2025_h1)` | 1777 | -0.0787 | [-0.0861, -0.0708] | -0.0160 | 0.0003 | [0.0000, 0.0011] | -0.0151 | -0.2258 | -0.2142 | 0.3671 | -0.0902 | -0.5836 | 0.7727 | 0.5286 | 0.0061 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_ShuffledVar` | 1777 | -0.0382 | [-0.0432, -0.0332] | +0.0505 | 0.0026 | [0.0017, 0.0035] | +0.0477 | -0.1449 | -0.1374 | 0.2783 | -0.0414 | -0.4938 | 0.7780 | 0.5287 | 0.0065 |
| Confirmatory | nadig_hepg2 | `Hybrid_on_TMVS` | 1777 | +0.0487 | [+0.0430, +0.0546] | +0.0354 | 0.0012 | [0.0004, 0.0026] | +0.0334 | +0.1770 | +0.1679 | 0.2898 | +0.0530 | +0.5794 | 0.7515 | 0.5286 | 0.0066 |
| Confirmatory | nadig_hepg2 | `FCT` | 1777 | +0.0408 | [+0.0377, +0.0438] | +0.2163 | 0.0468 | [0.0404, 0.0534] | +0.2041 | +0.2238 | +0.2123 | 0.1923 | +0.0423 | +1.1039 | 0.7770 | 0.5621 | 0.0065 |
| Confirmatory | nadig_jurkat | `Source_Grand_Mean` | 1779 | -0.0364 | [-0.0548, -0.0176] | +0.0181 | 0.0003 | [0.0000, 0.0016] | +0.0170 | -0.0587 | -0.0556 | 0.6547 | -0.0634 | -0.0849 | 0.5001 | 0.2621 | 0.0061 |
| Confirmatory | nadig_jurkat | `Transfer_Mean` | 1779 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.7927 | 0.5098 | 0.0069 |
| Confirmatory | nadig_jurkat | `CIPHER_Linear` | 1779 | -0.0433 | [-0.0642, -0.0214] | +0.0213 | 0.0005 | [0.0000, 0.0020] | +0.0201 | -0.0592 | -0.0561 | 0.7731 | -0.1065 | -0.0725 | 0.4936 | 0.0621 | 0.0062 |
| Confirmatory | nadig_jurkat | `Hybrid_Raw` | 1779 | +0.0105 | [+0.0079, +0.0134] | +0.1100 | 0.0121 | [0.0076, 0.0175] | +0.1038 | +0.1119 | +0.1060 | 0.0986 | +0.0106 | +1.0750 | 0.7917 | 0.5133 | 0.0069 |
| Confirmatory | nadig_jurkat | `Hybrid_ProjShared` | 1779 | +0.0040 | [+0.0025, +0.0056] | +0.0609 | 0.0037 | [0.0019, 0.0062] | +0.0575 | +0.0552 | +0.0523 | 0.0764 | +0.0040 | +0.6843 | 0.7851 | 0.5099 | 0.0070 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_MismatchedCov(nadig_hepg2)` | 1779 | -0.0034 | [-0.0044, -0.0025] | -0.0519 | 0.0027 | [0.0015, 0.0043] | -0.0490 | -0.0576 | -0.0546 | 0.0626 | -0.0034 | -0.8721 | 0.7878 | 0.5008 | 0.0070 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_MismatchedCov(nourreddine_kolf21j)` | 1779 | +0.0025 | [+0.0011, +0.0040] | +0.0231 | 0.0005 | [0.0001, 0.0014] | +0.0218 | +0.0314 | +0.0297 | 0.0856 | +0.0026 | +0.3469 | 0.7896 | 0.5043 | 0.0071 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 1779 | +0.0015 | [+0.0006, +0.0024] | +0.0216 | 0.0005 | [0.0001, 0.0012] | +0.0204 | +0.0272 | +0.0258 | 0.0575 | +0.0015 | +0.4486 | 0.7923 | 0.5082 | 0.0070 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 1779 | +0.0016 | [+0.0008, +0.0024] | +0.0357 | 0.0013 | [0.0005, 0.0024] | +0.0337 | +0.0318 | +0.0302 | 0.0534 | +0.0016 | +0.5648 | 0.7933 | 0.5091 | 0.0070 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_MismatchedCov(vcc2025_h1)` | 1779 | +0.0064 | [+0.0041, +0.0087] | +0.0370 | 0.0014 | [0.0005, 0.0026] | +0.0349 | +0.0447 | +0.0423 | 0.1509 | +0.0065 | +0.2804 | 0.7849 | 0.4961 | 0.0073 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_ShuffledGenes` | 1779 | +0.0039 | [+0.0023, +0.0056] | +0.0270 | 0.0007 | [0.0002, 0.0017] | +0.0255 | +0.0429 | +0.0407 | 0.0964 | +0.0040 | +0.4219 | 0.7890 | 0.5029 | 0.0072 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_DiagCov` | 1779 | +0.0037 | [+0.0021, +0.0053] | +0.0305 | 0.0009 | [0.0002, 0.0022] | +0.0288 | +0.0473 | +0.0448 | 0.0818 | +0.0037 | +0.5475 | 0.7904 | 0.5063 | 0.0070 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_SrcCorrTgtVar` | 1779 | +0.0019 | [+0.0013, +0.0025] | +0.0448 | 0.0020 | [0.0010, 0.0033] | +0.0423 | +0.0472 | +0.0448 | 0.0432 | +0.0019 | +1.0359 | 0.7907 | 0.5100 | 0.0070 |
| Confirmatory | nadig_jurkat | `Transfer_Mean_VarScaled` | 1779 | +0.0152 | [+0.0114, +0.0190] | +0.0874 | 0.0076 | [0.0059, 0.0096] | +0.0824 | +0.0856 | +0.0811 | 0.1879 | +0.0158 | +0.4318 | 0.7926 | 0.5108 | 0.0077 |
| Confirmatory | nadig_jurkat | `Transfer_Mean_GlobalScale` | 1779 | +0.0066 | [+0.0032, +0.0098] | +nan | nan | [nan, nan] | +nan | +0.0576 | +0.0546 | 0.1205 | +0.0067 | +0.4528 | 0.7927 | 0.5098 | 0.0078 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_MismatchedVar(nadig_hepg2)` | 1779 | -0.0267 | [-0.0306, -0.0228] | -0.1656 | 0.0274 | [0.0234, 0.0319] | -0.1563 | -0.1586 | -0.1503 | 0.1774 | -0.0275 | -0.8471 | 0.7893 | 0.4858 | 0.0071 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_MismatchedVar(nourreddine_kolf21j)` | 1779 | +0.0365 | [+0.0330, +0.0403] | +0.1458 | 0.0213 | [0.0180, 0.0247] | +0.1376 | +0.1914 | +0.1814 | 0.2013 | +0.0380 | +0.9011 | 0.7978 | 0.5194 | 0.0078 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 1779 | +0.0073 | [+0.0040, +0.0105] | +0.0982 | 0.0096 | [0.0072, 0.0125] | +0.0927 | +0.0573 | +0.0543 | 0.1337 | +0.0074 | +0.4061 | 0.7943 | 0.5127 | 0.0069 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 1779 | -0.0199 | [-0.0238, -0.0160] | -0.1219 | 0.0149 | [0.0115, 0.0190] | -0.1150 | -0.1318 | -0.1249 | 0.1596 | -0.0204 | -0.7824 | 0.7862 | 0.4907 | 0.0077 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_MismatchedVar(vcc2025_h1)` | 1779 | -0.0026 | [-0.0114, +0.0061] | +0.0337 | 0.0011 | [0.0003, 0.0027] | +0.0318 | -0.0072 | -0.0069 | 0.3847 | -0.0031 | -0.0178 | 0.7943 | 0.4821 | 0.0057 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_ShuffledVar` | 1779 | +0.0287 | [+0.0246, +0.0328] | +0.1144 | 0.0131 | [0.0110, 0.0153] | +0.1079 | +0.1183 | +0.1121 | 0.2561 | +0.0307 | +0.4376 | 0.7989 | 0.5058 | 0.0073 |
| Confirmatory | nadig_jurkat | `Hybrid_on_TMVS` | 1779 | +0.0192 | [+0.0153, +0.0231] | +0.1026 | 0.0105 | [0.0081, 0.0134] | +0.0968 | +0.0983 | +0.0932 | 0.2065 | +0.0201 | +0.4512 | 0.7863 | 0.5114 | 0.0078 |
| Confirmatory | nadig_jurkat | `FCT` | 1779 | +0.0149 | [+0.0117, +0.0181] | +0.0853 | 0.0073 | [0.0052, 0.0097] | +0.0805 | +0.0826 | +0.0782 | 0.1902 | +0.0154 | +0.4113 | 0.7878 | 0.5047 | 0.0078 |
| Secondary | nourreddine_kolf21j | `Source_Grand_Mean` | 28 | +0.1322 | [+0.0107, +0.2694] | +0.0515 | 0.0026 | [0.0000, 0.0221] | +0.0485 | +0.2137 | +0.2024 | 0.6533 | +0.2195 | +0.3098 | 0.5000 | 0.2327 | 0.0027 |
| Secondary | nourreddine_kolf21j | `Transfer_Mean` | 28 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.7130 | 0.3458 | 0.0036 |
| Secondary | nourreddine_kolf21j | `CIPHER_Linear` | 28 | +0.1340 | [-0.0044, +0.2848] | -0.0985 | 0.0097 | [0.0000, 0.0437] | -0.0928 | +0.1899 | +0.1798 | 0.7449 | +0.2778 | +0.2414 | 0.5251 | -0.0927 | 0.0020 |
| Secondary | nourreddine_kolf21j | `Hybrid_Raw` | 28 | -0.0046 | [-0.0280, +0.0173] | +0.0096 | 0.0001 | [0.0000, 0.0115] | +0.0091 | -0.0250 | -0.0237 | 0.1953 | -0.0048 | -0.1212 | 0.7183 | 0.3447 | 0.0037 |
| Secondary | nourreddine_kolf21j | `Hybrid_ProjShared` | 28 | -0.0072 | [-0.0310, +0.0147] | +0.0131 | 0.0002 | [0.0000, 0.0137] | +0.0124 | -0.0457 | -0.0433 | 0.1672 | -0.0075 | -0.2590 | 0.7235 | 0.3470 | 0.0038 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_MismatchedCov(nadig_hepg2)` | 28 | -0.0059 | [-0.0313, +0.0180] | +0.0278 | 0.0008 | [0.0000, 0.0173] | +0.0262 | -0.0393 | -0.0372 | 0.1596 | -0.0061 | -0.2333 | 0.7050 | 0.3509 | 0.0038 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_MismatchedCov(nadig_jurkat)` | 28 | -0.0070 | [-0.0331, +0.0151] | +0.0226 | 0.0005 | [0.0000, 0.0152] | +0.0213 | -0.0491 | -0.0465 | 0.1505 | -0.0072 | -0.3090 | 0.6997 | 0.3487 | 0.0038 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 28 | -0.0079 | [-0.0374, +0.0185] | +0.0287 | 0.0008 | [0.0000, 0.0195] | +0.0271 | -0.0512 | -0.0485 | 0.1627 | -0.0081 | -0.2978 | 0.7077 | 0.3541 | 0.0038 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 28 | -0.0138 | [-0.0432, +0.0072] | -0.0170 | 0.0003 | [0.0000, 0.0091] | -0.0160 | -0.0929 | -0.0880 | 0.1571 | -0.0142 | -0.5600 | 0.6958 | 0.3397 | 0.0038 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_MismatchedCov(vcc2025_h1)` | 28 | +0.0093 | [-0.0260, +0.0400] | +0.0318 | 0.0010 | [0.0000, 0.0425] | +0.0299 | +0.0476 | +0.0451 | 0.2066 | +0.0097 | +0.2182 | 0.7183 | 0.3427 | 0.0039 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_ShuffledGenes` | 28 | -0.0043 | [-0.0334, +0.0233] | +0.0479 | 0.0023 | [0.0000, 0.0224] | +0.0451 | -0.0247 | -0.0234 | 0.1841 | -0.0045 | -0.1270 | 0.7130 | 0.3590 | 0.0039 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_DiagCov` | 28 | -0.0044 | [-0.0340, +0.0234] | +0.0548 | 0.0030 | [0.0000, 0.0276] | +0.0517 | -0.0273 | -0.0259 | 0.1708 | -0.0045 | -0.1515 | 0.7037 | 0.3605 | 0.0038 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_SrcCorrTgtVar` | 28 | -0.0027 | [-0.0298, +0.0204] | +0.0639 | 0.0041 | [0.0000, 0.0393] | +0.0602 | -0.0216 | -0.0205 | 0.1304 | -0.0027 | -0.1570 | 0.7116 | 0.3572 | 0.0037 |
| Secondary | nourreddine_kolf21j | `Transfer_Mean_VarScaled` | 28 | -0.0205 | [-0.0660, +0.0238] | +0.1012 | 0.0102 | [0.0000, 0.0473] | +0.0953 | -0.0742 | -0.0703 | 0.2923 | -0.0225 | -0.2405 | 0.7130 | 0.3523 | 0.0046 |
| Secondary | nourreddine_kolf21j | `Transfer_Mean_GlobalScale` | 28 | -0.0339 | [-0.0702, -0.0018] | +nan | nan | [nan, nan] | +nan | -0.2009 | -0.1902 | 0.1785 | -0.0350 | -1.0657 | 0.7130 | 0.3458 | 0.0048 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_MismatchedVar(nadig_hepg2)` | 28 | -0.0006 | [-0.0104, +0.0141] | -0.0425 | 0.0018 | [0.0000, 0.0211] | -0.0400 | -0.0052 | -0.0050 | 0.1130 | -0.0006 | -0.0439 | 0.7156 | 0.3361 | 0.0037 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_MismatchedVar(nadig_jurkat)` | 28 | -0.0210 | [-0.0419, +0.0024] | +0.0351 | 0.0012 | [0.0000, 0.0088] | +0.0330 | -0.1439 | -0.1363 | 0.1540 | -0.0215 | -0.8846 | 0.7037 | 0.3473 | 0.0044 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 28 | +0.0060 | [-0.0030, +0.0151] | +0.0969 | 0.0094 | [0.0018, 0.0231] | +0.0914 | +0.0602 | +0.0570 | 0.1049 | +0.0060 | +0.5434 | 0.7130 | 0.3505 | 0.0037 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 28 | +0.0042 | [-0.0087, +0.0235] | +0.0121 | 0.0001 | [0.0000, 0.0105] | +0.0114 | +0.0311 | +0.0295 | 0.1409 | +0.0042 | +0.2092 | 0.7116 | 0.3451 | 0.0041 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_MismatchedVar(vcc2025_h1)` | 28 | +0.0293 | [-0.0096, +0.0628] | -0.0841 | 0.0071 | [0.0000, 0.0262] | -0.0793 | +0.0971 | +0.0919 | 0.3190 | +0.0326 | +0.2882 | 0.7011 | 0.3145 | 0.0024 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_ShuffledVar` | 28 | +0.0110 | [-0.0110, +0.0343] | +0.0274 | 0.0007 | [0.0000, 0.0095] | +0.0258 | +0.0451 | +0.0427 | 0.2587 | +0.0118 | +0.1651 | 0.7090 | 0.3351 | 0.0042 |
| Secondary | nourreddine_kolf21j | `Hybrid_on_TMVS` | 28 | -0.0278 | [-0.0867, +0.0264] | +0.0876 | 0.0077 | [0.0000, 0.0379] | +0.0825 | -0.0834 | -0.0789 | 0.3521 | -0.0317 | -0.2242 | 0.7302 | 0.3546 | 0.0048 |
| Secondary | nourreddine_kolf21j | `FCT` | 28 | +0.0304 | [-0.0035, +0.0699] | +0.1101 | 0.0121 | [0.0000, 0.0503] | +0.1037 | +0.1160 | +0.1099 | 0.2764 | +0.0328 | +0.3974 | 0.7130 | 0.3418 | 0.0040 |
| Secondary | replogle_k562_essential | `Source_Grand_Mean` | 1550 | +0.1509 | [+0.1195, +0.1845] | -0.0486 | 0.0024 | [0.0009, 0.0044] | -0.0465 | +0.1668 | +0.1603 | 0.9413 | +0.6913 | +0.1703 | 0.5001 | 0.2778 | 0.0041 |
| Secondary | replogle_k562_essential | `Transfer_Mean` | 1550 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.7920 | 0.5364 | 0.0053 |
| Secondary | replogle_k562_essential | `CIPHER_Linear` | 1550 | +0.1963 | [+0.1611, +0.2338] | +0.0168 | 0.0003 | [0.0000, 0.0014] | +0.0161 | +0.1848 | +0.1777 | 1.1051 | +1.7114 | +0.1608 | 0.4583 | 0.0479 | 0.0038 |
| Secondary | replogle_k562_essential | `Hybrid_Raw` | 1550 | +0.0140 | [+0.0113, +0.0169] | +0.0847 | 0.0072 | [0.0044, 0.0105] | +0.0810 | +0.1082 | +0.1040 | 0.1350 | +0.0143 | +0.7705 | 0.8017 | 0.5367 | 0.0052 |
| Secondary | replogle_k562_essential | `Hybrid_ProjShared` | 1550 | +0.0029 | [+0.0013, +0.0044] | +0.0372 | 0.0014 | [0.0004, 0.0029] | +0.0356 | +0.0320 | +0.0308 | 0.0932 | +0.0029 | +0.3302 | 0.7898 | 0.5351 | 0.0054 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(nadig_hepg2)` | 1550 | -0.0055 | [-0.0066, -0.0045] | -0.0695 | 0.0048 | [0.0033, 0.0068] | -0.0665 | -0.0720 | -0.0692 | 0.0797 | -0.0056 | -0.8683 | 0.7909 | 0.5295 | 0.0054 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(nadig_jurkat)` | 1550 | +0.0010 | [+0.0001, +0.0019] | +0.0127 | 0.0002 | [0.0000, 0.0006] | +0.0121 | +0.0157 | +0.0151 | 0.0672 | +0.0010 | +0.2247 | 0.7929 | 0.5342 | 0.0054 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(nourreddine_kolf21j)` | 1550 | -0.0008 | [-0.0024, +0.0009] | +0.0045 | 0.0000 | [0.0000, 0.0003] | +0.0043 | -0.0069 | -0.0067 | 0.1165 | -0.0008 | -0.0571 | 0.7940 | 0.5290 | 0.0055 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 1550 | +0.0040 | [+0.0030, +0.0049] | +0.0668 | 0.0045 | [0.0028, 0.0065] | +0.0639 | +0.0624 | +0.0600 | 0.0661 | +0.0040 | +0.9070 | 0.7899 | 0.5362 | 0.0054 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(vcc2025_h1)` | 1550 | +0.0058 | [+0.0029, +0.0086] | +0.0222 | 0.0005 | [0.0001, 0.0013] | +0.0213 | +0.0282 | +0.0271 | 0.2136 | +0.0061 | +0.1269 | 0.7847 | 0.5205 | 0.0057 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_ShuffledGenes` | 1550 | -0.0021 | [-0.0038, -0.0002] | +0.0024 | 0.0000 | [0.0000, 0.0003] | +0.0023 | -0.0171 | -0.0164 | 0.1254 | -0.0021 | -0.1309 | 0.7948 | 0.5290 | 0.0056 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_DiagCov` | 1550 | -0.0020 | [-0.0037, -0.0003] | +0.0022 | 0.0000 | [0.0000, 0.0003] | +0.0021 | -0.0186 | -0.0179 | 0.1119 | -0.0020 | -0.1599 | 0.7947 | 0.5307 | 0.0054 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_SrcCorrTgtVar` | 1550 | +0.0019 | [+0.0011, +0.0027] | +0.0300 | 0.0009 | [0.0003, 0.0019] | +0.0287 | +0.0343 | +0.0329 | 0.0568 | +0.0019 | +0.5802 | 0.7910 | 0.5356 | 0.0053 |
| Secondary | replogle_k562_essential | `Transfer_Mean_VarScaled` | 1550 | +0.0451 | [+0.0383, +0.0529] | +0.1711 | 0.0293 | [0.0237, 0.0357] | +0.1637 | +0.2174 | +0.2090 | 0.2156 | +0.0472 | +0.9691 | 0.7928 | 0.5454 | 0.0051 |
| Secondary | replogle_k562_essential | `Transfer_Mean_GlobalScale` | 1550 | +0.0016 | [+0.0013, +0.0019] | +nan | nan | [nan, nan] | +nan | +0.1851 | +0.1779 | 0.0091 | +0.0016 | +19.5735 | 0.7920 | 0.5364 | 0.0053 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(nadig_hepg2)` | 1550 | -0.0310 | [-0.0394, -0.0221] | -0.1571 | 0.0247 | [0.0180, 0.0315] | -0.1503 | -0.1440 | -0.1384 | 0.2241 | -0.0326 | -0.6175 | 0.7891 | 0.5162 | 0.0053 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(nadig_jurkat)` | 1550 | -0.0067 | [-0.0111, -0.0027] | +0.0949 | 0.0090 | [0.0067, 0.0118] | +0.0908 | -0.0357 | -0.0343 | 0.1957 | -0.0070 | -0.1751 | 0.7926 | 0.5383 | 0.0062 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(nourreddine_kolf21j)` | 1550 | +0.0209 | [+0.0134, +0.0276] | +0.1430 | 0.0205 | [0.0160, 0.0253] | +0.1368 | +0.0832 | +0.0800 | 0.2608 | +0.0224 | +0.3068 | 0.7948 | 0.5440 | 0.0062 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 1550 | -0.0136 | [-0.0214, -0.0044] | -0.0781 | 0.0061 | [0.0029, 0.0098] | -0.0747 | -0.0712 | -0.0685 | 0.1981 | -0.0141 | -0.3455 | 0.7870 | 0.5230 | 0.0059 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(vcc2025_h1)` | 1550 | +0.0748 | [+0.0595, +0.0907] | -0.0090 | 0.0001 | [0.0000, 0.0007] | -0.0086 | +0.1382 | +0.1328 | 0.5631 | +0.1082 | +0.2358 | 0.7966 | 0.5053 | 0.0036 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_ShuffledVar` | 1550 | +0.0620 | [+0.0534, +0.0710] | +0.0564 | 0.0032 | [0.0020, 0.0046] | +0.0540 | +0.1722 | +0.1655 | 0.3748 | +0.0718 | +0.4415 | 0.7958 | 0.5210 | 0.0048 |
| Secondary | replogle_k562_essential | `Hybrid_on_TMVS` | 1550 | +0.0479 | [+0.0408, +0.0560] | +0.1683 | 0.0283 | [0.0227, 0.0349] | +0.1610 | +0.2114 | +0.2031 | 0.2359 | +0.0506 | +0.8610 | 0.7905 | 0.5438 | 0.0052 |
| Secondary | replogle_k562_essential | `FCT` | 1550 | +0.0285 | [+0.0246, +0.0324] | +0.1468 | 0.0216 | [0.0175, 0.0261] | +0.1405 | +0.1198 | +0.1152 | 0.2475 | +0.0303 | +0.4654 | 0.7847 | 0.5397 | 0.0057 |
| Secondary | replogle_rpe1 | `Source_Grand_Mean` | 1779 | -0.0936 | [-0.1052, -0.0826] | +0.5109 | 0.2610 | [0.2437, 0.2777] | +0.4957 | -0.2389 | -0.2326 | 0.4025 | -0.1104 | -0.5780 | 0.5003 | 0.5336 | 0.0101 |
| Secondary | replogle_rpe1 | `Transfer_Mean` | 1779 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.8612 | 0.5467 | 0.0097 |
| Secondary | replogle_rpe1 | `CIPHER_Linear` | 1779 | -0.1493 | [-0.1615, -0.1377] | +0.3315 | 0.1099 | [0.0929, 0.1276] | +0.3217 | -0.3401 | -0.3312 | 0.4508 | -0.1812 | -0.7348 | 0.5219 | 0.2579 | 0.0112 |
| Secondary | replogle_rpe1 | `Hybrid_Raw` | 1779 | +0.0149 | [+0.0132, +0.0166] | +0.2586 | 0.0669 | [0.0533, 0.0812] | +0.2509 | +0.2572 | +0.2505 | 0.0596 | +0.0150 | +4.2048 | 0.8493 | 0.5662 | 0.0095 |
| Secondary | replogle_rpe1 | `Hybrid_ProjShared` | 1779 | +0.0071 | [+0.0063, +0.0079] | +0.1811 | 0.0328 | [0.0263, 0.0396] | +0.1758 | +0.1613 | +0.1571 | 0.0451 | +0.0071 | +3.4817 | 0.8537 | 0.5563 | 0.0096 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(nadig_hepg2)` | 1779 | +0.0048 | [+0.0042, +0.0054] | +0.1628 | 0.0265 | [0.0207, 0.0328] | +0.1579 | +0.1614 | +0.1572 | 0.0305 | +0.0048 | +5.1585 | 0.8598 | 0.5535 | 0.0096 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(nadig_jurkat)` | 1779 | -0.0021 | [-0.0025, -0.0018] | -0.0793 | 0.0063 | [0.0038, 0.0094] | -0.0770 | -0.0987 | -0.0961 | 0.0219 | -0.0021 | -4.3790 | 0.8606 | 0.5414 | 0.0097 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(nourreddine_kolf21j)` | 1779 | -0.0049 | [-0.0056, -0.0043] | -0.1360 | 0.0185 | [0.0147, 0.0229] | -0.1320 | -0.1253 | -0.1221 | 0.0404 | -0.0049 | -3.0221 | 0.8616 | 0.5321 | 0.0098 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 1779 | -0.0026 | [-0.0030, -0.0022] | -0.1293 | 0.0167 | [0.0125, 0.0215] | -0.1254 | -0.1111 | -0.1082 | 0.0240 | -0.0026 | -4.5120 | 0.8603 | 0.5387 | 0.0098 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(vcc2025_h1)` | 1779 | -0.0059 | [-0.0067, -0.0052] | -0.1063 | 0.0113 | [0.0089, 0.0141] | -0.1032 | -0.1028 | -0.1001 | 0.0591 | -0.0059 | -1.6932 | 0.8590 | 0.5249 | 0.0099 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_ShuffledGenes` | 1779 | -0.0054 | [-0.0062, -0.0048] | -0.1406 | 0.0198 | [0.0160, 0.0240] | -0.1364 | -0.1182 | -0.1152 | 0.0473 | -0.0055 | -2.4360 | 0.8614 | 0.5286 | 0.0099 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_DiagCov` | 1779 | -0.0046 | [-0.0053, -0.0040] | -0.1471 | 0.0216 | [0.0172, 0.0264] | -0.1427 | -0.1165 | -0.1134 | 0.0410 | -0.0047 | -2.7686 | 0.8612 | 0.5317 | 0.0098 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 1779 | -0.0009 | [-0.0012, -0.0007] | -0.0324 | 0.0011 | [0.0003, 0.0023] | -0.0315 | -0.0487 | -0.0474 | 0.0200 | -0.0009 | -2.3727 | 0.8601 | 0.5445 | 0.0097 |
| Secondary | replogle_rpe1 | `Transfer_Mean_VarScaled` | 1779 | +0.0311 | [+0.0285, +0.0339] | +0.0992 | 0.0098 | [0.0073, 0.0128] | +0.0963 | +0.2700 | +0.2630 | 0.1183 | +0.0315 | +2.2223 | 0.8556 | 0.5476 | 0.0096 |
| Secondary | replogle_rpe1 | `Transfer_Mean_GlobalScale` | 1779 | +0.0077 | [+0.0073, +0.0082] | +nan | nan | [nan, nan] | +nan | +0.4067 | +0.3961 | 0.0196 | +0.0078 | +20.2518 | 0.8612 | 0.5467 | 0.0097 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(nadig_hepg2)` | 1779 | +0.0237 | [+0.0216, +0.0258] | +0.0183 | 0.0003 | [0.0000, 0.0010] | +0.0177 | +0.2091 | +0.2037 | 0.1162 | +0.0240 | +1.7519 | 0.8566 | 0.5438 | 0.0096 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(nadig_jurkat)` | 1779 | +0.0205 | [+0.0190, +0.0221] | -0.0398 | 0.0016 | [0.0009, 0.0025] | -0.0386 | +0.2408 | +0.2345 | 0.0876 | +0.0207 | +2.6779 | 0.8602 | 0.5353 | 0.0099 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(nourreddine_kolf21j)` | 1779 | +0.0097 | [+0.0085, +0.0110] | +0.0146 | 0.0002 | [0.0000, 0.0007] | +0.0142 | +0.0992 | +0.0966 | 0.1005 | +0.0098 | +0.9619 | 0.8636 | 0.5336 | 0.0100 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 1779 | -0.0122 | [-0.0136, -0.0108] | -0.0248 | 0.0006 | [0.0002, 0.0012] | -0.0240 | -0.1650 | -0.1607 | 0.0758 | -0.0122 | -2.1202 | 0.8604 | 0.5317 | 0.0098 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(vcc2025_h1)` | 1779 | -0.0770 | [-0.0818, -0.0721] | +0.0205 | 0.0004 | [0.0001, 0.0011] | +0.0199 | -0.3566 | -0.3473 | 0.2216 | -0.0804 | -1.5669 | 0.8629 | 0.5323 | 0.0100 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_ShuffledVar` | 1779 | -0.0353 | [-0.0378, -0.0329] | +0.0256 | 0.0007 | [0.0003, 0.0012] | +0.0248 | -0.2317 | -0.2257 | 0.1563 | -0.0361 | -1.4439 | 0.8669 | 0.5173 | 0.0102 |
| Secondary | replogle_rpe1 | `Hybrid_on_TMVS` | 1779 | +0.0382 | [+0.0354, +0.0412] | +0.1609 | 0.0259 | [0.0212, 0.0311] | +0.1561 | +0.2975 | +0.2897 | 0.1319 | +0.0388 | +2.1964 | 0.8492 | 0.5544 | 0.0094 |
| Secondary | replogle_rpe1 | `FCT` | 1779 | +0.0183 | [+0.0167, +0.0201] | +0.1295 | 0.0168 | [0.0135, 0.0203] | +0.1257 | +0.1684 | +0.1640 | 0.1119 | +0.0186 | +1.4660 | 0.8625 | 0.5463 | 0.0096 |
| Confirmatory | vcc2025_h1 | `Source_Grand_Mean` | 95 | +0.0836 | [-0.0669, +0.2472] | +0.0826 | 0.0068 | [0.0006, 0.0200] | +0.0824 | +0.1526 | +0.1523 | 0.5487 | +0.1179 | +0.2776 | 0.4999 | 0.2164 | 0.0038 |
| Confirmatory | vcc2025_h1 | `Transfer_Mean` | 95 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.8917 | 0.4842 | 0.0040 |
| Confirmatory | vcc2025_h1 | `CIPHER_Linear` | 95 | +0.0465 | [-0.1124, +0.2189] | -0.1236 | 0.0153 | [0.0001, 0.0724] | -0.1233 | +0.0788 | +0.0786 | 0.5917 | +0.0712 | +0.1329 | 0.4614 | -0.0083 | 0.0037 |
| Confirmatory | vcc2025_h1 | `Hybrid_Raw` | 95 | -0.0136 | [-0.0309, +0.0009] | -0.1011 | 0.0102 | [0.0000, 0.0488] | -0.1009 | -0.0923 | -0.0921 | 0.1480 | -0.0139 | -0.6223 | 0.8433 | 0.4255 | 0.0044 |
| Confirmatory | vcc2025_h1 | `Hybrid_ProjShared` | 95 | -0.0141 | [-0.0316, +0.0002] | -0.1016 | 0.0103 | [0.0000, 0.0519] | -0.1013 | -0.0986 | -0.0984 | 0.1432 | -0.0144 | -0.6870 | 0.8419 | 0.4307 | 0.0044 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_MismatchedCov(nadig_hepg2)` | 95 | -0.0055 | [-0.0132, -0.0000] | -0.1065 | 0.0114 | [0.0001, 0.0572] | -0.1062 | -0.0921 | -0.0919 | 0.0595 | -0.0055 | -1.5432 | 0.8856 | 0.4692 | 0.0041 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_MismatchedCov(nadig_jurkat)` | 95 | -0.0002 | [-0.0025, +0.0020] | -0.0053 | 0.0000 | [0.0000, 0.0035] | -0.0053 | -0.0039 | -0.0039 | 0.0429 | -0.0002 | -0.0903 | 0.8795 | 0.4765 | 0.0040 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_MismatchedCov(nourreddine_kolf21j)` | 95 | -0.0014 | [-0.0046, +0.0021] | -0.0092 | 0.0001 | [0.0000, 0.0043] | -0.0092 | -0.0246 | -0.0245 | 0.0569 | -0.0014 | -0.4313 | 0.8797 | 0.4744 | 0.0041 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 95 | +0.0026 | [-0.0008, +0.0066] | +0.0747 | 0.0056 | [0.0000, 0.0274] | +0.0745 | +0.0601 | +0.0600 | 0.0429 | +0.0026 | +1.3995 | 0.8841 | 0.4821 | 0.0040 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 95 | +0.0033 | [-0.0012, +0.0098] | +0.0925 | 0.0086 | [0.0000, 0.0512] | +0.0922 | +0.0682 | +0.0680 | 0.0488 | +0.0033 | +1.3947 | 0.8940 | 0.4861 | 0.0040 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_ShuffledGenes` | 95 | +0.0027 | [-0.0031, +0.0101] | +0.0703 | 0.0049 | [0.0000, 0.0228] | +0.0701 | +0.0305 | +0.0304 | 0.0873 | +0.0027 | +0.3483 | 0.8885 | 0.4749 | 0.0044 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_DiagCov` | 95 | +0.0026 | [-0.0024, +0.0092] | +0.1012 | 0.0102 | [0.0001, 0.0440] | +0.1009 | +0.0454 | +0.0453 | 0.0565 | +0.0026 | +0.8016 | 0.8877 | 0.4835 | 0.0040 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 95 | -0.0004 | [-0.0079, +0.0046] | -0.0580 | 0.0034 | [0.0000, 0.0467] | -0.0579 | -0.0080 | -0.0080 | 0.0445 | -0.0004 | -0.1805 | 0.8832 | 0.4760 | 0.0040 |
| Confirmatory | vcc2025_h1 | `Transfer_Mean_VarScaled` | 95 | -0.1303 | [-0.2156, -0.0544] | -0.0474 | 0.0022 | [0.0000, 0.0246] | -0.0473 | -0.4112 | -0.4103 | 0.3175 | -0.1420 | -1.2922 | 0.8910 | 0.4832 | 0.0034 |
| Confirmatory | vcc2025_h1 | `Transfer_Mean_GlobalScale` | 95 | +0.0257 | [-0.0394, +0.0960] | +nan | nan | [nan, nan] | +nan | +0.1097 | +0.1094 | 0.2353 | +0.0272 | +0.4650 | 0.8917 | 0.4842 | 0.0032 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_MismatchedVar(nadig_hepg2)` | 95 | +0.0203 | [-0.0049, +0.0471] | +0.0436 | 0.0019 | [0.0000, 0.0121] | +0.0434 | +0.2057 | +0.2052 | 0.0987 | +0.0205 | +2.0789 | 0.8927 | 0.4824 | 0.0037 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_MismatchedVar(nadig_jurkat)` | 95 | -0.0260 | [-0.0551, +0.0003] | -0.0563 | 0.0032 | [0.0001, 0.0098] | -0.0562 | -0.2593 | -0.2587 | 0.1005 | -0.0262 | -2.5752 | 0.8901 | 0.4777 | 0.0042 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_MismatchedVar(nourreddine_kolf21j)` | 95 | -0.0254 | [-0.0638, +0.0107] | +0.0206 | 0.0004 | [0.0000, 0.0063] | +0.0205 | -0.1783 | -0.1779 | 0.1428 | -0.0259 | -1.2457 | 0.8897 | 0.4835 | 0.0043 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 95 | +0.0199 | [-0.0047, +0.0465] | +0.0385 | 0.0015 | [0.0000, 0.0065] | +0.0384 | +0.2116 | +0.2111 | 0.0940 | +0.0200 | +2.2450 | 0.8905 | 0.4811 | 0.0037 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 95 | +0.0227 | [-0.0064, +0.0542] | -0.0472 | 0.0022 | [0.0000, 0.0128] | -0.0471 | +0.2070 | +0.2065 | 0.1101 | +0.0230 | +1.8763 | 0.8931 | 0.4711 | 0.0039 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_ShuffledVar` | 95 | +0.0140 | [-0.0502, +0.0830] | +0.0302 | 0.0009 | [0.0000, 0.0071] | +0.0301 | +0.0514 | +0.0513 | 0.2730 | +0.0151 | +0.1878 | 0.8952 | 0.4811 | 0.0032 |
| Confirmatory | vcc2025_h1 | `Hybrid_on_TMVS` | 95 | -0.1444 | [-0.2311, -0.0689] | -0.1037 | 0.0108 | [0.0000, 0.0464] | -0.1035 | -0.4062 | -0.4053 | 0.3562 | -0.1611 | -1.1378 | 0.8152 | 0.3920 | 0.0038 |
| Confirmatory | vcc2025_h1 | `FCT` | 95 | +0.2017 | [+0.1201, +0.2910] | +0.3618 | 0.1309 | [0.0806, 0.1899] | +0.3608 | +0.5279 | +0.5267 | 0.3829 | +0.2245 | +1.3755 | 0.9047 | 0.5381 | 0.0049 |

## R18: TMVS Shuffled-Variance 100-Permutation Null vs H-B3 (97.5% CI Threshold)

| Target Line | TMVS r_orth_cal | TMVS CI lo | Perm100 p2.5 | Perm100 median | Perm100 p97.5 | Perm100 max | Max Mismatch CI hi | H-B3 Pass |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| nadig_hepg2 | +0.0324 | +0.0185 | +0.0059 | +0.0445 | +0.0746 | +0.0918 | +0.0998 | False |
| nadig_jurkat | +0.0874 | +0.0768 | +0.0521 | +0.1011 | +0.1382 | +0.1573 | +0.1571 | False |
| nourreddine_kolf21j | +0.1012 | -0.0274 | +0.0170 | +0.0696 | +0.1264 | +0.1611 | +0.1520 | False |
| replogle_k562_essential | +0.1711 | +0.1541 | +0.0238 | +0.0754 | +0.1172 | +0.1254 | +0.1590 | False |
| replogle_rpe1 | +0.0992 | +0.0855 | -0.0191 | +0.0169 | +0.0462 | +0.0509 | +0.0333 | True |
| vcc2025_h1 | -0.0474 | -0.1568 | -0.0185 | +0.0256 | +0.0510 | +0.0785 | +0.1101 | False |

## R19: Descriptive 3-Layer Residual Variance Decomposition (`R = Delta_t - Delta_src`, PREREG_A v4 §5.4)

| Target Line | Mask | N | Elig Genes | beta_ols_int [95% CI] | beta_mult | L1 Line Gain [95% CI] | L2 Per-Gene Gain [95% CI] | L3 Remaining [95% CI] |
|---|---|---:|---:|---|---:|---|---|---|
| nadig_hepg2 | `resp_top50` | 1777 | 2096 (C_dist>=5) | 1.2708 [1.2316, 1.3098] | 1.2923 | 0.0419 [0.0320, 0.0530] | 0.1227 [0.1115, 0.1370] | 0.8354 [0.8192, 0.8490] |
| nadig_hepg2 | `all_nontarget` | 1777 | 6148 (C_dist>=5) | 0.5650 [0.5445, 0.5857] | 0.5705 | 0.1070 [0.0966, 0.1178] | 0.0518 [0.0476, 0.0559] | 0.8412 [0.8330, 0.8493] |
| nadig_jurkat | `resp_top50` | 1779 | 2256 (C_dist>=5) | 1.0702 [1.0380, 1.1013] | 1.1404 | 0.0030 [0.0009, 0.0061] | 0.1299 [0.1202, 0.1416] | 0.8671 [0.8554, 0.8766] |
| nadig_jurkat | `all_nontarget` | 1779 | 6148 (C_dist>=5) | 0.4031 [0.3852, 0.4211] | 0.4143 | 0.2105 [0.1974, 0.2240] | 0.0386 [0.0354, 0.0419] | 0.7509 [0.7390, 0.7627] |
| nourreddine_kolf21j | `resp_top50` | 28 | 25 (C_dist>=5) | 0.7390 [0.5864, 0.9558] | 0.6499 | 0.0362 [0.0011, 0.0994] | 0.0244 [0.0000, 0.0525] | 0.9394 [0.8729, 0.9903] |
| nourreddine_kolf21j | `all_nontarget` | 28 | 6148 (C_dist>=5) | 0.1581 [0.1016, 0.2267] | 0.1434 | 0.4723 [0.3774, 0.5496] | 0.0084 [0.0063, 0.0225] | 0.5193 [0.4342, 0.6079] |
| replogle_k562_essential | `resp_top50` | 1550 | 2271 (C_dist>=5) | 0.8395 [0.8176, 0.8608] | 0.8388 | 0.0316 [0.0232, 0.0416] | 0.1500 [0.1403, 0.1632] | 0.8183 [0.8015, 0.8301] |
| replogle_k562_essential | `all_nontarget` | 1550 | 6148 (C_dist>=5) | 0.3372 [0.3235, 0.3506] | 0.3382 | 0.3897 [0.3749, 0.4053] | 0.0336 [0.0307, 0.0365] | 0.5767 [0.5628, 0.5903] |
| replogle_rpe1 | `resp_top50` | 1779 | 1564 (C_dist>=5) | 1.8719 [1.8143, 1.9279] | 2.1842 | 0.1569 [0.1422, 0.1724] | 0.0900 [0.0815, 0.0989] | 0.7531 [0.7366, 0.7697] |
| replogle_rpe1 | `all_nontarget` | 1779 | 6148 (C_dist>=5) | 0.8362 [0.8120, 0.8598] | 0.8639 | 0.0082 [0.0059, 0.0109] | 0.1154 [0.1098, 0.1210] | 0.8764 [0.8719, 0.8808] |
| vcc2025_h1 | `resp_top50` | 95 | 246 (C_dist>=5) | 0.8066 [0.5049, 1.3667] | 0.7907 | 0.0120 [0.0000, 0.1026] | 0.2325 [0.1243, 0.2577] | 0.7555 [0.6893, 0.8608] |
| vcc2025_h1 | `all_nontarget` | 95 | 6148 (C_dist>=5) | 0.3503 [0.2664, 0.4420] | 0.3529 | 0.2247 [0.1502, 0.3172] | 0.0458 [0.0329, 0.0641] | 0.7295 [0.6423, 0.7931] |

## S6: Split-Half Noise Ceiling (`resp_top50` vs Unselected `all_nontarget`, Amplitude-Orthogonalized & Raw)

| Target Line | N half-valid | rho_AB_orth (top50) | Ceil r_orth (top50) | rho_AB (top50) | Rep Slope (top50) | rho_AB_orth (all) | rho_AB (all) | TMVS r_orth_cal (all) | Hybrid r_orth_cal (all) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| nadig_hepg2 | 1777 | 0.8019 | 0.9434 | 0.8182 | 0.8182 | 0.2737 | 0.3551 | +0.1840 | +0.0353 |
| nadig_jurkat | 1779 | 0.8030 | 0.9438 | 0.8146 | 0.8146 | 0.2434 | 0.3768 | +0.1143 | +0.0411 |
| nourreddine_kolf21j | 28 | 0.7985 | 0.9423 | 0.8125 | 0.8126 | 0.1272 | 0.4294 | +0.0882 | -0.0247 |
| replogle_k562_essential | 1550 | 0.8441 | 0.9568 | 0.8583 | 0.8583 | 0.3150 | 0.5409 | +0.1232 | +0.0185 |
| replogle_rpe1 | 1779 | 0.8896 | 0.9703 | 0.9022 | 0.9022 | 0.5218 | 0.5440 | +0.1739 | +0.1625 |
| vcc2025_h1 | 95 | 0.9889 | 0.9972 | 0.9911 | 0.9911 | 0.8874 | 0.9187 | +0.0438 | -0.0021 |

## S5: Expressed Genes Sensitivity (Drop Bottom 25% Target Control Variance Genes)

| Target Line | Genes Kept | rho_AB_orth (S5) | Method | Slope OLS (S5) | r_orth_cal (S5) | r_orth_cal^2 (S5) | Joint 97.5% CI | r_cal raw (S5) |
|---|---:|---:|---|---:|---:|---:|---|---:|
| nadig_hepg2 | 4611/6148 | 0.8018 | `Hybrid_ProjShared` | +0.0062 | +0.0191 | 0.0004 | [-0.0002, +0.0387] | +0.0571 |
| nadig_hepg2 | 4611/6148 | 0.8018 | `Transfer_Mean_VarScaled` | +0.0426 | +0.0330 | 0.0011 | [+0.0191, +0.0476] | +0.1843 |
| nadig_hepg2 | 4611/6148 | 0.8018 | `Transfer_Mean_GlobalScale` | -0.0117 | +nan | nan | [nan, nan] | -0.2148 |
| nadig_hepg2 | 4611/6148 | 0.8018 | `Hybrid_on_TMVS` | +0.0488 | +0.0356 | 0.0013 | [+0.0197, +0.0516] | +0.1774 |
| nadig_hepg2 | 4611/6148 | 0.8018 | `FCT` | +0.0411 | +0.2173 | 0.0472 | [+0.2020, +0.2321] | +0.2259 |
| nadig_hepg2 | 4611/6148 | 0.8018 | `NegCtrl_HybridProj_SrcCorrTgtVar` | +0.0030 | +0.0335 | 0.0011 | [+0.0229, +0.0438] | +0.0641 |
| nadig_hepg2 | 4611/6148 | 0.8018 | `NegCtrl_TMVS_ShuffledVar` | -0.0380 | +0.0508 | 0.0026 | [+0.0413, +0.0597] | -0.1440 |
| nadig_jurkat | 4611/6148 | 0.8032 | `Hybrid_ProjShared` | +0.0040 | +0.0612 | 0.0037 | [+0.0437, +0.0792] | +0.0552 |
| nadig_jurkat | 4611/6148 | 0.8032 | `Transfer_Mean_VarScaled` | +0.0151 | +0.0885 | 0.0078 | [+0.0780, +0.0989] | +0.0849 |
| nadig_jurkat | 4611/6148 | 0.8032 | `Transfer_Mean_GlobalScale` | +0.0063 | +nan | nan | [nan, nan] | +0.0552 |
| nadig_jurkat | 4611/6148 | 0.8032 | `Hybrid_on_TMVS` | +0.0191 | +0.1036 | 0.0107 | [+0.0909, +0.1167] | +0.0977 |
| nadig_jurkat | 4611/6148 | 0.8032 | `FCT` | +0.0148 | +0.0852 | 0.0073 | [+0.0718, +0.0985] | +0.0826 |
| nadig_jurkat | 4611/6148 | 0.8032 | `NegCtrl_HybridProj_SrcCorrTgtVar` | +0.0019 | +0.0444 | 0.0020 | [+0.0313, +0.0570] | +0.0469 |
| nadig_jurkat | 4611/6148 | 0.8032 | `NegCtrl_TMVS_ShuffledVar` | +0.0285 | +0.1137 | 0.0129 | [+0.1045, +0.1227] | +0.1175 |
| nourreddine_kolf21j | 4611/6148 | 0.7987 | `Hybrid_ProjShared` | -0.0072 | +0.0143 | 0.0002 | [-0.0932, +0.1095] | -0.0450 |
| nourreddine_kolf21j | 4611/6148 | 0.7987 | `Transfer_Mean_VarScaled` | -0.0205 | +0.1000 | 0.0100 | [-0.0311, +0.2165] | -0.0740 |
| nourreddine_kolf21j | 4611/6148 | 0.7987 | `Transfer_Mean_GlobalScale` | -0.0337 | +nan | nan | [nan, nan] | -0.1997 |
| nourreddine_kolf21j | 4611/6148 | 0.7987 | `Hybrid_on_TMVS` | -0.0277 | +0.0873 | 0.0076 | [-0.0395, +0.1947] | -0.0828 |
| nourreddine_kolf21j | 4611/6148 | 0.7987 | `FCT` | +0.0301 | +0.1097 | 0.0120 | [-0.0134, +0.2239] | +0.1152 |
| nourreddine_kolf21j | 4611/6148 | 0.7987 | `NegCtrl_HybridProj_SrcCorrTgtVar` | -0.0026 | +0.0649 | 0.0042 | [-0.0785, +0.1987] | -0.0209 |
| nourreddine_kolf21j | 4611/6148 | 0.7987 | `NegCtrl_TMVS_ShuffledVar` | +0.0115 | +0.0276 | 0.0008 | [-0.0190, +0.0971] | +0.0472 |
| replogle_k562_essential | 4611/6148 | 0.8432 | `Hybrid_ProjShared` | +0.0029 | +0.0375 | 0.0014 | [+0.0209, +0.0538] | +0.0320 |
| replogle_k562_essential | 4611/6148 | 0.8432 | `Transfer_Mean_VarScaled` | +0.0420 | +0.1750 | 0.0306 | [+0.1578, +0.1923] | +0.2131 |
| replogle_k562_essential | 4611/6148 | 0.8432 | `Transfer_Mean_GlobalScale` | +0.0016 | +nan | nan | [nan, nan] | +0.1797 |
| replogle_k562_essential | 4611/6148 | 0.8432 | `Hybrid_on_TMVS` | +0.0448 | +0.1715 | 0.0294 | [+0.1541, +0.1891] | +0.2063 |
| replogle_k562_essential | 4611/6148 | 0.8432 | `FCT` | +0.0291 | +0.1492 | 0.0223 | [+0.1346, +0.1639] | +0.1221 |
| replogle_k562_essential | 4611/6148 | 0.8432 | `NegCtrl_HybridProj_SrcCorrTgtVar` | +0.0019 | +0.0299 | 0.0009 | [+0.0160, +0.0431] | +0.0352 |
| replogle_k562_essential | 4611/6148 | 0.8432 | `NegCtrl_TMVS_ShuffledVar` | +0.0650 | +0.0621 | 0.0039 | [+0.0522, +0.0722] | +0.1808 |
| replogle_rpe1 | 4611/6148 | 0.8892 | `Hybrid_ProjShared` | +0.0072 | +0.1824 | 0.0333 | [+0.1632, +0.2005] | +0.1635 |
| replogle_rpe1 | 4611/6148 | 0.8892 | `Transfer_Mean_VarScaled` | +0.0335 | +0.0984 | 0.0097 | [+0.0854, +0.1122] | +0.3121 |
| replogle_rpe1 | 4611/6148 | 0.8892 | `Transfer_Mean_GlobalScale` | +0.0076 | +nan | nan | [nan, nan] | +0.4106 |
| replogle_rpe1 | 4611/6148 | 0.8892 | `Hybrid_on_TMVS` | +0.0407 | +0.1636 | 0.0268 | [+0.1492, +0.1790] | +0.3355 |
| replogle_rpe1 | 4611/6148 | 0.8892 | `FCT` | +0.0199 | +0.1305 | 0.0170 | [+0.1174, +0.1438] | +0.1877 |
| replogle_rpe1 | 4611/6148 | 0.8892 | `NegCtrl_HybridProj_SrcCorrTgtVar` | -0.0010 | -0.0319 | 0.0010 | [-0.0477, -0.0167] | -0.0487 |
| replogle_rpe1 | 4611/6148 | 0.8892 | `NegCtrl_TMVS_ShuffledVar` | -0.0358 | +0.0272 | 0.0007 | [+0.0192, +0.0351] | -0.2370 |
| vcc2025_h1 | 4611/6148 | 0.9888 | `Hybrid_ProjShared` | -0.0145 | -0.0989 | 0.0098 | [-0.2258, +0.0079] | -0.1013 |
| vcc2025_h1 | 4611/6148 | 0.9888 | `Transfer_Mean_VarScaled` | -0.1304 | -0.0524 | 0.0027 | [-0.1601, +0.0561] | -0.4188 |
| vcc2025_h1 | 4611/6148 | 0.9888 | `Transfer_Mean_GlobalScale` | +0.0260 | +nan | nan | [nan, nan] | +0.1124 |
| vcc2025_h1 | 4611/6148 | 0.9888 | `Hybrid_on_TMVS` | -0.1448 | -0.1054 | 0.0111 | [-0.2160, +0.0035] | -0.4137 |
| vcc2025_h1 | 4611/6148 | 0.9888 | `FCT` | +0.2012 | +0.3608 | 0.1302 | [+0.2824, +0.4370] | +0.5275 |
| vcc2025_h1 | 4611/6148 | 0.9888 | `NegCtrl_HybridProj_SrcCorrTgtVar` | -0.0003 | -0.0624 | 0.0039 | [-0.2215, +0.0630] | -0.0078 |
| vcc2025_h1 | 4611/6148 | 0.9888 | `NegCtrl_TMVS_ShuffledVar` | +0.0143 | +0.0308 | 0.0009 | [-0.0075, +0.0853] | +0.0531 |

## Appendix: Alternative H-A4/H-B4 Discrimination & Centroid Conventions (R5)

| Target Line | Method | Ctrl disc_raw (A-T2) | Ctrl disc_raw_panel (A-T3) | Cent disc_raw | Cent disc_raw_panel | Cent de50_r | pds_cosine_panel | Slope OLS (All) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| nadig_hepg2 | `Source_Grand_Mean` | 0.5001 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.2247 |
| nadig_hepg2 | `Transfer_Mean` | 0.7706 | 0.7810 | 0.8774 | 0.8815 | 0.5340 | 0.7794 | +0.0000 |
| nadig_hepg2 | `CIPHER_Linear` | 0.5001 | 0.5005 | 0.4939 | 0.4974 | 0.0950 | 0.4985 | +0.2457 |
| nadig_hepg2 | `Hybrid_Raw` | 0.7449 | 0.7560 | 0.8606 | 0.8642 | 0.5143 | 0.7531 | +0.0042 |
| nadig_hepg2 | `Hybrid_ProjShared` | 0.7577 | 0.7694 | 0.8710 | 0.8747 | 0.5269 | 0.7689 | +0.0040 |
| nadig_hepg2 | `NegCtrl_HybridProj_ShuffledGenes` | 0.7681 | 0.7775 | 0.8751 | 0.8784 | 0.5180 | 0.7763 | +0.0016 |
| nadig_hepg2 | `NegCtrl_HybridProj_DiagCov` | 0.7685 | 0.7784 | 0.8755 | 0.8794 | 0.5255 | 0.7772 | +0.0015 |
| nadig_hepg2 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.7677 | 0.7778 | 0.8769 | 0.8808 | 0.5313 | 0.7765 | +0.0012 |
| nadig_hepg2 | `Transfer_Mean_VarScaled` | 0.7645 | 0.7770 | 0.8729 | 0.8788 | 0.5144 | 0.7757 | +0.0266 |
| nadig_hepg2 | `Transfer_Mean_GlobalScale` | 0.7706 | 0.7810 | 0.8795 | 0.8831 | 0.5395 | 0.7794 | +0.0188 |
| nadig_hepg2 | `NegCtrl_TMVS_ShuffledVar` | 0.7780 | 0.7890 | 0.8739 | 0.8778 | 0.5183 | 0.7877 | +0.0432 |
| nadig_hepg2 | `Hybrid_on_TMVS` | 0.7515 | 0.7654 | 0.8658 | 0.8718 | 0.5042 | 0.7651 | +0.0306 |
| nadig_hepg2 | `FCT` | 0.7770 | 0.7873 | 0.8799 | 0.8835 | 0.5475 | 0.7824 | +0.0385 |
| nadig_jurkat | `Source_Grand_Mean` | 0.5001 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.2701 |
| nadig_jurkat | `Transfer_Mean` | 0.7927 | 0.7978 | 0.8730 | 0.8766 | 0.4933 | 0.7900 | +0.0000 |
| nadig_jurkat | `CIPHER_Linear` | 0.4936 | 0.4954 | 0.4977 | 0.4986 | 0.3708 | 0.4996 | +0.3450 |
| nadig_jurkat | `Hybrid_Raw` | 0.7917 | 0.7966 | 0.8724 | 0.8762 | 0.5158 | 0.7877 | +0.0130 |
| nadig_jurkat | `Hybrid_ProjShared` | 0.7851 | 0.7905 | 0.8716 | 0.8753 | 0.4948 | 0.7837 | +0.0027 |
| nadig_jurkat | `NegCtrl_HybridProj_ShuffledGenes` | 0.7890 | 0.7940 | 0.8699 | 0.8740 | 0.4864 | 0.7864 | -0.0011 |
| nadig_jurkat | `NegCtrl_HybridProj_DiagCov` | 0.7904 | 0.7961 | 0.8702 | 0.8745 | 0.4886 | 0.7881 | -0.0012 |
| nadig_jurkat | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.7907 | 0.7959 | 0.8727 | 0.8764 | 0.4961 | 0.7887 | +0.0005 |
| nadig_jurkat | `Transfer_Mean_VarScaled` | 0.7926 | 0.7970 | 0.8705 | 0.8735 | 0.4679 | 0.7886 | -0.0427 |
| nadig_jurkat | `Transfer_Mean_GlobalScale` | 0.7927 | 0.7978 | 0.8714 | 0.8751 | 0.4610 | 0.7900 | -0.0547 |
| nadig_jurkat | `NegCtrl_TMVS_ShuffledVar` | 0.7989 | 0.8036 | 0.8762 | 0.8801 | 0.5473 | 0.7969 | +0.0038 |
| nadig_jurkat | `Hybrid_on_TMVS` | 0.7863 | 0.7910 | 0.8691 | 0.8717 | 0.4704 | 0.7836 | -0.0400 |
| nadig_jurkat | `FCT` | 0.7878 | 0.7924 | 0.8627 | 0.8670 | 0.4759 | 0.7755 | -0.0157 |
| nourreddine_kolf21j | `Source_Grand_Mean` | 0.5000 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.4161 |
| nourreddine_kolf21j | `Transfer_Mean` | 0.7130 | 0.7156 | 0.7447 | 0.7474 | 0.3766 | 0.6905 | +0.0000 |
| nourreddine_kolf21j | `CIPHER_Linear` | 0.5251 | 0.5251 | 0.5053 | 0.5053 | 0.5821 | 0.4960 | +0.5562 |
| nourreddine_kolf21j | `Hybrid_Raw` | 0.7183 | 0.7235 | 0.7579 | 0.7579 | 0.4289 | 0.6852 | +0.0238 |
| nourreddine_kolf21j | `Hybrid_ProjShared` | 0.7235 | 0.7222 | 0.7698 | 0.7672 | 0.3712 | 0.7011 | -0.0024 |
| nourreddine_kolf21j | `NegCtrl_HybridProj_ShuffledGenes` | 0.7130 | 0.7130 | 0.7579 | 0.7579 | 0.3582 | 0.6892 | +0.0009 |
| nourreddine_kolf21j | `NegCtrl_HybridProj_DiagCov` | 0.7037 | 0.7050 | 0.7619 | 0.7606 | 0.3626 | 0.6892 | +0.0012 |
| nourreddine_kolf21j | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.7116 | 0.7116 | 0.7381 | 0.7394 | 0.3863 | 0.6839 | -0.0002 |
| nourreddine_kolf21j | `Transfer_Mean_VarScaled` | 0.7130 | 0.7143 | 0.7460 | 0.7460 | 0.4081 | 0.6944 | -0.0980 |
| nourreddine_kolf21j | `Transfer_Mean_GlobalScale` | 0.7130 | 0.7156 | 0.7354 | 0.7354 | 0.2977 | 0.6905 | -0.1374 |
| nourreddine_kolf21j | `NegCtrl_TMVS_ShuffledVar` | 0.7090 | 0.7143 | 0.7500 | 0.7487 | 0.4662 | 0.6905 | -0.0335 |
| nourreddine_kolf21j | `Hybrid_on_TMVS` | 0.7302 | 0.7302 | 0.7632 | 0.7619 | 0.4013 | 0.7063 | -0.1004 |
| nourreddine_kolf21j | `FCT` | 0.7130 | 0.7130 | 0.7566 | 0.7500 | 0.4213 | 0.6799 | +0.0392 |
| replogle_k562_essential | `Source_Grand_Mean` | 0.5001 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.4245 |
| replogle_k562_essential | `Transfer_Mean` | 0.7920 | 0.8035 | 0.8951 | 0.8999 | 0.5037 | 0.7989 | +0.0000 |
| replogle_k562_essential | `CIPHER_Linear` | 0.4583 | 0.4628 | 0.4943 | 0.4948 | 0.4685 | 0.4746 | +0.5778 |
| replogle_k562_essential | `Hybrid_Raw` | 0.8017 | 0.8104 | 0.8989 | 0.9019 | 0.5402 | 0.8091 | +0.0256 |
| replogle_k562_essential | `Hybrid_ProjShared` | 0.7898 | 0.7998 | 0.8947 | 0.8989 | 0.4976 | 0.7955 | +0.0011 |
| replogle_k562_essential | `NegCtrl_HybridProj_ShuffledGenes` | 0.7948 | 0.8051 | 0.8967 | 0.9000 | 0.4851 | 0.8007 | -0.0019 |
| replogle_k562_essential | `NegCtrl_HybridProj_DiagCov` | 0.7947 | 0.8048 | 0.8966 | 0.9002 | 0.4872 | 0.8003 | -0.0020 |
| replogle_k562_essential | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.7910 | 0.8017 | 0.8952 | 0.8999 | 0.5064 | 0.7977 | +0.0002 |
| replogle_k562_essential | `Transfer_Mean_VarScaled` | 0.7928 | 0.8038 | 0.8985 | 0.9028 | 0.5585 | 0.7991 | +0.0316 |
| replogle_k562_essential | `Transfer_Mean_GlobalScale` | 0.7920 | 0.8035 | 0.8955 | 0.9003 | 0.5064 | 0.7989 | +0.0048 |
| replogle_k562_essential | `NegCtrl_TMVS_ShuffledVar` | 0.7958 | 0.8085 | 0.9040 | 0.9076 | 0.6061 | 0.8057 | +0.0887 |
| replogle_k562_essential | `Hybrid_on_TMVS` | 0.7905 | 0.8002 | 0.8983 | 0.9018 | 0.5501 | 0.7956 | +0.0327 |
| replogle_k562_essential | `FCT` | 0.7847 | 0.7987 | 0.8864 | 0.8968 | 0.5127 | 0.7849 | +0.0241 |
| replogle_rpe1 | `Source_Grand_Mean` | 0.5003 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.1076 |
| replogle_rpe1 | `Transfer_Mean` | 0.8612 | 0.8659 | 0.8770 | 0.8823 | 0.4361 | 0.8623 | +0.0000 |
| replogle_rpe1 | `CIPHER_Linear` | 0.5219 | 0.5242 | 0.5084 | 0.5115 | -0.1883 | 0.5250 | +0.0753 |
| replogle_rpe1 | `Hybrid_Raw` | 0.8493 | 0.8535 | 0.8888 | 0.8954 | 0.4728 | 0.8515 | +0.0131 |
| replogle_rpe1 | `Hybrid_ProjShared` | 0.8537 | 0.8578 | 0.8846 | 0.8919 | 0.4587 | 0.8553 | +0.0086 |
| replogle_rpe1 | `NegCtrl_HybridProj_ShuffledGenes` | 0.8614 | 0.8648 | 0.8707 | 0.8761 | 0.4165 | 0.8608 | -0.0050 |
| replogle_rpe1 | `NegCtrl_HybridProj_DiagCov` | 0.8612 | 0.8639 | 0.8707 | 0.8759 | 0.4192 | 0.8601 | -0.0049 |
| replogle_rpe1 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.8601 | 0.8640 | 0.8755 | 0.8808 | 0.4334 | 0.8601 | -0.0008 |
| replogle_rpe1 | `Transfer_Mean_VarScaled` | 0.8556 | 0.8608 | 0.8806 | 0.8862 | 0.4617 | 0.8577 | +0.0113 |
| replogle_rpe1 | `Transfer_Mean_GlobalScale` | 0.8612 | 0.8659 | 0.8795 | 0.8845 | 0.4421 | 0.8623 | -0.0022 |
| replogle_rpe1 | `NegCtrl_TMVS_ShuffledVar` | 0.8669 | 0.8712 | 0.8668 | 0.8738 | 0.3651 | 0.8672 | -0.0089 |
| replogle_rpe1 | `Hybrid_on_TMVS` | 0.8492 | 0.8539 | 0.8878 | 0.8948 | 0.4785 | 0.8517 | +0.0200 |
| replogle_rpe1 | `FCT` | 0.8625 | 0.8669 | 0.8765 | 0.8836 | 0.4447 | 0.8600 | +0.0237 |
| vcc2025_h1 | `Source_Grand_Mean` | 0.4999 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.2723 |
| vcc2025_h1 | `Transfer_Mean` | 0.8917 | 0.8907 | 0.8518 | 0.8506 | 0.5101 | 0.8758 | +0.0000 |
| vcc2025_h1 | `CIPHER_Linear` | 0.4614 | 0.4635 | 0.4834 | 0.4865 | 0.2829 | 0.4582 | +0.3362 |
| vcc2025_h1 | `Hybrid_Raw` | 0.8433 | 0.8458 | 0.8264 | 0.8263 | 0.4672 | 0.8396 | +0.0009 |
| vcc2025_h1 | `Hybrid_ProjShared` | 0.8419 | 0.8439 | 0.8356 | 0.8351 | 0.4587 | 0.8368 | -0.0018 |
| vcc2025_h1 | `NegCtrl_HybridProj_ShuffledGenes` | 0.8885 | 0.8867 | 0.8543 | 0.8538 | 0.4960 | 0.8729 | +0.0002 |
| vcc2025_h1 | `NegCtrl_HybridProj_DiagCov` | 0.8877 | 0.8871 | 0.8535 | 0.8529 | 0.5054 | 0.8740 | +0.0004 |
| vcc2025_h1 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.8832 | 0.8828 | 0.8506 | 0.8489 | 0.5011 | 0.8693 | +0.0001 |
| vcc2025_h1 | `Transfer_Mean_VarScaled` | 0.8910 | 0.8898 | 0.7382 | 0.7358 | 0.5551 | 0.8801 | +0.1390 |
| vcc2025_h1 | `Transfer_Mean_GlobalScale` | 0.8917 | 0.8907 | 0.7483 | 0.7475 | 0.5448 | 0.8758 | +0.1439 |
| vcc2025_h1 | `NegCtrl_TMVS_ShuffledVar` | 0.8952 | 0.8945 | 0.7046 | 0.7052 | 0.5407 | 0.8829 | +0.1672 |
| vcc2025_h1 | `Hybrid_on_TMVS` | 0.8152 | 0.8170 | 0.7225 | 0.7222 | 0.4842 | 0.8119 | +0.1372 |
| vcc2025_h1 | `FCT` | 0.9047 | 0.9043 | 0.8944 | 0.8935 | 0.5875 | 0.8923 | +0.0435 |
