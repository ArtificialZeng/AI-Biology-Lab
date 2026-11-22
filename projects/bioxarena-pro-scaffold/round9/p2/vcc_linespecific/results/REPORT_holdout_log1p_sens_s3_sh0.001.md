# Line-Specific Zero-Shot Perturbation Evaluation Report (`holdout_log1p_sens_s3_sh0.001`)

- **Code Freeze Commit (C1)**: `1ca9f60f49adf1ac7cce8266e40533273411db30` | **Bootstrap Seed (C2)**: `20261121`
- **Confirmatory Lines (3)**: nadig_hepg2, nadig_jurkat, vcc2025_h1
- **Secondary Lines (3)**: nourreddine_kolf21j, replogle_k562_essential, replogle_rpe1
- **Dual Confirmatory Summary (R13, 97.5% Bonferroni CIs)**: **Group A (Hybrid_ProjShared): 未确认 | Group B (Transfer_Mean_VarScaled): 未确认 (Bonferroni 97.5% CIs)**
- **Group A Primary (`Hybrid_ProjShared`)**: **未确认**
- **Group B Primary (`Transfer_Mean_VarScaled`)**: **未确认**
- **R17 Primary Correlation Metric**: `ls_resid_corr_orth_cal_top50` (amplitude-orthogonalized; raw `r_cal` reported as secondary)
- **R18 Global-Scale Scalar Baseline**: 单个标量 Transfer_Mean_GlobalScale 即可得到原始残差相关 [-0.126, +0.250]
- **R4 Co-fluctuation Structure Claim Allowed**: `目标线对照统计量` (`all beat SrcCorrTgtVar disjoint=False`)
- **R4 Marginal Variance Note**: `增益可由边际方差解释`
- **R14 Mechanism Statement (`Hybrid_on_TMVS` vs `TMVS` on `r_orth_cal`, 95% CI)**: `在本数据上，共波动结构没有超出边际方差的可检测增量` (`95%: 2/3`, `97.5%: 2/3`)
- **R20 Prior-Art Novelty Scope**: FCT is equivalent to VCC 2026 CellularAI (mu_t * 2^log2FC); CASTER performed supervised control-variance scaling. TMVS is framed strictly as a parameter-free variance-ratio variant of expression-aware transfer (not a new method); contribution is the pre-registered 3-layer decomposition of recoverable zero-shot line-specific residuals into (1) line-wide amplitude, (2) gene marginal variance, and (3) off-diagonal co-fluctuations.

## Hypothesis Adjudication Summary (Confirmatory Cohort, 97.5% CIs)

| Hypothesis | Group | Pass / Confirmatory | Threshold | Verdict |
|---|---|---:|---|---|
| **H-A1** (Slope 97.5% CI lo > 0.13) | A (`Hybrid_ProjShared`) | 0/3 | >= 3 | **未确认** |
| **H-A2** (Joint-boot `r_orth_cal` 97.5% CI lo > 0, valid>=95%, rho_AB_orth>0) | A (`Hybrid_ProjShared`) | 2/3 (fixed95: 2/3) | >= 3 | **未确认** |
| **H-A3** (Matched > ALL source mismatched cov on slope AND `r_orth_cal`) | A (`Hybrid_ProjShared`) | 0/3 | >= 3 | **未确认** |
| **H-A4** (Strict AND non-inferiority across 6 conventions) | A (`Hybrid_ProjShared`) | 3/3 inf (ctrl_disc_raw=3/3, ctrl_de50_r=2/3, ctrl_disc_raw_panel=3/3, cent_disc_raw=2/3, cent_de50_r=1/3, cent_disc_raw_panel=2/3; union=3) | <= 0 inf | **未确认** |
| **H-B1** (Slope 97.5% CI lo > 0.13) | B (`Transfer_Mean_VarScaled`) | 0/3 | >= 3 | **未确认** |
| **H-B2** (Joint-boot `r_orth_cal` 97.5% CI lo > 0, valid>=95%, rho_AB_orth>0) | B (`Transfer_Mean_VarScaled`) | 2/3 | >= 3 | **未确认** |
| **H-B3** (Matched `r_orth_cal` lo > max(ALL mismatched vars CI hi, 100-perm p97.5)) | B (`Transfer_Mean_VarScaled`) | 0/3 | >= 3 | **未确认** |
| **H-B4** (Strict AND non-inferiority across 6 conventions) | B (`Transfer_Mean_VarScaled`) | 3/3 inf (ctrl_disc_raw=2/3, ctrl_de50_r=1/3, ctrl_disc_raw_panel=2/3, cent_disc_raw=3/3, cent_de50_r=1/3, cent_disc_raw_panel=3/3; union=3) | <= 0 inf | **未确认** |

## Main Table: Primary View (A-T2 Own-Target-Zeroed Control Reference, 97.5% CIs, R5/R6/R11/R13/R17/R18)

| Cohort | Target Line | Method | N | Slope OLS (Top50) | 97.5% CI | r_orth_cal | r_orth_cal^2 | r_orth_cal^2 Joint 97.5% CI | Resid r_orth | r_cal (raw) | Resid r (raw) | sd_pred/sd_obs | Deming (delta=1) | Obs~Pred | Ctrl disc_raw | Ctrl de50_r | MSE |
|---|---|---|---:|---:|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Confirmatory | nadig_hepg2 | `Source_Grand_Mean` | 799 | +0.0194 | [-0.0105, +0.0513] | +0.3127 | 0.0978 | [0.0750, 0.1232] | +0.3020 | +0.0290 | +0.0281 | 0.6924 | +0.0373 | +0.0405 | 0.5003 | 0.6148 | 0.0072 |
| Confirmatory | nadig_hepg2 | `Transfer_Mean` | 799 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.8811 | 0.7053 | 0.0078 |
| Confirmatory | nadig_hepg2 | `CIPHER_Linear` | 799 | -0.0670 | [-0.1014, -0.0320] | +0.1174 | 0.0138 | [0.0072, 0.0226] | +0.1134 | -0.0830 | -0.0802 | 0.8348 | -0.2112 | -0.0961 | 0.4832 | 0.2594 | 0.0087 |
| Confirmatory | nadig_hepg2 | `Hybrid_Raw` | 799 | +0.0352 | [+0.0286, +0.0421] | +0.1227 | 0.0150 | [0.0079, 0.0246] | +0.1185 | +0.1832 | +0.1771 | 0.1987 | +0.0366 | +0.8915 | 0.8597 | 0.7050 | 0.0081 |
| Confirmatory | nadig_hepg2 | `Hybrid_ProjShared` | 799 | +0.0097 | [+0.0061, +0.0135] | +0.0490 | 0.0024 | [0.0006, 0.0055] | +0.0473 | +0.0720 | +0.0696 | 0.1396 | +0.0099 | +0.4988 | 0.8688 | 0.6982 | 0.0081 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_MismatchedCov(nadig_jurkat)` | 799 | -0.0012 | [-0.0023, -0.0002] | -0.0101 | 0.0001 | [0.0000, 0.0006] | -0.0097 | -0.0176 | -0.0170 | 0.0704 | -0.0012 | -0.2415 | 0.8806 | 0.7004 | 0.0079 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_MismatchedCov(nourreddine_kolf21j)` | 799 | +0.0064 | [+0.0047, +0.0082] | +0.0446 | 0.0020 | [0.0010, 0.0033] | +0.0430 | +0.0604 | +0.0584 | 0.1094 | +0.0065 | +0.5343 | 0.8788 | 0.7007 | 0.0081 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 799 | +0.0028 | [+0.0017, +0.0038] | +0.0269 | 0.0007 | [0.0002, 0.0016] | +0.0260 | +0.0398 | +0.0385 | 0.0723 | +0.0028 | +0.5322 | 0.8796 | 0.7042 | 0.0079 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 799 | -0.0046 | [-0.0056, -0.0035] | -0.0403 | 0.0016 | [0.0007, 0.0030] | -0.0390 | -0.0644 | -0.0622 | 0.0733 | -0.0046 | -0.8493 | 0.8820 | 0.6957 | 0.0080 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_MismatchedCov(vcc2025_h1)` | 799 | +0.0024 | [+0.0000, +0.0048] | +0.0202 | 0.0004 | [0.0000, 0.0011] | +0.0195 | +0.0141 | +0.0137 | 0.1787 | +0.0025 | +0.0765 | 0.8785 | 0.6842 | 0.0086 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_ShuffledGenes` | 799 | +0.0057 | [+0.0035, +0.0082] | +0.0221 | 0.0005 | [0.0001, 0.0012] | +0.0213 | +0.0411 | +0.0397 | 0.1441 | +0.0058 | +0.2757 | 0.8772 | 0.6922 | 0.0086 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_DiagCov` | 799 | +0.0065 | [+0.0045, +0.0087] | +0.0314 | 0.0010 | [0.0003, 0.0021] | +0.0303 | +0.0620 | +0.0599 | 0.1084 | +0.0066 | +0.5526 | 0.8779 | 0.7007 | 0.0080 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 799 | +0.0034 | [+0.0024, +0.0044] | +0.0405 | 0.0016 | [0.0008, 0.0029] | +0.0391 | +0.0518 | +0.0501 | 0.0681 | +0.0034 | +0.7349 | 0.8808 | 0.7056 | 0.0079 |
| Confirmatory | nadig_hepg2 | `Transfer_Mean_VarScaled` | 799 | +0.0357 | [+0.0273, +0.0441] | +0.0418 | 0.0017 | [0.0006, 0.0036] | +0.0404 | +0.1224 | +0.1183 | 0.3015 | +0.0392 | +0.3925 | 0.8745 | 0.6900 | 0.0076 |
| Confirmatory | nadig_hepg2 | `Transfer_Mean_GlobalScale` | 799 | -0.0085 | [-0.0111, -0.0057] | +nan | nan | [nan, nan] | +nan | -0.1306 | -0.1263 | 0.0670 | -0.0085 | -1.8867 | 0.8811 | 0.7053 | 0.0075 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_MismatchedVar(nadig_jurkat)` | 799 | +0.0076 | [+0.0038, +0.0112] | -0.0619 | 0.0038 | [0.0025, 0.0055] | -0.0598 | +0.0503 | +0.0487 | 0.1559 | +0.0078 | +0.3122 | 0.8803 | 0.6929 | 0.0087 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_MismatchedVar(nourreddine_kolf21j)` | 799 | +0.0088 | [+0.0049, +0.0121] | +0.0472 | 0.0022 | [0.0011, 0.0037] | +0.0456 | +0.0495 | +0.0479 | 0.1835 | +0.0091 | +0.2608 | 0.8859 | 0.6991 | 0.0088 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 799 | +0.0085 | [+0.0043, +0.0136] | +0.0951 | 0.0091 | [0.0059, 0.0130] | +0.0919 | +0.0599 | +0.0579 | 0.1465 | +0.0087 | +0.3953 | 0.8831 | 0.7043 | 0.0078 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 799 | +0.0150 | [+0.0071, +0.0227] | -0.0401 | 0.0016 | [0.0003, 0.0039] | -0.0388 | +0.0694 | +0.0671 | 0.2231 | +0.0158 | +0.3010 | 0.8731 | 0.6888 | 0.0083 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_MismatchedVar(vcc2025_h1)` | 799 | -0.0579 | [-0.0725, -0.0429] | -0.0181 | 0.0003 | [0.0000, 0.0018] | -0.0175 | -0.1403 | -0.1356 | 0.4270 | -0.0705 | -0.3176 | 0.8823 | 0.6729 | 0.0070 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_ShuffledVar` | 799 | -0.0196 | [-0.0289, -0.0099] | +0.0589 | 0.0035 | [0.0023, 0.0049] | +0.0569 | -0.0621 | -0.0601 | 0.3266 | -0.0220 | -0.1839 | 0.8902 | 0.6707 | 0.0077 |
| Confirmatory | nadig_hepg2 | `Hybrid_on_TMVS` | 799 | +0.0454 | [+0.0356, +0.0556] | +0.0549 | 0.0030 | [0.0012, 0.0058] | +0.0530 | +0.1288 | +0.1246 | 0.3644 | +0.0522 | +0.3420 | 0.8611 | 0.6800 | 0.0079 |
| Confirmatory | nadig_hepg2 | `FCT` | 799 | +0.0588 | [+0.0536, +0.0640] | +0.2611 | 0.0682 | [0.0585, 0.0788] | +0.2522 | +0.2618 | +0.2532 | 0.2322 | +0.0619 | +1.0905 | 0.8871 | 0.7269 | 0.0076 |
| Confirmatory | nadig_jurkat | `Source_Grand_Mean` | 934 | -0.0036 | [-0.0349, +0.0310] | +0.0153 | 0.0002 | [0.0000, 0.0020] | +0.0147 | -0.0049 | -0.0047 | 0.7616 | -0.0086 | -0.0062 | 0.5002 | 0.3335 | 0.0068 |
| Confirmatory | nadig_jurkat | `Transfer_Mean` | 934 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.8717 | 0.6172 | 0.0076 |
| Confirmatory | nadig_jurkat | `CIPHER_Linear` | 934 | +0.0067 | [-0.0288, +0.0460] | +0.0440 | 0.0019 | [0.0001, 0.0060] | +0.0421 | +0.0076 | +0.0073 | 0.9149 | +0.0408 | +0.0080 | 0.4832 | 0.0743 | 0.0070 |
| Confirmatory | nadig_jurkat | `Hybrid_Raw` | 934 | +0.0200 | [+0.0141, +0.0267] | +0.1290 | 0.0166 | [0.0097, 0.0258] | +0.1235 | +0.1478 | +0.1419 | 0.1408 | +0.0204 | +1.0082 | 0.8743 | 0.6225 | 0.0076 |
| Confirmatory | nadig_jurkat | `Hybrid_ProjShared` | 934 | +0.0058 | [+0.0028, +0.0089] | +0.0652 | 0.0043 | [0.0017, 0.0081] | +0.0624 | +0.0571 | +0.0549 | 0.1059 | +0.0059 | +0.5182 | 0.8685 | 0.6159 | 0.0078 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_MismatchedCov(nadig_hepg2)` | 934 | -0.0054 | [-0.0073, -0.0034] | -0.0484 | 0.0023 | [0.0008, 0.0047] | -0.0464 | -0.0569 | -0.0547 | 0.0987 | -0.0054 | -0.5537 | 0.8677 | 0.6050 | 0.0078 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_MismatchedCov(nourreddine_kolf21j)` | 934 | -0.0001 | [-0.0024, +0.0022] | +0.0055 | 0.0000 | [0.0000, 0.0007] | +0.0053 | -0.0006 | -0.0005 | 0.1129 | -0.0001 | -0.0048 | 0.8726 | 0.6078 | 0.0079 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 934 | +0.0006 | [-0.0010, +0.0023] | +0.0129 | 0.0002 | [0.0000, 0.0009] | +0.0123 | +0.0077 | +0.0074 | 0.0826 | +0.0006 | +0.0892 | 0.8720 | 0.6132 | 0.0078 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 934 | +0.0024 | [+0.0010, +0.0038] | +0.0288 | 0.0008 | [0.0001, 0.0023] | +0.0276 | +0.0350 | +0.0336 | 0.0712 | +0.0024 | +0.4722 | 0.8708 | 0.6151 | 0.0077 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_MismatchedCov(vcc2025_h1)` | 934 | +0.0076 | [+0.0028, +0.0124] | +0.0354 | 0.0013 | [0.0002, 0.0032] | +0.0339 | +0.0399 | +0.0383 | 0.1983 | +0.0079 | +0.1930 | 0.8702 | 0.5983 | 0.0084 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_ShuffledGenes` | 934 | +0.0025 | [-0.0002, +0.0052] | +0.0159 | 0.0003 | [0.0000, 0.0013] | +0.0152 | +0.0205 | +0.0197 | 0.1262 | +0.0025 | +0.1561 | 0.8706 | 0.6079 | 0.0081 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_DiagCov` | 934 | +0.0022 | [-0.0004, +0.0048] | +0.0171 | 0.0003 | [0.0000, 0.0017] | +0.0163 | +0.0211 | +0.0202 | 0.1072 | +0.0022 | +0.1887 | 0.8712 | 0.6112 | 0.0078 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_SrcCorrTgtVar` | 934 | +0.0015 | [+0.0004, +0.0027] | +0.0283 | 0.0008 | [0.0001, 0.0020] | +0.0271 | +0.0242 | +0.0232 | 0.0662 | +0.0015 | +0.3508 | 0.8709 | 0.6150 | 0.0077 |
| Confirmatory | nadig_jurkat | `Transfer_Mean_VarScaled` | 934 | +0.0095 | [+0.0028, +0.0157] | +0.0947 | 0.0090 | [0.0064, 0.0119] | +0.0906 | +0.0439 | +0.0422 | 0.2245 | +0.0100 | +0.1878 | 0.8709 | 0.6187 | 0.0087 |
| Confirmatory | nadig_jurkat | `Transfer_Mean_GlobalScale` | 934 | -0.0008 | [-0.0068, +0.0047] | +nan | nan | [nan, nan] | +nan | -0.0058 | -0.0056 | 0.1426 | -0.0008 | -0.0393 | 0.8717 | 0.6172 | 0.0088 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_MismatchedVar(nadig_hepg2)` | 934 | -0.0370 | [-0.0438, -0.0300] | -0.1779 | 0.0317 | [0.0259, 0.0383] | -0.1704 | -0.1784 | -0.1714 | 0.2160 | -0.0388 | -0.7933 | 0.8671 | 0.5856 | 0.0079 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_MismatchedVar(nourreddine_kolf21j)` | 934 | +0.0417 | [+0.0332, +0.0488] | +0.1588 | 0.0252 | [0.0208, 0.0300] | +0.1521 | +0.1724 | +0.1656 | 0.2516 | +0.0444 | +0.6581 | 0.8797 | 0.6284 | 0.0088 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 934 | +0.0141 | [+0.0076, +0.0211] | +0.1062 | 0.0113 | [0.0077, 0.0153] | +0.1017 | +0.0893 | +0.0858 | 0.1643 | +0.0145 | +0.5221 | 0.8745 | 0.6215 | 0.0076 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 934 | -0.0309 | [-0.0381, -0.0234] | -0.1368 | 0.0187 | [0.0136, 0.0252] | -0.1310 | -0.1676 | -0.1609 | 0.1918 | -0.0320 | -0.8389 | 0.8629 | 0.5928 | 0.0086 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_MismatchedVar(vcc2025_h1)` | 934 | +0.0187 | [+0.0056, +0.0318] | +0.0365 | 0.0013 | [0.0002, 0.0039] | +0.0349 | +0.0423 | +0.0406 | 0.4603 | +0.0237 | +0.0883 | 0.8728 | 0.5810 | 0.0062 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_ShuffledVar` | 934 | +0.0440 | [+0.0379, +0.0503] | +0.1291 | 0.0167 | [0.0136, 0.0202] | +0.1236 | +0.1496 | +0.1437 | 0.3060 | +0.0484 | +0.4696 | 0.8813 | 0.6108 | 0.0082 |
| Confirmatory | nadig_jurkat | `Hybrid_on_TMVS` | 934 | +0.0153 | [+0.0084, +0.0218] | +0.1111 | 0.0124 | [0.0085, 0.0168] | +0.1064 | +0.0631 | +0.0606 | 0.2520 | +0.0163 | +0.2406 | 0.8683 | 0.6179 | 0.0088 |
| Confirmatory | nadig_jurkat | `FCT` | 934 | +0.0200 | [+0.0149, +0.0253] | +0.1043 | 0.0109 | [0.0072, 0.0153] | +0.0999 | +0.0915 | +0.0879 | 0.2274 | +0.0211 | +0.3865 | 0.8603 | 0.6115 | 0.0085 |
| Secondary | nourreddine_kolf21j | `Source_Grand_Mean` | 18 | +0.2723 | [+0.0345, +0.5563] | +0.0672 | 0.0045 | [0.0000, 0.0395] | +0.0638 | +0.3484 | +0.3333 | 0.8168 | +0.5608 | +0.4081 | 0.5000 | 0.2906 | 0.0034 |
| Secondary | nourreddine_kolf21j | `Transfer_Mean` | 18 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.8137 | 0.4335 | 0.0042 |
| Secondary | nourreddine_kolf21j | `CIPHER_Linear` | 18 | +0.2738 | [+0.0234, +0.5581] | -0.1096 | 0.0120 | [0.0000, 0.0715] | -0.1040 | +0.3147 | +0.3010 | 0.9096 | +0.7333 | +0.3309 | 0.5817 | -0.1068 | 0.0023 |
| Secondary | nourreddine_kolf21j | `Hybrid_Raw` | 18 | -0.0129 | [-0.0552, +0.0227] | +0.0058 | 0.0000 | [0.0000, 0.0201] | +0.0055 | -0.0507 | -0.0485 | 0.2666 | -0.0139 | -0.1818 | 0.7810 | 0.4346 | 0.0045 |
| Secondary | nourreddine_kolf21j | `Hybrid_ProjShared` | 18 | -0.0139 | [-0.0573, +0.0250] | +0.0082 | 0.0001 | [0.0000, 0.0224] | +0.0078 | -0.0610 | -0.0584 | 0.2373 | -0.0147 | -0.2460 | 0.8007 | 0.4300 | 0.0046 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_MismatchedCov(nadig_hepg2)` | 18 | -0.0011 | [-0.0499, +0.0434] | +0.0821 | 0.0067 | [0.0000, 0.0404] | +0.0779 | -0.0047 | -0.0045 | 0.2404 | -0.0012 | -0.0189 | 0.8137 | 0.4504 | 0.0046 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_MismatchedCov(nadig_jurkat)` | 18 | -0.0089 | [-0.0592, +0.0332] | +0.0405 | 0.0016 | [0.0000, 0.0279] | +0.0385 | -0.0405 | -0.0388 | 0.2294 | -0.0094 | -0.1690 | 0.7810 | 0.4393 | 0.0046 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 18 | -0.0103 | [-0.0665, +0.0383] | +0.0448 | 0.0020 | [0.0000, 0.0320] | +0.0426 | -0.0441 | -0.0422 | 0.2441 | -0.0109 | -0.1729 | 0.7974 | 0.4455 | 0.0045 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 18 | -0.0172 | [-0.0742, +0.0222] | +0.0036 | 0.0000 | [0.0000, 0.0128] | +0.0034 | -0.0761 | -0.0728 | 0.2367 | -0.0182 | -0.3074 | 0.7908 | 0.4300 | 0.0046 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_MismatchedCov(vcc2025_h1)` | 18 | +0.0139 | [-0.0185, +0.0612] | +0.0771 | 0.0059 | [0.0000, 0.0563] | +0.0732 | +0.0684 | +0.0654 | 0.2122 | +0.0145 | +0.3082 | 0.8170 | 0.4335 | 0.0045 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_ShuffledGenes` | 18 | -0.0099 | [-0.0615, +0.0391] | +0.0527 | 0.0028 | [0.0000, 0.0354] | +0.0500 | -0.0389 | -0.0372 | 0.2674 | -0.0107 | -0.1391 | 0.7614 | 0.4404 | 0.0047 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_DiagCov` | 18 | -0.0091 | [-0.0621, +0.0388] | +0.0577 | 0.0033 | [0.0000, 0.0402] | +0.0547 | -0.0381 | -0.0365 | 0.2492 | -0.0097 | -0.1464 | 0.7647 | 0.4448 | 0.0045 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_SrcCorrTgtVar` | 18 | -0.0021 | [-0.0505, +0.0402] | +0.0987 | 0.0097 | [0.0000, 0.0653] | +0.0937 | -0.0113 | -0.0108 | 0.1965 | -0.0022 | -0.0552 | 0.8007 | 0.4587 | 0.0044 |
| Secondary | nourreddine_kolf21j | `Transfer_Mean_VarScaled` | 18 | -0.0133 | [-0.0926, +0.0531] | +0.1385 | 0.0192 | [0.0000, 0.0824] | +0.1314 | -0.0392 | -0.0375 | 0.3552 | -0.0152 | -0.1056 | 0.8105 | 0.4380 | 0.0055 |
| Secondary | nourreddine_kolf21j | `Transfer_Mean_GlobalScale` | 18 | -0.0686 | [-0.1382, -0.0098] | +nan | nan | [nan, nan] | +nan | -0.3304 | -0.3161 | 0.2170 | -0.0716 | -1.4568 | 0.8137 | 0.4335 | 0.0057 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_MismatchedVar(nadig_hepg2)` | 18 | +0.0172 | [-0.0138, +0.0686] | +0.0174 | 0.0003 | [0.0000, 0.0273] | +0.0165 | +0.1186 | +0.1135 | 0.1513 | +0.0176 | +0.7498 | 0.8203 | 0.4238 | 0.0042 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_MismatchedVar(nadig_jurkat)` | 18 | -0.0405 | [-0.0734, -0.0031] | +0.0487 | 0.0024 | [0.0000, 0.0124] | +0.0462 | -0.2243 | -0.2146 | 0.1889 | -0.0419 | -1.1361 | 0.8039 | 0.4343 | 0.0052 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 18 | +0.0090 | [-0.0028, +0.0225] | +0.1200 | 0.0144 | [0.0029, 0.0310] | +0.1139 | +0.0744 | +0.0712 | 0.1269 | +0.0092 | +0.5610 | 0.8105 | 0.4385 | 0.0043 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 18 | +0.0052 | [-0.0188, +0.0311] | +0.0278 | 0.0008 | [0.0000, 0.0206] | +0.0264 | +0.0334 | +0.0319 | 0.1623 | +0.0053 | +0.1969 | 0.8137 | 0.4309 | 0.0048 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_MismatchedVar(vcc2025_h1)` | 18 | +0.0488 | [+0.0029, +0.0906] | -0.1371 | 0.0188 | [0.0000, 0.0728] | -0.1301 | +0.1333 | +0.1275 | 0.3827 | +0.0570 | +0.3331 | 0.8007 | 0.3961 | 0.0027 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_ShuffledVar` | 18 | +0.0332 | [-0.0146, +0.0883] | +0.0389 | 0.0015 | [0.0000, 0.0115] | +0.0370 | +0.1165 | +0.1114 | 0.2975 | +0.0363 | +0.3745 | 0.8105 | 0.4234 | 0.0049 |
| Secondary | nourreddine_kolf21j | `Hybrid_on_TMVS` | 18 | -0.0272 | [-0.1282, +0.0572] | +0.1149 | 0.0132 | [0.0000, 0.0734] | +0.1091 | -0.0636 | -0.0608 | 0.4468 | -0.0339 | -0.1361 | 0.8007 | 0.4411 | 0.0059 |
| Secondary | nourreddine_kolf21j | `FCT` | 18 | +0.0747 | [+0.0070, +0.1592] | +0.1839 | 0.0338 | [0.0006, 0.1256] | +0.1746 | +0.2547 | +0.2436 | 0.3066 | +0.0819 | +0.7947 | 0.8170 | 0.4414 | 0.0044 |
| Secondary | replogle_k562_essential | `Source_Grand_Mean` | 859 | +0.2588 | [+0.2086, +0.3093] | -0.0643 | 0.0041 | [0.0017, 0.0076] | -0.0623 | +0.2554 | +0.2487 | 1.0408 | +1.1736 | +0.2390 | 0.5001 | 0.3941 | 0.0052 |
| Secondary | replogle_k562_essential | `Transfer_Mean` | 859 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.8868 | 0.6788 | 0.0073 |
| Secondary | replogle_k562_essential | `CIPHER_Linear` | 859 | +0.3499 | [+0.2931, +0.4075] | +0.0256 | 0.0007 | [0.0000, 0.0027] | +0.0248 | +0.2843 | +0.2768 | 1.2641 | +2.1697 | +0.2190 | 0.4382 | 0.0712 | 0.0049 |
| Secondary | replogle_k562_essential | `Hybrid_Raw` | 859 | +0.0312 | [+0.0254, +0.0373] | +0.1045 | 0.0109 | [0.0056, 0.0172] | +0.1013 | +0.1666 | +0.1622 | 0.1925 | +0.0324 | +0.8423 | 0.8931 | 0.6811 | 0.0069 |
| Secondary | replogle_k562_essential | `Hybrid_ProjShared` | 859 | +0.0018 | [-0.0013, +0.0047] | +0.0330 | 0.0011 | [0.0001, 0.0033] | +0.0319 | +0.0140 | +0.0136 | 0.1295 | +0.0018 | +0.1053 | 0.8840 | 0.6763 | 0.0075 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(nadig_hepg2)` | 859 | -0.0095 | [-0.0121, -0.0072] | -0.0929 | 0.0086 | [0.0060, 0.0117] | -0.0901 | -0.0818 | -0.0796 | 0.1191 | -0.0096 | -0.6682 | 0.8854 | 0.6666 | 0.0075 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(nadig_jurkat)` | 859 | +0.0002 | [-0.0021, +0.0023] | +0.0018 | 0.0000 | [0.0000, 0.0006] | +0.0017 | +0.0022 | +0.0022 | 0.1049 | +0.0002 | +0.0207 | 0.8862 | 0.6749 | 0.0074 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(nourreddine_kolf21j)` | 859 | -0.0057 | [-0.0087, -0.0030] | -0.0117 | 0.0001 | [0.0000, 0.0009] | -0.0114 | -0.0381 | -0.0371 | 0.1542 | -0.0059 | -0.2407 | 0.8924 | 0.6671 | 0.0076 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 859 | +0.0049 | [+0.0030, +0.0066] | +0.0735 | 0.0054 | [0.0022, 0.0096] | +0.0713 | +0.0566 | +0.0551 | 0.0882 | +0.0049 | +0.6250 | 0.8865 | 0.6789 | 0.0074 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(vcc2025_h1)` | 859 | +0.0061 | [+0.0017, +0.0106] | +0.0242 | 0.0006 | [0.0000, 0.0018] | +0.0235 | +0.0228 | +0.0222 | 0.2742 | +0.0066 | +0.0809 | 0.8829 | 0.6570 | 0.0081 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_ShuffledGenes` | 859 | -0.0084 | [-0.0118, -0.0052] | -0.0142 | 0.0002 | [0.0000, 0.0010] | -0.0138 | -0.0521 | -0.0508 | 0.1660 | -0.0087 | -0.3059 | 0.8930 | 0.6671 | 0.0078 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_DiagCov` | 859 | -0.0083 | [-0.0116, -0.0053] | -0.0155 | 0.0002 | [0.0000, 0.0012] | -0.0150 | -0.0571 | -0.0556 | 0.1501 | -0.0085 | -0.3704 | 0.8931 | 0.6691 | 0.0076 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_SrcCorrTgtVar` | 859 | +0.0016 | [-0.0001, +0.0031] | +0.0219 | 0.0005 | [0.0000, 0.0017] | +0.0212 | +0.0188 | +0.0183 | 0.0865 | +0.0016 | +0.2121 | 0.8870 | 0.6762 | 0.0074 |
| Secondary | replogle_k562_essential | `Transfer_Mean_VarScaled` | 859 | +0.0605 | [+0.0510, +0.0708] | +0.1822 | 0.0332 | [0.0258, 0.0418] | +0.1766 | +0.2450 | +0.2385 | 0.2536 | +0.0644 | +0.9406 | 0.8885 | 0.6919 | 0.0069 |
| Secondary | replogle_k562_essential | `Transfer_Mean_GlobalScale` | 859 | +0.0029 | [+0.0024, +0.0034] | +nan | nan | [nan, nan] | +nan | +0.2849 | +0.2774 | 0.0104 | +0.0029 | +26.6524 | 0.8868 | 0.6788 | 0.0072 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(nadig_hepg2)` | 859 | -0.0443 | [-0.0559, -0.0318] | -0.1781 | 0.0317 | [0.0233, 0.0407] | -0.1726 | -0.1728 | -0.1682 | 0.2633 | -0.0475 | -0.6387 | 0.8838 | 0.6487 | 0.0073 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(nadig_jurkat)` | 859 | -0.0148 | [-0.0209, -0.0083] | +0.1090 | 0.0119 | [0.0089, 0.0153] | +0.1056 | -0.0674 | -0.0656 | 0.2252 | -0.0156 | -0.2912 | 0.8873 | 0.6829 | 0.0087 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(nourreddine_kolf21j)` | 859 | +0.0193 | [+0.0089, +0.0291] | +0.1581 | 0.0250 | [0.0191, 0.0315] | +0.1533 | +0.0652 | +0.0635 | 0.3042 | +0.0213 | +0.2087 | 0.8903 | 0.6885 | 0.0087 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 859 | -0.0273 | [-0.0393, -0.0139] | -0.0816 | 0.0067 | [0.0029, 0.0112] | -0.0791 | -0.1225 | -0.1193 | 0.2291 | -0.0288 | -0.5206 | 0.8814 | 0.6614 | 0.0083 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(vcc2025_h1)` | 859 | +0.1440 | [+0.1173, +0.1729] | -0.0118 | 0.0001 | [0.0000, 0.0012] | -0.0115 | +0.2279 | +0.2219 | 0.6492 | +0.2352 | +0.3418 | 0.8908 | 0.6370 | 0.0044 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_ShuffledVar` | 859 | +0.0994 | [+0.0857, +0.1132] | +0.0614 | 0.0038 | [0.0023, 0.0056] | +0.0595 | +0.2382 | +0.2319 | 0.4286 | +0.1200 | +0.5412 | 0.8921 | 0.6575 | 0.0064 |
| Secondary | replogle_k562_essential | `Hybrid_on_TMVS` | 859 | +0.0623 | [+0.0522, +0.0736] | +0.1725 | 0.0298 | [0.0223, 0.0384] | +0.1672 | +0.2236 | +0.2177 | 0.2860 | +0.0675 | +0.7610 | 0.8856 | 0.6886 | 0.0071 |
| Secondary | replogle_k562_essential | `FCT` | 859 | +0.0419 | [+0.0355, +0.0489] | +0.1938 | 0.0376 | [0.0306, 0.0462] | +0.1879 | +0.1468 | +0.1429 | 0.2933 | +0.0458 | +0.4871 | 0.8763 | 0.6883 | 0.0076 |
| Secondary | replogle_rpe1 | `Source_Grand_Mean` | 1070 | -0.0729 | [-0.0905, -0.0547] | +0.5224 | 0.2729 | [0.2516, 0.2938] | +0.5104 | -0.1541 | -0.1508 | 0.4837 | -0.0944 | -0.3117 | 0.5005 | 0.6129 | 0.0106 |
| Secondary | replogle_rpe1 | `Transfer_Mean` | 1070 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.9301 | 0.6109 | 0.0105 |
| Secondary | replogle_rpe1 | `CIPHER_Linear` | 1070 | -0.1486 | [-0.1676, -0.1291] | +0.3427 | 0.1174 | [0.0962, 0.1405] | +0.3348 | -0.2751 | -0.2691 | 0.5523 | -0.2049 | -0.4873 | 0.5120 | 0.3076 | 0.0123 |
| Secondary | replogle_rpe1 | `Hybrid_Raw` | 1070 | +0.0240 | [+0.0204, +0.0276] | +0.2602 | 0.0677 | [0.0499, 0.0878] | +0.2542 | +0.2350 | +0.2299 | 0.1042 | +0.0242 | +2.2057 | 0.9200 | 0.6376 | 0.0102 |
| Secondary | replogle_rpe1 | `Hybrid_ProjShared` | 1070 | +0.0132 | [+0.0113, +0.0152] | +0.1964 | 0.0386 | [0.0292, 0.0491] | +0.1919 | +0.1714 | +0.1677 | 0.0790 | +0.0133 | +2.1227 | 0.9248 | 0.6249 | 0.0103 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(nadig_hepg2)` | 1070 | +0.0092 | [+0.0078, +0.0107] | +0.1746 | 0.0305 | [0.0229, 0.0394] | +0.1706 | +0.1544 | +0.1511 | 0.0609 | +0.0092 | +2.4805 | 0.9288 | 0.6208 | 0.0104 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(nadig_jurkat)` | 1070 | -0.0035 | [-0.0045, -0.0026] | -0.0491 | 0.0024 | [0.0007, 0.0049] | -0.0480 | -0.0817 | -0.0800 | 0.0443 | -0.0035 | -1.8053 | 0.9296 | 0.6037 | 0.0106 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(nourreddine_kolf21j)` | 1070 | -0.0074 | [-0.0087, -0.0061] | -0.1179 | 0.0139 | [0.0100, 0.0184] | -0.1152 | -0.1134 | -0.1110 | 0.0665 | -0.0074 | -1.6682 | 0.9278 | 0.5924 | 0.0108 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 1070 | -0.0042 | [-0.0053, -0.0032] | -0.0924 | 0.0085 | [0.0050, 0.0128] | -0.0903 | -0.0958 | -0.0937 | 0.0452 | -0.0042 | -2.0715 | 0.9293 | 0.6015 | 0.0106 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(vcc2025_h1)` | 1070 | -0.0109 | [-0.0125, -0.0092] | -0.1096 | 0.0120 | [0.0087, 0.0157] | -0.1071 | -0.1156 | -0.1131 | 0.0962 | -0.0110 | -1.1762 | 0.9290 | 0.5798 | 0.0109 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_ShuffledGenes` | 1070 | -0.0087 | [-0.0103, -0.0073] | -0.1342 | 0.0180 | [0.0137, 0.0230] | -0.1312 | -0.1172 | -0.1147 | 0.0761 | -0.0088 | -1.5060 | 0.9283 | 0.5856 | 0.0109 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_DiagCov` | 1070 | -0.0074 | [-0.0089, -0.0060] | -0.1367 | 0.0187 | [0.0140, 0.0241] | -0.1336 | -0.1122 | -0.1098 | 0.0672 | -0.0074 | -1.6327 | 0.9277 | 0.5904 | 0.0107 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 1070 | -0.0014 | [-0.0021, -0.0006] | -0.0002 | 0.0000 | [0.0000, 0.0006] | -0.0002 | -0.0350 | -0.0343 | 0.0400 | -0.0014 | -0.8561 | 0.9293 | 0.6085 | 0.0105 |
| Secondary | replogle_rpe1 | `Transfer_Mean_VarScaled` | 1070 | +0.0323 | [+0.0284, +0.0360] | +0.0920 | 0.0085 | [0.0053, 0.0123] | +0.0899 | +0.2294 | +0.2245 | 0.1439 | +0.0329 | +1.5600 | 0.9254 | 0.6139 | 0.0103 |
| Secondary | replogle_rpe1 | `Transfer_Mean_GlobalScale` | 1070 | +0.0084 | [+0.0077, +0.0092] | +nan | nan | [nan, nan] | +nan | +0.3598 | +0.3521 | 0.0240 | +0.0084 | +14.6828 | 0.9301 | 0.6109 | 0.0106 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(nadig_hepg2)` | 1070 | +0.0255 | [+0.0224, +0.0283] | +0.0264 | 0.0007 | [0.0001, 0.0018] | +0.0258 | +0.1851 | +0.1811 | 0.1408 | +0.0260 | +1.2858 | 0.9267 | 0.6101 | 0.0103 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(nadig_jurkat)` | 1070 | +0.0205 | [+0.0180, +0.0229] | -0.0451 | 0.0020 | [0.0011, 0.0033] | -0.0440 | +0.1930 | +0.1888 | 0.1083 | +0.0207 | +1.7433 | 0.9293 | 0.5974 | 0.0108 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(nourreddine_kolf21j)` | 1070 | +0.0103 | [+0.0082, +0.0125] | +0.0086 | 0.0001 | [0.0000, 0.0005] | +0.0084 | +0.0836 | +0.0818 | 0.1259 | +0.0105 | +0.6497 | 0.9318 | 0.5937 | 0.0110 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 1070 | -0.0124 | [-0.0145, -0.0104] | -0.0273 | 0.0007 | [0.0002, 0.0015] | -0.0267 | -0.1345 | -0.1316 | 0.0943 | -0.0125 | -1.3960 | 0.9295 | 0.5937 | 0.0106 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(vcc2025_h1)` | 1070 | -0.0823 | [-0.0902, -0.0745] | +0.0222 | 0.0005 | [0.0000, 0.0014] | +0.0217 | -0.3083 | -0.3017 | 0.2727 | -0.0882 | -1.1063 | 0.9313 | 0.5925 | 0.0109 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_ShuffledVar` | 1070 | -0.0349 | [-0.0384, -0.0314] | +0.0259 | 0.0007 | [0.0002, 0.0013] | +0.0253 | -0.1885 | -0.1844 | 0.1894 | -0.0362 | -0.9734 | 0.9344 | 0.5723 | 0.0112 |
| Secondary | replogle_rpe1 | `Hybrid_on_TMVS` | 1070 | +0.0455 | [+0.0410, +0.0498] | +0.1760 | 0.0310 | [0.0240, 0.0385] | +0.1720 | +0.2683 | +0.2625 | 0.1735 | +0.0469 | +1.5132 | 0.9208 | 0.6236 | 0.0101 |
| Secondary | replogle_rpe1 | `FCT` | 1070 | +0.0228 | [+0.0202, +0.0253] | +0.1296 | 0.0168 | [0.0124, 0.0217] | +0.1267 | +0.1592 | +0.1558 | 0.1463 | +0.0233 | +1.0647 | 0.9308 | 0.6151 | 0.0103 |
| Confirmatory | vcc2025_h1 | `Source_Grand_Mean` | 80 | +0.1897 | [-0.0224, +0.3795] | +0.1087 | 0.0118 | [0.0022, 0.0303] | +0.1084 | +0.2991 | +0.2985 | 0.6354 | +0.2911 | +0.4698 | 0.5003 | 0.2596 | 0.0048 |
| Confirmatory | vcc2025_h1 | `Transfer_Mean` | 80 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.8932 | 0.4900 | 0.0054 |
| Confirmatory | vcc2025_h1 | `CIPHER_Linear` | 80 | +0.1589 | [-0.0636, +0.3568] | -0.0468 | 0.0022 | [0.0000, 0.0189] | -0.0467 | +0.2317 | +0.2312 | 0.6870 | +0.2777 | +0.3365 | 0.4725 | 0.0066 | 0.0045 |
| Confirmatory | vcc2025_h1 | `Hybrid_Raw` | 80 | -0.0080 | [-0.0240, +0.0053] | -0.0384 | 0.0015 | [0.0000, 0.0144] | -0.0383 | -0.0429 | -0.0428 | 0.1878 | -0.0083 | -0.2280 | 0.8579 | 0.4384 | 0.0063 |
| Confirmatory | vcc2025_h1 | `Hybrid_ProjShared` | 80 | -0.0088 | [-0.0245, +0.0045] | -0.0398 | 0.0016 | [0.0000, 0.0149] | -0.0397 | -0.0494 | -0.0493 | 0.1792 | -0.0091 | -0.2752 | 0.8513 | 0.4420 | 0.0062 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_MismatchedCov(nadig_hepg2)` | 80 | -0.0022 | [-0.0070, +0.0021] | -0.0431 | 0.0019 | [0.0000, 0.0077] | -0.0429 | -0.0255 | -0.0255 | 0.0870 | -0.0022 | -0.2927 | 0.8802 | 0.4691 | 0.0057 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_MismatchedCov(nadig_jurkat)` | 80 | +0.0003 | [-0.0032, +0.0034] | +0.0001 | 0.0000 | [0.0000, 0.0038] | +0.0001 | +0.0044 | +0.0044 | 0.0638 | +0.0003 | +0.0683 | 0.8759 | 0.4757 | 0.0056 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_MismatchedCov(nourreddine_kolf21j)` | 80 | -0.0020 | [-0.0076, +0.0035] | -0.0005 | 0.0000 | [0.0000, 0.0025] | -0.0005 | -0.0262 | -0.0261 | 0.0774 | -0.0020 | -0.3370 | 0.8837 | 0.4778 | 0.0056 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 80 | +0.0006 | [-0.0024, +0.0046] | +0.0293 | 0.0009 | [0.0000, 0.0090] | +0.0292 | +0.0093 | +0.0093 | 0.0604 | +0.0006 | +0.1544 | 0.8807 | 0.4816 | 0.0056 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 80 | -0.0002 | [-0.0042, +0.0046] | +0.0012 | 0.0000 | [0.0000, 0.0059] | +0.0012 | -0.0025 | -0.0025 | 0.0607 | -0.0002 | -0.0415 | 0.8896 | 0.4836 | 0.0055 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_ShuffledGenes` | 80 | -0.0022 | [-0.0084, +0.0046] | +0.0280 | 0.0008 | [0.0000, 0.0056] | +0.0279 | -0.0186 | -0.0186 | 0.1192 | -0.0023 | -0.1561 | 0.8828 | 0.4699 | 0.0063 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_DiagCov` | 80 | -0.0020 | [-0.0079, +0.0045] | +0.0400 | 0.0016 | [0.0000, 0.0110] | +0.0399 | -0.0255 | -0.0254 | 0.0801 | -0.0021 | -0.3173 | 0.8842 | 0.4844 | 0.0056 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 80 | +0.0036 | [-0.0007, +0.0074] | +0.0008 | 0.0000 | [0.0000, 0.0043] | +0.0008 | +0.0593 | +0.0592 | 0.0609 | +0.0036 | +0.9730 | 0.8813 | 0.4821 | 0.0055 |
| Confirmatory | vcc2025_h1 | `Transfer_Mean_VarScaled` | 80 | -0.1644 | [-0.2633, -0.0561] | -0.0772 | 0.0060 | [0.0000, 0.0378] | -0.0770 | -0.4446 | -0.4436 | 0.3707 | -0.1842 | -1.1968 | 0.8900 | 0.4810 | 0.0042 |
| Confirmatory | vcc2025_h1 | `Transfer_Mean_GlobalScale` | 80 | +0.0683 | [-0.0229, +0.1504] | +nan | nan | [nan, nan] | +nan | +0.2506 | +0.2501 | 0.2730 | +0.0734 | +0.9160 | 0.8932 | 0.4900 | 0.0039 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_MismatchedVar(nadig_hepg2)` | 80 | +0.0353 | [+0.0006, +0.0658] | +0.0629 | 0.0040 | [0.0000, 0.0169] | +0.0627 | +0.3153 | +0.3147 | 0.1122 | +0.0357 | +2.8055 | 0.8905 | 0.4878 | 0.0050 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_MismatchedVar(nadig_jurkat)` | 80 | -0.0416 | [-0.0747, -0.0049] | -0.0690 | 0.0048 | [0.0006, 0.0141] | -0.0688 | -0.3536 | -0.3528 | 0.1179 | -0.0421 | -2.9934 | 0.8908 | 0.4831 | 0.0058 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_MismatchedVar(nourreddine_kolf21j)` | 80 | -0.0471 | [-0.0917, +0.0038] | -0.0074 | 0.0001 | [0.0000, 0.0041] | -0.0074 | -0.2907 | -0.2901 | 0.1623 | -0.0483 | -1.7869 | 0.8927 | 0.4864 | 0.0060 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 80 | +0.0362 | [+0.0015, +0.0667] | +0.0601 | 0.0036 | [0.0003, 0.0106] | +0.0599 | +0.3323 | +0.3316 | 0.1093 | +0.0366 | +3.0338 | 0.8934 | 0.4895 | 0.0049 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 80 | +0.0384 | [-0.0024, +0.0755] | -0.0222 | 0.0005 | [0.0000, 0.0082] | -0.0221 | +0.3025 | +0.3019 | 0.1271 | +0.0389 | +2.3751 | 0.8924 | 0.4802 | 0.0054 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_ShuffledVar` | 80 | +0.0602 | [-0.0300, +0.1412] | +0.0343 | 0.0012 | [0.0000, 0.0074] | +0.0342 | +0.1892 | +0.1888 | 0.3187 | +0.0667 | +0.5924 | 0.8930 | 0.4876 | 0.0038 |
| Confirmatory | vcc2025_h1 | `Hybrid_on_TMVS` | 80 | -0.1733 | [-0.2722, -0.0623] | -0.0816 | 0.0067 | [0.0000, 0.0345] | -0.0813 | -0.4149 | -0.4140 | 0.4186 | -0.2016 | -0.9891 | 0.8247 | 0.3963 | 0.0051 |
| Confirmatory | vcc2025_h1 | `FCT` | 80 | +0.2503 | [+0.1339, +0.3542] | +0.3505 | 0.1228 | [0.0775, 0.1802] | +0.3495 | +0.5489 | +0.5478 | 0.4569 | +0.2898 | +1.1988 | 0.9134 | 0.5552 | 0.0069 |

## R18: TMVS Shuffled-Variance 100-Permutation Null vs H-B3 (97.5% CI Threshold)

| Target Line | TMVS r_orth_cal | TMVS CI lo | Perm100 p2.5 | Perm100 median | Perm100 p97.5 | Perm100 max | Max Mismatch CI hi | H-B3 Pass |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| nadig_hepg2 | +0.0418 | +0.0246 | +0.0122 | +0.0538 | +0.0872 | +0.1121 | +0.1140 | False |
| nadig_jurkat | +0.0947 | +0.0800 | +0.0576 | +0.1110 | +0.1529 | +0.1704 | +0.1731 | False |
| nourreddine_kolf21j | +0.1385 | -0.0508 | +0.0234 | +0.0810 | +0.1457 | +0.1761 | +0.1760 | False |
| replogle_k562_essential | +0.1822 | +0.1605 | +0.0198 | +0.0772 | +0.1247 | +0.1349 | +0.1775 | False |
| replogle_rpe1 | +0.0920 | +0.0725 | -0.0189 | +0.0176 | +0.0479 | +0.0560 | +0.0423 | True |
| vcc2025_h1 | -0.0772 | -0.1943 | -0.0223 | +0.0205 | +0.0534 | +0.0819 | +0.1298 | False |

## R19: Descriptive 3-Layer Residual Variance Decomposition (`R = Delta_t - Delta_src`, PREREG_A v4 §5.4)

| Target Line | Mask | N | Elig Genes | beta_ols_int [95% CI] | beta_mult | L1 Line Gain [95% CI] | L2 Per-Gene Gain [95% CI] | L3 Remaining [95% CI] |
|---|---|---:|---:|---|---:|---|---|---|
| nadig_hepg2 | `resp_top50` | 799 | 859 (C_dist>=5) | 1.1439 [1.0990, 1.1893] | 1.2080 | 0.0160 [0.0078, 0.0266] | 0.1413 [0.1259, 0.1559] | 0.8428 [0.8259, 0.8605] |
| nadig_hepg2 | `all_nontarget` | 799 | 6148 (C_dist>=5) | 0.5973 [0.5726, 0.6230] | 0.6055 | 0.1310 [0.1137, 0.1492] | 0.0731 [0.0671, 0.0793] | 0.7958 [0.7817, 0.8089] |
| nadig_jurkat | `resp_top50` | 934 | 1216 (C_dist>=5) | 0.9939 [0.9548, 1.0314] | 1.1063 | 0.0000 [0.0000, 0.0018] | 0.1488 [0.1283, 0.1661] | 0.8512 [0.8336, 0.8712] |
| nadig_jurkat | `all_nontarget` | 934 | 6148 (C_dist>=5) | 0.4469 [0.4272, 0.4669] | 0.4638 | 0.2301 [0.2124, 0.2479] | 0.0586 [0.0541, 0.0631] | 0.7113 [0.6961, 0.7265] |
| nourreddine_kolf21j | `resp_top50` | 18 | 12 (C_dist>=5) | 0.6432 [0.5067, 0.8921] | 0.6164 | 0.0999 [0.0069, 0.2496] | 0.0168 [0.0000, 0.0486] | 0.8833 [0.7428, 0.9893] |
| nourreddine_kolf21j | `all_nontarget` | 18 | 6148 (C_dist>=5) | 0.1832 [0.1212, 0.2530] | 0.1765 | 0.5129 [0.4073, 0.5919] | 0.0110 [0.0082, 0.0328] | 0.4761 [0.3865, 0.5721] |
| replogle_k562_essential | `resp_top50` | 859 | 1253 (C_dist>=5) | 0.7814 [0.7574, 0.8063] | 0.7813 | 0.0769 [0.0589, 0.0969] | 0.1790 [0.1647, 0.1915] | 0.7440 [0.7221, 0.7661] |
| replogle_k562_essential | `all_nontarget` | 859 | 6148 (C_dist>=5) | 0.3469 [0.3333, 0.3614] | 0.3475 | 0.4627 [0.4421, 0.4827] | 0.0439 [0.0404, 0.0477] | 0.4933 [0.4753, 0.5121] |
| replogle_rpe1 | `resp_top50` | 1070 | 823 (C_dist>=5) | 1.6321 [1.5747, 1.6920] | 1.9497 | 0.1239 [0.1058, 0.1428] | 0.1035 [0.0899, 0.1139] | 0.7725 [0.7545, 0.7942] |
| replogle_rpe1 | `all_nontarget` | 1070 | 6148 (C_dist>=5) | 0.8038 [0.7757, 0.8306] | 0.8349 | 0.0162 [0.0119, 0.0216] | 0.1300 [0.1227, 0.1368] | 0.8538 [0.8481, 0.8595] |
| vcc2025_h1 | `resp_top50` | 80 | 212 (C_dist>=5) | 0.6190 [0.3722, 1.1112] | 0.6144 | 0.0625 [0.0002, 0.2073] | 0.1886 [0.0958, 0.2190] | 0.7489 [0.6462, 0.8932] |
| vcc2025_h1 | `all_nontarget` | 80 | 6148 (C_dist>=5) | 0.2691 [0.2006, 0.3557] | 0.2713 | 0.3477 [0.2407, 0.4584] | 0.0367 [0.0254, 0.0535] | 0.6156 [0.5140, 0.7093] |

## S6: Split-Half Noise Ceiling (`resp_top50` vs Unselected `all_nontarget`, Amplitude-Orthogonalized & Raw)

| Target Line | N half-valid | rho_AB_orth (top50) | Ceil r_orth (top50) | rho_AB (top50) | Rep Slope (top50) | rho_AB_orth (all) | rho_AB (all) | TMVS r_orth_cal (all) | Hybrid r_orth_cal (all) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| nadig_hepg2 | 799 | 0.8745 | 0.9659 | 0.8784 | 0.8784 | 0.4747 | 0.5524 | +0.1824 | +0.0482 |
| nadig_jurkat | 934 | 0.8465 | 0.9575 | 0.8562 | 0.8562 | 0.3765 | 0.5096 | +0.1220 | +0.0464 |
| nourreddine_kolf21j | 18 | 0.8195 | 0.9491 | 0.8436 | 0.8438 | 0.1627 | 0.4900 | +0.0983 | -0.0148 |
| replogle_k562_essential | 859 | 0.8861 | 0.9693 | 0.9008 | 0.9008 | 0.4547 | 0.6830 | +0.1387 | +0.0213 |
| replogle_rpe1 | 1070 | 0.9134 | 0.9771 | 0.9178 | 0.9178 | 0.6139 | 0.6386 | +0.1762 | +0.1778 |
| vcc2025_h1 | 80 | 0.9885 | 0.9971 | 0.9915 | 0.9915 | 0.8886 | 0.9329 | +0.0323 | +0.0064 |

## S5: Expressed Genes Sensitivity (Drop Bottom 25% Target Control Variance Genes)

| Target Line | Genes Kept | rho_AB_orth (S5) | Method | Slope OLS (S5) | r_orth_cal (S5) | r_orth_cal^2 (S5) | Joint 97.5% CI | r_cal raw (S5) |
|---|---:|---:|---|---:|---:|---:|---|---:|
| nadig_hepg2 | 4611/6148 | 0.8744 | `Hybrid_ProjShared` | +0.0097 | +0.0489 | 0.0024 | [+0.0235, +0.0737] | +0.0717 |
| nadig_hepg2 | 4611/6148 | 0.8744 | `Transfer_Mean_VarScaled` | +0.0358 | +0.0423 | 0.0018 | [+0.0251, +0.0610] | +0.1229 |
| nadig_hepg2 | 4611/6148 | 0.8744 | `Transfer_Mean_GlobalScale` | -0.0084 | +nan | nan | [nan, nan] | -0.1295 |
| nadig_hepg2 | 4611/6148 | 0.8744 | `Hybrid_on_TMVS` | +0.0455 | +0.0553 | 0.0031 | [+0.0354, +0.0766] | +0.1291 |
| nadig_hepg2 | 4611/6148 | 0.8744 | `FCT` | +0.0591 | +0.2624 | 0.0688 | [+0.2432, +0.2820] | +0.2637 |
| nadig_hepg2 | 4611/6148 | 0.8744 | `NegCtrl_HybridProj_SrcCorrTgtVar` | +0.0035 | +0.0423 | 0.0018 | [+0.0292, +0.0559] | +0.0533 |
| nadig_hepg2 | 4611/6148 | 0.8744 | `NegCtrl_TMVS_ShuffledVar` | -0.0193 | +0.0590 | 0.0035 | [+0.0479, +0.0698] | -0.0612 |
| nadig_jurkat | 4611/6148 | 0.8465 | `Hybrid_ProjShared` | +0.0059 | +0.0657 | 0.0043 | [+0.0410, +0.0908] | +0.0574 |
| nadig_jurkat | 4611/6148 | 0.8465 | `Transfer_Mean_VarScaled` | +0.0093 | +0.0966 | 0.0093 | [+0.0825, +0.1112] | +0.0430 |
| nadig_jurkat | 4611/6148 | 0.8465 | `Transfer_Mean_GlobalScale` | -0.0013 | +nan | nan | [nan, nan] | -0.0094 |
| nadig_jurkat | 4611/6148 | 0.8465 | `Hybrid_on_TMVS` | +0.0151 | +0.1129 | 0.0127 | [+0.0941, +0.1314] | +0.0624 |
| nadig_jurkat | 4611/6148 | 0.8465 | `FCT` | +0.0198 | +0.1043 | 0.0109 | [+0.0845, +0.1239] | +0.0914 |
| nadig_jurkat | 4611/6148 | 0.8465 | `NegCtrl_HybridProj_SrcCorrTgtVar` | +0.0015 | +0.0281 | 0.0008 | [+0.0117, +0.0447] | +0.0240 |
| nadig_jurkat | 4611/6148 | 0.8465 | `NegCtrl_TMVS_ShuffledVar` | +0.0438 | +0.1285 | 0.0165 | [+0.1162, +0.1411] | +0.1490 |
| nourreddine_kolf21j | 4611/6148 | 0.8198 | `Hybrid_ProjShared` | -0.0138 | +0.0096 | 0.0001 | [-0.1361, +0.1355] | -0.0605 |
| nourreddine_kolf21j | 4611/6148 | 0.8198 | `Transfer_Mean_VarScaled` | -0.0133 | +0.1379 | 0.0190 | [-0.0542, +0.2869] | -0.0391 |
| nourreddine_kolf21j | 4611/6148 | 0.8198 | `Transfer_Mean_GlobalScale` | -0.0684 | +nan | nan | [nan, nan] | -0.3297 |
| nourreddine_kolf21j | 4611/6148 | 0.8198 | `Hybrid_on_TMVS` | -0.0271 | +0.1152 | 0.0133 | [-0.0518, +0.2709] | -0.0632 |
| nourreddine_kolf21j | 4611/6148 | 0.8198 | `FCT` | +0.0744 | +0.1835 | 0.0337 | [+0.0212, +0.3537] | +0.2538 |
| nourreddine_kolf21j | 4611/6148 | 0.8198 | `NegCtrl_HybridProj_SrcCorrTgtVar` | -0.0021 | +0.0996 | 0.0099 | [-0.0832, +0.2562] | -0.0109 |
| nourreddine_kolf21j | 4611/6148 | 0.8198 | `NegCtrl_TMVS_ShuffledVar` | +0.0336 | +0.0387 | 0.0015 | [-0.0080, +0.1066] | +0.1186 |
| replogle_k562_essential | 4611/6148 | 0.8847 | `Hybrid_ProjShared` | +0.0019 | +0.0339 | 0.0011 | [+0.0082, +0.0582] | +0.0150 |
| replogle_k562_essential | 4611/6148 | 0.8847 | `Transfer_Mean_VarScaled` | +0.0562 | +0.1868 | 0.0349 | [+0.1655, +0.2091] | +0.2383 |
| replogle_k562_essential | 4611/6148 | 0.8847 | `Transfer_Mean_GlobalScale` | +0.0028 | +nan | nan | [nan, nan] | +0.2804 |
| replogle_k562_essential | 4611/6148 | 0.8847 | `Hybrid_on_TMVS` | +0.0581 | +0.1764 | 0.0311 | [+0.1534, +0.1995] | +0.2161 |
| replogle_k562_essential | 4611/6148 | 0.8847 | `FCT` | +0.0424 | +0.1966 | 0.0386 | [+0.1778, +0.2170] | +0.1482 |
| replogle_k562_essential | 4611/6148 | 0.8847 | `NegCtrl_HybridProj_SrcCorrTgtVar` | +0.0018 | +0.0223 | 0.0005 | [+0.0021, +0.0409] | +0.0207 |
| replogle_k562_essential | 4611/6148 | 0.8847 | `NegCtrl_TMVS_ShuffledVar` | +0.1034 | +0.0656 | 0.0043 | [+0.0546, +0.0767] | +0.2483 |
| replogle_rpe1 | 4611/6148 | 0.9128 | `Hybrid_ProjShared` | +0.0135 | +0.1979 | 0.0391 | [+0.1722, +0.2232] | +0.1739 |
| replogle_rpe1 | 4611/6148 | 0.9128 | `Transfer_Mean_VarScaled` | +0.0351 | +0.0903 | 0.0082 | [+0.0724, +0.1077] | +0.2710 |
| replogle_rpe1 | 4611/6148 | 0.9128 | `Transfer_Mean_GlobalScale` | +0.0083 | +nan | nan | [nan, nan] | +0.3631 |
| replogle_rpe1 | 4611/6148 | 0.9128 | `Hybrid_on_TMVS` | +0.0486 | +0.1789 | 0.0320 | [+0.1582, +0.1994] | +0.3033 |
| replogle_rpe1 | 4611/6148 | 0.9128 | `FCT` | +0.0253 | +0.1318 | 0.0174 | [+0.1144, +0.1498] | +0.1815 |
| replogle_rpe1 | 4611/6148 | 0.9128 | `NegCtrl_HybridProj_SrcCorrTgtVar` | -0.0014 | +0.0001 | 0.0000 | [-0.0211, +0.0222] | -0.0350 |
| replogle_rpe1 | 4611/6148 | 0.9128 | `NegCtrl_TMVS_ShuffledVar` | -0.0356 | +0.0274 | 0.0007 | [+0.0178, +0.0373] | -0.1939 |
| vcc2025_h1 | 4611/6148 | 0.9885 | `Hybrid_ProjShared` | -0.0095 | -0.0387 | 0.0015 | [-0.1212, +0.0380] | -0.0530 |
| vcc2025_h1 | 4611/6148 | 0.9885 | `Transfer_Mean_VarScaled` | -0.1639 | -0.0791 | 0.0063 | [-0.1943, +0.0276] | -0.4486 |
| vcc2025_h1 | 4611/6148 | 0.9885 | `Transfer_Mean_GlobalScale` | +0.0688 | +nan | nan | [nan, nan] | +0.2556 |
| vcc2025_h1 | 4611/6148 | 0.9885 | `Hybrid_on_TMVS` | -0.1734 | -0.0820 | 0.0067 | [-0.1842, +0.0135] | -0.4197 |
| vcc2025_h1 | 4611/6148 | 0.9885 | `FCT` | +0.2509 | +0.3495 | 0.1221 | [+0.2751, +0.4252] | +0.5497 |
| vcc2025_h1 | 4611/6148 | 0.9885 | `NegCtrl_HybridProj_SrcCorrTgtVar` | +0.0035 | -0.0024 | 0.0000 | [-0.0560, +0.0614] | +0.0580 |
| vcc2025_h1 | 4611/6148 | 0.9885 | `NegCtrl_TMVS_ShuffledVar` | +0.0611 | +0.0324 | 0.0010 | [-0.0082, +0.0848] | +0.1935 |

## Appendix: Alternative H-A4/H-B4 Discrimination & Centroid Conventions (R5)

| Target Line | Method | Ctrl disc_raw (A-T2) | Ctrl disc_raw_panel (A-T3) | Cent disc_raw | Cent disc_raw_panel | Cent de50_r | pds_cosine_panel | Slope OLS (All) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| nadig_hepg2 | `Source_Grand_Mean` | 0.5003 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.3059 |
| nadig_hepg2 | `Transfer_Mean` | 0.8811 | 0.8811 | 0.9026 | 0.9028 | 0.5411 | 0.8790 | +0.0000 |
| nadig_hepg2 | `CIPHER_Linear` | 0.4832 | 0.4858 | 0.4893 | 0.4929 | -0.0750 | 0.4848 | +0.3257 |
| nadig_hepg2 | `Hybrid_Raw` | 0.8597 | 0.8606 | 0.8927 | 0.8932 | 0.5530 | 0.8575 | +0.0124 |
| nadig_hepg2 | `Hybrid_ProjShared` | 0.8688 | 0.8695 | 0.8975 | 0.8977 | 0.5452 | 0.8675 | +0.0080 |
| nadig_hepg2 | `NegCtrl_HybridProj_ShuffledGenes` | 0.8772 | 0.8761 | 0.8990 | 0.8995 | 0.5292 | 0.8748 | +0.0023 |
| nadig_hepg2 | `NegCtrl_HybridProj_DiagCov` | 0.8779 | 0.8771 | 0.8994 | 0.9000 | 0.5375 | 0.8756 | +0.0021 |
| nadig_hepg2 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.8808 | 0.8806 | 0.9023 | 0.9024 | 0.5409 | 0.8782 | +0.0020 |
| nadig_hepg2 | `Transfer_Mean_VarScaled` | 0.8745 | 0.8755 | 0.8998 | 0.9007 | 0.5398 | 0.8736 | +0.0329 |
| nadig_hepg2 | `Transfer_Mean_GlobalScale` | 0.8811 | 0.8811 | 0.9039 | 0.9039 | 0.5361 | 0.8790 | +0.0248 |
| nadig_hepg2 | `NegCtrl_TMVS_ShuffledVar` | 0.8902 | 0.8900 | 0.9011 | 0.9029 | 0.4723 | 0.8881 | +0.0625 |
| nadig_hepg2 | `Hybrid_on_TMVS` | 0.8611 | 0.8627 | 0.8938 | 0.8944 | 0.5367 | 0.8613 | +0.0410 |
| nadig_hepg2 | `FCT` | 0.8871 | 0.8883 | 0.9080 | 0.9094 | 0.5682 | 0.8852 | +0.0576 |
| nadig_jurkat | `Source_Grand_Mean` | 0.5002 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.3056 |
| nadig_jurkat | `Transfer_Mean` | 0.8717 | 0.8637 | 0.9218 | 0.9182 | 0.5323 | 0.8530 | +0.0000 |
| nadig_jurkat | `CIPHER_Linear` | 0.4832 | 0.4874 | 0.4936 | 0.4951 | 0.4025 | 0.4947 | +0.4094 |
| nadig_jurkat | `Hybrid_Raw` | 0.8743 | 0.8673 | 0.9228 | 0.9191 | 0.5654 | 0.8560 | +0.0209 |
| nadig_jurkat | `Hybrid_ProjShared` | 0.8685 | 0.8609 | 0.9191 | 0.9156 | 0.5355 | 0.8507 | +0.0051 |
| nadig_jurkat | `NegCtrl_HybridProj_ShuffledGenes` | 0.8706 | 0.8619 | 0.9197 | 0.9161 | 0.5215 | 0.8516 | -0.0009 |
| nadig_jurkat | `NegCtrl_HybridProj_DiagCov` | 0.8712 | 0.8624 | 0.9202 | 0.9167 | 0.5237 | 0.8520 | -0.0011 |
| nadig_jurkat | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.8709 | 0.8629 | 0.9217 | 0.9182 | 0.5333 | 0.8524 | +0.0014 |
| nadig_jurkat | `Transfer_Mean_VarScaled` | 0.8709 | 0.8629 | 0.9135 | 0.9100 | 0.5053 | 0.8514 | -0.0493 |
| nadig_jurkat | `Transfer_Mean_GlobalScale` | 0.8717 | 0.8637 | 0.9141 | 0.9104 | 0.4919 | 0.8530 | -0.0645 |
| nadig_jurkat | `NegCtrl_TMVS_ShuffledVar` | 0.8813 | 0.8756 | 0.9308 | 0.9283 | 0.5930 | 0.8661 | +0.0074 |
| nadig_jurkat | `Hybrid_on_TMVS` | 0.8683 | 0.8607 | 0.9115 | 0.9076 | 0.5091 | 0.8495 | -0.0442 |
| nadig_jurkat | `FCT` | 0.8603 | 0.8519 | 0.9042 | 0.9006 | 0.5139 | 0.8359 | -0.0071 |
| nourreddine_kolf21j | `Source_Grand_Mean` | 0.5000 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.4560 |
| nourreddine_kolf21j | `Transfer_Mean` | 0.8137 | 0.8105 | 0.7778 | 0.7778 | 0.4715 | 0.8072 | +0.0000 |
| nourreddine_kolf21j | `CIPHER_Linear` | 0.5817 | 0.5882 | 0.5196 | 0.5196 | 0.6784 | 0.5098 | +0.6231 |
| nourreddine_kolf21j | `Hybrid_Raw` | 0.7810 | 0.7810 | 0.7516 | 0.7516 | 0.5053 | 0.7451 | +0.0212 |
| nourreddine_kolf21j | `Hybrid_ProjShared` | 0.8007 | 0.7974 | 0.7451 | 0.7451 | 0.4517 | 0.7941 | -0.0018 |
| nourreddine_kolf21j | `NegCtrl_HybridProj_ShuffledGenes` | 0.7614 | 0.7614 | 0.7222 | 0.7255 | 0.4352 | 0.7680 | +0.0024 |
| nourreddine_kolf21j | `NegCtrl_HybridProj_DiagCov` | 0.7647 | 0.7680 | 0.7255 | 0.7255 | 0.4470 | 0.7712 | +0.0030 |
| nourreddine_kolf21j | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.8007 | 0.8039 | 0.7451 | 0.7451 | 0.4716 | 0.7941 | +0.0027 |
| nourreddine_kolf21j | `Transfer_Mean_VarScaled` | 0.8105 | 0.8105 | 0.7647 | 0.7647 | 0.5042 | 0.7941 | -0.1109 |
| nourreddine_kolf21j | `Transfer_Mean_GlobalScale` | 0.8137 | 0.8105 | 0.7745 | 0.7745 | 0.3784 | 0.8072 | -0.1538 |
| nourreddine_kolf21j | `NegCtrl_TMVS_ShuffledVar` | 0.8105 | 0.8105 | 0.7712 | 0.7712 | 0.5577 | 0.8007 | -0.0326 |
| nourreddine_kolf21j | `Hybrid_on_TMVS` | 0.8007 | 0.7974 | 0.7614 | 0.7614 | 0.4906 | 0.7908 | -0.1127 |
| nourreddine_kolf21j | `FCT` | 0.8170 | 0.8170 | 0.7712 | 0.7778 | 0.5354 | 0.7810 | +0.0620 |
| replogle_k562_essential | `Source_Grand_Mean` | 0.5001 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.4825 |
| replogle_k562_essential | `Transfer_Mean` | 0.8868 | 0.8927 | 0.9374 | 0.9390 | 0.4539 | 0.8788 | +0.0000 |
| replogle_k562_essential | `CIPHER_Linear` | 0.4382 | 0.4424 | 0.4878 | 0.4894 | 0.4603 | 0.4612 | +0.6932 |
| replogle_k562_essential | `Hybrid_Raw` | 0.8931 | 0.8974 | 0.9480 | 0.9487 | 0.5117 | 0.8892 | +0.0465 |
| replogle_k562_essential | `Hybrid_ProjShared` | 0.8840 | 0.8888 | 0.9364 | 0.9374 | 0.4488 | 0.8750 | +0.0019 |
| replogle_k562_essential | `NegCtrl_HybridProj_ShuffledGenes` | 0.8930 | 0.8957 | 0.9375 | 0.9379 | 0.4326 | 0.8829 | -0.0032 |
| replogle_k562_essential | `NegCtrl_HybridProj_DiagCov` | 0.8931 | 0.8960 | 0.9376 | 0.9382 | 0.4340 | 0.8831 | -0.0034 |
| replogle_k562_essential | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.8870 | 0.8924 | 0.9372 | 0.9389 | 0.4572 | 0.8787 | +0.0006 |
| replogle_k562_essential | `Transfer_Mean_VarScaled` | 0.8885 | 0.8948 | 0.9425 | 0.9449 | 0.5267 | 0.8801 | +0.0408 |
| replogle_k562_essential | `Transfer_Mean_GlobalScale` | 0.8868 | 0.8927 | 0.9380 | 0.9396 | 0.4569 | 0.8788 | +0.0058 |
| replogle_k562_essential | `NegCtrl_TMVS_ShuffledVar` | 0.8921 | 0.8995 | 0.9541 | 0.9552 | 0.5777 | 0.8880 | +0.1132 |
| replogle_k562_essential | `Hybrid_on_TMVS` | 0.8856 | 0.8909 | 0.9414 | 0.9432 | 0.5183 | 0.8762 | +0.0427 |
| replogle_k562_essential | `FCT` | 0.8763 | 0.8846 | 0.9309 | 0.9368 | 0.4754 | 0.8610 | +0.0440 |
| replogle_rpe1 | `Source_Grand_Mean` | 0.5005 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.1628 |
| replogle_rpe1 | `Transfer_Mean` | 0.9301 | 0.9346 | 0.9173 | 0.9202 | 0.4525 | 0.9287 | +0.0000 |
| replogle_rpe1 | `CIPHER_Linear` | 0.5120 | 0.5157 | 0.5123 | 0.5144 | -0.2142 | 0.5164 | +0.1242 |
| replogle_rpe1 | `Hybrid_Raw` | 0.9200 | 0.9235 | 0.9316 | 0.9341 | 0.5099 | 0.9179 | +0.0225 |
| replogle_rpe1 | `Hybrid_ProjShared` | 0.9248 | 0.9291 | 0.9277 | 0.9306 | 0.4921 | 0.9237 | +0.0161 |
| replogle_rpe1 | `NegCtrl_HybridProj_ShuffledGenes` | 0.9283 | 0.9322 | 0.9087 | 0.9114 | 0.4214 | 0.9265 | -0.0082 |
| replogle_rpe1 | `NegCtrl_HybridProj_DiagCov` | 0.9277 | 0.9316 | 0.9085 | 0.9111 | 0.4260 | 0.9261 | -0.0080 |
| replogle_rpe1 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.9293 | 0.9333 | 0.9157 | 0.9186 | 0.4492 | 0.9275 | -0.0009 |
| replogle_rpe1 | `Transfer_Mean_VarScaled` | 0.9254 | 0.9302 | 0.9238 | 0.9256 | 0.4894 | 0.9251 | +0.0132 |
| replogle_rpe1 | `Transfer_Mean_GlobalScale` | 0.9301 | 0.9346 | 0.9208 | 0.9234 | 0.4611 | 0.9287 | -0.0036 |
| replogle_rpe1 | `NegCtrl_TMVS_ShuffledVar` | 0.9344 | 0.9387 | 0.9018 | 0.9074 | 0.3555 | 0.9327 | -0.0071 |
| replogle_rpe1 | `Hybrid_on_TMVS` | 0.9208 | 0.9258 | 0.9313 | 0.9338 | 0.5189 | 0.9207 | +0.0293 |
| replogle_rpe1 | `FCT` | 0.9308 | 0.9347 | 0.9181 | 0.9220 | 0.4648 | 0.9290 | +0.0358 |
| vcc2025_h1 | `Source_Grand_Mean` | 0.5003 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.3627 |
| vcc2025_h1 | `Transfer_Mean` | 0.8932 | 0.8921 | 0.8354 | 0.8362 | 0.5056 | 0.8758 | +0.0000 |
| vcc2025_h1 | `CIPHER_Linear` | 0.4725 | 0.4718 | 0.4826 | 0.4826 | 0.3693 | 0.4663 | +0.4676 |
| vcc2025_h1 | `Hybrid_Raw` | 0.8579 | 0.8590 | 0.8247 | 0.8266 | 0.4682 | 0.8538 | +0.0000 |
| vcc2025_h1 | `Hybrid_ProjShared` | 0.8513 | 0.8511 | 0.8293 | 0.8294 | 0.4572 | 0.8402 | -0.0001 |
| vcc2025_h1 | `NegCtrl_HybridProj_ShuffledGenes` | 0.8828 | 0.8812 | 0.8342 | 0.8334 | 0.4814 | 0.8698 | -0.0004 |
| vcc2025_h1 | `NegCtrl_HybridProj_DiagCov` | 0.8842 | 0.8834 | 0.8331 | 0.8335 | 0.4963 | 0.8709 | -0.0001 |
| vcc2025_h1 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.8813 | 0.8778 | 0.8345 | 0.8356 | 0.4958 | 0.8649 | +0.0006 |
| vcc2025_h1 | `Transfer_Mean_VarScaled` | 0.8900 | 0.8886 | 0.7253 | 0.7228 | 0.5599 | 0.8750 | +0.1924 |
| vcc2025_h1 | `Transfer_Mean_GlobalScale` | 0.8932 | 0.8921 | 0.7328 | 0.7334 | 0.5563 | 0.8758 | +0.1979 |
| vcc2025_h1 | `NegCtrl_TMVS_ShuffledVar` | 0.8930 | 0.8911 | 0.6888 | 0.6892 | 0.5681 | 0.8771 | +0.2344 |
| vcc2025_h1 | `Hybrid_on_TMVS` | 0.8247 | 0.8245 | 0.6954 | 0.6964 | 0.4871 | 0.8152 | +0.1923 |
| vcc2025_h1 | `FCT` | 0.9134 | 0.9123 | 0.8905 | 0.8905 | 0.6073 | 0.8998 | +0.0427 |
