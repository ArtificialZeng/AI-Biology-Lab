# Line-Specific Zero-Shot Perturbation Evaluation Report (`holdout_log1p_sens_s2b_rank1q1`)

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
| Confirmatory | nadig_hepg2 | `CIPHER_Linear` | 799 | -0.0670 | [-0.1014, -0.0320] | +0.1174 | 0.0138 | [0.0072, 0.0226] | +0.1134 | -0.0830 | -0.0803 | 0.8348 | -0.2112 | -0.0961 | 0.4832 | 0.2594 | 0.0088 |
| Confirmatory | nadig_hepg2 | `Hybrid_Raw` | 799 | +0.0352 | [+0.0286, +0.0421] | +0.1226 | 0.0150 | [0.0079, 0.0246] | +0.1185 | +0.1831 | +0.1771 | 0.1987 | +0.0366 | +0.8913 | 0.8597 | 0.7050 | 0.0081 |
| Confirmatory | nadig_hepg2 | `Hybrid_ProjShared` | 799 | +0.0089 | [+0.0054, +0.0127] | +0.0468 | 0.0022 | [0.0004, 0.0052] | +0.0452 | +0.0659 | +0.0637 | 0.1401 | +0.0091 | +0.4547 | 0.8696 | 0.6971 | 0.0081 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_MismatchedCov(nadig_jurkat)` | 799 | -0.0024 | [-0.0036, -0.0011] | -0.0127 | 0.0002 | [0.0000, 0.0008] | -0.0122 | -0.0325 | -0.0314 | 0.0754 | -0.0024 | -0.4172 | 0.8818 | 0.6991 | 0.0080 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_MismatchedCov(nourreddine_kolf21j)` | 799 | +0.0085 | [+0.0064, +0.0108] | +0.0573 | 0.0033 | [0.0019, 0.0052] | +0.0553 | +0.0759 | +0.0734 | 0.1160 | +0.0086 | +0.6329 | 0.8811 | 0.7015 | 0.0081 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 799 | +0.0047 | [+0.0034, +0.0059] | +0.0481 | 0.0023 | [0.0012, 0.0038] | +0.0464 | +0.0634 | +0.0613 | 0.0760 | +0.0047 | +0.8070 | 0.8785 | 0.7056 | 0.0079 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 799 | -0.0024 | [-0.0036, -0.0012] | -0.0192 | 0.0004 | [0.0000, 0.0012] | -0.0185 | -0.0326 | -0.0315 | 0.0771 | -0.0024 | -0.4087 | 0.8829 | 0.6970 | 0.0080 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_MismatchedCov(vcc2025_h1)` | 799 | +0.0070 | [+0.0043, +0.0096] | +0.0370 | 0.0014 | [0.0005, 0.0026] | +0.0357 | +0.0398 | +0.0385 | 0.1830 | +0.0073 | +0.2104 | 0.8741 | 0.6866 | 0.0087 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_ShuffledGenes` | 799 | +0.0031 | [+0.0003, +0.0061] | +0.0166 | 0.0003 | [0.0000, 0.0010] | +0.0160 | +0.0212 | +0.0205 | 0.1507 | +0.0032 | +0.1360 | 0.8816 | 0.6889 | 0.0087 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_DiagCov` | 799 | +0.0141 | [+0.0116, +0.0168] | +0.0771 | 0.0059 | [0.0036, 0.0091] | +0.0745 | +0.1231 | +0.1191 | 0.1186 | +0.0143 | +1.0038 | 0.8742 | 0.7052 | 0.0080 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 799 | +0.0051 | [+0.0039, +0.0062] | +0.0586 | 0.0034 | [0.0020, 0.0054] | +0.0566 | +0.0753 | +0.0729 | 0.0699 | +0.0051 | +1.0430 | 0.8810 | 0.7067 | 0.0079 |
| Confirmatory | nadig_hepg2 | `Transfer_Mean_VarScaled` | 799 | +0.0357 | [+0.0273, +0.0441] | +0.0419 | 0.0018 | [0.0006, 0.0037] | +0.0404 | +0.1224 | +0.1184 | 0.3017 | +0.0392 | +0.3924 | 0.8745 | 0.6900 | 0.0076 |
| Confirmatory | nadig_hepg2 | `Transfer_Mean_GlobalScale` | 799 | -0.0085 | [-0.0111, -0.0057] | +nan | nan | [nan, nan] | +nan | -0.1306 | -0.1263 | 0.0670 | -0.0085 | -1.8867 | 0.8811 | 0.7053 | 0.0075 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_MismatchedVar(nadig_jurkat)` | 799 | +0.0076 | [+0.0038, +0.0112] | -0.0620 | 0.0038 | [0.0025, 0.0055] | -0.0598 | +0.0503 | +0.0486 | 0.1559 | +0.0078 | +0.3117 | 0.8803 | 0.6929 | 0.0087 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_MismatchedVar(nourreddine_kolf21j)` | 799 | +0.0088 | [+0.0049, +0.0121] | +0.0471 | 0.0022 | [0.0011, 0.0036] | +0.0455 | +0.0493 | +0.0477 | 0.1837 | +0.0091 | +0.2597 | 0.8859 | 0.6990 | 0.0088 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 799 | +0.0085 | [+0.0044, +0.0136] | +0.0951 | 0.0091 | [0.0059, 0.0130] | +0.0919 | +0.0599 | +0.0580 | 0.1466 | +0.0087 | +0.3952 | 0.8831 | 0.7043 | 0.0078 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 799 | +0.0150 | [+0.0071, +0.0227] | -0.0401 | 0.0016 | [0.0003, 0.0039] | -0.0387 | +0.0694 | +0.0672 | 0.2232 | +0.0158 | +0.3009 | 0.8731 | 0.6888 | 0.0083 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_MismatchedVar(vcc2025_h1)` | 799 | -0.0579 | [-0.0725, -0.0430] | -0.0181 | 0.0003 | [0.0000, 0.0018] | -0.0175 | -0.1402 | -0.1356 | 0.4272 | -0.0705 | -0.3175 | 0.8823 | 0.6728 | 0.0070 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_ShuffledVar` | 799 | -0.0196 | [-0.0289, -0.0099] | +0.0589 | 0.0035 | [0.0023, 0.0049] | +0.0569 | -0.0621 | -0.0601 | 0.3268 | -0.0220 | -0.1838 | 0.8902 | 0.6706 | 0.0077 |
| Confirmatory | nadig_hepg2 | `Hybrid_on_TMVS` | 799 | +0.0446 | [+0.0349, +0.0547] | +0.0544 | 0.0030 | [0.0012, 0.0057] | +0.0526 | +0.1278 | +0.1236 | 0.3613 | +0.0512 | +0.3420 | 0.8618 | 0.6794 | 0.0079 |
| Confirmatory | nadig_hepg2 | `FCT` | 799 | +0.0588 | [+0.0536, +0.0640] | +0.2611 | 0.0682 | [0.0585, 0.0788] | +0.2522 | +0.2618 | +0.2532 | 0.2322 | +0.0619 | +1.0905 | 0.8871 | 0.7269 | 0.0076 |
| Confirmatory | nadig_jurkat | `Source_Grand_Mean` | 934 | -0.0036 | [-0.0349, +0.0310] | +0.0153 | 0.0002 | [0.0000, 0.0020] | +0.0147 | -0.0049 | -0.0047 | 0.7616 | -0.0086 | -0.0062 | 0.5002 | 0.3335 | 0.0068 |
| Confirmatory | nadig_jurkat | `Transfer_Mean` | 934 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.8717 | 0.6172 | 0.0076 |
| Confirmatory | nadig_jurkat | `CIPHER_Linear` | 934 | +0.0067 | [-0.0288, +0.0460] | +0.0440 | 0.0019 | [0.0001, 0.0060] | +0.0421 | +0.0076 | +0.0073 | 0.9149 | +0.0407 | +0.0079 | 0.4832 | 0.0743 | 0.0070 |
| Confirmatory | nadig_jurkat | `Hybrid_Raw` | 934 | +0.0200 | [+0.0141, +0.0267] | +0.1290 | 0.0166 | [0.0097, 0.0258] | +0.1235 | +0.1478 | +0.1419 | 0.1408 | +0.0204 | +1.0084 | 0.8743 | 0.6225 | 0.0076 |
| Confirmatory | nadig_jurkat | `Hybrid_ProjShared` | 934 | +0.0102 | [+0.0070, +0.0136] | +0.1040 | 0.0108 | [0.0061, 0.0168] | +0.0996 | +0.0950 | +0.0912 | 0.1113 | +0.0103 | +0.8193 | 0.8679 | 0.6193 | 0.0077 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_MismatchedCov(nadig_hepg2)` | 934 | -0.0042 | [-0.0061, -0.0023] | -0.0300 | 0.0009 | [0.0001, 0.0024] | -0.0287 | -0.0432 | -0.0415 | 0.1003 | -0.0042 | -0.4141 | 0.8682 | 0.6063 | 0.0078 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_MismatchedCov(nourreddine_kolf21j)` | 934 | +0.0014 | [-0.0011, +0.0038] | +0.0151 | 0.0002 | [0.0000, 0.0014] | +0.0145 | +0.0119 | +0.0114 | 0.1200 | +0.0014 | +0.0949 | 0.8717 | 0.6081 | 0.0079 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 934 | +0.0043 | [+0.0024, +0.0063] | +0.0529 | 0.0028 | [0.0011, 0.0051] | +0.0507 | +0.0507 | +0.0487 | 0.0891 | +0.0044 | +0.5468 | 0.8689 | 0.6163 | 0.0078 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 934 | -0.0016 | [-0.0031, +0.0001] | -0.0272 | 0.0007 | [0.0001, 0.0022] | -0.0261 | -0.0217 | -0.0209 | 0.0745 | -0.0016 | -0.2798 | 0.8714 | 0.6108 | 0.0078 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_MismatchedCov(vcc2025_h1)` | 934 | +0.0102 | [+0.0052, +0.0153] | +0.0464 | 0.0022 | [0.0006, 0.0046] | +0.0445 | +0.0528 | +0.0507 | 0.2008 | +0.0106 | +0.2524 | 0.8658 | 0.6010 | 0.0084 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_ShuffledGenes` | 934 | +0.0039 | [+0.0011, +0.0068] | +0.0263 | 0.0007 | [0.0000, 0.0023] | +0.0252 | +0.0305 | +0.0293 | 0.1319 | +0.0039 | +0.2219 | 0.8717 | 0.6080 | 0.0081 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_DiagCov` | 934 | +0.0017 | [-0.0012, +0.0045] | -0.0014 | 0.0000 | [0.0000, 0.0008] | -0.0014 | +0.0150 | +0.0144 | 0.1148 | +0.0017 | +0.1258 | 0.8660 | 0.6102 | 0.0079 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_SrcCorrTgtVar` | 934 | +0.0026 | [+0.0015, +0.0039] | +0.0433 | 0.0019 | [0.0007, 0.0037] | +0.0414 | +0.0402 | +0.0386 | 0.0685 | +0.0027 | +0.5633 | 0.8696 | 0.6160 | 0.0077 |
| Confirmatory | nadig_jurkat | `Transfer_Mean_VarScaled` | 934 | +0.0095 | [+0.0028, +0.0157] | +0.0947 | 0.0090 | [0.0064, 0.0119] | +0.0906 | +0.0439 | +0.0422 | 0.2246 | +0.0100 | +0.1878 | 0.8709 | 0.6187 | 0.0087 |
| Confirmatory | nadig_jurkat | `Transfer_Mean_GlobalScale` | 934 | -0.0008 | [-0.0068, +0.0047] | +nan | nan | [nan, nan] | +nan | -0.0058 | -0.0056 | 0.1426 | -0.0008 | -0.0393 | 0.8717 | 0.6172 | 0.0088 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_MismatchedVar(nadig_hepg2)` | 934 | -0.0370 | [-0.0438, -0.0300] | -0.1779 | 0.0316 | [0.0258, 0.0383] | -0.1703 | -0.1784 | -0.1713 | 0.2161 | -0.0388 | -0.7928 | 0.8671 | 0.5856 | 0.0079 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_MismatchedVar(nourreddine_kolf21j)` | 934 | +0.0417 | [+0.0332, +0.0488] | +0.1587 | 0.0252 | [0.0208, 0.0299] | +0.1520 | +0.1724 | +0.1655 | 0.2518 | +0.0444 | +0.6575 | 0.8798 | 0.6284 | 0.0088 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 934 | +0.0141 | [+0.0076, +0.0211] | +0.1061 | 0.0113 | [0.0077, 0.0153] | +0.1016 | +0.0893 | +0.0857 | 0.1644 | +0.0145 | +0.5215 | 0.8745 | 0.6215 | 0.0076 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 934 | -0.0309 | [-0.0382, -0.0234] | -0.1368 | 0.0187 | [0.0136, 0.0252] | -0.1310 | -0.1675 | -0.1609 | 0.1920 | -0.0320 | -0.8384 | 0.8629 | 0.5927 | 0.0086 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_MismatchedVar(vcc2025_h1)` | 934 | +0.0187 | [+0.0055, +0.0318] | +0.0364 | 0.0013 | [0.0002, 0.0039] | +0.0349 | +0.0423 | +0.0406 | 0.4605 | +0.0237 | +0.0882 | 0.8728 | 0.5809 | 0.0062 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_ShuffledVar` | 934 | +0.0440 | [+0.0380, +0.0503] | +0.1290 | 0.0166 | [0.0136, 0.0202] | +0.1235 | +0.1495 | +0.1436 | 0.3063 | +0.0484 | +0.4690 | 0.8814 | 0.6107 | 0.0082 |
| Confirmatory | nadig_jurkat | `Hybrid_on_TMVS` | 934 | +0.0196 | [+0.0126, +0.0263] | +0.1321 | 0.0175 | [0.0128, 0.0228] | +0.1265 | +0.0803 | +0.0772 | 0.2544 | +0.0210 | +0.3033 | 0.8678 | 0.6208 | 0.0088 |
| Confirmatory | nadig_jurkat | `FCT` | 934 | +0.0200 | [+0.0149, +0.0253] | +0.1043 | 0.0109 | [0.0072, 0.0153] | +0.0999 | +0.0915 | +0.0879 | 0.2274 | +0.0211 | +0.3865 | 0.8603 | 0.6115 | 0.0085 |
| Secondary | nourreddine_kolf21j | `Source_Grand_Mean` | 18 | +0.2723 | [+0.0345, +0.5563] | +0.0672 | 0.0045 | [0.0000, 0.0395] | +0.0638 | +0.3484 | +0.3333 | 0.8168 | +0.5608 | +0.4081 | 0.5000 | 0.2906 | 0.0034 |
| Secondary | nourreddine_kolf21j | `Transfer_Mean` | 18 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.8137 | 0.4335 | 0.0042 |
| Secondary | nourreddine_kolf21j | `CIPHER_Linear` | 18 | +0.2738 | [+0.0234, +0.5581] | -0.1096 | 0.0120 | [0.0000, 0.0714] | -0.1040 | +0.3147 | +0.3010 | 0.9096 | +0.7333 | +0.3309 | 0.5817 | -0.1068 | 0.0023 |
| Secondary | nourreddine_kolf21j | `Hybrid_Raw` | 18 | -0.0129 | [-0.0552, +0.0227] | +0.0058 | 0.0000 | [0.0000, 0.0201] | +0.0055 | -0.0506 | -0.0484 | 0.2666 | -0.0139 | -0.1817 | 0.7810 | 0.4346 | 0.0045 |
| Secondary | nourreddine_kolf21j | `Hybrid_ProjShared` | 18 | -0.0164 | [-0.0579, +0.0190] | -0.0130 | 0.0002 | [0.0000, 0.0234] | -0.0123 | -0.0710 | -0.0679 | 0.2412 | -0.0174 | -0.2816 | 0.7941 | 0.4299 | 0.0047 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_MismatchedCov(nadig_hepg2)` | 18 | +0.0027 | [-0.0447, +0.0465] | +0.0912 | 0.0083 | [0.0000, 0.0449] | +0.0865 | +0.0116 | +0.0111 | 0.2430 | +0.0029 | +0.0457 | 0.8105 | 0.4552 | 0.0047 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_MismatchedCov(nadig_jurkat)` | 18 | -0.0089 | [-0.0571, +0.0298] | +0.0277 | 0.0008 | [0.0000, 0.0232] | +0.0263 | -0.0407 | -0.0389 | 0.2289 | -0.0094 | -0.1701 | 0.7810 | 0.4360 | 0.0046 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 18 | -0.0089 | [-0.0624, +0.0373] | +0.0465 | 0.0022 | [0.0000, 0.0323] | +0.0441 | -0.0383 | -0.0366 | 0.2425 | -0.0094 | -0.1511 | 0.7843 | 0.4482 | 0.0045 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 18 | -0.0181 | [-0.0740, +0.0208] | -0.0105 | 0.0001 | [0.0000, 0.0159] | -0.0100 | -0.0800 | -0.0766 | 0.2371 | -0.0192 | -0.3229 | 0.7974 | 0.4271 | 0.0046 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_MismatchedCov(vcc2025_h1)` | 18 | +0.0118 | [-0.0219, +0.0586] | +0.0662 | 0.0044 | [0.0000, 0.0526] | +0.0628 | +0.0568 | +0.0543 | 0.2178 | +0.0124 | +0.2493 | 0.8137 | 0.4331 | 0.0046 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_ShuffledGenes` | 18 | -0.0101 | [-0.0600, +0.0373] | +0.0440 | 0.0019 | [0.0000, 0.0338] | +0.0417 | -0.0396 | -0.0379 | 0.2670 | -0.0109 | -0.1421 | 0.7484 | 0.4421 | 0.0048 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_DiagCov` | 18 | -0.0113 | [-0.0624, +0.0363] | +0.0517 | 0.0027 | [0.0000, 0.0418] | +0.0491 | -0.0471 | -0.0450 | 0.2499 | -0.0120 | -0.1801 | 0.7549 | 0.4481 | 0.0045 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_SrcCorrTgtVar` | 18 | -0.0017 | [-0.0482, +0.0394] | +0.0970 | 0.0094 | [0.0000, 0.0649] | +0.0921 | -0.0089 | -0.0085 | 0.1965 | -0.0017 | -0.0431 | 0.7974 | 0.4609 | 0.0044 |
| Secondary | nourreddine_kolf21j | `Transfer_Mean_VarScaled` | 18 | -0.0133 | [-0.0927, +0.0532] | +0.1385 | 0.0192 | [0.0000, 0.0824] | +0.1314 | -0.0392 | -0.0375 | 0.3554 | -0.0152 | -0.1054 | 0.8105 | 0.4380 | 0.0055 |
| Secondary | nourreddine_kolf21j | `Transfer_Mean_GlobalScale` | 18 | -0.0686 | [-0.1382, -0.0098] | +nan | nan | [nan, nan] | +nan | -0.3304 | -0.3161 | 0.2170 | -0.0716 | -1.4568 | 0.8137 | 0.4335 | 0.0057 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_MismatchedVar(nadig_hepg2)` | 18 | +0.0172 | [-0.0138, +0.0687] | +0.0173 | 0.0003 | [0.0000, 0.0273] | +0.0165 | +0.1185 | +0.1134 | 0.1514 | +0.0176 | +0.7489 | 0.8203 | 0.4238 | 0.0042 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_MismatchedVar(nadig_jurkat)` | 18 | -0.0405 | [-0.0734, -0.0031] | +0.0486 | 0.0024 | [0.0000, 0.0124] | +0.0461 | -0.2242 | -0.2145 | 0.1889 | -0.0420 | -1.1355 | 0.8039 | 0.4343 | 0.0052 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 18 | +0.0091 | [-0.0028, +0.0225] | +0.1200 | 0.0144 | [0.0029, 0.0310] | +0.1139 | +0.0745 | +0.0713 | 0.1270 | +0.0092 | +0.5611 | 0.8105 | 0.4385 | 0.0043 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 18 | +0.0052 | [-0.0188, +0.0311] | +0.0278 | 0.0008 | [0.0000, 0.0206] | +0.0263 | +0.0335 | +0.0320 | 0.1624 | +0.0053 | +0.1973 | 0.8137 | 0.4309 | 0.0048 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_MismatchedVar(vcc2025_h1)` | 18 | +0.0488 | [+0.0028, +0.0907] | -0.1371 | 0.0188 | [0.0000, 0.0728] | -0.1301 | +0.1332 | +0.1274 | 0.3828 | +0.0570 | +0.3329 | 0.8007 | 0.3960 | 0.0027 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_ShuffledVar` | 18 | +0.0332 | [-0.0145, +0.0885] | +0.0388 | 0.0015 | [0.0000, 0.0114] | +0.0369 | +0.1166 | +0.1115 | 0.2979 | +0.0364 | +0.3744 | 0.8105 | 0.4233 | 0.0049 |
| Secondary | nourreddine_kolf21j | `Hybrid_on_TMVS` | 18 | -0.0297 | [-0.1264, +0.0509] | +0.1015 | 0.0103 | [0.0000, 0.0665] | +0.0963 | -0.0700 | -0.0670 | 0.4433 | -0.0369 | -0.1511 | 0.7974 | 0.4410 | 0.0060 |
| Secondary | nourreddine_kolf21j | `FCT` | 18 | +0.0747 | [+0.0070, +0.1592] | +0.1839 | 0.0338 | [0.0006, 0.1256] | +0.1746 | +0.2547 | +0.2436 | 0.3066 | +0.0819 | +0.7947 | 0.8170 | 0.4414 | 0.0044 |
| Secondary | replogle_k562_essential | `Source_Grand_Mean` | 859 | +0.2588 | [+0.2086, +0.3093] | -0.0643 | 0.0041 | [0.0017, 0.0076] | -0.0623 | +0.2554 | +0.2487 | 1.0408 | +1.1736 | +0.2390 | 0.5001 | 0.3941 | 0.0052 |
| Secondary | replogle_k562_essential | `Transfer_Mean` | 859 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.8868 | 0.6788 | 0.0073 |
| Secondary | replogle_k562_essential | `CIPHER_Linear` | 859 | +0.3499 | [+0.2931, +0.4075] | +0.0256 | 0.0007 | [0.0000, 0.0027] | +0.0248 | +0.2843 | +0.2768 | 1.2641 | +2.1697 | +0.2190 | 0.4382 | 0.0712 | 0.0049 |
| Secondary | replogle_k562_essential | `Hybrid_Raw` | 859 | +0.0312 | [+0.0254, +0.0372] | +0.1044 | 0.0109 | [0.0056, 0.0172] | +0.1012 | +0.1665 | +0.1621 | 0.1925 | +0.0324 | +0.8423 | 0.8931 | 0.6811 | 0.0069 |
| Secondary | replogle_k562_essential | `Hybrid_ProjShared` | 859 | +0.0016 | [-0.0018, +0.0048] | +0.0435 | 0.0019 | [0.0002, 0.0048] | +0.0421 | +0.0121 | +0.0118 | 0.1385 | +0.0017 | +0.0850 | 0.8806 | 0.6767 | 0.0075 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(nadig_hepg2)` | 859 | -0.0084 | [-0.0112, -0.0060] | -0.0738 | 0.0054 | [0.0033, 0.0081] | -0.0715 | -0.0710 | -0.0692 | 0.1212 | -0.0085 | -0.5707 | 0.8847 | 0.6680 | 0.0075 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(nadig_jurkat)` | 859 | +0.0039 | [+0.0013, +0.0062] | +0.0498 | 0.0025 | [0.0006, 0.0054] | +0.0483 | +0.0358 | +0.0349 | 0.1122 | +0.0040 | +0.3112 | 0.8835 | 0.6777 | 0.0074 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(nourreddine_kolf21j)` | 859 | -0.0046 | [-0.0079, -0.0015] | -0.0049 | 0.0000 | [0.0000, 0.0006] | -0.0048 | -0.0287 | -0.0279 | 0.1646 | -0.0047 | -0.1696 | 0.8923 | 0.6668 | 0.0076 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 859 | +0.0017 | [-0.0002, +0.0035] | +0.0192 | 0.0004 | [0.0000, 0.0016] | +0.0186 | +0.0184 | +0.0179 | 0.0934 | +0.0017 | +0.1917 | 0.8895 | 0.6752 | 0.0074 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(vcc2025_h1)` | 859 | +0.0051 | [+0.0005, +0.0099] | +0.0308 | 0.0009 | [0.0002, 0.0023] | +0.0298 | +0.0190 | +0.0185 | 0.2770 | +0.0056 | +0.0669 | 0.8780 | 0.6584 | 0.0081 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_ShuffledGenes` | 859 | -0.0087 | [-0.0125, -0.0052] | -0.0176 | 0.0003 | [0.0000, 0.0013] | -0.0171 | -0.0513 | -0.0499 | 0.1747 | -0.0090 | -0.2859 | 0.8933 | 0.6660 | 0.0078 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_DiagCov` | 859 | -0.0123 | [-0.0160, -0.0087] | -0.0283 | 0.0008 | [0.0001, 0.0021] | -0.0274 | -0.0783 | -0.0762 | 0.1615 | -0.0126 | -0.4719 | 0.8903 | 0.6670 | 0.0076 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_SrcCorrTgtVar` | 859 | +0.0021 | [+0.0003, +0.0037] | +0.0354 | 0.0013 | [0.0002, 0.0030] | +0.0343 | +0.0236 | +0.0230 | 0.0902 | +0.0021 | +0.2548 | 0.8859 | 0.6768 | 0.0074 |
| Secondary | replogle_k562_essential | `Transfer_Mean_VarScaled` | 859 | +0.0605 | [+0.0510, +0.0709] | +0.1821 | 0.0332 | [0.0257, 0.0417] | +0.1765 | +0.2450 | +0.2385 | 0.2538 | +0.0644 | +0.9395 | 0.8885 | 0.6919 | 0.0069 |
| Secondary | replogle_k562_essential | `Transfer_Mean_GlobalScale` | 859 | +0.0029 | [+0.0024, +0.0034] | +nan | nan | [nan, nan] | +nan | +0.2849 | +0.2774 | 0.0104 | +0.0029 | +26.6524 | 0.8868 | 0.6788 | 0.0072 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(nadig_hepg2)` | 859 | -0.0443 | [-0.0559, -0.0318] | -0.1780 | 0.0317 | [0.0233, 0.0407] | -0.1726 | -0.1727 | -0.1682 | 0.2635 | -0.0475 | -0.6382 | 0.8838 | 0.6487 | 0.0073 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(nadig_jurkat)` | 859 | -0.0147 | [-0.0208, -0.0083] | +0.1090 | 0.0119 | [0.0089, 0.0153] | +0.1056 | -0.0672 | -0.0655 | 0.2253 | -0.0155 | -0.2905 | 0.8873 | 0.6829 | 0.0087 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(nourreddine_kolf21j)` | 859 | +0.0194 | [+0.0090, +0.0291] | +0.1580 | 0.0250 | [0.0191, 0.0315] | +0.1532 | +0.0654 | +0.0636 | 0.3044 | +0.0213 | +0.2090 | 0.8903 | 0.6885 | 0.0087 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 859 | -0.0273 | [-0.0393, -0.0139] | -0.0816 | 0.0067 | [0.0029, 0.0112] | -0.0791 | -0.1224 | -0.1191 | 0.2293 | -0.0288 | -0.5196 | 0.8814 | 0.6613 | 0.0083 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(vcc2025_h1)` | 859 | +0.1441 | [+0.1173, +0.1729] | -0.0118 | 0.0001 | [0.0000, 0.0012] | -0.0115 | +0.2278 | +0.2218 | 0.6495 | +0.2353 | +0.3415 | 0.8908 | 0.6369 | 0.0044 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_ShuffledVar` | 859 | +0.0995 | [+0.0858, +0.1133] | +0.0613 | 0.0038 | [0.0023, 0.0056] | +0.0594 | +0.2382 | +0.2319 | 0.4289 | +0.1201 | +0.5407 | 0.8921 | 0.6574 | 0.0064 |
| Secondary | replogle_k562_essential | `Hybrid_on_TMVS` | 859 | +0.0622 | [+0.0521, +0.0732] | +0.1757 | 0.0309 | [0.0231, 0.0398] | +0.1703 | +0.2203 | +0.2145 | 0.2898 | +0.0676 | +0.7401 | 0.8821 | 0.6889 | 0.0071 |
| Secondary | replogle_k562_essential | `FCT` | 859 | +0.0419 | [+0.0355, +0.0489] | +0.1938 | 0.0376 | [0.0306, 0.0462] | +0.1879 | +0.1468 | +0.1429 | 0.2933 | +0.0458 | +0.4871 | 0.8763 | 0.6883 | 0.0076 |
| Secondary | replogle_rpe1 | `Source_Grand_Mean` | 1070 | -0.0729 | [-0.0905, -0.0547] | +0.5224 | 0.2729 | [0.2516, 0.2938] | +0.5104 | -0.1541 | -0.1508 | 0.4837 | -0.0944 | -0.3117 | 0.5005 | 0.6129 | 0.0106 |
| Secondary | replogle_rpe1 | `Transfer_Mean` | 1070 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.9301 | 0.6109 | 0.0105 |
| Secondary | replogle_rpe1 | `CIPHER_Linear` | 1070 | -0.1486 | [-0.1676, -0.1291] | +0.3427 | 0.1174 | [0.0962, 0.1406] | +0.3349 | -0.2751 | -0.2691 | 0.5523 | -0.2049 | -0.4873 | 0.5120 | 0.3076 | 0.0123 |
| Secondary | replogle_rpe1 | `Hybrid_Raw` | 1070 | +0.0239 | [+0.0204, +0.0276] | +0.2602 | 0.0677 | [0.0499, 0.0878] | +0.2542 | +0.2349 | +0.2299 | 0.1042 | +0.0242 | +2.2062 | 0.9200 | 0.6376 | 0.0102 |
| Secondary | replogle_rpe1 | `Hybrid_ProjShared` | 1070 | +0.0140 | [+0.0120, +0.0160] | +0.2037 | 0.0415 | [0.0319, 0.0523] | +0.1990 | +0.1816 | +0.1777 | 0.0788 | +0.0141 | +2.2555 | 0.9253 | 0.6264 | 0.0103 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(nadig_hepg2)` | 1070 | +0.0084 | [+0.0068, +0.0101] | +0.1624 | 0.0264 | [0.0190, 0.0355] | +0.1587 | +0.1351 | +0.1321 | 0.0637 | +0.0085 | +2.0745 | 0.9297 | 0.6206 | 0.0105 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(nadig_jurkat)` | 1070 | -0.0054 | [-0.0065, -0.0043] | -0.0703 | 0.0049 | [0.0022, 0.0086] | -0.0687 | -0.1164 | -0.1138 | 0.0471 | -0.0054 | -2.4189 | 0.9302 | 0.6008 | 0.0106 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(nourreddine_kolf21j)` | 1070 | -0.0036 | [-0.0051, -0.0021] | -0.0748 | 0.0056 | [0.0030, 0.0089] | -0.0731 | -0.0532 | -0.0520 | 0.0684 | -0.0036 | -0.7610 | 0.9299 | 0.5966 | 0.0107 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 1070 | -0.0023 | [-0.0034, -0.0012] | -0.0589 | 0.0035 | [0.0012, 0.0065] | -0.0576 | -0.0498 | -0.0488 | 0.0468 | -0.0023 | -1.0431 | 0.9283 | 0.6036 | 0.0106 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(vcc2025_h1)` | 1070 | -0.0079 | [-0.0097, -0.0061] | -0.0901 | 0.0081 | [0.0051, 0.0115] | -0.0881 | -0.0823 | -0.0805 | 0.0983 | -0.0080 | -0.8190 | 0.9263 | 0.5832 | 0.0109 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_ShuffledGenes` | 1070 | -0.0058 | [-0.0076, -0.0042] | -0.0995 | 0.0099 | [0.0064, 0.0143] | -0.0972 | -0.0754 | -0.0738 | 0.0784 | -0.0058 | -0.9417 | 0.9291 | 0.5887 | 0.0109 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_DiagCov` | 1070 | +0.0029 | [+0.0013, +0.0043] | -0.0236 | 0.0006 | [0.0000, 0.0021] | -0.0230 | +0.0425 | +0.0416 | 0.0691 | +0.0029 | +0.6019 | 0.9256 | 0.6044 | 0.0107 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 1070 | -0.0003 | [-0.0012, +0.0006] | +0.0185 | 0.0003 | [0.0000, 0.0018] | +0.0181 | -0.0070 | -0.0069 | 0.0417 | -0.0003 | -0.1655 | 0.9293 | 0.6097 | 0.0105 |
| Secondary | replogle_rpe1 | `Transfer_Mean_VarScaled` | 1070 | +0.0323 | [+0.0284, +0.0360] | +0.0920 | 0.0085 | [0.0053, 0.0123] | +0.0899 | +0.2294 | +0.2244 | 0.1440 | +0.0330 | +1.5585 | 0.9254 | 0.6139 | 0.0103 |
| Secondary | replogle_rpe1 | `Transfer_Mean_GlobalScale` | 1070 | +0.0084 | [+0.0077, +0.0092] | +nan | nan | [nan, nan] | +nan | +0.3598 | +0.3521 | 0.0240 | +0.0084 | +14.6828 | 0.9301 | 0.6109 | 0.0106 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(nadig_hepg2)` | 1070 | +0.0255 | [+0.0224, +0.0283] | +0.0264 | 0.0007 | [0.0001, 0.0018] | +0.0258 | +0.1851 | +0.1811 | 0.1409 | +0.0260 | +1.2852 | 0.9267 | 0.6101 | 0.0103 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(nadig_jurkat)` | 1070 | +0.0204 | [+0.0180, +0.0229] | -0.0451 | 0.0020 | [0.0011, 0.0033] | -0.0440 | +0.1929 | +0.1887 | 0.1084 | +0.0207 | +1.7414 | 0.9293 | 0.5974 | 0.0108 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(nourreddine_kolf21j)` | 1070 | +0.0103 | [+0.0081, +0.0125] | +0.0086 | 0.0001 | [0.0000, 0.0005] | +0.0084 | +0.0833 | +0.0815 | 0.1260 | +0.0104 | +0.6472 | 0.9318 | 0.5937 | 0.0110 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 1070 | -0.0124 | [-0.0145, -0.0104] | -0.0273 | 0.0007 | [0.0002, 0.0015] | -0.0267 | -0.1345 | -0.1316 | 0.0944 | -0.0125 | -1.3947 | 0.9295 | 0.5936 | 0.0106 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(vcc2025_h1)` | 1070 | -0.0823 | [-0.0903, -0.0745] | +0.0222 | 0.0005 | [0.0000, 0.0014] | +0.0217 | -0.3083 | -0.3016 | 0.2728 | -0.0882 | -1.1056 | 0.9313 | 0.5925 | 0.0109 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_ShuffledVar` | 1070 | -0.0350 | [-0.0385, -0.0314] | +0.0259 | 0.0007 | [0.0002, 0.0013] | +0.0253 | -0.1885 | -0.1844 | 0.1896 | -0.0362 | -0.9727 | 0.9344 | 0.5723 | 0.0112 |
| Secondary | replogle_rpe1 | `Hybrid_on_TMVS` | 1070 | +0.0463 | [+0.0418, +0.0506] | +0.1797 | 0.0323 | [0.0253, 0.0401] | +0.1756 | +0.2733 | +0.2674 | 0.1732 | +0.0476 | +1.5438 | 0.9213 | 0.6248 | 0.0102 |
| Secondary | replogle_rpe1 | `FCT` | 1070 | +0.0228 | [+0.0202, +0.0253] | +0.1296 | 0.0168 | [0.0124, 0.0217] | +0.1267 | +0.1592 | +0.1558 | 0.1463 | +0.0233 | +1.0647 | 0.9308 | 0.6151 | 0.0103 |
| Confirmatory | vcc2025_h1 | `Source_Grand_Mean` | 80 | +0.1897 | [-0.0224, +0.3795] | +0.1087 | 0.0118 | [0.0022, 0.0303] | +0.1084 | +0.2991 | +0.2985 | 0.6354 | +0.2911 | +0.4698 | 0.5003 | 0.2596 | 0.0048 |
| Confirmatory | vcc2025_h1 | `Transfer_Mean` | 80 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.8932 | 0.4900 | 0.0054 |
| Confirmatory | vcc2025_h1 | `CIPHER_Linear` | 80 | +0.1589 | [-0.0636, +0.3568] | -0.0468 | 0.0022 | [0.0000, 0.0189] | -0.0467 | +0.2317 | +0.2312 | 0.6870 | +0.2777 | +0.3366 | 0.4725 | 0.0066 | 0.0045 |
| Confirmatory | vcc2025_h1 | `Hybrid_Raw` | 80 | -0.0080 | [-0.0240, +0.0053] | -0.0384 | 0.0015 | [0.0000, 0.0144] | -0.0383 | -0.0429 | -0.0428 | 0.1878 | -0.0083 | -0.2279 | 0.8579 | 0.4385 | 0.0063 |
| Confirmatory | vcc2025_h1 | `Hybrid_ProjShared` | 80 | -0.0064 | [-0.0216, +0.0063] | -0.0299 | 0.0009 | [0.0000, 0.0125] | -0.0298 | -0.0355 | -0.0355 | 0.1811 | -0.0066 | -0.1959 | 0.8468 | 0.4420 | 0.0063 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_MismatchedCov(nadig_hepg2)` | 80 | -0.0026 | [-0.0072, +0.0012] | -0.0471 | 0.0022 | [0.0001, 0.0080] | -0.0469 | -0.0297 | -0.0296 | 0.0873 | -0.0026 | -0.3392 | 0.8813 | 0.4699 | 0.0057 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_MismatchedCov(nadig_jurkat)` | 80 | -0.0003 | [-0.0048, +0.0030] | -0.0036 | 0.0000 | [0.0000, 0.0046] | -0.0036 | -0.0042 | -0.0042 | 0.0663 | -0.0003 | -0.0629 | 0.8753 | 0.4722 | 0.0056 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_MismatchedCov(nourreddine_kolf21j)` | 80 | -0.0015 | [-0.0069, +0.0045] | +0.0006 | 0.0000 | [0.0000, 0.0024] | +0.0006 | -0.0184 | -0.0184 | 0.0812 | -0.0015 | -0.2265 | 0.8813 | 0.4775 | 0.0057 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 80 | +0.0016 | [-0.0015, +0.0058] | +0.0459 | 0.0021 | [0.0000, 0.0123] | +0.0458 | +0.0243 | +0.0242 | 0.0645 | +0.0016 | +0.3755 | 0.8756 | 0.4823 | 0.0056 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 80 | -0.0003 | [-0.0044, +0.0047] | -0.0159 | 0.0003 | [0.0000, 0.0060] | -0.0159 | -0.0053 | -0.0053 | 0.0630 | -0.0003 | -0.0845 | 0.8894 | 0.4826 | 0.0056 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_ShuffledGenes` | 80 | +0.0003 | [-0.0058, +0.0072] | +0.0387 | 0.0015 | [0.0000, 0.0070] | +0.0386 | +0.0021 | +0.0021 | 0.1205 | +0.0003 | +0.0173 | 0.8862 | 0.4723 | 0.0064 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_DiagCov` | 80 | +0.0004 | [-0.0051, +0.0076] | +0.0583 | 0.0034 | [0.0000, 0.0151] | +0.0581 | +0.0054 | +0.0054 | 0.0823 | +0.0004 | +0.0651 | 0.8832 | 0.4872 | 0.0057 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 80 | +0.0030 | [-0.0004, +0.0063] | +0.0126 | 0.0002 | [0.0000, 0.0039] | +0.0126 | +0.0498 | +0.0496 | 0.0601 | +0.0030 | +0.8256 | 0.8790 | 0.4825 | 0.0055 |
| Confirmatory | vcc2025_h1 | `Transfer_Mean_VarScaled` | 80 | -0.1647 | [-0.2637, -0.0561] | -0.0772 | 0.0060 | [0.0000, 0.0378] | -0.0770 | -0.4447 | -0.4437 | 0.3711 | -0.1845 | -1.1958 | 0.8900 | 0.4810 | 0.0042 |
| Confirmatory | vcc2025_h1 | `Transfer_Mean_GlobalScale` | 80 | +0.0683 | [-0.0229, +0.1504] | +nan | nan | [nan, nan] | +nan | +0.2506 | +0.2501 | 0.2730 | +0.0734 | +0.9160 | 0.8932 | 0.4900 | 0.0039 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_MismatchedVar(nadig_hepg2)` | 80 | +0.0353 | [+0.0006, +0.0659] | +0.0629 | 0.0040 | [0.0000, 0.0169] | +0.0627 | +0.3153 | +0.3147 | 0.1122 | +0.0357 | +2.8041 | 0.8905 | 0.4878 | 0.0050 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_MismatchedVar(nadig_jurkat)` | 80 | -0.0416 | [-0.0748, -0.0049] | -0.0690 | 0.0048 | [0.0006, 0.0141] | -0.0688 | -0.3536 | -0.3529 | 0.1180 | -0.0421 | -2.9913 | 0.8908 | 0.4831 | 0.0058 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_MismatchedVar(nourreddine_kolf21j)` | 80 | -0.0471 | [-0.0918, +0.0038] | -0.0074 | 0.0001 | [0.0000, 0.0041] | -0.0074 | -0.2907 | -0.2901 | 0.1624 | -0.0483 | -1.7858 | 0.8927 | 0.4864 | 0.0060 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 80 | +0.0363 | [+0.0015, +0.0668] | +0.0601 | 0.0036 | [0.0003, 0.0106] | +0.0599 | +0.3323 | +0.3316 | 0.1094 | +0.0367 | +3.0317 | 0.8934 | 0.4895 | 0.0049 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 80 | +0.0384 | [-0.0024, +0.0756] | -0.0222 | 0.0005 | [0.0000, 0.0082] | -0.0221 | +0.3026 | +0.3020 | 0.1273 | +0.0390 | +2.3728 | 0.8924 | 0.4801 | 0.0054 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_ShuffledVar` | 80 | +0.0602 | [-0.0300, +0.1412] | +0.0342 | 0.0012 | [0.0000, 0.0074] | +0.0341 | +0.1891 | +0.1887 | 0.3188 | +0.0667 | +0.5921 | 0.8930 | 0.4876 | 0.0038 |
| Confirmatory | vcc2025_h1 | `Hybrid_on_TMVS` | 80 | -0.1711 | [-0.2709, -0.0606] | -0.0746 | 0.0056 | [0.0000, 0.0319] | -0.0743 | -0.4109 | -0.4100 | 0.4173 | -0.1990 | -0.9825 | 0.8128 | 0.3967 | 0.0051 |
| Confirmatory | vcc2025_h1 | `FCT` | 80 | +0.2503 | [+0.1339, +0.3542] | +0.3505 | 0.1228 | [0.0775, 0.1802] | +0.3495 | +0.5489 | +0.5478 | 0.4569 | +0.2898 | +1.1988 | 0.9134 | 0.5552 | 0.0069 |

## R18: TMVS Shuffled-Variance 100-Permutation Null vs H-B3 (97.5% CI Threshold)

| Target Line | TMVS r_orth_cal | TMVS CI lo | Perm100 p2.5 | Perm100 median | Perm100 p97.5 | Perm100 max | Max Mismatch CI hi | H-B3 Pass |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| nadig_hepg2 | +0.0419 | +0.0247 | +0.0122 | +0.0538 | +0.0872 | +0.1120 | +0.1140 | False |
| nadig_jurkat | +0.0947 | +0.0800 | +0.0575 | +0.1109 | +0.1528 | +0.1704 | +0.1730 | False |
| nourreddine_kolf21j | +0.1385 | -0.0508 | +0.0234 | +0.0809 | +0.1456 | +0.1760 | +0.1760 | False |
| replogle_k562_essential | +0.1821 | +0.1605 | +0.0198 | +0.0771 | +0.1247 | +0.1349 | +0.1774 | False |
| replogle_rpe1 | +0.0920 | +0.0725 | -0.0189 | +0.0176 | +0.0479 | +0.0560 | +0.0423 | True |
| vcc2025_h1 | -0.0772 | -0.1943 | -0.0224 | +0.0205 | +0.0534 | +0.0818 | +0.1299 | False |

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
| nadig_hepg2 | 799 | 0.8745 | 0.9659 | 0.8784 | 0.8784 | 0.4747 | 0.5524 | +0.1824 | +0.0451 |
| nadig_jurkat | 934 | 0.8465 | 0.9575 | 0.8562 | 0.8562 | 0.3765 | 0.5096 | +0.1220 | +0.0676 |
| nourreddine_kolf21j | 18 | 0.8195 | 0.9491 | 0.8436 | 0.8438 | 0.1627 | 0.4900 | +0.0983 | -0.0446 |
| replogle_k562_essential | 859 | 0.8861 | 0.9693 | 0.9008 | 0.9008 | 0.4547 | 0.6830 | +0.1387 | +0.0267 |
| replogle_rpe1 | 1070 | 0.9134 | 0.9771 | 0.9178 | 0.9178 | 0.6139 | 0.6386 | +0.1761 | +0.1680 |
| vcc2025_h1 | 80 | 0.9885 | 0.9971 | 0.9915 | 0.9915 | 0.8886 | 0.9329 | +0.0323 | +0.0065 |

## S5: Expressed Genes Sensitivity (Drop Bottom 25% Target Control Variance Genes)

| Target Line | Genes Kept | rho_AB_orth (S5) | Method | Slope OLS (S5) | r_orth_cal (S5) | r_orth_cal^2 (S5) | Joint 97.5% CI | r_cal raw (S5) |
|---|---:|---:|---|---:|---:|---:|---|---:|
| nadig_hepg2 | 4611/6148 | 0.8744 | `Hybrid_ProjShared` | +0.0089 | +0.0467 | 0.0022 | [+0.0211, +0.0723] | +0.0655 |
| nadig_hepg2 | 4611/6148 | 0.8744 | `Transfer_Mean_VarScaled` | +0.0359 | +0.0424 | 0.0018 | [+0.0251, +0.0610] | +0.1230 |
| nadig_hepg2 | 4611/6148 | 0.8744 | `Transfer_Mean_GlobalScale` | -0.0084 | +nan | nan | [nan, nan] | -0.1295 |
| nadig_hepg2 | 4611/6148 | 0.8744 | `Hybrid_on_TMVS` | +0.0448 | +0.0548 | 0.0030 | [+0.0349, +0.0763] | +0.1281 |
| nadig_hepg2 | 4611/6148 | 0.8744 | `FCT` | +0.0591 | +0.2624 | 0.0688 | [+0.2432, +0.2820] | +0.2637 |
| nadig_hepg2 | 4611/6148 | 0.8744 | `NegCtrl_HybridProj_SrcCorrTgtVar` | +0.0052 | +0.0604 | 0.0037 | [+0.0462, +0.0752] | +0.0769 |
| nadig_hepg2 | 4611/6148 | 0.8744 | `NegCtrl_TMVS_ShuffledVar` | -0.0193 | +0.0589 | 0.0035 | [+0.0479, +0.0698] | -0.0612 |
| nadig_jurkat | 4611/6148 | 0.8465 | `Hybrid_ProjShared` | +0.0102 | +0.1046 | 0.0109 | [+0.0789, +0.1301] | +0.0951 |
| nadig_jurkat | 4611/6148 | 0.8465 | `Transfer_Mean_VarScaled` | +0.0093 | +0.0966 | 0.0093 | [+0.0825, +0.1112] | +0.0430 |
| nadig_jurkat | 4611/6148 | 0.8465 | `Transfer_Mean_GlobalScale` | -0.0013 | +nan | nan | [nan, nan] | -0.0094 |
| nadig_jurkat | 4611/6148 | 0.8465 | `Hybrid_on_TMVS` | +0.0195 | +0.1341 | 0.0180 | [+0.1149, +0.1528] | +0.0797 |
| nadig_jurkat | 4611/6148 | 0.8465 | `FCT` | +0.0198 | +0.1043 | 0.0109 | [+0.0845, +0.1239] | +0.0914 |
| nadig_jurkat | 4611/6148 | 0.8465 | `NegCtrl_HybridProj_SrcCorrTgtVar` | +0.0026 | +0.0429 | 0.0018 | [+0.0257, +0.0601] | +0.0396 |
| nadig_jurkat | 4611/6148 | 0.8465 | `NegCtrl_TMVS_ShuffledVar` | +0.0438 | +0.1285 | 0.0165 | [+0.1161, +0.1410] | +0.1490 |
| nourreddine_kolf21j | 4611/6148 | 0.8198 | `Hybrid_ProjShared` | -0.0163 | -0.0116 | 0.0001 | [-0.1515, +0.1118] | -0.0704 |
| nourreddine_kolf21j | 4611/6148 | 0.8198 | `Transfer_Mean_VarScaled` | -0.0133 | +0.1378 | 0.0190 | [-0.0542, +0.2869] | -0.0390 |
| nourreddine_kolf21j | 4611/6148 | 0.8198 | `Transfer_Mean_GlobalScale` | -0.0684 | +nan | nan | [nan, nan] | -0.3297 |
| nourreddine_kolf21j | 4611/6148 | 0.8198 | `Hybrid_on_TMVS` | -0.0296 | +0.1016 | 0.0103 | [-0.0626, +0.2578] | -0.0696 |
| nourreddine_kolf21j | 4611/6148 | 0.8198 | `FCT` | +0.0744 | +0.1835 | 0.0337 | [+0.0212, +0.3537] | +0.2538 |
| nourreddine_kolf21j | 4611/6148 | 0.8198 | `NegCtrl_HybridProj_SrcCorrTgtVar` | -0.0016 | +0.0977 | 0.0096 | [-0.0894, +0.2553] | -0.0084 |
| nourreddine_kolf21j | 4611/6148 | 0.8198 | `NegCtrl_TMVS_ShuffledVar` | +0.0337 | +0.0386 | 0.0015 | [-0.0081, +0.1065] | +0.1187 |
| replogle_k562_essential | 4611/6148 | 0.8847 | `Hybrid_ProjShared` | +0.0018 | +0.0439 | 0.0019 | [+0.0163, +0.0691] | +0.0130 |
| replogle_k562_essential | 4611/6148 | 0.8847 | `Transfer_Mean_VarScaled` | +0.0562 | +0.1868 | 0.0349 | [+0.1654, +0.2091] | +0.2383 |
| replogle_k562_essential | 4611/6148 | 0.8847 | `Transfer_Mean_GlobalScale` | +0.0028 | +nan | nan | [nan, nan] | +0.2804 |
| replogle_k562_essential | 4611/6148 | 0.8847 | `Hybrid_on_TMVS` | +0.0580 | +0.1793 | 0.0321 | [+0.1560, +0.2029] | +0.2125 |
| replogle_k562_essential | 4611/6148 | 0.8847 | `FCT` | +0.0424 | +0.1966 | 0.0386 | [+0.1778, +0.2170] | +0.1482 |
| replogle_k562_essential | 4611/6148 | 0.8847 | `NegCtrl_HybridProj_SrcCorrTgtVar` | +0.0022 | +0.0360 | 0.0013 | [+0.0154, +0.0545] | +0.0253 |
| replogle_k562_essential | 4611/6148 | 0.8847 | `NegCtrl_TMVS_ShuffledVar` | +0.1035 | +0.0656 | 0.0043 | [+0.0545, +0.0767] | +0.2483 |
| replogle_rpe1 | 4611/6148 | 0.9128 | `Hybrid_ProjShared` | +0.0143 | +0.2051 | 0.0421 | [+0.1799, +0.2303] | +0.1838 |
| replogle_rpe1 | 4611/6148 | 0.9128 | `Transfer_Mean_VarScaled` | +0.0351 | +0.0903 | 0.0082 | [+0.0724, +0.1077] | +0.2710 |
| replogle_rpe1 | 4611/6148 | 0.9128 | `Transfer_Mean_GlobalScale` | +0.0083 | +nan | nan | [nan, nan] | +0.3631 |
| replogle_rpe1 | 4611/6148 | 0.9128 | `Hybrid_on_TMVS` | +0.0494 | +0.1828 | 0.0334 | [+0.1622, +0.2036] | +0.3079 |
| replogle_rpe1 | 4611/6148 | 0.9128 | `FCT` | +0.0253 | +0.1318 | 0.0174 | [+0.1144, +0.1498] | +0.1815 |
| replogle_rpe1 | 4611/6148 | 0.9128 | `NegCtrl_HybridProj_SrcCorrTgtVar` | -0.0003 | +0.0188 | 0.0004 | [-0.0047, +0.0433] | -0.0069 |
| replogle_rpe1 | 4611/6148 | 0.9128 | `NegCtrl_TMVS_ShuffledVar` | -0.0356 | +0.0274 | 0.0007 | [+0.0178, +0.0373] | -0.1940 |
| vcc2025_h1 | 4611/6148 | 0.9885 | `Hybrid_ProjShared` | -0.0071 | -0.0295 | 0.0009 | [-0.1132, +0.0444] | -0.0394 |
| vcc2025_h1 | 4611/6148 | 0.9885 | `Transfer_Mean_VarScaled` | -0.1641 | -0.0791 | 0.0063 | [-0.1943, +0.0276] | -0.4487 |
| vcc2025_h1 | 4611/6148 | 0.9885 | `Transfer_Mean_GlobalScale` | +0.0688 | +nan | nan | [nan, nan] | +0.2556 |
| vcc2025_h1 | 4611/6148 | 0.9885 | `Hybrid_on_TMVS` | -0.1712 | -0.0755 | 0.0057 | [-0.1781, +0.0190] | -0.4156 |
| vcc2025_h1 | 4611/6148 | 0.9885 | `FCT` | +0.2509 | +0.3495 | 0.1221 | [+0.2751, +0.4252] | +0.5497 |
| vcc2025_h1 | 4611/6148 | 0.9885 | `NegCtrl_HybridProj_SrcCorrTgtVar` | +0.0028 | +0.0083 | 0.0001 | [-0.0368, +0.0599] | +0.0473 |
| vcc2025_h1 | 4611/6148 | 0.9885 | `NegCtrl_TMVS_ShuffledVar` | +0.0611 | +0.0324 | 0.0010 | [-0.0082, +0.0848] | +0.1935 |

## Appendix: Alternative H-A4/H-B4 Discrimination & Centroid Conventions (R5)

| Target Line | Method | Ctrl disc_raw (A-T2) | Ctrl disc_raw_panel (A-T3) | Cent disc_raw | Cent disc_raw_panel | Cent de50_r | pds_cosine_panel | Slope OLS (All) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| nadig_hepg2 | `Source_Grand_Mean` | 0.5003 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.3059 |
| nadig_hepg2 | `Transfer_Mean` | 0.8811 | 0.8811 | 0.9026 | 0.9028 | 0.5411 | 0.8790 | +0.0000 |
| nadig_hepg2 | `CIPHER_Linear` | 0.4832 | 0.4858 | 0.4893 | 0.4930 | -0.0750 | 0.4848 | +0.3257 |
| nadig_hepg2 | `Hybrid_Raw` | 0.8597 | 0.8606 | 0.8927 | 0.8932 | 0.5530 | 0.8575 | +0.0124 |
| nadig_hepg2 | `Hybrid_ProjShared` | 0.8696 | 0.8702 | 0.8977 | 0.8980 | 0.5449 | 0.8672 | +0.0084 |
| nadig_hepg2 | `NegCtrl_HybridProj_ShuffledGenes` | 0.8816 | 0.8804 | 0.9008 | 0.9019 | 0.5283 | 0.8760 | +0.0041 |
| nadig_hepg2 | `NegCtrl_HybridProj_DiagCov` | 0.8742 | 0.8732 | 0.8990 | 0.8999 | 0.5405 | 0.8738 | +0.0031 |
| nadig_hepg2 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.8810 | 0.8808 | 0.9027 | 0.9029 | 0.5414 | 0.8786 | +0.0023 |
| nadig_hepg2 | `Transfer_Mean_VarScaled` | 0.8745 | 0.8754 | 0.8998 | 0.9007 | 0.5398 | 0.8736 | +0.0329 |
| nadig_hepg2 | `Transfer_Mean_GlobalScale` | 0.8811 | 0.8811 | 0.9039 | 0.9039 | 0.5361 | 0.8790 | +0.0248 |
| nadig_hepg2 | `NegCtrl_TMVS_ShuffledVar` | 0.8902 | 0.8900 | 0.9011 | 0.9029 | 0.4722 | 0.8881 | +0.0626 |
| nadig_hepg2 | `Hybrid_on_TMVS` | 0.8618 | 0.8634 | 0.8939 | 0.8947 | 0.5370 | 0.8609 | +0.0413 |
| nadig_hepg2 | `FCT` | 0.8871 | 0.8883 | 0.9080 | 0.9094 | 0.5682 | 0.8852 | +0.0576 |
| nadig_jurkat | `Source_Grand_Mean` | 0.5002 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.3056 |
| nadig_jurkat | `Transfer_Mean` | 0.8717 | 0.8637 | 0.9218 | 0.9182 | 0.5323 | 0.8530 | +0.0000 |
| nadig_jurkat | `CIPHER_Linear` | 0.4832 | 0.4874 | 0.4936 | 0.4951 | 0.4026 | 0.4947 | +0.4094 |
| nadig_jurkat | `Hybrid_Raw` | 0.8743 | 0.8673 | 0.9228 | 0.9191 | 0.5654 | 0.8560 | +0.0209 |
| nadig_jurkat | `Hybrid_ProjShared` | 0.8679 | 0.8605 | 0.9187 | 0.9154 | 0.5421 | 0.8482 | +0.0089 |
| nadig_jurkat | `NegCtrl_HybridProj_ShuffledGenes` | 0.8717 | 0.8630 | 0.9193 | 0.9155 | 0.5279 | 0.8467 | +0.0031 |
| nadig_jurkat | `NegCtrl_HybridProj_DiagCov` | 0.8660 | 0.8570 | 0.9184 | 0.9147 | 0.5132 | 0.8558 | -0.0036 |
| nadig_jurkat | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.8696 | 0.8616 | 0.9210 | 0.9174 | 0.5335 | 0.8523 | +0.0022 |
| nadig_jurkat | `Transfer_Mean_VarScaled` | 0.8709 | 0.8629 | 0.9135 | 0.9100 | 0.5053 | 0.8514 | -0.0493 |
| nadig_jurkat | `Transfer_Mean_GlobalScale` | 0.8717 | 0.8637 | 0.9141 | 0.9104 | 0.4919 | 0.8530 | -0.0645 |
| nadig_jurkat | `NegCtrl_TMVS_ShuffledVar` | 0.8814 | 0.8756 | 0.9307 | 0.9283 | 0.5930 | 0.8661 | +0.0075 |
| nadig_jurkat | `Hybrid_on_TMVS` | 0.8678 | 0.8604 | 0.9109 | 0.9072 | 0.5148 | 0.8475 | -0.0404 |
| nadig_jurkat | `FCT` | 0.8603 | 0.8519 | 0.9042 | 0.9006 | 0.5139 | 0.8359 | -0.0071 |
| nourreddine_kolf21j | `Source_Grand_Mean` | 0.5000 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.4560 |
| nourreddine_kolf21j | `Transfer_Mean` | 0.8137 | 0.8105 | 0.7778 | 0.7778 | 0.4715 | 0.8072 | +0.0000 |
| nourreddine_kolf21j | `CIPHER_Linear` | 0.5817 | 0.5882 | 0.5196 | 0.5196 | 0.6784 | 0.5098 | +0.6231 |
| nourreddine_kolf21j | `Hybrid_Raw` | 0.7810 | 0.7810 | 0.7516 | 0.7516 | 0.5053 | 0.7451 | +0.0212 |
| nourreddine_kolf21j | `Hybrid_ProjShared` | 0.7941 | 0.7908 | 0.7549 | 0.7582 | 0.4540 | 0.7484 | -0.0000 |
| nourreddine_kolf21j | `NegCtrl_HybridProj_ShuffledGenes` | 0.7484 | 0.7484 | 0.7222 | 0.7255 | 0.4431 | 0.7288 | +0.0036 |
| nourreddine_kolf21j | `NegCtrl_HybridProj_DiagCov` | 0.7549 | 0.7582 | 0.7255 | 0.7222 | 0.4283 | 0.7582 | -0.0034 |
| nourreddine_kolf21j | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.7974 | 0.7974 | 0.7353 | 0.7386 | 0.4673 | 0.7745 | +0.0027 |
| nourreddine_kolf21j | `Transfer_Mean_VarScaled` | 0.8105 | 0.8105 | 0.7647 | 0.7647 | 0.5042 | 0.7941 | -0.1108 |
| nourreddine_kolf21j | `Transfer_Mean_GlobalScale` | 0.8137 | 0.8105 | 0.7745 | 0.7745 | 0.3784 | 0.8072 | -0.1538 |
| nourreddine_kolf21j | `NegCtrl_TMVS_ShuffledVar` | 0.8105 | 0.8105 | 0.7712 | 0.7712 | 0.5578 | 0.8007 | -0.0325 |
| nourreddine_kolf21j | `Hybrid_on_TMVS` | 0.7974 | 0.7941 | 0.7647 | 0.7647 | 0.4928 | 0.7516 | -0.1109 |
| nourreddine_kolf21j | `FCT` | 0.8170 | 0.8170 | 0.7712 | 0.7778 | 0.5354 | 0.7810 | +0.0620 |
| replogle_k562_essential | `Source_Grand_Mean` | 0.5001 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.4825 |
| replogle_k562_essential | `Transfer_Mean` | 0.8868 | 0.8927 | 0.9374 | 0.9390 | 0.4539 | 0.8788 | +0.0000 |
| replogle_k562_essential | `CIPHER_Linear` | 0.4382 | 0.4424 | 0.4878 | 0.4894 | 0.4603 | 0.4612 | +0.6932 |
| replogle_k562_essential | `Hybrid_Raw` | 0.8931 | 0.8974 | 0.9480 | 0.9487 | 0.5117 | 0.8892 | +0.0465 |
| replogle_k562_essential | `Hybrid_ProjShared` | 0.8806 | 0.8850 | 0.9335 | 0.9344 | 0.4439 | 0.8754 | +0.0005 |
| replogle_k562_essential | `NegCtrl_HybridProj_ShuffledGenes` | 0.8933 | 0.8957 | 0.9365 | 0.9366 | 0.4283 | 0.8857 | -0.0030 |
| replogle_k562_essential | `NegCtrl_HybridProj_DiagCov` | 0.8903 | 0.8930 | 0.9342 | 0.9346 | 0.4177 | 0.8896 | -0.0098 |
| replogle_k562_essential | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.8859 | 0.8911 | 0.9360 | 0.9377 | 0.4556 | 0.8796 | +0.0005 |
| replogle_k562_essential | `Transfer_Mean_VarScaled` | 0.8885 | 0.8948 | 0.9425 | 0.9449 | 0.5267 | 0.8801 | +0.0408 |
| replogle_k562_essential | `Transfer_Mean_GlobalScale` | 0.8868 | 0.8927 | 0.9380 | 0.9396 | 0.4569 | 0.8788 | +0.0058 |
| replogle_k562_essential | `NegCtrl_TMVS_ShuffledVar` | 0.8921 | 0.8995 | 0.9541 | 0.9552 | 0.5777 | 0.8880 | +0.1133 |
| replogle_k562_essential | `Hybrid_on_TMVS` | 0.8821 | 0.8870 | 0.9384 | 0.9401 | 0.5137 | 0.8767 | +0.0413 |
| replogle_k562_essential | `FCT` | 0.8763 | 0.8846 | 0.9309 | 0.9368 | 0.4754 | 0.8610 | +0.0440 |
| replogle_rpe1 | `Source_Grand_Mean` | 0.5005 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.1628 |
| replogle_rpe1 | `Transfer_Mean` | 0.9301 | 0.9346 | 0.9173 | 0.9202 | 0.4525 | 0.9287 | +0.0000 |
| replogle_rpe1 | `CIPHER_Linear` | 0.5120 | 0.5157 | 0.5123 | 0.5143 | -0.2142 | 0.5164 | +0.1241 |
| replogle_rpe1 | `Hybrid_Raw` | 0.9200 | 0.9235 | 0.9316 | 0.9341 | 0.5099 | 0.9179 | +0.0225 |
| replogle_rpe1 | `Hybrid_ProjShared` | 0.9253 | 0.9295 | 0.9262 | 0.9296 | 0.4930 | 0.9242 | +0.0156 |
| replogle_rpe1 | `NegCtrl_HybridProj_ShuffledGenes` | 0.9291 | 0.9330 | 0.9092 | 0.9125 | 0.4287 | 0.9259 | -0.0064 |
| replogle_rpe1 | `NegCtrl_HybridProj_DiagCov` | 0.9256 | 0.9294 | 0.9146 | 0.9176 | 0.4472 | 0.9269 | -0.0012 |
| replogle_rpe1 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.9293 | 0.9334 | 0.9166 | 0.9196 | 0.4519 | 0.9275 | -0.0001 |
| replogle_rpe1 | `Transfer_Mean_VarScaled` | 0.9254 | 0.9302 | 0.9238 | 0.9256 | 0.4894 | 0.9251 | +0.0132 |
| replogle_rpe1 | `Transfer_Mean_GlobalScale` | 0.9301 | 0.9346 | 0.9208 | 0.9234 | 0.4611 | 0.9287 | -0.0036 |
| replogle_rpe1 | `NegCtrl_TMVS_ShuffledVar` | 0.9344 | 0.9387 | 0.9018 | 0.9073 | 0.3553 | 0.9327 | -0.0071 |
| replogle_rpe1 | `Hybrid_on_TMVS` | 0.9213 | 0.9260 | 0.9302 | 0.9328 | 0.5197 | 0.9211 | +0.0288 |
| replogle_rpe1 | `FCT` | 0.9308 | 0.9347 | 0.9181 | 0.9220 | 0.4648 | 0.9290 | +0.0358 |
| vcc2025_h1 | `Source_Grand_Mean` | 0.5003 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.3627 |
| vcc2025_h1 | `Transfer_Mean` | 0.8932 | 0.8921 | 0.8354 | 0.8362 | 0.5056 | 0.8758 | +0.0000 |
| vcc2025_h1 | `CIPHER_Linear` | 0.4725 | 0.4718 | 0.4826 | 0.4824 | 0.3693 | 0.4663 | +0.4676 |
| vcc2025_h1 | `Hybrid_Raw` | 0.8579 | 0.8589 | 0.8247 | 0.8266 | 0.4682 | 0.8538 | +0.0000 |
| vcc2025_h1 | `Hybrid_ProjShared` | 0.8468 | 0.8479 | 0.8353 | 0.8348 | 0.4530 | 0.8427 | -0.0061 |
| vcc2025_h1 | `NegCtrl_HybridProj_ShuffledGenes` | 0.8862 | 0.8848 | 0.8399 | 0.8394 | 0.4793 | 0.8772 | -0.0023 |
| vcc2025_h1 | `NegCtrl_HybridProj_DiagCov` | 0.8832 | 0.8842 | 0.8419 | 0.8430 | 0.4910 | 0.8801 | -0.0038 |
| vcc2025_h1 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.8790 | 0.8759 | 0.8339 | 0.8342 | 0.4979 | 0.8661 | +0.0013 |
| vcc2025_h1 | `Transfer_Mean_VarScaled` | 0.8900 | 0.8886 | 0.7252 | 0.7228 | 0.5599 | 0.8750 | +0.1924 |
| vcc2025_h1 | `Transfer_Mean_GlobalScale` | 0.8932 | 0.8921 | 0.7328 | 0.7334 | 0.5563 | 0.8758 | +0.1979 |
| vcc2025_h1 | `NegCtrl_TMVS_ShuffledVar` | 0.8930 | 0.8911 | 0.6888 | 0.6892 | 0.5681 | 0.8771 | +0.2344 |
| vcc2025_h1 | `Hybrid_on_TMVS` | 0.8128 | 0.8127 | 0.7070 | 0.7085 | 0.4816 | 0.8130 | +0.1863 |
| vcc2025_h1 | `FCT` | 0.9134 | 0.9123 | 0.8905 | 0.8905 | 0.6073 | 0.8998 | +0.0427 |
