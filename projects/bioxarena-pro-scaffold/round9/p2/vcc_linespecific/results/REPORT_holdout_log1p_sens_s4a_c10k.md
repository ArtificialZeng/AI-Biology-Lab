# Line-Specific Zero-Shot Perturbation Evaluation Report (`holdout_log1p_sens_s4a_c10k`)

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
| Confirmatory | nadig_hepg2 | `CIPHER_Linear` | 799 | -0.0679 | [-0.1022, -0.0330] | +0.1174 | 0.0138 | [0.0072, 0.0226] | +0.1134 | -0.0841 | -0.0814 | 0.8347 | -0.2137 | -0.0975 | 0.4832 | 0.2594 | 0.0087 |
| Confirmatory | nadig_hepg2 | `Hybrid_Raw` | 799 | +0.0342 | [+0.0278, +0.0410] | +0.1213 | 0.0147 | [0.0076, 0.0242] | +0.1172 | +0.1817 | +0.1757 | 0.1949 | +0.0355 | +0.9016 | 0.8603 | 0.7051 | 0.0081 |
| Confirmatory | nadig_hepg2 | `Hybrid_ProjShared` | 799 | +0.0094 | [+0.0059, +0.0131] | +0.0478 | 0.0023 | [0.0005, 0.0053] | +0.0462 | +0.0710 | +0.0686 | 0.1370 | +0.0096 | +0.5010 | 0.8693 | 0.6984 | 0.0081 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_MismatchedCov(nadig_jurkat)` | 799 | -0.0012 | [-0.0022, -0.0002] | -0.0093 | 0.0001 | [0.0000, 0.0005] | -0.0090 | -0.0172 | -0.0166 | 0.0695 | -0.0012 | -0.2393 | 0.8806 | 0.7006 | 0.0080 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_MismatchedCov(nourreddine_kolf21j)` | 799 | +0.0061 | [+0.0044, +0.0079] | +0.0425 | 0.0018 | [0.0009, 0.0030] | +0.0411 | +0.0574 | +0.0555 | 0.1098 | +0.0062 | +0.5060 | 0.8790 | 0.7003 | 0.0081 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 799 | +0.0026 | [+0.0015, +0.0037] | +0.0242 | 0.0006 | [0.0001, 0.0014] | +0.0234 | +0.0379 | +0.0367 | 0.0715 | +0.0026 | +0.5135 | 0.8797 | 0.7040 | 0.0079 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 799 | -0.0043 | [-0.0053, -0.0033] | -0.0387 | 0.0015 | [0.0006, 0.0028] | -0.0374 | -0.0616 | -0.0596 | 0.0725 | -0.0043 | -0.8214 | 0.8820 | 0.6959 | 0.0080 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_MismatchedCov(vcc2025_h1)` | 799 | +0.0023 | [-0.0001, +0.0047] | +0.0189 | 0.0004 | [0.0000, 0.0010] | +0.0182 | +0.0135 | +0.0131 | 0.1766 | +0.0024 | +0.0742 | 0.8783 | 0.6845 | 0.0086 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_ShuffledGenes` | 799 | +0.0055 | [+0.0033, +0.0079] | +0.0211 | 0.0004 | [0.0001, 0.0012] | +0.0203 | +0.0402 | +0.0389 | 0.1412 | +0.0056 | +0.2753 | 0.8774 | 0.6926 | 0.0086 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_DiagCov` | 799 | +0.0062 | [+0.0043, +0.0084] | +0.0300 | 0.0009 | [0.0003, 0.0019] | +0.0290 | +0.0608 | +0.0588 | 0.1063 | +0.0063 | +0.5531 | 0.8781 | 0.7008 | 0.0080 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 799 | +0.0034 | [+0.0023, +0.0044] | +0.0402 | 0.0016 | [0.0007, 0.0029] | +0.0388 | +0.0516 | +0.0499 | 0.0674 | +0.0034 | +0.7402 | 0.8808 | 0.7055 | 0.0079 |
| Confirmatory | nadig_hepg2 | `Transfer_Mean_VarScaled` | 799 | +0.0357 | [+0.0273, +0.0441] | +0.0421 | 0.0018 | [0.0006, 0.0037] | +0.0406 | +0.1224 | +0.1184 | 0.3018 | +0.0392 | +0.3922 | 0.8745 | 0.6900 | 0.0076 |
| Confirmatory | nadig_hepg2 | `Transfer_Mean_GlobalScale` | 799 | -0.0084 | [-0.0110, -0.0057] | +nan | nan | [nan, nan] | +nan | -0.1306 | -0.1263 | 0.0668 | -0.0085 | -1.8925 | 0.8811 | 0.7053 | 0.0075 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_MismatchedVar(nadig_jurkat)` | 799 | +0.0075 | [+0.0037, +0.0112] | -0.0632 | 0.0040 | [0.0026, 0.0057] | -0.0611 | +0.0497 | +0.0480 | 0.1562 | +0.0077 | +0.3074 | 0.8803 | 0.6929 | 0.0087 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_MismatchedVar(nourreddine_kolf21j)` | 799 | +0.0089 | [+0.0051, +0.0122] | +0.0478 | 0.0023 | [0.0012, 0.0037] | +0.0461 | +0.0499 | +0.0483 | 0.1839 | +0.0092 | +0.2623 | 0.8859 | 0.6990 | 0.0088 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 799 | +0.0087 | [+0.0044, +0.0138] | +0.0963 | 0.0093 | [0.0061, 0.0133] | +0.0930 | +0.0607 | +0.0587 | 0.1475 | +0.0088 | +0.3977 | 0.8831 | 0.7044 | 0.0078 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 799 | +0.0150 | [+0.0071, +0.0227] | -0.0404 | 0.0016 | [0.0003, 0.0040] | -0.0390 | +0.0693 | +0.0670 | 0.2236 | +0.0158 | +0.2997 | 0.8731 | 0.6887 | 0.0083 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_MismatchedVar(vcc2025_h1)` | 799 | -0.0581 | [-0.0727, -0.0431] | -0.0183 | 0.0003 | [0.0000, 0.0018] | -0.0177 | -0.1404 | -0.1358 | 0.4281 | -0.0708 | -0.3171 | 0.8823 | 0.6724 | 0.0070 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_ShuffledVar` | 799 | -0.0196 | [-0.0289, -0.0099] | +0.0591 | 0.0035 | [0.0023, 0.0049] | +0.0571 | -0.0620 | -0.0600 | 0.3269 | -0.0219 | -0.1836 | 0.8902 | 0.6707 | 0.0077 |
| Confirmatory | nadig_hepg2 | `Hybrid_on_TMVS` | 799 | +0.0451 | [+0.0353, +0.0553] | +0.0545 | 0.0030 | [0.0012, 0.0057] | +0.0526 | +0.1285 | +0.1243 | 0.3631 | +0.0518 | +0.3424 | 0.8615 | 0.6803 | 0.0079 |
| Confirmatory | nadig_hepg2 | `FCT` | 799 | +0.0588 | [+0.0536, +0.0640] | +0.2611 | 0.0682 | [0.0585, 0.0788] | +0.2522 | +0.2618 | +0.2532 | 0.2322 | +0.0619 | +1.0905 | 0.8871 | 0.7269 | 0.0076 |
| Confirmatory | nadig_jurkat | `Source_Grand_Mean` | 934 | -0.0036 | [-0.0349, +0.0310] | +0.0153 | 0.0002 | [0.0000, 0.0020] | +0.0147 | -0.0049 | -0.0047 | 0.7616 | -0.0086 | -0.0062 | 0.5002 | 0.3335 | 0.0068 |
| Confirmatory | nadig_jurkat | `Transfer_Mean` | 934 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.8717 | 0.6172 | 0.0076 |
| Confirmatory | nadig_jurkat | `CIPHER_Linear` | 934 | +0.0070 | [-0.0284, +0.0465] | +0.0465 | 0.0022 | [0.0002, 0.0063] | +0.0445 | +0.0080 | +0.0077 | 0.9150 | +0.0432 | +0.0084 | 0.4831 | 0.0742 | 0.0071 |
| Confirmatory | nadig_jurkat | `Hybrid_Raw` | 934 | +0.0201 | [+0.0143, +0.0267] | +0.1298 | 0.0169 | [0.0100, 0.0259] | +0.1243 | +0.1488 | +0.1429 | 0.1407 | +0.0205 | +1.0157 | 0.8744 | 0.6223 | 0.0076 |
| Confirmatory | nadig_jurkat | `Hybrid_ProjShared` | 934 | +0.0060 | [+0.0030, +0.0091] | +0.0658 | 0.0043 | [0.0017, 0.0082] | +0.0630 | +0.0582 | +0.0559 | 0.1066 | +0.0060 | +0.5243 | 0.8685 | 0.6157 | 0.0078 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_MismatchedCov(nadig_hepg2)` | 934 | -0.0054 | [-0.0073, -0.0035] | -0.0494 | 0.0024 | [0.0008, 0.0048] | -0.0473 | -0.0575 | -0.0552 | 0.0982 | -0.0055 | -0.5622 | 0.8678 | 0.6050 | 0.0078 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_MismatchedCov(nourreddine_kolf21j)` | 934 | -0.0003 | [-0.0025, +0.0020] | +0.0037 | 0.0000 | [0.0000, 0.0006] | +0.0036 | -0.0023 | -0.0022 | 0.1153 | -0.0003 | -0.0190 | 0.8726 | 0.6072 | 0.0079 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 934 | +0.0005 | [-0.0012, +0.0021] | +0.0110 | 0.0001 | [0.0000, 0.0008] | +0.0106 | +0.0062 | +0.0059 | 0.0825 | +0.0005 | +0.0719 | 0.8720 | 0.6131 | 0.0078 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 934 | +0.0025 | [+0.0011, +0.0039] | +0.0305 | 0.0009 | [0.0001, 0.0025] | +0.0292 | +0.0363 | +0.0349 | 0.0712 | +0.0025 | +0.4898 | 0.8708 | 0.6152 | 0.0077 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_MismatchedCov(vcc2025_h1)` | 934 | +0.0075 | [+0.0029, +0.0122] | +0.0349 | 0.0012 | [0.0002, 0.0030] | +0.0334 | +0.0397 | +0.0381 | 0.1972 | +0.0078 | +0.1934 | 0.8702 | 0.5981 | 0.0084 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_ShuffledGenes` | 934 | +0.0024 | [-0.0003, +0.0052] | +0.0153 | 0.0002 | [0.0000, 0.0013] | +0.0147 | +0.0199 | +0.0191 | 0.1259 | +0.0024 | +0.1515 | 0.8707 | 0.6078 | 0.0081 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_DiagCov` | 934 | +0.0021 | [-0.0005, +0.0047] | +0.0158 | 0.0003 | [0.0000, 0.0016] | +0.0152 | +0.0201 | +0.0193 | 0.1063 | +0.0021 | +0.1818 | 0.8713 | 0.6111 | 0.0078 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_SrcCorrTgtVar` | 934 | +0.0015 | [+0.0003, +0.0027] | +0.0270 | 0.0007 | [0.0001, 0.0019] | +0.0259 | +0.0233 | +0.0224 | 0.0662 | +0.0015 | +0.3379 | 0.8709 | 0.6148 | 0.0077 |
| Confirmatory | nadig_jurkat | `Transfer_Mean_VarScaled` | 934 | +0.0099 | [+0.0031, +0.0162] | +0.0967 | 0.0093 | [0.0068, 0.0124] | +0.0926 | +0.0456 | +0.0438 | 0.2257 | +0.0104 | +0.1940 | 0.8710 | 0.6189 | 0.0087 |
| Confirmatory | nadig_jurkat | `Transfer_Mean_GlobalScale` | 934 | -0.0008 | [-0.0068, +0.0048] | +nan | nan | [nan, nan] | +nan | -0.0058 | -0.0056 | 0.1432 | -0.0008 | -0.0391 | 0.8717 | 0.6172 | 0.0088 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_MismatchedVar(nadig_hepg2)` | 934 | -0.0370 | [-0.0438, -0.0301] | -0.1781 | 0.0317 | [0.0259, 0.0384] | -0.1705 | -0.1786 | -0.1716 | 0.2159 | -0.0388 | -0.7945 | 0.8671 | 0.5856 | 0.0079 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_MismatchedVar(nourreddine_kolf21j)` | 934 | +0.0418 | [+0.0334, +0.0489] | +0.1589 | 0.0252 | [0.0208, 0.0299] | +0.1521 | +0.1726 | +0.1658 | 0.2523 | +0.0446 | +0.6570 | 0.8798 | 0.6284 | 0.0088 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 934 | +0.0141 | [+0.0075, +0.0212] | +0.1061 | 0.0112 | [0.0076, 0.0153] | +0.1016 | +0.0888 | +0.0853 | 0.1653 | +0.0145 | +0.5157 | 0.8745 | 0.6215 | 0.0076 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 934 | -0.0310 | [-0.0383, -0.0235] | -0.1378 | 0.0190 | [0.0139, 0.0256] | -0.1319 | -0.1683 | -0.1616 | 0.1920 | -0.0322 | -0.8418 | 0.8629 | 0.5926 | 0.0086 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_MismatchedVar(vcc2025_h1)` | 934 | +0.0188 | [+0.0055, +0.0319] | +0.0364 | 0.0013 | [0.0002, 0.0039] | +0.0349 | +0.0424 | +0.0407 | 0.4614 | +0.0238 | +0.0882 | 0.8729 | 0.5807 | 0.0062 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_ShuffledVar` | 934 | +0.0439 | [+0.0378, +0.0502] | +0.1287 | 0.0166 | [0.0135, 0.0201] | +0.1232 | +0.1494 | +0.1435 | 0.3058 | +0.0483 | +0.4693 | 0.8814 | 0.6107 | 0.0082 |
| Confirmatory | nadig_jurkat | `Hybrid_on_TMVS` | 934 | +0.0158 | [+0.0088, +0.0225] | +0.1131 | 0.0128 | [0.0089, 0.0173] | +0.1083 | +0.0651 | +0.0625 | 0.2533 | +0.0169 | +0.2468 | 0.8683 | 0.6180 | 0.0088 |
| Confirmatory | nadig_jurkat | `FCT` | 934 | +0.0200 | [+0.0149, +0.0253] | +0.1043 | 0.0109 | [0.0072, 0.0153] | +0.0999 | +0.0915 | +0.0879 | 0.2274 | +0.0211 | +0.3865 | 0.8603 | 0.6115 | 0.0085 |
| Secondary | nourreddine_kolf21j | `Source_Grand_Mean` | 18 | +0.2723 | [+0.0345, +0.5563] | +0.0672 | 0.0045 | [0.0000, 0.0395] | +0.0638 | +0.3484 | +0.3333 | 0.8168 | +0.5608 | +0.4081 | 0.5000 | 0.2906 | 0.0034 |
| Secondary | nourreddine_kolf21j | `Transfer_Mean` | 18 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.8137 | 0.4335 | 0.0042 |
| Secondary | nourreddine_kolf21j | `CIPHER_Linear` | 18 | +0.2731 | [+0.0230, +0.5569] | -0.1123 | 0.0126 | [0.0000, 0.0698] | -0.1066 | +0.3138 | +0.3002 | 0.9096 | +0.7327 | +0.3300 | 0.5915 | -0.1037 | 0.0023 |
| Secondary | nourreddine_kolf21j | `Hybrid_Raw` | 18 | -0.0137 | [-0.0568, +0.0236] | -0.0054 | 0.0000 | [0.0000, 0.0212] | -0.0051 | -0.0530 | -0.0507 | 0.2691 | -0.0147 | -0.1885 | 0.7843 | 0.4329 | 0.0046 |
| Secondary | nourreddine_kolf21j | `Hybrid_ProjShared` | 18 | -0.0146 | [-0.0591, +0.0248] | -0.0028 | 0.0000 | [0.0000, 0.0242] | -0.0027 | -0.0637 | -0.0609 | 0.2402 | -0.0155 | -0.2537 | 0.8072 | 0.4282 | 0.0046 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_MismatchedCov(nadig_hepg2)` | 18 | -0.0013 | [-0.0496, +0.0432] | +0.0751 | 0.0056 | [0.0000, 0.0400] | +0.0713 | -0.0057 | -0.0054 | 0.2371 | -0.0014 | -0.0229 | 0.8170 | 0.4508 | 0.0046 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_MismatchedCov(nadig_jurkat)` | 18 | -0.0086 | [-0.0579, +0.0337] | +0.0346 | 0.0012 | [0.0000, 0.0287] | +0.0328 | -0.0395 | -0.0378 | 0.2277 | -0.0091 | -0.1660 | 0.7876 | 0.4405 | 0.0046 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 18 | -0.0103 | [-0.0664, +0.0380] | +0.0392 | 0.0015 | [0.0000, 0.0326] | +0.0372 | -0.0446 | -0.0427 | 0.2411 | -0.0109 | -0.1769 | 0.7941 | 0.4461 | 0.0045 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 18 | -0.0178 | [-0.0726, +0.0210] | -0.0070 | 0.0000 | [0.0000, 0.0128] | -0.0067 | -0.0797 | -0.0762 | 0.2333 | -0.0188 | -0.3267 | 0.7908 | 0.4292 | 0.0046 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_MismatchedCov(vcc2025_h1)` | 18 | +0.0121 | [-0.0210, +0.0598] | +0.0721 | 0.0052 | [0.0000, 0.0549] | +0.0684 | +0.0589 | +0.0564 | 0.2141 | +0.0126 | +0.2633 | 0.8235 | 0.4313 | 0.0046 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_ShuffledGenes` | 18 | -0.0107 | [-0.0636, +0.0389] | +0.0418 | 0.0018 | [0.0000, 0.0355] | +0.0397 | -0.0417 | -0.0399 | 0.2670 | -0.0115 | -0.1495 | 0.7712 | 0.4398 | 0.0048 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_DiagCov` | 18 | -0.0092 | [-0.0617, +0.0389] | +0.0503 | 0.0025 | [0.0000, 0.0403] | +0.0477 | -0.0391 | -0.0374 | 0.2458 | -0.0098 | -0.1520 | 0.7712 | 0.4458 | 0.0045 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_SrcCorrTgtVar` | 18 | -0.0026 | [-0.0497, +0.0396] | +0.0895 | 0.0080 | [0.0000, 0.0620] | +0.0849 | -0.0143 | -0.0136 | 0.1940 | -0.0027 | -0.0703 | 0.8072 | 0.4590 | 0.0044 |
| Secondary | nourreddine_kolf21j | `Transfer_Mean_VarScaled` | 18 | -0.0135 | [-0.0925, +0.0527] | +0.1386 | 0.0192 | [0.0000, 0.0821] | +0.1315 | -0.0398 | -0.0380 | 0.3541 | -0.0154 | -0.1074 | 0.8105 | 0.4380 | 0.0055 |
| Secondary | nourreddine_kolf21j | `Transfer_Mean_GlobalScale` | 18 | -0.0683 | [-0.1376, -0.0097] | +nan | nan | [nan, nan] | +nan | -0.3304 | -0.3161 | 0.2161 | -0.0713 | -1.4628 | 0.8137 | 0.4335 | 0.0057 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_MismatchedVar(nadig_hepg2)` | 18 | +0.0173 | [-0.0138, +0.0689] | +0.0179 | 0.0003 | [0.0000, 0.0273] | +0.0170 | +0.1191 | +0.1140 | 0.1514 | +0.0177 | +0.7528 | 0.8203 | 0.4239 | 0.0042 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_MismatchedVar(nadig_jurkat)` | 18 | -0.0412 | [-0.0748, -0.0033] | +0.0467 | 0.0022 | [0.0000, 0.0120] | +0.0443 | -0.2264 | -0.2165 | 0.1901 | -0.0426 | -1.1390 | 0.8039 | 0.4342 | 0.0052 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 18 | +0.0091 | [-0.0027, +0.0226] | +0.1209 | 0.0146 | [0.0030, 0.0312] | +0.1147 | +0.0753 | +0.0721 | 0.1269 | +0.0093 | +0.5678 | 0.8105 | 0.4386 | 0.0043 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 18 | +0.0054 | [-0.0186, +0.0313] | +0.0288 | 0.0008 | [0.0000, 0.0212] | +0.0274 | +0.0349 | +0.0334 | 0.1619 | +0.0055 | +0.2061 | 0.8137 | 0.4309 | 0.0048 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_MismatchedVar(vcc2025_h1)` | 18 | +0.0491 | [+0.0028, +0.0912] | -0.1377 | 0.0189 | [0.0000, 0.0731] | -0.1307 | +0.1339 | +0.1281 | 0.3837 | +0.0574 | +0.3338 | 0.8007 | 0.3959 | 0.0027 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_ShuffledVar` | 18 | +0.0328 | [-0.0147, +0.0874] | +0.0390 | 0.0015 | [0.0000, 0.0115] | +0.0370 | +0.1151 | +0.1101 | 0.2980 | +0.0360 | +0.3695 | 0.8105 | 0.4234 | 0.0049 |
| Secondary | nourreddine_kolf21j | `Hybrid_on_TMVS` | 18 | -0.0281 | [-0.1308, +0.0579] | +0.1071 | 0.0115 | [0.0000, 0.0716] | +0.1017 | -0.0658 | -0.0630 | 0.4464 | -0.0351 | -0.1410 | 0.8039 | 0.4405 | 0.0059 |
| Secondary | nourreddine_kolf21j | `FCT` | 18 | +0.0747 | [+0.0070, +0.1592] | +0.1839 | 0.0338 | [0.0006, 0.1256] | +0.1746 | +0.2547 | +0.2436 | 0.3066 | +0.0819 | +0.7947 | 0.8170 | 0.4414 | 0.0044 |
| Secondary | replogle_k562_essential | `Source_Grand_Mean` | 859 | +0.2588 | [+0.2086, +0.3093] | -0.0643 | 0.0041 | [0.0017, 0.0076] | -0.0623 | +0.2554 | +0.2487 | 1.0408 | +1.1736 | +0.2390 | 0.5001 | 0.3941 | 0.0052 |
| Secondary | replogle_k562_essential | `Transfer_Mean` | 859 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.8868 | 0.6788 | 0.0073 |
| Secondary | replogle_k562_essential | `CIPHER_Linear` | 859 | +0.3500 | [+0.2933, +0.4075] | +0.0267 | 0.0007 | [0.0000, 0.0028] | +0.0259 | +0.2844 | +0.2769 | 1.2641 | +2.1692 | +0.2191 | 0.4395 | 0.0717 | 0.0049 |
| Secondary | replogle_k562_essential | `Hybrid_Raw` | 859 | +0.0309 | [+0.0252, +0.0369] | +0.1055 | 0.0111 | [0.0057, 0.0174] | +0.1022 | +0.1666 | +0.1622 | 0.1905 | +0.0320 | +0.8516 | 0.8930 | 0.6812 | 0.0069 |
| Secondary | replogle_k562_essential | `Hybrid_ProjShared` | 859 | +0.0020 | [-0.0010, +0.0049] | +0.0351 | 0.0012 | [0.0001, 0.0035] | +0.0341 | +0.0162 | +0.0158 | 0.1288 | +0.0021 | +0.1224 | 0.8841 | 0.6765 | 0.0075 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(nadig_hepg2)` | 859 | -0.0092 | [-0.0119, -0.0069] | -0.0914 | 0.0084 | [0.0058, 0.0115] | -0.0886 | -0.0803 | -0.0782 | 0.1178 | -0.0093 | -0.6637 | 0.8854 | 0.6668 | 0.0075 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(nadig_jurkat)` | 859 | +0.0007 | [-0.0016, +0.0028] | +0.0063 | 0.0000 | [0.0000, 0.0007] | +0.0061 | +0.0069 | +0.0068 | 0.1046 | +0.0007 | +0.0646 | 0.8863 | 0.6752 | 0.0074 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(nourreddine_kolf21j)` | 859 | -0.0053 | [-0.0082, -0.0026] | -0.0094 | 0.0001 | [0.0000, 0.0007] | -0.0091 | -0.0351 | -0.0342 | 0.1556 | -0.0054 | -0.2195 | 0.8922 | 0.6671 | 0.0077 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 859 | +0.0048 | [+0.0030, +0.0065] | +0.0715 | 0.0051 | [0.0021, 0.0090] | +0.0693 | +0.0557 | +0.0542 | 0.0879 | +0.0048 | +0.6170 | 0.8865 | 0.6787 | 0.0074 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(vcc2025_h1)` | 859 | +0.0060 | [+0.0017, +0.0103] | +0.0241 | 0.0006 | [0.0001, 0.0017] | +0.0234 | +0.0227 | +0.0221 | 0.2711 | +0.0065 | +0.0814 | 0.8831 | 0.6573 | 0.0081 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_ShuffledGenes` | 859 | -0.0081 | [-0.0115, -0.0049] | -0.0129 | 0.0002 | [0.0000, 0.0009] | -0.0125 | -0.0508 | -0.0495 | 0.1642 | -0.0084 | -0.3015 | 0.8930 | 0.6673 | 0.0078 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_DiagCov` | 859 | -0.0081 | [-0.0114, -0.0050] | -0.0145 | 0.0002 | [0.0000, 0.0011] | -0.0141 | -0.0561 | -0.0546 | 0.1483 | -0.0083 | -0.3680 | 0.8930 | 0.6693 | 0.0076 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_SrcCorrTgtVar` | 859 | +0.0018 | [+0.0001, +0.0034] | +0.0247 | 0.0006 | [0.0000, 0.0019] | +0.0240 | +0.0219 | +0.0214 | 0.0861 | +0.0019 | +0.2479 | 0.8870 | 0.6763 | 0.0074 |
| Secondary | replogle_k562_essential | `Transfer_Mean_VarScaled` | 859 | +0.0611 | [+0.0515, +0.0717] | +0.1833 | 0.0336 | [0.0262, 0.0422] | +0.1777 | +0.2462 | +0.2396 | 0.2551 | +0.0651 | +0.9393 | 0.8885 | 0.6920 | 0.0069 |
| Secondary | replogle_k562_essential | `Transfer_Mean_GlobalScale` | 859 | +0.0028 | [+0.0024, +0.0033] | +nan | nan | [nan, nan] | +nan | +0.2849 | +0.2774 | 0.0101 | +0.0028 | +27.3970 | 0.8868 | 0.6788 | 0.0072 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(nadig_hepg2)` | 859 | -0.0441 | [-0.0557, -0.0316] | -0.1777 | 0.0316 | [0.0231, 0.0406] | -0.1723 | -0.1719 | -0.1673 | 0.2635 | -0.0473 | -0.6351 | 0.8838 | 0.6488 | 0.0073 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(nadig_jurkat)` | 859 | -0.0153 | [-0.0214, -0.0088] | +0.1094 | 0.0120 | [0.0089, 0.0155] | +0.1061 | -0.0695 | -0.0677 | 0.2260 | -0.0161 | -0.2995 | 0.8873 | 0.6830 | 0.0087 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(nourreddine_kolf21j)` | 859 | +0.0199 | [+0.0097, +0.0296] | +0.1589 | 0.0253 | [0.0193, 0.0318] | +0.1541 | +0.0673 | +0.0655 | 0.3043 | +0.0220 | +0.2153 | 0.8903 | 0.6886 | 0.0087 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 859 | -0.0275 | [-0.0395, -0.0141] | -0.0830 | 0.0069 | [0.0031, 0.0115] | -0.0804 | -0.1232 | -0.1200 | 0.2294 | -0.0290 | -0.5231 | 0.8814 | 0.6612 | 0.0083 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(vcc2025_h1)` | 859 | +0.1442 | [+0.1174, +0.1732] | -0.0123 | 0.0002 | [0.0000, 0.0012] | -0.0119 | +0.2277 | +0.2217 | 0.6506 | +0.2361 | +0.3408 | 0.8908 | 0.6367 | 0.0044 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_ShuffledVar` | 859 | +0.0998 | [+0.0861, +0.1136] | +0.0618 | 0.0038 | [0.0023, 0.0056] | +0.0599 | +0.2389 | +0.2326 | 0.4291 | +0.1205 | +0.5420 | 0.8921 | 0.6575 | 0.0064 |
| Secondary | replogle_k562_essential | `Hybrid_on_TMVS` | 859 | +0.0632 | [+0.0530, +0.0747] | +0.1748 | 0.0306 | [0.0230, 0.0392] | +0.1694 | +0.2260 | +0.2200 | 0.2871 | +0.0685 | +0.7662 | 0.8856 | 0.6889 | 0.0071 |
| Secondary | replogle_k562_essential | `FCT` | 859 | +0.0419 | [+0.0355, +0.0489] | +0.1938 | 0.0376 | [0.0306, 0.0462] | +0.1879 | +0.1468 | +0.1429 | 0.2933 | +0.0458 | +0.4871 | 0.8763 | 0.6883 | 0.0076 |
| Secondary | replogle_rpe1 | `Source_Grand_Mean` | 1070 | -0.0729 | [-0.0905, -0.0547] | +0.5224 | 0.2729 | [0.2516, 0.2938] | +0.5104 | -0.1541 | -0.1508 | 0.4837 | -0.0944 | -0.3117 | 0.5005 | 0.6129 | 0.0106 |
| Secondary | replogle_rpe1 | `Transfer_Mean` | 1070 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.9301 | 0.6109 | 0.0105 |
| Secondary | replogle_rpe1 | `CIPHER_Linear` | 1070 | -0.1494 | [-0.1682, -0.1300] | +0.3418 | 0.1169 | [0.0959, 0.1398] | +0.3340 | -0.2767 | -0.2707 | 0.5520 | -0.2058 | -0.4903 | 0.5123 | 0.3073 | 0.0124 |
| Secondary | replogle_rpe1 | `Hybrid_Raw` | 1070 | +0.0236 | [+0.0201, +0.0272] | +0.2592 | 0.0672 | [0.0494, 0.0870] | +0.2532 | +0.2341 | +0.2290 | 0.1029 | +0.0238 | +2.2245 | 0.9200 | 0.6375 | 0.0102 |
| Secondary | replogle_rpe1 | `Hybrid_ProjShared` | 1070 | +0.0130 | [+0.0112, +0.0150] | +0.1950 | 0.0380 | [0.0288, 0.0483] | +0.1905 | +0.1705 | +0.1668 | 0.0781 | +0.0131 | +2.1353 | 0.9249 | 0.6248 | 0.0103 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(nadig_hepg2)` | 1070 | +0.0091 | [+0.0077, +0.0106] | +0.1752 | 0.0307 | [0.0231, 0.0397] | +0.1712 | +0.1550 | +0.1516 | 0.0602 | +0.0092 | +2.5209 | 0.9288 | 0.6208 | 0.0104 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(nadig_jurkat)` | 1070 | -0.0035 | [-0.0045, -0.0026] | -0.0489 | 0.0024 | [0.0007, 0.0049] | -0.0478 | -0.0823 | -0.0805 | 0.0439 | -0.0035 | -1.8334 | 0.9295 | 0.6037 | 0.0106 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(nourreddine_kolf21j)` | 1070 | -0.0072 | [-0.0086, -0.0060] | -0.1155 | 0.0133 | [0.0095, 0.0177] | -0.1128 | -0.1107 | -0.1083 | 0.0666 | -0.0072 | -1.6249 | 0.9276 | 0.5923 | 0.0108 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 1070 | -0.0041 | [-0.0052, -0.0031] | -0.0915 | 0.0084 | [0.0049, 0.0125] | -0.0894 | -0.0946 | -0.0925 | 0.0448 | -0.0042 | -2.0638 | 0.9294 | 0.6016 | 0.0106 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(vcc2025_h1)` | 1070 | -0.0108 | [-0.0125, -0.0091] | -0.1101 | 0.0121 | [0.0088, 0.0157] | -0.1076 | -0.1161 | -0.1136 | 0.0953 | -0.0109 | -1.1920 | 0.9290 | 0.5803 | 0.0109 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_ShuffledGenes` | 1070 | -0.0085 | [-0.0101, -0.0071] | -0.1328 | 0.0176 | [0.0134, 0.0226] | -0.1297 | -0.1155 | -0.1130 | 0.0752 | -0.0085 | -1.5025 | 0.9284 | 0.5861 | 0.0109 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_DiagCov` | 1070 | -0.0072 | [-0.0087, -0.0059] | -0.1359 | 0.0185 | [0.0138, 0.0239] | -0.1328 | -0.1115 | -0.1091 | 0.0663 | -0.0073 | -1.6459 | 0.9277 | 0.5909 | 0.0107 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 1070 | -0.0013 | [-0.0021, -0.0005] | -0.0004 | 0.0000 | [0.0000, 0.0005] | -0.0004 | -0.0346 | -0.0339 | 0.0398 | -0.0013 | -0.8517 | 0.9293 | 0.6085 | 0.0105 |
| Secondary | replogle_rpe1 | `Transfer_Mean_VarScaled` | 1070 | +0.0323 | [+0.0284, +0.0360] | +0.0915 | 0.0084 | [0.0052, 0.0121] | +0.0894 | +0.2294 | +0.2244 | 0.1440 | +0.0330 | +1.5589 | 0.9254 | 0.6139 | 0.0103 |
| Secondary | replogle_rpe1 | `Transfer_Mean_GlobalScale` | 1070 | +0.0085 | [+0.0077, +0.0092] | +nan | nan | [nan, nan] | +nan | +0.3598 | +0.3521 | 0.0241 | +0.0085 | +14.5959 | 0.9301 | 0.6109 | 0.0106 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(nadig_hepg2)` | 1070 | +0.0255 | [+0.0224, +0.0283] | +0.0265 | 0.0007 | [0.0001, 0.0018] | +0.0259 | +0.1847 | +0.1807 | 0.1411 | +0.0260 | +1.2803 | 0.9267 | 0.6101 | 0.0103 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(nadig_jurkat)` | 1070 | +0.0207 | [+0.0183, +0.0232] | -0.0445 | 0.0020 | [0.0011, 0.0032] | -0.0434 | +0.1949 | +0.1907 | 0.1088 | +0.0210 | +1.7531 | 0.9293 | 0.5975 | 0.0108 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(nourreddine_kolf21j)` | 1070 | +0.0102 | [+0.0081, +0.0123] | +0.0085 | 0.0001 | [0.0000, 0.0005] | +0.0083 | +0.0826 | +0.0808 | 0.1259 | +0.0103 | +0.6419 | 0.9318 | 0.5935 | 0.0110 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 1070 | -0.0126 | [-0.0147, -0.0106] | -0.0280 | 0.0008 | [0.0002, 0.0016] | -0.0273 | -0.1353 | -0.1323 | 0.0950 | -0.0127 | -1.3927 | 0.9294 | 0.5935 | 0.0106 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(vcc2025_h1)` | 1070 | -0.0824 | [-0.0904, -0.0747] | +0.0222 | 0.0005 | [0.0000, 0.0014] | +0.0217 | -0.3082 | -0.3015 | 0.2734 | -0.0884 | -1.1031 | 0.9313 | 0.5923 | 0.0109 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_ShuffledVar` | 1070 | -0.0349 | [-0.0385, -0.0314] | +0.0262 | 0.0007 | [0.0002, 0.0014] | +0.0256 | -0.1889 | -0.1848 | 0.1891 | -0.0362 | -0.9769 | 0.9344 | 0.5727 | 0.0112 |
| Secondary | replogle_rpe1 | `Hybrid_on_TMVS` | 1070 | +0.0453 | [+0.0408, +0.0496] | +0.1746 | 0.0305 | [0.0236, 0.0380] | +0.1706 | +0.2681 | +0.2623 | 0.1729 | +0.0466 | +1.5169 | 0.9208 | 0.6236 | 0.0102 |
| Secondary | replogle_rpe1 | `FCT` | 1070 | +0.0228 | [+0.0202, +0.0253] | +0.1296 | 0.0168 | [0.0124, 0.0217] | +0.1267 | +0.1592 | +0.1558 | 0.1463 | +0.0233 | +1.0647 | 0.9308 | 0.6151 | 0.0103 |
| Confirmatory | vcc2025_h1 | `Source_Grand_Mean` | 80 | +0.1897 | [-0.0224, +0.3795] | +0.1087 | 0.0118 | [0.0022, 0.0303] | +0.1084 | +0.2991 | +0.2985 | 0.6354 | +0.2911 | +0.4698 | 0.5003 | 0.2596 | 0.0048 |
| Confirmatory | vcc2025_h1 | `Transfer_Mean` | 80 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.8932 | 0.4900 | 0.0054 |
| Confirmatory | vcc2025_h1 | `CIPHER_Linear` | 80 | +0.1604 | [-0.0631, +0.3596] | -0.0440 | 0.0019 | [0.0000, 0.0184] | -0.0439 | +0.2332 | +0.2327 | 0.6891 | +0.2812 | +0.3377 | 0.4672 | 0.0069 | 0.0045 |
| Confirmatory | vcc2025_h1 | `Hybrid_Raw` | 80 | -0.0065 | [-0.0218, +0.0064] | -0.0352 | 0.0012 | [0.0000, 0.0135] | -0.0351 | -0.0349 | -0.0349 | 0.1863 | -0.0067 | -0.1870 | 0.8595 | 0.4407 | 0.0063 |
| Confirmatory | vcc2025_h1 | `Hybrid_ProjShared` | 80 | -0.0075 | [-0.0228, +0.0053] | -0.0376 | 0.0014 | [0.0000, 0.0143] | -0.0375 | -0.0421 | -0.0421 | 0.1779 | -0.0077 | -0.2365 | 0.8535 | 0.4440 | 0.0063 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_MismatchedCov(nadig_hepg2)` | 80 | -0.0023 | [-0.0070, +0.0020] | -0.0428 | 0.0018 | [0.0000, 0.0076] | -0.0427 | -0.0260 | -0.0260 | 0.0868 | -0.0023 | -0.2991 | 0.8810 | 0.4691 | 0.0057 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_MismatchedCov(nadig_jurkat)` | 80 | +0.0005 | [-0.0029, +0.0035] | +0.0058 | 0.0000 | [0.0000, 0.0042] | +0.0057 | +0.0078 | +0.0078 | 0.0644 | +0.0005 | +0.1211 | 0.8772 | 0.4769 | 0.0056 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_MismatchedCov(nourreddine_kolf21j)` | 80 | -0.0019 | [-0.0075, +0.0037] | -0.0056 | 0.0000 | [0.0000, 0.0029] | -0.0056 | -0.0246 | -0.0246 | 0.0792 | -0.0020 | -0.3104 | 0.8839 | 0.4776 | 0.0057 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 80 | +0.0006 | [-0.0024, +0.0046] | +0.0277 | 0.0008 | [0.0000, 0.0087] | +0.0276 | +0.0102 | +0.0102 | 0.0607 | +0.0006 | +0.1676 | 0.8813 | 0.4816 | 0.0056 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 80 | -0.0001 | [-0.0041, +0.0045] | +0.0015 | 0.0000 | [0.0000, 0.0057] | +0.0015 | -0.0024 | -0.0024 | 0.0608 | -0.0001 | -0.0389 | 0.8878 | 0.4837 | 0.0056 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_ShuffledGenes` | 80 | -0.0023 | [-0.0083, +0.0044] | +0.0279 | 0.0008 | [0.0000, 0.0056] | +0.0278 | -0.0191 | -0.0190 | 0.1191 | -0.0023 | -0.1597 | 0.8824 | 0.4692 | 0.0064 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_DiagCov` | 80 | -0.0021 | [-0.0079, +0.0044] | +0.0390 | 0.0015 | [0.0000, 0.0106] | +0.0389 | -0.0260 | -0.0260 | 0.0800 | -0.0021 | -0.3248 | 0.8842 | 0.4841 | 0.0056 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 80 | +0.0038 | [-0.0004, +0.0077] | +0.0019 | 0.0000 | [0.0000, 0.0041] | +0.0019 | +0.0631 | +0.0630 | 0.0609 | +0.0039 | +1.0349 | 0.8820 | 0.4827 | 0.0055 |
| Confirmatory | vcc2025_h1 | `Transfer_Mean_VarScaled` | 80 | -0.1655 | [-0.2647, -0.0565] | -0.0782 | 0.0061 | [0.0000, 0.0384] | -0.0780 | -0.4454 | -0.4444 | 0.3723 | -0.1855 | -1.1937 | 0.8896 | 0.4808 | 0.0042 |
| Confirmatory | vcc2025_h1 | `Transfer_Mean_GlobalScale` | 80 | +0.0684 | [-0.0229, +0.1506] | +nan | nan | [nan, nan] | +nan | +0.2506 | +0.2501 | 0.2734 | +0.0735 | +0.9147 | 0.8932 | 0.4900 | 0.0039 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_MismatchedVar(nadig_hepg2)` | 80 | +0.0354 | [+0.0006, +0.0661] | +0.0628 | 0.0039 | [0.0000, 0.0169] | +0.0627 | +0.3159 | +0.3152 | 0.1124 | +0.0358 | +2.8039 | 0.8905 | 0.4877 | 0.0050 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_MismatchedVar(nadig_jurkat)` | 80 | -0.0435 | [-0.0783, -0.0051] | -0.0714 | 0.0051 | [0.0007, 0.0150] | -0.0712 | -0.3597 | -0.3589 | 0.1212 | -0.0441 | -2.9614 | 0.8908 | 0.4830 | 0.0058 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_MismatchedVar(nourreddine_kolf21j)` | 80 | -0.0464 | [-0.0907, +0.0040] | -0.0066 | 0.0000 | [0.0000, 0.0040] | -0.0066 | -0.2875 | -0.2869 | 0.1618 | -0.0476 | -1.7733 | 0.8929 | 0.4865 | 0.0060 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 80 | +0.0381 | [+0.0019, +0.0699] | +0.0632 | 0.0040 | [0.0003, 0.0116] | +0.0631 | +0.3404 | +0.3396 | 0.1121 | +0.0385 | +3.0295 | 0.8938 | 0.4896 | 0.0049 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 80 | +0.0382 | [-0.0024, +0.0753] | -0.0217 | 0.0005 | [0.0000, 0.0079] | -0.0217 | +0.3015 | +0.3009 | 0.1271 | +0.0388 | +2.3673 | 0.8922 | 0.4802 | 0.0054 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_ShuffledVar` | 80 | +0.0604 | [-0.0301, +0.1417] | +0.0340 | 0.0012 | [0.0000, 0.0073] | +0.0339 | +0.1894 | +0.1890 | 0.3195 | +0.0669 | +0.5914 | 0.8929 | 0.4875 | 0.0038 |
| Confirmatory | vcc2025_h1 | `Hybrid_on_TMVS` | 80 | -0.1730 | [-0.2719, -0.0621] | -0.0812 | 0.0066 | [0.0000, 0.0345] | -0.0809 | -0.4155 | -0.4146 | 0.4172 | -0.2009 | -0.9937 | 0.8271 | 0.3979 | 0.0051 |
| Confirmatory | vcc2025_h1 | `FCT` | 80 | +0.2503 | [+0.1339, +0.3542] | +0.3505 | 0.1228 | [0.0775, 0.1802] | +0.3495 | +0.5489 | +0.5478 | 0.4569 | +0.2898 | +1.1988 | 0.9134 | 0.5552 | 0.0069 |

## R18: TMVS Shuffled-Variance 100-Permutation Null vs H-B3 (97.5% CI Threshold)

| Target Line | TMVS r_orth_cal | TMVS CI lo | Perm100 p2.5 | Perm100 median | Perm100 p97.5 | Perm100 max | Max Mismatch CI hi | H-B3 Pass |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| nadig_hepg2 | +0.0421 | +0.0249 | +0.0123 | +0.0541 | +0.0874 | +0.1123 | +0.1154 | False |
| nadig_jurkat | +0.0967 | +0.0822 | +0.0575 | +0.1109 | +0.1529 | +0.1705 | +0.1730 | False |
| nourreddine_kolf21j | +0.1386 | -0.0513 | +0.0233 | +0.0814 | +0.1462 | +0.1759 | +0.1767 | False |
| replogle_k562_essential | +0.1833 | +0.1617 | +0.0202 | +0.0777 | +0.1252 | +0.1352 | +0.1783 | False |
| replogle_rpe1 | +0.0915 | +0.0721 | -0.0187 | +0.0177 | +0.0479 | +0.0560 | +0.0424 | True |
| vcc2025_h1 | -0.0782 | -0.1960 | -0.0220 | +0.0210 | +0.0534 | +0.0827 | +0.1299 | False |

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
| nadig_hepg2 | 799 | 0.8745 | 0.9659 | 0.8784 | 0.8784 | 0.4747 | 0.5524 | +0.1823 | +0.0478 |
| nadig_jurkat | 934 | 0.8465 | 0.9575 | 0.8562 | 0.8562 | 0.3765 | 0.5096 | +0.1229 | +0.0452 |
| nourreddine_kolf21j | 18 | 0.8195 | 0.9491 | 0.8436 | 0.8438 | 0.1627 | 0.4900 | +0.0981 | -0.0139 |
| replogle_k562_essential | 859 | 0.8861 | 0.9693 | 0.9008 | 0.9008 | 0.4547 | 0.6830 | +0.1390 | +0.0221 |
| replogle_rpe1 | 1070 | 0.9134 | 0.9771 | 0.9178 | 0.9178 | 0.6139 | 0.6386 | +0.1764 | +0.1734 |
| vcc2025_h1 | 80 | 0.9885 | 0.9971 | 0.9915 | 0.9915 | 0.8886 | 0.9329 | +0.0319 | +0.0057 |

## S5: Expressed Genes Sensitivity (Drop Bottom 25% Target Control Variance Genes)

| Target Line | Genes Kept | rho_AB_orth (S5) | Method | Slope OLS (S5) | r_orth_cal (S5) | r_orth_cal^2 (S5) | Joint 97.5% CI | r_cal raw (S5) |
|---|---:|---:|---|---:|---:|---:|---|---:|
| nadig_hepg2 | 4611/6148 | 0.8744 | `Hybrid_ProjShared` | +0.0094 | +0.0477 | 0.0023 | [+0.0224, +0.0728] | +0.0706 |
| nadig_hepg2 | 4611/6148 | 0.8744 | `Transfer_Mean_VarScaled` | +0.0359 | +0.0426 | 0.0018 | [+0.0252, +0.0612] | +0.1229 |
| nadig_hepg2 | 4611/6148 | 0.8744 | `Transfer_Mean_GlobalScale` | -0.0084 | +nan | nan | [nan, nan] | -0.1295 |
| nadig_hepg2 | 4611/6148 | 0.8744 | `Hybrid_on_TMVS` | +0.0453 | +0.0548 | 0.0030 | [+0.0350, +0.0760] | +0.1288 |
| nadig_hepg2 | 4611/6148 | 0.8744 | `FCT` | +0.0591 | +0.2624 | 0.0688 | [+0.2432, +0.2820] | +0.2637 |
| nadig_hepg2 | 4611/6148 | 0.8744 | `NegCtrl_HybridProj_SrcCorrTgtVar` | +0.0035 | +0.0420 | 0.0018 | [+0.0290, +0.0554] | +0.0531 |
| nadig_hepg2 | 4611/6148 | 0.8744 | `NegCtrl_TMVS_ShuffledVar` | -0.0193 | +0.0592 | 0.0035 | [+0.0481, +0.0701] | -0.0612 |
| nadig_jurkat | 4611/6148 | 0.8465 | `Hybrid_ProjShared` | +0.0060 | +0.0662 | 0.0044 | [+0.0420, +0.0907] | +0.0584 |
| nadig_jurkat | 4611/6148 | 0.8465 | `Transfer_Mean_VarScaled` | +0.0097 | +0.0987 | 0.0097 | [+0.0847, +0.1131] | +0.0446 |
| nadig_jurkat | 4611/6148 | 0.8465 | `Transfer_Mean_GlobalScale` | -0.0013 | +nan | nan | [nan, nan] | -0.0094 |
| nadig_jurkat | 4611/6148 | 0.8465 | `Hybrid_on_TMVS` | +0.0157 | +0.1148 | 0.0132 | [+0.0964, +0.1331] | +0.0644 |
| nadig_jurkat | 4611/6148 | 0.8465 | `FCT` | +0.0198 | +0.1043 | 0.0109 | [+0.0845, +0.1239] | +0.0914 |
| nadig_jurkat | 4611/6148 | 0.8465 | `NegCtrl_HybridProj_SrcCorrTgtVar` | +0.0015 | +0.0268 | 0.0007 | [+0.0108, +0.0434] | +0.0231 |
| nadig_jurkat | 4611/6148 | 0.8465 | `NegCtrl_TMVS_ShuffledVar` | +0.0437 | +0.1282 | 0.0164 | [+0.1158, +0.1407] | +0.1489 |
| nourreddine_kolf21j | 4611/6148 | 0.8198 | `Hybrid_ProjShared` | -0.0147 | -0.0017 | 0.0000 | [-0.1466, +0.1282] | -0.0635 |
| nourreddine_kolf21j | 4611/6148 | 0.8198 | `Transfer_Mean_VarScaled` | -0.0134 | +0.1379 | 0.0190 | [-0.0547, +0.2864] | -0.0396 |
| nourreddine_kolf21j | 4611/6148 | 0.8198 | `Transfer_Mean_GlobalScale` | -0.0681 | +nan | nan | [nan, nan] | -0.3297 |
| nourreddine_kolf21j | 4611/6148 | 0.8198 | `Hybrid_on_TMVS` | -0.0281 | +0.1072 | 0.0115 | [-0.0654, +0.2673] | -0.0656 |
| nourreddine_kolf21j | 4611/6148 | 0.8198 | `FCT` | +0.0744 | +0.1835 | 0.0337 | [+0.0212, +0.3537] | +0.2538 |
| nourreddine_kolf21j | 4611/6148 | 0.8198 | `NegCtrl_HybridProj_SrcCorrTgtVar` | -0.0026 | +0.0903 | 0.0082 | [-0.0922, +0.2489] | -0.0142 |
| nourreddine_kolf21j | 4611/6148 | 0.8198 | `NegCtrl_TMVS_ShuffledVar` | +0.0333 | +0.0388 | 0.0015 | [-0.0073, +0.1065] | +0.1172 |
| replogle_k562_essential | 4611/6148 | 0.8847 | `Hybrid_ProjShared` | +0.0021 | +0.0357 | 0.0013 | [+0.0104, +0.0596] | +0.0166 |
| replogle_k562_essential | 4611/6148 | 0.8847 | `Transfer_Mean_VarScaled` | +0.0568 | +0.1881 | 0.0354 | [+0.1666, +0.2102] | +0.2396 |
| replogle_k562_essential | 4611/6148 | 0.8847 | `Transfer_Mean_GlobalScale` | +0.0027 | +nan | nan | [nan, nan] | +0.2804 |
| replogle_k562_essential | 4611/6148 | 0.8847 | `Hybrid_on_TMVS` | +0.0589 | +0.1786 | 0.0319 | [+0.1559, +0.2015] | +0.2183 |
| replogle_k562_essential | 4611/6148 | 0.8847 | `FCT` | +0.0424 | +0.1966 | 0.0386 | [+0.1778, +0.2170] | +0.1482 |
| replogle_k562_essential | 4611/6148 | 0.8847 | `NegCtrl_HybridProj_SrcCorrTgtVar` | +0.0019 | +0.0247 | 0.0006 | [+0.0050, +0.0430] | +0.0231 |
| replogle_k562_essential | 4611/6148 | 0.8847 | `NegCtrl_TMVS_ShuffledVar` | +0.1038 | +0.0661 | 0.0044 | [+0.0551, +0.0771] | +0.2489 |
| replogle_rpe1 | 4611/6148 | 0.9128 | `Hybrid_ProjShared` | +0.0133 | +0.1963 | 0.0385 | [+0.1707, +0.2213] | +0.1728 |
| replogle_rpe1 | 4611/6148 | 0.9128 | `Transfer_Mean_VarScaled` | +0.0351 | +0.0898 | 0.0081 | [+0.0719, +0.1072] | +0.2711 |
| replogle_rpe1 | 4611/6148 | 0.9128 | `Transfer_Mean_GlobalScale` | +0.0083 | +nan | nan | [nan, nan] | +0.3631 |
| replogle_rpe1 | 4611/6148 | 0.9128 | `Hybrid_on_TMVS` | +0.0484 | +0.1775 | 0.0315 | [+0.1567, +0.1978] | +0.3032 |
| replogle_rpe1 | 4611/6148 | 0.9128 | `FCT` | +0.0253 | +0.1318 | 0.0174 | [+0.1144, +0.1498] | +0.1815 |
| replogle_rpe1 | 4611/6148 | 0.9128 | `NegCtrl_HybridProj_SrcCorrTgtVar` | -0.0014 | -0.0001 | 0.0000 | [-0.0213, +0.0223] | -0.0346 |
| replogle_rpe1 | 4611/6148 | 0.9128 | `NegCtrl_TMVS_ShuffledVar` | -0.0356 | +0.0276 | 0.0008 | [+0.0180, +0.0376] | -0.1942 |
| vcc2025_h1 | 4611/6148 | 0.9885 | `Hybrid_ProjShared` | -0.0081 | -0.0369 | 0.0014 | [-0.1186, +0.0376] | -0.0459 |
| vcc2025_h1 | 4611/6148 | 0.9885 | `Transfer_Mean_VarScaled` | -0.1649 | -0.0801 | 0.0064 | [-0.1960, +0.0274] | -0.4494 |
| vcc2025_h1 | 4611/6148 | 0.9885 | `Transfer_Mean_GlobalScale` | +0.0689 | +nan | nan | [nan, nan] | +0.2556 |
| vcc2025_h1 | 4611/6148 | 0.9885 | `Hybrid_on_TMVS` | -0.1731 | -0.0819 | 0.0067 | [-0.1857, +0.0126] | -0.4205 |
| vcc2025_h1 | 4611/6148 | 0.9885 | `FCT` | +0.2509 | +0.3495 | 0.1222 | [+0.2753, +0.4252] | +0.5498 |
| vcc2025_h1 | 4611/6148 | 0.9885 | `NegCtrl_HybridProj_SrcCorrTgtVar` | +0.0037 | -0.0010 | 0.0000 | [-0.0538, +0.0620] | +0.0617 |
| vcc2025_h1 | 4611/6148 | 0.9885 | `NegCtrl_TMVS_ShuffledVar` | +0.0613 | +0.0322 | 0.0010 | [-0.0082, +0.0847] | +0.1937 |

## Appendix: Alternative H-A4/H-B4 Discrimination & Centroid Conventions (R5)

| Target Line | Method | Ctrl disc_raw (A-T2) | Ctrl disc_raw_panel (A-T3) | Cent disc_raw | Cent disc_raw_panel | Cent de50_r | pds_cosine_panel | Slope OLS (All) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| nadig_hepg2 | `Source_Grand_Mean` | 0.5003 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.3059 |
| nadig_hepg2 | `Transfer_Mean` | 0.8811 | 0.8811 | 0.9026 | 0.9028 | 0.5411 | 0.8790 | +0.0000 |
| nadig_hepg2 | `CIPHER_Linear` | 0.4832 | 0.4858 | 0.4889 | 0.4927 | -0.0772 | 0.4848 | +0.3257 |
| nadig_hepg2 | `Hybrid_Raw` | 0.8603 | 0.8612 | 0.8931 | 0.8936 | 0.5530 | 0.8581 | +0.0121 |
| nadig_hepg2 | `Hybrid_ProjShared` | 0.8693 | 0.8699 | 0.8976 | 0.8979 | 0.5452 | 0.8679 | +0.0079 |
| nadig_hepg2 | `NegCtrl_HybridProj_ShuffledGenes` | 0.8774 | 0.8763 | 0.8991 | 0.8995 | 0.5295 | 0.8749 | +0.0023 |
| nadig_hepg2 | `NegCtrl_HybridProj_DiagCov` | 0.8781 | 0.8772 | 0.8994 | 0.9001 | 0.5374 | 0.8757 | +0.0021 |
| nadig_hepg2 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.8808 | 0.8805 | 0.9023 | 0.9024 | 0.5408 | 0.8782 | +0.0020 |
| nadig_hepg2 | `Transfer_Mean_VarScaled` | 0.8745 | 0.8754 | 0.8998 | 0.9007 | 0.5398 | 0.8736 | +0.0329 |
| nadig_hepg2 | `Transfer_Mean_GlobalScale` | 0.8811 | 0.8811 | 0.9039 | 0.9039 | 0.5361 | 0.8790 | +0.0247 |
| nadig_hepg2 | `NegCtrl_TMVS_ShuffledVar` | 0.8902 | 0.8900 | 0.9011 | 0.9029 | 0.4722 | 0.8881 | +0.0625 |
| nadig_hepg2 | `Hybrid_on_TMVS` | 0.8615 | 0.8631 | 0.8940 | 0.8946 | 0.5368 | 0.8617 | +0.0408 |
| nadig_hepg2 | `FCT` | 0.8871 | 0.8883 | 0.9080 | 0.9094 | 0.5682 | 0.8852 | +0.0576 |
| nadig_jurkat | `Source_Grand_Mean` | 0.5002 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.3056 |
| nadig_jurkat | `Transfer_Mean` | 0.8717 | 0.8637 | 0.9218 | 0.9182 | 0.5323 | 0.8530 | +0.0000 |
| nadig_jurkat | `CIPHER_Linear` | 0.4831 | 0.4873 | 0.4936 | 0.4950 | 0.4030 | 0.4953 | +0.4096 |
| nadig_jurkat | `Hybrid_Raw` | 0.8744 | 0.8675 | 0.9232 | 0.9197 | 0.5649 | 0.8558 | +0.0207 |
| nadig_jurkat | `Hybrid_ProjShared` | 0.8685 | 0.8610 | 0.9193 | 0.9159 | 0.5353 | 0.8506 | +0.0051 |
| nadig_jurkat | `NegCtrl_HybridProj_ShuffledGenes` | 0.8707 | 0.8620 | 0.9199 | 0.9163 | 0.5214 | 0.8517 | -0.0009 |
| nadig_jurkat | `NegCtrl_HybridProj_DiagCov` | 0.8713 | 0.8625 | 0.9203 | 0.9169 | 0.5237 | 0.8521 | -0.0011 |
| nadig_jurkat | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.8709 | 0.8629 | 0.9218 | 0.9182 | 0.5331 | 0.8524 | +0.0014 |
| nadig_jurkat | `Transfer_Mean_VarScaled` | 0.8710 | 0.8629 | 0.9134 | 0.9100 | 0.5051 | 0.8514 | -0.0496 |
| nadig_jurkat | `Transfer_Mean_GlobalScale` | 0.8717 | 0.8637 | 0.9141 | 0.9104 | 0.4918 | 0.8530 | -0.0648 |
| nadig_jurkat | `NegCtrl_TMVS_ShuffledVar` | 0.8814 | 0.8756 | 0.9308 | 0.9282 | 0.5928 | 0.8661 | +0.0072 |
| nadig_jurkat | `Hybrid_on_TMVS` | 0.8683 | 0.8607 | 0.9115 | 0.9077 | 0.5087 | 0.8496 | -0.0445 |
| nadig_jurkat | `FCT` | 0.8603 | 0.8519 | 0.9042 | 0.9006 | 0.5139 | 0.8359 | -0.0071 |
| nourreddine_kolf21j | `Source_Grand_Mean` | 0.5000 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.4560 |
| nourreddine_kolf21j | `Transfer_Mean` | 0.8137 | 0.8105 | 0.7778 | 0.7778 | 0.4715 | 0.8072 | +0.0000 |
| nourreddine_kolf21j | `CIPHER_Linear` | 0.5915 | 0.5850 | 0.5294 | 0.5261 | 0.6786 | 0.5065 | +0.6232 |
| nourreddine_kolf21j | `Hybrid_Raw` | 0.7843 | 0.7843 | 0.7582 | 0.7582 | 0.5071 | 0.7484 | +0.0206 |
| nourreddine_kolf21j | `Hybrid_ProjShared` | 0.8072 | 0.8072 | 0.7484 | 0.7516 | 0.4529 | 0.8007 | -0.0020 |
| nourreddine_kolf21j | `NegCtrl_HybridProj_ShuffledGenes` | 0.7712 | 0.7712 | 0.7190 | 0.7190 | 0.4344 | 0.7712 | +0.0020 |
| nourreddine_kolf21j | `NegCtrl_HybridProj_DiagCov` | 0.7712 | 0.7712 | 0.7288 | 0.7288 | 0.4472 | 0.7778 | +0.0027 |
| nourreddine_kolf21j | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.8072 | 0.8072 | 0.7418 | 0.7451 | 0.4720 | 0.8007 | +0.0026 |
| nourreddine_kolf21j | `Transfer_Mean_VarScaled` | 0.8105 | 0.8105 | 0.7647 | 0.7647 | 0.5043 | 0.7941 | -0.1104 |
| nourreddine_kolf21j | `Transfer_Mean_GlobalScale` | 0.8137 | 0.8105 | 0.7745 | 0.7745 | 0.3788 | 0.8072 | -0.1532 |
| nourreddine_kolf21j | `NegCtrl_TMVS_ShuffledVar` | 0.8105 | 0.8105 | 0.7712 | 0.7712 | 0.5581 | 0.8007 | -0.0321 |
| nourreddine_kolf21j | `Hybrid_on_TMVS` | 0.8039 | 0.8007 | 0.7647 | 0.7647 | 0.4917 | 0.8007 | -0.1123 |
| nourreddine_kolf21j | `FCT` | 0.8170 | 0.8170 | 0.7712 | 0.7778 | 0.5354 | 0.7810 | +0.0620 |
| replogle_k562_essential | `Source_Grand_Mean` | 0.5001 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.4825 |
| replogle_k562_essential | `Transfer_Mean` | 0.8868 | 0.8927 | 0.9374 | 0.9390 | 0.4539 | 0.8788 | +0.0000 |
| replogle_k562_essential | `CIPHER_Linear` | 0.4395 | 0.4437 | 0.4881 | 0.4897 | 0.4609 | 0.4619 | +0.6934 |
| replogle_k562_essential | `Hybrid_Raw` | 0.8930 | 0.8975 | 0.9478 | 0.9486 | 0.5109 | 0.8891 | +0.0458 |
| replogle_k562_essential | `Hybrid_ProjShared` | 0.8841 | 0.8889 | 0.9365 | 0.9375 | 0.4490 | 0.8752 | +0.0019 |
| replogle_k562_essential | `NegCtrl_HybridProj_ShuffledGenes` | 0.8930 | 0.8958 | 0.9375 | 0.9380 | 0.4328 | 0.8829 | -0.0032 |
| replogle_k562_essential | `NegCtrl_HybridProj_DiagCov` | 0.8930 | 0.8960 | 0.9377 | 0.9383 | 0.4342 | 0.8831 | -0.0034 |
| replogle_k562_essential | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.8870 | 0.8924 | 0.9373 | 0.9390 | 0.4571 | 0.8787 | +0.0006 |
| replogle_k562_essential | `Transfer_Mean_VarScaled` | 0.8885 | 0.8949 | 0.9425 | 0.9449 | 0.5272 | 0.8801 | +0.0409 |
| replogle_k562_essential | `Transfer_Mean_GlobalScale` | 0.8868 | 0.8927 | 0.9380 | 0.9396 | 0.4568 | 0.8788 | +0.0057 |
| replogle_k562_essential | `NegCtrl_TMVS_ShuffledVar` | 0.8921 | 0.8995 | 0.9541 | 0.9552 | 0.5779 | 0.8880 | +0.1133 |
| replogle_k562_essential | `Hybrid_on_TMVS` | 0.8856 | 0.8910 | 0.9414 | 0.9433 | 0.5191 | 0.8764 | +0.0428 |
| replogle_k562_essential | `FCT` | 0.8763 | 0.8846 | 0.9309 | 0.9368 | 0.4754 | 0.8610 | +0.0440 |
| replogle_rpe1 | `Source_Grand_Mean` | 0.5005 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.1628 |
| replogle_rpe1 | `Transfer_Mean` | 0.9301 | 0.9346 | 0.9173 | 0.9202 | 0.4525 | 0.9287 | +0.0000 |
| replogle_rpe1 | `CIPHER_Linear` | 0.5123 | 0.5157 | 0.5119 | 0.5139 | -0.2174 | 0.5157 | +0.1235 |
| replogle_rpe1 | `Hybrid_Raw` | 0.9200 | 0.9235 | 0.9315 | 0.9340 | 0.5093 | 0.9181 | +0.0222 |
| replogle_rpe1 | `Hybrid_ProjShared` | 0.9249 | 0.9292 | 0.9275 | 0.9305 | 0.4916 | 0.9237 | +0.0159 |
| replogle_rpe1 | `NegCtrl_HybridProj_ShuffledGenes` | 0.9284 | 0.9322 | 0.9089 | 0.9116 | 0.4221 | 0.9266 | -0.0081 |
| replogle_rpe1 | `NegCtrl_HybridProj_DiagCov` | 0.9277 | 0.9316 | 0.9086 | 0.9112 | 0.4265 | 0.9261 | -0.0079 |
| replogle_rpe1 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.9293 | 0.9333 | 0.9156 | 0.9186 | 0.4494 | 0.9275 | -0.0009 |
| replogle_rpe1 | `Transfer_Mean_VarScaled` | 0.9254 | 0.9301 | 0.9238 | 0.9256 | 0.4894 | 0.9251 | +0.0132 |
| replogle_rpe1 | `Transfer_Mean_GlobalScale` | 0.9301 | 0.9346 | 0.9208 | 0.9234 | 0.4611 | 0.9287 | -0.0036 |
| replogle_rpe1 | `NegCtrl_TMVS_ShuffledVar` | 0.9344 | 0.9387 | 0.9018 | 0.9073 | 0.3558 | 0.9327 | -0.0071 |
| replogle_rpe1 | `Hybrid_on_TMVS` | 0.9208 | 0.9257 | 0.9312 | 0.9338 | 0.5186 | 0.9207 | +0.0291 |
| replogle_rpe1 | `FCT` | 0.9308 | 0.9347 | 0.9181 | 0.9220 | 0.4648 | 0.9290 | +0.0358 |
| vcc2025_h1 | `Source_Grand_Mean` | 0.5003 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.3627 |
| vcc2025_h1 | `Transfer_Mean` | 0.8932 | 0.8921 | 0.8354 | 0.8362 | 0.5056 | 0.8758 | +0.0000 |
| vcc2025_h1 | `CIPHER_Linear` | 0.4672 | 0.4701 | 0.4813 | 0.4828 | 0.3693 | 0.4639 | +0.4680 |
| vcc2025_h1 | `Hybrid_Raw` | 0.8595 | 0.8604 | 0.8266 | 0.8269 | 0.4691 | 0.8565 | +0.0006 |
| vcc2025_h1 | `Hybrid_ProjShared` | 0.8535 | 0.8527 | 0.8297 | 0.8293 | 0.4576 | 0.8426 | -0.0002 |
| vcc2025_h1 | `NegCtrl_HybridProj_ShuffledGenes` | 0.8824 | 0.8820 | 0.8340 | 0.8339 | 0.4811 | 0.8699 | -0.0006 |
| vcc2025_h1 | `NegCtrl_HybridProj_DiagCov` | 0.8842 | 0.8832 | 0.8329 | 0.8343 | 0.4960 | 0.8712 | -0.0003 |
| vcc2025_h1 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.8820 | 0.8788 | 0.8348 | 0.8358 | 0.4955 | 0.8650 | +0.0006 |
| vcc2025_h1 | `Transfer_Mean_VarScaled` | 0.8896 | 0.8888 | 0.7248 | 0.7220 | 0.5596 | 0.8750 | +0.1925 |
| vcc2025_h1 | `Transfer_Mean_GlobalScale` | 0.8932 | 0.8921 | 0.7321 | 0.7331 | 0.5563 | 0.8758 | +0.1982 |
| vcc2025_h1 | `NegCtrl_TMVS_ShuffledVar` | 0.8929 | 0.8908 | 0.6888 | 0.6886 | 0.5679 | 0.8771 | +0.2347 |
| vcc2025_h1 | `Hybrid_on_TMVS` | 0.8271 | 0.8266 | 0.6954 | 0.6948 | 0.4874 | 0.8180 | +0.1924 |
| vcc2025_h1 | `FCT` | 0.9134 | 0.9123 | 0.8905 | 0.8905 | 0.6073 | 0.8998 | +0.0427 |
