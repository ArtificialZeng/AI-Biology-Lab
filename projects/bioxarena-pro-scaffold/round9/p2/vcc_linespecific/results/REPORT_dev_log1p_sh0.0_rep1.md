# Line-Specific Zero-Shot Perturbation Evaluation Report (`dev_log1p_sh0.0_rep1`)

- **Confirmatory Lines (2)**: replogle_k562_essential, replogle_rpe1
- **Secondary Lines (0)**: None
- **Dual Confirmatory Summary (R13, 97.5% Bonferroni CIs)**: **Group A (Hybrid_ProjShared): 未确认 | Group B (Transfer_Mean_VarScaled): 未确认 (Bonferroni 97.5% CIs)**
- **Group A Primary (`Hybrid_ProjShared`)**: **未确认**
- **Group B Primary (`Transfer_Mean_VarScaled`)**: **未确认**
- **R4 Co-fluctuation Structure Claim Allowed**: `目标线对照统计量` (`all beat SrcCorrTgtVar disjoint=False`)
- **R4 Marginal Variance Note**: `增益可由边际方差解释`
- **R14 Mechanism Statement (`Hybrid_on_TMVS` vs `TMVS`)**: `在本数据上，共波动结构没有超出边际方差的可检测增量` (`0/2`)

## Hypothesis Adjudication Summary (Confirmatory Cohort, 97.5% CIs)

| Hypothesis | Group | Pass / Confirmatory | Threshold | Verdict |
|---|---|---:|---|---|
| **H-A1** (Slope 97.5% CI lo > 0.13) | A (`Hybrid_ProjShared`) | 0/2 | >= 2 | **未确认** |
| **H-A2** (Joint-boot `r_cal` 97.5% CI lo > 0, valid>=95%) | A (`Hybrid_ProjShared`) | 2/2 (fixed95: 2/2) | >= 2 | **CONFIRMED** |
| **H-A3** (Matched > ALL source mismatched cov) | A (`Hybrid_ProjShared`) | 2/2 | >= 2 | **CONFIRMED** |
| **H-A4** (Strict AND non-inferiority across 6 conventions) | A (`Hybrid_ProjShared`) | 2/2 inf (ctrl_disc_raw=0/2, ctrl_de50_r=0/2, ctrl_disc_raw_panel=2/2, cent_disc_raw=1/2, cent_de50_r=0/2, cent_disc_raw_panel=1/2; union=2) | <= 0 inf | **未确认** |
| **H-B1** (Slope 97.5% CI lo > 0.13) | B (`Transfer_Mean_VarScaled`) | 0/2 | >= 2 | **未确认** |
| **H-B2** (Joint-boot `r_cal` 97.5% CI lo > 0, valid>=95%) | B (`Transfer_Mean_VarScaled`) | 2/2 | >= 2 | **CONFIRMED** |
| **H-B3** (Matched var > ALL mismatched + shuffled var on slope & `r_cal`) | B (`Transfer_Mean_VarScaled`) | 1/2 | >= 2 | **未确认** |
| **H-B4** (Strict AND non-inferiority across 6 conventions) | B (`Transfer_Mean_VarScaled`) | 0/2 inf (ctrl_disc_raw=0/2, ctrl_de50_r=0/2, ctrl_disc_raw_panel=0/2, cent_disc_raw=0/2, cent_de50_r=0/2, cent_disc_raw_panel=0/2; union=0) | <= 0 inf | **CONFIRMED** |

## Main Table: Primary View (A-T2 Own-Target-Zeroed Control Reference, 97.5% CIs, R5/R6/R11/R13)

| Cohort | Target Line | Method | N | Slope OLS (Top50) | 97.5% CI | Calib Resid r | r_cal^2 | r_cal^2 Joint 97.5% CI | Resid r | sd_pred/sd_obs | Deming (delta=1) | Obs~Pred | Ctrl disc_raw | Ctrl de50_r | MSE |
|---|---|---|---:|---:|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|
| Confirmatory | replogle_k562_essential | `Source_Grand_Mean` | 832 | +0.4606 | [+0.4124, +0.5100] | +0.4695 | 0.2204 | [0.1865, 0.2574] | +0.4634 | 0.9938 | +0.9866 | +0.4663 | 0.5002 | 0.3439 | 0.0079 |
| Confirmatory | replogle_k562_essential | `Transfer_Mean` | 832 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | 0.0000 | +0.0000 | +nan | 0.8201 | 0.5996 | 0.0150 |
| Confirmatory | replogle_k562_essential | `CIPHER_Linear` | 832 | +0.6998 | [+0.6517, +0.7487] | +0.5631 | 0.3171 | [0.2828, 0.3538] | +0.5559 | 1.2590 | +1.5018 | +0.4415 | 0.4811 | 0.0683 | 0.0051 |
| Confirmatory | replogle_k562_essential | `Hybrid_Raw` | 832 | +0.0200 | [+0.0150, +0.0248] | +0.1379 | 0.0190 | [0.0101, 0.0307] | +0.1361 | 0.1467 | +0.0204 | +0.9281 | 0.8221 | 0.6022 | 0.0149 |
| Confirmatory | replogle_k562_essential | `Hybrid_ProjShared` | 832 | +0.0050 | [+0.0029, +0.0071] | +0.0523 | 0.0027 | [0.0009, 0.0054] | +0.0516 | 0.0964 | +0.0050 | +0.5356 | 0.8190 | 0.6021 | 0.0152 |
| Confirmatory | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 832 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | 0.0000 | +0.0000 | +nan | 0.8201 | 0.5996 | 0.0150 |
| Confirmatory | replogle_k562_essential | `NegCtrl_HybridProj_ShuffledGenes` | 832 | -0.0020 | [-0.0046, +0.0006] | -0.0160 | 0.0003 | [0.0000, 0.0013] | -0.0158 | 0.1274 | -0.0020 | -0.1243 | 0.8220 | 0.5951 | 0.0154 |
| Confirmatory | replogle_k562_essential | `NegCtrl_HybridProj_DiagCov` | 832 | -0.0023 | [-0.0049, +0.0003] | -0.0191 | 0.0004 | [0.0000, 0.0017] | -0.0188 | 0.1218 | -0.0023 | -0.1545 | 0.8221 | 0.5954 | 0.0153 |
| Confirmatory | replogle_k562_essential | `NegCtrl_HybridProj_SrcCorrTgtVar` | 832 | +0.0011 | [+0.0004, +0.0018] | +0.0292 | 0.0009 | [0.0001, 0.0023] | +0.0288 | 0.0380 | +0.0011 | +0.7581 | 0.8206 | 0.5994 | 0.0150 |
| Confirmatory | replogle_k562_essential | `Transfer_Mean_VarScaled` | 832 | +0.0737 | [+0.0629, +0.0844] | +0.2198 | 0.0483 | [0.0344, 0.0646] | +0.2170 | 0.3396 | +0.0827 | +0.6390 | 0.8287 | 0.6116 | 0.0126 |
| Confirmatory | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 832 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | 0.0000 | +0.0000 | +nan | 0.8201 | 0.5996 | 0.0150 |
| Confirmatory | replogle_k562_essential | `NegCtrl_TMVS_ShuffledVar` | 832 | +0.1715 | [+0.1500, +0.1923] | +0.3169 | 0.1004 | [0.0741, 0.1307] | +0.3128 | 0.5483 | +0.2320 | +0.5705 | 0.8343 | 0.5883 | 0.0114 |
| Confirmatory | replogle_k562_essential | `Hybrid_on_TMVS` | 832 | +0.0787 | [+0.0676, +0.0896] | +0.2254 | 0.0508 | [0.0365, 0.0669] | +0.2225 | 0.3536 | +0.0892 | +0.6293 | 0.8272 | 0.6136 | 0.0127 |
| Confirmatory | replogle_k562_essential | `FCT` | 832 | +0.0462 | [+0.0378, +0.0539] | +0.1394 | 0.0194 | [0.0128, 0.0268] | +0.1376 | 0.3357 | +0.0519 | +0.4099 | 0.8117 | 0.6222 | 0.0153 |
| Confirmatory | replogle_rpe1 | `Source_Grand_Mean` | 832 | -0.0500 | [-0.0620, -0.0379] | -0.1384 | 0.0192 | [0.0112, 0.0286] | -0.1364 | 0.3665 | -0.0575 | -0.3722 | 0.5005 | 0.6210 | 0.0143 |
| Confirmatory | replogle_rpe1 | `Transfer_Mean` | 832 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | 0.0000 | +0.0000 | +nan | 0.9274 | 0.5548 | 0.0150 |
| Confirmatory | replogle_rpe1 | `CIPHER_Linear` | 832 | -0.1107 | [-0.1241, -0.0969] | -0.2652 | 0.0703 | [0.0550, 0.0871] | -0.2613 | 0.4234 | -0.1325 | -0.6172 | 0.5113 | 0.3256 | 0.0160 |
| Confirmatory | replogle_rpe1 | `Hybrid_Raw` | 832 | +0.0091 | [+0.0039, +0.0141] | +0.0951 | 0.0091 | [0.0015, 0.0237] | +0.0938 | 0.0974 | +0.0092 | +0.9625 | 0.9219 | 0.5785 | 0.0149 |
| Confirmatory | replogle_rpe1 | `Hybrid_ProjShared` | 832 | +0.0085 | [+0.0044, +0.0125] | +0.1054 | 0.0111 | [0.0026, 0.0258] | +0.1038 | 0.0823 | +0.0086 | +1.2622 | 0.9226 | 0.5788 | 0.0149 |
| Confirmatory | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 832 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | 0.0000 | +0.0000 | +nan | 0.9274 | 0.5548 | 0.0150 |
| Confirmatory | replogle_rpe1 | `NegCtrl_HybridProj_ShuffledGenes` | 832 | -0.0011 | [-0.0018, -0.0004] | -0.0288 | 0.0008 | [0.0001, 0.0022] | -0.0283 | 0.0395 | -0.0011 | -0.7171 | 0.9259 | 0.5435 | 0.0153 |
| Confirmatory | replogle_rpe1 | `NegCtrl_HybridProj_DiagCov` | 832 | -0.0010 | [-0.0016, -0.0003] | -0.0286 | 0.0008 | [0.0001, 0.0023] | -0.0282 | 0.0348 | -0.0010 | -0.8089 | 0.9258 | 0.5455 | 0.0151 |
| Confirmatory | replogle_rpe1 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 832 | +0.0018 | [+0.0014, +0.0021] | +0.1144 | 0.0131 | [0.0081, 0.0193] | +0.1128 | 0.0159 | +0.0018 | +7.1078 | 0.9271 | 0.5597 | 0.0150 |
| Confirmatory | replogle_rpe1 | `Transfer_Mean_VarScaled` | 832 | +0.0267 | [+0.0243, +0.0290] | +0.2145 | 0.0460 | [0.0386, 0.0540] | +0.2114 | 0.1262 | +0.0271 | +1.6753 | 0.9275 | 0.5544 | 0.0149 |
| Confirmatory | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 832 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | 0.0000 | +0.0000 | +nan | 0.9274 | 0.5548 | 0.0150 |
| Confirmatory | replogle_rpe1 | `NegCtrl_TMVS_ShuffledVar` | 832 | -0.0099 | [-0.0123, -0.0073] | -0.0599 | 0.0036 | [0.0019, 0.0058] | -0.0591 | 0.1673 | -0.0102 | -0.3531 | 0.9287 | 0.5355 | 0.0155 |
| Confirmatory | replogle_rpe1 | `Hybrid_on_TMVS` | 832 | +0.0352 | [+0.0310, +0.0393] | +0.2351 | 0.0553 | [0.0428, 0.0684] | +0.2317 | 0.1521 | +0.0360 | +1.5236 | 0.9226 | 0.5753 | 0.0147 |
| Confirmatory | replogle_rpe1 | `FCT` | 832 | +0.0111 | [+0.0089, +0.0135] | +0.0903 | 0.0082 | [0.0053, 0.0118] | +0.0890 | 0.1250 | +0.0113 | +0.7122 | 0.9292 | 0.5462 | 0.0154 |

## S6: Split-Half Noise Ceiling (`resp_top50` vs Unselected `all_nontarget`, Ceiling Selection Bias)

| Target Line | N half-valid | rho_AB (top50) | Ceil r (top50) | Rep Slope (top50) | rho_AB (all_nontarget) | Ceil r (all_nontarget) | TMVS r_cal (all) | Hybrid r_cal (all) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| replogle_k562_essential | 832 | 0.9503 | 0.9872 | 0.9503 | 0.8421 | 0.9562 | +0.4659 | +0.0286 |
| replogle_rpe1 | 832 | 0.9441 | 0.9855 | 0.9441 | 0.7160 | 0.9135 | +0.1214 | +0.1206 |

## S5: Expressed Genes Sensitivity (Drop Bottom 25% Target Control Variance Genes)

| Target Line | Genes Kept | rho_AB (S5) | Method | Slope OLS (S5) | Calib Resid r (S5) | r_cal^2 (S5) | Joint 97.5% CI |
|---|---:|---:|---|---:|---:|---:|---|
| replogle_k562_essential | 5419/7225 | 0.9495 | `Hybrid_ProjShared` | +0.0052 | +0.0549 | 0.0030 | [+0.0330, +0.0765] |
| replogle_k562_essential | 5419/7225 | 0.9495 | `Transfer_Mean_VarScaled` | +0.0777 | +0.2304 | 0.0531 | [+0.1960, +0.2643] |
| replogle_k562_essential | 5419/7225 | 0.9495 | `Hybrid_on_TMVS` | +0.0829 | +0.2360 | 0.0557 | [+0.2021, +0.2695] |
| replogle_k562_essential | 5419/7225 | 0.9495 | `FCT` | +0.0478 | +0.1436 | 0.0206 | [+0.1180, +0.1675] |
| replogle_k562_essential | 5419/7225 | 0.9495 | `NegCtrl_HybridProj_SrcCorrTgtVar` | +0.0010 | +0.0256 | 0.0007 | [+0.0069, +0.0450] |
| replogle_k562_essential | 5419/7225 | 0.9495 | `NegCtrl_TMVS_ShuffledVar` | +0.1982 | +0.4005 | 0.1604 | [+0.3769, +0.4241] |
| replogle_rpe1 | 5419/7225 | 0.9434 | `Hybrid_ProjShared` | +0.0087 | +0.1063 | 0.0113 | [+0.0519, +0.1624] |
| replogle_rpe1 | 5419/7225 | 0.9434 | `Transfer_Mean_VarScaled` | +0.0285 | +0.2327 | 0.0542 | [+0.2140, +0.2508] |
| replogle_rpe1 | 5419/7225 | 0.9434 | `Hybrid_on_TMVS` | +0.0372 | +0.2501 | 0.0626 | [+0.2205, +0.2779] |
| replogle_rpe1 | 5419/7225 | 0.9434 | `FCT` | +0.0133 | +0.1079 | 0.0116 | [+0.0913, +0.1254] |
| replogle_rpe1 | 5419/7225 | 0.9434 | `NegCtrl_HybridProj_SrcCorrTgtVar` | +0.0019 | +0.1198 | 0.0144 | [+0.0958, +0.1440] |
| replogle_rpe1 | 5419/7225 | 0.9434 | `NegCtrl_TMVS_ShuffledVar` | -0.0123 | -0.0853 | 0.0073 | [-0.0996, -0.0708] |

## Appendix: Alternative H-A4/H-B4 Discrimination & Centroid Conventions (R5)

| Target Line | Method | Ctrl disc_raw (A-T2) | Ctrl disc_raw_panel (A-T3) | Cent disc_raw | Cent disc_raw_panel | Cent de50_r | pds_cosine_panel | Slope OLS (All) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| replogle_k562_essential | `Source_Grand_Mean` | 0.5002 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.5329 |
| replogle_k562_essential | `Transfer_Mean` | 0.8201 | 0.8237 | 0.8784 | 0.8770 | 0.1708 | 0.8074 | +0.0000 |
| replogle_k562_essential | `CIPHER_Linear` | 0.4811 | 0.4734 | 0.5003 | 0.4974 | 0.7678 | 0.5138 | +0.8633 |
| replogle_k562_essential | `Hybrid_Raw` | 0.8221 | 0.8250 | 0.8877 | 0.8856 | 0.2364 | 0.8074 | +0.0179 |
| replogle_k562_essential | `Hybrid_ProjShared` | 0.8190 | 0.8213 | 0.8754 | 0.8732 | 0.1772 | 0.8061 | +0.0034 |
| replogle_k562_essential | `NegCtrl_HybridProj_ShuffledGenes` | 0.8220 | 0.8264 | 0.8794 | 0.8778 | 0.1624 | 0.8099 | +0.0010 |
| replogle_k562_essential | `NegCtrl_HybridProj_DiagCov` | 0.8221 | 0.8265 | 0.8796 | 0.8778 | 0.1621 | 0.8099 | +0.0009 |
| replogle_k562_essential | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.8206 | 0.8242 | 0.8782 | 0.8768 | 0.1780 | 0.8081 | +0.0009 |
| replogle_k562_essential | `Transfer_Mean_VarScaled` | 0.8287 | 0.8307 | 0.8995 | 0.8976 | 0.4137 | 0.8134 | +0.1119 |
| replogle_k562_essential | `NegCtrl_TMVS_ShuffledVar` | 0.8343 | 0.8318 | 0.9185 | 0.9120 | 0.5835 | 0.8194 | +0.2045 |
| replogle_k562_essential | `Hybrid_on_TMVS` | 0.8272 | 0.8283 | 0.8962 | 0.8938 | 0.4124 | 0.8121 | +0.1154 |
| replogle_k562_essential | `FCT` | 0.8117 | 0.8170 | 0.8701 | 0.8728 | 0.2991 | 0.7973 | +0.0655 |
| replogle_rpe1 | `Source_Grand_Mean` | 0.5005 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.1715 |
| replogle_rpe1 | `Transfer_Mean` | 0.9274 | 0.9313 | 0.9112 | 0.9130 | 0.3708 | 0.9212 | +0.0000 |
| replogle_rpe1 | `CIPHER_Linear` | 0.5113 | 0.5052 | 0.5202 | 0.5159 | -0.2263 | 0.5227 | +0.1445 |
| replogle_rpe1 | `Hybrid_Raw` | 0.9219 | 0.9243 | 0.9178 | 0.9190 | 0.4093 | 0.9150 | +0.0103 |
| replogle_rpe1 | `Hybrid_ProjShared` | 0.9226 | 0.9253 | 0.9182 | 0.9194 | 0.4081 | 0.9166 | +0.0101 |
| replogle_rpe1 | `NegCtrl_HybridProj_ShuffledGenes` | 0.9259 | 0.9291 | 0.9077 | 0.9086 | 0.3601 | 0.9198 | -0.0022 |
| replogle_rpe1 | `NegCtrl_HybridProj_DiagCov` | 0.9258 | 0.9290 | 0.9080 | 0.9088 | 0.3617 | 0.9198 | -0.0019 |
| replogle_rpe1 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.9271 | 0.9311 | 0.9121 | 0.9141 | 0.3775 | 0.9212 | +0.0011 |
| replogle_rpe1 | `Transfer_Mean_VarScaled` | 0.9275 | 0.9317 | 0.9187 | 0.9220 | 0.4093 | 0.9228 | +0.0144 |
| replogle_rpe1 | `NegCtrl_TMVS_ShuffledVar` | 0.9287 | 0.9330 | 0.9032 | 0.9041 | 0.3192 | 0.9228 | +0.0102 |
| replogle_rpe1 | `Hybrid_on_TMVS` | 0.9226 | 0.9256 | 0.9248 | 0.9274 | 0.4399 | 0.9178 | +0.0245 |
| replogle_rpe1 | `FCT` | 0.9292 | 0.9324 | 0.9139 | 0.9164 | 0.3733 | 0.9246 | +0.0227 |
