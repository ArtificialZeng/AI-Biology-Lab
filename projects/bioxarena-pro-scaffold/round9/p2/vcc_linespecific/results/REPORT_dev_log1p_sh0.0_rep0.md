# Line-Specific Zero-Shot Perturbation Evaluation Report (`dev_log1p_sh0.0_rep0`)

- **Confirmatory Lines (2)**: replogle_k562_essential, replogle_rpe1
- **Secondary Lines (0)**: None
- **Dual Confirmatory Summary (R13, 97.5% Bonferroni CIs)**: **Group A (Hybrid_ProjShared): 未确认 | Group B (Transfer_Mean_VarScaled): 未确认 (Bonferroni 97.5% CIs)**
- **Group A Primary (`Hybrid_ProjShared`)**: **未确认**
- **Group B Primary (`Transfer_Mean_VarScaled`)**: **未确认**
- **R4 Co-fluctuation Structure Claim Allowed**: `目标线对照统计量` (`all beat SrcCorrTgtVar disjoint=False`)
- **R4 Marginal Variance Note**: `增益可由边际方差解释`
- **R14 Mechanism Statement (`Hybrid_on_TMVS` vs `TMVS`)**: `共波动结构在边际方差之外提供了额外信息` (`2/2`)

## Hypothesis Adjudication Summary (Confirmatory Cohort, 97.5% CIs)

| Hypothesis | Group | Pass / Confirmatory | Threshold | Verdict |
|---|---|---:|---|---|
| **H-A1** (Slope 97.5% CI lo > 0.13) | A (`Hybrid_ProjShared`) | 0/2 | >= 2 | **未确认** |
| **H-A2** (Joint-boot `r_cal` 97.5% CI lo > 0, valid>=95%) | A (`Hybrid_ProjShared`) | 2/2 (fixed95: 2/2) | >= 2 | **CONFIRMED** |
| **H-A3** (Matched > ALL source mismatched cov) | A (`Hybrid_ProjShared`) | 2/2 | >= 2 | **CONFIRMED** |
| **H-A4** (Strict AND non-inferiority across 6 conventions) | A (`Hybrid_ProjShared`) | 2/2 inf (ctrl_disc_raw=2/2, ctrl_de50_r=0/2, ctrl_disc_raw_panel=2/2, cent_disc_raw=0/2, cent_de50_r=0/2, cent_disc_raw_panel=0/2; union=2) | <= 0 inf | **未确认** |
| **H-B1** (Slope 97.5% CI lo > 0.13) | B (`Transfer_Mean_VarScaled`) | 0/2 | >= 2 | **未确认** |
| **H-B2** (Joint-boot `r_cal` 97.5% CI lo > 0, valid>=95%) | B (`Transfer_Mean_VarScaled`) | 2/2 | >= 2 | **CONFIRMED** |
| **H-B3** (Matched var > ALL mismatched + shuffled var on slope & `r_cal`) | B (`Transfer_Mean_VarScaled`) | 1/2 | >= 2 | **未确认** |
| **H-B4** (Strict AND non-inferiority across 6 conventions) | B (`Transfer_Mean_VarScaled`) | 1/2 inf (ctrl_disc_raw=1/2, ctrl_de50_r=1/2, ctrl_disc_raw_panel=1/2, cent_disc_raw=0/2, cent_de50_r=0/2, cent_disc_raw_panel=0/2; union=1) | <= 0 inf | **未确认** |

## Main Table: Primary View (A-T2 Own-Target-Zeroed Control Reference, 97.5% CIs, R5/R6/R11/R13)

| Cohort | Target Line | Method | N | Slope OLS (Top50) | 97.5% CI | Calib Resid r | r_cal^2 | r_cal^2 Joint 97.5% CI | Resid r | sd_pred/sd_obs | Deming (delta=1) | Obs~Pred | Ctrl disc_raw | Ctrl de50_r | MSE |
|---|---|---|---:|---:|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|
| Confirmatory | replogle_k562_essential | `Source_Grand_Mean` | 1761 | +0.4566 | [+0.4232, +0.4910] | +0.4838 | 0.2341 | [0.2092, 0.2609] | +0.4749 | 0.9615 | +0.9207 | +0.4939 | 0.5001 | 0.2398 | 0.0066 |
| Confirmatory | replogle_k562_essential | `Transfer_Mean` | 1761 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | 0.0000 | +0.0000 | +nan | 0.7565 | 0.4557 | 0.0127 |
| Confirmatory | replogle_k562_essential | `CIPHER_Linear` | 1761 | +0.6424 | [+0.6079, +0.6786] | +0.5628 | 0.3168 | [0.2911, 0.3439] | +0.5524 | 1.1628 | +1.3110 | +0.4751 | 0.4857 | 0.0520 | 0.0040 |
| Confirmatory | replogle_k562_essential | `Hybrid_Raw` | 1761 | +0.0193 | [+0.0158, +0.0227] | +0.1409 | 0.0199 | [0.0128, 0.0277] | +0.1383 | 0.1399 | +0.0197 | +0.9891 | 0.7577 | 0.4567 | 0.0125 |
| Confirmatory | replogle_k562_essential | `Hybrid_ProjShared` | 1761 | +0.0047 | [+0.0032, +0.0062] | +0.0513 | 0.0026 | [0.0013, 0.0045] | +0.0504 | 0.0926 | +0.0047 | +0.5439 | 0.7500 | 0.4584 | 0.0128 |
| Confirmatory | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 1761 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | 0.0000 | +0.0000 | +nan | 0.7565 | 0.4557 | 0.0127 |
| Confirmatory | replogle_k562_essential | `NegCtrl_HybridProj_ShuffledGenes` | 1761 | -0.0003 | [-0.0021, +0.0015] | -0.0022 | 0.0000 | [0.0000, 0.0003] | -0.0021 | 0.1183 | -0.0003 | -0.0180 | 0.7529 | 0.4518 | 0.0130 |
| Confirmatory | replogle_k562_essential | `NegCtrl_HybridProj_DiagCov` | 1761 | -0.0005 | [-0.0023, +0.0013] | -0.0046 | 0.0000 | [0.0000, 0.0005] | -0.0045 | 0.1116 | -0.0005 | -0.0404 | 0.7532 | 0.4522 | 0.0129 |
| Confirmatory | replogle_k562_essential | `NegCtrl_HybridProj_SrcCorrTgtVar` | 1761 | +0.0007 | [+0.0002, +0.0012] | +0.0201 | 0.0004 | [0.0000, 0.0012] | +0.0197 | 0.0359 | +0.0007 | +0.5487 | 0.7554 | 0.4556 | 0.0127 |
| Confirmatory | replogle_k562_essential | `Transfer_Mean_VarScaled` | 1761 | +0.0621 | [+0.0541, +0.0702] | +0.2054 | 0.0422 | [0.0313, 0.0546] | +0.2016 | 0.3083 | +0.0683 | +0.6539 | 0.7637 | 0.4635 | 0.0108 |
| Confirmatory | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 1761 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | 0.0000 | +0.0000 | +nan | 0.7565 | 0.4557 | 0.0127 |
| Confirmatory | replogle_k562_essential | `NegCtrl_TMVS_ShuffledVar` | 1761 | +0.1664 | [+0.1520, +0.1808] | +0.3381 | 0.1143 | [0.0934, 0.1367] | +0.3319 | 0.5013 | +0.2122 | +0.6621 | 0.7692 | 0.4495 | 0.0098 |
| Confirmatory | replogle_k562_essential | `Hybrid_on_TMVS` | 1761 | +0.0668 | [+0.0584, +0.0752] | +0.2109 | 0.0445 | [0.0335, 0.0570] | +0.2070 | 0.3227 | +0.0742 | +0.6413 | 0.7566 | 0.4660 | 0.0109 |
| Confirmatory | replogle_k562_essential | `FCT` | 1761 | +0.0382 | [+0.0319, +0.0447] | +0.1231 | 0.0152 | [0.0105, 0.0205] | +0.1208 | 0.3163 | +0.0424 | +0.3820 | 0.7507 | 0.4702 | 0.0132 |
| Confirmatory | replogle_rpe1 | `Source_Grand_Mean` | 1761 | -0.0456 | [-0.0535, -0.0374] | -0.1384 | 0.0192 | [0.0132, 0.0261] | -0.1359 | 0.3356 | -0.0513 | -0.4050 | 0.5002 | 0.4963 | 0.0116 |
| Confirmatory | replogle_rpe1 | `Transfer_Mean` | 1761 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | 0.0000 | +0.0000 | +nan | 0.8434 | 0.4321 | 0.0127 |
| Confirmatory | replogle_rpe1 | `CIPHER_Linear` | 1761 | -0.0818 | [-0.0905, -0.0730] | -0.2264 | 0.0512 | [0.0415, 0.0616] | -0.2222 | 0.3682 | -0.0938 | -0.6036 | 0.5150 | 0.2506 | 0.0124 |
| Confirmatory | replogle_rpe1 | `Hybrid_Raw` | 1761 | +0.0096 | [+0.0069, +0.0122] | +0.1458 | 0.0213 | [0.0101, 0.0362] | +0.1432 | 0.0673 | +0.0097 | +2.1271 | 0.8267 | 0.4610 | 0.0125 |
| Confirmatory | replogle_rpe1 | `Hybrid_ProjShared` | 1761 | +0.0083 | [+0.0061, +0.0104] | +0.1466 | 0.0215 | [0.0111, 0.0352] | +0.1440 | 0.0577 | +0.0083 | +2.4929 | 0.8293 | 0.4585 | 0.0125 |
| Confirmatory | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 1761 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | 0.0000 | +0.0000 | +nan | 0.8434 | 0.4321 | 0.0127 |
| Confirmatory | replogle_rpe1 | `NegCtrl_HybridProj_ShuffledGenes` | 1761 | -0.0011 | [-0.0015, -0.0007] | -0.0398 | 0.0016 | [0.0007, 0.0029] | -0.0391 | 0.0291 | -0.0011 | -1.3433 | 0.8439 | 0.4239 | 0.0128 |
| Confirmatory | replogle_rpe1 | `NegCtrl_HybridProj_DiagCov` | 1761 | -0.0011 | [-0.0015, -0.0008] | -0.0462 | 0.0021 | [0.0010, 0.0038] | -0.0454 | 0.0253 | -0.0011 | -1.7902 | 0.8435 | 0.4249 | 0.0127 |
| Confirmatory | replogle_rpe1 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 1761 | +0.0015 | [+0.0013, +0.0017] | +0.1337 | 0.0179 | [0.0130, 0.0234] | +0.1312 | 0.0113 | +0.0015 | +11.6556 | 0.8409 | 0.4367 | 0.0126 |
| Confirmatory | replogle_rpe1 | `Transfer_Mean_VarScaled` | 1761 | +0.0222 | [+0.0204, +0.0240] | +0.2033 | 0.0413 | [0.0353, 0.0474] | +0.1996 | 0.1112 | +0.0225 | +1.7955 | 0.8386 | 0.4280 | 0.0126 |
| Confirmatory | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 1761 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | 0.0000 | +0.0000 | +nan | 0.8434 | 0.4321 | 0.0127 |
| Confirmatory | replogle_rpe1 | `NegCtrl_TMVS_ShuffledVar` | 1761 | -0.0075 | [-0.0091, -0.0059] | -0.0534 | 0.0029 | [0.0017, 0.0044] | -0.0525 | 0.1437 | -0.0077 | -0.3651 | 0.8494 | 0.4228 | 0.0130 |
| Confirmatory | replogle_rpe1 | `Hybrid_on_TMVS` | 1761 | +0.0305 | [+0.0281, +0.0330] | +0.2467 | 0.0608 | [0.0520, 0.0702] | +0.2422 | 0.1260 | +0.0310 | +1.9223 | 0.8270 | 0.4503 | 0.0125 |
| Confirmatory | replogle_rpe1 | `FCT` | 1761 | +0.0078 | [+0.0064, +0.0092] | +0.0684 | 0.0047 | [0.0032, 0.0065] | +0.0672 | 0.1158 | +0.0079 | +0.5804 | 0.8443 | 0.4189 | 0.0130 |

## S6: Split-Half Noise Ceiling (`resp_top50` vs Unselected `all_nontarget`, Ceiling Selection Bias)

| Target Line | N half-valid | rho_AB (top50) | Ceil r (top50) | Rep Slope (top50) | rho_AB (all_nontarget) | Ceil r (all_nontarget) | TMVS r_cal (all) | Hybrid r_cal (all) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| replogle_k562_essential | 1761 | 0.9294 | 0.9815 | 0.9294 | 0.7870 | 0.9385 | +0.4596 | +0.0233 |
| replogle_rpe1 | 1761 | 0.9304 | 0.9818 | 0.9304 | 0.6349 | 0.8813 | +0.0831 | +0.1289 |

## S5: Expressed Genes Sensitivity (Drop Bottom 25% Target Control Variance Genes)

| Target Line | Genes Kept | rho_AB (S5) | Method | Slope OLS (S5) | Calib Resid r (S5) | r_cal^2 (S5) | Joint 97.5% CI |
|---|---:|---:|---|---:|---:|---:|---|
| replogle_k562_essential | 5419/7225 | 0.9287 | `Hybrid_ProjShared` | +0.0047 | +0.0511 | 0.0026 | [+0.0352, +0.0673] |
| replogle_k562_essential | 5419/7225 | 0.9287 | `Transfer_Mean_VarScaled` | +0.0658 | +0.2160 | 0.0467 | [+0.1871, +0.2443] |
| replogle_k562_essential | 5419/7225 | 0.9287 | `Hybrid_on_TMVS` | +0.0705 | +0.2209 | 0.0488 | [+0.1930, +0.2486] |
| replogle_k562_essential | 5419/7225 | 0.9287 | `FCT` | +0.0396 | +0.1269 | 0.0161 | [+0.1064, +0.1469] |
| replogle_k562_essential | 5419/7225 | 0.9287 | `NegCtrl_HybridProj_SrcCorrTgtVar` | +0.0006 | +0.0178 | 0.0003 | [+0.0039, +0.0321] |
| replogle_k562_essential | 5419/7225 | 0.9287 | `NegCtrl_TMVS_ShuffledVar` | +0.1877 | +0.4115 | 0.1693 | [+0.3933, +0.4308] |
| replogle_rpe1 | 5419/7225 | 0.9296 | `Hybrid_ProjShared` | +0.0085 | +0.1483 | 0.0220 | [+0.1065, +0.1895] |
| replogle_rpe1 | 5419/7225 | 0.9296 | `Transfer_Mean_VarScaled` | +0.0234 | +0.2177 | 0.0474 | [+0.2020, +0.2325] |
| replogle_rpe1 | 5419/7225 | 0.9296 | `Hybrid_on_TMVS` | +0.0319 | +0.2604 | 0.0678 | [+0.2410, +0.2792] |
| replogle_rpe1 | 5419/7225 | 0.9296 | `FCT` | +0.0090 | +0.0791 | 0.0063 | [+0.0678, +0.0910] |
| replogle_rpe1 | 5419/7225 | 0.9296 | `NegCtrl_HybridProj_SrcCorrTgtVar` | +0.0015 | +0.1376 | 0.0189 | [+0.1180, +0.1568] |
| replogle_rpe1 | 5419/7225 | 0.9296 | `NegCtrl_TMVS_ShuffledVar` | -0.0093 | -0.0745 | 0.0056 | [-0.0851, -0.0635] |

## Appendix: Alternative H-A4/H-B4 Discrimination & Centroid Conventions (R5)

| Target Line | Method | Ctrl disc_raw (A-T2) | Ctrl disc_raw_panel (A-T3) | Cent disc_raw | Cent disc_raw_panel | Cent de50_r | pds_cosine_panel | Slope OLS (All) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| replogle_k562_essential | `Source_Grand_Mean` | 0.5001 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.5681 |
| replogle_k562_essential | `Transfer_Mean` | 0.7565 | 0.7664 | 0.8205 | 0.8226 | 0.2794 | 0.7614 | +0.0000 |
| replogle_k562_essential | `CIPHER_Linear` | 0.4857 | 0.4811 | 0.5001 | 0.4984 | 0.8039 | 0.5078 | +0.8342 |
| replogle_k562_essential | `Hybrid_Raw` | 0.7577 | 0.7668 | 0.8272 | 0.8283 | 0.3256 | 0.7614 | +0.0197 |
| replogle_k562_essential | `Hybrid_ProjShared` | 0.7500 | 0.7597 | 0.8192 | 0.8214 | 0.2794 | 0.7553 | +0.0027 |
| replogle_k562_essential | `NegCtrl_HybridProj_ShuffledGenes` | 0.7529 | 0.7651 | 0.8213 | 0.8239 | 0.2691 | 0.7604 | +0.0010 |
| replogle_k562_essential | `NegCtrl_HybridProj_DiagCov` | 0.7532 | 0.7650 | 0.8213 | 0.8239 | 0.2691 | 0.7603 | +0.0009 |
| replogle_k562_essential | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.7554 | 0.7656 | 0.8203 | 0.8229 | 0.2834 | 0.7607 | +0.0006 |
| replogle_k562_essential | `Transfer_Mean_VarScaled` | 0.7637 | 0.7728 | 0.8359 | 0.8379 | 0.4622 | 0.7665 | +0.1026 |
| replogle_k562_essential | `NegCtrl_TMVS_ShuffledVar` | 0.7692 | 0.7746 | 0.8476 | 0.8448 | 0.6117 | 0.7714 | +0.1890 |
| replogle_k562_essential | `Hybrid_on_TMVS` | 0.7566 | 0.7652 | 0.8343 | 0.8357 | 0.4577 | 0.7598 | +0.1052 |
| replogle_k562_essential | `FCT` | 0.7507 | 0.7630 | 0.8163 | 0.8242 | 0.3762 | 0.7548 | +0.0493 |
| replogle_rpe1 | `Source_Grand_Mean` | 0.5002 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.1875 |
| replogle_rpe1 | `Transfer_Mean` | 0.8434 | 0.8556 | 0.8470 | 0.8536 | 0.3400 | 0.8470 | +0.0000 |
| replogle_rpe1 | `CIPHER_Linear` | 0.5150 | 0.5132 | 0.5137 | 0.5134 | -0.1145 | 0.5215 | +0.1709 |
| replogle_rpe1 | `Hybrid_Raw` | 0.8267 | 0.8363 | 0.8603 | 0.8680 | 0.3795 | 0.8289 | +0.0086 |
| replogle_rpe1 | `Hybrid_ProjShared` | 0.8293 | 0.8392 | 0.8579 | 0.8663 | 0.3760 | 0.8325 | +0.0077 |
| replogle_rpe1 | `NegCtrl_HybridProj_ShuffledGenes` | 0.8439 | 0.8554 | 0.8440 | 0.8502 | 0.3317 | 0.8470 | -0.0018 |
| replogle_rpe1 | `NegCtrl_HybridProj_DiagCov` | 0.8435 | 0.8549 | 0.8441 | 0.8504 | 0.3329 | 0.8466 | -0.0016 |
| replogle_rpe1 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.8409 | 0.8538 | 0.8478 | 0.8545 | 0.3450 | 0.8451 | +0.0008 |
| replogle_rpe1 | `Transfer_Mean_VarScaled` | 0.8386 | 0.8517 | 0.8523 | 0.8598 | 0.3597 | 0.8443 | +0.0091 |
| replogle_rpe1 | `NegCtrl_TMVS_ShuffledVar` | 0.8494 | 0.8595 | 0.8391 | 0.8442 | 0.3111 | 0.8517 | +0.0096 |
| replogle_rpe1 | `Hybrid_on_TMVS` | 0.8270 | 0.8373 | 0.8615 | 0.8704 | 0.3890 | 0.8315 | +0.0167 |
| replogle_rpe1 | `FCT` | 0.8443 | 0.8531 | 0.8499 | 0.8568 | 0.3385 | 0.8471 | +0.0178 |
