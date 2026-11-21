# Line-Specific Zero-Shot Perturbation Evaluation Report (`dev_log1p_sh0.0_rep1`)

- **Code Freeze Commit (C1)**: `UNFROZEN_DEV` | **Bootstrap Seed (C2)**: `20261121`
- **Confirmatory Lines (2)**: replogle_k562_essential, replogle_rpe1
- **Secondary Lines (0)**: None
- **Dual Confirmatory Summary (R13, 97.5% Bonferroni CIs)**: **Group A (Hybrid_ProjShared): 未确认 | Group B (Transfer_Mean_VarScaled): 未确认 (Bonferroni 97.5% CIs)**
- **Group A Primary (`Hybrid_ProjShared`)**: **未确认**
- **Group B Primary (`Transfer_Mean_VarScaled`)**: **未确认**
- **R17 Primary Correlation Metric**: `ls_resid_corr_orth_cal_top50` (amplitude-orthogonalized; raw `r_cal` reported as secondary)
- **R18 Global-Scale Scalar Baseline**: 单个标量 Transfer_Mean_GlobalScale 即可得到原始残差相关 [+0.303, +0.555]
- **R4 Co-fluctuation Structure Claim Allowed**: `目标线对照统计量` (`all beat SrcCorrTgtVar disjoint=False`)
- **R4 Marginal Variance Note**: `增益可由边际方差解释`
- **R14 Mechanism Statement (`Hybrid_on_TMVS` vs `TMVS` on `r_orth_cal`, 95% CI)**: `共波动结构在边际方差之外提供了额外信息` (`95%: 2/2`, `97.5%: 2/2`)
- **R20 Prior-Art Novelty Scope**: FCT is equivalent to VCC 2026 CellularAI (mu_t * 2^log2FC); CASTER performed supervised control-variance scaling. TMVS is framed strictly as a parameter-free variance-ratio variant of expression-aware transfer (not a new method); contribution is the pre-registered 3-layer decomposition of recoverable zero-shot line-specific residuals into (1) line-wide amplitude, (2) gene marginal variance, and (3) off-diagonal co-fluctuations.

## Hypothesis Adjudication Summary (Confirmatory Cohort, 97.5% CIs)

| Hypothesis | Group | Pass / Confirmatory | Threshold | Verdict |
|---|---|---:|---|---|
| **H-A1** (Slope 97.5% CI lo > 0.13) | A (`Hybrid_ProjShared`) | 0/2 | >= 4 | **未确认** |
| **H-A2** (Joint-boot `r_orth_cal` 97.5% CI lo > 0, valid>=95%, rho_AB_orth>0) | A (`Hybrid_ProjShared`) | 2/2 (fixed95: 2/2) | >= 4 | **未确认** |
| **H-A3** (Matched > ALL source mismatched cov on slope AND `r_orth_cal`) | A (`Hybrid_ProjShared`) | 2/2 | >= 4 | **未确认** |
| **H-A4** (Strict AND non-inferiority across 6 conventions) | A (`Hybrid_ProjShared`) | 1/2 inf (ctrl_disc_raw=1/2, ctrl_de50_r=0/2, ctrl_disc_raw_panel=1/2, cent_disc_raw=1/2, cent_de50_r=0/2, cent_disc_raw_panel=1/2; union=2) | <= 1 inf | **CONFIRMED** |
| **H-B1** (Slope 97.5% CI lo > 0.13) | B (`Transfer_Mean_VarScaled`) | 0/2 | >= 4 | **未确认** |
| **H-B2** (Joint-boot `r_orth_cal` 97.5% CI lo > 0, valid>=95%, rho_AB_orth>0) | B (`Transfer_Mean_VarScaled`) | 2/2 | >= 4 | **未确认** |
| **H-B3** (Matched `r_orth_cal` lo > max(ALL mismatched vars CI hi, 100-perm p97.5)) | B (`Transfer_Mean_VarScaled`) | 1/2 | >= 4 | **未确认** |
| **H-B4** (Strict AND non-inferiority across 6 conventions) | B (`Transfer_Mean_VarScaled`) | 0/2 inf (ctrl_disc_raw=0/2, ctrl_de50_r=0/2, ctrl_disc_raw_panel=0/2, cent_disc_raw=0/2, cent_de50_r=0/2, cent_disc_raw_panel=0/2; union=0) | <= 1 inf | **CONFIRMED** |

## Main Table: Primary View (A-T2 Own-Target-Zeroed Control Reference, 97.5% CIs, R5/R6/R11/R13/R17/R18)

| Cohort | Target Line | Method | N | Slope OLS (Top50) | 97.5% CI | r_orth_cal | r_orth_cal^2 | r_orth_cal^2 Joint 97.5% CI | Resid r_orth | r_cal (raw) | Resid r (raw) | sd_pred/sd_obs | Deming (delta=1) | Obs~Pred | Ctrl disc_raw | Ctrl de50_r | MSE |
|---|---|---|---:|---:|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Confirmatory | replogle_k562_essential | `Source_Grand_Mean` | 832 | +0.4606 | [+0.4102, +0.5087] | -0.1212 | 0.0147 | [0.0093, 0.0210] | -0.1189 | +0.4695 | +0.4634 | 0.9938 | +0.9866 | +0.4663 | 0.5002 | 0.3439 | 0.0079 |
| Confirmatory | replogle_k562_essential | `Transfer_Mean` | 832 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.8201 | 0.5996 | 0.0150 |
| Confirmatory | replogle_k562_essential | `CIPHER_Linear` | 832 | +0.6998 | [+0.6477, +0.7486] | +0.0576 | 0.0033 | [0.0007, 0.0080] | +0.0565 | +0.5631 | +0.5559 | 1.2590 | +1.5018 | +0.4415 | 0.4811 | 0.0683 | 0.0051 |
| Confirmatory | replogle_k562_essential | `Hybrid_Raw` | 832 | +0.0200 | [+0.0152, +0.0247] | +0.0619 | 0.0038 | [0.0013, 0.0077] | +0.0607 | +0.1379 | +0.1361 | 0.1467 | +0.0204 | +0.9281 | 0.8221 | 0.6022 | 0.0149 |
| Confirmatory | replogle_k562_essential | `Hybrid_ProjShared` | 832 | +0.0050 | [+0.0028, +0.0072] | +0.0838 | 0.0070 | [0.0035, 0.0115] | +0.0822 | +0.0523 | +0.0516 | 0.0964 | +0.0050 | +0.5356 | 0.8190 | 0.6021 | 0.0152 |
| Confirmatory | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 832 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.8201 | 0.5996 | 0.0150 |
| Confirmatory | replogle_k562_essential | `NegCtrl_HybridProj_ShuffledGenes` | 832 | -0.0020 | [-0.0045, +0.0006] | -0.0018 | 0.0000 | [0.0000, 0.0007] | -0.0018 | -0.0160 | -0.0158 | 0.1274 | -0.0020 | -0.1243 | 0.8220 | 0.5951 | 0.0154 |
| Confirmatory | replogle_k562_essential | `NegCtrl_HybridProj_DiagCov` | 832 | -0.0023 | [-0.0048, +0.0003] | -0.0049 | 0.0000 | [0.0000, 0.0008] | -0.0048 | -0.0191 | -0.0188 | 0.1218 | -0.0023 | -0.1545 | 0.8221 | 0.5954 | 0.0153 |
| Confirmatory | replogle_k562_essential | `NegCtrl_HybridProj_SrcCorrTgtVar` | 832 | +0.0011 | [+0.0004, +0.0018] | +0.0235 | 0.0006 | [0.0000, 0.0020] | +0.0231 | +0.0292 | +0.0288 | 0.0380 | +0.0011 | +0.7581 | 0.8206 | 0.5994 | 0.0150 |
| Confirmatory | replogle_k562_essential | `Transfer_Mean_VarScaled` | 832 | +0.0737 | [+0.0626, +0.0848] | +0.1325 | 0.0176 | [0.0112, 0.0249] | +0.1300 | +0.2198 | +0.2170 | 0.3396 | +0.0827 | +0.6390 | 0.8287 | 0.6116 | 0.0126 |
| Confirmatory | replogle_k562_essential | `Transfer_Mean_GlobalScale` | 832 | +0.0308 | [+0.0286, +0.0330] | +nan | nan | [nan, nan] | +nan | +0.5618 | +0.5546 | 0.0556 | +0.0309 | +9.9793 | 0.8201 | 0.5996 | 0.0139 |
| Confirmatory | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 832 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.8201 | 0.5996 | 0.0150 |
| Confirmatory | replogle_k562_essential | `NegCtrl_TMVS_ShuffledVar` | 832 | +0.1715 | [+0.1508, +0.1916] | +0.0843 | 0.0071 | [0.0046, 0.0102] | +0.0827 | +0.3169 | +0.3128 | 0.5483 | +0.2320 | +0.5705 | 0.8343 | 0.5883 | 0.0114 |
| Confirmatory | replogle_k562_essential | `Hybrid_on_TMVS` | 832 | +0.0787 | [+0.0671, +0.0899] | +0.1500 | 0.0225 | [0.0148, 0.0316] | +0.1471 | +0.2254 | +0.2225 | 0.3536 | +0.0892 | +0.6293 | 0.8272 | 0.6136 | 0.0127 |
| Confirmatory | replogle_k562_essential | `FCT` | 832 | +0.0462 | [+0.0379, +0.0540] | +0.2357 | 0.0556 | [0.0483, 0.0634] | +0.2313 | +0.1394 | +0.1376 | 0.3357 | +0.0519 | +0.4099 | 0.8117 | 0.6222 | 0.0153 |
| Confirmatory | replogle_rpe1 | `Source_Grand_Mean` | 832 | -0.0500 | [-0.0625, -0.0370] | +0.5072 | 0.2572 | [0.2388, 0.2753] | +0.4999 | -0.1384 | -0.1364 | 0.3665 | -0.0575 | -0.3722 | 0.5005 | 0.6210 | 0.0143 |
| Confirmatory | replogle_rpe1 | `Transfer_Mean` | 832 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.9274 | 0.5548 | 0.0150 |
| Confirmatory | replogle_rpe1 | `CIPHER_Linear` | 832 | -0.1107 | [-0.1242, -0.0969] | +0.1929 | 0.0372 | [0.0196, 0.0603] | +0.1901 | -0.2652 | -0.2613 | 0.4234 | -0.1325 | -0.6172 | 0.5113 | 0.3256 | 0.0160 |
| Confirmatory | replogle_rpe1 | `Hybrid_Raw` | 832 | +0.0091 | [+0.0040, +0.0140] | +0.1584 | 0.0251 | [0.0105, 0.0462] | +0.1561 | +0.0951 | +0.0938 | 0.0974 | +0.0092 | +0.9625 | 0.9219 | 0.5785 | 0.0149 |
| Confirmatory | replogle_rpe1 | `Hybrid_ProjShared` | 832 | +0.0085 | [+0.0045, +0.0124] | +0.1527 | 0.0233 | [0.0102, 0.0421] | +0.1505 | +0.1054 | +0.1038 | 0.0823 | +0.0086 | +1.2622 | 0.9226 | 0.5788 | 0.0149 |
| Confirmatory | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 832 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.9274 | 0.5548 | 0.0150 |
| Confirmatory | replogle_rpe1 | `NegCtrl_HybridProj_ShuffledGenes` | 832 | -0.0011 | [-0.0018, -0.0004] | -0.0656 | 0.0043 | [0.0022, 0.0069] | -0.0647 | -0.0288 | -0.0283 | 0.0395 | -0.0011 | -0.7171 | 0.9259 | 0.5435 | 0.0153 |
| Confirmatory | replogle_rpe1 | `NegCtrl_HybridProj_DiagCov` | 832 | -0.0010 | [-0.0017, -0.0003] | -0.0713 | 0.0051 | [0.0027, 0.0081] | -0.0703 | -0.0286 | -0.0282 | 0.0348 | -0.0010 | -0.8089 | 0.9258 | 0.5455 | 0.0151 |
| Confirmatory | replogle_rpe1 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 832 | +0.0018 | [+0.0015, +0.0021] | +0.1211 | 0.0147 | [0.0090, 0.0215] | +0.1193 | +0.1144 | +0.1128 | 0.0159 | +0.0018 | +7.1078 | 0.9271 | 0.5597 | 0.0150 |
| Confirmatory | replogle_rpe1 | `Transfer_Mean_VarScaled` | 832 | +0.0267 | [+0.0244, +0.0291] | +0.0979 | 0.0096 | [0.0070, 0.0128] | +0.0965 | +0.2145 | +0.2114 | 0.1262 | +0.0271 | +1.6753 | 0.9275 | 0.5544 | 0.0149 |
| Confirmatory | replogle_rpe1 | `Transfer_Mean_GlobalScale` | 832 | +0.0058 | [+0.0052, +0.0063] | +nan | nan | [nan, nan] | +nan | +0.3078 | +0.3034 | 0.0190 | +0.0058 | +15.9564 | 0.9274 | 0.5548 | 0.0152 |
| Confirmatory | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 832 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.9274 | 0.5548 | 0.0150 |
| Confirmatory | replogle_rpe1 | `NegCtrl_TMVS_ShuffledVar` | 832 | -0.0099 | [-0.0122, -0.0074] | +0.0706 | 0.0050 | [0.0034, 0.0068] | +0.0696 | -0.0599 | -0.0591 | 0.1673 | -0.0102 | -0.3531 | 0.9287 | 0.5355 | 0.0155 |
| Confirmatory | replogle_rpe1 | `Hybrid_on_TMVS` | 832 | +0.0352 | [+0.0310, +0.0393] | +0.1664 | 0.0277 | [0.0189, 0.0379] | +0.1640 | +0.2351 | +0.2317 | 0.1521 | +0.0360 | +1.5236 | 0.9226 | 0.5753 | 0.0147 |
| Confirmatory | replogle_rpe1 | `FCT` | 832 | +0.0111 | [+0.0088, +0.0134] | +0.0532 | 0.0028 | [0.0012, 0.0051] | +0.0524 | +0.0903 | +0.0890 | 0.1250 | +0.0113 | +0.7122 | 0.9292 | 0.5462 | 0.0154 |

## R18: TMVS Shuffled-Variance 100-Permutation Null vs H-B3 (97.5% CI Threshold)

| Target Line | TMVS r_orth_cal | TMVS CI lo | Perm100 p2.5 | Perm100 median | Perm100 p97.5 | Perm100 max | Max Mismatch CI hi | H-B3 Pass |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| replogle_k562_essential | +0.1325 | +0.1058 | +0.0495 | +0.0913 | +0.1266 | +0.1296 | +0.0000 | False |
| replogle_rpe1 | +0.0979 | +0.0837 | -0.0048 | +0.0324 | +0.0571 | +0.0655 | +0.0000 | True |

## R19: Descriptive 3-Layer Residual Variance Decomposition (`R = Delta_t - Delta_src`, PREREG_A v4 §5.4)

| Target Line | Mask | N | Elig Genes | beta_ols_int [95% CI] | beta_mult | L1 Line Gain [95% CI] | L2 Per-Gene Gain [95% CI] | L3 Remaining [95% CI] |
|---|---|---:|---:|---|---:|---|---|---|
| replogle_k562_essential | `resp_top50` | 832 | 1217 (C_dist>=5) | 0.5611 [0.5404, 0.5822] | 0.5789 | 0.3076 [0.2778, 0.3375] | 0.1425 [0.1308, 0.1520] | 0.5500 [0.5241, 0.5785] |
| replogle_k562_essential | `all_nontarget` | 832 | 7225 (C_dist>=5) | 0.1989 [0.1903, 0.2075] | 0.1994 | 0.6964 [0.6838, 0.7090] | 0.0259 [0.0239, 0.0279] | 0.2777 [0.2663, 0.2890] |
| replogle_rpe1 | `resp_top50` | 832 | 629 (C_dist>=5) | 1.7340 [1.6604, 1.8133] | 2.0716 | 0.0920 [0.0782, 0.1068] | 0.1195 [0.1059, 0.1294] | 0.7885 [0.7727, 0.8082] |
| replogle_rpe1 | `all_nontarget` | 832 | 7225 (C_dist>=5) | 0.6230 [0.5993, 0.6450] | 0.6397 | 0.0493 [0.0429, 0.0563] | 0.1103 [0.1043, 0.1166] | 0.8404 [0.8332, 0.8473] |

## S6: Split-Half Noise Ceiling (`resp_top50` vs Unselected `all_nontarget`, Amplitude-Orthogonalized & Raw)

| Target Line | N half-valid | rho_AB_orth (top50) | Ceil r_orth (top50) | rho_AB (top50) | Rep Slope (top50) | rho_AB_orth (all) | rho_AB (all) | TMVS r_orth_cal (all) | Hybrid r_orth_cal (all) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| replogle_k562_essential | 832 | 0.9273 | 0.9810 | 0.9503 | 0.9503 | 0.5418 | 0.8421 | +0.1193 | +0.0509 |
| replogle_rpe1 | 832 | 0.9447 | 0.9857 | 0.9441 | 0.9441 | 0.6965 | 0.7160 | +0.1502 | +0.1268 |

## S5: Expressed Genes Sensitivity (Drop Bottom 25% Target Control Variance Genes)

| Target Line | Genes Kept | rho_AB_orth (S5) | Method | Slope OLS (S5) | r_orth_cal (S5) | r_orth_cal^2 (S5) | Joint 97.5% CI | r_cal raw (S5) |
|---|---:|---:|---|---:|---:|---:|---|---:|
| replogle_k562_essential | 5419/7225 | 0.9264 | `Hybrid_ProjShared` | +0.0052 | +0.0861 | 0.0074 | [+0.0619, +0.1096] | +0.0549 |
| replogle_k562_essential | 5419/7225 | 0.9264 | `Transfer_Mean_VarScaled` | +0.0777 | +0.1335 | 0.0178 | [+0.1067, +0.1591] | +0.2304 |
| replogle_k562_essential | 5419/7225 | 0.9264 | `Transfer_Mean_GlobalScale` | +0.0306 | +nan | nan | [nan, nan] | +0.5621 |
| replogle_k562_essential | 5419/7225 | 0.9264 | `Hybrid_on_TMVS` | +0.0829 | +0.1514 | 0.0229 | [+0.1230, +0.1792] | +0.2360 |
| replogle_k562_essential | 5419/7225 | 0.9264 | `FCT` | +0.0478 | +0.2371 | 0.0562 | [+0.2210, +0.2527] | +0.1436 |
| replogle_k562_essential | 5419/7225 | 0.9264 | `NegCtrl_HybridProj_SrcCorrTgtVar` | +0.0010 | +0.0203 | 0.0004 | [+0.0021, +0.0402] | +0.0256 |
| replogle_k562_essential | 5419/7225 | 0.9264 | `NegCtrl_TMVS_ShuffledVar` | +0.1982 | +0.0952 | 0.0091 | [+0.0807, +0.1099] | +0.4005 |
| replogle_rpe1 | 5419/7225 | 0.9441 | `Hybrid_ProjShared` | +0.0087 | +0.1569 | 0.0246 | [+0.1046, +0.2089] | +0.1063 |
| replogle_rpe1 | 5419/7225 | 0.9441 | `Transfer_Mean_VarScaled` | +0.0285 | +0.0968 | 0.0094 | [+0.0834, +0.1105] | +0.2327 |
| replogle_rpe1 | 5419/7225 | 0.9441 | `Transfer_Mean_GlobalScale` | +0.0057 | +nan | nan | [nan, nan] | +0.3132 |
| replogle_rpe1 | 5419/7225 | 0.9441 | `Hybrid_on_TMVS` | +0.0372 | +0.1695 | 0.0287 | [+0.1395, +0.1991] | +0.2501 |
| replogle_rpe1 | 5419/7225 | 0.9441 | `FCT` | +0.0133 | +0.0581 | 0.0034 | [+0.0404, +0.0756] | +0.1079 |
| replogle_rpe1 | 5419/7225 | 0.9441 | `NegCtrl_HybridProj_SrcCorrTgtVar` | +0.0019 | +0.1289 | 0.0166 | [+0.1032, +0.1534] | +0.1198 |
| replogle_rpe1 | 5419/7225 | 0.9441 | `NegCtrl_TMVS_ShuffledVar` | -0.0123 | +0.0815 | 0.0066 | [+0.0702, +0.0922] | -0.0853 |

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
| replogle_k562_essential | `Transfer_Mean_GlobalScale` | 0.8201 | 0.8237 | 0.8842 | 0.8828 | 0.2046 | 0.8074 | +0.0382 |
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
| replogle_rpe1 | `Transfer_Mean_GlobalScale` | 0.9274 | 0.9313 | 0.9141 | 0.9158 | 0.3789 | 0.9212 | -0.0060 |
| replogle_rpe1 | `NegCtrl_TMVS_ShuffledVar` | 0.9287 | 0.9330 | 0.9032 | 0.9041 | 0.3192 | 0.9228 | +0.0102 |
| replogle_rpe1 | `Hybrid_on_TMVS` | 0.9226 | 0.9256 | 0.9248 | 0.9274 | 0.4399 | 0.9178 | +0.0245 |
| replogle_rpe1 | `FCT` | 0.9292 | 0.9324 | 0.9139 | 0.9164 | 0.3733 | 0.9246 | +0.0227 |
