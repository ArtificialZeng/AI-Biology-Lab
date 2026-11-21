# Line-Specific Zero-Shot Perturbation Evaluation Report (`dev_log1p_sh0.0_rep0`)

- **Code Freeze Commit (C1)**: `UNFROZEN_DEV` | **Bootstrap Seed (C2)**: `20261121`
- **Confirmatory Lines (2)**: replogle_k562_essential, replogle_rpe1
- **Secondary Lines (0)**: None
- **Dual Confirmatory Summary (R13, 97.5% Bonferroni CIs)**: **Group A (Hybrid_ProjShared): 未确认 | Group B (Transfer_Mean_VarScaled): 未确认 (Bonferroni 97.5% CIs)**
- **Group A Primary (`Hybrid_ProjShared`)**: **未确认**
- **Group B Primary (`Transfer_Mean_VarScaled`)**: **未确认**
- **R17 Primary Correlation Metric**: `ls_resid_corr_orth_cal_top50` (amplitude-orthogonalized; raw `r_cal` reported as secondary)
- **R18 Global-Scale Scalar Baseline**: 单个标量 Transfer_Mean_GlobalScale 即可得到原始残差相关 [+0.263, +0.552]
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
| **H-A4** (Strict AND non-inferiority across 6 conventions) | A (`Hybrid_ProjShared`) | 2/2 inf (ctrl_disc_raw=2/2, ctrl_de50_r=0/2, ctrl_disc_raw_panel=2/2, cent_disc_raw=0/2, cent_de50_r=0/2, cent_disc_raw_panel=0/2; union=2) | <= 1 inf | **未确认** |
| **H-B1** (Slope 97.5% CI lo > 0.13) | B (`Transfer_Mean_VarScaled`) | 0/2 | >= 4 | **未确认** |
| **H-B2** (Joint-boot `r_orth_cal` 97.5% CI lo > 0, valid>=95%, rho_AB_orth>0) | B (`Transfer_Mean_VarScaled`) | 2/2 | >= 4 | **未确认** |
| **H-B3** (Matched `r_orth_cal` lo > max(ALL mismatched vars CI hi, 100-perm p97.5)) | B (`Transfer_Mean_VarScaled`) | 1/2 | >= 4 | **未确认** |
| **H-B4** (Strict AND non-inferiority across 6 conventions) | B (`Transfer_Mean_VarScaled`) | 1/2 inf (ctrl_disc_raw=1/2, ctrl_de50_r=1/2, ctrl_disc_raw_panel=1/2, cent_disc_raw=0/2, cent_de50_r=0/2, cent_disc_raw_panel=0/2; union=1) | <= 1 inf | **CONFIRMED** |

## Main Table: Primary View (A-T2 Own-Target-Zeroed Control Reference, 97.5% CIs, R5/R6/R11/R13/R17/R18)

| Cohort | Target Line | Method | N | Slope OLS (Top50) | 97.5% CI | r_orth_cal | r_orth_cal^2 | r_orth_cal^2 Joint 97.5% CI | Resid r_orth | r_cal (raw) | Resid r (raw) | sd_pred/sd_obs | Deming (delta=1) | Obs~Pred | Ctrl disc_raw | Ctrl de50_r | MSE |
|---|---|---|---:|---:|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Confirmatory | replogle_k562_essential | `Source_Grand_Mean` | 1761 | +0.4566 | [+0.4237, +0.4907] | -0.0697 | 0.0049 | [0.0026, 0.0075] | -0.0678 | +0.4838 | +0.4749 | 0.9615 | +0.9207 | +0.4939 | 0.5001 | 0.2398 | 0.0066 |
| Confirmatory | replogle_k562_essential | `Transfer_Mean` | 1761 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.7565 | 0.4557 | 0.0127 |
| Confirmatory | replogle_k562_essential | `CIPHER_Linear` | 1761 | +0.6424 | [+0.6067, +0.6783] | +0.0536 | 0.0029 | [0.0008, 0.0062] | +0.0521 | +0.5628 | +0.5524 | 1.1628 | +1.3110 | +0.4751 | 0.4857 | 0.0520 | 0.0040 |
| Confirmatory | replogle_k562_essential | `Hybrid_Raw` | 1761 | +0.0193 | [+0.0159, +0.0228] | +0.0571 | 0.0033 | [0.0016, 0.0057] | +0.0554 | +0.1409 | +0.1383 | 0.1399 | +0.0197 | +0.9891 | 0.7577 | 0.4567 | 0.0125 |
| Confirmatory | replogle_k562_essential | `Hybrid_ProjShared` | 1761 | +0.0047 | [+0.0032, +0.0062] | +0.0752 | 0.0057 | [0.0034, 0.0087] | +0.0731 | +0.0513 | +0.0504 | 0.0926 | +0.0047 | +0.5439 | 0.7500 | 0.4584 | 0.0128 |
| Confirmatory | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 1761 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.7565 | 0.4557 | 0.0127 |
| Confirmatory | replogle_k562_essential | `NegCtrl_HybridProj_ShuffledGenes` | 1761 | -0.0003 | [-0.0020, +0.0015] | +0.0090 | 0.0001 | [0.0000, 0.0007] | +0.0087 | -0.0022 | -0.0021 | 0.1183 | -0.0003 | -0.0180 | 0.7529 | 0.4518 | 0.0130 |
| Confirmatory | replogle_k562_essential | `NegCtrl_HybridProj_DiagCov` | 1761 | -0.0005 | [-0.0023, +0.0012] | +0.0066 | 0.0000 | [0.0000, 0.0006] | +0.0064 | -0.0046 | -0.0045 | 0.1116 | -0.0005 | -0.0404 | 0.7532 | 0.4522 | 0.0129 |
| Confirmatory | replogle_k562_essential | `NegCtrl_HybridProj_SrcCorrTgtVar` | 1761 | +0.0007 | [+0.0002, +0.0012] | +0.0176 | 0.0003 | [0.0000, 0.0010] | +0.0171 | +0.0201 | +0.0197 | 0.0359 | +0.0007 | +0.5487 | 0.7554 | 0.4556 | 0.0127 |
| Confirmatory | replogle_k562_essential | `Transfer_Mean_VarScaled` | 1761 | +0.0621 | [+0.0541, +0.0702] | +0.1114 | 0.0124 | [0.0081, 0.0177] | +0.1082 | +0.2054 | +0.2016 | 0.3083 | +0.0683 | +0.6539 | 0.7637 | 0.4635 | 0.0108 |
| Confirmatory | replogle_k562_essential | `Transfer_Mean_GlobalScale` | 1761 | +0.0283 | [+0.0267, +0.0299] | +nan | nan | [nan, nan] | +nan | +0.5620 | +0.5517 | 0.0513 | +0.0284 | +10.7521 | 0.7565 | 0.4557 | 0.0117 |
| Confirmatory | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 1761 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.7565 | 0.4557 | 0.0127 |
| Confirmatory | replogle_k562_essential | `NegCtrl_TMVS_ShuffledVar` | 1761 | +0.1664 | [+0.1518, +0.1798] | +0.0792 | 0.0063 | [0.0045, 0.0085] | +0.0769 | +0.3381 | +0.3319 | 0.5013 | +0.2122 | +0.6621 | 0.7692 | 0.4495 | 0.0098 |
| Confirmatory | replogle_k562_essential | `Hybrid_on_TMVS` | 1761 | +0.0668 | [+0.0584, +0.0751] | +0.1281 | 0.0164 | [0.0112, 0.0227] | +0.1245 | +0.2109 | +0.2070 | 0.3227 | +0.0742 | +0.6413 | 0.7566 | 0.4660 | 0.0109 |
| Confirmatory | replogle_k562_essential | `FCT` | 1761 | +0.0382 | [+0.0317, +0.0443] | +0.1968 | 0.0387 | [0.0337, 0.0439] | +0.1911 | +0.1231 | +0.1208 | 0.3163 | +0.0424 | +0.3820 | 0.7507 | 0.4702 | 0.0132 |
| Confirmatory | replogle_rpe1 | `Source_Grand_Mean` | 1761 | -0.0456 | [-0.0538, -0.0374] | +0.4969 | 0.2469 | [0.2338, 0.2609] | +0.4876 | -0.1384 | -0.1359 | 0.3356 | -0.0513 | -0.4050 | 0.5002 | 0.4963 | 0.0116 |
| Confirmatory | replogle_rpe1 | `Transfer_Mean` | 1761 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.8434 | 0.4321 | 0.0127 |
| Confirmatory | replogle_rpe1 | `CIPHER_Linear` | 1761 | -0.0818 | [-0.0906, -0.0729] | +0.2189 | 0.0479 | [0.0330, 0.0656] | +0.2148 | -0.2264 | -0.2222 | 0.3682 | -0.0938 | -0.6036 | 0.5150 | 0.2506 | 0.0124 |
| Confirmatory | replogle_rpe1 | `Hybrid_Raw` | 1761 | +0.0096 | [+0.0070, +0.0121] | +0.1906 | 0.0363 | [0.0223, 0.0533] | +0.1871 | +0.1458 | +0.1432 | 0.0673 | +0.0097 | +2.1271 | 0.8267 | 0.4610 | 0.0125 |
| Confirmatory | replogle_rpe1 | `Hybrid_ProjShared` | 1761 | +0.0083 | [+0.0063, +0.0103] | +0.1798 | 0.0323 | [0.0204, 0.0467] | +0.1764 | +0.1466 | +0.1440 | 0.0577 | +0.0083 | +2.4929 | 0.8293 | 0.4585 | 0.0125 |
| Confirmatory | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 1761 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.8434 | 0.4321 | 0.0127 |
| Confirmatory | replogle_rpe1 | `NegCtrl_HybridProj_ShuffledGenes` | 1761 | -0.0011 | [-0.0015, -0.0007] | -0.0675 | 0.0046 | [0.0029, 0.0064] | -0.0662 | -0.0398 | -0.0391 | 0.0291 | -0.0011 | -1.3433 | 0.8439 | 0.4239 | 0.0128 |
| Confirmatory | replogle_rpe1 | `NegCtrl_HybridProj_DiagCov` | 1761 | -0.0011 | [-0.0015, -0.0008] | -0.0791 | 0.0063 | [0.0041, 0.0087] | -0.0777 | -0.0462 | -0.0454 | 0.0253 | -0.0011 | -1.7902 | 0.8435 | 0.4249 | 0.0127 |
| Confirmatory | replogle_rpe1 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 1761 | +0.0015 | [+0.0013, +0.0017] | +0.1385 | 0.0192 | [0.0143, 0.0249] | +0.1359 | +0.1337 | +0.1312 | 0.0113 | +0.0015 | +11.6556 | 0.8409 | 0.4367 | 0.0126 |
| Confirmatory | replogle_rpe1 | `Transfer_Mean_VarScaled` | 1761 | +0.0222 | [+0.0204, +0.0240] | +0.0881 | 0.0078 | [0.0059, 0.0099] | +0.0864 | +0.2033 | +0.1996 | 0.1112 | +0.0225 | +1.7955 | 0.8386 | 0.4280 | 0.0126 |
| Confirmatory | replogle_rpe1 | `Transfer_Mean_GlobalScale` | 1761 | +0.0044 | [+0.0040, +0.0048] | +nan | nan | [nan, nan] | +nan | +0.2683 | +0.2635 | 0.0167 | +0.0044 | +15.8083 | 0.8434 | 0.4321 | 0.0129 |
| Confirmatory | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 1761 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.8434 | 0.4321 | 0.0127 |
| Confirmatory | replogle_rpe1 | `NegCtrl_TMVS_ShuffledVar` | 1761 | -0.0075 | [-0.0091, -0.0060] | +0.0659 | 0.0043 | [0.0032, 0.0057] | +0.0646 | -0.0534 | -0.0525 | 0.1437 | -0.0077 | -0.3651 | 0.8494 | 0.4228 | 0.0130 |
| Confirmatory | replogle_rpe1 | `Hybrid_on_TMVS` | 1761 | +0.0305 | [+0.0281, +0.0330] | +0.1660 | 0.0275 | [0.0216, 0.0339] | +0.1629 | +0.2467 | +0.2422 | 0.1260 | +0.0310 | +1.9223 | 0.8270 | 0.4503 | 0.0125 |
| Confirmatory | replogle_rpe1 | `FCT` | 1761 | +0.0078 | [+0.0064, +0.0092] | +0.0311 | 0.0010 | [0.0003, 0.0019] | +0.0305 | +0.0684 | +0.0672 | 0.1158 | +0.0079 | +0.5804 | 0.8443 | 0.4189 | 0.0130 |

## R18: TMVS Shuffled-Variance 100-Permutation Null vs H-B3 (97.5% CI Threshold)

| Target Line | TMVS r_orth_cal | TMVS CI lo | Perm100 p2.5 | Perm100 median | Perm100 p97.5 | Perm100 max | Max Mismatch CI hi | H-B3 Pass |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| replogle_k562_essential | +0.1114 | +0.0899 | +0.0461 | +0.0851 | +0.1159 | +0.1211 | +0.0000 | False |
| replogle_rpe1 | +0.0881 | +0.0769 | +0.0036 | +0.0392 | +0.0578 | +0.0616 | +0.0000 | True |

## R19: Descriptive 3-Layer Residual Variance Decomposition (`R = Delta_t - Delta_src`, PREREG_A v4 §5.4)

| Target Line | Mask | N | Elig Genes | beta_ols_int [95% CI] | beta_mult | L1 Line Gain [95% CI] | L2 Per-Gene Gain [95% CI] | L3 Remaining [95% CI] |
|---|---|---:|---:|---|---:|---|---|---|
| replogle_k562_essential | `resp_top50` | 1761 | 2466 (C_dist>=5) | 0.5272 [0.5101, 0.5442] | 0.5422 | 0.3043 [0.2816, 0.3281] | 0.1112 [0.1043, 0.1200] | 0.5845 [0.5621, 0.6047] |
| replogle_k562_essential | `all_nontarget` | 1761 | 7225 (C_dist>=5) | 0.1553 [0.1481, 0.1623] | 0.1560 | 0.7099 [0.6996, 0.7202] | 0.0157 [0.0144, 0.0170] | 0.2744 [0.2649, 0.2841] |
| replogle_rpe1 | `resp_top50` | 1761 | 1630 (C_dist>=5) | 1.7272 [1.6678, 1.7928] | 2.0007 | 0.0694 [0.0605, 0.0792] | 0.0856 [0.0772, 0.0928] | 0.8449 [0.8329, 0.8580] |
| replogle_rpe1 | `all_nontarget` | 1761 | 7225 (C_dist>=5) | 0.4920 [0.4719, 0.5110] | 0.5030 | 0.0811 [0.0742, 0.0884] | 0.0770 [0.0727, 0.0815] | 0.8419 [0.8358, 0.8478] |

## S6: Split-Half Noise Ceiling (`resp_top50` vs Unselected `all_nontarget`, Amplitude-Orthogonalized & Raw)

| Target Line | N half-valid | rho_AB_orth (top50) | Ceil r_orth (top50) | rho_AB (top50) | Rep Slope (top50) | rho_AB_orth (all) | rho_AB (all) | TMVS r_orth_cal (all) | Hybrid r_orth_cal (all) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| replogle_k562_essential | 1761 | 0.8935 | 0.9715 | 0.9294 | 0.9294 | 0.3957 | 0.7870 | +0.0970 | +0.0420 |
| replogle_rpe1 | 1761 | 0.9287 | 0.9813 | 0.9304 | 0.9304 | 0.5975 | 0.6349 | +0.1215 | +0.1394 |

## S5: Expressed Genes Sensitivity (Drop Bottom 25% Target Control Variance Genes)

| Target Line | Genes Kept | rho_AB_orth (S5) | Method | Slope OLS (S5) | r_orth_cal (S5) | r_orth_cal^2 (S5) | Joint 97.5% CI | r_cal raw (S5) |
|---|---:|---:|---|---:|---:|---:|---|---:|
| replogle_k562_essential | 5419/7225 | 0.8925 | `Hybrid_ProjShared` | +0.0047 | +0.0753 | 0.0057 | [+0.0583, +0.0927] | +0.0511 |
| replogle_k562_essential | 5419/7225 | 0.8925 | `Transfer_Mean_VarScaled` | +0.0658 | +0.1119 | 0.0125 | [+0.0902, +0.1334] | +0.2160 |
| replogle_k562_essential | 5419/7225 | 0.8925 | `Transfer_Mean_GlobalScale` | +0.0282 | +nan | nan | [nan, nan] | +0.5651 |
| replogle_k562_essential | 5419/7225 | 0.8925 | `Hybrid_on_TMVS` | +0.0705 | +0.1285 | 0.0165 | [+0.1060, +0.1513] | +0.2209 |
| replogle_k562_essential | 5419/7225 | 0.8925 | `FCT` | +0.0396 | +0.1984 | 0.0394 | [+0.1854, +0.2110] | +0.1269 |
| replogle_k562_essential | 5419/7225 | 0.8925 | `NegCtrl_HybridProj_SrcCorrTgtVar` | +0.0006 | +0.0158 | 0.0003 | [+0.0032, +0.0293] | +0.0178 |
| replogle_k562_essential | 5419/7225 | 0.8925 | `NegCtrl_TMVS_ShuffledVar` | +0.1877 | +0.0860 | 0.0074 | [+0.0749, +0.0970] | +0.4115 |
| replogle_rpe1 | 5419/7225 | 0.9280 | `Hybrid_ProjShared` | +0.0085 | +0.1835 | 0.0337 | [+0.1463, +0.2196] | +0.1483 |
| replogle_rpe1 | 5419/7225 | 0.9280 | `Transfer_Mean_VarScaled` | +0.0234 | +0.0874 | 0.0076 | [+0.0774, +0.0980] | +0.2177 |
| replogle_rpe1 | 5419/7225 | 0.9280 | `Transfer_Mean_GlobalScale` | +0.0043 | +nan | nan | [nan, nan] | +0.2718 |
| replogle_rpe1 | 5419/7225 | 0.9280 | `Hybrid_on_TMVS` | +0.0319 | +0.1700 | 0.0289 | [+0.1503, +0.1886] | +0.2604 |
| replogle_rpe1 | 5419/7225 | 0.9280 | `FCT` | +0.0090 | +0.0336 | 0.0011 | [+0.0206, +0.0455] | +0.0791 |
| replogle_rpe1 | 5419/7225 | 0.9280 | `NegCtrl_HybridProj_SrcCorrTgtVar` | +0.0015 | +0.1440 | 0.0207 | [+0.1253, +0.1629] | +0.1376 |
| replogle_rpe1 | 5419/7225 | 0.9280 | `NegCtrl_TMVS_ShuffledVar` | -0.0093 | +0.0747 | 0.0056 | [+0.0663, +0.0834] | -0.0745 |

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
| replogle_k562_essential | `Transfer_Mean_GlobalScale` | 0.7565 | 0.7664 | 0.8242 | 0.8261 | 0.3043 | 0.7614 | +0.0370 |
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
| replogle_rpe1 | `Transfer_Mean_GlobalScale` | 0.8434 | 0.8556 | 0.8496 | 0.8565 | 0.3441 | 0.8470 | -0.0073 |
| replogle_rpe1 | `NegCtrl_TMVS_ShuffledVar` | 0.8494 | 0.8595 | 0.8391 | 0.8442 | 0.3111 | 0.8517 | +0.0096 |
| replogle_rpe1 | `Hybrid_on_TMVS` | 0.8270 | 0.8373 | 0.8615 | 0.8704 | 0.3890 | 0.8315 | +0.0167 |
| replogle_rpe1 | `FCT` | 0.8443 | 0.8531 | 0.8499 | 0.8568 | 0.3385 | 0.8471 | +0.0178 |
