# Line-Specific Zero-Shot Perturbation Evaluation Report (`holdout_zscore_sens_s4b`)

- **Code Freeze Commit (C1)**: `1ca9f60f49adf1ac7cce8266e40533273411db30` | **Bootstrap Seed (C2)**: `20261121`
- **Confirmatory Lines (3)**: nadig_hepg2, nadig_jurkat, vcc2025_h1
- **Secondary Lines (3)**: nourreddine_kolf21j, replogle_k562_essential, replogle_rpe1
- **Dual Confirmatory Summary (R13, 97.5% Bonferroni CIs)**: **Group A (Hybrid_ProjShared): 未确认 | Group B (Transfer_Mean_VarScaled): 未确认 (Bonferroni 97.5% CIs)**
- **Group A Primary (`Hybrid_ProjShared`)**: **未确认**
- **Group B Primary (`Transfer_Mean_VarScaled`)**: **未确认**
- **R17 Primary Correlation Metric**: `ls_resid_corr_orth_cal_top50` (amplitude-orthogonalized; raw `r_cal` reported as secondary)
- **R18 Global-Scale Scalar Baseline**: 单个标量 Transfer_Mean_GlobalScale 即可得到原始残差相关 [-0.156, +0.280]
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
| **H-A4** (Strict AND non-inferiority across 6 conventions) | A (`Hybrid_ProjShared`) | 2/3 inf (ctrl_disc_raw=2/3, ctrl_de50_r=1/3, ctrl_disc_raw_panel=2/3, cent_disc_raw=0/3, cent_de50_r=2/3, cent_disc_raw_panel=0/3; union=3) | <= 0 inf | **未确认** |
| **H-B1** (Slope 97.5% CI lo > 0.13) | B (`Transfer_Mean_VarScaled`) | 0/3 | >= 3 | **未确认** |
| **H-B2** (Joint-boot `r_orth_cal` 97.5% CI lo > 0, valid>=95%, rho_AB_orth>0) | B (`Transfer_Mean_VarScaled`) | 0/3 | >= 3 | **未确认** |
| **H-B3** (Matched `r_orth_cal` lo > max(ALL mismatched vars CI hi, 100-perm p97.5)) | B (`Transfer_Mean_VarScaled`) | 0/3 | >= 3 | **未确认** |
| **H-B4** (Strict AND non-inferiority across 6 conventions) | B (`Transfer_Mean_VarScaled`) | 1/3 inf (ctrl_disc_raw=0/3, ctrl_de50_r=0/3, ctrl_disc_raw_panel=0/3, cent_disc_raw=0/3, cent_de50_r=1/3, cent_disc_raw_panel=0/3; union=1) | <= 0 inf | **未确认** |

## Main Table: Primary View (A-T2 Own-Target-Zeroed Control Reference, 97.5% CIs, R5/R6/R11/R13/R17/R18)

| Cohort | Target Line | Method | N | Slope OLS (Top50) | 97.5% CI | r_orth_cal | r_orth_cal^2 | r_orth_cal^2 Joint 97.5% CI | Resid r_orth | r_cal (raw) | Resid r (raw) | sd_pred/sd_obs | Deming (delta=1) | Obs~Pred | Ctrl disc_raw | Ctrl de50_r | MSE |
|---|---|---|---:|---:|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Confirmatory | nadig_hepg2 | `Source_Grand_Mean` | 799 | +0.1225 | [+0.0275, +0.2494] | +0.2672 | 0.0714 | [0.0565, 0.0881] | +0.2392 | +0.1710 | +0.1549 | 0.7906 | +0.2977 | +0.1960 | 0.5002 | 0.5125 | 0.0646 |
| Confirmatory | nadig_hepg2 | `Transfer_Mean` | 799 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.9050 | 0.6444 | 0.0656 |
| Confirmatory | nadig_hepg2 | `CIPHER_Linear` | 799 | +0.0763 | [-0.0258, +0.2105] | +0.1024 | 0.0105 | [0.0046, 0.0189] | +0.0916 | +0.0966 | +0.0875 | 0.8722 | +0.2917 | +0.1003 | 0.4791 | 0.1487 | 0.0740 |
| Confirmatory | nadig_hepg2 | `Hybrid_Raw` | 799 | +0.0066 | [+0.0047, +0.0087] | +0.0984 | 0.0097 | [0.0045, 0.0171] | +0.0881 | +0.1256 | +0.1138 | 0.0582 | +0.0066 | +1.9551 | 0.8954 | 0.6492 | 0.0662 |
| Confirmatory | nadig_hepg2 | `Hybrid_ProjShared` | 799 | +0.0013 | [+0.0000, +0.0026] | +0.0485 | 0.0024 | [0.0005, 0.0057] | +0.0435 | +0.0271 | +0.0245 | 0.0537 | +0.0013 | +0.4567 | 0.8995 | 0.6449 | 0.0664 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_MismatchedCov(nadig_jurkat)` | 799 | -0.0003 | [-0.0011, +0.0005] | -0.0288 | 0.0008 | [0.0000, 0.0027] | -0.0257 | -0.0121 | -0.0110 | 0.0284 | -0.0003 | -0.3859 | 0.9046 | 0.6429 | 0.0660 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_MismatchedCov(nourreddine_kolf21j)` | 799 | +0.0001 | [-0.0016, +0.0016] | +0.0174 | 0.0003 | [0.0000, 0.0012] | +0.0156 | +0.0042 | +0.0038 | 0.0392 | +0.0001 | +0.0960 | 0.9040 | 0.6429 | 0.0661 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 799 | -0.0003 | [-0.0010, +0.0003] | -0.0164 | 0.0003 | [0.0000, 0.0011] | -0.0146 | -0.0125 | -0.0113 | 0.0292 | -0.0003 | -0.3876 | 0.9046 | 0.6425 | 0.0660 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 799 | -0.0008 | [-0.0014, -0.0002] | -0.0180 | 0.0003 | [0.0000, 0.0015] | -0.0162 | -0.0312 | -0.0283 | 0.0286 | -0.0008 | -0.9871 | 0.9047 | 0.6426 | 0.0660 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_MismatchedCov(vcc2025_h1)` | 799 | +0.0004 | [-0.0009, +0.0017] | +0.0195 | 0.0004 | [0.0000, 0.0013] | +0.0174 | +0.0057 | +0.0052 | 0.0708 | +0.0004 | +0.0728 | 0.9042 | 0.6403 | 0.0677 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_ShuffledGenes` | 799 | +0.0006 | [-0.0018, +0.0025] | +0.0178 | 0.0003 | [0.0000, 0.0011] | +0.0159 | +0.0119 | +0.0108 | 0.0565 | +0.0006 | +0.1917 | 0.9039 | 0.6401 | 0.0672 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_DiagCov` | 799 | +0.0005 | [-0.0019, +0.0023] | +0.0265 | 0.0007 | [0.0000, 0.0022] | +0.0237 | +0.0148 | +0.0134 | 0.0369 | +0.0005 | +0.3633 | 0.9037 | 0.6432 | 0.0659 |
| Confirmatory | nadig_hepg2 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 799 | -0.0002 | [-0.0008, +0.0004] | +0.0023 | 0.0000 | [0.0000, 0.0004] | +0.0021 | -0.0099 | -0.0090 | 0.0211 | -0.0002 | -0.4251 | 0.9049 | 0.6442 | 0.0658 |
| Confirmatory | nadig_hepg2 | `Transfer_Mean_VarScaled` | 799 | +0.0002 | [-0.0001, +0.0007] | +0.0108 | 0.0001 | [0.0000, 0.0011] | +0.0097 | +0.0782 | +0.0709 | 0.0035 | +0.0002 | +20.3016 | 0.9050 | 0.6444 | 0.0655 |
| Confirmatory | nadig_hepg2 | `Transfer_Mean_GlobalScale` | 799 | +0.0002 | [-0.0001, +0.0007] | +nan | nan | [nan, nan] | +nan | +0.0817 | +0.0740 | 0.0030 | +0.0002 | +24.7864 | 0.9050 | 0.6444 | 0.0655 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_MismatchedVar(nadig_jurkat)` | 799 | +0.0000 | [-0.0000, +0.0001] | +0.0039 | 0.0000 | [0.0000, 0.0005] | +0.0035 | +0.0379 | +0.0343 | 0.0011 | +0.0000 | +32.6732 | 0.9050 | 0.6444 | 0.0656 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_MismatchedVar(nourreddine_kolf21j)` | 799 | -0.0001 | [-0.0003, +0.0001] | +0.0039 | 0.0000 | [0.0000, 0.0005] | +0.0035 | -0.0671 | -0.0608 | 0.0016 | -0.0001 | -38.0871 | 0.9050 | 0.6444 | 0.0657 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 799 | +0.0000 | [-0.0000, +0.0001] | +0.0039 | 0.0000 | [0.0000, 0.0005] | +0.0035 | +0.0352 | +0.0319 | 0.0010 | +0.0000 | +30.9234 | 0.9050 | 0.6444 | 0.0656 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 799 | +0.0000 | [-0.0000, +0.0001] | +0.0039 | 0.0000 | [0.0000, 0.0005] | +0.0035 | +0.0450 | +0.0408 | 0.0011 | +0.0000 | +36.5167 | 0.9050 | 0.6444 | 0.0656 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_MismatchedVar(vcc2025_h1)` | 799 | -0.0000 | [-0.0001, +0.0001] | -0.0039 | 0.0000 | [0.0000, 0.0005] | -0.0035 | -0.0050 | -0.0045 | 0.0037 | -0.0000 | -1.2145 | 0.9050 | 0.6443 | 0.0656 |
| Confirmatory | nadig_hepg2 | `NegCtrl_TMVS_ShuffledVar` | 799 | +0.0003 | [-0.0001, +0.0009] | +0.0387 | 0.0015 | [0.0005, 0.0031] | +0.0346 | +0.0985 | +0.0893 | 0.0036 | +0.0003 | +24.5079 | 0.9050 | 0.6444 | 0.0655 |
| Confirmatory | nadig_hepg2 | `Hybrid_on_TMVS` | 799 | +0.0016 | [+0.0001, +0.0030] | +0.0489 | 0.0024 | [0.0005, 0.0057] | +0.0438 | +0.0320 | +0.0290 | 0.0539 | +0.0016 | +0.5373 | 0.8994 | 0.6449 | 0.0663 |
| Confirmatory | nadig_jurkat | `Source_Grand_Mean` | 934 | +0.2172 | [+0.1462, +0.2920] | +0.2337 | 0.0546 | [0.0440, 0.0665] | +0.2107 | +0.2502 | +0.2296 | 0.9461 | +0.7875 | +0.2427 | 0.5001 | 0.4774 | 0.0415 |
| Confirmatory | nadig_jurkat | `Transfer_Mean` | 934 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.8998 | 0.6786 | 0.0507 |
| Confirmatory | nadig_jurkat | `CIPHER_Linear` | 934 | +0.1715 | [+0.0949, +0.2528] | +0.0538 | 0.0029 | [0.0005, 0.0073] | +0.0485 | +0.1756 | +0.1612 | 1.0642 | +1.4587 | +0.1514 | 0.4615 | 0.0374 | 0.0413 |
| Confirmatory | nadig_jurkat | `Hybrid_Raw` | 934 | +0.0017 | [-0.0002, +0.0038] | +0.0403 | 0.0016 | [0.0001, 0.0052] | +0.0363 | +0.0278 | +0.0255 | 0.0683 | +0.0017 | +0.3739 | 0.9032 | 0.6776 | 0.0504 |
| Confirmatory | nadig_jurkat | `Hybrid_ProjShared` | 934 | +0.0016 | [-0.0001, +0.0032] | +0.0307 | 0.0009 | [0.0000, 0.0034] | +0.0276 | +0.0296 | +0.0271 | 0.0581 | +0.0016 | +0.4670 | 0.8989 | 0.6773 | 0.0512 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_MismatchedCov(nadig_hepg2)` | 934 | +0.0016 | [+0.0005, +0.0028] | +0.0404 | 0.0016 | [0.0004, 0.0039] | +0.0364 | +0.0330 | +0.0303 | 0.0543 | +0.0016 | +0.5585 | 0.8979 | 0.6784 | 0.0513 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_MismatchedCov(nourreddine_kolf21j)` | 934 | -0.0013 | [-0.0027, +0.0001] | -0.0080 | 0.0001 | [0.0000, 0.0009] | -0.0072 | -0.0258 | -0.0237 | 0.0568 | -0.0013 | -0.4169 | 0.9031 | 0.6764 | 0.0514 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 934 | -0.0001 | [-0.0011, +0.0009] | -0.0115 | 0.0001 | [0.0000, 0.0010] | -0.0104 | -0.0025 | -0.0023 | 0.0470 | -0.0001 | -0.0487 | 0.8996 | 0.6766 | 0.0512 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 934 | -0.0016 | [-0.0025, -0.0004] | -0.0547 | 0.0030 | [0.0011, 0.0058] | -0.0493 | -0.0439 | -0.0403 | 0.0390 | -0.0016 | -1.0322 | 0.8995 | 0.6752 | 0.0511 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_MismatchedCov(vcc2025_h1)` | 934 | +0.0037 | [+0.0018, +0.0056] | +0.0377 | 0.0014 | [0.0005, 0.0029] | +0.0340 | +0.0385 | +0.0353 | 0.1040 | +0.0037 | +0.3395 | 0.8994 | 0.6746 | 0.0527 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_ShuffledGenes` | 934 | -0.0021 | [-0.0041, -0.0003] | -0.0216 | 0.0005 | [0.0000, 0.0016] | -0.0194 | -0.0368 | -0.0338 | 0.0632 | -0.0021 | -0.5349 | 0.9012 | 0.6752 | 0.0518 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_DiagCov` | 934 | -0.0022 | [-0.0040, -0.0005] | -0.0303 | 0.0009 | [0.0001, 0.0027] | -0.0273 | -0.0476 | -0.0437 | 0.0501 | -0.0022 | -0.8721 | 0.9014 | 0.6765 | 0.0512 |
| Confirmatory | nadig_jurkat | `NegCtrl_HybridProj_SrcCorrTgtVar` | 934 | +0.0005 | [-0.0001, +0.0011] | +0.0187 | 0.0004 | [0.0000, 0.0013] | +0.0169 | +0.0154 | +0.0142 | 0.0324 | +0.0005 | +0.4362 | 0.9005 | 0.6786 | 0.0509 |
| Confirmatory | nadig_jurkat | `Transfer_Mean_VarScaled` | 934 | +0.0000 | [-0.0000, +0.0001] | +0.0118 | 0.0001 | [0.0000, 0.0008] | +0.0107 | +0.0157 | +0.0144 | 0.0011 | +0.0000 | +12.6211 | 0.8998 | 0.6786 | 0.0507 |
| Confirmatory | nadig_jurkat | `Transfer_Mean_GlobalScale` | 934 | -0.0000 | [-0.0000, -0.0000] | +nan | nan | [nan, nan] | +nan | -0.1703 | -0.1563 | 0.0001 | -0.0000 | -1408.5332 | 0.8998 | 0.6786 | 0.0507 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_MismatchedVar(nadig_hepg2)` | 934 | +0.0006 | [+0.0003, +0.0008] | +0.0309 | 0.0010 | [0.0002, 0.0021] | +0.0279 | +0.1748 | +0.1605 | 0.0035 | +0.0006 | +45.8952 | 0.8998 | 0.6786 | 0.0506 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_MismatchedVar(nourreddine_kolf21j)` | 934 | -0.0003 | [-0.0005, -0.0002] | +0.0118 | 0.0001 | [0.0000, 0.0008] | +0.0107 | -0.1418 | -0.1301 | 0.0025 | -0.0003 | -51.9897 | 0.8998 | 0.6786 | 0.0508 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 934 | +0.0000 | [-0.0000, +0.0001] | +0.0118 | 0.0001 | [0.0000, 0.0008] | +0.0107 | +0.0079 | +0.0072 | 0.0011 | +0.0000 | +6.2878 | 0.8998 | 0.6786 | 0.0507 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 934 | +0.0000 | [+0.0000, +0.0001] | +0.0118 | 0.0001 | [0.0000, 0.0008] | +0.0107 | +0.0389 | +0.0357 | 0.0012 | +0.0000 | +31.0310 | 0.8998 | 0.6786 | 0.0507 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_MismatchedVar(vcc2025_h1)` | 934 | -0.0003 | [-0.0005, -0.0001] | -0.0196 | 0.0004 | [0.0000, 0.0013] | -0.0176 | -0.0678 | -0.0623 | 0.0045 | -0.0003 | -13.7877 | 0.8998 | 0.6785 | 0.0507 |
| Confirmatory | nadig_jurkat | `NegCtrl_TMVS_ShuffledVar` | 934 | +0.0000 | [-0.0000, +0.0001] | +0.0118 | 0.0001 | [0.0000, 0.0008] | +0.0107 | +0.0157 | +0.0144 | 0.0011 | +0.0000 | +12.6213 | 0.8998 | 0.6786 | 0.0507 |
| Confirmatory | nadig_jurkat | `Hybrid_on_TMVS` | 934 | +0.0016 | [-0.0001, +0.0032] | +0.0309 | 0.0010 | [0.0000, 0.0035] | +0.0278 | +0.0298 | +0.0274 | 0.0582 | +0.0016 | +0.4709 | 0.8989 | 0.6773 | 0.0512 |
| Secondary | nourreddine_kolf21j | `Source_Grand_Mean` | 18 | +0.4954 | [+0.1030, +0.7975] | +0.2229 | 0.0497 | [0.0108, 0.1228] | +0.2070 | +0.5627 | +0.5345 | 0.9270 | +0.8680 | +0.5766 | 0.5000 | 0.3846 | 0.0223 |
| Secondary | nourreddine_kolf21j | `Transfer_Mean` | 18 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.7778 | 0.4579 | 0.0331 |
| Secondary | nourreddine_kolf21j | `CIPHER_Linear` | 18 | +0.4852 | [+0.0433, +0.8318] | -0.0350 | 0.0012 | [0.0000, 0.0615] | -0.0325 | +0.5116 | +0.4860 | 0.9984 | +0.9967 | +0.4868 | 0.5784 | -0.0432 | 0.0177 |
| Secondary | nourreddine_kolf21j | `Hybrid_Raw` | 18 | -0.0015 | [-0.0135, +0.0132] | +0.0302 | 0.0009 | [0.0000, 0.0295] | +0.0281 | -0.0161 | -0.0153 | 0.1008 | -0.0016 | -0.1517 | 0.7745 | 0.4702 | 0.0340 |
| Secondary | nourreddine_kolf21j | `Hybrid_ProjShared` | 18 | +0.0003 | [-0.0105, +0.0140] | +0.0484 | 0.0023 | [0.0000, 0.0290] | +0.0449 | +0.0030 | +0.0028 | 0.0914 | +0.0003 | +0.0310 | 0.7745 | 0.4714 | 0.0341 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_MismatchedCov(nadig_hepg2)` | 18 | +0.0054 | [-0.0049, +0.0192] | +0.0902 | 0.0081 | [0.0000, 0.0382] | +0.0837 | +0.0578 | +0.0549 | 0.0991 | +0.0055 | +0.5538 | 0.7745 | 0.4812 | 0.0342 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_MismatchedCov(nadig_jurkat)` | 18 | +0.0034 | [-0.0061, +0.0158] | +0.0770 | 0.0059 | [0.0000, 0.0310] | +0.0715 | +0.0382 | +0.0363 | 0.0926 | +0.0034 | +0.3920 | 0.7549 | 0.4768 | 0.0342 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 18 | +0.0025 | [-0.0078, +0.0164] | +0.0442 | 0.0020 | [0.0000, 0.0267] | +0.0411 | +0.0277 | +0.0263 | 0.0949 | +0.0025 | +0.2771 | 0.7614 | 0.4732 | 0.0342 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 18 | +0.0045 | [-0.0046, +0.0150] | +0.0632 | 0.0040 | [0.0000, 0.0244] | +0.0587 | +0.0535 | +0.0508 | 0.0894 | +0.0046 | +0.5685 | 0.7255 | 0.4734 | 0.0342 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_MismatchedCov(vcc2025_h1)` | 18 | +0.0129 | [-0.0082, +0.0454] | +0.1318 | 0.0174 | [0.0000, 0.0982] | +0.1223 | +0.1189 | +0.1130 | 0.1142 | +0.0131 | +0.9892 | 0.7745 | 0.4616 | 0.0341 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_ShuffledGenes` | 18 | +0.0027 | [-0.0077, +0.0185] | +0.0528 | 0.0028 | [0.0000, 0.0290] | +0.0490 | +0.0288 | +0.0273 | 0.0997 | +0.0027 | +0.2740 | 0.7255 | 0.4709 | 0.0345 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_DiagCov` | 18 | +0.0039 | [-0.0071, +0.0182] | +0.0658 | 0.0043 | [0.0000, 0.0353] | +0.0611 | +0.0453 | +0.0430 | 0.0904 | +0.0039 | +0.4755 | 0.7288 | 0.4771 | 0.0342 |
| Secondary | nourreddine_kolf21j | `NegCtrl_HybridProj_SrcCorrTgtVar` | 18 | +0.0057 | [-0.0024, +0.0164] | +0.1157 | 0.0134 | [0.0000, 0.0544] | +0.1074 | +0.0858 | +0.0815 | 0.0704 | +0.0058 | +1.1574 | 0.7680 | 0.4808 | 0.0336 |
| Secondary | nourreddine_kolf21j | `Transfer_Mean_VarScaled` | 18 | -0.0013 | [-0.0022, -0.0001] | -0.0168 | 0.0003 | [0.0000, 0.0067] | -0.0156 | -0.4794 | -0.4553 | 0.0029 | -0.0013 | -159.7697 | 0.7778 | 0.4578 | 0.0332 |
| Secondary | nourreddine_kolf21j | `Transfer_Mean_GlobalScale` | 18 | -0.0013 | [-0.0022, -0.0001] | +nan | nan | [nan, nan] | +nan | -0.5146 | -0.4888 | 0.0026 | -0.0013 | -189.7535 | 0.7778 | 0.4579 | 0.0332 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_MismatchedVar(nadig_hepg2)` | 18 | +0.0011 | [+0.0001, +0.0018] | -0.0470 | 0.0022 | [0.0000, 0.0091] | -0.0437 | +0.4138 | +0.3931 | 0.0027 | +0.0011 | +145.1916 | 0.7778 | 0.4578 | 0.0330 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_MismatchedVar(nadig_jurkat)` | 18 | -0.0003 | [-0.0005, -0.0000] | -0.0168 | 0.0003 | [0.0000, 0.0067] | -0.0156 | -0.2663 | -0.2529 | 0.0012 | -0.0003 | -219.7586 | 0.7778 | 0.4578 | 0.0331 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 18 | -0.0003 | [-0.0005, -0.0000] | -0.0168 | 0.0003 | [0.0000, 0.0067] | -0.0156 | -0.2817 | -0.2676 | 0.0012 | -0.0003 | -227.0941 | 0.7778 | 0.4578 | 0.0331 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 18 | -0.0002 | [-0.0004, -0.0000] | -0.0168 | 0.0003 | [0.0000, 0.0067] | -0.0156 | -0.2137 | -0.2030 | 0.0011 | -0.0002 | -188.2389 | 0.7778 | 0.4578 | 0.0331 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_MismatchedVar(vcc2025_h1)` | 18 | -0.0002 | [-0.0005, +0.0001] | +0.0290 | 0.0008 | [0.0000, 0.0060] | +0.0269 | -0.0610 | -0.0580 | 0.0040 | -0.0002 | -14.5004 | 0.7778 | 0.4580 | 0.0331 |
| Secondary | nourreddine_kolf21j | `NegCtrl_TMVS_ShuffledVar` | 18 | -0.0013 | [-0.0022, -0.0001] | -0.0168 | 0.0003 | [0.0000, 0.0067] | -0.0156 | -0.4794 | -0.4553 | 0.0029 | -0.0013 | -159.7696 | 0.7778 | 0.4578 | 0.0332 |
| Secondary | nourreddine_kolf21j | `Hybrid_on_TMVS` | 18 | -0.0010 | [-0.0117, +0.0135] | +0.0482 | 0.0023 | [0.0000, 0.0290] | +0.0447 | -0.0119 | -0.0113 | 0.0915 | -0.0010 | -0.1240 | 0.7745 | 0.4714 | 0.0342 |
| Secondary | replogle_k562_essential | `Source_Grand_Mean` | 859 | +0.5389 | [+0.4050, +0.6657] | +0.0785 | 0.0062 | [0.0031, 0.0101] | +0.0747 | +0.5118 | +0.4962 | 1.0862 | +1.1806 | +0.4568 | 0.5001 | 0.4247 | 0.0361 |
| Secondary | replogle_k562_essential | `Transfer_Mean` | 859 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.9095 | 0.7048 | 0.0520 |
| Secondary | replogle_k562_essential | `CIPHER_Linear` | 859 | +0.5655 | [+0.4258, +0.6980] | -0.0117 | 0.0001 | [0.0000, 0.0018] | -0.0111 | +0.4833 | +0.4685 | 1.2071 | +1.4827 | +0.3881 | 0.4170 | -0.0021 | 0.0358 |
| Secondary | replogle_k562_essential | `Hybrid_Raw` | 859 | +0.0044 | [+0.0028, +0.0061] | +0.0109 | 0.0001 | [0.0000, 0.0012] | +0.0104 | +0.0628 | +0.0609 | 0.0725 | +0.0044 | +0.8395 | 0.9143 | 0.7039 | 0.0511 |
| Secondary | replogle_k562_essential | `Hybrid_ProjShared` | 859 | -0.0009 | [-0.0022, +0.0002] | -0.0011 | 0.0000 | [0.0000, 0.0006] | -0.0010 | -0.0179 | -0.0174 | 0.0546 | -0.0010 | -0.3188 | 0.9087 | 0.7038 | 0.0528 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(nadig_hepg2)` | 859 | +0.0008 | [+0.0000, +0.0015] | -0.0105 | 0.0001 | [0.0000, 0.0007] | -0.0100 | +0.0169 | +0.0164 | 0.0493 | +0.0008 | +0.3328 | 0.9091 | 0.7028 | 0.0526 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(nadig_jurkat)` | 859 | -0.0001 | [-0.0011, +0.0008] | -0.0250 | 0.0006 | [0.0001, 0.0017] | -0.0238 | -0.0018 | -0.0017 | 0.0429 | -0.0001 | -0.0406 | 0.9087 | 0.7025 | 0.0526 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(nourreddine_kolf21j)` | 859 | -0.0034 | [-0.0053, -0.0018] | +0.0066 | 0.0000 | [0.0000, 0.0007] | +0.0063 | -0.0621 | -0.0602 | 0.0566 | -0.0034 | -1.0625 | 0.9140 | 0.7031 | 0.0529 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 859 | -0.0009 | [-0.0017, -0.0002] | +0.0125 | 0.0002 | [0.0000, 0.0009] | +0.0119 | -0.0238 | -0.0230 | 0.0374 | -0.0009 | -0.6153 | 0.9096 | 0.7038 | 0.0524 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov(vcc2025_h1)` | 859 | +0.0014 | [-0.0005, +0.0033] | +0.0201 | 0.0004 | [0.0000, 0.0012] | +0.0192 | +0.0146 | +0.0142 | 0.0958 | +0.0014 | +0.1482 | 0.9065 | 0.6996 | 0.0544 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_ShuffledGenes` | 859 | -0.0057 | [-0.0080, -0.0036] | +0.0070 | 0.0000 | [0.0000, 0.0009] | +0.0066 | -0.0907 | -0.0880 | 0.0645 | -0.0057 | -1.3634 | 0.9148 | 0.7024 | 0.0534 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_DiagCov` | 859 | -0.0055 | [-0.0078, -0.0035] | +0.0064 | 0.0000 | [0.0000, 0.0010] | +0.0061 | -0.1042 | -0.1010 | 0.0549 | -0.0056 | -1.8398 | 0.9147 | 0.7035 | 0.0528 |
| Secondary | replogle_k562_essential | `NegCtrl_HybridProj_SrcCorrTgtVar` | 859 | -0.0004 | [-0.0011, +0.0002] | +0.0076 | 0.0001 | [0.0000, 0.0006] | +0.0072 | -0.0150 | -0.0145 | 0.0303 | -0.0004 | -0.4781 | 0.9102 | 0.7047 | 0.0523 |
| Secondary | replogle_k562_essential | `Transfer_Mean_VarScaled` | 859 | +0.0000 | [-0.0001, +0.0001] | +0.0569 | 0.0032 | [0.0011, 0.0063] | +0.0541 | +0.0118 | +0.0115 | 0.0013 | +0.0000 | +9.0827 | 0.9095 | 0.7048 | 0.0520 |
| Secondary | replogle_k562_essential | `Transfer_Mean_GlobalScale` | 859 | -0.0001 | [-0.0001, -0.0001] | +nan | nan | [nan, nan] | +nan | -0.4816 | -0.4669 | 0.0002 | -0.0001 | -2370.0588 | 0.9095 | 0.7048 | 0.0520 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(nadig_hepg2)` | 859 | +0.0016 | [+0.0012, +0.0021] | +0.0412 | 0.0017 | [0.0005, 0.0038] | +0.0392 | +0.4295 | +0.4164 | 0.0039 | +0.0016 | +106.2713 | 0.9095 | 0.7048 | 0.0518 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(nadig_jurkat)` | 859 | +0.0000 | [-0.0000, +0.0001] | +0.0569 | 0.0032 | [0.0011, 0.0063] | +0.0541 | +0.0347 | +0.0337 | 0.0013 | +0.0000 | +26.7655 | 0.9095 | 0.7048 | 0.0520 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(nourreddine_kolf21j)` | 859 | -0.0011 | [-0.0014, -0.0008] | +0.0569 | 0.0032 | [0.0011, 0.0063] | +0.0541 | -0.4067 | -0.3943 | 0.0028 | -0.0011 | -139.6838 | 0.9095 | 0.7048 | 0.0522 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 859 | +0.0001 | [+0.0001, +0.0002] | +0.0569 | 0.0032 | [0.0011, 0.0063] | +0.0541 | +0.1022 | +0.0991 | 0.0013 | +0.0001 | +78.5119 | 0.9095 | 0.7048 | 0.0520 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_MismatchedVar(vcc2025_h1)` | 859 | -0.0007 | [-0.0009, -0.0004] | -0.0572 | 0.0033 | [0.0013, 0.0062] | -0.0544 | -0.1342 | -0.1301 | 0.0052 | -0.0007 | -25.0697 | 0.9095 | 0.7047 | 0.0521 |
| Secondary | replogle_k562_essential | `NegCtrl_TMVS_ShuffledVar` | 859 | +0.0000 | [-0.0001, +0.0001] | +0.0569 | 0.0032 | [0.0011, 0.0063] | +0.0541 | +0.0118 | +0.0115 | 0.0013 | +0.0000 | +9.0836 | 0.9095 | 0.7048 | 0.0520 |
| Secondary | replogle_k562_essential | `Hybrid_on_TMVS` | 859 | -0.0009 | [-0.0022, +0.0003] | +0.0002 | 0.0000 | [0.0000, 0.0006] | +0.0002 | -0.0177 | -0.0171 | 0.0547 | -0.0009 | -0.3131 | 0.9087 | 0.7038 | 0.0528 |
| Secondary | replogle_rpe1 | `Source_Grand_Mean` | 1070 | -0.1062 | [-0.1360, -0.0780] | +0.4644 | 0.2157 | [0.1911, 0.2415] | +0.4351 | -0.2450 | -0.2330 | 0.4560 | -0.1318 | -0.5110 | 0.5003 | 0.6149 | 0.0706 |
| Secondary | replogle_rpe1 | `Transfer_Mean` | 1070 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.9418 | 0.6652 | 0.0703 |
| Secondary | replogle_rpe1 | `CIPHER_Linear` | 1070 | -0.1821 | [-0.2053, -0.1577] | +0.2790 | 0.0778 | [0.0603, 0.0980] | +0.2614 | -0.3803 | -0.3617 | 0.5035 | -0.2310 | -0.7185 | 0.5107 | 0.3370 | 0.0799 |
| Secondary | replogle_rpe1 | `Hybrid_Raw` | 1070 | +0.0057 | [+0.0048, +0.0068] | +0.2377 | 0.0565 | [0.0428, 0.0718] | +0.2227 | +0.2849 | +0.2710 | 0.0212 | +0.0057 | +12.8069 | 0.9387 | 0.6752 | 0.0695 |
| Secondary | replogle_rpe1 | `Hybrid_ProjShared` | 1070 | +0.0029 | [+0.0023, +0.0035] | +0.1610 | 0.0259 | [0.0187, 0.0347] | +0.1508 | +0.1727 | +0.1642 | 0.0176 | +0.0029 | +9.3289 | 0.9402 | 0.6713 | 0.0697 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(nadig_hepg2)` | 1070 | +0.0005 | [+0.0003, +0.0008] | +0.0467 | 0.0022 | [0.0009, 0.0040] | +0.0437 | +0.0411 | +0.0391 | 0.0140 | +0.0005 | +2.8010 | 0.9420 | 0.6670 | 0.0701 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(nadig_jurkat)` | 1070 | -0.0002 | [-0.0004, -0.0000] | -0.0160 | 0.0003 | [0.0000, 0.0010] | -0.0150 | -0.0199 | -0.0189 | 0.0103 | -0.0002 | -1.8308 | 0.9418 | 0.6646 | 0.0704 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(nourreddine_kolf21j)` | 1070 | +0.0002 | [-0.0003, +0.0007] | -0.0164 | 0.0003 | [0.0000, 0.0011] | -0.0153 | +0.0113 | +0.0107 | 0.0144 | +0.0002 | +0.7446 | 0.9421 | 0.6634 | 0.0705 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 1070 | -0.0003 | [-0.0005, -0.0001] | -0.0234 | 0.0005 | [0.0000, 0.0016] | -0.0219 | -0.0268 | -0.0255 | 0.0109 | -0.0003 | -2.3372 | 0.9416 | 0.6643 | 0.0705 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov(vcc2025_h1)` | 1070 | -0.0008 | [-0.0013, -0.0003] | -0.0323 | 0.0010 | [0.0003, 0.0023] | -0.0303 | -0.0365 | -0.0347 | 0.0237 | -0.0008 | -1.4617 | 0.9420 | 0.6611 | 0.0708 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_ShuffledGenes` | 1070 | +0.0001 | [-0.0005, +0.0008] | -0.0369 | 0.0014 | [0.0005, 0.0027] | -0.0346 | +0.0038 | +0.0036 | 0.0169 | +0.0001 | +0.2121 | 0.9419 | 0.6621 | 0.0707 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_DiagCov` | 1070 | +0.0003 | [-0.0002, +0.0010] | -0.0267 | 0.0007 | [0.0001, 0.0020] | -0.0250 | +0.0267 | +0.0254 | 0.0129 | +0.0003 | +1.9737 | 0.9418 | 0.6631 | 0.0705 |
| Secondary | replogle_rpe1 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 1070 | -0.0001 | [-0.0003, +0.0001] | -0.0186 | 0.0003 | [0.0000, 0.0012] | -0.0174 | -0.0150 | -0.0143 | 0.0083 | -0.0001 | -1.7252 | 0.9421 | 0.6647 | 0.0703 |
| Secondary | replogle_rpe1 | `Transfer_Mean_VarScaled` | 1070 | -0.0000 | [-0.0000, -0.0000] | +0.0211 | 0.0004 | [0.0000, 0.0020] | +0.0198 | -0.0310 | -0.0295 | 0.0006 | -0.0000 | -53.1147 | 0.9418 | 0.6653 | 0.0703 |
| Secondary | replogle_rpe1 | `Transfer_Mean_GlobalScale` | 1070 | -0.0000 | [-0.0000, -0.0000] | +nan | nan | [nan, nan] | +nan | -0.3940 | -0.3748 | 0.0000 | -0.0000 | -10137.4726 | 0.9418 | 0.6652 | 0.0703 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(nadig_hepg2)` | 1070 | -0.0006 | [-0.0006, -0.0005] | -0.0060 | 0.0000 | [0.0000, 0.0011] | -0.0057 | -0.3522 | -0.3350 | 0.0017 | -0.0006 | -200.3922 | 0.9418 | 0.6652 | 0.0702 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(nadig_jurkat)` | 1070 | +0.0000 | [+0.0000, +0.0000] | +0.0211 | 0.0004 | [0.0000, 0.0020] | +0.0198 | +0.0222 | +0.0211 | 0.0006 | +0.0000 | +38.0461 | 0.9418 | 0.6653 | 0.0703 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(nourreddine_kolf21j)` | 1070 | +0.0004 | [+0.0004, +0.0005] | +0.0211 | 0.0004 | [0.0000, 0.0020] | +0.0198 | +0.3523 | +0.3351 | 0.0012 | +0.0004 | +278.1028 | 0.9418 | 0.6653 | 0.0703 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 1070 | +0.0000 | [+0.0000, +0.0000] | +0.0211 | 0.0004 | [0.0000, 0.0020] | +0.0198 | +0.0400 | +0.0380 | 0.0006 | +0.0000 | +68.4294 | 0.9418 | 0.6653 | 0.0703 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_MismatchedVar(vcc2025_h1)` | 1070 | +0.0001 | [+0.0001, +0.0002] | -0.0136 | 0.0002 | [0.0000, 0.0015] | -0.0128 | +0.0574 | +0.0546 | 0.0022 | +0.0001 | +24.8952 | 0.9418 | 0.6652 | 0.0703 |
| Secondary | replogle_rpe1 | `NegCtrl_TMVS_ShuffledVar` | 1070 | -0.0000 | [-0.0000, -0.0000] | +0.0211 | 0.0004 | [0.0000, 0.0020] | +0.0198 | -0.0310 | -0.0295 | 0.0006 | -0.0000 | -53.1154 | 0.9418 | 0.6653 | 0.0703 |
| Secondary | replogle_rpe1 | `Hybrid_on_TMVS` | 1070 | +0.0029 | [+0.0023, +0.0035] | +0.1615 | 0.0261 | [0.0188, 0.0349] | +0.1513 | +0.1715 | +0.1631 | 0.0176 | +0.0029 | +9.2524 | 0.9402 | 0.6713 | 0.0697 |
| Confirmatory | vcc2025_h1 | `Source_Grand_Mean` | 80 | -0.0619 | [-0.0930, -0.0277] | +0.1447 | 0.0209 | [0.0074, 0.0422] | +0.1437 | -0.1998 | -0.1987 | 0.3117 | -0.0683 | -0.6374 | 0.5000 | 0.2579 | 0.0693 |
| Confirmatory | vcc2025_h1 | `Transfer_Mean` | 80 | +0.0000 | [+0.0000, +0.0000] | +nan | nan | [nan, nan] | +nan | +nan | +nan | 0.0000 | +0.0000 | +nan | 0.9109 | 0.5156 | 0.0697 |
| Confirmatory | vcc2025_h1 | `CIPHER_Linear` | 80 | -0.0901 | [-0.1201, -0.0555] | -0.0438 | 0.0019 | [0.0000, 0.0138] | -0.0435 | -0.2637 | -0.2622 | 0.3435 | -0.1011 | -0.7634 | 0.4546 | 0.0062 | 0.0711 |
| Confirmatory | vcc2025_h1 | `Hybrid_Raw` | 80 | -0.0020 | [-0.0080, +0.0026] | -0.0486 | 0.0024 | [0.0000, 0.0149] | -0.0483 | -0.0258 | -0.0257 | 0.0785 | -0.0020 | -0.3270 | 0.8731 | 0.4790 | 0.0729 |
| Confirmatory | vcc2025_h1 | `Hybrid_ProjShared` | 80 | -0.0019 | [-0.0075, +0.0021] | -0.0495 | 0.0024 | [0.0000, 0.0150] | -0.0492 | -0.0254 | -0.0253 | 0.0767 | -0.0019 | -0.3297 | 0.8731 | 0.4822 | 0.0724 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_MismatchedCov(nadig_hepg2)` | 80 | -0.0009 | [-0.0042, +0.0014] | -0.0470 | 0.0022 | [0.0000, 0.0113] | -0.0467 | -0.0229 | -0.0228 | 0.0416 | -0.0009 | -0.5477 | 0.9092 | 0.5054 | 0.0708 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_MismatchedCov(nadig_jurkat)` | 80 | +0.0007 | [-0.0006, +0.0021] | +0.0044 | 0.0000 | [0.0000, 0.0031] | +0.0044 | +0.0241 | +0.0239 | 0.0297 | +0.0007 | +0.8057 | 0.9047 | 0.5091 | 0.0702 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_MismatchedCov(nourreddine_kolf21j)` | 80 | -0.0014 | [-0.0042, +0.0006] | -0.0449 | 0.0020 | [0.0000, 0.0100] | -0.0447 | -0.0448 | -0.0445 | 0.0312 | -0.0014 | -1.4287 | 0.9065 | 0.5067 | 0.0703 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_MismatchedCov(replogle_k562_essential)` | 80 | +0.0001 | [-0.0014, +0.0019] | +0.0232 | 0.0005 | [0.0000, 0.0059] | +0.0230 | +0.0024 | +0.0024 | 0.0292 | +0.0001 | +0.0815 | 0.9051 | 0.5121 | 0.0702 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_MismatchedCov(replogle_rpe1)` | 80 | -0.0008 | [-0.0037, +0.0014] | -0.0048 | 0.0000 | [0.0000, 0.0079] | -0.0048 | -0.0270 | -0.0269 | 0.0298 | -0.0008 | -0.9027 | 0.9109 | 0.5135 | 0.0700 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_ShuffledGenes` | 80 | -0.0016 | [-0.0058, +0.0018] | +0.0019 | 0.0000 | [0.0000, 0.0033] | +0.0019 | -0.0281 | -0.0279 | 0.0582 | -0.0016 | -0.4798 | 0.9082 | 0.5062 | 0.0722 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_DiagCov` | 80 | -0.0009 | [-0.0044, +0.0018] | +0.0066 | 0.0000 | [0.0000, 0.0069] | +0.0066 | -0.0267 | -0.0266 | 0.0357 | -0.0010 | -0.7446 | 0.9089 | 0.5160 | 0.0700 |
| Confirmatory | vcc2025_h1 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 80 | -0.0005 | [-0.0022, +0.0007] | -0.0325 | 0.0011 | [0.0000, 0.0100] | -0.0323 | -0.0298 | -0.0296 | 0.0161 | -0.0005 | -1.8405 | 0.9074 | 0.5128 | 0.0699 |
| Confirmatory | vcc2025_h1 | `Transfer_Mean_VarScaled` | 80 | -0.0000 | [-0.0001, +0.0001] | +0.0196 | 0.0004 | [0.0000, 0.0030] | +0.0194 | -0.0016 | -0.0016 | 0.0017 | -0.0000 | -0.9425 | 0.9109 | 0.5157 | 0.0697 |
| Confirmatory | vcc2025_h1 | `Transfer_Mean_GlobalScale` | 80 | +0.0001 | [+0.0000, +0.0001] | +nan | nan | [nan, nan] | +nan | +0.2819 | +0.2803 | 0.0002 | +0.0001 | +1366.8696 | 0.9109 | 0.5156 | 0.0697 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_MismatchedVar(nadig_hepg2)` | 80 | -0.0003 | [-0.0003, -0.0002] | +0.0064 | 0.0000 | [0.0000, 0.0024] | +0.0064 | -0.2704 | -0.2689 | 0.0009 | -0.0003 | -289.1933 | 0.9109 | 0.5156 | 0.0696 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_MismatchedVar(nadig_jurkat)` | 80 | +0.0000 | [+0.0000, +0.0000] | -0.0065 | 0.0000 | [0.0000, 0.0024] | -0.0065 | +0.1975 | +0.1964 | 0.0001 | +0.0000 | +2109.5014 | 0.9109 | 0.5156 | 0.0697 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_MismatchedVar(nourreddine_kolf21j)` | 80 | +0.0002 | [+0.0001, +0.0003] | -0.0065 | 0.0000 | [0.0000, 0.0024] | -0.0065 | +0.2808 | +0.2793 | 0.0007 | +0.0002 | +383.5548 | 0.9109 | 0.5156 | 0.0697 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_MismatchedVar(replogle_k562_essential)` | 80 | +0.0000 | [+0.0000, +0.0000] | -0.0065 | 0.0000 | [0.0000, 0.0024] | -0.0065 | +0.2182 | +0.2170 | 0.0001 | +0.0000 | +2070.5117 | 0.9109 | 0.5156 | 0.0697 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_MismatchedVar(replogle_rpe1)` | 80 | +0.0000 | [+0.0000, +0.0000] | -0.0065 | 0.0000 | [0.0000, 0.0024] | -0.0065 | +0.0732 | +0.0728 | 0.0001 | +0.0000 | +1053.7848 | 0.9109 | 0.5156 | 0.0697 |
| Confirmatory | vcc2025_h1 | `NegCtrl_TMVS_ShuffledVar` | 80 | +0.0001 | [+0.0000, +0.0002] | +0.0204 | 0.0004 | [0.0000, 0.0041] | +0.0203 | +0.0849 | +0.0845 | 0.0015 | +0.0001 | +54.8305 | 0.9109 | 0.5156 | 0.0697 |
| Confirmatory | vcc2025_h1 | `Hybrid_on_TMVS` | 80 | -0.0019 | [-0.0075, +0.0021] | -0.0490 | 0.0024 | [0.0000, 0.0149] | -0.0487 | -0.0254 | -0.0253 | 0.0768 | -0.0020 | -0.3293 | 0.8733 | 0.4822 | 0.0724 |

## R18: TMVS Shuffled-Variance 100-Permutation Null vs H-B3 (97.5% CI Threshold)

| Target Line | TMVS r_orth_cal | TMVS CI lo | Perm100 p2.5 | Perm100 median | Perm100 p97.5 | Perm100 max | Max Mismatch CI hi | H-B3 Pass |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| nadig_hepg2 | +0.0108 | -0.0101 | -0.0604 | +0.0085 | +0.0441 | +0.0651 | +0.0230 | False |
| nadig_jurkat | +0.0118 | -0.0047 | +0.0118 | +0.0118 | +0.0118 | +0.0118 | +0.0463 | False |
| nourreddine_kolf21j | -0.0168 | -0.0820 | -0.0168 | -0.0168 | -0.0168 | -0.0168 | +0.0771 | False |
| replogle_k562_essential | +0.0569 | +0.0337 | +0.0569 | +0.0569 | +0.0569 | +0.0569 | +0.0795 | False |
| replogle_rpe1 | +0.0211 | -0.0058 | +0.0211 | +0.0211 | +0.0211 | +0.0211 | +0.0444 | False |
| vcc2025_h1 | +0.0196 | -0.0188 | -0.0460 | -0.0029 | +0.0516 | +0.0542 | +0.0470 | False |

## R19: Descriptive 3-Layer Residual Variance Decomposition (`R = Delta_t - Delta_src`, PREREG_A v4 §5.4)

| Target Line | Mask | N | Elig Genes | beta_ols_int [95% CI] | beta_mult | L1 Line Gain [95% CI] | L2 Per-Gene Gain [95% CI] | L3 Remaining [95% CI] |
|---|---|---:|---:|---|---:|---|---|---|
| nadig_hepg2 | `resp_top50` | 799 | 1359 (C_dist>=5) | 0.9153 [0.7971, 1.0443] | 1.0863 | 0.0055 [0.0000, 0.0363] | 0.2404 [0.2034, 0.2753] | 0.7541 [0.6967, 0.7953] |
| nadig_hepg2 | `all_nontarget` | 799 | 6148 (C_dist>=5) | 0.6136 [0.5850, 0.6437] | 0.6137 | 0.0920 [0.0765, 0.1090] | 0.0505 [0.0461, 0.0561] | 0.8574 [0.8426, 0.8698] |
| nadig_jurkat | `resp_top50` | 934 | 1810 (C_dist>=5) | 0.8531 [0.7974, 0.9090] | 0.9184 | 0.0244 [0.0092, 0.0476] | 0.2212 [0.1923, 0.2503] | 0.7543 [0.7149, 0.7923] |
| nadig_jurkat | `all_nontarget` | 934 | 6148 (C_dist>=5) | 0.3754 [0.3589, 0.3920] | 0.3727 | 0.2968 [0.2797, 0.3143] | 0.0356 [0.0331, 0.0388] | 0.6676 [0.6514, 0.6830] |
| nourreddine_kolf21j | `resp_top50` | 18 | 8 (C_dist>=5) | 0.5109 [0.4035, 0.8377] | 0.5210 | 0.2389 [0.0117, 0.4628] | 0.0052 [0.0000, 0.0047] | 0.7559 [0.5358, 0.9873] |
| nourreddine_kolf21j | `all_nontarget` | 18 | 6148 (C_dist>=5) | 0.1791 [0.1193, 0.2219] | 0.1800 | 0.4997 [0.3999, 0.5761] | 0.0108 [0.0076, 0.0353] | 0.4895 [0.4010, 0.5846] |
| replogle_k562_essential | `resp_top50` | 859 | 1450 (C_dist>=5) | 0.6115 [0.5601, 0.6759] | 0.6593 | 0.2180 [0.1427, 0.2952] | 0.2079 [0.1848, 0.2346] | 0.5741 [0.4941, 0.6550] |
| replogle_k562_essential | `all_nontarget` | 859 | 6148 (C_dist>=5) | 0.3396 [0.3253, 0.3549] | 0.3386 | 0.4399 [0.4152, 0.4633] | 0.0323 [0.0299, 0.0350] | 0.5278 [0.5052, 0.5512] |
| replogle_rpe1 | `resp_top50` | 1070 | 1096 (C_dist>=5) | 1.7377 [1.6270, 1.8593] | 1.8955 | 0.1404 [0.1082, 0.1771] | 0.1382 [0.1109, 0.1721] | 0.7213 [0.6634, 0.7668] |
| replogle_rpe1 | `all_nontarget` | 1070 | 6148 (C_dist>=5) | 0.6951 [0.6631, 0.7286] | 0.6870 | 0.0416 [0.0322, 0.0517] | 0.1030 [0.0921, 0.1177] | 0.8554 [0.8465, 0.8622] |
| vcc2025_h1 | `resp_top50` | 80 | 213 (C_dist>=5) | 1.8745 [1.4799, 2.3359] | 2.1315 | 0.0786 [0.0294, 0.1522] | 0.0449 [0.0066, 0.0687] | 0.8766 [0.8056, 0.9424] |
| vcc2025_h1 | `all_nontarget` | 80 | 6148 (C_dist>=5) | 0.4713 [0.3478, 0.6241] | 0.4713 | 0.0889 [0.0362, 0.1667] | 0.0424 [0.0300, 0.0672] | 0.8687 [0.7958, 0.9086] |

## S6: Split-Half Noise Ceiling (`resp_top50` vs Unselected `all_nontarget`, Amplitude-Orthogonalized & Raw)

| Target Line | N half-valid | rho_AB_orth (top50) | Ceil r_orth (top50) | rho_AB (top50) | Rep Slope (top50) | rho_AB_orth (all) | rho_AB (all) | TMVS r_orth_cal (all) | Hybrid r_orth_cal (all) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| nadig_hepg2 | 799 | 0.6687 | 0.8953 | 0.6964 | 0.6966 | 0.3903 | 0.4478 | -0.0108 | +0.0270 |
| nadig_jurkat | 934 | 0.6841 | 0.9013 | 0.7275 | 0.7276 | 0.3024 | 0.4651 | -0.0057 | +0.0218 |
| nourreddine_kolf21j | 18 | 0.7575 | 0.9285 | 0.8219 | 0.8222 | 0.1331 | 0.4509 | +0.0094 | +0.0161 |
| replogle_k562_essential | 859 | 0.8266 | 0.9514 | 0.8865 | 0.8865 | 0.4015 | 0.6239 | -0.0003 | +0.0021 |
| replogle_rpe1 | 1070 | 0.7819 | 0.9368 | 0.8258 | 0.8259 | 0.5306 | 0.5664 | -0.0046 | +0.1363 |
| vcc2025_h1 | 80 | 0.9741 | 0.9934 | 0.9785 | 0.9786 | 0.8734 | 0.8943 | -0.0169 | -0.0037 |

## S5: Expressed Genes Sensitivity (Drop Bottom 25% Target Control Variance Genes)

| Target Line | Genes Kept | rho_AB_orth (S5) | Method | Slope OLS (S5) | r_orth_cal (S5) | r_orth_cal^2 (S5) | Joint 97.5% CI | r_cal raw (S5) |
|---|---:|---:|---|---:|---:|---:|---|---:|
| nadig_hepg2 | 4611/6148 | 0.6995 | `Hybrid_ProjShared` | +0.0022 | +0.0677 | 0.0046 | [+0.0396, +0.0966] | +0.0420 |
| nadig_hepg2 | 4611/6148 | 0.6995 | `Transfer_Mean_VarScaled` | -0.0001 | -0.0406 | 0.0017 | [-0.0592, -0.0205] | -0.0279 |
| nadig_hepg2 | 4611/6148 | 0.6995 | `Transfer_Mean_GlobalScale` | -0.0000 | +nan | nan | [nan, nan] | -0.0023 |
| nadig_hepg2 | 4611/6148 | 0.6995 | `Hybrid_on_TMVS` | +0.0021 | +0.0666 | 0.0044 | [+0.0385, +0.0956] | +0.0404 |
| nadig_hepg2 | 4611/6148 | 0.6995 | `NegCtrl_HybridProj_SrcCorrTgtVar` | -0.0000 | +0.0026 | 0.0000 | [-0.0140, +0.0197] | -0.0008 |
| nadig_hepg2 | 4611/6148 | 0.6995 | `NegCtrl_TMVS_ShuffledVar` | +0.0000 | +0.0026 | 0.0000 | [-0.0108, +0.0163] | +0.0014 |
| nadig_jurkat | 4611/6148 | 0.7181 | `Hybrid_ProjShared` | +0.0018 | +0.0334 | 0.0011 | [+0.0063, +0.0617] | +0.0321 |
| nadig_jurkat | 4611/6148 | 0.7181 | `Transfer_Mean_VarScaled` | -0.0000 | +0.0060 | 0.0000 | [-0.0096, +0.0212] | -0.0055 |
| nadig_jurkat | 4611/6148 | 0.7181 | `Transfer_Mean_GlobalScale` | -0.0000 | +nan | nan | [nan, nan] | -0.1695 |
| nadig_jurkat | 4611/6148 | 0.7181 | `Hybrid_on_TMVS` | +0.0018 | +0.0335 | 0.0011 | [+0.0065, +0.0618] | +0.0320 |
| nadig_jurkat | 4611/6148 | 0.7181 | `NegCtrl_HybridProj_SrcCorrTgtVar` | +0.0003 | +0.0113 | 0.0001 | [-0.0067, +0.0288] | +0.0108 |
| nadig_jurkat | 4611/6148 | 0.7181 | `NegCtrl_TMVS_ShuffledVar` | -0.0000 | +0.0060 | 0.0000 | [-0.0096, +0.0212] | -0.0055 |
| nourreddine_kolf21j | 4611/6148 | 0.7546 | `Hybrid_ProjShared` | -0.0082 | +0.0073 | 0.0001 | [-0.1378, +0.1290] | -0.0890 |
| nourreddine_kolf21j | 4611/6148 | 0.7546 | `Transfer_Mean_VarScaled` | -0.0014 | -0.0195 | 0.0004 | [-0.1135, +0.0449] | -0.5157 |
| nourreddine_kolf21j | 4611/6148 | 0.7546 | `Transfer_Mean_GlobalScale` | -0.0013 | +nan | nan | [nan, nan] | -0.5452 |
| nourreddine_kolf21j | 4611/6148 | 0.7546 | `Hybrid_on_TMVS` | -0.0096 | +0.0071 | 0.0001 | [-0.1375, +0.1287] | -0.1038 |
| nourreddine_kolf21j | 4611/6148 | 0.7546 | `NegCtrl_HybridProj_SrcCorrTgtVar` | -0.0006 | +0.1117 | 0.0125 | [-0.0215, +0.2439] | -0.0084 |
| nourreddine_kolf21j | 4611/6148 | 0.7546 | `NegCtrl_TMVS_ShuffledVar` | -0.0014 | -0.0195 | 0.0004 | [-0.1135, +0.0450] | -0.5157 |
| replogle_k562_essential | 4611/6148 | 0.8136 | `Hybrid_ProjShared` | -0.0010 | +0.0117 | 0.0001 | [-0.0090, +0.0335] | -0.0187 |
| replogle_k562_essential | 4611/6148 | 0.8136 | `Transfer_Mean_VarScaled` | +0.0000 | +0.0797 | 0.0064 | [+0.0567, +0.1033] | +0.0043 |
| replogle_k562_essential | 4611/6148 | 0.8136 | `Transfer_Mean_GlobalScale` | -0.0001 | +nan | nan | [nan, nan] | -0.5008 |
| replogle_k562_essential | 4611/6148 | 0.8136 | `Hybrid_on_TMVS` | -0.0010 | +0.0134 | 0.0002 | [-0.0075, +0.0353] | -0.0186 |
| replogle_k562_essential | 4611/6148 | 0.8136 | `NegCtrl_HybridProj_SrcCorrTgtVar` | -0.0004 | +0.0146 | 0.0002 | [-0.0019, +0.0317] | -0.0136 |
| replogle_k562_essential | 4611/6148 | 0.8136 | `NegCtrl_TMVS_ShuffledVar` | +0.0000 | +0.0797 | 0.0064 | [+0.0567, +0.1033] | +0.0043 |
| replogle_rpe1 | 4611/6148 | 0.8028 | `Hybrid_ProjShared` | +0.0036 | +0.1796 | 0.0322 | [+0.1551, +0.2060] | +0.1972 |
| replogle_rpe1 | 4611/6148 | 0.8028 | `Transfer_Mean_VarScaled` | -0.0000 | +0.0056 | 0.0000 | [-0.0125, +0.0226] | -0.0554 |
| replogle_rpe1 | 4611/6148 | 0.8028 | `Transfer_Mean_GlobalScale` | -0.0000 | +nan | nan | [nan, nan] | -0.4030 |
| replogle_rpe1 | 4611/6148 | 0.8028 | `Hybrid_on_TMVS` | +0.0036 | +0.1796 | 0.0323 | [+0.1553, +0.2060] | +0.1954 |
| replogle_rpe1 | 4611/6148 | 0.8028 | `NegCtrl_HybridProj_SrcCorrTgtVar` | -0.0002 | -0.0192 | 0.0004 | [-0.0354, -0.0029] | -0.0205 |
| replogle_rpe1 | 4611/6148 | 0.8028 | `NegCtrl_TMVS_ShuffledVar` | -0.0000 | +0.0056 | 0.0000 | [-0.0125, +0.0226] | -0.0555 |
| vcc2025_h1 | 4611/6148 | 0.9748 | `Hybrid_ProjShared` | -0.0006 | -0.0381 | 0.0015 | [-0.0923, +0.0156] | -0.0088 |
| vcc2025_h1 | 4611/6148 | 0.9748 | `Transfer_Mean_VarScaled` | +0.0002 | -0.0056 | 0.0000 | [-0.0661, +0.0738] | +0.1201 |
| vcc2025_h1 | 4611/6148 | 0.9748 | `Transfer_Mean_GlobalScale` | +0.0000 | +nan | nan | [nan, nan] | +0.1793 |
| vcc2025_h1 | 4611/6148 | 0.9748 | `Hybrid_on_TMVS` | -0.0005 | -0.0382 | 0.0015 | [-0.0925, +0.0154] | -0.0065 |
| vcc2025_h1 | 4611/6148 | 0.9748 | `NegCtrl_HybridProj_SrcCorrTgtVar` | -0.0005 | -0.0223 | 0.0005 | [-0.0837, +0.0355] | -0.0323 |
| vcc2025_h1 | 4611/6148 | 0.9748 | `NegCtrl_TMVS_ShuffledVar` | -0.0000 | -0.0123 | 0.0002 | [-0.0560, +0.0350] | -0.0218 |

## Appendix: Alternative H-A4/H-B4 Discrimination & Centroid Conventions (R5)

| Target Line | Method | Ctrl disc_raw (A-T2) | Ctrl disc_raw_panel (A-T3) | Cent disc_raw | Cent disc_raw_panel | Cent de50_r | pds_cosine_panel | Slope OLS (All) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| nadig_hepg2 | `Source_Grand_Mean` | 0.5002 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.2402 |
| nadig_hepg2 | `Transfer_Mean` | 0.9050 | 0.9046 | 0.9257 | 0.9259 | 0.5200 | 0.9015 | +0.0000 |
| nadig_hepg2 | `CIPHER_Linear` | 0.4791 | 0.4823 | 0.4890 | 0.4908 | 0.1000 | 0.4801 | +0.2410 |
| nadig_hepg2 | `Hybrid_Raw` | 0.8954 | 0.8959 | 0.9222 | 0.9225 | 0.5271 | 0.8951 | +0.0050 |
| nadig_hepg2 | `Hybrid_ProjShared` | 0.8995 | 0.8996 | 0.9236 | 0.9234 | 0.5273 | 0.8975 | +0.0025 |
| nadig_hepg2 | `NegCtrl_HybridProj_ShuffledGenes` | 0.9039 | 0.9029 | 0.9243 | 0.9248 | 0.5142 | 0.9004 | +0.0006 |
| nadig_hepg2 | `NegCtrl_HybridProj_DiagCov` | 0.9037 | 0.9029 | 0.9244 | 0.9248 | 0.5177 | 0.9005 | +0.0005 |
| nadig_hepg2 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.9049 | 0.9043 | 0.9257 | 0.9259 | 0.5195 | 0.9013 | +0.0001 |
| nadig_hepg2 | `Transfer_Mean_VarScaled` | 0.9050 | 0.9046 | 0.9257 | 0.9259 | 0.5201 | 0.9015 | +0.0008 |
| nadig_hepg2 | `Transfer_Mean_GlobalScale` | 0.9050 | 0.9046 | 0.9257 | 0.9259 | 0.5201 | 0.9015 | +0.0008 |
| nadig_hepg2 | `NegCtrl_TMVS_ShuffledVar` | 0.9050 | 0.9046 | 0.9257 | 0.9259 | 0.5201 | 0.9015 | +0.0008 |
| nadig_hepg2 | `Hybrid_on_TMVS` | 0.8994 | 0.8996 | 0.9236 | 0.9234 | 0.5274 | 0.8974 | +0.0033 |
| nadig_jurkat | `Source_Grand_Mean` | 0.5001 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.3812 |
| nadig_jurkat | `Transfer_Mean` | 0.8998 | 0.8964 | 0.9410 | 0.9391 | 0.5676 | 0.8960 | +0.0000 |
| nadig_jurkat | `CIPHER_Linear` | 0.4615 | 0.4652 | 0.4952 | 0.4959 | 0.4900 | 0.4746 | +0.4735 |
| nadig_jurkat | `Hybrid_Raw` | 0.9032 | 0.8999 | 0.9423 | 0.9405 | 0.5829 | 0.9002 | +0.0097 |
| nadig_jurkat | `Hybrid_ProjShared` | 0.8989 | 0.8957 | 0.9399 | 0.9383 | 0.5641 | 0.8953 | +0.0017 |
| nadig_jurkat | `NegCtrl_HybridProj_ShuffledGenes` | 0.9012 | 0.8971 | 0.9402 | 0.9385 | 0.5571 | 0.8968 | -0.0014 |
| nadig_jurkat | `NegCtrl_HybridProj_DiagCov` | 0.9014 | 0.8973 | 0.9404 | 0.9388 | 0.5579 | 0.8970 | -0.0013 |
| nadig_jurkat | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.9005 | 0.8971 | 0.9412 | 0.9394 | 0.5669 | 0.8966 | +0.0005 |
| nadig_jurkat | `Transfer_Mean_VarScaled` | 0.8998 | 0.8964 | 0.9410 | 0.9391 | 0.5676 | 0.8960 | -0.0001 |
| nadig_jurkat | `Transfer_Mean_GlobalScale` | 0.8998 | 0.8964 | 0.9410 | 0.9391 | 0.5676 | 0.8960 | -0.0000 |
| nadig_jurkat | `NegCtrl_TMVS_ShuffledVar` | 0.8998 | 0.8964 | 0.9410 | 0.9391 | 0.5676 | 0.8960 | -0.0001 |
| nadig_jurkat | `Hybrid_on_TMVS` | 0.8989 | 0.8957 | 0.9399 | 0.9383 | 0.5640 | 0.8953 | +0.0016 |
| nourreddine_kolf21j | `Source_Grand_Mean` | 0.5000 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.5116 |
| nourreddine_kolf21j | `Transfer_Mean` | 0.7778 | 0.7778 | 0.7680 | 0.7680 | 0.5377 | 0.7876 | +0.0000 |
| nourreddine_kolf21j | `CIPHER_Linear` | 0.5784 | 0.5752 | 0.5131 | 0.5131 | 0.7895 | 0.5752 | +0.6110 |
| nourreddine_kolf21j | `Hybrid_Raw` | 0.7745 | 0.7745 | 0.7386 | 0.7418 | 0.5444 | 0.7810 | +0.0027 |
| nourreddine_kolf21j | `Hybrid_ProjShared` | 0.7745 | 0.7778 | 0.7386 | 0.7418 | 0.5225 | 0.7843 | -0.0013 |
| nourreddine_kolf21j | `NegCtrl_HybridProj_ShuffledGenes` | 0.7255 | 0.7255 | 0.7549 | 0.7549 | 0.5239 | 0.7451 | -0.0035 |
| nourreddine_kolf21j | `NegCtrl_HybridProj_DiagCov` | 0.7288 | 0.7288 | 0.7647 | 0.7614 | 0.5261 | 0.7451 | -0.0030 |
| nourreddine_kolf21j | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.7680 | 0.7680 | 0.7451 | 0.7451 | 0.5346 | 0.7745 | -0.0007 |
| nourreddine_kolf21j | `Transfer_Mean_VarScaled` | 0.7778 | 0.7778 | 0.7680 | 0.7680 | 0.5363 | 0.7876 | -0.0016 |
| nourreddine_kolf21j | `Transfer_Mean_GlobalScale` | 0.7778 | 0.7778 | 0.7680 | 0.7680 | 0.5365 | 0.7876 | -0.0016 |
| nourreddine_kolf21j | `NegCtrl_TMVS_ShuffledVar` | 0.7778 | 0.7778 | 0.7680 | 0.7680 | 0.5363 | 0.7876 | -0.0016 |
| nourreddine_kolf21j | `Hybrid_on_TMVS` | 0.7745 | 0.7778 | 0.7386 | 0.7418 | 0.5211 | 0.7843 | -0.0029 |
| replogle_k562_essential | `Source_Grand_Mean` | 0.5001 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.5228 |
| replogle_k562_essential | `Transfer_Mean` | 0.9095 | 0.9143 | 0.9513 | 0.9521 | 0.5629 | 0.9230 | +0.0000 |
| replogle_k562_essential | `CIPHER_Linear` | 0.4170 | 0.4209 | 0.4899 | 0.4907 | 0.4917 | 0.4296 | +0.6634 |
| replogle_k562_essential | `Hybrid_Raw` | 0.9143 | 0.9180 | 0.9560 | 0.9565 | 0.5829 | 0.9259 | +0.0186 |
| replogle_k562_essential | `Hybrid_ProjShared` | 0.9087 | 0.9128 | 0.9512 | 0.9520 | 0.5530 | 0.9217 | -0.0003 |
| replogle_k562_essential | `NegCtrl_HybridProj_ShuffledGenes` | 0.9148 | 0.9176 | 0.9518 | 0.9524 | 0.5496 | 0.9259 | -0.0025 |
| replogle_k562_essential | `NegCtrl_HybridProj_DiagCov` | 0.9147 | 0.9177 | 0.9518 | 0.9525 | 0.5503 | 0.9259 | -0.0026 |
| replogle_k562_essential | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.9102 | 0.9148 | 0.9514 | 0.9523 | 0.5601 | 0.9235 | -0.0002 |
| replogle_k562_essential | `Transfer_Mean_VarScaled` | 0.9095 | 0.9143 | 0.9512 | 0.9521 | 0.5628 | 0.9230 | -0.0001 |
| replogle_k562_essential | `Transfer_Mean_GlobalScale` | 0.9095 | 0.9143 | 0.9512 | 0.9521 | 0.5628 | 0.9230 | -0.0001 |
| replogle_k562_essential | `NegCtrl_TMVS_ShuffledVar` | 0.9095 | 0.9143 | 0.9512 | 0.9521 | 0.5628 | 0.9230 | -0.0001 |
| replogle_k562_essential | `Hybrid_on_TMVS` | 0.9087 | 0.9128 | 0.9512 | 0.9520 | 0.5529 | 0.9217 | -0.0004 |
| replogle_rpe1 | `Source_Grand_Mean` | 0.5003 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.1874 |
| replogle_rpe1 | `Transfer_Mean` | 0.9418 | 0.9452 | 0.9322 | 0.9345 | 0.4852 | 0.9396 | +0.0000 |
| replogle_rpe1 | `CIPHER_Linear` | 0.5107 | 0.5132 | 0.5017 | 0.5022 | -0.2128 | 0.5137 | +0.1462 |
| replogle_rpe1 | `Hybrid_Raw` | 0.9387 | 0.9426 | 0.9398 | 0.9420 | 0.5043 | 0.9389 | +0.0065 |
| replogle_rpe1 | `Hybrid_ProjShared` | 0.9402 | 0.9444 | 0.9372 | 0.9397 | 0.4991 | 0.9401 | +0.0047 |
| replogle_rpe1 | `NegCtrl_HybridProj_ShuffledGenes` | 0.9419 | 0.9450 | 0.9304 | 0.9330 | 0.4809 | 0.9390 | -0.0015 |
| replogle_rpe1 | `NegCtrl_HybridProj_DiagCov` | 0.9418 | 0.9450 | 0.9304 | 0.9329 | 0.4818 | 0.9390 | -0.0014 |
| replogle_rpe1 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.9421 | 0.9454 | 0.9321 | 0.9344 | 0.4845 | 0.9398 | -0.0004 |
| replogle_rpe1 | `Transfer_Mean_VarScaled` | 0.9418 | 0.9452 | 0.9321 | 0.9345 | 0.4852 | 0.9396 | +0.0000 |
| replogle_rpe1 | `Transfer_Mean_GlobalScale` | 0.9418 | 0.9452 | 0.9321 | 0.9345 | 0.4852 | 0.9396 | +0.0000 |
| replogle_rpe1 | `NegCtrl_TMVS_ShuffledVar` | 0.9418 | 0.9452 | 0.9321 | 0.9345 | 0.4852 | 0.9396 | +0.0000 |
| replogle_rpe1 | `Hybrid_on_TMVS` | 0.9402 | 0.9444 | 0.9372 | 0.9397 | 0.4991 | 0.9401 | +0.0047 |
| vcc2025_h1 | `Source_Grand_Mean` | 0.5000 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.5000 | +0.1550 |
| vcc2025_h1 | `Transfer_Mean` | 0.9109 | 0.9093 | 0.8620 | 0.8609 | 0.4794 | 0.9036 | +0.0000 |
| vcc2025_h1 | `CIPHER_Linear` | 0.4546 | 0.4554 | 0.4777 | 0.4777 | 0.1770 | 0.4472 | +0.1653 |
| vcc2025_h1 | `Hybrid_Raw` | 0.8731 | 0.8728 | 0.8519 | 0.8511 | 0.4608 | 0.8680 | -0.0032 |
| vcc2025_h1 | `Hybrid_ProjShared` | 0.8731 | 0.8728 | 0.8563 | 0.8557 | 0.4571 | 0.8690 | -0.0010 |
| vcc2025_h1 | `NegCtrl_HybridProj_ShuffledGenes` | 0.9082 | 0.9073 | 0.8598 | 0.8590 | 0.4673 | 0.9035 | +0.0002 |
| vcc2025_h1 | `NegCtrl_HybridProj_DiagCov` | 0.9089 | 0.9087 | 0.8619 | 0.8611 | 0.4778 | 0.9040 | +0.0003 |
| vcc2025_h1 | `NegCtrl_HybridProj_SrcCorrTgtVar` | 0.9074 | 0.9066 | 0.8606 | 0.8600 | 0.4779 | 0.9014 | -0.0003 |
| vcc2025_h1 | `Transfer_Mean_VarScaled` | 0.9109 | 0.9095 | 0.8622 | 0.8608 | 0.4794 | 0.9038 | -0.0002 |
| vcc2025_h1 | `Transfer_Mean_GlobalScale` | 0.9109 | 0.9093 | 0.8622 | 0.8609 | 0.4794 | 0.9036 | -0.0001 |
| vcc2025_h1 | `NegCtrl_TMVS_ShuffledVar` | 0.9109 | 0.9092 | 0.8622 | 0.8609 | 0.4793 | 0.9038 | -0.0002 |
| vcc2025_h1 | `Hybrid_on_TMVS` | 0.8733 | 0.8728 | 0.8568 | 0.8555 | 0.4571 | 0.8688 | -0.0012 |
