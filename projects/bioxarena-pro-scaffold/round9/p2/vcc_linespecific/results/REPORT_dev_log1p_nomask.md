# Line-Specific Zero-Shot Perturbation Evaluation Report (`dev_log1p_nomask`)

- **Evaluated Lines (2)**: replogle_k562_essential, replogle_rpe1
- **Primary Method**: `Hybrid_ProjShared`

## Hypothesis Adjudication Summary

| Hypothesis | Pass / Total | Threshold | Verdict |
|---|---:|---|---|
| **H-A1** (Line-specific slope 95% CI > 0.13) | 0/2 | >= 2 | **未确认** |
| **H-A2** (Noise-calibrated resid corr CI > 0) | 2/2 | >= 2 | **CONFIRMED** |
| **H-A3** (Matched > Mismatched cov, disjoint CI) | 2/2 | >= 2 | **CONFIRMED** |
| **H-A4** (Non-inferiority on disc_raw & de50_r; worst of 4) | 2/2 inf (cent_disc_raw=0/2, cent_de50_r=0/2, ctrl_disc_raw=2/2, ctrl_de50_r=0/2) | <= 1 inf | **未确认** |

## Per-Line Summary by Method

| Target Line | Method | N | Slope OLS (Top50) | 95% CI | Slope NoInt | Slope PerPert | Calib Resid r | Cent disc_raw | Cent de50_r | Ctrl disc_raw | Ctrl de50_r | MSE |
|---|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| replogle_k562_essential | `Source_Grand_Mean` | 1761 | +0.4566 | [+0.4269, +0.4865] | +0.4626 | +0.3726 | +0.4838 | 0.5000 | 0.0000 | 0.5001 | 0.2398 | 0.0066 |
| replogle_k562_essential | `Transfer_Mean` | 1761 | +0.0000 | [+0.0000, +0.0000] | +0.0000 | +0.0000 | +nan | 0.8205 | 0.2794 | 0.7565 | 0.4557 | 0.0127 |
| replogle_k562_essential | `CIPHER_Linear` | 1761 | +0.6424 | [+0.6114, +0.6738] | +0.6615 | +0.5641 | +0.5629 | 0.5001 | 0.8041 | 0.4857 | 0.0520 | 0.0040 |
| replogle_k562_essential | `Hybrid_Raw` | 1761 | +0.0191 | [+0.0161, +0.0220] | +0.0200 | +0.0221 | +0.1409 | 0.8272 | 0.3251 | 0.7577 | 0.4567 | 0.0125 |
| replogle_k562_essential | `Hybrid_ProjShared` | 1761 | +0.0046 | [+0.0033, +0.0059] | +0.0040 | +0.0050 | +0.0513 | 0.8193 | 0.2794 | 0.7501 | 0.4584 | 0.0128 |
| replogle_k562_essential | `NegCtrl_CIPHER_MismatchedCov` | 1761 | +0.6233 | [+0.5924, +0.6548] | +0.6415 | +0.5421 | +0.5526 | 0.4979 | 0.7835 | 0.4949 | 0.0755 | 0.0042 |
| replogle_k562_essential | `NegCtrl_Hybrid_MismatchedCov` | 1761 | +0.0000 | [+0.0000, +0.0000] | +0.0000 | +0.0000 | +nan | 0.8205 | 0.2794 | 0.7565 | 0.4557 | 0.0127 |
| replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov` | 1761 | +0.0000 | [+0.0000, +0.0000] | +0.0000 | +0.0000 | +nan | 0.8205 | 0.2794 | 0.7565 | 0.4557 | 0.0127 |
| replogle_k562_essential | `NegCtrl_CIPHER_ShuffledGenes` | 1761 | +0.6440 | [+0.6129, +0.6753] | +0.6644 | +0.5657 | +0.5620 | 0.5003 | 0.8096 | 0.5033 | -0.0012 | 0.0040 |
| replogle_k562_essential | `NegCtrl_Hybrid_ShuffledGenes` | 1761 | +0.0207 | [+0.0169, +0.0244] | +0.0229 | +0.0236 | +0.1268 | 0.8316 | 0.3331 | 0.7651 | 0.4494 | 0.0126 |
| replogle_k562_essential | `NegCtrl_HybridProj_ShuffledGenes` | 1761 | -0.0002 | [-0.0018, +0.0013] | -0.0002 | +0.0001 | -0.0022 | 0.8213 | 0.2692 | 0.7530 | 0.4519 | 0.0130 |
| replogle_k562_essential | `NegCtrl_CIPHER_DiagCov` | 1761 | +0.6436 | [+0.6124, +0.6749] | +0.6648 | +0.5652 | +0.5620 | 0.5003 | 0.8102 | 0.5000 | 0.0000 | 0.0039 |
| replogle_k562_essential | `NegCtrl_Hybrid_DiagCov` | 1761 | +0.0203 | [+0.0165, +0.0240] | +0.0233 | +0.0232 | +0.1281 | 0.8316 | 0.3328 | 0.7652 | 0.4498 | 0.0124 |
| replogle_k562_essential | `NegCtrl_HybridProj_DiagCov` | 1761 | -0.0005 | [-0.0020, +0.0010] | -0.0004 | -0.0002 | -0.0046 | 0.8214 | 0.2692 | 0.7533 | 0.4523 | 0.0128 |
| replogle_rpe1 | `Source_Grand_Mean` | 1761 | -0.0456 | [-0.0527, -0.0383] | -0.0435 | -0.0001 | -0.1384 | 0.5000 | 0.0000 | 0.5002 | 0.4963 | 0.0116 |
| replogle_rpe1 | `Transfer_Mean` | 1761 | +0.0000 | [+0.0000, +0.0000] | +0.0000 | +0.0000 | +nan | 0.8470 | 0.3400 | 0.8434 | 0.4321 | 0.0127 |
| replogle_rpe1 | `CIPHER_Linear` | 1761 | -0.0788 | [-0.0868, -0.0707] | -0.0840 | -0.0283 | -0.2160 | 0.5187 | -0.0816 | 0.5150 | 0.2506 | 0.0123 |
| replogle_rpe1 | `Hybrid_Raw` | 1761 | +0.0118 | [+0.0088, +0.0145] | +0.0163 | +0.0159 | +0.1458 | 0.8605 | 0.3854 | 0.8215 | 0.4638 | 0.0125 |
| replogle_rpe1 | `Hybrid_ProjShared` | 1761 | +0.0101 | [+0.0078, +0.0124] | +0.0134 | +0.0139 | +0.1466 | 0.8588 | 0.3814 | 0.8252 | 0.4614 | 0.0125 |
| replogle_rpe1 | `NegCtrl_CIPHER_MismatchedCov` | 1761 | -0.0906 | [-0.0976, -0.0834] | -0.1003 | -0.0442 | -0.2546 | 0.5062 | -0.2318 | 0.4954 | 0.1282 | 0.0125 |
| replogle_rpe1 | `NegCtrl_Hybrid_MismatchedCov` | 1761 | +0.0000 | [+0.0000, +0.0000] | +0.0000 | +0.0000 | +nan | 0.8470 | 0.3400 | 0.8434 | 0.4321 | 0.0127 |
| replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov` | 1761 | +0.0000 | [+0.0000, +0.0000] | +0.0000 | +0.0000 | +nan | 0.8470 | 0.3400 | 0.8434 | 0.4321 | 0.0127 |
| replogle_rpe1 | `NegCtrl_CIPHER_ShuffledGenes` | 1761 | -0.0954 | [-0.1023, -0.0884] | -0.1065 | -0.0485 | -0.2679 | 0.4998 | -0.2643 | 0.4868 | -0.0044 | 0.0127 |
| replogle_rpe1 | `NegCtrl_Hybrid_ShuffledGenes` | 1761 | -0.0048 | [-0.0054, -0.0042] | -0.0062 | -0.0043 | -0.1239 | 0.8380 | 0.3245 | 0.8485 | 0.4131 | 0.0129 |
| replogle_rpe1 | `NegCtrl_HybridProj_ShuffledGenes` | 1761 | -0.0014 | [-0.0018, -0.0010] | -0.0029 | -0.0017 | -0.0398 | 0.8431 | 0.3297 | 0.8437 | 0.4216 | 0.0128 |
| replogle_rpe1 | `NegCtrl_CIPHER_DiagCov` | 1761 | -0.0955 | [-0.1023, -0.0884] | -0.1087 | -0.0486 | -0.2683 | 0.4999 | -0.2665 | 0.5000 | 0.0000 | 0.0127 |
| replogle_rpe1 | `NegCtrl_Hybrid_DiagCov` | 1761 | -0.0049 | [-0.0055, -0.0043] | -0.0084 | -0.0043 | -0.1414 | 0.8382 | 0.3258 | 0.8487 | 0.4140 | 0.0128 |
| replogle_rpe1 | `NegCtrl_HybridProj_DiagCov` | 1761 | -0.0014 | [-0.0018, -0.0010] | -0.0030 | -0.0016 | -0.0462 | 0.8434 | 0.3313 | 0.8432 | 0.4230 | 0.0127 |
