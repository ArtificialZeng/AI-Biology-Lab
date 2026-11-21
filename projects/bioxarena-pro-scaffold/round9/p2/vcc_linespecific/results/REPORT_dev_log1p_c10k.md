# Line-Specific Zero-Shot Perturbation Evaluation Report (`dev_log1p_c10k`)

- **Evaluated Lines (2)**: replogle_k562_essential, replogle_rpe1
- **Primary Method**: `Hybrid_ProjShared`

## Hypothesis Adjudication Summary

| Hypothesis | Pass / Total | Threshold | Verdict |
|---|---:|---|---|
| **H-A1** (Line-specific slope 95% CI > 0.13) | 0/2 | >= 2 | **未确认** |
| **H-A2** (Noise-calibrated resid corr CI > 0) | 2/2 | >= 2 | **CONFIRMED** |
| **H-A3** (Matched > Mismatched cov, disjoint CI) | 2/2 | >= 2 | **CONFIRMED** |
| **H-A4** (Non-inferiority on disc_raw & de50_r; worst of 4) | 2/2 inf (cent_disc_raw=1/2, cent_de50_r=0/2, ctrl_disc_raw=2/2, ctrl_de50_r=0/2) | <= 1 inf | **未确认** |

## Per-Line Summary by Method

| Target Line | Method | N | Slope OLS (Top50) | 95% CI | Slope NoInt | Slope PerPert | Calib Resid r | Cent disc_raw | Cent de50_r | Ctrl disc_raw | Ctrl de50_r | MSE |
|---|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| replogle_k562_essential | `Source_Grand_Mean` | 1761 | +0.4566 | [+0.4269, +0.4865] | +0.4626 | +0.3726 | +0.4838 | 0.5000 | 0.0000 | 0.5001 | 0.2398 | 0.0066 |
| replogle_k562_essential | `Transfer_Mean` | 1761 | +0.0000 | [+0.0000, +0.0000] | +0.0000 | +0.0000 | +nan | 0.8205 | 0.2794 | 0.7565 | 0.4557 | 0.0127 |
| replogle_k562_essential | `CIPHER_Linear` | 1761 | +0.6424 | [+0.6114, +0.6738] | +0.6615 | +0.5642 | +0.5629 | 0.5002 | 0.8040 | 0.4869 | 0.0526 | 0.0040 |
| replogle_k562_essential | `Hybrid_Raw` | 1761 | +0.0192 | [+0.0162, +0.0221] | +0.0202 | +0.0223 | +0.1409 | 0.8271 | 0.3251 | 0.7576 | 0.4569 | 0.0125 |
| replogle_k562_essential | `Hybrid_ProjShared` | 1761 | +0.0047 | [+0.0035, +0.0061] | +0.0041 | +0.0052 | +0.0524 | 0.8192 | 0.2794 | 0.7502 | 0.4586 | 0.0128 |
| replogle_k562_essential | `NegCtrl_CIPHER_MismatchedCov` | 1761 | +0.6232 | [+0.5923, +0.6547] | +0.6414 | +0.5419 | +0.5526 | 0.4979 | 0.7833 | 0.4953 | 0.0762 | 0.0042 |
| replogle_k562_essential | `NegCtrl_Hybrid_MismatchedCov` | 1761 | +0.0000 | [+0.0000, +0.0000] | +0.0000 | +0.0000 | +nan | 0.8205 | 0.2794 | 0.7565 | 0.4557 | 0.0127 |
| replogle_k562_essential | `NegCtrl_HybridProj_MismatchedCov` | 1761 | +0.0000 | [+0.0000, +0.0000] | +0.0000 | +0.0000 | +nan | 0.8205 | 0.2794 | 0.7565 | 0.4557 | 0.0127 |
| replogle_k562_essential | `NegCtrl_CIPHER_ShuffledGenes` | 1761 | +0.6440 | [+0.6129, +0.6753] | +0.6644 | +0.5657 | +0.5620 | 0.5003 | 0.8096 | 0.5017 | -0.0006 | 0.0041 |
| replogle_k562_essential | `NegCtrl_Hybrid_ShuffledGenes` | 1761 | +0.0209 | [+0.0170, +0.0246] | +0.0230 | +0.0238 | +0.1270 | 0.8316 | 0.3332 | 0.7651 | 0.4494 | 0.0126 |
| replogle_k562_essential | `NegCtrl_HybridProj_ShuffledGenes` | 1761 | -0.0002 | [-0.0017, +0.0014] | -0.0001 | +0.0002 | -0.0017 | 0.8213 | 0.2691 | 0.7530 | 0.4519 | 0.0130 |
| replogle_k562_essential | `NegCtrl_CIPHER_DiagCov` | 1761 | +0.6436 | [+0.6124, +0.6749] | +0.6648 | +0.5652 | +0.5620 | 0.5003 | 0.8102 | 0.5000 | 0.0000 | 0.0039 |
| replogle_k562_essential | `NegCtrl_Hybrid_DiagCov` | 1761 | +0.0204 | [+0.0167, +0.0241] | +0.0234 | +0.0233 | +0.1283 | 0.8316 | 0.3328 | 0.7653 | 0.4498 | 0.0124 |
| replogle_k562_essential | `NegCtrl_HybridProj_DiagCov` | 1761 | -0.0004 | [-0.0019, +0.0011] | -0.0003 | -0.0001 | -0.0041 | 0.8214 | 0.2691 | 0.7532 | 0.4523 | 0.0129 |
| replogle_rpe1 | `Source_Grand_Mean` | 1761 | -0.0456 | [-0.0527, -0.0383] | -0.0435 | -0.0001 | -0.1384 | 0.5000 | 0.0000 | 0.5002 | 0.4963 | 0.0116 |
| replogle_rpe1 | `Transfer_Mean` | 1761 | +0.0000 | [+0.0000, +0.0000] | +0.0000 | +0.0000 | +nan | 0.8470 | 0.3400 | 0.8434 | 0.4321 | 0.0127 |
| replogle_rpe1 | `CIPHER_Linear` | 1761 | -0.0819 | [-0.0895, -0.0740] | -0.0885 | -0.0320 | -0.2265 | 0.5137 | -0.1148 | 0.5154 | 0.2513 | 0.0124 |
| replogle_rpe1 | `Hybrid_Raw` | 1761 | +0.0095 | [+0.0071, +0.0118] | +0.0133 | +0.0130 | +0.1443 | 0.8603 | 0.3794 | 0.8265 | 0.4609 | 0.0125 |
| replogle_rpe1 | `Hybrid_ProjShared` | 1761 | +0.0083 | [+0.0063, +0.0101] | +0.0110 | +0.0114 | +0.1457 | 0.8579 | 0.3759 | 0.8293 | 0.4584 | 0.0125 |
| replogle_rpe1 | `NegCtrl_CIPHER_MismatchedCov` | 1761 | -0.0914 | [-0.0984, -0.0843] | -0.1018 | -0.0450 | -0.2571 | 0.5052 | -0.2403 | 0.4967 | 0.1281 | 0.0126 |
| replogle_rpe1 | `NegCtrl_Hybrid_MismatchedCov` | 1761 | +0.0000 | [+0.0000, +0.0000] | +0.0000 | +0.0000 | +nan | 0.8470 | 0.3400 | 0.8434 | 0.4321 | 0.0127 |
| replogle_rpe1 | `NegCtrl_HybridProj_MismatchedCov` | 1761 | +0.0000 | [+0.0000, +0.0000] | +0.0000 | +0.0000 | +nan | 0.8470 | 0.3400 | 0.8434 | 0.4321 | 0.0127 |
| replogle_rpe1 | `NegCtrl_CIPHER_ShuffledGenes` | 1761 | -0.0954 | [-0.1023, -0.0884] | -0.1069 | -0.0485 | -0.2680 | 0.4998 | -0.2653 | 0.4859 | -0.0039 | 0.0127 |
| replogle_rpe1 | `NegCtrl_Hybrid_ShuffledGenes` | 1761 | -0.0040 | [-0.0045, -0.0035] | -0.0051 | -0.0035 | -0.1253 | 0.8400 | 0.3274 | 0.8482 | 0.4169 | 0.0128 |
| replogle_rpe1 | `NegCtrl_HybridProj_ShuffledGenes` | 1761 | -0.0011 | [-0.0015, -0.0008] | -0.0024 | -0.0014 | -0.0397 | 0.8439 | 0.3318 | 0.8438 | 0.4239 | 0.0128 |
| replogle_rpe1 | `NegCtrl_CIPHER_DiagCov` | 1761 | -0.0955 | [-0.1023, -0.0884] | -0.1087 | -0.0486 | -0.2683 | 0.4999 | -0.2665 | 0.5000 | 0.0000 | 0.0127 |
| replogle_rpe1 | `NegCtrl_Hybrid_DiagCov` | 1761 | -0.0040 | [-0.0045, -0.0035] | -0.0070 | -0.0036 | -0.1438 | 0.8401 | 0.3284 | 0.8482 | 0.4176 | 0.0128 |
| replogle_rpe1 | `NegCtrl_HybridProj_DiagCov` | 1761 | -0.0012 | [-0.0015, -0.0008] | -0.0025 | -0.0013 | -0.0467 | 0.8441 | 0.3329 | 0.8434 | 0.4249 | 0.0127 |
