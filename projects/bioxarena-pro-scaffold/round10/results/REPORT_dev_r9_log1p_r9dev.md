# Round 10 Confirmatory Group C (`FCT`) Report (`dev_r9_log1p_r9dev`)

- **Phase**: `dev_r9` | **PREREG_C Commit**: `` | **Code Freeze**: ``
- **Bootstrap**: `np.random.default_rng(20261122)` x `4000` (Bonferroni two-sided `97.5%` CIs)
- **Cohort Disclosure**: **Insufficient confirmatory cell lines (K_conf < 2); confirmatory hypothesis testing was not performed (确认性线不足，未做确认性检验).** (`[]`)
- **Secondary Lines**: `['nadig_hepg2', 'nadig_jurkat', 'nourreddine_kolf21j', 'replogle_k562_essential', 'replogle_rpe1', 'vcc2025_h1']`
- **Primary Group C (`FCT`) Verdict**: **NOT_EVALUATED_INSUFFICIENT_LINES**
- **Parallel Practical Utility (`H-C3` Non-Inferiority)**: **NOT_EVALUATED_INSUFFICIENT_LINES**

## Per-Line Results (`resp_top50`, Joint-Bootstrap 97.5% CIs)

| Target Line | Role | Method | N | r_orth_cal | 97.5% CI | Slope OLS | 97.5% CI | r_cal (raw) | rho_AB_orth | pds_cos |
|---|---|---|---:|---:|---|---:|---|---:|---:|---:|
| nadig_hepg2 | Secondary | `Source_Grand_Mean` | 797 | +0.3108 | [+0.2732, +0.3482] | +0.0200 | [-0.0099, +0.0528] | +0.0299 | +0.8742 | 0.5000 |
| nadig_hepg2 | Secondary | `Transfer_Mean` | 797 | +nan | [+nan, +nan] | +0.0000 | [+0.0000, +0.0000] | +nan | +0.8742 | 0.8790 |
| nadig_hepg2 | Secondary | `FCT` | 797 | +0.2614 | [+0.2424, +0.2812] | +0.0587 | [+0.0537, +0.0640] | +0.2620 | +0.8742 | 0.8852 |
| nadig_hepg2 | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 797 | +0.0457 | [+0.0330, +0.0595] | +0.0084 | [+0.0055, +0.0117] | +0.0498 | +0.8742 | 0.8719 |
| nadig_hepg2 | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 797 | +0.0552 | [+0.0415, +0.0687] | +0.0025 | [-0.0029, +0.0085] | +0.0102 | +0.8742 | 0.8750 |
| nadig_hepg2 | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 797 | +0.0665 | [+0.0527, +0.0801] | +0.0087 | [+0.0053, +0.0124] | +0.0518 | +0.8742 | 0.8740 |
| nadig_hepg2 | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 797 | -0.0044 | [-0.0233, +0.0135] | +0.0042 | [+0.0019, +0.0065] | +0.0349 | +0.8742 | 0.8767 |
| nadig_hepg2 | Secondary | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 797 | +0.1187 | [+0.1060, +0.1317] | +0.0281 | [+0.0149, +0.0435] | +0.0654 | +0.8742 | 0.8776 |
| nadig_hepg2 | Secondary | `NegCtrl_FCT_Perm0` | 797 | +0.1614 | [+0.1471, +0.1746] | -0.0145 | [-0.0363, +0.0087] | -0.0245 | +0.8742 | 0.8930 |
| nadig_hepg2 | Secondary | `Transfer_Mean_GlobalScale` | 797 | +nan | [+nan, +nan] | -0.0084 | [-0.0110, -0.0055] | -0.1300 | +0.8742 | 0.8790 |
| nadig_hepg2 | Secondary | `Transfer_Mean_VarScaled` | 797 | +0.0419 | [+0.0244, +0.0603] | +0.0356 | [+0.0268, +0.0440] | +0.1221 | +0.8742 | 0.8736 |
| nadig_hepg2 | Secondary | `CIPHER_Linear` | 797 | +0.1171 | [+0.0849, +0.1498] | -0.0664 | [-0.0994, -0.0295] | -0.0823 | +0.8742 | 0.4840 |
| nadig_hepg2 | Secondary | `Hybrid_ProjShared` | 797 | +0.0484 | [+0.0226, +0.0753] | +0.0099 | [+0.0061, +0.0137] | +0.0727 | +0.8742 | 0.8675 |
| nadig_hepg2 | Secondary | `Hybrid_on_TMVS` | 797 | +0.0548 | [+0.0342, +0.0761] | +0.0454 | [+0.0347, +0.0559] | +0.1288 | +0.8742 | 0.8614 |
| nadig_jurkat | Secondary | `Source_Grand_Mean` | 931 | +0.0114 | [-0.0160, +0.0410] | -0.0048 | [-0.0367, +0.0309] | -0.0066 | +0.8467 | 0.5000 |
| nadig_jurkat | Secondary | `Transfer_Mean` | 931 | +nan | [+nan, +nan] | +0.0000 | [+0.0000, +0.0000] | +nan | +0.8467 | 0.8534 |
| nadig_jurkat | Secondary | `FCT` | 931 | +0.1047 | [+0.0860, +0.1243] | +0.0201 | [+0.0148, +0.0256] | +0.0919 | +0.8467 | 0.8362 |
| nadig_jurkat | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 931 | +0.0644 | [+0.0498, +0.0791] | +0.0140 | [+0.0105, +0.0178] | +0.0806 | +0.8467 | 0.8559 |
| nadig_jurkat | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 931 | +0.0514 | [+0.0390, +0.0642] | +0.0170 | [+0.0119, +0.0222] | +0.0667 | +0.8467 | 0.8440 |
| nadig_jurkat | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 931 | +0.0703 | [+0.0557, +0.0846] | +0.0155 | [+0.0119, +0.0191] | +0.0889 | +0.8467 | 0.8413 |
| nadig_jurkat | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 931 | -0.0492 | [-0.0646, -0.0341] | -0.0101 | [-0.0129, -0.0075] | -0.0736 | +0.8467 | 0.8434 |
| nadig_jurkat | Secondary | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 931 | +0.0518 | [+0.0354, +0.0678] | +0.0145 | [-0.0048, +0.0369] | +0.0307 | +0.8467 | 0.8326 |
| nadig_jurkat | Secondary | `NegCtrl_FCT_Perm0` | 931 | +0.0733 | [+0.0598, +0.0870] | +0.0280 | [+0.0048, +0.0532] | +0.0488 | +0.8467 | 0.8504 |
| nadig_jurkat | Secondary | `Transfer_Mean_GlobalScale` | 931 | +nan | [+nan, +nan] | -0.0009 | [-0.0070, +0.0048] | -0.0063 | +0.8467 | 0.8534 |
| nadig_jurkat | Secondary | `Transfer_Mean_VarScaled` | 931 | +0.0941 | [+0.0800, +0.1092] | +0.0093 | [+0.0028, +0.0155] | +0.0432 | +0.8467 | 0.8518 |
| nadig_jurkat | Secondary | `CIPHER_Linear` | 931 | +0.0444 | [+0.0122, +0.0772] | +0.0071 | [-0.0316, +0.0467] | +0.0081 | +0.8467 | 0.4948 |
| nadig_jurkat | Secondary | `Hybrid_ProjShared` | 931 | +0.0719 | [+0.0460, +0.0971] | +0.0069 | [+0.0039, +0.0101] | +0.0663 | +0.8467 | 0.8503 |
| nadig_jurkat | Secondary | `Hybrid_on_TMVS` | 931 | +0.1140 | [+0.0957, +0.1329] | +0.0162 | [+0.0095, +0.0227] | +0.0667 | +0.8467 | 0.8493 |
| nourreddine_kolf21j | Secondary | `Source_Grand_Mean` | 18 | +0.0639 | [-0.0611, +0.1917] | +0.2720 | [+0.0372, +0.5522] | +0.3472 | +0.8195 | 0.5000 |
| nourreddine_kolf21j | Secondary | `Transfer_Mean` | 18 | +nan | [+nan, +nan] | +0.0000 | [+0.0000, +0.0000] | +nan | +0.8195 | 0.8072 |
| nourreddine_kolf21j | Secondary | `FCT` | 18 | +0.1839 | [+0.0183, +0.3537] | +0.0747 | [+0.0053, +0.1597] | +0.2547 | +0.8195 | 0.7810 |
| nourreddine_kolf21j | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 18 | +0.1081 | [+0.0320, +0.1671] | +0.0656 | [+0.0172, +0.1186] | +0.2460 | +0.8195 | 0.7876 |
| nourreddine_kolf21j | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 18 | +0.0633 | [-0.0556, +0.1439] | +0.0359 | [-0.0088, +0.0862] | +0.1300 | +0.8195 | 0.7745 |
| nourreddine_kolf21j | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 18 | +0.0870 | [+0.0104, +0.1472] | +0.0443 | [+0.0056, +0.1055] | +0.1609 | +0.8195 | 0.7712 |
| nourreddine_kolf21j | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 18 | +0.0528 | [-0.0635, +0.1160] | +0.0605 | [+0.0065, +0.1139] | +0.2007 | +0.8195 | 0.7810 |
| nourreddine_kolf21j | Secondary | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 18 | +0.0708 | [-0.0218, +0.1810] | +0.0300 | [+0.0024, +0.0597] | +0.0889 | +0.8195 | 0.7843 |
| nourreddine_kolf21j | Secondary | `NegCtrl_FCT_Perm0` | 18 | +0.0240 | [-0.0448, +0.1102] | +0.1273 | [+0.0103, +0.2504] | +0.2454 | +0.8195 | 0.7745 |
| nourreddine_kolf21j | Secondary | `Transfer_Mean_GlobalScale` | 18 | +nan | [+nan, +nan] | -0.0686 | [-0.1361, -0.0089] | -0.3304 | +0.8195 | 0.8072 |
| nourreddine_kolf21j | Secondary | `Transfer_Mean_VarScaled` | 18 | +0.1385 | [-0.0493, +0.2910] | -0.0133 | [-0.0922, +0.0565] | -0.0392 | +0.8195 | 0.7941 |
| nourreddine_kolf21j | Secondary | `CIPHER_Linear` | 18 | -0.1096 | [-0.2640, +0.0507] | +0.2738 | [+0.0192, +0.5483] | +0.3147 | +0.8195 | 0.5098 |
| nourreddine_kolf21j | Secondary | `Hybrid_ProjShared` | 18 | +0.0093 | [-0.1312, +0.1342] | -0.0137 | [-0.0577, +0.0258] | -0.0605 | +0.8195 | 0.7974 |
| nourreddine_kolf21j | Secondary | `Hybrid_on_TMVS` | 18 | +0.1155 | [-0.0524, +0.2727] | -0.0270 | [-0.1245, +0.0585] | -0.0633 | +0.8195 | 0.7876 |
| replogle_k562_essential | Secondary | `Source_Grand_Mean` | 858 | -0.0665 | [-0.0904, -0.0439] | +0.2575 | [+0.2086, +0.3107] | +0.2541 | +0.8859 | 0.5000 |
| replogle_k562_essential | Secondary | `Transfer_Mean` | 858 | +nan | [+nan, +nan] | +0.0000 | [+0.0000, +0.0000] | +nan | +0.8859 | 0.8791 |
| replogle_k562_essential | Secondary | `FCT` | 858 | +0.1939 | [+0.1739, +0.2147] | +0.0419 | [+0.0354, +0.0489] | +0.1469 | +0.8859 | 0.8612 |
| replogle_k562_essential | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 858 | +0.0487 | [+0.0303, +0.0741] | +0.0084 | [+0.0034, +0.0146] | +0.0406 | +0.8859 | 0.8739 |
| replogle_k562_essential | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 858 | +0.0493 | [+0.0302, +0.0727] | +0.0079 | [+0.0012, +0.0151] | +0.0363 | +0.8859 | 0.8599 |
| replogle_k562_essential | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 858 | -0.0100 | [-0.0265, +0.0081] | +0.0247 | [+0.0154, +0.0349] | +0.0737 | +0.8859 | 0.8648 |
| replogle_k562_essential | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 858 | -0.0509 | [-0.0769, -0.0139] | -0.0089 | [-0.0136, -0.0027] | -0.0528 | +0.8859 | 0.8697 |
| replogle_k562_essential | Secondary | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 858 | +0.1086 | [+0.0921, +0.1257] | +0.1509 | [+0.1236, +0.1800] | +0.2357 | +0.8859 | 0.8581 |
| replogle_k562_essential | Secondary | `NegCtrl_FCT_Perm0` | 858 | +0.1458 | [+0.1313, +0.1603] | +0.2535 | [+0.2194, +0.2899] | +0.3118 | +0.8859 | 0.8735 |
| replogle_k562_essential | Secondary | `Transfer_Mean_GlobalScale` | 858 | +nan | [+nan, +nan] | +0.0029 | [+0.0024, +0.0034] | +0.2852 | +0.8859 | 0.8791 |
| replogle_k562_essential | Secondary | `Transfer_Mean_VarScaled` | 858 | +0.1822 | [+0.1603, +0.2050] | +0.0605 | [+0.0513, +0.0708] | +0.2452 | +0.8859 | 0.8804 |
| replogle_k562_essential | Secondary | `CIPHER_Linear` | 858 | +0.0252 | [-0.0005, +0.0513] | +0.3501 | [+0.2936, +0.4117] | +0.2846 | +0.8859 | 0.4610 |
| replogle_k562_essential | Secondary | `Hybrid_ProjShared` | 858 | +0.0353 | [+0.0098, +0.0613] | +0.0021 | [-0.0009, +0.0049] | +0.0169 | +0.8859 | 0.8750 |
| replogle_k562_essential | Secondary | `Hybrid_on_TMVS` | 858 | +0.1733 | [+0.1491, +0.1977] | +0.0627 | [+0.0528, +0.0732] | +0.2246 | +0.8859 | 0.8761 |
| replogle_rpe1 | Secondary | `Source_Grand_Mean` | 1068 | +0.5247 | [+0.5041, +0.5449] | -0.0689 | [-0.0864, -0.0515] | -0.1456 | +0.9134 | 0.5000 |
| replogle_rpe1 | Secondary | `Transfer_Mean` | 1068 | +nan | [+nan, +nan] | +0.0000 | [+0.0000, +0.0000] | +nan | +0.9134 | 0.9286 |
| replogle_rpe1 | Secondary | `FCT` | 1068 | +0.1295 | [+0.1116, +0.1480] | +0.0228 | [+0.0201, +0.0255] | +0.1595 | +0.9134 | 0.9290 |
| replogle_rpe1 | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 1068 | +0.0885 | [+0.0747, +0.1019] | +0.0128 | [+0.0111, +0.0145] | +0.1265 | +0.9134 | 0.9310 |
| replogle_rpe1 | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 1068 | +0.0057 | [-0.0101, +0.0214] | +0.0045 | [+0.0024, +0.0065] | +0.0407 | +0.9134 | 0.9256 |
| replogle_rpe1 | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 1068 | +0.0093 | [-0.0030, +0.0219] | -0.0126 | [-0.0155, -0.0100] | -0.0891 | +0.9134 | 0.9284 |
| replogle_rpe1 | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 1068 | -0.0051 | [-0.0175, +0.0068] | -0.0031 | [-0.0047, -0.0014] | -0.0318 | +0.9134 | 0.9256 |
| replogle_rpe1 | Secondary | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 1068 | +0.0615 | [+0.0493, +0.0740] | -0.0264 | [-0.0326, -0.0207] | -0.0979 | +0.9134 | 0.9234 |
| replogle_rpe1 | Secondary | `NegCtrl_FCT_Perm0` | 1068 | +0.0381 | [+0.0267, +0.0489] | -0.1021 | [-0.1115, -0.0926] | -0.2814 | +0.9134 | 0.9303 |
| replogle_rpe1 | Secondary | `Transfer_Mean_GlobalScale` | 1068 | +nan | [+nan, +nan] | +0.0084 | [+0.0077, +0.0092] | +0.3596 | +0.9134 | 0.9286 |
| replogle_rpe1 | Secondary | `Transfer_Mean_VarScaled` | 1068 | +0.0922 | [+0.0734, +0.1111] | +0.0323 | [+0.0286, +0.0361] | +0.2294 | +0.9134 | 0.9250 |
| replogle_rpe1 | Secondary | `CIPHER_Linear` | 1068 | +0.3425 | [+0.3111, +0.3729] | -0.1484 | [-0.1673, -0.1297] | -0.2750 | +0.9134 | 0.5161 |
| replogle_rpe1 | Secondary | `Hybrid_ProjShared` | 1068 | +0.2088 | [+0.1814, +0.2357] | +0.0147 | [+0.0126, +0.0168] | +0.1842 | +0.9134 | 0.9231 |
| replogle_rpe1 | Secondary | `Hybrid_on_TMVS` | 1068 | +0.1840 | [+0.1620, +0.2061] | +0.0470 | [+0.0426, +0.0515] | +0.2740 | +0.9134 | 0.9202 |
| vcc2025_h1 | Secondary | `Source_Grand_Mean` | 80 | +0.1086 | [+0.0504, +0.1699] | +0.1892 | [-0.0125, +0.3778] | +0.2985 | +0.9885 | 0.5000 |
| vcc2025_h1 | Secondary | `Transfer_Mean` | 80 | +nan | [+nan, +nan] | +0.0000 | [+0.0000, +0.0000] | +nan | +0.9885 | 0.8758 |
| vcc2025_h1 | Secondary | `FCT` | 80 | +0.3505 | [+0.2765, +0.4225] | +0.2503 | [+0.1371, +0.3541] | +0.5489 | +0.9885 | 0.8998 |
| vcc2025_h1 | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 80 | +0.1240 | [+0.0788, +0.1661] | +0.0095 | [+0.0003, +0.0202] | +0.0603 | +0.9885 | 0.8848 |
| vcc2025_h1 | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 80 | +0.0550 | [+0.0115, +0.0974] | +0.0076 | [+0.0019, +0.0136] | +0.0565 | +0.9885 | 0.8691 |
| vcc2025_h1 | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 80 | +0.1453 | [+0.0900, +0.2022] | +0.0374 | [+0.0201, +0.0532] | +0.1905 | +0.9885 | 0.8826 |
| vcc2025_h1 | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 80 | +0.0536 | [+0.0136, +0.0921] | +0.0013 | [-0.0067, +0.0108] | +0.0104 | +0.9885 | 0.8687 |
| vcc2025_h1 | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 80 | -0.0541 | [-0.1092, +0.0006] | -0.0086 | [-0.0194, -0.0001] | -0.0700 | +0.9885 | 0.8786 |
| vcc2025_h1 | Secondary | `NegCtrl_FCT_Perm0` | 80 | +0.1802 | [+0.1328, +0.2279] | +0.1780 | [+0.0547, +0.2926] | +0.4000 | +0.9885 | 0.8965 |
| vcc2025_h1 | Secondary | `Transfer_Mean_GlobalScale` | 80 | +nan | [+nan, +nan] | +0.0683 | [-0.0199, +0.1513] | +0.2506 | +0.9885 | 0.8758 |
| vcc2025_h1 | Secondary | `Transfer_Mean_VarScaled` | 80 | -0.0772 | [-0.1890, +0.0299] | -0.1647 | [-0.2639, -0.0582] | -0.4447 | +0.9885 | 0.8750 |
| vcc2025_h1 | Secondary | `CIPHER_Linear` | 80 | -0.0468 | [-0.1398, +0.0367] | +0.1589 | [-0.0552, +0.3601] | +0.2317 | +0.9885 | 0.4663 |
| vcc2025_h1 | Secondary | `Hybrid_ProjShared` | 80 | -0.0403 | [-0.1280, +0.0377] | -0.0089 | [-0.0249, +0.0048] | -0.0496 | +0.9885 | 0.8410 |
| vcc2025_h1 | Secondary | `Hybrid_on_TMVS` | 80 | -0.0819 | [-0.1820, +0.0140] | -0.1735 | [-0.2755, -0.0654] | -0.4152 | +0.9885 | 0.8146 |

## 100-Perm `NegCtrl_FCT_Perm` (`k=0..99`) & Mismatched Mean Summary

- **nadig_hepg2**: perm100 `r_orth_cal` median=+0.1264, p97.5=+0.1696, p98.75=+0.1744; max mismatched hi=+0.1317
- **nadig_jurkat**: perm100 `r_orth_cal` median=+0.0685, p97.5=+0.1203, p98.75=+0.1285; max mismatched hi=+0.0846
- **nourreddine_kolf21j**: perm100 `r_orth_cal` median=+0.0577, p97.5=+0.1246, p98.75=+0.1381; max mismatched hi=+0.1810
- **replogle_k562_essential**: perm100 `r_orth_cal` median=+0.1157, p97.5=+0.1563, p98.75=+0.1576; max mismatched hi=+0.1257
- **replogle_rpe1**: perm100 `r_orth_cal` median=+0.0758, p97.5=+0.1064, p98.75=+0.1083; max mismatched hi=+0.1019
- **vcc2025_h1**: perm100 `r_orth_cal` median=+0.2101, p97.5=+0.2445, p98.75=+0.2521; max mismatched hi=+0.2022

## `SA-3` Round 9 Anchors Notice (`P7`) & `SA-6 / S7` Cell-Count Controls (`P6`)

- *Benchmark Anchor Note (`P7`):* `Hybrid_ProjShared` follows frozen Round 9 commit `1ca9f60f` 3D `span{1, b_S, a_S}` projection (whereas Round 9 `PREREG_A v4` prose described rank-1 `q_1`).
- **nadig_hepg2 S7**: Spearman(`r_orth_p`, ln n_cells)=+0.2084 [+0.1324, +0.2857], Spearman(`r_orth_p`, |k_p|)=-0.0154 [-0.0970, +0.0647], OLS-adjusted intercept=+0.2943 [+0.2763, +0.3113]
  * `T1_low` (N=266, n_cells=[25..48]): FCT `r_orth_cal`=+0.2117 [+0.1774, +0.2450], max_mis_hi=+0.1025, perm100_p975=+0.1474, HC1=True, HC2=True
  * `T2_mid` (N=266, n_cells=[48..75]): FCT `r_orth_cal`=+0.2589 [+0.2264, +0.2928], max_mis_hi=+0.1465, perm100_p975=+0.1847, HC1=True, HC2=True
  * `T3_high` (N=265, n_cells=[75..1213]): FCT `r_orth_cal`=+0.3333 [+0.3059, +0.3613], max_mis_hi=+0.1842, perm100_p975=+0.2027, HC1=True, HC2=True
- **nadig_jurkat S7**: Spearman(`r_orth_p`, ln n_cells)=-0.0683 [-0.1360, +0.0039], Spearman(`r_orth_p`, |k_p|)=-0.0950 [-0.1669, -0.0234], OLS-adjusted intercept=+0.1174 [+0.1021, +0.1331]
  * `T1_low` (N=311, n_cells=[25..69]): FCT `r_orth_cal`=+0.1144 [+0.0860, +0.1432], max_mis_hi=+0.0823, perm100_p975=+0.1194, HC1=True, HC2=False
  * `T2_mid` (N=310, n_cells=[69..114]): FCT `r_orth_cal`=+0.1263 [+0.0944, +0.1586], max_mis_hi=+0.1229, perm100_p975=+0.1146, HC1=True, HC2=False
  * `T3_high` (N=310, n_cells=[114..2555]): FCT `r_orth_cal`=+0.0674 [+0.0326, +0.1023], max_mis_hi=+0.0997, perm100_p975=+0.1314, HC1=True, HC2=False
- **nourreddine_kolf21j S7**: Spearman(`r_orth_p`, ln n_cells)=-0.1280 [-0.6819, +0.5794], Spearman(`r_orth_p`, |k_p|)=-0.4262 [-0.7634, +0.1441], OLS-adjusted intercept=+0.1490 [+0.0128, +0.2742]
  * `T1_low` (N=6, n_cells=[46..66]): FCT `r_orth_cal`=+0.2785 [+0.0304, +0.5342], max_mis_hi=+0.3162, perm100_p975=+0.1818, HC1=True, HC2=False
  * `T2_mid` (N=6, n_cells=[71..87]): FCT `r_orth_cal`=+0.0730 [-0.0360, +0.2265], max_mis_hi=+0.1493, perm100_p975=+0.0902, HC1=False, HC2=False
  * `T3_high` (N=6, n_cells=[89..148]): FCT `r_orth_cal`=+0.1593 [-0.1218, +0.4864], max_mis_hi=+0.3102, perm100_p975=+0.1767, HC1=False, HC2=False
- **replogle_k562_essential S7**: Spearman(`r_orth_p`, ln n_cells)=+0.1927 [+0.1194, +0.2646], Spearman(`r_orth_p`, |k_p|)=-0.0643 [-0.1410, +0.0105], OLS-adjusted intercept=+0.1943 [+0.1800, +0.2085]
  * `T1_low` (N=286, n_cells=[25..98]): FCT `r_orth_cal`=+0.1682 [+0.1348, +0.2007], max_mis_hi=+0.1151, perm100_p975=+0.1450, HC1=True, HC2=False
  * `T2_mid` (N=286, n_cells=[98..169]): FCT `r_orth_cal`=+0.1889 [+0.1472, +0.2393], max_mis_hi=+0.1333, perm100_p975=+0.1607, HC1=True, HC2=False
  * `T3_high` (N=286, n_cells=[169..992]): FCT `r_orth_cal`=+0.2235 [+0.1971, +0.2502], max_mis_hi=+0.1614, perm100_p975=+0.1696, HC1=True, HC2=True
- **replogle_rpe1 S7**: Spearman(`r_orth_p`, ln n_cells)=-0.0516 [-0.1203, +0.0171], Spearman(`r_orth_p`, |k_p|)=+0.0184 [-0.0518, +0.0842], OLS-adjusted intercept=+0.1887 [+0.1739, +0.2040]
  * `T1_low` (N=356, n_cells=[25..51]): FCT `r_orth_cal`=+0.1303 [+0.1055, +0.1556], max_mis_hi=+0.1215, perm100_p975=+0.1163, HC1=True, HC2=False
  * `T2_mid` (N=356, n_cells=[51..86]): FCT `r_orth_cal`=+0.1331 [+0.0961, +0.1704], max_mis_hi=+0.1201, perm100_p975=+0.1124, HC1=True, HC2=False
  * `T3_high` (N=356, n_cells=[86..3580]): FCT `r_orth_cal`=+0.1266 [+0.0964, +0.1595], max_mis_hi=+0.0964, perm100_p975=+0.0916, HC1=True, HC2=False
- **vcc2025_h1 S7**: Spearman(`r_orth_p`, ln n_cells)=-0.0843 [-0.3351, +0.1858], Spearman(`r_orth_p`, |k_p|)=+0.0414 [-0.2020, +0.2815], OLS-adjusted intercept=+0.3196 [+0.2587, +0.3818]
  * `T1_low` (N=27, n_cells=[33..407]): FCT `r_orth_cal`=+0.3537 [+0.2101, +0.4819], max_mis_hi=+0.1886, perm100_p975=+0.2849, HC1=True, HC2=False
  * `T2_mid` (N=27, n_cells=[450..959]): FCT `r_orth_cal`=+0.3867 [+0.2724, +0.4970], max_mis_hi=+0.3210, perm100_p975=+0.2449, HC1=True, HC2=False
  * `T3_high` (N=26, n_cells=[972..2428]): FCT `r_orth_cal`=+0.3040 [+0.2080, +0.4025], max_mis_hi=+0.2138, perm100_p975=+0.2407, HC1=True, HC2=False
