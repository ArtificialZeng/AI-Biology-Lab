# 探索性分析（未预注册 / Post-Hoc Exploratory）：FCT 特异性对照检验

> **重要声明（PREREG_A v4 合规声明）**：**本文件为未预注册的事后探索性分析（Exploratory / Not Pre-Registered），不改变任何预注册确认性假设（Group A 与 Group B 均维持 `未确认`）及主表结果。**

- **Code Freeze Commit**: `1ca9f60f49adf1ac7cce8266e40533273411db30` | **Bootstrap**: `np.random.default_rng(20261121)` x 4,000 | **Permutations**: `np.random.default_rng(k)` (`k = 0..99`, exact `COMMON_GENES = 6148` order)
- **分析目的**：主评测（R15 次要基线）观察到 Fold-Change Transfer (`FCT`) 在 3 条确认性线上取得较高 `r_orth_cal`（`+0.2611`, `+0.1043`, `+0.3505`）和斜率（`+0.0588`, `+0.0200`, `+0.2503`）。此处为其补充与 **H-B3** 完全同构的事后特异性对照（5 个错配源线对照均值 `ctrl_mean_libnorm(s')` 代替 `mu_t` + 100 次基因标签置换 `mu_t[perm_k]`），报告 **3 条确认性线**（`nadig_hepg2`, `nadig_jurkat`, `vcc2025_h1`）的 `r_orth_cal`、OLS 斜率及 97.5% 双侧 CI（外加次要线附录）。

## 1. 确认性线（K_conf = 3）探索性 H-B3 类比判定汇总（97.5% CI）

| Target Line | N | FCT Matched `r_orth_cal` [97.5% CI] | Max Mismatch Mean `r_orth_cal` CI hi | Perm100 `r_orth_cal` p2.5 / med / p97.5 | FCT `r_orth_cal` Beats All 6 Nulls? | FCT Matched `slope_ols` [97.5% CI] | Max Mismatch Mean `slope` CI hi | Perm100 `slope` p97.5 | FCT `slope` Beats All 6 Nulls? |
|---|---:|---|---:|---|---|---|---:|---:|---|
| **nadig_hepg2** | 799 | **+0.2611** [+0.2420, +0.2807] | +0.1319 | +0.0884 / +0.1265 / +0.1696 | **True** | **+0.0588** [+0.0536, +0.0640] | +0.0434 | -0.0043 | **True** |
| **nadig_jurkat** | 934 | **+0.1043** [+0.0849, +0.1238] | +0.0847 | +0.0215 / +0.0680 / +0.1201 | **False** | **+0.0200** [+0.0149, +0.0253] | +0.0367 | +0.0480 | **False** |
| **vcc2025_h1** | 80 | **+0.3505** [+0.2784, +0.4245] | +0.2022 | +0.1725 / +0.2101 / +0.2445 | **True** | **+0.2503** [+0.1339, +0.3542] | +0.0520 | +0.2409 | **False** |

## 2. 全 6 条线 100 次基因置换分布（次要线含附录）

| Cohort | Target Line | N | Matched `r_orth_cal` [97.5% CI] | Max Mismatch `r_orth_cal` CI hi | Perm100 `r_orth_cal` [p2.5, med, p97.5, max] | Beats All 6 (`r_orth_cal`) | Matched `slope` [97.5% CI] | Max Mismatch `slope` CI hi | Perm100 `slope` p97.5 | Beats All 6 (`slope`) |
|---|---|---:|---|---:|---|---|---|---:|---:|---|
| Confirmatory | nadig_hepg2 | 799 | +0.2611 [+0.2420, +0.2807] | +0.1319 | [+0.0884, +0.1265, +0.1696, +0.1822] | True | +0.0588 [+0.0536, +0.0640] | +0.0434 | -0.0043 | True |
| Confirmatory | nadig_jurkat | 934 | +0.1043 [+0.0849, +0.1238] | +0.0847 | [+0.0215, +0.0680, +0.1201, +0.1294] | False | +0.0200 [+0.0149, +0.0253] | +0.0367 | +0.0480 | False |
| Confirmatory | vcc2025_h1 | 80 | +0.3505 [+0.2784, +0.4245] | +0.2022 | [+0.1725, +0.2101, +0.2445, +0.2671] | True | +0.2503 [+0.1339, +0.3542] | +0.0520 | +0.2409 | False |
| Secondary | nourreddine_kolf21j | 18 | +0.1839 [+0.0214, +0.3545] | +0.1752 | [-0.0017, +0.0577, +0.1246, +0.1411] | False | +0.0747 [+0.0070, +0.1592] | +0.1202 | +0.1873 | False |
| Secondary | replogle_k562_essential | 859 | +0.1938 [+0.1748, +0.2150] | +0.1258 | [+0.0582, +0.1157, +0.1563, +0.1620] | True | +0.0419 [+0.0355, +0.0489] | +0.1807 | +0.2674 | False |
| Secondary | replogle_rpe1 | 1070 | +0.1296 [+0.1115, +0.1475] | +0.1023 | [+0.0339, +0.0757, +0.1065, +0.1100] | True | +0.0228 [+0.0202, +0.0253] | +0.0144 | -0.0773 | True |

## 3. 逐控制项明细（匹配 FCT vs 全部 5 个错配源线对照均值 `ctrl_mean_libnorm`）

| Cohort | Target Line | Control | N | `r_orth_cal` | 97.5% CI | `slope_ols` | 97.5% CI | `r_cal` (raw) | 97.5% CI |
|---|---|---|---:|---:|---|---:|---|---:|---|
| Confirmatory | nadig_hepg2 | `FCT_Matched` | 799 | +0.2611 | [+0.2420, +0.2807] | +0.0588 | [+0.0536, +0.0640] | +0.2618 | [+0.2425, +0.2807] |
| Confirmatory | nadig_hepg2 | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 799 | +0.0453 | [+0.0320, +0.0592] | +0.0084 | [+0.0054, +0.0114] | +0.0494 | [+0.0325, +0.0664] |
| Confirmatory | nadig_hepg2 | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 799 | +0.0554 | [+0.0415, +0.0690] | +0.0024 | [-0.0034, +0.0080] | +0.0099 | [-0.0144, +0.0330] |
| Confirmatory | nadig_hepg2 | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 799 | +0.0665 | [+0.0532, +0.0803] | +0.0087 | [+0.0052, +0.0123] | +0.0518 | [+0.0319, +0.0712] |
| Confirmatory | nadig_hepg2 | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 799 | -0.0048 | [-0.0225, +0.0129] | +0.0041 | [+0.0018, +0.0065] | +0.0342 | [+0.0146, +0.0537] |
| Confirmatory | nadig_hepg2 | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 799 | +0.1186 | [+0.1058, +0.1319] | +0.0279 | [+0.0155, +0.0434] | +0.0649 | [+0.0375, +0.0959] |
| Confirmatory | nadig_jurkat | `FCT_Matched` | 934 | +0.1043 | [+0.0849, +0.1238] | +0.0200 | [+0.0149, +0.0253] | +0.0915 | [+0.0693, +0.1144] |
| Confirmatory | nadig_jurkat | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 934 | +0.0640 | [+0.0503, +0.0788] | +0.0140 | [+0.0106, +0.0179] | +0.0807 | [+0.0619, +0.0997] |
| Confirmatory | nadig_jurkat | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 934 | +0.0511 | [+0.0389, +0.0633] | +0.0168 | [+0.0120, +0.0221] | +0.0659 | [+0.0480, +0.0851] |
| Confirmatory | nadig_jurkat | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 934 | +0.0700 | [+0.0555, +0.0847] | +0.0154 | [+0.0118, +0.0192] | +0.0885 | [+0.0693, +0.1074] |
| Confirmatory | nadig_jurkat | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 934 | -0.0491 | [-0.0644, -0.0343] | -0.0102 | [-0.0129, -0.0076] | -0.0741 | [-0.0942, -0.0544] |
| Confirmatory | nadig_jurkat | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 934 | +0.0520 | [+0.0352, +0.0684] | +0.0145 | [-0.0046, +0.0367] | +0.0307 | [-0.0100, +0.0763] |
| Confirmatory | vcc2025_h1 | `FCT_Matched` | 80 | +0.3505 | [+0.2784, +0.4245] | +0.2503 | [+0.1339, +0.3542] | +0.5489 | [+0.3515, +0.6808] |
| Confirmatory | vcc2025_h1 | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 80 | +0.1240 | [+0.0815, +0.1664] | +0.0095 | [+0.0003, +0.0205] | +0.0603 | [+0.0022, +0.1230] |
| Confirmatory | vcc2025_h1 | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 80 | +0.0550 | [+0.0114, +0.0978] | +0.0076 | [+0.0021, +0.0136] | +0.0565 | [+0.0148, +0.0993] |
| Confirmatory | vcc2025_h1 | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 80 | +0.1453 | [+0.0890, +0.2022] | +0.0374 | [+0.0194, +0.0520] | +0.1905 | [+0.0913, +0.2815] |
| Confirmatory | vcc2025_h1 | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 80 | +0.0536 | [+0.0131, +0.0916] | +0.0013 | [-0.0065, +0.0108] | +0.0104 | [-0.0575, +0.0734] |
| Confirmatory | vcc2025_h1 | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 80 | -0.0541 | [-0.1091, +0.0042] | -0.0086 | [-0.0194, -0.0001] | -0.0700 | [-0.1454, -0.0007] |
| Secondary | nourreddine_kolf21j | `FCT_Matched` | 18 | +0.1839 | [+0.0214, +0.3545] | +0.0747 | [+0.0070, +0.1592] | +0.2547 | [+0.0269, +0.4664] |
| Secondary | nourreddine_kolf21j | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 18 | +0.1081 | [+0.0347, +0.1678] | +0.0656 | [+0.0167, +0.1202] | +0.2460 | [+0.0877, +0.3851] |
| Secondary | nourreddine_kolf21j | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 18 | +0.0633 | [-0.0521, +0.1447] | +0.0359 | [-0.0099, +0.0889] | +0.1300 | [-0.0508, +0.2390] |
| Secondary | nourreddine_kolf21j | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 18 | +0.0870 | [+0.0088, +0.1439] | +0.0443 | [+0.0054, +0.1075] | +0.1609 | [+0.0433, +0.2526] |
| Secondary | nourreddine_kolf21j | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 18 | +0.0528 | [-0.0564, +0.1183] | +0.0605 | [+0.0060, +0.1158] | +0.2007 | [+0.0354, +0.2998] |
| Secondary | nourreddine_kolf21j | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 18 | +0.0708 | [-0.0223, +0.1752] | +0.0300 | [+0.0028, +0.0605] | +0.0889 | [+0.0087, +0.1699] |
| Secondary | replogle_k562_essential | `FCT_Matched` | 859 | +0.1938 | [+0.1748, +0.2150] | +0.0419 | [+0.0355, +0.0489] | +0.1468 | [+0.1236, +0.1708] |
| Secondary | replogle_k562_essential | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 859 | +0.0479 | [+0.0292, +0.0728] | +0.0083 | [+0.0032, +0.0145] | +0.0399 | [+0.0150, +0.0708] |
| Secondary | replogle_k562_essential | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 859 | +0.0493 | [+0.0292, +0.0738] | +0.0079 | [+0.0012, +0.0153] | +0.0362 | [+0.0056, +0.0680] |
| Secondary | replogle_k562_essential | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 859 | -0.0097 | [-0.0264, +0.0091] | +0.0248 | [+0.0151, +0.0345] | +0.0738 | [+0.0456, +0.1016] |
| Secondary | replogle_k562_essential | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 859 | -0.0506 | [-0.0766, -0.0143] | -0.0088 | [-0.0135, -0.0029] | -0.0524 | [-0.0828, -0.0168] |
| Secondary | replogle_k562_essential | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 859 | +0.1090 | [+0.0928, +0.1258] | +0.1511 | [+0.1238, +0.1807] | +0.2360 | [+0.2010, +0.2696] |
| Secondary | replogle_rpe1 | `FCT_Matched` | 1070 | +0.1296 | [+0.1115, +0.1475] | +0.0228 | [+0.0202, +0.0253] | +0.1592 | [+0.1397, +0.1783] |
| Secondary | replogle_rpe1 | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 1070 | +0.0883 | [+0.0740, +0.1023] | +0.0128 | [+0.0111, +0.0144] | +0.1265 | [+0.1095, +0.1424] |
| Secondary | replogle_rpe1 | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 1070 | +0.0056 | [-0.0098, +0.0208] | +0.0044 | [+0.0024, +0.0064] | +0.0404 | [+0.0219, +0.0589] |
| Secondary | replogle_rpe1 | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 1070 | +0.0096 | [-0.0027, +0.0224] | -0.0127 | [-0.0155, -0.0100] | -0.0893 | [-0.1084, -0.0702] |
| Secondary | replogle_rpe1 | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 1070 | -0.0050 | [-0.0170, +0.0064] | -0.0031 | [-0.0047, -0.0014] | -0.0316 | [-0.0491, -0.0149] |
| Secondary | replogle_rpe1 | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 1070 | +0.0613 | [+0.0490, +0.0732] | -0.0265 | [-0.0325, -0.0208] | -0.0983 | [-0.1198, -0.0774] |
