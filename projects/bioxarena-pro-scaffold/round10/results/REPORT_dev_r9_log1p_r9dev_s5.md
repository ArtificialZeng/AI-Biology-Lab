# Round 10 Confirmatory Group C (`FCT`) Report (`dev_r9_log1p_r9dev_s5`)

- **Phase**: `dev_r9` | **PREREG_C Commit**: `` | **Code Freeze**: ``
- **Bootstrap**: `np.random.default_rng(20261122)` x `500` (Bonferroni two-sided `97.5%` CIs)
- **Cohort Disclosure**: **Insufficient confirmatory cell lines (K_conf < 2); confirmatory hypothesis testing was not performed (确认性线不足，未做确认性检验).** (`[]`)
- **Secondary Lines**: `['nadig_hepg2', 'vcc2025_h1']`
- **Primary Group C (`FCT`) Verdict**: **NOT_EVALUATED_INSUFFICIENT_LINES**
- **Parallel Practical Utility (`H-C3` Non-Inferiority)**: **NOT_EVALUATED_INSUFFICIENT_LINES**

## Per-Line Results (`resp_top50`, Joint-Bootstrap 97.5% CIs)

| Target Line | Role | Method | N | r_orth_cal | 97.5% CI | Slope OLS | 97.5% CI | r_cal (raw) | rho_AB_orth | pds_cos |
|---|---|---|---:|---:|---|---:|---|---:|---:|---:|
| nadig_hepg2 | Secondary | `Source_Grand_Mean` | 797 | +0.3119 | [+0.2769, +0.3471] | +0.0211 | [-0.0090, +0.0568] | +0.0315 | +0.8741 | 0.5000 |
| nadig_hepg2 | Secondary | `Transfer_Mean` | 797 | +nan | [+nan, +nan] | +0.0000 | [+0.0000, +0.0000] | +nan | +0.8741 | 0.8790 |
| nadig_hepg2 | Secondary | `FCT` | 797 | +0.2627 | [+0.2451, +0.2839] | +0.0590 | [+0.0540, +0.0645] | +0.2639 | +0.8741 | 0.8852 |
| nadig_hepg2 | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 797 | +0.0422 | [+0.0285, +0.0538] | +0.0077 | [+0.0048, +0.0105] | +0.0458 | +0.8741 | 0.8719 |
| nadig_hepg2 | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 797 | +0.0520 | [+0.0382, +0.0664] | +0.0016 | [-0.0038, +0.0085] | +0.0066 | +0.8741 | 0.8750 |
| nadig_hepg2 | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 797 | +0.0651 | [+0.0516, +0.0777] | +0.0082 | [+0.0047, +0.0113] | +0.0502 | +0.8741 | 0.8740 |
| nadig_hepg2 | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 797 | -0.0109 | [-0.0294, +0.0074] | +0.0034 | [+0.0012, +0.0059] | +0.0285 | +0.8741 | 0.8767 |
| nadig_hepg2 | Secondary | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 797 | +0.1193 | [+0.1062, +0.1315] | +0.0285 | [+0.0155, +0.0448] | +0.0663 | +0.8741 | 0.8776 |
| nadig_hepg2 | Secondary | `NegCtrl_FCT_Perm0` | 797 | +0.1602 | [+0.1450, +0.1733] | -0.0144 | [-0.0370, +0.0129] | -0.0243 | +0.8741 | 0.8930 |
| nadig_hepg2 | Secondary | `Transfer_Mean_GlobalScale` | 797 | +nan | [+nan, +nan] | -0.0083 | [-0.0111, -0.0052] | -0.1288 | +0.8741 | 0.8790 |
| nadig_hepg2 | Secondary | `Transfer_Mean_VarScaled` | 797 | +0.0424 | [+0.0237, +0.0609] | +0.0357 | [+0.0276, +0.0442] | +0.1226 | +0.8741 | 0.8736 |
| nadig_hepg2 | Secondary | `CIPHER_Linear` | 797 | +0.1169 | [+0.0856, +0.1498] | -0.0657 | [-0.0995, -0.0252] | -0.0814 | +0.8741 | 0.4840 |
| nadig_hepg2 | Secondary | `Hybrid_ProjShared` | 797 | +0.0483 | [+0.0232, +0.0757] | +0.0098 | [+0.0064, +0.0137] | +0.0723 | +0.8741 | 0.8675 |
| nadig_hepg2 | Secondary | `Hybrid_on_TMVS` | 797 | +0.0551 | [+0.0342, +0.0782] | +0.0456 | [+0.0351, +0.0567] | +0.1291 | +0.8741 | 0.8614 |
| vcc2025_h1 | Secondary | `Source_Grand_Mean` | 80 | +0.1063 | [+0.0395, +0.1649] | +0.1896 | [-0.0072, +0.3800] | +0.3021 | +0.9885 | 0.5000 |
| vcc2025_h1 | Secondary | `Transfer_Mean` | 80 | +nan | [+nan, +nan] | +0.0000 | [+0.0000, +0.0000] | +nan | +0.9885 | 0.8758 |
| vcc2025_h1 | Secondary | `FCT` | 80 | +0.3495 | [+0.2799, +0.4170] | +0.2509 | [+0.1435, +0.3545] | +0.5497 | +0.9885 | 0.8998 |
| vcc2025_h1 | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 80 | +0.1265 | [+0.0838, +0.1721] | +0.0094 | [+0.0008, +0.0190] | +0.0625 | +0.9885 | 0.8848 |
| vcc2025_h1 | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 80 | +0.0616 | [+0.0218, +0.0962] | +0.0077 | [+0.0024, +0.0126] | +0.0589 | +0.9885 | 0.8691 |
| vcc2025_h1 | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 80 | +0.1523 | [+0.1016, +0.2088] | +0.0385 | [+0.0234, +0.0529] | +0.1979 | +0.9885 | 0.8826 |
| vcc2025_h1 | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 80 | +0.0475 | [+0.0079, +0.0881] | +0.0008 | [-0.0068, +0.0099] | +0.0063 | +0.9885 | 0.8687 |
| vcc2025_h1 | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 80 | -0.0469 | [-0.1052, +0.0020] | -0.0079 | [-0.0173, +0.0002] | -0.0669 | +0.9885 | 0.8786 |
| vcc2025_h1 | Secondary | `NegCtrl_FCT_Perm0` | 80 | +0.1813 | [+0.1378, +0.2275] | +0.1780 | [+0.0545, +0.2933] | +0.4077 | +0.9885 | 0.8965 |
| vcc2025_h1 | Secondary | `Transfer_Mean_GlobalScale` | 80 | +nan | [+nan, +nan] | +0.0688 | [-0.0165, +0.1533] | +0.2556 | +0.9885 | 0.8758 |
| vcc2025_h1 | Secondary | `Transfer_Mean_VarScaled` | 80 | -0.0791 | [-0.1917, +0.0214] | -0.1641 | [-0.2575, -0.0601] | -0.4487 | +0.9885 | 0.8750 |
| vcc2025_h1 | Secondary | `CIPHER_Linear` | 80 | -0.0451 | [-0.1381, +0.0394] | +0.1596 | [-0.0531, +0.3629] | +0.2353 | +0.9885 | 0.4663 |
| vcc2025_h1 | Secondary | `Hybrid_ProjShared` | 80 | -0.0394 | [-0.1273, +0.0383] | -0.0095 | [-0.0249, +0.0035] | -0.0532 | +0.9885 | 0.8410 |
| vcc2025_h1 | Secondary | `Hybrid_on_TMVS` | 80 | -0.0825 | [-0.1780, +0.0133] | -0.1736 | [-0.2725, -0.0715] | -0.4200 | +0.9885 | 0.8146 |

## 100-Perm `NegCtrl_FCT_Perm` (`k=0..99`) & Mismatched Mean Summary

- **nadig_hepg2**: perm100 `r_orth_cal` median=+0.1265, p97.5=+0.1692, p98.75=+0.1738; max mismatched hi=+0.1315
- **vcc2025_h1**: perm100 `r_orth_cal` median=+0.2108, p97.5=+0.2466, p98.75=+0.2555; max mismatched hi=+0.2088

## `SA-3` Round 9 Anchors Notice (`P7`) & `SA-6 / S7` Cell-Count Controls (`P6`)

- *Benchmark Anchor Note (`P7`):* `Hybrid_ProjShared` follows frozen Round 9 commit `1ca9f60f` 3D `span{1, b_S, a_S}` projection (whereas Round 9 `PREREG_A v4` prose described rank-1 `q_1`).
- **nadig_hepg2 S7**: Spearman(`r_orth_p`, ln n_cells)=+0.2084 [+0.1221, +0.2869], Spearman(`r_orth_p`, |k_p|)=-0.0154 [-0.0979, +0.0615], OLS-adjusted intercept=+0.2943 [+0.2763, +0.3119]
  * `T1_low` (N=266, n_cells=[25..48]): FCT `r_orth_cal`=+0.2117 [+0.1788, +0.2439], max_mis_hi=+0.1010, perm100_p975=+0.1474, HC1=True, HC2=True
  * `T2_mid` (N=266, n_cells=[48..75]): FCT `r_orth_cal`=+0.2589 [+0.2278, +0.2917], max_mis_hi=+0.1456, perm100_p975=+0.1847, HC1=True, HC2=True
  * `T3_high` (N=265, n_cells=[75..1213]): FCT `r_orth_cal`=+0.3333 [+0.3075, +0.3646], max_mis_hi=+0.1840, perm100_p975=+0.2027, HC1=True, HC2=True
- **vcc2025_h1 S7**: Spearman(`r_orth_p`, ln n_cells)=-0.0843 [-0.2967, +0.1786], Spearman(`r_orth_p`, |k_p|)=+0.0414 [-0.1823, +0.2529], OLS-adjusted intercept=+0.3196 [+0.2639, +0.3780]
  * `T1_low` (N=27, n_cells=[33..407]): FCT `r_orth_cal`=+0.3537 [+0.2043, +0.4879], max_mis_hi=+0.1935, perm100_p975=+0.2849, HC1=True, HC2=False
  * `T2_mid` (N=27, n_cells=[450..959]): FCT `r_orth_cal`=+0.3867 [+0.2742, +0.4873], max_mis_hi=+0.3166, perm100_p975=+0.2449, HC1=True, HC2=False
  * `T3_high` (N=26, n_cells=[972..2428]): FCT `r_orth_cal`=+0.3040 [+0.2157, +0.4035], max_mis_hi=+0.2100, perm100_p975=+0.2407, HC1=True, HC2=False
