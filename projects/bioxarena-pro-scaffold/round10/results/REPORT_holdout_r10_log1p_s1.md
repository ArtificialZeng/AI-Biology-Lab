# Round 10 Confirmatory Group C (`FCT`) Report (`holdout_r10_log1p_s1`)

- **Phase**: `holdout_r10` | **PREREG_C Commit**: `PREREG_C_COMMIT（主控核验通过）` | **Code Freeze**: `f8d2ebfc238d83a0c3a6c86696950ed292bf90e9`
- **Bootstrap**: `np.random.default_rng(20261122)` x `4000` (Bonferroni two-sided `97.5%` CIs)
- **Cohort Disclosure**: **K_conf = 2 confirmatory cell lines** (`['dld1', 'cd4t_rest']`)
- **Secondary Lines**: `['calu3', 'wtc11_i3neuron']`
- **Primary Group C (`FCT`) Verdict**: **未确认**
- **Parallel Practical Utility (`H-C3` Non-Inferiority)**: **未确认**

## Confirmatory Hypothesis Adjudication (97.5% Bonferroni CIs)

| Hypothesis | Pass / Total | Criterion | Verdict |
|---|---:|---|---|
| **H-C1** (`FCT` `r_orth_cal` 97.5% CI lo > 0) | 2/2 | Unanimous K_conf | **CONFIRMED** |
| **H-C2** (`FCT` `r_orth_cal` lo > max 6 mismatched hi & perm100 p97.5) | 0/2 | Unanimous K_conf | **未确认** |
| **H-C3** (Non-inferiority vs `Transfer_Mean` across 6 conventions) | 8/12 non-inf (`1` strict sup) | 0 inferior pairs | **未确认** |

## Per-Line Results (`resp_top50`, Joint-Bootstrap 97.5% CIs)

| Target Line | Role | Method | N | r_orth_cal | 97.5% CI | Slope OLS | 97.5% CI | r_cal (raw) | rho_AB_orth | pds_cos |
|---|---|---|---:|---:|---|---:|---|---:|---:|---:|
| calu3 | Secondary | `Source_Grand_Mean` | 53 | -0.0291 | [-0.1509, +0.0942] | +0.0373 | [-0.0022, +0.0801] | +0.0969 | +0.7505 | 0.5000 |
| calu3 | Secondary | `Transfer_Mean` | 53 | +nan | [+nan, +nan] | +0.0000 | [+0.0000, +0.0000] | +nan | +0.7505 | 0.6687 |
| calu3 | Secondary | `FCT` | 53 | +0.0523 | [-0.0185, +0.1255] | +0.0043 | [-0.0062, +0.0161] | +0.0270 | +0.7505 | 0.6676 |
| calu3 | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 53 | +0.0031 | [-0.0421, +0.0464] | -0.0003 | [-0.0052, +0.0043] | -0.0037 | +0.7505 | 0.6698 |
| calu3 | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 53 | -0.0723 | [-0.1196, -0.0283] | -0.0052 | [-0.0095, -0.0006] | -0.0550 | +0.7505 | 0.6560 |
| calu3 | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 53 | -0.0162 | [-0.0663, +0.0367] | +0.0017 | [-0.0036, +0.0075] | +0.0160 | +0.7505 | 0.6738 |
| calu3 | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 53 | -0.0347 | [-0.0861, +0.0190] | -0.0011 | [-0.0063, +0.0045] | -0.0115 | +0.7505 | 0.6560 |
| calu3 | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 53 | +0.0079 | [-0.0432, +0.0628] | +0.0003 | [-0.0048, +0.0063] | +0.0032 | +0.7505 | 0.6571 |
| calu3 | Secondary | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 53 | +0.0411 | [-0.0349, +0.1103] | +0.0110 | [-0.0024, +0.0254] | +0.0605 | +0.7505 | 0.6687 |
| calu3 | Secondary | `NegCtrl_FCT_Perm0` | 53 | +0.0408 | [-0.0272, +0.1122] | +0.0343 | [+0.0041, +0.0689] | +0.1252 | +0.7505 | 0.6821 |
| calu3 | Secondary | `FCT_FullNTMean` | 53 | +0.0523 | [-0.0185, +0.1255] | +0.0043 | [-0.0062, +0.0161] | +0.0270 | +0.7505 | 0.6676 |
| calu3 | Secondary | `Transfer_Mean_GlobalScale` | 53 | +nan | [+nan, +nan] | -0.0083 | [-0.0165, -0.0012] | -0.1264 | +0.7505 | 0.6687 |
| calu3 | Secondary | `Transfer_Mean_VarScaled` | 53 | +0.0098 | [-0.0558, +0.0759] | -0.0381 | [-0.0918, +0.0075] | -0.0996 | +0.7505 | 0.6676 |
| calu3 | Secondary | `CIPHER_Linear` | 53 | -0.0089 | [-0.2574, +0.2268] | +0.0403 | [-0.0679, +0.1393] | +0.0709 | +0.7505 | 0.5269 |
| calu3 | Secondary | `Hybrid_ProjShared` | 53 | -0.0287 | [-0.2516, +0.1790] | -0.0135 | [-0.0979, +0.0641] | -0.0354 | +0.7505 | 0.5980 |
| calu3 | Secondary | `Hybrid_on_TMVS` | 53 | -0.0195 | [-0.2252, +0.1644] | -0.0516 | [-0.1623, +0.0522] | -0.0943 | +0.7505 | 0.6089 |
| cd4t_rest | Confirmatory | `Source_Grand_Mean` | 1422 | +0.0229 | [-0.0113, +0.0559] | +0.0407 | [+0.0289, +0.0560] | +0.1623 | +0.6350 | 0.5000 |
| cd4t_rest | Confirmatory | `Transfer_Mean` | 1422 | +nan | [+nan, +nan] | +0.0000 | [+0.0000, +0.0000] | +nan | +0.6350 | 0.5880 |
| cd4t_rest | Confirmatory | `FCT` | 1422 | +0.0314 | [+0.0124, +0.0502] | +0.0058 | [+0.0002, +0.0133] | +0.0405 | +0.6350 | 0.5864 |
| cd4t_rest | Confirmatory | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 1422 | -0.0007 | [-0.0165, +0.0149] | +0.0017 | [+0.0006, +0.0029] | +0.0284 | +0.6350 | 0.5878 |
| cd4t_rest | Confirmatory | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 1422 | -0.0024 | [-0.0217, +0.0155] | -0.0039 | [-0.0062, -0.0017] | -0.0534 | +0.6350 | 0.5798 |
| cd4t_rest | Confirmatory | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 1422 | -0.0089 | [-0.0220, +0.0045] | +0.0040 | [+0.0023, +0.0061] | +0.0491 | +0.6350 | 0.5816 |
| cd4t_rest | Confirmatory | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 1422 | -0.0084 | [-0.0259, +0.0088] | -0.0006 | [-0.0018, +0.0008] | -0.0086 | +0.6350 | 0.5808 |
| cd4t_rest | Confirmatory | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 1422 | -0.0036 | [-0.0218, +0.0141] | -0.0013 | [-0.0029, +0.0001] | -0.0236 | +0.6350 | 0.5841 |
| cd4t_rest | Confirmatory | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 1422 | +0.0034 | [-0.0171, +0.0232] | +0.0087 | [+0.0008, +0.0202] | +0.0567 | +0.6350 | 0.5805 |
| cd4t_rest | Confirmatory | `NegCtrl_FCT_Perm0` | 1422 | +0.0215 | [+0.0076, +0.0353] | +0.0301 | [+0.0216, +0.0415] | +0.1636 | +0.6350 | 0.5909 |
| cd4t_rest | Confirmatory | `FCT_FullNTMean` | 1422 | +0.0314 | [+0.0124, +0.0502] | +0.0058 | [+0.0002, +0.0133] | +0.0405 | +0.6350 | 0.5864 |
| cd4t_rest | Confirmatory | `Transfer_Mean_GlobalScale` | 1422 | +nan | [+nan, +nan] | -0.0231 | [-0.0300, -0.0173] | -0.1995 | +0.6350 | 0.5880 |
| cd4t_rest | Confirmatory | `Transfer_Mean_VarScaled` | 1422 | +0.0164 | [+0.0003, +0.0328] | -0.0451 | [-0.0630, -0.0312] | -0.1777 | +0.6350 | 0.5880 |
| cd4t_rest | Confirmatory | `CIPHER_Linear` | 1422 | -0.0058 | [-0.0452, +0.0348] | +0.0451 | [+0.0291, +0.0635] | +0.1500 | +0.6350 | 0.4986 |
| cd4t_rest | Confirmatory | `Hybrid_ProjShared` | 1422 | -0.0095 | [-0.0391, +0.0208] | -0.0030 | [-0.0071, +0.0009] | -0.0278 | +0.6350 | 0.5657 |
| cd4t_rest | Confirmatory | `Hybrid_on_TMVS` | 1422 | +0.0049 | [-0.0178, +0.0285] | -0.0481 | [-0.0663, -0.0335] | -0.1684 | +0.6350 | 0.5723 |
| dld1 | Confirmatory | `Source_Grand_Mean` | 500 | +0.1427 | [+0.1057, +0.1810] | -0.0139 | [-0.0308, +0.0031] | -0.0350 | +0.7639 | 0.5000 |
| dld1 | Confirmatory | `Transfer_Mean` | 500 | +nan | [+nan, +nan] | +0.0000 | [+0.0000, +0.0000] | +nan | +0.7639 | 0.7049 |
| dld1 | Confirmatory | `FCT` | 500 | +0.0717 | [+0.0462, +0.0992] | +0.0107 | [+0.0077, +0.0139] | +0.0997 | +0.7639 | 0.7001 |
| dld1 | Confirmatory | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 500 | -0.0005 | [-0.0190, +0.0187] | +0.0002 | [-0.0013, +0.0017] | +0.0025 | +0.7639 | 0.7057 |
| dld1 | Confirmatory | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 500 | -0.0267 | [-0.0454, -0.0077] | -0.0019 | [-0.0036, -0.0002] | -0.0212 | +0.7639 | 0.6975 |
| dld1 | Confirmatory | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 500 | +0.0055 | [-0.0108, +0.0219] | -0.0018 | [-0.0040, +0.0005] | -0.0152 | +0.7639 | 0.7001 |
| dld1 | Confirmatory | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 500 | -0.0106 | [-0.0295, +0.0095] | -0.0023 | [-0.0043, -0.0005] | -0.0269 | +0.7639 | 0.6989 |
| dld1 | Confirmatory | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 500 | +0.0067 | [-0.0138, +0.0270] | +0.0014 | [-0.0002, +0.0031] | +0.0180 | +0.7639 | 0.7007 |
| dld1 | Confirmatory | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 500 | +0.0228 | [+0.0057, +0.0406] | -0.0036 | [-0.0105, +0.0033] | -0.0158 | +0.7639 | 0.6963 |
| dld1 | Confirmatory | `NegCtrl_FCT_Perm0` | 500 | +0.0574 | [+0.0387, +0.0764] | -0.0160 | [-0.0269, -0.0049] | -0.0527 | +0.7639 | 0.7054 |
| dld1 | Confirmatory | `FCT_FullNTMean` | 500 | +0.0717 | [+0.0462, +0.0992] | +0.0107 | [+0.0077, +0.0139] | +0.0997 | +0.7639 | 0.7001 |
| dld1 | Confirmatory | `Transfer_Mean_GlobalScale` | 500 | +nan | [+nan, +nan] | +0.0112 | [+0.0065, +0.0157] | +0.1026 | +0.7639 | 0.7049 |
| dld1 | Confirmatory | `Transfer_Mean_VarScaled` | 500 | +0.0456 | [+0.0117, +0.0791] | +0.0168 | [+0.0084, +0.0247] | +0.0840 | +0.7639 | 0.7072 |
| dld1 | Confirmatory | `CIPHER_Linear` | 500 | +0.0019 | [-0.0366, +0.0406] | -0.0458 | [-0.0653, -0.0267] | -0.0996 | +0.7639 | 0.4777 |
| dld1 | Confirmatory | `Hybrid_ProjShared` | 500 | +0.0193 | [-0.0258, +0.0603] | +0.0014 | [-0.0016, +0.0040] | +0.0196 | +0.7639 | 0.6958 |
| dld1 | Confirmatory | `Hybrid_on_TMVS` | 500 | +0.0478 | [+0.0099, +0.0839] | +0.0182 | [+0.0094, +0.0267] | +0.0850 | +0.7639 | 0.7007 |
| wtc11_i3neuron | Secondary | `Source_Grand_Mean` | 67 | +0.0422 | [-0.0156, +0.0981] | +0.2041 | [+0.1121, +0.3404] | +0.4594 | +0.7387 | 0.5000 |
| wtc11_i3neuron | Secondary | `Transfer_Mean` | 67 | +nan | [+nan, +nan] | +0.0000 | [+0.0000, +0.0000] | +nan | +0.7387 | 0.5407 |
| wtc11_i3neuron | Secondary | `FCT` | 67 | +0.0402 | [-0.0129, +0.0959] | -0.0119 | [-0.0228, -0.0012] | -0.0647 | +0.7387 | 0.5529 |
| wtc11_i3neuron | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 67 | +0.0081 | [-0.0310, +0.0565] | -0.0004 | [-0.0115, +0.0081] | -0.0040 | +0.7387 | 0.5382 |
| wtc11_i3neuron | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 67 | -0.0095 | [-0.0463, +0.0334] | -0.0088 | [-0.0147, -0.0023] | -0.0782 | +0.7387 | 0.5328 |
| wtc11_i3neuron | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 67 | -0.0109 | [-0.0510, +0.0299] | +0.0093 | [+0.0003, +0.0196] | +0.0707 | +0.7387 | 0.5423 |
| wtc11_i3neuron | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 67 | +0.0025 | [-0.0485, +0.0532] | -0.0021 | [-0.0129, +0.0072] | -0.0231 | +0.7387 | 0.5355 |
| wtc11_i3neuron | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 67 | +0.0148 | [-0.0349, +0.0638] | +0.0009 | [-0.0057, +0.0076] | +0.0089 | +0.7387 | 0.5319 |
| wtc11_i3neuron | Secondary | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 67 | +0.0116 | [-0.0336, +0.0568] | +0.0720 | [+0.0062, +0.1984] | +0.2431 | +0.7387 | 0.5353 |
| wtc11_i3neuron | Secondary | `NegCtrl_FCT_Perm0` | 67 | -0.0134 | [-0.0552, +0.0271] | +0.1261 | [+0.0546, +0.2454] | +0.3607 | +0.7387 | 0.5455 |
| wtc11_i3neuron | Secondary | `FCT_FullNTMean` | 67 | +0.0402 | [-0.0129, +0.0959] | -0.0119 | [-0.0228, -0.0012] | -0.0647 | +0.7387 | 0.5529 |
| wtc11_i3neuron | Secondary | `Transfer_Mean_GlobalScale` | 67 | +nan | [+nan, +nan] | -0.0066 | [-0.0103, -0.0039] | -0.4907 | +0.7387 | 0.5407 |
| wtc11_i3neuron | Secondary | `Transfer_Mean_VarScaled` | 67 | +0.0218 | [-0.0144, +0.0627] | -0.0225 | [-0.0303, -0.0154] | -0.2326 | +0.7387 | 0.5486 |
| wtc11_i3neuron | Secondary | `CIPHER_Linear` | 67 | +0.0016 | [-0.0866, +0.0819] | +0.2380 | [+0.1383, +0.3756] | +0.4374 | +0.7387 | 0.5084 |
| wtc11_i3neuron | Secondary | `Hybrid_ProjShared` | 67 | -0.0121 | [-0.0830, +0.0520] | +0.0003 | [-0.0138, +0.0133] | +0.0012 | +0.7387 | 0.5678 |
| wtc11_i3neuron | Secondary | `Hybrid_on_TMVS` | 67 | -0.0048 | [-0.0718, +0.0586] | -0.0222 | [-0.0393, -0.0066] | -0.0858 | +0.7387 | 0.5719 |

## 100-Perm `NegCtrl_FCT_Perm` (`k=0..99`) & Mismatched Mean Summary

- **calu3**: perm100 `r_orth_cal` median=+0.0301, p97.5=+0.0670, p98.75=+0.0712; max mismatched hi=+0.1103
- **cd4t_rest**: perm100 `r_orth_cal` median=+0.0065, p97.5=+0.0183, p98.75=+0.0210; max mismatched hi=+0.0232
- **dld1**: perm100 `r_orth_cal` median=+0.0369, p97.5=+0.0570, p98.75=+0.0584; max mismatched hi=+0.0406
- **wtc11_i3neuron**: perm100 `r_orth_cal` median=+0.0075, p97.5=+0.0422, p98.75=+0.0467; max mismatched hi=+0.0638

## `SA-3` Round 9 Anchors Notice (`P7`) & `SA-6 / S7` Cell-Count Controls (`P6`)

- *Benchmark Anchor Note (`P7`):* `Hybrid_ProjShared` follows frozen Round 9 commit `1ca9f60f` 3D `span{1, b_S, a_S}` projection (whereas Round 9 `PREREG_A v4` prose described rank-1 `q_1`).
- **calu3 S7**: Spearman(`r_orth_p`, ln n_cells)=-0.1363 [-0.4361, +0.1745], Spearman(`r_orth_p`, |k_p|)=-0.0123 [-0.2976, +0.2811], OLS-adjusted intercept=-0.0087 [-0.0822, +0.0650]
  * `T1_low` (N=18, n_cells=[17..54]): FCT `r_orth_cal`=+0.0743 [-0.0525, +0.1987], max_mis_hi=+0.1701, perm100_p975=+0.0754, HC1=False, HC2=False
  * `T2_mid` (N=18, n_cells=[58..117]): FCT `r_orth_cal`=+0.0053 [-0.1192, +0.1326], max_mis_hi=+0.1134, perm100_p975=+0.0388, HC1=False, HC2=False
  * `T3_high` (N=17, n_cells=[122..280]): FCT `r_orth_cal`=+0.0277 [-0.1107, +0.1512], max_mis_hi=+0.2376, perm100_p975=+0.1082, HC1=False, HC2=False
- **cd4t_rest S7**: Spearman(`r_orth_p`, ln n_cells)=-0.0100 [-0.0673, +0.0526], Spearman(`r_orth_p`, |k_p|)=+0.0840 [+0.0245, +0.1448], OLS-adjusted intercept=+0.0257 [+0.0130, +0.0382]
  * `T1_low` (N=474, n_cells=[10..41]): FCT `r_orth_cal`=+0.0214 [-0.0112, +0.0523], max_mis_hi=+0.0242, perm100_p975=+0.0085, HC1=False, HC2=False
  * `T2_mid` (N=474, n_cells=[41..113]): FCT `r_orth_cal`=+0.0497 [+0.0222, +0.0776], max_mis_hi=+0.0525, perm100_p975=+0.0330, HC1=True, HC2=False
  * `T3_high` (N=474, n_cells=[113..1786]): FCT `r_orth_cal`=+0.0376 [+0.0077, +0.0669], max_mis_hi=+0.0638, perm100_p975=+0.0333, HC1=True, HC2=False
- **dld1 S7**: Spearman(`r_orth_p`, ln n_cells)=-0.0571 [-0.1641, +0.0475], Spearman(`r_orth_p`, |k_p|)=+0.0924 [-0.0098, +0.1931], OLS-adjusted intercept=+0.0508 [+0.0333, +0.0685]
  * `T1_low` (N=167, n_cells=[10..30]): FCT `r_orth_cal`=+0.1031 [+0.0640, +0.1475], max_mis_hi=+0.0557, perm100_p975=+0.0681, HC1=True, HC2=False
  * `T2_mid` (N=167, n_cells=[30..46]): FCT `r_orth_cal`=+0.0550 [+0.0185, +0.0909], max_mis_hi=+0.0544, perm100_p975=+0.0707, HC1=True, HC2=False
  * `T3_high` (N=166, n_cells=[46..122]): FCT `r_orth_cal`=+0.0395 [-0.0056, +0.0883], max_mis_hi=+0.0546, perm100_p975=+0.0384, HC1=False, HC2=False
- **wtc11_i3neuron S7**: Spearman(`r_orth_p`, ln n_cells)=+0.0639 [-0.2038, +0.3303], Spearman(`r_orth_p`, |k_p|)=-0.0652 [-0.3423, +0.2135], OLS-adjusted intercept=+0.0799 [+0.0404, +0.1189]
  * `T1_low` (N=23, n_cells=[36..75]): FCT `r_orth_cal`=+0.0463 [-0.0352, +0.1296], max_mis_hi=+0.0898, perm100_p975=+0.0636, HC1=False, HC2=False
  * `T2_mid` (N=22, n_cells=[75..96]): FCT `r_orth_cal`=+0.0668 [-0.0111, +0.1495], max_mis_hi=+0.1184, perm100_p975=+0.0664, HC1=False, HC2=False
  * `T3_high` (N=22, n_cells=[97..196]): FCT `r_orth_cal`=+0.0064 [-0.0727, +0.1070], max_mis_hi=+0.1347, perm100_p975=+0.0310, HC1=False, HC2=False
