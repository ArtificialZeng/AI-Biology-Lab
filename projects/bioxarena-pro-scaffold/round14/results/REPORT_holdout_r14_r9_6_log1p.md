# Round 14 Prospective Same-Line Independent Replication Report (`PREREG_G` — `holdout_r14_r9_6_log1p`)

- **Phase**: `holdout_r14` | **Source Mode**: `r9_6` (`6` sources) | **PREREG_E Commit**: `d34560d7ae568b299dca15c766c2e356b483c39f` | **Code Freeze**: `1409333c29103ed1b0ddfaf01c641e756eb5a109`
- **Bootstrap**: `np.random.default_rng(20261224)` x `4000` (Bonferroni 2-claim + IUT two-sided `97.5%` CIs `[1.25%, 98.75%]`)
- **Cohort Disclosure**: **K_conf = 0 (count gate failed); descriptive reporting only.**
- **Secondary / Non-Confirmatory Lines (Reporting-Only)**: `['a549_myllia']`
- **Confirmatory Claim 1 (`H-E1` Target-Line Specificity)**: **NOT_EVALUATED_NO_CONFIRMATORY_LINES**
- **Confirmatory Claim 2 (`H-E2` Persistent `SA-8` Constant-Baseline Inflation)**: **NOT_EVALUATED_NO_CONFIRMATORY_LINES**

## `H-E1b`: Paired Specificity Contrasts `FCT − NegCtrl_FCT_MismatchedMean(s)` (`all_nontarget`, 97.5% CI)

| Target Line | Role | Mismatched Source | Estimate | 97.5% CI | `lo > 0` |
|---|---|---|---:|---|---|
| a549_myllia | Secondary | `nadig_hepg2` | +0.0657 | [+0.0390, +0.0923] | **PASS** |
| a549_myllia | Secondary | `nadig_jurkat` | +0.0739 | [+0.0456, +0.1014] | **PASS** |
| a549_myllia | Secondary | `nourreddine_kolf21j` | +0.0984 | [+0.0631, +0.1333] | **PASS** |
| a549_myllia | Secondary | `replogle_k562_essential` | +0.0766 | [+0.0436, +0.1091] | **PASS** |
| a549_myllia | Secondary | `replogle_rpe1` | +0.0698 | [+0.0429, +0.0966] | **PASS** |
| a549_myllia | Secondary | `vcc2025_h1` | +0.0979 | [+0.0615, +0.1388] | **PASS** |

## `H-E2` / `SA-8`: Split-Half Top-50 (`resp_top50_split`) & Joint Cross-Mask `E1(split) − E1(all)`

| Target Line | Role | N_common | `E1(all_common)` | `E1(split)` [97.5% CI] | Joint `E1(split) − E1(all)` [97.5% CI] | `FCT(split)` [97.5% CI] | `Perm100 split` Q97.5 |
|---|---|---:|---:|---|---|---|---:|
| a549_myllia | Secondary | 32 | -0.0037 | +0.1105 [+0.0454, +0.1854] | +0.1143 [+0.0471, +0.1918] | +0.0955 [+0.0172, +0.1777] | +0.1324 |

## Reporting-Only Legacy Contrasts (Unpaired `H-D2` & Paired `FCT − E1_ConstMed`)

| Target Line | Role | Unpaired `FCT_lo − max(Mis_hi)` | `FCT − E1 (all_nt)` [97.5% CI] | `FCT − E1 (resp_top50)` [97.5% CI] |
|---|---|---:|---|---|
| a549_myllia | Secondary | +0.0096 (PASS) | +0.1114 [+0.0690, +0.1543] | -0.0456 [-0.0780, -0.0093] |

## Reporting-Only Full Line Metrics (`all_nontarget` & Naive `resp_top50`, 97.5% CIs)

| Target Line | Role | Method | N | `all_nontarget r_orth_cal` | 97.5% CI | Naive `resp_top50 r_orth_cal` | `rho_AB_orth` | `valid_frac` |
|---|---|---|---:|---:|---|---:|---:|---:|
| a549_myllia | Secondary | `Source_Grand_Mean` | 32 | +0.2834 | [+0.1727, +0.3828] | +0.3908 | +0.3875 | 1.0000 |
| a549_myllia | Secondary | `Transfer_Mean` | 32 | +nan | [+nan, +nan] | +nan | +0.3875 | 0.0000 |
| a549_myllia | Secondary | `FCT` | 32 | +0.1077 | [+0.0721, +0.1419] | +0.0827 | +0.3875 | 1.0000 |
| a549_myllia | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 32 | +0.0420 | [+0.0218, +0.0615] | +0.0312 | +0.3875 | 1.0000 |
| a549_myllia | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 32 | +0.0338 | [+0.0106, +0.0553] | +0.0514 | +0.3875 | 1.0000 |
| a549_myllia | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 32 | +0.0092 | [-0.0095, +0.0262] | +0.0594 | +0.3875 | 1.0000 |
| a549_myllia | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 32 | +0.0311 | [+0.0079, +0.0519] | -0.0129 | +0.3875 | 1.0000 |
| a549_myllia | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 32 | +0.0379 | [+0.0162, +0.0625] | +0.0099 | +0.3875 | 1.0000 |
| a549_myllia | Secondary | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 32 | +0.0097 | [-0.0144, +0.0364] | +0.1561 | +0.3875 | 1.0000 |
| a549_myllia | Secondary | `NegCtrl_FCT_Perm0` | 32 | -0.0017 | [-0.0181, +0.0169] | +0.1231 | +0.3875 | 1.0000 |
| a549_myllia | Secondary | `FCT_FullNTMean` | 32 | +0.1077 | [+0.0721, +0.1419] | +0.0827 | +0.3875 | 1.0000 |
| a549_myllia | Secondary | `Transfer_Mean_GlobalScale` | 32 | +nan | [+nan, +nan] | +nan | +0.3875 | 0.0000 |
| a549_myllia | Secondary | `Transfer_Mean_VarScaled` | 32 | +0.1029 | [+0.0734, +0.1299] | +0.1494 | +0.3875 | 1.0000 |
| a549_myllia | Secondary | `CIPHER_Linear` | 32 | +0.1488 | [+0.0509, +0.2379] | +0.2408 | +0.3875 | 1.0000 |
| a549_myllia | Secondary | `Hybrid_ProjShared` | 32 | +0.0853 | [+0.0384, +0.1297] | +0.1493 | +0.3875 | 1.0000 |
| a549_myllia | Secondary | `Hybrid_on_TMVS` | 32 | +0.1174 | [+0.0661, +0.1633] | +0.2006 | +0.3875 | 1.0000 |
| a549_myllia | Secondary | `E1_ConstMed` | 32 | -0.0037 | [-0.0255, +0.0208] | +0.1283 | +0.3875 | 1.0000 |
| a549_myllia | Secondary | `E2_ConstP10` | 32 | -0.0023 | [-0.0241, +0.0224] | +0.1221 | +0.3875 | 1.0000 |
| a549_myllia | Secondary | `E2_ConstP90` | 32 | -0.0068 | [-0.0288, +0.0172] | +0.1297 | +0.3875 | 1.0000 |
| a549_myllia | Secondary | `E3_SrcMeanCtrl` | 32 | +0.0310 | [+0.0114, +0.0506] | +0.0830 | +0.3875 | 1.0000 |
| a549_myllia | Secondary | `E5_RankMatchedMu` | 32 | +0.0375 | [+0.0167, +0.0580] | +0.0830 | +0.3875 | 1.0000 |
| a549_myllia | Secondary | `E6_LogFCBar` | 32 | -0.0097 | [-0.0286, +0.0096] | +0.1143 | +0.3875 | 1.0000 |
| a549_myllia | Secondary | `E7_FCBarMinus1` | 32 | -0.0015 | [-0.0231, +0.0224] | +0.1152 | +0.3875 | 1.0000 |

## Reporting-Only `N_5` Null-Cone Analysis — `all_nontarget` (`r_cone_sp`, 97.5% CIs)

*Note: `Transfer_Mean`, `Transfer_Mean_GlobalScale`, and `Source_Grand_Mean` are **zero by exact mathematical construction** under `N_5`, and `E1_ConstMed`, `E2_ConstP10`, `E2_ConstP90`, `E6_LogFCBar`, `E7_FCBarMinus1` are **suppressed by cone construction** (`R2_N5 >= 0.99`); only `FCT`, `Mis_s`, single-draw `Perm0..99`, `E3`/`E5`, and non-cone target models are empirically informative.*

| Target Line | Method | `N_5 r_cone_sp` | 97.5% CI | `eta` | `R2_N5` | Structural Role |
|---|---|---:|---|---:|---:|---|
| a549_myllia | `Source_Grand_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_myllia | `Transfer_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_myllia | `FCT` | +0.0839 | [+0.0521, +0.1236] | 0.6824 | 0.534310 | `INFORMATIVE` |
| a549_myllia | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | +0.0352 | [+0.0185, +0.0517] | 0.8344 | 0.303855 | `INFORMATIVE` |
| a549_myllia | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | +0.0137 | [+0.0026, +0.0246] | 0.7798 | 0.391953 | `INFORMATIVE` |
| a549_myllia | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | +0.0077 | [-0.0002, +0.0163] | 0.8667 | 0.248754 | `INFORMATIVE` |
| a549_myllia | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | +0.0169 | [+0.0022, +0.0301] | 0.8333 | 0.305563 | `INFORMATIVE` |
| a549_myllia | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | +0.0259 | [+0.0068, +0.0476] | 0.8312 | 0.309112 | `INFORMATIVE` |
| a549_myllia | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | +0.0031 | [-0.0075, +0.0117] | 0.7138 | 0.490456 | `INFORMATIVE` |
| a549_myllia | `NegCtrl_FCT_Perm0` | +0.0015 | [-0.0019, +0.0052] | 0.6406 | 0.589664 | `INFORMATIVE` |
| a549_myllia | `FCT_FullNTMean` | +0.0839 | [+0.0521, +0.1236] | 0.6824 | 0.534310 | `INFORMATIVE` |
| a549_myllia | `Transfer_Mean_GlobalScale` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_myllia | `Transfer_Mean_VarScaled` | +0.0878 | [+0.0624, +0.1144] | 0.8721 | 0.239386 | `INFORMATIVE` |
| a549_myllia | `CIPHER_Linear` | +0.0492 | [+0.0168, +0.0835] | 0.4319 | 0.813480 | `INFORMATIVE` |
| a549_myllia | `Hybrid_ProjShared` | +0.0964 | [+0.0412, +0.1526] | 0.9973 | 0.005297 | `INFORMATIVE` |
| a549_myllia | `Hybrid_on_TMVS` | +0.1233 | [+0.0675, +0.1771] | 0.9784 | 0.042734 | `INFORMATIVE` |
| a549_myllia | `E1_ConstMed` | -0.0000 | [-0.0001, +0.0000] | 0.0026 | 0.999993 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_myllia | `E2_ConstP10` | -0.0000 | [-0.0000, +0.0000] | 0.0012 | 0.999999 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_myllia | `E2_ConstP90` | +0.0000 | [-0.0000, +0.0001] | 0.0042 | 0.999983 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_myllia | `E3_SrcMeanCtrl` | +0.0240 | [+0.0138, +0.0350] | 0.7506 | 0.436673 | `INFORMATIVE` |
| a549_myllia | `E5_RankMatchedMu` | +0.0281 | [+0.0169, +0.0398] | 0.8295 | 0.311883 | `INFORMATIVE` |
| a549_myllia | `E6_LogFCBar` | -0.0002 | [-0.0009, +0.0006] | 0.0927 | 0.991404 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_myllia | `E7_FCBarMinus1` | +0.0000 | [-0.0000, +0.0001] | 0.0033 | 0.999989 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |

## Reporting-Only `N_5` Null-Cone Analysis — `resp_top50_split` (`r_cone_sp`, 97.5% CIs)

*Note: `Transfer_Mean`, `Transfer_Mean_GlobalScale`, and `Source_Grand_Mean` are **zero by exact mathematical construction** under `N_5`, and `E1_ConstMed`, `E2_ConstP10`, `E2_ConstP90`, `E6_LogFCBar`, `E7_FCBarMinus1` are **suppressed by cone construction** (`R2_N5 >= 0.99`); only `FCT`, `Mis_s`, single-draw `Perm0..99`, `E3`/`E5`, and non-cone target models are empirically informative.*

| Target Line | Method | `N_5 r_cone_sp` | 97.5% CI | `eta` | `R2_N5` | Structural Role |
|---|---|---:|---|---:|---:|---|
| a549_myllia | `Source_Grand_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_myllia | `Transfer_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_myllia | `FCT` | -0.0150 | [-0.0497, +0.0139] | 0.2556 | 0.931847 | `INFORMATIVE` |
| a549_myllia | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | -0.0051 | [-0.0367, +0.0352] | 0.5549 | 0.690593 | `INFORMATIVE` |
| a549_myllia | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | +0.0145 | [-0.0270, +0.0567] | 0.5992 | 0.637909 | `INFORMATIVE` |
| a549_myllia | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | -0.0040 | [-0.0550, +0.0485] | 0.6457 | 0.574874 | `INFORMATIVE` |
| a549_myllia | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | -0.0149 | [-0.0485, +0.0253] | 0.6085 | 0.629223 | `INFORMATIVE` |
| a549_myllia | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | -0.0023 | [-0.0412, +0.0414] | 0.5569 | 0.689837 | `INFORMATIVE` |
| a549_myllia | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | -0.0039 | [-0.0339, +0.0270] | 0.5384 | 0.696060 | `INFORMATIVE` |
| a549_myllia | `NegCtrl_FCT_Perm0` | +0.0069 | [-0.0238, +0.0332] | 0.4252 | 0.816229 | `INFORMATIVE` |
| a549_myllia | `FCT_FullNTMean` | -0.0150 | [-0.0497, +0.0139] | 0.2556 | 0.931847 | `INFORMATIVE` |
| a549_myllia | `Transfer_Mean_GlobalScale` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_myllia | `Transfer_Mean_VarScaled` | +0.0437 | [+0.0021, +0.0836] | 0.5829 | 0.660198 | `INFORMATIVE` |
| a549_myllia | `CIPHER_Linear` | +0.0510 | [+0.0096, +0.0963] | 0.3299 | 0.889870 | `INFORMATIVE` |
| a549_myllia | `Hybrid_ProjShared` | +0.1317 | [+0.0318, +0.2280] | 0.9120 | 0.168162 | `INFORMATIVE` |
| a549_myllia | `Hybrid_on_TMVS` | +0.1305 | [+0.0427, +0.2136] | 0.8064 | 0.349722 | `INFORMATIVE` |
| a549_myllia | `E1_ConstMed` | +0.0000 | [-0.0000, +0.0001] | 0.0005 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_myllia | `E2_ConstP10` | -0.0000 | [-0.0000, +0.0000] | 0.0002 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_myllia | `E2_ConstP90` | -0.0000 | [-0.0001, +0.0001] | 0.0012 | 0.999999 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_myllia | `E3_SrcMeanCtrl` | +0.0130 | [-0.0210, +0.0562] | 0.4983 | 0.744642 | `INFORMATIVE` |
| a549_myllia | `E5_RankMatchedMu` | +0.0123 | [-0.0199, +0.0516] | 0.5072 | 0.735911 | `INFORMATIVE` |
| a549_myllia | `E6_LogFCBar` | -0.0000 | [-0.0003, +0.0001] | 0.0015 | 0.999997 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_myllia | `E7_FCBarMinus1` | +0.0000 | [-0.0000, +0.0000] | 0.0011 | 0.999998 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |

## 100-Perm Empirical Null (`all_nontarget`, `SA-8`, & Single-Draw `N_5`)

- **a549_myllia**: `all_nontarget` med=-0.0036, Q97.5=+0.0035; `SA-8 split` med=+0.0851, Q97.5=+0.1324; `N_5` single-draw med=-0.0003, med|r|=+0.0018, Q97.5=+0.0075 (med_eta=0.6361)

## Mandatory Governance, Provenance & Prior-Art Disclosures (`PREREG_E §0.2, §5, §7.2, §8`)

- **Post-Hoc Provenance**: `H-E1` (paired `FCT - Mis_s` on `all_nontarget`) and `H-E2` (`SA-8 resp_top50_split` constant-baseline inflation) were specified in `PREREG_E` after unblinding Rounds 9-11 (`N=12` lines) and evaluated zero-shot on unseen Round 12 lines.
- **Claim 2 Non-Mechanistic Scope**: Confirmation of `H-E2` establishes an **empirical evaluation artifact strictly under `r_orth_cal`** (`resp_top50_split` vs `all_nontarget`) — namely that single-line constant baseline `E1_ConstMed` attains inflated disattenuated residual correlation under magnitude-selected subsets — and **must NOT be extrapolated** as a mechanistic biological claim or applied directly to raw DEG discovery benchmarks without independent validation.
- **Modality Boundary & Single-Line Generalizability**: Source cohorts (`r9_6`, `r9_r10_8`) comprise **CRISPRi** perturbations whereas `a549_ko` (`GSE345058`) is **CRISPRko**. Failure on `a549_ko` is confounded with CRISPRi->CRISPRko modality shift rather than purely same-modality failure; passage constitutes cross-modality generalization evidence. Under `K_conf = 1`, population generalizability across additional lineages remains explicitly limited.
- **Prior-Art & Data Attribution**: Baseline-aware perturbation transfer context acknowledges **Kawsar et al. (2026), `CanPert/CASTER` (`doi:10.64898/2026.11.15.766815`)**; `a549_ko` CRISPRko dataset derives from **Liu et al. (2026), `GSE345058`**; `hek293_xucao` CRISPRi derives from **Xu & Cao (2023), PerturbSci-Kinetics**.
- **Reporting-Only §5 Scope**: Legacy unpaired `H-D2`, paired `FCT - E1_ConstMed`, `N_5` null-cone suites (`all_nontarget` & `resp_top50_split`), `SA-Src8` (`--source-mode r9_r10_8`), `SA-Cas9` (`a549_ko_cas9pos`), naive `resp_top50`, and `H-D_util` are strictly descriptive and carry zero multiplicity weight.

## Pre-Registered Descriptive Subgroup Analysis (`PREREG_G` Item 7 — Reporting-Only)

*Descriptive-only stratification of evaluated targets by overlap with Round 12 `a549_ko` (`in_r12_a549_overlap` vs `novel_not_in_r12_a549`); zero confirmatory adjudication weight.*

| Target Line | Subgroup | N_targets | `FCT` (`all_nontarget` `r_orth_cal`) [97.5% CI] | `rho_AB_orth` | `valid_frac` | Min Paired `FCT − Mis(s')` CI97.5_lo |
|---|---|---:|---|---:|---:|---:|
| a549_myllia | In Round 12 `a549_ko` Overlap | 4 | +0.1058 [+0.0609, +0.1338] | +0.5386 | 1.0000 | +0.0034 |
| a549_myllia | Novel (Not Tested in Round 12 `a549_ko`) | 28 | +0.1079 [+0.0704, +0.1436] | +0.3755 | 1.0000 | +0.0397 |

| Target Line | Subgroup | N_targets | Mismatched Source | Paired `FCT − Mis(s')` (`all_nontarget`) [97.5% CI] | `valid_frac` |
|---|---|---:|---|---|---:|
| a549_myllia | Overlap (`R12 a549_ko`) | 4 | `nadig_hepg2` | +0.0538 [+0.0034, +0.1029] | 1.0000 |
| a549_myllia | Overlap (`R12 a549_ko`) | 4 | `nadig_jurkat` | +0.0521 [+0.0461, +0.0711] | 1.0000 |
| a549_myllia | Overlap (`R12 a549_ko`) | 4 | `nourreddine_kolf21j` | +0.1023 [+0.0582, +0.1351] | 1.0000 |
| a549_myllia | Overlap (`R12 a549_ko`) | 4 | `replogle_k562_essential` | +0.0452 [+0.0269, +0.0706] | 1.0000 |
| a549_myllia | Overlap (`R12 a549_ko`) | 4 | `replogle_rpe1` | +0.0641 [+0.0056, +0.1248] | 0.9960 |
| a549_myllia | Overlap (`R12 a549_ko`) | 4 | `vcc2025_h1` | +0.1075 [+0.0802, +0.1500] | 1.0000 |
| a549_myllia | Novel (Not in `R12 a549_ko`) | 28 | `nadig_hepg2` | +0.0669 [+0.0397, +0.0958] | 1.0000 |
| a549_myllia | Novel (Not in `R12 a549_ko`) | 28 | `nadig_jurkat` | +0.0773 [+0.0467, +0.1061] | 1.0000 |
| a549_myllia | Novel (Not in `R12 a549_ko`) | 28 | `nourreddine_kolf21j` | +0.0991 [+0.0594, +0.1351] | 1.0000 |
| a549_myllia | Novel (Not in `R12 a549_ko`) | 28 | `replogle_k562_essential` | +0.0815 [+0.0458, +0.1151] | 1.0000 |
| a549_myllia | Novel (Not in `R12 a549_ko`) | 28 | `replogle_rpe1` | +0.0701 [+0.0421, +0.0976] | 1.0000 |
| a549_myllia | Novel (Not in `R12 a549_ko`) | 28 | `vcc2025_h1` | +0.0983 [+0.0576, +0.1412] | 1.0000 |

## Round 14 (`PREREG_G`) Governance & Interpretation Rules

- **Protocol Inheritance**: Claims 1–2 (`H-E1`, `H-E2`), `r_orth_cal`, `all_nontarget` / `resp_top50_split`, `r9_6`, 2-claim Bonferroni + IUT (`97.5%` CIs), `4,000` bootstraps, and `100` permutations are inherited identically from locked `PREREG_E`; seed=`20261224`.
- **Same-Cell-Line Independent-Dataset Scope ONLY**: Round 14 tests `a549_myllia` (`GSE337804`, Adduri et al. *Cell* 2026 / Myllia) independently under `K_conf <= 1` (`ch_1 -> Half A`, `ch_2 -> Half B`). Confirmation constitutes replication across two independent A549 datasets in two separately pre-registered prospective rounds (NEVER pooled with Round 12 `a549_ko`, and NEVER claimed as cross-cell-line generalization).
- **Chemistry & Technical Context**: 10x `5' HT v2` + Jumpcode `DepleteX`, `SpCas9` CRISPRn (`3` guides/gene), `20` `non-targeting` NTCs (`ctrl_T` `gene_desert` excluded), RNA QC `umi >= 300 & ngenes >= 200`, `max_pct_mt=None`.
