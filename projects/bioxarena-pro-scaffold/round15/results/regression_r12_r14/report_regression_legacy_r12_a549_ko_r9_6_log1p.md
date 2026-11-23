# Round 12 Prospective Group E Report (`regression_legacy_r12_a549_ko_r9_6_log1p`)

- **Phase**: `regression_legacy` | **Source Mode**: `r9_6` (`6` sources) | **PREREG_E Commit**: `` | **Code Freeze**: ``
- **Bootstrap**: `np.random.default_rng(20261201)` x `4000` (Bonferroni 2-claim + IUT two-sided `97.5%` CIs `[1.25%, 98.75%]`)
- **Cohort Disclosure**: **K_conf = 1 confirmatory line** (`['a549_ko']`) — confirmatory inference permitted under `PREREG_E §1`; **mandatory disclosure: single-line generalizability is limited.**
- **Confirmatory Claim 1 (`H-E1` Target-Line Specificity)**: **CONFIRMED**
- **Confirmatory Claim 2 (`H-E2` Persistent `SA-8` Constant-Baseline Inflation)**: **CONFIRMED**

## Confirmatory Hypothesis Adjudication (`97.5%` Two-Sided CIs via IUT + 2-Claim Bonferroni)

| Claim / Hypothesis | Pass / Total | Criterion | Verdict |
|---|---:|---|---|
| **Claim 1 (`H-E1`)** | 1/1 | `H-E1a ∩ H-E1b ∩ H-E1c` (Unanimous) | **CONFIRMED** |
| ├─ `H-E1a` (`FCT` all_nt 97.5% lo > 0, valid >= 0.95, rho_orth > 0) | 1/1 | Addendum A1 | **CONFIRMED** |
| ├─ `H-E1b` (Paired `FCT − Mis_s` 97.5% lo > 0, ALL sources) | 1/1 | Source-wise IUT | **CONFIRMED** |
| └─ `H-E1c` (`FCT` 97.5% lo > Perm100 `all_nontarget` Q0.975) | 1/1 | Perm null | **CONFIRMED** |
| **Claim 2 (`H-E2`)** | 1/1 | `H-E2a ∩ H-E2b` (Unanimous) | **CONFIRMED** |
| ├─ `H-E2a` (`E1_ConstMed` under `resp_top50_split` 97.5% lo > 0) | 1/1 | SA-8 inflation | **CONFIRMED** |
| └─ `H-E2b` (Joint `E1(resp_top50_split) − E1(all_nontarget)` 97.5% lo > 0) | 1/1 | Common-pert IUT | **CONFIRMED** |

## `H-E1b`: Paired Specificity Contrasts `FCT − NegCtrl_FCT_MismatchedMean(s)` (`all_nontarget`, 97.5% CI)

| Target Line | Role | Mismatched Source | Estimate | 97.5% CI | `lo > 0` |
|---|---|---|---:|---|---|
| a549_ko | Confirmatory | `nadig_hepg2` | +0.1869 | [+0.1693, +0.2049] | **PASS** |
| a549_ko | Confirmatory | `nadig_jurkat` | +0.0796 | [+0.0693, +0.0900] | **PASS** |
| a549_ko | Confirmatory | `nourreddine_kolf21j` | +0.1069 | [+0.0985, +0.1156] | **PASS** |
| a549_ko | Confirmatory | `replogle_k562_essential` | +0.1279 | [+0.1150, +0.1410] | **PASS** |
| a549_ko | Confirmatory | `replogle_rpe1` | +0.0804 | [+0.0682, +0.0917] | **PASS** |
| a549_ko | Confirmatory | `vcc2025_h1` | +0.1505 | [+0.1402, +0.1605] | **PASS** |

## `H-E2` / `SA-8`: Split-Half Top-50 (`resp_top50_split`) & Joint Cross-Mask `E1(split) − E1(all)`

| Target Line | Role | N_common | `E1(all_common)` | `E1(split)` [97.5% CI] | Joint `E1(split) − E1(all)` [97.5% CI] | `FCT(split)` [97.5% CI] | `Perm100 split` Q97.5 |
|---|---|---:|---:|---|---|---|---:|
| a549_ko | Confirmatory | 445 | -0.0201 | +0.2565 [+0.2195, +0.2981] | +0.2766 [+0.2412, +0.3164] | +0.2965 [+0.2633, +0.3310] | +0.2108 |

## Reporting-Only Legacy Contrasts (Unpaired `H-D2` & Paired `FCT − E1_ConstMed`)

| Target Line | Role | Unpaired `FCT_lo − max(Mis_hi)` | `FCT − E1 (all_nt)` [97.5% CI] | `FCT − E1 (resp_top50)` [97.5% CI] |
|---|---|---:|---|---|
| a549_ko | Confirmatory | +0.0504 (PASS) | +0.2099 [+0.1975, +0.2223] | +0.0207 [-0.0078, +0.0516] |

## Reporting-Only Full Line Metrics (`all_nontarget` & Naive `resp_top50`, 97.5% CIs)

| Target Line | Role | Method | N | `all_nontarget r_orth_cal` | 97.5% CI | Naive `resp_top50 r_orth_cal` | `rho_AB_orth` | `valid_frac` |
|---|---|---|---:|---:|---|---:|---:|---:|
| a549_ko | Confirmatory | `Source_Grand_Mean` | 445 | +0.4051 | [+0.3814, +0.4287] | +0.4520 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `Transfer_Mean` | 445 | +nan | [+nan, +nan] | +nan | +0.7887 | 0.0000 |
| a549_ko | Confirmatory | `FCT` | 445 | +0.1898 | [+0.1742, +0.2055] | +0.3073 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 445 | +0.0029 | [-0.0053, +0.0114] | +0.0608 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 445 | +0.1102 | [+0.0977, +0.1237] | +0.0957 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 445 | +0.0829 | [+0.0734, +0.0927] | +0.1139 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 445 | +0.0619 | [+0.0519, +0.0723] | -0.0184 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 445 | +0.1094 | [+0.0966, +0.1231] | +0.1284 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | 445 | +0.0393 | [+0.0303, +0.0483] | +0.2244 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `NegCtrl_FCT_Perm0` | 445 | -0.0012 | [-0.0080, +0.0057] | +0.2068 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `FCT_FullNTMean` | 445 | +0.1898 | [+0.1742, +0.2055] | +0.3073 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `Transfer_Mean_GlobalScale` | 445 | +nan | [+nan, +nan] | +nan | +0.7887 | 0.0000 |
| a549_ko | Confirmatory | `Transfer_Mean_VarScaled` | 445 | +0.2217 | [+0.2071, +0.2359] | +0.1871 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `CIPHER_Linear` | 445 | +0.2817 | [+0.2408, +0.3238] | +0.3128 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `Hybrid_ProjShared` | 445 | +0.1554 | [+0.1303, +0.1820] | +0.2307 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `Hybrid_on_TMVS` | 445 | +0.2373 | [+0.2127, +0.2625] | +0.2767 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `E1_ConstMed` | 445 | -0.0201 | [-0.0293, -0.0108] | +0.2866 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `E2_ConstP10` | 445 | -0.0229 | [-0.0318, -0.0137] | +0.2412 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `E2_ConstP90` | 445 | -0.0141 | [-0.0239, -0.0042] | +0.3257 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `E3_SrcMeanCtrl` | 445 | +0.0794 | [+0.0674, +0.0916] | +0.1845 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `E5_RankMatchedMu` | 445 | +0.1022 | [+0.0912, +0.1142] | +0.1224 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `E6_LogFCBar` | 445 | -0.0056 | [-0.0163, +0.0048] | +0.3167 | +0.7887 | 1.0000 |
| a549_ko | Confirmatory | `E7_FCBarMinus1` | 445 | -0.0241 | [-0.0329, -0.0150] | +0.2009 | +0.7887 | 1.0000 |

## Reporting-Only `N_5` Null-Cone Analysis — `all_nontarget` (`r_cone_sp`, 97.5% CIs)

*Note: `Transfer_Mean`, `Transfer_Mean_GlobalScale`, and `Source_Grand_Mean` are **zero by exact mathematical construction** under `N_5`, and `E1_ConstMed`, `E2_ConstP10`, `E2_ConstP90`, `E6_LogFCBar`, `E7_FCBarMinus1` are **suppressed by cone construction** (`R2_N5 >= 0.99`); only `FCT`, `Mis_s`, single-draw `Perm0..99`, `E3`/`E5`, and non-cone target models are empirically informative.*

| Target Line | Method | `N_5 r_cone_sp` | 97.5% CI | `eta` | `R2_N5` | Structural Role |
|---|---|---:|---|---:|---:|---|
| a549_ko | `Source_Grand_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko | `Transfer_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko | `FCT` | +0.1992 | [+0.1850, +0.2133] | 0.9350 | 0.125788 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | -0.0114 | [-0.0167, -0.0057] | 0.9336 | 0.128458 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | +0.0528 | [+0.0474, +0.0586] | 0.8569 | 0.265742 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | +0.0638 | [+0.0581, +0.0698] | 0.8980 | 0.193679 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | +0.0175 | [+0.0121, +0.0228] | 0.9276 | 0.139574 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | +0.0716 | [+0.0637, +0.0799] | 0.9355 | 0.124827 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | +0.0236 | [+0.0214, +0.0260] | 0.7829 | 0.387123 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_Perm0` | +0.0091 | [+0.0078, +0.0103] | 0.7029 | 0.505865 | `INFORMATIVE` |
| a549_ko | `FCT_FullNTMean` | +0.1992 | [+0.1850, +0.2133] | 0.9350 | 0.125788 | `INFORMATIVE` |
| a549_ko | `Transfer_Mean_GlobalScale` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko | `Transfer_Mean_VarScaled` | +0.1273 | [+0.1190, +0.1354] | 0.5880 | 0.654241 | `INFORMATIVE` |
| a549_ko | `CIPHER_Linear` | +0.0706 | [+0.0581, +0.0840] | 0.3199 | 0.897687 | `INFORMATIVE` |
| a549_ko | `Hybrid_ProjShared` | +0.1709 | [+0.1427, +0.2018] | 0.9926 | 0.014715 | `INFORMATIVE` |
| a549_ko | `Hybrid_on_TMVS` | +0.2049 | [+0.1821, +0.2282] | 0.8251 | 0.319223 | `INFORMATIVE` |
| a549_ko | `E1_ConstMed` | +0.0001 | [+0.0000, +0.0001] | 0.0028 | 0.999992 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko | `E2_ConstP10` | +0.0000 | [+0.0000, +0.0000] | 0.0007 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko | `E2_ConstP90` | -0.0000 | [-0.0000, -0.0000] | 0.0020 | 0.999996 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko | `E3_SrcMeanCtrl` | +0.0440 | [+0.0396, +0.0487] | 0.7678 | 0.410512 | `INFORMATIVE` |
| a549_ko | `E5_RankMatchedMu` | +0.0555 | [+0.0500, +0.0613] | 0.9400 | 0.116366 | `INFORMATIVE` |
| a549_ko | `E6_LogFCBar` | +0.0000 | [-0.0000, +0.0000] | 0.0123 | 0.999849 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko | `E7_FCBarMinus1` | -0.0001 | [-0.0002, -0.0001] | 0.0052 | 0.999973 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |

## Reporting-Only `N_5` Null-Cone Analysis — `resp_top50_split` (`r_cone_sp`, 97.5% CIs)

*Note: `Transfer_Mean`, `Transfer_Mean_GlobalScale`, and `Source_Grand_Mean` are **zero by exact mathematical construction** under `N_5`, and `E1_ConstMed`, `E2_ConstP10`, `E2_ConstP90`, `E6_LogFCBar`, `E7_FCBarMinus1` are **suppressed by cone construction** (`R2_N5 >= 0.99`); only `FCT`, `Mis_s`, single-draw `Perm0..99`, `E3`/`E5`, and non-cone target models are empirically informative.*

| Target Line | Method | `N_5 r_cone_sp` | 97.5% CI | `eta` | `R2_N5` | Structural Role |
|---|---|---:|---|---:|---:|---|
| a549_ko | `Source_Grand_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko | `Transfer_Mean` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko | `FCT` | +0.1301 | [+0.1099, +0.1518] | 0.6330 | 0.599272 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | -0.0374 | [-0.0509, -0.0237] | 0.7452 | 0.444621 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | +0.0170 | [-0.0006, +0.0330] | 0.8191 | 0.328820 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | +0.0495 | [+0.0359, +0.0634] | 0.8178 | 0.331200 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | +0.0023 | [-0.0127, +0.0169] | 0.8053 | 0.351417 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | +0.0041 | [-0.0133, +0.0216] | 0.7021 | 0.507010 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` | +0.0438 | [+0.0324, +0.0550] | 0.7157 | 0.487785 | `INFORMATIVE` |
| a549_ko | `NegCtrl_FCT_Perm0` | +0.0432 | [+0.0354, +0.0515] | 0.5730 | 0.671689 | `INFORMATIVE` |
| a549_ko | `FCT_FullNTMean` | +0.1301 | [+0.1099, +0.1518] | 0.6330 | 0.599272 | `INFORMATIVE` |
| a549_ko | `Transfer_Mean_GlobalScale` | +0.0000 | [+0.0000, +0.0000] | 0.0000 | 1.000000 | `ZERO_BY_EXACT_CONSTRUCTION` |
| a549_ko | `Transfer_Mean_VarScaled` | +0.1458 | [+0.1228, +0.1719] | 0.7164 | 0.486631 | `INFORMATIVE` |
| a549_ko | `CIPHER_Linear` | +0.0599 | [+0.0462, +0.0753] | 0.2838 | 0.919480 | `INFORMATIVE` |
| a549_ko | `Hybrid_ProjShared` | +0.1649 | [+0.1275, +0.2024] | 0.8701 | 0.242934 | `INFORMATIVE` |
| a549_ko | `Hybrid_on_TMVS` | +0.2070 | [+0.1726, +0.2424] | 0.8067 | 0.349245 | `INFORMATIVE` |
| a549_ko | `E1_ConstMed` | -0.0001 | [-0.0002, -0.0001] | 0.0015 | 0.999998 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko | `E2_ConstP10` | -0.0000 | [-0.0001, -0.0000] | 0.0004 | 1.000000 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko | `E2_ConstP90` | +0.0002 | [+0.0001, +0.0003] | 0.0023 | 0.999995 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko | `E3_SrcMeanCtrl` | +0.0078 | [-0.0065, +0.0220] | 0.6963 | 0.515032 | `INFORMATIVE` |
| a549_ko | `E5_RankMatchedMu` | +0.0122 | [-0.0041, +0.0282] | 0.7856 | 0.382769 | `INFORMATIVE` |
| a549_ko | `E6_LogFCBar` | -0.0001 | [-0.0002, +0.0000] | 0.0026 | 0.999993 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |
| a549_ko | `E7_FCBarMinus1` | +0.0004 | [+0.0002, +0.0007] | 0.0073 | 0.999946 | `SUPPRESSED_BY_CONE_CONSTRUCTION` |

## 100-Perm Empirical Null (`all_nontarget`, `SA-8`, & Single-Draw `N_5`)

- **a549_ko**: `all_nontarget` med=-0.0106, Q97.5=+0.0018; `SA-8 split` med=+0.1347, Q97.5=+0.2108; `N_5` single-draw med=+0.0004, med|r|=+0.0036, Q97.5=+0.0112 (med_eta=0.6975)

## Mandatory Governance, Provenance & Prior-Art Disclosures (`PREREG_E §0.2, §5, §7.2, §8`)

- **Post-Hoc Provenance**: `H-E1` (paired `FCT - Mis_s` on `all_nontarget`) and `H-E2` (`SA-8 resp_top50_split` constant-baseline inflation) were specified in `PREREG_E` after unblinding Rounds 9-11 (`N=12` lines) and evaluated zero-shot on unseen Round 12 lines.
- **Claim 2 Non-Mechanistic Scope**: Confirmation of `H-E2` establishes an **empirical evaluation artifact strictly under `r_orth_cal`** (`resp_top50_split` vs `all_nontarget`) — namely that single-line constant baseline `E1_ConstMed` attains inflated disattenuated residual correlation under magnitude-selected subsets — and **must NOT be extrapolated** as a mechanistic biological claim or applied directly to raw DEG discovery benchmarks without independent validation.
- **Modality Boundary & Single-Line Generalizability**: Source cohorts (`r9_6`, `r9_r10_8`) comprise **CRISPRi** perturbations whereas `a549_ko` (`GSE345058`) is **CRISPRko**. Failure on `a549_ko` is confounded with CRISPRi->CRISPRko modality shift rather than purely same-modality failure; passage constitutes cross-modality generalization evidence. Under `K_conf = 1`, population generalizability across additional lineages remains explicitly limited.
- **Prior-Art & Data Attribution**: Baseline-aware perturbation transfer context acknowledges **Kawsar et al. (2026), `CanPert/CASTER` (`doi:10.64898/2026.11.15.766815`)**; `a549_ko` CRISPRko dataset derives from **Liu et al. (2026), `GSE345058`**; `hek293_xucao` CRISPRi derives from **Xu & Cao (2023), PerturbSci-Kinetics**.
- **Reporting-Only §5 Scope**: Legacy unpaired `H-D2`, paired `FCT - E1_ConstMed`, `N_5` null-cone suites (`all_nontarget` & `resp_top50_split`), `SA-Src8` (`--source-mode r9_r10_8`), `SA-Cas9` (`a549_ko_cas9pos`), naive `resp_top50`, and `H-D_util` are strictly descriptive and carry zero multiplicity weight.

---

## Round 15 (`PREREG_H` v2) Governance, Multi-Line, & Wave Disclosures

- **Protocol**: `PREREG_H` v2 (`seed = 20261201`, `n_boot = 4000`, `n_perm = 100`)
- **Release Wave**: `r12_a549_ko` (`interim, more lines pending`) | **Pending Candidate Lines**: `['cal33', 'hela', 'ac16', 'molm13', 'u87mg']` | **Excluded `jurkat`**: `True`
- **Mandatory Scope Notice (`§8.2`)**: K_conf = 1 single-line confirmation only; cross-line generalizability is limited.
- **Arc Challenge Curation Disclosure (`§4.2` / `§10.3`)**: Competition challenge panels filter targets by `>= 80%` median knockdown, select the stronger of 2 sgRNAs per gene, and uniformly downsample to `400` cells/perturbation.
- **Prior-Art Notice (`§10.6`)**: Inspect NeurIPS 2026 `OpenReview lnpzb1HvwX` and `CanPert/CASTER` boundaries (`PRIOR_ART_H.md` / `PRIOR_ART_H_CONTROLLER_FAN2026.md`) cleanly.

## `SA-H1` (Leave-`vcc2025_h1`-Out Source Sensitivity — `5` Remaining Sources, Strictly Descriptive)

| Target Line | `n_targets` | `FCT r_orth_cal` [97.5% CI] | `rho_half_ab_orth` | Paired `FCT - Mis_s` Ranges (Min/Max 97.5% CI lo) |
|---|---:|---|---:|---|
| `a549_ko` | 439 | +0.1898 [+0.1741, +0.2059] | +0.7895 | [+0.0691, +0.1705] (`5` sources) |

## `SA-Panel` (`129` Shared Val-Test Targets vs Remainder — Strictly Descriptive)

- **Status**: `PANEL_OVERLAP_LIST_NOT_AVAILABLE` (`path=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round15/shared/vcc2026_val_test_overlap_129.txt`).
