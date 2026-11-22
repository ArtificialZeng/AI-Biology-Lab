# Round 11 Confirmatory Group D (`FCT`, Primary Mask `all_nontarget`) Report (`dev_r9_r9_6_log1p_smoke`)

- **Phase**: `dev_r9` | **Source Mode**: `r9_6` (`6` sources) | **PREREG_D Commit**: `` | **Code Freeze**: ``
- **Bootstrap**: `np.random.default_rng(20261123)` x `4000` (Bonferroni co-primary `98.33%` CIs `[0.8333%, 99.1667%]`)
- **Cohort Disclosure**: **Insufficient confirmatory cell lines (K_conf < 2); confirmatory hypothesis testing was not performed (确认性线不足，未做确认性检验).** (`[]`)
- **Secondary Lines**: `['vcc2025_h1']`
- **Primary Group D (`FCT` 全基因线特异信号) Verdict**: **NOT_EVALUATED_INSUFFICIENT_LINES**
- **Secondary Methodological Hypothesis (`H-D4`)**: **NOT_EVALUATED_INSUFFICIENT_LINES** (all lines descriptive: `CONFIRMED`)
- **Parallel Practical Utility (`H-C3`-style Non-Inferiority)**: **NOT_EVALUATED_INSUFFICIENT_LINES**

## Per-Line Primary Mask (`all_nontarget`, Joint-Bootstrap 98.33% CIs) & Secondary (`resp_top50`)

| Target Line | Role | Method | N | Primary `r_orth_cal` | Primary 98.33% CI | `resp_top50` `r_orth_cal` | `resp_top50` 98.33% CI | `rho_AB_orth` |
|---|---|---|---:|---:|---|---:|---|---:|
| vcc2025_h1 | Secondary | `Source_Grand_Mean` | 80 | +0.0446 | [+0.0131, +0.0768] | +0.1086 | [+0.0395, +0.1733] | +0.8886 |
| vcc2025_h1 | Secondary | `Transfer_Mean` | 80 | +nan | [+nan, +nan] | +nan | [+nan, +nan] | +0.8886 |
| vcc2025_h1 | Secondary | `FCT` | 80 | +0.1427 | [+0.1155, +0.1733] | +0.3505 | [+0.2754, +0.4246] | +0.8886 |
| vcc2025_h1 | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_hepg2)` | 80 | +0.0451 | [+0.0353, +0.0564] | +0.1240 | [+0.0811, +0.1716] | +0.8886 |
| vcc2025_h1 | Secondary | `NegCtrl_FCT_MismatchedMean(nadig_jurkat)` | 80 | +0.0180 | [+0.0073, +0.0294] | +0.0550 | [+0.0100, +0.0994] | +0.8886 |
| vcc2025_h1 | Secondary | `NegCtrl_FCT_MismatchedMean(nourreddine_kolf21j)` | 80 | +0.0537 | [+0.0411, +0.0680] | +0.1453 | [+0.0873, +0.2033] | +0.8886 |
| vcc2025_h1 | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_k562_essential)` | 80 | +0.0281 | [+0.0176, +0.0387] | +0.0536 | [+0.0100, +0.0936] | +0.8886 |
| vcc2025_h1 | Secondary | `NegCtrl_FCT_MismatchedMean(replogle_rpe1)` | 80 | -0.0134 | [-0.0298, +0.0003] | -0.0541 | [-0.1122, +0.0088] | +0.8886 |
| vcc2025_h1 | Secondary | `NegCtrl_FCT_Perm0` | 80 | +0.0318 | [+0.0193, +0.0465] | +0.1802 | [+0.1283, +0.2299] | +0.8886 |
| vcc2025_h1 | Secondary | `Transfer_Mean_GlobalScale` | 80 | +nan | [+nan, +nan] | +nan | [+nan, +nan] | +0.8886 |
| vcc2025_h1 | Secondary | `Transfer_Mean_VarScaled` | 80 | +0.0323 | [-0.0078, +0.0735] | -0.0772 | [-0.1905, +0.0410] | +0.8886 |
| vcc2025_h1 | Secondary | `CIPHER_Linear` | 80 | +0.0022 | [-0.0271, +0.0303] | -0.0468 | [-0.1487, +0.0452] | +0.8886 |
| vcc2025_h1 | Secondary | `Hybrid_ProjShared` | 80 | +0.0067 | [-0.0169, +0.0292] | -0.0403 | [-0.1358, +0.0451] | +0.8886 |
| vcc2025_h1 | Secondary | `Hybrid_on_TMVS` | 80 | +0.0208 | [-0.0068, +0.0468] | -0.0819 | [-0.1862, +0.0162] | +0.8886 |
| vcc2025_h1 | Secondary | `E1_ConstMed` | 80 | +0.0437 | [+0.0272, +0.0632] | +0.3076 | [+0.2383, +0.3791] | +0.8886 |

## `H-D3` (`FCT − E1`) & `H-D4` (`E1(resp_top50) − E1(all_nontarget)`) Paired Comparisons

| Target Line | Comparison | Estimate | 95% CI | 98.33% CI |
|---|---|---:|---|---|
| vcc2025_h1 | `H-D3_FCT_minus_E1_all_nontarget` | +0.0989 | [+0.0852, +0.1141] | [+0.0819, +0.1186] |
| vcc2025_h1 | `diag_FCT_minus_E1_resp_top50` | +0.0429 | [+0.0156, +0.0725] | [+0.0105, +0.0787] |
| vcc2025_h1 | `H-D4_E1_resp_top50_minus_all_nontarget` | +0.2639 | [+0.2120, +0.3147] | [+0.2029, +0.3242] |

## 100-Perm Null (`all_nontarget` Primary & `resp_top50` Secondary)

- **vcc2025_h1**: `all_nontarget` median=+0.0322, p97.5=+0.0386; `resp_top50` median=+0.2101, p97.5=+0.2445
