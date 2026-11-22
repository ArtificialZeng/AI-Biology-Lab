# Round 12 Independent Pre-Registration Audit Report (`dev_r11`)

- **Source Mode:** `r9_6` | **Seed:** `20261201` | **Bootstraps:** `400` (`97.5%` two-sided Bonferroni CIs) | **Perms:** `20`
- **Adjudication Status:** `DEV_DRY_RUN_DEV_R11` (`K_conf = 0`, confirmatory lines = `[]`)
- **Claim 1 (`H-E1`):** `DEV_ONLY_NO_ADJUDICATION` | **Claim 2 (`H-E2`):** `DEV_ONLY_NO_ADJUDICATION`

## Mandatory Pre-Registered Disclosures (`PREREG_E §0.1, §0.2, §7, §8`)
- **Post-Hoc Provenance (`E6`):** Both `H-E1` (paired `FCT - Mis_s`) and `H-E2` (`SA-8`) hypotheses — and the interpretation that prior T-cell unpaired `H-D2` failures reflected low unpaired power — were formulated after inspecting Rounds 9–11 data.
- **Claim 2 Empirical & Metric Scope (`E6`):** `SA-8` is strictly an empirical finding on residual metric `r_orth_cal` (`resp_top50_split`), NOT mechanistic proof (Round 11 Check 4 noiseless simulation predicted the wrong sign on `8/12` lines under `all_nontarget`), and does NOT extrapolate to Nicol et al. (2026) DEG metrics.
- **Claim 1 Novelty Scope (`§8`):** Must cite Kawsar et al. 2026 (`CanPert/CASTER`, `doi:10.64898/2026.11.15.766815`) for prior per-target paired matched-vs-mismatched control comparison; Claim 1 tests prospectively whether baseline-only `FCT` outperforms all `6` mismatched-line baselines on residual `r_orth_cal` in an unseen lineage and modality.
- **Scope & Modality:** Pre-freeze verification strictly on unblinded historical R9-R11 lines.
- **CRISPRi -> CRISPRko Modality Boundary:** DEV_R9_TO_R11_VALIDATION

## Confirmatory Endpoints (`H-E1` & `H-E2`, `r_orth_cal`)
| Line (`N`) | `FCT` `[97.5% CI]` (`vf`, `rho`) | `H-E1a` | `H-E1b` Paired Pass | `H-E1c` (`> Perm97.5`) | `H-E1` | `E1(split)` `[97.5% CI]` | `E1(split)-E1(all)` `[97.5% CI]` (`N_c`) | `H-E2` |
|---|---|:---:|:---:|:---:|:---:|---|---|:---:|
| `feng_ipsc` (`93`) | `+0.1212` `[+0.1061, +0.1382]` (`1.000`, `+0.557`) | `True` | `6/6` (`True`) | `True` (`+0.0434`) | **`True`** | `+0.2468` `[+0.1994, +0.2988]` | `+0.1938` `[+0.1523, +0.2416]` (`93`) | **`True`** |

## `H-E1b` Source-by-Source Paired `FCT - Mis_s` (`all_nontarget`, `r_orth_cal`)
| Target | Source `s'` | `Mis_s` `[97.5% CI]` | Paired `FCT - Mis_s` `[97.5% CI]` | Pass |
|---|---|---|---|---|
| `feng_ipsc` | `nadig_hepg2` | `+0.0437` `[+0.0321, +0.0560]` | `+0.0775` `[+0.0655, +0.0910]` | `True` |
| `feng_ipsc` | `nadig_jurkat` | `-0.0024` `[-0.0132, +0.0083]` | `+0.1236` `[+0.1049, +0.1404]` | `True` |
| `feng_ipsc` | `nourreddine_kolf21j` | `+0.0562` `[+0.0362, +0.0766]` | `+0.0650` `[+0.0537, +0.0794]` | `True` |
| `feng_ipsc` | `replogle_k562_essential` | `+0.0069` `[-0.0031, +0.0174]` | `+0.1143` `[+0.0966, +0.1328]` | `True` |
| `feng_ipsc` | `replogle_rpe1` | `-0.0234` `[-0.0350, -0.0114]` | `+0.1446` `[+0.1256, +0.1623]` | `True` |
| `feng_ipsc` | `vcc2025_h1` | `+0.0669` `[+0.0516, +0.0837]` | `+0.0543` `[+0.0417, +0.0667]` | `True` |

## §5 Non-Confirmatory Descriptive Summary (`只报告、不判定`)
*Explicit notice (`PREREG_E §5.3`): Within `N_5`, `Transfer_Mean`, `Transfer_Mean_GlobalScale` (`eta=0`), and `Source_Grand_Mean` (`eta=0`) equal zero by exact algebraic construction, and `E1_ConstMed`, `E2_ConstP10/P90`, `E6_LogFCBar`, `E7_FCBarMinus1` are suppressed by cone construction — none constitute empirical validation evidence.*
- **`feng_ipsc`**: Legacy unpaired `H-D2`=`True`; `FCT - E1 (all)`=`+0.0682` `[+0.0412, +0.0896]`; `N_5(all) FCT`=`+0.1325` (`eta=0.906`), `GS`=`+0.0000` (`eta=1.8e-15`), `E1`=`-0.0000` (`eta=0.0020`); `H-D_util` pass=`False` (`R^2(FCT orth/all)` approx via unprojected `r_cal^2` reported in JSON).
