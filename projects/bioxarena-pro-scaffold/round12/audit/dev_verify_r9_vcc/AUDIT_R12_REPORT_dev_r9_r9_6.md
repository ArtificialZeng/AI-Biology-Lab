# Round 12 Independent Pre-Registration Audit Report (`dev_r9`)

- **Source Mode:** `r9_6` | **Seed:** `20261201` | **Bootstraps:** `400` (`97.5%` two-sided Bonferroni CIs) | **Perms:** `20`
- **Adjudication Status:** `DEV_DRY_RUN_DEV_R9` (`K_conf = 0`, confirmatory lines = `[]`)
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
| `vcc2025_h1` (`80`) | `+0.1427` `[+0.1200, +0.1712]` (`1.000`, `+0.889`) | `True` | `5/5` (`True`) | `True` (`+0.0363`) | **`True`** | `+0.3100` `[+0.2552, +0.3757]` | `+0.2663` `[+0.2107, +0.3224]` (`80`) | **`True`** |

## `H-E1b` Source-by-Source Paired `FCT - Mis_s` (`all_nontarget`, `r_orth_cal`)
| Target | Source `s'` | `Mis_s` `[97.5% CI]` | Paired `FCT - Mis_s` `[97.5% CI]` | Pass |
|---|---|---|---|---|
| `vcc2025_h1` | `nadig_hepg2` | `+0.0451` `[+0.0357, +0.0554]` | `+0.0976` `[+0.0753, +0.1272]` | `True` |
| `vcc2025_h1` | `nadig_jurkat` | `+0.0180` `[+0.0077, +0.0278]` | `+0.1247` `[+0.1025, +0.1504]` | `True` |
| `vcc2025_h1` | `nourreddine_kolf21j` | `+0.0537` `[+0.0422, +0.0654]` | `+0.0890` `[+0.0744, +0.1096]` | `True` |
| `vcc2025_h1` | `replogle_k562_essential` | `+0.0281` `[+0.0189, +0.0368]` | `+0.1145` `[+0.0900, +0.1414]` | `True` |
| `vcc2025_h1` | `replogle_rpe1` | `-0.0134` `[-0.0269, -0.0003]` | `+0.1560` `[+0.1219, +0.1902]` | `True` |

## §5 Non-Confirmatory Descriptive Summary (`只报告、不判定`)
*Explicit notice (`PREREG_E §5.3`): Within `N_5`, `Transfer_Mean`, `Transfer_Mean_GlobalScale` (`eta=0`), and `Source_Grand_Mean` (`eta=0`) equal zero by exact algebraic construction, and `E1_ConstMed`, `E2_ConstP10/P90`, `E6_LogFCBar`, `E7_FCBarMinus1` are suppressed by cone construction — none constitute empirical validation evidence.*
- **`vcc2025_h1`**: Legacy unpaired `H-D2`=`True`; `FCT - E1 (all)`=`+0.0989` `[+0.0851, +0.1164]`; `N_5(all) FCT`=`+0.1088` (`eta=0.741`), `GS`=`+0.0000` (`eta=1.9e-15`), `E1`=`+0.0000` (`eta=0.0009`); `H-D_util` pass=`True` (`R^2(FCT orth/all)` approx via unprojected `r_cal^2` reported in JSON).
