# Round 12 Independent Pre-Registration Audit Report (`dev_r10`)

- **Source Mode:** `r9_6` | **Seed:** `20261201` | **Bootstraps:** `4000` (`97.5%` two-sided Bonferroni CIs) | **Perms:** `100`
- **Adjudication Status:** `DEV_DRY_RUN_DEV_R10` (`K_conf = 0`, confirmatory lines = `[]`)
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
| `cd4t_rest` (`302`) | `+0.0556` `[+0.0383, +0.0741]` (`1.000`, `+0.147`) | `True` | `6/6` (`True`) | `True` (`+0.0131`) | **`True`** | `+0.0532` `[+0.0087, +0.1005]` | `+0.0383` `[-0.0067, +0.0863]` (`284`) | **`False`** |

## `H-E1b` Source-by-Source Paired `FCT - Mis_s` (`all_nontarget`, `r_orth_cal`)
| Target | Source `s'` | `Mis_s` `[97.5% CI]` | Paired `FCT - Mis_s` `[97.5% CI]` | Pass |
|---|---|---|---|---|
| `cd4t_rest` | `nadig_hepg2` | `+0.0070` `[-0.0061, +0.0214]` | `+0.0486` `[+0.0259, +0.0722]` | `True` |
| `cd4t_rest` | `nadig_jurkat` | `+0.0163` `[-0.0042, +0.0376]` | `+0.0392` `[+0.0249, +0.0527]` | `True` |
| `cd4t_rest` | `nourreddine_kolf21j` | `+0.0017` `[-0.0095, +0.0137]` | `+0.0538` `[+0.0375, +0.0706]` | `True` |
| `cd4t_rest` | `replogle_k562_essential` | `+0.0042` `[-0.0115, +0.0204]` | `+0.0514` `[+0.0357, +0.0679]` | `True` |
| `cd4t_rest` | `replogle_rpe1` | `+0.0158` `[+0.0001, +0.0335]` | `+0.0398` `[+0.0192, +0.0596]` | `True` |
| `cd4t_rest` | `vcc2025_h1` | `+0.0208` `[+0.0024, +0.0400]` | `+0.0347` `[+0.0260, +0.0431]` | `True` |

## §5 Non-Confirmatory Descriptive Summary (`只报告、不判定`)
*Explicit notice (`PREREG_E §5.3`): Within `N_5`, `Transfer_Mean`, `Transfer_Mean_GlobalScale` (`eta=0`), and `Source_Grand_Mean` (`eta=0`) equal zero by exact algebraic construction, and `E1_ConstMed`, `E2_ConstP10/P90`, `E6_LogFCBar`, `E7_FCBarMinus1` are suppressed by cone construction — none constitute empirical validation evidence.*
- **`cd4t_rest`**: Legacy unpaired `H-D2`=`False`; `FCT - E1 (all)`=`+0.0426` `[+0.0244, +0.0608]`; `N_5(all) FCT`=`+0.0553` (`eta=0.745`), `GS`=`+0.0000` (`eta=2.8e-15`), `E1`=`+0.0000` (`eta=0.0012`); `H-D_util` pass=`True` (`R^2(FCT orth/all)` approx via unprojected `r_cal^2` reported in JSON).
