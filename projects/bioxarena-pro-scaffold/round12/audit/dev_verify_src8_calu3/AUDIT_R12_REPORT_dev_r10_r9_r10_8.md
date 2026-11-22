# Round 12 Independent Pre-Registration Audit Report (`dev_r10`)

- **Source Mode:** `r9_r10_8` | **Seed:** `20261201` | **Bootstraps:** `400` (`97.5%` two-sided Bonferroni CIs) | **Perms:** `20`
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
| `calu3` (`9`) | `+0.0597` `[-0.0003, +0.1453]` (`1.000`, `+0.362`) | `False` | `6/8` (`False`) | `False` (`+0.0190`) | **`False`** | `+0.0895` `[-0.0630, +0.2503]` | `+0.0685` `[-0.1009, +0.2468]` (`9`) | **`False`** |

## `H-E1b` Source-by-Source Paired `FCT - Mis_s` (`all_nontarget`, `r_orth_cal`)
| Target | Source `s'` | `Mis_s` `[97.5% CI]` | Paired `FCT - Mis_s` `[97.5% CI]` | Pass |
|---|---|---|---|---|
| `calu3` | `nadig_hepg2` | `-0.0029` `[-0.0386, +0.0397]` | `+0.0625` `[+0.0200, +0.1211]` | `True` |
| `calu3` | `nadig_jurkat` | `-0.0441` `[-0.0782, +0.0098]` | `+0.1038` `[+0.0633, +0.1632]` | `True` |
| `calu3` | `nourreddine_kolf21j` | `-0.0238` `[-0.0476, +0.0168]` | `+0.0835` `[+0.0344, +0.1426]` | `True` |
| `calu3` | `replogle_k562_essential` | `-0.0424` `[-0.0641, +0.0047]` | `+0.1021` `[+0.0630, +0.1513]` | `True` |
| `calu3` | `replogle_rpe1` | `+0.0012` `[-0.0507, +0.0851]` | `+0.0585` `[+0.0252, +0.0979]` | `True` |
| `calu3` | `vcc2025_h1` | `-0.0040` `[-0.0356, +0.0392]` | `+0.0637` `[-0.0039, +0.1497]` | `False` |
| `calu3` | `dld1` | `-0.0082` `[-0.0473, +0.0655]` | `+0.0678` `[+0.0438, +0.0930]` | `True` |
| `calu3` | `cd4t_rest` | `+0.0286` `[-0.0120, +0.0905]` | `+0.0311` `[-0.0462, +0.1147]` | `False` |

## §5 Non-Confirmatory Descriptive Summary (`只报告、不判定`)
*Explicit notice (`PREREG_E §5.3`): Within `N_5`, `Transfer_Mean`, `Transfer_Mean_GlobalScale` (`eta=0`), and `Source_Grand_Mean` (`eta=0`) equal zero by exact algebraic construction, and `E1_ConstMed`, `E2_ConstP10/P90`, `E6_LogFCBar`, `E7_FCBarMinus1` are suppressed by cone construction — none constitute empirical validation evidence.*
- **`calu3`**: Legacy unpaired `H-D2`=`False`; `FCT - E1 (all)`=`+0.0387` `[-0.0336, +0.1286]`; `N_5(all) FCT`=`+0.0551` (`eta=0.866`), `GS`=`+0.0000` (`eta=1.7e-15`), `E1`=`-0.0000` (`eta=0.0024`); `H-D_util` pass=`True` (`R^2(FCT orth/all)` approx via unprojected `r_cal^2` reported in JSON).
