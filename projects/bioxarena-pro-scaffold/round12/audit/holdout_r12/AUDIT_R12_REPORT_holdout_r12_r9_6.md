# Round 12 Independent Pre-Registration Audit Report (`holdout_r12`)

- **Source Mode:** `r9_6` | **Seed:** `20261201` | **Bootstraps:** `4000` (`97.5%` two-sided Bonferroni CIs) | **Perms:** `100`
- **Adjudication Status:** `CONFIRMATORY_SINGLE_LINE_K1` (`K_conf = 1`, confirmatory lines = `['a549_ko']`)
- **Claim 1 (`H-E1`):** `CONFIRMED` | **Claim 2 (`H-E2`):** `CONFIRMED`

## Mandatory Pre-Registered Disclosures (`PREREG_E §0.1, §0.2, §7, §8`)
- **Post-Hoc Provenance (`E6`):** Both `H-E1` (paired `FCT - Mis_s`) and `H-E2` (`SA-8`) hypotheses — and the interpretation that prior T-cell unpaired `H-D2` failures reflected low unpaired power — were formulated after inspecting Rounds 9–11 data.
- **Claim 2 Empirical & Metric Scope (`E6`):** `SA-8` is strictly an empirical finding on residual metric `r_orth_cal` (`resp_top50_split`), NOT mechanistic proof (Round 11 Check 4 noiseless simulation predicted the wrong sign on `8/12` lines under `all_nontarget`), and does NOT extrapolate to Nicol et al. (2026) DEG metrics.
- **Claim 1 Novelty Scope (`§8`):** Must cite Kawsar et al. 2026 (`CanPert/CASTER`, `doi:10.64898/2026.11.15.766815`) for prior per-target paired matched-vs-mismatched control comparison; Claim 1 tests prospectively whether baseline-only `FCT` outperforms all `6` mismatched-line baselines on residual `r_orth_cal` in an unseen lineage and modality.
- **Scope & Modality:** MANDATORY SCOPE NOTICE (K_conf = 1, 推广性有限): Confirmatory adjudication relies on a single line (a549_ko, CRISPRko); hek293_xucao status = SECONDARY_ONLY. Cross-line generalizability is limited.
- **CRISPRi -> CRISPRko Modality Boundary:** STRENGTHENS_CROSS_LINE_AND_CROSS_MODALITY_GENERALIZATION (CRISPRi -> CRISPRko A549)

## Confirmatory Endpoints (`H-E1` & `H-E2`, `r_orth_cal`)
| Line (`N`) | `FCT` `[97.5% CI]` (`vf`, `rho`) | `H-E1a` | `H-E1b` Paired Pass | `H-E1c` (`> Perm97.5`) | `H-E1` | `E1(split)` `[97.5% CI]` | `E1(split)-E1(all)` `[97.5% CI]` (`N_c`) | `H-E2` |
|---|---|:---:|:---:|:---:|:---:|---|---|:---:|
| `a549_ko` (`445`) | `+0.1898` `[+0.1742, +0.2055]` (`1.000`, `+0.789`) | `True` | `6/6` (`True`) | `True` (`+0.0018`) | **`True`** | `+0.2565` `[+0.2195, +0.2981]` | `+0.2766` `[+0.2412, +0.3164]` (`445`) | **`True`** |
| `hek293_xucao` (`16`) | `+0.1561` `[+0.0595, +0.2593]` (`1.000`, `+0.301`) | `True` | `6/6` (`True`) | `True` (`+0.0015`) | **`True`** | `+0.0157` `[-0.0988, +0.1350]` | `+0.0339` `[-0.0868, +0.1797]` (`16`) | **`False`** |
| `a549_ko_cas9pos` (`445`) | `+0.1950` `[+0.1797, +0.2108]` (`1.000`, `+0.739`) | `True` | `6/6` (`True`) | `True` (`+0.0004`) | **`True`** | `+0.2523` `[+0.2172, +0.2918]` | `+0.2735` `[+0.2398, +0.3112]` (`445`) | **`True`** |

## `H-E1b` Source-by-Source Paired `FCT - Mis_s` (`all_nontarget`, `r_orth_cal`)
| Target | Source `s'` | `Mis_s` `[97.5% CI]` | Paired `FCT - Mis_s` `[97.5% CI]` | Pass |
|---|---|---|---|---|
| `a549_ko` | `nadig_hepg2` | `+0.0029` `[-0.0053, +0.0114]` | `+0.1869` `[+0.1693, +0.2049]` | `True` |
| `a549_ko` | `nadig_jurkat` | `+0.1102` `[+0.0977, +0.1237]` | `+0.0796` `[+0.0693, +0.0900]` | `True` |
| `a549_ko` | `nourreddine_kolf21j` | `+0.0829` `[+0.0734, +0.0927]` | `+0.1069` `[+0.0985, +0.1156]` | `True` |
| `a549_ko` | `replogle_k562_essential` | `+0.0619` `[+0.0519, +0.0723]` | `+0.1279` `[+0.1150, +0.1410]` | `True` |
| `a549_ko` | `replogle_rpe1` | `+0.1094` `[+0.0966, +0.1231]` | `+0.0804` `[+0.0682, +0.0917]` | `True` |
| `a549_ko` | `vcc2025_h1` | `+0.0393` `[+0.0303, +0.0483]` | `+0.1505` `[+0.1402, +0.1605]` | `True` |
| `hek293_xucao` | `nadig_hepg2` | `-0.0354` `[-0.0630, -0.0036]` | `+0.1914` `[+0.1013, +0.2945]` | `True` |
| `hek293_xucao` | `nadig_jurkat` | `-0.0354` `[-0.0796, +0.0117]` | `+0.1914` `[+0.1223, +0.2607]` | `True` |
| `hek293_xucao` | `nourreddine_kolf21j` | `+0.0450` `[-0.0139, +0.1083]` | `+0.1110` `[+0.0666, +0.1588]` | `True` |
| `hek293_xucao` | `replogle_k562_essential` | `-0.0050` `[-0.0532, +0.0457]` | `+0.1611` `[+0.1047, +0.2234]` | `True` |
| `hek293_xucao` | `replogle_rpe1` | `-0.0437` `[-0.0657, -0.0099]` | `+0.1997` `[+0.1139, +0.2923]` | `True` |
| `hek293_xucao` | `vcc2025_h1` | `-0.0237` `[-0.0769, +0.0199]` | `+0.1798` `[+0.1043, +0.2601]` | `True` |
| `a549_ko_cas9pos` | `nadig_hepg2` | `+0.0043` `[-0.0039, +0.0129]` | `+0.1907` `[+0.1737, +0.2096]` | `True` |
| `a549_ko_cas9pos` | `nadig_jurkat` | `+0.1102` `[+0.0976, +0.1225]` | `+0.0849` `[+0.0739, +0.0957]` | `True` |
| `a549_ko_cas9pos` | `nourreddine_kolf21j` | `+0.0838` `[+0.0745, +0.0933]` | `+0.1112` `[+0.1028, +0.1200]` | `True` |
| `a549_ko_cas9pos` | `replogle_k562_essential` | `+0.0620` `[+0.0523, +0.0722]` | `+0.1330` `[+0.1193, +0.1464]` | `True` |
| `a549_ko_cas9pos` | `replogle_rpe1` | `+0.1115` `[+0.0987, +0.1249]` | `+0.0835` `[+0.0719, +0.0959]` | `True` |
| `a549_ko_cas9pos` | `vcc2025_h1` | `+0.0387` `[+0.0298, +0.0476]` | `+0.1563` `[+0.1465, +0.1665]` | `True` |

## §5 Non-Confirmatory Descriptive Summary (`只报告、不判定`)
*Explicit notice (`PREREG_E §5.3`): Within `N_5`, `Transfer_Mean`, `Transfer_Mean_GlobalScale` (`eta=0`), and `Source_Grand_Mean` (`eta=0`) equal zero by exact algebraic construction, and `E1_ConstMed`, `E2_ConstP10/P90`, `E6_LogFCBar`, `E7_FCBarMinus1` are suppressed by cone construction — none constitute empirical validation evidence.*
- **`a549_ko`**: Legacy unpaired `H-D2`=`True`; `FCT - E1 (all)`=`+0.2099` `[+0.1975, +0.2223]`; `N_5(all) FCT`=`+0.1992` (`eta=0.935`), `GS`=`+0.0000` (`eta=1.8e-15`), `E1`=`+0.0001` (`eta=0.0028`); `H-D_util` pass=`True` (`R^2(FCT orth/all)` approx via unprojected `r_cal^2` reported in JSON).
- **`hek293_xucao`**: Legacy unpaired `H-D2`=`False`; `FCT - E1 (all)`=`+0.1743` `[+0.0853, +0.2653]`; `N_5(all) FCT`=`+0.1966` (`eta=0.852`), `GS`=`+0.0000` (`eta=1.7e-15`), `E1`=`+0.0000` (`eta=0.0016`); `H-D_util` pass=`True` (`R^2(FCT orth/all)` approx via unprojected `r_cal^2` reported in JSON).
- **`a549_ko_cas9pos`**: Legacy unpaired `H-D2`=`True`; `FCT - E1 (all)`=`+0.2163` `[+0.2045, +0.2283]`; `N_5(all) FCT`=`+0.2049` (`eta=0.934`), `GS`=`+0.0000` (`eta=1.8e-15`), `E1`=`+0.0001` (`eta=0.0028`); `H-D_util` pass=`True` (`R^2(FCT orth/all)` approx via unprojected `r_cal^2` reported in JSON).
