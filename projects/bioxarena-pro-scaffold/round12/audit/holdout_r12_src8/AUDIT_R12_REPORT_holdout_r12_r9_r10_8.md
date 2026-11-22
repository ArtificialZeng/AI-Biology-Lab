# Round 12 Independent Pre-Registration Audit Report (`holdout_r12`)

- **Source Mode:** `r9_r10_8` | **Seed:** `20261201` | **Bootstraps:** `4000` (`97.5%` two-sided Bonferroni CIs) | **Perms:** `100`
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
| `a549_ko` (`528`) | `+0.1563` `[+0.1414, +0.1721]` (`1.000`, `+0.778`) | `True` | `8/8` (`True`) | `True` (`+0.0023`) | **`True`** | `+0.2302` `[+0.1955, +0.2677]` | `+0.2462` `[+0.2132, +0.2812]` (`528`) | **`True`** |
| `hek293_xucao` (`16`) | `+0.1544` `[+0.0600, +0.2509]` (`1.000`, `+0.299`) | `True` | `8/8` (`True`) | `True` (`+0.0002`) | **`True`** | `+0.0115` `[-0.0960, +0.1468]` | `+0.0305` `[-0.0798, +0.1882]` (`16`) | **`False`** |
| `a549_ko_cas9pos` (`526`) | `+0.1616` `[+0.1464, +0.1773]` (`1.000`, `+0.727`) | `True` | `8/8` (`True`) | `True` (`+0.0016`) | **`True`** | `+0.2291` `[+0.1961, +0.2643]` | `+0.2461` `[+0.2149, +0.2782]` (`526`) | **`True`** |

## `H-E1b` Source-by-Source Paired `FCT - Mis_s` (`all_nontarget`, `r_orth_cal`)
| Target | Source `s'` | `Mis_s` `[97.5% CI]` | Paired `FCT - Mis_s` `[97.5% CI]` | Pass |
|---|---|---|---|---|
| `a549_ko` | `nadig_hepg2` | `+0.0119` `[+0.0054, +0.0183]` | `+0.1445` `[+0.1283, +0.1611]` | `True` |
| `a549_ko` | `nadig_jurkat` | `+0.0965` `[+0.0840, +0.1088]` | `+0.0598` `[+0.0506, +0.0695]` | `True` |
| `a549_ko` | `nourreddine_kolf21j` | `+0.0748` `[+0.0657, +0.0840]` | `+0.0816` `[+0.0734, +0.0901]` | `True` |
| `a549_ko` | `replogle_k562_essential` | `+0.0567` `[+0.0471, +0.0662]` | `+0.0996` `[+0.0873, +0.1122]` | `True` |
| `a549_ko` | `replogle_rpe1` | `+0.0819` `[+0.0715, +0.0926]` | `+0.0745` `[+0.0657, +0.0839]` | `True` |
| `a549_ko` | `vcc2025_h1` | `+0.0388` `[+0.0303, +0.0476]` | `+0.1176` `[+0.1068, +0.1290]` | `True` |
| `a549_ko` | `dld1` | `+0.1140` `[+0.1013, +0.1265]` | `+0.0423` `[+0.0372, +0.0480]` | `True` |
| `a549_ko` | `cd4t_rest` | `+0.0112` `[+0.0026, +0.0201]` | `+0.1451` `[+0.1332, +0.1575]` | `True` |
| `hek293_xucao` | `nadig_hepg2` | `-0.0239` `[-0.0533, +0.0157]` | `+0.1783` `[+0.0914, +0.2733]` | `True` |
| `hek293_xucao` | `nadig_jurkat` | `-0.0250` `[-0.0633, +0.0185]` | `+0.1794` `[+0.1122, +0.2471]` | `True` |
| `hek293_xucao` | `nourreddine_kolf21j` | `+0.0484` `[-0.0090, +0.1079]` | `+0.1060` `[+0.0629, +0.1524]` | `True` |
| `hek293_xucao` | `replogle_k562_essential` | `+0.0043` `[-0.0431, +0.0540]` | `+0.1501` `[+0.0980, +0.2075]` | `True` |
| `hek293_xucao` | `replogle_rpe1` | `-0.0353` `[-0.0624, +0.0081]` | `+0.1898` `[+0.1052, +0.2814]` | `True` |
| `hek293_xucao` | `vcc2025_h1` | `-0.0222` `[-0.0682, +0.0173]` | `+0.1767` `[+0.1016, +0.2549]` | `True` |
| `hek293_xucao` | `dld1` | `+0.0114` `[-0.0296, +0.0558]` | `+0.1430` `[+0.0821, +0.2091]` | `True` |
| `hek293_xucao` | `cd4t_rest` | `-0.0239` `[-0.0492, +0.0008]` | `+0.1784` `[+0.0879, +0.2706]` | `True` |
| `a549_ko_cas9pos` | `nadig_hepg2` | `+0.0128` `[+0.0061, +0.0196]` | `+0.1489` `[+0.1323, +0.1658]` | `True` |
| `a549_ko_cas9pos` | `nadig_jurkat` | `+0.0971` `[+0.0855, +0.1090]` | `+0.0645` `[+0.0550, +0.0741]` | `True` |
| `a549_ko_cas9pos` | `nourreddine_kolf21j` | `+0.0760` `[+0.0673, +0.0850]` | `+0.0856` `[+0.0766, +0.0946]` | `True` |
| `a549_ko_cas9pos` | `replogle_k562_essential` | `+0.0572` `[+0.0483, +0.0663]` | `+0.1045` `[+0.0915, +0.1168]` | `True` |
| `a549_ko_cas9pos` | `replogle_rpe1` | `+0.0840` `[+0.0740, +0.0949]` | `+0.0777` `[+0.0680, +0.0874]` | `True` |
| `a549_ko_cas9pos` | `vcc2025_h1` | `+0.0382` `[+0.0301, +0.0468]` | `+0.1234` `[+0.1119, +0.1349]` | `True` |
| `a549_ko_cas9pos` | `dld1` | `+0.1160` `[+0.1044, +0.1282]` | `+0.0457` `[+0.0399, +0.0515]` | `True` |
| `a549_ko_cas9pos` | `cd4t_rest` | `+0.0113` `[+0.0027, +0.0199]` | `+0.1503` `[+0.1374, +0.1629]` | `True` |

## §5 Non-Confirmatory Descriptive Summary (`只报告、不判定`)
*Explicit notice (`PREREG_E §5.3`): Within `N_5`, `Transfer_Mean`, `Transfer_Mean_GlobalScale` (`eta=0`), and `Source_Grand_Mean` (`eta=0`) equal zero by exact algebraic construction, and `E1_ConstMed`, `E2_ConstP10/P90`, `E6_LogFCBar`, `E7_FCBarMinus1` are suppressed by cone construction — none constitute empirical validation evidence.*
- **`a549_ko`**: Legacy unpaired `H-D2`=`True`; `FCT - E1 (all)`=`+0.1723` `[+0.1591, +0.1858]`; `N_5(all) FCT`=`+0.1556` (`eta=0.854`), `GS`=`+0.0000` (`eta=1.8e-15`), `E1`=`+0.0000` (`eta=0.0037`); `H-D_util` pass=`True` (`R^2(FCT orth/all)` approx via unprojected `r_cal^2` reported in JSON).
- **`hek293_xucao`**: Legacy unpaired `H-D2`=`False`; `FCT - E1 (all)`=`+0.1734` `[+0.0826, +0.2641]`; `N_5(all) FCT`=`+0.1962` (`eta=0.850`), `GS`=`+0.0000` (`eta=1.4e-15`), `E1`=`+0.0001` (`eta=0.0028`); `H-D_util` pass=`True` (`R^2(FCT orth/all)` approx via unprojected `r_cal^2` reported in JSON).
- **`a549_ko_cas9pos`**: Legacy unpaired `H-D2`=`True`; `FCT - E1 (all)`=`+0.1787` `[+0.1647, +0.1924]`; `N_5(all) FCT`=`+0.1614` (`eta=0.855`), `GS`=`+0.0000` (`eta=1.8e-15`), `E1`=`+0.0001` (`eta=0.0037`); `H-D_util` pass=`True` (`R^2(FCT orth/all)` approx via unprojected `r_cal^2` reported in JSON).
