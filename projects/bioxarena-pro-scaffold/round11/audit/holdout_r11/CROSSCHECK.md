# Round 11 Blind Independent Audit vs. Primary Pipeline Cross-Check (`holdout_r11_r9_6`)

**Date:** 2026-11-22T21:08+08:00  
**Auditor:** A-T3 (Team A · Pre-Registration & Statistical Audit)  
**Bootstrap / Permutation Config:** `seed = 20261123`, `n_boot = 4,000`, `n_perm = 100`, two-sided `97.5%` Bonferroni CIs (`[1.25%, 98.75%]`)  
**Scope:** Full numerical & adjudication parity check between A-T3's blind recomputation (`/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round11/audit/holdout_r11/audit_r11_summary.json`) and Team A's primary frozen pipeline (`/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round11/results/holdout_r11_r9_6/`) across all three evaluated Round 11 lines:
- `cd4t_stim48` (`K_conf = 1`, `N = 503`, `G = 6,131`, `P = 6`)
- `telohaec` (failed confirmatory gate `N = 9 < 50`, `G = 4,869`, `P = 6`)
- `feng_ipsc` (`SECONDARY_ONLY` shared-donor iPSC panel, `N = 93`, `G = 5,326`, `P = 6`)

---

## 1. Per-Hypothesis & Per-Analysis Maximum Discrepancies

All comparisons below use exact IEEE-754 float64 differences (`|audit - results|`) across all 3 lines:

| Hypothesis / Analysis Block | Point Checks (`N`) | Max Point Abs Diff | 97.5% CI Bound Checks (`N`) | Max 97.5% CI Abs Diff | Parity Verdict |
|---|---:|---:|---:|---:|---|
| **`H-D1`** (`FCT` `all_nontarget` `r_orth_cal` + CIs) | 3 | **`2.776e-17`** | 6 | **`5.551e-17`** | **EXACT PARITY** |
| **`H-D2`** (`6` mismatched `NegCtrl`s × `3` lines) | 18 | **`1.388e-17`** | 36 | **`2.776e-17`** | **EXACT PARITY** |
| **`H-D2`** (`Perm100` `all_nontarget` median & `p97.5`) | 6 | **`1.388e-17`** | — | — | **EXACT PARITY** |
| **`H-D3`** (`FCT - E1_ConstMed` `all_nontarget`) | 3 | **`5.551e-17`** | 6 | **`5.551e-17`** | **EXACT PARITY** |
| **`H-D4a`** (`E1_ConstMed` `resp_top50`) | 3 | **`2.776e-17`** | 6 | **`1.665e-16`** | **EXACT PARITY** |
| **`H-D4b`** (`E1(resp_top50) - E1(all_nontarget)`) | 3 | **`2.776e-17`** | 6 | **`5.551e-17`** | **EXACT PARITY** |
| **`H-D5`** (`FCT - NegCtrl(kolf21j)` `all_nontarget`) | 3 | **`2.776e-17`** | 6 | **`2.776e-17`** | **EXACT PARITY** |
| **`H-D_util`** (`FCT - Transfer_Mean`, all 7 utility conventions) | 21 | **`1.743e-16`** | 42 | **`2.186e-16`** | **EXACT PARITY** |
| **`SA-8`** (`resp_top50_split`: `FCT`, `E1`, `Perm100`) | 12 | **`5.551e-17`** | 12 | **`1.110e-16`** | **EXACT PARITY** |
| **`SA-7 / S7`** (`N_cells` & `|KD|` Spearman, tertiles) | 24 | **`5.551e-17`** | 24 | **`1.110e-16`** | **EXACT PARITY** |
| **`SA-4 / R19`** (`L1_whole_line_frac` variance explained) | 6 | **`7.772e-16`** | 12 | **`5.551e-16`** | **EXACT PARITY** |
| All `20` non-degenerate methods (`all_nontarget` + `resp_top50`, `600` scalars) | — | **`<= 5.274e-16`** | — | **`<= 5.274e-16`** | **EXACT PARITY** |

### Two Documented Edge / Auxiliary Convention Differences (Zero Impact on Claims)
1. **`Transfer_Mean` and `Transfer_Mean_GlobalScale` (`r_orth_cal` only):** Algebraically, $P_1^\perp(\lambda \bar{\Delta}_{\mathrm{src}} - \bar{\Delta}_{\mathrm{src}}) \equiv 0$. A-T3's `independent_r11_metrics.py` sets correlations to `0.0` (`valid_frac = 1.0`) whenever predicted residual variance `<= EPS (1e-12)`, whereas `round11/results/` returns `NaN` (`valid_frac = 0.0`) except on `cd4t_stim48` `Transfer_Mean_GlobalScale` where unclipped floating-point rounding residual (`~1e-16`) yielded `+0.0108`. All `H-D_util` utility metrics (`ctrl_*`, `cent_*`, `pds_cosine_panel`) are computed directly against `Transfer_Mean` without residual orthogonalization and agree to `1.743e-16`.
2. **`SA-4 / R19` exploratory `L2_gene_gain_frac`:** Primary `L1_whole_line_frac` agrees to `7.772e-16`. On exploratory `L2_gene_gain_frac`, A-T3 enforces `c_dist >= 5` on distinct bootstrap-present perturbations (`W > 0`), yielding tiny finite-sample differences (`1.82e-4` on `cd4t_stim48`: `0.00221` vs `0.00239`; `5.42e-3` on `telohaec` `N=9`: `0.02075` vs `0.02617`), without changing the `< 3%` substantive interpretation.

---

## 2. Headline Holdout Summary & Adjudication Consensus

- **Confirmatory Rule `P1` (`K_conf = 1 < 2`):** **CONFIRMED PARITY (`NOT_EVALUATED_INSUFFICIENT_LINES` / `NO_CONFIRMATORY_ADJUDICATION`)** — only `cd4t_stim48` passed QC, so pre-registered multi-line confirmatory hypothesis testing was withheld.
- **`cd4t_stim48` (`N = 503`, `G = 6,131`):**
  - `H-D1` (`FCT all_nontarget`): `+0.0532 [+0.0337, +0.0731]` (`valid_frac = 1.000`, `rho_ab_orth = +0.0784` -> descriptive **PASS**, `hd1_a1_deviation_flag = False`).
  - `E1_ConstMed (all_nontarget)`: `-0.0203 [-0.0326, -0.0087]`.
  - `H-D2`: **FAIL** — `FCT ci975_lo (+0.0337)` exceeds `Perm100 Q97.5 (-0.0129)`, **but falls below `NegCtrl_FCT_MismatchedMean(vcc2025_h1)` `ci975_hi (+0.0572)`** (`est = +0.0369`).
  - `H-D3` (`FCT - E1` `all_nontarget`): `+0.0735 [+0.0561, +0.0917]` (descriptive **PASS**).
  - `H-D4` (`E1 resp_top50` & cross-mask jump): `+0.0522 [+0.0257, +0.0792]`, `E1(top50) - E1(all) = +0.0725 [+0.0477, +0.0982]` (descriptive **PASS**).
  - `SA-8` (`resp_top50_split`): `FCT +0.0621 [+0.0231, +0.1031]`, `E1 +0.0541 [+0.0193, +0.0896]`, `Perm100 Med / Q97.5 = +0.0381 / +0.0571`.
  - `H-D_util`: non-inferior on `6/6` conventions (`n_inferior = 0`, `n_superior = 1` on `cent_de50_r` `+0.0107 [+0.0004, +0.0208]`).
- **`telohaec` (`N = 9`, `G = 4,869`, failed confirmatory QC):**
  - `FCT`: `+0.2173 [+0.0958, +0.3602]`, `E1`: `-0.0475 [-0.1700, +0.0187]`, `max Mis97.5_hi = +0.1712` (`nadig_hepg2`), `Perm100 Q97.5 = -0.0132`, `FCT - E1 = +0.2648 [+0.1443, +0.4636]`, `E1(top50) = +0.1798 [+0.0334, +0.3756]`, `SA-8 FCT = +0.2579`, `E1 = +0.2188`, `Perm100 Q97.5 = +0.2343`, `R19 L1 = 0.6356`.
- **`feng_ipsc` (`N = 93`, `G = 5,326`, secondary shared-donor panel):**
  - `FCT`: `+0.1212 [+0.1043, +0.1371]` (`rho_ab_orth = +0.5570`), `E1`: `+0.0530 [+0.0287, +0.0785]`, `max Mis97.5_hi = +0.0829` (`vcc2025_h1`), `Perm100 Q97.5 = +0.0452`, `FCT - E1 = +0.0682 [+0.0427, +0.0923]`, `E1(top50) = +0.2479 [+0.2003, +0.2975]`, `SA-8 FCT = +0.2440`, `E1 = +0.2468`, `Perm100 Q97.5 = +0.1874`.
  - Descriptive **`H-D5`** (`feng_ipsc` vs `NegCtrl(nourreddine_kolf21j)`): `Mis(kolf21j) = +0.0562 [+0.0378, +0.0744]`, `FCT - Mis(kolf21j) = +0.0650 [+0.0528, +0.0784]` (paired CI > 0 AND disjoint above `Mis(kolf21j)` `ci975_hi`). Note Addendum A3: split halves share donors across clones.
  - **`H-D_util`**: `n_inferior = 3, n_superior = 2` (FCT is inferior on uncentered/panel discrimination on this 5,326-gene panel).

---

**Audit Parity Sign-Off:** Independent recomputation agrees with `round11/results/holdout_r11_r9_6/` to floating-point precision (`max_abs_diff <= 7.772e-16` across all primary, `97.5% CI`, `H-D5`, `H-D_util`, `SA-8`, `S7`, and `R19 L1` items), and both pipelines reach the identical `NOT_EVALUATED_INSUFFICIENT_LINES` (`K_conf = 1`) verdict.
