# Independent Exploratory Audit Report (`EXPLORATORY_REVIEW_R11.md` §B & §C)

**Date:** 2026-11-22T21:42+08:00  
**Auditor:** A-T3 (Team A · Pre-Registration & Statistical Audit)  
**Environment & Parameters:** `python3 -I`, `seed = 20261124`, `n_boot = 4,000`, `n_perm = 100`, two-sided Bonferroni `97.5%` CIs (`[1.25%, 98.75%]`)  
**Scope & Isolation:** Implemented cleanly under `/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round11/audit/exploratory_check/` reading ONLY standardized caches (`round9/p2/shared/lines/`, `round10/shared/lines/`, `round11/shared/lines/`) and frozen `round11/audit/` modules. Zero `.py` or raw result `.json`/`.csv` files under `exploratory/null_cone_theorem/` or `round11/src/` were read.

---

## Executive Verdict Across the 4 Independent Checks

| Audit Check | Outcome | Key Audit Findings |
|---|---|---|
| **1. `SA-8` (`12` Lines, Frozen Spec)** | **Substantively CONFIRMED (`9/12` strict CI > 0 vs Section 2's `10/12`), `N` Discrepancy Exposed** | Locked `P5` target-self-gene filter reduces `N` on `11/12` lines (`cd4t_stim48` `656 -> 503`, `cd4t_rest` `398 -> 302`, `feng_ipsc` `111 -> 93`). Joint `E1(split) - E1(all)` stays `+` on `12/12` lines and `97.5% CI > 0` on **`9/12`** (`cd4t_rest` becomes non-significant: `+0.0383 [-0.0061, +0.0864]`). |
| **2. `N_5` Independent Audit & Control Validity** | **`GlobalScale`/`SGM`/`E1..E7` Non-Evidentiary (By Construction); Real Controls Mixed-to-Pass** | `GlobalScale` & `Source_Grand_Mean` have exact $\eta < 10^{-14}$ ($R^2_{\mathcal{N}_5}=1$); `E1/E2/E6/E7` are cone-suppressed by construction ($\eta \le 0.085$, $R^2_{\mathcal{N}_5} \ge 0.9928$). True informative tests show: single-draw `Perm0..99` $\approx 0$ (`12/12` `FCT > Perm97.5`), **`7/12`** pass unpaired `MaxMis_hi` (`cd4t_rest` `-0.0009` & `cd4t_stim48` `-0.0188` BOTH fail unpaired vs `Mis(vcc2025_h1)`), and **`11/12`** pass paired joint bootstrap (`cd4t_rest` `+0.0345 [+0.0261, +0.0432]`, `cd4t_stim48` `+0.0141 [+0.0064, +0.0218]`). |
| **3. `Q5a` Blocked & `LOO-line` Test + `Bonferroni x3`** | **Naive `p=2.33e-10` Overstated; **`expr05` `FCT > log1p` SURVIVES (`Bonf3 p = 0.0059`), `all` FAILS, `Cal logFC > log1p` Borderline (`0.0410` on `expr05`)** | The 33 pairs share 10 active lines (`cd4t_stim48` in 9 edges; max vertex-disjoint matching = **`5`**). After replacing naive Wilcoxon with exact $2^{10}$ line-blocked / target-block sign-flips, binomial sign tests, dyadic CGM cluster t-tests, LOO-line drops & `Bonferroni x3`: `Cal FCT > log1p` robustly survives ONLY on `expr05` (`10/10` lines, `Bonf3 p = 0.0059`, CGM `0.0002`), is mixed on `expr01`, and fails non-parametric tests on `all` (`Bonf3 p = 1.000`). `Cal logFC > log1p` fails `all`/`expr01` line-perm & reaches `Bonf3 p = 0.0410` (`expr05`). |
| **4. Noiseless Real-`mu` Theorem 1 Simulation** | **`SA-8` Tail Prediction CONFIRMED (`12/12`); `all_nontarget` FALSIFIED (`4/12`, `+0.6084` vs real `-0.0203` on `cd4t_stim48`)** | Using real $(\mu_t, \mu_s, \delta_p)$ at zero noise, `E1(split)` is `+` (`+0.564` to `+0.869`) on `12/12` lines (`12/12` match). However noiseless `E1(all)` ALSO predicts `+` on `12/12` lines (`+0.6084` on `cd4t_stim48`, larger than toy `+0.227`), missing the real negative `E1(all)` on **`8/12`** lines (`cd4t_stim48 -0.0203`, `rpe1 -0.0934`) caused by real finite-sample source noise + cross-line rewiring (`mid_mus < 0.5` carries `65%` of `E1` variance & negative covariance `-96.2` on `cd4t_stim48`). |

---

## Check 1: `SA-8` Across `12` Cell Lines (`seed = 20261124`, `n_boot = 4,000`, `n_perm = 100`)

Computed via frozen `round11/audit/` `r_orth_cal` (`evaluate_joint_residuals`, `evaluate_sa8_resp_top50_split`, `evaluate_sa8_perm_null`) plus joint cross-mask bootstrap `evaluate_cross_mask_split_minus_all` (`sa8_12lines.csv`, `sa8_12lines.json`):

| Line (`Grp`) | `N` (Us / Sec.2) | `E1(all)` | `E1(top50)` | `E1(SA-8 split)` `[97.5% CI]` | `E1(split) − E1(all)` `[97.5% CI]` | Sec.2 `dE1` | `FCT(split)` `[97.5% CI]` | `Perm100(split) Med / Q97.5` |
|---|---:|---:|---:|---|---|---:|---|---|
| `nadig_hepg2` (R9) | `797 / 985` | `-0.0470` | `+0.2259` | **`+0.1977`** `[+0.1768, +0.2193]` | **`+0.2447`** `[+0.2241, +0.2660]` | `+0.2689` | `+0.2589` `[+0.2393, +0.2795]` | `+0.1109 / +0.1553` |
| `nadig_jurkat` (R9) | `931 / 1180` | `-0.0782` | `+0.1209` | **`+0.0751`** `[+0.0549, +0.0949]` | **`+0.1532`** `[+0.1356, +0.1708]` | `+0.1449` | `+0.0879` `[+0.0682, +0.1082]` | `+0.0403 / +0.0963` |
| `nourreddine_kolf21j` (R9) | `18 / 20` | `-0.0172` | `+0.0962` | `+0.1039` `[-0.0496, +0.2548]` | `+0.1211` `[-0.0158, +0.2564]` *(NS)* | `+0.1945*` | `+0.2133` `[+0.0556, +0.3732]` | `+0.0689 / +0.1593` |
| `replogle_k562_essential` (R9) | `858 / 1044` | `-0.0377` | `+0.1900` | **`+0.1758`** `[+0.1521, +0.2002]` | **`+0.2135`** `[+0.1902, +0.2368]` | `+0.2172` | `+0.2005` `[+0.1797, +0.2224]` | `+0.1045 / +0.1492` |
| `replogle_rpe1` (R9) | `1068 / 1300` | `-0.0934` | `+0.1415` | **`+0.1217`** `[+0.1062, +0.1373]` | **`+0.2151`** `[+0.1998, +0.2306]` | `+0.2298` | `+0.1358` `[+0.1179, +0.1536]` | `+0.0625 / +0.0964` |
| `vcc2025_h1` (R9) | `80 / 91` | `+0.0437` | `+0.3076` | **`+0.3100`** `[+0.2455, +0.3711]` | **`+0.2663`** `[+0.2097, +0.3202]` | `+0.2417` | `+0.3533` `[+0.2793, +0.4202]` | `+0.2130 / +0.2472` |
| `dld1` (R10) | `86 / 93` | `-0.0614` | `+0.1567` | **`+0.1294`** `[+0.0903, +0.1706]` | **`+0.1907`** `[+0.1500, +0.2356]` | `+0.1610` | `+0.1407` `[+0.0794, +0.2065]` | `+0.0689 / +0.1008` |
| `cd4t_rest` (R10) | `302 / 398` | `+0.0130` | `+0.0495` | **`+0.0532`** `[+0.0097, +0.1010]` | `+0.0383` `[-0.0061, +0.0864]` *(NS; n_c=284)* | `+0.0936*` | `+0.0648` `[+0.0134, +0.1163]` | `+0.0384 / +0.0633` |
| `calu3` (R10) | `8 / 9` | `+0.0338` | `+0.1321` | `+0.1385` `[-0.0350, +0.2988]` | `+0.1048` `[-0.0658, +0.2672]` *(NS)* | `+0.0818` | `+0.0925` `[-0.0685, +0.2418]` | `+0.0786 / +0.1622` |
| `cd4t_stim48` (R11) | `503 / 656` | `-0.0203` | `+0.0522` | **`+0.0541`** `[+0.0198, +0.0941]` | **`+0.0730`** `[+0.0384, +0.1098]` *(n_c=422)* | `+0.0805` | `+0.0621` `[+0.0216, +0.1058]` | `+0.0381 / +0.0571` |
| `telohaec` (R11) | `9 / 9` | `-0.0475` | `+0.1798` | `+0.2188` `[-0.0338, +0.4615]` | **`+0.2663`** `[+0.0196, +0.5715]` | `+0.2550` *(NS)* | `+0.2579` `[-0.0680, +0.6051]` | `+0.1070 / +0.2343` |
| `feng_ipsc` (R11) | `93 / 111` | `+0.0530` | `+0.2479` | **`+0.2468`** `[+0.1991, +0.2947]` | **`+0.1938`** `[+0.1480, +0.2405]` | `+0.1808` | `+0.2440` `[+0.1906, +0.2956]` | `+0.1482 / +0.1874` |

### Key Cell-by-Cell Reconciliation Findings
1. **Cohort `N` discrepancy on `11/12` lines:** Locked `PREREG_D` / `round11/audit/` enforces `P5` target self-gene presence inside `G_eval(t)` (`target_indices >= 0`), giving `cd4t_stim48 N=503`, `cd4t_rest N=302`, `dld1 N=86`, `feng_ipsc N=93`. `RESULTS_13LINE_NULL_CONE_AUDIT.md` Section 2 & 3 omitted that self-gene filter (`656`, `398`, `93`, `111`), creating cell differences up to `0.0734` (`nourreddine_kolf21j dE1` `+0.1211` vs `+0.1945`) and `0.0553` (`cd4t_rest dE1` `+0.0383` vs `+0.0936`).
2. **Exact alignment with locked Round 11:** Our `SA-8` point estimates (`E1(all)`, `E1(top50)`, `E1(split)`, `FCT(split)`, `Perm100`) match `round11/audit/holdout_r11/` and `round11/results/holdout_r11_r9_6/` identically on all Round 11 lines (`cd4t_stim48`, `telohaec`, `feng_ipsc`), plus `dld1` & `vcc2025_h1` (`A-T2.txt`).
3. **Pass-count shift:** Under locked cohorts, joint `E1(split) - E1(all)` passes `97.5% CI > 0` on **`9/12`** lines (`cd4t_rest` fails `[-0.0061, +0.0864]` on `n_common=284` split-valid perturbations, and `nourreddine_kolf21j` fails `[-0.0158, +0.2564]`; `telohaec` passes `[+0.0196, +0.5715]`), vs `10/12` in Section 2. Note also `cd4t_rest` (`302 -> 284`) and `cd4t_stim48` (`503 -> 422`) have `half_valid=False` rows so joint cross-mask evaluates common perturbations (`all_est = +0.0149` and `-0.0190`).

---

## Check 2: Independent `N_5` Implementation (`r_cone_sp`) & Basis Validity Audit

Independently implemented $\mathcal{N}_5(p) = \mathrm{span}\{\mathbf{1}, \Delta_{\mathrm{TMean},p}, b_S, g_{0.05}(\overline{\mathrm{FC}}_p), g_{1.0}(\overline{\mathrm{FC}}_p), g_{20.0}(\overline{\mathrm{FC}}_p)\}$ via SVD (`rank = 6` on all perturbations across 12 lines) and semi-partial `r_cone_sp` (`n5_audit.json`, `n5_methods_12lines.csv`, `n5_perm100_12lines.csv`, `n5_spotlight_vs_vcc2025_h1.csv`).

### 2A. Structural Audit: Controls Zeroed/Suppressed BY CONSTRUCTION (Non-Evidentiary) vs. Informative Controls

Across all `12` lines (`all_nontarget`):

| Predictor / Control | `r_cone_sp` Range (`12` lines) | Residual Fraction $\eta$ Range | Variance in Cone $R^2_{\mathcal{N}_5}$ | Audit Classification (Can it count as empirical validation?) |
|---|---|---|---:|---|
| `Transfer_Mean_GlobalScale` | `+0.0000` .. `+0.0000` | `8.92e-16` .. `7.30e-15` | `>= 1.000000` | **NO — ZERO BY EXACT CONSTRUCTION** ($(\lambda-1)\Delta_{\mathrm{TMean},p} \in \mathcal{N}_5$) |
| `Source_Grand_Mean` ($b_S$) | `+0.0000` .. `+0.0000` | `8.36e-16` .. `5.54e-15` | `>= 1.000000` | **NO — ZERO BY EXACT CONSTRUCTION** ($b_S - \Delta_{\mathrm{TMean},p} \in \mathcal{N}_5$) |
| `E1_ConstMed` ($c_{\mathrm{med}}$) | `-0.0001` .. `+0.0000` | `8.79e-04` .. `3.88e-03` | `>= 0.999985` | **NO — SUPPRESSED BY CONE CONSTRUCTION** ($g_c \in \mathrm{conv}\{g_{0.05},g_{1.0},g_{20.0}\}$) |
| `E2_ConstP10` | `-0.0000` .. `+0.0000` | `2.90e-04` .. `1.73e-03` | `>= 0.999997` | **NO — SUPPRESSED BY CONE CONSTRUCTION** |
| `E2_ConstP90` | `-0.0000` .. `+0.0001` | `1.40e-03` .. `3.05e-03` | `>= 0.999991` | **NO — SUPPRESSED BY CONE CONSTRUCTION** |
| `E6_LogFCBar` ($c \to \infty$) | `-0.0001` .. `+0.0001` | `1.37e-03` .. `8.49e-02` | `>= 0.992784` | **NO — SUPPRESSED BY CONE BOUNDARY CONSTRUCTION** |
| `E7_FCBarMinus1` ($c \to 0^+$) | `-0.0001` .. `+0.0001` | `1.89e-03` .. `1.06e-02` | `>= 0.999888` | **NO — SUPPRESSED BY CONE BOUNDARY CONSTRUCTION** |
| **`Single-draw Perm0..99`** | `med -0.0010 .. +0.0007` | `med 0.631 .. 0.713` | `~0.49 .. 0.60` | **YES — INFORMATIVE NULL** ($\sim 68\%$ non-cone energy; non-expected draws) |
| **`Mismatched FCT (Mis_s)`** | `-0.0344` .. `+0.0955` | `~0.68 .. 0.95` | non-cone | **YES — INFORMATIVE TARGET-BASELINE TEST** |
| `E3_SrcMeanCtrl` / `E5_RankMatchedMu` | `+0.0054..+0.0544` / `+0.0066..+0.0637` | `0.754..0.836` / `0.732..0.943` | `~0.11 .. 0.46` | **YES — REFERENCE BASELINES** |
| `FCT` | `+0.0553` .. `+0.1899` | `0.741` .. `0.934` | `0.128 .. 0.451` | **PRIMARY MODEL UNDER TEST** |
| `Hybrid_ProjShared` / `TMVS` | `-0.0184..+0.2400` / `+0.0220..+0.1535` | `0.984..0.996` / `0.425..0.924` | `0.008..0.032` | **YES — NON-CONE TARGET MODELS** (`cd4t_stim48` Hybrid = `+0.2199 [+0.1497,+0.2867]`) |

### 2B. Single-Draw `Perm0..99` & All `6` Mismatched Baselines `Mis_s` Across `12` Lines

- **`Perm0..99` single-draw empirical distributions (`non-expected`):**
  Every single permutation draw retains high non-cone residual energy ($\mathrm{med}\,\eta \in [0.6308, 0.7133]$), yet random gene-baseline permutations center tightly around `0` across all `12` lines (`median |r_cone_sp|` $\in [0.0011, 0.0058]$, `Q97.5` $\in [+0.0032, +0.0156]$):
  - `nadig_hepg2`: `min=-0.0089, p025=-0.0078, med=-0.0007, p975=+0.0075, max=+0.0085 (med_eta=0.6813)`
  - `nadig_jurkat`: `min=-0.0097, p025=-0.0084, med=+0.0007, p975=+0.0083, max=+0.0109 (med_eta=0.6829)`
  - `nourreddine_kolf21j`: `min=-0.0108, p025=-0.0091, med=+0.0001, p975=+0.0103, max=+0.0155 (med_eta=0.6804)`
  - `replogle_k562_essential`: `min=-0.0095, p025=-0.0071, med=+0.0000, p975=+0.0075, max=+0.0105 (med_eta=0.6857)`
  - `replogle_rpe1`: `min=-0.0097, p025=-0.0087, med=-0.0010, p975=+0.0070, max=+0.0078 (med_eta=0.7037)`
  - `vcc2025_h1`: `min=-0.0084, p025=-0.0044, med=-0.0003, p975=+0.0043, max=+0.0070 (med_eta=0.6470)`
  - `dld1`: `min=-0.0098, p025=-0.0066, med=-0.0001, p975=+0.0072, max=+0.0090 (med_eta=0.7023)`
  - `cd4t_rest`: `min=-0.0045, p025=-0.0032, med=+0.0001, p975=+0.0032, max=+0.0040 (med_eta=0.6317)`
  - `calu3`: `min=-0.0137, p025=-0.0110, med=+0.0000, p975=+0.0102, max=+0.0135 (med_eta=0.7133)`
  - `cd4t_stim48`: `min=-0.0055, p025=-0.0042, med=-0.0001, p975=+0.0037, max=+0.0046 (med_eta=0.6372)`
  - `telohaec`: `min=-0.0197, p025=-0.0144, med=-0.0003, p975=+0.0156, max=+0.0170 (med_eta=0.6308)`
  - `feng_ipsc`: `min=-0.0104, p025=-0.0082, med=+0.0006, p975=+0.0068, max=+0.0092 (med_eta=0.6880)`
  - `FCT ci975_lo > Perm100 Q97.5` holds on **`12/12`** lines.

- **All `6` Mismatched Baselines `Mis_s` vs `FCT` (`N_5`):**
  Full point/CI and paired `FCT - Mis_s` contrasts are stored in `n5_audit.json`. Under locked `P5` cohorts:
  - Unpaired strict criterion `FCT ci975_lo > max_s Mis_s ci975_hi`: passes **`7/12`** lines (`nadig_hepg2` `+0.1829 > +0.0599`, `nadig_jurkat` `+0.1266 > +0.1007`, `replogle_k562_essential` `+0.1796 > +0.0685`, `replogle_rpe1` `+0.1563 > +0.0840`, `vcc2025_h1` `+0.0907 > +0.0444`, `dld1` `+0.1289 > +0.0704`, `feng_ipsc` `+0.1161 > +0.0800`); **FAILS `5/12`** (`nourreddine_kolf21j` `+0.0440 <= +0.0864`, `cd4t_rest` `+0.0351 <= +0.0360`, `calu3` `+0.0244 <= +0.0806`, `cd4t_stim48` `+0.0559 <= +0.0747`, `telohaec` `+0.0811 <= +0.1589`).
  - Paired joint bootstrap `FCT - Mis_s ci975_lo > 0` across **all** sources: passes **`11/12`** lines (fails only `nourreddine_kolf21j - Mis(k562)` `+0.0525 [-0.0004, +0.0955]`).

### 2C. Spotlight Audit: `cd4t_rest` & `cd4t_stim48` Relative to `Mis(vcc2025_h1)` Under `N_5`

| Target Line (`N`) | `N5 FCT` `[97.5% CI]` | `N5 Mis(vcc2025_h1)` `[97.5% CI]` | Unpaired `FCT_lo - Mis_hi` | Paired `FCT - Mis(vcc2025_h1)` `[97.5% CI]` | `N5 MaxMis` (All 6 Sources) |
|---|---:|---:|---:|---:|---|
| **`cd4t_rest` (`302`)** | `+0.0553` `[+0.0351, +0.0745]` | `+0.0207` `[+0.0046, +0.0360]` | **`-0.0009` (`FAIL`)** | **`+0.0345` `[+0.0261, +0.0432]` (`PASS`)** | `hepg2 -0.0047 (hi +0.0056)`, `jurkat +0.0099 (hi +0.0194)`, `kolf21j -0.0060 (hi +0.0011)`, `k562 -0.0050 (hi +0.0026)`, `rpe1 +0.0052 (hi +0.0133)`, `vcc +0.0207 (hi +0.0360)` |
| **`cd4t_stim48` (`503`)** | `+0.0725` `[+0.0559, +0.0893]` | `+0.0584` `[+0.0416, +0.0747]` | **`-0.0188` (`FAIL`)** | **`+0.0141` `[+0.0064, +0.0218]` (`PASS`)** | `hepg2 +0.0059 (hi +0.0132)`, `jurkat +0.0155 (hi +0.0233)`, `kolf21j +0.0073 (hi +0.0135)`, `k562 +0.0140 (hi +0.0237)`, `rpe1 -0.0016 (hi +0.0056)`, `vcc +0.0584 (hi +0.0747)` |

*Crucial takeaway:* On locked `P5` cohorts, **not only `cd4t_stim48`, but `cd4t_rest` as well fails the strict unpaired `FCT_lo > Mis(vcc2025_h1)_hi` margin under `N_5`** (`-0.0009`, overturning `RESULTS_13LINE_NULL_CONE_AUDIT.md` Section 3's non-`P5` `+0.0394 > +0.0378` `PASS`). Both pass when compared via **paired joint bootstrap** (`+0.0345` and `+0.0141`, because `FCT` and `Mis(vcc2025_h1)` co-vary positively across perturbations), though `cd4t_stim48`'s paired margin over `vcc2025_h1` under `N_5` is tiny (`+0.0141`), well below `Hybrid_ProjShared` (`+0.2199`).

---

## Check 3: `Q5a` Blocked Permutation & `LOO-line` Sign Test Across All `3` Gene Universes (`Bonferroni x3`)

Graph audit of the **`33` undirected pairs** (`>= 25` shared reproducible perturbations): edges span only **`10` active lines** (`cd4t_stim48` appears in `9` edges; `nadig_jurkat`, `replogle_k562_essential`, `replogle_rpe1`, `nadig_hepg2` in `7–8` edges each; `nourreddine_kolf21j` in `1` edge), and the **exact maximum number of vertex-disjoint independent pairs is `5`**. Treating `33` edges as i.i.d. in a Wilcoxon test is invalid.

| Gene Universe | Comparison | Mean Diff (`33` pairs) | Naive Edge Wins & Wilcoxon `p2` | Active Line Positives (`N=10`) & Binomial `p2` (`Bonf3`) | Exact `2^10` Line Sign-Flip `p2` (`Bonf3`) | Directed Block `2^10` `p2` (`Tgt / Src Bonf3`) | Dyadic CGM Cluster `t_9` & `p2` (`Bonf3`) | `LOO-Line` Min Remaining Mean (`All > 0?`) |
|---|---|---:|---|---|---|---|---|---|
| **`all`** | `Cal logFC − log1p` | `+0.0140` | `20/33` (`2.79e-02`) | `7/10` (`0.3438` -> **`1.0000`**) | `0.9941` -> **`1.0000`** | — | `t=1.65, p=0.1344` (**`0.4031`**) | `cd4t_rest (+0.0049)` (`True`) |
| **`all`** | `Cal FCT − log1p` | `+0.0093` | `30/33` (`3.41e-06`) | `9/10` (`0.0215` -> **`0.0645`**) | `0.6875` -> **`1.0000`** | `Tgt 1.0000 / Src 1.0000` | `t=2.94, p=0.0166` (**`0.0498`**) | `vcc2025_h1 (+0.0073)` (`True`) |
| **`all`** | `Raw FCT − log1p` | `+0.0068` | `30/33` (`3.45e-07`) | `9/10` (`0.0215` -> **`0.0645`**) | `0.0176` -> **`0.0527`** | `Tgt 0.0703 / Src 0.1113` | `t=2.86, p=0.0188` (**`0.0565`**) | `vcc2025_h1 (+0.0051)` (`True`) |
| **`expr01`** | `Cal logFC − log1p` | `+0.0269` | `26/33` (`4.17e-05`) | `9/10` (`0.0215` -> **`0.0645`**) | `0.3027` -> **`0.9082`** | — | `t=3.86, p=0.0038` (**`0.0115`**) | `cd4t_rest (+0.0190)` (`True`) |
| **`expr01`** | `Cal FCT − log1p` | `+0.0149` | `31/33` (`4.82e-08`) | `9/10` (`0.0215` -> **`0.0645`**) | `0.0078` -> **`0.0234`** | `Tgt 0.0234 / Src 0.8906` | `t=6.70, p=0.0001` (**`0.0003`**) | `vcc2025_h1 (+0.0129)` (`True`) |
| **`expr01`** | `Raw FCT − log1p` | `+0.0090` | `31/33` (`1.63e-09`) | `9/10` (`0.0215` -> **`0.0645`**) | `0.0039` -> **`0.0117`** | `Tgt 0.0234 / Src 0.0293` | `t=3.66, p=0.0052` (**`0.0157`**) | `vcc2025_h1 (+0.0069)` (`True`) |
| **`expr05`** | `Cal logFC − log1p` | `+0.0202` | `26/33` (`9.08e-06`) | `9/10` (`0.0215` -> **`0.0645`**) | `0.0137` -> **`0.0410`** | — | `t=3.08, p=0.0132` (**`0.0397`**) | `cd4t_rest (+0.0146)` (`True`) |
| **`expr05`** | **`Cal FCT − log1p`** | `+0.0115` | `32/33` (`4.66e-10`) | **`10/10`** (`0.0020` -> **`0.0059`**) | **`0.0020` -> `0.0059`** | **`Tgt 0.0059 / Src 0.0410`** | `t=6.76, p=0.0001` (**`0.0002`**) | `vcc2025_h1 (+0.0099)` (`True`) |
| **`expr05`** | `Raw FCT − log1p` | `+0.0072` | `32/33` (`4.66e-10`) | **`10/10`** (`0.0020` -> **`0.0059`**) | **`0.0020` -> `0.0059`** | `Tgt 0.0059 / Src 0.0117` | `t=3.78, p=0.0044` (**`0.0131`**) | `vcc2025_h1 (+0.0055)` (`True`) |

*(Raw uncalibrated means across `all`/`expr01`/`expr05`: `Raw logFC` loses to `Raw log1p` on `33/33`, `32/33`, `27/33` edges -> `0/10`, `0/10`, `1/10` lines, confirming shot-noise reversal in raw correlations. Pairwise outlier driving non-parametric permutation insensitivity in `all`/`expr01`: `nourreddine_kolf21j <-> cd4t_stim48` [`N=36`], the sole edge for `kolf21j`, gave `Cal logFC - log1p = -0.1843` and `Cal FCT - log1p = -0.0563` on `all`.)*

---

## Check 4: Noiseless Real-Distribution Simulation for Theorem 1 (`(mu_t, mu_s, logFC)`)

Simulated $\Delta_{c,p,g}^{\mathrm{sim}} = \log(1 + \mu_{c,g}\exp(\delta_{p,g})) - \log(1 + \mu_{c,g})$ using real $\mu_t$, real $\mu_s$, and real source geometric-mean $\delta_{p,g} = \log\overline{\mathrm{FC}}_{p,g}$ (`zero sampling noise, zero biological rewiring`), comparing both unclipped `logFC` and `|logFC| <= 2` (`noiseless_sim.csv`, `noiseless_sim.json`):

| Line (`Grp`, `N`) | Mean $w_t$ / $w_S$ (`r(w)`) | `real E1(all)` | `sim E1(all)` (`clip2`) | `Thm1 a_all` | `all` Sign Match? | `real E1(split)` | `sim E1(split)` (`clip2`) | `Thm1 a_split` | `split` Sign Match? |
|---|---|---:|---:|---:|:---:|---:|---:|---:|:---:|
| `nadig_hepg2` (R9, 797) | `0.321 / 0.361 (+0.829)` | `-0.0470` | `+0.1587` (`+0.1578`) | `+0.0072` | **False** | `+0.1977` | `+0.7111` (`+0.7172`) | `+0.0478` | **True** |
| `nadig_jurkat` (R9, 931) | `0.362 / 0.353 (+0.840)` | `-0.0782` | `+0.3486` (`+0.3487`) | `+0.0239` | **False** | `+0.0751` | `+0.6044` (`+0.6119`) | `+0.0325` | **True** |
| `nourreddine_kolf21j` (R9, 18) | `0.351 / 0.356 (+0.751)` | `-0.0172` | `+0.3103` (`+0.3103`) | `+0.0265` | **False** | `+0.1039` | `+0.7663` (`+0.7663`) | `+0.0763` | **True** |
| `replogle_k562_essential` (R9, 858) | `0.338 / 0.358 (+0.859)` | `-0.0377` | `+0.2149` (`+0.2136`) | `+0.0129` | **False** | `+0.1758` | `+0.6641` (`+0.6639`) | `+0.0309` | **True** |
| `replogle_rpe1` (R9, 1068) | `0.326 / 0.361 (+0.837)` | `-0.0934` | `+0.1738` (`+0.1737`) | `+0.0098` | **False** | `+0.1217` | `+0.7010` (`+0.7082`) | `+0.0646` | **True** |
| `vcc2025_h1` (R9, 80) | `0.431 / 0.339 (+0.488)` | `+0.0437` | `+0.5927` (`+0.5916`) | `+0.0902` | **True** | `+0.3100` | `+0.8635` (`+0.8643`) | `+0.2030` | **True** |
| `dld1` (R10, 86) | `0.346 / 0.355 (+0.861)` | `-0.0614` | `+0.2017` (`+0.2016`) | `+0.0122` | **False** | `+0.1294` | `+0.6161` (`+0.6181`) | `+0.0384` | **True** |
| `cd4t_rest` (R10, 302) | `0.421 / 0.355 (+0.491)` | `+0.0130` | `+0.6290` (`+0.6289`) | `+0.0885` | **True** | `+0.0532` | `+0.8594` (`+0.8608`) | `+0.2399` | **True** |
| `calu3` (R10, 8) | `0.309 / 0.355 (+0.735)` | `+0.0338` | `+0.2007` (`+0.2007`) | `+0.0122` | **True** | `+0.1385` | `+0.7147` (`+0.7147`) | `+0.0765` | **True** |
| **`cd4t_stim48` (R11, 503)** | `0.430 / 0.355 (+0.567)` | **`-0.0203`** | **`+0.6084`** (`+0.6084`) | `+0.0839` | **False** | **`+0.0541`** | **`+0.8691`** (`+0.8731`) | `+0.2182` | **True** |
| `telohaec` (R11, 9) | `0.320 / 0.381 (+0.788)` | `-0.0475` | `+0.0887` (`+0.0904`) | `+0.0013` | **False** | `+0.2188` | `+0.5642` (`+0.5668`) | `+0.0328` | **True** |
| `feng_ipsc` (R11, 93) | `0.377 / 0.381 (+0.904)` | `+0.0530` | `+0.2837` (`+0.2836`) | `+0.0164` | **True** | `+0.2468` | `+0.6964` (`+0.6962`) | `+0.0411` | **True** |

### Why Noiseless Theorem 1 Matches `SA-8 split` (`12/12`) Yet Fails `all_nontarget` (`4/12`)
1. **`resp_top50_split` (`12/12` positive):** True strong signal-tail conditioning overwhelmingly dominates noise/rewiring ($a_{\mathrm{top50},\delta^2} \in [+0.0309, +0.2399]$), validating Theorem 1's qualitative prediction that out-of-sample split-half selection still produces spurious positive `E1` under `N_1`.
2. **`all_nontarget` (`4/12` match; `+0.6084` noiseless vs `-0.0203` real on `cd4t_stim48`):** Under zero noise and zero rewiring ($\delta_t = \delta_s$), regression-to-the-mean between imperfectly correlated $w_t, w_S$ leaves $a_{\mathrm{all},\delta^2} > 0$ on ALL `12/12` lines — even stronger (`+0.6084`) on `cd4t_stim48` than the `+0.227` Gaussian toy table. In real data, **`E1 = g_{c_med}(FC_bar)` amplifies low/mid-expression source variance (`mid_mus < 0.5` contributes `65.0%` of `E1` residual variance on `cd4t_stim48` and `62.6%` on `replogle_rpe1`)**, while cross-line rewiring flattens empirical `Delta_obs ~ TMean` slopes across expression strata (`cd4t_stim48`: `b_obs = 0.1005` mid vs `0.0992` high, whereas `b_E1` jumps `1.8133` mid vs `0.7621` high; `replogle_rpe1`: `b_obs` rises `0.582 -> 0.852` while `b_E1` drops `1.139 -> 0.493`), forcing negative pooled covariance (`-96.2` on `cd4t_stim48` `mid_mus`; `-783.1` & `-622.7` on `replogle_rpe1`) — explaining why `8/12` real lines have `E1(all) < 0` despite noiseless Theorem 1 predicting `+`.

---

## Verification & Freeze Integrity
- `python3 -I round11/audit/exploratory_check/test_exploratory_check.py -v`: **4/4 PASSED** (`0.089s`).
- Output artifacts in `/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round11/audit/exploratory_check/outputs/`: `sa8_12lines.{json,csv}`, `n5_audit.json`, `n5_methods_12lines.csv`, `n5_perm100_12lines.csv`, `n5_spotlight_vs_vcc2025_h1.csv`, `q5a_blocked.json`, `q5a_33pairs_all_universes.csv`, `noiseless_sim.{json,csv}`, `paired_rorth_12lines.{csv,json}`.

---

## Check 5 (Post-Audit Addendum): Paired `FCT − Mis_s` Under Primary `r_orth_cal` (`all_nontarget`, 12 Lines)

**Added:** 2026-11-22T22:51+08:00 (per Master Controller `[主控 22:50]` verification request separating Check 2B's `N_5` `r_cone_sp` metric from `H-E1b`'s `r_orth_cal` metric)  
**Script & Outputs:** [`run_paired_rorth_12lines.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round11/audit/exploratory_check/run_paired_rorth_12lines.py) -> [`paired_rorth_12lines.csv`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round11/audit/exploratory_check/outputs/paired_rorth_12lines.csv), [`paired_rorth_12lines.json`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round11/audit/exploratory_check/outputs/paired_rorth_12lines.json)  
**Method & Parameters:** Frozen `round11/audit/independent_r11_metrics.py` (`evaluate_paired_cal_diff`, `evaluate_joint_residuals`, `evaluate_fct_perm_null`) under `all_nontarget`, locked `P5` eligibility (`only_reproducible=True`), R9 leave-one-out (`5` non-self R9 sources for R9 target lines; all `6` R9 sources for R10/R11 lines), **`seed = 20261124`**, `n_boot = 4,000`, `n_perm = 100`, two-sided Bonferroni **`97.5%` CIs (`[1.25%, 98.75%]`)**.

### Summary Counts Under `r_orth_cal` vs Check 2B (`N_5 r_cone_sp`)
- **All Mismatched Baselines `97.5% CI_lo > 0` (`H-E1b`-type paired test):** **`11/12` lines** (`6/6` unseen R10+R11 lines pass all `6/6` sources; `5/6` R9 LOO lines pass all `5/5` sources). Only `nourreddine_kolf21j` (`N = 18`, R9 LOO) fails, missing `2/5` pairs under `r_orth_cal` (`vs nadig_jurkat +0.0568 [-0.0008, +0.1162]` and `vs replogle_k562_essential +0.0436 [-0.0105, +0.0917]`; vs `1/5` failing pairs under `N_5 r_cone_sp`).
- **Legacy Unpaired `H-D2` (`FCT CI97.5_lo > max(Mis_hi, Perm100 Q97.5)`):** **`7/12` lines** (`nadig_hepg2`, `nadig_jurkat`, `replogle_k562_essential`, `replogle_rpe1`, `vcc2025_h1`, `dld1`, `feng_ipsc` pass; `nourreddine_kolf21j`, `cd4t_rest`, `calu3`, `cd4t_stim48`, `telohaec` fail). Thus paired evaluation flips `cd4t_rest` (`302`), `cd4t_stim48` (`503`), `calu3` (`8`), and `telohaec` (`9`) from unpaired `FAIL` to `6/6 PASS`.
- **`H-E1a` + `H-E1b` + `H-E1c` Jointly Under `r_orth_cal`:** **`11/12` lines** (`12/12` clear `H-E1a` with `valid_frac = 1.0`, `rho_ab_orth > 0` and `H-E1c` `FCT CI97.5_lo > Perm100 Q97.5`).

### Full Line-by-Line Results (All 6 Source Contrasts per Line, `all_nontarget`, `r_orth_cal`)

| Line (`Grp`, `N`) | `FCT r_orth_cal` `[97.5% CI]` | `Perm Q97.5` | Unpaired `H-D2` | Paired Pass | `vs nadig_hepg2` | `vs nadig_jurkat` | `vs kolf21j` | `vs k562` | `vs rpe1` | `vs vcc2025_h1` |
|---|---|---:|:---:|:---:|---|---|---|---|---|---|
| `nadig_hepg2` (R9, `797`) | `+0.1848` `[+0.1772, +0.1929]` | `-0.0224` | `PASS` | **`5/5`** | *LOO self* | `+0.1217` `[+0.1147, +0.1293]` | `+0.1578` `[+0.1503, +0.1657]` | **`+0.1179`** `[+0.1108, +0.1255]` | `+0.1212` `[+0.1126, +0.1301]` | `+0.1825` `[+0.1763, +0.1893]` |
| `nadig_jurkat` (R9, `931`) | `+0.1081` `[+0.0975, +0.1185]` | `-0.0475` | `PASS` | **`5/5`** | `+0.0839` `[+0.0726, +0.0954]` | *LOO self* | `+0.0962` `[+0.0875, +0.1049]` | **`+0.0285`** `[+0.0231, +0.0342]` | `+0.1251` `[+0.1158, +0.1348]` | `+0.1420` `[+0.1336, +0.1509]` |
| `nourreddine_kolf21j` (R9, `18`) | `+0.1142` `[+0.0399, +0.1809]` | `-0.0025` | `FAIL` | **`3/5`** | `+0.0665` `[+0.0243, +0.1061]` | `+0.0568` **`[-0.0008, +0.1162]`** | *LOO self* | `+0.0436` **`[-0.0105, +0.0917]`** | `+0.0893` `[+0.0326, +0.1493]` | `+0.0922` `[+0.0348, +0.1546]` |
| `replogle_k562_essential` (R9, `858`) | `+0.1614` `[+0.1506, +0.1720]` | `-0.0186` | `PASS` | **`5/5`** | **`+0.1191`** `[+0.1103, +0.1277]` | `+0.1239` `[+0.1183, +0.1298]` | `+0.1456` `[+0.1381, +0.1531]` | *LOO self* | `+0.1878` `[+0.1768, +0.1987]` | `+0.1717` `[+0.1637, +0.1797]` |
| `replogle_rpe1` (R9, `1068`) | `+0.1510` `[+0.1429, +0.1590]` | `-0.0529` | `PASS` | **`5/5`** | `+0.1061` `[+0.0988, +0.1134]` | **`+0.0668`** `[+0.0614, +0.0723]` | `+0.1176` `[+0.1121, +0.1231]` | `+0.0852` `[+0.0787, +0.0918]` | *LOO self* | `+0.1812` `[+0.1743, +0.1884]` |
| `vcc2025_h1` (R9, `80`) | `+0.1427` `[+0.1168, +0.1707]` | `+0.0386` | `PASS` | **`5/5`** | `+0.0976` `[+0.0724, +0.1252]` | `+0.1247` `[+0.1012, +0.1494]` | **`+0.0890`** `[+0.0717, +0.1077]` | `+0.1145` `[+0.0903, +0.1403]` | `+0.1560` `[+0.1221, +0.1908]` | *LOO self* |
| `dld1` (R10, `86`) | `+0.1355` `[+0.1083, +0.1639]` | `-0.0352` | `PASS` | **`6/6`** | `+0.0966` `[+0.0719, +0.1238]` | `+0.1292` `[+0.1073, +0.1523]` | `+0.1140` `[+0.0947, +0.1341]` | `+0.1159` `[+0.0922, +0.1415]` | **`+0.0802`** `[+0.0626, +0.0994]` | `+0.1734` `[+0.1489, +0.1984]` |
| `cd4t_rest` (R10, `302`) | `+0.0556` `[+0.0373, +0.0736]` | `+0.0131` | `FAIL` | **`6/6`** | `+0.0486` `[+0.0253, +0.0716]` | `+0.0392` `[+0.0251, +0.0521]` | `+0.0538` `[+0.0368, +0.0706]` | `+0.0514` `[+0.0344, +0.0679]` | `+0.0398` `[+0.0188, +0.0591]` | **`+0.0347`** `[+0.0262, +0.0431]` |
| `calu3` (R10, `8`) | `+0.0820` `[+0.0381, +0.1591]` | `+0.0320` | `FAIL` | **`6/6`** | `+0.0628` `[+0.0115, +0.1388]` | `+0.1164` `[+0.0545, +0.1949]` | `+0.0906` `[+0.0403, +0.1539]` | `+0.1030` `[+0.0463, +0.1700]` | `+0.0451` `[+0.0046, +0.0945]` | **`+0.0716`** `[+0.0045, +0.1468]` |
| `cd4t_stim48` (R11, `503`) | `+0.0532` `[+0.0346, +0.0725]` | `-0.0129` | `FAIL` | **`6/6`** | `+0.0456` `[+0.0252, +0.0664]` | `+0.0374` `[+0.0256, +0.0496]` | `+0.0516` `[+0.0369, +0.0676]` | `+0.0356` `[+0.0207, +0.0508]` | `+0.0517` `[+0.0325, +0.0721]` | **`+0.0165`** `[+0.0086, +0.0247]` |
| `telohaec` (R11, `9`) | `+0.2173` `[+0.1047, +0.3559]` | `-0.0132` | `FAIL` | **`6/6`** | **`+0.1335`** `[+0.0378, +0.2270]` | `+0.2608` `[+0.1188, +0.4411]` | `+0.2130` `[+0.1008, +0.3432]` | `+0.2276` `[+0.0842, +0.3803]` | `+0.1403` `[+0.0483, +0.2352]` | `+0.2851` `[+0.1615, +0.4703]` |
| `feng_ipsc` (R11, `93`) | `+0.1212` `[+0.1037, +0.1377]` | `+0.0452` | `PASS` | **`6/6`** | `+0.0775` `[+0.0636, +0.0912]` | `+0.1236` `[+0.1052, +0.1422]` | `+0.0650` `[+0.0525, +0.0778]` | `+0.1143` `[+0.0962, +0.1319]` | `+0.1446` `[+0.1245, +0.1630]` | **`+0.0543`** `[+0.0407, +0.0674]` |

*Note:* Bold deltas highlight each passing line's tightest contrast (`min_pair` by `CI975_lo`, or `min_pair` by point estimate on `cd4t_rest` where `vcc2025_h1` gives min diff `+0.0347` and `replogle_rpe1` gives min `CI975_lo` `+0.0188`), plus the two failing pairs (`CI975_lo <= 0`) on `nourreddine_kolf21j`. Independent replication by Master Controller (`controller/checks/r12_paired_check.py`, `seed = 777`, `1,000` bootstraps) on `dld1`, `cd4t_rest`, `cd4t_stim48`, and `feng_ipsc` matches identically (`6/6` pairs `> 0` on all four lines).
