# Round 10 Independent Group C Statistical Audit & `C3` Cross-Check Report (`A-T3`)

- **Pre-Registration Protocol:** [`round10/prereg/PREREG_C.md`](../../prereg/PREREG_C.md) (`commit dc4b4a5`, `sha256: 0fb48aeb08aa6befe9583030f33ec668155e9e937616d291ea0c211e0d1bf4c6`)
- **Code Freeze & QC Lock:** `CODE_FREEZE_R10` (`commit f8d2ebfc238d83a0c3a6c86696950ed292bf90e9`), `QC_LOCK` (`commit 86872ed741c3aa681a7a9091e60d03bd9922038e`), `FREEZE_MANIFEST_R10.txt` (`22/22` `.py` hashes verified `OK` via `shasum -a 256 -c`)
- **Execution Mode:** Strict isolated Python (`python3 -I` via `runpy.run_path` isolating `sys.path` to `round10/audit/`, zero `.py` modifications, `n_boot = 4,000`, `n_perm = 100`, `seed = 20261122`, Bonferroni two-sided `97.5%` CIs)
- **Evaluated Cohort (`CONFIRMATORY_LINES.txt`):** `K_conf = 2` (`dld1` `N = 86`, `cd4t_rest` `N = 302` [`284` half-valid across 2 pre-locked donors `CE0006864`/`CE0008162`]); Secondary lines: `calu3` (`N = 8`), `wtc11_i3neuron` (`N = 2`). (`h9_neuron` excluded pre-download per `P2`.)

---

## 1. Confirmatory & Secondary Adjudication (`A-T3` Independent Verdict)

- **Primary Scientific Claim (`H-C1` & `H-C2`, `K_conf = 2` Unanimous Rule):** **`未确认 (NOT CONFIRMED)`**
  - **`H-C1` (`FCT r_orth_cal 97.5% CI_lo > 0`):** **`2/2` PASS** (`dld1`, `cd4t_rest`)
  - **`H-C2` (`FCT CI97.5_lo > max(6 Mismatched Source Means CI975_hi, Perm100 p97.5)`):** **`0/2` FAIL** (`dld1`: `+0.0870 <= max(+0.1092, +0.1179)`; `cd4t_rest`: `+0.0171 <= max(+0.0761, +0.0554)`)
- **Parallel Practical Non-Inferiority (`H-C3`, `6` Conventions vs. `Transfer_Mean`):** **`2/2` PASS (`非劣效通过 / H-C3 PASS`)** (`0/6` inferior conventions on both `dld1` and `cd4t_rest`; `1` strictly superior on `cd4t_rest` `cent_de50_r` `+0.0305 [+0.0177, +0.0435]`).

### 1.1 Primary Endpoint, Controls, & `SA-1`–`SA-5b` Summary

| Line (Role) | `N` (`\|G_eval\|`) | `FCT r_orth_cal [97.5% CI]` | Max Mismatch `CI975_hi` | `Perm100` `med / p97.5` | `H-C1 / H-C2 / H-C3` | `FCT slope_ols [97.5% CI]` | `all_nontarget` `r_orth_cal` ($\rho^{\perp}_{AB}$ `top50 -> all`) | `top75` / `FullNT` |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`dld1`** (Conf.) | `86` (`6,124`) | `+0.1401 [+0.0870, +0.1917]` | `+0.1092` | `+0.0892 / +0.1179` | `T / F / T` | `+0.0290 [+0.0194, +0.0396]` | `+0.1355` (`0.8277 -> 0.2533`) | `+0.1433 / +0.1401` |
| **`cd4t_rest`** (Conf.) | `302` (`6,123`) | `+0.0560 [+0.0171, +0.0938]` | `+0.0761` | `+0.0347 / +0.0554` | `T / F / T` | `+0.0156 [+0.0026, +0.0313]` | `+0.0556` (`0.6479 -> 0.1468`) | `+0.0543 / +0.0560` |
| `calu3` (Sec.) | `8` (`6,134`) | `+0.0879 [-0.0689, +0.2743]` | `+0.2702` | `+0.0776 / +0.1519` | `F / F / T` | `+0.0068 [-0.0181, +0.0400]` | `+0.0820` (`0.8776 -> 0.3631`) | `+0.0897 / +0.0879` |
| `wtc11_i3neuron` (Sec.) | `2` (`5,995`) | `-0.0313 [-0.0334, +0.2649]` | `+0.3224` | `-0.0219 / +0.1681` | `F / F / T` | `-0.0329 [-0.0684, +0.0235]` | `-0.0122` (`0.8562 -> 0.2223`) | `-0.0313 / -0.0313` |

- **6 Source Mismatched `ctrl_mean_libnorm` Controls (`r_orth_cal [97.5% CI]`):**
  - `dld1`: `nadig_hepg2 +0.0294 [-0.0174, +0.0726]`, `nadig_jurkat -0.0225 [-0.0739, +0.0297]`, `nourreddine_kolf21j +0.0300 [-0.0138, +0.0759]`, `replogle_k562_essential -0.0058 [-0.0557, +0.0447]`, `replogle_rpe1 +0.0463 [-0.0196, +0.1032]`, `vcc2025_h1 +0.0707 [+0.0331, +0.1092]`.
  - `cd4t_rest`: `nadig_hepg2 +0.0238 [hi +0.0599]`, `nadig_jurkat +0.0398 [hi +0.0761]`, `nourreddine_kolf21j +0.0213 [hi +0.0548]`, `replogle_k562_essential +0.0348 [hi +0.0691]`, `replogle_rpe1 +0.0353 [hi +0.0714]`, `vcc2025_h1 +0.0294 [hi +0.0680]`.
- **`SA-2` & `SA-3` Benchmark Context (`r_orth_cal` / `r_cal` / `slope_ols`):**
  - `dld1`: `GlobalScale` `0.0000` / `+0.2404` / `+0.0369`; `Hybrid_ProjShared (3D 1ca9f60f)` `-0.0429` / `-0.0461` / `-0.0053`; `TMVS` `+0.0888` / `+0.2149` / `+0.0536`; `Hybrid_on_TMVS` `+0.0373` / `+0.1726` / `+0.0483`; `CIPHER_Linear` `-0.0314` / `-0.2441` / `-0.1590`.
  - `cd4t_rest`: `GlobalScale` `0.0000` / `-0.2616` / `-0.0452`; `Hybrid_ProjShared` `-0.0244` / `-0.0286` / `-0.0047`; `TMVS` `+0.0198` / `-0.2445` / `-0.0896`; `Hybrid_on_TMVS` `-0.0040` / `-0.2267` / `-0.0943`; `CIPHER_Linear` `-0.0324` / `+0.2182` / `+0.0954`.

---

## 2. `H-C3` Six Conventions & `SA-4` / `R19` Decomposition

- **`H-C3` Confirmatory Non-Inferiority (`FCT - Transfer_Mean`, `97.5% CI`):**
  - `dld1`: `ctrl_disc_panel -0.0011 [-0.0059, +0.0038]`, `ctrl_disc_ownmasked -0.0018 [-0.0073, +0.0033]`, `ctrl_de50_r +0.0065 [-0.0018, +0.0150]`, `cent_disc_panel -0.0030 [-0.0120, +0.0040]`, `cent_disc_ownmasked -0.0030 [-0.0116, +0.0038]`, `cent_de50_r +0.0053 [-0.0071, +0.0179]`, descriptive `pds_cosine_panel -0.0055 [-0.0114, -0.0003]`.
  - `cd4t_rest`: `ctrl_disc_panel -0.0001 [-0.0084, +0.0079]`, `ctrl_disc_ownmasked -0.0005 [-0.0088, +0.0073]`, `ctrl_de50_r -0.0011 [-0.0117, +0.0093]`, `cent_disc_panel +0.0130 [-0.0025, +0.0283]`, `cent_disc_ownmasked +0.0124 [-0.0028, +0.0272]`, `cent_de50_r +0.0305 [+0.0177, +0.0435]` (`sup=True`), descriptive `pds_cosine_panel +0.0001 [-0.0074, +0.0072]`.
- **`R19` `all_nontarget` Three-Tier Decomposition (`L1 / L2 / L3`):**
  - `dld1`: `beta_t = -0.4872`, `0.1133 / 0.0271 / 0.8595`.
  - `cd4t_rest`: `beta_t = -0.8980`, `0.3098 / 0.0031 / 0.6872`.

---

## 3. `S7` (`P6`) Cell-Count & Knockdown-Depth Confounder Controls (`SA-6`)

- **`dld1` (`n_complete = 86`):**
  - `Spearman(r_orth_p, ln n_cells) = -0.0005 [-0.2382, +0.2399]`; `Spearman(r_orth_p, |k_p|) = +0.1322 [-0.1074, +0.3603]`.
  - Centered OLS `[1, ln(n_cells), |k_p|]`: intercept **`+0.1400 [+0.0979, +0.1826]`** (`slope_ln_ncells +0.0064`, `slope_abs_kp +0.0217`).
  - `n_cells` Tertiles (`T1=29, T2=29, T3=28`): `T1` (`25–33` cells) `+0.1403 [+0.0654]` (`H-C1=T, H-C2=F`), `T2` `+0.1456 [+0.0434]`, `T3` `+0.1290 [+0.0171]`.
- **`cd4t_rest` (`n_complete = 302`):**
  - `Spearman(r_orth_p, ln n_cells) = -0.1428 [-0.2637, -0.0167]`; `Spearman(r_orth_p, |k_p|) = +0.0828 [-0.0431, +0.2131]`.
  - Centered OLS `[1, ln(n_cells), |k_p|]`: intercept **`+0.0607 [+0.0281, +0.0934]`** (`slope_ln_ncells -0.0391`, `slope_abs_kp +0.0772`).
  - `n_cells` Tertiles (`T1=101, T2=101, T3=100`): `T1` (`25–86` cells) `+0.0988 [+0.0443]` (`H-C1=T, H-C2=F` vs `maxmis +0.1203` / `perm +0.0745`), `T2` `-0.0030 [-0.0967]`, `T3` `+0.0369 [-0.0310]`.

---

## 4. Dual-Team `C3` Cross-Check vs. `A-T2` (`round10/results/line_metrics_holdout_r10_log1p_primary.csv`)

- **Core `run_cross_check` Suite (`349` checks across all 4 Round 10 lines, stored in [`cross_check_r10.json`](cross_check_r10.json)):**
  - **`max_point_abs_diff = 6.661e-16`** (`< 1e-6`, `1.5e9`x below tolerance)
  - **`max_ci_abs_diff = 1.069e-05`** (`<= 0.01`, `935`x below tolerance)
  - **`passes_c3_tolerance = True`**
- **Extended `H-C3` Point + `97.5%` CI Parity (`36` confirmatory checks):** **`max_abs_diff = 2.351e-16`**.
- **Hypothesis & Governance Consensus (`100%` identical):**
  - `K_conf = 2` (`dld1`, `cd4t_rest`), `H-C1 = 2/2 PASS` (`valid_frac = 1.0`, `rho_ab_orth > 0` -> zero `D-R10-2` divergence), `H-C2 = 0/2 FAIL` (**`未确认`**), `H-C3 = 2/2 PASS` (`0` inferior pairs / `1` strictly superior).
  - Official CLI cross-check (`run_r10_audit.py --cross-check-csv .../line_metrics_holdout_r10_log1p_primary.csv`) is stored verbatim in [`audit_r10_summary.json`](audit_r10_summary.json) (`"cross_check"` block: `n_checks = 349`, `max_point_abs_diff = 6.661e-16`, `max_ci_abs_diff = 1.069e-05`, `passes_c3_tolerance = True`).
  - Note on `python3 -I` isolated invocation: frozen `round10/audit/run_r10_audit.py` and `round10/src/run_eval_r10.py` omit `sys.path.insert` inside the script, so both `A-T2` and `A-T3` executed strictly isolated `python3 -I` via `runpy.run_path` with `sys.path` set exclusively to the team's own frozen directory (leaving `22/22` manifest `.py` hashes untouched).

---

## 5. Sensitivity Runs (`S1 --all-targets` & `S5 --s5-top75`) Dual-Team Verification

- **`S1` (`--all-targets`, [`audit_r10_s1_alltargets_summary.json`](audit_r10_s1_alltargets_summary.json) vs. `A-T2` `line_metrics_holdout_r10_log1p_s1.csv`):**
  - **`364` `run_cross_check` checks across all 4 lines:** **`max_point_abs_diff = 1.978e-16`**, **`max_ci_abs_diff = 3.886e-16`**, **`passes_c3_tolerance = True`**.
  - Confirmatory lines: `dld1` (`N = 500`, `FCT r_orth_cal = +0.0717 [+0.0462, +0.0992]`, `max_mis_hi = +0.0406`, `perm100_p975 = +0.0570` -> `H-C1=T, H-C2=F`); `cd4t_rest` (`N = 1,422`, `+0.0314 [+0.0124, +0.0502]`, `max_mis_hi = +0.0232`, `perm100_p975 = +0.0183` -> `H-C1=T, H-C2=F`).
  - Identical `S1` verdicts: `H-C1` `2/2` PASS, `H-C2` `0/2` FAIL (**`未确认`**), `H-C3` `0/2` FAIL (**`未通过非劣效`**, `4` inferior pairs in `cent_*` once non-reproducible targets are included).
- **`S5` (`--s5-top75`, `audit_r10_summary.json` `sa5_top75_fct` vs. `A-T2` `line_metrics_holdout_r10_log1p_s5.csv`):**
  - `dld1`: `FCT r_orth_cal = +0.143259734837` (point diff `8.327e-17`); `cd4t_rest`: `+0.054313591518` (point diff `5.551e-17`). Both confirmatory runs preserve primary verdicts (`H-C1` `2/2`, `H-C2` `0/2` **`未确认`**, `H-C3` `2/2` **`非劣效通过`**).
