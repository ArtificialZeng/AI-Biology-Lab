# Round 11 (`PREREG_D` + Addendum A) Independent Statistical Audit (`dev_r9r10`)

- **Parameters:** `n_boot = 4000`, `n_perm = 100`, `seed = 20261122`, Bonferroni two-sided `97.5%` CIs (`[1.25%, 98.75%]`), `--source-mode r9_6`.
- **Adjudication Status:** `NO_CONFIRMATORY_ADJUDICATION` (`K_conf = 0`)
- **P1 Notice:** `Insufficient confirmatory cell lines (K_conf < 2); confirmatory hypothesis testing was not performed (确认性线不足，未做确认性检验).`

## Table 1: Primary Whole-Transcriptome Endpoints (`all_nontarget`, `97.5% CI`)

| Line (`N`) | `FCT` `r_orth_cal [97.5% CI]` | `valid_frac` / `rho_AB` | `E1_ConstMed` | `FCT - E1 [97.5% CI]` | `max Mis97.5_hi` | `Perm100 Q975` | `H-D1 / D2 / D3` |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `dld1` (`N=86`) | `+0.1355 [+0.1082, +0.1649]` | `1.000` / `+0.2533` | `-0.0614 [-0.0773, -0.0457]` | `+0.1969 [+0.1673, +0.2253]` | `+0.0793` | `-0.0352` | `True / True / True` |

## Table 2: Selection Trap (`H-D4` & `SA-8 resp_top50_split`) and Dual-Mask Permutations

| Line | `E1(resp_top50) [97.5% CI]` | `E1(top50 - all) [97.5% CI]` | `SA-8 FCT split` | `SA-8 E1 split` | `SA-8 Perm100 Med / Q975` | `Perm100 top50 Q975` | `H-D4a / D4b` |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `dld1` | `+0.1567 [+0.1068, +0.2060]` | `+0.2180 [+0.1637, +0.2719]` | `+0.1407 [+0.0756, +0.2058]` | `+0.1294 [+0.0900, +0.1713]` | `+0.0689 / +0.1008` | `+0.1179` | `True / True` |

## Table 3: `H-D5` (`all_nontarget` FCT vs `Mis(nourreddine_kolf21j)`) & Addendum Notes

| Line | `FCT` (`all_nontarget`) | `Mis(kolf21j)` (`all_nontarget`) | `FCT - Mis(kolf21j) [97.5% CI]` | `Diff CI > 0` |
| :--- | :---: | :---: | :---: | :---: |
| `dld1` | `+0.1355 [+0.1082, +0.1649]` | `+0.0215 [+0.0078, +0.0362]` | `+0.1140 [+0.0942, +0.1343]` | `True` |
