# Round 11 (`PREREG_D` + Addendum A) Independent Statistical Audit (`dev_r9r10`)

- **Parameters:** `n_boot = 4000`, `n_perm = 100`, `seed = 20261123`, Bonferroni two-sided `97.5%` CIs (`[1.25%, 98.75%]`), `--source-mode r9_6`.
- **Adjudication Status:** `NO_CONFIRMATORY_ADJUDICATION` (`K_conf = 0`)
- **P1 Notice:** `Insufficient confirmatory cell lines (K_conf < 2); confirmatory hypothesis testing was not performed (确认性线不足，未做确认性检验).`

## Table 1: Primary Whole-Transcriptome Endpoints (`all_nontarget`, `97.5% CI`)

| Line (`N`) | `FCT` `r_orth_cal [97.5% CI]` | `valid_frac` / `rho_AB` | `E1_ConstMed` | `FCT - E1 [97.5% CI]` | `max Mis97.5_hi` | `Perm100 Q975` | `H-D1 / D2 / D3` |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `dld1` (`N=86`) | `+0.1355 [+0.1079, +0.1620]` | `1.000` / `+0.2533` | `-0.0614 [-0.0774, -0.0467]` | `+0.1969 [+0.1672, +0.2251]` | `+0.0781` | `-0.0352` | `True / True / True` |
| `vcc2025_h1` (`N=80`) | `+0.1427 [+0.1172, +0.1712]` | `1.000` / `+0.8886` | `+0.0437 [+0.0282, +0.0624]` | `+0.0989 [+0.0832, +0.1168]` | `+0.0669` | `+0.0386` | `True / True / True` |

## Table 2: Selection Trap (`H-D4` & `SA-8 resp_top50_split`) and Dual-Mask Permutations

| Line | `E1(resp_top50) [97.5% CI]` | `E1(top50 - all) [97.5% CI]` | `SA-8 FCT split` | `SA-8 E1 split` | `SA-8 Perm100 Med / Q975` | `Perm100 top50 Q975` | `H-D4a / D4b` |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `dld1` | `+0.1567 [+0.1076, +0.2033]` | `+0.2180 [+0.1639, +0.2697]` | `+0.1407 [+0.0722, +0.2014]` | `+0.1294 [+0.0887, +0.1692]` | `+0.0689 / +0.1008` | `+0.1179` | `True / True` |
| `vcc2025_h1` | `+0.3076 [+0.2427, +0.3742]` | `+0.2639 [+0.2052, +0.3210]` | `+0.3533 [+0.2815, +0.4215]` | `+0.3100 [+0.2464, +0.3752]` | `+0.2130 / +0.2472` | `+0.2445` | `True / True` |
