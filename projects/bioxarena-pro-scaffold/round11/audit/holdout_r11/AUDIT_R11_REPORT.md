# Round 11 (`PREREG_D` + Addendum A) Independent Statistical Audit (`holdout_r11`)

- **Parameters:** `n_boot = 4000`, `n_perm = 100`, `seed = 20261123`, Bonferroni two-sided `97.5%` CIs (`[1.25%, 98.75%]`), `--source-mode r9_6`.
- **Adjudication Status:** `NO_CONFIRMATORY_ADJUDICATION` (`K_conf = 1`)
- **P1 Notice:** `Insufficient confirmatory cell lines (K_conf < 2); confirmatory hypothesis testing was not performed (确认性线不足，未做确认性检验).`

## Table 1: Primary Whole-Transcriptome Endpoints (`all_nontarget`, `97.5% CI`)

| Line (`N`) | `FCT` `r_orth_cal [97.5% CI]` | `valid_frac` / `rho_AB` | `E1_ConstMed` | `FCT - E1 [97.5% CI]` | `max Mis97.5_hi` | `Perm100 Q975` | `H-D1 / D2 / D3` |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `cd4t_stim48` (`N=503`) | `+0.0532 [+0.0337, +0.0731]` | `1.000` / `+0.0784` | `-0.0203 [-0.0326, -0.0087]` | `+0.0735 [+0.0561, +0.0917]` | `+0.0572` | `-0.0129` | `True / False / True` |
| `telohaec` (`N=9`) | `+0.2173 [+0.0958, +0.3602]` | `1.000` / `+0.1248` | `-0.0475 [-0.1700, +0.0187]` | `+0.2648 [+0.1443, +0.4636]` | `+0.1712` | `-0.0132` | `True / False / True` |
| `feng_ipsc` (`N=93`) | `+0.1212 [+0.1043, +0.1371]` | `1.000` / `+0.5570` | `+0.0530 [+0.0287, +0.0785]` | `+0.0682 [+0.0427, +0.0923]` | `+0.0829` | `+0.0452` | `True / True / True` |

## Table 2: Selection Trap (`H-D4` & `SA-8 resp_top50_split`) and Dual-Mask Permutations

| Line | `E1(resp_top50) [97.5% CI]` | `E1(top50 - all) [97.5% CI]` | `SA-8 FCT split` | `SA-8 E1 split` | `SA-8 Perm100 Med / Q975` | `Perm100 top50 Q975` | `H-D4a / D4b` |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `cd4t_stim48` | `+0.0522 [+0.0257, +0.0792]` | `+0.0725 [+0.0477, +0.0982]` | `+0.0621 [+0.0231, +0.1031]` | `+0.0541 [+0.0193, +0.0896]` | `+0.0381 / +0.0571` | `+0.0481` | `True / True` |
| `telohaec` | `+0.1798 [+0.0334, +0.3756]` | `+0.2273 [+0.0659, +0.4838]` | `+0.2579 [-0.0858, +0.6027]` | `+0.2188 [-0.0366, +0.4558]` | `+0.1070 / +0.2343` | `+0.1723` | `True / True` |
| `feng_ipsc` | `+0.2479 [+0.2003, +0.2975]` | `+0.1949 [+0.1457, +0.2429]` | `+0.2440 [+0.1909, +0.2936]` | `+0.2468 [+0.1982, +0.2981]` | `+0.1482 / +0.1874` | `+0.1876` | `True / True` |

## Table 3: `H-D5` (`all_nontarget` FCT vs `Mis(nourreddine_kolf21j)`) & Addendum Notes

| Line | `FCT` (`all_nontarget`) | `Mis(kolf21j)` (`all_nontarget`) | `FCT - Mis(kolf21j) [97.5% CI]` | `Diff CI > 0` |
| :--- | :---: | :---: | :---: | :---: |
| `cd4t_stim48` | `+0.0532 [+0.0337, +0.0731]` | `+0.0016 [-0.0089, +0.0121]` | `+0.0516 [+0.0360, +0.0672]` | `True` |
- **Disclosure (`cd4t_stim48`):** Addendum A2: halves CE0008678 and CE0010866 share the 24 sequencing lanes GSM9393875..98 (channel-NT subtracted, residual shared channel noise elevates ceiling / lowers r_orth_cal conservatively); PREREG_D Q2 weaker independence: shares CD4+ alpha-beta T lineage, library, & Marson lab with R10 cd4t_rest.
| `telohaec` | `+0.2173 [+0.0958, +0.3602]` | `+0.0042 [-0.0253, +0.0441]` | `+0.2130 [+0.0953, +0.3476]` | `True` |
| `feng_ipsc` | `+0.1212 [+0.1043, +0.1371]` | `+0.0562 [+0.0378, +0.0744]` | `+0.0650 [+0.0528, +0.0784]` | `True` |
- **Disclosure (`feng_ipsc`):** Addendum A3: sorted 17-line parity halves split clonal replicates across the same HipSci donors; split-half ceiling captures clonal + channel rather than inter-donor variance. Secondary-only (H-D5).
