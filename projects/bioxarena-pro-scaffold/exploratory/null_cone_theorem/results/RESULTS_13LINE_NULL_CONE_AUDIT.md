# Universal Target-Free Null-Cone Audit Across 13 Human Perturb-seq Atlases (>10.33M QC Single Cells)

- **Audit Parameters**: `n_boot = 4000`, `n_perm = 100`, `seed = 20261124`, Bonferroni two-sided `97.5%` CIs (`[1.25%, 98.75%]`), universal target-free `log1p` scales `c ∈ [0.05, 1.0, 20.0]`.
- **Scope**: `12` eval-valid cell lines (`N >= 9` reproducible shared perturbations; `wtc11_i3neuron` `N = 2` skipped) across Round 9 (6 lines LOLO), Round 10 (3 lines, R9 6-source pool), and Round 11 (3 lines, R9 6-source pool).
- **Demarcation**: All Round 9 (`PREREG_A v4`), Round 10 (`PREREG_C`), and Round 11 (`PREREG_D`) locked confirmatory verdicts remain strictly unchanged; the `N_4`/`N_5` null-cone cure, cross-benchmark audit, and `Q5a` conservation analysis are explicitly demarcated as **post-hoc / mechanistic methodological discoveries** built on the unblinded 13-line corpus.

## 1. Pure Noiseless Cross-Benchmark Audit (Theorem 1 & Theorem 2)

| Metric / Benchmark | `TMean` | `GlobalScale` | `SrcGrandMean` (`b_S`) | `E1` (zero-tgt) | `E6` (`logFC`) | `Perm` (scrambled) | `FCT` (true `μ_t`) |
|---|---:|---:|---:|---:|---:|---:|---:|
| PearsonDelta (all) | +0.9492 | +0.9492 | +0.2155 | +0.8586 | +0.8483 | +0.7368 | +1.0000 |
| Systema/SBB PearsonD50 (split) | +0.9693 | +0.9693 | +0.2665 | +0.9437 | +0.9676 | +0.8207 | +1.0000 |
| Miller/cell-eval2 top100 DEG r | +0.9684 | +0.9684 | +0.2855 | +0.9438 | +0.9648 | +0.8214 | +1.0000 |
| Miller/SBB DEG-weighted R2_w | +0.9359 | +0.8710 | -0.0776 | +0.8413 | +0.7987 | +0.6977 | +1.0000 |
| CGC/Molina Resid r (all) | +0.0000 | +0.1698 | -0.1660 | +0.1480 | +0.2876 | +0.1167 | +1.0000 |
| CGC/Molina Resid r (top50 split) | +0.0000 | +0.4183 | -0.4163 | -0.0233 | +0.7176 | -0.0060 | +1.0000 |
| N1 1D r_orth (all) | +0.0000 | +0.0000 | -0.0030 | +0.2270 | +0.2359 | +0.1620 | +1.0000 |
| N1 1D r_orth (top50 split) | +0.0000 | +0.0000 | +0.0142 | +0.4641 | +0.6556 | +0.2710 | +1.0000 |
| OURS N5 Cone Semi-Partial (all) | +0.0000 | +0.0000 | +0.0000 | -0.0000 | -0.0001 | +0.0090 | +0.9500 |
| OURS N5 Cone Semi-Partial (top50) | +0.0000 | +0.0000 | -0.0000 | -0.0016 | -0.0059 | -0.0095 | +0.5700 |

- **Universal Cone Low-Rank Proof**: Over 6 orders of magnitude `c ∈ [10^-3, 10^3]` and `log(FC) ∈ [-3, 3]`, the top-3 SVD components capture **`99.987%`** of normalized cone variance, and the 3 fixed target-free nodes `c ∈ {0.05, 1.0, 20.0}` achieve **`min R^2 = 0.99919`** (`median R^2 = 0.99985`).

## 2. Out-of-Sample Split-Half Selection (`SA-8`) Fails on 12/12 Human Cell Lines Under `N_1` — and Collapses to Zero Under `N_5`

| Cell Line (`Grp`, `N`) | `N1 E1(all)` | `N1 E1(top50_obs)` | `N1 E1(SA-8 split)` `[97.5% CI]` | `N1 ΔE1(SA-8 − all)` `[97.5% CI]` | `N1 Perm100 SA-8 (Med / Q97.5)` | `N5 E1(SA-8)` | `N5 Perm100 SA-8 (Med / Q97.5)` | `N5 ΔE1(SA-8 − all)` |
|---|---:|---:|---|---|---|---:|---|---:|
| `nadig_hepg2` (R9, 985) | -0.0455 | +0.2537 | +0.2234 [+0.2047, +0.2440] | +0.2689 [+0.2497, +0.2892] | +0.1161 / +0.1770 | -0.0001 | -0.0023 / +0.0375 | -0.0002 |
| `nadig_jurkat` (R9, 1180) | -0.0789 | +0.1248 | +0.0660 [+0.0467, +0.0858] | +0.1449 [+0.1276, +0.1633] | +0.0289 / +0.0942 | -0.0000 | +0.0006 / +0.0305 | -0.0000 |
| `nourreddine_kolf21j` (R9, 20) | -0.0381 | +0.1097 | +0.1564 [-0.0168, +0.3279] | +0.1945 [+0.0316, +0.3633] | +0.0932 / +0.1875 | +0.0000 | -0.0027 / +0.0321 | +0.0000 |
| `replogle_k562_essential` (R9, 1044) | -0.0436 | +0.1904 | +0.1736 [+0.1530, +0.1960] | +0.2172 [+0.1962, +0.2387] | +0.0986 / +0.1441 | -0.0002 | -0.0004 / +0.0118 | -0.0001 |
| `replogle_rpe1` (R9, 1300) | -0.0932 | +0.1602 | +0.1366 [+0.1215, +0.1527] | +0.2298 [+0.2153, +0.2460] | +0.0663 / +0.1071 | -0.0001 | +0.0003 / +0.0295 | -0.0001 |
| `vcc2025_h1` (R9, 91) | +0.0441 | +0.2866 | +0.2859 [+0.2317, +0.3420] | +0.2417 [+0.1933, +0.2909] | +0.1916 / +0.2295 | -0.0000 | -0.0011 / +0.0150 | -0.0000 |
| `dld1` (R10, 93) | -0.0652 | +0.1388 | +0.0957 [+0.0553, +0.1394] | +0.1610 [+0.1187, +0.2084] | +0.0511 / +0.0837 | -0.0000 | +0.0009 / +0.0175 | -0.0001 |
| `cd4t_rest` (R10, 398) | +0.0013 | +0.0875 | +0.0949 [+0.0629, +0.1291] | +0.0936 [+0.0607, +0.1277] | +0.0673 / +0.0916 | +0.0000 | +0.0009 / +0.0195 | +0.0000 |
| `calu3` (R10, 9) | +0.0201 | +0.0707 | +0.1019 [-0.0182, +0.1785] | +0.0818 [-0.0497, +0.1559] | +0.0522 / +0.1600 | -0.0000 | +0.0023 / +0.0661 | -0.0000 |
| `cd4t_stim48` (R11, 656) | -0.0276 | +0.0520 | +0.0529 [+0.0223, +0.0842] | +0.0805 [+0.0500, +0.1123] | +0.0382 / +0.0650 | +0.0000 | +0.0011 / +0.0233 | -0.0000 |
| `feng_ipsc` (R11, 111) | +0.0479 | +0.2316 | +0.2287 [+0.1775, +0.2817] | +0.1808 [+0.1370, +0.2272] | +0.1338 / +0.1688 | -0.0000 | -0.0005 / +0.0261 | +0.0000 |
| `telohaec` (R11, 9) | -0.0528 | +0.1726 | +0.2022 [-0.0440, +0.4675] | +0.2550 [-0.0016, +0.5973] | +0.0925 / +0.2340 | -0.0000 | -0.0007 / +0.0497 | -0.0000 |

- **Summary (`SA-8` Split-Half Audit)**: Under legacy `N_1`, out-of-sample split-half selection yields positive `E1(SA-8) − E1(all)` on **`12/12`** lines (`97.5% CI > 0` on **`10/12`** lines; `E1(SA-8) >= E1(top50_obs) − 0.003` on **`7/12`** lines). Under our target-free universal cone **`N_5`**, `max |E1(SA-8)| = 0.0002` across all `12` lines!

## 3. Universal Target-Free Null-Cone Cure (`N_5`) on `all_nontarget` Across All 12 Human Cell Lines

| Cell Line (`Grp`, `N`) | `N5 FCT` `[97.5% CI]` | `η(FCT)` | `N5 SGM (b_S)` | `N5 E1` | `N5 Perm100 (Med / Q97.5)` | `N5 MaxMis CI97.5_hi` (`Line`) | `N5 FCT − E1` `[97.5% CI]` | `FCT_lo > MaxMis_hi & Perm975` |
|---|---|---:|---:|---:|---|---|---|---|
| `nadig_hepg2` (R9, 985) | +0.1899 [+0.1840, +0.1964] | 0.934 | -0.0000 | +0.0000 | -0.0004 / +0.0076 | +0.0616 (`replogle_k562_essential`) | +0.1899 [+0.1840, +0.1964] | PASS |
| `nadig_jurkat` (R9, 1180) | +0.1314 [+0.1240, +0.1391] | 0.850 | -0.0000 | +0.0000 | +0.0005 / +0.0088 | +0.0995 (`replogle_k562_essential`) | +0.1314 [+0.1239, +0.1391] | PASS |
| `nourreddine_kolf21j` (R9, 20) | +0.1116 [+0.0498, +0.1691] | 0.882 | -0.0000 | -0.0000 | +0.0001 / +0.0091 | +0.0837 (`replogle_k562_essential`) | +0.1117 [+0.0498, +0.1691] | PERM_ONLY |
| `replogle_k562_essential` (R9, 1044) | +0.1832 [+0.1755, +0.1911] | 0.925 | -0.0000 | -0.0001 | +0.0003 / +0.0081 | +0.0648 (`nadig_jurkat`) | +0.1833 [+0.1756, +0.1911] | PASS |
| `replogle_rpe1` (R9, 1300) | +0.1674 [+0.1606, +0.1746] | 0.928 | +0.0000 | +0.0000 | -0.0008 / +0.0074 | +0.0850 (`nadig_jurkat`) | +0.1674 [+0.1606, +0.1745] | PASS |
| `vcc2025_h1` (R9, 91) | +0.1040 [+0.0855, +0.1239] | 0.746 | +0.0000 | +0.0000 | -0.0004 / +0.0040 | +0.0420 (`nourreddine_kolf21j`) | +0.1040 [+0.0855, +0.1239] | PASS |
| `dld1` (R10, 93) | +0.1480 [+0.1270, +0.1706] | 0.923 | -0.0000 | +0.0000 | -0.0001 / +0.0071 | +0.0670 (`replogle_rpe1`) | +0.1479 [+0.1270, +0.1705] | PASS |
| `cd4t_rest` (R10, 398) | +0.0567 [+0.0394, +0.0740] | 0.745 | +0.0000 | +0.0000 | +0.0005 / +0.0035 | +0.0378 (`vcc2025_h1`) | +0.0567 [+0.0394, +0.0740] | PASS |
| `calu3` (R10, 9) | +0.0854 [+0.0289, +0.1619] | 0.911 | -0.0000 | +0.0000 | -0.0006 / +0.0101 | +0.0780 (`replogle_rpe1`) | +0.0854 [+0.0289, +0.1619] | PERM_ONLY |
| `cd4t_stim48` (R11, 656) | +0.0734 [+0.0596, +0.0879] | 0.748 | +0.0000 | +0.0000 | -0.0001 / +0.0035 | +0.0766 (`vcc2025_h1`) | +0.0734 [+0.0596, +0.0879] | PERM_ONLY |
| `feng_ipsc` (R11, 111) | +0.1298 [+0.1153, +0.1446] | 0.902 | -0.0000 | -0.0000 | +0.0006 / +0.0070 | +0.0778 (`nourreddine_kolf21j`) | +0.1298 [+0.1153, +0.1446] | PASS |
| `telohaec` (R11, 9) | +0.1773 [+0.0811, +0.2819] | 0.801 | +0.0000 | -0.0000 | -0.0003 / +0.0156 | +0.1589 (`nadig_hepg2`) | +0.1773 [+0.0811, +0.2819] | PERM_ONLY |

- **Summary (`N_5`, `all_nontarget`)**: `FCT` exceeds `E1`, `Perm100`, and `MaxMis` point estimates on **`12/12`** lines; paired `FCT − E1` `97.5% CI > 0` on **`12/12`** lines; `FCT CI97.5_lo > Perm100 Q97.5` on **`12/12`** lines; and strict `FCT CI97.5_lo > max(MaxMis CI97.5_hi, Perm100 Q97.5)` holds on **`8/12`** lines.

## 4. Complete 48-Evaluation Artifact Annihilation Summary (`12 lines × 4 masks`)

| Artifact Class | Legacy `N_1` Range (`min .. max`) | `N_4` Semi-Partial Range | **`N_5` Semi-Partial Range** | `N_5` Median `η` |
|---|---|---|---|---:|
| `GlobalScale` (`R18`) | `-0.0007` .. `+0.0004` | `-0.0000` .. `+0.0000` | **`-0.0000` .. `+0.0000`** | `0.0000` |
| `Source_Grand_Mean` (`b_S`) | `-0.2795` .. `+0.4478` | `-0.0798` .. `+0.1865` | **`-0.0000` .. `+0.0000`** | `0.0000` |
| `E1_ConstMed` | `-0.0948` .. `+0.2866` | `-0.0002` .. `+0.0002` | **`-0.0002` .. `+0.0001`** | `0.0012` |
| `E2_ConstP10` | `-0.0957` .. `+0.2587` | `-0.0001` .. `+0.0000` | **`-0.0001` .. `+0.0000`** | `0.0006` |
| `E2_ConstP90` | `-0.0995` .. `+0.3004` | `-0.0004` .. `+0.0002` | **`-0.0003` .. `+0.0002`** | `0.0025` |
| `E6_LogFCBar` | `-0.1033` .. `+0.2986` | `-0.0003` .. `+0.0011` | **`-0.0003` .. `+0.0011`** | `0.0054` |
| `Perm100 Median` | `-0.0616` .. `+0.1916` | `-0.0062` .. `+0.0054` | **`-0.0043` .. `+0.0056`** | `0.4731` |

## 5. Direct Q5a Cross-Line Conservation Test (`33` Human Cell-Line Pairs, `>=25` Shared Reproducible Targets)

| Gene Universe | Pairs | Raw `log1p` | Raw `logFC` | Raw `FCT` (`FCT>log1p`, Wilcoxon `p`) | Cal `log1p` | Cal `logFC` (`logFC>log1p`, Wilcoxon `p`) | Cal `FCT` (`FCT>log1p`, Wilcoxon `p`) |
|---|---:|---:|---:|---|---:|---|---|
| `all` | 33 | +0.1915 | +0.1394 | +0.1979 (`28/33`, `p=3.41e-06`) | +0.3220 | +0.3216 (`16/33`, `p=9.44e-01`) | +0.3305 (`28/33`, `p=1.28e-05`) |
| `expr01` | 33 | +0.1938 | +0.1503 | +0.2024 (`32/33`, `p=2.33e-09`) | +0.3257 | +0.3333 (`21/33`, `p=9.33e-02`) | +0.3396 (`32/33`, `p=2.94e-07`) |
| `expr05` | 33 | +0.2221 | +0.2092 | +0.2290 (`33/33`, `p=2.33e-10`) | +0.3466 | +0.3566 (`25/33`, `p=1.47e-03`) | +0.3577 (`33/33`, `p=2.33e-10`) |

