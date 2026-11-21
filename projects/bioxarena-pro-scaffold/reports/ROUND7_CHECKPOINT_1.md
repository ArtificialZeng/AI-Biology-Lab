# Round 7 Checkpoint 1: Pre-Registered Unseen External Confirmation (Nadig et al. 2025, HepG2 + Jurkat)

- **Date**: 2026-11-21
- **Pre-Registration Chain (Strictly Prior to Reading Any Expression Matrix Values)**:
  - `c5294c0`: Initial `R7-0` pre-registration (`reports/ROUND7_PREREGISTRATION.md`)
  - `87f4bc2`: Addendum A errata (code-authoritative `SEEDS_10`, `make_strict_union_split`, `S3` naming, metadata-only lock rule)
  - `b9fe2bb`: Addendum B (official GEO `GSE264667` provenance/SHA-256, exact Replogle 2022 gem_group normalization, locked Anderson-Darling + Energy-distance `reproducible` criteria, top-5,000 `control_mean` readout selection)
- **Frozen Solver**: `src/xatlas_panhuman_solver.py` @ `3609811` (zero constants modified; external dataset strictly evaluation-only).

---

## 1. Pre-Registered Hypothesis Adjudication (`reports/r71_nadig_confirmation_summary.json`)

| Hypothesis | Priority | Metric & Comparison | Pre-Registered Decision Criterion | Observed Result (10 Strict Union Splits) | Verdict |
|---|---|---|---|---|---|
| **H1** | **Primary** | `disc_raw` (`S5 − AE`) | Mean `Δ > 0` **AND** `≥ 7/10` splits row-level 95% CI `> 0` | Mean `Δ = +0.0487` (`0.7152` vs `0.6665`), mean 95% CI `[+0.0161, +0.0812]`, **`9/10` splits CI `> 0`** | **确认 (CONFIRMED)** |
| **H2** | **Primary** | `delta_r` (`S5 − AE`) | Mean `Δ < 0` (metric inversion after perturbation centroid removal) | Mean `Δ = −0.0195` (`0.1502` vs `0.1698`), mean 95% CI `[−0.0446, +0.0056]`, `9/10` splits `< 0`, `4/10` splits CI `< 0`, `0/10` CI `> 0` | **确认 (CONFIRMED)** |
| **Main Thesis** | **Joint Primary** | **H1 AND H2 simultaneously** | **BOTH H1 AND H2 must pass** | **H1 = Pass AND H2 = Pass** | **确认 (CONFIRMED)** |
| **H3** | **Secondary** | `v5_orth_r` (`S5 − AE`) | `≥ 7/10` splits row-level 95% CI `> 0` | Mean `Δ = +0.0595` (`0.1741` vs `0.1146`), mean 95% CI `[+0.0437, +0.0753]`, **`10/10` splits CI `> 0`** | **确认 (CONFIRMED)** |

---

## 2. Official Data Provenance, Preprocessing Audit & Cache SHA-256

- **Official Release**: NCBI GEO `GSE264667` (Nadig, Replogle, O'Connor, Weissman et al., *Nature Genetics* 2025, PMID `40259084`; public NIH GEO open-access release). No author-published `normalized_bulk_01.h5ad` exists on GEO, Figshare, or Zenodo; preprocessing strictly reproduces Replogle et al. (*Cell* 2022) per-`gem_group` non-targeting control z-scoring via `src/prepare_nadig2025_cache.py`.
- **Unit Check**: Vectorized Scholz-Stephens two-sample Anderson-Darling test verified against `scipy.stats.anderson_ksamp` (`max_abs_diff < 1e-5`).
- **20-Perturbation Speed Benchmarks**: HepG2 `6.84s / 20 perts` (actual total `902.9s`); Jurkat `9.31s / 20 perts` (actual total `1,270.7s` using contiguous sequential chunk HDF5 reads).

| Cell Line (Slot) | Official Input File & SHA-256 | Raw Shape (`cells × genes`) | Deduplicated ENSG Targets (`rep / nonrep`) | Cache File SHA-256 (`cache_manifest.json`) |
|---|---|---|---|---|
| **HepG2** (`HCT116`) | `GSE264667_hepg2_raw_singlecell_01.h5ad`<br>`e1ad7c3c5a201c861a207a858aa7e59f5e6ac1955674c415f7de0d1dadadb52e` | `145,473 × 9,624` (`56` gems, `4,976` NT) | `2,343` (`1,000` / `1,343`) | `HepG2_X.npy`: `97e164ed92e1486b86001f9152c8978942121af8030fc1b60bec85e18783812f`<br>`HepG2_meta.csv`: `51786d722f1c3eb0cbc4a7403398c2ce8b5c302a3edec67b953b68f8e7eaf923`<br>`HepG2_genes.csv`: `af43ce86c0d1093ad0274acf0286080388319eeb77ed66c7d15243ea221fd85c` |
| **Jurkat** (`HEK293T`) | `GSE264667_jurkat_raw_singlecell_01.h5ad`<br>`ffbe15f2c8f7ffcfd7b0ba9e6937d4ebc2d03b0179fa8234648a59bcb82c04a3` | `262,956 × 8,882` (`55` gems, `12,013` NT) | `2,350` (`1,262` / `1,088`) | `Jurkat_X.npy`: `4a85657883033953923b440e58a795d599777fb428d056ac1bdabeefa6729d32`<br>`Jurkat_meta.csv`: `fa828370d93b1fea0bf18260c1af2b75fac07e8f1a6d3b6bc37181b1cc797509`<br>`Jurkat_genes.csv`: `200dca8696f89107da4354c29eb881a04041e420cf8827ba79aa129627c27304` |

- Shared selected readout genes across HepG2 & Jurkat: **`3,644` / `5,000`** Ensembl genes.
- Per-split strict simultaneous union holdout size (`SEEDS_10`): **`147` HepG2 + `191` Jurkat = `338` held-out reproducible rows** (`305–322` union targets excluded simultaneously from both lines' training sets).

---

## 3. Complete Main Benchmark Table Across All Methods & All Metrics (`reports/r71_nadig_5metric_combined.csv`)

*(Mean ± SD across `SEEDS_10 = [42, 0, 1, 2, 3, 7, 11, 99, 2024, 31337]`; macro-averaged `0.5 × HepG2 + 0.5 × Jurkat`; **all unfavorable metrics reported in main table**.)*

| Method | Role | `v5` (Orth. Residual) | `delta_r` (= Train-Centroid r) | `de50_r` (Top-50 DE r) | `de50_sign` (Top-50 Sign) | `disc_raw` (Primary H1) | `disc_v5` | `discrim_delta_pct` (Target-Masked) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| `Constant_Mean` | Trivial | `0.0000 ± 0.0000` | `-0.0001 ± 0.0013` | `0.0057 ± 0.0079` | `0.5000 ± 0.0043` | `0.5000 ± 0.0000` | `0.5000 ± 0.0000` | `0.5000 ± 0.0002` |
| `kNN_30` | Non-param | `0.1251 ± 0.0100` | `0.1616 ± 0.0100` | `0.2021 ± 0.0172` | `0.6195 ± 0.0087` | `0.6502 ± 0.0126` | `0.6512 ± 0.0165` | `0.6503 ± 0.0126` |
| **`AE`** | **Linear Comparator** | `0.1146 ± 0.0110` | **`0.1698 ± 0.0079`** | `0.2492 ± 0.0198` | **`0.6457 ± 0.0097`** | `0.6665 ± 0.0149` | `0.6509 ± 0.0145` | `0.6662 ± 0.0150` |
| `S0` | Linear Ridge | `0.1042 ± 0.0098` | `0.1610 ± 0.0067` | **`0.2511 ± 0.0188`** | `0.6451 ± 0.0094` | `0.6709 ± 0.0146` | `0.6550 ± 0.0130` | `0.6706 ± 0.0147` |
| `Simplified_ResMLP` | Simplified DL | `0.1292 ± 0.0107` | `0.1100 ± 0.0090` | `0.1817 ± 0.0253` | `0.6019 ± 0.0085` | `0.6622 ± 0.0146` | `0.6591 ± 0.0146` | `0.6627 ± 0.0144` |
| `Simplified_GNN` | Simplified DL | `0.1432 ± 0.0077` | `0.1208 ± 0.0059` | `0.1878 ± 0.0200` | `0.6090 ± 0.0073` | `0.6579 ± 0.0155` | `0.6564 ± 0.0154` | `0.6584 ± 0.0154` |
| `S1` (HybridTOM) | Ablation | `0.1506 ± 0.0079` | `0.1295 ± 0.0065` | `0.1985 ± 0.0137` | `0.6155 ± 0.0082` | `0.6991 ± 0.0100` | `0.6998 ± 0.0097` | `0.7001 ± 0.0100` |
| `S2` (+Deflation) | Ablation | `0.1444 ± 0.0080` | `0.1243 ± 0.0063` | `0.1987 ± 0.0120` | `0.6093 ± 0.0098` | `0.6963 ± 0.0091` | `0.6965 ± 0.0092` | `0.6972 ± 0.0091` |
| `S3` (+Cross-Graph) | Ablation | `0.1544 ± 0.0101` | `0.1330 ± 0.0084` | `0.2140 ± 0.0105` | `0.6180 ± 0.0103` | `0.7091 ± 0.0114` | `0.7095 ± 0.0118` | `0.7099 ± 0.0115` |
| `S4` (+NonRep) | Ablation | `0.1584 ± 0.0098` | `0.1367 ± 0.0082` | `0.2231 ± 0.0107` | `0.6229 ± 0.0092` | `0.7130 ± 0.0113` | `0.7139 ± 0.0121` | `0.7139 ± 0.0114` |
| **`S5`** | **PRIMARY (Frozen)** | **`0.1741 ± 0.0121`** | `0.1502 ± 0.0105` | `0.2265 ± 0.0175` | `0.6293 ± 0.0128` | **`0.7152 ± 0.0129`** | **`0.7157 ± 0.0132`** | **`0.7160 ± 0.0128`** |
| `S6` (+Reciprocal) | Ablation | `0.1756 ± 0.0125` | `0.1517 ± 0.0110` | `0.2384 ± 0.0195` | `0.6307 ± 0.0123` | **`0.7184 ± 0.0135`** | **`0.7190 ± 0.0135`** | **`0.7205 ± 0.0133`** |
| **`S7`** | **SECONDARY (Frozen)** | **`0.1871 ± 0.0150`** | `0.1617 ± 0.0128` | `0.2242 ± 0.0212` | `0.6277 ± 0.0134` | `0.6992 ± 0.0156` | `0.6993 ± 0.0160` | `0.7008 ± 0.0155` |
| `S7_no_S6_exploratory` | **Exploratory Only** | `0.1901 ± 0.0150` | `0.1644 ± 0.0127` | `0.2164 ± 0.0210` | `0.6298 ± 0.0142` | `0.6983 ± 0.0147` | `0.6984 ± 0.0154` | `0.6992 ± 0.0147` |

---

## 4. Row-Level Stratified Paired Bootstrap 95% CIs (`reports/r71_nadig_row_ci.csv`, 4,000 resamples)

| Comparison | `v5` Mean Δ [Avg 95% CI] (`CI>0`) | `delta_r` Mean Δ [Avg 95% CI] (`CI<0`) | `disc_raw` Mean Δ [Avg 95% CI] (`CI>0`) | `disc_v5` Mean Δ [Avg 95% CI] (`CI>0`) |
|---|---|---|---|---|
| **`S5 − AE` (Primary)** | **`+0.0595`** `[+0.0437, +0.0753]` (**10/10**) | **`−0.0195`** `[−0.0446, +0.0056]` (4/10 `<0`) | **`+0.0487`** `[+0.0161, +0.0812]` (**9/10**) | **`+0.0647`** `[+0.0382, +0.0917]` (**10/10**) |
| **`S7 − AE` (Secondary)** | **`+0.0725`** `[+0.0523, +0.0928]` (**10/10**) | `−0.0080` `[−0.0355, +0.0199]` (1/10 `<0`) | `+0.0327` `[−0.0005, +0.0661]` (4/10) | `+0.0484` `[+0.0207, +0.0765]` (**9/10**) |
| `S7_no_S6 − AE` (Expl.) | `+0.0755` `[+0.0556, +0.0951]` (**10/10**) | `−0.0054` `[−0.0326, +0.0222]` (0/10 `<0`) | `+0.0319` `[−0.0012, +0.0649]` (3/10) | `+0.0475` `[+0.0203, +0.0751]` (**9/10**) |
| `S5 − S0` | `+0.0699` `[+0.0545, +0.0852]` (**10/10**) | `−0.0107` `[−0.0341, +0.0127]` (0/10 `<0`) | `+0.0443` `[+0.0120, +0.0768]` (6/10) | `+0.0606` `[+0.0341, +0.0875]` (**10/10**) |
| `S6 − S5` | `+0.0014` `[−0.0038, +0.0069]` (0/10) | `+0.0015` `[−0.0031, +0.0063]` (0/10) | `+0.0032` `[−0.0048, +0.0116]` (0/10) | `+0.0034` `[−0.0046, +0.0118]` (0/10) |
| **`S5 − Simplified_GNN`** | **`+0.0309`** `[+0.0130, +0.0490]` (**9/10**) | **`+0.0295`** `[+0.0133, +0.0459]` (**9/10**) | **`+0.0572`** `[+0.0295, +0.0850]` (**10/10**) | **`+0.0593`** `[+0.0323, +0.0865]` (**10/10**) |
| **`S7 − Simplified_GNN`** | **`+0.0438`** `[+0.0225, +0.0657]` (**10/10**) | **`+0.0410`** `[+0.0217, +0.0606]` (**10/10**) | **`+0.0412`** `[+0.0123, +0.0701]` (**8/10**) | **`+0.0429`** `[+0.0148, +0.0712]` (**8/10**) |
| `S5 − Simplified_ResMLP` | `+0.0449` `[+0.0288, +0.0613]` (**10/10**) | `+0.0402` `[+0.0257, +0.0551]` (**10/10**) | `+0.0529` `[+0.0260, +0.0801]` (**10/10**) | `+0.0565` `[+0.0302, +0.0834]` (**10/10**) |

### Key Replication Notes (Unfavorable & Structural Findings Disclosed)
1. **Metric Inversion Confirmed on a Third Independent Atlas**: Linear baselines (`AE`, `S0`) outperform `S5` and `S7` on `delta_r` (which is algebraically identical to training-centroid-centered Systema Pearson), `de50_r`, and `de50_sign`, while `S5` decisively outperforms `AE`, `S0`, `Simplified_ResMLP`, and `Simplified_GNN` on `v5` (`10/10 CI > 0` vs `AE`) and raw-space perturbation discrimination `disc_raw` (`9/10 CI > 0` vs `AE`; `10/10 CI > 0` vs `Simplified_GNN`).
2. **S6 & S7 Trade-off Replicated**: Just as on xatlas and Replogle 2022, `S6` achieves the highest `disc_raw` (`0.7184`) across all methods (with near-neutral `S6 − S5 = +0.0014` on `v5`, 0/10 CI > 0), while `S7` Wiener/manifold snapping boosts `v5` (`0.1741 → 0.1871`) at the cost of lowering `disc_raw` (`0.7184 → 0.6992`, reducing `S7 − AE` `disc_raw` to `4/10 CI > 0`).
3. **Unfavorable Detail**: Unlike `S5 − AE` (`9/10 CI > 0` on `disc_raw`), `S5 − S0` on `disc_raw` reaches `6/10` splits with CI `> 0` (`+0.0443`, `[+0.0120, +0.0768]`), though `disc_v5` is `10/10 CI > 0` (`+0.0606`), and `S5` trails `AE` on `de50_r` (`0.2265` vs `0.2492`) and `de50_sign` (`0.6293` vs `0.6457`).

CHECKPOINT_1_DONE
