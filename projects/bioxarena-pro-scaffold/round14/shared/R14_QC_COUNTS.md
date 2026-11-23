# Round 14 QC & Eligibility Counts (`PREREG_G` §4.1, §6, §7)

- **Blinding Declaration**: Strictly **counts-only** audit executed after [PREREG_G.md](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round14/prereg/PREREG_G.md) (`commit d34560d`), [R14_GO.txt](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round14/status/R14_GO.txt), and [R14_LINES_AUTHORIZED.txt](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round14/status/R14_LINES_AUTHORIZED.txt), using unmodified [dataset_readers_r14.py](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round14/src/dataset_readers_r14.py) (`sha256: e0b0c461e33950ee22897e8c36c69fa15af4fb8e46e68e6d0fe91c73d53d4104`) and [verify_r10_lines.py](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round10/shared/scripts/verify_r10_lines.py) via `python3 -I`. Zero transfer predictions, correlations, `r_orth_cal`, slopes, residuals, `FCT`, `E1_ConstMed`, `N_5`, or subgroup metric evaluations were computed, and zero code files were modified.
- **Line Cache Gate**: [round14/shared/lines/a549_myllia](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round14/shared/lines/a549_myllia) verified complete (`ctrl_cells.npz`, `genes.txt`, `line_metadata.json`, `pseudobulk_delta.npz`, `targets.txt`), zero `_tmp_*` files, and all contract checks passed in `verify_and_count_line`.

---

## 1. `PREREG_G` §4.1 Feasibility, Integer Raw Slot & Per-Channel Cell QC

| Line | Role (`PREREG_G` §2) | Feasibility (`§4.1`) | Slot Used | Guide Call Mode | Raw Cells (`obs`) | `gene_desert` Excluded (`ctrl_T`) | RNA QC / Multi-Zero / Unmapped Excluded | Kept Single-sgRNA+NTC Cells | Targeting Cells | Total `ctrl_NT` NTC (`20` guides) | NTC Half A (`ch_1`) / Half B (`ch_2`) | Genes (`raw -> fin -> final`) |
|---|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| **`a549_myllia`** | Sole confirmatory candidate (`max K_conf = 1`) | **`PASS`** | **`/X`** (`== /raw/X`, `csr_nonneg_int`) | `obs_author_single_guide` | `53,511` | **`2,063`** | `0` / `0` / `0` | **`51,448`** | **`49,355`** | **`2,093`** | **`1,033` / `1,060`** | `18480 -> 15456 -> 15456` |

### Per-Channel Breakdown (`PREREG_G` §5 Disjoint Split-Half & NTC Subtraction Units)

| Channel Tag | File | Slot & Median Raw UMI | Assigned Half & Batch | Raw `obs` | `gene_desert` (`ctrl_T`) Excluded | Kept Cells | Targeting Cells | `ctrl_NT` NTC Kept | `< 15` NTC Fallback? |
|---|---|---|---|---:|---:|---:|---:|---:|---|
| **`a549_myllia_ch_1`** | `GSE337804_A549_ch_1.RDS.h5ad` | `/X` `csr_matrix` `[25702, 18711]`, `float64`, `median=13,578.0`, `n_nonint=0` | `Half A (1)` / `ch_1` | `25,702` | **`988`** | `24,714` | `23,681` | **`1,033`** | NO (`1,033 >= 15`) |
| **`a549_myllia_ch_2`** | `GSE337804_A549_ch_2.RDS.h5ad` | `/X` `csr_matrix` `[27809, 18808]`, `float64`, `median=13,717.0`, `n_nonint=0` | `Half B (2)` / `ch_2` | `27,809` | **`1,075`** | `26,734` | `25,674` | **`1,060`** | NO (`1,060 >= 15`) |
| **Total** | `2 .RDS.h5ad` | Exact match w/ `GSE337804_sgRNA_lib.csv.gz` | `2` batches (`ch_1`, `ch_2`) | **`53,511`** | **`2,063`** | **`51,448`** | **`49,355`** | **`2,093`** | `0` fallback cells |

- **Gene Standardizing Counts (`line_metadata.json`)**: `n_genes_multifile_dropped = 231`, `n_genes_unmapped_to_approved_hgnc = 379`, `n_genes_raw = 18,480`, `n_genes_finite_all_batches = 15,456`, `n_genes_empty_symbol_dropped = 0`, `n_genes_duplicate_symbol_dropped = 0`, `n_genes_low_expr_cap_dropped = 0`, `n_genes_final = 15,456`, **`|G_eval(a549_myllia)| = 6,126`** (`5,966` raw) **`>= 3,000`** (`PASS`).

---

## 2. `PREREG_G` §6 Confirmatory Count Gate & Addendum-B `P5` Attrition (`K_CONF = 0`)

- **Target Attrition Funnel (`218` -> `32`)**:
  - Library canonical targets (`n_raw_pert_groups`): **`218`**
  - Retained deduplicated targets (`n >= 10` cells): **`212`** (`6` dropped at `n < 10`), all **`212`** `half_valid` (`>= 2` cells in both `ch_1` and `ch_2`); `206` targets had `n >= 25` cells and were tested under Addendum-B.
  - **(1) Reproducible on `a549_myllia`** (`n >= 25`, `AD >= 50`, `Energy p < 0.001`): **`81`** (`81` raw)
  - **(2) Shared reproducible w/ `r9_6`** (`u6_rep = 2,117`): **`34`** (`34` raw)
  - **(3) `n_eval_targets`** (`self-gene in G_eval & half_valid`): **`32`** (`32` raw; `0` HGNC-bridged additions; the `2` `r9_6`-shared reproducible targets absent from `G_eval` are `LDB1` and `UBE3D`, which are in `genes.txt` but not in `R9_COMMON_6148`).

| `PREREG_G` §6 Gate Criterion | Required Threshold | Observed (`a549_myllia`) | Status |
|---|---|---:|---|
| 1. Technical feasibility & `/X` raw integers (§4.1) | `PASS` | `PASS` (`/X == /raw/X`, `n_nonint = 0`) | **`PASS`** |
| 2. **`n_eval_targets >= 50`** | `>= 50` | **`32`** (`32 < 50`, `29` of `61` dropped) | **`FAIL`** |
| 3. **`\|G_eval(a549_myllia)\| >= 3,000`** | `>= 3,000` | **`6,126`** (`5,966` raw) | **`PASS`** |
| 4. **`per-half NTC >= 15`** (`ctrl_NT` only) | `ch_1 >= 15` AND `ch_2 >= 15` | **`1,033`** / **`1,060`** | **`PASS`** |

- **Mechanical Verdict (`round14/status/CONFIRMATORY_LINES_R14.txt`)**: **`K_CONF = 0`** — **`a549_myllia` classified `SECONDARY_ONLY`**. Per `PREREG_G §3.2` and `§6`, Round 14 executes **zero confirmatory adjudication** (`不做确认性判定`); all downstream evaluation must be reported strictly descriptively under `PREREG_E §5` and `PREREG_G §7`, leaving Round 12 `a549_ko` conclusions unchanged.

---

## 3. `PREREG_G` §7 Pre-Registered Subgroups (`R12-Overlap` vs `Novel`, Counts Only)

Partitioning the **`32`** `eval_targets` strictly against pre-locked `PREREG_G §7` candidate lists and [round12/shared/lines/a549_ko/targets.txt](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/shared/lines/a549_ko/targets.txt):

| Subgroup (`PREREG_G` §7) | Candidate Pool (`u6_rep ∩ R9C`) | Surviving `eval_targets` | Survival Rate | Exact Surviving Symbol List (Alphabetical) |
|---|---:|---:|---:|---|
| **`overlap_r12_a549`** (present in R12 `a549_ko`) | `15` | **`4`** | `26.7%` (`4 / 15`) | `CREBBP`, `HDAC3`, `KDM6A`, `PTEN` |
| **`novel_not_in_r12_a549`** (absent from R12 `a549_ko`) | `46` | **`28`** | `60.9%` (`28 / 46`) | `BRPF1`, `CBX3`, `DHX29`, `DHX36`, `EHMT1`, `EP300`, `HDAC1`, `HSD17B10`, `INTS14`, `IRF2BP2`, `KAT6A`, `KDM1A`, `KDM3B`, `MSL1`, `MTPAP`, `PIK3R4`, `PRORP`, `RAF1`, `RPUSD3`, `RPUSD4`, `SHOC2`, `SMARCA5`, `SRFBP1`, `SUPT20H`, `SUPV3L1`, `TFPT`, `TRMT10C`, `UBA5` |
| **Total `eval_targets`** | **`61`** | **`32`** | `52.5%` (`32 / 61`) | Disjoint partition (`4 + 28 = 32`) |

### Counts-Only Attrition Summary of Dropped Candidates (`29` of `61`)

- **`overlap_r12_a549` dropped (`11` of `15`)**:
  - `n < 10` (absent from `targets.txt`): `BMS1`, `DHDDS`, `SEC61A1` (`3`)
  - Kept in `targets.txt` with `n < 25` (below Addendum-B min cells): `GMPPB` (`n=11`, `7/4`), `CMTR1` (`n=24`, `15/9`) (`2`)
  - `n >= 25` failing `AD >= 50` and/or `Energy p < 0.001`: `EIF4G1` (`n=26`, `ad=7`, `p=0.00333`), `RIOK2` (`n=29`, `ad=2`, `p=0.0891`), `SPCS3` (`n=33`, `ad=11`, `p=0.0199`), `SLC35B1` (`n=76`, `ad=6`, `p=0.000500`), `SMC3` (`n=130`, `ad=74`, `p=0.0396`), `CNOT10` (`n=199`, `ad=135`, `p=0.0297`) (`6`)
- **`novel_not_in_r12_a549` dropped (`18` of `46`)**:
  - `n < 10` (absent from `targets.txt`): `DOT1L`, `UTP6` (`2`)
  - Kept with `n < 25`: `UBA3` (`n=12`, `5/7`), `RRP9` (`n=14`, `6/8`), `DDX1` (`n=18`, `9/9`) (`3`)
  - `n >= 25` failing `AD >= 50` and/or `Energy p < 0.001`: `TFAM` (`n=29`, `ad=23`, `p=0.000500`), `CTCF` (`n=40`, `ad=3`, `p=0.0990`), `RPP14` (`n=41`, `ad=0`, `p=0.9604`), `FASTKD5` (`n=48`, `ad=20`, `p=0.000500`), `SIN3B` (`n=75`, `ad=49`, `p=0.00997`), `ZNF598` (`n=98`, `ad=0`, `p=0.4158`), `SRF` (`n=99`, `ad=2`, `p=0.7030`), `EWSR1` (`n=145`, `ad=8`, `p=0.00998`), `CALR` (`n=165`, `ad=15`, `p=0.0594`), `KDM2B` (`n=202`, `ad=67`, `p=0.00235`), `NEDD4` (`n=218`, `ad=0`, `p=0.8218`), `KEAP1` (`n=254`, `ad=1`, `p=0.8416`), `GSK3B` (`n=283`, `ad=2`, `p=0.5743`) (`13`)

---

## 4. Adapter Audit (`round14/shared/lines/a549_myllia/line_metadata.json`)

- **`a549_myllia` (`ch_1_HalfA_vs_ch_2_HalfB`, `x_layer="X"`)**: `{"line": "a549_myllia", "n_h5ad_files": 2, "channel_slots": {"ch_1": {"x_slot_used": "X", "selected_median_total_umi": 13578.0, "int_slots": ["X", "raw/X"]}, "ch_2": {"x_slot_used": "X", "selected_median_total_umi": 13717.0, "int_slots": ["X", "raw/X"]}}, "sgrna_lib_source": "sgrna_lib:GSE337804_sgRNA_lib.csv.gz", "guide_call_mode": "obs_author_single_guide", "min_umi": 300, "min_ngenes": 200, "max_pct_mt": null, "n_rna_qc_excluded": 0, "n_multi_or_zero_guide_excluded": 0, "n_safe_targeting_excluded": 2063, "n_unmapped_guide_excluded": 0, "n_qc_kept_cells": 51448, "n_kept_ntc": 2093, "n_kept_half_a_ntc": 1033, "n_kept_half_b_ntc": 1060, "n_low_ntc_channel_fallback_cells": 0, "half_rule_used": "ch_1_HalfA_vs_ch_2_HalfB", "channels_used": ["ch_1", "ch_2"], "median_guide_threshold": null}`
- **`CONFIRMATORY_LINES_R14.txt` (`K_CONF = 0`)**: `[]` (no passing confirmatory line; `SECONDARY_ONLY` descriptive reporting only).
