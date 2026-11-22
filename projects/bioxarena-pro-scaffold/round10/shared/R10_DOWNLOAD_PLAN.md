# Round 10 (`R10-1`) Download, Standardization & QC Plan (`R10_DOWNLOAD_PLAN.md`)

**Agent**: A-T1 (Team A · Data Engineering)  
**Date**: 2026-11-22  
**Compliance**: Aligned verbatim with `PREREG_C_AUDIT_R1.md` (`P1–P7`) and `【Claude 主控 · 下载方案约束】`.  
**Hard Lock Rule**: **Zero expression data files will be downloaded until `/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round10/status/R10_GO.txt` exists.** After `R10_GO.txt`, ingestion runs Addendum-B QC and **reports counts only** before code freeze.  
**Host Verified Environment**: Python `3.13.5` (`h5py`, `scipy`, `pandas`, `pickle` installed; `Rscript=None`, `rpy2=No`, `rds2py=No`), **`36 GiB` (`38,654,705,664 B`) system RAM**, `2.5 TiB` free disk.

---

## 1. Pre-Download Audit & `60.000 GB` Budget Decision (`P2`–`P4`)

### 1.1 Active Downloading Portfolio (`34.080 GB <= 60.000 GB`, `25.920 GB` Margin)

| Role (`PREREG_C`) | Dataset / Line | Exact Files Selected | Download Size | Filtering / Single-sgRNA Rule (`P4`) |
|---|---|---|---|---|
| **Confirmatory 1** | **DLD-1** (`GSE337988`) | **Low-MOI ONLY**: `24` GSMs `GSM9864227..GSM9864250` (`SAM24482910_1..24`): `*.filtered_feature_bc_matrix.h5` (`5.788 GB`) + `*.crispr.umi_correct.h5` (`0.526 GB`) + `*.hashing.umi_correct.h5` (`0.182 GB`) + `NGS6475` guide CSV (`0.0002 GB`) *(pure `h5py`, avoids `.rds`)* | **`6.496 GB`** (`6,495,678,716 B`) | Low-MOI channels ONLY (`MOI ~0.3`); MULTI-seq singlets + **strictly 1 sgRNA** assigned (1 targeting or 1 NTC). |
| **Confirmatory 2** | **Primary CD4+ αβ T Cells** (`GSE314342` `Rest`) | **Lexicographically first 2 donors (`P3`)**:<br>• **`CE0006864` (`D4`, batch `CD4i_R2`)**: `24` `GSM9393851..74_*_D4_Rest_*.h5` = **`12.113 GB`** (`12,113,175,779 B`; `2,693,903` cells, `76,636` NTC single-sgRNA)<br>• **`CE0008162` (`D1`, batch `CD4i_R1`)**: `23` `GSM9393828..50_*_D1_Rest_*.h5` = **`14.428 GB`** (`14,427,676,272 B`; `3,074,496` cells, `76,543` NTC single-sgRNA)<br>• Paired S3 CSV metadata (`sgrna_library_metadata.suppl_table.csv` + `sample_metadata` + `QC_summaries`): **`0.010 GB`** | **`26.551 GB`** (`26,540,852,051 B` `.h5` + `0.010 GB` CSV) | `Rest` condition ONLY; **strictly 1 sgRNA** per cell (`guide_type` targeting/non-targeting, excluding `multi sgRNA` and `low_quality`). Combined: `5,768,399` cells (`153,179` NTC; `3,347,949` targeting). |
| **Secondary 1** | **WTC11 `i3Neurons`** (Tian et al. 2021) | Zenodo `13350497`: `TianKampmann2021_CRISPRi.h5ad` *(CRISPRi only; fallback 4 GEO `.h5` lanes `GSM4632022/23/25/26` = `0.316 GB` + `MOESM6` if `.h5ad` omits lane labels)* | **`0.289 GB`** (`289,389,479 B`) | Single sgRNA (`CRISPRi` only; exclude `CRISPRa`). |
| **Secondary 2** | **Calu-3** (Sunshine et al. 2023) | Zenodo `13350497`: `SunshineHein2023.h5ad` | **`0.744 GB`** (`743,672,579 B`) | **Single-gene element ONLY** (`195` single-gene elements / NTC pairs; exclude `22` dual-gene pairs). |
| **Total Active Download** | **4 Lines (2 Confirmatory + 2 Secondary)** | — | **`34.080 GB`** | Strictly within `<= 60.000 GB` (`25.920 GB` spare) and `< 4 GB` streaming RAM. |

*(Audit of all 4 donor IDs in `GSE314342` `sample_metadata.suppl_table.csv` sorted lexicographically: **`CE0006864` (`D4`), `CE0008162` (`D1`), `CE0008678` (`D3`), `CE0010866` (`D2`)** -> first two are **`CE0006864` (`D4`)** and **`CE0008162` (`D1`)**, matching `P3`.)*

---

### 1.2 Mandatory Pre-Download Exclusion Proof for H9 Cortical Neurons (`IGVFDS9825FOBP`, `P2`)

**Decision**: **H9 (`IGVFDS9825FOBP`) is EXCLUDED BEFORE DOWNLOAD** — zero H9 `.h5` / `.pkl` / `.h5ad` files will be downloaded. Two independent verified reasons mandate exclusion under `P2`:

1. **IGVF Small Metadata Audit + Pre-QC Single-sgRNA Scarcity (`P2` / `P4`)**:
   - **IGVF small-file inventory**: On IGVF `IGVFDS9825FOBP`, the ONLY non-matrix guide file is `IGVFFI7085MUTI.csv.gz` (`244,119 B`), which contains **only** the `10,691` sgRNA sequences and target `ENSG` IDs — **zero per-cell or per-target single-sgRNA cell counts**. All `72` `input_file_sets` (`LW558..LW581`) contain only raw `.fastq.gz` reads (`>300 GB`) and `.yaml` seqspec metadata.
   - **Exact conservation bound from `PMC13629240 Table S2` (`24` libraries `LW558..LW581`)**:
     - Across all `24` HTO-superloaded libraries (`MOI ~4`, mean `3.8` sgRNAs/cell): `Number of cells detected by 10X` = `1,689,060`; `Exactly 1 HTO` = `974,093`; `Exactly 1 HTO and > 0 sgRNA` = `958,228`; **`Exactly 1 HTO and exactly 1 sgRNA` = `53,037`** (`5.53%`).
     - Post-QC final dataset (`PMC13629240` Methods): `549,727` QC cells (`59,714` NTC cells = `10.86%`).
     - Even assuming proportional QC retention (`549,727 / 958,228 = 57.37%`) and zero sgRNA-outlier loss (`609` sgRNAs were removed as outliers), total post-QC single-sgRNA HTO-singlets across the **entire** dataset are only `~30,427` cells (`~27,100` targeting cells after `~10.9%` controls) spread across `1,408` targeted genes -> **mean `~19.2` single-sgRNA cells per target gene** (below the `>= 25` cell gate for typical targets). Even the unobtainable pre-QC count (`~47,200` targeting cells across `1,408` genes ≈ `33.5` pre-filter cells/gene) cannot verify without unpacking the `46.074 GB` `.pkl` that `>= 50` `union6`-shared targets survive `>= 25` post-QC single-sgRNA cells.
2. **Hard Budget & `.pkl` Requirement Exclusion (`> 60.000 GB`)**:
   - **Why `46.074 GB` `.pkl` is mandatory**: On IGVF, `IGVFFI4730IQLW.h5ad` (`87.425 GB` > 60 GB alone) is explicitly `derived_from: [IGVFFI6916RLYF, IGVFFI9448NBQC]`. `IGVFFI6916RLYF.h5` (`4.925 GB`) holds cell × gene RNA counts only; **`IGVFFI9448NBQC.pkl` (`46.074 GB`, `aggr_combined_df_full_clones_removed.pkl`) is the sole public guide-assignment file on IGVF**, making minimum H9 download = **`50.999 GB`**.
   - **Budget breach**: Active 4 lines (`CD4+ T CE0006864+CE0008162` `26.551 GB` + `DLD-1` `6.496 GB` + `i3Neurons` `0.289 GB` + `Calu-3` `0.744 GB` = **`34.080 GB`**) + `H9` (`50.999 GB`) = **`85.079 GB`**, breaching the **`60.000 GB`** cumulative limit by **`25.079 GB`** (and unpacking `46.074 GB` pandas pickle risks OOM on `36 GiB` host RAM).
   - **Per `P2` & `R10-1`**: Exclude `H9` prior to download (**document in manuscript as `"excluded pre-download due to high-MOI (~4) single-sgRNA scarcity (~5.5% single-guide, mean ~19 cells/gene) and 60 GB download budget (51.0 GB paired with 34.1 GB = 85.1 GB)"`**), leaving **`K_conf = 2` (`DLD-1`, `CD4+ T Rest`)**, both required to pass `H-C1 AND H-C2` (`2/2` confirmatory rule).

*(Also check: why don't we pull `Low_assays.h5` + `Low_se.rds` [`5.291 GB`] for DLD-1? Neither `Rscript` nor `rds2py` is installed; per-GSM HDF5 trios [`6.496 GB`] read natively in Python `h5py` with zero R dependency.)*

---

## 2. Deterministic Disjoint Split-Half Rules (`PREREG_C` Locked)

Within each half independently, per-batch NTC mean `mu_lib`/`mu_log` and SD `sd_lib` are computed **strictly from single-sgRNA NTC cells within that half** (falling back to global half `sd_lib` only if a sub-batch has `< 10` NTC cells or zero variance). `half_valid = (n_cells_half_a >= 1) & (n_cells_half_b >= 1)`.

| Line | Disjoint Split Axis | Exact Partition | Biological / Batch Variance Note |
|---|---|---|---|
| **DLD-1 (`GSE337988`)** | **10x Chromium channel (`SAM24482910_1..24`)** | **Half A (12 odd channels)**: `SAM24482910_{1,3,5,7,9,11,13,15,17,19,21,23}` (`GSM9864227,29,31,33,35,37,39,41,43,45,47,49`)<br>**Half B (12 even channels)**: `SAM24482910_{2,4,6,8,10,12,14,16,18,20,22,24}` (`GSM9864228,30,32,34,36,38,40,42,44,46,48,50`) | Interleaved low-MOI captures (`12 vs 12`); within-channel NTC batch normalization cleanly isolates capture noise. |
| **CD4+ T (`GSE314342` `Rest`)** | **Donor (`donor_id`, `P3`)** | **Half A**: Donor **`CE0006864` (`D4`)** (`24` lanes `R2L01..24` in run `CD4i_R2`; male Asian age 34; `2.69M` cells, `76,636` NTCs)<br>**Half B**: Donor **`CE0008162` (`D1`)** (`23` lanes `R1L01..23` in run `CD4i_R1`; female Hispanic/Latin American age 29; `3.07M` cells, `76,543` NTCs) | **Critical Disclosure (`P3`)**: Inter-donor biological variation (different sex/age/ancestry) **and** processing run variation (`CD4i_R2` vs `CD4i_R1`) are absorbed into the split-half noise ceiling `r_split` (each donor uses its own lane-level NTC baseline). |
| **WTC11 `i3Neurons` (`Tian 2021`)** | **10x lane (`lib1..4`)** | **Half A**: `{lib1, lib3}` (`GSM4632022, GSM4632025`)<br>**Half B**: `{lib2, lib4}` (`GSM4632023, GSM4632026`) *(fallback if scPerturb `.h5ad` omitted lane column: use the 4 GEO `.h5` lanes `0.316 GB` + `MOESM6_ESM.xlsx`)* | Interleaved 10x capture lanes (`2 vs 2`). |
| **Calu-3 (`Sunshine 2023`)** | **Batch/lane in `.obs`** | If `.obs` contains multiple batch/lane IDs, split sorted batch IDs by even/odd parity; if single merged batch in scPerturb `.h5ad`, split 50/50 stratified by `infection_state` (`bystander`/`infected`/`borderline`, `seed=42`) and record exact branch in `line_metadata.json`. | Ensures infection-state proportions are balanced between halves. |

*(H9 pre-download excluded under `P2`; if ever re-evaluated on a higher-RAM/expanded-budget run, split parity would be `LW558,560,...,580` [12 even] vs `LW559,561,...,581` [12 odd].)*

---

## 3. Target & Readout-Gene Mapping to HGNC (`/Users/mac/.cache/r9_p2_raw/hgnc_complete_set.txt`)

### 3.1 Strict Deterministic Resolution Precedence
Using the exact local HGNC release cached in Round 9 (`hgnc_complete_set.txt`, filtered to `status == "Approved"`), every target perturbation and expression readout `.var` gene is mapped as follows:
1. **Exact Approved Symbol**: `sym in approved_symbols` -> keep `sym`.
2. **Ensembl Gene ID (`ENSG...`)**: unambiguous `ensembl_gene_id -> approved_symbol` (resolves `95` legacy names in `CD4+ T` e.g. `ARMT1 -> DCPH1`, `C15orf48 -> COXFA4L3`).
3. **Unambiguous Previous HGNC Symbol (`prev_symbol`)**: resolve unique `prev_symbol -> approved_symbol` (resolves `176` symbols in `DLD-1` e.g. `ADAL -> MAPDA`, `AES -> TLE5`; `8` in `i3Neurons` e.g. `ATP5C1 -> ATP5F1C`, `GARS -> GARS1`, `YARS -> YARS1`; `13` in `Calu-3` e.g. `ADCK4 -> COQ8B`, `DDX58 -> RIGI`).
4. **Unambiguous Alias Symbol (`alias_symbol`)**: resolve only when unambiguous 1-to-1.
5. **Line-Specific Parsing & Filtering (`P4` Strictly Applied)**:
   - **DLD-1**: Inspect `.h5` `matrix/features` + `.crispr.umi_correct.h5` + `.hashing.umi_correct.h5`. Keep barcodes passing MULTI-seq singlet threshold with **strictly 1 assigned sgRNA**. Strip isoform/promoter suffixes `^(.+?)_(?:P\d|ENST|\d)` (`BRD4_ENST00000263377(W)_1 -> BRD4`; `NTC_*` -> NTC). Out of `4,954` target prefixes: `4,715` Approved + `176` `prev_symbol` (`4,891` canonical HGNC); `63` unmapped clone/lncRNA prefixes (`AC002094.3`, etc.) retained only if unambiguous HGNC or dropped with loss count audited in `line_metadata.json`.
   - **CD4+ T**: From CellRanger Flex `.h5`, join assigned guide call against `sgrna_library_metadata.suppl_table.csv` `(target_gene_name, target_gene_id)` (`12,680` Approved + `95` via `ENSG` = `12,775/12,778`; `3` non-HGNC dropped). Keep **strictly 1 sgRNA per cell** (matching single-sgRNA UMI threshold from `QC_summaries_per_sample_lane.csv` / CellRanger `protospacer_calls_per_cell`, zero `multi sgRNA`).
   - **i3Neurons**: Strictly 1 sgRNA (`CRISPRi` only; `176` Approved + `8` `prev_symbol` = `184/184`; `'non-targeting'` -> NTC).
   - **Calu-3**: Strictly single-gene library elements (`170` Approved + `13` `prev_symbol` = `183/183`; discard `22` dual-gene `_` elements; `non_targeting_pair_*` -> NTC).

### 3.2 Cross-Round HGNC Nomenclature Bridge (`R9 Frozen` ↔ `R10 New Lines`)
Inspection of Round 9 frozen files revealed **`164` legacy `prev_symbol` genes inside `round9/p2/shared/COMMON_GENES.txt`** (`6,148` genes, e.g., `AARS`, `ADAL`, `ARMT1`, `ATP5MD`, `C11orf54`, `GARS`, `YARS`) and `~75` legacy symbols in `replogle_*` / `nadig_*` `targets.txt`.
- **Zero Modification of R9 Files**: Frozen `round9/p2/shared/` artifacts remain untouched.
- **Loader-Level Canonical Projection (`round10/src/`)**: When matching source lines (`R9`) against target lines (`R10`), **both** source and target symbol namespaces are projected through the identical read-only `hgnc_complete_set.txt` canonicalizer (breaking rare 2-to-1 collisions by highest baseline control expression / first sorted index). This aligns all `164` legacy `COMMON_GENES` without silent loss and lifts baseline target overlaps cleanly (`DLD-1`: `668 -> 669`; `CD4+ T Rest DE`: `1,452 -> 1,490`; `i3Neurons`: `70 -> 72`; `Calu-3`: `63 -> 65`).

---

## 4. Standardized Output Schema & Count-Only QC Verification (After `R10_GO.txt`)

Each dataset outputs `round10/shared/lines/<line>/` containing:
- `genes.txt`, `targets.txt`, `ctrl_cells.npz` (`<= 20,000` batch-stratified single-sgRNA NTC cells, plus `half_a_ctrl_mean_*`, `half_b_ctrl_mean_*`), `line_metadata.json`,
- `pseudobulk_delta.npz` with **Addendum-B** arrays (`delta_{log1p,zscore,libnorm}`, `half_{a,b}_delta_*`, `half_valid`, `reproducible` [`n_cells >= 25` AND `ad_counts >= 50` AND `energy_test_p_value < 0.001`], `ad_counts`, `energy_test_p_value`) **plus `S7` covariate arrays (`n_cells`, `knockdown_depth_log1p`, `knockdown_fraction`, `half_a_n_cells`, `half_b_n_cells`, `target_raw_symbols`, `target_ensembl_ids`, `target_hgnc_status`)**.
- **Blind Count-Only Gate (`P1` & `R10-1`)**: Once `R10_GO.txt` is issued and ingestion completes, A-T1 reports **only** `raw_cells`, `single_sgRNA_qc_cells`, `n_ctrl_total`, `G_final`, `targets_dedup`, `targets_reproducible`, `targets_rep_shared_canon_source_rep` (intersecting `reproducible` targets in each R10 line with targets reproducible in `>= 1` R9 source line + expressed in `G_eval(t)` per `P5`), confirming which candidate confirmatory lines clear **`>= 50`** (`K_conf`). **Zero downstream evaluation or FCT correlation metrics are run by A-T1.**
