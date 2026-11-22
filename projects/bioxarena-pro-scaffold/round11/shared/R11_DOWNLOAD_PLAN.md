# Round 11 Download & Pre-Unblinding Ingestion Plan (`R11_DOWNLOAD_PLAN.md`)

- **Prepared by:** `A-T1` (Team A · Data Engineering)
- **Date:** `2026-11-22`
- **Hard Execution Lock:** **NO expression matrices have been downloaded or opened.** All figures below come exclusively from public API JSON/TSV metadata (`api.data.igvf.org`, `api.figshare.com`, `www.ncbi.nlm.nih.gov/geo`, `ftp.ncbi.nlm.nih.gov`), format specification PDFs (`igvf_HDF5_file_format.pdf`), published analysis scripts (`Feng 2026 scripts.zip`), small guide/cell metadata tables (`IGVFFI9754AGFB`, `IGVFFI4575UMXX`, `IGVFFI9546UECF`, `GSE210681_Dialout_guides_per_CBC.txt.gz`, `TargetedScreen_Cell-Metadata.tsv.gz`), and frozen Round 10 CD4i CSVs. **Zero matrix downloads may start before `/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round11/status/R11_GO.txt` is issued by Master Controller (`Claude`).**

---

## 1. System Constraints, Budget Lock, & Pre-Download Compliance (`Q3`, `Q6`)

- **Host Compute & Storage:** macOS (`sysctl hw.memsize = 38,654,705,664 B` = **`36.0 GiB` RAM**), `h5py 3.12.1` (`HDF5 1.14.5`), local disk free `~2.5 TiB` (`> 30x` budget). `round11/.gitignore` already ignores `shared/raw/`, `shared/lines/`, `*.h5`, `*.h5ad`, `*.npz`, `*.pkl`, `*.rds`, `*.tar`, `*.mtx.gz`.
- **Master Controller Download Budget Cap:** **`80.000 GB`** (`80,000,000,000 B` decimal).
- **Round 9 / Round 10 Overlap Prohibition (`Q6`):**
  - Strictly excludes all `6` R9 source lines and all `4` R10 lines (`dld1`, `cd4t_rest`, `wtc11_i3neuron`, `calu3`).
  - `cd4t_stim48` uses donors **`CE0008678` (`D3`)** and **`CE0010866` (`D2`)** (`Stim48hr` condition, `CD4i_R2`), completely donor-disjoint from R10 `cd4t_rest` (`CE0008162` `D1` & `CE0006864` `D4`, `Rest`), though from the same `GSE314342` Marson lab CD4+ T cohort (weaker lineage independence explicitly disclosed).
  - `feng_ipsc` (`TargetedScreen`) strictly excludes **`kolf_2`** (`64,009` cells) and **`kolf_3`** (`69,622` cells) because R9 source `nourreddine_kolf21j` derives from the KOLF2 lineage.
  - **Prohibited Cohorts:** `X-Atlas/Orion` and `Cas9` double-strand-break (`CRISPR-KO`) excluded.

### Cumulative Storage Budget Table (`All <= 80.000 GB`)

| Tier | Dataset Key | Required Files (Exact Count) | Exact Bytes (Core Required) | Decimal GB (Core) | Binary GiB (Core) | Optional / Already-Cached Reference | Max Bytes (Incl. Optional Verification) |
|---|---|---|---:|---:|---:|---|---:|
| Confirmatory 1 | `telohaec` | `IGVFFI8645YMGX.h5ad` + `IGVFFI9754AGFB.tsv.gz` + `GSE210681_Dialout_guides_per_CBC.txt.gz` (3 files) | `1,486,778,253` | `1.487` | `1.385` | — | `1,486,778,253` |
| Confirmatory 2 | `hcasmc` | `IGVFFI5999NXAF.h5` + `IGVFFI4575UMXX.tsv.gz` + `IGVFFI9546UECF.tsv.gz` (3 files) | `42,885,287,884` | `42.885` | `39.940` | Excluded: `IGVFFI9007DBZO.h5ad` (`87.29 GB`), `IGVFFI7271ZXXZ.h5` (`3.68 GB` pilot) | `42,885,287,884` |
| Confirmatory 3 | `cd4t_stim48` | `48` GEO `.h5` files (`GSM9393875..98`, `D2_Stim48hr` & `D3_Stim48hr`, `CD4i_R2L25..48`) | `27,387,332,602` | `27.387` | `25.506` | `3` shared CSVs in `round10/shared/raw/cd4t_rest/` (`9,998,730 B`) | `27,397,331,332` |
| Secondary (`H-D5`) | `feng_ipsc` | `TargetedScreen_RNA-UMI-Counts.csv.gz` + `TargetedScreen_Cell-Metadata.tsv.gz` (2 files) | `4,494,461,026` | `4.494` | `4.186` | `TargetedScreen_Guide-UMI-Counts.csv.gz` (`73,778,450 B`) | `4,568,239,476` |
| **TOTAL** | **All 4** | **56 core files (+ optional 4)** | **`76,253,859,765`** | **`76.254`** | **`71.017`** | **Remaining vs `80.000 GB`: `3.662–3.746 GB`** | **`76,337,636,945` (`76.338 GB`)** |

---

## 2. Exact Master Control, `P4` Single-sgRNA, & Split-Half Summary (Pre-Locked for `PREREG_D`)

Per Master Controller (`PREREG_D_AUDIT_R1.md` `Q2`), **only true `non-targeting` (NTC) cells** are used as primary control cells (consistent with R9 & R10); **`safe-targeting` guides are excluded from both primary NTC controls and perturbed targets**. Also, **"unassigned"** cells are strictly excluded even when original authors pooled them.

| Dataset Key | Primary NTC Definition & Guide Count | `safe-targeting` Distinction & Action | Single-sgRNA (`P4`) Filtering Rule | Split-Half Rule (`PREREG_D_AUDIT_R1 Q2`) | Per-Unit NTC Adequacy (`>= 15` Check) |
|---|---|---|---|---|---|
| `telohaec` | **`type == "non-targeting"`** (`400` guides, `guide_id` prefix `negative_control__`, `Dialout` `gene_symbol == "negative-control"`); `2,656` pre-QC dialout cells | **`type == "safe-targeting"`** (`600` guides, prefix `safe_targeting__`, `Dialout` `gene_symbol == "safe-targeting"`; **NOTE:** `538/600` have `targeting == True`, `62` `False`); `3,912` pre-QC cells -> **EXCLUDE** | **`max_umi_ct >= 4` & `X2nd_place_umi_ct < 4`** (cross-verified between `.h5ad` CRISPR layer and `Dialout` `671,432` CBCs) + single-gene plain symbols resolvable via frozen `HGNC` (exclude `199` `-TSS*` alt-promoters & `168` `Enhancer-*`/`9p21-*`/`*-and-*`) | **Sorted 10x channel/lane ID parity** (`1..20` matching `11AMDox_1..6`, `2PMDox1_1..7`, `2PMDox2_1..7`): `Half A` (`10` odd lanes: `1,364` NTC / `130,519` tgt pre-QC), `Half B` (`10` even lanes: `1,292` NTC / `120,899` tgt pre-QC) | **20/20 lanes pass** (`min = 31` NTC in Lane 12; `78..179` across other 19 lanes) |
| `hcasmc` | **`type == "non-targeting"`** (`125` guides, `guide_id` prefix `non-targeting_`, `targeting == False`) | **`type == "safe-targeting"`** (`125` guides, `guide_id` prefix `0Safe_`, **ALSO `targeting == False`** — never filter by `targeting == False` alone!) -> **EXCLUDE** | Extract guide rows (`feature_type != "Gene Expression"`) per QC-passing cell; apply Poisson-Gaussian floor + single-sgRNA dominance (`n_pos == 1` & `t1 >= 5.0` & `(t2 < 3.0 \| t1 / max(t2, 0.5) >= 3.0)`, matching `prepare_r10_lines.py`) + RNA QC (`MT <= 5–10%`, `ngenes >= 500` pre-locked in `build_hcasmc_parts`) | **Sorted `14` Parse sublibrary ID parity** (`3A..3H, 4A..4F` from `/matrix/barcodes` suffix `__<sublib>`): `Half A` (`3A,3C,3E,3G,4A,4C,4E` = 7), `Half B` (`3B,3D,3F,3H,4B,4D,4F` = 7). **Locked fallback:** if `/matrix/barcodes` has `< 2` parseable sublibrary suffixes, use `seed = 42` stratified 50/50 split | Audited post-download on guide counts only before `CONFIRMATORY_LINES_R11.txt`; fallback to half pooled NTC if any sublibrary has `< 15` NTC |
| `cd4t_stim48` | **`sgRNA.startswith("NTC-")`** (`992` guides, blank `target_gene_name`, `flag=False`) | **`0` safe-targeting guides** (`315` flagged mismatch/distal guides with blank `target_gene_name` excluded cleanly by `target_gene_name` non-empty check) | **Identical to frozen R10 `build_cd4t_rest_parts`:** `MT <= 5.0%` & `n_genes >= 500`; per-donor pooled Poisson-Gaussian guide threshold (`min_thr=5, max_thr=20, default_thr=7`); `n_pos == 1` & `t1 >= 5.0` & `(t2 < 3.0 \| t1 / max(t2, 0.5) >= 3.0)`; `CD4_EXCEL_GUIDE_REMAP = {"1-Jun": "JUN-1", "2-Jun": "JUN-2"}` | **Disjoint Donor Split (lexicographic per `Q2`):** `Half A = CE0008678` (`D3_Stim48hr`, `24` lanes `R2L25..48`, `77,334` NTC / `1,808,899` tgt), `Half B = CE0010866` (`D2_Stim48hr`, `24` lanes `R2L25..48`, `87,747` NTC / `1,965,311` tgt) | **48/48 donor-lanes pass** (`min = 3,027` NTC/lane for `D3`; `min = 3,381` NTC/lane for `D2`) |
| `feng_ipsc` | **`Guide_Call.startswith("NonTarget_")`** (`20` guides; **`7,266`** non-KOLF cells) | **`0` safe-targeting guides**; **`470,561` `"unassigned"` non-KOLF cells MUST BE STRICTLY EXCLUDED** (authors pooled `"unassigned"` + `"NonTarget"`) | Already called by `4a_add_guides.R` (`freq_max_guide > th_freq` & `max_UMI > 2` [`>= 3` UMIs]; `0` multi-guide delimiters in `Guide_Call`). Exclude `kolf_2`/`kolf_3` (`133,631`) & `"unassigned"` (`470,561`); keep **`7,266` NTC + `550,406` single-sgRNA targeting** across `17` non-KOLF lines (`444` target genes) | **Sorted `17` non-KOLF `Cell_Line` ID parity (`1`-based odd vs even per `Q2`):** `Half A` (`9` odd lines, `9` donors: `eipl_1, fiaj_1, iudw_1, jejf_2, oikd_2, paab_3, pipw_4, tolg_4, zapk_3` -> `4,091` NTC / `311,764` tgt), `Half B` (`8` even lines, `8` donors: `eipl_3, fiaj_3, iudw_4, jejf_3, oikd_5, paab_4, pipw_5, tolg_6` -> `3,175` NTC / `238,642` tgt) | **17/17 lines pass** (`min = 86` NTC in `fiaj_3`; `>= 86` per line). `142` `(Pool x Inlet)` units avg `~51` NTC -> center by `Cell_Line x Inlet` if `>= 15` NTC, else line/half pooled NTC (`PREREG_D §2.3`) |

> **Important `feng_ipsc` Correction vs Earlier Notes:** `Feng 2026 TargetedScreen` (`Pica`) has **`19` HipSci lines (`10` donors)** total, leaving **`17` non-KOLF lines (`9` donors)** after excluding `kolf_2` & `kolf_3` (`43` non-KOLF lines / `30` donors spanned all `3` screens combined, where `Fitness` has only `32` NTCs and `NonFitness` has only `10` NTCs globally and were already excluded). All `444` target genes in these `17` `TargetedScreen` lines have `>= 25` total cells AND `>= 10` cells in both `Half A` and `Half B` (`333` in `R9_COMMON`).

---

## 3. Dataset-by-Dataset File Manifests, Protocols, & Verification Checks

### 3.1 Confirmatory 1 — `telohaec` (`TeloHAEC_2285`, `IGVFDS7943IIWZ` / `GSE210681`)

- **Biological State & Independence:** Immortalized human aortic endothelial cells (TeloHAEC, resting DOX-inducible KRAB-dCas9, arterial endothelium — lineage unseen in R9 & R10). **Target dir:** `round11/shared/raw/telohaec/`.
- **Exact Files (`1,486,778,253` B = `1.487 GB`):**
  1. `IGVFFI8645YMGX.h5ad` (`filtered cell by gene and guide matrix`, `cellranger-v8.0.0`):
     - **URL:** `https://api.data.igvf.org/matrix-files/IGVFFI8645YMGX/@@download/IGVFFI8645YMGX.h5ad` (S3: `s3://igvf-public/2025/06/02/19c2e2af-7ec5-4d90-9829-f5f005c1f8ea/IGVFFI8645YMGX.h5ad`)
     - **Size / MD5:** `1,473,190,952` B / `167d157d6ffd318a415739113fdea649`
  2. `IGVFFI9754AGFB.tsv.gz` (`GRCh38 guide RNA sequences`, `37,637` guides):
     - **URL:** `https://api.data.igvf.org/tabular-files/IGVFFI9754AGFB/@@download/IGVFFI9754AGFB.tsv.gz`
     - **Size / MD5:** `729,687` B / `10de9535f3914a1a7bcfa94b0ce11e18`
  3. `GSE210681_Dialout_guides_per_CBC.txt.gz` (`671,432` cell barcodes, PCR dialout):
     - **URL:** `ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE210nnn/GSE210681/suppl/GSE210681_Dialout_guides_per_CBC.txt.gz` (HTTPS mirror supported)
     - **Size / MD5:** `12,857,614` B / `96a6511e4e4021124c994f332f232c61`
- **Channel / Lane Metadata (`20` 10x Channels `-1`..`-20`):**
  - Paired IGVF `MeasurementSet` / `AuxiliarySet` (`40` input file sets + `IGVFDS9204JREI` guide design) and GEO GSMs (`GSM6435403..GSM6435422` scRNA-seq + `GSM6435423..GSM6435442` Dialout): `11AMDox_1..6` (`6`), `2PMDox1_1..7` (`7`), `2PMDox2_1..7` (`7`).
  - `Dialout` P4 counts (`max_umi_ct >= 4` & `X2nd_place_umi_ct < 4`) per lane `1..20` (NTC / SAFE / TGT):
    `L1: 117/176/10,757`, `L2: 100/116/8,959`, `L3: 103/131/8,260`, `L4: 78/112/8,096`, `L5: 100/147/9,075`, `L6: 112/147/9,190`, `L7: 147/224/14,013`, `L8: 132/201/12,899`, `L9: 130/226/13,845`, `L10: 175/244/16,263`, `L11: 146/233/14,936`, `L12: 31/74/3,593`, `L13: 125/199/12,692`, `L14: 166/231/15,607`, `L15: 149/228/14,500`, `L16: 164/247/15,421`, `L17: 179/255/17,545`, `L18: 169/255/15,692`, `L19: 168/229/14,896`, `L20: 165/237/15,179`.
- **Target Resolution Caution (`-TSS2` vs Plain vs Enhancer):**
  - Out of `2,618` targeting `gene_symbol` values present in `Dialout P4`, `2,251` are plain gene symbols (`2,133` directly approved HGNC, `1,064` in `R9_COMMON`), `199` end in `-TSS1/-TSS2/-TSS3` (of which **`192` overlap plain gene names as distinct alternative promoters, e.g. `AMPD2` vs `AMPD2-TSS2` — MUST NOT strip `-TSS*` and merge**), and `168` are enhancers/multi-locus (`Enhancer-*`, `9p21-*`, `*-and-*`). Strictly retain HGNC-resolvable single-gene TSS elements only (`plain` symbols resolving cleanly via `HGNCLookup`, excluding `-TSS*`, `Enhancer-*`, `9p21-*`, `*-and-*`).

### 3.2 Confirmatory 2 — `hcasmc` (`HCASMC-hTERT`, `IGVFDS6248NDUW`)

- **Biological State & Independence:** Immortalized human coronary artery smooth muscle cells (`Sp-dCas9-ZIM3-KRAB`, `Parse Biosciences Split-Pool v1.2.1`, MOI `0.1–0.2`, smooth muscle lineage unseen in R9 & R10). **Target dir:** `round11/shared/raw/hcasmc/`.
- **Exact Files (`42,885,287,884` B = `42.885 GB`):**
  1. `IGVFFI5999NXAF.h5` (`.../large_analysis/filtered_feature_barcode_matrix.hdf5`):
     - **URL:** `https://api.data.igvf.org/matrix-files/IGVFFI5999NXAF/@@download/IGVFFI5999NXAF.h5` (S3: `s3://igvf-public/2025/09/16/cfe59b71-8720-4d6f-9ece-12ffd4f25fc1/IGVFFI5999NXAF.h5`)
     - **Size / MD5:** `42,885,130,097` B / `17257b6ee904e06404e4d18e52bd5850`
     - **Provenance:** `analysis_step_version /052ad56e-ca4d-4b6f-abd6-a6d9d5055716/` (`parse-biosciences-v1.2.1`, workflow `IGVFWF9728QHLL`); `derived_from` `/tabular-files/IGVFFI4575UMXX/`.
  2. `IGVFFI4575UMXX.tsv.gz` (`GRCh38 GENCODE 43 guide RNA sequences`, direct `derived_from` table of `IGVFFI5999NXAF.h5`):
     - **URL:** `https://api.data.igvf.org/tabular-files/IGVFFI4575UMXX/@@download/IGVFFI4575UMXX.tsv.gz`
     - **Size / MD5:** `80,945` B / `38d9d80dfa4957cf429e4d990399543e`
  3. `IGVFFI9546UECF.tsv.gz` (companion curated guide TSV in `IGVFDS9377AZJL`; identical `3,163` `guide_id`/`spacer` set, plus genomic coordinates on `safe-targeting`):
     - **URL:** `https://api.data.igvf.org/tabular-files/IGVFFI9546UECF/@@download/IGVFFI9546UECF.tsv.gz`
     - **Size / MD5:** `76,842` B / `1b444eb41eaec5162e7bf9d72c8ee171`
- **Parse Sublibrary Unit Metadata (`14` Sublibraries = `14` `MeasurementSet` / `IGVFSM*` Samples):**
  - `Half A` (`7` odd sorted sublibraries): `3A` (`IGVFDS2119CSBN` / `IGVFSM9407HJID`), `3C` (`IGVFDS5571KMQO` / `IGVFSM7066JLBS`), `3E` (`IGVFDS5825NUBC` / `IGVFSM7595THKR`), `3G` (`IGVFDS1025EOLK` / `IGVFSM0725XZPH`), `4A` (`IGVFDS6062YTII` / `IGVFSM7917PCZW`), `4C` (`IGVFDS5293PXDV` / `IGVFSM0766JREB`), `4E` (`IGVFDS9889ACMT` / `IGVFSM8554KQDX`).
  - `Half B` (`7` even sorted sublibraries): `3B` (`IGVFDS2137AEUK` / `IGVFSM2944FBQP`), `3D` (`IGVFDS1805GZHW` / `IGVFSM5360BYEX`), `3F` (`IGVFDS3435YNUD` / `IGVFSM5106OTOD`), `3H` (`IGVFDS8137QEAT` / `IGVFSM6710JFWC`), `4B` (`IGVFDS1372ZQKQ` / `IGVFSM3202EVMF`), `4D` (`IGVFDS1798VANO` / `IGVFSM1489CASX`), `4F` (`IGVFDS9901WEMX` / `IGVFSM2781BJZV`).
- **Chunked / Out-of-Core Streaming Feasibility Analysis for `IGVFFI5999NXAF.h5` (`42.885 GB` vs `36.0 GiB` Host RAM) — ZERO Expression Values Inspected:**
  - **Specification Verified:** IGVF links `IGVFFI5999NXAF.h5` (`submitted_file_name: .../filtered_feature_barcode_matrix.hdf5`) directly to document `/documents/f8806e2e-92e8-45e1-9178-72d7fdaaf70c/` (**10x/Cell Ranger HDF5 Feature-Barcode Matrix Format**). As verified in `igvf_HDF5_file_format.pdf`, `/matrix` stores a **Compressed Sparse Column (CSC)** matrix with `data` (`uint32`), `indices` (`uint32`), `indptr` (`uint32`/`uint64` of length `N_cells + 1`), `shape` `(n_features, n_cells)`, `barcodes` (UTF-8 strings), and `features/{_all_tag_keys, derivation, feature_type, genome, id, name}`.
  - **Why Naive Loading Fails (`OOM` Risk):** Even in float32 sparse CSR/CSC or Scanpy `read_10x_h5()` whole-matrix load, a `42.89 GB` HDF5 file can expand past our **`38,654,705,664 B` (`36.0 GiB`)** physical RAM when combined with gzip decompression / dense projections.
  - **Safe Streaming Path (Already Proven in `prepare_r10_lines.py`):**
    1. Open with `h5py.File(..., "r")` (`h5py 3.12.1`) and read ONLY structural vectors first: `/matrix/shape`, `/matrix/features/{id,name,feature_type}`, `/matrix/barcodes`, and `/matrix/indptr` (`~8–80 MB` even at multi-million columns).
    2. Stream RNA QC (`umi, ngenes, pct_mt`) and extract guide-capture rows (`feature_type != "Gene Expression"`) in fixed `batch_cols = 20,000` column windows via `indptr[cs:ce+1]` -> `data[p0:p1], indices[p0:p1]` (`compute_10x_csc_rna_qc`, `< 300 MB` transient slice RAM).
    3. Compute Pass-2 top-5000 HVG / Anderson-Darling / 20-PC Energy Distance sequentially via `make_10x_csc_subset_readers` & `prepare_r10_lines.py` out-of-core target-group sharding (`<= 4–6 GB` peak RSS).
  - **Pre-Unblinding Feasibility Exclusion Gate (`Q3`):** Immediately after `R11_GO.txt` download and **strictly before computing any prediction or correlation metrics**, verify `/matrix/{shape,features/feature_type,barcodes,indptr,indices,data}` chunked sliceability, presence of `non-targeting_` guides, and `1`–`2` parseable sublibrary suffixes. If `IGVFFI5999NXAF.h5` is unchunked/corrupted or lacks extractable single-sgRNA guide rows, log a pre-`CONFIRMATORY_LINES_R11.txt` feasibility exclusion in the QC report.

### 3.3 Confirmatory 3 — `cd4t_stim48` (`GSE314342` `Stim48hr` Donors `CE0008678` & `CE0010866`)

- **Biological State & Independence Notice (`Q2`, `Q6`):** Primary human CD4+ alpha-beta T cells stimulated for **48 hr** (`CD4i_R2`, donors `D3 = CE0008678` [23F African American] & `D2 = CE0010866` [22F African American]). Donors and stimulation state (resting vs 48h TCR-activated) are completely new relative to R10 `cd4t_rest` (`CE0008162` D1 & `CE0006864` D4), **but cell type (primary CD4+ T) and originating lab (Marson lab) match R10 `cd4t_rest`** -> weaker biological independence explicitly flagged. **Target dir:** `round11/shared/raw/cd4t_stim48/`.
- **Exact Files (`48` GEO `.h5` files + `3` shared CSVs):**
  - `24` files `GSM9393875..GSM9393898` (`CD4i_R2L25..CD4i_R2L48`) for **`D3_Stim48hr` (`CE0008678`, `Half A`)**: **`13,590,537,285` B** (`13.591 GB`), `2,607,532` total cells, `77,334` single-NTC (`3,027–3,484` NTC/lane), `1,808,899` single-sgRNA targeting.
  - `24` files `GSM9393875..GSM9393898` (`CD4i_R2L25..CD4i_R2L48`) for **`D2_Stim48hr` (`CE0010866`, `Half B`)**: **`13,796,795,317` B** (`13.797 GB`), `2,863,571` total cells, `87,747` single-NTC (`3,381–3,879` NTC/lane), `1,965,311` single-sgRNA targeting.
  - Combined `48` `.h5` files = **`27,387,332,602` B** (`27.387 GB`; see Appendix A for all `48` exact GSM IDs, filenames, FTP/HTTPS URLs, and verified FTP byte sizes). Post-download integrity verified via exact FTP byte match + HDF5 `h5py` structural open across all `48` `.h5` files (GEO GSM text summaries omit MD5 hashes).
  - Shared metadata CSVs (`9,998,730` B total; **already cached locally** in `round10/shared/raw/cd4t_rest/` and reusable read-only/symlinked without duplicate download):
    - `sgrna_library_metadata.suppl_table.csv` (`9,935,249` B, `md5 = 6736271405e9ce94be0b5f8e8412a2c5`)
    - `QC_summaries_per_sample_lane.csv` (`60,515` B, `md5 = 0c5b3a3fbfefb9f8e736a1678228ad6d`)
    - `sample_metadata.suppl_table.csv` (`2,966` B, `md5 = 29ca85689789eaff35d255bc0f0db0ad`)

### 3.4 Secondary (`SA-KOLF` / `H-D5`, Non-Confirmatory) — `feng_ipsc` (`Feng 2026 TargetedScreen`)

- **Biological Role (`PREREG_D_AUDIT_R1 Q2`):** Undifferentiated HipSci human iPSCs (`TargetedScreen` / `Pica`, MIT License, Figshare `27989294 v2`, `DOI: 10.6084/m9.figshare.27989294.v2`). **Never enters confirmatory `K_conf`** because its cell type matches R9 source `nourreddine_kolf21j`; retained for `H-D5` (`all_nontarget` comparison of FCT vs KOLF2.1J mismatched baseline). **Target dir:** `round11/shared/raw/feng_ipsc/`.
- **Exact Files (`4,494,461,026` B core; `4,568,239,476` B incl. optional `Guide-UMI-Counts`):**
  1. `TargetedScreen_RNA-UMI-Counts.csv.gz` (`id=51072902`): `https://ndownloader.figshare.com/files/51072902` — **`4,480,969,391` B** / `md5 = f0d613b1ea8147425501baddba034158`
  2. `TargetedScreen_Cell-Metadata.tsv.gz` (`id=51051809`): `https://ndownloader.figshare.com/files/51051809` — **`13,491,635` B** / `md5 = 980a440a8132a20567c4c3610e6a5b20`
  3. *(Optional guide-UMI audit table)* `TargetedScreen_Guide-UMI-Counts.csv.gz` (`id=51051815`): `https://ndownloader.figshare.com/files/51051815` — **`73,778,450` B** / `md5 = 96342f750a3cd983d33bfea64aa2559f`
- **Line, Donor, & Inlet Composition (Observed in `TargetedScreen_Cell-Metadata.tsv.gz`):**
  - `19` lines total (`1,161,864` cells); exclude `kolf_2` (`64,009`) & `kolf_3` (`69,622`) -> **`17` non-KOLF lines across `9` donors** (`1,028,233` total = `470,561` `"unassigned"` excluded + `7,266` `NonTarget_*` NTC + `550,406` single-sgRNA targeting across `444` genes):
    - **`Half A` (`9` odd lines, `9` donors, `4,091` NTC / `311,764` tgt):** `eipl_1` (`496`/`38,203`), `fiaj_1` (`549`/`41,249`), `iudw_1` (`575`/`43,894`), `jejf_2` (`605`/`46,953`), `oikd_2` (`177`/`13,950`), `paab_3` (`397`/`28,521`), `pipw_4` (`488`/`36,690`), `tolg_4` (`182`/`16,550`), `zapk_3` (`622`/`45,754`).
    - **`Half B` (`8` even lines, `8` donors, `3,175` NTC / `238,642` tgt):** `eipl_3` (`501`/`36,679`), `fiaj_3` (`86`/`8,168`), `iudw_4` (`433`/`34,269`), `jejf_3` (`484`/`34,367`), `oikd_5` (`668`/`44,839`), `paab_4` (`543`/`41,964`), `pipw_5` (`92`/`9,233`), `tolg_6` (`368`/`29,123`).
  - `Batch` units span `7` `Pool`s (`Pool1..Pool7`, all `Day3`) across `142` 10x `Inlet`s (`1..144` minus missing `[63, 68]`).
- **Mandatory Streaming Ingest Constraint (`36 GiB` Host RAM):**
  - `TargetedScreen_RNA-UMI-Counts.csv.gz` is transposed (`rows = genes` [`ENSG:symbol:Gene-Expression`], `columns = 1,161,864 Cell_IDs`). Never invoke `pd.read_csv()` on the full table (naive dense load `> 300 GB`). Instead, stream `gzip.open(..., "rt")` row-by-row after `R11_GO.txt`, parse the header `Cell_ID` list to retain **only** the `~557,672` kept column indices (`non-KOLF & Guide_Call != "unassigned"`), discard `:CRISPR-Guide-Capture` rows if present, and accumulate sparse CSR/CSC integer blocks into an on-disk `.h5ad`/`.h5` cache (`< 3 GB` peak RSS).

---

## Appendix A: Exact `48` GEO Manifest for `cd4t_stim48` (`GSE314342` `CD4i_R2L25..R2L48`)

All URLs below are `ftp://ftp.ncbi.nlm.nih.gov/...` (or identical `https://ftp.ncbi.nlm.nih.gov/...`); sizes queried directly via `ftplib.FTP.size()` (`27,387,332,602` B exact sum: `D2 = 13,796,795,317` B, `D3 = 13,590,537,285` B).

| GSM | Donor | Filename | Exact Bytes | FTP URL |
|---|---|---|---:|---|
| `GSM9393875` | `D2_CE0010866` | `GSM9393875_CD4i_R2L25_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `541,464,459` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393875/suppl/GSM9393875_CD4i_R2L25_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393875` | `D3_CE0008678` | `GSM9393875_CD4i_R2L25_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `535,766,249` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393875/suppl/GSM9393875_CD4i_R2L25_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393876` | `D2_CE0010866` | `GSM9393876_CD4i_R2L26_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `559,751,620` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393876/suppl/GSM9393876_CD4i_R2L26_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393876` | `D3_CE0008678` | `GSM9393876_CD4i_R2L26_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `554,367,686` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393876/suppl/GSM9393876_CD4i_R2L26_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393877` | `D2_CE0010866` | `GSM9393877_CD4i_R2L27_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `570,825,846` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393877/suppl/GSM9393877_CD4i_R2L27_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393877` | `D3_CE0008678` | `GSM9393877_CD4i_R2L27_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `556,904,595` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393877/suppl/GSM9393877_CD4i_R2L27_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393878` | `D2_CE0010866` | `GSM9393878_CD4i_R2L28_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `571,245,020` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393878/suppl/GSM9393878_CD4i_R2L28_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393878` | `D3_CE0008678` | `GSM9393878_CD4i_R2L28_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `561,314,616` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393878/suppl/GSM9393878_CD4i_R2L28_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393879` | `D2_CE0010866` | `GSM9393879_CD4i_R2L29_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `625,758,561` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393879/suppl/GSM9393879_CD4i_R2L29_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393879` | `D3_CE0008678` | `GSM9393879_CD4i_R2L29_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `615,267,418` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393879/suppl/GSM9393879_CD4i_R2L29_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393880` | `D2_CE0010866` | `GSM9393880_CD4i_R2L30_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `568,375,214` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393880/suppl/GSM9393880_CD4i_R2L30_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393880` | `D3_CE0008678` | `GSM9393880_CD4i_R2L30_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `557,985,492` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393880/suppl/GSM9393880_CD4i_R2L30_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393881` | `D2_CE0010866` | `GSM9393881_CD4i_R2L31_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `593,264,701` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393881/suppl/GSM9393881_CD4i_R2L31_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393881` | `D3_CE0008678` | `GSM9393881_CD4i_R2L31_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `582,515,837` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393881/suppl/GSM9393881_CD4i_R2L31_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393882` | `D2_CE0010866` | `GSM9393882_CD4i_R2L32_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `588,018,510` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393882/suppl/GSM9393882_CD4i_R2L32_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393882` | `D3_CE0008678` | `GSM9393882_CD4i_R2L32_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `580,799,811` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393882/suppl/GSM9393882_CD4i_R2L32_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393883` | `D2_CE0010866` | `GSM9393883_CD4i_R2L33_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `622,523,435` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393883/suppl/GSM9393883_CD4i_R2L33_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393883` | `D3_CE0008678` | `GSM9393883_CD4i_R2L33_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `614,199,302` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393883/suppl/GSM9393883_CD4i_R2L33_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393884` | `D2_CE0010866` | `GSM9393884_CD4i_R2L34_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `587,131,134` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393884/suppl/GSM9393884_CD4i_R2L34_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393884` | `D3_CE0008678` | `GSM9393884_CD4i_R2L34_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `579,412,433` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393884/suppl/GSM9393884_CD4i_R2L34_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393885` | `D2_CE0010866` | `GSM9393885_CD4i_R2L35_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `587,878,471` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393885/suppl/GSM9393885_CD4i_R2L35_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393885` | `D3_CE0008678` | `GSM9393885_CD4i_R2L35_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `579,504,815` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393885/suppl/GSM9393885_CD4i_R2L35_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393886` | `D2_CE0010866` | `GSM9393886_CD4i_R2L36_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `582,783,193` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393886/suppl/GSM9393886_CD4i_R2L36_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393886` | `D3_CE0008678` | `GSM9393886_CD4i_R2L36_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `575,748,051` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393886/suppl/GSM9393886_CD4i_R2L36_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393887` | `D2_CE0010866` | `GSM9393887_CD4i_R2L37_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `603,286,635` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393887/suppl/GSM9393887_CD4i_R2L37_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393887` | `D3_CE0008678` | `GSM9393887_CD4i_R2L37_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `590,757,077` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393887/suppl/GSM9393887_CD4i_R2L37_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393888` | `D2_CE0010866` | `GSM9393888_CD4i_R2L38_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `583,119,601` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393888/suppl/GSM9393888_CD4i_R2L38_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393888` | `D3_CE0008678` | `GSM9393888_CD4i_R2L38_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `569,284,608` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393888/suppl/GSM9393888_CD4i_R2L38_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393889` | `D2_CE0010866` | `GSM9393889_CD4i_R2L39_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `556,291,146` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393889/suppl/GSM9393889_CD4i_R2L39_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393889` | `D3_CE0008678` | `GSM9393889_CD4i_R2L39_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `550,032,451` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393889/suppl/GSM9393889_CD4i_R2L39_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393890` | `D2_CE0010866` | `GSM9393890_CD4i_R2L40_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `572,369,966` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393890/suppl/GSM9393890_CD4i_R2L40_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393890` | `D3_CE0008678` | `GSM9393890_CD4i_R2L40_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `564,154,107` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393890/suppl/GSM9393890_CD4i_R2L40_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393891` | `D2_CE0010866` | `GSM9393891_CD4i_R2L41_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `595,658,476` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393891/suppl/GSM9393891_CD4i_R2L41_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393891` | `D3_CE0008678` | `GSM9393891_CD4i_R2L41_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `583,579,744` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393891/suppl/GSM9393891_CD4i_R2L41_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393892` | `D2_CE0010866` | `GSM9393892_CD4i_R2L42_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `603,037,243` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393892/suppl/GSM9393892_CD4i_R2L42_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393892` | `D3_CE0008678` | `GSM9393892_CD4i_R2L42_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `591,409,318` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393892/suppl/GSM9393892_CD4i_R2L42_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393893` | `D2_CE0010866` | `GSM9393893_CD4i_R2L43_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `548,210,994` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393893/suppl/GSM9393893_CD4i_R2L43_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393893` | `D3_CE0008678` | `GSM9393893_CD4i_R2L43_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `541,956,987` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393893/suppl/GSM9393893_CD4i_R2L43_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393894` | `D2_CE0010866` | `GSM9393894_CD4i_R2L44_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `581,247,205` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393894/suppl/GSM9393894_CD4i_R2L44_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393894` | `D3_CE0008678` | `GSM9393894_CD4i_R2L44_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `568,665,092` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393894/suppl/GSM9393894_CD4i_R2L44_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393895` | `D2_CE0010866` | `GSM9393895_CD4i_R2L45_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `523,469,509` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393895/suppl/GSM9393895_CD4i_R2L45_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393895` | `D3_CE0008678` | `GSM9393895_CD4i_R2L45_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `522,331,797` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393895/suppl/GSM9393895_CD4i_R2L45_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393896` | `D2_CE0010866` | `GSM9393896_CD4i_R2L46_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `571,845,946` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393896/suppl/GSM9393896_CD4i_R2L46_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393896` | `D3_CE0008678` | `GSM9393896_CD4i_R2L46_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `562,083,137` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393896/suppl/GSM9393896_CD4i_R2L46_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393897` | `D2_CE0010866` | `GSM9393897_CD4i_R2L47_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `527,199,793` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393897/suppl/GSM9393897_CD4i_R2L47_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393897` | `D3_CE0008678` | `GSM9393897_CD4i_R2L47_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `523,649,181` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393897/suppl/GSM9393897_CD4i_R2L47_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393898` | `D2_CE0010866` | `GSM9393898_CD4i_R2L48_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `532,038,639` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393898/suppl/GSM9393898_CD4i_R2L48_D2_Stim48hr_sample_filtered_feature_bc_matrix.h5` |
| `GSM9393898` | `D3_CE0008678` | `GSM9393898_CD4i_R2L48_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` | `528,847,481` | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9393nnn/GSM9393898/suppl/GSM9393898_CD4i_R2L48_D3_Stim48hr_sample_filtered_feature_bc_matrix.h5` |

R11_PLAN_DONE
