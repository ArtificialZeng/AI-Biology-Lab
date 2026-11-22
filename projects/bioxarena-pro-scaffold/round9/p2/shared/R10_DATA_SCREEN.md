# Round 10 (`R10-0`) New Dataset Screening Report

**Agent**: A-T1 (Team A · Data Engineering)  
**Date**: 2026-11-22  
**Audit Method**: API/FTP metadata, target/sgRNA tables, and peer-reviewed / preprint supplementary tables ONLY (**zero** expression matrices downloaded; **zero** evaluations run; **zero** X-Atlas/Orion data).

---

## 1. Round 9 Baseline & Round 10 Screening Gates

- **Existing 6 Source Lines (`round9/p2/shared/lines/`)**:
  `nadig_hepg2` (`2,343` targets), `nadig_jurkat` (`2,350`), `nourreddine_kolf21j` (`95`), `replogle_k562_essential` (`2,049`), `replogle_rpe1` (`2,351`), `vcc2025_h1` (`300`).
  - **`union6` Unique Target Count**: **`2,655`** unique HGNC symbols.
  - **Round 9 Cumulative Download Footprint**: `56,597,332,700 B` (`56.597333 GB`).
  - **Round 10 Additional Download Budget**: `<= 60.000 GB`.
- **Mandatory Acceptance Gates**:
  1. Never-used human cell line (distinct from `HepG2`, `Jurkat`, `KOLF2.1J`, `K562`, `RPE1`, `H1`).
  2. CRISPRi (or clearly comparable knockdown) single-cell perturbation sequencing (`Perturb-seq` / `CROP-seq`).
  3. Non-targeting control (`NTC`) cells present.
  4. **`>= 50` shared targets** with existing `union6` (`2,655` targets).
  5. Open public license (`CC BY`, `MIT`, or open public GEO / IGVF).
  6. Download footprint **`<= 60 GB`** (within Round 10 ceiling).
  7. **Strictly NO X-Atlas / Orion** or derivatives.

---

## 2. Mandatory Infrastructure Checks (Arc VCC 2026 & Nourreddine 2026)

### 2.1 Arc VCC 2026 Public Release Check
- **GCS API Check**: `https://storage.googleapis.com/storage/v1/b/arc-institute-virtual-cell-atlas/o?prefix=virtual-cell-challenge/&delimiter=/`
  - Returned prefixes: `['virtual-cell-challenge/2025/']` only.
  - `https://storage.googleapis.com/storage/v1/b/arc-institute-virtual-cell-atlas/o?prefix=virtual-cell-challenge/2026/` returned `items: []` (**empty**) as of 2026-11-22.
- **Public Competition Release Status**: Arc Virtual Cell Challenge 2026 (`Contexts A, B, C` validation on 2026-08-20; `D, E, F` test on 2026-10-22) is structured as a zero-shot generalization benchmark providing only unperturbed NTC controls (~18,400 cells/context) and 300 target gene symbols; **no labeled Perturb-seq perturbation response training/validation expression matrices for 2026 contexts exist in the public GCS bucket**.

### 2.2 Nourreddine et al. (2026, *Nat. Biotechnol.*) KOLF2.1J Pseudobulk Check
- **Publication / DOI**: *Nature Biotechnology* (2026), `10.1038/s41587-026-03199-w`, `PMID:42386990`.
- **Figshare+ Deposit (`https://api.figshare.com/v2/articles/27261219`, License: `CC BY 4.0`)**:
  - Contains **only** 6 single-cell `.h5ad` files and **zero standalone compact pseudobulk CSV/summary files**:
    - `KOLF_Pan_Genome_QC_Filtered.h5ad`: `189,393,177,972 B` (`189.3932 GB`, exceeds `60 GB`)
    - `KOLF_Strong_Perturbations.h5ad`: `46,718,752,086 B` (`46.7188 GB`, fits `<= 60 GB`, **same `KOLF2.1J` cell line** — useful only if expanding Round 9 `KOLF2.1J`, not a new line)
    - `KOLF_Chromatin_Modifiers_QC_Filtered.h5ad`: `4.6356 GB`; `Strong`: `1.5289 GB` (already ingested in R9)
    - `KOLF_Metabolic_Enzymes_QC_Filtered.h5ad`: `6.9844 GB`; `Strong`: `1.3263 GB` (already ingested in R9)
- **GitHub Repository (`https://github.com/y-doctor/KOLF2.1J_Perturbation_Cell_Atlas`)**:
  - Inspected `main` tree (`.gitignore`, `LICENSE`, `README.md`, `docs/`, `environment.yml`, `psp/`); `input_files/` referenced in `README.md` is not uploaded to `main`, and `psp/` provides only analysis/DESeq2 pseudobulk code, **no precomputed pseudobulk table**.

---

## 3. Tier A — Qualifying New Human Cell Lines (PASS All Gates)

### Summary Matrix (`>= 50` Shared Targets, CRISPRi scRNA-seq, Open License, `<= 60 GB` Option Available)

| Rank | Cell Line (Tissue / Type) | Modality | Dataset / Citation | License | Download Size (`<= 60 GB` Option) | Cells / NTC Cells | Unique Targets | `shared_union6` | Batch Info (Split-Half) |
|---|---|---|---|---|---|---|---|---|---|
| **1** | **DLD-1** (Colorectal adenocarcinoma, epithelial) | `Zim3-dCas9` CRISPRi | **GEO `GSE337988`** (Yeung et al., *bioRxiv* Jul 2026) | Open GEO (`10.64898/2026.07.10.737863` `CC BY`) | **`4.621 GB`** (`Low_assays.h5`; +`0.670 GB` `.rds`; `7.672 GB` Low+Med+High `.h5`) | **`1,338,060`** est. GEX (`Low`), `2,152,003` (`Low+Med+High`); **`504` NTC guides** | **`4,954`** | **`668`** | **Yes** (`24` Low-MOI 10x channels + MULTI-seq) |
| **2** | **Primary CD4+ αβ T Cells** (Healthy donors, immune) | `ZIM3KRAB-dCas9` CRISPRi | **GEO `GSE314342`** + S3 `marson2025_data` (Zhu et al., *Cell* Oct 2026) | **`MIT License`** (S3/VCM) + Open GEO (`PMC13614707`) | **`54.582 GB`** (`Rest` 94 `.h5` files; `2-donor` ~`27 GB`; or **`44.567 GB`** `pseudobulk_merged.h5ad`) | `Rest`: **`11,487,117`** total, **`313,989` NTC single-sgRNA**, `7,072,306` tgt single-sgRNA | **`12,778`** lib (`11,287` `Rest` DE-tested) | **`2,482`** lib / **`1,452`** (`Rest` DE) | **Yes** (`4` donors × `2` 10x Flex runs × `94` lanes) |
| **3** | **WTC11 i3Neurons** (hiPSC-derived excitatory neurons) | `Sp-dCas9-KRAB` CRISPRi (CROP-seq) | **Zenodo `13350497`** (`TianKampmann2021_CRISPRi.h5ad`) / `GSE152988` (Tian et al., *Nat Neurosci* 2021) | **`CC BY 4.0`** | **`289.39 MB`** (`0.289 GB` `.h5ad`; or `315.54 MB` 4× `.h5`) | **`27,189`** assigned cells (`368` NTC cells across `6` NTC sgRNAs) | **`184`** | **`70`** | **Yes** (`4` 10x Chromium lanes `lib1..lib4`) |
| **4** | **H9 Cortical Neurons** (WA09 female hESC Day-6 cortical neurons) | `Sp-dCas9` CRISPRi | **IGVF `IGVFDS9825FOBP`** (Gogate et al., *bioRxiv* Sep 2026, `PMC13629240`) | Open IGVF External Policy (preprint `CC-BY-NC 4.0`) | **`50.999 GB`** (`4.925 GB` `.h5` + `46.074 GB` `.pkl` guide + `0.24 MB` CSV; `.h5ad` is `87.43 GB` > 60 GB) | **`549,727`** QC cells; **`59,714` NTC cells** (`600` NTC sgRNAs) | **`1,408`** (`1,394` mapped targeting / `1,508` w/ pos ctrl) | **`244`** tgt (`251` w/ pos ctrl) | **Yes** (`24` HTO-superloaded libraries `LW558..LW581`) |
| **5** | **Calu-3** (Human lung adenocarcinoma, epithelial) | `dCas9-KRAB` CRISPRi | **Zenodo `13350497`** (`SunshineHein2023.h5ad`) / `GSE208240` (Sunshine et al., *Nat Commun* 2023) | **`CC BY 4.0`** | **`743.67 MB`** (`0.744 GB` `.h5ad`; `1.2 GB` GEO `.tar.gz`) | **`27,882`** single-element QC cells (`2,500` NTC cells, `22,619` single-gene cells) | **`195`** single-gene (`183` unique single genes + `22` dual pairs) | **`63`** | **Verify on ingestion** (check 10x channel/bystander `.obs` cols) |

---

### Detailed Candidate Profiles (Tier A)

#### 1. DLD-1 (`GSE337988`) — **Top Cancer-Epithelial Candidate**
- **Cell Line / Modality**: Human colorectal adenocarcinoma **DLD-1** engineered with doxycycline/Shield-1-inducible `Zim3-dCas9` CRISPRi (`Day 5` induction, MULTI-seq CMO cell hashing).
- **Citation**: Yeung J et al., *"Joint analysis of multiply perturbed cells improves statistical power and cost efficiency in Perturb-seq"*, *bioRxiv* 2026 (`doi:10.64898/2026.07.10.737863`, EuropePMC `PPR1277619`).
- **License**: Preprint `CC BY` (`api.biorxiv.org`); open public NCBI GEO deposit (`PRJNA1492064`).
- **URLs & Sizes**:
  - GEO Series: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE337988`
  - FTP root: `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE337nnn/GSE337988/suppl/`
  - Recommended Low-MOI processed RNA/guide HDF5: `GSE337988_sublib2_processed_objects_Low_assays.h5` = **`4,621,149,400 B` (`4.621 GB`)** (+ `GSE337988_sublib2_processed_objects_Low_se.rds` = `670,285,340 B` / `0.670 GB`).
  - Optional Medium/High MOI objects: `Med_assays.h5` = `1.839 GB` (`Med_se.rds` = `1.017 GB`), `High_assays.h5` = `1.212 GB` (`High_se.rds` = `3.116 GB`); complete `GSE337988_RAW.tar` = `26,495,651,840 B` (`26.496 GB`).
  - Guide reference (`0.22 MB`): `GSE337988_NGS6475_crispr_feature_barcode_MENGNTCUPDATE.csv.gz` (`23,480` sgRNAs).
- **Cells & NTC**: From all 48 GSM `metrics_summary.csv.gz` reports:
  - `Low` MOI (`~0.3`, `24` 10x channels): **`1,338,060` estimated GEX cells** (low-MOI regime minimizes multi-guide confounding!).
  - `Med` (`~3`, `12` channels): `477,773`; `High` (`~10`, `8` channels): `336,170`; `Pilot` (`4` channels): `125,616`.
  - Feature library contains **`504` non-targeting (`NTC_*`) sgRNAs** (`~2.15%` of library). Exact QC-filtered NTC cell count inside `Low_assays.h5` will be measured during ingestion.
- **Target Overlap (`GSE337988_NGS6475_crispr_feature_barcode_MENGNTCUPDATE.csv.gz`)**:
  - **`4,954`** unique target genes (`504` NTCs excluded).
  - **`shared_union6 = 668`**: `nadig_hepg2`: **`586`**, `nadig_jurkat`: **`586`**, `replogle_rpe1`: **`589`**, `replogle_k562_essential`: **`512`**, `vcc2025_h1`: **`85`**, `nourreddine_kolf21j`: **`26`**.
- **Batch Info**: **Yes** — `44` sublib2 10x Chromium channels (`24` Low MOI across sample pools `SAM24482910`/`11`/`12` + `96`-barcode `GSE337988_MULTIseq_hashing.csv.gz`). Clean channel-level or MULTI-seq-hashed 50/50 split-half validation is directly supported.

#### 2. Primary Human CD4+ αβ T Cells (`GSE314342` + S3 `marson2025_data`) — **Top Primary Normal Immune Candidate**
- **Cell Line / Modality**: Primary human CD4+ αβ T cells from **4 healthy donors** (`CE0008162`, `CE0010866`, `CE0008678`, `CE0006864`), genome-scale `ZIM3KRAB-dCas9` CRISPRi (`10x GEM-X Flex v2`, 3 conditions: `Rest`, `Stim8hr`, `Stim48hr`).
- **Citation**: Zhu X et al., *Cell* (2026-10-29), `doi:10.1016/j.cell.2026.08.002`, `PMID:42664972`, `PMC13614707`.
- **License**: **`https://mit-license.org/` (`MIT License`)** verified in S3 Croissant metadata (`marson_D1_Rest.assigned_guide.jsonld`) hosted on CZ Virtual Cell Models (`https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq`); paired open GEO `GSE314342`.
- **URLs & `<= 60 GB` Options**:
  - GEO Series: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314342` (`Rest`: `94` `.h5` lane matrices = **`54.582 GB`** [`50.834 GiB`]; `Stim8hr`: `94` `.h5` = `58.324 GB`; `Stim48hr`: `96` `.h5` = `58.247 GB`; full `RAW.tar` = `171.154 GB` exceeds budget).
  - S3 Root: `https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/`
    - Pre-aggregated pseudobulk across all 3 conditions: `GWCD4i.pseudobulk_merged.h5ad` = **`44,566,657,140 B` (`44.567 GB <= 60 GB`)** (guide × donor × condition pseudobulk UMI sums). Note: per-donor raw `.assigned_guide.h5ad` files on S3 are `118.610–172.796 GB` each (exceed `60 GB`).
    - Essential S3 metadata tables (`~10 MB` total): `suppl_tables/sgrna_library_metadata.suppl_table.csv` (`9.94 MB`), `suppl_tables/sample_metadata.suppl_table.csv` (`2,966 B`), `suppl_tables/QC_summaries_per_sample_lane.csv` (`60.5 KB`).
  - **Practical Budget Note**: Downloading a **2-donor subset of `Rest` `.h5` lanes** (~`46–48` `.h5` files ≈ `27.3 GB`) leaves >`30 GB` headroom to co-ingest `DLD-1` + `WTC11 i3Neurons` + `Calu-3` inside `60 GB`, OR all 4 `Rest` donors (`54.582 GB`) + `DLD-1` (`4.621 GB`) + `WTC11 i3Neurons` (`0.289 GB`) = `59.492 GB <= 60 GB`.
- **Cells & NTC (`QC_summaries_per_sample_lane.csv`, `Rest` = `94` lanes)**:
  - `Rest` total: **`11,487,117` cells** (`7,072,306` targeting single-sgRNA, **`313,989` NTC single-sgRNA**, `2,131,212` multi-sgRNA, `25,752` low-quality).
  - Per-donor `Rest` NTC / targeting single-sgRNA: `CE0008162` (`23` lanes: `76,543` NTC / `1,675,494` tgt), `CE0010866` (`23` lanes: `83,987` NTC / `1,899,470` tgt), `CE0008678` (`24` lanes: `76,823` NTC / `1,824,887` tgt), `CE0006864` (`24` lanes: `76,636` NTC / `1,672,455` tgt).
- **Target Overlap**:
  - `sgrna_library_metadata.suppl_table.csv` (`26,504` sgRNAs): **`12,778`** unique target genes + `994` NTC sgRNAs -> **`shared_union6 = 2,482`** (`nadig_hepg2`: `2,189`, `nadig_jurkat`: `2,195`, `replogle_rpe1`: `2,196`, `replogle_k562_essential`: `1,908`, `vcc2025_h1`: `285`, `nourreddine_kolf21j`: `90`).
  - `DE_stats.suppl_table.csv` (`Rest` DE-tested, after T-cell essential gene dropout/coverage filtering): **`11,287`** unique `Rest` targets -> **`shared_union6 = 1,452`** (`nadig_hepg2`: `1,207`, `nadig_jurkat`: `1,215`, `replogle_rpe1`: `1,209`, `replogle_k562_essential`: `939`, `vcc2025_h1`: `262`, `nourreddine_kolf21j`: `79`).
- **Batch Info**: **Yes** — 4 unrelated human donors (`D1..D4`), 2 10x processing batches (`CD4i_R1`, `CD4i_R2`), and 94 lanes (`lane_id`), ideal for donor-stratified and lane-stratified split-half checks.

#### 3. WTC11 hiPSC-Derived `i3Neurons` (`TianKampmann2021_CRISPRi.h5ad` / `GSE152988`) — **Top Compact Neuron Line**
- **Cell Line / Modality**: Differentiated excitatory human neurons (`i3Neurons`) derived from **WTC11** hiPSCs (different genetic background & terminally differentiated neuronal state vs `KOLF2.1J` undifferentiated/CM/ME and `H1` hESCs), `Sp-dCas9-KRAB` CRISPRi CROP-seq (`MOI 0.1–0.2`). *(Use CRISPRi subset strictly; exclude `CRISPRa`.)*
- **Citation**: Tian R, Kampmann M et al., *Nature Neuroscience* (2021), `doi:10.1038/s41593-021-00862-0`, `PMID:34031600`, `PMC8254803`.
- **License**: **`CC BY 4.0`** (Zenodo `13350497` scPerturb harmonization); raw 10x `.h5` open on GEO `GSE152988`.
- **URLs & Sizes**:
  - Zenodo `.h5ad` (recommended): `https://zenodo.org/api/records/13350497/files/TianKampmann2021_CRISPRi.h5ad/content` = **`289,389,479 B` (`289.39 MB` / `0.289 GB`)**.
  - GEO 4 CRISPRi 10x lanes (`GSM4632022/23/25/26`): `79.94 MB + 77.81 MB + 69.49 MB + 88.29 MB` = **`315.54 MB`** (`0.316 GB`).
- **Cells & NTC**: From `MOESM5_ESM.xlsx` (`CRISPRi_CROPseq_cell_counts`, `372` sgRNA entries): **`27,189` single-guide QC neurons**, **`368` NTC neurons** across `6` NTC sgRNAs (`~98,000` loaded neurons into `4` 10x lanes).
- **Target Overlap (`MOESM6_ESM.xlsx`, `CRISPRi_CROP-seq`, excluding `'non-targeting'`)**:
  - **`184`** target genes (`2` sgRNAs/gene + `6` NTCs = `374` guides).
  - **`shared_union6 = 70 >= 50`**: `nadig_hepg2`: **`60`**, `nadig_jurkat`: **`60`**, `replogle_rpe1`: **`61`**, `replogle_k562_essential`: **`53`**, `vcc2025_h1`: **`11`**, `nourreddine_kolf21j`: **`3`**.
- **Batch Info**: **Yes** — **`4` 10x Chromium lanes** (`GSM4632022`, `GSM4632023`, `GSM4632025`, `GSM4632026`, `lib1..lib4`).

#### 4. H9 Female hESC-Derived Immature Cortical Neurons (`IGVFDS9825FOBP`)
- **Cell Line / Modality**: **H9 (`WA09`)** female hESC-derived Day-6 immature cortical neurons (distinct donor/sex vs `H1` male hESC and neuronal lineage vs pluripotent `H1`/`KOLF2.1J`), `Sp-dCas9` CRISPRi targeting promoters (`genomic_element = promoter`).
- **Citation**: Gogate A, Chahrour M, Hon G et al., *bioRxiv* (2026-09-02), `doi:10.64898/2026.09.02.748925`, `PMID:42827575`, `PMC13629240`.
- **License & Access**: IGVF Data Portal (`controlled_access: False/None` — **open public** under IGVF External Resource Sharing Policy permitting free download, analysis, and academic publication); note *bioRxiv* preprint carries `CC-BY-NC 4.0` non-commercial license.
- **URLs & Sizes**:
  - AnalysisSet: `https://api.data.igvf.org/analysis-sets/IGVFDS9825FOBP/?format=json`
  - **Fits `<= 60 GB` via split HDF5 + PKL + CSV**:
    - RNA matrix (`IGVFFI6916RLYF.h5`): `https://api.data.igvf.org/matrix-files/IGVFFI6916RLYF/@@download/IGVFFI6916RLYF.h5` = **`4.925 GB`**
    - Guide matrix (`IGVFFI9448NBQC.pkl`): `https://api.data.igvf.org/matrix-files/IGVFFI9448NBQC/@@download/IGVFFI9448NBQC.pkl` = **`46.074 GB`**
    - Guide table (`IGVFFI7085MUTI.csv.gz`): `https://api.data.igvf.org/tabular-files/IGVFFI7085MUTI/@@download/IGVFFI7085MUTI.csv.gz` = `244,119 B` (`0.24 MB`)
    - Combined split bundle = **`50.999 GB` (`<= 60 GB`)**. *(Note: single `IGVFFI4730IQLW.h5ad` is **`87.425 GB > 60 GB`**.)*
- **Cells & NTC**: **`549,727` QC cells** (`18,085` NPCs, `28,493` intermediate, `13,136` immature neurons among NTCs = **`59,714` NTC cells** across `600` NTC sgRNAs + `588` negative-control guides).
- **Target Overlap (`IGVFFI7085MUTI.csv.gz`, mapped `ENSG -> HGNC`)**:
  - `1,408` ASD risk genes (`1,518` unique `ENSG` incl. pos ctrls -> `1,508` mapped HGNC; `1,394` mapped `type=='targeting'` symbols).
  - **`shared_union6 = 244` (`targeting`-only) / `251` (`targeting + positive control`)**:
    - `targeting`-only: `nadig_hepg2`: **`191`**, `nadig_jurkat`: **`193`**, `replogle_rpe1`: **`192`**, `replogle_k562_essential`: **`158`**, `vcc2025_h1`: **`41`**, `nourreddine_kolf21j`: **`28`**.
- **Batch Info & Technical Caveats**:
  - **Batch**: **Yes** — **`24` HTO-superloaded 10x libraries (`LW558..LW581`)**.
  - **High-MOI & Format Caveats**: Infected at `MOI ~4` (**mean `3.8` sgRNAs/cell**); single-sgRNA-per-cell fraction is small (~`5%`), requiring multi-guide linear/additive deconvolution during pseudobulk construction unless paired with `DLD-1` `High/Med` validation. Also consumes `51.00 GB` (`85%` of `60 GB` budget) and relies on Python `.pkl` for guide assignments.

#### 5. Calu-3 (`SunshineHein2023.h5ad` / `GSE208240`) — **Compact Reserve Epithelial Line (Condition Caveat)**
- **Cell Line / Modality**: Human lung adenocarcinoma epithelial line (**Calu-3**), `dCas9-KRAB` CRISPRi direct-capture Perturb-seq (`183` single host factor genes across `195` single-gene elements, `22` dual-gene pairs, `22` NTC pairs = `239` elements).
- **Citation**: Sunshine S, Hein MY, Weissman JS et al., *Nature Communications* (2023), `doi:10.1038/s41467-023-41788-4`, `PMID:37803001`, `PMC10558542`.
- **License**: **`CC BY 4.0`** (both *Nat. Commun.* PMC10558542 and Zenodo `13350497`).
- **URLs & Size**:
  - Zenodo `.h5ad`: `https://zenodo.org/api/records/13350497/files/SunshineHein2023.h5ad/content` = **`743,672,579 B` (`743.67 MB` / `0.744 GB`)**.
  - GEO: `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE208nnn/GSE208240/suppl/GSE208240_CRISPRi_perturbseq_sarscov2_filtered.tar.gz` (`1.2 GB`).
- **Cells & NTC (`MOESM4_ESM.xlsx`)**: **`27,882` QC single-element cells** (`22,619` single-gene cells; **`2,500` NTC cells** across `22` `non_targeting_pair_*` controls; rest dual-gene).
- **Target Overlap (`MOESM4_ESM.xlsx` strict single-gene elements)**:
  - **`183` unique single-gene targets** (`195` elements; excl. `22` dual-gene pairs).
  - **`shared_union6 = 63 >= 50`**: `nadig_hepg2`: **`56`**, `nadig_jurkat`: **`56`**, `replogle_rpe1`: **`58`**, `replogle_k562_essential`: **`51`**, `vcc2025_h1`: **`7`**, `nourreddine_kolf21j`: **`1`**.
- **Caveats (Record Transparently)**:
  1. **Viral Challenge State**: All cells were challenged with SARS-CoV-2 for 24h (`MOI 4`) and split post-hoc via viral RNA classifier into `infected` vs `uninfected bystander` states (plus borderline cluster `T`). Thus even `bystander` cells experienced paracrine/cytokine exposure rather than naive baseline culture.
  2. **Batch Column Check Needed**: The `.obs` columns of `SunshineHein2023.h5ad` should be verified for 10x channel/replicate labels (`GSE208240` combined `gex`/`guide` series entries; verify lane/GEM well metadata inside `.h5ad` before split-half).

---

## 4. Tier B — Auxiliary / Near-Miss Candidates (Useful with Specific Limitations)

| Dataset / Line | Modality & Size | Targets / `shared_union6` | Cells / NTC | Why Tier B (Specific Limitation) |
|---|---|---|---|---|
| **IGVF `IGVFDS3530HHEU`** (`GM25256 WTC-11` Day-1 **Endothelial Cells**, Engreitz Lab) | `Sp-dCas9-KOX1-KRAB` CC-Perturb-seq; `IGVFFI4783YGLZ.h5ad` = **`143.57 MB` (`0.144 GB`)**, open IGVF | **`298`** promoter genes; **`shared_union6 = 156`** (`HepG2: 136, Jurkat: 139, K562: 123, RPE1: 142, H1: 24, KOLF2.1J: 5`) | `16,925` QC cells (`d1_rep1: 5,856`, `d1_rep2: 11,069`); `15,962` guide-assigned | **High-MOI NTC Bottleneck**: ~6.2 guides/cell leaves ONLY **`64` pure-NTC cells** (`62` single-NTC, `75` NTC+safe-targeting) and `1,930` single-guide cells. High NTC baseline noise if not modeled via multi-guide regression. |
| **IGVF `IGVFDS6332VCTO`** (`WTC-11` Day-12 **Cardiomyocytes**, Sivakumar 2025, `PMC11785042`) | `Sp-dCas9` CRISPRi (`promoter` TF screen); Open IGVF | **`1,951`** promoter targets (`1,962` w/ pos ctrl) + `600` NTCs; **`shared_union6 = 331`** (`335` w/ pos ctrl) | `28` HTO-superloaded libraries (`LW310..`) | **Exceeds `60 GB` Limit**: `IGVFFI7637STPX.h5ad` = **`112.804 GB`**; split `.h5` (`4.637 GB`) + `.pkl` (`71.136 GB`) = **`75.773 GB > 60 GB`**. (Smaller released subpools `IGVFDS4003HZAB` `20.57 GB`, `IGVFDS9192YMGI` `18.63 GB`, and H9 CM `IGVFDS2266YDVM` `18.83 GB` / `IGVFDS6924DJAZ` `15.54 GB` use `CM_Pilot02` / `Cardiac_Pilot02` mostly targeting enhancers/193 pilot loci.) |
| **Feng et al. 2026** (`45 HipSci` Human **iPSC** Lines, `PMC12903452`, `Figshare 27989294` + `26819743`) | `dCas9-KRAB-MeCP2` CRISPRi; **`MIT License`**; **`7.93 GB`** (`7.676 GB` RNA UMI CSVs + `0.229 GB` guide + `0.021 GB` metadata) | **`7,234`** targets (`34–45` lines); **`shared_union6 = 2,375`** (`HepG2: 2,116, Jurkat: 2,119, K562: 1,907, RPE1: 2,121, H1: 270, KOLF2.1J: 73`) | `1,881,068` raw rows -> **`854,276` assigned cells** (`845,924` tgt + **`8,352` NTC** across `286` batches) | **Cell-State / Donor Overlap**: Undifferentiated **iPSCs** (same pluripotent cell state as `nourreddine_kolf21j`), and explicitly contains **`kolf_2` (`120,699` raw cells)** & **`kolf_3` (`118,138` raw cells)**. Suitable only as a multi-donor iPSC robustness benchmark after holding out `kolf_2/3`, NOT as a novel differentiated tissue. |

---

## 5. Recommended Round 10 Ingestion Portfolios (`<= 60 GB` Ceiling)

To maximize biological diversity (**primary immune + colorectal carcinoma + differentiated neuron [+ lung epithelium]**) cleanly within **`60 GB`** while minimizing multi-guide deconvolution risk:

1. **Portfolio A (All-Low/Single-Guide Clean Core + Pre-Aggregated or 2-Donor T Cells — `<= 50.3 GB`)**:
   - **`DLD-1` (`GSE337988` `Low_assays.h5` + `.rds` + guide CSV)**: **`5.29 GB`** (`668` shared targets, low-MOI colon carcinoma)
   - **`WTC11 i3Neurons` (`TianKampmann2021_CRISPRi.h5ad`, Zenodo)**: **`0.29 GB`** (`70` shared targets, low-MOI differentiated neuron)
   - **`Calu-3` (`SunshineHein2023.h5ad`, Zenodo, optional w/ infection covariate)**: **`0.74 GB`** (`63` shared targets, low-MOI lung epithelial)
   - **`Primary CD4+ T Cells` (`GSE314342` `Rest`)**: **2-donor subset (`46–48` `.h5` lanes ≈ `27.3 GB`)** OR S3 `GWCD4i.pseudobulk_merged.h5ad` (**`44.57 GB`**) (`1,452` `Rest` DE / `2,482` lib shared targets, primary normal immune lineage)
   - **Total**: **`~33.6 GB`** (2-donor T `Rest` scRNA + DLD-1 + i3Neurons + Calu-3) or **`~50.2 GB`** (S3 pseudobulk + DLD-1 + i3Neurons + Calu-3) — strictly `<= 60 GB` with comfortable staging margin.
2. **Portfolio B (Maximum Single-Cell Sample Size — `59.49 GB <= 60 GB`)**:
   - `Primary CD4+ T Cells` ALL 4 donors (`GSE314342` `Rest`, 94 `.h5` files = `54.582 GB`) + `DLD-1` (`Low_assays.h5` = `4.621 GB`) + `WTC11 i3Neurons` (`0.289 GB`) = **`59.492 GB`** (requires sequential donor streaming to disk so peak raw + processed stays well below quota).

---

## 6. Tier C — Verified Excluded Candidates (Do Not Re-Screen)

| Candidate / Accession | Cell Line / Modality | Verified Verification Metric | Exact Reason for Exclusion |
|---|---|---|---|
| **X-Atlas / Orion** | HCT116, HEK293T, etc. | — | **Explicitly prohibited** by project policy |
| **`GSE337804`** | A549, HCT116 (`219` targets, `shared=75`) | GEO Sample title typo `"CRISPIRi"`; `GSE337804_2cl_description.pdf` states *"low-MOI CRISPRn (= CRISPR knockout)"* | **Wrong modality** (Cas9 KO, NOT CRISPRi) |
| **`GSE273677`** (*Nat Metab* 2026) | `EndoC-betaH1` pancreatic β-cells (`61` T2D + `40` RQC) | `MOESM3_ESM.xlsx` union (`126` symbols ST1/2/14/15) -> **`shared_union6 = 20`** | **`< 50` shared targets** |
| **`GSE182308`** (Leng et al. 2022) | hiPSC-derived `iAstrocytes` (`485 MB`) | `GSM5526487_sgRNAcounts.csv.gz` -> **`31` targets**, `2` NTCs -> **`shared_union6 = 6`** | **`< 50` targets & `< 50` shared** |
| **`GSE289235`** (Boggess et al. 2025) | `WTC11` Day-21 Neurons (`1.01 GB`) | FACS screen had `1,343` genes, but `PMC11844385` CROP-seq targeted ONLY **`45` genes** | **`< 50` targets** |
| **`GSE178317`** (Dräger et al. 2022) | `WTC11` iTF-Microglia | Survival screen had `2,325` genes, but `PMC9448678` CROP-seq targeted ONLY **`39` genes** | **`< 50` targets** |
| **`GSE335887`** (2026) | `iMG` / `iTF` Microglia | GEO `feature_reference.csv.gz` -> **`30` targets** | **`< 50` targets** |
| **`GSE284197`** (2025) | Cortical Radial Glia | **`44` TF targets** | **`< 50` targets** |
| **`GSE272093`** (2024) | Cortical iAssembloids | GEO `sgRNA_assignments.txt.gz` -> `68` targets -> **`shared_union6 = 11`** | **`< 50` shared targets** |
| **`GSE238219`** (`PMC11865357`, 2025) | `HepaRG` differentiated hepatocytes (`644 MB`) | Table 3 Perturb-seq panel targeted ONLY **`11` genes** | **`< 50` targets** |
| **`GSE217812`** (`PMC10635647`, 2023) | `SGBS` human adipocytes (`KRAB-dCas9`) | CROP-seq targeted **`6` candidate + `6` positive control = `12` genes** | **`< 50` targets** |
| **`PMC10030154`** (Wessels et al. 2023) | `THP-1` AML (`CaRPool-seq` Cas13d, `219 MB`) | Bulk FACS was `439` genes, but scRNA-seq `CaRPool-seq` had **`28` single perturbations** (`26` + `2` neg ctrls) | **`< 50` targets** (plus combinatorial focus) |
| **`GSE293514`** (`PMC12808753`, 2026) | Human myoblasts (`MyoCRISPR-KO Lib`) | Split-toxin screen uses **Cas9 knockout (`KO`)**, NOT CRISPRi | **Wrong modality** (CRISPR KO) |
| **`GSE343369`** (UCSF Study 2) | `KOLF2.2J` iPSCs & Embryoid Bodies | GEO `perturbed_gene_list.csv.gz` -> `223` targets -> **`shared_union6 = 34`** (+ parental line of `KOLF2.1J`) | **`< 50` shared targets** & donor overlap |
| **`IGVFDS9629VXQK`** (Engreitz Lab) | `WTC-11` Endothelial GW MACS screen | `assay = protein abundance / gRNA sequencing` (bulk FACS/MACS sort, NOT scRNA-seq) | **Not scRNA-seq Perturb-seq** |
| **`CZ CellxGene dde06e0f...`** | PBMCs (OneK1K, cited in Zhu 2026) | Population sc-eQTL cohort (Yazar et al.), zero CRISPR perturbations | **Not Perturb-seq** |
