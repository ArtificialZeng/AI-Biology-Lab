# Round 12 Prospective Validation Data Screening (`A-T1`, Metadata-Only)

**Cutoff Date**: `2026-11-22` (strictly metadata-only audit: **zero** `.h5`, `.h5ad`, `.mtx`, large `.rds`, or `.tar` expression matrix files downloaded; **zero** account registrations, terms-clickthroughs, `X-Atlas/Orion` datasets, `Arc Virtual Cell Challenge [VCC]` datasets, or HTTP 401/403 bypasses).

---

## 1. Screening Criteria `(a)–(e)` & Yield Calibration Basis

Every candidate was evaluated against Round 12 eligibility criteria **(a)–(e)**:
- **(a)** Public single-cell **CRISPRi or CRISPRko** Perturb-seq in human cells with **1 targeting sgRNA per cell** (explicitly filtered or determinable via low MOI / guide-purity metadata) and **non-targeting control (NTC)** cells.
- **(b)** Cell line or primary/differentiated cell type **NOT used in Rounds 9–11** (Excluded lines/types: `K562`, `RPE1`, `HepG2`, `Jurkat`, `H1`, `KOLF2.1J`, `DLD1`, `Calu3`, `CD4+ T [resting & stimulated]`, `HipSci iPSC`, `TeloHAEC`, `HCASMC`, `WTC11 i3Neuron / iPSC`).
- **(c)** Overlap with frozen Round 9 reproducible target union (**`union6_rep = 2,117`** canonical HGNC symbols; **`R9_COMMON = 6,148`**; **`union6_rep ∩ R9_COMMON = 1,698`**) sufficient to yield **`n_eval_targets >= 50`** (plus `|G_eval(t)| >= 3,000` measured genes in `R9_COMMON`).
  - **Empirical Yield Calibration (R10/R11)**: Across evaluated lines (`round11/shared/R11_QC_COUNTS.md`), only **`~20%–40%`** of shared `u6_rep` targets survive single-line reproducibility QC + self-gene presence in `G_eval(t)`:
    - R11 `cd4t_stim48` (`PASS`): `503` `n_eval_targets` out of `2,117` `u6_rep` (**`23.8%`** of `u6_rep`; **`29.6%`** of `1,698` `u6_rep ∩ R9_COMMON`).
    - R11 `feng_ipsc` (`SECONDARY_ONLY`): `93` `n_eval_targets` out of `444` deduplicated library targets (**`20.9%`** of total targets).
    - R11 `telohaec` (`FAIL`): `2,198` deduplicated targets yielded only `17` line-reproducible targets and **`9`** `n_eval_targets` (`< 50`) despite `|G_eval| = 4,869`.
    - R11 `hcasmc` (`EXCLUDED_FEASIBILITY_Q3`): Parse flat-root `.h5` contained `0` per-cell guide assignments despite `955,600` cells and separate guide-library TSVs — highlighting the mandatory need to verify **per-cell** guide linkage before committing to any matrix download.
  - **Calibrated Screening Gate**: Under the `20%–40%` survival rule, a dataset needs **`u6_rep` intersection >= `125–250` (`u6_rep ∩ R9_COMMON >= 125–250`)** to reliably clear `n_eval_targets >= 50`, and **`u6_rep <= 50`** is mathematically impossible even at `100%` survival.
- **(d)** Open anonymous download (`HTTP 200 OK`), NO `X-Atlas/Orion`, NO terms/login/account-gated data (`VCC`, Broad SCP API `401/500`, `dbGaP`).
- **(e)** Publicly released prior to **`2026-11-22`**.

---

## 2. Primary & Borderline Candidate Summary Table

| Tier | Accession / Source | Cell Type / Line | Modality | Release Date | # Targets (`Raw / Canon`) | # Cells (`Total / QC`) | Guide Assignment & NTCs | `u6_rep` (`∩ R9C`) | Est. `n_eval_targets` (`20–40%` rule) | File Size(s) & Direct Access | Verdict |
|---|---|---|---|---:|---:|---:|---|---:|---:|---|---|
| **Tier 1 (Strong Qualifier)** | **`GSE345058`** (bioRxiv `10.64898/2026.06.01.728087`) | **`A549`** (`A549-doxCas9`) | **CRISPRko** (`CROP-seq`, Dox-Cas9) | `2026-06-02` (preprint) / `2026-08-26` (GEO) | **`1,000`** (`1,000` canon; `4` sgRNAs/gene; `1,000/1,000 >=50` cells & `half_valid`) | **`606,075` QC singlets** (`560,755` Perturbed + **`45,320` NTC**; `med 565`, `min 73` cells/gene) | **YES (exact)**: `GSE345058_cell_annotations.csv.gz` (`606,075 / 606,075` single guide, `0` multi-guide delimiters; `MOI=0.1` + Scrublet); **`211` NTC sgRNAs** (`10,451–11,919` NTC/batch across `4` batches) | **`530`** (**`452`** in `R9_COMMON`) | **`106–212`** (`90–181` in `R9C`; R11-rate `~126–134`) | `.h5ad` `23.3 GiB` (`25.05 GB`) or `.mtx.gz` `10G` + `norm.mtx.gz` `12G` + metadata CSVs `38M` (`HTTP 200` open) | **PASS `(a)–(e)` — Top Primary R12 Line** |
| **Tier 2 (Borderline / Uncertainty Risks)** | **`Zenodo 13350497: XuCao2023.h5ad`** / **`GSE218566`** (`PMID 36778497`, `PMC10961254`) | **`HEK293`** (`HEK293-idCas9`) | **CRISPRi** (`dCas9-KRAB-MeCP2`, `4sU` PerturbSci-Kinetics) | `2023-08-10` (GEO) / `2024-08-26` (Zenodo `v`) | **`228`** (`228` canon; `3` sgRNAs/gene; `181` targets in author DEG ST10 w/ `>=50` cells & `>40%` KD) | **`126,271` sgRNA singlets** (`161,966` matched; **`98,315`** retained post `>40%` KD filter, `med 484` cells/gene) | **YES (explicit)**: `MOI 0.1–0.2`, singlet rule (`top >=60% & >=3x 2nd sgRNA`); **`15` NTC sgRNAs** (~2.15% of library, exact cell count unverified w/o opening `.h5ad`) | **`164`** (**`143`** in `R9_COMMON`; ST10 `121` / `109`) | **`33–66`** (`29–57` in `R9C`) — **straddles `50` floor** | `XuCao2023.h5ad` **`336.64 MB`** (`CC-BY-4.0`, `HTTP 200` open); GEO `RAW.tar` `516M` | **BORDERLINE `(c)` + Shallow Depth Risk** (`|G_eval|>=3000` uncertain) |
| **Tier 2 (Likely Underpowered `(c)`)** | **`GSE327644`** (`PMID 42608566`, `10.1038/s41587-026-03258-2`) | **Primary Human Macrophages** (`CD14+` monocyte-derived, `Ctrl` & `LPS 5h`, donor `RG5581`) | **CRISPRko** (`Cas9-eVLP` + `lentiVPX` CROP-seq; *mislabelled `CRISPRi` in GEO tag*) | `2026-04-10` (GEO) / `2026-08-17` (Nat Biotechnol) | **`569` canon** (`2,332` targeting guides; text cites `583` including aliases/controls) | **`~2.7M Ctrl + ~3.3M LPS` FACS GFP+ loaded** across `8` Flex lanes (`Mac_L01..L08`, `16` `.h5` files; exact post-QC singlets unverified) | **In-matrix (`Cell Ranger v9.0.1`)**, *single-sgRNA MOI/purity NOT printed in GEO protocol*; **`593` Control sgRNAs** (~20% of library) | **`108`** (**`92`** in `R9_COMMON`) | **`22–43`** (`18–37` in `R9C`) — **likely `< 50`** | `GSE327644_RAW.tar` **`6.58 GB`** (`16` filtered `.h5` files, `567M–1.05G` each; `HTTP 200` open) | **HIGH RISK FAIL `(c)` + Flex Probe + Unverified MOI** |

---

## 3. Detailed Profiles of Primary & Borderline Candidates

### 3.1 `GSE345058` — CZ Biohub `A549` Perturbation Atlas (Strong Primary Qualifier)
- **Citation / Accession**: Liu et al., *"A multimodal perturbation atlas defines the phenotypic resolution of cellular morphology"*, bioRxiv `10.64898/2026.06.01.728087` (`2026-06-02`, preprint text `cc_by_nc_nd`), GEO Series **`GSE345058`** (Public **`2026-08-26`**, unrestricted GEO FTP/HTTPS).
- **Cell Type & Modality**: Human **`A549`** lung adenocarcinoma stably expressing doxycycline-inducible Cas9 (**`A549-doxCas9`**), transduced with CROP-seq lentiviral sgRNA library at low **MOI = 0.1** (**CRISPRko**). Fresh/new lineage — `A549` has never been used in Rounds 9–11.
- **Scale & Cell/Target Metadata (Directly Verified from `GSE345058_cell_annotations.csv.gz` & `GSE345058_sgRNA.csv.gz`)**:
  - **Total QC Cells**: **`606,075`** across `4` batches (`batch1`: `163,833`; `batch2`: `145,974`; `batch3`: `154,614`; `batch4`: `141,654`) and `32` 10x lanes.
  - **Targeting vs. Control Cells**: **`560,755`** `Perturbed` cells + **`45,320`** `Nontargeting` NTC cells (`batch1=11,839`, `batch2=11,111`, `batch3=11,919`, `batch4=10,451`), enabling clean disjoint Half A (`batch1+3`: `23,758` NTC) vs Half B (`batch2+4`: `21,562` NTC) partitioning with `>=10,000` NTCs per half.
  - **Single-sgRNA Determinability**: Every single one of `606,075` rows in `cell_annotations.csv.gz` carries a single clean `sgRNA` sequence and single `Perturbation` symbol (`0` multi-guide delimiters `,;|+ `), consistent with low `MOI = 0.1` + Scrublet/multi-sgRNA filtering; guide UMI `barcode_counts` spans `min=3.0, p5=4.0, med=12.0, max=163.0`.
  - **Targets & Cell Depth**: **`1,000`** protein-coding genes (`4` sgRNAs/gene = `4,000` targeting + `211` NTC = `4,211` guides). Perturbed cell counts range from `min=73` (`p5=267`, `p10=348.7`, `med=565.0`, `max=1,196`), with **`1,000 / 1,000` (`100%`) >= 50 cells** and **`1,000 / 1,000` (`100%`) `half_valid`** (`>=25` cells in both `batch1+3` and `batch2+4`).
  - **Transcriptome Depth & Gene Coverage**: `n_genes_by_counts`: `min=3,007`, `p5=4,435`, **`med=6,083`**, `max=11,451`; `total_counts`: `min=8,129`, `p5=14,911`, **`med=26,103`**; `pct_counts_mt`: `med=3.5%`, `max=10.0%`. `GSE345058_gene_annotations.csv.gz` (`29,892` features) covers **`6,133 / 6,148` (`99.76%`)** of `R9_COMMON` (`2,112 / 2,117` of `u6_rep`), guaranteeing `|G_eval(A549)| >> 3,000`.
- **Target Overlap & Estimated `n_eval_targets`**:
  - **`u6_rep = 530 / 1,000` (`53.0%`)**; **`u6_rep ∩ R9_COMMON = 452`** (`452/452 >= 50` cells & `half_valid`).
  - Calibrated at `20%–40%`: **`106–212` `eval_targets`** (`90–181` in `R9_COMMON`; at R11's `cd4t_stim48` survival rate of `23.8%` of `u6_rep` / `29.6%` of `u6_rep ∩ R9_COMMON`, expected yield is **`~126–134` `eval_targets`**).
- **Files & Sizes**:
  - `GSE345058_SC_raw_normalized_counts.h5ad` (`25.05 GB` / `23G`), OR raw sparse triplet `GSE345058_raw_counts.mtx.gz` (`10G`) + `GSE345058_normalized_counts.mtx.gz` (`12G`) + `GSE345058_cell_annotations.csv.gz` (`38M`) + `GSE345058_gene_annotations.csv.gz` (`370K`) + `GSE345058_sgRNA.csv.gz` (`67K`).
- **Explicit Uncertainties / Caveats**:
  1. **`CAS9` Transcript Zeroes**: In `cell_annotations.csv.gz`, **`161,953 / 606,075` (`26.7%`)** cells (`11,635 / 45,320` = `25.7%` of NTCs) have `total_counts_CAS9 == 0`, reflecting either dox-Cas9 expression variegation or mRNA capture dropout; if biological silencing occurred in a fraction of cells, effect magnitudes (and thus AD/Energy detectability for weaker hits) could be diluted, though `565` median cells/gene provides strong statistical power.
  2. **Preprint License Metadata**: BioRxiv preprint metadata lists **`cc_by_nc_nd`** for the manuscript; GEO hosts the data without explicit restrictions, though non-commercial terms on any preprint-bundled supplementary artifacts should be noted.
  3. **Memory / Storage Footprint**: Downloading `GSE345058_raw_counts.mtx.gz` (`10G`) + `GSE345058_cell_annotations.csv.gz` (`38M`) + `GSE345058_gene_annotations.csv.gz` (`370K`) + `GSE345058_sgRNA.csv.gz` (`67K`) avoids redundant `normalized_counts` in the `25.05 GB` (`23G`) `.h5ad`, while still requiring chunked sparse CSR streaming upon ingest.

---

### 3.2 `XuCao2023.h5ad` / `GSE218566` — `HEK293` PerturbSci-Kinetics (Borderline `(c)` + Shallow Depth Risk)
- **Citation / Accession**: Xu et al., *Nature Biotechnology* 2023 (`PMID 36778497`, `PMC10961254`, `10.1038/s41587-023-01948-9`), GEO **`GSE218566`** (`2023-08-10`), harmonized `scPerturb` release **`Zenodo 13350497: XuCao2023.h5ad`** (`2024-08-26`, `CC-BY-4.0`, **`336.64 MB`**, `HTTP 200 OK`).
- **Cell Type & Modality**: Human **`HEK293`** (`HEK293-idCas9`) expressing dox-inducible `dCas9-KRAB-MeCP2` (**CRISPRi**, `7d` dox + `2h` `200 uM 4sU` metabolic RNA labeling). `HEK293` is distinct from all R9–R11 lines.
- **Scale, Guide Assignment, & NTCs (`PMC10961254` & `41587_2023_1948_MOESM3_ESM.xlsx`)**:
  - Transduced across `4` replicates (`2` at `MOI = 0.1`, `2` at `MOI = 0.2`) with **`699` sgRNAs** (`684` targeting **`228` genes** at `3`/gene + **`15` NTCs**).
  - **`161,966` cells** recovered with matched sgRNAs (`88%`); explicit singlet rule (**top sgRNA `>= 60%` and `>= 3x` 2nd sgRNA**, `>= 10` sgRNA UMIs, `>= 300` UMIs, `>= 200` genes) yielded **`126,271` sgRNA singlets**; after authors' post-hoc exclusion of low-efficacy guides (`<= 40%` target KD vs NTC), **`98,315` cells** remained (`median 484 cells/gene`).
- **Target Overlap & Estimated `n_eval_targets`**:
  - Full library (`228` genes): **`u6_rep = 164`** (**`143`** in `R9_COMMON`, `87.2%`!).
  - Author DEG table (`ST10`, targets with `>= 50` cells & `> 40%` KD showing `>= 1` DEG): `181` targets -> **`u6_rep = 121` (`109` in `R9_COMMON`)**.
  - Estimated `n_eval_targets` (`20%–40%` on `u6_rep = 164` / `R9C = 143`): **`33–66`** (`29–57` in `R9C`). Because library selection is strongly enriched for shared RNA-processing/transcription/mitochondrial genes (`71.9%` of targets in `u6_rep`!) and `109–121` already showed transcriptional effects in `HEK293`, survival *might* approach or slightly clear `50` — OR fall just short at `~35–48`.
- **Explicit Uncertainties / Risks**:
  1. **Shallow Sequencing Depth (`|G_eval| >= 3,000` Risk)**: Combinatorial indexing yielded only **`~8,000` reads/cell**, **`median 2,155 UMIs/cell`** (QC floor `300 UMIs / 200 genes`). Low UMI counts per cell can cause excessive gene-detection filtering loss (`|G_eval|` falling near or below `3,000`) and attenuate per-gene perturbation AD/Energy reproducibility.
  2. **Uncertain `.h5ad` Snapshot & NTC Count**: Paper XML does not print the exact NTC cell count (`15 / 699` sgRNAs ~ `2.15%`, implying `~2,100–2,700` NTC cells if unfiltered — close to the `1,000 + 500` disjoint Half A/B NTC budget). Further, `scPerturb`'s `XuCao2023.h5ad` could contain `126,271` singlets or pre-filtered `98,315` high-KD cells, and replicate columns (`MOI 0.1` vs `0.2`) require `.obs` inspection upon `Q3`.
  3. **`4sU` Metabolic Labeling**: `2h` `4sU` pulse + chemical conversion prior to harvest could subtly shift baseline stress/mitochondrial transcripts relative to standard scRNA-seq, though matched within-experiment NTC controls normalize batch shifts.

---

### 3.3 `GSE327644` — Marson Lab Primary Human Macrophages (Borderline / Likely Underpowered `(c)`)
- **Citation / Accession**: *Nature Biotechnology* (`2026-08-17`, `PMID 42608566`, `10.1038/s41587-026-03258-2`), GEO **`GSE327644`** (Public **`2026-04-10`**, updated `2026-11-19`, open `HTTP 200 OK`).
- **Cell Type & Modality**: **Primary human monocyte-derived macrophages** (`donor RG5581`, M-CSF 11d, `Mac_Ctrl` untreated vs `Mac_LPS` `5h 50 ng/ml LPS`, `8` lanes `Mac_L01..L08`, `10x GEM-X Flex v2` probe set `v1.1.0 2024-A`). **Important Correction**: GEO characteristic field reads `"CRISPR interference"`, whereas `PMID 42608566` abstract and `!Sample_extract_protocol_ch1` prove this screen uses **Cas9-eVLP RNP + lentiVPX CROP-seq sgRNA (CRISPRko)**.
- **Scale & Library (`GSE327644_PRJNA1450051_marson_mac_perturbseq_guide_template.txt.gz`)**:
  - `2,925` guides = `2,332` targeting guides (**`569` unique canonical HGNC genes**; GEO protocol states `583` including unmapped/controls) + **`593` `Control` sgRNAs** (`20.3%` of library).
  - Sorted `~2.7M` unstimulated + `~3.3M` LPS GFP+ cells pooled for `8` Flex channels; `GSE327644_RAW.tar` (`6.58 GB`) holds `16` Cell Ranger v9.0.1 `.h5` matrices (`Mac_Ctrl` and `Mac_LPS` per lane, `567 MB – 1.05 GB` each).
- **Target Overlap & Estimated `n_eval_targets`**:
  - **`u6_rep = 108`**, **`u6_rep ∩ R9_COMMON = 92`**.
  - At `20%–40%` survival: **`22–43` `eval_targets`** (**`18–37` in `R9_COMMON`**). At R11 `cd4t_stim48` survival (`29.6%` of `R9C`), expected yield is `~27` (`< 50`). Clearing `50` requires `>= 54.3%` survival across QC + line-reproducibility + `G_eval`.
- **Explicit Uncertainties / Risks**:
  1. **Severe Underpowering Risk under `(c)`**: Expected `n_eval_targets` (`18–37`) sits well below `50`.
  2. **Unverified Per-Cell Singlet Purity / MOI**: Protocol specifies lentiviral volume (`4 ul / 1 ml`) without stating target MOI or `%` single-sgRNA cells; Flex `.h5` guide-assignment matrix and cell count per lane cannot be verified without opening `.h5` (`Q3` feasibility check mandatory if ever pursued).
  3. **10x Flex Targeted Probe Chemistry & Single Donor**: Uses fixed probe hybridization (`~18k` targeted genes) rather than 3'/5' poly-A capture, and all `8` lanes derive from a single donor (`RG5581`), forcing lane-split (`L01–L04 vs L05–L08`) rather than donor-split Half A/B and splitting cells 50/50 between `Ctrl` and `LPS`.

---

## 4. Quantitative Near-Miss Audit (Excluded Primarily by `(c)` `< 50` `n_eval_targets`)

| Dataset | Cell Type / Modality | Release Date | # Targets (`Canon`) | # Cells | `u6_rep` (`∩ R9C`) | Est. `n_eval` / Theoretical Max | Primary Reason Excluded |
|---|---|---|---|---:|---:|---:|---|
| **`GSE337804`** (Myllia Biotech) | `A549` CRISPRko (`MOI ~6%`, DepleteX 5' HT) | `2026-07-10` | `217` (`654` sgRNAs + `20` `ctrl_NT` + `20` `ctrl_T`) | `~72,600` singlets (`2` channels; `1.12 GB`) | **`65`** (`61`) | `13–26` (max `61`) | Fails `(c)` (`u6_rep=65 << ~125`); redundant w/ `GSE345058` (`A549`, `u6_rep=530`) |
| **`GSE281048` / `Zenodo 14518762`** (`Mixscale`, *Nat Cell Biol* `PMC12083445`) | `A549, MCF7, HT29, HAP1, BxPC3, K562` CRISPRi (`Parse Evercode Mega` + `14` NT) | `2024-11-04` / `2025-02-27` (`CC-BY-4.0`) | `IFNG: 59`, `IFNB: 61`, `TGFB: 52`, `TNFA: 55`, `INS: 44`; **Union `ST1f`: `218`** | `~2.6M` total (`6` lines x `5` stimulations; `2.4–5.2G` `.rds.gz` GEO / `2.64–5.60 GB` `.rds` Zenodo) | Per-path `9–20`; **Union `55`** (`46`) | `11–22` (max `46` in R9C even pooled!) | Fails `(c)` (`u6_rep=9–20` per pathway; `46` R9C even across all 5 pathways), PLUS **no unstimulated control** (`24h` cytokine stim in every arm) |
| **`Zenodo 13350497: FrangiehIzar2021_RNA.h5ad`** (`PMC8376399`, `SCP1064` `HTTP/2 500` on files API) | Melanoma `2686` CRISPRko (`MOI=0.1`, CITE-seq) | `2021` / `2024-08-26` (`CC-BY-4.0`) | **`249`** ST1 (`245` ST2a; only `141` Ctrl / `183` IFNγ / `181` co-culture analyzed) | `218,331` (`57,627` Ctrl, `87,590` IFNγ, `73,114` TIL; `1.46 GB`) | **`51`** (`44`) | `10–20` (max `44` in R9C) | Fails `(c)` (`u6_rep=51`, `R9C=44` strictly `< 50` even at 100% survival) |
| **`GSE291147`** (`PerturbFate`, *Nature* Apr 2026) | `A375` melanoma CRISPRi (`MOI~0.05`, dual-sgRNA + `DMSO/PLX`) | `2026-02-06` (GEO) / `Apr 2026` | **`144`** (`143` gene pairs + NT) | `~340k` (`Dual` + `Tri`, `302–663M` `.RDS`) | **`39`** (`34`) | `8–16` (max `34`) | Fails `(c)` (`u6_rep=39 < 50`) |
| **`GSE263747`** (*Cell Rep* 2024, `PMID 38691456`) | `SNU-761` HCC CRISPRko (`CROP-seq`, D5 & D10) | `2024-04-29` | **`108`** TF genes (`450` sgRNAs + NTC) | `556 MB` tar (`D5+D10+HTO`) | Unlisted (`108` total TFs `< ~125`) | `<< 50` | Fails `(c)` (`108` total TF library cannot reach `50` at `20–40%` survival) |
| **`GSE273677`** (*Nat Metab* 2026, `PMID 41482566`) | `EndoC-βH1` pancreatic β-cell CRISPRi (`MOI<20%`, 3' DC) | `2025-09-22` | **`101`** (`61` T2D + `40` RQC; `45+25` survived KD/GSEA) | `1.4 GB` tar | **`15`** (`14`) (`1` T2D + `14` RQC) | `3–6` (max `14`) | Fails `(c)` (`u6_rep=15 << 50`) |
| **`GSE283574`** (*Nat Commun* 2026, `PMC13392366`) | `H358` Cas12a-dTag CRISPRi (EpiGenFactor + Surface) | `2026-04-07` | `55` single + `17` dual | `299 MB` tar | **`29`** (`28`) | `6–12` (max `28`) | Fails `(c)` (`u6_rep=29 < 50`) |
| **`GSE280506`** (`PMC12354852`) | Human gastric organoids CRISPR | `2025-06-26` | `78` single-sgRNA genes | `—` | **`41`** (`34`) | `8–16` (max `34`) | Fails `(c)` (`u6_rep=41 < 50`) |
| **`Zenodo 20722064: A549_Day7.Rds`** (Drepanos et al. `10.64898/2026.07.04.736492`) | `A549` (& `K562`) CRISPRko + CRISPRi head-to-head | `2026-06-23` / `2026-07-04` (`CC-BY-4.0`) | **`25`** benchmark genes (`98` KO + `98` i + `14` neg/either) | `234,346` A549 singlets (`1.55 GB`) | **`13`** (`10`) | `3–5` (max `10`) | Fails `(c)` (`25` targets `<< 50`) |
| **`GSE284197`** (*Nature* 2026, `PMC12999477`) | Primary human cortical RG CRISPRi (`D0, D7`) | `2025-10-13` | `44` TFs (`164` sgRNAs + `20` NT) | `—` | `<= 44` | `<< 50` | Fails `(c)` (`44` genes `< 50`) |
| **`GSE319096`** | hiPSC SNaPs NPCs CRISPRko | `2026-03-01` | **`29`** (`115` guides) | `—` | **`8`** (`8`) | `< 8` | Fails `(c)` (`29` genes `< 50`) |
| **`Gupta et al. 2026`** (`10.1101/2026.08.26.747303`) | Human CRISPRko vs CRISPRi benchmark | `2026-08/10` | `29` UPR genes (`87` sgRNAs) | `—` | `<= 29` | `< 29` | Fails `(c)` (`29` genes `< 50`) + no open repo matrix |
| **`GSE200201`** (*Mol Cell* 2023, `PMC12593249`) | `MOLM-13` mSWI/SNF CRISPRko (`single` + `combo`) | `2023-04-05` | `28` BAF subunits (`+ 5` NTC) | `1.2 GB` tar | `<= 28` | `< 28` | Fails `(c)` (`28` genes `< 50`) |
| **`GSE271853`** (*Nature* 2026, `PMC12999519`) | Cortical organoids / NPCs CRISPRi | `2025-11-12` | `26` ASD/M5 genes (`96` sgRNAs + `18` NTC) | `—` | `<= 26` | `< 26` | Fails `(c)` (`26` genes `< 50`) + raw `dbGaP` controlled access |
| **`GSE305037`** (*Blood* 2026, `PMC13386094`) | `MOLM-13` AML (`Parse Evercode`, `DMSO/A-485`) | `2026-02-27` | `17` acetyl-interactome genes (`186` sgRNAs) | `1.1G` + `305M` | `<= 17` | `< 17` | Fails `(c)` (`17` genes `<< 50`) |
| **`GSE354698` / `Zenodo 20615641`** (`10.64898/2026.11.16.767105`) | `TF-iMG` microglia arrayed CRISPRko (`Parse`, `Novartis`) | `2026-11-18/19` | Arrayed: `4` genes (`Parse2`) / `9` genes (`Parse4`) | `451M–1.7G` `.rds` | `< 15` | `< 15` | Fails `(c)` (`< 15` targets `<< 50`) |

---

## 5. Exhaustive Portal-by-Portal Audit & Exclusions

### 5.1 IGVF Data Portal (`api.data.igvf.org`, `157` Released `Perturb` `AnalysisSet`s)
All released `Perturb`, `CRISPR`, and `CROP` `AnalysisSet`s (`157` + `174` + `1`) and raw `MeasurementSet` groups (`405` Perturb-seq sets) were audited:
- **Previously Screened / Excluded R9–R11 Lines**: `TeloHAEC` (`IGVFDS7943IIWZ`, `IGVFDS1612GWIZ`, etc. — failed R11 confirmatory `n_eval_targets=9`), `HCASMC` (`IGVFDS6546QFJY`, etc. — `Q3`-excluded in R11 due to `0` per-cell guide assignments in `.h5`), `K562` (`IGVFDS9078ZWQH`, etc.), `RPE-1` (`IGVFDS1890GPOK`), `WTC-11 / GM25256` (undifferentiated & differentiated EC/CM/NPC/MHC-targeted `AnalysisSet`s; high MOI `5.5` EC, `71 GB` `.pkl` CM / enhancers), `H1-hESC` (CRISPRa `IGVFDS4969LWNK` & dexamethasone multiome organoids `IGVFDS8286GRBO`).
- **Wrong Perturbation Target / Modality / Species / Size**:
  - **Enhancers (`accessible genome regions`, non-gene targets)**: `HFF-1` (`IGVFDS3405NGXF`), `H9` CM (`IGVFDS4762GRIW`, enhancers), `H9` neurons (`MOI=4`, enhancers/TFs), `WTC-11` enhancers (`IGVFDS4389OUWU`, `IGVFDS7340YDHF`).
  - **CRISPRa (not CRISPRi/ko)**: `Hs27` (`IGVFDS2001NDKP`, `Sp-dCas9-VPH` GW), `H1` (`IGVFDS4969LWNK`), `Astrocyte` (`IGVFDS0886UBIB` `90` genes `0` `AnalysisSet`s), `CD8+ T` (`IGVFDS5239SBBT`).
  - **Too Few Genes / Targeted Assay Only / Non-Perturb**: `CD8+` memory T CRISPRi (`IGVFDS9445RJOU`, **`24` genes** `< 50`), `B cell` (`5` genes; `IGVFDS2177GGWA`/`7845HWEA` unperturbed plasmablast `scRNA-seq`), `HCT116` (`IGVFDS7853IEOL`/`5308EMCG`, **`TAP-seq`** targeting `30` loci with differential TSVs only), `XM001` (`IGVFDS1028VSJQ`/`0835XTFR`, unperturbed BSA/palmitate/glucose `scRNA-seq`), `Myh11` (`IGVFDS8039PPZU`, *Mus musculus*).
- **Raw-Only `MeasurementSet`s (`0` Released Processed `AnalysisSet` Matrices)**: Verified `derived_analysis_sets == 0` for Huangfu `HUES8` (`IGVFDS0956UUQN`, `2026-11-05`, FASTQs only AND **`MOI = 12`** violating `(a)`!), Huangfu `WTC-11` (`IGVFDS5642SPLX`), Gersbach `WTC-11->hepatocyte` (`IGVFDS5284BXGZ`, `MOI = 3`, FASTQs only), and Gersbach `WTC-11` (`IGVFDS2050UJNL`). **Zero eligible ready-to-use IGVF `AnalysisSet`s remain for R12.**

### 5.2 GEO (`2016–2026/11`, `733` Human 2024–2026 Query Hits Screened)
All candidate human scCRISPR Perturb-seq / CROP-seq series were screened:
- Beyond Tier 1/2/3 entries above (`GSE345058`, `GSE327644`, `GSE218566`, `GSE337804`, `GSE281048`, `GSE291147`, `GSE263747`, `GSE273677`, `GSE283574`, `GSE280506`, `GSE284197`, `GSE319096`, `GSE200201`, `GSE271853`, `GSE305037`, `GSE354698`), remaining series fail `(a)`, `(b)`, or `(c)`:
  - `GSE326978` (`HAP1` & `MDA-MB-231`): uses **`RfxCas13d`** RNA-knockdown (`not` CRISPRi/ko) targeting lncRNAs (`0` coding `u6_rep`).
  - `GSE150062` (`WTC-11` iPSC neural induction, used in R10 + lncRNA), `GSE217812` (`SGBS` small candidate panel), `GSE159519` (`A549`, `31` genes, `u6_rep=10`), `GSE221578` (`AML`, `40` genes, `u6_rep=22`), `GSE303158` (`hPSC`, `41` genes, `u6_rep=11`), `GSE343371`/`GSE343366` (hPSC islet, `30`/`40` genes).

### 5.3 `scPerturb` (`Zenodo 13350497`, Updated `2024-08-26`, `54` Files, `CC-BY-4.0`)
All `54` `.h5ad` files were audited completely:
- **Evaluated Above**: `XuCao2023.h5ad` (`HEK293`, Tier 2 borderline `u6_rep=164`); `FrangiehIzar2021_RNA/protein.h5ad` (Melanoma `2686`, excluded `(c)` `u6_rep=51 / R9C=44`).
- **Excluded (R9–R11 Used Lines)**: `ReplogleWeissman2022` (`K562`, `RPE1`), `NadigOConner2024` (`HepG2`, `Jurkat`), `AdamsonWeissman2016` (`K562`), `DixitRegev2016` (`K562`), `GasperiniShendure2019` (`K562`), `NormanWeissman2019` (`K562` + `CRISPRa`), `XieHon2017` (`K562`), `SchraivogelSteinmetz2020` (`K562`), `DatlingerBock2017/2021` (`Jurkat`), `TianKampmann2019/2021` (`WTC11` iPSC/i3Neuron), `SunshineHein2023` (`Calu-3`), `JoungZhang2023` (`H1` + ORF overexpression).
- **Excluded (Small Targets `<50` / Wrong Modality / Non-Human)**: `PapalexiSatija2021` (`THP-1`, `25` genes), `WesselsSatija2023` (`THP-1` `Cas13d`, `28` genes), `ShifrutMarson2018` (`CD8+ T`, `20` genes); mouse (`LiangWang2023`, `SantinhaPlatt2023`, `LaraAstiasoHuntly2023`, `CuiHacohen2023`, `SchiebingerLander2019`, `WeinrebKlein2020`); small-molecule/cytokine chemical screens (`SrivatsanTrapnell2020`, `AissaBenevolenskaya2021`, `ChangYe2021`, `GehringPachter2019`, `LotfollahiTheis2023`, `McFarlandTsherniak2020`, `ZhaoSims2021`).

### 5.4 Arc Virtual Cell Atlas, Commercial Releases, Broad SCP, & Forbidden Sources
- **Arc Virtual Cell Atlas (`gs://arc-institute-virtual-cell-atlas`)**: Contains ONLY `scbasecount/` (unperturbed SRA observational scRNA-seq), `tahoe100M/` (small-molecule drugs, not CRISPR), and `virtual-cell-challenge/` (`2025 H1 hESC` `300` genes + `2026` benchmark — strictly **forbidden** under `(d)` and `H1` already used in R9). Zero eligible open CRISPR datasets.
- **Parse Biosciences / Nomic Bio `Cascade 100M` (announced `2026-10-29`)**: Cytokine/chemokine/small-molecule/drug + secretome proteomics screen (**NOT CRISPRi/ko**) and **commercial licensable** dataset (**fails `(a)` & `(d)`**).
- **Broad SCP (`SCP1064`, etc.) & `X-Atlas/Orion`**: `SCP1064` `/files/all` API returned `HTTP/2 500` (open mirror via `Zenodo 13350497` fails `(c)` anyway). `X-Atlas/Orion` excluded per `(d)`.

---

## 6. Round 12 Portfolio Recommendation & Pre-Eval (`Q3`) Action Plan

1. **Primary Confirmatory Target — `GSE345058` (`A549` CRISPRko)**:
   - Unambiguously **the only publicly available unscreened dataset** that cleanly satisfies **(a)–(e)** with massive headroom (`606,075` single-sgRNA cells, `45,320` NTCs, `1,000` genes, `u6_rep = 530`, `u6_rep ∩ R9_COMMON = 452`, `6,133 / 6,148` R9 common genes present, est. `106–212` `eval_targets` across `4` clean biological batches).
   - **Pre-download download optimization**: Prefer `GSE345058_raw_counts.mtx.gz` (`10 GB`) + `GSE345058_cell_annotations.csv.gz` (`38 MB`) + `GSE345058_gene_annotations.csv.gz` (`370 KB`) + `GSE345058_sgRNA.csv.gz` (`67 KB`) rather than the `25.05 GB` `.h5ad` (saves `13 GB` transfer & avoids double-storing normalized counts).
   - **Mandatory `Q3`/Pre-QC Checks upon Approval**: Verify CSR matrix dimensions align `606,075 x 29,892`, check whether `total_counts_CAS9 == 0` cells (`26.7%`) reduce on-target knockout/reproducibility effect or whether `530` `u6_rep` keeps high survival unmodified, and partition `batch1+3` (`23,758` NTC) vs `batch2+4` (`21,562` NTC).
2. **Optional Secondary Candidate (`K_conf` Margin Check) — `XuCao2023.h5ad` (`HEK293` CRISPRi, `336.64 MB`)**:
   - Extremely lightweight (`337 MB`) and open (`CC-BY-4.0`); worth running **`Q3` + counts-only `n_eval_targets` / `|G_eval|` gate** alongside `A549` **only if `K_conf >= 2` is required in R12**, while recognizing beforehand that it may sit near/below `50` `eval_targets` (`u6_rep = 164` / `R9C = 143` -> est. `33–66`) or fail `|G_eval| >= 3,000` due to shallow `2,155` median UMIs/cell.
3. **Do NOT Invest Download Bandwidth in `GSE327644` (`Macrophages`) Without Stronger Precedent**:
   - `u6_rep = 108` (`R9C = 92`) predicts `18–37 << 50` `eval_targets`, compounded by Flex targeted probe chemistry, single donor, and unverified per-cell MOI/singlet purity inside `6.58 GB` `.h5` tar. Keep solely as a documented contingency if `>= 55%` target survival were ever demonstrated.
4. **Portfolio Constraint Finding for Main Control**:
   - Across **all** major public portals (IGVF, GEO through `2026-11-22`, `scPerturb`, `Arc Virtual Cell Atlas`) after excluding the `13` R9–R11 cell types/lines, **`GSE345058` (`A549`) is the sole high-confidence `n_eval_targets >= 50` dataset** (and **`XuCao2023` `HEK293`** is the only plausible borderline runner-up; `HUES8` on IGVF is raw-only + `MOI=12`). If Round 12 requires `K_conf >= 2` independent confirmatory lines, either `XuCao2023` (`HEK293`) must pass `>=50` targets & `>=3000` `G_eval` in QC, OR PREREG criteria would need to evaluate `GSE327644` / revisit CRISPRa / previously tested lineage-distinct cell types.
