# Round 11 Human CRISPRi Perturb-seq Candidate Screen (`R11-0`, A-T1)

**Date**: 2026-11-22  
**Author**: A-T1 (Team A · Data Engineering)  
**Compliance with `R11_第十一轮行动指令.md` (`A-T1`)**:
- **Metadata ONLY**: Inspected API JSON/XML, GEO brief headers, published supplementary stats, and small guide/QC metadata files (`<= 13.5 MB`); **zero** expression matrices (`.h5`, `.h5ad`, `.pkl`, `.tar`, UMI count CSVs) were downloaded. (`Zenodo API HTTP 403` was encountered and NOT bypassed; `scPerturb` was cross-checked via `PMC12220817` + GEO.)
- **Hard Eligibility Filter (`未使用过的细胞` + `预计 n_eval_targets >= 50` + `<= 60 GB`)**: Only candidates on **unused human cells/donors** passing CRISPRi/knockdown modality, NT controls, public license, `<= 60 GB`, and **predicted `n_eval_targets >= 50` under the R10 empirical `20%–40%` reproducible-share rule** are admitted to the qualified candidate list (§2). Every candidate named in the R11 directive that falls short (`< 50`, `> 60 GB`, high-MOI `P4` collapse, wrong modality, or R9/R10-used lines) is explicitly documented in §3.

---

## 1. Existing R9 + R10 Baseline & Attrition Calibration Rule

- **Already-Used Cells (`K = 10`, excluded from "new cell" eligibility)**:
  - **R9 (`6` source lines)**: `K562`, `RPE1`, `HepG2`, `Jurkat`, `KOLF2.1J`, `H1`.
  - **R10 (`4` lines, already inspected)**: `DLD-1` (`GSE337988` Low-MOI -> **`86`** `eval_targets`), `Primary CD4+ T Rest` (`GSE314342` `D1=CE0008162 + D4=CE0006864` -> **`302`**), `Calu-3` (`8`), `WTC-11 i3Neuron` (`2`).
- **Reference Source Target Unions (`HGNC`-canonicalized)**:
  - **Main `PREREG_D` Source (`R9` 6 lines)**: `union6_all = 2,655`, **`union6_rep = 2,117`** (`union6_rep ∩ R9_COMMON [6,148] = 1,657`). *(A target must be in `union6_rep ∩ G_eval` and pass target-line QC/reproducibility to count toward `n_eval_targets >= 50`.)*
  - **Sensitivity `PREREG_D` Source (`R9 + dld1 + cd4t_rest`, 8 lines)**: **`union8_rep = 5,482`** (`∩ R9_COMMON = 2,996`); all 10 R9+R10 lines **`union10_rep = 5,491`**.
- **Calibration Basis (`20%–40%` rule)**:
  - In R10, raw `union6` overlap shrank to **`12.2%–12.9%`** when including low-coverage library tails (`dld1`: `~668 -> 86`; `cd4t_rest`: `~2,482 -> 302`; `calu3`: `63 -> 8`), and collapsed to **`2.9%`** (`70 -> 2`) when NTC count dropped to **`423`** (`wtc11_i3neuron`). Among well-covered shared `R9`-reproducible candidates (`~398` DE+R9-rep in `cd4t_rest` -> `302`, or `~20.8%` of `1,452` DE-shared), reproducible survival is **`~20%–40%`** of shared `union6_rep` (given sufficient single-sgRNA cells `>= 25` and strong NTC power `>> 500–1,000`). Accordingly, clearing **`n_eval_targets >= 50`** reliably requires **`>= ~150–250` shared `union6_rep` (`>= ~125–200` inside `R9_COMMON`)** **AND** strict `P4` single-sgRNA depth + adequate Half-A/Half-B NTC counts.

---

## 2. Qualifying R11 Candidates (`Unused Cells` · `CRISPRi` · `<= 60 GB` · `Est. n_eval_targets >= 50`)

| Rank | Dataset & Accession | Unused Cell / Donor / State | Modality, MOI & Library Scope | Best Download (`<=60 GB`) | Exact Single-sgRNA (`P4`) & NTC Basis | Shared Targets (`union6_rep` / `R9∩COMMON` / `u8_rep`) | Est. R11 `n_eval_targets` (`20%–40%`) & Basis | Independence / Risk Assessment |
|---|---|---|---|---|---|---|---|---|
| **1 (Primary Independent New Lineage)** | **`TeloHAEC_2285`** (`IGVFDS7943IIWZ` / `GSE210681`, `PMC10921916` *Nature* 2024, `PMID:38326615`) | **`TeloHAEC`** (human aortic endothelial, female, **resting** `treatments=None`) — **100% unused lineage & lab** | CRISPRi (`Sp-dCas9-KOX1-KRAB`), **`2,241`** canonical TSS genes (`2,285` promoters) + `1,000` controls (`400` NTC + `600` safe) | **`1.47 GB`** (`1.37 GiB`, `.h5ad` `IGVFFI8645YMGX`; IGVF `no certificate` open) | `GSE210681_Dialout` (`671,432` CBCs) + `PMC10921916`: at `UMI>=4` (`0`-guide=`352,686`), **`257,986` strict `1`-guide CBCs (`214,449` post-QC)**, `med=102`/gene pre-QC (`~85` post-QC), **`2,656` NTC + `3,912` safe** pre-QC (~2.2k NTC post-QC) across 20 10x lanes | `u6_rep=391` (`u6_all=475`, `R9∩C=299`) · **`UMI>=4 & >=25`: `u6_rep=373` (`288` in `R9∩C`)** / **`u8_rep=862`** | **`75–149`** (**`58–115`** in `R9∩C`; raw floor `12.9%×475=61`; vs R9+R10-8: **`172–345`**) | **Best overall choice**: completely independent arterial endothelial lineage, low-MOI singlets, >2k NTCs, 2.5% of storage budget. *(Verify IGVFFI8645YMGX built-in guide call threshold on prep vs raw `gpc==1` `272 u6_rep`.)* |
| **2 (Unused Donors + Unused State; Weak Line Independence)** | **`GSE314342` `Stim48hr (D2+D3)`** (`CE0010866 + CE0008678`, Zhu *Cell* 2026, `PMID:42664972`, `MIT License`) | Primary `CD4+ αβ T` (**48h TCR-stimulated**, **unused donors `D2` & `D3`**, both in batch `CD4i_R2`) | CRISPRi (`GEMX Flex`), **`12,650`** canonical GW targets + NTC (`~5–10%` controls) | **`27.39 GB`** (`48` `.h5` files: `D2 13.80 + D3 13.59 GB`; also `Stim8hr D2+D3 = 30.83 GB`, `Rest D2+D3 = 28.04 GB`) | Marson `QC_summaries`: **`3,774,210` single-sgRNA targeting (`296`/gene)** + **`165,081` single-sgRNA NTC (`87.7k` `D2` vs `77.3k` `D3`)** — +12.7% deeper than R10 `cd4t_rest` (`3.35M` tgt / `153k` NTC -> `302`) | **`u6_rep=2,081`** (**`1,649`** in `R9∩C`) / GW overlap across `u8_rep` (`~4,700+`) | **`~280–400`** (`>= 50` virtually certain, calibrated directly by R10 `cd4t_rest = 302` at lower cell count!) | **DISCLOSURE**: Same `Primary CD4+ T` cell type, sgRNA library & lab as R10 `cd4t_rest` (`D1/D4`), so **biological/lab independence is weak**; best used as an **in-lineage stimulated-state / independent-donor sensitivity arm**, not a substitute for an unrelated lineage. |
| **3 (Multi-Donor Unused iPSCs; Same State as R9 hPSC)** | **`Feng 2026 TargetedScreen`** excl. `kolf_2/3` (`PMC12903452`, `PMID:41330380` *Cell Genom* 2026, Figshare `27989294`, `MIT License`) | **`43` unused HipSci iPSC lines (`30` unused genetic donors)** (`kolf_2` & `kolf_3` excluded) | CRISPRi (`KOX1-KRAB`), **`444`** targeted genes (`TargetedScreen` only; single-sgRNA) | **`4.48 GB`** (`TargetedScreen_RNA-UMI-Counts.csv.gz` + `21.3 MB` metadata) | `TargetedScreen` alone holds **`142/142` (`100%`) NTC-containing inlets**: **`550,406` non-KOLF pert (`~1,240`/gene)** + **`7,266` non-KOLF NTCs** (`3,633`/half); `100%` of `444` targets `>= 25` | **`u6_rep=292`** (**`242`** in `R9∩C`) / **`u8_rep=337`** | Pooled excl `kolf`: **`58–117`** (**`48–97`** in `R9∩C`; depth `~1,240`/gene & `7,266` NTCs support **~`65–110`**) | **Caveats**: (a) Undifferentiated **iPSC state** already represented in R9 (`H1`, `KOLF2.1J`), so value is **cross-donor** rather than cross-tissue; (b) **`Fitness` (`32` NTC) & `NonFitness` (`10` NTC) fail same-batch NTC centering and are excluded**; (c) single-donor splits (`368–668` NTC/line) risk falling below `50` after half-splitting. |
| **4 (Borderline `>= 50` · New Vascular Lineage)** | **`HCASMC-hTERT` 971-gene** (`IGVFDS6248NDUW`, Quertermous Lab, IGVF `2025-10`) | **`HCASMC-hTERT`** (human coronary artery smooth muscle, female) — **100% unused lineage** | CRISPRi (`Sp-dCas9-ZIM3-KRAB`), Parse SPLIT-POOL, **`971`** CAD TSS genes + `250` ctrls | **`42.89 GB`** (`39.94 GiB`, `.hdf5` `IGVFFI5999NXAF`; `.h5ad` is `87.29 GB > 60`) | `hcasmc_workflow_large.pdf`: lentiviral **`MOI 0.1–0.2`** + 7d GFP FACS; `IGVFDS1435ERVG` (serum subset) = **`509,608` guide+ cells (`>500`/gene**, `~7.9%` ctrl -> `>20k` ctrls) | **`u6_rep=194`** (`u6_all=242`, **`161`** in `R9∩C`) / **`u8_rep=430`** (`308` in C) | **`39–78`** (`R9∩C`: **`32–64`**; vs `u8_rep`: **`86–172`**) — **BORDERLINE (`>= 50` requires `~26–31%` survival)** | **Borderline Warning**: Under strict `20%–40%`, the low edge (`39` / `32` in `R9∩C`) falls **below 50**, only clearing `50` if `~500+`/gene single-sgRNA depth + `ZIM3-KRAB` reach the upper quartile (`~30–40%`, midpoint **`~58`**). Preprint (`pubs=None`), consumes `42.9 GB` of budget, Parse `.hdf5` needs inspection. |

---

## 3. Audit of Named & Screened Candidates Excluded from §2 (Verified Failure Modes)

Per the updated directive ("只列预计 n_eval_targets ≥ 50 的候选"), the following datasets were audited against metadata and **excluded from §2** because they fail `est. >= 50` under `20%–40%`, exceed `60 GB`, lack strict single-sgRNA power, use already-evaluated cell lines, or use non-CRISPRi modalities:

### 3.1 IGVF `WTC-11` Endothelial & Cardiomyocyte Candidates
1. **`WTC-11` Endothelial `CC-Perturb-seq` (`IGVFDS3530HHEU` d1 `0.14 GB`, `IGVFDS7380DUFB` d2 `0.12 GB`, `IGVFDS6826YKCZ` d3 `0.23 GB`; `10.64898/2026.11.04.763740`) — EXCLUDED (`MOI=5.5` & Low-NTC Failure -> `est. 0–20 << 50`)**:
   - Direct inspection of `per-cell_quality_report` + `guide_quantifications` (`0.37–1.70 MB`) proves **`MOI = 5.5`**:
     - `IGVFDS3530HHEU` (d1 EC): `16,925` total cells -> `2,235` `guides==1` (`12,471` `>=2` guides), **`75` pure NTC (`71` single-guide NTC)** -> **`0` genes with `>= 25` single-sgRNA cells (`est. 0`)**.
     - `IGVFDS7380DUFB` (d2 EC): `21,980` total cells -> `2,978` `guides==1`, `93` pure NTC (`88` single) -> `2` genes `>= 25` (`0` in `u6_rep`).
     - `IGVFDS6826YKCZ` (d3 EC): `44,102` total cells -> `7,857` `guides==1`, `378` pure NTC (`363` single) -> `54` `u6_rep` genes `>= 25` (`est. 11–22 << 50`, plus `363` NTCs).
     - Pooled `d1+d2+d3` (`0.49 GB`): only `112` `u6_rep` genes reach `>= 25` single-sgRNA cells across mixed differentiation days with just **`522` single-guide NTCs (`546` pure)** -> `20%–40%` gives `22–45`, and Half-A/B splitting (`~261` NTCs/half) risks `wtc11_i3neuron`-like collapse (`423` NTC -> `2` targets). Also `WTC-11` was used in R10. Superseded by **`TeloHAEC_2285`**.
2. **`WTC-11` Cardiomyocytes Day 12 (`IGVFDS6332VCTO`, `PMC11785042` Sivakumar 2025) — EXCLUDED (`> 60 GB` & Borderline `u6_rep=181`)**:
   - Low MOI (`~1 sgRNA/cell`, BFP+ sorted), `IGVFFI0830FXFI` (`1.88 MB`) shows `1,863` canonical promoter/TF targets (`median = 571` cells/promoter, `1,848 >= 25`) + `30` enhancers -> **`u6_rep = 181` (`171` with `num_cell >= 25`; `R9∩C = 137`; `u8_rep = 681`)**. *(Correction vs R10_DATA_SCREEN: `331` reflected downstream `gene_names`; perturbed `intended_target_name` overlap is **`181`**.)*
   - Expected `20%–40%` on `171` = **`34–68`** (`27–55` in `R9_COMMON`, straddling `< 50`) **AND** minimum download is **`75.78 GB` (`.h5 4.64 + .pkl 71.14`)** or **`112.80 GB` (`.h5ad`)** with zero split download option on IGVF (`> 60 GB` + `> 75 GB` pickle RAM). Sister WTC-11/H9 CM optimization sets (`IGVFDS4003HZAB`, `9192YMGI`, `6924DJAZ`, `2873EPIU`) target only `~72–189` enhancers/promoters (`<< 50`); `IGVFDS1245HFTX` (`121.67–148.14 GB`) targets enhancers (`u6_rep=34` promoters).

### 3.2 All Screened iPSC-Derived Microglia & Astrocyte CRISPRi Datasets (`est. << 50` Targets)
An exhaustive GEO (`25+25` series) + EuropePMC (`500+` hits) scan confirmed **zero** public human iPSC-derived microglia or astrocyte CRISPRi scRNA-seq screens currently have enough shared gene perturbations (`u6_rep >= ~125–250`) to yield `>= 50` `eval_targets`:
- **`GSE335887`** (`PMID:42526428` *Neuron* 2026, `iTF-Microglia`): Perturbs only **`31` genes** (`est. < 15`).
- **`GSE306200`** (`PMID:42214129` 2026, `AICS-0090 WTC-11` microglia d0–d45): `184` genes in `sgRNAsequences.csv.gz` -> **`u6_rep = 16`** (`est. ~3–6`).
- **`GSE289721`** (`PMID:41340152` 2025, `iMGL` AD GWAS): `119 + 52 = 171` genes total (`73` guides / `~24` genes in GEO `feature_ref4`), far too small to yield `>= 50` R9-shared targets.
- **`PMC12339388`** (`PMID:40188316` 2025, `iTF-Microglia` ASD): `53` genes (`12` assessed due to low cell counts).
- **`GSE178317`** (`PMID:35953545`, Dräger 2022 *Nat. Neurosci.*, `iTF-Microglia`): `81 sgRNAs` (`~38` genes) CROP-seq.
- **`GSE354698`** (2026 arrayed Parse `TF-iMG`): `49`-row well map (`<= 48` arrayed wells).
- **`GSE311359`** (2025 `iMGL` MS risk loci): `9` noncoding MS GWAS loci (`5` cis genes identified).
- **`GSE182308`** (`PMID:36303069`, Leng 2022 *Nat. Neurosci.*, `hiPSC-iAstrocytes`): Only **`30` regulators** CROP-seq.
- **`GSE236057` / `GSE255009`** (`PMID:41413662` *Nat. Neurosci.* 2026, primary `NHAs` *AstroREG*): Targets **`979` candidate enhancers**, not gene TSSs.

### 3.3 High-MOI `P4`-Incompatible, Wrong-Modality & Used-Line Exclusions
- **High-MOI (`P4` collapse `< 50`)**:
  - **`H9` Cortical Neurons (`IGVFDS9825FOBP`, `PMC13629240`)**: `50.99 GB` `.h5+.pkl` (`.pkl 46.07 GB` vs `36 GiB` RAM; `.h5ad 87.43 GB > 60 GB`), `u6_rep = 207`, **enriched `MOI 2–4` (mean `3.8` sgRNAs/cell)** -> `53,037` pre-QC (`5.53%`) / `~19.2` post-QC cells/gene -> **`est. 0–15` strict `P4` targets** (already excluded in R10).
  - **`TeloHAEC_GW` (`IGVFDS6017INUS`, `44.51 GB` resting / `5371SLCS` `35.17 GB` IL1β)**: Genome-wide `19,273` genes (`u6_rep = 2,107`), **but `MOI = 5`** (`~3.4%` singlets across `19.3k` genes -> `<5` singlets/gene, **`est. ~0` strict `P4` targets**).
- **Used Cell Backgrounds (`!= 未使用细胞`)**: `GSE295214` (`H1 hESC` `7.8 GB`), `GSE343369` (`Kolf2.2J` EB `225` genes), `GSE339461` (`K562`), `Nourreddine 2026 KOLF_Strong` (`46.7 GB`). *(Note: `GSE314342 Stim48hr` and `Feng 2026 non-KOLF` also carry explicit within-lineage/within-state disclosures in §2 above.)*
- **Non-CRISPRi / Too Small (`2024–2026` scan)**: `GSE345058` (`A549` 1,000 genes **Cas9 CRISPR-KO**), `GSE283574` (`H358` **enAsCas12a KO** + `<65` genes), `GSE327644` (**CRISPR-KO** myeloid), `GSE327727` (**chemical** `D-SPIN` PBMC), `GSE309511` (`WI-38` **CRISPRa**), `GSE284197` (`~55` TFs), `GSE347104` (`SMCPilot` `22` samples), `IGVFDS3405NGXF` (`HFF-1` `87` enhancers), `GSE253842` (`SH-SY5Y` `3`-UTR loci), `GSE303901` (`DC-TAP-seq` targeted readout), `GSE279287` (`RCC` high-MOI GWAS variants), `GSE352893` (mixed `WTC11` `MOI=5` + `K562`), `GSE349847` (`OCT4/NANOG` dosage), and `scPerturb` (`PMC12220817`).

---

## 4. Concrete Round 11 Portfolio Recommendations (`K_conf >= 2`, Total `<= 60 GB`)

- **Core Mandatory Pick (`1.47 GB`)**: **`TeloHAEC_2285` (`IGVFDS7943IIWZ`)** — the only completely unused human non-T, non-hPSC line that simultaneously passes `<= 60 GB` (`1.47 GB`), strict low-MOI single-sgRNA QC (`214,449` singlets, `>2k` NTCs), peer review (*Nature* 2024), and reliable `n_eval >= 50` (**est. `75–149` vs R9**).
- **Pairing to Guarantee `K_conf >= 2` (choose based on whether `PREREG_D` prioritizes tissue novelty vs statistical headroom)**:
  - **Portfolio A (`Highest Biological Diversity`, `48.84 GB <= 60 GB`)**: `TeloHAEC_2285` (`1.47 GB`, est. `75–149`) + `Feng 2026 TargetedScreen` (`4.48 GB`, `43` non-KOLF iPSCs, est. `58–117`) + `HCASMC_971` (`42.89 GB`, est. `39–78` borderline) -> **zero Marson-CD4 / zero DLD1 overlap**, `2` somatic non-T/non-cancer vascular lines + `1` non-KOLF iPSC donor panel, with `TeloHAEC + Feng` safeguarding `K_conf >= 2` even if `HCASMC` finishes slightly under `50`.
  - **Portfolio B (`Maximum Statistical Power + 1 Brand-New Somatic Line`, `28.86 GB <= 60 GB`)**: `TeloHAEC_2285` (`1.47 GB`, est. `75–149`) + `GSE314342 Stim48hr (D2+D3)` (`27.39 GB`, est. `~280–400`) *(explicitly pre-registering `Stim48hr D2+D3` as an independent-donor TCR-activated T-cell confirmatory/sensitivity test given shared `CD4+ T` lineage & Marson lab with R10 `cd4t_rest`)*.

R11_SCREEN_DONE
