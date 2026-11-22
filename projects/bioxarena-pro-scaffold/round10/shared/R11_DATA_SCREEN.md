# Round 11 Human CRISPRi Perturb-seq Candidate Screen (`R11-0`, A-T1)

**Date**: 2026-11-22  
**Author**: A-T1 (Team A · Data Engineering)  
**Compliance**: Metadata, lightweight API JSON/XML, and small guide/QC/DE tables ONLY (`<= 13.5 MB`); **zero** expression matrices (`.h5`, `.h5ad`, `.pkl`, `.tar`, UMI count CSVs) were downloaded. All accessions, file sizes, licenses, and overlaps were API-verified in this session (with `Zenodo API HTTP 403` documented cleanly under no-bypass rules; `scPerturb` was cross-checked via `PMC12220817` + GEO instead).

---

## 1. Existing R9 + R10 Baseline & Empirical Attrition Calibration

### 1.1 Already-Used Human Cell Lines (`K = 10` total — DO NOT count as "new" cell lines)
| Round | Datasets (`lines/`) | Human Cell Line / State | Target Pool (`n_targets` / `n_rep`) | `eval_targets` vs R9-rep |
|---|---|---|---|---|
| **R9 (`K=6`)** | `replogle_k562_essential`, `replogle_rpe1`, `nadig_hepg2`, `nadig_jurkat`, `nourreddine_kolf21j`, `vcc2025_h1` | `K562`, `RPE1`, `HepG2`, `Jurkat`, `KOLF2.1J` (hiPSC/diff), `H1` (hESC) | `union6_all = 2,655` / **`union6_rep = 2,117`** (`∩ R9_COMMON = 1,657`) | Anchor reference (`K=6`) |
| **R10 (`K=4`)** | `dld1`, `cd4t_rest`, `wtc11_i3neuron`, `calu3` | `DLD-1` (colorectal, Low-MOI), `Primary CD4+ T` (Rest), `WTC-11 i3Neuron`, `Calu-3` | Combined `R9+R10` **`union10_rep = 5,491`** (`∩ R9_COMMON = 3,001`) | `dld1`: **86**, `cd4t_rest`: **302**, `calu3`: **8**, `wtc11_i3neuron`: **2** (`K_CONF = 2`) |

### 1.2 Observed R10 Empirical QC Attrition (`Addendum-B` strict single-sgRNA `P4` + dual-half NTC QC + `n_cells >= 25`)
Observed R10 survival rates calibrate all Round 11 estimates below (`shared_R9` = overlap with `union6_rep = 2,117`, and `R9_COMMON` subset `1,657`; `shared_R10` = overlap with `union10_rep = 5,491`):
- **`dld1` (`GSE337988` Low-MOI)**: `~668` raw library `union6` overlap (`93` passing canonical R9-rep pre-half) -> **`86` `eval_targets`** (**`12.9%` of raw `union6`**; **`92.5%` of high-coverage pre-half candidates**).
- **`cd4t_rest` (`GSE314342` GW)**: `~2,482` library `union6` overlap (`1,452` DE; `398` canonical R9-rep pre-half) -> **`302` `eval_targets`** (**`12.2%` of lib `union6`**, **`20.8%` of DE `union6`**, **`75.9%` of `398`**).
- **`calu3` (`SunshineHein2023`)**: `63` raw `union6` -> **`8` `eval_targets`** (**`12.7%` of raw `union6`**).
- **`wtc11_i3neuron` (`TianKampmann2021`)**: `70` raw `union6` -> **`2` `eval_targets`** (**`2.9%`** due to **`423`-NTC statistical power bottleneck** under Half-A/Half-B splitting).
- **Planning rule used below**: **`~20%–40%` of shared `R9`-reproducible (`union6_rep`) targets with adequate NTC & `>= 25` single-sgRNA cells** (or **`~12%–15%` of raw `union6_all`** when low-coverage tails are large), **with an explicit severe penalty when `NTC < 500–1,000` or `MOI > 1` strips single-sgRNA cells**.

---

## 2. Master Candidate Summary Table (Ranked by Round 11 Readiness)

| Rank / Tier | Dataset & Accession | Cell Line & State (New Lineage?) | Modality & Scale | Size (Best Format) & `<=60 GB`? | Single-sgRNA (`P4`) & NTC Power | Exact `shared_u6_rep` (`u6_all`) / `u10_rep` | Est. R11 `eval_targets` vs R9 / (vs R9+R10) | Key Verdict / Caveats |
|---|---|---|---|---|---|---|---|---|
| **1 (Tier A · Top Pick)** | **`TeloHAEC_2285`** (`IGVFDS7943IIWZ` / `GSE210681`, `PMC10921916` *Nature* 2024) | **`TeloHAEC`** aortic endothelial, **resting** (**NEW lineage**, female) | CRISPRi (`KOX1-KRAB`), **2,241** canonical TSS genes (`2,285` total promoters) + 1,000 ctrl sgRNAs | **`1.47 GB`** (`1.37 GiB`, `.h5ad` `IGVFFI8645YMGX`) — **YES (2.5%)** | **`214,449` post-QC `1`-guide singlets** (`257,986` UMI>=4 pre-QC, `med=102`/gene); **`2,656` NTC + `3,912` safe** pre-QC (~2.2k NTC post-QC) | `u6_rep=391` (`u6_all=475`, `R9∩C=299`) · **`UMI>=4 & >=25`: `u6_rep=373` (`288` in C)** / **`u10_rep=863`** | **`75–150`** (`60–115` in R9_C; ~59 floor) / **`(170–345)`** | **Best overall R11 candidate**: tiny size, peer-reviewed, clean low-MOI singlets, resting baseline, rich `>2k` NTCs |
| **2 (Tier A · New Lineage)** | **`HCASMC_971`** (`IGVFDS6248NDUW`, Quertermous Lab, IGVF `2025-10`) | **`HCASMC-hTERT`** coronary artery smooth muscle (**NEW lineage**, female) | CRISPRi (`Sp-dCas9-ZIM3-KRAB`), Parse SPLIT-POOL, **971** CAD TSS genes + 250 ctrls | **`42.89 GB`** (`.hdf5` `IGVFFI5999NXAF`; `.h5ad` `87.29 GB`) — **YES (`.hdf5`)** | **`MOI 0.1–0.2`** (`workflow_large.pdf`), `509,608` guide+ serum cells (`IGVFDS1435ERVG`, ~500+/gene, est. >20k ctrls) | **`u6_rep=194`** (`u6_all=242`, `R9∩C=161`) / **`u10_rep=430`** (`308` in C) | **`40–80`** (`32–65` in R9_C; ~30 floor) / **`(85–170)`** | Brand-new VSMC lineage + ZIM3-KRAB + low MOI; `.hdf5` `42.89 GB` fits `<=60 GB`, unreleased preprint (`pubs=None`), smaller `u6_rep` (`194`) |
| **3 (Tier B · Multi-Donor / Same hPSC State)** | **`Feng 2026 TargetedScreen`** (`PMC12903452`, `PMID:41330380` *Cell Genom* 2026, Figshare `27989294`) | **`43` non-KOLF** HipSci iPSC lines (`30` non-KOLF donors; **iPSC state**, excl `kolf_2/3`) | CRISPRi (`KOX1-KRAB`), **444** targeted genes (single-sgRNA; MIT license) | **`4.48 GB`** (`TargetedScreen` RNA `.csv.gz`) / `7.68 GB` all 3 — **YES** | **All `8,241` Targeted NTCs (`7,266` non-kolf)** in 142/142 inlets! `13` single lines each have `368–668` NTCs & `~444` targets `>=25` | Excl `kolf`: **`u6_rep=292`** (`R9∩C=242`) / **`u10_rep=338`** (Full pool: `1,072` `u6_rep` `>=25`) | Pooled Targeted: **`60–115`** (`50–95` R9_C); Per-donor: `~25–80`; (`u10`: `70–135`) | **`Fitness` (`32` NTC) & `NonFitness` (`10` NTC) lack batch NTCs**; **iPSC state already in R9 (`H1`, `KOLF2.1J`)**, useful mainly as donor-generalization control |
| **4 (Tier A2 · High-MOI / Budget >60 GB)** | **`WTC-11 CM Day12`** (`IGVFDS6332VCTO`, `PMC11785042` Sivakumar 2025) | `WTC-11` Day-12 cardiomyocyte (**new contractile lineage**, `WTC-11` donor used in R10 i3N) | CRISPRi (`Sp-dCas9`), **`1,863` canonical TF/promoter genes** + 30 enhancers | **`75.78 GB`** (`.h5 4.64 + .pkl 71.14`) or `112.80 GB` `.h5ad` — **NO (`>60 GB`)** | **Low MOI (`~1 sgRNA/cell`, BFP+ FACS)**; `IGVFFI0830FXFI` **`median=571` cells/promoter** (`1,848 >= 25`) | **`u6_rep=181`** (`171` with `num_cell>=25`; `R9∩C=137`) / **`u10_rep=682`** (`670 >= 25`) | **`35–75`** (`27–55` R9_C) / **`(135–270)`** *(requires expanded budget)* | No subpool split on IGVF; `.pkl` `71.14 GB` needs `>75 GB` RAM; earlier R10 cite of `331` counted downstream `gene_names`, true `intended_target` overlap is **`181`** |
| **5 (Tier A2 · High-MOI + Heavy RAM)** | **`H9 Cortical Neurons`** (`IGVFDS9825FOBP`, Gogate 2026 `PMC13629240`, `PMID:42827575`) | **`H9` hESC -> Day-21 Cortical Neurons** (new H9 donor; immature neuron state overlaps R10 i3N) | CRISPRi (`Sp-dCas9`), **`1,432` canonical genes** (`1,408` ASD + 600 NTC sgRNAs) | **`50.99 GB`** (`.h5 4.92 + .pkl 46.07`) or `87.43 GB` `.h5ad` — **Edge (`<=60 GB` on disk only)** | **Enriched MOI `2–4`, mean `3.8` sgRNAs/cell** (`549,727` QC cells): **strict `P4` infeasible (`~18–33` cells/gene -> most <25)** | **`u6_rep=207`** (`u6_all=244`, `R9∩C=181`) / **`u10_rep=439`** (`286` in C) | Strict `P4`: **`0–15`**; Multi-guide deconv: **`40–85`** / **`(90–175)`** | Strict `P4` fails (`~90%+` cells multi-guide); `.pkl 46.07 GB` risks host OOM unless `.h5ad` (`87.4 GB`) or streaming used |
| **6 (Tier B2 · High-MOI Exploratory)** | **`TeloHAEC_GW` Resting** (`IGVFDS6017INUS` `44.51 GB`; also IL1b `5371SLCS` `35.17 GB`) | `TeloHAEC` resting (`MOI=5` `CC-Perturb-seq`, `2026-10-12`) | CRISPRi (`KOX1-KRAB`), **19,273** genome-wide genes (`114,885` guides, `1,149` NTC) | **`44.51 GB`** (`41.45 GiB` `.h5ad`) — **YES** | **`MOI = 5`** (`~3.4%` singlets across 19.3k genes -> `<5` singlets/gene!): **strict `P4` fails**, needs deconvolution | **`u6_rep=2,107`** (`u6_all=2,639`, `R9∩C=1,652`) / **`u10_rep=5,460`** | Strict `P4`: **~0**; Multi-guide deconv: unknown pending cell-depth check | Covers `99.5%` of `union6_rep` within `<=60 GB`, **but strict single-sgRNA cannot work at MOI=5**; secondary to `TeloHAEC_2285` |
| **7 (Tier C · Low-Power EC)** | **`WTC-11 EC CC-Perturb`** (`IGVFDS3530HHEU` d1, `7380DUFB` d2, `6826YKCZ` d3; `10.64898/2026.11.04.763740`) | `WTC-11` -> Endothelial (`d1`/`d2`/`d3`, `WTC-11` donor used in R10 i3N) | CRISPRi (`KOX1-KRAB`), **298** TSS genes + 120 ctrl guides | **`0.14 GB`** (`d1`), **`0.49 GB`** (`d1+d2+d3` `.h5ad`) — **YES** | **`MOI = 5.5`!** `d1` (`3530`): **`2,235` single-guide, `75` pure NTC (`71` single), `0` `P4>=25`!** Pooled `d1–3`: `546` pure NTC (`522` single), `239` `P4>=25` | `u6_rep=151` (`u6_all=158`, `R9∩C=132`) / `u10_rep=191` (`d1–3` `P4>=25`: `112` `u6_rep`) | `3530` (d1) strict `P4`: **`0`** (deconv: `0–5`); Pooled `d1–3`: **`5–20`** (`P4`) / `15–40` (deconv) | **Subsumed by `TeloHAEC_2285`**: `MOI=5.5` + `75` (`d1`) / `546` (`d1–3`) pure NTCs mirrors `wtc11_i3neuron`'s `423`-NTC collapse (`2` targets) |

---

## 3. Detailed Technical Audit of Priority & Named Candidates

### 3.1 `TeloHAEC_2285` (`IGVFDS7943IIWZ` / `GSE210681`, Schnitzler et al. 2024 *Nature*, `PMC10921916`, `PMID:38326615`) — **#1 Recommended R11 Primary Line**
- **Biological Novelty**: Telomerase-immortalized human aortic endothelial cells (`TeloHAEC`, female) in **unstimulated resting state** (verified `treatments=None` on all IGVF `IGVFDS7943IIWZ` samples; cytokine stimulations in `PMC10921916` were separate bulk/validation assays). Completely unused endothelial lineage.
- **Exact Guide & Cell Verification (`IGVFFI9754AGFB` + `GSE210681_Dialout_guides_per_CBC.txt.gz`, `671,432` CBCs across 20 10x lanes)**:
  - Library targets `2,285` promoters (`2,241` resolved canonical protein-coding HGNC symbols + `106` non-coding/lncRNA/alt-locus IDs + `400` non-targeting + `600` safe-targeting control sgRNAs).
  - Matching `PMC10921916` (`n = 352,686` zero-guide cells): applying `UMI >= 4` per guide yields **`352,686` zero-guide CBCs**, **`257,986` strict single-guide (`1` guide) CBCs**, `60,760` multi-guide CBCs, **`2,656` NTC + `3,912` safe-targeting singlets**, and **`2,180` genes with `>= 25` single-sgRNA cells (`median = 102` cells/gene)**. Further GEX/doublet QC in `PMC10921916` retains **`214,449` single-guide cells** (`mean 929k` transcript UMIs/promoter, **`~85` cells/gene average**, **`~2,200+` NTC singlets**).
  - **Even under raw unthresholded `guides_per_cbc == 1` (`UMI >= 1`)**, `91,177` CBCs survive with `1,859` genes `>= 25` (`272` in `u6_rep`) and `1,006` NTC + `1,297` safe singlets (`IGVFFI8645YMGX` `.h5ad` should be checked on load for IGVF's binary guide layer vs dialout threshold).
- **Overlap & Expected Yield**: `373` `u6_rep` (`288` in `R9_COMMON`, `UMI>=4 >=25`), `863` `u10_rep`. Estimated **`75–150` `eval_targets` vs R9** (**`170–345` vs `R9+R10`**) in just **`1.47 GB` (`1.37 GiB`)**, open access (`no certificate`).

### 3.2 `HCASMC-hTERT` 971-gene Parse Perturb-seq (`IGVFDS6248NDUW`, Quertermous Lab) — **#2 Recommended New-Lineage Addition (`<= 60 GB`)**
- **Biological Novelty**: Human coronary artery smooth muscle cells (`HCASMC-hTERT`, female Hispanic) expressing `Sp-dCas9-ZIM3-KRAB` (stronger knockdown than KOX1-KRAB) — completely unused vascular smooth muscle lineage.
- **MOI & Design (`hcasmc_workflow_large.pdf` + `IGVFDS1435ERVG`)**: Lentiviral transduction at **`MOI 0.1–0.2`** (7-day selection, GFP+ FACS sorting, `14` Parse SPLIT-POOL libraries) across `971` CAD locus genes (`2,913` guides) + `125` NTC + `125` safe-targeting controls (`IGVFFI9546UECF`). Serum-fed subset `IGVFDS1435ERVG` alone has **`509,608` guide-positive cells** (`~500+` cells/gene).
- **Storage & Caveats**: Choose **`IGVFFI5999NXAF` `.hdf5` (`42.89 GB` = `39.94 GiB`)** rather than `.h5ad` (`87.29 GB`). Fits under `60 GB`, no `.pkl` RAM bottleneck, `u6_rep = 194` (`161` in `R9_COMMON`), `u10_rep = 430` -> estimated **`40–80` vs R9 (`85–170` vs R9+R10)**. Main caveat: consumes `71.5%` of a `60 GB` budget for half as many `u6_rep` targets as `TeloHAEC_2285`, and `Parse` `.hdf5` internal schema has not yet been inspected.

### 3.3 `WTC-11` Endothelial `CC-Perturb-seq` (`IGVFDS3530HHEU`, `IGVFDS7380DUFB`, `IGVFDS6826YKCZ`; Engreitz `2026-10-18`) — **Why It Is Demoted**
- Verified across **all three endothelial differentiation timepoints (`d1`, `d2`, `d3`)** + `d0` undifferentiated (`IGVFDS5257JGQO`, `0.31 GB`):
  - `d1` (`IGVFDS3530HHEU`, `0.14 GB`): `16,925` cells (`guides==1`: `2,235`; `guides>=2`: `12,471`), **`75` pure NTC (`71` single-guide NTC)** -> **`0` genes with `>= 25` single-sgRNA cells**.
  - `d2` (`IGVFDS7380DUFB`, `0.12 GB`): `21,980` cells (`guides==1`: `2,978`), `93` pure NTC (`88` single) -> `2` genes `>= 25` (`0` in `u6_rep`).
  - `d3` (`IGVFDS6826YKCZ`, `0.23 GB`): `44,102` cells (`guides==1`: `7,857`), `378` pure NTC (`363` single) -> `133` genes `>= 25` (`54` in `u6_rep`).
  - Pooled `d1+d2+d3` (`0.49 GB`, `83,007` cells): `546` pure NTC (`522` single), `239` single-guide `>= 25` (`112` `u6_rep`), OR `298` (`151` `u6_rep`) under multi-guide assignment. Even pooled, `522–546` NTCs (~`273`/half) risks reproducing `wtc11_i3neuron`'s `423`-NTC power failure across mixed d1–d3 cell states.

### 3.4 `Feng et al. 2026` 45 HipSci iPSCs (`PMC12903452`, *Cell Genomics* 2026, MIT License) — **Multi-Donor `TargetedScreen` Audit**
- **Donors**: `45` HipSci clones = **`31` genetic donors** (`kolf_2` `62,288 pert / 523 NTC` and `kolf_3` `49,465 pert / 458 NTC` share the R9 `KOLF2.1J` genetic background). Excluding `kolf*` leaves **`43` lines (`30` distinct non-KOLF donors)**: `734,234` perturbed singlets + `7,308` NTCs.
- **Sub-screen NTC Split (`Fitness` / `NonFitness` Disqualified for Strict Same-Batch NTC Centering)**:
  - `TargetedScreen` (`4.48 GB` RNA): **`550,406` non-kolf pert + `7,266` non-kolf NTCs across `142/142` inlets (`100%`)**, `444` targets (`100%` `>= 25`), **`u6_rep = 292`** (`242` in `R9_COMMON`), **`u10_rep = 338`**. Each of the **top `13` non-KOLF lines alone** (`oikd_5` `668` NTC, `zapk_3` `622`, `jejf_2` `605`, `iudw_1` `575`, `fiaj_1` `549`, etc.) covers all `444` targets (`288–291` in `u6_rep >= 25`) with `368–668` NTCs!
  - `FitnessGenes` (`1.45 GB`): `93,824` non-kolf pert (`937` `u6_rep >= 25`) but **only `32` non-kolf NTCs (`36` total incl `kolf`)** across `24/56` inlets.
  - `NonFitnessGenes` (`1.75 GB`): `90,004` non-kolf pert (`104` `u6_rep >= 25`) but **only `10` non-kolf NTCs (`12` total incl `kolf`)** across `9/88` inlets (`79` inlets have `0` NTCs).
- **Verdict**: Download ONLY `TargetedScreen_RNA-UMI-Counts.csv.gz` (`4.48 GB`, or `7.68 GB` if exploratory) if R11 needs an **out-of-donor iPSC (`H1`/`KOLF` already in R9) generalizability test**, expecting **`60–115` pooled `eval_targets` (`25–80` per top individual donor)**. Does NOT contribute a novel somatic differentiated tissue lineage.

### 3.5 `WTC-11` Cardiomyocytes (`IGVFDS6332VCTO`) & `H9` Cortical Neurons (`IGVFDS9825FOBP`) — **Budget & MOI Constraints**
- **`IGVFDS6332VCTO` (`WTC-11` Day-12 CMs, `PMC11785042`)**: Low MOI (`~1 sgRNA/cell`), high depth (`median = 571` cells/promoter, `p25 = 407`), **`1,863` canonical promoter/TF targets** (`u6_rep = 181`, `171` with `num_cell >= 25`; `u10_rep = 682`), plus `30` enhancers. *(Correction vs prior R10 screen: `331` reflected downstream `gene_names`; true perturbed `intended_target_name` overlap is **`181`**).* Requires **`>= 75.78 GB` (`4.64 GB .hdf5 + 71.14 GB .pkl`) or `112.80 GB` (`.h5ad`)** with zero split download option and `>75 GB` host RAM for `.pkl`; pilot subpools (`IGVFDS4003HZAB`, `9192YMGI`, `6924DJAZ`) only perturb `~72–189` enhancers/promoters.
- **`IGVFDS9825FOBP` (`H9` D21 Cortical Neurons, Gogate 2026 `PMC13629240`)**: **`1,432` canonical genes** (`u6_rep = 207`, `u10_rep = 439`), `549,727` QC cells, **mean `3.8` sgRNAs/cell** (`40th–90th` BFP FACS enriched for `MOI 2–4`). Strict `P4` yields `<~18–33` cells/gene (`0–15` eval targets). Multi-guide linear deconvolution yields **`~40–85` vs R9 (`90–175` vs R9+R10)** on **`50.99 GB` (`.h5 4.92 + .pkl 46.07`)** (requires `>= 64–96 GB` RAM to unpickle `46.07 GB`, else `87.43 GB .h5ad`).

---

## 4. Broad 2024–2026 Scan, Same-Line Expansions & Exclusions

- **Same-Line / New-State Expansions (Already-Used R9/R10 Cell Backgrounds — do NOT count as new lines)**:
  - `GSE295214` (`H1` hESC, `85` GSMs, `7.8 GB` `RAW.tar`): Genome-scale CRISPRi in **H1 hESC** (same background as R9 `vcc2025_h1`).
  - `GSE314342` `Primary CD4+ T Stim48hr` (TCR-activated state of R10 `cd4t_rest`, shares `11,903`-gene GW library; expected `~250–350` eval targets vs `u6_rep`).
  - `GSE337988` `DLD-1 Med/High MOI` (sister pools of R10 `dld1`, same `4,934` library).
  - `Nourreddine 2026 KOLF_Strong` (`46.7 GB`, same `KOLF2.1J` donor as R9).
- **Audited Exclusions (Tier C / Wrong Modality / Non-Gene Targets)**:
  - `GSE345058` (`A549`, `1,000` genes, `23 GB`): **Cas9 CRISPR-KO** (`Brunello`), NOT CRISPRi.
  - `GSE283574` (`H358`, `299 MB`, `PMID:42230636`): **`enAsCas12a` gene knockout + small `<65` gene panels**.
  - `GSE327644` (Primary myeloid VLP, `PMID:42608566`): **CRISPR-KO**, 4 GSMs.
  - `GSE327727` (`PBMC` `D-SPIN`, `PMID:42127893`): **Chemical/cytokine perturbations**, NOT CRISPRi.
  - `GSE306200` (`AICS-0090` WTC-11 iPSC-Microglia, `PMID:42214129`): CRISPRi `184` genes, only **`16` in `u6_rep`** (`40` in `u10_rep`) + WTC-11 donor.
  - `GSE289721` (iPSC-Microglia, `PMID:41340152`): `171` genes (`73` rows in GEO `feature_ref4`), small focused set.
  - `GSE347104` (`SMCPilot`, `116,062` cells, `778 MB`): only **`22` samples** (mostly CREs + few genes).
  - `GSE236057` / `GSE255009` (Primary astrocytes, `PMID:41413662`): `979` **candidate enhancers**, not gene TSSs.
  - `GSE284197` (Cortical fate `Nature` 2026, `PMID:41565813`): **`~55` TFs** (`165` guides).
  - `IGVFDS3405NGXF` (`HFF-1` fibroblasts `3.09 GB`): **`87` enhancers** only.
  - `GSE339461` (`K562`, already used in R9); `GSE253842` (`SH-SY5Y`, 3'-UTR variant screen, `PMID:42374457`); `GSE309511` (`WI-38`, `CRISPRa` 9p21 locus, `PMID:42520801`); `GSE349847` (`OCT4/NANOG` dosage).
  - `scPerturb` (`PMC12220817`, `Peidli 2024`): cross-checked GEO constituents (`Replogle K562/RPE1`, `Tian WTC-11`, `Gasperini/Schraivogel/Morris K562`, `Papalexi THP-1 KO`, `Datlinger Jurkat`); no missed `>100`-gene non-R9/R10 human CRISPRi whole-transcriptome screen.

---

## 5. Final Concrete Recommendations for Round 11 Pre-Registration

1. **Primary Low-MOI New-Lineage Selection (`<= 60 GB` budget — Strongest Choice)**:  
   **`TeloHAEC_2285` (`IGVFDS7943IIWZ` / `GSE210681`, `1.47 GB`)**: brand-new arterial endothelial lineage in resting baseline, strict low-MOI single-sgRNA (`214,449` QC singlets, `>2,200` NTCs), **`373–391` `u6_rep` (`863–896` `u10_rep`)**, predicted **`75–150` R11 `eval_targets` vs R9** (surpassing DLD-1's `86` in 2.5% of the storage budget).
2. **Secondary Low-MOI New-Lineage Addition (if total `~45 GB` disk budget fits)**:  
   **`HCASMC-hTERT` (`IGVFDS6248NDUW`, `IGVFFI5999NXAF.hdf5` `42.89 GB`)**: brand-new smooth muscle lineage with `ZIM3-KRAB`, `MOI 0.1–0.2`, `194` `u6_rep` (`430` `u10_rep`) -> predicted **`40–80` R11 `eval_targets` vs R9**, bringing R11 total to `~115–230` targets across **2 brand-new vascular lineages (`Endothelial + VSMC`)** under `44.4 GB` total!
3. **Optional iPSC Donor-Generalization Arm (`+4.48 GB`)**:  
   **`Feng 2026 TargetedScreen` (`Figshare 51072902`, `4.48 GB`)** excluding `kolf_2/3` (`43` non-KOLF lines, `7,266` NTCs, `292` `u6_rep` -> **`60–115` pooled targets**) ONLY if cross-donor robustness in an existing hPSC state is part of the R11 question (skip `Fitness`/`NonFitness` unless global non-batch NTCs are allowed).
4. **Skip under `<= 60 GB` Strict-`P4` Rule**:  
   Do **NOT** select `IGVFDS3530HHEU` (`WTC-11` EC `d1`: `0` `P4>=25` targets, `75` pure NTCs), `IGVFDS9825FOBP` (`H9` Neurons: `MOI~3.8` kills `P4` & `46.07 GB` `.pkl` stresses RAM), `IGVFDS6017INUS` (`TeloHAEC_GW`: `MOI=5` kills `P4`), or `IGVFDS6332VCTO` (`WTC-11` CMs: `75.78–112.80 GB` > 60 GB & `181` `u6_rep`), unless Round 11 explicitly expands both the disk/RAM budget (`>= 80–120 GB`) and adopts high-MOI linear deconvolution.

R11_SCREEN_DONE
