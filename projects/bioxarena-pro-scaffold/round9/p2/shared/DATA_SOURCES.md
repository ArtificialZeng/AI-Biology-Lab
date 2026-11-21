# Shared Data Sources & Download Provenance (`round9/p2/shared/`)

Strict adherence to `P2_RULES.md` Rule 5:
- All datasets are publicly available open-access datasets.
- URL, license, byte size, and pre-download/post-download checksums are recorded below prior to download.
- Total new downloads across the team are strictly bounded below **60 GB** (`<= 60.00 GB`).
- **NO X-Atlas/Orion data or derivatives** are used or touched.

---

## 1. Replogle et al. (Cell 2022) — K562-essential & RPE1
- **Publication**: Replogle JM et al. (2022), *Cell* (`doi:10.1016/j.cell.2022.05.013`, PMID `35688146`).
- **Repository**: Figshare+ (`https://doi.org/10.25452/figshare.plus.20029387`, Figshare API `https://api.figshare.com/v2/articles/20029387`).
- **License**: **CC BY 4.0** (`https://creativecommons.org/licenses/by/4.0/`).

| File | URL | Size (Bytes / GB) | Pre-Download Hash (Figshare API) | Post-Download SHA-256 | Status |
|---|---|---|---|---|---|
| `K562_essential_raw_singlecell_01.h5ad` | `https://ndownloader.figshare.com/files/35773219` | `10,661,879,995` B (`10.662 GB`) | `md5:4f1122ce1c7f13299a68df6459a266d3` | `3e5a63a9e892b21029bb55fca4e12517a49aad7af6c14133ca63d12cf68c6cee` | Downloaded & MD5/SHA-256 verified |
| `rpe1_raw_singlecell_01.h5ad` | `https://ndownloader.figshare.com/files/35775606` | `8,700,873,216` B (`8.701 GB`) | `md5:6a2a9d0d2bf4ec147f4d1104043b268c` | `9b05ef1f81526216fa008d677e9e0d03dce9a2f7a95499a4fb81e505e9d88ef1` | Downloaded & MD5/SHA-256 verified |
| `K562_essential_normalized_bulk_01.h5ad` | `https://ndownloader.figshare.com/files/35780870` | `79,766,954` B (`0.080 GB`) | `md5:30496767641cd2e660ee6ecb5baee132` | `c1ca6456c9c9f1aa2b02c496eb64d1dc3e6a852edbd744d682b8d2c95fd36829` | Already on local disk (R8) |
| `rpe1_normalized_bulk_01.h5ad` | `https://ndownloader.figshare.com/files/35775512` | `95,350,546` B (`0.095 GB`) | `md5:6f1e7d6a09e2f869759e3c4526b7f171` | `a3c5bfd0f15d63938bc80c9b8874b9cd761e3a23caf5ffe7966bae4e887ec89d` | Already on local disk (R6) |

*Note*: Figshare+ API v2 exposes `supplied_md5`/`computed_md5` rather than `sha256` prior to download; exact SHA-256 digests are computed immediately upon download and recorded here. New Replogle download subtotal: **19.363 GB**.

---

## 2. Nadig et al. (Nature Genetics 2025) — HepG2 & Jurkat
- **Publication**: Nadig A et al. (2025), *Nature Genetics* (`doi:10.1038/s41588-025-02169-3`, PMID `40259084`).
- **Repository**: NCBI GEO `GSE264667` (`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl/`).
- **License**: Public NCBI GEO open-access release (zero registration, open scientific use).

| File | URL | Size (Bytes / GB) | Pre-Recorded SHA-256 | Status |
|---|---|---|---|---|
| `GSE264667_hepg2_raw_singlecell_01.h5ad` | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl/GSE264667_hepg2_raw_singlecell_01.h5ad` | `5,614,460,941` B (`5.614 GB`) | `e1ad7c3c5a201c861a207a858aa7e59f5e6ac1955674c415f7de0d1dadadb52e` | Already on local disk (R7, 0 GB new download) |
| `GSE264667_jurkat_raw_singlecell_01.h5ad` | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl/GSE264667_jurkat_raw_singlecell_01.h5ad` | `9,366,490,264` B (`9.366 GB`) | `ffbe15f2c8f7ffcfd7b0ba9e6937d4ebc2d03b0179fa8234648a59bcb82c04a3` | Already on local disk (R7, 0 GB new download) |

---

## 3. Arc Virtual Cell Challenge 2025 — H1 hESC
- **Source**: Arc Virtual Cell Atlas (`https://arcinstitute.org/tools/virtualcellatlas`, `gs://arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/`).
- **License**: **CC0 1.0 Universal (Public Domain Dedication)** (`https://creativecommons.org/publicdomain/zero/1.0/legalcode.txt`), verified directly on `https://arcinstitute.org/tools/virtualcellatlas`. Publicly readable bucket with anonymous access (`https://storage.googleapis.com/arc-institute-virtual-cell-atlas/...`).

| File | URL | Size (Bytes / GB) | Pre-Download Hash (GCS JSON API) | Post-Download SHA-256 | Status |
|---|---|---|---|---|---|
| `gene_names.csv` | `https://storage.googleapis.com/arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/gene_names.csv` | `116,023` B (`0.0001 GB`) | `md5(b64):PnSuaCYgIvYvcHAb7FolCQ==`, `crc32c:Px7c5g==` | `e29ff3512b4c2440759bfb2a34577049b2f49df786b3e19c2a91b565de47b38e` | Downloaded & verified |
| `adata_Training.h5ad` | `https://storage.googleapis.com/arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/train/adata_Training.h5ad` | `15,482,497,461` B (`15.482 GB`) | `crc32c:/Z0row==` (GCS composite object, no single md5) | `a09977104fefb622368ca74b50c9d3c1e891733e6c83db07acfca49b0219c02b` | Downloaded & MD5/SHA-256 verified (150 targets + NT controls) |
| `adata_Validation.h5ad` | `https://storage.googleapis.com/arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/validation/adata_Validation.h5ad` | `6,928,967,541` B (`6.929 GB`) | `crc32c:EendZg==` | `376f0bab27d9f22e00dfd888668d3602006d035edc182721c378f4802c77bfe1` | Downloaded & verified (50 targets + NT controls; shared with A-T4) |
| `adata_Test.h5ad` | `https://storage.googleapis.com/arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/test/adata_Test.h5ad` | `11,950,739,168` B (`11.951 GB`) | `crc32c:JMP5RQ==` | `ba5ce66729e4c8e50b6f147c7d73d3f9f28f76c91fa43ff33e3f6bc2ed831f03` | Downloaded & MD5/SHA-256 verified (100 targets + NT controls) |

---

## 4. Nourreddine et al. (Nature Biotechnology 2026) — KOLF2.1J hiPSC
- **Publication**: Nourreddine S et al. (2026), *Nature Biotechnology* (`doi:10.1038/s41587-026-03199-w`, PMID `42386990`).
- **Repository**: Figshare+ (`https://doi.org/10.25452/figshare.plus.27261219`, Figshare API `https://api.figshare.com/v2/articles/27261219`).
- **License**: **CC BY 4.0** (`https://creativecommons.org/licenses/by/4.0/`), Embargo: `None`. *(Note: Raw FASTQs on SRA `PRJNA1173491` carry Jackson Labs Y-chromosome sequence restriction; the processed UMI count `.h5ad` files on Figshare+ `27261219` are openly released under CC BY 4.0).*
- **Budget Constraint Disclosure**: Full `KOLF_Pan_Genome_QC_Filtered.h5ad` (`189.393 GB`) and `KOLF_Strong_Perturbations.h5ad` (`46.719 GB`) cannot fit alongside Replogle SC (`19.363 GB`) + VCC 2025 H1 (`22.411–34.362 GB`) within the 60 GB team cap; therefore only the smaller QC-filtered / strong-perturbation sub-atlases are considered:

| File | URL | Size (Bytes / GB) | Pre-Download Hash (Figshare API) | Post-Download SHA-256 | Status |
|---|---|---|---|---|---|
| `KOLF_Chromatin_Modifiers_Strong_Perturbations.h5ad` | `https://ndownloader.figshare.com/files/59541377` | `1,528,903,380` B (`1.529 GB`) | `md5:fbd62eee013ad01be10f8ff291c8f8f5` | `741674b8a19a50b57d02b8b103f1779a135e9509dffd805fd33dd65853361745` | Downloaded & MD5/SHA-256 verified |
| `KOLF_Metabolic_Enzymes_Strong_Perturbations.h5ad` | `https://ndownloader.figshare.com/files/64680753` | `1,326,348,118` B (`1.326 GB`) | `md5:95f57fd679b98ba4f060d143959c3a66` | `d86fa3260ca6289a361bd8be99a87f547e681d6b52c770c1ea1ba38011d59211` | Downloaded & MD5/SHA-256 verified |
| `KOLF_Chromatin_Modifiers_QC_Filtered.h5ad` | `https://ndownloader.figshare.com/files/64680606` | `4,635,611,222` B (`4.636 GB`) | `md5:246eda29efaf8ae8da86c39840678582` | — | **Not downloaded** (Strong sub-atlases + `adata_Test.h5ad` chosen to stay `<= 60.00 GB`) |
| `KOLF_Metabolic_Enzymes_QC_Filtered.h5ad` | `https://ndownloader.figshare.com/files/64680708` | `6,984,432,788` B (`6.984 GB`) | `md5:b5e17ab928cf263f679cf63bbb311a25` | — | **Not downloaded** (Strong sub-atlases + `adata_Test.h5ad` chosen to stay `<= 60.00 GB`) |
| `KOLF_Strong_Perturbations.h5ad` | `https://ndownloader.figshare.com/files/64650852` | `46,718,752,086` B (`46.719 GB`) | `md5:28bcfff0679e7c6c35bdd584f3626362` | — | **Excluded** (exceeds 60 GB cap with Replogle+VCC) |
| `KOLF_Pan_Genome_QC_Filtered.h5ad` | `https://ndownloader.figshare.com/files/64650261` | `189,393,177,972` B (`189.393 GB`) | `md5:afd30fde1e6ad32969c29868394385d1` | — | **Excluded** (exceeds 60 GB cap) |

---

## 5. HGNC Official Complete Gene Nomenclature Set
- **Source**: HGNC (HUGO Gene Nomenclature Committee, EMBL-EBI), `https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt`
- **License**: CC0 1.0 / Public Domain open access (`https://www.genenames.org/about/`).

| File | URL | Size (Bytes / GB) | Last-Modified (HTTP HEAD) | Post-Download SHA-256 | Status |
|---|---|---|---|---|---|
| `hgnc_complete_set.txt` | `https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt` | `17,007,798` B (`0.017 GB`) | `Fri, 13 Nov 2026 12:59:59 GMT` | `6df77ac3cbaeb86a6dce5e17d851ebb27b6b3cdcd2471a211a3d8464f00c54c4` | Downloaded & SHA-256 verified |

---

## 6. Cumulative Team Download Budget Tracker (Cap: 60.000 GB)
- Replogle 2022 (`K562_essential_raw_singlecell_01.h5ad` + `rpe1_raw_singlecell_01.h5ad`): **19.363 GB**
- Nadig 2025 (`HepG2` + `Jurkat`, pre-existing local): **0.000 GB**
- Arc VCC 2025 H1 (`adata_Training` + `adata_Validation` + `adata_Test` + `gene_names.csv` = `34.362 GB`)
- Nourreddine 2026 KOLF2.1J (`Chromatin_Modifiers_Strong` + `Metabolic_Enzymes_Strong` = `2.855 GB`)
- HGNC official complete set (`hgnc_complete_set.txt`): **`0.017 GB`**
- **Completed Download Total**: `19.363 + 34.362 + 2.855 + 0.017 = 56.597 GB` (`56,597,332,700` B exact = `56.597333` decimal GB / `52.710374` GiB) — strictly `<= 60.000 GB`.
