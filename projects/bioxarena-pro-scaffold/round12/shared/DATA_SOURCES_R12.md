# Round 12 Data Sources, Checksums & Download Provenance (`round12/shared/`)

**Gating Check**: `round12/status/PREREG_E_COMMIT.txt` (`4c606f4732ea080768e27757c689d40d7892cd3d`) & `round12/status/R12_GO.txt` verified present prior to download.
**Cumulative Locked Manifest**: **`11,333,944,184` B (`11.333944 GB <= 15.000000 GB`)** across **5** entries under `round12/shared/raw/` (outside `scripts/`).
**Checksum Verification File**: `round12/shared/R12_CHECKSUMS.sha256` (`shasum -a 256 -c shared/R12_CHECKSUMS.sha256` from `round12/`).

---

## 1. A549 Cas9 KO Perturb-seq (`GSE345058`, Confirmatory Primary `a549_ko` + Sensitivity `a549_ko_cas9pos`)
- **Publication / Source**: NCBI GEO `GSE345058` (A549 lung carcinoma single-gene Cas9 KO Perturb-seq, 1,000 targeted genes + 211 non-targeting sgRNAs across 32 lanes / 4 batches).
- **Files**: `4` locked files totaling **`10,997,306,975` B (`10.997307 GB`)**.

| Filename | URL | Size (Bytes) | MD5 | Post-Download SHA-256 | Integrity |
|---|---|---:|---|---|---|
| `GSE345058_cell_annotations.csv.gz` | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE345nnn/GSE345058/suppl/GSE345058_cell_annotations.csv.gz` | `40,152,708` | `68515d8b700fb282e4d01a5a2d6db7f2` | `5561c19d1020b6bcb7bf7e9a6f7a875df00097d9674723063eb0928cbe849fdf` | `gzip -t OK` |
| `GSE345058_gene_annotations.csv.gz` | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE345nnn/GSE345058/suppl/GSE345058_gene_annotations.csv.gz` | `378,455` | `f73ec036e1a8d47ead3e6df163d5e292` | `7cbed598dcda0125bcc7aad1dbb79c23ac17f296368d2d4f2873535acf54d0e9` | `gzip -t OK` |
| `GSE345058_raw_counts.mtx.gz` | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE345nnn/GSE345058/suppl/GSE345058_raw_counts.mtx.gz` | `10,956,707,618` | `2878ebec6e8aee3be7f6cfe2f9244ebb` | `abbcdd70f7de6f726c5ef6e0f63d28104f4a383fcda7502ff06d481287352189` | `gzip -t OK` |
| `GSE345058_sgRNA.csv.gz` | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE345nnn/GSE345058/suppl/GSE345058_sgRNA.csv.gz` | `68,194` | `a67bcf0d2949175da649a0b16ae4a19c` | `ec973d729a75407a091a5b363ebb8b286ef5acebfe421ede7dc042a469b4a72d` | `gzip -t OK` |

---

## 2. HEK293 dCas9-KRAB-MeCP2 (`Zenodo 13350497 / XuCao2023.h5ad`, Confirmatory Candidate 2 Subject to `E3`)
- **Publication / Source**: Xu, Cao et al. (`Zenodo 13350497`, `XuCao2023.h5ad`).
- **Files**: `1` locked file totaling **`336,637,209` B (`0.336637 GB`)**.

| Filename | URL | Size (Bytes) | Expected & Verified MD5 | Post-Download SHA-256 | Integrity |
|---|---|---:|---|---|---|
| `XuCao2023.h5ad` | `https://zenodo.org/api/records/13350497/files/XuCao2023.h5ad/content` | `336,637,209` | `3efb25a2927c5c228f04f783f18004d8` | `0065ac9524f103fcb7903125b4d2d22e8f2c29e427db31b7cefd5a8f1ea846e0` | `h5py` readable (`/X`, `/layers`, `/obs`, `/var`, `/uns`) |
