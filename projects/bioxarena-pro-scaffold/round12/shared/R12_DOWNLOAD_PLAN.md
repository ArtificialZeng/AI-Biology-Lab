# Round 12 Download, Verification & Pre-Unblinding Ingestion Plan (`R12_DOWNLOAD_PLAN.md`)

- **Prepared by:** `A-T1` (Team A · Data Engineering)
- **Date:** `2026-11-22`
- **Hard Execution Lock (`PREREG_E` Gate):** **ZERO expression matrices (`.mtx.gz`, `.h5ad`) have been downloaded or inspected.** Every specification below comes strictly from public GEO (`GSE345058`, `GSE218566`), Zenodo (`13350497`), EuropePMC (`PMC10961254`, `PMC12220817`), and bioRxiv (`10.64898/2026.06.01.728087`) metadata plus the previously audited `GSE345058` annotation CSVs (`cell_annotations`, `gene_annotations`, `sgRNA`). **No matrix download may begin until Master Controller commits `PREREG_E` (`PREREG_E_COMMIT`) and issues `/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/status/R12_GO.txt`.**

---

## 1. Host Compute, Storage Cap, & Compliance Rules

- **Host Environment:** macOS (`hw.memsize = 38,654,705,664 B` = **`36.0 GiB` RAM**), `/Users/mac/AI-Biology-Lab` free disk `2.4 TiB`, Python `3.13.5`, `numpy 2.2.3`, `pandas 2.2.3`, `scipy 1.15.3`, `h5py 3.12.1` (**NOTE:** `anndata` is **NOT installed** — `.h5ad` must be read directly via `h5py` CSR/CSC primitives).
- **Download Budget Cap (`Directive §3`):** **`15.000 GB`** (`15,000,000,000 B` decimal).
- **Git Hygiene Note:** `round12/` does not yet contain a `.gitignore` file; Master Controller should ensure raw/derived matrix artifacts (`round12/shared/raw/`, `round12/shared/lines/`, `*.mtx.gz`, `*.h5ad`, `*.npz`, `*.h5`) are gitignored when locking `PREREG_E`.
- **Mandatory Exclusions:** Strictly NO `GSE345058_SC_raw_normalized_counts.h5ad` (`25.05 GB`), NO `GSE345058_normalized_counts.mtx.gz` (`12 GB`), NO `X-Atlas/Orion`, NO `Arc Virtual Cell Challenge [VCC]`, NO account/terms-gated data, NO fallback substitutions if any manifest check fails.

### Cumulative Storage Budget Table (`Total = 11.33 GB <= 15.00 GB`)

| Dataset Key | Role in R12 (`Directive §0, §3`) | Cell Type & Modality | Target Directory | Files (`Count`) | Exact Bytes (HTTP `Content-Length`) | Decimal `GB` | Binary `GiB` |
|---|---|---|---|---|---:|---:|---:|
| **`a549_ko`** | **Primary Confirmatory** (`PASS >= 50 targets & >= 3,000 \|G_eval\|` required) | `A549-doxCas9` (`GSE345058`, **CRISPRko**) | `round12/shared/raw/a549_ko/` | `4` files (`raw_counts.mtx.gz` + `3` annotation CSVs) | **`10,997,306,975`** | **`11.00`** (`10.997`) | **`10.242`** |
| **`hek293_xucao`** | **Confirmatory-Gated / Else Secondary** (`EXCLUDED_FEASIBILITY_Q3` if non-integer counts or missing guide/NTC labels) | `HEK293-idCas9` (`XuCao2023.h5ad`, **CRISPRi**) | `round12/shared/raw/hek293_xucao/` | `1` file (`XuCao2023.h5ad`) | **`336,637,209`** | **`0.34`** (`0.337`) | **`0.314`** |
| **TOTAL** | **`a549_ko` + `hek293_xucao`** | — | `round12/shared/raw/` | **`5` files** | **`11,333,944,184`** | **`11.33`** (`11.334`) | **`10.556`** |

> **Headroom vs `15.000 GB` Cap:** **`3,666,055,816 B` (`3.67 GB`) remaining**.

---

## 2. Exact File Manifests & SHA-256 Acquisition Protocol

Because upstream repositories differ in published hashes, **every downloaded file MUST pass upstream validation AND generate a local SHA-256 digest** recorded immediately post-lock into `round12/shared/R12_CHECKSUMS.sha256` and `round12/shared/DATA_SOURCES_R12.md`.

### 2.1 Dataset 1 — `a549_ko` (`GSE345058`, CZ Biohub)

- **Provenance & Release:** Liu et al., *"A multimodal perturbation atlas defines the phenotypic resolution of cellular morphology"*, bioRxiv `10.64898/2026.06.01.728087` (`2026-06-02`, preprint text `cc_by_nc_nd`), GEO **`GSE345058`** (Public **`2026-08-26`**, anonymous `HTTP 200 OK`).
- **Exact Files (`4` files, `10,997,306,975` B):**
  1. `GSE345058_raw_counts.mtx.gz`
     - **URL:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE345nnn/GSE345058/suppl/GSE345058_raw_counts.mtx.gz`
     - **HTTP `Content-Length` / `Last-Modified`:** **`10,956,707,618` B** (`10.957 GB` / `10.204 GiB`), `Tue, 25 Aug 2026 22:19:51 GMT`
  2. `GSE345058_cell_annotations.csv.gz`
     - **URL:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE345nnn/GSE345058/suppl/GSE345058_cell_annotations.csv.gz`
     - **HTTP `Content-Length` / `Last-Modified`:** **`40,152,708` B** (`38.29 MiB`), `Wed, 26 Aug 2026 09:55:08 GMT`
  3. `GSE345058_gene_annotations.csv.gz`
     - **URL:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE345nnn/GSE345058/suppl/GSE345058_gene_annotations.csv.gz`
     - **HTTP `Content-Length` / `Last-Modified`:** **`378,455` B** (`369.58 KiB`), `Wed, 26 Aug 2026 09:55:08 GMT`
  4. `GSE345058_sgRNA.csv.gz`
     - **URL:** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE345nnn/GSE345058/suppl/GSE345058_sgRNA.csv.gz`
     - **HTTP `Content-Length` / `Last-Modified`:** **`68,194` B** (`66.60 KiB`), `Wed, 26 Aug 2026 17:50:35 GMT`
- **SHA-256 Acquisition & Verification Method (`a549_ko`):**
  - GEO `GSE345058/suppl/` does **not** publish pre-computed SHA-256 or MD5 files (`filelist.txt` returns `HTTP/1.1 403 Forbidden`), while individual file `HEAD` requests expose exact `Content-Length` and `Last-Modified`.
  - **Pre-unblinding download gate:**
    1. Download via `curl -4 -fSL` with resume support (`-C -`) only after `round12/status/R12_GO.txt` exists.
    2. Assert `os.path.getsize(f) == expected_content_length` for all `4` files.
    3. Verify complete stream decompression via `gzip -t` on all `4` `.gz` files.
    4. Compute local SHA-256 (`hashlib.sha256` / `shasum -a 256`) of every `.gz` file immediately after download and freeze into `round12/shared/R12_CHECKSUMS.sha256`.

### 2.2 Dataset 2 — `hek293_xucao` (`XuCao2023.h5ad`, `Zenodo 13350497`)

- **Provenance & Release:** Xu et al., *Nature Biotechnology* 2023 (`PMID 36778497`, `PMC10961254`, `10.1038/s41587-023-01948-9`, GEO `GSE218566` `2023-08-10`), harmonized `scPerturb` v1.4 release **`Zenodo 13350497`** (`DOI: 10.5281/zenodo.13350497`, updated `2024-08-26`, `CC-BY-4.0`, `HTTP 200 OK`).
- **Exact File (`1` file, `336,637,209` B):**
  1. `XuCao2023.h5ad`
     - **URL:** `https://zenodo.org/api/records/13350497/files/XuCao2023.h5ad/content`
     - **Zenodo API Exact Size:** **`336,637,209` B** (`0.337 GB` / `321.04 MiB`)
     - **Zenodo API Published MD5:** **`md5:3efb25a2927c5c228f04f783f18004d8`**
- **SHA-256 Acquisition & Verification Method (`hek293_xucao`):**
  - Zenodo REST API (`https://zenodo.org/api/records/13350497`) provides an authoritative upstream **`md5` (`3efb25a2927c5c228f04f783f18004d8`)**, not SHA-256.
  - **Pre-unblinding download gate:**
    1. Download via `curl -4 -fSL` only after `round12/status/R12_GO.txt` exists.
    2. Assert `os.path.getsize == 336,637,209` B AND computed `hashlib.md5 == "3efb25a2927c5c228f04f783f18004d8"`.
    3. Verify HDF5 readability via `h5py.File("XuCao2023.h5ad", "r")` (structural check of `/X`, `/obs`, `/var` keys only).
    4. Compute local SHA-256 (`hashlib.sha256` / `shasum -a 256`) alongside the verified Zenodo MD5 and freeze into `round12/shared/R12_CHECKSUMS.sha256`.

---

## 3. Pre-Locked NTC Controls, Single-sgRNA Filtering, & Split-Half Rules (`PREREG_E`)

All cell filtering, NTC definitions, target normalization, and split-half partitions below are **locked prior to `PREREG_E_COMMIT`** (`Directive §3`).

| Rule Dimension | **`a549_ko` (`GSE345058`)** | **`hek293_xucao` (`XuCao2023.h5ad`)** |
|---|---|---|
| **Control (NTC) Rule** | Strictly **`Condition == "Nontargeting"`** (`Perturbation == "Nontargeting"`): `211` NTC sgRNAs (**`45,320` cells** total; `0` sgRNA sequence overlap with the `4,000` targeting guides; `0` safe-targeting guides exist). | Strictly true **non-targeting control (`NTC`)** single-sgRNA cells per harmonized `.obs` (`perturbation` / `guide_id` labeled `control` / `NTC` / `non-targeting` per `scPerturb` schema; exclude any `"unassigned"`, empty, or non-singlet label). (`15` NTC sgRNAs in library.) |
| **Single-sgRNA Rule** | Join strictly by unique **`cell_id`** (`606,075` unique values; **NEVER** join by `cell_barcode` which collides across `32` lanes -> `559,858` unique). Authors already enforced: `>= 2` guide reads, multiplet removal (any cell with multiple barcodes dropped unless top > `10x` 2nd barcode), Scrublet + top-0.5% count doublet removal, `n_genes > 3000`, `UMIs > 8000`, `mt < 10%`. All `606,075/606,075` rows have `1` clean `sgRNA` & `Perturbation` (`0` delimiters). | Use **authors' single-sgRNA rule (`nperts == 1` / top sgRNA `>= 60%` & `>= 3x` 2nd sgRNA, `>= 10` sgRNA UMIs, `>= 300` UMIs, `>= 200` genes)**; **DO NOT apply authors' post-hoc "`knockdown > 40%`" filter** (`Directive §3`: retain **ALL** single-sgRNA cells to avoid selecting samples on perturbation outcome). *(Note: during `Q3` QC, log verbatim into `R12_QC_COUNTS.md` whether upstream `XuCao2023.h5ad` holds the `126,271` all-singlet pool or was upstream-exported after authors' `98,315` filter.)* |
| **Split-Half Rule** | **Biological replicate batch split (`Batch`):** **`Half A = batch1 + batch3`** (`318,447` cells: `23,758` NTC + `294,689` tgt; `16` lanes) vs **`Half B = batch2 + batch4`** (`287,628` cells: `21,562` NTC + `266,066` tgt; `16` lanes). | **Pre-locked two-stage deterministic rule:** **(Stage 1)** Inspect `.obs` replicate columns (`replicate`, `rep`, `sample`, `batch`); if `>= 2` distinct biological/experimental replicate labels exist (`e.g. rep1..rep4`), sort labels lexicographically and assign **`1`-based odd -> `Half A`**, **`even -> `Half B`** (so `rep1, rep3` vs `rep2, rep4` balances `MOI=0.1` [`1,2`] and `MOI=0.2` [`3,4`]). **(Stage 2 — Locked Fallback)** Since GEO `GSE218566` pooled replicates under single `GSM6752591`, **if `.obs` contains `< 2` distinct replicate IDs**, use **stratified `50/50` split per canonical perturbation label (NTC control stratum split `50/50` separately) with fixed `seed = 42`** (on cell-barcode-sorted indices). *(Record which stage triggered into `R12_QC_COUNTS.md` before unblinding.)* |
| **NTC Centering Unit & Adequacy (`>= 15`, `E2`)** | Center strictly by **`Lane`** (`32` lanes: `Lane1-1..Lane8-4`), all nested strictly inside `Batch` / `Half` (`E2` aligned w/ R10–R11 per-channel centering). Verified `100%` adequacy: **`32/32` lanes have `1,244–1,577` NTCs** (`>> 15`). *(In `CAS9 > 0` sensitivity arm: **`32/32` lanes have `607–1,539` CAS9>0 NTCs** `>> 15`.)* | Center by `(replicate)` / `(Half)` unit if `>= 2` replicates exist and each has `>= 15` NTCs; under Stage-2 stratified `50/50` split (`seed = 42`), center strictly within `Half A` NTCs and `Half B` NTCs respectively (`>= 15` NTCs required per half). |
| **Feature & Target Standardization** | `gene_annotations.csv.gz` (`29,892` rows): strip Ensembl `.version`, map via frozen `HGNCLookup`, and **strictly exclude synthetic non-human vector row `gene_symbol == "Cas9"` (`genome == "CAS9"`, `ensembl_id == "Cas9"`)** leaving `29,891` `GRCh38` features (`0` symbol duplicates). Resolve all `1,000` `Perturbation` targets via `HGNCLookup` (`1,000/1,000` canonical, `530` in `u6_rep`, `452` in `R9_COMMON`). | Standardize `/var` human symbols/Ensembl IDs via frozen `HGNCLookup` (excluding non-human spike-in e.g. 3T3 if present, and excluding synthetic dCas9-KRAB-MeCP2 lentiviral feature per `GSE218566` custom `hg38+lentiviral` ref). Map `perturbation` column via `HGNCLookup` (`228` genes; `164` in `u6_rep`, `143` in `R9_COMMON`). |
| **Sensitivity Arm (`Directive §2`)** | Parallel derived arm **`a549_ko_cas9pos`**: exact same pipeline restricted to cells with **`total_counts_CAS9 > 0`** in `cell_annotations.csv.gz` (`444,122` cells total = `231,892` `Half A` [`17,435` NTC] + `212,230` `Half B` [`16,250` NTC]; `1,000/1,000` perturbed targets have `>= 60` cells, median `410`). Reported descriptively (`不判定`). | None (no Cas9 transcript split in CRISPRi). |

### Exact Pre-Verified `a549_ko` Batch & Lane NTC Breakdown (`GSE345058_cell_annotations.csv.gz`)

- **`Half A` (`batch1 + batch3`, `16` lanes):** `318,447` cells (`23,758` NTC + `294,689` Perturbed):
  - `batch1` (`8` lanes `Lane1-1..Lane8-1`): `163,833` cells (`11,839` NTC, per-lane NTC `1,396–1,577`; `CAS9 > 0`: `81,785` cells, `5,824` NTC).
  - `batch3` (`8` lanes `Lane1-3..Lane8-3`): `154,614` cells (`11,919` NTC, per-lane NTC `1,405–1,570`; `CAS9 > 0`: `150,107` cells, `11,611` NTC).
- **`Half B` (`batch2 + batch4`, `16` lanes):** `287,628` cells (`21,562` NTC + `266,066` Perturbed):
  - `batch2` (`8` lanes `Lane1-2..Lane8-2`): `145,974` cells (`11,111` NTC, per-lane NTC `1,322–1,446`; `CAS9 > 0`: `114,625` cells, `8,942` NTC).
  - `batch4` (`8` lanes `Lane1-4..Lane8-4`): `141,654` cells (`10,451` NTC, per-lane NTC `1,244–1,353`; `CAS9 > 0`: `97,605` cells, `7,308` NTC).

---

## 4. Host Memory & Matrix Ingest Safeguards (For `A-T1` / `A-T2`, Post-Lock)

1. **`a549_ko` (`GSE345058_raw_counts.mtx.gz`, `10.957 GB` compressed) — Orientation & Streaming:**
   - **Programmatic Header Check (Step 0 post-lock):** Before reading triplets into arrays, inspect ONLY lines 1–3 of `gzip.open("GSE345058_raw_counts.mtx.gz", "rt")` to capture `%%MatrixMarket` header and `n_rows n_cols nnz`:
     - Verify whether `shape == (606075, 29892)` (**cell-by-gene**, matching AnnData CSVs: `"Observation (cells; rows) and var (gene; columns)"`) OR `(29892, 606075)` (**gene-by-cell**, CellRanger order); halt cleanly on any dimension mismatch.
     - Verify **row/column order alignment** against `cell_annotations.csv.gz` (`606,075` rows indexed by `cell_id`) and `gene_annotations.csv.gz` (`29,892` rows, excluding `genome == "CAS9"`).
   - **Out-of-Core Streaming within `36.0 GiB` Host RAM:** Never materialize a `606,075 x 29,891` dense array (~`72 GB` float32). Stream the `mtx.gz` triplets in bounded row/coordinate chunks (`batch_nnz ~ 25–50M` entries -> temporary HDF5 CSC/CSR or lane-sharded sparse `.npz` on `2.4 TiB` disk, keeping peak RSS `< 6–8 GiB`), reusing R10/R11 chunked QC (`compute_10x_csc_rna_qc` / `prepare_r10_lines.py`) for Pass-1 & Pass-2 HVG/AD/Energy tests.
2. **`hek293_xucao` (`XuCao2023.h5ad`, `336.64 MB`) — Pre-Eval Feasibility Exclusion Gate (`E3`):**
   - `anndata` is absent on this host, so inspect and stream `/X`, `/layers`, `/obs`, `/var` directly with **`h5py`** (supporting both CSC/CSR groups and HDF5 categorical `__categories` / `categories`+`codes` encodings).
   - **Mandatory `E3` Feasibility Exclusion Rule (Post-Download, Strictly Pre-Eval):** Because both `FCT` and `μ` require pre-normalization raw integer UMI counts, immediately verify on `XuCao2023.h5ad`:
     1. Whether `/X` or any `/layers/*` contains **non-negative integer raw UMI counts** (`integer` dtype or exact `round(x) == x` across streamed chunks);
     2. Whether `/obs` contains **per-cell `sgRNA` and true `NTC` labels**.
     - If **EITHER** raw integer UMI counts are absent from both `/X` and `/layers` **OR** per-cell `sgRNA`/`NTC` labels are absent, mark `hek293_xucao` **`EXCLUDED_FEASIBILITY_Q3`**, exclude it from **both** confirmatory adjudication AND secondary reporting, log the exact missing field diagnostic into `round12/shared/R12_QC_COUNTS.md` & `round12/status/CONFIRMATORY_LINES_R12.txt`, and compute **zero** downstream metrics on `hek293_xucao`.
3. **Strict Counts-Only Confirmatory Gate (`CONFIRMATORY_LINES_R12.txt`) Precedes `CODE_FREEZE_R12.txt`:**
   - Only after `PREREG_E_COMMIT` and `R12_GO.txt`: run counts-only QC + feasibility (`E3`) + reproducibility (`n_eval_targets >= 50` and `|G_eval(t)| >= 3,000`, `round12/shared/R12_QC_COUNTS.md`, zero evaluation correlations/FCT/sliders), emit `round12/status/CONFIRMATORY_LINES_R12.txt` (`a549_ko` alone if `hek293_xucao` fails `E3` feasibility or counts gate; both if `hek293_xucao` passes; abort confirmatory claims if `a549_ko` fails), and only then proceed to `round12` scoring.

R12_PLAN_V2_READY
