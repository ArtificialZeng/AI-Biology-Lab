# Round 10 QC & Eligibility Counts (`PREREG_C` §2.2–2.4)

**Blinding Declaration**: Strictly **counts-only** audit executed prior to `CODE_FREEZE_R10.txt`. Zero transfer predictions, correlations, orthogonal calibration (`r_orth_cal`), slopes, `p_ds`, or `FCT` values were computed on any Round 10 dataset.
**Verification Gate**: All 4 prepared directories (`dld1`, `cd4t_rest`, `wtc11_i3neuron`, `calu3`) verified complete (`genes.txt`, `targets.txt`, `ctrl_cells.npz`, `pseudobulk_delta.npz`, `line_metadata.json`) with zero leftover temporary shards (`_tmp_*`) via `round10/shared/scripts/verify_r10_lines.py` and `test_prepare_r10.py`.

---

## 1. Per-Line Single-sgRNA QC & Split-Half Cell Counts

| Line | Status (`PREREG_C` §2.2) | Split-Half Partition (`Half A` vs `Half B`) | Kept Cells (`single-sgRNA + NTC`) | Targeting Single-sgRNA Cells | Total NTC Cells (`n_ctrl_total`) | NTC Half A / Half B | Sampled NTC (`ctrl_cells.npz`) | Finite Genes (`raw -> final`) |
|---|---|---|---:|---:|---:|---:|---:|---:|
| **`dld1`** | **Confirmatory (`PASS`)** | `12` odd (`ch_01..23`) vs `12` even (`ch_02..24`) Low-MOI channels (excl. `A1` & `NTC_MENG*`) | `701,121` | `685,008` | `16,113` | `8,032` / `8,081` | `16,113` | `36,601 -> 13,915 -> 13,901` |
| **`cd4t_rest`** | **Confirmatory (`PASS`)** | Donor `CE0006864` (`D4`, `24` lanes) vs `CE0008162` (`D1`, `23` lanes), `<=5%` MT & `>=500` genes | `3,502,875` | `3,346,822` | `156,053` | `78,174` / `77,879` | `20,000` | `18,130 -> 11,756 -> 11,756` |
| **`wtc11_i3neuron`** | Secondary-only | Lanes `lib_1, lib_3` vs `lib_2, lib_4` | `31,407` | `30,984` | `423` | `194` / `229` | `423` | `33,538 -> 13,541 -> 13,534` |
| **`calu3`** | Secondary-only | Even vs odd `gem_group` (`1,3,5,7` vs `2,4,6,8`); single-gene elements only (`2,667` dual-gene cells dropped) | `29,370` | `26,276` | `3,094` | `1,486` / `1,608` | `3,094` | `33,538 -> 13,704 -> 13,692` |
| **`h9`** (`IGVFDS9825FOBP`) | **Excluded pre-download (`P2`)** | Excluded under P2 (too few single-sgRNA cells per target; would exceed the 60 GB budget; see R10_DOWNLOAD_PLAN.md) — reason corrected by controller 12:15; zero files downloaded | `0` | `0` | `0` | `0` / `0` | `0` | `—` |

> **Notes on Ingestion Audits**:
> - `dld1`: `1,338,060` raw filtered RNA barcodes -> `1,099,062` Low-MOI hashing singlets excluding `A1` (`A2..B7`) -> `701,121` single-sgRNA + pure `NTC_501..1004` cells (`NTC_MENG1..4` excluded per author GEO note).
> - `cd4t_rest`: Pooled per-donor Poisson-Gaussian guide caller (`>= 5` UMI floor, single-sgRNA ratio `>= 3`) yielded `78,174` NTC + `1,673,859` targeting single-sgRNA on `CE0006864_D4` (vs S3 QC summary `76,636` / `1,672,455`) and `77,879` NTC + `1,672,963` targeting single-sgRNA on `CE0008162_D1` (vs S3 QC summary `76,543` / `1,675,494`), matching `QC_summaries_per_sample_lane.csv` within `~0.1%` on targeting single-sgRNA totals across all 47 lanes. Excel date-mangled CRISPR probe names (`1-Jun`, `2-Jun`) were restored cleanly to `JUN-1`, `JUN-2`.

---

## 2. `PREREG_C` §2.2 (`P5`) Step-by-Step Target Attrition & `|G_eval(t)|`

Criterion `P5` requires each evaluated perturbation on R10 target line `t` to pass three nested gates:
- **(Step 1)** Reproducible on target line `t` (`n_cells >= 25`, `ad_counts >= 50` at BH FDR `< 0.05` across top-5,000 finite genes vs `1,000` NTCs, and 20-PC Energy Distance `p < 0.001` over `2,000` perms vs `500` NTCs).
- **(Step 2)** Also reproducible on `>= 1` frozen Round 9 source line (canonicalized across R9 legacy symbols via `hgnc_complete_set.txt`; unbridged raw overlap shown alongside for comparison).
- **(Step 3 / `n_eval_targets`)** Target self-gene lies within the aligned evaluation gene universe `G_eval(t) = canon(R9_COMMON_6148) ∩ genes(t)` (confirmatory floor: `|G_eval(t)| >= 3,000` and `n_eval_targets >= 50`).

| Line | Deduplicated Targets (`n >= 10`) | `half_valid` Targets (`>=2` in A & B) | **(1)** Reproducible on `t` (`half_valid` rep) | **(2)** Shared Reproducible w/ R9 (`Canon / Raw`) | **(3) `n_eval_targets`** in `G_eval(t)` (`Canon / Raw`) | `half_valid` subset of `n_eval_targets` (`Canon`) | **`\|G_eval(t)\|`** (`Canon / Raw`) | Confirmatory Gate (`>=50` targets & `>=3000` genes) |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| **`dld1`** | `4,934` | `4,934` | `202` (`202`) | **`93`** / `90` | **`86`** / `83` | `86` | **`6,124`** / `5,962` | **`PASS` (`86 >= 50`, `6124 >= 3000`)** |
| **`cd4t_rest`** | `11,903` | `11,531` | `3,675` (`3,604`) | **`398`** / `396` | **`302`** / `301` | `284` | **`6,123`** / `5,962` | **`PASS` (`302 >= 50`, `6123 >= 3000`)** |
| **`wtc11_i3neuron`** | `184` | `184` | `4` (`4`) | **`2`** / `2` | **`2`** / `2` | `2` | **`5,995`** / `5,835` | `SECONDARY_ONLY` |
| **`calu3`** | `178` | `178` | `22` (`22`) | **`9`** / `9` | **`8`** / `8` | `8` | **`6,134`** / `5,972` | `SECONDARY_ONLY` |

---

## 3. Raw-Symbol vs HGNC-Bridged Differences (`dld1`: `83 -> 86`, `cd4t_rest`: `301 -> 302`)

Using frozen `hgnc_complete_set.txt` (`approved -> ENSG -> unambiguous prev_symbol -> unambiguous alias_symbol`) bridges outdated HGNC symbols present in Round 9 lines (`R9 COMMON_GENES.txt` contains legacy `prev_symbol`s) to modern approved HGNC symbols used in Round 10 annotations (`GRCh38-2020-A` / Flex v1 / Cell 2026):
1. **`dld1` (`83` raw string -> `86` canonical HGNC, `+3` targets)**:
   - Unbridged raw overlap is **`83`** (`90` before restricting to `G_eval(dld1)`), matching independent raw-string verification (`~83`).
   - Canonical HGNC bridging recovers **3 exact R9 `prev_symbol` targets** whose modern HGNC symbols are both reproducible in DLD-1 and expressed inside `G_eval(dld1)`:
     - `DARS` (`nadig_hepg2`, `nadig_jurkat`, `replogle_k562_essential`, `replogle_rpe1`) → **`DARS1`**
     - `NARS` (`nadig_hepg2`, `nadig_jurkat`, `replogle_k562_essential`, `replogle_rpe1`) → **`NARS1`**
     - `TAZ` (`nadig_jurkat`, `vcc2025_h1`) → **`TAFAZZIN`**
   - For readout genes, `|G_eval(dld1)|` expands from `5,962` raw strings to **`6,124`** canonical HGNC genes (`+162` bridged legacy R9 `COMMON_GENES` symbols).
2. **`cd4t_rest` (`301` raw string -> `302` canonical HGNC, `+1` target)**:
   - Unbridged raw overlap is **`301`** (`396` before `G_eval(cd4t_rest)`), and canonical bridging recovers **1 exact R9 `prev_symbol` target**:
     - `H2AFZ` (`nadig_hepg2`, `replogle_rpe1`) → **`H2AZ1`**
   - For readout genes, `|G_eval(cd4t_rest)|` expands from `5,962` raw strings to **`6,123`** canonical HGNC genes (`+161` bridged symbols).
   - Note that across `cd4t_rest`'s `302` evaluated targets, **`284`** satisfy `half_valid == True` (`>= 2` cells in both `CE0006864_D4` and `CE0008162_D1`), while `18` pass full-cohort reproducibility with `< 2` cells in one of the two donors (flagged cleanly via `half_valid=False` for `r_orth_cal` noise-ceiling discipline per `PREREG_C` §2.3).

---

## 4. Confirmatory Gate Outcome (`round10/status/CONFIRMATORY_LINES.txt`)

- **`CONFIRMATORY_LINES`**: `dld1`, `cd4t_rest` (**`K_conf = 2 >= 2`**)
- **`SECONDARY_LINES`**: `wtc11_i3neuron`, `calu3`
- **`EXCLUDED_PRE_DOWNLOAD`**: `h9` (`IGVFDS9825FOBP`, excluded under `PREREG_C` §2.2 `P2`)
