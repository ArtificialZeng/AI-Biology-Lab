# Round 11 QC & Eligibility Counts (`PREREG_D` §2.2–2.3 & Addendum A)

**Blinding Declaration**: Strictly **counts-only** audit executed prior to `CODE_FREEZE_R11.txt`. Zero transfer predictions, correlations, orthogonal calibration (`r_orth_cal`), slopes, residuals, or `FCT` values were computed on any Round 11 dataset.
**Verification Gate**: All prepared directories under `round11/shared/lines/` verified complete (`genes.txt`, `targets.txt`, `ctrl_cells.npz`, `pseudobulk_delta.npz`, `line_metadata.json`) with zero leftover temporary shards (`_tmp_*`).

---

## 1. Pre-Eval Feasibility (`Q3`), Overlap Audit (`Q6`), & Per-Line Split-Half NTC Counts

- **`hcasmc` `Q3` Pre-Eval Feasibility Gate (`IGVFFI5999NXAF.h5`)**: **`EXCLUDED_FEASIBILITY_Q3`** — RuntimeError: Parse flat-root CSR matrix (955,600 cells x 62,710 ENSG RNA features, nnz=2,679,007,357) contains ZERO CRISPR guide-capture features (guide_matches=0, non_ensg=['gene_id']); manifest guide TSVs (3163 rows) are guide library tables lacking per-cell barcodes/assignments. (shape=[955600, 62710], chunked `max_rss_gb=0.306` GB, barcodes us_suffix={'s13': 75136, 's1': 73515, 's10': 73067, 's9': 72728, 's2': 72285, 's5': 72271, 's7': 72086, 's12': 71968, 's3': 71907, 's6': 71745, 's11': 71715, 's8': 70764, 's4': 70600, 's14': 15813}).
- **`Q6` Pre-QC Disjointness Enforcement**:
  - `cd4t_stim48`: verified 0 files from `D1`/`D4`/`Rest`/`Stim8hr` in `round11/shared/raw/cd4t_stim48/`; strictly donors `CE0008678` (`D3_Stim48hr`) and `CE0010866` (`D2_Stim48hr`).
  - `feng_ipsc`: strictly excluded `kolf_2` & `kolf_3` pre-QC (`n_kolf_excluded=133,631`, 0 kept); strictly excluded all `unassigned` cells.

| Line | Status (`PREREG_D` §2.2–2.3) | Split-Half Partition (`Half A` vs `Half B`) | Kept Cells (`single-sgRNA + NTC`) | Targeting Cells | Total NTC (`n_ctrl_total`) | NTC Half A / Half B | Sampled NTC | Finite Genes (`raw -> final`) |
|---|---|---|---:|---:|---:|---:|---:|---:|
| **`telohaec`** | **Fail Confirmatory (`<50` or `<3000`)** | odd_10x_lanes_1..19_vs_even_10x_lanes_2..20 | `180,574` | `178,359` | `2,215` | `1,142` / `1,073` | `2,215` | `17472 -> 7324 -> 7324` |
| **`hcasmc`** | **Excluded (`Q3` Feasibility)** | `—` | `0` | `0` | `0` | `0` / `0` | `0` | `—` |
| **`cd4t_stim48`** | **Confirmatory (`PASS`)** | HalfA=CE0008678_D3_Stim48hr, HalfB=CE0010866_D2_Stim48hr (lexicographic) | `3,911,606` | `3,744,735` | `166,871` | `78,168` / `88,703` | `20,000` | `18130 -> 12458 -> 12458` |
| **`feng_ipsc`** | Secondary-only (`H-D5` / never `K_conf`) | 22 vs 30 batches | `557,672` | `550,406` | `7,266` | `4,091` / `3,175` | `7,266` | `36518 -> 7595 -> 7595` |

---

## 2. `PREREG_D` §2.3 Step-by-Step Target Attrition & `|G_eval(t)|`

Criterion `PREREG_D` §2.3 filters perturbations sequentially: **(1)** reproducible on target line `t` (`n_cells >= 25`, top-5000 AD `BH < 0.05 >= 50` vs `1,000` NTC, 20-PC Energy `p < 0.001` vs `500` NTC; `1` row per symbol); **(2)** reproducible on `>= 1` frozen Round 9 source line (Canon / Raw); **(3) `n_eval_targets`**: target self-gene lies in `G_eval(t) = canon(R9_COMMON_6148) ∩ genes(t)`.

| Line | Deduplicated Targets (`n >= 10`) | `half_valid` (`>=2` in A & B) | **(1)** Reproducible on `t` (`half_valid` rep) | **(2)** Shared Reproducible w/ R9 (`Canon / Raw`) | **(3) `n_eval_targets`** in `G_eval(t)` (`Canon / Raw`) | `half_valid` in `n_eval_targets` | **`\|G_eval(t)\|`** (`Canon / Raw`) | Confirmatory Gate (`>=50` targets & `>=3000` genes) |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| **`telohaec`** | `2,198` | `2,198` | `17` (`17`) | **`9`** / `8` | **`9`** / `8` | `9` | **`4,869`** / `4,741` | **`FAIL` (`n_eval=9`, `|G_eval|=4869`)** |
| **`cd4t_stim48`** | `11,730` | `10,155` | `7,301` (`6,489`) | **`655`** / `642` | **`503`** / `492` | `422` | **`6,131`** / `5,970` | **`PASS` (`503 >= 50`, `6131 >= 3000`)** |
| **`feng_ipsc`** | `444` | `444` | `135` (`135`) | **`111`** / `111` | **`93`** / `93` | `93` | **`5,326`** / `5,189` | SECONDARY_ONLY (`H-D5`) |

---

## 3. Adapter Audits, HGNC Bridges, & Confirmatory Gate Outcome (`round11/status/CONFIRMATORY_LINES_R11.txt`)

- **`telohaec`**: `{"n_h5ad_raw_cells": 214449, "n_failed_dialout_p4": 528, "n_safe_targeting_excluded": 3273, "n_tss_or_cis_excluded": 26955, "n_unmapped_excluded": 3119, "n_kept_cells": 180574, "n_kept_ntc": 2215, "lane_ntc_counts": {"1": 112, "2": 97, "3": 101, "4": 69, "5": 94, "6": 108, "7": 109, "8": 103, "9": 106, "10": 132, "11": 119, "12": 20, "13": 106, "14": 142, "15": 110, "16": 134, "17": 144, "18": 131, "19": 141, "20": 137}, "target_pooling_rule": "P4 single-sgRNA cells grouped by plain canonical HGNC symbol (15-guide/TSS library)", "half_rule": "odd_10x_lanes_1..19_vs_even_10x_lanes_2..20"}` | HGNC-bridged eval symbols: `['PHB1']`
- **`cd4t_stim48`**: `{"donors": {"CE0008678_D3": {"n_lanes": 24, "n_single_ntc": 78168, "n_single_targeting": 1796022, "median_guide_threshold": 5.0}, "CE0010866_D2": {"n_lanes": 24, "n_single_ntc": 88703, "n_single_targeting": 1948713, "median_guide_threshold": 5.0}}, "half_rule": "HalfA=CE0008678_D3_Stim48hr, HalfB=CE0010866_D2_Stim48hr (lexicographic)", "shared_lineage_lab_disclosure": "Shares CD4+ alpha-beta T cell type, library, and Marson lab with R10 cd4t_rest"}` | HGNC-bridged eval symbols: `['AFG2B', 'H2AZ1', 'HEMK2', 'IARS1', 'KRABD5', 'NUCLEOLIN', 'SKIC8', 'TAFAZZIN', 'TIMCC', 'WARS1', 'YARS1']`
- **`feng_ipsc`**: `{"n_metadata_total": 1161864, "n_kolf_excluded": 133631, "n_unassigned_excluded": 470561, "n_qc_kept_cells": 557672, "n_kept_ntc": 7266, "n_non_kolf_lines": 17, "half_a_lines": ["eipl_1", "fiaj_1", "iudw_1", "jejf_2", "oikd_2", "paab_3", "pipw_4", "tolg_4", "zapk_3"], "half_b_lines": ["eipl_3", "fiaj_3", "iudw_4", "jejf_3", "oikd_5", "paab_4", "pipw_5", "tolg_6"]}` | HGNC-bridged eval symbols: `[]`
- **Final `CONFIRMATORY_LINES_R11.txt` (`K_conf=1`)**: `['cd4t_stim48']`
- **Mandatory `Q2`/`Q7` Disclosure (`cd4t_stim48`)**: Shares primary `CD4+ alpha-beta T` cell type, sgRNA library, and Marson lab with Round 10 `cd4t_rest` (disjoint donors & `Stim48hr` state; shared R2L25..48 sequencing lanes across Half A/B per `Addendum A2`).
- **Mandatory `Addendum A3` Note (`feng_ipsc`)**: Sorted cell-line parity places sibling clones from `8` donors across Half A & Half B (clone/inlet noise ceiling, secondary-only `H-D5`).
