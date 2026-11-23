# Round 13 QC & Eligibility Counts (`PREREG_F` §2–§5)

**Blinding Declaration**: Strictly **counts-only** audit executed after `PREREG_F` lock (`5b51e5b`), `R13_GO.txt`, and `R13_LINES_AUTHORIZED.txt` (`dataset_readers_r13.py` `sha256: 22dfc34d9d432baa5c6b4fa7ea815b89c5f9c8ace51842ab0daec3ee609202a5`), before `CODE_FREEZE_R13.txt` / held-out evaluation. Zero transfer predictions, correlations, `r_orth_cal`, slopes, residuals, `FCT`, or `N_5` metrics were computed.
**Verification Gate**: Standardized directories under `round13/shared/lines/` verified complete (`genes.txt`, `targets.txt`, `ctrl_cells.npz`, `pseudobulk_delta.npz`, `line_metadata.json`) with zero `_tmp_*` files.

---

## 1. Pre-Metric Feasibility (`§5.1`) & Line Overview

| Line | Role (`PREREG_F` §2) | Feasibility (`§5.1`) | Channels Detected | Guide Call Source | NTC Rule | Raw Cells (Sum `.h5`) | Kept Cells (`single-sgRNA + NTC`) | Targeting Cells | Total NTC (`n_ctrl_total`) | NTC Half A / Half B | Sampled NTC | Finite Genes (`raw -> fin -> final`) |
|---|---|---|---|---|---|---:|---:|---:|---:|---:|---:|---|
| **`mac_ctrl`** | Confirmatory candidate | **`PASS`** | `L01,L02,L03,L04` (`4 .h5`) | `h5_crispr_umi_em` | `all_controls_no_safetargeting_distinction` (`593` NTC guides) | `574,620` | `122,354` | `97,675` | `24,679` | `13,354 / 11,325` | `20,000` | `18129 -> 13454 -> 13454` |
| **`mac_lps`** | Descriptive secondary only (`NEVER K_conf`) | **`PASS`** | `L01,L02,L03,L04` (`4 .h5`) | `h5_crispr_umi_em` | `all_controls_no_safetargeting_distinction` (`593` NTC guides) | `678,467` | `134,241` | `107,036` | `27,205` | `14,605 / 12,600` | `20,000` | `18129 -> 13436 -> 13436` |

---

## 2. Per-Channel Retention, NTC Counts & `< 15` NTC Fallback Audit (`§4` / `§9`)

| Line | Channel Tag | `.h5` Filename | Assigned Half | Kept Cells | Targeting Cells | Kept NTC Cells | Centering Batch Assigned | `< 15` NTC Fallback Triggered? |
|---|---|---|---|---:|---:|---:|---|---|
| **`mac_ctrl`** | `mac_ctrl_L01` | `GSM9662558_Mac_L01_Mac_Ctrl_sample_filtered_feature_bc_matrix.h5` | Half A (1) | `33,088` | `26,437` | `6,651` | `L01` (33,088) | NO |
| **`mac_ctrl`** | `mac_ctrl_L02` | `GSM9662559_Mac_L02_Mac_Ctrl_sample_filtered_feature_bc_matrix.h5` | Half B (2) | `35,621` | `28,448` | `7,173` | `L02` (35,621) | NO |
| **`mac_ctrl`** | `mac_ctrl_L03` | `GSM9662560_Mac_L03_Mac_Ctrl_sample_filtered_feature_bc_matrix.h5` | Half A (1) | `33,254` | `26,551` | `6,703` | `L03` (33,254) | NO |
| **`mac_ctrl`** | `mac_ctrl_L04` | `GSM9662561_Mac_L04_Mac_Ctrl_sample_filtered_feature_bc_matrix.h5` | Half B (2) | `20,391` | `16,239` | `4,152` | `L04` (20,391) | NO |
| **`mac_lps`** | `mac_lps_L01` | `GSM9662558_Mac_L01_Mac_LPS_sample_filtered_feature_bc_matrix.h5` | Half A (1) | `35,900` | `28,751` | `7,149` | `L01` (35,900) | NO |
| **`mac_lps`** | `mac_lps_L02` | `GSM9662559_Mac_L02_Mac_LPS_sample_filtered_feature_bc_matrix.h5` | Half B (2) | `39,334` | `31,349` | `7,985` | `L02` (39,334) | NO |
| **`mac_lps`** | `mac_lps_L03` | `GSM9662560_Mac_L03_Mac_LPS_sample_filtered_feature_bc_matrix.h5` | Half A (1) | `36,685` | `29,229` | `7,456` | `L03` (36,685) | NO |
| **`mac_lps`** | `mac_lps_L04` | `GSM9662561_Mac_L04_Mac_LPS_sample_filtered_feature_bc_matrix.h5` | Half B (2) | `22,322` | `17,707` | `4,615` | `L04` (22,322) | NO |

---

## 3. `PREREG_F` `P4`/`P5` Target Attrition, `|G_eval(t)|`, & §5 Confirmatory Gate

| Line | Deduplicated Targets (`n >= 10`) | `half_valid` (`>=2` in A & B) | **(1)** Reproducible on `t` (`half_valid`) | **(2)** Shared w/ `r9_6` (`Canon / Raw`) | **(3) `n_eval_targets`** in `G_eval` (`Canon / Raw`) | `half_valid` in `n_eval` | **`\|G_eval(t)\|`** (`Canon / Raw`) | §5 Confirmatory Gate (`>=50` targets & `>=3000` genes & `>=15` NTC/half) |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| **`mac_ctrl`** | `582` | `582` | `15` (`15`) | **`3`** / `3` | **`3`** / `3` | `3` | **`6,132`** / `5,971` | **`SECONDARY_ONLY` (`K_CONF=0`: `n_eval=3`, `|G_eval|=6132`, NTC `13354/11325`)** |
| **`mac_lps`** | `582` | `582` | `26` (`26`) | **`6`** / `6` | **`5`** / `5` | `5` | **`6,132`** / `5,971` | **`SECONDARY_ONLY`** (`PREREG_F §2` descriptive only; `n_eval=5`, `|G_eval|=6132`, NTC `14605/12600`) |

---

## 4. Adapter Audit JSON, HGNC Bridges & `CONFIRMATORY_LINES_R13.txt`

- **`mac_ctrl`** (`channel_parity_L01L03_vs_L02L04`): `{"line": "mac_ctrl", "n_h5_files": 4, "guide_call_mode": "h5_crispr_umi_em", "guide_meta_source": "template:GSE327644_PRJNA1450051_marson_mac_perturbseq_guide_template.txt.gz", "has_safe_targeting": false, "ntc_rule_used": "all_controls_no_safetargeting_distinction", "min_umi": 300, "min_ngenes": 200, "max_pct_mt": null, "n_rna_qc_excluded": 0, "n_multi_or_zero_guide_excluded": 452123, "n_safe_targeting_excluded": 0, "n_unmapped_guide_excluded": 143, "n_qc_kept_cells": 122354, "n_kept_ntc": 24679, "n_kept_half_a_ntc": 13354, "n_kept_half_b_ntc": 11325, "n_low_ntc_channel_fallback_cells": 0, "used_seed42_fallback": false, "half_rule_used": "channel_parity_L01L03_vs_L02L04", "channels_used": ["L01", "L02", "L03", "L04"], "median_guide_threshold": 5.0}` | HGNC-bridged eval symbols: `[]`
- **`mac_lps`** (`channel_parity_L01L03_vs_L02L04`): `{"line": "mac_lps", "n_h5_files": 4, "guide_call_mode": "h5_crispr_umi_em", "guide_meta_source": "template:GSE327644_PRJNA1450051_marson_mac_perturbseq_guide_template.txt.gz", "has_safe_targeting": false, "ntc_rule_used": "all_controls_no_safetargeting_distinction", "min_umi": 300, "min_ngenes": 200, "max_pct_mt": null, "n_rna_qc_excluded": 0, "n_multi_or_zero_guide_excluded": 544084, "n_safe_targeting_excluded": 0, "n_unmapped_guide_excluded": 142, "n_qc_kept_cells": 134241, "n_kept_ntc": 27205, "n_kept_half_a_ntc": 14605, "n_kept_half_b_ntc": 12600, "n_low_ntc_channel_fallback_cells": 0, "used_seed42_fallback": false, "half_rule_used": "channel_parity_L01L03_vs_L02L04", "channels_used": ["L01", "L02", "L03", "L04"], "median_guide_threshold": 5.0}` | HGNC-bridged eval symbols: `[]`
- **`CONFIRMATORY_LINES_R13.txt` (`K_CONF=0`)**: `[]`
