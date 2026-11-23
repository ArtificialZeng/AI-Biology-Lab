# Round 13 Clean-Room Reader Cross-Check Report (`READER_XCHECK_R13`)

- **Agent:** `A-T3` (Team A · Pre-Registration & Statistical Audit)
- **Scope:** Clean-room independent reconstruction of `mac_ctrl` and `mac_lps` directly from raw `.h5` files in `round13/shared/raw/gse327644/` via [`independent_r13_loader.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round13/audit/independent_r13_loader.py) (`run_r13_reader_xcheck.py`) into `round13/audit/reader_xcheck/`, compared item-by-item against `round13/shared/lines/` ([`reader_xcheck_comparison.json`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round13/audit/reader_xcheck/reader_xcheck_comparison.json)).
- **Blinding & Isolation Compliance:** Strictly zero reads/imports of `round13/src/`, strictly zero edits to `round13/shared/lines/`, and strictly zero prediction models, correlations, or evaluation metrics computed.
- **Gate Verified:** `PREREG_F_COMMIT = 5b51e5b079c2c01b44408b9ed9f8d6027b4ae4f5`, `R13_LINES_AUTHORIZED.txt`.
- **Overall Cross-Check Verdict:** **`R13_READER_XCHECK_MISMATCH`**

> **Executive Summary of Findings:** Both independent pipelines agree **bit-for-bit (`max_abs_diff = 0.0`)** on guide-calling mode (`H5_CRISPR_POISSON_GAUSSIAN_LINE_FIT` / `h5_crispr_umi_em`), **`median_guide_threshold = 5.0`**, every single NTC count per channel and per split half (`mac_ctrl`: `13,354 / 11,325`, `mac_lps`: `14,605 / 12,600`), all `ctrl_cells.npz` arrays (`X_log1p`, `X_libnorm`, `X_zscore`, `batch_labels`, and all `ctrl_mean_*`), final gene lists (`13,454` & `13,436`), and **`|G_eval| = 6,132`**. Two clean-room reader divergences produce the remaining differences:
> 1. **Per-Guide (`pert_key = guide_id`, `xcheck`) vs Gene-Pooled (`pert_key = canon_sym`, `shared`) Aggregation (`PRIMARY CAUSE`):** `shared` pooled all guides per canonical gene (`n_raw_pert_groups = 582`), whereas `xcheck` followed the Round 10/11 Flex/multi-guide standard (`cd4t_rest`, `cd4t_stim48`, `dld1`, `feng_ipsc`, `vcc2025_h1`), passing `guide_id` (`2,318` / `2,319` raw groups) and letting frozen `process_r10_line` select the top guide per gene (`best_by_sym`). This causes differences in reproducible target counts (`11 vs 15` on `mac_ctrl`; `21 vs 26` on `mac_lps`), `mac_lps` eval target set (`4` vs `5`, missing `MYC`), and `pseudobulk_delta` arrays (`delta_log1p` max diff `0.7297` / `0.8218`). Under both conventions, confirmatory candidate **`mac_ctrl` yields the exact same `3` eval targets (`['AKIRIN2', 'KAT5', 'ZC3H13']` `<< 50`)**, confirming `K_CONF = 0`.
> 2. **Unmapped `CRIPAK` Retention (`xcheck`: retained `143`/`142` cells, `583` targets) vs Exclusion (`shared`: excluded, `582` targets):** Mirroring Round 12 `a549_ko` `QARS` (`RULING_R12_HEK293_UNITS §3`), `xcheck` retained `CRIPAK` per `PREREG_F §3.4`, while `shared` filtered targeting guides to `hgnc.approved` only (`582`). This accounts for **100% of retained cell count differences** across all 8 `.h5` channels (`45+31+47+20 = 143` in `mac_ctrl`; `34+41+41+26 = 142` in `mac_lps`).

---

## 1. Line-Level Comparison Summary

### 1.1 `mac_ctrl` (`Mac_Ctrl`, Confirmatory Candidate Line -> Gate Failed, `K_CONF = 0`)

| Item | Independent `xcheck` (`A-T3`) | `round13/shared/lines/mac_ctrl` ([`R13_QC_COUNTS.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round13/shared/R13_QC_COUNTS.md)) | Match? | Exact Difference / Root Cause |
|---|---:|---:|:---:|---|
| **Channels Detected / Parity** | `4` (`L01,L03 -> Half A`; `L02,L04 -> Half B`) | `4` (`L01,L03 -> Half A`; `L02,L04 -> Half B`) | **PASS** | Exact agreement; `0` `<15` NTC fallbacks |
| **Guide Calling Mode & Median Thr** | `H5_CRISPR_POISSON_GAUSSIAN_LINE_FIT`, `median = 5.0` | `h5_crispr_umi_em`, `median = 5.0` | **PASS** | Exact agreement |
| **Raw Cells (4 `.h5`)** | `574,620` | `574,620` | **PASS** | Exact agreement |
| **Retained Cells Total** | **`122,497`** | **`122,354`** | **MISMATCH (`+143`)** | Exact difference = `CRIPAK` unmapped cells (`n_unmapped_guide_excluded = 143`) |
| **Total NTC Cells / Sampled** | **`24,679`** (`20,000` sampled) | **`24,679`** (`20,000` sampled) | **PASS** | `diff = 0` (`ctrl_batch_labels_equal = True`) |
| **Two-Half NTC Counts (`Half A / Half B`)** | **`13,354 / 11,325`** | **`13,354 / 11,325`** | **PASS** | Exact equality across Half A and Half B |
| **`ctrl_cells.npz` (`X_log1p`, `X_libnorm`, `X_zscore`, all `ctrl_mean_*`)** | `20,000 x 13,454` | `20,000 x 13,454` | **PASS (`0.0`)** | **Exact `0.0` bit-for-bit equality** across every control vector and `20,000`-cell matrix |
| **`n_genes_final` (`genes.txt`) / `|G_eval|`** | **`13,454` / `6,132`** | **`13,454` / `6,132`** | **PASS** | `genes_exact_equal = True` |
| **Raw Pert Groups -> Deduplicated Targets (`targets.txt`)** | `2,318` guides (`2,221` `>=10`) -> **`583`** (`incl. CRIPAK`) | `582` pooled genes (`582` `>=10`) -> **`582`** (`excl. CRIPAK`) | **MISMATCH** | Guide-key (`2,318`) vs gene-pooled (`582`) + `CRIPAK` (`only_in_xcheck = ['CRIPAK']`) |
| **Reproducible Targets (`targets[reproducible]`)** | **`11 / 583`** | **`15 / 582`** | **MISMATCH (`-4`)** | `11` shared; `4` only in `shared` (`AHR, IRF2, PRDM1, TRIB1`) gain power when 4 guides pool (`AHR` best guide `39` cells `ad=16` vs pooled `128` cells `ad=68`; `IRF2` `52` `ad=46` vs `183` `ad=122`) |
| **Evaluation Targets (`r9_6` `u6_rep`, rep, `G_eval`, `half_valid`)** | **`3`** (`['AKIRIN2', 'KAT5', 'ZC3H13']`) | **`3`** (`['AKIRIN2', 'KAT5', 'ZC3H13']`) | **PASS** | Exact set equality (`eval_targets_exact_equal = True`); `< 50 -> SECONDARY_ONLY (K_CONF=0)` |
| **`pseudobulk_delta` Max Abs Diff (`582 x 13,454` aligned)** | — | — | **MISMATCH (`> 1e-6`)** | `n_cells`: `205.0`; `half_a_n_cells`: `115.0`; `half_b_n_cells`: `93.0`; `ad_counts`: `313.0`; `delta_log1p`: **`0.729744`**; `half_a_delta_log1p`: `1.022624`; `half_b_delta_log1p`: `1.009811`; `knockdown_depth_log1p`: `0.540938` (caused by single best guide vs 4-guide gene pool) |

### 1.2 `mac_lps` (`Mac_LPS`, Descriptive Secondary Only)

| Item | Independent `xcheck` (`A-T3`) | `round13/shared/lines/mac_lps` ([`R13_QC_COUNTS.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round13/shared/R13_QC_COUNTS.md)) | Match? | Exact Difference / Root Cause |
|---|---:|---:|:---:|---|
| **Channels Detected / Parity** | `4` (`L01,L03 -> Half A`; `L02,L04 -> Half B`) | `4` (`L01,L03 -> Half A`; `L02,L04 -> Half B`) | **PASS** | Exact agreement; `0` `<15` NTC fallbacks |
| **Guide Calling Mode & Median Thr** | `H5_CRISPR_POISSON_GAUSSIAN_LINE_FIT`, `median = 5.0` | `h5_crispr_umi_em`, `median = 5.0` | **PASS** | Exact agreement |
| **Raw Cells (4 `.h5`)** | `678,467` | `678,467` | **PASS** | Exact agreement |
| **Retained Cells Total** | **`134,383`** | **`134,241`** | **MISMATCH (`+142`)** | Exact difference = `CRIPAK` unmapped cells (`n_unmapped_guide_excluded = 142`) |
| **Total NTC Cells / Sampled** | **`27,205`** (`20,000` sampled) | **`27,205`** (`20,000` sampled) | **PASS** | `diff = 0` (`ctrl_batch_labels_equal = True`) |
| **Two-Half NTC Counts (`Half A / Half B`)** | **`14,605 / 12,600`** | **`14,605 / 12,600`** | **PASS** | Exact equality across Half A and Half B |
| **`ctrl_cells.npz` (`X_log1p`, `X_libnorm`, `X_zscore`, all `ctrl_mean_*`)** | `20,000 x 13,436` | `20,000 x 13,436` | **PASS (`0.0`)** | **Exact `0.0` bit-for-bit equality** |
| **`n_genes_final` (`genes.txt`) / `|G_eval|`** | **`13,436` / `6,132`** | **`13,436` / `6,132`** | **PASS** | `genes_exact_equal = True` |
| **Raw Pert Groups -> Deduplicated Targets (`targets.txt`)** | `2,319` guides (`2,241` `>=10`) -> **`583`** (`incl. CRIPAK`) | `582` pooled genes (`582` `>=10`) -> **`582`** (`excl. CRIPAK`) | **MISMATCH** | Guide-key (`2,319`) vs gene-pooled (`582`) + `CRIPAK` |
| **Reproducible Targets (`targets[reproducible]`)** | **`21 / 583`** | **`26 / 582`** | **MISMATCH** | `19` shared; `7` only in `shared` (`CEBPB, DUSP1, FLI1, IRF1, MYC, NFKB1, TLR4`); **`2` only in `xcheck` (`CD14`, `TRAF3` — active single guide `ad=72`/`151` diluted to `ad=25`/`40` under 4-guide pooling)** |
| **Evaluation Targets (`r9_6` `u6_rep`, rep, `G_eval`, `half_valid`)** | **`4`** (`['AKIRIN2', 'KAT5', 'PPP1R10', 'ZC3H13']`) | **`5`** (`['AKIRIN2', 'KAT5', 'MYC', 'PPP1R10', 'ZC3H13']`) | **MISMATCH (`-1`)** | `MYC` (`xcheck` best guide `55` cells `ad=17, ep=0.0060`; `shared` pooled `178` cells `ad=65, ep=0.0005`) |
| **`pseudobulk_delta` Max Abs Diff (`582 x 13,436` aligned)** | — | — | **MISMATCH (`> 1e-6`)** | `n_cells`: `230.0`; `half_a_n_cells`: `119.0`; `half_b_n_cells`: `115.0`; `ad_counts`: `570.0`; `delta_log1p`: **`0.821834`**; `half_a_delta_log1p`: `1.074454`; `half_b_delta_log1p`: `1.106389`; `knockdown_depth_log1p`: `0.417769` |

---

## 2. Per-Channel Kept Cell & NTC Reconciliation (`L01–L04`)

| Line | Channel (Half) | `xcheck` Kept (`Targeting / NTC`) | `shared` Kept (`Targeting / NTC`) | Exact Diff (`Targeting / NTC`) | Matches `CRIPAK` Exclusion? |
|---|---|---:|---:|---:|---|
| `mac_ctrl` | `L01` (`Half A`) | `33,133` (`26,482 / 6,651`) | `33,088` (`26,437 / 6,651`) | **`+45 / 0`** | Part of `45+31+47+20 = 143` |
| `mac_ctrl` | `L02` (`Half B`) | `35,652` (`28,479 / 7,173`) | `35,621` (`28,448 / 7,173`) | **`+31 / 0`** | Exact `NTC` match (`0` diff) |
| `mac_ctrl` | `L03` (`Half A`) | `33,301` (`26,598 / 6,703`) | `33,254` (`26,551 / 6,703`) | **`+47 / 0`** | Exact `NTC` match (`0` diff) |
| `mac_ctrl` | `L04` (`Half B`) | `20,411` (`16,259 / 4,152`) | `20,391` (`16,239 / 4,152`) | **`+20 / 0`** | Total kept diff `= +143` (`CRIPAK`) |
| `mac_lps` | `L01` (`Half A`) | `35,934` (`28,785 / 7,149`) | `35,900` (`28,751 / 7,149`) | **`+34 / 0`** | Part of `34+41+41+26 = 142` |
| `mac_lps` | `L02` (`Half B`) | `39,375` (`31,390 / 7,985`) | `39,334` (`31,349 / 7,985`) | **`+41 / 0`** | Exact `NTC` match (`0` diff) |
| `mac_lps` | `L03` (`Half A`) | `36,726` (`29,270 / 7,456`) | `36,685` (`29,229 / 7,456`) | **`+41 / 0`** | Exact `NTC` match (`0` diff) |
| `mac_lps` | `L04` (`Half B`) | `22,348` (`17,733 / 4,615`) | `22,322` (`17,707 / 4,615`) | **`+26 / 0`** | Total kept diff `= +142` (`CRIPAK`) |

---

## 3. Root Causes & Governance Note for Master Controller Ruling

Both differences reflect pre-data architectural choices not uniquely pinned down in `PREREG_F.md` (neither depends on any downstream predictor or correlation):
1. **Guide Aggregation Unit (`guide_id` vs `canon_sym` pooling):** `PREREG_F §3.2` cited `prepare_r10_lines.py::build_cd4t_rest_parts`, which uses `pert_key = guide_id` (per-guide Addendum-B then `best_by_sym` selection, matching `7/7` prior guide-level readers in Rounds 9–11 including Marson Flex `cd4t_rest` and `cd4t_stim48`), whereas `round13/shared/lines/` pooled all `4` guides per canonical target symbol (`pert_key = canon_sym`, matching `a549_ko` and `hek293_xucao` in Round 12). Notably, pooling raises `n_eval_targets` on `mac_lps` (`4 -> 5`, adding `MYC`) while losing `CD14` & `TRAF3` from `reproducible` due to inactive-guide dilution, and leaves confirmatory candidate `mac_ctrl` unchanged at **`n_eval_targets = 3 << 50` (`K_CONF = 0`)**.
2. **`CRIPAK` (`hgnc.resolve` `status == 'unmapped'`):** `xcheck` retained `CRIPAK` (`143` / `142` cells; non-reproducible; absent from `r9_6` and `G_eval`) per `PREREG_F §3.4`, whereas `shared` excluded unapproved HGNC targets (`582` targets, matching `a549_ko` `QARS` in Round 12 `RULING_R12_HEK293_UNITS §3`).

No files under `round13/shared/lines/` or `round13/audit/` were altered post-comparison; reporting **`R13_READER_XCHECK_MISMATCH`** honestly awaits Master Controller ruling before downstream descriptive evaluation and freeze.
