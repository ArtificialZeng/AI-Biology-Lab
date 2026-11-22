# Round 12 Clean-Room Reader Cross-Check Report (`READER_XCHECK_R12`)

**Agent:** `A-T3` (Team A · Pre-Registration & Statistical Audit)  
**Scope:** Independent raw ingestion & standardization from `round12/shared/raw/` via `independent_r12_loader.py::standardize_raw_cohort_to_audit_line` into [`round12/audit/reader_xcheck/`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/audit/reader_xcheck), compared item-by-item against `round12/shared/lines/` (`round12/results/` and `round12/shared/lines/` strictly untouched; zero FCT, correlations, or evaluation metrics computed).  
**Status Summary:** **`R12_READER_XCHECK_MISMATCH`** (both `hek293_xucao` and `a549_ko` completed; `hek293_xucao` shows critical batch/split discrepancies collapsing `n_eval_targets` from `62` to `19`; `a549_ko` matches `100%` on `445` eval targets / `6,121` `G_eval` / two-half NTCs / `ctrl_mean` `0.0` and differs only on `QARS` `591` cells + `3.58e-7` float32 ULP in `delta_log1p`).

---

## 1. `hek293_xucao` (`XuCao2023.h5ad`) Cross-Check — **MISMATCH**

### 1.1 Pre-Fix Clean-Room Finding (`sha256: 4fbd6da...`)
An initial cross-check run (`task-3108`, `119.2s`) using the pre-check `independent_r12_loader.py` exposed a layer-selection bug in `check_hek293_xucao_feasibility`: `XuCao2023.h5ad` stores full cellular counts in `/X` (`sparse_int_int64`, sum `250,678,179` vs `obs.ncounts` sum `250,838,663`) AND a 4sU newly-synthesized subfraction in `/layers/nascent_counts` (`sparse_int_int64`, sum `45,564,947` vs `obs.nascent_UMI_counts` sum `45,571,050`). Because substring `"count"` matched `/layers/nascent_counts`, Run 1 erroneously chose `/layers/nascent_counts` over `/X`, keeping only `59,775` cells, `1,662` NTCs (`831 / 831`), `198` targets (`missing RBM8A, RUVBL1`), `15,432` genes (`|G_eval| = 6,015`), `0` reproducible targets, and `0` eval targets.

We patched [`check_hek293_xucao_feasibility`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/audit/independent_r12_loader.py#L370-L407) to prefer integer `/X` over subfraction layers (`nascent`, `unspliced`, `spliced`, `ambiguous`, `velocity`, `intron`, `exon`), added unit test Case G in [`test_r12_audit.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/audit/test_r12_audit.py) (`pyflakes` OK, `4/4 OK`), and re-ran `hek293_xucao` (`task-3127`, `chosen_counts_path = "/X"`, `281.4s`).

### 1.2 Post-Fix Item-by-Item Comparison (`/X` Run vs `round12/shared/lines/hek293_xucao/`)

| Item | Independent `reader_xcheck` (`A-T3`, `/X`) | `round12/shared/lines/hek293_xucao` | Match? | Exact Difference / Notes |
|---|---:|---:|:---:|---|
| **Retained Cells** | **`98,274`** (of `98,315` raw; `41` dropped `<300` human UMI) | **`98,274`** | **PASS** | `diff = 0` (`max_pct_mt = None`, `apply_author_kd40 = False`) |
| **Total NTC Cells** | **`2,757`** (`2,757` sampled) | **`2,757`** (`2,757` sampled) | **PASS** | `diff = 0` |
| **Two-Half NTC Counts (`Half A / Half B`)** | **`1,379 / 1,378`** | **`1,378 / 1,379`** | **MISMATCH** | `Half A <-> Half B` parity swapped (`1,379/1,378` vs `1,378/1,379`); every odd target's `half_a_n_cells` in `xcheck` equals `shared`'s `half_b_n_cells` (`max diff = 0.0`). Split cells also differ cell-by-cell (`half_a_ctrl_mean_libnorm` diff `2.025` vs `shared_ha` and `6.686` vs `shared_hb`) |
| **`batch_labels` / `n_batches`** | **`2` (`["half_1", "half_2"]`)** | **`1` (`["batch_0"]`)** | **MISMATCH** | `A-T3` passed `half_1`/`half_2` as batches (`disjoint_batches = True` in `process_r10_line`); `shared` used single `batch_0` (`disjoint_batches = False`) |
| **`n_genes_final` (`genes.txt`)** | **`22,114`** | **`29,194`** | **MISMATCH** | `n_batches = 2` required finite SD in each half separately, dropping `7,080` genes |
| **`\|G_eval\|` (with `r9_6`)** | **`6,126`** | **`6,127`** | **MISMATCH** | `1` gene missing in `xcheck`: **`HINT2`** (in `only_in_shared`) |
| **`targets.txt` (`n_targets_deduplicated`)** | **`200`** | **`200`** | **PASS** | Exact ordered equality (`xcheck == shared`, `n_cells` array diff `0`) |
| **`n_targets_reproducible`** | **`23`** | **`96`** | **MISMATCH** | `diff = -73`; `ad_counts` quantiles `[0, 1, 2, 12, 2963]` (`>=10`: `53`) vs `shared` `[0, 22.75, 51, 145.25, 3595]` (`>=10`: `182`), while `energy_test_p_value <= 0.05` is identical (`189 == 189`) |
| **`n_eval_targets` (`r9_6` shared, rep, in `G_eval`)** | **`19`** (`< 50` confirmatory floor!) | **`62`** (`>= 50` confirmatory floor) | **MISMATCH** | `18` shared; `1` only in `xcheck` (`PCNA`); `44` only in `shared` (`ARID1A`, `BRD7`, `CASC3`, `CPSF2`, `CSTF3`, `CTCF`, `CYC1`, `DHX15`, `DHX29`, `DICER1`, `DIS3`, `ERCC3`, `ETF1`, `EXOSC2`, `EXOSC5`, `EXOSC6`, `GTF2A1`, `GTF2H4`, `INO80`, `KDM1A`, `METTL16`, `METTL3`, `MNAT1`, `NCBP1`, `NUDT21`, `PCF11`, `PELO`, `PGK1`, `PKM`, `POLD1`, `POLD2`, `POLD3`, `POLR1B`, `POLR1D`, `RFC5`, `SMARCA4`, `SMARCB1`, `STAG2`, `TAF5`, `UQCRFS1`, `WDR12`, `XRN2`, `YY1`, `ZC3H13`) |
| **`ctrl_mean_libnorm` / `_full` / `_log1p` Max Abs Diff** | `22,114` common / `6,126` `G_eval` | — | **PASS** | **`0.0`** (`max_abs_diff = 0.0`, exact agreement across all `22,114` & `6,126` genes) |
| **`half_a/b_ctrl_mean_libnorm` Max Abs Diff** | `22,114` common | — | **MISMATCH** | `half_a`: `2.024597`, `half_b`: `2.029434` (`6.6856` / `6.6904` if swapped) |
| **`pseudobulk_delta` Max Abs Diff (`delta_log1p` / `delta_libnorm`)** | `200 x 22,114` (`6,126` `G_eval`) | — | **MISMATCH** | `delta_log1p`: **`0.003620`** (`eval_space`: **`0.002769`**); `delta_libnorm`: **`0.254744`**; `half_a_delta_log1p`: **`1.419790`**; `half_b_delta_log1p`: **`1.419778`** |

### 1.3 Root Causes Identified for `hek293_xucao` Discrepancies
1. **Batch-Unit & Anderson–Darling Tie-Breaking Collapse (`batches=["half_1","half_2"]` vs `["batch_0"]`, `n_eval_targets = 19` vs `62`):**
   - [`independent_r12_loader.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/audit/independent_r12_loader.py) set `batches = ["half_1", "half_2"]` (`n_batches = 2`, `disjoint_batches = True`), whereas `round12/shared/lines/hek293_xucao` used single `batch_0` (`n_batches = 1`, `disjoint_batches = False`).
   - Inside `process_r10_line` (`prepare_r10_lines.py`), `n_batches = 2`: (i) dropped `7,080` genes finite only in the full `2,757`-NTC pool (losing `HINT2` from `G_eval`); (ii) centered full `delta_log1p` using separate half NTC means rather than pooled `2,757` NTCs (`diff = 0.00362`); and (iii) normalized `Z_5000` by `mu_top[half_1]` vs `mu_top[half_2]`, breaking discrete raw-zero expression ties across the two halves and collapsing Anderson–Darling FDR hits (`ad_counts` median `51 -> 2`, `n_targets_reproducible = 96 -> 23`, `n_eval_targets = 62 -> 19 < 50`), while `energy_test_p_value <= 0.05` stayed identical (`189 == 189`).
   - Crucially, inside `process_r10_line` when `batches = ["batch_0"]`, `half_a_delta_*` and `half_b_delta_*` **still** strictly subtract each half's own NTC mean (`ha_mu_log` / `hb_mu_log`), so single `batch_0` preserves half-isolated split deltas (`ha_mu`, `hb_mu`) without fragmenting `finite_mask` or breaking discrete zero ties in Addendum-B AD QC.
2. **Stage-2 `seed=42` Half-Split Difference:**
   - Even independently of `batch_0` vs `half_1/half_2`, `xcheck` (`perm[0::2] -> Half A`, `perm[1::2] -> Half B` on canonical-target barcode-sorted cells) produced `1,379 / 1,378` NTCs and `half_a_n_cells == shared.half_b_n_cells` (`max diff = 0.0`), whereas `shared` had `1,378 / 1,379` AND non-identical half control/delta arrays (`half_a_ctrl_mean_libnorm` diff `2.025` vs `shared_ha`, `6.686` vs `shared_hb`), indicating different cell-level split assignment / permutation ordering under `seed=42`.

---

## 2. `a549_ko` (`GSE345058`) Cross-Check — **EXACT ON `eval_targets`/`G_eval`/`ctrl_mean`, `QARS` (`591` Cells) UNAPPROVED-FILTER DIFFERENCE**

Before running `a549_ko`, an initial `task-3159` run exposed a column-case bug in `build_a549_ko_audit_parts` (`KeyError: 'batch'` because `GSE345058_cell_annotations.csv.gz` uses capitalized `"Batch"`). We fixed [`independent_r12_loader.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/audit/independent_r12_loader.py) with case-insensitive column resolution validated before 10 GB `.mtx.gz` streaming, deducted `Cas9` vector counts from per-cell UMIs (`umi_human = umi_all - cas9_mtx_umi`), switched MTX ingestion to `pd.read_csv(engine="c")`, and updated [`test_r12_audit.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/audit/test_r12_audit.py) `test_2` (`pyflakes` OK, `4/4 OK`).

### 2.1 Item-by-Item Comparison (`round12/audit/reader_xcheck/a549_ko` vs `round12/shared/lines/a549_ko/`)

| Item | Independent `reader_xcheck` (`A-T3`) | `round12/shared/lines/a549_ko` | Match? | Exact Difference / Notes |
|---|---:|---:|:---:|---|
| **Retained Cells** | **`606,075`** (`0` excluded) | **`605,484`** (`591` excluded) | **MISMATCH (`+591`)** | Difference equals **`QARS`** (`591` cells; `hgnc.resolve("QARS", "") == ("QARS", "unmapped")`), excluded in `shared` (`"n_nonplain_or_unapproved_excluded": 591`) vs retained per `PREREG_E §3.1` (`606,075 / 606,075`) in `xcheck` |
| **Total NTC Cells** | **`45,320`** (`20,000` sampled) | **`45,320`** (`20,000` sampled) | **PASS** | `diff = 0` (`32` Lanes nested inside `batch1..4`) |
| **Two-Half NTC Counts (`Half A / Half B`)** | **`23,758 / 21,562`** | **`23,758 / 21,562`** | **PASS** | Exact equality (`Half A = batch1+batch3`, `Half B = batch2+batch4`) |
| **`n_genes_final` (`genes.txt`)** | **`17,337`** | **`17,337`** | **PASS** | `genes_txt_exact_equal = true` (same `17,337` symbols in exact order) |
| **`\|G_eval\|` (with `r9_6`)** | **`6,121`** | **`6,121`** | **PASS** | `exact_equal = true` |
| **`targets.txt` (`n_targets_deduplicated`)** | **`1,000`** | **`999`** | **MISMATCH (`+1`)** | `only_in_xcheck = ["QARS"]`; all `999` other symbols identical |
| **`n_targets_reproducible`** | **`972`** | **`971`** | **MISMATCH (`+1`)** | `QARS` passes Addendum-B QC (`n_cells=591`, `ad_counts=2582`, `ep=0.0005`); all `999` shared targets match identically (`971 == 971`) |
| **`n_eval_targets` (`r9_6` shared, rep, in `G_eval`)** | **`445`** | **`445`** | **PASS** | `exact_equal = true` (`100%` identical ordered `445`-target evaluation set; `QARS` is absent from `r9_6` and `G_eval`) |
| **`ctrl_mean_libnorm` / `_full` / `_log1p` / `half_a` / `half_b` Max Abs Diff** | `17,337` & `6,121` `G_eval` | — | **PASS** | **`0.0`** across every control mean array (`ctrl_mean_libnorm`, `_full`, `_log1p`, `half_a`, `half_b`, `eval_space_ctrl_mean_libnorm`) |
| **`pseudobulk_delta` Max Abs Diff (aligned `999 x 17,337` & `445 x 6,121`)** | `999 x 17,337` (`6,121` `G_eval`) | — | **FLOAT32 EPSILON DIFF** | `aligned_delta_log1p`: **`3.576e-7`** (`allclose_1e_6 = true`); `eval_space_delta_log1p`: **`3.576e-7`** (`allclose_1e_6 = true`); `half_a_delta_log1p`: **`4.172e-7`**; `half_b_delta_log1p`: **`4.768e-7`**; `aligned_delta_libnorm`: **`1.526e-5`** (due to `cas9_mtx_umi` deduction from per-cell library size / float32 ULP) |

---

## 3. Pre-Freeze Summary & Governance Impact

- **`a549_ko` (`GSE345058`, Mandatory Primary Line):** Independently verified to **`445 == 445`** eval targets, **`6,121 == 6,121`** `G_eval` genes, **`23,758 / 21,562`** two-half NTCs, **`0.0`** control mean diff, and **`<= 3.58e-7`** `delta_log1p` agreement across all `999` shared targets (`diff` on `QARS` `591` cells & `cas9_mtx_umi` library-size exclusion has zero impact on `445` evaluation targets).
- **`hek293_xucao` (`XuCao2023.h5ad`, Conditional Line):** **Critical discrepancy (`19 < 50` vs `62 >= 50` eval targets, plus Stage-2 split cell differences and `batch_0` vs `half_1`/`half_2` AD tie-breaking)** requires a Master Controller ruling before `CODE_FREEZE_R12` / holdout evaluation (specifically: whether Stage-2 HEK293 passes `"batch_0"` to `process_r10_line` per `shared/lines/hek293_xucao` — where `ha_delta`/`hb_delta` remain strictly half-NTC-centered while pooled `Z_5000` preserves discrete zero ties — and aligning Stage-2 `seed=42` barcode parity `1,378/1,379` vs `1,379/1,378`).

