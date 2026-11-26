# DEVIATION_R16_D7 — A-T3 SE-600M gene mapping: duplicate raw symbols mapped twice

Written by Claude (master controller), 2026-11-26 22:45 +0800.

## How it was found
- The main line finished Groups 1–2 at 22:07 (`R16_MAIN_REAL_DONE`).
- The controller then ran Addendum A §11 item 5, staged verification step 1: compare X_state embeddings first. The script is `controller/checks/r16_compare_embeddings.py`.
- It compared main `round16/results/embeddings/x_state_*.npz` with A-T3 `round16/audit/results/se600m_embeddings_r16.npz`. That file was written 20:35 by the run that later stopped on D6; its SE-600M code is the same as in the run started at 20:47.

| Pool | Ordered barcodes equal | max \|X_main − X_AT3\| |
|---|---|---|
| a549_ko, a549_ko resample | yes | 0 (bitwise) |
| cd4t_rest, cd4t_rest resample | yes | 0 (bitwise) |
| dld1, dld1 resample | yes | 0.049, 0.042 |
| feng_ipsc | yes | 0.052 |
| nadig_hepg2, nadig_hepg2 resample | yes | 0.044, 0.048 |
| nadig_jurkat, nadig_jurkat resample | yes | 0.094, 0.095 |

Tolerance is 1e-5. The four failing lines are exactly the lines whose raw gene lists contain duplicated symbols. From the controller dry run: dld1 and feng_ipsc have 10 each; each Nadig line has 1.

## Cause
- arc-state 0.11.1 `FilteredGenesCounts` (`state/emb/data/loader.py:223-226`) maps genes to SE-600M global indices by `adata.var_names`. It falls back to `gene_column` only if no var_name maps.
- Addendum A §11 item 3 locks the mapping:
  - set `adata.var_names` to the raw symbols;
  - then call `adata.var_names_make_unique()`;
  - so the first occurrence of each duplicated symbol is kept and no SE-600M global ID appears twice.
- **Main line (`round16/src/model_adapters_r16.py:715-718`):** follows this rule.
- **A-T3 (`round16/audit/audit_r16_state_adapter.py::encode_with_loaded_se_model`):** set `adata.var_names = gene_syms` without make_unique, so every copy of a duplicated symbol mapped to the same global ID.
- **Mapping check without the model:** `controller/checks/r16_d7_mapping_test.py` builds the `FilteredGenesCounts` map both ways on the dry-run gene lists:

| Line | SE-600M genes mapped (main / A-T3) | Duplicate global IDs (main / A-T3) | Columns mapped differently |
|---|---|---|---|
| a549_ko | 18,098 / 18,098 | 0 / 0 | 0 |
| cd4t_rest | 17,958 / 17,958 | 0 / 0 | 0 |
| dld1 | 19,324 / 19,332 | 0 / 8 | 8: TBCE, CYB561D2, MATR3, HSPA14, GOLGA8M, GGT1, ARMCX5-GPRASP2, TMSB15B |
| feng_ipsc | 19,279 / 19,287 | 0 / 8 | same 8 |
| nadig_hepg2 | 9,192 / 9,193 | 0 / 1 | 1: HSPA14 |
| nadig_jurkat | 8,441 / 8,442 | 0 / 1 | 1: HSPA14 |

The lines with identical maps are exactly the lines with bitwise-identical embeddings.

## Ruling and change
- 22:35: HOLD and ruling sent to A-T3 (`controller/msgs/r16_D7_AT3.txt`).
- A-T3 stopped its pass-1 run (pid 86771) and deleted the stale embedding file.
- The change is in A-T3's own code: one line, `adata.var_names_make_unique()`, added right after `adata.var_names = gene_syms`.
  - Controller checked it with `git diff`: 1 insertion, nothing else.
  - The other frozen files are unchanged and verify against FREEZE_MANIFEST_R16_D3.
  - New sha256: bd47c9f01922c958a48d1de051d8edd9fd16aac24bd78f49101e7e85def98e65.
- A-T3 verification:
  - py_compile passes.
  - test_audit_r16: 13/13 OK.
  - Duplicate-symbol self-check `D7_DUPLICATE_SYMBOL_SELFCHECK_OK`: only the first copy maps, there are no duplicate global IDs, and the mapped count equals the number of unique mapped symbols. It also checked an end-to-end `encode_cells` run.
- A-T3 restarted pass 1 from scratch at about 22:43 (pid 34240), recomputing all 11 embeddings. Pass 2 follows.

## What does not change
- Nothing changes in the preregistration, Addendum A, data, models, or the main line. The main line is compliant.
- No metric was involved in finding or fixing this. A-T3 was not shown any main-line metric or status. The ruling message names only the agreement and disagreement of embeddings by line, and the gene symbols involved.
- The A-T3 rerun will be checked again with step 1 (embeddings ≤ 1e-5) before predicted Δ (≤ 1e-4) and metrics (≤ 1e-5).
