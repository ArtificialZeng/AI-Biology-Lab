# DEVIATION_R17_D2 — control pools must come from the builder that made the line's observed data (telohaec)

Written by Claude (master controller), 2026-11-27 11:01 +0800, during the D1 controls-only pre-checks. No Round 17 array, metric or summary exists from either party.

## What happened
- **A-T3's pre-check, 10:40:** the telohaec control extraction raised `ValueError: telohaec NTC pool excluded: N_NTC=0 < floor 15`.
- **Both parties' frozen code** gets telohaec parts from the Round 11 library reader, `round11/src/dataset_readers_r11.py:build_telohaec_parts`. On the raw files that reader resolves no NTC cells. A-T3's diagnosis: it takes the first dial-out column (`gene_symbol` instead of the cell-barcode column), misses `guide_w_max`, and matches against stringified `obs["_index"]` tuples.
- **The line's observed data** (`round11/shared/lines/telohaec`: `n_nt_total` 2,215, `n_batches` 20) was built by the Round 11 script, `round11/shared/scripts/prepare_r11_lines.py:build_telohaec_parts`. A-T3 verified that this script reproduces 180,574 cells and 2,215 NTCs across 20 lanes.
- The main line's pre-check would fail on the same pool, because both parties call the same reader.

## Ruling (binding on both parties)
1. **General rule.** For every Round 17 line, the sampled control pool must be drawn from the same NTC pool as the line's standardized observed data. The total NTC count and the number of batches must equal the line's `line_metadata.json` (`n_nt_total`, `n_batches`), and the batch partition must be the same. Each party's pre-check must report `n_ctrl_total` and the number of batches per pool, so the controller can check this for all 17 pools.
2. **telohaec.** Both parties obtain the parts from `round11/shared/scripts/prepare_r11_lines.py:build_telohaec_parts(raw_dir, hgnc)`, the builder of `round11/shared/lines/telohaec`. Sampling (seed 20261701, stratification, QC) is otherwise unchanged. This changes only the telohaec branch of each party's control extraction.
3. **Barcode strings.** If a party cannot resolve real cell barcodes for a pool (for example, A-T3's resolver falls back to `part_row_<idx>` on hek293_xucao and calu3, whose raw `X` is CSC), then cell identity for verification step 1 is established by these, all exactly equal between the parties:
   - the raw-count matrix sha256;
   - the per-cell totals;
   - the batch labels;
   - n_ctrl_total.

   Barcode strings are still compared wherever both parties resolve them.
4. **Any other pool** whose n_ctrl_total or batch count differs from its line metadata is a deviation to report, not to work around.

## What does not change
PREREG_J, C1–C3, D1, the line roles (telohaec stays descriptive micro-N), and every other line.

## D2.5 — per-cell normalization totals (2026-11-27 12:50 +0800, controller; from the pre-check comparison, before any re-freeze)
- **Pre-check result** (`controller/checks/r17_compare_precheck.py`): on all 17 pools the two parties agree on n_ctrl_total (equal to every line's metadata), the batch labels, the gene axis and the raw-count sha256. Barcodes are equal wherever both parties resolve them. The one difference is telohaec, where 13 of 512 cells have per-cell totals differing by ≤ 2 UMI.
- **Cause.** A-T3's stored `cell_total_umis` comes from the parts-level totals, which include features the exported matrix does not carry. The main line's normalization uses the frozen R16 preprocessing, which takes row sums of the exported matrix.
- **Ruling.** For every pool, the per-cell totals used for normalization (`official_median_log1p` totals and their median, and `cell_total_umis`) are the row sums of the exported raw-count matrix: the matrix both parties feed to SE-600M and to HVG preprocessing. On the other 16 pools this is already true for both parties. A-T3 sets `cell_total_umis = raw_counts.sum(1)` (float64) and recomputes the median from it. The main line needs no change.
