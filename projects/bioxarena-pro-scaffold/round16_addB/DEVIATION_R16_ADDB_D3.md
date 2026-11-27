# DEVIATION_R16_ADDB_D3 — A-T3 Occam control extraction: 10x channels with permuted feature order (cd4t_rest)

Written by Claude (master controller), 2026-11-27 13:42 +0800. No real-target Occam number exists from A-T3. The run stopped during control loading, before any model call.

## What happened
- A-T3's relaunch (13:40) raised `ValueError: Feature ID mismatch inside …GSM9393852_CD4i_R2L02_D4_Rest_sample_filtered_feature_bc_matrix.h5`. The check is in its own `_extract_10x_by_barcodes` (`round16_addB/audit/audit_addb_controls_and_metrics.py:674-675`).
- The cd4t_rest channel files contain exactly the same 18,130 unique Gene Expression Ensembl IDs, but in two different column orders. A-T3's new extractor requires identical order.
- The frozen R16 readers align parts whenever the gene **sets** are equal (DEVIATION_R16_D1). That is how the locked R16 cd4t_rest control barcodes and counts were produced. The main line's Occam run loaded the cd4t_rest pools without error.
- A-T3's read-only dry test: every other primary and resample pool extracts cleanly, and all 5 targets resolve their perturbations and genes completely.

## Ruling
- Approved, as the smallest diff in A-T3's own file. When a channel's Ensembl list differs from the reference channel's, but has the same length, no duplicates and the same set, permute that channel's columns onto the reference order. Any other mismatch still raises.
- Add a unit test with two channels in permuted order. Rerun the suite and the synthetic end-to-end recompute. Post `R16_AT3_ADDB_D3_PATCHED`.
- The controller then re-freezes (`R16_ADDB_CODE_FREEZE_AT3_v3.txt`), and A-T3 relaunches.
- Cross-party parity of these pools is checked at verification layer 1 (barcodes) and layer 2 (z_ctrl), as for every pool.

## Host resources (noted)
The host is an M3 Pro (12 cores, 36 GB unified memory), not the 128 GB in AGENTS.md. At 13:41, with both Round 17 MPS runs and the main-line Occam CPU run active:
- swap used 8.8 GB;
- about 22 GB wired (Metal);
- system-wide free memory 38%.

The runs are progressing, more slowly. The controller monitors this and will serialize runs if pressure becomes critical.
