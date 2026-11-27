# DEVIATION_R17_D1 — multi-file lines with different gene axes (a549_myllia)

Written by Claude (master controller), 2026-11-27 09:54 +0800. No Round 17 array, metric or summary existed from either party when this was written.

## What happened
- **Main-line run, 09:52.** The run crashed after 1,442 s, while extracting controls for the third pool, a549_myllia: `ValueError: Gene axis mismatch across parts: a549_myllia_ch_2 != a549_myllia_ch_1`, raised by the frozen R16 `extract_ntc_bundle_from_parts_r16`. It wrote nothing except the log. The main line halted as required.
- **Cause.** GSE337804 ships two Seurat `.RDS.h5ad` channel exports, and each dropped low-detection genes independently:
  - ch_1: 18,711 genes;
  - ch_2: 18,808 genes;
  - shared: 18,480; ch_1 only: 231; ch_2 only: 328.
- **How the observed side was built.** The Round 14 standardized line (`round14/shared/lines/a549_myllia`) was built by frozen R10 `process_r10_line`. That code takes the common-key intersection across parts: the key is the gene name, or the gene id if there is no name, in part-0 order. Its metadata: `raw_shape [51448, 18480]`, `n_genes_multifile_dropped 231`.
- **A-T3's frozen code** (not yet reached in its run) uses a third rule: ch_1's full axis, zero-filled for ch_2 cells, with ch_2-only genes dropped. The two parties would therefore have produced different a549_myllia control matrices.

## Ruling (binding on both parties)
1. For any line whose parts (files or channels) have different gene axes, the sampled control raw-count matrix uses the R10 common-key intersection: the key is the stripped gene name, or the gene id when the name is empty, in part-0 order. This is exactly the axis the line's observed side was built on (`process_r10_line`). Of the Round 17 lines, only a549_myllia is affected.
2. Per-cell totals for that line are row sums over the intersected axis: `cell_total_umis`, the `official_median_log1p` normalization totals and its median. This applies the PREREG_I §2.1 / Addendum A rule (the full measured transcriptome; DEVIATION_R16_D8) to the transcriptome both channels measured.
3. Control sampling itself is unchanged: the same QC, stratification, seed 20261701, and cell order from the frozen parts.
4. **Main line:** patch only the control extraction path in `round17/src/ntc_extractors_r17.py`, to restrict parts to the common axis before calling the frozen R16 extractor (or an equivalent that gives identical output). Do not modify `round16/`.
5. **A-T3:** patch only `_extract_from_parts_with_barcode_lookup` (or its caller) in `round17/audit/audit_r17_controls.py`, to use the intersected axis and intersected row sums whenever the part axes differ.
6. **Before either party restarts Round 17,** each runs a controls-only pre-check: no SE-600M, no STATE, no metric. It calls its own frozen-plus-patched control extraction for all 17 pools (11 lines plus the 6 R16 mismatched pools) at seed 20261701, and writes per pool:
   - the ordered barcodes;
   - the batch labels;
   - the per-cell totals;
   - the gene axis (names and length);
   - the sha256 of the raw-count matrix.

   The controller compares the two parties, after the C2.4 prefix normalization. Both parties restart only after a full match and a re-freeze.

## What does not change
PREREG_J, C1–C3, every hypothesis and threshold, the line roles (a549_myllia stays secondary descriptive), and every other line's processing.
