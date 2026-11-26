# PREREG_J (Round 17): controller lock record

Claude (master controller), 2026-11-27 06:22 +0800.

## What is locked
The following files, byte for byte, are PREREG_J from this commit on. Their content stays uncommitted until the controller has verified A-T3's Round 16 pass 2, because they quote Round 16 discovery values. These SHA-256 values are the binding commitment, and the files will be committed unchanged.

| File | SHA-256 | Bytes |
|---|---|---|
| `round17/prereg/PREREG_J_DRAFT_v2.md` (the protocol) | `7033337d1f42588ee9277804c6a803eb940312fe95a4ee31e2c000555b0fb85b` | 67,941 |
| `round17/prereg/PRIOR_ART_J.md` | `3f9210df9c4a42dc17e68f08d2133c130f5265693ec2d0bd24558b5840180b75` | 21,687 |
| `round17/status/LINE_INVENTORY_J.json` | `102faba3d6df2b6c77b89dfb74c77c96de294cdfe589f69a00434ce3562e81b2` | 30,100 |
| `round17/src/build_line_inventory_j.py` | `3dce1a2f6dfa840179cd1cbef4a4d430aa70a68212b5191f6e043ed55ad96a39` | 31,884 |
| `round17/src/pfizer_telohaec_reader_r17.py` | `35e41ab80d513242e8710a5c0eb52e9ffd1a16b4710975927a145df6ef6683c3` | 41,756 |
| `round17/src/test_pfizer_telohaec_reader_r17.py` | `05edefceec7a517a2a32eb485fba6490b64adebe734e3dc766e58fbf04f89941` | 18,541 |

History, not binding: `round17/prereg/PREREG_J_DRAFT.md` (v1), `3f31bce74be8003d6412af32fd9674b86272bee6f7cbb1e7b28cba6c6146e90b`.

## Review trail
- v1 was reviewed by the controller (review R1, items J1–J7 plus three rulings). v2 resolves all of them.
- Controller checks on v2:
  - the TeloHAEC reader's synthetic tests rerun by the controller: 4/4 pass;
  - the hashes above equal the main line's `R17_PREREG_V2_READY` report.

## State at lock (no Round 17 inference of any kind has run)
- **Confirmatory units:** cd4t_stim48 (local), gse295214_h1 and pfizer_telohaec_basal. Each enters only under PREREG_J §2.1 and the PREREG_I floors.
  - gse295214_h1: R16 counts-only QC PASS_GROUP3_PRIMARY (`round16/shared/R16_GSE295214_QC_COUNTS.md`; data freeze `round16/status/FREEZE_MANIFEST_R16_DATA_G3.txt`). There has been no STATE or Occam inference on it. The R16 Group 3 GO is issued only after this lock.
  - pfizer_telohaec_basal: `/Users/mac/.cache/r17_telohaec_pfizer/imru_full.h5ad`. Size 53,290,123,849 B; md5 `0594ddd25e224e5e9cc22414e71daf48` (equal to Zenodo); SHA-256 `4b22ad11de666f7a9b3710740995da8efc916a73991a3615720df87017550565`. It has not been opened except for the structure check below. Counts-only QC (`prepare_pfizer_telohaec_r17`) runs only after this lock.
- **TeloHAEC structure-only check** before the lock (`controller/checks/r17/r17_telohaec_struct_check.py`; output `r17_telohaec_struct_out.json`, sha256 `faa38f5d…`). It read keys, shapes, dtypes, column names, condition categories and control-label flags; no matrix values and no per-target counts.
  - Every column the reader resolves exists: `condition` (categories `Treated`, `Untreated`), `gene_target`, `guide`, `flask` (2 levels), `welltag` (47 levels). There is no guide-count column, which the reader treats as optional.
  - The only non-gene labels in `gene_target` are `NO-TARGET` (the reader's NTC rule) and `SAFE_TARGET` (the reader's exclusion rule). There are no multi-guide or unassigned labels.
  - `X` is CSR float32, 864,115 × 38,606, with no `layers` and no `raw`. The integer-count check is the reader's post-lock fail-closed slot rule.
  - `var` has only its index, which mixes HGNC symbols with Ensembl IDs (12,962) where no symbol exists. The frozen R10 `HGNCLookup.resolve` maps approved symbols directly and Ensembl IDs through its Ensembl table, so the reader needs no change.
  - **Result: the reader spec matches the file; v2 is locked unchanged.**

## Sequencing (binding)
1. This lock commit is pushed.
2. The R16 Group 3 GO can follow.
3. Model-free, after the lock:
   - TeloHAEC counts-only QC, then its data freeze;
   - cd4t_stim48 and the descriptive lines, standardized as already done.
4. Round 17 code freeze and A-T3 release (after R16 closes), then `R17_GO`. No Round 17 STATE inference runs before `R17_GO`.
5. The full round17/ text is committed unchanged once A-T3's R16 pass 2 is verified. Its hashes must equal the table above.
