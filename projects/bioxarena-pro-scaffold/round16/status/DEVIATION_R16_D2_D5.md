# DEVIATION_R16_D2 – D5: input-format fixes found by a controls-only dry run (before any Round 16 result)

Written by Claude (master controller), 2026-11-26 17:10 +0800. These follow DEVIATION_R16_D1. No Round 16 prediction, metric or array exists yet. Both runs stopped before inference: A-T3 at 15:56, the main line at 16:02 on a controller hold.

## How these were found
After A-T3 stopped at feng_ipsc, the controller ran a dry run of raw-control extraction only (no model, no prediction) for all six lines, with both parties' code: controller/checks/r16_dryrun_controls.py and r16_dryrun_compare.py. Outputs were written to the controller scratch folder, not the repo. The controller also checked the two inputs that are only read after inference: Arc's archive file and the checkpoint's batch map.

## D2 — feng_ipsc cell IDs (both parties)
- **Problem:** metadata Cell_ID looks like `PC-P4-D3_I73_<bc>`; RNA-UMI-Counts headers look like `P4_I73_<bc>`. Exact matching found 0 of 1,161,864.
- **Rule:** the mapping round11 prepare_r11_lines.py L46-49 used to build the feng_ipsc line data, `re.sub(r"^PC-(P\d+)-D\d+_(I\d+_.*)$", r"\1_\2", Cell_ID)`. The controller verified it is one-to-one: 1,161,864 IDs give 1,161,864 distinct headers, and all are present.
- **Unchanged:** enumeration stays in CSV header order. The exported barcode is the CSV header string, as locked Addendum A section 7 requires. A collision raises.
- **Correction on record:** the controller's first D2 message to A-T3 wrongly said to export the metadata form. That was corrected within the hour (controller/msgs/r16_D2_AT3_correction.txt), before any run.
- **Implementations:** main line by the controller (run_eval_r16.py); A-T3 by A-T3 in its own code (audit_r16_raw_controls.py).

## D3 — Nadig legacy AnnData categoricals (main line)
- **Problem:** in both GSE264667 files, obs `gene` and var `gene_name` are legacy categoricals: integer codes plus `__categories`. The main-line reader returned the raw codes, so it found 0 non-targeting cells. Gene names would also have come through as integers.
- **Fix:** add the legacy branch exactly as round9 prepare_lines.read_obs_or_var_col does (the round9 code built these lines and is used by A-T3). Modern categoricals and plain columns read exactly as before.

## D4 — Nadig exported barcodes (A-T3)
- **Problem:** A-T3 exported synthetic row ids (`nadig_hepg2_row_31629`). Addendum A section 7 requires the AnnData obs index resolved via obs.attrs["_index"], which is `cell_barcode` here (unique: 145,473 HepG2 / 262,956 Jurkat).
- **Same cells already:** before the fix, the controller showed A-T3's sampled cells were already identical to the main line's. Raw counts were array_equal, and the row ids map one-to-one onto the main line's barcodes, for the primary and resample sets.
- **Fix:** A-T3 checks obs.attrs["_index"] first.

## D5 — Arc archive gene index (main line; Addendum A section 3 check)
- **Problem:** the var index of Arc's archived adata_pred.h5ad is positional ('0'..'1999'). The check matched genes by name, so it would have found 0 common genes, failed, and halted the run after nadig_hepg2.
- **Evidence:** under the checkpoint's ordered G_model (var_dims.pkl gene_names), Arc's predicted control profile correlates with observed HepG2 NTC controls at r = 0.997 across 2,000 genes, against r = 0.017 for alphabetical order and -0.031 for random order.
- **Fix:** if the archive var index is exactly '0'..'n-1', map it by position onto var_dims.pkl gene_names; a length mismatch raises.
- **Unchanged:** the pass thresholds (median >= 0.90, pooled >= 0.95) and the comparison protocol.
- **Test (controller/checks/r16_d5_archive_test.py, Arc's own deltas fed back through the patched check):**
  - correct names pass with median and pooled Pearson 1.0 on 60 perturbations x 500 genes;
  - shuffled names fail (r about 0.01).

## Also checked, no change needed
- The checkpoint batch map uses names '1'..'56', the same strings as Nadig gem_group labels (native-batch arm of the section 3 check).

## What does not change
Hypotheses, metrics, seeds, thresholds, cohorts, sampling rules, data files, models, and adjudication rules all stay as locked. Each fix only makes the code read the files as the locked protocol already specifies (Addendum A section 7 identifiers; PREREG_I section 1.3 check).

## Verification before release
- Dry-run cross-party comparison. Required: exact ordered barcodes, identical batches, identical raw counts, identical resample sets, for all 6 lines. Results are appended below when complete.
- Tests:
  - main-line unit suite test_unit_r16: 18/18 OK (rerun on final code e2a2f2cf...);
  - controller D1 test: pass;
  - controller D5 test: pass;
  - A-T3 test_audit_r16: 13/13 OK (A-T3 report).

## Dry-run results (controller, 2026-11-26 17:25 +0800; final code)
Main line (run_eval_r16.py e2a2f2cf...) vs A-T3 (audit_r16_raw_controls.py aede3b19...), primary seed 20261601 and resample seed 20261602:

| line | n | genes | barcodes (ordered) | batches | raw counts | resample barcodes | resample counts |
|---|---|---|---|---|---|---|---|
| a549_ko | 512 | 29,891 | identical | identical | array_equal | identical | array_equal |
| cd4t_rest | 512 | 18,130 | identical | identical | array_equal | identical | array_equal |
| dld1 | 512 | 36,601 | identical | identical | array_equal | identical | array_equal |
| feng_ipsc | 512 | 36,518 | identical | identical | array_equal | identical | array_equal |
| nadig_hepg2 | 512 | 9,624 | identical | identical | array_equal | identical | array_equal |
| nadig_jurkat | 512 | 8,882 | identical | identical | array_equal | identical | array_equal |

Final-code tests: main-line test_unit_r16 18/18 OK; controller r16_d1_reorder_test PASS; controller r16_d5_archive_test PASS; A-T3 test_audit_r16 13/13 OK.
