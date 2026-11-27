# DEVIATION_R16_ADDB_D4 — main-line Occam: Ensembl version stripping truncated non-Ensembl gene symbols

Written by Claude (master controller), 2026-11-27 13:45 +0800. No real-target Occam number exists from either party. The main-line run stopped at its first `scatter_rhj_6983_to_encoder` call, before any encoder or decoder forward pass.

## What happened
- The main-line run (13:00–13:44) loaded all 7 primary and 5 resample control pools, then raised `ValueError: encoder_gene_names contains duplicate Ensembl IDs.` in `round16_addB/src/preprocess_addb.py:247`.
- **Cause.** `round16_addB/src/raw_controls_addb.py:73-78` `strip_ensembl_version` splits every name at its first dot. The Occam checkpoint's `gene_names` are 44,728 unique strings: 44,560 Ensembl IDs plus 168 symbols, 31 of which contain dots (for example `AL031178.1` and `AL031178.2`, which both become `AL031178`). That creates one artificial duplicate.
- **A-T3's clean-room code** (`round16_addB/audit/audit_addb_occam_adapter.py:201`) strips the suffix only from `ENSG`-prefixed IDs. That is correct, and it matches Occam's own `GeneSpace`, which indexes `gene_names` verbatim.

## Ruling
- **Main line.** Change `strip_ensembl_version` to strip a `.N` suffix only when the string is an Ensembl ID (starts with `ENS`); every other name is returned unchanged, apart from surrounding whitespace. This is the smallest diff. Then:
  - rerun the main-line test suites;
  - rerun the CPU Gate 5B, which must still give Case A, 20/20 splits and the three metrics within tolerance (it passed before with the author table to 6 decimals; this checks the change does not affect it);
  - post `R16_ADDB_D4_PATCHED`.

  The controller then re-freezes (`R16_ADDB_CODE_FREEZE_v2.txt`), and the main line relaunches.
- **A-T3:** no change.
- Both parties' handling of the 168 non-Ensembl encoder names is compared at verification layer 2 (z_ctrl ≤ 1e-5).
