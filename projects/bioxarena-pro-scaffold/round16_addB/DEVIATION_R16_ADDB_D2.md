# DEVIATION_R16_ADDB_D2 — A-T3 Occam runner: startup bugs found before any real-target inference

Written by Claude (master controller), 2026-11-27 13:05 +0800.

## What happened
A-T3 ran the exact `R16_ADDB_GO_AT3.txt` command at 13:04. It exited with code 1 inside its own freeze gate, before any model call: `AddBBlindingGateError: Addendum B draft v2 SHA-256 mismatch`. A-T3 then audited the real-data path it had not exercised in its synthetic tests, and found three defects, all in its own frozen code (`round16_addB/audit/`):
1. **The crash.** `audit_addb_controls_and_metrics.py:71`: `EXPECTED_ADDB_V2_SHA256` contains a one-character typo (`…178ad7d7…` instead of `…178ad5d7…`). The real hash of `PREREG_I_ADDENDUM_B_DRAFT_v2.md` is `537467db51b00b8e178ad5d7638e25b5d809f350c8feae229d86cdb96acb281e`, the value in the lock record and in `R16_ADDB_CODE_FREEZE_AT3.txt`. The unit test that should have caught this skipped the hash branch.
2. **A latent crash.** `run_audit_addb.py:1467` calls `load_rhj_gene_axis(DEFAULT_GENE_AXIS_CSV, verify_sha256=True)` and unpacks three values. The frozen adapter's signature is `load_rhj_gene_axis(gene_axis_csv, expected_content_sha256=…) -> (ensembl_ids, symbols)`. The unused fallback at lines 1438–1439 has the same kind of mismatch, against the `run_stage_5b_cpu` signature and return key.
3. **A wrong descriptive field.** `audit_addb_controls_and_metrics.py:342, 429` read `n_kept_ntc` at the top level of `round16/shared/lines/gse295214_h1/line_metadata.json`. There the top-level key is `n_nt_total` (58,153); `n_kept_ntc` is nested under `adapter_audit`. As a result, `n_ctrl_total` and `ctrl_sampling_fraction` for gse295214_h1 fall back to 512 and 1.0. These are descriptive fields, not inputs to any metric.

## Ruling
- Approved, as the smallest diff, in A-T3's own files only:
  - item 1: correct the hash constant;
  - item 2: call `load_rhj_gene_axis(DEFAULT_GENE_AXIS_CSV)` with the frozen signature and unpack two values, and align the fallback with `run_stage_5b_cpu(adapter, batch_size=…)` / `overall_cpu_gates_pass`;
  - item 3: read `n_nt_total` (top level), falling back to `adapter_audit.n_kept_ntc`.
- Add a unit test that runs the freeze gate against the real `round16_addB/` files, so the hash branch is exercised. Rerun the full suite and the synthetic end-to-end recompute.
- The controller re-freezes as `R16_ADDB_CODE_FREEZE_AT3_v2.txt` before the relaunch. No real-target Occam number existed from either party when this was written.
- **Main line, for parity:** confirm how your code sets `n_ctrl_total` for gse295214_h1. It must be 58,153 (line metadata `n_nt_total`).
