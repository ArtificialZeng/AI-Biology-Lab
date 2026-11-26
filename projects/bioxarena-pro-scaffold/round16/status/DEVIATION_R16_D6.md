# DEVIATION_R16_D6 — A-T3 result export: JSON serialization of numpy values

Written by Claude (master controller), 2026-11-26 20:48 +0800.

## What happened
- **Crash:** 20:37 A-T3's pass-1 rerun stopped with `TypeError: Object of type ndarray is not JSON serializable` in `run_audit_r16.py::export_line_arm_npz`, at `json.dumps(arm_metadata, sort_keys=True)`. Real-data arm metadata (sampling audit, gene mapping, centroid diagnostics) contains numpy arrays.
- **State at the crash:** all 6 raw control pools were extracted, all 11 SE-600M embeddings computed, and STATE_SE inference run on nadig_hepg2. No metric had been computed; the crash came before `evaluate_and_adjudicate_line`. A-T3 remains blind to results.
- **Main line:** unaffected, still running.

## Change (A-T3, its own code: round16/audit/run_audit_r16.py, 20 insertions, 4 deletions)
- **The converter, `_json_default_r16`:**
  - np.ndarray → list;
  - np.generic → Python scalar;
  - Path → str;
  - sets and tuples → list (sorted if all elements are strings);
  - anything else raises TypeError.
- **Where it applies:** passed as `default=` to the five result-writing json.dumps calls (arm_metadata_json, both control-barcode registries, *_metrics.json, manifest_r16_audit.json).
- **Unchanged:** no value, key, sort order, computation or sampling code.
- Controller reviewed the diff. New sha256: 6ebf9992916ad63fded30be089326f908a33c61ecef5e486a3b02c1d31d4332f.

## Verification
- A-T3 test_audit_r16: 13/13 OK.
- A-T3 serialization self-check D6_SERIALIZATION_CHECK_OK: numpy arrays and scalars, Path, set/tuple, TypeError on an unknown object, and an export round-trip.
- A-T3 reran pass 1 at 20:47 with the identical command.

## What does not change
Nothing in the preregistration, Addendum A, data, models, or the main line. Only the JSON metadata files A-T3 writes become serializable.
