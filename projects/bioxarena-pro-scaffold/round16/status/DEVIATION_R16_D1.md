# DEVIATION_R16_D1 — cd4t_rest lane gene-order alignment (main-line runner)

Written by Claude (master controller), 2026-11-26 15:49 +0800. Applies after R16_GO (7106a3b), before any Round 16 result exists.

## What happened
- 15:05 the main line launched the dispatched Groups 1-2 command.
- About 15:42 it stopped with `ValueError: Gene axis mismatch across parts: CD4i_R2L02 != CD4i_R2L01`, raised by `run_eval_r16.py::extract_ntc_bundle_from_parts_r16` while extracting cd4t_rest raw controls (round16/results/run_groups12.log).
- The main line recorded R16_MAIN_DEVIATION, changed nothing and waited, as instructed.
- **State at the stop:**
  - round16/results/ held only the log.
  - No prediction, metric or array had been written.
  - The only line processed before the stop was a549_ko: raw controls extracted and SE-600M gene mapping done, all in memory.
  - Nobody has seen any Round 16 result.

## Diagnosis (controller, from the raw files)
- round10/shared/raw/cd4t_rest has 47 lanes: 24 D4 lanes (R2L01-24) and 23 D1 lanes (R1L01-23).
- **Same genes everywhere:**
  - Every lane has exactly 18,130 Gene Expression features.
  - Gene names are unique (0 duplicates).
  - The set of (Ensembl ID, name) pairs is identical in all 47 lanes.
- **Only the order differs.** There are two orders, one an exact permutation of the other:
  - 34 lanes in Ensembl-ID order, starting TSPAN6, TNMD, DPM1;
  - 13 lanes in genomic order, starting SAMD11, NOC2L, KLHL17: R2L01, R2L10, R2L13-R2L23.
- **Cause:**
  - The part builder puts D4 lanes first, so parts[0] = CD4i_R2L01, which is in genomic order. R2L02 is in ID order.
  - The frozen main-line code required an identical ordered gene list. Two other frozen code paths instead align lanes by name onto parts[0]:
    - round10/shared/scripts/prepare_r10_lines.py, which built the cd4t_rest line data;
    - round16/audit/audit_r16_raw_controls.py, A-T3's runner.
- **Other inputs checked:**
  - dld1: 24 raw files, one gene order.
  - feng_ipsc: one CSV.
  - nadig_hepg2 and nadig_jurkat: one h5ad each.
  - a549_ko: passed the check in the run.
  - cd4t_rest is the only affected line.

## Change (round16/src/run_eval_r16.py only; 23 insertions, 4 deletions)
- **Rule:** in `extract_ntc_bundle_from_parts_r16`, a part whose gene list is the same set of unique names as parts[0], in a different order, gets a column permutation onto the parts[0] order. Every other mismatch still raises.
- **Where it's applied:** `_read_part_rows_by_indices` applies the permutation to each part's rows, for the primary sample and the resample alike.
- **Record:** the reordered lane tags are written to `ctrl_audit["deviation_r16_d1_reordered_part_tags"]`.
- **Canonical gene axis:** parts[0] = CD4i_R2L01, the same axis Round 10 and A-T3 use.
- **File hash:** run_eval_r16.py went from cc73f5bc...71fb5 (frozen 848bdac) to b158462f...0581d.

## Why nothing preregistered changes
- No hypothesis, metric, seed, threshold, cohort, data file, model, or adjudication rule changes.
- **Which cells are sampled:** sampling runs over cells (ntc_refs) before any gene column is read, so the sampled barcodes for every line are exactly as before.
- **What the change touches:** only the column order of the raw-count rows from 13 cd4t_rest lanes. The counts of every cell are unchanged, just reindexed by gene name. Lines whose lanes already agree run through code that gives bitwise identical output.
- **A-T3:** unaffected. Its own code already aligns lanes the same way and was not changed.

## Verification
- controller/checks/r16_d1_reorder_test.py, run against the frozen 848bdac version loaded from git. R16_D1_TEST_PASS on all four cases:
  - A: identical order gives bitwise identical output.
  - B: permuted parts reproduce the unpermuted reference exactly.
  - C: a different gene set still raises.
  - D: a duplicated name still raises.
- The main line's own unit-test suite (round16/src/test_unit_r16.py) was run against the patched file: `Ran 18 tests in 109.306s`, `OK`.
- End to end, the planned cross-party comparisons check the result: barcodes exact, SE-600M embeddings within 1e-5, predicted delta within 1e-4.

## Re-freeze and re-run
- FREEZE_MANIFEST_R16_D1.txt = FREEZE_MANIFEST_R16 with the run_eval_r16.py hash replaced, plus this file and the D1 test. FREEZE_MANIFEST_R16.txt is kept unchanged as the historical record.
- R16_GO.txt stays in force with unchanged scope. The code reference becomes FREEZE_MANIFEST_R16_D1.
- The main line reruns Groups 1-2 from scratch with the identical command, then the sensitivity cohort.
