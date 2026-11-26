# DEVIATION_R16_D8 — A-T3 STATE_HVG normalization on the Nadig lines used obs UMI_count instead of the sum over the file's features

Written by Claude (master controller), 2026-11-27 02:50 +0800.

## How it was found
- A-T3 finished pass 1 at 02:46 (`R16_AT3_REAL_DONE`).
- The controller ran staged verification step 2 (`controller/checks/r16_step23.sh` → `r16_compare_arrays.py`) on each A-T3 file `round16/audit/results/{model}_{ckpt}__{line}.npz` against `round16/results/arrays/<line>/`.

| Model | Result |
|---|---|
| STATE_SE (primary, fixed_index_zero) | All 5 lines: predicted Δ for matched, every mismatched source and resample is **bitwise identical** (max \|diff\| 0). Observed Δ, halves, half_valid, transfer, grand mean and FCT are within 1e-5. Targets, genes and target indices are exactly equal. |
| STATE_HVG (secondary, official_median_log1p) | Every arm fed by **Nadig** control cells differs by 3.7e-4 to 2.8e-3 (tolerance 1e-4): the nadig_hepg2 and nadig_jurkat matched and resample arms, and every mismatched arm whose source is nadig_hepg2 or nadig_jurkat. Every arm fed by a549_ko, dld1, cd4t_rest or feng_ipsc controls agrees. On nadig_hepg2 matched, the correlation between the parties is 0.9999998 and the slope 0.99991. |

- Step 3 (controller `r16_check.py` on A-T3's STATE_SE arrays vs the main line's `summary_r16.json`): **R16_MAIN_VS_CONTROLLER_MATCH**, worst \|diff\| 5.6e-16. The primary STATE_SE result therefore agrees across all three parties.

## Cause
- `round16/audit/audit_r16_raw_controls.py`, Nadig branch (`_slice_nadig`): `sel_umis = row_sum_umis if nt_umis is None else nt_umis[sel_local]`. `nt_umis` is `obs["UMI_count"]` (line 265). So A-T3's `cell_total_umis` and `median_ntc_umi` for Nadig pools come from the dataset's UMI_count column.
- What the protocol locks:
  - PREREG_I §2.1 (line 83): `sc.pp.normalize_total` "on the full measured transcriptome of the sampled control pool".
  - Addendum A (line 59): normalization factors "computed over ALL raw Gene Expression features in the sampled control pool".
  - `normalize_total` sums the features in X: 9,624 (HepG2) and 8,882 (Jurkat).
- **Main line** (`round16/src/model_adapters_r16.py::preprocess_hvg_counts_r16`): uses the row sums over the file's features. It is compliant.
- **Controller data check:** for 40 sampled HepG2 controls, UMI_count / row sum is 1.0067–1.0141 (median 1.0091). Median UMI_count is 15,586.5 against a median row sum of 15,436. This per-cell scale difference matches the size and the Nadig-only pattern of the HVG gaps.
- The 10x lines have no UMI_count path (`nt_umis is None`), so they were not affected.

## Effect before the fix (descriptive; secondary model only)
- The A-T3 STATE_HVG metrics on Nadig-fed arms differ from the main line by up to about 5e-5 (tolerance 1e-5). Values are withheld here until the verdict, to keep A-T3 blind.
- All STATE_HVG line calls agree: H-I1 on all 5 lines, and the direction of STATE_HVG − FCT.
- STATE_SE (primary) is unaffected.

## Ruling and change
- 02:50: HOLD A-T3 pass 2 (pid 38323) and send the ruling (`controller/msgs/r16_D8_AT3.txt`).
- One-line fix in A-T3's own code: in the Nadig branch, `sel_umis = row_sum_umis`. QC, barcodes, sampling and gene mapping are unchanged.
- A-T3 then runs tests and a self-check that the totals equal the row sums and the barcodes are unchanged, posts `R16_AT3_D8_PATCHED`, reruns pass 1 for STATE_HVG only into `round16/audit/results/hvg_d8/` (`R16_AT3_HVG_D8_DONE`), and then runs pass 2 from scratch (`R16_AT3_SECONDARY_DONE`).

## What does not change
- The preregistration, Addendum A, the data, the models, the main line, and A-T3's STATE_SE outputs.
- No metric was involved in finding the cause. A-T3 was shown no main-line number.
