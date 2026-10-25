# Source and transformation provenance

The score archive URL, ETag, server length and Last-Modified header are in `download_manifest.json`. Each of its 217 members records the compressed/uncompressed lengths, CRC-32, original CSV SHA-256 and derived numeric-cache SHA-256. The final run reused cached numeric inputs: `range_bytes_this_run` and elapsed time describe that run, not the complete earlier acquisition cost. No multi-GB cache or original sequence-bearing CSV is shipped here.

Official source commit: `144fe22b07dfaeec2b366f2346203a9838a55b4c`.

- [Evaluator](https://github.com/OATML-Markslab/ProteinGym/blob/144fe22b07dfaeec2b366f2346203a9838a55b4c/proteingym/performance_DMS_benchmarks.py).
- [Configuration](https://github.com/OATML-Markslab/ProteinGym/blob/144fe22b07dfaeec2b366f2346203a9838a55b4c/config.json).
- [Constants](https://github.com/OATML-Markslab/ProteinGym/blob/144fe22b07dfaeec2b366f2346203a9838a55b4c/proteingym/constants.json).
- [Official assay-level NDCG table](https://github.com/OATML-Markslab/ProteinGym/blob/144fe22b07dfaeec2b366f2346203a9838a55b4c/benchmarks/DMS_zero_shot/substitutions/NDCG/DMS_substitutions_NDCG_DMS_level.csv).

`assay_groups.csv` is a projection of the official DMS substitution reference metadata, retaining only DMS_id, UniProt_ID and coarse_selection_type for hierarchical aggregation. It excludes target sequences and assay outcome values. The GAL4 fixture projection is described separately in `fixtures/manifest.json`. Hashes of shipped files are collected in `file_hashes.json`.

`official_ndcg_assay.csv` preserves the public rounded reference table; the audit's unrounded metrics are recomputed from released prediction files. These are distinct artifacts and should not be treated as identical precision measurements.

`priority_check.json` records a bounded issue/PR search and prior sources. It does not clear priority comprehensively or constitute author endorsement.
