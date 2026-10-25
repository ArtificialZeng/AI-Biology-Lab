# Protein Metric Auditor

A reproducible, numeric-only audit of row-order dependence in ProteinGym's **secondary NDCG@10% metric**. This project studies the frozen public v1.3 zero-shot substitution score archive; it does not run protein models or design sequences.

The official evaluator assigns distinct ordinal ranks to tied predictions using NumPy sorting. Reordering otherwise identical aligned records can therefore change NDCG. This handling was already described in [VenusREM-Harness, Appendix B.1 (2026)](https://arxiv.org/html/2609.34654v1#A2.SS1). Our contribution is a released-score census and reproducible witnesses, not discovery of tie handling or a new ranking formula. Priority is not established.

## Results and scope

| Frozen census quantity | Result |
| --- | ---: |
| Assays / model configurations | 217 / 95 |
| Rows across distinct assay tables | 2,465,767 |
| Finite assay/configuration cells | 20,615 |
| Cells with a nonzero tie-order envelope, threshold `1e-10` | 5,006 |
| Cells changed by original/reverse/seed-1 common row orders | 4,656 |
| Largest assay/configuration change among those three orders | 0.0422548 NDCG |
| Largest aggregate rank movement under the two common reorderings | 2 positions |
| Largest aggregate change after standard tie averaging | 0.000342989 NDCG |

The top five configurations are unchanged under reversal, seed-1 shuffle and tie averaging. No change to the primary **Spearman** ranking has been established. These are descriptive results for a fixed released archive, without a significance, new-protein generalization or biological-mechanism claim. This is not an audit of every 2026 model: `AIDO.Protein-RAG-16B` and `Protriever` occur in the pinned configuration but are absent from the score archive.

Within each assay, reversal and shuffle apply the **same row permutation to every model**, keeping labels and predictions aligned. In contrast, each model's minimum/maximum over ties can require a different permutation. Those separate envelopes—and the `tie_min_rank`/`tie_max_rank` columns in the archived aggregate table—**do not form one jointly attainable leaderboard**.

Aggregation follows the official hierarchy: mean assays within each UniProt/function cell, mean proteins within each function category, then equal-weight mean over categories. Plain assay averaging answers a different question.

## Offline reproduction

Use Python 3.11–3.13 and an isolated environment. Dependencies are pinned to the versions used for the census; unstable ordinal tie order may differ across NumPy versions or hardware. The archived census was generated on macOS arm64, Python 3.13.5, NumPy 2.2.3 and pandas 2.2.3.

```sh
python -m venv /path/to/local-venv
/path/to/local-venv/bin/python -m pip install -r requirements.txt
/path/to/local-venv/bin/python -B scripts/reproduce_sample.py
/path/to/local-venv/bin/python -B scripts/verify_results.py
```

Run these commands from this project directory. The sample contains 1,195 GAL4 yeast assay rows and six numeric columns: measured score plus ESCOTT, GEMME, ESMC-300M, ESMC-600M and VenusREM predictions. It contains **no sequence or mutation columns**. The script extracts the real pinned official `calc_ndcg` function via AST, compares original/reverse/seed-1 results with an independent implementation, checks standard tie-averaged invariance, and strictly matches the archived order-free statistics. Original/reverse/shuffle ordinal values are compared with the archive and every difference is reported; `--strict-reference` additionally makes those differences an error. This distinction permits a correct local evaluator on another platform without concealing differences in ordinal output. It rejects nonfinite or constant experimental values and invalid normalized gains; nonfinite model predictions are explicitly marked rather than silently dropped.

Outputs are written to the repository-root `artifacts/protein-metric-auditor/`, never into this tracked project. `verify_results.py` checks file hashes, census coverage and the hierarchical aggregates from the shipped cell table; it does not regenerate all predictions.

## Optional complete numeric replay

The complete source score ZIP is about 1.91 GB compressed; its [official URL](https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/zero_shot_substitutions_scores.zip), ETag, length, and all 217 member CRC/SHA-256 values are preserved in [the manifest](provenance/download_manifest.json). No download occurs automatically. Given that local ZIP, run:

```sh
python -B scripts/replay_census.py --score-zip /path/to/zero_shot_substitutions_scores.zip
```

Alternatively, use `--numeric-cache /path/to/numeric_cache` with the original lossless cache matching the recorded hashes. The multi-GB cache is intentionally omitted. `--assay GAL4_YEAST_Kitzman_2015` performs a limited replay. Only configured numeric columns are retained in memory and output; source sequences are not written. The complete replay is a strict archived-reference comparison, including ordinal outputs; its exact tie-order match may require the recorded platform and versions. A full replay reads all assays in the frozen release but performs no inference, training, optimization of protein function, or clinical interpretation.

## Files and provenance

- [Census summary](results/census_summary.json), [cell metrics](results/cell_metrics.csv), [aggregate metrics](results/aggregate_metrics.csv).
- [Numeric sample and transformation manifest](fixtures/manifest.json).
- [Pinned official evaluator](vendor/proteingym/performance_DMS_benchmarks.py), commit `144fe22b07dfaeec2b366f2346203a9838a55b4c`.
- [Third-party credit and licenses](THIRD_PARTY.md); [Chinese brief](docs/BRIEF.zh-CN.md).
- [Archived independent formula check](validation/archived_formula_review.json). This is a separate computational/AI-assisted check of the recorded earlier formula implementation, not external human peer review or certification of every line in this packaged project.

Standard tie averaging and extremal rearrangement bounds are established methods. The [September 2026 RewireBio audit](https://rewirebio.io/blog/what-our-biology-benchmark-scores-actually-counted/) already examined a different ProteinGym constant-prediction Top-recall issue; it is not presented here as a new result. No upstream issue or PR has been submitted, no maintainer endorsement obtained, and no media outcome is promised.
