# Result-table interpretation

`cell_metrics.csv` contains one row for each assay/model configuration; `n` counts measured rows within that assay. It is not a count of globally unique variants. All 20,615 cells in the frozen census are finite and have status `ok`.

`original`, `reversed` and `shuffle_seed_1` use actual common row orders across models within each assay. `tie_average` is standard expected NDCG under uniform random breaking of exact prediction-score ties. `tie_min` and `tie_max` are per-model extrema with all between-score rankings held fixed. A different permutation may be needed for each model's extrema; neither the cell nor aggregate extrema define a simultaneous leaderboard. The accompanying `tie_min_rank` and `tie_max_rank` columns only rank those separate envelopes.

`aggregate_metrics.csv` follows UniProt/function-cell means, then protein means within functional category, then equal category weighting. `census_summary.json` is the full archived descriptive summary. The scope is a secondary NDCG metric in a fixed public release, without a biological conclusion or an established effect on primary Spearman scores.
