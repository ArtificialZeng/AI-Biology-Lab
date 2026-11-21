# [Bug] Leaked `CPA` entries remain in all four `Results/Genetic_perturbation/*_performance_*.csv` files despite README notice

**Commit checked:** `63fd9b95071b894c4f39c8f7eb2bd42ca9a70a77` (`main`)

## Summary

The README notice (*"Notice on CPA results in the genetic perturbation benchmark"*, `README.md:103-107`) states:

> *"We have removed CPA entries from the latest four genetic perturbation performance summary files."*

While `CPA` was removed from `Results/Genetic_perturbation/all_dataset_genetic_seed1.csv`, `seed2.csv`, and `seed3.csv` (0 rows each), **all four summary files on `main` still contain `CPA`**:

- `Results/Genetic_perturbation/genetic_single_performance_top100.csv`: **2,181 `CPA` rows** — `CPA` still ranks **#1 overall** by both dataset mean (`ac_score = 0.5508`, ahead of `scGPT` `0.5356` and `scouter` `0.5207`) and perturbation mean (`0.6642`, ahead of `GenePert` `0.5098`).
- `Results/Genetic_perturbation/genetic_single_performance_top5000.csv`: **2,181 `CPA` rows** (`ac_score = 0.5787` dataset mean, #9; `0.6855` perturbation mean, #5).
- `Results/Genetic_perturbation/genetic_combo_performance_top100.csv`: **654 `CPA` rows** (`ac_score = 0.4962` dataset mean, #11).
- `Results/Genetic_perturbation/genetic_combo_performance_top5000.csv`: **654 `CPA` rows** (`ac_score = 0.5871` dataset mean, #9).

In addition, because per-perturbation normalized scores (`edistance_score`, `mse_score`, `was_score`, `sym_score`) and ranks (`Rank_*`, `ac_rank`) in those tables were computed **with `CPA` in the min-max pool**, `CPA` sets the `1.0` upper bound on **1,254/2,181** (`mse_score`) and **1,548/2,181** (`was_score`) rows in `genetic_single_performance_top100.csv` (and **1,638/2,181** `mse_score` rows in `top5000`). Simply dropping `CPA` rows changes the top rankings (`scGPT` #2 → #1 and `scouter` #3 → #2 on `genetic_single_performance_top100.csv` by dataset mean), and regenerating the four summary files from `all_dataset_genetic_seed1-3.csv` is needed so min-max bounds and `Rank_*`/`ac_rank` are no longer anchored by leaked `CPA`.

## Minimal Reproduction

```python
import pandas as pd

for s in [1, 2, 3]:
    raw = pd.read_csv(f"Results/Genetic_perturbation/all_dataset_genetic_seed{s}.csv")
    print(f"all_dataset_genetic_seed{s}.csv CPA rows:", (raw["method"] == "CPA").sum())

for fn in [
    "genetic_single_performance_top100.csv",
    "genetic_single_performance_top5000.csv",
    "genetic_combo_performance_top100.csv",
    "genetic_combo_performance_top5000.csv",
]:
    df = pd.read_csv(f"Results/Genetic_perturbation/{fn}")
    ds_means = df.groupby(["method", "DataSet"])["ac_score"].mean().groupby("method").mean().sort_values(ascending=False)
    pert_means = df.groupby("method")["ac_score"].mean().sort_values(ascending=False)
    print(
        f"{fn}: CPA_rows={(df['method'] == 'CPA').sum()}, "
        f"ds_mean={ds_means['CPA']:.4f} (#{list(ds_means.index).index('CPA') + 1}), "
        f"pert_mean={pert_means['CPA']:.4f} (#{list(pert_means.index).index('CPA') + 1})"
    )
```

## Suggested Fix

1. Regenerate `genetic_single_performance_top100.csv`, `genetic_single_performance_top5000.csv`, `genetic_combo_performance_top100.csv`, and `genetic_combo_performance_top5000.csv` from `all_dataset_genetic_seed{1,2,3}.csv` (excluding `CPA` *before* min-max normalization, `ac_score`, and `Rank_*`/`ac_rank` calculation), or rerun `CPA` with the corrected counterfactual script (`Perturbation_generalization/Genetic/myCPA.py`) first.
2. Until regenerated, update `README.md` so users know the four checked-in summary CSVs still include the leaked `CPA` rows.
