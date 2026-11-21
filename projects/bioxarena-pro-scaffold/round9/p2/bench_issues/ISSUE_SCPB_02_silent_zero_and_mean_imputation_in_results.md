# [Bug] Missing method–dataset runs are silently zero-filled (or mean-imputed) in `Results/*_performance_top{100,5000}.csv`

**Commit checked:** `63fd9b95071b894c4f39c8f7eb2bd42ca9a70a77` (`main`)

## Summary

Cross-checking the raw per-perturbation metric files (`all_dataset_genetic_seed1-3.csv`, `all_dataset_chemical.csv`, `cellular_iid_distance0.csv`) against the aggregated `*_performance_top100.csv` and `*_performance_top5000.csv` tables reveals that missing runs are silently filled with `0.0000` (when an entire `(DataSet, seed)` is absent) or constant mean values (when partially absent), rather than marked `NaN` or excluded from dataset averages:

1. **Genetic Single (`scFoundation`):** In `all_dataset_genetic_seed1-3.csv`, `scFoundation` is missing `Replogle_K562essential` (323×3 = 969), `Replogle_RPE1essential` (266×3 = 798), `Replogle_exp8` seeds 2–3 (24), and `Replogle_exp7` seed 3 (20) — **1,811 / 2,181 (83.0%)** single-perturbation evaluations. In `genetic_single_performance_top100.csv` and `top5000.csv`, all 1,811 missing rows are filled with exact `0.0000` across all metrics, dragging `scFoundation`'s `ac_score` from `0.3202` (top100) / `0.2355` (top5000) on the 7 covered datasets down to `0.2491` / `0.1832` (last overall).
2. **Chemical Single (`cycleCDR`, `PRnet`, `biolord`, `chemCPA`):**
   - `cycleCDR` has **0/152** rows for `sciplex3_K562` seed 3 in `all_dataset_chemical.csv`, yet appears in `chemical_single_performance_top100.csv` and `top5000.csv` with **152 all-zero rows** (reducing `cycleCDR` on `K562` to `0.3408` vs `0.5951` on `A549` and `0.6227` on `MCF7`).
   - Partially missing `(DataSet, seed)` runs are filled with identical repeated values in `chemical_single_performance_top100/5000.csv`: **120** repeated rows per cell line for `PRnet` seed 3 (see separate `myPRnet.py` bug report), and **2** repeated rows each for `biolord` and `chemCPA` in `sciplex3_MCF7` seed 1 (34/36 present in `all_dataset_chemical.csv`).
3. **Cellular Context IID (`biolord`):** In `Results/Cellular_context_iid/cellular_iid_distance0.csv`, `bioLord` **does** have evaluated scores on all 12 datasets—including `Afriat` (`pearson_distance = 0.9908`), `McFarland` (`0.6802`), and `crossSpecies` (**#1** among all 14 methods at `0.9965`). Yet in `cellular_iid_performance_top100.csv` and `top5000.csv`, `biolord` has `0.0000` for `cor`, `edistance_score`, `mse_score`, `was_score`, and `sym_score` on those 3 datasets (only `DEG_score` is non-zero in `top100`), dropping `biolord` from **#3** (`0.7177` top100, `0.8079` top5000 across the 9 non-zeroed datasets) to **#9** (`0.5389` top100, `0.6060` top5000 across 12 datasets).
4. **Related discrepancy (`pearson_distance` clipping):** `all_dataset_genetic_seed1-3.csv` (`5,585/11,340 = 49.3%` of `baseControl` rows `== 0.0`) and `all_dataset_chemical.csv` (`2,212/4,268 = 51.8%` `== 0.0`) clip negative PCC-delta at `0.0`, whereas `*_performance_top100/5000.csv` retain signed negative `cor` (`1,061/2,181 = 48.6%` < 0, min `-0.5707`).

## Minimal Reproduction

```python
import pandas as pd

# 1. scFoundation zero-filled rows in genetic_single_performance_top100.csv
g100 = pd.read_csv("Results/Genetic_perturbation/genetic_single_performance_top100.csv")
sf_zero = g100[(g100["method"] == "scFoundation") & (g100["ac_score"] == 0)]
print("scFoundation zero rows:", len(sf_zero), sf_zero.groupby("DataSet").size().to_dict())

# 2. cycleCDR zero-fill & PRnet/biolord/chemCPA constant-fill in chemical_single_performance_top100.csv
c100 = pd.read_csv("Results/Chemical_perturbation/chemical_single_performance_top100.csv")
cyc_zero = c100[(c100["method"] == "cycleCDR") & (c100[["cor", "edistance_score", "mse_score", "ac_score"]] == 0).all(axis=1)]
print("cycleCDR all-zero rows:", len(cyc_zero), cyc_zero["DataSet"].value_counts().to_dict())
for m in ["PRnet", "biolord", "chemCPA"]:
    vc = c100[c100["method"] == m].groupby(["DataSet", "cor", "edistance_score", "mse_score"]).size()
    print(m, "repeated rows:", vc[vc > 1].to_dict())

# 3. biolord present in cellular_iid_distance0.csv yet zeroed in cellular_iid_performance_top100.csv
d_iid = pd.read_csv("Results/Cellular_context_iid/cellular_iid_distance0.csv")
p100 = pd.read_csv("Results/Cellular_context_iid/cellular_iid_performance_top100.csv")
for ds in ["Afriat", "McFarland", "crossSpecies"]:
    d_pcc = d_iid[(d_iid["method"] == "bioLord") & (d_iid["DataSet"] == ds) & (d_iid["metric"] == "pearson_distance") & (d_iid["DEG"] == 100)]["performance"].mean()
    p_cor = p100[(p100["method"] == "biolord") & (p100["DataSet"] == ds)]["cor"].mean()
    print(f"Cellular IID {ds}: distance0 pcc={d_pcc:.4f} vs top100 cor={p_cor:.4f}")
```

## Suggested Fix

1. Do **not** impute missing runs with `0.0` or constant means when generating `*_performance_top{100,5000}.csv`; represent missing evaluations as `NaN` and report per-method dataset/perturbation coverage alongside matched-subset rankings.
2. Re-aggregate `cellular_iid_performance_top{100,5000}.csv` for `biolord` (`bioLord` in `cellular_iid_distance0.csv`) so its existing `Afriat`, `McFarland`, and `crossSpecies` evaluations are not zeroed out (and note `Cellular_context_generalization/i.i.d./` is also missing `mybiolord.py`).
3. Keep `pearson_distance` / `cor` signed consistently in `all_dataset_*.csv` or document the `0.0` floor.
