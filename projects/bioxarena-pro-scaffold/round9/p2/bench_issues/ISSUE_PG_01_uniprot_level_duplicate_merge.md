# [Bug] `uniprot_function_lookup` missing `.drop_duplicates()` and post-merge `Average` inflate and reweight `DMS_substitutions_*_Uniprot_level.csv`

**Commit checked:** `144fe22b07dfaeec2b366f2346203a9838a55b4c` (`main`)

## Summary

In `proteingym/performance_DMS_benchmarks.py`, line 136 omits `.drop_duplicates()` when building `uniprot_function_lookup` (unlike the adjacent lookups on lines 138 and 140), and `uniprot_level_average` (line 305) is computed **after** left-joining three metadata tables onto `uniprot_metric_performance` (lines 300–302).

Because `uniprot_function_lookup` retains all 217 rows—and `uniprot_Neff_lookup` retains 188 rows because two proteins (`DYR_ECOLI` and `SPG1_STRSG`) have two `MSA_Neff_L_category` values across assays—the merges expand the 186 unique `UniProt_ID`s to **221 data rows** (+ 1 `Average` row = 222 rows) in all five published `DMS_substitutions_*_Uniprot_level.csv` files, weighting multi-assay proteins up to **4×** (`BLAT_ECOLX`, `P53_HUMAN`, `DYR_ECOLI`, `SPG1_STRSG`) in the `Average` row.

**Scope note:** This does **not** affect the main leaderboard `Summary_performance_DMS_substitutions_*.csv` (`Average_<metric>` / `final_average`), because line 298 groups directly by `['UniProt_ID', 'Selection Type']` without using `uniprot_function_lookup`. It also does not affect `--indel_mode` outputs (no `_Uniprot_level.csv` is written).

## Minimal Reproduction

```python
import json
import pandas as pd

ref = pd.read_csv("reference_files/DMS_substitutions.csv")
with open("proteingym/constants.json") as f:
    clean_names = json.load(f)["clean_names"]
rev_clean = {v: k for k, v in clean_names.items()}

# Lines 136-141 of proteingym/performance_DMS_benchmarks.py
print("Unique UniProt_IDs:", ref["UniProt_ID"].nunique())  # 186
print("uniprot_function_lookup (missing .drop_duplicates()):", len(ref[["UniProt_ID", "coarse_selection_type"]]))  # 217
print("uniprot_Neff_lookup (.drop_duplicates()):", len(ref[["UniProt_ID", "MSA_Neff_L_category"]].drop_duplicates()))  # 188 (DYR_ECOLI, SPG1_STRSG)

for metric in ["Spearman", "AUC", "MCC", "NDCG", "Top_recall"]:
    dms = pd.read_csv(f"benchmarks/DMS_zero_shot/substitutions/{metric}/DMS_substitutions_{metric}_DMS_level.csv")
    uni = pd.read_csv(f"benchmarks/DMS_zero_shot/substitutions/{metric}/DMS_substitutions_{metric}_Uniprot_level.csv")
    meta = {"UniProt_ID", "MSA_Neff_L_category", "Taxon", "Selection Type"}
    models = [c for c in uni.columns if c not in meta]

    pub_data = uni[uni["UniProt_ID"].notna()]
    pub_mean = pub_data[models].mean(axis=0)
    true_186_mean = dms.rename(columns=rev_clean).groupby("UniProt_ID")[models].mean(numeric_only=True).mean(axis=0)

    r_pub = pub_mean.rank(ascending=False, method="min")
    r_true = true_186_mean.rank(ascending=False, method="min")
    print(
        f"{metric:10s}: rows={len(pub_data)}->186 | "
        f"rank_changes={int((r_pub != r_true).sum())}/{len(models)} | "
        f"ESM3={pub_mean['ESM3']:.3f}(#{int(r_pub['ESM3'])})->{true_186_mean['ESM3']:.3f}(#{int(r_true['ESM3'])})"
    )
```

Observed output on `main` (`144fe22b`):
- `Spearman`: `221 -> 186` rows, **49/97** models change rank (`58/97` if comparing rounded 3-decimal `Average` rows); `ESM3` moves `0.480 (#15) -> 0.493 (#10)`.
- `AUC`: **53/97** (`59/97` rounded) models change rank; `ESM3` moves `0.766 (#16) -> 0.771 (#11)`.
- `MCC`: **52/97** (`58/97` rounded) models change rank; `ESM3` moves `0.380 (#14) -> 0.392 (#11)`.
- `NDCG`: **66/97** (`71/97` rounded) models change rank.
- `Top_recall`: **61/97** (`63/97` rounded) models change rank; `ESM3` moves `0.247 (#4) -> 0.255 (#3)`.

## Suggested Fix

Note that adding `.drop_duplicates()` at line 136 alone still leaves **202 rows** (due to 14 dual-function UniProt IDs plus `DYR_ECOLI` and `SPG1_STRSG` having two `MSA_Neff_L_category` values). Computing `uniprot_level_average` before the metadata merges and collapsing metadata per `UniProt_ID` keeps `_Uniprot_level.csv` at one row per unique `UniProt_ID` (186 + `Average`):

```diff
--- a/proteingym/performance_DMS_benchmarks.py
+++ b/proteingym/performance_DMS_benchmarks.py
@@ -133,12 +133,13 @@ def main():
         config = json.load(f)
     with open(f"{os.path.dirname(os.path.realpath(__file__))}/constants.json") as f:
         constants = json.load(f)
-    uniprot_function_lookup = mapping_protein_seq_DMS[["UniProt_ID","coarse_selection_type"]]
+    def _collapse_meta(df, col):
+        return df[["UniProt_ID", col]].drop_duplicates().groupby("UniProt_ID", as_index=False)[col].agg("; ".join)
+
+    uniprot_function_lookup = _collapse_meta(mapping_protein_seq_DMS, "coarse_selection_type")
     uniprot_function_lookup.columns = ["UniProt_ID", "Selection Type"]
-    uniprot_Neff_lookup = mapping_protein_seq_DMS[['UniProt_ID','MSA_Neff_L_category']].drop_duplicates()
+    uniprot_Neff_lookup = _collapse_meta(mapping_protein_seq_DMS, 'MSA_Neff_L_category')
     uniprot_Neff_lookup.columns=['UniProt_ID','MSA_Neff_L_category']
-    uniprot_taxon_lookup = mapping_protein_seq_DMS[['UniProt_ID','taxon']].drop_duplicates()
+    uniprot_taxon_lookup = _collapse_meta(mapping_protein_seq_DMS, 'taxon')
     uniprot_taxon_lookup.columns=['UniProt_ID','Taxon']
     if args.indel_mode:
         args.performance_by_depth = False
@@ -297,12 +298,12 @@ def main():
             uniprot_metric_performance = performance_all_DMS[metric].groupby(['UniProt_ID']).mean(numeric_only=True)
             uniprot_function_metric_performance = performance_all_DMS[metric].groupby(['UniProt_ID',"Selection Type"]).mean(numeric_only=True)
+            uniprot_level_average = uniprot_metric_performance.drop(columns=['number_mutants']).mean(numeric_only=True)
             uniprot_metric_performance = uniprot_metric_performance.reset_index()
             uniprot_metric_performance = pd.merge(uniprot_metric_performance,uniprot_Neff_lookup,on='UniProt_ID', how='left')
             uniprot_metric_performance = pd.merge(uniprot_metric_performance,uniprot_taxon_lookup,on='UniProt_ID', how='left')
             uniprot_metric_performance = pd.merge(uniprot_metric_performance,uniprot_function_lookup,on="UniProt_ID",how="left")
             del uniprot_metric_performance['number_mutants']
             del uniprot_function_metric_performance["number_mutants"]
-            uniprot_level_average = uniprot_metric_performance.mean(numeric_only=True)
             uniprot_function_level_average = uniprot_function_metric_performance.groupby("Selection Type").mean(numeric_only=True)
```

Regenerating `benchmarks/DMS_zero_shot/substitutions/{Spearman,AUC,MCC,NDCG,Top_recall}/DMS_substitutions_*_Uniprot_level.csv` will then yield 186 protein rows + 1 `Average` row without changing `Summary_performance_*.csv`.
