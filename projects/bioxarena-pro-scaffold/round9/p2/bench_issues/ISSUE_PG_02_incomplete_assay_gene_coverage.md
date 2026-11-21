# [Benchmark / Reporting] Unequal assay/gene coverage in DMS (`Protriever` 200/217) and Clinical (`1,702–2,525/2,525`) summaries without coverage columns or matched-subset reporting

**Commit checked:** `144fe22b07dfaeec2b366f2346203a9838a55b4c` (`main`)

## Summary

Both `proteingym/performance_DMS_benchmarks.py` and `proteingym/performance_clinical_benchmarks.py` compute leaderboard averages (`mean()`) over each model's available assays/genes without reporting assay/gene coverage in the checked-in summary tables or providing a matched complete-case comparison:

1. **DMS Zero-Shot Substitutions (`Protriever`):** Among all 97 models, `Protriever` is missing **17/217** assays (**12** unique `UniProt_ID`s). Those 17 assays are substantially larger (`mean = 41,639.2` vs `8,789.5` mutants) and lower-scoring for complete-coverage models (`VenusREM`: `0.4648` on the 17 vs `0.5418` on the 200; `ProSST (K=512)`: `0.3590` vs `0.4998`). On the common **200-assay intersection**:
   - Spearman: `ESCOTT` rises `0.4759 (#9) -> 0.4812 (#8)` and `ProSST (K=512)` rises `0.4712 (#10) -> 0.4804 (#9)`, dropping `Protriever` (`0.4789`) from **#8 to #10**.
   - AUC: `ESCOTT` rises `0.7613 (#9) -> 0.7630 (#8)` ahead of `Protriever` (`0.7617`, `#8 -> #9`).
   - Top_recall: `Protriever` drops from **#16 to #20**.
2. **Clinical Substitutions (Zero-Shot & Supervised):**
   - **Schema mismatch in checked-in zero-shot summary:** `benchmarks/clinical_zero_shot/substitutions/AUC/Summary_performance_clinical_substitutions_AUC.csv` has a stray empty `Model_rank.1` column, lacks the `Num_genes` column that `performance_clinical_benchmarks.py:188` generates, and `config.json` lists `SiteRM` even though `SiteRM` is absent from `zero_shot_clinical_substitutions_scores.zip`.
   - **Unequal gene coverage:** 8/13 zero-shot models miss 76–790 genes (`MutPred`: `1,735/2,525`; `LIST-S2`: `1,893/2,525`), and all 17 supervised models miss 25–823 genes (`MPC`: `1,702/2,525`; `gMVP`: `1,772/2,525`). Missing genes are systematically smaller (`MutPred` missing: median 7 variants vs 15 present).
   - **Rank shifts on complete cases:** On the zero-shot intersection (`N = 1,129`), `MutPred` (`#10 -> #9`) swaps with `MutationAssessor` (`#9 -> #10`), and point estimates for `PoET (200M)` (`0.9200 #1 -> 0.9168 #3`) and `GEMME` (`0.9191 #3 -> 0.9192 #1`) reverse—though `GEMME - PoET = +0.0024` is within statistical noise (`95% bootstrap CI [-0.0039, +0.0088]`, Wilcoxon `p = 0.685`). On the supervised intersection (`N = 1,374`), `VARITY (R)` vs `BayesDel (noAF)` swap `#6/#7` and `PolyPhen2 (HDIV)` vs `MPC` swap `#13/#14`.
   - **Tiny-gene equal weighting & pooling sensitivity:** `883/2,525` (`35.0%`) clinical genes have `min(n_pos, n_neg) == 1`. Weighting by variant count or filtering `min_cls >= 5` (`N = 717`) moves `MutPred` from `0.8754 (#10)` to `0.8889` / `0.9051`; in pooled variant-level AUC (`62,727` variants), `MutPred` is #1 (`0.9173`) while `TranceptEVE L` drops from `#2` (`0.9198`) to `#10` (`0.8757`) due to uncalibrated per-gene score scales (`0.7772` after per-gene $z$-scoring).

## Minimal Reproduction

```python
import pandas as pd

# 1. DMS Protriever 200-assay intersection
dms_sp = pd.read_csv("benchmarks/DMS_zero_shot/substitutions/Spearman/DMS_substitutions_Spearman_DMS_level.csv")
meta = {"DMS ID", "Number of Mutants", "Selection Type", "UniProt ID", "MSA_Neff_L_category", "Taxon"}
models = [c for c in dms_sp.columns if c not in meta]

def cat_avg(df):
    return df.groupby(["UniProt ID", "Selection Type"])[models].mean(numeric_only=True).groupby("Selection Type").mean().mean()

full, c200 = cat_avg(dms_sp), cat_avg(dms_sp[dms_sp["Protriever"].notna()])
for m in ["Protriever", "ESCOTT", "ProSST (K=512)"]:
    print(m, "full:", round(full[m], 4), int(full.rank(ascending=False)[m]),
          "c200:", round(c200[m], 4), int(c200.rank(ascending=False)[m]))

# 2. Clinical complete-case intersections
for track in ["clinical_zero_shot", "clinical_supervised"]:
    df = pd.read_csv(f"benchmarks/{track}/substitutions/AUC/clinical_substitutions_AUC_DMS_level.csv")
    mc = [c for c in df.columns if c != "RefSeq ID"]
    cc = df.dropna(subset=mc)
    rf, rc = df[mc].mean().rank(ascending=False).astype(int), cc[mc].mean().rank(ascending=False).astype(int)
    print(track, f"complete_N={len(cc)}/{len(df)} changed:", [m for m in mc if rf[m] != rc[m]])
```

## Suggested Improvements

1. Regenerate `benchmarks/clinical_zero_shot/substitutions/AUC/Summary_performance_clinical_substitutions_AUC.csv` with `performance_clinical_benchmarks.py` (removing `Model_rank.1` and exposing `Num_genes`), and align `config.json` (`SiteRM`) with `zero_shot_clinical_substitutions_scores.zip`.
2. Add a `Num_assays` column to `Summary_performance_DMS_*.csv` (and `Num_genes` to clinical supervised/zero-shot summaries) so readers immediately see incomplete coverage (`200/217`, `1,735/2,525`, `1,702/2,525`).
3. Consider adding an optional `--complete_cases_only` flag in `performance_DMS_benchmarks.py` and `performance_clinical_benchmarks.py` (plus reporting `min_cls >= 5` / variant-weighted AUC sensitivity in clinical substitutions).
