# B1 Verification & Issue Triage Summary (ProteinGym & scPerturBench)

**Status:** Drafts only — **not posted**. All reproduction scripts were run locally against snapshots byte-identical to GitHub `main` HEAD on 2026-11-21.

## 1. Repository Snapshots & Existing Issue/PR Deduplication

| Repository | Verified `main` HEAD OID | Existing Issues & PRs Checked | Overlap with Drafts |
| :--- | :--- | :--- | :--- |
| `OATML-Markslab/ProteinGym` | `144fe22b07dfaeec2b366f2346203a9838a55b4c` | All 91 issues (#1–#133: 28 open, 63 closed) and 42 PRs (#7–#130) | **None.** Closest issues (#4/#15 on aggregation intent; #125/PR #128 on supervised MSA/taxon column swaps; #131 on v1.3 merge layout; #132 on cutoff sign/indel taxon; #133 on `calc_ndcg` ties) do **not** report `uniprot_function_lookup` duplication, `Protriever` missing 17/217 assays, or clinical unequal gene coverage. |
| `bm2-lab/scPerturBench` | `63fd9b95071b894c4f39c8f7eb2bd42ca9a70a77` | All 13 issues (#1–#13) and 1 PR (#14) | **None for our drafts.** Issue #1 reported CPA test-cell reconstruction leakage in genetic perturbation (leading to the README disclaimer), #3 asked about `pearson_distance`, #7 fixed identical IID/OOD DEG scores, #13/PR #14 fixed cellular-context `baseReg`/`baseMLP` paired direction. None report CPA still present in all four genetic summary CSVs, silent zero-/mean-imputation, `myPRnet.py` split-3 overwrite, `myCPA_chemical.py` unpatched view/key bug, or `generateExp` vs `edistance` asymmetry. |

## 2. Verification of Task-Card Items & S4 Findings (with Honest Corrections/Caveats)

### A. ProteinGym (`OATML-Markslab/ProteinGym`)

1. **`uniprot_function_lookup` missing `.drop_duplicates()` + post-merge `Average` (`proteingym/performance_DMS_benchmarks.py:136, 300-305`) — CONFIRMED (with scope refinement)**
   - **Repro:** `repro/repro_pg_uniprot_lookup.py`
   - **Exact mechanism:** `uniprot_function_lookup` (line 136) lacks `.drop_duplicates()`, leaving 217 rows instead of 200 (unlike lines 138 and 140 which call `.drop_duplicates()`). In addition, `uniprot_Neff_lookup` (line 138) has 188 rows because `DYR_ECOLI` (`High`/`Medium`) and `SPG1_STRSG` (`Low`/`Medium`) span two `MSA_Neff_L_category` values across assays. Merging both before computing `uniprot_level_average` (line 305) produces a Cartesian 4× expansion on `DYR_ECOLI`/`SPG1_STRSG` and duplicates multi-assay proteins (`BLAT_ECOLX` 4×, `P53_HUMAN` 4×, `HSP82_YEAST` 3×, etc.), inflating `DMS_substitutions_*_Uniprot_level.csv` from **186** unique UniProts to **221 data rows** (+ `Average` = 222). Adding `.drop_duplicates()` on line 136 alone still leaves **202 rows** unless `uniprot_level_average` is computed *before* metadata merges (or metadata lookups are collapsed per `UniProt_ID`).
   - **Quantified impact:**
     - **Affected files:** All five `benchmarks/DMS_zero_shot/substitutions/{Spearman,AUC,MCC,NDCG,Top_recall}/DMS_substitutions_*_Uniprot_level.csv` files. Comparing unrounded means over the published 221-row tables vs the true unweighted 186-UniProt mean changes ranks for **49/97** (Spearman), **53/97** (AUC), **52/97** (MCC), **66/97** (NDCG), and **61/97** (Top_recall) models (58/59/58/71/63 if comparing rounded 3-decimal published `Average` vs rounded 186-UniProt mean). E.g., `ESM3` moves in `_Uniprot_level.csv` Spearman `0.480 (#15) -> 0.493 (#10)`, AUC `0.766 (#16) -> 0.771 (#11)`, MCC `0.380 (#14) -> 0.392 (#11)`, Top_recall `0.247 (#4) -> 0.255 (#3)`.
     - **Unaffected files (explicitly verified):** `Summary_performance_DMS_substitutions_*.csv` (`Average_*` / `final_average`, max diff = `0.000`), because line 298 groups directly by `['UniProt_ID', 'Selection Type']` without using `uniprot_function_lookup`. `--indel_mode` outputs are also unaffected (no `_Uniprot_level.csv` written).
   - **Caveat on S4 (`line 298`):** S4 claimed line 298 double-counts 14 dual-function UniProt IDs in `Summary_performance`. Per ProteinGym's `README.md` and maintainer explanation in Issue #4, grouping by `(UniProt_ID, Selection Type)` before averaging across the 5 functional categories is the **intentional** v1.0+ hierarchical design, not a bug. Excluded from the bug issue draft.
   - **Draft:** `ISSUE_PG_01_uniprot_level_duplicate_merge.md`

2. **Incomplete assay/gene coverage (`Protriever` 17/217 in DMS; 8/13 zero-shot & 17/17 supervised in Clinical) — CONFIRMED (with statistical caveat)**
   - **Repro:** `repro/repro_pg_missing_coverage.py`
   - **DMS Zero-Shot (`Protriever`):** Only `Protriever` has missing DMS substitution assays: **17/217** across **12** unique UniProt IDs (note: S4 wrote 11 UniProts; actual count in `DMS_substitutions_Spearman_DMS_level.csv` is 12). The 17 missing assays are much larger (`mean = 41,639.2` mutants vs `8,789.5` on the 200 evaluated assays) and harder (`VenusREM`: `0.5418` on the 200 vs `0.4648` on the 17; `ProSST (K=512)`: `0.4998` vs `0.3590`). On the matched **200-assay intersection**, `ESCOTT` rises from `0.4759 (#9)` to `0.4812 (#8)` and `ProSST (K=512)` from `0.4712 (#10)` to `0.4804 (#9)`, dropping `Protriever` (`0.4789`) **#8 -> #10** (Spearman), **#8 -> #9** (AUC), and **#16 -> #20** (Top_recall).
   - **Clinical Substitutions:**
     - `benchmarks/clinical_zero_shot/substitutions/AUC/Summary_performance_clinical_substitutions_AUC.csv` contains a stray empty `Model_rank.1` column, omits `Num_genes` (which `performance_clinical_benchmarks.py:188` outputs), and `config.json` lists `SiteRM` which is absent from `zero_shot_clinical_substitutions_scores.zip`.
     - Coverage ranges from `1,735` (`MutPred`) to `2,525` in zero-shot (complete-case `N = 1,129`) and `1,702` (`MPC`) to `2,500` in supervised (complete-case `N = 1,374`). On `N = 1,129`, point-estimate Top-3 reverses (`PoET (200M)`: `0.9200 (#1) -> 0.9168 (#3)`; `GEMME`: `0.9191 (#3) -> 0.9192 (#1)`), and `MutPred` swaps `#10 <-> #9` with `MutationAssessor`. In supervised `N = 1,374`, `BayesDel (noAF)` vs `VARITY (R)` swap `#6/#7` and `MPC` vs `PolyPhen2 (HDIV)` swap `#13/#14`.
     - **Honest statistical caveat (Rule 8):** On `N = 1,129`, `GEMME - PoET (200M) = +0.0024` has a 95% bootstrap CI of `[-0.0039, +0.0088]` and paired Wilcoxon `p = 0.6852` (statistically tied). By contrast, `MutPred`'s missingness (`790` smaller genes, median 7 variants vs 15 present), equal per-gene weighting across `883/2,525` (`35.0%`) tiny genes with `min(n_pos, n_neg) == 1`, and uncalibrated cross-gene pooling (`TranceptEVE L` pooled AUC `0.8757` vs `0.9198` per-gene mean; `MutPred` `#1` at `0.9173` pooled) reflect genuine evaluation sensitivity.
   - **Draft:** `ISSUE_PG_02_incomplete_assay_gene_coverage.md`

---

### B. scPerturBench (`bm2-lab/scPerturBench`)

3. **Leaked `CPA` still present in all four genetic summary CSVs despite `README.md:107` — CONFIRMED**
   - **Repro:** `repro/repro_scpb_cpa_retained.py`
   - **Verified:** `all_dataset_genetic_seed1-3.csv` have `0` CPA rows, matching `README.md` lines 103–107. Yet all four summary CSVs in `Results/Genetic_perturbation/` still contain `CPA`:
     - `genetic_single_performance_top100.csv`: `2,181` CPA rows, **#1** by dataset mean (`0.5508` vs `scGPT` `0.5356`, `scouter` `0.5207`) and **#1** by perturbation mean (`0.6642` vs `GenePert` `0.5098`).
     - `genetic_single_performance_top5000.csv`: `2,181` CPA rows (`0.5787` ds mean #9; `0.6855` pert mean #5).
     - `genetic_combo_performance_top100.csv`: `654` CPA rows (`0.4962` ds mean #11).
     - `genetic_combo_performance_top5000.csv`: `654` CPA rows (`0.5871` ds mean #9).
   - **Extra finding:** Because `CPA` hit the per-perturbation min-max boundary `1.0` on **1,254/2,181** (`mse_score`) and **1,548/2,181** (`was_score`) rows in `genetic_single_performance_top100.csv` (and **1,638/2,181** `mse_score` rows in `top5000`), simply deleting `CPA` rows leaves remaining methods' normalized scores and `Rank_*`/`ac_rank` anchored to leaked `CPA`. Summary tables must be regenerated from `all_dataset_genetic_seed1-3.csv`.
   - **Draft:** `ISSUE_SCPB_01_cpa_still_in_genetic_summary_csvs.md`

4. **Silent zero-filling and constant mean-imputation of missing method-dataset runs in `Results/` (+ `pearson_distance` clipping) — CONFIRMED (with crucial detail on `biolord` IID & `PRnet`)**
   - **Repro:** `repro/repro_scpb_zero_and_mean_imputation.py`
   - **Verified:**
     - **Genetic single (`scFoundation`):** Missing **589/727** single perturbations in seed 1, **601/727** in seed 2, **621/727** in seed 3 (**1,811/2,181 = 83.0%** across seeds 1–3; note S4 quoted `589/958 = 61.5%` which mixed single missing count with single+combo denominator in seed 1) in `all_dataset_genetic_seed1-3.csv`. All 1,811 missing rows are filled with `0.0000` in `genetic_single_performance_top100.csv` and `top5000.csv` (dragging `scFoundation` to `0.2491`/`0.1832` vs `0.3202`/`0.2355` on the 7 covered datasets).
     - **Chemical single (`cycleCDR`, `PRnet`, `biolord`, `chemCPA`):** `cycleCDR` has `0/152` rows in `sciplex3_K562` seed 3 in `all_dataset_chemical.csv`, silently zero-filled (`152` all-zero rows in `top100`/`top5000`, dropping `K562` `ac_score` to `0.3408` vs `0.5951`/`0.6227` on `A549`/`MCF7`). Partially missing `(DataSet, seed)` runs are **mean-imputed** with identical repeated rows: `PRnet` (`120` repeated rows per cell line), `biolord` (`2` in `MCF7` seed 1), `chemCPA` (`2` in `MCF7` seed 1).
     - **Cellular IID (`biolord`):** Contrary to S4's assumption that `biolord` lacked IID evaluations on `Afriat`, `McFarland`, and `crossSpecies`, `cellular_iid_distance0.csv` **actually contains** `bioLord` results on all 12 datasets (including `pcc = 0.9908` on `Afriat`, `0.6802` on `McFarland`, and #1 `0.9965` on `crossSpecies`!). Yet `cellular_iid_performance_top100.csv` and `top5000.csv` zeroed out `cor`, `edistance_score`, `mse_score`, `was_score`, and `sym_score` on those 3 datasets, dropping `biolord` from **#3 (`0.7177` top100, `0.8079` top5000)** on the 9 non-zeroed datasets to **#9 (`0.5389` top100, `0.6060` top5000)**.
     - **`pearson_distance` clipping:** Raw CSVs (`all_dataset_genetic_seed1-3.csv`, `all_dataset_chemical.csv`) clip negative PCC-delta at `0.0` (`5,585/11,340 = 49.3%` of `baseControl` rows = `0.0`), while summary CSVs (`*_performance_top100/5000.csv`) keep signed negative `cor` (`1,061/2,181 = 48.6%` < 0).
   - **Draft:** `ISSUE_SCPB_02_silent_zero_and_mean_imputation_in_results.md`

5. **Chemical code bugs (`myPRnet.py` overwriting split 3; `myCPA_chemical.py` unpatched view/key bug) — CONFIRMED**
   - **Repro:** `repro/repro_scpb_chemical_code_bugs.py`
   - **Verified:**
     - `Perturbation_generalization/Chemical/myPRnet.py:30-31, 37-38, 64-65` copy `split_ood_multi_task2` onto `split_ood_multi_task3`. In `all_dataset_chemical.csv`, `PRnet` seed 3 has **exact set equality** `PRnet_seed3 == seed2_OOD ∩ seed3_OOD`: **32/152** perturbations (**78.9% dropped**) on `sciplex3_A549`, `sciplex3_K562`, `sciplex3_MCF7`, and **11/17** (**35.3% dropped**) on `sciplex3_comb`.
     - `Perturbation_generalization/Chemical/myCPA_chemical.py:103-106` still indexes `adata.obs['split']` (instead of `split_key = f'split_ood_multi_task{seed}'`) and performs an AnnData slice view assignment (`adata[...].X = tmp.values`), whereas genetic `myCPA.py:159-174` was already patched after Issue #1.
   - **Draft:** `ISSUE_SCPB_03_chemical_myPRnet_split3_and_myCPA_view_bug.md`

6. **Evaluation asymmetry (`generateExp(..., control_std)` variance injection vs deterministic mean predictions under Euclidean `edistance`) — CONFIRMED**
   - **Repro:** `repro/repro_scpb_generateExp_edistance.py`
   - **Verified:** Linear/simple baselines (`myGenePert.py:91-95,116`, `linearModel.py:79-83,104`, `linearModel_mean.py:113-117,137`, `genetic_baseMLP.py:215-219,240`, `genetic_baseReg.py:153-157,179`, `chemical_baseMLP.py:180-184,207`, `chemical_baseReg.py:100-104,127`, `linearModel_mean_chemical.py:99-103,126`, plus `Cellular_context_generalization` `baseControl.py`/`trainMean.py`) sample synthetic cell dispersion via `generateExp(cellNum, means, control_std)`. Deep/foundation models (`myscGPT1.py`, `myscouter.py`, `myscFoundation.py`, `myChemCPA.py`, `myCycleCDR_pre.py`, etc.) output denoised means or low-variance predictions. Under Euclidean energy distance ($2\mathbb{E}\|X-Y\|_2 - \mathbb{E}\|X-X'\|_2 - \mathbb{E}\|Y-Y'\|_2$), a zero-variance oracle predicting the exact perturbed mean gets `edistance = 14.812` at `D=5000` (13× worse than `baseControl` `1.134`), vs `0.237` with `generateExp`. Empirically in `DEG=5000` dataset means: `scGPT` has `mse = 0.003803` (`mse_score = 0.9531`, #1) and `pcc-delta = 0.3288` (#2) yet `edistance = 11.19` (vs `GenePert` `0.77`, `baseControl` `0.91`), placing #11/#12 behind `baseControl`; `chemCPA` has #1 `mse = 0.003777` and #1 `pcc-delta = 0.3976` yet `edistance = 14.83` (vs `baseControl` `1.03`), placing #5 behind #3 `baseControl`.
   - **Literature context (API-verified):** Unlike Sigliano (2026, `10.21203/rs.3.rs-9955756/v1`, single-dataset `pearson_distance` sign replay) and Li et al. (2026, *Sci. Adv.* `10.1126/sciadv.aed3414`, general variance compression across models on separate benchmarks), this documents the exact code-level asymmetry (`generateExp` applied asymmetrically to baseline scripts) within `scPerturBench` itself.
   - **Draft:** `ISSUE_SCPB_04_generateExp_edistance_variance_asymmetry.md`
