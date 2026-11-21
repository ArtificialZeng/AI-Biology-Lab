# [Bug] Chemical scripts: `myPRnet.py` overwrites `split_ood_multi_task3` with `split_ood_multi_task2`, and `myCPA_chemical.py` keeps unpatched view/key bug

**Commit checked:** `63fd9b95071b894c4f39c8f7eb2bd42ca9a70a77` (`main`)

## Summary

Two code bugs remain in `Perturbation_generalization/Chemical/`:

1. **`myPRnet.py` overwrites `split_ood_multi_task3` with `split_ood_multi_task2` (lines 30–31, 37–38, and 64–65):**
   In both `preData4PRnet_single` and `preData4PRnet_comb`, copy-paste lines assign `adata.obs["split_ood_multi_task2"]` to `adata.obs["split_ood_multi_task3"]`. For `seed = 3`, `PRnet` therefore trains and tests on **seed 2's OOD perturbations** (148 single / 20 combo), while `trainMean/savedModels3/result.h5ad` (line 144) supplies **seed 3's ground truth** (152 single / 17 combo). As a result, `calPerformance_chemical.py` evaluates only their intersection (`seed2_OOD ∩ seed3_OOD`): **32/152 (78.9% dropped)** on `sciplex3_A549`, `sciplex3_K562`, and `sciplex3_MCF7`, and **11/17 (35.3% dropped)** on `sciplex3_comb` in `Results/Chemical_perturbation/all_dataset_chemical.csv`.
2. **`myCPA_chemical.py` still uses `adata.obs['split']` and an AnnData slice view assignment (lines 103–106):**
   While `Perturbation_generalization/Genetic/myCPA.py:159-174` was patched after Issue #1 to use boolean NumPy masks (`adata.X[test_mask, :] = tmp.values`), `Perturbation_generalization/Chemical/myCPA_chemical.py:103-106` still indexes `'split'` (instead of `split_key = 'split_ood_multi_task{}'.format(seed)` from lines 74, 98, 109) and assigns to a slice view (`adata[...].X = tmp.values`), which does not update `adata.X` in place.

## Minimal Reproduction

```python
import pandas as pd

chem = pd.read_csv("Results/Chemical_perturbation/all_dataset_chemical.csv")
for ds in ["sciplex3_A549", "sciplex3_K562", "sciplex3_MCF7", "sciplex3_comb"]:
    s2 = set(chem[(chem["method"] == "trainMean") & (chem["DataSet"] == ds) & (chem["seed"] == 2)]["perturb"])
    s3 = set(chem[(chem["method"] == "trainMean") & (chem["DataSet"] == ds) & (chem["seed"] == 3)]["perturb"])
    pr3 = set(chem[(chem["method"] == "PRnet") & (chem["DataSet"] == ds) & (chem["seed"] == 3)]["perturb"])
    print(f"{ds}: |seed2|={len(s2)}, |seed3|={len(s3)}, |s2&s3|={len(s2 & s3)}, |PRnet_seed3|={len(pr3)}, exact_match={pr3 == (s2 & s3)}")
```

Output on `main`:
```
sciplex3_A549: |seed2|=148, |seed3|=152, |s2&s3|=32, |PRnet_seed3|=32, exact_match=True
sciplex3_K562: |seed2|=148, |seed3|=152, |s2&s3|=32, |PRnet_seed3|=32, exact_match=True
sciplex3_MCF7: |seed2|=148, |seed3|=152, |s2&s3|=32, |PRnet_seed3|=32, exact_match=True
sciplex3_comb: |seed2|=20, |seed3|=17, |s2&s3|=11, |PRnet_seed3|=11, exact_match=True
```

## Suggested Patches

```diff
--- a/Perturbation_generalization/Chemical/myPRnet.py
+++ b/Perturbation_generalization/Chemical/myPRnet.py
@@ -27,15 +27,15 @@ def preData4PRnet_single(adata):
     adata.obs["split_ood_multi_task1"] = adata.obs["split_ood_multi_task1"].astype("category")
     adata.obs["split_ood_multi_task2"] = (adata.obs["split_ood_multi_task2"].astype(str).replace("test", "valid"))
     adata.obs["split_ood_multi_task2"] = adata.obs["split_ood_multi_task2"].astype("category")
-    adata.obs["split_ood_multi_task3"] = (adata.obs["split_ood_multi_task2"].astype(str).replace("test", "valid"))
-    adata.obs["split_ood_multi_task3"] = adata.obs["split_ood_multi_task2"].astype("category")
+    adata.obs["split_ood_multi_task3"] = (adata.obs["split_ood_multi_task3"].astype(str).replace("test", "valid"))
+    adata.obs["split_ood_multi_task3"] = adata.obs["split_ood_multi_task3"].astype("category")
 
     adata.obs["split_ood_multi_task1"] = (adata.obs["split_ood_multi_task1"].astype(str).replace("ood", "test"))
     adata.obs["split_ood_multi_task1"] = adata.obs["split_ood_multi_task1"].astype("category")
     adata.obs["split_ood_multi_task2"] = (adata.obs["split_ood_multi_task2"].astype(str).replace("ood", "test"))
     adata.obs["split_ood_multi_task2"] = adata.obs["split_ood_multi_task2"].astype("category")
-    adata.obs["split_ood_multi_task3"] = (adata.obs["split_ood_multi_task2"].astype(str).replace("ood", "test"))
-    adata.obs["split_ood_multi_task3"] = adata.obs["split_ood_multi_task2"].astype("category")
+    adata.obs["split_ood_multi_task3"] = (adata.obs["split_ood_multi_task3"].astype(str).replace("ood", "test"))
+    adata.obs["split_ood_multi_task3"] = adata.obs["split_ood_multi_task3"].astype("category")
@@ -61,8 +61,8 @@ def preData4PRnet_comb(adata):
     adata.obs["split_ood_multi_task1"] = adata.obs["split_ood_multi_task1"].astype("category")
     adata.obs["split_ood_multi_task2"] = (adata.obs["split_ood_multi_task2"].astype(str).replace("val", "valid"))
     adata.obs["split_ood_multi_task2"] = adata.obs["split_ood_multi_task2"].astype("category")
-    adata.obs["split_ood_multi_task3"] = (adata.obs["split_ood_multi_task2"].astype(str).replace("val", "valid"))
-    adata.obs["split_ood_multi_task3"] = adata.obs["split_ood_multi_task2"].astype("category")
+    adata.obs["split_ood_multi_task3"] = (adata.obs["split_ood_multi_task3"].astype(str).replace("val", "valid"))
+    adata.obs["split_ood_multi_task3"] = adata.obs["split_ood_multi_task3"].astype("category")
```

```diff
--- a/Perturbation_generalization/Chemical/myCPA_chemical.py
+++ b/Perturbation_generalization/Chemical/myCPA_chemical.py
@@ -100,10 +100,12 @@ def getPredict_chemical(DataSet, seed, ood = 'ood'):
     adata = sc.read_h5ad('{}/cpa.h5ad'.format(modeldir))
     model = cpa.CPA.load(modeldir, adata = adata, use_gpu = True)
 
-    a = adata[adata.obs['split'] == 'ood'].shape[0]
-    tmp = adata[adata.obs['perturbation'] == 'control']
+    ood_mask = (adata.obs[split_key] == ood).to_numpy()
+    ctrl_mask = (adata.obs['perturbation'] == 'control').to_numpy()
+    a = int(ood_mask.sum())
+    tmp = adata[ctrl_mask].copy()
     tmp = tmp.to_df().sample(n=a, random_state=42, replace=True)
-    adata[adata.obs['split'] == 'ood'].X = tmp.values
+    adata.X[ood_mask, :] = tmp.values
 
     model.predict(adata, batch_size=2048)
```
