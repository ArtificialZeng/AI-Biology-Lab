# Round 6 — Phase A Completion Report (T0, T1, T3, T7)

**Date**: 2026-11-21  
**Controller**: Claude (Round 6 Action Directive) | **Executor**: Antigravity  
**Current Repository Scope**: `AI-Biology-Lab/projects/bioxarena-pro-scaffold`  
**Core Solver Status**: Hyperparameters 100% FROZEN (`src/xatlas_panhuman_solver.py` modified ONLY on path resolution line 35 via `BIOXARENA_DIR`; zero algorithm constants or logic modified).

---

## 1. Phase A Task Completion Summary

### T0 — Frozen Codabench Submission Pack (Standalone Commit `2a8b2e0`)
- **Commit**: `2a8b2e0` (`chore(T0): freeze Codabench submission pack and SHA-256 manifest for xatlas-orion-perturb`)
- **Artifacts**:
  - `outputs/perturbation_baselines/xatlas-orion-perturb/CODABENCH_FREEZE_MANIFEST.json`
  - `outputs/perturbation_baselines/xatlas-orion-perturb/submission.csv.sha256`
- **Verified Digests & Preflight**:
  - Generator commit: `c479bc6` (solver unchanged at `a9ab830`), command: `python3 projects/bioxarena-pro-scaffold/src/xatlas_panhuman_solver.py`
  - Official `training.submission_preflight.validate_public_submission`: **PASS (`100% compliant`)**
  - `submission.csv` SHA-256: `e2613a9b6616368ba5ca29edf132fe023761f0f72072b25176724207daa6a207`
  - `solution.py` SHA-256: `792c577408e1146c30e9d087796deed11c3ab95c1f30b164cad2d7c73ccae6e4`
  - 455 `predictions/*.npz` aggregate SHA-256: `42532c396e1c71b588ab4925aa885625bbcafd8a64d372f3861a9dc94a7430d7`
- **Official Codabench Single-Submission Result (executed by Master Controller under PI authorization)**:
  - Submission ID: `1062221`, configuration `panperturb-bridge__human-ai-collab`, status `Finished`
  - Official Private Test Score (`perturbation/xatlas-orion-perturb`, metric `per_source_mean_projected_residual_effect_pearson`): **`0.0998`** (`status ok`; archive `gemini/Codabench官方结果_submission1062221.html`)
  - Pre-registered interpretation: `+0.0173` above best autonomous 22-agent run (`0.0825`, `mlevolve__claude-opus-5-5`), inside the `0.0825–0.1025` noise band (`< 0.02` official noise margin) -> **leading within noise**; NOT a statistically significant SOTA breakthrough. Measured internal-to-test shift: `0.0998 − 0.1099 = −0.0101`. Human-AI multi-round iteration (~25 tuned constants) is NOT the same comparison category as single-shot 2-hour autonomous agents. Frozen submission pack will NOT be modified or re-submitted.

### T7 — Reproducibility, Path Generalization, Cache Verification & Dependency Lock (Commit `0e7d3a7`)
- **Commit**: `0e7d3a7` (`chore(T7): remove hardcoded paths, add prepare_replogle2022_cache.py and requirements-lock.txt`)
- **Hardcoded Path Removal Audit**:
  - `grep -rnE "/Volumes/|/Users/mac/" projects/bioxarena-pro-scaffold/src/*.py` -> **0 matches (`CLEAN_ZERO_HARDCODED_PATHS`)**.
  - `src/xatlas_panhuman_solver.py` line 35 updated ONLY for `BIOXARENA_DIR` environment variable + relative fallbacks; zero algorithm constants or logic modified.
  - **Bit-for-bit parity verification after path edit**:
    - Seed 42 `S0 = 0.070151` (exact match to `reports/xatlas_10split_ablation.csv`: `0.070151`)
    - Seed 42 `S7 = 0.126020` (exact match to `reports/xatlas_10split_ablation.csv`: `0.126020`)
    - Regenerated all 455 test `.npz` vectors in memory vs frozen T0 `predictions/*.npz`: **`max_abs_diff = 0.0`**!
- **Replogle 2022 Cache Script & SHA-256 Reproducibility**:
  - Script: `src/prepare_replogle2022_cache.py`
  - Official Figshare+ DOI `https://doi.org/10.25452/figshare.plus.20029387` verified SHA-256:
    - `rpe1_normalized_bulk_01.h5ad` (`files/35775512`): `a3c5bfd0f15d63938bc80c9b8874b9cd761e3a23caf5ffe7966bae4e887ec89d`
    - `K562_gwps_normalized_bulk_01.h5ad` (`files/35773217`): `37e48c474d8b5dead4151f96ea8f5fe7bbe6beb10eeea48685b740c3f74490a2`
  - Gene selection (`5,000` genes): drop non-finite columns (2 in RPE1, 73 in K562), select top 5,000 by `var/clean_mean` (pre-computed mean expression across clean cells in published h5ad; overlap with `var/mean` top-5000 is `4,991/5,000` in RPE1 and `4,996/5,000` in K562).
  - `reproducible` definition: `anderson_darling_counts >= 50` AND `energy_test_p_value < 0.001` (rows failing only `energy_test_p_value < 0.001`: 17 in RPE1, 29 in K562).
  - **Multi-guide deduplication disclosure**: When multiple guides map to the same Ensembl target gene, the script selects the single guide maximizing `(reproducible, ad_counts, num_cells_filtered)`. Because this uses perturbation-level summary statistics across the screen, it is explicitly disclosed in `docs/PAPER_DRAFT.md` Methods.
  - **Temp-directory cache regeneration SHA-256 check (all 6 files exact match)**:
    - `RPE1_X.npy`: `e6f4ba2db5208174b8c446821a54009ea9fa366d9a2acc3bf380f6ea05a40679` (`MATCH=True`)
    - `K562_X.npy`: `3c2a646d26b5fa61ffadd1e252e9205aceb7ca0c7088604687994ab2debbc25b` (`MATCH=True`)
    - `RPE1_genes.csv`: `01559ba49a1c80134a2e09afb0e6f1dc25634c808af3d723c24bc273844515f0` (`MATCH=True`)
    - `K562_genes.csv`: `ee1f3757cef19da019e39d96ed9b68ac4c239f97927aab7535521345867d7f17` (`MATCH=True`)
    - `RPE1_meta.csv`: `4d42f0f3dc762c361f7984f7fe1f28d8870a4b3e3ca8fb8f8911a93c74955092` (`MATCH=True`)
    - `K562_meta.csv`: `c8c4c20ac229d029186090730a056fec83336044cbbdd8ad6d7ed1126d5a3020` (`MATCH=True`)
- **Pinned Dependency Lockfile**: `requirements-lock.txt` (`Python 3.13.5`, `numpy==2.2.3`, `scipy==1.15.3`, `pandas==2.2.3`, `scikit-learn==1.6.1`, `torch==2.15.0`, `h5py==3.12.1`, `matplotlib==3.10.0`, `seaborn==0.13.2`, `statsmodels==0.14.4`).

### T3 — Row-Level Stratified Paired Bootstrap 95% CIs (4,000 Resamples per Split)
- **Script**: `src/compute_row_bootstrap_ci.py`
- **Outputs**: `reports/row_bootstrap_ci_xatlas.csv`, `reports/row_bootstrap_ci_replogle.csv`, `reports/row_bootstrap_ci_summary.json`

| Dataset & Protocol | Comparison | Mean Δ | Median Per-Split 95% CI | Mean Boot SE | Splits CI > 0 | Splits CI < 0 |
|:---|:---|:---:|:---:|:---:|:---:|:---:|
| xatlas BenchA (10 dev splits) | `S7 vs S0` | +0.0485 | `[+0.0318, +0.0677]` | 0.0091 | **10/10** | 0/10 |
| xatlas BenchA (10 dev splits) | `S7 vs AE` | +0.0427 | `[+0.0256, +0.0606]` | 0.0090 | **10/10** | 0/10 |
| xatlas BenchA (10 dev splits) | `S5 vs AE` | +0.0215 | `[+0.0091, +0.0317]` | 0.0058 | **10/10** | 0/10 |
| xatlas BenchA (10 dev splits) | `S6 vs S5` | +0.0140 | `[+0.0050, +0.0244]` | 0.0050 | 9/10 | 0/10 |
| xatlas BenchA (10 dev splits) | `S1 vs S0` | +0.0063 | `[-0.0015, +0.0174]` | 0.0045 | 4/10 | 0/10 |
| xatlas BenchB Strict Union+Refit | `S7 vs S0` | +0.0457 | `[+0.0276, +0.0635]` | 0.0089 | **10/10** | 0/10 |
| xatlas BenchB Strict Union+Refit | `S7 vs AE` | +0.0398 | `[+0.0215, +0.0569]` | 0.0089 | **10/10** | 0/10 |
| Replogle 2022 Strict Union (10 splits) | `S1 vs S0` | +0.0356 | `[+0.0256, +0.0450]` | 0.0050 | **10/10** | 0/10 |
| Replogle 2022 Strict Union (10 splits) | `S5 vs AE` | **+0.0434** | **`[+0.0323, +0.0549]`** | 0.0058 | **10/10** | 0/10 |
| Replogle 2022 Strict Union (10 splits) | `S7 vs AE` | +0.0392 | `[+0.0240, +0.0537]` | 0.0072 | **10/10** | 0/10 |
| Replogle 2022 Strict Union (10 splits) | `S7 vs S0` | +0.0476 | `[+0.0326, +0.0610]` | 0.0071 | **10/10** | 0/10 |
| Replogle 2022 Strict Union (10 splits) | `S5 vs GNN` | +0.0168 | `[+0.0038, +0.0302]` | 0.0067 | 7/10 | 0/10 |
| Replogle 2022 Strict Union (10 splits) | **`S7 vs GNN`** | +0.0126 | **`[-0.0041, +0.0286]`** | 0.0082 | **3/10** | 0/10 |
| Replogle 2022 Strict Union (10 splits) | **`S6 vs S5`** | **−0.0109** | **`[-0.0169, -0.0050]`** | 0.0030 | 0/10 | **9/10** |

### T1 — Complete P0 Document Remediation & Exact T3 Backfill
- Remediated all 5 target documents (`docs/PRESS_RELEASE_CN.md`, `docs/PAPER_DRAFT.md`, `README.md`, `ACCEPTANCE.md`, `COMMUNICATION_LOG.md`) with exact T3 row-level 95% CIs backfilled (zero `[T3-CI 待填]` placeholders remaining).
- Verified via exact acceptance grep:
  `grep -nE "p ?=|untouched|从未触碰|never-seen|\+0\.010|Causal Law|因果律|Pan-Human|全人源|全面击败|全面领跑|线性壁垒" docs/PRESS_RELEASE_CN.md docs/PAPER_DRAFT.md README.md ACCEPTANCE.md COMMUNICATION_LOG.md`
  - Only matches are the exact numeric S2 ablation delta `+0.0109` (`0.0722 − 0.0613`, from `reports/xatlas_10split_ablation.csv`), which starts with `+0.010` as a 4-digit prefix (zero validation-to-test `+0.010` speculation; replaced with measured Codabench shift `−0.0101`).

---

## 2. Reproduction Commands (Zero Hardcoded Paths)

```bash
# 1. Reproduce xatlas main 10-split table & frozen test predictions
BIOXARENA_DIR=/path/to/BioXArena-Pro python3 src/xatlas_panhuman_solver.py

# 2. Reproduce Replogle 2022 cache (with SHA-256 check) & 10-split external table
REPLOGLE_DATA_DIR=/path/to/replogle2022 python3 src/prepare_replogle2022_cache.py --verify-sha256
REPLOGLE_DATA_DIR=/path/to/replogle2022 python3 src/replogle_external_validation.py --mode pan
REPLOGLE_DATA_DIR=/path/to/replogle2022 python3 src/replogle_external_validation.py --mode dl

# 3. Reproduce T3 row-level stratified paired bootstrap 95% CIs
BIOXARENA_DIR=/path/to/BioXArena-Pro REPLOGLE_DATA_DIR=/path/to/replogle2022 python3 src/compute_row_bootstrap_ci.py --mode xatlas
BIOXARENA_DIR=/path/to/BioXArena-Pro REPLOGLE_DATA_DIR=/path/to/replogle2022 python3 src/compute_row_bootstrap_ci.py --mode replogle_gnn
BIOXARENA_DIR=/path/to/BioXArena-Pro REPLOGLE_DATA_DIR=/path/to/replogle2022 python3 src/compute_row_bootstrap_ci.py --mode replogle
```

PHASE_A_DONE
