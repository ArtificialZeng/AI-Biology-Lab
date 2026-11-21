# Round 7 Execution Report — Pre-Registered External Confirmation, Official Deep-Learning & Foundation-Model Comparisons, Sensitivity Analysis, Systema Audit, and Preprint Package

**Project**: BioXArena-Pro / Pan-Human Cross-Line Perturbation Bridge (`PanPerturbBridgeModel`, frozen commit `3609811`)  
**Executor**: Antigravity  
**Master Controller**: Claude (`给Antigravity的第七轮行动指令.md`, `审核报告_第七轮_检查点1.md`)  
**Status**: All tasks `R7-0` through `R7-5` completed. Frozen solver (`src/xatlas_panhuman_solver.py`) 100% untouched (`git diff 3609811` = 0 lines).

---

## 1. Executive Summary & Revised Core Scientific Takeaway

Round 7 executed a strictly pre-registered confirmatory evaluation on a third, previously unseen CRISPRi perturbation atlas (**Nadig et al., *Nature Genetics* 2025, `GSE264667`**, HepG2 + Jurkat), added official published deep-learning (`cell-gears==0.1.2`) and foundation-model (`scGPT_human` frozen gene-embedding ridge probe) baselines on the exact same 10 strict simultaneous union-held-out splits (`SEEDS_10`), stress-tested all 10 hyperparameters of frozen `S5` by $\pm 50\%$, resolved the Systema perturbation-centroid evaluation algebra across all three atlases, and produced a complete submission-ready preprint draft (`docs/PREPRINT_DRAFT.md`) with 4 publication-grade 300-DPI figures (`PDF + PNG`).

**Verified Core Narrative (incorporating all favorable and unfavorable results):**
1. **Pre-registered confirmation on unseen Nadig 2025 (`R7-1`, independently verified by Claude in `CHECKPOINT_1_DONE`)**:
   - **H1 (`disc_raw`, `S5 − AE`)**: `+0.0487` (`0.7152` vs `0.6665`, mean 95% CI `[+0.0161, +0.0812]`, **`9/10` splits CI `> 0`**) → **CONFIRMED** (caveat: lower bounds modest `+0.0011`..`+0.0038`, one split `−0.0080`).
   - **H2 (`delta_r` metric inversion, `S5 − AE`)**: `−0.0195` (`0.1502` vs `0.1698`, `9/10` splits negative, `4/10` splits CI `< 0`) → **CONFIRMED** (caveat: weaker split-level significance than Replogle's `8/10` CI `< 0`).
   - **Main Thesis (`H1 ∧ H2`)**: **CONFIRMED**.
   - **H3 (`v5` orthogonal residual, `S5 − AE`)**: `+0.0595` (`0.1741` vs `0.1146`, **`10/10` splits CI `> 0`**) → **CONFIRMED**.
2. **Official `cell-gears==0.1.2` vs frozen `S5` (`R7-2`)**:
   - Across all 10 splits on both Replogle 2022 (`99.83%`/`99.38%` GO test coverage) and Nadig 2025 (`99.66%`/`99.79%` GO test coverage), official single-line `GEARS` underperforms both `S5` (`10/10` splits CI `> 0` on `v5`, `disc_raw`, `disc_v5` on both datasets; `10/10` Replogle and `9/10` Nadig on `delta_r`) and linear `AE` (`10/10` CI `< 0` on `v5` and `delta_r` on both datasets; `10/10` Replogle and `7/10` Nadig on `disc_raw`).
3. **Strongest Baseline — `scGPT_human` Frozen Gene-Embedding Ridge Probe (`R7-2`, critical boundary finding)**:
   - A simple linear ridge probe over frozen 33M-cell pretrained `scGPT_human` gene embeddings (`scGPT_FrozenEmb_Ridge`, MIT license) outperforms linear `AE` across **every** metric on both datasets (`Replogle`: `delta_r = 0.1802`, `v5 = 0.1127`, `disc_raw = 0.6859`; `Nadig`: `delta_r = 0.2129`, `v5 = 0.1512`, `disc_raw = 0.7109`).
   - Against `scGPT_FrozenEmb_Ridge`, zero-pretraining `S5`:
     - **Wins on `v5` orthogonal residuals**: `+0.0274` (**`10/10` CI `> 0`**) on Replogle; `+0.0230` (**`7/10` CI `> 0`**) on Nadig.
     - **Statistically ties on `disc_raw` perturbation identification**: `+0.0025` (`0/10` significant either way) on Replogle; `+0.0043` (`0/10` significant either way) on Nadig.
     - **Trails substantially on `delta_r` bulk correlation**: `−0.0609` (**`10/10` CI `< 0`**) on Replogle; `−0.0627` (**`10/10` CI `< 0`**) on Nadig.
   - **Pre-registerable future hypothesis (deferred, NOT engineered in Round 7)**: Pretrained foundation-model gene representations (strong on shared/bulk directions) and zero-pretraining cross-line HybridTOM/bridge kernels (strong on orthogonal perturbation-specific residuals) exhibit complementary strengths; testing a principled combination is explicitly deferred to a pre-registered 4th dataset.
4. **Systema perturbation-centroid evaluation (`R7-4`)**:
   - Proven algebraically and numerically across all three atlases (`xatlas`, `Replogle`, `Nadig`, including `||context_baseline − mu_rep||_2 = 0.0000` on `xatlas`) that `delta_r` is already identical to training-centroid-centered Systema Pearson (`centroid_train_r`). Even after subtracting the mean perturbation centroid (`centroid_train_r` or `centroid_val_r`), linear models (`AE`, `scGPT_Probe`) retain higher Pearson `r` than `S5` because they predict amplitude along the dominant shared inter-perturbation variance axis (PC1 of training residuals); removing that shared axis (`v5`) or testing perturbation identification (`disc_raw`) inverts rankings relative to `AE`.
5. **Sensitivity & spectral boundary (`R7-3`)**:
   - `18/20` (Replogle) and `20/20` (Nadig) $\pm 50\%$ perturbations of `S5`'s 10 hyperparameters preserve `v5` `10/10` CI `> 0` over `AE` (`18/20` Replogle and `17/20` Nadig preserve `disc_raw >= 7/10`).
   - `trans_pk -50%` (`0.75 < 1`) causes a sharp collapse on Replogle (`v5 = 0.0602`, `disc_raw = 0.6052`, `0/10`). Spectral audit diagnoses the exact mathematical mechanism: concave power `p_k = 0.75 < 1` breaks positive semi-definiteness (`min_eig(K) = -6.177`, 6 negative eigenvalues in `K + 2I`, `min_eig = -4.177`, `cond = 4,473.5` vs `17.2` at `p_k = 1.5`).

---

## 2. Task-by-Task Audit & Deliverables

### R7-0: Pre-Registration & Zero-Look Governance
- **Commits**:
  - `c5294c0`: `reports/ROUND7_PREREGISTRATION.md` (locked primary model `S5`, secondary `S7`, exploratory `S7_no_S6`, H1/H2/H3 pass/fail criteria).
  - `87f4bc2`: Addendum A (corrected 3 errata caught by Claude before data acquisition: code-authoritative `SEEDS_10 = [42, 0, 1, 2, 3, 7, 11, 99, 2024, 31337]`, exact `make_strict_union_split`, and `S3` Pan-Human Cross-Line Gene Graph definition).
  - `b9fe2bb`: Addendum B (locked GEO `GSE264667` accession provenance/SHA-256, Replogle 2022 per-`gem_group` control z-score protocol, locked Scholz-Stephens AD `ad_counts >= 50` at BH FDR `< 0.05` + 20-PC 2,000-perm energy distance `p < 0.001` + `n_cells >= 25`, top-5,000 `control_mean` readout selection, and Ensembl gene ID harmonization strictly before reading any expression matrix `X`).

### R7-1: Unseen External Confirmation — Nadig et al. (2025) `GSE264667`
- **Commit**: `f49df42` (`CHECKPOINT_1_DONE`, independently verified bit-for-bit on seeds `42, 0, 1` by Claude in `审核报告_第七轮_检查点1.md`).
- **Scripts**: `src/prepare_nadig2025_cache.py`, `src/round7_r71_nadig_confirmation.py`.
- **Reports**:
  - `reports/ROUND7_CHECKPOINT_1.md`
  - `reports/r71_nadig_confirmation_summary.json`
  - `reports/r71_nadig_5metric_pan.csv`, `reports/r71_nadig_5metric_dl.csv`, `reports/r71_nadig_5metric_combined.csv`
  - `reports/r71_nadig_row_ci.csv`
- **Cache Provenance (stored outside git under `$NADIG_DATA_DIR/processed_cache`)**:
  - HepG2 (`2,343` ENSG targets: `1,000` reproducible, `1,343` non-reproducible): `HepG2_X.npy` SHA-256 `97e164ed...`, `HepG2_meta.csv` `51786d72...`, `HepG2_genes.csv` `af43ce86...` (`902.9s`).
  - Jurkat (`2,350` ENSG targets: `1,262` reproducible, `1,088` non-reproducible): `Jurkat_X.npy` SHA-256 `4a856578...`, `Jurkat_meta.csv` `fa828370...`, `Jurkat_genes.csv` `200dca86...` (`1,270.7s`).
  - Shared readout genes: `3,644`.
- **Headline Results (`reports/r71_nadig_confirmation_summary.json`)**:

| Method | `v5` (mean ± sd) | `delta_r` (mean ± sd) | `de50_r` (mean ± sd) | `de50_sign` (mean ± sd) | `disc_raw` (mean ± sd) | `disc_v5` (mean ± sd) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| `Constant_Mean` | `0.0000 ± 0.0000` | `-0.0001 ± 0.0013` | `0.0057 ± 0.0079` | `0.5000 ± 0.0043` | `0.5000 ± 0.0000` | `0.5000 ± 0.0000` |
| `AE` (Linear Autoencoder) | `0.1146 ± 0.0078` | `0.1698 ± 0.0079` | `0.2286 ± 0.0235` | `0.6457 ± 0.0091` | `0.6665 ± 0.0110` | `0.6509 ± 0.0116` |
| `S0` (Cross-Line Ridge) | `0.1045 ± 0.0067` | `0.1610 ± 0.0067` | `0.2194 ± 0.0209` | `0.6356 ± 0.0079` | `0.6499 ± 0.0133` | `0.6373 ± 0.0133` |
| `S5` (**Primary Frozen**) | **`0.1741 ± 0.0121`** | `0.1502 ± 0.0105` | `0.2407 ± 0.0256` | `0.6459 ± 0.0080` | **`0.7152 ± 0.0129`** | **`0.7157 ± 0.0119`** |
| `S7` (Secondary Frozen) | `0.1871 ± 0.0137` | `0.1617 ± 0.0124` | `0.2584 ± 0.0294` | `0.6562 ± 0.0076` | `0.6992 ± 0.0139` | `0.6993 ± 0.0112` |

### R7-2: Official Published Deep-Learning & Foundation-Model Comparisons
- **Scripts**:
  - `src/round7_r72_official_gears.py` (official `cell-gears==0.1.2`, default hyperparams `hidden_size=64, num_go_gnn_layers=1, num_gene_gnn_layers=1, decoder_hidden_size=16, epochs=20, lr=1e-3, weight_decay=5e-4, batch_size=32`).
  - `src/round7_r72_official_scgpt.py` (official `scGPT_human` frozen gene-embedding ridge probe, MIT license; checkpoint SHA-256 `ad0252a1...`, vocab SHA-256 `ee2b2c90...`; plus hardware timing audit of full `TransformerGenerator` fine-tuning).
- **Reports**:
  - `reports/r72_official_gears_replogle_{metrics.csv,row_ci.csv,summary.json}`
  - `reports/r72_official_gears_nadig_{metrics.csv,row_ci.csv,summary.json}`
  - `reports/r72_scgpt_probe_replogle_{metrics.csv,row_ci.csv}` (committed in `19ab299`)
  - `reports/r72_scgpt_probe_nadig_{metrics.csv,row_ci.csv}` (committed in `19ab299`)
- **Summary Table (All 10 Splits `SEEDS_10`, Exact Paired Row-Level Bootstrap CIs)**:

| Dataset | Model | `v5` | `delta_r` | `de50_r` | `de50_sign` | `disc_raw` | `disc_v5` |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Replogle 2022** | `Official_GEARS` (`cell-gears==0.1.2`) | `0.0304 ± 0.0168` | `0.0486 ± 0.0253` | `0.1077 ± 0.0332` | `0.6009 ± 0.0191` | `0.5621 ± 0.0208` | `0.5446 ± 0.0172` |
| **Replogle 2022** | `scGPT_FrozenEmb_Ridge` | `0.1127 ± 0.0067` | **`0.1802 ± 0.0060`** | **`0.2646 ± 0.0140`** | **`0.6536 ± 0.0064`** | `0.6859 ± 0.0080` | `0.6807 ± 0.0085` |
| **Replogle 2022** | `S5` (Frozen) | **`0.1401 ± 0.0060`** | `0.1193 ± 0.0056` | `0.2145 ± 0.0185` | `0.6256 ± 0.0079` | **`0.6884 ± 0.0065`** | **`0.6885 ± 0.0062`** |
| **Nadig 2025** | `Official_GEARS` (`cell-gears==0.1.2`) | `0.0577 ± 0.0143` | `0.1077 ± 0.0174` | `0.1973 ± 0.0246` | `0.6462 ± 0.0095` | `0.6095 ± 0.0149` | `0.5834 ± 0.0273` |
| **Nadig 2025** | `scGPT_FrozenEmb_Ridge` | `0.1512 ± 0.0113` | **`0.2129 ± 0.0119`** | **`0.3064 ± 0.0279`** | **`0.6792 ± 0.0120`** | `0.7109 ± 0.0141` | `0.7100 ± 0.0126` |
| **Nadig 2025** | `S5` (Frozen) | **`0.1741 ± 0.0121`** | `0.1502 ± 0.0105` | `0.2407 ± 0.0256` | `0.6459 ± 0.0080` | **`0.7152 ± 0.0129`** | **`0.7157 ± 0.0119`** |

- **Key Pairwise Row-Level 95% Bootstrap CIs (10 splits)**:
  - `S5 − Official_GEARS` (Replogle): `v5 = +0.1096` (`10/10` CI `> 0`), `delta_r = +0.0706` (`10/10`), `disc_raw = +0.1263` (`10/10`), `disc_v5 = +0.1439` (`10/10`). GO-covered-only scores match full scores (`v5 = 0.0305`, `disc_raw = 0.5620`).
  - `S5 − Official_GEARS` (Nadig): `v5 = +0.1164` (`10/10`), `delta_r = +0.0425` (`9/10`), `disc_raw = +0.1057` (`10/10`), `disc_v5 = +0.1323` (`10/10`). GO-covered-only scores match full scores (`v5 = 0.0579`, `disc_raw = 0.6102`).
  - `S5 − scGPT_FrozenEmb_Ridge` (Replogle): `v5 = +0.0274` (**`10/10` CI `> 0`**), `disc_raw = +0.0025` (`0/10` pos, `0/10` neg — tie), `disc_v5 = +0.0078` (`0/10` pos, `0/10` neg — tie), `delta_r = -0.0609` (**`10/10` CI `< 0`**).
  - `S5 − scGPT_FrozenEmb_Ridge` (Nadig): `v5 = +0.0230` (**`7/10` CI `> 0`**), `disc_raw = +0.0043` (`0/10` pos, `0/10` neg — tie), `disc_v5 = +0.0057` (`0/10` pos, `0/10` neg — tie), `delta_r = -0.0627` (**`10/10` CI `< 0`**).
  - `scGPT_FrozenEmb_Ridge − AE`: Beats `AE` on ALL metrics — Replogle `delta_r = +0.0332` (`8/10`), `disc_raw = +0.0392` (`8/10`), `disc_v5 = +0.0420` (`10/10`), `v5 = +0.0160` (`6/10`); Nadig `delta_r = +0.0432` (`10/10`), `disc_raw = +0.0444` (`8/10`), `disc_v5 = +0.0591` (`10/10`), `v5 = +0.0365` (`10/10`).
- **Transparent Hardware & Compatibility Audits (no approximation or under-training)**:
  - **Full `scGPT_human` `TransformerGenerator` Fine-Tuning**: Loaded all 167 tensors on Apple MPS (fixed PyTorch MPS nested-tensor bug via `model.transformer_encoder.use_nested_tensor = False`); single-epoch wall-clock timed at `606s` (~`10.1 min/epoch/cell line` at `ctx=1536, batch=4` without CUDA `flash-attn`), projecting to `~100 hours` for 20 splits × 15 official epochs. Transparently omitted on local workstation hardware; CELLxGENE 33M-cell pretraining corpus overlap caveat (unperturbed RPE1/K562/HepG2/Jurkat transcriptional states vs CRISPRi holdout) explicitly documented in `docs/PREPRINT_DRAFT.md`.
  - **Arc Institute `STATE` (`arc-state==0.11.1`)**: Installed and audited; `ST` defaults to `pert_rep=onehot` (incompatible with simultaneous union holdout of unseen gene targets without non-commercial `SE-600M` embeddings and NVIDIA CUDA cluster training). Transparently documented in `docs/PREPRINT_DRAFT.md`.

### R7-3: Frozen `S5` $\pm 50\%$ Hyperparameter Sensitivity Analysis
- **Script**: `src/round7_r73_sensitivity.py`.
- **Reports**:
  - `reports/r73_sensitivity_replogle.{csv,json}` (summary: `r73_sensitivity_replogle_summary.json`, committed `f49df42`)
  - `reports/r73_sensitivity_nadig.{csv,json}` (summary: `r73_sensitivity_nadig_summary.json`, committed `19ab299`)
- **Key Findings**:
  - Exact `0.000e+00` numerical parity verified between baseline config `S5_Frozen` and frozen `PanPerturbBridgeModel` S5 on both datasets.
  - **Replogle 2022**: `18/20` perturbed configs preserve `v5` `10/10` CI `> 0` and `disc_raw` `8–10/10` CI `> 0` over `AE`. Exceptions: `soft_p -50%` (`0.75`) drops `disc_raw` to `5/10` (still `v5` `10/10`), and `trans_pk -50%` (`0.75`) collapses (`v5 = 0.0602`, `disc_raw = 0.6052`, `0/10`).
  - **Nadig 2025**: `20/20` perturbed configs preserve `v5` `10/10` CI `> 0` over `AE`; `17/20` preserve `disc_raw >= 7/10` (`soft_p -50%`, `krr_alpha -50%`, and `trans_pk -50%` drop to `6/10`).
  - **Mathematical Diagnosis of `trans_pk = 0.75 < 1`**: Entrywise concave exponentiation $|K_{ij}|^{0.75}$ violates positive semi-definiteness (`min_eig(K) = -6.177`, 6 negative eigenvalues in `K + 2I` with `min_eig = -4.177`, condition number `4,473.5` vs `17.2` at `p_k = 1.5`). Consequently, `p_k >= 1` is a structural PSD requirement rather than a finely tuned hyperparameter.

### R7-4: Systema Perturbation-Centroid-Centered Evaluation (All 3 Datasets)
- **Script**: `src/round7_r74_systema_eval.py`.
- **Reports**:
  - `reports/r74_systema_xatlas_{scores.csv,row_ci.csv,summary.json}` (`f49df42`)
  - `reports/r74_systema_replogle_{scores.csv,row_ci.csv,summary.json}` (`f49df42`)
  - `reports/r74_systema_nadig_{scores.csv,row_ci.csv,summary.json}` (`19ab299`)
- **Key Findings**:
  - Across all 3 atlases — including `xatlas`, where `||context_baseline − mu_rep||_2 = 0.0000` (`r = 1.000000`) — `delta_r` is numerically identical to training-centroid-centered Systema Pearson (`centroid_train_r`).
  - Companion `centroid_val_r` (centered by validation-set mean perturbation centroid) reproduces the exact same ranking (`AE > S5`):
    - `xatlas`: `S5 = 0.1098 ± 0.0164` vs `AE = 0.1308 ± 0.0182` (`delta = -0.0210`, `2/10` CI `< 0`).
    - `Replogle`: `S5 = 0.1382 ± 0.0122` vs `AE = 0.1581 ± 0.0102` (`delta = -0.0199`, `6/10` CI `< 0`).
    - `Nadig`: `S5 = 0.1596 ± 0.0173` vs `AE = 0.1787 ± 0.0110` (`delta = -0.0191`, `3/10` CI `< 0`).

### R7-5: Preprint Draft & 4 Publication-Grade 300-DPI Figures
- **Manuscript**: `docs/PREPRINT_DRAFT.md`
  - Complete IMRaD preprint with Abstract, Introduction, Results (spanning xatlas, Replogle 2022, Nadig 2025, Official GEARS, `scGPT_FrozenEmb_Ridge`, LOMO ablation, $\pm 50\%$ sensitivity + PSD diagnosis, Systema audit, and Reactome/cis-leakage controls), Discussion (including CELLxGENE overlap caveat, hardware/compatibility audit for full scGPT & STATE, and explicitly deferred pre-registerable complementarity hypothesis), Methods, Data/Code Availability, and transparent AI-assistance attribution.
  - Every reported metric cites its exact committed file in `reports/`.
- **Figure Generator**: `src/round7_r75_generate_figures.py` (reads exclusively from committed `reports/*.{csv,json}`, zero hardcoded numbers).
- **4 Visually Verified 300-DPI Figures (`PDF + PNG`)**:
  1. `docs/figures/fig1_metric_inversion_three_datasets.{pdf,png}` — 3-dataset metric inversion (`delta_r` vs `v5` vs `disc_raw`), incorporating 10-split `Official_GEARS` and 10-split `scGPT_FrozenEmb_Ridge`.
  2. `docs/figures/fig2_preregistered_confirmation_forest.{pdf,png}` — 10-split row-level paired bootstrap 95% CI forest plot for pre-registered Nadig 2025 (`H1`, `H2`, `H3` + `S5 vs scGPT_Probe` & `S5 vs Official_GEARS`).
  3. `docs/figures/fig3_lomo_ablation_and_s5_sensitivity.{pdf,png}` — Dual-dataset LOMO component ablation + $\pm 50\%$ `S5` sensitivity & `trans_pk = 0.75` PSD collapse.
  4. `docs/figures/fig4_biology_reactome_and_cis_audit.{pdf,png}` — Reactome pathway enrichment & cis-target leakage / `S6` deflation audit.

---

## 3. Reproducibility & Environment Commands

All Round 7 scripts use environment variables (`BIOXARENA_DIR`, `REPLOGLE_DATA_DIR`, `NADIG_DATA_DIR`, `SCGPT_WEIGHTS_DIR`) with zero hardcoded local paths in `src/*.py`.

```bash
export BIOXARENA_DIR=/path/to/BioXArena-Pro
export REPLOGLE_DATA_DIR=/path/to/external_data/replogle2022
export NADIG_DATA_DIR=/path/to/external_data/nadig2025
export SCGPT_WEIGHTS_DIR=/path/to/scgpt_human_weights

# 1. R7-1: Build Nadig 2025 cache & run pre-registered 10-split confirmation
python3 projects/bioxarena-pro-scaffold/src/prepare_nadig2025_cache.py
python3 projects/bioxarena-pro-scaffold/src/round7_r71_nadig_confirmation.py

# 2. R7-2: Official cell-gears==0.1.2 (in isolated GEARS venv) & scGPT_human frozen probe
$GEARS_VENV/bin/python projects/bioxarena-pro-scaffold/src/round7_r72_official_gears.py --dataset replogle
$GEARS_VENV/bin/python projects/bioxarena-pro-scaffold/src/round7_r72_official_gears.py --dataset nadig
python3 projects/bioxarena-pro-scaffold/src/round7_r72_official_scgpt.py --dataset replogle --mode probe
python3 projects/bioxarena-pro-scaffold/src/round7_r72_official_scgpt.py --dataset nadig --mode probe

# 3. R7-3: S5 +/-50% sensitivity analysis
python3 projects/bioxarena-pro-scaffold/src/round7_r73_sensitivity.py --dataset replogle
python3 projects/bioxarena-pro-scaffold/src/round7_r73_sensitivity.py --dataset nadig

# 4. R7-4: Systema perturbation-centroid evaluation across all 3 datasets
python3 projects/bioxarena-pro-scaffold/src/round7_r74_systema_eval.py --dataset xatlas
python3 projects/bioxarena-pro-scaffold/src/round7_r74_systema_eval.py --dataset replogle
python3 projects/bioxarena-pro-scaffold/src/round7_r74_systema_eval.py --dataset nadig

# 5. R7-5: Regenerate all 4 publication figures (300 DPI PDF + PNG)
python3 projects/bioxarena-pro-scaffold/src/round7_r75_generate_figures.py
```

---

## 4. Governance & Release Boundary Checklist

- [x] Frozen solver `src/xatlas_panhuman_solver.py` untouched from `3609811` (`git diff 3609811` = 0 lines).
- [x] Zero post-hoc `scGPT + S5` hybrid engineered in Round 7 (strictly deferred to a pre-registered 4th dataset).
- [x] Zero hardcoded local paths in `src/*.py` or `docs/PREPRINT_DRAFT.md`.
- [x] All temporary scratch artifacts (`reports/_tmp_nadig_dl_rows/`) removed prior to final commit.
- [x] Zero external actions taken (no bioRxiv upload, no Codabench resubmission, no media/researcher outreach — sole authority reserved for Human PI).
- [x] Both git repositories (local SSD clone and Thunderbolt Sync Volume clone) pushed, synced, and verified clean at identical HEAD commit.

ROUND7_DONE
