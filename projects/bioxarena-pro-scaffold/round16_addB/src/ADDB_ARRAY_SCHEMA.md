# Round 16 Addendum B (`Occam`) — Auditable Array, Latent, Barcode & Summary Export Schema (`ADDB_ARRAY_SCHEMA.md`)

**Protocol:** `PREREG_I_ADDENDUM_B_DRAFT_v2.md` (locked in `PREREG_I_ADDENDUM_B_LOCK.md`, Lock `L1`/`L2`, `DEVIATION_R16_ADDB_D1`, and Gate 5A Item 3 `CPU float32` device lock).  
**Model & Checkpoint:** `Occam` (`rhj_fewshot/all_sources/seed_42`, `best_model.pth` SHA-256 `11d5ba5515e79710b3a952532c60aa29515046117292528c59d452a0937685ff`).  
**Locked Compute Device:** `cpu` (`CPU float32`, `torch.set_num_threads(4)`, deterministic).  
**Purpose:** Binding export contract shared by Main Line (`round16_addB/results/`) and Clean-Room Audit (`A-T3`) so Master Controller can verify three-party parity in strict staged order (`PREREG_I_ADDENDUM_B_DRAFT_v2.md §6.2`):
1. **Layer 1 (Control Barcodes):** Exact `100%` ordered string equality against locked Round 16 `control_barcodes_r16.json`;
2. **Layer 2 (`z_ctrl` Encoder Latents):** `max|Z_ctrl^Main - Z_ctrl^AT3| <= 1e-5` checked **before** predicted deltas;
3. **Layer 3 (Predicted Delta Arrays):** `max|delta_pred^Main - delta_pred^AT3| <= 1e-4` across every axis, baseline, and arm;
4. **Layer 4 (Same-Array Metric Recomputation):** `|metric^Controller - metric^Party| <= 1e-10` directly from exported `.npz` files;
5. **Layer 5 (Cross-Party Metrics & Verdicts):** `|metric^Main - metric^AT3| <= 1e-5` AND `100%` verdict parity (`H-I1_Occam`, `H-I2_Occam`, `H-I3a`, `H-I3b`, `pairwise_common_axis`).

---

## 1. Directory Layout (Mirrors Round 16 `round16/results/`)

```text
<results_dir>/                                      # Main Line: round16_addB/results/
├── control_barcodes_addb.json                      # Exact locked R16 barcodes reused across all pools
├── summary_addb.json                               # Full metrics, CIs, H-I1_Occam, H-I2_Occam, H-I3a, H-I3b
├── z_ctrl/
│   ├── z_ctrl_manifest_sha256.json                 # {filename: sha256} for all z_ctrl_*.npz
│   └── z_ctrl_<pool_key>.npz                       # float32 (512, 512) encoder latents per unique pool/axis/mask
└── arrays/
    ├── a549_ko/                                    # Group 2 Primary — CRISPRko (N=403, G=1837)
    ├── dld1/                                       # Group 2 Primary — CRISPRi  (N=75,  G=1832)
    ├── cd4t_rest/                                  # Group 2 Primary — CRISPRi  (N=245, G=1830)
    ├── gse295214_h1/                               # Group 3 Holdout — CRISPRi  (N=406, G=1837)
    └── feng_ipsc/                                  # Sensitivity Only (Q4)      (N=88,  G=1575)
        ├── observed.npz                            # Exact Round 16 observed-side keys & values
        ├── preds_SENS_FCT_full_r9_6_incl_vcc2025_h1.npz   # gse295214_h1 ONLY (copied from R16 Group 3)
        ├── preds_<arm>_{matched,mis_<s>,resample}.npz     # All Occam delta arrays (Section 4)
        └── array_manifest_sha256.json              # {filename: sha256} for observed.npz + all preds_*.npz
```

> **Group 1 Exclusion (`§2.1`):** `nadig_hepg2` and `nadig_jurkat` were observed during supervised training of `rhj_fewshot/all_sources/seed_42`. They are **NEVER** evaluated as targets and have no directory under `arrays/`; they enter strictly as `seen_by_checkpoint` mismatched control pools (`z_ctrl/` and `mis_nadig_hepg2` / `mis_nadig_jurkat` arms).

---

## 2. Target Lines, Exact `(N_eval^R16, |G_eval^R16|)` Dimensions, and `§2.3` Mismatched Strata

Every evaluated target line uses the **exact locked target list (`target_names`, canonical HGNC sorted ascending)** and **exact locked evaluation gene list (`gene_names`, canonical HGNC sorted ascending)** from frozen Round 16 `observed.npz` with **zero dropped targets and zero dropped genes** (`100.0%` matrix parity with `STATE_SE`):

| Target Line | `group_key` | `role` | Modality | `N_eval^R16` (Tier 1 / Tier 2) | `|G_eval^R16|` (`RHJ_6983` cov.) | Stratum A: `seen_by_checkpoint` | Stratum B: `unseen_by_checkpoint` | Frozen R16 Reference Dir (for `observed.npz` & `STATE_SE`) |
|---|---|---|---|---:|---:|---|---|---|
| `a549_ko` | `group2_crisprko` | `primary` | `CRISPRko` | `403` (`402` T1 + `1` T2 `RBM14`) | `1,837` (`100%`) | `nadig_hepg2`, `nadig_jurkat` | `cd4t_rest`, `dld1`, `feng_ipsc` | `round16/results/arrays/a549_ko/` |
| `dld1` | `group2_crispri` | `primary` | `CRISPRi` | `75` (`75` T1 + `0` T2) | `1,832` (`100%`) | `nadig_hepg2`, `nadig_jurkat` | `a549_ko`, `cd4t_rest`, `feng_ipsc` | `round16/results/arrays/dld1/` |
| `cd4t_rest` | `group2_crispri` | `primary` | `CRISPRi` | `245` (`245` T1 + `0` T2) | `1,830` (`100%`) | `nadig_hepg2`, `nadig_jurkat` | `a549_ko`, `dld1`, `feng_ipsc` | `round16/results/arrays/cd4t_rest/` |
| `gse295214_h1` | `group3_crispri` | `prospective_group3` | `CRISPRi` | `406` (`406` T1 + `0` T2) | `1,837` (`100%`) | `nadig_hepg2`, `nadig_jurkat` | `a549_ko`, `cd4t_rest`, `dld1`, `feng_ipsc` | `round16/results/group3/arrays/gse295214_h1/` |
| `feng_ipsc` | `sensitivity_only` | `sensitivity_only` | `CRISPRi` | `88` (`88` T1 + `0` T2) | `1,575` (`100%`) | `nadig_hepg2`, `nadig_jurkat` | `a549_ko`, `cd4t_rest`, `dld1` | `round16/results/sensitivity/arrays/feng_ipsc/` |

### Target & Gene Ordering Rules
1. **Target order (`target_names`, shape `(P,)`):** Exact ordered 1D Unicode array from frozen Round 16 `<r16_line_dir>/observed.npz["target_names"]` (strictly sorted canonical HGNC symbols). Perturbations are resolved via the locked Two-Tier `pert2id` rule (`§3.5`): Tier 1 unique canonical match in `benchmarks/rhj_fewshot/perturbations.csv` (`eligible == 1`), falling back to Tier 2 direct `ckpt["pert2id"]` (`RBM14` on `a549_ko`). Both `resolved_pert_keys` and `resolved_pert_ids` are exported inside every `preds_*.npz` in the identical row order `0..P-1`.
2. **Gene order (`gene_names`, shape `(G,)`):** Exact ordered 1D Unicode array from frozen Round 16 `<r16_line_dir>/observed.npz["gene_names"]` (strictly sorted canonical HGNC symbols). Output columns are projected from Occam's `6,983` RHJ Ensembl outputs canonicalized through frozen Round 10 `HGNCBridge` (`0` collisions, `100%` coverage across all 5 target lines).

---

## 3. Observed-Side File (`arrays/<line>/observed.npz`)

Identical to Round 16's `observed.npz` schema (all 18 keys preserved verbatim in exact name, shape, dtype, and value from frozen Round 16 so existing Round 16 readers and Controller checkers work with zero translation):

| Key | Shape | Dtype | Description |
|---|---:|---|---|
| `delta_obs` | `(P, G)` | `float64` | Target observed perturbation response matrix (`log1p(norm_1e4)`) |
| `delta_half_a` | `(P, G)` | `float64` | Half-A observed perturbation response matrix |
| `delta_half_b` | `(P, G)` | `float64` | Half-B observed perturbation response matrix |
| `half_valid` | `(P,)` | `bool` | Per-perturbation split-half validity mask |
| `delta_transfer` | `(P, G)` | `float64` | `Transfer_Mean` orthogonalization reference (`source_pool` mean) |
| `source_grand_mean` | `(G,)` | `float64` | Grand mean across source pool (`r_dorth_cal`) |
| `fct_delta_pred` | `(P, G)` | `float64` | Locked `FCT` multiplicative baseline predictions |
| `e1_constmed_delta_pred` | `(P, G)` | `float64` | `E1` constant-median baseline predictions |
| `target_indices` | `(P,)` | `int64` | Index of each target gene inside `gene_names` (`-1` if absent; self-masking) |
| `target_names` | `(P,)` | `str` | Canonical HGNC target symbols (exact Round 16 order) |
| `gene_names` | `(G,)` | `str` | Canonical HGNC evaluation gene symbols (exact Round 16 order) |
| `mu_target` | `(G,)` | `float64` | Target baseline expression vector |
| `mu_sources` | `(S_src, G)` | `float64` | Source-pool baseline expression matrix (`S_src = 6`, or `5` on `gse295214_h1`) |
| `matched_ctrl_barcodes` | `(512,)` | `str` | Exact locked Round 16 primary control barcodes (`seed = 20261601`) |
| `resample_ctrl_barcodes` | `(512,)` | `str` | Exact locked Round 16 disjoint resample control barcodes (`seed = 20261602`) |
| `n_ctrl_total` | `()` | `int64` | Total available raw NTC cells in target line |
| `n_ctrl_sampled` | `()` | `int64` | Sampled primary control count (`512`) |
| `ctrl_sampling_fraction` | `()` | `float64` | `n_ctrl_sampled / n_ctrl_total` |

Additionally, on `gse295214_h1` only, `arrays/gse295214_h1/preds_SENS_FCT_full_r9_6_incl_vcc2025_h1.npz` is copied byte-for-byte from `round16/results/group3/arrays/gse295214_h1/preds_SENS_FCT_full_r9_6_incl_vcc2025_h1.npz` (keys: `delta_pred`, `delta_transfer`, `source_grand_mean`, `target_names`, `gene_names`, `model_key`) for the fixed-universe full-`r9_6` `FCT` sensitivity check.

---

## 4. Predicted Delta Arrays (`arrays/<line>/preds_<arm>_{matched,mis_<s>,resample}.npz`)

### 4.1 Per-File Internal Keys (Identical Across All Axes, Baselines, and Arms)
Every `preds_*.npz` file is written by `export_addb.py::build_occam_delta_payload` (`np.savez_compressed`, `allow_pickle=False`) with the exact 15 keys:

| Key | Shape | Dtype | Description |
|---|---:|---|---|
| `delta_pred` | `(P, G)` | `float64` | Predicted perturbation delta matrix (exact `target_names` × `gene_names` order) |
| `target_names` | `(P,)` | `str` | Canonical HGNC target symbols (`== observed.npz["target_names"]`) |
| `gene_names` | `(G,)` | `str` | Canonical HGNC gene symbols (`== observed.npz["gene_names"]`) |
| `resolved_pert_keys` | `(P,)` | `str` | Occam `pert2id` string keys selected by Two-Tier rule (`§3.5`) |
| `resolved_pert_ids` | `(P,)` | `int64` | Occam integer `pert2id` node IDs (`0..18215`) |
| `model_key` | `()` | `str` | `"Occam"` |
| `checkpoint_id` | `()` | `str` | `"rhj_fewshot/all_sources/seed_42"` |
| `input_axis_mode` | `()` | `str` | One of `"rhj_6983_axis"`, `"pairwise_common_axis"`, `"full_44728_api"` |
| `baseline_mode` | `()` | `str` | One of `"reconstructed"` (primary) or `"real"` (sensitivity) |
| `control_source` | `()` | `str` | Control pool line name (e.g. `"a549_ko"`, `"nadig_hepg2"`, `"a549_ko_within_line_resample_seed20261602"`) |
| `control_barcodes` | `(512,)` | `str` | Exact ordered `512` control-cell barcodes used to generate `delta_pred` |
| `pair_counterpart` | `()` | `str` | `""` for `rhj_6983_axis`/`full_44728_api`; counterpart line name `s` (for matched) or `t` (for mis) on `pairwise_common_axis` |
| `n_active_rhj_genes` | `()` | `int64` | Number of active RHJ-axis input genes populated (`<= 6983`; for `full_44728_api`, total matched encoder genes populated) |
| `axis_mask_sha256` | `()` | `str` | SHA-256 of the `uint8` boolean RHJ mask (length `6,983`); `""` for `full_44728_api` |
| `delta_sha256` | `()` | `str` | SHA-256 of contiguous `float64` `delta_pred.tobytes()` |

### 4.2 Complete Arm & Filename Specification Per Scoring Axis

Let $\mathcal{S}(t) = \text{seen\_by\_checkpoint}(t) \cup \text{unseen\_by\_checkpoint}(t)$ denote target $t$'s locked mismatched sources from Section 2 (`5` sources for `a549_ko`, `dld1`, `cd4t_rest`, `feng_ipsc`; `6` sources for `gse295214_h1`).

#### Axis 1 — `rhj_6983_axis` (Primary Input Axis, `§3.3.3`)
1. **Primary baseline (`baseline="reconstructed"`, `§3.4`):**
   - `preds_Occam_rhj_6983_axis_reconstructed_matched.npz`
   - `preds_Occam_rhj_6983_axis_reconstructed_mis_<s>.npz` for every $s \in \mathcal{S}(t)$
   - `preds_Occam_rhj_6983_axis_reconstructed_resample.npz`
2. **Sensitivity baseline (`baseline="real"`, `§3.4`, zero extra forward passes):**
   - `preds_Occam_rhj_6983_axis_real_matched.npz`
   - `preds_Occam_rhj_6983_axis_real_mis_<s>.npz` for every $s \in \mathcal{S}(t)$
   - `preds_Occam_rhj_6983_axis_real_resample.npz`

*(14 files per 5-source target; 16 files for `gse295214_h1`.)*

#### Axis 2 — `pairwise_common_axis` (Mandatory `P4` Gene-Coverage Sensitivity, `§3.3.4`)
- Baseline: primary `baseline="reconstructed"`.
- For every mismatched source $s \in \mathcal{S}(t)$ whose pairwise common input mask $\mathcal{A}_{\mathrm{common}}(t, s) = \mathcal{A}_{\mathrm{RHJ}} \cap \mathcal{G}_{\mathrm{raw}}(t) \cap \mathcal{G}_{\mathrm{raw}}(s)$ has **`n_union_absent > 0`** (`n_common < 6,983`):
  - `preds_Occam_pairwise_common_axis_pair_<s>_reconstructed_matched.npz` (target `t` primary controls masked to $\mathcal{A}_{\mathrm{common}}(t, s)$, `control_source = t`, `pair_counterpart = s`)
  - `preds_Occam_pairwise_common_axis_pair_<s>_reconstructed_mis_<s>.npz` (source `s` primary controls masked to $\mathcal{A}_{\mathrm{common}}(t, s)$, `control_source = s`, `pair_counterpart = t`)
- **Pairs with `n_union_absent == 0` (`n_common == 6,983`):** Strictly among full-coverage lines `{dld1, feng_ipsc, nadig_hepg2, nadig_jurkat}` (3 pairs on `dld1`: `feng_ipsc`, `nadig_hepg2`, `nadig_jurkat`; 3 pairs on `feng_ipsc`: `dld1`, `nadig_hepg2`, `nadig_jurkat`). Their common axis is identical to `rhj_6983_axis`; no duplicate `.npz` files are written, and `summary_addb.json` evaluates those pairs directly from `preds_Occam_rhj_6983_axis_reconstructed_{matched,mis_<s>}.npz` with `"reused_primary_rhj_6983": true`.
- **Within-line `resample` on `pairwise_common_axis`:** Because target $t$ and its disjoint within-line resample share the exact same raw feature set $\mathcal{G}_{\mathrm{raw}}(t)$, $\mathcal{A}_{\mathrm{common}}(t, t) = \mathcal{A}_{\mathrm{RHJ}} \cap \mathcal{G}_{\mathrm{raw}}(t)$, which is identical to `preds_Occam_rhj_6983_axis_reconstructed_resample.npz` (no separate resample files).
- **Exact Pair Counts Requiring `pairwise_common_axis` Files:**
  - `a549_ko` (`6,950/6,983`): all `5` sources $\implies$ `10` files (`6,557` with `cd4t_rest`; `6,950` with `dld1`, `feng_ipsc`, `nadig_hepg2`, `nadig_jurkat`).
  - `cd4t_rest` (`6,570/6,983`): all `5` sources $\implies$ `10` files (`6,557` with `a549_ko`; `6,570` with `dld1`, `feng_ipsc`, `nadig_hepg2`, `nadig_jurkat`).
  - `dld1` (`6,983/6,983`): `2` sources (`a549_ko` `6,950`, `cd4t_rest` `6,570`) $\implies$ `4` files (`3` remaining sources reuse primary `6,983`).
  - `feng_ipsc` (`6,983/6,983`): `2` sources (`a549_ko` `6,950`, `cd4t_rest` `6,570`) $\implies$ `4` files (`3` remaining sources reuse primary `6,983`).
  - `gse295214_h1` (`6,978/6,983`, `5` RHJ genes absent): all `6` sources $\implies$ `12` files (`n_common <= 6,978` for every pair).

#### Axis 3 — `full_44728_api` (Secondary Input-Axis Sensitivity, `§3.3.5`)
- Baseline: primary `baseline="reconstructed"`.
- Row sums computed strictly over the `44,728`-matched Ensembl columns (`GeneSpace.prepare` convention, `preprocess_addb.py::scatter_full_44728_api`).
- Exact filenames per target line $t$:
  - `preds_Occam_full_44728_api_reconstructed_matched.npz`
  - `preds_Occam_full_44728_api_reconstructed_mis_<s>.npz` for every $s \in \mathcal{S}(t)$
  - `preds_Occam_full_44728_api_reconstructed_resample.npz`

*(7 files per 5-source target; 8 files for `gse295214_h1`.)*

---

## 5. Encoder Latent Export (`z_ctrl/z_ctrl_<pool_key>.npz` & `z_ctrl_manifest_sha256.json`)

Every unique encoded control pool `(input_axis_mode, pool_line, pool_kind, active_mask)` is saved ONCE via `export_addb.py::save_occam_z_ctrl_latents` for Layer-2 verification (`max|Z_ctrl^Main - Z_ctrl^AT3| <= 1e-5`) **prior to** comparing predicted deltas.

### 5.1 Per-File Keys
| Key | Shape | Dtype | Description |
|---|---:|---|---|
| `Z_ctrl` | `(512, 512)` | `float32` | Contiguous `float32` Occam encoder latents (`encoder(X_ctrl)`) |
| `control_barcodes` | `(512,)` | `str` | Exact ordered `512` control barcodes |
| `seed` | `()` | `int64` | `20261601` (`primary`) or `20261602` (`resample`) |
| `pool_name` | `()` | `str` | Pool line name (e.g. `"a549_ko"`) |
| `input_axis_mode` | `()` | `str` | `"rhj_6983_axis"`, `"pairwise_common_axis"`, or `"full_44728_api"` |
| `n_active_rhj_genes` | `()` | `int64` | Populated RHJ input genes (or matched encoder genes for `full_44728_api`) |
| `axis_mask_sha256` | `()` | `str` | SHA-256 of `uint8` RHJ mask (`""` for `full_44728_api`) |
| `z_ctrl_sha256` | `()` | `str` | SHA-256 of contiguous `float32` `Z_ctrl.tobytes()` |
| `audit_json` | `()` | `str` | JSON string of preprocessing audit (`n_present_rhj`, `n_absent_rhj`, `n_present_but_zeroed_by_mask`, `n_zero_libsize_cells`, `z_ctrl_sha256`) |

### 5.2 Deterministic `<pool_key>` Naming
1. **`rhj_6983_axis` (12 files):**
   - `rhj_6983_axis__<line>__primary` for `<line>` in `{a549_ko, cd4t_rest, dld1, feng_ipsc, gse295214_h1, nadig_hepg2, nadig_jurkat}` (`seed = 20261601`);
   - `rhj_6983_axis__<line>__resample` for `<line>` in `{a549_ko, cd4t_rest, dld1, feng_ipsc, gse295214_h1}` (`seed = 20261602`).
2. **`pairwise_common_axis` (genuinely masked primary pools where `n_present_but_zeroed_by_mask > 0`):**
   - `pairwise_common_axis__<pool_line>__pair_<line_a>__<line_b>__primary`, where `line_a < line_b` lexicographically is the unordered pair defining $\mathcal{A}_{\mathrm{common}}(\text{line\_a}, \text{line\_b})$ and `<pool_line>` $\in \{\text{line\_a}, \text{line\_b}\}$.
   - Note: Whenever `<pool_line>` has `n_present_but_zeroed_by_mask == 0` under $\mathcal{A}_{\mathrm{common}}(\text{line\_a}, \text{line\_b})$ (for example `a549_ko` paired with `dld1`, where all `33` absent genes were already absent in `a549_ko`), its input tensor and `Z_ctrl` are bit-for-bit identical to `rhj_6983_axis__<pool_line>__primary`, so the primary latent is reused without duplicating a redundant `z_ctrl_*.npz` file.
3. **`full_44728_api` (12 files):**
   - `full_44728_api__<line>__primary` for all 7 lines (`seed = 20261601`);
   - `full_44728_api__<line>__resample` for all 5 evaluated target lines (`seed = 20261602`).

---

## 6. Control Barcode File (`control_barcodes_addb.json`)

Exact JSON consolidating all locked Round 16 control pools reused in Addendum B (zero new control sampling):
- **Sources:**
  - `round16/results/control_barcodes_r16.json`: `a549_ko` (primary + resample), `dld1` (primary + resample), `cd4t_rest` (primary + resample), `nadig_hepg2` (primary), `nadig_jurkat` (primary), `feng_ipsc` (primary);
  - `round16/results/sensitivity/control_barcodes_r16.json`: `feng_ipsc` (resample, verifying primary matches `control_barcodes_r16.json` byte-for-byte);
  - `round16/results/group3/control_barcodes_r16.json`: `gse295214_h1` (primary + resample, verifying shared mismatched pools match `control_barcodes_r16.json` byte-for-byte).
- **Schema:**
  ```json
  {
    "provenance": {
      "primary_seed": 20261601,
      "resample_seed": 20261602,
      "r16_sources_sha256": {
        "round16/results/control_barcodes_r16.json": "<sha256>",
        "round16/results/sensitivity/control_barcodes_r16.json": "<sha256>",
        "round16/results/group3/control_barcodes_r16.json": "<sha256>"
      }
    },
    "pools": {
      "<line>": {
        "primary_seed": 20261601,
        "n_ctrl_total": 4976,
        "n_ctrl_sampled": 512,
        "barcodes": ["...512 ordered strings..."],
        "barcodes_sha256": "<sha256 of newline-joined barcodes>",
        "within_line_resample": {
          "eligible": true,
          "seed": 20261602,
          "n_resample": 512,
          "overlap_with_primary": 0,
          "regime": "full_512_disjoint",
          "barcodes": ["...512 ordered strings..."],
          "barcodes_sha256": "<sha256>"
        }
      }
    }
  }
  ```

---

## 7. `summary_addb.json` Schema (`H-I1_Occam`, `H-I2_Occam`, `H-I3a`, `H-I3b`)

All bootstrap metrics use `n_boot = 4,000`, `seed = 20261601`, `delta_tost = 0.04`, two-sided `97.5%` CI (`[1.25%, 98.75%]`) for superiority/inferiority, and `95.0%` CI (`[2.5%, 97.5%]`) inside `(-0.04, +0.04)` for TOST equivalence (`valid_frac >= 0.95`), sharing the **exact identical** bootstrap draw matrix `np.random.default_rng(20261601).integers(0, N_eval, size=(4000, N_eval))` across `Occam`, `STATE_SE`, and `FCT`.

### 7.1 Top-Level Keys
- `protocol`: `"PREREG_I_ADDENDUM_B_v2"`
- `model_key`: `"Occam"`
- `checkpoint_id`: `"rhj_fewshot/all_sources/seed_42"`
- `ckpt_sha256`: `"11d5ba5515e79710b3a952532c60aa29515046117292528c59d452a0937685ff"`
- `device`: `"cpu"`
- `torch_num_threads`: `4`
- `seed`: `20261601`, `resample_seed`: `20261602`, `n_boot`: `4000`, `delta_tost`: `0.04`
- `gates_verified`: `{code_freeze_path, go_path, lock_path, human_approval_path}`
- `control_barcodes_file`: `"control_barcodes_addb.json"`
- `z_ctrl_dir`: `"z_ctrl"`
- `z_ctrl_manifest_sha256`: `{...}`
- `r16_state_se_reference_sha256`: `{<line>: {<arm_key>: "<sha256>"}}`
- `group_adjudications`: strictly separated (never pooled across groups or modalities):
  - `group2_crisprko`: `{K: 1, lines: ["a549_ko"], H_I1_Occam_verdict_counts, H_I1_Occam_k_over_K, H_I2_Occam_verdict_counts, mandatory_disclosures}`
  - `group2_crispri`: `{K: 2, lines: ["dld1", "cd4t_rest"], H_I1_Occam_verdict_counts, H_I1_Occam_k_over_K, H_I2_Occam_verdict_counts, mandatory_disclosures}`
  - `group3_crispri`: `{K: 1, lines: ["gse295214_h1"], H_I1_Occam_verdict_counts, H_I1_Occam_k_over_K, H_I2_Occam_verdict_counts, mandatory_disclosures}`
  - `sensitivity_only`: `{K: 1, lines: ["feng_ipsc"], H_I1_Occam_verdict_counts, H_I1_Occam_k_over_K, H_I2_Occam_verdict_counts, mandatory_disclosures}`
- `lines`: `{<line>: <LineEntry>}`

### 7.2 Per-Line Keys (`summary_addb.json["lines"][<line>]`)

```text
lines[<line>]
├── line_name, group_key, role, modality, checkpoint_id
├── n_eval_targets, n_eval_genes, n_half_valid, n_self_masked
├── strata_sources: {seen_by_checkpoint: [...], unseen_by_checkpoint: [...]}
├── pert_resolution_audit: {n_tier1, n_tier2, tier2_targets, bridged_pairs, unresolved_targets}
├── output_gene_projection_audit: {n_target_genes, n_mapped, n_missing, rhj_hgnc_collisions}
├── mandatory_disclosures: [...]
├── array_manifest_sha256: {...}
├── primary                                      # rhj_6983_axis, baseline="reconstructed", r_orth_cal
│   ├── Occam_matched_orth                       # {estimate, ci95_lo, ci95_hi, ci975_lo, ci975_hi, boot_se, valid_frac, n_items, n_boot, raw_r, rho_ab_orth, ceiling, mode, double_orth}
│   ├── Occam_paired_vs_mismatched               # {n_rows, matched_summary, rho_ab_orth, ceiling, double_orth, per_source: {<s>: {...}}}
│   ├── H_I1_Occam                               # Exact adjudicate_h_i1_r16 dict: {line_summary, counts, all_within_delta, any_within_delta, matched_arm_valid, matched_positive_ci975, rho_ab_orth, delta, strata: {seen_by_checkpoint, unseen_by_checkpoint}, per_source: {<s>: {category, within_delta, exact_invariance, max_abs_diff, estimate, ci975_lo, ci975_hi, ci95_lo, ci95_hi, valid_frac, stratum}}}
│   ├── FCT_matched_orth                         # Exact FCT r_orth_cal summary (must match R16 to <=1e-12)
│   ├── Occam_vs_FCT                             # Paired Occam_matched - FCT summary
│   ├── H_I2_Occam                               # Exact adjudicate_h_i2_r16 dict: {classification, estimate, ci975_lo, ci975_hi, valid_frac} + anti_masking_m62_triggered
│   └── within_line_resample_reference           # Exact compute_paired_with_within_line_resample_r16 + summarize_within_line_resample_r16 dict (delta_resample, excess_vs_resample per source, within_line_resample_status, mandatory_narrative_rule)
├── paired_vs_STATE_SE                           # §4.2 pre-registered contrasts against committed R16 STATE_SE (fixed_index_zero)
│   ├── H_I3a                                    # Matched accuracy contrast: Occam_matched - STATE_SE_matched
│   │   ├── r_orth_cal                           # {classification [OCCAM_SUPERIOR_TO_STATE_SE | STATE_SE_SUPERIOR_TO_OCCAM | INDISTINGUISHABLE], estimate, ci975_lo, ci975_hi, ci95_lo, ci95_hi, boot_se, valid_frac, n_items, n_boot, occam_summary, state_se_summary}
│   │   └── r_dorth_cal                          # Descriptive double-orthogonalized mirror (identical subkeys)
│   ├── H_I3b                                    # Control-specificity contrast: Delta_pair_Occam(s') - Delta_pair_STATE_SE(s')
│   │   ├── r_orth_cal
│   │   │   ├── per_source[<s>]                  # {stratum, delta_spec: {classification [OCCAM_MORE_CONTROL_SPECIFIC | STATE_SE_MORE_CONTROL_SPECIFIC | INDISTINGUISHABLE], estimate, ci975_lo, ci975_hi, ci95_lo, ci95_hi, boot_se, valid_frac, occam_delta_pair, state_se_delta_pair}, delta_spec_excess: {classification, estimate, ci975_lo, ci975_hi, ci95_lo, ci95_hi, boot_se, valid_frac, occam_delta_excess, state_se_delta_excess}}
│   │   │   └── per_stratum_mean[seen|unseen]    # {sources, delta_spec: {classification, estimate, ci975_lo, ci975_hi, ci95_lo, ci95_hi, boot_se, valid_frac, occam_stratum_mean, state_se_stratum_mean}, delta_spec_excess: {...}}
│   │   └── r_dorth_cal                          # Descriptive double-orthogonalized mirror (identical subkeys)
│   ├── STATE_SE_reference_calls                 # {H_I1_line_summary, H_I2_classification, matched_r_orth_cal} from R16 arrays
│   └── advance_interpretation_case              # Pre-locked §4.3 outcome pattern tag (Case_1..Case_6)
└── secondary
    ├── r_dorth_cal                              # Full Occam matched / paired / H_I1_Occam_dorth / Occam_vs_FCT / H_I2_Occam_dorth / within_line_resample
    ├── pairwise_common_axis                     # Mandatory P4 arm (§3.3.4): {per_source[<s>]: {n_common, n_union_absent, reused_primary_rhj_6983, mask_sha256, matched_orth, mismatched_orth, paired_delta, pair_classification, primary_pair_classification, verdict_changed_vs_primary}, H_I1_Occam_common_axis, H_I3b_common_axis_crosscheck, any_verdict_changed_vs_primary}
    ├── baseline_real                            # rhj_6983_axis, baseline="real": {Occam_matched_orth, Occam_paired_vs_mismatched, H_I1_Occam_real, Occam_vs_FCT, H_I2_Occam_real, within_line_resample_reference}
    ├── full_44728_api                           # full_44728_api, baseline="reconstructed": {Occam_matched_orth, Occam_paired_vs_mismatched, H_I1_Occam_full44728, Occam_vs_FCT, H_I2_Occam_full44728, within_line_resample_reference}
    ├── practical_accuracy_Occam                 # Raw Pearson Delta, centered DE-50 Pearson, cosine similarity, M2 cross-mask shift (exact R16 definitions)
    ├── tier1_trained_only_sensitivity           # Present on a549_ko (N=402 excluding Tier-2 RBM14); null on lines with 0 Tier-2 targets
    └── extra_fct_sensitivity                    # Present on gse295214_h1 (FCT_full_r9_6_incl_vcc2025_h1 fixed to N=406); null on other lines
```

---

## 8. Stage 4 Independent From-Arrays Recomputation Contract

Given `<results_dir>/arrays/` (containing `observed.npz` and `preds_*.npz`) plus frozen Round 16 `STATE_SE` arrays (`round16/results/arrays/<line>/preds_STATE_SE_fixed_index_zero_*.npz`, and corresponding `group3/` / `sensitivity/` paths), `recompute_addb_from_exported_arrays(...)` recomputes EVERY statistic, bootstrap CI, `H-I1_Occam`, `H-I2_Occam`, `H-I3a`, `H-I3b`, `pairwise_common_axis`, `baseline_real`, `full_44728_api`, `tier1_trained_only_sensitivity`, and `extra_fct_sensitivity` without loading Occam or touching raw `.h5`/`.h5ad`/`.mtx`/`.csv` files, achieving `max_abs_diff <= 1e-12` (well inside the locked `<= 1e-10` Layer-4 tolerance).
