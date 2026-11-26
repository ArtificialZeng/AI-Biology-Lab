# Round 17 Exported Array & Embedding Schema (`R17_ARRAY_SCHEMA.md`)

This document specifies the exact directory layout, filenames, array keys, shapes, dtypes, and invariants exported by `round17/src/run_eval_r17.py` so that Master Controller and `A-T3` can independently verify every stage of the Round 17 three-party tolerance ladder (`PREREG_J §5.2`) and recompute **every** reported statistic, confidence interval, Holm–Bonferroni call, `J5a` `G4` decomposition, and `J5b` rank correlation strictly from exported `.npz` files alone (`recompute_cohort_from_exported_arrays_r17`, Stage 4 tolerance `<= 1e-10`).

---

## 1. Directory Layout (`round17/results/`)

```text
round17/results/
├── control_barcodes_r17.json                      # Stage 1 exact ordered control barcodes & batch labels per pool
├── embeddings/                                    # Stage 2 cached SE-600M X_state (2058-dim float32)
│   ├── x_state_<pool>.npz
│   └── embeddings_manifest_sha256.json
├── arrays/                                        # Stage 3 & Stage 4 auditable per-line arrays
│   └── <line>/
│       ├── observed.npz
│       ├── preds_STATE_SE_sweep56.npz
│       ├── preds_STATE_SE_g4_mismatched.npz
│       ├── preds_STATE_HVG_sweep56.npz
│       ├── preds_STATE_HVG_g4_mismatched.npz
│       ├── SENS_FCT_full_r9_6_incl_vcc2025_h1.npz # gse295214_h1 ONLY
│       └── array_manifest_sha256.json
└── summary_r17.json                               # Complete machine-readable results & adjudications
```

Evaluated `<line>` keys (`11` lines; `wtc11_i3neuron` omitted per Ruling 2; all 6 R16 discovery lines quarantined as mismatched control pools only):
- Confirmatory (`K_conf = 3`): `cd4t_stim48`, `gse295214_h1`, `pfizer_telohaec_basal`
- Within-family replication (`role = "within_family_replication"`): `pfizer_telohaec_cytokine`
- Secondary fresh (`role = "secondary_fresh"`): `a549_myllia`, `hek293_xucao`
- Descriptive micro-N (`role = "micro_n_descriptive"`): `telohaec`, `calu3`, `mac_lps`, `mac_ctrl`
- R16 subset sensitivity (`role = "r16_subset_sensitivity"`): `a549_ko_cas9pos`

---

## 2. Stage 1 & Stage 2 Artifacts

### 2.1 `control_barcodes_r17.json`
JSON dictionary keyed by `<pool>` (every evaluated line plus all 6 R16 mismatched pools: `nadig_hepg2`, `nadig_jurkat`, `a549_ko`, `dld1`, `cd4t_rest`, `feng_ipsc`):
- `ctrl_seed` (`int`): RNG seed used for stratified Hamilton control sampling (`20261701` default; configurable via `--ctrl-seed`).
- `n_ctrl_total` (`int`): Total QC-passing NTC cells available in `<pool>`.
- `n_ctrl_sampled` (`int`): Sampled control count (`512`, partitioned into `16 x 32` sets).
- `barcodes` (`list[str]`): Exact ordered cell barcode list matching rows `0..n_ctrl_sampled-1` of `X_state` and `X_hvg`.
- `batch_labels` (`list[str]`): Exact ordered experimental batch label per sampled cell.

### 2.2 `embeddings/x_state_<safe_pool>.npz`
Identical schema to frozen Round 16 `save_control_embeddings_r16`:
| Key | Shape | Dtype | Description |
|---|---|---|---|
| `X_state` | `(n_ctrl, 2058)` | `float32` | `SE-600M` control embeddings (`2048` cell + `10` dataset embedding), C-contiguous |
| `control_barcodes` | `(n_ctrl,)` | `str` | Exact ordered cell barcodes matching `control_barcodes_r17.json` |
| `seed` | `()` | `int64` | Per-pool RNG seed applied immediately before `DataLoader` iteration |
| `dataset_name` | `()` | `str` | Pool label passed to `FilteredGenesCounts` |
| `se_audit_json` | `()` | `str` | JSON audit metadata (`device`, `dtype`, `n_mapped_to_pe19790`, `x_state_sha256`, etc.) |

---

## 3. Per-Line `.npz` Schemas (`arrays/<line>/`)

Let `P = N_eval^R17` (canonically sorted ascending by HGNC target symbol) and `G = |G_eval^STATE|`. All `.npz` files are written via `np.savez_compressed` and must be readable with `allow_pickle=False`.

### 3.1 `observed.npz`
Contains all observed ground-truth deltas, split halves, `FCT` & `Transfer_Mean` baselines, unlabeled observed NTC control profile for `L2_ntc_recon_match`, locked split `s = 0` indices, and metadata:

| Key | Shape | Dtype | Description |
|---|---|---|---|
| `line_name` | `()` | `str` | Canonical line key |
| `role` | `()` | `str` | One of `confirmatory`, `within_family_replication`, `secondary_fresh`, `micro_n_descriptive`, `r16_subset_sensitivity` |
| `modality` | `()` | `str` | `CRISPRi` or `CRISPRko` |
| `checkpoint_id` | `()` | `str` | Always `zeroshot/ensemble4` in Round 17 |
| `delta_obs` | `(P, G)` | `float64` | Observed target perturbation pseudobulk log1p shifts |
| `delta_half_a` | `(P, G)` | `float64` | Observed Half A perturbation shifts |
| `delta_half_b` | `(P, G)` | `float64` | Observed Half B perturbation shifts |
| `half_valid` | `(P,)` | `bool` | Per-perturbation split-half validity indicator |
| `delta_transfer` | `(P, G)` | `float64` | Per-perturbation `Transfer_Mean` fitted on locked source pool |
| `source_grand_mean` | `(G,)` | `float64` | Source grand-mean vector for secondary `r_dorth_cal` |
| `fct_delta_pred` | `(P, G)` | `float64` | Reference `FCT` predicted shifts on `(T_eval, G_eval)` |
| `e1_constmed_delta_pred` | `(P, G)` | `float64` | Reference `E1_ConstMed` predicted shifts |
| `target_indices` | `(P,)` | `int64` | Column index of target gene inside `gene_names` (`-1` if absent from `G_eval`) |
| `target_names` | `(P,)` | `str` | Canonically sorted evaluated target symbols |
| `gene_names` | `(G,)` | `str` | Ordered evaluated canonical gene symbols (`G_eval^STATE`) |
| `mu_target` | `(G,)` | `float64` | Target control mean libnorm expression on `G_eval` |
| `mu_sources` | `(S_fct, G)` | `float64` | Source line control mean libnorm expressions (`S_fct = 5` for `gse295214_h1`, `6` otherwise) |
| `obs_ntc_mean_official_median_log1p` | `(G,)` | `float64` | Mean of sampled `n_ctrl` unlabeled NTC cells normalized via `official_median_log1p` on `G_eval` (sole observed input for `L2_ntc_recon_match`) |
| `split_s0_p1_indices` | `(floor(P/2),)` | `int64` | Locked split `s = 0` (`seed = 20261701`) `P1` row indices into `0..P-1` |
| `split_s0_p2_indices` | `(ceil(P/2),)` | `int64` | Locked split `s = 0` (`seed = 20261701`) `P2` row indices into `0..P-1` |
| `seen_sources` | `(S_seen,)` | `str` | `J5a` `seen_by_checkpoint` mismatched pool names |
| `unseen_sources` | `(S_unseen,)` | `str` | `J5a` `unseen_by_checkpoint` foreign mismatched pool names |
| `same_lineage_sources` | `(S_same,)` | `str` | `J5a` same-lineage comparator names (`["cd4t_rest"]` on `cd4t_stim48`, empty elsewhere) |
| `matched_ctrl_barcodes` | `(n_ctrl,)` | `str` | Ordered sampled control cell barcodes |
| `matched_ctrl_batch_labels` | `(n_ctrl,)` | `str` | Ordered sampled control batch labels |
| `n_ctrl_total` | `()` | `int64` | Total QC-passing NTC cells before sampling |
| `n_ctrl_sampled` | `()` | `int64` | Number of sampled control cells (`512`) |
| `ctrl_sampling_fraction` | `()` | `float64` | `n_ctrl_sampled / n_ctrl_total` |

### 3.2 `preds_<MODEL>_sweep56.npz` (`<MODEL>` in `{STATE_SE, STATE_HVG}`)
Contains all 56 selection-free matched predictions (`b = 0..55`, never dropped per `J4`), 4-fold and per-fold `'non-targeting'` pseudobulks for Ruling 3 `L2_ntc_recon_match`, label-free rule outputs, secondary neutral-batch comparators, and per-fold 56-index `r_orth_cal` summaries:

| Key | Shape | Dtype | Description |
|---|---|---|---|
| `delta_pred_56` | `(56, P, G)` | `float64` | Exact 4-fold averaged predicted deltas `Y_hat(p, b) - Y_hat('non-targeting', b)` for every `b = 0..55` |
| `batch_indices` | `(56,)` | `int64` | `[0, 1, ..., 55]` |
| `ntc_pb_ensemble_56` | `(56, G)` | `float64` | 4-fold average predicted `'non-targeting'` expression on `G_eval` in clipped `[0, 14]` space (for `L2_ntc_recon_match` primary & negative-L2 criteria) |
| `ntc_pb_per_fold_56` | `(4, 56, G)` | `float64` | Per-fold predicted `'non-targeting'` expression on `G_eval` (for Ruling 3 per-fold `L2` selection diagnostic) |
| `fold_names` | `(4,)` | `str` | `["zeroshot/hepg2", "zeroshot/jurkat", "zeroshot/k562", "zeroshot/rpe1"]` |
| `per_fold_r_orth_cal_56` | `(4, 56)` | `float64` | Secondary per-fold 56-index point `r_orth_cal` vector (`PREREG_J §1.2`) |
| `delta_pred_L1_ens56` | `(P, G)` | `float64` | Exact uniform mean `delta_pred_56.mean(axis=0)` |
| `delta_pred_L2_ntc_recon_match` | `(P, G)` | `float64` | Exact `delta_pred_56[l2_selected_batch_idx]` |
| `l2_selected_batch_idx` | `()` | `int64` | Index `b_NTC in 0..55` selected by primary Pearson `C_NTC(b)` |
| `delta_pred_bypass_zero` | `(P, G)` | `float64` | Secondary neutral-batch comparator (`batch_encoder.weight[0] = 0`) |
| `delta_pred_mean_batch_emb` | `(P, G)` | `float64` | Secondary neutral-batch comparator (`batch_encoder.weight[0] = mean_b W[b]`) |
| `delta_pred_multi_batch_ens8` | `(P, G)` | `float64` | Exact mean over `b in {0, 7, 14, 21, 28, 35, 42, 49}` |
| `delta_pred_medoid_idx` | `(P, G)` | `float64` | Exact `delta_pred_56[medoid_batch_idx]` |
| `medoid_batch_idx` | `()` | `int64` | Weight-only Euclidean geometric medoid index of 4-fold mean `batch_encoder.weight` |
| `target_names` | `(P,)` | `str` | Target names matching `observed.npz` |
| `gene_names` | `(G,)` | `str` | Gene names matching `observed.npz` |
| `model_key` | `()` | `str` | `STATE_SE` or `STATE_HVG` |
| `checkpoint_id` | `()` | `str` | `zeroshot/ensemble4` |
| `hvg_norm_mode` | `()` | `str` | `official_median_log1p` for `STATE_HVG`, `""` for `STATE_SE` |
| `control_source` | `()` | `str` | Target line key |
| `control_barcodes` | `(n_ctrl,)` | `str` | Ordered sampled control cell barcodes |

### 3.3 `preds_<MODEL>_g4_mismatched.npz` (`<MODEL>` in `{STATE_SE, STATE_HVG}`)
Contains all mismatched-control predicted deltas required for `J5a` (`PREREG_J §3.5`), covering both cross-fitted pairs (`(b*_1, b_dag_1)` evaluated on `P2`, plus optional `(b*_2, b_dag_2)` on `P1`) and pre-fixed continuity pairs (`(0, 17)` and `(0, 35)` on all `N_eval`):

| Key | Shape | Dtype | Description |
|---|---|---|---|
| `mismatched_sources` | `(S_mis,)` | `str` | Sorted union of `seen_sources`, `unseen_sources`, and `same_lineage_sources` |
| `seen_sources` | `(S_seen,)` | `str` | `seen_by_checkpoint` pool subset |
| `unseen_sources` | `(S_unseen,)` | `str` | `unseen_by_checkpoint` pool subset |
| `same_lineage_sources` | `(S_same,)` | `str` | Same-lineage comparator subset |
| `g4_batch_indices` | `(B_g4,)` | `int64` | Sorted unique batch indices (`sorted({0, 17, 35, b*_1, b_dag_1, [b*_2, b_dag_2]})`) |
| `delta_pred_mismatched` | `(S_mis, B_g4, P, G)` | `float64` | 4-fold mean predicted delta for mismatched source `s` at batch index `g4_batch_indices[k]` |
| `b_star_1` | `()` | `int64` | Best index selected strictly on `P1` of locked split `s = 0` |
| `b_dag_1` | `()` | `int64` | Worst index selected strictly on `P1` of locked split `s = 0` |
| `b_star_2` | `()` | `int64` | Best index selected strictly on `P2` of locked split `s = 0` |
| `b_dag_2` | `()` | `int64` | Worst index selected strictly on `P2` of locked split `s = 0` |
| `target_names` | `(P,)` | `str` | Target names matching `observed.npz` |
| `gene_names` | `(G,)` | `str` | Gene names matching `observed.npz` |
| `model_key` | `()` | `str` | `STATE_SE` or `STATE_HVG` |
| `checkpoint_id` | `()` | `str` | `zeroshot/ensemble4` |

### 3.4 `SENS_FCT_full_r9_6_incl_vcc2025_h1.npz` (`gse295214_h1` ONLY)
Mirrors Round 16 Addendum A Item 1 (fixed on primary `N_eval` and `G_eval`, updating `FCT`, `Transfer_Mean`, and `Source_Grand_Mean` using full `r9_6`):
- `fct_delta_pred` `(P, G)` `float64`
- `delta_transfer` `(P, G)` `float64`
- `source_grand_mean` `(G,)` `float64`
- `target_names` `(P,)` `str`
- `gene_names` `(G,)` `str`

### 3.5 `array_manifest_sha256.json`
JSON map `{filename: sha256_hex}` covering every `.npz` file in `arrays/<line>/`.

---

## 4. Independent Recomputation Guarantee (Stage 4)

Given ONLY `round17/results/arrays/`, running:
```bash
round16/.venv/bin/python round17/src/run_eval_r17.py --recompute-from-arrays round17/results/arrays
```
(or calling `recompute_cohort_from_exported_arrays_r17(Path("round17/results/arrays"))`) recomputes every `H-J1`, `H-J2` (including normal-calibrated & 20,000-bootstrap Bonferroni-56 companions), `H-J3`, `L2_ntc_recon_match` (including negative-L2 and per-fold selections), `J5a` cross-fitted & continuity `G4` 2×2 decompositions, `J5b` Spearman + 10,000-permutation null + BH FDR, separate-model Holm–Bonferroni, and pooled-model Holm–Bonferroni sensitivity with `max_abs_diff <= 1e-12` against `summary_r17.json` (well inside Stage 4's `1e-10` bound).
