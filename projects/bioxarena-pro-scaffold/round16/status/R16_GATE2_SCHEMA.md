# Round 16 Gate 2 Weight-Structure & Schema Audit (`R16_GATE2_SCHEMA.md`)

- **Date:** `2026-11-25 (+08:00)`
- **Author:** Main Line (`Antigravity`)
- **Directive:** [`controller/msgs/r16_MAIN_gate2.txt`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/controller/msgs/r16_MAIN_gate2.txt) + Master Controller supplemental directives (`2026-11-25 23:20` conservative `ci975 == 0.0` tie-precedence ruling; `2026-11-25 23:28` `SE-600M` 77/77 SHA-256 completion notice) + Ledger `G4` (`batch` vs. control-cell 2×2 decomposition draft)
- **Audited Artifacts (`round16/shared/models/`, `77` files, `R16_MODEL_CHECKSUMS.sha256` `77/77 OK`):**
  - `arcinstitute/ST-SE-Replogle` (`zeroshot/{hepg2,jurkat,k562,rpe1}`, 4 checkpoints, `202,210,763 B` each)
  - `arcinstitute/ST-HVG-Replogle` (`zeroshot/{hepg2,jurkat,k562,rpe1}`, 4 checkpoints, `164,071,711 B` each)
  - `arcinstitute/SE-600M` (`model.safetensors` `2,860,276,816 B`, `protein_embeddings.pt` `410,886,729 B`, `config.yaml`)
- **Isolated Environment:** `/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/.venv/bin/python` (`Python 3.11`, `arc-state 0.11.1`, `cell_load`, PyTorch MPS; executed with `WANDB_MODE=disabled HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1`).
- **Strict Blinding, Guard, & Non-Modification Attestation:**
  - **ZERO real cells inferred** (including zero real control cells from any `.h5ad`, `.h5`, `.mtx.gz`, or `.npz`).
  - **ZERO `.npz`, `.h5ad`, or Arc archived `adata_pred.h5ad` / `adata_real.h5ad` files opened** (all gene/target intersections computed strictly from plain-text `genes.txt`, `targets.txt`, `COMMON_GENES.txt`, and frozen `hgnc_complete_set.txt`).
  - **ZERO modifications** to `round16/shared/models/`, `round16/.venv/`, `round16/prereg/PREREG_I.md`, `round16/audit/`, `paper/`, Rounds 9–15 frozen files, or git `HEAD` (`80f7c17`).

---

## §1. `G_model`: Output Gene Vocabulary & `ST-SE` vs. `ST-HVG` Decoding Mechanism

### §1.1 Exact Output Gene Table (`|G_model| = 2,000`, Identical Across All 8 Checkpoints)
Direct read-only inspection across all 8 downloaded zero-shot checkpoints (`ST-SE-Replogle/zeroshot/{hepg2,jurkat,k562,rpe1}` and `ST-HVG-Replogle/zeroshot/{hepg2,jurkat,k562,rpe1}`) verified that the output gene list and exact ordering are **100% identical across all 8 checkpoints**:

- **Count:** **`|G_model| = 2,000`** unique gene symbols (`2,000` raw symbols; `2,000` unique canonical symbols under frozen R10 `HGNCBridge`).
- **First 5 genes (indices `0..4`):** `['HES4', 'ISG15', 'MIB2', 'CEP104', 'ACOT7']`
- **Last 5 genes (indices `1995..1999`):** `['MT-CO3', 'MT-ND3', 'MT-ND4L', 'MT-ND5', 'MT-ND6']`
- **Primary Source Files & Exact Locations:**
  - `round16/shared/models/ST-SE-Replogle/zeroshot/<line>/var_dims.pkl` (key `gene_names`, length `2,000`; scalar keys `input_dim=2058`, `output_dim=2058`, `gene_dim=6546`, `hvg_dim=2000`, `pert_dim=2024`, `batch_dim=56`)
  - `round16/shared/models/ST-HVG-Replogle/zeroshot/<line>/var_dims.pkl` (key `gene_names`, length `2,000`; scalar keys `input_dim=2000`, `output_dim=2000`, `gene_dim=6546`, `hvg_dim=2000`, `pert_dim=2024`, `batch_dim=56`)
  - Identically mirrored in `checkpoints/final.ckpt["hyper_parameters"]["gene_names"]` and `version_0/hparams.yaml` (`lines 12–2012`).
- **Why `2,000` instead of `6,546`:** Training config (`zeroshot/<line>/config.yaml:6`) sets `output_space: gene`, which selects `hvg_dim = 2000` rather than `gene_dim = 6546` (`state/_cli/_tx/_train.py:133–134`).

### §1.2 `ST-SE` Output Space & Exact Decoding Path
- **`ST-SE` backbone output (`preds`) lives in `2,058`-dim `X_state` embedding space**, NOT gene space:
  - `StateTransitionPerturbationModel.forward()` (`state/tx/models/state_transition.py:396–523`) encodes `2,058`-dim `ctrl_cell_emb` via `basal_encoder` (`Linear(2058,128)->...`), `2,024`-dim `pert_emb` via `pert_encoder`, and `batch` index via `batch_encoder` (`Embedding(56,128)`), runs an 8-layer bidirectional Llama backbone (`hidden_size=128`, `intermediate_size=3072`, `num_attention_heads=2`, `cell_set_len=32`, `use_rotary_embeddings=false`, `zeroshot/<line>/config.yaml:58–77`), and projects `res_pred + control_cells` (`state_transition.py:508–509`) through `project_out` to `[B*S, 2058]`.
- **Decoding to `2,000`-gene expression space (`LatentToGeneDecoder`):**
  - Because `finetune_vci_decoder: false` (`config.yaml:47`), checkpoint loading instantiates `self.gene_decoder = LatentToGeneDecoder(**decoder_cfg)` (`state/tx/models/base.py:15–116, 257–260`) with `decoder_cfg = {'latent_dim': 2058, 'gene_dim': 2000, 'hidden_dims': [1024, 1024, 512], 'dropout': 0.1, 'residual_decoder': False}`.
  - Architecture (`base.py:60–75`):
    `Linear(2058,1024) -> LayerNorm(1024) -> GELU -> Dropout(0.1) -> Linear(1024,1024) -> LayerNorm(1024) -> GELU -> Dropout(0.1) -> Linear(1024,512) -> LayerNorm(512) -> GELU -> Dropout(0.1) -> Linear(512,2000) -> ReLU`.
  - During `predict_step` (`state_transition.py:791–794`), `ST-SE` calls `self.gene_decoder(preds)` to produce `pert_cell_counts_preds` of shape `[N, 2000]`, which `_predict.py:526–537` and `_infer.py:838–849, 920–936` save as the gene-space prediction and clip to `[0.0, 14.0]`.
  - Weight-audit verification: all `14` `gene_decoder.*` tensors in `ST-SE-Replogle` are actively trained and fold-specific (`max_abs_diff` across folds `= 0.550–0.751`; LayerNorm weights deviate up to `0.333` from `1.0`, biases up to `0.417` from `0.0`).

### §1.3 Critical Gate-2 Discovery on `ST-HVG`: `gene_decoder` Is Untrained Initialization — Read `preds` Directly
- Although `ST-HVG-Replogle`'s `final.ckpt` `state_dict` also contains `14` `gene_decoder.*` tensors (`LatentToGeneDecoder(2000 -> 2000)`), **direct 4-fold weight inspection proves `ST-HVG`'s `gene_decoder` NEVER received gradients during training**:
  - `max_abs_diff(hepg2, jurkat) = 0.0`, `max_abs_diff(hepg2, k562) = 0.0`, `max_abs_diff(hepg2, rpe1) = 0.0` (exact bit-for-bit identical untouched seed-42 random initialization across all 4 folds);
  - LayerNorm weights are identically `1.0` (`max_dev = 0.0`) and biases identically `0.0` (`max_dev = 0.0`).
- **Source-code reason:** `PerturbationDataModule` (`cell_load/data_modules/perturbation_dataloader.py:204–210`), `_predict.py:350–354`, and `_infer.py:720–722, 838–849` all set **`store_raw_expression = False`** whenever `embed_key == "X_hvg"` and `output_space == "gene"`. Consequently, training batches omit `pert_cell_counts`, `gene_decoder` is never trained, and official inference uses `output_dict["preds"]` directly (which already has `ReLU` applied at `state_transition.py:514–516` and is clipped to `[0.0, 14.0]`).
- **Locked Adapter Rule:** `ST-SE` reads `clip(model.gene_decoder(preds), 0.0, 14.0)` (`2058 -> 2000`); `ST-HVG` reads `clip(preds, 0.0, 14.0)` directly (`2000`), **NEVER** `model.gene_decoder(preds)`.

---

## §2. `P_model`, Exact Intersections `N_eval^R16(t)`, and `|G_eval^R16(t)|` (Plain-Text Metadata Only; Zero `.npz` Opened)

### §2.1 Perturbation Vocabulary (`pert_onehot_map.pt`)
Verified 100% identical across all 8 checkpoints (`round16/shared/models/ST-{SE,HVG}-Replogle/zeroshot/<line>/pert_onehot_map.pt` and `var_dims.pkl["pert_names"]`):
- **Total keys:** **`2,024`** = **`2,023` non-control gene perturbations** + **`1` control token `'non-targeting'`** (`pert_dim = 2024`, float32 one-hot vectors).
- **First 5 keys (lexicographic):** `['AAAS', 'AAMP', 'AAR2', 'AARS', 'AARS2']`
- **Last 5 keys (lexicographic):** `['ZRANB2', 'ZRSR2', 'ZW10', 'ZWINT', 'non-targeting']`
- **Overlap with Round 9 source pools (from `targets.txt` only):** All `2,023 / 2,023` (`100.0%`) non-control keys in `P_model` appear inside `Replogle-4` (`2,391` union targets -> `2,023` in `P_model`) and `r9_6` (`2,655` union targets -> `2,023` in `P_model`).

### §2.2 Line-by-Line Gene & Target Intersections (Raw Exact Match vs. Frozen R10 `HGNCBridge`)
Computed strictly from `round*/shared/lines/<line>/{genes.txt,targets.txt}`, `round9/p2/shared/COMMON_GENES.txt` (`6,148` genes), `var_dims.pkl["gene_names"]` (`2,000`), `pert_onehot_map.pt` (`2,023` non-control), and frozen `/Users/mac/.cache/r9_p2_raw/hgnc_complete_set.txt` (`round10/src/data_loader_r10.py:64–154`). **Zero `.npz` or `.h5ad` files were opened.**

| Scope / Role | Line Key (`dir`) | `\|genes.txt\|` | `\|G ∩ G_model\|` (Raw / Canon) | **`\|G_eval^R16(t)\|`** `= \|G_R9 ∩ G_t ∩ G_model\|` (Raw / **Canon**) | `\|targets.txt\|` | **`\|targets.txt ∩ P_model\|`** (Raw / **Canon**) | `\|T ∩ P_model ∩ S_FCT(>=1)\|` (Raw / **Canon**) | HGNC-Bridged Target Examples (`+ΔT`) | Historical R9–R12 `P5` `N_eval` Reference (`R10/11/12_QC_COUNTS.md`, `PREREG_E`) |
|---|---|---:|---:|---:|---:|---:|---:|---|---:|
| **Group 1 Primary (CRISPRi, LOO)** | `nadig_hepg2` (`round9/p2/...`) | `9,432` | `1,989` / `1,989` | `1,843` / **`1,843`** | `2,343` | `1,996` / **`1,996`** | `1,996` / **`1,996`** | Same legacy namespace as `P_model` (`0` net gain) | `797` |
| | `nadig_jurkat` (`round9/p2/...`) | `8,826` | `1,995` / `1,995` | `1,843` / **`1,843`** | `2,350` | `1,999` / **`1,999`** | `1,999` / **`1,999`** | Same legacy namespace as `P_model` (`0` net gain) | `931` |
| **Group 2 Primary (CRISPRko)** | `a549_ko` (`round12/...`) | `17,337` | `1,939` / `1,985` | `1,797` / **`1,837`** | `999` | `498` / **`515`** | `498` / **`515`** | `+17`: `AARS1->AARS`, `CARS1`, `DARS1`, `EPRS1`, `GARS1`, `H2AX`, `HARS1`, `LARS1`, `MACIR`, ... | `445` |
| **Group 2 Primary (CRISPRi)** | `dld1` (`round10/...`) | `13,901` | `1,934` / `1,980` | `1,792` / **`1,832`** | `4,934` | `481` / **`490`** | `481` / **`490`** | `+9`: `COXFA4`, `DARS1`, `H2AX`, `LARS1`, `MACIR`, `NARS1`, `POLR1F`, `RMP64`, `TAFAZZIN` | `86` |
| | `cd4t_rest` (`round10/...`) | `11,756` | `1,814` / `1,854` | `1,790` / **`1,830`** | `11,903` | `1,387` / **`1,423`** | `1,387` / **`1,423`** | `+36`: `AFG2A`, `AFG2B`, `BBLN`, `CEP43`, `COXFA4`, `DARS1`, `H2AZ1`, ... | `302` |
| **Sensitivity Target / Retained Unseen `S_mis` (`Q4`)** | `feng_ipsc` (`round11/...`) | `7,595` | `1,642` / `1,681` | `1,540` / **`1,575`** | `444` | `282` / **`286`** | `282` / **`286`** | `+4`: `AARS1`, `GARS1`, `LARS1`, `WARS1` | `93` |
| **Sensitivity Only** | `cd4t_stim48` (`round11/...`) | `12,458` | `1,820` / `1,860` | `1,796` / `1,836` | `11,730` | `1,289` / `1,328` | `1,289` / `1,328` | `+39` | `503` |
| | `a549_ko_cas9pos` (`round12/...`) | `16,721` | `1,938` / `1,984` | `1,797` / `1,837` | `999` | `498` / `515` | `498` / `515` | `+17` | `445` |
| | `a549_myllia` (`round14/...`) | `15,456` | `1,886` / `1,927` | `1,796` / `1,836` | `212` | `50` / `50` | `50` / `50` | `0` | `32` |
| | `hek293_xucao` (`round12/...`) | `22,155` | `1,940` / `1,985` | `1,799` / `1,838` | `200` | `124` / `124` | `124` / `124` | `0` | `16` |
| **Optional Extensions (Not Local Raw)** | `replogle_k562_essential` (`round9/p2/...`) | `8,561` | `2,000` / `2,000` | `1,843` / `1,843` | `2,049` | `1,730` / `1,730` | `1,726` / `1,726` | Legacy namespace | `858` |
| | `replogle_rpe1` (`round9/p2/...`) | `8,746` | `2,000` / `2,000` | `1,843` / `1,843` | `2,351` | `1,999` / `1,999` | `1,999` / `1,999` | Legacy namespace | `1,068` |
| | `vcc2025_h1` (`round9/p2/...`) | `14,656` | `1,857` / `1,857` | `1,843` / `1,843` | `300` | `93` / `93` | `93` / `93` | `+1` (`TAFAZZIN`) | `80` |
| **Excluded (No Local Raw)** | `nourreddine_kolf21j` (`round9/p2/...`) | `15,000` | `2,000` / `2,000` | `1,843` / `1,843` | `95` | `23` / `23` | `23` / `23` | `0` | `18` |

### §2.3 Mandatory `PREREG_I v5` Disambiguations Revealed by Gate 2
1. **HGNC Canonical Bridging (Locked Same as Rounds 10–12):** `ST-SE`/`ST-HVG` (`G_model`, `P_model`) use legacy Replogle-4 HGNC symbols (`AARS`, `DARS`, `TAZ`, `H2AFZ`), identical to Round 9 `COMMON_GENES.txt`, whereas modern Group 2 lines (`a549_ko`, `dld1`, `cd4t_rest`, `feng_ipsc`) use updated HGNC symbols (`AARS1`, `DARS1`, `TAFAZZIN`, `H2AZ1`). Applying frozen R10 `HGNCBridge` (`round10/src/data_loader_r10.py:64–154`, already used by `FCT` in Rounds 10–12) harmonizes `STATE` and `FCT` to the exact same canonical symbol space (`|G_eval^R16| = 1,830–1,843` on all 5 primary lines; `+9` to `+36` targets recovered on Group 2).
2. **Exact Definition of $N_{\mathrm{eval}}^{\mathrm{R16}}(t)$ (Preventing Unintended Self-Gene Collapse):**
   - In Rounds 10–12 (`R10/R11/R12_QC_COUNTS.md`), `P5` target eligibility required: (1) reproducible on target line $t$ (`repr_mask` inside `pseudobulk_delta.npz`, not opened at Gate 2), (2) reproducible on $\ge 1$ source line in $\mathcal{S}_{\mathrm{FCT}}(t)$, and (3) self-gene $g_p \in G_{\mathrm{eval}}^{\mathrm{R9}}(t)$ (`~6,124` genes).
   - **Critical Gate-2 Catch:** Only `~18–22%` of perturbed genes happen to be themselves among Arc's `2,000` output HVGs (`G_model`) — e.g., on `a549_ko`, only `91` of `515` `P_model`-predictable targets have $g_p \in G_{\mathrm{eval}}^{\mathrm{R16}}$, and on `feng_ipsc` only `47 / 286`. Requiring $g_p \in G_{\mathrm{eval}}^{\mathrm{R16}}(t)$ on top of `P5` reproducibility would collapse `dld1` from `~73–86` down to `~15` targets!
   - **Recommended `PREREG_I v5` Lock:** Define $N_{\mathrm{eval}}^{\mathrm{R16}}(t)$ strictly as the intersection of the **frozen R9–R12 `P5`-qualified reproducible target set** (whose exact count is `797 / 931 / 445 / 86 / 302` on primary lines and `93` on `feng_ipsc`) with **`canon(P_model)`**, masking out $g_p$ from $G_{\mathrm{eval}}^{\mathrm{R16}}(t)$ whenever $g_p \in G_{\mathrm{eval}}^{\mathrm{R16}}(t)$ (`all_nontarget` = $G_{\mathrm{eval}}^{\mathrm{R16}}(t) \setminus \{g_p\}$, `target_indices[i] = -1` when $g_p \notin G_{\mathrm{eval}}^{\mathrm{R16}}(t)$). Because `P_model` covers `2,023 / 2,391` (`84.6%`) of Replogle-4 targets (`|targets.txt ∩ P_model|` = `1,996 / 1,999 / 515 / 490 / 1,423`), all 3 Group 2 primary lines (`a549_ko`, `dld1`, `cd4t_rest`) remain cleanly above `>= 50` ($\mathcal{G}_2^{\mathrm{qual}} = \{\text{a549\_ko}, \text{dld1}, \text{cd4t\_rest}\}$, $K_2 = 3 \ge 2$, degenerate Group-2 safeguard not triggered). Exact post-`P5` counts $|P5(t) \cap \mathrm{canon}(P_{\mathrm{model}})|$ are logged mechanically at Gate 5 when `.npz` is first unlocked.

---

## §3. `batch` & `cell_type` Covariate Mechanisms, Symmetric Neutral Scheme, & `G4` 2×2 Decomposition Draft

### §3.1 Exact Source-Code Behavior (`arc-state 0.11.1` & `cell_load`)
1. **`cell_type_onehot_map.torch` (`{'hepg2': 0, 'jurkat': 1, 'k562': 2, 'rpe1': 3}`):**
   - **NOT required by and NEVER passed into `StateTransitionPerturbationModel`!** Confirmed in `state/tx/models/state_transition.py:160–524` and all 8 `final.ckpt` `state_dict`s (`0` `cell_type_encoder` parameters).
   - Used exclusively by `PerturbationDataModule` (`cell_load/data_modules/perturbation_dataloader.py:687–713`) and `_infer.py:739–754` to partition cells when sampling control pools and aggregating `cell-eval` outputs.
2. **`batch_onehot_map.torch` & `batch_encoder = nn.Embedding(56, 128)` (`batch_encoder: true`, `use_batch_token: false`, `batch_col: gem_group`):**
   - **MANDATORY at inference:** `state/tx/models/state_transition.py:192–197, 425–441` requires `batch["batch"]` whenever `self.batch_encoder is not None`:
     ```python
     # state_transition.py:422-441
     seq_input = pert_embedding + control_cells
     if self.batch_encoder is not None:
         batch_indices = batch["batch"]             # KeyError if omitted!
         ...
         batch_embeddings = self.batch_encoder(batch_indices.long())  # [B, S, 128]
         seq_input = seq_input + batch_embeddings
     ```
   - `batch_onehot_map.torch` (identical across all 8 checkpoints) maps `56` string `gem_group` labels `'1'`..`'56'` sorted lexicographically (`'1'->0`, `'10'->1`, `'11'->2`, ..., `'56'->51`, `'6'->52`, `'7'->53`, `'8'->54`, `'9'->55`).
   - 4-fold weight inspection confirms **all 56 rows received training gradients across every fold** (`0` identical rows across any pair of folds, `min_row_diff = 0.043–0.061`, `0` fold-specific frozen rows), because Replogle/Nadig datasets each number their 10x channels starting from `'1'..'K'`.
3. **How Official `arc-state 0.11.1` Handles Unseen Cell Lines / Unseen Batches (`state/_cli/_tx/_infer.py:656–698, 822–825`):**
   - **Case A (a batch column exists, e.g. `batch` / `donor` / `lane`, with unseen labels):** `_infer.py:683–693` silently falls back to **`idxs[i] = 0`** (`batch_onehot_map` key `'1'`, `batch_encoder.weight[0]`) with a printed warning.
   - **Case B (no batch column found in `adata.obs`):** `_infer.py:694–697, 822–825` sets `uses_batch_encoder = False`, passes `bi = None` to `prepare_batch`, omits `batch["batch"]`, and **crashes with `KeyError: 'batch'` at `state_transition.py:427`!**
   - **Confounding Hazard Without Explicit Neutralization:** If native `.obs` batch strings were passed during matched-vs-mismatched testing, `nadig_hepg2`/`nadig_jurkat` (or any dataset with numeric lane IDs `'1'..'56'`) would index into `batch_encoder.weight[0..55]` while datasets with string batch IDs (e.g., `CE0006864_D4`, `batch1`) would all map to index `0` — violating matched=mismatched covariate symmetry!

### §3.2 Quantitative Weight & Synthetic Sensitivity Discovery (`batch_encoder` Dominates `seq_input`!)
Weight-norm audit + end-to-end pure-synthetic `SE-600M -> ST-SE` smoke test (`§6`) uncovered a major architectural mechanism directly relevant to Ledger `G4` and Zhang et al. (`arXiv:2611.05849`):
- Because `nn.Embedding(56, 128)` is initialized from $\mathcal{N}(0, 1)$ (expected L2 norm $\sqrt{128} = 11.31$), **every trained row of `batch_encoder.weight` has L2 norm `9.77–12.43` (median `11.09`; row `0` = `10.8311` in `ST-SE/hepg2`, `10.7867` in `ST-HVG/hepg2`)**, whereas:
  - `||basal_encoder(X_state)||_2 = 0.38–1.20` (on real `SE-600M`-encoded synthetic cells);
  - `||pert_encoder(pert)||_2 = 0.7371`;
  - `||mean_b(batch_encoder.weight)||_2 = 1.5453` (`ST-SE`) / `1.4661` (`ST-HVG`).
- Because `seq_input = pert_embedding + control_cells + batch_embeddings` (`state_transition.py:422–441`) performs unnormalized addition prior to the 8-layer Llama backbone, **the `batch_encoder` vector accounts for `>98%` of input L2 energy**!
- **Empirical End-to-End Synthetic 2×2 Sensitivity (`SE-600M -> ST-SE`, `20` perturbations + `'non-targeting'`, `512` synthetic cells, MPS `fp32`):**
  - **Swapping synthetic `X_ctrl` ONLY (`X_A` vs. `X_B`, fixed `b = 0`):** `l2_abs = 5.49`, `l2_delta_ntc_diff = 4.25`, `l2_pert_centered_diff = 2.63`, **`corr_delta_ntc = +0.9265`**, **`corr_pert_centered = +0.9283`** (`M5.3` PASS: strictly non-zero control sensitivity).
  - **Swapping ONLY `batch_idx` (`b = 0` vs. `b = 1, 17, 35`, fixed `X_A`):** `l2_abs = 19.48–22.18` (**3.6×–5.7× larger**), `l2_delta_ntc_diff = 13.65–17.92`, `l2_pert_centered_diff = 9.94–12.17`, and **`corr_delta_ntc = +0.037 to +0.116`**, **`corr_pert_centered = +0.049 to +0.123` (near-orthogonal!)**.
  - **Replacing `b = 0` with `zero` (`0 ∈ R^128`) or `mean` (`norm 1.55`):** Shifts `seq_input` `~10` L2 units out of the training norm shell (`[9.77, 12.43]`), producing `corr_delta_ntc(b=0 vs zero) = 0.146` and `corr_delta_ntc(b=0 vs mean) = 0.166`, while still passing `M5.3` (`ctrl_swap` `l2_abs = 6.07–6.88`, `corr_delta_ntc = 0.739–0.845`).

### §3.3 Locked Completely Symmetric Neutral Covariate Scheme (Matched = Mismatched)
- **Primary Scheme for `H-I1` & `H-I2` (Locked Identically Across Matched & Mismatched Arms in Both Group 1 and Group 2):**
  - **`fixed_index_zero` (`batch_idx = 0`, corresponding to `batch_onehot_map['1']`)** passed identically to every cell in both the matched control arm and every mismatched control arm across all Group 1 and Group 2 lines, with zero modification to frozen `Arc STATE` weights.
  - **Rationale:** (1) Matches official `arc-state 0.11.1` (`state/_cli/_tx/_infer.py:688`) unseen-batch fallback verbatim; (2) keeps `seq_input` strictly inside the trained Llama input norm shell (`||e_0||_2 = 10.83` vs. training median `11.09`), unlike `bypass_zero` (`norm 0`) or `mean_batch_emb` (`norm 1.55`); (3) preserves **100% covariate symmetry** (`b_matched == b_mismatched == 0`) with zero batch-identity leakage; (4) verified to pass `M5.3` non-zero `X_ctrl` sensitivity on real `SE-600M -> ST-SE` (`l2_abs = 5.49`, `max_abs = 0.973`, `l2_delta_ntc_diff = 4.25`) and `ST-HVG` (`l2_abs = 8.92`, `l2_delta_ntc_diff = 6.81`).

### §3.4 Pre-Registered Secondary Analysis Draft for `PREREG_I v5 §4.4` (Ledger `G4`: `batch` vs. Control-Cell Decomposition)
Because `batch_encoder` exerts a multi-fold larger input perturbation than `basal_encoder` in `Arc STATE`, `PREREG_I v5` adds a pre-registered **secondary descriptive analysis (`§4.4`, non-confirmatory)** with two components:
1. **Neutral-Batch Robustness Panel (Addressing Single-Index Dependence of `b = 0`):** Re-evaluate `H-I1` and `H-I2` under three alternative completely symmetric neutral schemes:
   - `(a) bypass_zero`: override `batch_embeddings = 0 ∈ R^128` identically for matched & mismatched arms (matches `A-T3`'s `bypass_zero` variant in `R16_CONTROLLER_NOTES.md`);
   - `(b) mean_batch_emb`: override `batch_embeddings = mean_{b=0..55}(batch_encoder.weight[b]) ∈ R^128`;
   - `(c) multi_batch_ens8`: run inference independently under `K = 8` fixed trained indices `b ∈ {0, 7, 14, 21, 28, 35, 42, 49}` (each individually in-distribution with norm `~11`) and average predictions in clipped decoded gene space before computing $\hat{\Delta}_{\mathrm{pred}}(p)$.
2. **2×2 Factorial `G4` Decomposition ("只换 batch 不换对照" vs. "只换对照不换 batch"):**
   Hold frozen `STATE_SE` weights and `pert_emb` fixed, and evaluate all 4 cells of the `(X_ctrl, batch)` factorial on:
   - **Group 1 (`nadig_hepg2`, `nadig_jurkat`):** Using each sampled control cell's native Replogle/Nadig `gem_group` index (`b_matched` from target line $t$'s `gem_group`s in `'1'..'56'`; `b_mismatched` from counterpart line $s'$'s `gem_group`s in `'1'..'56'`), noting that cell `(M, B_M)` matches native `state tx predict` for §1.3 `adata_pred.h5ad` parity verification.
   - **Group 2 (`a549_ko`, `dld1`, `cd4t_rest`, where native labels are unseen):** Using fixed paired trained indices `(b_1 = 0, b_2 = 17)` (and reporting spread across `{0, 1, 17, 35}`) alongside matched vs. mismatched `X_ctrl`.
   - **Exact Decomposition Metrics (reported under both `r_orth_cal` and `r_dorth_cal`, with `4,000` joint perturbation bootstrap `97.5%` CIs):**
     - `(M, B_M)`: `(X_ctrl = matched, b = b_matched)`
     - `(M, B_S)`: `(X_ctrl = matched, b = b_mismatched)` -> **Pure batch-swap effect ("只换 batch 不换对照"):** $\Delta_{\mathrm{batch}} = r(M, B_M) - r(M, B_S)$
     - `(S, B_M)`: `(X_ctrl = mismatched, b = b_matched)` -> **Pure control-swap effect ("只换对照不换 batch"):** $\Delta_{\mathrm{ctrl}} = r(M, B_M) - r(S, B_M)$
     - `(S, B_S)`: `(X_ctrl = mismatched, b = b_mismatched)` -> **Joint swap & interaction:** $\Delta_{\mathrm{both}} = r(M, B_M) - r(S, B_S)$, $\Delta_{\mathrm{int}} = \Delta_{\mathrm{both}} - \Delta_{\mathrm{ctrl}} - \Delta_{\mathrm{batch}}$, plus ratio $|\Delta_{\mathrm{batch}}| / (|\Delta_{\mathrm{ctrl}}| + |\Delta_{\mathrm{batch}}|)$.
   - Directly links Round 16's input-control audit to Zhang et al. (`arXiv:2611.05849` §C.5: `control + batch` ridge probe PDC `0.27–0.46`) by isolating how much of `STATE`'s context specificity comes from `X_ctrl` vs. `batch_encoder`.

---

## §4. Input/Output Numerical Spaces, `predict_residual=true`, `SE-600M` Gene Coverage, & Conversion to `delta_obs`

### §4.1 `SE-600M` Input & Output Space (`model.safetensors`, `protein_embeddings.pt`, `config.yaml`)
- **Required Input:** **Raw non-negative integer UMI counts** on the full measured transcriptome (`adata.X` CSR matrix, gene symbols in `var_names` / `var["gene_name"]`), mapped by exact gene-symbol lookup to `SE-600M/protein_embeddings.pt` (`19,790` genes × `5,120`-dim `float32` ESM2 embeddings; first 5 `['A1BG', 'A1CF', 'A2M', 'A2ML1', 'A3GALT2']`, last 5 `['ZYG11A', 'ZYG11B', 'ZYX', 'ZZEF1', 'ZZZ3']`).
- **Internal Normalization (`state/emb/data/loader.py:17–19, 217–240, 431–473, 524–554` & `state/emb/nn/model.py:117–119, 295–332`):**
  - `is_raw_integer_counts()` detects raw counts (`max > 35` or `expm1.sum > 5,000,000`, `loader.py:17–19, 431–454`) and applies **`torch.log1p(raw_counts)` directly** (WITHOUT `libnorm(1e4)`), ranks the top `pad_length - 1 = 2,047` expressed mapped genes (`SE-600M/config.yaml:27`), computes relative expression `100 * (counts / sum(counts))` (`loader.py:549`) fed into `count_encoder` + `bin_encoder` (`model.py:127–133, 296–320`), and runs a 16-layer Transformer (`emsize=2048`, `nhead=16`, `d_hid=2048`).
- **Output (`X_state`, shape `[N, 2058]`, `float32`):**
  - Concatenation `[emb_2048 || ds_emb_10]` (`state/emb/inference.py:178–184, 249–254`):
    - `emb_2048`: L2-normalized (`||emb_2048||_2 = 1.0000`, verified std `0.00000`) CLS cell embedding (`model.py:324–331`);
    - `ds_emb_10`: 10-dim linear projection `dataset_embedder(ds_emb)` (`model.py:118, 155–157`; verified norm `~0.44`, range `[-0.319, +0.207]` on synthetic raw counts).
- **Critical `arc-state 0.11.1` Safetensors & Apple Silicon MPS Findings:**
  1. Official CLI `state emb transform` (`state/_cli/_emb/_transform.py:89`, `state/emb/inference.py:125`) globs `*.ckpt` and calls `StateEmbeddingModel.load_from_checkpoint()` (`torch.load`), which **crashes with `_pickle.UnpicklingError`** on HuggingFace's `SE-600M/model.safetensors` (`converted_utc: 2026-02-24`, `238` `F32` tensors, `715,062,442` parameters).
  2. Direct instantiation `StateEmbeddingModel(...)` (`state/emb/nn/model.py:61–172`) + `safetensors.torch.load_file("model.safetensors")` loads cleanly with **`missing_keys = []`** (`0` missing keys!) and sole `unexpected_key = ['pe_embedding.weight']`, which we verified is **bit-for-bit identical (`max_abs_diff = 0.0`)** to `torch.vstack(list(protein_embeddings.pt.values()))` (`[19790, 5120]`).
  3. `state/emb/inference.py:130, 174, 208` hardcodes `cuda` vs. `cpu` (ignoring Apple MPS); on Apple M3 Pro, both `MPS bfloat16` and `MPS float16` abort inside macOS Metal (`MPSNDArrayMatrixMultiplication.mm:5799: Destination NDArray and Accumulator NDArray cannot have different datatype`), whereas **`MPS float32`** executes cleanly (`all_finite = True`, **`4.0×` faster than CPU `fp32`**, mean cosine parity vs. CPU `fp32` `= 0.999416`, `max_abs = 0.0123`).
- **`SE-600M` `19,790`-Gene Mapping Coverage (Lower Bound from `genes.txt` Only; Raw Full-Transcriptome Files Have Even Higher Coverage, Verified at Gate 5):**

| Line / Set | `\|genes.txt\|` | Exact Raw Match in `protein_embeddings.pt` (`19,790`) | `%` of `19,790` / `%` of Line | Canonical HGNC Match |
|---|---:|---:|---:|---:|
| `G_model` (`2,000` HVGs) | `2,000` | `1,911` | `9.66%` / `95.55%` | `1,935` |
| `R9_COMMON_6148` | `6,148` | `6,068` | `30.66%` / `98.70%` | `6,141` |
| `nadig_hepg2` (Group 1) | `9,432` | **`9,022`** | `45.59%` / `95.65%` | `9,122` |
| `nadig_jurkat` (Group 1) | `8,826` | **`8,391`** | `42.40%` / `95.07%` | `8,493` |
| `a549_ko` (Group 2) | `17,337` | **`13,766`** | `69.56%` / `79.40%` | `14,000` |
| `dld1` (Group 2) | `13,901` | **`11,714`** | `59.19%` / `84.27%` | `11,902` |
| `cd4t_rest` (Group 2) | `11,756` | **`11,535`** | `58.29%` / `98.12%` | `11,723` |
| `feng_ipsc` (`Q4` Sens / `S_mis`) | `7,595` | **`7,360`** | `37.19%` / `96.91%` | `7,452` |
| `cd4t_stim48` / `a549_ko_cas9pos` | `12,458` / `16,721` | `12,222` / `13,450` | `61.76%` / `67.96%` | `12,423` / `13,675` |
| `a549_myllia` / `hek293_xucao` | `15,456` / `22,155` | `14,801` / `14,419` | `74.79%` / `72.86%` | `15,045` / `14,653` |

Every primary and mismatched-source line maps `7,360–13,766` genes into `SE-600M`'s `19,790`-gene ESM2 table (far exceeding `pad_length - 1 = 2,047` top genes per cell), confirming zero gene-ID mapping collapse.

### §4.2 `ST-HVG` Input Space, `predict_residual=true` Semantics, & Exact Conversion to `delta_obs`
- **`ST-HVG` Input (`X_hvg`, shape `[N, 2000]`):** Aligned strictly to the ordered `2,000` `G_model` genes (`data_module.torch`: `is_log1p=True`, `exp_counts=False`, `normalize_counts=False`; `state/_cli/_tx/_preprocess_train.py:45–55` applies `sc.pp.normalize_total(adata)` + `sc.pp.log1p(adata)`). Missing `G_model` genes in a given line (e.g., `15–146` genes out of `2,000`, `§2.2`) are zero-filled on input and excluded from evaluation since evaluation restricts strictly to $G_{\mathrm{eval}}^{\mathrm{R16}}(t) = G_{\mathrm{R9\_6148}} \cap G_t^{\mathrm{fin}} \cap G_{\mathrm{model}}$.
  - *Library-size normalization note:* `_preprocess_train.py:46` uses `sc.pp.normalize_total(adata)` (default `target_sum=None`, per-dataset median total count) rather than `target_sum=1e4`. For `ST-HVG` primary input, we lock `libnorm(1e4) -> log1p` (exact consistency across lines with differing sequencing depths and matching `delta_obs`), and verify sensitivity against `median-norm -> log1p` during the §1.3 Group 1 `adata_pred.h5ad` check.
- **Exact `predict_residual=true` Semantics (`state/tx/models/state_transition.py:501–516`):**
  - When `output_space == "gene"` (true for ALL 8 `ST-SE` and `ST-HVG` checkpoints), `predict_residual=true` executes **line 509** (NOT line 506!):
    ```python
    # state_transition.py:508-509
    elif self.predict_residual:
        out_pred = self.project_out(res_pred + control_cells)
    ```
    where `control_cells = self.basal_encoder(basal)` is the `128`-dim hidden encoding (NOT raw input `basal`), and `res_pred` is the `[B, 32, 128]` Llama backbone output. Thus the residual skip connection occurs in `128`-dim hidden space prior to `project_out`.
  - **Batch-window parity (`state_transition.py:409–438`):** Passing `512` cells with `padded=True` (reshaping into `(16, 32, D)`) is bit-for-bit identical (`max_abs = 0.0` on CPU, `1.85e-6` on MPS) to looping over 16 disjoint `32`-cell windows with `padded=False` (as in `_infer.py:815–835`). **CRITICAL:** Passing `512` cells with `padded=False` silently creates a single length-512 sequence (`(1, 512, D)`) instead of `16 × 32`; our adapter MUST lock `padded=True` on exact multiples of `32` (`16 × 32`).
- **Output Space & Exact Conversion to `delta_obs` (`libnorm(1e4) -> log1p` shift):**
  - Both `ST-SE` (`clip(gene_decoder(preds), 0, 14)`) and `ST-HVG` (`clip(preds, 0, 14)`) directly output **`log1p`-normalized gene expression** (`is_log1p=True`, `_predict.py:535, 584`, `_infer.py:920–936`) — NO additional `log1p` or `expm1` transform is applied!
  - **Locked Perturbation Shift Formula (`delta_pred(p)` on $G_{\mathrm{eval}}^{\mathrm{R16}}(t)$):**
    $$\hat{\Delta}_{\mathrm{pred}}(p) = \frac{1}{n_{\mathrm{ctrl}}}\sum_{i=1}^{n_{\mathrm{ctrl}}} \operatorname{clip}_{[0,14]}\!\big(\hat{Y}(x_i, p, b)\big) - \frac{1}{n_{\mathrm{ctrl}}}\sum_{i=1}^{n_{\mathrm{ctrl}}} \operatorname{clip}_{[0,14]}\!\big(\hat{Y}(x_i, \text{"non-targeting"}, b)\big)$$
    projected onto the exact ordered $G_{\mathrm{eval}}^{\mathrm{R16}}(t)$ gene columns.
  - **Why subtracting predicted `"non-targeting"` (and NEVER raw $\overline{X}_{\mathrm{ctrl}}$) is mathematically required:**
    1. For `ST-SE`, `X_ctrl` lives in `2,058`-dim `X_state` space while decoded outputs live in `2,000`-gene space; subtracting predicted `"non-targeting"` keeps both terms inside the exact same `LatentToGeneDecoder` output space and matches official `state tx predict` / `cell-eval` (`_predict.py:576–597`, which contrasts predicted perturbations against predicted `control_pert = "non-targeting"`).
    2. For `ST-HVG`, `ST-HVG`'s `project_out` outputs Replogle-scale `log1p` baselines (`l2(pred_ntc - X_ctrl) = 58.2` vs. `l2(pred_p - pred_ntc) = 6.7` on synthetic counts); subtracting `pred("non-targeting")` cancels the decoder's fixed baseline offset cleanly.
    3. Moreover, `r_orth_cal` (pooled Pearson correlation after orthogonalization to $\bar{\Delta}^{\mathrm{src}}_p$) is invariant to any uniform positive rescaling between median-normalized and `1e4`-normalized `log1p` units.

---

## §5. `SE-600M` Pretraining Data Audit & Disclosure Plan

### §5.1 Verified Facts From Local Weights, Configs, & Source Code
1. **`round16/shared/models/SE-600M/config.yaml`:**
   - `dataset.current: scbasecamp-cellxgene-tahoe-filtered` (`lines 20, 36–42`), `num_datasets: 14420`, `experiment.name: vci_1.5.0_600M_basecount_tahoe_cxg` (`line 76`) — combining **scBaseCamp (`basecount`), CELLxGENE, and Tahoe-100M**.
   - `dataset.num_cells: 36238464` (`line 22`): used in `state/emb/train/trainer.py:41–42` solely to compute `EPOCH_LENGTH` / `warmup_steps` scheduling (Adduri et al., *Cell* 2026 abstract reports ~167M observational cells; `PRIOR_ART_I.md §4.1`).
   - Validation paths explicitly referenced in `config.yaml`: `rpe1_top5000_variable.h5ad` (`lines 26, 128`) and `replogle_perturbation.h5ad` (`line 138`, `validations.*.enable: false`).
   - `embeddings.current: esm2-cellxgene-basecamp-tahoe` (`lines 45, 52–58`): `Homo_sapiens.GRCh38.gene_symbol_to_embedding_ESM2.pt` (`num: 19790`, `size: 5120`).
2. **`round16/shared/models/SE-600M/model.safetensors` header metadata:**
   - `source_checkpoint: /home/aadduri/SE-600M/se600m_epoch16.ckpt`, `converted_utc: 2026-02-24T01:55:29.880997+00:00` (HuggingFace repo created `2025-06-23`, original `se600m_epoch16.ckpt` uploaded `2025-07-09`; `PRIOR_ART_I.md §4.1`).
3. **`round16/shared/models/SE-600M/README.md` (`lines 1–37`):** Contains CLI usage and license references only; lists no dataset manifest.

### §5.2 Items That Cannot Be Determined (`未能确定`) & Disclosure / Handling Plan
- **`未能确定` (Cannot Be Determined):** The exact sample/study manifest CSVs referenced in `SE-600M/config.yaml:41–42` (`/large_storage/ctc/userspace/aadduri/data/auxillary/esm_basecount_tahoe_cellxgene_{train,val}_filtered.csv`) reside on Arc's internal cluster and are NOT bundled in `SE-600M` or `arc-state`. Therefore, exact inclusion/exclusion of any dataset published prior to `2025-06` in `SE-600M`'s `14,420`-dataset training set is **`未能确定`**.
- **Dataset-by-Dataset Exposure & Locked Handling Plan (Integrating `PRIOR_ART_I.md §4.2` GEO/ENA Dates):**
  - **Group 1 (`nadig_hepg2`, `nadig_jurkat`, `GSE264667`, GEO public `2024-05-05`):** Predates `SE-600M` (`2025-06`); exact membership in `14,420` CSV is **`未能确定`** -> **Explicitly disclosed as potentially present in `SE-600M` pretraining** (while strictly held out during `ST-SE-Replogle` `zeroshot/{hepg2,jurkat}` transition training).
  - **`feng_ipsc` (`ENA PRJEB81502`, `first_public = 2024-11-08`):** Predates `SE-600M` (`2025-06`); exact membership is **`未能确定`** -> **Remains demoted to `SENSITIVITY_ONLY` as a target line** (`Q4`), while retained as an `unseen_by_checkpoint` mismatched-control source on unrelated targets.
  - **All 3 Group 2 Primary Target Lines — Strictly Chronologically Isolated From `SE-600M` (`2025-06/07`) and `ST-*-Replogle`:**
    - `cd4t_rest`: `GSE314342` (GEO public **`2026-02-06`**) — strictly post-dates `SE-600M` (zero pretraining overlap possible).
    - `dld1`: `GSE337988` (GEO public **`2026-07-10`**) — strictly post-dates `SE-600M` (zero pretraining overlap possible).
    - `a549_ko`: `GSE345058` (GEO public **`2026-08-26`**, ENA `2026-07-03`) — strictly post-dates `SE-600M` (zero pretraining overlap possible).
  - **Conclusion:** Gate 2 adds **zero new pretraining disclosures** beyond `PREREG_I v4` (`Q4`), and confirms all 3 Group 2 primary target lines (`a549_ko`, `dld1`, `cd4t_rest`) strictly post-date `SE-600M`.

---

## §6. Complete End-to-End Pure-Synthetic Smoke Test (Zero Real Cells) & Wall-Time Estimate

### §6.1 Synthetic Setup (100% Synthetic Poisson Counts, `seed = 20261601`)
- **Zero real cells used.** Generated two distinct synthetic raw integer UMI count matrices (`raw_A`, `raw_B`, `n_cells = 512`, `mean_UMI ≈ 14,041–14,814` across `4,000` genes drawn from `protein_embeddings.pt` for `SE-600M -> ST-SE`; and `2,000` genes for `ST-HVG`, with distinct expression profiles between Profile `A` and Profile `B`).
- Executed full pipelines on **Apple M3 Pro (`36 GB` unified RAM, `MPS` + `CPU` parity checks)**:
  1. **End-to-End `SE-600M (model.safetensors)` -> `ST-SE-Replogle (zeroshot/hepg2)` -> `LatentToGeneDecoder`** (`20` perturbations + `'non-targeting'`, `16 × 32` `padded=True`, across `b ∈ {0, 1, 17, 35, zero, mean}`);
  2. **End-to-End `ST-HVG-Replogle (zeroshot/hepg2)`** (`libnorm(1e4) -> log1p` -> `16 × 32` `padded=True`, verifying bit-exact `0.0` parity vs. 16 loops of `32`-cell `padded=False` windows).

### §6.2 Measured Performance, Memory, Shapes, Parity, & `M5.3` / `G4` Smoke Results

| Stage / Model | Device & Precision | Load Time | Inference Latency | Output Tensor Shapes | Peak RSS | Key Verification Result |
|---|---|---:|---:|---|---:|---|
| **`SE-600M`** (`715M` params, `model.safetensors`) | **MPS `float32`** (`bs=16`) | `1.68 s` (+ `0.55 s` PE) | **`330.6–331.7 s` / 512 cells** (`646.8 ms/cell`) | `X_state`: **`(512, 2058)`** (`[2048 || 10]`) | **`5,695.6 MB`** (`5.56 GiB`) | `all_finite=True`; `||emb2048||_2 = 1.0000` (std `0.0`); `||ds10||_2 = 0.4409`; `cos(A,B) = 0.3637`, `L2(A,B) = 1.1390` |
| **`SE-600M`** (Parity Benchmark) | CPU `float32` (`bs=16`) | — | `20.66 s` / 8 cells (`2,583 ms/cell`, est. `22.0 min` / 512 cells) | `(8, 2058)` | — | MPS `fp32` vs. CPU `fp32`: **`mean_cos = 0.999416`**, `max_abs = 0.01228` (`4.0×` MPS speedup) |
| **`SE-600M`** (Precision Guard) | MPS `bf16` & `fp16` | — | **Aborts (`rc=-6`)** | — | — | macOS Metal `MPSNDArrayMatrixMultiplication.mm:5799` assertion -> **Lock `MPS float32`!** |
| **`ST-SE` + `LatentToGeneDecoder`** (`zeroshot/hepg2`) | **MPS `float32`** (`16×32`, `padded=True`) | `0.31 s` | **`0.266 s` / 21 perts** (**`12.65 ms/pert`** for 512 cells) | Latent `(512, 2058)` -> Decoded **`(512, 2000)`** | `2,348.5 MB` (standalone) / `5,695.6 MB` (e2e) | `out ∈ [0.0, 5.91]`, `frac>14 = 0.0`; `padded=True` vs. 16×`32` windows `max_abs = 1.85e-6` (`0.0` CPU); **`M5.3` PASS** (`ctrl_swap b=0`: `l2_abs=5.49`, `max_abs=0.973`, `l2_delta=4.25`) |
| **`ST-SE` + `LatentToGeneDecoder`** | CPU `float32` (`16×32`, `padded=True`) | `0.30 s` | `38.09 ms/pert` (512 cells) | `(512, 2058)` -> `(512, 2000)` | `2,348.5 MB` | `padded=True` vs. 16×`32` windows `max_abs = 0.0` (exact) |
| **`ST-HVG`** (`zeroshot/hepg2`, direct `preds`) | **MPS `float32`** (`16×32`, `padded=True`) | `0.34 s` | **`10.21 ms/pert`** (512 cells) | **`(512, 2000)`** | `2,348.5 MB` | `out ∈ [0.0, 7.26]`; **`M5.3` PASS** (`ctrl_swap b=0`: `l2_abs=8.92`, `l2_delta=6.81`, `corr_pert_cen=0.789`) |
| **`ST-HVG`** | CPU `float32` (`16×32`, `padded=True`) | `0.33 s` | `31.56 ms/pert` (512 cells) | `(512, 2000)` | `2,348.5 MB` | `padded=True` vs. 16×`32` windows `max_abs = 0.0` (exact) |

### §6.3 End-to-End Synthetic `M5.3` & `G4` Contrast Summary (`SE-600M -> ST-SE`, `MPS fp32`)

| Contrast (`20` perturbations + `'non-targeting'`, `512` synthetic cells) | `L2` Pseudobulk Diff | Max Abs Diff | `L2` `Δ_ntc` Diff | `L2` Pert-Centered Diff | Corr `Δ_ntc` | Corr Pert-Centered |
|---|---:|---:|---:|---:|---:|---:|
| **Control swap (`A` vs. `B`) under locked primary `b = 0`** (`M5.3`) | **`5.492`** | **`0.973`** | **`4.255`** | **`2.627`** | `+0.9265` | `+0.9283` |
| Control swap (`A` vs. `B`) under `b = 17` | `3.871` | `0.500` | `3.291` | `2.253` | `+0.9488` | `+0.9446` |
| Control swap (`A` vs. `B`) under `mean_batch_emb` | `6.069` | `1.492` | `6.828` | `5.065` | `+0.8449` | `+0.8799` |
| Control swap (`A` vs. `B`) under `bypass_zero` (`0 ∈ R^128`) | `6.882` | `2.070` | `11.429` | `5.895` | `+0.7388` | `+0.8569` |
| **Batch swap ONLY (`b = 0` vs. `b = 1`, fixed `X_A`)** (`G4`) | **`21.192`** | `2.935` | `14.396` | `10.879` | **`+0.1161`** | **`+0.1014`** |
| **Batch swap ONLY (`b = 0` vs. `b = 17`, fixed `X_A`)** (`G4`) | **`20.957`** | `2.660` | `15.570` | `9.943` | **`+0.0368`** | **`+0.0492`** |
| **Batch swap ONLY (`b = 0` vs. `b = 35`, fixed `X_A`)** (`G4`) | **`22.182`** | `2.412` | `17.917` | `12.166` | **`+0.0771`** | **`+0.0646`** |
| Joint swap (`(A, b=0)` vs. `(B, b=17)`) | `20.189` | `2.659` | `14.994` | `9.546` | `+0.0551` | `+0.0698` |

### §6.4 Extrapolated Total Real-Cell Inference Wall Time (Apple M3 Pro, `MPS float32`)
1. **Stage 1 — `SE-600M` Control-Cell Embedding (Cached Once Per Control Pool):**
   - `6` local raw-count lines (`nadig_hepg2`, `nadig_jurkat`, `a549_ko`, `dld1`, `cd4t_rest`, `feng_ipsc`) × `512` NTC cells = `3,072` cells.
   - Measured `331 s` / `512` cells -> **`~33.1 minutes`** GPU compute + `~10–15 min` sparse disk I/O to extract NTC cells from raw `.h5ad`/`.h5`/`.mtx.gz` -> **`~45–50 minutes`**. (Adding the 4 optional sensitivity lines `cd4t_stim48`, `a549_ko_cas9pos`, `a549_myllia`, `hek293_xucao` adds `~22 min`.)
2. **Stage 2 — `ST-SE` Primary Inference (`b = 0`, `12.65 ms / pert` on MPS `fp32`):**
   - Group 1 (`2` lines × `1` ckpt × `6` control arms × `~800–931` P5 targets [or `~2,000` full]): `~2.1–5.1 min`.
   - Group 2 (`3` primary lines + `feng_ipsc` × `4` ckpts × `6` control arms = `96` arms × `~86–445` P5 targets [or `~286–1,423` full]): `~4.5–14.5 min`.
   - **Total Stage 2 (`ST-SE` Primary):** **`~7–20 minutes`**.
3. **Stage 3 — `ST-HVG` Secondary + `Perm100` Null (`10.21 ms / pert`):** **`~5–15 minutes`** (base arms) + **`~15–30 minutes`** (`100` gene-label permutations on `ST-HVG`).
4. **Stage 4 — Secondary `G4` 2×2 `batch` vs. `ctrl` Decomposition + Neutral-Batch Robustness (`bypass_zero`, `mean`, `K=8` ensemble):** **`~15–25 minutes`**.
5. **Grand Total Full Gate 5 Real-Cell Execution:** **`~1.5–2.2 hours`** end-to-end on local Apple M3 Pro (`peak RSS ~5.7 GB << 36 GB` RAM; zero cloud GPU required).

---

## §7. Concrete `PREREG_I v5` Edit Recommendations (`PREREG_I.md` NOT Modified Yet)

All edits below are queued for `PREREG_I v5` ONLY after Master Controller signs off on `R16_GATE2_SCHEMA.md`. *(Note: the potential addition of `Occam` as a second open-weight model is held untouched pending human download approval and a separate controller addendum, per `2026-11-25 23:28` directive.)*

### §7.1 Line-by-Line Replacement of Every `TBD_GATE2` Slot in `round16/prereg/PREREG_I.md`

| `PREREG_I.md` Line(s) | v4 Placeholder | Exact `v5` Constant Replacement |
|---|---|---|
| **Line 1, 3, 5** | `v4-DRAFT`, Gate-2 pending disclaimer | Update header to `v5-FINAL-PRE-LOCK`; record `round16/shared/R16_MODEL_CHECKSUMS.sha256` (`77/77 OK`), `R16_GATE2_SCHEMA.md` completion, and zero real cells inferred prior to lock. |
| **Line 62 (`G_model`)** | `TBD_GATE2` | **`|G_model| = 2,000` ordered HVG symbols** (identical across all 8 `ST-SE` & `ST-HVG` checkpoints; `var_dims.pkl["gene_names"]`; first 5 `['HES4','ISG15','MIB2','CEP104','ACOT7']`, last 5 `['MT-CO3','MT-ND3','MT-ND4L','MT-ND5','MT-ND6']`). Harmonized with frozen R10 `HGNCBridge` (`hgnc_complete_set.txt`), yielding fixed `|G_eval^R16(t)| = |canon(G_R9_6148) ∩ canon(G_t^fin) ∩ canon(G_model)|`: `nadig_hepg2 = 1,843`, `nadig_jurkat = 1,843`, `a549_ko = 1,837` (`1,797` raw), `dld1 = 1,832` (`1,792` raw), `cd4t_rest = 1,830` (`1,790` raw), `feng_ipsc = 1,575` (`1,540` raw). |
| **Line 63 (`P_model`)** | `TBD_GATE2` | **`2,024` keys** in `pert_onehot_map.pt` = **`2,023` non-control perturbations** + `'non-targeting'` (identical across all 8 checkpoints; first 5 `['AAAS','AAMP','AAR2','AARS','AARS2']`, last 5 `['ZRANB2','ZRSR2','ZW10','ZWINT','non-targeting']`). Exact `|targets.txt ∩ P_model|` (canon / raw): `nadig_hepg2 = 1,996/1,996`, `nadig_jurkat = 1,999/1,999`, `a549_ko = 515/498`, `dld1 = 490/481`, `cd4t_rest = 1,423/1,387`, `feng_ipsc = 286/282`. **Locked evaluation target rule:** $N_{\mathrm{eval}}^{\mathrm{R16}}(t) = |\text{R9–R12 P5-qualified reproducible targets}(t) \cap \mathrm{canon}(P_{\mathrm{model}})|$ (masking self-gene $g_p$ from $G_{\mathrm{eval}}^{\mathrm{R16}}(t)$ when present; NOT dropping targets whose $g_p \notin G_{\mathrm{model}}$). Locks $\mathcal{G}_2^{\mathrm{qual}} = \{\text{a549\_ko}, \text{dld1}, \text{cd4t\_rest}\}$ ($K_2 = 3 \ge 2$). |
| **Line 64 (`batch` / `cell_type`, `M5.3`)** | `TBD_GATE2` | `cell_type_onehot_map.torch` NEVER enters `StateTransitionPerturbationModel` (`state_transition.py`). `batch_encoder = nn.Embedding(56, 128)` (`batch_col=gem_group`, `'1'..'56'`) is mandatory; **Primary neutral scheme locked:** **`fixed_index_zero` (`batch_idx = 0`)** identically for matched and mismatched arms on Group 1 and Group 2 (matches official `_infer.py:688` unseen-batch fallback; keeps inputs in-distribution with `||e_0||_2 = 10.83`; zero weight edits). Gate-2 synthetic smoke test verified `M5.3` PASS (`SE-600M -> ST-SE` `l2_abs = 5.49`, `max_abs = 0.973`, `l2_delta = 4.25`; `ST-HVG` `l2_abs = 8.92`). Secondary sensitivity (`bypass_zero`, `mean_batch_emb`, `K=8` ensemble) + 2×2 `G4` `batch`-vs-`ctrl` decomposition locked in `§4.4`. |
| **Line 65 (Input/Output & `predict_residual`)** | `TBD_GATE2` | `SE-600M` (`model.safetensors`, loaded in `MPS float32` via `StateEmbeddingModel` + `safetensors.torch.load_file`, `0` missing keys; `bf16`/`fp16` disallowed on MPS due to Metal assertion) takes raw integer UMI counts mapped by exact symbol to `19,790`-gene `protein_embeddings.pt` (`5,120`-dim ESM2), applies internal `log1p(raw_counts)` + top-`2,047` gene ranking, and outputs `2,058`-dim `X_state = [L2-norm 2048 emb || 10 ds_emb]`. `ST-SE` runs `padded=True` `16 × 32` sets with hidden-space residual `project_out(res_pred + control_cells)` (`state_transition.py:509`) and decodes `2,058 -> 2,000` via trained `LatentToGeneDecoder` (`base.py:15–116`). `ST-HVG` takes `libnorm(1e4) -> log1p` on ordered `2,000` `G_model` genes and reads `preds` directly (`gene_decoder` verified 100% untrained init across all 4 folds). Both output `log1p`-normalized expression clipped to `[0.0, 14.0]`; $\hat{\Delta}_{\mathrm{pred}}(p) = \overline{\hat{Y}}(p) - \overline{\hat{Y}}(\text{"non-targeting"})$. |
| **Line 66 (`SE-600M` Gene Coverage)** | `TBD_GATE2` | Exact `genes.txt` -> `protein_embeddings.pt` (`19,790`) lower-bound coverage locked (`nadig_hepg2: 9,022`, `nadig_jurkat: 8,391`, `a549_ko: 13,766`, `dld1: 11,714`, `cd4t_rest: 11,535`, `feng_ipsc: 7,360`; `>> 2,047` top-gene sentence length). Exact full-raw-file mapped gene counts and matched-vs-mismatched `SE-600M` centroid cosine/Euclidean distances logged at Gate 5. |
| **Line 67 (`SE-600M` Pretraining Overlap)** | `TBD_GATE2` | `SE-600M/config.yaml:20,36–42,76` confirms `scbasecamp-cellxgene-tahoe-filtered` (`14,420` datasets); internal CSV manifests are `未能确定` (unshipped). Zero additional pretraining disclosures needed beyond `Q4`: `a549_ko` (`GSE345058`, `2026-08-26`), `dld1` (`GSE337988`, `2026-07-10`), and `cd4t_rest` (`GSE314342`, `2026-02-06`) all strictly post-date `SE-600M` (`2025-06/07`). |
| **Lines 84–92 (§2.3 Table)** | `TBD_GATE2` cells (pretraining & counts) | Fill exact `|G_eval^R16|` (canon/raw), `|targets.txt ∩ P_model|` (canon/raw) + P5 baseline (`N_eval^R16` locked rule), and exact GEO/ENA post-`SE-600M` dates (`GSE337988` `2026-07-10` for `dld1`, `GSE314342` `2026-02-06` for `cd4t_rest`, `GSE345058` `2026-08-26` for `a549_ko`). Also update Line 95 runtime estimate from `~30–60 min` to measured **`~1.5–2.2 hours`** (`SE-600M` MPS `fp32` `~33 min` compute + I/O + `ST-SE`/`ST-HVG`/`G4`). |

### §7.2 Mandatory `PREREG_I v5 §3.1` Conservative `ci975 == 0.0` Tie-Precedence Edit (`2026-11-25 23:20` Ruling)
Replace lines 117–122 of `PREREG_I.md` with the unambiguous 6-line precedence order (explicitly placing `ci975 == 0.0` on its own standalone top-precedence rule before `MATCHED_BETTER` / `MISMATCHED_BETTER` / `EQUIVALENT`):
1. **`INCONCLUSIVE`** if `valid_frac < 0.95` or `97.5%` CI is non-finite (`within_delta = False`);
2. **`INCONCLUSIVE` (Conservative Zero-Boundary Tie Precedence):** if either endpoint of the `97.5%` CI equals `0.0` (`ci975_lo == 0.0` or `ci975_hi == 0.0`), classify as **`INCONCLUSIVE`** regardless of `within_delta` (taking strict precedence over `MATCHED_BETTER`, `MISMATCHED_BETTER`, and `EQUIVALENT`);
3. **`MATCHED_BETTER`** if `ci975_lo > 0.0` *(accompanied by two-sided `within_delta` flag)*;
4. **`MISMATCHED_BETTER`** if `ci975_hi < 0.0` *(accompanied by two-sided `within_delta` flag)*;
5. **`EQUIVALENT`** if NOT `MATCHED_BETTER`, NOT `MISMATCHED_BETTER`, neither `97.5%` CI endpoint equals `0.0`, AND `within_delta == True`;
6. **`INCONCLUSIVE`** for all remaining cases (including exact `±0.04` boundary ties `ci95_lo == -0.04` or `ci95_hi == +0.04`).

- **Main-Line Code & Unit Tests Already Updated & Verified (`2026-11-25`):**
  - [`round16/src/adjudicate_r16.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/src/adjudicate_r16.py) (`SHA-256: 3263aec59995e24e7fa79c5b38702260b2a0d915fda278f27addc16b938896c6`): `if not valid_boot or lo975 == 0.0 or hi975 == 0.0: cat = "INCONCLUSIVE"`.
  - [`round16/src/test_unit_r16.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/src/test_unit_r16.py) (`SHA-256: eb9414004a2e230a442bef32bb2e1bd6b4916debd96b821e7dfc9a3107d93b4c`): added explicit pair-level and line-level tests for `ci975_lo == 0.0`, `ci975_hi == 0.0`, and `[0.0, 0.0]` with `ci95 ⊂ (-0.04, +0.04)` (`within_delta=True` yet `category=="INCONCLUSIVE"` and line returns `MIXED_OR_INCONCLUSIVE`).
  - Verification: `py_compile` rc=0, `pyflakes` rc=0, `8/8` unit tests `OK` under both `/opt/anaconda3/bin/python3 -I` (`1.45s`) and `round16/.venv/bin/python` (`1.31s`).

### §7.3 Addition of `PREREG_I v5 §4.4` (Secondary Descriptive `G4` Decomposition & Neutral-Batch Panel)
Append `§4.4` to `PREREG_I.md §4` verbatim from `§3.4` above (non-confirmatory secondary reporting of `(a) bypass_zero`, `(b) mean_batch_emb`, `(c) K=8` multi-batch ensemble, and the 2×2 `(X_ctrl, batch)` factorial decomposition $\Delta_{\mathrm{batch}}, \Delta_{\mathrm{ctrl}}, \Delta_{\mathrm{both}}, \Delta_{\mathrm{int}}$).
