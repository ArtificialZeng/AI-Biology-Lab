# Round 17 Pre-Registration Protocol Draft v2 (`PREREG_J_DRAFT_v2.md`)

> **STATUS: DRAFT v2 FOR MASTER CONTROLLER AUDIT & LOCK — ADDRESSES CONTROLLER REVIEW R1 ITEMS `J1`–`J7` & `§6` RULINGS 1–3 — ZERO INFERENCE EXECUTED.**
> **Hard-Guard, Blinding, & Release Attestation (`J6`):**
> - **v1 Preserved Untouched:** [`round17/prereg/PREREG_J_DRAFT.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round17/prereg/PREREG_J_DRAFT.md) remains unmodified as historical draft v1.
> - **Zero inference executed:** Zero `STATE` (`ST-SE-Replogle`, `ST-HVG-Replogle`) or `SE-600M` forward passes have been run on any Round 17 candidate line (`0` real cells inferred; `0` MPS/GPU compute).
> - **Model-free local inventory & pre-download reader only:** All local line counts in [`round17/status/LINE_INVENTORY_J.json`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round17/status/LINE_INVENTORY_J.json) come strictly from [`round17/src/build_line_inventory_j.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round17/src/build_line_inventory_j.py) (zero `torch`/`state` imports, zero correlation/metric calculations on target expression deltas). The Pfizer TeloHAEC reader [`round17/src/pfizer_telohaec_reader_r17.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round17/src/pfizer_telohaec_reader_r17.py) is verified pre-download on synthetic fixtures only ([`round17/src/test_pfizer_telohaec_reader_r17.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round17/src/test_pfizer_telohaec_reader_r17.py)); prior to `PREREG_J` lock, only HDF5 structural metadata (keys, shapes, dtypes, `obs`/`var` column names, condition categories) may ever be inspected on `/Users/mac/.cache/r17_telohaec_pfizer/imru_full.h5ad` — zero `X` values, zero perturbation-label lists, and zero counts before lock.
> - **Round 16 `A-T3` Blinding & Release Lock (`J6`):** Because `round17/` cites Round 16 main-line discovery numbers (`§0.2`, `§4.1`), **`round17/` stays uncommitted and strictly hidden from `A-T3` until Master Controller has verified `A-T3`'s Round 16 pass 2.** `A-T3` joins Round 17 **only** after Master Controller formally releases it from Round 16, and then re-implements strictly from the locked `PREREG_J` alone.
> - **Strict governance boundaries:** Zero reads of `round16/audit/` or `A-T3.txt`; zero git operations; zero media material; zero priority claims; always refer to the human strictly as `"Human project lead"`.
> - **Discovery-to-Confirmatory Quarantine:** All 6 Round 16 lines (`nadig_hepg2`, `nadig_jurkat`, `a549_ko`, `dld1`, `cd4t_rest`, `feng_ipsc`) are strictly quarantined as **Discovery Context Only** and are NEVER included in Round 17 confirmatory testing.

- **Project:** AI-Biology-Lab Round 17 — Pre-Registered Confirmatory Audit of Batch-Token (`batch_encoder`) Index Dependence, `STATE - FCT` Verdict Instability, and Label-Free Index Selection in Zero-Shot `Arc STATE` (`ST-SE-Replogle` Primary; `ST-HVG-Replogle` Co-Evaluated Second Model)
- **Protocol Author:** Main Line (`Antigravity`), revised strictly to close Master Controller Review R1 (`controller/msgs/r17_REVIEW_R1_MAIN.txt`, Items `J1`–`J7` & `§6` Rulings 1–3)
- **Companion Deliverables:**
  - Prior-Art Review: [`round17/prereg/PRIOR_ART_J.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round17/prereg/PRIOR_ART_J.md)
  - Model-Free Local Line Inventory: [`round17/status/LINE_INVENTORY_J.json`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round17/status/LINE_INVENTORY_J.json) & [`round17/src/build_line_inventory_j.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round17/src/build_line_inventory_j.py)
  - Pre-Download Pfizer TeloHAEC Reader & Synthetic Test: [`round17/src/pfizer_telohaec_reader_r17.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round17/src/pfizer_telohaec_reader_r17.py) & [`round17/src/test_pfizer_telohaec_reader_r17.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round17/src/test_pfizer_telohaec_reader_r17.py)
  - Parallel Candidate Screening (Swarm, Model-Free, `J7`): `round17/status/R17_CANDIDATES.md` (`A-T1`, may follow separately without blocking `PREREG_J` v2)
- **Parent Protocols (Locked, Unmodified):** [`round16/prereg/PREREG_I.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md) (`v5`, `c1fc7c6`), [`round16/prereg/PREREG_I_ADDENDUM_A.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I_ADDENDUM_A.md) (`Items 1–11`).
- **Fresh Unshared Round 17 RNG Seeds:**
  - Primary control sampling, `SE-600M` per-pool encoding, perturbation half-splits (`s = 0..19` via `20261701 + s`), joint perturbation bootstrapping (`n_boot = 4,000`), and `J5` index-label permutation null (`n_perm = 10,000`): **`20261701`** (`np.random.default_rng(20261701)`).
  - Secondary within-line disjoint control resample (when applicable): **`20261702`**.

---

## §0. Scientific Rationale, Architectural Mechanism, & Round 16 Discovery Context

### §0.1 Architectural Mechanism Verified Directly From Shipped Weights & `arc-state 0.11.1` Source
Read-only inspection of `arc-state 0.11.1` source code and all 8 shipped zero-shot checkpoints (`round16/shared/models/ST-{SE,HVG}-Replogle/zeroshot/{hepg2,jurkat,k562,rpe1}`, verified `77/77` SHA-256 OK, zero inference) establishes four exact facts:

1. **Additive `batch_encoder` Dominates Input Energy:** Both `ST-SE-Replogle` and `ST-HVG-Replogle` instantiate `self.batch_encoder = nn.Embedding(56, 128)` (`batch_dim = 56`, `hidden_dim = 128`, `batch_col = gem_group`, `use_batch_token = false`, `batch_predictor = false`; `state/tx/models/state_transition.py:192–197`) and add the retrieved vector unnormalized into the sequence input (`seq_input = pert_embedding + control_cells + self.batch_encoder(batch_indices.long())`, `state_transition.py:422–441`). Across all 8 checkpoints, row L2 norms are `9.77–12.43` (median `11.09`) in `ST-SE` and `8.81–12.52` (median `10.88`) in `ST-HVG` (near PyTorch $\mathcal{N}(0, 1)^{56 \times 128}$ initialization norm $\sqrt{128} = 11.31$), compared to `||basal_encoder(X_ctrl)||_2 ≈ 0.38–1.20` and `||pert_encoder(pert)||_2 ≈ 0.74`.
2. **Nearly Orthogonal Rows Shared Across All 4 Zero-Shot Folds:** Within each fold, the 56 rows are nearly mutually orthogonal (off-diagonal row cosine `mean = 0.0014, sd = 0.0919` in `ST-SE`; `0.0006 ± 0.0896` in `ST-HVG`), so averaging rows in embedding space (`mean_batch_emb`) collapses vector norm by $\approx \sqrt{56} \approx 7.5\times$ (`||e_mean||_2 = 1.55` in `ST-SE`, `1.47` in `ST-HVG`). Across the 4 zero-shot folds (`hepg2`, `jurkat`, `k562`, `rpe1`), identical row indices remain tightly aligned (`mean_cos = 0.9995–0.9997`, `min_cos >= 0.9991` in both `ST-SE` and `ST-HVG`), confirming that applying a single index $b \in \{0..55\}$ consistently across the 4-fold ensemble (`zeroshot/ensemble4`) is geometrically coherent.
3. **Silent Unseen-Batch Fallback to `batch_idx = 0` (`J4`, Verified at Gate 2):** At zero-shot inference on unseen cell lines, official `arc-state 0.11.1` (`state/_cli/_tx/_infer.py:678–697`, specifically **`state/_cli/_tx/_infer.py:688`**, verified in `round16/status/R16_GATE2_SCHEMA.md`) looks up `.obs[batch_col]` in `batch_onehot_map.torch` (`{'1'->0, '10'->1, ..., '56'->51, '6'->52, ..., '9'->55}`) and silently sets **`idxs[i] = 0  # fallback to zero`** (`fixed_index_zero`) for every unseen batch label (and crashes with `KeyError: 'batch'` at `state_transition.py:427` if no batch column is present).
4. **Unshipped Batch Provenance (`未能确定`):** Exhaustive read-only inspection of all 77 files under `round16/shared/models/` confirms that **neither per-`gem_group` training-cell embeddings nor a `gem_group`-to-training-cell-line lookup table are shipped**.

### §0.2 Round 16 Discovery Observations (Strictly Discovery Context Only — Excluded from Round 17 Confirmatory Family)
In Round 16 (`round16/results/summary_r16.json` and `round16/results/sensitivity/summary_r16.json`), secondary neutral-batch and 4-index (`{0, 1, 17, 35}`) diagnostics on the 6 Round 16 lines surfaced pronounced, line-varying batch-token dependence that motivates Round 17's confirmatory study:

- On **`a549_ko`** (`N = 403`, `FCT r_orth_cal = 0.2134`): `STATE_SE` shifted across just 4 trained indices from `b=17: 0.0863` and `b=0 (fixed_index_zero): 0.1143` (`STATE_SE - FCT = -0.0991`, `FCT_SUPERIOR_TO_MODEL`) up to `b=35: 0.2426` (`4`-index range **`0.1563`**, `SD = 0.0702`), while `multi_batch_ens8` reached `0.2874` (`STATE_SE - FCT = +0.0740`, `97.5% CI [+0.0412, +0.1090]`, flipping the `H-I2` verdict to `MODEL_SUPERIOR_TO_FCT`).
- On **`nadig_hepg2`** (`N = 740`, `FCT = 0.2015`): `fixed_index_zero` gave `0.1053` (`FCT_SUPERIOR_TO_MODEL`), whereas `bypass_zero` reached `0.2341` (`STATE_SE - FCT = +0.0326 [+0.0139, +0.0508]`, flipping to `MODEL_SUPERIOR_TO_FCT`) and `multi_batch_ens8` / `mean_batch_emb` reached `0.2008` / `0.2022` (`INDISTINGUISHABLE`).
- On **`dld1`** (`N = 75`), **`feng_ipsc`** (`N = 88`), and **`cd4t_rest`** (`N = 245`): the 4-index range spanned **`0.0598`**, **`0.0539`**, and **`0.0374`**, respectively (with `cd4t_rest` exhibiting very weak overall `STATE_SE` signal), and **different indices ranked highest on different cell lines** (`b=35` on `a549_ko`; `b=0` on `dld1` & `cd4t_rest`; `b=1` on `feng_ipsc`), while `G4` 2×2 `b0 vs b17` median `batch_share` across unseen mismatched controls reached `57.7%–94.4%`.

Because these patterns were discovered on the 6 Round 16 lines, **none of `{nadig_hepg2, nadig_jurkat, a549_ko, dld1, cd4t_rest, feng_ipsc}` may enter Round 17 confirmatory testing**. Round 17 pre-registers a full 56-index confirmatory audit on fresh, never-STATE-run lines.

### §0.3 Prior-Art Positioning (Zero Priority Claims)
Full primary-source verification is logged in [`round17/prereg/PRIOR_ART_J.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round17/prereg/PRIOR_ART_J.md). Round 17 explicitly credits **`CELLAUDIT`** (`arXiv:2609.27234`, frozen-checkpoint input replacement & factorial auditing), **`Zhang et al. (2026)`** (`arXiv:2611.05849` §C.5, `(control, batch / gem_group)` ridge probes reaching PDC `0.27–0.46` and internal `L0` direction swapping), **`Molina & Zhang (2026)`** (`10.64898/2026.07.24.740459`, four-component response decomposition on Replogle/Nadig), **`OCOO-T`** (`arXiv:2606.12838`), and **`VCBench`** (`arXiv:2604.27646`). **No priority claims are made.**

---

## §1. Audited Models, Selection-Free Index Set (`J4`), Fixed Inference Specification, & Endpoint Parity

### §1.1 Explicit Inclusion of Both `ST-SE` (Primary) and `ST-HVG` (Co-Evaluated Second Model)
Both zero-shot `Arc STATE` families shipped under `round16/shared/models/` carry the exact same `batch_encoder = nn.Embedding(56, 128)` (`batch_dim = 56`) architecture and cross-fold geometry (`mean_cos = 0.9997`). Therefore, **`ST-HVG-Replogle` IS EXPLICITLY INCLUDED** as the co-evaluated second model under the **exact same `H-J1`, `H-J2`, `H-J3`, and secondary `J5` (`G4` 2×2 + cross-line rank-correlation) rules**:

| Role | Model Key | Checkpoints (All 4/4 Unseen on All Fresh Lines) | Control Input & Decoding Lock (Identical to `PREREG_I §1.2 & §2.1`) | Batch-Encoder Architecture |
|---|---|---|---|---|
| **Primary Model** | **`STATE_SE`** | `ST-SE-Replogle/zeroshot/{hepg2,jurkat,k562,rpe1}` equal-weight 4-fold mean (`zeroshot/ensemble4`) | Raw integer UMI NTC counts $\to$ `SE-600M` (`MPS float32`, official `FilteredGenesCounts` + `VCIDatasetSentenceCollator(training=False)`, per-pool seed `20261701`, `bs=16`, `shuffle=False`, `num_workers=0`) $\to$ `X_state` `(n_ctrl, 2058)` $\to$ `ST-SE` (`padded=True`, `16×32` sets) $\to$ `clip(gene_decoder(preds), 0, 14)` | `nn.Embedding(56, 128)`, row norms `9.77–12.43` (med `11.09`), cross-fold `cos >= 0.9991` |
| **Co-Evaluated Second Model** | **`STATE_HVG`** | `ST-HVG-Replogle/zeroshot/{hepg2,jurkat,k562,rpe1}` equal-weight 4-fold mean (`zeroshot/ensemble4`) | Raw integer UMI NTC counts $\to$ official `arc-state` `official_median_log1p` (`sc.pp.normalize_total(target_sum=None) -> log1p`, `HGNCBridge` onto ordered `2,000` `G_model` columns) $\to$ `ST-HVG` (`padded=True`, `16×32` sets) $\to$ `clip(preds, 0, 14)` directly (**NEVER** untrained `gene_decoder`) | `nn.Embedding(56, 128)`, row norms `8.81–12.52` (med `10.88`), cross-fold `cos >= 0.9993` |
| **Reference Baseline** | **`FCT`** | Frozen `round10/src/fct_pipeline.py` on `r9_6` (`6` source lines; **`r9_6 \ {vcc2025_h1}`** primary for `gse295214_h1` per `PREREG_I §2.4`, full `r9_6` secondary) | Exact `PREREG_I §1.2 & §2.4` `FCT` evaluated on identical `(T_eval(t), G_eval(t))` | N/A (closed-form target baseline mean scaling) |

*Scientific rationale for including `STATE_HVG`:* Because `STATE_HVG` bypasses `SE-600M` and `LatentToGeneDecoder` completely while sharing `StateTransitionPerturbationModel`'s Llama backbone and `batch_encoder = nn.Embedding(56, 128)`, running identical `H-J1–H-J3` and `J5` tests on `STATE_SE` and `STATE_HVG` isolates whether batch-index instability is intrinsic to the transition backbone + `batch_encoder` rather than an artifact of `SE-600M` embeddings or the latent-to-gene decoder.

### §1.2 Selection-Free 56-Index Sweep Lock (`J4`), Control-Cell Sampling, & Per-Index Prediction Formula
- **Selection-Free Index Set Lock (`J4`):**
  - Round 17 evaluates **all 56 trained batch-embedding indices $b \in \{0, 1, \dots, 55\}$** (`batch_dim = 56`) exhaustively on every included line and model.
  - Row **`0`** (`fixed_index_zero`) is explicitly anchored as the official **`arc-state 0.11.1` unseen-batch fallback** (`state/_cli/_tx/_infer.py:688`, verified at Round 16 Gate 2 in `round16/status/R16_GATE2_SCHEMA.md`).
  - **No-Drop Lock (`J4`):** **No index in `{0, 1, ..., 55}` is EVER dropped, filtered, winsorized, or excluded after inference** under any circumstance (all 56 indices always enter `H-J1`, `H-J2`, `H-J3`, and `J5` verbatim).
- **Control Sampling (`PREREG_I §2.2` with Fresh Seed `20261701`):** For each line $t$, sample `n_ctrl = 512` QC-passing NTC cells without replacement proportionally stratified across experimental batches using exact integer Hamilton largest-remainder quotas (`seed = 20261701`) and partition into `16 × 32` sets (`padded=True`). All 56 batch indices and all perturbations (including `'non-targeting'`) share the **exact same sampled control cells and `16 × 32` ordering** (`control_barcodes_r17.json`).
- **Per-Index Prediction Matrix ($\hat{\Delta}^{(m, b)}_t$):** For each model $m \in \{\texttt{STATE\_SE}, \texttt{STATE\_HVG}\}$ and each trained batch index $b \in \{0, 1, \dots, 55\}$ applied identically to all `n_ctrl` cells (`batch_indices = torch.full((n_ctrl,), b, dtype=torch.long)`), compute per-fold decoded pseudobulk expression in clipped `[0, 14]` space, average across all 4 zero-shot folds (`zeroshot/{hepg2,jurkat,k562,rpe1}`), and subtract the 4-fold average predicted `'non-targeting'` expression under that **same** batch index $b$:
  $$\overline{\hat{Y}}^{(m, b)}_t(p) = \frac{1}{4}\sum_{f=1}^{4} \frac{1}{n_{\mathrm{ctrl}}}\sum_{i=1}^{n_{\mathrm{ctrl}}} \operatorname{clip}_{[0, 14]}\!\big(\hat{Y}^{(m, f)}(x_i, p, b)\big), \qquad \hat{\Delta}^{(m, b)}_t(p) = \overline{\hat{Y}}^{(m, b)}_t(p) - \overline{\hat{Y}}^{(m, b)}_t(\text{"non-targeting"}).$$
  *(Per-fold 56-index summaries are also saved as secondary diagnostics to verify that batch-index effects are not driven by a single fold.)*
- **Primary Endpoint (`r_orth_cal`) & Secondary Endpoint (`r_dorth_cal`):** Identical to `PREREG_I §1.1 & §4.1` (`round16/src/metrics_r16.py`): unselected `all_nontarget` ($G_{\mathrm{eval}}(t) \setminus \{g_p\}$), orthogonalized to per-perturbation `Transfer_Mean` $\bar{\Delta}^{\mathrm{src}}_p$ (fitted on `r9_6`, or `r9_6 \ {vcc2025_h1}` primary for `gse295214_h1`), Spearman–Brown disattenuated by orthogonal split-half correlation explicitly keyed **`rho_ab_orth`**, clipped to `[-1, 1]`.

---

## §2. Confirmatory Cohort (`J1`, Track 1 + Track 2 -> $K_{\mathrm{conf}} \in \{2, 3\}$), Pfizer TeloHAEC Reader Spec, & Complete Line Inventory

### §2.1 Locked Confirmatory Cohort Architecture (`J1`, Ruling 1) & Sequencing Rule
Per Master Controller Review R1 (`J1`, Ruling 1), Round 17 combines **Track 1 + Track 2** into a multi-context confirmatory cohort reaching **$K_{\mathrm{conf}} \in \{2, 3\}$ independent confirmatory units** (avoiding reliance on a single confirmatory line sharing a study/lineage with discovery):

| Confirmatory Unit # | Line Key | Biological Context / Modality | Study / Release (`SE-600M` `2025-06` Chronology) | Approval & Qualifying Gate (`PREREG_I` Floors: $N_{\mathrm{eval}} \ge 50$, $\|G_{\mathrm{eval}}^{\mathrm{R9}}\| \ge 3{,}000$, per-half $\mathrm{NTC} \ge 15$) | Mandatory Disclosures & Independence Status |
|---|---|---|---|---|---|
| **Unit 1 (Local, Track 1)** | **`cd4t_stim48`** | Primary human CD4$^+$ T cells (48h stimulated), `CRISPRi` (`10x Flex`) | `GSE314342`, `2026-02-06` (**Post-`SE-600M`**) | **QUALIFIED** (`N_eval = 401`, `\|G_eval^R9\| = 6,131`, `\|G_eval^STATE\| = 1,836`, `NTC = 166,871`, `78,168 / 88,703`) | Disjoint donors (`D2/D3` vs. `D1/D4`), disjoint lanes (`R2L25–48`), and stimulated vs. resting state relative to R16 discovery `cd4t_rest`, **but shares primary CD4$^+$ T lineage, sgRNA library, and Marson lab (`GSE314342`)** — mandatory shared-study disclosure. |
| **Unit 2 (Prospective, Track 2A)** | **`gse295214_h1`** | Human `H1` hESC, `CRISPRi` (`10x 3'`, 85 `.h5` lanes) | `GSE295214`, `2026-05-05` (**Post-`SE-600M`**) | **Human-approved** (`round16/status/R16_GSE295214_HUMAN_APPROVAL.txt`); **CONFIRMATORY ONLY IF** frozen R16 Task 1 counts-only QC (`prepare_gse295214_line_r16`, `PREREG_I §2.4`) passes unchanged (`PRIMARY` status). | Strictly study- and lineage-independent of all 6 R16 discovery lines; disclose `H1` lineage & Replogle/Weissman lab provenance (`PREREG_I §2.4`); `FCT` primary uses `r9_6 \ {vcc2025_h1}` (full `r9_6` sensitivity). |
| **Unit 3 (Prospective, Track 2B)** | **`pfizer_telohaec_basal`** *(Pfizer TeloHAEC `Untreated`)* | Immortalized human aortic endothelial (`TeloHAEC`, basal `Untreated`), `CRISPRi` (`Perturb-seq`) | **Zenodo `18792213`**, `imru_full.h5ad` (`53,290,123,849 B`, `md5 0594ddd25e224e5e9cc22414e71daf48`, **`CC BY 4.0`**, released **`2026-03-04`**); Wong et al., *bioRxiv* `10.64898/2026.03.04.709607` (**Post-`SE-600M`**) | **Human download approved** (`controller/approvals/R17_TELOHAEC_HUMAN_APPROVAL.txt`; structure-only pre-lock); **CONFIRMATORY ONLY AFTER** post-lock counts-only QC (`prepare_pfizer_telohaec_r17`) passes all `PREREG_I` floors unchanged. | Strictly study- and lineage-independent of all 6 R16 discovery lines and of `cd4t_stim48`/`gse295214_h1`. **Mandatory Disclosures:** (i) `Untreated` (basal) is the sole primary confirmatory condition for this unit; `Treated` (`IL-1β + TNF-α`, `pfizer_telohaec_treated`) is the same cells under a second treatment condition and is reported strictly as a **secondary within-family replication**, NEVER counted as a second independent confirmatory line; (ii) the TeloHAEC lineage also underlies the local micro-N `telohaec` line (`GSE210681` / `IGVFDS7943IIWZ`; `N = 9`; never run through `STATE`). |

#### Mandatory Prospective Sequencing Rule for `gse295214_h1` (`J1`)
- **`PREREG_J` IS LOCKED BEFORE ANY `STATE` INFERENCE ON `gse295214_h1`.**
- Master Controller holds the Round 16 Group 3 execution GO (`R16_G3_GO`) until `PREREG_J` is formally locked (`PREREG_J_COMMIT.txt`).
- Consequently, both Round 16's Group 3 index-0 / `{0, 1, 17, 35}` evaluation and Round 17's full 56-index confirmatory evaluation on `gse295214_h1` are **strictly prospective** (zero prior `STATE` inference on `gse295214_h1` before lock).
- Likewise, for `pfizer_telohaec_basal` (`imru_full.h5ad`), **zero `X` values, zero perturbation-label lists, and zero cell counts are inspected before `PREREG_J` lock** (HDF5 structural keys/shapes/dtypes/column names/condition labels only if needed pre-lock; counts-only QC executes strictly after `PREREG_J` lock, and inference executes strictly after code freeze and `R17_GO`).

#### Fixed Cohort-Level Claim Rule (`J1`)
- Report **`k / K_conf` per model (`STATE_SE` and `STATE_HVG` separately)** with **ZERO pooling of perturbations across lines**.
- A **general cross-context claim** (*"STATE zero-shot accuracy depends on the batch index on unseen contexts"*) is permitted **IF AND ONLY IF** a model achieves **`BATCH_DEPENDENCE_PRACTICALLY_MEANINGFUL` (`H-J1`) on at least `2` independent confirmatory units** (`k >= 2` out of `K_conf` independent units; noting explicitly that `pfizer_telohaec_treated` is NOT an independent unit and NEVER increments `k` or `K_conf`).
- If `k < 2` (or if `K_conf == 1`), **ONLY per-line statements** are made — no general cross-context claim is permitted.

---

### §2.2 Complete Local Never-STATE-Run Inventory & Resolved Scope Rulings (`§6` Ruling 2)
All 9 local never-STATE-run lines audited model-free in [`round17/status/LINE_INVENTORY_J.json`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round17/status/LINE_INVENTORY_J.json), updated per Master Controller Ruling 2:

| Round | Line Key | Modality / Assay | Study / Release (`SE-600M` Chronology) | Raw Int UMI Verified | Prior `STATE` Exposure | `NTC` Total (`Half A / B`; `samp`) | `\|G_eval^R9\|` / **`\|G_eval^STATE\|`** | **`N_eval` (`PREREG_I`)** (`half_valid`) | Locked R17 Role (Rulings 1–2) |
|---|---|---|---|---|---|---|---:|---:|---|
| `r11` | **`cd4t_stim48`** | `CRISPRi` (`10x Flex`) | `GSE314342`, `2026-02-06` (Post-`SE-600M`) | Yes (`27.39 GB`) | **`NONE`** | `166,871` (`78,168 / 88,703`; `512`) | `6,131` / **`1,836`** | **`401`** (`333`) | **CONFIRMATORY UNIT 1** *(shared-study disclosure vs. R16 `cd4t_rest`)* |
| `r14` | `a549_myllia` | `CRISPRko` (`clonal Cas9`) | `GSE337804`, `2026` (Post-`SE-600M`) | Yes (`1.17 GB`) | **`NONE`** | `2,093` (`1,033 / 1,060`; `512`) | `6,126` / **`1,836`** | **`20`** (`20`) | **SECONDARY DESCRIPTIVE** *(underpowered `N = 20`)* |
| `r12` | `hek293_xucao` | `CRISPRi` (`PerturbSci-Kinetics`) | `GSE218566`, `2023` (Predates `SE-600M`) | Yes (`0.34 GB`) | **`NONE`** | `2,757` (`1,378 / 1,379`; `512`) | `6,126` / **`1,838`** | **`16`** (`16`) | **SECONDARY DESCRIPTIVE** *(underpowered `N = 16`)* |
| `r11` | `telohaec` | `CRISPRi` (`10x 3'`) | `IGVFDS7943IIWZ` / `GSE210681`, `2023–24` | Yes (`1.47 GB`) | **`NONE`** | `2,215` (`1,142 / 1,073`; `512`) | `4,869` / **`1,475`** | **`9`** (`9`) | **DESCRIPTIVE MICRO-N** *(run; shares TeloHAEC lineage with Pfizer Unit 3)* |
| `r10` | `calu3` | `CRISPRi` (`10x 3'`) | `GSE208240`, `2023` | Yes (`0.74 GB`) | **`NONE`** | `3,094` (`1,486 / 1,608`; `512`) | `6,134` / **`1,840`** | **`8`** (`8`) | **DESCRIPTIVE MICRO-N** *(run)* |
| `r13` | `mac_lps` | `CRISPRko` (`10x Flex`) | `GSE327644`, `2026-08` (Post-`SE-600M`) | Yes (`3.57 GB`) | **`NONE`** | `27,205` (`14,605 / 12,600`; `512`) | `6,132` / **`1,837`** | **`4`** (`4`) | **DESCRIPTIVE MICRO-N** *(run; paired with `mac_ctrl`)* |
| `r13` | `mac_ctrl` | `CRISPRko` (`10x Flex`) | `GSE327644`, `2026-08` (Post-`SE-600M`) | Yes (`3.01 GB`) | **`NONE`** | `24,679` (`13,354 / 11,325`; `512`) | `6,132` / **`1,838`** | **`3`** (`3`) | **DESCRIPTIVE MICRO-N** *(run; paired with `mac_lps`)* |
| `r12` | `a549_ko_cas9pos` | `CRISPRko` (`CAS9>0` subset) | `GSE345058`, `2026-08-26` (Post-`SE-600M`) | Yes (`10.96 GB`) | **`NONE`** | `33,685` (`17,435 / 16,250`; `512`) | `6,121` / **`1,837`** | **`403`** (`403`) | **R16 SUBSET SENSITIVITY ONLY** *(run; strictly excluded from confirmatory & `K_conf`; 73.3% subset of R16 `a549_ko`)* |
| `r10` | `wtc11_i3neuron` | `CRISPRi` (`10x 3'`) | `GSE152988`, `2021` | Yes (`0.29 GB`) | **`NONE`** | `423` (`194 / 229`) | `5,995` / **`1,746`** | **`1`** (`1`) | **INFERENCE OMITTED** *(Ruling 2: `N = 1`, degenerate bootstrap & half-split; zero inference run)* |

*(Additional parallel candidates screened model-free by `A-T1` under `J7` — `round17/status/R17_CANDIDATES.md`, public post-`2025-06` CRISPR Perturb-seq screens absent from STATE-Replogle training with plausible $N_{\mathrm{eval}} \ge 50$, small HTTP range checks `< 100 MB` total, zero full downloads, no VCC 2026 participant-channel data — may be added to Track 2 ONLY after separate Human project lead download approval and post-lock counts-only QC under identical `PREREG_I` floors and `PREREG_J` rules.)*

---

### §2.3 Pre-Download Locked Reader Specification & Synthetic Verification for Pfizer TeloHAEC (`Zenodo 18792213`, `imru_full.h5ad`)
Implemented and verified pre-download in [`round17/src/pfizer_telohaec_reader_r17.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round17/src/pfizer_telohaec_reader_r17.py) with synthetic fixture suite [`round17/src/test_pfizer_telohaec_reader_r17.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round17/src/test_pfizer_telohaec_reader_r17.py):

1. **Source Artifact & Checksum Pin:**
   - File: `/Users/mac/.cache/r17_telohaec_pfizer/imru_full.h5ad` (local internal SSD cache, never sync share); Zenodo record `18792213`, exact size `53,290,123,849 B`, `md5 0594ddd25e224e5e9cc22414e71daf48`, licence `CC BY 4.0`, released `2026-03-04` (Wong et al., *bioRxiv* `10.64898/2026.03.04.709607`; schema documented in `pfizer-opensource/phenotype_reversion`).
2. **Bounded-RAM HDF5 Ingestion & Raw Integer Count Slot Enforcement:**
   - Never loads full `53.3 GB` `.h5ad` into memory; streams `obs`/`var` metadata first (decoding both modern AnnData categoricals `codes/categories` and legacy `__categories/<col>` tables — never using raw integer codes as string labels), then streams strictly selected rows of the verified raw integer UMI slot in fixed row chunks (`chunk_size = 4,096`) across `CSR`, `CSC`, or `Dense` HDF5 encodings.
   - Resolves the exact integer UMI matrix via frozen Round 14 `_resolve_h5ad_integer_slot` (`layers/counts` $\to$ `layers/raw_counts` $\to$ `raw/X` $\to$ `X`, verifying `frac_integer >= 0.9999`, `max >= 5`, `min >= 0`); **fails closed** (`RuntimeError("SLOT_AMBIGUOUS")`) if no verified integer count slot exists.
3. **Strict Condition Stratification (Zero Mixing Across Conditions):**
   - Reads `obs["condition"]` (fallback `cell_treat` if `condition` absent) and filters strictly to exact normalized label:
     - Primary confirmatory unit (`pfizer_telohaec_basal`): `condition == "Untreated"` (case-insensitive exact match to `untreated` / `basal` / `control` / `unstimulated`).
     - Secondary within-family replication (`pfizer_telohaec_treated`): `condition == "Treated"` (`IL-1β + TNF-α`, case-insensitive exact match to `treated` / `stimulated` / `il1b_tnfa` / `cytokine`).
   - Cells with unknown/ambiguous condition labels are excluded and counted in the reader audit JSON.
4. **Strict Perturbation & NTC Classification (Zero `SAFE_TARGET` Contamination):**
   - Reads `obs["gene_target"]` (with `obs["guide"]` auxiliary check; optional `n_guides == 1` filter if present):
     - **NTC (`ntc`):** Strict match to `NO-TARGET` / `NO_TARGET` / `NON-TARGETING` / `NTC` (case-insensitive).
     - **Excluded (`excluded_safe_target`):** `SAFE_TARGET` / `SAFE-TARGET` / `SAFE_TARGETING` / `INTERGENIC` / `CONTROL` (non-NTC) / multi-gene guides (`+`, `,`, `;`, `|`, `/`) / unassigned / `NA` / empty — **strictly excluded from BOTH NTC controls and target perturbations**.
     - **Single-gene targeting (`target`):** Single clean gene symbol canonicalized via frozen Round 10 `HGNCBridge`.
5. **Gene Filtering, Cell QC, & Deterministic Split-Half Assignment:**
   - Retains strictly `Gene Expression` / `RNA` features from `var` (excluding `CRISPR Guide Capture` / non-human transgenes), stripping Ensembl `.N` versions and applying frozen Round 10 `select_canonical_genes` (leading-occurrence deduplication) for gene names.
   - Cell QC computed strictly on RNA features: `n_umis >= 300`, `n_genes >= 200`, `max_pct_mt = None` (matching `PREREG_I §2.4` / `R10–R16`).
   - Split-half assignment (`PREREG_H §6` deterministic hierarchy):
     1. **Tier 1 (`flask`):** Lexicographical 1-based odd/even parity across unique `obs["flask"]` strata when `>= 2` flasks exist and both halves contain `>= 15` NTCs;
     2. **Tier 2 (`welltag`):** Lexicographical 1-based odd/even parity across unique `obs["welltag"]` strata if `flask` fails the `>= 15` NTC floor;
     3. **Tier 3 (`seed=20261701`):** Deterministic 50/50 within-stratum random split (`np.random.default_rng(20261701)`) only if neither biological replicate stratum yields `>= 15` per-half NTCs.
6. **Post-Lock Counts-Only Qualifying Gate (`prepare_pfizer_telohaec_r17`):**
   - Runs frozen Round 10 `process_r10_line` (`n_cells >= 25`, `ad_counts >= 50`, `energy_test_p_value < 0.001`) strictly on counts/metadata (zero model predictions, zero correlations) and checks exact `PREREG_I` floors (`N_eval >= 50` vs. `r9_6 ∩ canon(P_model)`, `|G_eval^R9| >= 3,000`, `per-half NTC >= 15`), writing `round17/shared/R17_PFIZER_TELOHAEC_QC_COUNTS.md`. If `pfizer_telohaec_basal` fails any floor, it is automatically downgraded to `EXCLUDED_FEASIBILITY` (or `SENSITIVITY_ONLY` if `20 <= N_eval < 50`) with zero protocol changes.
7. **Pre-Download Synthetic Test Verification:**
   - [`round17/src/test_pfizer_telohaec_reader_r17.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round17/src/test_pfizer_telohaec_reader_r17.py) verifies all 4 core invariants on synthetic `.h5ad` fixtures (`py_compile = 0`, `pyflakes = 0`, `4/4` tests PASS): (i) bit-exact `CSR`/`CSC`/`Dense` + legacy `__categories` parity; (ii) strict `Untreated` vs `Treated` separation & `SAFE_TARGET` / multi-gene exclusion; (iii) NTC raw-count extractor + end-to-end `process_r10_line` execution; (iv) fail-closed `EXCLUDED_FEASIBILITY` and `SLOT_AMBIGUOUS` guards.

---

## §3. Pre-Registered Hypotheses, Winner's-Curse & Floor-Effect Controls (`J2`), Mechanism Tie-In (`J5`), & Multiplicity Control

Every confirmatory unit $t \in \mathcal{C}_{\mathrm{conf}}$ (`cd4t_stim48`, plus `gse295214_h1` if R16 counts-only QC passes, plus `pfizer_telohaec_basal` if post-lock counts-only QC passes) is evaluated on both primary **`STATE_SE`** and co-evaluated second model **`STATE_HVG`** under the identical rules below. `pfizer_telohaec_treated` (within-family replication), secondary fresh lines (`a549_myllia`, `hek293_xucao`), micro-N descriptive lines (`telohaec`, `calu3`, `mac_ctrl`, `mac_lps`), and R16 subset sensitivity (`a549_ko_cas9pos`) are scored identically but reported strictly outside the confirmatory family (and `wtc11_i3neuron` is omitted per Ruling 2).

### §3.1 Multiplicity Control & Joint Bootstrap Architecture
- **Shared Joint Perturbation Bootstrap:** On each line $t$, generate `n_boot = 4,000` perturbation-resampling index draws (`seed = 20261701`, drawing $N_t$ perturbations with replacement per draw) shared identically across all 56 batch indices $b \in \{0..55\}$, `FCT`, `L1_ens56`, `L2_ntc_recon_match`, and secondary comparators, producing the exact `4000 × 56` bootstrap matrix $R^{\mathrm{boot}}_{t,m}$. Validity requires `valid_frac >= 0.95` and finite `rho_ab_orth > 0` on the observed sample.
- **Family-Wise Error Control ($\mathrm{FWER} \le 0.05$):**
  - Within each confirmatory unit $t$ and model $m$, there are **3 co-primary inferential hypothesis tests**:
    1. **`H-J1`**: Batch-index landscape heterogeneity & winner's-curse-free cross-fitted spread ($\widehat{\Delta}_{\mathrm{CF}, 0}$ and variance-component test $p_{\mathrm{var}}$);
    2. **`H-J3a`**: Label-free uniform 56-index ensemble `L1_ens56` vs. `fixed_index_zero` (`b = 0`);
    3. **`H-J3b`**: Label-free NTC control-reconstruction index `L2_ntc_recon_match` vs. `fixed_index_zero` (`b = 0`).
  - *(Note: `H-J2` is an audit quantification of verdict-call instability using the exact unchanged `PREREG_I §3.2` `97.5%` CI rule at each index $b$, paired with a Bonferroni-56 hardened flip check.)*
  - **Primary Family (`STATE_SE`, size $M_1 = 3 K_{\mathrm{conf}}$, i.e. `6` tests for $K_{\mathrm{conf}} = 2$ or `9` tests for $K_{\mathrm{conf}} = 3$):**
    - Within any single line, the base 3-test Bonferroni level is $\alpha_{\mathrm{line}} = 0.05 / 3 = 0.016667$ (two-sided **`98.333%` percentile bootstrap CI**, quantiles `[0.8333%, 99.1667%]`, $z = 2.3940$).
    - Across all $M_1 = 3 K_{\mathrm{conf}}$ confirmatory tests in the `STATE_SE` family, apply **Holm–Bonferroni step-down FWER control** ($\alpha_{(j)} = 0.05 / (3 K_{\mathrm{conf}} - j + 1)$ for ordered $p_{(1)} \le \dots \le p_{(M_1)}$, verifying equivalence with the corresponding two-sided `(1 - \alpha_{(j)})` percentile bootstrap CI; worst-case step-1 level is `0.008333` [`99.167%` CI, $z = 2.6383$] at $K_{\mathrm{conf}} = 2$ and `0.005556` [`99.444%` CI, $z = 2.7729$] at $K_{\mathrm{conf}} = 3$).
  - **Co-Evaluated Second Model Family (`STATE_HVG`, size $M_2 = 3 K_{\mathrm{conf}}$):** Identical independent Holm–Bonferroni FWER $\le 0.05$ control, PLUS a mandatory **pooled-model sensitivity** controlling FWER across all $6 K_{\mathrm{conf}}$ confirmatory tests (`STATE_SE + STATE_HVG` jointly).

---

### §3.2 Hypothesis `H-J1` — Batch-Index Landscape, Winner's-Curse Control, & Floor-Effect Guard (`J2`, `J4`)

**Scientific Question:** Holding matched target control cells $X_{\mathrm{ctrl}}$ and frozen weights fixed across **all 56 trained indices (`0..55`, never dropped post-inference per `J4`)**, how much does zero-shot `r_orth_cal` vary, and does true between-index spread exceed both perturbation-sampling noise and Round 16's practical threshold $\delta = 0.04$?

1. **Full Descriptive Landscape Outputs:** Report all 56 point estimates $\{r_{t,m}(b)\}_{b=0}^{55}$ and individual bootstrap SEs, `min`, `q10`, `q25`, `median`, `mean`, `q75`, `q90`, `max` ($\max_b r(b)$), `IQR`, `L1_ens56` score ($r(\texttt{L1\_ens56})$), cross-fitted oracle score ($\widehat{r}_{\mathrm{CF\text{-}Oracle}, 0}$), observed between-index SD $s_{\mathrm{obs}} = \sqrt{\frac{1}{55}\sum_{b=0}^{55}(r(b)-\bar{r})^2}$, raw range $R_{\mathrm{obs}} = \max_b r(b) - \min_b r(b)$, and exact rank & percentile of official fallback `b = 0` (`fixed_index_zero`, `_infer.py:688`).
2. **Winner's-Curse-Free Co-Estimator 1 — Perturbation Half-Split Cross-Fitted Spread ($\widehat{\Delta}_{\mathrm{CF}, 0}$):**
   - **Locked Primary Split (`s = 0`, deterministic seed `20261701`):** Sort the $N_t$ evaluated perturbations canonically by target symbol, permute using `np.random.default_rng(20261701)` stratified by `half_valid` status whenever any `half_valid == False` perturbations exist (exact alternating allocation after within-stratum permutation), and partition into disjoint halves $\mathcal{P}_1$ ($\lfloor N_t / 2 \rfloor$) and $\mathcal{P}_2$ ($\lceil N_t / 2 \rceil$).
   - **Cross-Fitting Step:**
     - Select best and worst indices strictly on $\mathcal{P}_1$: $\hat{b}^*_{1} = \arg\max_{b \in 0..55} r_{\mathcal{P}_1}(b)$, $\hat{b}^{\dagger}_{1} = \arg\min_{b \in 0..55} r_{\mathcal{P}_1}(b)$ (lowest-index tie-break). Evaluate strictly out-of-sample on $\mathcal{P}_2$: $\Delta_{1\to 2} = r_{\mathcal{P}_2}(\hat{b}^*_{1}) - r_{\mathcal{P}_2}(\hat{b}^{\dagger}_{1})$.
     - Symmetrically select $(\hat{b}^*_{2}, \hat{b}^{\dagger}_{2})$ strictly on $\mathcal{P}_2$ and evaluate strictly out-of-sample on $\mathcal{P}_1$: $\Delta_{2\to 1} = r_{\mathcal{P}_1}(\hat{b}^*_{2}) - r_{\mathcal{P}_1}(\hat{b}^{\dagger}_{2})$.
     - Primary winner's-curse-free spread & cross-fitted oracle score:
       $$\widehat{\Delta}_{\mathrm{CF}, 0}(t, m) = \frac{1}{2}\!\left(\Delta_{1\to 2} + \Delta_{2\to 1}\right), \qquad \widehat{r}_{\mathrm{CF\text{-}Oracle}, 0}(t, m) = \frac{1}{2}\!\left(r_{\mathcal{P}_2}(\hat{b}^*_1) + r_{\mathcal{P}_1}(\hat{b}^*_2)\right).$$
   - **Exact Within-Half Bootstrap Inference (`n_boot = 4,000`):** Conditional on $(\hat{b}^*_{1}, \hat{b}^{\dagger}_{1})$ fixed from $\mathcal{P}_1$, draw `4,000` bootstraps inside $\mathcal{P}_2$ (`seed = 20261701`); conditional on $(\hat{b}^*_{2}, \hat{b}^{\dagger}_{2})$ fixed from $\mathcal{P}_2$, draw `4,000` independent bootstraps inside $\mathcal{P}_1$ (`seed = 20261701 + 10000`); form exact paired half-averaged bootstrap draws for both $\widehat{\Delta}_{\mathrm{CF}, 0}^{(u)}$ and $\widehat{r}_{\mathrm{CF\text{-}Oracle}, 0}^{(u)}$ to obtain `98.333%`, Holm-adjusted, and `95.0%` percentile CIs and empirical two-sided $p$-values.
   - **Multi-Split Stability Check (`R_split = 20`, Secondary Diagnostic):** Repeat point estimation across `seed = 20261701 + s` (`s = 0..19`) and report mean, median, min, max, and IQR of $\widehat{\Delta}_{\mathrm{CF}, s}$ (only locked `s = 0` enters primary hypothesis testing).
3. **Winner's-Curse-Free Co-Estimator 2 — Noise-Debiased Between-Index SD ($\hat{\tau}_{\mathrm{batch}}$) & Centered-Bootstrap Homogeneity Test ($p_{\mathrm{var}}$):**
   - From the `4000 × 56` joint bootstrap matrix $R^{\mathrm{boot}}_{t,m}$ (all $N_t$ perturbations), with $C_{56} = I_{56} - \frac{1}{56}\mathbf{1}\mathbf{1}^\top$ and $S_{\mathrm{obs}}^2 = \frac{1}{55}\|C_{56}\hat{r}\|_2^2$:
     $$V_{\mathrm{noise}} = \frac{1}{55}\sum_{b=0}^{55} \operatorname{Var}_{\mathrm{boot}}\!\big(R^{\mathrm{boot}}_{u, b} - \bar{R}^{\mathrm{boot}}_{u, \cdot}\big), \qquad \hat{\tau}_{\mathrm{batch}}(t, m) = \sqrt{\max\!\big(0,\; S_{\mathrm{obs}}^2 - V_{\mathrm{noise}}\big)},$$
     $$Z^{(u)} = R^{\mathrm{boot}}_{u, \cdot} - \bar{R}^{\mathrm{boot}}_{\cdot}, \qquad p_{\mathrm{var}}(t, m) = \frac{1 + \sum_{u=1}^{4000} \mathbb{I}\!\big(\frac{1}{55}\|C_{56} Z^{(u)}\|_2^2 \ge S_{\mathrm{obs}}^2\big)}{4001}.$$
4. **Mandatory Floor-Effect Detectable-Signal Predicate (`J2`):**
   - Low-signal lines (such as R16 `cd4t_rest`) can look batch-invariant simply because no batch index produces positive signal above `Transfer_Mean`. Define the locked detectable-signal predicate:
     $$\mathrm{has\_detectable\_signal}(t, m) \;\iff\; \Big(\mathrm{CI}_{98.333\%}^{\mathrm{lo}}\big(\widehat{r}_{\mathrm{CF\text{-}Oracle}, 0}(t, m)\big) > 0.0\Big) \;\lor\; \Big(\mathrm{CI}_{98.333\%}^{\mathrm{lo}}\big(r_{t,m}(\texttt{L1\_ens56})\big) > 0.0\Big)$$
     (requiring `valid_frac >= 0.95` on the qualifying estimator).
   - **Mandatory Reporting (`J2`):** Report **$\max_b r_{t,m}(b)$, $r_{t,m}(\texttt{L1\_ens56})$ (with `98.333%` CI), and $\widehat{r}_{\mathrm{CF\text{-}Oracle}, 0}$ (with `98.333%` CI)** alongside **every** `H-J1` call.
5. **Locked 5-Way Decision Rule for `H-J1` (`J2` Updated; Anchored to $\delta = 0.04$):**
   Let $\alpha_{\mathrm{adj}}$ denote the locked Holm–Bonferroni multiplicity-adjusted level (`0.01667` for single-line base check, step-down $\alpha_{(j)}$ across $3 K_{\mathrm{conf}}$ confirmatory tests):
   1. **`INCONCLUSIVE` (Validity Fail):** `valid_frac < 0.95` on any of the 56 indices or on either cross-fit half, OR `rho_ab_orth <= 0` or non-finite CIs;
   2. **`BATCH_DEPENDENCE_PRACTICALLY_MEANINGFUL` (H-J1 Confirmed):** ALL of:
      - $p_{\mathrm{var}} < \alpha_{\mathrm{adj}}$ **AND** $\mathrm{CI}_{1-\alpha_{\mathrm{adj}}}^{\mathrm{lo}}(\widehat{\Delta}_{\mathrm{CF}, 0}) > 0.0$;
      - $\widehat{\Delta}_{\mathrm{CF}, 0} > 0.04$ **AND** $2\,\hat{\tau}_{\mathrm{batch}} > 0.04$ (i.e., $\hat{\tau}_{\mathrm{batch}} > 0.02$);
   3. **`BATCH_DEPENDENCE_DETECTED_SUB_DELTA`:** $p_{\mathrm{var}} < \alpha_{\mathrm{adj}}$ AND $\mathrm{CI}_{1-\alpha_{\mathrm{adj}}}^{\mathrm{lo}}(\widehat{\Delta}_{\mathrm{CF}, 0}) > 0.0$, yet $\widehat{\Delta}_{\mathrm{CF}, 0} \le 0.04$ OR $\hat{\tau}_{\mathrm{batch}} \le 0.02$ *(mandatory wording: detected above noise, but smaller than practical threshold $\delta = 0.04$)*;
   4. **`INCONCLUSIVE_LOW_SIGNAL` (`J2` Floor-Effect Guard):** NOT statistically detected ($p_{\mathrm{var}} \ge \alpha_{\mathrm{adj}}$ OR $\mathrm{CI}_{1-\alpha_{\mathrm{adj}}}^{\mathrm{lo}}(\widehat{\Delta}_{\mathrm{CF}, 0}) \le 0.0$), **AND** $\mathrm{has\_detectable\_signal}(t, m) == \mathrm{False}$ *(neither cross-fitted oracle nor `L1_ens56` has `98.333%` CI lower bound $> 0$; invariance CANNOT be awarded on a zero-signal floor)*;
   5. **`BATCH_INVARIANT_WITHIN_DELTA`:** NOT statistically detected, **AND** $\mathrm{has\_detectable\_signal}(t, m) == \mathrm{True}$, **AND** two-sided `95.0%` upper bound $\mathrm{CI}_{95\%}^{\mathrm{hi}}(\widehat{\Delta}_{\mathrm{CF}, 0}) < 0.04$ AND $\hat{\tau}_{\mathrm{batch}} \le 0.02$;
   6. **`INCONCLUSIVE` (Insufficient Precision):** All remaining cases (detectable signal present, heterogeneity not proven, yet `95.0%` upper bound $\ge 0.04$ or $\hat{\tau}_{\mathrm{batch}} > 0.02$; interpreted strictly as insufficient precision, NEVER as batch invariance proven).

---

### §3.3 Hypothesis `H-J2` — `STATE - FCT` Benchmark Call Instability Across All 56 Batch Indices (with `J2` Floor-Effect Qualifier)

**Scientific Question:** Using the exact `PREREG_I §3.2` `H-I2` decision rule (`STATE - FCT`, two-sided `97.5%` joint perturbation bootstrap CI), what fraction of all 56 trained batch indices (`0..55`, never dropped per `J4`) produce a **different** benchmark verdict from official fallback `fixed_index_zero` (`b = 0`, `_infer.py:688`)?

1. **Exact Per-Index `H-I2` Call:** For every $b \in \{0, \dots, 55\}$, compute $\Delta^{\mathrm{vsFCT}}_{t,m}(b) = r_{\mathrm{orth\_cal}}(\mathrm{model}_{m, b}) - r_{\mathrm{orth\_cal}}(\mathrm{FCT})$ using the locked `4,000` joint perturbation bootstraps (`seed = 20261701`):
   - `MODEL_SUPERIOR_TO_FCT` iff `valid_frac >= 0.95` and `ci975_lo > 0.0`;
   - `FCT_SUPERIOR_TO_MODEL` iff `valid_frac >= 0.95` and `ci975_hi < 0.0`;
   - `INDISTINGUISHABLE` otherwise (including exact `0.0` ties or `valid_frac < 0.95`).
2. **Primary `H-J2` Quantity & Companion Bonferroni-56 Guard:**
   - Baseline call under official fallback: $c_0 = \mathrm{call}_{t,m}(0)$.
   - **Call-instability fraction:** $f_{\mathrm{flip}}(t, m) = \frac{1}{56}\sum_{b=0}^{55} \mathbb{I}\!\big(\mathrm{call}_{t,m}(b) \neq c_0\big)$, reported with `(n_MODEL_SUPERIOR, n_INDISTINGUISHABLE, n_FCT_SUPERIOR)`, index lists, and `has_full_reversal = (n_MODEL_SUPERIOR >= 1) AND (n_FCT_SUPERIOR >= 1)`.
   - **Companion Bonferroni-56 Hardened Stability Check:** Report $f_{\mathrm{flip}}^{\mathrm{Bonf56}}(t, m)$ under both `z = 3.5117` normal-calibrated `99.9554%` CIs ($\alpha = 0.025 / 56$) AND `n_boot = 20,000` empirical quantile CIs side-by-side.
3. **Locked 5-Way Audit Classification for `H-J2` (`J2` Floor-Effect Qualifier Applied to `CALL_STABLE`):**
   1. **`CALL_UNSTABLE_SIGN_REVERSAL`:** `has_full_reversal == True`;
   2. **`CALL_UNSTABLE_SUBSTANTIAL`:** NOT full reversal, AND $f_{\mathrm{flip}}(t, m) \ge 5/56$ (`>= 8.93%`);
   3. **`CALL_UNSTABLE_MINOR`:** $1/56 \le f_{\mathrm{flip}}(t, m) < 5/56$ (`1–4` of 56 indices disagree with `b = 0`);
   4. **`CALL_STABLE`:** $f_{\mathrm{flip}}(t, m) = 0$ (all `56/56` indices agree with `b = 0`) **AND** $\mathrm{has\_detectable\_signal}(t, m) == \mathrm{True}$ (`CF-Oracle` OR `L1_ens56` `98.333%` CI lower bound $> 0$);
   5. **`INCONCLUSIVE_LOW_SIGNAL` (`J2`):** $f_{\mathrm{flip}}(t, m) = 0$ **AND** $\mathrm{has\_detectable\_signal}(t, m) == \mathrm{False}$ *(zero flips occurred only because the model has no detectable signal at any index; stability on a null floor cannot be interpreted as genuine call robustness)*.
   - **Mandatory Reporting:** Always report `FCT r_orth_cal` (`PREREG_I §3.2` Anti-Masking Rule), $\max_b r_{t,m}(b)$, and $r_{t,m}(\texttt{L1\_ens56})$ alongside every `H-J2` call.

---

### §3.4 Hypothesis `H-J3` — Pre-Specified Label-Free Batch-Index Rules vs. `fixed_index_zero` & Cross-Fitted Oracle (Ruling 3 Locked)

Because shipped `arc-state` checkpoints contain no per-`gem_group` training embeddings (`§0.1`), `H-J3` locks two strictly self-contained, zero-perturbation-label rules computable solely from shipped weights + the target line's own unlabeled NTC control cells:

1. **Rule `L1_ens56` (Uniform 56-Index Decoded-Space Ensemble — Co-Primary `H-J3a`):**
   - Uniform average of decoded predictions across all 56 trained indices (`0..55`) in clipped `[0, 14]` space before subtracting `'non-targeting'` ($\hat{\Delta}^{(m, \mathrm{L1})}_t(p) = \frac{1}{56}\sum_{b=0}^{55} \hat{\Delta}^{(m, b)}_t(p)$). Zero perturbation labels, zero free parameters, strictly inside the trained norm shell $\|e_b\|_2 \approx 11.1$.
2. **Rule `L2_ntc_recon_match` (Unlabeled NTC Control-Profile Reconstruction Index — Co-Primary `H-J3b`, Ruling 3):**
   - **Primary Selection Criterion (Locked):** Compute observed unlabeled NTC mean expression $\bar{Y}^{\mathrm{obs}}_{\mathrm{ctrl}}(t)$ on $G_{\mathrm{eval}}(t)$ using `official_median_log1p` over the sampled `n_ctrl` unperturbed NTC cells (zero perturbed cells, zero perturbation labels). For every $b \in \{0, \dots, 55\}$, compute the criterion:
     $$C_{\mathrm{NTC}}(b) = \operatorname{Pearson}_{g \in G_{\mathrm{eval}}(t)}\!\Big(\overline{\hat{Y}}^{(m, b)}_t(\text{"non-targeting"}),\; \bar{Y}^{\mathrm{obs}}_{\mathrm{ctrl}}(t)\Big), \qquad \hat{b}_{\mathrm{NTC}}(t, m) = \arg\max_{b \in \{0,\dots,55\}} C_{\mathrm{NTC}}(b)$$
     (lowest integer index tie-break).
   - **Mandatory Top-1 vs. Top-2 Selection Margin Disclosure (Ruling 3):** To expose near-degenerate selections transparently, **ALWAYS report:**
     - `top1_idx` ($\hat{b}_{\mathrm{NTC}}$) and `criterion_top1` ($C_{\mathrm{NTC}}^{(1)}$);
     - `top2_idx` (runner-up index) and `criterion_top2` ($C_{\mathrm{NTC}}^{(2)}$);
     - **Exact margin `criterion_margin_top1_minus_top2 = criterion_top1 - criterion_top2`** plus full 56-index criterion `min / median / IQR / max` and exact rank of `b = 0`.
   - **Secondary `L2` Diagnostics (Ruling 3):** Log (i) negative L2 distance criterion $-\|\overline{\hat{Y}}^{(m, b)}_t(\text{"non-targeting"}) - \bar{Y}^{\mathrm{obs}}_{\mathrm{ctrl}}(t)\|_2$ (with its own `top1_idx`, `top2_idx`, and `top1 - top2` gap) and (ii) per-fold index selection $\hat{b}_{\mathrm{NTC}}^{(f)}$ (`f = 1..4`) and agreement across folds.
3. **Locked Paired Contrasts & Decision Rules for Each Rule $L \in \{\texttt{L1\_ens56}, \texttt{L2\_ntc\_recon\_match}\}$:**
   - **Primary Contrast vs. `fixed_index_zero` (`b = 0`):** $\Delta_{L \text{ vs } 0}(t, m) = r_{\mathrm{orth\_cal}}(L) - r_{\mathrm{orth\_cal}}(b = 0)$ under the locked Holm–Bonferroni multiplicity-adjusted level $\alpha_{\mathrm{adj}}$ (`4,000` joint bootstraps, `seed = 20261701`):
     - `RULE_SUPERIOR_TO_ZERO` iff `valid_frac >= 0.95` AND `ci_adj_lo > 0.0` (annotated `> 0.04` practical vs. `<= 0.04` sub-delta);
     - `ZERO_SUPERIOR_TO_RULE` iff `valid_frac >= 0.95` AND `ci_adj_hi < 0.0`;
     - `EXACT_TIE_SELECTED_ZERO` if $L = \texttt{L2\_ntc\_recon\_match}$ and $\hat{b}_{\mathrm{NTC}} == 0$;
     - `INDISTINGUISHABLE` otherwise.
   - **Out-of-Sample Contrast vs. Cross-Fitted Oracle:** On locked split `s = 0` ($\mathcal{P}_1, \mathcal{P}_2$), compute $\Delta_{L \text{ vs CF-Oracle}} = \frac{1}{2}(r_{\mathcal{P}_2}(L) + r_{\mathcal{P}_1}(L)) - \widehat{r}_{\mathrm{CF\text{-}Oracle}, 0}$ with within-half joint bootstrap CI.
   - **Benchmark Impact:** Report `L vs FCT` (`PREREG_I §3.2` `97.5%` CI) alongside `b = 0 vs FCT`.
4. **Secondary Neutral-Batch Comparators (`PREREG_I §4.4` Continuity):** `bypass_zero`, `mean_batch_emb`, `multi_batch_ens8` (`{0,7,14,21,28,35,42,49}`), and weight-only `medoid_idx`.

---

### §3.5 Pre-Registered Secondary Mechanism Tie-In (`J5`) — Cross-Fitted `G4` 2×2 Decomposition & Cross-Line 56-Index Rank Correlation

Per Master Controller Review R1 (`J5`), Round 17 locks two secondary mechanism analyses on every fresh line (for BOTH `STATE_SE` and `STATE_HVG`):

1. **`J5a` — Per-Index `G4` 2×2 (`X_ctrl × batch`) Decomposition at Cross-Fitted Best/Worst Indices (Plus Pre-Fixed Continuity Pairs):**
   - **Cross-Fitted Best vs. Worst Pair (Zero Selection Leakage, `J5` Primary Mechanism Pair):** Using locked split `s = 0` ($\mathcal{P}_1, \mathcal{P}_2$ from `§3.2`), select the best and worst indices strictly on $\mathcal{P}_1$: $(B_{\mathrm{best}} = \hat{b}^*_1,\; B_{\mathrm{worst}} = \hat{b}^\dagger_1)$. Then run the exact R16-style `G4` 2×2 factorial (`(M, B_best)`, `(M, B_worst)`, `(S, B_best)`, `(S, B_worst)` AND symmetric reversed `(B_worst, B_best)`) **evaluated strictly on disjoint $\mathcal{P}_2$ ONLY** (reporting $\Delta_{\mathrm{batch}}, \Delta_{\mathrm{ctrl}}, \Delta_{\mathrm{both}}, \Delta_{\mathrm{int}}$, and `batch_share = |delta_batch| / (|delta_ctrl| + |delta_batch|)` under `r_orth_cal` and `r_dorth_cal`, `4,000` bootstraps strictly inside $\mathcal{P}_2$). *(A symmetric $\mathcal{P}_2 \to \mathcal{P}_1$ check may be logged alongside as an auxiliary diagnostic, with $(\hat{b}^*_1, \hat{b}^\dagger_1)$ evaluated on $\mathcal{P}_2$ as the locked `J5` readout.)*
   - **Pre-Fixed Continuity Pairs (Evaluated on All $N_t$ Perturbations):** Also retain `(B_M = 0, B_S = 17)` (with reversed `(17, 0)`) for `PREREG_I §4.4` continuity and `(B_M = 0, B_S = 35)` (with reversed `(35, 0)`) to test transfer of R16 `a549_ko`'s `b0 vs b35` separation.
   - **Mismatched Control Strata (`X_ctrl = S`):**
     - `seen_by_checkpoint` ($n = 2$): `{nadig_hepg2, nadig_jurkat}` (`3/4`-checkpoint overlap disclosed);
     - `unseen_by_checkpoint` ($n = 3$–`4`): `{a549_ko, dld1, cd4t_rest, feng_ipsc}` excluding any pool sharing `canonical_family` with target $t$ (`cd4t_rest` reported separately as `same_lineage_resting_comparator` on `cd4t_stim48`; `a549_ko` excluded on `a549_myllia` and `a549_ko_cas9pos`; all 4 pools foreign on `gse295214_h1` and `pfizer_telohaec_basal`).
2. **`J5b` — Cross-Line Spearman Rank Correlation of the 56 Per-Index Scores with Index-Label Permutation Null:**
   - **Definition:** For each model $m \in \{\texttt{STATE\_SE}, \texttt{STATE\_HVG}\}$ and every unordered pair of fresh lines $(t_1, t_2)$ (computed separately for: [i] all pairs among confirmatory units `C_conf`; [ii] `pfizer_telohaec_basal` vs. `pfizer_telohaec_treated` within-family pair; [iii] all fresh lines with $N \ge 16$; and [iv] cross-model `STATE_SE` vs. `STATE_HVG` on the identical line $t$), compute the exact Spearman rank correlation over the full selection-free 56-index vector:
     $$\rho_{\mathrm{Spearman}}(t_1, t_2; m) = \operatorname{Spearman}_{b \in \{0..55\}}\!\Big(r_{t_1, m}(b),\; r_{t_2, m}(b)\Big).$$
   - **Exact Index-Label Permutation Null:** Draw `n_perm = 10,000` independent uniform permutations $\pi_k \in S_{56}$ (`np.random.default_rng(20261701)`), compute $\rho^{(k)} = \operatorname{Spearman}(r_{t_1, m}(b),\; r_{t_2, m}(\pi_k(b)))$, and report the exact two-sided permutation $p$-value:
     $$p_{\mathrm{perm}}(t_1, t_2; m) = \frac{1 + \sum_{k=1}^{10000} \mathbb{I}\!\big(|\rho^{(k)}| \ge |\rho_{\mathrm{Spearman}}(t_1, t_2; m)|\big)}{10001}$$
     (plus Benjamini–Hochberg FDR across pairs within each model, and the top-5 / bottom-5 shared index overlap counts).
   - **Locked Mechanistic Interpretation (`J5`):**
     - **Strong positive shared ordering across independent lineages** (e.g., $\rho_{\mathrm{Spearman}} > 0$, $p_{\mathrm{perm}} < 0.05$ consistently across independent confirmatory units) points to a **line-independent property of the batch embeddings themselves** (e.g., norm, training-sample quality, or global alignment), indicating that **a fixed better default batch index may exist**.
     - **Near-zero or idiosyncratic ordering across independent lineages** points to **line-specific `(control × batch)` interaction**, indicating that **no universal fixed single-index default exists** (favoring ensembling such as `L1_ens56` or context-adaptive selection over replacing `b = 0` with another fixed index).

---

## §4. Power and Precision Statement (`J3` — Anchored to R16 SEs + Parametric $N_{\mathrm{eval}}$ Functions for Prospective Lines)

### §4.1 Exact Round 16 Discovery Planning Inputs & Heterogeneity Caveat (`J3`)
All planning standard errors below are read directly from Round 16 executed outputs (`round16/results/summary_r16.json` and `round16/results/sensitivity/summary_r16.json`, `n_boot = 4,000`):

| R16 Discovery Line | $N_{\mathrm{eval}}$ | `rho_ab_orth` | `b=0` `boot_se` | Paired `b0 - b17` `boot_se` ($\mathrm{SE}_{\mathrm{pair}}$) | Scaled Constant $c = \mathrm{SE}_{\mathrm{pair}}\sqrt{N}$ | `STATE_SE - FCT` `boot_se` | Observed 4-Index Range / SD |
|---|---:|---:|---:|---:|---:|---:|---:|
| `a549_ko` | `403` | `0.8271` | `0.00735` | **`0.00748`** | **`0.1502`** *(high-$\rho^\perp$ anchor)* | `0.01172` | **`0.1563`** / `0.0702` |
| `dld1` | `75` | `0.2931` | `0.00962` | **`0.00889`** | `0.0770` | `0.01778` | **`0.0598`** / `0.0254` |
| `cd4t_rest` | `245` | `0.2062` | `0.01174` | **`0.01118`** | **`0.1750`** *(low-$\rho^\perp$ conservative anchor)* | `0.01223` | **`0.0374`** / `0.0161` |
| `feng_ipsc` | `88` | `0.6177` | `0.01297` | **`0.01466`** | `0.1375` | `0.01795` | **`0.0539`** / `0.0247` |

> **Explicit Conditional-Power & Discovery-Heterogeneity Disclosure (`J3`):** All power figures below are strictly **conditional on a true cross-fitted spread of $\Delta = 0.04$ or $\Delta = 0.06$** (or $\tau_{\mathrm{batch}} = 0.02$). **Crucially, Round 16 discovery 4-index ranges differed widely across cell lines — from `0.1563` (`a549_ko`) down to `0.0598` (`dld1`), `0.0539` (`feng_ipsc`), and `0.0374` (`cd4t_rest`) — so this design guards against false positives rather than against low-signal lines; the `J2` floor-effect rule (`INCONCLUSIVE_LOW_SIGNAL`) handles low-signal lines explicitly.**

### §4.2 Precision & Power for Local Confirmatory Unit 1 (`cd4t_stim48`, $N_{\mathrm{eval}} = 401$) & Secondary Fresh Lines
As derived in v1 (`§4.2`), algebraic cancellation across the two independent half-splits yields $\mathrm{SE}(\widehat{\Delta}_{\mathrm{CF}, 0}) \approx \mathrm{SE}_{\mathrm{pair}}$, and for the 56-index variance-component test (`df = 55`), $55 S_{\mathrm{obs}}^2 / \sigma_e^2 \sim (1 + \tau_{\mathrm{batch}}^2 / \sigma_e^2)\,\chi^2_{55}$ with $\sigma_e = \mathrm{SE}_{\mathrm{pair}} / \sqrt{2}$.

- **`cd4t_stim48` ($N = 401$):**
  - Scaled from `cd4t_rest` ($\mathrm{SE}_{\mathrm{pair}} = 0.00874$): `98.33%` CI half-width `±0.0209` (`K=1`, $z=2.3940$) / worst-case `K=3` Holm `99.44%` half-width **`±0.0242`** ($z=2.7729$); conditional power at $\Delta = 0.04$: **`98.5%`** (`K=1`) / **`96.4%`** (`K=3` Holm step-1); conditional power at $\Delta = 0.06$: **`>99.9%`**; variance-test power at $\tau = 0.02$: **`>99.9%`**.
  - Conservative unscaled `cd4t_rest` ($\mathrm{SE}_{\mathrm{pair}} = 0.01118$): `K=1` half-width `±0.0268` (`88.2%` power at $\Delta=0.04$) / `K=3` Holm step-1 half-width **`±0.0310`** (`79.0%` power at $\Delta=0.04$, `99.5%` at $\Delta=0.06$, `>99.9%` variance-test power).
- **Secondary Fresh Lines (`a549_myllia` `N = 20`, `hek293_xucao` `N = 16`):** $\mathrm{SE}_{\mathrm{pair}} \approx 0.030$ (`98.33%` half-width `±0.0718`, `14.4%` power at $\Delta = 0.04$), confirming strictly descriptive secondary status.

### §4.3 Parametric Precision & Power Functions as a Function of $N_{\mathrm{eval}}$ for Prospective Units (`gse295214_h1` & `pfizer_telohaec_basal`, `J3`)
Because counts-only QC for `gse295214_h1` (Task 1) and `pfizer_telohaec_basal` (post-lock) occurs prior to any model inference, precision and conditional power are pre-registered **mechanically as closed-form functions of realized $N_{\mathrm{eval}}$**, bracketed by the two empirical R16 scaling envelopes:
$$\mathrm{SE}_{\mathrm{pair}}^{\mathrm{hi}}(N) = \frac{0.1502}{\sqrt{N}} \quad (\text{high-}\rho^\perp\ \texttt{a549\_ko}\text{-like}), \qquad \mathrm{SE}_{\mathrm{pair}}^{\mathrm{lo}}(N) = \frac{0.1750}{\sqrt{N}} \quad (\text{conservative low-}\rho^\perp\ \texttt{cd4t\_rest}\text{-like}),$$
$$\mathrm{HW}_{K}(N) = z_{K}\,\mathrm{SE}_{\mathrm{pair}}(N) \quad (z_1 = 2.3940,\; z_2 = 2.6383,\; z_3 = 2.7729), \qquad \mathrm{MDE}_{80\%, K}(N) = (z_K + 0.8416)\,\mathrm{SE}_{\mathrm{pair}}(N).$$

| Realized $N_{\mathrm{eval}}$ (From Counts-Only QC) | $\mathrm{SE}_{\mathrm{pair}}(N)$ `[hi .. lo]` | `98.33%` (`K=1`) Half-Width `[hi .. lo]` | Worst-Case `K=3` Holm (`99.44%`) Half-Width `[hi .. lo]` | `80%`-Power MDE (`K=3` Holm) `[hi .. lo]` | Cond. Power at $\Delta = 0.04$ (`K=3` Holm) `[hi .. lo]` | Cond. Power at $\Delta = 0.06$ (`K=3` Holm) `[hi .. lo]` | Variance-Test Power at $\tau = 0.02$ (`K=3` Holm, conservative `lo`) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| **`50`** *(floor)* | `0.02124 .. 0.02475` | `±0.0508 .. ±0.0592` | `±0.0589 .. ±0.0686` | `0.0768 .. 0.0895` | `18.7% .. 12.4%` | `52.1% .. 36.4%` | `97.1%` |
| **`75`** | `0.01734 .. 0.02021` | `±0.0415 .. ±0.0484` | `±0.0481 .. ±0.0560` | `0.0627 .. 0.0730` | `32.1% .. 21.4%` | `75.4% .. 57.8%` | `99.9%` |
| **`100`** | `0.01502 .. 0.01750` | `±0.0359 .. ±0.0419` | `±0.0416 .. ±0.0485` | `0.0543 .. 0.0633` | `45.7% .. 31.3%` | `88.9% .. 74.4%` | `>99.9%` |
| **`150`** | `0.01226 .. 0.01429` | `±0.0294 .. ±0.0342` | `±0.0340 .. ±0.0396` | `0.0443 .. 0.0516` | `68.8% .. 51.1%` | `98.3% .. 92.3%` | `>99.9%` |
| **`200`** | `0.01062 .. 0.01237` | `±0.0254 .. ±0.0296` | `±0.0294 .. ±0.0343` | `0.0384 .. 0.0447` | `84.0% .. 67.7%` | `99.8% .. 98.1%` | `>99.9%` |
| **`300`** | `0.00867 .. 0.01010` | `±0.0208 .. ±0.0242` | `±0.0240 .. ±0.0280` | `0.0313 .. 0.0365` | `96.7% .. 88.2%` | `>99.9% .. 99.9%` | `>99.9%` |
| **`400`** | `0.00751 .. 0.00875` | `±0.0180 .. ±0.0209` | `±0.0208 .. ±0.0243` | `0.0271 .. 0.0316` | `99.5% .. 96.4%` | `>99.9% .. >99.9%` | `>99.9%` |
| **`600`** | `0.00613 .. 0.00714` | `±0.0147 .. ±0.0171` | `±0.0170 .. ±0.0198` | `0.0222 .. 0.0258` | `>99.9% .. 99.8%` | `>99.9% .. >99.9%` | `>99.9%` |

- **Mechanical Pre-Inference Fill-In Slots (Populated Strictly From Counts-Only QC Before Any Model Inference):**
  - `gse295214_h1`: realized $N_{\mathrm{eval}} = \text{[TBD\_POST\_R16\_G3\_QC]} \implies$ plug directly into formulas above to record pre-inference $\mathrm{SE}_{\mathrm{pair}}$, $\mathrm{HW}_{K}$, $\mathrm{MDE}_{80\%}$, and conditional power at $\Delta \in \{0.04, 0.06\}$.
  - `pfizer_telohaec_basal`: realized $N_{\mathrm{eval}} = \text{[TBD\_POST\_LOCK\_QC]} \implies$ plug directly into formulas above before any model inference.
  - Note: Under worst-case `K = 3` Holm–Bonferroni, **`>= 80%` conditional power for $\Delta = 0.04$** is achieved at $N_{\mathrm{eval}} \ge 185$ (`hi`) to $250$ (`lo`), while the 56-index variance-component test exceeds **`97%` power** across the entire qualifying range $N_{\mathrm{eval}} \ge 50$.

---

## §5. Compute Plan & Three-Party Blinding/Verification Plan (`J6`, `A-T3` Post-R16)

### §5.1 Local Compute & Memory Budget (Apple Silicon Mac, `128 GB` Unified RAM, `MPS float32`)
Anchored to Round 16 verified `MPS float32` latencies (`SE-600M = 5.5 min` per `512`-cell pool, peak RSS `5.70 GB`; `ST-SE = 12.65 ms/pert/fold`; `ST-HVG = 10.21 ms/pert/fold`, peak RSS `2.35 GB`). Per **Ruling 2**, **`wtc11_i3neuron` (`N = 1`) is omitted from inference**:

| Line Key (Locked Role) | $N_{\mathrm{eval}}$ (+1 NTC) | `SE-600M` (`MPS fp32`) | `ST-SE` (`56 idx × 4 folds` + `G4`) | `ST-HVG` (`56 idx × 4 folds` + `G4`) | Total Compute |
|---|---:|---:|---:|---:|---:|
| **`cd4t_stim48`** (Confirmatory Unit 1) | `401` (`402`) | `5.5 min` | `19.0 + 7.0 = 26.0 min` | `15.3 + 5.7 = 21.0 min` | **`52.5 min`** |
| **`gse295214_h1`** (Confirmatory Unit 2, Track 2A) | $N_{\mathrm{H1}}$ (est. `~100–600`) | `5.5 min` *(reusable if R16 G3 uses same pool? No — R17 uses fresh seed `20261701`!)* | `0.065 × N_H1 min` | `0.052 × N_H1 min` | `5.5 + 0.117 × N_H1 min` |
| **`pfizer_telohaec_basal`** (Confirmatory Unit 3, Track 2B) | $N_{\mathrm{TH}}$ (est. `~100–500`) | `5.5 min` | `0.065 × N_TH min` | `0.052 × N_TH min` | `5.5 + 0.117 × N_TH min` |
| `pfizer_telohaec_treated` (Secondary Within-Family Replication) | $N_{\mathrm{TH,tr}}$ | `5.5 min` | `0.065 × N_TH,tr min` | `0.052 × N_TH,tr min` | `5.5 + 0.117 × N_TH,tr min` |
| `a549_myllia` (Secondary Fresh) | `20` (`21`) | `5.5 min` | `1.3 min` | `1.1 min` | `7.9 min` |
| `hek293_xucao` (Secondary Fresh) | `16` (`17`) | `5.5 min` | `1.1 min` | `0.9 min` | `7.5 min` |
| `telohaec`, `calu3`, `mac_ctrl`, `mac_lps` (4 Micro-N, Ruling 2) | `9, 8, 3, 4` | `22.0 min` | `1.9 min` | `1.6 min` | `25.5 min` |
| `a549_ko_cas9pos` (R16 Subset Sensitivity, Ruling 2) | `403` (`404`) | `5.5 min` | `26.1 min` | `21.1 min` | `52.7 min` |
| `wtc11_i3neuron` (`N = 1`, **OMITTED per Ruling 2**) | — | `0 min` | `0 min` | `0 min` | **`0 min`** *(omitted)* |
| **Local Subtotal (All 8 Executed Local Lines)** | — | **`44.0 min`** | **`56.4 min`** | **`45.7 min`** | **`2.44 hrs`** |
| **Prospective Track 2 (`gse295214_h1` + Pfizer Basal & Treated)** | — | **`16.5 min`** | `~1.0–2.0 hrs` | `~0.8–1.6 hrs` | **`~2.1–3.9 hrs`** |

- **Strict Operational Guardrails:** (1) Peak RSS `~5.70 GB` (`< 4.5%` of `128 GB`); (2) **Zero concurrent MPS execution** across Main Line and `A-T3`; (3) **10-minute MPS stall fallback** switches both parties identically to `CPU float32`; (4) `MPS float32` ONLY (`bf16`/`fp16` prohibited due to Metal assertion verified at R16 Gate 2).

### §5.2 Blinding & Release Gate (`J6`) and Staged Three-Party Verification Plan (`A-T3` Post-R16)
- **Strict R16-to-R17 Blinding & Release Sequence (`J6`):**
  1. `round17/` contains Round 16 main-line discovery numbers (`§0.2`, `§4.1`). Therefore, **`round17/` MUST remain uncommitted and hidden from `A-T3` until Master Controller has verified `A-T3`'s Round 16 pass 2.**
  2. `A-T3` joins Round 17 **only after Master Controller formally releases it from Round 16**; `A-T3` then independently implements `round17/audit/` strictly from the locked `PREREG_J` alone (zero shared Round 17 analysis code).
  3. Main Line's execution engine (`round17/src/`) enforces a **triple hard gate** refusing real-data execution unless `PREREG_J_COMMIT.txt`, `CODE_FREEZE_R17.txt`, and `R17_GO.txt` all exist and are non-empty.
- **Staged Three-Party Tolerance Ladder (Checked Strictly in Order by Master Controller):**

| Verification Stage | Artifacts Compared (Main Line vs. `A-T3` / Controller) | Locked Acceptance Tolerance |
|---|---|---|
| **Stage 1: Control-Cell Barcodes** | `control_barcodes_r17.json` & `.npz` `control_barcodes` (exact Hamilton largest-remainder stratified sample, `seed = 20261701`) | **`100%` exact ordered string equality** |
| **Stage 2: `SE-600M` `X_state` Embeddings** | `round17/results/embeddings/x_state_<line>.npz` (`float32`, official `FilteredGenesCounts` + `VCIDatasetSentenceCollator`, per-pool `seed = 20261701`, checked **before** transition passes) | **$\max_{i, d} \left|X_{\mathrm{state}}^{\mathrm{Main}} - X_{\mathrm{state}}^{\mathrm{AT3}}\right| \le 10^{-5}$** |
| **Stage 3: Predicted $\Delta$ Arrays (`56` Indices + `L1`/`L2`/`G4`)** | Exported per-line, per-model `.npz` arrays (`preds_STATE_{SE,HVG}_b{00..55}.npz`, `preds_*_L1_ens56.npz`, `preds_*_G4_*.npz`, plus `'non-targeting'` pseudobulks for `L2_ntc_recon_match`) | **$\max_{p, g} \left|\hat{\Delta}^{\mathrm{Main}} - \hat{\Delta}^{\mathrm{AT3}}\right| \le 10^{-4}$** |
| **Stage 4: Controller Same-Array Recomputation** | Master Controller recomputing every `H-J1`, `H-J2`, `H-J3`, and `J5` statistic directly from exported `.npz` arrays (`seed = 20261701`, `n_boot = 4,000`, `n_perm = 10,000`) | **$\left|\mathrm{metric}^{\mathrm{Controller}} - \mathrm{metric}^{\mathrm{Party}}\right| \le 10^{-10}$** |
| **Stage 5: Cross-Party Metrics & Verdicts** | Main Line vs. `A-T3` across all 56 indices, $\widehat{\Delta}_{\mathrm{CF}, 0}$, $\hat{\tau}_{\mathrm{batch}}$, $p_{\mathrm{var}}$, `has_detectable_signal`, $f_{\mathrm{flip}}$, $\hat{b}_{\mathrm{NTC}}$ + `top1-top2` gap, `L1`/`L2`, `J5` cross-fitted `G4`, & 56-index Spearman + permutation $p$ | **$\left|\mathrm{metric}^{\mathrm{Main}} - \mathrm{metric}^{\mathrm{AT3}}\right| \le 10^{-5}$** **AND `100%` verdict parity** (`H-J1`, `H-J2`, `H-J3`, exact $\hat{b}_{\mathrm{NTC}}$ index match) |

---

## §6. Summary of Resolved Controller Review R1 Items (`J1`–`J7`) & Rulings 1–3

| Review Item / Ruling | Resolution Locked in `PREREG_J_DRAFT_v2.md` | Section(s) |
|---|---|---|
| **`J1` & Ruling 1** (Cohort = Track 1 + Track 2 -> $K_{\mathrm{conf}} \in \{2, 3\}$) | Locked `cd4t_stim48` + `gse295214_h1` (prospective sequencing: `PREREG_J` locked before ANY `STATE` inference on `gse295214_h1`) + `pfizer_telohaec_basal` (`Zenodo 18792213`, `imru_full.h5ad`, `CC BY 4.0`, `Untreated` primary, `Treated` secondary within-family replication, TeloHAEC lineage disclosure vs. local `telohaec`, pre-download reader & synthetic tests verified). Locked `k/K` per model, zero pooling, general claim requires `>= 2` independent confirmatory units. | `§2.1`, `§2.3`, `§3.1` |
| **`J2`** (Floor effect) | Added detectable-signal predicate (`CF-Oracle` OR `L1_ens56` `98.333%` CI lower bound `> 0`) required for both `H-J1` `BATCH_INVARIANT_WITHIN_DELTA` and `H-J2` `CALL_STABLE`; otherwise classifies as `INCONCLUSIVE_LOW_SIGNAL`. Always reports `max_b r(b)` and `L1_ens56` score with every `H-J1` and `H-J2` call. | `§3.2`, `§3.3` |
| **`J3`** (Power statement) | Added explicit caveat that R16 discovery 4-index ranges vary widely (`0.1563`, `0.0598`, `0.0539`, `0.0374`), so design guards against false positives while `J2` handles low-signal lines. Added closed-form parametric precision/MDE/power table as a function of $N_{\mathrm{eval}} \in \{50..600\}$ for `gse295214_h1` and `pfizer_telohaec_basal`, filled mechanically from counts-only QC before inference. | `§4.1`, `§4.3` |
| **`J4`** (Selection-free index set) | Explicitly locked all 56 trained indices (`0..55`), identified row `0` (`fixed_index_zero`) as `arc-state 0.11.1`'s unseen-batch fallback (`state/_cli/_tx/_infer.py:688`, verified at Gate 2), and locked that **no index is ever dropped after inference**. | `§0.1`, `§1.2`, `§3.2` |
| **`J5`** (Mechanism tie-in) | Added (a) R16-style `G4` 2×2 decomposition at cross-fitted best/worst indices $(\hat{b}^*_1, \hat{b}^\dagger_1)$ selected on $\mathcal{P}_1$ and evaluated strictly on $\mathcal{P}_2$ only; (b) pairwise Spearman rank correlation of the 56 per-index scores across lines with exact `10,000`-permutation null over index labels (`seed = 20261701`) distinguishing shared batch-embedding property vs. line-specific interaction. | `§3.5` |
| **`J6`** (Blinding & release) | Recorded explicitly that `round17/` contains R16 main-line numbers (`§0.2`, `§4.1`), stays uncommitted and hidden from `A-T3` until Controller verifies `A-T3`'s R16 pass 2, and `A-T3` joins R17 only after release from R16, re-implementing strictly from locked `PREREG_J`. | Header, `§5.2` |
| **`J7`** (Additional candidates) | Documented `A-T1`'s parallel model-free screening deliverable (`round17/status/R17_CANDIDATES.md`, `< 100 MB` HTTP range reads, no full downloads, no VCC 2026 participant-channel data), non-blocking for `PREREG_J` v2 lock. | Header, `§2.2` |
| **Ruling 2** (Micro-N & `wtc11_i3neuron`) | Run 4 micro-N lines (`telohaec`, `calu3`, `mac_ctrl`, `mac_lps`) and `a549_ko_cas9pos` as descriptive / subset-sensitivity only; **omit inference on `wtc11_i3neuron` (`N = 1`)**. | `§2.2`, `§5.1` |
| **Ruling 3** (`L2_ntc_recon_match`) | Locked Pearson primary, negative L2 distance & per-fold selection secondary, AND mandatory reporting of `criterion_top1 - criterion_top2` selection margin to expose near-degenerate selections. | `§3.4` |
