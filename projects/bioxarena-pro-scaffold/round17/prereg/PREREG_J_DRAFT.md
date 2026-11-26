# Round 17 Pre-Registration Protocol Draft (`PREREG_J_DRAFT.md`)

> **STATUS: DRAFT FOR MASTER CONTROLLER AUDIT & HUMAN PROJECT LEAD REVIEW — MODEL-FREE INVENTORY ONLY — ZERO INFERENCE EXECUTED.**
> **Hard-Guard & Blinding Attestation:**
> - **Zero inference executed:** Zero `STATE` (`ST-SE-Replogle`, `ST-HVG-Replogle`) or `SE-600M` forward passes have been run on any Round 17 candidate line (`0` real cells inferred; `0` MPS/GPU compute).
> - **Model-free inventory only:** All line counts and eligibility flags in [`round17/status/LINE_INVENTORY_J.json`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round17/status/LINE_INVENTORY_J.json) were generated strictly by [`round17/src/build_line_inventory_j.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round17/src/build_line_inventory_j.py) (zero `torch`/`state` imports, zero correlation/metric calculations on target expression deltas).
> - **Strict workspace & governance boundaries:** Zero modifications under `round16/`, `round16_addB/`, or Rounds 9–15 frozen files; zero reads of `round16/audit/` or `A-T3.txt`; zero downloads; zero package installs; zero git operations; zero priority claims (`"first"`).
> - **Discovery-to-Confirmatory Quarantine:** All 6 Round 16 lines (`nadig_hepg2`, `nadig_jurkat`, `a549_ko`, `dld1`, `cd4t_rest`, `feng_ipsc`) are strictly quarantined as **Discovery Context Only** and are NEVER included in Round 17 confirmatory testing.

- **Project:** AI-Biology-Lab Round 17 — Pre-Registered Confirmatory Audit of Batch-Token (`batch_encoder`) Index Dependence, `STATE - FCT` Verdict Instability, and Label-Free Index Selection in Zero-Shot `Arc STATE` (`ST-SE-Replogle` Primary; `ST-HVG-Replogle` Co-Evaluated Second Model)
- **Protocol Author:** Main Line (`Antigravity`)
- **Companion Deliverables:**
  - Prior-Art Review: [`round17/prereg/PRIOR_ART_J.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round17/prereg/PRIOR_ART_J.md)
  - Model-Free Line Inventory: [`round17/status/LINE_INVENTORY_J.json`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round17/status/LINE_INVENTORY_J.json)
  - Inventory Builder Script: [`round17/src/build_line_inventory_j.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round17/src/build_line_inventory_j.py)
- **Parent Protocols (Locked, Unmodified):** [`round16/prereg/PREREG_I.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md) (`v5`, `c1fc7c6`), [`round16/prereg/PREREG_I_ADDENDUM_A.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I_ADDENDUM_A.md) (`Items 1–11`).
- **Fresh Unshared Round 17 RNG Seeds:**
  - Primary control sampling, `SE-600M` per-pool encoding, perturbation half-splits (`s = 0..19` via `20261701 + s`), and joint perturbation bootstrapping (`n_boot = 4,000`): **`20261701`** (`np.random.default_rng(20261701)`).
  - Secondary within-line disjoint control resample (when applicable): **`20261702`**.

---

## §0. Scientific Rationale, Architectural Mechanism, & Round 16 Discovery Context

### §0.1 Architectural Mechanism Verified Directly From Shipped Weights & `arc-state 0.11.1` Source
Read-only inspection of `arc-state 0.11.1` source code and all 8 shipped zero-shot checkpoints (`round16/shared/models/ST-{SE,HVG}-Replogle/zeroshot/{hepg2,jurkat,k562,rpe1}`, verified `77/77` SHA-256 OK, zero inference) establishes four exact facts:

1. **Additive `batch_encoder` Dominates Input Energy:** Both `ST-SE-Replogle` and `ST-HVG-Replogle` instantiate `self.batch_encoder = nn.Embedding(56, 128)` (`batch_dim = 56`, `hidden_dim = 128`, `batch_col = gem_group`, `use_batch_token = false`, `batch_predictor = false`; `state/tx/models/state_transition.py:192–197`) and add the retrieved vector unnormalized into the sequence input (`seq_input = pert_embedding + control_cells + self.batch_encoder(batch_indices.long())`, `state_transition.py:422–441`). Across all 8 checkpoints, row L2 norms are `9.77–12.43` (median `11.09`) in `ST-SE` and `8.81–12.52` (median `10.88`) in `ST-HVG` (near PyTorch $\mathcal{N}(0, 1)^{56 \times 128}$ initialization norm $\sqrt{128} = 11.31$), compared to `||basal_encoder(X_ctrl)||_2 ≈ 0.38–1.20` and `||pert_encoder(pert)||_2 ≈ 0.74`.
2. **Nearly Orthogonal Rows Shared Across All 4 Zero-Shot Folds:** Within each fold, the 56 rows are nearly mutually orthogonal (off-diagonal row cosine `mean = 0.0014, sd = 0.0919` in `ST-SE`; `0.0006 ± 0.0896` in `ST-HVG`), so averaging rows in embedding space (`mean_batch_emb`) collapses vector norm by $\approx \sqrt{56} \approx 7.5\times$ (`||e_mean||_2 = 1.55` in `ST-SE`, `1.47` in `ST-HVG`). Across the 4 zero-shot folds (`hepg2`, `jurkat`, `k562`, `rpe1`), identical row indices remain tightly aligned (`mean_cos = 0.9995–0.9997`, `min_cos >= 0.9991` in both `ST-SE` and `ST-HVG`), confirming that applying a single index $b \in \{0..55\}$ consistently across the 4-fold ensemble (`zeroshot/ensemble4`) is geometrically coherent.
3. **Silent Unseen-Batch Fallback to `batch_idx = 0` (`_infer.py:688`):** At zero-shot inference on unseen cell lines, official `arc-state 0.11.1` (`state/_cli/_tx/_infer.py:678–697`) looks up `.obs[batch_col]` in `batch_onehot_map.torch` (`{'1'->0, '10'->1, ..., '56'->51, '6'->52, ..., '9'->55}`) and silently sets **`idxs[i] = 0  # fallback to zero`** (`_infer.py:688`) for every unseen label (and crashes with `KeyError: 'batch'` at `state_transition.py:427` if no batch column is present).
4. **Unshipped Batch Provenance (`未能确定`):** Exhaustive read-only inspection of all 77 files under `round16/shared/models/` confirms that **neither per-`gem_group` training-cell embeddings nor a `gem_group`-to-training-cell-line lookup table are shipped**.

### §0.2 Round 16 Discovery Observations (Strictly Discovery Context Only — Excluded from Round 17 Confirmatory Family)
In Round 16 (`round16/results/summary_r16.json` and `round16/results/sensitivity/summary_r16.json`), secondary neutral-batch and 4-index (`{0, 1, 17, 35}`) diagnostics on the 6 Round 16 lines surfaced pronounced batch-token dependence that motivates Round 17's confirmatory study:

- On **`a549_ko`** (`N = 403`, `FCT r_orth_cal = 0.2134`): `STATE_SE` shifted across just 4 trained indices from `b=17: 0.0863` and `b=0 (fixed_index_zero): 0.1143` (`STATE_SE - FCT = -0.0991`, `FCT_SUPERIOR_TO_MODEL`) up to `b=35: 0.2426` (`4`-index range `0.1563`, `SD = 0.0702`), while `multi_batch_ens8` reached `0.2874` (`STATE_SE - FCT = +0.0740`, `97.5% CI [+0.0412, +0.1090]`, flipping the H-I2 verdict to `MODEL_SUPERIOR_TO_FCT`).
- On **`nadig_hepg2`** (`N = 740`, `FCT = 0.2015`): `fixed_index_zero` gave `0.1053` (`FCT_SUPERIOR_TO_MODEL`), whereas `bypass_zero` reached `0.2341` (`STATE_SE - FCT = +0.0326 [+0.0139, +0.0508]`, flipping to `MODEL_SUPERIOR_TO_FCT`) and `multi_batch_ens8` / `mean_batch_emb` reached `0.2008` / `0.2022` (`INDISTINGUISHABLE`).
- On **`dld1`** (`N = 75`), **`cd4t_rest`** (`N = 245`), and **`feng_ipsc`** (`N = 88`): the 4-index range spanned `0.0598`, `0.0374`, and `0.0539`, and **different indices ranked highest on different cell lines** (`b=35` on `a549_ko`; `b=0` on `dld1` & `cd4t_rest`; `b=1` on `feng_ipsc`), while `G4` 2×2 `b0 vs b17` median `batch_share` across unseen mismatched controls reached `57.7%–94.4%`.

Because these patterns were discovered on the 6 Round 16 lines, **none of `{nadig_hepg2, nadig_jurkat, a549_ko, dld1, cd4t_rest, feng_ipsc}` may enter Round 17 confirmatory testing**. Round 17 pre-registers a full 56-index confirmatory audit on fresh, never-STATE-run lines.

### §0.3 Prior-Art Positioning (Zero Priority Claims)
Full primary-source verification is logged in [`round17/prereg/PRIOR_ART_J.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round17/prereg/PRIOR_ART_J.md). Round 17 explicitly credits **`CELLAUDIT`** (`arXiv:2609.27234`, frozen-checkpoint input replacement & factorial auditing), **`Zhang et al. (2026)`** (`arXiv:2611.05849` §C.5, `(control, batch / gem_group)` ridge probes reaching PDC `0.27–0.46` and internal `L0` direction swapping), **`Molina & Zhang (2026)`** (`10.64898/2026.07.24.740459`, four-component response decomposition on Replogle/Nadig), **`OCOO-T`** (`arXiv:2606.12838`), and **`VCBench`** (`arXiv:2604.27646`). **No priority claims (`"first"`) are made.**

---

## §1. Audited Models, Fixed Inference Specification, & Endpoint Parity

### §1.1 Explicit Inclusion of Both `ST-SE` (Primary) and `ST-HVG` (Co-Evaluated Second Model)
Both zero-shot `Arc STATE` families shipped under `round16/shared/models/` carry the exact same `batch_encoder = nn.Embedding(56, 128)` (`batch_dim = 56`) architecture and cross-fold geometry (`mean_cos = 0.9997`). Therefore, **`ST-HVG-Replogle` IS EXPLICITLY INCLUDED** as the co-evaluated second model under the **exact same `H-J1`, `H-J2`, `H-J3`, and secondary `G4` 2×2 rules**:

| Role | Model Key | Checkpoints (All 4/4 Unseen on All Fresh Lines) | Control Input & Decoding Lock (Identical to `PREREG_I §1.2 & §2.1`) | Batch-Encoder Architecture |
|---|---|---|---|---|
| **Primary Model** | **`STATE_SE`** | `ST-SE-Replogle/zeroshot/{hepg2,jurkat,k562,rpe1}` equal-weight 4-fold mean (`zeroshot/ensemble4`) | Raw integer UMI NTC counts $\to$ `SE-600M` (`MPS float32`, official `FilteredGenesCounts` + `VCIDatasetSentenceCollator(training=False)`, per-pool seed `20261701`, `bs=16`, `shuffle=False`, `num_workers=0`) $\to$ `X_state` `(n_ctrl, 2058)` $\to$ `ST-SE` (`padded=True`, `16×32` sets) $\to$ `clip(gene_decoder(preds), 0, 14)` | `nn.Embedding(56, 128)`, row norms `9.77–12.43` (med `11.09`), cross-fold `cos >= 0.9991` |
| **Co-Evaluated Second Model** | **`STATE_HVG`** | `ST-HVG-Replogle/zeroshot/{hepg2,jurkat,k562,rpe1}` equal-weight 4-fold mean (`zeroshot/ensemble4`) | Raw integer UMI NTC counts $\to$ official `arc-state` `official_median_log1p` (`sc.pp.normalize_total(target_sum=None) -> log1p`, `HGNCBridge` onto ordered `2,000` `G_model` columns) $\to$ `ST-HVG` (`padded=True`, `16×32` sets) $\to$ `clip(preds, 0, 14)` directly (**NEVER** untrained `gene_decoder`) | `nn.Embedding(56, 128)`, row norms `8.81–12.52` (med `10.88`), cross-fold `cos >= 0.9993` |
| **Reference Baseline** | **`FCT`** | Frozen `round10/src/fct_pipeline.py` on `r9_6` (`6` source lines) | Exact `PREREG_I §1.2` `FCT` evaluated on identical `(T_eval(t), G_eval(t))` | N/A (closed-form target baseline mean scaling) |

*Scientific rationale for including `STATE_HVG`:* Because `STATE_HVG` bypasses `SE-600M` and `LatentToGeneDecoder` completely while sharing `StateTransitionPerturbationModel`'s Llama backbone and `batch_encoder = nn.Embedding(56, 128)`, running identical `H-J1–H-J3` tests on `STATE_SE` and `STATE_HVG` isolates whether batch-index instability is intrinsic to the transition backbone + `batch_encoder` rather than an artifact of `SE-600M` embeddings or the latent-to-gene decoder.

### §1.2 Control-Cell Sampling & Per-Index Prediction Formula
- **Control Sampling (`PREREG_I §2.2` with Fresh Seed `20261701`):** For each line $t$, sample `n_ctrl = 512` QC-passing NTC cells without replacement proportionally stratified across experimental batches using exact integer Hamilton largest-remainder quotas (`seed = 20261701`) and partition into `16 × 32` sets (`padded=True`). When `128 <= N_NTC < 512` (applies to `wtc11_i3neuron`, `423 -> 416 = 13 × 32`), truncate to `(N_NTC // 32) * 32`. All 56 batch indices and all perturbations (including `'non-targeting'`) share the **exact same sampled control cells and `16 × 32` ordering** (`control_barcodes_r17.json`).
- **Per-Index Prediction Matrix ($\hat{\Delta}^{(m, b)}_t$):** For each model $m \in \{\texttt{STATE\_SE}, \texttt{STATE\_HVG}\}$ and each trained batch index $b \in \{0, 1, \dots, 55\}$ applied identically to all `n_ctrl` cells (`batch_indices = torch.full((n_ctrl,), b, dtype=torch.long)`), compute per-fold decoded pseudobulk expression in clipped `[0, 14]` space, average across all 4 zero-shot folds (`zeroshot/{hepg2,jurkat,k562,rpe1}`), and subtract the 4-fold average predicted `'non-targeting'` expression under that **same** batch index $b$:
  $$\overline{\hat{Y}}^{(m, b)}_t(p) = \frac{1}{4}\sum_{f=1}^{4} \frac{1}{n_{\mathrm{ctrl}}}\sum_{i=1}^{n_{\mathrm{ctrl}}} \operatorname{clip}_{[0, 14]}\!\big(\hat{Y}^{(m, f)}(x_i, p, b)\big), \qquad \hat{\Delta}^{(m, b)}_t(p) = \overline{\hat{Y}}^{(m, b)}_t(p) - \overline{\hat{Y}}^{(m, b)}_t(\text{"non-targeting"}).$$
  *(Per-fold 56-index summaries are also saved as secondary diagnostics to verify that batch-index effects are not driven by a single fold.)*
- **Primary Endpoint (`r_orth_cal`) & Secondary Endpoint (`r_dorth_cal`):** Identical to `PREREG_I §1.1 & §4.1` (`round16/src/metrics_r16.py`): unselected `all_nontarget` ($G_{\mathrm{eval}}(t) \setminus \{g_p\}$), orthogonalized to per-perturbation `Transfer_Mean` $\bar{\Delta}^{\mathrm{src}}_p$ (fitted on `r9_6`), Spearman–Brown disattenuated by orthogonal split-half correlation explicitly keyed **`rho_ab_orth`**, clipped to `[-1, 1]`.

---

## §2. Model-Free Candidate Line Inventory, Family Structure, & Confirmatory Governance

### §2.1 Exhaustive Model-Free Inventory of All 9 Local Never-STATE-Run Lines
Computed strictly from local raw/standardized files and `var_dims.pkl` vocabulary metadata by [`round17/src/build_line_inventory_j.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round17/src/build_line_inventory_j.py) $\to$ [`round17/status/LINE_INVENTORY_J.json`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round17/status/LINE_INVENTORY_J.json) (zero `torch`/`state` imports, zero model outputs, zero target-expression correlations):

| Round | Line Key | Modality / Assay | Study / Release (`SE-600M` `2025-06` Chronology) | Raw Int UMI Verified (Size) | Prior `STATE` Exposure | `NTC` Total (`Half A / B`; `samp`) | `\|G_eval^R9\|` / **`\|G_eval^STATE\|`** | **`N_eval` (`PREREG_I`)** (`half_valid`; `R9_P5`) | `PREREG_I` Eligible? | Family Structure & Independence vs. R16 Discovery (`{nadig_hepg2, nadig_jurkat, a549_ko, dld1, cd4t_rest, feng_ipsc}`) | Locked R17 Role |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `r11` | **`cd4t_stim48`** | `CRISPRi` (`10x Flex`) | `GSE314342`, `2026-02-06` (**Post-`SE-600M` isolated**) | **Yes** (`27.39 GB`, 48 `.h5`) | **`NONE`** | `166,871` (`78,168 / 88,703`; `512`) | `6,131` / **`1,836`** | **`401`** (`333`; `503`) | **PASS** | Disjoint donors (`D2/D3` vs. `D1/D4`), disjoint lanes (`R2L25–48`), & stimulated (48h) vs. resting relative to R16 `cd4t_rest`, **but shares primary CD4$^+$ T lineage, sgRNA library, & Marson lab (`GSE314342`)** | **CONFIRMATORY CANDIDATE (Mandatory Shared-Study Disclosure)** |
| `r12` | `a549_ko_cas9pos` | `CRISPRko` (`CROP-seq`, `CAS9>0`) | `GSE345058`, `2026-08-26` (Post-`SE-600M`) | **Yes** (`10.96 GB`, `.mtx.gz`) | **`NONE`** | `33,685` (`17,435 / 16,250`; `512`) | `6,121` / **`1,837`** | **`403`** (`403`; `445`) | **PASS (counts only)** | **NON-INDEPENDENT SUBSET:** strict `CAS9>0` cell subset (`73.3%` of cells; `74.3%` of NTCs) of R16 discovery line `a549_ko` from identical raw `.mtx.gz` | **EXCLUDED FROM CONFIRMATORY** *(R16 subset sensitivity only)* |
| `r14` | `a549_myllia` | `CRISPRko` (`clonal Cas9`) | `GSE337804`, `2026` (**Post-`SE-600M` isolated**) | **Yes** (`1.17 GB`, 2 `.h5ad`) | **`NONE`** | `2,093` (`1,033 / 1,060`; `512`) | `6,126` / **`1,836`** | **`20`** (`20`; `32`) | FAIL (`< 50`) | Independent lab/chemistry/dataset (`Myllia`) from R16 `a549_ko`, same canonical cell line (`A549`, CRISPRko) | **SECONDARY DESCRIPTIVE** *(underpowered `N = 20`)* |
| `r12` | `hek293_xucao` | `CRISPRi` (`PerturbSci-Kinetics`) | `GSE218566`, `2023` (Predates `SE-600M`) | **Yes** (`0.34 GB`, 1 `.h5ad`) | **`NONE`** | `2,757` (`1,378 / 1,379`; `512`) | `6,126` / **`1,838`** | **`16`** (`16`; `16`) | FAIL (`< 50`) | Independent lineage (`HEK293`) & study from all R16 lines; `seed42` 50/50 split-half fallback | **SECONDARY DESCRIPTIVE** *(underpowered `N = 16`)* |
| `r11` | `telohaec` | `CRISPRi` (`10x 3'`) | `IGVFDS7943IIWZ` / `GSE210681`, `2023–24` (Predates `SE-600M`) | **Yes** (`1.47 GB`, 1 `.h5ad`) | **`NONE`** | `2,215` (`1,142 / 1,073`; `512`) | `4,869` / **`1,475`** | **`9`** (`9`; `9`) | FAIL (`< 50`) | Independent aortic endothelial lineage & study | **EXPLORATORY MICRO-N** (`N = 9`, descriptive only) |
| `r10` | `calu3` | `CRISPRi` (`10x 3'` bystander) | `GSE208240`, `2023` (Predates `SE-600M`) | **Yes** (`0.74 GB`, 1 `.h5ad`) | **`NONE`** | `3,094` (`1,486 / 1,608`; `512`) | `6,134` / **`1,840`** | **`8`** (`8`; `8`) | FAIL (`< 50`) | Independent lung epithelial (`Calu-3`) lineage & study | **EXPLORATORY MICRO-N** (`N = 8`, descriptive only) |
| `r13` | `mac_lps` | `CRISPRko` (`10x Flex`) | `GSE327644`, `2026-08` (**Post-`SE-600M` isolated**) | **Yes** (`3.57 GB`, 4 `.h5`) | **`NONE`** | `27,205` (`14,605 / 12,600`; `512`) | `6,132` / **`1,837`** | **`4`** (`4`; `5`) | FAIL (`< 50`) | Independent primary macrophage lineage vs. R16; **PAIRED with `mac_ctrl`** (same cells/donor/library under LPS stimulation) | **EXPLORATORY MICRO-N** (`N = 4`, paired with `mac_ctrl`) |
| `r13` | `mac_ctrl` | `CRISPRko` (`10x Flex`) | `GSE327644`, `2026-08` (**Post-`SE-600M` isolated**) | **Yes** (`3.01 GB`, 4 `.h5`) | **`NONE`** | `24,679` (`13,354 / 11,325`; `512`) | `6,132` / **`1,838`** | **`3`** (`3`; `3`) | FAIL (`< 50`) | Independent primary macrophage lineage vs. R16; **PAIRED with `mac_lps`** (same cells/donor/library unstimulated) | **EXPLORATORY MICRO-N** (`N = 3`, paired with `mac_lps`) |
| `r10` | `wtc11_i3neuron` | `CRISPRi` (`10x 3'`) | `GSE152988`, `2021` (Predates `SE-600M`) | **Yes** (`0.29 GB`, 1 `.h5ad`) | **`NONE`** | `423` (`194 / 229`; **`416`** `13×32`) | `5,995` / **`1,746`** | **`1`** (`1`; `2`) | FAIL (`< 50`) | Independent neuronal lineage (`WTC11 i3Neuron`); truncated control sample (`416 < 512`) | **EXCLUDED FROM INFERENCE** (`N = 1`: degenerate bootstrap & half-split; descriptive point check only if run) |

### §2.2 Critical Confirmatory Governance Finding & Decision Point for Master Controller / Human Project Lead
The model-free inventory reveals a decisive structural constraint across the 9 local never-STATE-run lines that MUST be settled explicitly before locking `PREREG_J`:

1. **Local Qualifying Census:**
   - **Zero (`0/9`)** local never-STATE-run lines are simultaneously (i) completely independent in study/lineage from all 6 R16 discovery lines **and** (ii) satisfy `PREREG_I`'s `N_eval >= 50` floor.
   - **Exactly ONE (`1/9`)** local non-subset line satisfies all `PREREG_I` eligibility gates: **`cd4t_stim48`** (`N_eval = 401`, `|G_eval^STATE| = 1,836`, `NTC = 166,871`, strictly post-dates `SE-600M`, zero prior STATE exposure, disjoint donors `D2/D3` & 48h stimulated biological state vs. R16 discovery `cd4t_rest`, with mandatory shared-study/lineage/library disclosure).
   - **`a549_ko_cas9pos`** (`N_eval = 403`) is a strict `CAS9>0` cell subset (`73.3%`) of R16 discovery `a549_ko` and is **strictly disqualified from confirmatory claims** (retained only as a subset-consistency check).
   - Among the remaining 7 truly non-R16-study local lines, the two largest are **`a549_myllia`** (`N = 20`, independent lab/dataset `GSE337804`, same `A549` line) and **`hek293_xucao`** (`N = 16`, independent `HEK293` lineage & study), followed by 4 micro-N lines (`telohaec = 9`, `calu3 = 8`, `mac_lps = 4`, `mac_ctrl = 3`) and degenerate `wtc11_i3neuron` (`N = 1`).
2. **Pre-Lock Cohort Architecture (Two Clean Options for Master Controller & Human Project Lead):**
   - **Track 1 (Zero-Download Local Execution, Ready Immediately Post-Lock):**
     - **Confirmatory Primary ($K_{\mathrm{conf}} = 1$):** **`cd4t_stim48`** (`N = 401`, CRISPRi), explicitly labeled **"Single-Line Confirmatory Test ($K_{\mathrm{conf}} = 1$) with Shared-Study/Lineage Disclosure (`GSE314342` vs. R16 `cd4t_rest`; disjoint donors `D2/D3` & stimulated 48h state)"** — mirroring Round 12's transparent `CONFIRMATORY_SINGLE_LINE_K1` reporting discipline (`PREREG_E`).
     - **Pre-Registered Secondary Fresh Lines (Non-Confirmatory, Direction & Landscape Replication):** **`a549_myllia`** (`N = 20`, CRISPRko) and **`hek293_xucao`** (`N = 16`, CRISPRi), plus micro-N descriptive table (`telohaec` `N = 9`, `calu3` `N = 8`, `mac_ctrl` `N = 3`, `mac_lps` `N = 4`, explicitly noting the `mac_ctrl`/`mac_lps` pair is not mutually independent; `wtc11_i3neuron` `N = 1` excluded from inferential statistics).
     - **Non-Independent R16 Subset Consistency Check (Strictly Separated):** `a549_ko_cas9pos` (`N = 403`), reported strictly in a separate sensitivity table to test whether removing `CAS9 == 0` cells alters R16's `a549_ko` batch-index pattern — NEVER counted toward $K_{\mathrm{conf}}$.
   - **Track 2 (Optional Prospective Cohort Expansion — Requires Separate Human Project Lead Download Approval):**
     - If a strictly study- and lineage-independent $K_{\mathrm{conf}} \ge 2$ confirmatory cohort is required, `PREREG_J` applies **verbatim without modifying any hypothesis, estimator, or threshold** to any separately approved prospective dataset meeting `N_eval >= 50, |G_eval^R9| >= 3,000, NTC_half >= 15` (e.g., `GSE295214` `H1 hESC` already specified in `PREREG_I §2.4`, or future public VCC 2026 lines once released on GCS; zero downloads performed in this draft stage).
3. **Quarantine of Round 16 Discovery Lines:** `{nadig_hepg2, nadig_jurkat, a549_ko, dld1, cd4t_rest, feng_ipsc}` are strictly excluded from all Round 17 confirmatory claims and $K_{\mathrm{conf}}$ counts. (If a full 56-index descriptive retrospective sweep is later added on R16 discovery lines for visualization only, it MUST be quarantined in an explicit `"Retrospective Discovery Context (Non-Confirmatory)"` appendix.)

---

## §3. Pre-Registered Hypotheses, Winner's-Curse Controls, Decision Rules, & Multiplicity Control

Every confirmatory line $t \in \mathcal{C}_{\mathrm{conf}}$ (`cd4t_stim48` under Track 1; plus any approved prospective lines under Track 2) is evaluated on both primary **`STATE_SE`** and co-evaluated second model **`STATE_HVG`** under the identical rules below. Secondary fresh lines (`a549_myllia`, `hek293_xucao`, micro-N lines) are scored identically but reported strictly as non-confirmatory descriptive results.

### §3.1 Multiplicity Control & Joint Bootstrap Architecture
- **Shared Joint Perturbation Bootstrap:** On each line $t$, generate `n_boot = 4,000` perturbation-resampling index draws (`seed = 20261701`, drawing $N_t$ perturbations with replacement per draw) shared identically across all 56 batch indices $b \in \{0..55\}$, `FCT`, `L1_ens56`, `L2_ntc_recon_match`, and secondary comparators, producing the exact `4000 × 56` bootstrap matrix $R^{\mathrm{boot}}_{t,m}$. Validity requires `valid_frac >= 0.95` (proportion of bootstrap draws yielding finite `rho_ab_orth > 0` and finite `r_orth_cal`) and finite `rho_ab_orth > 0` on the observed sample.
- **Family-Wise Error Control ($\mathrm{FWER} \le 0.05$):**
  - Within each confirmatory line $t$ and model $m$, there are **3 co-primary inferential hypothesis tests**:
    1. **`H-J1`**: Batch-index landscape heterogeneity & winner's-curse-free cross-fitted spread ($\widehat{\Delta}_{\mathrm{CF}, 0}$ and variance-component test $p_{\mathrm{var}}$);
    2. **`H-J3a`**: Label-free uniform 56-index ensemble `L1_ens56` vs. `fixed_index_zero` (`b = 0`);
    3. **`H-J3b`**: Label-free NTC control-reconstruction index `L2_ntc_recon_match` vs. `fixed_index_zero` (`b = 0`).
  - *(Note: `H-J2` is an audit quantification of verdict-call instability using the exact unchanged `PREREG_I §3.2` `97.5%` CI rule at each index $b$, paired with a Bonferroni-56 hardened flip check; it is NOT an additional directional hypothesis test competing for $\alpha$.)*
  - **Primary Family (`STATE_SE`, size $M_1 = 3 K_{\mathrm{conf}}$):** Bonferroni per-test level within a single line is $\alpha_{\mathrm{test}} = 0.05 / 3 = 0.016667$ (two-sided **`98.333%` percentile bootstrap CI**, quantiles `[0.8333%, 99.1667%]`, critical $z = 2.3940$). When $K_{\mathrm{conf}} > 1$, apply **Holm–Bonferroni step-down FWER control** across all $3 K_{\mathrm{conf}}$ tests (equivalently verifying that tests passing Holm also pass the corresponding Bonferroni-adjusted percentile CI at $\alpha / (3 K_{\mathrm{conf}} - j + 1)$).
  - **Co-Evaluated Second Model Family (`STATE_HVG`, size $M_2 = 3 K_{\mathrm{conf}}$):** Identical independent FWER $\le 0.05$ control ($\alpha_{\mathrm{test}} = 0.05 / 3$ within line, Holm–Bonferroni across $3 K_{\mathrm{conf}}$), PLUS a mandatory **pooled-model sensitivity** controlling FWER across all $6 K_{\mathrm{conf}}$ tests (`STATE_SE + STATE_HVG` jointly, $\alpha_{\mathrm{pool}} = 0.05 / 6 = 0.008333$ when $K_{\mathrm{conf}} = 1$, two-sided `99.167%` CI).

---

### §3.2 Hypothesis `H-J1` — Batch-Index Landscape & Winner's-Curse-Free Spread Across All 56 Indices (`0..55`)

**Scientific Question:** Holding matched target control cells $X_{\mathrm{ctrl}}$ and frozen weights fixed, how much does zero-shot `r_orth_cal` vary across all 56 trained batch indices $b \in \{0, 1, \dots, 55\}$, and does the true between-index spread exceed both perturbation-sampling noise and Round 16's practical-effect threshold $\delta = 0.04$?

1. **Full Descriptive Landscape Outputs:** Report all 56 point estimates $\{r_{t,m}(b)\}_{b=0}^{55}$ and their individual bootstrap SEs, empirical `min`, `q10`, `q25`, `median`, `mean`, `q75`, `q90`, `max`, `IQR`, observed between-index SD $s_{\mathrm{obs}} = \sqrt{\frac{1}{55}\sum_{b=0}^{55}(r(b)-\bar{r})^2}$, raw range $R_{\mathrm{obs}} = \max_b r(b) - \min_b r(b)$, and exact rank & percentile of official fallback `b = 0` (`fixed_index_zero`).
2. **Winner's-Curse-Free Co-Estimator 1 — Perturbation Half-Split Cross-Fitted Spread ($\widehat{\Delta}_{\mathrm{CF}, 0}$):**
   - **Why required:** Taking $\max_{b} \hat{r}(b) - \min_{b} \hat{r}(b)$ on the same perturbations inflates spread via selection bias (winner's curse), even under a sharp null $\theta_0 = \dots = \theta_{55}$.
   - **Locked Primary Split (`s = 0`, deterministic seed `20261701`):** Randomly permute the $N_t$ evaluated perturbations (sorted canonically by target symbol before permutation for cross-process determinism) using `np.random.default_rng(20261701)` and split into two disjoint halves $\mathcal{P}_1$ (size $\lfloor N_t / 2 \rfloor$) and $\mathcal{P}_2$ (size $\lceil N_t / 2 \rceil$), requiring `half_valid` counts to remain balanced (stratified by `half_valid` status when any `half_valid == False` perturbations exist, as in `cd4t_stim48` where `333/401` are `half_valid`, using exact deterministic alternating allocation after within-stratum permutation).
   - **Cross-Fitting Step:**
     - Select best and worst indices strictly on $\mathcal{P}_1$: $\hat{b}^*_{1} = \arg\max_{b \in 0..55} r_{\mathcal{P}_1}(b)$, $\hat{b}^{\dagger}_{1} = \arg\min_{b \in 0..55} r_{\mathcal{P}_1}(b)$ (lowest-index tie-break). Evaluate their paired contrast strictly out-of-sample on disjoint $\mathcal{P}_2$: $\Delta_{1\to 2} = r_{\mathcal{P}_2}(\hat{b}^*_{1}) - r_{\mathcal{P}_2}(\hat{b}^{\dagger}_{1})$.
     - Symmetrically select $(\hat{b}^*_{2}, \hat{b}^{\dagger}_{2})$ strictly on $\mathcal{P}_2$ and evaluate strictly out-of-sample on $\mathcal{P}_1$: $\Delta_{2\to 1} = r_{\mathcal{P}_1}(\hat{b}^*_{2}) - r_{\mathcal{P}_1}(\hat{b}^{\dagger}_{2})$.
     - Primary winner's-curse-free spread estimator:
       $$\widehat{\Delta}_{\mathrm{CF}, 0}(t, m) = \frac{1}{2}\!\left(\Delta_{1\to 2} + \Delta_{2\to 1}\right).$$
     - Under the sharp null $\theta_0 = \dots = \theta_{55}$ with independent perturbation errors across $\mathcal{P}_1 \perp \mathcal{P}_2$, $\mathbb{E}[\Delta_{1\to 2} \mid \mathcal{P}_1] = 0$ and $\mathbb{E}[\Delta_{2\to 1} \mid \mathcal{P}_2] = 0$, so $\mathbb{E}[\widehat{\Delta}_{\mathrm{CF}, 0}] = 0$ (zero winner's curse!).
   - **Exact Within-Half Bootstrap Inference (`n_boot = 4,000`):** Conditional on $(\hat{b}^*_{1}, \hat{b}^{\dagger}_{1})$ selected on $\mathcal{P}_1$, draw `4,000` perturbation bootstraps strictly inside $\mathcal{P}_2$ (`seed = 20261701`) to obtain $\Delta_{1\to 2}^{(u)}$; independently draw `4,000` perturbation bootstraps strictly inside $\mathcal{P}_1$ (`seed = 20261701 + 10000`) for fixed $(\hat{b}^*_{2}, \hat{b}^{\dagger}_{2})$ to obtain $\Delta_{2\to 1}^{(u)}$; form the exact paired half-averaged bootstrap distribution $\widehat{\Delta}_{\mathrm{CF}, 0}^{(u)} = \frac{1}{2}(\Delta_{1\to 2}^{(u)} + \Delta_{2\to 1}^{(u)})$ and extract its two-sided `(1 - \alpha_test)` (`98.333%`) and `95.0%` percentile CIs and empirical two-sided $p$-value.
   - **Multi-Split Stability Check (`R_split = 20`, Secondary Diagnostic):** Repeat the cross-fitting point estimate across 20 independent deterministic splits (`seed = 20261701 + s`, `s = 0..19`) and report mean, median, min, max, and IQR of $\widehat{\Delta}_{\mathrm{CF}, s}$ to verify `s = 0` is not a split outlier (only locked `s = 0` enters primary hypothesis testing).
3. **Winner's-Curse-Free Co-Estimator 2 — Noise-Debiased Between-Index SD ($\hat{\tau}_{\mathrm{batch}}$) & Centered-Bootstrap Homogeneity Test ($p_{\mathrm{var}}$):**
   - From the `4000 × 56` joint bootstrap matrix $R^{\mathrm{boot}}_{t,m}$ (evaluated on all $N_t$ perturbations), define the centering matrix $C_{56} = I_{56} - \frac{1}{56}\mathbf{1}\mathbf{1}^\top$, observed between-index variance $S_{\mathrm{obs}}^2 = \frac{1}{55}\|C_{56} \hat{r}\|_2^2$, and exact covariance-aware perturbation-sampling noise floor:
     $$V_{\mathrm{noise}} = \frac{1}{55}\sum_{b=0}^{55} \operatorname{Var}_{\mathrm{boot}}\!\big(R^{\mathrm{boot}}_{u, b} - \bar{R}^{\mathrm{boot}}_{u, \cdot}\big), \qquad \hat{\tau}_{\mathrm{batch}}(t, m) = \sqrt{\max\!\big(0,\; S_{\mathrm{obs}}^2 - V_{\mathrm{noise}}\big)}.$$
   - **Global Null Test (`H0: tau_batch = 0`):** Form centered null draws $Z^{(u)} = R^{\mathrm{boot}}_{u, \cdot} - \bar{R}^{\mathrm{boot}}_{\cdot}$ (which preserve the exact $56 \times 56$ perturbation-resampling covariance $\hat{\Sigma}_{\mathrm{noise}}$ while imposing true between-index variance `0`), compute $S_{\mathrm{null}}^2(u) = \frac{1}{55}\|C_{56} Z^{(u)}\|_2^2$, and define exact one-sided homogeneity $p$-value:
     $$p_{\mathrm{var}}(t, m) = \frac{1 + \sum_{u=1}^{4000} \mathbb{I}\!\big(S_{\mathrm{null}}^2(u) \ge S_{\mathrm{obs}}^2\big)}{4001}.$$
   - *(Note: Because $4,000$ draws give minimum resolution $p_{\mathrm{var}} \ge 1/4001 \approx 2.5 \times 10^{-4} \ll 0.01667$, this cleanly resolves Bonferroni and Holm–Bonferroni thresholds.)*
4. **Locked 4-Way Decision Rule for `H-J1` (Anchored to R16 Practical Threshold $\delta = 0.04$):**
   Let $\alpha_{\mathrm{adj}}$ denote the locked multiplicity-adjusted level (`0.01667` when $K_{\mathrm{conf}} = 1$, Holm–Bonferroni when $K_{\mathrm{conf}} > 1$):
   1. **`INCONCLUSIVE` (Validity Fail):** `valid_frac < 0.95` on any of the 56 indices or on either cross-fit half, OR `rho_ab_orth <= 0` or non-finite CIs;
   2. **`BATCH_DEPENDENCE_PRACTICALLY_MEANINGFUL` (H-J1 Confirmed):** ALL of:
      - $p_{\mathrm{var}} < \alpha_{\mathrm{adj}}$ **AND** $\mathrm{CI}_{1-\alpha_{\mathrm{adj}}}^{\mathrm{lo}}(\widehat{\Delta}_{\mathrm{CF}, 0}) > 0.0$ *(statistically significant non-zero batch-index variation, winner's-curse-free)*;
      - $\widehat{\Delta}_{\mathrm{CF}, 0} > 0.04$ **AND** $2\,\hat{\tau}_{\mathrm{batch}} > 0.04$ (i.e., $\hat{\tau}_{\mathrm{batch}} > 0.02$) *(practical spread exceeds Round 16's $\delta = 0.04$ meaningful-effect threshold)*;
   3. **`BATCH_DEPENDENCE_DETECTED_SUB_DELTA`:** $p_{\mathrm{var}} < \alpha_{\mathrm{adj}}$ AND $\mathrm{CI}_{1-\alpha_{\mathrm{adj}}}^{\mathrm{lo}}(\widehat{\Delta}_{\mathrm{CF}, 0}) > 0.0$, yet $\widehat{\Delta}_{\mathrm{CF}, 0} \le 0.04$ OR $\hat{\tau}_{\mathrm{batch}} \le 0.02$ *(mandatory wording: detected above noise but smaller than practical threshold $\delta = 0.04$)*;
   4. **`BATCH_INVARIANT_WITHIN_DELTA`:** NOT statistically detected ($p_{\mathrm{var}} \ge \alpha_{\mathrm{adj}}$ OR $\mathrm{CI}_{1-\alpha_{\mathrm{adj}}}^{\mathrm{lo}}(\widehat{\Delta}_{\mathrm{CF}, 0}) \le 0.0$), AND two-sided `95.0%` upper bound $\mathrm{CI}_{95\%}^{\mathrm{hi}}(\widehat{\Delta}_{\mathrm{CF}, 0}) < 0.04$ AND $\hat{\tau}_{\mathrm{batch}} \le 0.02$;
   5. **`INCONCLUSIVE` (Insufficient Precision):** All remaining cases (explicitly interpreted as insufficient statistical precision, NEVER as batch invariance proven).

---

### §3.3 Hypothesis `H-J2` — `STATE - FCT` Benchmark Call Instability Across All 56 Batch Indices

**Scientific Question:** Using the exact `PREREG_I §3.2` `H-I2` decision rule (`STATE - FCT`, two-sided `97.5%` joint perturbation bootstrap CI), what fraction of the 56 trained batch indices produce a **different** benchmark verdict from `arc-state 0.11.1`'s default `fixed_index_zero` (`b = 0`) fallback?

1. **Exact Per-Index `H-I2` Call (Identical to `PREREG_I §3.2`):** For every $b \in \{0, \dots, 55\}$, compute the paired contrast $\Delta^{\mathrm{vsFCT}}_{t,m}(b) = r_{\mathrm{orth\_cal}}(\mathrm{model}_{m, b}) - r_{\mathrm{orth\_cal}}(\mathrm{FCT})$ using the locked `4,000` joint perturbation bootstrap draws (`seed = 20261701`) and classify strictly as:
   - `MODEL_SUPERIOR_TO_FCT` iff `valid_frac >= 0.95` and `ci975_lo > 0.0`;
   - `FCT_SUPERIOR_TO_MODEL` iff `valid_frac >= 0.95` and `ci975_hi < 0.0`;
   - `INDISTINGUISHABLE` otherwise (including exact `0.0` ties or `valid_frac < 0.95`).
2. **Primary `H-J2` Quantity & Companion Multiplicity-Hardened Audit:**
   - Baseline call under official fallback: $c_0 = \mathrm{call}_{t,m}(0)$.
   - **Call-instability fraction:**
     $$f_{\mathrm{flip}}(t, m) = \frac{1}{56}\sum_{b=0}^{55} \mathbb{I}\!\big(\mathrm{call}_{t,m}(b) \neq c_0\big),$$
     reported with exact 3-category counts `(n_MODEL_SUPERIOR, n_INDISTINGUISHABLE, n_FCT_SUPERIOR)`, exact index lists for each category, and strict sign-reversal indicator `has_full_reversal = (n_MODEL_SUPERIOR >= 1) AND (n_FCT_SUPERIOR >= 1)`.
   - **Companion Bonferroni-56 Hardened Stability Check (Anti-Borderline-Noise Guard):** To distinguish genuine batch-driven verdict flips from borderline CIs grazing zero across 56 evaluations, also compute $f_{\mathrm{flip}}^{\mathrm{Bonf56}}(t, m)$ using two-sided `99.9554%` CIs ($\alpha = 0.025 / 56$ per index, via `n_boot = 20,000` nested/extended draws strictly for this companion check or normal-calibrated `z = 3.5117` using `boot_se` — **locked to exact `z = 3.5117` normal-calibrated check AND `n_boot = 20,000` empirical quantile check**, reported side-by-side!).
3. **Locked 4-Way Audit Classification for `H-J2`:**
   1. **`CALL_UNSTABLE_SIGN_REVERSAL`:** `has_full_reversal == True` (both `MODEL_SUPERIOR_TO_FCT` and `FCT_SUPERIOR_TO_MODEL` appear among `0..55`);
   2. **`CALL_UNSTABLE_SUBSTANTIAL`:** NOT full reversal, AND $f_{\mathrm{flip}}(t, m) \ge 5/56$ (`>= 8.93%`, at least 5 of 56 indices disagree with `fixed_index_zero`);
   3. **`CALL_UNSTABLE_MINOR`:** $1/56 \le f_{\mathrm{flip}}(t, m) < 5/56$ (`1–4` of 56 indices disagree with `b = 0`);
   4. **`CALL_STABLE`:** $f_{\mathrm{flip}}(t, m) = 0$ (all `56/56` indices match `fixed_index_zero`).
   *(Anti-Masking Rule from `PREREG_I §3.2`: Always report `FCT` `r_orth_cal` alongside `H-J2`, noting explicitly that on low-`FCT` lines such as T-cell screens where `FCT` is modest, transitions between `INDISTINGUISHABLE` and `FCT_SUPERIOR_TO_MODEL` require narrower shifts than on high-`FCT` lines.)*

---

### §3.4 Hypothesis `H-J3` — Pre-Specified Label-Free Batch-Index Rules vs. `fixed_index_zero` & Cross-Fitted Oracle

**Why *"Nearest Training-Batch Embedding"* Cannot Be Defined From Shipped Files:** CPU read-only inspection of all 77 files under `round16/shared/models/` verified that zero per-`gem_group` training-cell embeddings or `gem_group`-to-cell-line mappings are shipped (`§0.1`). Defining a training-batch centroid rule would require downloading and encoding raw Replogle-4 training datasets outside shipped model artifacts. Therefore, `H-J3` locks **two strictly self-contained, zero-perturbation-label rules computable solely from shipped checkpoints + the target line's own unlabeled NTC control cells**:

1. **Rule `L1_ens56` (Uniform 56-Index Decoded-Space Ensemble — Co-Primary `H-J3a`):**
   - Average decoded predictions uniformly over all 56 trained batch indices in clipped `[0, 14]` gene space before subtracting `'non-targeting'` (equivalently $\hat{\Delta}^{(m, \mathrm{L1})}_t(p) = \frac{1}{56}\sum_{b=0}^{55} \hat{\Delta}^{(m, b)}_t(p)$).
   - *Properties:* Zero perturbation labels, zero free parameters, zero external data; every single forward pass remains strictly inside the trained norm shell $\|e_b\|_2 \approx 11.1$ (unlike `mean_batch_emb` whose norm collapses to `1.55` or `bypass_zero` at norm `0`).
2. **Rule `L2_ntc_recon_match` (Unlabeled NTC Control-Profile Reconstruction Index — Co-Primary `H-J3b`):**
   - Using strictly the target line's sampled `n_ctrl` unlabeled unperturbed NTC cells (zero perturbed cells, zero perturbation labels), compute the observed unlabeled NTC mean `log1p`-normalized expression $\bar{Y}^{\mathrm{obs}}_{\mathrm{ctrl}}(t)$ on $G_{\mathrm{eval}}(t)$ (using `libnorm(1e4) -> log1p`, or `official_median_log1p` — note Pearson correlation across genes is invariant to additive shifts and positive scaling, and both normalizations yield monotonic profiles; **locked strictly to `official_median_log1p` mean of the sampled `n_ctrl` NTC cells on $G_{\mathrm{eval}}(t)$** matching `arc-state`'s training space).
   - For each $b \in \{0, \dots, 55\}$, compute the 4-fold mean predicted `'non-targeting'` pseudobulk $\overline{\hat{Y}}^{(m, b)}_t(\text{"non-targeting"})$ restricted to $G_{\mathrm{eval}}(t)$ and select:
     $$\hat{b}_{\mathrm{NTC}}(t, m) = \arg\max_{b \in \{0,\dots,55\}} \operatorname{Pearson}_{g \in G_{\mathrm{eval}}(t)}\!\Big(\overline{\hat{Y}}^{(m, b)}_t(\text{"non-targeting"}),\; \bar{Y}^{\mathrm{obs}}_{\mathrm{ctrl}}(t)\Big)$$
     (breaking exact ties by smallest integer index $b$).
   - *Properties:* Zero perturbation labels, zero external data; directly tests whether matching a training batch token's unperturbed baseline reconstruction to the target line's unperturbed transcriptome recovers a favorable batch index for perturbation responses. *(Per-fold $\hat{b}_{\mathrm{NTC}}^{(f)}$ selection is reported as a secondary diagnostic.)*
3. **Locked Paired Contrasts & Decision Rules for Each Rule $L \in \{\texttt{L1\_ens56}, \texttt{L2\_ntc\_recon\_match}\}$:**
   - **Primary Contrast vs. Official Fallback `fixed_index_zero` (`b = 0`):**
     $$\Delta_{L \text{ vs } 0}(t, m) = r_{\mathrm{orth\_cal}}(L) - r_{\mathrm{orth\_cal}}(b = 0),$$
     using the shared `4,000` joint perturbation bootstraps (`seed = 20261701`) and locked multiplicity-adjusted confidence level `1 - \alpha_adj` (`98.333%` when $K_{\mathrm{conf}} = 1$, Holm–Bonferroni when $K_{\mathrm{conf}} > 1$):
     - **`RULE_SUPERIOR_TO_ZERO`** iff `valid_frac >= 0.95` AND `ci_adj_lo > 0.0` (annotated `> 0.04` practical vs. `<= 0.04` sub-delta);
     - **`ZERO_SUPERIOR_TO_RULE`** iff `valid_frac >= 0.95` AND `ci_adj_hi < 0.0`;
     - **`INDISTINGUISHABLE`** otherwise (including exact `0.0` ties; note for `L2_ntc_recon_match`, if $\hat{b}_{\mathrm{NTC}} == 0$, report explicitly as `EXACT_TIE_SELECTED_ZERO` rather than `INDISTINGUISHABLE`!).
   - **Paired Gap vs. Winner's-Curse-Free Cross-Fitted Oracle (Mandatory Benchmark, Evaluated Strictly Out-of-Sample):**
     On locked split `s = 0` ($\mathcal{P}_1, \mathcal{P}_2$), define the cross-fitted oracle score $\widehat{r}_{\mathrm{CF\text{-}Oracle}, 0} = \frac{1}{2}\big(r_{\mathcal{P}_2}(\hat{b}^*_1) + r_{\mathcal{P}_1}(\hat{b}^*_2)\big)$ and the identically split-averaged rule score $\widehat{r}_{L, 0} = \frac{1}{2}\big(r_{\mathcal{P}_2}(L) + r_{\mathcal{P}_1}(L)\big)$, reporting their exact paired out-of-sample contrast $\Delta_{L \text{ vs CF-Oracle}} = \widehat{r}_{L, 0} - \widehat{r}_{\mathrm{CF\text{-}Oracle}, 0}$ with within-half joint bootstrap `98.333%` CI (so `L1_ens56` can even validly **outperform** the best single index if ensembling reduces single-index variance!).
   - **Benchmark Impact (`H-I2` Call Under Rule $L$):** Report `L vs FCT` under exact `PREREG_I §3.2` `97.5%` CI alongside `b = 0 vs FCT`.
4. **Secondary Non-Primary Neutral-Batch Comparators (Descriptive Continuity with `PREREG_I §4.4`):**
   - `bypass_zero` ($e = 0 \in \mathbb{R}^{128}$), `mean_batch_emb` ($e = \bar{e} \in \mathbb{R}^{128}$), `multi_batch_ens8` ($b \in \{0, 7, 14, 21, 28, 35, 42, 49\}$), and `medoid_idx` (weight-only index $\arg\min_b \sum_{b'}\|e_b - e_{b'}\|_2$ on the 4-fold mean Gram matrix).

---

### §3.5 Pre-Registered Secondary Mechanism Analysis — `G4` 2×2 (`X_ctrl × batch`) Factorial on Fresh Lines

To replicate and mechanistically decompose Round 16's `G4` finding on fresh lines without data-snooping index pairs, evaluate the exact `PREREG_I §4.4` / Addendum A §6 2×2 factorial for BOTH `STATE_SE` and `STATE_HVG` using **two strictly pre-fixed index pairs** (NEVER chosen from fresh-line outcomes):

1. **Pre-Fixed Pair 1 (Primary Continuity with `PREREG_I §4.4`):** `(B_M = 0, B_S = 17)` AND reversed symmetry check `(B_M = 17, B_S = 0)`.
2. **Pre-Fixed Pair 2 (Fixed Strictly From R16 Discovery `a549_ko`):** `(B_M = 0, B_S = 35)` AND reversed symmetry check `(B_M = 35, B_S = 0)` — explicitly testing whether the R16 `b0 vs b35` separation transfers to unseen lines or is cell-line-idiosyncratic.
3. **Locked Mismatched Control Pools (`X_ctrl = S`, Two-Stratum Reporting):**
   - **`seen_by_checkpoint` ($n = 2$):** `{nadig_hepg2, nadig_jurkat}` (carrying `PREREG_I` `M5.2` `3/4`-checkpoint overlap disclosure).
   - **`unseen_by_checkpoint` ($n = 3$–`4`, excluding same canonical family):** Local R16 control pools `{a549_ko, dld1, cd4t_rest, feng_ipsc}` excluding any pool sharing `canonical_family` with target line $t$ (specifically: on `cd4t_stim48`, foreign-lineage unseen pool = `{a549_ko, dld1, feng_ipsc}`, while `cd4t_rest` is reported separately as a **`same_lineage_resting_comparator`**; on `a549_myllia` and `a549_ko_cas9pos`, foreign-lineage unseen pool = `{dld1, cd4t_rest, feng_ipsc}`).
4. **Reported Quantities:** $\Delta_{\mathrm{batch}} = r(M, B_M) - r(M, B_S)$, $\Delta_{\mathrm{ctrl}} = r(M, B_M) - r(S, B_M)$, $\Delta_{\mathrm{both}} = r(M, B_M) - r(S, B_S)$, interaction $\Delta_{\mathrm{int}} = \Delta_{\mathrm{both}} - \Delta_{\mathrm{ctrl}} - \Delta_{\mathrm{batch}}$, and $|\Delta_{\mathrm{batch}}| / (|\Delta_{\mathrm{ctrl}}| + |\Delta_{\mathrm{batch}}|)$ under both `r_orth_cal` and `r_dorth_cal` (`4,000` joint bootstraps).

---

## §4. Power and Precision Statement (Grounded Directly in Round 16 Empirical Bootstrap SEs)

### §4.1 Exact Round 16 Planning Inputs Cited From `summary_r16.json`
All planning standard errors below are read directly from Round 16's executed outputs (`round16/results/summary_r16.json` and `round16/results/sensitivity/summary_r16.json`, `n_boot = 4,000`, exact JSON paths cited):

| R16 Discovery Line (`JSON source`) | $N_{\mathrm{eval}}$ | `rho_ab_orth` | `primary.STATE_SE_matched_orth.boot_se` (`b=0`) | Mean `boot_se` across `{0,1,17,35}` (`g4_single_index_spread`) | Paired `b0 - b17` `boot_se` (`g4_factorial_decomposition.b0_vs_b17*.r_orth_cal.delta_batch.boot_se`) | `primary.STATE_SE_vs_FCT.boot_se` | Observed 4-Index Range / SD (`g4_single_index_spread`) |
|---|---:|---:|---:|---:|---:|---:|---:|
| `nadig_hepg2` (`summary_r16.json`) | `740` | `0.5918` | `0.00332` | — | — *(native gem_group)* | `0.00585` | — |
| `nadig_jurkat` (`summary_r16.json`) | `873` | `0.4532` | `0.00261` | — | — *(native gem_group)* | `0.00713` | — |
| `a549_ko` (`summary_r16.json`) | `403` | `0.8271` | `0.00735` | `0.00740` | **`0.00748`** | `0.01172` | `0.1563` / `0.0702` |
| `dld1` (`summary_r16.json`) | `75` | `0.2931` | `0.00962` | `0.01203` | **`0.00889`** | `0.01778` | `0.0598` / `0.0254` |
| `cd4t_rest` (`summary_r16.json`) | `245` | `0.2062` | `0.01174` | `0.01153` | **`0.01118`** | `0.01223` | `0.0374` / `0.0161` |
| `feng_ipsc` (`sensitivity/summary_r16.json`) | `88` | `0.6177` | `0.01297` | `0.01048` | **`0.01466`** | `0.01795` | `0.0539` / `0.0247` |

### §4.2 Analytical Precision & Power for `H-J1` ($\widehat{\Delta}_{\mathrm{CF}, 0}$ and $\hat{\tau}_{\mathrm{batch}}$) and `H-J3`
For the cross-fitted spread $\widehat{\Delta}_{\mathrm{CF}, 0} = \frac{1}{2}(\Delta_{1\to 2} + \Delta_{2\to 1})$, each half contains $N/2$ perturbations so $\mathrm{SE}(\Delta_{1\to 2}) \approx \sqrt{2}\,\mathrm{SE}_{\mathrm{pair}}$, while averaging the two independent halves divides SE by $\sqrt{2}$, yielding exact algebraic cancellation: $\mathrm{SE}(\widehat{\Delta}_{\mathrm{CF}, 0}) \approx \mathrm{SE}_{\mathrm{pair}}$ (identical order of magnitude to a full-sample paired batch contrast `b_1 - b_2`, and identical to `H-J3`'s paired rule contrast $\Delta_{L \text{ vs } 0}$). For the global variance-component test (`df = 55`), approximating equicorrelated perturbation-sampling error $\sigma_e = \mathrm{SE}_{\mathrm{pair}} / \sqrt{2}$ yields $55 S_{\mathrm{obs}}^2 / \sigma_e^2 \sim (1 + \tau_{\mathrm{batch}}^2 / \sigma_e^2)\,\chi^2_{55}$ under $H_1$. Evaluating at our locked Bonferroni 3-test critical value $\alpha_{\mathrm{test}} = 0.05 / 3 = 0.01667$ (two-sided $z = 2.3940$):

| Target / Planning Scenario | Assumed $\mathrm{SE}_{\mathrm{pair}}$ (Basis in R16) | `98.33%` CI Half-Width | MDE (`80%` / `90%` Power) | Power at $\Delta = 0.04$ (`R16 δ`) | Power at $\Delta = 0.06$ | Variance Test Power ($\tau = 0.02$, i.e. $2\tau = 0.04$) | Variance Test Power ($\tau = 0.015$) | Interpretation & Power Assessment |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| **`cd4t_stim48` (`N = 401`) — scaled from `cd4t_rest`** | `0.00874` (`0.01118 × √(245/401)`) | **`±0.0209`** | `0.0283` / `0.0321` | **`98.5%`** | **`>99.9%`** | **`>99.9%`** | **`>99.9%`** | **Well-powered** for both $\delta = 0.04$ spread & `H-J3` gains |
| **`cd4t_stim48` (`N = 401`) — conservative unscaled `cd4t_rest`** | `0.01118` (unscaled `cd4t_rest` `b0-b17` SE) | **`±0.0268`** | `0.0362` / `0.0411` | **`88.2%`** | **`99.9%`** | **`>99.9%`** | **`>99.9%`** | **Powered (`>88%`)** even without `N=401` variance gain |
| High-N / High-$\rho^\perp_{AB}$ Prospective (`N ≈ 400`, `a549_ko`-like) | `0.00748` (`a549_ko` `b0-b17` SE) | **`±0.0179`** | `0.0242` / `0.0275` | **`99.8%`** | **`>99.9%`** | **`>99.9%`** | **`>99.9%`** | **Very high power** (`>99%`) |
| Moderate-N (`N ≈ 75–88`, `dld1`–`feng_ipsc` range) | `0.00889 – 0.01466` | `±0.0213 – ±0.0351` | `0.0288 – 0.0474` / `0.0327 – 0.0539` | `63.1% – 98.2%` | `95.5% – >99.9%` | `>99.9%` | `>99.9%` | Adequate for variance test & large spreads (`>= 0.06`); borderline at `0.04` if SE near `0.015` |
| **Secondary Fresh Local (`a549_myllia` `N=20`, `hek293_xucao` `N=16`)** | `~0.0300` (scaled $\sim 3\times$ for `N=16–20`) | **`±0.0718`** | `0.0971` / `0.1103` | **`14.4%`** | **`34.7%`** | `89.8%` | `54.9%` | **Severely underpowered for $\delta = 0.04$ cross-fitted spread**; strictly descriptive secondary only |

*Key Power Takeaway:* **`cd4t_stim48` (`N = 401`) achieves `88.2%–98.5%` power** at $\Delta = 0.04$ (and **`>99.9%` power** for the 56-index variance-component test at $\tau_{\mathrm{batch}} = 0.02$) under Bonferroni 3-test adjustment (`98.33%` CIs), whereas **`a549_myllia` (`N = 20`) and `hek293_xucao` (`N = 16`) have only `~14%` power** at $\Delta = 0.04$ — confirming quantitatively why they MUST remain descriptive secondary cohorts.

---

## §5. Compute Plan & Three-Party Independent Verification Plan (`A-T3` Post-R16)

### §5.1 Local Compute & Memory Budget (Apple Silicon Mac, `128 GB` Unified RAM, `MPS float32`)
All latencies below are anchored directly to Round 16's verified `MPS float32` measurements (`R16_GATE2_SCHEMA.md §6.2` & `round16/results/summary_r16.json`: `SE-600M = 331 s [5.5 min]` per `512`-cell pool, peak RSS `5.70 GB`; `ST-SE = 12.65 ms / pert` per fold on `16×32`; `ST-HVG = 10.21 ms / pert` per fold on `16×32`, peak RSS `2.35 GB`):

| Line Key (Role) | $N_{\mathrm{eval}}$ (+1 NTC) | `SE-600M` (`MPS fp32`) | `ST-SE` (`56 idx × 4 folds` + `G4`) | `ST-HVG` (`56 idx × 4 folds` + `G4`) | Total Compute |
|---|---:|---:|---:|---:|---:|
| **`cd4t_stim48`** (Confirmatory Candidate) | `401` (`402`) | `5.5 min` | `19.0 + 4.7 = 23.7 min` | `15.3 + 3.8 = 19.1 min` | **`48.3 min`** |
| `a549_myllia` (Secondary Fresh) | `20` (`21`) | `5.5 min` | `1.0 + 0.2 = 1.2 min` | `0.8 + 0.2 = 1.0 min` | `7.7 min` |
| `hek293_xucao` (Secondary Fresh) | `16` (`17`) | `5.5 min` | `0.8 + 0.2 = 1.0 min` | `0.6 + 0.2 = 0.8 min` | `7.3 min` |
| `telohaec`, `calu3`, `mac_ctrl`, `mac_lps` (4 Micro-N) | `9, 8, 3, 4` | `22.0 min` | `1.6 min` | `1.4 min` | `25.0 min` |
| *Subtotal: 7 Non-Subset Fresh Lines (`N >= 3`)* | — | *`38.5 min`* | *`27.5 min`* | *`22.3 min`* | ***`1.47 hrs`*** |
| `a549_ko_cas9pos` (Optional R16 Subset Sensitivity) | `403` (`404`) | `5.5 min` | `19.1 + 4.8 = 23.9 min` | `15.4 + 3.8 = 19.2 min` | `48.6 min` |
| `wtc11_i3neuron` (`N = 1`, Excluded from Inference) | `1` (`2`) | `4.5 min` | `0.1 min` | `0.1 min` | `4.7 min` |
| **Grand Total (All 9 Local Lines, If All Executed)** | — | **`48.6 min`** | **`51.6 min`** | **`41.6 min`** | **`2.36 hrs`** *(+ `~10–15 min` raw disk I/O & `~3 min` 4,000 bootstraps $\approx$ **`2.6 hrs`**)* |

- **Peak RAM & Hardware Safety:** Peak RSS during `SE-600M` (`MPS float32`, `bs=16`) is **`~5.70 GB` (`< 4.5%` of `128 GB` unified memory)**; peak RSS during `ST-SE` / `ST-HVG` is **`~2.35 GB`**.
- **Strict Operational Guardrails:**
  1. **Zero Concurrent MPS Execution:** Main Line and `A-T3` NEVER run MPS jobs concurrently.
  2. **10-Minute MPS Stall Fallback (Matching `PREREG_I Addendum B`):** If any `MPS float32` step stalls without progress for `> 10 min`, abort cleanly and switch BOTH parties identically to `CPU float32` (logging the deviation transparently).
  3. **Array Storage Footprint:** Storing all 56 predicted delta matrices `(N_eval, G_eval)` in `float32` `.npz` for `cd4t_stim48` (`401 × 1,836 × 56 × 4 B ≈ 165 MB` uncompressed per model, `~330 MB` for `ST-SE + ST-HVG`) plus secondary lines (`< 50 MB`) totals **`< 0.5 GB`** (`< 0.85 GB` if `a549_ko_cas9pos` is also exported) — trivial local disk footprint with full bit-level auditability.
  4. **Precision Lock:** `MPS float32` ONLY (`bf16`/`fp16` prohibited on Apple Silicon MPS due to Metal `MPSNDArrayMatrixMultiplication.mm:5799` assertion verified in R16 Gate 2).

### §5.2 Independent Re-Implementation & Staged Three-Party Verification Plan (`A-T3` After R16 Closes)
**Sequencing Lock:** `A-T3` will NOT begin Round 17 code implementation or execution until **Round 16 is formally closed** by Master Controller. Once Round 16 closes and `PREREG_J` passes lock (`PREREG_J_COMMIT.txt`), code freeze (`CODE_FREEZE_R17.txt`), and execution release (`R17_GO.txt`), Round 17 enforces the exact three-party verification architecture established in `PREREG_I Addendum A §7 & §11`:

1. **Clean-Room Dual Implementation:**
   - **Main Line** implements `round17/src/` (`run_eval_r17.py`, `metrics_r17.py`, `adjudicate_r17.py`, `test_unit_r17.py`) with a **triple hard gate** refusing real-data execution unless `PREREG_J_COMMIT.txt`, `CODE_FREEZE_R17.txt`, and `R17_GO.txt` all exist and are non-empty.
   - **`A-T3`** independently implements `round17/audit/` strictly from `PREREG_J.md` (sharing only the frozen upstream `arc-state 0.11.1` library collator/model loader and frozen `round10/src/data_loader_r10.py` `HGNCBridge`, zero shared Round 17 analysis code), verified first on pure-synthetic smoke fixtures (`PYTHONHASHSEED=1` vs. `99991` cross-process determinism: `0.0` on CPU, `<= 1e-6` on MPS).
2. **Staged Three-Party Tolerance Ladder (Checked Strictly in Order by Master Controller):**

| Verification Stage | Artifacts Compared (Main Line vs. `A-T3` / Controller) | Locked Acceptance Tolerance |
|---|---|---|
| **Stage 1: Control-Cell Barcodes** | `control_barcodes_r17.json` & `.npz` `control_barcodes` (exact Hamilton largest-remainder stratified sample, `seed = 20261701`) | **`100%` exact ordered string equality** |
| **Stage 2: `SE-600M` `X_state` Embeddings** | `round17/results/embeddings/x_state_<line>.npz` (`float32`, official `FilteredGenesCounts` + `VCIDatasetSentenceCollator`, per-pool `seed = 20261701`, checked **before** any transition pass) | **$\max_{i, d} \left|X_{\mathrm{state}}^{\mathrm{Main}} - X_{\mathrm{state}}^{\mathrm{AT3}}\right| \le 10^{-5}$** |
| **Stage 3: Predicted $\Delta$ Arrays (`56` Indices + `L1`/`L2`/`G4`)** | Exported per-line, per-model, per-index `.npz` arrays (`preds_STATE_{SE,HVG}_b{00..55}.npz`, `preds_*_L1_ens56.npz`, `preds_*_G4_*.npz`, plus `'non-targeting'` pseudobulks for `L2_ntc_recon_match`) | **$\max_{p, g} \left|\hat{\Delta}^{\mathrm{Main}} - \hat{\Delta}^{\mathrm{AT3}}\right| \le 10^{-4}$** *(independent `MPS float32` passes)* |
| **Stage 4: Controller Same-Array Recomputation** | Master Controller recomputing every `H-J1`, `H-J2`, `H-J3`, and `G4` statistic directly from exported `.npz` arrays (`seed = 20261701`, `n_boot = 4,000`) | **$\left|\mathrm{metric}^{\mathrm{Controller}} - \mathrm{metric}^{\mathrm{Party}}\right| \le 10^{-10}$** |
| **Stage 5: Cross-Party Metrics & Verdicts** | Main Line vs. `A-T3` metrics across all 56 indices, $\widehat{\Delta}_{\mathrm{CF}, 0}$, $\hat{\tau}_{\mathrm{batch}}$, $p_{\mathrm{var}}$, $f_{\mathrm{flip}}$, $\hat{b}_{\mathrm{NTC}}$, `L1`/`L2` contrasts, and `G4` 2×2 | **$\left|\mathrm{metric}^{\mathrm{Main}} - \mathrm{metric}^{\mathrm{AT3}}\right| \le 10^{-5}$** **AND `100%` verdict parity** (`H-J1`, `H-J2`, `H-J3`, exact $\hat{b}_{\mathrm{NTC}}$ index match) |

*(Any discrepancy at any stage halts adjudication and must be reported item-by-item with full root-cause disclosure; selective omission is strictly prohibited.)*

---

## §6. Open Review Items Queued for Master Controller & Human Project Lead Prior to Lock

1. **Confirmatory Cohort Scope (`§2.2`):** Confirm whether Round 17 locks **Track 1** (zero-download local `K_conf = 1` [`cd4t_stim48`, `N = 401`, mandatory shared-study disclosure with R16 `cd4t_rest`] + secondary fresh lines `a549_myllia` `N = 20`, `hek293_xucao` `N = 16`, micro-N descriptive lines, and separated `a549_ko_cas9pos` subset sensitivity) OR **Track 2** (adding a separately human-approved prospective holdout such as `GSE295214` `H1 hESC` or future public VCC 2026 cohorts to reach strictly independent $K_{\mathrm{conf}} \ge 2$).
2. **Micro-N & `wtc11_i3neuron` Execution Scope (`§5.1`):** Confirm whether to run the full 56-index sweep on the 4 micro-N lines (`telohaec` `N=9`, `calu3` `N=8`, `mac_lps` `N=4`, `mac_ctrl` `N=3`; adds only `~25 min`) and `a549_ko_cas9pos` subset sensitivity (`+48.6 min`), while omitting inference on `wtc11_i3neuron` (`N = 1`).
3. **Rule `L2_ntc_recon_match` Reference Metric (`§3.4`):** Confirm Pearson correlation on `G_eval(t)` between 4-fold mean predicted `'non-targeting'` and observed `official_median_log1p` NTC mean as primary `L2`, with negative L2 distance and per-fold index selection logged as secondary diagnostics.
