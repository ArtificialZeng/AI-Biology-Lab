# Round 17 Prior-Art Review: Batch / Covariate / Context Tokens at Zero-Shot Inference in Perturbation-Response Models (`PRIOR_ART_J.md`)

- **Protocol Served:** [`round17/prereg/PREREG_J_DRAFT.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round17/prereg/PREREG_J_DRAFT.md) (Confirmatory Audit of `batch_encoder` Index Dependence in Zero-Shot `Arc STATE`)
- **Date:** `2026-11-27 (+08:00)`
- **Author:** Main Line (`Antigravity`)
- **Upstream Companion Prior-Art Logs (Unmodified):** [`round16/status/PRIOR_ART_I.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/status/PRIOR_ART_I.md), [`round16/status/PRIOR_ART_I_CONTROLLER_CHECK.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/status/PRIOR_ART_I_CONTROLLER_CHECK.md), [`round16_addB/PRIOR_ART_I_ADD2_v2.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16_addB/PRIOR_ART_I_ADD2_v2.md).
- **Evidence, Blinding, & Non-Priority Standard:**
  1. Every cited work below carries an explicit resolving primary URL/DOI **and** an explicit reading status (`VERIFIED_PRIMARY` vs. `UNREAD` / `未能读取` whenever a full text is behind Cloudflare `HTTP 403`, paywall `HTTP 429`, or was not retrieved in full). Unread full texts are NEVER reconstructed from memory.
  2. Public Virtual Cell Challenge 2026 (VCC 2026) references are drawn **strictly from public webpages, bioRxiv abstracts, and public GitHub repositories** — zero VCC participant-channel data or private communications are used.
  3. **Zero Priority Claims:** This review and `PREREG_J_DRAFT.md` make **NO priority claims** (never claiming to be the earliest or sole audit of batch/covariate conditioning or input sensitivity). Direct methodological and probing precedents (`CELLAUDIT`, `Zhang et al. 2026`, `Molina & Zhang 2026`, `VCBench`, `OCOO-T`) are explicitly credited.
  4. **Hard-Guard Compliance:** Zero inference, zero MPS/GPU work, zero downloads, zero reads of `round16/audit/` or `A-T3.txt`, zero git commands.

---

## §1. Executive Summary of the Target Mechanism & Prior-Art Landscape

Round 17 investigates a specific architectural and operational mechanism in released zero-shot `Arc STATE` checkpoints (`ST-SE-Replogle` and `ST-HVG-Replogle`, `arc-state 0.11.1`):

- Every released canonical zero-shot checkpoint enables `batch_encoder = nn.Embedding(56, 128)` (`batch_col = gem_group`, `use_batch_token = false`, `batch_predictor = false`), adding `self.batch_encoder(batch_indices.long())` directly to `pert_embedding + control_cells` prior to the 8-layer Llama backbone (`state/tx/models/state_transition.py:192–197, 422–441`).
- Read-only weight inspection across all 8 checkpoints (`round16/shared/models/ST-{SE,HVG}-Replogle/zeroshot/{hepg2,jurkat,k562,rpe1}`) shows that the 56 rows retain near-initialization $\mathcal{N}(0, 1)^{56 \times 128}$ norms (`9.77–12.43`, median `11.09` in `ST-SE`; `8.81–12.52`, median `10.88` in `ST-HVG`), remain nearly mutually orthogonal within each fold (off-diagonal row cosine `mean = 0.0014, sd = 0.0919` in `ST-SE`; `0.0006 ± 0.0896` in `ST-HVG`), and remain tightly aligned across all 4 zero-shot folds (`mean_cos = 0.9995–0.9997`), while carrying `>98%` of input L2 energy relative to `basal_encoder` (`~0.38–1.20`) and `pert_encoder` (`~0.74`).
- At zero-shot inference on unseen cell lines, official `arc-state 0.11.1` (`state/_cli/_tx/_infer.py:678–697`) silently falls back to **`batch_idx = 0`** (`gem_group = '1'`) whenever a batch column exists with unseen labels (`idxs[i] = 0  # fallback to zero`, line `688`), and raises `KeyError: 'batch'` (`state_transition.py:427`) if no batch column is present.

Across all verified primary sources reviewed below:
1. **Batch/covariate conditioning and probing in perturbation models are well-established topics** (`Zhang et al. 2026` §C.5 `control + batch` ridge probes; `VCBench` `DecoderOnly` covariate baseline; `OCOO-T` discrete cell-line vs. mean-control comparison; `CPA` and `LPM` discrete context embeddings; `CELLAUDIT` frozen-checkpoint input replacement).
2. **None of the verified primary sources** systematically sweep all 56 trained `batch_encoder` indices (`0..55`) of frozen zero-shot `Arc STATE` (`ST-SE` and `ST-HVG`) on fresh unseen cell lines, quantify `STATE - FCT` decision-call instability across the 56 indices under a pre-registered confidence-interval rule, or benchmark pre-specified label-free index-selection / index-ensembling rules against `fixed_index_zero` and a cross-fitted oracle. *(Any unreadable full texts listed in §4 are explicitly flagged as potential partial overlap rather than assumed absent.)*

---

## §2. Core Primary Sources (Directly Verified vs. Explicitly Marked `UNREAD`)

### §2.1 `Arc STATE` Paper, Documentation, and Official `arc-state 0.11.1` Source Code

| Item | Resolving Primary Link / Local Verification Path | Reading Status | Verified Technical Content |
|---|---|---|---|
| **Published paper (`Cell`, 2026 Sep 17)** | <https://doi.org/10.1016/j.cell.2026.07.052> (PubMed: <https://pubmed.ncbi.nlm.nih.gov/42673963/>; preprint <https://doi.org/10.1101/2025.06.26.661135>) | **Abstract: `VERIFIED_PRIMARY`**; **Full text: `UNREAD` (`未能读取`: Cloudflare `HTTP 403` / Elsevier `HTTP 429`)** | Formal STATE citation (Adduri et al., *Cell* 2026). Full paper prose on `batch_encoder` / unseen-batch inference is `UNREAD`; we rely strictly on the directly verified installed `arc-state 0.11.1` code and shipped HuggingFace configs below. |
| **Arc release blog & GitHub repo** | <https://arcinstitute.org/news/blog/state> ; <https://github.com/ArcInstitute/state> | **`VERIFIED_PRIMARY`** (blog & local installed `arc-state 0.11.1` package) | Public `ArcInstitute/state` documentation and `cell-load` tutorials specify `data.kwargs.batch_col=gem_group` and `model.kwargs.batch_encoder=true` when training multi-batch transition models. |
| **Shipped HF model repos (`ST-SE-Replogle`, `ST-HVG-Replogle`)** | <https://huggingface.co/arcinstitute/ST-SE-Replogle> ; <https://huggingface.co/arcinstitute/ST-HVG-Replogle> ; local verified checksums (`77/77 OK`) in `round16/shared/models/` | **`VERIFIED_PRIMARY`** | All 8 canonical zero-shot checkpoints (`hepg2`, `jurkat`, `k562`, `rpe1` × `ST-SE`, `ST-HVG`) set `batch_encoder: true`, `use_batch_token: false`, `batch_col: gem_group`, `batch_dim: 56` (`batch_onehot_map.torch` maps `'1'..'56'` lexicographically `'1'->0, '10'->1, ..., '9'->55`). No per-`gem_group` training-cell embeddings or `gem_group`-to-cell-line provenance table are shipped (`未能确定`). |
| **Exact unseen-batch fallback (`arc-state 0.11.1`)** | Installed `round16/.venv/lib/python3.11/site-packages/state/_cli/_tx/_infer.py` (lines `656–698, 822–825`) and `state/tx/models/state_transition.py` (lines `188–230, 422–448`) | **`VERIFIED_PRIMARY`** (exact source code read) | `_infer.py:683–693`: if `adata.obs` contains a batch column, each cell's string label is looked up in `label_to_idx`; any unseen batch label prints a warning and executes `idxs[i] = 0  # fallback to zero` (line `688`), indexing `batch_encoder.weight[0]` (`gem_group = '1'`). If no batch column exists (`_infer.py:694–697`), `uses_batch_encoder = False` omits `batch["batch"]`, raising `KeyError: 'batch'` at `state_transition.py:427`. `state_transition.py:440–441` adds `batch_embeddings` directly to `pert_embedding + control_cells`. |

---

### §2.2 `Zhang et al. (2026)` — `arXiv:2611.05849` (§C.5 `control + batch` Ridge Probe & `L0` Decoder Swap)

- **Citation:** Zhang et al. (2026), *"A low-rank decoder bottleneck bounds reliability in foundation-model perturbation prediction"* (NeurIPS 2026).
- **Resolving Link:** <https://arxiv.org/abs/2611.05849>
- **Reading Status:** **`VERIFIED_PRIMARY` (relevant sections & Table 7 verified from arXiv HTML via `PRIOR_ART_I_CONTROLLER_CHECK.md` & `PRIOR_ART_I.md`; secondary full-appendix re-download not performed under Round 17 zero-download rule).**
- **Directly Verified Verbatim Quotes:**
  1. Context/batch conditioning architecture:
     > *"the decoder receives 64 control cells from the same batch (gem group) as context for each perturbation"*
  2. §C.5 auxiliary linear probe on control and batch inputs:
     > *"A ridge on the control and batch inputs alone reaches PDC 0.27–0.46"*
  3. Table 7 cross-cell-line internal direction swap caption:
     > *"Change in mean PDC when the top L0 direction of each perturbation is replaced by the one the same encoder gives it in a decoder trained on another cell line."*
- **Relevance & Precise Boundary vs. Round 17 (`PREREG_J`):**
  - Zhang et al. §C.5 establish directly that `(control, batch / gem_group)` inputs alone predict substantial perturbation-discrimination signal in `STATE` decoders (and also evaluate `perturbation + batch` and shuffled-perturbation ridge probes).
  - However, Zhang et al. evaluate linear probes trained on `(control, batch)` and swap internal layer-0 (`L0`) SVD directions across separate per-cell-line decoders on the Replogle 4 lines (`K562, RPE1, HepG2, Jurkat`). They do **not** sweep the frozen `batch_encoder = nn.Embedding(56, 128)` lookup index `b ∈ {0..55}` on fresh unseen cell lines, do **not** quantify `STATE - FCT` ranking flips across batch indices (`H-J2`), and do **not** evaluate label-free unseen-batch inference rules (`H-J3`) or compare `ST-SE` vs. `ST-HVG` under the same 56-index protocol.

---

### §2.3 `Molina & Zhang (2026)` — Perturbation Response Decomposition

- **Citation:** Alexis Molina & Xinyi Zhang (2026), *"Perturbation response decomposition enables biologically aligned generalization to unseen perturbations and cellular contexts"*, `bioRxiv 10.64898/2026.07.24.740459` (Europe PMC `PPR1286577`).
- **Resolving Links:** <https://doi.org/10.64898/2026.07.24.740459> ; official MIT repo <https://github.com/xinyizhanglab/perturbation-decomposition>
- **Reading Status:** **bioRxiv REST API abstract & GitHub source (`README.md`, `decomposition/anova.py`, `evaluation/metrics.py`, `models/state/predict.py`): `VERIFIED_PRIMARY`**; **bioRxiv preprint full text/PDF: `UNREAD` (`未能读取`, Cloudflare `HTTP 403`)**.
- **Directly Verified Verbatim Quotes (from `PRIOR_ART_I_ADD2_v2.md`):**
  > *"The global response component is low-dimensional, reflects recurrent proliferation and stress response programs, and can be inferred from control gene expression. In contrast, the perturbation and cell-line specific components are high-dimensional and cannot be recovered from control expression."*
  > *"STATE | Source-only transfer"* / *"STATE | `arc-state==0.10.2`; Python 3.11, PyTorch 2.10.0+cu128"*
- **Relevance & Boundary vs. Round 17:** Evaluates `STATE` (`arc-state==0.10.2`) via four-component observational ANOVA on the 4 Replogle/Nadig lines. Because the preprint full text is `UNREAD`, we note explicitly that any appendix discussion of `STATE` covariate settings in the preprint body cannot be verified from the abstract alone; their public GitHub code contains no 56-index `batch_encoder` sweep, no `STATE - FCT` call-instability test, and no label-free batch-index selection benchmark.

---

### §2.4 `CELLAUDIT` (`Li et al. 2026`) — Direct Methodological Precedent for Frozen-Checkpoint Input Audits

- **Citation:** Mengran Li, Bo Li, Chengyang Zhang, Yang Yan, Jinfeng Xu, & Zhenchao Tang (2026), *"Discover, Falsify, Revise: Auditing Input-Use Claims from Source Code to Predictive Contribution in Agent-Discovered Cell Models"*, `arXiv:2609.27234v1`.
- **Resolving Link:** <https://arxiv.org/abs/2609.27234>
- **Reading Status:** **`VERIFIED_PRIMARY` (full `230,078`-char arXiv HTML read in `PRIOR_ART_I_ADD2_v2.md`)**.
- **Directly Verified Verbatim Quotes:**
  > *"To distinguish these cases, CellAudit combines source inspection with matched input replacement (Fisher et al., 2019; Chamma et al., 2023) on frozen checkpoints, testing both whole-model dependence and target-relevant contribution."*
  > *"The plate-family audit crosses the observed or permuted compound representation with the observed or an alternative control profile from the same plate family."*
- **Relevance & Boundary vs. Round 17:** Direct methodological precedent for pairing source-code inspection with frozen-checkpoint input/covariate replacement and factorial crossing (`X_ctrl × batch` in our secondary `G4` mechanism test is directly analogous in spirit to CELLAUDIT's crossed factorial audits). Full-text verification confirms `0` mentions of `Arc STATE` or discrete `batch_encoder` index sweeps; we credit CELLAUDIT explicitly and claim no novelty for frozen-checkpoint input-perturbation auditing itself.

---

### §2.5 Public Virtual Cell Challenge 2026 (VCC 2026) Write-Ups (Zero Participant-Channel Data)

| Source | Resolving Public Link | Reading Status | Verified Findings on Batch / Covariate / Context Handling |
|---|---|---|---|
| **Arc Institute `VCC 2026 Wrap-Up` (`2026-11-19`)** | <https://arcinstitute.org/news/virtual-cell-challenge-2026-wrap-up> | **`VERIFIED_PRIMARY`** | Summarizes 2026 zero-shot multi-context challenge (`CellularAI` #1, `flowstate` #2, `X-Cell Transfer` #3); reports small cell-line-specific slopes (`0.02–0.13` vs. `1.02` replicate) and notes covariate-adjusted log-fold-change statistical transfer outperformed end-to-end neural models. Contains **no** audit of `STATE`'s discrete `batch_encoder` lookup table (`0..55`) or `_infer.py:688` unseen-batch fallback. |
| **Dobin et al. (2026) — `cell-eval2` / VCC 2026 metrics** | <https://doi.org/10.64898/2026.11.17.767214> | **API abstract & `vcc2026-metrics-brief.md`: `VERIFIED_PRIMARY`**; **bioRxiv full text: `UNREAD` (`未能读取`, Cloudflare `HTTP 403`)** | Defines 6 replicate-calibrated scoring metrics anchored on target control cells; does not evaluate `STATE`'s batch-index dependence. |
| **`Illumina/CellularAI` public repo & blog** | <https://github.com/Illumina/CellularAI> ; <https://www.illumina.com/science/genomics-research/articles/introducing-the-cellularai-algorithm--predicting-single-cell-res.html> | **`VERIFIED_PRIMARY`** (`README.md`, `docs/input-files.md`, public blog) | Statistical transfer pipeline using control-context pooling, covariate-adjusted summary statistics, and control-conditioned residuals (proprietary Billion Cell Atlas summary bundle `UNREAD` / access-controlled). Does not use discrete training-batch lookup tokens at inference. |
| **`Occam` (Bendidi et al., 2026)** | <https://arxiv.org/abs/2611.21959> ; <https://github.com/IhabBendidi/Occam> | **`VERIFIED_PRIMARY`** (arXiv HTML & repo) | Conditions zero-shot predictions directly on unlabeled target control cells without a discrete training-batch embedding table (`batch_encoder`) at inference. |

*Note on Public VCC 2026 Ecosystem:* Web search across public VCC 2026 write-ups (`arcinstitute.org`, `virtualcellchallenge.org`, `github.com/ArcInstitute/state`, `github.com/Illumina/CellularAI`, `bioRxiv`) surfaced no public VCC 2026 post-mortem or competitor blog analyzing `Arc STATE`'s 56-index `batch_encoder` fallback sensitivity. *(Private Discord/Slack/participant channels were strictly never accessed.)*

---

## §3. Broader Perturbation-Model Context, Batch, & Covariate Token Literature

| Work & Citation | Resolving Primary Link | Reading Status | Role of Batch / Context / Covariate Tokens & Distinction from Round 17 |
|---|---|---|---|
| **`VCBench`** — Mao et al. (2026) | <https://arxiv.org/abs/2604.27646> | **`VERIFIED_PRIMARY`** (§4.3 in `PRIOR_ART_I.md`) | Includes `DecoderOnly` baseline predicting solely from perturbation + covariate identities without control cells, establishing that covariate identities alone can drive non-trivial scores; does not sweep `STATE`'s `batch_encoder` indices `0..55`. |
| **`OCOO-T`** — Jiang et al. (2026) | <https://arxiv.org/abs/2606.12838> | **`VERIFIED_PRIMARY`** (§4.4 in `PRIOR_ART_I.md`) | Explicitly documents the trade-off between discrete learned cell-line embeddings (`S1`, strong on seen contexts, inapplicable zero-shot) and mean control-cell profile conditioning (`S2–S7`, transferable to unseen contexts). Directly parallels Round 17's question of how discrete `batch_encoder` tokens behave when transferred zero-shot to unseen lines. |
| **`LPM`** — Miladinovic et al. (2025) | <https://arxiv.org/abs/2503.23535> | **`VERIFIED_PRIMARY`** (arXiv HTML in `PRIOR_ART_I_ADD2_v2.md`) | Disentangles discrete context embeddings during training and explicitly notes *"inability to extrapolate to unseen contexts in a zero-shot manner"*. |
| **`CPA`** — Lotfollahi et al. (*Mol. Syst. Biol.* 2023) | <https://doi.org/10.15252/msb.202211517> | **`UNREAD` (`未能读取` in this session; canonical published reference only)** | Foundational compositional perturbation autoencoder with discrete additive covariate embeddings (cell type / batch / dose); known architectural predecessor for additive covariate embeddings, whose discrete lookup tables do not natively generalize to unseen covariate categories without retraining/adaptation. |
| **`TxPert`** — Wenkel et al. (2025) | <https://arxiv.org/abs/2505.14919> | **`VERIFIED_PRIMARY`** (§4.1.3 arXiv HTML) | Training-time ablation showing removing basal state encoder outperformed an MLP basal encoder on unseen K562; does not evaluate inference-time batch-index variation on frozen `STATE`. |
| **`Fan et al. (2026)`** — *Zero-Shot Cross-Cell-Line Perturbation Prediction Needs Controls* | <https://github.com/xyfan22/crossline-controls> ; <https://openreview.net/forum?id=lnpzb1HvwX> | **GitHub `README.md`: `VERIFIED_PRIMARY`**; **OpenReview PDF: `UNREAD` (`未能读取`, `HTTP 403`)** | Evaluates `target-control-copy` and `shifted-control` baselines against `CellFlow` and scFM flow models; does not sweep `STATE`'s `batch_encoder` index. |
| **`COMPASS`** — Liang & Singh (2026) | <https://doi.org/10.64898/2026.08.03.742643> | **API abstract: `VERIFIED_PRIMARY`**; **bioRxiv full text: `UNREAD` (`未能读取`, `HTTP 403`)** | Two-component linear perturbation decomposition outperforming `STATE`, `CPA`, `scGPT`, `GEARS`; does not audit inference-time `batch_encoder` fallback indices. |
| **`Systema`** — Viñas Torné et al. (*Nat. Biotechnol.* 2025) | <https://doi.org/10.1038/s41587-025-02777-8> | **`VERIFIED_PRIMARY`** (PMC13271886 in `PRIOR_ART_I.md`) | Origin-shift perturbation evaluation framework (`O_pert`); provides context for `r_dorth_cal` sensitivity metric (no novelty claimed for `r_dorth_cal`). |
| **`Nicol et al. (2026)`** | <https://doi.org/10.64898/2026.05.07.723486> | **`VERIFIED_PRIMARY`** (PMC13192888 in `PRIOR_ART_I.md`) | Establishes control-splitting to eliminate shared-control correlation bias; informs split-half noise-ceiling design. |

---

## §4. Explicit Unreadable / Gated Primary Source Register (`UNREAD` / `未能读取`)

To maintain strict epistemic transparency, every primary document whose full text could NOT be read in this audit is listed here; Round 17 makes no claim about details inside their unread bodies:

1. **Adduri et al. (*Cell* 2026, `doi:10.1016/j.cell.2026.07.052`; bioRxiv `10.1101/2025.06.26.661135`) full text:** Cloudflare `HTTP 403` / Elsevier `HTTP 429` (`UNREAD`; abstract, Arc blog, HF configs, and installed `arc-state 0.11.1` code verified).
2. **Molina & Zhang (2026, `bioRxiv 10.64898/2026.07.24.740459`) full text:** Cloudflare `HTTP 403` (`UNREAD`; bioRxiv API abstract & public GitHub repo verified).
3. **Dobin et al. (2026, `bioRxiv 10.64898/2026.11.17.767214`) full text:** Cloudflare `HTTP 403` (`UNREAD`; bioRxiv API abstract & public `vcc2026-metrics-brief.md` verified).
4. **COMPASS (Liang & Singh 2026, `bioRxiv 10.64898/2026.08.03.742643`) full text:** Cloudflare `HTTP 403` (`UNREAD`; bioRxiv API abstract verified).
5. **Fan et al. (2026, OpenReview `lnpzb1HvwX`) PDF/forum:** OpenReview `HTTP 403` (`UNREAD`; GitHub `xyfan22/crossline-controls/README.md` verified).
6. **CPA (Lotfollahi et al. 2023, `doi:10.15252/msb.202211517`) full text:** Not fetched under Round 17 zero-download boundary (`UNREAD` in this session; cited strictly for canonical discrete additive covariate-embedding architecture).
7. **Arc `ST-*-Replogle` training `gem_group` `'1'..'56'` provenance mapping & `SE-600M` `14,420`-dataset CSV manifest:** Not shipped in any HuggingFace repository (`未能确定`).

---

## §5. Summary Positioning for `PREREG_J`

1. **Credited Precedents:** `CELLAUDIT` (`arXiv:2609.27234`) establishes frozen-checkpoint input/covariate replacement auditing; `Zhang et al. (2026)` (`arXiv:2611.05849` §C.5) establishes that `(control, batch / gem_group)` inputs carry major linear predictive power in `STATE` decoders; `OCOO-T` (`arXiv:2606.12838`), `VCBench` (`arXiv:2604.27646`), and `LPM` (`arXiv:2503.23535`) document discrete covariate/context-embedding limitations under zero-shot transfer; `Molina & Zhang (2026)` (`10.64898/2026.07.24.740459`) decompose `STATE` responses on Replogle/Nadig.
2. **Targeted Gap Tested in Round 17 (`PREREG_J`):** No verified primary source performs a pre-registered, winner's-curse-corrected full-table sweep of all 56 `batch_encoder` indices (`0..55`) across both `ST-SE-Replogle` and `ST-HVG-Replogle` on fresh never-STATE-run lines to quantify (a) true cross-index variation in `r_orth_cal` (`H-J1`), (b) `STATE - FCT` decision-call instability relative to `arc-state 0.11.1`'s default `fixed_index_zero` fallback (`H-J2`), and (c) whether pre-specified label-free index rules outperform `fixed_index_zero` (`H-J3`). Zero priority or novelty-first claims are made.
