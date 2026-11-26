# Round 16 Pre-Lock Addendum B Draft (`v2`) — Secondary Open-Weight Model: `Occam` (`PREREG_I_ADDENDUM_B_DRAFT_v2.md`)

> **STATUS: `v2` PRE-LOCK DRAFT FOR MASTER CONTROLLER (`Claude`) AUDIT & EXPLICIT HUMAN LICENSE/DOWNLOAD APPROVAL (`R16_ADDB_HUMAN_APPROVAL.txt`). DO NOT LOCK WITHOUT CONTROLLER SIGN-OFF.**
> **Addresses Master Controller Review R1 ([`ADDB_CONTROLLER_REVIEW_R1.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16_addB/ADDB_CONTROLLER_REVIEW_R1.md), Items `P1–P6`):**
> - **`P1`:** Integrates `Molina & Zhang 2026` (`bioRxiv 10.64898/2026.07.24.740459`, `xinyizhanglab/perturbation-decomposition`), `COMPASS` (`Liang & Singh 2026`, `bioRxiv 10.64898/2026.08.03.742643`), `CELLAUDIT` (`Li et al. 2026`, `arXiv:2609.27234v1` — explicitly credited as direct methodological precedent for frozen-checkpoint input replacement), and `Stańczak & Niazi 2026` (*Cells*, `doi:10.3390/cells15181711`, `PMC13605750`), and removes all unverified secondary-search OpenReview citations; drops any wording implying the input-replacement design itself is novel (`§7.2` & companion `PRIOR_ART_I_ADD2_v2.md`).
> - **`P2`:** Prohibits modifying frozen `round16/.venv`; pre-specifies isolated `round16_addB/.venv` from Occam's `pyproject.toml` pins with exact macOS `arm64` fallbacks + `CROSS_VERSION_ATOL/RTOL`; locks Occam source archive of commit `9ddb890e...` (code only, no LFS, SHA-256 recorded) under human approval alongside weights/assets; specifies shared unchanged upstream model classes + independent clean-room adapters (`§1`, `§6.1–§6.2`).
> - **`P3`:** Reconciles `621` held-out test perturbations vs. `611` split-half cohort (`min_cells = 10` dropping `10` shallow targets); documents exactly what `splits_manifest.json` and `cohort.json` contain vs. `splits.npz` (`600,792 B`); reads HepG2 strictly from local `/Users/mac/.cache/r9_p2_raw/nadig2025/GSE264667_hepg2_raw_singlecell_01.h5ad` (SHA-256 `e1ad7c3c...`, not sync share); locks dual Gate 5B tolerances (`±0.005 / ±0.002` if `20/20` split SHA-256s match `splits_manifest.json`; `±0.015 / ±0.005` only if fresh deterministic splits are required) (`§2.1`, `§5.2`).
> - **`P4`:** Pre-registers mandatory **`pairwise_common_axis`** sensitivity arm zeroing union-absent RHJ axis columns in both matched and mismatched inputs for every pair involving `cd4t_rest` (`413` absent) or `a549_ko` (`33` absent), with verified exact common-axis counts (`6,557`, `6,570`, `6,950`) (`§3.3`, `§4.1`).
> - **`P5`:** Corrects Eq. 13 notation to explicit **per-cell encoding** `z_i = E_theta(x_i)` followed by averaging per-cell decoded outputs (`y_hat_{i,p}`), matching `src/occam/vcell/model.py` and `src/occam/api.py`, and states explicit inference determinism (zero sampling, softplus mean output, dropout off) (`§3.1`, `§3.4`).
> - **`P6`:** Retains strict no-concurrent-MPS rule and pre-locks the **10-minute no-progress MPS stall rule** switching both parties to `CPU float32` for the entire Addendum B run (`§6.1`).
>
> **Parent Protocols (Strictly Unmodified):**
> - [`round16/prereg/PREREG_I.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md) (`v5-FINAL-PRE-LOCK`, commit `c1fc7c6d24dd128bc5e2d969d8c78224a2cbc118`, SHA-256 `43d242e460fd64db2e4af1e628e886cf1d94fb4ccd1e95925928a99cb96ee4a8`);
> - [`round16/prereg/PREREG_I_ADDENDUM_A.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I_ADDENDUM_A.md) (Items 1–11, SHA-256 `bdb16818c6e2499d1eb10f8690f3905dd19ce1d7af4c423c9f5afe5455ed3a3b`).
> **Companion Prior-Art Refresh (`v2`):** [`round16_addB/PRIOR_ART_I_ADD2_v2.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16_addB/PRIOR_ART_I_ADD2_v2.md).
> **Supersedes:** [`PREREG_I_ADDENDUM_B_DRAFT.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16_addB/PREREG_I_ADDENDUM_B_DRAFT.md) (`v1`, retained untouched for audit trail).
>
> **Pre-Download, Blinding, & Freeze Compliance Attestation:**
> 1. Zero Occam weights (`.pth`), feature tables (`.pkl`), binary panels (`.npz`), code archives, or `GSE295214` files downloaded; zero Occam code cloned or executed; zero GPU/MPS compute run.
> 2. Zero modifications inside `round16/` (except appending `R16_ADDB_V2_READY` to `round16/status/MAIN.txt`); zero reads of `round16/audit/`, `A-T3` logs, or `A-T3.txt`; zero git operations; zero priority claims or outreach text.

- **Secondary Audited Model:** **Occam** (Bendidi, Kenyon-Dean, Denton, Genovesio, & Noutahi, `arXiv:2611.21959`, 2026-11-16) — `rhj_fewshot` all-source virtual-cell predictor (`ABL_N_ALL`, `seed_42`).
- **Scientific Purpose:** Extend Round 16's matched-vs-mismatched target-control audit (`PREREG_I.md`) beyond `Arc STATE` (`STATE_SE` / `STATE_HVG`) to an architecturally independent open-weight perturbation model (MLP denoising autoencoder + multiplicative mean-free FiLM + GO/STRING/GenePT/ESM2 perturbation encoder; zero discrete cell-line or batch IDs), enabling (a) standalone `H-I1_Occam` and `H-I2_Occam` adjudication on unseen Group 2 targets (`{a549_ko, dld1, cd4t_rest}` and optional Group 3 `gse295214_h1`), and (b) a pre-registered paired **`STATE_SE` vs. `Occam`** contrast on identical `(targets, genes, control cells)`.
- **Locked RNG Seeds (Identical to Parent Protocol):** Primary `20261601`; within-line disjoint resample `20261602`.

---

## §1. Exact Artifacts to Hash at Download, Code Archive Pin, and License Disclosure (Item 1; `P2`)

### §1.1 Pinned Upstream Revisions & Code Archive (`P2`)
- **Source Code Archive (Code Only, Zero LFS Binaries):**
  - Repository: `https://github.com/IhabBendidi/Occam`
  - **Pinned Commit SHA:** **`9ddb890e65ca1ec2f81882c78311002a75b52f67`** (verified `HEAD` `2026-11-26`).
  - **Exact Download URL (after approval only):** `https://github.com/IhabBendidi/Occam/archive/9ddb890e65ca1ec2f81882c78311002a75b52f67.tar.gz` (source archive only — never `git clone`, zero Git LFS binaries). Its exact downloaded `.tar.gz` SHA-256 and extracted tree checksum MUST be recorded in `round16_addB/R16_ADDB_CHECKSUMS.sha256`.
  - **Code License:** **`MIT License`** (`LICENSE`, `Copyright (c) 2026 Occam contributors`, covering `"this software and associated documentation files (the 'Software')"`).
- **Model & Benchmark Hub Repository:** `https://huggingface.co/bendidiihab/Occam` (`type: models`, `gated: false`, `tags: ["arxiv:2611.21959", "region:us"]`)
  - **Pinned Hub Revision SHA:** **`acd080b788d028911807b87e11a56498f984f279`** (`lastModified: 2026-11-18T12:25:03.000Z`).

### §1.2 Exact Weights, Assets, and Code-Archive Human Approval Scope (`P2`)
Primary inspection of HuggingFace API metadata, `MODEL_CARD.md`, and `assets/README.md` reveals **three explicit licensing gaps on weights and derived embedding tables**:
1. **HuggingFace Hub metadata (`bendidiihab/Occam`):** `cardData` is `null` — no structured `license` field is declared.
2. **Model card (`IhabBendidi/Occam/MODEL_CARD.md:2-3`):** YAML frontmatter explicitly declares `license: other` / `license_name: see-repository`, while the repository root ships only the `MIT License` covering `"Software"`, with **no standalone weight-license text**.
3. **Derived feature tables (`bendidiihab/Occam/assets/README.md:26-33`):** Verbatim statement:
   > *"GenePT: the GenePT v2 release (133,736 gene-symbol keys, text-embedding-3-large, 3072-d) ... ESM2: precomputed ESM2 (1280-d, esm2_t33_650M) per-gene protein embeddings from the Hugging Face dataset `genbio-ai/foundation-models-perturbation` @ `d67ba210` (`ESM2_(D=1280).h5ad`) ... Licensing of both sources is not recorded here; the Hub repo is private (see the README)."*

**Mandatory Human Approval Gate (`round16_addB/R16_ADDB_HUMAN_APPROVAL.txt`):** Prior to downloading ANY artifact or building `round16_addB/.venv`, explicit human approval MUST cover **BOTH**:
- **(a) Occam model weights and derived feature assets** (`bendidiihab/Occam` @ `acd080b7...`, `license: other / see-repository`, unrecorded GenePT/ESM2 licensing); **AND**
- **(b) Occam pinned source code archive** (`IhabBendidi/Occam` @ `9ddb890e65ca1ec2f81882c78311002a75b52f67`, `MIT License`, code-only `.tar.gz`).

### §1.3 Exact Artifacts to Download and SHA-256 Hash-Lock
All 4 binary artifacts below reside under `bendidiihab/Occam` @ `acd080b788d028911807b87e11a56498f984f279` and total **`2,929,553,709` bytes (`2.93 GB` / `2.73 GiB`)**. At download time, local SHA-256 checksums MUST match the pre-locked upstream HuggingFace LFS SHA-256 digests below bit-for-bit (`round16_addB/R16_ADDB_CHECKSUMS.sha256`):

| # | Artifact / Hub Path (`bendidiihab/Occam` @ `acd080b7` or GitHub @ `9ddb890e`) | Exact Size (Bytes) | Pre-Locked Upstream SHA-256 | Role & Runtime Requirement |
|---|---|---|---|---|
| 0 | GitHub archive `IhabBendidi/Occam/archive/9ddb890e65ca1ec2f81882c78311002a75b52f67.tar.gz` | Recorded at download | Recorded in `R16_ADDB_CHECKSUMS.sha256` at download | **REQUIRED (`P2`) — Pinned source archive (code only, no LFS).** Provides unchanged upstream model classes (`AutoEncoder`, `HybridNormGlobalSkipDecoder`, `VirtualCellModel`, `FiLMInteraction`, `ConcatKnowledgeModel`) imported identically by Main Line and `A-T3`. |
| 1 | `checkpoints/occam_seeds/rhj_fewshot/all_sources/seed_42/best_model.pth` | **`2,550,066,734`** (`2.55 GB`) | `11d5ba5515e79710b3a952532c60aa29515046117292528c59d452a0937685ff` | **REQUIRED — Primary model checkpoint (`ABL_N_ALL`, `seed_42`).** Verified in `src/occam/vcell/pipeline.py::from_vcell_checkpoint:903-960` to be **100% self-contained for inference**: encoder, decoder, FiLM, `gene_names` (`44,728`), `pert2id` (`18,216`), adjacency buffers, and energy-scaled GenePT/ESM2 tables are baked into `gnn_state_dict`. |
| 2 | `assets/genept_full_occam_rhj_fewshot.pkl` | **`220,037,958`** (`220.0 MB`) | `a93e1131beb9722228e32d847ebc6f19b0e6d77c3a0a8374b0552320624181c8` | **Training provenance / optional full-bundle verification** (`17,843/18,216` keys, 3072-d GenePT v2; `resolved_config.yaml:97`; NOT re-read at inference because scaled vectors are baked into #1). |
| 3 | `assets/esm2_t33_650M_UR50D_occam_rhj_fewshot.pkl` | **`93,705,169`** (`93.7 MB`) | `cb551710feaca63c77f4e0e008064a3f7a9e645a8900c45c082ab81867ed6a16` | **Training provenance / optional full-bundle verification** (`18,149/18,216` keys, 1280-d ESM2; `resolved_config.yaml:98`; NOT re-read at inference). |
| 4 | `assets/genept_occam_rhj_fewshot.pkl` | **`65,743,848`** (`65.7 MB`) | `595068eacae03f6256003953558f67bbf0e1f03e86422f0801f40debb32df668` | **Provenance audit only** (`2,376` keys; per `assets/README.md:18`, used by PerturbDiff competitor, `max_abs_diff = 0.0` vs. #2 on shared keys; NOT read by Occam). |

- **Download Scope Option for Human Approval:** Because Artifact #1 (`best_model.pth`, `2.55 GB`) is self-contained at inference, human approval may authorize either **(a) Minimal Inference Bundle** (Artifact #0 code archive + Artifact #1 `best_model.pth`, `2.55 GB`, plus tiny text sidecars below) or **(b) Complete `rhj_fewshot` Artifact Bundle** (Artifact #0 + Artifacts #1–#4, `2.93 GB`, plus tiny text sidecars below).
- **Required Tiny Text/JSON/CSV Provenance Sidecars (`< 400 KB` total, non-LFS):**
  - `checkpoints/occam_seeds/rhj_fewshot/all_sources/seed_42/resolved_config.yaml` (`4,435 B`, git OID `51cff88a31fbcdcce6c8ba4b87983c0d6783596f`);
  - `benchmarks/rhj_fewshot/gene_axis.csv` (`156,021 B`, SHA-256 `ef64f2ba7ff7816ed947796066cd8a68014884c6176619d750ef79462c78d4c0`);
  - `benchmarks/rhj_fewshot/perturbations.csv` (`68,294 B`);
  - `benchmarks/rhj_fewshot/splithalf/cohort.json` (`101,500 B` — frozen `611`-perturbation HepG2 cohort for §5.2);
  - `benchmarks/rhj_fewshot/splithalf/splits_manifest.json` (`4,894 B` — exact `20` split SHA-256 digests and `blake2b` seed recipe for §5.2);
  - `assets/*.pkl.json` (`30,037 B` total) if Artifacts #2–#4 are included.
- **Explicitly Excluded Downloads (Never Downloaded in Addendum B):** `seed_43`–`seed_46` (`10.2 GB`), `occam_small_encoder_decoder.pth` (`893 MB`), `benchmarks/rhj_fewshot/splithalf/panel.npz` (`171 MB`, rebuilt locally from existing local HepG2 in §5.2), `benchmarks/rhj_fewshot/splithalf/splits.npz` (`600,792 B`, regenerated deterministically via `draw_splits` and verified against `splits_manifest.json`'s `20` SHA-256 hashes in §5.2), and all competitor checkpoints.

---

## §2. Evaluation Targets, Two-Stratum Mismatched Control Pools, and Pretraining Disclosures (Item 2; `P3`)

### §2.1 Training Scope of Occam `rhj_fewshot/all_sources/seed_42` (`ABL_N_ALL`) & `621` vs. `611` Reconciliation (`P3`)
Verified directly from `resolved_config.yaml:44-54`, `benchmarks/rhj_fewshot/CARD.md:12-20`, `benchmarks/rhj_fewshot/splithalf/CARD.md:34, 41-42`, and `splits_manifest.json:20-22, 90-101`:
- **Supervised perturbation training lines (3 lines):** `Replogle RPE1` (`2,381` eligible perts), `Nadig Jurkat` (`2,381` eligible perts), and `Nadig HepG2` (`1,657` train / `103` validation / **`621`** held-out eligible test perturbations), all projected onto the **`6,983`-gene Ensembl RHJ axis**.
- **Exact `621` vs. `611` Reconciliation (`P3`):** Of the `621` held-out eligible HepG2 test perturbations in `rhj_fewshot`, Occam's split-half panel builder (`scripts/build_splithalf_panel.py`, `min_cells = 10` so each half has `m >= 5`) excludes **`10` shallow test perturbations** (`CHMP3`, `COG3`, `CPSF2`, `EIF3F`, `HAUS1`, `PRPF38B`, `RFFL`, `RPF2`, `TCP1`, `ZNF676`), yielding the frozen **`611`-perturbation** split-half evaluation cohort (`cohort.json`) used for Occam's published `rhj_fewshot` scores and our Gate 5B reproduction check (§5.2).
- **Consequence for Round 16 Cohorts:**
  - **Group 1 (`nadig_hepg2`, `nadig_jurkat`) is EXCLUDED as an evaluation target for `H-I1_Occam` and `H-I2_Occam`** because both lines were directly observed during supervised perturbation training of `seed_42`. Both lines enter Addendum B strictly as (i) the **`seen_by_checkpoint` mismatched-control stratum** (`1/1` seen for this single checkpoint), and (ii) `nadig_hepg2`'s `611`-perturbation held-out cohort strictly for the pre-evaluation **Adapter Reproduction Gate (§5.2)**.
  - **Primary Evaluated Targets are Group 2 ONLY (`{a549_ko, dld1, cd4t_rest}`)**, plus **Group 3 (`gse295214_h1`)** strictly contingent on human download approval and passing `PREREG_I §2.4`'s counts-only gate, plus **`feng_ipsc`** as a **`SENSITIVITY_ONLY`** target (identical to `PREREG_I §0.2` / `Q4`).

### §2.2 Mandatory `DiverseRNA-1.4M` Encoder Pretraining Disclosure
- Occam's frozen encoder (`44,728 -> 2048 -> 1024 -> 512`) was pretrained on `DiverseRNA` (`1,434,299` profiles across 9 datasets; Kenyon-Dean et al., `arXiv:2605.31562`), with the released softplus+NB decoder finalized `2026-09-16` (`MODEL_CARD.md:23`).
- **Critical Provenance Caveat (`MODEL_CARD.md:183-186`):**
  > *"Pretraining-corpus caveats: 43.6% of rows are non-integer, 2.4% are bulk RNA-seq, and 7 of 9 constituent datasets are unidentified. Because of the last point, no claim in this card can assert that any evaluation cell line is an unseen context."*
- **Date-Based Isolation Audit for Addendum B Targets:**
  - **`a549_ko` (`GSE345058`, GEO `2026-08-26`, ENA `2026-07-03`)** and **`dld1` (`GSE337988`, GEO `2026-07-10`)** strictly post-date `DiverseRNA` (`arXiv:2605.31562`, May 2026).
  - **`cd4t_rest` (`GSE314342`, GEO `2026-02-06`)** strictly post-dates `SE-600M` (`2025-06/07`), **but precedes `arXiv:2605.31562` by ~3.5 months**. Because 7 of 9 `DiverseRNA` blocks are unattributed (`未能确定`), **encoder-pretraining absence for `cd4t_rest` cannot be formally guaranteed**, even though `cd4t_rest` was 100% unseen in Occam's supervised perturbation training (`rhj_fewshot`). This disclosure MUST accompany every `cd4t_rest` Occam table and narrative.
  - **`gse295214_h1` (`GSE295214`, GEO `2026-05-05`, ENA `2026-06-19`)** overlaps the May 2026 `DiverseRNA` publication month (supervised unseen; encoder pretraining absence `未能确定`, alongside `H1`/Replogle-Weissman disclosures from `PREREG_I §2.4.5`).
  - **`feng_ipsc` (`ENA PRJEB81502`, `2024-11-08`)** predates BOTH `SE-600M` and `DiverseRNA` $\implies$ kept strictly **`SENSITIVITY_ONLY`** as a target line, while **retained as an `unseen_by_checkpoint` mismatched-control source** (never in `rhj_fewshot` supervised training).

### §2.3 Locked Target Cohort, Stratum Assignments, and Verified Exact `(N_eval^R16, |G_eval^R16|)` Table

| Scope / Group | Target Line $t$ | Modality & Lineage | Mismatched Stratum A: `seen_by_checkpoint` ($n=2$, strictly `1/1` seen) | Mismatched Stratum B: `unseen_by_checkpoint` | Shared `R16` Locked Targets $N_{\mathrm{eval}}^{\mathrm{R16}}$ (`rhj_2381` trained + graph-extrap.) | Shared `R16` Locked Genes $\|G_{\mathrm{eval}}^{\mathrm{R16}}\|$ (`%` inside `RHJ_6983`) | Native Full `RHJ_6983 ∩ G_R9` Sensitivity Genes | Mandatory Disclosures |
|---|---|---|---|---|---|---|---|---|
| **Group 2 Primary — CRISPRko ($K_{2,\mathrm{ko}}=1$)** | `a549_ko` | Lung (**CRISPRko**) | `{nadig_hepg2, nadig_jurkat}` | `{dld1, cd4t_rest, feng_ipsc}` ($n=3$) | **`403`** (`403` hv): **`402`** trained in `rhj_2381` + **`1`** (`RBM14`) graph-extrapolated | **`1,837`** (**100.0%**, `0` missing, `0` collisions) | `6,876` | CRISPRko-vs-CRISPRi shift (`I11d`); `33` RHJ input genes absent (`P4` pairwise common axis required); `RBM14` untrained node (`402`-target sensitivity); post-dates DiverseRNA |
| **Group 2 Primary — CRISPRi ($K_{2,\mathrm{i}}=2$)** | `dld1` | Colorectal (CRISPRi) | `{nadig_hepg2, nadig_jurkat}` | `{a549_ko, cd4t_rest, feng_ipsc}` ($n=3$) | **`75`** (`75` hv): **`75/75`** trained in `rhj_2381` | **`1,832`** (**100.0%**, `0` missing, `0` collisions) | `6,883` | Low/borderline TOST power at `N = 75` (`M4`); `P4` pairwise common axis required for pairs with `a549_ko` / `cd4t_rest`; post-dates DiverseRNA |
| **Group 2 Primary — CRISPRi ($K_{2,\mathrm{i}}=2$)** | `cd4t_rest` | Primary CD4$^+$ T (CRISPRi) | `{nadig_hepg2, nadig_jurkat}` | `{a549_ko, dld1, feng_ipsc}` ($n=3$) | **`245`** (`230` hv): **`245/245`** trained in `rhj_2381` | **`1,830`** (**100.0%**, `0` missing, `0` collisions) | `6,489` | `413/6,983` RHJ input genes absent (`P4` pairwise common axis required for all pairs); `2026-02-06` precedes DiverseRNA `2026-05` (`7/9` unattributed); fragile relative $\delta$ (`M6.1`) |
| **Group 3 Prospective Holdout ($K_3 \in \{0,1\}$)** | `gse295214_h1` *(ONLY if approved & §2.4 gate passes)* | `H1 hESC` (CRISPRi `10x 5'`) | `{nadig_hepg2, nadig_jurkat}` | `{a549_ko, dld1, cd4t_rest, feng_ipsc}` ($n=4$) | Determined post-approval by `PREREG_I §2.4` $\cap \mathrm{canon}(P_{\mathrm{Occam}})$ (`>= 50` required) | `G_eval^R16(gse295214_h1) ∩ RHJ_6983` | `RHJ_6983 ∩ G_R9` | Reported strictly separately; `P4` pairwise common axis required; shared `H1` lineage (`vcc2025_h1`), Replogle/Weissman provenance, DiverseRNA May 2026 timing |
| **Sensitivity Only (`Q4`)** | `feng_ipsc` | Human iPSC pool (CRISPRi) | `{nadig_hepg2, nadig_jurkat}` | `{a549_ko, dld1, cd4t_rest}` ($n=3$) | **`88`** (`88` hv): **`88/88`** trained in `rhj_2381` | **`1,575`** (**100.0%**, `0` missing, `0` collisions) | `5,863` | `SENSITIVITY_ONLY` target (`ENA 2024-11-08` predates SE-600M & DiverseRNA); **retained as `unseen_by_checkpoint` mismatched source**; `P4` pairwise common axis required for pairs with `a549_ko`/`cd4t_rest` |

*Exact Matrix Parity:* On all three Group 2 primary targets (`a549_ko`, `dld1`, `cd4t_rest`) and on `feng_ipsc`, **100.0% of `STATE_SE`'s locked targets (`N_eval^R16`) and 100.0% of `STATE_SE`'s locked evaluation genes (`G_eval^R16`) are covered by Occam with zero dropped targets and zero dropped genes**.

---

## §3. Exact Occam Forward Definition, Control-Swap Arms, Vocabulary & Perturbation Mapping, and `P4` Pairwise Common Axis (Item 3; `P4`, `P5`)

### §3.1 Exact Per-Cell Forward Definition, Determinism, and Arm Structure (`P5`)
Verified directly from `src/occam/vcell/model.py:104-114, 293-347, 396-497, 523-574, 795-804`, `src/occam/api.py:454-522`, and `resolved_config.yaml:99-105`:
1. **Per-Cell Encoding & Decoding (Official API Convention — NO Pre-Averaging of Inputs):**
   For each individual preprocessed control cell $\mathbf{x}_i \in \mathbb{R}_{\ge 0}^{44,728}$ ($i = 1, \dots, n_{\mathrm{ctrl}}$; note that `arXiv:2611.21959 Appendix B.1` wrote $\bar{x}_{ij} = \log(1 + y_{ij})$ for per-cell log-normalized expression, **not** a sample mean across cells; we write $\mathbf{x}_i$ here to remove all ambiguity):
   $$\mathbf{z}_i = E_\theta(\mathbf{x}_i) \in \mathbb{R}^{512}, \qquad \hat{\mathbf{z}}_{i, p} = \mathbf{z}_i + \big[\sigma(g(\mathbf{e}_p)) - \sigma(g(\bar{\mathbf{e}}))\big] \odot \mathbf{z}_i, \qquad \hat{\mathbf{y}}_{i, p} = D_\omega(\hat{\mathbf{z}}_{i, p}) = \operatorname{softplus}(r_{i,p}) \in \mathbb{R}_{>0}^{44,728},$$
   where:
   - $\mathbf{e}_p \in \mathbb{R}^{128}$ is row $\mathrm{pid}(p)$ of the cached deterministic perturbation table $\mathbf{Z}_{\mathrm{graph}} = \mathrm{GNN}(\cdot) \in \mathbb{R}^{18,216 \times 128}$;
   - $\bar{\mathbf{e}} = \frac{1}{18,216}\sum_{q=0}^{18,215} \mathbf{e}_q$ is the fixed mean across **all `18,216` graph nodes** (`model.py:300-301, 566` — never a batch statistic);
   - $g(\cdot)$ is the 2-layer `gamma_net` MLP (`128 -> 256 -> 512`, `film_activation="none"`, `film_mode="multiplicative_only"`, `film_scale=1.0`, `film_mean_free=True`);
   - For the unperturbed reference (`[-1]`, `pert_mask = 0`), `_slot_embed` yields $\mathbf{e}_0 = \mathbf{0}, m_0 = 0$, and `_interact` computes $\delta \mathbf{z}_{i, -1} = \mathrm{interaction}(\mathbf{z}_i, \mathbf{0})$ without mean-free subtraction (`model.py:560-573, 763-770, 795-804`), giving $\hat{\mathbf{y}}_{i, -1} = D_\omega(\mathbf{z}_i + \delta \mathbf{z}_{i, -1})$ identically to `api.py:519`.
   - Pseudobulk predicted expression for perturbation $p$ is the **mean of the per-cell decoded predictions** across all $n_{\mathrm{ctrl}}$ control cells: $\overline{\hat{\mathbf{y}}}_p = \frac{1}{n_{\mathrm{ctrl}}}\sum_{i=1}^{n_{\mathrm{ctrl}}} \hat{\mathbf{y}}_{i, p}$ (`api.py:507`, never $D_\omega(E_\theta(\bar{\mathbf{x}}))$!).
2. **Exact Inference Determinism (`P5`):**
   Occam inference under this protocol is **100% deterministic**:
   - `model.eval()` and `torch.no_grad()` disable every dropout module (`adj_dropout=0.2`, decoder `dropout=0.3`, encoder `dropout=0.0`);
   - `source_emb = None` and `modality_emb = None` (`source_dropout` is inert; zero discrete cell-line, batch, or source IDs enter `forward`, verified in Gate 5A);
   - The decoder emits the deterministic continuous `softplus` log-expression mean directly (**zero count/Poisson/NB sampling**);
   - No stochastic operation exists anywhere inside Occam's forward pass (the only random seeds in Addendum B are the already-locked `control_barcodes_r16.json` cell draws [`20261601` / `20261602`], the bootstrap draws [`20261601`], and Stage 5B's deterministic `blake2b` split seed `42`).
3. **Single Pure Control-Swap Arm Definition:**
   - **Matched arm (`Occam_matched(t)`):** Feed frozen `seed_42` target line $t$'s `512` locked control cells (`seed = 20261601`).
   - **Mismatched arm (`Occam_mis(s')(t)`):** Hold frozen `seed_42` weights and perturbation indices fixed; replace target line $t$'s `512` control cells with mismatched source $s'$'s `512` locked control cells (`seed = 20261601`).
   - **Within-line disjoint resample arm (`Occam_resample(t)`):** Hold frozen `seed_42` weights and perturbation indices fixed; replace target line $t$'s primary `512` controls with target line $t$'s disjoint `512` resample controls (`seed = 20261602`, `0` overlap).

### §3.2 Exact Reuse of Round 16 Locked Control Barcodes
- **Primary `512`-cell pools (`seed = 20261601`):** Load exact ordered barcode lists verbatim from `round16/results/control_barcodes_r16.json` for `{a549_ko, dld1, cd4t_rest, nadig_hepg2, nadig_jurkat, feng_ipsc}` (all `512` cells each).
- **Disjoint `512`-cell within-line resample pools (`seed = 20261602`, `overlap_with_primary == 0`):** Load exact ordered barcode lists verbatim from `round16/results/control_barcodes_r16.json` for `{a549_ko, dld1, cd4t_rest}` and `round16/results/sensitivity/control_barcodes_r16.json` for `feng_ipsc`.
- **Group 3 (`gse295214_h1`, if approved & qualified):** Reuse the exact `PREREG_I §2.2` (`20261601`) and `Addendum A §8` (`20261602`) ordered barcode lists generated for `STATE_SE`.
- **Local Raw Count Paths:** Read `nadig_hepg2` and `nadig_jurkat` strictly from local `/Users/mac/.cache/r9_p2_raw/nadig2025/GSE264667_{hepg2,jurkat}_raw_singlecell_01.h5ad` (not from the mounted sync share, `P3`), and Group 2 / `feng_ipsc` from local SSD (`round10–12/shared/raw/`).
- **Batching:** Feed cells in exact locked order `0..511`, `batch_size = 256`, `shuffle = False`.

### §3.3 Input Gene Mapping, Duplicate Guard, Normalization, and Mandatory `P4` `pairwise_common_axis` Sensitivity
1. **Native Ensembl Gene ID Extraction (`0` Duplicates Across All 6 Local Datasets):**
   Strip `.N` version suffixes (`id.split(".")[0]`) to match Occam's `44,728`-gene Ensembl vocabulary:
   - `dld1` (10x `features/id`), `feng_ipsc` (Ensembl token before `:`), `nadig_hepg2` (`var["gene_id"]`), `nadig_jurkat` (`var["gene_id"]`): **`6,983 / 6,983`** RHJ axis genes present (`0` absent);
   - `a549_ko` (`ensembl_id`): **`6,950 / 6,983`** RHJ axis genes present (**`33` absent**, zero-filled);
   - `cd4t_rest` (10x `features/id`): **`6,570 / 6,983`** RHJ axis genes present (**`413` absent**, zero-filled);
   - `gse295214_h1` (if approved): 10x `features/id` stripped of `.N`.
2. **Duplicate-ID / Duplicate-Symbol Guard (`DEVIATION_R16_D7` Mirror):**
   If any dataset (e.g., `gse295214_h1`) has duplicate version-stripped Ensembl IDs (or duplicate symbols under a symbol map), **strictly keep the leading occurrence in original raw file column order** and drop subsequent duplicates before scattering into the `44,728` encoder input.
3. **Locked Primary Input Protocol (`rhj_6983_axis` — Matches Official `rhj_fewshot` Training & `k562_seen`/`h1_seen` Evaluation):**
   Compute per-cell library size $L_i = \sum_{g \in \mathcal{G}_{\mathrm{raw\_file}}} x_{i, g}$ over **all raw measured `Gene Expression` features** (`build_splithalf_panel.py:100-150`, `splithalf/CARD.md:32`), compute `log1p(x_raw / L_i * 2000.0)` in `float32` (`0.0` if $L_i = 0$), and scatter **strictly the `6,983` RHJ axis Ensembl IDs** (`gene_axis.csv`, SHA-256 `ef64f2ba7ff7816ed947796066cd8a68014884c6176619d750ef79462c78d4c0`) into their exact columns of the `44,728`-gene encoder tensor, zeroing all `37,745` non-axis columns and any absent axis genes.
4. **Mandatory Pre-Registered Gene-Coverage Sensitivity Arm: `pairwise_common_axis` (`P4`):**
   - *Confound addressed:* Because `cd4t_rest` is missing `413` RHJ axis genes and `a549_ko` is missing `33` (both zero-filled), swapping controls between pools of unequal axis coverage alters **both** input gene coverage/sparsity **and** biological cell-line state, which could artificially degrade or improve a mismatched arm independently of true context use.
   - *Exact specification:* While primary `H-I1_Occam` / `H-I2_Occam` stays on `rhj_6983_axis` as drafted, for **every `(target t, mismatched source s')` pair that involves `cd4t_rest` or `a549_ko`** (as target OR as mismatched source, plus any pair involving `gse295214_h1` if it lacks RHJ axis genes):
     1. Define the pairwise common input axis $\mathcal{A}_{\mathrm{common}}(t, s') = \mathcal{A}_{\mathrm{RHJ}} \cap \mathcal{G}_{\mathrm{raw}}(t) \cap \mathcal{G}_{\mathrm{raw}}(s')$;
     2. Compute per-cell library sizes over each cell's full raw file and apply `target_sum = 2000.0 -> log1p` identically to primary `rhj_6983_axis`;
     3. **Zero in BOTH the matched target input (`t`) AND the mismatched source input (`s'`) every RHJ-axis input gene absent from $\mathcal{A}_{\mathrm{common}}(t, s')$**, then run frozen `Occam` on both equally masked pools and score on the identical `G_eval^R16(t)` output genes.
   - *Pre-Locked Verified Common-Axis Input Gene Counts:*
     - `(a549_ko, cd4t_rest)` in either target/source role: **`6,557 / 6,983`** retained (`426` union-absent genes zeroed in both arms);
     - `(cd4t_rest, s)` for $s \in \{\text{dld1}, \text{feng\_ipsc}, \text{nadig\_hepg2}, \text{nadig\_jurkat}\}$ in either role: **`6,570 / 6,983`** retained (`413` zeroed in both arms);
     - `(a549_ko, s)` for $s \in \{\text{dld1}, \text{feng\_ipsc}, \text{nadig\_hepg2}, \text{nadig\_jurkat}\}$ in either role: **`6,950 / 6,983`** retained (`33` zeroed in both arms);
     - Pairs strictly among `{dld1, feng_ipsc, nadig_hepg2, nadig_jurkat}` already share `6,983 / 6,983` genes (`0` absent $\implies$ identical to primary).
   - *Mandatory Reporting:* Report `pairwise_common_axis` paired `r_orth_cal` ($\Delta^{\mathrm{pair}}_{t,\mathrm{common}}(s')$, `97.5%` & `95.0%` TOST CIs, and classification) **side-by-side next to primary `rhj_6983_axis`** for every pair involving `cd4t_rest` or `a549_ko` (and `gse295214_h1` if qualified). Any change in pair/line verdict between primary and `pairwise_common_axis` MUST be highlighted explicitly as gene-coverage sensitivity.
5. **Secondary Input-Axis Sensitivity (`full_44728_api`):** Also report Occam's `GeneSpace.prepare` default (row-sum over all `44,728`-matched Ensembl genes, populating all matched columns) descriptively.

### §3.4 Locked Output-to-Delta Conversion (`baseline="reconstructed"` Primary vs. `baseline="real"` Sensitivity)
- **Primary (`baseline="reconstructed"`, Occam `api.py:486, 510-521` Default & `STATE` `PREREG_I §2.1.4` Analogue):**
  $$\hat{\Delta}_{\mathrm{pred}}^{\mathrm{Occam}}(p) = \frac{1}{n_{\mathrm{ctrl}}}\sum_{i=1}^{n_{\mathrm{ctrl}}} \hat{\mathbf{y}}_{i, p} - \frac{1}{n_{\mathrm{ctrl}}}\sum_{i=1}^{n_{\mathrm{ctrl}}} \hat{\mathbf{y}}_{i, -1},$$
  where $\hat{\mathbf{y}}_{i, -1} = D_\omega(\mathbf{z}_i + \mathrm{interaction}(\mathbf{z}_i, \mathbf{0}))$ is the exact per-cell zero-perturbation forward (`api.py:519`, `model.py:795-804`), canceling the `~5×` decoder control-reconstruction bias (`F-recon1`, `docs/api.md:96-99`, `model.py:824-825`).
- **Secondary Output-Baseline Sensitivity (`baseline="real"` / `injctrl`):** Subtract mean preprocessed real controls $\frac{1}{n_{\mathrm{ctrl}}}\sum_{i=1}^{n_{\mathrm{ctrl}}} \mathbf{x}_i$ (matching Occam's `official_tables.md` `nocal` convention), reported descriptively without mixing conventions (`docs/api.md:98-99`).

### §3.5 Critical Perturbation-Key Resolution Lock (Two-Tier `pert2id` HGNC Rule)
Occam's `18,216`-gene `master_gene_set_sorted.csv` contains **`317` duplicate HGNC pairs** where both a legacy symbol (e.g., `AARS`, `DARS`, `TAZ`, `H2AFZ`) and the current approved HGNC symbol (`AARS1`, `DARS1`, `TAFAZZIN`, `H2AZ1`) exist as separate nodes, whereas `rhj_fewshot`'s `2,381` trained perturbations (`0` HGNC collisions) trained **only** the legacy node. To prevent routing `18` Group 2 targets onto untrained duplicate nodes:
1. **Tier 1 (Trained `rhj_fewshot` Node Priority):** Given canonical HGNC target $p \in \mathcal{T}_{\mathrm{eval}}^{\mathrm{R16}}(t)$, resolve to the unique key $k$ in `benchmarks/rhj_fewshot/perturbations.csv` (`eligible == 1`) with `HGNCBridge.canonicalize(k) == p`. Resolves **`402/403` `a549_ko`** (14 bridged pairs: `AARS1->AARS`, `CARS1->CARS`, `DARS1->DARS`, `EPRS1->EPRS`, `GARS1->GARS`, `HARS1->HARS`, `LARS1->LARS`, `MARS1->MARS`, `NOPCHAP1->C12orf45`, `PHB1->PHB`, `POLR1G->CD3EAP`, `POLR1H->ZNRD1`, `RARS1->RARS`, `TARS1->TARS`), **`75/75` `dld1`** (3 bridged: `DARS1->DARS`, `NARS1->NARS`, `TAFAZZIN->TAZ`), **`245/245` `cd4t_rest`** (1 bridged: `H2AZ1->H2AFZ`), and **`88/88` `feng_ipsc`** (`0` bridged).
2. **Tier 2 (Direct `pert2id` Fallback):** If absent from `rhj_fewshot`'s `2,381` eligible keys, resolve $p$ directly in `ckpt["pert2id"]` (resolving **`RBM14`** in `a549_ko` $\to$ `403/403` total; descriptive `402`-target trained-only sensitivity also reported).
- **Output Gene Projection:** Gather the `6,983` RHJ Ensembl outputs, canonicalize via frozen R10 `HGNCBridge` (`0` collisions), and select exact ordered `G_eval^R16(t)` columns (`1,837`, `1,832`, `1,830`, `1,575`); also report descriptive `Occam`-native `RHJ_6983 ∩ G_eval^R9` (`6,876`, `6,883`, `6,489`, `5,863` genes).

---

## §4. Metrics, Decision Cascades (`H-I1_Occam`, `H-I2_Occam`), and Paired `STATE_SE` vs. `Occam` Contrasts (Item 4)

### §4.1 Metrics & Decision Rules Copied Verbatim from `PREREG_I.md` + `PREREG_I_ADDENDUM_A.md`
All primary endpoints, parameters, and decision rules are inherited **verbatim**:
- **Primary Endpoint:** `r_orth_cal` (uncentered orthogonalization to per-perturbation `Transfer_Mean` via frozen `round10/src/metrics_r10.py::_orth_pair`, pooled moment correlation on unselected **`all_nontarget`**, Spearman–Brown disattenuation via **`rho_ab_orth`** clipped to `[-1, 1]`).
- **Secondary Endpoints:** `r_dorth_cal` (`PREREG_I §4.1`), practical metrics (raw Pearson $\Delta$, centered DE-50 Pearson, cosine similarity, `M2` cross-mask shift), `pairwise_common_axis` (`§3.3.4`), `full_44728_api` (`§3.3.5`), and `baseline="real"` (`§3.4`).
- **Bootstrap & Multiplicity:** `n_boot = 4,000` joint perturbation bootstraps (`seed = 20261601`, identical resample indices across Occam, `STATE_SE`, and `FCT`); Bonferroni $\alpha_{\mathrm{claim}} = 0.05 / 2 = 0.025$ across `H-I1_Occam` and `H-I2_Occam`:
  - Two-sided **`97.5%` percentile CI (`[1.25%, 98.75%]`)** for superiority/inferiority;
  - Two-sided **`95.0%` percentile CI (`[2.5%, 97.5%]`)** inside open **$(-\delta, +\delta)$ with $\delta = 0.04$** (`valid_frac >= 0.95`) for TOST equivalence.
- **Per-Pair Priority Cascade (`PREREG_I §3.1` + `Addendum A §10`):**
  1. `INCONCLUSIVE` (`valid_frac < 0.95` or non-finite `97.5%` CI; `within_delta = False`);
  2. `EXACT_INVARIANCE` (both delta matrices finite AND $d_{\max} \le 10^{-12}$; counted identically as `EQUIVALENT` via `n_equiv_effective = counts["EQUIVALENT"] + counts["EXACT_INVARIANCE"]`);
  3. `INCONCLUSIVE` (conservative zero-tie precedence: `ci975_lo == 0.0` OR `ci975_hi == 0.0`);
  4. `MATCHED_BETTER` (`ci975_lo > 0.0`, with `within_delta` flag);
  5. `MISMATCHED_BETTER` (`ci975_hi < 0.0`, with `within_delta` flag);
  6. `EQUIVALENT` (`within_delta == True`);
  7. `INCONCLUSIVE` (all remaining cases including exact `±0.04` boundary ties).
- **Line-Level `H-I1_Occam` 8-State Cascade (`PREREG_I §3.1`):** `MISSING_STRATUM` $\to$ `INCONCLUSIVE` $\to$ `LINE_MATCHED_SPECIFIC` $\to$ `LINE_MATCHED_SPECIFIC_SUB_DELTA` (**MUST** follow `Addendum A §9` mandatory wording *"能检测到，但小于事先规定的最小有意义效应 δ=0.04"* jointly with §8 within-line resample, NEVER standalone *"模型用到了对照"*) $\to$ `LINE_EQUIVALENT_NONZERO_SIGNAL` $\to$ `LINE_EQUIVALENT_NULL_SIGNAL` $\to$ `LINE_MISMATCHED_BETTER` $\to$ `MIXED_OR_INCONCLUSIVE` (any `INCONCLUSIVE` reported strictly as `检验力不足`).
- **Within-Line Disjoint Resample Arm (`Addendum A §8`):** Compute $\Delta^{\mathrm{pair}}(\mathrm{resample}) = r(\mathrm{matched}) - r(\mathrm{resample})$ and per-source $\Delta^{\mathrm{excess}}(s') = \Delta^{\mathrm{pair}}(s') - \Delta^{\mathrm{pair}}(\mathrm{resample})$ using locked `seed = 20261602` controls and `Addendum A §8.3` narrative rules.
- **Claim `H-I2_Occam` (`Occam` vs. `FCT`, `PREREG_I §3.2`):** Paired $\Delta^{\mathrm{vsFCT}}_{\mathrm{Occam}, t} = r^{\mathrm{orth\_cal}}(\mathrm{Occam}_{\mathrm{matched}}(t)) - r^{\mathrm{orth\_cal}}(\mathrm{FCT}(t))$ (`97.5%` CI $\to$ `MODEL_SUPERIOR_TO_FCT` / `FCT_SUPERIOR_TO_MODEL` / `INDISTINGUISHABLE`), subject to the `M6.2` Anti-Masking Rule.
- **Group-Level Reporting:** Report Group 2 CRISPRko (`a549_ko`), Group 2 CRISPRi (`{dld1, cd4t_rest}`), optional Group 3 (`gse295214_h1`), and `SENSITIVITY_ONLY` (`feng_ipsc`) strictly separately.

### §4.2 Pre-Registered Paired `STATE_SE` vs. `Occam` Contrasts (`H-I3a` & `H-I3b`)
Using the identical `(N_eval^R16, G_eval^R16)` evaluation matrices, identical locked control pools, and identical `n_boot = 4,000` joint perturbation bootstrap draws (`seed = 20261601`) on `r_orth_cal` (`all_nontarget`; mirrored descriptively on `r_dorth_cal`), comparing `STATE_SE` (`fixed_index_zero`) vs. `Occam` (`rhj_6983_axis`, `baseline="reconstructed"`, plus `pairwise_common_axis` cross-check):
1. **`H-I3a` (Matched-Control Accuracy Contrast):** $\Delta^{\mathrm{acc}}_{t} = r^{\mathrm{orth\_cal}}(\mathrm{Occam}_{\mathrm{matched}}(t)) - r^{\mathrm{orth\_cal}}(\mathrm{STATE\_SE}_{\mathrm{matched}}(t))$, classified via `97.5%` CI (`valid_frac >= 0.95`, zero-tie $\to$ `INDISTINGUISHABLE`) into `OCCAM_SUPERIOR_TO_STATE_SE` / `STATE_SE_SUPERIOR_TO_OCCAM` / `INDISTINGUISHABLE`.
2. **`H-I3b` (Matched-vs-Mismatched Control-Specificity Contrast):** Per mismatched source $s'$ AND per stratum mean (`seen_by_checkpoint`, `unseen_by_checkpoint` separately), $\Delta^{\mathrm{spec}}_{t}(s') = \Delta^{\mathrm{pair}}_{\mathrm{Occam}, t}(s') - \Delta^{\mathrm{pair}}_{\mathrm{STATE\_SE}, t}(s')$ (plus resample-adjusted $\Delta^{\mathrm{spec,excess}}_{t}(s')$), classified via `97.5%` CI into `OCCAM_MORE_CONTROL_SPECIFIC` / `STATE_SE_MORE_CONTROL_SPECIFIC` / `INDISTINGUISHABLE`.

*Scope Note:* `H-I3a` and `H-I3b` are pre-registered comparative secondary contrasts reported per line without modifying Round 16's locked primary `STATE_SE` Bonferroni family.

### §4.3 Pre-Locked Advance Outcome Interpretation Matrix (`STATE_SE` vs. `Occam`)

| Joint Outcome Pattern on Target Line $t$ | Pre-Locked Scientific Interpretation (Fixed Before Inference) |
|---|---|
| **Case 1:** Both `STATE_SE` and `Occam` fail target-control specificity (`LINE_EQUIVALENT_*` or `SUB_DELTA` not exceeding within-line resample) while `FCT` succeeds | **Cross-architecture replication of target-control under-utilization:** Failure to exploit target-line controls generalizes across both Transformer-embedding (`STATE_SE`) and multiplicative-FiLM MLP (`Occam`), ruling out a `STATE`-only artifact (such as `batch_encoder` dominance). |
| **Case 2:** `Occam` achieves `LINE_MATCHED_SPECIFIC` (and exceeds within-line resample loss, verified under `pairwise_common_axis`) while `STATE_SE` does not (`OCCAM_MORE_CONTROL_SPECIFIC`) | **Architecture-specific control sensitivity:** Multiplicative FiLM modulation on per-cell control latents (`Occam`) extracts target-line-specific control signal where `STATE_SE` does not — falsifying any universal claim that all neural virtual-cell models ignore target controls (read jointly with `H-I2_Occam` vs. `FCT`). |
| **Case 3:** `STATE_SE` achieves `LINE_MATCHED_SPECIFIC` while `Occam` does not (`STATE_SE_MORE_CONTROL_SPECIFIC`) | **Scale/representation advantage for `STATE_SE`:** `SE-600M`'s 167M-cell embeddings preserve target-line control specificity that Occam's `DiverseRNA-1.4M` + `6,983`-axis encoder fails to maintain under cross-study shift (check `pairwise_common_axis` for `cd4t_rest`/`a549_ko` coverage effects). |
| **Case 4:** `Occam` wins matched accuracy (`OCCAM_SUPERIOR_TO_STATE_SE`) YET fails `H-I1_Occam` (`LINE_EQUIVALENT_NONZERO_SIGNAL` or `<=` within-line resample) | **Control-agnostic perturbation-prior gain:** Occam's higher `r_orth_cal` stems from stronger target-invariant perturbation representations (GO/STRING/GenePT/ESM2) rather than target-cell-line control conditioning — confirming that matched accuracy alone cannot prove context specificity. |
| **Case 5:** Both `STATE_SE` and `Occam` achieve `LINE_MATCHED_SPECIFIC` AND exceed within-line resample loss (holding under `pairwise_common_axis`) | **Cross-architecture target-control utilization:** Both neural models genuinely utilize target-line controls beyond `Transfer_Mean`; `H-I2` and `H-I3a`/`H-I3b` quantify relative standing against closed-form `FCT`. |
| **Case 6:** Either/both return `MIXED_OR_INCONCLUSIVE` (`dld1` `N=75` risk) or `LINE_MISMATCHED_BETTER` | **Power-limited or OOD-degradation outcome:** `INCONCLUSIVE` is reported strictly as `检验力不足` (never "equivalence ruled out"); `LINE_MISMATCHED_BETTER` indicates OOD degradation under target controls relative to familiar training-line (`seen`) or cleaner (`unseen`) controls (cross-check `pairwise_common_axis`). |

---

## §5. Pre-Locked Two-Stage Adapter Reproduction Gate (Item 5; `P3`)

Both stages below MUST pass **BEFORE ANY Group 2, Group 3, or `feng_ipsc` prediction is scored**.

### §5.1 Stage 5A — Structural Integrity, Pure-Synthetic Control Sensitivity, and MPS-vs-CPU Parity (0 Real Cells)
1. **Checksum & Checkpoint Assertions:** `best_model.pth` SHA-256 == `11d5ba5515e79710b3a952532c60aa29515046117292528c59d452a0937685ff`; `is_virtual_cell_checkpoint is True`; `target_sum == 2000.0`; `output_activation == "softplus"`; `interaction_type == "film"`, `film_mode == "multiplicative_only"`, `film_mean_free is True`, `film_hidden_dim == 256`; `len(gene_names) == 44728`, `len(pert2id) == 18216`; `gene_axis_sha256 == "ef64f2ba7ff7816ed947796066cd8a68014884c6176619d750ef79462c78d4c0"`; `source_emb is None`, `modality_emb is None`, zero `batch_encoder` parameters; `strict=True` load cleanly.
2. **Pure-Synthetic Smoke & Determinism Assertions (`n = 32` synthetic cells, `seed = 20261601`):**
   - Strictly positive finite `softplus` outputs (`dead_gene_frac < 0.01`);
   - Exact forward determinism: two repeated calls under `model.eval()` produce `max_abs == 0.0` (CPU) / `<= 1e-6` (MPS); cross-process `PYTHONHASHSEED` (`1` vs. `99991`) yields `max_abs == 0.0` (CPU) / `<= 1e-6` (MPS);
   - Swapping synthetic `X_ctrl_A != X_ctrl_B` produces `max_abs_delta_diff > 1e-4`;
   - Single-perturbation identity invariant: `forward_composed(..., "decode_seq_mean")` equals `forward_composed(..., "sum")` within `<= 1e-6`.
3. **MPS-vs-CPU `float32` Device Parity Check (and PyTorch Cross-Version Check):**
   - Compare identical synthetic forward outputs on `CPU float32` vs. `MPS float32`: lock **`MPS float32`** iff `max_abs(z_ctrl_mps - z_ctrl_cpu) <= 1e-5` AND `max_abs(delta_mps - delta_cpu) <= 1e-4`; otherwise lock **`CPU float32`** for both Main Line and `A-T3`.
   - If running on Fallback Tier 2 (`torch != 2.4.1`, §6.1), run Occam's `scripts/verify_checkpoint.py` structural checks and verify cross-version numerical agreement within Occam's locked `CROSS_VERSION_ATOL = 1e-4, CROSS_VERSION_RTOL = 1e-3`.

### §5.2 Stage 5B — Quantitative Author-Benchmark Reproduction Gate on Local `nadig_hepg2` (`P3`)
- **Local File Path (Strictly Local Cache, NOT Sync Share):**
  - Read HepG2 strictly from **`/Users/mac/.cache/r9_p2_raw/nadig2025/GSE264667_hepg2_raw_singlecell_01.h5ad`** (controller-verified SHA-256 **`e1ad7c3c5a201c861a207a858aa7e59f5e6ac1955674c415f7de0d1dadadb52e`**; `145,473 cells × 9,624 genes`, `4,976` controls, `6,983/6,983` RHJ axis Ensembl IDs in `var["gene_id"]`, perturbation labels in `obs["gene"]`, batches in `obs["gem_group"]`, dense HDF5 `X` dataset — whereas Occam's repackaged `nadig_2024_hepg2.h5ad` has SHA-256 `f31a85d8...`, CSR `X`, `var["ensembl_gene_id"]`, `obs["gene_name"]`, and `obs["dataset_batch_num"]`).
- **Explicit Audit of Author Split Artifacts (`P3`):**
  - Having read `benchmarks/rhj_fewshot/splithalf/splits_manifest.json`, `cohort.json`, `scripts/build_splithalf_panel.py`, and `src/occam/vcell/panel/splithalf.py`: **neither `splits_manifest.json` nor `cohort.json` stores cell-level split index arrays**. `cohort.json` stores the frozen `611`-perturbation list and per-perturbation cell counts (`min_cells = 10` subset of the `621` held-out HepG2 test perturbations, §2.1); `splits_manifest.json` stores the deterministic `blake2b(SPLIT_SALT | 42 | split_index | group)` recipe (`SPLIT_SALT = "occam.vcell.panel.splithalf/v1"`) and the `20` exact SHA-256 digests (`pert_A_sha256`, `pert_B_sha256`, `ctrl_A_sha256`, `ctrl_B_sha256` for splits `0..4`); the binary row-index arrays themselves live in `splithalf/splits.npz` (`600,792 B`, indexing Occam's `f31a85d8...` file).
- **Exact Gate 5B Protocol:**
  1. From local `/Users/mac/.cache/r9_p2_raw/nadig2025/GSE264667_hepg2_raw_singlecell_01.h5ad`, filter the `42,402` panel rows (`37,426` cells across the `611` `cohort.json` perturbations + `4,976` `non-targeting` controls), normalizing over all `9,624` file genes (`target_sum = 2000.0 -> log1p`) and projecting onto the `6,983`-gene RHJ axis.
  2. Run Occam's exact deterministic `draw_splits(labels, batches, cohort, n_splits=5, base_seed=42)` (`src/occam/vcell/panel/splithalf.py:219-252`, mapping `obs["gene"] -> labels`, `obs["gem_group"] -> batches`).
  3. Compute the `20` SHA-256 digests of `pert_A`, `pert_B`, `ctrl_A`, `ctrl_B` (splits `0..4`) and check them directly against `splits_manifest.json` (`32–87`), verifying `max |A − B| <= 1` per control batch across all `5` splits.
  4. Score `OccamAdapter` (`seed_42`, `rhj_6983_axis`, **`injctrl` / `baseline="real"` `nocal`** convention matching Occam's `official_tables.md`) against author-archived HF `official_comparisons/selected_scores_by_seed_split.csv` (`rhj_fewshot | all | ABL_N_ALL | nocal`, `seed_42`, `611` perts, `5` splits):
     - Author `ABL_N_ALL` `seed_42` `nocal` **`pearson_delta_mean = 0.394178`** (`5`-split range `0.392278–0.397797`);
     - Author `ABL_N_ALL` `seed_42` `nocal` **`mae_mean = 0.030157`**;
     - Author `experimental` replicate `nocal` **`pearson_delta_mean = 0.372174`** (`5`-split range `0.370873–0.373463`).
- **Pre-Locked Dual Pass Criteria (`P3`, ALL 3 Metrics Required; Determined Automatically by Split-Hash Parity):**
  - **Case A — Exact Author Split Parity (`20/20` SHA-256 match against `splits_manifest.json`):**
    1. `abs(adapter_mean_pearson_delta - 0.394178) <= 0.005`;
    2. `abs(adapter_mean_mae - 0.030157) <= 0.002`;
    3. `abs(local_experimental_pearson_delta - 0.372174) <= 0.002`.
  - **Case B — Fresh Deterministic `blake2b` Splits (`< 20/20` hash match solely because local `e1ad7c3c...` cell/batch ordering differs from Occam's repackaged `f31a85d8...`, logged explicitly):**
    1. `abs(adapter_mean_pearson_delta - 0.394178) <= 0.015`;
    2. `abs(adapter_mean_mae - 0.030157) <= 0.005`;
    3. `abs(local_experimental_pearson_delta - 0.372174) <= 0.015`.
- **Failure Protocol:** On any failure in Stage 5A or Stage 5B, **halt immediately prior to any Group 2/3/sensitivity inference** and escalate to Master Controller. No adjudication rule or threshold may be altered.

---

## §6. Isolated Environment Plan, Compute/MPS Stall Rule, and Clean-Room Independent Audit Route (Item 6; `P2`, `P6`)

### §6.1 Isolated `round16_addB/.venv`, Compute Budget, and Pre-Locked 10-Minute MPS Stall Rule (`P2`, `P6`)
- **Strict Environment Isolation (`P2`):**
  - **NEVER install, upgrade, or modify any package in `round16/.venv`** (strictly frozen for Round 16).
  - Build a dedicated **`round16_addB/.venv`** strictly after human approval (`R16_ADDB_HUMAN_APPROVAL.txt`):
    - **Primary target (Occam `pyproject.toml` pins):** `Python 3.8.19` + `torch==2.4.1`, `numpy==1.24.4` (`>=1.22,<1.25`), `scipy==1.9.3`, `scikit-learn==1.1.3`, `pandas==2.0.3`, `anndata==0.9.2`, `h5py==3.9.0`, `PyYAML>=6.0`, `tqdm>=4.66`, `psutil>=5.9` (PyPI verified to serve native `macosx_11_0_arm64`/`macosx_12_0_arm64` `cp38` wheels for `torch==2.4.1`, `scipy==1.9.3`, `scikit-learn==1.1.3`, and `numpy==1.24.4`; `torch-geometric` is unnecessary at inference for `adjacency_map`).
    - **Pre-Specified macOS `arm64` Fallback Chain (Locked Pre-Approval):** Because no `python3.8` interpreter is currently installed on this host (installed: `Python 3.11.15` at `/Users/mac/.local/bin/python3.11`, `3.13.5`, `3.14.7`), if `Python 3.8.19` or legacy `anndata==0.9.2`/`h5py==3.9.0` cannot be built cleanly on macOS `arm64`, fall back in strict order:
      1. *Fallback Tier 1 (Exact PyTorch Pin):* Fresh `round16_addB/.venv` on **`Python 3.11.15`** + **`torch==2.4.1`** (`torch-2.4.1-cp311-none-macosx_11_0_arm64.whl`) + compatible `numpy`/`scipy`/`anndata`/`h5py`/`pandas`/`PyYAML`;
      2. *Fallback Tier 2 (Verified Host Stack):* Fresh `round16_addB/.venv` on **`Python 3.11.15`** + **`torch==2.15.0`**, `numpy==2.4.6`, `scipy==1.17.1`, `anndata==0.12.19`, `h5py==3.16.0`, `pandas==2.3.3`, `PyYAML==6.0.4` (separate venv matching `R16_VENV_FREEZE.txt` versions), strictly gated by Occam's cross-version tolerance **`CROSS_VERSION_ATOL = 1e-4, CROSS_VERSION_RTOL = 1e-3`** (`scripts/verify_checkpoint.py:106-107`) plus Gate 5A & Gate 5B.
    - Freeze exact installed packages to `round16_addB/R16_ADDB_VENV_FREEZE.txt` prior to code freeze; both Main Line and `A-T3` run against the identical frozen `round16_addB/.venv` specification.
- **Memory & Runtime Budget on Apple M3 Pro (`36 GB` Unified RAM):**
  - Pure PyTorch MLP encoder (`94.3M` params) + 4-block SE residual decoder (`128.5M` params) + baked `ConcatKnowledgeModel` & FiLM (`~250M` params/buffers) $\implies$ **`~3.5–5.5 GB` peak RSS** in `float32` (`<< 36 GB`).
  - Caching `z_ctrl = encoder(X_ctrl)` once per control pool (`10` primary/resample pools + `pairwise_common_axis` masked pools $\implies$ `< 20 s` total) and `z_graph` once yields estimated total runtime **`~5–15 minutes` on `MPS float32`** or **`~20–45 minutes` on `CPU float32`** (including Stage 5B + `pairwise_common_axis` + `full_44728_api` sensitivities).
- **Strict No-Concurrent-MPS & 10-Minute Stall Fallback Rule (`P6`):**
  1. **Zero GPU/MPS execution** until Master Controller explicitly confirms `A-T3` has finished its Round 16 run, locks Addendum B, records `R16_ADDB_HUMAN_APPROVAL.txt`, freezes Addendum B code, and issues `R16_ADDB_GO.txt`.
  2. **Pre-Locked 10-Minute MPS Stall Rule (`P6`):** Every execution script emits per-pool/per-line progress logs with wall-clock timestamps and cumulative process CPU time (`os.times()`). **If an MPS run shows no progress for 10 consecutive minutes (no new log line AND zero incremental process CPU time), immediately abort the MPS run and switch BOTH parties (Main Line and `A-T3`) to `CPU float32` for the entire Addendum B run.** This rule is fixed pre-lock and never adjusted post-hoc.

### §6.2 Independent Clean-Room Audit Route (`A-T3` + Master Controller; `P2`)
- **Shared Upstream Model Classes vs. Independent Clean-Room Adapters (`P2`):**
  - Both Main Line (`round16_addB/src/`) and `A-T3` (independent audit directory) import **Occam's own unchanged model classes** directly from the hash-locked pinned source archive of commit `9ddb890e65ca1ec2f81882c78311002a75b52f67` (zero patches to Occam model math).
  - **Each party independently writes its own adapter layer** from this locked protocol: raw HDF5/10x count loading, Ensembl extraction & duplicate handling, `rhj_6983_axis` / `pairwise_common_axis` / `full_44728_api` preprocessing, Two-Tier `pert2id` resolution, `baseline="reconstructed"`/`"real"` delta computation, Stage 5B HepG2 reproduction harness, and `.npz` export.
- **Staged Three-Party Export & Acceptance Tolerances (Mirroring `Addendum A §7` & `§11`):**

| Audit Layer | Parties Compared | Locked Acceptance Criterion |
|---|---|---|
| **1. Control-cell barcodes** | Main Line vs. `A-T3` vs. locked `control_barcodes_r16.json` | **Exact `100%` ordered string equality** (primary `20261601` & resample `20261602`) |
| **2. Occam `z_ctrl` encoder latents** | Main Line vs. `A-T3` (exported `z_ctrl_<pool>.npz`, checked **prior to** deltas) | $\max_{i, d} |z_{\mathrm{ctrl}}^{\mathrm{Main}}[i, d] - z_{\mathrm{ctrl}}^{\mathrm{AT3}}[i, d]| \le 10^{-5}$ |
| **3. Predicted delta arrays** | Main Line vs. `A-T3` (primary `rhj_6983_axis`, `pairwise_common_axis`, `full_44728_api`, `baseline="real"`) | $\max_{p, g} |\hat{\Delta}^{\mathrm{Main}}_{p, g} - \hat{\Delta}^{\mathrm{AT3}}_{p, g}| \le 10^{-4}$ |
| **4. Same-array metric recomputation** | Master Controller recomputing from exported `.npz` arrays (`seed = 20261601`, `n_boot = 4,000`) | $|\mathrm{metric}^{\mathrm{Controller}} - \mathrm{metric}^{\mathrm{Party}}| \le 10^{-10}$ |
| **5. Cross-party metrics & verdicts** | Main Line vs. `A-T3` | $|\mathrm{metric}^{\mathrm{Main}} - \mathrm{metric}^{\mathrm{AT3}}| \le 10^{-5}$ **AND `100%` verdict parity** (`H-I1_Occam`, `H-I2_Occam`, `H-I3a`, `H-I3b`, `pairwise_common_axis`) |

---

## §7. Limitations and Exact Overlap with Occam's `"Transfer"` Baseline & Close Prior Art (Item 7; `P1`)

### §7.1 Pre-Specified Limitations (Disclosed in All Reports)
1. **Unattributed `DiverseRNA-1.4M` Pretraining Blocks:** 7 of 9 encoder pretraining datasets are unidentified (`MODEL_CARD.md:184`). While `a549_ko` (`2026-07/08`) and `dld1` (`2026-07`) strictly post-date `arXiv:2605.31562`, `cd4t_rest` (`2026-02-06`) and `gse295214_h1` (`2026-05-05`) cannot be guaranteed absent from unsupervised encoder pretraining (supervised perturbation absence is 100% certain for all Group 2/3 lines).
2. **Input Gene-Axis Coverage Asymmetry (`cd4t_rest` `413` absent, `a549_ko` `33` absent):** Directly controlled via the pre-registered `pairwise_common_axis` sensitivity arm (§3.3.4, `P4`), reported alongside primary `rhj_6983_axis` for every affected pair.
3. **`a549_ko` Modality Shift & `RBM14`:** `a549_ko` is CRISPRko whereas `rhj_fewshot` is CRISPRi, and `1/403` target (`RBM14`) was not in `rhj_2381` (addressed via the `402`-target sensitivity).
4. **Single Training Seed (`seed_42`):** Single-seed checkpoint (`2.55 GB` vs. `12.75 GB` for 5 seeds); inter-seed variance across `42–46` is not estimated here.
5. **Decoder Control-Reconstruction Bias (`F-recon1`) & Library-Scale Convention:** Occam predicts in `log1p(norm_2000)` vs. Round 16 `delta_obs` in `log1p(norm_1e4)`; both `baseline="reconstructed"` (primary, canceling `F-recon1`'s `~5×` bias) and `baseline="real"` (sensitivity) are reported transparently.
6. **Retrospective Epistemic Scope & Statistical Power:** Group 2 cohorts were previously observed in Rounds 10–12 (`预注册的回顾性模型审计`); `dld1` (`N = 75`) has low/borderline TOST power (`PREREG_I §3.3`), and `cd4t_rest` has a wide relative $\delta = 0.04$ (`71.9%` of `FCT`).

### §7.2 Exact Overlap and Scientific Boundary vs. Occam's `"Transfer"` Baseline, `CELLAUDIT`, and `Molina & Zhang 2026` (`P1`)
Round 16 + Addendum B builds directly on three closely related precedents and claims **no novelty for the input-replacement intervention paradigm itself**:
1. **Direct Methodological Precedent — `CELLAUDIT` (`Li et al. 2026`, `arXiv:2609.27234v1`):** Already established frozen-checkpoint matched input replacement (including compound, dose, and within-plate-family control-profile replacement plus disjoint control-well refits) to test whether fitted predictions depend on an input and gain predictive accuracy from it (`"predictive generalization need not imply generalization of input-use claims"`). Round 16 + Addendum B applies that established frozen-checkpoint input-replacement principle to a distinct domain and hypothesis class: **cross-cell-line** mismatched control-cell swapping (`seen_by_checkpoint` vs. `unseen_by_checkpoint` + within-line disjoint resampling) on open-weight CRISPR virtual-cell models (`Arc STATE_SE` and `Occam`) under `Transfer_Mean`-orthogonalized `r_orth_cal` and TOST equivalence ($\delta = 0.04$) with multiplicative `FCT` as a positive control.
2. **Closest Observational Precedent — `Molina & Zhang 2026` (`bioRxiv 10.64898/2026.07.24.740459`, `xinyizhanglab/perturbation-decomposition`):** Decomposes measured CRISPRi responses across 4 Replogle/Nadig lines into global, cell-line, perturbation, and interaction components ($\delta_{c,p} = \mu + \alpha_c + \beta_p + \gamma_{c,p}$), concludes observationally that `"the perturbation and cell-line specific components are high-dimensional and cannot be recovered from control expression"`, and evaluates `STATE` (`arc-state==0.10.2`) source-only transfer via component recovery (full bioRxiv text `未能读取`, treated as a potentially substantial overlap; full analysis & reconciliation with Round 12 `K1` `FCT` in [`PRIOR_ART_I_ADD2_v2.md §2`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16_addB/PRIOR_ART_I_ADD2_v2.md)). Round 16 + Addendum B complements that observational decomposition with a causal input-replacement intervention on frozen `STATE_SE` and `Occam` checkpoints across post-2025 screens (`a549_ko`, `dld1`, `cd4t_rest`) using `FCT` to test whether neural models capture even the baseline-expression-gated multiplicative residual recoverable from target controls.
3. **Occam's Own `"Transfer"` Baseline (`Bendidi et al. 2026`, `arXiv:2611.21959`):**
   - Verbatim (`Appendix E.2`): *"Transfer. Add the same perturbation’s measured effect from other training cell lines to each sampled B control from the evaluation cell line. The source effects are weighted by their cell counts. This baseline uses perturbation identity but does not learn how the effect changes between cell lines."*
   - Verbatim (`Appendix A`): *"Transferring measured effects from other training cell lines gives higher Pearson Δ than Occam in held-out cell lines, and no model leads on every metric."*
   - *Exact boundary:* While Occam's paper already shows additive `"Transfer"` beats Occam on raw `Pearson Δ` in held-out lines (`h1_seen`, `k562_seen`), it never swaps Occam's own `control_expression` across mismatched cell lines, never removes `"Transfer"` via `r_orth_cal` (where `"Transfer"` equals `0.0` by construction), never runs TOST $\delta=0.04$ or within-line disjoint resampling, and never runs a paired `STATE_SE` vs. `Occam` vs. multiplicative `FCT` comparison on post-2025 Group 2 screens.
