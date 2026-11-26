# Round 16 Pre-Lock Addendum B Draft — Secondary Open-Weight Model: `Occam` (`PREREG_I_ADDENDUM_B_DRAFT.md`)

> **STATUS: PRE-LOCK DRAFT FOR MASTER CONTROLLER (`Claude`) AUDIT & EXPLICIT HUMAN LICENSE/DOWNLOAD APPROVAL (`R16_ADDB_HUMAN_APPROVAL.txt`).**
> **Parent Protocols (Strictly Unmodified):**
> - [`round16/prereg/PREREG_I.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md) (`v5-FINAL-PRE-LOCK`, commit `c1fc7c6d24dd128bc5e2d969d8c78224a2cbc118`, SHA-256 `43d242e460fd64db2e4af1e628e886cf1d94fb4ccd1e95925928a99cb96ee4a8`);
> - [`round16/prereg/PREREG_I_ADDENDUM_A.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I_ADDENDUM_A.md) (Items 1–11, SHA-256 `bdb16818c6e2499d1eb10f8690f3905dd19ce1d7af4c423c9f5afe5455ed3a3b`).
> **Companion Prior-Art Refresh:** [`round16_addB/PRIOR_ART_I_ADD2.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16_addB/PRIOR_ART_I_ADD2.md).
>
> **Pre-Download, Blinding, & Freeze Compliance Attestation:**
> 1. Drafted under `PREREG_I.md §0.4` (`"Any addition of a second open-weight control-conditioned perturbation model (e.g., Occam) requires separate human download approval and a standalone Master Controller addendum"`) and Master Controller directive `2026-11-26 16:18:49Z`.
> 2. **Zero Occam weight files (`.pth`), feature tables (`.pkl`), binary panels (`.npz`), or `GSE295214` files have been downloaded.** Zero Occam code has been cloned or executed. Zero GPU/MPS compute was run while drafting.
> 3. **Zero modifications** were made to any file inside `round16/` (except appending the `R16_ADDB_DRAFT_READY` pointer line to `round16/status/MAIN.txt`). **Zero reads** of `round16/audit/`, `A-T3` logs, or `A-T3.txt` occurred. **Zero git commands** were executed.
> 4. All exact gene-overlap and target-count audits in §2–§3 were computed strictly from frozen local target/gene masks and Occam's tiny public text CSVs (`benchmarks/rhj_fewshot/{gene_axis.csv,perturbations.csv}`, `data/graphs/master_gene_set_sorted.csv`), with zero Occam inference and zero Round 16 prediction metrics consulted.
> 5. No priority or outreach claims are made anywhere in this document.

- **Secondary Audited Model:** **Occam** (Bendidi, Kenyon-Dean, Denton, Genovesio, & Noutahi, `arXiv:2611.21959`, 2026-11-16) — `rhj_fewshot` all-source virtual-cell predictor (`ABL_N_ALL`, `seed_42`).
- **Scientific Purpose:** Extend Round 16's matched-vs-mismatched target-control audit (`PREREG_I.md`) beyond `Arc STATE` (`STATE_SE` / `STATE_HVG`) to an architecturally independent open-weight perturbation model (MLP denoising autoencoder + multiplicative mean-free FiLM + GO/STRING/GenePT/ESM2 perturbation encoder; zero discrete cell-line or batch IDs), enabling both (a) standalone `H-I1_Occam` and `H-I2_Occam` adjudication on unseen Group 2 targets (`{a549_ko, dld1, cd4t_rest}` and optional Group 3 `gse295214_h1`), and (b) a pre-registered paired **`STATE_SE` vs. `Occam`** head-to-head contrast on identical `(targets, genes, control cells)`.
- **Locked RNG Seeds (Identical to Parent Protocol):** Primary `20261601`; within-line disjoint resample `20261602`.

---

## §1. Exact Artifacts to Hash at Download, Repository Pins, and License Disclosure (Item 1)

### §1.1 Pinned Upstream Revisions
- **Source Code Repository:** `https://github.com/IhabBendidi/Occam`
  - **Pinned Commit SHA:** **`9ddb890e65ca1ec2f81882c78311002a75b52f67`** (verified `HEAD` at `2026-11-26`).
  - **Code License:** **`MIT License`** (`LICENSE`, `Copyright (c) 2026 Occam contributors`, granting rights over `"this software and associated documentation files (the 'Software')"`).
- **Model & Benchmark Hub Repository:** `https://huggingface.co/bendidiihab/Occam` (`type: models`, `gated: false`, `tags: ["arxiv:2611.21959", "region:us"]`)
  - **Pinned Hub Revision SHA:** **`acd080b788d028911807b87e11a56498f984f279`** (`lastModified: 2026-11-18T12:25:03.000Z`).

### §1.2 Exact Weights & Assets License Situation (Mandatory Human Approval Required)
Contrary to informal shorthand describing the entire Occam release as MIT, primary inspection of HuggingFace API metadata, `MODEL_CARD.md`, and `assets/README.md` reveals **three distinct licensing gaps on weights and derived embedding tables**:
1. **HuggingFace Hub metadata (`bendidiihab/Occam`):** `cardData` is `null` — no structured `license` field is declared on the Hub repository.
2. **Model card (`IhabBendidi/Occam/MODEL_CARD.md:2-3`):** YAML frontmatter explicitly declares:
   ```yaml
   license: other
   license_name: see-repository
   ```
   while the repository root ships only the `MIT License` covering `"Software"`, with **no standalone weight-license text**.
3. **Derived feature tables (`bendidiihab/Occam/assets/README.md:26-33`):** Verbatim provenance and licensing statement:
   > *"GenePT: the GenePT v2 release (133,736 gene-symbol keys, text-embedding-3-large, 3072-d) ... ESM2: precomputed ESM2 (1280-d, esm2_t33_650M) per-gene protein embeddings from the Hugging Face dataset `genbio-ai/foundation-models-perturbation` @ `d67ba210` (`ESM2_(D=1280).h5ad`) ... Licensing of both sources is not recorded here; the Hub repo is private (see the README)."*

**Governance Requirement:** Because weights and baked/external GenePT + ESM2 tables lack an explicit standalone open weight license and exceed `100 MB`, **NO file below may be downloaded until explicit human license/research-use approval is recorded in `round16_addB/R16_ADDB_HUMAN_APPROVAL.txt` (analogous to `round16/status/R16_HUMAN_APPROVAL.txt`) AND Master Controller locks Addendum B.**

### §1.3 Exact Artifacts to Download and SHA-256 Hash-Lock
All 4 binary artifacts below reside under `bendidiihab/Occam` @ `acd080b788d028911807b87e11a56498f984f279` and total **`2,929,553,709` bytes (`2.93 GB` / `2.73 GiB`)**. At download time, local SHA-256 checksums MUST match the upstream HuggingFace LFS SHA-256 digests listed below bit-for-bit (`round16_addB/R16_ADDB_CHECKSUMS.sha256`):

| # | Hub Relative Path (`bendidiihab/Occam` @ `acd080b7`) | Exact Size (Bytes) | Pre-Locked Upstream HF LFS SHA-256 | Role & Runtime Requirement |
|---|---|---|---|---|
| 1 | `checkpoints/occam_seeds/rhj_fewshot/all_sources/seed_42/best_model.pth` | **`2,550,066,734`** (`2.55 GB`) | `11d5ba5515e79710b3a952532c60aa29515046117292528c59d452a0937685ff` | **REQUIRED — Primary model checkpoint (`ABL_N_ALL`, `seed_42`).** Verified in `src/occam/vcell/pipeline.py::from_vcell_checkpoint:903-960` and `checkpoints/occam_seeds/rhj_fewshot/all_sources/CARD.md` to be **100% self-contained for inference**: encoder, decoder, FiLM, `gene_names` (`44,728`), `pert2id` (`18,216`), adjacency buffers, and energy-scaled GenePT/ESM2 tables are baked into `gnn_state_dict`. |
| 2 | `assets/genept_full_occam_rhj_fewshot.pkl` | **`220,037,958`** (`220.0 MB`) | `a93e1131beb9722228e32d847ebc6f19b0e6d77c3a0a8374b0552320624181c8` | **Training provenance / optional full-bundle verification** (`17,843/18,216` keys, 3072-d GenePT v2; `resolved_config.yaml:97`; NOT re-read at inference because scaled vectors are baked into #1). |
| 3 | `assets/esm2_t33_650M_UR50D_occam_rhj_fewshot.pkl` | **`93,705,169`** (`93.7 MB`) | `cb551710feaca63c77f4e0e008064a3f7a9e645a8900c45c082ab81867ed6a16` | **Training provenance / optional full-bundle verification** (`18,149/18,216` keys, 1280-d ESM2; `resolved_config.yaml:98`; NOT re-read at inference). |
| 4 | `assets/genept_occam_rhj_fewshot.pkl` | **`65,743,848`** (`65.7 MB`) | `595068eacae03f6256003953558f67bbf0e1f03e86422f0801f40debb32df668` | **Provenance audit only** (`2,376` keys; per `assets/README.md:18`, used by the PerturbDiff competitor, `max_abs_diff = 0.0` vs. #2 on shared keys; NOT read by Occam). |

- **Download Scope Option for Human Approval:** Because `best_model.pth` (#1, `2.55 GB`) is strictly self-contained at inference, human approval may explicitly choose either **(a) Minimal Inference Bundle** (Artifact #1 only, `2.55 GB`, plus tiny text sidecars below) or **(b) Complete `rhj_fewshot` Artifact Bundle** (Artifacts #1–#4, `2.93 GB`, plus tiny text sidecars below). Whichever scope is approved is frozen in `R16_ADDB_HUMAN_APPROVAL.txt` and `R16_ADDB_CHECKSUMS.sha256` prior to execution.
- **Required Tiny Text/JSON/CSV Provenance Sidecars (`< 400 KB` total, non-LFS, verified against HF git blob OIDs):**
  - `checkpoints/occam_seeds/rhj_fewshot/all_sources/seed_42/resolved_config.yaml` (`4,435 B`, git OID `51cff88a31fbcdcce6c8ba4b87983c0d6783596f`);
  - `benchmarks/rhj_fewshot/gene_axis.csv` (`156,021 B`, SHA-256 `ef64f2ba7ff7816ed947796066cd8a68014884c6176619d750ef79462c78d4c0` — exact `6,983` RHJ Ensembl IDs + symbols);
  - `benchmarks/rhj_fewshot/perturbations.csv` (`68,294 B` — exact `2,381` eligible trained perturbations + `pert2id`);
  - `benchmarks/rhj_fewshot/splithalf/cohort.json` (`101,500 B` — exact `611` held-out HepG2 test perturbations for §5 Gate 5B);
  - `benchmarks/rhj_fewshot/splithalf/splits_manifest.json` (`4,894 B`);
  - `assets/*.pkl.json` (`30,037 B` total) if Artifacts #2–#4 are included.
- **Explicitly Excluded Downloads (Never Downloaded in Addendum B):** `seed_43`–`seed_46` (`10.2 GB`), `occam_small_encoder_decoder.pth` (`893 MB`), `benchmarks/rhj_fewshot/splithalf/panel.npz` (`171 MB`, rebuilt locally from existing `GSE264667_hepg2_raw_singlecell_01.h5ad` in §5), and all competitor checkpoints.

---

## §2. Evaluation Targets, Two-Stratum Mismatched Control Pools, and Pretraining Disclosures (Item 2)

### §2.1 Training Scope of Occam `rhj_fewshot/all_sources/seed_42` (`ABL_N_ALL`)
Verified directly from `resolved_config.yaml:44-54` and `benchmarks/rhj_fewshot/CARD.md:12-20`:
- **Supervised perturbation training lines (3 lines):** `Replogle RPE1` (`2,381` eligible perts), `Nadig Jurkat` (`2,381` eligible perts), and `Nadig HepG2` (`1,657` train / `103` val / `621` held-out test perts), all projected onto the **`6,983`-gene Ensembl RHJ axis**.
- **Consequence for Round 16 Cohorts:**
  - **Group 1 (`nadig_hepg2`, `nadig_jurkat`) is EXCLUDED as an evaluation target for `H-I1_Occam` and `H-I2_Occam`** because both lines were directly observed during supervised perturbation training of `seed_42`. Both lines enter Addendum B strictly in two non-target roles: (i) as the **`seen_by_checkpoint` mismatched-control stratum** (where each line is **100% `1/1` seen** by the single `seed_42` checkpoint, sharper than `STATE_SE`'s `3/4`-fold ensemble overlap), and (ii) `nadig_hepg2`'s `611` held-out test perturbations serve strictly for the pre-evaluation **Adapter Reproduction Gate (§5B)**.
  - **Primary Evaluated Targets are Group 2 ONLY (`{a549_ko, dld1, cd4t_rest}`)**, plus **Group 3 (`gse295214_h1`)** strictly contingent on human download approval and passing `PREREG_I §2.4`'s counts-only gate, plus **`feng_ipsc`** as a **`SENSITIVITY_ONLY`** target (identical to `PREREG_I §0.2` / `Q4`).

### §2.2 Mandatory `DiverseRNA-1.4M` Encoder Pretraining Disclosure
- Occam's frozen encoder (`44,728 -> 2048 -> 1024 -> 512`) was pretrained on `DiverseRNA` (`1,434,299` profiles across 9 datasets; Kenyon-Dean et al., `arXiv:2605.31562`), with the released softplus+NB decoder finalized `2026-09-16` (`MODEL_CARD.md:23`).
- **Critical Provenance Caveat:** Occam's own `MODEL_CARD.md:183-186` states verbatim:
  > *"Pretraining-corpus caveats: 43.6% of rows are non-integer, 2.4% are bulk RNA-seq, and 7 of 9 constituent datasets are unidentified. Because of the last point, no claim in this card can assert that any evaluation cell line is an unseen context."*
  And `benchmarks/rhj_fewshot/CARD.md:103-104` reiterates: *"7 of 9 DiverseRNA blocks are unattributed, so it cannot be asserted that any of these three cell lines is absent from pretraining."*
- **Date-Based Isolation Audit for Addendum B Targets:**
  - **`a549_ko` (`GSE345058`, GEO `2026-08-26`, ENA `2026-07-03`)** and **`dld1` (`GSE337988`, GEO `2026-07-10`)** strictly post-date `DiverseRNA` (`arXiv:2605.31562`, May 2026) and Occam's frozen encoder.
  - **`cd4t_rest` (`GSE314342`, GEO `2026-02-06`)** strictly post-dates `SE-600M` (`2025-06/07`), **but precedes `arXiv:2605.31562` by ~3.5 months**. Because 7 of 9 `DiverseRNA` blocks are unattributed (`未能确定`), **encoder-pretraining absence for `cd4t_rest` cannot be formally guaranteed**, even though `cd4t_rest` was 100% unseen in Occam's supervised perturbation training (`rhj_fewshot`) and uses a `probe-based` assay (`18,130` genes). This disclosure MUST accompany every `cd4t_rest` Occam table and narrative.
  - **`gse295214_h1` (`GSE295214`, GEO `2026-05-05`, ENA `2026-06-19`)** overlaps the May 2026 `DiverseRNA` publication month (supervised unseen; encoder pretraining absence `未能确定`, alongside `H1`/Replogle-Weissman disclosures from `PREREG_I §2.4.5`).
  - **`feng_ipsc` (`ENA PRJEB81502`, `2024-11-08`)** predates BOTH `SE-600M` and `DiverseRNA` by `> 1.5` years $\implies$ kept strictly **`SENSITIVITY_ONLY`** as a target line, while **retained as an `unseen_by_checkpoint` mismatched-control source** (never in `rhj_fewshot` supervised training).

### §2.3 Locked Target Cohort, Stratum Assignments, and Verified Exact `(N_eval^R16, |G_eval^R16|)` Table
Every evaluated target count and gene-universe count below was verified directly against `round16` frozen reproducibility masks, `benchmarks/rhj_fewshot/gene_axis.csv` (`6,983` genes), `benchmarks/rhj_fewshot/perturbations.csv` (`2,381` eligible trained targets), and `master_gene_set_sorted.csv` (`18,216` perturbation graph nodes):

| Scope / Group | Target Line $t$ | Modality & Lineage | Mismatched Stratum A: `seen_by_checkpoint` ($n=2$, strictly `1/1` seen) | Mismatched Stratum B: `unseen_by_checkpoint` | Shared `R16` Locked Targets $N_{\mathrm{eval}}^{\mathrm{R16}}$ (`rhj_2381` trained + graph-extrap.) | Shared `R16` Locked Genes $\|G_{\mathrm{eval}}^{\mathrm{R16}}\|$ (`%` inside `RHJ_6983`) | Native Full `RHJ_6983 ∩ G_R9` Sensitivity Genes | Mandatory Disclosures |
|---|---|---|---|---|---|---|---|---|
| **Group 2 Primary — CRISPRko ($K_{2,\mathrm{ko}}=1$)** | `a549_ko` | Lung (**CRISPRko**) | `{nadig_hepg2, nadig_jurkat}` | `{dld1, cd4t_rest, feng_ipsc}` ($n=3$) | **`403`** (`403` hv): **`402`** trained in `rhj_2381` + **`1`** (`RBM14`) graph-extrapolated | **`1,837`** (**100.0%**, `0` missing, `0` collisions) | `6,876` | CRISPRko-vs-CRISPRi modality shift (`I11d`); `RBM14` untrained node (`402`-target sensitivity reported); post-dates DiverseRNA |
| **Group 2 Primary — CRISPRi ($K_{2,\mathrm{i}}=2$)** | `dld1` | Colorectal (CRISPRi) | `{nadig_hepg2, nadig_jurkat}` | `{a549_ko, cd4t_rest, feng_ipsc}` ($n=3$) | **`75`** (`75` hv): **`75/75`** trained in `rhj_2381` | **`1,832`** (**100.0%**, `0` missing, `0` collisions) | `6,883` | Low/borderline TOST power at `N = 75` (`M4`); post-dates DiverseRNA |
| **Group 2 Primary — CRISPRi ($K_{2,\mathrm{i}}=2$)** | `cd4t_rest` | Primary CD4$^+$ T (CRISPRi) | `{nadig_hepg2, nadig_jurkat}` | `{a549_ko, dld1, feng_ipsc}` ($n=3$) | **`245`** (`230` hv): **`245/245`** trained in `rhj_2381` | **`1,830`** (**100.0%**, `0` missing, `0` collisions) | `6,489` | `413/6,983` RHJ input genes absent (zero-filled); `2026-02-06` precedes DiverseRNA `2026-05` (`7/9` unattributed); fragile relative $\delta$ (`M6.1`) |
| **Group 3 Prospective Holdout ($K_3 \in \{0,1\}$)** | `gse295214_h1` *(ONLY if approved & §2.4 gate passes)* | `H1 hESC` (CRISPRi `10x 5'`) | `{nadig_hepg2, nadig_jurkat}` | `{a549_ko, dld1, cd4t_rest, feng_ipsc}` ($n=4$) | Determined post-approval by `PREREG_I §2.4` $\cap \mathrm{canon}(P_{\mathrm{Occam}})$ (`>= 50` required) | `G_eval^R16(gse295214_h1) ∩ RHJ_6983` | `RHJ_6983 ∩ G_R9` | Reported strictly separately; shared `H1` lineage (`vcc2025_h1`), Replogle/Weissman provenance, DiverseRNA May 2026 timing |
| **Sensitivity Only (`Q4`)** | `feng_ipsc` | Human iPSC pool (CRISPRi) | `{nadig_hepg2, nadig_jurkat}` | `{a549_ko, dld1, cd4t_rest}` ($n=3$) | **`88`** (`88` hv): **`88/88`** trained in `rhj_2381` | **`1,575`** (**100.0%**, `0` missing, `0` collisions) | `5,863` | `SENSITIVITY_ONLY` target (`ENA 2024-11-08` predates SE-600M & DiverseRNA); **retained as `unseen_by_checkpoint` mismatched source** |

*Critical Parity Guarantee:* On all three Group 2 primary targets (`a549_ko`, `dld1`, `cd4t_rest`) and on `feng_ipsc`, **100.0% of `STATE_SE`'s locked targets (`N_eval^R16`) and 100.0% of `STATE_SE`'s locked evaluation genes (`G_eval^R16`) are covered by Occam with zero dropped targets and zero dropped genes**. Primary evaluation of Occam and the paired `STATE_SE` vs. `Occam` contrast therefore operate on the **exact identical `(N_eval^R16, G_eval^R16)` evaluation matrices** and the **exact identical observed pseudobulk / split-half / `Transfer_Mean` / `FCT` arrays** locked in Round 16.

---

## §3. Exact Occam Inputs, Mismatched-Control Arms, Vocabulary & Perturbation Mapping, Normalization, and Control-Cell Reuse (Item 3)

### §3.1 Exact Context Inputs and Single Pure Control-Swap Arm Definition
Verified directly from `resolved_config.yaml:99-105`, `src/occam/vcell/model.py`, `src/occam/vcell/pipeline.py`, and `arXiv:2611.21959 Appendix C.1`:
- **Model forward signature:** `VirtualCellModel.forward(control_expression, pert_indices)` receives **ONLY** (i) single-cell preprocessed log-expression vectors `control_expression ∈ R^{n_ctrl × 44,728}` and (ii) perturbation node indices `pert_indices`.
- **Latent modulation (`Appendix C.1 Eq. 13`, `resolved_config.yaml:99-105`):**
  $$\mathbf{z}_{\mathrm{ctrl}} = E_\theta(\bar{\mathbf{x}}_{\mathrm{ctrl}}) \in \mathbb{R}^{512}, \qquad \hat{\mathbf{z}} = \mathbf{z}_{\mathrm{ctrl}} + \big[\sigma(g(\mathbf{e}_p)) - \sigma(g(\bar{\mathbf{e}}))\big] \odot \mathbf{z}_{\mathrm{ctrl}}, \qquad \hat{\mathbf{y}} = D_\omega(\hat{\mathbf{z}}) = \operatorname{softplus}(r) \in \mathbb{R}_{>0}^{44,728},$$
  with `interaction_type = "film"`, `film_mode = "multiplicative_only"`, `film_scale = 1.0`, `film_hidden_dim = 256`, `film_activation = "none"`, `film_mean_free = true`, and 4-block SE residual `HybridNormGlobalSkipDecoder` (`dropout = 0.3` disabled in `eval()` mode).
- **Zero Discrete Context/Batch IDs:** Occam has **NO cell-line ID embedding, NO batch ID embedding, and NO assay/source embedding** (`source_emb = None`, `modality_emb = None`, verified programmatically in Gate 5A §5.1). Therefore, unlike `STATE_SE` (which required `fixed_index_zero` plus a 2×2 `G4` batch decomposition), Occam's matched-vs-mismatched audit requires a **single pure control-swap arm per source $s'$**:
  - **Matched arm (`Occam_matched(t)`):** Feed frozen `seed_42` target line $t$'s own `512` locked control cells (`seed = 20261601`) and target perturbation list $\mathcal{T}_{\mathrm{eval}}^{\mathrm{R16}}(t)$.
  - **Mismatched arm (`Occam_mis(s')(t)`):** Hold frozen `seed_42` weights and perturbation indices $\mathcal{T}_{\mathrm{eval}}^{\mathrm{R16}}(t)$ completely fixed; replace target line $t$'s `512` control cells with mismatched source $s'$'s `512` locked control cells (`seed = 20261601`).
  - **Within-line disjoint resample arm (`Occam_resample(t)`):** Hold frozen `seed_42` weights and perturbation indices fixed; replace target line $t$'s primary `512` control cells (`seed = 20261601`) with target line $t$'s disjoint `512` resample control cells (`seed = 20261602`, `0` overlap with primary).

### §3.2 Exact Reuse of Round 16 Locked Control Barcodes
To eliminate any control-sampling divergence between `STATE_SE` and `Occam`, Addendum B **reuses 100% identically** the ordered `512`-cell control barcode lists already locked by Round 16 (`seed = 20261601` primary; `seed = 20261602` disjoint resample):
- **Primary `512`-cell pools (`seed = 20261601`):** Load exact ordered barcode lists from `round16/results/control_barcodes_r16.json` for `{a549_ko, dld1, cd4t_rest, nadig_hepg2, nadig_jurkat, feng_ipsc}` (verified: `512` cells each, drawn from `45,320`, `4,462`, `89,759`, `4,976`, `12,013`, `55,050` QC-passing NTCs respectively).
- **Disjoint `512`-cell within-line resample pools (`seed = 20261602`, `overlap_with_primary == 0`):** Load exact ordered barcode lists from `round16/results/control_barcodes_r16.json` for `{a549_ko, dld1, cd4t_rest}` and `round16/results/sensitivity/control_barcodes_r16.json` for `feng_ipsc`.
- **Group 3 (`gse295214_h1`, if approved & qualified):** Reuse the exact same `PREREG_I §2.2` (`seed = 20261601`) and `Addendum A §8` (`seed = 20261602`) ordered barcode lists generated for `STATE_SE`.
- **Execution order:** Because Occam is strictly single-cell independent (`eval()` mode, `torch.no_grad()`, no cross-cell attention or batch norm), control cells are passed in exact locked index order `0..511` with `batch_size = 256`, `shuffle = False`.

### §3.3 Input Gene Vocabulary Mapping, Duplicate-ID Rule (`DEVIATION_R16_D7` Mirror), and Normalization
1. **Native Ensembl Gene ID Extraction (Verified 0 Duplicate Ensembl IDs Across All 6 Local Datasets):**
   Occam's `44,728`-gene input vocabulary (`ckpt["gene_names"]`) consists of Ensembl gene IDs (`44,560` `ENSG...` IDs, versionless, `docs/api.md:122-124`). Every local raw dataset carries native Ensembl IDs extracted deterministically as follows (stripping any trailing `.N` version suffix via `id.split(".")[0]`):
   - `a549_ko`: `ensembl_id` column (`29,892` raw, `0` duplicates, **`6,950 / 6,983`** RHJ axis genes present, `33` absent zero-filled);
   - `dld1`: 10x `matrix/features/id` (`36,601` raw, `0` duplicates, **`6,983 / 6,983`** RHJ axis genes present);
   - `cd4t_rest`: 10x `matrix/features/id` (`18,130` raw, `0` duplicates, **`6,570 / 6,983`** RHJ axis genes present; **`413` absent genes zero-filled** — explicitly disclosed alongside `h1_seen`'s `448` absent genes in `docs/h1_seen.md:37`);
   - `feng_ipsc`: Ensembl prefix before `:` (`36,518` raw, `0` duplicates, **`6,983 / 6,983`** RHJ axis genes present);
   - `nadig_hepg2`: `var["gene_id"]` (`9,624` raw, `0` duplicates, **`6,983 / 6,983`** RHJ axis genes present);
   - `nadig_jurkat`: `var["gene_id"]` (`8,882` raw, `0` duplicates, **`6,983 / 6,983`** RHJ axis genes present);
   - `gse295214_h1` (if approved): 10x `matrix/features/id` stripped of `.N` version suffixes.
2. **Duplicate-ID / Duplicate-Symbol Guard (Mirroring `DEVIATION_R16_D7` & `Addendum A §11.3`):**
   If any future dataset (e.g., `gse295214_h1`) contains duplicate version-stripped Ensembl IDs (or duplicate raw symbols if a symbol map is ever invoked), **strictly keep the first occurrence in original raw file column order** and drop subsequent duplicates before scattering into the `44,728` vocabulary — ensuring zero ambiguous collisions in `GeneSpace.align` and `100%` deterministic column mapping between Main Line and `A-T3`.
3. **Locked Input Axis and Normalization (Primary vs. Secondary Sensitivity):**
   - **Primary Input Protocol (`rhj_6983_axis` — Matches Occam Official `rhj_fewshot` Training & `k562_seen`/`h1_seen` Evaluation):**
     1. Compute each control cell's library size $L_i = \sum_{g \in \mathcal{G}_{\mathrm{raw\_file}}} x_{i, g}$ over **all raw measured `Gene Expression` features** in that dataset's raw count matrix (matching `stream_v2.load_chunk_sparse` and `benchmarks/rhj_fewshot/splithalf/CARD.md:32`: *"full-panel row sum over all file genes $\to \times 2{,}000 \to \log1p$"*; if $L_i = 0$, leave row at `0.0`).
     2. Apply `x_norm = np.log1p(x_raw / L_i * 2000.0)` in `float32`.
     3. Scatter **strictly the `6,983` RHJ axis Ensembl IDs** (`benchmarks/rhj_fewshot/gene_axis.csv`, SHA-256 `ef64f2ba7ff7816ed947796066cd8a68014884c6176619d750ef79462c78d4c0`) into their exact columns of the `44,728`-gene encoder tensor, setting **all `37,745` non-axis columns and any unmeasured axis genes strictly to `0.0`**.
     - *Rationale (directly from `benchmarks/rhj_fewshot/CARD.md:45-48`):* Occam's `rhj_fewshot` perturbation decoder and FiLM layers were trained exclusively on inputs with non-zero support on the `6,983`-gene axis (`axis_scope: input_and_loss`, `resolved_config.yaml:40`); feeding `37,745` non-axis genes causes an out-of-training-axis distribution shift warned against by the authors as *"silent and severe"*.
   - **Pre-Registered Secondary Input Sensitivity (`full_44728_api`):**
     Also evaluate Occam using its `occam.api.GeneSpace.prepare` default (library size over all Ensembl genes matching the `44,728`-gene encoder vocabulary, scattering all matched columns rather than zeroing the `37,745` non-axis columns), reported strictly as a descriptive sensitivity alongside primary `rhj_6983_axis`.

### §3.4 Locked Output-to-Delta Conversion (`baseline="reconstructed"` Primary vs. `baseline="real"` Sensitivity)
- **Primary Output-to-Delta Convention (`baseline="reconstructed"`, Occam `api.py:486, 510-521` Default & Exact `STATE` `PREREG_I §2.1.4` Analogue):**
  For every perturbation $p \in \mathcal{T}_{\mathrm{eval}}^{\mathrm{R16}}(t)$ and control pool $X_{\mathrm{ctrl}} = \{x_1, \dots, x_{n_{\mathrm{ctrl}}}\}$:
  $$\hat{\Delta}_{\mathrm{pred}}^{\mathrm{Occam}}(p) = \frac{1}{n_{\mathrm{ctrl}}}\sum_{i=1}^{n_{\mathrm{ctrl}}} \hat{Y}_{\mathrm{Occam}}\big(x_i, [\mathrm{pid}(p)]\big) - \frac{1}{n_{\mathrm{ctrl}}}\sum_{i=1}^{n_{\mathrm{ctrl}}} \hat{Y}_{\mathrm{Occam}}\big(x_i, [-1]\big),$$
  where `[-1]` invokes the exact zero-perturbation forward pass `D(z_ctrl + f(z_ctrl, 0))` (`src/occam/vcell/model.py:763-804`, `src/occam/api.py:519`), NOT raw `D(z_ctrl)`.
  - *Rationale:* Occam's own `docs/api.md:96-99` (`F-recon1`) and `src/occam/vcell/model.py:824-825` document that the decoder's control-reconstruction error $\|D(E(c)) - c\| = 15.34$ is **`~5×` larger than the perturbation effect (`2.88`)**, and `docs/h1_seen.md:94-99` documents `~9×` MAE inflation on H1 when subtracting real controls (`injctrl`). Using `baseline="reconstructed"` cancels this constant decoder offset and mirrors Round 16's locked `STATE_SE` delta formula ($\overline{\hat{Y}}(p) - \overline{\hat{Y}}(\text{"non-targeting"})$, `PREREG_I §2.1.4`).
- **Pre-Registered Secondary Output-Baseline Sensitivity (`baseline="real"` / `injctrl`):**
  Additionally compute $\hat{\Delta}_{\mathrm{pred}}^{\mathrm{Occam,real}}(p) = \frac{1}{n_{\mathrm{ctrl}}}\sum_{i=1}^{n_{\mathrm{ctrl}}} \hat{Y}_{\mathrm{Occam}}(x_i, [\mathrm{pid}(p)]) - \frac{1}{n_{\mathrm{ctrl}}}\sum_{i=1}^{n_{\mathrm{ctrl}}} \bar{x}_{i,\mathrm{ctrl}}$ (matching Occam's `official_tables.md` `nocal` evaluation convention), reported strictly as a descriptive sensitivity so readers can inspect both conventions side by side without mixing them (`docs/api.md:98-99`).

### §3.5 Critical Perturbation-Key Resolution Lock (Two-Tier `pert2id` HGNC Rule)
Primary inspection of Occam's `18,216`-gene perturbation graph (`master_gene_set_sorted.csv`) uncovered a critical alias trap:
- **`master_gene_set_sorted.csv` contains `317` duplicate HGNC pairs** where BOTH an older Replogle/Nadig symbol (e.g., `AARS`, `DARS`, `TAZ`, `H2AFZ`) AND the updated approved HGNC symbol (`AARS1`, `DARS1`, `TAFAZZIN`, `H2AZ1`) exist as **two separate graph nodes**!
- **`rhj_fewshot`'s `2,381` eligible trained perturbations (`benchmarks/rhj_fewshot/perturbations.csv`) have `0` HGNC collisions** — Occam's supervised training updated ONLY the legacy node! Naively looking up canonical HGNC symbols directly in `pert2id` would silently route `18` Group 2 targets to **untrained duplicate nodes**.
- **Locked Two-Tier Resolution Rule (`resolve_occam_pert_key_r16`):**
  1. **Tier 1 (Trained `rhj_fewshot` Node Priority):** Given canonical HGNC target $p \in \mathcal{T}_{\mathrm{eval}}^{\mathrm{R16}}(t)$, search the `2,381` eligible keys of `benchmarks/rhj_fewshot/perturbations.csv` for the unique key $k$ satisfying `HGNCBridge.canonicalize(k) == p` (`0` collisions verified). This resolves **`402/403` `a549_ko`** (bridging 14 pairs: `AARS1->AARS`, `CARS1->CARS`, `DARS1->DARS`, `EPRS1->EPRS`, `GARS1->GARS`, `HARS1->HARS`, `LARS1->LARS`, `MARS1->MARS`, `NOPCHAP1->C12orf45`, `PHB1->PHB`, `POLR1G->CD3EAP`, `POLR1H->ZNRD1`, `RARS1->RARS`, `TARS1->TARS`), **`75/75` `dld1`** (bridging 3 pairs: `DARS1->DARS`, `NARS1->NARS`, `TAFAZZIN->TAZ`), **`245/245` `cd4t_rest`** (bridging 1 pair: `H2AZ1->H2AFZ`), and **`88/88` `feng_ipsc`** (`0` bridged).
  2. **Tier 2 (Direct `pert2id` Fallback for Untrained Graph Nodes):** If $p$ has no canonical match in `rhj_fewshot`'s `2,381` eligible keys, resolve $p$ directly in `ckpt["pert2id"]` (resolving the sole remaining target **`RBM14`** in `a549_ko` to `pert2id["RBM14"]`, yielding **`403/403`** full parity with `STATE_SE`).
  3. **Descriptive Sensitivity on `a549_ko`:** Report both the primary `403`-target set (exact `STATE_SE` parity) and the `402`-target `rhj_2381`-trained-only subset (excluding `RBM14`).
- **Output Gene Projection:** Map Occam's `44,728`-gene output to the `6,983` RHJ Ensembl IDs (`gene_axis.csv`), canonicalize via frozen R10 `HGNCBridge` (`0` collisions), and gather the exact ordered `G_eval^R16(t)` columns (`1,837` `a549_ko`, `1,832` `dld1`, `1,830` `cd4t_rest`, `1,575` `feng_ipsc`). Additionally export a descriptive `Occam`-native full `RHJ_6983 ∩ G_eval^R9` sensitivity (`6,876`, `6,883`, `6,489`, `5,863` genes).

---

## §4. Metrics, Decision Cascades (`H-I1_Occam`, `H-I2_Occam`), and Paired `STATE_SE` vs. `Occam` Contrasts (Item 4)

### §4.1 Metrics & Decision Rules Copied Verbatim from `PREREG_I.md` + `PREREG_I_ADDENDUM_A.md`
All primary endpoints, statistical parameters, and decision cascades are inherited **verbatim** without modification:
- **Primary Endpoint:** `r_orth_cal` (uncentered orthogonalization to per-perturbation `Transfer_Mean` $\bar{\Delta}^{\mathrm{src}}_p$ via frozen `round10/src/metrics_r10.py::_orth_pair`, pooled moment correlation across perturbations on unselected **`all_nontarget`**, Spearman–Brown disattenuation using **`rho_ab_orth`** clipped to `[-1, 1]`).
- **Secondary Endpoints:** `r_dorth_cal` (Systema-inspired projection onto $\operatorname{span}\{\mathbf{1}, \bar{\Delta}^{\mathrm{src}}_p, \bar{b}_S\}^\perp$, `PREREG_I §4.1`), plus practical utility metrics (raw Pearson $\Delta$, centered DE-50 Pearson, cosine similarity, and `M2` `resp_top50_split - all_nontarget` shift).
- **Bootstrap & Multiple-Testing Parameters:** `n_boot = 4,000` joint perturbation bootstraps, `seed = 20261601` (identical resample indices across matched, mismatched, `FCT`, and `STATE_SE` for exact pairing); Bonferroni $\alpha_{\mathrm{claim}} = 0.05 / 2 = 0.025$ across the two Occam claims (`H-I1_Occam` and `H-I2_Occam`):
  - Two-sided **`97.5%` percentile CI (`[1.25%, 98.75%]`)** for superiority/inferiority;
  - Two-sided **`95.0%` percentile CI (`[2.5%, 97.5%]`)** vs. open interval **$(-\delta, +\delta)$ with $\delta = 0.04$** for TOST equivalence (`valid_frac >= 0.95` required).
- **Per-Pair 5-Way Priority Cascade (`PREREG_I §3.1` + `Addendum A §10`):**
  1. `INCONCLUSIVE` if `valid_frac < 0.95` or `97.5%` CI non-finite (`within_delta = False`);
  2. `EXACT_INVARIANCE` if both predicted delta matrices are finite AND $d_{\max} = \max_{p,g} |\hat{\Delta}^{\mathrm{matched}}_{p,g} - \hat{\Delta}^{\mathrm{mis}(s')}_{p,g}| \le 10^{-12}$ (counted as equivalent in stratum/line rules: `n_equiv_effective = counts["EQUIVALENT"] + counts["EXACT_INVARIANCE"]`);
  3. `INCONCLUSIVE` (conservative zero-tie precedence) if `ci975_lo == 0.0` OR `ci975_hi == 0.0`;
  4. `MATCHED_BETTER` if `ci975_lo > 0.0` (with `within_delta` flag);
  5. `MISMATCHED_BETTER` if `ci975_hi < 0.0` (with `within_delta` flag);
  6. `EQUIVALENT` if `within_delta == True`;
  7. `INCONCLUSIVE` otherwise (including exact `±0.04` boundary ties).
- **Line-Level `H-I1_Occam` 8-State Priority Cascade (`PREREG_I §3.1`):**
  1. `MISSING_STRATUM` (`n_seen == 0` or `n_unseen == 0`);
  2. `INCONCLUSIVE` (`NOT matched_arm_valid`);
  3. `LINE_MATCHED_SPECIFIC` (`matched_positive_ci975 == True`, ALL pairs in BOTH strata `MATCHED_BETTER`, `NOT all_within_delta`);
  4. `LINE_MATCHED_SPECIFIC_SUB_DELTA` (`matched_positive_ci975 == True`, ALL pairs in BOTH strata `MATCHED_BETTER`, `all_within_delta == True` — **MUST** follow `Addendum A §9` mandatory wording: *"能检测到，但小于事先规定的最小有意义效应 δ=0.04"*, reported jointly with §8 within-line resample, NEVER reported standalone as *"模型用到了对照"*);
  5. `LINE_EQUIVALENT_NONZERO_SIGNAL` (`matched_positive_ci975 == True` AND ALL pairs in BOTH strata `EQUIVALENT`/`EXACT_INVARIANCE`);
  6. `LINE_EQUIVALENT_NULL_SIGNAL` (`matched_positive_ci975 == False` AND ALL pairs in BOTH strata `EQUIVALENT`/`EXACT_INVARIANCE`);
  7. `LINE_MISMATCHED_BETTER` (ALL pairs in BOTH strata `MISMATCHED_BETTER`);
  8. `MIXED_OR_INCONCLUSIVE` (all remaining cases; any `INCONCLUSIVE` interpreted as insufficient power `检验力不足`, never as equivalence ruled out).
- **Within-Line Disjoint Resample Arm (`Addendum A §8`):** Compute $\Delta^{\mathrm{pair}}(\mathrm{resample}) = r(\mathrm{matched}) - r(\mathrm{resample})$ and per-source excess loss $\Delta^{\mathrm{excess}}(s') = \Delta^{\mathrm{pair}}(s') - \Delta^{\mathrm{pair}}(\mathrm{resample})$ using the locked `seed = 20261602` disjoint controls (`0` overlap) and identical `Addendum A §8.3` narrative rules.
- **Claim `H-I2_Occam` (`Occam` vs. Closed-Form `FCT`, `PREREG_I §3.2`):** Compute $\Delta^{\mathrm{vsFCT}}_{\mathrm{Occam}, t} = r^{\mathrm{orth\_cal}}(\mathrm{Occam}_{\mathrm{matched}}(t)) - r^{\mathrm{orth\_cal}}(\mathrm{FCT}(t))$ with `97.5%` CI (`MODEL_SUPERIOR_TO_FCT` / `FCT_SUPERIOR_TO_MODEL` / `INDISTINGUISHABLE`), subject to the identical **Anti-Masking Rule (`M6.2`)** (line-by-line reporting; winning/tying on weak-`FCT` `cd4t_rest` NEVER masks inferiority on strong-`FCT` `a549_ko` or target-control failures).
- **Group-Level Reporting:** Report Group 2 CRISPRko (`a549_ko`, $K_{2,\mathrm{ko}}=1$), Group 2 CRISPRi (`{dld1, cd4t_rest}`, $K_{2,\mathrm{i}}=2$), optional Group 3 (`gse295214_h1`, $K_3 \in \{0,1\}$), and `SENSITIVITY_ONLY` (`feng_ipsc`) strictly separately; never pool CRISPRko with CRISPRi or merge Group 3 with Group 2.

### §4.2 Pre-Registered Paired `STATE_SE` vs. `Occam` Contrasts on Shared Group 2 Targets (and Group 3 if Qualified)
Because Occam and `STATE_SE` share **100.0% identical targets `N_eval^R16(t)`, evaluation genes `G_eval^R16(t)`, observed/split-half arrays, `Transfer_Mean`, and locked `512`-cell control pools**, Addendum B pre-registers two paired head-to-head contrasts using the exact same `n_boot = 4,000` joint perturbation bootstrap draws (`seed = 20261601`) on `r_orth_cal` (`all_nontarget`; mirrored descriptively on `r_dorth_cal`):

1. **Contrast `H-I3a` — Matched-Control Accuracy Contrast (`Occam` vs. `STATE_SE`):**
   $$\Delta^{\mathrm{acc}}_{t} = r^{\mathrm{orth\_cal}}\!\big(\mathrm{Occam}_{\mathrm{matched}}(t)\big) - r^{\mathrm{orth\_cal}}\!\big(\mathrm{STATE\_SE}_{\mathrm{matched}}(t)\big),$$
   using `STATE_SE`'s locked primary `fixed_index_zero` matched arm and `Occam`'s locked primary `rhj_6983_axis` + `baseline="reconstructed"` matched arm. Classified via two-sided **`97.5%` CI (`[1.25%, 98.75%]`)** with `valid_frac >= 0.95` and conservative zero-tie rule (`ci975_lo == 0.0 or ci975_hi == 0.0 -> INDISTINGUISHABLE`):
   - `OCCAM_SUPERIOR_TO_STATE_SE` (`ci975_lo > 0.0`);
   - `STATE_SE_SUPERIOR_TO_OCCAM` (`ci975_hi < 0.0`);
   - `INDISTINGUISHABLE` (otherwise).
2. **Contrast `H-I3b` — Matched-vs-Mismatched Control-Specificity Contrast (`Occam` vs. `STATE_SE`):**
   For every mismatched source $s' \in \mathcal{S}_{\mathrm{mis}}^{\mathrm{seen}}(t) \cup \mathcal{S}_{\mathrm{mis}}^{\mathrm{unseen}}(t)$ AND for each stratum mean ($\overline{\Delta}^{\mathrm{pair}}_{\mathrm{seen}}$, $\overline{\Delta}^{\mathrm{pair}}_{\mathrm{unseen}}$ separately):
   $$\Delta^{\mathrm{spec}}_{t}(s') = \Delta^{\mathrm{pair}}_{\mathrm{Occam}, t}(s') - \Delta^{\mathrm{pair}}_{\mathrm{STATE\_SE}, t}(s') = \Big[r(\mathrm{Occam}_{\mathrm{matched}}) - r(\mathrm{Occam}_{\mathrm{mis}(s')})\Big] - \Big[r(\mathrm{STATE\_SE}_{\mathrm{matched}}) - r(\mathrm{STATE\_SE}_{\mathrm{mis}(s')})\Big],$$
   plus the analogous resample-adjusted excess contrast $\Delta^{\mathrm{spec,excess}}_t(s') = \Delta^{\mathrm{excess}}_{\mathrm{Occam}, t}(s') - \Delta^{\mathrm{excess}}_{\mathrm{STATE\_SE}, t}(s')$. Classified via the same `97.5%` CI rule (`OCCAM_MORE_CONTROL_SPECIFIC` / `STATE_SE_MORE_CONTROL_SPECIFIC` / `INDISTINGUISHABLE`).
   *Note on Multiple-Testing Scope:* `H-I3a` and `H-I3b` are pre-registered **comparative secondary contrasts** reported per line alongside `H-I1_Occam` and `H-I2_Occam` without re-spending or altering Round 16's locked primary `STATE_SE` Bonferroni family (`PREREG_I §3`).

### §4.3 Pre-Locked Advance Interpretation Matrix for `STATE_SE` vs. `Occam` Outcomes
Every major joint outcome on a target line $t$ is assigned a **fixed advance interpretation** before any Occam inference is run:

| Joint Outcome Pattern on Target Line $t$ | Pre-Locked Scientific Interpretation (No Post-Hoc Reframing Permitted) |
|---|---|
| **Case 1:** Both `STATE_SE` and `Occam` fail target-control specificity (`LINE_EQUIVALENT_NONZERO_SIGNAL`, `LINE_EQUIVALENT_NULL_SIGNAL`, or `SUB_DELTA` not exceeding within-line resample) while `FCT` succeeds | **Cross-architecture replication of target-control under-utilization:** Failure to exploit target-line control cells in unseen contexts generalizes across both Transformer-embedding (`STATE_SE`) and multiplicative-FiLM MLP (`Occam`) virtual-cell architectures, ruling out a `STATE`-only artifact (such as `batch_encoder` dominance). |
| **Case 2:** `Occam` achieves `LINE_MATCHED_SPECIFIC` (and exceeds within-line resample loss) while `STATE_SE` does not (`OCCAM_MORE_CONTROL_SPECIFIC` in `H-I3b`) | **Architecture-specific control sensitivity:** Multiplicative FiLM modulation on raw control-cell latents (`Occam`) extracts target-line-specific control signal where `STATE_SE` does not — falsifying any universal claim that all neural virtual-cell models ignore target controls. (Must still be read alongside `H-I2_Occam` vs. `FCT`.) |
| **Case 3:** `STATE_SE` achieves `LINE_MATCHED_SPECIFIC` while `Occam` does not (`STATE_SE_MORE_CONTROL_SPECIFIC` in `H-I3b`) | **Scale/representation advantage for `STATE_SE`:** `SE-600M`'s 167M-cell pretrained embeddings capture target-line control specificity that Occam's `DiverseRNA-1.4M` + `6,983`-axis encoder fails to preserve under cross-study shift. |
| **Case 4:** `Occam` wins matched accuracy (`OCCAM_SUPERIOR_TO_STATE_SE` in `H-I3a`) YET fails `H-I1_Occam` (`LINE_EQUIVALENT_NONZERO_SIGNAL` or `<=` within-line resample) | **Control-agnostic perturbation-prior gain:** Occam's higher `r_orth_cal` arises from stronger target-invariant perturbation features (GO/STRING/GenePT/ESM2) or cross-line residual alignment rather than target-cell-line control conditioning — demonstrating why matched-control accuracy alone cannot establish context specificity. |
| **Case 5:** Both `STATE_SE` and `Occam` achieve `LINE_MATCHED_SPECIFIC` AND exceed within-line resample loss | **Cross-architecture target-control utilization:** Both neural architectures genuinely utilize target-line control profiles beyond `Transfer_Mean`; `H-I2` (`model - FCT`) and `H-I3a`/`H-I3b` then quantify whether either surpasses or falls short of closed-form `FCT`. |
| **Case 6:** Either or both return `MIXED_OR_INCONCLUSIVE` (especially `dld1` at `N=75`) or `LINE_MISMATCHED_BETTER` | **Bounded / power-limited or harmful-conditioning outcome:** `INCONCLUSIVE` is strictly reported as `检验力不足` (never "equivalence ruled out"); `LINE_MISMATCHED_BETTER` indicates out-of-distribution degradation when conditioning on target controls relative to familiar training-line (`seen`) or cleaner (`unseen`) controls, reported stratum-by-stratum. |

---

## §5. Pre-Locked Two-Stage Adapter Reproduction Gate (Item 5)

Analogous to `PREREG_I §1.3` and `Addendum A §3`, our `OccamAdapter` MUST pass **both stages** below **BEFORE ANY Group 2, Group 3, or `feng_ipsc` prediction is scored**. All pass criteria are locked prior to download.

### §5.1 Stage 5A — Structural Integrity, Pure-Synthetic Control Sensitivity, and MPS-vs-CPU Parity (0 Real Cells)
1. **Checksum & Checkpoint Contract Assertions:**
   - Local `best_model.pth` SHA-256 equals `11d5ba5515e79710b3a952532c60aa29515046117292528c59d452a0937685ff`;
   - `ckpt["is_virtual_cell_checkpoint"] is True`;
   - `ckpt["target_sum"] == 2000.0`;
   - `ckpt["decoder_config"]["output_activation"] == "softplus"`;
   - `ckpt["interaction_config"]` matches `interaction_type="film"`, `film_mode="multiplicative_only"`, `film_mean_free=True`, `film_hidden_dim=256`;
   - `len(ckpt["gene_names"]) == 44728`, `len(ckpt["pert2id"]) == 18216`;
   - `(ckpt.get("benchmark") or {}).get("gene_axis_sha256") == "ef64f2ba7ff7816ed947796066cd8a68014884c6176619d750ef79462c78d4c0"`;
   - Model has **zero** discrete cell-line/batch/source embedding parameters (`source_emb is None`, `batch_encoder` absent); encoder + decoder + FiLM + GNN load cleanly with `strict=True` (no missing/unexpected keys).
2. **Pure-Synthetic Smoke & Control-Sensitivity Assertions (`n = 32` synthetic cells, `seed = 20261601`):**
   - Outputs are finite and strictly positive (`softplus` contract; `dead_gene_frac < 0.01`, ruling out `rectified_sigmoid` mismatch);
   - Swapping synthetic control matrix `X_ctrl_A != X_ctrl_B` (with perturbation index fixed) produces strictly non-zero output shift (`max_abs_delta_diff > 1e-4`);
   - Single-perturbation identity invariant holds: `forward_composed(..., "decode_seq_mean")` equals `forward_composed(..., "sum")` within `<= 1e-6`.
3. **MPS-vs-CPU `float32` Device Parity Check:**
   - Compare identical synthetic batch forward passes on `CPU float32` vs. `MPS float32`:
     - If `max_abs(z_ctrl_mps - z_ctrl_cpu) <= 1e-5` AND `max_abs(delta_mps - delta_cpu) <= 1e-4`, lock **`MPS float32`** as execution device (also verify cross-process `PYTHONHASHSEED` determinism `1` vs. `99991`: `max_abs == 0.0` CPU, `<= 1e-6` MPS);
     - Otherwise automatically fall back and lock **`CPU float32`** for both Main Line and `A-T3`.

### §5.2 Stage 5B — Quantitative Author-Benchmark Reproduction Gate on Local `nadig_hepg2` (Zero Extra Dataset Downloads)
- **Dataset Parity Verified Pre-Lock:** Local `GSE264667_hepg2_raw_singlecell_01.h5ad` (already on mounted sync drive; SHA-256 `e1ad7c3c5a201c861a207a858aa7e59f5e6ac1955674c415f7de0d1dadadb52e`) has exact shape **`145,473 cells × 9,624 genes`**, `4,976` non-targeting controls, and **`6,983 / 6,983`** RHJ axis Ensembl IDs (`var["gene_id"]`) — matching Occam's `nadig_2024_hepg2.h5ad` (`145,473 × 9,624`, `4,976` controls, `6,983` axis genes; repackaged with `var["ensembl_gene_id"]` and `obs["dataset_batch_num"]`, SHA-256 `f31a85d8...`, `benchmarks/rhj_fewshot/splithalf/CARD.md:30`).
- **Exact Reproduction Protocol:**
  1. Using local `GSE264667_hepg2_raw_singlecell_01.h5ad`, construct the `611`-perturbation held-out HepG2 test cohort defined in `benchmarks/rhj_fewshot/splithalf/cohort.json` (`42,402` cells = `37,426` perturbed + `4,976` controls), normalizing each cell over all `9,624` file genes to `target_sum = 2000.0 -> log1p` and projecting onto the `6,983`-gene RHJ axis.
  2. Perform `5` batch-stratified (by `gem_group`) disjoint A/B splits (`m_p = floor(n_p / 2)`, `2,488` controls per half, max per-batch imbalance `<= 1`).
  3. Run `OccamAdapter` (`seed_42`, `rhj_6983_axis`, **`injctrl` / `baseline="real"` `nocal`** convention matching Occam's `official_tables.md` & `selected_scores_by_seed_split.csv`: predict `m_p` cells per perturbation from sampled half-B controls, compute predicted delta against full half-B control mean, score mean `pearson_delta` and `mae` vs. half-A reference across the `611` perturbations and `5` splits).
- **Pre-Locked Reference Values (`bendidiihab/Occam/official_comparisons/selected_scores_by_seed_split.csv`, `rhj_fewshot | all | ABL_N_ALL | nocal`, `seed_42`, `n_perts = 611`, `5` splits):**
  - Author `ABL_N_ALL` `seed_42` `nocal` **`pearson_delta_mean = 0.394178`** (5-split values: `0.392850, 0.397797, 0.395306, 0.392278, 0.392659`);
  - Author `ABL_N_ALL` `seed_42` `nocal` **`mae_mean = 0.030157`** (`mae_delta_mean = 0.030284`);
  - Author `experimental` replicate `nocal` **`pearson_delta_mean = 0.372174`** (5-split values: `0.371669, 0.370873, 0.373463, 0.372176, 0.372690`).
- **Locked Pass Criteria (ALL 3 Required; Accounts for Local Split-Seed & Cross-Version PyTorch Tolerance `±0.015` / `±0.005`):**
  1. `abs(adapter_mean_pearson_delta - 0.394178) <= 0.015` (i.e. `[0.3792, 0.4092]`);
  2. `abs(adapter_mean_mae - 0.030157) <= 0.005` (i.e. `[0.0252, 0.0352]`);
  3. `abs(local_experimental_pearson_delta - 0.372174) <= 0.015` (confirming local panel/split construction fidelity).
- **Failure Protocol:** If ANY Stage 5A or Stage 5B check fails, **halt immediately before running any Group 2/3/sensitivity inference** and escalate to Master Controller to diagnose the adapter. No adjudication rule or threshold may be altered.

---

## §6. Compute Budget and Clean-Room Independent Audit Route (Item 6)

### §6.1 Compute Plan on Apple M3 Pro (`36 GB` Unified Memory)
- **Hardware & Environment:** `Apple M3 Pro` (`36 GB` unified RAM), isolated `round16/.venv` (`Python 3.11`, PyTorch `2.x`). Note that Occam's `pyproject.toml` pins `Python 3.8` + `torch==2.4.1`, and `scripts/verify_checkpoint.py:64-68, 106-107` explicitly documents cross-PyTorch-version float drift tolerance `CROSS_VERSION_ATOL = 1e-4, CROSS_VERSION_RTOL = 1e-3`; our Gate 5A + 5B checks empirically validate numerical fidelity in `round16/.venv` without installing a separate legacy Python 3.8 stack unless Gate 5B fails.
- **Memory Footprint:** Pure PyTorch MLP encoder (`94.3M` params) + 4-block SE residual decoder (`128.5M` params) + baked `ConcatKnowledgeModel` perturbation table & FiLM (`~250M` params/buffers) — zero FlashAttention or custom CUDA kernels. Loaded in `float32`, model + batch activations occupy **`~3.5–5.5 GB` peak RSS** (`<< 36 GB`).
- **Caching & Runtime Optimization:**
  - `z_ctrl = encoder(X_ctrl)` (`512 × 512`) is computed **strictly ONCE per control pool** (`6` primary pools + `4` within-line resample pools = `10` pools = `5,120` cells total $\implies$ **`< 10 seconds`**) and cached alongside `z_graph = pert_model()` (`18,216 × 128`, computed once).
  - For each `(target line t, control arm)` only the lightweight FiLM modulation (`512 × 512`) and decoder pass (`512 × 44,728`) run per perturbation (`403 + 75 + 245 + 88 = 811` target perts $\times$ `7` arms [`matched` + `resample` + `5` mismatched] + `baseline` `[-1]` $\approx 5,700$ perturbation evaluations $\times 512$ cells):
    - **Primary (`MPS float32`, if Gate 5A parity passes):** estimated **`~3–10 minutes`** total (plus `~2–4 minutes` for Stage 5B HepG2 reproduction gate and `~3–8 minutes` for `full_44728_api` sensitivity);
    - **Fallback (`CPU float32`):** estimated **`~15–35 minutes`** total.
- **Strict GPU/MPS Contention Rule:** **ZERO GPU/MPS execution** is permitted until Master Controller explicitly confirms that `A-T3` has completed its parallel Round 16 run, locks Addendum B, records `R16_ADDB_HUMAN_APPROVAL.txt`, freezes Addendum B code, and issues `R16_ADDB_GO.txt`.

### §6.2 Independent Clean-Room Audit Route (`A-T3` + Master Controller)
After Addendum B lock (`PREREG_I_ADDENDUM_B_COMMIT.txt`), code freeze (`CODE_FREEZE_R16_ADDB.txt`), and release (`R16_ADDB_GO.txt`), execution follows the exact staged three-party verification protocol of `PREREG_I_ADDENDUM_A §7` & `§11`:
1. **Clean-Room Dual Implementation:** Main Line (`round16_addB/src/`) and `A-T3` (independent clean-room adapter in its audit workspace) independently implement `OccamAdapter` strictly from this locked addendum and Occam's pinned commit `9ddb890e65ca1ec2f81882c78311002a75b52f67`, passing identical synthetic parity (`<= 1e-10`) and cross-process `PYTHONHASHSEED` (`1` vs. `99991`) determinism tests before `CODE_FREEZE_R16_ADDB.txt`.
2. **Staged Three-Party Export & Acceptance Tolerances:**

| Audit Layer | Parties Compared | Locked Acceptance Criterion |
|---|---|---|
| **1. Control-cell barcodes** | Main Line vs. `A-T3` vs. locked `control_barcodes_r16.json` | **Exact `100%` ordered string equality** (primary `20261601` & resample `20261602`) |
| **2. Occam `z_ctrl` encoder latents** | Main Line vs. `A-T3` (exported `z_ctrl_<pool>.npz`, `512 × 512`, checked **first**) | $\max_{i, d} |z_{\mathrm{ctrl}}^{\mathrm{Main}}[i, d] - z_{\mathrm{ctrl}}^{\mathrm{AT3}}[i, d]| \le 10^{-5}$ |
| **3. Predicted delta arrays** | Main Line vs. `A-T3` (exported `preds_Occam_<arm>.npz`, primary & sensitivities) | $\max_{p, g} |\hat{\Delta}^{\mathrm{Main}}_{p, g} - \hat{\Delta}^{\mathrm{AT3}}_{p, g}| \le 10^{-4}$ |
| **4. Same-array metric recomputation** | Master Controller recomputing from exported `.npz` arrays (`seed = 20261601`, `n_boot = 4,000`) | $|\mathrm{metric}^{\mathrm{Controller}} - \mathrm{metric}^{\mathrm{Party}}| \le 10^{-10}$ |
| **5. Cross-party metrics & verdicts** | Main Line vs. `A-T3` | $|\mathrm{metric}^{\mathrm{Main}} - \mathrm{metric}^{\mathrm{AT3}}| \le 10^{-5}$ **AND `100%` verdict parity** (`H-I1_Occam`, `H-I2_Occam`, `H-I3a`, `H-I3b`) |

---

## §7. Limitations and Exact Boundary vs. Occam's Own `"Transfer"` Baseline (Item 7)

### §7.1 Pre-Specified Limitations (Disclosed in All Reports)
1. **Unattributed `DiverseRNA-1.4M` Pretraining Blocks:** 7 of 9 encoder pretraining datasets are unidentified (`MODEL_CARD.md:184`). While `a549_ko` (`2026-07/08`) and `dld1` (`2026-07`) strictly post-date `arXiv:2605.31562`, `cd4t_rest` (`2026-02-06`) and `gse295214_h1` (`2026-05-05`) cannot be guaranteed absent from unsupervised encoder pretraining (supervised perturbation absence is 100% certain for all Group 2/3 lines).
2. **`cd4t_rest` Input Gene-Axis Coverage Deficit:** `413 / 6,983` (`5.9%`) of Occam's RHJ input axis genes are absent from `cd4t_rest`'s `18,130`-gene probe panel and must be zero-filled at input (analogous to `h1_seen`'s `448` zero-filled genes in `docs/h1_seen.md:37`). Because this asymmetry affects only `cd4t_rest`'s matched arm (and `cd4t_rest` when used as a mismatched source), any matched-arm penalty on `cd4t_rest` MUST be interpreted jointly with its `413`-gene zero-fill rate and cross-checked against `full_44728_api`.
3. **`a549_ko` Modality Shift & `RBM14`:** `a549_ko` is a CRISPRko (`SpCas9`) screen whereas `rhj_fewshot` is CRISPRi (`dCas9-KRAB`), and `1/403` target (`RBM14`) was not in `rhj_2381` (addressed via the `402`-target sensitivity).
4. **Single Training Seed (`seed_42`):** Addendum B evaluates `seed_42` (consistent with Occam's default release checkpoint and disk budget `2.55 GB` vs. `12.75 GB` for all 5 seeds). Across Occam's published 5 seeds (`42–46`), `seed_42` sits closely aligned with the 5-seed mean (`rhj_fewshot` `0.3942`; `h1_seen` `0.1162` vs. 5-seed `0.1071`; `k562_seen` `0.2302` vs. 5-seed `0.2115`), yet single-seed variance is not estimated here.
5. **Decoder Control-Reconstruction Bias (`F-recon1`) & Library-Scale Convention:** Occam predicts in `log1p(norm_2000)` space whereas Round 16 observed deltas (`delta_obs`) use `log1p(norm_1e4)`. While Pearson correlation (`r_orth_cal`) is invariant to global affine rescaling, nonlinear `log1p` curvature differs modestly across `2,000` vs. `10,000` library scaling, and `baseline="reconstructed"` vs. `baseline="real"` removes vs. retains the `~5×` `F-recon1` offset; reporting both `baseline="reconstructed"` (primary) and `baseline="real"` (sensitivity) keeps this trade-off completely transparent.
6. **Retrospective Epistemic Scope & Statistical Power:** Group 2 cohorts were previously observed in Rounds 10–12 (`预注册的回顾性模型审计`); `dld1` (`N = 75`) has low/borderline TOST power (`PREREG_I §3.3`), and `cd4t_rest` has a wide relative $\delta = 0.04$ (`71.9%` of `FCT`).

### §7.2 Exact Overlap and Scientific Distinction vs. Occam's Own `"Transfer"` Baseline (`arXiv:2611.21959`)
Occam's paper (`arXiv:2611.21959`) already includes a `"Transfer"` baseline and candidly reports that `"Transfer"` beats Occam on raw `Pearson Δ` in held-out cell lines:
- **Verbatim definition (`Appendix E.2`):**
  > *"Transfer. Add the same perturbation’s measured effect from other training cell lines to each sampled B control from the evaluation cell line. The source effects are weighted by their cell counts. This baseline uses perturbation identity but does not learn how the effect changes between cell lines."*
- **Verbatim limitation (`Appendix A`):**
  > *"Transferring measured effects from other training cell lines gives higher Pearson Δ than Occam in held-out cell lines, and no model leads on every metric."*

**Where Round 16 + Addendum B overlaps and where it is strictly distinct:**
1. **Acknowledged Overlap:** Both Occam's paper and Round 16 recognize that copying a perturbation's cross-line average delta onto target controls is a formidable baseline on unorthogonalized correlation (and Occam itself loses to `Transfer` on `h1_seen` [`0.116` vs. `0.233` for `seed_42`] and `k562_seen` [`0.230` vs. `0.271`]).
2. **Four Concrete Differences (Why `"Transfer"` Does NOT Answer Addendum B's Question):**
   - **(a) Exact Removal of `"Transfer"` via Orthogonalization (`r_orth_cal`):** Under `r_orth_cal`, the additive cross-line transfer vector (`Transfer_Mean`) is projected out of every predicted and observed delta, so `"Transfer"` evaluates identically to **`0.0` by construction**. Addendum B directly isolates the **cell-line-specific residual** that `"Transfer"` misses.
   - **(b) Direct Input Control-Swap Intervention on Frozen Occam Weights:** Occam's paper **never swaps Occam's own input `control_expression`** across matched vs. mismatched (`seen_by_checkpoint` and `unseen_by_checkpoint`) cell lines, never runs TOST equivalence testing ($\delta = 0.04$), and never compares mismatched-control loss against a within-line disjoint control-resampling null (`Addendum A §8`).
   - **(c) Comparison Against Target-Control-Aware Multiplicative `FCT` (`H-I2_Occam`):** Unlike Occam's additive `"Transfer"` baseline (which ignores target baseline expression levels when computing deltas), `FCT` explicitly modulates cross-line fold changes by the target line's own control expression.
   - **(d) Paired `STATE_SE` vs. `Occam` Contrast on Unseen Post-Pretraining Group 2 Screens:** Addendum B tests `Occam` and `STATE_SE` head-to-head on identical `(targets, genes, control cells)` across independent post-2025 screens (`a549_ko`, `dld1`, `cd4t_rest`) that never appeared in `rhj_fewshot` or `ST-SE-Replogle`.
