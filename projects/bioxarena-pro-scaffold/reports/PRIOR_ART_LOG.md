# Prior-Art Verification Log (`R8-8`)

- **Compiled at**: 2026-11-21 16:08+08:00
- **Search window**: Past 24 months (`2024-11` through `2026-11`), plus foundational 2024 anchors where relevant.
- **Databases / indexes queried**: bioRxiv, arXiv, PubMed / PMC (`nih.gov`), Nature / Cell / Springer Nature family, ICML / NeurIPS / ICLR proceedings, and Google Scholar web index via `search_web`.
- **Governance Rule (Mandatory Gate)**: No novelty claim may enter `docs/PREPRINT_DRAFT.md`, `docs/preprint/*.pdf`, or any external material without an entry in this log. Words such as **"first"**, **"novel discovery"**, **"首次"**, or **"首个"** are strictly prohibited unless explicitly supported by a verified `新` verdict below.

---

## Executive Verdict Table

| Claim ID | Claim Being Tested | Verdict | Can Claim "First / 首次"? | Permitted Framing |
| :--- | :--- | :--- | :--- | :--- |
| **PA-1** | `S5` cross-cell-line bridging for single-gene perturbation response prediction | **已有 (概念) / 部分重叠 (具体形式)** | **NO** | Closed-form, zero-neural-parameter, source-centered dual-ridge kernel transfer baseline evaluated zero-shot across 3 human CRISPRi atlases |
| **PA-2** | Explaining linear/frozen-probe baseline advantage on $\Delta$-Pearson via PC1 shared-axis decomposition | **部分重叠** | **NO** | Pre-registered rank-1 PC1 off-axis projection test ($S = \operatorname{span}\{1, b, a\}$) extending Ahlmann-Eltze (2025), Systema (2025/2026), and Miller (2026) across 3 CRISPRi atlases |
| **PA-3** | Parameter-free orthogonal fusion of foundation-model embeddings and kernel transfer models (`OrthoFusion`) | **部分重叠** (and **UNCONFIRMED** on 4th dataset) | **NO** | Exploratory parameter-free rank-1 subspace projection combining scGPT probe and `S5`; confirmatory 4th-dataset validation remains `PENDING` |
| **PA-4** | Pre-registered, commit-hashed zero-look confirmation with hard-stop dataset rejection in perturbation benchmarking | **部分重叠** | **NO** | Operational implementation of commit-hashed pre-registration and negative-gate enforcement in single-cell perturbation benchmarking |
| **PA-5** | Dual-AI autonomous research workflow (Executor AI + Cross-Vendor Auditor AI + Git audit ledger + Pre-registration) | **已有 (核心范式已有) / 部分重叠 (生物学全流程公开错误台账实证)** | **NO** | Concrete empirical case study of a cross-vendor audited autonomous research loop in computational biology with a public 34-entry unvarnished error ledger (explicitly situated alongside `CrossAudit`, *Virtual Lab*, *A Society of AI Scientists*, PCF, and `arXiv:2610.22465`) |

> **Summary Gate Finding**: **NONE** of the five claims qualifies as unconditionally `新` (first-in-field). All five have direct prior art or substantial overlap within the 2024-11 to 2026-11 window—most critically including **`CrossAudit` (`arXiv:2608.28631`)** for git-native cross-vendor AI auditing and **Miller et al. (*Nature Biotechnology* 2026 / `bioRxiv 10.1101/2025.10.20.683304`)** for perturbation metric calibration and ranking inversion. All "first / 首次 / 首个" language is banned across the manuscript and external materials.

---

## Detailed Search Records

### PA-1: Cross-Cell-Line Transfer / Bridging for Perturbation Response Prediction (`S5`)

- **Databases queried**: bioRxiv, arXiv, Nature/Cell family, GitHub/Google Scholar via `search_web`.
- **Exact search queries**:
  1. `"cross-cell-line" OR "cross-context" perturbation response prediction transfer learning 2024 2025 2026 bioRxiv Nature`
  2. `"scLAMBDA" OR "DECANT" OR "StateXDiff" OR "HyperMap" OR "GRAFT" OR "TranScouter" perturbation prediction bioRxiv 2025 2026`
  3. `"STATE" "Adduri" Arc Institute 2025 perturbation prediction "scGPT" Cui Nature Methods 2024 "GEARS" Roohani Nature Biotechnology 2024`
- **Time window**: 2024-11 to 2026-11 (plus 2024 anchors `GEARS`, `scGPT`).
- **Matched Prior Art**:
  1. **Adduri, Roohani, Hsu et al. (2025), *STATE*** (Arc Institute; NeurIPS 2025 / preprint): Transformer architecture trained on >100M perturbed and 167M observational cells explicitly targeting cross-context generalization across unobserved cell lines and perturbation types.
  2. **HyperMap (bioRxiv, April 2026)**: *"HyperMap: An Efficient Framework for Transferring Perturbation Responses Across Diverse Biological Contexts"* — meta-learning framework mapping existing perturbation atlases to new biological contexts (distinct cell lines, iPSC donors, small molecules) using small sets of shared perturbation "seeds", plus `HyperMapDB`.
  3. **TranScouter (bioRxiv, August 2026)**: *"A structured study of cross-condition prediction of transcriptional responses to gene perturbations"* — lightweight encoder-decoder for cross-condition / cross-cell-line prediction under seen and unseen perturbation splits.
  4. **GRAFT (bioRxiv, November 2026)**: *"Cross-modal transfer enables zero-shot prediction of drug responses across cell types"* — cross-modal and cross-cell-type transfer anchoring chemical responses on genetic perturbation networks.
  5. **Virtual Cell Challenge benchmarks (2025–2026, bioRxiv)**: Explicitly evaluate cross-cell-line CRISPRi transfer (e.g., H1 embryonic stem cells and KOLF variants), including **source-centering** and differential-expression kernel baselines.
  6. **Wilson et al. (bioRxiv, October 2026)**: Diagnostic analysis of *Context Geometry Compression (CGC)* in cross-context virtual cell models, showing models capture shared perturbation axes while collapsing cell-context specificity.
  7. **Roohani et al. (2024), *GEARS*** (*Nature Biotechnology*) & **Cui et al. (2024), *scGPT*** (*Nature Methods*): Foundational graph-based and transformer baselines for single-cell perturbation prediction.
- **Verdict**: **已有 (概念) / 部分重叠 (具体实现)**.
- **Impact on Manuscript**: Cross-cell-line bridging via shared perturbations is an established active subfield (STATE, HyperMap, TranScouter, Virtual Cell Challenge). `S5` must be framed strictly as a **lightweight, closed-form, source-centered dual-ridge kernel baseline** with zero neural training parameters whose fixed constants transfer without retuning across three human CRISPRi atlases.

---

### PA-2: Explaining Linear / Frozen-Probe Baseline Advantage via PC1 Shared-Axis Decomposition

- **Databases queried**: bioRxiv, arXiv, PubMed/PMC, Nature Methods, Nature Biotechnology via `search_web`.
- **Exact search queries**:
  1. `perturbation prediction linear baseline PC1 shared axis mean shift Systema Ahlmann-Eltze Miller 2025 2026`
  2. `"Miller" "Deep learning perturbation models can outperform baselines on calibrated metrics" Nature Biotechnology 2026 "2604.27646"`
  3. `"Ahlmann-Eltze" "Deep-learning-based gene perturbation effect prediction" Nature Methods 2025 "Systema" Viñas Torné Nature Biotechnology 2025 2026`
- **Time window**: 2024-11 to 2026-11.
- **Matched Prior Art**:
  1. **Ahlmann-Eltze, Huber, & Anders (*Nature Methods*, August 2025)**: *"Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines"* — demonstrated that low-rank PCA linear baselines match or beat complex deep learning models, pointing to shared low-dimensional structure in perturbation responses.
  2. **Viñas Torné, Wiatrak, Brbić et al. (*Nature Biotechnology*, August 2025 / June 2026), *Systema***: *"Systema: a framework for evaluating genetic perturbation response prediction beyond systematic variation"* — showed standard metrics (MSE, global Pearson) are inflated by systematic mean shifts and shared confounding variation, and proposed centroid-removed evaluation.
  3. **Miller et al. (*Nature Biotechnology*, 2026; `bioRxiv 10.1101/2025.10.20.683304`; `arXiv:2604.27646`)**: *"Deep learning perturbation models can outperform baselines on calibrated metrics"* — evaluated 14 datasets using positive ("interpolated duplicate") and negative controls and Dynamic Range Fraction (DRF), proving that unweighted $\Delta$-Pearson and MSE are poorly calibrated, model rankings invert across metrics, and frozen foundation-model embedding probes are strong competitors.
  4. **Stable-Shift / scDFM (arXiv, 2025–2026)**: Explicitly separate shared latent axes (PC1 / global shift) from gene-specific interaction effects.
- **Verdict**: **部分重叠**.
- **Impact on Manuscript**: Metric miscalibration, ranking inversion, and systematic mean-shift inflation are **already established prior art** (Ahlmann-Eltze 2025, Systema 2025/2026, Miller 2026). Our contribution in `R8-2` must be framed narrowly and accurately as a **pre-registered rank-1 PC1 off-axis decomposition** ($S = \operatorname{span}\{1, b, a\}$) quantifying how much of the apparent $\Delta$-Pearson advantage of PCA-linear (`AE`) and frozen `scGPT` probes over `S5` resides specifically along the single dominant post-centroid variation axis across three CRISPRi atlases (`67%–273%` collapse off-axis).

---

### PA-3: Parameter-Free Orthogonal Fusion of Foundation-Model Embeddings and Kernel / Domain Models (`OrthoFusion`)

- **Databases queried**: bioRxiv, arXiv, PubMed/PMC, OUP, CVF via `search_web`.
- **Exact search queries**:
  1. `"orthogonal" fusion OR projection foundation model linear baseline single-cell perturbation prediction scGPT 2024 2025 2026`
  2. `"scLAMBDA" OR "DECANT" OR "StateXDiff" OR "HyperMap" OR "GRAFT" OR "TranScouter" perturbation prediction bioRxiv 2025 2026`
- **Time window**: 2024-11 to 2026-11.
- **Matched Prior Art**:
  1. **Wang et al. (bioRxiv Dec 2024 / 2025–2026), *scLAMBDA* & *scLAMBDA-Program***: *"Modeling and predicting single-cell multi-gene perturbation responses with scLAMBDA"* — integrates LLM gene embeddings with transcriptomic representations via explicit disentanglement of basal cell state and perturbation-specific components.
  2. **Orthogonal transformation / Cayley-projection disentanglement models (2025–2026, arXiv / OUP)** and **scGenePT**: Enforce mathematical orthogonality between basal state, transcriptomic tokens, and external LLM/prior-knowledge embeddings to prevent modality collapse.
  3. **StateXDiff (arXiv, May 2026)**: Multimodal latent diffusion combining transcriptomic and protein-level features with mechanism-aware templates.
- **Verdict**: **部分重叠** (and **empirically UNCONFIRMED** on an unseen 4th dataset because both Feng 2025 and Jiang 2025 failed pre-registered sample-size gates).
- **Impact on Manuscript**: Orthogonal disentanglement/fusion is an active existing design pattern. Furthermore, `OrthoFusion_scGPT_S5` only passed exploratory `H4` (parity/gain on $\Delta$-r and $v_5$) and failed `H5` ($\ge 6/10$ splits CI $> 0$ on `disc_raw` over both parents) on seen datasets, with 4th-dataset confirmation `PENDING`. It must be presented strictly as an **exploratory, unconfirmed geometric test**, never as a headline novelty or proven SOTA method.

---

### PA-4: Pre-Registered Zero-Look Confirmation in Perturbation Prediction Benchmarks

- **Databases queried**: bioRxiv, arXiv, Nature Methods, ICML, CZI Virtual Cells, Google Scholar via `search_web`.
- **Exact search queries**:
  1. `"pre-registered" OR "preregistered" benchmark single-cell perturbation prediction CRISPR 2024 2025 2026`
  2. `"cross-vendor" OR "cross-model" AI scientist audit "pre-registration" OR "preregistration" autonomous research 2025 2026 arXiv`
  3. `"2608.28631" "CrossAudit" OR "2605.03571" "Preregistration for Experiments with AI Agents"`
- **Time window**: 2024-11 to 2026-11.
- **Matched Prior Art**:
  1. **scPerturBench (*Nature Methods*, 2026)**: Standardized benchmark of 27 single-cell perturbation methods across 29 datasets for unseen contexts and perturbations.
  2. **PerturBench (arXiv, late 2025)**, **OP3 (NeurIPS / 2024)**, **PertEval-scFM (ICML 2025)**, and **CZ Virtual Cells Platform Benchmarks (2026)**: Standardized public benchmarking platforms for single-cell perturbation prediction.
  3. **Vaccaro (arXiv:2606.11217, June 2026)**: *"Preregistration for Experiments with AI Agents"* — formalizes preregistration protocols specifically to curb researcher degrees of freedom and outcome-contingent redesigns in *in silico* AI-agent experiments.
  4. **OSF-preregistered audits of AI Scientist-v2 (NeurIPS 2025/2026 evaluation tracks)**: Use formal pre-registration prior to autonomous execution phases.
- **Verdict**: **部分重叠**.
- **Impact on Manuscript**: While commit-hashed pre-registration with hard-stop dataset rejection (`Feng 2025` / `Jiang 2025` stopped prior to model evaluation) appears uncommon in published single-cell perturbation benchmarks, standardized blind benchmarks (Virtual Cell Challenge, OP3, scPerturBench) and formal AI-agent preregistration (`arXiv:2606.11217`) already exist. We therefore describe our pre-registration as a **strict methodological discipline applied in this study**, NOT as the "first pre-registered perturbation benchmark".

---

### PA-5: Dual-AI Autonomous Research Workflow (Executor AI + Cross-Vendor Auditor AI + Pre-Registration + Public Audit Ledger)

- **Databases queried**: Nature, arXiv, bioRxiv, ICML, NeurIPS, ICLR, Stanford / alphaXiv / Sakana / Google DeepMind / FutureHouse publications via `search_web`.
- **Exact search queries**:
  1. `"Virtual Lab" Swanson Nature 2025 "Scientific Critic" AI scientist`
  2. `"A Society of AI Scientists" OR "Propose-Critique-Falsify" ICML 2026 autonomous research audit`
  3. `"2608.05179" OR "2610.22465" autonomous research agents verification gap reproducibility auditing`
  4. `"Sakana" "AI Scientist" v2 "AI co-scientist" Gottweis FutureHouse "Robin" "Kosmos" "Agents4Science" 2025 2026`
  5. `"A Society of AI Scientists" 2026 arXiv preprint peer review audit`
  6. `"Towards an AI co-scientist" Gottweis arXiv 2025 "Kosmos" Edison Scientific arXiv 2025`
  7. `"Agents4Science" 2025 conference AI scientists Stanford James Zou`
  8. `"cross-vendor" OR "cross-model" AI scientist audit "pre-registration" OR "preregistration" autonomous research 2025 2026 arXiv`
  9. `"2608.28631" "CrossAudit" OR "2605.03571" "Preregistration for Experiments with AI Agents"`
- **Time window**: 2024-11 to 2026-11 (plus 2024 Sakana v1 anchor).
- **Mandatory 9-Reference + Newly Discovered Exact-Overlap Comparison Matrix**:

| # | Prior Work | Citation / Date | Core Mechanism & Scale | Direct Overlap with Our System | Concrete Distinction (Without Overclaiming) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | **Virtual Lab** | Swanson, Zou, Pak et al., *Nature* (Jul/Oct 2025) | Multi-agent lab with AI PI, specialist agents, **Scientific Critic agent**, and human-in-the-loop; wet-lab validated 92 SARS-CoV-2 nanobodies. | Critic agent challenging assumptions and plans inside an AI biology workflow. | Virtual Lab uses single-vendor in-meeting critique + wet-lab validation; our system uses an external cross-vendor auditor re-executing code and blocking claims across 8 git-tracked rounds (computational only, no wet-lab). |
| 2 | **A Society of AI Scientists** | alphaXiv / Stanford / MIT / ETH Zurich et al. (Nov 2026) | 936 autonomous agents across 35 labs over 12 weeks; self-organized journal with **executable code peer review** overturning 77/156 (49.4%) of headline conclusions. | Independent AI agents auditing code/data programmatically and overturning inflated claims before publication. | Vastly larger ecosystem scale (936 preprints); our work is a single deep longitudinal case study (1 project, 8 rounds, 34 logged interventions) showing how repeated cross-vendor re-execution turns an initially overclaimed L1 draft into a pre-registered multi-atlas benchmark—and even catches an auditor's own missing literature search. |
| 3 | **Propose-Critique-Falsify (PCF)** | Bansal & Gangwani, ICML AI4Science Workshop (2026) | Benchmark for multi-agent scientific claim verification (`Propose → Critique → Falsify`); identifies "conservatism-as-abstention" and "self-refinement trap". | Multi-agent adversarial falsification of scientific claims. | PCF benchmarks text/synthetic verification pipelines; our case documents real code-executing falsification (e.g., independent seed sweeps, bit-for-bit reproduction, blind Codabench verification) and confirms that *human/auditor oversight* is still needed to catch upstream literature blind spots. |
| 4 | **Autonomous Research Agents Survey** | Ding, Nannapaneni, Liu, Zhang, `arXiv:2608.05179` (Aug 2026) | Surveys runnable AI scientist systems across 7 audit dimensions; defines the **verification gap** (83% release code vs. only 38% release traces, seeds, or novelty-verification logs). | Directly names the exact verification gap our public audit ledger + `PRIOR_ART_LOG.md` + `AUTONOMY_EVIDENCE.md` address. | Survey/taxonomy; our repository acts as a concrete compliance case study supplying all 7 audited artifacts (including negative gates and retracted claims). |
| 5 | **Agentic Reproducibility Auditing** | Dong, Zhang, Zheng, `arXiv:2610.22465` (Oct 2026) | *"Author-in-the-loop Agentic Reproducibility Auditing for a Machine Learning Conference"* — GitHub two-stage AI audit across 36 conference papers (133 tables/figures), 20 substantive flags resolved. | GitHub-tracked agentic code reproducibility auditing with structured author-auditor feedback loops. | Post-submission conference reproducibility audit of human-authored ML papers; ours is an *in-loop, pre-submission* dual-AI research cycle on a biological benchmark where the executor itself is an AI agent. |
| 6 | **Sakana AI Scientist & v2** | Lu et al. (2024); Yamada et al., `arXiv:2504.08066` (Apr 2025) | End-to-end automated ML research (v2 uses Best-First Tree Search; ICLR 2025 workshop paper). | End-to-end autonomous execution from code to manuscript. | Single-system self-review without cross-vendor adversarial auditing or public retraction ledger; focused on ML experiments rather than external biological atlases. |
| 7 | **Google AI Co-Scientist** | Gottweis et al., `arXiv:2502.18864` (Feb 2025) | Gemini 2.0 multi-agent "generate, debate, evolve" tournament for biomedical hypothesis generation + wet-lab validation (AML, liver fibrosis). | Multi-agent biomedical debate/ranking. | Hypothesis-generation system (human wet-lab execution); does not execute longitudinal multi-round computational auditing against external codebases and blind leaderboards. |
| 8 | **FutureHouse Robin & Edison Kosmos** | Robin (`arXiv:2505.13400`, May 2025); Kosmos (`arXiv:2511.02824`, Nov 2025) | Robin: multi-agent drug discovery (ripasudil for dAMD). Kosmos: 12-hour world-model autonomous data analysis + literature synthesis. | Autonomous multi-hour computational biology data analysis + literature synthesis. | World-model exploration without independent cross-vendor adversarial audit or commit-hashed preregistration failure ledger. |
| 9 | **Agents4Science 2025** | Zou, Queen, Thakkar, Sun, Bianchi et al., Stanford (Oct 2025) | First AI-authored, AI-reviewed conference (47/300 accepted); documented severe hallucination, context drift, and need for guardrails. | Transparent empirical evaluation of AI scientists + AI reviewers. | Conference-level meta-evaluation; our case documents the *micro-dynamics* of why unguarded AI executors overclaim (Rounds 1–2) and what concrete guardrails (independent execution, prior-art gates, hard stops) actually prevent those failures. |
| 10 | **CrossAudit** *(Critical Additional Hit)* | `arXiv:2608.28631` (Aug 2026) | *"CrossAudit: A Git-Native, Cross-Vendor Audit Loop for Agentic Science"* — Generator agent audited at every git increment by an independent cross-vendor Auditor agent against a human-authored rulebook, with deterministic code checks as backstop. | **EXACT ARCHITECTURAL OVERLAP**: Cross-vendor Executor + Auditor + Git tamper-evident ledger + deterministic code checks! | Because `CrossAudit` (`arXiv:2608.28631`) already formalizes the git-native cross-vendor audit architecture, we **CANNOT** claim the architecture itself is novel! Our contribution is strictly an **empirical computational-biology case study** applying a cross-vendor audited loop end-to-end (including blind public Codabench evaluation, three human CRISPRi atlases, a 34-entry unvarnished failure/retraction ledger, and an honest record of a shared auditor-executor failure mode: missing prior-art screening in Round 7). |
| 11 | **AI-Agent Preregistration** *(Additional Hit)* | Vaccaro, `arXiv:2606.11217` (Jun 2026) | *"Preregistration for Experiments with AI Agents"* — formal preregistration templates to prevent hidden iteration and outcome-contingent redesigns. | Directly covers preregistration of AI-agent experimental protocols. | Methodological proposal; our study provides a concrete git-hashed execution trace including two hard-stop dataset rejections (`Feng 2025`, `Jiang 2025`). |

- **Verdict**: **已有 (核心双 AI 跨厂商 Git 审计架构与 AI 预注册范式均已有在先文献)** / **部分重叠 (公开全量错误台账的端到端单细胞扰动生物学案例实证)**.
- **Impact on Manuscript & External Communications**:
  - **STRICTLY PROHIBITED**: "First dual-AI autonomous research system", "first cross-vendor AI audit loop", "首个双 AI 自主科研系统", "首次提出跨厂商 AI 审计", "metric inversion is a novel discovery", "first pre-registered perturbation benchmark".
  - **PERMITTED**: "An empirical case study of a cross-vendor audited autonomous research workflow (following the verification-gap agenda of `arXiv:2608.05179` and cross-vendor git-audit principles of `arXiv:2608.28631`) applied end-to-end to single-cell CRISPRi perturbation modeling, accompanied by a complete 34-entry public audit ledger documenting both executor overclaims and an auditor prior-art oversight."

---

## Prohibited vs. Permitted Phrasing Audit Checklist (Enforced for `R8-10` PDF v1 & All Subsequent Outputs)

1. **Banned Phrases**:
   - "first", "novel discovery", "首次", "首个", "首创" (unless referring solely to a verifiable dataset-specific fact explicitly bounded to *this repository*)
   - "we discover that perturbation metrics invert model rankings" → **REPLACED** with explicit citation to Ahlmann-Eltze (2025), Systema (2025/2026), and Miller et al. (*Nature Biotechnology* 2026; `arXiv:2604.27646`)
   - "OrthoFusion beats all baselines" → **REPLACED** with "exploratory on 2 seen atlases (`H4` met, `H5` not met); 4th-dataset confirmation `PENDING`"
   - "cross-vendor AI auditing is a novel paradigm" → **REPLACED** with explicit citation to `CrossAudit` (`arXiv:2608.28631`), *Virtual Lab* (*Nature* 2025), *A Society of AI Scientists* (2026), PCF (ICML 2026), `arXiv:2608.05179`, and `arXiv:2610.22465`.
2. **Required Disclosures in Manuscript**:
   - Explicitly cite Miller et al. (2026), Ahlmann-Eltze et al. (2025), Systema (2025/2026), STATE (2025), HyperMap (2026), TranScouter (2026), and the autonomous-science prior art (`CrossAudit`, *Virtual Lab*, *A Society of AI Scientists*, `arXiv:2608.05179`, `arXiv:2610.22465`, `arXiv:2606.11217`).
   - Explicitly state in `R8-10` PDF v1 that `R8-1` (native single-cell GEARS) and `R8-6` (Miller DRF calibration) are **in progress** and excluded from v1 claims/abstract until completed.
