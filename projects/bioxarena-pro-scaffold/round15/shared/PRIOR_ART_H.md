# Round 15 Prior-Art Supplement (`PRIOR_ART_H`) — Post-2026-11-16 & VCC 2026 Verification

- **Agent**: A-T5 (Team A · Literature & Official Metrics)
- **Verification Date**: `2026-11-23` (CST)
- **Baseline**: [`round13/status/PRIOR_ART_CONCLUSION_R13.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round13/status/PRIOR_ART_CONCLUSION_R13.md)
- **Sources Swept**: `bioRxiv` / `medRxiv` API (`2026-11-16..2026-11-23`), `arXiv` API (`q-bio.*` & perturbation queries through `2026-11-23`), `EuropePMC` & `PubMed` E-utilities, `Crossref` / `OpenAlex` / `DataCite` / `doi.org`, `Arc Institute` (`arcinstitute.org`, `virtualcellchallenge.org`, `ArcInstitute/cell-eval2`), `GitHub` (`Illumina/CellularAI`, `Xaira-Therapeutics/X-Cell`, and VCC 2026 team repositories), `Ion Genomics`, and `TIME`.

---

## Summary Verdicts

| Question | Focus | Verdict | Key Scope Boundary |
|---|---|---|---|
| **Q1** (vs Claim 1) | VCC 2026 6 lines (`Cal-33`, `HeLa`, `Jurkat`, `AC16`, `MOLM13`, `U-87MG`): target-line controls vs mismatched-line controls | **部分做过** | 仅见：(a) 非 VCC 2026 公开 3' 屏上 `CASTER_wrongctrl` 消融（仅同名覆盖 Nadig `Jurkat`，未测其余 5 株 Arc 2026 Flex 系）；(b) 部分参赛者仅在 **PC 去噪轴/局部上下文特征/同系批次**层面对比他系对照 PC vs 目标系对照 PC；**尚未找到**在 VCC 2026 六株细胞系上做整套预测的"目标系自身对照 vs 错配细胞系对照"系统消融。 |
| **Q2** (vs Claim 2) | Post-`2026-11-16`: "even after independent split-half gene selection, a predictor without the target line's baseline still scores high" | **没有找到** | `2026-11-16` 之后公开的 `OCCAM`、`mixing-vs-selection` v1.3.1、`PerturbJEPA`、`ReliFlow` 等均未报告"独立 split-half 基因选择后无目标系基线预测器仍获高分"。 |
| **Q3** (Dobin Preprint) | Availability of `doi:10.64898/2026.11.17.767214`, definition of `line-specific slope (0.02–0.13 vs 1.02)`, and evaluation of `FCT`-class methods | **未能读取**（原文不可读）；companion 文档与官方博文对公式及 `FCT` 检验均**没有找到** | `doi:10.64898/2026.11.17.767214` 在 Crossref / bioRxiv API / OpenAlex / DataCite / EuropePMC / `doi.org` 均返回 404 或无记录，bioRxiv 网页返回 HTTP 403（**`未能读取`**）。已可读的 `ArcInstitute/cell-eval2` (`vcc2026-metrics.md`) 与 Arc 官方 Wrap-Up 博文未公开该斜率公式，也未专门检验 `FCT` 类方法。 |

---

## Q1: Has Anyone Directly Compared "Target Line's Own Controls" vs "Mismatched Cell-Line Controls" on the VCC 2026 Cell Lines (`Cal-33`, `HeLa`, `Jurkat`, `AC16`, `MOLM13`, `U-87MG`)?

**Verdict: 部分做过（仅限于非 VCC 2026 公开数据集的单一模型对照替换，或 VCC 参赛方案中针对控制主成分去噪轴/特征向量的局部对比；未在 VCC 2026 六株细胞系上做完整预测消融）。**

### Detailed Primary Evidence

1. **Arc Institute VCC 2026 Official Wrap-Up — Winners & Overall Challenge Design**
   - **Verifiable URL**: `https://arcinstitute.org/news/virtual-cell-challenge-2026-wrap-up` (published `2026-11-19`)
   - **Verbatim quote**:
     > "Each team first estimated the effect of every knockdown on pseudobulk expression from measured responses to the same knockdown in other cell types, adapted these estimates to the new cell line using only its control cells, and then generated the 400 cells by modifying the control cells of the line itself."
   - **Finding**: Validation lines were `Cal-33`, `HeLa`, `Jurkat` and test lines were `AC16`, `MOLM13`, `U-87MG`. All three winning teams (`1st BioInsight AI Lab / CellularAI`, `2nd FlowState`, `3rd X-team / X-Cell Transfer`) used each target line's own supplied control cells; Arc explicitly noted that future ablation studies are needed to isolate component contributions and did not compare target-line controls against mismatched cell-line controls.

2. **1st–7th Place Winning Teams (`CellularAI`, `FlowState`, `X-Cell Transfer`, `Obradovic Lab`, `C-team`, `Aginglab.com`, `TAO`)**
   - **Verifiable URLs**:
     - `https://github.com/Illumina/CellularAI` (`docs/evaluation.md`)
     - `https://github.com/Xaira-Therapeutics/X-Cell` (`README.md`)
     - `https://time.com/article/2026/11/20/scientists-are-pursuing-the-dream-of-a-virtual-cell/`
     - `https://www.iongenomics.bio/p/illumina-team-wins-2026-virtual-cell` (`2026-11-20`)
   - **Verbatim quotes**:
     - `CellularAI` (`docs/evaluation.md`):
       > "CellularAI ships one public evaluation harness, **`vcc2025`**, built from the Virtual Cell Challenge 2025 H1 hESC CRISPRi data."
     - `X-Cell` (`README.md`):
       > "Trained on **X-Atlas/Pisces** (25.6M perturbed single cells, 7 CRISPRi Perturb-seq screens), X-Cell integrates multi-modal biological priors through cross-attention and generalizes zero-shot to unseen cell types and perturbations."
     - `TIME` (`2026-11-20`, reporting on FlowState / winners):
       > "While the models did well at predicting the general picture, they did less well at predicting the individual quirks of the different cell types, a key part of a virtual cell model."
   - **Finding**: `CellularAI`'s only open evaluation harness evaluates on `vcc2025` (`H1 hESC`). Neither `Illumina/CellularAI`, `Xaira-Therapeutics/X-Cell`, `FlowState` (Khodaveisi & Bagheri), nor public coverage/listings for ranks 4–7 (`Obradovic Lab`, `C-team / The_Way57`, `Aginglab.com / VirtualOmics`, `TAO`) contains a mismatched-control ablation on `Cal-33`, `HeLa`, `Jurkat`, `AC16`, `MOLM13`, or `U-87MG`. No standalone post-`2026-11-16` repositories or preprints from `FlowState` or ranks 4–7 were indexed across GitHub, bioRxiv, arXiv, EuropePMC, or PubMed.

3. **Kawsar et al. (`CanPert` / `CASTER`, bioRxiv `2026-11-17`) — `CASTER_wrongctrl` Ablation**
   - **Verifiable DOI / URL**: `https://doi.org/10.64898/2026.11.15.766815` (bioRxiv API `date: 2026-11-17`, v1; HTML/PDF HTTP 403 **`未能读取`**) & companion archive `https://doi.org/10.5281/zenodo.24239385` (`CanPert` code repository).
   - **Verbatim quote** (bioRxiv API abstract):
     > "CASTER led only at 2.5% screened (0.62 versus 0.47) and fell below the transferred mean without a screen."
   - **What they actually did (`CanPert` `src/canpert/results_v1.py` & `calibrated_metrics.csv`)**: Tested `CASTER_wrongctrl` ("CASTER, wrong-line controls") only on the four public essential-gene 3' screens (`K562`, `RPE1` from Replogle et al. 2022; `HepG2`, `Jurkat` from Nadig et al. 2025 `GSE264667`). Only `Jurkat` overlaps by cell-line name with VCC 2026 (from an earlier 3' screen, not Arc 2026 Flex); `Cal-33`, `HeLa`, `AC16`, `MOLM13`, and `U-87MG` were **not** evaluated.

4. **VCC 2026 Open Competitor Writeups: `jam-sudo/vcc2026-writeup` (`2026-11-16`) & `1minute99/Virtual-cell-challenge-2026` (`2026-11-23`)**
   - **Verifiable URLs**:
     - `https://github.com/jam-sudo/vcc2026-writeup/blob/main/README.md` (`updated 2026-11-16`, validation `0.1945` on `A/B/C`)
     - `https://github.com/1minute99/Virtual-cell-challenge-2026/blob/master/docs/05-why-it-does-not-work.md` (`updated 2026-11-23`, validation `0.2545`, rank `129/1,635`)
     - `https://github.com/1minute99/Virtual-cell-challenge-2026/blob/master/docs/01-challenge-and-final-system.md`
   - **Verbatim quotes**:
     - `jam-sudo` (`README.md` §9):
       > "I had tried this with the source lines' control cells and it had failed."
     - `1minute99` (`docs/05-why-it-does-not-work.md` §A5):
       > "The same number of PCs from other lines' controls held 0.0659, from gene-permuted controls 0.0060, and a random subspace 0.0014 (≈ k/G)."
     - `1minute99` (`docs/01-challenge-and-final-system.md` §4):
       > "Its rebuilt inputs and budget choice were sound; the one gap found (−0.030) came from feeding the rehearsal a different draw of control cells, which led to the rule "emit from the controls the organisers ship for that context"."
   - **What they actually did**: Both teams only compared **which line's PCA subspace to project out** during residual control-PC deflation (`jam-sudo`: source lines' control PCs vs target line's control PCs on VCC 2026 validation contexts; `1minute99`: destination control PCs vs other lines' control PCs across 7 public leave-one-line-out panels, ΔPDS `+0.0058 [0.0039, 0.0076]`, plus comparing different draws of H1 controls and scrambled-control bugs in S2/S3). Neither systematically replaced target-line controls with mismatched-line controls as the prediction input/baseline across the six VCC 2026 lines (`Cal-33`, `HeLa`, `Jurkat`, `AC16`, `MOLM13`, `U-87MG`).

5. **Other VCC 2026 Public Participant Repos (`saarufer`, `ZehaoLi666`, `MMatinGerami`, `zhuzhenghao-max`) & NeurIPS 2026 Controls Paper**
   - **Verifiable URLs**:
     - `https://github.com/saarufer/vcc2026-atlas-transfer` (`rank 48/806`, `0.2529` final)
     - `https://openreview.net/forum?id=lnpzb1HvwX` (Fan et al., NeurIPS 2026, *Zero-Shot Cross-Cell-Line Perturbation Prediction Needs Controls*; OpenReview Turnstile gated **`未能读取`**; indexed in `https://github.com/Boom5426/Awesome-Virtual-Cell`)
   - **Verbatim quotes**:
     - `saarufer` (`README.md`):
       > "The only target-side input is the line's own non-targeting control cells."
     - `Awesome-Virtual-Cell` (`README.md` catalog entry for `lnpzb1HvwX`):
       > "- **[Cross-Cell-Line Controls]** `[Perturbation]` `[Evaluation & Measurement]` `[Benchmark]` `[Foundation Model]` Cross-Cell-Line Perturbation Prediction Needs Controls (**NeurIPS 2026**) [[paper](https://openreview.net/forum?id=lnpzb1HvwX)]"
   - **Finding**: All public VCC 2026 participant pipelines took each context's supplied control cells as fixed input without mismatched-control swaps. Based on search metadata, Fan et al. evaluates zero-shot cross-line stress tests on Replogle/Nadig (`K562/RPE1/HepG2/Jurkat`) with target-control copy/shift baselines (`未能读取` full text due to Cloudflare Turnstile), not VCC 2026 lines.

---

## Q2: Any Post-`2026-11-16` Works Reporting "Even After Independent Split-Half Gene Selection, a Predictor Without the Target Line's Baseline Still Scores High"?

**Verdict: 没有找到。**

All post-`2026-11-16` papers and public VCC 2026 analyses audited below either (a) evaluate predictors that explicitly **condition on the target line's own controls**, or (b) expose metric inflation / ceiling failures rather than claiming an unanchored predictor scores high after independent split-half gene selection:

1. **`OCCAM` — Bendidi et al. (`arXiv:2611.21959v1`, `2026-11-16`)**
   - **Verifiable URL**: `https://arxiv.org/abs/2611.21959`
   - **Verbatim quote** (Appendix E.2):
     > "We use separate A and B controls to avoid this source of spurious correlation, while reusing B controls between model evaluation and experimental reproducibility."
   - **Why it does not match Claim 2**: OCCAM introduces an A/B split-half protocol with disjoint A and B controls to prevent shared-control inflation of Pearson $\Delta$ (citing Nicol et al. 2026 and Fukutani et al. 2026) and evaluates on HepG2, K562, and 80 VCC-2025 H1 perturbations. Crucially, **every** model and simple baseline (`Control cell`, `Control mean`, `Mean perturbation`, `Transfer`, `Occam`, `STATE`, etc.) **receives and uses the evaluation line's own `B` controls**; no target-baseline-free predictor is tested, nor is independent split-half gene selection combined with `r_orth_cal` shown to preserve high unanchored scores.

2. **`mixing-vs-selection` v1.3.1 — Zhu & Cui (`2026-11-21`)**
   - **Verifiable DOI / URL**: `https://doi.org/10.5281/zenodo.24223291` / `https://github.com/zhuzhenghao-max/mixing-vs-selection`
   - **Verbatim quote** (`README.md`):
     > "Gene-shuffled truth still scored 0.655 in both tiers, and truth plus noise at four times the response norm scored 1.000."
   - **Why it does not match Claim 2**: Audits leave-one-line-out PDS across Replogle/Nadig/X-Atlas/H1/Tahoe-100M, showing PDS depends heavily on pool composition and target selection rather than claiming high unanchored scores after independent split-half gene selection.

3. **`PerturbJEPA` — Pang et al. (`arXiv:2611.20575v1`, `2026-11-16`)**
   - **Verifiable URL**: `https://arxiv.org/abs/2611.20575`
   - **Verbatim quote** (arXiv abstract):
     > "Trained with separate drug and dose ranking objectives, the model predicts candidate responses from a common control population and matches them to the observed response without cell pairing or expression reconstruction."
   - **Why it does not match Claim 2**: Evaluates drug-dose retrieval within cellular backgrounds on Tahoe-100M from a shared control population; does not test target-baseline-free cross-line scoring or split-half gene selection.

4. **`ReliFlow` — David & Peng (`arXiv:2611.23999v1`, `2026-11-18`)**
   - **Verifiable URL**: `https://arxiv.org/abs/2611.23999`
   - **Verbatim quote** (arXiv abstract):
     > "We introduce a conditional flow-matching model that estimates the reliability of single-gene effects from training-only cross-fitting and uses those estimates to calibrate the additive source distribution."
   - **Why it does not match Claim 2**: Addresses combinatorial perturbation response prediction (Norman double genes, ComboSciPlex drug pairs), not unanchored cross-line scoring under independent split-half gene selection.

5. **`Tarikere (2026-11-18)` — Representation-Aware OOD Model Selection**
   - **Verifiable DOI**: `https://doi.org/10.21203/rs.3.rs-11768707/v1` (Research Square)
   - **Verbatim quote** (EuropePMC abstract):
     > "These results provide exploratory evidence that representation coverage can be useful for characterizing OOD behavior, but they do not establish universal superiority of adaptive selection."
   - **Why it does not match Claim 2**: Single-dataset Adamson OOD model-selection analysis; no split-half gene selection or target-baseline-free evaluation.

6. **`He et al. (`2026-11-17`)` — Cell-State Conditioning**
   - **Verifiable DOI**: `https://doi.org/10.64898/2026.11.16.767041` (bioRxiv webpage HTTP 403 **`未能读取`**, abstract readable via EuropePMC)
   - **Verbatim quote** (EuropePMC abstract):
     > "Two incidental findings apply beyond this model: credible-set membership mislabels over half of putatively condition-specific loci, and gating a ratio statistic on its larger term manufactures signal."
   - **Why it does not match Claim 2**: Studies IFN-$\beta$ monocyte eQTL/expression conditioning and ratio-gating artifacts, not target-baseline-free perturbation predictors.

7. **Post-`2026-11-16` bioRxiv/arXiv Updates & Challenge Commentary (`Cupaioli & Merelli` v3 `2026-11-19`; `UPM` v2 `2026-11-18`; `MitoWorld` `2026-11-18`; `CellChain` v2 `2026-11-16`; `Ion Genomics` `2026-11-20`)**
   - **Verifiable DOI / URLs**:
     - `https://doi.org/10.64898/2026.11.03.763968` (bioRxiv v3 `2026-11-19`, full text HTTP 403 **`未能读取`**, abstract via bioRxiv API)
     - `https://doi.org/10.64898/2026.10.25.761625` (bioRxiv v2 `2026-11-18`, full text HTTP 403 **`未能读取`**, abstract via bioRxiv API)
     - `https://arxiv.org/abs/2611.23953` (`2026-11-18`)
     - `https://www.iongenomics.bio/p/illumina-team-wins-2026-virtual-cell` (`2026-11-20`)
   - **Verbatim quotes**:
     - `Cupaioli & Merelli` (`10.64898/2026.11.03.763968` v3 abstract):
       > "Fine-tuning changed classification accuracy and the size of predicted effects but not their ranking, gene networks or the advantage of a training-mean baseline."
     - `Ion Genomics` (`2026-11-20`):
       > "Today’s models add little beyond what a measurement in another cell line already tells you."
   - **Why they do not match Claim 2**: None examines or reports high calibrated/orthogonalized scoring of target-baseline-free predictors under independent split-half gene selection.

---

## Q3: Is Alexander Dobin's VCC 2026 Metrics Preprint (`doi:10.64898/2026.11.17.767214`) Readable? Does It Define the Line-Specific Slope (`0.02–0.13` vs `1.02`) and Test `FCT`-Class Methods?

**Verdict: `未能读取`（`doi:10.64898/2026.11.17.767214` 原预印本全文与元数据目前均不可读）；在已可读的 Arc 官方 `cell-eval2` 文档及赛后总结中，对 `line-specific slope` 的精确回归定义与 `FCT` 类方法专项检验均`没有找到`。**

### 1. Accessibility Verification of `doi:10.64898/2026.11.17.767214` (`2026-11-23`)
- Checked official APIs and resolver on `2026-11-23 16:33–16:40 CST`:
  - `https://doi.org/10.64898/2026.11.17.767214` & `https://doi.org/api/handles/10.64898/2026.11.17.767214`: `HTTP 404 Not Found`
  - Crossref API (`https://api.crossref.org/works/10.64898/2026.11.17.767214`): `HTTP 404 Not Found`
  - bioRxiv API (`https://api.biorxiv.org/details/biorxiv/10.64898/2026.11.17.767214`) & medRxiv API: `{"messages":[{"status":"no posts found"}]}`
  - OpenAlex API & DataCite API: `HTTP 404 Not Found`
  - EuropePMC (`DOI:10.64898/2026.11.17.767214`): `hitCount: 0`
  - Direct bioRxiv URLs (`https://www.biorxiv.org/content/10.64898/2026.11.17.767214v1` and `.full.pdf`): `HTTP 403 Forbidden` (Cloudflare challenge, **`未能读取`**)
  - PubMed author query (`Dobin A[Author] AND 2026[PDAT]`): returns only 4 earlier 2026 papers (`Tahoe-100M` `10.1016/j.cell.2026.08.035`, `scBaseCount` `10.1016/j.cell.2026.08.025`, `State` `10.1016/j.cell.2026.07.052`, `BINSEQ` `10.1371/journal.pcbi.1014181`).
  - `https://dobin.org/` (`HTTP 200`): contains no link to `767214` or VCC 2026 preprint PDF.
- **Status**: Preprint `doi:10.64898/2026.11.17.767214` itself is **`未能读取`**; neither its mathematical formula for `line-specific slope` nor whether it evaluates `FCT`-class predictors can be verified from the preprint itself.

### 2. What Readable Official Arc / Dobin Sources Actually Show
1. **Arc Institute Winners & Wrap-Up (`2026-11-19`)**
   - **Verifiable URL**: `https://arcinstitute.org/news/virtual-cell-challenge-2026-wrap-up`
   - **Verbatim quote**:
     > "Among the genes that respond significantly to a knockdown in a given line, the cosine similarity between the predicted responses to a target in different lines ranged from 0.77 to 0.89, compared with 0.55 for the measured responses and 0.91 to 0.95 for two independent samples of cells from the same line, and the line-specific component of the response that the models did predict was small (a slope of 0.02 to 0.13, compared with 1.02 for a replicate)."
   - **Assessment**: This summary paragraph is the only readable public statement mentioning the `0.02 to 0.13 vs 1.02` slope (for top-10 entries on validation lines, restricted to genes responding significantly in a given line), but it gives no mathematical specification of how the line-specific component/slope is computed, and does not report evaluating unanchored `FCT`-class predictors.

2. **`ArcInstitute/cell-eval2` Official VCC 2026 Metric Specification (`2026/08/19`, maintained by Alexander Dobin)**
   - **Verifiable URL**: `https://github.com/ArcInstitute/cell-eval2/blob/main/docs/vcc2026_metrics/vcc2026-metrics-brief.md` (and companion `docs/vcc2026_metrics/vcc2026-metrics.md`, `src/cell_eval2/configs/vcc2026.yaml`, `src/cell_eval2/baseline.py`)
   - **Verbatim quote** (`vcc2026-metrics-brief.md`):
     > "The baseline $b$ is the context's mean perturbation response: an equal-weight average, over the constructs passing a cell-count and knockdown-efficiency filter, of the per-construct mean count vector, assigned identically to every perturbation and scored as an ordinary submission."
   - **Assessment**:
     - **Line-specific slope**: Not defined anywhere in `vcc2026-metrics.md`, `vcc2026-metrics-brief.md`, or `ArcInstitute/cell-eval2` code (full-repository search for `0.02`, `0.13`, `767214`, and `line-specific` returned zero hits; `slope` appears only for the piecewise-linear scoring tail derivative in `docs/metrics.md`, `metrics/delta.py`, and `tests/test_score.py`).
     - **`FCT`-class methods**: `vcc2026-metrics.md` only defines the 6 official competition metrics (`pds_cosine`, `expr_mse_unbiased_capped_norm`, `de_wilcoxon_direction_fidelity_yield_raw`, `de_wilcoxon_direction_reach_raw`, `de_wilcoxon_sig_jaccard`, `de_wilcoxon_lfc_nmae`), the target-context mean-response baseline $b$, and the 5-split replicate anchor $r$; it does **not** test or score cross-line `FCT`-class predictors without target-line controls.
