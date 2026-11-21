# Competition & Preprint Monitor (`bioRxiv` + `arXiv` + `EuropePMC`)

**Maintainer**: Agent A-T5 (Literature & Official Metrics)  
**Baseline Scan Completed**: `2026-11-21T19:00+08:00` (`2026-11-21T11:00Z`)  
**Current Alert Status**: **`CLEAR` (No active scoop on Team A's zero-shot control co-fluctuation linear-response hypothesis for line-specific residuals)**  
**Recurring Cadence**: Every 2 hours (`0 */2 * * *`, daemon schedule)

---

## 1. Scoop Alert Criterion

Trigger `ALERT` immediately at the top of `round9/p2/status/A-T5.txt` and log here if any newly indexed `bioRxiv`, `medRxiv`, `arXiv`, or `EuropePMC` preprint satisfies **all or core overlap** with Team A's falsifiable question:
1. Uses **unperturbed target-line single-cell control gene–gene co-fluctuations / covariance** ($C_t$ or $\Sigma_t$) in **zero-shot cross-cell-line mode** (0 target perturbation labels used), **AND**
2. Predicts or evaluates the **cell-line-specific perturbation residual component** (after removing cross-line average / shared perturbation response) or directly targets the VCC 2026 Wrap-Up `0.02–0.13` line-specific slope ceiling.

Secondary trigger (`WATCH_UPDATE` / metric check):
- Indexing of **Alexander Dobin's 2026 metrics preprint (`10.64898/2026.11.17.767214`)** so we can verify whether Arc's published `line-specific slope` equation matches §2 of [LINE_SPECIFIC_METRIC.md](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2/vcc_linespecific/audit/LIT/LINE_SPECIFIC_METRIC.md) and flag any discrepancy.

---

## 2. Baseline Scan Coverage (Through `2026-11-21`)

| Source | Scope & Query | Records Inspected | Scoop Risk |
|---|---|---|---|
| **`api.biorxiv.org` (Full, all categories)** | Date window `2026-11-16` to `2026-11-21` (pages 1–53; `bioinformatics`: 186, `genomics`: 71, `systems_biology`: 32, unfiltered all categories: `1,572`); plus `2026-11-20..2026-11-21` total check (`263`) | 1,572 total (27 keyword-matched inspected in detail) | **None** |
| **`export.arxiv.org/api/query`** | (a) Nov 2026 `(perturbation OR "virtual cell") AND ("cross-cell" OR "cell line" OR "line-specific" OR "unseen cell" OR "co-fluctuation" OR covariance OR transcriptomic)` (30 hits); (b) Nov 2026 `("Virtual Cell Challenge" OR "VCC 2026" OR "single-cell perturbation" OR "Perturb-seq")` (2 hits); (c) Nov 2026 `(cat:q-bio.* OR "single-cell" OR "Perturb-seq") AND (perturbation OR "virtual cell")` (30 hits); (d) 2026-08..11 targeted cross-line/line-specific/VCC query (1 hit); (e) All-time `"Virtual Cell Challenge" OR "cell-eval2" OR (au:Dobin AND cell)` (5 hits) and `"co-fluctuation" AND (perturbation OR "gene expression" OR CRISPR)` (1 hit) | 69 query hits deduplicated | **None** |
| **EuropePMC REST API (`PPR`, `PMC`, `MED`)** | Oct–Nov 2026 queries for `"Virtual Cell Challenge" OR "cell-eval2"`, `"line-specific" OR "cross-cell-line"`, `"co-fluctuation" OR "CIPHER"`, and `AUTH:"Dobin"`; plus Crossref/OpenAlex/PubMed check on `10.64898/2026.11.17.767214` | 46 records + 1 DOI watch | **None** (1 new Oct 2026 VCC transfer preprint flagged below) |

---

## 3. Monitored Recent Preprints & Relevance Triage

### 3.1 Priority Watch Item (Metric Formula Confirmation)
- **`bioRxiv 10.64898/2026.11.17.767214`** — Alexander Dobin (*VCC 2026 metrics preprint*, cited by Arc Wrap-Up 2026-11-19).
  - **Baseline Status (`2026-11-21`)**: Unindexed in `api.biorxiv.org` (`{"status":"no posts found"}`), Crossref (`Resource not found`), OpenAlex, EuropePMC (`hitCount: 0`), and PubMed. `www.biorxiv.org` HTML gated by HTTP 403 Cloudflare (not bypassed). Tracked every 2 hours.

### 3.2 Recent Cross-Cell-Line / Perturbation Preprints (Oct–Nov 2026, Verified Non-Scoops)

| Date | Identifier | Title & Authors | Method Summary | Overlap / Risk Assessment |
|---|---|---|---|---|
| `2026-11-19` | `bioRxiv 10.64898/2026.11.03.763968` | *Where single-cell foundation models help: a public-data benchmark...* (Cupaioli & Merelli; `fmaudit`) | Benchmarks 5 FMs (`scGPT`, `Geneformer V1/V2`, `Nicheformer`, `scFoundation`, `CellOracle`, `STATE`) vs simple baselines | **Supporting context**: Zero-shot FMs fail 46/46 knockout tests vs expression-matched controls. **No scoop.** |
| `2026-11-18` | `bioRxiv 10.64898/2026.10.25.761625` | *How Much Can AI Virtual Cells Gain from Unified Perturbational Modeling?* (`UPM`; Wei et al., `youlab2025/UPM`) | >1B-param VAE trained on 89.7M cells across 71 contexts with gene-channel alignment & shared perturbation encoding; evaluated on 5 held-out Tahoe-100M lines vs `State`/`Stack` | **Adjacent competitor**: Massive supervised multi-context VAE; does not use target-line control co-fluctuation linear response $C_t$. **No scoop.** |
| `2026-11-18` | `bioRxiv 10.64898/2026.10.20.760537` | *Averaging in gene perturbation benchmarks hides perturbation-strength variation* (Sharma et al.) | Audits perturbation-strength variation via guide abundance in X-Atlas/Orion CRISPRi (`HCT116`, `HEK293T`) and retrieval in `STATE` / `Cell2Sentence` | **Benchmark audit**: Evaluates within-perturbation heterogeneity; not cross-line prediction. **No scoop.** |
| `2026-11-18` | `arXiv:2611.23999v1` | *ReliFlow: Reliability-Calibrated Residual Transport for Combinatorial Single-Cell Perturbation Prediction* (David & Peng) | Cross-fitted single-effect reliability + conditional flow matching residual transport on Norman / ComboSciPlex | **Combinatorial only**: Not cross-cell-line or control covariance. **No scoop.** |
| `2026-11-17` | `bioRxiv 10.64898/2026.11.15.766815` | *Transferring CRISPRi knockdown responses to new cell lines requires screening part of each line* (`CASTER`; Kawsar et al.) | Set Transformer over 32 control cells + `PerturbMean` prior across 4 lines | **Audited in `PRIOR_ART_A.md`**: Fails below `PerturbMean` at `0%` screened; no gene–gene covariance $C_t$. **No scoop.** |
| `2026-11-16` | `arXiv:2611.21959v1` | *OCCAM: Improved Single-Cell Representation Learning and Perturbation Effect Prediction* (Bendidi et al.) | Denoising MLP + gene interaction + text + protein sequence features; evaluated on cross-context few-shot & zero-shot held-out lines | **Adjacent competitor**: Neural MLP + external priors; no unperturbed covariance $C_t$ or residual slope decomposition. **No scoop.** |
| `2026-11-16` | `bioRxiv 10.64898/2026.09.24.754259` | *PerturbBridge: Conditional Latent Schrödinger Bridge...* (Liu et al., `echQoQ/PerturbBridge`) | Response-anchored autoencoder + Schrödinger Bridge transformer on Norman, Sci-Plex3, Replogle RPE1 | **Generative transport**: Not zero-shot cross-line control covariance. **No scoop.** |
| `2026-11-13` | `arXiv:2611.17705v1` | *S1Pert: Scalable Target-Specific Gene Response Prediction with System-One Intelligence* (Sheng et al., `liujiawen-jpg/S1Pert`) | Frozen System-One encoder + cached pairwise interactions + hierarchical cross-attention decoder | **Supervised architecture**: No $C_t$ co-fluctuation linear response. **No scoop.** |
| `2026-11-08` | `arXiv:2611.11217v1` | *scBIG++: Learning Across Genes and Cells...* (Ruan et al., `ttruan2426-dot/scBIG`) | Gene-Relation Clustering + Dynamic SoftSub + Set-GCAE + latent flow across unseen targets, combinations, & cross-cell-line transfer | **Deep generative flow**: Does not test CIPHER control co-fluctuations on line-specific residuals. **No scoop.** |
| `2026-11-07` | `arXiv:2611.10658v1` | *The Delta Gap: Cellular World Models for Perturbation Response Prediction Are Far Harder Than Global Correlation Implies* (Junhui Gao) | Quantifies global vs delta Pearson (`GEARS` `0.977` vs `0.076` on Replogle K562); adds post-hoc RWR on co-expression graphs + kNN; evaluates K562 $\to$ Norman A549 (741 genes) | **Closely adjacent methodology**: Uses co-expression RWR on perturbation geometry / delta-aligned training rather than unperturbed target-line control covariance $C_t$ for line-specific residuals. **No scoop (monitor citation).** |
| `2026-11-05` | `arXiv:2611.14559v1` | *Predicting transcriptional interactions with GrnPFN* (Haihui Zhang) | Prior-data fitted network pretrained on simulated GRNs for combinatorial interactions | **Combinatorial only**: **No scoop.** |
| `2026-11-04` | `arXiv:2611.05849v1` | *A low-rank decoder bottleneck bounds reliability in foundation-model perturbation prediction* (Zhang et al.) | Audits `STATE` ST decoder (`SE-600M`, `scGPT`, `Tahoe-x1`, 4 CRISPRi lines): effective rank 4–32 bounds OOD accuracy via NN training distance | **Supporting mechanistic diagnosis**: Explains why FM decoders collapse toward training-average responses; no $C_t$ predictor. **No scoop.** |
| `2026-10-25` | `bioRxiv 10.64898/2026.10.20.760514` | *Centered cross-context perturbation transfer on the public 2025 Virtual Cell Challenge benchmark* (Luan, Chen, Yu & Yang; SJTU) | Combines public perturbation summaries (`KOLF2.1J` centered on `H1` training targets) + DE kernel fitted to `H1` training targets + isotonic regression (`TEST100` composite `0.341344` vs champion `0.340052`) | **Adjacent VCC cross-context transfer baseline**: Fits DE kernel on target-line (`H1`) training targets and transfers centered `KOLF2.1J` ranks; does **not** use target control co-fluctuation $C_t$ zero-shot. **No scoop.** |

*(Other non-perturbation / domain-unrelated hits among the 1,572 Nov 16–21 bioRxiv records—e.g. `Cryo-EM`, `MetaboRx`, `Infernal-GPU`, `Elflow`, `ARC-Impute`, `Drug-Prot`, `ImmunoTwin`, `CellVADIS`, `DRAFT`, `Circulating/clinical`—were reviewed and excluded.)*

---

## 4. Incremental Monitor Run Log

| Run Timestamp (UTC+8) | Sources Checked | New Matches Evaluated | Alert Status | Notes |
|---|---|---|---|---|
| `2026-11-21 19:00` | `bioRxiv` (`2026-11-16..2026-11-21`, 1,572 preprints), `arXiv` (`2026-11-01..2026-11-21` + targeted all-time), `EuropePMC` (Oct–Nov 2026), `Dobin DOI` watch | 13 perturbation/VCC records triaged (see §3) | **`CLEAR`** | Baseline scan complete. `10.64898/2026.11.17.767214` unindexed as of `19:00+08:00`. |
| `2026-11-21 20:17` | `bioRxiv` (`2026-11-19..2026-11-21`, 564 preprints), `arXiv` recent 30, `Dobin DOI` watch | 0 genuinely new (`10.64898/2026.11.18.767605` ImmunoTwin already triaged) | **`CLEAR`** | Cron iteration 1 (`run_monitor_scan.py 2`). `10.64898/2026.11.17.767214` still unindexed across bioRxiv/Crossref/EuropePMC. |
| `2026-11-21 22:18` | `bioRxiv` (`2026-11-19..2026-11-21`, 566 preprints), `arXiv` recent 30, `Dobin DOI` watch | 0 new matches (`new_biorxiv_matches: []`, `new_arxiv_matches: []`) | **`CLEAR`** | Cron iteration 2 (`run_monitor_scan.py 2`). `10.64898/2026.11.17.767214` still unindexed across bioRxiv/Crossref/EuropePMC/OpenAlex/PubMed. |
