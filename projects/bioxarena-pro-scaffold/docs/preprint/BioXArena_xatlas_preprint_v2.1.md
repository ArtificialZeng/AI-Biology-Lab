---
title: "Pre-Registered Three-Atlas Confirmation, Principal-Axis Decomposition, and Cross-Line Kernel Modeling Under Metric-Dependent Perturbation Rankings"
author:
  - "Zijian Zeng"
institute: "AI Biology Lab (github.com/ArtificialZeng)"
date: "November 21, 2026 (Preprint v2.1)"
geometry: "margin=0.85in"
fontsize: 10.5pt
linestretch: 1.15
colorlinks: true
linkcolor: NavyBlue
citecolor: NavyBlue
urlcolor: NavyBlue
header-includes:
  - \usepackage{booktabs}
  - \usepackage{float}
  - \usepackage{caption}
  - \captionsetup{font=small,labelfont=bf}
  - \usepackage{amsmath,amssymb}
  - \usepackage{microtype}
  - \usepackage{titling}
  - \pretitle{\begin{center}\Large\bfseries}
  - \posttitle{\end{center}\vspace{0.2em}}
  - \preauthor{\begin{center}\normalsize}
  - \postauthor{\\ \small AI Biology Lab (\texttt{github.com/ArtificialZeng})\end{center}\vspace{-0.4em}}
  - \predate{\begin{center}\small}
  - \postdate{\end{center}\vspace{0.5em}}
  - \let\origfigure\figure
  - \let\endorigfigure\endfigure
  - \renewenvironment{figure}[1][2]{\expandafter\origfigure\expandafter[H]}{\endorigfigure}
---

## Abstract

Recent perturbation-modeling benchmarks established that uncalibrated correlation and mean squared error (MSE) metrics poorly separate biological signal from mean baselines, that model rankings shift across metrics, and that frozen foundation-model probes are strong competitors (Miller et al., 2026; Mao et al., 2026). Across three human CRISPRi atlases (xatlas-orion-perturb, Replogle et al. 2022, and unseen Nadig et al. 2025), we provide: (1) pre-registered zero-look confirmation of ranking reversal under simultaneous union holdout; (2) rank-1 principal-axis (PC1) decomposition explaining linear models' post-centroid advantage; and (3) a frozen cross-line kernel baseline (S5) evaluated under Miller Dynamic Range Fraction ($\mathrm{DRF}_{\mathrm{raw}}$) calibration. On xatlas-orion-perturb, one frozen Codabench submission scored 0.0998 (+0.0173 above the best 2-hour agent at 0.0825, within the 0.02 noise margin; produced by multi-round iteration of a dual-AI executor–auditor system under a human-set goal, not comparable to single-shot 2-hour autonomous agents). On unseen Nadig 2025, H1 ($\mathrm{disc}_{\mathrm{raw}}$ +0.0487, 9/10 CI > 0), H2 ($\Delta r$ −0.0195, 9/10 negative, 4/10 CI < 0), and H3 ($v_5$ +0.0595, 10/10 CI > 0) confirmed. Removing training PC1 cuts AE−S5 on $\Delta r$ by 67%–273% across all three atlases. Under uniform $\mathrm{DRF}_{\mathrm{raw}}$ calibration, MSE is poorly calibrated (0.035–0.063), $v_5$ is only moderately calibrated (0.38–0.42, below $\Delta r$ at 0.48–0.51), and $\mathrm{de50\_}r$, $\mathrm{de50\_sign}$, $\mathrm{disc}_{\mathrm{raw}}$, and $\mathrm{disc}_{v_5}$ are best calibrated (0.81–0.98). Consequently, S5's advantage over the 33M-cell frozen scGPT probe is confined to a moderately calibrated metric ($v_5$); on the best-calibrated metrics, standalone S5 ties ($\mathrm{disc}_{\mathrm{raw}}$) or trails ($\mathrm{de50\_sign}$, $\mathrm{de50\_}r$) the scGPT probe.

---

## 1. Introduction

Predicting transcriptome-wide responses to unseen single-gene CRISPR perturbations is a central benchmark for virtual-cell modeling. Systematic evaluations showed that complex deep learning architectures often failed to outperform simple linear or perturbation-mean baselines under global correlation metrics (Ahlmann-Eltze et al., 2025; Viñas Torné et al., 2025). In parallel, cross-context transfer frameworks have addressed translating perturbation signatures across cell types and donors: STATE trains a large multi-context transformer across observational and perturbed single-cell atlases (Adduri et al., 2025); TranScouter studies lightweight cross-condition encoder-decoders (Zhu and Li, 2026); GRAFT anchors cross-cell-type drug prediction on genetic perturbations (Luo et al., 2026); and—most closely related to our cross-line bridging method S5—**HyperMap** introduces a parameter-efficient meta-learning framework that translates existing perturbation atlases to unseen cell lines, induced pluripotent stem cell (iPSC) donors, and compounds using a sparse set of shared perturbation "seeds" (Dhaka, Gao, and Ideker, 2026).

Two recent benchmark studies further demonstrated that evaluation metric choice fundamentally governs model rankings. Miller et al. (2026) introduced a positive-control (duplicate/interpolated duplicate) and negative-control (mean) calibration framework quantified by the Dynamic Range Fraction (DRF), demonstrating across 14 datasets that MSE and unweighted $\Delta$-Pearson are poorly calibrated, that switching to calibrated metrics alters model rankings, and that linear probes on frozen foundation-model embeddings perform strongly. Independently, Mao et al. (2026) documented across virtual-cell benchmarks that different evaluation metrics yield substantially different model rankings. Accordingly, we do not claim metric-dependent ranking reversal or seed-based cross-cell-line transfer as new discoveries. Building directly on Miller et al. (2026), Viñas Torné et al. (2025), Ahlmann-Eltze et al. (2025), Mao et al. (2026), and Dhaka et al. (2026), this study provides four concrete contributions:

1. **Pre-registered, zero-look three-atlas confirmation**: Locking hypotheses, preprocessing, and frozen solver constants in git before inspecting unseen Nadig et al. (2025) expression values, we test metric-dependent ranking reversal across three laboratories under strict simultaneous union holdout and 4,000-resample cell-line-stratified row-level paired bootstrap inference.
2. **Direct rank-1 PC1 off-axis decomposition beyond centroid subtraction**: Systema (Viñas Torné et al., 2025) showed that subtracting the mean perturbation centroid removes the constant shared shift. We prove algebraically that training-centroid subtraction alone ($\Delta r \equiv \mathrm{centroid\_train\_}r$) still leaves linear baselines ahead of cross-line residual kernels, and confirm via pre-registered orthogonal projection across all three atlases that this residual linear advantage resides predominantly along the single dominant shared training residual axis (PC1).
3. **Closed-form cross-line kernel baseline (S5), Miller DRF calibration, and exploratory orthogonal fusion**: Unlike HyperMap's multi-atlas neural meta-learner over 18 contexts (Dhaka et al., 2026), S5 is a closed-form, zero-neural-parameter, source-centered dual-ridge kernel translator trained strictly on a single paired cell-line atlas with frozen constants transferred zero-shot across laboratories. Using uniform $\mathrm{DRF}_{\mathrm{raw}}$ calibration, we show transparently that S5's sole clear advantage over a 33M-cell frozen scGPT probe is confined to a moderately calibrated metric ($v_5$, $\mathrm{DRF}_{\mathrm{raw}}$ 0.38–0.42), whereas on the highest-DRF metrics (0.81–0.98) standalone S5 ties ($\mathrm{disc}_{\mathrm{raw}}$) or trails ($\mathrm{de50\_sign}$, $\mathrm{de50\_}r$) the scGPT probe. We also report an exploratory parameter-free orthogonal subspace fusion alongside two pre-registered negative screening gates on underpowered fourth datasets.
4. **Audited empirical case study of dual-AI research governance**: Recent metascience works have defined the verification gap in autonomous research agents (Ding et al., 2026), formalized git-native cross-vendor auditing (Dong and Chen, 2026), and evaluated multi-agent critique, reproducibility auditing, and AI-agent preregistration (Swanson et al., 2025; Ahmad et al., 2026; Dong, Zhang, and Zheng, 2026; Vaccaro, 2026; ICML, 2026). This repository supplies an end-to-end computational biology trace accompanied by a 34-entry unvarnished audit ledger, an API-verified reference audit, and an autonomy/human-intervention log.

We evaluate single-gene perturbation prediction across paired human cell lines $(A, B)$ across seven metrics spanning reconstruction error ($\mathrm{MSE}$), centroid-centered and top-differential-expression reconstruction ($\Delta r \equiv \mathrm{centroid\_train\_}r$, $\mathrm{centroid\_val\_}r$, $\mathrm{de50\_}r$, $\mathrm{de50\_sign}$), orthogonal-residual correlation ($v_5$, Pearson correlation on non-target genes after QR-projecting orthogonal to $\mathrm{span}\{\mathbf{1}, \mathbf{b}_c, \mathbf{a}_c\}$), and tie-aware perturbation discrimination percentiles ($\mathrm{disc}_{\mathrm{raw}}$ in unprojected space and $\mathrm{disc}_{v_5}$ in projected space):
$$\mathrm{disc}_{\mathrm{raw}}(i) = 1 - \frac{1}{M - 1}\!\left[\sum_{j=1}^M \mathbb{I}(c_{ij} > c_{ii}) + \tfrac{1}{2}\!\left(\sum_{j=1}^M \mathbb{I}(c_{ij} = c_{ii}) - 1\right)\right] \in [0, 1],$$
where $c_{ij} = \mathrm{corr}(\hat{y}_i - \mathbf{b}_c,\, y_j - \mathbf{b}_c)$, 0.50 denotes random retrieval (including constant-mean ties), and 1.00 denotes exact top-1 retrieval.

---

## 2. Methods

### 2.1 Core Architecture: Soft-Power HybridTOM Kernel and Cross-Line Perturb-Bridge (S5, Primary Method)

All solver constants were permanently frozen at git commit 3609811 on BioXArena-Pro xatlas-orion-perturb prior to external atlas evaluation. Given paired perturbation matrices across two cell lines $(A, B)$ sharing readout genes $\mathcal{G}_{\mathrm{sh}}$, primary method S5 executes in closed form without gradient descent:

1. **Leak-free training subspace**: Within each split's training set for cell line $c \in \{A, B\}$, set baseline $\mathbf{b}_c$ to the mean of reproducible training perturbations and shared axis $\mathbf{a}_c$ to the unit-norm leading right singular vector of centered reproducible training residuals. Project training targets orthogonal to $\mathrm{span}\{\mathbf{1}, \mathbf{b}_c, \mathbf{a}_c\}$.
2. **Self-masked, precision-weighted Soft-Power HybridTOM kernel (S1)**: Form baseline-subtracted, $\sqrt{n_{\mathrm{cells}}}$-weighted training matrices with each perturbed target gene's self-knockdown entry masked to zero. Compute gene-gene cosine similarity $G^{(0)}$ across training perturbations, then construct the signed soft-power HybridTOM kernel ($p = 1.5$):
$$K_{\mathrm{HybridTOM}}(G; p) = \frac{1}{2}\!\left(\mathrm{sign}(G)|G|^p + \mathrm{TOM}\big(\mathrm{sign}(G)|G|^p\big)\right),$$
combining direct gene-gene co-variation with higher-order topological neighborhood overlap.
3. **Multi-scale spectral deflation (S2) and pan-human cross-line gene graph (S3)**: Deflate the top $k \in \{0, 1, 2, 5, 10\}$ principal components of each line's training matrix (weights 1.0, 0.8, 0.5, 0.3, 0.2) and fuse the partner cell line's shared-gene correlation pyramid (weights 0.6, 0.8, 0.8, 0.6, 0.5).
4. **Empirical-SNR non-reproducible anchors (S4)**: Retain non-reproducible measured training perturbations with down-weighted signal-to-noise weights $w_{\mathrm{nr}} = 0.25 \cdot \mathrm{clip}\big((z_{\mathrm{eff}} - 0.9\,\mathrm{med})/(0.3\,\mathrm{med}),\, 0.1,\, 1.5\big)$.
5. **GeneSlope-CosKRR cross-line translator and symmetric bridge (S5)**: Fit a bidirectional nonlinear cross-line translator (skip weight 0.8, ridge $\alpha = 2.0$, kernel power $p_k = 1.5$) on shared training perturbations measured reproducibly in both lines, combining gene-wise shrunk slope skips on shared genes with signed-power cosine kernel ridge regression on 5,000-gene residual profiles. The translator projects partner-only reproducible training perturbations into the target cell line (weight 0.6) and symmetrically fuses translated partner-line KRR predictions ($\alpha = 0.5$) of held-out queries (bridge weight 0.9).

### 2.2 Secondary Extensions (S6, S7) and Exploratory Parameter-Free Orthogonal Fusion

- **Stage S6 (Reciprocal Response Borrowing)**: Borrows residual responses from training perturbations that strongly perturb the query gene ($|z| > 4$, quadratic ReLU weight). Although S6 attains the highest $\mathrm{disc}_{\mathrm{raw}}$ across all three atlases, its orthogonal-residual behavior is dataset-specific (+0.0130 on xatlas vs. −0.0153 on Replogle) and cannot rule out CRISPRi *cis*-proximity confounding on anonymized xatlas.
- **Stage S7 (Output SVD Wiener Denoising and Softmax Manifold Snapping)**: Applies spectral Wiener shrinkage ($\tau = 0.020$) and cosine-similarity manifold snapping (weight 0.45, temperature 18.0, top 40 neighbors). S7 is the exact frozen Codabench submission pipeline; it consistently increases $v_5$ correlation while reducing $\mathrm{disc}_{\mathrm{raw}}$ and inflating MSE.
- **Exploratory Parameter-Free Orthogonal Fusion (OrthoFusion)**: Let $Q_c \in \mathbb{R}^{G \times k}$ ($k \le 3$) be the orthonormal QR basis of $S_c = \mathrm{span}\{\mathbf{1}, \mathbf{b}_c, \mathbf{a}_c\}$, with $P_{S_c}(\mathbf{u}) = Q_c(Q_c^\top \mathbf{u})$ and $P_{\perp,c}(\mathbf{u}) = \mathbf{u} - P_{S_c}(\mathbf{u})$. With zero tunable weights:
$$\hat{\mathbf{y}}_{\mathrm{fusion}} = \mathbf{b}_c + P_{S_c}\!\big(\hat{\mathbf{y}}_{\mathrm{scGPT}} - \mathbf{b}_c\big) + P_{\perp,c}\!\big(\hat{\mathbf{y}}_{\mathrm{S5}} - \mathbf{b}_c\big),$$
assigning the shared subspace to the frozen scGPT embedding ridge probe and the orthogonal complement to frozen S5.

### 2.3 Datasets, Pre-Registered Preprocessing, Baselines, DRF Calibration, and Inference

1. **Development atlas (BioXArena-Pro xatlas-orion-perturb)**: HCT116 and HEK293T cell lines (5,000 readout genes per line, 4,040 shared; anonymized gene identifiers never de-anonymized), evaluated on 10 development splits (seeds 42, 0, 1, 2, 3, 7, 11, 99, 2024, 31337), 10 unseen resplits (seeds 101–110), and one frozen Codabench private test submission (455 targets).
2. **External zero-shot atlas 1 (Replogle et al., 2022)**: RPE1 and K562 genome-wide CRISPRi bulk normalized datasets (Figshare+ DOI 10.25452/figshare.plus.20029387), top 5,000 finite genes by variance-to-mean ratio, reproducibility defined by Anderson-Darling count $\ge 50$ and energy-test $p < 0.001$, deduplicated to one row per Ensembl target.
3. **Pre-registered unseen confirmatory atlas 2 (Nadig et al., 2025)**: HepG2 (145,473 cells $\times$ 9,624 genes) and Jurkat (262,956 cells $\times$ 8,882 genes) raw single-cell Perturb-seq matrices (NCBI GEO GSE264667, PMID 40259084). Preprocessing was locked before inspecting expression values: per-GEM-group library-size normalization and z-scoring against non-targeting controls matching Replogle et al. (2022), top 5,000 finite genes by control mean (3,644 shared Ensembl genes), vectorized Scholz-Stephens two-sample Anderson-Darling test (count $\ge 50$ at BH FDR $< 0.05$, unit-verified within $10^{-5}$ vs. SciPy) plus 20-PC 2,000-permutation energy-distance test ($p < 0.001$) with $n_{\mathrm{cells}} \ge 25$, yielding 2,343 HepG2 (1,000 reproducible) and 2,350 Jurkat (1,262 reproducible) deduplicated Ensembl targets.
4. **GEARS (official package, pseudobulk-adapted)**: Official PyPI package cell-gears 0.1.2 run with default hyperparameters (hidden size 64, 1 GO GNN layer, 1 gene GNN layer, decoder hidden size 16, 20 epochs, learning rate $10^{-3}$, weight decay $5 \times 10^{-4}$, batch size 32) independently per cell line on each split's training set (90%/10% internal early stopping, zero test leakage). Because this completed evaluation supplied one pseudobulk profile per perturbation, 10 jittered copies ($\mathcal{N}(0, 10^{-4})$) of the training control centroid as control cells, and train-only pre-populated top-20 differential gene lists with differential-expression calculation skipped rather than native single-cell count distributions, all outputs are labeled **GEARS (official package, pseudobulk-adapted)**.
5. **Frozen scGPT gene-embedding probe**: Official scGPT human checkpoint (Cui et al., 2024; 33M normal human cells, MIT License) 512-dimensional frozen gene token embeddings paired with ridge regression ($\alpha = 10.0$), evaluated across all 10 splits on Replogle 2022 and Nadig 2025.
6. **Miller Dynamic Range Fraction (DRF) metric calibration**: Following Miller et al. (2026), for each metric we compute negative control $m_{\mathrm{neg}}$ (training reproducible mean baseline) and positive control $m_{\mathrm{pos}}$ under both raw replicates ($\mathrm{DRF}_{\mathrm{raw}}$, primary uniform estimator reported throughout prose) and Miller interpolated duplicates ($\mathrm{DRF}_{\mathrm{interp}}$, blending replicate profiles toward the training mean by gene-wise BH-FDR differential-expression weight $w_g = \mathrm{clip}(1 - q_g, 0, 1)$, reported alongside $\mathrm{DRF}_{\mathrm{raw}}$ in Table 3). For higher-is-better metrics (ceiling 1.0), $\mathrm{DRF} = (m_{\mathrm{pos}} - m_{\mathrm{neg}}) / (1 - m_{\mathrm{neg}})$; for MSE (floor 0.0), $\mathrm{DRF} = (m_{\mathrm{neg}} - m_{\mathrm{pos}}) / m_{\mathrm{neg}}$. On Replogle 2022, positive controls use the independent experimental Day 6 K562 essential replicate screen predicting the Day 8 K562 genome-wide screen (RPE1 has no second normalized bulk replicate screen locally). On Nadig 2025, positive controls use exact Addendum-B normalized 50/50 single-cell split-half pseudobulks (Half A predicting Half B) on both HepG2 and Jurkat.
7. **Strict simultaneous union holdout and row-level paired bootstrap**: For each seed, 20% of unique reproducible in-readout targets are sampled independently per cell line, their union is taken, and all union targets are excluded simultaneously from both lines' training sets. Metrics are macro-averaged equally across cell lines. Because overlapping splits invalidate split-level paired $t$-tests, all paired comparisons report 4,000-resample cell-line-stratified row-level paired bootstrap 95% confidence intervals per split and counts of splits (out of 10) excluding zero.

### 2.4 AI-Assistance and Governance Statement

Zijian Zeng directed this study, set the research goal, made all repository governance and public-release decisions, and bears sole authorship responsibility. Two AI systems assisted under a cross-vendor separation of duties (following Dong and Chen, 2026; Ding et al., 2026): Antigravity (Gemini) implemented modeling code, executed experiments, generated figures, and drafted manuscript text; Claude served as an independent cross-vendor technical auditor reviewing pre-registration protocols, enforcing leakage and statistical boundaries, and re-executing scripts and recomputing metrics across audit checkpoints. During Rounds 1–3, the human manually forwarded Claude's audit critiques to Antigravity (15 direct human messages to Antigravity total, all prior to 2026-11-20 22:36 UTC+8); from Round 6 onward, forwarding switched to an automated system relay (26 relays, zero direct human messages to Antigravity during Rounds 5–8). The single Codabench private test submission (Submission 1062221) was executed by the auditor on the human's behalf under human authorization. Every intervention, retracted claim, and scope downgrade across Rounds 1–8—including a shared Round 7 oversight omitting literature screening prior to locking a draft thesis—is archived unvarnished in the public audit ledger, prior-art verification log, reference check table, and autonomy evidence pack (Section 6). Neither AI system is listed as an author.

---

## 3. Results

### 3.1 Metric-Dependent Ranking Reversal Across Three Perturb-Seq Atlases

Consistent with Miller et al. (2026) and Mao et al. (2026), method rankings invert systematically across all three human CRISPRi atlases between centroid/differential-expression metrics ($\Delta r$, $\mathrm{de50\_}r$, $\mathrm{de50\_sign}$) and perturbation-specificity metrics ($v_5$, $\mathrm{disc}_{\mathrm{raw}}$) (**Table 1**, **Figure 1**).

| Dataset | Method | $\Delta r$ ($=\mathrm{cent\_tr}$) | $\mathrm{cent\_val}\ r$ | $\mathrm{de50}\ r$ | $\mathrm{de50\_sign}$ | $v_5$ (Orth.) | $\mathrm{disc}_{\mathrm{raw}}$ |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **xatlas** | Constant Mean | 0.0608 (−0.0019) | 0.0519 | 0.0941 | 0.5420 | 0.0090 | 0.5000* |
| (Dev 10) | AE (Linear) | **0.1367** (**0.1404**) | **0.1308** | 0.2039 | **0.6048** | 0.0672 | 0.5856* |
| | S0 (Ridge) | 0.1327 (0.1363) | 0.1258 | **0.2050** | 0.6046 | 0.0613 | 0.5891* |
| | Simplified GNN | 0.0432 | — | 0.1025 | 0.5464 | 0.0543 | 0.5900* |
| | **S5 (Primary)** | 0.0813 (0.0871) | 0.1098 | 0.1576 | 0.5768 | 0.0887 | 0.6262* |
| | S6 (+Reciprocal) | 0.0937 (0.0990) | 0.1185 | 0.1885 | 0.5875 | 0.1027 | **0.6453*** |
| | **S7 (Full)** | 0.0995 (0.1051) | 0.1257 | 0.1789 | 0.5857 | **0.1099** | 0.6286* |
| **Replogle 2022** | Constant Mean | −0.0322 (−0.0028) | 0.0520 | −0.0531 | 0.4992 | 0.0009 | 0.5000 |
| (Strict 10) | AE (Linear) | 0.1470 | 0.1581 | 0.2026 | 0.6255 | 0.0967 | 0.6468 |
| | S0 (Ridge) | 0.1408 | 0.1517 | 0.2119 | 0.6256 | 0.0883 | 0.6511 |
| | Simplified GNN | 0.0971 | — | 0.1892 | 0.5999 | 0.1233 | 0.6403* |
| | GEARS (official package, pseudobulk-adapted) | 0.0486 | — | 0.1077 | 0.6009 | 0.0304 | 0.5621 |
| | **scGPT FrozenEmb Ridge** | 0.1802 | — | 0.2646 | **0.6536** | 0.1127 | 0.6859 |
| | **S5 (Primary)** | 0.1193 | 0.1382 | 0.2170 | 0.6064 | **0.1401** | 0.6884 |
| | OrthoFusion (Exploratory) | **0.1906** | — | **0.2829** | **0.6628** | **0.1401** | **0.6988** |
| | S6 (+Reciprocal) | 0.1106 | 0.1273 | 0.2155 | 0.6020 | 0.1292 | 0.6993* |
| | **S7 (Full)** | 0.1162 | 0.1346 | 0.2013 | 0.5962 | 0.1359 | 0.6736 |
| | S7 no S6 (Exploratory) | 0.1288 | 0.1496 | — | — | 0.1512 | — |
| **Nadig 2025** | Constant Mean | −0.0001 | 0.0292 | 0.0057 | 0.5000 | 0.0000 | 0.5000 |
| (Pre-Reg Unseen) | AE (Linear) | 0.1698 | 0.1787 | 0.2492 | 0.6457 | 0.1146 | 0.6665 |
| | S0 (Ridge) | 0.1610 | 0.1698 | 0.2511 | 0.6451 | 0.1042 | 0.6709 |
| | Simplified GNN | 0.1208 | — | 0.1878 | 0.6090 | 0.1432 | 0.6579 |
| | GEARS (official package, pseudobulk-adapted) | 0.1077 | — | 0.1973 | 0.6462 | 0.0577 | 0.6095 |
| | **scGPT FrozenEmb Ridge** | 0.2129 | — | 0.3064 | 0.6792 | 0.1512 | 0.7109 |
| | **S5 (Primary)** | 0.1502 | 0.1596 | 0.2265 | 0.6293 | 0.1741 | 0.7152 |
| | OrthoFusion (Exploratory) | **0.2144** | — | **0.3166** | **0.6808** | 0.1741 | **0.7324** |
| | S6 (+Reciprocal) | 0.1517 | 0.1612 | 0.2384 | 0.6307 | 0.1756 | 0.7184 |
| | **S7 (Full)** | 0.1617 | 0.1721 | 0.2242 | 0.6277 | 0.1871 | 0.6992 |
| | S7 no S6 (Exploratory) | 0.1644 | 0.1745 | 0.2164 | 0.6298 | **0.1901** | 0.6983 |

: Multi-metric comparison across three independent human CRISPRi Perturb-seq atlases (10-split macro averages). On xatlas, values in parentheses denote strict training-centroid subtraction ($\mathrm{centroid\_train\_}r$). Values marked with an asterisk (*) denote target-gene-masked discrimination percentiles, agreeing with unmasked $\mathrm{disc}_{\mathrm{raw}}$ within 0.002. OrthoFusion $\mathrm{de50}\ r$ and $\mathrm{de50\_sign}$ reflect the unified 7-metric DRF evaluation protocol (Section 6).

![Metric-dependent ranking reversal across three independent human Perturb-seq atlases (xatlas-orion-perturb, Replogle et al. 2022, and pre-registered unseen Nadig et al. 2025). Linear baselines (AE, S0) and frozen scGPT embedding probes lead on centroid-centered reconstruction ($\Delta r$), whereas frozen cross-line kernel S5 and S7 lead on orthogonal residual correlation ($v_5$) and raw perturbation discrimination ($\mathrm{disc}_{\mathrm{raw}}$).](../figures/fig1_metric_inversion_three_datasets.png){width=96%}

### 3.2 Mechanistic PC1 Axis Decomposition Beyond Centroid Subtraction

Direct audit verified that baseline vector $\mathbf{b}_c$ equals the reproducible training perturbation centroid across all three datasets—including xatlas ($\|\mathrm{context\_baseline} - \mu_{\mathrm{rep}}\|_2 = 0.0000$, $r = 1.000000$)—establishing exact identity $\Delta r \equiv \mathrm{centroid\_train\_}r$. Even after subtracting the training or validation centroid, linear models (AE, S0) still outperform S5 on $\mathrm{centroid\_val\_}r$ (S5 minus AE: −0.0210 on xatlas, −0.0199 on Replogle, −0.0191 on Nadig).

Testing our pre-registered mechanistic hypothesis that this residual linear advantage arises from capturing scalar amplitude along the single dominant training residual axis $\mathbf{a}_c$ (PC1; pass rule $\mathrm{off\_axis\_}r(\mathrm{AE} - \mathrm{S5}) \le \frac{1}{3}\,\Delta r(\mathrm{AE} - \mathrm{S5})$) confirmed the mechanism across all 3/3 atlases (**Table 2**):

| Atlas | Comparison | Full $\Delta r$ Gap | Off-Axis $r$ Gap | Retention Ratio | Axis $r$ (Ref vs. S5) | Verdict |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **xatlas** | AE − S5 | +0.0533 | +0.0174 | **0.327** ($\le 1/3$) | 0.1870 vs. −0.0974 | **Confirmed** |
| **Replogle 2022** | AE − S5 | +0.0277 | −0.0017 | **−0.063** (sign flip) | 0.3379 vs. 0.0770 | **Confirmed** |
| **Replogle 2022** | scGPT Probe − S5 | +0.0609 | +0.0240 | 0.395 | 0.4270 vs. 0.0770 | Supporting |
| **Nadig 2025** | AE − S5 | +0.0195 | −0.0337 | **−1.725** (sign flip) | 0.4206 vs. −0.0320 | **Confirmed** |
| **Nadig 2025** | scGPT Probe − S5 | +0.0627 | +0.0071 | 0.113 | 0.4774 vs. −0.0320 | Supporting |

: Pre-registered PC1 shared-axis vs. off-axis decomposition across three human CRISPRi atlases (10-split means). Projecting out the single training PC1 axis cuts AE's $\Delta r$ lead over S5 by 67%–273% (flipping its sign on Replogle and Nadig) and cuts the scGPT probe's lead by 61% on Replogle and 89% on Nadig.

### 3.3 Miller DRF Metric Calibration and Rankings Across Calibration Tiers

Applying the Miller et al. (2026) positive/negative-control DRF calibration framework to Replogle 2022 (K562 independent Day 6 essential replicate vs. Day 8 genome-wide screen) and Nadig 2025 (HepG2 + Jurkat 50/50 single-cell split-half controls) yields highly consistent calibration tiers across both laboratories (**Table 3**, ordered from lowest to highest $\mathrm{DRF}_{\mathrm{raw}}$):

| Metric | Nadig / Replogle $\mathrm{DRF}_{\mathrm{raw}}$ | Nadig / Replogle $\mathrm{DRF}_{\mathrm{interp}}$ | Calibration Tier | Replogle Ranking (Best → Worst among Core Models) | Nadig Ranking (Best → Worst among Core Models) |
| :--- | :---: | :---: | :--- | :--- | :--- |
| **MSE** *(lower better)* | **0.035 / 0.063** | 0.223 / 0.202 | **Poor** | OrthoF (0.0629) < scGPT (0.0642) < S5 (0.0660) < AE (0.0669) < Mean (0.0693) < **S7 (0.0859)** | scGPT (0.0641) < OrthoF (0.0652) < S5 (0.0690) < AE (0.0699) < Mean (0.0707) < **S7 (0.0813)** |
| **$v_5$** | **0.377 / 0.420** | 0.379 / 0.400 | **Moderate** | S5 = OrthoF (0.1401) > S7 (0.1359) > scGPT (0.1127) > AE (0.0967) | S7 (0.1871) > S5 = OrthoF (0.1741) > scGPT (0.1512) > AE (0.1146) |
| **$\Delta r$** | 0.484 / 0.513 | 0.463 / 0.475 | Moderate | OrthoF (0.1906) > scGPT (0.1802) > AE (0.1470) > S5 (0.1193) > S7 (0.1162) | OrthoF (0.2144) > scGPT (0.2129) > AE (0.1698) > S7 (0.1617) > S5 (0.1502) |
| **$\mathrm{de50\_}r$** | **0.808 / 0.843** | 0.756 / 0.703 | **High** | OrthoF (0.2829) > scGPT (0.2646) > **S5 (0.2170)** > AE (0.2026) > S7 (0.2013) | OrthoF (0.3166) > scGPT (0.3064) > AE (0.2492) > **S5 (0.2265)** > S7 (0.2242) |
| **$\mathrm{de50\_sign}$** | **0.849 / 0.902** | 0.849 / 0.902 | **Very High** | OrthoF (0.6628) > scGPT (0.6536) > AE (0.6255) > **S5 (0.6064)** > S7 (0.5962) | OrthoF (0.6808) > scGPT (0.6792) > AE (0.6457) > **S5 (0.6293)** > S7 (0.6277) |
| **$\mathrm{disc}_{\mathrm{raw}}$** | **0.973 / 0.976** | 0.969 / 0.987 | **Highest** | OrthoF (0.6988) > S5 (0.6884) $\approx$ scGPT (0.6859) > S7 (0.6736) > AE (0.6468) | OrthoF (0.7324) > S5 (0.7152) $\approx$ scGPT (0.7109) > S7 (0.6992) > AE (0.6665) |
| **$\mathrm{disc}_{v_5}$** | **0.977 / 0.983** | 0.979 / 0.989 | **Highest** | S5 = OrthoF (0.6886) > scGPT (0.6807) > S7 (0.6744) > AE (0.6387) | S5 = OrthoF (0.7157) > scGPT (0.7100) > S7 (0.6993) > AE (0.6509) |

: Full Miller et al. (2026) Dynamic Range Fraction (DRF) calibration table and model rankings across Nadig 2025 (HepG2 + Jurkat split-half) and Replogle 2022 (K562 replicate), ordered by ascending $\mathrm{DRF}_{\mathrm{raw}}$ (10-split means). OrthoF = exploratory OrthoFusion; scGPT = scGPT FrozenEmb Ridge.

Four conclusions follow directly from **Table 3** under uniform $\mathrm{DRF}_{\mathrm{raw}}$ calibration—exposing critical limitations of standalone S5 and S7:

1. **Confirmation of Miller et al. (2026) calibration diagnosis**: Standard MSE exhibits near-zero raw calibration ($\mathrm{DRF}_{\mathrm{raw}}$ 0.035 on Nadig, 0.063 on Replogle), whereas top-50 differential-expression metrics ($\mathrm{de50\_}r$ 0.808–0.843, $\mathrm{de50\_sign}$ 0.849–0.902) and perturbation discrimination percentiles ($\mathrm{disc}_{\mathrm{raw}}$ 0.973–0.976, $\mathrm{disc}_{v_5}$ 0.977–0.983) are well calibrated.
2. **S5's advantage over the frozen foundation-model probe is confined to a moderately calibrated metric**: Orthogonal residual correlation $v_5$—the sole metric on which zero-pretraining S5 clearly outperforms the 33M-cell frozen scGPT probe—achieves only moderate calibration ($\mathrm{DRF}_{\mathrm{raw}}$ 0.377 on Nadig, 0.420 on Replogle), **lower than even $\Delta r$** ($\mathrm{DRF}_{\mathrm{raw}}$ 0.484 on Nadig, 0.513 on Replogle).
3. **Unfavorable performance of standalone S5 and S7 on the best-calibrated metrics and MSE**: Across the four best-calibrated metrics ($\mathrm{DRF}_{\mathrm{raw}}$ 0.808–0.983), standalone S5 statistically **ties** the scGPT probe on $\mathrm{disc}_{\mathrm{raw}}$ (0/10 splits significant), edges slightly ahead numerically on $\mathrm{disc}_{v_5}$, and **trails substantially on both $\mathrm{de50\_sign}$** (Replogle 0.6064 vs. 0.6536 scGPT and 0.6255 AE; Nadig 0.6293 vs. 0.6792 scGPT and 0.6457 AE) and **$\mathrm{de50\_}r$** (Replogle 0.2170 vs. 0.2646 scGPT; Nadig 0.2265 vs. 0.3064 scGPT and 0.2492 AE). On MSE, S5 trails the scGPT probe (Replogle 4/10 CI > 0; Nadig 10/10 CI > 0), and S7's manifold snapping degrades MSE worse than Constant Mean on both atlases (0.0859 vs. 0.0693; 0.0813 vs. 0.0707).
4. **Exploratory OrthoFusion reconciles the high-DRF metrics on seen atlases**: Without tuned weights, exploratory OrthoFusion ranks #1 or tied #1 across all four best-calibrated metrics simultaneously on both seen atlases and lowers MSE relative to S5 on 10/10 splits on both atlases, though unseen fourth-dataset confirmation remains pending.

### 3.4 Pre-Registered Unseen External Confirmation on Nadig et al. (2025)

Prior to reading any Nadig et al. (2025) expression value, we locked primary method S5, preprocessing rules, and decision criteria in git. All three hypotheses were confirmed (**Figure 2**; seeds 42, 0, and 1 independently replicated verbatim by the auditor):

- **H1 (Primary — $\mathrm{disc}_{\mathrm{raw}}$)**: S5 minus AE = +0.0487 (0.7152 vs. 0.6665), mean 95% CI [+0.0161, +0.0812], **9/10 splits strictly > 0** (confirmed; caveat: several split lower bounds lie near zero, +0.0011 to +0.0038, plus −0.0080 on seed 2024).
- **H2 (Primary — $\Delta r$ inversion)**: S5 minus AE = −0.0195 (0.1502 vs. 0.1698), **9/10 splits negative, 4/10 CI < 0**, 0/10 CI > 0 (confirmed; weaker per-split significance than Replogle's 8/10).
- **Joint primary thesis (H1 and H2)**: Confirmed.
- **H3 (Secondary — $v_5$)**: S5 minus AE = +0.0595 (0.1741 vs. 0.1146), mean 95% CI [+0.0437, +0.0753], **10/10 splits strictly > 0** (confirmed).

S5 also beat Simplified GNN ($v_5$ +0.0309 [9/10], $\mathrm{disc}_{\mathrm{raw}}$ +0.0572 [10/10], $\Delta r$ +0.0295 [9/10]) and Simplified ResMLP ($v_5$ +0.0449 [10/10], $\mathrm{disc}_{\mathrm{raw}}$ +0.0529 [10/10]). S6 gave the highest $\mathrm{disc}_{\mathrm{raw}}$ (0.7184), while S7 raised $v_5$ (0.1871, 10/10 vs. AE) at the expense of $\mathrm{disc}_{\mathrm{raw}}$ (0.6992, 4/10 vs. AE).

![Pre-registered zero-look confirmatory forest plot on unseen Nadig et al. (2025; HepG2 and Jurkat) showing split-level 4,000-resample stratified row-level paired bootstrap 95% confidence intervals for S5 minus AE across $\mathrm{disc}_{\mathrm{raw}}$ (H1, 9/10 CI > 0), $\Delta r$ (H2, 9/10 negative, 4/10 CI < 0), and $v_5$ (H3, 10/10 CI > 0) (A), and cross-comparator effect sizes including scGPT Probe and GEARS (official package, pseudobulk-adapted) (B).](../figures/fig2_preregistered_confirmation_forest.png){width=96%}

### 3.5 External Deep Learning, Frozen scGPT Probe, and Exploratory Orthogonal Fusion

1. **GEARS (official package, pseudobulk-adapted)**: Across 10 splits, official cell-gears 0.1.2 under our pseudobulk-adapted input protocol achieved $v_5 = 0.0304 \pm 0.0168$, $\Delta r = 0.0486 \pm 0.0253$, and $\mathrm{disc}_{\mathrm{raw}} = 0.5621 \pm 0.0208$ on Replogle 2022 (GO coverage 99.83% RPE1 / 99.38% K562), and $v_5 = 0.0577 \pm 0.0143$, $\Delta r = 0.1077 \pm 0.0174$, $\mathrm{disc}_{\mathrm{raw}} = 0.6095 \pm 0.0149$, and $\mathrm{de50\_sign} = 0.6462 \pm 0.0095$ on Nadig 2025 (GO coverage 99.66% HepG2 / 99.79% Jurkat). Under pseudobulk adaptation, GEARS trailed both AE and S5 on every metric; a native single-cell GEARS evaluation on Nadig 2025 is currently running (Section 5).
2. **Frozen scGPT gene-embedding probe**: Consistent with Miller et al. (2026), scGPT FrozenEmb Ridge is our strongest single baseline, beating linear AE across every metric on both Replogle ($\Delta r$ +0.0332 [8/10], $\mathrm{disc}_{\mathrm{raw}}$ +0.0392 [8/10], $\mathrm{disc}_{v_5}$ +0.0420 [10/10], $v_5$ +0.0160 [6/10]) and Nadig ($\Delta r$ +0.0432 [10/10], $\mathrm{disc}_{\mathrm{raw}}$ +0.0444 [8/10], $\mathrm{disc}_{v_5}$ +0.0591 [10/10], $v_5$ +0.0365 [10/10]). Directly comparing zero-pretraining S5 vs. scGPT FrozenEmb Ridge shows S5 wins **only on moderately calibrated $v_5$** (+0.0274 [10/10] Replogle; +0.0230 [7/10] Nadig; $\mathrm{DRF}_{\mathrm{raw}}$ 0.377–0.420), **statistically ties on $\mathrm{disc}_{\mathrm{raw}}$** (+0.0025 / +0.0043, 0/10 significant), and **trails substantially on $\Delta r$** (−0.0609 / −0.0627, 10/10 CI < 0) as well as best-calibrated $\mathrm{de50\_sign}$ and $\mathrm{de50\_}r$ (**Table 3**).
3. **Exploratory OrthoFusion and fourth-dataset hard stop**: On seen Replogle and Nadig atlases, zero-weight OrthoFusion passes exploratory non-inferiority (H4: 0/10 splits CI < 0 vs. scGPT probe on $\Delta r$, exact algebraic parity with S5 on $v_5$ within $10^{-12}$) while not meeting dual-sided synergy significance (H5: $\ge 7/10$ CI > 0 on $\mathrm{disc}_{\mathrm{raw}}$ vs. both parents; observed Replogle 5/10 vs. scGPT and 0/10 vs. S5; Nadig 7/10 vs. scGPT and 3/10 vs. S5). Screening independent non-Orion public multi-line CRISPRi studies under our locked statistical power gate ($\ge 100$ reproducible targets per line, $\ge 80$ shared reproducible targets, $\ge 2{,}000$ shared readout genes), Candidate 1 (Feng et al., 2026; Figshare 27989294 TargetedScreen, top-2 donor iPSC lines jejf\_2 and zapk\_3) yielded only 38 (27 in readout) and 86 (62 in readout) reproducible targets (37 shared) at Day 3 harvest, and Candidate 2 (Jiang et al., 2025; Zenodo 14518762) contained only 44–61 targets per matched cytokine context (IFNB 61, IFNG 59, TNFA 55, TGFB1 52, INS 44). We enforced the pre-registered hard stop before running any model and record fourth-dataset confirmation as **pending**.

### 3.6 Leave-One-Module-Out Transferability and ±50% Sensitivity Analysis

Dual-dataset LOMO ablation (**Figure 3A**) separates transferable modules (positive on both xatlas and Replogle: S1 HybridTOM +0.0109 [7/10] / +0.0103 [10/10]; S7 Wiener snapping +0.0072 / +0.0067 on $v_5$, hurting $\mathrm{disc}_{\mathrm{raw}}$ and MSE; S5 Cross-Line Bridge +0.0020 / +0.0061 [9/10]; S2 spectral deflation +0.0108 / +0.0035 [1/10]) from dataset-specific modules (S6 Reciprocal Borrowing +0.0130 [8/10] on xatlas vs. **−0.0153 [9/10 CI < 0]** on Replogle; S4 non-reproducible anchors +0.0043 vs. −0.0006).

One-factor-at-a-time ±50% perturbation of all 10 frozen S5 constants (**Figure 3B**) retains $v_5$ 10/10 splits CI > 0 vs. AE on 18/20 perturbed configurations on Replogle and 20/20 on Nadig, while exposing a sharp **concave-power boundary ($p < 1$)**: setting translator kernel power $p_k = 0.75$ (−50%) collapses Replogle ($v_5 = 0.0602$, shift −0.0798; $\mathrm{disc}_{\mathrm{raw}} = 0.6052$, 0/10 CI > 0) and causes the largest Nadig drop ($v_5 = 0.1571$, −0.0170). Spectral diagnostics on Replogle seed 42 show concave power $\mathrm{sign}(G)|G|^{0.75}$ destroys positive semidefiniteness ($\lambda_{\min}(K) = -6.177$), introducing six negative eigenvalues in $K + 2I$ ($\lambda_{\min} = -4.177$) and inflating the condition number from 17.2 ($p_k = 1.5$) to 4,473.5 ($p_k = 0.75$), whereas every convex power $p_k \ge 1.0$ preserves strict positive definiteness.

![Leave-one-module-out (LOMO) transferability ablation across xatlas and Replogle 2022 (A) and one-factor-at-a-time ±50% sensitivity sweep of all 10 frozen S5 constants on Replogle 2022 and Nadig 2025 (B), pinpointing the concave-power ($p_k < 1$) non-PSD instability boundary.](../figures/fig3_lomo_ablation_and_s5_sensitivity.png){width=96%}

### 3.7 Biological Pathway Enrichment, *Cis*-Proximity Audit, and Official Codabench Test

Evaluating local Reactome enrichment over the GMT-annotated evaluated background ($M_{\mathrm{bg}} = 1{,}452$, **Figure 4A**): Top-150 S7 minus S0 gain genes (128 annotated) enrich **mRNA Splicing – Major Pathway** (38/145, 2.97-fold, BH FDR = $2.90 \times 10^{-8}$; 0/50 random controls), and Top-150 S5 minus S2 cross-line bridge gain genes (131 annotated) enrich **Mitochondrial Translation** (20/82, 2.70-fold, BH FDR = $6.98 \times 10^{-3}$; 0/50 random controls). Full-universe sensitivity ($M_{\mathrm{bg}} = 1{,}784$) is consistent (3.12-fold, FDR = $1.18 \times 10^{-8}$; 2.90-fold, FDR = $2.76 \times 10^{-3}$). Auditing S6 reciprocal donors against GRCh38 GENCODE v46 coordinates (10,530/10,539 Replogle genes mapped, **Figure 4B**) shows negative-$z$ pairs ($z < -4$) are enriched 474-fold for $< 10\ \mathrm{kb}$ proximity (1.422% vs. 0.003% null), yet $< 10\ \mathrm{kb}$ pairs contribute only 0.255% of active $|z| > 4$ weight and removing them retains 97.41% of Top-40 own-reciprocal gains on Replogle; on anonymized xatlas (78% negative-$z$ top donors), local CRISPRi *cis*-repression cannot be ruled out.

On the blind official BioXArena-Pro xatlas-orion-perturb private test set, a single frozen submission (Submission 1062221, panperturb-bridge human-ai-collab, status ok) scored **0.0998** (+0.0173 above the best 2-hour autonomous agent at 0.0825, within the 0.02 noise margin; validation-to-test shift $0.0998 - 0.1099 = -0.0101$). Because this pipeline reflects multi-round iteration by a dual-AI executor–auditor system under a human-set goal (~25 fixed constants), it is not directly comparable to single-shot 2-hour autonomous agent runs.

![Reactome pathway enrichment of S7 minus S0 (mRNA Splicing – Major Pathway) and S5 minus S2 (Mitochondrial Translation) gain genes over the GMT-annotated evaluated background ($M_{\mathrm{bg}} = 1{,}452$) (A), and GRCh38 GENCODE v46 genomic-distance audit of S6 reciprocal donors on Replogle 2022 (B).](../figures/fig4_biology_reactome_and_cis_audit.png){width=96%}

---

## 4. Discussion

Rather than claiming metric-dependent ranking reversal or seed-based cross-cell-line transfer as new discoveries—both of which are established in prior literature (Miller et al., 2026; Mao et al., 2026; Dhaka, Gao, and Ideker, 2026; Adduri et al., 2025; Zhu and Li, 2026)—this study establishes four bounded empirical conclusions. First, under commit-hashed zero-look pre-registration and stratified row-level bootstrap CIs, ranking reversal reproduces on unseen Nadig et al. (2025) data. Second, our PC1 decomposition confirms across 3/3 atlases that even after Systema training-centroid centering, linear models' advantage over cross-line residual kernels resides predominantly along the single shared training PC1 variation axis (removing 67%–273% of AE's $\Delta r$ lead). Third, our **uniform $\mathrm{DRF}_{\mathrm{raw}}$ calibration audit** directly tempers any strong claim for standalone S5: while MSE is poorly calibrated ($\mathrm{DRF}_{\mathrm{raw}}$ 0.035–0.063) and top-50 DE and discrimination metrics are well calibrated ($\mathrm{DRF}_{\mathrm{raw}}$ 0.808–0.983), **$v_5$ is only moderately calibrated ($\mathrm{DRF}_{\mathrm{raw}}$ 0.377–0.420, below $\Delta r$ at 0.484–0.513)**. Because $v_5$ is the sole metric on which S5 clearly beats the 33M-cell frozen scGPT probe, **S5's advantage over the frozen foundation-model probe is confined to a moderately calibrated metric**; on the best-calibrated metrics ($\mathrm{DRF}_{\mathrm{raw}}$ 0.808–0.983), standalone S5 ties ($\mathrm{disc}_{\mathrm{raw}}$) or trails ($\mathrm{de50\_sign}$, $\mathrm{de50\_}r$) the scGPT probe. Fourth, positioning **S5 directly against its closest prior work, HyperMap (Dhaka, Gao, and Ideker, 2026)**: HyperMap uses sparse perturbation "seeds" within a parameter-efficient neural meta-learner trained across an 18-context compendium (HyperMapDB) to generalize across cell lines, iPSC donors, and small molecules, whereas S5 operates strictly in closed form with zero neural training parameters on a single paired atlas ($A \leftrightarrow B$) via a source-centered signed-power cosine KRR translator whose constants transfer zero-shot across laboratories. S5 thus serves as an interpretable, convex ($p_k \ge 1$) non-neural reference baseline rather than a foundation-model replacement, while exploratory OrthoFusion illustrates how combining foundation-model shared-subspace accuracy with S5's orthogonal residuals reconciles all four best-calibrated metrics numerically pending fourth-dataset confirmation.

From a metascience standpoint, our eight-round git trace provides a concrete empirical case study alongside CrossAudit (Dong and Chen, 2026), Virtual Lab (Swanson et al., 2025), A Society of AI Scientists (Ahmad et al., 2026), and conference reproducibility auditing (Dong, Zhang, and Zheng, 2026). Independent code re-execution by a cross-vendor auditor caught 34 distinct issues (retracting 13 inflated early claims and enforcing negative-gate hard stops on two underpowered datasets), while also exposing a shared auditor–executor failure mode—omitting systematic prior-art screening before Round 7—that motivated our mandatory prior-art verification and API-checked reference gates.

---

## 5. Limitations and In-Progress Audits

1. **Prior art on core themes**: Metric calibration and ranking reversal were established by Miller et al. (2026) and Mao et al. (2026); seed-based cross-context perturbation transfer exists in HyperMap (Dhaka et al., 2026), STATE (Adduri et al., 2025), TranScouter (Zhu and Li, 2026), and GRAFT (Luo et al., 2026); and git-native cross-vendor AI auditing is formalized in CrossAudit (Dong and Chen, 2026). No priority claim is made for these paradigms.
2. **Narrow calibration boundary and unfavorable performance of standalone S5 and S7**: Under uniform $\mathrm{DRF}_{\mathrm{raw}}$ calibration, S5's sole clear advantage over scGPT FrozenEmb Ridge ($v_5$) lies on a moderately calibrated metric ($\mathrm{DRF}_{\mathrm{raw}}$ 0.377–0.420, lower than $\Delta r$ at 0.484–0.513). Across the best-calibrated metrics ($\mathrm{DRF}_{\mathrm{raw}}$ 0.808–0.983), standalone S5 ties scGPT on $\mathrm{disc}_{\mathrm{raw}}$ and trails both scGPT and AE on $\mathrm{de50\_sign}$ (and scGPT on $\mathrm{de50\_}r$). Stage S7 degrades MSE worse than Constant Mean on both Replogle and Nadig.
3. **Unconfirmed fourth-dataset validation for OrthoFusion**: OrthoFusion results on Replogle 2022 and Nadig 2025 are exploratory (seen datasets) and do not meet dual-sided H5 significance; fourth-dataset confirmation remains pending after Feng et al. (2026) and Jiang et al. (2025) failed locked power gates.
4. **In-progress native single-cell GEARS audit**: Completed GEARS evaluations used pseudobulk adaptation with differential-expression calculation skipped (GEARS, official package, pseudobulk-adapted). A native single-cell GEARS evaluation on Nadig 2025 (native differential expression enabled, seeds 42, 0, and 1) is currently running and excluded from Preprint v2.1 claims.
5. **Compute and architecture constraints for full foundation-model fine-tuning**: Full 10-split scGPT TransformerGenerator fine-tuning (~100 hours on Apple Silicon MPS without CUDA FlashAttention; an exploratory 2-seed Nadig run was discontinued after 2h 41m) and Arc Institute STATE (defaulting to one-hot perturbation representations without non-commercial SE-600M embeddings and multi-GPU CUDA) were not completed.
6. **Potential scGPT pretraining corpus overlap**: scGPT human was pretrained on 33M normal human cells from CELLxGENE (2023); while Nadig et al. (2025; GEO release May 2024) post-dates scGPT pretraining, Replogle et al. (2022) pre-dates 2023 and broad co-expression overlap cannot be ruled out.
7. **Other empirical caveats**: The Codabench score (0.0998) sits within 0.02 of the best 2-hour autonomous agent (0.0825) and belongs to a different development category; several Nadig H1 split lower bounds lie near zero; S6 cannot rule out CRISPRi *cis*-repression on anonymized xatlas; and concave kernel power ($p_k < 1$) violates positive semidefiniteness.

---

## 6. Data, Code, and Audit Availability

All code, frozen manifests, pre-registration protocols, prior-art logs, reference API verification tables, and audit artifacts are publicly archived at `https://github.com/ArtificialZeng/AI-Biology-Lab` under `projects/bioxarena-pro-scaffold/`, with zero hardcoded local file paths (`BIOXARENA_DIR`, `REPLOGLE_DATA_DIR`, `NADIG_DATA_DIR`, `FENG_DATA_DIR`, `SCGPT_WEIGHTS_DIR`) and a locked dependency specification (`requirements-lock.txt`):

- **Frozen Solver & Codabench Pack**: `src/xatlas_panhuman_solver.py` (commit `3609811`), `outputs/perturbation_baselines/xatlas-orion-perturb/CODABENCH_FREEZE_MANIFEST.json` (`submission.csv` SHA-256 `e2613a9b6616368ba5ca29edf132fe023761f0f72072b25176724207daa6a207`), and the official Codabench Submission 1062221 archive.
- **Governance, Pre-Registration & Audit Artifacts**: `reports/ROUND7_PREREGISTRATION.md` (commits `c5294c0`, `87f4bc2`, `b9fe2bb`), `reports/ROUND8_PREREGISTRATION.md` (commits `ce82040`, `840b502`), `reports/PRIOR_ART_LOG.md`, `reports/REFERENCE_CHECK.csv`, `reports/AUDIT_LEDGER.csv`, `reports/AUDIT_LEDGER_SUMMARY.md`, `reports/AUTONOMY_EVIDENCE.md`, and `reports/autonomy_timeline.csv`.
- **Result Summary Files**: Multi-metric summaries (`reports/t2_all_5metric_summary.json`, `reports/t2_discrim_delta_row_ci_summary.json`, `reports/r74_systema_*_summary.json`), Nadig pre-registered confirmation (`reports/r71_nadig_confirmation_summary.json`), PC1 axis decomposition (`reports/r82_axis_decomposition_*_summary.json`), Miller DRF calibration (`reports/r86_miller_drf_{replogle,nadig}_summary.json`, `reports/r86_model_7metrics_*.csv`, `reports/r86_mse_row_ci_*.csv`), external baselines (`reports/r72_official_gears_*_summary.json`, `reports/r72_scgpt_probe_*`), exploratory fusion & fourth-dataset gate (`reports/r83_exploratory_orthofusion_*_summary.json`, `reports/r83_confirmatory_feng2025_summary.json`), LOMO & sensitivity (`reports/t4_lomo_summary.json`, `reports/r73_sensitivity_*_summary.json`), and biology/cis audits (`reports/t5_local_reactome_enrichment_summary.json`, `reports/t6_replogle_crispri_proximity_summary.json`).
- **Public External Data & Checkpoints**: Replogle et al. (2022) Figshare+ DOI `10.25452/figshare.plus.20029387` (`rpe1` SHA-256 `a3c5bfd0...`, `K562_gwps` `37e48c47...`, `K562_essential` `c1ca6456...`); Nadig et al. (2025) NCBI GEO `GSE264667` (`hepg2` `e1ad7c3c...`, `jurkat` `ffbe15f2...`); Feng et al. (2026, halted at QC gate) Figshare DOI `10.6084/m9.figshare.27989294.v2`; Jiang et al. (2025, halted at target-count audit) Zenodo `14518762`; `scGPT_human` weights (MIT License, `ad0252a1...`) and vocab (`ee2b2c90...`); Reactome Pathways GMT (CC BY 4.0, `a9feaa3a...`); GENCODE v46 Basic GTF (`d620d548...`).

---

## References

1. Adduri, A. K., Gautam, D., Bevilacqua, B., Imran, A. et al. (2025). Predicting cellular responses to perturbation across diverse contexts with STATE. *bioRxiv*, doi:10.1101/2025.06.26.661135.
2. Ahlmann-Eltze, C., Huber, W. & Anders, S. (2025). Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines. *Nature Methods*, doi:10.1038/s41592-025-02772-6.
3. Ahmad, R., Anderson, M., Qu, A. et al. (2026). A Society of AI Scientists. *alphaXiv*, `https://www.alphaxiv.org/abs/2611.society-ai-scientists`.
4. Cui, H., Wang, C., Maan, H., Pang, K. et al. (2024). scGPT: toward building a foundation model for single-cell multi-omics using generative AI. *Nature Methods*, doi:10.1038/s41592-024-02201-0.
5. Dhaka, B., Gao, J. & Ideker, T. (2026). HyperMap: An Efficient Framework for Transferring Perturbation Responses Across Diverse Biological Contexts. *bioRxiv*, doi:10.64898/2026.04.23.720505.
6. Ding, T., Nannapaneni, A., Liu, B. & Zhang, L. (2026). Autonomous Research Agents: A Survey of AI Scientists and the Verification Gap. *arXiv:2608.05179*.
7. Dong, Y., Zhang, L. & Zheng, T. (2026). Author-in-the-loop Agentic Reproducibility Auditing for a Machine Learning Conference. *arXiv:2610.22465*.
8. Dong, Z. & Chen, Y. (2026). CrossAudit: A Git-Native, Cross-Vendor Audit Loop for Agentic Science. *arXiv:2608.28631*.
9. Feng, C., Peets, E. M., Zhou, Y., Crepaldi, L. et al. (2026). A genome-scale single-cell CRISPRi map of trans gene regulation across human pluripotent stem cell lines. *Cell Genomics*, doi:10.1016/j.xgen.2025.101076 (PMID 41330380).
10. Jiang, L., Dalgarno, C., Papalexi, E., Mascio, I. et al. (2025). Systematic reconstruction of molecular pathway signatures using scalable single-cell perturbation screens. *Nature Cell Biology*, doi:10.1038/s41556-025-01622-z (PMID 40011560).
11. Luo, C., Hu, M., Wu, Y., Guo, D. et al. (2026). Cross-modal transfer enables zero-shot prediction of drug responses across cell types. *bioRxiv*, doi:10.64898/2026.10.27.762107.
12. Mao, X., Zhang, S., Wen, Q., Wen, X. et al. (2026). Benchmarking virtual cell models for in-the-wild perturbation response. *arXiv:2604.27646*.
13. Miller, H. E., Mejia, G. M., Leblanc, F. J. A., Swain, B., Wang, B. & de Lima Camillo, L. P. (2026). Deep learning perturbation models can outperform baselines on calibrated metrics. *Nature Biotechnology*, doi:10.1038/s41587-026-03307-w; *bioRxiv*, doi:10.1101/2025.10.20.683304.
14. Nadig, A., Replogle, J. M., Pogson, A. N., Murthy, M. et al. (2025). Transcriptome-wide analysis of differential expression in perturbation atlases. *Nature Genetics*, doi:10.1038/s41588-025-02169-3 (PMID 40259084; GEO GSE264667).
15. Propose, Critique, Falsify: Benchmarking Self-Verifying AI Scientists. *ICML 2026* (`https://icml.cc/virtual/2026/73543`).
16. Replogle, J. M., Saunders, R. A., Pogson, A. N., Hussmann, J. A. et al. (2022). Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq. *Cell*, doi:10.1016/j.cell.2022.05.013.
17. Roohani, Y., Huang, K. & Leskovec, J. (2024). Predicting transcriptional outcomes of novel multigene perturbations with GEARS. *Nature Biotechnology*, doi:10.1038/s41587-023-01905-6.
18. Swanson, K., Wu, W., Bulaong, N. L., Pak, J. E. & Zou, J. (2025). The Virtual Lab of AI agents designs new SARS-CoV-2 nanobodies. *Nature*, doi:10.1038/s41586-025-09442-9.
19. Vaccaro, M. (2026). Preregistration for Experiments with AI Agents. *arXiv:2606.11217*.
20. Viñas Torné, R., Wiatrak, M., Piran, Z., Fan, S. et al. (2025). Systema: a framework for evaluating genetic perturbation response prediction beyond systematic variation. *Nature Biotechnology*, doi:10.1038/s41587-025-02777-8.
21. Zhu, O. & Li, J. (2026). A structured study of cross-condition prediction of transcriptional responses to gene perturbations. *bioRxiv*, doi:10.64898/2026.07.30.741892.
