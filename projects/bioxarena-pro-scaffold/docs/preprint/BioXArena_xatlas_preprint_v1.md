---
title: "Pre-Registered Three-Atlas Confirmation, Principal-Axis Decomposition, and Cross-Line Kernel Modeling Under Metric-Dependent Perturbation Rankings"
author:
  - "Zijian Zeng"
institute: "AI Biology Lab (github.com/ArtificialZeng)"
date: "November 21, 2026 (Preprint v1)"
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

Recent benchmarks reported that simple linear regression matches or outperforms complex deep learning models in predicting single-gene CRISPR perturbation transcriptomes (Ahlmann-Eltze et al., 2025; Viñas Torné et al., 2025), while subsequent calibration studies established that uncalibrated correlation and mean squared error (MSE) metrics poorly separate perturbation-specific signal from mean baselines, that changing evaluation metrics substantially alters model rankings, and that linear probes on frozen foundation-model embeddings perform strongly (Miller et al., 2026; arXiv:2604.27646). Building directly on those established findings, this study provides an audited, pre-registered empirical evaluation across three independent human Perturb-seq atlases (BioXArena-Pro xatlas-orion-perturb, Replogle et al. 2022, and unseen Nadig et al. 2025) with three completed contributions: (1) a pre-registered, zero-look confirmation of metric-dependent ranking reversal across laboratories under strict simultaneous union holdout and cell-line-stratified row-level paired bootstrap inference; (2) a direct mechanistic decomposition demonstrating that even after Systema training-centroid subtraction ($\Delta r \equiv \mathrm{centroid\_train\_}r$), linear models' residual advantage over cross-line residual kernels resides primarily along the single dominant shared training principal component (PC1)—removing that axis reduces the gap between the linear baseline (AE) and our cross-line kernel model (S5) to at most one-third or reverses its sign across all three atlases (xatlas: +0.0533 to +0.0174, ratio 0.327; Replogle: +0.0277 to −0.0017, ratio −0.063; Nadig: +0.0195 to −0.0337, ratio −1.725); and (3) a frozen, zero-pretraining cross-line kernel baseline (S5, Cross-Line Perturb-Bridge) alongside an exploratory parameter-free orthogonal decomposition (OrthoFusion) combining frozen scGPT embeddings on the shared subspace with S5 on the orthogonal complement.

After freezing all solver constants on xatlas-orion-perturb—where a single frozen Codabench private test submission scored 0.0998 (+0.0173 above the best 2-hour autonomous agent score of 0.0825, within the 0.02 noise margin and explicitly categorized as multi-round human-assisted engineering rather than single-shot autonomous execution)—we evaluated frozen S5 zero-shot across two external Perturb-seq atlases. In a pre-registered, zero-look confirmatory test locked in git prior to expression inspection on Nadig et al. (2025; HepG2 and Jurkat), all three pre-declared hypotheses were confirmed: primary raw-space perturbation discrimination H1 ($\mathrm{disc}_{\mathrm{raw}}$ S5 minus AE = +0.0487, 9/10 splits with 95% CI strictly above zero, though with several lower bounds near zero), primary metric-inversion hypothesis H2 ($\Delta r$ S5 minus AE = −0.0195, 9/10 splits negative and 4/10 splits with 95% CI strictly below zero), and secondary orthogonal-residual hypothesis H3 ($v_5$ S5 minus AE = +0.0595, 10/10 splits with 95% CI strictly above zero). Consistent with Miller et al. (2026), a linear ridge probe on frozen 33M-cell pretrained scGPT gene embeddings (scGPT FrozenEmb Ridge) is our strongest single baseline: zero-pretraining S5 outperforms it only on orthogonal residuals $v_5$ (+0.0274, 10/10 CI > 0 on Replogle; +0.0230, 7/10 CI > 0 on Nadig), statistically ties on $\mathrm{disc}_{\mathrm{raw}}$ (+0.0025 and +0.0043, 0/10 significant), and trails substantially on $\Delta r$ (−0.0609 and −0.0627, 10/10 CI < 0). Parameter-free orthogonal fusion (OrthoFusion) achieves exploratory numerical Pareto improvements on both seen atlases (Replogle $\Delta r = 0.1906$, $v_5 = 0.1401$, $\mathrm{disc}_{\mathrm{raw}} = 0.6988$; Nadig $\Delta r = 0.2144$, $v_5 = 0.1741$, $\mathrm{disc}_{\mathrm{raw}} = 0.7324$), satisfying exploratory non-inferiority (H4) without meeting dual-sided synergy significance (H5), while prospective fourth-dataset confirmation remains pending after two screened public multi-line CRISPRi studies (Feng et al., 2025; Jiang et al., 2025) failed pre-registered statistical power gates. Pathway enrichment over the GMT-annotated evaluated background ($M_{\mathrm{bg}} = 1{,}452$) shows S7 minus S0 and S5 minus S2 gains enrich mRNA Splicing – Major Pathway (2.97-fold, BH FDR = $2.90 \times 10^{-8}$) and Mitochondrial Translation (2.70-fold, BH FDR = $6.98 \times 10^{-3}$; 0/50 random controls). The entire eight-round project was conducted under a cross-vendor AI execution-and-auditing workflow with commit-hashed pre-registration, a formal prior-art verification log, and an unvarnished 34-entry public audit ledger documenting every corrected error and retracted claim.

---

## 1. Introduction

Predicting transcriptome-wide responses to unseen single-gene CRISPR perturbations is a central testbed for virtual-cell modeling. Early systematic evaluations (Ahlmann-Eltze et al., 2025; Viñas Torné et al., 2025) demonstrated that complex deep learning models frequently failed to outperform simple linear regression or perturbation-mean baselines under standard Pearson correlation metrics. In parallel, cross-context transfer frameworks—including large multi-context transformers such as STATE (Adduri et al., 2025) as well as meta-learning and cross-condition architectures such as HyperMap (2026), TranScouter (2026), and GRAFT (2026)—have actively explored transferring perturbation signatures across distinct cellular backgrounds.

Two pivotal benchmark studies subsequently showed that evaluation metric choice fundamentally governs model rankings in perturbation prediction. Miller et al. (2026; bioRxiv 10.1101/2025.10.20.683304) introduced a positive-control (interpolated duplicate) and negative-control (mean) calibration framework quantified by the Dynamic Range Fraction (DRF), proving across 14 datasets that standard $\Delta$-Pearson and MSE metrics are poorly calibrated, that switching to calibrated metrics alters method rankings, and that linear probes on frozen pretrained embeddings are highly competitive. Independently, arXiv:2604.27646 documented across virtual-cell benchmarks that different evaluation metrics yield substantially different model rankings. Accordingly, we do not claim metric-dependent ranking reversal or cross-cell-line transfer as new discoveries. Rather, building directly on Miller et al. (2026), Viñas Torné et al. (2025), Ahlmann-Eltze et al. (2025), and arXiv:2604.27646, this work contributes three completed empirical analyses and a transparent cross-vendor audit case study:

1. **Pre-registered, zero-look three-atlas confirmation**: Under strict simultaneous union holdout and 4,000-resample cell-line-stratified row-level paired bootstrap confidence intervals, we lock hypotheses, preprocessing, and frozen solver constants in git prior to inspecting unseen expression data on Nadig et al. (2025), confirming consistent metric-dependent ranking reversal across three independent laboratories.
2. **Direct rank-1 PC1 off-axis decomposition beyond centroid subtraction**: Systema (Viñas Torné et al., 2025) established that subtracting the mean perturbation centroid removes the constant shared shift. Here we verify algebraically that centroid subtraction alone ($\Delta r \equiv \mathrm{centroid\_train\_}r$) still leaves linear models (AE, S0) ahead of cross-line residual kernels, and directly test via pre-registered orthogonal projection that this residual linear advantage originates predominantly from predicting scalar amplitude along the single dominant shared training residual axis (PC1).
3. **Frozen zero-pretraining cross-line kernel modeling (S5) and exploratory parameter-free orthogonal fusion**: We evaluate a closed-form, source-centered dual-ridge kernel baseline (S5) whose constants remain frozen across all three human CRISPRi atlases, sharply localize its domain of advantage over a 33M-cell scGPT gene-embedding probe to the orthogonal perturbation-specific residual subspace ($v_5$), and report an exploratory zero-parameter orthogonal subspace fusion alongside two pre-registered negative screening gates on candidate fourth datasets.
4. **Empirical case study of cross-vendor audited autonomous research**: Recent metascience studies have formalized the verification gap in autonomous research agents (Ding et al., 2026), cross-vendor git-native auditing protocols (CrossAudit, 2026), multi-agent critique and falsification benchmarks (Swanson et al., 2025; Bansal and Gangwani, 2026; A Society of AI Scientists, 2026; Dong et al., 2026), and AI-agent preregistration (Vaccaro, 2026). This repository provides a concrete, end-to-end computational biology trace containing a 34-entry public audit ledger of every executor exaggeration, bug, and auditor oversight encountered across eight rounds.

We formalize single-gene perturbation evaluation across paired human cell lines $(A, B)$ around two complementary metric families:

- **Centroid-centered and top-differential-expression reconstruction** ($\Delta r \equiv \mathrm{centroid\_train\_}r$, $\mathrm{centroid\_val\_}r$, $\mathrm{de50\_}r$, $\mathrm{de50\_sign}$): Reconstruction accuracy of perturbation profiles relative to the mean perturbation centroid.
- **Perturbation-specific signal recovery**: (i) orthogonal-residual Pearson correlation ($v_5$), computed on held-out genes after QR-projecting both true and predicted profiles orthogonal to $\mathrm{span}\{\mathbf{1}, \mathbf{b}_c, \mathbf{a}_c\}$, where $\mathbf{b}_c$ is the reproducible training perturbation centroid and $\mathbf{a}_c$ is the dominant shared stress axis (PC1 of training residuals), excluding the perturbed target gene; and (ii) raw-space perturbation discrimination ($\mathrm{disc}_{\mathrm{raw}}$, plus $\mathrm{disc}_{v_5}$), measuring in unprojected space how highly the true target $y_i$ ranks among all $M$ held-out perturbations $\{y_j\}_{j=1}^M$ when retrieved by Pearson correlation with $\hat{y}_i$:
$$\mathrm{disc}_{\mathrm{raw}}(i) = 1 - \frac{1}{M - 1} \sum_{j=1}^M \mathbb{I}\!\left[\mathrm{corr}(\hat{y}_i - \mathbf{b}_c,\, y_j - \mathbf{b}_c) > \mathrm{corr}(\hat{y}_i - \mathbf{b}_c,\, y_i - \mathbf{b}_c)\right] \in [0, 1],$$
where 0.50 corresponds to random retrieval and 1.00 corresponds to exact top-1 retrieval.

---

## 2. Methods

### 2.1 Core Architecture: Soft-Power HybridTOM Kernel and Cross-Line Perturb-Bridge (S5, Primary Method)

All solver constants were permanently frozen at git commit 3609811 on BioXArena-Pro xatlas-orion-perturb prior to any external atlas evaluation. Given paired perturbation matrices across two cell lines $(A, B)$ sharing readout genes $\mathcal{G}_{\mathrm{sh}}$, primary method S5 runs in closed form without gradient descent:

1. **Leak-free training subspace**: Within each split's training set for cell line $c \in \{A, B\}$, set baseline $\mathbf{b}_c$ to the mean of reproducible training perturbations and shared axis $\mathbf{a}_c$ to the unit-norm first right singular vector of centered reproducible training residuals. Project training targets orthogonal to $\mathrm{span}\{\mathbf{1}, \mathbf{b}_c, \mathbf{a}_c\}$.
2. **Self-masked, precision-weighted Soft-Power HybridTOM kernel (S1)**: Construct baseline-subtracted, $\sqrt{n_{\mathrm{cells}}}$-weighted training matrices with each perturbed target gene's self-knockdown entry masked to zero. Compute gene-gene cosine similarity $G^{(0)}$ across training perturbations, then form the signed soft-power HybridTOM kernel ($p = 1.5$):
$$K_{\mathrm{HybridTOM}}(G; p) = \frac{1}{2}\!\left(\mathrm{sign}(G)|G|^p + \mathrm{TOM}\big(\mathrm{sign}(G)|G|^p\big)\right),$$
combining direct gene-gene co-variation with higher-order topological neighborhood overlap.
3. **Multi-scale spectral deflation (S2) and pan-human cross-line gene graph (S3)**: Deflate the top $k \in \{0, 1, 2, 5, 10\}$ principal components of each line's training matrix (weights 1.0, 0.8, 0.5, 0.3, 0.2) and fuse the partner cell line's shared-gene correlation pyramid (weights 0.6, 0.8, 0.8, 0.6, 0.5).
4. **Empirical-SNR non-reproducible anchors (S4)**: Retain non-reproducible measured training perturbations with down-weighted signal-to-noise weights $w_{\mathrm{nr}} = 0.25 \cdot \mathrm{clip}\big((z_{\mathrm{eff}} - 0.9\,\mathrm{med})/(0.3\,\mathrm{med}),\, 0.1,\, 1.5\big)$.
5. **GeneSlope-CosKRR cross-line translator and symmetric bridge (S5)**: Fit a bidirectional nonlinear cross-line translator (skip weight 0.8, ridge regularization $\alpha = 2.0$, kernel power $p_k = 1.5$) on shared training perturbations measured reproducibly in both lines, combining gene-wise shrunk slope skips on shared genes with signed-power cosine kernel ridge regression on 5,000-gene residual profiles. The translator projects partner-only reproducible training perturbations into the target cell line (weight 0.6) and symmetrically fuses translated partner-line KRR predictions ($\alpha = 0.5$) of held-out queries (bridge weight 0.9).

### 2.2 Secondary Extensions (S6, S7) and Exploratory Parameter-Free Orthogonal Fusion

- **Stage S6 (Reciprocal Response Borrowing)**: Borrows residual responses from training perturbations that strongly alter the query gene ($|z| > 4$, quadratic ReLU weight). Although S6 attains the highest raw-space discrimination across all three atlases, its orthogonal-residual behavior is dataset-specific (+0.0130 on xatlas versus −0.0153 on Replogle) and cannot rule out CRISPRi *cis*-proximity confounding on anonymized xatlas.
- **Stage S7 (Output SVD Wiener Denoising and Softmax Manifold Snapping)**: Applies spectral Wiener shrinkage ($\tau = 0.020$) and cosine-similarity manifold snapping (weight 0.45, temperature 18.0, top 40 neighbors). S7 is the exact frozen Codabench submission pipeline; it consistently increases $v_5$ correlation while lowering $\mathrm{disc}_{\mathrm{raw}}$.
- **Exploratory Parameter-Free Orthogonal Decomposition Fusion (OrthoFusion)**: Let $Q_c \in \mathbb{R}^{G \times k}$ ($k \le 3$) denote the orthonormal QR basis of the training shared subspace $S_c = \mathrm{span}\{\mathbf{1}, \mathbf{b}_c, \mathbf{a}_c\}$, with $P_{S_c}(\mathbf{u}) = Q_c(Q_c^\top \mathbf{u})$ and $P_{\perp,c}(\mathbf{u}) = \mathbf{u} - P_{S_c}(\mathbf{u})$. Without any tunable weight,
$$\hat{\mathbf{y}}_{\mathrm{fusion}} = \mathbf{b}_c + P_{S_c}\!\big(\hat{\mathbf{y}}_{\mathrm{scGPT}} - \mathbf{b}_c\big) + P_{\perp,c}\!\big(\hat{\mathbf{y}}_{\mathrm{S5}} - \mathbf{b}_c\big),$$
assigning the shared subspace to the frozen scGPT embedding ridge probe and the orthogonal complement to frozen S5.

### 2.3 Datasets, Pre-Registered Preprocessing, Baseline Protocols, and Statistical Inference

1. **Development atlas (BioXArena-Pro xatlas-orion-perturb)**: HCT116 and HEK293T cell lines (5,000 readout genes per line, 4,040 shared; anonymized gene identifiers never de-anonymized), evaluated on 10 development splits (seeds 42, 0, 1, 2, 3, 7, 11, 99, 2024, 31337), 10 unseen resplits (seeds 101–110), and one frozen Codabench private test submission (455 targets).
2. **External zero-shot atlas 1 (Replogle et al., 2022)**: RPE1 and K562 genome-wide CRISPRi bulk normalized datasets (Figshare+ DOI 10.25452/figshare.plus.20029387), selecting the top 5,000 finite genes by variance-to-mean ratio, defining reproducibility by Anderson-Darling count $\ge 50$ and energy-test $p < 0.001$, and deduplicating to one row per Ensembl target.
3. **Pre-registered unseen confirmatory atlas 2 (Nadig et al., 2025)**: HepG2 (145,473 cells $\times$ 9,624 genes) and Jurkat (262,956 cells $\times$ 8,882 genes) raw single-cell Perturb-seq matrices (NCBI GEO GSE264667, PMID 40259084). Preprocessing was locked prior to inspecting any expression value: per-GEM-group library-size normalization and z-scoring against non-targeting control cells matching Replogle et al. (2022), top 5,000 finite genes by control mean (3,644 shared Ensembl genes), vectorized Scholz-Stephens two-sample Anderson-Darling test (count $\ge 50$ at BH FDR $< 0.05$, unit-verified within $10^{-5}$ against SciPy) plus 20-PC 2,000-permutation energy-distance test ($p < 0.001$) with $n_{\mathrm{cells}} \ge 25$, yielding 2,343 HepG2 (1,000 reproducible) and 2,350 Jurkat (1,262 reproducible) deduplicated Ensembl targets.
4. **Pseudobulk-adapted official GEARS protocol**: Official package cell-gears 0.1.2 was run with default hyperparameters (hidden size 64, 1 GO GNN layer, 1 gene GNN layer, decoder hidden size 16, 20 epochs, learning rate $10^{-3}$, weight decay $5 \times 10^{-4}$, batch size 32) independently per cell line on each split's training set (90%/10% internal early stopping, zero test leakage). Importantly, this completed evaluation supplied one pseudobulk profile per perturbation, 10 jittered copies ($\mathcal{N}(0, 10^{-4})$) of the training control centroid as control cells, and pre-populated train-only top-20 differential gene lists with differential-expression calculation skipped (skip\_calc\_de enabled), rather than native single-cell count distributions. All results from this protocol are therefore labeled **GEARS (official package, pseudobulk-adapted)**.
5. **Frozen scGPT gene-embedding probe**: Official scGPT human checkpoint (33M normal human cells, MIT License) 512-dimensional frozen gene token embeddings paired with ridge regression ($\alpha = 10.0$), evaluated across all 10 splits on Replogle 2022 and Nadig 2025.
6. **Strict simultaneous union holdout and row-level paired bootstrap**: For each seed, 20% of unique reproducible in-readout targets are sampled independently per cell line, their union is formed, and all union targets are excluded simultaneously from both cell lines' training sets (zero cross-line leakage). All metrics are macro-averaged equally across the two cell lines. Because split-level paired $t$-tests are invalid under overlapping training splits, every paired comparison reports 4,000-resample cell-line-stratified row-level paired bootstrap 95% confidence intervals per split and the count of splits (out of 10) whose 95% CI excludes zero.

### 2.4 AI-Assistance and Independent Verification Statement

This study was directed by Zijian Zeng, who holds sole authorship responsibility for scientific oversight, integrity, and all external release decisions. Two AI systems assisted under an explicit separation of duties inspired by cross-vendor audit protocols (CrossAudit, 2026; Ding et al., 2026): Antigravity (Gemini) implemented modeling code, executed experiments, rendered figures from committed reports, and drafted manuscript text; Claude served as an independent cross-vendor technical auditor reviewing pre-registration protocols, checking data-leakage and statistical boundaries, and independently re-executing scripts and recomputing numbers from raw output tables across seven audit checkpoints. Every intervention, retracted claim, and scope downgrade across Rounds 1–8—including an auditor oversight in Round 7 where literature screening was omitted prior to locking a draft thesis—is recorded without omission in the public audit ledger and prior-art verification log (Section 6). Neither AI system is listed as an author, in accordance with ICMJE and journal editorial policies.

---

## 3. Results

### 3.1 Metric-Dependent Ranking Reversal Across Three Perturb-Seq Atlases

Consistent with Miller et al. (2026) and arXiv:2604.27646, method rankings invert systematically across all three human CRISPRi atlases between centroid/differential-expression reconstruction metrics ($\Delta r$, $\mathrm{de50\_}r$, $\mathrm{de50\_sign}$) and perturbation-specificity metrics ($v_5$, $\mathrm{disc}_{\mathrm{raw}}$) (**Table 1**, **Figure 1**).

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
| | GEARS (pseudobulk-adapted) | 0.0486 | — | 0.1077 | 0.6009 | 0.0304 | 0.5621 |
| | **scGPT FrozenEmb Ridge** | 0.1802 | — | 0.2646 | **0.6536** | 0.1127 | 0.6859 |
| | **S5 (Primary)** | 0.1193 | 0.1382 | 0.2170 | 0.6064 | **0.1401** | 0.6884 |
| | OrthoFusion (Exploratory) | **0.1906** | — | **0.2738** | 0.6476 | **0.1401** | **0.6988** |
| | S6 (+Reciprocal) | 0.1106 | 0.1273 | 0.2155 | 0.6020 | 0.1292 | 0.6993* |
| | **S7 (Full)** | 0.1162 | 0.1346 | 0.2013 | 0.5962 | 0.1359 | 0.6736 |
| | S7 no S6 (Exploratory) | 0.1288 | 0.1496 | — | — | 0.1512 | — |
| **Nadig 2025** | Constant Mean | −0.0001 | 0.0292 | 0.0057 | 0.5000 | 0.0000 | 0.5000 |
| (Pre-Reg Unseen) | AE (Linear) | 0.1698 | 0.1787 | 0.2492 | 0.6457 | 0.1146 | 0.6665 |
| | S0 (Ridge) | 0.1610 | 0.1698 | 0.2511 | 0.6451 | 0.1042 | 0.6709 |
| | Simplified GNN | 0.1208 | — | 0.1878 | 0.6090 | 0.1432 | 0.6579 |
| | GEARS (pseudobulk-adapted) | 0.1077 | — | 0.1973 | 0.6462 | 0.0577 | 0.6095 |
| | **scGPT FrozenEmb Ridge** | 0.2129 | — | **0.3064** | **0.6792** | 0.1512 | 0.7109 |
| | **S5 (Primary)** | 0.1502 | 0.1596 | 0.2265 | 0.6293 | 0.1741 | 0.7152 |
| | OrthoFusion (Exploratory) | **0.2144** | — | 0.2899 | 0.6639 | 0.1741 | **0.7324** |
| | S6 (+Reciprocal) | 0.1517 | 0.1612 | 0.2384 | 0.6307 | 0.1756 | 0.7184 |
| | **S7 (Full)** | 0.1617 | 0.1721 | 0.2242 | 0.6277 | 0.1871 | 0.6992 |
| | S7 no S6 (Exploratory) | 0.1644 | 0.1745 | 0.2164 | 0.6298 | **0.1901** | 0.6983 |

: Multi-metric comparison across three independent human CRISPRi Perturb-seq atlases (10-split macro averages). On xatlas, values in parentheses denote strict training-centroid subtraction ($\mathrm{centroid\_train\_}r$). Values marked `*` denote target-gene-masked discrimination percentiles, agreeing with unmasked $\mathrm{disc}_{\mathrm{raw}}$ within 0.002.

![Metric-dependent ranking reversal across three independent human Perturb-seq atlases (xatlas-orion-perturb, Replogle et al. 2022, and pre-registered unseen Nadig et al. 2025). Linear baselines (AE, S0) and frozen scGPT embedding probes lead on centroid-centered reconstruction ($\Delta r$), whereas frozen cross-line kernel S5 and S7 lead on orthogonal residual correlation ($v_5$) and raw perturbation discrimination ($\mathrm{disc}_{\mathrm{raw}}$).](../figures/fig1_metric_inversion_three_datasets.png){width=96%}

### 3.2 Mechanistic PC1 Axis Decomposition Beyond Centroid Subtraction

Direct algebraic audit verified that the baseline vector $\mathbf{b}_c$ equals the reproducible training perturbation centroid across all three datasets—including on xatlas ($\|\mathrm{context\_baseline} - \mu_{\mathrm{rep}}\|_2 = 0.0000$, $r = 1.000000$)—establishing exact identity between $\Delta r$ and Systema training-centroid correlation ($\mathrm{centroid\_train\_}r$). Notably, even after subtracting the training or validation perturbation centroid, linear models (AE, S0) still outperform S5 on $\mathrm{centroid\_val\_}r$ (S5 minus AE: −0.0210 on xatlas, −0.0199 on Replogle, −0.0191 on Nadig).

To test our pre-registered mechanistic hypothesis that this residual linear advantage arises from capturing scalar amplitude along the single dominant shared training residual axis $\mathbf{a}_c$ (PC1), we executed the pre-declared PC1 axis versus off-axis decomposition (pass rule: $\mathrm{off\_axis\_}r(\mathrm{AE} - \mathrm{S5}) \le \frac{1}{3}\,\Delta r(\mathrm{AE} - \mathrm{S5})$) across all three atlases (**Table 2**):

| Atlas | Comparison | Full $\Delta r$ Gap | Off-Axis $r$ Gap | Retention Ratio | Axis $r$ (Ref vs. S5) | Verdict |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **xatlas** | AE − S5 | +0.0533 | +0.0174 | **0.327** ($\le 1/3$) | 0.1870 vs. −0.0974 | **Confirmed** |
| **Replogle 2022** | AE − S5 | +0.0277 | −0.0017 | **−0.063** (sign flip) | 0.3379 vs. 0.0770 | **Confirmed** |
| **Replogle 2022** | scGPT Probe − S5 | +0.0609 | +0.0240 | 0.395 | 0.4270 vs. 0.0770 | Supporting |
| **Nadig 2025** | AE − S5 | +0.0195 | −0.0337 | **−1.725** (sign flip) | 0.4206 vs. −0.0320 | **Confirmed** |
| **Nadig 2025** | scGPT Probe − S5 | +0.0627 | +0.0071 | 0.113 | 0.4774 vs. −0.0320 | Supporting |

: Pre-registered PC1 shared-axis vs. off-axis decomposition across three human CRISPRi atlases (10-split means). Removing the single dominant training residual axis $\mathbf{a}_c$ eliminates $\ge 67\%$ (and over 100% on Replogle and Nadig) of linear AE's $\Delta r$ advantage over S5.

Across all three atlases (3/3), projecting out the single training PC1 axis cuts AE's $\Delta r$ advantage by 67% to 273% (reversing its sign on both Replogle and Nadig) and cuts the frozen scGPT probe's $\Delta r$ advantage by 61% on Replogle and 89% on Nadig.

### 3.3 Pre-Registered Unseen External Confirmation on Nadig et al. (2025)

Prior to reading any Nadig et al. (2025) expression values, we locked the frozen primary method (S5), preprocessing rules, and quantitative pass criteria in git. All three pre-registered hypotheses were confirmed on unseen HepG2 and Jurkat Perturb-seq profiles (**Figure 2**; seeds 42, 0, and 1 were additionally re-executed verbatim by the independent auditor):

- **H1 (Primary — Raw perturbation discrimination $\mathrm{disc}_{\mathrm{raw}}$)**: S5 minus AE = +0.0487 (0.7152 vs. 0.6665), mean 95% CI [+0.0161, +0.0812], with **9/10 splits strictly above zero** (confirmed). *Caveat*: Several split-level lower bounds lie close to zero (+0.0011, +0.0024, +0.0038, and −0.0080 on seed 2024), indicating a consistent but modest raw-space margin.
- **H2 (Primary — Metric inversion on $\Delta r$)**: S5 minus AE = −0.0195 (0.1502 vs. 0.1698), with **9/10 splits negative**, **4/10 splits with 95% CI strictly below zero**, and 0/10 above zero (confirmed). *Caveat*: Direction is consistent, though per-split significance (4/10) is weaker than on Replogle 2022 (8/10).
- **Joint primary thesis (H1 and H2)**: Confirmed.
- **H3 (Secondary — Orthogonal residual $v_5$)**: S5 minus AE = +0.0595 (0.1741 vs. 0.1146), mean 95% CI [+0.0437, +0.0753], with **10/10 splits strictly above zero** (confirmed).

Frozen S5 also outperformed simplified neural baselines on Nadig 2025 (vs. Simplified GNN: $v_5$ +0.0309 [9/10], $\mathrm{disc}_{\mathrm{raw}}$ +0.0572 [10/10], $\Delta r$ +0.0295 [9/10]; vs. Simplified ResMLP: $v_5$ +0.0449 [10/10], $\mathrm{disc}_{\mathrm{raw}}$ +0.0529 [10/10]). As observed on xatlas and Replogle, S6 attained the highest $\mathrm{disc}_{\mathrm{raw}}$ (0.7184), whereas S7 raised $v_5$ (0.1871, 10/10 vs. AE) at the cost of lower $\mathrm{disc}_{\mathrm{raw}}$ (0.6992, reducing S7 minus AE to 4/10 splits CI > 0).

![Pre-registered zero-look confirmatory forest plot on unseen Nadig et al. (2025; HepG2 and Jurkat) showing split-level 4,000-resample stratified row-level paired bootstrap 95% confidence intervals for S5 minus AE across $\mathrm{disc}_{\mathrm{raw}}$ (H1, 9/10 CI > 0), $\Delta r$ (H2, 9/10 negative, 4/10 CI < 0), and $v_5$ (H3, 10/10 CI > 0).](../figures/fig2_preregistered_confirmation_forest.png){width=96%}

### 3.4 External Deep Learning, Foundation-Model Probe, and Exploratory Orthogonal Fusion

1. **GEARS (official package, pseudobulk-adapted)**: Across all 10 splits, official cell-gears 0.1.2 trained under our pseudobulk-adapted protocol achieved $v_5 = 0.0304 \pm 0.0168$, $\Delta r = 0.0486 \pm 0.0253$, and $\mathrm{disc}_{\mathrm{raw}} = 0.5621 \pm 0.0208$ on Replogle 2022 (GO-graph coverage 99.83% RPE1 / 99.38% K562), and $v_5 = 0.0577 \pm 0.0143$, $\Delta r = 0.1077 \pm 0.0174$, $\mathrm{disc}_{\mathrm{raw}} = 0.6095 \pm 0.0149$, and $\mathrm{de50\_sign} = 0.6462 \pm 0.0095$ on Nadig 2025 (GO-graph coverage 99.66% HepG2 / 99.79% Jurkat). Under this pseudobulk-adapted input regime, GEARS trailed both AE and S5 on every metric; because pseudobulk adaptation omits native single-cell dispersion and differential-expression selection, a native single-cell GEARS audit on Nadig 2025 is in progress and excluded from current claims (Section 5).
2. **Frozen scGPT gene-embedding probe**: Consistent with Miller et al. (2026), scGPT FrozenEmb Ridge (512-dim frozen token embeddings with ridge regression $\alpha = 10.0$) is our strongest single baseline, significantly outperforming linear AE on every metric across both Replogle ($\Delta r$ +0.0332 [8/10], $\mathrm{disc}_{\mathrm{raw}}$ +0.0392 [8/10], $\mathrm{disc}_{v_5}$ +0.0420 [10/10], $v_5$ +0.0160 [6/10]) and Nadig ($\Delta r$ +0.0432 [10/10], $\mathrm{disc}_{\mathrm{raw}}$ +0.0444 [8/10], $\mathrm{disc}_{v_5}$ +0.0591 [10/10], $v_5$ +0.0365 [10/10]). Directly comparing zero-pretraining S5 against this 33M-cell foundation-model probe reveals a sharp geometric boundary: **S5 significantly outperforms scGPT FrozenEmb Ridge only on orthogonal residuals $v_5$** (+0.0274, 10/10 CI > 0 on Replogle; +0.0230, 7/10 CI > 0 on Nadig), **statistically ties on $\mathrm{disc}_{\mathrm{raw}}$** (+0.0025 Replogle, +0.0043 Nadig, 0/10 splits significant in either direction), and **trails substantially on $\Delta r$** (−0.0609 Replogle, −0.0627 Nadig, 10/10 CI < 0).
3. **Exploratory parameter-free orthogonal fusion and pre-registered fourth-dataset hard stop**: Applying zero-weight OrthoFusion strictly as an exploratory test on the already-seen Replogle and Nadig atlases yields numerical Pareto gains over both parent models: on Replogle 2022, $\Delta r = 0.1906$ (+0.0104 vs. scGPT probe, 0/10 CI < 0), $v_5 = 0.1401$ (exact algebraic parity with S5, maximum absolute deviation $< 10^{-12}$), and $\mathrm{disc}_{\mathrm{raw}} = 0.6988$ (+0.0128 vs. scGPT probe [5/10 CI > 0]; +0.0103 vs. S5 [0/10 CI > 0]); on Nadig 2025, $\Delta r = 0.2144$ (+0.0015 vs. scGPT probe, 0/10 CI < 0), $v_5 = 0.1741$ (exact parity with S5), and $\mathrm{disc}_{\mathrm{raw}} = 0.7324$ (+0.0216 vs. scGPT probe [7/10 CI > 0]; +0.0172 vs. S5 [3/10 CI > 0]). This satisfies exploratory non-inferiority (H4) on both atlases but **does not meet** our pre-declared dual-sided synergy criterion (H5: $\ge 7/10$ splits CI > 0 against both parents simultaneously). Furthermore, when screening independent non-Orion public multi-line CRISPRi datasets under our locked statistical power gate ($\ge 100$ reproducible targets per line, $\ge 80$ shared reproducible targets, $\ge 2{,}000$ shared readout genes), Candidate 1 (Feng et al., 2025; Figshare 27989294 TargetedScreen, top-2 individual donor iPSC lines jejf\_2 and zapk\_3) yielded only 38 (27 in readout) and 86 (62 in readout) reproducible targets (37 shared) at Day 3 harvest, and Candidate 2 (Jiang et al., 2025; Zenodo 14518762) contained only 44–61 targets per matched cytokine context (IFNB 61, IFNG 59, TNFA 55, TGFB1 52, INS 44). Enforcing strict pre-registration discipline, we halted prior to running any predictive model and record fourth-dataset confirmation of OrthoFusion as **pending**.

### 3.5 Leave-One-Module-Out Transferability and ±50% Sensitivity Analysis

Dual-dataset leave-one-module-out (LOMO) ablation (**Figure 3A**) separates transferable components from dataset-specific components:

- **Transferable modules (positive on both xatlas and Replogle)**: S1 Soft-Power HybridTOM (+0.0109 on xatlas [7/10 CI > 0] / +0.0103 on Replogle [10/10 CI > 0]), S7 Wiener snapping (+0.0072 / +0.0067 on $v_5$, while reducing $\mathrm{disc}_{\mathrm{raw}}$), S5 Cross-Line Bridge (+0.0020 / +0.0061, Replogle 9/10 CI > 0), and S2 multi-scale spectral deflation (+0.0108 / +0.0035, Replogle 1/10 CI > 0).
- **Dataset-specific modules**: S6 Reciprocal Borrowing (+0.0130 [8/10 CI > 0] on xatlas versus **−0.0153 [9/10 CI < 0]** on Replogle; omitting S6 raises Replogle $v_5$ to 0.1512) and S4 non-reproducible anchors (+0.0043 vs. −0.0006).

One-factor-at-a-time ±50% perturbation of all 10 frozen constants in primary method S5 (**Figure 3B**) demonstrates broad stability across both external atlases: 18/20 perturbed configurations on Replogle and 20/20 on Nadig retain $v_5$ 10/10 splits CI > 0 against AE, with 17–18/20 retaining $\mathrm{disc}_{\mathrm{raw}}$ superiority. Crucially, the sweep exposed a sharp **concave-power regime boundary ($p < 1$)**: reducing the translator kernel power $p_k$ by 50% from 1.5 to 0.75 collapses Replogle performance ($v_5 = 0.0602$, shift −0.0798; $\mathrm{disc}_{\mathrm{raw}} = 0.6052$, 0/10 CI > 0) and causes the largest Nadig drop ($v_5 = 0.1571$, −0.0170). Spectral diagnostics on Replogle seed 42 pinpoint the exact mathematical cause: concave element-wise power $\mathrm{sign}(G)|G|^{0.75}$ breaks kernel positive semidefiniteness ($\lambda_{\min}(K) = -6.177$), introducing six negative eigenvalues in $K + 2I$ ($\lambda_{\min} = -4.177$, $\min|\lambda| = 0.0267$) and inflating the linear-system condition number from 17.2 ($p_k = 1.5$) to 4,473.5 ($p_k = 0.75$). Every convex power $p_k \ge 1.0$ preserves strict positive definiteness (zero negative eigenvalues; condition number 37.5 at 1.0, 17.2 at 1.5, 7.2 at 2.25). Per protocol, all frozen constants remained untouched.

![Leave-one-module-out (LOMO) transferability ablation across xatlas and Replogle 2022 (A) and one-factor-at-a-time ±50% hyperparameter sensitivity sweep of all 10 frozen S5 constants on Replogle 2022 and Nadig 2025 (B), exposing the concave-power ($p_k < 1$) non-PSD instability boundary.](../figures/fig3_lomo_ablation_and_s5_sensitivity.png){width=96%}

### 3.6 Biological Pathway Enrichment and CRISPRi *Cis*-Proximity Audit

Evaluating local Reactome pathway enrichment strictly over the GMT-annotated evaluated gene background ($M_{\mathrm{bg}} = 1{,}452$, **Figure 4A**):

- Top-150 S7 minus S0 gain genes (128 annotated in Reactome GMT) significantly enrich **mRNA Splicing – Major Pathway** (38/145 genes, 2.97-fold enrichment, BH FDR = $2.90 \times 10^{-8}$; 0/50 random 150-gene controls passing FDR $< 0.05$).
- Top-150 S5 minus S2 cross-line bridge gain genes (131 annotated in GMT) selectively enrich **Mitochondrial Translation** (20/82 genes, 2.70-fold enrichment, BH FDR = $6.98 \times 10^{-3}$; 0/50 random controls).
- Sensitivity analysis on the full evaluated gene universe ($M_{\mathrm{bg}} = 1{,}784$) yields consistent enrichments (mRNA Splicing: 3.12-fold, FDR = $1.18 \times 10^{-8}$; Mitochondrial Translation: 2.90-fold, FDR = $2.76 \times 10^{-3}$).

Auditing S6 reciprocal donor pairs against GRCh38 GENCODE v46 coordinates (10,530/10,539 Replogle genes mapped, **Figure 4B**) shows that negative-$z$ pairs ($z < -4$) are enriched 474-fold for $< 10\ \mathrm{kb}$ genomic proximity (1.422% vs. 0.003% null; 3.32% of negative-$z$ weight), yet $< 10\ \mathrm{kb}$ pairs constitute only 0.255% of total active $|z| > 4$ weight, and removing all $< 10\ \mathrm{kb}$ pairs retains 97.41% of Top-40 own-reciprocal gains on Replogle. On anonymized xatlas—where 78% of Top-40 S6 winners draw from negative-$z$ top donors—local CRISPRi *cis*-repression cannot be ruled out.

![Reactome pathway enrichment of S7 minus S0 (mRNA Splicing – Major Pathway) and S5 minus S2 (Mitochondrial Translation) gain genes over the GMT-annotated evaluated background ($M_{\mathrm{bg}} = 1{,}452$) (A), and GRCh38 GENCODE v46 genomic-distance audit of S6 reciprocal donors on Replogle 2022 (B).](../figures/fig4_biology_reactome_and_cis_audit.png){width=96%}

### 3.7 Official Codabench Private Test Evaluation

A single frozen submission on BioXArena-Pro xatlas-orion-perturb (Submission 1062221, panperturb-bridge human-ai-collab, status ok) scored **0.0998** on the blind official private test set—+0.0173 above the best 2-hour autonomous agent score of 0.0825 (within the 0.02 benchmark noise margin; validation-to-test shift $0.0998 - 0.1099 = -0.0101$). Because the frozen pipeline reflects multi-round iterative human-AI development (~25 fixed constants) rather than single-shot 2-hour autonomous execution, it is reported strictly as a separate development category.

---

## 4. Discussion

Rather than claiming metric-dependent ranking reversal or cross-cell-line transfer as new discoveries—both of which are established in recent literature (Miller et al., 2026; arXiv:2604.27646; Adduri et al., 2025; HyperMap, 2026; TranScouter, 2026)—this study establishes three narrowly bounded empirical results across three independent Perturb-seq laboratories. First, under commit-hashed zero-look pre-registration and cell-line-stratified row-level bootstrap inference, the ranking reversal between centroid-reconstruction metrics and perturbation-specificity metrics reproduces on unseen Nadig et al. (2025) data. Second, our pre-registered PC1 decomposition demonstrates that even after Systema training-centroid centering, linear models' advantage over cross-line residual kernels resides predominantly along the single shared training PC1 variation axis: projecting out that axis removes 67%–273% of AE's $\Delta r$ lead on all three atlases. Third, comparing a zero-pretraining cross-line kernel (S5) against a 33M-cell frozen scGPT embedding probe isolates a clean geometric complementarity—the foundation-model probe dominates shared-subspace reconstruction ($\Delta r$), whereas S5 significantly outperforms it on orthogonal perturbation-specific residuals ($v_5$) while tying on raw discrimination ($\mathrm{disc}_{\mathrm{raw}}$). Exploratory parameter-free orthogonal fusion (OrthoFusion) reconciles both regimes numerically on seen atlases without tuned weights, but prospective fourth-dataset confirmation awaits a publicly available non-Orion multi-cell-line CRISPRi atlas satisfying our locked statistical power gate.

From a metascience perspective, our eight-round execution trace illustrates both the value and the boundaries of cross-vendor AI auditing (CrossAudit, 2026; Ding et al., 2026; A Society of AI Scientists, 2026; Dong et al., 2026). Independent code re-execution by a cross-vendor auditor caught 34 distinct failure modes—retracting 13 inflated early claims (such as lucky-seed selection, cross-dataset comparisons, and invalid split-level $t$-tests) and enforcing negative-gate hard stops on two underpowered fourth datasets. Equally importantly, our public audit ledger records a shared blind spot: omitting systematic prior-art search before Round 7 allowed an already-published metric-inversion narrative (Miller et al., 2026) to survive until Round 8, motivating the formal prior-art verification gate enforced here.

---

## 5. Limitations and In-Progress Audits

1. **Prior art on core themes**: Metric-dependent ranking reversal and frozen embedding-probe strength were previously established by Miller et al. (2026) and arXiv:2604.27646; cross-context perturbation transfer exists in STATE (2025), HyperMap (2026), and TranScouter (2026); and cross-vendor git-native AI auditing is formalized in CrossAudit (2026) and related metascience studies. No priority ("first") claim is made for these paradigms.
2. **Performance against frozen foundation-model probes and traditional metrics**: Standalone frozen S5 trails AE, S0, and scGPT FrozenEmb Ridge on $\Delta r$, $\mathrm{centroid\_val\_}r$, $\mathrm{de50\_}r$, and $\mathrm{de50\_sign}$. Against scGPT FrozenEmb Ridge, S5 wins only on $v_5$, statistically ties on $\mathrm{disc}_{\mathrm{raw}}$ (0/10 splits significant), and trails substantially on $\Delta r$ (10/10 CI < 0).
3. **Unconfirmed fourth-dataset validation for OrthoFusion**: OrthoFusion results on Replogle 2022 and Nadig 2025 are exploratory (seen datasets) and do not meet dual-sided H5 significance ($\ge 7/10$ CI > 0 against both parents); fourth-dataset confirmation remains pending after Feng et al. (2025) and Jiang et al. (2025) failed locked sample-size gates.
4. **In-progress compute audits (excluded from Preprint v1 claims)**: Two follow-up audits are currently running and intentionally excluded from v1 results: (i) a native single-cell GEARS evaluation on Nadig 2025 with native differential-expression calculation enabled to replace the pseudobulk-adapted GEARS caveat, and (ii) Miller et al. (2026) split-half Dynamic Range Fraction (DRF) metric calibration and model rankings across Replogle 2022 and Nadig 2025. Both will be incorporated verbatim—regardless of outcome—in Preprint v2.
5. **Compute and architecture constraints for full foundation-model fine-tuning**: Full 10-split scGPT TransformerGenerator fine-tuning (~100 hours on Apple Silicon MPS without CUDA FlashAttention; an exploratory 2-seed Nadig fine-tuning run was discontinued after 2h 41m under compute prioritization) and Arc Institute STATE (defaulting to one-hot perturbation representations incompatible with union holdout of unseen targets without non-commercial SE-600M embeddings and multi-GPU CUDA) were not completed.
6. **Potential scGPT pretraining corpus overlap**: scGPT human was pretrained on 33M normal human cells from CELLxGENE (2023). Although Nadig et al. (2025; GEO release May 2024) post-dates scGPT pretraining, Replogle et al. (2022) pre-dates 2023 and broad co-expression overlap cannot be ruled out.
7. **Other empirical caveats**: The Codabench private test score (0.0998) sits within 0.02 of the best 2-hour autonomous agent (0.0825) and belongs to a different development category; several Nadig H1 split lower bounds lie near zero; S6 cannot rule out CRISPRi *cis*-repression on anonymized xatlas; and concave kernel power ($p_k < 1$) violates positive semidefiniteness.

---

## 6. Data, Code, and Audit Availability

All code, frozen manifests, pre-registration protocols, prior-art logs, and audit artifacts are publicly archived at `https://github.com/ArtificialZeng/AI-Biology-Lab` under `projects/bioxarena-pro-scaffold/`, with zero hardcoded local file paths (configured via `BIOXARENA_DIR`, `REPLOGLE_DATA_DIR`, `NADIG_DATA_DIR`, `FENG_DATA_DIR`, and `SCGPT_WEIGHTS_DIR`) and a locked dependency specification (`requirements-lock.txt`):

- **Frozen Solver & Codabench Pack**: `src/xatlas_panhuman_solver.py` (commit `3609811`), `outputs/perturbation_baselines/xatlas-orion-perturb/CODABENCH_FREEZE_MANIFEST.json` (`submission.csv` SHA-256 `e2613a9b6616368ba5ca29edf132fe023761f0f72072b25176724207daa6a207`), and the official Codabench Submission 1062221 result archive.
- **Pre-Registration, Prior-Art Gate & Audit Ledger**: `reports/ROUND7_PREREGISTRATION.md` (commits `c5294c0`, `87f4bc2`, `b9fe2bb`), `reports/ROUND8_PREREGISTRATION.md` (commits `ce82040`, `840b502`), `reports/PRIOR_ART_LOG.md`, `reports/AUDIT_LEDGER.csv`, and `reports/AUDIT_LEDGER_SUMMARY.md`.
- **Result Summary Files**: Multi-metric summaries (`reports/t2_all_5metric_summary.json`, `reports/t2_discrim_delta_row_ci_summary.json`, `reports/r74_systema_*_summary.json`), Nadig pre-registered confirmation (`reports/r71_nadig_confirmation_summary.json`), PC1 axis decomposition (`reports/r82_axis_decomposition_*_summary.json`), external baselines (`reports/r72_official_gears_*_summary.json`, `reports/r72_scgpt_probe_*`), exploratory fusion & fourth-dataset gate (`reports/r83_exploratory_orthofusion_*_summary.json`, `reports/r83_confirmatory_feng2025_summary.json`), LOMO & sensitivity (`reports/t4_lomo_summary.json`, `reports/r73_sensitivity_*_summary.json`), and biology/cis audits (`reports/t5_local_reactome_enrichment_summary.json`, `reports/t6_replogle_crispri_proximity_summary.json`).
- **Public External Data & Checkpoints**: Replogle et al. (2022) via Figshare+ DOI `10.25452/figshare.plus.20029387` (`rpe1_normalized_bulk_01.h5ad` SHA-256 `a3c5bfd0...`, `K562_gwps_normalized_bulk_01.h5ad` SHA-256 `37e48c47...`); Nadig et al. (2025) via NCBI GEO `GSE264667` (`hepg2` SHA-256 `e1ad7c3c...`, `jurkat` SHA-256 `ffbe15f2...`); Feng et al. (2025, halted at QC gate) via Figshare DOI `10.6084/m9.figshare.27989294.v2`; Jiang et al. (2025, halted at target-count audit) via Zenodo `14518762`; `scGPT_human` checkpoint (MIT License, SHA-256 `ad0252a1...`) and vocabulary (`ee2b2c90...`); Reactome Pathways GMT (CC BY 4.0, SHA-256 `a9feaa3a...`); GENCODE v46 Basic GTF (SHA-256 `d620d548...`).

---

## References

1. Adduri, A. K., Roohani, Y., Hsu, P. et al. (2025). STATE: Predicting cellular perturbation responses across contexts. *Arc Institute / NeurIPS 2025*.
2. Ahlmann-Eltze, C., Huber, W. & Anders, S. (2025). Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines. *Nature Methods*, August 2025.
3. Bansal, A. & Gangwani, I. (2026). Propose-Critique-Falsify: Benchmarking multi-agent scientific claim verification. *ICML AI4Science Workshop 2026*.
4. CrossAudit Authors (2026). CrossAudit: A Git-Native, Cross-Vendor Audit Loop for Agentic Science. *arXiv:2608.28631*.
5. Cui, H., Wang, C., Maan, H. et al. (2024). scGPT: toward building a foundation model for single-cell multi-omics using generative AI. *Nature Methods* 21, 1470–1480.
6. Ding, T., Nannapaneni, A., Liu, B. & Zhang, L. (2026). Autonomous Research Agents: A Survey of AI Scientists and the Verification Gap. *arXiv:2608.05179*.
7. Dong, Y., Zhang, L. & Zheng, T. (2026). Author-in-the-loop Agentic Reproducibility Auditing for a Machine Learning Conference. *arXiv:2610.22465*.
8. Feng, C. et al. (2025). Targeted Perturb-seq across human iPSC donor lines. *Cell Genomics* (2025); PMID 41330380; Figshare DOI `10.6084/m9.figshare.27989294.v2`.
9. GRAFT Authors (2026). Cross-modal transfer enables zero-shot prediction of drug responses across cell types. *bioRxiv*, November 2026.
10. HyperMap Authors (2026). HyperMap: An Efficient Framework for Transferring Perturbation Responses Across Diverse Biological Contexts. *bioRxiv*, April 2026.
11. Jiang, L., Dalgarno, C. et al. (2025). Perturb-seq across cytokine stimulation contexts. *Nature Cell Biology* (2025); GEO `GSE281048`; Zenodo `14518762`.
12. Miller, H. et al. (2026). Deep learning perturbation models can outperform baselines on calibrated metrics. *Nature Biotechnology* (2026); *bioRxiv* `10.1101/2025.10.20.683304`; *arXiv:2604.27646*.
13. Nadig, A. et al. (2025). Transcriptome-wide characterization of genetic perturbation responses across human cell lines. *Nature Genetics* (2025); GEO `GSE264667`; PMID 40259084.
14. Replogle, J. M., Saunders, R. A., Pogson, A. N. et al. (2022). Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq. *Cell* 185, 2559–2575.
15. Roohani, Y., Huang, K. & Leskovec, J. (2024). Predicting transcriptional outcomes of novel multigene perturbations with GEARS. *Nature Biotechnology* 42, 927–935.
16. Society of AI Scientists Consortium (2026). A Society of AI Scientists. *alphaXiv / Stanford / MIT / ETH Zurich*, November 2026.
17. Swanson, K., Zou, J., Pak, J. E. et al. (2025). The Virtual Lab of AI agents designs new SARS-CoV-2 nanobodies. *Nature* (2025).
18. TranScouter Authors (2026). A structured study of cross-condition prediction of transcriptional responses to gene perturbations. *bioRxiv*, August 2026.
19. Vaccaro, M. (2026). Preregistration for Experiments with AI Agents. *arXiv:2606.11217*.
20. Viñas Torné, R., Wiatrak, M., Brbić, M. et al. (2025/2026). Systema: a framework for evaluating genetic perturbation response prediction beyond systematic variation. *Nature Biotechnology* (2025/2026).
