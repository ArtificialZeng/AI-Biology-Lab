# PanPerturb-Bridge: Multi-Scale Spectral Deflation, Cross-Line Causal Bridging, and Reciprocal Perturbation Law for Unseen CRISPRi Response Prediction

**Working Research Manuscript (Level 2 Confirmed Milestone)**  
**Target Benchmark Task**: BioXArena-Pro `perturbation/xatlas-orion-perturb`  
**Primary Executable Implementation**: `src/xatlas_panhuman_solver.py` (exact code parity with `outputs/perturbation_baselines/xatlas-orion-perturb/solution.py`)

---

## Abstract

Predicting genome-wide transcriptional responses to unseen CRISPRi gene knockdowns remains one of the hardest open challenges in computational systems biology, particularly when generic cellular stress responses are orthogonally projected out (as in the BioXArena-Pro `xatlas-orion-perturb` v5 evaluation protocol). Standard spectral baselines (`100-SVD + Ridge`) suffer from four structural blind spots: (i) training-only self-knockdown spikes (`6.8–8.6σ`) corrupt co-expression geometry; (ii) dominant generic stress modes mask fine-grained pathway regulons; (iii) cell-line-specific models discard both the 4,040 shared readout genes across human cell lines (covering 84.6% of unseen test targets) and 2,500+ non-reproducible knockdown measurements carrying recoverable signal-to-noise; and (iv) all standard regression models completely ignore 4,136 reproducible genome-wide knockdowns whose target genes lie outside the 5,000-gene readout panel.

We introduce **PanPerturb-Bridge**, a unified, cell-line-agnostic framework combining seven complementary, biologically grounded components: (1) self-spike masking, baseline subtraction, $\sqrt{n_{\text{cells}}}$ inverse-noise weighting, and soft-thresholded Hybrid 1-hop/2-hop Topological Overlap (`HybridTOM`, $p=1.5$); (2) a Multi-Scale Spectral Deflation Pyramid (`d0, d1, d2, d5, d10`); (3) Pan-Human cross-line gene-graph denoising; (4) Empirical-SNR weighted un-normalized Kernel Ridge Regression incorporating non-reproducible measured targets; (5) a `GeneSlope + Soft-Power CosKRR` cross-line perturbation translator (5-fold CV `0.3783 / 0.3551` on 781 paired reproducible knockdowns) powering a symmetric bidirectional bridge; (6) a Thresholded Reciprocal Causal Perturbation Law ($\mathrm{ReLU}(|z_{k,g}| - 4.0)^2$) unlocking all 6,018 reproducible knockdowns—including the 4,136 unmeasured targets; and (7) output SVD Wiener shrinkage with empirical perturbation manifold snapping.

Across 10 independent random splits under the official BioXArena-Pro `strict_group_mean_from_row_scores` grader (including 5 untouched held-out confirmation splits never accessed during development), **PanPerturb-Bridge** raises mean projected Pearson correlation from `0.0613 ± 0.0040` (`100-SVD + Ridge`) to **`0.1099 ± 0.0087`** (**`+0.0485 ± 0.0085`**, `+79.3%` relative, `10/10` positive splits, paired $t = 17.08$, $p = 3.63 \times 10^{-8}$, `2.43×` the official `0.02` noise margin), with every single stage monotonically improving both HCT116 (`0.0767 → 0.1462`), HEK293T (`0.0460 → 0.0736`), Dev-5 (`0.0606 → 0.1134`), and Untouched Held-5 (`0.0621 → 0.1064`). Under an ultra-strict audit that simultaneously removes the union of held-out validation targets (`~360` genes) from *both* cell lines and refits `context_baseline` and `residual_axis` strictly on remaining training rows, **PanPerturb-Bridge** achieves **`0.1115 ± 0.0089`** vs `0.0658 ± 0.0033` (**`+0.0457 ± 0.0081`**, `10/10`, $t = 17.01$, $p = 3.77 \times 10^{-8}$; Untouched Held-5 `0.1097`).

---

## 1. Problem Formulation & Evaluation Protocol

For each human cell line $c \in \{\text{HCT116}, \text{HEK293T}\}$, the dataset provides pseudo-bulk transcriptional log-fold-change vectors $y_k \in \mathbb{R}^{5000}$ across 5,000 readout genes for $N_c$ CRISPRi knockdowns (`reproducible` $\in \{0,1\}$, cell count $n_k$), along with a line-specific scoring gene subset $\mathcal{S}_c$, context baseline $b_c$, and dominant residual axis $a_c$ (the first right-singular vector of unit-normalized, $\mathrm{span}\{\mathbf{1}, b_c\}$-projected reproducible training perturbations).

Under the official `xatlas-orion-perturb` v5 grader (`tasks/perturbation/xatlas-orion-perturb/grade.py`), for each held-out target gene $g$ (excluding $g$ and any private cis exclusions from $\mathcal{S}_c$), let $Q$ be an orthonormal basis for $\mathrm{span}\{\mathbf{1}, b_c, a_c\}$ on the active scoring indices, and $P_Q v = v - Q Q^\top v$. The row score is:

$$r(y, \hat{y}) = \mathrm{clip}\!\left(\frac{\langle P_Q(y - b_c),\, P_Q(\hat{y} - b_c) \rangle}{\|P_Q(y - b_c)\|_2 \, \|P_Q(\hat{y} - b_c)\|_2},\, -1,\, 1\right)$$

with a non-degeneracy gate returning $0$ whenever $\|P_Q(\hat{y}-b_c)\|_2 \le 10^{-4}(\|\hat{y}\|_2 + \|b_c\|_2) + 10^{-12}$. Row scores are aggregated via official `tasks.utils.strict_group_mean_from_row_scores` (mean within each cell line, then equal-weighted mean across HCT116 and HEK293T).

**Zero Test-Target Leakage Property**: In the public benchmark, 4,040 of 5,000 readout genes are shared between HCT116 and HEK293T, covering 385 of 455 test targets (84.6%), while **0 of 455 test targets** appear in `train.csv` for either cell line. All validation protocols strictly enforce this invariant.

---

## 2. Methodology: PanPerturb-Bridge

All hyperparameters are strictly cell-line-agnostic (identical across HCT116 and HEK293T).

### Stage 1: Self-Spike Masking, Inverse-Noise Weighting, and HybridTOM ($p=1.5$)
In `train_y`, direct CRISPRi knockdown creates a massive negative spike on the targeted gene itself (`-8.6σ` in HCT116, `-6.8σ` in HEK293T) that never occurs for unseen test target columns. We zero out `Y[i, idx(p_i)] = 0`, subtract $b_c$, and weight each perturbation row $i$ by $w_i = \sqrt{n_i} / \overline{\sqrt{n}}$ (inverse shot-noise standard deviation) to form $\tilde{Y}_c$. Let $G_c^{(0)}$ be the column-normalized cosine correlation matrix of $\tilde{Y}_c$. Retaining the unit diagonal inside the signed soft-thresholded adjacency $A = \mathrm{sign}(G) |G|^{p}$ ($p=1.5$) and taking row-normalized inner products $K = \hat{A}\hat{A}^\top$ yields a **Hybrid 1-hop + 2-hop Topological Overlap Matrix (`HybridTOM`)** preserving direct heterodimer edges alongside shared-neighbor co-regulation.

### Stage 2: Multi-Scale Spectral Deflation Pyramid (`d0, d1, d2, d5, d10`)
A single generic stress mode dominates transcriptional variance (`57.1%` in HEK293T), causing standard correlation kernels to cluster perturbations along stress magnitude rather than specific pathway identity—precisely the direction projected out by the v5 grader. Using the top-10 SVD of $\tilde{Y}_c = U \Sigma V^\top$, we construct deflated residual matrices $\tilde{Y}_c^{(k)} = \tilde{Y}_c - U_{:,:k} \Sigma_{:k} V_{:k,:}^\top$ for $k \in \{0, 1, 2, 5, 10\}$ and combine their `HybridTOM` kernels with decaying weights $\{1.0, 0.8, 0.5, 0.3, 0.2\}$.

### Stage 3: Pan-Human Cross-Line Gene Graph Fusion
For the 4,040 shared readout genes, we fuse cell-line-specific and complementary-line deflated correlations *prior* to nonlinear soft-thresholding:

$$G_{\text{pan}, c}^{(k)}[\mathcal{I}_{\text{sh}}, \mathcal{I}_{\text{sh}}] = G_c^{(k)}[\mathcal{I}_{\text{sh}}, \mathcal{I}_{\text{sh}}] + w_{\text{pan}}^{(k)} G_{\bar{c}}^{(k)}[\mathcal{I}'_{\text{sh}}, \mathcal{I}'_{\text{sh}}]$$

with $w_{\text{pan}} = \{0.6, 0.8, 0.8, 0.6, 0.5\}$. Fusing before soft-power suppression reinforces conserved human gene modules while cancelling cell-line-specific measurement noise.

### Stage 4: Empirical-SNR Weighted Un-Normalized Kernel Ridge Regression
Over 2,500 measured targets are labeled `reproducible=0` due to modest cell count or weaker effect size, yet many exhibit substantial empirical signal-to-noise $Z_{\text{eff}}(i) = \|P_Q(y_i - b_c)\|_2 \sqrt{n_i}$. We incorporate them via weighted centered KRR with **un-normalized sample weights**:

$$\hat{Y}_{\text{val}} = \bar{y}_w + K_{\text{val},c} \left(K_{\text{tr},c} + \mathrm{diag}\!\left(\frac{\alpha}{w_i}\right)\right)^{-1} (Y_{\text{tr}} - \bar{y}_w)$$

Keeping reproducible weights fixed at $w_i = 1.0$ (never normalizing $w$ by $\bar{w}$) guarantees that reproducible anchors retain exact regularization $\alpha = 0.5$, while noisy non-reproducible anchors receive heavier shrinkage $\alpha / w_i$.

### Stage 5: GeneSlope-CosKRR Cross-Line Translator & Symmetric Bidirectional Bridge
Using paired reproducible perturbations observed in both cell lines (strictly excluding validation targets; $N = 781$ on full train), we learn bidirectional response translators $\mathcal{T}_{\bar{c} \to c}$ combining an analytical gene-wise ridge slope skip connection on the 4,040 shared genes with a soft-power ($p=1.5$) cosine KRR residual over all 5,000 genes. On 5-fold CV over the 781 paired reproducible knockdowns (`src/xatlas_translator_cv.py`), $\mathcal{T}$ achieves `0.3783 ± 0.2635` (`HEK293T→HCT116`, vs `0.2376` raw copy) and `0.3551 ± 0.2431` (`HCT116→HEK293T`, vs `0.2420` raw copy). This translator powers (a) augmenting each line's training targets with translated reproducible perturbations from the complementary line ($w=0.6$), and (b) a **Symmetric Bidirectional Bridge** that predicts each shared query target in $\bar{c}$, translates $\mathcal{T}_{\bar{c}\to c}(\hat{y}_{\bar{c}})$ into line $c$, and fuses it with weight `0.9`.

### Stage 6: Thresholded Reciprocal Causal Perturbation Law (Unlocking 4,136 Unmeasured Knockdowns)
A critical structural feature of `xatlas-orion-perturb` is that **4,136 reproducible knockdowns** (`668` in HCT116, `3,468` in HEK293T) target genes *outside* the 5,000-gene readout panel and cannot serve as labeled regression inputs in feature-space models. Yet by biological reciprocity within protein complexes and transcriptional feedback loops, if knocking down gene $k$ strongly perturbs target gene $g$ beyond noise ($|z_{k,g}| > 4.0$, where $Z$ is the self-masked, QR-projected, column-standardized effect matrix across all reproducible knockdowns), then knocking down $g$ frequently induces a correlated downstream signature $y_k$. Weighting all reproducible knockdowns by:

$$w_{g,k}^{\text{rec}} = \mathrm{ReLU}(|z_{k,g}| - 4.0)^2, \qquad \hat{y}_g^{\text{rec}} = \frac{\sum_k w_{g,k}^{\text{rec}} \, P_Q(y_k - b_c)}{\sum_k w_{g,k}^{\text{rec}}}$$

with confidence-gated norm-matched fusion (`rec_w = 0.55`) operates in both the own cell line and inside the complementary bridge *before* translation, delivering the single largest single-stage gain (`+0.0140` across 10 splits).

### Stage 7: Output SVD Wiener Denoising & Empirical Perturbation Manifold Snapping
Finally, we learn the top-250 principal perturbation subspace $(V_{\text{out}}, S_{\text{out}})$ from all reproducible training knockdowns in line $c$, apply spectral Wiener shrinkage $s_j^2 / (s_j^2 + \tau s_1^2)$ ($\tau = 0.020$), perform top-$40$ temperature-scaled cosine + reciprocal-z attention snapping onto the empirical reproducible perturbation manifold (`snap_w = 0.45, temp = 18.0`), and apply a light final Wiener pass ($\tau = 0.015$).

---

## 3. Experimental Results

All experiments execute the exact same `PanPerturbBridgeModel` class (`src/xatlas_panhuman_solver.py`) evaluated via official `tasks.utils.strict_group_mean_from_row_scores`. Ten independent random splits are partitioned into **5 development splits** (`42, 0, 1, 2, 3`) and **5 untouched held-out confirmation splits** (`7, 11, 99, 2024, 31337`) never inspected during method design.

### 3.1 Benchmark A: 10-Split Progressive Ablation (`reports/xatlas_10split_ablation.csv`)

| Stage | Component Added | 10-Split Mean ± SD | HCT116 | HEK293T | Dev-5 | **Untouched Held-5** | Δ vs S0 (Mean ± SD) | Pos. Splits | Paired $p$-value |
|:---:|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **S0** | `100-SVD + Ridge (α=10)` (Official Ref) | 0.0613 ± 0.0040 | 0.0767 | 0.0460 | 0.0606 | 0.0621 | — | — | — |
| **S1** | `+ SelfMask + BaseSub + √n + HybridTOM` | 0.0677 ± 0.0068 | 0.0866 | 0.0487 | 0.0704 | 0.0649 | +0.0063 ± 0.0057 | 8/10 | $8.87 \times 10^{-3}$ |
| **S2** | `+ Multi-Scale Deflation Pyramid (d0..10)` | 0.0722 ± 0.0061 | 0.0935 | 0.0509 | 0.0737 | 0.0707 | +0.0109 ± 0.0046 | 10/10 | $5.69 \times 10^{-5}$ |
| **S3** | `+ Pan-Human Cross-Line Gene Graph` | 0.0748 ± 0.0053 | 0.0941 | 0.0555 | 0.0767 | 0.0729 | +0.0135 ± 0.0048 | 10/10 | $1.46 \times 10^{-5}$ |
| **S4** | `+ Empirical-SNR Non-Reproducible Anchors` | 0.0811 ± 0.0038 | 0.1033 | 0.0589 | 0.0829 | 0.0793 | +0.0197 ± 0.0040 | 10/10 | $1.14 \times 10^{-7}$ |
| **S5** | `+ GeneSlope-CosKRR Symmetric Bridge` | 0.0887 ± 0.0046 | 0.1138 | 0.0635 | 0.0906 | 0.0868 | +0.0273 ± 0.0040 | 10/10 | $7.03 \times 10^{-9}$ |
| **S6** | `+ Thresholded Reciprocal Causal Law` | 0.1027 ± 0.0065 | 0.1371 | 0.0682 | 0.1039 | 0.1014 | +0.0413 ± 0.0067 | 10/10 | $1.85 \times 10^{-8}$ |
| **S7** | **`+ SVD Wiener + Manifold Snapping (FULL)`** | **0.1099 ± 0.0087** | **0.1462** | **0.0736** | **0.1134** | **0.1064** | **+0.0485 ± 0.0085** | **10/10** | **$3.63 \times 10^{-8}$** |

### 3.2 Benchmark B: Ultra-Strict Simultaneous Union-Held-Out + Refitted Baseline & Axis (`reports/xatlas_strict_union_refit_audit.csv`)

To rule out any subtle leakage from publicly shipped `context_baseline.npz` / `residual_axis.npz` or indirect cross-line target overlap across splits, we simultaneously remove the union of HCT116 and HEK293T validation targets (`355–365` unique genes per split) from *both* cell lines and refit `context_baseline` and `residual_axis` strictly on remaining reproducible training rows:

| Model (10 Splits, Strict Union + Refitted Axis) | 10-Split Mean ± SD | Dev-5 | Untouched Held-5 | Paired Δ ± SD | Pos. Splits | Paired $t$ / $p$ |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|
| **S0** (`100-SVD + Ridge`) | 0.0658 ± 0.0033 | 0.0653 | 0.0663 | — | — | — |
| **S7** (`PanPerturb-Bridge`) | **0.1115 ± 0.0089** | **0.1132** | **0.1097** | **+0.0457 ± 0.0081** | **10/10** | $t=17.01,\; p=3.77 \times 10^{-8}$ |

---

## 4. Limitations & Strict Comparability Statement (D7)

1. **Internal Validation vs. Private Codabench Test**: All comparative claims above are made strictly on identical internal 10-split benchmarks with identical official scoring and zero target leakage. Because our internal 10-split `S0` baseline (`0.0613 ± 0.0040`) sits `~0.010` above the official `description.md` SVD reference (`0.051`), direct numerical comparison with the private Codabench leaderboard best (`0.0825`) is withheld pending live Codabench evaluation.
2. **Scope of Generalization**: `PanPerturb-Bridge` demonstrates a confirmed **Level 2 breakthrough** (`+0.0485`, `2.43×` official `0.02` noise threshold) on `xatlas-orion-perturb`; extending multi-scale causal-bridge modeling to other low-ceiling BioXArena-Pro tasks and evaluating on Codabench private test remain active next steps toward Level 3/4.
