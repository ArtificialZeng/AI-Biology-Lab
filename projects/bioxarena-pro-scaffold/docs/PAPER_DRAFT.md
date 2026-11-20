# Residual-Space Cross-Line Perturbation Modeling for Unseen CRISPRi Knockdowns: Metric Trade-Offs, Leave-One-Out Transferability, and External Validation

**Working Manuscript Draft (Round 6 Complete — All T0–T7 Integrated)**  
**Primary Benchmark**: BioXArena-Pro `perturbation/xatlas-orion-perturb` (HCT116, HEK293T)  
**External Independent Dataset**: Replogle et al. (*Cell* 2022) Perturb-seq (K562 GWPS, RPE1 Essential; Figshare+ 20029387)  
**Official Codabench Private Test (Single Frozen Submission `1062221`)**: `0.0998` (leading within the `0.02` noise band; human-AI multi-round iteration, NOT comparable to single-shot 2h autonomous agents)

---

## Abstract

Evaluating transcriptional response predictors for unseen CRISPRi gene knockdowns exposes a fundamental metric trade-off between **unprojected global correlation** (which retains dominant shared stress variation across perturbations) and **residualized or retrieval-based specificity** (which evaluates whether a model distinguishes one gene knockdown from another). Under standard unprojected $\Delta$-expression Pearson correlation, simple low-rank linear baselines outperform many deep models (Ahlmann-Eltze et al., *Nature Methods* 2025). Here we systematically evaluate a frozen, cell-line-agnostic residual-space kernel and cross-line bridging pipeline (`Cross-Line Perturb-Bridge`, `PanPerturbBridgeModel`) across two independent human Perturb-seq atlases—BioXArena-Pro `xatlas-orion-perturb` (HCT116, HEK293T) and Replogle et al. (*Cell* 2022; K562, RPE1)—using **per-split row-level stratified paired bootstrap 95% confidence intervals (4,000 resamples)** and **five complementary metrics** (v5 orthogonal-residual Pearson, unprojected $\Delta$-expression Pearson, Top-50 DE Pearson, Top-50 sign accuracy, and v5-independent raw-space perturbation discrimination rank percentile `disc_raw`).

Our findings establish four concrete, boundary-aware conclusions:
1. **Metric trade-off across both atlases**: On `xatlas-orion-perturb` (10 development splits), full model S7 improves v5 residual Pearson over `100-SVD + Ridge` (`0.0613 → 0.1099`, median 95% CI `[+0.0318, +0.0677]`, 10/10 CI > 0) and the Ahlmann-Eltze bilinear baseline (`0.0672 → 0.1099`, `[+0.0256, +0.0606]`, 10/10 CI > 0); on a single official Codabench private test submission (`Submission ID 1062221`), S7 achieved **`0.0998`** (measured validation-to-test shift `−0.0101`), `+0.0173` above the best 2-hour autonomous agent run (`0.0825`), which lies within the official `0.02` noise threshold and represents a different development category (human-AI multi-round iteration with ~25 tuned constants). Evaluated **zero-shot with frozen hyperparameters** on Replogle 2022 (10 strict simultaneous union-held-out splits), co-primary variants **S5** (`0.1401`) and **S7** (`0.1359`) both exceed linear baselines on v5 (`AE 0.0967`: S5 `[+0.0323, +0.0549]`, 10/10 CI > 0; S7 `[+0.0240, +0.0537]`, 10/10 CI > 0), and **S5** also exceeds linear baselines on **v5-independent raw-space perturbation discrimination (`disc_raw`)** (`0.6884` vs `AE 0.6468`, `+0.0417`, `[+0.0198, +0.0721]`, 9/10 CI > 0; S7 `+0.0268`, 4/10 CI > 0). Conversely, against a self-implemented untuned simplified GNN (`0.1233`), **S7** CI excludes 0 on only **3/10 splits** (`[-0.0041, +0.0286]`), while **S5** excludes 0 on **7/10 splits** (`[+0.0038, +0.0302]`). Critically, on **unprojected $\Delta$-expression Pearson, Top-50 DE Pearson, and Top-50 sign accuracy**, both S5 and S7 **lose to linear baselines** (`Replogle delta_r S7−AE = −0.0308`, 8/10 CI < 0; `S5−AE = −0.0277`, 8/10 CI < 0) because orthogonal residual projection deliberately forfeits the shared stress axis.
2. **Leave-one-module-out (LOMO) transferability**: Across both atlases, only **four modules** are consistently positive in LOMO: **soft-power `HybridTOM(p=1.5)`** (`+0.0109` xatlas / `+0.0103` Replogle, 10/10 CI > 0 on Replogle), **output SVD Wiener + manifold snapping** (`+0.0072` / `+0.0067` on v5, but reduces `disc_raw` on both datasets), **cross-line bridging (`S5`)** (`+0.0020` / `+0.0061`, 9/10 CI > 0 on Replogle), and **multi-scale spectral deflation (`S2`)** (`+0.0108` / `+0.0035`, Replogle 1/10 CI > 0). In contrast, **Reciprocal Response Borrowing (`S6`)** is strongly dataset-specific (`+0.0130`, 8/10 CI > 0 on xatlas, but **`−0.0153`, 9/10 CI < 0** on Replogle; omitting `S6` raises Replogle v5 to `0.1512`), and non-reproducible anchors (`S4`) benefit xatlas only (`+0.0043` vs `−0.0006`).
3. **Evaluated-background pathway enrichment & CRISPRi cis-proximity audit**: Local `ReactomePathways.gmt` hypergeometric testing strictly within the 1,784 evaluated Replogle genes confirms that top-predicted sets of S7 and S0 overlap heavily (`113/164`, `68.9%`), whereas gain-ranked cohorts enrich **mRNA Splicing** (`S7−S0` Top-150: `3.12×`, `FDR = 1.18e-08`) and **Mitochondrial translation** (`S5−S2` Top-150: `2.90×`, `FDR = 2.76e-03`). On Replogle GRCh38 TSS coordinates (GENCODE v46), negative-z (`z < −4`) reciprocal donor pairs show a **474×** enrichment for `< 10 kb` genomic proximity (`1.422%` vs `0.003%` null, `3.32%` weight), though removing `< 10 kb` donors retains `97.41%` of Replogle's top-40 positive gains. Because xatlas gene identities are anonymized and cannot be inspected under benchmark rules, CRISPRi local cis-repression remains an **unexcluded alternative explanation** for xatlas's S6 gain (where 78% of top-40 gain genes have strong negative-z donors).

---

## 1. Main Quantitative Results: Residual Metrics vs. Traditional Metrics (T2 & T3)

### 1.1 Five-Metric External Benchmark on Replogle et al. (*Cell* 2022, 10 Strict Union Splits, Zero Retuning)

**Preprocessing Disclosure (`src/prepare_replogle2022_cache.py`)**: Top 5,000 finite readout genes per line were selected via `var/clean_mean` (`4,991/5,000` RPE1 and `4,996/5,000` K562 overlap with `var/mean`); `reproducible = 1` iff `ad_counts >= 50` AND `energy_test_p_value < 0.001`. When multiple guides targeted the same Ensembl gene, a single row was selected maximizing `(reproducible, ad_counts, num_cells_filtered)`, which uses perturbation-level summary statistics across the screen.

| Method (Frozen / Untuned; **S5 & S7 Co-Primary**) | v5 Residual $r$ | Δ-Expr $r$ | Top-50 DE $r$ | Top-50 Sign | `Discrim(Delta)` (raw) | `Discrim(v5)` | Row-Level Median 95% CI (v5 / `disc_raw`) |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---|
| Constant Mean | 0.0009 | −0.0322 | −0.0531 | 0.4992 | 0.5000 | 0.4999 | — |
| `30-kNN` | 0.0971 | **0.1479** | 0.1822 | 0.6095 | 0.6365 | 0.6329 | — |
| Ahlmann-Eltze Bilinear (AE) | 0.0967 | 0.1470 | 0.2026 | 0.6255 | 0.6461 | 0.6375 | — |
| `100-SVD + Ridge (S0)` | 0.0883 | 0.1408 | 0.2119 | **0.6256** | 0.6505 | 0.6404 | — |
| Self-implemented untuned simplified ResMLP | 0.1091 | 0.0876 | 0.1896 | 0.5901 | 0.6447 | 0.6529 | — |
| Self-implemented untuned simplified GNN | 0.1233 | 0.0971 | 0.1892 | 0.5999 | 0.6403 | 0.6499 | — |
| **S5** (`+ Cross-Line Bridge`, **external peak**) | **0.1401** | 0.1193 | **0.2170** | 0.6064 | 0.6890 | 0.6893 | v5 vs AE `[+0.0323, +0.0549]` (10/10); vs GNN `[+0.0038, +0.0302]` (7/10); `disc_raw` vs AE `[+0.0198, +0.0721]` (9/10) |
| S6 (`+ Reciprocal Borrowing`, **v5 negative**) | 0.1292 | 0.1106 | 0.2155 | 0.6020 | **0.6993** | **0.7001** | v5 vs S5 `[-0.0169, -0.0050]` (9/10 CI < 0) |
| **S7** (`+ Wiener + Manifold Snap`, full frozen) | 0.1359 | 0.1162 | 0.2013 | 0.5962 | 0.6762 | 0.6771 | v5 vs AE `[+0.0240, +0.0537]` (10/10); vs GNN `[-0.0041, +0.0286]` (3/10); `disc_raw` vs AE `[-0.0011, +0.0546]` (4/10) |

### 1.2 Five-Metric Benchmark & Official Private Test on `xatlas-orion-perturb`

| Method / Protocol | v5 (`dev10`) | Δ-Expr $r$ | Top-50 DE $r$ | Top-50 Sign | `Discrim(Delta)` (raw) | `Discrim(v5)` | `resplit10` v5 / `Discrim(Delta)` |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| Constant Mean | 0.0090 | 0.0608 | 0.0941 | 0.5420 | 0.5000 | 0.4999 | 0.0073 / 0.5000 |
| `30-kNN` | 0.0410 | 0.0832 | 0.1327 | 0.5669 | 0.5698 | 0.5664 | 0.0391 / 0.5639 |
| AE Bilinear Baseline | 0.0672 | **0.1367** | 0.2039 | **0.6048** | 0.5856 | 0.5952 | 0.0646 / 0.5829 |
| `100-SVD + Ridge (S0)` | 0.0613 | 0.1327 | **0.2050** | 0.6046 | 0.5891 | 0.5993 | 0.0596 / 0.5870 |
| Self-implemented untuned ResMLP / GNN | `0.0486 / 0.0543` | `0.0390 / 0.0432` | `0.0964 / 0.1025` | `0.5425 / 0.5464` | `0.5951 / 0.5900` | `0.5962 / 0.5922` | — |
| **S5** (`+ Cross-Line Bridge`) | 0.0887 | 0.0813 | 0.1576 | 0.5768 | 0.6262 | 0.6262 | 0.0878 / 0.6237 |
| S6 (`+ Reciprocal Borrowing`) | 0.1027 | 0.0937 | 0.1885 | 0.5875 | **0.6453** | **0.6471** | 0.1006 / **0.6439** |
| **S7** (Full Frozen) | **0.1099** | 0.0995 | 0.1789 | 0.5857 | 0.6286 | 0.6313 | **0.1073** / 0.6263 |
| **Codabench Official Private Test (`1062221`)** | **`0.0998`** | — | — | — | — | — | `+0.0173` vs best 2h agent `0.0825` (within `0.02` noise; human-AI multi-round) |

**Why Residual Methods Lose on Traditional Metrics (`DeltaR`, `Top-50 DE`, `Sign`) While Gaining on Perturbation Discrimination**: Shared systemic stress variation across CRISPRi perturbations accounts for a large fraction of unprojected variance and disproportionately populates each knockdown's Top-50 differentially expressed genes (as evidenced by `Constant Mean` alone reaching `DE50 = 0.0941` and `Sign = 0.5420` on xatlas). Because `S1–S7` explicitly project training targets orthogonal to `residual_axis` ($P_Q$), they discard the shared stress component, systematically lowering unprojected `DeltaR`, `Top-50 DE`, and `Sign Accuracy` relative to linear models that predict along the stress axis. In contrast, cross-perturbation retrieval (`disc_raw`) cannot exploit a constant stress direction shared across perturbations, revealing a genuine specificity advantage for `S5` (`+0.0417`, 9/10 CI > 0 on Replogle) and exposing `S7` manifold snapping's bias-variance trade-off (higher average v5 cosine at the expense of reduced individual distinctiveness).

---

## 2. Leave-One-Module-Out (LOMO) Ablation Across Both Atlases (T4)

| Module Left Out from Full `S7` | xatlas Contrib (`Full − LOO`, med CI) | Replogle Contrib (`Full − LOO`, med CI) | Classification |
|:---|:---:|:---:|:---|
| **`S1_HybridTOM` ($p=1.5$)** | **`+0.0109`** (`[+0.0027, +0.0193]`, 7/10 CI>0) | **`+0.0103`** (`[+0.0038, +0.0155]`, **10/10 CI>0**) | **Transferable core contribution** |
| **`S7_WienerSnap`** | **`+0.0072`** (`[+0.0008, +0.0142]`, 6/10 CI>0) | **`+0.0067`** (`[+0.0006, +0.0114]`, **8/10 CI>0**) | **Transferable for v5** (hurts `disc_raw`) |
| **`S5_CrossLineBridge`** | `+0.0020` (`[-0.0058, +0.0098]`, 8/10 pos) | **`+0.0061`** (`[+0.0022, +0.0103]`, **9/10 CI>0**) | **Transferable core contribution** |
| **`S2_Deflation` (`d0..d10`)** | **`+0.0108`** (`[+0.0004, +0.0198]`, 5/10 CI>0) | `+0.0035` (`[-0.0051, +0.0113]`, 9/10 pos, 1/10 CI>0) | **Transferable (weaker on Replogle)** |
| `S1_SqrtNCells` | `+0.0011` (`[-0.0020, +0.0050]`) | **`+0.0068`** (`[+0.0009, +0.0123]`, 6/10 CI>0) | Replogle-strong, xatlas-mild |
| `S1_SelfMask` / `S1_BaseSub` / `S3_CrossLineGraph` | `+0.0004` / `−0.0007` / `−0.0002` | `+0.0006` / `+0.0007` / `+0.0016` | Minor / redundant in full system |
| `S4_NonRepAnchors` | **`+0.0043`** (`[+0.0006, +0.0073]`, 6/10 CI>0) | `−0.0006` (`[-0.0018, +0.0008]`, 2/10 CI<0) | **Dataset-specific (xatlas only)** |
| **`S6_Reciprocal`** | **`+0.0130`** (`[+0.0016, +0.0253]`, **8/10 CI>0**) | **`−0.0153`** (`[-0.0240, -0.0080]`, **9/10 CI<0**) | **Dataset-specific / externally negative** (LOO `v5 = 0.1512`) |

---

## 3. Pathway Enrichment (T5) & CRISPRi Proximity Confounding Audit (T6)

- **T5 Local `ReactomePathways.gmt` (`M_bg = 1,784` evaluated unique Replogle genes)**: Top-accuracy genes of `S7` and `S0` overlap `113/164` (`68.9%`) and enrich the same housekeeping pathways. Gain-ranked sets against the evaluated background yield **mRNA Splicing - Major Pathway** (`S7−S0` Top-150: `38/145`, `3.12×`, `FDR = 1.18e-08`; `1,452` GMT-annotated background `FDR = 1.7e-08`; random 150-gene control `0/50`) and **Mitochondrial translation** (`S5−S2` Top-150: `20/82`, `2.90×`, `FDR = 2.76e-03`).
- **T6 CRISPRi Proximity Audit (Replogle GRCh38 GENCODE v46, `10,530/10,539` mapped)**: Negative-z donor pairs (`z < −4`, `n = 3,939`) exhibit a **474×** enrichment for `< 10 kb` TSS proximity (`1.422%` vs `0.003%` null, `3.32%` weight; `< 1 kb` `1.371%`, 685×). On Replogle, removing `< 10 kb` donors retains `97.41%` of top-40 ownRec positive gains (`0.2536 → 0.2471`), showing that selected Replogle qualitative cases (`DICER1–DGCR8/DROSHA`, `NUDT21–CPSF6`, `MYC–TELO2/TTI1/TRRAP`) are trans rather than cis. However, on anonymized `xatlas`, 78% of top-40 `S6`-gain genes have strong negative-z top donors, and because benchmark rules forbid inspecting xatlas gene identities or coordinates, **CRISPRi local cis-repression cannot be excluded as a contributor to xatlas's S6 gain**.
