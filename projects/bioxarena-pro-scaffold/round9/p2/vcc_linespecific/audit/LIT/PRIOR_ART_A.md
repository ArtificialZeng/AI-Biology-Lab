# Prior Art Audit: VCC 2026 Top-3, CASTER, CGC, HyperMap & Co-Fluctuation Transfer

**Author**: Agent A-T5 (Literature & Official Metrics)  
**Date**: 2026-11-21  
**Core Question**: *Has any published or public work already used target-line unperturbed gene–gene co-fluctuations ($C_t$) in zero-shot mode to recover the cell-line-specific residual component missed by cross-line average baselines?*

---

## 1. Executive Verdict

**Within our verified literature and code audit (through 2026-11-21), NO prior work has tested target-line gene–gene co-fluctuation linear response ($C_{t,:,g}/C_{t,g,g}$ or CIPHER-style covariance propagation) in zero-shot transfer under a cross-line shared + line-specific residual decomposition.**

Key nuanced findings across the 6 primary targets + 3 adjacent works:
1. **BioInsight/Illumina `CellularAI` (VCC 2026 #1, `0.4399`)** uses **mean/summed** target control logCPM in its `ControlResidual` model (`src/cellularai/prediction/residual.py`), and uses top-$k$ target control PCs (`control_pcs` in `totals.py`) **only for post-hoc smoothing/reweighting** (`D' = D + lam * V_t^T (V_t D)`) of the pooled consensus profile—not gene-specific CIPHER linear response. Empirically, VCC 2026 Wrap-Up showed top-10 entries (including BioInsight) still only reached line-specific slope **`0.02–0.13`** (vs. `1.02` replicate).
2. **`FlowState` (VCC 2026 #2, `0.3348`)** and **Xaira `X-team` / `X-Cell Transfer` (VCC 2026 #3, `0.3263`)** condition neural predictors or reference-pooling weights on target control cells/profiles without explicit covariance propagation $C_t$, and both likewise collapsed toward cross-line average responses (`0.77–0.89` cross-line cosine, `0.02–0.13` slope).
3. **`CASTER` (Kawsar et al., 2026-11-17)** uses a Set Transformer over $S=32$ target control cells + `PerturbMean` prior (`CanPert/src/canpert/model.py`) with **no gene–gene covariance matrix**. In zero-shot transfer (`0%` screened), CASTER **failed pre-registered H2** and scored **below `PerturbMean`**.
4. **`CGC` (Huang et al., 2026-10-01)** diagnoses *Context Geometry Compression* across contexts without proposing a co-fluctuation predictor.
5. **`HyperMap` (Dhaka et al., 2026-04-27)** is **not zero-shot** (`n_adapt_genes = 10/20` measured target perturbation seeds + inner-loop SGD in `hypermap/trainer.py`), operates on single control cell vectors without gene–gene covariance $C_t$, and does not evaluate cross-line shared + line-specific residuals.
6. **`CIPHER v2` (Kuznets-Speck et al., 2026-09-03)** uses control covariance $\Sigma$, **but always fits its per-perturbation scalar $\hat{a}_g$ on target-dataset perturbation responses $\Delta x$** (even in `Fig. 4` cross-dataset transfer `run_fig4.py:102`) and was never tested zero-shot on unseen cell lines or on cross-line residual slopes.
7. **Important Skeptical Prior (`Molina & Zhang 2026`, `CASTER 2026`, `CGC 2026`)**: Three 2026 studies explicitly report failure to recover gene-level line-specific components from target control expression alone (via first-moment / neural set representations). Whether **second-order gene–gene co-fluctuation structure $C_t$** breaks this `0.13` slope ceiling—or confirms the negative bound under strict zero-shot preregistration—is the exact falsifiable gap Team A tests.

---

## 2. Detailed Method Audit Table

| Work / Model | Source Verified | How Target Control Cells Are Used | Uses Target Gene–Gene Co-Fluctuations $C_t$? | Zero-Shot (0 Target Perts)? | Isolates Line-Specific Residual? |
|---|---|---|---|---|---|
| **BioInsight / Illumina `CellularAI`** (VCC 2026 #1, `0.4399`) | Wrap-Up (2026-11-19); `Illumina/CellularAI` repo (`method.md`, `residual.py`, `totals.py`, `weights.py`, `configs/example.yaml`) | (1) `sim_factor`: Spearman corr of target control log1p(CPM) (2k HVGs) for random-effects meta-analysis; (2) `ControlResidual`: **summed** control log1p(CPM) + 5 global PCA scores ridge-regressed onto low-rank residuals (`preserve_signs`, `[0.5, 2.0]` cap); (3) `totals.py` `control_pcs`: top-100 control cell PCs for post-hoc `CP`/`PW` profile smoothing & jackknife spread matching | **Partial (post-hoc PC smoothing only, NOT linear response $C_{t,:,g}$)**: `ControlResidual` uses 1st-moment summed controls; `cp_design` smooths already-pooled effects along control PCs (`D + lam * V^T V D`), not targeted gene-$g$ covariance column propagation | **Yes** | **No** (still `0.02–0.13` slope in Wrap-Up) |
| **`FlowState`** (VCC 2026 #2, `0.3348`) | Wrap-Up (2026-11-19) | Untreated control cells of new line passed as conditioning prompt to 3-head effect model (rank, direction, magnitude), trained on 222 public datasets; cell generator adjusts control counts | **No** (neural context prompt; no explicit covariance propagation documented) | **Yes** | **No** (`JACCARD -0.0018`; top-10 slope `0.02–0.13`) |
| **Xaira `X-team` / `X-Cell Transfer`** (VCC 2026 #3, `0.3263`) | Wrap-Up (2026-11-19); bioRxiv `10.64898/2026.03.18.712807` (2026-03-20); `Xaira-Therapeutics/X-Cell` docs | `X-Cell` masked diffusion over cell sets with Flamingo cross-attention to 6 priors; VCC `X-Cell Transfer` pools reference screens by **control similarity**, gates with `X-Cell` agreement, and redistributes counts over 400 control cells | **No** (set-level masked diffusion + reference control similarity; no $C_t$ linear response) | **Yes** | **No** (top-10 slope `0.02–0.13`) |
| **Kawsar et al. `CASTER`** (`10.64898/2026.11.15.766815`, 2026-11-17) | bioRxiv API; Zenodo `24239385` (`CanPert/src/canpert/{model,train,results_text}.py`) | Set Transformer (`S=32` control cells) -> cell MLP `f_cell`, mean-pooled `ctx` -> FiLM `(gamma, beta)` + gene scaling `g_scale(ctx)` modulating `PerturbMean` prior | **No** (attention is across $S=32$ cells in embedding space, no $G \times G$ gene covariance $C_t$) | Tested **both** 0% & partial (`2.5–30%`) | **Failed at 0%**: fell **below `PerturbMean`** (`H2` not supported) |
| **Huang, Wang, Li & Wilson `CGC`** (`10.64898/2026.09.30.755756`, 2026-10-01) | bioRxiv API | Diagnostic benchmark of *Context Geometry Compression* across cell contexts & drug matrices | **No** (evaluates recovery at pathway vs. gene resolution & compact baseline representations) | Diagnostic | **Diagnoses collapse** of context-specific effects |
| **Dhaka, Gao & Ideker `HyperMap`** (`10.64898/2026.04.23.720505`, 2026-04-27) | bioRxiv API; `bhavya1929/HyperMap` (`hypermap/{model,trainer}.py`) | Single control cell vectors `control_exp` fed to `CellEncoder`; hypernetwork generates predictor weights, adapted at test time via inner-loop SGD | **No** (single-cell MLP encoder, no gene–gene covariance $C_t$) | **No** (requires `n_adapt_genes = 10–20` target seeds) | **No** (no shared vs. residual decomposition) |
| **Kuznets-Speck et al. `CIPHER v2`** (`10.1101/2025.06.27.661814` v2, 2026-09-03) | bioRxiv API; `GoyalLab/CIPHER` (`core.py`, `run_fig4.py`, `run_fig6.py`) | Computes control covariance $\Sigma$; predicts $\hat{\Delta}_g = \hat{a}_g \Sigma_{:,g}$ | **Yes**, on same dataset (`Fig. 3`) or source $\Sigma$ (`Fig. 4`) | **No** ($\hat{a}_g$ always fit on target $\Delta x$) | **No** (`Fig. 6` projects out control mean $\mu$, not cross-line shared response) |
| **Molina & Zhang** (`10.64898/2026.07.24.740459`, 2026-07-27) | bioRxiv API; `xinyizhanglab/perturbation-decomposition` (`anova.py`, `controls.py`) | Exports control `mean` and `controls.npz` for model inference; tests whether control expression predicts ANOVA components $(\mu, \alpha_c, \beta_p, \gamma_{c,p})$ | **No** (code audit confirms no gene–gene covariance $C_t$ predictor) | **Yes** (unseen cell lines via priors) | **Yes (ANOVA)** — finds control expression predicts $\mu$/global only, **not** line-specific $\gamma_{c,p}$ |

---

## 3. Per-Work Technical Notes (Verified Against Source Text & Code)

### 3.1 BioInsight AI Lab / Illumina `CellularAI` (VCC 2026 1st Place, Score `0.4399`)
- **Components & Scores** (Wrap-Up, 2026-11-19): Overall `0.4399` (`PDS 0.8313`, `MSE 0.6842`, `NMAE 0.4399`, `FIDELITY 0.1908`, `REACH 0.4288`, `JACCARD 0.0642`); called ~500 DEGs vs. real `22–44`. Trained on proprietary **Illumina Billion Cell Atlas** + Nourreddine 2026 + Feng 2026 + Zhu 2026.
- **Exact Control Usage Verified in `Illumina/CellularAI` Source Code**:
  1. **Meta-analysis pooling (`pooling/weights.py:62–72`)**: Computes Spearman correlation between target control `log1p(CPM)` and each training dataset's control `log1p(CPM)` over top 2,000 HVGs (`sim_hvg: 2000`), weighting random-effects meta-analysis (`sim: 0.25`). Uses **mean profile only**.
  2. **Residual model (`prediction/residual.py:64–192`)**: Fits a low-rank SVD basis `U` (`residual_rank: 160`) across training-line residuals (`own-target` masked by off-diagonal row mean). Predicts residual coordinates via per-gene ridge regression (`residual_ridge: 10.0`) on standardized **summed** control `log1p(CPM)` + 5 global control PCA scores (`residual_global: 5`). Applied with shrinkage `residual_alpha: 0.25`, `bounded_ratio` (`[0.5, 2.0]`), and **`preserve_signs`** (cannot reverse direction!). Own targets are left untouched. Uses **first-moment summed controls**, NOT cell–cell covariance.
  3. **Control-program smoothing & spread matching (`prediction/totals.py:123–190`, `prediction/build.py:47`)**: `control_pcs()` extracts top `pw_k = 100` right singular vectors $V_t$ (`sv**2`) of `centred log1p(TS * x / lib)` over `pcs_cells = 6000` target control cells. `cp_design()` and `pw_design()` then smooth/reweight the already-pooled + residual-adjusted log2 fold-change profile $D$ along $V_t$ (`D' = D + lam * g_p * V_t^T (V_t D)` for top `33%` / `25%` perturbations by $|z| \ge 5.0$), rescaling amplitude to match delete-one jackknife variation credit (`amp_share = 0.9`).
- **Implication for Team A**: Even `CellularAI`'s `cp_design`/`pw_design` only smooths an existing consensus effect vector $D$ along dominant control PCs—it never propagates a perturbation from target gene $g$ via $C_{t,:,g}/C_{t,g,g}$, and explicitly forbids sign flips in `ControlResidual`. Despite proprietary billion-cell training data and `control_pcs` smoothing, top-10 models including BioInsight still recovered only `0.02–0.13` line-specific slope.

### 3.2 `FlowState` (VCC 2026 2nd Place, Score `0.3348`)
- **Components & Scores** (Wrap-Up, 2026-11-19): Overall `0.3348` (`PDS 0.8363`, `MSE 0.3632`, `NMAE 0.3286`, `FIDELITY 0.1695`, `REACH 0.3128`, `JACCARD -0.0018`). Built by independent researchers Sina Khodaveisi & Navid Bagheri using exclusively 222 public perturbation datasets (128 studies, 86 contexts: single-cell + bulk ENCODE/LINCS/PerturbAtlas/KnockTF/ARCHS4 across CRISPRi/KO/siRNA/shRNA).
- **Control Usage**: Feeds target untreated control cells as a prompt into a 3-head effect model (response rank, direction, magnitude), trained by holding out contexts; cell generator adjusts supplied control counts. No public codebase or evidence of target control covariance $C_t$ propagation.

### 3.3 Xaira `X-team` / `X-Cell Transfer` (VCC 2026 3rd Place, Score `0.3263`)
- **Components & Scores** (Wrap-Up, 2026-11-19): Overall `0.3263` (`PDS 0.7346`, `MSE 0.4681`, `NMAE 0.3038`, `FIDELITY 0.0659`, `REACH 0.3186`, `JACCARD 0.0667`); called `122–136` DEGs vs. real `22–44`.
- **Architecture & Control Usage** (`10.64898/2026.03.18.712807`, `X-Cell/model/`): `X-Cell` is a set-level masked diffusion transformer (55M Mini to 4.9B Ultra) trained on 25.6M–37M+ cells including `X-Atlas/Pisces` (excluded by Rule 5), with Flamingo cross-attention to ESM-2, STRING, GenePT, DepMap, and JUMP-CP priors. For VCC 2026, `X-Cell Transfer` anchored predictions to reference CRISPRi screens weighted by cell count and **control-cell similarity** to the target line, gated by `X-Cell` agreement/confidence, plus a GBDT DEG classifier and count redistribution across 400 sampled control cells. No gene–gene co-fluctuation linear response $C_t$ is used.

### 3.4 Kawsar et al. `CASTER` (`10.64898/2026.11.15.766815`, 2026-11-17)
- **Evaluated lines**: K562, RPE1 (Replogle 2022) and HepG2, Jurkat (Nadig 2025), using 2,000 control-selected genes, excluding target genes (`CanPert/src/canpert/results_text.py:142–157`).
- **Architecture (`model.py:60–111`)**: Encodes $S=32$ control cells via 4-layer MLP `f_cell`, mean-pools context vector `ctx = h_cell.mean(dim=1)`, applies FiLM conditioning to `PerturbMean` transfer prior `d_src`, runs cell-set Transformer `f_st`, and predicts $\hat{X} = X_{\mathrm{ctrl}} + (1 + g(c)) \odot d_{\mathrm{src}} + f_{\mathrm{out}}(O)$. No gene–gene covariance matrix $C_t$ is computed.
- **Key Negative Zero-Shot Finding (`results_text.py:377–406`, `abstract`)**:
  - With `2.5%`–`30%` of target perturbations screened (`H1`), CASTER beats `PerturbMean` (`0.62` vs. `0.47` noise-corrected at `2.5%`; `0.80` vs. `0.44` at `30%`, behind a `0.84` linear model).
  - With **`0%` target perturbations (zero-shot, `H2`)**, CASTER fell **significantly below `PerturbMean`** on all three primary metrics (discrimination dropped below `PerturbMean`, pulling toward generic response). Note CASTER also found that predicted DEG-size differentials tracked DepMap Chronos line-specific dependency in K562 (`H4b`, `results_text.py:280–282`), though zero-shot control-only adaptation failed.

### 3.5 Huang, Wang, Li & Wilson `CGC` (`10.64898/2026.09.30.755756`, 2026-10-01)
- Identifies **Context Geometry Compression (CGC)**: models preserve basal cell identity and generic shared perturbation vectors while compressing differences across contexts (`gamma`). Demonstrates that even in large drug-response matrices with near-complete perturbation coverage, gene-level context-specific differences remain difficult to recover relative to coarse pathway summaries. Does not propose or test target control covariance $C_t$ for recovering line-specific residuals.

### 3.6 Dhaka, Gao & Ideker `HyperMap` (`10.64898/2026.04.23.720505`, 2026-04-27)
- Meta-learning hypernetwork (`bhavya1929/HyperMap/hypermap/model.py`, `trainer.py`): `CellEncoder` + `PerturbationEncoder` + `ContextFusion` + `HyperNetwork` generating predictor weights from `(cell_emb || context_vector)`.
- **Why distinct from Team A**: (1) Requires **target-line perturbation "seeds"** (`n_adapt_genes = 10` or `20`) and `higher` inner-loop SGD steps (`n_adapt_steps = 1..3`) at test time—not zero-shot; (2) encodes individual control cell vectors `control_exp` after subtracting control mean (`_preprocess_data`), never computing gene–gene covariance $C_t$; (3) does not isolate or score the cross-line shared vs. line-specific residual slope.

---

## 4. Verification Sources
- Arc VCC 2026 Wrap-Up (2026-11-19): `https://arcinstitute.org/news/virtual-cell-challenge-2026-wrap-up`
- `Illumina/CellularAI`: `https://github.com/Illumina/CellularAI` (`README.md`, `docs/method.md`, `configs/example.yaml`, `src/cellularai/`)
- `Xaira-Therapeutics/X-Cell`: `https://xaira-therapeutics.github.io/X-Cell/`, `https://github.com/Xaira-Therapeutics/X-Cell`, bioRxiv `10.64898/2026.03.18.712807` (verified via `api.biorxiv.org`)
- Kawsar et al. `CASTER`: bioRxiv `10.64898/2026.11.15.766815` (verified via `api.biorxiv.org`), Zenodo `10.5281/zenodo.24239385` (`CanPert.zip`)
- Huang et al. `CGC`: bioRxiv `10.64898/2026.09.30.755756` (verified via `api.biorxiv.org`)
- Dhaka et al. `HyperMap`: bioRxiv `10.64898/2026.04.23.720505` (verified via `api.biorxiv.org`), `https://github.com/bhavya1929/HyperMap`
- Molina & Zhang: bioRxiv `10.64898/2026.07.24.740459` (verified via `api.biorxiv.org`), `https://github.com/xinyizhanglab/perturbation-decomposition`
- Kuznets-Speck et al. `CIPHER v2`: bioRxiv `10.1101/2025.06.27.661814` v2 (verified via `api.biorxiv.org`), `https://github.com/GoyalLab/CIPHER`
- Verma et al. `PIE`: bioRxiv `10.64898/2026.10.02.756297` (verified via `api.biorxiv.org`)
