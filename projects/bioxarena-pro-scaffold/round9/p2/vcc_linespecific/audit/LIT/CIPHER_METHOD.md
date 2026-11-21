# CIPHER v2 Method Specification: Prediction Formula, Single-Parameter Fitting, and Preprocessing

**Author**: Agent A-T5 (Literature & Official Metrics)  
**Date**: 2026-11-21  
**Primary Preprint**: Kuznets-Speck, B., Kumari, N., Senthilkumar, I., ..., Vaikuntanathan, S. & Goyal, Y. (2026). *Co-fluctuations of genes in single cells predict transcriptome-wide outcomes to perturbations* (`bioRxiv` DOI: `10.1101/2025.06.27.661814`, **v2 posted 2026-09-03**; v1 posted 2025-07-01, `PMID:40631127`).

---

## 1. Source Verification & Access Provenance

- **bioRxiv API (`api.biorxiv.org`)**: Verified v1 (2025-07-01, 9 authors) and v2 (2026-09-03, 19 authors; title updated to *"Co-fluctuations of genes in single cells predict transcriptome-wide outcomes to perturbations"*).
- **PubMed E-utilities**: Verified `PMID:40631127` (v1 record).
- **Direct bioRxiv XML/PDF access note**: As of 2026-11-21, `biorxiv.org` full-text XML (`2025.06.27.661814.source.xml`) returns HTTP 403 Cloudflare challenge (not bypassed per rules), and the preprint is not in EuropePMC OA (`hitCount: 0`).
- **Complete Official v2 Code & Companion Material Verified Directly**:
  1. Official `GoyalLab/CIPHER` GitHub repository (updated 2026-08-03, matching PyPI `cipher-perturb==0.1.0`): `cipher/{core,covariance,normalize,data,preprocess,forward}.py`, `notebooks/src/{run_fig3,run_fig4,run_fig6,suppl_forward6}.py`, `benchmarks/{README.md,parameters.csv,run_files/run_CIPHER.py}`.
  2. Official paper figures & captions in repo: `resources/paper/fig{1,2,3,4,6}.pdf`.
  3. Official companion documentation website: `https://goyallab.github.io/CIPHERWebsite/` (`CIPHER.html`, `Results.html`, `packageDocumentation.html`).
  *(Citations below reference exact paper Figure/Results sections, website pages, and source file line numbers; manuscript prose section numbers inside the Cloudflare-gated PDF are not guessed.)*

---

## 2. Theoretical Formulation & Forward Prediction Formula

### 2.1 Linear Response Relation (`CIPHERWebsite/CIPHER.html`, `Fig. 1`, `Fig. 2e`, `cipher/core.py:1–14`)
Drawing on linear response theory in statistical physics, CIPHER models the mean post-perturbation transcriptional shift $\Delta \mathbf{x} = \langle \mathbf{x} \rangle_u - \langle \mathbf{x} \rangle_0 \in \mathbb{R}^{G}$ as linear in the **unperturbed (control) cell gene–gene covariance matrix** $\Sigma \in \mathbb{R}^{G \times G}$:

$$\Delta \mathbf{x} \approx \Sigma \mathbf{u}, \qquad \Sigma_{ij} = \langle \delta x_i \, \delta x_j \rangle_0 = \frac{1}{N_{\mathrm{ctrl}} - 1} \sum_{c \in \mathcal{C}_0} (x_{ci} - \bar{x}_{0,i})(x_{cj} - \bar{x}_{0,j}),$$

where $\mathcal{C}_0$ is the unperturbed control cell population, $\bar{\mathbf{x}}_0 = \langle \mathbf{x} \rangle_0$ is the control mean vector, and $\mathbf{u} \in \mathbb{R}^G$ is a sparse perturbation driver vector.

### 2.2 Single-Gene Forward Prediction (`Fig. 3`, `cipher/core.py:10–11, 51–61`)
For a known single-gene perturbation targeting gene $g \in \{0, \dots, G-1\}$, CIPHER restricts $\mathbf{u}$ to support $\{g\}$ ($\mathbf{u} = a_g \mathbf{e}_g$), yielding a rank-1 projection onto column $g$ of $\Sigma$:

$$\widehat{\Delta \mathbf{x}}_g = \hat{a}_g \, \Sigma_{:, g}.$$

---

## 3. Single-Parameter Fitting Procedure (`cipher/core.py:33–61`, `Fig. 3f–i`, `benchmarks/run_CIPHER.py`)

### 3.1 CIPHER v2 Native Fitting (`cipher/core.py:forward_fit`, lines 33–48)
For each perturbation $g$, CIPHER fits **one scalar amplitude $\hat{a}_g$ per perturbation** (`benchmarks/README.md` §5: *100 fitted parameters for 100 perturbations; closed-form least squares*) via ordinary least squares through the origin over a gene mask $M_{\mathrm{fit}} \subseteq \{0,\dots,G-1\}$:

$$\hat{a}_g = \arg\min_{a \in \mathbb{R}} \sum_{i \in M_{\mathrm{fit}}} \left(\Delta x^{\mathrm{obs}}_{g, i} - a \, \Sigma_{i, g}\right)^2 = \frac{\sum_{i \in M_{\mathrm{fit}}} \Sigma_{i, g} \, \Delta x^{\mathrm{obs}}_{g, i}}{\sum_{i \in M_{\mathrm{fit}}} \Sigma_{i, g}^2},$$

returning `NaN` (coerced to `0.0` in `forward_predict`) if $\sum_{i \in M_{\mathrm{fit}}} \Sigma_{i, g}^2 \le 10^{-8}$ (`EPS = 1e-8`).

**How $M_{\mathrm{fit}}$ and $\Delta \mathbf{x}^{\mathrm{obs}}_g$ are constructed across CIPHER v2 experiments**:
1. **In-sample / full-transcriptome (`Fig. 3c–e`, `run_fig3.py:85–95`)**: $M_{\mathrm{fit}} = \{0, \dots, G-1\}$ (all filtered genes, target included).
2. **Gene-holdout validation (`Fig. 3f–i`, `run_fig3.py:106–130`, `cipher/core.py:gene_holdout_masks`)**: 50% of genes randomly assigned to $M_{\mathrm{fit}}$ and 50% held out for evaluation (`CONV_HOLDOUT_FRAC = 0.5`, `exclude_target_fit=True`, `exclude_target_eval=True`), verifying convergence of $\hat{a}_g$ vs. $K = |M_{\mathrm{fit}}|$ genes (`Fig. 3h–i`).
3. **External DL benchmark (`Fig. 3j–n`, `benchmarks/run_files/run_CIPHER.py:128–138`)**: Perturbed cells of gene $i$ are split 50/50 at the **cell level** (`train_idx` vs. `test_idx`), and $\hat{a}_g$ is fit on $\Delta \mathbf{x}^{\mathrm{train}}_g = \bar{\mathbf{x}}^{\mathrm{train}}_g - \bar{\mathbf{x}}^{\mathrm{train}}_{\mathrm{ctrl}}$ over all genes.
4. **Cross-dataset / CellxGene atlas covariance transfer (`Fig. 4`, `notebooks/src/run_fig4.py:91–106`)**: **Crucial distinction**: Even in `Fig. 4`, CIPHER only replaces $\Sigma$ with $\Sigma_{\mathrm{source}}$ (another study or CellxGene atlas) while **still fitting $\hat{a}_g$ on the host dataset's own observed perturbation vector $\Delta \mathbf{x}_{\mathrm{host}}$** (`forward_predict(Sigma_src, dx_host, tgi)` at `run_fig4.py:102`).

### 3.2 CRITICAL AUDIT NOTE: CIPHER v2 Native vs. Zero-Shot Cross-Line Predictor (`A-T2.md` §2)
In our **zero-shot cross-line** problem, **no perturbation cells exist for target line $t$**, so CIPHER's native $\hat{a}_g$ fit on target-line $\Delta \mathbf{x}_{g, t}$ is **infeasible**. Our A-T2 predictor (`Δ ≈ λ · C_t[:, g] / C_t[g, g]`) is a **zero-shot adaptation**:

1. **Why divide by $C_{t, g, g}$?** In linear fluctuation / conditional Gaussian theory, $\mathbb{E}[\delta x_i \mid \delta x_g] = \frac{\Sigma_{i,g}}{\Sigma_{g,g}} \, \delta x_g$. Dividing column $C_{t, :, g}$ by target gene variance $C_{t, g, g}$ normalizes the $g$-th coordinate to $1$ (representing unit knockdown/shift at target $g$, removing gene- and line-specific baseline variance scale), whereas CIPHER v2's unnormalized $\Sigma_{:, g}$ absorbs scale into the perturbation-specific $\hat{a}_g$.
2. **How $\lambda$ is fit strictly on source lines** (zero leakage from target line $t$):
   - **Global scalar $\lambda$ (1 parameter across all source perturbations & lines)**: Fit by pooled OLS over source-line pairs $(s, g)$ where $s \neq t$ (excluding target gene $g$ or including per locked config):
     $$\hat{\lambda} = \frac{\sum_{s \neq t} \sum_g \sum_{i \neq g} \left(\frac{C_{s, i, g}}{C_{s, g, g}}\right) \Delta^{\mathrm{obs}}_{s, g, i}}{\sum_{s \neq t} \sum_g \sum_{i \neq g} \left(\frac{C_{s, i, g}}{C_{s, g, g}}\right)^2}$$
     *(or fit on residual space if predicting line-specific residuals; see `A-T2.md` item 2 vs. item 3)*.
   - **Numerical stability guard**: Require $C_{t, g, g} > \epsilon$ (e.g., $\epsilon = 10^{-8}$) in target line $t$; if target gene $g$ has zero/near-zero variance in $t$'s control cells, set the CIPHER covariance prediction for $(t, g)$ to $\mathbf{0}$ (cannot propagate fluctuations from an unexpressed/zero-variance gene).

---

## 4. Data Preprocessing, Filtering, Covariance Estimation & Normalization

### 4.1 Gene and Perturbation Filtering (`cipher/data.py:load_dataset`, lines 202–300; `cipher/preprocess.py:67–88`)
- **Gene expression cutoff**: Keep genes whose raw count mean in control cells (`cutoff_source="control"`) meets `expression_threshold = 1.0` (`run_fig3.py`, `run_fig4.py`, `run_fig6.py`; `figS7` also tests `0.1`).
- **Force-include targets**: `force_include_targets=True` retains perturbed target genes even if control raw mean $< 1.0$.
- **Perturbation filtering**: Single-gene perturbations (`only_single_gene=True`, parsing labels via `strip_pert_label` / `infer_target_gene`), target present in `var_names` (`require_target_in_var=True`), minimum `min_samples = 100` perturbed cells (`50` in benchmark evaluation `run_fig3_benchmarks.py:38`).
- **Cross-dataset gene alignment (`run_fig4.py:109–137`)**: Intersect gene sets first, then apply the expression threshold on the **host** control mean (`ctrl_mean >= 1.0`) plus host targets (`run_fig4.py:112–125`) so shallow source datasets do not collapse the gene axis.

### 4.2 Control Subsampling & Covariance Estimation (`cipher/preprocess.py:368–410`, `cipher/covariance.py:17–52`)
- **Control mean**: Computed over **all** control cells (`control_raw`).
- **Control covariance $\Sigma$**: Computed on up to `cov_max_cells = 10,000` randomly subsampled control cells (without replacement, deterministic seed `stable_seed(seed, dataset_name)`), via unbiased sample covariance $\Sigma = \frac{1}{n-1}(X^\top X - n \bar{\mathbf{x}}\bar{\mathbf{x}}^\top)$ (`float64` accumulation, stored `float32`).
- **Regularization**: Forward prediction in `run_fig3.py` / `run_fig4.py` / `run_CIPHER.py` uses **unregularized** `compute_covariance` (`ridge_abs=0, ridge_rel=0, shrink=0`). Precompute (`PreprocessConfig`) adds tiny `ridge_abs = 1e-8` to the diagonal only to stabilize inverse solvers (`pinv`/`ridge`). Reverse inference optionally uses convex Ledoit-Wolf-style shrinkage `(1 - s)*Sigma + s*mean(diag)*I` with `s = 1e-3` (`cipher/covariance.py:38`).
- **Mean-field null models** (`cipher/covariance.py:67–129`, `run_fig3.py:65–81`, `Fig. 3b–g`):
  - Analytic diagonal mean-field (`Fig. 3c–e`): $\Sigma_{\mathrm{mf}} = \mathrm{diag}(\mathrm{diag}(\Sigma))$ (destroys all off-diagonal covariance).
  - Empirical column-shuffled mean-field (`Fig. 3g–i`): independently permute each gene column across control cells (`meanfield_covariance`).
  - Fully shuffled (`Fig. 4c`): permute every entry across cells and genes (`shuffled_covariance`) or `N(0,1)` noise covariance (`run_fig4.py:81–87`).

### 4.3 Six Normalization Modes & First-Order Raw Pushback (`cipher/normalize.py:28–35`, `suppl_forward6.py:850–906`, `Fig. 6j–l`, `Fig. S5`, `Fig. S16`)

| Mode | Formula ($x =$ raw count, $\ell_c = \sum_g x_{cg}$) | First-Order Raw Pushback ($\Delta x_{\mathrm{raw}} \approx s \odot \Delta y$) |
|---|---|---|
| `raw` (**paper default**) | $y = x$ | $s = \mathbf{1}$ |
| `log1p` | $y = \log(1 + x)$ | $s = 1 + \langle x_0 \rangle$ |
| `frequency` | $y = x / \ell_c$ | $s = L_0 = \langle \ell_c \rangle_0$ |
| `libsize10k` | $y = 10^4 \, x / \ell_c$ | *(not used in `suppl_forward6` pushback)* |
| `log1CP10k` | $y = \log(1 + 10^4 \, x / \ell_c)$ | $s = \langle x_0 \rangle + L_0 / 10^4$ |
| `pflog` *(recommended variance-stabilizing)* | $y_{cg} = \log(x_{cg} + \mathrm{pc}) - \frac{1}{G}\sum_{g'} \log(x_{cg'} + \mathrm{pc})$ | $s = \langle x_0 \rangle + \mathrm{pc}$ |

where `pflog` fits negative-binomial dispersion $\alpha$ across **all** un-filtered control genes ($\mu_g > 0, v_g > 0$) via OLS through the origin of $(v_g - \mu_g)$ on $\mu_g^2$, clamped to $[10^{-8}, 10^8]$, with $\mathrm{pc} = 1/(4\alpha)$ (`cipher/normalize.py:159–238`).

**Why CIPHER defaults to `raw` in single-dataset figures (`Fig. 3`, `Fig. 6l`) vs. our cross-line Addendum-B choice**:
- `Fig. 6l` (`run_fig6.py` / `Fig. 6k–l`) shows that on Replogle RPE1 raw counts, library-size normalization reduces forward Pearson vs. `raw` (`0.119 ± 0.060` vs. `0.445 ± 0.048`, Holm-adjusted $q = 0.00032$) because dividing by variable cell total UMI induces spurious negative correlations.
- **However**, in **zero-shot cross-cell-line transfer** across heterogeneous atlases (Replogle normalized bulk vs. Nadig raw scRNA-seq vs. Nourreddine/VCC), raw UMI depth differs dramatically across laboratories and cell lines, and A-T1's shared pipeline (`p2/shared/lines/`) provides Addendum-B (`log1p` library-size normalized + control z-scored) pseudobulk $\Delta$ alongside control cells (`ctrl_cells.npz`). Therefore, A-T2 must compute $C_t$ from target line $t$'s control cells in the **exact same normalization/z-score space** as $\Delta^{\mathrm{obs}}$ so that $C_{t,:,g} / C_{t,g,g}$ and $\Delta^{\mathrm{obs}}$ share identical units and scale (while preserving raw/non-z-scored variants only if explicitly tested in A-T1/A-T2).

---

## 5. CIPHER v2's `μ`-Axis vs. Orthogonal Decomposition (`Fig. 6c–f`, `run_fig6.py:244–326`)

In `Fig. 6c–f`, CIPHER v2 decomposes a dataset's perturbation responses along the **within-line control mean vector** $\hat{\boldsymbol{\mu}} = \bar{\mathbf{x}}_0 / \|\bar{\mathbf{x}}_0\|_2$ (interpreted as global total-RNA / growth scaling $c = 1 - \Delta x / \mu$):
- `dx_perp = dx - <dx, mu_hat> * mu_hat`, refitting CIPHER (`forward_fit(Sigma[:, g], dx_perp)`) on the $\mu$-orthogonal residual.
- **Key Distinction from Team A's Question**: CIPHER v2's `Fig. 6` decomposition removes the **within-line control mean axis $\boldsymbol{\mu}$** while still fitting $\hat{a}_g$ on target-line $\Delta \mathbf{x}$. CIPHER v2 **never** tests whether target-line control covariance $C_t$ recovers **zero-shot cross-line residuals** ($r_{p,t} = \Delta^{\mathrm{obs}}_{p,t} - \bar{\Delta}^{\mathrm{src}}_{p,(-t)}$) after removing the cross-line shared perturbation response and/or training PC1!
