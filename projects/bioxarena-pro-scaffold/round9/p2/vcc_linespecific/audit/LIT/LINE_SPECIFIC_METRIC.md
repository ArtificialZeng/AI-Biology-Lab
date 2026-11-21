# VCC 2026 Line-Specific Slope & Official `cell-eval2` Metrics Specification

**Author**: Agent A-T5 (Literature & Official Metrics)  
**Date**: 2026-11-21  
**Target Consumers**: A-T2 (`src/metrics.py`), A-T3 (`prereg/PREREG_A.md` & `audit/`), A-T4 (`prospective/`)

---

## 1. Primary Source Audit & Epistemic Boundary (What Is Official vs. Underspecified)

### 1.1 Verified Official Sources
1. **Arc VCC 2026 Winners and Wrap-Up** (`https://arcinstitute.org/news/virtual-cell-challenge-2026-wrap-up`, published 2026-11-19, fetched HTTP 200):
   Verbatim excerpt from *"Where models fell short"*:
   > "On the validation lines, the top ten entries predicted nearly the same response to a given knockdown in every line. Among the genes that respond significantly to a knockdown in a given line, the cosine similarity between the predicted responses to a target in different lines ranged from 0.77 to 0.89, compared with 0.55 for the measured responses and 0.91 to 0.95 for two independent samples of cells from the same line, and the line-specific component of the response that the models did predict was small (a slope of 0.02 to 0.13, compared with 1.02 for a replicate). On the same genes, the predicted responses agreed with the measured ones at 0.38 to 0.52, compared with 0.55 for the response to the same knockdown measured in a different validation line. Together, these observations indicate that current models recover roughly what a cross-line average would predict, whereas the part of the response that is specific to a cell line, which is precisely what a virtual cell model would need in order to improve on an experiment in another model system, remains largely missing."

2. **`ArcInstitute/cell-eval2` Official Repository & Docs** (`v0.22.0`, released 2026-11-06; `docs/vcc2026_metrics/vcc2026-metrics.md` and `vcc2026-metrics-brief.md` dated 2026/08/19; `docs/metrics.md`; `src/cell_eval2/configs/vcc2026.yaml`; maintainer: Alexander Dobin, `alexander.dobin@arcinstitute.org`):
   - Directly linked by `https://virtualcellchallenge.org/evaluation` as the canonical specification for the **six scored competition metrics** (`pds_cosine`, `expr_mse_unbiased_capped_norm`, `de_wilcoxon_direction_fidelity_yield_raw`, `de_wilcoxon_direction_reach_raw`, `de_wilcoxon_sig_jaccard`, `de_wilcoxon_lfc_nmae`).
   - **Negative Verification Fact**: Exhaustive inspection of `cell-eval2` (`docs/`, `src/cell_eval2/`, `tests/`, `CHANGELOG.md`) confirms that `cell-eval2` implements the 6 competition metrics (§1–§6) and 3 unscored within-perturbation single-cell covariance diagnostics (`§8`: `cell_cov_corr`, `cell_cov_eig_ratio`, `cell_cov_ccc`), and does **NOT** implement or document the Wrap-Up blog's cross-line *"line-specific component ... (a slope of 0.02 to 0.13, compared with 1.02 for a replicate)"* or cross-line cosine similarity.

3. **Alexander Dobin's 2026 Metrics Preprint** (`https://www.biorxiv.org/content/10.64898/2026.11.17.767214`):
   - Linked in the Wrap-Up under *"The six metrics"* (`"Full definitions are available on the evaluation page and a new preprint prepared by the lead architect of the Challenge metrics, Alexander Dobin, Arc's Director of Bioinformatics."`) and *"Looking ahead"* (`"2026 metrics preprint"`).
   - **API Verification Audit (2026-11-21)**: Checked `api.biorxiv.org` (both `biorxiv` and `medrxiv`, including full 2026-11-16..2026-11-21 `bioinformatics` [186 records] and `genomics` [71 records] dumps), `api.crossref.org`, `api.openalex.org`, EuropePMC OA API, and PubMed (`Dobin[Author] AND 2026[PDAT]`). As of 2026-11-21, `10.64898/2026.11.17.767214` is **not yet indexed** (`no posts found` / `404`), and `www.biorxiv.org` HTML returns HTTP 403 Cloudflare challenge (not bypassed per rules).

### 1.2 Explicit Gap Statement for the "Line-Specific Slope"
Neither the Wrap-Up HTML nor `cell-eval2` prints an algebraic equation or code link for the Wrap-Up's post-hoc `0.02–0.13` (replicate `1.02`) line-specific slope. What Arc's prose **explicitly fixes** vs. what Team A **must lock by preregistration** is:

| Component | Explicitly Fixed by Arc Wrap-Up / `cell-eval2` | Underspecified in Arc Prose (Locked Below for A-T2/A-T3) |
|---|---|---|
| **Core Quantity** | Regression slope of **predicted** line-specific component on **measured** line-specific component (so cross-line-mean predictor $\to 0$, attenuated predictions $\to 0.02\text{--}0.13$, biological replicate $\to \approx 1.0$) | Intercept included ($\hat{\beta}_{\mathrm{OLS}}$) vs. forced zero intercept ($\hat{\beta}_0$); pooled over $(p, g)$ vs. perturbation-averaged $\bar{\beta}_{\mathrm{pert}}$ |
| **Line-Specific Decomposition** | Response minus cross-line average response (*"recover roughly what a cross-line average would predict, whereas the part ... specific to a cell line ... remains largely missing"*) | **Leave-target-out** source mean $\bar{\Delta}^{\mathrm{src}}_p = \frac{1}{|\mathcal{S}_p|}\sum_{s \in \mathcal{S}_p}\Delta^{\mathrm{obs}}_{p,s}$ (leak-free; sets Transferred Mean baseline strictly to $0$) vs. all-lines mean |
| **Gene Subset** | *"Among the genes that respond significantly to a knockdown in a given line"*; all `cell-eval2` metrics exclude the perturbed target gene $p$ itself | In `p2/shared/lines/<line>/pseudobulk_delta.npz` (where single-cell perturbed counts are not stored uniformly across Replogle/Nadig/etc.), define responding genes deterministically from $\Delta^{\mathrm{obs}}_{p,t}$ alone (e.g., Top-50 $|\Delta^{\mathrm{obs}}_{p,t,g}|$, $g \neq p$), plus all shared non-target genes $\mathcal{G}\setminus\{p\}$ |
| **Expression Space** | `cell-eval2` uses `bulk_lognorm` (`TS = 5e4`) for expression/PDS and `log2FC` for DE | `p2/shared/lines/` uses Addendum-B standardized pseudobulk $\Delta$ (`log1p` normalized + control z-scored, matching Replogle 2022 & Nadig 2025) |

---

## 2. Locked Executable Specification for Team A (`A-T2` & `A-T3`)

Let $t$ be the held-out target cell line (zero-shot: only unperturbed controls and target list are observable for $t$), and let $\mathcal{S}_p$ be the set of source cell lines in which perturbation $p$ is measured. Let $\mathcal{G}$ be the shared readout gene panel (`COMMON_GENES.txt`).

For each evaluated perturbation $p \in \mathcal{P}_t$ in target line $t$:
- $\Delta^{\mathrm{obs}}_{p, t} \in \mathbb{R}^{|\mathcal{G}|}$: measured pseudobulk response in target line $t$ relative to target line $t$'s controls.
- $\hat{\Delta}_{p, t} \in \mathbb{R}^{|\mathcal{G}|}$: model-predicted response in target line $t$.
- **Transferred cross-line mean baseline** (leave-target-line-out, zero leakage):
  $$\bar{\Delta}^{\mathrm{src}}_{p, (-t)} = \frac{1}{|\mathcal{S}_p|} \sum_{s \in \mathcal{S}_p} \Delta^{\mathrm{obs}}_{p, s}.$$

### 2.1 Line-Specific Residual Decomposition
Define the **true** and **predicted** line-specific residuals for $(p, t)$ as:
$$r_{p, t} = \Delta^{\mathrm{obs}}_{p, t} - \bar{\Delta}^{\mathrm{src}}_{p, (-t)}, \qquad \hat{r}_{p, t} = \hat{\Delta}_{p, t} - \bar{\Delta}^{\mathrm{src}}_{p, (-t)}.$$

*Boundary Check*: For the Transferred Mean baseline $\hat{\Delta}_{p,t} = \bar{\Delta}^{\mathrm{src}}_{p, (-t)}$, we have $\hat{r}_{p,t} \equiv \mathbf{0}$ identically, so its line-specific slope is exact `0.0000`.

### 2.2 Gene Mask $\mathcal{M}_{p,t}$ (Mandatory Target-Gene Exclusion)
For every perturbation $p$ in target line $t$:
1. Exclude the perturbed target gene $p$ itself ($g \neq p$) and any non-finite gene entries (`~isfinite`).
2. **Primary mask (`resp_top50`, Wrap-Up-aligned responding genes)**: Let $\mathcal{M}^{\mathrm{top50}}_{p,t}$ be the top $K = 50$ non-target genes ranked by $|\Delta^{\mathrm{obs}}_{p, t, g}|$ in target line $t$ (matching the median DEG breadth `22–44` on VCC validation and `25–102` on VCC test, and matching `de50` in Preprint v2.1; if $|\mathcal{G} \setminus \{p\}| < 50$, take all valid non-target genes).
3. **Full-panel sensitivity mask (`all_nontarget`)**: $\mathcal{M}^{\mathrm{all}}_{p,t} = \{g \in \mathcal{G} : g \neq p,\ \text{finite}\}$.
*(If a dataset also provides per-gene FDR $< 0.05$ from single-cell DE, `sig_fdr05` may be logged as a secondary diagnostic, provided the cohort size meets `min_gate_size >= 10`.)*

### 2.3 Line-Specific Slope Estimators
Collect all masked pairs $(x_k, y_k) = (r_{p,t,g},\; \hat{r}_{p,t,g})$ across $p \in \mathcal{P}_t$ and $g \in \mathcal{M}_{p,t}$, for $k = 1, \dots, N$ ($N = \sum_{p \in \mathcal{P}_t} |\mathcal{M}_{p,t}|$):

1. **Primary Metric (`line_specific_slope`, pooled OLS slope with intercept on `resp_top50`)**:
   $$\hat{\beta}_{\mathrm{OLS}} = \frac{\sum_{k=1}^{N} (x_k - \bar{x})(y_k - \bar{y})}{\sum_{k=1}^{N} (x_k - \bar{x})^2}, \qquad \bar{x} = \frac{1}{N}\sum_{k=1}^N x_k,\ \bar{y} = \frac{1}{N}\sum_{k=1}^N y_k.$$
   *(Return `NaN` if $\sum (x_k - \bar{x})^2 \le 10^{-12}$.)*

2. **Mandatory Sensitivity Variants (log alongside primary, never substitute post-hoc)**:
   - **Zero-intercept pooled OLS slope** ($\hat{\beta}_0$):
     $$\hat{\beta}_0 = \frac{\sum_{k=1}^N x_k y_k}{\sum_{k=1}^N x_k^2}.$$
   - **Per-perturbation mean OLS slope** ($\bar{\beta}_{\mathrm{pert}}$):
     $$\bar{\beta}_{\mathrm{pert}} = \frac{1}{|\mathcal{P}_t|} \sum_{p \in \mathcal{P}_t} \frac{\sum_{g \in \mathcal{M}_{p,t}} (r_{p,t,g} - \bar{r}_{p,t})(\hat{r}_{p,t,g} - \bar{\hat{r}}_{p,t})}{\sum_{g \in \mathcal{M}_{p,t}} (r_{p,t,g} - \bar{r}_{p,t})^2}.$$
   - Same three quantities on `all_nontarget` ($\mathcal{M}^{\mathrm{all}}_{p,t}$).

### 2.4 Line-Specific Noise-Calibrated Correlation (for Hypothesis H-A2)
1. **Raw line-specific Pearson correlation**:
   $$r_{\mathrm{resid}} = \frac{\sum_{k=1}^N (x_k - \bar{x})(y_k - \bar{y})}{\sqrt{\sum_{k=1}^N (x_k - \bar{x})^2 \sum_{k=1}^N (y_k - \bar{y})^2}}.$$
   *(Also compute perturbation-averaged $\bar{r}_{\mathrm{resid,pert}} = \frac{1}{|\mathcal{P}_t|}\sum_p \mathrm{corr}_g(r_{p,t},\, \hat{r}_{p,t})$ on $\mathcal{M}_{p,t}$.)*

2. **Split-Half Noise Ceiling (Spearman–Brown, aligned with `cell-eval2` §6b)**:
   When target line $t$ has replicate screens (e.g., Replogle K562 Day 6 vs. Day 8) or single-cell split-halves $(A, B)$ with **disjoint control cells** per half (per `cell-eval2` §6b so control sampling noise is not shared across halves):
   - Form $r^{(A)}_{p,t} = \Delta^{(A)}_{p,t} - \bar{\Delta}^{\mathrm{src}}_{p, (-t)}$ and $r^{(B)}_{p,t} = \Delta^{(B)}_{p,t} - \bar{\Delta}^{\mathrm{src}}_{p, (-t)}$ using the **independent** source mean $\bar{\Delta}^{\mathrm{src}}_{p, (-t)}$.
   - Compute half-split correlation $\rho_{AB} = \mathrm{corr}(r^{(A)}, r^{(B)})$ on $\mathcal{M}_{p,t}$ (evaluated at the true-full mask or each half's mask; record exact convention).
   - Following `cell-eval2` §6b, **only if $\rho_{AB} > 0$**, apply Spearman–Brown doubling to estimate full-depth reliability:
     $$\rho_{\mathrm{SB}} = \frac{2\,\rho_{AB}}{1 + \rho_{AB}}, \qquad \text{ceiling}_{r} = \sqrt{\rho_{\mathrm{SB}}}.$$
     *(If $\rho_{AB} \le 0$, report ceiling as `NaN`: no positive repeatability evidence for the line-specific residual at this depth.)*
   - **Noise-calibrated line-specific correlation**:
     $$r_{\mathrm{resid,cal}} = \frac{r_{\mathrm{resid}}}{\sqrt{\rho_{\mathrm{SB}}}}.$$
     *(Note for A-T3: Because $\sqrt{\rho_{\mathrm{SB}}} > 0$ is a positive line-level scalar, if $\rho_{\mathrm{SB}}$ is held fixed as a line constant during row bootstrap, `CI(r_resid_cal) > 0` is algebraically equivalent to `CI(r_resid) > 0`. If $\rho_{\mathrm{SB}}$ is also resampled inside each perturbation bootstrap draw when paired half-splits are available per perturbation, the CI reflects ceiling uncertainty too — PREREG_A.md must state which of the two is locked.)*
   - **Split-half replicate slope benchmark**: Also report $\hat{\beta}_{\mathrm{rep}} = \mathrm{OLS}(r^{(B)} \sim r^{(A)})$ (and averaged $A \leftrightarrow B$); note that OLS attenuation $\mathrm{Var}(r^\star)/(\mathrm{Var}(r^\star)+\sigma_{\mathrm{half}}^2)$ typically places split-half OLS slopes below $1.0$, whereas full-vs-replicate or Deming slope sits near $1.0$.

### 2.5 Row-Level (Perturbation-Level) Paired Bootstrap Inference
- **Unit of resampling**: Perturbation $p \in \mathcal{P}_t$ (rows), NOT individual $(p,g)$ entries, because genes within a perturbation share strong co-expression correlation.
- **Resamples**: $B = 4{,}000$ deterministic draws of $M = |\mathcal{P}_t|$ perturbations with replacement (matching Preprint v2.1 `compute_row_bootstrap_ci.py`).
- In each bootstrap draw $b$, concatenate the pre-masked vectors $(r_{p,t}[\mathcal{M}_{p,t}],\; \hat{r}_{p,t}[\mathcal{M}_{p,t}])$ for the $M$ sampled perturbations and recompute $\hat{\beta}_{\mathrm{OLS}}^{(b)}$, $r_{\mathrm{resid}}^{(b)}$, and pairwise model differences $\Delta \theta^{(b)} = \theta_{\mathrm{model1}}^{(b)} - \theta_{\mathrm{model2}}^{(b)}$. Report 2.5th and 97.5th percentiles.

### 2.6 Wrap-Up Companion Diagnostics (Cosine Agreement & Cross-Line Similarity)
For completeness with the Wrap-Up paragraph:
- **Predicted-vs-measured cosine agreement** on $\mathcal{M}_{p,t}$: $\frac{1}{|\mathcal{P}_t|}\sum_{p \in \mathcal{P}_t} \cos(\hat{\Delta}_{p,t}[\mathcal{M}_{p,t}],\; \Delta^{\mathrm{obs}}_{p,t}[\mathcal{M}_{p,t}])$.
- **Predicted cross-line cosine similarity** (across pairs of held-out lines $t \neq t'$ on shared targets): $\cos(\hat{\Delta}_{p,t},\; \hat{\Delta}_{p,t'})$ vs. measured $\cos(\Delta^{\mathrm{obs}}_{p,t},\; \Delta^{\mathrm{obs}}_{p,t'})$ (Arc observed `0.77–0.89` for top-10 models vs. `0.55` measured).

---

## 3. Official `cell-eval2` `vcc2026` 6-Metric Suite (`v0.22.0`, `configs/vcc2026.yaml`)

For A-T4 prospective validation and reference parity against `ArcInstitute/cell-eval2` (`rule_version: 3`, `docs/vcc2026_metrics/vcc2026-metrics.md` & `vcc2026-metrics-brief.md`):

### 3.1 Input Contract & Shared Parameters (`vcc2026.yaml`)
- **Submission contract**: Non-negative **raw integer counts** (`input_type: counts`, `max_counts_per_cell: 1e6`), **400 cells per perturbation** per context (`vcc prep` upload rule; non-targeting rows forbidden on upload and supplied from real controls at scoring via `control_source: real`), `18,533` genes matching `gene_names.csv`.
- **Expression comparator (`bulk_lognorm`)**: Group-sum counts $P_{p,g} = \sum_{c \in \mathcal{C}_p} y_{c,g}$, normalized to `bulk_target_sum` $\mathrm{TS} = 5 \times 10^4$:
  $$b_{p,g} = \log\!\left(1 + \mathrm{TS}\,\frac{P_{p,g}}{\sum_{g'} P_{p,g'}}\right), \qquad \delta_q = b_q - b_{\mathrm{ctrl}},\qquad \hat{\delta}_p = \hat{b}_p - b_{\mathrm{ctrl}}.$$
- **DE table (`de_wilcoxon_*`)**: Two-sided Wilcoxon rank-sum test (normal approximation) on per-cell `1e6`-normalized (CPM) non-logged counts vs. real control (`control_source: real`), filtered to genes with **control mean > 5 CPM** (`filter_gene_min_cpm_cell: 5.0`), per-perturbation BH FDR `< 0.05` (`fdr_scope: per_pert`), fold change $\mathrm{lfc}_{p,g} = \log_2 \frac{m^{\mathrm{pert}}_{p,g} + 10^{-9}}{m^{\mathrm{ctrl}}_g + 10^{-9}}$ (`mean_calc: arithmetic`, `epsilon: 1e-9`).

### 3.2 The Six Scored Metrics & Official Reference Bounds

| Metric ID | Formula Summary | Target Exclusion | $b$ (Context Mean) | $r$ (5-Split Half) | Clamping / Cohort |
|---|---|---|---|---|---|
| `pds_cosine` | $1 - \frac{k_p}{n-1}$, midrank tie on cosine distance $d_{p,q}$ | **All 300 panel targets** (`panel`) | `0.500` | `0.927 – 0.984` | Unclamped; $n=300$ |
| `expr_mse_unbiased_capped_norm` | $\frac{\sum_p N_p}{\sum_p D_p}$ with delete-one jackknife $C_p, \hat{C}_p$ and spread cap $\rho \le 1$ | Own target ($G_p$) | `0.986 – 0.992` | `0.028 – 0.045` | **Clamped `[0, 1]`**; panel ratio |
| `de_wilcoxon_direction_fidelity_yield_raw` | $\frac{k}{\max(n_{\mathrm{pred}}, n_{\mathrm{real}})}$ over adjudicable calls | Own target | `0.505 – 0.522` | `0.795 – 0.832` | Unclamped; `295–300` |
| `de_wilcoxon_direction_reach_raw` | $\frac{k^\star}{n_{\mathrm{real}}}$, $k^\star = \max\{k : P(k) \ge 0.9\}$ (sig first, $p_{\mathrm{adj}} \uparrow$, $p \uparrow$, $|\mathrm{lfc}| \downarrow$) | Own target | `0.047 – 0.097` | `0.958 – 0.978` | Unclamped; `290–300` |
| `de_wilcoxon_sig_jaccard` | $\frac{|R_p \cap \hat{R}_p|}{|R_p \cup \hat{R}_p|}$ (empty union $= 1$) | Own target | `0.021 – 0.037` | `0.375 – 0.423` | Unclamped; `300` |
| `de_wilcoxon_lfc_nmae` | $\frac{\sum_{g \in S_p} |\widehat{\mathrm{lfc}}_{p,g} - \mathrm{lfc}_{p,g}|}{\sum_{g \in S_p} |\mathrm{lfc}_{p,g}|}$ (missing/non-finite $\to 0$) | Own target | `1.0009 – 1.0017` | `0.369 – 0.431` | **Floor `-6`**; $|S_p| \ge 10$ (`209–261`) |

- **Per-context scaled score**: $s = \frac{u - b}{r - b}$ (unweighted mean of the 6 scaled metrics across the 3 contexts).

---

## 4. Executable Reference Implementation (for A-T2 & A-T3 Cross-Check)

```python
import numpy as np

def compute_line_specific_metrics(
    delta_obs: np.ndarray,       # (P, G) target-line measured pseudobulk deltas
    delta_pred: np.ndarray,      # (P, G) predicted deltas for target line
    delta_src_mean: np.ndarray,  # (P, G) leave-target-out mean across source lines
    target_gene_idx: np.ndarray, # (P,) int index into G (-1 if target not in readout)
    top_k: int = 50,
    n_boot: int = 4000,
    seed: int = 42,
) -> dict:
    """Reference implementation of the VCC 2026 Wrap-Up line-specific slope & residual metrics."""
    P, G = delta_obs.shape
    r_true = delta_obs - delta_src_mean
    r_pred = delta_pred - delta_src_mean

    def _masked_pairs(mode: str):
        xs_list, ys_list, per_pert_slopes, per_pert_corrs = [], [], [], []
        for p in range(P):
            valid = np.isfinite(r_true[p]) & np.isfinite(r_pred[p]) & np.isfinite(delta_obs[p])
            t_idx = int(target_gene_idx[p])
            if 0 <= t_idx < G:
                valid[t_idx] = False
            idx = np.flatnonzero(valid)
            if idx.size == 0:
                continue
            if mode == "resp_top50" and idx.size > top_k:
                order = np.argsort(np.abs(delta_obs[p, idx]))[-top_k:]
                idx = idx[order]
            xp = r_true[p, idx].astype(np.float64)
            yp = r_pred[p, idx].astype(np.float64)
            xs_list.append(xp)
            ys_list.append(yp)
            xv = xp - xp.mean()
            yv = yp - yp.mean()
            denom = float(np.dot(xv, xv))
            if denom > 1e-12:
                per_pert_slopes.append(float(np.dot(xv, yv) / denom))
            yn = float(np.dot(yv, yv))
            if denom > 1e-12 and yn > 1e-12:
                per_pert_corrs.append(float(np.dot(xv, yv) / np.sqrt(denom * yn)))
        return xs_list, ys_list, per_pert_slopes, per_pert_corrs

    def _pooled_stats(xs_list, ys_list):
        x = np.concatenate(xs_list)
        y = np.concatenate(ys_list)
        xc, yc = x - x.mean(), y - y.mean()
        vx = float(np.dot(xc, xc))
        vy = float(np.dot(yc, yc))
        sxx = float(np.dot(x, x))
        slope_ols = float(np.dot(xc, yc) / vx) if vx > 1e-12 else float("nan")
        slope_no_int = float(np.dot(x, y) / sxx) if sxx > 1e-12 else float("nan")
        corr = float(np.dot(xc, yc) / np.sqrt(vx * vy)) if (vx > 1e-12 and vy > 1e-12) else float("nan")
        return slope_ols, slope_no_int, corr

    out = {}
    rng = np.random.default_rng(seed)
    for mode in ("resp_top50", "all_nontarget"):
        xs_list, ys_list, p_slopes, p_corrs = _masked_pairs(mode)
        slope_ols, slope_no_int, corr = _pooled_stats(xs_list, ys_list)
        m = len(xs_list)
        boot_ols, boot_noint, boot_corr = np.empty(n_boot), np.empty(n_boot), np.empty(n_boot)
        for b in range(n_boot):
            draw = rng.integers(0, m, size=m)
            b_ols, b_noint, b_c = _pooled_stats([xs_list[i] for i in draw], [ys_list[i] for i in draw])
            boot_ols[b], boot_noint[b], boot_corr[b] = b_ols, b_noint, b_c
        out[mode] = {
            "slope_ols": slope_ols,
            "slope_ols_ci95": np.nanpercentile(boot_ols, [2.5, 97.5]).tolist(),
            "slope_no_intercept": slope_no_int,
            "slope_no_intercept_ci95": np.nanpercentile(boot_noint, [2.5, 97.5]).tolist(),
            "per_pert_slope_mean": float(np.mean(p_slopes)) if p_slopes else float("nan"),
            "resid_corr": corr,
            "resid_corr_ci95": np.nanpercentile(boot_corr, [2.5, 97.5]).tolist(),
            "per_pert_corr_mean": float(np.mean(p_corrs)) if p_corrs else float("nan"),
        }
    return out
```

---

## 5. Verification Log (URLs & APIs Checked on 2026-11-21)
- Arc VCC 2026 Wrap-Up: `https://arcinstitute.org/news/virtual-cell-challenge-2026-wrap-up` (HTTP 200; full text archived & verified).
- Arc VCC 2026 Evaluation Page: `https://virtualcellchallenge.org/evaluation` (Next.js static bundle inspected; links `ArcInstitute/cell-eval2/blob/main/docs/vcc2026_metrics/vcc2026-metrics.pdf`).
- `ArcInstitute/cell-eval2` (`v0.22.0`): `https://github.com/ArcInstitute/cell-eval2/archive/refs/heads/main.tar.gz` (`docs/vcc2026_metrics/vcc2026-metrics.md`, `vcc2026-metrics-brief.md`, `docs/metrics.md`, `src/cell_eval2/configs/vcc2026.yaml`, `src/cell_eval2/metrics/{discrimination,delta,de}.py`).
- Alexander Dobin Preprint (`doi:10.64898/2026.11.17.767214`): checked `api.biorxiv.org`, `api.crossref.org`, `api.openalex.org`, EuropePMC OA API, PubMed E-utilities (`Dobin[Author] AND 2026[PDAT]`) — currently unindexed (`no posts found` / `404`); `www.biorxiv.org` HTTP 403 Cloudflare challenge (not bypassed).
