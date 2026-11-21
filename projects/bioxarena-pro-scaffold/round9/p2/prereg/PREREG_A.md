# Pre-Registration A (`PREREG_A`): Zero-Shot Cross-Cell-Line CRISPRi Perturbation Prediction via Target-Line Unperturbed Control Co-Fluctuations

**Study Team**: BioXArena Team A (`A-T1` Data, `A-T2` Solver/Eval, `A-T3` Pre-Registration & Statistical Audit, `A-T4` Prospective, `A-T5` Literature & Metrics)  
**Principal Investigator**: Zijian Zeng (Round 9 Phase 2, AI-Biology-Lab)  
**Date Locked**: 2026-11-21 (Pre-Evaluation Lock; awaiting master commit recorded in `round9/p2/status/PREREG_COMMIT.txt`)  
**Upstream Reference Rules**: `round9/p2/P2_RULES.md`, `round9/p2/cards/A-T3.md`  
**Prior-Art & Official Metric Specs**: `round9/p2/vcc_linespecific/audit/LIT/{LINE_SPECIFIC_METRIC.md,CIPHER_METHOD.md,PRIOR_ART_A.md,COMPETITION_MONITOR.md}`

---

## 1. Falsifiable Research Question & Scientific Hypothesis

In the Arc Virtual Cell Challenge 2026 Wrap-Up (2026-11-19), top-10 challenge models achieved high bulk/global agreement yet recovered a **cell-line-specific perturbation slope of only `0.02–0.13`** (versus **`1.02`** for experimental replicates). This indicates that existing cross-context perturbation predictors predominantly capture context-invariant shared responses while missing cell-line-specific downstream propagation.

**Core Falsifiable Question**: Under strict **zero-shot cross-cell-line generalization**—using **zero perturbation labels or perturbed cells from the held-out target cell line**—can the target cell line's unperturbed single-cell control gene–gene co-fluctuation covariance $C_t$ recover the cell-line-specific residual component $\Delta_{p,t} - \bar{\Delta}^{\mathrm{src}}_{p,(-t)}$ missed by cross-line mean transfer (`Transfer_Mean`) without degrading overall perturbation discrimination (`disc_raw`) or differentially expressed gene correlation (`de50_r`)?

---

## 2. Rule 3 Zero-Look & Rule 5 Clean-Data Guarantee

1. **Zero-Look Prior to Commit Lock**: Prior to the creation of `round9/p2/status/PREREG_COMMIT.txt`, **no team member or script evaluated, tuned, inspected, or loaded expression matrices (`.npz`) from any held-out cell line** (`nadig_hepg2`, `nadig_jurkat`, `nourreddine_kolf21j`, `vcc2025_h1`). Preprocessing sanity checks, numerical cross-checks (`A-T2` vs `A-T3`), and locked predictor/representation choices were performed **exclusively on the Replogle et al. (2022) `K562 <-> RPE1` development pair** (`dev_sweep.py` / `--dev-only`).
2. **Open-Access Clean Data & Zero X-Atlas**: Every dataset is open-access (`CC BY 4.0`, `CC0 1.0`, or public NCBI GEO), verified by URL, byte size, and SHA-256 digest in `round9/p2/shared/DATA_SOURCES.md`. Total newly downloaded volume across the team is **`56.597333` GB** (`56,597,332,700` B, strictly `<= 60.000` GB cap). Zero X-Atlas / Orion datasets, derivatives, or gene panels are used.

---

## 3. Pre-Registered Dataset Inventory & Standardized Preprocessing

### 3.1 Evaluated Cell Lines ($K = 6$)

| Cell Line (`p2/shared/lines/`) | Source & License | Raw Input SHA-256 | Raw Cells × Genes | NT Ctrl (Sampled / Total) | Aligned Final Genes | Dedup Targets (All / Rep) | In `COMMON_GENES` (All / Rep) | LOLO Targets >=1 Src (All / Rep-Target-in-CG) |
|---|---|---|---|---:|---:|---:|---:|---:|
| `replogle_k562_essential` *(dev pair)* | Replogle et al. (2022), Figshare+ 20029387 (`CC BY 4.0`) | `3e5a63a9e892b21029bb55fca4e12517a49aad7af6c14133ca63d12cf68c6cee` | 310,385 × 8,563 | 10,691 / 10,691 | 8,561 | 2,049 / 1,113 | 1,558 / 910 | 2,039 / 859 |
| `replogle_rpe1` *(dev pair)* | Replogle et al. (2022), Figshare+ 20029387 (`CC BY 4.0`) | `9b05ef1f81526216fa008d677e9e0d03dce9a2f7a95499a4fb81e505e9d88ef1` | 247,914 × 8,749 | 11,485 / 11,485 | 8,746 | 2,351 / 1,575 | 1,780 / 1,290 | 2,350 / 1,070 |
| `nadig_hepg2` *(held-out)* | Nadig et al. (2025), GEO `GSE264667` (Open GEO) | `e1ad7c3c5a201c861a207a858aa7e59f5e6ac1955674c415f7de0d1dadadb52e` | 145,473 × 9,624 | 4,976 / 4,976 | 9,432 | 2,343 / 1,000 | 1,777 / 811 | 2,343 / 799 |
| `nadig_jurkat` *(held-out)* | Nadig et al. (2025), GEO `GSE264667` (Open GEO) | `ffbe15f2c8f7ffcfd7b0ba9e6937d4ebc2d03b0179fa8234648a59bcb82c04a3` | 262,956 × 8,882 | 12,013 / 12,013 | 8,826 | 2,350 / 1,262 | 1,779 / 995 | 2,350 / 934 |
| `nourreddine_kolf21j` *(held-out)* | Nourreddine et al. (2026), Figshare+ 27261219 (`CC BY 4.0`; `CM+ME Strong`) | `741674b8a1...3361745` (`CM`), `d86fa3260c...1d59211` (`ME`) *(exact digests in `DATA_SOURCES.md`)* | 38,322 × 34,583 | 17,895 / 17,895 | 15,000 | 95 / 67 | 69 / 51 | 32 / 18 |
| `vcc2025_h1` *(held-out)* | Arc VCC 2025 H1 hESC (`Train+Val+Test`, `CC0 1.0`) | `a09977104f...19c02b` (`Train`), `376f0bab27...77bfe1` (`Val`), `ba5ce66729...31f03` (`Test`) | 414,694 × 18,080 *(76,352 dup Val/Test NT clones removed)* | 20,000 / 38,176 | 14,656 | 300 / 204 | 228 / 168 | 106 / 80 |

- **Aligned Common Readout Genes (`round9/p2/shared/COMMON_GENES.txt`)**: **$G = 6,148$** HGNC symbols shared across all 6 lines (after HGNC mapping via `hgnc_complete_set.txt` `sha256=6df77ac3cbaeb86a6dce5e17d851ebb27b6b3cdcd2471a211a3d8464f00c54c4`, dropping empty/duplicate/non-finite gene columns and capping `nourreddine_kolf21j` at 15,000 highest-expressed genes).
- **Budget Disclosure (`nourreddine_kolf21j`)**: `KOLF_Strong_Perturbations.h5ad` (`46.719 GB`) and `KOLF_Pan_Genome_QC_Filtered.h5ad` (`189.393 GB`) exceed the `60.00 GB` team cap alongside `Replogle` (`19.363 GB`) + `VCC 2025 H1` (`34.362 GB`); therefore `nourreddine_kolf21j` uses the two Strong sub-atlases (`2.855 GB`), sharing 32 all / 18 reproducible targets in `COMMON_GENES`.
- **Disjoint Deduplication & Split-Halves**: Multiple sgRNAs per target gene are aggregated per cell line (`min >= 10` cells). Within each batch, perturbed cells and NT control cells are deterministically partitioned 50/50 into disjoint halves (`half_a`, `half_b`; zero shared control or perturbed cells between Half A and Half B; `half_valid = True` on 100% of deduplicated targets across all 6 lines).
- **Reproducibility Criterion (`reproducible`)**: Perturbations satisfying `n_cells >= 25`, `ad_counts >= 50`, and `energy_test_p_value < 0.001`.

### 3.2 Locked Primary Configuration & Sensitivity Suite (Locked from Replogle `K562 <-> RPE1` Dev Only)

- **Primary Representation (`--repr` / `--space`)**: **`log1p`** (`log1p(1e4 * UMI / libsize)` minus per-batch NT control mean)
- **Primary Covariance Shrinkage (`--shrink`)**: **`0.0`** (diagonal variance floor `max(diag(C), 1e-8)` without matrix inversion)
- **Control Sample Scope (`cov_max_cells`)**: All NT control cells stored in `ctrl_cells.npz` (`4,976–20,000` per line)
- **Self-Target Masking (`mask_self`)**: **`True`** ($v_\ell(p)_{g_p} = 0$ on target gene $g_p$ so self-knockdown depth never inflates slope or discrimination)
- **Primary Evaluation Target Filter (`only_reproducible`)**: **`True`** — evaluate targets $p \in \mathcal{P}_t$ in `COMMON_GENES` that pass `reproducible == True` on target line $t$ and have at least one `reproducible == True` source line $s \ne t$ (averaging $\bar{\Delta}^{\mathrm{src}}_{p,(-t)}$ and $\bar{v}_{\mathcal{S}_p}(p)$ over reproducible source lines $\mathcal{S}_p$).
- **Primary Zero-Shot Model**: **`Hybrid_ProjShared`** (reporting `Hybrid_Raw`, `CIPHER_Linear`, `Transfer_Mean`, and `Source_Grand_Mean` alongside).
- **Mandatory Pre-Registered Sensitivity Checks *(secondary, non-confirmatory)***:
  - **S1**: `only_reproducible = False` (all deduplicated shared targets, expanding `nourreddine_kolf21j` and `vcc2025_h1` counts)
  - **S2**: `Hybrid_Raw` vs `Hybrid_ProjShared` comparison
  - **S3**: `shrink = 1e-3` covariance regularization check
  - **S4**: `cov_max_cells = 10,000` subsample check & `zscore` representation sensitivity

**Dev Verification Rationale (`Replogle K562 <-> RPE1` Only)**:  
On the 7-config `Replogle` dev sweep (`7,225` common genes, `1,761` shared all / `832` shared reproducible targets):
1. `log1p` outperformed `zscore` (which suppressed `Hybrid_ProjShared` line-specific slope $>3\times$ to `0.0011–0.0016` and lowered `K562` `cent_de50_r` by `-0.011`) and `libnorm` (which exploded MSE $>100\times$ to `1.46–1.72` and induced a `120\times` $\lambda$ instability `-2.80` vs `-0.023`).
2. `Hybrid_ProjShared` eliminated unprojected `Hybrid_Raw`'s $-\bar{v}_{\mathcal{S}_p}$ shared-axis leak (`NegCtrl_Hybrid_DiagCov` `+0.0206` vs `NegCtrl_HybridProj_DiagCov` `-0.0005` on `K562` `all`), passing **H-A2** (`cal_r` `+0.0523 [0.0333, 0.0715]` on `K562`, `+0.1054 [0.0580, 0.1526]` on `RPE1` for `rep1`), **H-A3** (`2/2` disjoint CIs vs `MismatchedCov`, `ShuffledGenes`, and `DiagCov`), and **H-A4** (`4/4` non-inferiority sub-metrics confirmed on `rep1`, improving `de50_r` by `+0.0026` to `+0.0373`), while showing transparently on dev that **H-A1** (`slope` `0.0050–0.0085`) stays far below `0.13` because unscaled zero-shot covariance residual spread is $\sim 15\times$ smaller than observed residual spread (`rule 8` respected: no post-hoc rescaling).

---

## 4. Exact LOO-CL Zero-Shot Predictors & Negative Controls

For each held-out target cell line $t \in \{1, \dots, 6\}$ and remaining source lines $\mathcal{S} = \{1, \dots, 6\} \setminus \{t\}$:

### 4.1 Zero-Shot Information Boundary
Target line $t$ exposes **strictly** its unperturbed non-targeting control matrix $X_t^{\mathrm{ctrl}} \in \mathbb{R}^{N_t^{\mathrm{ctrl}} \times G}$ and target gene names $\mathcal{P}_t$. **Zero** perturbed cells or perturbation deltas $\Delta_{p,t}$ from line $t$ enter covariance estimation, shared-axis SVD, or $\lambda_{(-t)}$ estimation.

### 4.2 Mathematical Formulations

1. **Single-Line Unperturbed Covariance Linear Response** (`CIPHER_METHOD.md`, adapted from Kuznets-Speck et al. 2025/2026 bioRxiv `10.1101/2025.06.27.661814` v2; CIPHER v2 fits native per-perturbation `a_hat` on perturbed target data, whereas Team A tests strictly zero-shot normalized control transfer):
   $$v_\ell(p) = \frac{C_{\ell,:, g_p}}{\max(C_{\ell, g_p, g_p},\, 10^{-8})}, \qquad v_\ell(p)_{g_p} \leftarrow 0.$$

2. **Source-Only Signed Scalar Calibration $\hat{\lambda}_{(-t)}$**:
   Fitted across source lines $s \in \mathcal{S}$ on non-target genes (`target` excluded):
   $$\hat{\lambda}_{(-t)} = \mathrm{clip}\!\left(\frac{\sum_{s \in \mathcal{S}} \sum_{p \in \mathcal{P}_s} \langle \Delta_{p,s},\, v_s(p) \rangle}{\max\!\left(\sum_{s \in \mathcal{S}} \sum_{p \in \mathcal{P}_s} \|v_s(p)\|_2^2,\, 10^{-12}\right)},\; -10.0,\; 10.0\right).$$
   Under CRISPRi knockdown ($\Delta_{g_p} < 0$), positively co-varying genes in $v_s(p)$ also decrease ($\Delta_k < 0$), yielding $\hat{\lambda}_{(-t)} < 0$ (observed `-0.25` to `-0.60` on Replogle `log1p`).

3. **Source-Only Shared-Axis Projection (`P_perp`)**:
   Leading right singular vector $q_1$ (rank-$1$ PC1 across source perturbation vectors $\{\Delta_{p,s}\}_{s \in \mathcal{S}, p \in \mathcal{P}_s}$) defines $P_\perp = I_G - q_1 q_1^T$.

4. **Candidate & Primary Predictors for Query $p \in \mathcal{P}_t$**:
   - **`Transfer_Mean`**: $\hat{\Delta}^{\mathrm{TM}}_{p,t} = \bar{\Delta}^{\mathrm{src}}_{p,(-t)} = \frac{1}{|\mathcal{S}_p|} \sum_{s \in \mathcal{S}_p} \Delta_{p,s}$ (cross-line mean).
   - **`CIPHER_Linear`**: $\hat{\Delta}^{\mathrm{CIPHER}}_{p,t} = \hat{\lambda}_{(-t)}\, v_t(p)$ (pure control covariance, no source perturbation mean).
   - **`Hybrid_Raw`**: $\hat{\Delta}^{\mathrm{HybRaw}}_{p,t} = \bar{\Delta}^{\mathrm{src}}_{p,(-t)} + \hat{\lambda}_{(-t)}\big(v_t(p) - \bar{v}_{\mathcal{S}_p}(p)\big)$.
   - **`Hybrid_ProjShared` *(Primary Locked)***: $\hat{\Delta}^{\mathrm{HybProj}}_{p,t} = \bar{\Delta}^{\mathrm{src}}_{p,(-t)} + \hat{\lambda}_{(-t)}\, P_\perp\big(v_t(p) - \bar{v}_{\mathcal{S}_p}(p)\big)$.

5. **Pre-Registered Specificity & Negative Controls *(mandatory)***:
   - **`NegCtrl_HybridProj_MismatchedCov`** *(Primary H-A3 Comparator)*: Replace target $C_t$ with another line's $C_{s'}$ ($s' \in \mathcal{S}$), keeping $\bar{\Delta}^{\mathrm{src}}_{p,(-t)}$, $\hat{\lambda}_{(-t)}$, $P_\perp$, and $\bar{v}_{\mathcal{S}_p}(p)$ fixed.
   - **`NegCtrl_HybridProj_ShuffledGenes`**: Randomly permute gene indices of target $C_t$ only ($\Pi C_t \Pi^T$).
   - **`NegCtrl_HybridProj_ShuffledAllCov`** *(A-T3 Audit)*: Permute gene indices on **both** target $C_t$ and source covariances $C_s$ so even minor $-\bar{v}_{\mathcal{S}_p}$ residuals are symmetrized.
   - **`NegCtrl_HybridProj_DiagCov` / Mean-Field**: Zero all off-diagonals of $C_t$ ($C_t \leftarrow \mathrm{diag}(C_t)$).

---

## 5. Exact Evaluation Metrics, Noise Ceiling, & Bootstrap Protocol

All residual definitions match `round9/p2/vcc_linespecific/audit/LIT/LINE_SPECIFIC_METRIC.md`:

### 5.1 Line-Specific Residual Decomposition & Slope
For each target line $t$ and perturbation $p \in \mathcal{P}_t$:
$$r^{\mathrm{obs}}_{p,t} = \Delta^{\mathrm{obs}}_{p,t} - \bar{\Delta}^{\mathrm{src}}_{p,(-t)}, \qquad r^{\mathrm{pred}}_{p,t} = \hat{\Delta}_{p,t} - \bar{\Delta}^{\mathrm{src}}_{p,(-t)}.$$
- **Primary Gene Mask (`resp_top50`)**: Mask out own target gene $g_p$ (if in `COMMON_GENES`), then select the top 50 genes by $|\Delta^{\mathrm{obs}}_{p,g}|$ (or all finite non-target genes if $<50$). Sensitivity mask `all_nontarget` includes all $G-1$ non-target genes.
- **Primary Line-Specific Slope (`ls_slope_ols_top50`)**: Pooled OLS slope with intercept across stacked $(x, y) = (r^{\mathrm{obs}}_{p,g}, r^{\mathrm{pred}}_{p,g})$ in `resp_top50` over $p \in \mathcal{P}_t$: $\hat{\beta}_{\mathrm{OLS}} = \frac{\sum (x-\bar{x})(y-\bar{y})}{\sum (x-\bar{x})^2}$. Companion sensitivity slopes (`slope_no_intercept`, `per_pert_slope_mean`, `all_nontarget`) and Wrap-Up profile cosine agreement (`cos_agreement_top50`) are logged alongside.

### 5.2 Split-Half Noise Ceiling & Joint-Bootstrap Disattenuated Correlation
Using A-T1's disjoint 50/50 split-halves (`half_a`, `half_b`, `half_valid`):
1. Compute split-half line-specific residuals $r^{(A)}_{p,t} = \Delta^{(A)}_{p,t} - \bar{\Delta}^{\mathrm{src}}_{p,(-t)}$ and $r^{(B)}_{p,t} = \Delta^{(B)}_{p,t} - \bar{\Delta}^{\mathrm{src}}_{p,(-t)}$ over the identical `resp_top50` gene mask.
2. Pooled half-split correlation $\rho_{AB} = \mathrm{corr}(r^{(A)}, r^{(B)})$. If $\rho_{AB} \le 0$ on full sample, the line has **zero measurable ceiling** and **H-A2 fails** automatically.
3. For $\rho_{AB} > 0$, Spearman-Brown reliability and correlation noise ceiling:
   $$\rho_{\mathrm{SB}} = \frac{2\rho_{AB}}{1 + \rho_{AB}}, \qquad c_{\mathrm{ceil}} = \sqrt{\rho_{\mathrm{SB}}}, \qquad r_{\mathrm{cal}} = \frac{\mathrm{corr}(r^{\mathrm{obs}}, r^{\mathrm{pred}})}{c_{\mathrm{ceil}}}.$$
4. **Joint Perturbation-Level Bootstrap (`n_boot = 4000`)**: Resample perturbations with replacement ($b=1,\dots,4000$). In every draw $b$, recompute BOTH $r_{\mathrm{resid}}^{(b)}$ AND $\rho_{AB}^{(b)} \to c_{\mathrm{ceil}}^{(b)}$ on the **identical resampled perturbation indices** $I^{(b)}$ so noise-ceiling estimation error is jointly propagated. Require $\ge 95\%$ finite draws (if $>5\%$ of draws give $\rho_{AB}^{(b)} \le 0$, CI fails as `NaN` rather than censoring zero-ceiling draws). Fixed-ceiling CI (`r_resid_ci95 / c_ceil`) is recorded strictly for secondary comparison. *(Analytical control-cell split-half variance ceiling retained solely as fallback if perturbation halves were absent.)*

### 5.3 Overall Perturbation Discrimination & DE50 Correlation (H-A4)
- **`disc_raw`**: Tie-aware mid-rank discrimination percentile $\mathrm{disc}(i) = 1 - \frac{\#\{j:c_{ij}>c_{ii}\} + 0.5\#\{j\ne i: c_{ij}=c_{ii}\}}{|\mathcal{P}_t|-1}$ (comparing both A-T2's own-target-zeroed convention and A-T3's all-panel-targets-masked convention + `cell-eval2` uncentered `pds_cosine_panel`).
- **`de50_r`**: Per-perturbation Pearson correlation on top 50 $|\Delta^{\mathrm{obs}}_{p,g} - \text{ref}_g|$ non-target genes.
- Both metrics are evaluated relative to unperturbed control (`ctrl_*`, $\text{ref}=0$) **and** `Source_Grand_Mean` centroid (`cent_*`, $\text{ref} = \bar{\Delta}^{\mathrm{grand}}_{(-t)}$).

---

## 6. Locked Hypotheses & Hard-Coded Adjudication Rule ($K = 6$)

With $K = 6$ evaluated cell lines, pre-registered thresholds enforce $N_{\mathrm{pass}} = \max(4, \lceil 0.8 K \rceil) = \mathbf{5/6}$ (`83.3%`, strictly $\ge 4/5$) and $N_{\mathrm{inf}} = \lfloor 0.2 K \rfloor = \mathbf{1/6}$ (`16.7%`, strictly $\le 1/5$):

1. **H-A1 (Line-Specific Recovery Above VCC 2026 Ceiling)**:  
   On **$\ge 5/6$ held-out cell lines**, the locked primary predictor **`Hybrid_ProjShared`** achieves `ls_slope_ols_top50` **95% perturbation-bootstrap CI lower bound $> 0.13$**.  
   *(Methodological Lock)*: As proven on Replogle dev, pure `CIPHER_Linear` exhibits large apparent OLS slopes (`+0.64`) indistinguishable from diagonal-only / mismatched covariance nulls (`NegCtrl_CIPHER_DiagCov` `+0.64`) purely from $-\bar{\Delta}^{\mathrm{src}}$ magnitude mismatch, and unprojected `Hybrid_Raw` exhibits $-\bar{v}_{\mathcal{S}_p}$ shared-axis leakage (`NegCtrl_Hybrid_DiagCov` `+0.021`). Therefore H-A1 is adjudicated strictly on **`Hybrid_ProjShared`**; `CIPHER_Linear` and `Hybrid_Raw` cannot substitute for a failed `Hybrid_ProjShared`.
2. **H-A2 (Noise-Calibrated Line-Specific Correlation)**:  
   On **$\ge 5/6$ held-out cell lines**, full-sample $\rho_{AB} > 0$, finite joint-bootstrap draws $\ge 95\%$, and noise-calibrated residual correlation $r_{\mathrm{cal}}$ has **joint-bootstrap 95% CI lower bound $> 0$**.
3. **H-A3 (Target-Covariance Specificity)**:  
   On **$\ge 5/6$ held-out cell lines**, matched target covariance (`Hybrid_ProjShared`) achieves strictly higher `ls_slope_ols_top50` than mismatched covariance (`NegCtrl_HybridProj_MismatchedCov`) with **non-overlapping 95% CIs** (`mismatch_ci95_hi < matched_ci95_lo`).
4. **H-A4 (Non-Inferiority on Overall Perturbation Accuracy)**:  
   Relative to `Transfer_Mean`, `Hybrid_ProjShared` is non-inferior on **BOTH `disc_raw` and `de50_r`** (satisfied across both `ctrl_*` and `cent_*` conventions): at most **$\le 1/6$ held-out lines** exhibit paired bootstrap 95% CI strictly below zero ($\text{diff\_ci95\_hi} < 0$ for `Hybrid_ProjShared - Transfer_Mean`).

### Hard-Coded Verdict Rule (Zero Flexibility)
- **`CONFIRMED`**: Returned **only** if **ALL FOUR** criteria (**H-A1 AND H-A2 AND H-A3 AND H-A4**) pass on the locked primary configuration.
- **`未确认` (NOT CONFIRMED)**: If **ANY** of **H-A1, H-A2, H-A3, or H-A4** fails its threshold—including if signal correlation (`H-A2`) and specificity (`H-A3`) pass while `ls_slope_ols_top50` remains below `0.13` (`H-A1`) or `disc_raw`/`de50_r` degrades (`H-A4`)—the primary study conclusion **MUST state `未确认`** explicitly in all summaries and audit reports. Post-hoc threshold relaxation, model swapping, metric substitution, or line removal is prohibited; diagnostic analyses (`S1–S4`) are confined to a secondary section.

---

## 7. Prospective VCC 2026 Generalization Protocol (`AC16`, `MOLM13`, `U-87MG`)

Aligned with `round9/p2/vcc_linespecific/src/prospective/{HUMAN_STEPS.md,CELL_EVAL2_REPRO.md,vcc2026_format.py,prospective_predict.py}` and `P2_RULES.md` Rule 6:
1. **Compliance & Data Access Boundary**: Agents do **not** create Challenge accounts, complete KYC verification, download gated `vcc_2026_controls_final.zip`, or push to external repositories. Human PI (`Zijian Zeng`) either obtains written permission from `help@virtualcellchallenge.org` to run offline inference on `controls-final` post-deadline or awaits CC0 1.0 release on `gs://arc-institute-virtual-cell-atlas/virtual-cell-challenge/2026/`.
2. **Pre-GT Prediction Generation & Public Hash Timestamping**: Once compliant unperturbed `controls-final` inputs are obtained **before** test ground truth becomes public, `prospective_predict.py` runs the exact locked primary predictor (`log1p`, `shrink=0.0`, `mask_self=True`, `Hybrid_ProjShared`) fitted strictly on the 6-line source reference and target `non-targeting` control matrices (`18,533` gene panel, `400` cells/pert/context, dispersed CSR integer counts validated via `vcc2026_format.py`). Prior to public GT release, only `SHA256SUMS.txt`, `prediction_manifest.json`, `TIMESTAMP_DECLARATION.md`, and code (never restricted input `.h5ad`/CSV files) are committed to public GitHub (+ optional OpenTimestamps `.ots`).
3. **Post-GT Prospective Adjudication**: Once VCC 2026 GT is publicly released under CC0 1.0, score unmodified hashed predictions using:
   - Official **`ArcInstitute/cell-eval2 v0.22.0`** (`--preset vcc2026`: 6 scored metrics `de_wilcoxon_lfc_nmae`, `expr_mse_unbiased_capped_norm`, `de_wilcoxon_direction_fidelity_yield_raw`, `de_wilcoxon_direction_reach_raw`, `de_wilcoxon_sig_jaccard`, `pds_cosine`, and `avg_score` vs `prep-real-bundle`).
   - Exact same **H-A1..H-A4** line-specific slope, noise-calibrated correlation, mismatched-covariance control (`AC16/MOLM13/U-87MG`), and non-inferiority adjudication defined in §§5–6, reporting negative prospective verdicts verbatim as `未确认`.

---

## 8. Independent Statistical Audit Plan (`A-T3`)

Upon `round9/p2/status/PREREG_COMMIT.txt` creation and `A-T2` completion of `round9/p2/vcc_linespecific/results/`:
1. Execute standalone `round9/p2/vcc_linespecific/audit/audit_runner.py` (**zero imports from `src/`**) under locked primary settings (`--space log1p --shrink 0.0 --cov-max-cells 0 --project-shared-axis --status-dir round9/p2/status`) and run sensitivity checks (`--all-targets` S1, raw unprojected S2, `--shrink 1e-3` S3).
2. Verify zero-leakage invariance (corrupting held-out $\Delta_t$ changes predicted $\hat{\Delta}_t$ by `0.0`), run `audit/cross_check.py` against `A-T2`'s `line_metrics`, `bootstrap_ci`, and `adjudication` outputs (`tol = 1e-4`), and document all matches and discrepancies in `round9/p2/vcc_linespecific/audit/AUDIT_REPORT.md`.

---

## 9. API-Verified Bibliography

All DOIs/arXiv IDs verified via Crossref, bioRxiv, or arXiv APIs on 2026-11-21:
1. **Kuznets-Speck B et al. (19 authors)** (2025/2026). *Co-fluctuations of genes in single cells predict transcriptome-wide outcomes to perturbations* (`CIPHER`). bioRxiv `doi:10.1101/2025.06.27.661814` (v1: 2025-07-01; v2: 2026-09-03), PMID `40631127`.
2. **Replogle JM et al. (17 authors)** (2022). *Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq*. *Cell*, `doi:10.1016/j.cell.2022.05.013`, PMID `35688146`.
3. **Nadig A et al. (8 authors)** (2025). *Transcriptome-wide analysis of differential expression in perturbation atlases*. *Nature Genetics*, `doi:10.1038/s41588-025-02169-3`, PMID `40259084`.
4. **Nourreddine S et al. (22 authors)** (2026). *A genome-scale CRISPRi perturbation atlas of human induced pluripotent stem cells*. *Nature Biotechnology*, `doi:10.1038/s41587-026-03199-w`, PMID `42386990`.
5. **Miller HE et al. (6 authors)** (2026). *Deep learning perturbation models can outperform baselines on calibrated metrics*. *Nature Biotechnology*, `doi:10.1038/s41587-026-03307-w`.
6. **Adduri A et al. (81 authors)** (2026). *Virtual Cell Challenge 2026: Benchmarking zero-shot generalization across cellular contexts*. *Cell*, `doi:10.1016/j.cell.2026.08.004`.
7. **Adduri AK et al. (32 authors)** (2026). *Predicting cellular responses to perturbation across diverse contexts with State*. *Cell*, `doi:10.1016/j.cell.2026.07.052`.
8. **Viñas Torné R et al. (8 authors)** (2025). *Systema: a framework for evaluating genetic perturbation response prediction beyond systematic variation*. *Nature Biotechnology*, `doi:10.1038/s41587-025-02777-8`.
9. **Kawsar M et al. (6 authors)** (2026). *Transferring CRISPRi knockdown responses to new cell lines requires screening part of each line* (`CASTER`). bioRxiv `doi:10.64898/2026.11.15.766815` (2026-11-17).
10. **Huang Y et al. (4 authors)** (2026). *The recoverable resolution of cellular perturbation-response prediction* (`CGC`). bioRxiv `doi:10.64898/2026.09.30.755756` (2026-10-01).
11. **Dhaka B et al. (3 authors)** (2026). *HyperMap: An Efficient Framework for Transferring Perturbation Responses Across Diverse Biological Contexts*. bioRxiv `doi:10.64898/2026.04.23.720505` (2026-04-27).
12. **Luan et al. (4 authors)** (2026). *Centered cross-context perturbation transfer on the public 2025 Virtual Cell Challenge benchmark*. bioRxiv `doi:10.64898/2026.10.20.760514` (2026-10-25).
13. **Molina et al. (2 authors)** (2026). *Perturbation response decomposition enables biologically aligned generalization to unseen perturbations and cellular contexts*. bioRxiv `doi:10.64898/2026.07.24.740459` (2026-07-27).
14. **Gao J** (2026). *The Delta Gap: Cellular World Models for Perturbation Response Prediction Are Far Harder Than Global Correlation Implies*. `arXiv:2611.10658v1` (2026-11-07).
15. **Zhang Z, Huang Y, Feng F et al. (4 authors)** (2026). *A low-rank decoder bottleneck bounds reliability in foundation-model perturbation prediction*. `arXiv:2611.05849v1` (2026-11-04).
16. **Arc Institute (`cell-eval2` & VCC 2026 Wrap-Up)**: `ArcInstitute/cell-eval2` `v0.22.0` (`--preset vcc2026`, verified directly via repository tests) and Arc VCC 2026 Wrap-Up (`2026-11-19`, top-10 line-specific slope `0.02–0.13`, replicate `1.02`). *Unindexed DOI Disclosure*: Dobin metrics preprint (`10.64898/2026.11.17.767214`) returned HTTP 404 / 0 results across Crossref, bioRxiv, OpenAlex, EuropePMC, and PubMed APIs on 2026-11-21 (`www.biorxiv.org` HTML Cloudflare-gated, not bypassed); exact locked equations are explicitly fixed in §§4–6 above.
