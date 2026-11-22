# Pre-Registration Deviations Log (`round9/p2/vcc_linespecific`)

## D1. Shared-Subspace Projector `P_perp`: Rank-1 `q_1` (`PREREG_A v4 §4.2`) vs. Rank-3 `span{1, b_S, a_S}` (`CODE_FREEZE 1ca9f60f`)

- **Discovery Stage:** Identified during Master Controller review of the Round 9 preprint draft (`round9/p2/status/PAPER_REVIEW_R1.md`, Blocking Item `B1`, `2026-11-22`).
- **Pre-Freeze Origin:** Present prior to `CODE_FREEZE` (`1ca9f60f49adf1ac7cce8266e40533273411db30`, `PREREG_A.md` commit `4cc89ebae67f0ed0748b80f994d26dcc9c2b95ef`).
- **Exact Discrepancy:**
  1. **Pre-registration prose (`PREREG_A.md` v4 §4.2):** Wrote the rank-1 orthogonal projector $P_\perp = I_G - q_1 q_1^\top$, where $q_1$ denotes the first principal component of source-line perturbation vectors.
  2. **Frozen primary codebase (`1ca9f60f`, `src/predictors.py::fit_shared_subspace`) & independent audit codebase (`audit/independent_predictors.py::fit_source_shared_axis_basis`):** Both removed the **3-dimensional orthonormal subspace** $\operatorname{span}\{\mathbf{1}, \bar{b}_S, a_S\}$, where $\bar{b}_S$ is the source perturbation centroid and $a_S$ is the top right singular vector of row-normalized source residuals orthogonal to $\operatorname{span}\{\mathbf{1}, \bar{b}_S\}$ (matching `src/xatlas_panhuman_solver.py::fit_baseline_and_axis`).
  3. **Code independence note:** Because `A-T3`'s independent audit implementation (`audit/independent_predictors.py`) replicated the same rank-3 $\operatorname{span}\{\mathbf{1}, \bar{b}_S, a_S\}$ projection as `src/predictors.py`, the two implementations were not independent with respect to projector dimension in `1ca9f60f`.
- **Affected Predictors & Controls:**
  - `Hybrid_ProjShared` (Group A Primary),
  - `Hybrid_on_TMVS` (`R14` secondary increment),
  - `NegCtrl_HybridProj_MismatchedCov(s')` (all 5 source lines),
  - `NegCtrl_HybridProj_SrcCorrTgtVar`, `NegCtrl_HybridProj_ShuffledGenes`, `NegCtrl_HybridProj_DiagCov`.
  - *(Unprojected `Hybrid_Raw` (`S2`), `Transfer_Mean`, `CIPHER_Linear`, `Transfer_Mean_VarScaled` (Group B), `Transfer_Mean_GlobalScale` (`R18`), and `FCT` (`R15`) do not use $P_\perp$ and are strictly invariant.)*
- **Governance & Post-Freeze Sensitivity Switch (`S2b`):**
  - Primary confirmatory adjudication remains anchored to the frozen `1ca9f60f` rank-3 implementation (`proj_mode="span3"` kept as default, bit-identical to `1ca9f60f`).
  - Added an **opt-in, non-default** switch (`--proj-mode rank1_q1` in `src/run_eval.py` / `proj_mode="rank1_q1"` in `src/predictors.py`) to execute sensitivity **`S2b`** matching literal `PREREG_A v4 §4.2` rank-1 centered PCA PC1 (`q_1 = Vt[0]` of `Y - mean(Y, axis=0)`, $P_\perp = I_G - q_1 q_1^\top$) on held-out (`results/*_holdout_log1p_sens_s2b_rank1q1.*`).
- **Observed `S2b` (`rank1_q1`) Result Summary (`results/*_holdout_log1p_sens_s2b_rank1q1.*`, 597.37s):**
  - **Confirmatory Verdicts Unchanged:** `Group A (Hybrid_ProjShared): 未确认 | Group B (Transfer_Mean_VarScaled): 未确认`. Zero individual Group A hypotheses changed verdict:
    - `H-A1` (`slope 97.5% CI lo > 0.13`): **未确认** (`0/3` pass: `nadig_hepg2 +0.0089 [+0.0054, +0.0127]`, `nadig_jurkat +0.0102 [+0.0070, +0.0136]`, `vcc2025_h1 -0.0064 [-0.0216, +0.0063]`).
    - `H-A2` (`r_orth_cal 97.5% CI lo > 0`): **未确认** (`2/3` pass: `nadig_hepg2 +0.0468 [+0.0211, +0.0724]`, `nadig_jurkat +0.1040 [+0.0783, +0.1296]`, `vcc2025_h1 -0.0299 [-0.1118, +0.0445]` FAIL).
    - `H-A3` (`matched > ALL 5 mismatched cov on BOTH slope AND r_orth_cal`): **未确认** (`0/3` pass). Notably, under `rank1_q1`, `nadig_hepg2` `NegCtrl_HybridProj_DiagCov` (`+0.0771`, slope `+0.0141`) also exceeded `Hybrid_ProjShared` (`+0.0468`, slope `+0.0089`), whereas frozen `span3` suppressed this marginal-variance bleed-through.
    - `H-A4` (`strict AND non-inferiority across 6 conventions`): **未确认** (`3/3` inferior; `ctrl_disc_raw=3/3`, `ctrl_de50_r=2/3`, `ctrl_disc_raw_panel=3/3`, `cent_disc_raw=2/3`, `cent_de50_r=1/3`, `cent_disc_raw_panel=2/3`).
  - **`R14` Co-Fluctuation Increment (`Hybrid_on_TMVS` - `TMVS` on `r_orth_cal`) Unchanged:** Fails unanimity at `2/3` (`confirmed=False`; `95%: 2/3`, `97.5%: 2/3`), retaining allowed statement `"在本数据上，共波动结构没有超出边际方差的可检测增量"` (`nadig_hepg2` `+0.0126` `ci975_lo=+0.0014`, `nadig_jurkat` `+0.0375` `ci975_lo=+0.0238`, `vcc2025_h1` `+0.0026` `ci95_lo=-0.0611`, `ci975_lo=-0.0711` FAIL).
