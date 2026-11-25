# Round 16 Proposal (`PROPOSAL_R16_FOUNDATION_AUDIT.md`): Auditing Open-Weight Virtual Cell Models & Double-Orthogonalized Specificity

- **Author:** Main Line (Antigravity) — **PROPOSAL ONLY (Awaiting Master Controller & Human Project Lead Approval per `controller/README.md`)**
- **Date:** `2026-11-24` (`10:42 +08:00`)
- **Governance Compliance:** Zero modifications to Rounds 9–15 frozen files; zero evaluations run on any held-out dataset; zero model weights downloaded; zero edits to `paper/` or `README*`; zero media/outreach text.

---

## 1. Scientific Gap Keeping `MANUSCRIPT_DRAFT_v5.md` Below *Nature Methods / Nature Biotechnology* Bar

As of `2026-11-24`:
1. **`Round 15` (`PREREG_H`, frozen `3dc383b`) is blocked on external release:** `storage.googleapis.com/storage/v1/b/arc-institute-virtual-cell-atlas/o?delimiter=/&prefix=virtual-cell-challenge%2F2026%2F` returned `0` items (`2026-11-24 10:37 +08:00`).
2. **Closed-form baselines only:** Rounds 9–14 audit closed-form baselines (`FCT`, `TMVS`, `CIPHER`, `Source_Grand_Mean`) and do not test whether published deep-learning Virtual Cell models (`Arc STATE`) or VCC 2026 winning algorithms (`Illumina CellularAI`) actually pass matched-vs-mismatched target-control specificity.
3. **Single confirmatory line (`K_conf = 1`):** `a549_ko` (`GSE345058`) is the sole confirmatory line after three pre-registered gate failures (`HEK293` 16, `GSE327644` macrophages 3, `GSE337804` Day-16 A549 32).
4. **Two unresolved methodological questions left open in `MANUSCRIPT_DRAFT_v5.md`:**
   - **§2.3 (`Source_Grand_Mean` `+0.4051` > `FCT` `+0.1898`):** Single-axis orthogonalization (`r_orth_cal` against $\bar{\Delta}^{\mathrm{src}}_p$ alone) leaves the cross-perturbation shared axis $\bar{b}_S$ un-removed, allowing a target-agnostic, perturbation-agnostic constant vector to outscore every target-conditioned predictor.
   - **Limitation 5 (`Claim 2` mechanism marked "unknown"):** Round 11's naive noiseless simulation predicted the wrong sign on `8/12` lines under `all_nontarget`, leaving the cause of `E1_ConstMed` jumping from `-0.020` (`all_nontarget`) to `+0.257` (`resp_top50_split`) unexplained.

---

## 2. Verified Local-SSD Feasibility (Zero-Download Metadata & Hardware Audit)

All checks executed on local SSD (`/Users/mac/AI-Biology-Lab`, `2026-11-24`):
- **Local Compute:** Apple M3 Pro, **36.0 GB unified memory**, **2.4 TB free local SSD**, `/opt/anaconda3/bin/python3` (`Python 3.13.5`) with **`torch=2.15.0` (`MPS=True`)**, `numpy=2.2.3`, `scipy=1.15.3`, `sklearn=1.6.1`, `h5py=3.12.1`, `transformers=5.22.2`, `safetensors=0.8.0`.
- **Open-Weight Models Verified (Public, Ungated, Zero Login):**
  1. **Arc Institute `STATE` (`ST-SE-Replogle` & `ST-HVG-Replogle` + `SE-600M`):**
     - Source: `https://huggingface.co/arcinstitute/ST-SE-Replogle`, `https://huggingface.co/arcinstitute/ST-HVG-Replogle`, `https://huggingface.co/arcinstitute/SE-600M` (code: `https://github.com/ArcInstitute/state`).
     - License: **Arc Research Institute State Model Non-Commercial License** (`2025-06-23`, permits non-commercial academic research use, execution, outputs, and dissemination).
     - Canonical **Zero-Shot Checkpoints** (~525 MB each; `SE-600M` `model.safetensors` = 2.73 GB):
       - `zeroshot/hepg2/checkpoints/final.ckpt` (trained on `jurkat`, `k562`, `rpe1`; `hepg2` held out)
       - `zeroshot/jurkat/checkpoints/final.ckpt` (trained on `hepg2`, `k562`, `rpe1`; `jurkat` held out)
       - `zeroshot/k562/checkpoints/final.ckpt` (trained on `hepg2`, `jurkat`, `rpe1`; `k562` held out)
       - `zeroshot/rpe1/checkpoints/final.ckpt` (trained on `hepg2`, `jurkat`, `k562`; `rpe1` held out)
     - **Zero-Shot Line Eligibility:**
       - Within Replogle-4: **4 native zero-shot held-out lines** (`nadig_hepg2`, `nadig_jurkat`, `replogle_k562_essential`, `replogle_rpe1`), each evaluated strictly with its corresponding leave-one-line-out `zeroshot/<line>/final.ckpt`!
       - Cross-study unseen lineages (strictly absent from all Replogle-4 training sets): `a549_ko` (`N=445`), `dld1` (`N=86`), `cd4t_rest` (`N=302`), `cd4t_stim48` (`N=503`), `vcc2025_h1` (`N=80`), `feng_ipsc` (`N=93`) — up to **`K_conf = 5–9`** independent cell contexts, all already harmonized on local SSD! *(Note: exact prospective vs. retrospective status of already-seen R9–R14 lines must be disclosed transparently per E6, while model predictions on them are strictly out-of-sample for the frozen `STATE` weights!)*
  2. **Illumina `CellularAI` (1st Place, VCC 2026):**
     - Source: `https://github.com/Illumina/CellularAI` (**Apache-2.0 License**).
     - Open code implements the complete 5-stage pipeline (`pool effects -> plan responses -> residual model -> design counts -> generate cells`). Its proprietary Illumina Billion Cell Atlas summary bundle is gated (`EXCLUDED_ACCESS`), **but the exact Apache-2.0 algorithm can be run by building the summary catalog strictly from our public `r9_6` source atlas** (`pool effects` weighted by control similarity + control-conditioned residual model), enabling a clean, license-compliant audit of VCC 2026's winning architecture under matched vs. mismatched controls!
  3. **Xaira Therapeutics `X-Cell` (3rd Place, VCC 2026):**
     - Verified `https://huggingface.co/Xaira-Therapeutics/X-Cell` (`cc-by-nc-sa-4.0`): README explicitly states *"Status: Model weights and inference code coming soon"* (`0` weight files). Classify `UNRELEASED` and exclude unless released before freeze.

---

## 3. Proposed Round 16 Scope (Two Coupled Pillars — Pre-Registered as `PREREG_I`)

### Pillar A: Matched-vs-Mismatched Control Audit of Open-Weight Virtual Cell Models (`Arc STATE` + Open `CellularAI`)
For each qualifying zero-shot cell line $t$ and each frozen model $M \in \{\text{STATE-SE}, \text{STATE-HVG}, \text{CellularAI-r9\_6}, \text{FCT}, \text{CIPHER\_Linear}\}$:
1. Run model $M$ with:
   - **Matched target controls** ($X^{\mathrm{ctrl}}_t$),
   - **Mismatched cell-line controls** ($X^{\mathrm{ctrl}}_{s'}$ for each source/other cell line $s' \neq t$), and
   - **Gene-permuted target controls** (when gene-space input applies, e.g. `ST-HVG`, `CellularAI`, `FCT`).
2. Score under **both**:
   - **Legacy `r_orth_cal`** (for full backward comparability with Rounds 9–15), and
   - **Double-orthogonalized `r_dorth_cal`** (Pillar B below, removing $\operatorname{span}\{\mathbf{1}, \bar{\Delta}^{\mathrm{src}}_p, \bar{b}_S\}$),
   across **both `all_nontarget` and `resp_top50_split`**, with paired 97.5% joint perturbation bootstraps (`n_boot = 4,000`, fresh RNG seed).
3. **Why this is top-journal decisive regardless of outcome:**
   - If `STATE` / `CellularAI` **fail** matched-vs-mismatched control specificity under `r_dorth_cal` (or collapse once `Source_Grand_Mean` / `GlobalScale` are projected out), it directly explains Arc's `0.02–0.13` line-specific slope finding and exposes a fundamental limitation of current virtual cell foundation models.
   - If `STATE` / `CellularAI` **pass** matched-vs-mismatched specificity AND exceed `FCT`, it provides pre-registered quantitative evidence of context-specific perturbation transfer beyond closed-form baseline scaling (novelty strictly contingent on `A-T5` `PRIOR_ART_I.md`).

### Pillar B: Closing Both Theoretical Gaps in `MANUSCRIPT_DRAFT_v5.md`
1. **Double-Orthogonalized Specificity Metric (`r_dorth_cal`, Resolving §2.3):**
   - Project observed residual $y_{t,p} = \Delta^{\mathrm{obs}}_{t,p} - \bar{\Delta}^{\mathrm{src}}_p$ and predicted residual $\hat{y}_{t,p} = \hat{\Delta}_{t,p} - \bar{\Delta}^{\mathrm{src}}_p$ onto $\operatorname{span}\{\mathbf{1}, \bar{\Delta}^{\mathrm{src}}_p, \bar{b}_S\}^\perp$.
   - **Pure-synthetic verification (`2026-11-24`, zero real data):**
     - `r_orth(all)`: `GlobalScale = 0.000000`, **`Source_Grand_Mean = +0.342431`**, `0.6*TM + 0.4*SGM = +0.342431`, `FCT = +0.226332`.
     - `r_dorth(all)`: `GlobalScale = 0.000000`, **`Source_Grand_Mean = +0.000000`**, **`0.6*TM + 0.4*SGM = +0.000000`**, `FCT = +0.225161`.
   - Theorem guarantees exact zero score for **every** target-agnostic affine combination $\alpha_0 \mathbf{1} + \alpha_1 \bar{\Delta}^{\mathrm{src}}_p + \alpha_2 \bar{b}_S$, eliminating both amplitude scaling (`R18`) and global stress-axis (`Source_Grand_Mean`) confounding simultaneously.
2. **Exact Mechanism of Claim 2 (`SA-8` Biological Winner's Curse Under Split-Half Selection, Resolving Limitation 5):**
   - Cell-level split-half (`Half A` selects top-50 $|\Delta^A_{t,p,g}|$, `Half B` scores) removes only independent measurement noise ($\epsilon^A \perp \epsilon^B$), **not** conditioning on true biological amplitude $|\Delta^{\mathrm{true}}_{t,p,g}|$.
   - Because $\Delta_{t,p,g} \approx w(\mu_{t,g})\delta_{p,g}$ with $w(\mu) = \mu / (1 + \mu)$, conditioning on top-50 $|\Delta^A_{t,p,g}|$ simultaneously (i) excludes silent/low-$\mu_{t,g}$ genes where noisy `FC_bar` estimates cancel `E1_ConstMed` on `all_nontarget`, and (ii) selects genes with disproportionately high target expression relative to source ($\mathbb{E}[w(\mu_t) - w(\mu_s) \mid \text{top-50}] \gg 0$; in pure-synthetic check, `w(mu_t)` rose `0.2207 -> 0.5404` and `w_t - w_s` rose `+0.0018 -> +0.0692`, driving an `E1` jump of `+0.197` to `+0.411` depending on source noise & cross-line $\mu$ divergence).
   - Proposed deliverable: formal proposition + exact cross-regime synthetic verification (fixing R11's oversimplified simulation) + pre-registered empirical diagnostic tests across `STATE`, `CellularAI`, and `FCT`.

---

## 4. Requested Decision From Master Controller & Human Project Lead

1. **Approve `Round 16` (`PREREG_I`)** covering **Pillar A (`Arc STATE` + Open `CellularAI` matched-vs-mismatched audit) + Pillar B (`r_dorth_cal` + Claim 2 analytical mechanism)** on local M3 Pro SSD?
2. **Designate blinding/protocol structure:**
   - Option 1 (Strictest): Treat already-observed R9–R14 lines as **retrospective model-audit cohorts** (explicitly labeled as pre-existing target lines, with frozen `STATE` weights never trained on unseen lines and fresh `PREREG_I` hypotheses locked before any model inference!), AND simultaneously add `STATE` + `CellularAI` + `r_dorth_cal` as an **independent pre-release companion addendum** for Arc VCC 2026 (without touching frozen `PREREG_H`!).
   - Option 2: Also screen one genuinely untouched public CRISPR dataset (not used in R9–R14 and not in Replogle-4) as a brand-new prospective `PREREG_I` holdout for `STATE` + `CellularAI` + `FCT`.
3. **Human 1-minute item (unchanged):** Save `Nicol et al. 2026 v2` and `Heidari et al. 2026` PDFs locally via Chrome (`Zheng` profile) so prior-art verification for `MANUSCRIPT_DRAFT_v5.md` can close cleanly.

*Standing by for Master Controller directive before downloading any model checkpoints or running any code on real data.*
