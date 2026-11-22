# Universal Target-Free Null-Cone Cure for Zero-Shot Perturbation Benchmarks (`exploratory/null_cone_theorem/`)

## 1. Scientific Overview & Relationship to Rounds 9–11

This self-contained methodological module builds directly on the lab's three commit-locked, pre-registered zero-shot Perturb-seq audits across 13 human cell lines (`10,334,369` QC single cells; Rounds 9–11). **All Round 9 (`PREREG_A v4`), Round 10 (`PREREG_C`), and Round 11 (`PREREG_D`) locked confirmatory verdicts remain strictly unchanged and untouched**; every discovery documented here is explicitly demarcated as a **mechanistic / post-hoc methodological contribution** derived from the unblinded 13-line corpus.

### Four Verified Contributions (Zero Unchecked Claims)

1. **Theorem 1 — Deterministic $w(\mu)$ Signal-Tail Conditioning Bias & Why Out-of-Sample Split-Half DEG Selection Fails (`12/12` Human Cell Lines):**
   - Nicol, Shivakumar & Irizarry (2026, `doi:10.64898/2026.05.07.723486`) demonstrated that in-sample top-$k$ DEG selection amplifies finite-sample shared control noise (`ctl-mean`) and proposed split-half DEG selection as the cure.
   - Across **all 12 eval-valid human Perturb-seq cell lines** (`10.33M` single cells) — including disjoint human donors (`cd4t_rest`, `cd4t_stim48`) — we prove mathematically (under zero sampling noise) and empirically (`n_boot = 4,000`, `n_perm = 100`) that **out-of-sample split-half top-50 selection (`SA-8 resp_top50_split`) fails to eliminate false-positive residual scores** for baseline-free transforms (`E1`, `E6`) and gene-permuted nulls (`Perm100`), and often **amplifies** them relative to in-sample `resp_top50`.
   - **Mechanism:** `log1p(CP10k)` applies a monotone gene-specific Jacobian $w(\mu_{c,g}) = \mu_{c,g}/(1+\mu_{c,g})$. Conditioning on large target effects $|\Delta_{t,g}| \approx w(\mu_{t,g})|\delta_g|$ deterministically enriches genes near $w_t \to 1$, compressing $\mathrm{Var}_M(w_t)$ while source weights $\bar{w}_S$ regress toward their mean, forcing a nonzero $\delta^2$-weighted intercept $a_{M,\delta^2} > 0$ that survives at infinite sample size.

2. **Theorem 2 — Universal Target-Free Null-Cone Semi-Partial Metric (`N_5`):**
   - Every constant-baseline predictor (`E1`, `E2(P10/P90)`), asymptotic fold-change limit (`E6 = \log\overline{\mathrm{FC}}`, `E7 = \overline{\mathrm{FC}}-1`), and expected gene-permutation null ($\mathbb{E}_\pi[\mathrm{FCT}(\mu_\pi)]$) lies inside the convex cone $\mathcal{C}(\overline{\mathrm{FC}}_p) = \overline{\mathrm{conv}}\{g_c(\overline{\mathrm{FC}}_p) = \log(1 + c\,\overline{\mathrm{FC}}_p) - \log(1 + c) : c > 0\}$.
   - Across 6 orders of magnitude ($c \in [10^{-3}, 10^3]$), 3 fixed universal **target-free** nodes $\{g_{0.05}, g_{1.0}, g_{20.0}\}$ capture **`99.987%`** of cone variance (`min R^2 = 0.99919`).
   - Defining the 5D **100% target-free** null subspace
     $$\mathcal{N}_5(p) = \mathrm{span}\Big\{\mathbf{1},\; \Delta_{\mathrm{TMean},p},\; b_S,\; g_{0.05}(\overline{\mathrm{FC}}_p),\; g_{1.0}(\overline{\mathrm{FC}}_p),\; g_{20.0}(\overline{\mathrm{FC}}_p)\Big\}$$
     and evaluating the **Noise-Calibrated Semi-Partial Null-Cone Correlation**
     $$r_{\mathrm{cone\text{-}sp}}(\Delta_{\mathrm{obs}}, \hat{\Delta}) = \frac{1}{c_{\mathcal{N}_5}} \frac{\langle P_{\mathcal{N}_5}^\perp r_{\mathrm{obs}},\; P_{\mathcal{N}_5}^\perp r_{\mathrm{pred}} \rangle}{\|P_{\mathcal{N}_5}^\perp r_{\mathrm{obs}}\|_2 \cdot \|(r_{\mathrm{pred}})_{\mathrm{centered}}\|_2}$$
     simultaneously annihilates all four benchmark hallucinations (`GlobalScale`, `Source_Grand_Mean` $b_S$, `E1/E2/E6/E7`, and `Perm100`) across all 4 masks (`all_nontarget`, `resp_top50`, `SA-8 resp_top50_split`, `top50_src`) while preserving orthogonal energy fraction $\eta \approx 0.75\text{–}0.996$ for genuinely target-specific models (`FCT`, `TMVS`, `Hybrid_ProjShared`, true biological rewiring).

3. **Direct Q5a Cross-Line Conservation Test (`33` Human Cell-Line Pairs; Zero Prior Art per `PRIOR_ART_D.md`):**
   - Uncalibrated raw `logFC` loses to `log1p` (`0/33` on `all` and `expr01`) due to low-expression shot noise.
   - Once disattenuated by split-half noise ceilings on expressed genes ($\mu \ge 0.5$), true conservation reverses (`logFC > log1p` on `25/33` pairs, paired Wilcoxon `p = 1.47e-3`).
   - Combining conserved `logFC` with target-baseline `log1p` compression (`FCT`) beats direct `log1p` transfer on `28/33` (`all`), `32/33` (`expr01`), and **`33/33`** (`expr05`, Wilcoxon `p = 2.33e-10`) human cell-line pairs.

4. **Honest Effect-Size & Scope Boundary:**
   - Although target-baseline `log1p` rescaling (`FCT`) carries genuine, reproducible gene-specific cross-line signal after removing all four benchmark artifacts under `N_5`, its calibrated residual correlation remains modest (`r_cone_sp ~ 0.05–0.20`, $R^2 \approx 0.3\%\text{–}4\%$), confirming that `>95%` of cell-type-specific transcriptional rewiring remains an open biological challenge.

---

## 2. Reproducibility Commands

All scripts use `python3 -I` and read standardized caches in `round9/p2/shared/lines/`, `round10/shared/lines/`, and `round11/shared/lines/` strictly read-only:

```bash
# 1. Run 5-test unit & adversarial verification suite (0.08s)
python3 -I exploratory/null_cone_theorem/test_null_cone.py -v

# 2. Run full 13-cell-line audit (N1/N4/N5, 100 perms, 4,000 joint bootstraps, Q5a, noiseless audit)
python3 -I exploratory/null_cone_theorem/run_13line_audit.py

# 3. Generate verified Markdown report directly from JSON (zero manual transcription)
python3 -I exploratory/null_cone_theorem/generate_report.py

# 4. Generate 4 publication-grade PNG + vector PDF figures directly from JSON
python3 -I exploratory/null_cone_theorem/generate_figures.py
```

## 3. Freeze & Governance Integrity

- Zero files modified in `src/xatlas_panhuman_solver.py` (`3609811`), `round9/`, `round10/`, or `round11/`.
- Zero X-Atlas/Orion data accessed; zero expression matrices committed to git.
