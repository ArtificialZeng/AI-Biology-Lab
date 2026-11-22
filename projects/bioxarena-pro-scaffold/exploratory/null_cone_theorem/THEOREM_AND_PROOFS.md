# Mathematical Appendix: Exact Proofs of Theorem 1, Theorem 2, and Proposition 3

> **Scope & Status:** Formal derivations supporting `exploratory/null_cone_theorem/`. All empirical constants and cross-cell-line figures are verified programmatically by `test_null_cone.py` (5/5 unit & adversarial checks) and `run_13line_audit.py` across 13 human Perturb-seq atlases (`10,334,369` QC single cells).

---

## 1. Setup and Nonlinear `log1p` Geometry

For cell line $c \in \mathcal{S} \cup \{t\}$, gene $g \in \{1,\dots,G\}$, and perturbation $p \in \{1,\dots,N\}$, let $\mu_{c,g} \ge 0$ denote the unperturbed library-normalized (`CP10k`) mean expression and $\mathrm{FC}_{c,p,g} = \exp(\delta_{c,p,g}) > 0$ denote the multiplicative expression ratio (perturbed over control). In standard `log1p(CP10k)` space, the perturbation response is
$$\Delta_{c,p,g} = \log\!\big(1 + \mu_{c,g}\exp(\delta_{c,p,g})\big) - \log(1 + \mu_{c,g}) \triangleq g_{\mu_{c,g}}\!\big(\mathrm{FC}_{c,p,g}\big).$$

**Lemma 1 (Exact Derivatives and Concave Saturation).** For every $\mu \ge 0$ and $\delta \in \mathbb{R}$,
$$\frac{\partial}{\partial \delta} g_\mu(e^\delta)\Big|_{\delta=0} = w(\mu) \triangleq \frac{\mu}{1+\mu} \in [0, 1), \qquad \frac{\partial^2}{\partial \delta^2} g_\mu(e^\delta)\Big|_{\delta=0} = w(\mu)\big(1 - w(\mu)\big) \in \left[0, \frac{1}{4}\right].$$
Hence the Taylor expansion around $\delta = 0$ is
$$\Delta_{c,p,g} = w(\mu_{c,g})\,\delta_{c,p,g} + \kappa(\mu_{c,g})\,\delta_{c,p,g}^2 + O(\delta_{c,p,g}^3), \qquad \kappa(\mu) = \frac{1}{2}w(\mu)\big(1 - w(\mu)\big).$$
Moreover, for large negative perturbations ($\delta \to -\infty$), $g_\mu(e^\delta)$ saturates at the finite floor $-\log(1 + \mu)$, whereas for large positive perturbations ($\delta \to +\infty$), $g_\mu(e^\delta) \sim \delta + \log w(\mu)$.

---

## 2. Theorem 1 (Deterministic Signal-Tail $w(\mu)$ Conditioning Bias and Split-Half Survival)

Consider the **zero-noise, zero-rewiring idealization** in which perturbation $p$ induces an identical multiplicative log-fold-change $\delta_g = \log\mathrm{FC}_{p,g}$ across all source lines $s \in \mathcal{S}$ and target line $t$ ($\sigma_{\mathrm{noise}} = 0$, $\delta_{t,g} = \delta_{s,g} = \delta_g$). Let
$$\Delta_{\mathrm{TMean},g} = \frac{1}{|\mathcal{S}|}\sum_{s\in\mathcal{S}} \Delta_{s,g} \approx \bar{w}_S(g)\,\delta_g, \qquad \bar{w}_S(g) \triangleq \frac{1}{|\mathcal{S}|}\sum_{s\in\mathcal{S}} w(\mu_{s,g}),$$
and $\Delta_{t,g} \approx w_t(g)\,\delta_g$ where $w_t(g) = w(\mu_{t,g})$.

Let $M \subseteq \{1,\dots,G\}$ be any evaluation gene mask, and let $P_1^\perp$ denote orthogonal projection onto $\mathrm{span}\{\mathbf{1}_M, \Delta_{\mathrm{TMean},M}\}^\perp$ (the legacy 1D residual metric `N_1` / `r_orth_cal`).

**Theorem 1 (First-Order Residual Leakage and Infinite-Sample Split-Half Survival).**
1. *(Exact first-order inner product with zero-target predictors).* Assuming zero-mean perturbation effects on $M$, the 1D-orthogonalized target residual satisfies
   $$P_1^\perp \Delta_{t,M} = \big(w_t - \beta_{1,M}\bar{w}_S\big)\,\delta_M + O(\|\delta\|^2), \qquad \beta_{1,M} = \frac{\sum_{g\in M}\delta_g^2\, w_t(g)\,\bar{w}_S(g)}{\sum_{g\in M}\delta_g^2\, \bar{w}_S(g)^2}.$$
   Its inner product with the **zero-target** log-fold-change vector $E6_M = \delta_M$ (and, up to positive scalar $w(c)$, with any constant-baseline predictor $E1_M(c) = g_c(\overline{\mathrm{FC}}_M) \approx w(c)\delta_M$) equals
   $$\langle P_1^\perp \Delta_{t,M},\; P_1^\perp \delta_M \rangle = \left(\sum_{g\in M}\delta_g^2\right)\cdot a_{M,\delta^2} + O(\|\delta\|^3),$$
   where the $\delta^2$-weighted regression intercept is
   $$a_{M,\delta^2} \triangleq \mathbb{E}_{M,\delta^2}[w_t] - \beta_{1,M}\,\mathbb{E}_{M,\delta^2}[\bar{w}_S], \qquad \mathbb{E}_{M,\delta^2}[u] \triangleq \frac{\sum_{g\in M}\delta_g^2 u_g}{\sum_{g\in M}\delta_g^2}.$$
2. *(Why out-of-sample split-half top-$k$ DEG selection fails — and can amplify bias).* Let $M_{\mathrm{top}k}^{(A)} = \mathrm{Top}_k(|\Delta_{t,g}^{(A)}|)$ be selected on an independent half/donor $A$ and evaluated on disjoint half/donor $B$ (`SA-8 resp_top50_split`). As sample size $n \to \infty$, $\Delta_t^{(A)} \to \Delta_t$ and $\Delta_t^{(B)} \to \Delta_t$ deterministically, so split-half selection converges to the noiseless signal-tail mask $M_{\mathrm{top}k}^* = \{g : w_t(g)|\delta_g| \ge \tau_k\}$. Whenever $0 \le \mathrm{Corr}(w_t, \bar{w}_S) < 1$ across cell lines, conditioning on $w_t(g)|\delta_g| \ge \tau_k$ truncates $w_t$ from below near its ceiling $1$ while $\bar{w}_S$ regresses toward its marginal mean: $\mathrm{Cov}_{M^*,\delta^2}(w_t, \bar{w}_S) / \mathrm{Var}_{M^*,\delta^2}(\bar{w}_S)$ shrinks while $\mathbb{E}_{M^*,\delta^2}[w_t] \to 1$, forcing $a_{M^*,\delta^2} > 0$.
   - *Purification Paradox:* At finite sample size, in-sample `resp_top50` mixes true high-$w_t$ signal genes with noise-selected genes (whose $w_t$ is closer to average). Split-half selection (`SA-8`) eliminates finite-sample shared control noise (Nicol, Shivakumar & Irizarry 2026) **and thereby purifies the high-$w_t$ true signal tail**, which can preserve or **increase** $a_{M,\delta^2} > 0$ and inflate `E1` / `Perm100` scores!
3. *(Whole-transcriptome drift & permutation null inheritance).*
   - On $M = \mathrm{all\_nontarget}$, whenever target line $t$ has higher overall library depth or $w_t$ shift relative to $\mathcal{S}$ (e.g., `vcc2025_h1`, `feng_ipsc`), $a_{\mathrm{all},\delta^2} > 0$ induces false-positive `E1`/`E6`/`Perm100` scores even without top-$k$ selection; conversely, when $a_{\mathrm{all},\delta^2} < 0$ (e.g., `rpe1`, `dld1`, `cd4t_stim48`), `E1`/`E6`/`Perm100` drift systematically negative.
   - For any random gene permutation $\pi$ of $\mu_t$, its expectation $\bar{f}(\overline{\mathrm{FC}}_g) = \mathbb{E}_\pi[g_{\mu_{t,\pi(g)}}(\overline{\mathrm{FC}}_g)] = \int g_u(\overline{\mathrm{FC}}_g)\,dF_{\mu_t}(u)$ is a purely baseline-free convex mixture of $\{g_u\}$, inheriting the identical $a_{M,\delta^2}$ drift.

---

## 3. Theorem 2 (Universal Finite-Rank Target-Free Null Cone & Semi-Partial De-Biasing)

Define the **source-only baseline-free transformation cone** for perturbation $p$:
$$\mathcal{C}(\overline{\mathrm{FC}}_p) \triangleq \overline{\mathrm{conv}}\Big\{ g_c(\overline{\mathrm{FC}}_p) = \log(1 + c\,\overline{\mathrm{FC}}_p) - \log(1 + c) \;:\; c \in (0, \infty) \Big\}.$$

**Theorem 2 (Universal Low-Rank Cone Coverage & Exact Semi-Partial Cure).**
1. *(Inclusion of all baseline-free cheats).* Every constant-baseline predictor (`E1`, `E2(P10)`, `E2(P90)`), asymptotic fold-change predictor (`E6 = \log\overline{\mathrm{FC}}`, `E7 = \overline{\mathrm{FC}} - 1`), and expected permutation null ($\mathbb{E}_\pi[\mathrm{FCT}(\mu_\pi)]$) belongs to $\mathcal{C}(\overline{\mathrm{FC}}_p)$ (with `E7` on the $c \to 0^+$ ray and `E6` on the $c \to +\infty$ ray).
2. *(Low-rank universality — zero target data).* On the normalized family $\{g_c(e^\delta)/\|g_c(e^\delta)\|_2 : c \in [10^{-3}, 10^3],\; \delta \in [-3, 3]\}$, the top-3 singular values capture **`99.987%`** of total variance, and three fixed universal target-free nodes $c \in \{0.05, 1.0, 20.0\}$ span every $g_c$ over 6 orders of magnitude of $c$ with **$\min_{c \in [10^{-3},10^3]} R^2(g_c) = 0.99919$** (median `0.99985`).
3. *(Why naive Pearson after projection fails — the $0/0$ curvature blow-up).* Let
   $$\mathcal{N}_5(p) \triangleq \mathrm{span}\Big\{\mathbf{1},\; \Delta_{\mathrm{TMean},p},\; b_S,\; g_{0.05}(\overline{\mathrm{FC}}_p),\; g_{1.0}(\overline{\mathrm{FC}}_p),\; g_{20.0}(\overline{\mathrm{FC}}_p)\Big\},$$
   where $b_S$ is the source-only grand mean perturbation vector. For any baseline-free predictor $r_{\mathrm{pred}} = \hat{\Delta}_p - \Delta_{\mathrm{TMean},p}$, the orthogonal residual energy fraction $\eta(r_{\mathrm{pred}}) = \|P_{\mathcal{N}_5}^\perp r_{\mathrm{pred}}\|_2 / \|(r_{\mathrm{pred}})_{\mathrm{centered}}\|_2 \le 0.004$. However, standard orthogonal Pearson correlation divides by $\|P_{\mathcal{N}_5}^\perp r_{\mathrm{pred}}\|_2$, re-normalizing a $10^{-3}$ higher-order Taylor remainder to unit length and re-creating spurious `0/0` correlation.
4. *(Exact semi-partial cure).* Defining the **Noise-Calibrated Semi-Partial Null-Cone Correlation**
   $$r_{\mathrm{cone\text{-}sp}}(\Delta_{\mathrm{obs}}, \hat{\Delta}) \triangleq \frac{1}{c_{\mathcal{N}_5}} \frac{\langle P_{\mathcal{N}_5}^\perp r_{\mathrm{obs}},\; P_{\mathcal{N}_5}^\perp r_{\mathrm{pred}} \rangle}{\|P_{\mathcal{N}_5}^\perp r_{\mathrm{obs}}\|_2 \cdot \|(r_{\mathrm{pred}})_{\mathrm{centered}}\|_2} = \eta(r_{\mathrm{pred}}) \cdot r_{\mathrm{cone\text{-}orth}}$$
   satisfies four simultaneous guarantees using **zero target-line information**:
   - **(i) Exact invariance to `GlobalScale` (`R18`) and `Source_Grand_Mean` (`b_S`):** $r_{\mathrm{cone\text{-}sp}} \equiv 0$.
   - **(ii) Uniform suppression of the entire baseline-free cone:** $|r_{\mathrm{cone\text{-}sp}}| \le \eta / c_{\mathcal{N}_5} \approx 0$ (`<= 0.0004` empirically across 48 real human evaluations!) for all `E1`, `E2(P10/P90)`, `E6`, `E7` across all masks (`all_nontarget`, `resp_top50`, `SA-8 resp_top50_split`, `top50_src`).
   - **(iii) Near-zero centering of gene-permuted nulls (`Perm100`):** Eliminates the mean drift $\mathbb{E}_\pi[\mathrm{Perm}_\pi] \in \mathcal{C}(\overline{\mathrm{FC}}_p)$ across all masks.
   - **(iv) Zero attenuation of genuinely target-specific models:** For any true biological rewiring residual $r_{\mathrm{pred}} \perp \mathcal{N}_5$, $\eta = 1 - O(5/G) > 0.995$ (verified in `test_null_cone.py::test_4` and on real `Hybrid_ProjShared` $\eta = 0.996$, `FCT` $\eta = 0.74\text{–}0.94$).

---

## 4. Proposition 3 (Shot-Noise Reversal in `logFC` vs `log1p` Cross-Line Conservation — Q5a)

Because $\hat{\delta}_{c,p,g} = \log(\hat{\mu}_{c,p,g} / \hat{\mu}_{c,g})$ amplifies sampling variance by $1/\mu_{c,g}^2$ on lowly expressed genes whereas $\hat{\Delta}_{c,p,g} = \log(1+\hat{\mu}_{c,p,g}) - \log(1+\hat{\mu}_{c,g})$ damps low-expression noise by $w(\mu_{c,g})^2 \approx \mu_{c,g}^2$:
1. **Raw uncalibrated cross-line correlation** systematically favors `log1p` over `logFC` even when true multiplicative fold-changes are more conserved across cell types.
2. **Noise-ceiling disattenuation** ($\hat{r}_{\mathrm{cal}}(A, B) = \hat{r}(A, B) / (c_{\mathrm{SB},A} c_{\mathrm{SB},B})$) removes differential measurement attenuation, revealing the true conservation hierarchy on expressed genes ($\mu \ge 0.5$: `logFC > log1p` on `25/33` pairs, Wilcoxon $p = 1.47 \times 10^{-3}$).
3. **Target-baseline `FCT` mapping** $\log(1 + \mu_B \exp(\hat{\delta}_A)) - \log(1 + \mu_B)$ simultaneously leverages multiplicative conservation AND restores `log1p` low-expression noise damping on the target baseline, explaining why `FCT` outperforms direct `log1p` transfer on `28/33` (`all`), `32/33` (`expr01`), and **`33/33`** (`expr05`, Wilcoxon $p = 2.33 \times 10^{-10}$) cell-line pairs.
