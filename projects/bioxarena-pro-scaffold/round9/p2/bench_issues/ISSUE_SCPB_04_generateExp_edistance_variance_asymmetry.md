# [Evaluation] Asymmetric `generateExp(..., control_std)` variance injection heavily penalizes deterministic models under `edistance` in 5,000-gene space

**Commit checked:** `63fd9b95071b894c4f39c8f7eb2bd42ca9a70a77` (`main`)

## Summary

In the perturbation generalization scripts, linear and simple baselines sample synthetic cell-level dispersion around predicted means via `generateExp(cellNum, means, control_std)`:
- Genetic: `myGenePert.py:91-95,116`, `linearModel.py:79-83,104`, `linearModel_mean.py:113-117,137`, `genetic_baseMLP.py:215-219,240`, `genetic_baseReg.py:153-157,179`
- Chemical: `chemical_baseMLP.py:180-184,207`, `chemical_baseReg.py:100-104,127`, `linearModel_mean_chemical.py:99-103,126`
- Cellular context: `i.i.d./baseControl.py`, `i.i.d./trainMean.py`, `o.o.d./trainMean.py`

In contrast, deep and foundation models (`myscGPT1.py`, `myscouter.py`, `myscFoundation.py`, `myChemCPA.py`, `myCycleCDR_pre.py`, etc.) output denoised conditional expectations or low-variance predictions without control dispersion injection.

Because `pertpy`'s Euclidean energy distance is $\mathrm{E\text{-}dist}(X, Y) = 2\mathbb{E}\|X - Y\|_2 - \mathbb{E}\|X - X'\|_2 - \mathbb{E}\|Y - Y'\|_2$, a zero-variance predictor suffers a missing $\mathbb{E}\|Y - Y'\|_2 \approx \sqrt{2\sum_{g=1}^D \sigma_g^2}$ term that grows as $\sqrt{D}$—penalizing even an **oracle mean predictor** by >10× in `D = 5000` gene space while rewarding `baseControl` (which has `cor ≈ 0` but full `control_std`).

**Empirical impact in published `DEG=5000` dataset means:**
- **Genetic Single (`top5000`):** `scGPT` ranks **#1 in `mse_score` (`0.9531`, raw `mse = 0.003803`)** and **#2 in `pearson_distance` (`0.3288`, `cor = 0.3247`)**, yet has `edistance = 11.19` (vs `GenePert` `0.77`, `trainMean` `0.80`, `baseControl` `0.91`), pulling its `edistance_score` down to `0.2359` (#14) and overall `ac_score` to `0.5355` (#12 with CPA / #11 without CPA), **below `baseControl` (`0.5676`, `cor = 0.0026`)**. `scouter` shows the same gap (`pcc-delta = 0.3207`, `edistance = 11.49`, `ac_score = 0.4535`).
- **Chemical Single (`top5000`):** `chemCPA` ranks **#1 in `mse` (`0.003777`, `mse_score = 0.9916`)** and **#1 in `cor` (`0.3981`)**, and `biolord` ranks **#2 in both (`mse = 0.004066`, `cor = 0.3594`)**, yet `edistance = 14.83` (`chemCPA`) and `16.41` (`biolord`) vs `0.97` (`trainMean`) and `1.03` (`baseControl`) push them to **#5 (`0.6760`)** and **#6 (`0.5261`)**, both **below `baseControl` (`#3`, `0.7326`, `cor = -0.0016`)**.

## Minimal Reproduction

```python
import numpy as np
from scipy.spatial.distance import cdist

def edistance(x, y):
    return 2.0 * cdist(x, y).mean() - cdist(x, x).mean() - cdist(y, y).mean()

rng = np.random.default_rng(42)
n_cells = 300
for d in [100, 5000]:
    sigma = rng.uniform(0.25, 0.45, size=d)
    mu_ctrl, mu_stim = np.zeros(d), np.zeros(d)
    mu_stim[: min(50, d)] = 0.8
    x_stim = rng.normal(mu_stim, sigma, size=(n_cells, d))

    y_oracle_mean = np.tile(mu_stim, (n_cells, 1))                # exact mean, zero variance
    y_oracle_genexp = rng.normal(mu_stim, sigma, (n_cells, d))    # exact mean + generateExp()
    y_base_control = rng.normal(mu_ctrl, sigma, (n_cells, d))     # wrong mean + generateExp()

    print(
        f"D={d:4d} | Oracle(mean)={edistance(x_stim, y_oracle_mean):6.3f} | "
        f"Oracle+generateExp={edistance(x_stim, y_oracle_genexp):5.3f} | "
        f"baseControl={edistance(x_stim, y_base_control):5.3f}"
    )
```

Output:
```
D= 100 | Oracle(mean)= 2.079 | Oracle+generateExp=0.029 | baseControl=5.152
D=5000 | Oracle(mean)=14.812 | Oracle+generateExp=0.237 | baseControl=1.134
```

## Suggested Improvements

1. Evaluate distributional distances (`edistance`, `wasserstein`, `sym_kldiv`) **only** under a symmetric post-processing protocol (either disabling `generateExp` across all models—leaving mean accuracy to `MSE`/`PCC-delta`—or applying uniform `control_std` dispersion calibration to all deterministic predictors).
2. Report mean-prediction metrics (`MSE`, `PCC-delta`, `DEG overlap`) separately from distributional spread metrics (`edistance`, `wasserstein`, `sym_kldiv`) so `baseControl` dispersion matching cannot override `cor ≈ 0` in composite `ac_score` rankings at `D = 5000`.
