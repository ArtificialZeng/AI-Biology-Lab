"""Minimal reproduction for scPerturBench issue SCPB-04:
Evaluation asymmetry: linear/simple baselines inject control_std Gaussian noise via
generateExp() whereas deep/foundation models output denoised means, incurring ~11-15x
higher Euclidean energy distance in 5,000-gene space even when MSE and PCC-delta are top-ranked.
"""
import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

SCPB_ROOT = os.environ.get(
    "SCPERTURB_ROOT",
    "/tmp/round9_S4/repos/scPerturBench/scPerturBench-main",
)

def edistance(x: np.ndarray, y: np.ndarray) -> float:
    """Standard Euclidean energy distance: 2 E||X-Y|| - E||X-X'|| - E||Y-Y'||."""
    d_xy = cdist(x, y, metric="euclidean").mean()
    d_xx = cdist(x, x, metric="euclidean").mean()
    d_yy = cdist(y, y, metric="euclidean").mean()
    return float(2.0 * d_xy - d_xx - d_yy)

print("=== 1. Self-contained numerical verification (D=5000 and D=100) ===")
rng = np.random.default_rng(42)
n_cells = 300
for d_genes in [100, 5000]:
    sigma = rng.uniform(0.25, 0.45, size=d_genes)
    mu_ctrl = np.zeros(d_genes)
    mu_stim = np.zeros(d_genes)
    mu_stim[: min(50, d_genes)] = 0.8  # true perturbation shift on 50 DEGs
    x_stim = rng.normal(loc=mu_stim, scale=sigma, size=(n_cells, d_genes))

    # Case A: Oracle mean predictor (zero MSE error vs mu_stim, zero intra-sample variance)
    y_oracle_mean = np.tile(mu_stim, (n_cells, 1))
    # Case B: Same oracle mean + generateExp(..., control_std=sigma)
    y_oracle_noisy = rng.normal(loc=mu_stim, scale=sigma, size=(n_cells, d_genes))
    # Case C: baseControl (wrong mean mu_ctrl, but has control_std=sigma)
    y_ctrl_noisy = rng.normal(loc=mu_ctrl, scale=sigma, size=(n_cells, d_genes))

    print(
        f"  D={d_genes:4d} | "
        f"Oracle(denoised mean) edist={edistance(x_stim, y_oracle_mean):7.3f} | "
        f"Oracle+generateExp edist={edistance(x_stim, y_oracle_noisy):6.3f} | "
        f"baseControl edist={edistance(x_stim, y_ctrl_noisy):6.3f}"
    )

print("\n=== 2. Empirical discrepancy in published scPerturBench CSVs (DEG=5000) ===")
gs = pd.concat(
    [pd.read_csv(f"{SCPB_ROOT}/Results/Genetic_perturbation/all_dataset_genetic_seed{s}.csv") for s in [1, 2, 3]],
    ignore_index=True,
)
single_ds = [
    "Adamson", "Frangieh", "TianActivation", "TianInhibition",
    "Replogle_exp7", "Replogle_exp8", "Papalexi",
    "Replogle_RPE1essential", "Replogle_K562essential",
]
gs_5k = gs[(gs["DataSet"].isin(single_ds)) & (gs["DEG"] == 5000)]
tab_g = gs_5k.groupby(["method", "DataSet", "metric"])["performance"].mean().groupby(["method", "metric"]).mean().unstack()
print("Genetic Single DEG=5000 (dataset mean) selected methods:")
print(tab_g.loc[["GenePert", "trainMean", "linearModel", "baseControl", "scGPT", "scouter"], ["edistance", "mse", "pearson_distance"]].to_string())

chem = pd.read_csv(f"{SCPB_ROOT}/Results/Chemical_perturbation/all_dataset_chemical.csv")
cs_5k = chem[(chem["DataSet"].isin(["sciplex3_A549", "sciplex3_K562", "sciplex3_MCF7"])) & (chem["DEG"] == 5000)]
tab_c = cs_5k.groupby(["method", "DataSet", "metric"])["performance"].mean().groupby(["method", "metric"]).mean().unstack()
print("\nChemical Single DEG=5000 (dataset mean) selected methods:")
print(tab_c.loc[["trainMean", "baseControl", "chemCPA", "bioLord", "cycleCDR"], ["edistance", "mse", "pearson_distance"]].to_string())
