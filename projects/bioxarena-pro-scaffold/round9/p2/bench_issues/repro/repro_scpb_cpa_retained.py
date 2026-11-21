"""Minimal reproduction for scPerturBench issue SCPB-01:
README.md (lines 103-107) states CPA was removed from the four genetic perturbation
performance summary files, yet all four summary CSVs on main still retain CPA
(ranking #1 in genetic_single_performance_top100.csv) and CPA anchors min-max bounds.
"""
import os
import pandas as pd

SCPB_ROOT = os.environ.get(
    "SCPERTURB_ROOT",
    "/tmp/round9_S4/repos/scPerturBench/scPerturBench-main",
)

for seed in [1, 2, 3]:
    raw = pd.read_csv(f"{SCPB_ROOT}/Results/Genetic_perturbation/all_dataset_genetic_seed{seed}.csv")
    print(f"all_dataset_genetic_seed{seed}.csv: CPA rows = {(raw['method'] == 'CPA').sum()}")

summary_files = [
    "genetic_single_performance_top100.csv",
    "genetic_single_performance_top5000.csv",
    "genetic_combo_performance_top100.csv",
    "genetic_combo_performance_top5000.csv",
]
for fn in summary_files:
    df = pd.read_csv(f"{SCPB_ROOT}/Results/Genetic_perturbation/{fn}")
    cpa_df = df[df["method"] == "CPA"]
    ds_means = df.groupby(["method", "DataSet"])["ac_score"].mean().groupby("method").mean().sort_values(ascending=False)
    pert_means = df.groupby("method")["ac_score"].mean().sort_values(ascending=False)
    ds_rank = int(ds_means.index.get_loc("CPA")) + 1
    pert_rank = int(pert_means.index.get_loc("CPA")) + 1
    print(
        f"\n{fn}: CPA rows={len(cpa_df)} | "
        f"ds_mean={ds_means['CPA']:.4f} (#{ds_rank}/{len(ds_means)}) | "
        f"pert_mean={pert_means['CPA']:.4f} (#{pert_rank}/{len(pert_means)})"
    )
    print("  Top 4 (with CPA):", {k: round(float(v), 4) for k, v in ds_means.head(4).items()})
    ds_no_cpa = ds_means.drop("CPA")
    print("  Top 4 (w/o CPA) :", {k: round(float(v), 4) for k, v in ds_no_cpa.head(4).items()})
    anchor_cols = [c for c in ["edistance_score", "mse_score", "was_score", "sym_score"] if c in cpa_df.columns]
    anchors = {c: (int((cpa_df[c] == 0).sum()), int((cpa_df[c] == 1).sum())) for c in anchor_cols}
    print("  CPA min-max boundary hits (0s, 1s):", anchors)
