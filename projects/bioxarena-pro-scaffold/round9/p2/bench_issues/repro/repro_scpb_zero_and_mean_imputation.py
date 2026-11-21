"""Minimal reproduction for scPerturBench issue SCPB-02:
Silent zero-filling and constant mean-imputation of missing method-dataset runs
in Results/*_performance_top{100,5000}.csv vs raw distance CSVs, plus pearson_distance
zero-clipping discrepancy.
"""
import os
import pandas as pd

SCPB_ROOT = os.environ.get(
    "SCPERTURB_ROOT",
    "/tmp/round9_S4/repos/scPerturBench/scPerturBench-main",
)

print("=== 1. Genetic Single: scFoundation silent zero-filling ===")
gs = pd.concat(
    [pd.read_csv(f"{SCPB_ROOT}/Results/Genetic_perturbation/all_dataset_genetic_seed{s}.csv") for s in [1, 2, 3]],
    ignore_index=True,
)
single_ds = [
    "Adamson", "Frangieh", "TianActivation", "TianInhibition",
    "Replogle_exp7", "Replogle_exp8", "Papalexi",
    "Replogle_RPE1essential", "Replogle_K562essential",
]
for s in [1, 2, 3]:
    sub = gs[(gs["seed"] == s) & (gs["DataSet"].isin(single_ds))]
    exp = sub.groupby("DataSet")["perturb"].nunique()
    sf = sub[sub["method"] == "scFoundation"].groupby("DataSet")["perturb"].nunique().reindex(exp.index, fill_value=0)
    miss = (exp - sf)[(exp - sf) > 0].to_dict()
    print(f"  Seed {s}: scFoundation missing single perturbs = {int((exp - sf).sum())}/{int(exp.sum())} -> {miss}")

for fn in ["genetic_single_performance_top100.csv", "genetic_single_performance_top5000.csv"]:
    df = pd.read_csv(f"{SCPB_ROOT}/Results/Genetic_perturbation/{fn}")
    sf_zero = df[(df["method"] == "scFoundation") & (df["ac_score"] == 0)]
    print(f"  {fn}: scFoundation zero rows={len(sf_zero)}/{(df['method']=='scFoundation').sum()} -> {sf_zero.groupby('DataSet').size().to_dict()}")

print("\n=== 2. Chemical Single: cycleCDR zero-fill & PRnet/biolord/chemCPA mean-imputation ===")
chem = pd.read_csv(f"{SCPB_ROOT}/Results/Chemical_perturbation/all_dataset_chemical.csv")
for m in ["cycleCDR", "PRnet", "bioLord", "chemCPA"]:
    cov = chem[(chem["method"] == m) & (chem["DataSet"].str.startswith("sciplex3_")) & (chem["DataSet"] != "sciplex3_comb")]
    tab = cov.groupby(["DataSet", "seed"])["perturb"].nunique().unstack(fill_value=0)
    print(f"  all_dataset_chemical.csv {m} coverage:\n{tab.to_string()}")

c100 = pd.read_csv(f"{SCPB_ROOT}/Results/Chemical_perturbation/chemical_single_performance_top100.csv")
cyc_zero = c100[(c100["method"] == "cycleCDR") & (c100[["cor", "edistance_score", "mse_score", "ac_score"]] == 0).all(axis=1)]
print(f"  chemical_single_performance_top100.csv cycleCDR all-zero rows: {len(cyc_zero)} ({cyc_zero['DataSet'].value_counts().to_dict()})")
for m in ["PRnet", "biolord", "chemCPA"]:
    sub = c100[c100["method"] == m]
    vc = sub.groupby(["DataSet", "cor", "edistance_score", "mse_score"]).size()
    print(f"  chemical_single_performance_top100.csv {m} exact repeated rows: {vc[vc > 1].to_dict()}")

print("\n=== 3. Cellular IID: biolord zeroed in top100/5000 despite present in distance0 ===")
d_iid = pd.read_csv(f"{SCPB_ROOT}/Results/Cellular_context_iid/cellular_iid_distance0.csv")
p100 = pd.read_csv(f"{SCPB_ROOT}/Results/Cellular_context_iid/cellular_iid_performance_top100.csv")
p5k = pd.read_csv(f"{SCPB_ROOT}/Results/Cellular_context_iid/cellular_iid_performance_top5000.csv")
for ds in ["Afriat", "McFarland", "crossSpecies"]:
    d_val = d_iid[(d_iid["method"] == "bioLord") & (d_iid["DataSet"] == ds) & (d_iid["metric"] == "pearson_distance") & (d_iid["DEG"] == 100)]["performance"].mean()
    p1_val = p100[(p100["method"] == "biolord") & (p100["DataSet"] == ds)]["cor"].mean()
    p5_val = p5k[(p5k["method"] == "biolord") & (p5k["DataSet"] == ds)]["cor"].mean()
    print(f"  {ds:12s}: distance0 pcc={d_val:.4f} vs top100 cor={p1_val:.4f}, top5000 cor={p5_val:.4f}")

ds12 = p100.groupby(["method", "DataSet"])["ac_score"].mean().groupby("method").mean().sort_values(ascending=False)
ds9 = p100[~p100["DataSet"].isin(["Afriat", "McFarland", "crossSpecies"])].groupby(["method", "DataSet"])["ac_score"].mean().groupby("method").mean().sort_values(ascending=False)
print(f"  biolord Cellular IID top100 rank: 12-dataset={ds12['biolord']:.4f} (#{int(ds12.index.get_loc('biolord'))+1}) vs 9-dataset={ds9['biolord']:.4f} (#{int(ds9.index.get_loc('biolord'))+1})")

print("\n=== 4. pearson_distance zero-clipping in raw CSVs vs signed cor in summary CSVs ===")
bc_raw = gs[(gs["method"] == "baseControl") & (gs["metric"] == "pearson_distance")]
g100 = pd.read_csv(f"{SCPB_ROOT}/Results/Genetic_perturbation/genetic_single_performance_top100.csv")
bc_sum = g100[g100["method"] == "baseControl"]
print(f"  Genetic raw baseControl pearson_distance == 0: {(bc_raw['performance']==0).sum()}/{len(bc_raw)} ({(bc_raw['performance']==0).mean()*100:.1f}%), min={bc_raw['performance'].min()}")
print(f"  Genetic top100 baseControl cor < 0: {(bc_sum['cor']<0).sum()}/{len(bc_sum)} ({(bc_sum['cor']<0).mean()*100:.1f}%), min={bc_sum['cor'].min()}")
