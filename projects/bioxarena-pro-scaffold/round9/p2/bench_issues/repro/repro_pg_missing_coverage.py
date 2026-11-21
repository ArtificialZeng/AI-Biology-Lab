"""Minimal reproduction for ProteinGym issue PG-02:
Incomplete assay coverage in DMS_zero_shot (Protriever: 200/217) and unequal
gene coverage in clinical_zero_shot (N=1,129/2,525) & clinical_supervised (N=1,374/2,525).
"""
import json
import os
import zipfile
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from sklearn.metrics import roc_auc_score

PG_ROOT = os.environ.get(
    "PROTEINGYM_ROOT",
    "/tmp/round9_S4/repos/ProteinGym/ProteinGym-main",
)
CLIN_ZIP = os.environ.get(
    "PROTEINGYM_CLIN_ZIP",
    "/tmp/round9_S4/pg_clinical/zero_shot_clinical_substitutions_scores.zip",
)

print("=== 1. DMS Zero-Shot Substitutions: Protriever missing 17/217 assays ===")
meta_cols = {"DMS ID", "Number of Mutants", "Selection Type", "UniProt ID", "MSA_Neff_L_category", "Taxon"}
for metric in ["Spearman", "AUC", "MCC", "NDCG", "Top_recall"]:
    dms_df = pd.read_csv(f"{PG_ROOT}/benchmarks/DMS_zero_shot/substitutions/{metric}/DMS_substitutions_{metric}_DMS_level.csv")
    sum_df = pd.read_csv(f"{PG_ROOT}/benchmarks/DMS_zero_shot/substitutions/{metric}/Summary_performance_DMS_substitutions_{metric}.csv")
    mcols = [c for c in dms_df.columns if c not in meta_cols]
    na_models = dms_df[mcols].isna().sum()
    na_models = na_models[na_models > 0].to_dict()

    def cat_avg(df):
        uf = df.groupby(["UniProt ID", "Selection Type"])[mcols].mean(numeric_only=True).reset_index()
        return uf.groupby("Selection Type")[mcols].mean(numeric_only=True).mean(axis=0)

    full_cat = cat_avg(dms_df)
    c200_cat = cat_avg(dms_df[dms_df["Protriever"].notna()])
    r_full = full_cat.rank(ascending=False, method="min")
    r_c200 = c200_cat.rank(ascending=False, method="min")
    if metric == "Spearman":
        miss = dms_df["Protriever"].isna()
        print(f"Missing assays={miss.sum()}/{len(dms_df)}, unique UniProt_IDs={dms_df.loc[miss, 'UniProt ID'].nunique()}")
        print(f"Mean mutants: missing17={dms_df.loc[miss, 'Number of Mutants'].mean():.1f}, present200={dms_df.loc[~miss, 'Number of Mutants'].mean():.1f}")
        print(f"Models with NaNs: {na_models}")
    print(
        f"[{metric:10s}] Protriever: {full_cat['Protriever']:.4f}(#{int(r_full['Protriever'])})->{c200_cat['Protriever']:.4f}(#{int(r_c200['Protriever'])}) | "
        f"ESCOTT: {full_cat['ESCOTT']:.4f}(#{int(r_full['ESCOTT'])})->{c200_cat['ESCOTT']:.4f}(#{int(r_c200['ESCOTT'])}) | "
        f"ProSST(K=512): {full_cat['ProSST (K=512)']:.4f}(#{int(r_full['ProSST (K=512)'])})->{c200_cat['ProSST (K=512)']:.4f}(#{int(r_c200['ProSST (K=512)'])})"
    )

print("\n=== 2. Clinical Substitutions (Zero-Shot & Supervised) ===")
for track in ["clinical_zero_shot", "clinical_supervised"]:
    dms = pd.read_csv(f"{PG_ROOT}/benchmarks/{track}/substitutions/AUC/clinical_substitutions_AUC_DMS_level.csv")
    sm = pd.read_csv(f"{PG_ROOT}/benchmarks/{track}/substitutions/AUC/Summary_performance_clinical_substitutions_AUC.csv")
    mcols = [c for c in dms.columns if c != "RefSeq ID"]
    cc = dms.dropna(subset=mcols)
    full_m = dms[mcols].mean()
    cc_m = cc[mcols].mean()
    rf = full_m.rank(ascending=False, method="min").astype(int)
    rc = cc_m.rank(ascending=False, method="min").astype(int)
    cov = dms[mcols].notna().sum()
    print(f"[{track}] summary_cols={list(sm.columns)} | total_genes={len(dms)}, complete_case_N={len(cc)}")
    print(f"  Min/Max gene coverage: {cov.idxmin()}={cov.min()} .. {cov.idxmax()}={cov.max()}")
    changed = [m for m in mcols if rf[m] != rc[m]]
    print(f"  Models with rank change (full vs N={len(cc)}): {changed}")
    for m in changed:
        print(f"    {m:18s}: N={cov[m]:4d}, full={full_m[m]:.4f}(#{rf[m]}) -> cc={cc_m[m]:.4f}(#{rc[m]})")

if os.path.exists(CLIN_ZIP):
    with open(f"{PG_ROOT}/config.json") as f:
        cfg = json.load(f)["model_list_zero_shot_substitutions_clinical"]
    with open(f"{PG_ROOT}/proteingym/constants.json") as f:
        clean_names = json.load(f)["clean_names"]
    with zipfile.ZipFile(CLIN_ZIP) as zf:
        names = [n for n in zf.namelist() if n.endswith(".csv")]
        pooled_rows, gstats = [], []
        for n in names:
            df = pd.read_csv(zf.open(n))
            gid = n.split("/")[-1].replace(".csv", "")
            y = np.array([1 if v == "Pathogenic" else (0 if v == "Benign" else v) for v in df["DMS_bin_score"].values])
            y = 1 - y
            gstats.append({"gene": gid, "n": len(y), "min_cls": int(min((y == 1).sum(), (y == 0).sum()))})
            df["_gene"] = gid
            df["_y"] = y
            pooled_rows.append(df)
    gs = pd.DataFrame(gstats)
    pooled = pd.concat(pooled_rows, ignore_index=True)
    missing_from_zip = [m for m in cfg if m not in pooled.columns]
    models = [m for m in cfg if m in pooled.columns]
    print(f"\nConfig models missing from clinical score zip: {missing_from_zip}")
    print(f"Genes with min_cls == 1: {(gs['min_cls']==1).sum()}/{len(gs)} ({(gs['min_cls']==1).mean()*100:.1f}%)")
    for m in ["MutPred", "TranceptEVE_L", "PoET", "GEMME"]:
        s = pooled[m].values * cfg[m]["directionality"]
        mask = ~np.isnan(s)
        auc_p = roc_auc_score(pooled["_y"].values[mask], s[mask])
        df_z = pd.DataFrame({"g": pooled["_gene"], "y": pooled["_y"], "s": s})[mask]
        df_z["sz"] = df_z.groupby("g")["s"].transform(lambda x: (x - x.mean()) / x.std(ddof=0) if x.std(ddof=0) > 0 else 0.0)
        auc_z = roc_auc_score(df_z["y"].values, df_z["sz"].values)
        print(f"  {clean_names.get(m,m):15s}: pooled_AUC={auc_p:.4f}, per_gene_zscore_pooled_AUC={auc_z:.4f}")

    # Paired bootstrap & Wilcoxon on N=1129
    czs = pd.read_csv(f"{PG_ROOT}/benchmarks/clinical_zero_shot/substitutions/AUC/clinical_substitutions_AUC_DMS_level.csv")
    cc = czs.dropna(subset=[c for c in czs.columns if c != "RefSeq ID"])
    rng = np.random.default_rng(42)
    diff = cc["GEMME"].values - cc["PoET (200M)"].values
    boot = rng.choice(diff, size=(10000, len(diff)), replace=True).mean(axis=1)
    _, pval = wilcoxon(diff)
    print(f"  Paired N={len(cc)} GEMME - PoET: diff={diff.mean():+.4f}, 95%CI=[{np.percentile(boot,2.5):+.4f}, {np.percentile(boot,97.5):+.4f}], Wilcoxon p={pval:.4f}")
