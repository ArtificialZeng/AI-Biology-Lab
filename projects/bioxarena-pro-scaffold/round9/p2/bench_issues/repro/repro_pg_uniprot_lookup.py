"""Minimal reproduction for ProteinGym issue PG-01:
uniprot_function_lookup missing .drop_duplicates() and multi-row merges
before uniprot_level_average in proteingym/performance_DMS_benchmarks.py.
"""
import json
import os
import pandas as pd

PG_ROOT = os.environ.get(
    "PROTEINGYM_ROOT",
    "/tmp/round9_S4/repos/ProteinGym/ProteinGym-main",
)

ref = pd.read_csv(f"{PG_ROOT}/reference_files/DMS_substitutions.csv")
with open(f"{PG_ROOT}/proteingym/constants.json") as f:
    clean_names = json.load(f)["clean_names"]
rev_clean = {v: k for k, v in clean_names.items()}

# Replicate lines 136-141 of proteingym/performance_DMS_benchmarks.py
uniprot_function_lookup_buggy = ref[["UniProt_ID", "coarse_selection_type"]].rename(
    columns={"coarse_selection_type": "Selection Type"}
)
uniprot_function_lookup_dedup = uniprot_function_lookup_buggy.drop_duplicates()
uniprot_Neff_lookup = ref[["UniProt_ID", "MSA_Neff_L_category"]].drop_duplicates()
uniprot_taxon_lookup = ref[["UniProt_ID", "taxon"]].drop_duplicates()

print(f"Reference assays: {len(ref)}, unique UniProt_IDs: {ref['UniProt_ID'].nunique()}")
print(f"uniprot_function_lookup rows: buggy={len(uniprot_function_lookup_buggy)}, dedup={len(uniprot_function_lookup_dedup)}")
print(f"uniprot_Neff_lookup rows: {len(uniprot_Neff_lookup)} (dups: {uniprot_Neff_lookup[uniprot_Neff_lookup.duplicated('UniProt_ID', keep=False)]['UniProt_ID'].unique().tolist()})")

metrics = ["Spearman", "AUC", "MCC", "NDCG", "Top_recall"]
for metric in metrics:
    dms_df = pd.read_csv(
        f"{PG_ROOT}/benchmarks/DMS_zero_shot/substitutions/{metric}/DMS_substitutions_{metric}_DMS_level.csv"
    )
    uni_pub = pd.read_csv(
        f"{PG_ROOT}/benchmarks/DMS_zero_shot/substitutions/{metric}/DMS_substitutions_{metric}_Uniprot_level.csv"
    )
    sum_pub = pd.read_csv(
        f"{PG_ROOT}/benchmarks/DMS_zero_shot/substitutions/{metric}/Summary_performance_DMS_substitutions_{metric}.csv"
    )

    meta_cols = {"UniProt_ID", "MSA_Neff_L_category", "Taxon", "Selection Type"}
    raw_models = [c for c in uni_pub.columns if c not in meta_cols]
    dms_raw = dms_df.rename(columns=rev_clean)

    # Published Uniprot_level rows & Average row
    pub_data = uni_pub[uni_pub["UniProt_ID"].notna()]
    pub_avg_row = uni_pub[uni_pub["UniProt_ID"].isna()][raw_models].iloc[0]
    pub_raw_mean = pub_data[raw_models].mean(axis=0)

    # True unweighted 186-UniProt average (before metadata merges)
    true_186 = dms_raw.groupby("UniProt_ID")[raw_models].mean(numeric_only=True)
    true_186_raw_mean = true_186.mean(axis=0)
    true_186_r3 = true_186_raw_mean.round(3)

    # Only adding .drop_duplicates() at line 136 still yields 202 rows
    m202 = (
        true_186.reset_index()
        .merge(uniprot_Neff_lookup, on="UniProt_ID", how="left")
        .merge(uniprot_taxon_lookup, on="UniProt_ID", how="left")
        .merge(uniprot_function_lookup_dedup, on="UniProt_ID", how="left")
    )

    r_pub_raw = pub_raw_mean.rank(ascending=False, method="min")
    r_186_raw = true_186_raw_mean.rank(ascending=False, method="min")
    r_pub_r3 = pub_avg_row.rank(ascending=False, method="min")
    r_186_r3 = true_186_r3.rank(ascending=False, method="min")

    # Verify Summary_performance final_average is unaffected
    uf_perf = dms_raw.groupby(["UniProt_ID", "Selection Type"])[raw_models].mean(numeric_only=True).reset_index()
    final_avg = uf_perf.groupby("Selection Type")[raw_models].mean(numeric_only=True).mean(axis=0).round(3)
    sum_map = dict(zip(sum_pub["Model_name"], sum_pub[f"Average_{metric}"]))
    max_sum_diff = max(abs(final_avg[m] - sum_map[clean_names.get(m, m)]) for m in raw_models)

    print(
        f"[{metric:10s}] pub_data_rows={len(pub_data)} (m202={len(m202)}, true={len(true_186)}) | "
        f"rank_changes_raw={int((r_pub_raw != r_186_raw).sum())}/{len(raw_models)} "
        f"(round3={int((r_pub_r3 != r_186_r3).sum())}/{len(raw_models)}) | "
        f"ESM3: {pub_raw_mean['ESM3']:.3f}(#{int(r_pub_raw['ESM3'])})->{true_186_raw_mean['ESM3']:.3f}(#{int(r_186_raw['ESM3'])}) | "
        f"max_Summary_diff={max_sum_diff:.3f}"
    )
