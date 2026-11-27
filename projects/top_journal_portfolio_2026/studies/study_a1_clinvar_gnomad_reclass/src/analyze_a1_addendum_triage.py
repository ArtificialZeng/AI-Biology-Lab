#!/opt/anaconda3/bin/python3
"""
Study A-1 Addendum: Clinical/Technical 4-Layer Triage of All 55 AD >=2-Star Alerts,
Actionable VUS-to-Benign Liberation Debt Under gnomAD v4.1, and Curation-Velocity Mediation Test.
"""
import json, os, warnings
import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")

BASE = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_a1_clinvar_gnomad_reclass"
CLEAN = os.path.join(BASE, "data", "clean")
RES = os.path.join(BASE, "results")

j = pd.read_parquet(os.path.join(CLEAN, "a1_clinvar_gnomad_exact_joined.parquet"))
cat = pd.read_csv(os.path.join(RES, "a1_ad_highconf_reclassification_catalog.csv"))

# 1. Exact 4-Layer Clinical/Technical Triage of all 55 AD >=2-star P/LP variants (FAF95_max8 >= 1e-4)
MAX_AN = 1614324.0  # 807,162 x 2 alleles in gnomAD v4.1
ESTABLISHED_CONTINENTAL_FOUNDERS = {
    "18-31598655-G-A",    # TTR p.Val142Ile (AFR cardiac amyloidosis founder)
    "2-21006288-C-T",     # APOB p.Arg3527Gln (NFE FH founder)
    "2-21006289-G-A",     # APOB p.Arg3527Trp (EAS/SAS FH founder)
    "11-47342698-G-A",    # MYBPC3 p.Arg502Trp (NFE most common HCM allele)
    "11-2768881-C-T",     # KCNQ1 p.Arg518Ter (NFE LQT1 founder)
    "10-43119548-G-A",    # RET p.Val804Met (reduced-penetrance MEN2A allele)
    "16-23621362-C-T",    # PALB2 p.Trp1038Ter (NFE founder)
    "7-6005918-C-A",      # PMS2 p.Ser46Ile (NFE/AMR low-penetrance Lynch founder)
    "19-38499993-G-A",    # RYR1 p.Gly2434Arg (NFE MH pharmacogenetic founder)
    "19-11120425-C-A",    # LDLR p.Cys681Ter (MID Lebanese FH founder)
    "13-32346841-C-T",    # BRCA2 p.Arg2318Ter (EAS Japanese founder)
    "13-32398209-TTGTA-T",# BRCA2 c.9699_9702del (AMR Spanish/Hispanic founder)
    "16-23629985-CAT-C",  # PALB2 c.2167_2168del (AMR Hispanic/Italian founder)
    "7-5977629-G-A",      # PMS2 p.Arg802Ter (exon 14 pseudogene homology region / SAS)
}

def triage_row(r):
    an_frac = float(r["an_joint"]) / MAX_AN
    if an_frac < 0.80:
        return "L1_LowCallability_TechnicalArtifact"
    if r["faf95_grpmax5"] < 1e-4 and r["faf95_pop8"] in ("fin", "asj"):
        return "L2_BottleneckFounder_FIN_ASJ"
    if r["variant_id"] in ESTABLISHED_CONTINENTAL_FOUNDERS:
        return "L3_EstablishedContinentalFounder_LowPen"
    return "L4_PriorityNonBottleneckAuditCandidate"

cat["an_frac"] = cat["an_joint"] / MAX_AN
cat["triage_layer"] = cat.apply(triage_row, axis=1)
cat["nfe_invisible"] = (cat["faf95_nfe"] < 1e-4)
cat.to_csv(os.path.join(RES, "a1_ad_highconf_triaged_55.csv"), index=False)

triage_summary = {}
for layer, g in cat.groupby("triage_layer"):
    triage_summary[layer] = {
        "n": int(len(g)),
        "share": round(len(g) / len(cat), 4),
        "n_3star": int((g["stars"] == 3).sum()),
        "n_2star": int((g["stars"] == 2).sum()),
        "n_nfe_invisible": int(g["nfe_invisible"].sum()),
        "n_with_hom": int((g["hom_joint"] >= 1).sum()),
        "domains": {str(k): int(v) for k, v in g["domain"].value_counts().items()},
        "pops": {str(k): int(v) for k, v in g["faf95_pop8"].value_counts().items()},
        "variants": g[["variant_id", "gene", "domain", "stars", "last_evaluated", "hgvs_name", "an_joint", "hom_joint", "faf95_nfe", "faf95_max8", "faf95_pop8"]].to_dict(orient="records"),
    }

# 2. Actionable VUS-to-Benign Liberation Debt across ALL 84 ACMG SF v3.3 genes (and 72 AD genes)!
vus_all = j[j["clnsig_cat"] == "VUS"].copy()
vus_ad = vus_all[vus_all["inheritance"] == "AD"].copy()
vus_cs = vus_all[vus_all["csq_group"].isin(["pLoF", "missense", "inframe_indel", "splice_other", "synonymous"])].copy()

def vus_debt_stats(df_v, label):
    res = {"label": label, "n_total_vus_in_gnomad": int(len(df_v)), "tiers": {}}
    # Restrict strictly to high-callability sites (AN >= 80%) so zero low-coverage artifacts pollute VUS liberation!
    hc = df_v[df_v["an_joint"] >= 0.80 * MAX_AN].copy()
    res["n_high_callability_vus"] = int(len(hc))
    for tname, thr in [
        ("faf95_grpmax5_ge_1e4", 1e-4),
        ("faf95_grpmax5_ge_5e4_BS1_AD", 5e-4),
        ("faf95_grpmax5_ge_1e3_BS1_strict", 1e-3),
        ("faf95_grpmax5_ge_5e3_BA1_near", 5e-3),
        ("af_grpmax5_ge_5e2_BA1_standalone", 5e-2),
    ]:
        if "af_grpmax5" in tname:
            hit = hc[hc["af_grpmax5"] >= thr]
            nfe_inv = hit[hit["af_nfe"] < thr]
        else:
            hit = hc[hc["faf95_grpmax5"] >= thr]
            nfe_inv = hit[hit["faf95_nfe"] < thr]
        res["tiers"][tname] = {
            "n_vus_liberable": int(len(hit)),
            "n_nfe_invisible": int(len(nfe_inv)),
            "nfe_invisible_share": round(len(nfe_inv) / len(hit), 4) if len(hit) else 0.0,
            "by_pop5": {str(k): int(v) for k, v in hit["faf95_pop5"].value_counts().items()},
            "by_domain": {str(k): int(v) for k, v in hit["domain"].value_counts().items()},
            "median_submitters": float(hit["n_submitters"].median()) if len(hit) else None,
        }
    return res

vus_lib = {
    "all_84_genes_all_vus": vus_debt_stats(vus_all, "All 84 ACMG SF v3.3 genes (all VUS)"),
    "ad_72_genes_all_vus": vus_debt_stats(vus_ad, "72 AD ACMG SF v3.3 genes (all VUS)"),
    "all_84_genes_coding_splice_vus": vus_debt_stats(vus_cs, "All 84 ACMG SF v3.3 genes (coding/splice VUS)"),
}

# 3. Curation Velocity Mediation Analysis (How much of Non-NFE VUS OR=7.01 is mediated by n_submitters?)
sub = j[j["csq_group"].isin(["pLoF", "missense", "inframe_indel", "splice_other", "synonymous"]) &
        j["clnsig_cat"].isin(["P_LP", "VUS", "CONFLICTING", "B_LB"]) &
        (j["an_joint"] >= 0.80 * MAX_AN)].copy()
is_nn = (sub["af_nonnfe5"] >= 1e-4) & (sub["af_nonnfe5"] / np.maximum(sub["af_nfe"], 1e-6) >= 5.0)
is_eu = (sub["af_nfe"] >= 1e-4) & (sub["af_nfe"] / np.maximum(sub["af_nonnfe5"], 1e-6) >= 5.0)
comp = sub[is_nn | is_eu].copy()
comp["non_nfe_enr"] = is_nn.loc[comp.index].astype(int)
comp["af_ref"] = np.where(comp["non_nfe_enr"] == 1, comp["af_nonnfe5"], comp["af_nfe"])
comp["log10_af"] = np.log10(comp["af_ref"].clip(lower=1e-6))
comp["log_subs"] = np.log1p(comp["n_submitters"].fillna(1).clip(lower=1))
comp["y_vus"] = (comp["clnsig_cat"] == "VUS").astype(int)
comp["y_conf"] = (comp["clnsig_cat"] == "CONFLICTING").astype(int)
comp["y_plp"] = (comp["clnsig_cat"] == "P_LP").astype(int)

# Base model (high-callability AN>=80% sensitivity!) vs Submitters-Adjusted Mediation model
# Drop zero-outcome-variance genes for logistic FE (exact conditional FE support)
def fit_logit_gene_fe(df_in, formula):
    ycol = formula.split("~")[0].strip()
    gv = df_in.groupby("gene")[ycol].nunique()
    valid_g = gv[gv == 2].index
    d_sub = df_in[df_in["gene"].isin(valid_g)].copy()
    return smf.logit(formula, data=d_sub).fit(disp=False, maxiter=150)

m_vus_base = fit_logit_gene_fe(comp, "y_vus ~ non_nfe_enr + log10_af + I(log10_af**2) + C(csq_group) + C(gene)")
m_vus_med  = fit_logit_gene_fe(comp, "y_vus ~ non_nfe_enr + log_subs + log10_af + I(log10_af**2) + C(csq_group) + C(gene)")
m_conf_base = fit_logit_gene_fe(comp, "y_conf ~ non_nfe_enr + log10_af + I(log10_af**2) + C(csq_group) + C(gene)")
m_conf_med  = fit_logit_gene_fe(comp, "y_conf ~ non_nfe_enr + log_subs + log10_af + I(log10_af**2) + C(csq_group) + C(gene)")

# Submitter count disparity controlling for AF, consequence, and gene!
m_subs = smf.ols("log_subs ~ non_nfe_enr + log10_af + I(log10_af**2) + C(csq_group) + C(gene)", data=comp).fit(cov_type="HC3")

# Submitter-count breakdown by individual ancestry
sub_by_anc = {}
for p in ["nfe", "afr", "amr", "eas", "sas", "mid", "asj", "fin"]:
    if p == "nfe":
        g = comp[comp["non_nfe_enr"] == 0]
    else:
        g = sub[(sub[f"af_{p}"] >= 1e-4) & (sub[f"af_{p}"] / np.maximum(sub["af_nfe"], 1e-6) >= 5.0) & (sub["pop_max8"] == p)]
    sub_by_anc[p] = {
        "n": int(len(g)),
        "median_submitters": float(g["n_submitters"].median()) if len(g) else None,
        "mean_submitters": round(float(g["n_submitters"].mean()), 2) if len(g) else None,
        "single_submitter_share": round(float((g["n_submitters"] == 1).mean()), 4) if len(g) else None,
    }

# Also compute Carrier Burden specifically for L4 Priority Non-Bottleneck Audit Candidates!
l4 = cat[cat["triage_layer"] == "L4_PriorityNonBottleneckAuditCandidate"]
l4_burden = {}
for p in ["nfe", "afr", "amr", "eas", "sas", "mid", "asj", "fin"]:
    afs = l4[f"af_{p}"].values.clip(0.0, 1.0)
    cp = 1.0 - np.prod(1.0 - (2.0 * afs * (1.0 - afs) + afs * afs))
    l4_burden[p] = {"per_1000": round(float(cp * 1000.0), 3), "carrier_prob": float(cp)}
for p in l4_burden:
    l4_burden[p]["ratio_vs_nfe"] = round(l4_burden[p]["carrier_prob"] / max(l4_burden["nfe"]["carrier_prob"], 1e-9), 2)

mediation = {
    "n_high_callability_comp": int(len(comp)),
    "log_submitters_gap_nonnfe_vs_nfe": {
        "coef": float(m_subs.params["non_nfe_enr"]),
        "se": float(m_subs.bse["non_nfe_enr"]),
        "p": float(m_subs.pvalues["non_nfe_enr"]),
        "pct_fewer_submitters": round((1.0 - float(np.exp(m_subs.params["non_nfe_enr"]))) * 100.0, 1),
    },
    "vus_base_or": float(np.exp(m_vus_base.params["non_nfe_enr"])),
    "vus_base_p": float(m_vus_base.pvalues["non_nfe_enr"]),
    "vus_submitter_adjusted_or": float(np.exp(m_vus_med.params["non_nfe_enr"])),
    "vus_submitter_adjusted_p": float(m_vus_med.pvalues["non_nfe_enr"]),
    "vus_log_or_attenuation_pct": round((1.0 - float(m_vus_med.params["non_nfe_enr"]) / float(m_vus_base.params["non_nfe_enr"])) * 100.0, 1),
    "conf_base_or": float(np.exp(m_conf_base.params["non_nfe_enr"])),
    "conf_base_p": float(m_conf_base.pvalues["non_nfe_enr"]),
    "conf_submitter_adjusted_or": float(np.exp(m_conf_med.params["non_nfe_enr"])),
    "conf_submitter_adjusted_p": float(m_conf_med.pvalues["non_nfe_enr"]),
    "conf_log_or_attenuation_pct": round((1.0 - float(m_conf_med.params["non_nfe_enr"]) / float(m_conf_base.params["non_nfe_enr"])) * 100.0, 1),
    "submitters_by_ancestry": sub_by_anc,
}

out = {
    "triage_4layer_ad_55": triage_summary,
    "l4_priority_carrier_burden": l4_burden,
    "vus_liberation_debt": vus_lib,
    "curation_velocity_mediation": mediation,
}
out_p = os.path.join(RES, "a1_addendum_triage_and_vus.json")
with open(out_p, "w") as f:
    json.dump(out, f, indent=2)

print("=== 1. 4-LAYER TRIAGE OF ALL 55 AD >=2-STAR ALERTS ===")
for k, v in triage_summary.items():
    print(f"  {k:40s}: n={v['n']:2d} ({v['share']*100:4.1f}%) | 3*={v['n_3star']:2d}, 2*={v['n_2star']:2d} | NFE-invis={v['n_nfe_invisible']:2d} | hom={v['n_with_hom']} | pops={v['pops']}")
print("\n=== 2. L4 PRIORITY AUDIT CANDIDATES (EXACT LIST) ===")
for r in triage_summary["L4_PriorityNonBottleneckAuditCandidate"]["variants"]:
    print(f"  {r['gene']:7s} | {r['variant_id']:20s} | {r['stars']}* | {r['last_evaluated']:12s} | {r['faf95_pop8']:3s} FAF95={r['faf95_max8']:.5f} (NFE={r['faf95_nfe']:.6f}) hom={r['hom_joint']} | {r['hgvs_name'][:52]}")
print("\n=== 3. L4 PRIORITY CARRIER BURDEN PER 1,000 INDIVIDUALS ===")
print(json.dumps(l4_burden, indent=2))
print("\n=== 4. ACTIONABLE VUS-TO-BENIGN LIBERATION DEBT (HIGH-CALLABILITY AN>=80%) ===")
for k, v in vus_lib.items():
    print(f"  [{k}] total_HC_VUS={v['n_high_callability_vus']:,}")
    for tn, tv in v["tiers"].items():
        print(f"    {tn:34s}: n={tv['n_vus_liberable']:5d} | NFE-invisible={tv['n_nfe_invisible']:5d} ({tv['nfe_invisible_share']*100:5.1f}%) | pops={tv['by_pop5']}")
print("\n=== 5. CURATION VELOCITY MEDIATION ===")
print(json.dumps(mediation, indent=2))
