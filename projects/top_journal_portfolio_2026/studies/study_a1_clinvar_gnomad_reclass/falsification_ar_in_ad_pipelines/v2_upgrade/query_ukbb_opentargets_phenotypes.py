#!/opt/anaconda3/bin/python3
import gzip, json, os, time, urllib.request
import pandas as pd

BASE = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_a1_clinvar_gnomad_reclass"
OUT_DIR = os.path.join(BASE, "falsification_ar_in_ad_pipelines/v2_upgrade")
os.makedirs(OUT_DIR, exist_ok=True)

df_j = pd.read_parquet(os.path.join(BASE, "data/clean/a1_clinvar_gnomad_exact_joined.parquet"))
def is_plp(s):
    if pd.isna(s): return False
    s = str(s)
    return ("Conflicting" not in s) and (("Pathogenic" in s) or ("Likely_pathogenic" in s))

ad_j = df_j[(df_j["inheritance"] == "AD") & (df_j["clnsig"].apply(is_plp))].copy()
vids_all = set(ad_j["variation_id"].astype(str))
vcf_mc = {}
with gzip.open(os.path.join(BASE, "data/raw/clinvar.vcf.gz"), "rt", errors="replace") as fh:
    for line in fh:
        if line.startswith("#"): continue
        p = line.rstrip("\n").split("\t")
        if p[2] in vids_all:
            info = dict(kv.split("=", 1) for kv in p[7].split(";") if "=" in kv)
            vcf_mc[p[2]] = info.get("MC", "")

ad_j["vcf_mc"] = ad_j["variation_id"].astype(str).map(lambda x: vcf_mc.get(x, ""))
ad_j["is_plof"] = ad_j["vcf_mc"].str.lower().apply(lambda s: any(k in s for k in ["nonsense", "frameshift", "splice_donor", "splice_acceptor"]))
ad_j["het_joint"] = (ad_j["ac_joint"] - 2 * ad_j["hom_joint"]).clip(lower=0)

# Build query target sets:
# 1. ALL 2-4* APOB pLoF (n=39)
apob_plof = ad_j[(ad_j["gene"] == "APOB") & (ad_j["stars"] >= 2) & (ad_j["is_plof"])].sort_values("het_joint", ascending=False)
# 2. ALL 2-4* APOB non-pLoF (true AD FH missense!)
apob_fh = ad_j[(ad_j["gene"] == "APOB") & (ad_j["stars"] >= 2) & (~ad_j["is_plof"])].sort_values("het_joint", ascending=False)
# 3. Top 20 2-4* LDLR true AD FH controls (by het count!)
ldlr_fh = ad_j[(ad_j["gene"] == "LDLR") & (ad_j["stars"] >= 2)].sort_values("het_joint", ascending=False).head(20)
# 4. ALL 2-4* PCSK9 pLoF (opposite-direction low-LDL!) and non-pLoF (true AD FH GoF!)
pcsk9_all = ad_j[(ad_j["gene"] == "PCSK9") & (ad_j["stars"] >= 2)].sort_values("het_joint", ascending=False)
# 5. All 10 3* RYR1 AR-only alleles
ryr1_ar = pd.read_csv(os.path.join(BASE, "falsification_ar_in_ad_pipelines/ryr1_10_3star_cm_vcep_ar_alleles_gnomad_v41.csv"))

targets = []
for _, r in apob_plof.iterrows():
    targets.append((r["variant_id"], "APOB", "APOB_pLoF_FHBL", r["hgvs_name"], int(r["het_joint"]), int(r["stars"])))
for _, r in apob_fh.iterrows():
    targets.append((r["variant_id"], "APOB", "APOB_Missense_TrueFH", r["hgvs_name"], int(r["het_joint"]), int(r["stars"])))
for _, r in ldlr_fh.iterrows():
    targets.append((r["variant_id"], "LDLR", "LDLR_TrueFH", r["hgvs_name"], int(r["het_joint"]), int(r["stars"])))
for _, r in pcsk9_all.iterrows():
    grp = "PCSK9_pLoF_LowLDL" if r["is_plof"] else "PCSK9_Missense_TrueFH"
    targets.append((r["variant_id"], "PCSK9", grp, r["hgvs_name"], int(r["het_joint"]), int(r["stars"])))
for _, r in ryr1_ar.iterrows():
    targets.append((r["variant_id"], "RYR1", "RYR1_3star_AR_Myopathy", r["hgvs_name"], int(r["het_joint"]), 3))

print(f"Total variants to query against Open Targets Platform GraphQL: {len(targets)}")

opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
ot_url = "https://api.platform.opentargets.org/api/v4/graphql"
UA = {"User-Agent": "A1Xv2PhenotypeAudit/1.0", "Content-Type": "application/json"}

q_ot = """
query VarCred($vid: String!) {
  variant(variantId: $vid) {
    id
    rsIds
    credibleSets(studyTypes: [gwas], page: {index: 0, size: 100}) {
      count
      rows {
        study { id traitFromSource nSamples }
        beta
        standardError
        pValueMantissa
        pValueExponent
        finemappingMethod
      }
    }
  }
}
"""

# Also extract specific LDL-C, Total Cholesterol, ApoB, Lipid-lowering/Statin, CAD/Heart disease hits!
results = []
for vid_dash, gene, grp, hgvs, het, stars in targets:
    vid_us = vid_dash.replace("-", "_")
    rec = {
        "variant_id": vid_dash, "ot_id": vid_us, "gene": gene, "group": grp,
        "hgvs": hgvs, "gnomad_het": het, "stars": stars, "n_credible_sets": 0,
        "ldl_hits": [], "chol_hits": [], "apob_hits": [], "cad_statin_hits": [], "other_top_hits": []
    }
    try:
        req = urllib.request.Request(ot_url, data=json.dumps({"query": q_ot, "variables": {"vid": vid_us}}).encode(), headers=UA)
        with opener.open(req, timeout=20) as r:
            d = json.loads(r.read().decode("utf-8"))
        v = (d.get("data") or {}).get("variant")
        if v and v.get("credibleSets"):
            cs = v["credibleSets"]
            rec["n_credible_sets"] = cs.get("count", 0)
            rec["rsids"] = v.get("rsIds", [])
            for row in cs.get("rows", []):
                st = row.get("study") or {}
                tr = str(st.get("traitFromSource") or "")
                tl = tr.lower()
                item = {
                    "study_id": st.get("id"),
                    "trait": tr,
                    "n_samples": st.get("nSamples"),
                    "beta": row.get("beta"),
                    "se": row.get("standardError"),
                    "p_mantissa": row.get("pValueMantissa"),
                    "p_exp": row.get("pValueExponent"),
                    "method": row.get("finemappingMethod"),
                }
                if any(k in tl for k in ["ldl", "low density lipoprotein", "low-density lipoprotein"]):
                    rec["ldl_hits"].append(item)
                elif "apolipoprotein b" in tl or "apob" in tl:
                    rec["apob_hits"].append(item)
                elif "cholesterol" in tl and "hdl" not in tl:
                    rec["chol_hits"].append(item)
                elif any(k in tl for k in ["coronary", "myocardial", "heart disease", "statin", "atorvastatin", "simvastatin", "cholesterol lowering", "lipid lowering", "hypercholesterolemia", "ischemic heart"]):
                    rec["cad_statin_hits"].append(item)
                elif len(rec["other_top_hits"]) < 5:
                    rec["other_top_hits"].append(item)
    except Exception as e:
        rec["error"] = str(e)
    results.append(rec)
    time.sleep(0.12)

out_json = os.path.join(OUT_DIR, "ukbb_opentargets_lipid_cad_phenotype_audit.json")
with open(out_json, "w") as f:
    json.dump(results, f, indent=2)

# Build clean summary CSV of every variant with >=1 LDL/Cholesterol/ApoB/CAD hit!
rows = []
for r in results:
    best_ldl = min(r["ldl_hits"], key=lambda x: (x["p_exp"] if x["p_exp"] is not None else 0, x["p_mantissa"] or 9)) if r["ldl_hits"] else None
    best_chol = min(r["chol_hits"], key=lambda x: (x["p_exp"] if x["p_exp"] is not None else 0, x["p_mantissa"] or 9)) if r["chol_hits"] else None
    best_apob = min(r["apob_hits"], key=lambda x: (x["p_exp"] if x["p_exp"] is not None else 0, x["p_mantissa"] or 9)) if r["apob_hits"] else None
    best_cad = min(r["cad_statin_hits"], key=lambda x: (x["p_exp"] if x["p_exp"] is not None else 0, x["p_mantissa"] or 9)) if r["cad_statin_hits"] else None
    if any([best_ldl, best_chol, best_apob, best_cad]) or r["n_credible_sets"] > 0:
        rows.append({
            "variant_id": r["variant_id"], "gene": r["gene"], "group": r["group"],
            "hgvs": r["hgvs"], "gnomad_het": r["gnomad_het"], "n_cs": r["n_credible_sets"],
            "ldl_beta": best_ldl["beta"] if best_ldl else None,
            "ldl_p": f"{best_ldl['p_mantissa']:.2f}e{best_ldl['p_exp']}" if best_ldl and best_ldl["p_mantissa"] else None,
            "ldl_n": best_ldl["n_samples"] if best_ldl else None,
            "ldl_trait": best_ldl["trait"] if best_ldl else None,
            "chol_beta": best_chol["beta"] if best_chol else None,
            "chol_p": f"{best_chol['p_mantissa']:.2f}e{best_chol['p_exp']}" if best_chol and best_chol["p_mantissa"] else None,
            "chol_n": best_chol["n_samples"] if best_chol else None,
            "chol_trait": best_chol["trait"] if best_chol else None,
            "apob_beta": best_apob["beta"] if best_apob else None,
            "apob_p": f"{best_apob['p_mantissa']:.2f}e{best_apob['p_exp']}" if best_apob and best_apob["p_mantissa"] else None,
            "cad_beta": best_cad["beta"] if best_cad else None,
            "cad_p": f"{best_cad['p_mantissa']:.2f}e{best_cad['p_exp']}" if best_cad and best_cad["p_mantissa"] else None,
            "cad_trait": best_cad["trait"] if best_cad else None,
        })
pd.DataFrame(rows).to_csv(os.path.join(OUT_DIR, "ukbb_opentargets_lipid_cad_summary.csv"), index=False)
print("DONE! Matched variants with credible sets:", len(rows))
