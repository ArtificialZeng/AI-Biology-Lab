#!/opt/anaconda3/bin/python3
import json, os, time, urllib.request, urllib.parse
import pandas as pd

OUT_DIR = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_a1_clinvar_gnomad_reclass/falsification_ar_in_ad_pipelines/v2_upgrade"
opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
UA = {"User-Agent": "A1Xv2UKBBPheno/1.0", "Content-Type": "application/json"}

# 1. Inspect ALL credible sets already downloaded for APOB pLoF, APOB TrueFH, LDLR TrueFH (we saved 100 per variant, let's check if LDL direct / ApoB / Cholesterol / Statin / CAD are in other_top_hits or need broader keyword match!)
with open(os.path.join(OUT_DIR, "ukbb_opentargets_lipid_cad_phenotype_audit.json")) as f:
    ot_all = json.load(f)

# Also query Open Targets for ALL 39 APOB pLoF variants + top 15 LDLR + APOB missense + PCSK9 using larger page size (200) and checking ALL trait names!
# First check what traits we already got in ukbb_opentargets_lipid_cad_phenotype_audit.json!
summary = []
for r in ot_all:
    if r["n_credible_sets"] > 0:
        all_h = r["ldl_hits"] + r["chol_hits"] + r["apob_hits"] + r["cad_statin_hits"] + r["other_top_hits"]
        summary.append({
            "variant_id": r["variant_id"], "gene": r["gene"], "group": r["group"], "hgvs": r["hgvs"],
            "gnomad_het": r["gnomad_het"], "n_cs": r["n_credible_sets"],
            "rsids": r.get("rsids", []),
            "sample_traits": [(h["trait"], h["beta"], f"{h['p_mantissa']:.2f}e{h['p_exp']}" if h["p_mantissa"] else None, h["n_samples"]) for h in all_h[:15]]
        })

# 2. Query Open Targets for MORE LDLR / APOB / PCSK9 variants & also query EBI GWAS Catalog by rsID for ALL 39 APOB pLoFs!
# Load all 39 APOB pLoFs + top APOB FH + LDLR FH + PCSK9 from parquet
BASE = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_a1_clinvar_gnomad_reclass"
df_j = pd.read_parquet(os.path.join(BASE, "data/clean/a1_clinvar_gnomad_exact_joined.parquet"))
def is_plp(s):
    if pd.isna(s): return False
    s = str(s)
    return ("Conflicting" not in s) and (("Pathogenic" in s) or ("Likely_pathogenic" in s))

# Query GWAS Catalog REST API for the top 15 APOB pLoF rsIDs (Backman et al. 2021 UKB 454k Exomes & Karczewski Genebass 2022 are indexed in GWAS Catalog!)
apob_rsids = []
for r in ot_all:
    if r["gene"] in ("APOB", "LDLR", "PCSK9") and r["gnomad_het"] >= 10:
        row_m = df_j[df_j["variant_id"] == r["variant_id"]]
        rs = str(row_m["rs"].iloc[0]) if len(row_m) > 0 and pd.notna(row_m["rs"].iloc[0]) else ""
        if rs.startswith("rs"):
            apob_rsids.append((r["variant_id"], r["gene"], r["group"], r["hgvs"], r["gnomad_het"], rs))

gwas_cat_hits = []
for vid, gene, grp, hgvs, het, rsid in apob_rsids[:35]:
    try:
        u = f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{rsid}/associations?projection=associationBySnp"
        req = urllib.request.Request(u, headers={"User-Agent": "A1Xv2GWASCat/1.0", "Accept": "application/json"})
        with opener.open(req, timeout=15) as resp:
            d = json.loads(resp.read().decode("utf-8"))
        assocs = (d.get("_embedded") or {}).get("associations", [])
        parsed = []
        for a in assocs:
            beta_num = a.get("betaNum")
            beta_dir = a.get("betaDirection")
            beta_unit = a.get("betaUnit")
            or_val = a.get("orPerCopyNum")
            pval = a.get("pvalue")
            efos = [t.get("trait") for t in (a.get("efoTraits") or [])]
            loci = a.get("loci") or []
            risk_a = [ra.get("riskAlleleName") for loc in loci for ra in (loc.get("strongestRiskAlleles") or [])]
            parsed.append({
                "traits": efos, "beta": beta_num, "dir": beta_dir, "unit": beta_unit,
                "or": or_val, "p": pval, "risk_alleles": risk_a,
            })
        if parsed:
            gwas_cat_hits.append({
                "variant_id": vid, "rsid": rsid, "gene": gene, "group": grp, "hgvs": hgvs,
                "gnomad_het": het, "n_assocs": len(parsed), "assocs": parsed[:25]
            })
    except Exception as e:
        pass
    time.sleep(0.15)

# 3. Also query Open Targets GraphQL for APOB, LDLR, PCSK9, RYR1 target->known gene-burden / rare-variant genetics (especially UKB 450k exome gene burden!)
q_gene = """
query GeneStudies($ensg: String!) {
  target(ensemblId: $ensg) {
    id
    approvedSymbol
    associatedDiseases(page: {index: 0, size: 15}) {
      rows {
        disease { id name }
        score
        datatypeScores { id score }
      }
    }
  }
}
"""
# Also query Open Targets direct `variant` -> `evidences` OR broader credibleSets (size 200) for top 5 APOB pLoFs & APOB Arg3527Gln & LDLR controls
# Actually check Genebass / Hail public API endpoint or Open Targets `search`!
# Wait: check exact UKB 450k Exome APOB pLoF variants via GWAS Catalog (Backman 2021 Nature GCST90081xxx / GCST90079xxx!)
out = {
    "opentargets_variants_with_cs": summary,
    "gwas_catalog_hits": gwas_cat_hits,
}
with open(os.path.join(OUT_DIR, "ukbb_gwas_catalog_and_ot_expanded.json"), "w") as f:
    json.dump(out, f, indent=2)
print(f"DONE! OT CS variants={len(summary)}, GWAS Catalog matched variants={len(gwas_cat_hits)}")
