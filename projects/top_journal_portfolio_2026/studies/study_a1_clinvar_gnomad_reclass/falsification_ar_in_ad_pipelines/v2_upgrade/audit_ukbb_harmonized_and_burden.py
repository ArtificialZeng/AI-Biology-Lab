#!/opt/anaconda3/bin/python3
import json, urllib.request, urllib.parse
import pandas as pd

v2 = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_a1_clinvar_gnomad_reclass/falsification_ar_in_ad_pipelines/v2_upgrade"
df = pd.read_csv(f"{v2}/unbiased_7427_ad_gold_labels.csv")

print("=== Breakdown of NOT_MEETING categories in 2-4star >=80% AN (n=186) ===")
sub = df[(df["stars"] >= 2) & (df["hi_call"]) & (df["gold_cat"] != "MEETS_ACMG_SF")]
print(f"Total non-SF in 2-4star >=80% AN: {len(sub)}")
print(sub.groupby(["gold_cat", "gene"]).size().sort_values(ascending=False).to_string())

print("\n=== Check PCSK9 in unbiased_7427 ===")
pcsk9 = df[df["gene"] == "PCSK9"]
print(pcsk9[["variant_id", "hgvs_name", "stars", "is_plof_mc", "gold_cat", "het_joint"]].head(20).to_string())

opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
ot_url = "https://api.platform.opentargets.org/api/v4/graphql"

gcsts = [
    "GCST90079024", "GCST90079033", "GCST90079018", "GCST90077883",
    "GCST90081493", "GCST90079766", "GCST90079765", "GCST90078178",
    "GCST90078226", "GCST90498120", "GCST90499022", "GCST90499437",
    "GCST90497920", "GCST90691895", "GCST007848", "GCST90692159"
]
q = """
query StudyInfo($id: String!) {
  study(studyId: $id) {
    id
    studyType
    traitFromSource
    publicationFirstAuthor
    publicationDate
    publicationJournal
    pubmedId
    nSamples
    nCases
    nControls
    cohorts
  }
}
"""
study_meta = {}
for sid in gcsts:
    req = urllib.request.Request(
        ot_url,
        data=json.dumps({"query": q, "variables": {"id": sid}}).encode(),
        headers={"Content-Type": "application/json", "User-Agent": "A1Xv2Audit/1.0"}
    )
    try:
        with opener.open(req, timeout=15) as r:
            d = json.loads(r.read().decode())
            study_meta[sid] = d.get("data", {}).get("study")
    except Exception as e:
        study_meta[sid] = {"error": str(e)}

print("\n=== Open Targets Study Metadata ===")
for sid, m in study_meta.items():
    if m and "error" not in m:
        print(f"{sid}: {m.get('publicationFirstAuthor')} ({m.get('publicationDate')}) {m.get('publicationJournal')} PMID={m.get('pubmedId')} N={m.get('nSamples')} cases={m.get('nCases')} | {m.get('traitFromSource')}")
    else:
        print(f"{sid}: {m}")

# Query Open Targets for ALL UKB PICS studies across ALL 6 variants PLUS APOB p.Arg3527Trp & PCSK9 pLoF (if any!)
q2 = """
query VarCS($vid: String!) {
  variant(variantId: $vid) {
    id
    rsIds
    alleleFrequencies { populationName alleleFrequency }
    credibleSets(page: {index: 0, size: 150}) {
      count
      rows {
        study { id traitFromSource publicationFirstAuthor publicationDate pubmedId nSamples nCases nControls }
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
targets = [
    ("2_21006019_CA_C", "APOB p.Gly3617fs", "APOB_pLoF_FHBL"),
    ("2_21032391_G_A", "APOB p.Arg439Ter", "APOB_pLoF_FHBL"),
    ("2_21006288_C_T", "APOB p.Arg3527Gln", "APOB_Missense_TrueFH"),
    ("19_11113620_G_A", "LDLR p.Asp482Asn", "LDLR_TrueFH"),
    ("19_11120436_C_T", "LDLR p.Pro685Leu", "LDLR_TrueFH"),
    ("19_11102774_G_A", "LDLR p.Glu101Lys", "LDLR_TrueFH"),
]
# Add PCSK9 pLoF alleles (if any in df)
for _, r in pcsk9.iterrows():
    if r["gold_cat"] != "MEETS_ACMG_SF" and int(r["het_joint"]) >= 5:
        ot_v = str(r["variant_id"]).replace("-", "_")
        targets.append((ot_v, f"PCSK9 {r['hgvs_name']}", r["gold_cat"]))

same_cohort_rows = []
for vid, label, grp in targets:
    req = urllib.request.Request(
        ot_url,
        data=json.dumps({"query": q2, "variables": {"vid": vid}}).encode(),
        headers={"Content-Type": "application/json", "User-Agent": "A1Xv2Audit/1.0"}
    )
    try:
        with opener.open(req, timeout=20) as r:
            vd = json.loads(r.read().decode())["data"]["variant"]
        if not vd:
            continue
        for row in vd["credibleSets"]["rows"]:
            sid = row["study"]["id"]
            if row["finemappingMethod"] == "PICS" or sid.startswith("GCST9007") or sid.startswith("GCST9008") or sid.startswith("GCST9049"):
                same_cohort_rows.append({
                    "variant_id": vid, "label": label, "group": grp,
                    "study_id": sid,
                    "author": row["study"]["publicationFirstAuthor"],
                    "pmid": row["study"]["pubmedId"],
                    "trait": row["study"]["traitFromSource"],
                    "n": row["study"]["nSamples"],
                    "n_cases": row["study"]["nCases"],
                    "beta": row["beta"], "se": row["standardError"],
                    "p": f"{row['pValueMantissa']:.2f}e{row['pValueExponent']}",
                    "p_exp": row["pValueExponent"],
                    "method": row["finemappingMethod"],
                })
    except Exception as e:
        print(f"Err {vid}: {e}")

print(f"\n=== Same-cohort UKB PICS rows: {len(same_cohort_rows)} ===")
sdf = pd.DataFrame(same_cohort_rows)
if len(sdf):
    print(sdf[["label", "study_id", "author", "pmid", "n", "beta", "se", "p", "trait"]].to_string(index=False))

# Check EBI GWAS Catalog REST API for exact Backman 2021 (PMID 34662886) and Nightingale NMR / UKB Exome studies!
# Also query Open Targets for APOB / LDLR / PCSK9 gene-level Exome rare-variant burden via target(ensemblId) -> associatedDiseases / geneticConstraint / etc.
# And check AZPheWAS & Genebass public endpoints!
az_probe = {}
for url in [
    "https://azphewas.com/api/v1/gene/APOB/phenotypes?release=470k&model=ptv",
    "https://azphewas.com/api/v1/genes/APOB",
    "https://app.genebass.org/api/gene/ENSG00000084674",
    "https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/rs982371659/associations",
    "https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/rs142066904/associations",
    "https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/rs5742904/associations",
    "https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/rs139624145/associations",
]:
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "A1Xv2Audit/1.0", "Accept": "application/json"})
        with opener.open(req, timeout=12) as r:
            body = r.read().decode()
            az_probe[url] = {"status": r.status, "len": len(body), "preview": body[:350]}
            print(f"PROBE OK [{r.status}] len={len(body)} {url}")
    except Exception as e:
        az_probe[url] = {"error": str(e)}
        print(f"PROBE ERR {url} -> {e}")

# Also query Open Targets for ALL 39 APOB pLoF variants in 2-4star AND top 15 APOB 1star pLoF + PCSK9 pLoF (all 95 APOB pLoF!)
apob_all_plof = df[(df["gene"] == "APOB") & (df["gold_cat"] == "NOT_MEETING_APOB_FHBL_LOWLDL")].sort_values("het_joint", ascending=False)
print(f"\nScanning ALL {len(apob_all_plof)} APOB pLoF alleles in Open Targets (het_joint >= 5)...")
apob_extra_hits = []
for _, r in apob_all_plof[apob_all_plof["het_joint"] >= 5].iterrows():
    ot_v = str(r["variant_id"]).replace("-", "_")
    req = urllib.request.Request(
        ot_url,
        data=json.dumps({"query": q2, "variables": {"vid": ot_v}}).encode(),
        headers={"Content-Type": "application/json", "User-Agent": "A1Xv2Audit/1.0"}
    )
    try:
        with opener.open(req, timeout=12) as resp:
            vd = json.loads(resp.read().decode())["data"]["variant"]
        if vd and vd["credibleSets"]["count"] > 0:
            for row in vd["credibleSets"]["rows"]:
                apob_extra_hits.append({
                    "variant_id": r["variant_id"], "stars": int(r["stars"]),
                    "gnomad_het": int(r["het_joint"]), "hgvs": r["hgvs_name"],
                    "study_id": row["study"]["id"], "author": row["study"]["publicationFirstAuthor"],
                    "pmid": row["study"]["pubmedId"], "trait": row["study"]["traitFromSource"],
                    "n": row["study"]["nSamples"], "beta": row["beta"], "se": row["standardError"],
                    "p": f"{row['pValueMantissa']:.2f}e{row['pValueExponent']}",
                    "method": row["finemappingMethod"]
                })
    except Exception as e:
        pass
print(f"Found {len(apob_extra_hits)} Open Targets credible-set hits across APOB pLoFs (het>=5)!")
if apob_extra_hits:
    edf = pd.DataFrame(apob_extra_hits)
    print(edf[["variant_id", "stars", "gnomad_het", "study_id", "author", "beta", "se", "p", "trait"]].to_string(index=False))

with open(f"{v2}/ukbb_study_meta_and_harmonized_hits.json", "w") as f:
    json.dump({
        "study_meta": study_meta,
        "same_cohort_rows": same_cohort_rows,
        "apob_extra_hits": apob_extra_hits,
        "az_probe": az_probe
    }, f, indent=2)
print("Saved ukbb_study_meta_and_harmonized_hits.json")
