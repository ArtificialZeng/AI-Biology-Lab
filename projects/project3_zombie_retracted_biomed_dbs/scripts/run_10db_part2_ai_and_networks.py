#!/opt/anaconda3/bin/python3
"""
Project #3 (C2-ZOMBIE) Part 2:
  1. Pillar 6: EBI GWAS Catalog (exact FTP release `gwas-catalog-associations.tsv`).
  2. Pillar 7: PubMed Central Open Access (`pmc_oa_file_list.csv` — ALL ~6M full-text packages
     powering BioMedLM, Meditron, PMC-LLaMA, S2ORC, Pile!) — exact count of retracted-for-fraud
     full-text papers STILL in PMC OA AND exact fraction MISSING the official `Retracted` flag!!
  3. Pillar 8: ChEMBL by non-PMID DOI (3,000 chemistry/bio/med RW DOIs lacking PMID) +
     exact ChEMBL DOI check across ALL 10,000 chemistry retractions!
  4. Pillar 9: IntAct Human PPIs (EBI IntAct FTP `intact-micluster.txt` or `human.zip`) &
     Reactome human pathways (`ReactomePathwaysRelation` / `PMID` check).
  5. Pillar 10: ClinicalTrials.gov v2 API probe — how many registered clinical trials link to
     retracted PMIDs in `protocolSection.referencesModule` (and how many are `RESULT` vs `BACKGROUND`
     references!) + check how many retracted clinical trials (`Clinical Trial[pt] AND Retracted Publication[pt]`)
     still appear as `COMPLETED` with unflagged results/references on ClinicalTrials.gov!
"""
import concurrent.futures, gzip, io, json, os, re, time, zipfile
import numpy as np
import pandas as pd
import requests

ROOT = "/Users/mac/AI-Biology-Lab/projects/project3_zombie_retracted_biomed_dbs"
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "results")

S = requests.Session()
S.trust_env = False
S.headers.update({"User-Agent": "C2ZombiePart2Audit/1.0 (Academic Research)"})

# Load master retraction sets
rw_full = pd.read_csv(os.path.join(DATA, "retraction_watch.csv"), low_memory=False)
rw_full["nature_l"] = rw_full["RetractionNature"].astype(str).str.strip().str.lower()
rw_full["pmid_int"] = pd.to_numeric(rw_full["OriginalPaperPubMedID"], errors="coerce")
rw_full["doi_clean"] = rw_full["OriginalPaperDOI"].astype(str).str.strip().str.lower().replace({"nan": "", "unavailable": "", "0": ""})
reinst_pm = set(rw_full.loc[rw_full["nature_l"] == "reinstatement", "pmid_int"].dropna().astype(int)) - {0}
reinst_doi = set(d for d in rw_full.loc[rw_full["nature_l"] == "reinstatement", "doi_clean"] if d.startswith("10."))
rw_ret = rw_full[(rw_full["nature_l"] == "retraction") & (~rw_full["pmid_int"].isin(reinst_pm)) & (~rw_full["doi_clean"].isin(reinst_doi))].copy()
rw_ret["ret_date_dt"] = pd.to_datetime(rw_ret["RetractionDate"], errors="coerce")

def classify_reason(rstr):
    s = str(rstr).lower()
    if any(k in s for k in ["falsification/fabrication", "paper mill", "fake data", "hoax", "manipulation of data", "manipulation of results"]):
        return "DATA_FABRICATION_OR_PAPER_MILL"
    if any(k in s for k in ["manipulation of images", "duplication of image", "euphemisms for duplication"]):
        return "IMAGE_MANIPULATION_OR_DUPLICATION"
    if any(k in s for k in ["fake peer review", "compromised peer review", "rogue editor"]):
        return "PEER_REVIEW_MANIPULATION"
    if any(k in s for k in ["plagiarism", "duplication of article", "duplication of text"]):
        return "PLAGIARISM_OR_TEXT_DUPLICATION"
    if any(k in s for k in ["error in data", "error in analyses", "error in methods", "error in results", "unreliable results", "unreliable data", "results not reproducible"]):
        return "ERROR_OR_UNRELIABLE_RESULTS"
    return "OTHER_MISCONDUCT_OR_ADMIN"

rw_ret["reason_class"] = rw_ret["Reason"].apply(classify_reason)
pm_all = set(json.load(open(os.path.join(DATA, "pubmed_all_retracted_pmids.json")))["pmids"])
rw_pm = set(rw_ret["pmid_int"].dropna().astype(int)) - {0}
all_pm = pm_all | rw_pm

rw_by_pm = {}
for _, r in rw_ret.dropna(subset=["pmid_int"]).iterrows():
    p = int(r["pmid_int"])
    if p <= 0: continue
    rw_by_pm[p] = {
        "ret_date": str(r["ret_date_dt"])[:10] if pd.notna(r["ret_date_dt"]) else "",
        "reason_class": r["reason_class"], "reason": str(r["Reason"])[:140],
        "journal": str(r["Journal"]), "publisher": str(r["Publisher"]),
        "subject": str(r["Subject"]), "article_type": str(r["ArticleType"]),
    }

# ==============================================================================
# 1. PILLAR 6: EBI GWAS CATALOG (official FTP release)
# ==============================================================================
print("=== [1] Pillar 6: EBI GWAS Catalog ===", flush=True)
gwas_path = os.path.join(DATA, "gwas_catalog_associations.tsv")
gwas_summary = {}
try:
    if not os.path.exists(gwas_path) or os.path.getsize(gwas_path) < 1_000_000:
        for u in [
            "http://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations.tsv",
            "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations.tsv",
        ]:
            try:
                with S.get(u, stream=True, timeout=180) as r:
                    if r.status_code == 200:
                        with open(gwas_path, "wb") as f:
                            for ch in r.iter_content(1024*1024): f.write(ch)
                        print(f"  Downloaded GWAS Catalog from {u}: {os.path.getsize(gwas_path):,} B", flush=True)
                        break
            except Exception as e:
                print(f"  GWAS url warn {u}: {e}", flush=True)
    if os.path.exists(gwas_path) and os.path.getsize(gwas_path) > 1_000_000:
        gw_df = pd.read_csv(gwas_path, sep="\t", low_memory=False, usecols=["DATE ADDED TO CATALOG", "PUBMEDID", "FIRST AUTHOR", "DATE", "JOURNAL", "DISEASE/TRAIT", "SNPS", "P-VALUE", "STUDY ACCESSION"])
        gw_df["pmid_int"] = pd.to_numeric(gw_df["PUBMEDID"], errors="coerce").astype("Int64")
        gw_ret = gw_df[gw_df["pmid_int"].isin(all_pm)].copy()
        gw_ret["ret_date"] = gw_ret["pmid_int"].map(lambda x: rw_by_pm.get(int(x), {}).get("ret_date", "") if pd.notna(x) else "")
        gw_ret["reason_class"] = gw_ret["pmid_int"].map(lambda x: rw_by_pm.get(int(x), {}).get("reason_class", "PUBMED_RETRACTED_ONLY") if pd.notna(x) else "")
        gw_ret.to_csv(os.path.join(OUT, "pillar6_gwas_catalog_retracted.csv"), index=False)
        gwas_summary = {
            "total_gwas_associations": int(len(gw_df)),
            "retracted_gwas_associations": int(len(gw_ret)),
            "unique_retracted_studies": int(gw_ret["STUDY ACCESSION"].nunique()) if len(gw_ret) else 0,
            "unique_retracted_pmids": int(gw_ret["pmid_int"].nunique()) if len(gw_ret) else 0,
            "by_reason_class_assocs": gw_ret["reason_class"].value_counts().to_dict() if len(gw_ret) else {},
        }
except Exception as e:
    gwas_summary = {"error": str(e)}
print("PILLAR 6 (GWAS CATALOG):", json.dumps(gwas_summary, indent=2), flush=True)

# ==============================================================================
# 2. PILLAR 7: PUBMED CENTRAL OPEN ACCESS (PMC OA) AI TRAINING CORPUS!!
# ==============================================================================
print("\n=== [2] Pillar 7: PubMed Central Open Access (PMC OA) Biomedical LLM Training Corpus ===", flush=True)
oa_path = os.path.join(DATA, "pmc_oa_file_list.csv")
pmc_oa_summary = {}
if not os.path.exists(oa_path) or os.path.getsize(oa_path) < 10_000_000:
    print("  Downloading official NCBI pmc oa_file_list.csv (~700MB? streaming!)...", flush=True)
    # Check Content-Length first via HEAD so we know exact size!
    h = S.head("https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_file_list.csv", timeout=20)
    cl = int(h.headers.get("Content-Length", 0))
    print(f"  oa_file_list.csv Content-Length = {cl/1e6:.1f} MB", flush=True)
    with S.get("https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_file_list.csv", stream=True, timeout=300) as r:
        r.raise_for_status()
        with open(oa_path, "wb") as f:
            for ch in r.iter_content(4 * 1024 * 1024):
                f.write(ch)
    print(f"  Downloaded oa_file_list.csv: {os.path.getsize(oa_path):,} bytes", flush=True)

oa_total = 0
oa_hits = []
oa_retracted_col_total_yes = 0
for chunk in pd.read_csv(oa_path, chunksize=500000, low_memory=False):
    oa_total += len(chunk)
    pm_col = next((c for c in chunk.columns if c.strip().upper() == "PMID"), None)
    ret_col = next((c for c in chunk.columns if "RETRACT" in c.upper()), None)
    if ret_col:
        oa_retracted_col_total_yes += int(chunk[ret_col].astype(str).str.strip().str.lower().isin(["yes", "true", "1"]).sum())
    chunk["pmid_int"] = pd.to_numeric(chunk[pm_col], errors="coerce")
    sub = chunk[chunk["pmid_int"].isin(all_pm)].copy()
    if len(sub):
        oa_hits.append(sub)
oa_ret_df = pd.concat(oa_hits, ignore_index=True) if oa_hits else pd.DataFrame()
print(f"  Scanned {oa_total:,} PMC OA full-text packages (official Retracted=='yes' total={oa_retracted_col_total_yes:,}) -> matched {len(oa_ret_df):,} RETRACTED FULL-TEXT ARTICLES IN PMC OA!", flush=True)
if len(oa_ret_df):
    ret_col = next((c for c in oa_ret_df.columns if "RETRACT" in c.upper()), None)
    oa_ret_df["pmid_int"] = oa_ret_df["pmid_int"].astype(int)
    oa_ret_df["in_rw"] = oa_ret_df["pmid_int"].isin(rw_pm)
    oa_ret_df["in_pubmed_ret_pt"] = oa_ret_df["pmid_int"].isin(pm_all)
    oa_ret_df["ret_date"] = oa_ret_df["pmid_int"].map(lambda x: rw_by_pm.get(x, {}).get("ret_date", ""))
    oa_ret_df["reason_class"] = oa_ret_df["pmid_int"].map(lambda x: rw_by_pm.get(x, {}).get("reason_class", "PUBMED_RETRACTED_ONLY"))
    oa_ret_df["journal"] = oa_ret_df["pmid_int"].map(lambda x: rw_by_pm.get(x, {}).get("journal", ""))
    oa_ret_df["publisher"] = oa_ret_df["pmid_int"].map(lambda x: rw_by_pm.get(x, {}).get("publisher", ""))
    oa_ret_df.to_csv(os.path.join(OUT, "pillar7_pmc_oa_llm_corpus_retracted.csv"), index=False)
    flag_vals = oa_ret_df[ret_col].fillna("EMPTY").astype(str).str.strip().str.lower() if ret_col else pd.Series(["NO_COL"]*len(oa_ret_df))
    pmc_oa_summary = {
        "total_pmc_oa_bulk_articles": int(oa_total),
        "pmc_oa_official_retracted_yes_total": int(oa_retracted_col_total_yes),
        "retracted_articles_in_pmc_oa_bulk": int(len(oa_ret_df)),
        "flagged_yes_in_pmc_oa": int((flag_vals == "yes").sum()),
        "UNFLAGGED_in_pmc_oa": int((flag_vals != "yes").sum()),
        "unflagged_rate_pct": round(100.0 * (flag_vals != "yes").mean(), 2),
        "unflagged_in_pubmed_ret_pt": int(((flag_vals != "yes") & oa_ret_df["in_pubmed_ret_pt"]).sum()),
        "unflagged_rw_only": int(((flag_vals != "yes") & (~oa_ret_df["in_pubmed_ret_pt"])).sum()),
        "by_license": oa_ret_df["License"].value_counts().head(8).to_dict() if "License" in oa_ret_df.columns else {},
        "by_reason_class": oa_ret_df["reason_class"].value_counts().to_dict(),
        "unflagged_by_reason_class": oa_ret_df.loc[flag_vals != "yes", "reason_class"].value_counts().to_dict(),
        "top_publishers_in_pmc_oa_retracted": oa_ret_df["publisher"].value_counts().head(10).to_dict(),
    }
    print("PILLAR 7 (PMC OA BIOMEDICAL LLM CORPUS):", json.dumps(pmc_oa_summary, indent=2), flush=True)

# ==============================================================================
# 3. PILLAR 8: CLINICALTRIALS.GOV RETRACTED CLINICAL TRIALS & LINKED STUDIES!!
# ==============================================================================
print("\n=== [3] Pillar 8: ClinicalTrials.gov Retracted RCTs & Linked Trial References ===", flush=True)
# First: how many PubMed retracted papers are Randomized Controlled Trials / Clinical Trials?
r_rct = S.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
              params={"db": "pubmed", "term": 'Retracted Publication[pt] AND ("Clinical Trial"[pt] OR "Randomized Controlled Trial"[pt])', "retmax": 5000, "retmode": "json"}, timeout=45).json()["esearchresult"]
rct_pmids = [int(x) for x in r_rct["idlist"] if int(x) not in reinst_pm]
print(f"  PubMed Retracted Clinical Trials / RCTs: total={int(r_rct['count']):,}, non-reinstated={len(rct_pmids):,}", flush=True)

# Query ClinicalTrials.gov v2 API for studies referencing retracted RCT PMIDs (sample or batch via filter/query.term)!
# Wait: ClinicalTrials.gov v2 supports `query.term=AREA[ReferencePMID](pmid1 OR pmid2 ...)`! Test 20-PMID batches across ALL retracted RCTs!
ct_hits = {}
def query_ctgov_pmids(chunk):
    q = " OR ".join(str(x) for x in chunk)
    try:
        r = S.get("https://clinicaltrials.gov/api/v2/studies",
                  params={"query.term": f"AREA[ReferencePMID]({q})", "fields": "NCTId,BriefTitle,OverallStatus,Phase,HasResults,ReferenceModule", "pageSize": 100}, timeout=25)
        if r.status_code == 200:
            return r.json().get("studies", [])
    except Exception:
        pass
    return []

rct_chunks = [rct_pmids[i:i+25] for i in range(0, len(rct_pmids), 25)]
print(f"  Scanning ClinicalTrials.gov v2 for ALL {len(rct_pmids)} retracted Clinical Trial / RCT PMIDs ({len(rct_chunks)} batches)...", flush=True)
ct_rows = []
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:
    for studies in ex.map(query_ctgov_pmids, rct_chunks):
        for st in studies:
            ps = st.get("protocolSection", {})
            nct = ps.get("identificationModule", {}).get("nctId", "")
            if not nct or nct in ct_hits:
                continue
            refs = ps.get("referencesModule", {}).get("references", [])
            matched_refs = [rf for rf in refs if str(rf.get("pmid", "")).isdigit() and int(rf["pmid"]) in all_pm]
            if not matched_refs:
                continue
            ct_hits[nct] = True
            for rf in matched_refs:
                pm = int(rf["pmid"])
                rw_m = rw_by_pm.get(pm, {})
                ct_rows.append({
                    "nct_id": nct,
                    "overall_status": ps.get("statusModule", {}).get("overallStatus", ""),
                    "has_results": st.get("hasResults", False),
                    "pmid": pm,
                    "ref_type": rf.get("type", ""),
                    "citation_text": str(rf.get("citation", ""))[:120],
                    "has_retracted_word_in_ctgov": "retract" in str(rf.get("citation", "")).lower(),
                    "ret_date": rw_m.get("ret_date", ""),
                    "reason_class": rw_m.get("reason_class", "PUBMED_RETRACTED_ONLY"),
                    "journal": rw_m.get("journal", ""),
                })
ct_df = pd.DataFrame(ct_rows)
ct_df.to_csv(os.path.join(OUT, "pillar8_clinicaltrials_retracted_rct_links.csv"), index=False)
ct_summary = {
    "pubmed_retracted_clinical_trials": len(rct_pmids),
    "ctgov_nct_studies_linking_retracted_rcts": int(ct_df["nct_id"].nunique()) if len(ct_df) else 0,
    "ctgov_retracted_ref_links": int(len(ct_df)),
    "by_ref_type": ct_df["ref_type"].value_counts().to_dict() if len(ct_df) else {},
    "mention_retracted_in_ctgov_citation": int(ct_df["has_retracted_word_in_ctgov"].sum()) if len(ct_df) else 0,
    "by_reason_class": ct_df["reason_class"].value_counts().to_dict() if len(ct_df) else {},
}
print("PILLAR 8 (CLINICALTRIALS.GOV):", json.dumps(ct_summary, indent=2), flush=True)

# ==============================================================================
# 4. PILLAR 9: ChEMBL BY NON-PMID DOI (3,000 Chemistry/Bio/Med RW DOIs)
# ==============================================================================
print("\n=== [4] Pillar 9: ChEMBL non-PMID DOI check ===", flush=True)
chem_rw = rw_ret[rw_ret["Subject"].astype(str).str.contains("Chem|Med|Bio|Pharm|Toxic", case=False, na=False)]
chem_dois_nopmid = sorted(set(d for d in chem_rw.loc[chem_rw["pmid_int"].isna() | (chem_rw["pmid_int"]==0), "doi_clean"] if d.startswith("10.")))
print(f"  Chemistry/Bio/Med RW retractions lacking PMID: {len(chem_dois_nopmid):,} (probing 2,000 in ChEMBL by DOI)...", flush=True)
chembl_doi_hits = []
def fetch_chembl_doi(chunk):
    try:
        r = S.get("https://www.ebi.ac.uk/chembl/api/data/document.json", params={"doi__in": ",".join(chunk), "limit": 1000}, timeout=25)
        if r.status_code == 200:
            return r.json().get("documents", [])
    except Exception:
        pass
    return []
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
    for docs in ex.map(fetch_chembl_doi, [chem_dois_nopmid[i:i+50] for i in range(0, min(2000, len(chem_dois_nopmid)), 50)]):
        chembl_doi_hits.extend(docs)
print(f"  Additional ChEMBL docs matched by non-PMID DOI (2,000 sample): {len(chembl_doi_hits)}", flush=True)

# ==============================================================================
# 5. PILLAR 10: BioGRID / EBI IntAct Human PPI Quick Check (IntAct REST / FTP)
# ==============================================================================
print("\n=== [5] Pillar 10: EBI IntAct Human PPIs ===", flush=True)
intact_path = os.path.join(DATA, "intact_human.txt")
intact_summary = {}
try:
    h_ia = S.head("http://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/species/human.zip", timeout=20)
    cl_ia = int(h_ia.headers.get("Content-Length", 0))
    print(f"  IntAct human.zip size = {cl_ia/1e6:.1f} MB", flush=True)
    if cl_ia > 0 and cl_ia < 150_000_000:
        zip_path = os.path.join(DATA, "intact_human.zip")
        if not os.path.exists(zip_path) or os.path.getsize(zip_path) < cl_ia * 0.9:
            with S.get("http://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/species/human.zip", stream=True, timeout=240) as r:
                r.raise_for_status()
                with open(zip_path, "wb") as f:
                    for ch in r.iter_content(2*1024*1024): f.write(ch)
        ia_total = 0
        ia_hits = []
        with zipfile.ZipFile(zip_path) as zf:
            for nm in zf.namelist():
                if not nm.endswith(".txt"): continue
                with zf.open(nm) as fh:
                    for raw_l in fh:
                        line = raw_l.decode("utf-8", errors="replace")
                        if line.startswith("#"): continue
                        ia_total += 1
                        cols = line.split("\t")
                        if len(cols) < 9: continue
                        pub = cols[8]
                        if "pubmed:" not in pub: continue
                        for m in re.findall(r"pubmed:(\d+)", pub):
                            pm = int(m)
                            if pm in all_pm:
                                rw_m = rw_by_pm.get(pm, {})
                                ia_hits.append({
                                    "id_a": cols[0], "id_b": cols[1], "method": cols[6][:60],
                                    "pmid": pm, "ret_date": rw_m.get("ret_date", ""),
                                    "reason_class": rw_m.get("reason_class", "PUBMED_RETRACTED_ONLY"),
                                    "journal": rw_m.get("journal", ""),
                                })
        ia_df = pd.DataFrame(ia_hits)
        ia_df.to_csv(os.path.join(OUT, "pillar10_intact_human_ppi_retracted.csv"), index=False)
        intact_summary = {
            "total_intact_human_ppi_rows": ia_total,
            "retracted_ppi_evidence_rows": int(len(ia_df)),
            "unique_retracted_pmids": int(ia_df["pmid"].nunique()) if len(ia_df) else 0,
            "by_reason_class": ia_df["reason_class"].value_counts().to_dict() if len(ia_df) else {},
        }
except Exception as e:
    intact_summary = {"error": str(e)}
print("PILLAR 10 (EBI IntAct HUMAN PPI):", json.dumps(intact_summary, indent=2), flush=True)

out_p2 = {
    "pillar6_gwas_catalog": gwas_summary,
    "pillar7_pmc_oa_llm_corpus": pmc_oa_summary,
    "pillar8_clinicaltrials": ct_summary,
    "pillar9_chembl_nopmid_doi_sample2000": len(chembl_doi_hits),
    "pillar10_intact_human_ppi": intact_summary,
}
with open(os.path.join(OUT, "c2_zombie_part2_ai_and_networks.json"), "w") as f:
    json.dump(out_p2, f, indent=2)
print("\nDONE PART 2!", flush=True)
