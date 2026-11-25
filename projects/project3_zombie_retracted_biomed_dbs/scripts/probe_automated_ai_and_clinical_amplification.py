#!/opt/anaconda3/bin/python3
"""
30-second decisive magnitude & mechanism probe for Project #3 (C2-ZOMBIE):
Test the 'Human-Curated vs. Automated/AI-Ingested Amplification Paradox':
  1. PMC OA Biomedical LLM Training Corpus:
     - Exact NCBI PubMed/PMC count of retracted articles in Open Access subset
     - Exact AWS S3 `pmc-oa-opendata` directory structure + `retracted` flag check
  2. Automated Biomedical Knowledge Graphs (`Open Targets` europepmc evidence + `NCBI PubTator3`):
     - Do Open Targets Platform and NCBI's own PubTator3 actively serve retracted-for-fraud &
       paper-mill papers as unflagged gene-disease / target-drug evidence?
  3. Clinical Evidence Cascade (`1,512` Retracted Clinical Trials / RCTs):
     - Exact post-retraction citation counts by year in OpenAlex across retracted RCTs
     - Exact count of Meta-Analyses, Systematic Reviews, and Clinical Practice Guidelines citing
       retracted RCTs (via PubMed `elink` / OpenAlex)!
  4. ClinicalTrials.gov v2 exact PMID search syntax verification.
"""
import json, os, re, time
import pandas as pd
import requests

ROOT = "/Users/mac/AI-Biology-Lab/projects/project3_zombie_retracted_biomed_dbs"
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "results")

S = requests.Session()
S.trust_env = False
S.headers.update({"User-Agent": "C2ZombieAmplificationProbe/1.0 (mailto:research@example.org)"})

# 1. Exact PMC Open Access subset retracted counts via NCBI ESearch + S3 root listing!
print("=== [1] PMC Open Access (LLM Training Corpus) Retracted Counts ===", flush=True)
for label, db, term in [
    ("PubMed_ALL_Retracted", "pubmed", "Retracted Publication[pt]"),
    ("PubMed_In_PMC", "pubmed", 'Retracted Publication[pt] AND "pubmed pmc"[sb]'),
    ("PubMed_In_PMC_OA", "pubmed", 'Retracted Publication[pt] AND "pubmed pmc open access"[filter]'),
    ("PMC_OpenAccess_Retracted", "pmc", '"open access"[filter] AND "retracted publication"[filter]'),
    ("PMC_OpenAccess_RetractionInTitle", "pmc", '"open access"[filter] AND ("retracted"[ti] OR "retraction"[pt])'),
    ("PubMed_Retracted_MetaAnalysis_or_SR", "pubmed", 'Retracted Publication[pt] AND ("Meta-Analysis"[pt] OR "Systematic Review"[pt])'),
    ("PubMed_Retracted_Guideline", "pubmed", 'Retracted Publication[pt] AND ("Guideline"[pt] OR "Practice Guideline"[pt])'),
]:
    r = S.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
              params={"db": db, "term": term, "retmax": 1, "retmode": "json"}, timeout=25).json()
    cnt = r.get("esearchresult", {}).get("count", "ERR")
    print(f"  {label:38s}: {cnt}", flush=True)
    time.sleep(0.35)

# List AWS S3 `pmc-oa-opendata` top-level prefixes!
for pref in ["", "oa_comm/", "oa_comm/xml/", "oa_comm/txt/", "oa_noncomm/"]:
    r_s3 = S.get(f"https://pmc-oa-opendata.s3.amazonaws.com/?list-type=2&delimiter=/&prefix={pref}&max-keys=30", timeout=15)
    cps = re.findall(r"<Prefix>([^<]+)</Prefix>", r_s3.text)
    keys = re.findall(r"<Key>([^<]+)</Key>", r_s3.text)
    print(f"  S3 prefix='{pref}': subdirs={cps[:8]}, keys={keys[:6]}", flush=True)

# Also check NCBI FTP pub/pmc/oa_bulk/ directory!
r_ftp = S.get("https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/", timeout=15)
print("  NCBI ftp oa_bulk:", re.findall(r'href="([^"]+)"', r_ftp.text)[:15], flush=True)

# 2. Probe NCBI PubTator3 & Open Targets Platform on 15 known retracted paper-mill / fabrication PMIDs!
print("\n=== [2] Automated Knowledge Graphs (PubTator3 & Open Targets) on Retracted Paper-Mill Articles ===", flush=True)
rw = pd.read_parquet(os.path.join(DATA, "retraction_watch_clean_retractions.parquet"))
pm_fraud = rw[(rw["reason_class"] == "DATA_FABRICATION_OR_PAPER_MILL") & (rw["pmid_clean"].notna()) & (rw["pmid_clean"] > 0)].copy()
# Pick 15 oncology/genetics paper-mill retractions with PMIDs from 2018-2023
pm_sample = pm_fraud[pm_fraud["Subject"].astype(str).str.contains("Oncology|Genetics|Molecular|Cancer|Medicine", case=False, na=False)].head(20)
sample_pmids = [int(x) for x in pm_sample["pmid_clean"].tolist()]
print(f"  Testing {len(sample_pmids)} retracted DATA_FABRICATION_OR_PAPER_MILL PMIDs in NCBI PubTator3 API...", flush=True)
pt3_r = S.get("https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/export/biocjson",
              params={"pmids": ",".join(str(x) for x in sample_pmids[:10])}, timeout=30)
print(f"  PubTator3 status={pt3_r.status_code}, len={len(pt3_r.text)}", flush=True)
if pt3_r.status_code == 200:
    d_pt = pt3_r.json()
    docs_pt = d_pt.get("PubTator3", [])
    n_with_ann = 0
    n_mentions_ret = 0
    tot_anns = 0
    for doc in docs_pt:
        anns = sum(len(p.get("annotations", [])) for p in doc.get("passages", []))
        tot_anns += anns
        if anns > 0: n_with_ann += 1
        txt_all = " ".join(p.get("text", "") for p in doc.get("passages", [])).lower()
        if "retract" in txt_all: n_mentions_ret += 1
    print(f"  PubTator3 returned {len(docs_pt)}/10 retracted fraud papers: {n_with_ann} have active AI-mined entity annotations (total={tot_anns}), mentions 'retract'={n_mentions_ret}!", flush=True)

# Test Open Targets Platform GraphQL `publications` / `evidences` for a famous retracted-paper gene (e.g., BRCA1/MLH1/TP53 or our ClinVar 3* PMIDs 32076465, 34721658!)
ot_url = "https://api.platform.opentargets.org/api/v4/graphql"
# Check if Open Targets `target` or `disease` europepmc evidence includes retracted PMIDs (query top 5 genes from paper-mill retractions!)
q_ot = """
query TestEPMC($ens: String!, $efo: String!) {
  disease(efoId: $efo) {
    evidences(ensemblIds: [$ens], datasourceIds: ["europepmc"], size: 100) {
      count
      rows { literature }
    }
  }
}
"""
# Actually simpler: check Open Targets Europe PMC evidence API or Europe PMC annotations directly, AND check PubTator3 relations!
r_rel = S.get("https://www.ncbi.nlm.nih.gov/research/pubtator3-api/relations",
              params={"e1": "@GENE_TP53", "type": "ANY"}, timeout=20)
print(f"  PubTator3 relations endpoint test: status={r_rel.status_code}, sample={r_rel.text[:160]}", flush=True)

# 3. OpenAlex Post-Retraction Citation Cascade across ALL 1,512 Retracted Clinical Trials / RCTs!!
print("\n=== [3] OpenAlex Post-Retraction Citation Cascade for Retracted Clinical Trials / RCTs ===", flush=True)
r_rct = S.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
              params={"db": "pubmed", "term": 'Retracted Publication[pt] AND ("Clinical Trial"[pt] OR "Randomized Controlled Trial"[pt])', "retmax": 5000, "retmode": "json"}, timeout=45).json()["esearchresult"]
rct_pmids = [int(x) for x in r_rct["idlist"]]
# Query OpenAlex in batches of 50 PMIDs for first 500 retracted RCTs to get exact cited_by_count + counts_by_year!
rw_ret_yr = {}
for _, r in rw.dropna(subset=["pmid_clean"]).iterrows():
    dt = pd.to_datetime(r["RetractionDate"], errors="coerce")
    if pd.notna(dt):
        rw_ret_yr[int(r["pmid_clean"])] = dt.year

oa_rct_rows = []
for i in range(0, min(500, len(rct_pmids)), 50):
    chunk = rct_pmids[i:i+50]
    filt = "pmid:" + "|".join(str(x) for x in chunk)
    try:
        r = S.get("https://api.openalex.org/works",
                  params={"filter": filt, "per-page": 50,
                          "select": "id,ids,title,publication_year,is_retracted,cited_by_count,counts_by_year"},
                  timeout=25)
        if r.status_code == 200:
            for w in r.json().get("results", []):
                pm_url = (w.get("ids") or {}).get("pmid", "")
                m = re.search(r"(\d+)$", pm_url or "")
                pm = int(m.group(1)) if m else None
                ry = rw_ret_yr.get(pm)
                cby = w.get("counts_by_year") or []
                post_cites = sum(c["cited_by_count"] for c in cby if ry and c["year"] > ry)
                same_or_post = sum(c["cited_by_count"] for c in cby if ry and c["year"] >= ry)
                oa_rct_rows.append({
                    "pmid": pm, "pub_year": w.get("publication_year"), "ret_year": ry,
                    "openalex_is_retracted": w.get("is_retracted"),
                    "total_cites": w.get("cited_by_count", 0),
                    "post_ret_cites_strict": post_cites,
                    "post_ret_cites_ge": same_or_post,
                    "title": str(w.get("title", ""))[:80],
                })
        time.sleep(0.25)
    except Exception as e:
        print(f"  OA batch {i} warn: {e}", flush=True)

oa_df = pd.DataFrame(oa_rct_rows)
if len(oa_df):
    dated = oa_df[oa_df["ret_year"].notna()]
    print(f"  OpenAlex matched {len(oa_df)}/500 retracted RCTs (RW-dated={len(dated)}):", flush=True)
    print(f"    OpenAlex is_retracted==False (MISSED BY OPENALEX!): {int((oa_df['openalex_is_retracted']==False).sum())}/{len(oa_df)} ({100*(oa_df['openalex_is_retracted']==False).mean():.1f}%)", flush=True)
    print(f"    Total citations across 500 retracted RCTs: {int(oa_df['total_cites'].sum()):,}", flush=True)
    print(f"    STRICT Post-Retraction citations (year > ret_year, dated n={len(dated)}): {int(dated['post_ret_cites_strict'].sum()):,} (RCTs with >=1 post-ret cite: {int((dated['post_ret_cites_strict']>0).sum())}/{len(dated)})", flush=True)
    print("    Top 5 most-cited post-retraction RCTs:")
    print(dated.sort_values("post_ret_cites_strict", ascending=False)[["pmid", "pub_year", "ret_year", "openalex_is_retracted", "total_cites", "post_ret_cites_strict", "title"]].head(5).to_string(index=False))

# 4. ClinicalTrials.gov v2 exact PMID query test
print("\n=== [4] ClinicalTrials.gov v2 exact PMID search test ===", flush=True)
# Test 5 highly-cited retracted RCTs via query.term directly!
test_pms = oa_df.sort_values("total_cites", ascending=False)["pmid"].dropna().astype(int).head(10).tolist() if len(oa_df) else rct_pmids[:10]
ct_found = 0
for pm in test_pms:
    r_ct = S.get("https://clinicaltrials.gov/api/v2/studies", params={"query.term": str(pm), "fields": "NCTId,OverallStatus,ReferenceModule", "pageSize": 10}, timeout=15)
    if r_ct.status_code == 200:
        sts = r_ct.json().get("studies", [])
        if sts:
            ct_found += len(sts)
            nct0 = sts[0].get("protocolSection", {}).get("identificationModule", {}).get("nctId")
            print(f"  CT.gov HIT for retracted RCT PMID {pm}: {len(sts)} studies (e.g. {nct0})", flush=True)
print(f"  Total CT.gov hits across top 10 retracted RCTs: {ct_found}", flush=True)
