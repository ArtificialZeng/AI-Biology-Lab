#!/opt/anaconda3/bin/python3
"""
Project #3 (C2-ZOMBIE): Persistent Retracted-for-Fraud Literature & 'Database Laundering'
Across Production Biomedical Knowledge Bases & AI Training Benchmarks
(ClinVar, ChEMBL, RCSB PDB, UniProtKB/Swiss-Prot).

Step 1: Ingest complete Crossref Retraction Watch official database (CSV) + NCBI PubMed
Retracted Publication[pt] PMIDs and build a clean, deduplicated retracted-literature master
registry classified by retraction mechanism (Fraud/Fabrication/PaperMill vs ImageDuplication
vs HonestError/Reproducibility vs Plagiarism/PeerReviewManipulation).
"""
import hashlib, json, os, time, urllib.request
import pandas as pd

ROOT = "/Users/mac/AI-Biology-Lab/projects/project3_zombie_retracted_biomed_dbs"
DATA = os.path.join(ROOT, "data")
os.makedirs(DATA, exist_ok=True)

opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
UA = {"User-Agent": "C2ZombieAudit/1.0 (Academic Research)"}

rw_csv = os.path.join(DATA, "retraction_watch.csv")
rw_url = "https://gitlab.com/crossref/retraction-watch-data/-/raw/main/retraction_watch.csv"

if not os.path.exists(rw_csv) or os.path.getsize(rw_csv) < 10_000_000:
    print(f"Downloading Crossref Retraction Watch CSV from {rw_url}...")
    req = urllib.request.Request(rw_url, headers=UA)
    with opener.open(req, timeout=180) as r, open(rw_csv, "wb") as f:
        total = 0
        while True:
            chunk = r.read(1024 * 1024)
            if not chunk:
                break
            f.write(chunk)
            total += len(chunk)
    print(f"Downloaded {total:,} bytes to {rw_csv}")

# Compute SHA-256 for reproducibility
h = hashlib.sha256()
with open(rw_csv, "rb") as f:
    for blk in iter(lambda: f.read(1024 * 1024), b""):
        h.update(blk)
rw_sha = h.hexdigest()

df = pd.read_csv(rw_csv, low_memory=False)
print(f"Retraction Watch rows: {len(df):,}, cols: {list(df.columns)}")
print("RetractionNature counts:\n", df["RetractionNature"].value_counts().head(10).to_string())

# Filter strictly to true Retractions (exclude pure Expressions of Concern / Corrections unless also retracted)
ret = df[df["RetractionNature"].astype(str).str.strip().str.lower() == "retraction"].copy()
ret["pmid_clean"] = pd.to_numeric(ret["OriginalPaperPubMedID"], errors="coerce").astype("Int64")
ret["doi_clean"] = ret["OriginalPaperDOI"].astype(str).str.strip().str.lower().replace({"nan": "", "unavailable": "", "0": ""})

# Classify retraction reasons
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

ret["reason_class"] = ret["Reason"].apply(classify_reason)
ret["orig_date"] = pd.to_datetime(ret["OriginalPaperDate"], errors="coerce")
ret["ret_date"] = pd.to_datetime(ret["RetractionDate"], errors="coerce")
ret["lag_days"] = (ret["ret_date"] - ret["orig_date"]).dt.days

# Also fetch ALL PubMed Retracted Publication[pt] PMIDs via ESearch webenv/chunking
pmids_pubmed = set()
esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=Retracted+Publication%5Bpt%5D&retmax=100000&retmode=json"
try:
    with opener.open(urllib.request.Request(esearch_url, headers=UA), timeout=60) as r:
        ed = json.loads(r.read().decode())["esearchresult"]
        cnt = int(ed["count"])
        pmids_pubmed.update(int(x) for x in ed["idlist"])
        print(f"PubMed Retracted Publication[pt] total={cnt:,}, fetched={len(pmids_pubmed):,}")
except Exception as e:
    print(f"PubMed esearch warn: {e}")

rw_pmids = set(ret["pmid_clean"].dropna().astype(int).tolist()) - {0}
rw_dois = set(d for d in ret["doi_clean"].tolist() if d.startswith("10."))
union_pmids = rw_pmids | pmids_pubmed

ret_slim = ret[[
    "Record ID", "Title", "Subject", "Journal", "Publisher", "Country",
    "OriginalPaperDate", "RetractionDate", "lag_days", "pmid_clean", "doi_clean",
    "RetractionDOI", "Reason", "reason_class", "Paywalled"
]].copy()
ret_slim.to_parquet(os.path.join(DATA, "retraction_watch_clean_retractions.parquet"), index=False)

with open(os.path.join(DATA, "pubmed_retracted_pmids.json"), "w") as f:
    json.dump({"pubmed_retracted_pt_count": len(pmids_pubmed), "rw_retraction_pmid_count": len(rw_pmids), "union_pmid_count": len(union_pmids), "pmids": sorted(union_pmids)}, f)

summary = {
    "rw_csv_sha256": rw_sha,
    "rw_total_rows": int(len(df)),
    "rw_true_retractions": int(len(ret)),
    "rw_unique_valid_pmids": int(len(rw_pmids)),
    "rw_unique_valid_dois": int(len(rw_dois)),
    "pubmed_retracted_pt_pmids": int(len(pmids_pubmed)),
    "union_unique_retracted_pmids": int(len(union_pmids)),
    "reason_class_counts": ret["reason_class"].value_counts().to_dict(),
    "median_retraction_lag_days": float(ret["lag_days"].median()) if ret["lag_days"].notna().any() else None,
}
with open(os.path.join(DATA, "retraction_master_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print("SUMMARY:", json.dumps(summary, indent=2))
