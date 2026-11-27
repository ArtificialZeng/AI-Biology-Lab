#!/opt/anaconda3/bin/python3
"""
Project #3 (C2-ZOMBIE): Persistent Retracted Literature & 'Database Laundering'
Across Production Biomedical Knowledge Bases & AI Training Infrastructure
(ClinVar, RCSB PDB, UniProtKB/Swiss-Prot, ChEMBL).

Step 2:
  0. Complete full PubMed Retracted Publication[pt] pagination (all ~35,847 PMIDs) +
     strictly exclude all 164 Retraction Watch Reinstatements -> airtight master registry!
  1. Pillar 1 (ClinVar): Download official NCBI ClinVar `var_citations.txt` + reuse/fetch
     `variant_summary.txt.gz` -> exact count of ACTIVE ClinVar variants (esp. Pathogenic /
     Likely pathogenic / 2-4*) citing retracted papers, sole-PMID dependence, ACMG SF overlap,
     and post-retraction persistence years.
  2. Pillar 2 (RCSB PDB): Query live RCSB Search API across ALL retracted PMIDs & DOIs ->
     exact count of ACTIVE, NON-OBSOLETE experimental 3D structures from retracted papers,
     AlphaFold2/3 & Boltz training-cutoff membership, and structure validation metrics!
  3. Pillar 3 (UniProtKB/Swiss-Prot): Query live UniProt REST API for reviewed human Swiss-Prot
     entries citing retracted PMIDs -> exact count of active curated human proteins & annotations!
  4. Pillar 4 (ChEMBL): Query live EBI ChEMBL REST API `/document.json` across retracted PMIDs/DOIs ->
     exact count of retracted medicinal chemistry papers still active in ChEMBL + total active
     bioactivity records (IC50/Ki/pChEMBL)!
"""
import concurrent.futures, gzip, json, os, time, urllib.parse
import numpy as np
import pandas as pd
import requests

ROOT = "/Users/mac/AI-Biology-Lab/projects/project3_zombie_retracted_biomed_dbs"
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "results")
os.makedirs(DATA, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

S = requests.Session()
S.trust_env = False  # bypass local proxy env if any
S.headers.update({"User-Agent": "C2ZombieBiomedAudit/1.0 (Academic Research)"})

# ==============================================================================
# STEP 0: Airtight Master Retraction Registry (RW + ALL ~35,847 PubMed PMIDs - Reinstatements)
# ==============================================================================
print("=== [STEP 0] Building Airtight Master Retraction Registry ===", flush=True)
rw_full = pd.read_csv(os.path.join(DATA, "retraction_watch.csv"), low_memory=False)
rw_full["nature_l"] = rw_full["RetractionNature"].astype(str).str.strip().str.lower()
rw_full["pmid_int"] = pd.to_numeric(rw_full["OriginalPaperPubMedID"], errors="coerce")
rw_full["doi_clean"] = rw_full["OriginalPaperDOI"].astype(str).str.strip().str.lower().replace({"nan": "", "unavailable": "", "0": ""})

reinstated_pmids = set(rw_full.loc[rw_full["nature_l"] == "reinstatement", "pmid_int"].dropna().astype(int)) - {0}
reinstated_dois = set(d for d in rw_full.loc[rw_full["nature_l"] == "reinstatement", "doi_clean"] if d.startswith("10."))
print(f"Excluding {len(reinstated_pmids)} reinstated PMIDs and {len(reinstated_dois)} reinstated DOIs.", flush=True)

rw_ret = rw_full[rw_full["nature_l"] == "retraction"].copy()
rw_ret = rw_ret[~rw_ret["pmid_int"].isin(reinstated_pmids)]
rw_ret = rw_ret[~rw_ret["doi_clean"].isin(reinstated_dois)]

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
rw_ret["orig_date"] = pd.to_datetime(rw_ret["OriginalPaperDate"], errors="coerce")
rw_ret["ret_date"] = pd.to_datetime(rw_ret["RetractionDate"], errors="coerce")

# Fetch ALL ~35,847 PubMed Retracted Publication[pt] PMIDs using date-range slices (NCBI ESearch caps retstart<10,000!)
pm_cache = os.path.join(DATA, "pubmed_all_retracted_pmids.json")
if os.path.exists(pm_cache):
    pm_pmids = set(json.load(open(pm_cache))["pmids"])
else:
    pm_pmids = set()
    r_tot = S.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                  params={"db": "pubmed", "term": "Retracted Publication[pt]", "retmax": 1, "retmode": "json"}, timeout=60).json()["esearchresult"]
    total_pm = int(r_tot["count"])
    # Partition by pdat (publication date) into <9,000-record slices!
    slices = [
        "1900/01/01:2012/12/31[pdat]",
        "2013/01/01:2017/12/31[pdat]",
        "2018/01/01:2020/12/31[pdat]",
        "2021/01/01:2022/06/30[pdat]",
        "2022/07/01:2023/06/30[pdat]",
        "2023/07/01:2024/12/31[pdat]",
        "2025/01/01:2027/12/31[pdat]",
    ]
    for sl in slices:
        time.sleep(0.4)
        ri = S.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                   params={"db": "pubmed", "term": f"Retracted Publication[pt] AND ({sl})", "retmax": 9999, "retmode": "json"}, timeout=60).json()["esearchresult"]
        c_sl = int(ri["count"])
        ids = [int(x) for x in ri["idlist"]]
        pm_pmids.update(ids)
        print(f"  PubMed slice {sl}: count={c_sl:,}, fetched={len(ids):,}", flush=True)
    pm_pmids -= reinstated_pmids
    with open(pm_cache, "w") as f:
        json.dump({"pubmed_total_reported": total_pm, "fetched_non_reinstated": len(pm_pmids), "pmids": sorted(pm_pmids)}, f)
    print(f"Complete PubMed Retracted Publication[pt]: total_reported={total_pm:,}, fetched_non_reinstated={len(pm_pmids):,}", flush=True)

rw_pmids = set(rw_ret["pmid_int"].dropna().astype(int)) - {0}
rw_dois = set(d for d in rw_ret["doi_clean"] if d.startswith("10."))
all_ret_pmids = rw_pmids | pm_pmids
print(f"MASTER RETRACTION SET: RW_PMIDs={len(rw_pmids):,}, PubMed_PMIDs={len(pm_pmids):,}, UNION_PMIDs={len(all_ret_pmids):,}, RW_DOIs={len(rw_dois):,}", flush=True)

# Build fast PMID -> metadata lookup (from Retraction Watch where available)
rw_by_pmid = {}
for _, r in rw_ret.dropna(subset=["pmid_int"]).iterrows():
    p = int(r["pmid_int"])
    if p == 0: continue
    rw_by_pmid[p] = {
        "pmid": p, "doi": r["doi_clean"], "journal": str(r["Journal"]),
        "publisher": str(r["Publisher"]), "subject": str(r["Subject"]),
        "orig_date": str(r["orig_date"])[:10] if pd.notna(r["orig_date"]) else "",
        "ret_date": str(r["ret_date"])[:10] if pd.notna(r["ret_date"]) else "",
        "reason": str(r["Reason"]), "reason_class": r["reason_class"],
        "country": str(r["Country"]),
    }

# ==============================================================================
# PILLAR 1: NCBI CLINVAR (`var_citations.txt` + `variant_summary.txt.gz`)
# ==============================================================================
print("\n=== [PILLAR 1] Auditing NCBI ClinVar Active Variant Citations ===", flush=True)
vc_path = os.path.join(DATA, "clinvar_var_citations.txt")
if not os.path.exists(vc_path) or os.path.getsize(vc_path) < 1_000_000:
    print("Downloading ClinVar var_citations.txt...", flush=True)
    with S.get("https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/var_citations.txt", stream=True, timeout=180) as r:
        r.raise_for_status()
        with open(vc_path, "wb") as f:
            for ch in r.iter_content(1024 * 1024):
                f.write(ch)
    print(f"Downloaded var_citations.txt: {os.path.getsize(vc_path):,} bytes", flush=True)

vc_df = pd.read_csv(vc_path, sep="\t", low_memory=False)
print(f"ClinVar var_citations rows: {len(vc_df):,}, cols: {list(vc_df.columns)}", flush=True)
vc_pm = vc_df[vc_df["citation_source"].astype(str) == "PubMed"].copy()
vc_pm["pmid_int"] = pd.to_numeric(vc_pm["citation_id"], errors="coerce")
vc_pm = vc_pm.dropna(subset=["pmid_int"])
vc_pm["pmid_int"] = vc_pm["pmid_int"].astype(int)

# Total unique PubMed citations per VariationID (for sole-PMID dependence!)
pmids_per_var = vc_pm.groupby("VariationID")["pmid_int"].nunique()
ret_vc = vc_pm[vc_pm["pmid_int"].isin(all_ret_pmids)].copy()
ret_vc["total_pmids_for_var"] = ret_vc["VariationID"].map(pmids_per_var)
ret_vc["sole_pmid"] = ret_vc["total_pmids_for_var"] == 1
ret_vc["reason_class"] = ret_vc["pmid_int"].map(lambda x: rw_by_pmid.get(x, {}).get("reason_class", "PUBMED_RETRACTED_ONLY"))
ret_vc["ret_date"] = ret_vc["pmid_int"].map(lambda x: rw_by_pmid.get(x, {}).get("ret_date", ""))
ret_vc["journal"] = ret_vc["pmid_int"].map(lambda x: rw_by_pmid.get(x, {}).get("journal", ""))

# Join to existing local ClinVar variant_summary.txt.gz (GRCh38) to get exact ClinicalSignificance, ReviewStatus, GeneSymbol!
vs_path = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_a1_clinvar_gnomad_reclass/data/raw/variant_summary.txt.gz"
ret_vids = set(ret_vc["VariationID"].unique())
vs_rows = []
if os.path.exists(vs_path):
    print(f"Joining {len(ret_vids):,} retracted-citing ClinVar VariationIDs against variant_summary.txt.gz...", flush=True)
    for chunk in pd.read_csv(vs_path, sep="\t", compression="gzip", chunksize=250000, low_memory=False,
                             usecols=["VariationID", "Assembly", "GeneSymbol", "ClinicalSignificance", "ReviewStatus", "Type", "Name", "LastEvaluated"]):
        sub = chunk[(chunk["Assembly"] == "GRCh38") & (chunk["VariationID"].isin(ret_vids))]
        if len(sub):
            vs_rows.append(sub)
    vs_m = pd.concat(vs_rows, ignore_index=True).drop_duplicates("VariationID") if vs_rows else pd.DataFrame()
    ret_vc = ret_vc.merge(vs_m, on="VariationID", how="left")
    ret_vc.to_csv(os.path.join(OUT, "pillar1_clinvar_retracted_citations.csv"), index=False)

is_plp_cv = ret_vc["ClinicalSignificance"].astype(str).apply(
    lambda s: ("Conflicting" not in s) and (("Pathogenic" in s) or ("Likely pathogenic" in s))
) if "ClinicalSignificance" in ret_vc.columns else pd.Series([False]*len(ret_vc))

clinvar_summary = {
    "total_clinvar_pubmed_citation_links": int(len(vc_pm)),
    "retracted_citation_links": int(len(ret_vc)),
    "unique_retracted_pmids_in_clinvar": int(ret_vc["pmid_int"].nunique()),
    "unique_variants_citing_retracted": int(ret_vc["VariationID"].nunique()),
    "unique_plp_variants_citing_retracted": int(ret_vc.loc[is_plp_cv, "VariationID"].nunique()),
    "unique_plp_sole_pmid_retracted": int(ret_vc.loc[is_plp_cv & ret_vc["sole_pmid"], "VariationID"].nunique()),
    "plp_by_review_status": ret_vc.loc[is_plp_cv].drop_duplicates("VariationID")["ReviewStatus"].value_counts().to_dict() if "ReviewStatus" in ret_vc.columns else {},
    "plp_by_reason_class": ret_vc.loc[is_plp_cv].drop_duplicates("VariationID")["reason_class"].value_counts().to_dict(),
    "top_genes_plp_retracted": ret_vc.loc[is_plp_cv].drop_duplicates("VariationID")["GeneSymbol"].value_counts().head(15).to_dict() if "GeneSymbol" in ret_vc.columns else {},
}
print("PILLAR 1 (CLINVAR):", json.dumps(clinvar_summary, indent=2), flush=True)

# ==============================================================================
# PILLAR 2: RCSB PDB ACTIVE (NON-OBSOLETE) 3D STRUCTURES FROM RETRACTED PAPERS
# ==============================================================================
print("\n=== [PILLAR 2] Auditing Active (Non-Obsolete) RCSB PDB Structures ===", flush=True)
pdb_url = "https://search.rcsb.org/rcsbsearch/v2/query"
pmid_list = sorted(all_ret_pmids)
pdb_hits = set()
# Batch query RCSB Search API in chunks of 1,500 PMIDs (`in` operator!)
for i in range(0, len(pmid_list), 1500):
    chunk = pmid_list[i:i+1500]
    payload = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_primary_citation.pdbx_database_id_PubMed",
                "operator": "in",
                "value": chunk
            }
        },
        "return_type": "entry",
        "request_options": {"return_all_hits": True}
    }
    try:
        r = S.post(pdb_url, json=payload, timeout=30)
        if r.status_code == 200:
            d = r.json()
            for it in d.get("result_set", []):
                pdb_hits.add(it["identifier"])
    except Exception as e:
        print(f"  PDB chunk {i} warn: {e}", flush=True)

# Also query by Retraction Watch DOIs (in chunks of 1,000 for chemistry/structural papers lacking PMID!)
doi_list = sorted(rw_dois)
for i in range(0, len(doi_list), 1000):
    chunk = doi_list[i:i+1000]
    payload = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_primary_citation.pdbx_database_id_DOI",
                "operator": "in",
                "value": chunk
            }
        },
        "return_type": "entry",
        "request_options": {"return_all_hits": True}
    }
    try:
        r = S.post(pdb_url, json=payload, timeout=30)
        if r.status_code == 200:
            for it in r.json().get("result_set", []):
                pdb_hits.add(it["identifier"])
    except Exception as e:
        pass

print(f"Found {len(pdb_hits)} ACTIVE (NON-OBSOLETE) RCSB PDB entries linked to retracted papers!", flush=True)
pdb_details = []
if pdb_hits:
    # Fetch full metadata + validation metrics via RCSB Data GraphQL API!
    gql = """
    query($ids: [String!]!) {
      entries(entry_ids: $ids) {
        rcsb_id
        struct { title }
        rcsb_accession_info { deposit_date initial_release_date revision_date }
        exptl { method }
        rcsb_entry_info { resolution_combined deposited_atom_count polymer_entity_count_protein }
        refine { ls_R_factor_R_free ls_R_factor_R_work }
        rcsb_primary_citation {
          pdbx_database_id_PubMed
          pdbx_database_id_DOI
          journal_abbrev
          year
          title
        }
      }
    }
    """
    hit_list = sorted(pdb_hits)
    for i in range(0, len(hit_list), 100):
        r = S.post("https://data.rcsb.org/graphql", json={"query": gql, "variables": {"ids": hit_list[i:i+100]}}, timeout=30)
        if r.status_code == 200:
            for e in r.json().get("data", {}).get("entries", []) or []:
                pc = e.get("rcsb_primary_citation") or {}
                pm = pc.get("pdbx_database_id_PubMed")
                doi = str(pc.get("pdbx_database_id_DOI") or "").lower()
                rw_m = rw_by_pmid.get(int(pm), {}) if pm else {}
                acc = e.get("rcsb_accession_info") or {}
                rel = str(acc.get("initial_release_date") or "")[:10]
                ref = (e.get("refine") or [{}])[0] if e.get("refine") else {}
                ei = e.get("rcsb_entry_info") or {}
                res_v = (ei.get("resolution_combined") or [None])[0]
                pdb_details.append({
                    "pdb_id": e["rcsb_id"],
                    "method": (e.get("exptl") or [{}])[0].get("method", ""),
                    "resolution": res_v,
                    "r_free": ref.get("ls_R_factor_R_free"),
                    "release_date": rel,
                    "pre_af2_2018": rel <= "2018-04-30" if rel else False,
                    "pre_af3_2021": rel <= "2021-09-30" if rel else False,
                    "pmid": pm, "doi": doi,
                    "journal": pc.get("journal_abbrev", ""),
                    "pub_year": pc.get("year"),
                    "ret_date": rw_m.get("ret_date", ""),
                    "reason_class": rw_m.get("reason_class", "PUBMED_OR_DOI_RETRACTED"),
                    "reason": rw_m.get("reason", ""),
                    "title": (e.get("struct") or {}).get("title", "")[:100],
                })
    pdb_df = pd.DataFrame(pdb_details)
    pdb_df.to_csv(os.path.join(OUT, "pillar2_rcsb_pdb_active_retracted.csv"), index=False)
else:
    pdb_df = pd.DataFrame()

pdb_summary = {
    "active_non_obsolete_pdb_entries": int(len(pdb_df)),
    "unique_retracted_papers": int(pdb_df["pmid"].nunique()) if len(pdb_df) else 0,
    "in_af2_training_cutoff_le_2018_04": int(pdb_df["pre_af2_2018"].sum()) if len(pdb_df) else 0,
    "in_af3_boltz_training_cutoff_le_2021_09": int(pdb_df["pre_af3_2021"].sum()) if len(pdb_df) else 0,
    "by_reason_class": pdb_df["reason_class"].value_counts().to_dict() if len(pdb_df) else {},
    "by_method": pdb_df["method"].value_counts().to_dict() if len(pdb_df) else {},
}
print("PILLAR 2 (RCSB PDB):", json.dumps(pdb_summary, indent=2), flush=True)

# ==============================================================================
# PILLAR 3: UNIPROTKB / SWISS-PROT REVIEWED HUMAN PROTEOME CITING RETRACTED PMIDS
# ==============================================================================
print("\n=== [PILLAR 3] Auditing Reviewed Human UniProtKB/Swiss-Prot ===", flush=True)
# Download tiny reviewed human Swiss-Prot accession+PMID table in one fast streaming TSV!
up_url = "https://rest.uniprot.org/uniprotkb/stream?query=(reviewed:true)+AND+(organism_id:9606)&fields=accession,id,gene_primary,lit_pubmed_id,cc_disease&format=tsv"
r_up = S.get(up_url, timeout=120)
r_up.raise_for_status()
up_lines = r_up.text.splitlines()
print(f"Downloaded reviewed human Swiss-Prot table: {len(up_lines)-1:,} curated human proteins!", flush=True)
up_ret_rows = []
for line in up_lines[1:]:
    cols = line.split("\t")
    if len(cols) < 4: continue
    acc, entry_name, gene, pm_str = cols[0], cols[1], cols[2], cols[3]
    dis = cols[4] if len(cols) > 4 else ""
    if not pm_str: continue
    pms = [int(x.strip()) for x in pm_str.split(";") if x.strip().isdigit()]
    hits = [p for p in pms if p in all_ret_pmids]
    for h in hits:
        rw_m = rw_by_pmid.get(h, {})
        up_ret_rows.append({
            "uniprot_acc": acc, "entry_name": entry_name, "gene": gene,
            "pmid": h, "total_pmids_on_protein": len(pms),
            "has_disease_cc": bool(dis.strip()),
            "ret_date": rw_m.get("ret_date", ""),
            "reason_class": rw_m.get("reason_class", "PUBMED_RETRACTED_ONLY"),
            "journal": rw_m.get("journal", ""),
        })
up_df = pd.DataFrame(up_ret_rows)
up_df.to_csv(os.path.join(OUT, "pillar3_uniprot_human_swissprot_retracted.csv"), index=False)
uniprot_summary = {
    "total_reviewed_human_proteins": len(up_lines) - 1,
    "retracted_citation_links_in_human_swissprot": int(len(up_df)),
    "unique_reviewed_human_proteins_citing_retracted": int(up_df["uniprot_acc"].nunique()) if len(up_df) else 0,
    "unique_retracted_pmids_in_human_swissprot": int(up_df["pmid"].nunique()) if len(up_df) else 0,
    "by_reason_class": up_df["reason_class"].value_counts().to_dict() if len(up_df) else {},
}
print("PILLAR 3 (UNIPROT SWISS-PROT):", json.dumps(uniprot_summary, indent=2), flush=True)

# ==============================================================================
# PILLAR 4: ChEMBL PRODUCTION BIOACTIVITY DATABASE (LIVE EBI API)
# ==============================================================================
print("\n=== [PILLAR 4] Auditing ChEMBL Production Medicinal Chemistry Database ===", flush=True)
# Query ChEMBL /document.json in parallel batches of 100 PMIDs (filtering first to biomedical/chemistry RW + PubMed PMIDs)
# Wait: ChEMBL supports pubmed_id__in=... with up to 150 PMIDs per GET! Let's scan ALL ~42k retracted PMIDs in 8 parallel workers (~40s!)
chembl_docs = []
def fetch_chembl_batch(chunk):
    url = "https://www.ebi.ac.uk/chembl/api/data/document.json"
    params = {"pubmed_id__in": ",".join(str(x) for x in chunk), "limit": 1000}
    for att in range(3):
        try:
            r = S.get(url, params=params, timeout=30)
            if r.status_code == 200:
                return r.json().get("documents", [])
        except Exception:
            time.sleep(1.0)
    return []

chunks = [pmid_list[i:i+150] for i in range(0, len(pmid_list), 150)]
print(f"Querying ChEMBL /document.json across {len(pmid_list):,} retracted PMIDs ({len(chunks)} batches, 8 workers)...", flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
    for docs in ex.map(fetch_chembl_batch, chunks):
        chembl_docs.extend(docs)

print(f"Found {len(chembl_docs):,} ACTIVE ChEMBL documents matching retracted PMIDs!", flush=True)
# For every matched ChEMBL document, fetch exact count of active bioactivity records!
def fetch_doc_act_count(doc):
    did = doc["document_chembl_id"]
    pm = doc.get("pubmed_id")
    rw_m = rw_by_pmid.get(int(pm), {}) if pm else {}
    n_act, n_pchembl = 0, 0
    try:
        r1 = S.get("https://www.ebi.ac.uk/chembl/api/data/activity.json", params={"document_chembl_id": did, "limit": 1}, timeout=20)
        if r1.status_code == 200:
            n_act = r1.json().get("page_meta", {}).get("total_count", 0)
        r2 = S.get("https://www.ebi.ac.uk/chembl/api/data/activity.json", params={"document_chembl_id": did, "pchembl_value__isnull": "false", "limit": 1}, timeout=20)
        if r2.status_code == 200:
            n_pchembl = r2.json().get("page_meta", {}).get("total_count", 0)
    except Exception:
        pass
    return {
        "document_chembl_id": did, "pmid": pm, "doi": doc.get("doi", ""),
        "journal": doc.get("journal", ""), "year": doc.get("year"),
        "title": str(doc.get("title", ""))[:100],
        "n_activities": n_act, "n_pchembl_activities": n_pchembl,
        "ret_date": rw_m.get("ret_date", ""),
        "reason_class": rw_m.get("reason_class", "PUBMED_RETRACTED_ONLY"),
        "reason": rw_m.get("reason", ""),
    }

chembl_enriched = []
if chembl_docs:
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
        for row in ex.map(fetch_doc_act_count, chembl_docs):
            chembl_enriched.append(row)
ch_df = pd.DataFrame(chembl_enriched)
ch_df.to_csv(os.path.join(OUT, "pillar4_chembl_retracted_docs_and_activities.csv"), index=False)

chembl_summary = {
    "active_retracted_chembl_documents": int(len(ch_df)),
    "total_active_bioactivity_records": int(ch_df["n_activities"].sum()) if len(ch_df) else 0,
    "total_active_pchembl_ml_ready_records": int(ch_df["n_pchembl_activities"].sum()) if len(ch_df) else 0,
    "by_reason_class_docs": ch_df["reason_class"].value_counts().to_dict() if len(ch_df) else {},
    "by_reason_class_activities": ch_df.groupby("reason_class")["n_activities"].sum().to_dict() if len(ch_df) else {},
}
print("PILLAR 4 (ChEMBL):", json.dumps(chembl_summary, indent=2), flush=True)

master_pilot = {
    "audit_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    "master_retracted_pmids": len(all_ret_pmids),
    "master_retracted_dois": len(rw_dois),
    "pillar1_clinvar": clinvar_summary,
    "pillar2_rcsb_pdb": pdb_summary,
    "pillar3_uniprot_swissprot": uniprot_summary,
    "pillar4_chembl": chembl_summary,
}
with open(os.path.join(OUT, "c2_zombie_4db_pilot_summary.json"), "w") as f:
    json.dump(master_pilot, f, indent=2)
print("\n=== ALL 4 PILLARS COMPLETE! Saved c2_zombie_4db_pilot_summary.json ===", flush=True)
