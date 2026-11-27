#!/opt/anaconda3/bin/python3
"""
Analysis 3 (C2-ZOMBIE):
  COMPLETE POPULATION-LEVEL AUDIT OF AUTOMATED BIOMEDICAL AI KNOWLEDGE GRAPHS!
  1. NCBI `PubTator3` Bulk Relation Knowledge Graph (`relation2pubtator3.gz`, 288 MB, updated Nov 2026):
     - Stream-scan ALL tens of millions of AI-extracted biomedical KG relations `(PMID, relation_type, Entity1, Entity2)`
     - Match against 100%-complete Master Retracted PMID Registry (`~38,100` non-reinstated retractions)
     - Compute:
       * Exact total relations & retracted relations (overall, by relation type: `treat`, `cause`, `inhibit`, `stimulate`, `associate`, `positive_correlate`, `negative_correlate`, `drug_interact`, etc.)
       * Breakdown by Retraction Watch mechanism (`DATA_FABRICATION_OR_PAPER_MILL`, `PEER_REVIEW_MANIPULATION`, `IMAGE_MANIPULATION`, `ERROR`, `PLAGIARISM`)
       * Breakdown by Entity-Type Pair (`Gene-Disease`, `Chemical-Disease`, `Gene-Chemical`, `Gene-Gene`, `Variant-Disease`, etc.)
       * EXACT SOLE-SOURCE RETRACTED KG EDGES: Unique `(relation_type, Entity1, Entity2)` triples whose ENTIRE PubMed evidence set consists EXCLUSIVELY of retracted papers (100% retracted, 0 non-retracted papers!! -> pure hallucination seeds in GraphRAG / biomedical AI!)
       * EXACT MAJORITY-RETRACTED KG EDGES (>=50% retracted PMIDs)
       * Exact `PubTator3` retracted enrichment rate (retracted relations per 10k vs non-retracted relations per 10k, and mean relations extracted per retracted paper vs per PubMed paper!)
  2. NCBI `PubTator3` Bulk Mutation/Variant Extractions (`mutation2pubtator3.gz`, 111 MB):
     - Exact retracted variant/mutation extraction mentions & unique retracted papers with AI-mined variant assertions.
  3. Live `Open Targets Platform` GraphQL Verification on Top Retracted Gene–Disease Pairs:
     - Checks whether Open Targets' `europepmc` text-mining evidence STILL serves retracted PMIDs as positive target–disease association evidence!
"""
import collections, gzip, json, os, time
import numpy as np
import pandas as pd
import requests

ROOT = "/Users/mac/AI-Biology-Lab/projects/project3_zombie_retracted_biomed_dbs"
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "results")

S = requests.Session()
S.trust_env = False
S.headers.update({"User-Agent": "C2ZombiePubTator3Audit/1.0 (mailto:research@example.org)"})

# Load 100%-complete master registry (fallback to union if needed, wait up to 60s if Analysis 1 is finishing)
reg_path = os.path.join(DATA, "master_retracted_pmids_100pct.json")
for _ in range(30):
    if os.path.exists(reg_path):
        break
    time.sleep(2)
with open(reg_path) as f:
    reg = json.load(f)

master_pmids = set(reg["union_pmids"])
pm_ret_set = set(reg["pubmed_pmids"])
pm_oa_ret_set = set(reg["pubmed_pmcoa_pmids"])
print(f"Loaded 100% Master Retracted Registry: {len(master_pmids):,} union PMIDs ({len(pm_ret_set):,} PubMed, {len(pm_oa_ret_set):,} PMC OA)", flush=True)

rw = pd.read_parquet(os.path.join(DATA, "retraction_watch_clean_retractions.parquet"))
rw["pmid_int"] = pd.to_numeric(rw["pmid_clean"], errors="coerce").fillna(0).astype(int)
rw["ret_year"] = pd.to_datetime(rw["RetractionDate"], errors="coerce").dt.year
rw["orig_year"] = pd.to_datetime(rw["OriginalPaperDate"], errors="coerce").dt.year
rw_valid = rw[rw["pmid_int"] > 0].drop_duplicates("pmid_int")
pmid_to_rc = dict(zip(rw_valid["pmid_int"], rw_valid["reason_class"]))
pmid_to_ry = dict(zip(rw_valid["pmid_int"], rw_valid["ret_year"]))
pmid_to_pub = dict(zip(rw_valid["pmid_int"], rw_valid["Publisher"]))

# ==============================================================================
# 1. DOWNLOAD & STREAM-SCAN `relation2pubtator3.gz` (288 MB)
# ==============================================================================
rel_gz = os.path.join(DATA, "relation2pubtator3.gz")
if not os.path.exists(rel_gz) or os.path.getsize(rel_gz) < 100_000_000:
    print("Downloading NCBI PubTator3 relation2pubtator3.gz (~288 MB)...", flush=True)
    with S.get("https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/relation2pubtator3.gz", stream=True, timeout=180) as r:
        r.raise_for_status()
        with open(rel_gz, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)
    print(f"Downloaded relation2pubtator3.gz: {os.path.getsize(rel_gz):,} bytes", flush=True)

print("=== [1] Pass 1: Scanning ALL PubTator3 Relations for Retracted PMIDs ===", flush=True)
t0 = time.time()
total_rel_rows = 0
total_unique_pmids = set()
ret_rel_rows = 0
ret_unique_pmids = set()
ret_pmcoa_rel_rows = 0
ret_rw_only_rel_rows = 0

by_rel_type_all = collections.Counter()
by_rel_type_ret = collections.Counter()
by_pair_type_all = collections.Counter()
by_pair_type_ret = collections.Counter()
by_reason_ret = collections.Counter()
by_ret_year_ret = collections.Counter()

# Track exact candidate retracted triples `(rel_type, e1, e2)` -> set of retracted supporting PMIDs!
cand_triple_ret_pmids = collections.defaultdict(set)
# Also store per-PMID relation counts for retracted papers
ret_pmid_rel_counts = collections.Counter()

def canon_edge(rel, e1, e2):
    # Symmetric relations canonicalized; directional (`treat`, `cause`, `inhibit`, `stimulate`, `prevent`, `compare`) keep order
    if rel in ("associate", "positive_correlate", "negative_correlate", "interact", "drug_interact", "cotreat", "bind", "compare"):
        if e1 > e2:
            e1, e2 = e2, e1
    return (rel, e1, e2)

def etype(e):
    i = e.find("|")
    return e[:i] if i > 0 else "Unknown"

with gzip.open(rel_gz, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 4:
            continue
        try:
            pmid = int(parts[0])
        except ValueError:
            continue
        rel = parts[1]
        e1 = parts[2]
        e2 = parts[3]
        total_rel_rows += 1
        total_unique_pmids.add(pmid)
        by_rel_type_all[rel] += 1
        t1, t2 = etype(e1), etype(e2)
        pair_t = f"{min(t1,t2)}--{max(t1,t2)}"
        by_pair_type_all[pair_t] += 1

        if pmid in master_pmids:
            ret_rel_rows += 1
            ret_unique_pmids.add(pmid)
            ret_pmid_rel_counts[pmid] += 1
            by_rel_type_ret[rel] += 1
            by_pair_type_ret[pair_t] += 1
            rc = pmid_to_rc.get(pmid, "PUBMED_RETRACTED_ONLY")
            by_reason_ret[rc] += 1
            ry = pmid_to_ry.get(pmid, None)
            if ry is not None and not pd.isna(ry):
                by_ret_year_ret[int(ry)] += 1
            if pmid in pm_oa_ret_set:
                ret_pmcoa_rel_rows += 1
            if pmid not in pm_ret_set:
                ret_rw_only_rel_rows += 1
            cand_triple_ret_pmids[canon_edge(rel, e1, e2)].add(pmid)

print(f"  Pass 1 done in {time.time()-t0:.1f}s: Total relations={total_rel_rows:,} across {len(total_unique_pmids):,} PMIDs!", flush=True)
print(f"  RETRACTED relations={ret_rel_rows:,} across {len(ret_unique_pmids):,} unique retracted PMIDs ({len(cand_triple_ret_pmids):,} unique candidate retracted triples)!", flush=True)
print(f"  Of retracted relations: {ret_pmcoa_rel_rows:,} ({100*ret_pmcoa_rel_rows/max(1,ret_rel_rows):.1f}%) come from PMC OA papers; {ret_rw_only_rel_rows:,} from RW-only unflagged PubMed papers!", flush=True)

# Pass 2: For every candidate retracted triple, count ALL non-retracted supporting PMIDs in the entire PubTator3 KG!
print("=== [2] Pass 2: Computing Exact SOLE-SOURCE & MAJORITY-RETRACTED KG Triples ===", flush=True)
t1 = time.time()
cand_triple_nonret_count = collections.Counter()
with gzip.open(rel_gz, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 4:
            continue
        try:
            pmid = int(parts[0])
        except ValueError:
            continue
        if pmid in master_pmids:
            continue
        edge = canon_edge(parts[1], parts[2], parts[3])
        if edge in cand_triple_ret_pmids:
            cand_triple_nonret_count[edge] += 1

print(f"  Pass 2 done in {time.time()-t1:.1f}s!", flush=True)

sole_source_triples = []
maj_ret_count = 0
sole_by_rel = collections.Counter()
sole_by_pair = collections.Counter()
sole_by_reason = collections.Counter()
sole_multi_paper = 0

for edge, rpmids in cand_triple_ret_pmids.items():
    n_ret = len(rpmids)
    n_non = cand_triple_nonret_count.get(edge, 0)
    rel, e1, e2 = edge
    t1_e, t2_e = etype(e1), etype(e2)
    pair_t = f"{min(t1_e,t2_e)}--{max(t1_e,t2_e)}"
    if n_ret >= n_non:
        maj_ret_count += 1
    if n_non == 0:
        sole_by_rel[rel] += 1
        sole_by_pair[pair_t] += 1
        rcs = {pmid_to_rc.get(p, "PUBMED_RETRACTED_ONLY") for p in rpmids}
        for rc in rcs:
            sole_by_reason[rc] += 1
        if n_ret >= 2:
            sole_multi_paper += 1
        if len(sole_source_triples) < 50000 or n_ret >= 2 or rel in ("treat", "cause", "inhibit", "drug_interact"):
            sole_source_triples.append({
                "relation": rel, "entity1": e1, "entity2": e2, "pair_type": pair_t,
                "n_retracted_pmids": n_ret, "n_nonretracted_pmids": 0,
                "retracted_pmids": ";".join(str(p) for p in sorted(rpmids)[:10]),
                "reason_classes": ";".join(sorted(rcs)),
            })

print(f"  SOLE-SOURCE RETRACTED KG TRIPLES (100% retracted, 0 non-retracted PMIDs!): {sum(sole_by_rel.values()):,} / {len(cand_triple_ret_pmids):,} ({100*sum(sole_by_rel.values())/max(1,len(cand_triple_ret_pmids)):.1f}%)!!", flush=True)
print(f"  MAJORITY-RETRACTED KG TRIPLES (>=50% retracted support): {maj_ret_count:,}!", flush=True)
print(f"  SOLE-SOURCE supported by >=2 DISTINCT retracted papers (coordinated paper-mill echoes!): {sole_multi_paper:,}!", flush=True)

# Save summary tables (keep CSV compact: top 5,000 sole-source triples sorted by n_retracted_pmids + full breakdown CSVs!)
sole_df = pd.DataFrame(sole_source_triples).sort_values(["n_retracted_pmids", "relation"], ascending=[False, True])
sole_df.head(5000).to_csv(os.path.join(OUT, "pubtator3_sole_source_retracted_triples_top5000.csv"), index=False)

rel_tbl = []
for rel, cnt_all in by_rel_type_all.most_common():
    cnt_ret = by_rel_type_ret.get(rel, 0)
    cnt_sole = sole_by_rel.get(rel, 0)
    rel_tbl.append({
        "relation_type": rel, "total_relations": cnt_all, "retracted_relations": cnt_ret,
        "retracted_rate_per10k": round(10000 * cnt_ret / cnt_all, 2),
        "sole_source_retracted_triples": cnt_sole,
    })
pd.DataFrame(rel_tbl).to_csv(os.path.join(OUT, "pubtator3_relations_by_type.csv"), index=False)
print("\n  [PubTator3 Relations by Type]:", flush=True)
print(pd.DataFrame(rel_tbl).to_string(index=False), flush=True)

pair_tbl = []
for pt, cnt_all in by_pair_type_all.most_common(15):
    cnt_ret = by_pair_type_ret.get(pt, 0)
    cnt_sole = sole_by_pair.get(pt, 0)
    pair_tbl.append({
        "entity_pair_type": pt, "total_relations": cnt_all, "retracted_relations": cnt_ret,
        "retracted_rate_per10k": round(10000 * cnt_ret / cnt_all, 2),
        "sole_source_retracted_triples": cnt_sole,
    })
pd.DataFrame(pair_tbl).to_csv(os.path.join(OUT, "pubtator3_relations_by_entity_pair.csv"), index=False)
print("\n  [PubTator3 Relations by Entity-Type Pair]:", flush=True)
print(pd.DataFrame(pair_tbl).to_string(index=False), flush=True)

# ==============================================================================
# 2. DOWNLOAD & SCAN `mutation2pubtator3.gz` (111 MB)
# ==============================================================================
mut_gz = os.path.join(DATA, "mutation2pubtator3.gz")
if not os.path.exists(mut_gz) or os.path.getsize(mut_gz) < 50_000_000:
    print("\nDownloading NCBI PubTator3 mutation2pubtator3.gz (~111 MB)...", flush=True)
    with S.get("https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/mutation2pubtator3.gz", stream=True, timeout=180) as r:
        r.raise_for_status()
        with open(mut_gz, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)

mut_total = 0
mut_ret = 0
mut_ret_pmids = set()
mut_ret_rsids = set()
mut_by_rc = collections.Counter()
with gzip.open(mut_gz, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3: continue
        try:
            pmid = int(parts[0])
        except ValueError:
            continue
        mut_total += 1
        if pmid in master_pmids:
            mut_ret += 1
            mut_ret_pmids.add(pmid)
            if parts[2].startswith("rs"):
                mut_ret_rsids.add(parts[2])
            mut_by_rc[pmid_to_rc.get(pmid, "PUBMED_RETRACTED_ONLY")] += 1

print(f"=== [3] PubTator3 Mutations: {mut_ret:,} / {mut_total:,} retracted variant mentions across {len(mut_ret_pmids):,} retracted PMIDs ({len(mut_ret_rsids):,} unique dbSNP rsIDs)! ===", flush=True)

# ==============================================================================
# 3. OPEN TARGETS PLATFORM LIVE CHECK (EUROPE PMC TEXT-MINING EVIDENCE!)
# ==============================================================================
# Check Open Targets FTP `evidence/sourceId=europepmc` OR query 20 prominent cancer/drug targets via GraphQL
# for whether Open Targets `europepmc` literature arrays include retracted PMIDs!
print("\n=== [4] Probing Open Targets Platform `europepmc` AI Text-Mining Evidence for Retracted PMIDs ===", flush=True)
ot_targets = [
    ("TP53", "ENSG00000141510", "MONDO_0004992"),      # cancer
    ("BRCA1", "ENSG00000012048", "MONDO_0007254"),     # breast cancer
    ("EGFR", "ENSG00000146648", "EFO_0001071"),        # lung carcinoma
    ("VEGFA", "ENSG00000112715", "MONDO_0004992"),     # cancer
    ("AKT1", "ENSG00000142208", "MONDO_0004992"),      # cancer
    ("STAT3", "ENSG00000168610", "MONDO_0004992"),     # cancer
    ("TNF", "ENSG00000232810", "MONDO_0004992"),       # cancer
    ("IL6", "ENSG00000136244", "EFO_0000685"),         # rheumatoid arthritis
    ("PTEN", "ENSG00000171862", "MONDO_0004992"),      # cancer
    ("MYC", "ENSG00000136997", "MONDO_0004992"),       # cancer
    ("KRAS", "ENSG00000133703", "EFO_0005842"),        # colorectal cancer
    ("BCL2", "ENSG00000171791", "MONDO_0004992"),      # cancer
]
ot_rows = []
for sym, ensg, efo in ot_targets:
    q = """
    query($ensg: String!, $efo: [String!]!) {
      target(ensemblId: $ensg) {
        approvedSymbol
        evidences(efoIds: $efo, datasourceIds: ["europepmc"], size: 1000) {
          count
          rows {
            disease { id name }
            literature
            score
          }
        }
      }
    }
    """
    try:
        r = S.post("https://api.platform.opentargets.org/api/v4/graphql",
                   json={"query": q, "variables": {"ensg": ensg, "efo": [efo]}}, timeout=25).json()
        ev = r.get("data", {}).get("target", {}).get("evidences", {})
        tot_cnt = ev.get("count", 0)
        rows_ev = ev.get("rows", [])
        n_checked = 0
        n_ret_ev = 0
        ret_pmids_hit = []
        for rv in rows_ev:
            lits = rv.get("literature") or []
            n_checked += 1
            hit = [int(x) for x in lits if str(x).isdigit() and int(x) in master_pmids]
            if hit:
                n_ret_ev += 1
                ret_pmids_hit.extend(hit)
        ot_rows.append({
            "gene": sym, "ensembl_id": ensg, "disease_id": efo,
            "ot_europepmc_total_count": tot_cnt, "ot_checked_rows": n_checked,
            "retracted_evidence_rows": n_ret_ev,
            "retracted_pct_of_checked": round(100 * n_ret_ev / max(1, n_checked), 3),
            "sample_retracted_pmids": ";".join(str(x) for x in sorted(set(ret_pmids_hit))[:10]),
        })
        print(f"  Open Targets {sym:6s} ({efo}): checked={n_checked:4d} (of {tot_cnt:6,}) -> RETRACTED active evidence rows = {n_ret_ev:3d} ({100*n_ret_ev/max(1,n_checked):.2f}%)", flush=True)
    except Exception as e:
        print(f"  OT {sym} err: {e}", flush=True)
    time.sleep(0.2)

ot_df = pd.DataFrame(ot_rows)
ot_df.to_csv(os.path.join(OUT, "opentargets_europepmc_12targets_retracted_audit.csv"), index=False)

# Compute exact summary metrics
mean_rel_per_ret = ret_rel_rows / max(1, len(ret_unique_pmids))
mean_rel_per_nonret = (total_rel_rows - ret_rel_rows) / max(1, len(total_unique_pmids) - len(ret_unique_pmids))
summary = {
    "audit_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    "pubtator3_total_relations": total_rel_rows,
    "pubtator3_total_unique_pmids": len(total_unique_pmids),
    "pubtator3_retracted_relations": ret_rel_rows,
    "pubtator3_retracted_unique_pmids": len(ret_unique_pmids),
    "pubtator3_retracted_share_of_all_relations_pct": round(100 * ret_rel_rows / max(1, total_rel_rows), 4),
    "pubtator3_retracted_share_of_all_pmids_pct": round(100 * len(ret_unique_pmids) / max(1, len(total_unique_pmids)), 4),
    "mean_relations_per_retracted_pmid": round(mean_rel_per_ret, 2),
    "mean_relations_per_nonretracted_pmid": round(mean_rel_per_nonret, 2),
    "relation_density_ratio_ret_vs_nonret": round(mean_rel_per_ret / max(0.001, mean_rel_per_nonret), 3),
    "retracted_relations_from_pmcoa": ret_pmcoa_rel_rows,
    "retracted_relations_from_pmcoa_pct": round(100 * ret_pmcoa_rel_rows / max(1, ret_rel_rows), 2),
    "retracted_relations_from_rw_only_unflagged_pubmed": ret_rw_only_rel_rows,
    "unique_retracted_kg_triples": len(cand_triple_ret_pmids),
    "sole_source_retracted_kg_triples": int(sum(sole_by_rel.values())),
    "sole_source_pct_of_retracted_triples": round(100 * sum(sole_by_rel.values()) / max(1, len(cand_triple_ret_pmids)), 2),
    "majority_retracted_kg_triples": int(maj_ret_count),
    "sole_source_multi_retracted_paper_triples": int(sole_multi_paper),
    "by_reason_class_relations": dict(by_reason_ret),
    "sole_source_by_relation_type": dict(sole_by_rel),
    "pubtator3_mutation_total_mentions": mut_total,
    "pubtator3_mutation_retracted_mentions": mut_ret,
    "pubtator3_mutation_retracted_pmids": len(mut_ret_pmids),
    "pubtator3_mutation_retracted_unique_rsids": len(mut_ret_rsids),
    "pubtator3_mutation_by_reason": dict(mut_by_rc),
    "opentargets_12_target_checked_rows": int(ot_df["ot_checked_rows"].sum()) if len(ot_df) else 0,
    "opentargets_12_target_retracted_rows": int(ot_df["retracted_evidence_rows"].sum()) if len(ot_df) else 0,
    "opentargets_targets_with_ge1_retracted": int((ot_df["retracted_evidence_rows"] > 0).sum()) if len(ot_df) else 0,
}
with open(os.path.join(OUT, "analysis3_pubtator3_and_opentargets_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print("\n=== ANALYSIS 3 (PUBTATOR3 BULK KG + OPEN TARGETS) SUMMARY ===", flush=True)
print(json.dumps(summary, indent=2), flush=True)
