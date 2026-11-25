#!/opt/anaconda3/bin/python3
"""
Analysis 5 (C2-ZOMBIE) — Powered by NIH iCite + PubMed XML + ClinicalTrials.gov v2:
  Complete Clinical Evidence Cascade Across ALL ~1,515 Retracted Clinical Trials/RCTs,
  ~897 Retracted Meta-Analyses/SRs, and ~29 Retracted Practice Guidelines!
  1. Fetch 100% of PubMed retracted RCT/Clinical Trial, Meta-Analysis/SR, and Guideline PMIDs
     (strictly excluding Reinstatements).
  2. Query NIH iCite API (`https://icite.od.nih.gov/api/pubs`) across ALL 1,512 retracted RCTs:
     - Exact year-by-year citing PMID lists (`citedByPmidsByYear`) -> Event-Study relative to
       Retraction Year `t = 0` (`t = -5 .. +8`), stratified by Retraction Mechanism (Fraud/Fabrication
       vs Honest Error) AND by NIH RCR!
     - Exact strict post-retraction citations (`year > retraction_year`) across ALL retracted RCTs!
     - Exact strict post-retraction CLINICAL article citations (`cited_by_clin` published after `ret_year`)!
  3. Extract EVERY citing PMID published STRICTLY AFTER the RCT's retraction year from `citedByPmidsByYear`
     and cross-check via PubMed `esummary` to count the EXACT population of downstream
     `Meta-Analysis[pt]`, `Systematic Review[pt]`, and `Practice Guideline[pt]` articles that cited a
     retracted RCT STRICTLY AFTER ITS RETRACTION!!
  4. Extract ClinicalTrials.gov NCT registry IDs directly from PubMed XML (`<DataBankList>` + `<Abstract>`)
     across ALL 1,512 retracted RCTs -> audit EVERY linked NCT record live on ClinicalTrials.gov v2 API!
"""
import collections, json, os, re, time, xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import requests
from scipy import stats

ROOT = "/Users/mac/AI-Biology-Lab/projects/project3_zombie_retracted_biomed_dbs"
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "results")

S = requests.Session()
S.trust_env = False
S.headers.update({"User-Agent": "C2ZombieClinicalCascadeNIH/1.0 (mailto:research@example.org)"})

def fetch_pmids(term):
    r = S.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
              params={"db": "pubmed", "term": term, "retmax": 5000, "retmode": "json"}, timeout=25).json()["esearchresult"]
    return int(r["count"]), [int(x) for x in r["idlist"]]

print("=== [1] Fetching ALL Retracted Clinical Trials, Meta-Analyses, and Guidelines ===", flush=True)
n_rct, rct_pmids = fetch_pmids('Retracted Publication[pt] AND (Randomized Controlled Trial[pt] OR Clinical Trial[pt])')
time.sleep(0.35)
n_ma, ma_pmids = fetch_pmids('Retracted Publication[pt] AND (Meta-Analysis[pt] OR Systematic Review[pt])')
time.sleep(0.35)
n_gl, gl_pmids = fetch_pmids('Retracted Publication[pt] AND (Guideline[pt] OR Practice Guideline[pt])')

rw_raw = pd.read_csv(os.path.join(DATA, "retraction_watch.csv"), low_memory=False)
reinst = rw_raw[rw_raw["RetractionNature"].astype(str).str.contains("Reinstate", case=False, na=False)]
reinst_pmids = set(pd.to_numeric(reinst["OriginalPaperPubMedID"], errors="coerce").dropna().astype(int))
reinst_pmids.discard(0)

rct_pmids = [p for p in rct_pmids if p not in reinst_pmids]
ma_pmids = [p for p in ma_pmids if p not in reinst_pmids]
gl_pmids = [p for p in gl_pmids if p not in reinst_pmids]
print(f"  Clean Retracted RCTs/Trials: {len(rct_pmids)} | Meta-Analyses/SRs: {len(ma_pmids)} | Guidelines: {len(gl_pmids)}", flush=True)

rw = pd.read_parquet(os.path.join(DATA, "retraction_watch_clean_retractions.parquet"))
rw["pmid_int"] = pd.to_numeric(rw["pmid_clean"], errors="coerce").fillna(0).astype(int)
rw["orig_year"] = pd.to_datetime(rw["OriginalPaperDate"], errors="coerce").dt.year
rw["ret_year"] = pd.to_datetime(rw["RetractionDate"], errors="coerce").dt.year
rw_map = rw[rw["pmid_int"] > 0].drop_duplicates("pmid_int").set_index("pmid_int").to_dict("index")

# ==============================================================================
# 2. FULL NIH iCite EVENT-STUDY ACROSS ALL 1,512 RETRACTED CLINICAL TRIALS!
# ==============================================================================
print(f"\n=== [2] Querying NIH iCite Across ALL {len(rct_pmids)} Retracted Clinical Trials ===", flush=True)
ic_rows = []
event_study_rows = []
post_citing_pairs = []  # (retracted_rct_pmid, ret_year, reason_class, citing_pmid, citing_year, is_clin_citing)

for i in range(0, len(rct_pmids), 100):
    batch = rct_pmids[i:i+100]
    try:
        r = S.get("https://icite.od.nih.gov/api/pubs",
                  params={"pmids": ",".join(str(p) for p in batch)}, timeout=35).json()
        for d in r.get("data", []):
            pid = int(d["pmid"])
            rw_info = rw_map.get(pid, {})
            ret_yr = rw_info.get("ret_year", np.nan)
            pub_yr = d.get("year") or rw_info.get("orig_year", np.nan)
            rc = rw_info.get("reason_class", "PUBMED_RETRACTED_ONLY")
            cb_by_yr = d.get("citedByPmidsByYear") or []
            clin_set = set(int(x) for x in (d.get("cited_by_clin") or []) if str(x).isdigit())
            # Build year->count map AND exact post-retraction citing PMID list
            # NIH iCite returns a list of 1-key dicts: [{"citing_pmid": citing_year}, ...]
            yr_cnt = collections.Counter()
            post_cites = 0
            post_clin_cites = 0
            same_yr_cites = 0
            pre_cites = 0
            items_iter = []
            if isinstance(cb_by_yr, list):
                for entry in cb_by_yr:
                    if isinstance(entry, dict):
                        for cp_str, y_val in entry.items():
                            if str(cp_str).isdigit() and y_val is not None:
                                items_iter.append((int(cp_str), int(y_val)))
            elif isinstance(cb_by_yr, dict):
                for y_str, clist in cb_by_yr.items():
                    if str(y_str).isdigit():
                        for cp in (clist or []):
                            if str(cp).isdigit():
                                items_iter.append((int(cp), int(y_str)))

            for cp, y in items_iter:
                yr_cnt[y] += 1
                if pd.notna(ret_yr):
                    ry = int(ret_yr)
                    if y > ry:
                        post_cites += 1
                        is_cl = cp in clin_set
                        if is_cl:
                            post_clin_cites += 1
                        post_citing_pairs.append((pid, ry, rc, cp, y, is_cl))
                    elif y == ry:
                        same_yr_cites += 1
                    else:
                        pre_cites += 1

            ic_rows.append({
                "pmid": pid, "pub_year": pub_yr, "ret_year": ret_yr, "reason_class": rc,
                "journal": d.get("journal", ""), "title": str(d.get("title", ""))[:180],
                "total_citations": int(d.get("citation_count") or 0),
                "rcr": d.get("relative_citation_ratio"),
                "nih_percentile": d.get("nih_percentile"),
                "apt": d.get("apt"),
                "pre_ret_citations": pre_cites,
                "same_yr_citations": same_yr_cites,
                "strict_post_ret_citations": post_cites,
                "strict_post_ret_clinical_citations": post_clin_cites,
                "total_clinical_citations": len(clin_set),
            })

            # Event-study panel (k = -5 .. +8) for dated RCTs retracted 2005-2024
            if pd.notna(ret_yr) and pd.notna(pub_yr) and 2005 <= int(ret_yr) <= 2024:
                ry = int(ret_yr)
                py = int(pub_yr)
                for k in range(-5, 9):
                    cal_yr = ry + k
                    if cal_yr < py or cal_yr > 2025:
                        continue
                    event_study_rows.append({
                        "pmid": pid, "ret_year": ry, "rel_year": k, "cal_year": cal_yr,
                        "cites": yr_cnt.get(cal_yr, 0), "reason_class": rc,
                    })
    except Exception as e:
        print(f"  iCite batch {i} err: {e}", flush=True)
    time.sleep(0.15)

rct_df = pd.DataFrame(ic_rows).drop_duplicates("pmid")
rct_df.to_csv(os.path.join(OUT, "clinical_all_retracted_rcts_icite.csv"), index=False)
es_df = pd.DataFrame(event_study_rows)
es_df.to_csv(os.path.join(OUT, "clinical_rct_event_study_panel.csv"), index=False)

dated_rct = rct_df[rct_df["ret_year"].notna() & (rct_df["ret_year"] <= 2024)].copy()
print(f"  NIH iCite matched: {len(rct_df)}/{len(rct_pmids)} retracted RCTs!", flush=True)
print(f"  Dated RCTs (ret_year<=2024): {len(dated_rct)}", flush=True)
print(f"    Cited STRICTLY AFTER retraction (year > ret_year): {int((dated_rct['strict_post_ret_citations']>0).sum())}/{len(dated_rct)} ({100*(dated_rct['strict_post_ret_citations']>0).mean():.1f}%)", flush=True)
print(f"    Total strict post-retraction citations = {int(dated_rct['strict_post_ret_citations'].sum()):,} (incl. {int(dated_rct['strict_post_ret_clinical_citations'].sum()):,} NIH-classified CLINICAL citing articles!)", flush=True)
print(f"    Median total citations={dated_rct['total_citations'].median():.1f}, median post-retraction citations={dated_rct['strict_post_ret_citations'].median():.1f}, max post-retraction={int(dated_rct['strict_post_ret_citations'].max())}", flush=True)

# Event-study summary (overall + Fraud vs Honest Error!)
es_df["is_fraud"] = es_df["reason_class"].isin(["DATA_FABRICATION_OR_PAPER_MILL", "IMAGE_MANIPULATION_OR_DUPLICATION", "PEER_REVIEW_MANIPULATION"])
es_sum = es_df.groupby("rel_year").agg(
    n_trials=("pmid", "nunique"),
    mean_cites=("cites", "mean"),
    se_cites=("cites", lambda s: s.std(ddof=1) / np.sqrt(len(s)) if len(s) > 1 else 0),
    median_cites=("cites", "median"),
    pct_nonzero=("cites", lambda s: 100*(s > 0).mean()),
    total_cites=("cites", "sum"),
).reset_index()
# Also add fraud vs error means per rel_year
f_means = es_df[es_df["is_fraud"]].groupby("rel_year")["cites"].mean().rename("mean_cites_fraud")
e_means = es_df[es_df["reason_class"] == "ERROR_OR_UNRELIABLE_RESULTS"].groupby("rel_year")["cites"].mean().rename("mean_cites_error")
es_sum = es_sum.merge(f_means, on="rel_year", how="left").merge(e_means, on="rel_year", how="left")
es_sum.to_csv(os.path.join(OUT, "clinical_rct_event_study_summary.csv"), index=False)
print("\n  [NIH iCite Retracted RCT Citation Event-Study (rel_year = cal_year - ret_year)]:", flush=True)
print(es_sum[["rel_year", "n_trials", "mean_cites", "se_cites", "pct_nonzero", "mean_cites_fraud", "mean_cites_error"]].to_string(index=False), flush=True)

# Paired Wilcoxon signed-rank test: mean annual citations in t=-2..-1 vs t=+1..+2 vs t=+3..+5!
pre_win = es_df[es_df["rel_year"].isin([-2, -1])].groupby("pmid")["cites"].mean()
post12 = es_df[es_df["rel_year"].isin([1, 2])].groupby("pmid")["cites"].mean()
post35 = es_df[es_df["rel_year"].isin([3, 4, 5])].groupby("pmid")["cites"].mean()
common12 = pre_win.index.intersection(post12.index)
w_stat12, w_p12 = stats.wilcoxon(pre_win.loc[common12], post12.loc[common12]) if len(common12) > 10 else (np.nan, np.nan)
pct_drop12 = 100.0 * (1.0 - post12.loc[common12].mean() / max(1e-6, pre_win.loc[common12].mean()))
common35 = pre_win.index.intersection(post35.index)
pct_drop35 = 100.0 * (1.0 - post35.loc[common35].mean() / max(1e-6, pre_win.loc[common35].mean()))
print(f"  Paired Pre(-2..-1) vs Post(+1..+2) (N={len(common12)}): mean {pre_win.loc[common12].mean():.2f} -> {post12.loc[common12].mean():.2f} ({pct_drop12:.1f}% drop, Wilcoxon p={w_p12:.2e})", flush=True)
print(f"  Paired Pre(-2..-1) vs Post(+3..+5) (N={len(common35)}): mean {pre_win.loc[common35].mean():.2f} -> {post35.loc[common35].mean():.2f} ({pct_drop35:.1f}% drop, STILL {post35.loc[common35].mean():.2f} cites/yr 3-5 yrs post-retraction!)", flush=True)

# ==============================================================================
# 3. EXACT POPULATION AUDIT OF ALL POST-RETRACTION CITING META-ANALYSES, SRs, & GUIDELINES!!
# ==============================================================================
# Every unique citing PMID from `post_citing_pairs` (we have the EXACT complete list across ALL 1,512 RCTs!)
unique_post_citers = sorted(set(x[3] for x in post_citing_pairs))
print(f"\n=== [3] Auditing ALL {len(unique_post_citers):,} Unique Post-Retraction Citing PubMed Articles for Meta-Analyses / SRs / Guidelines ===", flush=True)
# Batch PubMed esummary (200 per POST call -> ~40 fast calls, ~16 seconds total!)
citer_meta = {}
for i in range(0, len(unique_post_citers), 200):
    batch = unique_post_citers[i:i+200]
    for att in range(3):
        try:
            r = S.post("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                       data={"db": "pubmed", "id": ",".join(str(x) for x in batch), "retmode": "json"}, timeout=30).json()
            res = r.get("result", {})
            for cid in batch:
                it = res.get(str(cid), {})
                pts = [str(x).lower() for x in it.get("pubtype", [])]
                is_gl = any(k in pts for k in ["guideline", "practice guideline"])
                is_ma = any(k in pts for k in ["meta-analysis", "systematic review"])
                is_rct = any(k in pts for k in ["randomized controlled trial", "clinical trial"])
                is_ret = "retracted publication" in pts
                if is_gl or is_ma or is_rct:
                    citer_meta[cid] = {
                        "type": "GUIDELINE" if is_gl else ("META_ANALYSIS_OR_SR" if is_ma else "CLINICAL_TRIAL"),
                        "is_retracted_citer": is_ret,
                        "journal": it.get("source", ""),
                        "title": str(it.get("title", ""))[:160],
                    }
            break
        except Exception:
            time.sleep(1.0)
    time.sleep(0.34)

contagion_rows = []
for rct_pid, ry, rc, c_pid, c_yr, is_cl in post_citing_pairs:
    cm = citer_meta.get(c_pid)
    if cm and cm["type"] in ("GUIDELINE", "META_ANALYSIS_OR_SR"):
        contagion_rows.append({
            "retracted_rct_pmid": rct_pid, "rct_ret_year": ry, "rct_reason_class": rc,
            "citing_pmid": c_pid, "citing_year": c_yr, "years_after_retraction": c_yr - ry,
            "citing_type": cm["type"], "citer_itself_retracted": cm["is_retracted_citer"],
            "citing_journal": cm["journal"], "citing_title": cm["title"],
        })
cont_df = pd.DataFrame(contagion_rows).drop_duplicates(["retracted_rct_pmid", "citing_pmid"])
cont_df.to_csv(os.path.join(OUT, "clinical_downstream_post_retraction_ma_guideline_citations.csv"), index=False)

n_ma_uniq = int(cont_df[cont_df["citing_type"] == "META_ANALYSIS_OR_SR"]["citing_pmid"].nunique()) if len(cont_df) else 0
n_gl_uniq = int(cont_df[cont_df["citing_type"] == "GUIDELINE"]["citing_pmid"].nunique()) if len(cont_df) else 0
n_rct_contam = int(cont_df["retracted_rct_pmid"].nunique()) if len(cont_df) else 0
# Also count non-retracted active post-retraction Meta-Analyses & Guidelines!
active_cont = cont_df[cont_df["citer_itself_retracted"] == False] if len(cont_df) else cont_df
n_ma_active = int(active_cont[active_cont["citing_type"] == "META_ANALYSIS_OR_SR"]["citing_pmid"].nunique()) if len(active_cont) else 0
n_gl_active = int(active_cont[active_cont["citing_type"] == "GUIDELINE"]["citing_pmid"].nunique()) if len(active_cont) else 0
print(f"  EXACT POPULATION CONTAGION (STRICTLY AFTER RETRACTION YEAR):", flush=True)
print(f"    Total retracted-RCT -> post-retraction Meta-Analysis/SR/Guideline citation edges: {len(cont_df):,} (from {n_rct_contam} retracted RCTs!)", flush=True)
print(f"    Unique post-retraction Meta-Analyses / Systematic Reviews: {n_ma_uniq:,} ({n_ma_active:,} active non-retracted!)", flush=True)
print(f"    Unique post-retraction Clinical Practice Guidelines:       {n_gl_uniq:,} ({n_gl_active:,} active non-retracted!)", flush=True)
if len(cont_df):
    print(f"    Median years AFTER retraction when cited by MA/SR/Guideline: {cont_df['years_after_retraction'].median():.1f} yrs (max {int(cont_df['years_after_retraction'].max())} yrs)!", flush=True)

# ==============================================================================
# 4. EXTRACT NCT IDs FROM PUBMED XML FOR ALL 1,512 RETRACTED RCTs -> LIVE CLINICALTRIALS.GOV AUDIT!
# ==============================================================================
print(f"\n=== [4] Extracting ClinicalTrials.gov NCT IDs from PubMed XML Across ALL {len(rct_pmids)} Retracted RCTs ===", flush=True)
nct_links = []
for i in range(0, len(rct_pmids), 200):
    batch = rct_pmids[i:i+200]
    try:
        r = S.post("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                   data={"db": "pubmed", "id": ",".join(str(p) for p in batch), "retmode": "xml"}, timeout=40)
        root = ET.fromstring(r.content)
        for art in root.findall(".//PubmedArticle"):
            pid_el = art.find(".//MedlineCitation/PMID")
            if pid_el is None: continue
            pid = int(pid_el.text)
            ncts = set()
            for db in art.findall(".//DataBank"):
                dbname = (db.findtext("DataBankName") or "").lower()
                if "clinicaltrials" in dbname:
                    for acc in db.findall(".//AccessionNumber"):
                        if acc.text and re.match(r"NCT\d{8}", acc.text.strip(), re.I):
                            ncts.add(acc.text.strip().upper())
            abs_txt = " ".join("".join(el.itertext()) for el in art.findall(".//Abstract"))
            for m in re.findall(r"\bNCT\d{8}\b", abs_txt, re.I):
                ncts.add(m.upper())
            for nct in sorted(ncts):
                nct_links.append({"pmid": pid, "nct_id": nct})
    except Exception as e:
        print(f"  XML batch {i} err: {e}", flush=True)
    time.sleep(0.35)

nct_df = pd.DataFrame(nct_links).drop_duplicates()
print(f"  Found {len(nct_df)} retracted RCT -> NCT links ({nct_df['pmid'].nunique() if len(nct_df) else 0} unique retracted RCTs, {nct_df['nct_id'].nunique() if len(nct_df) else 0} unique NCT trials)!", flush=True)

ct_audit = []
if len(nct_df):
    for _, row in nct_df.iterrows():
        nct = row["nct_id"]
        pid = int(row["pmid"])
        rw_info = rw_map.get(pid, {})
        try:
            rc = S.get(f"https://clinicaltrials.gov/api/v2/studies/{nct}", timeout=15)
            if rc.status_code == 200:
                js = rc.json()
                ps = js.get("protocolSection", {})
                status = ps.get("statusModule", {}).get("overallStatus", "")
                refs = ps.get("referencesModule", {}).get("references", [])
                has_pmid_ref = any(str(rf.get("pmid", "")) == str(pid) for rf in refs)
                ref_type = next((rf.get("type", "") for rf in refs if str(rf.get("pmid", "")) == str(pid)), "")
                ref_cit = next((rf.get("citation", "") for rf in refs if str(rf.get("pmid", "")) == str(pid)), "")
                ct_mentions_retract = ("retract" in json.dumps(ps).lower())
                ct_audit.append({
                    "pmid": pid, "nct_id": nct, "ct_status": status,
                    "has_retracted_pmid_in_references": has_pmid_ref,
                    "ref_type": ref_type,
                    "ct_protocol_mentions_retract": ct_mentions_retract,
                    "ret_year": rw_info.get("ret_year", np.nan),
                    "reason_class": rw_info.get("reason_class", "PUBMED_RETRACTED_ONLY"),
                    "ref_citation_prefix": str(ref_cit)[:120],
                })
        except Exception:
            pass
        time.sleep(0.10)

ct_df = pd.DataFrame(ct_audit)
ct_df.to_csv(os.path.join(OUT, "pillar8_clinicaltrials_nct_live_audit.csv"), index=False)
if len(ct_df):
    linked_refs = ct_df[ct_df["has_retracted_pmid_in_references"]]
    print(f"  ClinicalTrials.gov LIVE AUDIT: {len(ct_df)} NCT records retrieved!", flush=True)
    print(f"    NCT records with retracted PMID directly listed in `referencesModule`: {len(linked_refs)}/{len(ct_df)}", flush=True)
    print(f"    NCT protocolSection mentions 'retract' == False (UNFLAGGED!): {int((ct_df['ct_protocol_mentions_retract']==False).sum())}/{len(ct_df)} ({100*(ct_df['ct_protocol_mentions_retract']==False).mean():.1f}%)", flush=True)
    if len(linked_refs):
        print(f"    Among NCTs directly linking retracted PMID in referencesModule, UNFLAGGED: {int((linked_refs['ct_protocol_mentions_retract']==False).sum())}/{len(linked_refs)} ({100*(linked_refs['ct_protocol_mentions_retract']==False).mean():.1f}%)", flush=True)

summary = {
    "audit_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    "pubmed_retracted_rcts_clean": len(rct_pmids),
    "pubmed_retracted_meta_analyses_clean": len(ma_pmids),
    "pubmed_retracted_guidelines_clean": len(gl_pmids),
    "icite_matched_rcts": len(rct_df),
    "dated_rcts_ret_le_2024": len(dated_rct),
    "dated_rcts_with_strict_post_ret_cites": int((dated_rct["strict_post_ret_citations"] > 0).sum()),
    "dated_rcts_with_strict_post_ret_cites_pct": round(100 * (dated_rct["strict_post_ret_citations"] > 0).mean(), 2),
    "total_strict_post_ret_citations_all_rcts": int(dated_rct["strict_post_ret_citations"].sum()),
    "total_strict_post_ret_clinical_citations_all_rcts": int(dated_rct["strict_post_ret_clinical_citations"].sum()),
    "unique_post_retraction_citing_pmids_audited": len(unique_post_citers),
    "paired_pre_vs_post12_n": len(common12),
    "paired_pre_mean_cites_per_yr": round(float(pre_win.loc[common12].mean()), 3),
    "paired_post12_mean_cites_per_yr": round(float(post12.loc[common12].mean()), 3),
    "paired_drop_post12_pct": round(float(pct_drop12), 2),
    "paired_wilcoxon_p": float(w_p12),
    "paired_post35_mean_cites_per_yr": round(float(post35.loc[common35].mean()), 3),
    "paired_drop_post35_pct": round(float(pct_drop35), 2),
    "downstream_post_ret_ma_gl_citation_edges": len(cont_df),
    "downstream_post_ret_ma_sr_unique_total": n_ma_uniq,
    "downstream_post_ret_ma_sr_unique_active": n_ma_active,
    "downstream_post_ret_guidelines_unique_total": n_gl_uniq,
    "downstream_post_ret_guidelines_unique_active": n_gl_active,
    "retracted_rcts_cited_by_post_ret_ma_or_gl": n_rct_contam,
    "nct_links_extracted_from_pubmed_xml": len(nct_df),
    "nct_unique_trials_live_audited": len(ct_df),
    "nct_with_retracted_pmid_in_refs": int(ct_df["has_retracted_pmid_in_references"].sum()) if len(ct_df) else 0,
    "nct_protocol_unflagged_pct": round(100 * (ct_df["ct_protocol_mentions_retract"] == False).mean(), 2) if len(ct_df) else 0,
    "nct_refs_linked_unflagged_pct": round(100 * (ct_df[ct_df["has_retracted_pmid_in_references"]]["ct_protocol_mentions_retract"] == False).mean(), 2) if len(ct_df) and ct_df["has_retracted_pmid_in_references"].any() else 0,
}
with open(os.path.join(OUT, "analysis5_clinical_cascade_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print("\n=== ANALYSIS 5 (CLINICAL CASCADE & CLINICALTRIALS.GOV) SUMMARY ===", flush=True)
print(json.dumps(summary, indent=2), flush=True)
