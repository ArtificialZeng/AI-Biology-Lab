#!/opt/anaconda3/bin/python3
"""
T00 Institutional & API v2 Forensic Audit for Study B-2:
Empirically audits ClinicalTrials.gov API v2 fields under the NLM 2020-01-29 QC Posting Policy
(https://www.nlm.nih.gov/pubs/techbull/jf20/jf20_clinicaltrials_qc.html):
  1. Lists all API v2 /studies/metadata fields related to results submission, QC, posting, and unposted tracking.
  2. Queries live counts and full record structures across 2026-04..2026-11, including hasResults=false checks.
  3. Decomposes resultsFirstSubmitDate vs resultsFirstSubmitQcDate vs resultsFirstPostDate in our local universe.
  4. Audits derivedSection.miscInfoModule.submissionTracking (and protocolSection.statusModule) for:
       (a) Overdue ACTs with hasResults=false updated since 2026-04-01 (do they hold unposted/in-QC submissions?!)
       (b) Posted 2026 ACTs with lag<=30d vs lag>60d (identifying posted-with-QC-comments vs post-QC!).
  5. Audits ALL overdue ACTs with hasResults=false in our universe for any unposted submissionTracking / dispFirstSubmitDate!
"""
import json, os, urllib.request, urllib.parse, time
import numpy as np
import pandas as pd

BASE = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_b2_fda_2026_sweep"
OUT = os.path.join(BASE, "audit", "b2_t00_nlm_qc_api_audit.json")
os.makedirs(os.path.dirname(OUT), exist_ok=True)


def get_json(url):
    req = urllib.request.Request(url, headers={"User-Agent": "AntigravityResearchAudit/1.0"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read().decode())


def main():
    t0 = time.time()
    res = {"executed_at": time.strftime("%Y-%m-%d %H:%M:%S CST")}

    # 1. Inspect API v2 /studies/metadata
    meta = get_json("https://clinicaltrials.gov/api/v2/studies/metadata")
    matched_fields = []
    def walk_meta(nodes, prefix=""):
        for n in nodes:
            nm = n.get("name", "")
            full = f"{prefix}.{nm}" if prefix else nm
            low = full.lower()
            if any(k in low for k in ["result", "qc", "submit", "post", "unpost", "track", "disp", "delay", "waiv", "mcp", "release", "violation"]):
                matched_fields.append({"path": full, "type": n.get("type"), "piece": n.get("piece"), "title": n.get("title")})
            if n.get("children"):
                walk_meta(n["children"], full)
    walk_meta(meta)
    res["metadata_matching_fields"] = matched_fields
    print(f"[1] Matched {len(matched_fields)} metadata fields", flush=True)

    # 2. Query live counts & sample structures
    queries = {
        "sub_2026_11": 'AREA[ResultsFirstSubmitDate]RANGE[2026-11-01, 2026-11-25]',
        "sub_2026_10": 'AREA[ResultsFirstSubmitDate]RANGE[2026-10-01, 2026-10-31]',
        "sub_2026_09": 'AREA[ResultsFirstSubmitDate]RANGE[2026-09-01, 2026-09-30]',
        "sub_2026_08": 'AREA[ResultsFirstSubmitDate]RANGE[2026-08-01, 2026-08-31]',
        "sub_2026_04_07": 'AREA[ResultsFirstSubmitDate]RANGE[2026-04-01, 2026-07-31]',
        "post_2026_11": 'AREA[ResultsFirstPostDate]RANGE[2026-11-01, 2026-11-25]',
        "post_2026_10": 'AREA[ResultsFirstPostDate]RANGE[2026-10-01, 2026-10-31]',
        "sub_hasResults_false_all": 'AREA[ResultsFirstSubmitDate]RANGE[2020-01-01, 2026-11-25] AND AREA[HasResults]false',
        "sub_hasResults_true_2026_11": 'AREA[ResultsFirstSubmitDate]RANGE[2026-11-01, 2026-11-25] AND AREA[HasResults]true',
        "qc_null_but_sub_2026": 'AREA[ResultsFirstSubmitDate]RANGE[2026-01-01, 2026-11-25] AND AREA[ResultsFirstSubmitQcDate]MISSING',
        "qc_null_but_post_2026": 'AREA[ResultsFirstPostDate]RANGE[2026-01-01, 2026-11-25] AND AREA[ResultsFirstSubmitQcDate]MISSING',
        "unposted_2026_10_11": 'AREA[HasResults]false AND AREA[LastUpdatePostDate]RANGE[2026-10-01, 2026-11-25] AND (AREA[IsFDARegulatedDrug]true OR AREA[IsFDARegulatedDevice]true) AND AREA[LocationCountry]"United States" AND (AREA[OverallStatus]COMPLETED OR AREA[OverallStatus]TERMINATED) AND AREA[PrimaryCompletionDate]RANGE[2017-01-18, 2025-09-30]',
    }
    q_summary = {}
    samples = {}
    for k, q in queries.items():
        try:
            params = urllib.parse.urlencode({
                "query.term": q,
                "countTotal": "true",
                "pageSize": 25 if k in ("sub_2026_11", "sub_2026_10", "post_2026_11", "sub_hasResults_false_all", "qc_null_but_sub_2026", "qc_null_but_post_2026", "unposted_2026_10_11") else 5,
            })
            d = get_json(f"https://clinicaltrials.gov/api/v2/studies?{params}")
            studies = d.get("studies", [])
            q_summary[k] = {"totalCount": d.get("totalCount"), "returned": len(studies)}
            if len(studies) and k in ("sub_2026_11", "sub_2026_10", "post_2026_11", "sub_hasResults_false_all", "qc_null_but_sub_2026", "qc_null_but_post_2026", "unposted_2026_10_11"):
                recs = []
                for s in studies:
                    ps = s.get("protocolSection", {})
                    st = ps.get("statusModule", {})
                    ds = s.get("derivedSection", {})
                    ms = ds.get("miscInfoModule", {})
                    rs = s.get("resultsSection", {})
                    recs.append({
                        "nctId": ps.get("identificationModule", {}).get("nctId"),
                        "hasResults": s.get("hasResults"),
                        "sponsorClass": (ps.get("sponsorCollaboratorsModule", {}).get("leadSponsor") or {}).get("class"),
                        "overallStatus": st.get("overallStatus"),
                        "pcd": (st.get("primaryCompletionDateStruct") or {}).get("date"),
                        "resultsFirstSubmitDate": st.get("resultsFirstSubmitDate"),
                        "resultsFirstSubmitQcDate": st.get("resultsFirstSubmitQcDate"),
                        "resultsFirstPostDate": (st.get("resultsFirstPostDateStruct") or {}).get("date"),
                        "resultsFirstPostType": (st.get("resultsFirstPostDateStruct") or {}).get("type"),
                        "dispFirstSubmitDate": st.get("dispFirstSubmitDate"),
                        "dispFirstPostDate": (st.get("dispFirstPostDateStruct") or {}).get("date"),
                        "delayedPosting": st.get("delayedPosting"),
                        "resultsWaived": st.get("resultsWaived"),
                        "lastUpdatePostDate": (st.get("lastUpdatePostDateStruct") or {}).get("date"),
                        "submissionTracking": ms.get("submissionTracking"),
                        "hasResultsSectionKeys": sorted(list(rs.keys())) if rs else [],
                    })
                samples[k] = recs
        except Exception as e:
            q_summary[k] = {"error": str(e)[:200]}
    res["api_query_counts"] = q_summary
    res["api_samples"] = samples
    print(f"[2] API query counts: {json.dumps(q_summary)}", flush=True)

    # 3. Decompose our local dataset
    df = pd.read_parquet(os.path.join(BASE, "data", "clean", "b2_trial_universe_steady.parquet"))
    act = df[df["arm"] == "us_fdareg_nonp1"].copy()
    for c in ["results_first_submit_date", "results_first_submit_qc_date", "results_first_post_date", "primary_completion_date", "last_update_post_date"]:
        act[c + "_dt"] = pd.to_datetime(act[c], errors="coerce")

    sub_any = act[act["results_first_submit_date_dt"].notna()].copy()
    sub_any["lag_post_minus_sub"] = (sub_any["results_first_post_date_dt"] - sub_any["results_first_submit_date_dt"]).dt.days
    sub_any["lag_qc_minus_sub"] = (sub_any["results_first_submit_qc_date_dt"] - sub_any["results_first_submit_date_dt"]).dt.days
    sub_any["lag_post_minus_qc"] = (sub_any["results_first_post_date_dt"] - sub_any["results_first_submit_qc_date_dt"]).dt.days
    sub_any["is_industry"] = (sub_any["lead_sponsor_class"] == "INDUSTRY").astype(int)

    recent = sub_any[sub_any["results_first_submit_date_dt"] >= "2024-01-01"].copy()
    recent["ym"] = recent["results_first_submit_date_dt"].dt.to_period("M").astype(str)
    recent["qc_missing"] = recent["results_first_submit_qc_date_dt"].isna()
    recent["post_missing"] = recent["results_first_post_date_dt"].isna()
    recent["post_before_qc"] = recent["results_first_post_date_dt"] < recent["results_first_submit_qc_date_dt"]
    recent["post_within_30d"] = recent["lag_post_minus_sub"] <= 30

    lag_stats = {}
    for grp_name, gdf in [
        ("ALL_2024_2026", recent),
        ("INDUSTRY_2024_2025", recent[(recent["is_industry"]==1) & (recent["results_first_submit_date_dt"]<"2026-01-01")]),
        ("ACADEMIC_2024_2025", recent[(recent["is_industry"]==0) & (recent["results_first_submit_date_dt"]<"2026-01-01")]),
        ("INDUSTRY_2026_04_09", recent[(recent["is_industry"]==1) & (recent["results_first_submit_date_dt"]>="2026-04-01") & (recent["results_first_submit_date_dt"]<"2026-10-01")]),
        ("ACADEMIC_2026_04_09", recent[(recent["is_industry"]==0) & (recent["results_first_submit_date_dt"]>="2026-04-01") & (recent["results_first_submit_date_dt"]<"2026-10-01")]),
    ]:
        lag_stats[grp_name] = {
            "n": int(len(gdf)),
            "has_results_true": int(gdf["has_results"].sum()),
            "qc_date_missing": int(gdf["qc_missing"].sum()),
            "post_date_missing": int(gdf["post_missing"].sum()),
            "post_before_qc_count": int(gdf["post_before_qc"].sum()),
            "post_within_30d_pct": round(100.0 * gdf["post_within_30d"].mean(), 2) if len(gdf) else None,
            "lag_post_minus_sub_quantiles": {str(q): float(gdf["lag_post_minus_sub"].quantile(q)) for q in [0.1, 0.25, 0.5, 0.75, 0.9]} if gdf["lag_post_minus_sub"].notna().any() else None,
            "lag_qc_minus_sub_quantiles": {str(q): float(gdf["lag_qc_minus_sub"].quantile(q)) for q in [0.1, 0.25, 0.5, 0.75, 0.9]} if gdf["lag_qc_minus_sub"].notna().any() else None,
            "lag_post_minus_qc_quantiles": {str(q): float(gdf["lag_post_minus_qc"].quantile(q)) for q in [0.1, 0.25, 0.5, 0.75, 0.9]} if gdf["lag_post_minus_qc"].notna().any() else None,
        }
    res["local_dataset_lag_decomposition"] = lag_stats

    m26 = []
    for ym in [f"2026-{m:02d}" for m in range(1, 12)]:
        g = recent[recent["ym"] == ym]
        m26.append({
            "ym": ym,
            "n_total": int(len(g)),
            "n_ind": int((g["is_industry"]==1).sum()),
            "n_acad": int((g["is_industry"]==0).sum()),
            "n_has_results_false": int((~g["has_results"]).sum()),
            "n_qc_missing": int(g["qc_missing"].sum()),
            "n_post_before_qc": int(g["post_before_qc"].sum()),
            "median_post_minus_sub": float(g["lag_post_minus_sub"].median()) if len(g) else None,
            "median_qc_minus_sub": float(g["lag_qc_minus_sub"].median()) if len(g) else None,
            "median_post_minus_qc": float(g["lag_post_minus_qc"].median()) if len(g) else None,
        })
    res["local_2026_monthly_qc_decomposition"] = m26

    # 4. CRITICAL CHECK: Fetch derivedSection.miscInfoModule.submissionTracking for ALL overdue ACTs with has_results=False
    #    updated since 2026-04-01 (and specifically count ALL of them in batches of 50!)
    overdue_no_res = act[(~act["has_results"]) & (act["results_first_submit_date"].isna()) & (act["last_update_post_date_dt"] >= "2026-04-01")].copy()
    res["overdue_no_res_updated_post_mar26_total"] = int(len(overdue_no_res))
    all_ncts = overdue_no_res.sort_values("last_update_post_date_dt", ascending=False)["nct_id"].tolist()
    found_tracking = []
    # Check up to 600 most recently updated overdue no-results ACTs!
    check_ncts = all_ncts[:600]
    for i in range(0, len(check_ncts), 50):
        chunk = check_ncts[i:i+50]
        params = urllib.parse.urlencode({
            "query.id": " OR ".join(chunk),
            "fields": "protocolSection.identificationModule.nctId,protocolSection.statusModule,protocolSection.sponsorCollaboratorsModule.leadSponsor,derivedSection.miscInfoModule,hasResults",
            "pageSize": 50,
        })
        d = get_json(f"https://clinicaltrials.gov/api/v2/studies?{params}")
        for s in d.get("studies", []):
            st = s.get("protocolSection", {}).get("statusModule", {})
            ms = s.get("derivedSection", {}).get("miscInfoModule", {})
            trk = ms.get("submissionTracking")
            if trk or st.get("resultsFirstSubmitDate") or st.get("dispFirstSubmitDate"):
                found_tracking.append({
                    "nctId": s.get("protocolSection", {}).get("identificationModule", {}).get("nctId"),
                    "sponsorClass": (s.get("protocolSection", {}).get("sponsorCollaboratorsModule", {}).get("leadSponsor") or {}).get("class"),
                    "hasResults": s.get("hasResults"),
                    "resultsFirstSubmitDate": st.get("resultsFirstSubmitDate"),
                    "resultsFirstSubmitQcDate": st.get("resultsFirstSubmitQcDate"),
                    "resultsFirstPostDate": (st.get("resultsFirstPostDateStruct") or {}).get("date"),
                    "dispFirstSubmitDate": st.get("dispFirstSubmitDate"),
                    "submissionTracking": trk,
                })
    # Summarize found_tracking by Industry vs Academic and by estimatedResultsFirstSubmitDate / firstMcpInfo / submissionInfos!
    n_ind_trk = sum(1 for x in found_tracking if x["sponsorClass"] == "INDUSTRY")
    n_aca_trk = sum(1 for x in found_tracking if x["sponsorClass"] != "INDUSTRY")
    n_with_est_sub = sum(1 for x in found_tracking if (x.get("submissionTracking") or {}).get("estimatedResultsFirstSubmitDate"))
    n_with_mcp = sum(1 for x in found_tracking if (x.get("submissionTracking") or {}).get("firstMcpInfo"))
    n_with_subinfos = sum(1 for x in found_tracking if (x.get("submissionTracking") or {}).get("submissionInfos"))
    res["overdue_no_res_submissionTracking_audit"] = {
        "checked_ncts": len(check_ncts),
        "n_with_tracking_or_submit": len(found_tracking),
        "n_industry": n_ind_trk,
        "n_academic": n_aca_trk,
        "n_with_estimatedResultsFirstSubmitDate": n_with_est_sub,
        "n_with_firstMcpInfo": n_with_mcp,
        "n_with_submissionInfos": n_with_subinfos,
        "examples": found_tracking[:25],
    }
    print(f"[4] Checked {len(check_ncts)} overdue no-res ACTs updated >=2026-04-01: {len(found_tracking)} have submissionTracking/submit! (Ind={n_ind_trk}, Acad={n_aca_trk}, estSub={n_with_est_sub}, mcp={n_with_mcp})", flush=True)

    # 5. Check submissionTracking on 50 posted 2026 ACTs (25 lag<=30d, 25 lag>60d)
    with_res_26 = recent[(recent["results_first_submit_date_dt"] >= "2026-04-01")].copy()
    fast_ncts = with_res_26[with_res_26["lag_post_minus_sub"] <= 30]["nct_id"].head(25).tolist()
    slow_ncts = with_res_26[with_res_26["lag_post_minus_sub"] > 60]["nct_id"].head(25).tolist()
    check_50 = fast_ncts + slow_ncts
    res_50_tracking = []
    if check_50:
        params = urllib.parse.urlencode({
            "query.id": " OR ".join(check_50),
            "fields": "protocolSection.identificationModule.nctId,protocolSection.statusModule,protocolSection.sponsorCollaboratorsModule.leadSponsor,derivedSection.miscInfoModule,hasResults",
            "pageSize": 50,
        })
        d = get_json(f"https://clinicaltrials.gov/api/v2/studies?{params}")
        for s in d.get("studies", []):
            st = s.get("protocolSection", {}).get("statusModule", {})
            ms = s.get("derivedSection", {}).get("miscInfoModule", {})
            res_50_tracking.append({
                "nctId": s.get("protocolSection", {}).get("identificationModule", {}).get("nctId"),
                "sponsorClass": (s.get("protocolSection", {}).get("sponsorCollaboratorsModule", {}).get("leadSponsor") or {}).get("class"),
                "hasResults": s.get("hasResults"),
                "resultsFirstSubmitDate": st.get("resultsFirstSubmitDate"),
                "resultsFirstSubmitQcDate": st.get("resultsFirstSubmitQcDate"),
                "resultsFirstPostDate": (st.get("resultsFirstPostDateStruct") or {}).get("date"),
                "submissionTracking": ms.get("submissionTracking"),
            })
    res["posted_2026_fast_vs_slow_submissionTracking"] = res_50_tracking

    # 6. Probe API search AREA names for unposted/submissionTracking
    unposted_field_probes = {}
    for fname in ["SubmissionTracking", "EstimatedResultsFirstSubmitDate", "FirstMcpPostDate", "ResultsFirstSubmitDate", "ResultsFirstSubmitQcDate", "DispFirstSubmitDate", "UnpostedEvent"]:
        try:
            q = f'AREA[{fname}]RANGE[2026-04-01, 2026-11-25]'
            params = urllib.parse.urlencode({"query.term": q, "countTotal": "true", "pageSize": 3})
            d = get_json(f"https://clinicaltrials.gov/api/v2/studies?{params}")
            unposted_field_probes[fname] = {"ok": True, "totalCount": d.get("totalCount")}
        except Exception as e:
            unposted_field_probes[fname] = {"ok": False, "error": str(e)[:150]}
    res["api_search_area_probes"] = unposted_field_probes
    res["elapsed_sec"] = round(time.time() - t0, 1)

    with open(OUT, "w") as fh:
        json.dump(res, fh, indent=2)
    print(f"DONE in {res['elapsed_sec']}s -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
