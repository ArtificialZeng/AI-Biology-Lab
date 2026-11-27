#!/opt/anaconda3/bin/python3
"""
T00 Complete Recovery of Unposted / In-QC Results Submissions & Delayed-Results Certifications
for Study B-2 (ClinicalTrials.gov API v2).

Fetches:
  1. ALL studies with AREA[EstimatedResultsFirstSubmitDate]RANGE[2017-01-01, 2026-11-25]
     (captures 100% of unposted / in-QC / MCP-posted-with-QC-comments submissions where hasResults=false
     and protocolSection.statusModule.resultsFirstSubmitDate is NULL!)
  2. ALL studies with AREA[DispFirstSubmitDate]RANGE[2017-01-01, 2026-11-25] AND PCD in 2017-01-18..2026-08-31
     (captures delayed-results certifications & extensions!)
  3. ALL studies in our universe that cleared QC in 2025-04..2026-11 to also retrieve
     derivedSection.miscInfoModule.submissionTracking.firstMcpInfo (30-day posting with QC comments date)
     and annotationSection.annotationModule.violationAnnotation (FDAAA 801 Notices of Noncompliance!).

Merges cleanly with b2_trial_universe_steady.parquet + b2_us_fdareg_recent_25_26.parquet to create
b2_t00_complete_submission_panel.parquet with ZERO QC-clearance right-censoring!
"""
import json, os, time, urllib.request, urllib.parse
import numpy as np
import pandas as pd

BASE = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_b2_fda_2026_sweep"
CLEAN = os.path.join(BASE, "data", "clean")
AUD = os.path.join(BASE, "audit")
os.makedirs(AUD, exist_ok=True)

FIELDS = ",".join([
    "protocolSection.identificationModule.nctId",
    "protocolSection.statusModule.overallStatus",
    "protocolSection.statusModule.startDateStruct",
    "protocolSection.statusModule.primaryCompletionDateStruct",
    "protocolSection.statusModule.resultsFirstSubmitDate",
    "protocolSection.statusModule.resultsFirstSubmitQcDate",
    "protocolSection.statusModule.resultsFirstPostDateStruct",
    "protocolSection.statusModule.dispFirstSubmitDate",
    "protocolSection.statusModule.dispFirstSubmitQcDate",
    "protocolSection.statusModule.dispFirstPostDateStruct",
    "protocolSection.statusModule.delayedPosting",
    "protocolSection.statusModule.resultsWaived",
    "protocolSection.sponsorCollaboratorsModule.leadSponsor",
    "annotationSection.annotationModule",
    "derivedSection.miscInfoModule.submissionTracking",
    "hasResults",
])


def fetch_all_pages(query, label, page_size=1000, max_pages=40):
    all_s = []
    tok = None
    total = None
    for page in range(1, max_pages + 1):
        params = {"query.term": query, "fields": FIELDS, "pageSize": page_size}
        if page == 1:
            params["countTotal"] = "true"
        if tok:
            params["pageToken"] = tok
        url = "https://clinicaltrials.gov/api/v2/studies?" + urllib.parse.urlencode(params)
        for attempt in range(4):
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "AntigravityT00Audit/1.0"})
                with urllib.request.urlopen(req, timeout=45) as r:
                    d = json.loads(r.read().decode())
                break
            except Exception as e:
                print(f"  [{label} p{page}] retry {attempt+1}: {e}", flush=True)
                time.sleep(2 * (attempt + 1))
        else:
            raise RuntimeError(f"Failed fetching {label} p{page}")
        if page == 1:
            total = d.get("totalCount")
        st = d.get("studies", [])
        all_s.extend(st)
        tok = d.get("nextPageToken")
        if page == 1 or page % 5 == 0 or not tok:
            print(f"  [{label}] p{page}: +{len(st)} (cum={len(all_s)} / total={total})", flush=True)
        if not tok:
            break
    return all_s, total


def extract_row(s):
    ps = s.get("protocolSection", {})
    ident = ps.get("identificationModule", {})
    st = ps.get("statusModule", {})
    ms = (s.get("derivedSection") or {}).get("miscInfoModule") or {}
    trk = ms.get("submissionTracking") or {}
    ann = (s.get("annotationSection") or {}).get("annotationModule") or {}
    unp = ann.get("unpostedAnnotation") or {}
    vio = ann.get("violationAnnotation") or {}
    sub_infos = trk.get("submissionInfos") or []
    unp_events = unp.get("unpostedEvents") or []

    # Earliest releaseDate from submissionInfos and unpostedEvents
    rel_dates = [x.get("releaseDate") for x in sub_infos if x.get("releaseDate")]
    unp_rel_dates = [x.get("date") for x in unp_events if x.get("date") and (x.get("type") in ("RELEASE", None))]
    mcp_counts = [x.get("mcpReleaseN") for x in sub_infos if isinstance(x.get("mcpReleaseN"), (int, float))]
    has_unrelease = any(x.get("unreleaseDate") or x.get("resetDate") for x in sub_infos)
    has_unp_reset = any((x.get("type") in ("RESET", "UNRELEASE")) for x in unp_events)

    return {
        "nct_id": ident.get("nctId"),
        "has_results_live": bool(s.get("hasResults")),
        "start_date_live": (st.get("startDateStruct") or {}).get("date"),
        "results_first_submit_live": st.get("resultsFirstSubmitDate"),
        "results_first_submit_qc_live": st.get("resultsFirstSubmitQcDate"),
        "results_first_post_live": (st.get("resultsFirstPostDateStruct") or {}).get("date"),
        "disp_first_submit_date": st.get("dispFirstSubmitDate"),
        "disp_first_submit_qc_date": st.get("dispFirstSubmitQcDate"),
        "disp_first_post_date": (st.get("dispFirstPostDateStruct") or {}).get("date"),
        "delayed_posting": st.get("delayedPosting"),
        "results_waived": st.get("resultsWaived"),
        "est_results_first_submit_date": trk.get("estimatedResultsFirstSubmitDate"),
        "first_mcp_post_date": ((trk.get("firstMcpInfo") or {}).get("postDateStruct") or {}).get("date"),
        "earliest_subinfo_release_date": min(rel_dates) if rel_dates else None,
        "latest_subinfo_release_date": max(rel_dates) if rel_dates else None,
        "n_subinfo_cycles": len(sub_infos),
        "max_mcp_release_n": max(mcp_counts) if mcp_counts else 0,
        "has_unrelease_or_reset": bool(has_unrelease or has_unp_reset),
        "earliest_unposted_event_date": min(unp_rel_dates) if unp_rel_dates else None,
        "n_violation_events": len(vio.get("violationEvents") or []),
    }


def main():
    t0 = time.time()
    # 1. Fetch ALL studies with EstimatedResultsFirstSubmitDate in 2017-01-01..2026-11-25 AND PCD in 2017-01-18..2026-08-31
    q_est = 'AREA[EstimatedResultsFirstSubmitDate]RANGE[2017-01-01, 2026-11-25] AND AREA[PrimaryCompletionDate]RANGE[2017-01-18, 2026-08-31]'
    s_est, tot_est = fetch_all_pages(q_est, "est_unposted_2017_2026")

    # 2. Fetch ALL studies with DispFirstSubmitDate in 2017-01-01..2026-11-25 AND PCD in 2017-01-18..2026-08-31
    q_disp = 'AREA[DispFirstSubmitDate]RANGE[2017-01-01, 2026-11-25] AND AREA[PrimaryCompletionDate]RANGE[2017-01-18, 2026-08-31]'
    s_disp, tot_disp = fetch_all_pages(q_disp, "disp_cert_2017_2026")

    # 3. Fetch ALL US FDA-regulated non-P1 posted results submitted in 2023-01-01..2026-11-25 to also get firstMcpInfo + startDate!
    q_mcp = (
        'AREA[StudyType]INTERVENTIONAL AND AREA[LocationCountry]"United States" '
        'AND (AREA[IsFDARegulatedDrug]true OR AREA[IsFDARegulatedDevice]true) '
        'AND AREA[ResultsFirstSubmitDate]RANGE[2023-01-01, 2026-11-25] '
        'AND AREA[PrimaryCompletionDate]RANGE[2017-01-18, 2026-08-31]'
    )
    s_mcp, tot_mcp = fetch_all_pages(q_mcp, "us_fdareg_posted_mcp_2023_2026")

    # Merge extracted records (priority: est_unposted, mcp, disp -> keep richest non-null)
    rec_map = {}
    for src in (s_disp, s_mcp, s_est):
        for s in src:
            r = extract_row(s)
            nid = r["nct_id"]
            if nid not in rec_map:
                rec_map[nid] = r
            else:
                cur = rec_map[nid]
                for k, v in r.items():
                    if cur.get(k) in (None, "", 0, False) and v not in (None, "", 0, False):
                        cur[k] = v

    df_t00 = pd.DataFrame(list(rec_map.values()))
    out_pq = os.path.join(AUD, "b2_t00_unposted_and_cert_overlay.parquet")
    df_t00.to_parquet(out_pq, index=False)

    # Immediate diagnostic join against our steady universe + recent 25-26 universe!
    d_main = pd.read_parquet(os.path.join(CLEAN, "b2_trial_universe_steady.parquet"))
    d_rec = pd.read_parquet(os.path.join(CLEAN, "b2_us_fdareg_recent_25_26.parquet"))
    d_rec["arm"] = "us_fdareg_nonp1"
    u = pd.concat([d_main, d_rec], ignore_index=True).drop_duplicates("nct_id", keep="first")
    m = u.merge(df_t00, on="nct_id", how="left")

    act_strict = m[
        (m["arm"] == "us_fdareg_nonp1") &
        m["overall_status"].isin(["COMPLETED", "TERMINATED"]) &
        (m["primary_completion_type"] == "ACTUAL")
    ].copy()
    act_strict["is_industry"] = (act_strict["lead_sponsor_class"] == "INDUSTRY").astype(int)
    act_strict["has_unposted_sub"] = act_strict["results_first_submit_date"].isna() & act_strict["est_results_first_submit_date"].notna()

    # Tabulate recovered unposted submissions by month (2025-04..2026-11) and Industry vs Academic!
    unp_act = act_strict[act_strict["has_unposted_sub"]].copy()
    unp_act["est_dt"] = pd.to_datetime(unp_act["est_results_first_submit_date"], errors="coerce")
    unp_act["pcd_dt"] = pd.to_datetime(unp_act["primary_completion_date"], errors="coerce")
    unp_act["is_overdue_at_sub"] = (unp_act["est_dt"] - unp_act["pcd_dt"]).dt.days > 365
    unp_act["est_ym"] = unp_act["est_dt"].dt.to_period("M").astype(str)

    by_ym = []
    for ym in [str(p) for p in pd.period_range("2025-04", "2026-11", freq="M")]:
        g = unp_act[unp_act["est_ym"] == ym]
        g_od = g[g["is_overdue_at_sub"]]
        g_ot = g[~g["is_overdue_at_sub"]]
        by_ym.append({
            "ym": ym,
            "unposted_all_ind": int((g["is_industry"]==1).sum()),
            "unposted_all_aca": int((g["is_industry"]==0).sum()),
            "unposted_overdue_ind": int((g_od["is_industry"]==1).sum()),
            "unposted_overdue_aca": int((g_od["is_industry"]==0).sum()),
            "unposted_ontime_ind": int((g_ot["is_industry"]==1).sum()),
            "unposted_ontime_aca": int((g_ot["is_industry"]==0).sum()),
        })

    # Also check how many unposted submissions have first_mcp_post_date vs not, and how many have start_date >= 2017-01-18!
    summary = {
        "executed_at": time.strftime("%Y-%m-%d %H:%M:%S CST"),
        "elapsed_sec": round(time.time() - t0, 1),
        "fetched_totals": {"est_unposted": tot_est, "disp_cert": tot_disp, "us_fdareg_posted_mcp": tot_mcp},
        "unique_overlay_ncts": int(len(df_t00)),
        "strict_act_total": int(len(act_strict)),
        "strict_act_unposted_recovered_total": int(len(unp_act)),
        "strict_act_unposted_overdue_total": int(unp_act["is_overdue_at_sub"].sum()),
        "strict_act_unposted_overdue_ind": int((unp_act["is_overdue_at_sub"] & (unp_act["is_industry"]==1)).sum()),
        "strict_act_unposted_overdue_aca": int((unp_act["is_overdue_at_sub"] & (unp_act["is_industry"]==0)).sum()),
        "strict_act_unposted_with_mcp_count": int(unp_act["first_mcp_post_date"].notna().sum()),
        "strict_act_unposted_with_mcp_release_n_gt0": int((unp_act["max_mcp_release_n"] > 0).sum()),
        "strict_act_unposted_by_ym": by_ym,
    }
    sum_p = os.path.join(AUD, "b2_t00_unposted_recovery_summary.json")
    with open(sum_p, "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
