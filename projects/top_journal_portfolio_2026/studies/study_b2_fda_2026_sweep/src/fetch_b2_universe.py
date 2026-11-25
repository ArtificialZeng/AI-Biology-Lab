#!/opt/anaconda3/bin/python3
"""
Fetch full ClinicalTrials.gov trial-level universe for Study B-2
using the official clinical-trials-database skill wrapper (clinical_trials_api.py).

Arms pulled (PrimaryCompletionDate 2017-01-18..2025-02-28, INTERVENTIONAL, COMPLETED|TERMINATED):
  1. us_fdareg_nonp1 : US location + (IsFDARegulatedDrug=true OR IsFDARegulatedDevice=true) + NOT Phase1-only
  2. us_notfdareg    : US location + IsFDARegulatedDrug=false + IsFDARegulatedDevice=false
  3. us_p1_only      : US location + Phase=PHASE1/EARLY_PHASE1 + NOT Phase2/3/4
  4. nonus_p23_nonfda: NOT US location + (Phase2 OR Phase3) + IsFDARegulatedDrug=false + IsFDARegulatedDevice=false

Plus full resultsSection (outcomeMeasuresModule + adverseEventsModule) for all us_fdareg_nonp1
trials with ResultsFirstSubmitDate in 2024-01-01..2026-10-31 (both overdue and on-time comparators)
for the H4 Lemons / Adverse-Selection analysis!
"""
import json
import os
import subprocess
import sys
import time
import hashlib
import pandas as pd

API = "/Users/mac/.gemini/config/plugins/scientific-literature/skills/clinical_trials_database/scripts/clinical_trials_api.py"
BASE = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_b2_fda_2026_sweep"
RAW_DIR = os.path.join(BASE, "data", "raw")
CLEAN_DIR = os.path.join(BASE, "data", "clean")
os.makedirs(RAW_DIR, exist_ok=True)
os.makedirs(CLEAN_DIR, exist_ok=True)

PROTOCOL_FIELDS = ",".join([
    "protocolSection.identificationModule.nctId",
    "protocolSection.identificationModule.briefTitle",
    "protocolSection.identificationModule.organization",
    "protocolSection.statusModule.overallStatus",
    "protocolSection.statusModule.whyStopped",
    "protocolSection.statusModule.startDateStruct",
    "protocolSection.statusModule.primaryCompletionDateStruct",
    "protocolSection.statusModule.completionDateStruct",
    "protocolSection.statusModule.resultsFirstSubmitDate",
    "protocolSection.statusModule.resultsFirstSubmitQcDate",
    "protocolSection.statusModule.resultsFirstPostDateStruct",
    "protocolSection.statusModule.studyFirstPostDateStruct",
    "protocolSection.statusModule.lastUpdatePostDateStruct",
    "protocolSection.sponsorCollaboratorsModule.leadSponsor",
    "protocolSection.sponsorCollaboratorsModule.responsibleParty",
    "protocolSection.oversightModule",
    "protocolSection.designModule.studyType",
    "protocolSection.designModule.phases",
    "protocolSection.designModule.enrollmentInfo",
    "protocolSection.designModule.designInfo",
    "protocolSection.conditionsModule.conditions",
    "protocolSection.armsInterventionsModule.interventions.type",
    "hasResults",
])

RESULTS_FIELDS = ",".join([
    "protocolSection.identificationModule.nctId",
    "protocolSection.statusModule.overallStatus",
    "protocolSection.statusModule.whyStopped",
    "protocolSection.statusModule.primaryCompletionDateStruct",
    "protocolSection.statusModule.resultsFirstSubmitDate",
    "protocolSection.statusModule.resultsFirstSubmitQcDate",
    "protocolSection.statusModule.resultsFirstPostDateStruct",
    "protocolSection.sponsorCollaboratorsModule.leadSponsor",
    "protocolSection.oversightModule",
    "protocolSection.designModule.phases",
    "protocolSection.designModule.enrollmentInfo",
    "protocolSection.conditionsModule.conditions",
    "protocolSection.armsInterventionsModule.interventions.type",
    "resultsSection.outcomeMeasuresModule.outcomeMeasures.type",
    "resultsSection.outcomeMeasuresModule.outcomeMeasures.analyses.pValue",
    "resultsSection.adverseEventsModule.eventGroups",
    "hasResults",
])

ARMS = {
    "us_fdareg_nonp1": (
        'AREA[StudyType]INTERVENTIONAL AND AREA[LocationCountry]"United States" '
        'AND (AREA[IsFDARegulatedDrug]true OR AREA[IsFDARegulatedDevice]true) '
        'AND NOT ((AREA[Phase]PHASE1 OR AREA[Phase]EARLY_PHASE1) AND NOT (AREA[Phase]PHASE2 OR AREA[Phase]PHASE3 OR AREA[Phase]PHASE4)) '
        'AND (AREA[OverallStatus]COMPLETED OR AREA[OverallStatus]TERMINATED) '
        'AND AREA[PrimaryCompletionDate]RANGE[2017-01-18, 2025-02-28]'
    ),
    "us_notfdareg": (
        'AREA[StudyType]INTERVENTIONAL AND AREA[LocationCountry]"United States" '
        'AND AREA[IsFDARegulatedDrug]false AND AREA[IsFDARegulatedDevice]false '
        'AND (AREA[OverallStatus]COMPLETED OR AREA[OverallStatus]TERMINATED) '
        'AND AREA[PrimaryCompletionDate]RANGE[2017-01-18, 2025-02-28]'
    ),
    "us_p1_only": (
        'AREA[StudyType]INTERVENTIONAL AND AREA[LocationCountry]"United States" '
        'AND (AREA[Phase]PHASE1 OR AREA[Phase]EARLY_PHASE1) '
        'AND NOT (AREA[Phase]PHASE2 OR AREA[Phase]PHASE3 OR AREA[Phase]PHASE4) '
        'AND (AREA[OverallStatus]COMPLETED OR AREA[OverallStatus]TERMINATED) '
        'AND AREA[PrimaryCompletionDate]RANGE[2017-01-18, 2025-02-28]'
    ),
    "nonus_p23_nonfda": (
        'AREA[StudyType]INTERVENTIONAL AND NOT AREA[LocationCountry]"United States" '
        'AND (AREA[Phase]PHASE2 OR AREA[Phase]PHASE3) '
        'AND AREA[IsFDARegulatedDrug]false AND AREA[IsFDARegulatedDevice]false '
        'AND (AREA[OverallStatus]COMPLETED OR AREA[OverallStatus]TERMINATED) '
        'AND AREA[PrimaryCompletionDate]RANGE[2017-01-18, 2025-02-28]'
    ),
}

LEMONS_QUERY = (
    'AREA[StudyType]INTERVENTIONAL AND AREA[LocationCountry]"United States" '
    'AND (AREA[IsFDARegulatedDrug]true OR AREA[IsFDARegulatedDevice]true) '
    'AND NOT ((AREA[Phase]PHASE1 OR AREA[Phase]EARLY_PHASE1) AND NOT (AREA[Phase]PHASE2 OR AREA[Phase]PHASE3 OR AREA[Phase]PHASE4)) '
    'AND (AREA[OverallStatus]COMPLETED OR AREA[OverallStatus]TERMINATED) '
    'AND AREA[PrimaryCompletionDate]RANGE[2017-01-18, 2026-09-30] '
    'AND AREA[ResultsFirstSubmitDate]RANGE[2024-01-01, 2026-10-31]'
)


def fetch_paginated(arm_name, query, fields, page_size=1000, max_pages=60):
    arm_dir = os.path.join(RAW_DIR, arm_name)
    os.makedirs(arm_dir, exist_ok=True)
    all_studies = []
    token = None
    total_expected = None
    page = 0
    while page < max_pages:
        page += 1
        out_file = os.path.join(arm_dir, f"page_{page:03d}.json")
        if os.path.exists(out_file) and os.path.getsize(out_file) > 100:
            with open(out_file) as f:
                data = json.load(f)
        else:
            cmd = [
                "uv", "run", API, "search",
                "--advanced", query,
                "--fields", fields,
                "--limit", str(page_size),
                "--output", out_file,
            ]
            if page == 1:
                cmd.append("--count-total")
            if token:
                cmd.extend(["--page-token", token])
            for attempt in range(4):
                p = subprocess.run(cmd, capture_output=True, text=True)
                if p.returncode == 0 and os.path.exists(out_file):
                    break
                print(f"[{arm_name} p{page}] retry {attempt+1}: {p.stderr.strip()[:160]}", flush=True)
                time.sleep(2.0 * (attempt + 1))
            else:
                raise RuntimeError(f"Failed fetching {arm_name} page {page}: {p.stderr}")
            with open(out_file) as f:
                data = json.load(f)
        if page == 1:
            total_expected = data.get("totalCount")
        studies = data.get("studies", [])
        all_studies.extend(studies)
        token = data.get("nextPageToken")
        if page == 1 or page % 5 == 0 or not token:
            print(f"  [{arm_name}] page {page}: +{len(studies)} (cum={len(all_studies)} / total={total_expected})", flush=True)
        if not token:
            break
    return all_studies, total_expected


def flatten_protocol_study(s, arm_name):
    ps = s.get("protocolSection", {})
    ident = ps.get("identificationModule", {})
    stat = ps.get("statusModule", {})
    spon = ps.get("sponsorCollaboratorsModule", {})
    lead = spon.get("leadSponsor", {})
    rp = spon.get("responsibleParty", {})
    over = ps.get("oversightModule", {})
    des = ps.get("designModule", {})
    enr = des.get("enrollmentInfo", {})
    conds = ps.get("conditionsModule", {}).get("conditions", [])
    intrs = ps.get("armsInterventionsModule", {}).get("interventions", [])
    intr_types = sorted(set(x.get("type", "") for x in intrs if x.get("type")))
    phases = des.get("phases", [])
    return {
        "arm": arm_name,
        "nct_id": ident.get("nctId"),
        "org_name": ident.get("organization", {}).get("fullName"),
        "org_class": ident.get("organization", {}).get("class"),
        "overall_status": stat.get("overallStatus"),
        "why_stopped": stat.get("whyStopped"),
        "start_date": (stat.get("startDateStruct") or {}).get("date"),
        "primary_completion_date": (stat.get("primaryCompletionDateStruct") or {}).get("date"),
        "primary_completion_type": (stat.get("primaryCompletionDateStruct") or {}).get("type"),
        "completion_date": (stat.get("completionDateStruct") or {}).get("date"),
        "results_first_submit_date": stat.get("resultsFirstSubmitDate"),
        "results_first_submit_qc_date": stat.get("resultsFirstSubmitQcDate"),
        "results_first_post_date": (stat.get("resultsFirstPostDateStruct") or {}).get("date"),
        "study_first_post_date": (stat.get("studyFirstPostDateStruct") or {}).get("date"),
        "last_update_post_date": (stat.get("lastUpdatePostDateStruct") or {}).get("date"),
        "has_results": bool(s.get("hasResults")),
        "lead_sponsor_name": lead.get("name"),
        "lead_sponsor_class": lead.get("class"),
        "responsible_party_type": rp.get("type"),
        "is_fda_drug": over.get("isFdaRegulatedDrug"),
        "is_fda_device": over.get("isFdaRegulatedDevice"),
        "is_us_export": over.get("isUsExport"),
        "is_unapproved_device": over.get("isUnapprovedDevice"),
        "fdaaa801_violation": over.get("fdaaa801Violation"),
        "has_dmc": over.get("oversightHasDmc"),
        "phases": "|".join(phases) if phases else "",
        "enrollment": enr.get("count"),
        "enrollment_type": enr.get("type"),
        "allocation": (des.get("designInfo") or {}).get("allocation"),
        "masking": ((des.get("designInfo") or {}).get("maskingInfo") or {}).get("masking"),
        "n_conditions": len(conds),
        "primary_condition": conds[0] if conds else "",
        "intervention_types": "|".join(intr_types),
    }


def parse_p_value(raw):
    if raw is None:
        return None
    s = str(raw).strip().lower()
    if not s or s in ("na", "n/a", "ns", "not applicable", "not calculated", "none"):
        return None
    s = s.replace("p", "").replace("=", "").replace(" ", "").replace(",", ".")
    lt = False
    gt = False
    if s.startswith("<="):
        lt = True
        s = s[2:]
    elif s.startswith("<"):
        lt = True
        s = s[1:]
    elif s.startswith(">="):
        gt = True
        s = s[2:]
    elif s.startswith(">"):
        gt = True
        s = s[1:]
    try:
        v = float(s)
        if v < 0 or v > 1.0:
            return None
        if lt and v == 0.05:
            return 0.049999
        if gt and v == 0.05:
            return 0.050001
        return v
    except Exception:
        return None


def flatten_lemons_study(s):
    row = flatten_protocol_study(s, "lemons_us_fdareg")
    rs = s.get("resultsSection", {}) or {}
    oms = (rs.get("outcomeMeasuresModule") or {}).get("outcomeMeasures", []) or []
    prim_pvals = []
    n_prim_outcomes = 0
    n_prim_analyses = 0
    for om in oms:
        if (om.get("type") or "").upper() == "PRIMARY":
            n_prim_outcomes += 1
            for an in (om.get("analyses") or []):
                n_prim_analyses += 1
                pv = parse_p_value(an.get("pValue"))
                if pv is not None:
                    prim_pvals.append(pv)
    row["n_primary_outcomes"] = n_prim_outcomes
    row["n_primary_analyses"] = n_prim_analyses
    row["n_primary_pvals_parsed"] = len(prim_pvals)
    row["has_primary_pval"] = len(prim_pvals) > 0
    row["min_primary_pval"] = min(prim_pvals) if prim_pvals else None
    row["any_sig_primary_p05"] = (min(prim_pvals) < 0.05) if prim_pvals else None
    row["all_nonsig_primary_p05"] = (min(prim_pvals) >= 0.05) if prim_pvals else None

    egs = (rs.get("adverseEventsModule") or {}).get("eventGroups", []) or []
    sae_aff = 0
    sae_risk = 0
    d_aff = 0
    d_risk = 0
    for eg in egs:
        sa = eg.get("seriousNumAffected")
        sr = eg.get("seriousNumAtRisk")
        if isinstance(sa, (int, float)) and isinstance(sr, (int, float)) and sr > 0:
            sae_aff += sa
            sae_risk += sr
        da = eg.get("deathsNumAffected")
        dr = eg.get("deathsNumAtRisk")
        if isinstance(da, (int, float)) and isinstance(dr, (int, float)) and dr > 0:
            d_aff += da
            d_risk += dr
    row["sae_affected"] = sae_aff
    row["sae_at_risk"] = sae_risk
    row["sae_rate"] = (sae_aff / sae_risk) if sae_risk > 0 else None
    row["deaths_affected"] = d_aff
    row["deaths_at_risk"] = d_risk
    row["mortality_rate"] = (d_aff / d_risk) if d_risk > 0 else None
    return row


def main():
    t0 = time.time()
    meta = {"arms": {}, "started_at": time.strftime("%Y-%m-%d %H:%M:%S")}
    all_rows = []
    for arm_name, q in ARMS.items():
        print(f"=== Fetching arm: {arm_name} ===", flush=True)
        studies, expected = fetch_paginated(arm_name, q, PROTOCOL_FIELDS, page_size=1000, max_pages=60)
        rows = [flatten_protocol_study(s, arm_name) for s in studies]
        all_rows.extend(rows)
        meta["arms"][arm_name] = {"expected": expected, "fetched": len(rows)}

    df = pd.DataFrame(all_rows)
    # Check overlap across arms (should be zero by construction)
    n_unique = df["nct_id"].nunique()
    print(f"Total universe rows: {len(df)}, unique NCTs: {n_unique}", flush=True)
    univ_path = os.path.join(CLEAN_DIR, "b2_trial_universe.parquet")
    df.to_parquet(univ_path, index=False)

    print("=== Fetching Lemons resultsSection cohort (2024-2026 US FDA-reg submitters) ===", flush=True)
    lem_studies, lem_expected = fetch_paginated("lemons_results", LEMONS_QUERY, RESULTS_FIELDS, page_size=500, max_pages=40)
    lem_rows = [flatten_lemons_study(s) for s in lem_studies]
    df_lem = pd.DataFrame(lem_rows)
    lem_path = os.path.join(CLEAN_DIR, "b2_lemons_results.parquet")
    df_lem.to_parquet(lem_path, index=False)
    meta["lemons"] = {"expected": lem_expected, "fetched": len(df_lem)}
    meta["universe_rows"] = len(df)
    meta["universe_unique_nct"] = int(n_unique)
    meta["elapsed_sec"] = round(time.time() - t0, 1)

    for p in (univ_path, lem_path):
        h = hashlib.sha256(open(p, "rb").read()).hexdigest()
        meta[os.path.basename(p) + "_sha256"] = h

    meta_path = os.path.join(CLEAN_DIR, "fetch_manifest.json")
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)
    print(json.dumps(meta, indent=2), flush=True)


if __name__ == "__main__":
    main()
