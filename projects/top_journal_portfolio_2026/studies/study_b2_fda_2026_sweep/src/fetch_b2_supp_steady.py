#!/opt/anaconda3/bin/python3
"""
Fetch supplementary rolling-entry cohorts (PrimaryCompletionDate 2025-03-01..2025-09-30)
AND ACTIVE_NOT_RECRUITING trials (PCD 2017-01-18..2025-09-30) across all 4 arms
so that Study B-2 has a 100% stationary rolling-entry delinquent risk set through 2026-09-30
(zero artificial cohort-aging drop!) AND captures active-but-primary-completed FDAAA trials.
"""
import json, os, subprocess, time
import pandas as pd
import sys
sys.path.insert(0, os.path.dirname(__file__))
from fetch_b2_universe import PROTOCOL_FIELDS, fetch_paginated, flatten_protocol_study, CLEAN_DIR

SUPP_ARMS = {
    "us_fdareg_nonp1_25": (
        'AREA[StudyType]INTERVENTIONAL AND AREA[LocationCountry]"United States" '
        'AND (AREA[IsFDARegulatedDrug]true OR AREA[IsFDARegulatedDevice]true) '
        'AND NOT ((AREA[Phase]PHASE1 OR AREA[Phase]EARLY_PHASE1) AND NOT (AREA[Phase]PHASE2 OR AREA[Phase]PHASE3 OR AREA[Phase]PHASE4)) '
        'AND (AREA[OverallStatus]COMPLETED OR AREA[OverallStatus]TERMINATED OR AREA[OverallStatus]ACTIVE_NOT_RECRUITING) '
        'AND AREA[PrimaryCompletionDate]RANGE[2025-03-01, 2025-09-30]'
    ),
    "us_fdareg_nonp1_anr": (
        'AREA[StudyType]INTERVENTIONAL AND AREA[LocationCountry]"United States" '
        'AND (AREA[IsFDARegulatedDrug]true OR AREA[IsFDARegulatedDevice]true) '
        'AND NOT ((AREA[Phase]PHASE1 OR AREA[Phase]EARLY_PHASE1) AND NOT (AREA[Phase]PHASE2 OR AREA[Phase]PHASE3 OR AREA[Phase]PHASE4)) '
        'AND AREA[OverallStatus]ACTIVE_NOT_RECRUITING '
        'AND AREA[PrimaryCompletionDate]RANGE[2017-01-18, 2025-02-28]'
    ),
    "us_notfdareg_25": (
        'AREA[StudyType]INTERVENTIONAL AND AREA[LocationCountry]"United States" '
        'AND AREA[IsFDARegulatedDrug]false AND AREA[IsFDARegulatedDevice]false '
        'AND (AREA[OverallStatus]COMPLETED OR AREA[OverallStatus]TERMINATED OR AREA[OverallStatus]ACTIVE_NOT_RECRUITING) '
        'AND AREA[PrimaryCompletionDate]RANGE[2025-03-01, 2025-09-30]'
    ),
    "us_p1_only_25": (
        'AREA[StudyType]INTERVENTIONAL AND AREA[LocationCountry]"United States" '
        'AND (AREA[Phase]PHASE1 OR AREA[Phase]EARLY_PHASE1) '
        'AND NOT (AREA[Phase]PHASE2 OR AREA[Phase]PHASE3 OR AREA[Phase]PHASE4) '
        'AND (AREA[OverallStatus]COMPLETED OR AREA[OverallStatus]TERMINATED OR AREA[OverallStatus]ACTIVE_NOT_RECRUITING) '
        'AND AREA[PrimaryCompletionDate]RANGE[2025-03-01, 2025-09-30]'
    ),
    "nonus_p23_nonfda_25": (
        'AREA[StudyType]INTERVENTIONAL AND NOT AREA[LocationCountry]"United States" '
        'AND (AREA[Phase]PHASE2 OR AREA[Phase]PHASE3) '
        'AND AREA[IsFDARegulatedDrug]false AND AREA[IsFDARegulatedDevice]false '
        'AND (AREA[OverallStatus]COMPLETED OR AREA[OverallStatus]TERMINATED OR AREA[OverallStatus]ACTIVE_NOT_RECRUITING) '
        'AND AREA[PrimaryCompletionDate]RANGE[2025-03-01, 2025-09-30]'
    ),
}

ARM_MAP = {
    "us_fdareg_nonp1_25": "us_fdareg_nonp1",
    "us_fdareg_nonp1_anr": "us_fdareg_nonp1",
    "us_notfdareg_25": "us_notfdareg",
    "us_p1_only_25": "us_p1_only",
    "nonus_p23_nonfda_25": "nonus_p23_nonfda",
}

def main():
    t0 = time.time()
    supp_rows = []
    meta = {}
    for k, q in SUPP_ARMS.items():
        studies, exp = fetch_paginated(k, q, PROTOCOL_FIELDS, page_size=1000, max_pages=20)
        rows = [flatten_protocol_study(s, ARM_MAP[k]) for s in studies]
        supp_rows.extend(rows)
        meta[k] = {"expected": exp, "fetched": len(rows)}
    df_base = pd.read_parquet(os.path.join(CLEAN_DIR, "b2_trial_universe.parquet"))
    df_supp = pd.DataFrame(supp_rows)
    df_all = pd.concat([df_base, df_supp], ignore_index=True)
    arm_prio = {"us_fdareg_nonp1": 0, "us_p1_only": 1, "us_notfdareg": 2, "nonus_p23_nonfda": 3}
    df_all["_p"] = df_all["arm"].map(arm_prio).fillna(9)
    df_all = df_all.sort_values("_p").drop_duplicates("nct_id", keep="first").drop(columns=["_p"]).reset_index(drop=True)
    out_pq = os.path.join(CLEAN_DIR, "b2_trial_universe_steady.parquet")
    df_all.to_parquet(out_pq, index=False)
    print(f"Merged steady-state universe: {len(df_all):,} unique trials in {time.time()-t0:.1f}s", flush=True)
    print(json.dumps(meta, indent=2), flush=True)

if __name__ == "__main__":
    main()
