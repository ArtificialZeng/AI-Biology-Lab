#!/opt/anaconda3/bin/python3
"""
Study B-2 Prospective Validation & Exact Risk-Set On-Time Denominator Pull:
  1. Pull US FDA-regulated non-P1 trials with PCD 2025-10-01..2026-08-31 so we have EXACT trial-month
     risk-set denominators for BOTH On-Time (<=365d post-PCD) and Overdue (>365d post-PCD) trials
     throughout 2025-04..2026-07 (eliminating any count-composition critique!).
  2. Prospectively locked UNSEEN legacy validation cohort (never inspected before):
     Pre-Final-Rule US Phase 2/3 Drug/Biologic/Device trials (PCD 2012-01-01..2017-01-17)
     to test our "Stale-Data Decay" prediction out-of-sample!
"""
import json, os, sys, time
import pandas as pd

BASE = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_b2_fda_2026_sweep"
CLEAN = os.path.join(BASE, "data", "clean")
sys.path.insert(0, os.path.join(BASE, "src"))
from fetch_b2_universe import PROTOCOL_FIELDS, fetch_paginated, flatten_protocol_study

Q_ONTIME_RECENT = (
    'AREA[StudyType]INTERVENTIONAL AND AREA[LocationCountry]"United States" '
    'AND (AREA[IsFDARegulatedDrug]true OR AREA[IsFDARegulatedDevice]true) '
    'AND NOT ((AREA[Phase]PHASE1 OR AREA[Phase]EARLY_PHASE1) AND NOT (AREA[Phase]PHASE2 OR AREA[Phase]PHASE3 OR AREA[Phase]PHASE4)) '
    'AND (AREA[OverallStatus]COMPLETED OR AREA[OverallStatus]TERMINATED) '
    'AND AREA[PrimaryCompletionDate]RANGE[2025-10-01, 2026-08-31]'
)

Q_UNSEEN_LEGACY = (
    'AREA[StudyType]INTERVENTIONAL AND AREA[LocationCountry]"United States" '
    'AND (AREA[Phase]PHASE2 OR AREA[Phase]PHASE3) '
    'AND (AREA[InterventionType]DRUG OR AREA[InterventionType]BIOLOGICAL OR AREA[InterventionType]DEVICE) '
    'AND (AREA[OverallStatus]COMPLETED OR AREA[OverallStatus]TERMINATED) '
    'AND AREA[PrimaryCompletionDate]RANGE[2012-01-01, 2017-01-17]'
)

def main():
    t0 = time.time()
    print("=== 1. Fetching 2025-10..2026-08 US FDA-reg completions for exact On-Time risk-set denominator ===", flush=True)
    s1, exp1 = fetch_paginated("us_fdareg_recent_25_26", Q_ONTIME_RECENT, PROTOCOL_FIELDS, page_size=1000, max_pages=10)
    df1 = pd.DataFrame([flatten_protocol_study(s, "us_fdareg_nonp1") for s in s1])
    df1.to_parquet(os.path.join(CLEAN, "b2_us_fdareg_recent_25_26.parquet"), index=False)

    print("=== 2. Fetching UNSEEN 2012-2017 pre-Final-Rule legacy FDAAA validation cohort ===", flush=True)
    s2, exp2 = fetch_paginated("unseen_legacy_2012_2017", Q_UNSEEN_LEGACY, PROTOCOL_FIELDS, page_size=1000, max_pages=25)
    df2 = pd.DataFrame([flatten_protocol_study(s, "unseen_legacy_12_17") for s in s2])
    df2.to_parquet(os.path.join(CLEAN, "b2_unseen_legacy_12_17.parquet"), index=False)

    print(json.dumps({
        "recent_25_26": {"expected": exp1, "fetched": len(df1)},
        "unseen_legacy_12_17": {"expected": exp2, "fetched": len(df2)},
        "elapsed_sec": round(time.time() - t0, 1)
    }, indent=2), flush=True)

if __name__ == "__main__":
    main()
