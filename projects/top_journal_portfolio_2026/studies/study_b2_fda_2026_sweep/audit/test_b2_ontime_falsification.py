#!/opt/anaconda3/bin/python3
"""
THE ULTIMATE SMOKING-GUN FALSIFICATION TEST FOR STUDY B-2:
  1. Compare ON-TIME (non-delinquent, delay <= 365d) vs OVERDUE (delinquent, delay > 365d)
     monthly submissions by Industry vs Academic across 2025-04..2026-09!
     Both pass through the EXACT same NLM QC queue and share the exact same 2026 academic/industry macro shocks,
     yet ONLY OVERDUE trials received the March 30, 2026 FDA delinquency sweep!
  2. Test whether Academic's 2026-06..08 drop also occurs in ON-TIME (non-overdue) submissions
     (which would prove differential NLM QC right-censoring or academic macro shock rather than sweep effect!).
  3. Also test PROSPECTIVE UNSEEN VALIDATION on the PRE-FINAL-RULE 2012-2016 overdue cohort (never touched in pilots!)
     AND test Industry vs Academic in exact 2026-04..05 (180+ days QC buffer) vs 2026-06..08!
"""
import json, os, subprocess, warnings
import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.api as sm
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")

BASE = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_b2_fda_2026_sweep"
CLEAN = os.path.join(BASE, "data", "clean")
AUD = os.path.join(BASE, "audit")

# Load b2_lemons_results.parquet (ALL 5,118 US FDA-regulated non-P1 submissions in 2024-01-01..2026-10-31,
# covering BOTH On-Time (<=365d) and Overdue (>365d) submissions!)
lem = pd.read_parquet(os.path.join(CLEAN, "b2_lemons_results.parquet"))
import sys
sys.path.insert(0, os.path.join(BASE, "src"))
from analyze_b2_v2_stationary import parse_ctgov_date

lem["pcd"] = parse_ctgov_date(lem["primary_completion_date"], conservative_end=True)
lem["sub_date"] = parse_ctgov_date(lem["results_first_submit_date"], conservative_end=False)
lem["post_date"] = parse_ctgov_date(lem["results_first_post_date"], conservative_end=False)
lem = lem[lem["pcd"].notna() & lem["sub_date"].notna()].copy()
lem["delay_days"] = (lem["sub_date"] - lem["pcd"]).dt.days
lem["post_lag_days"] = (lem["post_date"] - lem["sub_date"]).dt.days
lem["is_overdue"] = (lem["delay_days"] > 365).astype(int)
lem["is_industry"] = (lem["lead_sponsor_class"] == "INDUSTRY").astype(int)
lem["ym"] = lem["sub_date"].dt.to_period("M").astype(str)

# Restrict to 2025-04..2026-09
sub = lem[(lem["ym"] >= "2025-04") & (lem["ym"] <= "2026-09")].copy()
piv = sub.groupby(["ym", "is_overdue", "is_industry"])["nct_id"].count().unstack(fill_value=0)
print("=== 1. MONTHLY SUBMISSION COUNTS: ON-TIME (overdue=0) vs OVERDUE (overdue=1) x INDUSTRY/ACADEMIC ===")
piv_flat = sub.groupby(["ym", "is_overdue", "is_industry"])["nct_id"].count().reset_index()
piv_flat["col"] = piv_flat.apply(lambda r: f"{'OD' if r.is_overdue else 'OT'}_{'IND' if r.is_industry else 'ACA'}", axis=1)
pt = piv_flat.pivot(index="ym", columns="col", values="nct_id").fillna(0).astype(int)
pt["OT_IND_share"] = (pt["OT_IND"] / (pt["OT_IND"] + pt["OT_ACA"]) * 100).round(1)
pt["OD_IND_share"] = (pt["OD_IND"] / (pt["OD_IND"] + pt["OD_ACA"]) * 100).round(1)
print(pt.to_string())

# Compute exact DDD on counts: Overdue vs On-Time x Industry vs Academic x Post (Apr-Jul 2026, Apr-Jun 2026, Apr-Sep 2026)!
# PLUS exact H=120d & H=90d zero-censoring uniform horizon on Overdue vs On-Time!
print("\n=== 2. SMOKING-GUN DDD TEST: OVERDUE vs ON-TIME x INDUSTRY vs ACADEMIC x POST ===")
results = {"monthly_table": pt.to_dict(orient="index"), "tests": {}}

for tag, end_ym, h_cap in [
    ("raw_apr_jul_2026", "2026-07", None),
    ("raw_apr_jun_2026", "2026-06", None),
    ("raw_apr_sep_2026", "2026-09", None),
    ("H120_zero_cens_apr_jul", "2026-07", 120),
    ("H90_zero_cens_apr_aug", "2026-08", 90),
    ("H150_zero_cens_apr_jun", "2026-06", 150),
    ("raw_apr_may_2026", "2026-05", None),
]:
    s = lem[(lem["ym"] >= "2025-04") & (lem["ym"] <= end_ym)].copy()
    if h_cap is not None:
        s = s[s["post_lag_days"].notna() & (s["post_lag_days"] <= h_cap)].copy()
    s["post"] = (s["ym"] >= "2026-04").astype(int)
    # 2x2x2 counts
    c = s.groupby(["is_overdue", "post", "is_industry"])["nct_id"].count().unstack(fill_value=0)
    # Within On-Time (placebo): Industry/Academic ratio Post vs Pre!
    ot_pre_i, ot_pre_a = int(c.loc[(0, 0), 1]), int(c.loc[(0, 0), 0])
    ot_pos_i, ot_pos_a = int(c.loc[(0, 1), 1]), int(c.loc[(0, 1), 0])
    ot_or, ot_p = stats.fisher_exact([[ot_pos_i, ot_pos_a], [ot_pre_i, ot_pre_a]])
    # Within Overdue (treated): Industry/Academic ratio Post vs Pre!
    od_pre_i, od_pre_a = int(c.loc[(1, 0), 1]), int(c.loc[(1, 0), 0])
    od_pos_i, od_pos_a = int(c.loc[(1, 1), 1]), int(c.loc[(1, 1), 0])
    od_or, od_p = stats.fisher_exact([[od_pos_i, od_pos_a], [od_pre_i, od_pre_a]])
    # Exact Logistic DDD: is_industry ~ is_overdue * post (plus exact month FE in conditional/unconditional)
    m_ddd = smf.logit("is_industry ~ is_overdue * post", data=s).fit(disp=False)
    ddd_or = float(np.exp(m_ddd.params["is_overdue:post"]))
    ddd_p = float(m_ddd.pvalues["is_overdue:post"])
    # Also check raw Post/Pre monthly count rate ratio in On-Time vs Overdue separately for Industry & Academic!
    n_pre_m = 12.0
    n_pos_m = float(len(pd.period_range("2026-04", end_ym, freq="M")))
    ot_ind_rr = (ot_pos_i / n_pos_m) / (ot_pre_i / n_pre_m)
    ot_aca_rr = (ot_pos_a / n_pos_m) / (ot_pre_a / n_pre_m)
    od_ind_rr = (od_pos_i / n_pos_m) / (od_pre_i / n_pre_m)
    od_aca_rr = (od_pos_a / n_pos_m) / (od_pre_a / n_pre_m)

    results["tests"][tag] = {
        "on_time_2x2": {"pre_ind": ot_pre_i, "pre_aca": ot_pre_a, "post_ind": ot_pos_i, "post_aca": ot_pos_a, "ind_vs_aca_OR": round(float(ot_or), 3), "p": float(ot_p), "ind_mo_ratio": round(ot_ind_rr, 3), "aca_mo_ratio": round(ot_aca_rr, 3)},
        "overdue_2x2": {"pre_ind": od_pre_i, "pre_aca": od_pre_a, "post_ind": od_pos_i, "post_aca": od_pos_a, "ind_vs_aca_OR": round(float(od_or), 3), "p": float(od_p), "ind_mo_ratio": round(od_ind_rr, 3), "aca_mo_ratio": round(od_aca_rr, 3)},
        "ddd_overdue_vs_ontime_OR": round(ddd_or, 3),
        "ddd_p": ddd_p,
    }
    print(f"  [{tag:24s}] ON-TIME Ind/Aca OR={ot_or:.3f} (p={ot_p:.4f}) [Ind_RR={ot_ind_rr:.2f}, Aca_RR={ot_aca_rr:.2f}] | "
          f"OVERDUE Ind/Aca OR={od_or:.3f} (p={od_p:.4f}) [Ind_RR={od_ind_rr:.2f}, Aca_RR={od_aca_rr:.2f}] | "
          f"DDD(OD vs OT) OR={ddd_or:.3f} (p={ddd_p:.4f})")

# Also check Historical Placebo 2025 on the EXACT SAME On-Time vs Overdue DDD!
s25 = lem[(lem["ym"] >= "2024-04") & (lem["ym"] <= "2025-07")].copy()
s25["post"] = (s25["ym"] >= "2025-04").astype(int)
m25 = smf.logit("is_industry ~ is_overdue * post", data=s25).fit(disp=False)
print(f"  [PLACEBO 2025 Apr-Jul    ] DDD(OD vs OT) OR={np.exp(m25.params['is_overdue:post']):.3f} (p={m25.pvalues['is_overdue:post']:.4f})")
results["placebo_2025_ddd"] = {"or": float(np.exp(m25.params["is_overdue:post"])), "p": float(m25.pvalues["is_overdue:post"])}

# Check QC lag in 2026Q2 specifically for On-Time vs Overdue!
print("\n=== 3. QC LAG IN 2026-04..07: ON-TIME vs OVERDUE x INDUSTRY/ACADEMIC ===")
s_q = lem[(lem["ym"] >= "2026-04") & (lem["ym"] <= "2026-07") & lem["post_lag_days"].notna()]
print(s_q.groupby(["is_overdue", "is_industry"])["post_lag_days"].agg(["count", "median", lambda x: np.percentile(x, 75)]).rename(columns={"<lambda_0>": "p75"}).to_string())

with open(os.path.join(AUD, "b2_ontime_smoking_gun.json"), "w") as f:
    json.dump(results, f, indent=2)
print("\nSaved b2_ontime_smoking_gun.json")
