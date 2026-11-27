#!/opt/anaconda3/bin/python3
"""
Ruthless Forensic Falsification Audit for Study B-2:
  1. ACTUAL vs ESTIMATED PrimaryCompletionDate audit across all 4 arms x Industry/Academic
  2. Symmetric COMPLETED|TERMINATED + ACTUAL-only re-estimation (removing ACTIVE_NOT_RECRUITING arm asymmetry!)
  3. Exact DAILY & WEEKLY submission trajectory around March 30, 2026 (2026-02-01..2026-06-30)
     to verify whether Industry surged *after* March 30 vs pre-existing mid-March / Q1 seasonality!
  4. Calendar-month Q1/Q2 seasonality audit across 2022, 2023, 2024, 2025, 2026!
  5. NLM QC lag dynamics over time (2024-Q1..2026-Q3) + Uniform-Horizon sensitivity across H=60, 90, 120, 150, 180d
  6. Two-way clustering (Sponsor + Month), sponsor-name deduplication audit, and linear-pretrend-adjusted bounds!
"""
import json, os, warnings
import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.api as sm
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")

BASE = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_b2_fda_2026_sweep"
CLEAN = os.path.join(BASE, "data", "clean")
RES = os.path.join(BASE, "results")
AUD_DIR = os.path.join(BASE, "audit")
os.makedirs(AUD_DIR, exist_ok=True)

import sys
sys.path.insert(0, os.path.join(BASE, "src"))
from analyze_b2_v2_stationary import load_steady_universe, build_stationary_monthly_panel, run_within_act_governance_did, run_age_adjusted_ddd

df = load_steady_universe()

out = {}

# 1. Audit primary_completion_type (ACTUAL vs ESTIMATED) & overall_status across all 4 arms x Industry
print("=== 1. ACTUAL vs ESTIMATED PCD & STATUS AUDIT ===")
tab1 = df.groupby(["arm", "is_industry", "primary_completion_type"], dropna=False)["nct_id"].count().unstack(fill_value=0)
print(tab1.to_string())
tab_st = df.groupby(["arm", "overall_status"], dropna=False)["nct_id"].count().unstack(fill_value=0)
print("\nOverall status by arm:")
print(tab_st.to_string())
out["pcd_type_counts"] = {str(k): {str(c): int(v) for c, v in row.items()} for k, row in tab1.to_dict(orient="index").items()}
out["status_counts"] = {str(k): {str(c): int(v) for c, v in row.items()} for k, row in tab_st.to_dict(orient="index").items()}

# 2. Exact WEEKLY & DAILY trajectory around March 30, 2026 (strictly COMPLETED|TERMINATED + ACTUAL PCD!)
df_strict = df[(df["overall_status"].isin(["COMPLETED", "TERMINATED"])) & (df["primary_completion_type"] == "ACTUAL")].copy()
print(f"\nStrict symmetric universe (COMPLETED|TERMINATED & ACTUAL PCD): {len(df_strict):,} / {len(df):,} ({len(df_strict)/len(df)*100:.1f}%)")
out["strict_n_trials"] = int(len(df_strict))

act_s = df_strict[df_strict["arm"] == "us_fdareg_nonp1"].copy()
# Exact Monday-anchored weekly counts & hazards from week -12 (2026-01-05) to week +17 (2026-07-27)
t0 = pd.Timestamp("2026-03-30")
weekly_rows = []
for k in range(-16, 18):
    ws = t0 + pd.Timedelta(days=7 * k)
    we = ws + pd.Timedelta(days=6, hours=23, minutes=59, seconds=59)
    ar = (act_s["due_date"] < ws) & (act_s["sub_date"].isna() | (act_s["sub_date"] >= ws))
    sub = act_s.loc[ar]
    ev = (sub["sub_date"] >= ws) & (sub["sub_date"] <= we)
    row = {"rel_week": k, "week_start": ws.strftime("%Y-%m-%d")}
    for iv, nm in [(1, "ind"), (0, "aca")]:
        m = sub["is_industry"] == iv
        rk = int(m.sum())
        e = int((ev & m).sum())
        row[f"{nm}_risk"] = rk
        row[f"{nm}_events"] = e
        row[f"{nm}_haz_pct"] = round(e / rk * 100.0, 3) if rk > 0 else 0.0
    weekly_rows.append(row)

print("\n=== 2. EXACT MONDAY-ANCHORED WEEKLY OVERDUE SUBMISSIONS AROUND 2026-03-30 (STRICT ACTUAL) ===")
wdf = pd.DataFrame(weekly_rows)
print(wdf.to_string(index=False))
out["weekly_exact_around_mar30"] = weekly_rows

# Exact day-by-day overdue submissions in March 15 – April 20, 2026 (check Mar 1..29 vs Mar 30..Apr 20!)
mar_apr = act_s[(act_s["sub_date"] >= "2026-03-01") & (act_s["sub_date"] <= "2026-04-30") & (act_s["sub_date"] > act_s["due_date"])].copy()
mar_apr["day"] = mar_apr["sub_date"].dt.strftime("%Y-%m-%d")
mar_pre = mar_apr[mar_apr["sub_date"] < "2026-03-30"]
mar_30_31 = mar_apr[(mar_apr["sub_date"] >= "2026-03-30") & (mar_apr["sub_date"] <= "2026-03-31")]
apr_all = mar_apr[mar_apr["sub_date"] >= "2026-04-01"]
print("\n=== 3. MARCH 2026 INTRA-MONTH DECOMPOSITION (Mar 1–29 vs Mar 30–31 vs Apr 1–30) ===")
for lbl, sub_d, n_days in [("Mar 1-29 (Pre-sweep, 29d)", mar_pre, 29), ("Mar 30-31 (Sweep days, 2d)", mar_30_31, 2), ("Apr 1-30 (Post-sweep, 30d)", apr_all, 30)]:
    ni = int((sub_d["is_industry"] == 1).sum())
    na = int((sub_d["is_industry"] == 0).sum())
    print(f"  {lbl:30s}: Industry={ni:2d} ({ni/n_days:.2f}/day) | Academic={na:2d} ({na/n_days:.2f}/day)")
out["march_intramonth_split"] = {
    "mar_1_29_29d": {"ind": int((mar_pre["is_industry"] == 1).sum()), "aca": int((mar_pre["is_industry"] == 0).sum())},
    "mar_30_31_2d": {"ind": int((mar_30_31["is_industry"] == 1).sum()), "aca": int((mar_30_31["is_industry"] == 0).sum())},
    "apr_1_30_30d": {"ind": int((apr_all["is_industry"] == 1).sum()), "aca": int((apr_all["is_industry"] == 0).sum())},
}

# 4. Multi-year Q1/Q2 Monthly Seasonality (2023, 2024, 2025, 2026 Jan-Jul) in Strict Actual ACT_FDA!
print("\n=== 4. MULTI-YEAR MONTHLY SEASONALITY CHECK (ACT_FDA Industry vs Academic, Jan–Jul across 2023–2026) ===")
seas = []
for yr in (2023, 2024, 2025, 2026):
    for mo in range(1, 8):
        ms = pd.Timestamp(f"{yr}-{mo:02d}-01")
        me = ms + pd.offsets.MonthEnd(0)
        ar = (act_s["due_date"] < ms) & (act_s["sub_date"].isna() | (act_s["sub_date"] >= ms))
        sub = act_s.loc[ar]
        ev = (sub["sub_date"] >= ms) & (sub["sub_date"] <= me)
        hi = float((ev & (sub["is_industry"] == 1)).sum() / max((sub["is_industry"] == 1).sum(), 1) * 100)
        ha = float((ev & (sub["is_industry"] == 0)).sum() / max((sub["is_industry"] == 0).sum(), 1) * 100)
        seas.append({"year": yr, "month": mo, "ind_h": round(hi, 2), "aca_h": round(ha, 2), "ratio_i_a": round(hi / max(ha, 0.01), 2)})
sdf = pd.DataFrame(seas)
print("Industry hazard (%) by year x month:")
print(sdf.pivot(index="year", columns="month", values="ind_h").to_string())
print("Academic hazard (%) by year x month:")
print(sdf.pivot(index="year", columns="month", values="aca_h").to_string())
print("Industry/Academic Ratio by year x month:")
print(sdf.pivot(index="year", columns="month", values="ratio_i_a").to_string())
out["multiyear_seasonality"] = seas

# 5. Re-run Primary Models STRICTLY on Symmetric COMPLETED|TERMINATED + ACTUAL PCD universe!
print("\n=== 5. STRICT SYMMETRIC (COMPLETED|TERMINATED + ACTUAL PCD ONLY) MODELS ===")
p_strict_jul = build_stationary_monthly_panel(df_strict, "2025-04", "2026-07", "2026-04")
p_strict_jun = build_stationary_monthly_panel(df_strict, "2025-04", "2026-06", "2026-04")
p_strict_p25 = build_stationary_monthly_panel(df_strict, "2024-04", "2025-07", "2025-04")
p_strict_p24 = build_stationary_monthly_panel(df_strict, "2023-04", "2024-07", "2024-04")

out["strict_models"] = {
    "within_act_jul": run_within_act_governance_did(p_strict_jul, "Strict Symmetric Apr-Jul 2026"),
    "within_act_jun": run_within_act_governance_did(p_strict_jun, "Strict Symmetric Apr-Jun 2026"),
    "within_act_p25": run_within_act_governance_did(p_strict_p25, "Strict Symmetric Placebo 2025"),
    "within_act_p24": run_within_act_governance_did(p_strict_p24, "Strict Symmetric Placebo 2024"),
    "ddd_clean_exempt_jul": run_age_adjusted_ddd(p_strict_jul, ["us_notfdareg", "nonus_p23_nonfda"], "Strict Symmetric DDD vs Clean Exempt Apr-Jul 2026"),
    "ddd_us_notfdareg_jul": run_age_adjusted_ddd(p_strict_jul, ["us_notfdareg"], "Strict Symmetric DDD vs US Non-FDA Apr-Jul 2026"),
    "ddd_nonus_p23_jul": run_age_adjusted_ddd(p_strict_jul, ["nonus_p23_nonfda"], "Strict Symmetric DDD vs Non-US P23 Apr-Jul 2026"),
}
for k, v in out["strict_models"].items():
    if "poisson_age_fe" in v:
        print(f"  [{k:24s}] LPM={v['lpm_age_fe']['coef_pp']:+.3f}pp (p={v['lpm_age_fe']['p']:.3e}) | Pois_IRR={v['poisson_age_fe']['irr']:.3f} {v['poisson_age_fe']['irr_ci95']} (p={v['poisson_age_fe']['p']:.3e})")
    else:
        print(f"  [{k:24s}] LPM_DDD={v['lpm_ddd_age_fe']['ddd']['coef_pp']:+.3f}pp (p={v['lpm_ddd_age_fe']['ddd']['p']:.3e}) | Pois_DDD={v['poisson_ddd_age_fe']['ddd']['irr']:.3f} (p={v['poisson_ddd_age_fe']['ddd']['p']:.3e}) | Ind={v['stratified_did']['industry']['irr']:.3f} (p={v['stratified_did']['industry']['p']:.3e}) | Aca={v['stratified_did']['academic']['irr']:.3f} (p={v['stratified_did']['academic']['p']:.3e})")

# 6. NLM QC Lag Drift Over Time (2024-Q1..2026-Q3) + Exact Uniform-Horizon Grid H=60, 90, 120, 150, 180d on STRICT ACTUAL!
print("\n=== 6. NLM QC LAG DRIFT OVER TIME (2024-Q1..2026-Q3) ===")
s_lag = act_s[act_s["sub_date"].notna() & (act_s["sub_date"] >= "2024-01-01") & (act_s["sub_date"] <= "2026-09-30") & act_s["post_lag_days"].notna()].copy()
s_lag["qtr"] = s_lag["sub_date"].dt.to_period("Q").astype(str)
q_lag = s_lag.groupby(["qtr", "is_industry"])["post_lag_days"].agg(["count", "median", lambda x: np.percentile(x, 75)]).rename(columns={"<lambda_0>": "p75"})
print(q_lag.to_string())
out["qc_lag_by_quarter"] = {f"{k[0]}_I{k[1]}": {"n": int(r["count"]), "median": float(r["median"]), "p75": float(r["p75"])} for k, r in q_lag.iterrows()}

# Uniform-Horizon Grid on STRICT ACTUAL (using exact end_ym where end_ym + H <= 2026-11-25)
print("\n=== 7. UNIFORM-HORIZON ZERO-CENSORING GRID ON STRICT ACTUAL ===")
uh_grid = {}
for H, end_ym in [(60, "2026-09"), (90, "2026-08"), (120, "2026-07"), (150, "2026-06"), (180, "2026-05")]:
    months = pd.period_range("2025-04", end_ym, freq="M")
    recs = []
    for m in months:
        ms, me = m.start_time, m.end_time
        ar = (df_strict["due_date"] < ms) & (df_strict["sub_date"].isna() | (df_strict["sub_date"] >= ms))
        sub = df_strict.loc[ar, ["nct_id", "arm", "is_act", "is_industry", "lead_sponsor_clean", "due_date", "sub_date", "post_lag_days"]].copy()
        sub["ym"] = str(m)
        sub["post"] = int(str(m) >= "2026-04")
        sub["y"] = ((sub["sub_date"] >= ms) & (sub["sub_date"] <= me) & sub["post_lag_days"].notna() & (sub["post_lag_days"] <= H)).astype(int)
        od_m = ((ms - sub["due_date"]).dt.days / 30.4375).clip(lower=0.0)
        sub["age_bin"] = pd.cut(od_m, bins=[-0.01, 3.0, 6.0, 12.0, 24.0, 48.0, 999.0],
                                labels=["0_3m", "3_6m", "6_12m", "12_24m", "24_48m", "48m_plus"]).astype(str)
        recs.append(sub)
    pan = pd.concat(recs, ignore_index=True)
    act_p = pan[pan["arm"] == "us_fdareg_nonp1"].copy()
    cdf = act_p.groupby(["lead_sponsor_clean", "is_industry", "age_bin", "ym", "post"], as_index=False).agg(events=("y", "sum"), risk=("y", "count"))
    cdf["ind_post"] = cdf["is_industry"] * cdf["post"]
    cdf["cell_ia"] = cdf["is_industry"].astype(str) + "_" + cdf["age_bin"]
    mp = smf.glm("events ~ ind_post + C(cell_ia) + C(ym)", data=cdf, family=sm.families.Poisson(), offset=np.log(cdf["risk"].values)).fit(
        cov_type="cluster", cov_kwds={"groups": cdf["lead_sponsor_clean"]}
    )
    b, se = float(mp.params["ind_post"]), float(mp.bse["ind_post"])
    # Also Clean-Exempt DDD
    ce = pan[pan["arm"].isin(["us_fdareg_nonp1", "us_notfdareg", "nonus_p23_nonfda"])].copy()
    cd2 = ce.groupby(["lead_sponsor_clean", "arm", "is_act", "is_industry", "age_bin", "ym", "post"], as_index=False).agg(events=("y", "sum"), risk=("y", "count"))
    cd2["ddd"] = cd2["is_act"] * cd2["is_industry"] * cd2["post"]
    cd2["act_post"] = cd2["is_act"] * cd2["post"]
    cd2["ind_post"] = cd2["is_industry"] * cd2["post"]
    cd2["cell_aia"] = cd2["arm"].astype(str) + "_" + cd2["is_industry"].astype(str) + "_" + cd2["age_bin"]
    cd2 = cd2[cd2.groupby("cell_aia")["events"].transform("sum") > 0]
    mp2 = smf.glm("events ~ ddd + act_post + ind_post + C(cell_aia) + C(ym)", data=cd2, family=sm.families.Poisson(), offset=np.log(cd2["risk"].values)).fit(
        cov_type="cluster", cov_kwds={"groups": cd2["lead_sponsor_clean"]}
    )
    uh_grid[f"H{H}_{end_ym}"] = {
        "n_act_events": int(act_p["y"].sum()),
        "within_act_irr": float(np.exp(b)),
        "within_act_ci95": [float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se))],
        "within_act_p": float(mp.pvalues["ind_post"]),
        "ddd_irr": float(np.exp(mp2.params["ddd"])),
        "ddd_p": float(mp2.pvalues["ddd"]),
    }
    print(f"  H={H:3d}d (thru {end_ym}): events={act_p['y'].sum():3d} | Within-ACT IRR={np.exp(b):.3f} [{np.exp(b-1.96*se):.2f},{np.exp(b+1.96*se):.2f}] (p={mp.pvalues['ind_post']:.4f}) | DDD IRR={np.exp(mp2.params['ddd']):.3f} (p={mp2.pvalues['ddd']:.4f})")
out["uniform_horizon_grid_strict"] = uh_grid

# 8. Two-Way Clustering (Sponsor + Month) on Strict Within-ACT LPM & Poisson (Cameron-Gelbach-Miller formula)
print("\n=== 8. TWO-WAY CLUSTERING (SPONSOR + CALENDAR MONTH) CHECK ===")
act_j = p_strict_jul[p_strict_jul["arm"] == "us_fdareg_nonp1"].copy()
act_j["ind_post"] = act_j["is_industry"] * act_j["post"]
act_j["y_bp"] = act_j["y"] * 100.0
act_j["cell_ia"] = act_j["is_industry"].astype(str) + "_" + act_j["age_bin"]
act_j["sm_inter"] = act_j["lead_sponsor_clean"] + "_" + act_j["ym"]
m_s = smf.ols("y_bp ~ ind_post + C(cell_ia) + C(ym)", data=act_j).fit(cov_type="cluster", cov_kwds={"groups": act_j["lead_sponsor_clean"]})
m_m = smf.ols("y_bp ~ ind_post + C(cell_ia) + C(ym)", data=act_j).fit(cov_type="cluster", cov_kwds={"groups": act_j["ym"]})
m_sm = smf.ols("y_bp ~ ind_post + C(cell_ia) + C(ym)", data=act_j).fit(cov_type="cluster", cov_kwds={"groups": act_j["sm_inter"]})
var_2w = float(m_s.bse["ind_post"]**2 + m_m.bse["ind_post"]**2 - m_sm.bse["ind_post"]**2)
se_2w = float(np.sqrt(max(var_2w, m_s.bse["ind_post"]**2)))
b_2w = float(m_s.params["ind_post"])
p_2w = float(2 * stats.norm.sf(abs(b_2w / se_2w)))
print(f"  Strict LPM ind_post={b_2w:+.3f}pp | SE_sponsor={m_s.bse['ind_post']:.3f} | SE_month={m_m.bse['ind_post']:.3f} | SE_2way={se_2w:.3f} (p_2way={p_2w:.4f})")
out["twoway_cluster_strict"] = {"coef_pp": b_2w, "se_sponsor": float(m_s.bse["ind_post"]), "se_month": float(m_m.bse["ind_post"]), "se_2way": se_2w, "p_2way": p_2w}

# Also check exact Post-Sweep Start Sensitivity: March 30, 2026 weekly exact DiD (26 pre weeks vs 17 post weeks through Jul 26!)
print("\n=== 9. EXACT WEEKLY (DAY-0 = 2026-03-30) POISSON DiD ON STRICT ACTUAL ===")
w_recs = []
for k in range(-52, 18):  # 52 pre weeks (2025-03-31..2026-03-29) vs 18 post weeks (2026-03-30..2026-08-02)
    ws = t0 + pd.Timedelta(days=7 * k)
    we = ws + pd.Timedelta(days=6, hours=23, minutes=59, seconds=59)
    ar = (act_s["due_date"] < ws) & (act_s["sub_date"].isna() | (act_s["sub_date"] >= ws))
    sub = act_s.loc[ar, ["is_industry", "lead_sponsor_clean", "due_date", "sub_date"]].copy()
    sub["wk"] = k
    sub["post"] = int(k >= 0)
    sub["y"] = ((sub["sub_date"] >= ws) & (sub["sub_date"] <= we)).astype(int)
    od_m = ((ws - sub["due_date"]).dt.days / 30.4375).clip(lower=0.0)
    sub["age_bin"] = pd.cut(od_m, bins=[-0.01, 3.0, 6.0, 12.0, 24.0, 48.0, 999.0],
                            labels=["0_3m", "3_6m", "6_12m", "12_24m", "24_48m", "48m_plus"]).astype(str)
    w_recs.append(sub)
wp = pd.concat(w_recs, ignore_index=True)
wcdf = wp.groupby(["lead_sponsor_clean", "is_industry", "age_bin", "wk", "post"], as_index=False).agg(events=("y", "sum"), risk=("y", "count"))
wcdf["ind_post"] = wcdf["is_industry"] * wcdf["post"]
wcdf["cell_ia"] = wcdf["is_industry"].astype(str) + "_" + wcdf["age_bin"]
mw = smf.glm("events ~ ind_post + C(cell_ia) + C(wk)", data=wcdf, family=sm.families.Poisson(), offset=np.log(wcdf["risk"].values)).fit(
    cov_type="cluster", cov_kwds={"groups": wcdf["lead_sponsor_clean"]}
)
bw, sew = float(mw.params["ind_post"]), float(mw.bse["ind_post"])
print(f"  Exact Weekly (k=-52..+17, split at 2026-03-30): IRR={np.exp(bw):.3f} [{np.exp(bw-1.96*sew):.2f},{np.exp(bw+1.96*sew):.2f}] (p={mw.pvalues['ind_post']:.4e})")
out["weekly_exact_poisson_strict"] = {"irr": float(np.exp(bw)), "ci95": [float(np.exp(bw-1.96*sew)), float(np.exp(bw+1.96*sew))], "p": float(mw.pvalues["ind_post"])}

# Also check if excluding Nov 2025 (outlier) OR restricting pre-period to 2025-04..2026-03 vs 2024-04..2026-03 changes anything!
p_2yr = build_stationary_monthly_panel(df_strict, "2024-04", "2026-07", "2026-04")
res_2yr = run_within_act_governance_did(p_2yr, "24-month pre-period (2024-04..2026-07)")
p_nonov = build_stationary_monthly_panel(df_strict, "2025-04", "2026-07", "2026-04", drop_ym=["2025-11"])
res_nonov = run_within_act_governance_did(p_nonov, "Drop Nov 2025 outlier")
print(f"  24-month pre (2024-04..2026-07): Pois IRR={res_2yr['poisson_age_fe']['irr']:.3f} (p={res_2yr['poisson_age_fe']['p']:.4e}) | LPM={res_2yr['lpm_age_fe']['coef_pp']:+.3f}pp (p={res_2yr['lpm_age_fe']['p']:.4e})")
print(f"  Drop Nov 2025 outlier:           Pois IRR={res_nonov['poisson_age_fe']['irr']:.3f} (p={res_nonov['poisson_age_fe']['p']:.4e}) | LPM={res_nonov['lpm_age_fe']['coef_pp']:+.3f}pp (p={res_nonov['lpm_age_fe']['p']:.4e})")
out["sensitivity_2yr_pre"] = res_2yr
out["sensitivity_drop_nov25"] = res_nonov

with open(os.path.join(AUD_DIR, "b2_forensic_audit.json"), "w") as f:
    json.dump(out, f, indent=2)
print("\nSaved b2_forensic_audit.json")
