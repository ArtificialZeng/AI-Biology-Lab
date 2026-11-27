#!/opt/anaconda3/bin/python3
"""
Study B-2 — Attack 2D & 2E: Falsifying the Surviving Acute March 30–April 12 Discontinuity
==========================================================================================
CLASSIFICATION: EXPLICITLY POST-HOC EXPLORATORY FALSIFICATION (preserving all negative/null results!).

Tests:
  2D-1. Historical Q1 Quarter-End Placebos (2024-03-30..2024-04-12 and 2025-03-30..2025-04-12 vs 2026):
        Does ClinicalTrials.gov exhibit a generic Q1 quarter-end (Mar 30–31) submission spike in 2024/2025?
  2D-2. Multi-Quarter-End Bunching Diagnostic (Q1/Q2/Q3/Q4 final 2 business days in 2024–2025 vs non-QE days).
  2D-3. Q1-Donut Pure Private-Email Window (2026-04-01..2026-04-12):
        Drops Mar 30–31 (Q1-end) completely AND stops before Apr 13 (public FDA press release).
        Tests both plain Period FE and Sponsor x Period FE (explicitly noting Sponsor x Period FE absorbs only
        common sponsor-period shocks, NOT status-differential shocks, selection, or measurement error!).
  2D-4. Year-over-Year Same-Calendar-Window Comparison (2026 vs 2024–2025 pooled Mar 30–Apr 12 & Apr 1–12).
  2E.   Reporting-Proxy / Measurement-Error Sensitivity:
        Repeats acute window tests on posted-only `results_first_submit` vs `sub_combined_proxy`.
"""
import json, os, sys, warnings
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")

BASE = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_b2_fda_2026_sweep"
AUD = os.path.join(BASE, "audit")
OUT_DIR = os.path.join(BASE, "falsification_confounding_suite")
sys.path.insert(0, AUD)
sys.path.insert(0, OUT_DIR)

from audit_b2_t00_t06_t07_closure import load_clean_t00_universe
from run_b2_attack2_confounding_falsification import fit_ppml_fe

df = load_clean_t00_universe()
act = df[df["is_act"] == 1].copy()
print("Loaded ACTs:", len(act), "columns matching submit/date:", [c for c in act.columns if "sub" in c or "date" in c or "first" in c])

# Parse posted-only submit date columns cleanly as datetime
act["posted_first_submit_dt"] = pd.to_datetime(act["results_first_submit_date"], errors="coerce")
posted_col = "posted_first_submit_dt"
print("Posted-only submit col:", posted_col, "non-null:", int(act[posted_col].notna().sum()), "| sub_posted_qc non-null:", int(act["sub_posted_qc"].notna().sum()))

def daily_rates_for_range(act_df, sub_col, start_str, end_str):
    rows = []
    for d_ts in pd.date_range(start_str, end_str, freq="D"):
        ar = (act_df["pcd"] < d_ts) & (act_df[sub_col].isna() | (act_df[sub_col] >= d_ts))
        sub = act_df.loc[ar, ["is_industry", "due_date", sub_col]]
        is_od = (sub["due_date"] < d_ts)
        y_d = (sub[sub_col] == d_ts)
        for od in (0, 1):
            for ind in (0, 1):
                m = (is_od == (od == 1)) & (sub["is_industry"] == ind)
                rows.append({
                    "date": d_ts.strftime("%Y-%m-%d"),
                    "year": d_ts.year,
                    "md": d_ts.strftime("%m-%d"),
                    "dow": d_ts.day_name()[:3],
                    "is_overdue": od,
                    "is_industry": ind,
                    "events": int((y_d & m).sum()),
                    "risk": int(m.sum()),
                })
    return pd.DataFrame(rows)

def summarize_window(ddf, s_md, e_md, yr, lbl):
    g = ddf[(ddf["year"] == yr) & (ddf["md"] >= s_md) & (ddf["md"] <= e_md)]
    out = {"year": yr, "window": lbl, "start_md": s_md, "end_md": e_md, "days": int(g["date"].nunique())}
    for od, ind, k in [(1,1,"ind_od"), (0,1,"ind_ot"), (1,0,"acad_od"), (0,0,"acad_ot")]:
        c = g[(g["is_overdue"] == od) & (g["is_industry"] == ind)]
        ev = int(c["events"].sum())
        rk = int(c["risk"].sum())
        r30 = 100.0 * (ev / rk) * 30.4375 if rk > 0 else 0.0
        out[f"{k}_events"] = ev
        out[f"{k}_rate30_pct"] = round(r30, 4)
    out["ind_od_ot_ratio"] = round(out["ind_od_rate30_pct"] / out["ind_ot_rate30_pct"], 4) if out["ind_ot_rate30_pct"] > 0 else None
    out["acad_od_ot_ratio"] = round(out["acad_od_rate30_pct"] / out["acad_ot_rate30_pct"], 4) if out["acad_ot_rate30_pct"] > 0 else None
    return out

# 1. Historical Q1 Quarter-End Placebo across 2024, 2025, 2026!
sub_col = "sub_combined_proxy"
d_all = pd.concat([
    daily_rates_for_range(act, sub_col, "2024-02-01", "2024-05-15"),
    daily_rates_for_range(act, sub_col, "2025-02-01", "2025-05-15"),
    daily_rates_for_range(act, sub_col, "2026-02-01", "2026-05-15"),
], ignore_index=True)

q1_placebo_rows = []
for yr in [2024, 2025, 2026]:
    for lbl, smd, emd in [
        ("Feb 1-28", "02-01", "02-28"),
        ("Mar 1-29 (Pre-Mar30)", "03-01", "03-29"),
        ("Mar 30-31 (Q1 Final 2d)", "03-30", "03-31"),
        ("Apr 1-12 (Q1-Donut Pure Window)", "04-01", "04-12"),
        ("Mar 30-Apr 12 (Full 14d Window)", "03-30", "04-12"),
        ("Apr 13-30 (Post-Apr13 Window)", "04-13", "04-30"),
    ]:
        q1_placebo_rows.append(summarize_window(d_all, smd, emd, yr, lbl))

# 2. Other Quarter-Ends in 2024-2025 (Q2: Jun 29-30; Q3: Sep 29-30; Q4: Dec 30-31) vs non-QE month baseline
other_qe_rows = []
for lbl, s_dt, e_dt, base_s, base_e in [
    ("2024 Q2-end (Jun 28-30)", "2024-06-28", "2024-06-30", "2024-06-01", "2024-06-27"),
    ("2024 Q3-end (Sep 29-30)", "2024-09-29", "2024-09-30", "2024-09-01", "2024-09-28"),
    ("2024 Q4-end (Dec 30-31)", "2024-12-30", "2024-12-31", "2024-12-01", "2024-12-29"),
    ("2025 Q2-end (Jun 29-30)", "2025-06-29", "2025-06-30", "2025-06-01", "2025-06-28"),
    ("2025 Q3-end (Sep 29-30)", "2025-09-29", "2025-09-30", "2025-09-01", "2025-09-28"),
    ("2025 Q4-end (Dec 30-31)", "2025-12-30", "2025-12-31", "2025-12-01", "2025-12-29"),
]:
    dqe = daily_rates_for_range(act, sub_col, base_s, e_dt)
    g_qe = dqe[dqe["date"] >= s_dt]
    g_bs = dqe[dqe["date"] <= base_e]
    rec = {"quarter_end": lbl}
    for od, ind, k in [(1,1,"ind_od"), (0,1,"ind_ot"), (1,0,"acad_od")]:
        c_qe = g_qe[(g_qe["is_overdue"]==od) & (g_qe["is_industry"]==ind)]
        c_bs = g_bs[(g_bs["is_overdue"]==od) & (g_bs["is_industry"]==ind)]
        r_qe = 100.0 * (c_qe["events"].sum() / c_qe["risk"].sum()) * 30.4375
        r_bs = 100.0 * (c_bs["events"].sum() / c_bs["risk"].sum()) * 30.4375
        rec[f"{k}_qe_ev"] = int(c_qe["events"].sum())
        rec[f"{k}_qe_rate30"] = round(r_qe, 4)
        rec[f"{k}_base_rate30"] = round(r_bs, 4)
        rec[f"{k}_qe_to_base_ratio"] = round(r_qe / r_bs, 3) if r_bs > 0 else None
    other_qe_rows.append(rec)

# 3. Formal Q1-Donut Regression (splitting 2026-03-30..31 vs 2026-04-01..12 vs 2026-04-13+)
#    Run for BOTH `sub_combined_proxy` AND `posted_only`!
def run_q1_donut_reg(act_df, scol, label):
    regime_periods = []
    for m in pd.period_range("2025-04", "2026-02", freq="M"):
        regime_periods.append((str(m), m.start_time, m.end_time, "PRE_BASE"))
    regime_periods.append(("2026-03_PRE", pd.Timestamp("2026-03-01"), pd.Timestamp("2026-03-29"), "PRE_MAR"))
    regime_periods.append(("2026-0330_Q1END", pd.Timestamp("2026-03-30"), pd.Timestamp("2026-03-31"), "MAR30_31_Q1END"))
    regime_periods.append(("2026-0401_PUREPRIV", pd.Timestamp("2026-04-01"), pd.Timestamp("2026-04-12"), "APR01_12_PUREPRIV"))
    regime_periods.append(("2026-0413_PUB1", pd.Timestamp("2026-04-13"), pd.Timestamp("2026-05-12"), "PUBLIC_PRESS"))
    regime_periods.append(("2026-0513_PUB2", pd.Timestamp("2026-05-13"), pd.Timestamp("2026-06-30"), "PUBLIC_PRESS"))
    for m in pd.period_range("2026-07", "2026-09", freq="M"):
        regime_periods.append((str(m), m.start_time, m.end_time, "PUBLIC_PRESS"))

    recs = []
    for pname, ps, pe, regime in regime_periods:
        ndays = (pe - ps).days + 1
        ar = (act_df["pcd"] < ps) & (act_df[scol].isna() | (act_df[scol] >= ps))
        sub = act_df.loc[ar, ["lead_sponsor_clean", "is_industry", "pcd", "due_date", scol]].copy()
        sub["period"] = pname
        sub["regime"] = regime
        sub["is_overdue"] = (sub["due_date"] < ps).astype(int)
        sub["y"] = ((sub[scol] >= ps) & (sub[scol] <= pe)).astype(int)
        sub["exposure"] = ndays / 30.4375
        age_m = ((ps - sub["pcd"]).dt.days / 30.4375).clip(lower=0.0)
        sub["pcd_age_bin"] = pd.cut(
            age_m, bins=[-0.01, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 24.0, 36.0, 60.0, 999.0],
            labels=["0_3", "3_6", "6_9", "9_12", "12_15", "15_18", "18_24", "24_36", "36_60", "60p"]
        ).astype(str)
        recs.append(sub)
    rpan = pd.concat(recs, ignore_index=True)
    rpan["sponsor_period"] = rpan["lead_sponsor_clean"] + "__" + rpan["period"]

    rcdf = rpan.groupby(["lead_sponsor_clean", "is_overdue", "is_industry", "pcd_age_bin", "period", "regime", "sponsor_period"], as_index=False).agg(
        events=("y", "sum"), risk_exp=("exposure", "sum")
    )
    rcdf["cell"] = rcdf["is_overdue"].astype(str) + "_" + rcdf["is_industry"].astype(str) + "_" + rcdf["pcd_age_bin"]
    rcdf = rcdf[rcdf.groupby("cell")["events"].transform("sum") > 0].copy().reset_index(drop=True)
    rgs = ["PRE_MAR", "MAR30_31_Q1END", "APR01_12_PUREPRIV", "PUBLIC_PRESS"]
    for rg in rgs:
        rcdf[f"ddd_{rg}"] = ((rcdf["is_overdue"] == 1) & (rcdf["is_industry"] == 1) & (rcdf["regime"] == rg)).astype(float)
        rcdf[f"od_{rg}"] = ((rcdf["is_overdue"] == 1) & (rcdf["regime"] == rg)).astype(float)
        rcdf[f"in_{rg}"] = ((rcdf["is_industry"] == 1) & (rcdf["regime"] == rg)).astype(float)

    # 1. Plain Period FE (DDD & Ind OD/OT & Acad OD/OT)
    f_ddd = smf.glm(
        "events ~ " + " + ".join([f"ddd_{rg} + od_{rg} + in_{rg}" for rg in rgs]) + " + C(cell) + C(period)",
        data=rcdf, family=sm.families.Poisson(), offset=np.log(rcdf["risk_exp"].values)
    ).fit(cov_type="cluster", cov_kwds={"groups": rcdf["lead_sponsor_clean"]})

    rci = rcdf[rcdf["is_industry"] == 1].copy()
    rci = rci[rci.groupby("cell")["events"].transform("sum") > 0].copy().reset_index(drop=True)
    f_ind = smf.glm(
        "events ~ " + " + ".join([f"od_{rg}" for rg in rgs]) + " + C(cell) + C(period)",
        data=rci, family=sm.families.Poisson(), offset=np.log(rci["risk_exp"].values)
    ).fit(cov_type="cluster", cov_kwds={"groups": rci["lead_sponsor_clean"]})

    rca = rcdf[rcdf["is_industry"] == 0].copy()
    rca = rca[rca.groupby("cell")["events"].transform("sum") > 0].copy().reset_index(drop=True)
    f_acad = smf.glm(
        "events ~ " + " + ".join([f"od_{rg}" for rg in rgs]) + " + C(cell) + C(period)",
        data=rca, family=sm.families.Poisson(), offset=np.log(rca["risk_exp"].values)
    ).fit(cov_type="cluster", cov_kwds={"groups": rca["lead_sponsor_clean"]})

    # 2. Sponsor x Period FE (Ind OD/OT and Acad OD/OT)
    cdi = pd.get_dummies(rci["cell"], prefix="ci", drop_first=True, dtype=float)
    r_ind_sp = fit_ppml_fe(pd.concat([rci, cdi], axis=1), "events", "risk_exp", [f"od_{rg}" for rg in rgs] + list(cdi.columns), "sponsor_period", "lead_sponsor_clean")
    cda = pd.get_dummies(rca["cell"], prefix="ca", drop_first=True, dtype=float)
    r_acad_sp = fit_ppml_fe(pd.concat([rca, cda], axis=1), "events", "risk_exp", [f"od_{rg}" for rg in rgs] + list(cda.columns), "sponsor_period", "lead_sponsor_clean")

    # 3. Also run HISTORICAL PLACEBO regression on 2024-04..2025-09 (treating 2025-03-30..31 and 2025-04-01..12 as fake sweep!)
    res = {"label": label, "sub_col": scol, "period_fe": {}, "sponsor_period_fe": {}}
    for rg in rgs:
        def pack(fit, k):
            b, se, p = float(fit.params[k]), float(fit.bse[k]), float(fit.pvalues[k])
            return {"irr": round(np.exp(b), 4), "ci95": [round(np.exp(b-1.96*se), 4), round(np.exp(b+1.96*se), 4)], "p": p}
        res["period_fe"][rg] = {
            "ot_od_ddd": pack(f_ddd, f"ddd_{rg}"),
            "ind_od_ot": pack(f_ind, f"od_{rg}"),
            "acad_od_ot": pack(f_acad, f"od_{rg}"),
        }
        res["sponsor_period_fe"][rg] = {
            "ind_od_ot": {"irr": round(r_ind_sp["irr"][f"od_{rg}"], 4), "ci95": [round(x, 4) for x in r_ind_sp["ci95"][f"od_{rg}"]], "p": r_ind_sp["p"][f"od_{rg}"]},
            "acad_od_ot": {"irr": round(r_acad_sp["irr"][f"od_{rg}"], 4), "ci95": [round(x, 4) for x in r_acad_sp["ci95"][f"od_{rg}"]], "p": r_acad_sp["p"][f"od_{rg}"]},
        }
    return res

def run_historical_2025_placebo(act_df, scol):
    regime_periods = []
    for m in pd.period_range("2024-04", "2025-02", freq="M"):
        regime_periods.append((str(m), m.start_time, m.end_time, "PRE_BASE"))
    regime_periods.append(("2025-03_PRE", pd.Timestamp("2025-03-01"), pd.Timestamp("2025-03-29"), "PRE_MAR"))
    regime_periods.append(("2025-0330_Q1END", pd.Timestamp("2025-03-30"), pd.Timestamp("2025-03-31"), "MAR30_31_Q1END"))
    regime_periods.append(("2025-0401_PURE", pd.Timestamp("2025-04-01"), pd.Timestamp("2025-04-12"), "APR01_12_PURE"))
    regime_periods.append(("2025-0413_POST", pd.Timestamp("2025-04-13"), pd.Timestamp("2025-05-12"), "POST_PLACEBO"))
    recs = []
    for pname, ps, pe, regime in regime_periods:
        ndays = (pe - ps).days + 1
        ar = (act_df["pcd"] < ps) & (act_df[scol].isna() | (act_df[scol] >= ps))
        sub = act_df.loc[ar, ["lead_sponsor_clean", "is_industry", "pcd", "due_date", scol]].copy()
        sub["period"] = pname; sub["regime"] = regime
        sub["is_overdue"] = (sub["due_date"] < ps).astype(int)
        sub["y"] = ((sub[scol] >= ps) & (sub[scol] <= pe)).astype(int)
        sub["exposure"] = ndays / 30.4375
        age_m = ((ps - sub["pcd"]).dt.days / 30.4375).clip(lower=0.0)
        sub["pcd_age_bin"] = pd.cut(
            age_m, bins=[-0.01, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 24.0, 36.0, 60.0, 999.0],
            labels=["0_3", "3_6", "6_9", "9_12", "12_15", "15_18", "18_24", "24_36", "36_60", "60p"]
        ).astype(str)
        recs.append(sub)
    rpan = pd.concat(recs, ignore_index=True)
    rcdf = rpan.groupby(["lead_sponsor_clean", "is_overdue", "is_industry", "pcd_age_bin", "period", "regime"], as_index=False).agg(events=("y", "sum"), risk_exp=("exposure", "sum"))
    rcdf["cell"] = rcdf["is_overdue"].astype(str) + "_" + rcdf["is_industry"].astype(str) + "_" + rcdf["pcd_age_bin"]
    rcdf = rcdf[rcdf.groupby("cell")["events"].transform("sum") > 0].copy().reset_index(drop=True)
    rgs = ["PRE_MAR", "MAR30_31_Q1END", "APR01_12_PURE", "POST_PLACEBO"]
    for rg in rgs:
        rcdf[f"ddd_{rg}"] = ((rcdf["is_overdue"] == 1) & (rcdf["is_industry"] == 1) & (rcdf["regime"] == rg)).astype(float)
        rcdf[f"od_{rg}"] = ((rcdf["is_overdue"] == 1) & (rcdf["regime"] == rg)).astype(float)
        rcdf[f"in_{rg}"] = ((rcdf["is_industry"] == 1) & (rcdf["regime"] == rg)).astype(float)
    rci = rcdf[rcdf["is_industry"] == 1].copy()
    rci = rci[rci.groupby("cell")["events"].transform("sum") > 0].copy().reset_index(drop=True)
    f_ind = smf.glm("events ~ " + " + ".join([f"od_{rg}" for rg in rgs]) + " + C(cell) + C(period)", data=rci, family=sm.families.Poisson(), offset=np.log(rci["risk_exp"].values)).fit(cov_type="cluster", cov_kwds={"groups": rci["lead_sponsor_clean"]})
    rca = rcdf[rcdf["is_industry"] == 0].copy()
    rca = rca[rca.groupby("cell")["events"].transform("sum") > 0].copy().reset_index(drop=True)
    f_acad = smf.glm("events ~ " + " + ".join([f"od_{rg}" for rg in rgs]) + " + C(cell) + C(period)", data=rca, family=sm.families.Poisson(), offset=np.log(rca["risk_exp"].values)).fit(cov_type="cluster", cov_kwds={"groups": rca["lead_sponsor_clean"]})
    out = {}
    for rg in rgs:
        bi, sei, pi = float(f_ind.params[f"od_{rg}"]), float(f_ind.bse[f"od_{rg}"]), float(f_ind.pvalues[f"od_{rg}"])
        ba, sea, pa = float(f_acad.params[f"od_{rg}"]), float(f_acad.bse[f"od_{rg}"]), float(f_acad.pvalues[f"od_{rg}"])
        out[rg] = {
            "ind_od_ot": {"irr": round(np.exp(bi), 4), "ci95": [round(np.exp(bi-1.96*sei), 4), round(np.exp(bi+1.96*sei), 4)], "p": pi},
            "acad_od_ot": {"irr": round(np.exp(ba), 4), "ci95": [round(np.exp(ba-1.96*sea), 4), round(np.exp(ba+1.96*sea), 4)], "p": pa},
        }
    return out

reg_combined = run_q1_donut_reg(act, "sub_combined_proxy", "Combined Proxy (Posted + Unposted Tracking)")
reg_posted = run_q1_donut_reg(act, posted_col, f"Posted-Only ({posted_col})") if posted_col else None
placebo_2025 = run_historical_2025_placebo(act, "sub_combined_proxy")

out = {
    "classification": "POST_HOC_EXPLORATORY_FALSIFICATION",
    "sponsor_fe_causal_boundary": (
        "Sponsor x Month (or Sponsor x Period) fixed effects absorb ONLY shocks common across all status cells "
        "(Overdue and On-Time) within the same sponsor-period. Status-differential shocks within a sponsor-period, "
        "dynamic risk-set selection (into Overdue vs On-Time and into dual-status multi-trial firms), and reporting/eligibility "
        "measurement error remain unabsorbed; thus neither a significant nor a null Sponsor x Period FE estimate alone guarantees causal identification."
    ),
    "q1_calendar_placebo_2024_2025_vs_2026": q1_placebo_rows,
    "other_quarter_ends_2024_2025": other_qe_rows,
    "q1_donut_2026_combined_proxy": reg_combined,
    "q1_donut_2026_posted_only": reg_posted,
    "historical_placebo_2025_same_calendar_dates": placebo_2025,
}

out_path = os.path.join(OUT_DIR, "b2_attack2d_q1_placebo_and_measurement.json")
with open(out_path, "w") as f:
    json.dump(out, f, indent=2)
print("WROTE:", out_path)
