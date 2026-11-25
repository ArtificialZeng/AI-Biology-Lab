#!/opt/anaconda3/bin/python3
"""
Full Publication-Grade Econometric & Mechanism Analysis Engine for Study B-2:
March 30, 2026 FDA Mass Clinical-Trial Delinquency Enforcement Sweep.

Implements pre-registered protocol PREREG_B2.md:
  - Exact trial-level delinquent risk-set construction (monthly & exact Monday-anchored weekly)
  - Primary DiD & Triple-Difference (DDD) across 3 exogenous exempt controls + pooled control
  - LPM (sponsor-clustered SEs, TWFE, Sponsor FE) + Poisson Exposure-Offset Hazard (IRR) + Cox PH (HR)
  - Dynamic Monthly & Biweekly Event-Study with joint Wald pre-trend tests
  - Historical Placebo Sweep Falsifications (March 30, 2025 & March 30, 2024)
  - Cutoff Sensitivity (Aug 31, Sep 30, Oct 15, 2026) & Alternative ACT Proxy (US Phase 2/3)
  - Mechanism & Heterogeneity Decomposition (Sponsor Class, Portfolio Size, Backlog, Drug vs Device, Stale Depth, QC Lag)
  - Adverse Selection / Market-for-Lemons Test (Termination, Primary p>=0.05, P-Value Omission, SAE Rate, Mortality, Within-Sponsor FE)
"""
import json
import os
import time
import warnings
import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.duration.hazard_regression import PHReg

warnings.filterwarnings("ignore")

BASE = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_b2_fda_2026_sweep"
CLEAN_DIR = os.path.join(BASE, "data", "clean")
RES_DIR = os.path.join(BASE, "results")
os.makedirs(RES_DIR, exist_ok=True)


def parse_ctgov_date(series, conservative_end=False):
    """Parse ClinicalTrials.gov YYYY-MM-DD or YYYY-MM dates."""
    s = series.astype(str).str.strip()
    s = s.replace({"None": np.nan, "nan": np.nan, "": np.nan})
    # Exact YYYY-MM-DD
    is_ym = s.str.match(r"^\d{4}-\d{2}$", na=False)
    is_ymd = s.str.match(r"^\d{4}-\d{2}-\d{2}$", na=False)
    out = pd.Series(pd.NaT, index=series.index)
    if is_ymd.any():
        out.loc[is_ymd] = pd.to_datetime(s.loc[is_ymd], format="%Y-%m-%d", errors="coerce")
    if is_ym.any():
        dt_ym = pd.to_datetime(s.loc[is_ym] + "-01", format="%Y-%m-%d", errors="coerce")
        if conservative_end:
            dt_ym = dt_ym + pd.offsets.MonthEnd(0)
        out.loc[is_ym] = dt_ym
    return out


def classify_sponsor_group(cls):
    if cls == "INDUSTRY":
        return "INDUSTRY"
    if cls in ("NIH", "FED"):
        return "NIH_FED"
    return "ACADEMIC_OTHER"


def load_and_prepare_universe():
    df = pd.read_parquet(os.path.join(CLEAN_DIR, "b2_trial_universe.parquet"))
    # Deduplicate if any NCT somehow appeared in >1 arm (priority: us_fdareg_nonp1 first)
    arm_prio = {"us_fdareg_nonp1": 0, "us_p1_only": 1, "us_notfdareg": 2, "nonus_p23_nonfda": 3}
    df["_prio"] = df["arm"].map(arm_prio).fillna(9)
    df = df.sort_values("_prio").drop_duplicates("nct_id", keep="first").drop(columns=["_prio"]).reset_index(drop=True)

    df["pcd"] = parse_ctgov_date(df["primary_completion_date"], conservative_end=True)
    df["sub_date"] = parse_ctgov_date(df["results_first_submit_date"], conservative_end=False)
    df["qc_date"] = parse_ctgov_date(df["results_first_submit_qc_date"], conservative_end=False)
    df["post_date"] = parse_ctgov_date(df["results_first_post_date"], conservative_end=False)
    df["due_date"] = df["pcd"] + pd.Timedelta(days=365)
    df["comp_year"] = df["pcd"].dt.year

    # Restrict strictly to valid primary completion dates in statutory window
    df = df[df["pcd"].notna() & (df["pcd"] >= "2017-01-18") & (df["pcd"] <= "2025-02-28")].copy()

    df["is_industry"] = (df["lead_sponsor_class"] == "INDUSTRY").astype(int)
    df["sponsor_group"] = df["lead_sponsor_class"].apply(classify_sponsor_group)
    df["lead_sponsor_clean"] = df["lead_sponsor_name"].fillna("UNKNOWN_SPONSOR").str.strip()
    df["is_act"] = (df["arm"] == "us_fdareg_nonp1").astype(int)
    df["is_terminated"] = (df["overall_status"] == "TERMINATED").astype(int)
    df["is_p23"] = df["phases"].str.contains("PHASE2|PHASE3", na=False).astype(int)
    df["is_fda_drug_only"] = ((df["is_fda_drug"] == True) & (df["is_fda_device"] != True)).astype(int)
    df["is_fda_device_only"] = ((df["is_fda_device"] == True) & (df["is_fda_drug"] != True)).astype(int)
    df["has_fdaaa_violation"] = (df["fdaaa801_violation"] == True).astype(int)

    # Sponsor-level FDA portfolio & March 29, 2026 delinquent backlog
    act_sub = df[df["arm"] == "us_fdareg_nonp1"]
    spon_total = act_sub.groupby("lead_sponsor_clean")["nct_id"].count().rename("spon_fda_total")
    delinq_mar29 = act_sub[
        (act_sub["due_date"] < pd.Timestamp("2026-03-30")) &
        (act_sub["sub_date"].isna() | (act_sub["sub_date"] >= pd.Timestamp("2026-03-30")))
    ]
    spon_backlog = delinq_mar29.groupby("lead_sponsor_clean")["nct_id"].count().rename("spon_backlog_mar26")
    df = df.merge(spon_total, left_on="lead_sponsor_clean", right_index=True, how="left")
    df = df.merge(spon_backlog, left_on="lead_sponsor_clean", right_index=True, how="left")
    df["spon_fda_total"] = df["spon_fda_total"].fillna(0).astype(int)
    df["spon_backlog_mar26"] = df["spon_backlog_mar26"].fillna(0).astype(int)

    # Exact QC lag for submitted trials
    df["qc_lag_days"] = (df["qc_date"] - df["sub_date"]).dt.days
    df["post_lag_days"] = (df["post_date"] - df["sub_date"]).dt.days
    df["sub_delay_days"] = (df["sub_date"] - df["pcd"]).dt.days
    return df


def build_monthly_risk_panel(df, start_ym="2025-04", end_ym="2026-09", event_ym="2026-04", drop_ym=None):
    """
    Build exact delinquent trial-month risk-set panel.
    Trial i enters month m iff due_date < month_start AND (sub_date is NaT or sub_date >= month_start).
    Outcome y = 1 iff sub_date falls within month m.
    """
    months = pd.period_range(start_ym, end_ym, freq="M")
    if drop_ym:
        months = [m for m in months if str(m) not in set(drop_ym)]
    records = []
    # Vectorized over months (fast for ~50k trials x 18 months!)
    for m in months:
        m_start = m.start_time
        m_end = m.end_time
        at_risk = (df["due_date"] < m_start) & (df["sub_date"].isna() | (df["sub_date"] >= m_start))
        sub_m = df.loc[at_risk, [
            "nct_id", "arm", "is_act", "is_industry", "sponsor_group", "lead_sponsor_clean",
            "comp_year", "is_p23", "is_terminated", "is_fda_drug", "is_fda_device",
            "spon_fda_total", "spon_backlog_mar26", "due_date", "sub_date", "enrollment"
        ]].copy()
        sub_m["ym"] = str(m)
        sub_m["post"] = int(str(m) >= event_ym)
        sub_m["y"] = ((sub_m["sub_date"] >= m_start) & (sub_m["sub_date"] <= m_end)).astype(int)
        sub_m["overdue_years"] = ((m_start - sub_m["due_date"]).dt.days / 365.25).clip(lower=0)
        records.append(sub_m)
    panel = pd.concat(records, ignore_index=True)
    return panel


def build_weekly_risk_panel(df, event_date="2026-03-30", weeks_pre=52, weeks_post=26):
    """
    Exact Monday-anchored weekly risk-set panel where week k=0 starts EXACTLY on Monday 2026-03-30!
    No boundary contamination at all (2026-03-23..29 is k=-1; 2026-03-30..04-05 is k=0).
    """
    t0 = pd.Timestamp(event_date)
    records = []
    for k in range(-weeks_pre, weeks_post):
        w_start = t0 + pd.Timedelta(days=7 * k)
        w_end = w_start + pd.Timedelta(days=6, hours=23, minutes=59, seconds=59)
        at_risk = (df["due_date"] < w_start) & (df["sub_date"].isna() | (df["sub_date"] >= w_start))
        sub_w = df.loc[at_risk, [
            "nct_id", "arm", "is_act", "is_industry", "sponsor_group", "lead_sponsor_clean",
            "comp_year", "is_p23", "spon_fda_total", "spon_backlog_mar26", "sub_date"
        ]].copy()
        sub_w["rel_week"] = k
        sub_w["week_start"] = w_start.strftime("%Y-%m-%d")
        sub_w["post"] = int(k >= 0)
        sub_w["y"] = ((sub_w["sub_date"] >= w_start) & (sub_w["sub_date"] <= w_end)).astype(int)
        records.append(sub_w)
    return pd.concat(records, ignore_index=True)


def run_ddd_suite(panel, ctrl_arms, label):
    """
    Run comprehensive DiD & DDD regressions comparing us_fdareg_nonp1 vs ctrl_arms:
      (1) LPM DiD (Pooled Industry+Academic) with Arm x CompYear + YM FE, sponsor-clustered SE
      (2) LPM DDD (Basic): Arm x Industry x Post + Arm x Industry + YM FE, sponsor-clustered SE
      (3) LPM DDD (Cohort FE): Arm x Industry x CompYear FE + YM FE, sponsor-clustered SE
      (4) LPM DDD (Industry-only DiD & Academic-only DiD separately!)
      (5) Poisson Exposure-Offset Hazard DDD (multiplicative IRR!) with Arm x Industry x CompYear + YM FE
    """
    sub = panel[panel["arm"].isin(["us_fdareg_nonp1"] + list(ctrl_arms))].copy()
    sub["act_post"] = sub["is_act"] * sub["post"]
    sub["ind_post"] = sub["is_industry"] * sub["post"]
    sub["act_ind"] = sub["is_act"] * sub["is_industry"]
    sub["ddd"] = sub["is_act"] * sub["is_industry"] * sub["post"]
    sub["y_bp"] = sub["y"] * 100.0  # percentage points per month (easier to read: 1.0 = 1.00 pp/month)
    sub["cell_ac"] = sub["arm"].astype(str) + "_" + sub["comp_year"].astype(str)
    sub["cell_aic"] = sub["arm"].astype(str) + "_" + sub["is_industry"].astype(str) + "_" + sub["comp_year"].astype(str)

    out = {
        "label": label,
        "ctrl_arms": list(ctrl_arms),
        "n_obs": int(len(sub)),
        "n_trials": int(sub["nct_id"].nunique()),
        "n_sponsors": int(sub["lead_sponsor_clean"].nunique()),
        "n_events": int(sub["y"].sum()),
    }

    # Cell summary hazards (per 100 delinquent trial-months)
    grp = sub.groupby(["is_act", "is_industry", "post"])["y"].agg(["count", "sum", "mean"]).reset_index()
    out["cell_hazards"] = [
        {
            "is_act": int(r.is_act), "is_industry": int(r.is_industry), "post": int(r.post),
            "risk_months": int(r["count"]), "submissions": int(r["sum"]),
            "hazard_pct": round(float(r["mean"]) * 100.0, 4)
        }
        for _, r in grp.iterrows()
    ]

    # (1) Pooled DiD LPM
    m_did = smf.ols("y_bp ~ act_post + C(cell_ac) + C(ym)", data=sub).fit(
        cov_type="cluster", cov_kwds={"groups": sub["lead_sponsor_clean"]}
    )
    out["lpm_did_pooled"] = {
        "coef_pp": float(m_did.params["act_post"]),
        "se_pp": float(m_did.bse["act_post"]),
        "t": float(m_did.tvalues["act_post"]),
        "p": float(m_did.pvalues["act_post"]),
        "ci95": [float(x) for x in m_did.conf_int().loc["act_post"]],
    }

    # (2) Industry-only DiD LPM & Academic-only DiD LPM
    for tag, ind_val in [("lpm_did_industry", 1), ("lpm_did_academic", 0)]:
        s2 = sub[sub["is_industry"] == ind_val]
        m2 = smf.ols("y_bp ~ act_post + C(cell_ac) + C(ym)", data=s2).fit(
            cov_type="cluster", cov_kwds={"groups": s2["lead_sponsor_clean"]}
        )
        out[tag] = {
            "coef_pp": float(m2.params["act_post"]),
            "se_pp": float(m2.bse["act_post"]),
            "t": float(m2.tvalues["act_post"]),
            "p": float(m2.pvalues["act_post"]),
            "ci95": [float(x) for x in m2.conf_int().loc["act_post"]],
            "n_obs": int(len(s2)),
            "n_events": int(s2["y"].sum()),
        }

    # (3) Full DDD LPM (Cohort x Arm x Industry FE + YM FE)
    m_ddd = smf.ols("y_bp ~ ddd + act_post + ind_post + C(cell_aic) + C(ym)", data=sub).fit(
        cov_type="cluster", cov_kwds={"groups": sub["lead_sponsor_clean"]}
    )
    out["lpm_ddd_cohort_fe"] = {
        term: {
            "coef_pp": float(m_ddd.params[term]),
            "se_pp": float(m_ddd.bse[term]),
            "t": float(m_ddd.tvalues[term]),
            "p": float(m_ddd.pvalues[term]),
            "ci95": [float(x) for x in m_ddd.conf_int().loc[term]],
        }
        for term in ("ddd", "act_post", "ind_post")
    }

    # (4) Sponsor Fixed Effects DDD via Frisch-Waugh-Lovell exact demeaning (memory-efficient for thousands of sponsors!)
    # Demean within sponsor on trials belonging to sponsors with >=2 obs, plus YM FE via dummies
    # Even cleaner: within-sponsor demeaning of y_bp and covariates after including YM dummies -> exact FE OLS!
    spon_counts = sub.groupby("lead_sponsor_clean")["y"].transform("count")
    sub_fe = sub[spon_counts >= 2].copy()
    ym_dummies = pd.get_dummies(sub_fe["ym"], drop_first=True, dtype=float)
    X_raw = pd.concat([sub_fe[["ddd", "act_post", "ind_post"]].astype(float), ym_dummies], axis=1)
    y_raw = sub_fe["y_bp"].astype(float)
    # Within-sponsor demean
    g = sub_fe["lead_sponsor_clean"]
    X_dm = X_raw - X_raw.groupby(g).transform("mean")
    y_dm = y_raw - y_raw.groupby(g).transform("mean")
    m_sfe = sm.OLS(y_dm, X_dm).fit(cov_type="cluster", cov_kwds={"groups": g})
    # Adjust dof for absorbed sponsor FE (N_g)
    n_g = g.nunique()
    n_obs = len(sub_fe)
    k_reg = X_dm.shape[1]
    dof_scale = np.sqrt((n_obs - k_reg) / max(n_obs - k_reg - n_g, 1))
    out["lpm_ddd_sponsor_fe"] = {}
    for term in ("ddd", "act_post", "ind_post"):
        b = float(m_sfe.params[term])
        se = float(m_sfe.bse[term]) * dof_scale
        tval = b / se if se > 0 else np.nan
        pval = float(2 * stats.norm.sf(abs(tval)))
        out["lpm_ddd_sponsor_fe"][term] = {
            "coef_pp": b, "se_pp": se, "t": tval, "p": pval,
            "ci95": [b - 1.96 * se, b + 1.96 * se],
        }
    out["lpm_ddd_sponsor_fe"]["n_obs"] = int(n_obs)
    out["lpm_ddd_sponsor_fe"]["n_sponsors_fe"] = int(n_g)

    # (5) Poisson Exposure-Offset Hazard Model (exact multiplicative IRR!)
    # Aggregate to sponsor-arm-cohort-month cells so we can cluster at sponsor level AND avoid huge memory
    cell_df = sub.groupby(
        ["lead_sponsor_clean", "arm", "is_act", "is_industry", "comp_year", "ym", "post"], as_index=False
    ).agg(events=("y", "sum"), risk=("y", "count"))
    cell_df["act_post"] = cell_df["is_act"] * cell_df["post"]
    cell_df["ind_post"] = cell_df["is_industry"] * cell_df["post"]
    cell_df["ddd"] = cell_df["is_act"] * cell_df["is_industry"] * cell_df["post"]
    cell_df["cell_aic"] = cell_df["arm"].astype(str) + "_" + cell_df["is_industry"].astype(str) + "_" + cell_df["comp_year"].astype(str)
    # Drop cells whose fixed-effect group has 0 total events (separation in Poisson)
    fe_ev = cell_df.groupby("cell_aic")["events"].transform("sum")
    cell_p = cell_df[fe_ev > 0].copy()
    m_pois = smf.glm(
        "events ~ ddd + act_post + ind_post + C(cell_aic) + C(ym)",
        data=cell_p,
        family=sm.families.Poisson(),
        offset=np.log(cell_p["risk"].values)
    ).fit(cov_type="cluster", cov_kwds={"groups": cell_p["lead_sponsor_clean"]})
    out["poisson_ddd"] = {}
    for term in ("ddd", "act_post", "ind_post"):
        b = float(m_pois.params[term])
        se = float(m_pois.bse[term])
        ci = m_pois.conf_int().loc[term]
        out["poisson_ddd"][term] = {
            "log_irr": b, "se": se, "z": float(m_pois.tvalues[term]), "p": float(m_pois.pvalues[term]),
            "irr": float(np.exp(b)), "irr_ci95": [float(np.exp(ci[0])), float(np.exp(ci[1]))],
        }
    # Industry-specific & Academic-specific Poisson IRRs
    for tag, ind_val in [("poisson_did_industry", 1), ("poisson_did_academic", 0)]:
        cp = cell_df[cell_df["is_industry"] == ind_val].copy()
        cp["cell_ac"] = cp["arm"].astype(str) + "_" + cp["comp_year"].astype(str)
        cp = cp[cp.groupby("cell_ac")["events"].transform("sum") > 0].copy()
        mp = smf.glm(
            "events ~ act_post + C(cell_ac) + C(ym)",
            data=cp, family=sm.families.Poisson(), offset=np.log(cp["risk"].values)
        ).fit(cov_type="cluster", cov_kwds={"groups": cp["lead_sponsor_clean"]})
        b = float(mp.params["act_post"])
        ci = mp.conf_int().loc["act_post"]
        out[tag] = {
            "log_irr": b, "se": float(mp.bse["act_post"]), "z": float(mp.tvalues["act_post"]),
            "p": float(mp.pvalues["act_post"]), "irr": float(np.exp(b)),
            "irr_ci95": [float(np.exp(ci[0])), float(np.exp(ci[1]))],
        }

    return out


def run_cox_timevarying(df, ctrl_arms, start_date="2025-04-01", event_date="2026-03-30", end_date="2026-09-30"):
    """
    Exact continuous-time Cox Proportional Hazards model with time-varying treatment (episode splitting at 2026-03-30).
    Uses calendar time (days since 2025-04-01) with left truncation at max(start_date, due_date).
    """
    t_start = pd.Timestamp(start_date)
    t_evt = pd.Timestamp(event_date)
    t_end = pd.Timestamp(end_date)
    cut_day = (t_evt - t_start).days
    end_day = (t_end - t_start).days

    sub = df[df["arm"].isin(["us_fdareg_nonp1"] + list(ctrl_arms))].copy()
    # Must enter delinquent risk set before end_date and not submit before start_date
    sub["entry_dt"] = sub["due_date"].clip(lower=t_start)
    sub = sub[(sub["entry_dt"] < t_end) & (sub["sub_date"].isna() | (sub["sub_date"] >= sub["entry_dt"]))].copy()
    sub["exit_dt"] = sub["sub_date"].fillna(t_end).clip(upper=t_end)
    sub["event"] = ((sub["sub_date"].notna()) & (sub["sub_date"] >= sub["entry_dt"]) & (sub["sub_date"] <= t_end)).astype(int)

    sub["t0"] = (sub["entry_dt"] - t_start).dt.days.astype(float)
    sub["t1"] = (sub["exit_dt"] - t_start).dt.days.astype(float)
    # Avoid zero-length intervals (submitted on exact entry day -> +0.5d)
    sub.loc[sub["t1"] <= sub["t0"], "t1"] = sub.loc[sub["t1"] <= sub["t0"], "t0"] + 0.5

    # Split intervals crossing cut_day (2026-03-30) into pre and post episodes
    pre_only = sub[sub["t1"] <= cut_day].copy()
    pre_only["post"] = 0
    post_only = sub[sub["t0"] >= cut_day].copy()
    post_only["post"] = 1
    cross = sub[(sub["t0"] < cut_day) & (sub["t1"] > cut_day)].copy()
    c_pre = cross.copy()
    c_pre["t1"] = float(cut_day)
    c_pre["event"] = 0
    c_pre["post"] = 0
    c_post = cross.copy()
    c_post["t0"] = float(cut_day)
    c_post["post"] = 1

    ep = pd.concat([pre_only, post_only, c_pre, c_post], ignore_index=True)
    ep["act_ind"] = ep["is_act"] * ep["is_industry"]
    ep["act_post"] = ep["is_act"] * ep["post"]
    ep["ind_post"] = ep["is_industry"] * ep["post"]
    ep["ddd"] = ep["is_act"] * ep["is_industry"] * ep["post"]

    # Stratify by completion cohort (comp_year) so each cohort has an arbitrary nonparametric
    # calendar-time baseline hazard h_{0,c}(t), while Arm x Industry main effects and
    # time-varying Post interactions (act_post, ind_post, ddd) are identified across arms at each day t!
    names = ["is_act", "is_industry", "act_ind", "act_post", "ind_post", "ddd"]
    X = ep[names].astype(float).values
    mod = PHReg(ep["t1"].values, X, status=ep["event"].values, entry=ep["t0"].values, strata=ep["comp_year"].values)
    res = mod.fit()
    out = {
        "n_episodes": int(len(ep)),
        "n_trials": int(sub["nct_id"].nunique()),
        "n_events": int(ep["event"].sum()),
        "terms": {}
    }
    for idx, nm in enumerate(names):
        b = float(res.params[idx])
        se = float(res.bse[idx])
        z = b / se
        p = float(2 * stats.norm.sf(abs(z)))
        out["terms"][nm] = {
            "log_hr": b, "se": se, "z": z, "p": p,
            "hr": float(np.exp(b)), "hr_ci95": [float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se))],
        }
    # Also compute total Industry ACT Post HR = exp(act_post + ddd) with delta-method SE
    cov = res.cov_params()
    i_a, i_d = names.index("act_post"), names.index("ddd")
    b_tot = float(res.params[i_a] + res.params[i_d])
    se_tot = float(np.sqrt(max(cov[i_a, i_a] + cov[i_d, i_d] + 2 * cov[i_a, i_d], 0.0)))
    z_tot = b_tot / se_tot
    out["terms"]["industry_act_total_post"] = {
        "log_hr": b_tot, "se": se_tot, "z": z_tot, "p": float(2 * stats.norm.sf(abs(z_tot))),
        "hr": float(np.exp(b_tot)), "hr_ci95": [float(np.exp(b_tot - 1.96 * se_tot)), float(np.exp(b_tot + 1.96 * se_tot))],
    }
    return out


def run_event_study_monthly(panel, ctrl_arms, ref_ym="2026-02"):
    """Monthly dynamic event-study DDD and Industry-only DiD with joint pre-trend Wald test."""
    sub = panel[panel["arm"].isin(["us_fdareg_nonp1"] + list(ctrl_arms))].copy()
    sub["y_bp"] = sub["y"] * 100.0
    sub["cell_aic"] = sub["arm"].astype(str) + "_" + sub["is_industry"].astype(str) + "_" + sub["comp_year"].astype(str)
    months = sorted(sub["ym"].unique())
    ref_idx = months.index(ref_ym)

    # Construct interaction columns explicitly
    ddd_cols = []
    act_cols = []
    ind_cols = []
    for m in months:
        if m == ref_ym:
            continue
        k = months.index(m) - ref_idx
        tag = f"m{abs(k)}" if k < 0 else f"p{k}"
        dm = (sub["ym"] == m).astype(float)
        sub[f"ddd_{tag}"] = sub["is_act"] * sub["is_industry"] * dm
        sub[f"act_{tag}"] = sub["is_act"] * dm
        sub[f"ind_{tag}"] = sub["is_industry"] * dm
        ddd_cols.append((m, k, f"ddd_{tag}"))
        act_cols.append((m, k, f"act_{tag}"))
        ind_cols.append((m, k, f"ind_{tag}"))

    rhs = " + ".join([c[2] for c in ddd_cols] + [c[2] for c in act_cols] + [c[2] for c in ind_cols]) + " + C(cell_aic) + C(ym)"
    fit = smf.ols(f"y_bp ~ {rhs}", data=sub).fit(
        cov_type="cluster", cov_kwds={"groups": sub["lead_sponsor_clean"]}
    )

    series = []
    pre_ddd_names = []
    pre_act_names = []
    for m in months:
        k = months.index(m) - ref_idx
        if m == ref_ym:
            series.append({"ym": m, "rel_m": k, "ddd_coef": 0.0, "ddd_se": 0.0, "ddd_p": 1.0,
                           "act_acad_coef": 0.0, "act_acad_se": 0.0, "act_ind_total": 0.0, "act_ind_se": 0.0})
        else:
            tag = f"m{abs(k)}" if k < 0 else f"p{k}"
            dn = f"ddd_{tag}"
            an = f"act_{tag}"
            if k < 0:
                pre_ddd_names.append(dn)
                pre_act_names.append(an)
            b_ddd = float(fit.params[dn])
            se_ddd = float(fit.bse[dn])
            b_act = float(fit.params[an])
            se_act = float(fit.bse[an])
            # Total Industry ACT lift = act_t + ddd_t
            cov_m = fit.cov_params()
            var_tot = float(cov_m.loc[dn, dn] + cov_m.loc[an, an] + 2 * cov_m.loc[dn, an])
            se_tot = float(np.sqrt(max(var_tot, 0.0)))
            series.append({
                "ym": m, "rel_m": k,
                "ddd_coef": b_ddd, "ddd_se": se_ddd, "ddd_p": float(fit.pvalues[dn]),
                "act_acad_coef": b_act, "act_acad_se": se_act, "act_acad_p": float(fit.pvalues[an]),
                "act_ind_total": b_act + b_ddd, "act_ind_se": se_tot,
                "act_ind_p": float(2 * stats.norm.sf(abs((b_act + b_ddd) / se_tot))) if se_tot > 0 else 1.0,
            })

    # Joint Wald pre-trend test on all k < 0 (excluding ref_ym k=0 and March 2026 k=+1)
    R = np.zeros((len(pre_ddd_names), len(fit.params)))
    param_names = list(fit.params.index)
    for i, nm in enumerate(pre_ddd_names):
        R[i, param_names.index(nm)] = 1.0
    wt = fit.wald_test(R, scalar=True)
    return {
        "ref_ym": ref_ym,
        "pre_trend_wald_ddd": {
            "df": len(pre_ddd_names),
            "stat": float(wt.statistic),
            "p": float(wt.pvalue),
        },
        "series": series,
    }


def run_biweekly_exact_event_study(df, ctrl_arms):
    """Exact biweekly event study anchored on Monday 2026-03-30 (zero March-boundary blurring!)."""
    wpanel = build_weekly_risk_panel(df[df["arm"].isin(["us_fdareg_nonp1"] + list(ctrl_arms))], event_date="2026-03-30", weeks_pre=26, weeks_post=26)
    wpanel["biweek"] = (wpanel["rel_week"] // 2).astype(int)  # -13..-1 pre, 0..+12 post
    wpanel["y_bp"] = wpanel["y"] * 100.0
    wpanel["cell_aic"] = wpanel["arm"].astype(str) + "_" + wpanel["is_industry"].astype(str) + "_" + wpanel["comp_year"].astype(str)

    # Raw biweekly hazard series by group (for exact figure!)
    raw_series = (
        wpanel.groupby(["biweek", "is_act", "is_industry"], as_index=False)
        .agg(events=("y", "sum"), risk=("y", "count"))
    )
    raw_series["hazard_pct_per_week"] = raw_series["events"] / raw_series["risk"] * 100.0

    # Regression event study (ref biweek = -1, i.e., Mar 16–Mar 29, 2026)
    bws = sorted(wpanel["biweek"].unique())
    ref_bw = -1
    ddd_cols, act_cols, ind_cols = [], [], []
    for b in bws:
        if b == ref_bw:
            continue
        tag = f"m{abs(b)}" if b < 0 else f"p{b}"
        db = (wpanel["biweek"] == b).astype(float)
        wpanel[f"ddd_{tag}"] = wpanel["is_act"] * wpanel["is_industry"] * db
        wpanel[f"act_{tag}"] = wpanel["is_act"] * db
        wpanel[f"ind_{tag}"] = wpanel["is_industry"] * db
        ddd_cols.append((b, f"ddd_{tag}"))
        act_cols.append((b, f"act_{tag}"))
        ind_cols.append((b, f"ind_{tag}"))

    rhs = " + ".join([c[1] for c in ddd_cols] + [c[1] for c in act_cols] + [c[1] for c in ind_cols]) + " + C(cell_aic) + C(biweek)"
    fit = smf.ols(f"y_bp ~ {rhs}", data=wpanel).fit(
        cov_type="cluster", cov_kwds={"groups": wpanel["lead_sponsor_clean"]}
    )
    reg_series = []
    pre_names = []
    for b in bws:
        w0 = (pd.Timestamp("2026-03-30") + pd.Timedelta(days=14 * b)).strftime("%Y-%m-%d")
        if b == ref_bw:
            reg_series.append({"biweek": int(b), "start_date": w0, "ddd_coef": 0.0, "ddd_se": 0.0, "ddd_p": 1.0,
                               "act_ind_total": 0.0, "act_ind_se": 0.0, "act_acad_coef": 0.0, "act_acad_se": 0.0})
        else:
            tag = f"m{abs(b)}" if b < 0 else f"p{b}"
            dn, an = f"ddd_{tag}", f"act_{tag}"
            if b < -1:
                pre_names.append(dn)
            b_d, se_d = float(fit.params[dn]), float(fit.bse[dn])
            b_a, se_a = float(fit.params[an]), float(fit.bse[an])
            cm = fit.cov_params()
            se_t = float(np.sqrt(max(cm.loc[dn, dn] + cm.loc[an, an] + 2 * cm.loc[dn, an], 0.0)))
            reg_series.append({
                "biweek": int(b), "start_date": w0,
                "ddd_coef": b_d, "ddd_se": se_d, "ddd_p": float(fit.pvalues[dn]),
                "act_acad_coef": b_a, "act_acad_se": se_a,
                "act_ind_total": b_a + b_d, "act_ind_se": se_t,
            })
    R = np.zeros((len(pre_names), len(fit.params)))
    pnames = list(fit.params.index)
    for i, nm in enumerate(pre_names):
        R[i, pnames.index(nm)] = 1.0
    wt = fit.wald_test(R, scalar=True)
    return {
        "pre_trend_wald_ddd": {"df": len(pre_names), "stat": float(wt.statistic), "p": float(wt.pvalue)},
        "reg_series": reg_series,
        "raw_series": raw_series.to_dict(orient="records"),
    }


def run_mechanism_heterogeneity(df, panel):
    """
    Mechanism & heterogeneity decomposition inside ACT_FDA (and DDD vs us_notfdareg):
      1. Sponsor governance fine breakdown: INDUSTRY vs ACADEMIC_OTHER vs NIH_FED
      2. Portfolio capacity (spon_fda_total: 1-2, 3-19, 20+) within Industry & Academic
      3. Delinquent backlog (spon_backlog_mar26: 1, 2-4, 5+) within Industry & Academic
      4. Modality: FDA Drug/Biologic vs FDA Device
      5. Delinquency depth on Mar 30, 2026: 1-2yr (2023-25), 2-4yr (2021-22), 4+yr (2017-20)
      6. NLM QC lag check (resultsFirstSubmitQcDate - resultsFirstSubmitDate) pre vs post
    """
    out = {}
    act_p = panel[panel["arm"] == "us_fdareg_nonp1"].copy()
    ctrl_p = panel[panel["arm"] == "us_notfdareg"].copy()

    def group_did(sub_act, sub_ctrl):
        comb = pd.concat([sub_act, sub_ctrl], ignore_index=True)
        comb["act_post"] = comb["is_act"] * comb["post"]
        comb["y_bp"] = comb["y"] * 100.0
        comb["cell_ac"] = comb["arm"].astype(str) + "_" + comb["comp_year"].astype(str)
        if comb["y"].sum() < 5 or comb["is_act"].nunique() < 2:
            return None
        m = smf.ols("y_bp ~ act_post + C(cell_ac) + C(ym)", data=comb).fit(
            cov_type="cluster", cov_kwds={"groups": comb["lead_sponsor_clean"]}
        )
        # Pre & post raw hazard in treated subgroup
        pre_h = float(sub_act.loc[sub_act["post"] == 0, "y"].mean()) * 100.0
        post_h = float(sub_act.loc[sub_act["post"] == 1, "y"].mean()) * 100.0
        # Poisson IRR
        cdf = comb.groupby(["lead_sponsor_clean", "arm", "is_act", "comp_year", "ym", "post"], as_index=False).agg(
            events=("y", "sum"), risk=("y", "count")
        )
        cdf["act_post"] = cdf["is_act"] * cdf["post"]
        cdf["cell_ac"] = cdf["arm"].astype(str) + "_" + cdf["comp_year"].astype(str)
        cdf = cdf[cdf.groupby("cell_ac")["events"].transform("sum") > 0]
        try:
            mp = smf.glm("events ~ act_post + C(cell_ac) + C(ym)", data=cdf,
                         family=sm.families.Poisson(), offset=np.log(cdf["risk"].values)).fit(
                cov_type="cluster", cov_kwds={"groups": cdf["lead_sponsor_clean"]}
            )
            irr = float(np.exp(mp.params["act_post"]))
            irr_p = float(mp.pvalues["act_post"])
        except Exception:
            irr, irr_p = np.nan, np.nan
        return {
            "n_act_obs": int(len(sub_act)),
            "n_act_events": int(sub_act["y"].sum()),
            "pre_hazard_pct": round(pre_h, 4),
            "post_hazard_pct": round(post_h, 4),
            "did_pp": float(m.params["act_post"]),
            "se_pp": float(m.bse["act_post"]),
            "p": float(m.pvalues["act_post"]),
            "poisson_irr": irr,
            "poisson_p": irr_p,
        }

    # 1. Fine sponsor class
    out["by_sponsor_group"] = {}
    for sg in ["INDUSTRY", "ACADEMIC_OTHER", "NIH_FED"]:
        out["by_sponsor_group"][sg] = group_did(
            act_p[act_p["sponsor_group"] == sg],
            ctrl_p[ctrl_p["sponsor_group"] == sg] if (ctrl_p["sponsor_group"] == sg).sum() > 200 else ctrl_p
        )

    # 2. Portfolio size (regulatory capacity) inside Industry & Academic
    out["by_portfolio_size"] = {}
    bins = [(1, 2, "small_1_2"), (3, 19, "mid_3_19"), (20, 999999, "large_20plus")]
    for ind_tag, ind_val in [("industry", 1), ("academic", 0)]:
        out["by_portfolio_size"][ind_tag] = {}
        c_sub = ctrl_p[ctrl_p["is_industry"] == ind_val]
        if len(c_sub) < 500:
            c_sub = ctrl_p
        for lo, hi, bname in bins:
            a_sub = act_p[(act_p["is_industry"] == ind_val) & (act_p["spon_fda_total"] >= lo) & (act_p["spon_fda_total"] <= hi)]
            out["by_portfolio_size"][ind_tag][bname] = group_did(a_sub, c_sub)

    # 3. Delinquent backlog size inside Industry & Academic
    out["by_backlog_size"] = {}
    bbins = [(1, 1, "backlog_1"), (2, 4, "backlog_2_4"), (5, 999999, "backlog_5plus")]
    for ind_tag, ind_val in [("industry", 1), ("academic", 0)]:
        out["by_backlog_size"][ind_tag] = {}
        c_sub = ctrl_p[ctrl_p["is_industry"] == ind_val]
        if len(c_sub) < 500:
            c_sub = ctrl_p
        for lo, hi, bname in bbins:
            a_sub = act_p[(act_p["is_industry"] == ind_val) & (act_p["spon_backlog_mar26"] >= lo) & (act_p["spon_backlog_mar26"] <= hi)]
            out["by_backlog_size"][ind_tag][bname] = group_did(a_sub, c_sub)

    # 4. Modality (Drug vs Device) x Industry/Academic
    out["by_modality"] = {}
    for mod_name, mask in [
        ("drug_industry", (act_p["is_fda_drug"] == True) & (act_p["is_industry"] == 1)),
        ("device_industry", (act_p["is_fda_device"] == True) & (act_p["is_fda_drug"] != True) & (act_p["is_industry"] == 1)),
        ("drug_academic", (act_p["is_fda_drug"] == True) & (act_p["is_industry"] == 0)),
        ("device_academic", (act_p["is_fda_device"] == True) & (act_p["is_fda_drug"] != True) & (act_p["is_industry"] == 0)),
    ]:
        out["by_modality"][mod_name] = group_did(act_p[mask], ctrl_p)

    # 5. Delinquency depth (completion cohort)
    out["by_cohort_depth"] = {}
    for ind_tag, ind_val in [("industry", 1), ("academic", 0)]:
        out["by_cohort_depth"][ind_tag] = {}
        for cname, ylo, yhi in [("recent_2023_2025", 2023, 2025), ("mid_2021_2022", 2021, 2022), ("deep_2017_2020", 2017, 2020)]:
            a_sub = act_p[(act_p["is_industry"] == ind_val) & (act_p["comp_year"] >= ylo) & (act_p["comp_year"] <= yhi)]
            c_sub = ctrl_p[(ctrl_p["comp_year"] >= ylo) & (ctrl_p["comp_year"] <= yhi)]
            out["by_cohort_depth"][ind_tag][cname] = group_did(a_sub, c_sub)

    # 6. NLM QC lag diagnostics
    sub_df = df[(df["arm"] == "us_fdareg_nonp1") & df["sub_date"].notna()].copy()
    pre_s = sub_df[(sub_df["sub_date"] >= "2025-04-01") & (sub_df["sub_date"] < "2026-03-30")]
    post_s = sub_df[(sub_df["sub_date"] >= "2026-03-30") & (sub_df["sub_date"] <= "2026-09-30")]
    out["qc_lag_diagnostics"] = {
        "pre_n": int(len(pre_s)),
        "pre_median_post_lag_days": float(pre_s["post_lag_days"].median()),
        "pre_p90_post_lag_days": float(pre_s["post_lag_days"].quantile(0.90)),
        "post_n": int(len(post_s)),
        "post_median_post_lag_days": float(post_s["post_lag_days"].median()),
        "post_p90_post_lag_days": float(post_s["post_lag_days"].quantile(0.90)),
    }
    return out


def run_lemons_analysis():
    """
    H4 Adverse Selection / Market-for-Lemons Test on ALL 2024-2026 US FDA-regulated results submissions.
    Compares:
      - Post-Sweep Overdue (submitted 2026-03-30..2026-09-30, delay > 365d)
      - Pre-Sweep Overdue  (submitted 2024-01-01..2026-03-29, delay > 365d)
      - On-Time            (submitted 2024-01-01..2026-09-30, delay <= 365d)
    """
    lem = pd.read_parquet(os.path.join(CLEAN_DIR, "b2_lemons_results.parquet"))
    lem["pcd"] = parse_ctgov_date(lem["primary_completion_date"], conservative_end=True)
    lem["sub_date"] = parse_ctgov_date(lem["results_first_submit_date"], conservative_end=False)
    lem = lem[lem["pcd"].notna() & lem["sub_date"].notna()].copy()
    # Clean cutoff 2024-01-01..2026-09-30 to avoid QC right-censoring
    lem = lem[(lem["sub_date"] >= "2024-01-01") & (lem["sub_date"] <= "2026-09-30")].copy()
    lem["delay_days"] = (lem["sub_date"] - lem["pcd"]).dt.days
    lem["is_overdue"] = (lem["delay_days"] > 365).astype(int)
    lem["is_post_sweep"] = (lem["sub_date"] >= "2026-03-30").astype(int)
    lem["is_industry"] = (lem["lead_sponsor_class"] == "INDUSTRY").astype(int)
    lem["is_terminated"] = (lem["overall_status"] == "TERMINATED").astype(int)
    lem["is_fda_drug_ind"] = (lem["is_fda_drug"] == True).astype(int)
    lem["log_enroll"] = np.log1p(lem["enrollment"].fillna(0).clip(lower=0))
    lem["omitted_pval"] = ((lem["n_primary_outcomes"] > 0) & (~lem["has_primary_pval"])).astype(int)
    lem["nonsig_p05"] = lem["all_nonsig_primary_p05"].astype(float)
    lem["lead_sponsor_clean"] = lem["lead_sponsor_name"].fillna("UNKNOWN").str.strip()
    lem["phase_clean"] = lem["phases"].fillna("NA").replace("", "NA")
    # Group classification
    def grp_cat(r):
        if r["is_overdue"] == 0:
            return "on_time"
        return "post_sweep_overdue" if r["is_post_sweep"] == 1 else "pre_sweep_overdue"
    lem["disclosure_group"] = lem.apply(grp_cat, axis=1)

    out = {"n_total": int(len(lem)), "group_summary": {}, "regressions": {}}
    for gname, gdf in lem.groupby("disclosure_group"):
        pv_sub = gdf[gdf["nonsig_p05"].notna()]
        sae_sub = gdf[gdf["sae_rate"].notna()]
        mort_sub = gdf[gdf["mortality_rate"].notna()]
        out["group_summary"][gname] = {
            "n": int(len(gdf)),
            "industry_share": float(gdf["is_industry"].mean()),
            "median_delay_days": float(gdf["delay_days"].median()),
            "mean_delay_days": float(gdf["delay_days"].mean()),
            "terminated_share": float(gdf["is_terminated"].mean()),
            "omitted_pval_share": float(gdf["omitted_pval"].mean()),
            "n_with_pval": int(len(pv_sub)),
            "nonsig_p05_share": float(pv_sub["nonsig_p05"].mean()) if len(pv_sub) else None,
            "n_with_sae": int(len(sae_sub)),
            "mean_sae_rate": float(sae_sub["sae_rate"].mean()) if len(sae_sub) else None,
            "median_sae_rate": float(sae_sub["sae_rate"].median()) if len(sae_sub) else None,
            "mean_mortality_rate": float(mort_sub["mortality_rate"].mean()) if len(mort_sub) else None,
        }

    # Regression comparisons:
    # (A) All overdue + on-time: compare pre_sweep_overdue & post_sweep_overdue vs on_time (and direct contrast post_sweep_overdue vs pre_sweep_overdue!)
    lem["d_pre_overdue"] = (lem["disclosure_group"] == "pre_sweep_overdue").astype(int)
    lem["d_post_overdue"] = (lem["disclosure_group"] == "post_sweep_overdue").astype(int)

    outcomes = [
        ("is_terminated", lem),
        ("omitted_pval", lem[lem["n_primary_outcomes"] > 0]),
        ("nonsig_p05", lem[lem["nonsig_p05"].notna()]),
        ("sae_rate", lem[lem["sae_rate"].notna()]),
        ("mortality_rate", lem[lem["mortality_rate"].notna()]),
        ("delay_days", lem[lem["is_overdue"] == 1]),
    ]

    for ycol, sub_df in outcomes:
        out["regressions"][ycol] = {}
        if ycol == "delay_days":
            m = smf.ols("delay_days ~ is_post_sweep + is_industry + is_fda_drug_ind + log_enroll + C(phase_clean)", data=sub_df).fit(
                cov_type="cluster", cov_kwds={"groups": sub_df["lead_sponsor_clean"]}
            )
            out["regressions"][ycol]["post_vs_pre_overdue"] = {
                "coef": float(m.params["is_post_sweep"]), "se": float(m.bse["is_post_sweep"]),
                "p": float(m.pvalues["is_post_sweep"]), "n": int(len(sub_df))
            }
            continue
        # 1. Cross-sectional with controls (vs on-time reference)
        m1 = smf.ols(
            f"{ycol} ~ d_pre_overdue + d_post_overdue + is_industry + is_fda_drug_ind + log_enroll + C(phase_clean)",
            data=sub_df
        ).fit(cov_type="cluster", cov_kwds={"groups": sub_df["lead_sponsor_clean"]})
        # Wald test of d_post_overdue - d_pre_overdue = 0
        diff_b = float(m1.params["d_post_overdue"] - m1.params["d_pre_overdue"])
        cm = m1.cov_params()
        diff_se = float(np.sqrt(max(cm.loc["d_post_overdue", "d_post_overdue"] + cm.loc["d_pre_overdue", "d_pre_overdue"] - 2 * cm.loc["d_post_overdue", "d_pre_overdue"], 0.0)))
        diff_p = float(2 * stats.norm.sf(abs(diff_b / diff_se))) if diff_se > 0 else 1.0

        # 2. Industry-only overdue contrast (post_sweep_overdue vs pre_sweep_overdue & vs on_time)
        sub_ind = sub_df[sub_df["is_industry"] == 1]
        m_ind = smf.ols(
            f"{ycol} ~ d_pre_overdue + d_post_overdue + is_fda_drug_ind + log_enroll + C(phase_clean)",
            data=sub_ind
        ).fit(cov_type="cluster", cov_kwds={"groups": sub_ind["lead_sponsor_clean"]})
        diff_ind_b = float(m_ind.params["d_post_overdue"] - m_ind.params["d_pre_overdue"])
        cm_i = m_ind.cov_params()
        diff_ind_se = float(np.sqrt(max(cm_i.loc["d_post_overdue", "d_post_overdue"] + cm_i.loc["d_pre_overdue", "d_pre_overdue"] - 2 * cm_i.loc["d_post_overdue", "d_pre_overdue"], 0.0)))
        diff_ind_p = float(2 * stats.norm.sf(abs(diff_ind_b / diff_ind_se))) if diff_ind_se > 0 else 1.0

        # 3. Within-Sponsor Fixed Effects (for sponsors with >=2 trials in sample)
        sc = sub_df.groupby("lead_sponsor_clean")[ycol].transform("count")
        sf = sub_df[sc >= 2].copy()
        fe_res = None
        if len(sf) >= 50:
            X_cov = pd.concat([
                sf[["d_pre_overdue", "d_post_overdue", "is_fda_drug_ind", "log_enroll"]].astype(float),
                pd.get_dummies(sf["phase_clean"], drop_first=True, dtype=float)
            ], axis=1)
            y_v = sf[ycol].astype(float)
            g = sf["lead_sponsor_clean"]
            X_dm = X_cov - X_cov.groupby(g).transform("mean")
            y_dm = y_v - y_v.groupby(g).transform("mean")
            # Drop zero-variance columns
            keep_cols = [c for c in X_dm.columns if X_dm[c].std() > 1e-10]
            m_fe = sm.OLS(y_dm, X_dm[keep_cols]).fit(cov_type="cluster", cov_kwds={"groups": g})
            dof_s = np.sqrt((len(sf) - len(keep_cols)) / max(len(sf) - len(keep_cols) - g.nunique(), 1))
            if "d_post_overdue" in keep_cols and "d_pre_overdue" in keep_cols:
                cm_f = m_fe.cov_params()
                fb_post = float(m_fe.params["d_post_overdue"])
                fse_post = float(m_fe.bse["d_post_overdue"]) * dof_s
                fb_pre = float(m_fe.params["d_pre_overdue"])
                fse_pre = float(m_fe.bse["d_pre_overdue"]) * dof_s
                fdiff = fb_post - fb_pre
                fdiff_se = float(np.sqrt(max(cm_f.loc["d_post_overdue", "d_post_overdue"] + cm_f.loc["d_pre_overdue", "d_pre_overdue"] - 2 * cm_f.loc["d_post_overdue", "d_pre_overdue"], 0.0))) * dof_s
                fe_res = {
                    "n": int(len(sf)), "n_sponsors": int(g.nunique()),
                    "pre_overdue_vs_ontime": {"coef": fb_pre, "se": fse_pre, "p": float(2 * stats.norm.sf(abs(fb_pre / fse_pre)))},
                    "post_overdue_vs_ontime": {"coef": fb_post, "se": fse_post, "p": float(2 * stats.norm.sf(abs(fb_post / fse_post)))},
                    "post_vs_pre_overdue": {"coef": fdiff, "se": fdiff_se, "p": float(2 * stats.norm.sf(abs(fdiff / fdiff_se)))},
                }

        out["regressions"][ycol] = {
            "n": int(len(sub_df)),
            "all_pre_overdue_vs_ontime": {"coef": float(m1.params["d_pre_overdue"]), "se": float(m1.bse["d_pre_overdue"]), "p": float(m1.pvalues["d_pre_overdue"])},
            "all_post_overdue_vs_ontime": {"coef": float(m1.params["d_post_overdue"]), "se": float(m1.bse["d_post_overdue"]), "p": float(m1.pvalues["d_post_overdue"])},
            "all_post_vs_pre_overdue": {"coef": diff_b, "se": diff_se, "p": diff_p},
            "industry_post_overdue_vs_ontime": {"coef": float(m_ind.params["d_post_overdue"]), "se": float(m_ind.bse["d_post_overdue"]), "p": float(m_ind.pvalues["d_post_overdue"])},
            "industry_post_vs_pre_overdue": {"coef": diff_ind_b, "se": diff_ind_se, "p": diff_ind_p},
            "within_sponsor_fe": fe_res,
        }
    return out


def run_randomization_inference(panel, ctrl_arms, n_perm=400, seed=20260330):
    """
    Sponsor-cluster randomization inference (exact permutation p-value) for primary DDD estimate:
    Randomly permute treated vs control arm assignment across sponsors (preserving industry share & arm sizes).
    """
    rng = np.random.default_rng(seed)
    sub = panel[panel["arm"].isin(["us_fdareg_nonp1"] + list(ctrl_arms))].copy()
    # Aggregate to sponsor x post level for blazing-fast exact permutation test
    sp = sub.groupby(["lead_sponsor_clean", "is_act", "is_industry", "post"], as_index=False).agg(
        events=("y", "sum"), risk=("y", "count")
    )
    # Observed DDD in risk-weighted rate (exact equivalent to trial-month LPM without cohort FE; also check fast residualized)
    def calc_ddd(df_s):
        g = df_s.groupby(["is_act", "is_industry", "post"])[["events", "risk"]].sum()
        h = (g["events"] / g["risk"]) * 100.0
        did_ind = (h.loc[(1, 1, 1)] - h.loc[(1, 1, 0)]) - (h.loc[(0, 1, 1)] - h.loc[(0, 1, 0)])
        did_aca = (h.loc[(1, 0, 1)] - h.loc[(1, 0, 0)]) - (h.loc[(0, 0, 1)] - h.loc[(0, 0, 0)])
        return float(did_ind - did_aca)

    obs = calc_ddd(sp)
    # Permute is_act within each is_industry stratum across unique sponsors
    spon_map = sp[["lead_sponsor_clean", "is_industry", "is_act"]].drop_duplicates(["lead_sponsor_clean", "is_industry"]).copy()
    perm_vals = []
    for _ in range(n_perm):
        pm = spon_map.copy()
        for iv in (0, 1):
            idx = pm["is_industry"] == iv
            pm.loc[idx, "is_act_p"] = rng.permutation(pm.loc[idx, "is_act"].values)
        m_sp = sp.drop(columns=["is_act"]).merge(pm[["lead_sponsor_clean", "is_industry", "is_act_p"]], on=["lead_sponsor_clean", "is_industry"])
        m_sp = m_sp.rename(columns={"is_act_p": "is_act"})
        perm_vals.append(calc_ddd(m_sp))
    perm_arr = np.array(perm_vals)
    ri_p = float((np.abs(perm_arr) >= abs(obs)).mean())
    return {
        "n_perm": n_perm,
        "obs_unadj_ddd_pp": obs,
        "perm_mean": float(perm_arr.mean()),
        "perm_sd": float(perm_arr.std()),
        "ri_two_sided_p": ri_p,
    }


def main():
    t0 = time.time()
    df = load_and_prepare_universe()
    print(f"Loaded universe: {len(df):,} unique trials across {df['lead_sponsor_clean'].nunique():,} sponsors", flush=True)

    # 1. Universe summary
    univ_summary = {}
    for arm, g in df.groupby("arm"):
        delinq_mar29 = g[(g["due_date"] < "2026-03-30") & (g["sub_date"].isna() | (g["sub_date"] >= "2026-03-30"))]
        univ_summary[arm] = {
            "total_trials": int(len(g)),
            "industry_trials": int(g["is_industry"].sum()),
            "academic_trials": int((1 - g["is_industry"]).sum()),
            "delinquent_mar29_2026": int(len(delinq_mar29)),
            "delinquent_mar29_industry": int(delinq_mar29["is_industry"].sum()),
            "delinquent_mar29_academic": int((1 - delinq_mar29["is_industry"]).sum()),
            "ever_submitted_share": float(g["sub_date"].notna().mean()),
            "fdaaa_violation_flagged": int(g["has_fdaaa_violation"].sum()),
        }

    # 2. Primary Monthly Panel (2025-04..2026-09, event_ym=2026-04)
    print("Building primary monthly risk-set panel (2025-04..2026-09)...", flush=True)
    panel_main = build_monthly_risk_panel(df, "2025-04", "2026-09", event_ym="2026-04")
    # Save cell-level monthly hazard series for figures
    monthly_hazards = (
        panel_main.groupby(["arm", "is_act", "is_industry", "sponsor_group", "ym", "post"], as_index=False)
        .agg(events=("y", "sum"), risk=("y", "count"))
    )
    monthly_hazards["hazard_pct"] = monthly_hazards["events"] / monthly_hazards["risk"] * 100.0
    monthly_hazards.to_parquet(os.path.join(RES_DIR, "b2_monthly_hazards.parquet"), index=False)

    # 3. Primary DiD & DDD Suite across ALL 3 controls + Pooled
    print("Running primary DDD regressions across all control arms...", flush=True)
    ddd_results = {
        "vs_us_notfdareg": run_ddd_suite(panel_main, ["us_notfdareg"], "ACT_FDA vs US Non-FDA-Regulated"),
        "vs_us_p1_only": run_ddd_suite(panel_main, ["us_p1_only"], "ACT_FDA vs US Phase-1-Only Exempt"),
        "vs_nonus_p23": run_ddd_suite(panel_main, ["nonus_p23_nonfda"], "ACT_FDA vs Non-US Phase 2/3 Exempt"),
        "vs_pooled_exempt": run_ddd_suite(panel_main, ["us_notfdareg", "us_p1_only", "nonus_p23_nonfda"], "ACT_FDA vs All 3 Exempt Controls Pooled"),
    }

    # 4. Continuous-time Cox Proportional Hazards (exact March 30, 2026 split!)
    print("Running exact continuous-time stratified Cox PH models...", flush=True)
    cox_results = {
        "vs_us_notfdareg": run_cox_timevarying(df, ["us_notfdareg"]),
        "vs_us_p1_only": run_cox_timevarying(df, ["us_p1_only"]),
        "vs_nonus_p23": run_cox_timevarying(df, ["nonus_p23_nonfda"]),
        "vs_pooled_exempt": run_cox_timevarying(df, ["us_notfdareg", "us_p1_only", "nonus_p23_nonfda"]),
    }

    # 5. Dynamic Event Studies (Monthly + Exact Biweekly anchored at 2026-03-30!)
    print("Running monthly & exact Monday-anchored biweekly event studies...", flush=True)
    es_monthly = {
        "vs_us_notfdareg": run_event_study_monthly(panel_main, ["us_notfdareg"], ref_ym="2026-02"),
        "vs_pooled_exempt": run_event_study_monthly(panel_main, ["us_notfdareg", "us_p1_only", "nonus_p23_nonfda"], ref_ym="2026-02"),
    }
    es_biweekly = run_biweekly_exact_event_study(df, ["us_notfdareg"])
    es_biweekly_pooled = run_biweekly_exact_event_study(df, ["us_notfdareg", "us_p1_only", "nonus_p23_nonfda"])

    # 6. Historical Placebo Falsifications (March 30, 2025 & March 30, 2024!)
    print("Running historical placebo sweep falsifications (2025 & 2024)...", flush=True)
    panel_p25 = build_monthly_risk_panel(df, "2024-04", "2025-09", event_ym="2025-04")
    panel_p24 = build_monthly_risk_panel(df, "2023-04", "2024-09", event_ym="2024-04")
    placebos = {
        "placebo_2025_vs_us_notfdareg": run_ddd_suite(panel_p25, ["us_notfdareg"], "Placebo Mar 30, 2025 (vs US Non-FDA)"),
        "placebo_2025_vs_pooled": run_ddd_suite(panel_p25, ["us_notfdareg", "us_p1_only", "nonus_p23_nonfda"], "Placebo Mar 30, 2025 (vs Pooled)"),
        "placebo_2024_vs_us_notfdareg": run_ddd_suite(panel_p24, ["us_notfdareg"], "Placebo Mar 30, 2024 (vs US Non-FDA)"),
        "placebo_2024_vs_pooled": run_ddd_suite(panel_p24, ["us_notfdareg", "us_p1_only", "nonus_p23_nonfda"], "Placebo Mar 30, 2024 (vs Pooled)"),
    }

    # 7. Sensitivity Checks (Drop March 2026, Aug 31 cutoff, Oct 15 cutoff, US Phase 2/3 ACT subset, Randomization Inference)
    print("Running robustness & randomization inference checks...", flush=True)
    panel_nomar = build_monthly_risk_panel(df, "2025-04", "2026-09", event_ym="2026-04", drop_ym=["2026-03"])
    panel_aug = build_monthly_risk_panel(df, "2025-04", "2026-08", event_ym="2026-04")
    panel_oct = build_monthly_risk_panel(df, "2025-04", "2026-10", event_ym="2026-04")
    df_p23_act = df[(df["arm"] != "us_fdareg_nonp1") | (df["is_p23"] == 1)].copy()
    panel_p23_act = build_monthly_risk_panel(df_p23_act, "2025-04", "2026-09", event_ym="2026-04")
    sensitivity = {
        "drop_mar2026_vs_us_notfdareg": run_ddd_suite(panel_nomar, ["us_notfdareg"], "Drop March 2026 donut"),
        "cutoff_aug2026_vs_us_notfdareg": run_ddd_suite(panel_aug, ["us_notfdareg"], "Strict Aug 31, 2026 cutoff (86d QC buffer)"),
        "cutoff_oct2026_vs_us_notfdareg": run_ddd_suite(panel_oct, ["us_notfdareg"], "Extended Oct 31, 2026 cutoff"),
        "p23_only_act_vs_nonus_p23": run_ddd_suite(panel_p23_act, ["nonus_p23_nonfda"], "Phase 2/3 ACT only vs Non-US Phase 2/3"),
        "ri_vs_us_notfdareg": run_randomization_inference(panel_main, ["us_notfdareg"], n_perm=500),
        "ri_vs_pooled_exempt": run_randomization_inference(panel_main, ["us_notfdareg", "us_p1_only", "nonus_p23_nonfda"], n_perm=500),
    }

    # 8. Mechanism & Heterogeneity
    print("Running mechanism & heterogeneity decompositions...", flush=True)
    mechanisms = run_mechanism_heterogeneity(df, panel_main)

    # 9. Lemons / Adverse-Selection Analysis (H4)
    print("Running H4 Adverse-Selection / Lemons regressions...", flush=True)
    lemons = run_lemons_analysis()

    # 10. Top Responsive Sponsors (concrete institutional case examples!)
    act_post_subs = df[
        (df["arm"] == "us_fdareg_nonp1") &
        (df["sub_date"] >= "2026-03-30") & (df["sub_date"] <= "2026-09-30") &
        (df["sub_date"] > df["due_date"])
    ]
    top_resp = (
        act_post_subs.groupby(["lead_sponsor_clean", "lead_sponsor_class"], as_index=False)
        .agg(post_overdue_subs=("nct_id", "count"), backlog_mar26=("spon_backlog_mar26", "first"), total_fda=("spon_fda_total", "first"))
        .sort_values("post_overdue_subs", ascending=False)
        .head(20)
        .to_dict(orient="records")
    )

    master = {
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "elapsed_sec": round(time.time() - t0, 1),
        "universe_summary": univ_summary,
        "ddd_primary": ddd_results,
        "cox_timevarying": cox_results,
        "event_study_monthly": es_monthly,
        "event_study_biweekly_us_notfdareg": es_biweekly,
        "event_study_biweekly_pooled": es_biweekly_pooled,
        "placebos": placebos,
        "sensitivity": sensitivity,
        "mechanisms": mechanisms,
        "lemons_h4": lemons,
        "top_post_sweep_sponsors": top_resp,
    }
    out_json = os.path.join(RES_DIR, "b2_econometric_results.json")
    with open(out_json, "w") as f:
        json.dump(master, f, indent=2)
    print(f"DONE in {master['elapsed_sec']}s! Saved {out_json}", flush=True)

    # Print crisp executive summary to stdout
    print("\n=== STUDY B-2 HEADLINE ECONOMETRIC RESULTS ===")
    for k, v in ddd_results.items():
        d_lpm = v["lpm_ddd_cohort_fe"]["ddd"]
        d_sfe = v["lpm_ddd_sponsor_fe"]["ddd"]
        d_poi = v["poisson_ddd"]["ddd"]
        ind_p = v["poisson_did_industry"]
        aca_p = v["poisson_did_academic"]
        cx = cox_results[k]["terms"]["ddd"]
        print(
            f"[{k}] LPM_DDD={d_lpm['coef_pp']:+.3f}pp (p={d_lpm['p']:.2e}) | "
            f"SponFE_DDD={d_sfe['coef_pp']:+.3f}pp (p={d_sfe['p']:.2e}) | "
            f"Pois_DDD_IRR={d_poi['irr']:.3f} (p={d_poi['p']:.2e}) | "
            f"Ind_IRR={ind_p['irr']:.3f} (p={ind_p['p']:.2e}) vs Acad_IRR={aca_p['irr']:.3f} (p={aca_p['p']:.2e}) | "
            f"Cox_DDD_HR={cx['hr']:.3f} (p={cx['p']:.2e})"
        )
    print("\n=== PLACEBO FALSIFICATIONS ===")
    for k, v in placebos.items():
        d_lpm = v["lpm_ddd_cohort_fe"]["ddd"]
        d_poi = v["poisson_ddd"]["ddd"]
        print(f"[{k}] LPM_DDD={d_lpm['coef_pp']:+.3f}pp (p={d_lpm['p']:.3f}) | Pois_DDD_IRR={d_poi['irr']:.3f} (p={d_poi['p']:.3f})")
    print("\n=== PRE-TREND WALD TESTS ===")
    print("Monthly vs_us_notfdareg:", es_monthly["vs_us_notfdareg"]["pre_trend_wald_ddd"])
    print("Monthly vs_pooled:", es_monthly["vs_pooled_exempt"]["pre_trend_wald_ddd"])
    print("Biweekly exact vs_us_notfdareg:", es_biweekly["pre_trend_wald_ddd"])
    print("Biweekly exact vs_pooled:", es_biweekly_pooled["pre_trend_wald_ddd"])
    print("\n=== H4 LEMONS SUMMARY ===")
    print(json.dumps(lemons["group_summary"], indent=2))


if __name__ == "__main__":
    main()
