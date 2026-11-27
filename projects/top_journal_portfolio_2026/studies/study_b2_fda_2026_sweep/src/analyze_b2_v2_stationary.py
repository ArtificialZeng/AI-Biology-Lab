#!/opt/anaconda3/bin/python3
"""
Study B-2 v2 Stationary Rolling-Entry Econometric Engine (Management Science / JFE Standard).

Key econometric upgrades over fixed-cohort v1:
  1. Uses full 53,927-trial steady-state rolling-entry universe (PCD 2017-01-18..2025-09-30,
     including ACTIVE_NOT_RECRUITING trials past primary completion date) so fresh 12-month-overdue
     trials enter continuously in EVERY month 2025-04..2026-09 (zero cohort-aging cutoff distortion!).
  2. Controls for exact Overdue-Age Fixed Effects (months_since_due bins: 0-3m, 3-6m, 6-12m, 12-24m,
     24-48m, 48m+) in ALL LPM and Poisson regressions!
  3. Demonstrates & verifies Multiplicative QC-Censoring Invariance: checks NLM QC lag equivalence
     across Industry vs Academic & runs BOTH clean QC-buffer windows (Apr-Jun 2026 >=148d buffer,
     Apr-Jul 2026 >=117d buffer) AND full Apr-Sep 2026 window!
  4. Runs clean Within-ACT Governance DiD (ACT_FDA Industry vs ACT_FDA Academic) + Triple-Difference
     (DDD vs US Non-FDA, vs Non-US P2/3, vs Pooled Clean Exempt) + Cox PH (stratified by overdue-age/cohort).
  5. Direct Within-Firm Compliance Spillover Test on statutorily exempt US Phase-1 trials: compares
     Industry sponsors who ALSO had >=1 overdue FDA-regulated ACT on Mar 29, 2026 ("Sweep-Notified Firms")
     vs Industry sponsors with 0 overdue ACTs ("Unnotified Firms")!
  6. Exact Poisson Monthly Event-Study (with joint pre-trend Wald test) + 2025 & 2024 Placebos +
     Randomization Inference (500 sponsor-level permutations) + full H4 Lemons audit!
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
    s = series.astype(str).str.strip().replace({"None": np.nan, "nan": np.nan, "": np.nan})
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


def load_steady_universe():
    df = pd.read_parquet(os.path.join(CLEAN_DIR, "b2_trial_universe_steady.parquet"))
    df["pcd"] = parse_ctgov_date(df["primary_completion_date"], conservative_end=True)
    df["sub_date"] = parse_ctgov_date(df["results_first_submit_date"], conservative_end=False)
    df["qc_date"] = parse_ctgov_date(df["results_first_submit_qc_date"], conservative_end=False)
    df["post_date"] = parse_ctgov_date(df["results_first_post_date"], conservative_end=False)
    df["due_date"] = df["pcd"] + pd.Timedelta(days=365)
    df["comp_year"] = df["pcd"].dt.year

    df = df[df["pcd"].notna() & (df["pcd"] >= "2017-01-18") & (df["pcd"] <= "2025-09-30")].copy()
    df["is_industry"] = (df["lead_sponsor_class"] == "INDUSTRY").astype(int)
    df["sponsor_group"] = df["lead_sponsor_class"].apply(classify_sponsor_group)
    df["lead_sponsor_clean"] = df["lead_sponsor_name"].fillna("UNKNOWN_SPONSOR").str.strip()
    df["is_act"] = (df["arm"] == "us_fdareg_nonp1").astype(int)
    df["is_terminated"] = (df["overall_status"] == "TERMINATED").astype(int)
    df["is_anr"] = (df["overall_status"] == "ACTIVE_NOT_RECRUITING").astype(int)
    df["is_p23"] = df["phases"].str.contains("PHASE2|PHASE3", na=False).astype(int)
    df["post_lag_days"] = (df["post_date"] - df["sub_date"]).dt.days
    df["qc_lag_days"] = (df["qc_date"] - df["sub_date"]).dt.days

    # Sponsor-level FDA portfolio & March 29, 2026 delinquent backlog (incl ACTIVE_NOT_RECRUITING!)
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
    return df


def build_stationary_monthly_panel(df, start_ym="2025-04", end_ym="2026-07", event_ym="2026-04", drop_ym=None):
    """
    Stationary rolling-entry delinquent trial-month panel with exact overdue-age bins!
    Every month m includes ALL trials whose statutory due_date < m_start and unsubmitted before m_start.
    """
    months = pd.period_range(start_ym, end_ym, freq="M")
    if drop_ym:
        months = [m for m in months if str(m) not in set(drop_ym)]
    records = []
    for m in months:
        m_start = m.start_time
        m_end = m.end_time
        at_risk = (df["due_date"] < m_start) & (df["sub_date"].isna() | (df["sub_date"] >= m_start))
        sub_m = df.loc[at_risk, [
            "nct_id", "arm", "is_act", "is_industry", "sponsor_group", "lead_sponsor_clean",
            "comp_year", "is_p23", "is_terminated", "is_anr", "is_fda_drug", "is_fda_device",
            "spon_fda_total", "spon_backlog_mar26", "due_date", "sub_date"
        ]].copy()
        sub_m["ym"] = str(m)
        sub_m["post"] = int(str(m) >= event_ym)
        sub_m["y"] = ((sub_m["sub_date"] >= m_start) & (sub_m["sub_date"] <= m_end)).astype(int)
        od_m = ((m_start - sub_m["due_date"]).dt.days / 30.4375).clip(lower=0.0)
        sub_m["overdue_months"] = od_m
        sub_m["age_bin"] = pd.cut(
            od_m,
            bins=[-0.01, 3.0, 6.0, 12.0, 24.0, 48.0, 999.0],
            labels=["0_3m", "3_6m", "6_12m", "12_24m", "24_48m", "48m_plus"]
        ).astype(str)
        records.append(sub_m)
    return pd.concat(records, ignore_index=True)


def run_within_act_governance_did(panel, label):
    """
    DIRECT GOVERNANCE TEST INSIDE FDA-REGULATED ACTs (us_fdareg_nonp1):
    Industry vs Academic overdue trials, controlling for exact Overdue-Age x Industry FE + Calendar-Month FE!
    Completely immune to control-arm baseline differences AND multiplicative NLM QC censoring!
    """
    act = panel[panel["arm"] == "us_fdareg_nonp1"].copy()
    act["ind_post"] = act["is_industry"] * act["post"]
    act["y_bp"] = act["y"] * 100.0
    act["cell_ia"] = act["is_industry"].astype(str) + "_" + act["age_bin"]

    # Raw hazards
    g = act.groupby(["is_industry", "post"])["y"].agg(["count", "sum", "mean"]).reset_index()
    haz = {f"ind{int(r.is_industry)}_post{int(r.post)}": {"risk": int(r["count"]), "events": int(r["sum"]), "hazard_pct": round(float(r["mean"]) * 100, 4)} for _, r in g.iterrows()}

    # 1. LPM with Industry x Overdue-Age FE + Month FE
    m_lpm = smf.ols("y_bp ~ ind_post + C(cell_ia) + C(ym)", data=act).fit(
        cov_type="cluster", cov_kwds={"groups": act["lead_sponsor_clean"]}
    )

    # 2. LPM with Sponsor FE + Overdue-Age FE + Month FE
    sc = act.groupby("lead_sponsor_clean")["y"].transform("count")
    af = act[sc >= 2].copy()
    X_s = pd.concat([
        af[["ind_post"]].astype(float),
        pd.get_dummies(af["age_bin"], drop_first=True, dtype=float),
        pd.get_dummies(af["ym"], drop_first=True, dtype=float),
    ], axis=1)
    y_s = af["y_bp"].astype(float)
    grp = af["lead_sponsor_clean"]
    X_dm = X_s - X_s.groupby(grp).transform("mean")
    y_dm = y_s - y_s.groupby(grp).transform("mean")
    m_sfe = sm.OLS(y_dm, X_dm).fit(cov_type="cluster", cov_kwds={"groups": grp})
    dof_s = np.sqrt((len(af) - X_dm.shape[1]) / max(len(af) - X_dm.shape[1] - grp.nunique(), 1))
    b_sfe = float(m_sfe.params["ind_post"])
    se_sfe = float(m_sfe.bse["ind_post"]) * dof_s
    p_sfe = float(2 * stats.norm.sf(abs(b_sfe / se_sfe)))

    # 3. Poisson Exposure-Offset Hazard (EXACTLY invariant to arbitrary month-level QC censoring q(m)!)
    cdf = act.groupby(["lead_sponsor_clean", "is_industry", "age_bin", "ym", "post"], as_index=False).agg(
        events=("y", "sum"), risk=("y", "count")
    )
    cdf["ind_post"] = cdf["is_industry"] * cdf["post"]
    cdf["cell_ia"] = cdf["is_industry"].astype(str) + "_" + cdf["age_bin"]
    m_pois = smf.glm(
        "events ~ ind_post + C(cell_ia) + C(ym)", data=cdf,
        family=sm.families.Poisson(), offset=np.log(cdf["risk"].values)
    ).fit(cov_type="cluster", cov_kwds={"groups": cdf["lead_sponsor_clean"]})
    bp = float(m_pois.params["ind_post"])
    sep = float(m_pois.bse["ind_post"])

    # 4. Poisson with Sponsor-Group 3-way breakdown (INDUSTRY vs ACADEMIC_OTHER vs NIH_FED, ref = ACADEMIC_OTHER)
    cdf3 = act.groupby(["lead_sponsor_clean", "sponsor_group", "age_bin", "ym", "post"], as_index=False).agg(
        events=("y", "sum"), risk=("y", "count")
    )
    cdf3["ind_post"] = ((cdf3["sponsor_group"] == "INDUSTRY") & (cdf3["post"] == 1)).astype(int)
    cdf3["nih_post"] = ((cdf3["sponsor_group"] == "NIH_FED") & (cdf3["post"] == 1)).astype(int)
    cdf3["cell_sa"] = cdf3["sponsor_group"] + "_" + cdf3["age_bin"]
    cdf3 = cdf3[cdf3.groupby("cell_sa")["events"].transform("sum") > 0]
    m_p3 = smf.glm(
        "events ~ ind_post + nih_post + C(cell_sa) + C(ym)", data=cdf3,
        family=sm.families.Poisson(), offset=np.log(cdf3["risk"].values)
    ).fit(cov_type="cluster", cov_kwds={"groups": cdf3["lead_sponsor_clean"]})

    return {
        "label": label,
        "n_obs": int(len(act)),
        "n_trials": int(act["nct_id"].nunique()),
        "n_sponsors": int(act["lead_sponsor_clean"].nunique()),
        "n_events": int(act["y"].sum()),
        "raw_hazards": haz,
        "lpm_age_fe": {
            "coef_pp": float(m_lpm.params["ind_post"]), "se_pp": float(m_lpm.bse["ind_post"]),
            "t": float(m_lpm.tvalues["ind_post"]), "p": float(m_lpm.pvalues["ind_post"]),
            "ci95": [float(x) for x in m_lpm.conf_int().loc["ind_post"]],
        },
        "lpm_sponsor_age_fe": {
            "coef_pp": b_sfe, "se_pp": se_sfe, "z": b_sfe / se_sfe, "p": p_sfe,
            "ci95": [b_sfe - 1.96 * se_sfe, b_sfe + 1.96 * se_sfe],
            "n_sponsors_fe": int(grp.nunique()),
        },
        "poisson_age_fe": {
            "log_irr": bp, "se": sep, "z": float(m_pois.tvalues["ind_post"]), "p": float(m_pois.pvalues["ind_post"]),
            "irr": float(np.exp(bp)), "irr_ci95": [float(np.exp(bp - 1.96 * sep)), float(np.exp(bp + 1.96 * sep))],
        },
        "poisson_3group_vs_academic": {
            "industry_vs_academic": {
                "irr": float(np.exp(m_p3.params["ind_post"])),
                "p": float(m_p3.pvalues["ind_post"]),
                "ci95": [float(np.exp(x)) for x in m_p3.conf_int().loc["ind_post"]],
            },
            "nih_fed_vs_academic": {
                "irr": float(np.exp(m_p3.params["nih_post"])),
                "p": float(m_p3.pvalues["nih_post"]),
                "ci95": [float(np.exp(x)) for x in m_p3.conf_int().loc["nih_post"]],
            },
        },
    }


def run_age_adjusted_ddd(panel, ctrl_arms, label):
    """
    Stationary Overdue-Age-Adjusted DiD & DDD (LPM + Multiplicative Poisson + Sponsor FE).
    """
    sub = panel[panel["arm"].isin(["us_fdareg_nonp1"] + list(ctrl_arms))].copy()
    sub["act_post"] = sub["is_act"] * sub["post"]
    sub["ind_post"] = sub["is_industry"] * sub["post"]
    sub["ddd"] = sub["is_act"] * sub["is_industry"] * sub["post"]
    sub["y_bp"] = sub["y"] * 100.0
    sub["cell_aia"] = sub["arm"].astype(str) + "_" + sub["is_industry"].astype(str) + "_" + sub["age_bin"]

    # LPM DDD with Arm x Industry x Overdue-Age FE + Month FE
    m_lpm = smf.ols("y_bp ~ ddd + act_post + ind_post + C(cell_aia) + C(ym)", data=sub).fit(
        cov_type="cluster", cov_kwds={"groups": sub["lead_sponsor_clean"]}
    )

    # Poisson Exposure-Offset Hazard DDD with Arm x Industry x Overdue-Age FE + Month FE
    cdf = sub.groupby(["lead_sponsor_clean", "arm", "is_act", "is_industry", "age_bin", "ym", "post"], as_index=False).agg(
        events=("y", "sum"), risk=("y", "count")
    )
    cdf["act_post"] = cdf["is_act"] * cdf["post"]
    cdf["ind_post"] = cdf["is_industry"] * cdf["post"]
    cdf["ddd"] = cdf["is_act"] * cdf["is_industry"] * cdf["post"]
    cdf["cell_aia"] = cdf["arm"].astype(str) + "_" + cdf["is_industry"].astype(str) + "_" + cdf["age_bin"]
    cdf_p = cdf[cdf.groupby("cell_aia")["events"].transform("sum") > 0].copy()

    m_pois = smf.glm(
        "events ~ ddd + act_post + ind_post + C(cell_aia) + C(ym)",
        data=cdf_p, family=sm.families.Poisson(), offset=np.log(cdf_p["risk"].values)
    ).fit(cov_type="cluster", cov_kwds={"groups": cdf_p["lead_sponsor_clean"]})

    # Industry-only & Academic-only Poisson DiD (comparing ACT vs Exempt within same sponsor class & overdue age!)
    ind_aca = {}
    for tag, iv in [("industry", 1), ("academic", 0)]:
        s2 = cdf[cdf["is_industry"] == iv].copy()
        s2["cell_aa"] = s2["arm"].astype(str) + "_" + s2["age_bin"]
        s2 = s2[s2.groupby("cell_aa")["events"].transform("sum") > 0]
        mp = smf.glm(
            "events ~ act_post + C(cell_aa) + C(ym)", data=s2,
            family=sm.families.Poisson(), offset=np.log(s2["risk"].values)
        ).fit(cov_type="cluster", cov_kwds={"groups": s2["lead_sponsor_clean"]})
        b = float(mp.params["act_post"])
        se = float(mp.bse["act_post"])
        ind_aca[tag] = {
            "irr": float(np.exp(b)), "log_irr": b, "se": se, "z": float(mp.tvalues["act_post"]),
            "p": float(mp.pvalues["act_post"]),
            "ci95": [float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se))],
        }
        # Also LPM
        s2_l = sub[sub["is_industry"] == iv].copy()
        s2_l["cell_aa"] = s2_l["arm"].astype(str) + "_" + s2_l["age_bin"]
        ml = smf.ols("y_bp ~ act_post + C(cell_aa) + C(ym)", data=s2_l).fit(
            cov_type="cluster", cov_kwds={"groups": s2_l["lead_sponsor_clean"]}
        )
        ind_aca[tag + "_lpm"] = {
            "coef_pp": float(ml.params["act_post"]), "se_pp": float(ml.bse["act_post"]),
            "p": float(ml.pvalues["act_post"]),
        }

    # Cell summary hazards
    g = sub.groupby(["is_act", "is_industry", "post"])["y"].agg(["count", "sum", "mean"]).reset_index()
    haz = [
        {"is_act": int(r.is_act), "is_industry": int(r.is_industry), "post": int(r.post),
         "risk": int(r["count"]), "events": int(r["sum"]), "hazard_pct": round(float(r["mean"]) * 100, 4)}
        for _, r in g.iterrows()
    ]

    return {
        "label": label,
        "n_obs": int(len(sub)),
        "n_trials": int(sub["nct_id"].nunique()),
        "n_events": int(sub["y"].sum()),
        "cell_hazards": haz,
        "lpm_ddd_age_fe": {
            t: {"coef_pp": float(m_lpm.params[t]), "se_pp": float(m_lpm.bse[t]),
                "t": float(m_lpm.tvalues[t]), "p": float(m_lpm.pvalues[t]),
                "ci95": [float(x) for x in m_lpm.conf_int().loc[t]]}
            for t in ("ddd", "act_post", "ind_post")
        },
        "poisson_ddd_age_fe": {
            t: {"log_irr": float(m_pois.params[t]), "se": float(m_pois.bse[t]),
                "z": float(m_pois.tvalues[t]), "p": float(m_pois.pvalues[t]),
                "irr": float(np.exp(m_pois.params[t])),
                "irr_ci95": [float(np.exp(x)) for x in m_pois.conf_int().loc[t]]}
            for t in ("ddd", "act_post", "ind_post")
        },
        "stratified_did": ind_aca,
    }


def run_poisson_event_study_within_act(panel, ref_ym="2026-02"):
    """
    EXACT Multiplicative Poisson Event Study inside ACT_FDA (Industry vs Academic, with Industry x Age FE + Month FE)!
    100% invariant to arbitrary month-level NLM QC censoring q(m)!
    """
    act = panel[panel["arm"] == "us_fdareg_nonp1"].copy()
    cdf = act.groupby(["lead_sponsor_clean", "is_industry", "age_bin", "ym"], as_index=False).agg(
        events=("y", "sum"), risk=("y", "count")
    )
    cdf["cell_ia"] = cdf["is_industry"].astype(str) + "_" + cdf["age_bin"]
    months = sorted(cdf["ym"].unique())
    ref_i = months.index(ref_ym)
    terms = []
    for m in months:
        if m == ref_ym:
            continue
        k = months.index(m) - ref_i
        tag = f"m{abs(k)}" if k < 0 else f"p{k}"
        cdf[f"ip_{tag}"] = ((cdf["is_industry"] == 1) & (cdf["ym"] == m)).astype(float)
        terms.append((m, k, f"ip_{tag}"))

    rhs = " + ".join([t[2] for t in terms]) + " + C(cell_ia) + C(ym)"
    fit = smf.glm(
        f"events ~ {rhs}", data=cdf, family=sm.families.Poisson(), offset=np.log(cdf["risk"].values)
    ).fit(cov_type="cluster", cov_kwds={"groups": cdf["lead_sponsor_clean"]})

    # Also compute exact raw monthly hazards for ACT Industry, ACT Academic, US Non-FDA Academic, Non-US P23 Industry
    series = []
    pre_names = []
    for m in months:
        k = months.index(m) - ref_i
        if m == ref_ym:
            series.append({"ym": m, "rel_m": k, "log_irr": 0.0, "se": 0.0, "irr": 1.0, "irr_ci95": [1.0, 1.0], "p": 1.0})
        else:
            tag = f"m{abs(k)}" if k < 0 else f"p{k}"
            nm = f"ip_{tag}"
            if k < 0:
                pre_names.append(nm)
            b = float(fit.params[nm])
            se = float(fit.bse[nm])
            series.append({
                "ym": m, "rel_m": k, "log_irr": b, "se": se,
                "irr": float(np.exp(b)),
                "irr_ci95": [float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se))],
                "p": float(fit.pvalues[nm]),
            })
    R = np.zeros((len(pre_names), len(fit.params)))
    pnames = list(fit.params.index)
    for i, nm in enumerate(pre_names):
        R[i, pnames.index(nm)] = 1.0
    wt = fit.wald_test(R, scalar=True)

    # Also LPM Event Study inside ACT with Overdue-Age FE
    act["y_bp"] = act["y"] * 100.0
    act["cell_ia"] = act["is_industry"].astype(str) + "_" + act["age_bin"]
    for m, k, nm in terms:
        act[nm] = ((act["is_industry"] == 1) & (act["ym"] == m)).astype(float)
    fit_lpm = smf.ols(f"y_bp ~ {rhs}", data=act).fit(cov_type="cluster", cov_kwds={"groups": act["lead_sponsor_clean"]})
    R2 = np.zeros((len(pre_names), len(fit_lpm.params)))
    pnames2 = list(fit_lpm.params.index)
    for i, nm in enumerate(pre_names):
        R2[i, pnames2.index(nm)] = 1.0
    wt_lpm = fit_lpm.wald_test(R2, scalar=True)
    for row in series:
        m = row["ym"]
        if m == ref_ym:
            row["lpm_pp"], row["lpm_se"], row["lpm_p"] = 0.0, 0.0, 1.0
        else:
            k = row["rel_m"]
            nm = f"ip_{'m' + str(abs(k)) if k < 0 else 'p' + str(k)}"
            row["lpm_pp"] = float(fit_lpm.params[nm])
            row["lpm_se"] = float(fit_lpm.bse[nm])
            row["lpm_p"] = float(fit_lpm.pvalues[nm])

    return {
        "ref_ym": ref_ym,
        "poisson_pre_trend_wald": {"df": len(pre_names), "stat": float(wt.statistic), "p": float(wt.pvalue)},
        "lpm_pre_trend_wald": {"df": len(pre_names), "stat": float(wt_lpm.statistic), "p": float(wt_lpm.pvalue)},
        "series": series,
    }


def run_phase1_spillover_test(df, panel):
    """
    WITHIN-FIRM ORGANIZATIONAL COMPLIANCE SPILLOVER TEST:
    Inside statutorily exempt US Phase-1-Only Industry trials (us_p1_only, is_industry=1),
    compare:
      - Exposed Firms: Industry sponsors with >=1 overdue FDA-regulated ACT on March 29, 2026
        (i.e., received FDA's March 30, 2026 sweep email!)
      - Unexposed Firms: Industry sponsors with 0 overdue ACTs on March 29, 2026
    Test whether the post-sweep US Phase-1 submission lift occurs specifically among Sweep-Notified firms!
    """
    p1 = panel[(panel["arm"] == "us_p1_only") & (panel["is_industry"] == 1)].copy()
    p1["exposed_firm"] = (p1["spon_backlog_mar26"] >= 1).astype(int)
    p1["exp_post"] = p1["exposed_firm"] * p1["post"]
    p1["y_bp"] = p1["y"] * 100.0
    p1["cell_ea"] = p1["exposed_firm"].astype(str) + "_" + p1["age_bin"]

    # Raw hazards pre vs post (and April 2026 specifically!)
    g = p1.groupby(["exposed_firm", "post"])["y"].agg(["count", "sum", "mean"]).reset_index()
    haz = {f"exp{int(r.exposed_firm)}_post{int(r.post)}": {"risk": int(r["count"]), "events": int(r["sum"]), "hazard_pct": round(float(r["mean"]) * 100, 4)} for _, r in g.iterrows()}

    # April 2026 sharp spike test (2026-04 vs pre-sweep 2025-04..2026-03)
    p1_apr = p1[p1["ym"] <= "2026-04"].copy()
    p1_apr["is_apr26"] = (p1_apr["ym"] == "2026-04").astype(int)
    p1_apr["exp_apr"] = p1_apr["exposed_firm"] * p1_apr["is_apr26"]
    g_apr = p1_apr.groupby(["exposed_firm", "is_apr26"])["y"].agg(["count", "sum", "mean"]).reset_index()
    haz_apr = {f"exp{int(r.exposed_firm)}_apr{int(r.is_apr26)}": {"risk": int(r["count"]), "events": int(r["sum"]), "hazard_pct": round(float(r["mean"]) * 100, 4)} for _, r in g_apr.iterrows()}

    # Poisson DiD across full post window
    cdf = p1.groupby(["lead_sponsor_clean", "exposed_firm", "age_bin", "ym", "post"], as_index=False).agg(
        events=("y", "sum"), risk=("y", "count")
    )
    cdf["exp_post"] = cdf["exposed_firm"] * cdf["post"]
    cdf["cell_ea"] = cdf["exposed_firm"].astype(str) + "_" + cdf["age_bin"]
    cdf = cdf[cdf.groupby("cell_ea")["events"].transform("sum") > 0]
    mp = smf.glm("events ~ exp_post + C(cell_ea) + C(ym)", data=cdf,
                 family=sm.families.Poisson(), offset=np.log(cdf["risk"].values)).fit(
        cov_type="cluster", cov_kwds={"groups": cdf["lead_sponsor_clean"]}
    )
    ml = smf.ols("y_bp ~ exp_post + C(cell_ea) + C(ym)", data=p1).fit(
        cov_type="cluster", cov_kwds={"groups": p1["lead_sponsor_clean"]}
    )

    # Also April 2026 Poisson DiD
    cdf_a = p1_apr.groupby(["lead_sponsor_clean", "exposed_firm", "age_bin", "ym", "is_apr26"], as_index=False).agg(
        events=("y", "sum"), risk=("y", "count")
    )
    cdf_a["exp_apr"] = cdf_a["exposed_firm"] * cdf_a["is_apr26"]
    cdf_a["cell_ea"] = cdf_a["exposed_firm"].astype(str) + "_" + cdf_a["age_bin"]
    cdf_a = cdf_a[cdf_a.groupby("cell_ea")["events"].transform("sum") > 0]
    mp_a = smf.glm("events ~ exp_apr + C(cell_ea) + C(ym)", data=cdf_a,
                   family=sm.families.Poisson(), offset=np.log(cdf_a["risk"].values)).fit(
        cov_type="cluster", cov_kwds={"groups": cdf_a["lead_sponsor_clean"]}
    )
    return {
        "n_obs": int(len(p1)),
        "n_trials": int(p1["nct_id"].nunique()),
        "n_exposed_firms": int(p1.loc[p1["exposed_firm"] == 1, "lead_sponsor_clean"].nunique()),
        "n_unexposed_firms": int(p1.loc[p1["exposed_firm"] == 0, "lead_sponsor_clean"].nunique()),
        "hazards_full_post": haz,
        "hazards_apr2026": haz_apr,
        "poisson_full_post": {
            "irr": float(np.exp(mp.params["exp_post"])), "p": float(mp.pvalues["exp_post"]),
            "ci95": [float(np.exp(x)) for x in mp.conf_int().loc["exp_post"]],
        },
        "lpm_full_post": {
            "coef_pp": float(ml.params["exp_post"]), "se_pp": float(ml.bse["exp_post"]), "p": float(ml.pvalues["exp_post"]),
        },
        "poisson_apr2026_spike": {
            "irr": float(np.exp(mp_a.params["exp_apr"])), "p": float(mp_a.pvalues["exp_apr"]),
            "ci95": [float(np.exp(x)) for x in mp_a.conf_int().loc["exp_apr"]],
        },
    }


def run_qc_lag_symmetry_test(df):
    """Verify NLM QC review-to-posting lag symmetry between Industry & Academic submissions (proving multiplicative invariance)."""
    s = df[(df["arm"] == "us_fdareg_nonp1") & df["sub_date"].notna() & df["post_lag_days"].notna()].copy()
    s = s[(s["sub_date"] >= "2024-01-01") & (s["sub_date"] <= "2026-06-30") & (s["post_lag_days"] >= 0)]
    ind = s[s["is_industry"] == 1]["post_lag_days"]
    aca = s[s["is_industry"] == 0]["post_lag_days"]
    _, mw_p = stats.mannwhitneyu(ind, aca)
    # Also check post-sweep specifically (Apr-Jul 2026)
    s_post = df[(df["arm"] == "us_fdareg_nonp1") & (df["sub_date"] >= "2026-04-01") & (df["sub_date"] <= "2026-07-31") & df["post_lag_days"].notna()]
    ind_p = s_post[s_post["is_industry"] == 1]["post_lag_days"]
    aca_p = s_post[s_post["is_industry"] == 0]["post_lag_days"]
    _, mw_pp = stats.mannwhitneyu(ind_p, aca_p) if (len(ind_p) > 5 and len(aca_p) > 5) else (np.nan, np.nan)
    return {
        "overall_2024_2026q2": {
            "industry_n": int(len(ind)), "industry_median_days": float(ind.median()), "industry_iqr": [float(ind.quantile(0.25)), float(ind.quantile(0.75))],
            "academic_n": int(len(aca)), "academic_median_days": float(aca.median()), "academic_iqr": [float(aca.quantile(0.25)), float(aca.quantile(0.75))],
            "mann_whitney_p": float(mw_p),
        },
        "post_sweep_apr_jul_2026": {
            "industry_n": int(len(ind_p)), "industry_median_days": float(ind_p.median()) if len(ind_p) else None,
            "academic_n": int(len(aca_p)), "academic_median_days": float(aca_p.median()) if len(aca_p) else None,
            "mann_whitney_p": float(mw_pp),
        },
    }


def run_exact_biweekly_within_act(df):
    """Exact Monday-anchored biweekly Poisson & LPM event-study inside ACT_FDA (2026-03-30 = biweek 0)."""
    act_df = df[df["arm"] == "us_fdareg_nonp1"].copy()
    t0 = pd.Timestamp("2026-03-30")
    recs = []
    # 13 pre biweeks (-13..-1) and 9 post biweeks (0..+8, ending 2026-08-02 -> >=115d QC buffer!)
    for bw in range(-13, 9):
        b_start = t0 + pd.Timedelta(days=14 * bw)
        b_end = b_start + pd.Timedelta(days=13, hours=23, minutes=59, seconds=59)
        ar = (act_df["due_date"] < b_start) & (act_df["sub_date"].isna() | (act_df["sub_date"] >= b_start))
        sub = act_df.loc[ar, ["nct_id", "is_industry", "lead_sponsor_clean", "due_date", "sub_date"]].copy()
        sub["biweek"] = bw
        sub["y"] = ((sub["sub_date"] >= b_start) & (sub["sub_date"] <= b_end)).astype(int)
        od_m = ((b_start - sub["due_date"]).dt.days / 30.4375).clip(lower=0.0)
        sub["age_bin"] = pd.cut(od_m, bins=[-0.01, 3.0, 6.0, 12.0, 24.0, 48.0, 999.0],
                                labels=["0_3m", "3_6m", "6_12m", "12_24m", "24_48m", "48m_plus"]).astype(str)
        recs.append(sub)
    bp = pd.concat(recs, ignore_index=True)

    # Raw biweekly hazards
    raw = bp.groupby(["biweek", "is_industry"], as_index=False).agg(events=("y", "sum"), risk=("y", "count"))
    raw["hazard_pct"] = raw["events"] / raw["risk"] * 100.0
    raw["start_date"] = raw["biweek"].map(lambda b: (t0 + pd.Timedelta(days=14 * b)).strftime("%Y-%m-%d"))

    # Poisson biweekly event study (ref = -1, Mar 16-29, 2026)
    cdf = bp.groupby(["lead_sponsor_clean", "is_industry", "age_bin", "biweek"], as_index=False).agg(
        events=("y", "sum"), risk=("y", "count")
    )
    cdf["cell_ia"] = cdf["is_industry"].astype(str) + "_" + cdf["age_bin"]
    bws = sorted(cdf["biweek"].unique())
    ref_bw = -1
    terms = []
    for b in bws:
        if b == ref_bw:
            continue
        tag = f"m{abs(b)}" if b < 0 else f"p{b}"
        cdf[f"ip_{tag}"] = ((cdf["is_industry"] == 1) & (cdf["biweek"] == b)).astype(float)
        terms.append((b, f"ip_{tag}"))
    rhs = " + ".join([t[1] for t in terms]) + " + C(cell_ia) + C(biweek)"
    fit = smf.glm(f"events ~ {rhs}", data=cdf, family=sm.families.Poisson(), offset=np.log(cdf["risk"].values)).fit(
        cov_type="cluster", cov_kwds={"groups": cdf["lead_sponsor_clean"]}
    )
    series, pre_names = [], []
    for b in bws:
        sd = (t0 + pd.Timedelta(days=14 * b)).strftime("%Y-%m-%d")
        if b == ref_bw:
            series.append({"biweek": int(b), "start_date": sd, "log_irr": 0.0, "se": 0.0, "irr": 1.0, "irr_ci95": [1.0, 1.0], "p": 1.0})
        else:
            tag = f"m{abs(b)}" if b < 0 else f"p{b}"
            nm = f"ip_{tag}"
            if b < -1:
                pre_names.append(nm)
            val = float(fit.params[nm])
            se = float(fit.bse[nm])
            series.append({
                "biweek": int(b), "start_date": sd, "log_irr": val, "se": se,
                "irr": float(np.exp(val)), "irr_ci95": [float(np.exp(val - 1.96 * se)), float(np.exp(val + 1.96 * se))],
                "p": float(fit.pvalues[nm]),
            })
    R = np.zeros((len(pre_names), len(fit.params)))
    pnames = list(fit.params.index)
    for i, nm in enumerate(pre_names):
        R[i, pnames.index(nm)] = 1.0
    wt = fit.wald_test(R, scalar=True)
    return {
        "pre_trend_wald_poisson": {"df": len(pre_names), "stat": float(wt.statistic), "p": float(wt.pvalue)},
        "series": series,
        "raw_biweekly": raw.to_dict(orient="records"),
    }


def run_ri_within_act(panel, n_perm=500, seed=20260330):
    """Exact Sponsor-Cluster Randomization Inference for Within-ACT Industry vs Academic Governance DiD."""
    rng = np.random.default_rng(seed)
    act = panel[panel["arm"] == "us_fdareg_nonp1"].copy()
    sp = act.groupby(["lead_sponsor_clean", "is_industry", "post"], as_index=False).agg(events=("y", "sum"), risk=("y", "count"))
    def calc_log_irr_did(df_s):
        g = df_s.groupby(["is_industry", "post"])[["events", "risk"]].sum()
        h = g["events"] / g["risk"]
        return float(np.log((h.loc[(1, 1)] / h.loc[(1, 0)]) / (h.loc[(0, 1)] / h.loc[(0, 0)])))
    obs = calc_log_irr_did(sp)
    smap = sp[["lead_sponsor_clean", "is_industry"]].drop_duplicates("lead_sponsor_clean").copy()
    perms = []
    for _ in range(n_perm):
        pm = smap.copy()
        pm["ind_p"] = rng.permutation(pm["is_industry"].values)
        m = sp.drop(columns=["is_industry"]).merge(pm[["lead_sponsor_clean", "ind_p"]], on="lead_sponsor_clean").rename(columns={"ind_p": "is_industry"})
        perms.append(calc_log_irr_did(m))
    pa = np.array(perms)
    return {
        "n_perm": n_perm, "obs_log_irr": obs, "obs_irr": float(np.exp(obs)),
        "perm_mean": float(pa.mean()), "perm_sd": float(pa.std()),
        "ri_two_sided_p": float((np.abs(pa) >= abs(obs)).mean()),
    }


def main():
    t0 = time.time()
    df = load_steady_universe()
    print(f"Loaded steady-state universe: {len(df):,} unique trials across {df['lead_sponsor_clean'].nunique():,} sponsors", flush=True)

    # Universe stock on March 29, 2026 (including ACTIVE_NOT_RECRUITING!)
    stock = {}
    for arm, g in df.groupby("arm"):
        d29 = g[(g["due_date"] < "2026-03-30") & (g["sub_date"].isna() | (g["sub_date"] >= "2026-03-30"))]
        stock[arm] = {
            "total_trials": int(len(g)),
            "delinq_mar29_total": int(len(d29)),
            "delinq_mar29_industry": int(d29["is_industry"].sum()),
            "delinq_mar29_academic": int((1 - d29["is_industry"]).sum()),
            "delinq_mar29_unique_sponsors": int(d29["lead_sponsor_clean"].nunique()),
        }

    # Build primary clean QC-buffered steady-state panel: 2025-04..2026-07 (12m pre + 4m post, ALL >=117d QC buffer!)
    # Plus Q2 ultra-clean (2025-04..2026-06, >=148d buffer) AND full 6m post (2025-04..2026-09) AND donut (drop 2026-03)
    panel_jul = build_stationary_monthly_panel(df, "2025-04", "2026-07", "2026-04")
    panel_jun = build_stationary_monthly_panel(df, "2025-04", "2026-06", "2026-04")
    panel_sep = build_stationary_monthly_panel(df, "2025-04", "2026-09", "2026-04")
    panel_donut = build_stationary_monthly_panel(df, "2025-04", "2026-07", "2026-04", drop_ym=["2026-03"])

    # Historical Placebos (exact same 12m pre + 4m post stationary rolling-entry structure!)
    panel_p25 = build_stationary_monthly_panel(df, "2024-04", "2025-07", "2025-04")
    panel_p24 = build_stationary_monthly_panel(df, "2023-04", "2024-07", "2024-04")

    # Save stationary monthly hazards across full 2025-04..2026-09 for plotting!
    mh = (
        panel_sep.groupby(["arm", "is_act", "is_industry", "sponsor_group", "ym", "post"], as_index=False)
        .agg(events=("y", "sum"), risk=("y", "count"))
    )
    mh["hazard_pct"] = mh["events"] / mh["risk"] * 100.0
    mh.to_parquet(os.path.join(RES_DIR, "b2_v2_stationary_monthly_hazards.parquet"), index=False)

    print("1. Running Within-ACT Governance DiD (2026 Primary, Q2, Full 6m, Donut, Placebos 2025 & 2024)...", flush=True)
    within_act = {
        "primary_apr_jul_2026": run_within_act_governance_did(panel_jul, "Primary Apr-Jul 2026 (>=117d QC buffer)"),
        "ultraclean_apr_jun_2026": run_within_act_governance_did(panel_jun, "Q2 Apr-Jun 2026 (>=148d QC buffer)"),
        "full_apr_sep_2026": run_within_act_governance_did(panel_sep, "Full 6m Apr-Sep 2026"),
        "donut_drop_mar2026": run_within_act_governance_did(panel_donut, "Donut dropping March 2026"),
        "placebo_2025_apr_jul": run_within_act_governance_did(panel_p25, "Placebo March 30, 2025 (Apr-Jul 2025)"),
        "placebo_2024_apr_jul": run_within_act_governance_did(panel_p24, "Placebo March 30, 2024 (Apr-Jul 2024)"),
    }

    print("2. Running Stationary Overdue-Age-Adjusted DDD vs All Exempt Control Arms...", flush=True)
    ddd_v2 = {
        "vs_us_notfdareg_jul": run_age_adjusted_ddd(panel_jul, ["us_notfdareg"], "ACT vs US Non-FDA (Apr-Jul 2026)"),
        "vs_nonus_p23_jul": run_age_adjusted_ddd(panel_jul, ["nonus_p23_nonfda"], "ACT vs Non-US P2/3 (Apr-Jul 2026)"),
        "vs_clean_exempt_jul": run_age_adjusted_ddd(panel_jul, ["us_notfdareg", "nonus_p23_nonfda"], "ACT vs Clean Exempt (US Non-FDA + Non-US P2/3, Apr-Jul 2026)"),
        "vs_us_p1_jul": run_age_adjusted_ddd(panel_jul, ["us_p1_only"], "ACT vs US Phase-1 (Spillover Arm, Apr-Jul 2026)"),
        "vs_clean_exempt_sep": run_age_adjusted_ddd(panel_sep, ["us_notfdareg", "nonus_p23_nonfda"], "ACT vs Clean Exempt (Full Apr-Sep 2026)"),
        "placebo_2025_vs_clean_exempt": run_age_adjusted_ddd(panel_p25, ["us_notfdareg", "nonus_p23_nonfda"], "Placebo 2025 vs Clean Exempt"),
        "placebo_2024_vs_clean_exempt": run_age_adjusted_ddd(panel_p24, ["us_notfdareg", "nonus_p23_nonfda"], "Placebo 2024 vs Clean Exempt"),
    }

    print("3. Running Exact Multiplicative Poisson Event Studies (Monthly & Monday-Anchored Biweekly)...", flush=True)
    es_mon = run_poisson_event_study_within_act(panel_sep, ref_ym="2026-02")
    es_bw = run_exact_biweekly_within_act(df)

    print("4. Running Within-Firm Phase-1 Compliance Spillover Test & NLM QC Lag Symmetry Test...", flush=True)
    spillover = run_phase1_spillover_test(df, panel_jul)
    qc_sym = run_qc_lag_symmetry_test(df)
    ri = run_ri_within_act(panel_jul, n_perm=500)

    # 5. Heterogeneity inside stationary panel (Portfolio Size, Backlog, Cohort Depth, Modality, Completed vs ANR)
    print("5. Running Stationary Heterogeneity & Mechanism Decompositions...", flush=True)
    act_j = panel_jul[panel_jul["arm"] == "us_fdareg_nonp1"].copy()
    ctrl_j = panel_jul[panel_jul["arm"].isin(["us_notfdareg", "nonus_p23_nonfda"])].copy()
    def strat_pois(a_sub, c_sub):
        comb = pd.concat([a_sub, c_sub], ignore_index=True)
        cdf = comb.groupby(["lead_sponsor_clean", "arm", "is_act", "age_bin", "ym", "post"], as_index=False).agg(
            events=("y", "sum"), risk=("y", "count")
        )
        cdf["act_post"] = cdf["is_act"] * cdf["post"]
        cdf["cell_aa"] = cdf["arm"].astype(str) + "_" + cdf["age_bin"]
        cdf = cdf[cdf.groupby("cell_aa")["events"].transform("sum") > 0]
        if cdf["events"].sum() < 5 or cdf["is_act"].nunique() < 2:
            return None
        mp = smf.glm("events ~ act_post + C(cell_aa) + C(ym)", data=cdf,
                     family=sm.families.Poisson(), offset=np.log(cdf["risk"].values)).fit(
            cov_type="cluster", cov_kwds={"groups": cdf["lead_sponsor_clean"]}
        )
        b = float(mp.params["act_post"])
        se = float(mp.bse["act_post"])
        return {
            "n_act_events": int(a_sub["y"].sum()),
            "pre_h": round(float(a_sub.loc[a_sub["post"] == 0, "y"].mean()) * 100, 3),
            "post_h": round(float(a_sub.loc[a_sub["post"] == 1, "y"].mean()) * 100, 3),
            "raw_ratio": round(float(a_sub.loc[a_sub["post"] == 1, "y"].mean() / max(a_sub.loc[a_sub["post"] == 0, "y"].mean(), 1e-9)), 3),
            "irr": float(np.exp(b)), "ci95": [float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se))],
            "p": float(mp.pvalues["act_post"]),
        }

    het = {"by_backlog": {}, "by_portfolio": {}, "by_cohort": {}, "by_status": {}}
    for ind_tag, iv in [("industry", 1), ("academic", 0)]:
        het["by_backlog"][ind_tag] = {}
        for lo, hi, nm in [(1, 1, "backlog_1"), (2, 4, "backlog_2_4"), (5, 9999, "backlog_5plus")]:
            het["by_backlog"][ind_tag][nm] = strat_pois(
                act_j[(act_j["is_industry"] == iv) & (act_j["spon_backlog_mar26"] >= lo) & (act_j["spon_backlog_mar26"] <= hi)],
                ctrl_j
            )
        het["by_portfolio"][ind_tag] = {}
        for lo, hi, nm in [(1, 2, "small_1_2"), (3, 19, "mid_3_19"), (20, 9999, "large_20plus")]:
            het["by_portfolio"][ind_tag][nm] = strat_pois(
                act_j[(act_j["is_industry"] == iv) & (act_j["spon_fda_total"] >= lo) & (act_j["spon_fda_total"] <= hi)],
                ctrl_j
            )
        het["by_cohort"][ind_tag] = {}
        for lo, hi, nm in [(2023, 2025, "recent_2023_2025"), (2021, 2022, "mid_2021_2022"), (2017, 2020, "deep_2017_2020")]:
            het["by_cohort"][ind_tag][nm] = strat_pois(
                act_j[(act_j["is_industry"] == iv) & (act_j["comp_year"] >= lo) & (act_j["comp_year"] <= hi)],
                ctrl_j
            )
        het["by_status"][ind_tag] = {
            "completed_or_terminated": strat_pois(act_j[(act_j["is_industry"] == iv) & (act_j["is_anr"] == 0)], ctrl_j),
            "active_not_recruiting": strat_pois(act_j[(act_j["is_industry"] == iv) & (act_j["is_anr"] == 1)], ctrl_j),
        }

    # Include v1 Lemons H4 results directly for unified reference
    v1 = json.load(open(os.path.join(RES_DIR, "b2_econometric_results.json")))
    res = {
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "elapsed_sec": round(time.time() - t0, 1),
        "steady_universe_total": int(len(df)),
        "steady_universe_sponsors": int(df["lead_sponsor_clean"].nunique()),
        "stock_mar29_2026": stock,
        "within_act_governance_did": within_act,
        "ddd_stationary_age_adjusted": ddd_v2,
        "event_study_monthly_within_act": es_mon,
        "event_study_biweekly_within_act": es_bw,
        "phase1_within_firm_spillover": spillover,
        "qc_lag_symmetry": qc_sym,
        "randomization_inference_within_act": ri,
        "heterogeneity_stationary": het,
        "lemons_h4": v1["lemons_h4"],
        "top_post_sweep_sponsors": v1["top_post_sweep_sponsors"],
    }
    out_p = os.path.join(RES_DIR, "b2_v2_stationary_results.json")
    with open(out_p, "w") as f:
        json.dump(res, f, indent=2)

    print(f"\nDONE in {res['elapsed_sec']}s -> {out_p}")
    print("\n=== 1. WITHIN-ACT GOVERNANCE DiD (Industry vs Academic, Overdue-Age FE + Month FE) ===")
    for k, v in within_act.items():
        l = v["lpm_age_fe"]
        sf = v["lpm_sponsor_age_fe"]
        po = v["poisson_age_fe"]
        print(f"  [{k:26s}] LPM={l['coef_pp']:+.3f}pp (p={l['p']:.2e}) | SponFE={sf['coef_pp']:+.3f}pp (p={sf['p']:.2e}) | Pois_IRR={po['irr']:.3f} [{po['irr_ci95'][0]:.2f},{po['irr_ci95'][1]:.2f}] (p={po['p']:.2e})")
    print("\n=== 2. STATIONARY AGE-ADJUSTED DDD (ACT vs Exempt Controls) ===")
    for k, v in ddd_v2.items():
        ld = v["lpm_ddd_age_fe"]["ddd"]
        pd_ = v["poisson_ddd_age_fe"]["ddd"]
        ii = v["stratified_did"]["industry"]
        aa = v["stratified_did"]["academic"]
        print(f"  [{k:30s}] LPM_DDD={ld['coef_pp']:+.3f}pp (p={ld['p']:.2e}) | Pois_DDD_IRR={pd_['irr']:.3f} (p={pd_['p']:.2e}) | Ind_IRR={ii['irr']:.3f} (p={ii['p']:.2e}) vs Acad_IRR={aa['irr']:.3f} (p={aa['p']:.2e})")
    print("\n=== 3. PRE-TREND WALD TESTS (EXACT POISSON + LPM) ===")
    print("  Monthly Poisson pre-trend:", es_mon["poisson_pre_trend_wald"])
    print("  Monthly LPM pre-trend:    ", es_mon["lpm_pre_trend_wald"])
    print("  Biweekly Poisson pre-trend:", es_bw["pre_trend_wald_poisson"])
    print("\n=== 4. PHASE-1 WITHIN-FIRM SPILLOVER ===")
    print(json.dumps(spillover, indent=2))
    print("\n=== 5. QC LAG SYMMETRY & RANDOMIZATION INFERENCE ===")
    print("QC symmetry:", json.dumps(qc_sym, indent=2))
    print("RI within-ACT:", json.dumps(ri, indent=2))


if __name__ == "__main__":
    main()
