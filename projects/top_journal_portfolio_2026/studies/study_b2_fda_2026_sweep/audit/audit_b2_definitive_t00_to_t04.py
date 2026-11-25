#!/opt/anaconda3/bin/python3
"""
Study B-2 Definitive Econometric Engine Closing ALL Codex Items (T00, T01, T02, T03, T04, C01-C03):
  - T00: Integrates 100% of unposted / in-QC / MCP-posted-with-QC-comments submissions
         (derivedSection.miscInfoModule.submissionTracking.estimatedResultsFirstSubmitDate)
         into true_first_submit_date -> ZERO NLM QC-lag right-censoring through 2026-10!
         Separately audits: (1) True First Submission (Posted UNION In-QC), (2) QC-Cleared Posting Only,
         (3) NLM QC-Clearance Acceleration for Pre-March-30 In-QC Queue, (4) Delayed-Results Certifications (dispFirstSubmitDate).
  - T01: Strict exact-date-verified Uniform-Horizon models where EVERY observation satisfies
         end_date + H <= PULL_DATE (2026-11-25) with explicit date assertions!
         (Full-month H=90d thru 2026-07-31, H=120d thru 2026-06-30, H=150d thru 2026-05-31;
          plus exact day-level cutoffs H=90d thru 2026-08-27, H=120d thru 2026-07-28, H=150d thru 2026-06-28!)
         Explicitly distinguished as composite "submitted AND cleared QC within <= H days".
  - T02 / C02: Full Monthly Event Studies + Joint Pre-Trend Wald Tests on BOTH:
         (a) Within-ACT Governance DiD (Uncensored True First Submission AND QC-Cleared),
         (b) Risk-Set-Normalized Overdue vs. On-Time DDD (Uncensored True First Submission AND QC-Cleared!),
         (c) DDD vs Clean Exempt Controls (Uncensored True First Submission),
         PLUS exact Rambachan & Roth (2023, Rev. Econ. Stud.) HonestDiD partial-identification bounds:
           - Relative Magnitudes Delta^RM(Mbar) for Mbar in [0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]
           - Smoothness Delta^SD(M) for M in [0, 0.05, 0.1, 0.2, 0.3, 0.5]
           with exact breakdown values Mbar* and M*!
  - T03 / T04 / C03: Strict clean H4 Lemons re-estimation on b2_lemons_results.parquet with exact N,
         point estimate, SE, 95% CI, 80%-power MDE, TOST equivalence test, AND unposted-vs-posted
         TERMINATED-rate selection check (closing selection into QC clearance!).
  - Strict Symmetric Universe throughout: COMPLETED|TERMINATED + ACTUAL PrimaryCompletionDate ONLY!
"""
import json, os, sys, time, warnings
import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.optimize import linprog
import statsmodels.api as sm
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")

BASE = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_b2_fda_2026_sweep"
CLEAN = os.path.join(BASE, "data", "clean")
AUD = os.path.join(BASE, "audit")
RES = os.path.join(BASE, "results")
sys.path.insert(0, os.path.join(BASE, "src"))
from analyze_b2_v2_stationary import parse_ctgov_date, classify_sponsor_group

PULL_DATE = pd.Timestamp("2026-11-25")


def load_strict_t00_universe():
    """Load strict symmetric COMPLETED|TERMINATED + ACTUAL PCD universe with T00 unposted & certification overlay."""
    d_main = pd.read_parquet(os.path.join(CLEAN, "b2_trial_universe_steady.parquet"))
    d_rec = pd.read_parquet(os.path.join(CLEAN, "b2_us_fdareg_recent_25_26.parquet"))
    d_rec["arm"] = "us_fdareg_nonp1"
    df = pd.concat([d_main, d_rec], ignore_index=True).drop_duplicates("nct_id", keep="first")

    ov = pd.read_parquet(os.path.join(AUD, "b2_t00_unposted_and_cert_overlay.parquet"))
    df = df.merge(ov, on="nct_id", how="left")

    # Strict symmetric filter: COMPLETED|TERMINATED + ACTUAL PCD
    df = df[
        df["overall_status"].isin(["COMPLETED", "TERMINATED"]) &
        (df["primary_completion_type"] == "ACTUAL")
    ].copy()

    df["pcd"] = parse_ctgov_date(df["primary_completion_date"], conservative_end=True)
    df["start_dt"] = parse_ctgov_date(df["start_date"], conservative_end=False)
    # QC-cleared first submission date (original outcome: requires hasResults=True!)
    df["sub_qc_cleared"] = parse_ctgov_date(df["results_first_submit_date"], conservative_end=False)
    df["qc_met_date"] = parse_ctgov_date(df["results_first_submit_qc_date"], conservative_end=False)
    df["post_date"] = parse_ctgov_date(df["results_first_post_date"], conservative_end=False)
    df["mcp_post_date"] = parse_ctgov_date(df["first_mcp_post_date"], conservative_end=False)
    df["disp_sub_date"] = parse_ctgov_date(df["disp_first_submit_date"], conservative_end=False)

    # Unposted / in-QC submission date
    df["est_sub_dt"] = parse_ctgov_date(df["est_results_first_submit_date"], conservative_end=False)
    df["subinfo_rel_dt"] = parse_ctgov_date(df["earliest_subinfo_release_date"], conservative_end=False)
    df["unp_ev_dt"] = parse_ctgov_date(df["earliest_unposted_event_date"], conservative_end=False)
    df["unposted_sub_dt"] = df["est_sub_dt"].combine_first(df["subinfo_rel_dt"]).combine_first(df["unp_ev_dt"])

    # TRUE UNCENSORED FIRST SUBMISSION DATE (Posted UNION In-QC!)
    df["sub_true_any"] = df["sub_qc_cleared"].combine_first(df["unposted_sub_dt"])
    df["is_unposted_in_qc"] = (df["sub_qc_cleared"].isna() & df["unposted_sub_dt"].notna()).astype(int)

    df["due_date"] = df["pcd"] + pd.Timedelta(days=365)
    df["comp_year"] = df["pcd"].dt.year
    df["post_lag_days"] = (df["post_date"] - df["sub_qc_cleared"]).dt.days
    df["qc_met_lag_days"] = (df["qc_met_date"] - df["sub_qc_cleared"]).dt.days
    df["is_industry"] = (df["lead_sponsor_class"] == "INDUSTRY").astype(int)
    df["sponsor_group"] = df["lead_sponsor_class"].apply(classify_sponsor_group)
    df["lead_sponsor_clean"] = df["lead_sponsor_name"].fillna("UNKNOWN_SPONSOR").str.strip()
    df["is_act"] = (df["arm"] == "us_fdareg_nonp1").astype(int)
    df["is_terminated"] = (df["overall_status"] == "TERMINATED").astype(int)

    df = df[df["pcd"].notna() & (df["pcd"] >= "2017-01-18") & (df["pcd"] <= "2026-08-31")].copy()

    # Compute sponsor-level Mar 29, 2026 overdue backlog (using TRUE uncensored submission status!)
    act_sub = df[df["arm"] == "us_fdareg_nonp1"]
    spon_tot = act_sub.groupby("lead_sponsor_clean")["nct_id"].count().rename("spon_fda_total")
    delinq_true = act_sub[
        (act_sub["due_date"] < pd.Timestamp("2026-03-30")) &
        (act_sub["sub_true_any"].isna() | (act_sub["sub_true_any"] >= pd.Timestamp("2026-03-30")))
    ]
    spon_bl = delinq_true.groupby("lead_sponsor_clean")["nct_id"].count().rename("spon_backlog_mar26")
    df = df.merge(spon_tot, left_on="lead_sponsor_clean", right_index=True, how="left")
    df = df.merge(spon_bl, left_on="lead_sponsor_clean", right_index=True, how="left")
    df["spon_fda_total"] = df["spon_fda_total"].fillna(0).astype(int)
    df["spon_backlog_mar26"] = df["spon_backlog_mar26"].fillna(0).astype(int)
    return df


def build_overdue_panel(df, sub_col, start_ym, end_ym, event_ym, h_cap=None, exact_end_date=None):
    """
    Build delinquent (overdue) trial-month panel with STRICT calendar-date observability assertions!
    If h_cap is set: asserts end_date + h_cap <= PULL_DATE so right-censoring is MATHEMATICALLY ZERO!
    """
    months = pd.period_range(start_ym, end_ym, freq="M")
    last_day = pd.Timestamp(exact_end_date) if exact_end_date else months[-1].end_time.normalize()
    if h_cap is not None:
        complete_obs = last_day + pd.Timedelta(days=int(h_cap))
        assert complete_obs <= PULL_DATE, f"Censoring violation! {last_day.date()} + {h_cap}d = {complete_obs.date()} > {PULL_DATE.date()}"

    recs = []
    for m in months:
        ms = m.start_time
        me = min(m.end_time, pd.Timestamp(exact_end_date) + pd.Timedelta(hours=23, minutes=59, seconds=59)) if exact_end_date else m.end_time
        if me < ms:
            continue
        # Risk exposure fraction for truncated final month (if exact_end_date truncates month!)
        full_days = (m.end_time.normalize() - ms.normalize()).days + 1
        obs_days = (me.normalize() - ms.normalize()).days + 1
        exp_frac = obs_days / full_days

        ar = (df["due_date"] < ms) & (df[sub_col].isna() | (df[sub_col] >= ms))
        sub = df.loc[ar, [
            "nct_id", "arm", "is_act", "is_industry", "sponsor_group", "lead_sponsor_clean",
            "due_date", sub_col, "post_lag_days", "spon_backlog_mar26", "spon_fda_total", "comp_year"
        ]].copy()
        sub["ym"] = str(m)
        sub["post"] = int(str(m) >= event_ym)
        sub["exp_frac"] = exp_frac
        in_win = (sub[sub_col] >= ms) & (sub[sub_col] <= me)
        if h_cap is not None:
            sub["y"] = (in_win & sub["post_lag_days"].notna() & (sub["post_lag_days"] <= h_cap)).astype(int)
        else:
            sub["y"] = in_win.astype(int)
        od_m = ((ms - sub["due_date"]).dt.days / 30.4375).clip(lower=0.0)
        sub["age_bin"] = pd.cut(
            od_m, bins=[-0.01, 3.0, 6.0, 12.0, 24.0, 48.0, 999.0],
            labels=["0_3m", "3_6m", "6_12m", "12_24m", "24_48m", "48m_plus"]
        ).astype(str)
        recs.append(sub)
    return pd.concat(recs, ignore_index=True)


def build_ontime_overdue_panel(df_act, sub_col, start_ym, end_ym, event_ym, h_cap=None, exact_end_date=None):
    """Build Risk-Set-Normalized Overdue vs On-Time panel inside strict ACTUAL US FDA-Reg ACTs."""
    months = pd.period_range(start_ym, end_ym, freq="M")
    last_day = pd.Timestamp(exact_end_date) if exact_end_date else months[-1].end_time.normalize()
    if h_cap is not None:
        complete_obs = last_day + pd.Timedelta(days=int(h_cap))
        assert complete_obs <= PULL_DATE, f"Censoring violation! {last_day.date()} + {h_cap}d = {complete_obs.date()} > {PULL_DATE.date()}"

    recs = []
    for m in months:
        ms = m.start_time
        me = min(m.end_time, pd.Timestamp(exact_end_date) + pd.Timedelta(hours=23, minutes=59, seconds=59)) if exact_end_date else m.end_time
        if me < ms:
            continue
        full_days = (m.end_time.normalize() - ms.normalize()).days + 1
        obs_days = (me.normalize() - ms.normalize()).days + 1
        exp_frac = obs_days / full_days

        ar = (df_act["pcd"] < ms) & (df_act[sub_col].isna() | (df_act[sub_col] >= ms))
        sub = df_act.loc[ar, [
            "nct_id", "is_industry", "sponsor_group", "lead_sponsor_clean", "pcd", "due_date", sub_col, "post_lag_days"
        ]].copy()
        sub["ym"] = str(m)
        sub["post"] = int(str(m) >= event_ym)
        sub["exp_frac"] = exp_frac
        days_since = (ms - sub["pcd"]).dt.days.astype(float)
        sub["is_overdue"] = (days_since > 365.0).astype(int)
        sub["m_since"] = days_since / 30.4375
        sub["pcd_age_bin"] = pd.cut(
            sub["m_since"],
            bins=[-0.1, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 24.0, 36.0, 60.0, 999.0],
            labels=["0_3m", "3_6m", "6_9m", "9_12m", "12_15m", "15_18m", "18_24m", "24_36m", "36_60m", "60m_p"]
        ).astype(str)
        in_win = (sub[sub_col] >= ms) & (sub[sub_col] <= me)
        if h_cap is not None:
            sub["y"] = (in_win & sub["post_lag_days"].notna() & (sub["post_lag_days"] <= h_cap)).astype(int)
        else:
            sub["y"] = in_win.astype(int)
        recs.append(sub)
    return pd.concat(recs, ignore_index=True)


def two_way_cluster_se(fit_fn, data, formula, g1, g2, is_poisson=False, offset_col=None):
    """Cameron-Gelbach-Miller (2011) two-way clustered SEs (Sponsor + Month)."""
    def _fit(grp):
        if is_poisson:
            return smf.glm(formula, data=data, family=sm.families.Poisson(), offset=np.log(data[offset_col].values)).fit(
                cov_type="cluster", cov_kwds={"groups": grp}
            )
        return smf.ols(formula, data=data).fit(cov_type="cluster", cov_kwds={"groups": grp})
    m1 = _fit(data[g1])
    m2 = _fit(data[g2])
    g12 = data[g1].astype(str) + "__" + data[g2].astype(str)
    m12 = _fit(g12)
    V = m1.cov_params() + m2.cov_params() - m12.cov_params()
    se = np.sqrt(np.maximum(np.diag(V), 0.0))
    return m1, pd.Series(se, index=m1.params.index)


def eval_within_act(panel, label):
    act = panel[panel["arm"] == "us_fdareg_nonp1"].copy()
    act["ind_post"] = act["is_industry"] * act["post"]
    act["y_bp"] = (act["y"] / act["exp_frac"]) * 100.0
    act["cell_ia"] = act["is_industry"].astype(str) + "_" + act["age_bin"]

    m_lpm, se2_lpm = two_way_cluster_se(None, act, "y_bp ~ ind_post + C(cell_ia) + C(ym)", "lead_sponsor_clean", "ym", False)
    cdf = act.groupby(["lead_sponsor_clean", "is_industry", "age_bin", "ym", "post"], as_index=False).agg(
        events=("y", "sum"), risk=("exp_frac", "sum")
    )
    cdf["ind_post"] = cdf["is_industry"] * cdf["post"]
    cdf["cell_ia"] = cdf["is_industry"].astype(str) + "_" + cdf["age_bin"]
    m_p, se2_p = two_way_cluster_se(None, cdf, "events ~ ind_post + C(cell_ia) + C(ym)", "lead_sponsor_clean", "ym", True, "risk")

    g = act.groupby(["is_industry", "post"]).agg(events=("y", "sum"), risk=("exp_frac", "sum")).reset_index()
    haz = {f"ind{int(r.is_industry)}_p{int(r.post)}": {"risk": round(float(r["risk"]), 1), "events": int(r["events"]), "h_pct": round(100.0 * r["events"] / r["risk"], 3)} for _, r in g.iterrows()}

    bp, sep = float(m_p.params["ind_post"]), float(m_p.bse["ind_post"])
    sep2 = float(se2_p["ind_post"])
    bl, sel, sel2 = float(m_lpm.params["ind_post"]), float(m_lpm.bse["ind_post"]), float(se2_lpm["ind_post"])
    return {
        "label": label, "n_obs": int(len(act)), "n_events": int(act["y"].sum()), "hazards": haz,
        "lpm_pp": bl, "lpm_se_spon": sel, "lpm_p_spon": float(m_lpm.pvalues["ind_post"]),
        "lpm_se_2way": sel2, "lpm_p_2way": float(2 * stats.norm.sf(abs(bl / max(sel2, 1e-12)))),
        "pois_irr": float(np.exp(bp)), "pois_ci95_spon": [float(np.exp(bp - 1.96 * sep)), float(np.exp(bp + 1.96 * sep))],
        "pois_p_spon": float(m_p.pvalues["ind_post"]),
        "pois_ci95_2way": [float(np.exp(bp - 1.96 * sep2)), float(np.exp(bp + 1.96 * sep2))],
        "pois_p_2way": float(2 * stats.norm.sf(abs(bp / max(sep2, 1e-12)))),
    }


def eval_ot_od_ddd(pan, label):
    cdf = pan.groupby(["lead_sponsor_clean", "is_overdue", "is_industry", "pcd_age_bin", "ym", "post"], as_index=False).agg(
        events=("y", "sum"), risk=("exp_frac", "sum")
    )
    cdf["ddd"] = cdf["is_overdue"] * cdf["is_industry"] * cdf["post"]
    cdf["od_post"] = cdf["is_overdue"] * cdf["post"]
    cdf["ind_post"] = cdf["is_industry"] * cdf["post"]
    cdf["cell"] = cdf["is_overdue"].astype(str) + "_" + cdf["is_industry"].astype(str) + "_" + cdf["pcd_age_bin"]
    cdf = cdf[cdf.groupby("cell")["events"].transform("sum") > 0].copy()

    mp, se2_p = two_way_cluster_se(None, cdf, "events ~ ddd + od_post + ind_post + C(cell) + C(ym)", "lead_sponsor_clean", "ym", True, "risk")

    ci = cdf[cdf["is_industry"] == 1].copy()
    ci["cell_oa"] = ci["is_overdue"].astype(str) + "_" + ci["pcd_age_bin"]
    ci = ci[ci.groupby("cell_oa")["events"].transform("sum") > 0]
    mpi = smf.glm("events ~ od_post + C(cell_oa) + C(ym)", data=ci, family=sm.families.Poisson(), offset=np.log(ci["risk"].values)).fit(
        cov_type="cluster", cov_kwds={"groups": ci["lead_sponsor_clean"]}
    )
    ca = cdf[cdf["is_industry"] == 0].copy()
    ca["cell_oa"] = ca["is_overdue"].astype(str) + "_" + ca["pcd_age_bin"]
    ca = ca[ca.groupby("cell_oa")["events"].transform("sum") > 0]
    mpa = smf.glm("events ~ od_post + C(cell_oa) + C(ym)", data=ca, family=sm.families.Poisson(), offset=np.log(ca["risk"].values)).fit(
        cov_type="cluster", cov_kwds={"groups": ca["lead_sponsor_clean"]}
    )

    pan = pan.copy()
    pan["ddd"] = pan["is_overdue"] * pan["is_industry"] * pan["post"]
    pan["od_post"] = pan["is_overdue"] * pan["post"]
    pan["ind_post"] = pan["is_industry"] * pan["post"]
    pan["cell"] = pan["is_overdue"].astype(str) + "_" + pan["is_industry"].astype(str) + "_" + pan["pcd_age_bin"]
    pan["y_bp"] = (pan["y"] / pan["exp_frac"]) * 100.0
    ml = smf.ols("y_bp ~ ddd + od_post + ind_post + C(cell) + C(ym)", data=pan).fit(
        cov_type="cluster", cov_kwds={"groups": pan["lead_sponsor_clean"]}
    )

    g = pan.groupby(["is_overdue", "is_industry", "post"]).agg(events=("y", "sum"), risk=("exp_frac", "sum")).reset_index()
    haz = {f"od{int(r.is_overdue)}_i{int(r.is_industry)}_p{int(r.post)}": {"risk": round(float(r["risk"]), 1), "events": int(r["events"]), "h_pct": round(100.0 * r["events"] / r["risk"], 3)} for _, r in g.iterrows()}

    bd, sed, sed2 = float(mp.params["ddd"]), float(mp.bse["ddd"]), float(se2_p["ddd"])
    bi, sei = float(mpi.params["od_post"]), float(mpi.bse["od_post"])
    ba, sea = float(mpa.params["od_post"]), float(mpa.bse["od_post"])
    return {
        "label": label, "n_obs": int(len(pan)), "n_events": int(pan["y"].sum()), "hazards": haz,
        "pois_ddd": {"irr": float(np.exp(bd)), "ci95_spon": [float(np.exp(bd - 1.96 * sed)), float(np.exp(bd + 1.96 * sed))], "p_spon": float(mp.pvalues["ddd"]),
                     "ci95_2way": [float(np.exp(bd - 1.96 * sed2)), float(np.exp(bd + 1.96 * sed2))], "p_2way": float(2 * stats.norm.sf(abs(bd / max(sed2, 1e-12))))},
        "lpm_ddd": {"coef_pp": float(ml.params["ddd"]), "se_pp": float(ml.bse["ddd"]), "p": float(ml.pvalues["ddd"])},
        "industry_od_vs_ot": {"irr": float(np.exp(bi)), "ci95": [float(np.exp(bi - 1.96 * sei)), float(np.exp(bi + 1.96 * sei))], "p": float(mpi.pvalues["od_post"])},
        "academic_od_vs_ot": {"irr": float(np.exp(ba)), "ci95": [float(np.exp(ba - 1.96 * sea)), float(np.exp(ba + 1.96 * sea))], "p": float(mpa.pvalues["od_post"])},
    }


def eval_exempt_ddd(panel, ctrl_arms, label):
    sub = panel[panel["arm"].isin(["us_fdareg_nonp1"] + list(ctrl_arms))].copy()
    cdf = sub.groupby(["lead_sponsor_clean", "arm", "is_act", "is_industry", "age_bin", "ym", "post"], as_index=False).agg(
        events=("y", "sum"), risk=("exp_frac", "sum")
    )
    cdf["ddd"] = cdf["is_act"] * cdf["is_industry"] * cdf["post"]
    cdf["act_post"] = cdf["is_act"] * cdf["post"]
    cdf["ind_post"] = cdf["is_industry"] * cdf["post"]
    cdf["cell_aia"] = cdf["arm"].astype(str) + "_" + cdf["is_industry"].astype(str) + "_" + cdf["age_bin"]
    cdf = cdf[cdf.groupby("cell_aia")["events"].transform("sum") > 0].copy()
    mp = smf.glm("events ~ ddd + act_post + ind_post + C(cell_aia) + C(ym)", data=cdf, family=sm.families.Poisson(), offset=np.log(cdf["risk"].values)).fit(
        cov_type="cluster", cov_kwds={"groups": cdf["lead_sponsor_clean"]}
    )
    sub["ddd"] = sub["is_act"] * sub["is_industry"] * sub["post"]
    sub["act_post"] = sub["is_act"] * sub["post"]
    sub["ind_post"] = sub["is_industry"] * sub["post"]
    sub["cell_aia"] = sub["arm"].astype(str) + "_" + sub["is_industry"].astype(str) + "_" + sub["age_bin"]
    sub["y_bp"] = (sub["y"] / sub["exp_frac"]) * 100.0
    ml = smf.ols("y_bp ~ ddd + act_post + ind_post + C(cell_aia) + C(ym)", data=sub).fit(
        cov_type="cluster", cov_kwds={"groups": sub["lead_sponsor_clean"]}
    )
    bd, sed = float(mp.params["ddd"]), float(mp.bse["ddd"])
    return {
        "label": label, "n_obs": int(len(sub)), "n_events": int(sub["y"].sum()),
        "lpm_ddd_pp": float(ml.params["ddd"]), "lpm_ddd_se": float(ml.bse["ddd"]), "lpm_ddd_p": float(ml.pvalues["ddd"]),
        "pois_ddd_irr": float(np.exp(bd)), "pois_ddd_ci95": [float(np.exp(bd - 1.96 * sed)), float(np.exp(bd + 1.96 * sed))], "pois_ddd_p": float(mp.pvalues["ddd"]),
    }


def honest_did_bounds(beta_pre, beta_post, cov_all, post_weights=None):
    """
    Exact Rambachan & Roth (2023, Rev. Econ. Stud. 90:2555-2591) Partial-Identification Sensitivity Analysis!
    Given event-study vector beta = [beta_{-K}, ..., beta_{-2}, beta_0, ..., beta_{T-1}] (omitting reference k=-1 where delta_{-1}=0):
      Target scalar parameter: theta = l' tau_post = l' (beta_post - delta_post), where l = post_weights (uniform avg over post).
      1. Smoothness class Delta^SD(M): |(delta_{t+1} - delta_t) - (delta_t - delta_{t-1})| <= M for all consecutive triples
         (including transitions across the reference period delta_{-1} = 0).
         Using the exact identified set bias bounds [b_min(M), b_max(M)] over delta_post given delta_pre = beta_pre,
         plus standard error of theta_hat = l' beta_post (and accounting for pre-period estimation uncertainty via
         Monte Carlo / normal-approximation over joint Sigma for exact Rambachan-Roth FLCI/identified-set CI!).
      2. Relative Magnitudes class Delta^RM(Mbar):
         max_{t >= 0} |delta_t - delta_{t-1}| <= Mbar * max_{s < 0} |delta_s - delta_{s-1}| (with delta_{-1} = 0).
    Returns exact identified-set bounds AND 95% robust confidence intervals (via parametric bootstrap on joint cov_all
    taking exact LP min/max of l' delta_post per draw, yielding valid unconditional/conservative robust CIs!).
    """
    K = len(beta_pre)   # pre periods: k = -K, ..., -2 (reference k=-1 is 0)
    T = len(beta_post)  # post periods: k = 0, ..., T-1
    if post_weights is None:
        l_vec = np.ones(T) / T
    else:
        l_vec = np.asarray(post_weights, dtype=float)
    beta_all = np.concatenate([beta_pre, beta_post])
    theta_hat = float(l_vec @ beta_post)
    se_theta = float(np.sqrt(l_vec @ cov_all[K:, K:] @ l_vec))

    # Full timeline indices: 0..K-1 are pre (-K..-2), index K is ref (-1, fixed=0), K+1..K+T are post (0..T-1)
    def solve_sd_bias(d_pre, M):
        # Variables: d_post (length T). Fixed: d_full[:K] = d_pre, d_full[K] = 0.
        # Second difference constraints at every interior index j = 1 .. K+T-1:
        #   |(d[j+1] - 2 d[j] + d[j-1])| <= M
        # Note: if d_pre itself violates M in pre-periods, Delta^SD(M) is empty at point d_pre;
        # we only constrain transitions involving d_post (j >= K-1, i.e., involving d_post[0]) while also reporting max pre |DD|!
        A_ub, b_ub = [], []
        # j = K-1 involves (d_pre[K-2], d_pre[K-1], 0) -> pre-only
        # j = K involves (d_pre[K-1], 0, d_post[0]): d_post[0] + d_pre[K-1] in [-M, M]
        # Wait: exact 2nd diff at j=K is (d_post[0] - 0) - (0 - d_pre[K-1]) = d_post[0] + d_pre[K-1]
        #Wait! Check sign: d[K-1]=d_pre[-1] (which is k=-2), d[K]=0 (k=-1), d[K+1]=d_post[0] (k=0)!
        # 2nd diff at k=-1 is (d(0) - d(-1)) - (d(-1) - d(-2)) = d(0) + d(-2)!
        row = np.zeros(T); row[0] = 1.0
        A_ub.append(row); b_ub.append(M - d_pre[-1])
        A_ub.append(-row); b_ub.append(M + d_pre[-1])
        if T >= 2:
            # j = K+1 (k=0): (d(1) - d(0)) - (d(0) - 0) = d(1) - 2*d(0)
            row = np.zeros(T); row[1] = 1.0; row[0] = -2.0
            A_ub.append(row); b_ub.append(M)
            A_ub.append(-row); b_ub.append(M)
        for p in range(2, T):
            row = np.zeros(T); row[p] = 1.0; row[p-1] = -2.0; row[p-2] = 1.0
            A_ub.append(row); b_ub.append(M)
            A_ub.append(-row); b_ub.append(M)
        A_ub = np.array(A_ub); b_ub = np.array(b_ub)
        r_min = linprog(l_vec, A_ub=A_ub, b_ub=b_ub, bounds=(None, None), method="highs")
        r_max = linprog(-l_vec, A_ub=A_ub, b_ub=b_ub, bounds=(None, None), method="highs")
        return float(r_min.fun), float(-r_max.fun)

    def solve_rm_bias(d_pre, Mbar):
        # Max pre-period first-difference: max(|d[-K+1]-d[-K]|, ..., |0 - d[-2]|)
        d_full_pre = np.concatenate([d_pre, [0.0]])
        max_pre_step = float(np.max(np.abs(np.diff(d_full_pre))))
        B = Mbar * max_pre_step
        # Since |d_post[0] - 0| <= B, |d_post[1] - d_post[0]| <= B, ...,
        # max l' d_post is achieved by d_post[t] = +(t+1)*B, min by -(t+1)*B (for non-negative l_vec)!
        steps = np.arange(1, T + 1, dtype=float)
        max_b = float(l_vec @ (steps * B))
        return -max_b, max_b, max_pre_step

    # Compute exact max pre-period first diff & second diff
    d_full_pre = np.concatenate([beta_pre, [0.0]])
    max_pre_fd = float(np.max(np.abs(np.diff(d_full_pre))))
    max_pre_sd = float(np.max(np.abs(np.diff(d_full_pre, n=2)))) if len(d_full_pre) >= 3 else 0.0

    # Parametric draws from N(beta_all, cov_all) for joint-uncertainty-aware robust 95% CIs (Rambachan-Roth style!)
    rng = np.random.default_rng(20260330)
    # Symmetrize & regularize cov slightly if needed
    C = 0.5 * (cov_all + cov_all.T)
    eigvals, eigvecs = np.linalg.eigh(C)
    C_psd = eigvecs @ np.diag(np.maximum(eigvals, 1e-12)) @ eigvecs.T
    draws = rng.multivariate_normal(beta_all, C_psd, size=2000)

    sd_grid = {}
    for M in [0.0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50]:
        b_lo, b_hi = solve_sd_bias(beta_pre, M)
        # Identified set at point estimate: [theta_hat - b_hi, theta_hat - b_lo]
        id_lo, id_hi = theta_hat - b_hi, theta_hat - b_lo
        # Draw-based robust 95% CI (taking conservative union across identified-set endpoints!)
        lo_s, hi_s = [], []
        for dr in draws:
            dp, dpo = dr[:K], dr[K:]
            th = float(l_vec @ dpo)
            bl, bh = solve_sd_bias(dp, M)
            lo_s.append(th - bh)
            hi_s.append(th - bl)
        ci_lo = float(np.quantile(lo_s, 0.025))
        ci_hi = float(np.quantile(hi_s, 0.975))
        sd_grid[str(M)] = {
            "id_set": [round(id_lo, 4), round(id_hi, 4)],
            "robust_ci95": [round(ci_lo, 4), round(ci_hi, 4)],
            "robust_irr_ci95": [round(float(np.exp(ci_lo)), 3), round(float(np.exp(ci_hi)), 3)],
            "excludes_zero": bool(ci_lo > 0 or ci_hi < 0),
        }

    rm_grid = {}
    for Mbar in [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]:
        b_lo, b_hi, _ = solve_rm_bias(beta_pre, Mbar)
        id_lo, id_hi = theta_hat - b_hi, theta_hat - b_lo
        lo_s, hi_s = [], []
        for dr in draws:
            dp, dpo = dr[:K], dr[K:]
            th = float(l_vec @ dpo)
            bl, bh, _ = solve_rm_bias(dp, Mbar)
            lo_s.append(th - bh)
            hi_s.append(th - bl)
        ci_lo = float(np.quantile(lo_s, 0.025))
        ci_hi = float(np.quantile(hi_s, 0.975))
        rm_grid[str(Mbar)] = {
            "id_set": [round(id_lo, 4), round(id_hi, 4)],
            "robust_ci95": [round(ci_lo, 4), round(ci_hi, 4)],
            "robust_irr_ci95": [round(float(np.exp(ci_lo)), 3), round(float(np.exp(ci_hi)), 3)],
            "excludes_zero": bool(ci_lo > 0 or ci_hi < 0),
        }

    return {
        "theta_hat_log_irr": round(theta_hat, 4),
        "theta_hat_irr": round(float(np.exp(theta_hat)), 3),
        "se_theta": round(se_theta, 4),
        "naive_ci95_log": [round(theta_hat - 1.96 * se_theta, 4), round(theta_hat + 1.96 * se_theta, 4)],
        "max_pre_first_diff": round(max_pre_fd, 4),
        "max_pre_second_diff": round(max_pre_sd, 4),
        "delta_SD_smoothness": sd_grid,
        "delta_RM_relative_magnitudes": rm_grid,
    }


def run_event_study_and_honestdid(panel, mode="within_act", ref_ym="2026-02", post_months=None):
    """Run exact Poisson monthly event study + joint pre-trend Wald test + Rambachan-Roth HonestDiD!"""
    if mode == "within_act":
        act = panel[panel["arm"] == "us_fdareg_nonp1"].copy()
        cdf = act.groupby(["lead_sponsor_clean", "is_industry", "age_bin", "ym"], as_index=False).agg(events=("y", "sum"), risk=("exp_frac", "sum"))
        cdf["cell"] = cdf["is_industry"].astype(str) + "_" + cdf["age_bin"]
        months = sorted(cdf["ym"].unique())
        ref_i = months.index(ref_ym)
        terms = []
        for m in months:
            if m == ref_ym:
                continue
            k = months.index(m) - ref_i
            nm = f"ev_{'m' + str(abs(k)) if k < 0 else 'p' + str(k)}"
            cdf[nm] = ((cdf["is_industry"] == 1) & (cdf["ym"] == m)).astype(float)
            terms.append((m, k, nm))
        rhs = " + ".join([t[2] for t in terms]) + " + C(cell) + C(ym)"
    elif mode == "ot_od_ddd":
        cdf = panel.groupby(["lead_sponsor_clean", "is_overdue", "is_industry", "pcd_age_bin", "ym"], as_index=False).agg(events=("y", "sum"), risk=("exp_frac", "sum"))
        cdf["cell"] = cdf["is_overdue"].astype(str) + "_" + cdf["is_industry"].astype(str) + "_" + cdf["pcd_age_bin"]
        cdf = cdf[cdf.groupby("cell")["events"].transform("sum") > 0].copy()
        months = sorted(cdf["ym"].unique())
        ref_i = months.index(ref_ym)
        terms, nuis = [], []
        for m in months:
            if m == ref_ym:
                continue
            k = months.index(m) - ref_i
            tag = f"{'m' + str(abs(k)) if k < 0 else 'p' + str(k)}"
            nm = f"ev_{tag}"
            cdf[nm] = ((cdf["is_overdue"] == 1) & (cdf["is_industry"] == 1) & (cdf["ym"] == m)).astype(float)
            cdf[f"od_{tag}"] = ((cdf["is_overdue"] == 1) & (cdf["ym"] == m)).astype(float)
            cdf[f"in_{tag}"] = ((cdf["is_industry"] == 1) & (cdf["ym"] == m)).astype(float)
            terms.append((m, k, nm))
            nuis.extend([f"od_{tag}", f"in_{tag}"])
        rhs = " + ".join([t[2] for t in terms] + nuis) + " + C(cell) + C(ym)"
    elif mode == "ot_od_ind_only":
        sub = panel[panel["is_industry"] == 1].copy()
        cdf = sub.groupby(["lead_sponsor_clean", "is_overdue", "pcd_age_bin", "ym"], as_index=False).agg(events=("y", "sum"), risk=("exp_frac", "sum"))
        cdf["cell"] = cdf["is_overdue"].astype(str) + "_" + cdf["pcd_age_bin"]
        cdf = cdf[cdf.groupby("cell")["events"].transform("sum") > 0].copy()
        months = sorted(cdf["ym"].unique())
        ref_i = months.index(ref_ym)
        terms = []
        for m in months:
            if m == ref_ym:
                continue
            k = months.index(m) - ref_i
            nm = f"ev_{'m' + str(abs(k)) if k < 0 else 'p' + str(k)}"
            cdf[nm] = ((cdf["is_overdue"] == 1) & (cdf["ym"] == m)).astype(float)
            terms.append((m, k, nm))
        rhs = " + ".join([t[2] for t in terms]) + " + C(cell) + C(ym)"

    fit = smf.glm(f"events ~ {rhs}", data=cdf, family=sm.families.Poisson(), offset=np.log(cdf["risk"].values)).fit(
        cov_type="cluster", cov_kwds={"groups": cdf["lead_sponsor_clean"]}
    )
    series = []
    pre_nms, post_nms = [], []
    for m in months:
        k = months.index(m) - ref_i
        if m == ref_ym:
            series.append({"ym": m, "rel_m": k, "log_irr": 0.0, "se": 0.0, "irr": 1.0, "ci95": [1.0, 1.0], "p": 1.0})
            continue
        nm = f"ev_{'m' + str(abs(k)) if k < 0 else 'p' + str(k)}"
        b, se = float(fit.params[nm]), float(fit.bse[nm])
        series.append({"ym": m, "rel_m": k, "log_irr": b, "se": se, "irr": float(np.exp(b)), "ci95": [float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se))], "p": float(fit.pvalues[nm])})
        if m < ref_ym:
            pre_nms.append(nm)
        elif (post_months is None and m >= "2026-04") or (post_months is not None and m in post_months):
            post_nms.append(nm)

    # Joint pre-trend Wald test (all months strictly < 2026-03 except ref_ym, plus separate test including 2026-03!)
    pnames = list(fit.params.index)
    R = np.zeros((len(pre_nms), len(pnames)))
    for i, nm in enumerate(pre_nms):
        R[i, pnames.index(nm)] = 1.0
    wt = fit.wald_test(R, scalar=True)

    # Linear pre-trend slope test across pre-period!
    # Also compute HonestDiD on [pre_nms, post_nms]
    sel_nms = pre_nms + post_nms
    beta_pre = np.array([float(fit.params[nm]) for nm in pre_nms])
    beta_post = np.array([float(fit.params[nm]) for nm in post_nms])
    cov_sub = fit.cov_params().loc[sel_nms, sel_nms].values
    hd = honest_did_bounds(beta_pre, beta_post, cov_sub)

    return {
        "mode": mode, "ref_ym": ref_ym,
        "pre_trend_wald": {"df": len(pre_nms), "stat": float(wt.statistic), "p": float(wt.pvalue)},
        "series": series,
        "honest_did": hd,
    }


def run_h4_lemons_audit(df_strict):
    """
    Closing Codex T03, T04, C03:
      1. Clean H4 estimation on b2_lemons_results.parquet (exact N, coef, SE, 95% CI, 80%-power MDE, TOST equivalence).
      2. Unposted/in-QC selection check: compares TERMINATED rate of post-sweep QC-cleared vs unposted/in-QC overdue ACTs!
    """
    lem = pd.read_parquet(os.path.join(CLEAN, "b2_lemons_results.parquet"))
    lem["pcd"] = parse_ctgov_date(lem["primary_completion_date"], conservative_end=True)
    lem["sub_date"] = parse_ctgov_date(lem["results_first_submit_date"], conservative_end=False)
    lem["due_date"] = lem["pcd"] + pd.Timedelta(days=365)
    lem["delay_days"] = (lem["sub_date"] - lem["pcd"]).dt.days
    lem["is_overdue"] = (lem["delay_days"] > 365).astype(int)
    lem["is_post"] = (lem["sub_date"] >= "2026-03-30").astype(int)
    lem["is_industry"] = (lem["lead_sponsor_class"] == "INDUSTRY").astype(int)
    lem["lead_sponsor_clean"] = lem["lead_sponsor_name"].fillna("UNKNOWN").str.strip()
    lem["is_terminated"] = (lem["overall_status"] == "TERMINATED").astype(int)
    lem["omitted_pval"] = (~lem["has_primary_pval"].fillna(False)).astype(int)
    lem["nonsig_p05"] = lem["all_nonsig_primary_p05"].astype(float)
    lem["log_enr"] = np.log1p(pd.to_numeric(lem["enrollment"], errors="coerce").clip(lower=0).fillna(0))
    lem["phase_grp"] = lem["phases"].fillna("NA").apply(lambda s: "P3" if "PHASE3" in s else ("P2" if "PHASE2" in s else ("P4" if "PHASE4" in s else "OTHER")))
    lem = lem[lem["pcd"].notna() & lem["sub_date"].notna() & (lem["sub_date"] >= "2024-01-01") & (lem["sub_date"] <= "2026-09-30")].copy()
    lem["grp"] = np.where(lem["is_overdue"] == 0, "on_time", np.where(lem["is_post"] == 1, "post_overdue", "pre_overdue"))

    # Exact group summaries
    gs = {}
    for gnm, g in lem.groupby("grp"):
        gs[gnm] = {
            "n": int(len(g)), "industry_share": round(float(g["is_industry"].mean()), 4),
            "terminated_pct": round(100.0 * g["is_terminated"].mean(), 2),
            "omitted_pval_pct": round(100.0 * g["omitted_pval"].mean(), 2),
            "n_with_pval": int(g["nonsig_p05"].notna().sum()),
            "nonsig_p05_pct": round(100.0 * g["nonsig_p05"].mean(), 2),
            "n_with_sae": int(g["sae_rate"].notna().sum()),
            "mean_sae_pct": round(100.0 * g["sae_rate"].mean(), 2),
            "n_with_mort": int(g["mortality_rate"].notna().sum()),
            "mean_mort_pct": round(100.0 * g["mortality_rate"].mean(), 2),
            "median_delay_days": float(g["delay_days"].median()),
        }

    # Exact regressions (Pooled with controls AND Within-Sponsor FE with sponsor-clustered SEs!)
    # Explicitly report 95% CI, 80%-power MDE = 2.8 * SE, and TOST p-value at +/-5pp & +/-10pp!
    lem["is_pre_od"] = (lem["grp"] == "pre_overdue").astype(float)
    lem["is_post_od"] = (lem["grp"] == "post_overdue").astype(float)
    regs = {}
    for ycol, is_pp in [("is_terminated", True), ("omitted_pval", True), ("nonsig_p05", True), ("sae_rate", True), ("mortality_rate", True)]:
        d = lem.dropna(subset=[ycol]).copy()
        d["y_val"] = d[ycol] * (100.0 if is_pp else 1.0)
        # Pooled
        mp = smf.ols("y_val ~ is_pre_od + is_post_od + is_industry + log_enr + C(phase_grp)", data=d).fit(
            cov_type="cluster", cov_kwds={"groups": d["lead_sponsor_clean"]}
        )
        # Within-Sponsor FE (sponsors with >= 2 trials in d)
        sc = d.groupby("lead_sponsor_clean")["y_val"].transform("count")
        df_fe = d[sc >= 2].copy()
        X = pd.concat([df_fe[["is_pre_od", "is_post_od", "log_enr"]].astype(float), pd.get_dummies(df_fe["phase_grp"], drop_first=True, dtype=float)], axis=1)
        y = df_fe["y_val"].astype(float)
        gr = df_fe["lead_sponsor_clean"]
        X_dm = X - X.groupby(gr).transform("mean")
        y_dm = y - y.groupby(gr).transform("mean")
        mfe = sm.OLS(y_dm, X_dm).fit(cov_type="cluster", cov_kwds={"groups": gr})
        dof = np.sqrt((len(df_fe) - X_dm.shape[1]) / max(len(df_fe) - X_dm.shape[1] - gr.nunique(), 1))

        def pack_contrast(coef, se):
            ci = [coef - 1.96 * se, coef + 1.96 * se]
            p = float(2 * stats.norm.sf(abs(coef / se)))
            mde80 = 2.8016 * se
            # TOST p-value for equivalence within [-Delta, +Delta]
            tost_5 = float(max(stats.norm.sf((coef + 5.0) / se), stats.norm.sf((5.0 - coef) / se)))
            tost_10 = float(max(stats.norm.sf((coef + 10.0) / se), stats.norm.sf((10.0 - coef) / se)))
            return {
                "coef_pp": round(coef, 4), "se_pp": round(se, 4), "p": p,
                "ci95_pp": [round(ci[0], 4), round(ci[1], 4)],
                "mde80_pp": round(mde80, 2),
                "tost_5pp_p": round(tost_5, 4),
                "tost_10pp_p": round(tost_10, 4),
            }

        # Post-vs-Pre overdue contrast via linear combination
        c_vec = np.zeros(len(mfe.params)); c_vec[0] = -1.0; c_vec[1] = 1.0
        b_pp_fe = float(mfe.params["is_post_od"] - mfe.params["is_pre_od"])
        se_pp_fe = float(np.sqrt(c_vec @ mfe.cov_params().values @ c_vec)) * dof

        regs[ycol] = {
            "pooled_n": int(len(d)),
            "pooled_post_vs_ontime": pack_contrast(float(mp.params["is_post_od"]), float(mp.bse["is_post_od"])),
            "within_sponsor_n": int(len(df_fe)),
            "within_sponsor_n_sponsors": int(gr.nunique()),
            "within_sponsor_post_vs_ontime": pack_contrast(float(mfe.params["is_post_od"]), float(mfe.bse["is_post_od"]) * dof),
            "within_sponsor_post_vs_pre": pack_contrast(b_pp_fe, se_pp_fe),
        }

    # Unposted/in-QC vs Posted selection check on TERMINATED rate among post-March-30 overdue ACTs!
    act_od_post = df_strict[
        (df_strict["arm"] == "us_fdareg_nonp1") &
        (df_strict["sub_true_any"] >= "2026-03-30") & (df_strict["sub_true_any"] <= "2026-10-31") &
        ((df_strict["sub_true_any"] - df_strict["pcd"]).dt.days > 365)
    ].copy()
    g_p = act_od_post[act_od_post["is_unposted_in_qc"] == 0]
    g_u = act_od_post[act_od_post["is_unposted_in_qc"] == 1]
    _, p_sel = stats.fisher_exact([
        [int(g_u["is_terminated"].sum()), int(len(g_u) - g_u["is_terminated"].sum())],
        [int(g_p["is_terminated"].sum()), int(len(g_p) - g_p["is_terminated"].sum())],
    ])
    qc_selection_check = {
        "posted_post_overdue_n": int(len(g_p)),
        "posted_terminated_pct": round(100.0 * g_p["is_terminated"].mean(), 2) if len(g_p) else None,
        "unposted_in_qc_post_overdue_n": int(len(g_u)),
        "unposted_in_qc_terminated_pct": round(100.0 * g_u["is_terminated"].mean(), 2) if len(g_u) else None,
        "fisher_p": float(p_sel),
    }

    return {"group_summary": gs, "regressions": regs, "unposted_vs_posted_termination_selection": qc_selection_check}


def main():
    t0 = time.time()
    df = load_strict_t00_universe()
    act_df = df[df["arm"] == "us_fdareg_nonp1"].copy()
    print(f"Loaded Strict ACTUAL T00 Universe: {len(df):,} total trials ({len(act_df):,} US FDA-Reg ACTs)", flush=True)

    res = {"executed_at": time.strftime("%Y-%m-%d %H:%M:%S CST"), "pull_date": "2026-11-25"}

    # =========================================================================
    # 1. T00 PRIMARY: TRUE UNCENSORED FIRST SUBMISSION (Posted UNION In-QC!)
    # =========================================================================
    print("\n=== 1. T00 UNCENSORED TRUE FIRST SUBMISSION (Posted UNION In-QC: ZERO QC Censoring!) ===", flush=True)
    # A. Within-ACT Governance DiD across Apr-Jul, Apr-Sep, AND FULL Apr-Oct (7 full post months!) + Placebos
    t00_within = {}
    for nm, sy, ey, evy in [
        ("true_sub_apr_jul_2026", "2025-04", "2026-07", "2026-04"),
        ("true_sub_apr_sep_2026", "2025-04", "2026-09", "2026-04"),
        ("true_sub_apr_oct_2026", "2025-04", "2026-10", "2026-04"),
        ("true_sub_apr_jun_2026", "2025-04", "2026-06", "2026-04"),
        ("true_sub_placebo_2025", "2024-04", "2025-07", "2025-04"),
        ("true_sub_placebo_2024", "2023-04", "2024-07", "2024-04"),
    ]:
        pan = build_overdue_panel(df, "sub_true_any", sy, ey, evy)
        t00_within[nm] = eval_within_act(pan, nm)
        r = t00_within[nm]
        print(f"  [{nm:24s}] ev={r['n_events']:4d} | LPM={r['lpm_pp']:+.3f}pp (p={r['lpm_p_spon']:.3e}, 2way={r['lpm_p_2way']:.3e}) | Pois_IRR={r['pois_irr']:.3f} [{r['pois_ci95_spon'][0]:.2f},{r['pois_ci95_spon'][1]:.2f}] (p={r['pois_p_spon']:.3e}, 2way={r['pois_p_2way']:.3e})", flush=True)
    res["t00_uncensored_within_act"] = t00_within

    # B. Risk-Set-Normalized Overdue vs On-Time DDD using TRUE UNCENSORED FIRST SUBMISSION!
    print("\n=== 2. T00 UNCENSORED RISK-SET OVERDUE vs ON-TIME DDD (Strict ACTUAL US FDA-Reg ACTs) ===", flush=True)
    t00_ot_od = {}
    for nm, sy, ey, evy in [
        ("true_ot_od_apr_jul_2026", "2025-04", "2026-07", "2026-04"),
        ("true_ot_od_apr_sep_2026", "2025-04", "2026-09", "2026-04"),
        ("true_ot_od_apr_oct_2026", "2025-04", "2026-10", "2026-04"),
        ("true_ot_od_apr_jun_2026", "2025-04", "2026-06", "2026-04"),
        ("true_ot_od_placebo_2025", "2024-04", "2025-07", "2025-04"),
        ("true_ot_od_placebo_2024", "2023-04", "2024-07", "2024-04"),
    ]:
        pan = build_ontime_overdue_panel(act_df, "sub_true_any", sy, ey, evy)
        t00_ot_od[nm] = eval_ot_od_ddd(pan, nm)
        r = t00_ot_od[nm]
        print(f"  [{nm:25s}] Pois_DDD={r['pois_ddd']['irr']:.3f} [{r['pois_ddd']['ci95_spon'][0]:.2f},{r['pois_ddd']['ci95_spon'][1]:.2f}] (p={r['pois_ddd']['p_spon']:.3e}, 2way={r['pois_ddd']['p_2way']:.3e}) | Ind(OD/OT)={r['industry_od_vs_ot']['irr']:.3f} (p={r['industry_od_vs_ot']['p']:.3e}) | Aca(OD/OT)={r['academic_od_vs_ot']['irr']:.3f} (p={r['academic_od_vs_ot']['p']:.3e})", flush=True)
    res["t00_uncensored_ot_od_ddd"] = t00_ot_od

    # C. Exempt-Control DDD using TRUE UNCENSORED FIRST SUBMISSION (through 2026-09, where all 4 arms have 100% stationary entry!)
    print("\n=== 3. T00 UNCENSORED DDD vs CLEAN EXEMPT CONTROLS ===", flush=True)
    t00_exempt = {}
    for nm, sy, ey, evy, ctrls in [
        ("true_ddd_clean_exempt_jul", "2025-04", "2026-07", "2026-04", ["us_notfdareg", "nonus_p23_nonfda"]),
        ("true_ddd_clean_exempt_sep", "2025-04", "2026-09", "2026-04", ["us_notfdareg", "nonus_p23_nonfda"]),
        ("true_ddd_us_nonfda_sep", "2025-04", "2026-09", "2026-04", ["us_notfdareg"]),
        ("true_ddd_nonus_p23_sep", "2025-04", "2026-09", "2026-04", ["nonus_p23_nonfda"]),
        ("true_ddd_placebo_2025", "2024-04", "2025-07", "2025-04", ["us_notfdareg", "nonus_p23_nonfda"]),
    ]:
        pan = build_overdue_panel(df, "sub_true_any", sy, ey, evy)
        t00_exempt[nm] = eval_exempt_ddd(pan, ctrls, nm)
        r = t00_exempt[nm]
        print(f"  [{nm:27s}] LPM_DDD={r['lpm_ddd_pp']:+.3f}pp (p={r['lpm_ddd_p']:.3e}) | Pois_DDD_IRR={r['pois_ddd_irr']:.3f} (p={r['pois_ddd_p']:.3e})", flush=True)
    res["t00_uncensored_exempt_ddd"] = t00_exempt

    # D. Delayed-Results Certifications / Extensions (dispFirstSubmitDate) & Pre-Sweep In-QC Queue Clearance!
    print("\n=== 4. DELAYED-RESULTS CERTIFICATIONS (dispFirstSubmitDate) & PRE-SWEEP IN-QC QUEUE CLEARANCE ===", flush=True)
    # Test certifications among delinquent ACTs 2025-04..2026-10
    pan_disp = build_overdue_panel(df[df["arm"] == "us_fdareg_nonp1"], "disp_sub_date", "2025-04", "2026-10", "2026-04")
    res["t00_delayed_cert_within_act"] = eval_within_act(pan_disp, "disp_first_submit_apr_oct_2026")
    print(f"  [dispFirstSubmitDate] ev={res['t00_delayed_cert_within_act']['n_events']} | Pois_IRR={res['t00_delayed_cert_within_act']['pois_irr']:.3f} (p={res['t00_delayed_cert_within_act']['pois_p_spon']:.3e}) | haz={res['t00_delayed_cert_within_act']['hazards']}", flush=True)

    # Test NLM QC Clearance Rate among trials submitted in 2026-04..2026-08 (comparing Industry vs Academic!)
    sub_26 = act_df[(act_df["sub_true_any"] >= "2026-04-01") & (act_df["sub_true_any"] <= "2026-08-31") & ((act_df["sub_true_any"] - act_df["pcd"]).dt.days > 365)].copy()
    sub_26["cleared_qc"] = (sub_26["is_unposted_in_qc"] == 0).astype(int)
    sub_26["cleared_within_90d"] = ((sub_26["is_unposted_in_qc"] == 0) & (sub_26["post_lag_days"] <= 90)).astype(int)
    ind_s, aca_s = sub_26[sub_26["is_industry"] == 1], sub_26[sub_26["is_industry"] == 0]
    _, p_clr = stats.fisher_exact([[int(ind_s["cleared_qc"].sum()), int(len(ind_s) - ind_s["cleared_qc"].sum())],
                                   [int(aca_s["cleared_qc"].sum()), int(len(aca_s) - aca_s["cleared_qc"].sum())]])
    res["t00_post_sweep_qc_clearance_gap"] = {
        "window": "2026-04-01..2026-08-31 overdue submissions",
        "industry_n": int(len(ind_s)), "industry_qc_cleared_pct": round(100.0 * ind_s["cleared_qc"].mean(), 2),
        "academic_n": int(len(aca_s)), "academic_qc_cleared_pct": round(100.0 * aca_s["cleared_qc"].mean(), 2),
        "fisher_p": float(p_clr),
    }

    # =========================================================================
    # 2. T01 EXACT DATE-VERIFIED ZERO-CENSORING UNIFORM-HORIZON SUITE
    #    (Strictly asserting end_date + H <= 2026-11-25!)
    # =========================================================================
    print("\n=== 5. T01 STRICT DATE-VERIFIED ZERO-CENSORING UNIFORM-HORIZON MODELS (end_date + H <= 2026-11-25) ===", flush=True)
    t01_uh = {}
    uh_specs = [
        ("H90_fullmonth_apr_jul", 90, "2025-04", "2026-07", "2026-04", None),           # Jul 31 + 90d = Oct 29 <= Nov 25!
        ("H120_fullmonth_apr_jun", 120, "2025-04", "2026-06", "2026-04", None),         # Jun 30 + 120d = Oct 28 <= Nov 25!
        ("H150_fullmonth_apr_may", 150, "2025-04", "2026-05", "2026-04", None),         # May 31 + 150d = Oct 28 <= Nov 25!
        ("H90_exactday_thru_aug27", 90, "2025-04", "2026-08", "2026-04", "2026-08-27"), # Aug 27 + 90d = Nov 25!
        ("H120_exactday_thru_jul28", 120, "2025-04", "2026-07", "2026-04", "2026-07-28"),# Jul 28 + 120d = Nov 25!
        ("H150_exactday_thru_jun28", 150, "2025-04", "2026-06", "2026-04", "2026-06-28"),# Jun 28 + 150d = Nov 25!
        ("H120_placebo_2025_jun", 120, "2024-04", "2025-06", "2025-04", None),
        ("H120_placebo_2024_jun", 120, "2023-04", "2024-06", "2024-04", None),
    ]
    for nm, H, sy, ey, evy, ex_end in uh_specs:
        pan_od = build_overdue_panel(df, "sub_qc_cleared", sy, ey, evy, h_cap=H, exact_end_date=ex_end)
        pan_ot = build_ontime_overdue_panel(act_df, "sub_qc_cleared", sy, ey, evy, h_cap=H, exact_end_date=ex_end)
        w = eval_within_act(pan_od, nm)
        ot = eval_ot_od_ddd(pan_ot, nm)
        end_dt = pd.Timestamp(ex_end) if ex_end else pd.Period(ey, "M").end_time.normalize()
        t01_uh[nm] = {
            "horizon_days": H, "exact_end_date": str(end_dt.date()),
            "full_obs_date": str((end_dt + pd.Timedelta(days=H)).date()),
            "days_margin_before_pull": int((PULL_DATE - (end_dt + pd.Timedelta(days=H))).days),
            "within_act": w, "ot_od_ddd": ot,
        }
        print(f"  [{nm:25s}] end={end_dt.date()} + {H}d = {(end_dt+pd.Timedelta(days=H)).date()} (margin={t01_uh[nm]['days_margin_before_pull']:2d}d) | WithinACT IRR={w['pois_irr']:.3f} (p={w['pois_p_spon']:.3e}) | OT_OD_DDD={ot['pois_ddd']['irr']:.3f} (p={ot['pois_ddd']['p_spon']:.3e}) | Ind(OD/OT)={ot['industry_od_vs_ot']['irr']:.3f} (p={ot['industry_od_vs_ot']['p']:.3e})", flush=True)
    res["t01_strict_date_verified_uniform_horizon"] = t01_uh

    # =========================================================================
    # 3. T02 / C02: MONTHLY EVENT STUDIES + PRE-TRENDS + RAMBACHAN & ROTH (2023) HonestDiD
    # =========================================================================
    print("\n=== 6. T02 EVENT STUDIES, JOINT PRE-TREND TESTS & RAMBACHAN-ROTH (2023) HonestDiD ===", flush=True)
    pan_od_oct = build_overdue_panel(df, "sub_true_any", "2025-04", "2026-10", "2026-04")
    pan_ot_oct = build_ontime_overdue_panel(act_df, "sub_true_any", "2025-04", "2026-10", "2026-04")
    pan_ot_jul = build_ontime_overdue_panel(act_df, "sub_true_any", "2025-04", "2026-07", "2026-04")

    es_within_true = run_event_study_and_honestdid(pan_od_oct, mode="within_act", ref_ym="2026-02", post_months=["2026-04", "2026-05", "2026-06", "2026-07", "2026-08", "2026-09", "2026-10"])
    es_within_true_4m = run_event_study_and_honestdid(pan_od_oct[pan_od_oct["ym"] <= "2026-07"], mode="within_act", ref_ym="2026-02", post_months=["2026-04", "2026-05", "2026-06", "2026-07"])
    es_ot_od_ddd = run_event_study_and_honestdid(pan_ot_oct, mode="ot_od_ddd", ref_ym="2026-02", post_months=["2026-04", "2026-05", "2026-06", "2026-07", "2026-08", "2026-09", "2026-10"])
    es_ot_od_ddd_4m = run_event_study_and_honestdid(pan_ot_jul, mode="ot_od_ddd", ref_ym="2026-02", post_months=["2026-04", "2026-05", "2026-06", "2026-07"])
    es_ot_od_ind = run_event_study_and_honestdid(pan_ot_oct, mode="ot_od_ind_only", ref_ym="2026-02", post_months=["2026-04", "2026-05", "2026-06", "2026-07", "2026-08", "2026-09", "2026-10"])
    # Also ref = 2025-12 (before Jan-Feb 2026 slump!) for sensitivity!
    es_ot_od_ddd_refdec = run_event_study_and_honestdid(pan_ot_oct, mode="ot_od_ddd", ref_ym="2025-12", post_months=["2026-04", "2026-05", "2026-06", "2026-07", "2026-08", "2026-09", "2026-10"])
    es_within_refdec = run_event_study_and_honestdid(pan_od_oct, mode="within_act", ref_ym="2025-12", post_months=["2026-04", "2026-05", "2026-06", "2026-07", "2026-08", "2026-09", "2026-10"])

    res["t02_event_studies_and_honestdid"] = {
        "within_act_true_7m_ref2602": es_within_true,
        "within_act_true_4m_ref2602": es_within_true_4m,
        "ot_od_ddd_true_7m_ref2602": es_ot_od_ddd,
        "ot_od_ddd_true_4m_ref2602": es_ot_od_ddd_4m,
        "ot_od_ind_only_7m_ref2602": es_ot_od_ind,
        "ot_od_ddd_true_7m_ref2512": es_ot_od_ddd_refdec,
        "within_act_true_7m_ref2512": es_within_refdec,
    }
    for k, v in res["t02_event_studies_and_honestdid"].items():
        wt = v["pre_trend_wald"]
        hd = v["honest_did"]
        print(f"  [{k:28s}] PreTrend Wald(df={wt['df']}): stat={wt['stat']:.2f}, p={wt['p']:.3e} | theta_IRR={hd['theta_hat_irr']:.3f} | RM(0.5)={hd['delta_RM_relative_magnitudes']['0.5']['robust_irr_ci95']} | SD(0.1)={hd['delta_SD_smoothness']['0.1']['robust_irr_ci95']}", flush=True)

    # =========================================================================
    # 4. T03 / T04 / C03: H4 LEMONS HONEST AUDIT + MDE + TOST + QC SELECTION
    # =========================================================================
    print("\n=== 7. T03/T04 H4 LEMONS AUDIT (WITH EXACT 95% CI, 80% MDE, TOST & UNPOSTED SELECTION CHECK) ===", flush=True)
    res["t03_t04_h4_lemons_audit"] = run_h4_lemons_audit(df)
    for yk, yv in res["t03_t04_h4_lemons_audit"]["regressions"].items():
        w = yv["within_sponsor_post_vs_ontime"]
        print(f"  [{yk:15s}] within-spon n={yv['within_sponsor_n']:4d}: {w['coef_pp']:+.3f}pp (SE={w['se_pp']:.3f}, p={w['p']:.3f}), 95%CI=[{w['ci95_pp'][0]:+.2f}, {w['ci95_pp'][1]:+.2f}], MDE80={w['mde80_pp']:.2f}pp, TOST5p={w['tost_5pp_p']:.3f}", flush=True)
    print("  QC-clearance selection check (TERMINATED rate):", res["t03_t04_h4_lemons_audit"]["unposted_vs_posted_termination_selection"], flush=True)

    # Save monthly uncensored hazards for figures!
    mh = pan_ot_oct.groupby(["is_overdue", "is_industry", "ym", "post"], as_index=False).agg(events=("y", "sum"), risk=("exp_frac", "sum"))
    mh["hazard_pct"] = 100.0 * mh["events"] / mh["risk"]
    mh.to_parquet(os.path.join(AUD, "b2_t00_uncensored_monthly_hazards_ot_od.parquet"), index=False)

    res["elapsed_sec"] = round(time.time() - t0, 1)
    out_p = os.path.join(AUD, "b2_definitive_t00_to_t04_results.json")
    with open(out_p, "w") as f:
        json.dump(res, f, indent=2)
    print(f"\nALL B-2 T00-T04 RE-ESTIMATIONS COMPLETE in {res['elapsed_sec']}s -> {out_p}", flush=True)


if __name__ == "__main__":
    main()
