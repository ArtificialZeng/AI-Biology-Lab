#!/opt/anaconda3/bin/python3
"""
Study B-2 Definitive T00 + T06 + T07 Econometric Closure Engine:
  - T00: Uses STRICT combined results-submission proxy = resultsFirstSubmitDate UNION
         derivedSection.miscInfoModule.submissionTracking.estimatedResultsFirstSubmitDate
         (ZERO releaseDate fallback, ZERO unpostedEvents fallback, ZERO dispFirstSubmitDate mixing!).
         Also runs sensitivity excluding all trials with any delayed-results certification (dispFirstSubmitDate).
  - T07: Fixes the pre_nms bug cleanly: defines fixed calendar PRE_MONTHS = 2025-04..2026-02 (11 months),
         fixed transition month = 2026-03 (March 30 sweep), POST_MONTHS = 2026-04..2026-10 (7 months).
         Proves and asserts exact numerical invariance of the joint pre-trend Wald test (df=10) AND
         linear basis transformation (beta, full cov) across normalization baselines ref=2026-02 vs ref=2025-12!
  - T06: Exports complete contiguous chronological (betahat, full sigma, period labels, l_vec) payloads
         (NEVER dropping interior month 2026-03!) for BOTH monthly (ref=2026-02 with 2026-03 at t=+1 weight 0,
         and ref=2026-03) AND quarterly (2025Q2..2026Q3, ref=2026Q1) specifications to be evaluated directly
         by the OFFICIAL CRAN R package HonestDiD (Rambachan & Roth 2023)!
"""
import json, os, sys, time, warnings
import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.api as sm
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")

BASE = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_b2_fda_2026_sweep"
CLEAN = os.path.join(BASE, "data", "clean")
AUD = os.path.join(BASE, "audit")
sys.path.insert(0, os.path.join(BASE, "src"))
from analyze_b2_v2_stationary import parse_ctgov_date, classify_sponsor_group

PULL_DATE = pd.Timestamp("2026-11-25")
PRE_MONTHS_FIXED = [str(p) for p in pd.period_range("2025-04", "2026-02", freq="M")]  # 11 months
TRANSITION_MONTH = "2026-03"  # Sweep executed March 30, 2026
POST_MONTHS_FIXED = [str(p) for p in pd.period_range("2026-04", "2026-10", freq="M")]  # 7 months


def load_clean_t00_universe():
    """Load strict symmetric COMPLETED|TERMINATED + ACTUAL PCD universe with STRICT T00 proxy."""
    d_main = pd.read_parquet(os.path.join(CLEAN, "b2_trial_universe_steady.parquet"))
    d_rec = pd.read_parquet(os.path.join(CLEAN, "b2_us_fdareg_recent_25_26.parquet"))
    d_rec["arm"] = "us_fdareg_nonp1"
    df = pd.concat([d_main, d_rec], ignore_index=True).drop_duplicates("nct_id", keep="first")
    ov = pd.read_parquet(os.path.join(AUD, "b2_t00_unposted_and_cert_overlay.parquet"))
    df = df.merge(ov, on="nct_id", how="left")

    df = df[
        df["overall_status"].isin(["COMPLETED", "TERMINATED"]) &
        (df["primary_completion_type"] == "ACTUAL")
    ].copy()
    df["pcd"] = parse_ctgov_date(df["primary_completion_date"], conservative_end=True)
    df = df[df["pcd"].notna() & (df["pcd"] >= "2017-01-18") & (df["pcd"] <= "2026-08-31")].copy()

    # Strictly separate posted QC-cleared vs explicit estimatedResultsFirstSubmitDate (NO fallback!)
    df["sub_posted_qc"] = parse_ctgov_date(df["results_first_submit_date"], conservative_end=False)
    df["est_sub_strict"] = parse_ctgov_date(df["est_results_first_submit_date"], conservative_end=False)
    df["sub_combined_proxy"] = df["sub_posted_qc"].combine_first(df["est_sub_strict"])
    df["has_disp_cert"] = df["disp_first_submit_date"].notna()

    df["due_date"] = df["pcd"] + pd.Timedelta(days=365)
    df["is_industry"] = (df["lead_sponsor_class"] == "INDUSTRY").astype(int)
    df["lead_sponsor_clean"] = df["lead_sponsor_name"].fillna("UNKNOWN_SPONSOR").str.strip()
    df["is_act"] = (df["arm"] == "us_fdareg_nonp1").astype(int)
    return df


def build_overdue_panel(df, sub_col, start_ym, end_ym, event_ym):
    months = pd.period_range(start_ym, end_ym, freq="M")
    recs = []
    for m in months:
        ms, me = m.start_time, m.end_time
        ar = (df["due_date"] < ms) & (df[sub_col].isna() | (df[sub_col] >= ms))
        sub = df.loc[ar, ["nct_id", "arm", "is_act", "is_industry", "lead_sponsor_clean", "due_date", sub_col, "has_disp_cert"]].copy()
        sub["ym"] = str(m)
        sub["post"] = int(str(m) >= event_ym)
        sub["y"] = ((sub[sub_col] >= ms) & (sub[sub_col] <= me)).astype(int)
        od_m = ((ms - sub["due_date"]).dt.days / 30.4375).clip(lower=0.0)
        sub["age_bin"] = pd.cut(
            od_m, bins=[-0.01, 3.0, 6.0, 12.0, 24.0, 48.0, 999.0],
            labels=["0_3m", "3_6m", "6_12m", "12_24m", "24_48m", "48m_plus"]
        ).astype(str)
        recs.append(sub)
    return pd.concat(recs, ignore_index=True)


def build_ot_od_panel(act_df, sub_col, start_ym, end_ym, event_ym):
    months = pd.period_range(start_ym, end_ym, freq="M")
    recs = []
    for m in months:
        ms, me = m.start_time, m.end_time
        ar = (act_df["pcd"] < ms) & (act_df[sub_col].isna() | (act_df[sub_col] >= ms))
        sub = act_df.loc[ar, ["nct_id", "is_industry", "lead_sponsor_clean", "pcd", "due_date", sub_col, "has_disp_cert"]].copy()
        sub["ym"] = str(m)
        sub["post"] = int(str(m) >= event_ym)
        sub["is_overdue"] = (sub["due_date"] < ms).astype(int)
        sub["y"] = ((sub[sub_col] >= ms) & (sub[sub_col] <= me)).astype(int)
        age_m = ((ms - sub["pcd"]).dt.days / 30.4375).clip(lower=0.0)
        sub["pcd_age_bin"] = pd.cut(
            age_m, bins=[-0.01, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 24.0, 36.0, 60.0, 999.0],
            labels=["0_3", "3_6", "6_9", "9_12", "12_15", "15_18", "18_24", "24_36", "36_60", "60p"]
        ).astype(str)
        recs.append(sub)
    return pd.concat(recs, ignore_index=True)


def fit_monthly_es_fixed_pre(pan, mode, ref_ym):
    """
    T07-CLEAN Monthly Poisson Event Study:
      - Pre-treatment window is FIXED by real calendar date: PRE_MONTHS_FIXED = 2025-04..2026-02 (11 months),
        INDEPENDENT of ref_ym!
      - Joint pre-trend Wald test ALWAYS tests all 10 non-reference pre-months (df = 10), invariant to ref_ym!
      - Returns full chronological ordered (months_est, betahat, cov_matrix) including transition month 2026-03
        so no interior time step is ever dropped when passed to official HonestDiD!
    """
    if mode == "within_act":
        s = pan[pan["arm"] == "us_fdareg_nonp1"].copy()
        cdf = s.groupby(["lead_sponsor_clean", "is_industry", "age_bin", "ym"], as_index=False).agg(events=("y", "sum"), risk=("y", "count"))
        cdf["cell"] = cdf["is_industry"].astype(str) + "_" + cdf["age_bin"]
        months = sorted(cdf["ym"].unique())
        terms = []
        for m in months:
            if m == ref_ym: continue
            nm = "ev_" + m.replace("-", "_")
            cdf[nm] = ((cdf["is_industry"] == 1) & (cdf["ym"] == m)).astype(float)
            terms.append((m, nm))
        rhs = " + ".join([t[1] for t in terms]) + " + C(cell) + C(ym)"
    elif mode == "ot_od_ddd":
        cdf = pan.groupby(["lead_sponsor_clean", "is_overdue", "is_industry", "pcd_age_bin", "ym"], as_index=False).agg(events=("y", "sum"), risk=("y", "count"))
        cdf["cell"] = cdf["is_overdue"].astype(str) + "_" + cdf["is_industry"].astype(str) + "_" + cdf["pcd_age_bin"]
        cdf = cdf[cdf.groupby("cell")["events"].transform("sum") > 0].copy()
        months = sorted(cdf["ym"].unique())
        terms, nuis = [], []
        for m in months:
            if m == ref_ym: continue
            tag = m.replace("-", "_")
            nm = "ev_" + tag
            cdf[nm] = ((cdf["is_overdue"] == 1) & (cdf["is_industry"] == 1) & (cdf["ym"] == m)).astype(float)
            cdf["od_" + tag] = ((cdf["is_overdue"] == 1) & (cdf["ym"] == m)).astype(float)
            cdf["in_" + tag] = ((cdf["is_industry"] == 1) & (cdf["ym"] == m)).astype(float)
            terms.append((m, nm))
            nuis.extend(["od_" + tag, "in_" + tag])
        rhs = " + ".join([t[1] for t in terms] + nuis) + " + C(cell) + C(ym)"
    elif mode == "ot_od_ind_only":
        s = pan[pan["is_industry"] == 1].copy()
        cdf = s.groupby(["lead_sponsor_clean", "is_overdue", "pcd_age_bin", "ym"], as_index=False).agg(events=("y", "sum"), risk=("y", "count"))
        cdf["cell"] = cdf["is_overdue"].astype(str) + "_" + cdf["pcd_age_bin"]
        cdf = cdf[cdf.groupby("cell")["events"].transform("sum") > 0].copy()
        months = sorted(cdf["ym"].unique())
        terms = []
        for m in months:
            if m == ref_ym: continue
            nm = "ev_" + m.replace("-", "_")
            cdf[nm] = ((cdf["is_overdue"] == 1) & (cdf["ym"] == m)).astype(float)
            terms.append((m, nm))
        rhs = " + ".join([t[1] for t in terms]) + " + C(cell) + C(ym)"
    else:
        raise ValueError(mode)

    fit = smf.glm(f"events ~ {rhs}", data=cdf, family=sm.families.Poisson(), offset=np.log(cdf["risk"].values)).fit(
        cov_type="cluster", cov_kwds={"groups": cdf["lead_sponsor_clean"]}
    )

    # FIXED pre-trend test: ALL months in PRE_MONTHS_FIXED except ref_ym (ALWAYS df = 10 when ref_ym in PRE_MONTHS_FIXED,
    # or df = 10 pairwise equality when ref_ym == 2026-03 by testing beta(m) - beta(2026-02) = 0 for m in 2025-04..2026-01!)
    pnames = list(fit.params.index)
    pre_in_model = [(m, nm) for (m, nm) in terms if m in PRE_MONTHS_FIXED]
    if ref_ym in PRE_MONTHS_FIXED:
        assert len(pre_in_model) == 10, f"T07 invariant failed! Expected df=10, got {len(pre_in_model)}"
        R = np.zeros((len(pre_in_model), len(pnames)))
        for i, (_, nm) in enumerate(pre_in_model):
            R[i, pnames.index(nm)] = 1.0
    else:
        # ref_ym is 2026-03 (outside PRE_MONTHS_FIXED): test equality of all 11 pre-months (df = 10)!
        assert len(pre_in_model) == 11
        ref_pre_nm = dict(pre_in_model)["2026-02"]
        others = [(m, nm) for (m, nm) in pre_in_model if m != "2026-02"]
        R = np.zeros((len(others), len(pnames)))
        for i, (_, nm) in enumerate(others):
            R[i, pnames.index(nm)] = 1.0
            R[i, pnames.index(ref_pre_nm)] = -1.0
    wt = fit.wald_test(R, scalar=True)

    # Also test joint pre-trend including March 2026 (2025-04..2026-03, df=11) for complete transparency!
    pre_incl_mar = [(m, nm) for (m, nm) in terms if m <= "2026-03"]
    R_mar = np.zeros((len(pre_incl_mar), len(pnames)))
    for i, (_, nm) in enumerate(pre_incl_mar):
        R_mar[i, pnames.index(nm)] = 1.0
    wt_mar = fit.wald_test(R_mar, scalar=True)

    # Extract full chronologically ordered estimated vector & covariance
    ord_nms = [nm for (_, nm) in terms]
    ord_ms = [m for (m, _) in terms]
    betahat = np.array([float(fit.params[nm]) for nm in ord_nms])
    sigma = fit.cov_params().loc[ord_nms, ord_nms].values

    series = []
    for m in months:
        if m == ref_ym:
            series.append({"ym": m, "log_irr": 0.0, "se": 0.0, "irr": 1.0, "ci95": [1.0, 1.0], "p": 1.0})
        else:
            nm = "ev_" + m.replace("-", "_")
            b, se = float(fit.params[nm]), float(fit.bse[nm])
            series.append({
                "ym": m, "log_irr": round(b, 5), "se": round(se, 5), "irr": round(float(np.exp(b)), 4),
                "ci95": [round(float(np.exp(b - 1.96 * se)), 4), round(float(np.exp(b + 1.96 * se)), 4)],
                "p": float(fit.pvalues[nm]),
            })

    return {
        "mode": mode, "ref_ym": ref_ym,
        "pre_months_fixed": PRE_MONTHS_FIXED,
        "wald_fixed_pre_df10": {"df": int(R.shape[0]), "stat": float(wt.statistic), "p": float(wt.pvalue)},
        "wald_pre_incl_mar_df11": {"df": int(R_mar.shape[0]), "stat": float(wt_mar.statistic), "p": float(wt_mar.pvalue)},
        "ordered_months": ord_ms,
        "betahat": betahat.tolist(),
        "sigma": sigma.tolist(),
        "series": series,
    }


def fit_quarterly_es(pan, mode):
    """Contiguous quarterly event study: 2025Q2, 2025Q3, 2025Q4, 2026Q1 (ref), 2026Q2, 2026Q3."""
    p = pan[pan["ym"] <= "2026-09"].copy()
    p["qtr"] = pd.PeriodIndex(p["ym"], freq="M").asfreq("Q").astype(str)
    qtrs = sorted(p["qtr"].unique())
    ref_q = "2026Q1"
    if mode == "within_act":
        s = p[p["arm"] == "us_fdareg_nonp1"].copy()
        cdf = s.groupby(["lead_sponsor_clean", "is_industry", "age_bin", "qtr"], as_index=False).agg(events=("y", "sum"), risk=("y", "count"))
        cdf["cell"] = cdf["is_industry"].astype(str) + "_" + cdf["age_bin"]
        terms = [(q, "ev_" + q) for q in qtrs if q != ref_q]
        for q, nm in terms:
            cdf[nm] = ((cdf["is_industry"] == 1) & (cdf["qtr"] == q)).astype(float)
        rhs = " + ".join([t[1] for t in terms]) + " + C(cell) + C(qtr)"
    elif mode == "ot_od_ddd":
        cdf = p.groupby(["lead_sponsor_clean", "is_overdue", "is_industry", "pcd_age_bin", "qtr"], as_index=False).agg(events=("y", "sum"), risk=("y", "count"))
        cdf["cell"] = cdf["is_overdue"].astype(str) + "_" + cdf["is_industry"].astype(str) + "_" + cdf["pcd_age_bin"]
        cdf = cdf[cdf.groupby("cell")["events"].transform("sum") > 0].copy()
        terms, nuis = [], []
        for q in qtrs:
            if q == ref_q: continue
            nm = "ev_" + q
            cdf[nm] = ((cdf["is_overdue"] == 1) & (cdf["is_industry"] == 1) & (cdf["qtr"] == q)).astype(float)
            cdf["od_" + q] = ((cdf["is_overdue"] == 1) & (cdf["qtr"] == q)).astype(float)
            cdf["in_" + q] = ((cdf["is_industry"] == 1) & (cdf["qtr"] == q)).astype(float)
            terms.append((q, nm)); nuis.extend(["od_" + q, "in_" + q])
        rhs = " + ".join([t[1] for t in terms] + nuis) + " + C(cell) + C(qtr)"
    elif mode == "ot_od_ind_only":
        s = p[p["is_industry"] == 1].copy()
        cdf = s.groupby(["lead_sponsor_clean", "is_overdue", "pcd_age_bin", "qtr"], as_index=False).agg(events=("y", "sum"), risk=("y", "count"))
        cdf["cell"] = cdf["is_overdue"].astype(str) + "_" + cdf["pcd_age_bin"]
        cdf = cdf[cdf.groupby("cell")["events"].transform("sum") > 0].copy()
        terms = [(q, "ev_" + q) for q in qtrs if q != ref_q]
        for q, nm in terms:
            cdf[nm] = ((cdf["is_overdue"] == 1) & (cdf["qtr"] == q)).astype(float)
        rhs = " + ".join([t[1] for t in terms]) + " + C(cell) + C(qtr)"

    fit = smf.glm(f"events ~ {rhs}", data=cdf, family=sm.families.Poisson(), offset=np.log(cdf["risk"].values)).fit(
        cov_type="cluster", cov_kwds={"groups": cdf["lead_sponsor_clean"]}
    )
    pre_nms = [nm for (q, nm) in terms if q < ref_q]
    pnames = list(fit.params.index)
    R = np.zeros((len(pre_nms), len(pnames)))
    for i, nm in enumerate(pre_nms): R[i, pnames.index(nm)] = 1.0
    wt = fit.wald_test(R, scalar=True)
    ord_nms = [nm for (_, nm) in terms]
    ord_qs = [q for (q, _) in terms]
    series = []
    for q in qtrs:
        if q == ref_q:
            series.append({"qtr": q, "log_irr": 0.0, "se": 0.0, "irr": 1.0, "ci95": [1.0, 1.0], "p": 1.0})
        else:
            nm = "ev_" + q
            b, se = float(fit.params[nm]), float(fit.bse[nm])
            series.append({"qtr": q, "log_irr": round(b, 5), "se": round(se, 5), "irr": round(float(np.exp(b)), 4),
                           "ci95": [round(float(np.exp(b - 1.96 * se)), 4), round(float(np.exp(b + 1.96 * se)), 4)],
                           "p": float(fit.pvalues[nm])})
    return {
        "mode": mode, "ref_q": ref_q,
        "wald_pre_df3": {"df": 3, "stat": float(wt.statistic), "p": float(wt.pvalue)},
        "ordered_quarters": ord_qs,
        "numPrePeriods": 3, "numPostPeriods": 2,
        "betahat": [float(fit.params[nm]) for nm in ord_nms],
        "sigma": fit.cov_params().loc[ord_nms, ord_nms].values.tolist(),
        "series": series,
    }


def eval_disp_excl_sensitivity(act_df):
    """Sensitivity excluding ALL trials that ever filed a delayed-results certification (dispFirstSubmitDate)."""
    sub_df = act_df[~act_df["has_disp_cert"]].copy()
    pan = build_ot_od_panel(sub_df, "sub_combined_proxy", "2025-04", "2026-09", "2026-04")
    cdf = pan.groupby(["lead_sponsor_clean", "is_overdue", "is_industry", "pcd_age_bin", "ym", "post"], as_index=False).agg(events=("y", "sum"), risk=("y", "count"))
    cdf["ddd"] = cdf["is_overdue"] * cdf["is_industry"] * cdf["post"]
    cdf["od_post"] = cdf["is_overdue"] * cdf["post"]
    cdf["ind_post"] = cdf["is_industry"] * cdf["post"]
    cdf["cell"] = cdf["is_overdue"].astype(str) + "_" + cdf["is_industry"].astype(str) + "_" + cdf["pcd_age_bin"]
    cdf = cdf[cdf.groupby("cell")["events"].transform("sum") > 0]
    mp = smf.glm("events ~ ddd + od_post + ind_post + C(cell) + C(ym)", data=cdf, family=sm.families.Poisson(), offset=np.log(cdf["risk"].values)).fit(
        cov_type="cluster", cov_kwds={"groups": cdf["lead_sponsor_clean"]}
    )
    ci = cdf[cdf["is_industry"] == 1].copy()
    ci = ci[ci.groupby("cell")["events"].transform("sum") > 0]
    mpi = smf.glm("events ~ od_post + C(cell) + C(ym)", data=ci, family=sm.families.Poisson(), offset=np.log(ci["risk"].values)).fit(
        cov_type="cluster", cov_kwds={"groups": ci["lead_sponsor_clean"]}
    )
    bd, sed = float(mp.params["ddd"]), float(mp.bse["ddd"])
    bi, sei = float(mpi.params["od_post"]), float(mpi.bse["od_post"])
    return {
        "n_trials_excl_disp": int(len(sub_df)), "n_trials_dropped_disp": int(act_df["has_disp_cert"].sum()),
        "ot_od_ddd_irr": round(float(np.exp(bd)), 4),
        "ot_od_ddd_ci95": [round(float(np.exp(bd - 1.96 * sed)), 4), round(float(np.exp(bd + 1.96 * sed)), 4)],
        "ot_od_ddd_p": float(mp.pvalues["ddd"]),
        "ind_od_ot_irr": round(float(np.exp(bi)), 4),
        "ind_od_ot_ci95": [round(float(np.exp(bi - 1.96 * sei)), 4), round(float(np.exp(bi + 1.96 * sei)), 4)],
        "ind_od_ot_p": float(mpi.pvalues["od_post"]),
    }


def main():
    t0 = time.time()
    df = load_clean_t00_universe()
    act_df = df[df["arm"] == "us_fdareg_nonp1"].copy()

    pan_od = build_overdue_panel(df, "sub_combined_proxy", "2025-04", "2026-10", "2026-04")
    pan_ot = build_ot_od_panel(act_df, "sub_combined_proxy", "2025-04", "2026-10", "2026-04")

    # 1. Run Monthly Event Studies & Verify T07 Mathematical Base-Period Invariance (ref=2026-02 vs ref=2025-12 vs ref=2026-03)!
    es_monthly = {}
    invariance_checks = {}
    for mode, pan in [("ot_od_ind_only", pan_ot), ("ot_od_ddd", pan_ot), ("within_act", pan_od)]:
        r2602 = fit_monthly_es_fixed_pre(pan, mode, "2026-02")
        r2512 = fit_monthly_es_fixed_pre(pan, mode, "2025-12")
        r2603 = fit_monthly_es_fixed_pre(pan, mode, "2026-03")
        es_monthly[f"{mode}_ref2602"] = r2602
        es_monthly[f"{mode}_ref2512"] = r2512
        es_monthly[f"{mode}_ref2603"] = r2603

        # Assert exact Wald invariance across ref=2026-02, ref=2025-12, and ref=2026-03!
        s1 = r2602["wald_fixed_pre_df10"]["stat"]
        s2 = r2512["wald_fixed_pre_df10"]["stat"]
        s3 = r2603["wald_fixed_pre_df10"]["stat"]
        assert abs(s1 - s2) < 1e-4 and abs(s1 - s3) < 1e-4, f"T07 Wald invariance violated for {mode}: {s1} vs {s2} vs {s3}"

        # Assert exact linear basis transformation from ref=2025-12 -> ref=2026-02!
        # In r2512, every month m's coefficient relative to 2026-02 is beta_2512(m) - beta_2512(2026-02),
        # and for m=2025-12 it is 0 - beta_2512(2026-02)!
        ms_2512 = r2512["ordered_months"]
        ms_2602 = r2602["ordered_months"]
        idx_2602_in_2512 = ms_2512.index("2026-02")
        A = np.zeros((len(ms_2602), len(ms_2512)))
        for i, m in enumerate(ms_2602):
            if m != "2025-12":
                A[i, ms_2512.index(m)] = 1.0
            A[i, idx_2602_in_2512] = -1.0
        b_trans = A @ np.array(r2512["betahat"])
        S_trans = A @ np.array(r2512["sigma"]) @ A.T
        max_b_diff = float(np.max(np.abs(b_trans - np.array(r2602["betahat"]))))
        max_s_diff = float(np.max(np.abs(S_trans - np.array(r2602["sigma"]))))
        assert max_b_diff < 1e-6 and max_s_diff < 1e-6, f"T07 basis transform failed: {max_b_diff}, {max_s_diff}"
        invariance_checks[mode] = {
            "wald_df": 10,
            "wald_stat_ref2602": round(s1, 6), "wald_p_ref2602": r2602["wald_fixed_pre_df10"]["p"],
            "wald_stat_ref2512": round(s2, 6), "wald_p_ref2512": r2512["wald_fixed_pre_df10"]["p"],
            "wald_stat_ref2603": round(s3, 6), "wald_p_ref2603": r2603["wald_fixed_pre_df10"]["p"],
            "max_wald_stat_diff": float(max(abs(s1 - s2), abs(s1 - s3))),
            "basis_transform_max_beta_diff": max_b_diff,
            "basis_transform_max_sigma_diff": max_s_diff,
        }

    # 2. Contiguous Quarterly Event Studies (2025Q2..2026Q3, ref=2026Q1)
    es_quarterly = {
        "ot_od_ind_only": fit_quarterly_es(pan_ot, "ot_od_ind_only"),
        "ot_od_ddd": fit_quarterly_es(pan_ot, "ot_od_ddd"),
        "within_act": fit_quarterly_es(pan_od, "within_act"),
    }

    # 3. Certification-exclusion sensitivity check
    disp_sens = eval_disp_excl_sensitivity(act_df)

    # 4. Build exact payload for Official R HonestDiD package!
    #    Every payload has CONTIGUOUS calendar time steps with ZERO omitted interior periods!
    honest_payload = {"specs": {}}
    for mode in ("ot_od_ind_only", "ot_od_ddd", "within_act"):
        # Monthly spec 1: ref=2026-02 -> 10 pre (2025-04..2026-01), 8 post (2026-03..2026-10)
        #   - l_vec_apr_impact = [0, 1, 0, 0, 0, 0, 0, 0] (April 2026 immediate impact)
        #   - l_vec_apr_sep_avg = [0, 1/6, 1/6, 1/6, 1/6, 1/6, 1/6, 0] (Apr-Sep 6m avg, weight 0 on transition Mar)
        #   - l_vec_apr_oct_avg = [0, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7] (Apr-Oct 7m avg, weight 0 on transition Mar)
        r2602 = es_monthly[f"{mode}_ref2602"]
        assert r2602["ordered_months"][:10] == PRE_MONTHS_FIXED[:10]
        assert r2602["ordered_months"][10:] == ["2026-03"] + POST_MONTHS_FIXED
        honest_payload["specs"][f"monthly_{mode}_ref2602_apr_impact"] = {
            "description": f"Monthly {mode} (ref=2026-02, contiguous incl 2026-03 at t=+1): April 2026 on-impact",
            "periods": r2602["ordered_months"], "ref_period": "2026-02",
            "numPrePeriods": 10, "numPostPeriods": 8,
            "betahat": r2602["betahat"], "sigma": r2602["sigma"],
            "l_vec": [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        }
        honest_payload["specs"][f"monthly_{mode}_ref2602_apr_sep_avg"] = {
            "description": f"Monthly {mode} (ref=2026-02, contiguous incl 2026-03 at t=+1): Apr-Sep 6m avg",
            "periods": r2602["ordered_months"], "ref_period": "2026-02",
            "numPrePeriods": 10, "numPostPeriods": 8,
            "betahat": r2602["betahat"], "sigma": r2602["sigma"],
            "l_vec": [0.0, 1/6, 1/6, 1/6, 1/6, 1/6, 1/6, 0.0],
        }
        # Monthly spec 2: ref=2026-03 (29/31 days pre-sweep) -> 11 pre (2025-04..2026-02), 7 post (2026-04..2026-10)
        r2603 = es_monthly[f"{mode}_ref2603"]
        assert r2603["ordered_months"][:11] == PRE_MONTHS_FIXED
        assert r2603["ordered_months"][11:] == POST_MONTHS_FIXED
        honest_payload["specs"][f"monthly_{mode}_ref2603_apr_sep_avg"] = {
            "description": f"Monthly {mode} (ref=2026-03, contiguous 11 pre + 7 post): Apr-Sep 6m avg",
            "periods": r2603["ordered_months"], "ref_period": "2026-03",
            "numPrePeriods": 11, "numPostPeriods": 7,
            "betahat": r2603["betahat"], "sigma": r2603["sigma"],
            "l_vec": [1/6, 1/6, 1/6, 1/6, 1/6, 1/6, 0.0],
        }
        # Quarterly spec: ref=2026Q1 -> 3 pre (2025Q2..2025Q4), 2 post (2026Q2..2026Q3)
        rq = es_quarterly[mode]
        honest_payload["specs"][f"quarterly_{mode}_ref2026Q1_Q2_impact"] = {
            "description": f"Quarterly {mode} (ref=2026Q1, 3 pre + 2 post): 2026Q2 impact",
            "periods": rq["ordered_quarters"], "ref_period": "2026Q1",
            "numPrePeriods": 3, "numPostPeriods": 2,
            "betahat": rq["betahat"], "sigma": rq["sigma"],
            "l_vec": [1.0, 0.0],
        }
        honest_payload["specs"][f"quarterly_{mode}_ref2026Q1_Q2Q3_avg"] = {
            "description": f"Quarterly {mode} (ref=2026Q1, 3 pre + 2 post): 2026Q2-Q3 avg",
            "periods": rq["ordered_quarters"], "ref_period": "2026Q1",
            "numPrePeriods": 3, "numPostPeriods": 2,
            "betahat": rq["betahat"], "sigma": rq["sigma"],
            "l_vec": [0.5, 0.5],
        }

    with open(os.path.join(AUD, "b2_t06_honestdid_input_payload.json"), "w") as f:
        json.dump(honest_payload, f, indent=2)

    out = {
        "executed_at": time.strftime("%Y-%m-%d %H:%M:%S CST"),
        "elapsed_sec": round(time.time() - t0, 1),
        "t07_invariance_proof": invariance_checks,
        "t00_excl_disp_cert_sensitivity": disp_sens,
        "monthly_es_fixed_pre": {k: {kk: vv for kk, vv in v.items() if kk != "sigma"} for k, v in es_monthly.items()},
        "quarterly_es": {k: {kk: vv for kk, vv in v.items() if kk != "sigma"} for k, v in es_quarterly.items()},
    }
    with open(os.path.join(AUD, "b2_t07_invariance_and_es_clean.json"), "w") as f:
        json.dump(out, f, indent=2)

    print("=== T07 BASE-PERIOD INVARIANCE VERIFICATION (ALL df=10) ===")
    print(json.dumps(invariance_checks, indent=2))
    print("\n=== T00 EXCLUDING DELAYED-RESULTS CERTIFICATION SENSITIVITY ===")
    print(json.dumps(disp_sens, indent=2))
    print("\n=== QUARTERLY EVENT STUDIES (ref=2026Q1) ===")
    for k, v in es_quarterly.items():
        print(f"  {k:16s}: Wald(df=3) p={v['wald_pre_df3']['p']:.4f} | Q2={v['series'][-2]['irr']} (p={v['series'][-2]['p']:.4f}) | Q3={v['series'][-1]['irr']} (p={v['series'][-1]['p']:.4f})")


if __name__ == "__main__":
    main()
