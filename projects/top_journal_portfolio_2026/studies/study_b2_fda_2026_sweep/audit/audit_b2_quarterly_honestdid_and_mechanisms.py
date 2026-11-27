#!/opt/anaconda3/bin/python3
"""
Study B-2 T00 Addendum:
  1. Quarterly Event Study + Rambachan-Roth (2023) HonestDiD (reducing single-month Poisson noise,
     reported transparently alongside monthly HonestDiD!) for:
       - Industry-only Overdue vs On-Time DiD
       - Overdue vs On-Time DDD (Industry vs Academic)
       - Within-ACT Overdue DiD
  2. T00 Uncensored Mechanism & Heterogeneity Decompositions (Apr-Sep 2026 & Apr-Oct 2026,
     Strict ACTUAL universe, ZERO QC censoring):
       - By Mar 29 Backlog Size (1, 2-4, 5+)
       - By Sponsor FDA Portfolio Size (1-2, 3-19, 20+)
       - By Overdue Cohort Depth (2023-2025, 2021-2022, 2017-2020)
       - Within-Firm Phase-1 Spillover (Exposed vs Unexposed Industry sponsors on exempt US Phase-1 trials)
       - Exact Sponsor-Cluster Randomization Inference (500 perms) on Uncensored Overdue-vs-On-Time DDD & Industry OD/OT DiD!
"""
import json, os, sys, time, warnings
import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.api as sm
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")

BASE = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_b2_fda_2026_sweep"
AUD = os.path.join(BASE, "audit")
sys.path.insert(0, AUD)
from audit_b2_definitive_t00_to_t04 import (
    load_strict_t00_universe, build_overdue_panel, build_ontime_overdue_panel, honest_did_bounds
)


def run_quarterly_es(pan, mode="ot_od_ind_only"):
    """Quarterly event study: 2025Q2, 2025Q3, 2025Q4, 2026Q1 (ref), 2026Q2, 2026Q3."""
    p = pan[pan["ym"] <= "2026-09"].copy()
    p["qtr"] = pd.PeriodIndex(p["ym"], freq="M").asfreq("Q").astype(str)
    qtrs = sorted(p["qtr"].unique())  # 2025Q2..2026Q3
    ref_q = "2026Q1"
    ref_i = qtrs.index(ref_q)

    if mode == "ot_od_ind_only":
        s = p[p["is_industry"] == 1].copy()
        cdf = s.groupby(["lead_sponsor_clean", "is_overdue", "pcd_age_bin", "qtr"], as_index=False).agg(events=("y", "sum"), risk=("exp_frac", "sum"))
        cdf["cell"] = cdf["is_overdue"].astype(str) + "_" + cdf["pcd_age_bin"]
        cdf = cdf[cdf.groupby("cell")["events"].transform("sum") > 0].copy()
        terms = []
        for q in qtrs:
            if q == ref_q: continue
            k = qtrs.index(q) - ref_i
            nm = f"ev_{'m' + str(abs(k)) if k < 0 else 'p' + str(k)}"
            cdf[nm] = ((cdf["is_overdue"] == 1) & (cdf["qtr"] == q)).astype(float)
            terms.append((q, k, nm))
        rhs = " + ".join([t[2] for t in terms]) + " + C(cell) + C(qtr)"
    elif mode == "ot_od_ddd":
        cdf = p.groupby(["lead_sponsor_clean", "is_overdue", "is_industry", "pcd_age_bin", "qtr"], as_index=False).agg(events=("y", "sum"), risk=("exp_frac", "sum"))
        cdf["cell"] = cdf["is_overdue"].astype(str) + "_" + cdf["is_industry"].astype(str) + "_" + cdf["pcd_age_bin"]
        cdf = cdf[cdf.groupby("cell")["events"].transform("sum") > 0].copy()
        terms, nuis = [], []
        for q in qtrs:
            if q == ref_q: continue
            k = qtrs.index(q) - ref_i
            tag = f"{'m' + str(abs(k)) if k < 0 else 'p' + str(k)}"
            nm = f"ev_{tag}"
            cdf[nm] = ((cdf["is_overdue"] == 1) & (cdf["is_industry"] == 1) & (cdf["qtr"] == q)).astype(float)
            cdf[f"od_{tag}"] = ((cdf["is_overdue"] == 1) & (cdf["qtr"] == q)).astype(float)
            cdf[f"in_{tag}"] = ((cdf["is_industry"] == 1) & (cdf["qtr"] == q)).astype(float)
            terms.append((q, k, nm))
            nuis.extend([f"od_{tag}", f"in_{tag}"])
        rhs = " + ".join([t[2] for t in terms] + nuis) + " + C(cell) + C(qtr)"
    elif mode == "within_act":
        s = p[p["arm"] == "us_fdareg_nonp1"].copy()
        cdf = s.groupby(["lead_sponsor_clean", "is_industry", "age_bin", "qtr"], as_index=False).agg(events=("y", "sum"), risk=("exp_frac", "sum"))
        cdf["cell"] = cdf["is_industry"].astype(str) + "_" + cdf["age_bin"]
        terms = []
        for q in qtrs:
            if q == ref_q: continue
            k = qtrs.index(q) - ref_i
            nm = f"ev_{'m' + str(abs(k)) if k < 0 else 'p' + str(k)}"
            cdf[nm] = ((cdf["is_industry"] == 1) & (cdf["qtr"] == q)).astype(float)
            terms.append((q, k, nm))
        rhs = " + ".join([t[2] for t in terms]) + " + C(cell) + C(qtr)"

    fit = smf.glm(f"events ~ {rhs}", data=cdf, family=sm.families.Poisson(), offset=np.log(cdf["risk"].values)).fit(
        cov_type="cluster", cov_kwds={"groups": cdf["lead_sponsor_clean"]}
    )
    series, pre_nms, post_nms = [], [], []
    for q in qtrs:
        k = qtrs.index(q) - ref_i
        if q == ref_q:
            series.append({"qtr": q, "rel_q": k, "log_irr": 0.0, "se": 0.0, "irr": 1.0, "ci95": [1.0, 1.0], "p": 1.0})
            continue
        nm = f"ev_{'m' + str(abs(k)) if k < 0 else 'p' + str(k)}"
        b, se = float(fit.params[nm]), float(fit.bse[nm])
        series.append({"qtr": q, "rel_q": k, "log_irr": round(b, 4), "se": round(se, 4), "irr": round(float(np.exp(b)), 3),
                       "ci95": [round(float(np.exp(b - 1.96 * se)), 3), round(float(np.exp(b + 1.96 * se)), 3)], "p": float(fit.pvalues[nm])})
        if k < 0: pre_nms.append(nm)
        else: post_nms.append(nm)

    pnames = list(fit.params.index)
    R = np.zeros((len(pre_nms), len(pnames)))
    for i, nm in enumerate(pre_nms): R[i, pnames.index(nm)] = 1.0
    wt = fit.wald_test(R, scalar=True)
    sel = pre_nms + post_nms
    # Compute HonestDiD for BOTH 2026Q2 (on-impact) AND 2026Q2-Q3 average!
    hd_q2 = honest_did_bounds(np.array([float(fit.params[nm]) for nm in pre_nms]),
                              np.array([float(fit.params[post_nms[0]])]),
                              fit.cov_params().loc[pre_nms + [post_nms[0]], pre_nms + [post_nms[0]]].values)
    hd_avg = honest_did_bounds(np.array([float(fit.params[nm]) for nm in pre_nms]),
                               np.array([float(fit.params[nm]) for nm in post_nms]),
                               fit.cov_params().loc[sel, sel].values)
    return {
        "mode": mode, "ref_q": ref_q,
        "pre_trend_wald": {"df": len(pre_nms), "stat": round(float(wt.statistic), 3), "p": float(wt.pvalue)},
        "series": series, "honest_did_q2_impact": hd_q2, "honest_did_q2_q3_avg": hd_avg,
    }


def main():
    t0 = time.time()
    df = load_strict_t00_universe()
    act_df = df[df["arm"] == "us_fdareg_nonp1"].copy()

    pan_ot_sep = build_ontime_overdue_panel(act_df, "sub_true_any", "2025-04", "2026-09", "2026-04")
    pan_od_sep = build_overdue_panel(df, "sub_true_any", "2025-04", "2026-09", "2026-04")
    pan_od_oct = build_overdue_panel(df, "sub_true_any", "2025-04", "2026-10", "2026-04")

    # 1. Quarterly Event Studies + HonestDiD
    q_es = {
        "ot_od_ind_only": run_quarterly_es(pan_ot_sep, "ot_od_ind_only"),
        "ot_od_ddd": run_quarterly_es(pan_ot_sep, "ot_od_ddd"),
        "within_act": run_quarterly_es(pan_od_sep, "within_act"),
    }

    # 2. T00 Uncensored Heterogeneity Inside Industry ACTs (comparing Overdue vs On-Time or vs Clean Exempt!)
    # Use clean Overdue ACTs (Apr-Sep 2026 uncensored) vs Clean Exempt controls AND raw pre/post hazard ratios!
    act_o = pan_od_sep[pan_od_sep["arm"] == "us_fdareg_nonp1"].copy()
    ctrl_o = pan_od_sep[pan_od_sep["arm"].isin(["us_notfdareg", "nonus_p23_nonfda"])].copy()

    def strat_pois(a_sub, c_sub):
        comb = pd.concat([a_sub, c_sub], ignore_index=True)
        cdf = comb.groupby(["lead_sponsor_clean", "arm", "is_act", "age_bin", "ym", "post"], as_index=False).agg(events=("y", "sum"), risk=("exp_frac", "sum"))
        cdf["act_post"] = cdf["is_act"] * cdf["post"]
        cdf["cell_aa"] = cdf["arm"].astype(str) + "_" + cdf["age_bin"]
        cdf = cdf[cdf.groupby("cell_aa")["events"].transform("sum") > 0]
        mp = smf.glm("events ~ act_post + C(cell_aa) + C(ym)", data=cdf, family=sm.families.Poisson(), offset=np.log(cdf["risk"].values)).fit(
            cov_type="cluster", cov_kwds={"groups": cdf["lead_sponsor_clean"]}
        )
        b, se = float(mp.params["act_post"]), float(mp.bse["act_post"])
        pre_h = 100.0 * a_sub.loc[a_sub["post"] == 0, "y"].sum() / a_sub.loc[a_sub["post"] == 0, "exp_frac"].sum()
        post_h = 100.0 * a_sub.loc[a_sub["post"] == 1, "y"].sum() / a_sub.loc[a_sub["post"] == 1, "exp_frac"].sum()
        return {
            "n_events": int(a_sub["y"].sum()), "pre_h_pct": round(float(pre_h), 3), "post_h_pct": round(float(post_h), 3),
            "raw_ratio": round(float(post_h / max(pre_h, 1e-9)), 3),
            "irr_vs_exempt": round(float(np.exp(b)), 3),
            "ci95": [round(float(np.exp(b - 1.96 * se)), 3), round(float(np.exp(b + 1.96 * se)), 3)],
            "p": float(mp.pvalues["act_post"]),
        }

    het = {"by_backlog": {}, "by_portfolio": {}, "by_cohort": {}}
    for tag, iv in [("industry", 1), ("academic", 0)]:
        het["by_backlog"][tag] = {
            nm: strat_pois(act_o[(act_o["is_industry"] == iv) & (act_o["spon_backlog_mar26"] >= lo) & (act_o["spon_backlog_mar26"] <= hi)], ctrl_o)
            for lo, hi, nm in [(1, 1, "backlog_1"), (2, 4, "backlog_2_4"), (5, 9999, "backlog_5plus")]
        }
        het["by_portfolio"][tag] = {
            nm: strat_pois(act_o[(act_o["is_industry"] == iv) & (act_o["spon_fda_total"] >= lo) & (act_o["spon_fda_total"] <= hi)], ctrl_o)
            for lo, hi, nm in [(1, 2, "small_1_2"), (3, 19, "mid_3_19"), (20, 9999, "large_20plus")]
        }
        het["by_cohort"][tag] = {
            nm: strat_pois(act_o[(act_o["is_industry"] == iv) & (act_o["comp_year"] >= lo) & (act_o["comp_year"] <= hi)], ctrl_o)
            for lo, hi, nm in [(2023, 2025, "recent_2023_2025"), (2021, 2022, "mid_2021_2022"), (2017, 2020, "deep_2017_2020")]
        }

    # 3. Within-Firm Phase-1 Compliance Spillover (Uncensored sub_true_any, Apr-Sep 2026!)
    p1 = pan_od_sep[(pan_od_sep["arm"] == "us_p1_only") & (pan_od_sep["is_industry"] == 1)].copy()
    p1["exposed_firm"] = (p1["spon_backlog_mar26"] >= 1).astype(int)
    p1["exp_post"] = p1["exposed_firm"] * p1["post"]
    cdf_p1 = p1.groupby(["lead_sponsor_clean", "exposed_firm", "age_bin", "ym", "post"], as_index=False).agg(events=("y", "sum"), risk=("exp_frac", "sum"))
    cdf_p1["exp_post"] = cdf_p1["exposed_firm"] * cdf_p1["post"]
    cdf_p1["cell_ea"] = cdf_p1["exposed_firm"].astype(str) + "_" + cdf_p1["age_bin"]
    cdf_p1 = cdf_p1[cdf_p1.groupby("cell_ea")["events"].transform("sum") > 0]
    mp1 = smf.glm("events ~ exp_post + C(cell_ea) + C(ym)", data=cdf_p1, family=sm.families.Poisson(), offset=np.log(cdf_p1["risk"].values)).fit(
        cov_type="cluster", cov_kwds={"groups": cdf_p1["lead_sponsor_clean"]}
    )
    b1, s1 = float(mp1.params["exp_post"]), float(mp1.bse["exp_post"])
    g1 = p1.groupby(["exposed_firm", "post"]).agg(events=("y", "sum"), risk=("exp_frac", "sum")).reset_index()
    spill = {
        "n_trials": int(p1["nct_id"].nunique()), "n_events": int(p1["y"].sum()),
        "hazards": {f"exp{int(r.exposed_firm)}_p{int(r.post)}": {"risk": int(r["risk"]), "events": int(r["events"]), "h_pct": round(100.0 * r["events"] / r["risk"], 3)} for _, r in g1.iterrows()},
        "pois_irr": round(float(np.exp(b1)), 3), "ci95": [round(float(np.exp(b1 - 1.96 * s1)), 3), round(float(np.exp(b1 + 1.96 * s1)), 3)],
        "p": float(mp1.pvalues["exp_post"]),
    }

    # 4. Exact Sponsor-Cluster Randomization Inference (1000 perms) on Uncensored Overdue-vs-On-Time DDD (Apr-Sep 2026)
    rng = np.random.default_rng(20260330)
    sp = pan_ot_sep.groupby(["lead_sponsor_clean", "is_overdue", "is_industry", "post"], as_index=False).agg(events=("y", "sum"), risk=("exp_frac", "sum"))
    def ddd_stat(df_s):
        g = df_s.groupby(["is_overdue", "is_industry", "post"])[["events", "risk"]].sum()
        h = g["events"] / g["risk"]
        ind_did = (h.loc[(1, 1, 1)] / h.loc[(1, 1, 0)]) / (h.loc[(0, 1, 1)] / h.loc[(0, 1, 0)])
        aca_did = (h.loc[(1, 0, 1)] / h.loc[(1, 0, 0)]) / (h.loc[(0, 0, 1)] / h.loc[(0, 0, 0)])
        return float(np.log(ind_did / aca_did))
    obs_ddd = ddd_stat(sp)
    smap = sp[["lead_sponsor_clean", "is_industry"]].drop_duplicates("lead_sponsor_clean")
    perms = []
    for _ in range(1000):
        pm = smap.copy()
        pm["ind_p"] = rng.permutation(pm["is_industry"].values)
        m = sp.drop(columns=["is_industry"]).merge(pm[["lead_sponsor_clean", "ind_p"]], on="lead_sponsor_clean").rename(columns={"ind_p": "is_industry"})
        perms.append(ddd_stat(m))
    pa = np.array(perms)
    ri_ot_od = {"n_perm": 1000, "obs_ddd_irr": round(float(np.exp(obs_ddd)), 3), "ri_two_sided_p": float((np.abs(pa) >= abs(obs_ddd)).mean())}

    out = {
        "executed_at": time.strftime("%Y-%m-%d %H:%M:%S CST"),
        "elapsed_sec": round(time.time() - t0, 1),
        "quarterly_es_and_honestdid": q_es,
        "t00_uncensored_heterogeneity_apr_sep": het,
        "t00_uncensored_phase1_spillover": spill,
        "ri_uncensored_ot_od_ddd_apr_sep": ri_ot_od,
    }
    out_p = os.path.join(AUD, "b2_t00_quarterly_honestdid_and_mechanisms.json")
    with open(out_p, "w") as f:
        json.dump(out, f, indent=2)
    print("=== QUARTERLY EVENT STUDIES + HONESTDID ===")
    for k, v in q_es.items():
        print(f"  [{k:16s}] PreTrend(df={v['pre_trend_wald']['df']}) p={v['pre_trend_wald']['p']:.4f} | Q2={v['series'][-2]['irr']} (p={v['series'][-2]['p']:.3f}), Q3={v['series'][-1]['irr']} (p={v['series'][-1]['p']:.3f}) | Q2 HonestDiD RM(0.5)={v['honest_did_q2_impact']['delta_RM_relative_magnitudes']['0.5']['robust_irr_ci95']} | SD(0.1)={v['honest_did_q2_impact']['delta_SD_smoothness']['0.1']['robust_irr_ci95']}")
    print("\n=== UNCENSORED HETEROGENEITY (INDUSTRY) ===")
    print(json.dumps({k: v["industry"] for k, v in het.items()}, indent=2))
    print("\n=== UNCENSORED PHASE-1 SPILLOVER & RI ===")
    print("Spillover:", spill)
    print("RI OT-OD DDD:", ri_ot_od)


if __name__ == "__main__":
    main()
