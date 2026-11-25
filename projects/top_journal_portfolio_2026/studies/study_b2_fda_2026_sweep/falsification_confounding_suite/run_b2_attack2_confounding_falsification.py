#!/opt/anaconda3/bin/python3
"""
Study B-2 — Attack 2: Contemporaneous-Confounding & Exact-Timing Adversarial Falsification Suite
==============================================================================================
CLASSIFICATION: EXPLICITLY POST-HOC EXPLORATORY FALSIFICATION (Locked in NOVELTY_IDENTIFICATION_GATE.md
on 2026-11-25; NEVER presented as pre-registered). All negative results and attenuations preserved verbatim!

Tests:
  2A. Exact Daily (2026-03-01..2026-05-15) & Weekly (aligned to Monday 2026-03-30) Combined-Proxy Hazard
      Discontinuity separating:
        - Pre-sweep March (2026-03-01..2026-03-29) -> tests whether March's 1.86%->3.69% rebound preceded the sweep!
        - Private FDA Email Window (2026-03-30..2026-04-12, 2 weeks before public press release!)
        - Public FDA Press Release Window (2026-04-13..2026-06-30)
  2B. Exact High-Dimensional Fixed-Effects Poisson QMLE (exact PPMLHDFE concentrated likelihood with
      sponsor-clustered sandwich SEs, validated against statsmodels GLM on low-dim FE):
        - Sponsor × Month FE (compares the EXACT SAME sponsor's Overdue vs On-Time ACTs in the EXACT SAME month!)
        - Therapeutic-Area × Month FE
        - Phase × Intervention-Type × Month FE
  2C. On-Time Mechanical Deadline Cliff & Robustness Battery:
        - Cliff-Free On-Time Controls (0-9m post-PCD only; and 0-6m post-PCD only, removing the 10-12m rush!)
        - Fine 1-month PCD-age fixed effects (replacing coarse bins!)
        - Donut-March (2025-04..2026-02 vs 2026-04..2026-09, excluding 2026-03!)
        - Leave-One-Sponsor-Out across top 20 overdue sponsors!
"""

import json
import os
import re
import sys
import warnings
import numpy as np
import pandas as pd
import scipy.optimize as opt
import scipy.stats as stats
import statsmodels.api as sm
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")

BASE = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_b2_fda_2026_sweep"
AUD = os.path.join(BASE, "audit")
OUT_DIR = os.path.join(BASE, "falsification_confounding_suite")
os.makedirs(OUT_DIR, exist_ok=True)

sys.path.insert(0, AUD)
from audit_b2_t00_t06_t07_closure import load_clean_t00_universe, build_ot_od_panel


def fit_ppml_fe(df, y_col, risk_col, x_cols, fe_col, cluster_col):
    """
    Exact Poisson QMLE with arbitrary high-dimensional absorbed fixed effects `fe_col`
    (Hausman-Hall-Griliches 1984 / Correia-Guimaraes-Zylkin 2020 PPMLHDFE analytical concentration)
    and exact sponsor-clustered sandwich variance-covariance matrix (accounting for absorbed FE scores!).

    Math verification:
      For group g in fe_col, concentrated MLE satisfies:
        alpha_g(beta) = Y_g / sum_{i in g} R_i exp(x_i' beta)
      Dropping all-zero FE groups (Y_g == 0) is exact because alpha_g = 0 contributes 0 to score/Hessian.
      Concentrated log-likelihood (up to constant):
        L_c(beta) = sum_i y_i (x_i' beta) - sum_g Y_g log(sum_{i in g} R_i exp(x_i' beta))
      Score for observation i in group g:
        s_i(beta) = (y_i - mu_i) * (x_i - xbar_g),  where mu_i = alpha_g R_i exp(x_i' beta),
        and xbar_g = sum_{j in g} (mu_j / Y_g) x_j  (exact Frisch-Waugh-Lovell weighted projection!).
      Hessian (concentrated):
        H(beta) = sum_i mu_i (x_i - xbar_g)(x_i - xbar_g)'  (exact!).
      Clustered sandwich variance:
        V_cl = (G/(G-1)) * H^{-1} (sum_c S_c S_c') H^{-1}, where S_c = sum_{i in cluster c} s_i(beta).
    """
    # Drop all-zero FE groups (exact for Poisson FE!)
    g_sum = df.groupby(fe_col)[y_col].transform("sum")
    d = df[g_sum > 0].copy().reset_index(drop=True)

    y = d[y_col].values.astype(float)
    R = d[risk_col].values.astype(float)
    X = d[x_cols].values.astype(float)
    g_codes, g_uniques = pd.factorize(d[fe_col].values)
    c_codes, c_uniques = pd.factorize(d[cluster_col].values)
    n_g = len(g_uniques)
    n_c = len(c_uniques)
    Y_g = np.bincount(g_codes, weights=y, minlength=n_g)

    def neg_ll_and_grad(beta):
        eta = np.clip(X @ beta, -30.0, 30.0)
        w = R * np.exp(eta)
        denom = np.bincount(g_codes, weights=w, minlength=n_g)
        ll = np.dot(y, eta) - np.dot(Y_g, np.log(denom))
        # mu_i = Y_{g(i)} * w_i / denom_{g(i)}
        mu = (Y_g[g_codes] / denom[g_codes]) * w
        resid = y - mu
        grad = X.T @ resid
        return -ll, -grad

    def hess_fn(beta):
        eta = np.clip(X @ beta, -30.0, 30.0)
        w = R * np.exp(eta)
        denom = np.bincount(g_codes, weights=w, minlength=n_g)
        mu = (Y_g[g_codes] / denom[g_codes]) * w
        # xbar_g = (sum_{i in g} mu_i x_i) / Y_g
        k = X.shape[1]
        xbar = np.zeros((n_g, k))
        for j in range(k):
            xbar[:, j] = np.bincount(g_codes, weights=mu * X[:, j], minlength=n_g) / Y_g
        X_tilde = X - xbar[g_codes]
        H = (X_tilde * mu[:, None]).T @ X_tilde
        return H

    res = opt.minimize(
        neg_ll_and_grad, np.zeros(X.shape[1]), jac=True, hess=hess_fn,
        method="trust-ncg", options={"gtol": 1e-10, "maxiter": 200}
    )
    if not res.success or np.max(np.abs(res.jac)) > 1e-6:
        res = opt.minimize(neg_ll_and_grad, res.x, jac=True, method="BFGS", options={"gtol": 1e-10, "maxiter": 500})

    beta = res.x
    eta = np.clip(X @ beta, -30.0, 30.0)
    w = R * np.exp(eta)
    denom = np.bincount(g_codes, weights=w, minlength=n_g)
    mu = (Y_g[g_codes] / denom[g_codes]) * w
    k = X.shape[1]
    xbar = np.zeros((n_g, k))
    for j in range(k):
        xbar[:, j] = np.bincount(g_codes, weights=mu * X[:, j], minlength=n_g) / Y_g
    X_tilde = X - xbar[g_codes]
    H = (X_tilde * mu[:, None]).T @ X_tilde
    H_inv = np.linalg.pinv(H)

    # Exact observation scores projected orthogonal to FE: s_i = (y_i - mu_i) * X_tilde_i
    # Note: sum_{i in g} (y_i - mu_i) = 0, so (y_i - mu_i)*xbar_g sums to 0 within FE group g,
    # while when FE is nested inside cluster (Sponsor x Month inside Sponsor), X_tilde vs X gives IDENTICAL S_c;
    # when FE crosses clusters (e.g. TherapeuticArea x Month across Sponsors), X_tilde is the EXACT profile score!
    scores_obs = X_tilde * (y - mu)[:, None]
    S_c = np.zeros((n_c, k))
    for j in range(k):
        S_c[:, j] = np.bincount(c_codes, weights=scores_obs[:, j], minlength=n_c)
    meat = S_c.T @ S_c
    # Small-sample cluster adjustment G/(G-1)
    adj = n_c / max(1, n_c - 1)
    V = adj * (H_inv @ meat @ H_inv)
    se = np.sqrt(np.maximum(np.diag(V), 0.0))
    z = beta / np.where(se > 0, se, np.nan)
    pval = 2.0 * (1.0 - stats.norm.cdf(np.abs(z)))
    return {
        "converged": bool(res.success),
        "max_grad": float(np.max(np.abs(res.jac))),
        "n_obs_active": int(len(d)),
        "n_events_active": int(y.sum()),
        "n_fe_groups_active": int(n_g),
        "n_clusters_active": int(n_c),
        "params": {c: float(beta[i]) for i, c in enumerate(x_cols)},
        "se": {c: float(se[i]) for i, c in enumerate(x_cols)},
        "irr": {c: float(np.exp(beta[i])) for i, c in enumerate(x_cols)},
        "ci95": {c: [float(np.exp(beta[i] - 1.96 * se[i])), float(np.exp(beta[i] + 1.96 * se[i]))] for i, c in enumerate(x_cols)},
        "p": {c: float(pval[i]) for i, c in enumerate(x_cols)},
    }


def classify_therapeutic_area(cond):
    if pd.isna(cond):
        return "Other"
    s = str(cond).lower()
    if any(k in s for k in ["cancer", "carcinoma", "tumor", "tumour", "leukemia", "lymphoma", "melanoma", "myeloma", "sarcoma", "glioma", "neoplasm", "oncolog", "glioblastoma", "metastat"]):
        return "Oncology"
    if any(k in s for k in ["heart", "cardiac", "cardio", "coronary", "hypertension", "atrial", "arrhythm", "stroke", "vascular", "thromb", "atheroscler", "aortic", "myocardial", "lipoprotein", "cholesterol"]):
        return "Cardiovascular"
    if any(k in s for k in ["alzheimer", "parkinson", "epilep", "migraine", "depress", "schizo", "bipolar", "anxiety", "autism", "neuro", "brain", "dementia", "multiple sclerosis", "pain", "psych", "adhd", "insomnia", "opioid", "alcohol"]):
        return "Neuro_Psych"
    if any(k in s for k in ["covid", "hiv", "hepatitis", "influenza", "malaria", "tuberculosis", "infect", "virus", "viral", "bacter", "sepsis", "pneumonia", "vaccine", "rsv", "herpes"]):
        return "Infectious"
    if any(k in s for k in ["diabet", "obes", "thyroid", "metabol", "lipid", "renal", "kidney", "liver", "nash", "mash", "steato", "endocrin"]):
        return "Metabolic_Renal"
    if any(k in s for k in ["asthma", "copd", "arthritis", "lupus", "psoriasis", "crohn", "colitis", "dermatitis", "eczema", "autoimmune", "immun", "inflam", "fibrosis", "allergy"]):
        return "Immune_Resp"
    return "Other"


def main():
    print("Loading strict symmetric T00 universe...")
    df = load_clean_t00_universe()
    act = df[df["is_act"] == 1].copy()
    act["ta"] = act["primary_condition"].apply(classify_therapeutic_area)
    act["phase_clean"] = act["phases"].fillna("NA").astype(str).str.strip()
    act["intv_clean"] = act["intervention_types"].fillna("NA").astype(str).str.split("|").str[0].str.strip()
    print(f"Total clean symmetric ACTs: {len(act):,} (Industry={act['is_industry'].sum():,}, Academic={(act['is_industry']==0).sum():,})")

    # ==========================================================================
    # VERIFY PPMLHDFE ENGINE AGAINST STATSMODELS GLM ON BASELINE MONTHLY SPEC!
    # ==========================================================================
    pan_base = build_ot_od_panel(act, "sub_combined_proxy", "2025-04", "2026-09", "2026-04")
    # Build same aggregated panel as baseline OT/OD DDD
    cdf0 = pan_base.groupby(["lead_sponsor_clean", "is_overdue", "is_industry", "pcd_age_bin", "ym", "post"], as_index=False).agg(events=("y", "sum"), risk=("y", "count"))
    cdf0["ddd"] = cdf0["is_overdue"] * cdf0["is_industry"] * cdf0["post"]
    cdf0["od_post"] = cdf0["is_overdue"] * cdf0["post"]
    cdf0["ind_post"] = cdf0["is_industry"] * cdf0["post"]
    cdf0["cell"] = cdf0["is_overdue"].astype(str) + "_" + cdf0["is_industry"].astype(str) + "_" + cdf0["pcd_age_bin"]
    cdf0 = cdf0[cdf0.groupby("cell")["events"].transform("sum") > 0].copy()

    # Baseline statsmodels GLM
    fit_sm = smf.glm("events ~ ddd + od_post + ind_post + C(cell) + C(ym)", data=cdf0,
                     family=sm.families.Poisson(), offset=np.log(cdf0["risk"].values)).fit(
        cov_type="cluster", cov_kwds={"groups": cdf0["lead_sponsor_clean"]}
    )
    # Cross-check our fit_ppml_fe absorbing C(ym) with explicit C(cell) dummies!
    cell_dums = pd.get_dummies(cdf0["cell"], prefix="c", drop_first=True, dtype=float)
    cdf0_chk = pd.concat([cdf0.reset_index(drop=True), cell_dums.reset_index(drop=True)], axis=1)
    x_chk = ["ddd", "od_post", "ind_post"] + list(cell_dums.columns)
    chk_res = fit_ppml_fe(cdf0_chk, "events", "risk", x_chk, "ym", "lead_sponsor_clean")
    print(f"\n=== PPML ENGINE VALIDATION vs STATSMODELS GLM ===")
    print(f"  statsmodels baseline OT/OD DDD: beta={fit_sm.params['ddd']:.6f}, SE={fit_sm.bse['ddd']:.6f}, IRR={np.exp(fit_sm.params['ddd']):.4f}, p={fit_sm.pvalues['ddd']:.6g}")
    print(f"  fit_ppml_fe baseline OT/OD DDD: beta={chk_res['params']['ddd']:.6f}, SE={chk_res['se']['ddd']:.6f}, IRR={chk_res['irr']['ddd']:.4f}, p={chk_res['p']['ddd']:.6g}")
    assert abs(fit_sm.params["ddd"] - chk_res["params"]["ddd"]) < 1e-5, "PPML beta mismatch!"
    assert abs(fit_sm.bse["ddd"] - chk_res["se"]["ddd"]) < 5e-4, "PPML SE mismatch!"

    # ==========================================================================
    # TEST 2A: EXACT WEEKLY & DAILY DISCONTINUITY AROUND 2026-03-30 vs 2026-04-13
    # ==========================================================================
    print("\n=== TEST 2A: EXACT WEEKLY & DAILY TIMING AROUND 2026-03-30 (PRIVATE) vs 2026-04-13 (PUBLIC) ===")
    sub_col = "sub_combined_proxy"
    # Weekly bins aligned strictly to Monday 2026-03-30! (Week 0 = 2026-03-30..2026-04-05)
    t0 = pd.Timestamp("2026-03-30")
    weekly_rows = []
    for w in range(-16, 20):  # 16 weeks before (Dec 2025) through 19 weeks after (Aug 2026)
        ws = t0 + pd.Timedelta(days=7 * w)
        we = ws + pd.Timedelta(days=6)
        ar = (act["pcd"] < ws) & (act[sub_col].isna() | (act[sub_col] >= ws))
        sub = act.loc[ar, ["is_industry", "due_date", sub_col]].copy()
        sub["is_overdue"] = (sub["due_date"] < ws).astype(int)
        sub["y"] = ((sub[sub_col] >= ws) & (sub[sub_col] <= we)).astype(int)
        for od in (0, 1):
            for ind in (0, 1):
                g = sub[(sub["is_overdue"] == od) & (sub["is_industry"] == ind)]
                ev = int(g["y"].sum())
                rk = int(len(g))
                # Annualize or express per-30-days (= weekly_rate * 30/7 * 100%) for intuitive comparison to monthly %!
                weekly_rows.append({
                    "rel_week": w,
                    "week_start": ws.strftime("%Y-%m-%d"),
                    "week_end": we.strftime("%Y-%m-%d"),
                    "is_overdue": od,
                    "is_industry": ind,
                    "events": ev,
                    "risk": rk,
                    "weekly_pct": round(100.0 * ev / rk, 4) if rk > 0 else 0.0,
                    "equiv_30d_pct": round(100.0 * (ev / rk) * (30.4375 / 7.0), 4) if rk > 0 else 0.0,
                })
    wdf = pd.DataFrame(weekly_rows)
    wdf.to_csv(os.path.join(OUT_DIR, "b2_attack2_weekly_hazards.csv"), index=False)

    # Print Industry Overdue & Industry On-Time & Academic Overdue weekly around sweep (rel_week -8 to +10)
    wpiv = wdf[(wdf["rel_week"] >= -8) & (wdf["rel_week"] <= 10)].pivot(
        index=["rel_week", "week_start", "week_end"],
        columns=["is_overdue", "is_industry"],
        values=["events", "equiv_30d_pct"]
    )
    print("\nWeekly submissions & 30d-equiv hazard (%) around 2026-03-30 (Week 0) & 2026-04-13 (Week +2):")
    for idx, row in wpiv.iterrows():
        rw, ws, we = idx
        tag = " <-- PRIVATE EMAILS (Mar 30)" if rw == 0 else (" <-- PUBLIC PRESS (Apr 13)" if rw == 2 else "")
        print(
            f"  W{rw:+03d} ({ws}..{we}) | "
            f"Ind_OD: {int(row[('events',1,1)]):3d} ({row[('equiv_30d_pct',1,1)]:5.2f}%) | "
            f"Ind_OT: {int(row[('events',0,1)]):3d} ({row[('equiv_30d_pct',0,1)]:5.2f}%) | "
            f"Acad_OD: {int(row[('events',1,0)]):3d} ({row[('equiv_30d_pct',1,0)]:5.2f}%) | "
            f"Acad_OT: {int(row[('events',0,0)]):3d} ({row[('equiv_30d_pct',0,0)]:5.2f}%){tag}"
        )

    # Daily breakdown for March 1, 2026 .. May 15, 2026 to diagnose March 2026's 1.86%->3.69% rebound!
    daily_rows = []
    for d_ts in pd.date_range("2026-02-01", "2026-05-15", freq="D"):
        ar = (act["pcd"] < d_ts) & (act[sub_col].isna() | (act[sub_col] >= d_ts))
        sub = act.loc[ar, ["is_industry", "due_date", sub_col]]
        is_od = (sub["due_date"] < d_ts)
        y_d = (sub[sub_col] == d_ts)
        ind_od_ev = int((y_d & is_od & (sub["is_industry"] == 1)).sum())
        ind_od_rk = int((is_od & (sub["is_industry"] == 1)).sum())
        acad_od_ev = int((y_d & is_od & (sub["is_industry"] == 0)).sum())
        acad_od_rk = int((is_od & (sub["is_industry"] == 0)).sum())
        ind_ot_ev = int((y_d & (~is_od) & (sub["is_industry"] == 1)).sum())
        ind_ot_rk = int(((~is_od) & (sub["is_industry"] == 1)).sum())
        daily_rows.append({
            "date": d_ts.strftime("%Y-%m-%d"),
            "dow": d_ts.day_name()[:3],
            "ind_od_events": ind_od_ev, "ind_od_risk": ind_od_rk,
            "acad_od_events": acad_od_ev, "acad_od_risk": acad_od_rk,
            "ind_ot_events": ind_ot_ev, "ind_ot_risk": ind_ot_rk,
        })
    ddf = pd.DataFrame(daily_rows)
    ddf.to_csv(os.path.join(OUT_DIR, "b2_attack2_daily_feb_may2026.csv"), index=False)

    # Exact sub-period comparison (CRITICAL MARCH REBOUND TEST!)
    sub_windows = [
        ("Feb 1-28 (Pre-Sweep)", "2026-02-01", "2026-02-28"),
        ("Mar 1-15 (Pre-Sweep Early Mar)", "2026-03-01", "2026-03-15"),
        ("Mar 16-29 (Pre-Sweep Late Mar)", "2026-03-16", "2026-03-29"),
        ("Mar 1-29 (ALL Pre-Sweep Mar)", "2026-03-01", "2026-03-29"),
        ("Mar 30-31 (Email Days 1-2)", "2026-03-30", "2026-03-31"),
        ("Mar 30-Apr 12 (Private Email Only)", "2026-03-30", "2026-04-12"),
        ("Apr 13-30 (Post-Public Press Apr)", "2026-04-13", "2026-04-30"),
        ("May 1-15 (Post-Public Press May)", "2026-05-01", "2026-05-15"),
    ]
    print("\n=== SUB-PERIOD DAILY RATE DIAGNOSTIC (DIAGNOSING THE MARCH 2026 JUMP!) ===")
    sub_win_summary = []
    for lbl, s_d, e_d in sub_windows:
        g = ddf[(ddf["date"] >= s_d) & (ddf["date"] <= e_d)]
        nd = len(g)
        ind_od_ev = int(g["ind_od_events"].sum())
        ind_od_rate30 = 100.0 * (ind_od_ev / g["ind_od_risk"].sum()) * 30.4375
        acad_od_ev = int(g["acad_od_events"].sum())
        acad_od_rate30 = 100.0 * (acad_od_ev / g["acad_od_risk"].sum()) * 30.4375
        ind_ot_ev = int(g["ind_ot_events"].sum())
        ind_ot_rate30 = 100.0 * (ind_ot_ev / g["ind_ot_risk"].sum()) * 30.4375
        ratio_ind_acad = ind_od_rate30 / acad_od_rate30 if acad_od_rate30 > 0 else np.nan
        ratio_od_ot = ind_od_rate30 / ind_ot_rate30 if ind_ot_rate30 > 0 else np.nan
        print(f"  {lbl:36s} ({nd:2d}d) | Ind_OD={ind_od_ev:3d} ({ind_od_rate30:5.2f}%/30d) | Acad_OD={acad_od_ev:3d} ({acad_od_rate30:5.2f}%) | Ind_OT={ind_ot_ev:3d} ({ind_ot_rate30:5.2f}%) | Ind/Acad={ratio_ind_acad:.2f} | OD/OT={ratio_od_ot:.2f}")
        sub_win_summary.append({
            "window": lbl, "start": s_d, "end": e_d, "days": nd,
            "ind_od_events": ind_od_ev, "ind_od_rate30_pct": round(ind_od_rate30, 4),
            "acad_od_events": acad_od_ev, "acad_od_rate30_pct": round(acad_od_rate30, 4),
            "ind_ot_events": ind_ot_ev, "ind_ot_rate30_pct": round(ind_ot_rate30, 4),
            "ind_to_acad_od_ratio": round(float(ratio_ind_acad), 4) if pd.notna(ratio_ind_acad) else None,
            "ind_od_to_ot_ratio": round(float(ratio_od_ot), 4) if pd.notna(ratio_od_ot) else None,
        })

    # Formal weekly event-study Poisson regression (W-12..W+12, ref=W-1 = Mar 23-29!)
    # Check Private Window (W0..W+1 = Mar 30..Apr 12) vs Public Window (W+2..W+12) vs Pre-Sweep March (W-4..W-1)!
    print("\n=== FORMAL SUB-WINDOW POISSON DDD & IND OD/OT (SEPARATING MAR 1-29 vs MAR 30-APR 12 vs APR 13+) ===")
    # Construct trial-level biweekly/regime panel across 2025-04-01..2026-09-30 using exact date windows:
    # Pre-Baseline: 2025-04-01..2026-02-28 (11 monthly periods)
    # Pre-Sweep March: 2026-03-01..2026-03-29 (29 days!)
    # Private-Email Window: 2026-03-30..2026-04-12 (14 days!)
    # Public-Press Early: 2026-04-13..2026-05-31 (49 days!)
    # Public-Press Late: 2026-06-01..2026-09-30 (4 monthly periods)
    regime_periods = []
    for m in pd.period_range("2025-04", "2026-02", freq="M"):
        regime_periods.append((str(m), m.start_time, m.end_time, "PRE_BASE"))
    regime_periods.append(("2026-03_PRE", pd.Timestamp("2026-03-01"), pd.Timestamp("2026-03-29"), "PRE_MAR"))
    regime_periods.append(("2026-0330_PRIV", pd.Timestamp("2026-03-30"), pd.Timestamp("2026-04-12"), "PRIVATE_EMAIL"))
    regime_periods.append(("2026-0413_PUB1", pd.Timestamp("2026-04-13"), pd.Timestamp("2026-05-12"), "PUBLIC_PRESS"))
    regime_periods.append(("2026-0513_PUB2", pd.Timestamp("2026-05-13"), pd.Timestamp("2026-06-30"), "PUBLIC_PRESS"))
    for m in pd.period_range("2026-07", "2026-09", freq="M"):
        regime_periods.append((str(m), m.start_time, m.end_time, "PUBLIC_PRESS"))

    reg_recs = []
    for pname, ps, pe, regime in regime_periods:
        ndays = (pe - ps).days + 1
        ar = (act["pcd"] < ps) & (act[sub_col].isna() | (act[sub_col] >= ps))
        sub = act.loc[ar, ["lead_sponsor_clean", "is_industry", "pcd", "due_date", sub_col]].copy()
        sub["period"] = pname
        sub["regime"] = regime
        sub["is_overdue"] = (sub["due_date"] < ps).astype(int)
        sub["y"] = ((sub[sub_col] >= ps) & (sub[sub_col] <= pe)).astype(int)
        sub["exposure"] = ndays / 30.4375
        age_m = ((ps - sub["pcd"]).dt.days / 30.4375).clip(lower=0.0)
        sub["pcd_age_bin"] = pd.cut(
            age_m, bins=[-0.01, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 24.0, 36.0, 60.0, 999.0],
            labels=["0_3", "3_6", "6_9", "9_12", "12_15", "15_18", "18_24", "24_36", "36_60", "60p"]
        ).astype(str)
        reg_recs.append(sub)
    rpan = pd.concat(reg_recs, ignore_index=True)
    rcdf = rpan.groupby(["lead_sponsor_clean", "is_overdue", "is_industry", "pcd_age_bin", "period", "regime"], as_index=False).agg(
        events=("y", "sum"), risk_exp=("exposure", "sum")
    )
    rcdf["cell"] = rcdf["is_overdue"].astype(str) + "_" + rcdf["is_industry"].astype(str) + "_" + rcdf["pcd_age_bin"]
    rcdf = rcdf[rcdf.groupby("cell")["events"].transform("sum") > 0].copy()
    for rg in ["PRE_MAR", "PRIVATE_EMAIL", "PUBLIC_PRESS"]:
        rcdf[f"ddd_{rg}"] = ((rcdf["is_overdue"] == 1) & (rcdf["is_industry"] == 1) & (rcdf["regime"] == rg)).astype(float)
        rcdf[f"od_{rg}"] = ((rcdf["is_overdue"] == 1) & (rcdf["regime"] == rg)).astype(float)
        rcdf[f"in_{rg}"] = ((rcdf["is_industry"] == 1) & (rcdf["regime"] == rg)).astype(float)

    f_reg_ddd = smf.glm(
        "events ~ ddd_PRE_MAR + ddd_PRIVATE_EMAIL + ddd_PUBLIC_PRESS + od_PRE_MAR + od_PRIVATE_EMAIL + od_PUBLIC_PRESS + in_PRE_MAR + in_PRIVATE_EMAIL + in_PUBLIC_PRESS + C(cell) + C(period)",
        data=rcdf, family=sm.families.Poisson(), offset=np.log(rcdf["risk_exp"].values)
    ).fit(cov_type="cluster", cov_kwds={"groups": rcdf["lead_sponsor_clean"]})

    # Industry-only OD vs OT across the 3 exact regimes
    rci = rcdf[rcdf["is_industry"] == 1].copy()
    rci = rci[rci.groupby("cell")["events"].transform("sum") > 0].copy()
    f_reg_ind = smf.glm(
        "events ~ od_PRE_MAR + od_PRIVATE_EMAIL + od_PUBLIC_PRESS + C(cell) + C(period)",
        data=rci, family=sm.families.Poisson(), offset=np.log(rci["risk_exp"].values)
    ).fit(cov_type="cluster", cov_kwds={"groups": rci["lead_sponsor_clean"]})

    timing_tests = {"ot_od_ddd": {}, "ind_od_ot": {}}
    for rg in ["PRE_MAR", "PRIVATE_EMAIL", "PUBLIC_PRESS"]:
        b, se, pv = float(f_reg_ddd.params[f"ddd_{rg}"]), float(f_reg_ddd.bse[f"ddd_{rg}"]), float(f_reg_ddd.pvalues[f"ddd_{rg}"])
        timing_tests["ot_od_ddd"][rg] = {"irr": round(np.exp(b), 4), "ci95": [round(np.exp(b - 1.96 * se), 4), round(np.exp(b + 1.96 * se), 4)], "p": pv}
        bi, sei, pvi = float(f_reg_ind.params[f"od_{rg}"]), float(f_reg_ind.bse[f"od_{rg}"]), float(f_reg_ind.pvalues[f"od_{rg}"])
        timing_tests["ind_od_ot"][rg] = {"irr": round(np.exp(bi), 4), "ci95": [round(np.exp(bi - 1.96 * sei), 4), round(np.exp(bi + 1.96 * sei), 4)], "p": pvi}
        print(f"  {rg:15s} | OT/OD DDD: IRR={np.exp(b):.4f} [{np.exp(b-1.96*se):.4f}, {np.exp(b+1.96*se):.4f}] p={pv:.4g} | Ind OD/OT: IRR={np.exp(bi):.4f} [{np.exp(bi-1.96*sei):.4f}, {np.exp(bi+1.96*sei):.4f}] p={pvi:.4g}")

    # ==========================================================================
    # TEST 2B: ULTRA-TIGHT SPONSOR × MONTH FE, THERAPEUTIC-AREA × MONTH FE, PHASE×INTV×MONTH FE
    # ==========================================================================
    print("\n=== TEST 2B: HIGH-DIMENSIONAL FE CONFOUNDING TESTS (EXACT PPMLHDFE, SPONSOR-CLUSTERED) ===")
    # Enrich monthly panel with TA, phase, intervention
    meta_cols = act[["nct_id", "ta", "phase_clean", "intv_clean"]].drop_duplicates("nct_id")
    pan = pan_base.merge(meta_cols, on="nct_id", how="left")
    pan["sponsor_ym"] = pan["lead_sponsor_clean"] + "__" + pan["ym"]
    pan["ta_ym"] = pan["ta"] + "__" + pan["ym"]
    pan["phase_intv_ym"] = pan["phase_clean"] + "__" + pan["intv_clean"] + "__" + pan["ym"]

    # Check how many dual-status sponsor-months exist (sponsors holding BOTH overdue & on-time ACTs in same month!)
    sm_stat = pan.groupby("sponsor_ym").agg(has_od=("is_overdue", "max"), has_ot=("is_overdue", "min"), ev=("y", "sum"))
    n_dual_sm = int(((sm_stat["has_od"] == 1) & (sm_stat["has_ot"] == 0)).sum())
    n_dual_sm_ev = int(((sm_stat["has_od"] == 1) & (sm_stat["has_ot"] == 0) & (sm_stat["ev"] > 0)).sum())
    print(f"  Dual-status Sponsor×Month cells (both OD & OT trials in same firm-month): {n_dual_sm:,} total, {n_dual_sm_ev:,} with >=1 submission")

    fe_results = {}

    # Helper to run OT/OD DDD and Ind OD/OT under any absorbed FE column!
    def run_fe_spec(pan_in, fe_col, label):
        # 1. Full OT/OD DDD with absorbed fe_col
        # Note: when fe_col == "sponsor_ym", is_industry × ym is COMPLETELY ABSORBED by Sponsor × Month FE!
        agg = pan_in.groupby(["lead_sponsor_clean", fe_col, "is_overdue", "is_industry", "pcd_age_bin", "post"], as_index=False).agg(
            events=("y", "sum"), risk=("y", "count")
        )
        agg["ddd"] = (agg["is_overdue"] * agg["is_industry"] * agg["post"]).astype(float)
        agg["od_post"] = (agg["is_overdue"] * agg["post"]).astype(float)
        agg["ind_post"] = (agg["is_industry"] * agg["post"]).astype(float)
        agg["cell"] = agg["is_overdue"].astype(str) + "_" + agg["is_industry"].astype(str) + "_" + agg["pcd_age_bin"]
        # Drop cells with zero events
        agg = agg[agg.groupby("cell")["events"].transform("sum") > 0].copy().reset_index(drop=True)
        cd = pd.get_dummies(agg["cell"], prefix="c", drop_first=True, dtype=float)
        agg_x = pd.concat([agg, cd], axis=1)
        if fe_col == "sponsor_ym":
            # ind_post is collinear with sponsor_ym (since each sponsor is either industry or non-industry)!
            x_ddd = ["ddd", "od_post"] + list(cd.columns)
        else:
            x_ddd = ["ddd", "od_post", "ind_post"] + list(cd.columns)
        r_ddd = fit_ppml_fe(agg_x, "events", "risk", x_ddd, fe_col, "lead_sponsor_clean")

        # 2. Industry-Only OD vs OT under absorbed fe_col
        ai = agg[agg["is_industry"] == 1].copy().reset_index(drop=True)
        ai["cell_i"] = ai["is_overdue"].astype(str) + "_" + ai["pcd_age_bin"]
        ai = ai[ai.groupby("cell_i")["events"].transform("sum") > 0].copy().reset_index(drop=True)
        cdi = pd.get_dummies(ai["cell_i"], prefix="ci", drop_first=True, dtype=float)
        ai_x = pd.concat([ai, cdi], axis=1)
        r_ind = fit_ppml_fe(ai_x, "events", "risk", ["od_post"] + list(cdi.columns), fe_col, "lead_sponsor_clean")

        # 3. Academic-Only OD vs OT under absorbed fe_col
        aa = agg[agg["is_industry"] == 0].copy().reset_index(drop=True)
        aa["cell_a"] = aa["is_overdue"].astype(str) + "_" + aa["pcd_age_bin"]
        aa = aa[aa.groupby("cell_a")["events"].transform("sum") > 0].copy().reset_index(drop=True)
        cda = pd.get_dummies(aa["cell_a"], prefix="ca", drop_first=True, dtype=float)
        aa_x = pd.concat([aa, cda], axis=1)
        r_acad = fit_ppml_fe(aa_x, "events", "risk", ["od_post"] + list(cda.columns), fe_col, "lead_sponsor_clean")

        out = {
            "fe_col": fe_col,
            "ot_od_ddd": {
                "irr": round(r_ddd["irr"]["ddd"], 4), "ci95": [round(x, 4) for x in r_ddd["ci95"]["ddd"]],
                "p": r_ddd["p"]["ddd"], "n_events": r_ddd["n_events_active"], "n_fe": r_ddd["n_fe_groups_active"], "n_clusters": r_ddd["n_clusters_active"],
            },
            "ind_od_ot": {
                "irr": round(r_ind["irr"]["od_post"], 4), "ci95": [round(x, 4) for x in r_ind["ci95"]["od_post"]],
                "p": r_ind["p"]["od_post"], "n_events": r_ind["n_events_active"], "n_fe": r_ind["n_fe_groups_active"], "n_clusters": r_ind["n_clusters_active"],
            },
            "acad_od_ot": {
                "irr": round(r_acad["irr"]["od_post"], 4), "ci95": [round(x, 4) for x in r_acad["ci95"]["od_post"]],
                "p": r_acad["p"]["od_post"], "n_events": r_acad["n_events_active"], "n_fe": r_acad["n_fe_groups_active"], "n_clusters": r_acad["n_clusters_active"],
            },
        }
        print(
            f"  {label:30s} | DDD={out['ot_od_ddd']['irr']:.4f} {out['ot_od_ddd']['ci95']} (p={out['ot_od_ddd']['p']:.4g}) | "
            f"Ind OD/OT={out['ind_od_ot']['irr']:.4f} {out['ind_od_ot']['ci95']} (p={out['ind_od_ot']['p']:.4g}) | "
            f"Acad OD/OT={out['acad_od_ot']['irr']:.4f} {out['acad_od_ot']['ci95']} (p={out['acad_od_ot']['p']:.4g})"
        )
        return out

    fe_results["baseline_month_fe"] = run_fe_spec(pan, "ym", "Baseline Month FE")
    fe_results["sponsor_x_month_fe"] = run_fe_spec(pan, "sponsor_ym", "Sponsor × Month FE")
    fe_results["ta_x_month_fe"] = run_fe_spec(pan, "ta_ym", "Therapeutic-Area × Month FE")
    fe_results["phase_intv_x_month_fe"] = run_fe_spec(pan, "phase_intv_ym", "Phase × Intv × Month FE")

    # ==========================================================================
    # TEST 2C: ON-TIME MECHANICAL CLIFF EXCLUSION, FINE 1M AGE FE, DONUT-MARCH, LOSO
    # ==========================================================================
    print("\n=== TEST 2C: ON-TIME MECHANICAL CLIFF (0-9m / 0-6m), DONUT-MARCH, FINE AGE, LOSO ===")
    sens_results = {}

    # 1. Cliff-free On-Time: drop 9_12m On-Time bin (keep 0-9m On-Time + all Overdue!)
    pan_09 = pan[pan["pcd_age_bin"] != "9_12"].copy()
    sens_results["cliff_free_ot_0_9m_month_fe"] = run_fe_spec(pan_09, "ym", "Cliff-Free OT (0-9m) Month FE")
    sens_results["cliff_free_ot_0_9m_sponsor_ym_fe"] = run_fe_spec(pan_09, "sponsor_ym", "Cliff-Free OT (0-9m) Spon×Mo")

    # 2. Stricter Cliff-free On-Time: drop 6_9m and 9_12m (keep 0-6m On-Time + all Overdue!)
    pan_06 = pan[~pan["pcd_age_bin"].isin(["6_9", "9_12"])].copy()
    sens_results["cliff_free_ot_0_6m_month_fe"] = run_fe_spec(pan_06, "ym", "Cliff-Free OT (0-6m) Month FE")

    # 3. Donut-March (drop 2026-03 entirely!)
    pan_donut = pan[pan["ym"] != "2026-03"].copy()
    sens_results["donut_march_month_fe"] = run_fe_spec(pan_donut, "ym", "Donut-March (Drop 2026-03)")
    sens_results["donut_march_sponsor_ym_fe"] = run_fe_spec(pan_donut, "sponsor_ym", "Donut-March Sponsor×Month FE")

    # 4. Fine 1-month PCD-age bins (0,1,2,...,35, 36-47, 48-59, 60+)!
    pan_fine = pan.copy()
    # Recompute exact age_m from ym start - pcd
    ym_start = pd.to_datetime(pan_fine["ym"] + "-01")
    age_m_int = np.floor(((ym_start - pan_fine["pcd"]).dt.days / 30.4375).clip(lower=0.0)).astype(int)
    pan_fine["pcd_age_bin"] = np.where(
        age_m_int < 36, age_m_int.astype(str) + "m",
        np.where(age_m_int < 48, "36_47m", np.where(age_m_int < 60, "48_59m", "60p"))
    )
    sens_results["fine_1month_age_fe"] = run_fe_spec(pan_fine, "ym", "Fine 1-Month Age FE (0..35m+)")

    # 5. Leave-One-Sponsor-Out (Top 20 Industry & Academic sponsors by overdue risk-months!)
    top20_sponsors = (
        pan[pan["is_overdue"] == 1].groupby("lead_sponsor_clean").size()
        .sort_values(ascending=False).head(20).index.tolist()
    )
    loso_rows = []
    for sp in top20_sponsors:
        sub_p = pan[pan["lead_sponsor_clean"] != sp]
        agg = sub_p.groupby(["lead_sponsor_clean", "ym", "is_overdue", "is_industry", "pcd_age_bin", "post"], as_index=False).agg(
            events=("y", "sum"), risk=("y", "count")
        )
        agg["ddd"] = (agg["is_overdue"] * agg["is_industry"] * agg["post"]).astype(float)
        agg["od_post"] = (agg["is_overdue"] * agg["post"]).astype(float)
        agg["ind_post"] = (agg["is_industry"] * agg["post"]).astype(float)
        agg["cell"] = agg["is_overdue"].astype(str) + "_" + agg["is_industry"].astype(str) + "_" + agg["pcd_age_bin"]
        agg = agg[agg.groupby("cell")["events"].transform("sum") > 0].copy().reset_index(drop=True)
        cd = pd.get_dummies(agg["cell"], prefix="c", drop_first=True, dtype=float)
        r = fit_ppml_fe(pd.concat([agg, cd], axis=1), "events", "risk", ["ddd", "od_post", "ind_post"] + list(cd.columns), "ym", "lead_sponsor_clean")
        loso_rows.append({
            "dropped_sponsor": sp,
            "ddd_irr": round(r["irr"]["ddd"], 4),
            "ddd_p": round(r["p"]["ddd"], 6),
        })
    loso_df = pd.DataFrame(loso_rows)
    loso_df.to_csv(os.path.join(OUT_DIR, "b2_attack2_leave_one_sponsor_out_top20.csv"), index=False)
    print(f"\n  Leave-One-Sponsor-Out (Top 20 overdue sponsors): DDD IRR range = [{loso_df['ddd_irr'].min():.4f}, {loso_df['ddd_irr'].max():.4f}], max p = {loso_df['ddd_p'].max():.4g}")
    print(loso_df.head(8).to_string(index=False))

    # Save complete Attack 2 JSON!
    full_out = {
        "classification": "POST_HOC_EXPLORATORY_FALSIFICATION",
        "locked_in": "NOVELTY_IDENTIFICATION_GATE.md (2026-11-25)",
        "dual_status_sponsor_months_total": n_dual_sm,
        "dual_status_sponsor_months_with_events": n_dual_sm_ev,
        "test_2a_sub_window_daily_diagnostics": sub_win_summary,
        "test_2a_formal_regime_regressions": timing_tests,
        "test_2b_high_dim_fe_results": fe_results,
        "test_2c_cliff_donut_fineage_results": sens_results,
        "test_2c_loso_top20_summary": {
            "min_ddd_irr": float(loso_df["ddd_irr"].min()),
            "max_ddd_irr": float(loso_df["ddd_irr"].max()),
            "max_ddd_p": float(loso_df["ddd_p"].max()),
            "all_p_lt_005": bool((loso_df["ddd_p"] < 0.05).all()),
        },
    }
    out_json = os.path.join(OUT_DIR, "b2_attack2_confounding_falsification.json")
    with open(out_json, "w") as fh:
        json.dump(full_out, fh, indent=2)
    print(f"\nDONE! Saved {out_json}")


if __name__ == "__main__":
    main()
