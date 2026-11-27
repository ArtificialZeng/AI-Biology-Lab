#!/opt/anaconda3/bin/python3
"""
Generate 4 Publication-Grade Multi-Panel Figures for Study B-2 (Management Science / JFE style).
Every plotted point is read directly from verified econometric output files:
  - b2_v2_stationary_results.json
  - b2_v2_censoring_pretrend_audit.json
  - b2_v2_stationary_monthly_hazards.parquet
  - b2_trial_universe_steady.parquet (for exact RI permutation histogram)
"""
import json, os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

warnings.filterwarnings("ignore")

BASE = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_b2_fda_2026_sweep"
RES = os.path.join(BASE, "results")
FIG_DIR = os.path.join(BASE, "figures")
ART_DIR = "/Users/mac/.gemini/antigravity/brain/19f7514f-ce6a-4297-92cb-9803e648afd1"
os.makedirs(FIG_DIR, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9.5,
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
    "axes.labelsize": 9.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.9,
    "xtick.labelsize": 8.6,
    "ytick.labelsize": 8.6,
    "legend.fontsize": 8.5,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.facecolor": "white",
})

# Colorblind-safe Okabe-Ito / Nature palette
C_IND = "#0072B2"    # Deep Blue (Industry)
C_ACA = "#D55E00"    # Vermilion (Academic)
C_NIH = "#009E73"    # Teal Green (NIH/Fed)
C_CTRL = "#7F8C8D"   # Slate Gray (Exempt Controls)
C_PLAC = "#95A5A6"   # Light Gray (Placebos)
C_ACC = "#CC79A7"    # Muted Purple (Accent)
C_GOLD = "#E69F00"   # Amber Gold

v2 = json.load(open(os.path.join(RES, "b2_v2_stationary_results.json")))
aud = json.load(open(os.path.join(RES, "b2_v2_censoring_pretrend_audit.json")))
mh = pd.read_parquet(os.path.join(RES, "b2_v2_stationary_monthly_hazards.parquet"))


def save_both(fig, name):
    p1 = os.path.join(FIG_DIR, name)
    p2 = os.path.join(ART_DIR, name)
    fig.savefig(p1, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(p2, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved {p1} ({os.path.getsize(p1):,} B)")


# ==============================================================================
# FIGURE 1: Monthly Stationary Hazards & Dynamic Event-Study (Industry vs Academic)
# ==============================================================================
def make_fig1():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.0, 5.4), gridspec_kw={"width_ratios": [1.05, 1.05]})
    fig.subplots_adjust(wspace=0.28)

    # Panel A: Monthly submission hazard (% of delinquent stock per month)
    piv = mh.groupby(["ym", "arm", "is_industry"], as_index=False).agg(ev=("events", "sum"), rk=("risk", "sum"))
    piv["h"] = piv["ev"] / piv["rk"] * 100.0
    months = sorted(piv["ym"].unique())
    x = np.arange(len(months))

    def get_s(arm_list, ind_val=None):
        s = piv[piv["arm"].isin(arm_list)]
        if ind_val is not None:
            s = s[s["is_industry"] == ind_val]
        g = s.groupby("ym").agg(ev=("ev", "sum"), rk=("rk", "sum"))
        return np.array([g.loc[m, "ev"] / g.loc[m, "rk"] * 100.0 if m in g.index else np.nan for m in months])

    h_ind = get_s(["us_fdareg_nonp1"], 1)
    h_aca = get_s(["us_fdareg_nonp1"], 0)
    h_ctrl = get_s(["us_notfdareg", "nonus_p23_nonfda"], None)

    i_feb = months.index("2026-02")
    i_mar = months.index("2026-03")
    i_apr = months.index("2026-04")
    i_jul = months.index("2026-07")
    pct_surge = (h_ind[i_apr] / h_ind[i_feb] - 1.0) * 100.0

    ax1.axvspan(i_apr - 0.35, i_jul + 0.45, color="#EBF5FB", alpha=0.85, zorder=0, label="Primary Post-Sweep (Apr–Jul 2026, ≥117d buffer)")
    ax1.axvspan(i_jul + 0.45, len(months) - 0.5, color="#F4F6F7", alpha=0.85, zorder=0, label="Partial NLM QC Lag Zone (Aug–Sep 2026)")
    ax1.axvline(i_mar + 0.45, color="#2C3E50", linestyle="--", linewidth=1.4, zorder=2)

    ax1.plot(x, h_ind, marker="o", lw=2.3, ms=5.5, color=C_IND, label="FDA-Regulated ACT: Industry (Corporate)", zorder=4)
    ax1.plot(x, h_aca, marker="s", lw=2.1, ms=5.2, color=C_ACA, linestyle="-.", label="FDA-Regulated ACT: Academic / Non-Industry", zorder=3)
    ax1.plot(x, h_ctrl, marker="^", lw=1.7, ms=4.5, color=C_CTRL, linestyle=":", label="Exempt Controls (US Non-FDA + Non-US P2/3)", zorder=3)

    ax1.annotate("FDA Mass Sweep\n(Mar 30, 2026)", xy=(i_mar + 0.45, 5.35), xytext=(i_mar - 4.4, 5.50),
                 fontsize=8.3, fontweight="bold", color="#2C3E50",
                 arrowprops=dict(arrowstyle="->", color="#2C3E50", lw=1.1))
    ax1.annotate(f"Industry surges\n+{pct_surge:.0f}% vs Feb\n({h_ind[i_apr]:.2f}% in Apr)", xy=(i_apr, h_ind[i_apr]), xytext=(i_apr + 1.0, 4.45),
                 fontsize=8.0, color=C_IND, fontweight="bold",
                 arrowprops=dict(arrowstyle="->", color=C_IND, lw=1.0))

    ax1.set_xticks(x[::2])
    ax1.set_xticklabels([months[i] for i in range(0, len(months), 2)], rotation=35, ha="right")
    ax1.set_ylabel("Monthly Delinquent Submission Hazard (%)")
    ax1.set_ylim(-0.15, 6.25)
    ax1.set_title("A. Stationary Monthly Delinquent Hazard (N = 53,927)", pad=10)
    ax1.grid(axis="y", alpha=0.28, linestyle="-")
    ax1.legend(loc="upper center", bbox_to_anchor=(0.48, -0.19), ncol=2, frameon=True, facecolor="white", fontsize=7.8)

    # Panel B: Dynamic Event-Study (Poisson log-IRR normalized to pre-sweep mean 2025-04..2026-02)
    es = v2["event_study_monthly_within_act"]["series"]
    ym_es = [r["ym"] for r in es]
    pre_logs = [r["log_irr"] for r in es if r["ym"] <= "2026-02"]
    c_shift = float(np.mean(pre_logs))
    xe = np.arange(len(ym_es))
    vals = np.array([r["log_irr"] - c_shift for r in es])
    ses = np.array([r["se"] for r in es])
    irrs = np.exp(vals)
    lo = np.exp(vals - 1.96 * ses)
    hi = np.exp(vals + 1.96 * ses)

    ax2.axvspan(i_apr - 0.35, i_jul + 0.45, color="#EBF5FB", alpha=0.85, zorder=0)
    ax2.axvspan(i_jul + 0.45, len(ym_es) - 0.5, color="#F4F6F7", alpha=0.85, zorder=0)
    ax2.axhline(1.0, color="#7F8C8D", linestyle="--", lw=1.1, zorder=1)
    ax2.axvline(i_mar + 0.45, color="#2C3E50", linestyle="--", lw=1.4, zorder=2)

    for idx_m in range(len(ym_es)):
        is_post = ym_es[idx_m] >= "2026-04"
        is_mar = ym_es[idx_m] == "2026-03"
        col = C_IND if is_post else (C_GOLD if is_mar else "#566573")
        mk = "o" if is_post else ("D" if is_mar else "o")
        if ses[idx_m] > 0:
            ax2.errorbar(xe[idx_m], irrs[idx_m], yerr=[[irrs[idx_m] - lo[idx_m]], [hi[idx_m] - irrs[idx_m]]],
                         fmt=mk, color=col, ecolor=col, elinewidth=1.3, capsize=2.8, ms=5.8, zorder=4)
        else:
            ax2.plot(xe[idx_m], irrs[idx_m], marker="o", color=col, ms=5.8, zorder=4)

    ax2.set_yscale("log")
    ax2.set_yticks([0.25, 0.5, 1.0, 2.0, 4.0])
    ax2.set_yticklabels(["0.25×", "0.50×", "1.0×", "2.0×", "4.0×"])
    ax2.set_ylim(0.16, 6.4)
    ax2.set_xticks(xe[::2])
    ax2.set_xticklabels([ym_es[i] for i in range(0, len(ym_es), 2)], rotation=35, ha="right")
    ax2.set_ylabel("Industry / Academic Relative Hazard Ratio (Log Scale)")
    ax2.set_title("B. Within-ACT Poisson Event Study (Pre-Mean = 1.0)", pad=10)
    ax2.grid(axis="y", alpha=0.28)

    txt = (
        "Primary Post-Sweep (Apr–Jul 2026):\n"
        "  Poisson IRR = 1.564 (95% CI: 1.15–2.12, p = 0.0041)\n"
        "  Linear Pre-Trend Slope = −0.049/mo (p = 0.029)\n"
        "  (Negative pre-trend biases against post lift)"
    )
    ax2.text(0.03, 0.03, txt, transform=ax2.transAxes, fontsize=7.8,
             bbox=dict(boxstyle="round,pad=0.38", facecolor="white", edgecolor="#BDC3C7", alpha=0.95))
    ax2.annotate("Nov 2025\nAcad. dip", xy=(ym_es.index("2025-11"), irrs[ym_es.index("2025-11")]),
                 xytext=(ym_es.index("2025-11") - 3.3, 4.35), fontsize=7.6, color="#566573",
                 arrowprops=dict(arrowstyle="->", color="#7F8C8D", lw=0.9))

    save_both(fig, "fig1_b2_event_study_and_hazards.png")


# ==============================================================================
# FIGURE 2: Specification Curve, Zero-Censoring Audits, Placebos & Randomization Inference
# ==============================================================================
def make_fig2():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.6, 5.5), gridspec_kw={"width_ratios": [1.32, 0.85]})
    fig.subplots_adjust(wspace=0.28)

    wa = v2["within_act_governance_did"]
    uh = aud["uniform_horizon_zero_censoring"]
    rows = [
        ("Primary Apr–Jul 2026 (≥117d QC buffer)", wa["primary_apr_jul_2026"]["poisson_age_fe"]["irr"],
         wa["primary_apr_jul_2026"]["poisson_age_fe"]["irr_ci95"], wa["primary_apr_jul_2026"]["poisson_age_fe"]["p"], "main"),
        ("Ultra-Clean Q2 Apr–Jun 2026 (≥148d buffer)", wa["ultraclean_apr_jun_2026"]["poisson_age_fe"]["irr"],
         wa["ultraclean_apr_jun_2026"]["poisson_age_fe"]["irr_ci95"], wa["ultraclean_apr_jun_2026"]["poisson_age_fe"]["p"], "main"),
        ("Full 6-Month Post Window (Apr–Sep 2026)", wa["full_apr_sep_2026"]["poisson_age_fe"]["irr"],
         wa["full_apr_sep_2026"]["poisson_age_fe"]["irr_ci95"], wa["full_apr_sep_2026"]["poisson_age_fe"]["p"], "main"),
        ("Donut Spec (Excluding March 2026)", wa["donut_drop_mar2026"]["poisson_age_fe"]["irr"],
         wa["donut_drop_mar2026"]["poisson_age_fe"]["irr_ci95"], wa["donut_drop_mar2026"]["poisson_age_fe"]["p"], "main"),
        ("Zero-Censoring Horizon H=120d (Apr–Jul)", uh["H120_apr_jul_2026"]["within_act_irr"],
         uh["H120_apr_jul_2026"]["within_act_irr_ci95"], uh["H120_apr_jul_2026"]["within_act_irr_p"], "censor"),
        ("Zero-Censoring Horizon H=90d (Apr–Aug)", uh["H90_apr_aug_2026"]["within_act_irr"],
         uh["H90_apr_aug_2026"]["within_act_irr_ci95"], uh["H90_apr_aug_2026"]["within_act_irr_p"], "censor"),
        ("IPCW QC-CDF Adjusted (Apr–Jul 2026)", aud["ipcw_apr_jul_2026"]["ipcw_irr"],
         aud["ipcw_apr_jul_2026"]["ipcw_ci95"], aud["ipcw_apr_jul_2026"]["ipcw_p"], "censor"),
        ("IPCW QC-CDF Adjusted (Apr–Sep 2026)", aud["ipcw_apr_sep_2026"]["ipcw_irr"],
         aud["ipcw_apr_sep_2026"]["ipcw_ci95"], aud["ipcw_apr_sep_2026"]["ipcw_p"], "censor"),
        ("Placebo Sweep: March 30, 2025 (All)", wa["placebo_2025_apr_jul"]["poisson_age_fe"]["irr"],
         wa["placebo_2025_apr_jul"]["poisson_age_fe"]["irr_ci95"], wa["placebo_2025_apr_jul"]["poisson_age_fe"]["p"], "placebo"),
        ("Placebo Sweep: March 30, 2024 (All)", wa["placebo_2024_apr_jul"]["poisson_age_fe"]["irr"],
         wa["placebo_2024_apr_jul"]["poisson_age_fe"]["irr_ci95"], wa["placebo_2024_apr_jul"]["poisson_age_fe"]["p"], "placebo"),
        ("Placebo 2025 (Zero-Censoring H=120d)", uh["H120_placebo_2025"]["within_act_irr"],
         uh["H120_placebo_2025"]["within_act_irr_ci95"], uh["H120_placebo_2025"]["within_act_irr_p"], "placebo"),
        ("Placebo 2024 (Zero-Censoring H=120d)", uh["H120_placebo_2024"]["within_act_irr"],
         uh["H120_placebo_2024"]["within_act_irr_ci95"], uh["H120_placebo_2024"]["within_act_irr_p"], "placebo"),
    ]
    y_pos = np.arange(len(rows))[::-1]
    ax1.axvline(1.0, color="#2C3E50", linestyle="--", lw=1.1, zorder=1)
    ax1.axhline(3.5, color="#BDC3C7", linestyle="-", lw=0.9)

    col_map = {"main": C_IND, "censor": C_NIH, "placebo": C_ACA}
    for i, (lbl, irr, ci, pval, grp) in enumerate(rows):
        y = y_pos[i]
        c = col_map[grp]
        mk = "o" if grp == "main" else ("s" if grp == "censor" else "X")
        ax1.errorbar(irr, y, xerr=[[irr - ci[0]], [ci[1] - irr]], fmt=mk, color=c, ecolor=c,
                     elinewidth=1.6, capsize=3.2, ms=6.2, zorder=4)
        p_txt = f"p={pval:.4f}" if pval >= 0.001 else f"p={pval:.1e}"
        ax1.text(ci[1] + 0.05, y, f"{irr:.2f}× ({p_txt})", va="center", ha="left", fontsize=7.8, color=c, fontweight="bold" if grp != "placebo" else "normal")

    ax1.set_yticks(y_pos)
    ax1.set_yticklabels([r[0] for r in rows], fontsize=8.3)
    ax1.set_xlim(0.50, 3.35)
    ax1.set_xlabel("Industry vs. Academic Post-Sweep Hazard Ratio (Poisson IRR, 95% CI)")
    ax1.set_title("A. Specification Robustness, Zero-Censoring Audits & Placebos")
    ax1.grid(axis="x", alpha=0.25)

    handles = [
        mpatches.Patch(color=C_IND, label="Primary & Window Sensitivity (2026)"),
        mpatches.Patch(color=C_NIH, label="Zero-Censoring Horizon & IPCW (2026)"),
        mpatches.Patch(color=C_ACA, label="Historical Placebo Sweeps (2024 & 2025)"),
    ]
    ax1.legend(handles=handles, loc="lower right", frameon=True, facecolor="white", fontsize=7.7)

    # Panel B: Recompute exact 500-perm RI distribution (deterministic seed=20260330 matching analyze_b2_v2_stationary.py!)
    import sys
    sys.path.insert(0, os.path.join(BASE, "src"))
    from analyze_b2_v2_stationary import load_steady_universe, build_stationary_monthly_panel
    df = load_steady_universe()
    p_jul = build_stationary_monthly_panel(df, "2025-04", "2026-07", "2026-04")
    act = p_jul[p_jul["arm"] == "us_fdareg_nonp1"]
    sp = act.groupby(["lead_sponsor_clean", "is_industry", "post"], as_index=False).agg(events=("y", "sum"), risk=("y", "count"))
    def calc_lirr(d_s):
        g = d_s.groupby(["is_industry", "post"])[["events", "risk"]].sum()
        h = g["events"] / g["risk"]
        return float(np.log((h.loc[(1, 1)] / h.loc[(1, 0)]) / (h.loc[(0, 1)] / h.loc[(0, 0)])))
    obs = calc_lirr(sp)
    rng = np.random.default_rng(20260330)
    smap = sp[["lead_sponsor_clean", "is_industry"]].drop_duplicates("lead_sponsor_clean").copy()
    perms = []
    for _ in range(500):
        pm = smap.copy()
        pm["ind_p"] = rng.permutation(pm["is_industry"].values)
        m = sp.drop(columns=["is_industry"]).merge(pm[["lead_sponsor_clean", "ind_p"]], on="lead_sponsor_clean").rename(columns={"ind_p": "is_industry"})
        perms.append(calc_lirr(m))
    perms = np.array(perms)
    ri_p = float((np.abs(perms) >= abs(obs)).mean())

    ax2.hist(perms, bins=26, color="#BDC3C7", edgecolor="white", alpha=0.9, density=True, label="Null Permutations\n(500 sponsor shuffles)")
    ax2.axvline(obs, color=C_IND, lw=2.3, linestyle="-", label=f"Observed = {obs:+.3f}\n(IRR = {np.exp(obs):.2f}×)")
    ax2.axvline(-obs, color=C_IND, lw=1.2, linestyle=":", alpha=0.6)
    ax2.annotate(f"Exact RI p = {ri_p:.3f}\n(z = {(obs - perms.mean())/perms.std():.2f})",
                 xy=(obs, 1.35), xytext=(0.14, 1.98), fontsize=8.4, fontweight="bold", color=C_IND,
                 arrowprops=dict(arrowstyle="->", color=C_IND, lw=1.2),
                 bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=C_IND, alpha=0.95))
    ax2.set_xlabel("Permuted Within-ACT Governance Log-IRR")
    ax2.set_ylabel("Density")
    ax2.set_title("B. Sponsor-Cluster Randomization Inference")
    ax2.set_xlim(-0.72, 0.72)
    ax2.set_ylim(0, 2.75)
    ax2.legend(loc="upper left", frameon=True, facecolor="white", fontsize=8.0)
    ax2.grid(axis="y", alpha=0.25)

    save_both(fig, "fig2_b2_identification_and_placebos.png")


# ==============================================================================
# FIGURE 3: Governance Mechanisms — Hierarchy, Backlog, Stale-Data Decay & Spillover
# ==============================================================================
def make_fig3():
    fig, axes = plt.subplots(1, 3, figsize=(14.4, 4.7))
    fig.subplots_adjust(wspace=0.30)
    het = v2["heterogeneity_stationary"]
    wa_3 = v2["within_act_governance_did"]["primary_apr_jul_2026"]["poisson_3group_vs_academic"]

    # Panel A: Sponsor Backlog Gradient (Industry vs Academic Poisson DiD IRR vs Clean Exempt)
    ax = axes[0]
    cats = ["Backlog = 1\nTrial", "Backlog = 2–4\nTrials", "Backlog = 5+\nTrials"]
    keys = ["backlog_1", "backlog_2_4", "backlog_5plus"]
    x = np.arange(len(cats))
    w = 0.34
    ind_v = [het["by_backlog"]["industry"][k]["irr"] for k in keys]
    ind_lo = [het["by_backlog"]["industry"][k]["irr"] - het["by_backlog"]["industry"][k]["ci95"][0] for k in keys]
    ind_hi = [het["by_backlog"]["industry"][k]["ci95"][1] - het["by_backlog"]["industry"][k]["irr"] for k in keys]
    aca_v = [het["by_backlog"]["academic"][k]["irr"] for k in keys]
    aca_lo = [het["by_backlog"]["academic"][k]["irr"] - het["by_backlog"]["academic"][k]["ci95"][0] for k in keys]
    aca_hi = [het["by_backlog"]["academic"][k]["ci95"][1] - het["by_backlog"]["academic"][k]["irr"] for k in keys]

    ax.axhline(1.0, color="#7F8C8D", ls="--", lw=1.0)
    ax.bar(x - w/2, ind_v, w, color=C_IND, alpha=0.88, label="Industry (Corporate)", yerr=[ind_lo, ind_hi], capsize=3.5, error_kw={"lw": 1.2})
    ax.bar(x + w/2, aca_v, w, color=C_ACA, alpha=0.82, label="Academic / Non-Industry", yerr=[aca_lo, aca_hi], capsize=3.5, error_kw={"lw": 1.2})
    for i in range(len(cats)):
        ax.text(x[i] - w/2, ind_v[i] + ind_hi[i] + 0.14, f"{ind_v[i]:.2f}×", ha="center", va="bottom", fontsize=7.8, fontweight="bold", color=C_IND)
        ax.text(x[i] + w/2, aca_v[i] + aca_hi[i] + 0.14, f"{aca_v[i]:.2f}×", ha="center", va="bottom", fontsize=7.8, color=C_ACA)
    ax.set_xticks(x)
    ax.set_xticklabels(cats)
    ax.set_ylabel("Post-Sweep Hazard Lift (Poisson IRR vs. Exempt)")
    ax.set_ylim(0, 8.8)
    ax.set_title("A. Responsiveness by Delinquent Backlog")
    ax.legend(loc="upper right", frameon=True, facecolor="white", fontsize=8.0)
    ax.grid(axis="y", alpha=0.25)

    # Panel B: Stale-Data Decay Curve (Cohort Freshness) + Large-Sponsor Governance Divergence
    ax = axes[1]
    c_lbls = ["Recent\n(2023–25)", "Medium\n(2021–22)", "Deeply Stale\n(2017–20)", "Large Org\n(20+ Trials)"]
    ind_vals = [
        het["by_cohort"]["industry"]["recent_2023_2025"]["irr"],
        het["by_cohort"]["industry"]["mid_2021_2022"]["irr"],
        het["by_cohort"]["industry"]["deep_2017_2020"]["irr"],
        het["by_portfolio"]["industry"]["large_20plus"]["irr"],
    ]
    ind_err = [
        [het["by_cohort"]["industry"]["recent_2023_2025"]["irr"] - het["by_cohort"]["industry"]["recent_2023_2025"]["ci95"][0],
         het["by_cohort"]["industry"]["mid_2021_2022"]["irr"] - het["by_cohort"]["industry"]["mid_2021_2022"]["ci95"][0],
         het["by_cohort"]["industry"]["deep_2017_2020"]["irr"] - het["by_cohort"]["industry"]["deep_2017_2020"]["ci95"][0],
         het["by_portfolio"]["industry"]["large_20plus"]["irr"] - het["by_portfolio"]["industry"]["large_20plus"]["ci95"][0]],
        [het["by_cohort"]["industry"]["recent_2023_2025"]["ci95"][1] - het["by_cohort"]["industry"]["recent_2023_2025"]["irr"],
         het["by_cohort"]["industry"]["mid_2021_2022"]["ci95"][1] - het["by_cohort"]["industry"]["mid_2021_2022"]["irr"],
         het["by_cohort"]["industry"]["deep_2017_2020"]["ci95"][1] - het["by_cohort"]["industry"]["deep_2017_2020"]["irr"],
         het["by_portfolio"]["industry"]["large_20plus"]["ci95"][1] - het["by_portfolio"]["industry"]["large_20plus"]["irr"]],
    ]
    aca_vals = [
        het["by_cohort"]["academic"]["recent_2023_2025"]["irr"],
        het["by_cohort"]["academic"]["mid_2021_2022"]["irr"],
        het["by_cohort"]["academic"]["deep_2017_2020"]["irr"],
        het["by_portfolio"]["academic"]["large_20plus"]["irr"],
    ]
    aca_err = [
        [het["by_cohort"]["academic"]["recent_2023_2025"]["irr"] - het["by_cohort"]["academic"]["recent_2023_2025"]["ci95"][0],
         het["by_cohort"]["academic"]["mid_2021_2022"]["irr"] - het["by_cohort"]["academic"]["mid_2021_2022"]["ci95"][0],
         het["by_cohort"]["academic"]["deep_2017_2020"]["irr"] - het["by_cohort"]["academic"]["deep_2017_2020"]["ci95"][0],
         het["by_portfolio"]["academic"]["large_20plus"]["irr"] - het["by_portfolio"]["academic"]["large_20plus"]["ci95"][0]],
        [het["by_cohort"]["academic"]["recent_2023_2025"]["ci95"][1] - het["by_cohort"]["academic"]["recent_2023_2025"]["irr"],
         het["by_cohort"]["academic"]["mid_2021_2022"]["ci95"][1] - het["by_cohort"]["academic"]["mid_2021_2022"]["irr"],
         het["by_cohort"]["academic"]["deep_2017_2020"]["ci95"][1] - het["by_cohort"]["academic"]["deep_2017_2020"]["irr"],
         het["by_portfolio"]["academic"]["large_20plus"]["ci95"][1] - het["by_portfolio"]["academic"]["large_20plus"]["irr"]],
    ]
    x2 = np.arange(len(c_lbls))
    ax.axhline(1.0, color="#7F8C8D", ls="--", lw=1.0)
    ax.axvline(2.5, color="#BDC3C7", ls=":", lw=1.1)
    ax.errorbar(x2 - 0.14, ind_vals, yerr=ind_err, fmt="o", color=C_IND, ecolor=C_IND, capsize=3.5, ms=6.2, lw=1.5, label="Industry")
    ax.errorbar(x2 + 0.14, aca_vals, yerr=aca_err, fmt="s", color=C_ACA, ecolor=C_ACA, capsize=3.5, ms=5.8, lw=1.5, label="Academic")
    ax.annotate("p = 3.8e-5", xy=(0 - 0.14, ind_vals[0]), xytext=(0.05, 3.05), fontsize=7.8, fontweight="bold", color=C_IND,
                arrowprops=dict(arrowstyle="->", color=C_IND, lw=0.9))
    ax.annotate("p = 0.036", xy=(3 - 0.14, ind_vals[3]), xytext=(2.62, 2.95), fontsize=7.8, fontweight="bold", color=C_IND,
                arrowprops=dict(arrowstyle="->", color=C_IND, lw=0.9))
    ax.set_xticks(x2)
    ax.set_xticklabels(c_lbls, fontsize=8.2)
    ax.set_ylim(0, 4.1)
    ax.set_ylabel("Poisson IRR vs. Exempt Controls (95% CI)")
    ax.set_title("B. Stale-Data Decay & Institutional Scale")
    ax.legend(loc="upper right", frameon=True, facecolor="white", fontsize=8.0)
    ax.grid(axis="y", alpha=0.25)

    # Panel C: Centralized vs Decentralized Hierarchy + April 2026 Phase-1 Within-Firm Spillover
    ax = axes[2]
    sp = v2["phase1_within_firm_spillover"]["hazards_apr2026"]
    # Plot April 2026 Phase-1 Spillover bars cleanly (Pre-Sweep vs April 2026 for Notified vs Unnotified Industry Firms)
    g_lbls = ["Unnotified Firms\n(0 Overdue ACTs)", "Sweep-Notified Firms\n(≥1 Overdue ACT)"]
    pre_p1 = [sp["exp0_apr0"]["hazard_pct"], sp["exp1_apr0"]["hazard_pct"]]
    apr_p1 = [sp["exp0_apr1"]["hazard_pct"], sp["exp1_apr1"]["hazard_pct"]]
    x3 = np.arange(2)
    ax.bar(x3 - 0.18, pre_p1, 0.34, color="#95A5A6", label="Pre-Sweep Baseline (Apr 25–Mar 26)")
    ax.bar(x3 + 0.18, apr_p1, 0.34, color=C_GOLD, label="April 2026 Post-Sweep Month")
    for i in range(2):
        ax.text(x3[i] - 0.18, pre_p1[i] + 0.025, f"{pre_p1[i]:.2f}%", ha="center", va="bottom", fontsize=8.0, color="#566573")
        ax.text(x3[i] + 0.18, apr_p1[i] + 0.025, f"{apr_p1[i]:.2f}%", ha="center", va="bottom", fontsize=8.2, fontweight="bold", color="#B7791F")
    ax.annotate("3.54× surge in\nexempt Phase 1\n(12/1,561 vs 38/17,492)", xy=(1 + 0.18, apr_p1[1]), xytext=(0.04, 0.62),
                fontsize=7.8, fontweight="bold", color="#B7791F",
                arrowprops=dict(arrowstyle="->", color="#B7791F", lw=1.1),
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=C_GOLD, alpha=0.94))
    nih_irr = wa_3["nih_fed_vs_academic"]["irr"]
    nih_p = wa_3["nih_fed_vs_academic"]["p"]
    ind_irr = wa_3["industry_vs_academic"]["irr"]
    ind_p = wa_3["industry_vs_academic"]["p"]
    ax.text(0.03, 0.34, f"Within-ACT Hierarchy vs Academic:\n• Industry: IRR = {ind_irr:.2f}× (p = {ind_p:.3f})\n• NIH/Fed:  IRR = {nih_irr:.2f}× (p = {nih_p:.3f})",
            transform=ax.transAxes, fontsize=7.4,
            bbox=dict(boxstyle="round,pad=0.32", facecolor="#F8F9F9", edgecolor="#BDC3C7", alpha=0.95))
    ax.set_xticks(x3)
    ax.set_xticklabels(g_lbls, fontsize=8.3)
    ax.set_ylabel("Exempt US Phase-1 Monthly Hazard (%)")
    ax.set_ylim(0, 1.05)
    ax.set_title("C. Exempt Phase-1 Spillover & Hierarchy")
    ax.legend(loc="upper left", frameon=True, facecolor="white", fontsize=7.6)
    ax.grid(axis="y", alpha=0.25)

    save_both(fig, "fig3_b2_governance_mechanisms.png")


# ==============================================================================
# FIGURE 4: H4 Lemons (Strategic Concealment) vs Organizational Inattention Test
# ==============================================================================
def make_fig4():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.6, 5.0), gridspec_kw={"width_ratios": [1.05, 1.05]})
    fig.subplots_adjust(wspace=0.28)
    lem = v2["lemons_h4"]
    gs = lem["group_summary"]
    reg = lem["regressions"]

    # Panel A: Raw rates across On-Time (n=2,822), Pre-Sweep Overdue (n=1,930), Post-Sweep Overdue (n=312)
    metrics = [
        ("Early\nTerminated", "terminated_share"),
        ("Primary\np ≥ 0.05", "nonsig_p05_share"),
        ("P-Value\nOmitted", "omitted_pval_share"),
        ("Mean SAE\nRate", "mean_sae_rate"),
        ("All-Cause\nMortality", "mean_mortality_rate"),
    ]
    x = np.arange(len(metrics))
    w = 0.25
    v_on = [gs["on_time"][k] * 100 for _, k in metrics]
    v_pre = [gs["pre_sweep_overdue"][k] * 100 for _, k in metrics]
    v_post = [gs["post_sweep_overdue"][k] * 100 for _, k in metrics]

    ax1.bar(x - w, v_on, w, color="#85929E", label=f"On-Time ≤365d (n = {gs['on_time']['n']:,})")
    ax1.bar(x, v_pre, w, color=C_GOLD, label=f"Pre-Sweep Overdue (n = {gs['pre_sweep_overdue']['n']:,})")
    ax1.bar(x + w, v_post, w, color=C_IND, label=f"Post-Sweep Overdue (n = {gs['post_sweep_overdue']['n']:,})")
    for i in range(len(metrics)):
        ax1.text(x[i] - w, v_on[i] + 1.0, f"{v_on[i]:.1f}", ha="center", va="bottom", fontsize=7.0, color="#566573")
        ax1.text(x[i], v_pre[i] + 1.0, f"{v_pre[i]:.1f}", ha="center", va="bottom", fontsize=7.0, color="#9A7D0A")
        ax1.text(x[i] + w, v_post[i] + 1.0, f"{v_post[i]:.1f}", ha="center", va="bottom", fontsize=7.2, fontweight="bold", color=C_IND)

    ax1.set_xticks(x)
    ax1.set_xticklabels([m[0] for m in metrics], fontsize=8.5)
    ax1.set_ylabel("Share / Rate (%)")
    ax1.set_ylim(0, 96)
    ax1.set_title("A. Clinical Outcomes Across Disclosure Timeliness (N = 5,064)")
    ax1.legend(loc="upper right", frameon=True, facecolor="white", fontsize=8.0)
    ax1.grid(axis="y", alpha=0.25)

    # Panel B: Cross-Sectional vs Within-Sponsor FE Coefficient Plot (Post-Sweep Overdue vs On-Time, percentage points!)
    y_keys = [
        ("Early Terminated", "is_terminated"),
        ("Primary p ≥ 0.05 (Null)", "nonsig_p05"),
        ("P-Value Omitted", "omitted_pval"),
        ("Serious Adverse Event Rate", "sae_rate"),
        ("All-Cause Mortality Rate", "mortality_rate"),
    ]
    yp = np.arange(len(y_keys))[::-1]
    ax2.axvline(0.0, color="#2C3E50", ls="--", lw=1.1)

    for i, (lbl, k) in enumerate(y_keys):
        y = yp[i]
        cs = reg[k]["all_post_overdue_vs_ontime"]
        fe = reg[k]["within_sponsor_fe"]["post_overdue_vs_ontime"]
        ax2.errorbar(cs["coef"] * 100, y + 0.14, xerr=1.96 * cs["se"] * 100, fmt="o", color=C_ACA, ecolor=C_ACA,
                     capsize=3.2, ms=5.8, lw=1.5, label="Cross-Sectional Controls" if i == 0 else None)
        ax2.errorbar(fe["coef"] * 100, y - 0.14, xerr=1.96 * fe["se"] * 100, fmt="s", color=C_IND, ecolor=C_IND,
                     capsize=3.2, ms=6.0, lw=1.6, label="Within-Sponsor Fixed Effects" if i == 0 else None)
        ax2.text(16.8, y - 0.14, f"FE: {fe['coef']*100:+.1f}pp (p={fe['p']:.2f})", va="center", ha="left", fontsize=7.9, color=C_IND, fontweight="bold")

    ax2.set_yticks(yp)
    ax2.set_yticklabels([k[0] for k in y_keys], fontsize=8.6)
    ax2.set_xlim(-17.5, 34.5)
    ax2.set_ylim(-0.75, 4.85)
    ax2.set_xlabel("Post-Sweep Overdue Minus On-Time Difference (Percentage Points, 95% CI)")
    ax2.set_title("B. Within-Sponsor FE Falsifies 'Lemons Hiding' Hypothesis")
    ax2.legend(loc="upper center", bbox_to_anchor=(0.50, -0.15), ncol=2, frameon=True, facecolor="white", fontsize=8.2)
    ax2.grid(axis="x", alpha=0.25)

    save_both(fig, "fig4_b2_lemons_vs_inattention.png")


if __name__ == "__main__":
    make_fig1()
    make_fig2()
    make_fig3()
    make_fig4()
    print("All 4 Study B-2 publication figures generated successfully!")
