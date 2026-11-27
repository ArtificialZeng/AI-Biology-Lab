#!/opt/anaconda3/bin/python3
"""
Generate ALL 4 Revised Publication-Grade Figures (300 DPI) for Study B-2
strictly synchronized to T00/T01/T02/T03/T04/T06/T07 verified JSONs
(including OFFICIAL CRAN R HonestDiD v0.2.8 FLCI & C-LF confidence sets!).
Mirrors PNGs into BOTH studies/study_b2_fda_2026_sweep/figures/ AND the conversation artifact directory.
"""
import json, os, shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_b2_fda_2026_sweep"
AUD = os.path.join(BASE, "audit")
FIG = os.path.join(BASE, "figures")
ART = "/Users/mac/.gemini/antigravity/brain/19f7514f-ce6a-4297-92cb-9803e648afd1"
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9.3,
    "axes.titlesize": 10.3,
    "axes.titleweight": "bold",
    "axes.labelsize": 9.3,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

b_def = json.load(open(os.path.join(AUD, "b2_definitive_t00_to_t04_results.json")))
b_cmp = json.load(open(os.path.join(AUD, "b2_compact_es_and_hazards.json")))
b_qhm = json.load(open(os.path.join(AUD, "b2_t00_quarterly_honestdid_and_mechanisms.json")))
b_t07 = json.load(open(os.path.join(AUD, "b2_t07_invariance_and_es_clean.json")))
b_hd  = json.load(open(os.path.join(AUD, "b2_t06_official_r_honestdid_results.json")))

C_NAVY  = "#1B365D"
C_CORAL = "#D95F02"
C_TEAL  = "#008080"
C_BLUE  = "#2166AC"
C_RED   = "#B2182B"
C_SLATE = "#636363"
C_GOLD  = "#E6AB02"

def save_fig(fig, name):
    p1 = os.path.join(FIG, name)
    p2 = os.path.join(ART, name)
    fig.savefig(p1, dpi=300, bbox_inches="tight", facecolor="white")
    shutil.copy2(p1, p2)
    plt.close(fig)
    print(f"Saved {name}", flush=True)

# =============================================================================
# FIGURE 1: Combined Posted+Unposted-Tracking Monthly Hazards & T07 Event Studies
# =============================================================================
fig, axes = plt.subplots(1, 2, figsize=(14.4, 4.9), gridspec_kw={"width_ratios": [1.08, 1.05]})

# Panel A: Monthly Hazards (4 cells: Overdue vs On-Time x Industry vs Academic)
ax = axes[0]
hz = b_cmp["uncensored_monthly_hazards"]
yms = [r["ym"] for r in hz]
x = np.arange(len(yms))
ax.axvspan(11.5, len(yms) - 0.5, color="#FFF3E0", alpha=0.55, zorder=0)
ax.axvline(11.5, color=C_RED, linestyle="--", lw=1.3, zorder=2)
ax.plot(x, [r["od1_i1_h"] for r in hz], "o-", color=C_CORAL, lw=2.2, ms=5.5, label="Overdue Industry ACTs (Treated Target)")
ax.plot(x, [r["od0_i1_h"] for r in hz], "s--", color=C_NAVY, lw=1.6, ms=4.5, alpha=0.85, label="On-Time Industry ACTs (Within-Sector Control)")
ax.plot(x, [r["od1_i0_h"] for r in hz], "^-", color=C_TEAL, lw=1.8, ms=5.0, label="Overdue Academic ACTs")
ax.plot(x, [r["od0_i0_h"] for r in hz], "d:", color="#7570B3", lw=1.4, ms=4.2, alpha=0.8, label="On-Time Academic ACTs")
ax.text(11.75, 9.45, "Mar 30, 2026\nFDA Sweep", color=C_RED, fontsize=8.0, fontweight="bold", va="top")
ax.set_xticks(x[::2])
ax.set_xticklabels([yms[i][2:] for i in range(0, len(yms), 2)], rotation=30, ha="right", fontsize=8.0)
ax.set_ylim(0.0, 10.2)
ax.set_ylabel("Monthly Results-Submission Hazard (%)\n(Combined Posted + Unposted Tracking Proxy)")
ax.set_title("A. Strict-ACTUAL ACT Monthly Submission Hazards\n(Risk-Set Normalized Across All 4 Treatment Cells)")
ax.legend(frameon=True, facecolor="white", edgecolor="#DDDDDD", fontsize=7.4, loc="lower left", ncol=2)

# Panel B: Monthly Event Study (ref=2026-02; T07 Fixed Pre-Window df=10)
ax = axes[1]
es_ind = b_t07["monthly_es_fixed_pre"]["ot_od_ind_only_ref2602"]["series"]
es_ddd = b_t07["monthly_es_fixed_pre"]["ot_od_ddd_ref2602"]["series"]
ax.axvspan(11.5, len(yms) - 0.5, color="#FFF3E0", alpha=0.55, zorder=0)
ax.axvline(11.5, color=C_RED, linestyle="--", lw=1.3, zorder=2)
ax.axhline(0.0, color="#333333", linestyle="-", lw=0.9, zorder=1)
# Plot log-IRR (exact additive scale so error bars are symmetric!)
yi = [r["log_irr"] for r in es_ind]
ei = [1.96 * r["se"] for r in es_ind]
yd = [r["log_irr"] for r in es_ddd]
ed = [1.96 * r["se"] for r in es_ddd]
ax.errorbar(x - 0.14, yi, yerr=ei, fmt="o-", color=C_CORAL, capsize=2.8, lw=1.7, ms=5.0,
            label="Industry Overdue vs. On-Time DiD (Pre Wald df=10, p=0.429)")
ax.errorbar(x + 0.14, yd, yerr=ed, fmt="s--", color=C_NAVY, capsize=2.8, lw=1.4, ms=4.5, alpha=0.85,
            label="Overdue vs. On-Time DDD (Pre Wald df=10, p=0.059)")
ax.text(10.0, -0.78, "Ref\n26-02", ha="center", fontsize=7.3, color="#444444")
ax.text(11.0, -0.78, "Trans.\n26-03", ha="center", fontsize=7.3, color=C_RED)
ax.set_xticks(x[::2])
ax.set_xticklabels([yms[i][2:] for i in range(0, len(yms), 2)], rotation=30, ha="right", fontsize=8.0)
ax.set_ylim(-0.98, 2.85)
ax.set_ylabel("Event-Study Coefficient: log(IRR) vs. Feb 2026\n(Sponsor-Clustered 95% CI)")
ax.set_title("B. Monthly Event Studies (T07 Fixed Pre-Window df=10)\n(Sharp April 2026 Industry Jump; Within-ACT Pre-Trend Rejected)")
ax.legend(frameon=True, facecolor="white", edgecolor="#DDDDDD", fontsize=7.5, loc="upper left")

fig.tight_layout(w_pad=2.2)
save_fig(fig, "fig1_b2_event_study_and_hazards.png")

# =============================================================================
# FIGURE 2: Identification Forest Plot (T00+T01+Placebos) & Official CRAN HonestDiD v0.2.8
# =============================================================================
fig, axes = plt.subplots(1, 2, figsize=(14.8, 5.4), gridspec_kw={"width_ratios": [1.18, 1.0]})

# Panel A: Forest plot across T00 Combined Proxy, Cert-Excluded, T01 Exact Zero-Censoring, and Placebos
ax = axes[0]
ot = b_def["t00_uncensored_ot_od_ddd"]
t1 = b_def["t01_strict_date_verified_uniform_horizon"]
dp = b_t07["t00_excl_disp_cert_sensitivity"]

rows_f = [
    ("Combined Proxy (Apr–Jun, Ind OD/OT)",      ot["true_ot_od_apr_jun_2026"]["industry_od_vs_ot"]["irr"], ot["true_ot_od_apr_jun_2026"]["industry_od_vs_ot"]["ci95"], C_CORAL),
    ("Combined Proxy (Apr–Sep, Ind OD/OT)",      ot["true_ot_od_apr_sep_2026"]["industry_od_vs_ot"]["irr"], ot["true_ot_od_apr_sep_2026"]["industry_od_vs_ot"]["ci95"], C_CORAL),
    ("Combined Proxy (Apr–Oct, Ind OD/OT)",      ot["true_ot_od_apr_oct_2026"]["industry_od_vs_ot"]["irr"], ot["true_ot_od_apr_oct_2026"]["industry_od_vs_ot"]["ci95"], C_CORAL),
    ("Combined Proxy (Apr–Sep, OT/OD DDD)",      ot["true_ot_od_apr_sep_2026"]["pois_ddd"]["irr"],          ot["true_ot_od_apr_sep_2026"]["pois_ddd"]["ci95_spon"],     C_NAVY),
    ("Excl. Delayed Certs (Apr–Sep, OT/OD DDD)", dp["ot_od_ddd_irr"],                                       dp["ot_od_ddd_ci95"],                                        C_NAVY),
    ("T01 Zero-Cens H90 (thru Jul 31, DDD)",     t1["H90_fullmonth_apr_jul"]["ot_od_ddd"]["pois_ddd"]["irr"], t1["H90_fullmonth_apr_jul"]["ot_od_ddd"]["pois_ddd"]["ci95_spon"], C_TEAL),
    ("T01 Zero-Cens H120 (thru Jun 30, DDD)",    t1["H120_fullmonth_apr_jun"]["ot_od_ddd"]["pois_ddd"]["irr"], t1["H120_fullmonth_apr_jun"]["ot_od_ddd"]["pois_ddd"]["ci95_spon"], C_TEAL),
    ("T01 Zero-Cens H150 (thru Jun 28, DDD)",    t1["H150_exactday_thru_jun28"]["ot_od_ddd"]["pois_ddd"]["irr"], t1["H150_exactday_thru_jun28"]["ot_od_ddd"]["pois_ddd"]["ci95_spon"], C_TEAL),
    ("Clean-Exempt DDD (Apr–Sep, Proxy)",        b_def["t00_uncensored_exempt_ddd"]["true_ddd_clean_exempt_sep"]["pois_ddd_irr"], b_def["t00_uncensored_exempt_ddd"]["true_ddd_clean_exempt_sep"]["pois_ddd_ci95"], "#7570B3"),
    ("Placebo 2025 (Apr–Sep, OT/OD DDD)",        ot["true_ot_od_placebo_2025"]["pois_ddd"]["irr"],          ot["true_ot_od_placebo_2025"]["pois_ddd"]["ci95_spon"],      C_SLATE),
    ("Placebo 2024 (Apr–Sep, OT/OD DDD)",        ot["true_ot_od_placebo_2024"]["pois_ddd"]["irr"],          ot["true_ot_od_placebo_2024"]["pois_ddd"]["ci95_spon"],      C_SLATE),
]
yy = np.arange(len(rows_f))[::-1]
ax.axvline(1.0, color="#333333", linestyle="--", lw=1.1)
ax.axhline(1.5, color="#BBBBBB", linestyle=":", lw=1.0)
for y_i, (lbl, pt, ci, col) in zip(yy, rows_f):
    ax.plot([ci[0], ci[1]], [y_i, y_i], color=col, lw=2.1)
    ax.scatter([pt], [y_i], color=col, s=55, zorder=4, edgecolors="white", linewidths=0.5)
    ax.text(ci[1] * 1.05, y_i, f"{pt:.2f} [{ci[0]:.2f}, {ci[1]:.2f}]", va="center", fontsize=7.6, color="#222222")
ax.set_yticks(yy)
ax.set_yticklabels([r[0] for r in rows_f], fontsize=8.0)
ax.set_xscale("log")
ax.set_xlim(0.48, 9.8)
ax.set_xlabel("Poisson Incidence Rate Ratio (IRR, Log Scale; Sponsor-Clustered 95% CI)")
ax.set_title("A. Identification & Date-Verified Zero-Censoring Robustness\n(RI 1,000-Perm p = 0.010; Placebos 2024–25 Indistinguishable from 1.0)")

# Panel B: Official CRAN R HonestDiD v0.2.8 Sensitivity (Quarterly 2026Q2 Impact — Delta^RM C-LF & Delta^SD FLCI!)
ax = axes[1]
sp_ddd = b_hd["specs"]["quarterly_ot_od_ddd_ref2026Q1_Q2_impact"]
sp_ind = b_hd["specs"]["quarterly_ot_od_ind_only_ref2026Q1_Q2_impact"]

# Show BOTH Delta^RM (Mbar = 0, 0.25, 0.5, 0.75, 1.0) for Quarterly Q2 DDD & Ind-Only!
mbars = [0.0, 0.25, 0.5, 0.75, 1.0]
def get_rm(sp):
    d = {float(r["Mbar"]): (float(r["lb"]), float(r["ub"])) for r in sp["delta_RM_CLF"]}
    return [d[m] for m in mbars]
rm_ddd = get_rm(sp_ddd)
rm_ind = get_rm(sp_ind)
xm = np.arange(len(mbars))
ax.axhline(0.0, color=C_RED, linestyle="--", lw=1.2, label="Null Effect (log IRR = 0, IRR = 1)")
for i, m in enumerate(mbars):
    # DDD C-LF bar
    lb_d, ub_d = rm_ddd[i]
    ax.plot([i - 0.14, i - 0.14], [lb_d, ub_d], color=C_NAVY, lw=3.0)
    ax.scatter([i - 0.14], [sp_ddd["theta_hat_log_irr"]], color=C_NAVY, s=42, zorder=4)
    # Ind-Only C-LF bar
    lb_i, ub_i = rm_ind[i]
    ax.plot([i + 0.14, i + 0.14], [lb_i, ub_i], color=C_CORAL, lw=3.0)
    ax.scatter([i + 0.14], [sp_ind["theta_hat_log_irr"]], color=C_CORAL, s=42, zorder=4)
    if i <= 2:
        ax.text(i - 0.14, lb_d - 0.12, f"{lb_d:+.2f}", ha="center", va="top", fontsize=7.0, color=C_NAVY)
        ax.text(i + 0.14, ub_i + 0.08, f"{ub_i:+.2f}", ha="center", va="bottom", fontsize=7.0, color=C_CORAL)

ax.plot([], [], color=C_NAVY, lw=3.0, label="Q2 OT/OD DDD (C-LF; >0 at Mbar≤0.25, M≤0.10 FLCI)")
ax.plot([], [], color=C_CORAL, lw=3.0, label="Q2 Ind OD/OT (C-LF; >0 at Mbar≤0.25, M≤0.05 FLCI)")
ax.text(0.50, 0.04, "Note: Official HonestDiD v0.2.8 (gridPoints=1000, approx. grid endpoints).\nCIs cross 0 at Mbar=0.50 (DDD lb=-0.016, Ind lb=-0.028)\nand Q2–Q3 Ind avg includes 0 even at Mbar=0.",
        transform=ax.transAxes, ha="center", va="bottom", fontsize=7.2,
        bbox=dict(boxstyle="round,pad=0.35", fc="#F7F7F7", ec="#CCCCCC"))
ax.set_xticks(xm)
ax.set_xticklabels([f"M̄ = {m:g}" for m in mbars], fontsize=8.5)
ax.set_ylim(-1.15, 2.35)
ax.set_ylabel("Official HonestDiD v0.2.8 95% C-LF CS\n(log IRR Scale; 2026Q2 vs. 2026Q1)")
ax.set_title("B. Official CRAN HonestDiD v0.2.8 Sensitivity (Δᴿᴹ)\n(Transparent Partial-Identification Breakdown Bounds)")
ax.legend(frameon=True, facecolor="white", edgecolor="#DDDDDD", fontsize=7.3, loc="upper left")

fig.tight_layout(w_pad=2.2)
save_fig(fig, "fig2_b2_identification_and_placebos.png")

# =============================================================================
# FIGURE 3: Mechanisms & Heterogeneity (Combined Proxy, Apr–Sep 2026) + Null Phase-1 Spillover
# =============================================================================
fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.6))
het = b_qhm["t00_uncensored_heterogeneity_apr_sep"]
sp1 = b_qhm["t00_uncensored_phase1_spillover"]

def plot_het_panel(ax, items, title, xlabel):
    xx = np.arange(len(items))
    pts = [v["irr_vs_exempt"] for _, v, _ in items]
    lo  = [v["irr_vs_exempt"] - v["ci95"][0] for _, v, _ in items]
    hi  = [v["ci95"][1] - v["irr_vs_exempt"] for _, v, _ in items]
    cols = [c for _, _, c in items]
    ax.axhline(1.0, color="#333333", linestyle="--", lw=1.0)
    for i in range(len(items)):
        ax.errorbar(xx[i], pts[i], yerr=[[lo[i]], [hi[i]]], fmt="o", color=cols[i], capsize=4.5, lw=2.0, ms=7.5)
        ax.text(xx[i], pts[i] + hi[i] + 0.14, f"{pts[i]:.2f}\n(p={items[i][1]['p']:.3f})" if items[i][1]['p'] >= 0.001 else f"{pts[i]:.2f}\n(p<0.001)",
                ha="center", va="bottom", fontsize=7.7, fontweight="bold", color=cols[i])
    ax.set_xticks(xx)
    ax.set_xticklabels([lbl for lbl, _, _ in items], fontsize=8.3)
    ax.set_xlim(-0.55, len(items) - 0.45)
    ax.set_ylim(0.2, 5.45)
    ax.set_xlabel(xlabel, fontsize=8.6)
    ax.set_title(title)

plot_het_panel(axes[0], [
    ("Small\n(1–2 FDA)",  het["by_portfolio"]["industry"]["small_1_2"],   C_CORAL),
    ("Medium\n(3–19)",    het["by_portfolio"]["industry"]["mid_3_19"],    C_GOLD),
    ("Large\n(20+ FDA)",  het["by_portfolio"]["industry"]["large_20plus"], C_BLUE),
], "A. By Sponsor Portfolio Size\n(Large Pharma Already Complied: 8.2%/mo)", "Sponsor FDA-Regulated Trial Portfolio")
axes[0].set_ylabel("Industry Overdue IRR vs. Clean Exempt Controls\n(Combined Proxy, Apr–Sep 2026; 95% CI)")

plot_het_panel(axes[1], [
    ("Recent\n(2023–25)", het["by_cohort"]["industry"]["recent_2023_2025"], C_CORAL),
    ("Mid\n(2021–22)",    het["by_cohort"]["industry"]["mid_2021_2022"],    C_GOLD),
    ("Legacy\n(2017–20)", het["by_cohort"]["industry"]["deep_2017_2020"],   C_SLATE),
], "B. By Delinquency Cohort Depth\n(Response Concentrated in Recent Overdue)", "Primary Completion Year Cohort")

# Panel C: Backlog size + explicit NULL Within-Firm Phase-1 Spillover!
sp_item = {"irr_vs_exempt": sp1["pois_irr"], "ci95": sp1["ci95"], "p": sp1["p"]}
plot_het_panel(axes[2], [
    ("Backlog\n1 Trial",   het["by_backlog"]["industry"]["backlog_1"],     C_NAVY),
    ("Backlog\n2–4",       het["by_backlog"]["industry"]["backlog_2_4"],   C_TEAL),
    ("Backlog\n5+",        het["by_backlog"]["industry"]["backlog_5plus"], C_BLUE),
    ("Exempt\nPhase-1\nSpillover", sp_item, C_RED),
], "C. By Overdue Backlog & Phase-1 Spillover\n(No Detectable Phase-1 Spillover: IRR=0.89)", "Mar 29 Backlog / Within-Firm Spillover")

fig.tight_layout(w_pad=1.8)
save_fig(fig, "fig3_b2_governance_mechanisms.png")

# =============================================================================
# FIGURE 4: Honest H4 Lemons Audit (Within-Sponsor FE CIs + ±5pp TOST + Selection Check)
# =============================================================================
fig, axes = plt.subplots(1, 2, figsize=(13.8, 4.8), gridspec_kw={"width_ratios": [1.35, 0.95]})

# Panel A: Within-Sponsor FE point estimates + exact 95% CIs + ±5 pp TOST equivalence band
ax = axes[0]
reg = b_def["t03_t04_h4_lemons_audit"]["regressions"]
h4_rows = [
    ("Mortality Rate\n(n=4,180; MDE80=4.0pp)",     reg["mortality_rate"]["within_sponsor_post_vs_ontime"], C_TEAL),
    ("Serious Adv. Event Rate\n(n=4,179; MDE80=4.4pp)", reg["sae_rate"]["within_sponsor_post_vs_ontime"], C_TEAL),
    ("Early Termination\n(n=4,224; MDE80=10.1pp)", reg["is_terminated"]["within_sponsor_post_vs_ontime"], C_NAVY),
    ("Omitted Primary p-Value\n(n=4,224; MDE80=10.0pp)", reg["omitted_pval"]["within_sponsor_post_vs_ontime"], C_NAVY),
    ("Primary p ≥ 0.05 (Null)\n(n=1,035; MDE80=22.9pp!)", reg["nonsig_p05"]["within_sponsor_post_vs_ontime"], C_RED),
]
yh = np.arange(len(h4_rows))[::-1]
ax.axvspan(-5, 5, color="#E0F2F1", alpha=0.85, zorder=0, label="±5 pp Post-Hoc Exploratory TOST Margin (Safety p≤0.002)")
ax.axvline(0.0, color="#333333", linestyle="--", lw=1.0, zorder=1)
for y_i, (lbl, d, col) in zip(yh, h4_rows):
    pt = d["coef_pp"]
    lo, hi = d["ci95_pp"]
    ax.plot([lo, hi], [y_i, y_i], color=col, lw=2.6, zorder=3)
    ax.scatter([pt], [y_i], color=col, s=62, zorder=4, edgecolors="white", linewidths=0.6)
    t_str = f"{pt:+.2f} pp [{lo:+.1f}, {hi:+.1f}]"
    ax.text(hi + 0.7, y_i, t_str, va="center", fontsize=7.9, color=col if col == C_RED else "#222222", fontweight="bold" if col == C_RED else "normal")
ax.set_yticks(yh)
ax.set_yticklabels([r[0] for r in h4_rows], fontsize=8.3)
ax.set_xlim(-18.5, 30.5)
ax.set_ylim(-0.6, 4.95)
ax.set_xlabel("Within-Sponsor FE Difference: Post-Sweep Overdue minus On-Time (pp, 95% CI)")
ax.set_title("A. H4 Outcome Audit: Post-Hoc ±5 pp Safety Sensitivity\nvs. Underpowered Efficacy Test (Does Not Exclude Concealment)")
ax.legend(frameon=True, facecolor="white", edgecolor="#DDDDDD", fontsize=7.8, loc="upper right")

# Panel B: Selection Audit — Unposted-in-QC vs Posted Termination Rate & QC Clearance Rate
ax = axes[1]
sel = b_def["t03_t04_h4_lemons_audit"]["unposted_vs_posted_termination_selection"]
qcg = b_def["t00_post_sweep_qc_clearance_gap"]
cats = ["Post-Sweep Overdue\nTerminated Rate\n(Posted vs. In-QC)", "Apr–Aug Overdue\nQC Clearance Rate\n(Industry vs. Academic)"]
v_a = [sel["posted_terminated_pct"], qcg["industry_qc_cleared_pct"]]
v_b = [sel["unposted_in_qc_terminated_pct"], qcg["academic_qc_cleared_pct"]]
lbl_a = [f"QC-Cleared\n(n={sel['posted_post_overdue_n']})\n{v_a[0]:.1f}%", f"Industry\n(n={qcg['industry_n']})\n{v_a[1]:.1f}%"]
lbl_b = [f"In-QC Unposted\n(n={sel['unposted_in_qc_post_overdue_n']})\n{v_b[0]:.1f}%", f"Academic\n(n={qcg['academic_n']})\n{v_b[1]:.1f}%"]
ps  = [sel["fisher_p"], qcg["fisher_p"]]
xc = np.arange(2)
w = 0.34
ax.bar(xc - w/2, v_a, width=w, color=C_NAVY, edgecolor="white")
ax.bar(xc + w/2, v_b, width=w, color=C_CORAL, edgecolor="white")
for i in range(2):
    ax.text(xc[i] - w/2, v_a[i] + 1.8, lbl_a[i], ha="center", va="bottom", fontsize=7.5, color=C_NAVY, fontweight="bold")
    ax.text(xc[i] + w/2, v_b[i] + 1.8, lbl_b[i], ha="center", va="bottom", fontsize=7.5, color=C_CORAL, fontweight="bold")
    ax.text(xc[i], max(v_a[i], v_b[i]) + 17.0, f"Fisher p = {ps[i]:.3f}", ha="center", fontsize=8.2, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.25", fc="#F5F5F5", ec="#CCCCCC"))
ax.set_xticks(xc)
ax.set_xticklabels(cats, fontsize=8.3)
ax.set_ylim(0, 108)
ax.set_ylabel("Share of Trials (%)")
ax.set_title("B. QC-Clearance Selection & Backlog Audits\n(No Evidence of Selective Stalling of Terminated Trials)")

fig.tight_layout(w_pad=2.2)
save_fig(fig, "fig4_b2_lemons_vs_inattention.png")
print("ALL 4 REVISED B-2 FIGURES GENERATED!", flush=True)
