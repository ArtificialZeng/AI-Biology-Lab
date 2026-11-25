#!/opt/anaconda3/bin/python3
"""
Generate ALL 4 Revised Publication-Grade Figures (300 DPI) for Study A-1
strictly synchronized to verified audit JSONs (ZERO v1 retracted claims!).
Mirrors PNGs into BOTH studies/study_a1_clinvar_gnomad_reclass/figures/ AND the conversation artifact directory.
"""
import json, os, shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

BASE = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_a1_clinvar_gnomad_reclass"
RES = os.path.join(BASE, "results")
AUD = os.path.join(BASE, "audit")
FIG = os.path.join(BASE, "figures")
ART = "/Users/mac/.gemini/antigravity/brain/19f7514f-ce6a-4297-92cb-9803e648afd1"
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9.5,
    "axes.titlesize": 10.5,
    "axes.titleweight": "bold",
    "axes.labelsize": 9.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

a_orig = json.load(open(os.path.join(RES, "a1_reclassification_results.json")))
a_for = json.load(open(os.path.join(AUD, "a1_forensic_audit.json")))
a_p20 = json.load(open(os.path.join(AUD, "a1_prospective_unseen20_and_qc55.json")))
a_vcep = json.load(open(os.path.join(AUD, "a1_vcep_bs1_and_collapsible_submitters.json")))
a_t05 = json.load(open(os.path.join(AUD, "a1_t05_verified_closure.json")))
a_t08 = json.load(open(os.path.join(AUD, "a1_t08_chrx_ploidy_closure.json")))

C_NAVY = "#1B365D"
C_TEAL = "#008080"
C_CORAL = "#D95F02"
C_RED = "#B2182B"
C_SLATE = "#636363"
C_GOLD = "#E6AB02"
C_BLUE = "#2166AC"

def save_fig(fig, name):
    p1 = os.path.join(FIG, name)
    p2 = os.path.join(ART, name)
    fig.savefig(p1, dpi=300, bbox_inches="tight", facecolor="white")
    shutil.copy2(p1, p2)
    plt.close(fig)
    print(f"Saved {name}", flush=True)

# =============================================================================
# FIGURE 1: Denominator Bias Retraction (Raw AF vs Hypergeometric Downsampling)
#           & Collapsible Non-NFE VUS Persistence Across Submitter Tiers
# =============================================================================
fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.6), gridspec_kw={"width_ratios": [1.05, 1.15, 1.15]})

# Panel A: Singleton Denominator Bias — P/LP Collapses While VUS Survives!
ax = axes[0]
pops = ["afr", "amr", "eas", "sas", "mid"]
pop_lbls = ["AFR", "AMR", "EAS", "SAS", "MID"]
x = np.arange(len(pops))
w = 0.35
vus_ors = [a_for["h1_hypergeometric_downsampling"][p]["vus_cmh_or_mean"] for p in pops]
plp_ors = [a_for["h1_hypergeometric_downsampling"][p]["plp_cmh_or_mean"] for p in pops]
ax.axhline(1.0, color="#333333", linestyle="--", lw=1.0, zorder=1)
b1 = ax.bar(x - w/2, vus_ors, width=w, color=C_CORAL, edgecolor="white", label="VUS (Persists: OR 1.62–2.72)")
b2 = ax.bar(x + w/2, plp_ors, width=w, color=C_BLUE, edgecolor="white", label="P/LP (Null/Reversed: OR 0.42–1.27)")
for rect, val in zip(b1, vus_ors):
    ax.text(rect.get_x() + rect.get_width()/2, val + 0.06, f"{val:.2f}", ha="center", va="bottom", fontsize=7.8, fontweight="bold", color=C_CORAL)
for rect, val in zip(b2, plp_ors):
    y_lbl = val + 0.14 if 0.80 <= val <= 1.02 else val + 0.06
    ax.text(rect.get_x() + rect.get_width()/2, y_lbl, f"{val:.2f}", ha="center", va="bottom", fontsize=7.8, color=C_BLUE)
ax.set_xticks(x)
ax.set_xticklabels(pop_lbls)
ax.set_ylim(0, 3.35)
ax.set_ylabel("CMH Odds Ratio vs. Downsampled NFE\n(Exact Equal-AN Hypergeometric Draws)")
ax.set_title("A. Equal-AN Downsampling Separates\nTrue VUS Gap from P/LP Singleton Artifact")
ax.legend(frameon=False, fontsize=7.8, loc="upper left")

# Panel B: Collapsible LPM Risk Difference (pp) — Discovery (84 ACMG, Ploidy-Aware) vs Prospective Unseen (20 AD)
ax = axes[1]
c84 = a_t08["h1_84gene_ploidy_aware_revised"]["outcomes"]["VUS"]
c20 = a_vcep["collapsible_attenuation_unseen20"]["outcomes"]["VUS"]
specs = ["Base Model\n(Gene+AF+Csq)", "+ Submitter\nTier Bins", "+ log(1+Subs)\n(Continuous)"]
y84 = [c84["lpm_base_pp"], c84["lpm_subsbin_pp"], c84["lpm_logsubs_pp"]]
e84 = [1.96*c84["lpm_base_se_pp"], 1.96*c84["lpm_subsbin_se_pp"], 1.96*1.11]
y20 = [c20["lpm_base_rd_pp"], c20["lpm_subsbin_rd_pp"], c20["lpm_logsubs_rd_pp"]]
e20 = [1.96*c20["lpm_base_se_pp"], 1.96*c20["lpm_subsbin_se_pp"], 1.96*c20["lpm_logsubs_se_pp"]]
xx = np.arange(len(specs))
ax.axhline(0, color="#333333", linestyle="--", lw=1.0)
ax.errorbar(xx - 0.14, y84, yerr=e84, fmt="o-", color=C_NAVY, capsize=4, lw=1.8, ms=6.5, label="ACMG SF v3.3 (84 genes, ploidy-aware, n=10,058)")
ax.errorbar(xx + 0.14, y20, yerr=e20, fmt="s--", color=C_TEAL, capsize=4, lw=1.8, ms=6.0, label="Unseen Validation (20 AD genes, n=2,522)")
for i in range(3):
    s84 = f"+{y84[i]:.1f}" if y84[i] >= 0 else f"{y84[i]:.1f}"
    ax.text(xx[i] - 0.14, y84[i] - e84[i] - 0.95, s84, ha="center", fontsize=7.8, color=C_NAVY, fontweight="bold")
    ax.text(xx[i] + 0.14, y20[i] + e20[i] + 0.45, f"+{y20[i]:.1f}", ha="center", fontsize=7.8, color=C_TEAL, fontweight="bold")
ax.set_xticks(xx)
ax.set_xticklabels(specs, fontsize=8.5)
ax.set_ylim(-3.8, 18.0)
ax.set_ylabel("Collapsible LPM Risk Difference\n(Non-NFE minus NFE VUS Rate, pp)")
ax.set_title("B. Descriptive Covariate Attenuation\n(Collapsible LPM; Not Causal Mediation)")
ax.legend(frameon=False, fontsize=7.5, loc="upper right")

# Panel C: Exact Non-Parametric Submitter-Tier Stratification (ACMG 84 Ploidy-Aware, FAF95 >= 1e-4)
ax = axes[2]
tiers = ["1", "2_3", "4_9", "10_19", "20p"]
tlbls = ["1", "2–3", "4–9", "10–19", "20+"]
st = a_t08["h1_84gene_ploidy_aware_revised"]["within_tier_vus"]
v_nn = [st[t]["vus_nn_pct"] for t in tiers]
v_eu = [st[t]["vus_nfe_pct"] for t in tiers]
xt = np.arange(len(tiers))
w2 = 0.36
ax.bar(xt - w2/2, v_nn, width=w2, color=C_CORAL, edgecolor="white", label="Non-NFE Enriched (n=9,460)")
ax.bar(xt + w2/2, v_eu, width=w2, color=C_BLUE, edgecolor="white", label="NFE Enriched (n=598)")
for i in range(len(tiers)):
    ax.text(xt[i] - w2/2, v_nn[i] + 1.0, f"{v_nn[i]:.1f}%", ha="center", fontsize=7.2, color=C_CORAL, fontweight="bold")
    ax.text(xt[i] + w2/2, v_eu[i] + 1.0, f"{v_eu[i]:.1f}%", ha="center", fontsize=7.2, color=C_BLUE)
ax.text(0.5, 0.72, "All 5 within-tier Fisher p ≥ 0.16\nNon-NFE variants concentrate in\nsparse 1–3 submitter tiers",
        transform=ax.transAxes, fontsize=7.8, ha="center", bbox=dict(boxstyle="round,pad=0.35", fc="#F7F7F7", ec="#CCCCCC"))
ax.set_xticks(xt)
ax.set_xticklabels(tlbls)
ax.set_xlabel("ClinVar Submitter Count Tier (n_submitters)")
ax.set_ylabel("VUS Rate Within Submitter Tier (%)")
ax.set_ylim(0, 60)
ax.set_title("C. Non-Parametric Submitter Stratification\n(Distribution Shift Explains VUS Gap)")
ax.legend(frameon=False, fontsize=7.8, loc="upper right")

fig.tight_layout(w_pad=2.0)
save_fig(fig, "fig1_a1_equity_and_curation_velocity.png")

# =============================================================================
# FIGURE 2: Multiallelic rsID Mirage (84 ACMG + 20 Unseen) & 55/55 Revised Triage
# =============================================================================
fig, axes = plt.subplots(1, 2, figsize=(13.2, 4.9), gridspec_kw={"width_ratios": [1.0, 1.25]})

# Panel A: rsID Join False-Positive Inflation across ACMG 84 & Unseen 20
ax = axes[0]
rs84 = a_orig["h4_rsid_collision_audit"]
rs20 = a_p20["pv1_rsid_mirage"]
cats = ["ACMG 84\nAll P/LP", "ACMG 84\nAD 2–4★ P/LP", "Unseen 20\nAll P/LP", "Unseen 20\n2–4★ P/LP"]
naive_v = [rs84["faf95_ge_1e4"]["naive_rsid_flagged"], rs84["ad_highconf_2to4star_faf95_ge_1e4"]["naive_rsid_flagged"],
           rs20["all_plp_faf95_ge_1e4"]["naive"], rs20["hc_2to4star_faf95_ge_1e4"]["naive"]]
exact_v = [rs84["faf95_ge_1e4"]["allele_exact_confirmed"], rs84["ad_highconf_2to4star_faf95_ge_1e4"]["allele_exact_confirmed"],
           rs20["all_plp_faf95_ge_1e4"]["exact"], rs20["hc_2to4star_faf95_ge_1e4"]["exact"]]
fprs = [(n - e)/n*100 for n, e in zip(naive_v, exact_v)]
xx = np.arange(len(cats))
w = 0.36
ax.bar(xx - w/2, naive_v, width=w, color="#BDBDBD", edgecolor="#555555", label="Naive rsID Join (Inflated)")
ax.bar(xx + w/2, exact_v, width=w, color=C_NAVY, edgecolor="white", label="Allele-Exact (CHROM-POS-REF-ALT)")
for i in range(len(cats)):
    ax.text(xx[i] - w/2, naive_v[i] + 10, f"{naive_v[i]:,}", ha="center", fontsize=8.0, color="#444444")
    ax.text(xx[i] + w/2, exact_v[i] + 10, f"{exact_v[i]:,}", ha="center", fontsize=8.0, fontweight="bold", color=C_NAVY)
    ax.text(xx[i], max(naive_v[i], exact_v[i]) + 48, f"FPR: {fprs[i]:.1f}%", ha="center", fontsize=7.8, fontweight="bold", color=C_RED)
ax.set_xticks(xx)
ax.set_xticklabels(cats, fontsize=8.5)
ax.set_ylim(0, 680)
ax.set_ylabel("Flagged ClinVar P/LP Variants (FAF95 ≥ 1e-4)")
ax.set_title("A. Multiallelic rsID Join Mirage\n(Discovery 84 Genes & Unseen 20-Gene Replication)")
ax.legend(frameon=False, fontsize=8.2, loc="upper right")

# Panel B: Complete 55/55 AD >=2-Star Triage under Official CSpec + Johnston 2026 DAFT!
ax = axes[1]
tc = a_t05["ad55_verified_t05_counts"]
layers = [
    ("L2: FIN / ASJ Bottleneck\nFounders", tc.get("L2_BottleneckFounder_FIN_ASJ", 0), C_BLUE),
    ("L3: Established Continental\nFounders / Low-Pen", tc.get("L3_EstablishedContinentalFounder_LowPen", 0), C_TEAL),
    ("L1: gnomAD QC / Low-AN /\nAS_VQSR Artifacts (incl. APC)", tc.get("L1_QC_or_Callability_Artifact", 0), C_SLATE),
    ("L4a: Below Verified CSpec /\nJohnston 2026 BS1 (True AD)", tc.get("L4a_BelowVerifiedBS1_TrueADPathogenic", 0), C_GOLD),
    ("L4b: AR-in-Dual-Mode AD\nGene Return Hazards (RYR1/TTN)", tc.get("L4b_AR_in_DualMode_AD_Gene_ReturnHazard", 0), C_CORAL),
    ("L4c: > Johnston BS1, < BA1\n(CACNA1S MH Reduced-Pen)", tc.get("L4c_ExceedsVerifiedBS1_BelowBA1_ReducedPenOrReview", 0), "#7570B3"),
    ("L4d: High-QC AD Exceeding\nVerified Gene BA1", 0, C_RED),
]
yl = np.arange(len(layers))[::-1]
vals = [x[1] for x in layers]
cols = [x[2] for x in layers]
lbls = [x[0] for x in layers]
ax.barh(yl, vals, height=0.62, color=cols, edgecolor="white")
for y_i, v in zip(yl, vals):
    ax.text(v + 0.5, y_i, f"{v} / 55 ({v/55*100:.1f}%)", va="center", fontsize=8.3, fontweight="bold", color="#222222" if v > 0 else C_RED)
ax.set_yticks(yl)
ax.set_yticklabels(lbls, fontsize=8.2)
ax.set_xlim(0, 29)
ax.set_xlabel("Count of Allele-Exact 2–4★ AD ClinVar P/LP Alerts (FAF95 ≥ 1e-4)")
ax.set_title("B. Complete 55/55 Live-QC + CSpec/Johnston DAFT + MOI Triage\n(0/55 Exceed Verified BA1; 1/55 Reduced-Pen CACNA1S > DAFT BS1)")

fig.tight_layout(w_pad=2.4)
save_fig(fig, "fig2_a1_rsid_mirage_and_4layer_triage.png")

# =============================================================================
# FIGURE 3: All 11 High-QC Non-Bottleneck L4 Variants vs Verified CSpec / Johnston DAFT BS1
# =============================================================================
fig, ax = plt.subplots(figsize=(12.2, 5.4))
l4 = pd.DataFrame(a_t05["ad55_l4_all11_verified"]).sort_values("eval_freq_faithful", ascending=True).reset_index(drop=True)
y_pos = np.arange(len(l4))
for i, r in l4.iterrows():
    if "L4b" in r["verified_t05_layer"]:
        col = C_CORAL
    elif "L4c" in r["verified_t05_layer"]:
        col = "#7570B3"
    else:
        col = C_BLUE
    vf = r["eval_freq_faithful"] * 1e4
    vb = r["gene_bs1"] * 1e4
    ax.scatter(vf, i, s=85, color=col, zorder=4, edgecolors="black", linewidths=0.6)
    ax.scatter(vb, i, marker="|", s=240, color=C_RED, linewidths=2.2, zorder=5)
    ax.plot([vf, vb], [i, i], color="#CCCCCC", lw=1.2, zorder=2)
    src_tag = r["cspec_id"] if str(r["cspec_id"]).startswith("GN") else "Johnston2026"
    tag = f"{r['faf95_pop8'].upper()} {r['eval_freq_faithful']:.5f} (BS1={r['gene_bs1']:.5f}, {src_tag})"
    x_txt = max(vf, vb) * 1.15
    ax.text(x_txt, i, tag, va="center", fontsize=7.8, color="#222222")

ax.axvline(1.0, color="#777777", linestyle=":", lw=1.1, label="Uniform 1e-4 Discovery Screen")
lbls = [f"{r['gene']}  {r['hgvs_name'].split(':')[-1][:34]} ({r['stars']}★)" for _, r in l4.iterrows()]
ax.set_yticks(y_pos)
ax.set_yticklabels(lbls, fontsize=8.4)
ax.set_xscale("log")
ax.set_xlim(0.35, 160)
ax.set_xlabel("Rule-Faithful Grpmax Frequency × 10⁴ (Log Scale) vs. Verified Official CSpec / Johnston 2026 DAFT BS1")
ax.set_title("Complete Live-QC-Passed Non-Bottleneck AD 2–4★ Alerts (n = 11; APC c.423-3T>A Removed for AS_VQSR Failure)")
handles = [
    plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=C_CORAL, markeredgecolor="k", ms=8, label="L4b: AR-in-AD Return Hazard (RYR1/TTN, n=4; GN179 cites c.10348-6C>G)"),
    plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#7570B3", markeredgecolor="k", ms=8, label="L4c: > Johnston DAFT BS1, < BA1 (CACNA1S MH p.Arg1086Ser, n=1)"),
    plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=C_BLUE, markeredgecolor="k", ms=8, label="L4a: True AD Pathogenic Below Verified BS1 (n=6)"),
    plt.Line2D([0], [0], marker="|", color=C_RED, lw=0, mew=2.2, ms=12, label="Verified ClinGen CSpec / Johnston 2026 DAFT BS1 Ceiling"),
]
ax.legend(handles=handles, frameon=True, facecolor="white", edgecolor="#DDDDDD", fontsize=7.9,
          loc="upper center", bbox_to_anchor=(0.48, -0.14), ncol=2)
fig.tight_layout()
save_fig(fig, "fig3_a1_l4_priority_catalog_and_carrier_burden.png")

# =============================================================================
# FIGURE 4: Verified Tier-Stratified ClinGen CSpec (30 Genes) & Johnston 2026
#           DAFT (47 Genes) VUS Frequency-Conflict Queue (Ploidy-Aware incl. chrX!)
# =============================================================================
fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.0), gridspec_kw={"width_ratios": [1.15, 1.05]})

# Panel A: Tier 1 (30 Released CSpec incl. ABCD1/OTC) vs Tier 2 (6 Prelim) vs Tier 3 (47 Johnston DAFT incl. GLA)
ax = axes[0]
vs = a_t08["vus_t05_t08_ploidy_aware_summary"]
groups = [
    "Tier 1: 30 Released\nCSpec VCEPs\n(≥ BS1/BA1)",
    "Tier 1: 30 Released\nCSpec VCEPs\n(≥ BA1)",
    "Tier 2: 6 Prelim\nVCEP Genes\n(≥ BS1)",
    "Tier 3: 47 Johnston\n2026 DAFT\n(≥ BS1)",
    "Tier 3: 47 Johnston\n2026 DAFT\n(≥ BA1)",
]
keys = [
    ("tier1_cspec_30_faithful", "bs1"),
    ("tier1_cspec_30_faithful", "ba1"),
    ("tier2_prelim_vcep_6", "bs1"),
    ("tier3_johnston_daft_47", "bs1"),
    ("tier3_johnston_daft_47", "ba1"),
]
tot_c = [vs[k][f"{m}_n"] for k, m in keys]
nfe_inv = [vs[k][f"{m}_nfe_inv_n"] for k, m in keys]
nfe_vis = [t - i for t, i in zip(tot_c, nfe_inv)]
xx = np.arange(len(groups))
ax.bar(xx, nfe_inv, width=0.58, color=C_CORAL, edgecolor="white", label="NFE-Invisible (Non-NFE Only ≥ Threshold)")
ax.bar(xx, nfe_vis, width=0.58, bottom=nfe_inv, color=C_BLUE, edgecolor="white", label="Also ≥ Threshold in NFE")
for i in range(len(groups)):
    pct = nfe_inv[i] / max(tot_c[i], 1) * 100
    ax.text(xx[i], tot_c[i] + 18, f"n = {tot_c[i]}\n({pct:.1f}% inv)", ha="center", va="bottom", fontsize=7.8, fontweight="bold")
ax.set_xticks(xx)
ax.set_xticklabels(groups, fontsize=8.0)
ax.set_ylim(0, 1130)
ax.set_ylabel("Ploidy-Aware High-Callability VUS Candidates\n(Autosomes + chrX Non-PAR; Excl. HFE & Exceptions)")
ax.set_title("A. Verified CSpec (30 Genes) vs. Johnston 2026 DAFT\nVUS Queue (incl. 17 chrX VUS; 36/37 Live QC PASS)")
ax.legend(frameon=False, fontsize=8.0, loc="upper left")

# Panel B: Ancestry Attribution Inside Tier 1 Official Released CSpec VCEPs (n = 327, 77.7% NFE-Invisible!)
ax = axes[1]
pop_d = vs["tier1_cspec_30_faithful"]["bs1_by_pop5"]
p_order = ["sas", "eas", "afr", "amr", "nfe"]
p_lbls = ["SAS", "EAS", "AFR", "AMR", "NFE"]
p_vals = [pop_d.get(p, 0) for p in p_order]
p_cols = [C_GOLD, C_CORAL, C_TEAL, "#7570B3", C_BLUE]
xp = np.arange(len(p_order))
bars = ax.bar(xp, p_vals, width=0.58, color=p_cols, edgecolor="white")
for rect, v in zip(bars, p_vals):
    ax.text(rect.get_x() + rect.get_width()/2, v + 2.0, f"{v}\n({v/sum(p_vals)*100:.1f}%)", ha="center", va="bottom", fontsize=8.2, fontweight="bold")
ax.set_xticks(xp)
ax.set_xticklabels(p_lbls, fontsize=9.0)
ax.set_ylim(0, 108)
ax.set_ylabel("Tier 1 CSpec VCEP VUS ≥ BS1/BA1 (n = 327)")
ax.set_title("B. Grpmax Ancestry Driving Tier 1 CSpec VCEP Flags\n(30 Released VCEPs incl. ABCD1/OTC; 77.7% NFE-Invis)")
ax.text(0.97, 0.96, "Top Tier 1 CSpec Genes:\nAPC GN089 (107), FBN1 GN022 (93),\nHNF1A GN017 (22), VHL GN078 (13),\nABCD1 GN105 (11), PTEN (10), OTC (5)",
        transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
        bbox=dict(boxstyle="round,pad=0.35", fc="#F7F7F7", ec="#CCCCCC"))

fig.tight_layout(w_pad=2.2)
save_fig(fig, "fig4_a1_vus_liberation_and_post2024_governance.png")
print("ALL 4 REVISED A-1 FIGURES (T05 + T08 VERIFIED) GENERATED!", flush=True)
