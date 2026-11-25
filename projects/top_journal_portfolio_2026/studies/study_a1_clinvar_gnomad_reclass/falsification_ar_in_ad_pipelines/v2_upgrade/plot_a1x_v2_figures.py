#!/opt/anaconda3/bin/python3
"""
Polished zero-overlap publication figures for A1-X v2:
  1. fig1_v2_unbiased_8pipeline_and_1000g.png
  2. fig2_v2_ukbb_opposite_direction_lipids.png
"""
import json, os, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

V2 = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_a1_clinvar_gnomad_reclass/falsification_ar_in_ad_pipelines/v2_upgrade"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 8.5,
    "axes.titlesize": 9.5,
    "axes.labelsize": 8.6,
    "xtick.labelsize": 8.0,
    "ytick.labelsize": 8.0,
    "legend.fontsize": 7.4,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 300,
    "savefig.dpi": 300,
})

C_NAVY = "#1F4E79"
C_BLUE = "#3B82C4"
C_TEAL = "#1B9E77"
C_RED = "#C0392B"
C_ORANGE = "#D95F02"
C_AMBER = "#E69F00"
C_PURPLE = "#6A3D9A"
C_GRAY = "#6B7280"
C_LGRAY = "#E5E7EB"

# ==============================================================================
# FIGURE 1: Unbiased 4,799-Allele 8-Pipeline Audit + 1000G Real Genomes + 3* RYR1
# ==============================================================================
bdf = pd.read_csv(os.path.join(V2, "unbiased_7427_benchmark_v2_airtight.csv"))
b80 = bdf[bdf["subset"] == "HIGHCONF_2to4star_GE80AN"].copy()
g1k_exec = json.load(open(os.path.join(V2, "g1k_2504_real_pipeline_execution_results.json")))["canonical_2504"]

fig = plt.figure(figsize=(14.4, 9.2))
gs = fig.add_gridspec(2, 2, height_ratios=[1.18, 1.0], width_ratios=[1.12, 1.0], hspace=0.42, wspace=0.30)

# --- Panel A (Top Left): Unbiased 8-Pipeline Benchmark (2-4*, >=80% AN, n=4,799) ---
axA = fig.add_subplot(gs[0, 0])
order = [
    "Exomiser (98f4e0b)",
    "AutoGVP (054e2b7)",
    "CPHI-DRAGEN (974860d)",
    "Viking Genes (3ea7abe)",
    "CPSR v2.4.0 (7f68f9d)",
    "SFtool (0d72249)",
    "DNA (f459b22)",
    "OwnGenome (c324b61)",
    "Proposed Tier-4 Patch",
]
labels_pretty = [
    "Exomiser  [Tier A]",
    "AutoGVP  [Tier A]",
    "CPHI-DRAGEN  [Tier A]",
    "Viking Genes  [Tier A]",
    "CPSR v2.4.0  [Tier A]",
    "SFtool  [Tier A]",
    "DNA  [Tier B]",
    "OwnGenome  [Tier B]",
    "Tier-4 Patch  [Fix]",
]
b80 = b80.set_index("pipeline").loc[order].reset_index()
y = np.arange(len(order))[::-1]
h = 0.36

axA.barh(y + h/2, b80["carrier_fdr_joint_pct"], height=h, color=C_BLUE, label="Global Joint Carrier FDR (%)")
axA.barh(y - h/2, b80["carrier_fdr_sas_pct"], height=h, color=C_RED, label="South Asian (SAS) Carrier FDR (%)")

axA.axhline(2.5, color=C_GRAY, linestyle="--", linewidth=0.8, alpha=0.65)
axA.axhline(0.5, color=C_TEAL, linestyle="-", linewidth=1.0, alpha=0.85)

for i, row in b80.iterrows():
    yi = y[i]
    val_j = row["carrier_fdr_joint_pct"]
    val_s = row["carrier_fdr_sas_pct"]
    fp_n = int(row["FP"])
    spec = row["non_sf_spec_pct"]
    if row["pipeline"] == "Proposed Tier-4 Patch":
        axA.text(0.35, yi, "0.00% FDR  |  0/179 non-SF emitted  (100% recall, 100% spec.)",
                 va="center", ha="left", fontsize=7.6, fontweight="bold", color=C_TEAL)
    else:
        axA.text(val_j + 0.30, yi + h/2, f"{val_j:.2f}%", va="center", ha="left", fontsize=7.0, color=C_NAVY)
        axA.text(val_s + 0.30, yi - h/2, f"{val_s:.2f}%  ({fp_n}/179 FP; spec {spec:.1f}%)",
                 va="center", ha="left", fontsize=7.0, color=C_RED, fontweight="bold" if i in (0, 7) else "normal")

axA.set_yticks(y)
axA.set_yticklabels(labels_pretty, fontsize=8.0)
axA.set_xlim(0, 26.0)
axA.set_ylim(-0.65, 9.25)
axA.set_xlabel("Carrier-Weighted False Discovery Rate Among Emitted 2–4★ AD SF Calls (%)")
axA.set_title("A   Unbiased Benchmark (4,799 High-Callability 2–4★ AD Alleles)", loc="left", fontweight="bold", pad=22)
axA.legend(loc="lower center", bbox_to_anchor=(0.48, 1.005), ncol=2, frameon=True, facecolor="white", edgecolor=C_LGRAY, fontsize=7.4)
axA.grid(axis="x", linestyle=":", alpha=0.45)

# --- Panel B (Top Right): Real 1000G (N=2,504) vs gnomAD v4.1 (N=807,162) ---
axB = fig.add_subplot(gs[0, 1])
pops = ["EUR/NFE", "AMR", "EAS", "AFR", "SAS"]
pop_keys_gn = ["nfe", "amr", "eas", "afr", "sas"]
pop_keys_1k = ["EUR", "AMR", "EAS", "AFR", "SAS"]

row_exo = b80[b80["pipeline"] == "Exomiser (98f4e0b)"].iloc[0]
row_og = b80[b80["pipeline"] == "OwnGenome (c324b61)"].iloc[0]
gn_exo = [row_exo[f"carrier_fp_{k}"] for k in pop_keys_gn]
gn_og = [row_og[f"carrier_fp_{k}"] for k in pop_keys_gn]
g1k_hc24 = [g1k_exec["Exomiser_AutoGVP"]["by_superpop_primary"][k]["hc24_fp_per1k"] for k in pop_keys_1k]
g1k_hc24_n = [g1k_exec["Exomiser_AutoGVP"]["by_superpop_primary"][k]["hc24_fp"] for k in pop_keys_1k]
g1k_hc24_tot = [g1k_exec["Exomiser_AutoGVP"]["by_superpop_primary"][k]["hc24_emitted"] for k in pop_keys_1k]

x = np.arange(len(pops))
w = 0.25
axB.bar(x - w, gn_exo, width=w, color=C_NAVY, label="gnomAD v4.1 UB: Exomiser/AutoGVP")
axB.bar(x, gn_og, width=w, color=C_BLUE, label="gnomAD v4.1 UB: OwnGenome")
axB.bar(x + w, g1k_hc24, width=w, color=C_ORANGE, edgecolor="#7F2B00", linewidth=0.7, label="1000G Real (N=2,504 unrelated)")
axB.scatter([3 - w/2, 4 - w/2], [0.91, 3.34], marker="D", s=34, color="#FFD700", edgecolor="black", linewidth=0.8, zorder=6, label="Single-Variant Lower Bound (gnomAD)")

for xi, (ge, go, r1k, nfp, ntot) in enumerate(zip(gn_exo, gn_og, g1k_hc24, g1k_hc24_n, g1k_hc24_tot)):
    axB.text(xi - w, ge + 0.12, f"{ge:.2f}", ha="center", va="bottom", fontsize=6.7, color=C_NAVY)
    axB.text(xi, go + 0.12, f"{go:.2f}", ha="center", va="bottom", fontsize=6.7, color=C_BLUE)
    if nfp > 0:
        axB.text(xi + w, r1k + 0.14, f"{r1k:.2f}\n({nfp}/{ntot})", ha="center", va="bottom", fontsize=6.8, color="#8C2D04", fontweight="bold")
    else:
        axB.text(xi + w, 0.10, f"0/{ntot}", ha="center", va="bottom", fontsize=6.5, color=C_GRAY)

# Put gray info box in clean whitespace above EUR/AMR/EAS (y=3.1..5.3)
axB.text(-0.38, 3.25, "1000G N=2,504 unrelated (2–4★):\n• 0 multi-allele overlaps across 3,202\n• Primary AD: 4/49 calls (8.2% FDR)\n• EUR primary=0/503 (+2/503 if\n  NTHL1→TSC2 GENEINFO bleed)",
         fontsize=6.8, color="#374151", va="bottom",
         bbox=dict(boxstyle="round,pad=0.30", facecolor="#F9FAFB", edgecolor="#D1D5DB", linewidth=0.8))

# Put red SAS callout cleanly at top-middle (y=6.2..7.5) pointing right to SAS bar (4+w, 6.13)
axB.annotate("3/6 (50.0%) of SAS 2–4★ AD calls\nare 3★ AR RYR1 p.Arg3550Trp",
             xy=(4 + w - 0.05, 6.25), xytext=(0.40, 6.45),
             fontsize=7.2, fontweight="bold", color=C_RED,
             bbox=dict(boxstyle="round,pad=0.28", facecolor="#FFF5F5", edgecolor=C_RED, linewidth=0.8),
             arrowprops=dict(arrowstyle="->", color=C_RED, lw=1.0, connectionstyle="arc3,rad=-0.08"))

axB.set_xticks(x)
axB.set_xticklabels(pops, fontsize=8.2)
axB.set_xlim(-0.55, 4.65)
axB.set_ylim(0, 7.65)
axB.set_ylabel("Non-SF Calls per 1,000 Individuals / Diploid Equiv.")
axB.set_title("B   Ancestry Disparity: Real 1000G Genomes vs. gnomAD v4.1", loc="left", fontweight="bold", pad=30)
axB.legend(loc="lower center", bbox_to_anchor=(0.48, 1.005), ncol=2, frameon=True, facecolor="white", edgecolor=C_LGRAY, fontsize=6.8, columnspacing=0.9)
axB.grid(axis="y", linestyle=":", alpha=0.45)

# --- Panel C (Bottom Left): Within-Gene Discordant Share (2-4*, AN-normalized) ---
axC = fig.add_subplot(gs[1, 0])
genes_c = ["DES\n(AR null)", "CACNA1S\n(AR null)", "RYR1\n(AR myop.)", "APOB\n(Low-LDL pLoF)", "TTN\n(Low-PSI/nonTV)"]
shares_c = [85.6, 55.6, 45.2, 37.5, 15.0]
raw_lbls = ["85.6%  (48/59 raw inst.)", "55.6%  (47/101)", "45.2%  (835/1,882)", "37.5%  (429/1,213; opposite sign!)", "15.0%  (348/2,125)"]
cols_c = [C_PURPLE, C_PURPLE, C_RED, C_ORANGE, C_AMBER]
yc = np.arange(len(genes_c))[::-1]
axC.barh(yc, shares_c, height=0.58, color=cols_c, alpha=0.90)
for yi, val, lbl in zip(yc, shares_c, raw_lbls):
    axC.text(val + 1.5, yi, lbl, va="center", ha="left", fontsize=7.6, fontweight="bold" if "APOB" in lbl or "RYR1" in lbl else "normal")
axC.set_yticks(yc)
axC.set_yticklabels(genes_c, fontsize=8.1)
axC.set_xlim(0, 134)
axC.set_xlabel("Discordant Share of AN-Normalized 2–4★ ClinVar P/LP Rate (%)")
axC.set_title("C   Within-Gene Discordant Share in Dual-Mode AD SF Genes", loc="left", fontweight="bold", pad=10)
axC.grid(axis="x", linestyle=":", alpha=0.45)

# --- Panel D (Bottom Right): Exhaustive 3* ClinGen Expert-Panel RYR1 Ancestry Breakdown ---
axD = fig.add_subplot(gs[1, 1])
ryr_pops = ["SAS\n(n=184)", "AFR\n(n=68)", "Global\n(n=1,377)", "AMR\n(n=19)", "NFE\n(n=994)", "EAS\n(n=11)"]
ar_pct = [83.2, 57.4, 30.8, 21.1, 19.4, 18.2]
ad_pct = [100.0 - p for p in ar_pct]
xd = np.arange(len(ryr_pops))
axD.bar(xd, ar_pct, width=0.58, color=C_RED, label="Congenital Myopathies VCEP (AR-only; 10 alleles)")
axD.bar(xd, ad_pct, bottom=ar_pct, width=0.58, color=C_TEAL, alpha=0.82, label="MHS VCEP (AD Malignant Hyperthermia; 59 alleles)")
for xi, ap in enumerate(ar_pct):
    axD.text(xi, ap / 2, f"{ap:.1f}%", ha="center", va="center", fontsize=7.8, fontweight="bold", color="white")
axD.set_xticks(xd)
axD.set_xticklabels(ryr_pops, fontsize=7.9)
axD.set_ylim(0, 121)
axD.set_yticks([0, 20, 40, 60, 80, 100])
axD.set_ylabel("Share of All 3★ Expert-Panel RYR1 Instances (%)")
axD.set_title("D   Exhaustive Live XML Audit of All 69 3★ RYR1 Alleles", loc="left", fontweight="bold", pad=28)
axD.text(2.5, 103.5, "9/10 AR-only 3★ RYR1 alleles (418/424 instances = 98.6%)\ncarry polluted OMIM:145600 MH tags in clinvar.vcf.gz",
         ha="center", va="bottom", fontsize=7.2, fontweight="bold", color=C_RED,
         bbox=dict(boxstyle="round,pad=0.25", facecolor="#FFF5F5", edgecolor=C_RED, linewidth=0.75))
axD.legend(loc="lower center", bbox_to_anchor=(0.48, 1.005), ncol=1, frameon=True, facecolor="white", edgecolor=C_LGRAY, fontsize=7.1)

fig.savefig(os.path.join(V2, "fig1_v2_unbiased_8pipeline_and_1000g.png"), bbox_inches="tight")
plt.close(fig)
print("Saved polished fig1_v2_unbiased_8pipeline_and_1000g.png")

# ==============================================================================
# FIGURE 2: DIRECT UK BIOBANK PHENOTYPE PROOF (BACKMAN ET AL. NATURE 2021, N=411k)
# ==============================================================================
mdf = pd.read_csv(os.path.join(V2, "backman2021_ukbb_exact_apob_ldlr_clinvar_matches.csv"))
meta = json.load(open(os.path.join(V2, "backman2021_ukbb_meta_summary.json")))

ldl = mdf[(mdf["trait"] == "Direct_LDL_UKB30780") & (mdf["stars"] >= 2)].copy()
apob_plof = ldl[(ldl["gene"] == "APOB") & (ldl["gold_cat"] == "NOT_MEETING_APOB_FHBL_LOWLDL")].sort_values("alt_beta").copy()
apob_fh = ldl[(ldl["gene"] == "APOB") & (ldl["gold_cat"] == "MEETS_ACMG_SF")].sort_values("alt_beta", ascending=False).copy()
ldlr_top = ldl[ldl["gene"] == "LDLR"].sort_values("p_value").head(6).sort_values("alt_beta", ascending=False).copy()

def short_hgvs(gene, h):
    m = re.search(r"\((p\.[^)]+)\)", str(h))
    if m:
        return f"{gene} {m.group(1)}"
    m2 = re.search(r":(c\.[^\s(]+)", str(h))
    return f"{gene} {m2.group(1)}" if m2 else f"{gene} {h}"

fig2 = plt.figure(figsize=(14.4, 8.8))
gs2 = fig2.add_gridspec(2, 2, height_ratios=[1.12, 1.0], width_ratios=[1.32, 1.0], hspace=0.40, wspace=0.34)

# --- Panel A (Full Left Column): Variant-Level UKBB Direct LDL-C Forest Plot ---
ax1 = fig2.add_subplot(gs2[:, 0])
plot_items = []
for _, r in apob_plof.iterrows():
    plot_items.append((short_hgvs("APOB", r["hgvs"]), r["alt_beta"], r["se"], r["p_value"], int(r["gnomad_het"]), "APOB_PLOF"))
ivw_p = meta["Direct_LDL_UKB30780"]["HIGHCONF_2to4star"]["APOB_pLoF_FHBL"]
plot_items.append(("IVW Mean: All 13 APOB pLoF (2–4★)", ivw_p["ivw_beta_sd"], ivw_p["ivw_se"], float(ivw_p["ivw_p"]), -1, "IVW_PLOF"))
plot_items.append(("SEP", np.nan, np.nan, np.nan, 0, "SEP"))
for _, r in apob_fh.iterrows():
    plot_items.append((short_hgvs("APOB", r["hgvs"]) + " (AD FH)", r["alt_beta"], r["se"], r["p_value"], int(r["gnomad_het"]), "APOB_FH"))
for _, r in ldlr_top.iterrows():
    plot_items.append((short_hgvs("LDLR", r["hgvs"]) + " (AD FH)", r["alt_beta"], r["se"], r["p_value"], int(r["gnomad_het"]), "LDLR_FH"))
ivw_l = meta["Direct_LDL_UKB30780"]["HIGHCONF_2to4star"]["LDLR_TrueFH"]
plot_items.append(("IVW Mean: All 84 LDLR FH (2–4★)", ivw_l["ivw_beta_sd"], ivw_l["ivw_se"], float(ivw_l["ivw_p"]), -1, "IVW_FH"))

y1 = np.arange(len(plot_items))[::-1]
ax1.axvline(0, color="black", linestyle="--", linewidth=1.0, zorder=1)
ax1.axvspan(-4.1, 0, color="#EFF6FF", alpha=0.55, zorder=0)
ax1.axvspan(0, 2.7, color="#FEF2F2", alpha=0.50, zorder=0)

ytick_pos, ytick_lbl = [], []
for idx, (lbl, beta, se, pv, het, kind) in enumerate(plot_items):
    yi = y1[idx]
    if kind == "SEP":
        ax1.axhline(yi, color=C_GRAY, linestyle="-", linewidth=0.95, alpha=0.75)
        continue
    ytick_pos.append(yi)
    ytick_lbl.append(lbl)
    ci = 1.96 * se
    if kind == "APOB_PLOF":
        col = C_NAVY if pv < 0.05 else "#7AA6D2"
        ax1.errorbar(beta, yi, xerr=ci, fmt="o", color=col, ecolor=col, elinewidth=1.35, capsize=2.5, markersize=5.2, zorder=4)
        if pv < 1e-4:
            # Place p-value cleanly on the left without hitting frame
            ax1.text(beta - ci - 0.08, yi, f"p={pv:.1e}", va="center", ha="right", fontsize=6.6, color=C_NAVY, fontweight="bold")
    elif kind == "IVW_PLOF":
        ax1.errorbar(beta, yi, xerr=ci, fmt="D", color="#0B2545", ecolor="#0B2545", elinewidth=2.0, capsize=3.5, markersize=7.2, zorder=5)
        ax1.text(beta - ci - 0.10, yi, f"{beta:+.2f} SD (p=4.0e-194)", va="center", ha="right", fontsize=7.2, fontweight="bold", color="#0B2545")
    elif kind == "APOB_FH":
        ax1.errorbar(beta, yi, xerr=ci, fmt="s", color=C_RED, ecolor=C_RED, elinewidth=1.4, capsize=2.5, markersize=5.5, zorder=4)
        if pv < 1e-4:
            ax1.text(beta + ci + 0.08, yi, f"p={pv:.1e}", va="center", ha="left", fontsize=6.7, color=C_RED, fontweight="bold")
    elif kind == "LDLR_FH":
        ax1.errorbar(beta, yi, xerr=ci, fmt="^", color=C_ORANGE, ecolor=C_ORANGE, elinewidth=1.35, capsize=2.5, markersize=5.4, zorder=4)
        if pv < 1e-4:
            ax1.text(beta + ci + 0.08, yi, f"p={pv:.1e}", va="center", ha="left", fontsize=6.6, color="#9A3412")
    elif kind == "IVW_FH":
        ax1.errorbar(beta, yi, xerr=ci, fmt="D", color=C_RED, ecolor=C_RED, elinewidth=2.0, capsize=3.5, markersize=7.2, zorder=5)
        ax1.text(beta + ci + 0.09, yi, f"{beta:+.2f} SD (p=6.3e-37)*", va="center", ha="left", fontsize=7.2, fontweight="bold", color=C_RED)

ax1.set_yticks(ytick_pos)
ax1.set_yticklabels(ytick_lbl, fontsize=7.6)
for tick in ax1.get_yticklabels():
    if "IVW Mean" in tick.get_text():
        tick.set_fontweight("bold")

ax1.set_xlim(-4.05, 2.75)
ax1.set_ylim(-0.9, len(plot_items) + 1.55)
ax1.set_xlabel("Per-Allele Effect on UK Biobank Direct LDL-C (SD units, ±95% CI; N = 411,351)")
ax1.set_title("A   UK Biobank Direct LDL-C (Backman et al., Nature 2021)", loc="left", fontweight="bold", pad=10)

# Dedicated clean header band ABOVE top data point with zero cross-zero collision!
ax1.text(-3.95, len(plot_items) + 0.85, "← LOWER LDL-C (FHBL: 13 APOB pLoFs)",
         fontsize=7.3, fontweight="bold", color=C_NAVY, va="center")
ax1.text(0.10, len(plot_items) + 0.85, "HIGHER LDL-C (True AD FH) →",
         fontsize=7.3, fontweight="bold", color=C_RED, va="center")

# Put contrast callout box cleanly in the completely empty bottom-left whitespace (x=-3.85..-0.35, y=1.5..5.8)!
ax1.text(-3.82, 4.0,
         "Direct Within-APOB Contrast\n(Same UKB Cohort, N = 411,351):\n"
         "• Misreturned APOB pLoFs (n=13):\n   IVW = −1.73 ± 0.06 SD (p = 4.0e-194)\n"
         "• True AD FH APOB missense (n=2):\n   IVW = +0.98 ± 0.05 SD (p = 5.0e-79)\n"
         "• Opposite-sign gap: Δβ = −2.70 SD\n   (Z = −34.7, Mann-Whitney p = 7.3e-6)",
         fontsize=7.2, color="#1F2937", va="center", ha="left",
         bbox=dict(boxstyle="round,pad=0.38", facecolor="white", edgecolor=C_RED, linewidth=1.0, alpha=0.96))
ax1.grid(axis="x", linestyle=":", alpha=0.45)

# --- Panel B (Top Right): Three-Biomarker Harmonized UKBB IVW Contrast ---
ax2 = fig2.add_subplot(gs2[0, 1])
traits_b = ["Direct LDL-C\n(N = 411,351)", "Total Cholesterol\n(N = 412,075)", "Apolipoprotein B\n(N = 410,113)"]
tkeys_b = ["Direct_LDL_UKB30780", "Total_Cholesterol_UKB30690", "ApoB_UKB30640"]
b_plof = [meta[k]["HIGHCONF_2to4star"]["APOB_pLoF_FHBL"]["ivw_beta_sd"] for k in tkeys_b]
se_plof = [1.96 * meta[k]["HIGHCONF_2to4star"]["APOB_pLoF_FHBL"]["ivw_se"] for k in tkeys_b]
b_afh = [meta[k]["HIGHCONF_2to4star"]["APOB_Missense_TrueFH"]["ivw_beta_sd"] for k in tkeys_b]
se_afh = [1.96 * meta[k]["HIGHCONF_2to4star"]["APOB_Missense_TrueFH"]["ivw_se"] for k in tkeys_b]
b_lfh = [meta[k]["HIGHCONF_2to4star"]["LDLR_TrueFH"]["ivw_beta_sd"] for k in tkeys_b]
se_lfh = [1.96 * meta[k]["HIGHCONF_2to4star"]["LDLR_TrueFH"]["ivw_se"] for k in tkeys_b]

xb = np.arange(len(traits_b))
wb = 0.25
ax2.axhline(0, color="black", linewidth=0.9)
ax2.bar(xb - wb, b_plof, yerr=se_plof, width=wb, color=C_NAVY, capsize=3, label="Misreturned APOB pLoF (n=12–13)")
ax2.bar(xb, b_afh, yerr=se_afh, width=wb, color=C_RED, capsize=3, label="True FH: APOB Missense (n=2)")
ax2.bar(xb + wb, b_lfh, yerr=se_lfh, width=wb, color=C_ORANGE, capsize=3, label="True FH: LDLR P/LP (n=83–84)*")

for xi, (vp, va, vl) in enumerate(zip(b_plof, b_afh, b_lfh)):
    ax2.text(xi - wb, vp - 0.16, f"{vp:+.2f}", ha="center", va="top", fontsize=7.4, fontweight="bold", color=C_NAVY)
    ax2.text(xi, va + 0.14, f"{va:+.2f}", ha="center", va="bottom", fontsize=7.4, fontweight="bold", color=C_RED)
    ax2.text(xi + wb, vl + 0.09, f"{vl:+.2f}", ha="center", va="bottom", fontsize=7.2, color="#9A3412")

ax2.set_xticks(xb)
ax2.set_xticklabels(traits_b, fontsize=8.0)
ax2.set_ylim(-2.25, 1.50)
ax2.set_ylabel("IVW Mean Effect Size (SD, ±95% CI)")
ax2.set_title("B   Harmonized Opposite-Sign Effect Across 3 Biomarkers", loc="left", fontweight="bold", pad=34)
ax2.legend(loc="lower center", bbox_to_anchor=(0.48, 1.008), ncol=2, frameon=True, facecolor="white", edgecolor=C_LGRAY, fontsize=6.9)
ax2.grid(axis="y", linestyle=":", alpha=0.45)

# --- Panel C (Bottom Right): Clinical Outcomes, Medication Use & NMR Signature ---
ax3 = fig2.add_subplot(gs2[1, 1])
clin_labels = [
    "LDLR p.Glu101Lys: Ezetimibe",
    "LDLR p.Pro685Leu: ICD-10 E78.0",
    "LDLR p.Asp482Asn: Lipid meds",
    "APOB p.Arg3527Gln: High chol.",
    "APOB p.Arg3527Gln: Lipid meds",
    "APOB p.Arg3527Gln: Aortic sten.",
    "APOB p.Gly3617fs: Bilirubin (SD)**",
    "APOB p.Arg439Ter: HDL chol. (SD)**",
    "APOB p.Gly3617fs: VLDL lip. (SD)**",
]
clin_beta = [4.232, 2.952, 2.205, 2.181, 1.818, 1.209, 0.712, -1.533, -1.594]
clin_se   = [0.599, 0.393, 0.251, 0.135, 0.138, 0.156, 0.114, 0.133, 0.111]
clin_note = [
    "OR=68.8, p=1.6e-12",
    "OR=19.1, p=6.3e-14",
    "OR=9.1, p=1.8e-18",
    "OR=8.9, p=2.1e-58",
    "OR=6.2, p=8.4e-40",
    "OR=3.4, p=1.0e-14",
    "p=4.2e-10 (hepatic)",
    "p=7.6e-31 (NMR)",
    "p=1.5e-46 (NMR)",
]
clin_col  = [C_ORANGE, C_ORANGE, C_ORANGE, C_RED, C_RED, C_RED, C_PURPLE, C_NAVY, C_NAVY]
yc3 = np.arange(len(clin_labels))[::-1]
ax3.axvline(0, color="black", linestyle="--", linewidth=0.9)
ax3.axhline(2.5, color=C_GRAY, linestyle="-", linewidth=0.8, alpha=0.65)
for i in range(len(clin_labels)):
    ax3.errorbar(clin_beta[i], yc3[i], xerr=1.96*clin_se[i], fmt="o", color=clin_col[i], ecolor=clin_col[i], capsize=2.3, markersize=4.8, zorder=4)
    if clin_beta[i] >= 0:
        ax3.text(clin_beta[i] + 1.96*clin_se[i] + 0.14, yc3[i], clin_note[i], va="center", ha="left", fontsize=6.7, color=clin_col[i])
    else:
        # Place cleanly to the right of zero line with crisp halo so zero-line never strikes text!
        ax3.text(0.18, yc3[i], f"{clin_beta[i]:+.2f} SD, {clin_note[i]}", va="center", ha="left", fontsize=6.7, color=clin_col[i], fontweight="bold")
ax3.set_yticks(yc3)
ax3.set_yticklabels(clin_labels, fontsize=7.2)
ax3.set_xlim(-2.2, 8.2)
ax3.set_xlabel("Log-OR (Clinical Outcomes/Meds) or SD (**Biomarkers/NMR)")
ax3.set_title("C   Divergent Clinical Sequelae & NMR Biomarker Profiles", loc="left", fontweight="bold", pad=10)
ax3.grid(axis="x", linestyle=":", alpha=0.45)

fig2.text(0.01, -0.02,
          "*Unadjusted UKB LDLR LDL-C estimates are attenuated by high lipid-lowering treatment uptake among LDLR carriers (Panel C, OR 9–69×). "
          "**Bottom 3 rows of Panel C show SD units from Backman et al. (2021) and Zoodsma et al. (Nat Genet 2025, N=434,646); top 6 rows show log-OR (UKB / Small et al. 2026).",
          fontsize=7.2, color="#4B5563")

fig2.savefig(os.path.join(V2, "fig2_v2_ukbb_opposite_direction_lipids.png"), bbox_inches="tight")
plt.close(fig2)
print("Saved polished fig2_v2_ukbb_opposite_direction_lipids.png")
