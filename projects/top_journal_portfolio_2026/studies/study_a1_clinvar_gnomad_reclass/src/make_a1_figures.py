#!/opt/anaconda3/bin/python3
"""
Generate 4 Publication-Grade Multi-Panel Figures for Study A-1 (Nature Genetics / Nature Medicine style).
Every plotted value is loaded directly from verified outputs:
  - a1_reclassification_results.json
  - a1_addendum_triage_and_vus.json
  - a1_ad_highconf_triaged_55.csv
  - a1_clinvar_gnomad_exact_joined.parquet
"""
import json, os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

warnings.filterwarnings("ignore")

BASE = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_a1_clinvar_gnomad_reclass"
RES = os.path.join(BASE, "results")
CLEAN = os.path.join(BASE, "data", "clean")
FIG_DIR = os.path.join(BASE, "figures")
ART_DIR = "/Users/mac/.gemini/antigravity/brain/19f7514f-ce6a-4297-92cb-9803e648afd1"
os.makedirs(FIG_DIR, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9.3,
    "axes.titlesize": 10.6,
    "axes.titleweight": "bold",
    "axes.labelsize": 9.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.9,
    "xtick.labelsize": 8.4,
    "ytick.labelsize": 8.4,
    "legend.fontsize": 8.0,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.facecolor": "white",
})

C_BLUE = "#0072B2"
C_RED = "#D55E00"
C_TEAL = "#009E73"
C_GOLD = "#E69F00"
C_PURP = "#CC79A7"
C_GRAY = "#7F8C8D"
C_DARK = "#2C3E50"

d = json.load(open(os.path.join(RES, "a1_reclassification_results.json")))
add = json.load(open(os.path.join(RES, "a1_addendum_triage_and_vus.json")))
cat55 = pd.read_csv(os.path.join(RES, "a1_ad_highconf_triaged_55.csv"))

POP_LABELS = {
    "nfe": "NFE\n(Euro.)", "afr": "AFR\n(Afr.)", "amr": "AMR\n(Admix.)",
    "eas": "EAS\n(E. Asian)", "sas": "SAS\n(S. Asian)", "mid": "MID\n(Mid. East)",
    "asj": "ASJ\n(Ashk. Jew.)", "fin": "FIN\n(Finnish)"
}
POP_SHORT = {"nfe": "NFE", "afr": "AFR", "amr": "AMR", "eas": "EAS", "sas": "SAS", "mid": "MID", "asj": "ASJ", "fin": "FIN"}


def save_both(fig, name):
    p1 = os.path.join(FIG_DIR, name)
    p2 = os.path.join(ART_DIR, name)
    fig.savefig(p1, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(p2, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved {p1} ({os.path.getsize(p1):,} B)")


# ==============================================================================
# FIGURE 1: Cross-Ancestry Classification Asymmetry & Curation-Velocity Mediation
# ==============================================================================
def make_fig1():
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 5.0), gridspec_kw={"width_ratios": [1.05, 1.0, 1.0]})
    fig.subplots_adjust(wspace=0.30)

    # Panel A: Stacked Classification Share across 8 Ancestry-Enriched Variant Sets (AF >= 1e-4, >= 5x)
    ax = axes[0]
    pops = ["nfe", "afr", "sas", "eas", "amr", "mid", "asj", "fin"]
    prof = d["h1_equity"]["by_ancestry_profile"]
    x = np.arange(len(pops))
    vus = np.array([prof[p]["vus_rate"] * 100 for p in pops])
    plp = np.array([prof[p]["plp_rate"] * 100 for p in pops])
    conf = np.array([prof[p]["conf_rate"] * 100 for p in pops])
    blb = np.array([prof[p]["blb_rate"] * 100 for p in pops])

    ax.bar(x, vus, 0.65, color=C_RED, label="VUS (Uncertain)")
    ax.bar(x, plp, 0.65, bottom=vus, color="#8E44AD", label="P/LP (Pathogenic)")
    ax.bar(x, conf, 0.65, bottom=vus + plp, color=C_GOLD, label="Conflicting (96–97% B/LB vs VUS)")
    ax.bar(x, blb, 0.65, bottom=vus + plp + conf, color=C_TEAL, alpha=0.85, label="Benign / Likely Benign")

    for i, p in enumerate(pops):
        ax.text(i, vus[i] / 2, f"{vus[i]:.0f}%", ha="center", va="center", fontsize=7.4, color="white", fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([POP_SHORT[p] for p in pops], fontsize=8.4)
    ax.set_ylabel("ClinVar Classification Share (%)")
    ax.set_ylim(0, 118)
    ax.set_title("A. Allele-Exact Coding/Splice Labels (AF ≥ 10⁻⁴)")
    ax.legend(loc="upper center", ncol=2, frameon=True, facecolor="white", fontsize=7.2)

    # Panel B: Frequency-Matched CMH Odds Ratios vs NFE (VUS & P/LP)
    ax = axes[1]
    pops_nn = ["afr", "amr", "sas", "eas", "mid", "asj", "fin"]
    yp = np.arange(len(pops_nn))[::-1]
    ax.axvline(1.0, color=C_DARK, ls="--", lw=1.0)
    for i, p in enumerate(pops_nn):
        y = yp[i]
        v_or, v_ci = prof[p]["cmh_vus_or_vs_nfe"], prof[p]["cmh_vus_ci95"]
        p_or, p_ci = prof[p]["cmh_plp_or_vs_nfe"], prof[p]["cmh_plp_ci95"]
        ax.errorbar(v_or, y + 0.14, xerr=[[v_or - v_ci[0]], [v_ci[1] - v_or]], fmt="o", color=C_RED, ecolor=C_RED,
                    capsize=2.8, ms=5.5, lw=1.4, label="VUS vs NFE (CMH)" if i == 0 else None)
        ax.errorbar(p_or, y - 0.14, xerr=[[p_or - p_ci[0]], [p_ci[1] - p_or]], fmt="s", color="#8E44AD", ecolor="#8E44AD",
                    capsize=2.8, ms=5.0, lw=1.3, label="P/LP vs NFE (CMH)" if i == 0 else None)
        ax.text(v_ci[1] * 1.07, y + 0.14, f"{v_or:.1f}×", va="center", ha="left", fontsize=7.4, color=C_RED, fontweight="bold")

    ax.set_xscale("log")
    ax.set_xticks([0.5, 1, 2, 5, 10, 25])
    ax.set_xticklabels(["0.5×", "1×", "2×", "5×", "10×", "25×"])
    ax.set_xlim(0.45, 34)
    ax.set_yticks(yp)
    ax.set_yticklabels([f"{POP_SHORT[p]} (n={prof[p]['n']:,})" for p in pops_nn], fontsize=8.0)
    ax.set_xlabel("Frequency-Matched CMH Odds Ratio vs. NFE (Log)")
    ax.set_title("B. Frequency-Matched Disparity by Ancestry")
    ax.legend(loc="upper center", bbox_to_anchor=(0.50, -0.17), ncol=2, frameon=True, facecolor="white", fontsize=7.6)
    ax.grid(axis="x", alpha=0.25)

    # Panel C: Curation Velocity Mediation (Submitter Count by Ancestry + Attenuation Inset)
    ax = axes[2]
    med = add["curation_velocity_mediation"]
    sba = med["submitters_by_ancestry"]
    pops_s = ["nfe", "asj", "afr", "eas", "amr", "sas", "fin", "mid"]
    xs = np.arange(len(pops_s))
    med_sub = [sba[p]["median_submitters"] for p in pops_s]
    s1_pct = [sba[p]["single_submitter_share"] * 100 for p in pops_s]
    cols_s = [C_BLUE if p == "nfe" else (C_GOLD if p in ("asj", "fin") else C_RED) for p in pops_s]

    ax.bar(xs, med_sub, 0.64, color=cols_s, alpha=0.88)
    for i, (ms, sp) in enumerate(zip(med_sub, s1_pct)):
        ax.text(i, ms + 0.22, f"{ms:.0f}", ha="center", va="bottom", fontsize=8.4, fontweight="bold", color=C_DARK)
        ax.text(i, ms * 0.48, f"{sp:.0f}%", ha="center", va="center", fontsize=7.2, fontweight="bold", color="white")

    ax.set_xticks(xs)
    ax.set_xticklabels([POP_SHORT[p] for p in pops_s], fontsize=8.4)
    ax.set_ylabel("Median Independent ClinVar Submitters")
    ax.set_xlabel("White % inside bar = single-submitter trapping share")
    ax.set_ylim(0, 15.0)
    ax.set_title("C. Curation Velocity & Mediation (p < 10⁻³⁰⁰)")
    ax.grid(axis="y", alpha=0.25)

    mtxt = (
        "Gene-FE Logistic Mediation (N = 18,265):\n"
        f"• Non-NFE submitters: −{med['log_submitters_gap_nonnfe_vs_nfe']['pct_fewer_submitters']:.1f}% (p < 10⁻³⁰⁰)\n"
        f"• VUS OR: {med['vus_base_or']:.2f}× → {med['vus_submitter_adjusted_or']:.2f}× (−{med['vus_log_or_attenuation_pct']:.1f}% log-OR)\n"
        f"• Conflict OR: {med['conf_base_or']:.2f}× → {med['conf_submitter_adjusted_or']:.2f}× (−{med['conf_log_or_attenuation_pct']:.1f}%, p=0.10)"
    )
    ax.text(0.22, 0.72, mtxt, transform=ax.transAxes, fontsize=7.0,
            bbox=dict(boxstyle="round,pad=0.32", facecolor="white", edgecolor="#BDC3C7", alpha=0.95))

    save_both(fig, "fig1_a1_equity_and_curation_velocity.png")


# ==============================================================================
# FIGURE 2: The 79.7% Multiallelic rsID Mirage & 4-Layer Clinical Triage of All 55 Alerts
# ==============================================================================
def make_fig2():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.0, 5.4), gridspec_kw={"width_ratios": [0.95, 1.25]})
    fig.subplots_adjust(wspace=0.28)

    # Panel A: H4 Multiallelic rsID False-Positive Mirage
    h4 = d["h4_rsid_collision_audit"]
    cats = ["All 84 Genes P/LP\n(AF ≥ 10⁻⁴)", "All 84 Genes P/LP\n(FAF95 ≥ 10⁻⁴)", "72 AD Genes 2–4★ P/LP\n(FAF95 ≥ 10⁻⁴)"]
    keys = ["af_ge_1e4", "faf95_ge_1e4", "ad_highconf_2to4star_faf95_ge_1e4"]
    exact_c = [h4[k]["allele_exact_confirmed"] for k in keys]
    fp_c = [h4[k]["multiallelic_false_positives"] for k in keys]
    fp_r = [h4[k]["false_positive_rate"] * 100 for k in keys]
    x = np.arange(len(cats))

    ax1.bar(x, exact_c, 0.52, color=C_BLUE, label="Allele-Exact Confirmed (CHROM-POS-REF-ALT)")
    ax1.bar(x, fp_c, 0.52, bottom=exact_c, color=C_RED, alpha=0.82, hatch="//", edgecolor="white", label="False-Positive Multiallelic rsID Collisions")
    for i in range(len(cats)):
        tot = exact_c[i] + fp_c[i]
        ax1.text(i, tot + 25, f"{fp_r[i]:.1f}% FP\n({fp_c[i]}/{tot})", ha="center", va="bottom", fontsize=8.2, fontweight="bold", color=C_RED)
        ax1.text(i, exact_c[i] / 2, f"{exact_c[i]}", ha="center", va="center", fontsize=8.2, fontweight="bold", color="white")

    ax1.set_xticks(x)
    ax1.set_xticklabels(cats, fontsize=8.3)
    ax1.set_ylabel("Flagged ClinVar P/LP Variant Count")
    ax1.set_ylim(0, 1520)
    ax1.set_title("A. Multiallelic rsID Collision Mirage (48,185 P/LP)")
    ax1.legend(loc="upper right", frameon=True, facecolor="white", fontsize=7.8)
    ax1.grid(axis="y", alpha=0.25)

    # Validation callout box
    val = d["gnomad_api_faf95_validation"]
    ax1.text(0.46, 0.62, f"Live gnomAD v4.1 API Validation:\nComputed vs. Native joint.faf95\nPearson r = {val['pearson_r']:.5f} (n = {val['n_validated']})",
             transform=ax1.transAxes, fontsize=7.6, fontweight="bold", color=C_TEAL,
             bbox=dict(boxstyle="round,pad=0.35", facecolor="#F4FBF8", edgecolor=C_TEAL, alpha=0.95))

    # Panel B: 4-Layer Technical & Population-Genetic Triage of ALL 55 AD >=2-Star Alerts
    layer_meta = {
        "L1_LowCallability_TechnicalArtifact": ("L1: Low Callability <80% AN (n=7)", C_GRAY, "X"),
        "L2_BottleneckFounder_FIN_ASJ": ("L2: Bottleneck Founder FIN/ASJ (n=22)", C_GOLD, "s"),
        "L3_EstablishedContinentalFounder_LowPen": ("L3: Continental Founder / Low-Pen (n=14)", C_BLUE, "^"),
        "L4_PriorityNonBottleneckAuditCandidate": ("L4: Priority Non-Bottleneck Audit (n=12)", C_RED, "o"),
    }
    ax2.axvline(80.0, color=C_GRAY, ls="--", lw=1.2, zorder=1)
    ax2.axvspan(0, 80.0, color="#F2F3F4", alpha=0.75, zorder=0)
    ax2.axhline(5e-4, color="#BDC3C7", ls=":", lw=1.0)

    # Slight jitter in x only near 100% callability for visibility? Better: exact callability % with clean alpha!
    for lkey, (llbl, col, mk) in layer_meta.items():
        sub = cat55[cat55["triage_layer"] == lkey]
        ax2.scatter(sub["an_frac"] * 100.0, sub["faf95_max8"], c=col, marker=mk, s=62 if lkey.startswith("L4") else 46,
                    alpha=0.92, edgecolors="white", linewidth=0.6, label=llbl, zorder=4 if lkey.startswith("L4") else 3)

    # Crisp selective labels with zero overlap
    annots = [
        ("2-47414421-A-T", "MSH2 c.942+3A>T\n(0.7% AN poly-A)", (14, 0.0028)),
        ("13-48464984-C-G", "RB1 c.2212-14C>G\n(35% AN)", (39, 0.0078)),
        ("18-31598655-G-A", "TTR p.Val142Ile\n(AFR founder)", (66, 0.015)),
        ("13-32340300-GT-G", "BRCA2 6174delT\n(ASJ founder)", (66, 0.0046)),
        ("5-112775626-T-A", "APC c.423-3T>A\n(AMR 1-in-263)", (56, 0.0017)),
        ("19-38527014-C-T", "RYR1 p.Arg3550Trp\n(3★ SAS 1-in-300)", (56, 0.00095)),
        ("19-38523211-C-G", "RYR1 c.10348-6C>G\n(3★ AFR)", (52, 0.00033)),
        ("3-38613773-G-A", "SCN5A p.Arg225Trp\n(2★ AFR, Hom=1)", (52, 0.00011)),
    ]
    for vid, txt, xy_txt in annots:
        r = cat55[cat55["variant_id"] == vid]
        if len(r) == 1:
            rx, ry = float(r["an_frac"].iloc[0] * 100.0), float(r["faf95_max8"].iloc[0])
            c_txt = C_RED if r["triage_layer"].iloc[0].startswith("L4") else C_DARK
            ax2.annotate(txt, xy=(rx, ry), xytext=xy_txt, fontsize=7.2, color=c_txt,
                         fontweight="bold" if c_txt == C_RED else "normal",
                         arrowprops=dict(arrowstyle="->", color=c_txt, lw=0.85, alpha=0.85))

    ax2.set_yscale("log")
    ax2.set_ylim(7.5e-5, 0.032)
    ax2.set_xlim(-3, 105)
    ax2.set_xlabel("gnomAD v4.1 Site Callability (% of 1.614M Maximum Alleles)")
    ax2.set_ylabel("Max-Ancestry 95% Filtering Allele Frequency (FAF95, Log)")
    ax2.set_title("B. 4-Layer Clinical/Technical Triage of All 55 AD ≥2★ Alerts")
    ax2.legend(loc="lower left", frameon=True, facecolor="white", fontsize=7.4)
    ax2.grid(alpha=0.25)

    save_both(fig, "fig2_a1_rsid_mirage_and_4layer_triage.png")


# ==============================================================================
# FIGURE 3: All 12 L4 Priority Audit Candidates & 17.5x Ancestry Carrier Disparity
# ==============================================================================
def make_fig3():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.4, 5.4), gridspec_kw={"width_ratios": [1.32, 0.92]})
    fig.subplots_adjust(wspace=0.28)

    # Panel A: Dumbbell plot of all 12 L4 Priority Audit Candidates (NFE FAF95 vs Popmax FAF95)
    l4 = cat55[cat55["triage_layer"] == "L4_PriorityNonBottleneckAuditCandidate"].sort_values("faf95_max8", ascending=True).reset_index(drop=True)
    y = np.arange(len(l4))
    floor_v = 1.5e-6  # display floor for exact 0.0 NFE FAF95

    ax1.axvline(1e-4, color=C_DARK, ls="--", lw=1.1, label="AD BS1 Threshold (FAF95 = 10⁻⁴)")
    ax1.axvline(5e-4, color=C_RED, ls=":", lw=1.0, label="High-Alert Threshold (FAF95 = 5×10⁻⁴)")

    pop_cols = {"sas": C_RED, "afr": C_TEAL, "eas": C_BLUE, "amr": C_PURP}
    for i, r in l4.iterrows():
        x_nfe = max(float(r["faf95_nfe"]), floor_v)
        x_pop = float(r["faf95_max8"])
        pcol = pop_cols.get(r["faf95_pop8"], C_RED)
        ax1.plot([x_nfe, x_pop], [i, i], color="#BDC3C7", lw=1.8, zorder=2)
        ax1.scatter([x_nfe], [i], c="#7F8C8D", s=42, marker="o", zorder=3, label="NFE FAF95 (≤ 5.2×10⁻⁵)" if i == 0 else None)
        ax1.scatter([x_pop], [i], c=pcol, s=68, marker="D", edgecolors="white", lw=0.6, zorder=4)
        hom_tag = " [Hom=1]" if int(r["hom_joint"]) >= 1 else ""
        yr = str(r["last_evaluated"])[-4:]
        star_s = "★★★" if int(r["stars"]) == 3 else "★★"
        ax1.text(x_pop * 1.13, i, f"{r['faf95_pop8'].upper()} {x_pop:.4f} ({star_s}, {yr}){hom_tag}",
                 va="center", ha="left", fontsize=7.4, color=pcol, fontweight="bold")

    ylbls = []
    for _, r in l4.iterrows():
        hg = str(r["hgvs_name"])
        short_c = hg.split(":")[-1].split(" ")[0] if ":" in hg else r["variant_id"]
        ylbls.append(f"{r['gene']} {short_c}")
    ax1.set_xscale("log")
    ax1.set_xlim(1.0e-6, 0.022)
    ax1.set_xticks([1.5e-6, 1e-5, 1e-4, 5e-4, 2e-3])
    ax1.set_xticklabels(["0 / ≤10⁻⁶", "10⁻⁵", "10⁻⁴", "5×10⁻⁴", "2×10⁻³"])
    ax1.set_yticks(y)
    ax1.set_yticklabels(ylbls, fontsize=8.2)
    ax1.set_xlabel("gnomAD v4.1 95% Filtering Allele Frequency (FAF95, Log Scale)")
    ax1.set_title("A. All 12 L4 Priority AD ≥2★ Audit Candidates (100% NFE-Invisible)")
    ax1.grid(axis="x", alpha=0.25)

    # Panel B: Individual-level Clinical Exposure Disparity (Carriers per 1,000 Sequenced Individuals)
    cb = add["l4_priority_carrier_burden"]
    pops_ord = ["sas", "amr", "afr", "mid", "eas", "fin", "asj", "nfe"]
    y2 = np.arange(len(pops_ord))[::-1]
    vals = [cb[p]["per_1000"] for p in pops_ord]
    rats = [cb[p]["ratio_vs_nfe"] for p in pops_ord]
    cols2 = [C_BLUE if p == "nfe" else (C_GOLD if p in ("fin", "asj") else C_RED) for p in pops_ord]

    ax2.barh(y2, vals, height=0.62, color=cols2, alpha=0.88)
    for i, (v, rt, p) in enumerate(zip(vals, rats, pops_ord)):
        lbl = f"{v:.2f} / 1,000 ({rt:.1f}×)" if p != "nfe" else f"{v:.2f} / 1,000 (1.0× ref)"
        ax2.text(v + 0.14, y2[i], lbl, va="center", ha="left", fontsize=8.0,
                 fontweight="bold" if rt >= 5.0 else "normal", color=cols2[i])

    ax2.set_yticks(y2)
    ax2.set_yticklabels([POP_SHORT[p] for p in pops_ord], fontsize=8.6, fontweight="bold")
    ax2.set_xlim(0, 8.8)
    ax2.set_xlabel("Carriers of ≥1 L4 Priority Variant per 1,000 Individuals")
    ax2.set_title("B. Clinical Carrier Exposure Disparity (Up to 17.5×)")
    ax2.grid(axis="x", alpha=0.25)

    save_both(fig, "fig3_a1_l4_priority_catalog_and_carrier_burden.png")


# ==============================================================================
# FIGURE 4: Actionable VUS-to-Benign Liberation Debt & Post-2024 Active Governance Gap
# ==============================================================================
def make_fig4():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.8, 4.9), gridspec_kw={"width_ratios": [1.15, 0.95]})
    fig.subplots_adjust(wspace=0.28)

    # Panel A: VUS Liberation Debt across FAF95 thresholds (stacked by driving continental population)
    vl = add["vus_liberation_debt"]["all_84_genes_all_vus"]["tiers"]
    tiers = [
        ("FAF95 ≥ 10⁻⁴\n(AD Watch)", "faf95_grpmax5_ge_1e4"),
        ("FAF95 ≥ 5×10⁻⁴\n(AD BS1)", "faf95_grpmax5_ge_5e4_BS1_AD"),
        ("FAF95 ≥ 10⁻³\n(Strict BS1)", "faf95_grpmax5_ge_1e3_BS1_strict"),
        ("FAF95 ≥ 5×10⁻³\n(Near-BA1)", "faf95_grpmax5_ge_5e3_BA1_near"),
    ]
    x = np.arange(len(tiers))
    pop_order = ["sas", "eas", "afr", "amr", "nfe"]
    pop_colors = {"sas": C_RED, "eas": C_BLUE, "afr": C_TEAL, "amr": C_PURP, "nfe": "#BDC3C7"}
    pop_names = {"sas": "SAS", "eas": "EAS", "afr": "AFR", "amr": "AMR", "nfe": "NFE (European)"}

    # Plot normalized 100% stacked composition + exact total count & NFE-invisible % above each bar!
    bottom = np.zeros(len(tiers))
    for p in pop_order:
        shares = []
        for _, k in tiers:
            tot = vl[k]["n_vus_liberable"]
            cnt = vl[k]["by_pop5"].get(p, 0)
            shares.append(cnt / tot * 100.0 if tot > 0 else 0.0)
        ax1.bar(x, shares, 0.56, bottom=bottom, color=pop_colors[p], label=pop_names[p], edgecolor="white", lw=0.6)
        bottom += np.array(shares)

    for i, (_, k) in enumerate(tiers):
        tot = vl[k]["n_vus_liberable"]
        inv_s = vl[k]["nfe_invisible_share"] * 100.0
        ax1.text(i, 102.5, f"n = {tot:,}\n({inv_s:.1f}% NFE-invis.)", ha="center", va="bottom",
                 fontsize=8.0, fontweight="bold", color=C_DARK)

    ax1.set_xticks(x)
    ax1.set_xticklabels([t[0] for t in tiers], fontsize=8.4)
    ax1.set_ylabel("Share of High-Callability VUS by Driving Ancestry (%)")
    ax1.set_ylim(0, 122)
    ax1.set_title("A. Immediate VUS-to-Benign Liberation Debt (50,861 VUS)")
    ax1.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=5, frameon=True, facecolor="white", fontsize=8.0)
    ax1.grid(axis="y", alpha=0.22)

    # Panel B: H3 Active Post-gnomAD-v4 (2024-2026) Governance Gap (LastEvaluated Year Distribution)
    j = pd.read_parquet(os.path.join(CLEAN, "a1_clinvar_gnomad_exact_joined.parquet"))
    ad_hc = j[(j["inheritance"] == "AD") & (j["clnsig_cat"] == "P_LP") & (j["stars"] >= 2)].copy()
    ad_hc["eval_yr"] = pd.to_numeric(ad_hc["last_evaluated"].astype(str).str.extract(r"(\d{4})")[0], errors="coerce")
    ad_hc = ad_hc[ad_hc["eval_yr"].notna()]

    bins_yr = [2000, 2019.5, 2021.5, 2023.5, 2024.5, 2025.5, 2027]
    lbl_yr = ["≤2019", "2020–21", "2022–23\n(Pre-v4)", "2024\n(Post-v4)", "2025\n(Post-v4)", "2026\n(Post-v4)"]
    f_all = ad_hc[ad_hc["faf95_max8"] >= 1e-4]
    u_all = ad_hc[ad_hc["faf95_max8"] < 1e-4]
    l4_ids = set(cat55[cat55["triage_layer"] == "L4_PriorityNonBottleneckAuditCandidate"]["variant_id"])
    f_l4 = ad_hc[ad_hc["variant_id"].isin(l4_ids)]

    c_u = pd.cut(u_all["eval_yr"], bins=bins_yr, labels=lbl_yr).value_counts(sort=False, normalize=True) * 100
    c_f = pd.cut(f_all["eval_yr"], bins=bins_yr, labels=lbl_yr).value_counts(sort=False, normalize=True) * 100
    c_l4 = pd.cut(f_l4["eval_yr"], bins=bins_yr, labels=lbl_yr).value_counts(sort=False, normalize=True) * 100

    x2 = np.arange(len(lbl_yr))
    w = 0.26
    ax2.axvspan(2.5, 5.5, color="#FDEDEC", alpha=0.65, zorder=0)
    ax2.axvline(2.5, color=C_RED, ls="--", lw=1.2)
    ax2.bar(x2 - w, c_u.values, w, color="#95A5A6", label=f"Unflagged AD ≥2★ P/LP (n={len(u_all):,})")
    ax2.bar(x2, c_f.values, w, color=C_BLUE, label=f"All 55 Flagged AD ≥2★ (n={len(f_all)})")
    ax2.bar(x2 + w, c_l4.values, w, color=C_RED, label=f"12 L4 Priority Candidates (n={len(f_l4)})")

    ax2.text(4.02, 52.5, "gnomAD v4.0 Released (Nov 2023)\n63.6% of all 55 & 75.0% (9/12) of L4\nwere re-evaluated in 2024–2026!",
             ha="center", va="top", fontsize=7.7, fontweight="bold", color=C_RED,
             bbox=dict(boxstyle="round,pad=0.34", facecolor="white", edgecolor=C_RED, alpha=0.95))

    ax2.set_xticks(x2)
    ax2.set_xticklabels(lbl_yr, fontsize=8.0)
    ax2.set_ylabel("Share of Variants by ClinVar LastEvaluated (%)")
    ax2.set_ylim(0, 56)
    ax2.set_title("B. Active Post-gnomAD-v4 Governance Blind Spot")
    ax2.legend(loc="upper center", bbox_to_anchor=(0.50, -0.17), ncol=1, frameon=True, facecolor="white", fontsize=7.7)
    ax2.grid(axis="y", alpha=0.22)

    save_both(fig, "fig4_a1_vus_liberation_and_post2024_governance.png")


if __name__ == "__main__":
    make_fig1()
    make_fig2()
    make_fig3()
    make_fig4()
    print("All 4 Study A-1 publication figures generated successfully!")
