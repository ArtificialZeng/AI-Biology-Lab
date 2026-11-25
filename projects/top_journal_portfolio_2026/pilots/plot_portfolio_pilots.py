#!/usr/bin/env python3
"""Generate publication-grade 6-panel empirical pilot evidence figure for Top-Journal Portfolio 2026."""

from __future__ import annotations

import json
import shutil
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PORT = Path("/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026")
FIG_PATH = PORT / "figures" / "fig_portfolio_6pilots_evidence.png"
ART_PATH = Path("/Users/mac/.gemini/antigravity/brain/19f7514f-ce6a-4297-92cb-9803e648afd1/fig_portfolio_6pilots_evidence.png")

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9.5,
    "axes.titlesize": 10.5,
    "axes.titleweight": "bold",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 220,
})


def main() -> None:
    a1 = json.loads((PORT / "pilots" / "pilot_a1_summary.json").read_text())
    a2 = json.loads((PORT / "pilots" / "pilot_a2_summary.json").read_text())
    a5 = json.loads((PORT / "pilots" / "pilot_a5_summary.json").read_text())
    b1 = json.loads((PORT / "pilots" / "pilot_b1_summary.json").read_text())
    b2 = json.loads((PORT / "pilots" / "pilot_b2_summary.json").read_text())
    b2_sema = json.loads((PORT / "pilots" / "pilot_b2_semaglutide_addendum.json").read_text())
    b3_raw = json.loads((PORT / "data" / "pilot_b3" / "drug_shortages_all.json").read_text())["results"]

    fig, axes = plt.subplots(2, 3, figsize=(17.5, 10.2), constrained_layout=True)
    fig.suptitle(
        "Empirical Feasibility Pilots Across 6 Unstudied Top-Journal Research Problems (Live Public Data, Nov 2026)",
        fontsize=13.5, fontweight="bold", y=1.02,
    )

    # --- Panel A (Track A #1): Pilot A-2 ClinVar Reclassification Debt (gnomAD v4.1) ---
    ax = axes[0, 0]
    cats = ["Unresolved\nVUS", "Persistent\nP / LP", "Resolved\nBenign (B/LB)"]
    ft = a2["fisher_tests"]
    non_nfe = [ft["vus_among_clinvar"]["non_nfe_rate"] * 100,
               ft["plp_among_clinvar"]["non_nfe_rate"] * 100,
               ft["benign_resolved_among_clinvar"]["non_nfe_rate"] * 100]
    nfe = [ft["vus_among_clinvar"]["nfe_rate"] * 100,
           ft["plp_among_clinvar"]["nfe_rate"] * 100,
           ft["benign_resolved_among_clinvar"]["nfe_rate"] * 100]
    x = np.arange(len(cats))
    w = 0.36
    b1_bars = ax.bar(x - w/2, non_nfe, w, label="Non-European enriched (N=2,407)", color="#C0392B")
    b2_bars = ax.bar(x + w/2, nfe, w, label="European enriched (N=76)", color="#2E86C1")
    for bars in (b1_bars, b2_bars):
        for bar in bars:
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2, h + 1.0, f"{h:.1f}%", ha="center", va="bottom", fontsize=8.5)
    ax.set_xticks(x)
    ax.set_xticklabels(cats)
    ax.set_ylabel("% of Matched ClinVar Records")
    ax.set_ylim(0, 49)
    ax.set_title("A1 [Nat. Genet.]: ClinVar Reclassification Debt\n(8 ACMG Genes, gnomAD v4.1 Matched AF >=1e-4)")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.text(0.28, 0.45, "VUS OR=3.42 (p=1.2e-4)\nP/LP OR=7.39 (p=0.012)\nBenign OR=0.34 (p=2.5e-5)\n7/18 exact-allele P/LP",
            transform=ax.transAxes, fontsize=7.8, bbox=dict(boxstyle="round,pad=0.3", fc="#FDEDEC", ec="#C0392B", alpha=0.92))

    # --- Panel B (Track A #2): Pilot A-1 Modality-Decomposed Target LOEUF (mAb vs Small Molecule) ---
    ax = axes[0, 1]
    profs = a1["profiles"]
    groups = [
        ("SM\nBBW=1", [p["loeuf"] for p in profs if p["modality"] == "Small molecule" and p["black_box_warning"] == 1 and p["loeuf"] is not None], "#D35400"),
        ("SM\nBBW=0", [p["loeuf"] for p in profs if p["modality"] == "Small molecule" and p["black_box_warning"] == 0 and p["loeuf"] is not None], "#F5B041"),
        ("mAb\nBBW=1", [p["loeuf"] for p in profs if p["modality"] == "Antibody" and p["black_box_warning"] == 1 and p["loeuf"] is not None], "#1F618D"),
        ("mAb\nBBW=0", [p["loeuf"] for p in profs if p["modality"] == "Antibody" and p["black_box_warning"] == 0 and p["loeuf"] is not None], "#85C1E9"),
    ]
    for i, (lbl, vals, col) in enumerate(groups):
        bp = ax.boxplot([vals], positions=[i], widths=0.45, patch_artist=True,
                        boxprops=dict(facecolor=col, alpha=0.65), medianprops=dict(color="black", lw=2))
        jitter = np.linspace(-0.1, 0.1, len(vals))
        ax.scatter(np.full(len(vals), i) + jitter, vals, color=col, edgecolor="black", s=32, zorder=3)
    ax.set_xticks(range(4))
    ax.set_xticklabels([g[0] for g in groups])
    ax.set_ylabel("Primary Target gnomAD LOEUF\n(Lower = Stronger Constraint)")
    ax.set_title("A2 [Nat. Biotechnol.]: Modality-Constraint Asymmetry\n(LOEUF Predicts SM BBW [p=0.023], Null in mAbs [p=0.64])")
    ax.axvline(1.5, color="gray", ls="--", lw=1)
    ax.text(0.03, 0.86, "Small Molecules:\nMed 0.38 vs 1.06\np = 0.023", transform=ax.transAxes, fontsize=8,
            bbox=dict(boxstyle="round,pad=0.25", fc="#FEF5E7", ec="#D35400"))
    ax.text(0.55, 0.86, "Monoclonal Abs:\nMed 0.74 vs 0.71\np = 0.642 (Null)", transform=ax.transAxes, fontsize=8,
            bbox=dict(boxstyle="round,pad=0.25", fc="#EBF5FB", ec="#1F618D"))

    # --- Panel C (Track A #3): Pilot A-5 Trial-to-FAERS Geriatric Gap ---
    ax = axes[0, 2]
    rows_g = [r for r in a5["rows"] if r["trial_to_faers_65p_gap"] is not None]
    rows_g = sorted(rows_g, key=lambda r: r["trial_to_faers_65p_gap"])
    y = np.arange(len(rows_g))
    t_sh = [r["mean_trial_65plus_share"] * 100 for r in rows_g]
    f_sh = [r["faers_65plus_share"] * 100 for r in rows_g]
    for yi, ts, fs in zip(y, t_sh, f_sh):
        col = "#C0392B" if fs > ts else "#27AE60"
        ax.plot([ts, fs], [yi, yi], color=col, lw=2.2, alpha=0.8)
    ax.scatter(t_sh, y, color="#2E86C1", s=45, label="Phase III Trial Enrolled (>=65)", zorder=4)
    ax.scatter(f_sh, y, color="#C0392B", s=45, marker="D", label="Real-World FAERS (>=65)", zorder=4)
    ax.set_yticks(y)
    ax.set_yticklabels([r["drug"].title() for r in rows_g], fontsize=8.5)
    ax.set_xlabel("Share of Patients Aged >= 65 (%)")
    ax.set_xlim(0, 84)
    ax.set_ylim(-0.7, 9.7)
    ax.set_title("A3 [Nat. Med.]: Trial-to-Real-World Age Exclusion Gap\n(9/10 Drugs Under-Enroll >=65; Peak +33.3 pp)")
    ax.legend(frameon=False, fontsize=7.8, loc="upper left")

    # --- Panel D (Track B #1): Pilot B-2 Dual FAERS Distortion (Dilution + Lawyer Inflation) ---
    ax = axes[1, 0]
    labels = ["Tirzepatide\nPancreatitis\n(vs Sitagliptin)", "Semaglutide\nDeath Rate\n(MD vs Pooled)", "Tirzepatide\nDeath Rate\n(MD vs Pooled)"]
    md_vals = [1.665, 4.125, 2.781]
    pool_vals = [0.801, 1.799, 0.593]
    x = np.arange(len(labels))
    b_md = ax.bar(x - w/2, md_vals, w, label="Physician-Verified (qual=1)", color="#1A5276")
    b_pl = ax.bar(x + w/2, pool_vals, w, label="Pooled FAERS (Diluted by 75–91% Consumers)", color="#A9CCE3")
    for bars in (b_md, b_pl):
        for bar in bars:
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2, h + 0.09, f"{h:.2f}%", ha="center", va="bottom", fontsize=8.2)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_ylabel("Reporting Rate (%) in 2023–2026 FAERS")
    ax.set_ylim(0, 5.2)
    ax.set_title("B1 [Mktg Sci / ISR / MS]: Viral FAERS Distortion\n(2.1–4.7x Dilution + >90x Lawyer Litigation Inflation)")
    ax.legend(frameon=False, fontsize=7.8, loc="upper right")
    ax.text(0.02, 0.67, "Tirz vs Sita Pancreatitis:\n  MD: OR=3.06 (p=4.0e-6)\n  Pooled: OR=1.04 (p=0.78!)\nSema Lawyer IGE: 75.0%\n  vs MD <3.2% (OR>=90.5)",
            transform=ax.transAxes, fontsize=7.6, bbox=dict(boxstyle="round,pad=0.28", fc="#EBF5FB", ec="#1A5276", alpha=0.9))

    # --- Panel E (Track B #2): Pilot B-1 March 30, 2026 FDA Mass Delinquency Enforcement Sweep ---
    ax = axes[1, 1]
    ms = b1["monthly_series"]
    sel_months = ["2025-10", "2025-11", "2025-12", "2026-01", "2026-02", "2026-03", "2026-04", "2026-05", "2026-06", "2026-07", "2026-08"]
    def get_vals(arm):
        mp = {x["month"]: x["count"] for x in ms[arm]}
        return [mp.get(m, 0) for m in sel_months]
    v_ind = get_vals("treated_p23_us_industry")
    v_oth = get_vals("treated_p23_us_other")
    v_non = get_vals("ctrl_p23_nonus")
    xm = np.arange(len(sel_months))
    ax.plot(xm, v_ind, marker="o", lw=2.4, color="#B03A2E", label="US Industry P2/3 (FDAAA Treated)")
    ax.plot(xm, v_oth, marker="s", lw=2.0, ls="--", color="#2874A6", label="US Academic/Other P2/3 (Treated)")
    ax.plot(xm, v_non, marker="^", lw=2.0, ls=":", color="#1E8449", label="Non-US P2/3 (Exempt Control)")
    ax.axvline(5.0, color="black", ls="-.", lw=1.4)
    ax.set_ylim(0, 34)
    ax.text(4.85, 29.2, "Mar 30, 2026\nFDA Sweep\n(2,200+ sponsors)", fontsize=7.8, fontweight="bold", ha="right")
    ax.set_xticks(xm)
    ax.set_xticklabels([m[2:] for m in sel_months], rotation=40, fontsize=8)
    ax.set_ylabel("Overdue Trial Results Submitted / Mo")
    ax.set_title("B2 [MS / JFE / Res. Policy]: Mar 2026 FDA Sweep\n(Industry Surge OR=4.70 vs Non-US, p=8.6e-4)")
    ax.legend(frameon=False, fontsize=7.6, loc="upper right")

    # --- Panel F (Track B #3): Pilot B-3 openFDA Drug Shortages Duration & cGMP Recall Overlap ---
    ax = axes[1, 2]
    def parse_dt(s):
        for fmt in ("%m/%d/%Y", "%Y-%m-%d", "%Y%m%d"):
            try:
                return datetime.strptime(str(s).strip(), fmt)
            except ValueError:
                pass
        return None
    inj_y, oral_y = [], []
    for r in b3_raw:
        di, du = parse_dt(r.get("initial_posting_date") or ""), parse_dt(r.get("update_date") or "")
        if di and du and du >= di:
            yrs = (du - di).days / 365.25
            df = str(r.get("dosage_form") or "").lower()
            if any(k in df for k in ["inject", "solution", "infusion", "intravenous", "powder, for solution"]):
                inj_y.append(yrs)
            elif any(k in df for k in ["tablet", "capsule", "oral"]):
                oral_y.append(yrs)
    bins = np.linspace(0, 13, 27)
    ax.hist(inj_y, bins=bins, alpha=0.65, color="#8E44AD", density=True, label=f"Sterile Injectables (n={len(inj_y)}, med 6.5 yr)")
    ax.hist(oral_y, bins=bins, alpha=0.60, color="#16A085", density=True, label=f"Oral Solids (n={len(oral_y)}, med 3.3 yr)")
    ax.axvline(np.median(inj_y), color="#5B2C6F", lw=2, ls="--")
    ax.axvline(np.median(oral_y), color="#0E6251", lw=2, ls="--")
    ax.set_xlabel("Shortage Persistence Duration (Years)")
    ax.set_ylabel("Density")
    ax.set_title("B3 [MSOM / OR / MS]: Drug Shortage Persistence\n(Injectables 2x Longer p=8.6e-71; cGMP Overlap 3.84x)")
    ax.legend(frameon=False, fontsize=7.8, loc="upper right")
    ax.text(0.45, 0.56, "71.8% reasons UNSPECIFIED\n0.8% admit cGMP, yet 48.4%\nof shortage firms had\n2023–26 cGMP recalls\n(vs 5.1% labeling placebo)",
            transform=ax.transAxes, fontsize=7.6, bbox=dict(boxstyle="round,pad=0.28", fc="#F5EEF8", ec="#8E44AD", alpha=0.9))

    fig.savefig(FIG_PATH, dpi=220, bbox_inches="tight")
    shutil.copy2(FIG_PATH, ART_PATH)
    print(f"Saved figure to {FIG_PATH} ({FIG_PATH.stat().st_size} bytes) and copied to {ART_PATH}")


if __name__ == "__main__":
    main()
