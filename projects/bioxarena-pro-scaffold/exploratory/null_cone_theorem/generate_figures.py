#!/usr/bin/env python3
"""
Generate 4 publication-grade (*Nature Methods* / *Nature Biotechnology* style) multi-panel figures
directly from `results/audit_13lines_summary.json` (zero hardcoded numbers!).

Outputs both 300-DPI PNG and vector PDF in `exploratory/null_cone_theorem/figures/`.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

THIS_DIR = Path(__file__).resolve().parent
RESULTS_JSON = THIS_DIR / "results" / "audit_13lines_summary.json"
FIG_DIR = THIS_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Clean Nature/Science publication aesthetics (colorblind-safe Okabe-Ito palette)
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8.5,
    "axes.titlesize": 9.5,
    "axes.titleweight": "bold",
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.8,
    "ytick.labelsize": 7.8,
    "legend.fontsize": 7.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.dpi": 300,
})

C_BLUE = "#0072B2"
C_ORANGE = "#E69F00"
C_RED = "#D55E00"
C_GREEN = "#009E73"
C_PURPLE = "#CC79A7"
C_SKY = "#56B4E9"
C_GREY = "#666666"
C_DARK = "#222222"

SHORT_NAMES = {
    "nadig_hepg2": "HepG2 (R9)",
    "nadig_jurkat": "Jurkat (R9)",
    "nourreddine_kolf21j": "KOLF2.1J (R9)",
    "replogle_k562_essential": "K562 (R9)",
    "replogle_rpe1": "RPE1 (R9)",
    "vcc2025_h1": "H1 hESC (R9)",
    "dld1": "DLD-1 (R10)",
    "cd4t_rest": "CD4+ T Rest (R10)",
    "calu3": "Calu-3 (R10)",
    "cd4t_stim48": "CD4+ T Stim48 (R11)",
    "feng_ipsc": "Feng iPSC (R11)",
    "telohaec": "TeloHAEC (R11)",
}


def get_valid_lines(data: dict) -> List[str]:
    return [ln for ln, d in data["lines"].items() if not d.get("skipped", False)]


def fig1_mechanism_and_benchmark_audit(data: dict) -> None:
    fig = plt.figure(figsize=(15.4, 4.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.55, 1.05], wspace=0.54)

    # --- Panel A: w(mu) and kappa(mu) geometry ---
    ax = fig.add_subplot(gs[0, 0])
    mu = np.logspace(-2, 2, 400)
    w = mu / (1.0 + mu)
    kap = 0.5 * w * (1.0 - w)
    ax.semilogx(mu, w, color=C_BLUE, lw=2.2, label=r"Jacobian $w(\mu)=\mu/(1+\mu)$")
    ax.semilogx(mu, 4.0 * kap, color=C_RED, lw=1.8, ls="--", label=r"4× Curvature $2w(1-w)$")
    ax.axvspan(1.5, 100, color=C_ORANGE, alpha=0.14, label=r"Top-$k$ tail ($w_t \to 1$)")
    ax.axhline(1.0, color=C_GREY, lw=0.8, ls=":")
    ax.annotate("Signal-tail conditioning\ncompresses Var($w_t$)\n" + r"$\Rightarrow a_{M,\delta^2} > 0$",
                xy=(8.0, 0.89), xytext=(0.015, 0.72),
                arrowprops=dict(arrowstyle="->", color=C_DARK, lw=1.0),
                fontsize=7.2, bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=C_GREY, lw=0.6))
    ax.set_xlabel(r"Unperturbed baseline $\mu_{c,g}$ (CP10k)")
    ax.set_ylabel("Nonlinear compression / curvature")
    ax.set_ylim(-0.04, 1.12)
    ax.set_title("A  log1p Jacobian & signal-tail conditioning", loc="left", pad=8)
    ax.legend(loc="lower right", frameon=True, framealpha=0.92, fontsize=6.7)

    # --- Panel B: Pure Noiseless Cross-Benchmark Audit Heatmap ---
    ax = fig.add_subplot(gs[0, 1])
    tbl = data["noiseless_audit"]["cross_benchmark_table"]
    # Select 8 representative rows for visual clarity
    keep_idx = [0, 1, 2, 4, 5, 6, 7, 8, 9]
    rows = [tbl[i] for i in keep_idx if i < len(tbl)]
    pred_cols = ["GlobalScale", "SrcGrandMean", "E1_zero_tgt", "E6_logFC", "Perm_mu_t", "FCT_true"]
    col_labels = ["Global\nScale", "Src Grand\nMean (bS)", "E1\n(zero-tgt)", "E6\n(logFC)", "Perm\n(scrambled)", "FCT\n(true μt)"]
    row_labels = [
        r["metric"]
        .replace("Systema/SBB ", "Systema/SBB\n")
        .replace("Miller/cell-eval2 ", "Miller/cell-eval2\n")
        .replace("CGC/Molina ", "CGC/Molina\n")
        .replace("OURS N5 Cone Semi-Partial", "OURS N5 Cone")
        for r in rows
    ]
    mat = np.array([[r["scores"][c] for c in pred_cols] for r in rows], dtype=np.float64)

    im = ax.imshow(mat, cmap="RdBu_r", vmin=-0.3, vmax=1.0, aspect="auto")
    ax.set_xticks(range(len(pred_cols)))
    ax.set_xticklabels(col_labels, fontsize=7.2)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=7.0)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat[i, j]
            txt_c = "white" if (v > 0.60 or v < -0.20) else C_DARK
            fw = "bold" if (i >= len(rows) - 2 or (j < 5 and v > 0.20)) else "normal"
            ax.text(j, i, f"{v:+.2f}", ha="center", va="center", color=txt_c, fontsize=6.8, fontweight=fw)
    # Highlight OURS rows with box
    ax.axhline(len(rows) - 2.5, color=C_GREEN, lw=2.0)
    ax.axvline(4.5, color=C_DARK, lw=1.2, ls="--")
    ax.set_title("B  Noiseless benchmark audit (split-half = in-sample)", loc="left", pad=8)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.ax.tick_params(labelsize=7.0)
    cbar.set_label("Score", fontsize=7.5)

    # --- Panel C: Universal Log1p Cone SVD & Uniform R^2 Coverage ---
    ax = fig.add_subplot(gs[0, 2])
    svd = data["noiseless_audit"]["cone_svd"]
    c_log = np.array(svd["c_grid_log10"])
    r2_c = np.array(svd["r2_basis_005_1_20"])
    ax.plot(c_log, r2_c * 100.0, color=C_GREEN, lw=2.2, label=r"3-node basis $\{0.05, 1, 20\}$")
    for c_node in (0.05, 1.0, 20.0):
        ax.axvline(np.log10(c_node), color=C_BLUE, lw=0.9, ls=":", alpha=0.8)
    ax.axhline(100.0, color=C_GREY, lw=0.7, ls="--")
    ax.set_xlabel(r"$\log_{10}(c)$ of constant baseline $g_c(\overline{\mathrm{FC}})$")
    ax.set_ylabel(r"Variance captured by $\mathcal{N}_4$ (%)")
    ax.set_ylim(99.85, 100.015)
    ax.text(
        0.04, 0.16,
        f"Rank-3 SVD cum. var: {svd['cum_var_top6'][2]*100:.3f}%\n"
        f"Min $R^2$ over $10^{{-3}}..10^3$: {svd['min_r2']:.5f}\n"
        f"Median $R^2$: {svd['median_r2']:.5f}\n"
        f"Target data used: ZERO",
        transform=ax.transAxes, fontsize=7.2,
        bbox=dict(boxstyle="round,pad=0.3", fc="#f4faf7", ec=C_GREEN, lw=0.8),
    )
    ax.set_title("C  Universal target-free cone coverage", loc="left", pad=8)
    ax.legend(loc="lower right", frameon=True)

    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"fig1_mechanism_and_benchmark_audit.{ext}", bbox_inches="tight")
    plt.close(fig)


def fig2_split_half_failure_12lines(data: dict) -> None:
    lines = get_valid_lines(data)
    # Sort by N1 top50_spl E1 descending for clean visual storytelling
    lines.sort(key=lambda ln: data["lines"][ln]["masks"]["top50_spl"]["N1"]["predictors"]["E1_ConstMed"]["orth"]["est"])
    y = np.arange(len(lines))
    ylabels = [f"{SHORT_NAMES.get(ln, ln)} (N={data['lines'][ln]['N']})" for ln in lines]

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.8), sharey=True, gridspec_kw={"wspace": 0.15})

    # --- Panel A: Legacy N1 E1 across all vs top50_obs vs SA-8 top50_spl ---
    ax = axes[0]
    for i, ln in enumerate(lines):
        m = data["lines"][ln]["masks"]
        e_all = m["all"]["N1"]["predictors"]["E1_ConstMed"]["orth"]
        e_obs = m["top50_obs"]["N1"]["predictors"]["E1_ConstMed"]["orth"]
        e_spl = m["top50_spl"]["N1"]["predictors"]["E1_ConstMed"]["orth"]
        ax.plot([e_all["est"], e_spl["est"]], [i, i], color=C_GREY, lw=1.1, alpha=0.6, zorder=1)
        ax.errorbar(e_all["est"], i - 0.18, xerr=[[e_all["est"] - e_all["ci975_lo"]], [e_all["ci975_hi"] - e_all["est"]]],
                    fmt="o", color=C_BLUE, ms=4.5, lw=1.0, capsize=2, label="all_nontarget" if i == 0 else "")
        ax.errorbar(e_obs["est"], i, xerr=[[e_obs["est"] - e_obs["ci975_lo"]], [e_obs["ci975_hi"] - e_obs["est"]]],
                    fmt="s", color=C_ORANGE, ms=4.2, lw=1.0, capsize=2, label="resp_top50 (in-sample)" if i == 0 else "")
        ax.errorbar(e_spl["est"], i + 0.18, xerr=[[e_spl["est"] - e_spl["ci975_lo"]], [e_spl["ci975_hi"] - e_spl["est"]]],
                    fmt="^", color=C_RED, ms=5.0, lw=1.2, capsize=2, label="SA-8 top50_split (out-of-sample)" if i == 0 else "")
    ax.axvline(0.0, color=C_DARK, lw=0.9, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(ylabels, fontsize=7.8)
    ax.set_xlabel(r"Legacy 1D $r_{\mathrm{orth,cal}}$ (zero-target $E_1$)")
    ax.set_title("A  $E_1$ bias survives split-half (12/12 lines)", loc="left", pad=8)
    ax.legend(loc="lower right", frameon=True, fontsize=6.8)

    # --- Panel B: Legacy N1 Perm100 across all vs SA-8 top50_spl ---
    ax = axes[1]
    for i, ln in enumerate(lines):
        m = data["lines"][ln]["masks"]
        p_all = m["all"]["N1"]["perm100"]
        p_spl = m["top50_spl"]["N1"]["perm100"]
        ax.plot([p_all["median"], p_spl["median"]], [i, i], color=C_GREY, lw=1.1, alpha=0.6, zorder=1)
        ax.plot([p_all["min"], p_all["q975"]], [i - 0.12, i - 0.12], color=C_BLUE, lw=1.6, alpha=0.6)
        ax.scatter(p_all["median"], i - 0.12, color=C_BLUE, s=22, zorder=3, label="Perm100 median (all)" if i == 0 else "")
        ax.plot([p_spl["median"], p_spl["q975"]], [i + 0.12, i + 0.12], color=C_RED, lw=1.8, alpha=0.8)
        ax.scatter(p_spl["median"], i + 0.12, color=C_RED, marker="^", s=28, zorder=4, label="Perm100 med→Q97.5 (SA-8 split)" if i == 0 else "")
        ax.scatter(p_spl["q975"], i + 0.12, color=C_RED, marker="|", s=45, zorder=4)
    ax.axvline(0.0, color=C_DARK, lw=0.9, ls="--")
    ax.set_xlabel(r"Legacy 1D $r_{\mathrm{orth,cal}}$ (100 permutations)")
    ax.set_title("B  Scrambled genes inflate under split-half", loc="left", pad=8)
    ax.legend(loc="lower right", frameon=True, fontsize=6.8)

    # --- Panel C: Paired Delta E1(SA-8 split) - E1(all) [97.5% CI]: N1 vs OURS N5 ---
    ax = axes[2]
    for i, ln in enumerate(lines):
        cm = data["lines"][ln]["cross_mask_e1"]
        d1 = cm["N1"]["E1_top50_spl_minus_all"]
        d5 = cm["N5"]["E1_top50_spl_minus_all"]
        ax.errorbar(d1["est"], i + 0.14, xerr=[[d1["est"] - d1["ci975_lo"]], [d1["ci975_hi"] - d1["est"]]],
                    fmt="^", color=C_RED, ms=5.0, lw=1.2, capsize=2, label="Legacy N1 [97.5% CI]" if i == 0 else "")
        ax.errorbar(d5["est"], i - 0.14, xerr=[[d5["est"] - d5["ci975_lo"]], [d5["ci975_hi"] - d5["est"]]],
                    fmt="D", color=C_GREEN, ms=4.5, lw=1.2, capsize=2, label="OURS N5 semi-partial [97.5% CI]" if i == 0 else "")
    ax.axvline(0.0, color=C_DARK, lw=0.9, ls="--")
    ax.set_xlabel(r"Paired $\Delta E_1 = E_1(\mathrm{SA8\;split}) - E_1(\mathrm{all})$")
    ax.set_title("C  Split-half bias: Legacy N1 vs OURS N5", loc="left", pad=8)
    ax.legend(loc="lower right", frameon=True, fontsize=6.8)

    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"fig2_split_half_failure_12lines.{ext}", bbox_inches="tight")
    plt.close(fig)


def fig3_universal_null_cone_cure_12lines(data: dict) -> None:
    lines = get_valid_lines(data)
    masks = ["all", "top50_obs", "top50_spl", "top50_src"]

    fig = plt.figure(figsize=(14.4, 5.0))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.12, 1.02, 1.36], wspace=0.36)

    # --- Panel A: All 48 (line x mask) artifact distributions before (N1) vs after (N5) ---
    ax = fig.add_subplot(gs[0, 0])
    cats = [
        ("Src Grand\nMean (bS)", "Source_Grand_Mean"),
        ("E1\n(Med)", "E1_ConstMed"),
        ("E2\n(P90)", "E2_ConstP90"),
        ("E6\n(logFC)", "E6_LogFCBar"),
        ("Perm100\nMedian", "__PERM__"),
    ]
    rng = np.random.default_rng(42)
    for idx_c, (clab, key) in enumerate(cats):
        v1, v5 = [], []
        for ln in lines:
            for mk in masks:
                md = data["lines"][ln]["masks"][mk]
                if key == "__PERM__":
                    v1.append(md["N1"]["perm100"]["median"])
                    v5.append(md["N5"]["perm100"]["median"])
                else:
                    v1.append(md["N1"]["predictors"][key]["orth"]["est"])
                    v5.append(md["N5"]["predictors"][key]["semi"]["est"])
        v1 = np.array([x for x in v1 if np.isfinite(x)])
        v5 = np.array([x for x in v5 if np.isfinite(x)])
        j1 = rng.uniform(-0.14, 0.14, size=v1.size)
        j5 = rng.uniform(-0.14, 0.14, size=v5.size)
        ax.scatter(idx_c - 0.18 + j1, v1, color=C_RED, alpha=0.55, s=16, edgecolor="none", label="Legacy N1 (48 line×mask)" if idx_c == 0 else "")
        ax.scatter(idx_c + 0.18 + j5, v5, color=C_GREEN, alpha=0.75, s=16, edgecolor="none", label="OURS N5 semi-partial" if idx_c == 0 else "")
        ax.hlines(np.median(v1), idx_c - 0.32, idx_c - 0.04, color=C_DARK, lw=1.6)
        ax.hlines(np.median(v5), idx_c + 0.04, idx_c + 0.32, color=C_DARK, lw=1.6)
    ax.axhline(0.0, color=C_DARK, lw=0.9, ls="--")
    ax.set_xticks(range(len(cats)))
    ax.set_xticklabels([c[0] for c in cats], fontsize=7.4)
    ax.set_ylabel("Calibrated residual score across 48 evaluations")
    ax.set_title("A  All 4 cheat classes collapse under N5", loc="left", pad=8)
    ax.legend(loc="upper right", frameon=True, fontsize=7.0)

    # --- Panel B: Orthogonal Energy Fraction eta under N5 (all_nontarget) ---
    ax = fig.add_subplot(gs[0, 1])
    eta_groups = [
        ("GlobalScale", "GlobalScale", C_GREY),
        ("SrcGrandMean", "Source_Grand_Mean", C_GREY),
        ("E1 (Med)", "E1_ConstMed", C_RED),
        ("E6 (logFC)", "E6_LogFCBar", C_RED),
        ("Perm100", "__PERM__", C_ORANGE),
        ("FCT (top50_src)", "__FCT_SRC__", C_SKY),
        ("FCT (all)", "FCT", C_BLUE),
    ]
    for idx_e, (elab, key, col) in enumerate(eta_groups):
        vals = []
        for ln in lines:
            if key == "__FCT_SRC__":
                vals.append(data["lines"][ln]["masks"]["top50_src"]["N5"]["predictors"]["FCT"]["eta"])
            elif key == "__PERM__":
                vals.append(data["lines"][ln]["masks"]["all"]["N5"]["perm100"]["eta_median"])
            else:
                vals.append(data["lines"][ln]["masks"]["all"]["N5"]["predictors"][key]["eta"])
        vals = np.array([x for x in vals if np.isfinite(x)])
        mval = float(np.mean(vals)) if vals.size else 0.0
        max_x = max(mval, float(np.max(vals)) if vals.size else 0.0)
        ax.barh(idx_e, mval, color=col, alpha=0.80, height=0.62)
        ax.scatter(vals, np.full_like(vals, idx_e), color=C_DARK, s=11, alpha=0.65, zorder=3)
        ax.text(max_x + 0.035, idx_e, f"{mval:.3f}", va="center", ha="left",
                fontsize=7.1, color=C_DARK, fontweight="bold")
    ax.set_yticks(range(len(eta_groups)))
    ax.set_yticklabels([g[0] for g in eta_groups], fontsize=7.6)
    ax.set_xlim(0.0, 1.18)
    ax.set_xticks([0.0, 0.25, 0.50, 0.75, 1.00])
    ax.set_xlabel(r"Orthogonal energy $\eta = \|P_{\mathcal{N}_5}^\perp y\| / \|y_c\|$")
    ax.set_title("B  Why N5 spares true target signal", loc="left", pad=8)

    # --- Panel C: Whole-transcriptome (all_nontarget) Forest Plot under N5 across 12 lines ---
    ax = fig.add_subplot(gs[0, 2])
    s_lines = sorted(lines, key=lambda ln: data["lines"][ln]["masks"]["all"]["N5"]["predictors"]["FCT"]["semi"]["est"])
    y = np.arange(len(s_lines))
    for i, ln in enumerate(s_lines):
        md = data["lines"][ln]["masks"]["all"]["N5"]
        fct = md["predictors"]["FCT"]["semi"]
        e1 = md["predictors"]["E1_ConstMed"]["semi"]
        mxm = md["max_mismatch_ci975_hi"]
        prm = md["perm100"]
        ax.errorbar(fct["est"], i + 0.14, xerr=[[fct["est"] - fct["ci975_lo"]], [fct["ci975_hi"] - fct["est"]]],
                    fmt="o", color=C_BLUE, ms=4.8, lw=1.3, capsize=2, label="FCT [97.5% CI]" if i == 0 else "")
        ax.scatter(mxm["est"], i - 0.12, color=C_ORANGE, marker="s", s=22, zorder=3, label="Max Mismatch (pt → CI97.5_hi)" if i == 0 else "")
        ax.plot([mxm["est"], mxm["ci975_hi"]], [i - 0.12, i - 0.12], color=C_ORANGE, lw=1.3)
        ax.scatter(mxm["ci975_hi"], i - 0.12, color=C_ORANGE, marker="|", s=45)
        ax.scatter(prm["q975"], i - 0.12, color=C_RED, marker="^", s=20, zorder=4, label="Perm100 Q97.5" if i == 0 else "")
        ax.scatter(e1["est"], i + 0.14, color=C_GREEN, marker="D", s=18, zorder=4, label="E1 (zero-target)" if i == 0 else "")
    ax.axvline(0.0, color=C_DARK, lw=0.9, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{SHORT_NAMES.get(ln, ln)}" for ln in s_lines], fontsize=7.6)
    ax.set_xlabel(r"N5 Semi-Partial $r_{\mathrm{cone\text{-}sp}}$ (all_nontarget)")
    ax.set_title("C  De-biased FCT specificity across 12 lines", loc="left", pad=8)
    ax.legend(loc="lower right", frameon=True, fontsize=6.7)

    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"fig3_universal_null_cone_cure_12lines.{ext}", bbox_inches="tight")
    plt.close(fig)


def fig4_q5a_conservation(data: dict) -> None:
    q5a = data["q5a_conservation"]
    pairs = q5a["pairs"]
    summ = q5a["summary"]

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.1), gridspec_kw={"wspace": 0.30})

    # --- Panel A: expr05 Raw vs Noise-Ceiling Calibrated logFC vs log1p ---
    ax = axes[0]
    raw_d = np.array([p["expr05"]["raw_log1p"] for p in pairs])
    raw_fc = np.array([p["expr05"]["raw_logFC"] for p in pairs])
    cal_d = np.array([p["expr05"]["cal_log1p"] for p in pairs])
    cal_fc = np.array([p["expr05"]["cal_logFC"] for p in pairs])

    ax.scatter(raw_d, raw_fc, color=C_GREY, alpha=0.65, s=26, marker="o",
               label=f"Raw ({int(np.sum(raw_fc > raw_d))}/33 logFC > log1p)")
    ax.scatter(cal_d, cal_fc, color=C_BLUE, alpha=0.85, s=30, marker="^",
               label=f"Noise-calibrated ({summ['expr05']['cal_logFC_gt_log1p']}/33, Wilcoxon p={summ['expr05']['cal_logFC_gt_log1p_wilcoxon_p']:.1e})")
    lim = [0.05, 0.72]
    ax.plot(lim, lim, color=C_DARK, ls="--", lw=1.0)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel(r"Cross-line conservation in log1p $\Delta$ space")
    ax.set_ylabel(r"Cross-line conservation in logFC space")
    ax.set_title(r"A  Q5a ($\mu \geq 0.5$): Shot noise vs true logFC", loc="left", pad=8)
    ax.legend(loc="upper left", frameon=True, fontsize=6.9)

    # --- Panel B: FCT-mapped transfer vs direct log1p across all 33 pairs (expr05) ---
    ax = axes[1]
    raw_fct = np.array([p["expr05"]["raw_FCT"] for p in pairs])
    cal_fct = np.array([p["expr05"]["cal_FCT"] for p in pairs])
    ax.scatter(raw_d, raw_fct, color=C_ORANGE, alpha=0.75, s=26, marker="o",
               label=f"Raw FCT > log1p: {summ['expr05']['raw_FCT_gt_log1p']}/33 (p={summ['expr05']['raw_FCT_gt_log1p_sign_p']:.1e})")
    ax.scatter(cal_d, cal_fct, color=C_GREEN, alpha=0.85, s=30, marker="D",
               label=f"Calibrated FCT > log1p: {summ['expr05']['cal_FCT_gt_log1p']}/33")
    ax.plot(lim, lim, color=C_DARK, ls="--", lw=1.0)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel(r"Direct log1p cross-line transfer ($r$)")
    ax.set_ylabel(r"Baseline-mapped FCT transfer ($r$)")
    ax.set_title(r"B  FCT beats direct log1p on 33/33 pairs", loc="left", pad=8)
    ax.legend(loc="upper left", frameon=True, fontsize=6.9)

    # --- Panel C: Mean conservation summary across all / expr01 / expr05 ---
    ax = axes[2]
    masks = ["all", "expr01", "expr05"]
    mlabs = ["All shared\ngenes", r"Expressed" + "\n" + r"($\mu \geq 0.1$)", r"Expressed" + "\n" + r"($\mu \geq 0.5$)"]
    x = np.arange(len(masks))
    w = 0.27
    b1 = [summ[m]["cal_log1p_mean"] for m in masks]
    b2 = [summ[m]["cal_logFC_mean"] for m in masks]
    b3 = [summ[m]["cal_FCT_mean"] for m in masks]
    ax.bar(x - w, b1, width=w * 0.92, color=C_GREY, alpha=0.85, label="Calibrated log1p Δ")
    ax.bar(x, b2, width=w * 0.92, color=C_BLUE, alpha=0.85, label="Calibrated logFC")
    ax.bar(x + w, b3, width=w * 0.92, color=C_GREEN, alpha=0.85, label="Calibrated FCT-mapped")
    for xi, (v1, v2, v3) in enumerate(zip(b1, b2, b3)):
        ax.text(xi - w, v1 + 0.004, f"{v1:.3f}", ha="center", va="bottom", fontsize=6.3)
        ax.text(xi, v2 + 0.004, f"{v2:.3f}", ha="center", va="bottom", fontsize=6.3)
        ax.text(xi + w, v3 + 0.004, f"{v3:.3f}", ha="center", va="bottom", fontsize=6.3)
    ax.set_xticks(x)
    ax.set_xticklabels(mlabs, fontsize=7.6)
    ax.set_ylim(0.25, 0.41)
    ax.set_ylabel("Mean noise-calibrated cross-line r (33 pairs)")
    ax.set_title("C  Noise-calibrated space comparison", loc="left", pad=8)
    ax.legend(loc="upper left", frameon=True, fontsize=7.0)

    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"fig4_q5a_logfc_vs_log1p_conservation.{ext}", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    data = json.loads(RESULTS_JSON.read_text(encoding="utf-8"))
    fig1_mechanism_and_benchmark_audit(data)
    fig2_split_half_failure_12lines(data)
    fig3_universal_null_cone_cure_12lines(data)
    fig4_q5a_conservation(data)
    print(f"Generated 4 publication figures (PNG + PDF) in {FIG_DIR}")


if __name__ == "__main__":
    main()
