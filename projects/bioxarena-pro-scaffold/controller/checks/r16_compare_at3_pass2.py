"""Controller: staged verification of A-T3 pass 2 (secondary diagnostics) against the main line's summary_r16.json.

A-T3 pass 2 writes round16/audit/results/secondary/{model}_{ckpt}__{line}_metrics.json. This script checks:
- G4 2x2 (X_ctrl, batch) decomposition, STATE_SE (A-T3 results.g4_factorial):
    decompositions[s]            vs main secondary.g4_factorial_decomposition[native_vs_<s> | b0_vs_b17_mis_<s>]
    symmetry_check_b17_vs_b0[s]  vs main ...[b17_vs_b0_sym_mis_<s>]  (Group 2 only)
  delta_batch / delta_ctrl / delta_both (r_orth_cal) and delta_batch / delta_ctrl (r_dorth_cal): estimate, CI bounds,
  valid_frac; interaction and batch share point values; stratum medians.
- HVG gene-label Perm100 null, STATE_HVG (A-T3 results.hvg_perm100) vs main secondary.STATE_HVG.perm100_r_orth_cal and
  perm100_r_dorth_cal: the 100 permutation scores in order, median, mean, q025, q975, max, empirical p, and the call
  "matched 97.5% lower bound > null q975".
Tolerance |diff| <= 1e-5 (Addendum A cross-party metric tolerance); categorical items identical.
Usage: python3 -I r16_compare_at3_pass2.py <summary_r16.json> <at3_secondary_dir> [lines...]
"""
import json
import sys
from pathlib import Path

import numpy as np

TOL = 1e-5
S = json.load(open(sys.argv[1]))["lines"]
D = Path(sys.argv[2])
lines = sys.argv[3:] or ["nadig_hepg2", "nadig_jurkat", "a549_ko", "dld1", "cd4t_rest"]
BOUNDS = ["estimate", "ci95_lo", "ci95_hi", "ci975_lo", "ci975_hi", "valid_frac"]
worst, bad, n, notes = 0.0, [], 0, []


def num(tag, a, b):
    global worst, n
    n += 1
    a, b = float(a), float(b)
    if np.isnan(a) and np.isnan(b):
        return
    d = abs(a - b); worst = max(worst, d)
    if not d <= TOL:
        bad.append(f"{tag}: main={a} at3={b} diff={d:.3g}")


def eq(tag, a, b):
    global n
    n += 1
    if a != b:
        bad.append(f"{tag}: main={a!r} at3={b!r}")


def summ(tag, m, a):
    for k in BOUNDS:
        num(f"{tag}.{k}", m[k], a[k])


def load(model, ln):
    f = sorted(D.glob(f"{model}_*__{ln}_metrics.json"))
    if not f:
        bad.append(f"{ln}: no A-T3 {model} metrics file in {D}")
        return None
    return json.load(open(f[0]))["results"]


def g4_pair(tag, m, a):
    mo, md = m["r_orth_cal"], m["r_dorth_cal"]
    summ(f"{tag}.orth.delta_batch", mo["delta_batch"], a["delta_batch_orth"])
    summ(f"{tag}.orth.delta_ctrl", mo["delta_ctrl"], a["delta_ctrl_orth"])
    summ(f"{tag}.orth.delta_both", mo["delta_both"], a["delta_both_orth"])
    num(f"{tag}.orth.interaction", mo["delta_interaction_point"], a["delta_interaction_orth"])
    num(f"{tag}.orth.batch_share", mo["batch_share_point"], a["batch_share_orth"])
    summ(f"{tag}.dorth.delta_batch", md["delta_batch"], a["delta_batch_dorth"])
    summ(f"{tag}.dorth.delta_ctrl", md["delta_ctrl"], a["delta_ctrl_dorth"])


for ln in lines:
    sec = S[ln]["secondary"]
    # --- G4 (STATE_SE)
    a = load("STATE_SE", ln)
    if a is not None:
        g = a.get("g4_factorial")
        mg = sec["g4_factorial_decomposition"]
        if g is None:
            bad.append(f"{ln}: A-T3 STATE_SE has no g4_factorial")
        else:
            native = ln in ("nadig_hepg2", "nadig_jurkat")
            want = sorted(k.split("native_vs_", 1)[1] for k in mg if k.startswith("native_vs_")) if native else \
                sorted(k.split("b0_vs_b17_mis_", 1)[1] for k in mg if k.startswith("b0_vs_b17_mis_"))
            eq(f"{ln}.g4.sources", want, sorted(g["decompositions"]))
            for s in want:
                if s not in g["decompositions"]:
                    continue
                g4_pair(f"{ln}.g4.{s}", mg[f"native_vs_{s}" if native else f"b0_vs_b17_mis_{s}"], g["decompositions"][s])
                if not native:
                    sym = (g.get("symmetry_check_b17_vs_b0") or {}).get(s)
                    if sym is None:
                        bad.append(f"{ln}.g4.sym.{s}: missing in A-T3")
                    else:
                        g4_pair(f"{ln}.g4.sym.{s}", mg[f"b17_vs_b0_sym_mis_{s}"], sym)
            if not native:  # single-index spread b in {0, 1, 17, 35} (PREREG_I §4.4.3), r_orth_cal
                sp_m = sec["g4_single_index_spread"].get("per_index_r_orth_cal", {})
                sp_a = g.get("batch_spread") or {}
                eq(f"{ln}.g4.spread.indices", sorted(k[1:] for k in sp_m), sorted(sp_a))
                for b in sp_a:
                    if f"b{b}" in sp_m:
                        summ(f"{ln}.g4.spread.b{b}", sp_m[f"b{b}"], sp_a[b])
            ms_all = sec["g4_stratum_summary"]
            ms = ms_all.get("native" if native else "b0_vs_b17", {})
            for st in ("seen_by_checkpoint", "unseen_by_checkpoint"):
                at, mt = g["stratum_summaries"].get(st, {}), ms.get(st, {})
                eq(f"{ln}.g4.{st}.n_sources", mt.get("n_sources", 0), at.get("n_sources", 0))
                if at.get("n_sources", 0) and mt.get("n_sources", 0):
                    num(f"{ln}.g4.{st}.median_delta_batch", mt["r_orth_cal"]["delta_batch"]["median"], at["median_delta_batch_orth"])
                    num(f"{ln}.g4.{st}.median_delta_ctrl", mt["r_orth_cal"]["delta_ctrl"]["median"], at["median_delta_ctrl_orth"])
                    if "batch_share_point" in mt["r_orth_cal"]:
                        num(f"{ln}.g4.{st}.median_batch_share", mt["r_orth_cal"]["batch_share_point"]["median"], at["median_batch_share_orth"])
    # --- HVG Perm100 (STATE_HVG)
    a = load("STATE_HVG", ln)
    if a is not None:
        p = a.get("hvg_perm100")
        mp, mpd = sec["STATE_HVG"]["perm100_r_orth_cal"], sec["STATE_HVG"]["perm100_r_dorth_cal"]
        if p is None:
            bad.append(f"{ln}: A-T3 STATE_HVG has no hvg_perm100")
        else:
            pa, pm = np.asarray(p["perm_scores"], float), np.asarray(mp["perm_estimates"], float)
            eq(f"{ln}.perm.n_perm", int(mp["n_perm"]), int(p["n_perm"]))
            eq(f"{ln}.perm.n_finite", int(mp["n_finite_perms"]), int(p["n_finite"]))
            if pa.shape == pm.shape:
                dd = float(np.nanmax(np.abs(pa - pm)))
                num(f"{ln}.perm.scores_maxdiff", 0.0, dd)
                if dd > TOL:
                    notes.append(f"{ln}: perm scores differ element-wise (max {dd:.3g}); sorted max diff "
                                 f"{float(np.nanmax(np.abs(np.sort(pa) - np.sort(pm)))):.3g} (seed-rule check)")
            for km, ka in [("median", "median"), ("mean", "mean"), ("q025", "q025"), ("q975", "q975"), ("max", "max"),
                           ("empirical_p_ge", "p_emp"), ("observed_estimate", "matched_estimate")]:
                num(f"{ln}.perm.{km}", mp[km], p[ka])
            m_lo = sec["STATE_HVG"]["r_orth_cal_matched"]["ci975_lo"]
            eq(f"{ln}.perm.matched_gt_q975", bool(m_lo > mp["q975"]), bool(p["matched_gt_q975"]))
            ds = p.get("dorth_summary")
            if ds is None:
                bad.append(f"{ln}.perm.dorth: missing in A-T3")
            else:
                for k in ["median", "mean", "q025", "q975"]:
                    num(f"{ln}.perm.dorth.{k}", mpd[k], ds[k])

print(json.dumps({"lines": lines, "n_checks": n, "n_fail": len(bad), "worst_abs_diff": worst}))
for x in notes:
    print("NOTE", x)
for b in bad:
    print("FAIL", b)
print("R16_AT3_PASS2_MATCH" if not bad else "R16_AT3_PASS2_DIFF")
