"""Controller: compare r16_check.py recomputation with the main line's round16/results/summary_r16.json.

Addendum A tolerance for cross-party metrics: |diff| <= 1e-5 on every estimate and CI bound; categorical
calls (per-source categories, H-I1 line verdict, H-I2 call) must be identical.
Usage: python3 -I r16_compare_main.py <summary_r16.json> <controller_json_dir> [lines...]
"""
import json
import sys
from pathlib import Path

TOL = 1e-5
S = json.load(open(sys.argv[1]))["lines"]
C = Path(sys.argv[2])
lines = sys.argv[3:] or ["nadig_hepg2", "nadig_jurkat", "a549_ko", "dld1", "cd4t_rest"]
worst, bad = 0.0, []


def num(tag, a, b):
    global worst
    d = abs(float(a) - float(b)); worst = max(worst, d)
    if not d <= TOL:
        bad.append(f"{tag}: main={a} controller={b} diff={d:.3g}")


def summ(tag, m, c):  # main summary dict vs controller {point, ci975, ci95, valid_frac}
    num(tag + ".point", m["estimate"], c["point"])
    num(tag + ".ci975_lo", m["ci975_lo"], c["ci975"][0]); num(tag + ".ci975_hi", m["ci975_hi"], c["ci975"][1])
    num(tag + ".ci95_lo", m["ci95_lo"], c["ci95"][0]); num(tag + ".ci95_hi", m["ci95_hi"], c["ci95"][1])
    num(tag + ".valid_frac", m["valid_frac"], c["valid_frac"])


def same(tag, a, b):
    if a != b:
        bad.append(f"{tag}: main={a!r} controller={b!r}")


for ln in lines:
    p = S[ln]["primary"]; c = json.load(open(C / f"{ln}.json"))
    same(f"{ln}.N", S[ln]["n_eval_targets"], c["N"]); same(f"{ln}.G", S[ln]["n_eval_genes"], c["G"])
    num(f"{ln}.rho_ab_orth", p["STATE_SE_matched_orth"]["rho_ab_orth"], c["rho_ab_orth"])
    summ(f"{ln}.matched", p["STATE_SE_matched_orth"], c["arms"]["matched"])
    summ(f"{ln}.FCT", p["FCT_matched_orth"], c["arms"]["FCT"])
    ps = p["STATE_SE_paired_vs_mismatched"]["per_source"]
    same(f"{ln}.sources", sorted(ps), sorted(c["pairs"]))
    strata = p["H_I1_adjudication"]["strata"]
    cat = {s: v["category"] for st in strata.values() for s, v in st["per_source"].items()}
    for s in c["pairs"]:
        summ(f"{ln}.pair[{s}]", ps[s], c["pairs"][s])
        num(f"{ln}.mis[{s}]", ps[s]["mismatched_estimate"], c["arms"]["mis_" + s]["point"])
        same(f"{ln}.pair[{s}].category", cat.get(s), c["pairs"][s]["call"])
    hi1 = p["H_I1_adjudication"]
    same(f"{ln}.H_I1", hi1["line_summary"], c["H_I1_line"])
    summ(f"{ln}.SE_minus_FCT", p["STATE_SE_vs_FCT"], c["H_I2"])
    same(f"{ln}.H_I2", p["H_I2_adjudication"]["classification"], c["H_I2"]["call"])

print(f"lines={lines} worst_abs_diff={worst:.3g} tol={TOL}")
for b in bad:
    print("MISMATCH", b)
print("R16_MAIN_VS_CONTROLLER_MATCH" if not bad else f"R16_MAIN_VS_CONTROLLER_DIFF ({len(bad)})")
