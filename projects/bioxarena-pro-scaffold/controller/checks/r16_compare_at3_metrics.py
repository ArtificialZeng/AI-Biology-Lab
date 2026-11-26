"""Controller: staged verification step 3b. A-T3's own metric JSONs ({model}_{ckpt}__{line}_metrics.json, from A-T3's
independent scoring code) against the main line's round16/results/summary_r16.json, for one model.

Addendum A tolerance for cross-party metrics: |diff| <= 1e-5 on every estimate and CI bound; categorical calls
(per-source H-I1 categories, stratum summaries, H-I1 line verdict, H-I2 call) identical.
Main-line location: STATE_SE -> lines/<line>/primary/<...>; STATE_HVG -> lines/<line>/secondary/STATE_HVG/<...>.
The H-I2 call is not stored in the main summary for STATE_HVG, so it is derived from the main 97.5% CI with the
PREREG_I rule (lo > 0 MODEL_SUPERIOR_TO_FCT; hi < 0 FCT_SUPERIOR_TO_MODEL; else INDISTINGUISHABLE).
Usage: python3 -I r16_compare_at3_metrics.py <summary_r16.json> <at3_dir> <STATE_HVG|STATE_SE> [lines...]
"""
import json
import sys
from pathlib import Path

TOL = 1e-5
S = json.load(open(sys.argv[1]))["lines"]
D = Path(sys.argv[2])
model = sys.argv[3]
lines = sys.argv[4:] or ["nadig_hepg2", "nadig_jurkat", "a549_ko", "dld1", "cd4t_rest"]
worst, bad, n = 0.0, [], 0
BOUNDS = ["estimate", "ci95_lo", "ci95_hi", "ci975_lo", "ci975_hi", "valid_frac"]


def num(tag, a, b):
    global worst, n
    n += 1
    d = abs(float(a) - float(b)); worst = max(worst, d)
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


def h2_call(m):
    if m["ci975_lo"] > 0:
        return "MODEL_SUPERIOR_TO_FCT"
    if m["ci975_hi"] < 0:
        return "FCT_SUPERIOR_TO_MODEL"
    return "INDISTINGUISHABLE"


for ln in lines:
    f = sorted(D.glob(f"{model}_*__{ln}_metrics.json"))
    if not f:
        bad.append(f"{ln}: no A-T3 metrics file"); continue
    a = json.load(open(f[0]))["results"]
    if model == "STATE_HVG":
        m = S[ln]["secondary"]["STATE_HVG"]
        keys = {"matched": "r_orth_cal_matched", "dmatched": "r_dorth_cal_matched", "fct": "r_orth_cal_vs_fct",
                "paired": "r_orth_cal_paired_vs_mismatched", "h1": "r_orth_cal_descriptive_h_i1"}
        mres = m["within_line_resample_reference"]["r_orth_cal"]
    else:
        raise SystemExit("STATE_SE metric layout: use r16_compare_main.py on the controller recompute")
    summ(f"{ln}.matched", m[keys["matched"]], a["orth_matched"])
    summ(f"{ln}.dorth_matched", m[keys["dmatched"]], a["dorth_matched"])
    summ(f"{ln}.vs_fct", m[keys["fct"]], a["vs_fct_orth"])
    eq(f"{ln}.h_i2_call", h2_call(m[keys["fct"]]), a["h_i2"]["classification"])
    for k in ["estimate", "ci975_lo", "ci975_hi"]:
        num(f"{ln}.h_i2.{k}", m[keys["fct"]][k], a["h_i2"][k])
    mp, ap = m[keys["paired"]]["per_source"], a["orth_paired"]["per_source"]
    eq(f"{ln}.paired_sources", sorted(mp), sorted(ap))
    for s in mp:
        if s in ap:
            summ(f"{ln}.paired.{s}", mp[s], ap[s])
    mh, ah = m[keys["h1"]], a["h_i1"]
    eq(f"{ln}.h_i1.line_summary", mh["line_summary"], ah["line_summary"])
    eq(f"{ln}.h_i1.counts", mh["counts"], ah["counts"])
    for st in mh["strata"]:
        ms, as_ = mh["strata"][st], ah["strata"].get(st, {})
        for k in ["n_pairs", "sources", "counts", "all_within_delta", "any_within_delta"]:
            eq(f"{ln}.h_i1.{st}.{k}", ms.get(k), as_.get(k))
        for k in ms:
            if k.endswith("summary") or k in ("stratum_summary", "verdict", "category"):
                eq(f"{ln}.h_i1.{st}.{k}", ms[k], as_.get(k))
    for s in mh["per_source"]:
        mps, aps = mh["per_source"][s], ah["per_source"].get(s, {})
        eq(f"{ln}.h_i1.{s}.category", mps["category"], aps.get("category"))
        eq(f"{ln}.h_i1.{s}.within_delta", mps["within_delta"], aps.get("within_delta"))
        for k in ["estimate", "ci975_lo", "ci975_hi"]:
            num(f"{ln}.h_i1.{s}.{k}", mps[k], aps[k])
    summ(f"{ln}.resample_delta", mres["delta_pair_resample"], a["resample_orth"]["delta_pair_resample"])

print(json.dumps({"model": model, "lines": lines, "n_checks": n, "n_fail": len(bad), "worst_abs_diff": worst}))
for b in bad:
    print("FAIL", b)
print(f"R16_AT3_METRICS_MATCH {model}" if not bad else f"R16_AT3_METRICS_DIFF {model}")
