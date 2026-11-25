"""Controller (third-party) implementation of PREREG_I v4 §3 adjudication, written from the PREREG text only.

It never imports round16/src or round16/audit. Inputs are plain dicts of bootstrap summaries:
  pair / contrast: {"valid_frac": float, "ci975": (lo, hi), "ci95": (lo, hi)}
  matched arm:     {"valid_frac": float, "rho_ab_orth": float, "ci975": (lo, hi)}
Run `python3 -I r16_check_adjudicate.py` to execute the built-in branch tests.
"""
import math

DELTA = 0.04
MIN_VALID = 0.95


def _finite(ci):
    return ci is not None and len(ci) == 2 and all(isinstance(v, (int, float)) and math.isfinite(v) for v in ci)


def within_delta(p, delta=DELTA):
    """§3.1: valid_frac >= 0.95, both CIs finite, and two-sided ci95 strictly inside (-delta, +delta)."""
    return bool(p["valid_frac"] >= MIN_VALID and _finite(p["ci975"]) and _finite(p["ci95"])
                and p["ci95"][0] > -delta and p["ci95"][1] < delta)


def classify_pair(p, delta=DELTA):
    """§3.1 four-way call for one (line, mismatched source) pair."""
    if p["valid_frac"] < MIN_VALID or not _finite(p["ci975"]):
        return "INCONCLUSIVE"
    lo, hi = p["ci975"]
    if lo == 0.0 or hi == 0.0:  # §3.1 rule 5: exact ties at 0 are INCONCLUSIVE (controller reading; v5 must state it)
        return "INCONCLUSIVE"
    if lo > 0.0:
        return "MATCHED_BETTER"
    if hi < 0.0:
        return "MISMATCHED_BETTER"
    if within_delta(p, delta):
        return "EQUIVALENT"
    return "INCONCLUSIVE"


def matched_arm_valid(m):
    return bool(m["valid_frac"] >= MIN_VALID and isinstance(m["rho_ab_orth"], (int, float))
                and math.isfinite(m["rho_ab_orth"]) and m["rho_ab_orth"] > 0.0 and _finite(m["ci975"]))


def adjudicate_line_h_i1(matched, seen, unseen, delta=DELTA):
    """§3.1 line-level priority order. seen/unseen: {source: pair dict}."""
    if len(seen) == 0 or len(unseen) == 0:
        return "MISSING_STRATUM"
    if not matched_arm_valid(matched):
        return "INCONCLUSIVE"
    pos = matched["ci975"][0] > 0.0
    calls = [classify_pair(p, delta) for p in list(seen.values()) + list(unseen.values())]
    all_wd = all(within_delta(p, delta) for p in list(seen.values()) + list(unseen.values()))
    if pos and all(c == "MATCHED_BETTER" for c in calls):
        return "LINE_MATCHED_SPECIFIC_SUB_DELTA" if all_wd else "LINE_MATCHED_SPECIFIC"
    if all(c == "EQUIVALENT" for c in calls):
        return "LINE_EQUIVALENT_NONZERO_SIGNAL" if pos else "LINE_EQUIVALENT_NULL_SIGNAL"
    if all(c == "MISMATCHED_BETTER" for c in calls):
        return "LINE_MISMATCHED_BETTER"
    return "MIXED_OR_INCONCLUSIVE"


def classify_h_i2(c):
    """§3.2 model-vs-FCT three-way call."""
    if c["valid_frac"] >= MIN_VALID and _finite(c["ci975"]):
        if c["ci975"][0] > 0.0:
            return "MODEL_SUPERIOR_TO_FCT"
        if c["ci975"][1] < 0.0:
            return "FCT_SUPERIOR_TO_MODEL"
    return "INDISTINGUISHABLE"


def group_report(line_verdicts, group2_qualifying):
    """§3.1 group rule: k/K per group; Group 2 with < 2 qualifying primary lines is descriptive only."""
    out = {}
    for g, verdicts in line_verdicts.items():
        counts = {}
        for v in verdicts.values():
            counts[v] = counts.get(v, 0) + 1
        out[g] = {"K": len(verdicts), "counts": counts}
    if group2_qualifying < 2:
        for g in out:
            if g.startswith("group2"):
                out[g]["status"] = "SINGLE_STRATUM_DESCRIPTIVE_ONLY"
    return out


def _p(lo975, hi975, lo95=None, hi95=None, vf=1.0):
    lo95 = lo975 if lo95 is None else lo95
    hi95 = hi975 if hi95 is None else hi95
    return {"valid_frac": vf, "ci975": (lo975, hi975), "ci95": (lo95, hi95)}


def _tests():
    M_ok = {"valid_frac": 1.0, "rho_ab_orth": 0.5, "ci975": (0.05, 0.12)}
    M_null = {"valid_frac": 1.0, "rho_ab_orth": 0.5, "ci975": (-0.01, 0.03)}
    eq = _p(-0.02, 0.025, -0.015, 0.02)
    mb = _p(0.03, 0.09, 0.035, 0.085)
    mb_sub = _p(0.004, 0.026, 0.006, 0.024)
    mis = _p(-0.09, -0.03, -0.085, -0.035)
    wide = _p(-0.05, 0.06, -0.045, 0.055)
    T = [
        (classify_pair(eq), "EQUIVALENT"), (classify_pair(mb), "MATCHED_BETTER"), (classify_pair(mis), "MISMATCHED_BETTER"),
        (classify_pair(wide), "INCONCLUSIVE"), (classify_pair(_p(-0.02, 0.02, vf=0.9)), "INCONCLUSIVE"),
        (classify_pair(_p(float("nan"), 0.02)), "INCONCLUSIVE"), (classify_pair(_p(0.0, 0.03, -0.001, 0.03)), "INCONCLUSIVE"),
        (classify_pair(_p(-0.03, 0.04, -0.02, 0.04)), "INCONCLUSIVE"),  # open boundary at +delta
        (within_delta(_p(-0.02, 0.02, vf=0.9)), False), (within_delta(mb_sub), True),
        (adjudicate_line_h_i1(M_ok, {}, {"a": eq}), "MISSING_STRATUM"),
        (adjudicate_line_h_i1({"valid_frac": 1.0, "rho_ab_orth": 0.0, "ci975": (0.05, 0.1)}, {"s": eq}, {"u": eq}), "INCONCLUSIVE"),
        (adjudicate_line_h_i1({"valid_frac": 0.9, "rho_ab_orth": 0.5, "ci975": (0.05, 0.1)}, {"s": eq}, {"u": eq}), "INCONCLUSIVE"),
        (adjudicate_line_h_i1(M_ok, {"s": mb}, {"u": mb, "v": mb}), "LINE_MATCHED_SPECIFIC"),
        (adjudicate_line_h_i1(M_ok, {"s": mb_sub}, {"u": mb_sub}), "LINE_MATCHED_SPECIFIC_SUB_DELTA"),
        (adjudicate_line_h_i1(M_ok, {"s": mb_sub}, {"u": mb}), "LINE_MATCHED_SPECIFIC"),
        (adjudicate_line_h_i1(M_ok, {"s": eq}, {"u": eq}), "LINE_EQUIVALENT_NONZERO_SIGNAL"),
        (adjudicate_line_h_i1(M_null, {"s": eq}, {"u": eq}), "LINE_EQUIVALENT_NULL_SIGNAL"),
        (adjudicate_line_h_i1(M_null, {"s": mb}, {"u": mb}), "MIXED_OR_INCONCLUSIVE"),  # all MB but matched arm not > 0
        (adjudicate_line_h_i1(M_ok, {"s": mis}, {"u": mis}), "LINE_MISMATCHED_BETTER"),
        (adjudicate_line_h_i1(M_ok, {"s": mb}, {"u": eq}), "MIXED_OR_INCONCLUSIVE"),
        (adjudicate_line_h_i1(M_ok, {"s": eq}, {"u": wide}), "MIXED_OR_INCONCLUSIVE"),
        (classify_h_i2(_p(0.01, 0.05)), "MODEL_SUPERIOR_TO_FCT"), (classify_h_i2(_p(-0.05, -0.01)), "FCT_SUPERIOR_TO_MODEL"),
        (classify_h_i2(_p(-0.01, 0.05)), "INDISTINGUISHABLE"), (classify_h_i2(_p(0.01, 0.05, vf=0.9)), "INDISTINGUISHABLE"),
        (group_report({"group1": {"a": "X"}, "group2_crispri": {"b": "Y"}}, 1)["group2_crispri"].get("status"), "SINGLE_STRATUM_DESCRIPTIVE_ONLY"),
        (group_report({"group2_crispri": {"b": "Y", "c": "Y"}}, 2)["group2_crispri"].get("status"), None),
    ]
    bad = [(i, got, exp) for i, (got, exp) in enumerate(T) if got != exp]
    print(f"r16_check_adjudicate: {len(T) - len(bad)}/{len(T)} branch tests pass", "" if not bad else bad)
    return not bad


if __name__ == "__main__":
    raise SystemExit(0 if _tests() else 1)
