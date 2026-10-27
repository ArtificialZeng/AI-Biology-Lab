#!/usr/bin/env python3
"""Exact derivative algebra and finite, labelled floating-point replay.

The universal theorem is analytic; this program does not certify all real
parameters by enumeration. Output is deterministic and never writes the archive.
"""
import json
from fractions import Fraction as Q
from pathlib import Path

import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def exact_checks():
    a, b, u, A, B = sp.symbols("a b u A B")
    D = u * A + (1-u) * B
    G = a*u*A + b*(1-u)*B
    identity = sp.expand(G - (a-b)*u*(1-u)*(A-B) - (a*u+b*(1-u))*D)
    require(identity == 0, "Implicit derivative numerator identity failed")
    gap_identity = sp.expand(a*b*(a/b+b/a-2) - (a-b)**2)
    require(gap_identity == 0, "Endpoint displacement gap identity failed")
    count = strict = 0
    slopes = []
    for aa in range(1, 6):
        for bb in range(1, 6):
            for ui in range(1, 7):
                uu = Q(ui, 7)
                for v in [Q(6, 5), Q(3, 2), Q(2)]:
                    # Construct an exact valid transfer with tau=log(v).
                    AA, BB = v**aa, v**bb
                    DD = uu*AA + (1-uu)*BB
                    delta = 1/DD
                    vv = uu*AA/DD
                    GG = aa*uu*AA + bb*(1-uu)*BB
                    slope = (aa*uu+bb*(1-uu))*DD/GG
                    direct = 1-Q(aa-bb)*uu*(1-uu)*(AA-BB)/GG
                    require(0 < delta < 1 and 0 < vv < 1, "Invalid exact witness")
                    require(slope == direct and 0 < slope <= 1, "Slope witness failed")
                    require((AA-DD)*(aa-bb) >= 0 and (DD-BB)*(aa-bb) >= 0,
                            "Waiting-time endpoint ordering failed")
                    if aa != bb:
                        require(slope < 1, "Expected strict contraction")
                        strict += 1
                    else:
                        require(vv == uu and slope == 1, "Neutral transfer failed")
                    slopes.append(slope)
                    count += 1
    return {"symbolic_identity": str(identity), "endpoint_gap_identity": str(gap_identity),
            "exact_constructed_transfers": count, "strict_transfers": strict,
            "minimum_slope": str(min(slopes)), "maximum_slope": str(max(slopes)),
            "scope": "Two polynomial identities checked symbolically; 450 rational transfer witnesses, not an exhaustive real-parameter proof."}


def mf(q):
    return mp.mpf(q.numerator)/q.denominator


def transfer(z, a, b, delta):
    """Solve original density threshold in waiting time, without a map formula."""
    u = 1/(1+mp.exp(-z))
    L = -mp.log(delta)
    lo, hi = L/max(a, b), L/min(a, b)
    if a == b:
        return z
    for _ in range(112):
        mid = (lo+hi)/2
        density = delta*(u*mp.exp(a*mid)+(1-u)*mp.exp(b*mid))
        if density < 1:
            lo = mid
        else:
            hi = mid
    tau = (lo+hi)/2
    out = z+(a-b)*tau
    eps = mp.mpf("1e-27")
    require(L*(1-b/a)-eps <= out-z <= L*(a/b-1)+eps,
            "Floating replay violates analytic one-step envelope")
    return out


def classify(rates):
    upper = sum((a/b-1 for a, b in rates), Q(0))
    lower = sum((1-b/a for a, b in rates), Q(0))
    if all(a == b for a, b in rates):
        label = "neutral"
    elif upper <= 0:
        label = "species2_exclusion"
    elif lower >= 0:
        label = "species1_exclusion"
    else:
        require(lower < 0 < upper, "Classification gap")
        label = "interior_cycle"
    return label, lower, upper


def numerical_checks():
    mp.mp.dps = 60
    model = json.loads((ROOT/"model.json").read_text())
    delta_q = Q(model["dilution"])
    require(0 < delta_q < 1 and Q(model["threshold"]) == 1,
            "Invalid model dilution/threshold")
    delta = mf(delta_q)
    L = -mp.log(delta)
    eps = mp.mpf("1e-26")
    records = []
    transfer_count = 0
    for case in model["examples"]:
        rates_q = [(Q(a), Q(b)) for a, b in case["rates"]]
        require(rates_q and all(a > 0 and b > 0 for a, b in rates_q), "Rates must be positive")
        label, lower, upper = classify(rates_q)
        per_phase = []
        for phase in range(len(rates_q)):
            rates = rates_q[phase:]+rates_q[:phase]
            finals = []
            for initial in [-8, -1, 0, 3, 8]:
                z = mp.mpf(initial)
                previous_displacement = None
                for period in range(18):
                    old = z
                    for a, b in rates:
                        z = transfer(z, mf(a), mf(b), delta)
                        transfer_count += 1
                    displacement = z-old
                    require(L*mf(lower)-eps <= displacement <= L*mf(upper)+eps,
                            "Period envelope failed")
                    if label == "species2_exclusion":
                        require(displacement < eps, "Wrong exclusion direction")
                    elif label == "species1_exclusion":
                        require(displacement > -eps, "Wrong exclusion direction")
                    elif label == "neutral":
                        require(abs(displacement) <= eps, "Neutrality failed")
                    elif previous_displacement is not None:
                        require(abs(displacement) <= abs(previous_displacement)+eps,
                                "Interior-cycle displacement expanded")
                        require(displacement*previous_displacement >= -eps,
                                "Increasing period map overshot fixed point")
                    previous_displacement = displacement
                    require(mpf_envelope(initial, period+1, lower, upper, L, z, eps),
                            "Cumulative envelope failed")
                finals.append({"initial_log_odds": initial, "final_log_odds": mp.nstr(z, 20)})
            per_phase.append({"starting_phase": phase, "runs": finals})
        records.append({"name": case["name"], "classification": label,
                        "V_over_L": str(lower), "U_over_L": str(upper),
                        "phases": per_phase})
    return {"kind": "floating-point-observation", "decimal_precision": 60,
            "waiting_time_bisections": 112, "periods_per_run": 18,
            "transfers_replayed": transfer_count, "cases": records,
            "scope": "Seven finite examples, every phase and five initial log odds; not proof of infinite-time convergence or universal quantifiers."}


def mpf_envelope(initial, n, lower, upper, L, z, eps):
    return initial+n*L*mf(lower)-eps <= z <= initial+n*L*mf(upper)+eps


def compute():
    return {"schema": "serial-transfer-owner-certificate-v1",
            "exact": exact_checks(), "numerical": numerical_checks()}


if __name__ == "__main__":
    print(json.dumps(compute(), sort_keys=True, indent=2))
