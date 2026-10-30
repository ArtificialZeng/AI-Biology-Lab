#!/usr/bin/env python3
"""Exact reaction reconstruction, elimination and stability identities.

This supplies finite/topological and symbolic certificates; the analytic
all-positive/global-stability argument is stated separately in proof.md.
No private paths, network, floating-point inference or file writes.
"""
from itertools import product
import json
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def zero(value, message):
    require(s.cancel(value) == 0, message)


def compute():
    model = json.loads((ROOT/"model.json").read_text())
    require(model.get("schema") == "groebner-positive-regular-crn-v1"
            and model.get("species") == ["X", "Y"]
            and len(model.get("reactions", [])) == 6, "Expected the frozen six-reaction two-species model")
    x, y, t = s.symbols("x y t")
    a, b, c, d, p, q = rates = s.symbols("a b c d p q")
    names = dict(zip("abcdpq", rates))
    require({r["rate"] for r in model["reactions"]} == set(names), "Independent rate labels differ")
    columns, flows = [], []
    for reaction in model["reactions"]:
        before, after = reaction["reactant"], reaction["product"]
        require(len(before) == len(after) == 2 and all(isinstance(v, int) and v >= 0 for v in before+after), "Bad complex")
        require(sum(before) <= 2 and sum(after) <= 2, "Not bimolecular")
        columns.append(s.Matrix([v-u for u, v in zip(before, after)]))
        flows.append(names[reaction["rate"]]*x**before[0]*y**before[1])
    S = s.Matrix.hstack(*columns)
    field = s.simplify(S*s.Matrix(flows))
    A, B = c+p, d+q
    Delta = s.expand(A*B-p*q)
    u, v = (B*a+q*b)/Delta, (p*a+A*b)/Delta
    expected = s.Matrix([a-A*x*x+q*y*y, b+p*x*x-B*y*y])
    for actual, target in zip(field, expected):
        zero(actual-target, "Reaction reconstruction differs")
    require(S.rank() == 2 and len(S.T.nullspace()) == 0, "Conservation-law rank differs")
    intersections = []
    # Unique unordered nontrivial partitions: column0 always on the left.
    for mask in range((1 << 5)-1):
        left = [0]+[i+1 for i in range(5) if mask & (1 << i)]
        right = [i for i in range(6) if i not in left]
        r1 = s.Matrix.hstack(*(columns[i] for i in left)).rank()
        r2 = s.Matrix.hstack(*(columns[i] for i in right)).rank()
        intersection = r1+r2-S.rank()
        require(intersection > 0, "Independent stoichiometric decomposition found")
        intersections.append(intersection)
    require(len(intersections) == 31, "Incomplete partition domain")

    f1, f2 = field
    zero((-B*f1-q*f2)/Delta-(x*x-u), "First ideal recovery failed")
    zero((-p*f1-A*f2)/Delta-(y*y-v), "Second ideal recovery failed")
    domain = s.QQ.frac_field(*rates)
    bases = []
    for order in [(x, y), (y, x)]:
        gb = s.groebner([f1, f2], *order, order="lex", domain=domain)
        want = [order[0]**2-(u if order[0] == x else v),
                order[1]**2-(u if order[1] == x else v)]
        require(len(gb.polys) == 2, "Unexpected generic basis length")
        for got, target in zip(gb.polys, want):
            zero(got.as_expr()-target, "Generic native basis differs")
        bases.append({"order": [str(z) for z in order], "degrees": [2, 2],
                      "coordinate_separating": False})

    J = field.jacobian([x, y])
    zero(J.det()-4*Delta*x*y, "Jacobian determinant identity failed")
    saturation = 1-(x*y)*(x*y)/(u*v)+(y*y*(x*x-u)+u*(y*y-v))/(u*v)
    zero(saturation, "Torus inverse certificate failed")
    M = s.Matrix([[-A, q], [p, -B]])
    DM = s.diag(p, q)*M
    require(DM == DM.T, "Weighted matrix not symmetric")
    zero(DM.det()-p*q*Delta, "Lyapunov determinant differs")
    require(s.Poly(Delta, c, d, p, q).terms() ==
            s.Poly(c*d+c*q+p*d, c, d, p, q).terms(), "Positive Delta expansion differs")

    U, V, lam = s.symbols("U V lam")
    D0 = U-lam*lam*V
    g = t**4-2*(U+lam*lam*V)*t*t+D0*D0
    X = ((3*U+lam*lam*V)*t-t**3)/(2*D0)
    Y = (t**3-(U+3*lam*lam*V)*t)/(2*lam*D0)
    for value in [X*X-U, Y*Y-V]:
        numerator = s.cancel(value).as_numer_denom()[0]
        zero(s.rem(numerator, g, t), "Mixed-coordinate square identity failed")
    zero(X+lam*Y-t, "Mixed-coordinate sum failed")

    # All64 corners are supplementary rational cases, not the all-rate proof.
    cases = []
    inputs = list(product([s.Rational(1), s.Rational(2)], repeat=6))
    inputs += [(s.Rational(1, 7), s.Rational(11, 5), s.Rational(3), s.Rational(2, 3), s.Rational(5, 2), s.Rational(13, 11)),
               (s.Rational(43), s.Rational(7), s.Rational(1), s.Rational(1), s.Rational(1, 10), s.Rational(1, 10))]
    for values in inputs:
        sub = dict(zip(rates, values))
        uu, vv, dd = map(lambda z: s.cancel(z.subs(sub)), [u, v, Delta])
        require(uu > 0 and vv > 0 and dd > 0, "Positive equilibrium squares failed")
        for order in [(x, y), (y, x)]:
            gb = s.groebner([f1.subs(sub), f2.subs(sub)], *order, order="lex", domain=s.QQ)
            target = s.groebner([x*x-uu, y*y-vv], *order, order="lex", domain=s.QQ)
            require(gb == target, "Rational specialized ideal differs")
        cases.append({"rates": [str(z) for z in values], "u": str(uu), "v": str(vv), "Delta": str(dd)})
    unit = {r: 1 for r in rates}
    unit_J = J.subs(unit).subs({x: 1, y: 1})
    require(unit_J == s.Matrix([[-4, 2], [2, -4]]), "Unit witness Jacobian differs")
    require(unit_J.eigenvals() == {s.Integer(-2): 1, s.Integer(-6): 1}, "Unit witness stability differs")
    return {"schema": "groebner-positive-owner-certificate-v1",
            "stoichiometric_columns": [[int(v) for v in col] for col in columns],
            "stoichiometric_rank": 2, "conservation_laws": 0,
            "unordered_partitions_checked": 31,
            "partition_intersection_dimensions": intersections,
            "all_rate_symbolic_checks": ["ODE reconstructed", "two ideal recovery identities",
                 "both native reduced lex bases", "Jacobian determinant4Delta*x*y",
                 "torus inverse/saturation", "symmetric Lyapunov matrix determinantpqDelta",
                 "positive Delta=cd+cq+pd", "explicit mixed-coordinate identities"],
            "generic_native_bases": bases,
            "unit_witness_J": [[int(v) for v in row] for row in unit_J.tolist()],
            "unit_witness_eigenvalues": ["-6", "-2"],
            "exact_positive_cases": len(cases), "cases": cases,
            "scope": "Symbolic identities and a complete31-partition topological test, plus66finite rational cases. Global attraction and all-positive geometric claims require the analytic proof; no biological experiment or novelty certification."}


if __name__ == "__main__":
    print(json.dumps(compute(), sort_keys=True, indent=2))
