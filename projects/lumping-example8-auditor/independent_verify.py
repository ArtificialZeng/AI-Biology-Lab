#!/usr/bin/env python3
"""Read-only, exact independent checks for arXiv:2606.28895v2 Example 8.

Requires Python 3.10+ and SymPy. No network, coordinate-equation verifier, tolerance, subprocess,
or file input is used. The default only prints JSON. --output PATH explicitly
requests creation of a new JSON report, and refuses to overwrite any file.

The reaction model and coordinate definitions are due to Justin Eilertsen,
Valery G. Romanovski, Santiago Schnell and Sebastian Walcher; their model
follows Gijima and Peacock-Lopez. This code checks a source-specific correction,
not a new general lumping theorem. Analytic reasoning and limitations are in
PROOF.md. This is a symbolic certificate, not a proof-assistant formalization.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import sympy as sp


class CheckFailure(RuntimeError):
    """An exact identity or a certificate condition did not hold."""


def require(condition: bool, message: str) -> None:
    # Explicit exceptions continue to work under python -O.
    if not condition:
        raise CheckFailure(message)


def equal(actual, expected, label: str) -> None:
    difference = actual - expected
    if isinstance(difference, sp.MatrixBase):
        ok = all(sp.expand(entry) == 0 for entry in difference)
    else:
        ok = sp.expand(difference) == 0
    require(ok, label)


def certificate() -> dict:
    x = sp.Matrix(sp.symbols("x1:7"))
    x1, x2, x3, x4, x5, x6 = x
    rates = sp.symbols("k1 km1 k2 km2 k3 km3 k4 km4 k5 km5")
    k1, km1, k2, km2, k3, km3, k4, km4, k5, km5 = rates

    # A separate transcription of reactions (27)--(31), rather than the other verifier's
    # three coordinate equations: each pair is (reactant vector, product vector).
    pairs = [
        ((1, 0, 1, 0, 0, 0), (0, 0, 0, 1, 0, 0)),
        ((0, 1, 0, 1, 0, 0), (0, 0, 0, 0, 0, 1)),
        ((0, 1, 1, 0, 0, 0), (0, 0, 0, 0, 1, 0)),
        ((1, 0, 0, 0, 1, 0), (0, 0, 0, 0, 0, 1)),
        ((0, 0, 0, 0, 0, 1), (0, 0, 2, 0, 0, 0)),
    ]
    rhs = sp.zeros(6, 1)
    for j, (left, right) in enumerate(pairs):
        reactant, product = sp.Matrix(left), sp.Matrix(right)
        forward = rates[2 * j] * sp.prod(x[i] ** left[i] for i in range(6))
        reverse = rates[2 * j + 1] * sp.prod(x[i] ** right[i] for i in range(6))
        rhs += (product - reactant) * (forward - reverse)
    rhs = sp.expand(rhs)

    # Independent use of the full six-dimensional linear coordinate map.
    T = sp.Matrix(
        [
            [0, 0, 1, 1, 1, 1],
            [1, -1, 0, 1, -1, 0],
            [0, 1, 1, 1, 2, 2],
            [0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 1],
        ]
    )
    require(T.det() != 0, "coordinate map must be invertible")
    y = sp.Matrix(sp.symbols("y1:7"))
    transformed = sp.expand((T * rhs).subs(dict(zip(x, T.inv() * y)), simultaneous=True))
    c1, c2 = sp.symbols("c1 c2")
    reduced_rhs = transformed.subs({y[1]: c1 - c2, y[2]: c2, k5: 0, km5: 0})
    J = sp.simplify(
        reduced_rhs[3:6, 0].jacobian(y[3:6, 0]).subs({y[0]: 0, y[3]: 0, y[4]: 0, y[5]: 0})
    )
    expected = sp.Matrix(
        [
            [-k1 * c1 - k2 * c2 - km1, -k1 * c1, -k1 * c1 + km2],
            [-k3 * c2, -k3 * c2 - k4 * c1 - km3, -k3 * c2 + km4],
            [k2 * c2, k4 * c1, -km2 - km4],
        ]
    )
    equal(J, expected, "reaction-derived Jacobian must equal corrected matrix")
    printed = expected.copy()
    printed[1, 0] = -k2 * c2
    delta = sp.zeros(3)
    delta[1, 0] = (k3 - k2) * c2
    equal(printed - J, delta, "only printed (2,1) differs")

    witness = {c1: 1, c2: 1, k1: 7, km1: 1, k2: 7, km2: 7, k3: 1, km3: 1, k4: 1, km4: 1}
    require(all(v > 0 for v in witness.values()), "witness fast rates and totals positive")
    z = sp.Symbol("z")
    jp, ja = printed.subs(witness), J.subs(witness)
    cp = sp.factor(jp.charpoly(z).as_expr())
    ca = sp.factor(ja.charpoly(z).as_expr())
    equal(cp, (z + 8) * (z * z + 18 * z - 4), "printed witness polynomial")
    equal(ca, (z + 8) * (z * z + 18 * z + 38), "corrected witness polynomial")
    require(85 > 9**2, "-9+sqrt(85)>0 exactly")
    # At y=0 this witness represents x=(1,1,0,0,0,0), a physical state.
    equal(
        (T.inv() * y).subs({y[0]: 0, y[1]: 0, y[2]: 1, y[3]: 0, y[4]: 0, y[5]: 0}),
        sp.Matrix([1, 1, 0, 0, 0, 0]),
        "physical witness equilibrium",
    )

    # Alternate algebraic proof certificate, independent of generator spectral facts.
    v = sp.symbols("a b c d e f g h")
    a, b, c, d, e, f, g, h = v
    substitutions = {a: k1 * c1, b: k2 * c2, c: k3 * c2, d: k4 * c1, e: km1, f: km2, g: km3, h: km4}
    jabs = sp.Matrix([[-a - b - e, -a, -a + f], [-c, -c - d - g, -c + h], [b, d, -f - h]])
    equal(jabs.subs(substitutions), J, "abstract matrix parameter substitution")
    coeffs = sp.Poly(jabs.charpoly(z).as_expr(), z).all_coeffs()
    equal(coeffs[0], 1, "monic cubic")
    A1, A2, A3 = coeffs[1:]
    margins = [A1, A2, A3, sp.expand(A1 * A2 - A3)]
    names = ["A1", "A2", "A3", "A1*A2-A3"]
    boundaries = [
        e + f + g + h,
        e * f + e * g + e * h + f * g + g * h,
        e * g * (f + h),
        (e + g) * (e + f + h) * (f + g + h),
    ]
    positivity = {}
    for name, expr, base in zip(names, margins, boundaries):
        polynomial = sp.Poly(expr, *v)
        require(
            all(coefficient > 0 for coefficient in polynomial.coeffs()),
            name + " has a nonpositive coefficient",
        )
        equal(expr.subs({a: 0, b: 0, c: 0, d: 0}), base, name + " positive boundary polynomial")
        positivity[name] = {
            "monomials": len(polynomial.terms()),
            "minimum_coefficient": str(min(polynomial.coeffs())),
            "boundary": str(sp.factor(base)),
            "expanded": str(sp.expand(expr)),
        }
    # All a,b,c,d >= 0 and e,f,g,h > 0: nonnegative additional terms plus
    # the positive displayed boundary polynomials satisfy cubic Routh-Hurwitz.

    # Check the author's generator embedding as a second route.
    Q = sp.Matrix([[-a - c, e, g, 0], [a, -e - b, 0, f], [c, 0, -g - d, h], [0, b, d, -f - h]])
    H = sp.Matrix([[1, 1, 1, 1], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    qcoords = H * Q * H.inv()
    equal(sp.ones(1, 4) * Q, sp.zeros(1, 4), "generator column sums")
    equal(qcoords[0, :], sp.zeros(1, 4), "conserved total row")
    equal(qcoords[1:4, 1:4], jabs, "generator transverse block")
    equal(Q.charpoly(z).as_expr(), z * jabs.charpoly(z).as_expr(), "generator spectral factor")

    # Optional all-time example: independent reaction-derived identities only.
    eps = sp.Symbol("epsilon")
    slow = rhs.subs({**{r: 1 for r in rates[:8]}, k5: eps, km5: 0})
    x0 = sp.Matrix(
        [
            sp.Rational(1, 2),
            sp.Rational(1, 2),
            sp.Rational(13, 6),
            sp.Rational(13, 12),
            sp.Rational(13, 12),
            sp.Rational(13, 24),
        ]
    )
    equal(slow.subs(dict(zip(x, x0))).subs(eps, 0), sp.zeros(6, 1), "critical initial equilibrium")
    equal((T * slow)[:3, 0], sp.Matrix([eps * x6, 0, 0]), "lumped observable and integrals")
    equal((T * x0)[:3, 0], sp.Matrix([sp.Rational(39, 8), 0, 7]), "initial lumped values")
    mass = sp.Matrix([[1, 1, 2, 3, 3, 4]])
    equal((mass * slow)[0], 0, "positive mass conservation")
    equal((mass * x0)[0], 14, "initial mass")
    u, p, vv, w = sp.symbols("u p v w")
    symmetric = {x1: u, x2: u, x3: p, x4: vv, x5: vv, x6: w}
    equal((slow[0] - slow[1]).subs(symmetric), 0, "symmetry x1=x2 invariant")
    equal((slow[3] - slow[4]).subs(symmetric), 0, "symmetry x4=x5 invariant")
    symrhs = slow.subs(symmetric)
    f4 = sp.Matrix([symrhs[0], symrhs[2], symrhs[3], symrhs[5]])
    equal(
        f4,
        sp.Matrix(
            [
                -u * p + vv - u * vv + w,
                -2 * u * p + 2 * vv + 2 * eps * w,
                u * p - vv - u * vv + w,
                2 * u * vv - (2 + eps) * w,
            ]
        ),
        "four-dimensional symmetric system",
    )
    equal((sp.Matrix([[1, 1, 3, 2]]) * f4)[0], 0, "symmetric conserved mass")
    equal((sp.Matrix([[0, 1, 2, 1]]) * f4)[0], eps * w, "monotone h derivative")
    equal(7 - sp.Rational(39, 8), sp.Rational(17, 8), "exact supremum difference")
    # The implication from these identities to h(t)->7 uses the analytic
    # invariant-set argument in PROOF.md; no numerical integration claims it.

    return {
        "status": "passed",
        "sympy_version": sp.__version__,
        "scope": "source-specific displayed-Jacobian correction and exact algebra; analytic limits reviewed separately",
        "source": "https://arxiv.org/abs/2606.28895v2",
        "coordinate_equation_verifier_executed": False,
        "corrected_jacobian": str(J),
        "printed_witness_matrix": str(jp),
        "actual_witness_matrix": str(ja),
        "printed_characteristic_polynomial": str(cp),
        "actual_characteristic_polynomial": str(ca),
        "printed_positive_eigenvalue": "-9 + sqrt(85)",
        "routh_hurwitz_certificate": positivity,
        "generator_identity": "det(zI-Q) = z det(zI-J)",
        "parameter_domain": "a,b,c,d>=0; e,f,g,h>0, hence eight fast rates>0 and c1,c2>=0",
        "optional_slow_drift": {
            "h0": "39/8",
            "h_limit": "7, analytic invariant-set proof required",
            "supremum_error": "17/8 for each epsilon>0",
        },
        "formal_proof_assistant": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, help="explicitly create this new JSON file; never overwrite"
    )
    args = parser.parse_args()
    try:
        result = certificate()
        report = json.dumps(result, indent=2, ensure_ascii=False) + "\n"
        if args.output is not None:
            # Exclusive creation is an intentional write exception to read-only default.
            with args.output.open("x", encoding="utf-8") as handle:
                handle.write(report)
        sys.stdout.write(report)
        return 0
    except (CheckFailure, OSError) as error:
        sys.stderr.write(f"check failed: {error}\n")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
