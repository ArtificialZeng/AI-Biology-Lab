"""Exact SymPy checks for Example 8, arXiv:2606.28895v2, pp.38--39.

Run with Python + SymPy. The analytic spectral and limit arguments are in PROOF.md.
No numerical eigenvalue tolerance is used. Checks use explicit exceptions.
The default only prints one JSON document and writes no reports or archives.
Use --output PATH only to explicitly create a new report (never overwrite).
"""

from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
import sympy as s


class CheckFailure(RuntimeError):
    """An exact symbolic identity failed, including under Python -O."""


def certificate() -> dict:
    k1, k2, k3, k4, km1, km2, km3, km4 = s.symbols("k1 k2 k3 k4 km1 km2 km3 km4", positive=True)
    c1, c2 = s.symbols("c1 c2", nonnegative=True)
    y1, y4, y5, y6 = s.symbols("y1 y4 y5 y6")
    rhs = s.Matrix(
        [
            k1 * (-y1 + c1 - y4 - y6) * (y1 - y4 - y5 - y6)
            - k2 * (-y1 + c2 - y5 - y6) * y4
            - km1 * y4
            + km2 * y6,
            k3 * (-y1 + c2 - y5 - y6) * (y1 - y4 - y5 - y6)
            - k4 * (-y1 + c1 - y4 - y6) * y5
            - km3 * y5
            + km4 * y6,
            k2 * (-y1 + c2 - y5 - y6) * y4 + k4 * (-y1 + c1 - y4 - y6) * y5 - km2 * y6 - km4 * y6,
        ]
    )
    actual = rhs.jacobian([y4, y5, y6]).subs({y1: 0, y4: 0, y5: 0, y6: 0})
    printed = s.Matrix(
        [
            [-k1 * c1 - k2 * c2 - km1, -k1 * c1, -k1 * c1 + km2],
            [-k2 * c2, -k3 * c2 - k4 * c1 - km3, -k3 * c2 + km4],
            [k2 * c2, k4 * c1, -km2 - km4],
        ]
    )
    delta = s.simplify(printed - actual)
    expected_delta = s.zeros(3)
    expected_delta[1, 0] = (k3 - k2) * c2
    if not delta == expected_delta:
        raise CheckFailure("Exact identity failed: delta == expected_delta")
    witness = {c1: 1, c2: 1, k1: 7, k2: 7, km2: 7, k3: 1, k4: 1, km1: 1, km3: 1, km4: 1}
    jp, ja = (printed.subs(witness), actual.subs(witness))
    lam = s.Symbol("lambda")
    cp = s.factor(jp.charpoly(lam).as_expr())
    ca = s.factor(ja.charpoly(lam).as_expr())
    if not cp == (lam + 8) * (lam**2 + 18 * lam - 4):
        raise CheckFailure("Exact identity failed: cp == (lam + 8) * (lam ** 2 + 18 * lam - 4)")
    if not ca == (lam + 8) * (lam**2 + 18 * lam + 38):
        raise CheckFailure("Exact identity failed: ca == (lam + 8) * (lam ** 2 + 18 * lam + 38)")
    eig = -9 + s.sqrt(85)
    if not eig.is_positive is True:
        raise CheckFailure("Exact identity failed: eig.is_positive is True")
    a, b, c, d, e, f, g, h = s.symbols("a b c d e f g h", nonnegative=True)
    q = s.Matrix([[-a - c, e, g, 0], [a, -e - b, 0, f], [c, 0, -g - d, h], [0, b, d, -f - h]])
    hh = s.Matrix([[1, 1, 1, 1], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    transformed = s.simplify(hh * q * hh.inv())
    if not s.ones(1, 4) * q == s.zeros(1, 4):
        raise CheckFailure("Exact identity failed: s.ones(1, 4) * q == s.zeros(1, 4)")
    if not transformed[0, :] == s.zeros(1, 4):
        raise CheckFailure("Exact identity failed: transformed[0, :] == s.zeros(1, 4)")
    substitutions = {a: k1 * c1, b: k2 * c2, c: k3 * c2, d: k4 * c1, e: km1, f: km2, g: km3, h: km4}
    if not s.simplify(transformed[1:4, 1:4].subs(substitutions) - actual) == s.zeros(3):
        raise CheckFailure(
            "Exact identity failed: s.simplify(transformed[1:4, 1:4].subs(substitutions) - actual) == s.zeros(3)"
        )
    x1, x2, x3, x4, x5, x6, eps = s.symbols("x1 x2 x3 x4 x5 x6 eps")
    xx = s.Matrix([x1, x2, x3, x4, x5, x6])
    f6 = s.Matrix(
        [
            -x1 * x3 + x4 - x1 * x5 + x6,
            -x2 * x4 + x6 - x2 * x3 + x5,
            -x1 * x3 + x4 - x2 * x3 + x5 + 2 * eps * x6,
            x1 * x3 - x4 - x2 * x4 + x6,
            -x1 * x5 + x6 + x2 * x3 - x5,
            x2 * x4 + x1 * x5 - 2 * x6 - eps * x6,
        ]
    )
    t = s.Matrix([[0, 0, 1, 1, 1, 1], [1, -1, 0, 1, -1, 0], [0, 1, 1, 1, 2, 2]])
    if not t.rank() == 3:
        raise CheckFailure("Exact identity failed: t.rank() == 3")
    if not t * f6 == s.Matrix([eps * x6, 0, 0]):
        raise CheckFailure("Exact identity failed: t * f6 == s.Matrix([eps * x6, 0, 0])")
    x0 = s.Matrix(
        [
            s.Rational(1, 2),
            s.Rational(1, 2),
            s.Rational(13, 6),
            s.Rational(13, 12),
            s.Rational(13, 12),
            s.Rational(13, 24),
        ]
    )
    if not f6.subs(dict(zip(xx, x0))).subs(eps, 0) == s.zeros(6, 1):
        raise CheckFailure(
            "Exact identity failed: f6.subs(dict(zip(xx, x0))).subs(eps, 0) == s.zeros(6, 1)"
        )
    if not t * x0 == s.Matrix([s.Rational(39, 8), 0, 7]):
        raise CheckFailure("Exact identity failed: t * x0 == s.Matrix([s.Rational(39, 8), 0, 7])")
    u, v, p, w = s.symbols("u v p w")
    sym = {x1: u, x2: u, x3: p, x4: v, x5: v, x6: w}
    if not s.simplify((f6[0] - f6[1]).subs(sym)) == 0:
        raise CheckFailure("Exact identity failed: s.simplify((f6[0] - f6[1]).subs(sym)) == 0")
    if not s.simplify((f6[3] - f6[4]).subs(sym)) == 0:
        raise CheckFailure("Exact identity failed: s.simplify((f6[3] - f6[4]).subs(sym)) == 0")
    result = {
        "status": "passed",
        "sympy_version": s.__version__,
        "printed_minus_actual": str(delta),
        "witness_parameters": {str(k): str(v) for k, v in witness.items()},
        "printed_witness_matrix": str(jp),
        "actual_witness_matrix": str(ja),
        "printed_characteristic_polynomial": str(cp),
        "actual_characteristic_polynomial": str(ca),
        "printed_positive_eigenvalue": str(eig),
        "actual_general_jacobian": str(actual),
        "generator": str(q),
        "coordinate_change": str(hh),
        "conjugated_generator": str(transformed),
        "optional_initial_condition": str(x0),
        "optional_initial_lumped_values": str(t * x0),
        "optional_alltime_error": "17/8, proved analytically in PROOF.md",
    }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, help="Explicitly create a new JSON file; refuse overwrite"
    )
    args = parser.parse_args()
    try:
        result = certificate()
        report = json.dumps(result, indent=2, ensure_ascii=False) + "\n"
        if args.output is not None:
            with args.output.open("x", encoding="utf-8") as handle:
                handle.write(report)
        sys.stdout.write(report)
        return 0
    except (CheckFailure, OSError) as error:
        sys.stderr.write(f"check failed: {error}\n")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
