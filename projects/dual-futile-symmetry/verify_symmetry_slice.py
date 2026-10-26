#!/usr/bin/env python3
"""Exact, standard-library-only verification of a symmetry slice of arXiv:2608.27081v2.

This checks polynomial identities, not novelty or the unrestricted conjecture.
No floating-point arithmetic, external packages, sampling or eigensolver is used.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from itertools import permutations

NAMES = ("a", "b", "e", "f", "u", "w", "v", "z", "lambda")
ZERO_EXP = (0,) * len(NAMES)


class Poly:
    def __init__(self, value=0):
        if isinstance(value, Poly):
            self.terms = value.terms.copy()
        elif isinstance(value, dict):
            self.terms = {k: Fraction(v) for k, v in value.items() if v}
        else:
            self.terms = {ZERO_EXP: Fraction(value)} if value else {}

    def __add__(self, other):
        result = self.terms.copy()
        for powers, coefficient in Poly(other).terms.items():
            result[powers] = result.get(powers, 0) + coefficient
        return Poly(result)

    __radd__ = __add__

    def __neg__(self):
        return Poly({p: -c for p, c in self.terms.items()})

    def __sub__(self, other):
        return self + -Poly(other)

    def __rsub__(self, other):
        return Poly(other) + -self

    def __mul__(self, other):
        result = {}
        for p, c in self.terms.items():
            for q, d in Poly(other).terms.items():
                powers = tuple(x + y for x, y in zip(p, q))
                result[powers] = result.get(powers, 0) + c * d
        return Poly(result)

    __rmul__ = __mul__

    def __pow__(self, power):
        if not isinstance(power, int) or power < 0:
            raise ValueError("nonnegative integer power required")
        result = Poly(1)
        for _ in range(power):
            result = result * self
        return result

    def __eq__(self, other):
        return self.terms == Poly(other).terms

    def coefficient(self, variable, exponent):
        index = NAMES.index(variable)
        result = {}
        for powers, c in self.terms.items():
            if powers[index] == exponent:
                reduced = list(powers)
                reduced[index] = 0
                result[tuple(reduced)] = c
        return Poly(result)

    def evaluate(self, values):
        total = Fraction(0)
        for powers, c in self.terms.items():
            term = c
            for name, exponent in zip(NAMES, powers):
                if exponent:
                    term *= Fraction(values[name]) ** exponent
            total += term
        return total

    def positive_coefficients(self):
        return bool(self.terms) and all(c > 0 for c in self.terms.values())

    def digest(self):
        data = [[list(p), str(c)] for p, c in sorted(self.terms.items())]
        return hashlib.sha256(json.dumps(data, separators=(",", ":")).encode()).hexdigest()


def variable(index):
    powers = list(ZERO_EXP)
    powers[index] = 1
    return Poly({tuple(powers): 1})


a, b, e, f, u, w, v, z, lam = [variable(i) for i in range(len(NAMES))]


def matrix(entries):
    return [[Poly(x) for x in row] for row in entries]


def multiply(left, right):
    return [[sum((left[i][k] * right[k][j] for k in range(len(right))), Poly())
             for j in range(len(right[0]))] for i in range(len(left))]


def determinant(mat):
    """Direct Leibniz formula: independent of CAS determinant algorithms."""
    n = len(mat)
    result = Poly()
    for perm in permutations(range(n)):
        inversions = sum(perm[i] > perm[j] for i in range(n) for j in range(i + 1, n))
        term = Poly(-1 if inversions % 2 else 1)
        for i, j in enumerate(perm):
            term *= mat[i][j]
            if not term.terms:
                break
        result += term
    return result


def charpoly(mat):
    return determinant([[lam - mat[i][j] if i == j else -mat[i][j]
                         for j in range(len(mat))] for i in range(len(mat))])


def check(condition, message):
    # Deliberately not Python assert: these checks also execute under python -O.
    if not condition:
        raise ValueError(message)


def verify(corruption=None):
    h = a + e + u + w
    j = 2*b + f + v + z
    k = f + v + z
    G = matrix([
        [0, 0, -w, 0, z, 0], [0, 0, 0, z, 0, -w],
        [a, 0, -h, -e, 0, 0], [-b, -b, -f, -b-f-v-z, -b, 0],
        [-b, -b, 0, -b, -b-f-v-z, -f], [0, a, 0, 0, -e, -h],
    ])
    A = matrix([[0, -w, z], [a, -h, -e], [-2*b, -f, -j]])
    B = matrix([[0, -w, -z], [a, -h, -e], [0, -f, -k]])
    # Columns: Q1+Q2, U+X, V+Y, Q1-Q2, U-X, V-Y.
    P = matrix([[1,0,0,1,0,0], [1,0,0,-1,0,0],
                [0,1,0,0,1,0], [0,0,1,0,0,1],
                [0,0,1,0,0,-1], [0,1,0,0,-1,0]])
    if corruption == "matrix_entry":
        G[0][2] = w
    if corruption == "parity_block":
        B[0][2] = z
    D = [[(A[i][jj] if i < 3 and jj < 3 else
           B[i-3][jj-3] if i >= 3 and jj >= 3 else Poly())
          for jj in range(6)] for i in range(6)]
    check(determinant(P) != 0, "singular parity change of basis")
    check(multiply(G, P) == multiply(P, D), "G P = P diag(A,B) failed")
    pA, pB = charpoly(A), charpoly(B)
    check(charpoly(G) == pA*pB, "six-dimensional characteristic factorization failed")
    cA = [pA.coefficient("lambda", i) for i in (3,2,1,0)]
    cB = [pB.coefficient("lambda", i) for i in (3,2,1,0)]
    delta = w*k - f*z
    if corruption == "threshold_sign":
        delta = w*k + f*z
    expectedA = [Poly(1), h+j, h*j-e*f+a*w+2*b*z,
                 a*w*j+2*b*e*w+a*f*z+2*b*z*h]
    expectedB = [Poly(1), h+k, h*k-e*f+a*w, a*delta]
    check(cA == expectedA and cB == expectedB, "cubic coefficient identities failed")
    check(cA[1].positive_coefficients() and cA[2].positive_coefficients()
          and cA[3].positive_coefficients(), "even coefficient positivity failed")
    check(cB[1].positive_coefficients() and cB[2].positive_coefficients(),
          "odd c1,c2 positivity failed")
    HA = cA[1]*cA[2]-cA[3]
    HB = cB[1]*cB[2]-cB[3]
    check(HA.positive_coefficients() and HB.positive_coefficients(),
          "cubic Routh coefficient positivity failed")
    # Independently verify a compact positive decomposition supplied during review.
    alpha = a+u+w
    beta = 2*b+v+z
    grouped_even = (alpha**2*j + e*alpha*j + e**2*beta + e*j*beta
                    + a*w*h + 2*b*z*j
                    + e*((a+u)*beta+w*(v+z))
                    + (u+w)*j**2 + a*(j**2-f*z))
    grouped_odd = (h+k)*(h*k-e*f)+a*w*h+a*f*z
    check(HA == grouped_even and HB == grouped_odd, "grouped Routh identities failed")
    check((h*k-e*f).positive_coefficients() and (j**2-f*z).positive_coefficients(),
          "positive group premise failed")
    # These examples cover all three exact threshold signs, including a rational boundary.
    examples = []
    for expected_sign, values in [
        (1, dict(a=1,b=1,e=1,f=1,u=1,w=1,v=1,z=1)),
        (-1, dict(a=1,b=1,e=1,f=4,u=1,w=1,v=1,z=4)),
        (0, dict(a=1,b=1,e=1,f=4,u=1,w=Fraction(16,9),v=1,z=4)),
    ]:
        val = delta.evaluate(values)
        sign = (val > 0) - (val < 0)
        check(sign == expected_sign, "threshold fixture failed")
        # det(G)=cA3*cB3 since both factors are degree three.
        detvalue = determinant(G).evaluate(values)
        check(((detvalue > 0)-(detvalue < 0)) == expected_sign,
              "independent six-dimensional determinant sign failed")
        examples.append({"parameters": {x: str(y) for x,y in values.items()},
                         "delta": str(val), "detG": str(detvalue)})
    return {
        "source": "https://arxiv.org/html/2608.27081v2",
        "domain": "a,b,e,f,u,w,v,z strictly positive; paired-rate Z2-equivariant slice",
        "checks": ["invertible parity basis", "GP=PD", "six-dimensional charpoly factorization",
                   "cubic coefficient identities", "strict positive-coefficient Routh polynomials",
                   "compact positive-group identities", "three exact threshold regimes"],
        "basis_determinant": str(determinant(P).evaluate({})),
        "even_Routh": {"monomials": len(HA.terms), "min_coefficient": str(min(HA.terms.values())), "sha256": HA.digest()},
        "odd_Routh": {"monomials": len(HB.terms), "min_coefficient": str(min(HB.terms.values())), "sha256": HB.digest()},
        "examples": examples,
        "limitation": "Algebraic certificate for a symmetry subcase; no novelty or unrestricted-conjecture certification.",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true", help="reject three deliberately corrupted claims")
    args = parser.parse_args()
    receipt = verify()
    if args.self_test:
        rejected = []
        for corruption in ("matrix_entry", "parity_block", "threshold_sign"):
            try:
                verify(corruption)
            except ValueError:
                rejected.append(corruption)
            else:
                raise SystemExit("negative control was incorrectly accepted: " + corruption)
        receipt["rejected_negative_controls"] = rejected
    print(json.dumps(receipt, ensure_ascii=False, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
