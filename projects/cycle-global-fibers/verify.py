#!/usr/bin/env python3
"""Exact finite checks for the separately proved infinite two-leak fiber theorem.

Default: every placement, every attainable m, 3<=n<=7; rational arithmetic.
Independent symbolic matrix/cofactor verification: every placement, 3<=n<=5.
The finite enumeration checks the implementation, not the infinite proof.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as F
from itertools import permutations
import json
from math import factorial, prod
from pathlib import Path


def check(condition, message="exact verification failed"):
    if not condition:
        raise ValueError(message)


def polynomial(ds):
    """Ascending coefficients of product(s+d)."""
    coeff = [F(1)]
    for d in ds:
        nxt = [F(0)] * (len(coeff) + 1)
        for j, c in enumerate(coeff):
            nxt[j] += c * d
            nxt[j + 1] += c
        coeff = nxt
    return tuple(coeff)


def divide(a, b):
    a = list(a)
    quotient = [F(0)] * max(1, len(a) - len(b) + 1)
    while len(a) >= len(b):
        degree = len(a) - len(b)
        c = a[-1] / b[-1]
        quotient[degree] = c
        for j, v in enumerate(b):
            a[degree + j] -= c * v
        while a and a[-1] == 0:
            a.pop()
    return tuple(quotient), tuple(a)


def signature(k, leaks, p):
    d = [k[i] + leaks.get(i, F(0)) for i in range(len(k))]
    C = prod(k)
    kappa = prod(k[: p - 1])
    D = list(polynomial(d))
    D[0] -= C
    N = tuple(kappa * c for c in polynomial(d[p:]))
    return tuple(D), N


def inverse(d, r, t, p, kappa, C):
    """r,t are 0-based, p is the 1-based output label."""
    k = list(d)
    k[r] = kappa / prod(d[i] for i in range(p - 1) if i != r)
    k[t] = C / (kappa * prod(d[i] for i in range(p - 1, len(d)) if i != t))
    return tuple(k), {r: d[r] - k[r], t: d[t] - k[t]}


def fixture(n, p, r, t, m):
    alpha = tuple(F(i) for i in range(2, p + 2))
    beta = tuple(F(i) for i in range(p + 2, n + 2))
    A, B = prod(alpha), prod(beta)
    lower, upper = F(3, 2), F(2 * m + 3, 2)
    kappa = A / upper
    C = kappa * B * lower
    # Output diagonal = 2, guaranteeing a strictly positive source.
    d = alpha[1:] + alpha[:1] + beta
    k, leaks = inverse(d, r, t, p, kappa, C)
    check(all(v > 0 for v in k))
    check(all(v > 0 for v in leaks.values()))
    return alpha, beta, kappa, C, k, leaks


def validate_preimage(k, leaks, p, expected):
    check(signature(k, leaks, p) == expected, "branch fails coefficient reconstruction")


def check_fixture(n, p, r, t, m, keep=False):
    alpha, beta, kappa, C, original_k, original_leaks = fixture(n, p, r, t, m)
    D, N = expected = signature(original_k, original_leaks, p)
    check(N[-1] == kappa)
    T = tuple(c / kappa for c in N)
    U, remainder = divide(D, T)
    check(remainder == (-C,), remainder)
    # Adding C only changes the constant coefficient of the quotient numerator.
    check(U == polynomial(alpha))
    check(T == polynomial(beta))
    a, b = C / (kappa * T[0]), U[0] / kappa
    check(sum(a < v < b for v in alpha) == m)
    seen = set()
    positive = 0
    root_counts = {str(v): 0 for v in alpha}
    saved = []
    for front in permutations(alpha):
        for back in permutations(beta):
            d = front + back
            k, leaks = inverse(d, r, t, p, kappa, C)
            vector = k + (leaks[r], leaks[t])
            check(vector not in seen)
            seen.add(vector)
            validate_preimage(k, leaks, p, expected)
            check(all(v > 0 for v in k))
            check(k[r] / d[r] == front[-1] / b)
            check(k[t] / d[t] == a / front[-1])
            is_positive = all(v > 0 for v in leaks.values())
            check(is_positive == (a < front[-1] < b))
            positive += is_positive
            if is_positive:
                root_counts[str(front[-1])] += 1
            if keep:
                saved.append({
                    "edge_rates": [str(v) for v in k],
                    "leak_rates": {str(r + 1): str(leaks[r]), str(t + 1): str(leaks[t])},
                    "output_diagonal": str(front[-1]),
                    "strictly_positive": is_positive,
                })
    algebraic = factorial(p) * factorial(n - p)
    predicted_positive = m * factorial(p - 1) * factorial(n - p)
    check(len(seen) == algebraic)
    check(positive == predicted_positive)
    record = {
        "n": n, "p": p, "leaks": [r + 1, t + 1], "m": m,
        "complex_fiber_count": len(seen), "positive_fiber_count": positive,
        "interval": [str(a), str(b)], "eligible_root_branch_counts": root_counts,
    }
    if keep:
        record.update({
            "D_ascending": [str(v) for v in D],
            "N_ascending": [str(v) for v in N],
            "kappa": str(kappa), "C": str(C),
            "U_roots_negated": [str(v) for v in alpha],
            "T_roots_negated": [str(v) for v in beta],
            "all_preimages": saved,
        })
    return record


def negative_controls():
    """Reject two incorrect edge-product branches and a wrong root partition."""
    n, p, r, t = 4, 2, 0, 2
    alpha, beta, kappa, C, k, leaks = fixture(n, p, r, t, 1)
    expected = signature(k, leaks, p)
    bad = []
    for label, index in (("prefix_product", r), ("downstream_product", t)):
        wrong_k, wrong_leaks = list(k), dict(leaks)
        wrong_k[index] += k[index]
        wrong_leaks[index] -= k[index]  # preserve every diagonal, break product data
        bad.append((label, tuple(wrong_k), wrong_leaks))
    ds = list(alpha[1:] + alpha[:1] + beta)
    ds[0], ds[p] = ds[p], ds[0]
    wrong_k, wrong_leaks = inverse(tuple(ds), r, t, p, kappa, C)
    bad.append(("root_partition", wrong_k, wrong_leaks))
    rejected = []
    for label, wrong_k, wrong_leaks in bad:
        try:
            validate_preimage(wrong_k, wrong_leaks, p, expected)
        except ValueError:
            rejected.append(label)
        else:
            raise ValueError("invalid inverse branch was accepted: " + label)
    return rejected


def independent_symbolic_matrices(max_n=5):
    import sympy as sp

    s = sp.Symbol("s")
    checked = []
    for n in range(3, max_n + 1):
        for p in range(2, n):
            for r in range(p - 1):
                for t in range(p - 1, n):
                    k = sp.symbols(f"k1:{n+1}")
                    er, et = sp.symbols("er et")
                    ds = list(k)
                    ds[r] += er
                    ds[t] += et
                    M = sp.zeros(n)
                    for i in range(n):
                        M[i, i] = -ds[i]
                        M[(i + 1) % n, i] = k[i]
                    observed_D = M.charpoly(s).as_expr()
                    B = s * sp.eye(n) - M
                    observed_N = (-1) ** (p - 1) * B.minor_submatrix(0, p - 1).det()
                    predicted_D = sp.prod(s + v for v in ds) - sp.prod(k)
                    predicted_N = sp.prod(k[: p - 1]) * sp.prod(s + v for v in ds[p:])
                    check(sp.expand(observed_D - predicted_D) == 0)
                    check(sp.expand(observed_N - predicted_N) == 0)
                    checked.append({"n": n, "p": p, "leaks": [r + 1, t + 1]})
    return checked


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=7)
    parser.add_argument("--self-test", action="store_true", help="reject three invalid inverse branches")
    parser.add_argument("--output", type=Path, help="write a separate receipt; default stdout")
    args = parser.parse_args()
    check(args.max_n >= 3, "max-n must be at least 3")
    symbolic = independent_symbolic_matrices()
    records = []
    for n in range(3, args.max_n + 1):
        for p in range(2, n):
            for r in range(p - 1):
                for t in range(p - 1, n):
                    for m in range(1, p + 1):
                        records.append(check_fixture(n, p, r, t, m))
    examples = [check_fixture(4, 3, 0, t, m, True) for t in (2, 3) for m in (1, 2, 3)]
    examples += [check_fixture(3, 2, 0, 1, m, True) for m in (1, 2)]
    result = {
        "status": "PASS", "arithmetic": "exact rational / symbolic integer polynomials",
        "finite_scope": {"min_n": 3, "max_n": args.max_n,
                         "placements": "every 1<=r<p<=t<=n, 2<=p<n",
                         "m": "every 1<=m<=p"},
        "symbolic_matrix_count": len(symbolic),
        "symbolic_matrix_scope": symbolic,
        "rational_fixture_count": len(records),
        "preimage_vectors_checked": sum(row["complex_fiber_count"] for row in records),
        "fixtures": records, "exact_examples": examples,
        "limitations": "Completeness for arbitrary n follows from PROOF.md, not finite tests. No novelty or independent peer-review claim.",
    }
    if args.self_test:
        result["rejected_negative_controls"] = negative_controls()
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end="")



if __name__ == "__main__":
    main()
