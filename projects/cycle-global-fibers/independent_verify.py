#!/usr/bin/env python3
"""Independent standard-library exact verifier; never imports the candidate checker.

Run: python3 independent_checker.py --output verification-receipt.json
All algebra is in Z[s,k_1,...,k_n,e_r,e_t] or fractions.Fraction.
The infinite completeness and reducedness arguments are in PROOF.md.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as Q
from itertools import permutations
from math import factorial, isqrt, prod
from pathlib import Path
import hashlib
import json


def require(condition, message="independent exact verification failed"):
    if not condition:
        raise ValueError(message)


# Sparse multivariate polynomial ring. Variable 0 is s; the next n variables
# are edge rates, followed by the two independently labelled leak rates.
def const(c, v):
    return {} if c == 0 else {(0,) * v: c}


def variable(i, v):
    a = [0] * v
    a[i] = 1
    return {tuple(a): 1}


def add(a, b, sign=1):
    c = dict(a)
    for m, x in b.items():
        c[m] = c.get(m, 0) + sign * x
        if not c[m]:
            del c[m]
    return c


def mul(a, b):
    c = defaultdict(int)
    for ma, ca in a.items():
        for mb, cb in b.items():
            c[tuple(x + y for x, y in zip(ma, mb))] += ca * cb
    return {m: x for m, x in c.items() if x}


def times(polys, v):
    a = const(1, v)
    for b in polys:
        a = mul(a, b)
    return a


def det_sparse(matrix, v):
    """Generic Leibniz determinant, pruning zero entries only; no cycle formula."""
    n = len(matrix)
    total = {}
    nonzero_terms = 0

    def visit(row, cols, used, term):
        nonlocal total, nonzero_terms
        if row == n:
            inversions = sum(cols[i] > cols[j] for i in range(n) for j in range(i + 1, n))
            total = add(total, term, -1 if inversions % 2 else 1)
            nonzero_terms += 1
            return
        for col, entry in enumerate(matrix[row]):
            if entry and col not in used:
                visit(row + 1, cols + [col], used | {col}, mul(term, entry))

    visit(0, [], set(), const(1, v))
    return total, nonzero_terms


def ode_polynomials(n, p, r, t):
    """Build sI-A from dx_i=-d_i*x_i+k_(i-1)*x_(i-1)+delta_i1*u."""
    v = n + 3
    s = variable(0, v)
    edges = [variable(i + 1, v) for i in range(n)]
    ds = list(edges)
    ds[r - 1] = add(ds[r - 1], variable(n + 1, v))
    ds[t - 1] = add(ds[t - 1], variable(n + 2, v))
    matrix = [[{} for _ in range(n)] for _ in range(n)]
    for i in range(n):
        matrix[i][i] = add(s, ds[i])
        matrix[(i + 1) % n][i] = add({}, edges[i], -1)
    D, terms_d = det_sparse(matrix, v)
    minor = [[matrix[i][j] for j in range(n) if j != p - 1]
             for i in range(n) if i != 0]
    N, terms_n = det_sparse(minor, v)
    if (p + 1) % 2:
        N = add({}, N, -1)
    expected_d = add(times([add(s, d) for d in ds], v), times(edges, v), -1)
    expected_n = times(edges[:p - 1] + [add(s, d) for d in ds[p:]], v)
    require(D == expected_d and N == expected_n)
    require(terms_d == 2 and terms_n == 1)
    return D, N


def evaluate_s(poly, pars):
    degree = max((m[0] for m in poly), default=0)
    out = [Q(0)] * (degree + 1)
    for monomial, coefficient in poly.items():
        value = Q(coefficient)
        for x, e in zip(pars, monomial[1:]):
            if e:
                value *= x ** e
        out[monomial[0]] += value
    return trim(out)


def trim(a):
    a = list(a)
    while len(a) > 1 and not a[-1]:
        a.pop()
    return tuple(a)


def divrem(a, b):
    a = list(trim(a))
    b = trim(b)
    q = [Q(0)] * max(1, len(a) - len(b) + 1)
    while len(a) >= len(b) and any(a):
        shift = len(a) - len(b)
        c = a[-1] / b[-1]
        q[shift] += c
        for i, x in enumerate(b):
            a[i + shift] -= c * x
        a = list(trim(a))
    return trim(q), trim(a)


def root_values(poly):
    """Factor the actual recovered polynomial, without reading planted diagonals.

    Fixtures have integer negated roots. Exact division to degree zero certifies
    that no other complex root is omitted in these particular fixtures.
    """
    f = trim(poly)
    require(f[-1] == 1 and all(x.denominator == 1 for x in f))
    c = abs(int(f[0]))
    divisors = set()
    for x in range(1, isqrt(c) + 1):
        if c % x == 0:
            divisors.update((x, c // x))
    roots = []
    for alpha in sorted(divisors):
        while len(f) > 1:
            quotient, remainder = divrem(f, (Q(alpha), Q(1)))
            if any(remainder):
                break
            roots.append(Q(alpha))
            f = quotient
    require(f == (Q(1),), "fixture is not completely integer-split")
    require(len(roots) == len(set(roots)), "fixture has a repeated within-block root")
    return roots


def inverse_from_data(D, N, n, p, r, t):
    kappa = N[-1]
    T = tuple(x / kappa for x in N)
    U, remainder = divrem(D, T)
    require(len(remainder) == 1)
    C = -remainder[0]
    require(C and kappa)
    front, back = root_values(U), root_values(T)
    require(len(front) == p and len(back) == n - p)
    for first in permutations(front):
        for last in permutations(back):
            ds = first + last
            edges = list(ds)
            # Solve the independent path-product and complementary-product
            # constraints, each linear in its one unknown leaking edge.
            coefficient_r = prod(ds[i] for i in range(p - 1) if i != r - 1)
            coefficient_t = prod(ds[i] for i in range(p - 1, n) if i != t - 1)
            edges[r - 1] = kappa / coefficient_r
            edges[t - 1] = (C / kappa) / coefficient_t
            yield tuple(edges) + (ds[r - 1] - edges[r - 1], ds[t - 1] - edges[t - 1])


def rank_exact(rows):
    a = [list(row) for row in rows]
    rank = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(rank, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        d = a[rank][col]
        a[rank] = [x / d for x in a[rank]]
        for i in range(rank + 1, len(a)):
            d = a[i][col]
            if d:
                a[i] = [x - d * y for x, y in zip(a[i], a[rank])]
        rank += 1
        if rank == len(a):
            break
    return rank


def coefficient_jacobian(D, N, pars, n, p):
    rows = [[Q(0)] * (n + 2) for _ in range(n + n - p + 1)]
    for poly, offset, cap in [(D, 0, n - 1), (N, n, n - p)]:
        for monomial, coefficient in poly.items():
            power_s = monomial[0]
            if power_s > cap:
                continue
            for j, power_j in enumerate(monomial[1:]):
                if not power_j:
                    continue
                value = Q(coefficient * power_j)
                for k, (x, exponent) in enumerate(zip(pars, monomial[1:])):
                    exponent -= (k == j)
                    if exponent:
                        value *= x ** exponent
                rows[offset + power_s][j] += value
    return rows


def seed(n, p, r, t, m):
    front = [Q(5 + 3 * j) for j in range(p)]
    back = [Q(101 + 2 * j) for j in range(n - p)]
    ds = list(reversed(front[1:])) + [front[0]] + back
    a = Q(4)
    b = front[m - 1] + 1
    kappa = prod(front) / b
    C = kappa * prod(back) * a
    edges = list(ds)
    edges[r - 1] = kappa / prod(ds[i] for i in range(p - 1) if i != r - 1)
    edges[t - 1] = C / (kappa * prod(ds[i] for i in range(p - 1, n) if i != t - 1))
    pars = tuple(edges) + (ds[r - 1] - edges[r - 1], ds[t - 1] - edges[t - 1])
    require(all(x > 0 for x in pars))
    return pars


def check_fiber(Dsym, Nsym, pars, n, p, r, t, expected_m=None):
    D, N = evaluate_s(Dsym, pars), evaluate_s(Nsym, pars)
    points = list(inverse_from_data(D, N, n, p, r, t))
    require(len(set(points)) == len(points) == factorial(p) * factorial(n - p))
    kappa = N[-1]
    T = tuple(x / kappa for x in N)
    U, remainder = divrem(D, T)
    C = -remainder[0]
    a, b = C / (kappa * T[0]), U[0] / kappa
    m = sum(a < alpha < b for alpha in root_values(U))
    positive = sum(all(x > 0 for x in point) for point in points)
    require(positive == m * factorial(p - 1) * factorial(n - p))
    if expected_m is not None:
        require(m == expected_m)
    for point in points:
        require(all(x > 0 for x in point[:n]))
        require(evaluate_s(Dsym, point) == D)
        require(evaluate_s(Nsym, point) == N)
    require(rank_exact(coefficient_jacobian(Dsym, Nsym, pars, n, p)) == n + 2)
    return {"n": n, "p": p, "leaks": [r, t], "m": m,
            "complex_affine_count": len(points), "positive_count": positive,
            "all_parameters_nonzero_count": sum(all(x != 0 for x in point) for point in points),
            "interval": [str(a), str(b)], "D_ascending": list(map(str, D)),
            "N_ascending": list(map(str, N))}, points


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-n", type=int, default=6)
    ap.add_argument("--symbolic-max-n", type=int, default=8)
    ap.add_argument("--output", type=Path, help="write a separate receipt; default stdout")
    args = ap.parse_args()
    require(args.max_n >= 3 and args.symbolic_max_n >= 3, "both finite ceilings must be at least 3")
    symbolic_cases, fibers, branch_count = [], [], 0
    for n in range(3, max(args.max_n, args.symbolic_max_n) + 1):
        for p in range(2, n):
            for r in range(1, p):
                for t in range(p, n + 1):
                    Ds, Ns = ode_polynomials(n, p, r, t)
                    symbolic_cases.append([n, p, r, t])
                    if n <= args.max_n:
                        for m in range(1, p + 1):
                            result, points = check_fiber(Ds, Ns, seed(n, p, r, t, m), n, p, r, t, m)
                            fibers.append(result)
                            branch_count += len(points)

    special = []
    Ds, Ns = ode_polynomials(3, 2, 1, 2)
    for tag, pars in [
        ("upper_endpoint_zero_leak_in_other_preimage", (Q(2), Q(3, 2), Q(4), Q(1), Q(1, 2))),
        ("lower_endpoint_zero_leak_in_other_preimage", (Q(3, 2), Q(2), Q(4), Q(1, 2), Q(1))),
    ]:
        record, points = check_fiber(Ds, Ns, pars, 3, 2, 1, 2, 1)
        require(record["complex_affine_count"] == 2 and record["all_parameters_nonzero_count"] == 1)
        record["tag"] = tag
        record["all_preimages"] = [list(map(str, x)) for x in points]
        special.append(record)
    # Equal diagonal values across the two recovered factors are allowed.
    Ds, Ns = ode_polynomials(4, 2, 1, 3)
    record, _ = check_fiber(Ds, Ns, (Q(12, 7), Q(2), Q(9, 4), Q(5), Q(9, 7), Q(3, 4)), 4, 2, 1, 3, 2)
    record["tag"] = "cross_block_coincident_root"
    special.append(record)
    receipt = {"status": "PASS", "arithmetic": "integer symbolic polynomials and exact rational arithmetic only",
               "imports_candidate_checker": False,
               "symbolic_ode_cases": len(symbolic_cases), "symbolic_case_scope": symbolic_cases,
               "rational_fixture_count": len(fibers), "rational_branch_count": branch_count,
               "rational_scope": {"n": [3, args.max_n], "placements": "all allowed", "m": "all 1..p"},
               "exact_full_column_rank": "n+2 on every rational source fixture and all special cases",
               "symbolic_ode_scope": {"n": [3, max(args.max_n, args.symbolic_max_n)], "placements": "all allowed"},
               "special_cases": special, "fixtures": fibers,
               "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               "limitations": "Finite checks do not prove arbitrary-n completeness or novelty. This is not a runner release-gate receipt."}
    rendered = json.dumps(receipt, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end="")



if __name__ == "__main__":
    main()
