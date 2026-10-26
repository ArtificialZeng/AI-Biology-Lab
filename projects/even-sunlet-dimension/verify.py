"""Exact tropical exponent matrix for the source 3-sunlet map; stdout only."""
import flint
from itertools import product
from math import gcd, lcm
import json
from flint import fmpq_mat

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def finite_group(moduli):
    elements = tuple(product(*(range(m) for m in moduli)))
    index = {g: i for i, g in enumerate(elements)}
    zero = (0,) * len(moduli)

    def add(g, h):
        return tuple(((a + b) % m for a, b, m in zip(g, h, moduli)))

    def neg(g):
        return tuple((-a % m for a, m in zip(g, moduli)))

    def order(g):
        return lcm(*(m // gcd(a, m) for a, m in zip(g, moduli)))
    return (elements, index, zero, add, neg, order)

def matrix(moduli):
    elements, index, zero, add, neg, order = finite_group(moduli)
    n = len(elements)
    generator = next((g for g in elements if order(g) == n), None)
    if generator is not None:
        X = set()
        g = zero
        for _ in range(n // 2):
            g = add(g, generator)
            X.add(g)
    else:
        X = {g for g in elements if g != zero and (g == neg(g) or g < neg(g))}
    Y = set(elements) - X - {zero}
    need(all(((g in X) != (neg(g) in X) for g in elements if g != neg(g))), 'inverse-half')
    need(all((g in X for g in elements if g != zero and g == neg(g))), 'involutions in X')
    nu = {g: 4 if g == zero else 0 for g in elements}
    mu = {g: -1 if g == zero else 2 if g in X else -2 for g in elements}
    weights = [{g: 1 for g in elements} for _ in range(3)]
    weights += [{g: nu[g] + 1 for g in elements}, {g: 3 for g in elements}, {g: 3 + mu[g] for g in elements}]
    need(all((value > 0 for block in weights for value in block.values())), 'positive actual lambda')
    M = [[0] * (n * n) for _ in range(6 * n)]
    choices = []
    margins = []
    for i, g in enumerate(elements):
        for j, h in enumerate(elements):
            gh = add(g, h)
            leaf3 = neg(gh)
            delta = nu[gh] - nu[h] - mu[g]
            common = weights[0][g] + weights[1][h] + weights[2][leaf3]
            weight_first = common + weights[3][gh] + weights[4][g]
            weight_second = common + weights[3][h] + weights[5][g]
            need(weight_first - weight_second == delta, 'literal monomial weight comparison')
            need(delta != 0, 'strict tropical chamber')
            branch = 1 if delta < 0 else 2
            choices.append(branch)
            margins.append(abs(delta))
            expected = 2 if g == zero else (2 if h == neg(g) else 1) if g in X else 1 if h == zero else 2
            need(branch == expected, 'casewise monomial selection including involutions')
            labels = [(0, g), (1, h), (2, leaf3), (3, gh if branch == 1 else h), (4 if branch == 1 else 5, g)]
            for edge, label in labels:
                M[edge * n + index[label]][i * n + j] += 1
    return (M, (elements, index, zero, add, neg, order), X, Y, choices, min(margins))

def kernel_basis(info):
    elements, index, zero, add, neg, order = info
    n = len(elements)
    z = index[zero]
    out = []
    for i, g in enumerate(elements):
        v = [0] * (6 * n)
        v[i] = 1
        v[5 * n + i] = -1
        if g != zero:
            v[4 * n + i] = -1
        out.append(v)
    for edge in [1, 2, 3]:
        v = [0] * (6 * n)
        for i, g in enumerate(elements):
            v[edge * n + i] = 1
            v[5 * n + i] = -1
            if g != zero:
                v[4 * n + i] = -1
        out.append(v)
    v = [0] * (6 * n)
    v[4 * n + z] = 1
    out.append(v)
    return out
CASES = [(5,), (6,), (8,), (10,), (12,), (14,), (18,), (20,), (2, 4), (2, 2, 3), (3, 3), (2, 3, 3), (4, 4), (2, 2, 4), (2, 2, 2, 4), (2, 2), (2, 2, 2), (2, 2, 2, 2)]

def build_certificate():
    need(flint.__version__ == '0.9.0', 'python-flint 0.9.0 required')
    results = []
    controls = {(2, 2): 13, (2, 2, 2): 29, (2, 2, 2, 2): 61}
    for mods in CASES:
        rows, info, X, Y, choices, margin = matrix(mods)
        elements, index, zero, add, neg, order = info
        n = len(elements)
        A = fmpq_mat(rows)
        B = fmpq_mat(kernel_basis(info))
        rank = A.rank()
        brank = B.rank()
        need(B * A == fmpq_mat(n + 4, n * n), 'explicit kernel identities')
        need(brank == n + 4, 'kernel parameter independence')
        exponent = lcm(*mods)
        target = 5 * n - 4
        diagnostic = n >= 5 and exponent > 2
        member = diagnostic and n % 2 == 0
        if diagnostic:
            need(rank == target, 'target class or known odd control diagnostic rank')
        if mods in controls:
            need(rank == controls[mods], 'known source elementary2 negative control')
        results.append({'group': list(mods), 'order': n, 'exponent': exponent, 'cyclic': any((order(g) == n for g in elements)), 'rank_QQ': rank, 'parameter_upper_target': target, 'target_class': member, 'kernel_basis_rank': brank, 'strict_margin': margin, 'first_monomial_columns': choices.count(1), 'second_monomial_columns': choices.count(2), 'X': [list(g) for g in sorted(X)], 'Y': [list(g) for g in sorted(Y)]})
    return {'ok': True, 'method': 'exact rational ranks and explicit integer left-kernel vectors; global theorem supplied separately', 'finite_named_groups': len(CASES), 'cases': results, 'floating_point_used': False, 'finite_tests_not_infinite_proof': True, 'elementary2_not_resolved_by_this_argument': True}

from archive_tools import run
if __name__ == "__main__":
    raise SystemExit(run(build_certificate,"evidence/exact-certificate.json","exact tropical ranks and integer kernel basis"))
