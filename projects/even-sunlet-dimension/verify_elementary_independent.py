"""Independent, standalone, standard-library checker for the frozen R12 claim.

No owner code or private path is imported. This is finite exact evidence, not
the uniform proof. Run with Python 3.10+; stdout is deterministic JSON.
"""
from __future__ import annotations
import hashlib
import json
from math import gcd, isqrt
from functools import reduce
PRIME = 1000000007

def rank_mod(vectors, prime=PRIME):
    """Exact Gaussian elimination over F_prime; integer Q-rank lower bound."""
    pivots = {}
    for source in vectors:
        row = {i: a % prime for i, a in source.items() if a % prime}
        while row:
            col = min(row)
            if col not in pivots:
                scale = pow(row[col], -1, prime)
                pivots[col] = {i: a * scale % prime for i, a in row.items()}
                break
            factor = row[col]
            for i, a in pivots[col].items():
                value = (row.get(i, 0) - factor * a) % prime
                if value:
                    row[i] = value
                else:
                    row.pop(i, None)
    return len(pivots)

def rank_q(vectors):
    """Exact integer row elimination with primitive gcd normalization."""
    pivots = {}
    for source in vectors:
        row = {i: a for i, a in source.items() if a}
        while row:
            col = min(row)
            if col not in pivots:
                common = reduce(gcd, (abs(a) for a in row.values()))
                sign = 1 if row[col] > 0 else -1
                pivots[col] = {i: sign * a // common for i, a in row.items()}
                break
            pivot = pivots[col]
            a, b = (row[col], pivot[col])
            new = {i: b * row.get(i, 0) - a * pivot.get(i, 0) for i in row.keys() | pivot.keys()}
            row = {i: x for i, x in new.items() if x}
            if row:
                common = reduce(gcd, (abs(x) for x in row.values()))
                row = {i: x // common for i, x in row.items()}
    return len(pivots)

def literal_monomials(g, h, n):
    """Five edge-coordinate factors from each literal source summand."""
    k = g ^ h
    first = {g: 1, n + h: 1, 2 * n + k: 1, 3 * n + k: 1, 4 * n + g: 1}
    second = {g: 1, n + h: 1, 2 * n + k: 1, 3 * n + h: 1, 5 * n + g: 1}
    return (first, second)

def explicit_kernel(n):
    """The n a-directions, B/C/D directions, and unused e_0 direction."""
    vectors = []
    for g in range(n):
        v = {g: 1, 5 * n + g: -1}
        if g:
            v[4 * n + g] = -1
        vectors.append(v)
    for block in (1, 2, 3):
        v = {block * n + h: 1 for h in range(n)}
        v.update({5 * n + h: -1 for h in range(n)})
        v.update({4 * n + h: -1 for h in range(1, n)})
        vectors.append(v)
    vectors.append({4 * n: 1})
    return vectors

def scalar_constraints(r, hamming=True):
    """Eliminate t_g by comparing all oriented pairs with the same XOR."""
    n = 1 << r
    order = sorted(range(n), key=lambda x: (x.bit_count(), x)) if hamming else list(range(n))
    place = {x: i for i, x in enumerate(order)}
    equations = []
    for g in range(1, n):
        x0, y0 = (0, g) if place[0] < place[g] else (g, 0)
        for x in range(n):
            y = x ^ g
            if place[x] >= place[y]:
                continue
            row = {}
            for j, c in ((y, 1), (x, -1), (y0, -1), (x0, 1)):
                row[j] = row.get(j, 0) + c
            equations.append({j: c for j, c in row.items() if c})
    return equations

def one_case(r):
    n = 1 << r
    nu = [4 * (n * x.bit_count() + x) for x in range(n)]
    mu = [-1] + [1] * (n - 1)
    weights = [1] * (3 * n) + [v + 1 for v in nu] + [2] * n + [2 + v for v in mu]
    columns, decisions = ([], {})
    margins, direct_reduction_ok = ([], True)
    signature = hashlib.sha256()
    first_count = 0
    for g in range(n):
        for h in range(n):
            first, second = literal_monomials(g, h, n)
            w1 = sum((weights[i] * a for i, a in first.items()))
            w2 = sum((weights[i] * a for i, a in second.items()))
            if w1 == w2:
                raise AssertionError('tie')
            choose_first = w1 < w2
            decisions[g, h] = choose_first
            first_count += choose_first
            columns.append(first if choose_first else second)
            margins.append(abs(w1 - w2))
            direct_reduction_ok &= w1 - w2 == nu[g ^ h] - nu[h] - mu[g]
            signature.update(json.dumps([g, h, w1, w2, sorted(columns[-1].items())], separators=(',', ':')).encode())
    require(direct_reduction_ok, 'Independent invariant, original line 129')
    require(all((decisions[g, h] != decisions[g, h ^ g] for g in range(1, n) for h in range(n))), 'Independent invariant, original line 130')
    require(not any((decisions[0, h] for h in range(n))), 'Independent invariant, original line 132')
    basis = explicit_kernel(n)
    annihilates = all((sum((a * v.get(i, 0) for i, a in col.items())) == 0 for col in columns for v in basis))
    require(annihilates, 'Independent invariant, original line 136')
    lower = rank_mod(columns)
    basis_lower = rank_mod(basis)
    require(basis_lower == n + 4, 'Independent invariant, original line 139')
    upper = 6 * n - basis_lower
    qrank = rank_q(columns) if r <= 3 else None
    require(qrank is None or qrank >= lower, 'Independent invariant, original line 142')
    if r >= 3:
        require(lower == upper == 5 * n - 4, 'Independent invariant, original line 144')
    scalar = scalar_constraints(r)
    scalar_rank = rank_q(scalar) if r <= 3 else rank_mod(scalar)
    integer_scalar_rank = rank_mod(scalar_constraints(r, False))
    return {'r': r, 'n': n, 'matrix_shape': [6 * n, n * n], 'held_out_beyond_owner_cases_r_le_6': r == 7, 'positive_weights': min(weights) > 0, 'strict_min_margin': min(margins), 'literal_weight_identity': direct_reduction_ok, 'paired_choices_opposite': True, 'all_g0_choose_second': True, 'first_count': first_count, 'matrix_input_sha256': signature.hexdigest(), 'modulus': PRIME, 'rank_mod_prime': lower, 'rank_Q_direct_integer_elimination': qrank, 'explicit_integer_kernel_annihilates': annihilates, 'kernel_basis_count': len(basis), 'kernel_basis_rank_mod_prime': basis_lower, 'rank_Q_upper_bound_from_independent_integer_kernel': upper, 'rank_Q_certified_by_modular_kernel_sandwich': lower if lower == upper else None, 'scalar_constraint_rank': scalar_rank, 'scalar_rank_method': 'exact_Q_integer_elimination' if r <= 3 else 'exact_mod_prime_lower_bound', 'scalar_solution_dimension_upper_bound': n - scalar_rank, 'integer_order_scalar_solution_dimension_upper_bound_mod_prime': n - integer_scalar_rank, 'boundary_note': 'selected rank14 is not the known order4 variety dimension15' if r == 2 else None}

def build_certificate():
    require(PRIME > 2 and all((PRIME % divisor for divisor in range(2, isqrt(PRIME) + 1))), 'Independent invariant, original line 170')
    return {'schema': 'standalone-independent-sunlet-check-v1', 'independent_implementation': 'literal five-factor monomial construction; no owner imports', 'method': 'small exact Q rank; exact finite-field rank gives Q lower bound; explicit independent integer kernel gives Q upper bound', 'not_a_proof_of_uniform_quantifier': True, 'modulus_primality_checked_by_exhaustive_trial_division': True, 'cases': [one_case(r) for r in (1, 2, 3, 4, 7)]}
from archive_tools import require, run
if __name__ == '__main__':
    raise SystemExit(run(build_certificate, 'evidence/elementary-independent-certificate.json', 'independent literal monomials, exact modular/kernel rank sandwich'))
