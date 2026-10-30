#!/usr/bin/env python3
"""Standalone exact check from state-space matrices. Requires SymPy; emits JSON to stdout."""
import itertools
import json
import sympy as S


def require(test, label):
    if not test:
        raise ValueError(label)


def cycle_matrix(k, leaks):
    n = len(k)
    A = S.zeros(n)
    for j, rate in enumerate(k):
        A[j, j] -= rate + leaks.get(j + 1, 0)
        A[(j + 1) % n, j] += rate
    return A


def det_derivative(M, V):
    return S.expand(sum(
        (-1) ** (i + j) * M.minor_submatrix(i, j).det() * V[i, j]
        for i in range(M.rows) for j in range(M.cols) if V[i, j] != 0
    ))


def main():
    s = S.Symbol('s')
    k1, k2, k3, leak = S.symbols('k1 k2 k3 leak')
    M = s * S.eye(3) - cycle_matrix([k1, k2, k3], {2: leak})
    D = S.expand(M.det())
    # Cramer's rule: delete input row and output column, retaining the cofactor sign.
    N = S.expand(-M.minor_submatrix(0, 1).det())
    coeff = S.Matrix([D.coeff(s, q) for q in (2, 1, 0)] + [N.coeff(s, q) for q in (1, 0)])
    jac = coeff.jacobian([k1, k2, k3, leak])
    point = {k1: 2, k2: 1, k3: 3, leak: 1}
    require(k2.subs(point) + leak.subs(point) == k1.subs(point), 'hyperplane')
    minor = S.factor(jac.extract([0, 2, 3, 4], range(4)).det())
    require(minor == -k1**2 * k3, 'symbolic minor')
    require(minor.subs(point) == -12, 'minor value')
    require(jac.subs(point).rank() == 4, 'rank')
    require(S.gcd(D.subs(point), N.subs(point)) == 1, 'coprimality')
    require(D.subs(s, -k3) == -k1*k2*k3, 'independent cancellation test')
    # Independent held-out state-space tests: no formula for N or D is supplied.
    # These include p=n, leaks at p, ell=1, and >2-leak algebraic examples.
    cases = []
    for n in range(3, 7):
        for p in range(3, n + 1):
            for a, ell in ((1, p - 1), (p - 1, 1)):
                for leaks in (set([ell]), set([ell, p]), set(range(1, n + 1)) - {a}):
                    rates = list(map(S.Integer, range(2, n + 2)))
                    leak_values = {i: S.Integer(i + 1) for i in leaks}
                    rates[a - 1] = rates[ell - 1] + leak_values[ell]
                    dk = [S.Integer(0)] * n
                    dk[a - 1], dk[ell - 1] = rates[a - 1], -rates[ell - 1]
                    de = {ell: rates[ell - 1] - rates[a - 1]}
                    A = cycle_matrix(rates, leak_values)
                    dA = cycle_matrix(dk, de)
                    M0 = s*S.eye(n) - A
                    dM = -dA
                    dD = det_derivative(M0, dM)
                    minorM = M0.minor_submatrix(0, p - 1)
                    minorV = dM.minor_submatrix(0, p - 1)
                    dN = (-1)**(p + 1) * det_derivative(minorM, minorV)
                    require(dD == 0 and S.expand(dN) == 0, 'state-space derivative')
                    cases.append({'n': n, 'p': p, 'a': a, 'ell': ell, 'leaks': sorted(leaks)})
    print(json.dumps({
        'status': 'passed', 'endpoint_D': str(D), 'endpoint_N': str(N),
        'coefficient_values': list(map(int, coeff.subs(point))),
        'selected_minor': str(minor), 'selected_minor_value': int(minor.subs(point)),
        'rank': int(jac.subs(point).rank()),
        'all_maximal_minors_at_point': [int(jac.extract(rows, range(4)).det().subs(point)) for rows in itertools.combinations(range(5), 4)],
        'gcd_at_point': str(S.gcd(D.subs(point), N.subs(point))),
        'held_out_state_space_derivative_cases': len(cases),
        'unique_held_out_state_space_cases': len({(c['n'], c['p'], c['a'], c['ell'], tuple(c['leaks'])) for c in cases}),
        'held_out_scope': 'n=3..6; all p>=3; first/last interior positions interchanged; 1,2,and n-1 leak patterns; one positive integer point per case',
        'held_out_cases': cases,
        'limits': ['Finite state-space checks supplement the written all-n proof.', 'No candidate file was imported or overwritten.']
    }, indent=2))


if __name__ == '__main__':
    main()
