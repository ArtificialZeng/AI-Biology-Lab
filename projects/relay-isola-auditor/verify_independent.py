#!/usr/bin/env python3
"""Independent exact audit. Read-only: every run prints JSON and writes no files.

Uses integer monomial dictionaries rather than the owner's SymPy identities.
The topological argument is in proof.md, not claimed machine-formalized.
"""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OWNER = ROOT


def require(condition, message):
    if not condition:
        raise ValueError(message)


def clean(p):
    return {m: c for m, c in p.items() if c}


def add(p, q):
    out = dict(p)
    for m, c in q.items():
        out[m] = out.get(m, 0) + c
    return clean(out)


def scale(p, n):
    return clean({m: n*c for m, c in p.items()})


def derivative(p, i):
    out = {}
    for m, c in p.items():
        if m[i]:
            a = list(m)
            a[i] -= 1
            out[tuple(a)] = c*m[i]
    return clean(out)


def value(p, z):
    return sum(c * z[0]**m[0] * z[1]**m[1] * z[2]**m[2]
               for m, c in p.items())


def y_is_lambda(p):
    out = {}
    for (a, b, c), v in p.items():
        m = (a, 0, b+c)
        out[m] = out.get(m, 0) + v
    return clean(out)


def x_is_zero(p):
    return {m: c for m, c in p.items() if m[0] == 0}


def x_is_four(p):
    out = {}
    for (a, b, c), v in p.items():
        m = (0, b, c)
        out[m] = out.get(m, 0) + v*4**a
    return clean(out)


def main():
    data = json.loads((OWNER / 'evidence/exact-certificate.json').read_text())
    expected = [
        ((1,0),(0,0),'7'), ((2,0),(3,0),'4'),
        ((3,0),(2,0),'1'), ((1,1),(2,1),'4'),
        ((1,2),(0,0),'1'), ((1,2),(1,4),'1'),
        ((0,0),(0,1),'lambda'), ((0,1),(0,0),'1'),
    ]
    actual = [(tuple(r['source']), tuple(r['product']), r['constant'])
              for r in data['reactions']]
    require(actual == expected, 'The serialized reaction input changed')
    f = [{}, {}]
    consumed_columns = []
    for index, (u, v, k) in enumerate(actual, 1):
        k_degree = int(k == 'lambda')
        coefficient = 1 if k_degree else int(k)
        require(coefficient > 0, 'Nonpositive fixed coefficient')
        rate = {(u[0], u[1], k_degree): coefficient}
        negative = {i for i in (0,1) if v[i] < u[i]}
        dependency = {i for i in (0,1) if u[i] > 0}
        if negative:
            require(dependency <= negative, f'Net-consumption support failure in {index}')
            require(all(u[i] >= 1 for i in negative), f'Consumed factor absent in {index}')
            consumed_columns.append(index)
        for i in (0,1):
            f[i] = add(f[i], scale(rate, v[i]-u[i]))
    require(f[0] == {(1,0,0):-7, (2,0,0):4, (3,0,0):-1,
                     (1,1,0):4, (1,2,0):-1}, 'Wrong first field component')
    require(f[1] == {(0,0,1):1, (0,1,0):-1}, 'Wrong second field component')
    require(all(m[0] >= 1 for m in f[0]), 'Cannot divide first component by x')
    H = {(a-1,b,c): v for (a,b,c),v in f[0].items()}
    HL = y_is_lambda(H)
    require(HL == {(0,0,0):-7, (1,0,0):4, (2,0,0):-1,
                   (0,0,1):4, (0,0,2):-1}, 'Wrong reduced conic')
    # 1 - (x-2)^2 - (lambda-2)^2 has exactly these coefficients.
    J = [[derivative(f[i],j) for j in (0,1)] for i in (0,1)]
    boundary = [[y_is_lambda(x_is_zero(J[i][j])) for j in (0,1)] for i in (0,1)]
    require(boundary == [[{(0,0,0):-7,(0,0,1):4,(0,0,2):-1},{}],
                         [{},{(0,0,0):-1}]], 'Wrong boundary Jacobian')
    require(add(y_is_lambda(J[0][0]), scale(HL,-1)) ==
            {(1,0,0):4,(2,0,0):-2}, 'Wrong interior eigenvalue')
    siphons = []
    for mask in (1,2,3):
        subset = {i for i in (0,1) if mask & (1 << i)}
        if all(not any(v[i] > 0 for i in subset) or any(u[i] > 0 for i in subset)
               for u,v,k in actual):
            siphons.append(sorted(subset))
    require(siphons == [[0]], 'Unexpected siphon classification')
    fold_output = []
    for L in (1,3):
        z = (2,L,L)
        A = [[value(J[i][j],z) for j in (0,1)] for i in (0,1)]
        b = 4 if L == 1 else -4
        require(A == [[0,b],[0,-1]], 'Fold Jacobian mismatch')
        p = (1,0)
        q = (1,b)
        require([sum(A[i][j]*p[j] for j in (0,1)) for i in (0,1)] == [0,0],
                'Right nullvector mismatch')
        require([sum(q[i]*A[i][j] for i in (0,1)) for j in (0,1)] == [0,0],
                'Left nullvector mismatch')
        require(sum(p[i]*q[i] for i in (0,1)) == 1, 'Normalization mismatch')
        f_lambda = [value(derivative(f[i],2),z) for i in (0,1)]
        quadratic = [value(derivative(derivative(f[i],0),0),z) for i in (0,1)]
        a = sum(q[i]*f_lambda[i] for i in (0,1))
        c = sum(q[i]*quadratic[i] for i in (0,1))
        require(a == b and a != 0, 'Parameter transversality fails')
        require(c == -4 and c != 0, 'Quadratic nondegeneracy fails')
        fold_output.append({'lambda':L, 'J':A, 'q':q, 'q_f_lambda':a, 'q_f_xx':c})
    require(x_is_four(f[0]) == {(0,0,0):-28,(0,1,0):16,(0,2,0):-4},
            'Uniform x-bound polynomial fails')
    # Same polynomial is -12-4(y-2)^2; at y=4, f2=lambda-4.
    require(max(sum(u) for u,v,k in actual) == 3, 'Molecularity changed')
    # Independently check the reviewer's full-rate-space limitation path.
    # Substitute y=lambda=2 in the integer field, then replace X-loss 7 by a(t).
    univariate = {}
    for (power_x, power_y, power_lam), coefficient in f[0].items():
        univariate[power_x] = univariate.get(power_x, 0) + coefficient * 2**(power_y+power_lam)
    univariate = clean(univariate)
    require(univariate == {1:-3, 2:4, 3:-1}, 'Full-rate slice polynomial differs')
    rate_path = {0:4, 1:4, 2:-1}
    path_field = dict(univariate)
    path_field[1] = path_field.get(1, 0) + 7
    for power, coefficient in rate_path.items():
        path_field[power+1] = path_field.get(power+1, 0) - coefficient
    require(clean(path_field) == {}, 'Full-rate path is not equilibria')
    require(sum(c*3**i for i,c in rate_path.items()) == 7, 'Path does not start on original circle')
    require(rate_path[0] == 4, 'Path boundary limit rate changed')
    digest = {str(p.relative_to(OWNER)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in [OWNER/'source.md', OWNER/'proof.md', OWNER/'verify_independent.py',
                        OWNER/'evidence/exact-certificate.json']}
    print(json.dumps({'verified':True, 'method':'integer monomial dictionaries; no asserts',
                      'ordinary_operation':'read-only; stdout JSON only',
                      'consumed_columns':consumed_columns, 'siphons_zero_based':siphons,
                      'folds':fold_output, 'full_rate_path_identity_checked':True,
                      'reviewed_sha256':digest}, indent=2))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, KeyError, TypeError) as exc:
        print(json.dumps({'verified':False,'error':str(exc)}, sort_keys=True))
        raise SystemExit(1)
