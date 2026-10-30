"""Independent exact checks: standard library only; no original checker import.

Uses resolvents at positive rational frequencies and dR = R(dA)R, rather
than Markov recurrences. All numbers are fractions.Fraction. Finite checks
support, and do not prove, the arbitrary-tree theorem.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json


def require(value, message):
    if not value:
        raise ArithmeticError(message)


def inverse(a):
    n = len(a)
    b = [[F(x) for x in row] + [F(i == j) for j in range(n)]
         for i, row in enumerate(a)]
    for j in range(n):
        p = next(i for i in range(j, n) if b[i][j])
        b[j], b[p] = b[p], b[j]
        d = b[j][j]
        b[j] = [v / d for v in b[j]]
        for i in range(n):
            if i != j and b[i][j]:
                d = b[i][j]
                b[i] = [x - d*y for x, y in zip(b[i], b[j])]
    return [row[n:] for row in b]


def rank_pivots(a):
    b = [[F(x) for x in row] for row in a]
    m, n = len(b), len(b[0])
    q, cols = 0, []
    for j in range(n):
        p = next((i for i in range(q, m) if b[i][j]), None)
        if p is None:
            continue
        b[q], b[p] = b[p], b[q]
        d = b[q][j]
        b[q] = [x / d for x in b[q]]
        for i in range(q + 1, m):
            if b[i][j]:
                d = b[i][j]
                b[i] = [x - d*y for x, y in zip(b[i], b[q])]
        cols.append(j)
        q += 1
        if q == min(m, n):
            break
    return q, cols


def determinant(a):
    b = [[F(x) for x in row] for row in a]
    n, d = len(b), F(1)
    for j in range(n):
        p = next((i for i in range(j, n) if b[i][j]), None)
        if p is None:
            return F(0)
        if p != j:
            b[j], b[p] = b[p], b[j]
            d = -d
        d *= b[j][j]
        for i in range(j + 1, n):
            if b[i][j]:
                t = b[i][j] / b[j][j]
                for k in range(j + 1, n):
                    b[i][k] -= t*b[j][k]
                b[i][j] = F(0)
    return d


def model(n, edges, leaks):
    directed = [(u, v) for u, v in edges for u, v in ((u, v), (v, u))]
    params = directed + [(v, None) for v in leaks]
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]
    values = primes[:len(params)]
    a = [[F(0) for _ in range(n)] for _ in range(n)]
    for (source, target), value in zip(params, values):
        a[source][source] -= value
        if target is not None:
            a[target][source] += value
    return a, params, values


def path(n, edges, start, end):
    todo, parent = [start], {start: None}
    while todo:
        v = todo.pop()
        if v == end:
            break
        for x, y in edges:
            u = y if x == v else x if y == v else None
            if u is not None and u not in parent:
                parent[u] = v
                todo.append(u)
    out = [end]
    while out[-1] != start:
        out.append(parent[out[-1]])
    return out[::-1]


def path_edges(nodes):
    return {tuple(sorted((u, v))) for u, v in zip(nodes, nodes[1:])}


def derivatives(a, params, root, sensor):
    n = len(a)
    rows, rootrows, labels = [], [], []
    for z in range(1, 2*n + 1):
        r = inverse([[F(z*(i == j)) - a[i][j] for j in range(n)]
                     for i in range(n)])
        for v in (root, sensor):
            row = [(r[v][target] - r[v][source])*r[source][root]
                   if target is not None else -r[v][source]*r[source][root]
                   for source, target in params]
            rows.append(row)
            labels.append([z, v])
            if v == root:
                rootrows.append(row)
    return rows, rootrows, labels


def harmonic(a, leaks):
    n, (left, right) = len(a), leaks
    interior = [v for v in range(n) if v not in leaks]
    mat, rhs = [], []
    for source in interior:
        mat.append([F(sum(a[v][source] for v in range(n) if v != source))
                    if target == source else -a[target][source]
                    for target in interior])
        rhs.append(a[right][source])
    inv = inverse(mat) if interior else []
    h = [F(0) for _ in range(n)]
    h[right] = F(1)
    for i, v in enumerate(interior):
        h[v] = sum(inv[i][j]*rhs[j] for j in range(len(interior)))
    return h


def run_case(n, edges, leaks, root, sensor, include_minor=False):
    a, params, values = model(n, edges, leaks)
    j, jr, labels = derivatives(a, params, root, sensor)
    rank, cols = rank_pivots(j)
    rr, _ = rank_pivots(jr)
    shared = sorted(path_edges(path(n, edges, root, sensor)) &
                    path_edges(path(n, edges, *leaks)))
    expected = 2*n if shared else 2*n - 1
    require(rank == expected, 'joint resolvent derivative rank mismatch')
    require(rr == 2*n - 1, 'root resolvent derivative rank mismatch')
    h = harmonic(a, leaks)
    require(bool(h[root] != h[sensor]) == bool(shared), 'voltage criterion mismatch')
    result = dict(n=n, edges=edges, leaks=leaks, root=root, sensor=sensor,
                  shared_edges=shared, root_rank=rr, joint_rank=rank,
                  voltage=[str(x) for x in h], params=params, values=values)
    if include_minor:
        _, selected = rank_pivots(list(map(list, zip(*j))))
        minor = [[j[i][c] for c in cols] for i in selected]
        det = determinant(minor)
        require(det != 0, 'selected minor is zero')
        result['minor_row_labels'] = [labels[i] for i in selected]
        result['minor_parameter_columns'] = cols
        result['minor_determinant'] = str(det)
    if not shared:
        t = F(1, 97)
        scales = [1+t*(v-h[root]) for v in h]
        b = [[scales[i]*a[i][k]/scales[k] for k in range(n)] for i in range(n)]
        loss = [-sum(b[i][k] for i in range(n)) for k in range(n)]
        require(all(v > 0 for v in scales), 'nonpositive similarity scale')
        require(all(loss[v] > 0 if v in leaks else loss[v] == 0 for v in range(n)),
                'leak domain violated')
        require(a != b and scales[root] == scales[sensor] == 1, 'curve not a witness')
        # The exact identity BQ=QA, Qe_r=e_r, and c_v Q=c_v proves equality
        # at every frequency, independent of any sampled transfer evaluations.
        require(all(b[i][k]*scales[k] == scales[i]*a[i][k]
                    for i in range(n) for k in range(n)), 'similarity identity failed')
        result['positive_exact_ambiguity'] = dict(t=str(t), scales=[str(x) for x in scales],
                                                 leaks=[str(x) for x in loss])
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    # A seven-vertex branched tree differs from the original six-vertex control.
    edges = [(0,1), (1,2), (2,3), (1,4), (4,5), (4,6)]
    selected = [run_case(7, edges, (3,6), 0,5, True),
                run_case(7, edges, (3,6), 0,1, True),
                run_case(7, edges, (3,6), 3,6),
                run_case(7, edges, (0,6), 0,3),
                run_case(7, edges, (0,6), 3,0)]
    # Explicit smallest-tree and all placement checks on n=3,4 representatives.
    small = []
    for n, ee in [(2, [(0,1)]), (3,[(0,1),(1,2)]),
                  (4,[(0,1),(1,2),(2,3)]), (4,[(0,1),(0,2),(0,3)])]:
        counts = {'n': n, 'edges': ee, 'cases': 0, 'full': 0, 'corank_one': 0}
        for a in range(n):
            for b in range(a+1,n):
                for r in range(n):
                    for s in range(n):
                        if s == r:
                            continue
                        c = run_case(n,ee,(a,b),r,s)
                        counts['cases'] += 1
                        counts['full' if c['joint_rank'] == 2*n else 'corank_one'] += 1
        small.append(counts)
    result = dict(method='exact rational resolvent derivatives; standard library only',
                  generic_theorem_proved_by_computation=False,
                  seven_vertex_controls=selected, exhaustive_small_trees=small,
                  all_checks_passed=True,
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    output = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(output)
    else:
        print(output, end='')


if __name__ == '__main__':
    main()
