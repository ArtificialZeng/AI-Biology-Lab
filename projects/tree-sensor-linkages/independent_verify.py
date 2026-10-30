#!/usr/bin/env python3
"""Independent finite exact audit; standard library only; default writes nothing.

Uses rational elimination, exhaustive vertex separators (not path packing),
forest enumeration, and resolvent derivatives (not Markov recurrences).
This is finite evidence, not a proof of an all-positive infinite theorem.
"""
from fractions import Fraction as F
from itertools import combinations
import json
import sys


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def rref(rows, ncols=None):
    a = [[F(v) for v in row] for row in rows]
    width = len(a[0]) if a else (ncols or 0)
    pivots = []
    for j in range(width if ncols is None else ncols):
        p = next((i for i in range(len(pivots), len(a)) if a[i][j]), None)
        if p is None:
            continue
        i = len(pivots)
        a[i], a[p] = a[p], a[i]
        scale = a[i][j]
        a[i] = [v / scale for v in a[i]]
        for k in range(len(a)):
            if k != i and a[k][j]:
                scale = a[k][j]
                a[k] = [x - scale * y for x, y in zip(a[k], a[i])]
        pivots.append(j)
        if len(pivots) == len(a):
            break
    return a, pivots


def rank(rows):
    return len(rref(rows)[1])


def inverse(a):
    n = len(a)
    if not n:
        return []
    aug = [list(a[i]) + [int(i == j) for j in range(n)] for i in range(n)]
    red, pivots = rref(aug, n)
    require(len(pivots) == n, "singular exact matrix")
    return [row[n:] for row in red]


def det(a):
    a = [[F(v) for v in row] for row in a]
    out = F(1)
    for j in range(len(a)):
        p = next((i for i in range(j, len(a)) if a[i][j]), None)
        if p is None:
            return F(0)
        if p != j:
            a[j], a[p] = a[p], a[j]
            out = -out
        v = a[j][j]
        out *= v
        for i in range(j + 1, len(a)):
            c = a[i][j] / v
            for k in range(j + 1, len(a)):
                a[i][k] -= c * a[j][k]
    return out


def subsets(n, nonempty=False):
    return [tuple(v for v in range(n) if (mask >> v) & 1)
            for mask in range(int(nonempty), 1 << n)]


def adjacency(n, edges):
    adj = [set() for _ in range(n)]
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    return adj


def cut_number(n, edges, observations, boundary):
    """Menger min vertex separator, including source/target vertices."""
    adj = adjacency(n, edges)
    for size in range(n + 1):
        for cut in combinations(range(n), size):
            forbidden = set(cut)
            seen = set(observations) - forbidden
            todo = list(seen)
            for v in todo:
                for u in adj[v] - forbidden - seen:
                    seen.add(u)
                    todo.append(u)
            if not seen.intersection(boundary):
                return size
    raise RuntimeError("vertex separator enumeration failed")


def rates_for(edges, variant):
    arcs = sorted([(u, v) for u, v in edges] + [(v, u) for u, v in edges])
    return {arc: F((j + 2) ** 2 + 2 * variant + (j % 3),
                   j % 4 + 1) for j, arc in enumerate(arcs)}


def harmonic(n, edges, boundary, rates):
    b = list(boundary)
    inside = [j for j in range(n) if j not in b]
    adj = adjacency(n, edges)
    k = [[F(0) for _ in range(n)] for _ in range(n)]
    for j in range(n):
        k[j][j] = sum((rates[j, i] for i in adj[j]), F(0))
        for i in adj[j]:
            k[j][i] = -rates[j, i]
    m = [([F(int(j == i)) for i in range(n)] if j in b else k[j][:])
         for j in range(n)]
    mi = inverse(m)
    return [[mi[j][i] for i in b] for j in range(n)], k, m


def components(n, chosen_edges):
    adj = adjacency(n, chosen_edges)
    left = set(range(n))
    result = []
    while left:
        root = min(left)
        seen = {root}
        todo = [root]
        for v in todo:
            for u in adj[v] - seen:
                seen.add(u)
                todo.append(u)
        result.append(seen)
        left -= seen
    return result


def forest_check(n, edges, boundary, w, h, k, m, rates):
    """Check full-minor existence, unique pairing, and positive forest magnitude."""
    b = tuple(boundary)
    pairings = set()
    forest_sum = F(0)
    count = 0
    for chosen in combinations(edges, n - len(b)):
        comps = components(n, chosen)
        if len(comps) != len(b):
            continue
        if any(len(c.intersection(b)) != 1 or len(c.intersection(w)) != 1
               for c in comps):
            continue
        pairing = tuple(sorted((next(iter(c.intersection(w))),
                                next(iter(c.intersection(b)))) for c in comps))
        pairings.add(pairing)
        adj = adjacency(n, chosen)
        weight = F(1)
        for c in comps:
            root = next(iter(c.intersection(b)))
            seen = {root}
            todo = [root]
            for parent in todo:
                for child in adj[parent] - seen:
                    weight *= rates[child, parent]
                    seen.add(child)
                    todo.append(child)
        forest_sum += weight
        count += 1
    minor = det([[h[i][j] for j in range(len(b))] for i in w])
    numerator = det([[k[i][j] for j in range(n) if j not in w]
                     for i in range(n) if i not in b])
    require(len(pairings) <= 1, "more than one full linkage pairing on tree")
    require(abs(minor * det(m)) == abs(numerator), "Jacobi complement mismatch")
    require(abs(numerator) == forest_sum, "toward-boundary forest sum mismatch")
    require(bool(minor) == bool(count), "full minor/forest existence mismatch")
    return count


def physical_matrix(n, edges, boundary, rates, leaks):
    a = [[F(0) for _ in range(n)] for _ in range(n)]
    for (j, i), v in rates.items():
        a[i][j] += v
        a[j][j] -= v
    for b in boundary:
        a[b][b] -= leaks[b]
    return a


def resolvent(a, z):
    return inverse([[F(z) * (i == j) - a[i][j] for j in range(len(a))]
                    for i in range(len(a))])


def physical_rank(n, edges, boundary, root, observations, rates, leaks):
    a = physical_matrix(n, edges, boundary, rates, leaks)
    arcs = sorted(rates)
    rows = []
    for z in range(1, 2 * n):
        g = resolvent(a, z)
        for o in observations:
            rows.append([(g[o][i] - g[o][j]) * g[j][root] for j, i in arcs]
                        + [-g[o][b] * g[b][root] for b in boundary])
    return rank(rows)


def null_basis(rows, width):
    red, pivots = rref(rows, width)
    free = [j for j in range(width) if j not in pivots]
    result = []
    for f in free:
        v = [F(0)] * width
        v[f] = F(1)
        for i, p in enumerate(pivots):
            v[p] = -red[i][f]
        result.append(v)
    return result


def gauge_check(n, edges, boundary, root, observations, rates, leaks):
    h, _, _ = harmonic(n, edges, boundary, rates)
    bases = null_basis([h[o] for o in observations], len(boundary))
    a = physical_matrix(n, edges, boundary, rates, leaks)
    for bv in bases:
        hv = [sum((x * y for x, y in zip(row, bv)), F(0)) for row in h]
        require(hv[root] == 0, "root gauge not fixed")
        t = F(1)
        for _ in range(100):
            d = [1 + t * v for v in hv]
            if min(d) > 0:
                new_rates = {(j, i): v * d[i] / d[j] for (j, i), v in rates.items()}
                new_leaks = {b: -a[b][b] - sum((v for (j, i), v in new_rates.items()
                                               if j == b), F(0)) for b in boundary}
                if min(new_leaks.values()) > 0:
                    break
            t /= 2
        else:
            raise RuntimeError("failed to find positive rational gauge step")
        anew = physical_matrix(n, edges, boundary, new_rates, new_leaks)
        require(anew != a, "nonzero gauge did not change parameters")
        require(all(anew[i][j] == d[i] * a[i][j] / d[j]
                    for i in range(n) for j in range(n)), "gauge left physical model")
        for z in range(1, n + 2):
            g, gnew = resolvent(a, z), resolvent(anew, z)
            require(all(gnew[o][root] == g[o][root] for o in observations),
                    "observed transfer changed under exact gauge")
            require(all(gnew[o][root] == d[o] * g[o][root] for o in range(n)),
                    "resolvent similarity identity failed")
    return len(bases)


def run():
    branch = (6, ((0, 1), (1, 2), (1, 3), (3, 4), (3, 5)))
    harmonic_tests = 0
    forest_tests = 0
    forest_terms = 0
    sensor_constructions = 0
    for n, edges in [(1, ()), (4, ((0, 1), (1, 2), (2, 3))), branch]:
        for variant in [0, 1]:
            rates = rates_for(edges, variant)
            for b in subsets(n, True):
                h, k, m = harmonic(n, edges, b, rates)
                require(all(sum(row) == 1 and min(row) >= 0 for row in h),
                        "harmonic maximum principle failed")
                for o in subsets(n, True):
                    q = cut_number(n, edges, o, b)
                    require(rank([h[i] for i in o]) == q,
                            f"harmonic rank mismatch: {n}, {b}, {o}, {variant}")
                    harmonic_tests += 1
                for w in combinations(range(n), len(b)):
                    forest_terms += forest_check(n, edges, b, w, h, k, m, rates)
                    forest_tests += 1
                for root in range(n):
                    removed = next(j for j, value in zip(b, h[root]) if value > 0)
                    obs = tuple(sorted({root}.union(set(b) - {removed})))
                    require(len(obs) == len(b), "sensor count construction failed")
                    require(cut_number(n, edges, obs, b) == len(b),
                            "minimal sensor construction did not yield full linkage")
                    sensor_constructions += 1
    n, edges = branch
    physical_cases = [
        (1, (), (0,), 0, (0,), "one-vertex"),
        (n, edges, (0, 2, 4, 5), 1, (1,), "branch-root-only"),
        (n, edges, (0, 2, 4, 5), 1, (1, 3), "branch-bottleneck"),
        (n, edges, (0, 2, 4, 5), 1, (0, 1, 3), "boundary-intersection"),
        (n, edges, (0, 2, 4, 5), 1, (0, 1, 2, 4), "full-linkage"),
        (n, edges, (0, 1, 2, 3, 4, 5), 1, (1, 3), "all-leaks"),
        (n, edges, (1, 4, 5), 1, (1, 2, 4), "root-leak"),
        (4, ((0, 1), (1, 2), (2, 3)), (1, 3), 0, (0, 2), "first-leak-truncation"),
    ]
    records = []
    gauge_directions = 0
    for n, edges, b, root, obs, name in physical_cases:
        for variant in [0, 1]:
            rates = rates_for(edges, variant)
            leaks = {v: F(7 + 3 * j + variant, 2 + j % 2) for j, v in enumerate(b)}
            q = cut_number(n, edges, obs, b)
            got = physical_rank(n, edges, b, root, obs, rates, leaks)
            want = 2 * n - 2 + q
            require(got == want, f"physical rank mismatch: {name}, {variant}, {got}, {want}")
            ngauge = gauge_check(n, edges, b, root, obs, rates, leaks)
            require(ngauge == len(b) - q, f"gauge dimension mismatch: {name}")
            gauge_directions += ngauge
            records.append({"name": name, "variant": variant, "n": n,
                            "boundary": b, "root": root, "observations": obs,
                            "q": q, "physical_rank": got, "gauge_dimension": ngauge})
    edges = ((0, 1), (1, 2), (2, 3))
    h, _, _ = harmonic(4, edges, (1, 3), rates_for(edges, 0))
    require(h[0][1] == 0 and h[0][0] == 1, "first boundary absorption check failed")
    cycle = ((0, 1), (1, 2), (2, 3), (3, 0))
    unit_rates = {(j, i): F(1) for u, v in cycle for j, i in [(u, v), (v, u)]}
    hc, _, _ = harmonic(4, cycle, (0, 2), unit_rates)
    cycle_rank = rank([hc[o] for o in (1, 3)])
    cycle_q = cut_number(4, cycle, (1, 3), (0, 2))
    require(cycle_rank == 1 and cycle_q == 2, "non-tree cancellation control failed")
    symmetric_star = ((0, 1), (0, 2))
    symmetric_rates = {(j, i): F(1) for u, v in symmetric_star
                       for j, i in [(u, v), (v, u)]}
    nongeneric_rank = physical_rank(3, symmetric_star, (0,), 0, (0,),
                                   symmetric_rates, {0: F(1)})
    require(nongeneric_rank == 3, "positive nongeneric rank control failed")
    return {"status": "passed", "arithmetic": "exact fractions; Python standard library",
            "method": {"linkage": "exhaustive vertex separators",
                       "harmonic": "rational inverse of boundary-replaced Laplacian",
                       "forest": "enumerate tree edge subsets; orient toward boundary",
                       "physical": "direct resolvent derivative at 2n-1 positive integer z points",
                       "gauge": "exact rational positive similarity; direct physical matrix and resolvent comparison"},
            "harmonic_placements": harmonic_tests, "full_forest_minor_checks": forest_tests,
            "forest_terms": forest_terms, "minimal_sensor_constructions": sensor_constructions,
            "physical_cases": records,
            "nonzero_gauge_directions_checked": gauge_directions,
            "partial_minor_control": {"tree": "0-1-2-3", "boundary": [1, 3],
                                      "H_0_to_3": str(h[0][1]), "H_0_to_1": str(h[0][0])},
            "non_tree_control": {"cycle": "0-1-2-3-0", "boundary": [0, 2],
                                 "observations": [1, 3], "harmonic_rank": cycle_rank,
                                 "linkage": cycle_q},
            "nongeneric_positive_control": {"tree": "star 1-0-2", "root": 0,
                                            "boundary": [0], "all_rates_and_leak": 1,
                                            "physical_rank": nongeneric_rank,
                                            "generic_rank": 5},
            "limits": "Finite rational checks only. Infinite/all-positive and generic statements require the written proof."}


if __name__ == "__main__":
    try:
        print(json.dumps(run(), indent=2))
    except Exception as exc:
        print(json.dumps({"status": "failed", "error": str(exc)}), file=sys.stderr)
        raise SystemExit(1)
