"""Exact-arithmetic falsification checks for the frozen two-leak tree claim.

The finite-field ranks are exact lower-bound certificates over Q, not a proof
of the infinite theorem. The explicit positive rational similarity witness and
selected rational Jacobian minor are separate exact checks.
"""
import argparse
import itertools
import hashlib
import json
from pathlib import Path

import networkx as nx
import numpy as np
import sympy as sp

P = 1000003


def check(condition, message="exact verification failed"):
    if not condition:
        raise ValueError(message)


def rank_mod(mat, p=P):
    a = np.array(mat, dtype=np.int64).copy() % p
    rank = 0
    for col in range(a.shape[1]):
        pivots = np.flatnonzero(a[rank:, col])
        if not len(pivots):
            continue
        pivot = rank + int(pivots[0])
        a[[rank, pivot]] = a[[pivot, rank]]
        a[rank] = (a[rank] * pow(int(a[rank, col]), -1, p)) % p
        for row in range(rank + 1, a.shape[0]):
            if a[row, col]:
                a[row] = (a[row] - a[row, col] * a[rank]) % p
        rank += 1
        if rank == min(a.shape):
            break
    return rank


def model(g, leaks, seed=0, symbolic=False):
    n = len(g)
    directed = [(i, j) for i, j in sorted(g.edges()) for i, j in [(i, j), (j, i)]]
    params = directed + [(i, None) for i in leaks]
    values = [2 + int(hashlib.sha256(f"tree-sensor:{seed}:{j}".encode()).hexdigest()[:8], 16) % 997
              for j in range(len(params))]
    if symbolic:
        a = sp.zeros(n)
    else:
        a = np.zeros((n, n), dtype=np.int64)
    for (source, target), value in zip(params, values):
        a[source, source] -= value
        if target is not None:
            a[target, source] += value
    return a, params, values


def markov_jac(g, leaks, root, sensor, seed=0, rational=False):
    n = len(g)
    a, params, values = model(g, leaks, seed, symbolic=rational)
    if rational:
        v = sp.eye(n)[:, root]
        deriv = sp.zeros(n, len(params))
    else:
        v = np.eye(n, dtype=np.int64)[:, root]
        deriv = np.zeros((n, len(params)), dtype=np.int64)
    rows = []
    for k in range(2 * n):
        rows += [list(deriv[root, :]) if rational else deriv[root].copy(),
                 list(deriv[sensor, :]) if rational else deriv[sensor].copy()]
        next_deriv = a * deriv if rational else (a @ deriv) % P
        for j, (source, target) in enumerate(params):
            next_deriv[source, j] -= v[source]
            if target is not None:
                next_deriv[target, j] += v[source]
        deriv = next_deriv if rational else next_deriv % P
        v = a * v if rational else (a @ v) % P
    return sp.Matrix(rows) if rational else np.array(rows), a, params, values


def path_edges(g, u, v):
    nodes = nx.shortest_path(g, u, v)
    return {tuple(sorted(e)) for e in zip(nodes, nodes[1:])}


def finite_atlas():
    summary = []
    failures = []
    for n in range(2, 7):
        count = good = bad = 0
        for ti, g in enumerate(nx.nonisomorphic_trees(n)):
            g = nx.convert_node_labels_to_integers(g, ordering="sorted")
            for leaks in itertools.combinations(range(n), 2):
                lp = path_edges(g, *leaks)
                for root in range(n):
                    for sensor in range(n):
                        if root == sensor:
                            continue
                        predicted = 2 * n if lp & path_edges(g, root, sensor) else 2 * n - 1
                        jac, *_ = markov_jac(g, leaks, root, sensor)
                        observed = rank_mod(jac)
                        for seed in (1, 2):
                            if observed == predicted:
                                break
                            # An unlucky point does not refute a generic rank.
                            jac2, *_ = markov_jac(g, leaks, root, sensor, seed=seed)
                            observed = max(observed, rank_mod(jac2))
                        count += 1
                        good += int(predicted == 2 * n)
                        bad += int(predicted == 2 * n - 1)
                        if observed != predicted:
                            failures.append({"n": n, "tree_index": ti, "edges": list(g.edges()),
                                             "leaks": leaks, "root": root, "sensor": sensor,
                                             "predicted": predicted, "observed": observed})
        summary.append({"n": n, "cases": count, "full_rank_cases": good, "corank_one_cases": bad})
    return {"prime": P, "domain": "all nonisomorphic undirected trees, 2<=n<=6; every root, unordered two-leak set, and distinct additional sensor", "cases": summary, "mismatches": failures}


def rational_controls():
    # Branching at both 0 and 1, beyond the path and star special cases.
    g = nx.Graph([(0, 1), (1, 2), (1, 3), (0, 4), (0, 5)])
    jac, a, params, values = markov_jac(g, (2, 3), 0, 2, rational=True)
    _, independent_rows = jac.T.rref()
    rows = list(independent_rows)
    determinant = jac[rows, :].det()
    check(len(rows) == 12 and determinant != 0)

    bad_jac, *_ = markov_jac(g, (2, 3), 0, 4, rational=True)
    check(bad_jac.rank() == 11)
    # Original outgoing-rate harmonic voltage h_2=0, h_3=1.
    h = sp.symbols("h0:6")
    equations = [h[2], h[3] - 1]
    for source in (0, 1, 4, 5):
        equations.append(sum(a[target, source] * (h[target] - h[source]) for target in g.neighbors(source)))
    sol = sp.solve(equations, h)
    voltage = sp.Matrix([sol[x] for x in h])
    t = sp.Rational(1, 101)
    scales = sp.ones(6, 1) + t * (voltage - sp.ones(6, 1) * voltage[0])
    d = sp.diag(*scales)
    new_a = d * a * d.inv()
    new_leaks = [-sum(new_a[:, j]) for j in range(6)]
    check(all(scales[i] > 0 for i in range(6)))
    check(all(new_a[v, u] > 0 for u, v in params if v is not None))
    check(new_leaks[2] > 0 and new_leaks[3] > 0)
    check(all(new_leaks[j] == 0 for j in (0, 1, 4, 5)))
    check(new_a != a and scales[0] == scales[4] == 1)
    # All 2n Markov terms verify the exact transfer equality; similarity proves
    # equality at all orders independently of this finite regression check.
    v = sp.eye(6)[:, 0]
    new_v = v.copy()
    for _ in range(12):
        check(v[0] == new_v[0] and v[4] == new_v[4])
        v = a * v
        new_v = new_a * new_v
    # A sensor at leak 2 must detect this same nontrivial curve.
    check(scales[2] != 1)
    # Deliberately wrong criterion: intersection at a vertex alone is insufficient.
    vertex_jac, *_ = markov_jac(g, (2, 3), 0, 1, rational=True)
    check(vertex_jac.rank() == 11)
    check(set(nx.shortest_path(g, 2, 3)) & set(nx.shortest_path(g, 0, 1)) == {1})
    return {"edges": list(g.edges()), "root": 0, "leaks": [2, 3], "good_sensor": 2,
            "bad_sensor": 4, "parameters_source_target_then_leaks": params,
            "parameter_values": values, "good_markov_jacobian_minor_rows": rows,
            "good_markov_jacobian_minor_determinant": str(determinant),
            "bad_exact_rank": 11, "harmonic_voltage": [str(x) for x in voltage],
            "vertex_only_intersection_sensor": 1, "vertex_only_exact_rank": 11,
            "t": str(t), "diagonal_scales": [str(x) for x in scales],
            "transformed_leaks": [str(x) for x in new_leaks],
            "exact_positive_indistinguishable_witness": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="write a separate receipt; default stdout")
    args = parser.parse_args()
    result = {"finite_atlas": finite_atlas(), "rational_controls": rational_controls()}
    check(not result["finite_atlas"]["mismatches"], "finite atlas has mismatches")
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
