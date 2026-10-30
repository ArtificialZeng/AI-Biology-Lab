#!/usr/bin/env python3
"""Fresh referee: raw parent/target transitions, no owner scientific imports.

Run with /opt/anaconda3/bin/python raw_referee.py > evidence/raw-referee.json.
Fixtures are constructed here and serialized in the result. Exact generator
tests use integer scaling; exact fixation uses an independent Dirichlet solve.
"""
import json
import math
from fractions import Fraction as F

from sympy import ZZ
from sympy.polys.matrices import DomainMatrix


class Graph:
    """Minimal undirected graph from serialized edges; no fixture generator."""
    def __init__(self, n, edges):
        self.n = n
        self.edge_list = tuple(tuple(edge) for edge in edges)
        self.adj = [set() for _ in range(n)]
        for a, b in self.edge_list:
            assert 0 <= a < b < n and b not in self.adj[a]
            self.adj[a].add(b)
            self.adj[b].add(a)
    def __len__(self):
        return self.n
    def __iter__(self):
        return iter(range(self.n))
    def edges(self):
        return self.edge_list
    def neighbors(self, vertex):
        return sorted(self.adj[vertex])
    def degree(self, vertex):
        return len(self.adj[vertex])
    def connected(self):
        seen = {0}
        pending = [0]
        while pending:
            for vertex in self.adj[pending.pop()] - seen:
                seen.add(vertex)
                pending.append(vertex)
        return len(seen) == self.n


def fixtures():
    # Five connected noncomplete fixtures were checked against the frozen
    # owner's fixture list up to graph isomorphism in the separate input receipt.
    return [
        ('cycle8-heldout', Graph(8, [[0, 1], [0, 7], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]]), True, True),
        ('complement-cycle8-heldout', Graph(8, [[0, 2], [0, 3], [0, 4], [0, 5], [0, 6], [1, 3], [1, 4], [1, 5], [1, 6], [1, 7], [2, 4], [2, 5], [2, 6], [2, 7], [3, 5], [3, 6], [3, 7], [4, 6], [4, 7], [5, 7]]), True, True),
        ('complement-C3-union-C5-heldout', Graph(8, [[0, 3], [0, 4], [0, 5], [0, 6], [0, 7], [1, 3], [1, 4], [1, 5], [1, 6], [1, 7], [2, 3], [2, 4], [2, 5], [2, 6], [2, 7], [3, 5], [3, 6], [4, 6], [4, 7], [5, 7]]), True, True),
        ('seeded-cubic10-heldout', Graph(10, [[0, 1], [0, 2], [0, 4], [1, 4], [1, 7], [2, 5], [2, 9], [3, 4], [3, 7], [3, 8], [5, 6], [5, 9], [6, 8], [6, 9], [7, 8]]), True, False),
        ('seeded-quartic10-heldout', Graph(10, [[0, 1], [0, 4], [0, 7], [0, 8], [1, 2], [1, 5], [1, 7], [2, 3], [2, 5], [2, 6], [3, 5], [3, 7], [3, 9], [4, 5], [4, 8], [4, 9], [6, 7], [6, 8], [6, 9], [8, 9]]), True, False),
        ('complete5-control', Graph(5, [[0, 1], [0, 2], [0, 3], [0, 4], [1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]]), False, True),
        ('complete2-control', Graph(2, [[0, 1]]), False, True),
        ('two-triangles-outside-scope', Graph(6, [[0, 1], [0, 2], [1, 2], [3, 4], [3, 5], [4, 5]]), False, True),
    ]


def raw_rows(graph, reverse=False):
    n = len(graph)
    d = graph.degree(0)
    scale = math.lcm(d, n - 1)
    full = (1 << n) - 1
    neighbors = [tuple(graph.neighbors(i)) for i in range(n)]
    complete = [tuple(j for j in range(n) if j != i) for i in range(n)]
    rows = []
    for state in range(full + 1):
        row = {}
        for parent in range(n):
            parent_mutant = bool(state & (1 << parent))
            use_complete = parent_mutant if reverse else not parent_mutant
            targets = complete[parent] if use_complete else neighbors[parent]
            weight = scale // len(targets)
            for target in targets:
                if parent_mutant:
                    outcome = state | (1 << target)
                else:
                    outcome = state & ~(1 << target)
                row[outcome] = row.get(outcome, 0) + weight
        assert sum(row.values()) == n * scale
        rows.append(row)
    assert rows[0] == {0: n * scale}
    assert rows[full] == {full: n * scale}
    return scale, rows


def stats(graph, state):
    n = len(graph)
    d = graph.degree(0)
    k = state.bit_count()
    e = sum(bool(state & (1 << a)) and bool(state & (1 << b))
            for a, b in graph.edges())
    gs = [sum(bool(state & (1 << u)) for u in graph.neighbors(v))
          for v in graph if not (state & (1 << v))]
    # W = 2(N-1)V. E = (N-1)^2 D.
    w = d * k * (2*n - 1 + k) - 2*(n - 1)*e
    error = sum(((n - 1)*g - d*k)**2 for g in gs)
    return k, e, gs, w, error


def exact_solve(rows, scale, n):
    """Dirichlet equations from raw events, all labelled transient states."""
    size = (1 << n) - 2
    full = size + 1
    a = {}
    b = {}
    for state in range(1, full):
        i = state - 1
        a[i] = {i: n * scale}
        b[i] = {0: rows[state].get(full, 0)}
        for destination, weight in rows[state].items():
            if 0 < destination < full:
                j = destination - 1
                a[i][j] = a[i].get(j, 0) - weight
    matrix = DomainMatrix(a, (size, size), ZZ)
    rhs = DomainMatrix(b, (size, 1), ZZ)
    num, den = matrix.solve_den(rhs, method="rref")
    vector = num.to_Matrix()
    answer = [F(0)] + [F(int(vector[i, 0]), int(den)) for i in range(size)] + [F(1)]
    # Verify the returned rational solution against every raw labelled row.
    for state, row in enumerate(rows):
        assert n * scale * answer[state] == sum(w * answer[t] for t, w in row.items())
    return answer


def one_case(name, graph, heldout, solve):
    n = len(graph)
    degrees = {graph.degree(v) for v in graph}
    assert len(degrees) == 1
    d = degrees.pop()
    q = n - 1 - d
    full = (1 << n) - 1
    scale, forward = raw_rows(graph)
    _, reverse = raw_rows(graph, reverse=True)
    data = [stats(graph, state) for state in range(full + 1)]
    for state in range(full + 1):
        k, e, gs, w, error = data[state]
        cut = sum(gs)
        assert cut == d*k - 2*e
        row = forward[state]
        drift_k = sum(weight * (data[t][0] - k) for t, weight in row.items())
        drift_k2 = sum(weight * (data[t][0]**2 - k*k) for t, weight in row.items())
        drift_e = sum(weight * (data[t][1] - e) for t, weight in row.items())
        u = F(cut, d)
        loss = F(k*(n-k), n-1)
        assert F(drift_k, scale) == u-loss
        assert F(drift_k2, scale) == 2*k*(u-loss)+u+loss
        assert F(drift_e, scale) == F(sum(g*g for g in gs),d)-F(2*e*(n-k),n-1)
        drift_w = sum(weight*(data[t][3]-w) for t, weight in row.items())
        assert d*(n-1)*drift_w == -2*scale*error
        # Exact pathwise color reversal of the complete transition row.
        assert reverse[state] == {full ^ t: weight for t, weight in forward[full ^ state].items()}
    c = F(d*n*n, n-1)
    assert data[0][3] == 0
    assert data[full][3] == 2*(n-1)*c
    for vertex in graph:
        singleton = 1 << vertex
        near_full = full ^ singleton
        assert F(data[singleton][3],2*(n-1)) == c/n
        assert F(data[singleton][4],(n-1)**2) == F(d*q,n-1)
        assert F(data[near_full][3],2*(n-1)) == c*F(n-1,n)
        assert data[near_full][4] == 0
        expected_d1 = F(sum(weight*data[t][4] for t,weight in forward[near_full].items()),
                        n*scale*(n-1)**2)
        assert expected_d1 == F(2*d*q,n*(n-1)**2)
    result = {"name": name, "N": n, "d": d, "connected": graph.connected(),
              "held_out_from_owner_up_to_isomorphism": heldout,
              "edges": sorted([sorted(edge) for edge in graph.edges()]),
              "raw_states": full+1, "exact_generator_and_reversal": True}
    if solve:
        f = exact_solve(forward, scale, n)
        r = exact_solve(reverse, scale, n)
        assert all(f[s] == 1-r[full^s] for s in range(full+1))
        for v in graph:
            assert f[1<<v] <= F(1,n)-F(q,d*n**3)
            assert r[1<<v] >= F(1,n)+F(2*q,d*n**4*(n-1))
            if q == 0:
                assert f[1<<v] == r[1<<v] == F(1,n)
        result.update({"exact_fixation_both_orientations": True,
                       "forward_singletons": [str(f[1<<v]) for v in graph],
                       "reverse_singletons": [str(r[1<<v]) for v in graph],
                       "all_singleton_bounds_pass": True})
    else:
        result["exact_fixation_both_orientations"] = False
    return result


def main():
    tests = [one_case(name, graph, heldout, solve)
             for name, graph, heldout, solve in fixtures()]
    print(json.dumps({"schema":"fresh-raw-labelled-referee-v1",
                      "private_or_owner_input_required":False,
                      "candidate_code_imported":False, "test_count":len(tests),
                      "total_raw_states_per_orientation":sum(t["raw_states"] for t in tests),
                      "tests":tests}, indent=2))


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run without Python optimization so all exact assertions execute")
    main()
