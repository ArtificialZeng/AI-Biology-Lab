"""Finite graph-contract checker; no stochastic-image claim is encoded."""
from itertools import product
from collections import deque
import json

LEAVES = {"1", "2", "3", "4"}
BLOB = {"u", "v", "a", "b"}
RETICS = {"a", "b"}
DIRECTED = [("u", "a"), ("v", "a"), ("u", "b"), ("v", "b")]
UNDIRECTED = [("u", "v"), ("a", "s"), ("b", "t"),
              ("s", "1"), ("s", "2"), ("t", "3"), ("t", "4")]
EDGES = DIRECTED + UNDIRECTED
VERTICES = set(sum(([x, y] for x, y in EDGES), []))

def require(predicate, message):
    if not predicate:
        raise ValueError(message)

def neighbors(edges):
    adj = {v: set() for v in set(sum(([x, y] for x, y in edges), []))}
    for x, y in edges:
        adj[x].add(y)
        adj[y].add(x)
    return adj

def reachable(adj, start, avoid=None):
    if start == avoid:
        return set()
    found = {start}
    todo = [start]
    while todo:
        x = todo.pop()
        for y in adj[x]:
            if y != avoid and y not in found:
                found.add(y)
                todo.append(y)
    return found

def is_binary_rooted(edges):
    vertices = VERTICES | {"rho"}
    out = {v: [] for v in vertices}
    indegree = {v: 0 for v in vertices}
    for x, y in edges:
        out[x].append(y)
        indegree[y] += 1
    for v in vertices:
        wanted = ((0, 2) if v == "rho" else (1, 0) if v in LEAVES
                  else (2, 1) if v in RETICS else (1, 2))
        if (indegree[v], len(out[v])) != wanted:
            return False
    work = dict(indegree)
    todo = deque([v for v in vertices if work[v] == 0])
    seen = []
    while todo:
        x = todo.popleft()
        seen.append(x)
        for y in out[x]:
            work[y] -= 1
            if work[y] == 0:
                todo.append(y)
    if len(seen) != len(vertices) or reachable(out, "rho") != vertices:
        return False
    # A nonroot lies on all root-to-leaf paths iff deleting it leaves no leaf reachable.
    if any(not (reachable(out, "rho", v) & LEAVES)
           for v in vertices if v != "rho"):
        return False
    return True

def unique_path(edges, start, end):
    adj = neighbors(edges)
    todo = deque([(start, [start])])
    seen = {start}
    while todo:
        x, path = todo.popleft()
        if x == end:
            return path
        for y in sorted(adj[x]):
            if y not in seen:
                seen.add(y)
                todo.append((y, path + [y]))
    raise ValueError("Disconnected displayed tree")

def run():
    require(len({frozenset(e) for e in EDGES}) == len(EDGES), "Parallel edge")
    adj = neighbors(EDGES)
    require(len(reachable(adj, "u")) == len(VERTICES), "Disconnected graph")
    boundary = {v for v in BLOB if adj[v] - BLOB}
    require(boundary == {"a", "b"}, "Wrong blob boundary")
    internal = [e for e in EDGES if set(e) <= BLOB]
    require(len(internal) == 5, "Wrong diamond")
    require(reachable(neighbors(internal), "u") == BLOB, "Blob disconnected")
    bridges = []
    for i, edge in enumerate(EDGES):
        reduced = EDGES[:i] + EDGES[i+1:]
        if edge[1] not in reachable(neighbors(reduced), edge[0]):
            bridges.append(list(edge))
    require(not any(set(e) <= BLOB for e in bridges), "Blob contains cut edge")

    rootings = []
    failed_root_edges = []
    # Exhaust all edge-insertion rootings without hidden parallel/degree-two expansions.
    # An analytic cut/direction argument handles the source's more general rootings.
    for root_edge in EDGES:
        fixed = [e for e in DIRECTED if e != root_edge]
        remaining = [e for e in UNDIRECTED if e != root_edge]
        fixed += [("rho", root_edge[0]), ("rho", root_edge[1])]
        accepted = []
        for bits in product((0, 1), repeat=len(remaining)):
            edges = fixed + [(y, x) if bit else (x, y)
                             for (x, y), bit in zip(remaining, bits)]
            if not is_binary_rooted(edges):
                continue
            parent_sets = {v: sorted(x for x, y in edges if y == v) for v in RETICS}
            displayed = []
            for pa, pb in product(parent_sets["a"], parent_sets["b"]):
                chosen = {"a": pa, "b": pb}
                tree = [(x, y) for x, y in edges
                        if y not in chosen or chosen[y] == x]
                require(len(tree) == len(VERTICES), "Displayed edge count")
                require(len(reachable(neighbors(tree), "rho")) == len(VERTICES) + 1,
                        "Displayed graph not tree")
                path = unique_path(tree, "a", "b")
                require(len(path) - 1 >= 2, "Unexpected direct attachment edge")
                displayed.append({"a_parent": pa, "b_parent": pb,
                                  "attachment_path": path,
                                  "original_edge_count": len(path)-1})
            accepted.append({"root_edge": list(root_edge),
                             "directed_edges": [list(e) for e in sorted(edges)],
                             "displayed_trees": displayed})
        if accepted:
            rootings += accepted
        else:
            failed_root_edges.append(list(root_edge))
    require(len(rootings) == 5, "Expected exactly five simple edge rootings")
    require(all(set(r["root_edge"]) <= BLOB for r in rootings), "Outside rooting exists")
    require(len(failed_root_edges) == 6, "Expected six outside-root exclusions")
    collapsed = [("s", "t"), ("s", "1"), ("s", "2"), ("t", "3"), ("t", "4")]
    lengths = [t["original_edge_count"] for r in rootings for t in r["displayed_trees"]]
    return {"graph_kind": "own explicit level-2 binary semidirected construction",
            "source": "Brits et al. arXiv:2607.12919v3 Definitions2.1 and5.1",
            "vertices": sorted(VERTICES), "reticulations": sorted(RETICS),
            "directed_reticulation_edges": DIRECTED, "undirected_edges": UNDIRECTED,
            "blob": sorted(BLOB), "attachment_vertices": sorted(boundary),
            "cut_edges": bridges, "collapsed_quartet_edges": collapsed,
            "edge_insertion_rootings": rootings,
            "excluded_outside_root_edges": failed_root_edges,
            "rooting_count": len(rootings), "displayed_path_count": len(lengths),
            "minimum_displayed_attachment_path_edges": min(lengths),
            "maximum_displayed_attachment_path_edges": max(lengths),
            "claim_scope": "Finite graph/domain certificate only; no kernel or image nonmembership claim"}

if __name__ == "__main__":
    print(json.dumps(run(), sort_keys=True, indent=2))
