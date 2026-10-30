"""Independent exact Pearson oracle over explicit graph and allele partitions.
No tskit statistics, tree traversal, genotype decoder, variance or correlation routine
is used to form expected values. Run with a chosen native tskit on sys.path.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json
import math
import sys
import numpy as np
import tskit
import _tskit

def pearson_squared(g, w):
    n = len(w)
    a = n * sum((Q(x) * Q(y) for (x, y) in zip(g, w))) - sum(g) * sum(map(Q, w))
    b = n * sum((Q(x) ** 2 for x in g)) - sum(g) ** 2
    c = n * sum((Q(y) ** 2 for y in w)) - sum(map(Q, w)) ** 2
    if c <= 0:
        raise ValueError('trait variance is not positive')
    return Q(0) if b == 0 else a * a / (b * c)

def fixture(n):
    m = n // 2
    times = [Q(0)] * (n + 1) + [Q(2), Q(3), Q(5)]
    (left, right, root) = (n + 1, n + 2, n + 3)
    descendants = [set()] + [{j} for j in range(n)] + [set(range(m)), set(range(m, n)), set(range(n))]
    edges = [(left, j + 1) for j in range(m)] + [(right, j + 1) for j in range(m, n)] + [(root, left), (root, right)]
    alleles = [(Q(2), descendants[right]), (Q(8), descendants[n])]
    tab = tskit.TableCollection(sequence_length=10)
    for (j, t) in enumerate(times):
        tab.nodes.add_row(time=float(t), flags=tskit.NODE_IS_SAMPLE if 1 <= j <= n else 0)
    for (parent, child) in edges:
        tab.edges.add_row(0, 10, parent, child)
    for (pos, clade_node) in [(2, right), (8, n)]:
        site = tab.sites.add_row(pos, ancestral_state='A')
        tab.mutations.add_row(site=site, node=clade_node, derived_state='T')
    tab.sort()
    return (tab.tree_sequence(), times, descendants, edges, alleles)

def expected_values(n, traits, mode, windows, normalise, data):
    (times, desc, edges, alleles) = data
    inheritance = lambda s: [int(i in s) for i in range(n)]
    out = []
    for (lo, hi) in zip(windows[:-1], windows[1:]):
        span = Q(hi) - Q(lo)
        if mode == 'site':
            row = []
            for w in traits:
                val = Q(0)
                for (pos, s) in alleles:
                    if lo <= pos < hi:
                        val += (pearson_squared(inheritance(s), w) + pearson_squared(inheritance(set(range(n)) - s), w)) / 2
                row.append(val / span if normalise else val)
        elif mode == 'branch':
            row = []
            for w in traits:
                val = sum(((times[p] - times[c]) * pearson_squared(inheritance(desc[c]), w) for (p, c) in edges))
                row.append(val if normalise else val * span)
        elif mode == 'node':
            row = [[pearson_squared(inheritance(s), w) * (1 if normalise else span) for w in traits] for s in desc]
        else:
            raise ValueError(mode)
        out.append(row)
    return np.array(out, dtype=float)

def main():
    results = []
    for n in (2, 4, 6):
        (ts, *data) = fixture(n)
        for exponent in (None, 20, 26, 27, 30, 40):
            for sign in (1,) if exponent is None else (-1, 1):
                offset = 0 if exponent is None else sign * 2 ** exponent
                traits = [[offset + i for i in range(n)], [-offset + (i * i + 2 * i) for i in range(n)]]
                W = np.array(traits, dtype=float).T
                for mode in ('site', 'branch', 'node'):
                    for windows in ([0, 10], [0, 4, 10]):
                        for norm in (False, True):
                            expected = expected_values(n, traits, mode, windows, norm, data)
                            observed = ts.trait_correlation(W, windows=windows, mode=mode, span_normalise=norm)
                            ok = observed.shape == expected.shape and bool(np.allclose(observed, expected, rtol=2e-12, atol=2e-12))
                            results.append({'n': n, 'offset': offset, 'mode': mode, 'windows': windows, 'span_normalise': norm, 'matches': ok, 'observed': observed.tolist() if np.all(np.isfinite(observed)) else 'nonfinite', 'expected': expected.tolist()})
    receipt = {'tskit_version': tskit.__version__, 'tskit_module': tskit.__file__, 'extension': _tskit.__file__, 'extension_sha256': hashlib.sha256(Path(_tskit.__file__).read_bytes()).hexdigest(), 'cases': len(results), 'passes': sum((x['matches'] for x in results)), 'failures': sum((not x['matches'] for x in results)), 'four_sample_site_oracle': 'For g=(0,0,1,1),w=(0,1,2,3)+c, rational r_squared=4/5', 'results': results}
    print(json.dumps(receipt, indent=2, allow_nan=False))
if __name__ == '__main__':
    main()
