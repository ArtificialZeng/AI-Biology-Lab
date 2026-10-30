"""Exact Pybus--Harvey interval oracle; imports no TreeSwift code.

Uses lineage-count jumps at exact rational node depths rather than resolving
polytomies into binary nodes. The numerator below is an algebraic reordering of
the Pybus--Harvey cumulative-interval formula implemented by ape::gammaStat.
"""
from dataclasses import dataclass, field
from fractions import Fraction
from math import sqrt
import re

@dataclass
class Node:
    children: list = field(default_factory=list)
    label: str = ''
    edge: Fraction = Fraction(0)


def parse_newick(text):
    tokens = re.findall(r'[(),:;]|[^\s(),:;]+', text)
    cursor = 0
    def parse():
        nonlocal cursor
        n = Node()
        if tokens[cursor] == '(':
            cursor += 1
            n.children.append(parse())
            while tokens[cursor] == ',':
                cursor += 1
                n.children.append(parse())
            assert tokens[cursor] == ')'
            cursor += 1
        if cursor < len(tokens) and tokens[cursor] not in ',():;':
            n.label = tokens[cursor]
            cursor += 1
        if tokens[cursor] == ':':
            cursor += 1
            n.edge = Fraction(tokens[cursor])
            cursor += 1
        return n
    root = parse()
    assert tokens[cursor:] == [';']
    return root


def gamma_from_intervals(intervals):
    """Return exact coefficient c, radicand r, gamma=c*sqrt(r)."""
    g = [Fraction(x) for x in intervals]
    n = len(g) + 1
    if n < 3:
        raise ValueError('At least three extant tips required')
    total = sum(k*x for k,x in enumerate(g, 2))
    if total == 0:
        raise ZeroDivisionError('Gamma is undefined at zero total lineage time')
    # Each k*g_k occurs in n-k of the cumulative sums in the published formula.
    centered = sum(k*(n-k)*x for k,x in enumerate(g, 2))/Fraction(n-2) - total/2
    coefficient = centered/total
    radicand = 12*(n-2)
    return {'intervals': [str(x) for x in g], 'tip_count': n,
            'total_lineage_time': str(total),
            'coefficient': str(coefficient), 'radicand': radicand,
            'squared_value': str(coefficient*coefficient*radicand),
            'value': float(coefficient)*sqrt(radicand)}


def oracle_newick(text):
    root = parse_newick(text)
    while len(root.children) == 1:
        root = root.children[0]
    events = {}
    tips = []
    stack = [(root, Fraction(0))]
    while stack:
        node,depth = stack.pop()
        if node.children:
            if len(node.children) > 1:
                events[depth] = events.get(depth, 0) + len(node.children)-1
            stack.extend((child, depth+child.edge) for child in node.children)
        else:
            tips.append(depth)
    assert len(tips) >= 3
    assert len(set(tips)) == 1, 'Scientific oracle fixtures must be ultrametric'
    assert min(events) == 0
    horizon = tips[0]
    times = sorted(events)
    assert horizon >= times[-1]
    # Count extant lineages on each open interval. A degree-d branching creates
    # d-1 additional lineages at once; omitted intermediate intervals have zero length.
    intervals = {k:Fraction(0) for k in range(2,len(tips)+1)}
    count = 1
    for i,time in enumerate(times):
        count += events[time]
        stop = times[i+1] if i+1 < len(times) else horizon
        intervals[count] = stop-time
    assert count == len(tips)
    result = gamma_from_intervals([intervals[k] for k in range(2,len(tips)+1)])
    result['event_jumps'] = {str(k):v for k,v in sorted(events.items())}
    result['tip_depths'] = [str(t) for t in tips]
    return result
