"""Independent N=3 Fleming--Viot oracle using two-colour count states.

Only standard-library exact Fraction arithmetic. The 18-variable founder
operator and any symbolic current formula are deliberately not imported.
"""
from collections import defaultdict
from fractions import Fraction as F
from itertools import permutations
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def compositions(n):
    return [(a, b, n-a-b) for a in range(n+1) for b in range(n-a+1)]


def solve(a, b):
    """Exact Gauss--Jordan, returning the unique square-system solution."""
    n = len(a)
    require(len(b) == n and all(len(row) == n for row in a), 'square system')
    rows = [[F(x) for x in row] + [F(rhs)] for row, rhs in zip(a, b)]
    for col in range(n):
        pivot = next((i for i in range(col, n) if rows[i][col]), None)
        require(pivot is not None, 'singular exact system')
        rows[col], rows[pivot] = rows[pivot], rows[col]
        p = rows[col][col]
        rows[col] = [v / p for v in rows[col]]
        for i in range(n):
            if i != col and rows[i][col]:
                factor = rows[i][col]
                rows[i] = [x-factor*y for x, y in zip(rows[i], rows[col])]
    result = [row[-1] for row in rows]
    require(all(sum(F(v)*x for v, x in zip(row, result)) == rhs
                for row, rhs in zip(a, b)), 'solve residual')
    return result


POP = tuple(compositions(3))
MIXED = tuple(r+b for nred in (1, 2) for r in compositions(nred)
              for b in compositions(3-nred))
require(len(POP) == 10 and len(MIXED) == 36, 'state counts')


def coloured_jumps(state, killing):
    """Victim/donor multiplicities exclude the victim itself, even at same site."""
    jumps = defaultdict(F)
    for victim, count in enumerate(state):
        if not count:
            continue
        colour, site = divmod(victim, 3)
        # Natural motion is independent for every member of each colour.
        for destination in range(3):
            if destination != site:
                target = list(state)
                target[victim] -= 1
                target[3*colour+destination] += 1
                jumps[tuple(target)] += count
        # Each of the other two actual slots has donor rate k_site/2.
        for donor, donor_count in enumerate(state):
            available = donor_count - (donor == victim)
            if available <= 0:
                continue
            target = list(state)
            target[victim] -= 1
            target[donor] += 1
            target = tuple(target)
            if target != state:
                jumps[target] += F(count*available, 2) * killing[site]
    require(all(sum(t) == 3 and min(t) >= 0 and rate > 0
                for t, rate in jumps.items()), 'valid coloured jumps')
    return jumps


def colour_fixation(killing):
    index = {state:i for i, state in enumerate(MIXED)}
    a = [[F(0) for _ in MIXED] for _ in MIXED]
    b = [F(0) for _ in MIXED]
    for state, i in index.items():
        for target, rate in coloured_jumps(state, killing).items():
            a[i][i] += rate
            nred = sum(target[:3])
            if nred == 3:
                b[i] += rate
            elif nred == 0:
                pass
            else:
                a[i][index[target]] -= rate
    values = solve(a, b)
    require(all(0 < x < 1 for x in values), 'strict mixed fixation probabilities')
    result = dict(zip(MIXED, values))
    for state in MIXED:
        value = result[state]
        residual = F(0)
        for target, rate in coloured_jumps(state, killing).items():
            nred = sum(target[:3])
            after = F(1) if nred == 3 else F(0) if nred == 0 else result[target]
            residual += rate * (after-value)
        require(residual == 0, 'independent Dirichlet residual')
    return result


def population_stationary(killing):
    index = {state:i for i, state in enumerate(POP)}
    q = [[F(0) for _ in POP] for _ in POP]
    for state, row in index.items():
        for site, count in enumerate(state):
            if not count:
                continue
            for destination in range(3):
                if destination != site:
                    target = list(state)
                    target[site] -= 1
                    target[destination] += 1
                    # Natural moves, plus death and a donor at destination.
                    rate = count + F(count*state[destination], 2)*killing[site]
                    q[row][index[tuple(target)]] += rate
                    q[row][row] -= rate
    require(all(sum(row) == 0 for row in q), 'generator row sums')
    a = [[q[j][i] for j in range(len(POP))] for i in range(len(POP))]
    a[-1] = [F(1) for _ in POP]
    pi = solve(a, [F(0)]*(len(POP)-1)+[F(1)])
    require(sum(pi) == 1 and min(pi) > 0, 'stationary probability')
    require(all(sum(pi[j]*q[j][i] for j in range(len(POP))) == 0
                for i in range(len(POP))), 'full stationary residual')
    return dict(zip(POP, pi))


def compute(killing):
    killing = tuple(map(F, killing))
    require(len(killing) == 3 and min(killing) > 0, 'positive killing')
    fixation = colour_fixation(killing)
    pi = population_stationary(killing)
    def h(site, other):
        red = tuple(int(j == site) for j in range(3))
        return fixation[red+tuple(other)]
    for c in POP:
        normalization = F(0)
        for site, count in enumerate(c):
            if count:
                other = list(c); other[site] -= 1
                normalization += count*h(site, other)
        require(normalization == 1, 'founder normalization via two-colour oracle')
    transport = [[F(0) for _ in range(3)] for _ in range(3)]
    for source in range(3):
        for target in range(3):
            for c in POP:
                if c[source]:
                    other = list(c); other[source] -= 1
                    transport[source][target] += pi[c]*c[source]*h(target, other)
    current = transport[0][2]-transport[2][0]
    require(current == transport[2][1]-transport[1][2]
            == transport[1][0]-transport[0][1], 'cycle current conservation')
    return {'killing':list(map(str, killing)), 'delta':str(current),
            'population_states':len(POP), 'mixed_colour_states':len(MIXED),
            'population_pi':{str(c):str(pi[c]) for c in POP},
            'one_red_fixation':{str((i,d)):str(h(i,d)) for i in range(3)
                                for d in compositions(2)},
            'transport':[[str(x) for x in row] for row in transport],
            'all_exact_residuals_zero':True}


def main():
    points = list(dict.fromkeys(permutations((1, 2, 4))))
    points += [(1,1,4),(1,4,1),(4,1,1),(1,1,1),(2,2,2),
               (F(1,2),F(3,2),F(5,2)),(1,2,3),(2,4,6),
               (F(1,3),F(2,5),F(7,4))]
    results = [compute(point) for point in points]
    require(F(results[0]['delta']) == F(421413160,139117744582101),
            'source N3 baseline witness')
    baseline_pi = [103091,113786,157674,219448,138594,191246,264159,239661,329225,364330]
    require([F(results[0]['population_pi'][str(c)]) for c in POP]
            == [F(x,2121214) for x in baseline_pi], 'source stationary table')
    baseline_founders = [
        [138983140,121868982,108221552,108047078,96896756,87445362],
        [117160434,100452990,87445362,87483403,77057267,68542574],
        [87445362,72587826,61430106,61676473,52983701,46241930],
    ]
    require(all(F(results[0]['one_red_fixation'][str((i,d))])
                == F(baseline_founders[i][j],262336086)
                for i in range(3) for j,d in enumerate(compositions(2))),
            'source table of all eighteen founder probabilities')
    baseline_flux = {(0,1):113067935100535,(0,2):78610412314267,
                     (1,0):113068777926855,(1,2):64147123261485,
                     (2,0):78609569487947,(2,1):64147966087805}
    require(all(F(results[0]['transport'][i][j]) == F(value,278235489164202)
                for (i,j),value in baseline_flux.items()), 'source transport table')
    require(all(F(item['delta']) == 0 for item in results
                if len(set(item['killing'])) < 3), 'pair-equal finite controls')
    print(json.dumps({'method':'two-colour finite Dirichlet and separate uncoloured generator',
                      'exact_arithmetic':'standard-library Fraction',
                      'source_baseline_reproduced':True,'points':results},
                     indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
