"""Fresh referee verifier. Does not import or execute candidate code."""
import hashlib, itertools, json, math
from pathlib import Path
HERE = Path(__file__).resolve().parent

def exact_rank(integer_rows):
    rows = [r[:] for r in integer_rows]
    if not rows:
        return 0
    r = 0
    for col in range(len(rows[0])):
        p = next((i for i in range(r, len(rows)) if rows[i][col]), None)
        if p is None:
            continue
        rows[r], rows[p] = (rows[p], rows[r])
        pivot = rows[r][col]
        for i in range(r + 1, len(rows)):
            q = rows[i][col]
            if q:
                a = rows[i]
                for j in range(col, len(a)):
                    a[j] = pivot * a[j] - q * rows[r][j]
                divisor = math.gcd(*a)
                if divisor > 1:
                    rows[i] = [x // divisor for x in a]
        r += 1
        if r == len(rows):
            break
    return r

def check(moduli, orientation=0):
    elements = list(itertools.product(*(range(m) for m in moduli)))
    n = len(elements)
    idx = {g: i for i, g in enumerate(elements)}
    zero = elements[0]

    def add(g, h):
        return tuple(((a + b) % m for a, b, m in zip(g, h, moduli)))

    def neg(g):
        return tuple((-a % m for a, m in zip(g, moduli)))
    involutions = [g for g in elements[1:] if add(g, g) == zero]
    pairs = []
    seen = {zero, *involutions}
    for g in elements:
        if g not in seen:
            pairs.append((g, neg(g)))
            seen.update((g, neg(g)))
    if len(moduli) == 1:
        X = set(elements[1:n // 2 + 1])
    else:
        X = set(involutions)
        X.update((pair[orientation >> i & 1] for i, pair in enumerate(pairs)))
    Y = set(elements) - X - {zero}

    def generated(generators):
        reached = {zero}
        pending = [zero]
        while pending:
            h = pending.pop()
            for g in generators:
                k = add(h, g)
                if k not in reached:
                    reached.add(k)
                    pending.append(k)
        return reached
    weights = []
    for edge in range(6):
        for g in elements:
            if edge <= 2:
                v = 1
            elif edge == 3:
                v = 5 if g == zero else 1
            elif edge == 4:
                v = 3
            else:
                v = 2 if g == zero else 5 if g in X else 1
            weights.append(v)
    columns = []
    comparisons = []
    branch_ok = True
    for g, h in itertools.product(elements, repeat=2):
        gh = add(g, h)
        third = neg(gh)
        one = [0] * (6 * n)
        two = one[:]
        for edge, label in [(0, g), (1, h), (2, third), (3, gh), (4, g)]:
            one[edge * n + idx[label]] += 1
        for edge, label in [(0, g), (1, h), (2, third), (3, h), (5, g)]:
            two[edge * n + idx[label]] += 1
        weight_one = sum((a * b for a, b in zip(one, weights)))
        weight_two = sum((a * b for a, b in zip(two, weights)))
        if weight_one == weight_two:
            raise ValueError('Tie')
        first = weight_one < weight_two
        predicted = g in X and h != neg(g) or (g in Y and h == zero)
        branch_ok = branch_ok and first == predicted
        columns.append(one if first else two)
        comparisons.append({'g': g, 'h': h, 'first_weight': weight_one, 'second_weight': weight_two, 'choice': 1 if first else 2})
    matrix = [list(row) for row in zip(*columns)]
    rank = exact_rank(matrix)
    basis = []
    for free in range(n + 4):
        a = [int(free == i) for i in range(n)]
        B = int(free == n)
        C = int(free == n + 1)
        D = int(free == n + 2)
        e0 = int(free == n + 3)
        K = -B - C - D
        e = [K - x for x in a]
        e[0] = e0
        f = [K - x for x in a]
        basis.append(a + [B] * n + [C] * n + [D] * n + e + f)
    annihilation = all((sum((x * y for x, y in zip(v, col))) == 0 for v in basis for col in columns))
    record = {'moduli': moduli, 'orientation': orientation, 'n': n, 'shape': [6 * n, n * n], 'X': sorted(X), 'Y': sorted(Y), 'X_generates': len(generated(X)) == n, 'Y_generates': len(generated(Y)) == n, 'all_weights_positive': min(weights) > 0, 'minimum_absolute_margin': min((abs(x['first_weight'] - x['second_weight']) for x in comparisons)), 'no_ties': True, 'all_branch_predictions_match': branch_ok, 'first_count': sum((x['choice'] == 1 for x in comparisons)), 'second_count': sum((x['choice'] == 2 for x in comparisons)), 'rank_Q': rank, 'left_nullity': 6 * n - rank, 'expected_5n_minus_4': 5 * n - 4, 'proposed_kernel_basis_rank_Q': exact_rank(basis), 'proposed_basis_annihilates': annihilation, 'column_major_matrix_sha256': hashlib.sha256(json.dumps(columns, separators=(',', ':')).encode()).hexdigest()}
    return (record, {'elements': elements, 'weights': weights, 'columns': columns, 'comparisons': comparisons})

def build_certificate():
    cases = [([6], 0), ([2, 4], 0), ([2, 4], 1), ([2, 4], 2), ([2, 4], 3), ([2, 2, 3], 0), ([4, 4], 0), ([4], 0), ([2, 2], 0)]
    records = []
    exact_inputs = []
    for mods, orientation in cases:
        record, data = check(mods, orientation)
        records.append(record)
        exact_inputs.append({'moduli': mods, 'orientation': orientation, **data})
        expected = 15 if mods == [4] else 13 if mods == [2, 2] else 5 * record['n'] - 4
        if record['rank_Q'] != expected or record['proposed_kernel_basis_rank_Q'] != record['n'] + 4 or (not record['proposed_basis_annihilates']) or (not record['no_ties']) or (not record['all_weights_positive']) or (not record['all_branch_predictions_match']):
            raise RuntimeError('Independent exact diagnostic failed')
    receipt = {'method': 'integer row elimination with gcd normalization over Q; separate monomial construction from literal parameterization; no candidate code imported or executed', 'finite_diagnostics_only': True, 'cases': records}
    return receipt

from archive_tools import run
if __name__ == "__main__":
    raise SystemExit(run(build_certificate,"evidence/independent-certificate.json","independent literal monomials and integer elimination"))
