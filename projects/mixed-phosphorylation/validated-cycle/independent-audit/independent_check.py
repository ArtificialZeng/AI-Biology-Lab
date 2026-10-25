"""Read-only reconstruction and exact postprocessing for the limiting oscillator.

Uses only Python's standard library. This does not replace CAPD integration:
it independently checks the equation reconstruction and arithmetic consequences
of deliberately widened serialized CAPD derivative/residual enclosures.
Writes its report only into this dedicated fresh-audit directory.
"""
from fractions import Fraction as Q
from pathlib import Path
import ast
import hashlib
import itertools
import json
import re
import argparse

CANDIDATE = Path(__file__).resolve().parents[1]
ROOT = CANDIDATE.parents[2]
LOG = CANDIDATE / 'evidence/linux-interval-output.txt'
N = 7
ZERO = (0,) * N


class Polynomial:
    def __init__(self, value=0):
        self.terms = {ZERO: Q(value)} if value else {}

    @staticmethod
    def from_terms(terms):
        result = Polynomial()
        result.terms = {key: value for key, value in terms.items() if value}
        return result

    def __add__(self, other):
        other = other if isinstance(other, Polynomial) else Polynomial(other)
        terms = dict(self.terms)
        for key, value in other.terms.items():
            terms[key] = terms.get(key, 0) + value
        return self.from_terms(terms)

    __radd__ = __add__

    def __neg__(self):
        return self.from_terms({key: -value for key, value in self.terms.items()})

    def __sub__(self, other):
        return self + -(other if isinstance(other, Polynomial) else Polynomial(other))

    def __rsub__(self, other):
        return Polynomial(other) + -self

    def __mul__(self, other):
        other = other if isinstance(other, Polynomial) else Polynomial(other)
        terms = {}
        for key, value in self.terms.items():
            for key2, value2 in other.terms.items():
                product = tuple(key[index] + key2[index] for index in range(N))
                terms[product] = terms.get(product, 0) + value * value2
        return self.from_terms(terms)

    __rmul__ = __mul__

    def __eq__(self, other):
        other = other if isinstance(other, Polynomial) else Polynomial(other)
        return self.terms == other.terms

    def substitute(self, index, value):
        result = Polynomial()
        for key, coefficient in self.terms.items():
            reduced = list(key)
            exponent = reduced[index]
            reduced[index] = 0
            term = self.from_terms({tuple(reduced): coefficient})
            for unused in range(exponent):
                term = term * value
            result = result + term
        return result


def polynomial_checks():
    variables = []
    for index in range(N):
        exponent = list(ZERO)
        exponent[index] = 1
        variables.append(Polynomial.from_terms({tuple(exponent): Q(1)}))
    u, v, w, a, b, y, epsilon = variables
    p = 8 - a - b
    s = 40 - u - v - w - a - b
    # Independently transcribed from the published arrows, ordered species:
    # S0, K, S0K, S1K, S2, P, S2P, S1, S1P.
    before = [
        {0: 1, 1: 1}, {2: 1}, {2: 1}, {3: 1}, {4: 1, 5: 1},
        {6: 1}, {6: 1}, {7: 1, 5: 1}, {8: 1}, {8: 1},
    ]
    after = [
        {2: 1}, {0: 1, 1: 1}, {3: 1}, {4: 1, 1: 1}, {6: 1},
        {4: 1, 5: 1}, {7: 1, 5: 1}, {8: 1}, {7: 1, 5: 1}, {0: 1, 5: 1},
    ]
    # First binding flux is exactly (epsilon*y)*(epsilon^-1-u-v).
    fluxes = [
        y - epsilon*y*(u+v), u, u, v, 100*w*p, a, Q(9, 10)*a,
        3*(s-epsilon*y)*p, b, 100*b,
    ]
    full = [Polynomial() for unused in range(9)]
    for reactants, products, flux in zip(before, after, fluxes):
        for index, amount in reactants.items():
            full[index] = full[index] - amount*flux
        for index, amount in products.items():
            full[index] = full[index] + amount*flux
    slow = [full[index].substitute(6, Polynomial()).substitute(5, u+100*b)
            for index in [2, 3, 4, 6, 8]]
    source = (CANDIDATE / 'verify_cycle.cpp').read_text()
    expression = re.search(
        r'IMap vf\("par:r;var:u,v,w,a,b;fun:([^;]+);"\)', source).group(1)
    environment = dict(zip(['u', 'v', 'w', 'a', 'b'], variables))
    environment['r'] = Polynomial(Q(19, 10))

    def parse(node):
        if isinstance(node, ast.Name):
            return environment[node.id]
        if isinstance(node, ast.Constant):
            return Polynomial(node.value)
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            return -parse(node.operand)
        if isinstance(node, ast.BinOp):
            left, right = parse(node.left), parse(node.right)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
        raise ValueError(f'Unexpected expression node: {type(node).__name__}')

    encoded = [parse(ast.parse(part, mode='eval').body) for part in expression.split(',')]
    matches = [left == right for left, right in zip(slow, encoded)]
    fast = full[0] == -y + epsilon*(u+v)*y + u + 100*b
    p_identity = (-slow[3]-slow[4]) == -(100*w+3*s)*p + Q(19, 10)*a + 101*b
    s_identity = (-sum(slow, Polynomial())) == Q(9, 10)*a + b - 3*p*s
    assert all(matches) and fast and p_identity and s_identity
    return {'slow_rhs_equal_checker': matches, 'fast_identity': fast,
            'dependent_p_identity': p_identity, 'dependent_s_identity': s_identity}


class Interval:
    def __init__(self, lower, upper=None):
        self.lower = Q(lower)
        self.upper = Q(lower if upper is None else upper)

    def __add__(self, other):
        other = other if isinstance(other, Interval) else Interval(other)
        return Interval(self.lower+other.lower, self.upper+other.upper)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.upper, -self.lower)

    def __sub__(self, other):
        return self + -(other if isinstance(other, Interval) else Interval(other))

    def __rsub__(self, other):
        return Interval(other) + -self

    def __mul__(self, other):
        other = other if isinstance(other, Interval) else Interval(other)
        products = [left*right for left in [self.lower, self.upper]
                    for right in [other.lower, other.upper]]
        return Interval(min(products), max(products))

    __rmul__ = __mul__

    def abs_upper(self):
        return max(abs(self.lower), abs(self.upper))


def matmul(left, right):
    return [[sum((left[i][k]*right[k][j] for k in range(4)), Interval(0))
             for j in range(4)] for i in range(4)]


def subtract(left, right):
    return [[left[i][j]-right[i][j] for j in range(4)] for i in range(4)]


def norm_upper(matrix):
    return max(sum(entry.abs_upper() for entry in row) for row in matrix)


def widened_intervals(text):
    result = []
    for lower, upper in re.findall(r'\[([^,\[\]]+),\s*([^,\[\]]+)\]', text):
        lower, upper = Q(lower), Q(upper)
        slack = Q(1, 10**14)*max(1, abs(lower), abs(upper))
        result.append(Interval(lower-slack, upper+slack))
    return result


def interval_checks():
    seed = (CANDIDATE / 'seed.h').read_text()
    log = LOG.read_text()
    center = [Q(float(value)) for value in re.search(
        r'center_values\[4\] = \{([^}]+)\}', seed).group(1).split(',')]
    rows = ast.literal_eval(re.search(
        r'inverse_values\[4\]\[4\] = (\{.*?\});', seed, re.S).group(1)
        .replace('{', '[').replace('}', ']'))
    preconditioner = [[Interval(Q(float(value))) for value in row] for row in rows]
    entries = widened_intervals(log.split('\nDP {')[1].split('\nDETERMINANT_R')[0])
    assert len(entries) == 16
    derivative = [entries[index:index+4] for index in range(0, 16, 4)]
    residual = widened_intervals(log.split('\nRESIDUAL ')[1].split('\nDP ')[0])
    assert len(residual) == 4
    identity = [[Interval(i == j) for j in range(4)] for i in range(4)]
    correction = subtract(identity, matmul(preconditioner, subtract(derivative, identity)))
    radius = Q(float('1e-8'))
    krawczyk = [Interval(center[i]) - sum(
        (preconditioner[i][j]*residual[j] for j in range(4)), Interval(0)) + sum(
        (correction[i][j]*Interval(-radius, radius) for j in range(4)), Interval(0))
        for i in range(4)]
    inclusion = all(center[i]-radius < box.lower and box.upper < center[i]+radius
                    for i, box in enumerate(krawczyk))
    power = derivative
    for unused in range(3):
        power = matmul(power, derivative)
    determinant = Q(0)
    for permutation in itertools.permutations(range(4)):
        term = Q((-1)**sum(permutation[i] > permutation[j]
                          for i in range(4) for j in range(i+1, 4)))
        for i in range(4):
            term *= preconditioner[i][permutation[i]].lower
        determinant += term
    correction_norm, fourth_power_norm = norm_upper(correction), norm_upper(power)
    displacement = max(max(abs(box.lower-center[i]), abs(box.upper-center[i]))
                       for i, box in enumerate(krawczyk))
    assert inclusion and determinant != 0 and correction_norm < 1 and fourth_power_norm < 1
    rational_values = {
        'determinant_binary_preconditioner': determinant,
        'correction_norm_upper': correction_norm,
        'fourth_power_norm_upper': fourth_power_norm,
        'krawczyk_displacement_upper': displacement,
        'ideal_binary_radius': radius,
    }
    return {
        'strict_krawczyk_inclusion': inclusion,
        'all_logged_derivative_residual_intervals_widened_by':
            '10^-14 * max(1, abs(lower), abs(upper)) on each side',
        'exact_rational_values': {key: str(value) for key, value in rational_values.items()},
        'decimal_display_only': {key: float(value) for key, value in rational_values.items()},
    }


def main():
    global LOG
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--log', type=Path, default=LOG)
    parser.add_argument('--output', type=Path, default=ROOT/'artifacts/validated-cycle/independent-check.json')
    args = parser.parse_args()
    if not __debug__:
        parser.error('Assertions must be enabled')
    LOG = args.log.resolve()
    inputs = [CANDIDATE / name for name in [
        'verify_cycle.cpp', 'seed.h', 'CMakeLists.txt', 'build_and_run.sh',
        'candidate_lemma.md', 'finite_kinase_persistence.md']]
    inputs.append(LOG)
    report = {
        'status': 'passed',
        'scope': 'Independent exact polynomial reconstruction and exact rational '
                 'postprocessing of widened serialized Linux interval bounds. '
                 'Does not independently reimplement CAPD integration or establish priority.',
        'polynomial_checks': polynomial_checks(),
        'interval_checks': interval_checks(),
        'input_sha256': {str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                         for path in inputs},
        'implementation_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    destination = args.output
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
