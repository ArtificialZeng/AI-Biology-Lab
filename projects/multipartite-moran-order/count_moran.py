"""Exact native complete-multipartite mixed-Moran count oracle (stdlib only).

lambda_bd is the probability of Bd; the journal's delta is 1-lambda_bd.
The dense rational solver has a finite resource guard; the theorem in proof.md
has no size cap. This module does not mechanically certify that theorem.
"""
from fractions import Fraction
from itertools import product
from math import prod
import argparse
import json


def require(condition, message):
    if not condition:
        raise ValueError(message)


def rational(value):
    require(not isinstance(value, (float, bool)),
            "Use an integer, rational string, or Fraction; no floating input")
    return Fraction(value)


def parameters(part_sizes, fitness, lambda_bd):
    parts = tuple(part_sizes)
    require(len(parts) >= 2 and all(type(n) is int and n > 0 for n in parts),
            "At least two positive integer part sizes are required")
    r, lam = rational(fitness), rational(lambda_bd)
    require(r > 0 and 0 <= lam <= 1, "Positive fitness and Bd probability in [0,1] required")
    return parts, r, lam


def count_rates(part_sizes, counts, fitness, lambda_bd):
    parts, r, lam = parameters(part_sizes, fitness, lambda_bd)
    k = tuple(counts)
    require(len(k) == len(parts) and all(type(v) is int and 0 <= v <= n
                                       for v, n in zip(k, parts)), "Invalid count vector")
    N, K = sum(parts), sum(k)
    W = N + (r - 1) * K
    degrees = [N - n for n in parts]
    births, deaths = [], []
    for i, n in enumerate(parts):
        t = K - k[i]
        F = degrees[i] + (r - 1) * t
        parent_mut = sum((Fraction(k[j], degrees[j]) for j in range(len(parts)) if j != i),
                         Fraction(0))
        parent_res = sum((Fraction(parts[j] - k[j], degrees[j])
                          for j in range(len(parts)) if j != i), Fraction(0))
        b = (n - k[i]) * (lam * r * parent_mut / W
                         + (1 - lam) * r * t / (N * F))
        a = k[i] * (lam * parent_res / W
                    + (1 - lam) * (degrees[i] - t) / (N * F))
        births.append(b)
        deaths.append(a)
    return births, deaths


def native_count_kernel(part_sizes, fitness, lambda_bd, max_states=256):
    parts, r, lam = parameters(part_sizes, fitness, lambda_bd)
    require(type(max_states) is int and max_states > 0, "Invalid state resource guard")
    require(prod(n + 1 for n in parts) <= max_states,
            "Count-state resource guard exceeded; choose max_states deliberately")
    states = list(product(*(range(n + 1) for n in parts)))
    kernel = {}
    for k in states:
        b, a = count_rates(parts, k, r, lam)
        row = {}
        for i in range(len(parts)):
            for sign, rate in ((1, b[i]), (-1, a[i])):
                if rate:
                    target = list(k)
                    target[i] += sign
                    target = tuple(target)
                    require(all(0 <= v <= n for v, n in zip(target, parts)),
                            "Jump outside count domain")
                    row[target] = row.get(target, Fraction(0)) + rate
        diagonal = 1 - sum(row.values(), Fraction(0))
        require(diagonal >= 0, "Negative native holding probability")
        row[k] = row.get(k, Fraction(0)) + diagonal
        kernel[k] = row
    return kernel


def harmonic_values(kernel, zero, full):
    """Exact absorbing harmonic solution with explicit kernel/residual checks."""
    states = list(kernel)
    require(zero in kernel and full in kernel and zero != full, "Invalid boundaries")
    for s, row in kernel.items():
        require(all(t in kernel and isinstance(v, Fraction) and v >= 0
                    for t, v in row.items()), "Invalid transition row")
        require(sum(row.values(), Fraction(0)) == 1, "Transition row is not normalized")
    require(kernel[zero].get(zero) == 1 and kernel[full].get(full) == 1,
            "Boundaries are not absorbing")
    trans = [s for s in states if s not in (zero, full)]
    index = {s: i for i, s in enumerate(trans)}
    M = len(trans)
    aug = [[Fraction(int(i == j)) - kernel[s].get(t, Fraction(0))
            for j, t in enumerate(trans)] + [kernel[s].get(full, Fraction(0))]
           for i, s in enumerate(trans)]
    for col in range(M):
        pivot = next((i for i in range(col, M) if aug[i][col]), None)
        require(pivot is not None, "Singular transient system")
        aug[col], aug[pivot] = aug[pivot], aug[col]
        for row in range(col + 1, M):
            if not aug[row][col]:
                continue
            factor = aug[row][col] / aug[col][col]
            aug[row][col] = Fraction(0)
            for j in range(col + 1, M + 1):
                aug[row][j] -= factor * aug[col][j]
    x = [Fraction(0)] * M
    for i in range(M - 1, -1, -1):
        x[i] = (aug[i][M] - sum((aug[i][j] * x[j] for j in range(i + 1, M)),
                               Fraction(0))) / aug[i][i]
    values = {zero: Fraction(0), full: Fraction(1)}
    values.update({s: x[index[s]] for s in trans})
    for s, row in kernel.items():
        require(0 <= values[s] <= 1, "Harmonic probability outside [0,1]")
        require(values[s] == sum((v * values[t] for t, v in row.items()), Fraction(0)),
                "Nonzero exact harmonic residual")
    return values


def fixation_probabilities(part_sizes, fitness, lambda_bd, max_states=256):
    parts, r, lam = parameters(part_sizes, fitness, lambda_bd)
    kernel = native_count_kernel(parts, r, lam, max_states=max_states)
    return harmonic_values(kernel, tuple(0 for _ in parts), parts)


def universal_count_attractivity_criterion(part_sizes, lambda_bd):
    """Analytic criterion in proof.md; it is universal in r, not pointwise in r."""
    parts, _, lam = parameters(part_sizes, 1, lambda_bd)
    return lam == 0 or len(parts) == 2 or len(set(parts)) == 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parts", required=True, help="Comma-separated positive sizes")
    parser.add_argument("--fitness", required=True, help="Positive rational mutant fitness")
    parser.add_argument("--bd", required=True, help="Rational Bd probability in [0,1]")
    parser.add_argument("--initial", required=True, help="Comma-separated mutant counts")
    parser.add_argument("--max-states", type=int, default=256)
    args = parser.parse_args()
    parts = tuple(int(v) for v in args.parts.split(","))
    initial = tuple(int(v) for v in args.initial.split(","))
    parts, r, lam = parameters(parts, args.fitness, args.bd)
    count_rates(parts, initial, r, lam)
    values = fixation_probabilities(parts, r, lam, max_states=args.max_states)
    print(json.dumps({"part_sizes": parts, "initial_counts": initial,
                      "fitness": str(r), "lambda_Bd": str(lam),
                      "fixation_probability_exact": str(values[initial]),
                      "count_states": len(values), "binary_vertex_states": 2 ** sum(parts),
                      "all_fitness_native_count_attractivity_criterion":
                          universal_count_attractivity_criterion(parts, lam),
                      "theorem_mechanically_certified": False}, sort_keys=True))


if __name__ == "__main__":
    main()
