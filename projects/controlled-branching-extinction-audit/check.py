"""Exact algebra checks, deliberately separate from infinite-state analysis."""

from fractions import Fraction as F
import json
from pathlib import Path


def require(condition, message):
    if not condition:
        raise ValueError(message)


def multiply(p, q):
    out = [F(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return out


def main():
    data = json.loads(Path(__file__).with_name("witness.json").read_text())
    require(data["schema_version"] == 1, "unsupported schema")
    require(data["dimension"] == 2, "dimension")
    require(data["individual_rates"] == [1, 1], "rates")
    require(data["type2_offspring"] == [1, 0], "type 2 law")
    require(data["type1_support"] == [[0, 0], [2, 0]], "type 1 support")

    example = data["inaccessible_example"]
    states = [(1, 0), (2, 0), (3, 0)]
    require(example["closed_states"] == [list(s) for s in states], "states")
    require(example["carrying_capacity"] == 2, "capacity")
    require(example["type1_mean"] == "2 if 0<=z<=1; 3-z if 1<z<3; 0 if z>=3", "mean law")
    q = [[F(0) for _ in states] for _ in states]
    for row, (n, _) in enumerate(states):
        mean = F(2) if n == 1 else F(3 - n)
        for offspring, probability in [(0, 1 - mean / 2), (2, mean / 2)]:
            if probability:
                target = (n + offspring - 1, 0)
                require(target in states, "class not closed")
                col = states.index(target)
                q[row][col] += n * probability
                q[row][row] -= n * probability
    require(q == example["generator"], "generator reconstruction")
    require(example["limiting_matrix"] == [[-1, 0], [1, -1]], "Hurwitz limit")
    require(all(sum(row) == 0 for row in q), "row sums")
    require(all(q[i][j] >= 0 for i in range(3) for j in range(3) if i != j), "jump rates")

    example = data["critical_limit_example"]
    require(example["carrying_capacity"] == 1, "critical capacity")
    require(example["type1_mean"] == "1+(1-z)/(1+z)^3", "critical law")
    require(example["limiting_matrix"] == [[0, 0], [1, -1]], "critical limit")
    require(example["a_n"] == "(n-1)/(n+1)^3", "a law")
    require(F(example["sum_a_upper_bound"]) == F(1, 27) + F(1, 3), "series bound")
    require(F(example["all_finite_products_lower_bound"]) == 1 - 2 * F(10, 27), "product bound")
    require(example["occupation_lower_bound"] == "7/(27*n)", "occupation bound")

    cube_z_plus_one = multiply(multiply([1, 1], [1, 1]), [1, 1])
    require(cube_z_plus_one == [1, 3, 3, 1], "denominator expansion")
    # For z>=0, (2-m(z))*(1+z)^3 = z^3+3z^2+4z.
    upper_numerator = cube_z_plus_one[:]
    upper_numerator[0] -= 1
    upper_numerator[1] += 1
    require(upper_numerator == [0, 4, 3, 1], "upper probability certificate")
    # For z=t+1>=1, 4*(m(z)-3/4)*(1+z)^3
    # is (t+2)^3-4t = t^3+6t^2+8t+8, strictly positive.
    shifted_cube = multiply(multiply([2, 1], [2, 1]), [2, 1])
    shifted_cube[1] -= 4
    require(shifted_cube == [8, 8, 6, 1], "lower probability certificate")
    require(all(c > 0 for c in shifted_cube), "positive coefficients")
    # For 0<=a<=1/4, (1-a)/(1+a) - (1-2a) = 2*a^2/(1+a).
    product = multiply([1, -2], [1, 1])
    require([F(1) - product[0], F(-1) - product[1], -product[2]] == [0, 0, 2], "product inequality")

    # Finite regression checks of the occupation recurrence, not a proof of its
    # infinite sum or of eventual absorption. Those arguments are in proof.md.
    ratio_product = F(1)
    previous_upcrossings = None
    for n in range(1, 101):
        a = F(n - 1, (n + 1) ** 3)
        require(0 <= a <= F(1, 4), "offspring probability range")
        birth, death = F(n, 2) * (1 - a), F(n, 2) * (1 + a)
        occupation = ratio_product / death
        if n == 1:
            require(death * occupation == 1, "absorbing flux")
        else:
            require(death * occupation == previous_upcrossings, "edge crossing identity")
        require(occupation >= F(7, 27 * n), "finite occupation lower bound")
        previous_upcrossings = birth * occupation
        ratio_product *= birth / death
    print(json.dumps({
        "status": "exact-algebra-checks-passed",
        "closed_class_states": 3,
        "finite_occupation_regressions": 100,
        "infinite_state_theorems_verified_by_script": False,
        "requires_analytic_proof": "proof.md",
    }, indent=2))


if __name__ == "__main__":
    main()
