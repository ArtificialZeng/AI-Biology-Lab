#!/usr/bin/env python3
"""Exact finite checks for the frozen Conjecture 43 candidate; no probability proof.

Python 3 standard library only. Run from any working directory. Inputs are read
beside this script; the sole written output is evidence/certificate.json there.
No assertions, simulation, random sampling, floating point, or external imports.
"""

from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path


class CheckFailure(RuntimeError):
    pass


def require(condition, message):
    if not condition:
        raise CheckFailure(message)


def rational(value):
    require(type(value) in (int, str), "rates must be integers or rational strings")
    try:
        return F(value)
    except (ValueError, ZeroDivisionError) as error:
        raise CheckFailure("invalid rational rate") from error


def vector(value):
    require(isinstance(value, list) and len(value) == 3, "three coordinates required")
    require(all(type(x) is int and x >= 0 for x in value), "stoichiometries must be nonnegative integers")
    return value


def exact(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: exact(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [exact(item) for item in value]
    return value


def row_times_matrix(row, matrix):
    return [sum((row[i] * matrix[i][j] for i in range(3)), F(0)) for j in range(3)]


def matrix_times_column(matrix, column):
    return [sum((matrix[i][j] * column[j] for j in range(3)), F(0)) for i in range(3)]


def analyze(model):
    """Derive all finite facts from input arrows; fail if frozen properties change."""
    require(model.get("species") == ["A", "B", "C"], "species/order must be A,B,C")
    require(model.get("inflows") is False, "candidate has no inflows")
    require(vector(model.get("initial_state")) == [0, 1, 0], "initial population must be one B")
    require(model.get("initial_environment") == "either", "initial environment must be either state")
    require(model.get("environment") == "independent symmetric1<->2 ratekappa", "environment declaration changed")
    require(model.get("kappa") == "every strictly positive value", "kappa quantifier changed")
    reactions = model.get("reactions")
    require(isinstance(reactions, dict) and set(reactions) == {"1", "2"}, "two environments required")
    matrices, derived, lineage, propensity, upward = {}, {}, {}, {}, {}
    max_jump, total_reactions = 0, 0
    for environment in ("1", "2"):
        arrows = reactions[environment]
        require(isinstance(arrows, list) and len(arrows) == 3, "three arrows per environment required")
        matrix = [[F(0) for _ in range(3)] for _ in range(3)]
        total_coeff, up_coeff = [F(0)] * 3, [F(0)] * 3
        details, source_indices, b_arrows = [], [], []
        local_lineage = {}
        for arrow in arrows:
            source, target = vector(arrow.get("source")), vector(arrow.get("product"))
            require(sum(source) == 1, "every source must be unimolecular (no zero source)")
            rate = rational(arrow.get("rate"))
            require(rate > 0, "every existing reaction rate must be strictly positive")
            j = source.index(1)
            source_indices.append(j)
            delta = [target[i] - source[i] for i in range(3)]
            for i in range(3):
                matrix[i][j] += rate * delta[i]
            total_coeff[j] += rate
            jump = sum(target) - 1
            max_jump = max(max_jump, abs(jump))
            if jump > 0:
                up_coeff[j] += rate
            require(target[1] == 0, "B must never be produced")
            if j == 1:
                b_arrows.append((target, rate))
            else:
                require(all(target[i] == 0 for i in range(3) if i != j), "A/C offspring must retain lineage type")
                require(target[j] in (0, 2), "A/C events must be death or binary splitting")
                local_lineage[model["species"][j]] = {
                    "birth_rate": rate if target[j] == 2 else F(0),
                    "death_rate": rate if target[j] == 0 else F(0),
                    "net_per_capita_growth": rate * (target[j] - 1),
                }
            details.append({"source": source, "product": target, "rate": rate,
                            "source_species": model["species"][j], "delta": delta,
                            "drift_column_contribution": [rate * d for d in delta]})
        require(sorted(source_indices) == [0, 1, 2], "exactly one arrow for each source species per environment")
        require(len(b_arrows) == 1 and b_arrows[0][1] == 1, "sole B consumption rate must be one in each environment")
        expected_target = [1, 0, 0] if environment == "1" else [0, 0, 1]
        require(b_arrows[0][0] == expected_target, "B conversion product changed")
        matrices[environment], derived[environment] = matrix, details
        lineage[environment], propensity[environment], upward[environment] = local_lineage, total_coeff, up_coeff
        total_reactions += len(arrows)
    expected = {"1": [[1, 1, 0], [0, -1, 0], [0, 0, -2]],
                "2": [[-2, 0, 0], [0, -1, 0], [0, 1, 1]]}
    require(matrices == expected, "derived drift matrices differ from frozen candidate")
    weights = {"1": [2, 1, 0], "2": [0, 1, 2]}
    weighted, supports = {}, {}
    for environment, weight in weights.items():
        weighted[environment] = row_times_matrix(weight, matrices[environment])
        supports[environment] = [i for i, coordinate in enumerate(weight) if coordinate > 0]
        require(weighted[environment] == weight, "left-vector eigenvalue-one identity failed")
        require(all(weighted[environment][i] > 0 for i in supports[environment]), "required supported drift must be positive")
    shared = sorted(set(supports["1"]) & set(supports["2"]))
    require(shared == [1], "common support must be exactly B")
    initial_support = [i for i, x in enumerate(model["initial_state"]) if x > 0]
    require(bool(set(shared) & set(initial_support)), "initial support misses shared support")
    for species, good, bad in (("A", "1", "2"), ("C", "2", "1")):
        require(lineage[good][species]["birth_rate"] == 1 and lineage[good][species]["death_rate"] == 0,
                "matching environment must give birth rate one and death rate zero")
        require(lineage[bad][species]["birth_rate"] == 0 and lineage[bad][species]["death_rate"] == 2,
                "opposite environment must give birth rate zero and death rate two")
    q_coefficient = [[-1, 1], [1, -1]]
    stationary = [F(1, 2), F(1, 2)]
    require(all(sum(stationary[i] * q_coefficient[i][j] for i in range(2)) == 0 for j in range(2)),
            "uniform stationary distribution identity failed")
    averaged = {species: sum(stationary[int(env) - 1] * lineage[env][species]["net_per_capita_growth"]
                            for env in ("1", "2")) for species in ("A", "C")}
    require(all(value == F(-1, 2) for value in averaged.values()), "occupation-weighted exponent must be -1/2")
    require(max_jump <= 1, "population jump size must be at most one")
    require(all(0 <= c <= 2 for coeff in propensity.values() for c in coeff), "molecular total jump rate bound 2N failed")
    require(all(0 <= c <= 1 for coeff in upward.values() for c in coeff), "population-increase rate bound N failed")
    # Exhaustive exact small tuples cross-check reaction propensities against the
    # independently assembled matrix and both coefficientwise bounds.
    tuple_count = 0
    for state in product(range(5), repeat=3):
        for environment in ("1", "2"):
            drift = [F(0)] * 3
            total, increasing = F(0), F(0)
            for arrow in derived[environment]:
                j = model["species"].index(arrow["source_species"])
                hazard = arrow["rate"] * state[j]
                total += hazard
                if sum(arrow["delta"]) > 0:
                    increasing += hazard
                for i in range(3):
                    drift[i] += hazard * arrow["delta"][i]
            require(drift == matrix_times_column(matrices[environment], state), "exact tuple drift mismatch")
            require(total <= 2 * sum(state) and increasing <= sum(state), "exact tuple rate bound failed")
            tuple_count += 1
    return {
        "checks_passed": True, "reaction_count": total_reactions,
        "structural_checks": {"all_existing_rate_constants_strictly_positive": True,
                              "all_sources_unimolecular": True, "inflows_absent": True,
                              "initially_one_B_in_either_environment": True,
                              "all_A_and_C_offspring_retain_lineage_type": True},
        "convention": "column population x; drift M_i x; left weights v_i M_i",
        "arrows_and_derived_contributions": derived, "drift_matrices": matrices,
        "left_vectors": weights, "left_products": weighted,
        "support_species": {env: [model["species"][i] for i in indices] for env, indices in supports.items()},
        "shared_support": [model["species"][i] for i in shared],
        "initial_support": [model["species"][i] for i in initial_support],
        "B_structure": {"never_produced": True, "sole_consumption_rate_each_environment": "1",
                        "initially_only_particle": True, "no_other_particles_before_conversion": True,
                        "conversion_products": {"1": "A", "2": "C"}},
        "lineage_rates": lineage,
        "environment_algebra": {"generator": "kappa times Q_coefficient, kappa > 0",
                                "Q_coefficient": q_coefficient, "stationary_distribution": stationary,
                                "off_diagonal_coefficients_positive": True},
        "growth_algebra": {"formula": "1 * good_environment_time - 2 * bad_environment_time",
                           "ergodic_occupation_fractions_used_by_human_proof": stationary,
                           "occupation_weighted_per_capita_exponents": averaged,
                           "ergodic_theorem_proved_by_checker": False},
        "population_domination_algebra": {"total_molecular_jump_rate_coefficients": propensity,
                                          "total_molecular_jump_rate_bound": "2 * N",
                                          "including_environment_jump_rate_bound": "2 * N + kappa",
                                          "population_increase_rate_coefficients": upward,
                                          "population_increase_rate_bound": "N",
                                          "max_absolute_population_jump": max_jump,
                                          "human_argument": "linear Yule pure-birth domination; environment has constant rate kappa"},
        "exact_small_tuple_cross_checks": {"states": "{0,1,2,3,4}^3 in both environments",
                                          "state_environment_pairs": tuple_count, "all_passed": True},
    }


def negative_controls(model):
    cases = []
    mutations = (
        ("zero_existing_rate", "every existing reaction rate must be strictly positive", ("rate", "1", 0, 0)),
        ("altered_positive_birth_rate", "derived drift matrices differ from frozen candidate", ("rate", "1", 0, 2)),
        ("b_consumption_rate_changed", "sole B consumption rate must be one", ("rate", "1", 1, 2)),
        ("zero_source_inflow_arrow", "every source must be unimolecular", ("source", "1", 0, [0, 0, 0])),
        ("b_regeneration", "B must never be produced", ("product", "1", 0, [1, 1, 0])),
        ("lineage_type_changes", "A/C offspring must retain lineage type", ("product", "1", 0, [1, 0, 1])),
    )
    for name, expected, (field, environment, index, value) in mutations:
        altered = deepcopy(model)
        altered["reactions"][environment][index][field] = value
        try:
            analyze(altered)
        except CheckFailure as error:
            require(expected in str(error), "negative control failed for an unintended reason: " + name)
            cases.append({"name": name, "rejected": True, "reason": str(error)})
        else:
            raise CheckFailure("altered model unexpectedly passed: " + name)
    return cases


def main():
    directory = Path(__file__).resolve().parent
    raw = (directory / "model.json").read_bytes()
    model = json.loads(raw)
    certificate = {
        "schema": "conj43-owner-exact-finite-checks-v1",
        "status": "finite_algebra_and_structure_verified_only",
        "input_sha256": sha256(raw).hexdigest(),
        "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "dependencies": "Python 3 standard library only; fractions.Fraction exact arithmetic",
        "scope_limit": "This certificate does not prove any probability theorem, nonexplosion, ergodic theorem, extinction, source interpretation, novelty, or acceptance. Those require the human-readable proof and independent review.",
        "probability_theorem_proved": False,
        "formalization": False,
        "scientific_admission_or_acceptance": False,
        "finite_checks": analyze(model),
        "meaningful_negative_controls": negative_controls(model),
    }
    rendered = json.dumps(exact(certificate), indent=2, sort_keys=True) + "\n"
    evidence = directory / "evidence"
    evidence.mkdir(exist_ok=True)
    (evidence / "certificate.json").write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
