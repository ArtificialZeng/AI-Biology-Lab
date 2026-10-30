"""Exact finite checks supporting, not replacing, the analytical proof."""
from fractions import Fraction as Q
from itertools import combinations, product
from pathlib import Path
import json
import sys

H = ((1, 1, 1, 1), (1, 1, -1, -1),
     (1, -1, 1, -1), (1, -1, -1, 1))

def require(test, message):
    if not test:
        raise ValueError(message)

def inverse(eigen):
    return tuple(sum(Q(H[g][x]) * eigen[g] for g in range(4))/4
                 for x in range(4))

def transform(p):
    return tuple(sum(Q(H[g][x])*p[x] for x in range(4)) for g in range(4))

def convolve(p, q):
    return tuple(sum(p[x]*q[x ^ z] for x in range(4)) for z in range(4))

def jc(t):
    return inverse((Q(1), t, t, t))

def validate(p, model, positive=False):
    require(sum(p) == 1 and all(v >= 0 for v in p), "Not stochastic")
    if positive:
        require(all(v > 0 for v in p), "Not strictly positive")
    a = transform(p)
    require(a[0] == 1 and all(0 < v < 1 for v in a[1:]), "Spectrum outside domain")
    if model == "JC":
        require(a[1] == a[2] == a[3], "Wrong JC relation")
    if model == "K2P":
        require(a[1] == a[3], "Wrong K2P relation")
    return a

def direct_diamond_joint(T, J, da, db):
    # Direct ancestral-variable summation, without the displayed-path polynomial.
    out = {}
    for a, b in product(range(4), repeat=2):
        total = Q(0)
        for r, u, v in product(range(4), repeat=3):
            pa = da*T[u ^ a] + (1-da)*T[v ^ a]
            pb = db*J[u ^ b] + (1-db)*J[v ^ b]
            total += J[r ^ u]*J[r ^ v]*pa*pb/4
        out[a, b] = total
    return out

def serial(value):
    if isinstance(value, Q):
        return str(value)
    if isinstance(value, dict):
        return {str(k): serial(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [serial(v) for v in value]
    return value

def compute():
    support_sets = [set(s) for k in (3, 4) for s in combinations(range(4), k)]
    support_checks = 0
    for s, t in product(support_sets, repeat=2):
        require({x ^ y for x in s for y in t} == set(range(4)), "Support lemma failed")
        support_checks += 1

    nonidentity_supports = [set(s) for k in (2, 3, 4)
                           for s in combinations(range(4), k) if 0 in s]
    arbitrary_supports = [set(s) for k in (1, 2, 3, 4)
                          for s in combinations(range(4), k)]
    broader_support_checks = 0
    for s, t, u in product(nonidentity_supports, nonidentity_supports, arbitrary_supports):
        actual = {x ^ y ^ z for x in s for y in t for z in u}
        require(len(actual) in (2, 4), "Three-support attachment obstruction failed")
        broader_support_checks += 1

    boundary = (Q(3, 5), Q(1, 5), Q(0), Q(1, 5))
    eig = validate(boundary, "K2P")
    require(eig == (1, Q(3, 5), Q(1, 5), Q(3, 5)), "Wrong boundary spectrum")
    require(inverse(eig) == boundary, "Hadamard inversion")
    require(convolve(jc(Q(1, 2)), jc(Q(1, 2))) == jc(Q(1, 4)), "Rational root split")
    pendants = (jc(Q(1, 4)), jc(Q(1, 2)), jc(Q(1, 2)), jc(Q(1, 2)))
    law = {}
    for leaves in product(range(4), repeat=4):
        x1, x2, x3, x4 = leaves
        law[leaves] = sum(
            boundary[s ^ t]*pendants[0][s ^ x1]*pendants[1][s ^ x2]
            *pendants[2][t ^ x3]*pendants[3][t ^ x4]/4
            for s, t in product(range(4), repeat=2))
    require(sum(law.values()) == 1 and min(law.values()) > 0, "Invalid target distribution")
    recovered = []
    pair_records = []
    for g in range(1, 4):
        correlations = {}
        for i, j in combinations(range(4), 2):
            correlations[i+1, j+1] = sum(p*H[g][x[i] ^ x[j]] for x, p in law.items())
        require(correlations[1, 2] == Q(1, 8), "Within-left correlation")
        require(correlations[3, 4] == Q(1, 4), "Within-right correlation")
        square = correlations[1, 3]*correlations[2, 4]/(correlations[1, 2]*correlations[3, 4])
        require(square == eig[g]**2 and eig[g] > 0, "Internal-edge invariant")
        recovered.append(eig[g])
        pair_records.append({"character": g, "correlations": correlations,
                             "internal_eigenvalue_squared": square,
                             "unique_positive_root": eig[g]})
    require(inverse((Q(1), *recovered)) == boundary, "Recovered forbidden boundary")

    targets = [("JC", jc(Q(1, 5))),
               ("K2P", convolve(boundary, jc(Q(1, 2)))),
               ("K3P", inverse((Q(1), Q(1, 5), Q(1, 3), Q(2, 5))))]
    constructions = []
    for model, p in targets:
        a = validate(p, model, positive=True)
        for da, db in [(Q(1, 2), Q(1, 2)), (Q(1, 3), Q(2, 5)), (Q(4, 5), Q(1, 7))]:
            weights = (da*db, da*(1-db), (1-da)*db, (1-da)*(1-db))
            lengths = (2, 4, 4, 2)
            L = max(lengths)
            F = max(max(a[1:]), 1-4*min(p))
            t = 1-(1-F)/(2*L)
            f = sum(w*t**(ell-1) for w, ell in zip(weights, lengths))
            require(0 < F < f < 1 and 0 < t < 1, "Bad construction bound")
            T = tuple((v-(1-f)/4)/f for v in p)
            validate(T, model, positive=True)
            require(convolve(T, jc(f)) == p, "Compensating kernel identity")
            actual = direct_diamond_joint(T, jc(t), da, db)
            require(all(actual[x, y] == p[x ^ y]/4 for x, y in actual),
                    "Direct DAG summation disagrees")
            constructions.append({"model": model, "target": p, "inheritance": [da, db],
                                  "L": L, "F": F, "t": t, "f": f, "T": T,
                                  "direct_DAG_joint_cells_checked": len(actual)})
    return serial({
        "scope": "Exact finite kernel, target-membership, observable-invariant and explicit-construction checks; general statements depend on proof.md",
        "character_table": H,
        "support_pair_checks": support_checks,
        "broader_attachment_support_checks": broader_support_checks,
        "boundary_probability_row": boundary, "boundary_eigenvalues": eig,
        "target_rational_root_split": "J_(1/2) J_(1/2) = J_(1/4)",
        "quartet_target_probabilities": {"".join(map(str, x)): p for x, p in law.items()},
        "target_probability_sum": sum(law.values()),
        "target_minimum_leaf_probability": min(law.values()),
        "observable_internal_edge_reconstruction": pair_records,
        "rooted_joint_constructions": constructions,
        "number_of_constructed_cases": len(constructions),
        "no_claims": ["No finite search proves general nonmembership", "No human or trusted-runner certification", "No arbitrary-root-law or CTMC result"]})

def main():
    actual = compute()
    if sys.argv[1:] == ["compute"]:
        print(json.dumps(actual, indent=2, sort_keys=True))
    elif sys.argv[1:] in ([], ["verify"]):
        expected = json.loads((Path(__file__).parent/"certificate.json").read_text())
        require(actual == expected, "Certificate differs")
        print(json.dumps({"exact_verification": True, "target_leaf_patterns": 256,
                          "constructed_cases": 9, "support_pairs": 25,
                          "broader_attachment_support_cases": 735}, sort_keys=True))
    else:
        raise SystemExit("Use verify or compute")

if __name__ == "__main__":
    main()
