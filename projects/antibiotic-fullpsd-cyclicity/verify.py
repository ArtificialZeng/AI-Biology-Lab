"""Read-only exact identities and finite native-map corroboration.

The analytic all-parameter theorem is in proof.md. This script is not a
mechanical proof of that theorem. Normal and optimized Python run the same
explicit checks; no assert statements are used.
"""
from pathlib import Path
from fractions import Fraction as Q
import cmath
import hashlib
import itertools
import json
import math
import sympy as sp

HERE = Path(__file__).resolve().parent
MEMBERS = {
    "README.md", "proof.md", "sources.md", "examples.json",
    "certificate.json", "verify.py", "requirements.txt", "LICENSE",
}
WORDS = {"".join(t) for t in itertools.permutations("PSD")}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def archive_check():
    manifest = json.loads((HERE / "archive_manifest.json").read_text())
    require(manifest.get("schema") == 1, "manifest schema")
    require(set(manifest.get("files", {})) == MEMBERS, "manifest membership")
    observed = {p.name for p in HERE.iterdir() if p.is_file()}
    require(observed == MEMBERS | {"archive_manifest.json"}, "archive file closure")
    for name, expected in manifest["files"].items():
        path = HERE / name
        require(not path.is_symlink(), "archive symlink")
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        require(actual == expected, "archive hash: " + name)


def cyclic(a, b):
    return b in (a[1:] + a[:1], a[2:] + a[:2])


def classify(support):
    require(bool(support), "empty support")
    require(len(set(support)) == len(support), "duplicate support word")
    require(all(c in WORDS for c in support), "not a full-PSD column")
    require(not any(cyclic(a, b) for a, b in itertools.combinations(support, 2)),
            "input has a cyclic pair")
    fixed = [(i, support[0][i]) for i in range(3)
             if all(c[i] == support[0][i] for c in support)]
    require(len(support) <= 2 and bool(fixed), "support reduction")
    if len(support) == 2:
        require(len(fixed) == 1, "two-type fixed role")
    return fixed


def check_supports():
    words = sorted(WORDS)
    records = []
    for size in range(1, 7):
        for support in itertools.combinations(words, size):
            if not any(cyclic(a, b) for a, b in itertools.combinations(support, 2)):
                records.append((support, classify(support)))
    roles = {role: sum(len(s) == 2 and fixed[0][1] == role
                       for s, fixed in records) for role in "PSD"}
    counts = {
        "all_nonempty": 2 ** 6 - 1,
        "no_cyclic_pair": len(records),
        "singleton": sum(len(s) == 1 for s, _ in records),
        "two_type": sum(len(s) == 2 for s, _ in records),
        "two_type_common_roles": roles,
    }
    require(counts == {"all_nonempty": 63, "no_cyclic_pair": 15,
                      "singleton": 6, "two_type": 9,
                      "two_type_common_roles": {"P": 3, "S": 3, "D": 3}},
            "complete support counts")
    return counts, {tuple(s) for s, _ in records}


def check_symbolic_identities():
    # Start with constrained native channel survival, not a predeclared G.
    a, b, y1, y2, m, n = sp.symbols("a b y1 y2 m n", positive=True)
    x0 = 1 - y1 - y2
    h = sp.exp(-a * x0)
    k = sp.exp(-b * x0)
    h1, h2 = sp.exp(-a * y2), sp.exp(-a * y1)
    uP = [1 - sp.exp(-b * y2) * (1 - h),
          1 - sp.exp(-b * y1) * (1 - h)]
    uD = [1 - k * (1 - h1), 1 - k * (1 - h2)]
    expected = {
        "P": -m * n * sp.exp(-b * (y1 + y2)) * b * (1 - h)
             * (2 * a * h + b * (1 - h)) / (uP[0] * uP[1]),
        "D": -m * n * k ** 2 * a
             * (b * ((1 - h1) * h2 + (1 - h2) * h1) + a * h1 * h2)
             / (uD[0] * uD[1]),
    }
    for role, u in (("P", uP), ("D", uD)):
        G = sp.Matrix([[mult * sp.diff(sp.log(u[i]), y)
                        for y in (y1, y2)]
                       for i, mult in enumerate((m, n))])
        require(sp.simplify(G.det() - expected[role]) == 0,
                "native determinant identity " + role)
    B = sp.diag(y1, y2) - sp.Matrix([y1, y2]) * sp.Matrix([[y1, y2]])
    require(sp.expand(B.det() - x0 * y1 * y2) == 0, "chart determinant")
    A, Bv, C, D = sp.symbols("A B C D")
    require(sp.expand(A * D - (A + Bv) * (D + C)
                      + A * C + Bv * D + Bv * C) == 0,
            "strict-sign polynomial identity")
    return 4


def survival(x, columns, kp, kd):
    q = [1 + 0j, 1 + 0j, 1 + 0j]
    for col in columns:
        produced = sum(x[i] for i in range(3) if col[i] == "P")
        degraded = sum(x[i] for i in range(3) if col[i] == "D")
        hit = cmath.exp(-kd * degraded) * (1 - cmath.exp(-kp * produced))
        for i in range(3):
            if col[i] == "S":
                q[i] *= 1 - hit
    return q


def native_map(y, columns, kp, kd, growth):
    x = [1 - y[0] - y[1], y[0], y[1]]
    q = survival(x, columns, kp, kd)
    weighted = [x[i] * growth[i] * q[i] for i in range(3)]
    den = sum(weighted)
    return [weighted[1] / den, weighted[2] / den]


def direct_jacobian(y, columns, kp, kd, growth):
    J = [[0.0, 0.0], [0.0, 0.0]]
    for j in range(2):
        yy = list(map(complex, y))
        yy[j] += 1e-25j
        value = native_map(yy, columns, kp, kd, growth)
        for i in range(2):
            J[i][j] = value[i].imag / 1e-25
    return J


def multipliers(J):
    trace = J[0][0] + J[1][1]
    det = J[0][0] * J[1][1] - J[0][1] * J[1][0]
    root = cmath.sqrt(trace * trace - 4 * det)
    return [(trace + root) / 2, (trace - root) / 2]


def check_floating_fixtures(cases, all_supports):
    require(len(cases) == 60, "floating fixture count")
    represented = set()
    for case in cases:
        support = case["support"]
        fixed = classify(support)
        represented.add(tuple(sorted(support)))
        counts = case["multiplicities"]
        require(len(counts) == len(support), "multiplicity dimensions")
        require(all(type(n) is int and n > 0 for n in counts), "multiplicities")
        columns = [c for c, n in zip(support, counts) for _ in range(n)]
        x = case["x"]
        require(len(x) == 3 and all(math.isfinite(v) and v > 0 for v in x),
                "interior fixture")
        require(abs(sum(x) - 1) < 1e-12, "fixture simplex")
        kp, kd = case["KP"], case["KD"]
        require(all(math.isfinite(v) and v > 0 for v in (kp, kd)),
                "positive finite strengths")
        q = survival(x, columns, kp, kd)
        require(all(v.real > 0 and abs(v.imag) < 1e-15 for v in q),
                "positive survival")
        # Freeze rates at this point BEFORE differentiating the map.
        # These positive-rate fixtures are not all source-cost parameter fits.
        growth = [1 / v.real for v in q]
        got = native_map(x[1:], columns, kp, kd, growth)
        require(max(abs(got[i] - x[i + 1]) for i in range(2)) < 1e-12,
                "placed fixed point")
        eig = multipliers(direct_jacobian(x[1:], columns, kp, kd, growth))
        invariant = len(support) == 1 or fixed[0][1] == "S"
        if invariant:
            require(min(abs(v - 1) for v in eig) < 1e-8,
                    "immune-ratio unit multiplier")
        else:
            require(max(abs(v.imag) for v in eig) < 1e-8
                    and max(v.real for v in eig) > 1 + 1e-9,
                    "common-role unstable multiplier")
    require(represented == all_supports, "fixture support coverage")
    return len(cases)


def check_ordered_cost_examples(examples):
    require(len(examples) == 2, "ordered-cost example count")
    roles = set()
    numerical = []
    for example in examples:
        cols = example["columns"]
        fixed = classify(cols)
        role = example["common_role"]
        require(fixed == [(0, role)] and role in ("P", "D"),
                "ordered-cost common role")
        roles.add(role)
        g, cp, cd, cr = (Q(example[k]) for k in ("g", "cP", "cD", "cR"))
        require(g > 0 and 0 < cr < cd < cp, "strict cost order")
        growth = [g - cp * sum(c[i] == "P" for c in cols)
                  - cd * sum(c[i] == "D" for c in cols) for i in range(3)]
        require(min(growth) > 0, "positive cost-derived growth")
        parg = Q(example["KP_over_3_log_argument"])
        darg = Q(example["KD_over_3_log_argument"])
        require(parg > 1 and darg > 1, "positive logarithmic strengths")
        # x=(1/3,1/3,1/3), so exponentials are EXACT rationals.
        hp, hd = 1 / parg, 1 / darg
        q = [Q(1), Q(1), Q(1)]
        for col in cols:
            for i in range(3):
                if col[i] == "S":
                    q[i] *= 1 - hd * (1 - hp)
        f = [r * s for r, s in zip(growth, q)]
        require(f[0] == f[1] == f[2] == Q(example["expected_fitness"]),
                "exact ordered-cost fixed point")
        kp, kd = 3 * math.log(float(parg)), 3 * math.log(float(darg))
        y = [1 / 3, 1 / 3]
        numeric_growth = list(map(float, growth))
        got = native_map(y, cols, kp, kd, numeric_growth)
        require(max(abs(got[i] - y[i]) for i in range(2)) < 1e-12,
                "ordered-cost native replay")
        eig = multipliers(direct_jacobian(y, cols, kp, kd, numeric_growth))
        require(max(abs(v.imag) for v in eig) < 1e-8
                and max(v.real for v in eig) > 1 + 1e-9,
                "ordered-cost unstable multiplier")
        numerical.append({"common_role": role,
                          "multipliers_approximately": [v.real for v in eig]})
    require(roles == {"P", "D"}, "ordered-cost role coverage")
    return len(examples), numerical


def prior_comparison():
    # PP/SD edge of columns PSD,PDS: no degrader for the antibiotic
    # that the SD strain is sensitive to. KP=log2 gives survival1/2.
    g, cp, cd = Q(1), Q(1, 5), Q(1, 10)
    rpp, rsd = g - 2 * cp, g - cd
    invade_sd = rsd * Q(1, 2) / rpp
    invade_pp = rpp / rsd
    require(invade_sd == Q(3, 4) and invade_pp == Q(2, 3),
            "exact prior comparison")
    require(0 < invade_sd < 1 and 0 < invade_pp < 1,
            "strict pairwise bistability")
    return [str(invade_sd), str(invade_pp)]


def main():
    archive_check()
    cert = json.loads((HERE / "certificate.json").read_text())
    examples = json.loads((HERE / "examples.json").read_text())
    require(cert.get("schema") == examples.get("schema") == 1, "data schema")
    require(cert.get("global_theorem_mechanically_certified") is False,
            "theorem-certification scope")
    supports, all_supports = check_supports()
    identities = check_symbolic_identities()
    fixtures = check_floating_fixtures(examples["floating_positive_rate_fixtures"],
                                       all_supports)
    exact_examples, numerical = check_ordered_cost_examples(examples["ordered_cost_examples"])
    require(cert["support_counts"] == supports, "support certificate")
    require(cert["symbolic_identity_count"] == identities, "identity certificate")
    require(cert["floating_native_fixtures"] == fixtures, "fixture certificate")
    require(cert["exact_ordered_cost_equilibria"] == exact_examples,
            "ordered-cost certificate")
    require(cert["prior_invasion_multipliers"] == prior_comparison(),
            "prior comparison certificate")
    print(json.dumps({
        "status": "PASS",
        "support_counts": supports,
        "exact_symbolic_identities": identities,
        "exact_ordered_cost_equilibria": exact_examples,
        "floating_native_fixtures": fixtures,
        "ordered_cost_numerical_illustrations": numerical,
        "global_theorem_mechanically_certified": False,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
