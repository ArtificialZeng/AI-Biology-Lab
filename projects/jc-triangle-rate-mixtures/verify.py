"""Read-only finite exact checks for the JC triangle scalar-rate note.

Python standard library only. The arbitrary-law theorem is proved in proof.md;
these five examples verify the probability model and implementation exactly.
"""
import sys
sys.dont_write_bytecode = True
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
INPUTS = {"README.md", "proof.md", "sources.md", "witnesses.json",
          "certificate.json", "verify.py", "LICENSE"}
LETTERS = "AGCT"
STATES = range(4)
PATTERNS = list(product(STATES, repeat=3))
NAMES = ["".join(LETTERS[s] for s in pat) for pat in PATTERNS]
ZERO = F(0)

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def check_archive():
    manifest = json.loads((HERE / "archive_manifest.json").read_text())
    require(manifest["algorithm"] == "sha256", "archive algorithm")
    require(set(manifest["files"]) == INPUTS, "archive input scope")
    actual = set()
    for p in HERE.rglob("*"):
        require(not p.is_symlink(), "archive symlink")
        if p.is_file():
            actual.add(p.relative_to(HERE).as_posix())
    require(actual == INPUTS | {"archive_manifest.json"}, "unexpected archive file")
    for name, expected in manifest["files"].items():
        require(hashlib.sha256((HERE/name).read_bytes()).hexdigest() == expected,
                "archive hash differs: " + name)


def matrix(a):
    return [[(1 + 3 * a) / 4 if i == j else (1 - a) / 4 for j in STATES] for i in STATES]

def displayed_tree_probabilities(a):
    """Keep all three original internal vertices, even degree-two vertices.

    T4 deletes a6: H-L1:a1,U-L2:a2,V-L3:a3,U-H:a4,U-V:a5.
    T6 deletes a4: H-L1:a1,U-L2:a2,V-L3:a3,U-V:a5,V-H:a6.
    Root at U with probability 1/4. No contraction is used here.
    """
    m1, m2, m3, m4, m5, m6 = [matrix(v) for v in a]
    t4, t6 = ([], [])
    for s1, s2, s3 in PATTERNS:
        p4, p6 = (ZERO, ZERO)
        for h, u, v in product(STATES, repeat=3):
            pendant = F(1, 4) * m1[h][s1] * m2[u][s2] * m3[v][s3]
            p4 += pendant * m4[u][h] * m5[u][v]
            p6 += pendant * m5[u][v] * m6[v][h]
        t4.append(p4)
        t6.append(p6)
    return (t4, t6)

def star_probabilities(b):
    """Separate implementation used only to verify the contraction explanation."""
    m = [matrix(v) for v in b]
    return [sum((F(1, 4) * m[0][c][s[0]] * m[1][c][s[1]] * m[2][c][s[2]] for c in STATES), ZERO) for s in PATTERNS]

def character(k, s):
    return -1 if sum((bin(ki & si).count('1') for ki, si in zip(k, s))) % 2 else 1

def fourier(p):
    return [sum((c * v for c, v in zip(row, p)), ZERO) for row in SIGNS]

def inverse_fourier(q):
    return [sum((SIGNS[k][i] * q[k] for k in range(64)), ZERO) / 64 for i in range(64)]

def eq6(a, delta):
    a1, a2, a3, a4, a5, a6 = a
    return (a2 * a3 * a5, delta * a1 * a3 * a4 * a5 + (1 - delta) * a1 * a3 * a6, delta * a1 * a2 * a4 + (1 - delta) * a1 * a2 * a5 * a6, delta * a1 * a2 * a3 * a4 * a5 + (1 - delta) * a1 * a2 * a3 * a5 * a6)

def lift_four_coordinates(q):
    """Lift source's four classes, qAAA=1, and all 48 zero entries."""
    x, y, z, t = q
    answer = []
    for a, b, c in PATTERNS:
        if a ^ b ^ c:
            answer.append(ZERO)
        elif a == b == c == 0:
            answer.append(F(1))
        elif a == 0:
            answer.append(x)
        elif b == 0:
            answer.append(y)
        elif c == 0:
            answer.append(z)
        else:
            answer.append(t)
    return answer

def four_from_transform(q):
    return tuple((q[NAMES.index(s)] for s in ('ACC', 'CAC', 'CCA', 'CGT')))

def predicate_values(q, orientation):
    """Full source Theorem3.1 iff set, with direct leaf permutations."""
    x, y, z, t = q
    if orientation == 2:
        x, y = (y, x)
    elif orientation == 3:
        x, z = (z, x)
    values = {}
    for name, v in zip(('x', 'y', 'z', 't'), (x, y, z, t)):
        values[name] = v
        values['1-' + name] = 1 - v
    values.update({'x-t': x - t, 't-x*y': t - x * y, 'x*y*z-t^2': x * y * z - t * t, 't+x*(t-y-z)': t + x * (t - y - z)})
    return values

def memberships(q):
    ans = {}
    for k in (1, 2, 3):
        vals = predicate_values(q, k)
        ans['M' + str(k)] = {'strict_member': all((v > 0 for v in vals.values())), 'margins': vals, 'failed_conditions': [s for s, v in vals.items() if v <= 0]}
    return ans

def covariance(f, g, weights):
    ef = sum((w * v for w, v in zip(weights, f)), ZERO)
    eg = sum((w * v for w, v in zip(weights, g)), ZERO)
    direct = sum((w * u * v for w, u, v in zip(weights, f, g)), ZERO) - ef * eg
    pair = sum((wi * wj * (fi - fj) * (gi - gj) / 2 for wi, fi, gi in zip(weights, f, g) for wj, fj, gj in zip(weights, f, g)), ZERO)
    require(direct == pair, 'exact check failed: direct == pair')
    return direct

def preservation(qs, weights):
    xs, ys, zs, ts = map(list, zip(*qs))
    hs = [y + z - t for y, z, t in zip(ys, zs, ts)]
    mixed = [sum((w * q[i] for w, q in zip(weights, qs)), ZERO) for i in range(4)]
    x, y, z, t = mixed
    reports = {}
    for label, g, point, mixed_margin in (('T-X*Y', ys, [t - x * y for x, y, z, t in qs], t - x * y), ('T-X*Z', zs, [t - x * z for x, y, z, t in qs], t - x * z), ('Phi_X', hs, [t - x * (y + z - t) for x, y, z, t in qs], t - x * (y + z - t))):
        expectation = sum((w * v for w, v in zip(weights, point)), ZERO)
        cov = covariance(xs, g, weights)
        require(expectation + cov == mixed_margin, 'exact check failed: expectation + cov == mixed_margin')
        require(all((v > 0 for v in point)) and cov >= 0 and (mixed_margin > 0), 'exact check failed: all((v > 0 for v in point)) and cov >= 0 and (mixed_margin > 0)')
        reports[label] = {'pointwise_margins': point, 'expected_margin': expectation, 'covariance': cov, 'mixed_margin': mixed_margin, 'identity_residual': mixed_margin - expectation - cov}
    monotonic = {name: all((v > u for v, u in zip(vals, vals[1:]))) for name, vals in (('x', xs), ('y', ys), ('z', zs), ('h=y+z-t', hs))}
    require(all(monotonic.values()) and all((v > 0 for v in hs)), 'exact check failed: all(monotonic.values()) and all((v > 0 for v in hs))')
    reports['ordered_rates_strictly_decreasing_functions'] = monotonic
    reports['positive_h_values'] = hs
    return reports

def jsonable(obj):
    if isinstance(obj, F):
        return str(obj)
    if isinstance(obj, dict):
        return {str(k): jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (tuple, list)):
        return [jsonable(v) for v in obj]
    return obj

SIGNS = [[character(k, s) for s in PATTERNS] for k in PATTERNS]

def compute():
    raw = json.loads((HERE / "witnesses.json").read_text())
    require(set(raw) == {"cases"} and len(raw["cases"]) == 5, "witness scope")
    results = {}
    total_categories = 0
    for case in raw["cases"]:
        name = case["name"]
        base = list(map(F, case["a"]))
        delta = F(case["delta"])
        rates = case["rates"]
        weights = list(map(F, case["weights"]))
        require(len(base) == 6 and all(0 < v < 1 for v in base), "strict edges")
        require(0 < delta < 1, "strict generating inheritance")
        require(all(type(r) is int and r > 0 for r in rates)
                and rates == sorted(set(rates)), "positive integer example rates")
        require(len(rates) == len(weights) and all(w > 0 for w in weights)
                and sum(weights) == 1, "mixture weights")
        ps, qs, category_reports = [], [], []
        for rate in rates:
            edges = [v ** rate for v in base]
            p4, p6 = displayed_tree_probabilities(edges)
            b4 = [edges[0]*edges[3], edges[1], edges[2]*edges[4]]
            b6 = [edges[0]*edges[5], edges[1]*edges[4], edges[2]]
            require(p4 == star_probabilities(b4) and p6 == star_probabilities(b6),
                    "original internal-state sums and path contractions differ")
            require(sum(p4) == sum(p6) == 1 and min(p4) > 0 and min(p6) > 0,
                    "displayed tree is not strictly stochastic")
            prob = [delta*u + (1-delta)*v for u, v in zip(p4, p6)]
            transformed = fourier(prob)
            require(transformed == lift_four_coordinates(eq6(edges, delta)),
                    "full Fourier vector and source Eq6 differ")
            require(inverse_fourier(transformed) == prob, "inverse transform")
            coords = four_from_transform(transformed)
            require(memberships(coords)["M1"]["strict_member"], "source category membership")
            ps.append(prob)
            qs.append(coords)
            category_reports.append({"rate": rate, "coordinates": coords,
                                     "all64_native_and_Fourier_equal": True})
            total_categories += 1
        prob = [sum((w*p[i] for w, p in zip(weights, ps)), ZERO) for i in range(64)]
        transformed = fourier(prob)
        q = tuple(sum((w*v[i] for w, v in zip(weights, qs)), ZERO) for i in range(4))
        require(transformed == lift_four_coordinates(q), "mixed all64 Fourier coordinates")
        require(inverse_fourier(lift_four_coordinates(q)) == prob, "mixed all64 probabilities")
        require(sum(prob) == 1 and min(prob) > 0, "mixed strictly stochastic distribution")
        require(sum(v == 0 for v in transformed) == 48 and transformed[0] == 1,
                "JC Fourier zero/one classes")
        x, y, z, t = q
        determinant = x*y*z - t*t
        member = memberships(q)
        require(member["M1"]["strict_member"] == (determinant > 0), "example dichotomy")
        require(determinant > 0 or not any(v["strict_member"] for v in member.values()),
                "example wrong-orientation compatibility")
        orbit = {key: prob[NAMES.index(key)] for key in ("AAA", "ACC", "CAC", "CCA", "CGT")}
        counts = {key: 0 for key in orbit}
        for pat, value in zip(PATTERNS, prob):
            aa, bb, cc = pat
            key = ("AAA" if aa == bb == cc else "ACC" if bb == cc else
                   "CAC" if aa == cc else "CCA" if aa == bb else "CGT")
            require(value == orbit[key], "JC probability equality orbit")
            counts[key] += 1
        require(list(counts.values()) == [4, 12, 12, 12, 24], "orbit counts")
        results[name] = {"coordinates": q, "D": determinant,
                         "memberships": member, "probability_sum": sum(prob),
                         "minimum_probability": min(prob), "orbit_probabilities": orbit,
                         "orbit_counts": counts, "full_probabilities": prob,
                         "full_Fourier": transformed, "categories": category_reports,
                         "covariance_checks": preservation(qs, weights)}
    expected_q = (F(9,128), F(53,1024), F(53,1024), F(33,2048))
    example = results["equal_edges_two_rate_nonclosure"]
    require(example["coordinates"] == expected_q and example["D"] == F(-9567,134217728),
            "nonclosure witness exact values")
    require(results["asymmetric_three_rate_true_only"]["memberships"]["M1"]["strict_member"]
            and not results["asymmetric_three_rate_true_only"]["memberships"]["M2"]["strict_member"]
            and not results["asymmetric_three_rate_true_only"]["memberships"]["M3"]["strict_member"],
            "heldout true-only example")
    return jsonable({"schema": 1, "finite_examples": len(results),
                     "native_rate_categories": total_categories,
                     "universal_proof": "proof.md; finite checks are not that proof",
                     "cases": results})

def main():
    require(len(sys.argv) == 1, "usage: python verify.py")
    check_archive()
    actual = compute()
    expected = json.loads((HERE / "certificate.json").read_text())
    require(actual == expected, "exact static certificate differs")
    print(json.dumps({"status": "PASS", "finite_examples": actual["finite_examples"],
                      "native_rate_categories": actual["native_rate_categories"],
                      "all64_probability_and_Fourier_checks": True,
                      "archive_modified": False}, sort_keys=True))

if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, ValueError, KeyError, TypeError, OSError) as exc:
        print("FAIL: " + str(exc), file=sys.stderr)
        sys.exit(1)
