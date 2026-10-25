"""Fresh exact postprocessing of widened serialized Linux interval enclosures.

The CAPD integration is an input, not independently reimplemented here.
Both the logged DP and an independently section-corrected flow derivative are
used to solve the branch sensitivity problem with exact Fraction arithmetic.
"""
from pathlib import Path
from fractions import Fraction as Q
import ast
import hashlib
import itertools
import json
import re
import argparse
import sympy as sp

CANDIDATE = Path(__file__).resolve().parents[1]
ROOT = CANDIDATE.parents[3]
LOG = CANDIDATE / "evidence/linux-sensitivity-output.txt"


class I:
    def __init__(self, a, b=None):
        self.a, self.b = Q(a), Q(a if b is None else b)
        assert self.a <= self.b

    def __add__(self, y):
        y = y if isinstance(y, I) else I(y)
        return I(self.a + y.a, self.b + y.b)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.b, -self.a)

    def __sub__(self, y):
        return self + -(y if isinstance(y, I) else I(y))

    def __rsub__(self, y):
        return I(y) - self

    def __mul__(self, y):
        y = y if isinstance(y, I) else I(y)
        candidates = [x*z for x in (self.a, self.b) for z in (y.a, y.b)]
        return I(min(candidates), max(candidates))

    __rmul__ = __mul__

    def __truediv__(self, y):
        y = y if isinstance(y, I) else I(y)
        assert y.a*y.b > 0
        return self * I(1/y.b, 1/y.a)

    def absmax(self):
        return max(abs(self.a), abs(self.b))

    def intersection(self, y):
        return I(max(self.a, y.a), min(self.b, y.b))

    def report(self):
        return {"exact": [str(self.a), str(self.b)],
                "decimal_display_only": [float(self.a), float(self.b)]}


def intervals(text):
    out = []
    for a, b in re.findall(r"\[([^,\[\]]+),\s*([^,\[\]]+)\]", text):
        a, b = Q(a), Q(b)
        widening = Q(1, 10**14)*max(1, abs(a), abs(b))
        out.append(I(a-widening, b+widening))
    return out


def block(text, start, end):
    return intervals(text.split("\n"+start+" ", 1)[1].split("\n"+end+" ", 1)[0])


def matrix(entries, n):
    assert len(entries) == n*n
    return [entries[i*n:(i+1)*n] for i in range(n)]


def matmul(a, b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))), I(0))
             for j in range(len(b[0]))] for i in range(len(a))]


def matvec(a, b):
    return [sum((x*y for x, y in zip(row, b)), I(0)) for row in a]


def solve(dp, preconditioner):
    A = [[dp[i+1][j+1]-(i == j) for j in range(4)] for i in range(4)]
    RA = matmul(preconditioner, A)
    E = [[(i == j)-RA[i][j] for j in range(4)] for i in range(4)]
    eta = max(sum(x.absmax() for x in row) for row in E)
    assert eta < 1
    rhs = [-dp[i+1][5] for i in range(4)]
    rb = matvec(preconditioner, rhs)
    bound = max(x.absmax() for x in rb)/(1-eta)
    z = [I(-bound, bound) for unused in range(4)]
    for unused in range(12):
        ez = matvec(E, z)
        z = [old.intersection(r+x) for old, r, x in zip(z, rb, ez)]
    return z, eta, bound


def first_jet():
    u, v, w, a, b, eps = sp.symbols("u v w a b eps")
    q = [u, v, w, a, b]
    p, s, h0 = 8-a-b, 40-u-v-w-a-b, u+100*b
    f0 = sp.Matrix([100*b-u, u-v, v-100*w*p+a,
                    100*w*p-sp.Rational(19, 10)*a, 3*p*s-101*b])
    # The first row of Dh0 is (1,0,0,0,100).
    h1 = (u+v)*h0-f0[0]-100*f0[4]
    h = h0+eps*h1
    induced = sp.Matrix([h*(1-eps*(u+v))-2*u, u-v,
                         v-100*w*p+a, 100*w*p-sp.Rational(19, 10)*a,
                         3*p*(s-eps*h)-101*b])
    invariance = eps*(sp.Matrix([h]).jacobian(q)*induced)[0] + h-eps*(u+v)*h-h0
    assert all(sp.expand(invariance).coeff(eps, i) == 0 for i in (0, 1))
    f1 = sp.Matrix([u+10000*b-300*p*s, 0, 0, 0, -3*p*h0])
    assert all(sp.expand(induced[i]).coeff(eps, 1) == sp.expand(f1[i]) for i in range(5))
    source = (CANDIDATE/"verify_period_sensitivity.cpp").read_text()
    expression = re.search(r'IMap vf\("par:r;var:u,v,w,a,b,e;fun:([^;]+);"\)', source).group(1)
    env = dict(zip(["u", "v", "w", "a", "b"], q))
    env.update(e=eps, r=sp.Rational(19, 10))
    encoded = [sp.sympify(part, locals=env) for part in expression.split(",")]
    assert all(sp.expand(encoded[i]-f0[i]-eps*f1[i]) == 0 for i in range(5))
    assert encoded[5] == 0
    return {"invariance_orders_0_1": True,
            "induced_first_jet_equals_cpp": True,
            "h1_identity": "(u+v)*(u+100*b)+u+10000*b-300*p*s",
            "F1": [str(x) for x in f1]}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'artifacts/validated-cycle/independent-period-check.json')
    args=parser.parse_args()
    if not __debug__: parser.error('Assertions must be enabled')
    text = LOG.read_text()
    seed = (CANDIDATE/"seed.h").read_text()
    Rrows = ast.literal_eval(re.search(r'inverse_values\[4\]\[4\] = (\{.*?\});', seed, re.S)
                            .group(1).replace("{", "[").replace("}", "]"))
    R = [[I(Q(float(x))) for x in row] for row in Rrows]
    determinant = Q(0)
    for perm in itertools.permutations(range(4)):
        term = Q((-1)**sum(perm[i] > perm[j] for i in range(4) for j in range(i+1, 4)))
        for i in range(4):
            term *= R[i][perm[i]].a
        determinant += term
    assert determinant != 0
    M = matrix(block(text, "FLOW_DERIVATIVE", "RETURN_MAP_DERIVATIVE"), 6)
    DP = matrix(block(text, "RETURN_MAP_DERIVATIVE", "KRAWCZYK"), 6)
    ret = block(text, "RETURN_BOX", "FLOW_DERIVATIVE")
    normal = 100*ret[4]-ret[0]
    assert normal.a > 0
    tau = [-M[0][j]/normal for j in range(6)]
    # Compute the full return-map derivative afresh from the flow enclosure.
    u, v, w, a, b = ret[:5]
    p, s = 8-a-b, 40-u-v-w-a-b
    f0 = [normal, u-v, v-100*w*p+a, 100*w*p-Q(19, 10)*a,
          3*p*s-101*b, I(0)]
    DP2 = [[M[i][j]+f0[i]*tau[j] for j in range(6)] for i in range(6)]
    outcomes = {}
    for name, derivative in [("logged_return_map_derivative", DP),
                             ("fresh_section_correction_from_flow", DP2)]:
        z, eta, bound = solve(derivative, R)
        T1 = tau[5]+sum((tau[j+1]*z[j] for j in range(4)), I(0))
        assert T1.a > 0
        assert Q("7.65771717") < T1.a < T1.b < Q("7.68396573")
        outcomes[name] = {
            "preconditioned_remainder_norm_upper_exact": str(eta),
            "preconditioned_remainder_norm_upper_display": float(eta),
            "initial_neumann_bound_exact": str(bound),
            "initial_neumann_bound_display": float(bound),
            "z_e": [x.report() for x in z],
            "T1": T1.report(), "T1_strictly_positive": True,
        }
    c = [Q(float(x)) for x in re.search(r'center_values\[4\] = \{([^}]+)\}', seed).group(1).split(",")]
    radius = Q(float("1e-9"))
    K = block(text, "KRAWCZYK", "FIXED_POINT_INCLUSION")
    assert len(K) == 4 and all(c[i]-radius < x.a < x.b < c[i]+radius for i, x in enumerate(K))
    metadata_path = LOG.with_name("run-metadata.json")
    metadata = json.loads(metadata_path.read_text())
    assert metadata["status"] == "completed" and metadata["conclusion"] == "success"
    for name in ("verify_period_sensitivity.cpp", "seed.h", "CMakeLists.txt", "build_and_run.sh"):
        recorded = json.loads((CANDIDATE/'evidence/linux-run.json').read_text())['source_sha256'][name]
        assert hashlib.sha256((CANDIDATE/name).read_bytes()).hexdigest() == recorded
    files = [LOG, metadata_path, CANDIDATE/"verify_period_sensitivity.cpp", CANDIDATE/"seed.h",
             CANDIDATE/"proof_link.md", CANDIDATE/"verify_first_jet.py"]
    report = {
        "status": "passed",
        "scope": "Fresh exact rational postprocessing of widened Linux interval enclosures; CAPD integration itself is an input. Krawczyk inclusion is independently compared with X from logged K, not independently reconstructed because this log omits the point residual.",
        "first_jet": first_jet(),
        "logged_intervals_widened_on_each_side_by": "10^-14 * max(1, abs(lower), abs(upper))",
        "binary_preconditioner_determinant_exact": str(determinant),
        "linux_ci_provenance": {"metadata": metadata,
                                "numerical_inputs_match_archived_run_hashes": True},
        "logged_krawczyk_strictly_inside_exact_center_radius_box": True,
        "normal_velocity": normal.report(),
        "tau_z": [x.report() for x in tau[1:5]], "tau_e": tau[5].report(),
        "sensitivity_solves": outcomes,
        "simple_safe_exact_decimal_enclosure": ["7.65771717", "7.68396573"],
        "input_sha256": {str(x.relative_to(ROOT)): hashlib.sha256(x.read_bytes()).hexdigest() for x in files},
        "implementation_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({"status": "passed", "first_jet": report["first_jet"],
                      "coefficient_enclosures": {name: values["T1"]["decimal_display_only"]
                                                 for name, values in outcomes.items()}}, indent=2))


if __name__ == "__main__":
    main()
