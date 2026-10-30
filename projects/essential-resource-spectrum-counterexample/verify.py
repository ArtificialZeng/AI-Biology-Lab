"""Read-only exact reconstruction of the published two-pathway resource model."""
import sys
sys.dont_write_bytecode = True
import hashlib
import json
from pathlib import Path
import sympy as s

HERE = Path(__file__).resolve().parent
INPUTS = {"README.md", "proof.md", "sources.md", "witness.json",
          "certificate.json", "verify.py", "requirements.txt", "LICENSE"}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_archive():
    manifest = json.loads((HERE / "archive_manifest.json").read_text())
    require(manifest["algorithm"] == "sha256", "unsupported archive algorithm")
    require(set(manifest["files"]) == INPUTS, "archive input scope differs")
    actual = set()
    for path in HERE.rglob("*"):
        require(not path.is_symlink(), "archive contains a symlink")
        if path.is_file():
            actual.add(path.relative_to(HERE).as_posix())
    require(actual == INPUTS | {"archive_manifest.json"}, "unexpected archive file")
    for name in sorted(INPUTS):
        require(digest(HERE / name) == manifest["files"][name],
                "archive hash mismatch: " + name)


def matrix_strings(value):
    return [[str(x) for x in row] for row in value.tolist()]


def positive_polynomial(value, n, label):
    coefficients = s.Poly(s.expand(value), n).all_coeffs()
    require(all(c >= 0 for c in coefficients) and any(c > 0 for c in coefficients),
            label + " is not coefficientwise positive for n>0")


def routh_column(poly, x):
    """Independent ordinary Routh recursion; used away from the boundary."""
    coefficients = s.Poly(poly, x).all_coeffs()
    width = (len(coefficients) + 1) // 2
    row0 = coefficients[0::2] + [s.Integer(0)] * (width-len(coefficients[0::2]))
    row1 = coefficients[1::2] + [s.Integer(0)] * (width-len(coefficients[1::2]))
    rows = [row0, row1]
    while len(rows) < len(coefficients):
        a, b = rows[-2], rows[-1]
        require(b[0] != 0, "zero Routh pivot at test density")
        rows.append([s.cancel((b[0]*a[j+1]-a[0]*b[j+1])/b[0])
                     for j in range(width-1)] + [s.Integer(0)])
    return [int(s.sign(row[0])) for row in rows]


def reconstruct(witness):
    require(set(witness) == {"parameters"}, "unexpected witness schema")
    p = witness["parameters"]
    require(set(p) == {"C", "P", "Pt", "e", "theta", "r"}, "parameter schema")
    def mat(name):
        raw = p[name]
        require(len(raw) == 3 and all(len(row) == 3 for row in raw), "matrix shape")
        return s.Matrix([[s.Rational(x) for x in row] for row in raw])
    C, P, Pt = (mat(name) for name in ("C", "P", "Pt"))
    e, theta, r = (s.Rational(p[name]) for name in ("e", "theta", "r"))
    one = s.ones(3, 1)
    cd = C[0, 0]
    require(0 < e < 1 and 0 < theta < 1 and r > 0, "scalar physical domain")
    require(all(v > 0 for v in C) and all(C[i,i] == cd for i in range(3)),
            "positive common-diagonal consumption")
    require(all(v >= 0 for v in P) and all(v >= 0 for v in Pt), "recycling positivity")
    require(P*one == Pt*one == one, "recycling must be ROW stochastic")
    require(all(P[i,i] == Pt[i,i] == 0 for i in range(3)), "witness zero diagonals")
    eps = s.Matrix(3, 3, lambda i,j: e if i == j else 1)
    margins = [r*(C[j,k]-e*cd) for j in range(3) for k in range(3) if j != k]
    require(all(v > 0 for v in margins), "strict Liebig limiter fails")
    actual_min = [min(eps[j,k]*C[j,k]*r for j in range(3)) for k in range(3)]
    require(actual_min == [e*cd*r]*3, "literal original minimum differs")
    n = s.Symbol("n", positive=True)
    x = s.Symbol("x")
    R = s.Matrix(s.symbols("R0:3", positive=True))
    N = s.Matrix(s.symbols("N0:3", positive=True))
    rho = s.Matrix(s.symbols("rho0:3"))
    beta = (1-theta)*e*cd
    eta = beta*r
    # This branch is justified by the literal minimum and strict margins above.
    g = s.Matrix([eps[k,k]*C[k,k]*R[k] for k in range(3)])
    # Original Eq. (1) component sums: excess AND biomass conversion pathways.
    resources = s.Matrix([rho[i]-R[i]*sum(C[i,k]*N[k] for k in range(3))
        +sum(P[j,i]*sum((C[j,k]*R[j]-g[k]/eps[j,k])*N[k]
                        for k in range(3)) for j in range(3))
        +theta*sum(Pt[j,i]*sum(g[k]/eps[j,k]*N[k] for k in range(3))
                   for j in range(3)) for i in range(3)])
    consumers = s.Matrix([N[k]*((1-theta)*g[k]-eta) for k in range(3)])
    field = resources.col_join(consumers)
    eq = {**dict(zip(R,[r]*3)), **dict(zip(N,[n]*3))}
    nofeed = resources.subs(dict(zip(rho,[0]*3))).subs(eq)
    feed = s.simplify(-nofeed/n)
    require(all(v > 0 for v in feed), "not every inflow is positive")
    balance = {**eq, **dict(zip(rho,n*feed))}
    require(all(s.expand(v.subs(balance)) == 0 for v in field), "Eq1 equilibrium")
    # A second conservation oracle, before equilibrium substitution.
    weight = 1/e+2
    mass = sum(resources)+weight*sum(consumers)-sum(rho)+weight*eta*sum(N)
    require(s.expand(mass) == 0, "original biomass accounting fails")
    J = field.jacobian(list(R)+list(N)).subs(balance)
    A, B = s.simplify(J[:3,:3]/n), s.simplify(J[:3,3:]/r)
    Et = s.Matrix(3,3,lambda i,j: 1 if i == j else e)
    D = s.diag(*(C*one))
    Bformula = -C+P.T*(C-cd*Et)+theta*cd*Pt.T*Et
    Aformula = -D+P.T*D-cd*P.T*Et+theta*cd*Pt.T*Et
    require(A == Aformula and B == Bformula, "component/block derivatives disagree")
    block = (n*A).row_join(r*B).col_join((beta*n*s.eye(3)).row_join(s.zeros(3)))
    require(J == block, "full original Jacobian block mismatch")
    require(feed == -r*B*one, "feed and B normalization mismatch")
    kappa = (1-theta)*cd*(1+2*e)
    require(one.T*A == one.T*B == -kappa*one.T, "aggregate conservation block")
    require(A-B == -(s.eye(3)-P.T)*(D-C), "A/B decomposition")
    require(s.expand(sum(n*feed)-3*r*n*kappa) == 0, "total feed balance")
    # Full determinant, then a separate smaller quadratic matrix-pencil route.
    full = s.Poly(J.charpoly(x).as_expr(),x)
    pencil = (x*x*s.eye(3)-n*x*A-beta*n*r*B).det()
    require(s.expand(full.as_expr()-pencil) == 0, "independent pencil route differs")
    aggregate = s.Poly(x*x+n*kappa*x+beta*n*r*kappa,x)
    quartic, remainder = s.div(full,aggregate)
    require(remainder.is_zero and quartic.degree() == 4, "aggregate factor division")
    bshape, remainder = s.div(s.Poly(B.charpoly(x).as_expr(),x),s.Poly(x+kappa,x))
    require(remainder.is_zero and bshape.degree() == 2, "B aggregate factor")
    disc = s.discriminant(bshape.as_expr(),x)
    require(disc > 0 and all(v > 0 for v in bshape.all_coeffs()), "B is not real negative")
    coefficients = [s.factor(v) for v in quartic.all_coeffs()]
    require(coefficients[0] == 1, "quartic monicity")
    a1,a2,a3,a4 = coefficients[1:]
    d2 = s.factor(a1*a2-a3)
    d3 = s.factor(a1*a2*a3-a3*a3-a1*a1*a4)
    for label,value in zip(("a1","a2","a3","a4","delta2"),(a1,a2,a3,a4,d2)):
        positive_polynomial(value,n,label)
    threshold_poly = s.Poly(s.cancel(d3/n**4),n)
    require(threshold_poly.degree() == 1, "threshold is not linear")
    slope,intercept = threshold_poly.all_coeffs()
    require(slope > 0 and intercept < 0, "threshold sign classification fails")
    nc = s.factor(-intercept/slope)
    require(s.expand(d3-n**4*slope*(n-nc)) == 0, "all-n threshold identity")
    low,high = s.Rational(1,10000000),s.Integer(1)
    require(0 < low < nc < high, "positive low/high density placement")
    column = [s.Integer(1),a1,s.cancel(d2/a1),s.cancel(d3/d2),a4]
    low_signs = [int(s.sign(v.subs(n,low))) for v in column]
    high_signs = [int(s.sign(v.subs(n,high))) for v in column]
    require(low_signs == [1,1,1,-1,1] and high_signs == [1,1,1,1,1], "quartic signs")
    full_low = routh_column(full.as_expr().subs(n,low),x)
    full_high = routh_column(full.as_expr().subs(n,high),x)
    require(full_low == [1,1,1,1,1,-1,1] and full_high == [1]*7, "full Routh recursion")
    # Scoped S=2 algebra, including repeated B eigenvalues and self recycling.
    u,v,c12,c21,cd2 = s.symbols("u v c12 c21 cd2", real=True)
    P2 = s.Matrix([[1-u,u],[v,1-v]])
    C2 = s.Matrix([[cd2,c12],[c21,cd2]])
    L2 = (s.eye(2)-P2.T)*(s.diag(*(C2*s.ones(2,1)))-C2)
    w = s.Matrix([1,-1])
    ell = 2*(u*c12+v*c21)
    require(s.simplify(L2*w-ell*w) == s.zeros(2,1), "S2 identity")
    return {"parameters":p, "Cdiag":str(cd), "strict_limiter_min_margin":str(min(margins)),
        "eta":str(eta), "rho_per_n":[str(v) for v in feed],
        "biomass_weight":str(weight), "A":matrix_strings(A), "B":matrix_strings(B),
        "aggregate_kappa":str(kappa), "aggregate_polynomial":str(aggregate.as_expr()),
        "B_composition_poly":str(bshape.as_expr()), "B_composition_discriminant":str(disc),
        "J_composition_coefficients":[str(v) for v in coefficients],
        "routh_delta2":str(d2), "routh_delta3":str(d3), "n_critical":str(nc),
        "n_low":str(low), "low_routh_signs":low_signs,
        "n_high":str(high), "high_routh_signs":high_signs,
        "full_low_routh_signs":full_low, "full_high_routh_signs":full_high,
        "full_determinant_and_pencil_agree":True,
        "S2_identity":"ell=2*(u*c12+v*c21)"}


def main():
    require(s.__version__ == "1.13.3", "install the pinned sympy==1.13.3")
    require(len(sys.argv) == 1, "usage: python verify.py")
    check_archive()
    witness = json.loads((HERE/"witness.json").read_text())
    computed = reconstruct(witness)
    archived = json.loads((HERE/"certificate.json").read_text())
    require(computed == archived, "static certificate differs from exact reconstruction")
    print(json.dumps({"status":"PASS", "checks":"archive, original Eq1, exact spectrum and all-n Routh",
                      "n_critical":computed["n_critical"], "archive_modified":False}, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, ValueError, KeyError, TypeError, OSError) as error:
        print("FAIL: " + str(error), file=sys.stderr)
        sys.exit(1)
