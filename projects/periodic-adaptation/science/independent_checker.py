"""Portable exact algebra checker for the accompanying existing reaction model.

Requires Python 3 and SymPy. Reads model.json and proof.md in this directory.
Checks symbolic identities and records file hashes; DOES NOT verify the ODE proof,
novelty, literature, peer review or a formal release/provenance gate.
"""
from pathlib import Path
import hashlib
import json
import sympy as sp

HERE = Path(__file__).resolve().parent
EXPECTED_SHA256 = {
    "model.json": "bf8a492e3824eea0e3a0f0ef493810a079a18d541d9984bf832b02245f8556c5",
    "proof.md": "3084f8f5f0a0fc1c41fc7b1fb560ba5711e17e99f38b2b4f3bdef9851a96a05c",
}
verified_sha256 = {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                   for name in EXPECTED_SHA256}
assert verified_sha256 == EXPECTED_SHA256, "Companion file differs from reviewed bytes"

a, b, c, d, u, x, y = sp.symbols("a b c d u x y", positive=True)
checks = {}
def ident(label, expression):
    reduced = sp.cancel(expression)
    assert reduced == 0, (label, reduced)
    checks[label] = "zero rational function"

# Species rows are source labels (1=S, 2=Y, 3=X). Reaction columns are
# S->X, X->S, S->X+Y, X+Y->S, derived as product minus reactant vectors.
reactants = sp.Matrix([[1, 0, 1, 0], [0, 0, 0, 1], [0, 1, 0, 1]])
products = sp.Matrix([[0, 1, 0, 1], [0, 0, 1, 0], [1, 0, 1, 0]])
stoich = products - reactants
fluxes = sp.Matrix([a*u, b*x, c*u, d*x*y])
rhs = stoich * fluxes
fx, fy = rhs[2], rhs[1]
assert stoich.T.nullspace() == [sp.Matrix([1, 0, 1])]
ident("reaction-derived x", fx - ((a+c)*u-b*x-d*x*y))
ident("reaction-derived y", fy - (c*u-d*x*y))

h = (a+c)/c
yb = b*c/(a*d)
alpha = 1/h
q = a*d/c
ident("input-eliminated w dynamics", fy-alpha*fx - a*d/(a+c)*x*(yb-y))
eta = y-yb
z = x-h*eta
ident("z derivative", fx-h*fy-q*x*eta)
ident("log identity numerator", q*x*eta-q*z*eta-q*h*eta**2)
ident("constant equilibrium x", fx.subs({x:a*u/b, y:yb}))
ident("constant equilibrium y", fy.subs({x:a*u/b, y:yb}))

# Face identities use explicit nonnegative slack variables; their signs are
# then immediate for every allowed u<=U, y>=0, slack>=0, ell>=0.
U, slack, ell = sp.symbols("U slack ell", nonnegative=True)
M = (a+c)*U/b + slack
L = yb+ell
ident("upper-x boundary decomposition", fx.subs(x,M)-((a+c)*(u-U)-b*slack-d*M*y))
ident("slanted boundary decomposition", (fy-alpha*fx).subs(y,L+x/h) + a*d/(a+c)*x*(ell+x/h))

# Derive difference matrix directly by substitution of both independent states.
X, Y, p, v = sp.symbols("X Y p v", real=True)
f1 = sp.Matrix([fx, fy]).subs({x:X+p, y:Y-v}, simultaneous=True)
f2 = sp.Matrix([fx, fy]).subs({x:X, y:Y}, simultaneous=True)
df = f1-f2
A = sp.Matrix([[-b-d*(Y-v), d*X], [d*(Y-v), -d*X]])
ident("first difference equation", df[0]-(A*sp.Matrix([p,v]))[0])
ident("second difference equation", -df[1]-(A*sp.Matrix([p,v]))[1])
eps, boundY, m = sp.symbols("eps boundY m", positive=True)
ident("weighted first column", A[0,0]+(1+eps)*A[1,0]-(-b+eps*d*(Y-v)))
ident("weighted second column", A[0,1]+(1+eps)*A[1,1]+eps*d*X)
ident("first-column worst-case", (-b+eps*d*boundY).subs(eps,b/(2*d*boundY))+b/2)

# Independent harmonic calculation: solve the two complex amplitude equations.
u0, omega, s = sp.symbols("u0 omega s", positive=True)
xb = a*u0/b
J = sp.Matrix([fx,fy]).jacobian([x,y]).subs({x:xb,y:yb})
B = b+d*(xb+yb)
K = b*d*xb
ident("Hurwitz trace", sp.trace(J)+B)
ident("Hurwitz determinant", J.det()-K)
harmonic = (s*sp.eye(2)-J).inv()*sp.Matrix([a+c,c])
ident("y harmonic response", harmonic[1]-c*s/(s*s+B*s+K))
gain_squared = sp.simplify((harmonic[1].subs(s,sp.I*omega)*harmonic[1].subs(s,-sp.I*omega)))
denom = (K-omega**2)**2+B**2*omega**2
ident("sinusoidal squared gain", gain_squared-c*c*omega**2/denom)
vfreq = sp.symbols("vfreq", positive=True)
denv = (K-vfreq)**2+B**2*vfreq
ident("frequency derivative numerator", denv-vfreq*sp.diff(denv,vfreq)-(K*K-vfreq*vfreq))

result = {
    "job_id": None,
    "proposer_job_id": None,
    "source_binding_type": "fixed-reviewed-model-and-proof-byte-hashes-only",
    "verified_sha256": verified_sha256,
    "formal_freeze_passed": False,
    "symbolic_checks": checks,
    "symbolic_check_count": len(checks),
    "reaction_stoichiometry_source_species_1_2_3": [list(stoich.row(i)) for i in range(3)],
    "conservation_basis_source_species_1_2_3": [[1,0,1]],
    "scope": "Symbolic identities and exact file hashes only; analytic ODE steps require the written referee argument. No finite grid is used as evidence of a universal statement.",
}
print(json.dumps(result, sort_keys=True, indent=2, default=str))
