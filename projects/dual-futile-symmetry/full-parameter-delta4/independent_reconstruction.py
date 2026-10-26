"""Independent exact reconstruction from source reaction columns and reduction.

No proposer verifier, multiplier, cache, or numerical optimization is imported.
Run under Python -O: essential checks use explicit exceptions, not assert.
"""
from pathlib import Path
from hashlib import sha256
import importlib.util
import itertools
import json
import platform
import sys

import sympy as s
import flint
from flint import fmpz_mpoly_ctx

sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
INPUT = HERE
NAMES = ('a','b4','b10','c','e1','e4','f7','f10','u2','u3','v5','v6','x8','x9','y11','y12')

def require(truth, message):
    if not truth:
        raise RuntimeError(message)

from integrity import verify_manifest
verify_manifest(INPUT)

# Species order A,B,C,E,F,U,V,X,Y. Source (1.1)/(1.2) supplies the
# twelve elementary stoichiometric columns. Only r1,r4,r7,r10 have two
# reactant derivatives; source pp.1-2 labels all sixteen derivatives.
S = s.zeros(9,12)
reactions = [
    {0:-1,3:-1,5:1}, {0:1,3:1,5:-1}, {1:1,3:1,5:-1},
    {1:-1,3:-1,6:1}, {1:1,3:1,6:-1}, {2:1,3:1,6:-1},
    {2:-1,4:-1,7:1}, {2:1,4:1,7:-1}, {1:1,4:1,7:-1},
    {1:-1,4:-1,8:1}, {1:1,4:1,8:-1}, {0:1,4:1,8:-1}]
for j, col in enumerate(reactions):
    for i, value in col.items(): S[i,j] = value
symbols = s.symbols(' '.join(NAMES))
a,b4,b10,c,e1,e4,f7,f10,u2,u3,v5,v6,x8,x9,y11,y12 = symbols
Jrate = s.zeros(12,9)
for i, row in enumerate([{0:a,3:e1},{5:u2},{5:u3},{1:b4,3:e4},
        {6:v5},{6:v6},{2:c,4:f7},{7:x8},{7:x9},{1:b10,4:f10},{8:y11},{8:y12}]):
    for j, value in row.items(): Jrate[i,j] = value
J9 = S * Jrate

# q=(Q1,Q2,U,V,Y,X), Q1=A+U, Q2=C+X. On fixed totals,
# A=Q1-U, C=Q2-X, B=W3-Q1-Q2-V-Y, E=W1-U-V, F=W2-X-Y.
T = s.Matrix([[1,0,-1,0,0,0],[-1,-1,0,-1,-1,0],[0,1,0,0,0,-1],
              [0,0,-1,-1,0,0],[0,0,0,0,-1,-1],[0,0,1,0,0,0],
              [0,0,0,1,0,0],[0,0,0,0,0,1],[0,0,0,0,1,0]])
P = s.Matrix([[1,0,0,0,0,1,0,0,0],[0,0,1,0,0,0,0,1,0],
              [0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,1,0,0],
              [0,0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,1,0]])
C = s.Matrix([[0,0,0,1,0,1,1,0,0],[0,0,0,0,1,0,0,1,1],
              [1,1,1,0,0,1,1,1,1]])
require(P*T == s.eye(6), 'Projection/embedding mismatch')
require(C*S == s.zeros(3,12) and C*T == s.zeros(3,6), 'Conservation failure')
gred = P * J9 * T
require(J9*T == T*gred, 'Tangent subspace not invariant')

# Compare proposer matrix only after deriving ours from the primary model.
spec = importlib.util.spec_from_file_location('package_matrix', INPUT / 'matrix.py')
owner_matrix = importlib.util.module_from_spec(spec)
spec.loader.exec_module(owner_matrix)
require(tuple(owner_matrix.NAMES) == NAMES, 'Name order mismatch')
require(gred == owner_matrix.matrix(symbols, s.Matrix), 'Package matrix differs from source reduction')

ctx = fmpz_mpoly_ctx.get(NAMES)
zero, one = ctx.constant(0), ctx.constant(1)
def convert(expr):
    return ctx.from_dict({powers:int(value) for powers,value in s.Poly(expr,*symbols).terms()})
g = [[convert(gred[i,j]) for j in range(6)] for i in range(6)]
def matmul(A,B):
    return [[sum((A[i][k]*B[k][j] for k in range(6)),zero) for j in range(6)] for i in range(6)]

# Faddeev-LeVerrier uses traces and exact division rather than the owner's
# principal minors or discovery determinant. Each quotient is multiplied
# back and checked so an inexact integer quotient cannot slip through.
B = [[one if i==j else zero for j in range(6)] for i in range(6)]
cs = [one]
for k in range(1,7):
    AB = matmul(g,B)
    tr = sum((AB[i][i] for i in range(6)),zero)
    ck = (-tr)//k
    require(k*ck == -tr, 'Nonexact Newton division')
    cs.append(ck)
    B = [[AB[i][j]+(ck if i==j else zero) for j in range(6)] for i in range(6)]
require(all(v.is_zero() for row in B for v in row), 'Cayley-Hamilton residual')

# Derive the fourth Hurwitz polynomial using SymPy's generic determinant
# before evaluation; this does not trust a copied closed formula.
z = s.symbols('z0:7')
H4 = s.Matrix([[z[2*j-i+1] if 0 <= 2*j-i+1 <= 6 else 0 for j in range(4)] for i in range(4)])
generic = s.Poly(s.expand(H4.det().subs(z[0],1)),*z[1:])
d4 = zero
for exponents,value in generic.terms():
    term = ctx.constant(int(value))
    for ck,e in zip(cs[1:],exponents): term *= ck**e
    d4 += term

v = dict(zip(NAMES,ctx.gens()))
A = v['e4']+v['v5']+v['v6']
Bq = v['f10']+v['y11']+v['y12']
L = v['a']+v['c']+v['u2']+v['u3']+v['x8']+v['x9']
Q = (L+v['f7'])*A**3+(v['e1']*(v['v5']+v['v6'])+v['b4']*v['v6'])*A**2+\
    (L+v['e1'])*Bq**3+(v['f7']*(v['y11']+v['y12'])+v['b10']*v['y12'])*Bq**2+\
    A*Bq*(A+Bq)**2
cert = json.loads((INPUT/'certificate.json').read_text())
require(cert['names'] == list(NAMES) and cert['denominator'] == 1, 'Certificate convention mismatch')
require(all(len(e)==16 and sum(e)==4 and int(co)>0 for e,co in cert['Q']), 'Malformed multiplier')
certdict = {tuple(e):int(co) for e,co in cert['Q']}
require(len(certdict)==len(cert['Q']), 'Duplicate multiplier monomial')
require(Q == ctx.from_dict(certdict), 'Compressed multiplier differs from certificate')
require(len(Q)==266 and sum(Q.coeffs())==756, 'Multiplier statistics mismatch')

R = d4-Q*cs[6]
require(len(R)>0, 'Zero remainder cannot imply strictness')
min_coeff = min(R.coeffs())
require(min_coeff>0, 'Nonpositive remainder coefficient')
witness = tuple(1 if n=='u3' else 2 if n=='v6' else 3 if n=='x9' else 4 if n=='y12' else 0 for n in NAMES)
require(R[witness] == 1, 'Claimed strictness witness differs')
core = v['a']*v['c']*v['e4']*v['f10']*v['v6']*v['y12']
negative = ctx.from_dict({tuple(e):-int(co) for e,co in d4.terms() if co<0})
require(negative == core*Q, 'Negative-part interpretation false')
def stats(poly):
    coeffs = poly.coeffs()
    return {'terms':len(poly),'degree':int(poly.total_degree()),'negative_coefficients':sum(x<0 for x in coeffs),
            'minimum_coefficient':str(min(coeffs)),'coefficient_sum':str(sum(coeffs))}
receipt = {
    'verified': True, 'integrity_verified':True,
    'python':platform.python_version(),'sympy':s.__version__,'python_flint':flint.__version__,
    'matrix_route':'primary reaction columns S and partial derivatives Jrate -> full9x9 S*Jrate -> conservation embedding T -> P*J9*T',
    'characteristic_route':'Faddeev-LeVerrier traces with exact integer division and Cayley-Hamilton residual',
    'hurwitz_route':'generic 4x4 Hurwitz determinant generated by SymPy, evaluated in integer polynomial ring',
    'source_reduction_matches_frozen_matrix':True,'conservation_and_tangent_intertwining_checked':True,
    'coefficients':[stats(c) for c in cs],'Delta4':stats(d4),'Q':stats(Q),'R':stats(R),
    'strict_monomial':'u3*v6^2*x9^3*y12^4','strict_monomial_coefficient':int(R[witness]),
    'negative_part_equals_core_Q':True,'compressed_Q_equals_certificate':True,
    'generic_Delta4':str(generic.as_expr()),
    'optimization_mode_preserves_checks':True,
    'independent_script_sha256':sha256(Path(__file__).read_bytes()).hexdigest()
}
print(json.dumps(receipt,indent=2))
