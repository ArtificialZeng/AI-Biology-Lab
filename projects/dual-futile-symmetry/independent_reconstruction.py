#!/usr/bin/env python3
"""Independent reconstruction from stoichiometry and rate derivatives.

This file does not import the proposed verifier or symbolic derivation.
The reaction order and species names follow arXiv:2608.27081v2 (1.1)-(1.2).
"""
import json
from pathlib import Path
import sympy as s

species = ('A','B','C','E','F','U','V','X','Y')
vectors = [
    {'A':-1,'E':-1,'U':1}, {'A':1,'E':1,'U':-1},
    {'B':1,'E':1,'U':-1}, {'B':-1,'E':-1,'V':1},
    {'B':1,'E':1,'V':-1}, {'C':1,'E':1,'V':-1},
    {'C':-1,'F':-1,'X':1}, {'C':1,'F':1,'X':-1},
    {'B':1,'F':1,'X':-1}, {'B':-1,'F':-1,'Y':1},
    {'B':1,'F':1,'Y':-1}, {'A':1,'F':1,'Y':-1},
]
N = s.Matrix([[v.get(x,0) for v in vectors] for x in species])
names = 'a b4 b10 c e1 e4 f7 f10 u2 u3 v5 v6 x8 x9 y11 y12'
p = dict(zip(names.split(), s.symbols(names, positive=True)))
gradients = [
    {'A':p['a'],'E':p['e1']}, {'U':p['u2']}, {'U':p['u3']},
    {'B':p['b4'],'E':p['e4']}, {'V':p['v5']}, {'V':p['v6']},
    {'C':p['c'],'F':p['f7']}, {'X':p['x8']}, {'X':p['x9']},
    {'B':p['b10'],'F':p['f10']}, {'Y':p['y11']}, {'Y':p['y12']},
]
J = N*s.Matrix([[g.get(x,0) for x in species] for g in gradients])
# Inject a tangent vector of a compatibility class in order q1,q2,U,V,Y,X.
I = s.Matrix([
    [1,0,-1,0,0,0],[-1,-1,0,-1,-1,0],[0,1,0,0,0,-1],
    [0,0,-1,-1,0,0],[0,0,0,0,-1,-1],
    [0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,0,1],
    [0,0,0,0,1,0],
])
K = s.Matrix([
    [1,0,0,0,0,1,0,0,0],[0,0,1,0,0,0,0,1,0],
    [0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,1,0,0],
    [0,0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,1,0],
])
if K*I != s.eye(6):
    raise RuntimeError('bad compatibility coordinates')
Gfull = K*J*I
if J*I != I*Gfull:
    raise RuntimeError('stoichiometric class not invariant')

# Reduced involution, computed as the permutation q1<->q2,U<->X,V<->Y.
permutation = [1,0,5,4,3,2]
S = s.eye(6)[:,permutation]
equations = list(Gfull*S-S*Gfull)
constraints = s.linsolve(equations, list(p.values()))
equation_rank = s.linear_eq_to_matrix(equations, list(p.values()))[0].rank()
a,b,e,f,u,w,v,z = s.symbols('a b e f u w v z', positive=True)
sub = {p['a']:a,p['c']:a,p['b4']:b,p['b10']:b,
       p['e1']:e,p['f7']:e,p['e4']:f,p['f10']:f,
       p['u2']:u,p['x8']:u,p['u3']:w,p['x9']:w,
       p['v5']:v,p['y11']:v,p['v6']:z,p['y12']:z}
G = Gfull.subs(sub)
if G*S != S*G or equation_rank != 8:
    raise RuntimeError('equivariance restrictions')

# Choose invariant eigenvectors of S independently of the supplied P.
Peven = (S-s.eye(6)).nullspace()
Podd = (S+s.eye(6)).nullspace()
P = s.Matrix.hstack(*(Peven+Podd))
D = s.simplify(P.inv()*G*P)
if D[:3,3:] != s.zeros(3) or D[3:,:3] != s.zeros(3):
    raise RuntimeError('invariant subspaces failed')
even,odd = D[:3,:3],D[3:,3:]
lam = s.symbols('lam')
polys = [block.charpoly(lam).as_expr().expand() for block in (even,odd)]
h,j,k = a+e+u+w, 2*b+f+v+z, f+v+z
# Explicitly check the candidate's basis convention as well as the fresh nullspaces.
Pdisplay = s.Matrix([[1,0,0,1,0,0],[1,0,0,-1,0,0],
                     [0,1,0,0,1,0],[0,0,1,0,0,1],
                     [0,0,1,0,0,-1],[0,1,0,0,-1,0]])
Adisplay = s.Matrix([[0,-w,z],[a,-h,-e],[-2*b,-f,-j]])
Bdisplay = s.Matrix([[0,-w,-z],[a,-h,-e],[0,-f,-k]])
if G*Pdisplay != Pdisplay*s.diag(Adisplay,Bdisplay) or Pdisplay.det() == 0:
    raise RuntimeError('displayed GP=PD identity mismatch')
ca = [1,h+j,h*j-e*f+a*w+2*b*z,a*w*j+2*b*e*w+a*f*z+2*b*z*h]
cb = [1,h+k,h*k-e*f+a*w,a*(w*k-f*z)]
for poly, coefficients in zip(polys,(ca,cb)):
    if s.expand(poly-s.Poly.from_list(coefficients,lam).as_expr()) != 0:
        raise RuntimeError('cubic mismatch')
if s.expand(G.charpoly(lam).as_expr()-polys[0]*polys[1]) != 0:
    raise RuntimeError('six-dimensional factorization mismatch')
HA,HB = s.expand(ca[1]*ca[2]-ca[3]),s.expand(cb[1]*cb[2]-cb[3])
alpha,beta = a+u+w,2*b+v+z
# Derived by splitting alpha*j**2 and e*beta*h in the raw Routh expression.
HAgroup = (alpha*j*h +(u+w)*j**2 +a*(j**2-f*z)
           +e*beta*(a+u+e+j)+e*w*(v+z)+a*w*h+2*b*z*j)
HBgroup = (h+k)*(alpha*k+e*(v+z))+a*w*h+a*f*z
if s.expand(HA-HAgroup) or s.expand(HB-HBgroup):
    raise RuntimeError('fresh positive decomposition mismatch')
polyHA,polyHB = [s.Poly(x,a,b,e,f,u,w,v,z) for x in (HA,HB)]
if any(c<=0 for poly in (polyHA,polyHB) for c in poly.coeffs()):
    raise RuntimeError('Routh coefficient sign mismatch')
if s.expand(G.det()-ca[3]*cb[3]):
    raise RuntimeError('determinant mismatch')

# Positive threshold parameterization: w=f*z/(f+v+z), all other inputs positive.
Gzero=G.subs(w,f*z/(f+v+z))
# The cubic derivative at zero is c2B>0, checked symbolically above via its groups.
boundary_derivative = s.factor(s.diff(polys[1],lam).subs(lam,0).subs(w,f*z/k))
report = {
    'method':'N times formal rate-gradient matrix, then compatibility injection/projection',
    'source':'https://arxiv.org/html/2608.27081v2',
    'species_order':species,
    'full_reduced_matrix':str(Gfull),
    'equivariance_linear_rank':equation_rank,
    'equivariance_solution':str(constraints),
    'basis_determinant':str(P.det()),
    'candidate_basis_determinant':str(Pdisplay.det()),
    'parity_basis':str(P),
    'parity_blocks':str(D),
    'cubics':[str(s.factor(x)) for x in polys],
    'Routh_monomials':[len(x.terms()) for x in (polyHA,polyHB)],
    'minimum_Routh_coefficients':[str(min(x.coeffs())) for x in (polyHA,polyHB)],
    'fresh_even_positive_identity':str(HAgroup),
    'odd_positive_identity':str(HBgroup),
    'determinant_factorization':str(s.factor(ca[3]*cb[3])),
    'boundary_odd_derivative_at_zero':str(boundary_derivative),
    'checks_passed':True,
    'limitation':'Exact polynomial identities; analytical root count is in the referee report.',
}
print(json.dumps(report,indent=2,ensure_ascii=False))
