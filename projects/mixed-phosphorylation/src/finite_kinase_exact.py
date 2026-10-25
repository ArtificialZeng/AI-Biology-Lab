"""Exact finite-kinase stability reconstruction from ten reaction arrows.

Run only with python src/finite_kinase_exact.py.
No imports from the existing model or certificate implementation.
"""
from pathlib import Path
import json
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix

START = time.monotonic()
ROOT = Path(__file__).resolve().parents[1]
x = s.Matrix(s.symbols('x1:10'))
z, lam = s.symbols('z lambda')
K = s.Symbol('K')
rates = [s.Integer(1), 1, 1, 1, 100, 1, s.Rational(9, 10), 3, 1, 100]
# 0-based indices; reactants/products with their stoichiometric multiplicity.
arrows = [({0:1,1:1},{2:1}), ({2:1},{0:1,1:1}),
          ({2:1},{3:1}), ({3:1},{4:1,1:1}),
          ({4:1,5:1},{6:1}), ({6:1},{4:1,5:1}),
          ({6:1},{7:1,5:1}), ({7:1,5:1},{8:1}),
          ({8:1},{7:1,5:1}), ({8:1},{0:1,5:1})]
N = s.zeros(9,10)
v = s.zeros(10,1)
for j, (before,after) in enumerate(arrows):
    v[j] = rates[j] * s.prod(x[i]**power for i,power in before.items())
    for i,number in before.items(): N[i,j] -= number
    for i,number in after.items(): N[i,j] += number
C = s.Matrix([[0,1,1,1,0,0,0,0,0], [0,0,0,0,0,1,1,0,1],
              [1,0,1,1,1,0,1,1,1]])
assert C*N == s.zeros(3,10)
full = N*v
selected = [0,2,3,4,6,8]
q = s.Matrix([x[i] for i in selected])
eliminate = {x[1]: K-x[2]-x[3], x[5]:8-x[6]-x[8],
             x[7]:40-x[0]-x[2]-x[3]-x[4]-x[6]-x[8]}
reduced = full.extract(selected,[0]).subs(eliminate, simultaneous=True)
p = 8-s.Rational(1009,9)*z
w = s.Rational(19,9)*z/p
middle = s.Rational(101,3)*z/p
r = s.factor(40-s.Rational(2809,9)*z-s.Rational(322,9)*z/p)
free_kinase = s.factor(200*z/r)
K_of_z = s.factor(200*z+free_kinase)
eq9 = s.Matrix([r,free_kinase,100*z,100*z,w,p,s.Rational(1000,9)*z,middle,z])
assert all(s.cancel(t)==0 for t in full.subs(dict(zip(x,eq9)), simultaneous=True))
assert all(s.cancel(t)==0 for t in C*eq9-s.Matrix([K_of_z,8,40]))
eq = dict(zip(q,eq9.extract(selected,[0])))
J = reduced.jacobian(q).subs(eq, simultaneous=True).subs(K,K_of_z)
J = J.applyfunc(s.cancel)
print('reaction reconstruction and equilibrium identities exact', flush=True)
dm = DomainMatrix.from_Matrix(J).to_field()
coeffs = [s.factor(dm.domain.to_sympy(t)) for t in dm.charpoly()]
print('characteristic polynomial computed', time.monotonic()-START, flush=True)
out = {'method':'exact reconstruction from reaction arrows, rational-function characteristic polynomial',
       'species_order':[str(t) for t in x], 'reduced_species_indices_1based':[i+1 for i in selected],
       'stoichiometric_matrix':N.tolist(), 'conservation_matrix':C.tolist(),
       'free_S0':str(r), 'free_phosphatase':str(p), 'K_of_z':str(K_of_z),
       'equilibrium':[str(t) for t in eq9],
       'coefficients':[str(t) for t in coeffs], 'hurwitz':[]}
for size in range(1,7):
    H = s.Matrix(size,size,lambda i,j:coeffs[2*i-j+1] if 0<=2*i-j+1<=6 else 0)
    hd = DomainMatrix.from_Matrix(H).to_field()
    det = s.factor(hd.domain.to_sympy(hd.det()))
    out['hurwitz'].append({'index':size,'expression':str(det)})
    print('Delta',size,'computed',time.monotonic()-START,flush=True)
dest = ROOT/'evidence'/'finite_kinase_algebra.json'
dest.write_text(json.dumps(out,indent=2,default=str)+'\n')
print('saved',str(dest),flush=True)
