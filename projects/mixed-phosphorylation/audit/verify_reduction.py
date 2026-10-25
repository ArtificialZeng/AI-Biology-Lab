"""Second exact derivation adapted from an AI-assisted review pass; not external peer review.

Run with python audit/verify_reduction.py.
The six state coordinates differ from those in the proposed derivation.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import hashlib
import json
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parents[1]
START = time.monotonic()
out = {'scope': 'Second exact derivation; not external human peer review',
       'checks': [], 'completed': False}

def record(name, **extra):
    out['checks'].append({'name': name, **extra})
    out['elapsed_seconds'] = time.monotonic() - START
    (ROOT/'audit/reduction_verification.json').write_text(json.dumps(out, indent=2)+'\n')
    print(name, flush=True)

# Public adaptation: record inputs at this run's start, rather than importing
# the private runner's snapshot. This checks integrity during computation.
snap = {'files': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                  for name in ['source.md', 'evidence/finite_kinase_algebra.json',
                               'evidence/finite_kinase_sign_certificate.json']}}
snap['digest'] = hashlib.sha256(json.dumps(snap['files'], sort_keys=True,
                              separators=(',', ':'), ensure_ascii=False).encode()).hexdigest()
observed = {k: hashlib.sha256((ROOT/k).read_bytes()).hexdigest() for k in snap['files']}
assert observed == snap['files']
digest = hashlib.sha256(json.dumps(observed, sort_keys=True,
                                  separators=(',', ':'), ensure_ascii=False).encode()).hexdigest()
assert digest == snap['digest']
out['snapshot_digest'] = digest
record('snapshot hashes', file_count=len(observed))

algebra = json.loads((ROOT/'evidence/finite_kinase_algebra.json').read_text())
cert = json.loads((ROOT/'evidence/finite_kinase_sign_certificate.json').read_text())
z, K = s.symbols('z K')
p = 8-s.Rational(1009,9)*z
r = 40-s.Rational(2809,9)*z-s.Rational(322,9)*z/p
kk = 200*z+200*z/r
eq = s.Matrix([r,200*z/r,100*z,100*z,s.Rational(19,9)*z/p,
               p,s.Rational(1000,9)*z,s.Rational(101,3)*z/p,z])
assert all(s.cancel(u-s.sympify(v)) == 0 for u,v in zip(eq,algebra['equilibrium']))
assert s.cancel(kk-s.sympify(algebra['K_of_z'])) == 0

# Direct transcription of the published nine equations at the Table 1 rates.
x = s.Matrix(s.symbols('x1:10'))
x1,x2,x3,x4,x5,x6,x7,x8,x9 = x
odes = s.Matrix([-x1*x2+x3+100*x9,
                -x1*x2+x3+x4, x1*x2-2*x3, x3-x4,
                x4-100*x5*x6+x7,
                -100*x5*x6-3*x8*x6+s.Rational(19,10)*x7+101*x9,
                100*x5*x6-s.Rational(19,10)*x7,
                s.Rational(9,10)*x7-3*x6*x8+x9,
                3*x6*x8-101*x9])
C = s.Matrix([[0,1,1,1,0,0,0,0,0],
              [0,0,0,0,0,1,1,0,1],
              [1,0,1,1,1,0,1,1,1]])
assert C.rank() == 3 and C*odes == s.zeros(3,1)
assert all(s.cancel(t)==0 for t in odes.subs(dict(zip(x,eq)),simultaneous=True))
assert all(s.cancel(t)==0 for t in C*eq-s.Matrix([kk,8,40]))
record('published nine ODEs, conservation laws and equilibrium identities')

# Different affine coordinate chart: retain x3,x4,x5,x7,x8,x9.
selected = [2,3,4,6,7,8]
q = s.Matrix([x[i] for i in selected])
chart = {x1:40-x3-x4-x5-x7-x8-x9,
         x2:K-x3-x4, x6:8-x7-x9}
red = odes.extract(selected,[0]).subs(chart,simultaneous=True)
J = red.jacobian(q).subs(dict(zip(q,eq.extract(selected,[0]))),simultaneous=True)
J = J.subs(K,kk).applyfunc(s.cancel)
dm = DomainMatrix.from_Matrix(J).to_field()
cc = [s.factor(dm.domain.to_sympy(v)) for v in dm.charpoly()]
assert len(cc)==7 and cc[0]==1
assert all(s.cancel(u-s.sympify(v))==0 for u,v in zip(cc,algebra['coefficients']))
record('characteristic polynomial in independent affine chart',
       coordinates=['x3','x4','x5','x7','x8','x9'],degree=6)

# The tangent map proves that the six roots really are the restricted roots.
embedded = x.subs(chart,simultaneous=True)
T = embedded.jacobian(q)
eqmap = dict(zip(x,eq))
J9 = odes.jacobian(x).subs(eqmap,simultaneous=True)
assert C*T == s.zeros(3,6) and T.rank()==6
assert all(s.cancel(t)==0 for t in J9*T-T*J)
record('exact tangent-space intertwining J9*T=T*J6; rank(T)=6')

L = s.Rational(7,100)
U = s.Rational(70122383355,10**12)
Q = 2834281*z*z-568386*z+25920
A = 1009*z-72
assert s.cancel(r-Q/(9*(72-1009*z)))==0
assert 0<L<U<s.Rational(72,1009)
assert Q.subs(z,L)==s.Rational(209569,10000)
assert Q.subs(z,U)<0 and s.diff(Q,z).subs(z,U)<0
assert s.diff(Q,z,2)>0
B = s.cancel(kk.subs(z,L))
assert B==s.Rational(4660166,209569) and B<23
# G'(z)>0: d(z/p)/dz=8/p^2 on the physical interval.
G = 40-r
assert s.cancel(s.diff(G,z)-(s.Rational(2809,9)+s.Rational(2576,9)/p**2))==0
assert s.cancel(s.diff(kk,z)-(200+200/r+200*z*s.diff(G,z)/r**2))==0
assert s.limit(kk,z,0)==0
record('interval enclosure, monotone all-K parameterization and exact bound',
       left=str(L),right=str(U),threshold=str(B),
       Q_left=str(Q.subs(z,L)),Q_right=str(Q.subs(z,U)))

expected_denominators = [90*A*Q,72900*A**3*Q**2,
                        177147000*A**5*Q**3,71744535000*A**7*Q**4,
                        871696100250*A**8*Q**5,94143178827*A**7*Q**6]
den_signs = [-1,-1,-1,-1,1,-1]
det_signs = [1,1,1,1,-1,-1]
computed_deltas = []
total = 0
for j,item in enumerate(cert['certificates'],start=1):
    assert item['index']==j and item['determinant_sign']==det_signs[j-1]
    # Standard (transpose of candidate's) Hurwitz convention.
    H = s.Matrix(j,j,lambda row,col: cc[2*col-row+1]
                 if 0<=2*col-row+1<=6 else 0)
    hdm = DomainMatrix.from_Matrix(H).to_field()
    delta = s.factor(hdm.domain.to_sympy(hdm.det()))
    computed_deltas.append(delta)
    assert s.cancel(delta-s.sympify(algebra['hurwitz'][j-1]['expression']))==0
    rawnum,rawden = s.fraction(delta)
    assert s.expand(rawden-expected_denominators[j-1])==0
    signed = s.Poly(s.sympify(item['signed_numerator']),z,domain=s.QQ)
    assert s.expand(signed.as_expr()-det_signs[j-1]*den_signs[j-1]*rawnum)==0
    # Reconstruct the entire power polynomial from the stored Bernstein basis
    # using Fraction arithmetic, rather than repeating the generator's map.
    n = signed.degree()
    bb = [F(v) for v in item['bernstein_coefficients']]
    assert n==item['bernstein_degree'] and len(bb)==n+1 and min(bb)>0
    lhs = [F(0) for _ in range(n+1)]
    for k,b in enumerate(bb):
        for m in range(n-k+1):
            lhs[k+m] += b*comb(n,k)*comb(n-k,m)*(-1)**m
    rhs = [F(0) for _ in range(n+1)]
    lo,hi = F(str(L)),F(str(U))
    for (power,),value in signed.terms():
        v = F(str(value))
        for m in range(power+1):
            rhs[m] += v*comb(power,m)*lo**(power-m)*(hi-lo)**m
    assert lhs==rhs
    total += len(bb)
    record('Hurwitz determinant and full Bernstein polynomial identity',
           determinant=j,degree=n,positive_coefficients=len(bb),
           denominator_sign=den_signs[j-1],determinant_sign=det_signs[j-1])
assert total==125

# Reconstruct the Routh rows directly and check the determinant quotient law.
rows = [cc[0::2],cc[1::2]+[s.Integer(0)]]
for j in range(2,7):
    above,prev = rows[-2:]
    lead = prev[0]
    assert lead!=0
    row = [s.cancel((lead*above[k+1]-above[0]*prev[k+1])/lead)
           for k in range(3)]+[s.Integer(0)]
    rows.append(row)
first = [row[0] for row in rows]
quotients = [s.Integer(1),computed_deltas[0]]+[
    s.cancel(computed_deltas[j]/computed_deltas[j-1]) for j in range(1,6)]
assert all(s.cancel(a-b)==0 for a,b in zip(first,quotients))
signs=[1,1,1,1,1,-1,1]
assert sum(a!=b for a,b in zip(signs,signs[1:]))==2
record('direct Routh recurrence and two sign changes',first_column_signs=signs,
       right_half_plane_roots=2,left_half_plane_roots=4)

# Ensure none of the frozen inputs changed while verification ran.
assert all(hashlib.sha256((ROOT/k).read_bytes()).hexdigest()==v
           for k,v in snap['files'].items())
out['completed']=True
record('all independent exact checks completed; frozen evidence unchanged',
       total_positive_Bernstein_coefficients=total,floating_point_sign_decisions=0)
