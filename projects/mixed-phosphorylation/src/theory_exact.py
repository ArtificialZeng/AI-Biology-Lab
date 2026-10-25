"""Exact algebraic certificate for the K=infinity equilibrium, P=8,S=40.

This is an ordinary SymPy calculation over Q, not an interval ODE proof and
not a novelty claim.  No floating-point number decides any assertion below.
Run: python src/theory_exact.py
"""
from pathlib import Path
import json
import sympy as sp

z, lam = sp.symbols('z lambda')
Q = 2834281*z**2 - 568386*z + 25920
lo = sp.Rational(70122383354, 10**12)
hi = sp.Rational(70122383355, 10**12)
assert Q.subs(z, lo) > 0 > Q.subs(z, hi)
assert sp.diff(Q,z).subs(z, hi) < 0
assert 0 < lo < hi < sp.Rational(72,1009)

P, S = sp.Rational(8), sp.Rational(40)
p = P-sp.Rational(1009,9)*z
w = sp.Rational(19,9)*z/p
v = sp.Rational(101,3)*z/p

# Slow coordinates (a,b,w,u,z)=(x3,x4,x5,x7,x9).
a,b,W,u,Z = sp.symbols('a b W u Z')
pp = P-u-Z
vv = S-a-b-W-u-Z
F = sp.Matrix([-a+100*Z, a-b, b-100*W*pp+u,
               100*W*pp-sp.Rational(19,10)*u,
               3*pp*vv-101*Z])
eq = {a:100*z,b:100*z,W:w,u:sp.Rational(1000,9)*z,Z:z}

def field_reduce(expr):
    """Reduce a rational function exactly in Q[z]/(Q)."""
    num, den = sp.fraction(sp.cancel(expr))
    ni = sp.Poly(num,z,domain=sp.QQ).rem(sp.Poly(Q,z))
    di = sp.invert(sp.Poly(den,z,domain=sp.QQ),sp.Poly(Q,z))
    return ni.mul(di).rem(sp.Poly(Q,z)).as_expr()

assert all(field_reduce(f.subs(eq)) == 0 for f in F)
J = F.jacobian([a,b,W,u,Z]).subs(eq)
coeffs = J.charpoly(lam).all_coeffs()
root = min(sp.solve(Q,z))

def certify(expr, name, expected):
    reduced = field_reduce(expr)
    # Reduced polynomial is affine; its range on [lo,hi] is exactly between
    # its endpoint values. The unique Q root in that interval therefore has
    # the same certified sign.
    assert sp.Poly(reduced,z).degree() <= 1
    endpoints = [sp.factor(reduced.subs(z,t)) for t in (lo,hi)]
    if expected > 0:
        assert min(endpoints)>0
    else:
        assert max(endpoints)<0
    return {'name':name,'reduced_mod_Q':str(reduced),
            'interval_endpoints':[str(x) for x in endpoints],
            'certified_sign':expected,
            'display_only_approx':str(sp.N(expr.subs(z,root),18))}

certs = [certify(c,f'c{i}',1) for i,c in enumerate(coeffs[1:],1)]
for size in range(1,5):
    H = sp.Matrix([[coeffs[2*i-j+1] if 0<=2*i-j+1<=5 else 0
                    for j in range(size)] for i in range(size)])
    det = H.det(method='domain-ge')
    certs.append(certify(det,f'Delta{size}',-1 if size==4 else 1))

# Check independently the slow-manifold invariance expansion through O(eps).
eps,y=sp.symbols('eps y')
q=sp.Matrix([a,b,W,u,Z])
h0=a+100*Z
A=a+b
Lh0=(sp.Matrix([h0]).jacobian(q)*F)[0]
h1=sp.expand(A*h0-Lh0)
F1=sp.Matrix([-Lh0,0,0,0,-3*pp*h0])
raw=sp.Matrix([y*(1-eps*A)-2*a,a-b,b-100*W*pp+u,
               100*W*pp-sp.Rational(19,10)*u,
               3*pp*(vv-eps*y)-101*Z])
assert all(sp.expand(raw[i].subs(y,h0+eps*h1)-F[i]-eps*F1[i]).coeff(eps,1)==0
           for i in range(5))
resid=sp.expand(-(h0+eps*h1)+eps*A*(h0+eps*h1)+a+100*Z
                -eps*(sp.Matrix([h0+eps*h1]).jacobian(q)*raw.subs(y,h0+eps*h1))[0])
assert resid.coeff(eps,0)==0 and resid.coeff(eps,1)==0

out={'scope':'Table-1 rates; P=8,S=40; limit K=infinity',
     'Q':str(Q),'root_interval':[str(lo),str(hi)],
     'root_count_argument':'Q changes + to - and Q_prime<0 on interval',
     'certificates':certs,
     'routh_first_column_signs':[1,1,1,1,-1,1],
     'right_half_plane_root_count':2,
     'h0':str(h0),'h1':str(h1),'F1':[str(x) for x in F1],
     'caveat':'Continuity gives existence of K0; no explicit finite K0 or periodic-orbit certificate here.'}
dest=Path(__file__).resolve().parents[1]/'evidence'/'theory_exact_certificate.json'
dest.write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps({'root_interval':out['root_interval'],
 'signs':[(x['name'],x['certified_sign'],x['display_only_approx']) for x in certs],
 'RHP_roots':2,'slow_manifold_order1_identity':'exactly verified',
 'certificate':str(dest)},indent=2))
