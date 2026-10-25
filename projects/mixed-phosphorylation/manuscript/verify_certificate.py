"""Reconstruct and verify the manuscript's exact finite-kinase certificate.

Run from the project root with:
    python manuscript/verify_certificate.py
Writes manuscript/certificate_check.json; never edits archived evidence.
"""
from pathlib import Path
import json
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parents[1]
START = time.monotonic()
z = s.Symbol('z')
L = s.Rational(7, 100)
U = s.Rational(70122383355, 10**12)
p = 8-s.Rational(1009, 9)*z
r = 40-s.Rational(2809, 9)*z-s.Rational(322, 9)*z/p
u, w, v = 200*z/r, s.Rational(19, 9)*z/p, s.Rational(101, 3)*z/p
J = s.Matrix([
    [-u-r-2, -u-r, -u, -u, -u, -u],
    [1, -1, 0, 0, 0, 0],
    [0, 1, -100*p, 100*w+1, 0, 100*w],
    [0, 0, 100*p, -100*w-s.Rational(19, 10), 0, -100*w],
    [0, 0, 0, s.Rational(9, 10)+3*v, -3*p, 1+3*v],
    [0, 0, 0, -3*v, 3*p, -3*v-101]
]).applyfunc(s.cancel)

def exact_determinant(matrix):
    dm = DomainMatrix.from_Matrix(matrix).to_field()
    return dm.domain.to_sympy(dm.det())

dm = DomainMatrix.from_Matrix(J).to_field()
cc = [s.factor(dm.domain.to_sympy(c)) for c in dm.charpoly()]
algebra = json.loads((ROOT/'evidence/finite_kinase_algebra.json').read_text())
cert = json.loads((ROOT/'evidence/finite_kinase_sign_certificate.json').read_text())
assert all(s.cancel(c-s.sympify(a)) == 0
           for c, a in zip(cc, algebra['coefficients']))
A = 1009*z-72
Q = 2834281*z*z-568386*z+25920
D = [90*A*Q, 72900*A**3*Q**2, 177147000*A**5*Q**3,
     71744535000*A**7*Q**4, 871696100250*A**8*Q**5,
     94143178827*A**7*Q**6]
degrees = [4, 10, 17, 24, 30, 34]
exponents = [5, 14, 24, 34, 41, 47]
eta = [s.Rational(24, 1000), s.Rational(28, 10000),
       s.Rational(14, 10000), s.Rational(8, 10000),
       s.Rational(4, 10000), s.Rational(4, 10000)]
rows = []
for j in range(1, 7):
    H = s.Matrix(j, j, lambda a, b:
                 cc[2*b-a+1] if 0 <= 2*b-a+1 <= 6 else 0)
    delta = s.factor(exact_determinant(H))
    numerator, denominator = s.fraction(delta)
    assert s.expand(denominator-D[j-1]) == 0
    h = s.Poly((-1 if j <= 5 else 1)*numerator, z, domain=s.QQ)
    n = h.degree()
    assert n == degrees[j-1]
    aa = [s.diff(h.as_expr(), z, ell).subs(z, L)
          *(U-L)**ell/s.factorial(ell) for ell in range(n+1)]
    bb = [sum(aa[ell]*s.binomial(k, ell)/s.binomial(n, ell)
              for ell in range(k+1)) for k in range(n+1)]
    assert bb[0] > 10**exponents[j-1], (j, bb[0])
    assert all(b > eta[j-1]*bb[0] for b in bb)
    assert bb == [s.Rational(b) for b in cert['certificates'][j-1]['bernstein_coefficients']]
    rows.append({'index': j, 'degree': n, 'b0': str(bb[0]),
                 'b0_lower_bound_power_of_ten': exponents[j-1],
                 'minimum_ratio_exact': str(min(bb)/bb[0]),
                 'ratio_lower_bound': str(eta[j-1]),
                 'all_rational_checks_passed': True})
    print('exact table row', j, 'verified', flush=True)
out = {'method': 'reconstruction from manuscript six-by-six Jacobian; exact rational arithmetic',
       'sympy_version': s.__version__, 'rows': rows,
       'coefficient_count': sum(n+1 for n in degrees),
       'elapsed_seconds': time.monotonic()-START}
(ROOT/'manuscript/certificate_check.json').write_text(json.dumps(out, indent=2)+'\n')
print('125 exact Bernstein coefficients and all manuscript table bounds verified', flush=True)
