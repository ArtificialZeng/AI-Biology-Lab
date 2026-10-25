"""Rational Bernstein certificates for a continuum of finite kinase totals.

This script consumes the independently reconstructed algebra, not a numerical
fit. Positive Bernstein coefficients prove positivity on the whole interval.
"""
from pathlib import Path
import json
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
z=s.Symbol('z')
data=json.loads((ROOT/'evidence/finite_kinase_algebra.json').read_text())
left=s.Rational(7,100)
right=s.Rational(70122383355,10**12)
Q=2834281*z*z-568386*z+25920
assert Q.subs(z,left)>0 and Q.subs(z,right)<0
assert s.diff(Q,z).subs(z,right)<0
assert 0<left<right<s.Rational(72,1009)

def bernstein_coefficients(poly,lo,hi):
    # p(lo+(hi-lo)*t) = sum a_j t^j, with a_j from Taylor's formula.
    poly=s.Poly(poly,z,domain=s.QQ)
    n=poly.degree()
    derivative=poly
    aa=[]
    for j in range(n+1):
        aa.append(derivative.eval(lo)*(hi-lo)**j/s.factorial(j))
        derivative=derivative.diff()
    return [sum(aa[j]*s.binomial(k,j)/s.binomial(n,j) for j in range(k+1))
            for k in range(n+1)]

certificates=[]
for row in data['hurwitz']:
    numerator,denominator=s.fraction(s.sympify(row['expression']))
    # On [left,alpha), Q>0 and 1009z-72<0. The complete factorization
    # printed here makes the denominator sign separately inspectable.
    coeff,factors=s.factor_list(denominator,z)
    denominator_sign=s.sign(coeff)
    denominator_factors=[]
    for f,power in factors:
        if s.degree(f,z)==1:
            assert f.subs(z,left)*f.subs(z,right)>0
            sign=s.sign(f.subs(z,left))
        else:
            assert s.degree(f,z)==2 and s.rem(f,Q,z)==0
            sign=s.sign(s.Poly(f,z).LC())
        denominator_sign*=sign**power
        denominator_factors.append({'factor':str(f),'power':int(power),'sign_on_physical_interval':int(sign)})
    desired=1 if row['index']<=4 else -1
    scaled_numerator=s.expand(desired*denominator_sign*numerator)
    bb=bernstein_coefficients(scaled_numerator,left,right)
    assert all(b>0 for b in bb), ('sign not certified',row['index'])
    certificates.append({'index':row['index'],'determinant_sign':desired,
        'denominator_constant':str(coeff),'denominator_factors':denominator_factors,
        'signed_numerator':str(scaled_numerator),'bernstein_degree':len(bb)-1,
        'bernstein_coefficients':[str(b) for b in bb]})
    print('Delta',row['index'],'certified',desired,'degree',len(bb)-1,flush=True)

threshold=s.factor(s.sympify(data['K_of_z']).subs(z,left))
assert threshold<23
out={'scope':'Table-1 rates, P=8,S=40, all K>=K(left), in particular all K>=23',
     'z_interval_enclosure':[str(left),str(right)],
     'physical_interval':'left <= z < alpha; alpha is the first positive root of Q',
     'Q':str(Q),'threshold_exact':str(threshold),'threshold_display_only':str(s.N(threshold,18)),
     'certificates':certificates,'routh_first_column_signs':[1,1,1,1,1,-1,1],
     'right_half_plane_root_count':2,
     'caveat':'This certifies equilibrium instability only, not periodic-orbit existence.'}
(ROOT/'evidence/finite_kinase_sign_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('threshold',threshold,'<23; exactly two RHP roots',flush=True)
