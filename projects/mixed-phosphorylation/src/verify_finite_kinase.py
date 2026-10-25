"""A second exact formulation of the finite-kinase certificate.

Uses the published nine ODEs directly, instead of reaction-arrow construction.
Checks the continuous characteristic identity and every stored Bernstein
expansion. This is a reproducibility check in the researcher context, not a
fresh mathematical referee.
"""
from pathlib import Path
from fractions import Fraction
from math import comb
import hashlib
import json
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]
x=s.Matrix(s.symbols('x1:10'))
z=s.Symbol('z')
x1,x2,x3,x4,x5,x6,x7,x8,x9=x
F=s.Matrix([-x1*x2+x3+100*x9,
            -x1*x2+x3+x4,
            x1*x2-2*x3,
            x3-x4,
            x4-100*x5*x6+x7,
            -100*x5*x6-3*x8*x6+s.Rational(19,10)*x7+101*x9,
            100*x5*x6-s.Rational(19,10)*x7,
            s.Rational(9,10)*x7-3*x6*x8+x9,
            3*x6*x8-101*x9])
algebra=json.loads((ROOT/'evidence/finite_kinase_algebra.json').read_text())
certificate=json.loads((ROOT/'evidence/finite_kinase_sign_certificate.json').read_text())
eq=s.Matrix([s.sympify(t) for t in algebra['equilibrium']])
assert all(s.cancel(t)==0 for t in F.subs(dict(zip(x,eq)),simultaneous=True))
J9=F.jacobian(x).subs(dict(zip(x,eq)),simultaneous=True).applyfunc(s.cancel)
dm=DomainMatrix.from_Matrix(J9).to_field()
char9=[dm.domain.to_sympy(t) for t in dm.charpoly()]
char6=[s.sympify(t) for t in algebra['coefficients']]
assert len(char9)==10 and len(char6)==7
assert all(s.cancel(a-b)==0 for a,b in zip(char9,char6+[0,0,0]))
print('9D characteristic polynomial = lambda^3 * 6D polynomial (exact rational-function identity)',flush=True)

lo,hi=map(Fraction,certificate['z_interval_enclosure'])
results=[]
for item in certificate['certificates']:
    n=item['bernstein_degree']
    bb=list(map(Fraction,item['bernstein_coefficients']))
    pp=s.Poly(s.sympify(item['signed_numerator']),z,domain=s.QQ)
    assert n==pp.degree() and len(bb)==n+1 and min(bb)>0
    power=[Fraction(str(c)) for c in reversed(pp.all_coeffs())]
    # A degree-n polynomial is identified by n+1 distinct rational points.
    # These exact evaluations use stdlib Fraction and do not reuse the
    # generator's Taylor-to-Bernstein conversion.
    for j in range(n+1):
        t=Fraction(j,n) if n else Fraction(0)
        zz=lo+(hi-lo)*t
        lhs=sum(c*zz**i for i,c in enumerate(power))
        rhs=sum(bb[k]*comb(n,k)*t**k*(1-t)**(n-k) for k in range(n+1))
        assert lhs==rhs
    # Check the advertised rational function is the signed numerator over
    # exactly the serialized denominator, including its sign.
    row=algebra['hurwitz'][item['index']-1]
    determinant=s.sympify(row['expression'])
    denominator=s.Rational(item['denominator_constant'])*s.prod(
        s.sympify(v['factor'])**v['power'] for v in item['denominator_factors'])
    den_sign=s.sign(denominator.subs(z,s.Rational(str(lo))))
    assert s.cancel(determinant-item['determinant_sign']*den_sign*pp.as_expr()/denominator)==0
    results.append({'hurwitz_index':item['index'],'degree':n,
                    'all_Bernstein_coefficients_strictly_positive':True,
                    'exact_Fraction_identity_points':n+1})
    print('Bernstein identity verified for Delta',item['index'],flush=True)

source_hash=hashlib.sha256((ROOT/'source.md').read_bytes()).hexdigest()
# The public source overview is a paraphrase. Its current hash is provenance,
# not authentication of the private source snapshot; the root runner checks
# the public archive manifest before executing this file.
out={'scope':'Reconstruction consistency and finite exact certificate verification; not independent peer review',
     'nine_dimensional_characteristic_identity':'identically lambda^3 times stored finite polynomial',
     'Bernstein_identity_checks':results,'public_source_overview_sha256':source_hash,
     'floating_point_sign_decisions':0}
(ROOT/'evidence/finite_kinase_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print('all exact checks passed',flush=True)
