"""Exact quadratic Fourier normalization for the six-phase clock.

Outputs polynomial identities and explicit coefficient-majorant constants.
No finite computation is used to infer the global periodic-orbit question.
"""
import json
from pathlib import Path
from fractions import Fraction
import sympy as s

b=Path(__file__).resolve().parent
I=s.I; root=s.sqrt(3); zeta=(1+I*root)/2
c=s.symbols('c0:6'); A,B=s.symbols('A B')
lambdas=[s.simplify(zeta**(-k)-1-s.cos(2*k*s.pi/3)) for k in range(6)]
F=[]
for k in range(6):
    quadratic=sum(-s.cos(2*j*s.pi/3)*c[(k-j)%6]*c[j] for j in range(6))
    F.append(s.expand(lambdas[k]*c[k]+quadratic))
omega=root/2
H=s.simplify(s.Rational(1,2)/(2*lambdas[1]-lambdas[2]))
g=s.simplify(H-s.Rational(1,2))
assert s.expand_complex(H)==s.Rational(2,7)+I*root/7
assert s.re(g)==-s.Rational(3,14)

# Real local coordinates: A=ar+i*ai; stable correction Y=(u,v,br,bi).
ar,ai,u,v,br,bi=s.symbols('ar ai u v br bi',real=True)
AA=ar+I*ai; AB=ar-I*ai; yy=br+I*bi; yb=br-I*bi
sub={c[0]:AA*AB+u,c[1]:AA,c[2]:H*AA**2+yy,
     c[3]:v,c[4]:s.conjugate(H)*AB**2+yb,c[5]:AB}
actual=[s.expand(s.expand_complex(f.subs(sub))) for f in F]
ap=actual[1]; abp=s.conjugate(ap)
R1=s.expand(ap-lambdas[1]*AA-g*AA**2*AB)
yp=[s.expand(actual[0]-(ap*AB+AA*abp)),
    s.expand(actual[3]),
    s.expand(actual[2]-2*H*AA*ap)]
stable_linear=[-u,-3*v,lambdas[2]*yy]
R2=[s.expand(a-l) for a,l in zip(yp,stable_linear)]
real_R1=[s.expand(s.re(R1)),s.expand(s.im(R1))]
real_R2=[s.expand(s.re(R2[0])),s.expand(s.re(R2[1])),
         s.expand(s.re(R2[2])),s.expand(s.im(R2[2]))]
variables=[ar,ai,u,v,br,bi]

def term_certificate(expr,center):
    coeff_sum=s.Rational(0); counts={}
    for powers,coeff in s.Poly(expr,*variables).terms():
        da=sum(powers[:2]); dy=sum(powers[2:])
        if center:
            assert (da>=1 and dy>=1) or dy>=2,(powers,coeff)
        else:
            assert (dy==0 and da>=3) or (dy>=1 and da>=1) or dy>=2,(powers,coeff)
        qa=s.expand(coeff).coeff(root,0); qb=s.expand(coeff).coeff(root,1)
        assert s.simplify(coeff-qa-qb*root)==0
        assert qa.is_Rational and qb.is_Rational
        coeff_sum+=abs(qa)+2*abs(qb)
        key=f'center_degree_{da}_stable_degree_{dy}'
        counts[key]=counts.get(key,0)+1
    return coeff_sum,counts

cert1=[term_certificate(e,True) for e in real_R1]
cert2=[term_certificate(e,False) for e in real_R2]
C=max(s.Integer(1),s.ceiling(sum(x[0] for x in cert1+cert2)))
K=max(s.Integer(1),s.ceiling(28*C*C/s.Integer(3)))
assert K>=28*C*C/s.Integer(3)
delta=s.Rational(1,100000)
assert delta<=s.Rational(1,16)
assert delta<=K/(8*C*(1+2*K))
assert delta**2<=s.Rational(3)/(56*K*C*C)

# Exact counterexample to the unrelated simple entropy-inequality route.
w=list(map(Fraction,['1/50','1/50','1/25','26','5/4','1/50']))
S=sum(w);Q=sum(w[i]*w[(i+2)%6] for i in range(6));R=sum(w[(i-1)%6]/w[i] for i in range(6))
assert R*Q-S*S<0
scale=S/Q
entropy_derivative=S*S/Q-R

out={'schema':'six-clock-local-algebra-v1','kind':'exact-polynomial-and-method-obstruction-certificates',
     'Fourier_convention':'w_i=sum_k c_k*exp(i*pi*k*i/3); c_(6-k)=conj(c_k)',
     'linear_eigenvalues':[str(s.simplify(z)) for z in lambdas],
     'H_quadratic_mode2':str(H),'center_cubic_coefficient_g':str(g),
     'Re_g':str(s.re(g)),'Im_g':str(s.im(g)),
     'real_center_remainders':[str(z) for z in real_R1],
     'real_stable_remainders':[str(z) for z in real_R2],
     'center_term_degrees':[x[1] for x in cert1],
     'stable_term_degrees':[x[1] for x in cert2],
     'coefficient_majorant_C':int(C),'quadratic_Lyapunov_weight_K':int(K),
     'certified_small_coordinate_radius':str(delta),
     'delta_requirements':['delta<=1/16','delta<=K/(8*C*(1+2*K))','delta^2<=3/(56*K*C^2)'],
     'analytic_payoff':'Using W=|A|^2+K||Y||^2 gives Wprime<=-3|A|^4/14-K||Y||^2/4 within the stated small-coordinate neighborhood; proof must justify the estimates.',
     'entropy_candidate_failure':{'profile':[str(z) for z in w],'RQ_minus_S_squared':str(R*Q-S*S),
       'scale_for_positive_entropy_derivative':str(scale),'positive_entropy_derivative':str(entropy_derivative),
       'scope':'Refutes this proposed entropy route only; not the original periodic-orbit question.'},
     'global_H6_status':'unresolved; no stable closed orbit certified'}
(b/'local-algebra.json').write_text(json.dumps(out,sort_keys=True,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['H_quadratic_mode2','center_cubic_coefficient_g','coefficient_majorant_C','quadratic_Lyapunov_weight_K','global_H6_status']},indent=2))
