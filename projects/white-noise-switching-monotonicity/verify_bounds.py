"""Exact rational checks for the integral bounds proved in proof.md.
This is reproducibility bookkeeping, not a proof assistant.
"""
import sys
sys.dont_write_bytecode = True
from fractions import Fraction as Q
from math import factorial
from archive_tools import run, require


def build_certificate():
    mu1,mu2,sigma1,sigma2,rho=Q(1),Q(101,1000),Q(1),Q(1,10),Q(0)
    V=sigma1**2-2*rho*sigma1*sigma2+sigma2**2
    q=(mu2-mu1+sigma1**2-rho*sigma1*sigma2)/V
    phi1=mu1-sigma1**2/2;phi2=mu2-sigma2**2/2
    require(V==Q(101,100) and q==Q(1,10), 'Source optimum and noise difference')
    require(phi1>phi2 and 0<q<Q(1,4), 'Ordered source hedging branch')
    require(mu1>0 and mu2>0 and sigma1>0 and sigma2>0 and -1<=rho<=1,
            'Admissible source parameters')
    lam1,lam2=Q(1,1000),Q(1,500)
    omega1,omega2=V*lam1,V*lam2
    require(0<omega1<omega2, 'Ordered positive switching rates')
    nu=(1-2*q)*(1+2*lam2)
    A=2*lam2*(1-q);B=2*lam2*q
    require((nu,A,B)==(Q(502,625),Q(9,2500),Q(1,2500)), 'Exact odds density parameters')
    require(Q(4,5)<nu<1 and 0<1-nu<Q(1,5), 'Real-power comparison exponents')
    require((Q(25,8))**5>1/A, 'A^(nu-1) < 25/8 by rational fifth-power comparison')
    require((1/B)**4>Q(500)**5, 'B^(-nu) > 500 by rational fifth-power comparison')
    e_partial=sum((Q(1,factorial(k)) for k in range(4)),Q(0))
    require(e_partial==Q(8,3), 'e strictly exceeds its four-term Taylor sum 8/3')
    exp_minus_one_upper=Q(3,8)
    delta=Q(1,100)
    require(A*B/delta==Q(9,62500), 'Endpoint exponential uniform bound')
    N_upper=Q(25,8)*(1/(1-nu)+exp_minus_one_upper)
    D_lower=Q(500)*(1-A*B/delta)*(1-delta-exp_minus_one_upper)
    require(N_upper==Q(134225,7872), 'Numerator upper bound')
    require(D_lower==Q(7686393,25000), 'Denominator lower bound')
    require(D_lower>18*N_upper>0, 'Exact strict moment bound ER(lambda2)<1/18')
    ER_upper=N_upper/D_lower
    separator=phi1-omega1*q
    require(separator==Q(499899,1000000), 'Rational fitness separator')
    F2_upper=phi1+omega2*((1-q)*ER_upper-q)
    require(F2_upper<separator, 'Exact finite switching-rate decrease')
    stringify=lambda values:{k:str(v) for k,v in values.items()}
    return {'schema':'white-noise-fixed-target-finite-pair-v1',
      'source_parameters':stringify({'mu1':mu1,'mu2':mu2,'sigma1':sigma1,'sigma2':sigma2,'rho':rho,'V':V,'q_fixed_q_infinity':q,'phi1':phi1,'phi2':phi2}),
      'rates':stringify({'lambda1':lam1,'lambda2':lam2,'omega1':omega1,'omega2':omega2}),
      'odds_density_lambda2':stringify({'nu':nu,'A':A,'B':B}),
      'rational_integral_bounds':stringify({'N_strict_upper':N_upper,'D_strict_lower':D_lower,'D_lower_minus_18N_upper':D_lower-18*N_upper,'ER_strict_upper':ER_upper,'one_eighteenth':Q(1,18)}),
      'fitness_bounds':stringify({'F1_strict_lower_separator':separator,'F2_strict_upper':F2_upper,'separator_minus_F2_upper':separator-F2_upper}),
      'strict_order':'F(omega2) < F2_upper < 499899/1000000 < F(omega1)',
      'analytic_dependencies':['Original-number Ito generator and stationary-law reconstruction','Positive odds density and finite moments','Stationary identity F=phi1+omega*((1-q)E[R]-q)','Integral rescalings and elementary exponential bounds in proof.md'],
      'arithmetic':'fractions.Fraction, no floating computation or assert gates',
      'limitations':['No Julia execution','Not a formal proof-assistant certificate','Not finite-population extinction or survival','Not optimized-q envelope or full ambiguous source conjecture','No novelty or human-review certification']}


if __name__ == '__main__':
    raise SystemExit(run(build_certificate, 'evidence/bounds-certificate.json',
                        'exact rational integral-bound certificate'))
