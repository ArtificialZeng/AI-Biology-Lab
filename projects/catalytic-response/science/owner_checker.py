#!/usr/bin/env python3
"""Recompute general symbolic identities and one exact six-rate witness.

Numerical stress tests are supplementary and never certify the theorem.
Run with Python 3 and SymPy; writes evidence/verification.json.
"""
from pathlib import Path
import json
import math
import random
import sympy as s


def main():
    C, x = s.symbols('C x')
    a, b, c, d, e, f = s.symbols('a b c d e f', positive=True)
    D = a+c+d*C
    M = a+(b-a)*C
    N = c+(d-b-c)*C-d*C*C
    A, B = M/D, N/D
    T = a*d-(a+c)*(b+c)
    H = (e*C*D*x*x+(b*C*D-c*M)*x-f*N).expand()
    H_expanded = d*(b*x+e*x*x+f)*C*C + (
        a*(b+c)*x+(a+c)*e*x*x+f*(b+c-d))*C-c*(a*x+f)
    odeA = a*B+b*C-c*A-d*A*C
    odeB = d*A*C+e*x*C-a*B-f*B/x
    odeC = c*A+f*B/x-b*C-e*x*C
    S = T-2*d*(a+c)*C-d*d*C*C
    R = a*d+(a+c)*(a-b)
    identities = {
        'normalization': s.factor(A+B+C-1),
        'A_nullcline': s.factor(odeA),
        'H_expansion': s.factor(H-H_expanded),
        'H_is_minus_xD_C_ode': s.factor(H+x*D*odeC),
        'B_ode_follows': s.factor(odeA+odeB+odeC),
        'B_C_derivative': s.factor(s.diff(B,C)*D*D-S),
        'S_strict_decrease_identity': s.factor(s.diff(S,C)+2*d*D),
        'R_relation': s.factor(R-(a+c)**2-T),
        'B_division': s.factor(B-(1+(2*a+c-b)/d-(D+R/D)/d)),
        'A_division': s.factor(A-((b-a)/d+R/(d*D))),
    }
    assert all(v == 0 for v in identities.values())

    # Rational critical point, with all six rates strictly positive.
    rates = {a:s.Rational(1), b:s.Rational(1), c:s.Rational(1),
             d:s.Rational(9), e:s.Rational(1), f:s.Rational(3,5)}
    pt = {C:s.Rational(1,9), x:s.Rational(3)}
    exact_A = s.factor(A.subs(rates).subs(pt))
    exact_B = s.factor(B.subs(rates).subs(pt))
    assert exact_A == s.Rational(1,3) and exact_B == s.Rational(5,9)
    assert H.subs(rates).subs(pt) == 0
    assert S.subs(rates).subs(pt) == 0
    assert all(s.factor(v.subs(rates).subs(pt)) == 0 for v in (odeA,odeB,odeC))
    Hc = s.diff(H,C).subs(rates).subs(pt)
    Hx = s.diff(H,x).subs(rates).subs(pt)
    Cx = -Hx/Hc
    Bxx = s.diff(S,C).subs(rates).subs(pt)*Cx*Cx / D.subs(rates).subs(pt)**2
    assert Cx < 0 and Bxx < 0

    # Supplementary floating checks: random positive rates from 10^-3 to 10^3,
    # x from 10^-4 to 10^4. Direct ODE evaluation and Jacobian sensitivity.
    # These are discovery/implementation checks, not exhaustive evidence.
    rng = random.Random(20261027)
    max_residual = 0.0
    max_derivative_error = 0.0
    sign_checks = 0
    for _ in range(2000):
        av,bv,cv,dv,ev,fv = [10**rng.uniform(-3,3) for _ in range(6)]
        xv = 10**rng.uniform(-4,4)
        h2 = dv*(bv*xv+ev*xv*xv+fv)
        h1 = av*(bv+cv)*xv+(av+cv)*ev*xv*xv+fv*(bv+cv-dv)
        h0 = -cv*(av*xv+fv)
        disc = math.sqrt(h1*h1-4*h2*h0)
        Cn = -2*h0/(h1+disc) if h1 >= 0 else (-h1+disc)/(2*h2)
        Dn = av+cv+dv*Cn
        An = (av+(bv-av)*Cn)/Dn
        Bn = (cv+(dv-bv-cv)*Cn-dv*Cn*Cn)/Dn
        assert min(An,Bn,Cn) > 0 and abs(An+Bn+Cn-1)<1e-10
        fa = av*Bn+bv*Cn-cv*An-dv*An*Cn
        fb = dv*An*Cn+ev*xv*Cn-av*Bn-fv*Bn/xv
        fc = cv*An+fv*Bn/xv-bv*Cn-ev*xv*Cn
        scale = 1+av+bv+cv+dv+ev*xv+fv/xv
        max_residual=max(max_residual,abs(fa)/scale,abs(fb)/scale,abs(fc)/scale)
        # A,C-coordinate Jacobian, direct RHS partial wrt x.
        j11=-av-cv-dv*Cn; j12=bv-av-dv*An
        j21=cv-fv/xv; j22=-bv-ev*xv-fv/xv
        det=j11*j22-j12*j21
        gx=-fv*Bn/(xv*xv)-ev*Cn
        direct_Cx=-j11*gx/det
        direct_Ax=j12*gx/det
        direct_Bx=-direct_Ax-direct_Cx
        hc=2*h2*Cn+h1
        hx=dv*(bv+2*ev*xv)*Cn*Cn+(
            av*(bv+cv)+2*(av+cv)*ev*xv)*Cn-cv*av
        param_Cx=-hx/hc
        sn=av*dv-(av+cv)*(bv+cv)-2*dv*(av+cv)*Cn-dv*dv*Cn*Cn
        param_Bx=sn*param_Cx/(Dn*Dn)
        assert direct_Cx<0 and param_Cx<0 and det>0
        err=abs(direct_Bx-param_Bx)/(1+abs(direct_Bx)+abs(param_Bx))
        max_derivative_error=max(max_derivative_error,err)
        # Exclude cancellation-near-zero from supplementary sign agreement.
        if abs(param_Bx)>1e-12 and abs(direct_Bx)>1e-12:
            assert direct_Bx*param_Bx>0
            sign_checks+=1
    assert max_residual<1e-10
    assert max_derivative_error<1e-7
    result={
        'scope':'Exact identities in six symbolic rates; one exact positive rational critical witness; nonexhaustive supplementary floating checks.',
        'symbolic_identities':{k:str(v) for k,v in identities.items()},
        'exact_witness':{
            'rates':{str(k):str(v) for k,v in rates.items()},
            'x':'3','A':str(exact_A),'B':str(exact_B),'C':'1/9',
            'H_C':str(Hc),'H_x':str(Hx),'C_x':str(Cx),'B_xx':str(Bxx),
            'T':str(T.subs(rates)),
            'high_input_B':'1/2','overshoot':'1/18',
        },
        'numerical_checks':{
            'number':2000,'seed':20261027,
            'rate_range':'[10^-3,10^3] log-uniform independently',
            'x_range':'[10^-4,10^4] log-uniform',
            'max_scaled_ode_residual':max_residual,
            'max_scaled_sensitivity_discrepancy':max_derivative_error,
            'noncancelling_derivative_sign_agreements':sign_checks,
            'status':'floating-point evidence only, not proof or exhaustive classification',
        },
        'sympy_version':s.__version__,
    }
    out=Path(__file__).resolve().parent/'evidence'/'verification.json'
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
