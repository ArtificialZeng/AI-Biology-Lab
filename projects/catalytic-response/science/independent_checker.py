#!/usr/bin/env python3
"""Independent portable exact checker. Python 3 stdlib only; no private inputs.

Polynomial identities and exact ODE/Jacobian tests certify the algebraic steps,
not the positivity/quantifier proof or novelty. See the separate referee report.
"""
from fractions import Fraction as Q
from itertools import product
import json

NVAR = 8


class P:
    def __init__(self, terms=0):
        self.t = {k: Q(v) for k, v in terms.items() if v} if isinstance(terms, dict) else ({(0,) * NVAR: Q(terms)} if terms else {})

    def __add__(self, other):
        other = other if isinstance(other, P) else P(other)
        out = self.t.copy()
        for k, v in other.t.items():
            out[k] = out.get(k, 0) + v
        return P(out)

    __radd__ = __add__

    def __neg__(self):
        return P({k: -v for k, v in self.t.items()})

    def __sub__(self, other):
        return self + (-other if isinstance(other, P) else -Q(other))

    def __rsub__(self, other):
        return P(other) - self

    def __mul__(self, other):
        other = other if isinstance(other, P) else P(other)
        out = {}
        for k, v in self.t.items():
            for l, w in other.t.items():
                m = tuple(u + z for u, z in zip(k, l))
                out[m] = out.get(m, 0) + v * w
        return P(out)

    __rmul__ = __mul__

    def diff(self, i):
        out = {}
        for k, v in self.t.items():
            if k[i]:
                l = list(k)
                l[i] -= 1
                out[tuple(l)] = v * k[i]
        return P(out)


def var(i):
    k = [0] * NVAR
    k[i] = 1
    return P({tuple(k): 1})


def exact_point(rates, state, x):
    a, b, c, d, e, f = map(Q, rates)
    A, B, C = map(Q, state)
    x = Q(x)
    assert min(a, b, c, d, e, f, A, B, C, x) > 0
    assert A + B + C == 1
    rhs = (a*B+b*C-c*A-d*A*C,
           d*A*C+e*x*C-a*B-f*B/x,
           c*A+f*B/x-b*C-e*x*C)
    assert rhs == (0, 0, 0)
    # Independent effective-generator nullspace via its three principal cofactors.
    # Columns are source states; the catalytic A->B rate is d*C.
    wA = a*b+a*e*x+b*f/x
    wB = d*C*(b+e*x)+c*e*x
    wC = c*(a+f/x)+d*C*f/x
    Z = wA+wB+wC
    assert (A,B,C) == (wA/Z,wB/Z,wC/Z)
    # Differentiate the original ODE after eliminating B by conservation.
    j11, j12 = -a-c-d*C, b-a-d*A
    j21, j22 = c-f/x, -b-e*x-f/x
    det = j11*j22-j12*j21
    g = -f*B/(x*x)-e*C
    Ax, Cx = j12*g/det, -j11*g/det
    Bx = -Ax-Cx
    assert det > 0 and Cx < 0
    assert j11*Ax+j12*Cx == 0
    assert j21*Ax+j22*Cx+g == 0
    p, D = a+c, a+c+d*C
    T = a*d-(a+c)*(b+c)
    S = T-2*d*p*C-d*d*C*C
    assert Bx == S*Cx/(D*D)
    assert (Bx > 0) == (S < 0)
    assert (Bx < 0) == (S > 0)
    return dict(Ax=Ax, Bx=Bx, Cx=Cx, det=det, T=T, S=S)


def main():
    if not __debug__:
        raise RuntimeError('Run with assertions enabled; Python -O is unsupported')
    a,b,c,d,e,f,C,x = [var(i) for i in range(NVAR)]
    p = a+c
    D = p+d*C
    M = a+(b-a)*C
    N = c+(d-b-c)*C-d*C*C
    T = a*d-p*(b+c)
    S = T-2*d*p*C-d*d*C*C
    H = e*C*D*x*x+(b*C*D-c*M)*x-f*N
    # All identities originate from the published three ODEs, never a resultant.
    ids = {
        'mass': M+N+C*D-D,
        'A_ODE_numerator': a*N+b*C*D-c*M-d*M*C,
        'C_ODE_numerator': x*c*M+f*N-x*b*C*D-e*x*x*C*D+H,
        'B_ODE_numerator': x*d*M*C+e*x*x*C*D-a*x*N-f*N-H,
        'quadratic': H-(d*(b*x+e*x*x+f)*C*C+(a*(b+c)*x+p*e*x*x+f*(b+c-d))*C-c*(a*x+f)),
        'B_derivative': N.diff(6)*D-d*N-S,
        'S_derivative': S.diff(6)+2*d*D,
        'input_monotonicity': x*H.diff(7)-H-e*C*D*x*x-f*N,
        'jacobian_numerator': H.diff(6)*D-d*H-(D*D*(b*x+e*x*x+f)-((b-a)*D-d*M)*(c*x-f)),
        'sensitivity_sign': D*D+(b-a)*D-d*M+S,
        'high_input_difference': p*N-c*D-C*(T-p*d*C),
    }
    for name, poly in ids.items():
        assert not poly.t, (name, poly.t)

    # Held-out exact equilibrium constructed from species fractions and ODEs:
    # choose a,b,d,A,B,C, then solve for c and f, rather than solve candidate H.
    checks = signs = 0
    sign_counts = {'positive': 0, 'negative': 0, 'zero': 0}
    for i,j,k,av,bv,dv,xv in product(range(1,5),range(1,5),range(1,5),(1,3),(1,4),(1,7,23),(1,3)):
        A,B,Cv = Q(i,i+j+k),Q(j,i+j+k),Q(k,i+j+k)
        cv = (av*B+bv*Cv-dv*A*Cv)/A
        if cv <= 0:
            continue
        ev = Q(1) + cv*A/(xv*Cv)
        fv = xv*((bv+ev*xv)*Cv-cv*A)/B
        out = exact_point((av,bv,cv,dv,ev,fv),(A,B,Cv),xv)
        label = 'positive' if out['Bx'] > 0 else 'negative' if out['Bx'] < 0 else 'zero'
        sign_counts[label] += 1
        checks += 1
    assert all(sign_counts.values()), sign_counts

    extreme_checks = 0
    for i,j,k,av,bv,dv,xv in product((1,10**6),(1,10**6),(1,10**6),
                                     (Q(1,10**6),Q(10**6)),
                                     (Q(1,10**6),Q(10**6)),
                                     (Q(1,10**6),Q(10**6)),
                                     (Q(1,10**9),Q(10**9))):
        Av,Bv,Cv = Q(i,i+j+k),Q(j,i+j+k),Q(k,i+j+k)
        cv = (av*Bv+bv*Cv-dv*Av*Cv)/Av
        if cv <= 0:
            continue
        ev = Q(1) + cv*Av/(xv*Cv)
        fv = xv*((bv+ev*xv)*Cv-cv*Av)/Bv
        exact_point((av,bv,cv,dv,ev,fv),(Av,Bv,Cv),xv)
        extreme_checks += 1

    # Independent exact critical witness, distinct from the owner's witness.
    rates = (Q(2),Q(3),Q(1),Q(14),Q(1),Q(4,3))
    state = (Q(3,7),Q(3,7),Q(1,7))
    w = exact_point(rates,state,Q(2))
    assert w['Bx'] == 0 and w['T'] == 16
    p,D = Q(3),Q(5)
    S_C = -2*Q(14)*D
    Bxx = S_C*w['Cx']*w['Cx']/(D*D)
    assert Bxx < 0
    assert state[1]-Q(1,3) == Q(2,21)

    # Exact equality case: ad=(a+c)(b+c), at a finite input.
    eq = exact_point((2,3,1,6,1,Q(20,19)),(Q(5,9),Q(19,63),Q(1,7)),2)
    assert eq['T'] == 0 and eq['Bx'] > 0
    # Adversarial approach to threshold from either side, at exact equilibria.
    threshold_checks = 0
    for denom in (10,100,10000):
        for delta in (-Q(1,denom),Q(0),Q(1,denom)):
            dv=Q(6)+delta
            for Cv in (Q(1,100000),Q(1,10)):
                Av=(Q(2)+Cv)/(Q(3)+dv*Cv)
                Bv=1-Av-Cv
                xv=Q(2)
                ev=Q(1)+Av/(xv*Cv)
                fv=xv*((3+ev*xv)*Cv-Av)/Bv
                exact_point((2,3,1,dv,ev,fv),(Av,Bv,Cv),xv)
                threshold_checks += 1
    result = {
        'checker':'independent Python 3 stdlib exact algebra, no private inputs',
        'polynomial_identities':sorted(ids),
        'held_out_exact_ODE_Jacobian_points':checks,
        'held_out_effective_generator_nullspace_points':checks,
        'additional_extreme_exact_ODE_nullspace_sensitivity_points':extreme_checks,
        'held_out_B_derivative_sign_counts':sign_counts,
        'near_threshold_exact_points':threshold_checks,
        'equality_finite_input_Bx':str(eq['Bx']),
        'independent_peak':{'rates':list(map(str,rates)), 'A_B_C':list(map(str,state)),
                            'x':'2','T':'16','Cx':str(w['Cx']), 'Bxx':str(Bxx), 'overshoot':'2/21'},
        'status':'all checks passed; report supplies quantified proof and scope limits'
    }
    print(json.dumps(result,sort_keys=True,separators=(',',':')))


if __name__ == '__main__':
    main()
