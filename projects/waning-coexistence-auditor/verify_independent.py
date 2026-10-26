#!/usr/bin/env python3
"""Read-only referee checks. No import of the proposer checker or its formulas.

Pinned to SymPy 1.13.3. This script only emits JSON during ordinary verification.
"""
import hashlib
import json
from pathlib import Path
import sympy as S
import sys
sys.dont_write_bytecode = True
from archive_tools import run


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def same(a, b, message):
    check(S.cancel(a - b) == 0, message)


def build_certificate():
    check(S.__version__ == '1.13.3', 'Install pinned SymPy 1.13.3')
    base = Path(__file__).resolve().parent
    digest = '3d8848b2461d9f4f6d3c0bb4733d5cbe0714a31367c6d0db5161bc7019342af9'
    names = 's i1 i21 i2 i12 r1 r2 r12'.split()
    x = S.symbols(' '.join(names), positive=True)
    at = dict(zip(names, range(8)))
    s, i1, i21, i2, i12, r1, r2, r12 = x
    L, mu, g1, g2, b1, b2, e1, e2, sig1, sig2 = S.symbols(
        'Lambda mu gamma1 gamma2 beta1 beta2 eta1 eta2 sigma1 sigma2', positive=True)
    t1, t2, t12 = S.symbols('theta1 theta2 theta12', nonnegative=True)
    v1, v2 = g1 + mu, g2 + mu
    f = S.zeros(8, 1)
    reactions = []

    def add(substrates, products, rate):
        nonlocal f
        net = S.zeros(8, 1)
        for name in substrates:
            net[at[name]] -= 1
        for name in products:
            net[at[name]] += 1
        f += rate * net
        reactions.append((substrates, products))

    add([], ['s'], L)
    for name, value in zip(names, x):
        add([name], [], mu * value)
    # Primary and secondary infection, each separately catalyzed by primary/secondary cases.
    for catalyst, infectivity, transmission, susceptible, infected in [
        ('i1', 1, b1, 's', 'i1'), ('i21', e1, b1, 's', 'i1'),
        ('i2', 1, b2, 's', 'i2'), ('i12', e2, b2, 's', 'i2'),
        ('i1', 1, b1*sig1, 'r2', 'i21'), ('i21', e1, b1*sig1, 'r2', 'i21'),
        ('i2', 1, b2*sig2, 'r1', 'i12'), ('i12', e2, b2*sig2, 'r1', 'i12')]:
        add([susceptible, catalyst], [infected, catalyst],
            transmission * infectivity * x[at[susceptible]] * x[at[catalyst]])
    for start, end, coefficient in [('i1','r1',g1), ('i21','r12',g1),
                                  ('i2','r2',g2), ('i12','r12',g2),
                                  ('r1','s',t1), ('r2','s',t2), ('r12','s',t12)]:
        add([start], [end], coefficient * x[at[start]])
    check(len(reactions) == 24, 'Reaction count')
    same(sum(f), L - mu*sum(x), 'Mass balance')

    # Independently read Eq(6), rather than treating the reaction list as its own check.
    J1, J2 = i1+e1*i21, i2+e2*i12
    published = [L-mu*s-b1*J1*s-b2*J2*s+t1*r1+t2*r2+t12*r12,
                 b1*J1*s-v1*i1, b1*sig1*J1*r2-v1*i21,
                 b2*J2*s-v2*i2, b2*sig2*J2*r1-v2*i12,
                 g1*i1-(b2*sig2*J2+mu+t1)*r1,
                 g2*i2-(b1*sig1*J1+mu+t2)*r2,
                 g1*i21+g2*i12-(mu+t12)*r12]
    for a, b in zip(f, published):
        same(a, b, 'Source Eq(6) mismatch')

    # Derive amplitudes directly from the two recovery equations after the infected kernels.
    a, b = S.symbols('a b')
    kernel = {i1:a*s, i21:a*sig1*r2, i2:b*s, i12:b*sig2*r1,
              b1:v1/(s+e1*sig1*r2), b2:v2/(s+e2*sig2*r1)}
    fk = [S.cancel(z.subs(kernel, simultaneous=True)) for z in f]
    M, rhs = S.linear_eq_to_matrix([fk[5], fk[6]], [a, b])
    amp = M.inv() * rhs
    D = S.factor(M.det())
    same(D, g1*g2*s*s - sig1*sig2*v1*v2*r1*r2, 'Cramer determinant')
    for k in range(4):
        same(fk[1+k], 0, 'Infected kernel equation')
    rr12 = S.solve(fk[7].subs({a:amp[0], b:amp[1]}), r12)[0]
    reconstructed = {**kernel, a:amp[0], b:amp[1], r12:rr12}
    # Apply kernel first, then amplitudes, avoiding ambiguous simultaneous nested substitution.
    rest = {a:amp[0], b:amp[1], r12:rr12}
    for k in range(1, 8):
        same(fk[k].subs(rest, simultaneous=True), 0, 'Recovery balance')
    population = sum(z.subs(kernel, simultaneous=True).subs(rest, simultaneous=True) for z in x)
    same(fk[0].subs(rest, simultaneous=True), L-mu*population, 'Reduced population equation')

    p, q = S.symbols('p q', positive=True)
    curve = {r1:(q-s)/(e2*sig2), r2:(p-s)/(e1*sig1)}
    N = population.subs(curve, simultaneous=True)
    Np = S.cancel(N.subs(s,p))
    expected_Np = p+(v1+t1)*(q-p)/(g1*e2*sig2)
    same(Np, expected_Np, 'Right endpoint population')
    inv = (p+e2*sig2*g1*(L/mu-p)/(v1+t1))/q
    same(L-mu*Np, mu*(v1+t1)*q/(g1*e2*sig2)*(inv-1), 'Endpoint invasion sign')
    polynomial = S.cancel(D*(mu*population-L))
    same(polynomial.subs(curve, simultaneous=True).subs(s,p),
         -g2*p*p*mu*(v1+t1)*q/(e2*sig2)*(inv-1), 'Source PG endpoint factor')

    # Derive the strict-population lower bound as an identity modulo r1 balance.
    excess = sum(x)-s-(v1+t1)*r1/g1
    positive_remainder = i21+i2+i12+r2+r12 + b2*sig2*J2*r1/g1
    same(excess-positive_remainder, f[5]/g1, 'Converse identity')

    cert = json.loads((base/'evidence/exact-population-certificate.json').read_text())
    lookup = {str(z):z for z in [*x,L,mu,g1,g2,b1,b2,e1,e2,sig1,sig2,t1,t2,t12]}
    values = {lookup[key]:S.Rational(value)
              for block in ['parameters','state'] for key,value in cert['witness'][block].items()}
    check(all(values[z]>0 for z in x), 'Witness positivity')
    for z in f:
        same(z.subs(values), 0, 'Witness equilibrium')
    determinant = S.factor(f.jacobian(x).subs(values).det())
    same(determinant, S.Rational(416,243), 'Witness nondegeneracy')
    w_p = (v1/b1).subs(values)
    w_q = (v2/b2).subs(values)
    w_s0 = (L/mu).subs(values)
    w_inv = inv.subs(values).subs({p:w_p,q:w_q})
    same(w_inv,S.Rational(61,72), 'Witness invasion')
    same(sum(z.subs(values) for z in x), w_s0, 'Witness population')
    return {'ok':True,'snapshot_digest':digest,'sympy':S.__version__,
                      'reaction_count':len(reactions),'symbolic_checks':[
                      'independent reaction/source equation reconstruction','demographic mass identity',
                      'matrix-derived infected amplitudes and all seven balances',
                      'reduced mass equation','right endpoint and invasion factor',
                      'full-waning converse identity'],
                      'witness':{'p':str(w_p),'q':str(w_q),'s0':str(w_s0),
                                 'invasion':str(w_inv),'jacobian_determinant':str(determinant)},
                      'writes':'stdout only','limits':'Analytic positivity, IVT, quantifiers, and novelty need the accompanying manual report.'}


if __name__ == '__main__':
    raise SystemExit(run(build_certificate, 'evidence/independent-certificate.json', 'independent matrix-derived reconstruction'))
