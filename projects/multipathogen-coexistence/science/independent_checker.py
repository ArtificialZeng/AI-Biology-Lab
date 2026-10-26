#!/usr/bin/env python3
"""Read-only exact replay from a separately implemented source ODE.

Python 3 and SymPy required. No candidate import, file writes, network, or floats.
State ordering is per-host S,L_1,...,L_V,I_1,...,I_V,R.
"""
import json
import sympy as s


def check(ok, label):
    if not bool(ok):
        raise RuntimeError(label)


def parameters(contact):
    return dict(b=[s.Integer(1)] * 2, d=[s.Integer(1)] * 2,
                eps=s.ones(2, 3), gam=s.ones(2, 3), beta=contact)


def ode(par):
    hosts, strains = len(par['b']), len(par['beta'])
    S = s.symbols('x0:' + str(hosts))
    L = s.Matrix(hosts, strains, lambda p, v: s.Symbol(f'y{p}_{v}'))
    I = s.Matrix(hosts, strains, lambda p, v: s.Symbol(f'z{p}_{v}'))
    R = s.symbols('r0:' + str(hosts))
    states, rates, infection, infected = [], [], [], []
    for p in range(hosts):
        force = [sum(par['beta'][v][p, q] * I[q, v]
                     for q in range(hosts)) for v in range(strains)]
        states.extend([S[p]] + list(L.row(p)) + list(I.row(p)) + [R[p]])
        rates.extend([par['b'][p] - par['d'][p] * S[p] - S[p] * sum(force)])
        rates.extend([S[p] * force[v] - (par['eps'][p, v] + par['d'][p]) * L[p, v]
                      for v in range(strains)])
        rates.extend([par['eps'][p, v] * L[p, v] - (par['gam'][p, v] + par['d'][p]) * I[p, v]
                      for v in range(strains)])
        rates.extend([sum(par['gam'][p, v] * I[p, v] for v in range(strains)) - par['d'][p] * R[p]])
        infected.extend(list(L.row(p)) + list(I.row(p)))
        infection.extend([S[p] * force[v] for v in range(strains)] + [s.Integer(0)] * strains)
    rates = s.Matrix(rates)
    inf_rates = s.Matrix([rates[states.index(x)] for x in infected])
    zero = dict.fromkeys(states, s.Integer(0))
    zero.update({S[p]: par['b'][p] / par['d'][p] for p in range(hosts)})
    F = s.Matrix(infection).jacobian(infected).subs(zero)
    W = F - inf_rates.jacobian(infected).subs(zero)
    ngm = F * W.inv()
    return states, rates, S, L, I, R, ngm


def kernel(par, v):
    # The row-scaled equilibrium kernel and source column-scaled NGM are similar.
    c = s.diag(*[par['eps'][p, v] / ((par['eps'][p, v] + par['d'][p]) *
                 (par['gam'][p, v] + par['d'][p])) for p in range(2)])
    S0 = s.diag(*[par['b'][p] / par['d'][p] for p in range(2)])
    B = S0 * par['beta'][v] * c
    A = S0 * c * par['beta'][v]
    check(c.inv() * A * c == B, 'NGM equilibrium similarity')
    rho = s.factor((B.trace() + s.sqrt((B[0, 0] - B[1, 1]) ** 2 + 4 * B[0, 1] * B[1, 0])) / 2)
    return c * par['beta'][v], B, rho


def active_point(states, S, L, I, R, susceptible, infectious):
    point = dict.fromkeys(states, s.Integer(0))
    for p in range(2):
        point[S[p]] = susceptible
        for v in range(3):
            point[I[p, v]] = infectious[v]
            point[L[p, v]] = 2 * infectious[v]
        point[R[p]] = sum(infectious)
    return point


def replay():
    q = s.Rational
    t = s.Symbol('t')
    base = parameters([s.ones(2) * a for a in [1, 3, 4]])
    xs, f, S, L, I, R, ngm = ode(base)
    check(len(xs) == 16, 'state dimension')
    for p in range(2):
        block = xs[8*p:8*p+8]
        check(s.expand(sum(f[8*p:8*p+8]) - (1 - sum(block))) == 0, 'host conservation')
    ngchar = s.factor(ngm.charpoly(t).as_expr())
    check(s.expand(ngchar - t**9 * (t-q(1,2)) * (t-q(3,2)) * (t-2)) == 0, 'full infected next-generation spectrum')
    r0 = [kernel(base, v)[2] for v in range(3)]
    check(r0 == [q(1,2), q(3,2), s.Integer(2)], 'basic reproduction numbers')
    identities = []
    for v, lam in enumerate([1, 3, 4]):
        # Locate rows by named state; no use of candidate indexing or its helper.
        observed = sum(f[xs.index(L[p, v])] + 2*f[xs.index(I[p, v])] for p in range(2))
        expected = sum(I[p, v] for p in range(2)) * (lam*sum(S)-4)
        check(s.expand(observed-expected) == 0, 'strain balance identity')
        identities.append(str(s.factor(observed)))
    check(s.solve([3*t-4, 4*t-4], [t]) == [], 'incompatible susceptible sums')
    points = []
    for v in [1, 2]:
        ss = q(2, [1,3,4][v])
        weights = [s.Integer(0)]*3
        weights[v] = (1-ss)/4
        pt = active_point(xs, S, L, I, R, ss, weights)
        check(f.subs(pt) == s.zeros(16,1), 'exact singleton ODE residual')
        check(all(a>=0 for a in pt.values()), 'nonnegative singleton')
        points.append(pt)
    fullchar = s.factor(f.jacobian(xs).subs(points[1]).charpoly(t).as_expr())
    wanted = (t+1)**3*(t+2)**7*(t+3)*(t*t+4*t+1)*(t**3+6*t*t+8*t+4)
    check(s.expand(fullchar-wanted) == 0, 'full 16-dimensional Jacobian spectrum')
    # Positive coefficients and a1*a2>a3 are the cubic Hurwitz conditions.
    check(6*8>4, 'Hurwitz cubic')
    tie = parameters([s.ones(2)*a for a in [1,4,4]])
    tx, tf, ts, tl, ti, tr, _ = ode(tie)
    tied = active_point(tx, ts, tl, ti, tr, q(1,2), [0,q(1,16),q(1,16)])
    check(tf.subs(tied) == s.zeros(16,1), 'positive tied-strain control')
    fullrank = parameters([s.Matrix([[1,q(1,2)],[q(1,2),1]]),
                           s.Matrix([[4,1],[1,4]]), s.Matrix([[5,2],[2,5]])])
    _, _, ss, _, _, _, n2 = ode(fullrank)
    frho = [kernel(fullrank,v)[2] for v in range(3)]
    check(frho == [q(3,8),q(5,4),q(7,4)], 'full-rank reproduction numbers')
    check(all(B.det()!=0 and all(e>0 for e in B) for B in fullrank['beta']), 'full-rank positivity')
    d2 = s.factor((s.eye(2)-s.diag(*ss)*fullrank['beta'][1]/4).det())
    d3 = s.factor((s.eye(2)-s.diag(*ss)*fullrank['beta'][2]/4).det())
    T,U = s.symbols('T U')
    eliminant = s.solve([1-T+q(15,16)*U,1-q(5,4)*T+q(21,16)*U],[T,U])
    check(eliminant == {T:q(8,3), U:q(16,9)}, 'full-rank necessary equilibrium equations')
    check(eliminant[T] > 2, 'violates total susceptible feasibility')
    check(s.expand(d3-q(5,4)*d2 - (q(9,64)*ss[0]*ss[1]-q(1,4))) == 0, 'full-rank determinant elimination')
    # Reviewer-selected unequal rates, host sizes, and non-symmetric contacts.
    held = dict(b=[q(21,20),q(19,20)], d=[q(9,10),q(11,10)],
                eps=s.Matrix([[q(11,10),q(9,10),q(21,20)],
                              [q(19,20),q(23,20),q(17,20)]]),
                gam=s.Matrix([[q(101,100),q(99,100),q(103,100)],
                              [q(97,100),q(102,100),q(98,100)]]),
                beta=[s.Matrix([[1,q(51,100)],[q(49,100),q(99,100)]]),
                      s.Matrix([[q(401,100),q(99,100)],[q(101,100),q(399,100)]]),
                      s.Matrix([[q(499,100),q(203,100)],[q(197,100),q(501,100)]])])
    hs = [kernel(held,v) for v in range(3)]
    hrho = [a[2] for a in hs]
    margin = hs[2][0] - hs[1][0]
    check(all(a>0 for a in margin), 'heterogeneous strict effective-kernel dominance')
    check(hrho[0]<1 and hrho[1]>1 and hrho[2]>1, 'heterogeneous DFE inequalities')
    hx,hf,hS,hL,hI,hR,_ = ode(held)
    for p in range(2):
        check(s.expand(sum(hf[8*p:8*p+8])-held['b'][p]+held['d'][p]*sum(hx[8*p:8*p+8])) == 0,
              'heterogeneous host conservation')
    for v in range(3):
        c = [held['eps'][p,v]/((held['eps'][p,v]+held['d'][p])*(held['gam'][p,v]+held['d'][p])) for p in range(2)]
        for p in range(2):
            lzero = {(hL[p,v]): (held['gam'][p,v]+held['d'][p])*hI[p,v]/held['eps'][p,v]}
            reduced = s.factor(hf[hx.index(hL[p,v])].subs(lzero)*c[p])
            expected = hS[p]*sum(hs[v][0][p,j]*hI[j,v] for j in range(2))-hI[p,v]
            check(s.expand(reduced-expected)==0, 'general equilibrium kernel derived from ODE')
    return dict(schema='independent-source-replay-v1',checks_passed=True,
                arithmetic='exact rational and symbolic',state_dimension=16,
                simple_R0=list(map(str,r0)),full_infected_NGM_charpoly=str(ngchar),
                simple_balance_identities=identities,winning_16D_charpoly=str(fullchar),
                fullrank_R0=list(map(str,frho)),fullrank_determinants=[str(d2),str(d3)],
                necessary_fullrank_S_sum=str(eliminant[T]),
                heldout_parameters={k:([str(x) for x in v] if k in ['b','d'] else
                    [[str(x) for x in row] for row in v.tolist()] if k in ['eps','gam'] else
                    [[[str(x) for x in row] for row in b.tolist()] for b in v]) for k,v in held.items()},
                heldout_R0=list(map(str,hrho)),
                heldout_effective_kernel_margin=[[str(x) for x in row] for row in margin.tolist()],
                conclusion='No equilibrium can have both supercritical strains endemic; PF dominance applies on an open neighborhood.',
                provenance='Fresh AI review; no trusted runner job_id or human peer review asserted.')


if __name__ == '__main__':
    print(json.dumps(replay(), sort_keys=True, indent=2))
