"""Independent N4 victim/donor two-colour checker; imports no primary mathematical code."""
import itertools
import json
import pathlib
import sympy as s

def load(p):
    return json.loads(p.read_text())

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def comps(n, bins):
    if bins == 1:
        yield (n,)
    else:
        for x in range(n + 1):
            for tail in comps(n - x, bins - 1):
                yield ((x,) + tail)

def move(x, i, j):
    y = list(x)
    y[i] -= 1
    y[j] += 1
    return tuple(y)

def build_certificate():
    require(s.__version__ == '1.13.3', 'SymPy 1.13.3 required')
    N = 4
    dossier = pathlib.Path(__file__).resolve().parent
    t = s.Symbol('t')
    gens = (t,)
    kill = (1, t, t + 1)

    def poly(x):
        return s.Poly(s.sympify(x, locals={str(g): g for g in gens}), *gens, domain=s.QQ)
    zero = poly(0)
    order = load(dossier / 'evidence/four-particle/state-ordering.json')
    pop = [tuple(x) for x in order['population']]
    founders = [(x['i'], tuple(x['d'])) for x in order['founder']]
    require(set(pop) == set(comps(N, 3)) and len(pop) == len(set(pop)), 'Independent polynomial/state condition 41 failed')
    require(set(founders) == set(itertools.product(range(3), comps(N - 1, 3))), 'Independent polynomial/state condition 42 failed')
    require(len(founders) == len(set(founders)), 'Independent polynomial/state condition 43 failed')
    px = load(dossier / 'evidence/four-particle/population-solution.json')
    hx = load(dossier / 'evidence/four-particle/founder-solution.json')
    pn = list(map(poly, px['numerators']))
    pd = poly(px['denominator'])
    hn = list(map(poly, hx['numerators']))
    hd = poly(hx['denominator'])
    p = dict(zip(pop, pn))
    h = dict(zip(founders, hn))
    kk = list(map(poly, kill))
    require(sum(pn, zero) == pd, 'Independent polynomial/state condition 50 failed')
    stat = {c: zero for c in pop}
    population_events = 0
    for c in pop:
        for i in range(3):
            if not c[i]:
                continue
            for j in range(3):
                if i == j:
                    continue
                rate = poly(c[i]) + kk[i] * s.Rational(c[i] * c[j], N - 1)
                mass = p[c] * rate
                stat[c] -= mass
                stat[move(c, i, j)] += mass
                population_events += 1
    require(all((v.is_zero for v in stat.values())), 'Independent polynomial/state condition 62 failed')
    require(N * h[0, (N - 1, 0, 0)] == hd, 'Independent polynomial/state condition 63 failed')
    for c in pop:
        total = zero
        for i in range(3):
            if c[i]:
                d = list(c)
                d[i] -= 1
                total += h[i, tuple(d)] * c[i]
        require(total == hd, 'Independent polynomial/state condition 70 failed')
    for i, d in founders:
        out = -kk[i] * h[i, d]
        for j in range(3):
            if j != i:
                out += h[j, d] - h[i, d]
        for victim in range(3):
            if not d[victim]:
                continue
            for dest in range(3):
                if victim != dest:
                    out += (h[i, move(d, victim, dest)] - h[i, d]) * d[victim]
            contribution = h[i, move(d, victim, i)] * s.Rational(2, N - 1) - h[i, d]
            for donor in range(3):
                donors = d[donor] - (victim == donor)
                if donors:
                    contribution += h[i, move(d, victim, donor)] * s.Rational(donors, N - 1)
            out += contribution * kk[victim] * d[victim]
        require(out.is_zero, 'Independent polynomial/state condition 88 failed')
    coloured = list(comps(N, 6))
    mixed = [z for z in coloured if 0 < sum(z[:3]) < N]

    def H(z):
        if sum(z[:3]) == 0:
            return zero
        if sum(z[:3]) == N:
            return hd
        c = tuple((z[i] + z[i + 3] for i in range(3)))
        total = zero
        for i in range(3):
            if z[i]:
                d = list(c)
                d[i] -= 1
                total += h[i, tuple(d)] * z[i]
        return total
    Hcache = {z: H(z) for z in coloured}
    residuals = []
    motion_events = copy_events = same_site_cross = 0
    for z in mixed:
        out = zero
        for victim in range(6):
            if not z[victim]:
                continue
            ci, site = divmod(victim, 3)
            for destsite in range(3):
                if site == destsite:
                    continue
                dest = 3 * ci + destsite
                out += (Hcache[move(z, victim, dest)] - Hcache[z]) * z[victim]
                motion_events += 1
            for donor in range(6):
                donors = z[donor] - (victim == donor)
                if not donors:
                    continue
                copy_events += 1
                if donor == victim:
                    continue
                if donor % 3 == site and donor // 3 != ci:
                    same_site_cross += 1
                out += (Hcache[move(z, victim, donor)] - Hcache[z]) * kk[site] * s.Rational(z[victim] * donors, N - 1)
        residuals.append(out)
    require(all((v.is_zero for v in residuals)), 'Independent polynomial/state condition 126 failed')
    C = zero
    for c in pop:
        for i, j, sgn in [(0, 2, 1), (2, 0, -1)]:
            if c[i]:
                d = list(c)
                d[i] -= 1
                C += p[c] * h[j, tuple(d)] * (sgn * c[i])
    cx = load(dossier / 'evidence/four-particle/symbolic-current.json')
    A = poly(cx['numerator'])
    D = poly(cx['denominator'])
    require(C * D == A * pd * hd, 'Independent polynomial/state condition 137 failed')
    result = {'N': N, 'population_states': len(pop), 'population_transition_terms': population_events, 'population_residuals_zero': True, 'normalization_and_anchor': True, 'direct_founder_Eq3_6_residuals_zero': len(founders), 'mixed_two_colour_states': len(mixed), 'colour_residuals_zero': True, 'motion_event_terms': motion_events, 'copy_event_terms_including_null': copy_events, 'same_site_cross_colour_event_terms': same_site_cross, 'current_cross_product_zero': True, 'population_denominator_coefficients': len(pd.coeffs()), 'founder_denominator_coefficients': len(hd.coeffs())}
    t = gens[0]
    u = s.Symbol('u')
    pdshift = s.Poly(pd.as_expr().subs(t, u + 1), u)
    hdshift = s.Poly(hd.as_expr().subs(t, u + 1), u)
    require(all((x > 0 for x in pdshift.all_coeffs())), 'Independent polynomial/state condition 151 failed')
    require(all((x < 0 for x in hdshift.all_coeffs())), 'Independent polynomial/state condition 152 failed')
    require(all((x > 0 for x in D.all_coeffs())), 'Independent polynomial/state condition 153 failed')
    phase = load(dossier / 'evidence/four-particle/phase-polynomial.json')
    P = poly(phase['polynomial'])
    require(P.all_coeffs() == list(map(s.Integer, phase['descending_coefficients'])), 'Independent polynomial/state condition 156 failed')
    require(A == poly(2 * t * (t - 1) * (3 * t * t + 13 * t + 8)) * P, 'Independent polynomial/state condition 157 failed')
    signs = [s.sign(x) for x in P.all_coeffs() if x]
    require(sum((x != y for x, y in zip(signs, signs[1:]))) == 1, 'Independent polynomial/state condition 159 failed')
    lo = s.Rational(4769, 446)
    hi = s.Rational(3133, 293)
    require(1 < lo < hi and P.eval(1) < 0 and (P.eval(lo) < 0 < P.eval(hi)), 'Independent polynomial/state condition 161 failed')
    require(P.LC() > 0 and P.TC() < 0 and (s.gcd(P, P.diff()).degree() == 0), 'Independent polynomial/state condition 162 failed')
    vals = {str(v): str(s.cancel(A.eval(v) / D.eval(v))) for v in [2, 100]}
    require(s.Rational(vals['2']) < 0 < s.Rational(vals['100']), 'Independent polynomial/state condition 164 failed')
    result.update({'path_denominators_nonzero': True, 'D_all_coefficients_positive': True, 'phase_degree': P.degree(), 'Descartes_variations': 1, 'phase_simple': True, 'root_interval': [str(lo), str(hi)], 'root_is_above1': True, 'exact_points': vals, 'crossing_derivative_positive': True})
    return result
from archive_tools import run
if __name__ == '__main__':
    raise SystemExit(run(build_certificate, 'evidence/four-particle/independent-certificate.json', 'N4 full ninety-six-state two-colour polynomial Dirichlet residuals'))
