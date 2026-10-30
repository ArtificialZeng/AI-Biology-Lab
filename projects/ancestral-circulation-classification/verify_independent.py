"""Independent full two-colour polynomial checker; shared candidate data, separate event construction."""
import hashlib
import json
from collections import defaultdict
from pathlib import Path
import sympy as sp
from flint import fmpq, fmpq_mpoly_ctx
import flint

def build_certificate():
    P = Path(__file__).resolve().parent
    a, b, c = sp.symbols('a b c')
    gens = (a, b, c)
    ctx = fmpq_mpoly_ctx.get(['a', 'b', 'c'])
    Z = ctx.constant(0)
    O = ctx.constant(1)
    k = ctx.gens()

    def require(ok, message):
        if not ok:
            raise RuntimeError(message)
    require(sp.__version__ == '1.13.3', 'SymPy 1.13.3 required')
    require(flint.__version__ == '0.9.0', 'python-flint 0.9.0 required')

    def cp(n):
        return tuple(((x, y, n - x - y) for x in range(n + 1) for y in range(n - x + 1)))

    def pol(expr):
        q = sp.Poly(sp.sympify(expr), *gens, domain=sp.QQ)
        return ctx.from_dict({mon: fmpq(str(co)) for mon, co in q.terms()})

    def load(label):
        dat = json.loads((P / 'evidence' / (label + '-solution.json')).read_text())
        return (list(map(pol, dat['numerators'])), pol(dat['denominator']))
    pn, pd = load('population')
    hn, hd = load('founder')
    pop = cp(3)
    founders = tuple(((i, d) for i in range(3) for d in cp(2)))
    pi = dict(zip(pop, pn))
    h = dict(zip(founders, hn))
    order = json.loads((P / 'evidence/state-ordering.json').read_text())
    require(tuple(map(tuple, order['population'])) == pop, 'population labels')
    require(tuple(((v['i'], tuple(v['d'])) for v in order['founder'])) == founders, 'founder labels')
    require(len(pn) == 10 and len(hn) == 18, 'serialized sizes')

    def moved(state, origin, target):
        out = list(state)
        out[origin] -= 1
        out[target] += 1
        return tuple(out)

    def events(state, colours):
        """Every ordered victim/donor pair is counted with donor denominator N-1=2."""
        jumps = defaultdict(lambda: ctx.constant(0))
        for victim, number in enumerate(state):
            if not number:
                continue
            colour, site = divmod(victim, 3)
            for dest in range(3):
                if dest != site:
                    jumps[moved(state, victim, 3 * colour + dest)] += number
            for donor, donor_number in enumerate(state):
                available = donor_number - int(victim == donor)
                if available > 0 and donor != victim:
                    jumps[moved(state, victim, donor)] += k[site] * fmpq(number * available, 2)
        require(all((min(s) >= 0 and sum(s) == 3 for s in jumps)), 'valid event targets')
        return jumps
    residuals = {s: ctx.constant(0) for s in pop}
    population_edges = 0
    for state in pop:
        jumps = events(state, 1)
        population_edges += len(jumps)
        for target, rate in jumps.items():
            residuals[target] += pi[state] * rate
            residuals[state] -= pi[state] * rate
    require(all((v == 0 for v in residuals.values())), 'independent population stationarity')
    require(sum(pn, Z) == pd, 'population normalization')
    require(all((co > 0 for co in pd.to_dict().values())), 'population denominator positive')
    require(all((co < 0 for co in hd.to_dict().values())), 'founder denominator negative')
    normalizations = 0
    for state in pop:
        value = Z
        for site, count in enumerate(state):
            if count:
                other = list(state)
                other[site] -= 1
                value += count * h[site, tuple(other)]
        require(value == hd, 'all-founder normalization')
        normalizations += 1
    mixed = tuple((red + blue for nred in (1, 2) for red in cp(nred) for blue in cp(3 - nred)))
    require(len(mixed) == 36, 'mixed two-colour states')

    def red_win_numerator(state):
        nr = sum(state[:3])
        if nr == 0:
            return Z
        if nr == 3:
            return hd
        total = tuple((state[j] + state[3 + j] for j in range(3)))
        ans = Z
        for site, reds in enumerate(state[:3]):
            if reds:
                other = list(total)
                other[site] -= 1
                ans += reds * h[site, tuple(other)]
        return ans
    coloured_edges = 0
    same_site_opposite_colour_edges = 0
    for state in mixed:
        value = red_win_numerator(state)
        residual = Z
        jumps = events(state, 2)
        coloured_edges += len(jumps)
        for target, rate in jumps.items():
            residual += rate * (red_win_numerator(target) - value)
            if tuple((state[j] + state[3 + j] for j in range(3))) == tuple((target[j] + target[3 + j] for j in range(3))):
                same_site_opposite_colour_edges += 1
        require(residual == 0, 'full two-colour Dirichlet polynomial residual')
    require(3 * h[0, (2, 0, 0)] == hd, 'anchor')
    curr_num = Z
    for state in pop:
        for source, target, sign in ((0, 2, 1), (2, 0, -1)):
            if state[source]:
                other = list(state)
                other[source] -= 1
                curr_num += sign * pi[state] * state[source] * h[target, tuple(other)]
    dat = json.loads((P / 'evidence/symbolic-current.json').read_text())
    num, den = map(sp.sympify, (dat['numerator'], dat['denominator']))
    require(curr_num * pol(den) == pol(num) * pd * hd, 'independent current cross-product')
    V = (a - b) * (a - c) * (b - c)
    S = a + b + c + 9
    K = a * b * c + 2 * (a * b + a * c + b * c) + 3 * (a + b + c)
    Q, remainder = sp.div(-num, V * S * K, *gens)
    require(remainder == 0, 'exact Vandermonde quotient')
    require(pol(-V * S * K * Q) == pol(num), 'numerator sign identity')
    factor_scalar, denominator_factors = sp.factor_list(den)
    require(factor_scalar == 1, 'denominator overall scalar is +1')
    require(len(denominator_factors) == 2 and all((e == 1 for _, e in denominator_factors)), 'two simple denominator factors')
    denominator_terms = []
    for factor, exponent in denominator_factors:
        coeffs = sp.Poly(factor, *gens).coeffs()
        require(all((v > 0 for v in coeffs)), 'denominator coefficient positivity')
        denominator_terms.append(len(coeffs))
    for perm in ({a: b, b: a}, {b: c, c: b}):
        require(pol(Q.subs(perm, simultaneous=True)) == pol(Q), 'quotient symmetry')
    t, x, y = sp.symbols('t x y')
    ordered = sp.Poly(sp.expand(Q.subs({a: t, b: t + x, c: t + x + y}, simultaneous=True)), t, x, y)
    require(all((v > 0 for v in ordered.coeffs())), 'whole ordered quotient positive coefficients')
    groups = sp.Poly(ordered.as_expr(), x, y)
    C = 27 * (t + 2) * (t + 3) ** 2 * (3 * t ** 2 + 26 * t + 60)
    require(all((sp.expand(groups.coeff_monomial(m) - C) == 0 for m in (x * x, x * y, y * y))), 'strict positive quadratic bound')
    require(sp.expand(Q.subs({a: t, b: t, c: t})) == 0, 'full diagonal quotient zero')
    return {'ok': True, 'method': 'independent victim/donor event generator; full polynomial residuals', 'candidate_code_imported': False, 'candidate_data_read': ['state-ordering.json', 'population-solution.json', 'founder-solution.json', 'symbolic-current.json'], 'population_states': len(pop), 'population_directed_edges': population_edges, 'population_stationary_polynomial_residuals_zero': len(residuals), 'population_normalization': True, 'founder_normalization_polynomial_identities': normalizations, 'mixed_colour_states': len(mixed), 'mixed_colour_directed_edges': coloured_edges, 'same_site_colour_change_edges_retained': same_site_opposite_colour_edges, 'mixed_colour_Dirichlet_polynomial_residuals_zero': len(mixed), 'domain': 'all real a,b,c; physical conclusion on a,b,c>0', 'population_denominator_positive_terms': len(pd.to_dict()), 'founder_denominator_negative_terms': len(hd.to_dict()), 'current_cross_product_identity': True, 'denominator_factor_overall_scalar': str(factor_scalar), 'denominator_factor_positive_terms': sorted(denominator_terms), 'Q_total_degree': int(sp.total_degree(Q)), 'Q_symmetric': True, 'ordered_Q_coefficient_groups': len(groups.terms()), 'ordered_Q_positive_scalar_coefficients': len(ordered.terms()), 'positive_quadratic_bound': str(C), 'sympy_version': sp.__version__}

from archive_tools import run

if __name__ == "__main__":
    raise SystemExit(run(build_certificate, "evidence/independent-certificate.json", "full two-colour Dirichlet polynomial residuals"))
