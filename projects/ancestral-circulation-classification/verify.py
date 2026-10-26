"""Read-only exact certificate verification; no numeric linear algebra."""
import flint
import json
from pathlib import Path
import sympy as s
from flint import fmpq, fmpq_mpoly_ctx
import model
HERE = Path(__file__).resolve().parent
A, B, C = s.symbols('a b c')
GENS = (A, B, C)
CTX = fmpq_mpoly_ctx.get(['a', 'b', 'c'])

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def poly(expr):
    q = s.Poly(expr, *GENS, domain=s.QQ)
    return CTX.from_dict({mon: fmpq(str(co)) for mon, co in q.terms()})

def load_solution(label):
    j = json.loads((HERE / 'evidence' / (label + '-solution.json')).read_text())
    return ([poly(s.sympify(v)) for v in j['numerators']], poly(s.sympify(j['denominator'])))

def matvec(M, v):
    return [sum((poly(M[i, j]) * v[j] for j in range(M.cols)), CTX.constant(0)) for i in range(M.rows)]

def require_zero(values, message):
    require(all((v == 0 for v in values)), message)

def build_certificate():
    require(s.__version__ == '1.13.3', 'SymPy 1.13.3 required')
    require(flint.__version__ == '0.9.0', 'python-flint 0.9.0 required')
    pop, states, ci, si, L, T, P, pr, H, hr = model.normalized_systems(GENS)
    pn, pd = load_solution('population')
    hn, hd = load_solution('founder')
    require_zero(matvec(L.T, pn), 'Population stationarity residual')
    require(sum(pn, CTX.constant(0)) == pd, 'Population normalization')
    require_zero(matvec(T, hn), 'Founder equation residual')
    require(3 * hn[si[0, (2, 0, 0)]] == hd, 'Founder anchor')
    for c in pop:
        check = sum((c[i] * hn[si[i, tuple((c[j] - int(i == j) for j in range(3)))]] for i in range(3) if c[i]), CTX.constant(0))
        require(check == hd, 'Founder conservation normalization')
    require(all((co > 0 for co in pd.to_dict().values())), 'Population denominator positivity')
    require(all((co < 0 for co in hd.to_dict().values())), 'Founder denominator negative positivity certificate')
    data = json.loads((HERE / 'evidence/symbolic-current.json').read_text())
    num = s.sympify(data['numerator'])
    den = s.sympify(data['denominator'])
    nr = CTX.constant(0)
    for c in pop:
        if c[0]:
            d = list(c)
            d[0] -= 1
            nr += pn[ci[c]] * c[0] * hn[si[2, tuple(d)]]
        if c[2]:
            d = list(c)
            d[2] -= 1
            nr -= pn[ci[c]] * c[2] * hn[si[0, tuple(d)]]
    require(nr * poly(den) == poly(num) * pd * hd, 'Exact current cross-product identity')
    nf = s.factor_list(num)
    df = s.factor_list(den)
    require(nf[0] == -1, 'Numerator global sign')
    Q = max((f for f, m in nf[1]), key=s.total_degree)
    V = (A - B) * (A - C) * (B - C)
    S = A + B + C + 9
    K = A * B * C + 2 * (A * B + A * C + B * C) + 3 * (A + B + C)
    require(poly(num) == -poly(V * S * K * Q), 'Displayed numerator factorization')
    require(len(df[1]) == 2 and all((m == 1 for f, m in df[1])), 'Two denominator factors')
    require(df[0] > 0, 'Current denominator overall scalar positivity')
    denominator_counts = []
    for f, m in df[1]:
        cs = s.Poly(f, *GENS).coeffs()
        require(all((co > 0 for co in cs)), 'Current denominator coefficient positivity')
        denominator_counts.append(len(cs))
    for swap in [{A: B, B: A}, {B: C, C: B}]:
        require(poly(Q.subs(swap, simultaneous=True)) == poly(Q), 'Q symmetry')
    t, x, y = s.symbols('t x y')
    ordered = s.Poly(s.expand(Q.subs({A: t, B: t + x, C: t + x + y}, simultaneous=True)), t, x, y)
    require(all((co > 0 for co in ordered.coeffs())), 'Ordered Q positive coefficients')
    groups = s.Poly(ordered.as_expr(), x, y)
    require(len(groups.terms()) == 30 and len(ordered.terms()) == 94, 'Ordered Q coverage')
    minimum = 27 * (t + 2) * (t + 3) ** 2 * (3 * t * t + 26 * t + 60)
    for mon in [x * x, x * y, y * y]:
        require(s.expand(groups.coeff_monomial(mon) - minimum) == 0, 'Strict off-diagonal positivity term')
    require(Q.subs({A: t, B: t, C: t}) == 0, 'Q vanishes on full diagonal')
    eps = s.Symbol('eps')
    substituted_n = s.Poly(s.expand(num.subs({A: 1 - eps, B: 1, C: 1 + eps}, simultaneous=True)), eps)
    substituted_d = s.Poly(s.expand(den.subs({A: 1 - eps, B: 1, C: 1 + eps}, simultaneous=True)), eps)
    expected_n = s.Poly(-24 * eps ** 5 * (3 * eps ** 2 - 16) * (479 * eps ** 2 + 38448), eps)
    expected_d = s.Poly((5 * eps ** 6 + 279 * eps ** 4 + 4704 * eps ** 2 + 790272) * (143 * eps ** 6 + 840 * eps ** 4 - 137808 * eps ** 2 + 11757312), eps)
    require((substituted_n * expected_d - expected_n * substituted_d).is_zero, 'Exact perturbation-slice identity')
    require(all((substituted_n.coeff_monomial(eps ** j) == 0 for j in range(5))), 'No lower perturbation term')
    require(substituted_n.coeff_monomial(eps ** 5) / substituted_d.coeff_monomial(1) == s.Rational(89, 56010528), 'Exact quintic coefficient')
    value = s.cancel(num.subs({A: 1, B: 2, C: 4}) / den.subs({A: 1, B: 2, C: 4}))
    require(value == s.Rational(421413160, 139117744582101), 'Published rational source witness')
    e = s.symbols('e')
    substitution = {A: 1 - e, B: 1, C: 1 + e}
    slice_num = s.Poly(s.expand(num.subs(substitution)), e)
    slice_den = s.Poly(s.expand(den.subs(substitution)), e)
    first_power = next((j for j in range(slice_num.degree() + 1) if slice_num.nth(j)))
    quintic = slice_num.nth(first_power) / slice_den.nth(0)
    require(first_power == 5 and quintic == s.Rational(89, 56010528), 'Exact N3 quintic response')
    displayed_numerator = -24 * e ** 5 * (3 * e ** 2 - 16) * (479 * e ** 2 + 38448)
    displayed_denominator = (5 * e ** 6 + 279 * e ** 4 + 4704 * e ** 2 + 790272) * (143 * e ** 6 + 840 * e ** 4 - 137808 * e ** 2 + 11757312)
    require(s.expand(slice_num.as_expr() * displayed_denominator - displayed_numerator * slice_den.as_expr()) == 0, 'Exact displayed symmetric-perturbation identity')
    result = {'ok': True, 'kind': 'exact polynomial identities plus positive-coefficient certificate', 'population_states': 10, 'founder_unknowns': 18, 'all_positive_domain': True, 'ordered_Q_coefficient_groups': 30, 'ordered_Q_positive_scalar_coefficients': 94, 'denominator_factor_terms': sorted(denominator_counts), 'current_at_source_point': str(value), 'zero_locus': 'k0=k1 or k0=k2 or k1=k2', 'sign': 'sign(-(k0-k1)*(k0-k2)*(k1-k2))', 'archived_evidence_modified': False, 'proof_assistant_formalization': False, 'N3_symmetric_perturbation_first_power': first_power, 'N3_symmetric_perturbation_coefficient': str(quintic)}
    return result
from archive_tools import run
if __name__ == '__main__':
    raise SystemExit(run(build_certificate, 'evidence/exact-certificate.json', 'ten-state and eighteen-founder polynomial residuals'))
