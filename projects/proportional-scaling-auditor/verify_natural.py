"""Exact mortality-ordered witnesses from the population equations, without floats.

The rational dual calculation is a port of the frozen exact candidate calculation.
It independently checks the explicit quotient Jacobian. Verification is read-only;
archive_tools controls all optional external regeneration.
"""
from fractions import Fraction as Q
import sys
sys.dont_write_bytecode = True
from archive_tools import run, require


class Dual:
    def __init__(self, value, gradient=(0, 0)):
        self.value = Q(value)
        self.gradient = tuple(Q(x) for x in gradient)

    def __add__(self, other):
        other = as_dual(other)
        return Dual(self.value + other.value,
                    tuple(a + b for a, b in zip(self.gradient, other.gradient)))

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.value, tuple(-x for x in self.gradient))

    def __sub__(self, other):
        return self + -as_dual(other)

    def __rsub__(self, other):
        return as_dual(other) + -self

    def __mul__(self, other):
        other = as_dual(other)
        return Dual(self.value * other.value,
                    tuple(a * other.value + self.value * b
                          for a, b in zip(self.gradient, other.gradient)))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = as_dual(other)
        require(other.value != 0, 'Division by a zero dual value')
        return Dual(self.value / other.value,
                    tuple((a * other.value - self.value * b) / other.value**2
                          for a, b in zip(self.gradient, other.gradient)))

    def __rtruediv__(self, other):
        return as_dual(other) / self


def as_dual(value):
    return value if isinstance(value, Dual) else Dual(value)


def matrix(mu, transfers):
    return [[mu[i] + sum(transfers[k][i] for k in range(2)) if i == j
             else -transfers[i][j] for j in range(2)] for i in range(2)]


def quotient(nu, mu_s, mu, transfers, entry, beta, y):
    """Differentiate the original number equations along the N=1 fraction chart."""
    x = [Dual(y[0], (1, 0)), Dual(y[1], (0, 1))]
    susceptible = 1 - x[0] - x[1]
    total = susceptible + sum(x)
    force = sum(beta[i] * x[i] for i in range(2)) / total
    pit = matrix(mu, transfers)
    raw_s = nu * total - mu_s * susceptible - susceptible * force
    raw_x = [susceptible * force * entry[i]
             - sum(pit[i][j] * x[j] for j in range(2)) for i in range(2)]
    raw_total = raw_s + sum(raw_x)
    reduced = [(raw_x[i] * total - x[i] * raw_total) / (total * total)
               for i in range(2)]
    return reduced, raw_s, raw_x, raw_total, pit


def signs(jacobian):
    trace = jacobian[0][0] + jacobian[1][1]
    determinant = (jacobian[0][0] * jacobian[1][1]
                   - jacobian[0][1] * jacobian[1][0])
    return trace, determinant


def strings(values):
    return [str(x) for x in values]


def matrix_strings(values):
    return [strings(row) for row in values]


def point(name, nu, mu_s, mu, transfers, entry, beta, s, y, strict=False):
    require(nu > 0 and mu_s > 0, name + ': positive vital rates')
    require(all(m >= mu_s for m in mu), name + ': ordered mortalities')
    require(all(b >= 0 for b in beta) and any(b > 0 for b in beta),
            name + ': nonnegative, nonzero infectivity')
    require(all(p > 0 for p in entry) and sum(entry) == 1,
            name + ': entry probabilities directly reach both nodes')
    require(s > 0 and all(v > 0 for v in y) and s + sum(y) == 1,
            name + ': strictly positive normalized fractions')
    require(all(v >= 0 for row in transfers for v in row)
            and transfers[0][0] == transfers[1][1] == 0,
            name + ': source transfer convention')
    if strict:
        require(all(m > mu_s for m in mu) and all(b > 0 for b in beta)
                and transfers[0][1] > 0 and transfers[1][0] > 0,
                name + ': strict mortality, infectivity and transfers')
    reduced, raw_s, raw_x, raw_total, pit = quotient(
        nu, mu_s, mu, transfers, entry, beta, y)
    require(all(z.value == 0 for z in reduced), name + ': fraction stationarity')
    growth = raw_total.value
    require(raw_s.value == growth * s
            and all(raw_x[i].value == growth * y[i] for i in range(2)),
            name + ': number solution has the stated balanced growth')
    require(growth > 0, name + ': positive population growth')
    require(all(sum(pit[i][j] for i in range(2)) == mu[j] for j in range(2)),
            name + ': PIT mortality column sums')
    jacobian = [list(z.gradient) for z in reduced]
    delta = [m - mu_s for m in mu]
    force_value = sum(beta[i] * y[i] for i in range(2))
    explicit = [[entry[i] * (s * beta[j] - force_value) - pit[i][j]
                 - (growth if i == j else 0) + y[i] * delta[j]
                 for j in range(2)] for i in range(2)]
    require(jacobian == explicit, name + ': dual and explicit Jacobians agree')
    trace, determinant = signs(jacobian)
    return {'name': name, 'nu': str(nu), 'muS': str(mu_s), 'mu': strings(mu),
            'transfers': matrix_strings(transfers), 'M': matrix_strings(pit),
            'P': strings(entry), 'beta': strings(beta), 's': str(s), 'y': strings(y),
            'fraction_rhs': strings(z.value for z in reduced),
            'raw_numbers_balanced_growth_residual': ['0', '0', '0'],
            'g_numbers_growth': str(growth), 'Jacobian': matrix_strings(jacobian),
            'trace': str(trace), 'determinant': str(determinant),
            'characteristic_discriminant': str(trace**2 - 4 * determinant),
            'original_numbers_dual_agrees_explicit_Jacobian': True,
            'mortality_order_verified': True, 'all_nodes_directly_reachable': True,
            'no_return_to_S': True, 'no_infected_births': True,
            'same_per_capita_recruitment_into_S': True, 'all_classes_in_N': True,
            'strict_excess_beta_and_transfer_requested': strict}, trace, determinant


def disease_free(nu, mu_s, mu, transfers, entry, beta):
    reduced, raw_s, raw_x, raw_total, _ = quotient(
        nu, mu_s, mu, transfers, entry, beta, [Q(0), Q(0)])
    require(all(z.value == 0 for z in reduced), 'DFE fraction stationarity')
    require(raw_s.value == nu - mu_s and all(z.value == 0 for z in raw_x),
            'DFE number equation and growth')
    jacobian = [list(z.gradient) for z in reduced]
    trace, determinant = signs(jacobian)
    require(trace < 0 and determinant > 0, 'DFE is a hyperbolic sink')
    return {'Jacobian': matrix_strings(jacobian), 'trace': str(trace),
            'determinant': str(determinant), 'stable': True,
            'g_numbers_growth': str(raw_total.value)}


def build_certificate():
    nu, mu_s = Q(1), Q(1, 10)
    mu = [Q(1, 10), Q(91, 10)]
    transfers = [[Q(0), Q(0)], [Q(0), Q(0)]]
    entry, beta = [Q(1, 10), Q(9, 10)], [Q(135, 31), Q(0)]
    saddle, trace, determinant = point('simple_positive_growth_saddle', nu, mu_s,
        mu, transfers, entry, beta, Q(31, 54), [Q(37, 108), Q(1, 12)])
    require(trace == -Q(1239, 124) and determinant == -Q(37, 31),
            'Simple saddle exact trace and determinant')
    sink, trace, determinant = point('same_model_second_positive_sink', nu, mu_s,
        mu, transfers, entry, beta, Q(31, 66), [Q(175, 396), Q(35, 396)])
    require(trace < 0 and determinant == Q(525, 341), 'Positive hyperbolic sink')
    dfe = disease_free(nu, mu_s, mu, transfers, entry, beta)
    mu = [Q(11, 100), Q(911, 100)]
    transfers = [[Q(0), Q(1, 1000)], [Q(1, 1000), Q(0)]]
    entry = [Q(256241, 2504630), Q(2248389, 2504630)]
    beta = [Q(5008981, 1147000), Q(1, 1000)]
    strict, trace, determinant = point('strict_positive_growth_saddle', nu, mu_s,
        mu, transfers, entry, beta, Q(31, 54), [Q(37, 108), Q(1, 12)], strict=True)
    require(strict['g_numbers_growth'] == '787/5400', 'Strict growth value')
    require(trace == -Q(387768093698339, 38782943235000)
            and determinant == -Q(1155754208768357, 969573580875000),
            'Strict saddle exact trace and determinant')
    strict['DFE'] = disease_free(nu, mu_s, mu, transfers, entry, beta)
    strict['positive_sink_claimed'] = False
    return {'schema': 'natural-mortality-exact-certificate-v1',
            'arithmetic': 'fractions.Fraction rational dual differentiation of original population equations',
            'simple_saddle': saddle, 'same_model_positive_sink': sink,
            'same_model_DFE': dfe, 'strict_saddle': strict,
            'scope': 'Exact witnesses refute universal fraction GAS under ordered mortality. '
                     'Two sinks are certified only for the simple model. No oscillations, '
                     'new biological mechanism, priority or external human validation is claimed.'}


if __name__ == '__main__':
    raise SystemExit(run(build_certificate, 'evidence/natural-mortality-certificate.json',
                        'ordered-mortality original-number rational dual calculation'))
