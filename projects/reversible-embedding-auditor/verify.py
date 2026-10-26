"""Original exact construction, refactored for read-only certificate replay."""
import sys
sys.dont_write_bytecode = True
import sympy as s
from archive_tools import run

def build_certificate():
    if s.__version__ != '1.13.3':
        raise RuntimeError('Use pinned SymPy1.13.3')
    U = s.Matrix([[3, -3, -3, -3], [3, -1, 5, -1], [3, 5, -1, -1], [3, -1, -1, 5]]) / 6
    if not U.T * U == s.eye(4):
        raise RuntimeError('Exact identity failed (original line 16)')
    e, v, w, u = [U[:, i] for i in range(4)]
    E = e * e.T
    P = v * v.T + w * w.T
    R = u * u.T
    J = 2 * v * w.T - s.Rational(1, 2) * w * v.T
    A = -R - 3 * P
    B = 2 * J
    if not E + P + R == s.eye(4):
        raise RuntimeError('Exact identity failed (original line 24)')
    if not all((X * X == X for X in (E, P, R))):
        raise RuntimeError('Exact identity failed (original line 25)')
    if not all((X * Y == s.zeros(4) for X, Y in [(E, P), (P, E), (E, R), (R, E), (P, R), (R, P)])):
        raise RuntimeError('Exact identity failed (original line 26)')
    if not J * J == -P:
        raise RuntimeError('Exact identity failed (original line 28)')
    if not J * E == E * J == J * R == R * J == s.zeros(4):
        raise RuntimeError('Exact identity failed (original line 29)')
    if not J * P == P * J == J:
        raise RuntimeError('Exact identity failed (original line 30)')
    if not A * J == J * A:
        raise RuntimeError('Exact identity failed (original line 31)')
    ones = s.ones(4, 1)
    if not A * ones == B * ones == s.zeros(4, 1):
        raise RuntimeError('Exact identity failed (original line 33)')
    if not ones.T * A == ones.T * B == s.zeros(1, 4):
        raise RuntimeError('Exact identity failed (original line 34)')
    h, t = s.symbols('h t', positive=True)
    Q_form = h * A + t * B
    if not U.T * Q_form * U == s.Matrix([[0, 0, 0, 0], [0, -3 * h, 4 * t, 0], [0, -t, -3 * h, 0], [0, 0, 0, -h]]):
        raise RuntimeError('Exact identity failed (original line 39)')
    lower_bounds = {}
    for i in range(4):
        for j in range(4):
            if i == j:
                continue
            expr = s.expand(Q_form[i, j])
            ch, ct = (expr.coeff(h), expr.coeff(t))
            lb = ch * (6 if ch >= 0 else 9) + ct * (3 if ct >= 0 else s.Rational(22, 7))
            if not (lb >= 0 and (ch != 0 or ct != 0)):
                raise RuntimeError('Exact identity failed (original line 55)')
            lower_bounds[f'{i + 1},{j + 1}'] = str(lb)
    M = E + s.Rational(1, 512) * R + s.Rational(1, 512 ** 3) * P
    den = s.ilcm(*[x.q for x in M])
    if not den == 268435456:
        raise RuntimeError('Exact identity failed (original line 60)')
    if not (M == M.T and M * ones == ones):
        raise RuntimeError('Exact identity failed (original line 61)')
    if not all((x > 0 for x in M)):
        raise RuntimeError('Exact identity failed (original line 62)')
    if not U.T * M * U == s.diag(1, s.Rational(1, 512 ** 3), s.Rational(1, 512 ** 3), s.Rational(1, 512)):
        raise RuntimeError('Exact identity failed (original line 63)')
    if not A[0, 3] == -s.Rational(1, 12):
        raise RuntimeError('Exact identity failed (original line 65)')
    receipt = {'status': 'exact-identities-passed', 'scope': 'one rational symmetric 4-state Markov matrix; Thm3.16 (1)=>(3)', 'M_denominator': int(den), 'M_numerator': [[int(x) for x in row] for row in (den * M).tolist()], 'eigenvalues_M': ['1', '2^-9', '2^-27', '2^-27'], 'h': '9*log(2)', 'Q': 'h*A + pi*B', '36A': [[str(x) for x in row] for row in (36 * A).tolist()], '12B': [[str(x) for x in row] for row in (12 * B).tolist()], 'principal_log': 'h*A', 'principal_log_14': '-3*log(2)/4 < 0', 'Q_offdiagonal_strictly_positive': True, 'Q_offdiagonal_interval_infima': lower_bounds, 'analytic_exponential_identity': 'exp(Q)=E+exp(-h)*R+exp(-3*h)*P=M; exp(2*pi*J)=I', 'analytic_inputs': ['log(2)>2/3', 'log(2)<1', '3<pi<22/7'], 'floating_point_used_for_verification': False}
    return receipt

if __name__ == '__main__':
    raise SystemExit(run(build_certificate, 'evidence/exact-certificate.json', 'orthogonal-projector construction'))
