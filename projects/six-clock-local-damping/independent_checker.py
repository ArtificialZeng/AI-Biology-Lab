"""Fresh exact verifier: start with stoichiometric reactions, not Fourier RHS.

Public packet adaptation of the fresh reaction-based checker.
Run through verify.py for read-only verification in temporary child directories.
No owner module is imported or executed. Source polynomials are compared only
after the physical-vector-field reconstruction has finished.
"""
from pathlib import Path
import hashlib
import json
import platform
import sys
import time
import sympy as sp

if not __debug__:
    raise RuntimeError("Scientific assertions require enabled checks; use verify.py")

START = time.process_time()
HERE = Path(__file__).resolve().parent
CANDIDATE = HERE


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def zero(expr):
    return sp.expand(sp.expand_complex(expr)) == 0

# Rebuild alpha*T_i -> alpha*(T_i+T_(i+1)) production and
# kappa*T_i*T_(i+2) annihilation as six distinct unordered channels.
T = sp.symbols("T0:6", real=True)
alpha, kappa = sp.symbols("alpha kappa", positive=True)
physical = [sp.S.Zero for _ in range(6)]
pairs = set()
for j in range(6):
    physical[(j + 1) % 6] += alpha * T[j]
    pair = tuple(sorted((j, (j + 2) % 6)))
    assert pair not in pairs
    pairs.add(pair)
    flux = kappa * T[j] * T[(j + 2) % 6]
    for index in pair:
        physical[index] -= flux
assert len(pairs) == 6
for j in range(6):
    assert sp.expand(physical[j] - alpha*T[(j-1)%6]
                     + kappa*T[j]*(T[(j-2)%6]+T[(j+2)%6])) == 0

Z = sp.symbols("Z0:6", real=True)
scaled = [sp.expand(f.subs({T[j]:alpha*Z[j]/(2*kappa)
                          for j in range(6)})*2*kappa/alpha**2)
          for f in physical]
assert all(f.subs(dict.fromkeys(Z, 1)) == 0 for f in scaled)
J = sp.Matrix(scaled).jacobian(Z).subs(dict.fromkeys(Z, 1))
lam = sp.Symbol("lambda")
characteristic = sp.factor(J.charpoly(lam).as_expr())
expected_char = (lam+1)*(lam+3)*(4*lam**2+3)*(4*lam**2+8*lam+7)/16
assert sp.expand(characteristic-expected_char) == 0

# Direct real trigonometric reconstruction followed by projection of the
# physical ODE. The owner instead constructs the Fourier convolution first.
x, y, u, v, p, q = sp.symbols("ar ai u v br bi", real=True)
variables = (x, y, u, v, p, q)
rho = x*x+y*y
A = x+sp.I*y
B = p+sp.I*q
rt = sp.sqrt(3)
H = (2+sp.I*rt)/7
g = -sp.Rational(3,14)+sp.I*rt/7
root = (1+sp.I*rt)/2

def realexpand(expr):
    return sp.expand(sp.re(sp.expand_complex(expr)))

w = [sp.expand(u+rho+(-1)**j*v
               +2*realexpand(A*root**j)
               +2*realexpand((H*A*A+B)*root**(2*j))) for j in range(6)]
f = [sp.expand(F.subs({Z[j]:1+w[j] for j in range(6)})) for F in scaled]
projection = [sp.expand(sum(f[j]*sp.expand_complex(root**(-k*j))
                           for j in range(6))/6) for k in range(4)]
ap = projection[1]
linearA = -sp.I*rt*A/2
center_complex = sp.expand(ap-linearA-g*A*rho)
newY = [sp.expand(projection[0]-2*x*sp.re(ap)-2*y*sp.im(ap)),
        projection[3], sp.expand(projection[2]-2*H*A*ap)]
stable_remainder = [sp.expand(newY[0]+u),
                    sp.expand(newY[1]+3*v),
                    sp.expand(newY[2]-(-1-sp.I*rt/2)*B)]
own_center = [realexpand(center_complex),
              sp.expand(sp.im(sp.expand_complex(center_complex)))]
own_stable = [realexpand(stable_remainder[0]),
              realexpand(stable_remainder[1]),
              realexpand(stable_remainder[2]),
              sp.expand(sp.im(sp.expand_complex(stable_remainder[2])))]
claimed_R1 = -u*A/2+B*sp.conjugate(A)-v*(sp.conjugate(H)*sp.conjugate(A)**2+sp.conjugate(B))/2
assert zero(center_complex-claimed_R1)
assert zero(sp.Rational(1,2)/(2*(-sp.I*rt/2)-(-1-sp.I*rt/2))-H)
assert zero(H-sp.Rational(1,2)-g)

# Independent monomial coverage and rational coefficient majorants.
totals = []
degrees = []
for index, expr in enumerate(own_center+own_stable):
    tally = sp.S.Zero
    degree_counts = {}
    for exponent, coefficient in sp.Poly(expr, *variables).terms():
        a, b = sum(exponent[:2]), sum(exponent[2:])
        if index < 2:
            assert (a >= 1 and b >= 1) or b >= 2
        else:
            assert (b == 0 and a >= 3) or (a >= 1 and b >= 1) or b >= 2
        coefficient = sp.expand(coefficient)
        rational = coefficient.coeff(rt, 0)
        radical = coefficient.coeff(rt, 1)
        assert rational.is_Rational and radical.is_Rational
        assert sp.expand(coefficient-rational-rt*radical) == 0
        tally += abs(rational)+2*abs(radical)
        key = f"{a},{b}"
        degree_counts[key] = degree_counts.get(key, 0)+1
    totals.append(tally)
    degrees.append(degree_counts)
total_majorant = sum(totals)
C = sp.Integer(61)
K = sp.Integer(34730)
delta = sp.Rational(1,100000)
assert total_majorant <= C
assert K >= 28*C*C/3
assert delta <= sp.Rational(1,16)
assert delta <= K/(8*C*(1+2*K))
assert delta*delta <= 3/(56*K*C*C)

# Check the two Young inequalities as exact completed squares and the
# remaining scalar-budget deductions, not floating-point samples.
r, s = sp.symbols("r s", nonnegative=True)
assert sp.expand(sp.Rational(3,28)*r**4+28*C*C*s*s/3-2*C*r*r*s
                 -(sp.sqrt(sp.Rational(3,28))*r*r-sp.sqrt(28*C*C/3)*s)**2) == 0
assert sp.expand(K*s*s/2+2*K*C*C*r**6-2*K*C*r**3*s
                 -(sp.sqrt(K/2)*s-sp.sqrt(2*K*C*C)*r**3)**2) == 0
assert 2*C*delta*(1+2*K) <= K/4
assert 2*K*C*C*delta*delta <= sp.Rational(3,28)
assert -2*K+28*C*C/3+K/2+K/4 <= -K/4
assert -sp.Rational(3,7)+sp.Rational(3,28)+sp.Rational(3,28) == -sp.Rational(3,14)
# |H|<=1/2; stable inverse DFT contribution <=sqrt(6)||Y||<=3||Y||.
# Thus |w_i| <= 2r+3s+2r^2 <= 5delta+2delta^2 <= 6delta+delta^2 <1.
assert sp.expand(H*sp.conjugate(H)) == sp.Rational(1,7)
assert sp.Rational(1,7) <= sp.Rational(1,4)
assert 5*delta+2*delta*delta <= 6*delta+delta*delta < 1

# A comparison calculation for the inspected 2017 signed-cosine predecessor.
# This is our deduction from its displayed ODE, not an attributed new theorem.
fp, fm, zp, zm, fast = sp.symbols("fp fm zp zm fast", real=True)
prior_rhs = sp.Matrix([zp-fast*fm*fp, zm-fast*fm*fp,
                       fm-fast*zm*zp, fp-fast*zm*zp])
prior_energy = (fp-fm)**2+(zp-zm)**2
assert sp.expand((sp.Matrix([prior_energy]).jacobian((fp,fm,zp,zm))*prior_rhs)[0]) == 0

# Compare only now, after reconstructing all decisive identities.
record = json.loads((CANDIDATE/"local-algebra.json").read_text())
loc = {str(z):z for z in variables}
loc.update({"sqrt":sp.sqrt,"I":sp.I})
for actual, saved in zip(own_center, record["real_center_remainders"], strict=True):
    assert sp.expand(actual-sp.sympify(saved, locals=loc)) == 0
for actual, saved in zip(own_stable, record["real_stable_remainders"], strict=True):
    assert sp.expand(actual-sp.sympify(saved, locals=loc)) == 0
assert zero(g-sp.sympify(record["center_cubic_coefficient_g"], locals=loc))


out = {
    "schema":"independent-six-clock-exact-check-v1",
    "all_checks_passed":True,
    "method":"Stoichiometric physical ODE -> real inverse DFT -> physical RHS projection -> polynomial comparison; owner script not imported or executed",
    "sympy_version":sp.__version__,
    "physical_reaction_pairs":[list(z) for z in sorted(pairs)],
    "characteristic_polynomial":str(characteristic),
    "H":str(H),"g":str(g),
    "component_majorants":[str(z) for z in totals],
    "sum_component_majorants":str(total_majorant),
    "monomial_degrees_center_stable":degrees,
    "C":int(C),"K":int(K),"delta":str(delta),
    "strict_derivative_bound":"Vprime <= -3*r^4/14-K*s^2/4",
    "positivity_bound":"|w_i|<=2r+3s+2r^2<=5delta+2delta^2<=6delta+delta^2<1",
    "frozen_polynomial_match":True,
    "prior_four_species_signed_energy_derivative_zero":True,
    "analytic_items_not_machine_formalized":["local existence and continuation for polynomial ODE", "compact positively invariant smaller sublevels", "strict Lyapunov/omega-limit convergence", "contribution judgment"],
    "global_H6_status":"unresolved",
    "formal_math_acceptance":False,
}
(HERE/"independent-certificate.json").write_text(json.dumps(out, indent=2, sort_keys=True)+"\n")
print(json.dumps({"all_checks_passed":True,"C":int(C),"K":int(K),
                  "sum_component_majorants":str(total_majorant),
                  "python":platform.python_version(),"sympy":sp.__version__,
                  "scope":"local theorem only"}, sort_keys=True))
