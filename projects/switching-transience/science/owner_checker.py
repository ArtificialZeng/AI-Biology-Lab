from pathlib import Path
from fractions import Fraction
import itertools
import json
import math
import sympy as s

ROOT = Path(__file__).resolve().parent
MODEL = json.loads((ROOT / 'model.json').read_text())
alpha, kappa, t = s.symbols('alpha kappa t', positive=True)
x = s.Matrix(s.symbols('x1 x2 x3', nonnegative=True))
identities = {}

def equal(name, expression):
    result = s.simplify(expression)
    if result != 0:
        raise RuntimeError(f'{name} failed: {result}')
    identities[name] = str(result)

def drift(env):
    answer = s.zeros(3, 1)
    for reaction in MODEL['reactions'][str(env)]:
        source = reaction['source']
        if sum(source) not in (0, 1):
            raise RuntimeError('not at-most-unimolecular source')
        propensity = s.sympify(reaction['rate_constant'], locals={'alpha': alpha})
        for j in range(3):
            if source[j]:
                propensity *= x[j]
        answer += propensity * (s.Matrix(reaction['product']) - s.Matrix(source))
    return answer

D1, D2 = drift(1), drift(2)
M1, M2 = D1.jacobian(x), D2.jacobian(x)
b1, b2 = D1.subs(dict.fromkeys(x, 0)), D2.subs(dict.fromkeys(x, 0))
expected = [s.Matrix([[3, 1, 0], [1, -1, 0], [0, 0, -alpha]]),
            s.Matrix([[-alpha, 0, 0], [0, -1, 1], [0, 1, 3]])]
for env, M in enumerate([M1, M2], start=1):
    for i in range(3):
        for j in range(3):
            equal(f'env{env}_mean_entry_{i}_{j}', M[i,j] - expected[env-1][i,j])
for env, b, expected_b in [(1, b1, s.Matrix([0,0,1])), (2, b2, s.Matrix([1,0,0]))]:
    for i in range(3):
        equal(f'env{env}_inflow_{i}', b[i] - expected_b[i])
equal('env1_total_drift', sum(D1) - (4*x[0]-alpha*x[2]+1))
equal('env2_total_drift', sum(D2) - (4*x[2]-alpha*x[0]+1))
A = s.Matrix([[3,1],[1,-1]])
for env, M, indices in [(1,M1,[0,1]),(2,M2,[2,1])]:
    active = M.extract(indices, indices)
    for i in range(2):
        for j in range(2):
            equal(f'env{env}_active_entry_{i}_{j}', active[i,j]-A[i,j])
    for reaction in MODEL['reactions'][str(env)]:
        source = reaction['source']
        if sum(source) == 1 and any(source[j] for j in indices):
            if any(reaction['product'][j] for j in range(3) if j not in indices):
                raise RuntimeError('tracked active particle creates an inactive type')
            increment = sum(reaction['product'][j]-source[j] for j in indices)
            if increment not in (0,4):
                raise RuntimeError('active family can decrease')
r = s.sqrt(5)
lp, lm = 1+r, 1-r
cp, cm = (1-2/r)/2, (1+2/r)/2
vp, vm = s.Matrix([1/(2*r),cp]), s.Matrix([-1/(2*r),cm])
for label, eigenvalue, vector in [('plus',lp,vp),('minus',lm,vm)]:
    for i in range(2):
        equal(f'active_eigenvector_{label}_{i}', (A*vector-eigenvalue*vector)[i])
for i in range(2):
    equal(f'active_initial_condition_{i}', (vp+vm-s.Matrix([0,1]))[i])
equal('positive_coefficients_sum', cp+cm-1)
if not (r > s.Rational(20,9) and cp > s.Rational(1,20) and lp > 3):
    raise RuntimeError('strict algebraic inequalities failed')
exp3_partial = sum(Fraction(3**j, math.factorial(j)) for j in range(9))
if exp3_partial != Fraction(89641,4480) or exp3_partial <= 20:
    raise RuntimeError('strict exp3 positive Taylor lower bound failed')
meanU = s.exp(4*t)
secondU = 5*s.exp(8*t)-4*s.exp(4*t)
equal('Yule_plus4_first_moment_ODE', s.diff(meanU,t)-4*meanU)
equal('Yule_plus4_second_moment_ODE', s.diff(secondU,t)-8*secondU-16*meanU)
equal('Yule_plus4_first_initial', meanU.subs(t,0)-1)
equal('Yule_plus4_second_initial', secondU.subs(t,0)-1)
m,v,M,R = s.symbols('m v M R', positive=True)
next_R = (m*m*(R*M*M)+v*M)/(m*m*M*M)
equal('quenched_second_moment_recurrence',next_R-(R+(v/(m*m))/M))
delta_cases = []
for K in range(1,9):
    # Repeated integration by parts for Gamma(K+1), exact rational kappa points.
    for kap in [s.Rational(1,10),s.Rational(1,2),s.Rational(1),s.Rational(7,3)]:
        delta = kap*s.factorial(K)/(kap+K)**(K+1)
        if not (0 < delta <= 1):
            raise RuntimeError('finite-set extinction lower bound not a probability bound')
    delta_cases.append(K)

# The generic irreducibility words are proved in proof.md; these finite checks
# independently catch wrong directions, net stoichiometry and zero-reactant bugs.
def apply(state, env, reaction_index):
    reaction=MODEL['reactions'][str(env)][reaction_index]
    if any(state[j] < reaction['source'][j] for j in range(3)):
        raise RuntimeError('unavailable reactant on irreducibility path')
    return tuple(state[j]-reaction['source'][j]+reaction['product'][j] for j in range(3))

path_checks=0
for state in itertools.product(range(4), repeat=3):
    current=state
    for _ in range(state[1]):current=apply(current,1,1)
    for _ in range(state[2]):current=apply(current,1,2)
    for _ in range(current[0]):current=apply(current,2,0)
    if current != (0,0,0):raise RuntimeError('to-zero path failed')
    target=state
    current=apply((0,0,0),2,1)
    for _ in range(target[1]):current=apply(current,1,0)
    for _ in range(target[2]):current=apply(current,1,3)
    for _ in range(current[0]):current=apply(current,2,0)
    for _ in range(target[0]):current=apply(current,2,1)
    if current != target:raise RuntimeError('from-zero path failed')
    path_checks+=2

certificate={
    'schema':'switching-shared-species-exact-checks-v1',
    'scope':'Exact source Example44; proof candidate all alpha>0, 0<kappa<=1; generalConjecture43 unresolved',
    'source_reaction_count':8,
    'source_mean_matrices':[[[str(M[i,j]) for j in range(3)] for i in range(3)] for M in [M1,M2]],
    'source_inflows':[[str(b[i]) for i in range(3)] for b in [b1,b2]],
    'active_pair':'A->4A+B at rateA; B->A at rateB; active total increments4 or0',
    'identity_count':len(identities),
    'symbolic_identity_residuals':identities,
    'strict_exact_inequalities':{'sqrt5_gt_20_over9':True,'c_plus_gt_1_over20':True,'lambda_plus_gt3':True,'exp3_positive_partial_sum':str(exp3_partial),'partial_sum_minus20':str(exp3_partial-20),'implies_log20_lt3':True},
    'finite_path_checks':{'all_nonnegative_count_triples_in_0_to3':64,'two_paths_each':path_checks,'generic_all_count_proof':'proof.md Section6; this finite check alone is not a generic theorem'},
    'delta_positive_finite_K_checks':delta_cases,
    'probability_proof_status':'Not formalized by this script; conditional branching, SLLN/Borel-Cantelli, survival and transience rely on the complete argument in proof.md',
    'floating_tests':0,
    'external_BPRE_survival_theorem_required':False,
    'human_or_formal_review_asserted':False
}
output=ROOT/'evidence'/'certificate.json'
output.parent.mkdir(parents=True,exist_ok=True)
output.write_text(json.dumps(certificate,sort_keys=True,indent=2)+'\n')
print(json.dumps(certificate,sort_keys=True))
