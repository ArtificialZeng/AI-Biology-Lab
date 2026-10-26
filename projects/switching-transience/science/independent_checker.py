#!/usr/bin/env python3
"""Portable fresh finite exact checks; this is NOT a probability-proof certificate.

Run with Python 3 standard library only. No proposer code is imported or run.
Default paths resolve relative to this file; outputs stay in this review folder.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import random
import sys

if sys.flags.optimize != 0:
    raise RuntimeError("Run with PYTHONOPTIMIZE=0: exact checks require assertions enabled")

HERE = Path(__file__).resolve().parent
HISTORICAL_DIGEST = '930f8e2af59997fe098e9fdcd646c2c4d2f55d1a5f5f847615bb8574d3f404ce'
INPUT_HASHES = {
    'model.json': '2102a8b0c0f49abd6a8969a8bf697b34b19bf290e0acb7f773b129a37693e31a',
    'proof.md': 'c557f687601012c34a99d15da610f2c5594f53847524eca882686e5d03c8374f',
}
actual = {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in INPUT_HASHES}
assert actual == INPUT_HASHES

# Re-enter the eight arrows from the visually inspected ORIGINAL page 30.
# They are not reconstructed from a mean matrix or copied from proposer code.
arrows = {
    '1': [([1,0,0], [4,1,0], '1'), ([0,1,0], [1,0,0], '1'),
          ([0,0,1], [0,0,0], 'alpha'), ([0,0,0], [0,0,1], '1')],
    '2': [([1,0,0], [0,0,0], 'alpha'), ([0,0,0], [1,0,0], '1'),
          ([0,1,0], [0,0,1], '1'), ([0,0,1], [0,1,4], '1')],
}
model = json.loads((HERE / 'model.json').read_text())
for env in arrows:
    assert [(r['source'], r['product'], r['rate_constant']) for r in model['reactions'][env]] == arrows[env]
    active = (0,1) if env == '1' else (2,1)
    changes = []
    for src, dst, rate in arrows[env]:
        if sum(src) == 1 and any(src[j] for j in active):
            assert rate == '1'
            assert all(dst[j] == 0 for j in range(3) if j not in active)
            changes.append(sum(dst[j]-src[j] for j in active))
    assert sorted(changes) == [0,4]

# Q(sqrt(5)) exact arithmetic. Independent scalar ODE derivation:
# b'=a-b, a'=3a+b implies b''-2b'-4b=0, b(0)=1,b'(0)=-1.
def add(x,y): return (x[0]+y[0], x[1]+y[1])
def mul(x,y): return (x[0]*y[0]+5*x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def scale(c,x): return (c*x[0],c*x[1])
one, zero = (F(1),F(0)), (F(0),F(0))
lp, lm = (F(1),F(1)), (F(1),F(-1))
cp, cm = (F(1,2),F(-1,5)), (F(1,2),F(1,5))
for lam in (lp,lm):
    assert add(add(mul(lam,lam),scale(-2,lam)),(-F(4),F(0))) == zero
assert add(cp,cm) == one
assert add(mul(cp,lp),mul(cm,lm)) == (-F(1),F(0))
# Alternative rational-square proof: 2 < sqrt(5) < 9/4.
# Therefore lambda+>3 and c+=1/2-sqrt(5)/5>1/20.
assert F(4) < 5 < F(81,16)
exp3_lower = sum((F(3**j, math.factorial(j)) for j in range(11)), F(0))
assert exp3_lower > 20

# Pure-birth generator for jump +4: Lz=4z, Lz^2=8z^2+16z.
def expand_delta_square(jump): return (2*jump, jump*jump)
assert expand_delta_square(4) == (8,16)
# Coefficients of exponentials e^(8t),e^(4t): q=5e8-4e4.
assert (8*5,4*(-4)) == (8*5,8*(-4)+16)
assert 5-4 == 1

# Independent factorial-moment route to the normalized second-moment index.
# F'=m^2 F + f M, f=v+m^2-m; E[Y^2]=F+M.
# Finite exact rational substitutions supplement the universal derivation in report.
rng = random.Random(481930)
recurrences = []
for _ in range(37):
    M = F(rng.randrange(1,50),rng.randrange(1,18))
    m = F(rng.randrange(1,50),rng.randrange(1,18))
    v = F(rng.randrange(1,50),rng.randrange(1,18))
    fac = F(rng.randrange(0,50),rng.randrange(1,18))
    nxt = m*m*fac + (v+m*m-m)*M
    before = (fac+M)/(M*M)
    after = (nxt+m*M)/(m*m*M*M)
    assert after-before == v/(m*m*M)
    recurrences.append([str(M),str(m),str(v),str(fac)])

# Held-out paths outside the proposer's complete {0,1,2,3}^3 grid.
# Use S3 births to create S2, reversing the proposer S1 construction.
def fire(state,env,index):
    src,dst,_ = arrows[str(env)][index]
    assert all(state[j]>=src[j] for j in range(3))
    return tuple(state[j]-src[j]+dst[j] for j in range(3))
states = {(0,0,17),(9,0,0),(0,11,0),(17,19,23)}
while len(states)<104:
    states.add(tuple(rng.randrange(4,27) for _ in range(3)))
for state in sorted(states):
    x = state
    for _ in range(x[1]): x=fire(x,2,2)
    for _ in range(x[0]): x=fire(x,2,0)
    for _ in range(x[2]): x=fire(x,1,2)
    assert x == (0,0,0)
    x = fire((0,0,0),1,3)
    for _ in range(state[1]): x=fire(x,2,3)
    for _ in range(x[2]): x=fire(x,1,2)
    for _ in range(state[2]): x=fire(x,1,3)
    for _ in range(state[0]): x=fire(x,2,1)
    assert x == state

delta_checks=[]
for K in (9,11,17,25):
    for kap in (F(1,13),F(2,7),F(5,6),F(1)):
        d=kap*math.factorial(K)/(kap+K)**(K+1)
        assert 0<d<=1
        delta_checks.append({'K':K,'kappa':str(kap),'delta':str(d)})

result = {
    'schema':'fresh-finite-exact-checks-v1',
    'historical_snapshot_digest':HISTORICAL_DIGEST,
    'historical_snapshot_rechecked':False,
    'input_hashes_verified':actual,
    'input_scope':'Only same-directory model.json and proof.md; eight-file private snapshot is not rechecked by this portable variant',
    'source_arrows_reentered_and_matched':8,
    'independent_derivation':'scalar b second-order ODE in Q(sqrt(5)); factorial-moment recurrence; reversed S3 irreducibility paths',
    'exp3_taylor_terms':11,
    'exp3_exact_lower_bound':str(exp3_lower),
    'factorial_recurrence_exact_rational_instances':len(recurrences),
    'heldout_states': [list(x) for x in sorted(states)],
    'heldout_finite_paths':2*len(states),
    'delta_heldout_instances':delta_checks,
    'infinite_probability_argument_formalized':False,
    'finite_checks_prove_universal_claim':False,
    'proposer_code_imported_or_run':False,
    'formal_or_human_certification':False,
    'passed':True,
}
(HERE/'portable-certificate.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps(result,sort_keys=True))
