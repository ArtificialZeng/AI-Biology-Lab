#!/usr/bin/env python3
"""Independent finite checks for the frozen serial-transfer argument.

This is supplementary numerical/algebraic evidence, not a quantified proof.
Standalone: no arguments or input files; emits one JSON object on stdout.
Requires Python 3, mpmath, and sympy. It never reads a private research path.
"""
from fractions import Fraction as Q
import json
import random

import mpmath as mp
import sympy as sp

# Recompute the implicit derivative numerator as a formal polynomial.
a, b, u, A, B = sp.symbols('a b u A B')
G = a*u*A + b*(1-u)*B
D = u*A+(1-u)*B
remainder = sp.expand(G-(a-b)*u*(1-u)*(A-B)-(a*u+b*(1-u))*D)
assert remainder == 0

# Exact coexistence cycle: delta=1/8, (a,b)=(1,2) then (2,1),
# u=2/3 -> 1/3 -> 2/3, with exp(tau)=4 in each phase.
delta = Q(1, 8)
cycle = []
uc = Q(2, 3)
derivative_product = Q(1)
for ac, bc in [(1, 2), (2, 1)]:
    pre1, pre2 = delta*uc*4**ac, delta*(1-uc)*4**bc
    assert pre1+pre2 == 1
    vc = pre1
    deriv = (ac*uc+bc*(1-uc))/(ac*vc+bc*(1-vc))
    derivative_product *= deriv
    cycle.append({'a': ac, 'b': bc, 'u': str(uc), 'v': str(vc),
                  'exp_tau': '4', 'derivative': str(deriv)})
    uc = vc
assert uc == Q(2, 3)
assert derivative_product == Q(16, 25)

mp.mp.dps = 90

def phase(zz, aa, bb, dd):
    """Compute the hit time from the two original population equations."""
    aa, bb, dd, zz = map(mp.mpf, (aa, bb, dd, zz))
    freq = 1/(1+mp.exp(-zz))
    complement = 1/(1+mp.exp(zz))
    lo, hi = mp.log(1/dd)/max(aa, bb), mp.log(1/dd)/min(aa, bb)
    if aa == bb:
        tau = lo
    else:
        for _ in range(210):
            mid = (lo+hi)/2
            total = dd*(freq*mp.exp(aa*mid)+complement*mp.exp(bb*mid))
            if total < 1:
                lo = mid
            else:
                hi = mid
        tau = (lo+hi)/2
    pop1 = dd*freq*mp.exp(aa*tau)
    pop2 = dd*complement*mp.exp(bb*tau)
    out = mp.log(pop1/pop2)
    deriv = (aa*freq+bb*complement)/(aa*pop1+bb*pop2)
    return out, tau, deriv, abs(pop1+pop2-1)

ratepairs = [('1','1'), ('1','2'), ('2','1'), ('1','4'), ('3','2'),
             ('2','3'), ('0.01','30'), ('30','0.01'), ('7','11')]
dilutions = ['0.0009765625', '0.125', '0.5', '0.9']
logits = ['-60','-12','-2','0','1','9','60']
max_residual = mp.mpf(0)
max_identity_error = mp.mpf(0)
phase_checks = 0
for dd in dilutions:
    for aa, bb in ratepairs:
        for zz in logits:
            zz = mp.mpf(zz)
            out, tau, der, residual = phase(zz, aa, bb, dd)
            aa0, bb0, L = mp.mpf(aa), mp.mpf(bb), mp.log(1/mp.mpf(dd))
            gain = out-zz
            lower, upper = L*(1-bb0/aa0), L*(aa0/bb0-1)
            tol = mp.mpf('1e-42')*(1+abs(lower)+abs(upper))
            assert lower-tol <= gain <= upper+tol, (dd,aa,bb,str(zz),
                str(lower),str(gain),str(upper),str(tau))
            assert 0 < der <= 1+tol
            if aa != bb:
                assert der < 1
            max_residual = max(max_residual, residual)
            max_identity_error = max(max_identity_error,
                                     abs(gain-(aa0-bb0)*tau))
            phase_checks += 1

# Finite schedule attacks include phase reordering, mixed selection,
# critical U=0 and V=0, coexistence and neutral phases.
rng = random.Random(20261027)
schedules = [
    [('1','4'), ('3','2')],  # strict U<0 despite a favourable second phase
    [('1','2'), ('3','2')],  # U=0, nonneutral
    [('2','1'), ('2','3')],  # V=0, nonneutral
    [('1','2'), ('2','1')],  # interior cycle
    [('1','1'), ('2','2')],  # fully neutral
]
for _ in range(36):
    schedules.append([ratepairs[rng.randrange(len(ratepairs))]
                      for _ in range(1+rng.randrange(7))])
schedule_checks = 0
for sched in schedules:
    for shift in range(len(sched)):
        ordered = sched[shift:]+sched[:shift]
        L = mp.log(8)
        upper = L*sum(mp.mpf(aa)/mp.mpf(bb)-1 for aa,bb in ordered)
        lower = L*sum(1-mp.mpf(bb)/mp.mpf(aa) for aa,bb in ordered)
        for initial in ['-30','-1','0','4','30']:
            zz0, zz = mp.mpf(initial), mp.mpf(initial)
            product = mp.mpf(1)
            for aa, bb in ordered:
                zz, _, der, _ = phase(zz, aa, bb, '0.125')
                product *= der
            tol = mp.mpf('1e-40')*(1+abs(lower)+abs(upper))
            assert lower-tol <= zz-zz0 <= upper+tol
            if any(aa != bb for aa,bb in ordered):
                assert 0 < product < 1
            else:
                assert abs(zz-zz0) < tol
            if upper < -tol:
                assert zz < zz0
            if lower > tol:
                assert zz > zz0
            schedule_checks += 1

result = {
    'schema': 'serial-transfer-independent-finite-checks-v1',
    'dependencies': ['mpmath', 'sympy'],
    'symbolic_derivative_numerator_remainder': str(remainder),
    'exact_coexistence_cycle': cycle,
    'exact_return_derivative': str(derivative_product),
    'phase_checks': phase_checks,
    'schedule_phase_initial_checks': schedule_checks,
    'mpmath_decimal_digits': mp.mp.dps,
    'max_threshold_residual': mp.nstr(max_residual, 8),
    'max_odds_identity_error': mp.nstr(max_identity_error, 8),
    'status': 'finite checks completed; no counterexample found',
    'proof_status': 'supplementary only; no formal or quantified certification',
}
print(json.dumps(result, indent=2))
