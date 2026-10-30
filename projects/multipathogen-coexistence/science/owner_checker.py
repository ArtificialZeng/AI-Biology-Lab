#!/usr/bin/env python3
"""Exact verification of a counterexample to Djuikem--Arino Conjecture 5.

Uses the complete 16-variable source ODE, not numerical trajectories.
All checks are explicit and remain active under Python -O.
"""
from __future__ import annotations
import argparse
import itertools
import json
from pathlib import Path
import sympy as sp


def require(condition, message):
    if not condition:
        raise ValueError(message)


def model(P, lambdas):
    V = len(lambdas)
    S = sp.symbols('S1:' + str(P + 1))
    L = [[sp.Symbol(f'L{p+1}_{v+1}') for v in range(V)] for p in range(P)]
    I = [[sp.Symbol(f'I{p+1}_{v+1}') for v in range(V)] for p in range(P)]
    R = sp.symbols('R1:' + str(P + 1))
    totals = [sum(I[p][v] for p in range(P)) for v in range(V)]
    infection = [sum(lambdas[v] * totals[v] for v in range(V)) for p in range(P)]
    fS = [1 - S[p] * (1 + infection[p]) for p in range(P)]
    fL = [[lambdas[v] * S[p] * totals[v] - 2 * L[p][v] for v in range(V)] for p in range(P)]
    fI = [[L[p][v] - 2 * I[p][v] for v in range(V)] for p in range(P)]
    fR = [sum(I[p]) - R[p] for p in range(P)]
    xs = list(S) + list(itertools.chain.from_iterable(L)) + list(itertools.chain.from_iterable(I)) + list(R)
    fs = fS + list(itertools.chain.from_iterable(fL)) + list(itertools.chain.from_iterable(fI)) + fR
    return S, L, I, R, fs, xs, fS, fL, fI, fR


def equilibrium(P, lambdas, weights):
    S,L,I,R,fs,xs,*_ = model(P, lambdas)
    active = [v for v,t in enumerate(weights) if t > 0]
    if active:
        lam = lambdas[active[0]]
        s = sp.Rational(4, P) / lam
    else:
        s = sp.Integer(1)
    point = {x:s for x in S}
    point.update({L[p][v]:2*weights[v] for p in range(P) for v in range(len(lambdas))})
    point.update({I[p][v]:weights[v] for p in range(P) for v in range(len(lambdas))})
    point.update({R[p]:sum(weights) for p in range(P)})
    residual = [sp.factor(f.subs(point)) for f in fs]
    require(all(x == 0 for x in residual), 'proposed point is not an equilibrium')
    require(all(x >= 0 for x in point.values()), 'point outside nonnegative state space')
    require(all(sp.simplify(point[S[p]] + sum(point[z] for z in L[p]) + sum(point[z] for z in I[p]) + point[R[p]]) == 1 for p in range(P)), 'wrong host totals')
    return point, residual


def certificate():
    P = 2
    lambdas = [sp.Integer(1), sp.Integer(3), sp.Integer(4)]
    V = len(lambdas)
    S,L,I,R,fs,xs,fS,fL,fI,fR = model(P,lambdas)
    require(len(xs) == len(fs) == 16, 'full source state dimension mismatch')
    zero_point = {x:0 for x in xs}
    zero_point.update({x:1 for x in S})
    require(all(sp.expand(f.subs(zero_point)) == 0 for f in fs), 'DFE residual')
    spectra = []
    for lam in lambdas:
        B = sp.ones(P)*lam/4  # source equation(9), not mutant/within-host R0
        eigen = B.eigenvals()
        require(set(eigen) == {sp.Integer(0),lam/2}, 'next generation eigenvalue mismatch')
        spectra.append({'B':[[str(z) for z in row] for row in B.tolist()], 'eigenvalues':sorted(map(str,eigen)), 'R0':str(lam/2)})
    require(lambdas[0]/2 < 1 and all(lam/2 > 1 for lam in lambdas[1:]), 'Conjecture5 strict hypotheses not met')
    identities = []
    for v,lam in enumerate(lambdas):
        Jv = sum(I[p][v] for p in range(P))
        lhs = sum(fL[p][v] + 2*fI[p][v] for p in range(P))
        rhs = Jv*(lam*sum(S) - 4)
        require(sp.expand(lhs-rhs) == 0, 'infection equilibrium identity failed')
        identities.append(str(sp.factor(rhs)))
    total_residuals=[]
    for p in range(P):
        actual=fS[p]+sum(fL[p])+sum(fI[p])+fR[p]
        expected=1-S[p]-sum(L[p])-sum(I[p])-R[p]
        require(sp.expand(actual-expected)==0,'host population identity failed')
        total_residuals.append('0')
    # Independent exhaustive SUPPORT enumeration; continuous variables are resolved
    # analytically by the displayed identity, not by sampling a finite grid.
    supports=[]
    for bits in itertools.product([0,1],repeat=V):
        active=[v for v,x in enumerate(bits) if x]
        distinct={lambdas[v] for v in active}
        if len(distinct)>1:
            classification='impossible: incompatible required susceptible sums'
            target=None
        elif active and P*lambdas[active[0]]<=4:
            classification='impossible: required susceptible mass >= disease-free host total'
            target=str(4/lambdas[active[0]])
        elif active:
            classification='feasible singleton-strain equilibrium'
            target=str(4/lambdas[active[0]])
        else:
            classification='disease-free equilibrium'
            target=None
        supports.append({'active_strains':[v+1 for v in active],'required_S_total':target,'classification':classification})
    require(supports[3]['active_strains']==[2,3] and supports[3]['classification'].startswith('impossible'), 'wrong target support')
    # The source conclusion requires support {2,3}. It would force T=4/3 and T=1.
    T=sp.Symbol('T')
    bezout=sp.expand(4*(3*T-4)-3*(4*T-4))
    require(bezout == -4,'nonexistence contradiction identity')
    pts=[]
    for v in [1,2]:
        s=sp.Rational(4,P)/lambdas[v]
        weights=[sp.Integer(0)]*V
        weights[v]=(1-s)/4
        pt,res=equilibrium(P,lambdas,weights)
        pts.append({'active_strain':v+1,'point':{str(k):str(val) for k,val in pt.items()},'residuals':list(map(str,res))})
    winning_point,_=equilibrium(P,lambdas,[0,0,sp.Rational(1,8)])
    J=sp.Matrix(fs).jacobian(xs).subs(winning_point)
    z=sp.Symbol('z')
    charpoly=sp.factor(J.charpoly(z).as_expr())
    expected=(z+1)**3*(z+2)**7*(z+3)*(z*z+4*z+1)*(z**3+6*z*z+8*z+4)
    require(sp.expand(charpoly-expected)==0,'wrong full 16D characteristic polynomial')
    require(6>0 and 8>0 and 4>0 and 6*8>4,'cubic Hurwitz certificate failed')
    # Held controls: equal supercritical lambda permits coexisting positive strains;
    # below/equal threshold cannot give a positive endemic equilibrium.
    equal_point,_=equilibrium(P,[1,4,4],[0,sp.Rational(1,16),sp.Rational(1,16)])
    require(all(equal_point[I[p][v]]>0 for p in range(P) for v in [1,2]),'equal-rate coexistence control failed')
    # Different positive host labels / pathogen labels change no claim.
    for perm in itertools.permutations(range(V)):
        ls=[lambdas[v] for v in perm]
        weights=[sp.Rational(1,8) if v==2 else sp.Integer(0) for v in perm]
        equilibrium(P,ls,weights)
    for hosts,lam in [(1,8),(3,2),(4,3),(5,1)]:
        s=sp.Rational(4,hosts)/lam
        equilibrium(hosts,[sp.Integer(lam)],[(1-s)/4])
    fullrank_beta = [sp.Matrix([[1,sp.Rational(1,2)],[sp.Rational(1,2),1]]),
                     sp.Matrix([[4,1],[1,4]]),sp.Matrix([[5,2],[2,5]])]
    fullrank_R0 = [sp.Rational(3,8),sp.Rational(5,4),sp.Rational(7,4)]
    for B,rho in zip(fullrank_beta,fullrank_R0):
        require(B.det() != 0 and all(x>0 for x in B),'fullrank positive contact control')
        require(max(B.eigenvals())/4 == rho,'fullrank R0 error')
    require(all(x>0 for x in fullrank_beta[2]-fullrank_beta[1]),'strict dominance not satisfied')
    s1,s2=sp.symbols('s1 s2')
    T2=sp.diag(s1,s2)*fullrank_beta[1]/4
    T3=sp.diag(s1,s2)*fullrank_beta[2]/4
    d2=sp.factor((sp.eye(2)-T2).det())
    d3=sp.factor((sp.eye(2)-T3).det())
    require(sp.expand(d3-sp.Rational(5,4)*d2-(sp.Rational(9,64)*s1*s2-sp.Rational(1,4)))==0,'fullrank elimination identity')
    fullrank={'beta':[[[str(x) for x in row] for row in B.tolist()] for B in fullrank_beta],
              'R0':list(map(str,fullrank_R0)),
              'det_I_minus_T2':str(d2),'det_I_minus_T3':str(d3),
              'elimination':'d3-(5/4)*d2=(9/64)*S1*S2-1/4; d2=d3=0 force S1*S2=16/9 and S1+S2=8/3>2',
              'all_contact_matrices_full_rank':True,'strict_effective_kernel_dominance':True}
    return {'schema':'exact-multipathogen-counterexample-v1','source_model':'Djuikem-Arino2025 equation(1), sourceR0 equations(8)-(9), Conjecture5','scope':'P=2,V=3,k=1; all model parameters strictly positive; beta_pqv=lambda_v including cross-host contacts','state_dimension':16,'parameters':{'b_p':'1','d_p':'1','epsilon_pv':'1','gamma_pv':'1','lambda_v':['1','3','4'],'beta_pqv':'lambda_v for every p,q'},'next_generation':spectra,'fullrank_positive_witness':fullrank,'exact_polynomial_infection_identities':identities,'host_total_identity_residuals':total_residuals,'target_support':[2,3],'incompatible_susceptible_sums':['4/3','1'],'bezout_contradiction':'4*(3*T-4)-3*(4*T-4)=-4','equilibrium_support_classification':supports,'feasible_endemic_equilibria':pts,'winning_equilibrium_charpoly':str(charpoly),'winning_equilibrium_stability':'all16 eigenvalues have strictly negative real part, exact factor and Routh-Hurwitz proof','invasion_R_ratios':['1/4','3/4'],'controls':{'equal_lambda_coexistence':True,'all_six_pathogen_permutations':6,'additional_host_sizes':[1,3,4,5]},'evidence_kind':'exact identities and exact symbolic linear algebra; no floating point inference','limits':['Source-specific counterexample; competitive exclusion is established biology','General heterogeneous multi-pathogen dynamics not classified','No global convergence or wet-lab prediction','No author confirmation or universal priority claim']}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path);args=ap.parse_args()
    data=certificate();out=json.dumps(data,sort_keys=True,indent=2)+'\n'
    if args.output:
        with args.output.open('x',encoding='utf-8') as f:f.write(out)
    else:print(out,end='')

if __name__=='__main__':main()
