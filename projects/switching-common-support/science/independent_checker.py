#!/usr/bin/env python3
"""Independent exact finite verification. No probability theorem is machine-proved.

Run from any directory: python3 /path/to/checker.py
Inputs are sibling model.json and proof.md only; Python stdlib suffices.
This public checker verifies only the two input hashes; it does not verify the original six-file snapshot.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

BASE = Path(__file__).resolve().parent
MODEL_SHA = '8f1b8ecb4dc277f15ca5bb698deb600fcfb523bc7ee44efda84349164de99827'
PROOF_SHA = 'b978df02bb0b902489066d21cfff2b27106aacf557d79fca2ca18d4ef4747c68'
for filename, digest in [('model.json',MODEL_SHA),('proof.md',PROOF_SHA)]:
    assert sha256((BASE/filename).read_bytes()).hexdigest() == digest, filename+' changed'
model = json.loads((BASE/'model.json').read_text())
assert model['species'] == ['A','B','C']
assert model['initial_state'] == [0,1,0]
assert model['initial_environment'] == 'either'
assert model['inflows'] is False
assert model['kappa'] == 'every strictly positive value'
assert model['environment'] == 'independent symmetric1<->2 ratekappa'
assert set(model['reactions']) == {'1','2'}

matrices={}; bridge={}; types={}; checked_states=0
witnesses={'1':[2,1,0],'2':[0,1,2]}
for env,reactions in model['reactions'].items():
    matrix=[[F(0) for j in range(3)] for i in range(3)]
    bridge_reactions=[]
    for r in reactions:
        source,target,rate=r['source'],r['product'],F(r['rate'])
        assert len(source)==len(target)==3 and all(isinstance(z,int) and z>=0 for z in source+target)
        assert sum(source)==1 and rate>0
        j=source.index(1)
        delta=[target[i]-source[i] for i in range(3)]
        for i in range(3): matrix[i][j] += rate*delta[i]
        assert target[1] == 0, 'No reaction may produce B'
        if j==1:
            assert target in ([1,0,0],[0,0,1])
            bridge_reactions.append((rate,target))
    assert len(bridge_reactions)==1
    bridge[env]={'total_conversion_hazard':sum(r for r,t in bridge_reactions),'products':[t for r,t in bridge_reactions]}
    assert bridge[env]['total_conversion_hazard']==1
    matrices[env]=matrix
    # Exhaustive independent generator-on-coordinate check on a held-out cube.
    for x in product(range(6),repeat=3):
        drift=[F(0),F(0),F(0)]
        for r in reactions:
            j=r['source'].index(1)
            for i in range(3):
                drift[i] += F(r['rate'])*x[j]*(r['product'][i]-r['source'][i])
        assert drift == [sum(matrix[i][j]*x[j] for j in range(3)) for i in range(3)]
        checked_states+=1
    # After B has disappeared, each remaining one-type subspace is invariant.
    types[env]={}
    for j in [0,2]:
        rates={'birth':F(0),'death':F(0)}
        for r in reactions:
            if r['source'][j]:
                assert all(r['product'][i]==0 for i in range(3) if i!=j)
                if r['product'][j]==2: rates['birth']+=F(r['rate'])
                elif r['product'][j]==0: rates['death']+=F(r['rate'])
                else: raise AssertionError('Unsupported one-type reaction')
        types[env][model['species'][j]]=rates
expected={
 '1':[[1,1,0],[0,-1,0],[0,0,-2]],
 '2':[[-2,0,0],[0,-1,0],[0,1,1]]}
assert matrices==expected
left={e:[sum(F(witnesses[e][i])*m[i][j] for i in range(3)) for j in range(3)] for e,m in matrices.items()}
support={e:[j+1 for j,z in enumerate(v) if z] for e,v in witnesses.items()}
for e in matrices:
    assert left[e]==witnesses[e]
    assert all(left[e][j-1]>0 for j in support[e])
common=sorted(set.intersection(*(set(x) for x in support.values())))
assert common==[2]
# Exact stationary law solves wQ=0 for Q=[[-1,1],[1,-1]], w1+w2=1.
stationary=[F(1,2),F(1,2)]
Q=[[-1,1],[1,-1]]
assert [sum(stationary[i]*Q[i][j] for i in range(2)) for j in range(2)]==[0,0]
log_rates={s:sum(stationary[int(e)-1]*(types[e][s]['birth']-types[e][s]['death']) for e in types) for s in ['A','C']}
assert log_rates=={'A':F(-1,2),'C':F(-1,2)}
# At kappa=1, the annealed one-type first-moment matrix has negative determinant.
# It therefore has one positive eigenvalue despite the quenched log rate -1/2.
annealed_at_kappa1=[[F(0),F(1)],[F(1),F(-3)]]
annealed_det=annealed_at_kappa1[0][0]*annealed_at_kappa1[1][1]-annealed_at_kappa1[0][1]*annealed_at_kappa1[1][0]
assert annealed_det==-1

def clean(v):
    if isinstance(v,F): return str(v)
    if isinstance(v,dict): return {k:clean(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)): return [clean(x) for x in v]
    return v
certificate=clean({
 'ok':True,'method':'independent stdlib exact arithmetic; no probability proof encoded',
 'input_hashes':{'model.json':MODEL_SHA,'proof.md':PROOF_SHA},
 'matrices':matrices,'left_vectors':witnesses,'left_products':left,
 'supports_1_based':support,'common_support_1_based':common,
 'bridge':bridge,'post_conversion_rates':types,
 'stationary_distribution':stationary,'quenched_asymptotic_log_mean_per_time':log_rates,
 'generator_coordinate_checks':checked_states,
 'annealed_stress_at_kappa1':{'matrix':annealed_at_kappa1,'determinant':annealed_det,'positive_eigenvalue':True},
 'scope_limit':'Finite algebra checks hypotheses/invariants only. Almost-sure extinction requires the analytical argument in proof.md.'})
print(json.dumps(certificate,sort_keys=True,indent=2))
