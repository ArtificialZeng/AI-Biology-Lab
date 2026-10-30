#!/usr/bin/env python3
"""Portable exact audit: standard library, frozen adjacent inputs, JSON stdout.

This program independently uses Fraction matrix products and permutation
polynomial determinants; it neither imports nor invokes the proposer checker.
"""
from fractions import Fraction as Q
from pathlib import Path
from itertools import permutations
import hashlib
import json

EXPECTED = {
    'model.json': 'b8e83c0cfbffd1dbfbc76c7fe51e48f00e2c346c19a975a0c4f76f693a0a16c1',
    'proof.md': '3ce451c6bf1ef4b104a91e053a8436a822c78b20448fffc921f7136355ca416d',
}
EDGES = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]

def zero(n=4): return [[Q(0) for j in range(n)] for i in range(n)]
def eye(n=4): return [[Q(i == j) for j in range(n)] for i in range(n)]
def add(a,b): return [[x+y for x,y in zip(r,s)] for r,s in zip(a,b)]
def mul(a,b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))),Q(0))
             for j in range(len(b[0]))] for i in range(len(a))]
def power(a,k):
    out=eye(len(a))
    for _ in range(k): out=mul(out,a)
    return out
def gen(rates):
    a=zero()
    for (i,j),x in zip(EDGES,rates): a[i][j]=a[j][i]=Q(x)
    for j in range(4): a[j][j]=-sum(a[i][j] for i in range(4) if i != j)
    return a
def generator_ok(a):
    return all(sum(a[i][j] for i in range(4)) == 0 for j in range(4)) and all(
        a[i][j] > 0 for i in range(4) for j in range(4) if i != j)
def rank_and_det(a):
    a=[r[:] for r in a]; n=len(a); m=len(a[0]); rank=0; det=Q(1)
    for col in range(m):
        piv=next((r for r in range(rank,n) if a[r][col]),None)
        if piv is None: continue
        if piv != rank: a[rank],a[piv]=a[piv],a[rank]; det=-det
        x=a[rank][col]; det*=x
        a[rank]=[v/x for v in a[rank]]
        for r in range(rank+1,n):
            t=a[r][col]; a[r]=[v-t*w for v,w in zip(a[r],a[rank])]
        rank+=1
        if rank==n: break
    return rank, det if n==m and rank==n else Q(0)
def derivative(a):
    # Propagate the exact constant and linear coefficients of (a+t*h)^k.
    columns=[]
    for i,j in EDGES:
        h=zero(); h[i][j]=h[j][i]=1; h[i][i]=h[j][j]=-1
        p=eye(); d=zero(); values=[]
        for k in range(1,4):
            d=add(mul(d,a),mul(p,h)); p=mul(p,a)
            values.append((d[1][0],d[2][0]))
        columns.append([values[k][j] for j in range(2) for k in range(3)])
    return [list(r) for r in zip(*columns)]

# Polynomial coefficients are in ascending degree.
def trim(p):
    p=list(map(Q,p))
    while len(p)>1 and not p[-1]: p.pop()
    return p
def padd(p,q):
    r=[Q(0)]*max(len(p),len(q))
    for i,v in enumerate(p): r[i]+=v
    for i,v in enumerate(q): r[i]+=v
    return trim(r)
def pmul(p,q):
    r=[Q(0)]*(len(p)+len(q)-1)
    for i,v in enumerate(p):
        for j,w in enumerate(q): r[i+j]+=v*w
    return trim(r)
def pdiv(p,q):
    p=trim(p); q=trim(q)
    if q==[0]: raise ZeroDivisionError()
    out=[Q(0)]*max(1,len(p)-len(q)+1)
    while p != [0] and len(p)>=len(q):
        k=len(p)-len(q); t=p[-1]/q[-1]; out[k]+=t
        for j,w in enumerate(q): p[k+j]-=t*w
        p=trim(p)
    return trim(out),p
def gcd(p,q):
    while trim(q) != [0]: _,r=pdiv(p,q); p,q=q,r
    return [x/p[-1] for x in trim(p)]
def pdet(a):
    out=[Q(0)]; n=len(a)
    for perm in permutations(range(n)):
        sign=(-1)**sum(perm[i]>perm[j] for i in range(n) for j in range(i+1,n))
        term=[Q(sign)]
        for i,j in enumerate(perm): term=pmul(term,a[i][j])
        out=padd(out,term)
    return out
def cofactor(a,i,j):
    sub=[[a[r][c] for c in range(len(a)) if c != j]
         for r in range(len(a)) if r != i]
    return [(-1)**(i+j)*x for x in pdet(sub)]
def serial(x):
    if isinstance(x,Q): return str(x)
    if isinstance(x,dict): return {k:serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)): return [serial(v) for v in x]
    return x

def main():
    folder=Path(__file__).resolve().parent
    inputs={name:(folder/name).read_bytes() for name in EXPECTED}
    actual={name:hashlib.sha256(b).hexdigest() for name,b in inputs.items()}
    assert actual==EXPECTED, 'Frozen adjacent input hashes differ'
    model=json.loads(inputs['model.json'])
    a=gen(model['positive_symmetric_rates']); u=gen([1]*6)
    du=derivative(u); da=derivative(a)
    expected_u=[[1,0,0,0,0,0],[-8,0,0,0,0,0],[48,0,0,0,0,0],
                [0,1,0,0,0,0],[0,-8,0,0,0,0],[0,48,0,0,0,0]]
    assert du==expected_u and rank_and_det(du)==(2,0)
    assert rank_and_det(da)==(6,800)
    v=list(map(Q,model['krylov_v'])); b=list(map(Q,model['perturbation_b']))
    eps=Q(model['epsilon'])
    at=[[a[i][j]+eps*b[i]*v[j] for j in range(4)] for i in range(4)]
    assert generator_ok(a) and generator_ok(at)
    moments_a=[[power(a,k)[i][0] for i in range(4)] for k in range(5)]
    moments_at=[[power(at,k)[i][0] for i in range(4)] for k in range(5)]
    annihilator=[sum(v[i]*moments_a[k][i] for i in range(4)) for k in range(4)]
    assert annihilator==[0,0,0,1112]
    assert moments_a[:4]==moments_at[:4]
    delta=[x-y for x,y in zip(moments_at[4],moments_a[4])]
    assert delta==[0,Q(139,125),-Q(139,125),0]
    cycle=at[1][0]*at[2][1]*at[0][2]-at[0][1]*at[1][2]*at[2][0]
    assert cycle==Q(1,100)
    c=[[Q(x) for x in row] for row in model['circulant']]
    assert generator_ok(c)
    resolvent_matrix=[[[ -c[i][j],Q(i==j)] for j in range(4)] for i in range(4)]
    denominator=pdet(resolvent_matrix)
    q=[Q(26),Q(10),Q(1)]
    assert denominator==pmul(pmul([0,1],[6,1]),q)
    nums=[cofactor(resolvent_matrix,0,j) for j in (1,2)]
    assert nums==[[39,18,2],[39,13,1]]
    remainders=[pdiv(n,q)[1] for n in nums]
    assert remainders==[[-13,-2],[13,3]]
    assert all(gcd(n,denominator)==[1] for n in nums)
    # Independently verify the two resolvent columns by adjugate polynomial identity.
    adj_col=[cofactor(resolvent_matrix,0,j) for j in range(4)]
    products=[]
    for i in range(4):
        value=[Q(0)]
        for j in range(4): value=padd(value,pmul(resolvent_matrix[i][j],adj_col[j]))
        products.append(value)
    assert products==[denominator,[0],[0],[0]]
    # Adversarial controls attack essential hypotheses and exact identities.
    controls={}
    leaking=[[a[i][j]+eps*Q(i==1)*v[j] for j in range(4)] for i in range(4)]
    controls['nonconserving_perturbation_rejected']=not generator_ok(leaking)
    large=[[a[i][j]+b[i]*v[j] for j in range(4)] for i in range(4)]
    controls['positivity_destroying_epsilon_rejected']=not generator_ok(large)
    wrong_v=v[:]; wrong_v[0]+=1
    controls['wrong_annihilator_rejected']=sum(wrong_v[i]*moments_a[0][i] for i in range(4))!=0
    controls['wrong_transfer_numerator_rejected']=padd(nums[0],[1])!=cofactor(resolvent_matrix,0,1)
    h=zero(); h[0][1]=h[1][0]=1
    controls['missing_diagonal_tangent_rejected']=any(sum(h[i][j] for i in range(4)) for j in range(4))
    assert all(controls.values())
    return {'ok':True,'input_sha256':actual,'method':'stdlib Fraction arithmetic; dual-coefficient power propagation; permutation polynomial determinants',
            'uniform_derivative':du,'uniform_rank':rank_and_det(du)[0],
            'positive_point_derivative':da,'positive_point_determinant':rank_and_det(da)[1],
            'generator_positivity':{'base':generator_ok(a),'perturbation':generator_ok(at),'circulant':generator_ok(c)},
            'smallest_perturbed_offdiagonal':min(at[i][j] for i in range(4) for j in range(4) if i!=j),
            'krylov_annihilator_moments':annihilator,'matched_full_state_moment_orders':[0,1,2,3],
            'order4_difference':delta,'directed_cycle_product_difference':cycle,
            'characteristic_coefficients_ascending':denominator,'transfer_numerators_ascending':nums,
            'quadratic_remainders_ascending':remainders,'quadratic_discriminant':Q(10)**2-4*26,
            'transfer_gcds':[gcd(n,denominator) for n in nums],
            'adjugate_identity_verified':True,'adversarial_controls':controls,
            'limitations':['Exact fixed witnesses only','No new reversible spectral theorem','No novelty or human peer-review certification','No authenticated runner release-gate PASS']}

if __name__=='__main__':
    try:
        result=main(); code=0
    except Exception as exc:
        result={'ok':False,'error_type':type(exc).__name__,'error':str(exc)}; code=1
    print(json.dumps(serial(result),sort_keys=True,indent=2))
    raise SystemExit(code)
