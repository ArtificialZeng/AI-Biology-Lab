#!/usr/bin/env python3
"""Exact modular counterexample search; not a proof of an infinite theorem.

Uses exhaustive vertex-path set packing, harmonic linear solves, and a separately
computed physical Markov-parameter Jacobian. Mod-p nonzero minors certify
rational nonzero minors at the recorded positive integer parameter points.
"""
from __future__ import annotations
import argparse
from functools import lru_cache
import hashlib
import itertools
import json
from pathlib import Path
import random
import sys
import networkx as nx
import sympy as sp
P = 1000003
BASE = Path(__file__).resolve().parent

def rank(rows, p=P):
    a = [list(r) for r in rows]
    if not a: return 0
    col = 0
    for j in range(len(a[0])):
        z = next((i for i in range(col,len(a)) if a[i][j]%p),None)
        if z is None: continue
        a[col],a[z]=a[z],a[col]
        inv=pow(a[col][j]%p,-1,p)
        a[col]=[(x*inv)%p for x in a[col]]
        for i in range(col+1,len(a)):
            f=a[i][j]%p
            if f: a[i]=[(x-f*y)%p for x,y in zip(a[i],a[col])]
        col+=1
        if col==len(a): break
    return col

def solve(a,b,p=P):
    n=len(a)
    if not n: return []
    m=len(b[0]); c=[list(x)+list(y) for x,y in zip(a,b)]
    for j in range(n):
        z=next((i for i in range(j,n) if c[i][j]%p),None)
        if z is None: raise ArithmeticError('Dirichlet determinant zero modulo p')
        c[j],c[z]=c[z],c[j]
        v=pow(c[j][j]%p,-1,p); c[j]=[(x*v)%p for x in c[j]]
        for i in range(n):
            if i==j:continue
            v=c[i][j]%p
            if v:c[i]=[(x-v*y)%p for x,y in zip(c[i],c[j])]
    return [r[n:n+m] for r in c]

def verts(mask,n):return [i for i in range(n) if (mask>>i)&1]

def path_masks(t):
    n=len(t)
    return {(i,j):sum(1<<v for v in nx.shortest_path(t,i,j)) for i in range(n) for j in range(n)}

def linkage(obs,leaks,pm):
    obs=tuple(obs);leaks=tuple(leaks)
    @lru_cache(None)
    def rec(i,used):
        if i==len(obs):return 0
        ans=rec(i+1,used)
        for b in leaks:
            mask=pm[obs[i],b]
            if not mask&used:ans=max(ans,1+rec(i+1,used|mask))
        return ans
    return rec(0,0)

def parameters(t,leaks,seed):
    n=len(t); rng=random.Random(seed)
    arcs=sorted((j,i) for a,b in t.edges() for j,i in [(a,b),(b,a)])
    rates={a:rng.randrange(1,100) for a in arcs}
    ell={b:rng.randrange(1,100) for b in leaks}
    return arcs,rates,ell

def harmonic(t,leaks,rates):
    n=len(t); leaks=list(leaks); inside=[j for j in range(n) if j not in leaks]
    h=[[int(j==b) for b in leaks] for j in range(n)]
    a=[[sum(rates[j,i] for i in t[j]) if j==v else -rates.get((j,v),0) for v in inside] for j in inside]
    rhs=[[rates.get((j,b),0) for b in leaks] for j in inside]
    for j,row in zip(inside,solve(a,rhs)):h[j]=row
    return h

def jacobian(t,leaks,r,arcs,rates,ell,p=P):
    n=len(t); pars=arcs+[(b,b) for b in leaks]; m=len(pars)
    a=[[0]*n for _ in range(n)]
    for (j,i),rate in rates.items():a[i][j]+=rate;a[j][j]-=rate
    for b,x in ell.items():a[b][b]-=x
    v=[int(j==r) for j in range(n)]
    dv=[[0]*m for _ in range(n)]; rows=[[] for _ in range(n)]
    for power in range(1,2*n):
        nd=[[sum(a[i][j]*dv[j][k] for j in range(n))%p for k in range(m)] for i in range(n)]
        for k,(j,i) in enumerate(arcs):
            nd[i][k]=(nd[i][k]+v[j])%p;nd[j][k]=(nd[j][k]-v[j])%p
        for k,b in enumerate(leaks,len(arcs)):nd[b][k]=(nd[b][k]-v[b])%p
        v=[sum(a[i][j]*v[j] for j in range(n))%p for i in range(n)]
        dv=nd
        for o in range(n):rows[o].append(dv[o][:])
    return rows

def trees(n):
    if n==1:yield nx.empty_graph(1)
    else:yield from nx.nonisomorphic_trees(n)

def run(harmonic_n=7,jacobian_n=5,heldout=100):
    counts={'harmonic_placements':0,'physical_jacobian_placements':0,'heldout_cases':0,'trees_by_n':{},'physical_seed_retries':0}
    digest=hashlib.sha256()
    for n in range(1,max(harmonic_n,jacobian_n)+1):
        ts=list(trees(n));counts['trees_by_n'][str(n)]=len(ts)
        for ti,t in enumerate(ts):
            pm=path_masks(t)
            for bm in range(1,1<<n):
                b=verts(bm,n)
                arcs,rates,ell=parameters(t,b,1729+100000*n+1000*ti+bm)
                h=harmonic(t,b,rates)
                qs={om:linkage(verts(om,n),b,pm) for om in range(1,1<<n)}
                if n<=harmonic_n:
                    for om,q in qs.items():
                        got=rank([h[o] for o in verts(om,n)])
                        if not (got==q): raise AssertionError(('harmonic',n,ti,bm,om,got,q))
                        counts['harmonic_placements']+=1
                        digest.update(f'H,{n},{ti},{bm},{om},{got};'.encode())
                if n<=jacobian_n:
                    for r in range(n):
                        pending={om for om in qs if (om>>r)&1}
                        for trial in range(4):
                            if trial:
                                arcs,rates,ell=parameters(t,b,982451653+100000*n+1000*ti+bm+trial*10000000)
                                counts['physical_seed_retries']+=1
                            rows=jacobian(t,b,r,arcs,rates,ell)
                            for om in list(pending):
                                expected=2*n-2+qs[om]
                                got=rank([row for o in verts(om,n) for row in rows[o]])
                                if not (got<=expected): raise AssertionError(('physical upper bound',n,ti,bm,om,r,got,expected))
                                if got==expected:
                                    pending.remove(om);counts['physical_jacobian_placements']+=1
                                    digest.update(f'J,{n},{ti},{bm},{om},{r},{got};'.encode())
                            if not pending:break
                        if not (not pending): raise AssertionError(('physical deficient',n,ti,bm,r,sorted(pending)))
    rng=random.Random(20261026)
    for case in range(heldout):
        n=rng.randrange(6,17)
        t=nx.from_prufer_sequence([rng.randrange(n) for _ in range(n-2)])
        b=sorted(rng.sample(range(n),rng.randrange(1,n+1)))
        r=rng.randrange(n);obs=sorted(set([r]+rng.sample(range(n),rng.randrange(1,n+1))))
        pm=path_masks(t);q=linkage(obs,b,pm)
        arcs,rates,ell=parameters(t,b,47000+case)
        h=harmonic(t,b,rates)
        if not (rank([h[o] for o in obs])==q): raise AssertionError('Verification failed')
        rows=jacobian(t,b,r,arcs,rates,ell)
        got=rank([row for o in obs for row in rows[o]])
        if not (got==2*n-2+q): raise AssertionError(('heldout',case,n,b,obs,r,q,got))
        counts['heldout_cases']+=1
        digest.update(f'X,{case},{n},{b},{obs},{r},{got};'.encode())
    return {'status':'finite-checks-passed','prime':P,'arithmetic':'exact modular integer arithmetic; each rank equality is a finite observation; theorem proof is separate',
            'harmonic_max_n':harmonic_n,'physical_jacobian_max_n':jacobian_n,'heldout_n_range':[6,16],
            'coverage':'one representative of every unlabelled tree; every nonempty leak set and every nonempty observation set; every retained root for physical checks; labelled placements on each representative',
            'counts':counts,'result_digest_sha256':digest.hexdigest()}

def rational_witness():
    t=nx.Graph([(0,1),(1,2),(1,3),(3,4),(3,5),(5,6)])
    b=[0,2,4,6]; arcs,rates,ell=parameters(t,b,3571)
    n=len(t);inside=[v for v in range(n) if v not in b]
    a=sp.Matrix([[sum(rates[j,i] for i in t[j]) if j==v else -rates.get((j,v),0) for v in inside] for j in inside])
    rhs=sp.Matrix([[rates.get((j,l),0) for l in b] for j in inside])
    h=sp.zeros(n,len(b))
    for j,row in zip(inside,(a.inv()*rhs).tolist()):h[j,:]=sp.Matrix([row])
    for k,j in enumerate(b):h[j,k]=1
    pm=path_masks(t)
    obs_cases=[[1],[1,3],[1,3,5],[0,1,3,5],[0,1,2,4],[0,1,2,4,6]]
    result=[]
    for obs in obs_cases:
        sub=h[obs,:];q=linkage(obs,b,pm)
        if not (sub.rank()==q): raise AssertionError('Verification failed')
        result.append({'observations':obs,'harmonic_rank':q,'generic_physical_rank':2*n-2+q,'generic_ambiguity_dimension':len(b)-q})
    return {'n':n,'edges':sorted(sorted(e) for e in t.edges()),'root':1,'leaks':b,'rates':{f'{j}->{i}':v for (j,i),v in rates.items()},'leak_rates':ell,
            'harmonic_matrix':[[str(x) for x in row] for row in h.tolist()],'placements':result,
            'arithmetic':'exact rational harmonic ranks; physical rank labels follow theorem, not independently rational Jacobians'}

if __name__ == '__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--harmonic-max-n',type=int,default=7)
    ap.add_argument('--jacobian-max-n',type=int,default=5)
    ap.add_argument('--heldout',type=int,default=100)
    ap.add_argument('--output',type=Path,help='Explicit output directory; default writes only stdout')
    args=ap.parse_args()
    if not (1<=args.harmonic_max_n<=7 and 1<=args.jacobian_max_n<=7 and args.heldout>=0):
        ap.error('n must be between 1 and 7; heldout must be nonnegative')
    out={'finite_checks':run(args.harmonic_max_n,args.jacobian_max_n,args.heldout),'rational_witness':rational_witness()}
    if args.output:
        args.output.mkdir(parents=True,exist_ok=True)
        (args.output/'certificate.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
