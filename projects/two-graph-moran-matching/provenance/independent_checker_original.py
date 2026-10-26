#!/usr/bin/env python3
"""Fresh referee: raw parent-target derivation, with no owner-code imports.

Uses the four labelled deleted-edge endpoints and symbolic counts of remaining
mutants/residents; it does not start from the proposed six-state rate table.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib, json
import sympy as s

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'matching-fixation-research'
snapshot=json.loads((SOURCE/'snapshot.json').read_text())
def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()
assert all(digest(SOURCE/p)==h for p,h in snapshot['files'].items())
assert hashlib.sha256(json.dumps(snapshot['files'],sort_keys=True,separators=(',',':')).encode()).hexdigest()==snapshot['digest']

n,c,r=s.symbols('n c r')
def V(k,a): return (n-2)*k+a+2*k*(n-k)/(n*(n-1))

def raw_symbolic(bits):
    a=sum(bits[i]*bits[i+1] for i in (0,2))
    k=sum(bits)+c
    old=V(k,a)
    qv=0
    # Individually labelled endpoint parents and endpoint targets.
    for i,typ in enumerate(bits):
        degree=n-2 if typ else n-1
        for j,target in enumerate(bits):
            if j==i or (typ and j==(i^1)) or typ==target: continue
            delta=typ-target
            qv+=(V(k+delta,a+delta*bits[j^1])-old)/degree
        # The endpoint parent replaces an unpaired opposite-type target.
        targets=r if typ else c
        qv+=targets*(V(k+(1 if typ else -1),a)-old)/degree
    # Unpaired parents replace endpoints of the other type.
    for j,typ in enumerate(bits):
        if typ:
            qv+=r*(V(k-1,a-bits[j^1])-old)/(n-1)
        else:
            qv+=c*(V(k+1,a+bits[j^1])-old)/(n-1)
    # Unpaired parents replace unpaired targets of the other type.
    qv+=c*r*(V(k+1,a)+V(k-1,a)-2*old)/(n-1)
    return s.factor(qv)

representatives={(0,0):(0,0,0,0),(0,1):(1,0,0,0),(0,2):(1,0,1,0),
                 (1,0):(1,1,0,0),(1,1):(1,1,1,0),(2,0):(1,1,1,1)}
stored=json.loads((SOURCE/'evidence/six-polynomial-certificates.json').read_text())
expected={(z['a'],z['b']):sum(t['coefficient']*c**t['powers'][0]*r**t['powers'][1] for t in z['terms']) for z in stored['symbolic_identities']}
polys={}
for ab,bits in representatives.items():
    qv=raw_symbolic(bits)
    poly=s.Poly(s.cancel((-n*(n-2)*(n-1)**2*qv).subs(n,c+r+4)),c,r)
    assert s.expand(poly.as_expr()-expected[ab])==0
    assert all(coef>=0 for coef in poly.coeffs())
    polys[ab]=poly.as_expr()
assert s.cancel(V(0,0))==0
C=n*(n-2)+2
assert s.cancel(V(n,2)-C)==0
assert s.cancel(V(1,0)-C/n)==0
assert s.cancel(raw_symbolic((1,0,0,0)).subs({c:0,r:n-4})+4/(n*(n-1)))==0
assert s.cancel(raw_symbolic((0,0,0,0)).subs({c:1,r:n-5})+4/(n*(n-1)))==0

def vf(mask,N):
    k=bin(mask).count('1')
    a=sum(bool(mask&(1<<i)) and bool(mask&(1<<(i+1))) for i in (0,2))
    return F((N-2)*k+a)+F(2*k*(N-k),N*(N-1))

coverage=[]
for N in (9,10):
    values=[vf(mask,N) for mask in range(1<<N)]
    worst=F(-1000); zero=0; singleton=0
    for mask in range(1<<N):
        bits=[(mask>>i)&1 for i in range(N)]
        drift=F(0)
        for parent in range(N):
            typ=bits[parent]
            targets=[j for j in range(N) if j!=parent and not (typ and parent<4 and j==(parent^1))]
            for target in targets:
                after=(mask | (1<<target)) if typ else (mask & ~(1<<target))
                drift+=(values[after]-values[mask])/F(N*len(targets))
        a=sum(bits[i]*bits[i+1] for i in (0,2))
        b=sum(bits[i]!=bits[i+1] for i in (0,2))
        cv=sum(bits[4:]); rv=N-4-cv
        cert=int(polys[a,b].subs({c:cv,r:rv}))
        assert drift==-F(cert,N*N*(N-2)*(N-1)**2)
        assert drift<=0
        if bin(mask).count('1')==1:
            assert drift==-F(4,N*N*(N-1));singleton+=1
        worst=max(worst,drift)
        zero+=drift==0
    coverage.append({'n':N,'labelled_states':1<<N,'singleton_states':singleton,'maximum_drift':str(worst),'zero_drift_states':zero})

# Exact neutral identical-graph control: total mutant count must have zero drift.
for N in (4,9):
    for mask in range(1<<N):
        drift=F(0)
        for i in range(N):
            for j in range(N):
                if i!=j:drift+=F(((mask>>i)&1)-((mask>>j)&1),N*(N-1))
        assert drift==0

assert all(digest(SOURCE/p)==h for p,h in snapshot['files'].items())
out={'snapshot_digest':snapshot['digest'],'files_verified':len(snapshot['files']),
     'raw_generator':'Q=N(P-I); parent weights 1/d_type(parent)',
     'symbolic_scope':'all n>=4, exactly two vertex-disjoint deletions',
     'six_polynomials':{str(ab):str(p) for ab,p in polys.items()},
     'boundaries_and_singleton_drift':True,'heldout_labelled_audit':coverage,
     'identical_complete_graph_k_drift_control':'zero for every n=4,9 state',
     'frozen_inputs_unchanged_after':True,'sympy_version':s.__version__}
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
