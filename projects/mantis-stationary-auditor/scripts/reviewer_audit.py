"""Fresh reviewer oracle. No owner oracle or Mantis reachability/generator imported."""
from collections import deque
from fractions import Fraction as F
from math import factorial,prod
import inspect,json,numpy as np
from mantis import CRNetwork
from mantis.parsing import parse_reactions
from mantis.stochastic_stationary import _reachable_finite_class,StationaryDistribution

def check(ok,message):
    if not ok:raise RuntimeError(message)
def oracle(start,reactions):
    edges={};todo=deque([start]);seen={start}
    while todo:
        x=todo.popleft();out={}
        for a,b,k in reactions:
            if all(xi>=ai for xi,ai in zip(x,a)):
                y=tuple(xi-ai+bi for xi,ai,bi in zip(x,a,b));r=F(k)*prod(factorial(xi)//factorial(xi-ai) for xi,ai in zip(x,a))
                if r and y!=x:
                    out[y]=out.get(y,F(0))+r
                    if y not in seen:seen.add(y);todo.append(y)
        edges[x]=out
    states=sorted(seen);n=len(states);q=[[F(0)]*n for _ in states]
    for i,x in enumerate(states):
        for y,k in edges[x].items():q[i][states.index(y)]+=k;q[i][i]-=k
    a=[[q[j][i] for j in range(n)]+[F(0)] for i in range(n)];a[-1]=[F(1)]*n+[F(1)]
    for i in range(n):
        pivot=next(j for j in range(i,n) if a[j][i]);a[i],a[pivot]=a[pivot],a[i];p=a[i][i];a[i]=[v/p for v in a[i]]
        for j in range(n):
            if j!=i:
                p=a[j][i];a[j]=[v-p*w for v,w in zip(a[j],a[i])]
    pi=[r[-1] for r in a];check(sum(pi)==1 and all(p>0 for p in pi),'normalization/positivity')
    check(all(sum(pi[i]*q[i][j] for i in range(n))==0 for j in range(n)),'exact stationarity')
    return dict(zip(states,pi))

def compare(text,allowed,start,spec,reactions,means,chem=None):
    rx=parse_reactions([text] if isinstance(text,str) else text);rates={r.rate_key:float(t[2]) for r,t in zip(rx,reactions)}
    support=_reachable_finite_class(np.array(allowed),np.array(start),rx,spec,rates,chem or {})
    ref=oracle(start,reactions);check(set(map(tuple,support))==set(ref),'support mismatch')
    law=StationaryDistribution(spec,np.array(means,dtype=float),1.,True,True,_states=support)
    err=max(abs(law.probability(dict(zip(spec,x)))-float(ref.get(x,0))) for x in allowed)
    check(err<1e-10,'probability mismatch');check(law.probability(dict(zip(spec,[-1]*len(spec))))==0,'negative query')
    check(law.probability(dict(zip(spec,[sum(start)+1]*len(spec))))==0,'outside query')
    cache=law._support;law.probability(dict(zip(spec,start)));check(law._support is cache,'cache reuse')
    for j,s in enumerate(spec):
        mean=sum(F(x[j])*p for x,p in ref.items());check(abs(law.expected_counts()[s]-float(mean))<1e-10,'mean mismatch')
        vals,probs=law.marginal(s);check(abs(sum(probs)-1)<1e-10,'marginal mass')
    draws=law.sample(50,5);check(set(zip(*(draws[s] for s in spec)))<=set(ref),'sample support')
    return err
n=0;maxerr=0
for d in range(1,6):
    for total in range(0,11):
        allowed=[(i,total-i) for i in range(total+1)]
        for kf,kb in [(1,1),(2,3),(5,2)]:
            ratio=(kb/kf)**(1/d);means=[(total+1)*ratio/(1+ratio),(total+1)/(1+ratio)]
            for start in allowed:
                e=compare(f'{d}A <-> {d}B',allowed,start,['A','B'],[((d,0),(0,d),kf),((0,d),(d,0),kb)],means);maxerr=max(e,maxerr);n+=1
extra=0
for total in range(1,7):
    allowed=[(a,b,total-a-b) for a in range(total+1) for b in range(total-a+1)]
    for start in allowed:
        # Equal Poisson means satisfy detailed balance for every reversible pair.
        rr=[((2,0,0),(0,2,0),1),((0,2,0),(2,0,0),1),((0,2,0),(0,0,2),1),((0,0,2),(0,2,0),1)]
        compare(['2A <-> 2B','2B <-> 2C'],allowed,start,['A','B','C'],rr,[1,1,1]);extra+=1
for c in [0,1,2,3]:
    rr=[((1,0,2),(0,1,2),1),((0,1,2),(1,0,2),1)]
    compare('A + 2C <-> B + 2C',[(i,3-i,c) for i in range(4)],(3,0,c),['A','B','C'],rr,[1,1,1]);extra+=1
# Positive fixed reservoir factors must enable structural edges even when their squares underflow.
rx=parse_reactions(['A + 2C <-> B + 2C']);rv={r.rate_key:1. for r in rx};allowed=np.array([(0,2),(1,1),(2,0)])
for chem,expected in [({'C':1e-300},3),({'C':0.0},1)]:
    out=_reachable_finite_class(allowed,np.array([2,0]),rx,['A','B'],rv,chem);check(len(out)==expected,'fixed factor sign control')
# Rate 0 removes an edge; a positive subnormal rate remains structurally active.
rx=parse_reactions(['A <-> B']);rv={rx[0].rate_key:0.,rx[1].rate_key:1.};check(len(_reachable_finite_class(allowed,np.array([2,0]),rx,['A','B'],rv,{}))==1,'zero rate')
rv[rx[0].rate_key]=5e-324;check(len(_reachable_finite_class(allowed,np.array([2,0]),rx,['A','B'],rv,{}))==3,'subnormal rate')
# An incomplete purported finite class is rejected.
try:_reachable_finite_class(allowed[:2],np.array([2,0]),rx,['A','B'],rv,{})
except ValueError:pass
else:raise RuntimeError('incomplete class accepted')
# Constructor positional compatibility: support field must not add an init parameter.
check('_support' not in inspect.signature(StationaryDistribution).parameters,'dataclass constructor changed')
# Legacy plane mixture is stationary; it fails to describe the chosen initial component.
old=[F(1,4),F(1,2),F(1,4)];Q=[[-2,0,2],[0,0,0],[2,0,-2]]
check(all(sum(old[i]*Q[i][j] for i in range(3))==0 for j in range(3)),'legacy negative control')
print(json.dumps({'passed':True,'fresh_exact_two_species_cases':n,'fresh_multispecies_and_catalyst_cases':extra,'probability_comparison_tolerance':1e-10,'all_probability_comparisons_passed':True,'fixed_factors_zero_and_positive_underflow':True,'incomplete_class_rejected':True,'cache_and_constructor_compatibility':True,'legacy_mixture_exact_residual_zero':True},indent=2))
