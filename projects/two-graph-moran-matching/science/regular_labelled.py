"""Full labelled parent-target checks; independent of certificate.py and quotient rates."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path

def require(x,message):
    if not x:raise ArithmeticError(message)

def graph_data(rec):
    N=rec['N'];adj=[set() for _ in range(N)]
    for u,v in rec['edges']:
        require(0<=u<N and 0<=v<N and u!=v,'invalid edge')
        require(v not in adj[u],'duplicate edge');adj[u].add(v);adj[v].add(u)
    require(all(len(x)==rec['degree'] for x in adj),'irregular fixture')
    seen={0};todo=[0]
    while todo:
        for j in adj[todo.pop()]-seen:seen.add(j);todo.append(j)
    require(len(seen)==N,'disconnected fixture')
    return N,rec['degree'],adj

def values(N,d,adj,mask):
    M={i for i in range(N) if mask>>i&1};k=len(M)
    # Count internal edges twice, then divide using integer arithmetic.
    e=sum(j in M for i in M for j in adj[i])//2
    V=F(d*k*(2*N-1+k),2*(N-1))-e
    D=sum((F(len(adj[j]&M))-F(d*k,N-1))**2 for j in range(N) if j not in M)
    return V,D

def events(N,adj,mask,reverse=False):
    for parent in range(N):
        mutant=bool(mask>>parent&1)
        useG=mutant != reverse
        targets=adj[parent] if useG else set(range(N))-{parent}
        for target in sorted(targets):
            child=mask|(1<<target) if mutant else mask&~(1<<target)
            yield parent,target,child,F(1,N*len(targets))

def check_fixture(rec):
    N,d,adj=graph_data(rec);full=(1<<N)-1
    VD=[values(N,d,adj,m) for m in range(1<<N)];V=[x[0] for x in VD];D=[x[1] for x in VD];C=F(d*N*N,N-1);q=N-1-d
    require(V[0]==0 and V[full]==C,'boundary calibration')
    event_count=0;identity_checks=0;complement_checks=0
    for mask in range(1<<N):
        drift=F(0);mass=F(0)
        for parent,target,child,p in events(N,adj,mask):
            drift+=N*p*(V[child]-V[mask]);mass+=p;event_count+=1
        require(mass==1,f'{rec["name"]}: row mass')
        require(drift==-D[mask]/d,f'{rec["name"]} mask {mask}: squaredrift')
        require(drift<=0,'drift sign');identity_checks+=1
    for v in range(N):
        single=1<<v;near=full^single
        require(V[single]==C/N and D[single]==F(d*q,N-1),'singleton')
        require(V[near]/C==F(N-1,N) and D[near]==0,'nearfix')
        expectedD=sum(p*D[child] for _,_,child,p in events(N,adj,near))
        require(expectedD==F(2*d*q,N*(N-1)**2),'nearfix nextD')
        # Verify color-complement relation at the event level, not via shared rates.
        lhs={(par,targ,full^child,p) for par,targ,child,p in events(N,adj,single,True)}
        rhs=set(events(N,adj,near,False))
        require(lhs==rhs,'color-complement kernel');complement_checks+=1
    return {'name':rec['name'],'N':N,'d':d,'states':1<<N,'parent_target_events':event_count,'identity_checks':identity_checks,'color_complement_kernel_checks':complement_checks,'complete_control':q==0,'singleton_gap':str(F(q,d*N**3)),'reverse_gap':str(F(2*q,d*N**4*(N-1)))}

def run(fixture):
    records=json.loads(fixture.read_text())['graphs'];results=[check_fixture(x) for x in records]
    return {'schema':'regular-opponent-labelled-v1','graphs':results,'totals':{k:sum(x[k] for x in results) for k in ['states','parent_target_events','identity_checks','color_complement_kernel_checks']},'all_checks_passed':True,'scope':'finite raw labelled kernels only; uniform theorem is in proof.md; no fixation linear-system solve used'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--fixtures',type=Path,default=Path(__file__).parent/'evidence/graphs.json');p.add_argument('--output',type=Path);a=p.parse_args();d=run(a.fixtures);s=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if a.output:
        with a.output.open('x') as f:f.write(s)
    else:print(s,end='')
