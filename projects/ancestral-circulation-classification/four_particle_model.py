"""General finite triangle count/founder systems from source Eq3.2/3.6."""
import sympy as s
from sympy.polys.matrices import DomainMatrix

def require(ok,message):
    if not ok: raise RuntimeError(message)
def counts(n): return [(a,b,n-a-b) for a in range(n+1) for b in range(n-a+1)]
def moved(c,a,b):
    out=list(c);out[a]-=1;out[b]+=1;return tuple(out)
def systems(N,k):
    require(N>=2 and len(k)==3,'model parameters')
    pop=counts(N);states=[(i,d) for i in range(3) for d in counts(N-1)]
    ci={c:j for j,c in enumerate(pop)};si={x:j for j,x in enumerate(states)}
    L=s.zeros(len(pop));T=s.zeros(len(states));den=N-1
    for c in pop:
        row=ci[c]
        for i in range(3):
            if c[i]:
                for j in range(3):
                    if i!=j:
                        rate=c[i]+s.Rational(c[i]*c[j],den)*k[i]
                        L[row,ci[moved(c,i,j)]]+=rate;L[row,row]-=rate
    for i,d in states:
        row=si[i,d]
        for j in range(3):
            if i!=j:T[row,si[j,d]]+=1;T[row,row]-=1
        for a in range(3):
            if d[a]:
                for b in range(3):
                    if a!=b:T[row,si[i,moved(d,a,b)]]+=d[a];T[row,row]-=d[a]
        T[row,row]-=k[i]
        for a in range(3):
            if d[a]:
                rate=d[a]*k[a];T[row,si[i,moved(d,a,i)]]+=rate*s.Rational(2,den)
                for b in range(3):
                    donor=d[b]-int(a==b)
                    if donor:T[row,si[i,moved(d,a,b)]]+=rate*s.Rational(donor,den)
                T[row,row]-=rate
    P=L.T.copy();P[0,:]=s.ones(1,len(pop));pr=s.zeros(len(pop),1);pr[0]=1
    H=T.copy();anchor=si[0,(N-1,0,0)];H[anchor,:]=s.zeros(1,len(states));H[anchor,anchor]=1
    hr=s.zeros(len(states),1);hr[anchor]=s.Rational(1,N)
    return pop,states,ci,si,L,T,P,pr,H,hr

def exact_solve(M,b):
    dm=DomainMatrix.from_Matrix(M);db=DomainMatrix.from_Matrix(b).convert_to(dm.domain)
    num,den=dm.solve_den(db,method='rref')
    return num.to_Matrix(),dm.domain.to_sympy(den)
