"""Literal N=3 three-site Fleming-Viot count and founder systems."""
import sympy as s
from sympy.polys.matrices import DomainMatrix


def counts(n):
    return [(a,b,n-a-b) for a in range(n+1) for b in range(n-a+1)]


def moved(c,a,b):
    out=list(c);out[a]-=1;out[b]+=1
    return tuple(out)


def matrices(k):
    pop=counts(3); oth=counts(2)
    states=[(i,d) for i in range(3) for d in oth]
    ci={c:j for j,c in enumerate(pop)};si={x:j for j,x in enumerate(states)}
    L=s.zeros(10);T=s.zeros(18)
    for c in pop:
        row=ci[c]
        for i in range(3):
            if not c[i]:continue
            for j in range(3):
                if i==j:continue
                rate=c[i]+s.Rational(1,2)*c[i]*c[j]*k[i]
                L[row,ci[moved(c,i,j)]]+=rate;L[row,row]-=rate
    for i,d in states:
        row=si[i,d]
        for j in range(3):
            if i==j:continue
            T[row,si[j,d]]+=1;T[row,row]-=1
        for a in range(3):
            if not d[a]:continue
            for b in range(3):
                if a==b:continue
                T[row,si[i,moved(d,a,b)]]+=d[a];T[row,row]-=d[a]
        T[row,row]-=k[i]
        for a in range(3):
            if not d[a]:continue
            rate=d[a]*k[a]
            T[row,si[i,moved(d,a,i)]]+=rate
            for b in range(3):
                donor=d[b]-int(a==b)
                if donor:
                    T[row,si[i,moved(d,a,b)]]+=rate*s.Rational(donor,2)
            T[row,row]-=rate
    return pop,states,ci,si,L,T


def normalized_systems(k):
    pop,states,ci,si,L,T=matrices(k)
    P=L.T.copy();P[0,:]=s.ones(1,10);pr=s.zeros(10,1);pr[0]=1
    H=T.copy();anchor=si[0,(2,0,0)];H[anchor,:]=s.zeros(1,18);H[anchor,anchor]=1
    hr=s.zeros(18,1);hr[anchor]=s.Rational(1,3)
    return pop,states,ci,si,L,T,P,pr,H,hr


def exact_solve(A,b):
    dm=DomainMatrix.from_Matrix(A);db=DomainMatrix.from_Matrix(b).convert_to(dm.domain)
    numerator,denominator=dm.solve_den(db,method="rref")
    return numerator.to_Matrix(),dm.domain.to_sympy(denominator)


def current(k):
    pop,states,ci,si,L,T,P,pr,H,hr=normalized_systems(k)
    pn,pd=exact_solve(P,pr);hn,hd=exact_solve(H*3,hr*3)
    F=s.zeros(3)
    for i in range(3):
        for j in range(3):
            for c in pop:
                if c[i]:
                    d=list(c);d[i]-=1
                    F[i,j]+=pn[ci[c]]*c[i]*hn[si[j,tuple(d)]]
    return s.cancel((F[0,2]-F[2,0])/(pd*hd)),(pn,pd,hn,hd,F)
