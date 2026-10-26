#!/usr/bin/env python3
"""Exact full-model reconstruction and the population proof identities; stdout only."""
import json
import sympy as S
import sys
sys.dont_write_bytecode = True
from archive_tools import run


def require(test,why):
 if not test: raise ValueError(why)
def equal(a,b,why):
 require(S.cancel(a-b)==0,why)

def build_certificate():
 require(S.__version__ == "1.13.3", "Install pinned SymPy 1.13.3")
 s,i1,i21,i2,i12,r1,r2,r12=S.symbols('s i1 i21 i2 i12 r1 r2 r12')
 state=[s,i1,i21,i2,i12,r1,r2,r12]
 L,mu,g1,g2,b1,b2,e1,e2,h1,h2,t1,t2,t12=S.symbols('Lambda mu gamma1 gamma2 beta1 beta2 eta1 eta2 sigma1 sigma2 theta1 theta2 theta12',positive=True)
 v1,v2=g1+mu,g2+mu
 force1=i1+e1*i21;force2=i2+e2*i12
 f=S.Matrix([L-mu*s-b1*force1*s-b2*force2*s+t1*r1+t2*r2+t12*r12,
 b1*force1*s-v1*i1,b1*h1*force1*r2-v1*i21,
 b2*force2*s-v2*i2,b2*h2*force2*r1-v2*i12,
 g1*i1-(b2*h2*force2+mu+t1)*r1,
 g2*i2-(b1*h1*force1+mu+t2)*r2,
 g1*i21+g2*i12-(mu+t12)*r12])
 # Reconstruct every ODE term independently from standard mass-action reactions.
 reactions=[]
 def reaction(u,v,k):
  U=tuple(u.count(j) for j in range(8));V=tuple(v.count(j) for j in range(8))
  reactions.append((U,V,k))
 reaction([], [0],L)
 for j in range(8): reaction([j],[],mu)
 for sus,inf,new,k in [(0,1,1,b1),(0,2,1,b1*e1),(0,3,3,b2),(0,4,3,b2*e2),
                       (6,1,2,b1*h1),(6,2,2,b1*h1*e1),(5,3,4,b2*h2),(5,4,4,b2*h2*e2)]:
  reaction([sus,inf],[inf,new],k)
 for a,z,k in [(1,5,g1),(2,7,g1),(3,6,g2),(4,7,g2),(5,0,t1),(6,0,t2),(7,0,t12)]:
  reaction([a],[z],k)
 direct=S.zeros(8,1)
 for u,v,k in reactions:
  rate=k*S.prod(state[j]**u[j] for j in range(8))
  direct += S.Matrix([v[j]-u[j] for j in range(8)])*rate
 for a,b in zip(direct,f):equal(a,b,'Reaction reconstruction failed')
 equal(sum(f),L-mu*sum(state),'Full population identity failed')
 # Use free r1,r2 and rank-one degeneracy to avoid imposing arbitrary parameters.
 d1,d2,d12=mu+t1,mu+t2,mu+t12
 D=g1*g2*s*s-h1*h2*v1*v2*r1*r2
 A1=r1*(g2*s*d1+h2*v2*d2*r2)
 A2=r2*(g1*s*d2+h1*v1*d1*r1)
 c1,c2=A1/D,A2/D
 recovered12=(g1*h1*r2*c1+g2*h2*r1*c2)/d12
 reconstruction={i1:s*c1,i21:h1*r2*c1,i2:s*c2,i12:h2*r1*c2,r12:recovered12,
                 b1:v1/(s+e1*h1*r2),b2:v2/(s+e2*h2*r1)}
 reconstructed=[S.cancel(z.subs(reconstruction,simultaneous=True)) for z in f]
 for j in range(1,8):equal(reconstructed[j],0,'Non-susceptible balance failed: '+str(j))
 N=s+r1+r2+(s+h1*r2)*c1+(s+h2*r1)*c2+recovered12
 equal(reconstructed[0],L-mu*N,'Reduced susceptible equation is not Lambda-mu*N')
 W1=s+h1*r2*(1+g1/d12);W2=s+h2*r1*(1+g2/d12)
 P=mu*((s+r1+r2-L/mu)*D+A1*W1+A2*W2)
 equal(-D*reconstructed[0],P,'Denominator-cleared population polynomial failed')
 # Ordered right endpoint p<q; P(p) factors as minus a positive times invasion deficit.
 p,q=S.symbols('p q',positive=True)
 R1=(q-s)/(e2*h2);R2=(p-s)/(e1*h1)
 endpoint={s:p,r1:(q-p)/(e2*h2),r2:0}
 inv=(p+e2*h2*g1*(L/mu-p)/(v1+t1))/q
 equal(P.subs(endpoint,simultaneous=True),
       -g2*p*p*mu*(v1+t1)*q/(e2*h2)*(inv-1),'Ordered endpoint invasion identity failed')
 # Positive population lower-bound identity for the converse.
 # At coexistence r1 balance, N exceeds s+(v1+theta1)*r1/gamma1.
 excess=sum(state)-s-(v1+t1)/g1*r1
 expected=i21+i2+i12+r2+r12+b2*h2*force2*r1/g1
 equal(excess-expected,f[5]/g1,'Population lower-bound balance identity failed')
 # A rational witness to why the lower DFE threshold cannot be omitted.
 params={L:S.Rational(23,18),mu:1,g1:1,g2:1,b1:2,b2:1,e1:30,e2:10,h1:S.Rational(1,2),h2:S.Rational(1,2),t1:0,t2:0,t12:0}
 point={s:S.Rational(1,2),i1:S.Rational(1,3),i21:S.Rational(1,90),i2:S.Rational(1,18),i12:S.Rational(1,60),r1:S.Rational(3,10),r2:S.Rational(1,30),r12:S.Rational(1,36)}
 for z in point.values():require(z>0,'Witness coordinate nonpositive')
 for z in f:equal(z.subs(params).subs(point),0,'Rational equilibrium witness failed')
 inv_w=S.Rational(1,2)*(1+5*S.Rational(5,36))
 equal(inv_w,S.Rational(61,72),'Witness invasion value differs')
 det=S.factor(f.jacobian(state).subs(params).subs(point).det())
 require(det!=0,'Rational witness is degenerate')
 return {'verified':True,'arithmetic':'exact rational and symbolic identities','reaction_count':len(reactions),
 'full_population_identity':'sum(f)=Lambda-mu*N','seven_reconstructed_balances':'all exactly zero',
 'reduced_population_identity':'F_G(s)=Lambda-mu*N(s)',
 'denominator_cleared_identity':'P_G=mu*((s+r1+r2-Lambda/mu)*Delta+A1*W1+A2*W2)',
 'endpoint_identity':'P_G(p)=-gamma2*p^2*mu*(v1+theta1)*q/(eta2*sigma2)*(R2(E1)-1)',
 'converse_bound':'N>s+(v1+theta1)*r1/gamma1 at every strictly positive coexistence equilibrium',
 'witness':{'parameters':{str(k):str(v) for k,v in params.items()},'state':{str(k):str(v) for k,v in point.items()},'R1_DFE':'23/18','R2_DFE':'23/36','R2_E1':str(inv_w),'Jacobian_determinant':str(det)},
 'limitations':['Existence and conditional converse only; no uniqueness, stability or full Question2 root-count classification.','The general A1-A5 Conjecture1 remains unresolved.']}

if __name__ == '__main__':
 raise SystemExit(run(build_certificate, 'evidence/exact-population-certificate.json', 'original symbolic reconstruction'))
