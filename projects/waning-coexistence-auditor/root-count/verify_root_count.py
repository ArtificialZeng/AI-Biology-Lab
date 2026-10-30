#!/usr/bin/env python3
"""Exact full-model reconstruction and the population proof identities; stdout only."""
import json
import argparse
import hashlib
from pathlib import Path
import sympy as S

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

 # New pass2 identities: exact source PG, strict-convexity factors, and cubic transformation.
 curve={r1:R1,r2:R2}
 Nc=N.subs(curve,simultaneous=True)
 Pc=P.subs(curve,simultaneous=True)
 degree=S.Poly(Pc,s).degree()
 require(degree<=3,'Source PG exceeds cubic degree')
 K=h1*h2*v1*v2*R1*R2/(g1*g2*s*s)
 H=1/(1-K);T=R1*R2/s
 C=v1*v2/(g1*g2*e1*e2)
 equal(K,C*(p/s-1)*(q/s-1),'K factorization')
 equal(S.diff(K,s),C*((p+q)*s-2*p*q)/s**3,'K derivative')
 equal(S.diff(K,s,2),2*C*(3*p*q-(p+q)*s)/s**4,'K second derivative')
 equal(T,(p*q/s-p-q+s)/(e1*e2*h1*h2),'T factorization')
 equal((s*c1).subs(curve,simultaneous=True),H*(d1/g1*R1+h2*v2*d2/(g1*g2)*T),'Primary1 convex decomposition')
 equal((s*c2).subs(curve,simultaneous=True),H*(d2/g2*R2+h1*v1*d1/(g1*g2)*T),'Primary2 convex decomposition')
 Np=p+(v1+t1)*(q-p)/(g1*e2*h2)
 equal(Nc.subs(s,p),Np,'General endpoint population')
 tie_N=S.cancel(Nc.subs(q,p))
 tie_derivative=1-(v1+t1)/(g1*e2*h2)-(v2+t2)/(g2*e1*h1)
 equal(tie_N.subs(s,p),p,'Coincident endpoint population')
 equal(S.diff(tie_N,s).subs(s,p),tie_derivative,'Coincident endpoint derivative')
 # Universal cubic identity includes lower original degrees.
 u0,u1,u2,u3,z,w=S.symbols('u0 u1 u2 u3 z w')
 ell=p-z;G=u0+u1*s+u2*s*s+u3*s**3
 Gq=S.Poly(S.cancel((1+w)**3*G.subs(s,(z+p*w)/(1+w))),w)
 aa=G.subs(s,p);bb=3*aa-ell*S.diff(G,s).subs(s,p)
 zz=G.subs(s,z);cc=3*zz+ell*S.diff(G,s).subs(s,z)
 equal(Gq.as_expr(),aa*w**3+bb*w*w+cc*w+zz,'Mobius cubic coefficient identity')
 a,b,c,z0=S.symbols('a b c z0')
 discr=b*b*c*c-4*a*c**3-4*b**3*z0-27*a*a*z0*z0+18*a*b*c*z0
 equal(S.discriminant(a*w**3+b*w*w+c*w+z0,w),discr,'Cubic discriminant formula')
 coordinates=[s,s*c1,h1*r2*c1,s*c2,h2*r1*c2,r1,r2,recovered12]
 coordcurve=[S.cancel(x.subs(curve,simultaneous=True)) for x in coordinates]
 def example(name,parameters,roots,expected):
  pv=S.cancel((v1/b1).subs(parameters));qv=S.cancel((v2/b2).subs(parameters))
  values={**parameters,p:pv,q:qv}
  actual=S.factor(Pc.subs(values))
  equal(actual,expected,'Example PG factorization: '+name)
  rows=[]
  for root,multiplicity in roots:
   point_values=[S.cancel(x.subs(values).subs(s,root)) for x in coordcurve]
   require(all(x>0 for x in point_values),'Example coordinate positivity: '+name)
   point_map=dict(zip(state,point_values))
   for balance in f: equal(balance.subs(parameters).subs(point_map),0,'Eight-equation example balance: '+name)
   equal(sum(point_values),parameters[L]/parameters[mu],'Example population: '+name)
   for order in range(multiplicity):equal(S.diff(actual,s,order).subs(s,root),0,'Example root multiplicity')
   require(S.diff(actual,s,multiplicity).subs(s,root)!=0,'Example multiplicity not exact')
   rows.append({'s':str(root),'multiplicity':multiplicity,'state':{str(k):str(v) for k,v in point_map.items()}})
  return {'parameters':{str(k):str(v) for k,v in parameters.items()},'PG':str(actual),'degree':int(S.degree(actual,s)),'roots':rows}
 common={mu:S.Integer(1),g1:S.Integer(1),g2:S.Integer(1),b1:S.Integer(2),b2:S.Integer(1),e1:S.Integer(20),e2:S.Integer(10),h1:S.Integer(1),h2:S.Integer(1),t2:S.Integer(1),t12:S.Integer(1)}
 pair={**common,t1:S.Rational(535,216),L:S.Rational(2011,1440)}
 fold={**common,t1:S.Rational(535,262),L:S.Rational(872,655)}
 tie={**common,b2:S.Integer(2),e1:S.Integer(10),t1:S.Integer(1),L:S.Integer(1)}
 examples={
 'all_positive_waning_two_roots':example('pair',pair,[(S.Rational(2,3),1),(S.Rational(3,4),1)],(3*s-2)*(4*s-3)*(854*s+251)/21600),
 'all_positive_waning_double_root':example('fold',fold,[(S.Rational(2,3),2)],(3*s-2)**2*(595*s+169)/10480),
 'coincident_critical_threshold':example('tie',tie,[(S.Rational(11,31),1)],(s-1)*(4*s+1)*(31*s-11)/250)}
 foldvalues={**fold,p:S.Integer(1),q:S.Integer(2)}
 foldN=S.cancel(Nc.subs(foldvalues))
 foldsecond=S.cancel(S.diff(foldN,s,2).subs(s,S.Rational(2,3)))
 equal(foldsecond,S.Rational(229095,102704),'Fold second derivative')
 equal(foldN.subs(s,1)-fold[L],S.Rational(191,2620),'Fold endpoint gap')
 equal(S.diff(Pc,L),-D.subs(curve,simultaneous=True),'Population-level transversality')
 # Exact original degree collapse without an excluded parameter limit.
 lower={**common,e1:S.Integer(1),e2:S.Integer(2),t1:S.Integer(1),p:S.Integer(1),q:S.Integer(2)}
 degree2=S.factor(Pc.subs({**lower,L:S.Integer(3)}))
 degree1=S.factor(Pc.subs({**lower,L:S.Rational(5,2)}))
 equal(degree2,(s*s-34*s+32)/2,'Degree2 boundary family')
 equal(degree1,14*(1-s),'Degree1 boundary family')
 return {'schema':'waning-root-count-v1','verified':True,'sympy':S.__version__,
 'arithmetic':'exact symbolic and rational arithmetic','reaction_count':len(reactions),
 'generic_checks':['reaction and source ODE reconstruction','population identity','all seven reconstructed balances','source PG=-Delta*FG exactly','degree at most three','convex-factor identities','coincident endpoint identities','Mobius coefficient identities','cubic discriminant'],
 'source_PG_degree_bound':int(degree),'classification':{
 'a<0':'one simple','a=0 and b<0':'one simple','a=0 and b>=0':'none',
 'a>0 and (b<0 or c<0) and discriminant>0':'two simple',
 'a>0 and (b<0 or c<0) and discriminant=0':'one double','other a>0 cases':'none'},
 'examples':examples,'fold_population_second_derivative':str(foldsecond),
 'positive_parameter_degree_collapse':{'degree2':str(degree2),'degree1':str(degree1)},
 'limits':['Sign/convexity/real-root arguments are analytic proof, not formally encoded by symbolic identities.','Broad abstract Conjecture1 remains unresolved.','Scalar multiplicity makes no ODE stability, persistence or dynamical saddle-node claim.','Fresh review of this pass is pending; no formal runner identity or human acceptance.']}

def encode(obj):
 return (json.dumps(obj,indent=2,sort_keys=True)+'\n').encode()

if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('mode',choices=['compute','verify'])
 args=parser.parse_args()
 try:
  computed=encode(build_certificate())
  if args.mode=='compute': print(computed.decode(),end='')
  else:
   archive=Path(__file__).resolve().parent/'evidence/exact-root-certificate.json'
   require(computed==archive.read_bytes(),'Exact certificate differs from archive')
   print(json.dumps({'ok':True,'archive_modified':False,'certificate_sha256':hashlib.sha256(computed).hexdigest()},sort_keys=True))
 except Exception as exc:
  print(json.dumps({'ok':False,'error':str(exc)},sort_keys=True));raise SystemExit(1)
