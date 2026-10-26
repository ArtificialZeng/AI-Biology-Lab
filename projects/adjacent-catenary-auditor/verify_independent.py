#!/usr/bin/env python3
"""Independent exact rank builder: direct cofactor and suffix recurrence."""
import sys
sys.dont_write_bytecode = True
from itertools import product
import sympy as S
from archive_tools import require, run

z=S.Symbol('z')

def build(n,ell):
 p=S.symbols(f'p1:{n}');q=S.symbols(f'q1:{n}');e=S.Symbol('e');pars=p+q+(e,)
 M=S.zeros(n)
 for j in range(n):
  M[j,j]=z+(p[j] if j<n-1 else 0)+(q[j-1] if j>0 else 0)+(e if ell==j+1 else 0)
 for j in range(n-1):M[j+1,j]=-p[j];M[j,j+1]=-q[j]
 def suffix(j):
  f1=S.Integer(1);f2=S.Integer(0)
  for k in range(n-1,j-1,-1):
   f1,f2=S.expand(M[k,k]*f1-(p[k]*q[k]*f2 if k<n-1 else 0)),f1
  return f1
 D=suffix(0);B=suffix(1);R=suffix(2)
 N=S.expand((-1)**3*M.minor_submatrix(0,1).det(method='domain-ge'))
 require(S.expand(N-p[0]*R)==0, 'Direct cofactor and tail determinant disagree')
 c=S.Matrix(S.Poly(D,z).all_coeffs()[1:]+S.Poly(N,z).all_coeffs())
 J=c.jacobian(pars)
 return pars,p,q,e,M,D,B,R,N,J

def inspect(n,ell,vals,label):
 pars,p,q,e,M,D,B,R,N,J=build(n,ell)
 sub=dict(zip(pars,vals));a=p[0]+(e if ell==1 else 0)
 C=S.factor(R.subs(z,-a));T=S.cancel((B-q[0]*R)/z) if ell==1 else None
 K=S.expand(q[0]*R.subs(z,-a)-e*T.subs(z,-a)) if ell==1 else C
 det=S.factor(J.det(method='domain-ge')) if n<=4 else None
 rj=J.subs(sub).rank()
 # No division/cancellation is used in the rational differential.  Its exact
 # coefficient matrix is equivalent to delta H after multiplying by D^2.
 Dv=D.subs(sub);Nv=N.subs(sub)
 cols=[]
 for t in pars:
  poly=S.Poly(S.expand(Dv*S.diff(N,t).subs(sub)-Nv*S.diff(D,t).subs(sub)),z)
  cols.append([poly.nth(k) for k in range(2*n)])
 rh=S.Matrix.hstack(*map(S.Matrix,cols)).rank()
 kv=K.subs(sub);cv=C.subs(sub)
 expected_c=2*n-2 if kv==0 else 2*n-1
 expected_h=2*n-2 if kv==0 or cv==0 else 2*n-1
 require(rj==expected_c, 'Coefficient rank mismatch: '+str((n,ell,vals,rj,kv)))
 require(rh==expected_h, 'Rational transfer rank mismatch: '+str((n,ell,vals,rh,kv,cv)))
 gcd=S.gcd(S.Poly(Dv,z),S.Poly(Nv,z)).monic().as_expr()
 require(S.degree(gcd,z)==(1 if cv==0 else 0), 'Common-factor degree mismatch')
 return {'n':n,'ell':ell,'label':label,'parameters':[str(v) for v in vals],'K_or_R_value':str(kv),'cancellation_value':str(cv),'coefficient_rank':rj,'rational_transfer_differential_rank':rh,'gcd':str(gcd),'determinant_at_point':str(det.subs(sub)) if det is not None else None}

def build_certificate():
 require(S.__version__ == '1.13.3', 'Install pinned SymPy 1.13.3')
 results=[]
 witnesses={1:[1,1,1,1,1,1,1],2:[1,1,1,1,2,2,1],3:[1,1,1,1,1,2,1],4:[1,1,1,1,1,1,1]}
 for ell in range(1,5):
  results.append(inspect(4,ell,witnesses[ell],'positive coefficient-singular witness'))
  results.append(inspect(4,ell,[2,3,5,7,11,13,17],'fixed generic check'))
 results.append(inspect(4,1,[1,1,S.Rational(1,2),1,3,3,1],'cancellation with full coefficient rank'))
 # Held-out length5: independently selected exact positive points and singular
 # witnesses found by a bounded integer grid in tail rates.  This is finite
 # evidence supporting, not proving, the all-n theorem.
 for ell in range(1,6):
  results.append(inspect(5,ell,[2,3,5,7,11,13,17,19,23],'held-out fixed generic check'))
  pars,p,q,e,M,D,B,R,N,J=build(5,ell)
  found=None
  # Vary at most three small parameters; all other rates=1.
  for pp in product(range(1,5),repeat=3):
   vals=[1]*9;vals[0]=pp[0];vals[2]=pp[1];vals[7]=pp[2]
   sub=dict(zip(pars,vals));a=p[0]+(e if ell==1 else 0)
   if ell==1:
    # q1 occurs only in B and has been canceled from T.
    T=S.cancel((B-q[0]*R)/z);rv=R.subs(z,-a).subs(sub);tv=T.subs(z,-a).subs(sub)
    if rv and tv/rv>0:
     vals[4]=S.cancel(tv/rv);found=vals;break
   elif R.subs(z,-a).subs(sub)==0:
    found=vals;break
  if found is not None:results.append(inspect(5,ell,found,'held-out positive singular witness'))
 return {'arithmetic':'exact SymPy rationals','coverage':'n4 all leak positions, positive singular/generic witnesses and ell1 cancellation; n5 all positions generic and bounded witness search','results':results}

if __name__ == '__main__':
 raise SystemExit(run(build_certificate, 'evidence/independent-checks.json', 'direct cofactor and suffix recurrence'))
