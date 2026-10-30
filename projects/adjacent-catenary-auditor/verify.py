#!/usr/bin/env python3
"""Exact symbolic coefficient-Jacobian builder; no archive writes."""
import sys
sys.dont_write_bytecode = True
import sympy as s
from archive_tools import require, run

z=s.symbols('z')

def model(n,ell):
 p=s.symbols('p1:'+str(n));q=s.symbols('q1:'+str(n));e=s.symbols('e')
 A=s.zeros(n)
 for i in range(n-1):
  A[i+1,i]=p[i];A[i,i+1]=q[i]
  A[i,i]-=p[i];A[i+1,i+1]-=q[i]
 A[ell-1,ell-1]-=e
 M=z*s.eye(n)-A
 D=M.det(method='domain-ge');N=s.expand(p[0]*M[2:,2:].det(method='domain-ge'))
 c=s.Matrix(s.Poly(D,z).all_coeffs()[1:]+s.Poly(N,z).all_coeffs())
 return p,q,e,A,D,N,c

def build_certificate():
 require(s.__version__ == '1.13.3', 'Install pinned SymPy 1.13.3')
 out=[]
 for n in range(2,5):
  for ell in range(1,n+1):
   p,q,e,A,D,N,c=model(n,ell)
   pars=p+q+(e,);J=c.jacobian(pars)
   det=s.factor(J.det(method='domain-ge'))
   r={'n':n,'ell':ell,'parameters':[str(x) for x in pars],'A':str(A),'D':str(s.Poly(D,z).as_expr()),'N':str(s.factor(N)),'coefficient_map':[str(x) for x in c],'jacobian_determinant':str(det)}
   out.append(r)
 return out

if __name__ == '__main__':
 raise SystemExit(run(build_certificate, 'evidence/small-path-jacobians.json', 'direct characteristic determinant'))
