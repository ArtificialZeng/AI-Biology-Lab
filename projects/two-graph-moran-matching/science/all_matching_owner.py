"""Exact uniform drift certificate. Writes only this dossier's evidence."""
import json
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parent
n,m,a,b,c,d,r=s.symbols('n m a b c d r'); A,B,E=s.symbols('A B E')
k=2*a+b+c;T=(2*a+b)/(n-2)+c/(n-1);L=(n-k)/(n-1)
rates=[((a,b+1,c),2*(m-a-b)*T),((a+1,b-1,c),b*(T-1/(n-2))),((a-1,b+1,c),2*a*L),((a,b-1,c),b*L),((a,b,c+1),(n-2*m-c)*T),((a,b,c-1),c*L)]
W=(n-2+(m-2)/n)*k+a+m*k*(n-k)/(n*(n-1))
q=s.factor(sum(rate*(W.subs(dict(zip((a,b,c),target)),simultaneous=True)-W) for target,rate in rates))
D=n*(n-2)*(n-1)**2
P=s.Poly(s.expand(s.cancel(-D*q).subs(n,2*m+c+r).subs(m,a+b+d)),a,b,d,c,r)
a2=a*(a-1);a3=a2*(a-2);b2=b*(b-1);b3=b2*(b-2);b4=b3*(b-3);d2=d*(d-1);d3=d2*(d-2)
coeffs={
 (2,1):2*(a+b+d),
 (1,2):2*(a+b+d),
 (2,0):4*a*d+b2+6*b*d+4*d*d,
 (0,2):4*a+b*b+2*b*d+b,
 (1,1):2*(4*a*a+7*a*b+8*a*d+2*b*b+2*b2+9*b*d+2*d*d+2*d2),
 (1,0):16*a*a*d+4*a*b2+36*a*b*d+24*a*d*d+4*b3+5*b2+24*b2*d+28*b*d2+38*b*d+8*d3+16*d2,
 (0,1):8*a3+28*a2+8*a+12*a2*b+8*a*b2+24*a*b+8*a2*d+20*a*d+12*a*b*d+4*b3+9*b2+2*b+12*b2*d+8*b*d2+18*b*d,
 (0,0):d*(16*a3+72*a2+40*a)+b*d*(32*a2+84*a+24)+d2*(16*a2+24*a)+b2*(2+48*d)+b3*(6+2*a+14*d)+2*b4+28*a*b2*d+b*d2*(36+24*a)+20*b2*d2+8*b*d3,
}
assert s.expand(P.as_expr()-sum(co*c**i*r**j for (i,j),co in coeffs.items()))==0
cases=[]
for av,alabel in [(0,'a=0'),(1+A,'a=1+A')]:
 for bv,blabel in [(0,'b=0'),(1,'b=1'),(2+B,'b=2+B')]:
  for dv,dlabel in [(0,'d=0'),(1+E,'d=1+E')]:
   pol=s.Poly(s.expand(P.as_expr().subs({a:av,b:bv,d:dv})),A,B,E,c,r)
   assert all(co>=0 for co in pol.coeffs())
   cases.append({'case':[alabel,blabel,dlabel],'variables':['A','B','E','c','r'],'terms':[{'powers':list(pow),'coefficient':int(co)} for pow,co in pol.terms() if co]})
assert s.factor(q.subs({a:0,b:1,c:0}))==-2*m/(n*(n-1))
assert s.factor(q.subs({a:0,b:0,c:1}))==-2*m/(n*(n-1))
C=n*(n-2)+2*m-2
assert s.factor(W.subs({a:m,b:0,c:n-2*m})-C)==0
assert s.factor(W.subs({a:0,b:1,c:0})-C/n)==0
assert s.factor(W.subs({a:0,b:0,c:1})-C/n)==0
out={'schema':'exact-uniform-drift-certificate-v1','sympy':s.__version__,'variables':['a','b','d','c','r'],'P_terms':[{'powers':list(pow),'coefficient':int(co)} for pow,co in P.terms()], 'P_expression':str(P.as_expr()),'c_r_coefficient_identities':{f'{i},{j}':str(co) for (i,j),co in coeffs.items()},'nonnegative_integer_basis':'x_j=x(x-1)...(x-j+1), j=2,3,4; every displayed coefficient is a sum of products with nonnegative numerical coefficients','exhaustive_shift_cases':cases,'assertions':['exact generator identity','8 coefficient identities','12 exhaustive shifted coefficient certificates','both singleton drifts','fixation and singleton normalization'],'scope':'integers N>=3,1<=m<=floor(N/2),a,b,d,c,r>=0,m=a+b+d,N=2m+c+r'}
(ROOT/'evidence'/'uniform-certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':'exact symbolic identities passed','P_terms':len(P.terms()),'shift_cases':len(cases),'scope':out['scope']},indent=2))
