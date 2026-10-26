"""Symbolic all-N certificate, independent of finite harmonic solves."""
from pathlib import Path
import json
import sympy as s
N,a,b,c,r=s.symbols('N a b c r')
k=2*a+b+c
T=(2*a+b)/(N-2)+c/(N-1)
L=(N-k)/(N-1)
V=(N-2)*k+a+2*k*(N-k)/(N*(N-1))
transition=[((a,b+1,c),2*(2-a-b)*T),((a+1,b-1,c),b*(T-1/(N-2))),((a-1,b+1,c),2*a*L),((a,b-1,c),b*L),((a,b,c+1),(N-4-c)*T),((a,b,c-1),c*L)]
QV=s.factor(sum(rate*(V.subs({a:t[0],b:t[1],c:t[2]},simultaneous=True)-V)for t,rate in transition))
D=N*(N-2)*(N-1)**2
certificates={
 (0,0):4*c*(r+4)*(c+r+2),
 (0,1):2*(c+r+3)*(2*c*r+5*c+2*r+4),
 (0,2):2*(2*c**2*r+c**2+2*c*r**2+12*c*r+5*c+3*r**2+11*r+2),
 (1,0):4*(r+2)*(c**2+c*r+5*c+r+5),
 (1,1):2*r*(2*c**2+2*c*r+13*c+3*r+17),
 (2,0):4*r*(c**2+c*r+8*c+2*r+18),
}
records=[]
for (aa,bb),expected in certificates.items():
 actual=s.factor((-D*QV).subs({a:aa,b:bb,N:c+r+4}))
 assert s.expand(actual-expected)==0,(aa,bb,actual)
 pp=s.Poly(s.expand(expected),c,r)
 assert all(z>=0 for z in pp.coeffs())
 records.append({'a':aa,'b':bb,'factorized_polynomial':str(s.factor(expected)),'terms':[{'powers':list(m),'coefficient':int(v)}for m,v in pp.terms()]})
assert len(records)==6
C=N*(N-2)+2
assert s.simplify(V.subs({a:0,b:0,c:0}))==0
assert s.simplify(V.subs({a:2,b:0,c:N-4})-C)==0
for initial in [{a:0,b:1,c:0},{a:0,b:0,c:1}]:
 assert s.factor(V.subs(initial)-C/N)==0
 assert s.factor(QV.subs(initial)+4/(N*(N-1)))==0
out={'statement':'m=2,N>=4,all states c,r=N-4-c nonnegative; Q=N(P-I)','symbolic_identities':records,'singleton_actual_drift':'-4/[N^2(N-1)]','singleton_upper_bound':'1/N-4/[N^2(N-1)(N(N-2)+2)]','sympy_version':s.__version__}
p=Path(__file__).parent/'evidence';p.mkdir(exist_ok=True)
(p/'six-polynomial-certificates.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'verified_symbolic_cases':6,'nonnegative_integer_coefficient_certificates':True,'endpoints_and_singleton_identity':True,'proof_domain':'all integer N>=4,not bounded enumeration'}))
