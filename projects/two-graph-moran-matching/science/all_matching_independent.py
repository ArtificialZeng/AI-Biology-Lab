"""Fresh referee checker: event-count derivation, exact algebra, held-out labels.

No import or execution of the researcher's generator/certificate/checker.
Only the frozen hash maps and (after independent derivation) stored P terms are read.
"""
import hashlib, json, os, sys
from pathlib import Path
import sympy as S

ROOT = Path(__file__).resolve().parent
RESEARCH = ROOT



n,m,a,b,c,d,r = S.symbols('n m a b c d r')
k=2*a+b+c
# Q=N(P-I). Count events by parent type, and then by target type.
# The 2a+b matched mutant parents each see N-k resident targets, except that
# a mixed pair deletes exactly one possible mutant-parent/resident-target event.
birth=(2*a+b)*(n-k)/(n-2)-b/(n-2)+c*(n-k)/(n-1)
death=k*(n-k)/(n-1) # every resident parent sees every mutant target
# Each of b mixed pairs becomes full by any mutant parent except its partner.
# A full pair ceases to be full if either of its 2a mutants is overwritten.
fullpair=b*((2*a+b-1)/(n-2)+c/(n-1))-2*a*(n-k)/(n-1)
alpha=n-2+(m-2)/n
q=S.factor(alpha*(birth-death)+fullpair+
           m*((n-2*k)*(birth-death)-(birth+death))/(n*(n-1)))
D=n*(n-2)*(n-1)**2
poly=S.Poly(S.expand(S.cancel(-D*q).subs(n,2*m+c+r).subs(m,a+b+d)),a,b,d,c,r)
fall=lambda x,j:S.prod(x-i for i in range(j))
a2,a3,b2,b3,b4,d2,d3=fall(a,2),fall(a,3),fall(b,2),fall(b,3),fall(b,4),fall(d,2),fall(d,3)
F={
 (2,1):2*(a+b+d), (1,2):2*(a+b+d),
 (2,0):4*a*d+b2+6*b*d+4*d**2,
 (0,2):4*a+b**2+2*b*d+b,
 (1,1):2*(4*a*a+7*a*b+8*a*d+2*b*b+2*b2+9*b*d+2*d*d+2*d2),
 (1,0):16*a*a*d+4*a*b2+36*a*b*d+24*a*d*d+4*b3+5*b2+24*b2*d+28*b*d2+38*b*d+8*d3+16*d2,
 (0,1):8*a3+28*a2+8*a+12*a2*b+8*a*b2+24*a*b+8*a2*d+20*a*d+12*a*b*d+4*b3+9*b2+2*b+12*b2*d+8*b*d2+18*b*d,
 (0,0):d*(16*a3+72*a2+40*a)+b*d*(32*a2+84*a+24)+d2*(16*a2+24*a)+b2*(2+48*d)+b3*(6+2*a+14*d)+2*b4+28*a*b2*d+b*d2*(36+24*a)+20*b2*d2+8*b*d3,
}
assert S.expand(poly.as_expr()-sum(v*c**i*r**j for (i,j),v in F.items()))==0
coefficient_checks={f'{i},{j}':S.expand(S.Poly(poly.as_expr(),c,r).coeff_monomial(c**i*r**j)-v)==0 for (i,j),v in F.items()}
assert all(coefficient_checks.values())
A,B,E=S.symbols('A B E')
regions=[]
for av in (0,1+A):
 for bv in (0,1,2+B):
  for dv in (0,1+E):
   p=S.Poly(S.expand(poly.as_expr().subs({a:av,b:bv,d:dv})),A,B,E,c,r)
   assert all(co>=0 for co in p.coeffs())
   regions.append({'a':str(av),'b':str(bv),'d':str(dv),'terms':len(p.terms()),'minimum_coefficient':int(min(p.coeffs()))})
W=alpha*k+a+m*k*(n-k)/(n*(n-1))
C=n*(n-2)+2*m-2
assert S.factor(W.subs({a:m,b:0,c:n-2*m})-C)==0
for sub in ({a:0,b:1,c:0},{a:0,b:0,c:1}):
 assert S.factor(W.subs(sub)-C/n)==0
 assert S.factor(q.subs(sub)+2*m/(n*(n-1)))==0
assert S.factor(q.subs({m:0,a:0,b:0}))==0
# Verify the source m=1 exact hitting function for every integer pair occupancy.
source_m1_Q=S.factor((n-2)*(birth-death)+fullpair)
for av,bv in ((0,0),(0,1),(1,0)):
 assert S.factor(source_m1_Q.subs({a:av,b:bv,n:2+c+r}))==0
stored=json.loads((RESEARCH/'evidence'/'uniform-certificate.json').read_text())
stored_poly=sum(t['coefficient']*S.prod(x**j for x,j in zip((a,b,d,c,r),t['powers'])) for t in stored['P_terms'])
assert S.expand(stored_poly-poly.as_expr())==0
p_eval=S.lambdify((a,b,d,c,r),poly.as_expr(),'math')

def scaled_w(mask,n,m):
    k=mask.bit_count()
    full=sum(((mask>>(2*j))&3)==3 for j in range(m))
    return ((n*(n-2)+m-2)*(n-1))*k+n*(n-1)*full+m*k*(n-k)

def raw_case(n,m):
    values=[scaled_w(x,n,m) for x in range(1<<n)]
    events=zero=0
    for mask,old in enumerate(values):
        z=0
        for parent in range(n):
            mut=bool(mask&(1<<parent))
            matched=mut and parent<2*m
            for target in range(n):
                if target==parent or (matched and target==(parent^1)): continue
                new=mask|(1<<target) if mut else mask&~(1<<target)
                z+=(values[new]-old)*((n-1) if matched else (n-2))
                events+=1
        aa=sum(((mask>>(2*j))&3)==3 for j in range(m))
        bb=sum(((mask>>(2*j))&3) in (1,2) for j in range(m))
        cc=(mask>>(2*m)).bit_count()
        assert -z==p_eval(aa,bb,m-aa-bb,cc,n-2*m-cc),(n,m,mask,z)
        assert z<=0
        if mask.bit_count()==1: assert z==-2*m*(n-2)*(n-1)
        zero+=(z==0)
    return {'N':n,'m':m,'states':1<<n,'parent_target_events':events,'zero_drift_states':zero,'all_singletons_checked':n}

cases=[raw_case(nn,mm) for nn,mm in [(3,1),(4,2),(12,5),(12,6),(13,6),(3,0),(7,0)]]
out={'status':'all assertions passed','method':'independent birth/death/full-pair event counts, no owner transition-table import',
     'symbolic_QW':str(q),'independent_polynomial_terms':len(poly.terms()),'eight_coefficient_identities':coefficient_checks,
     'integer_regions':regions,'calibration_and_two_singleton_drifts':True,'m0_control_drift_zero':True,
     'source_m1_hitting_function_verified':True,
     'raw_labelled_cases':cases,'total_raw_labelled_states':sum(v['states'] for v in cases),
     'total_parent_target_events':sum(v['parent_target_events'] for v in cases),
     'scope':'all N>=3,1<=m<=floor(N/2) proved by general event-count identities and integer positivity; finite labels are supplementary',
     'MATH_JOB_ID_present':bool(os.environ.get('MATH_JOB_ID'))}
(ROOT/'results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ('status','independent_polynomial_terms','total_raw_labelled_states','total_parent_target_events')}))
