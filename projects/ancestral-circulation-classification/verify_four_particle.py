"""Read-only exact polynomial certificate for the four-particle path."""
from pathlib import Path
import json,sympy as s
from four_particle_model import systems,require
HERE=Path(__file__).resolve().parent
t=s.symbols('t');u=s.symbols('u')
def load(label):
    j=json.loads((HERE/'evidence'/'four-particle'/(label+'-solution.json')).read_text())
    return s.Matrix([s.sympify(v) for v in j['numerators']]),s.sympify(j['denominator'])
def zero(v):return s.Poly(s.expand(v),t).is_zero
def build_certificate():
    require(s.__version__ == '1.13.3', 'SymPy 1.13.3 required')
    pop,states,ci,si,L,T,P,pr,H,hr=systems(4,(s.Integer(1),t,t+1))
    pn,pd=load('population');hn,hd=load('founder')
    require(all(zero(v) for v in L.T*pn),'Stationary polynomial residuals')
    require(zero(sum(pn)-pd),'Stationary normalization')
    require(all(zero(v) for v in T*hn),'Founder polynomial residuals')
    require(zero(4*hn[si[0,(3,0,0)]]-hd),'Founder anchor')
    for c in pop:
        hsum=sum(c[i]*hn[si[i,tuple(c[j]-int(i==j) for j in range(3))]] for i in range(3) if c[i])
        require(zero(hsum-hd),'All founder normalizations')
    pshift=s.Poly(s.expand(pd.subs(t,u+1)),u);hshift=s.Poly(s.expand(hd.subs(t,u+1)),u)
    require(all(v>0 for v in pshift.coeffs()),'Population denominator positive on t>1')
    require(all(v<0 for v in hshift.coeffs()),'Founder denominator negative on t>1')
    j=json.loads((HERE/'evidence/four-particle/symbolic-current.json').read_text());num=s.sympify(j['numerator']);den=s.sympify(j['denominator'])
    raw=s.Integer(0)
    for c in pop:
        if c[0]:
            d=list(c);d[0]-=1;raw+=pn[ci[c]]*c[0]*hn[si[2,tuple(d)]]
        if c[2]:
            d=list(c);d[2]-=1;raw-=pn[ci[c]]*c[2]*hn[si[0,tuple(d)]]
    require(zero(raw*den-num*pd*hd),'Current crossproduct')
    phase,remainder=s.div(num,2*t*(t-1)*(3*t*t+13*t+8),t)
    require(remainder==0,'Exact numerator factorization')
    phase=s.Poly(phase,t);coeffs=phase.all_coeffs();signs=[int(s.sign(v)) for v in coeffs if v]
    variations=sum(x!=y for x,y in zip(signs,signs[1:]))
    require(phase.degree()==19 and variations==1,'Descartes one-positive-root certificate')
    require(sum(1 for v in coeffs if v>0)==7 and sum(1 for v in coeffs if v<0)==13,'Complete coefficient pattern')
    require(phase.eval(1)<0,'Root is above1')
    require(all(v>0 for v in s.Poly(den,t).coeffs()),'Current denominator coefficient positivity')
    left=s.Rational(4769,446);right=s.Rational(3133,293)
    require(1<left<right and phase.eval(left)<0<phase.eval(right),'Exact root isolator')
    expected={2:s.Rational(-7553807516755111581622,1096028058625091839860416349),100:s.Rational(458226106942025952304289497342294745709790415205,4105486436463376678337301093835196669201815807466961979)}
    for x,v in expected.items():require(num.subs(t,x)/den.subs(t,x)==v,'Exact current witness')
    result={'ok':True,'N':4,'domain':'t>1, k=(1,t,t+1), unit triangle motion','population_states':15,'founder_unknowns':30,'all_population_and_founder_polynomial_residuals_zero':True,'founder_normalizations':15,'phase_polynomial_degree':19,'positive_coefficients':7,'negative_coefficients':13,'Descartes_sign_changes':variations,'unique_simple_root_above1':True,'root_interval':[str(left),str(right)],'current_sign_below_root':-1,'current_sign_above_root':1,'witnesses':{str(k):str(v) for k,v in expected.items()},'archive_modified':False,'minimum_population_for_within_order_reversal':4,'minimality_dependency':'sourceN2balance plus companionN3fullorderingtheorem','all_positive_N4_classification_claimed':False}
    return result
from archive_tools import run
if __name__ == '__main__':
    raise SystemExit(run(build_certificate, 'evidence/four-particle/exact-certificate.json', 'N4 fifteen-state and thirty-founder polynomial residuals; exact path threshold'))
