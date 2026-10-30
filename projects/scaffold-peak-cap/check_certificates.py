"""Separate checker of serialized original-quadratic interval certificates.

Does not import the generator. No floating arithmetic or Python assert is used.
"""
from fractions import Fraction
from pathlib import Path
import json
import sys


def need(ok, message):
    if not ok:
        raise ValueError(message)


def rat(value):
    need(type(value) in (str,int), 'rationals must be strings or integers')
    return Fraction(value)


def poly(u,q,a):
    return u*u+(q-a-1)*u-q


def verify(case):
    p=case['p']
    need(type(p) is int and p>=2, 'p')
    need(case['m']==p+1 and case['occupied']==list(range(1,p+1))
         and case['vacant']==[p+1], 'occupancy')
    c=list(map(rat,case['centers']))
    need(c==[Fraction(1,10)]+[Fraction(2*j-3,2*j-1) for j in range(2,p+1)],'centers')
    points=[]
    for j in range(2,p+1):
        points.extend([Fraction(4*j-7,4*j-3),Fraction(4*j-5,4*j-1)])
    points.extend([Fraction(p,p+1),Fraction(3)])
    need(len(points)==2*p and all(a<b for a,b in zip(points,points[1:])), 'ordered witnesses')
    need(all(x!=1 and all(x!=t for t in c) for x in points), 'no limit singularities')
    E=[sum((x*t/(x-t)**2 for t in c),Fraction(0))+x*(1+x)/abs(x-1)**3 for x in points]
    emax=max(E)
    M=case['M']
    z=8*emax
    need(type(M) is int and M==1-((-z.numerator)//z.denominator), 'explicit M')
    need(rat(case['Emax'])==emax and M>8*emax,'M bound')
    K=list(map(rat,case['K'])); alpha=list(map(rat,case['alpha']))
    need(K==[M/t for t in c]+[Fraction(M)] and alpha==[Fraction(M-1)]*(p+1),'parameters')
    need(all(t>0 for t in K+alpha), 'positive finite parameters')
    witness=case['witnesses']
    need(len(witness)==2*p,'witness count')
    signs=[]
    for index,(x,e,w) in enumerate(zip(points,E,witness)):
        need(rat(w['X'])==x and rat(w['E'])==e and rat(w['error_bound'])==e/M,'witness arithmetic')
        count=sum(t<x for t in c)
        mu=x/(1-x) if x<1 else 1/(x-1)
        exact_limit=1-count+mu
        expected=Fraction(1,4) if index<2*p-2 and index%2==0 else (
            Fraction(-1,4) if index<2*p-2 else (Fraction(1) if index==2*p-2 else Fraction(3,2)-p))
        need(exact_limit==expected==rat(w['F_infty']), 'limit step counts')
        need(e/M<Fraction(1,8),'error margin')
        roots=w['U_brackets']
        need(len(roots)==p+1,'root count')
        tlo=[];thi=[];br=[]
        for k,a,pair in zip(K,alpha,roots):
            need(len(pair)==2,'interval shape')
            lo,hi=map(rat,pair);q=k*x
            # Product of roots is -q<0, so these signs isolate the unique positive root.
            need(1<lo<hi and poly(lo,q,a)<0<poly(hi,q,a),'U root isolation')
            br.append((lo,hi))
            tlo.append(q/(hi*hi+q));thi.append(q/(lo*lo+q))
        lo,hi=br[-1];q=K[-1]*x
        ml=(lo-1)*q/(hi*hi+q);mh=(hi-1)*q/(lo*lo+q)
        lower=1-sum(thi[:-1],Fraction(0))+ml
        upper=1-sum(tlo[:-1],Fraction(0))+mh
        need(lower<upper,'nonempty derivative interval')
        need(list(map(rat,w['F_bounds']))==[lower,upper],'serialized derivative bounds')
        sign=1 if index%2==0 else -1
        need(w['expected_sign']==sign,'alternating sign metadata')
        need((sign==1 and lower>0) or (sign==-1 and upper<0),'finite derivative sign')
        signs.append(sign)
    return {'p':p,'M':M,'exact_valid':True,'certified_signs':signs,
            'certified_disjoint_maximum_intervals':p}


def main():
    path=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('cases.json')
    data=json.loads(path.read_text())
    need(data['schema']=='one-vacancy-rational-radical-v1','schema')
    result=[verify(case) for case in data['cases']]
    need(len({r['p'] for r in result})==len(result),'distinct cases')
    print(json.dumps({'status':'exact finite certificates valid','cases':result},indent=2))


if __name__=='__main__':
    main()
