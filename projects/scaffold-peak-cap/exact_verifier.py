"""Exact rational certificate checker; floating arithmetic is never decisive."""
from fractions import Fraction as Q
import json,sys
from math import isqrt

def rational(x):
    if isinstance(x,(float,bool)):
        raise TypeError('Use rational strings/integers, not floats/bools')
    return Q(x)

def poly(u,q,alpha): return u*u+(q-alpha-1)*u-q

def sqrt_brackets(z,bits=120):
    """Strict dyadic rational brackets for a positive rational square root."""
    assert isinstance(z,Q) and z>0
    scale=1<<bits
    floor=isqrt(z.numerator*z.denominator*scale*scale)//z.denominator
    lo=Q(floor,scale); hi=Q(floor+1,scale)
    if lo*lo==z: lo-=Q(1,scale)
    assert lo>=0 and lo*lo<z<hi*hi
    return lo,hi

def radical_bounds(K,alpha,x,bits=120):
    out=[]
    for k,a in zip(K,alpha):
        q=k*x; aa=a+1; rad=(q-aa)**2+4*q
        sl,sh=sqrt_brackets(rad,bits)
        lo=(aa-q+sl)/2; hi=(aa-q+sh)/2
        assert 1<=lo<hi and poly(lo,q,a)<0<poly(hi,q,a)
        out.append((lo,hi))
    return out

def bounds(K,alpha,x,bits=100):
    out=[]
    for k,a in zip(K,alpha):
        q=k*x; lo=Q(1); hi=1+a
        for _ in range(bits):
            mid=(lo+hi)/2
            z=poly(mid,q,a)
            if z<0: lo=mid
            elif z>0: hi=mid
            else:
                lo=(lo+mid)/2; hi=(mid+hi)/2
        assert 1<=lo<hi and poly(lo,q,a)<0<poly(hi,q,a)
        out.append((lo,hi))
    return out

def sign_interval(K,alpha,x,brackets):
    t=[]
    for k,a,(lo,hi) in zip(K,alpha,brackets):
        q=k*x
        assert 1<=lo<hi and poly(lo,q,a)<0<poly(hi,q,a)
        t.append((q/(hi*hi+q),q/(lo*lo+q)))
    lo,hi=brackets[2]; q=K[2]*x
    blo=(lo-1)*q/(hi*hi+q); bhi=(hi-1)*q/(lo*lo+q)
    return 1-t[0][1]-t[1][1]+blo, 1-t[0][0]-t[1][0]+bhi

def make(K,alpha,points,expected,bits=100):
    K=list(map(rational,K)); alpha=list(map(rational,alpha)); points=list(map(rational,points))
    assert len(K)==len(alpha)==3 and all(a>0 for a in K+alpha)
    assert all(x>0 for x in points) and all(a<b for a,b in zip(points,points[1:]))
    assert len(points)==len(expected)
    evidence=[]
    for x,s in zip(points,expected):
        br=bounds(K,alpha,x,bits)
        lo,hi=sign_interval(K,alpha,x,br)
        assert (s==1 and lo>0) or (s==-1 and hi<0), (str(x),str(lo),str(hi),s)
        evidence.append({'X':str(x),'brackets':[[str(a),str(b)] for a,b in br], 'F_bounds':[str(lo),str(hi)],'sign':s})
    return {'K':list(map(str,K)),'alpha':list(map(str,alpha)),'bits':bits,'points':evidence,'method':'Fraction-only quadratic sign isolation and interval evaluation','certified':True}

def verify(doc):
    K=list(map(rational,doc['K']));alpha=list(map(rational,doc['alpha']))
    assert len(K)==len(alpha)==3 and all(a>0 for a in K+alpha)
    xs=[];sg=[]
    for z in doc['points']:
        x=rational(z['X']); xs.append(x); assert x>0
        br=[tuple(map(rational,b)) for b in z['brackets']]
        assert len(br)==3
        lo,hi=sign_interval(K,alpha,x,br)
        s=z['sign']; sg.append(s)
        assert (s==1 and lo>0) or (s==-1 and hi<0)
        assert [str(lo),str(hi)]==z['F_bounds']
    assert all(a<b for a,b in zip(xs,xs[1:]))
    return {'exact_valid':True,'signs':sg,'at_least_strict_maxima':sum(a==1 and b==-1 for a,b in zip(sg,sg[1:]))}

if __name__=='__main__':
    if len(sys.argv)>1:
        print(json.dumps(verify(json.load(open(sys.argv[1])))))
    else:
        doc=make(['1']*3,['1/10','1/10','100'],['1/10','10','100','1000'],[1,-1,1,-1])
        print(json.dumps(doc,indent=2))
