"""Reconstruct the exact third-derivative certificate. No floating arithmetic."""
from pathlib import Path
import sympy as S
import json,hashlib
ROOT=Path(__file__).resolve().parent

def derive():
    t,u,m,n=S.symbols('t u m n')
    x,y=S.symbols('x y')
    p=t*(1-t)*(1+2*m); q=u*(1-u)*(1+2*n)
    def D(z):
        return S.diff(z,t)*p+S.diff(z,u)*q+S.diff(z,m)*m*(1+m)*(1-2*t)+S.diff(z,n)*n*(1+n)*(1-2*u)
    g=t+u-1;h=g*(g+1)
    h1=S.expand(D(h));h2=S.expand(D(h1));h3=S.expand(D(h2))
    N=S.Poly(S.expand(h*h*(h1-h3)+6*h*h1*h2-6*h1**3),t,u,m,n)
    den_power=max(a+b for (a,b,c,e),co in N.terms())
    assert den_power==9
    d=1+x+y
    num=S.Poly(S.expand(sum(co*(1+y)**a*(1+x)**b*d**(den_power-a-b)*m**c*n**e for (a,b,c,e),co in N.terms())),x,y,m,n)
    terms=[{'powers':list(powers),'coefficient':int(co)} for powers,co in num.terms()]
    assert len(terms)==322
    assert all(z['coefficient']<0 for z in terms)
    assert all(z['coefficient']==int(z['coefficient']) for z in terms)
    return terms

if __name__=='__main__':
    terms=derive()
    expected=json.load(open(ROOT/'coefficient-certificate.json'))
    assert expected['variables']==['x','y','m','n']
    assert expected['denominator']=='(1+x+y)^9'
    assert expected['terms']==terms
    raw=(ROOT/'coefficient-certificate.json').read_bytes()
    print(json.dumps({'exact_identity_reconstructed':True,'integer_coefficients':322,'all_strictly_negative':True,'certificate_sha256':hashlib.sha256(raw).hexdigest(),'meaning':'Proves L third derivative <0 on G>0; see proof.md for parameterization and Rolle argument'}))
