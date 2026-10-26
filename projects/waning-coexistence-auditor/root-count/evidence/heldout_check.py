#!/usr/bin/env python3
"""Independent held-out rational cases; stdout only; no proposer checker imported."""
import json
import sympy as S

s=S.symbols('s')
def require(ok, message):
    if not ok: raise RuntimeError(message)
def sign(x):
    z=S.sign(S.simplify(x));require(z in (-1,0,1),'Sign not exactly resolved: '+str(x));return int(z)
def rat(x): return S.Rational(x)
def parameters(mu,g1,g2,p,q,d,eta1,sig1,sig2,theta1,theta2,theta12):
    mu,g1,g2,p,q,d,eta1,sig1,sig2,theta1,theta2,theta12=map(rat,[mu,g1,g2,p,q,d,eta1,sig1,sig2,theta1,theta2,theta12])
    eta2=(g1+mu)*(g2+mu)*(p-d)*(q-d)/(g1*g2*d*d*eta1)
    return dict(mu=mu,g1=g1,g2=g2,p=p,q=q,eta1=eta1,eta2=eta2,sig1=sig1,sig2=sig2,theta1=theta1,theta2=theta2,theta12=theta12)
def formulas(z):
    mu,g1,g2,p,q,e1,e2,h1,h2,t1,t2,t12=[z[k] for k in ['mu','g1','g2','p','q','eta1','eta2','sig1','sig2','theta1','theta2','theta12']]
    r1=(q-s)/(e2*h2);r2=(p-s)/(e1*h1);v1=g1+mu;v2=g2+mu
    D=S.Poly(g1*g2*s*s-h1*h2*v1*v2*r1*r2,s)
    A1=r1*(g2*s*(mu+t1)+h2*v2*(mu+t2)*r2)
    A2=r2*(g1*s*(mu+t2)+h1*v1*(mu+t1)*r1)
    W1=s+h1*r2*(1+g1/(mu+t12));W2=s+h2*r1*(1+g2/(mu+t12))
    numerator=S.expand((s+r1+r2)*D.as_expr()+A1*W1+A2*W2)
    N=S.cancel(numerator/D.as_expr())
    return D,numerator,N

def isolate_feasible(P,D,p):
    # Rational isolating intervals; Sturm checks exclude a D-pole from each interval.
    for exponent in range(3,18):
        roots=P.intervals(eps=S.Rational(1,10**exponent));accepted=[];pending=False
        for (lo,hi),mult in roots:
            if hi<=0 or lo>=p: continue
            if lo<=0 or hi>=p:
                pending=True;break
            if D.count_roots(lo,hi)!=0 or D.eval(lo)==0 or D.eval(hi)==0:
                pending=True;break
            state=sign(D.eval((lo+hi)/2))
            if state>0: accepted.append({'interval':[str(lo),str(hi)],'multiplicity':int(mult)})
        if not pending:return accepted
    raise RuntimeError('Could not separate polynomial roots from feasible endpoints')

def main():
    require(S.__version__=='1.13.3','Pinned SymPy 1.13.3 required')
    bases=[
        parameters(1,1,1,1,2,'1/10',30,'1/2','1/2',0,0,0),
        parameters(2,3,5,2,3,'6/5',7,3,2,0,4,0),
        parameters(1,2,4,'3/2','5/2','6/5',2,'1/3',5,3,0,2),
        parameters(1,1,2,1,1,'1/4',4,2,'1/3',1,1,1)]
    cases=[]
    for j,z in enumerate(bases):
        _,_,N=formulas(z);Np=S.cancel(N.subs(s,z['p']))
        for label,level in [('half',Np/2),('endpoint',Np),('above',Np+1),('three_quarters',3*Np/4)]:
            cases.append((f'rational-pole-{j}-{label}',z,level))
    standard=dict(mu=rat(1),g1=rat(1),g2=rat(1),p=rat(1),q=rat(2),eta1=rat(20),eta2=rat(10),sig1=rat(1),sig2=rat(1),theta1=rat('535/262'),theta2=rat(1),theta12=rat(1))
    for label,offset in [('fold',0),('fold-below',rat('-1/100')),('fold-above',rat('1/100'))]:
        cases.append((label,standard,rat('872/655')+offset))
    pair={**standard,'theta1':rat('535/216')};cases.append(('two-rational-roots',pair,rat('2011/1440')))
    endpoint_tangent={**standard,'theta1':rat('54/11')};cases.append(('endpoint-a-b-zero',endpoint_tangent,rat('93/55')))
    tie={**standard,'q':rat(1),'eta1':rat(10),'theta1':rat(1)};cases.append(('equal-critical-thresholds',tie,rat(1)))
    collapse={**standard,'eta1':rat(1),'eta2':rat(2),'theta1':rat(1)}
    cases.append(('P-degree-two',collapse,rat(3)))
    cases.append(('P-degree-one-endpoint',collapse,rat('5/2')))
    records=[]
    for name,z,s0 in cases:
        require(all(v>0 for k,v in z.items() if not k.startswith('theta')),'Positive parameter failure')
        require(all(z[k]>=0 for k in ['theta1','theta2','theta12']) and s0>0,'Domain failure')
        D,numerator,N=formulas(z);p=z['p'];mu=z['mu'];P=S.Poly(S.expand(mu*(numerator-s0*D.as_expr())),s)
        require(P.degree()<=3,'Population polynomial exceeds cubic')
        candidates=[root for root in S.solve(D.as_expr(),s) if sign(root)>0 and sign(p-root)>0]
        require(len(candidates)==1,'Feasible pole not unique');d=candidates[0];ell=p-d
        a=S.simplify(P.eval(p));b=S.simplify(3*a-ell*P.diff().eval(p));c=S.simplify(3*P.eval(d)+ell*P.diff().eval(d));D0=S.simplify(P.eval(d))
        require(sign(D0)>0,'Left feasible numerator is not positive')
        A,B,C=map(sign,[a,b,c])
        u3,u2,u1,u0=[P.nth(k) for k in [3,2,1,0]]
        binary_disc=18*u3*u2*u1*u0-4*u2**3*u0+u2*u2*u1*u1-4*u3*u1**3-27*u3*u3*u0*u0
        disc=sign(binary_disc)
        t=S.symbols('t');Q=a*t**3+b*t**2+c*t+D0
        require(S.simplify(S.cancel((1+t)**3*P.as_expr().subs(s,(d+p*t)/(1+t))-Q))==0,'Mobius coefficient identity failed')
        if A<0 or (A==0 and B<0):expected=[1]
        elif A==0:expected=[]
        elif (B<0 or C<0) and disc==1:expected=[1,1]
        elif (B<0 or C<0) and disc==0:expected=[2]
        else:expected=[]
        roots=isolate_feasible(P,D,p);actual=sorted(r['multiplicity'] for r in roots)
        require(actual==expected,'Classification mismatch: '+name+' '+str(actual)+' versus '+str(expected))
        if name=='fold': require(S.factor(P.as_expr())==S.factor((3*s-2)**2*(595*s+169)/10480),'Fold factor mismatch')
        if name=='two-rational-roots':require(P.eval(rat('2/3'))==0 and P.eval(rat('3/4'))==0,'Pair roots not exact')
        if name=='equal-critical-thresholds':require(P.eval(rat('11/31'))==0,'Tie root not exact')
        if name=='P-degree-two':require(P.degree()==2 and S.expand(P.as_expr()-(s*s-34*s+32)/2)==0,'Quadratic collapse identity failed')
        if name=='P-degree-one-endpoint':require(P.degree()==1 and S.expand(P.as_expr()-14*(1-s))==0,'Linear collapse identity failed')
        records.append({'name':name,'parameters':{k:str(v) for k,v in z.items()},'s0':str(s0),'D':str(D.as_expr()),'D_degree':D.degree(),'D_leading_sign':sign(D.LC()),'P':str(P.as_expr()),'P_degree':P.degree(),'feasible_pole':str(d),'endpoint_a_sign':A,'b_sign':B,'c_sign':C,'cubic_discriminant_sign':disc,'predicted_multiplicities':expected,'exact_isolated_feasible_roots':roots})
    print(json.dumps({'ok':True,'scope':'24 exact rational held-out parameter cases; no infinite theorem certification','sympy':S.__version__,'method':'exact rational polynomial root isolation and Sturm separation from D-poles; exact algebraic coefficient signs','cases':records},indent=2,sort_keys=True))
if __name__=='__main__':
    try:main()
    except Exception as exc:
        print(json.dumps({'ok':False,'error':str(exc)},sort_keys=True));raise SystemExit(1)
