from pathlib import Path
import json
import sympy as s

def need(condition, message):
    if not condition:
        raise ValueError(message)

def main():
    p=Path(__file__).resolve().parent
    w=json.loads((p/'witness.json').read_text())
    x,y,t=s.symbols('x y t')
    expected={(i,n-i) for n in range(4) for i in range(n+1)}
    vertices=[tuple(v) for v in w['complexes']]
    need(len(vertices)==10 and set(vertices)==expected,'complex domain')
    exceptions={}
    for e in w['exception_rates']:
        v,z=tuple(e['source']),tuple(e['target'])
        need(v in expected and z in expected and v!=z,'invalid reaction')
        need((v,z) not in exceptions,'duplicate replacement')
        exceptions[v,z]=s.Rational(e['rate'])
    base=s.Rational(w['base_rate'])
    need(base>0 and all(k>0 for k in exceptions.values()),'nonpositive rate')
    F=[s.Integer(0),s.Integer(0)];edges=0
    for v in vertices:
        for z in vertices:
            if v==z:
                continue
            k=exceptions.get((v,z),base)
            need(k>0,'rate')
            for j in range(2):
                F[j]+=k*x**v[0]*y**v[1]*(z[j]-v[j])
            edges+=1
    need(edges==90,'edge count')
    F=[s.expand(f) for f in F]
    need(isinstance(w['F'],list) and len(w['F'])==2,'stored field must have two components')
    need(all(s.expand(f-s.sympify(stored))==0 for f,stored in zip(F,w['F'])),'rate-to-field identity')
    need(all(s.Poly(f,x,y).total_degree()<=3 for f in F),'field degree')
    h=(x-2)**2+4*(y-5)**2-1
    need(s.expand(h-s.sympify(w['h']))==0,'ellipse')
    L=[sum(s.Rational(w['radial_coefficients'][3*j+i])*m for i,m in enumerate([1,x,y])) for j in range(2)]
    q=s.sympify(w['q']);K=s.sympify(w['cofactor'])
    need(s.Poly(q,x,y).total_degree()<=2,'q degree')
    need(s.expand(F[0]-h*L[0]+q*s.diff(h,y))==0,'first Darboux component')
    need(s.expand(F[1]-h*L[1]-q*s.diff(h,x))==0,'second Darboux component')
    need(s.expand(K-s.diff(h,x)*L[0]-s.diff(h,y)*L[1])==0,'cofactor')
    need(s.expand(s.diff(h,x)*F[0]+s.diff(h,y)*F[1]-K*h)==0,'invariance')
    C=(1-t*t)/(1+t*t);S=2*t/(1+t*t)
    qt=s.cancel(q.subs({x:2+C,y:5+S/2}))
    need(s.cancel(qt-s.sympify(w['q_rational_numerator'])/s.sympify(w['q_rational_denominator']))==0,'tangent rational identity')
    a,b=s.fraction(qt);pa,pb=s.Poly(a,t),s.Poly(b,t)
    need(pa.count_roots(-s.oo,s.oo)==0,'tangent zero on oval')
    need(pb.count_roots(-s.oo,s.oo)==0,'parameter denominator zero')
    need(qt.subs(t,0)>0 and q.subs({x:1,y:5})>0,'tangent sign or excluded point')
    rotations=[(C,S),(-S,C),(-C,-S),(S,-C)]
    parts=[s.cancel(K.subs({x:2+u,y:5+v/2})/(2*q.subs({x:2+u,y:5+v/2})*(1+t*t))) for u,v in rotations]
    G=s.cancel(sum(parts))
    need(s.cancel(G-s.sympify(w['fourquarter_integrand_numerator'])/s.sympify(w['fourquarter_integrand_denominator']))==0,'Floquet integrand identity')
    a,b=s.fraction(G)
    need(s.Poly(a,t).count_roots(0,1)==0,'Floquet numerator sign changes')
    need(s.Poly(b,t).count_roots(0,1)==0,'Floquet denominator pole')
    need(G.subs(t,0)<0,'Floquet integrand sign')
    result={'status':'exact-witness-checks-passed','directed_edges':edges,'source_and_target_degree_max':3,'all_rates_positive':True,'complete_reversible_graph':True,'Darboux_identities_exact':True,'ellipse_positive_and_regular':'analytic in proof.md','tangent_no_real_parameter_roots':True,'excluded_point_q':str(q.subs({x:1,y:5})),'fourquarter_Floquet_integrand_strictly_negative':True,'conclusion_with_analytic_argument':'one hyperbolic attracting algebraic limit cycle','finite_witness_not_global_classification':True}
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    main()
