"""Exact replay for the frozen proportional-scaling correction candidate.

All decisive operations use QQ polynomial identities, rational matrix tests,
and rational interval endpoints. No floating-point root is used.
"""
from pathlib import Path
import json
import sympy as S


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


u,v,N = S.symbols('u v N', positive=True)
nu_birth=S.Rational(34,5)
mu_s=S.Integer(5)
mortality=S.Matrix([1,17])
P=S.Matrix([S.Rational(1,10),S.Rational(9,10)])
beta=S.Matrix([24,1])
M=S.diag(*mortality)
s=1-u-v
y=S.Matrix([u,v])
g=nu_birth-mu_s*s-(mortality.T*y)[0]
lam=(beta.T*y)[0]
require(all(x>0 for x in [nu_birth,mu_s,*mortality,*P,*beta]), 'positive parameters')
require(sum(P)==1 and M==S.diag(1,17),'probability/mortality conditions')
# Both P entries are positive, hence both postinfection nodes are reached from S.
# Independent construction in numbers, including the nonlinear 1/N incidence.
state=N*S.Matrix([s,u,v])
force=(beta.T*state[1:3,0])[0]/sum(state)
F_numbers=S.Matrix([nu_birth*sum(state)-mu_s*state[0]-force*state[0],
                     P[0]*force*state[0]-mortality[0]*state[1],
                     P[1]*force*state[0]-mortality[1]*state[2]])
Ndot=S.factor(sum(F_numbers))
require(S.factor(Ndot/N-g)==0,'total population derivative')
F_quotient=S.Matrix([S.factor((F_numbers[j+1]-y[j]*Ndot)/N) for j in range(2)])
F_compact=s*lam*P-M*y-g*y
F_printed=S.Matrix([(-4*u+v-64*u*u+95*u*v-v*v)/10,
                    (216*u-179*v-216*u*u-265*u*v+111*v*v)/10])
require(all(S.factor(a-b)==0 for a,b in zip(F_quotient,F_compact)),'quotient rule')
require(all(S.factor(a-b)==0 for a,b in zip(F_quotient,F_printed)),'printed polynomial reconstruction')

point={u:S.Rational(1,10),v:S.Rational(1,10)}
require(F_quotient.subs(point)==S.zeros(2,1),'exact positive equilibrium of proportions')
require(g.subs(point)==1,'balanced-growth rate')
require(F_numbers.subs(point)==state.subs(point),'numbers derivative equals state, not zero')
J=F_printed.jacobian([u,v])
J_saddle=J.subs(point)
require(J_saddle==S.Matrix([[-73,103],[1463,-1833]])/100,'saddle Jacobian')
require(J_saddle.det()==-S.Rational(211,125),'negative saddle determinant')
require(J_saddle.det()<0,'real eigenvalues of opposite sign')

J_dfe=J.subs({u:0,v:0})
require(S.trace(J_dfe)==-S.Rational(183,10) and J_dfe.det()==5,'DFE strict Hurwitz test')
# A second, attracting positive equilibrium: alpha is the unique root in (13/100,7/50).
q=41464*v*v-72863*v+9180
lo,hi=S.Rational(13,100),S.Rational(7,50)
require(q.subs(v,lo)>0 and q.subs(v,hi)<0,'algebraic-root IVT bracket')
require(S.diff(q,v).subs(v,hi)<0,'strictly decreasing q on the bracket')
u_alpha=S.Rational(125,754)-S.Rational(37,377)*v
for poly in F_printed:
    require(S.rem(S.together(poly.subs(u,u_alpha)),q,v)==0,'algebraic equilibrium exact remainder')
require(u_alpha.subs(v,hi)>0,'u positive throughout bracket')
require((1-u_alpha-v).subs(v,hi)>0,'s positive throughout bracket')
trace_alpha=S.rem(S.together(S.trace(J).subs(u,u_alpha)),q,v)
det_alpha=S.rem(S.together(J.det().subs(u,u_alpha)),q,v)
require(trace_alpha==S.Rational(13405,377)*v-S.Rational(187107,7540),'algebraic trace reduction')
require(det_alpha==S.Rational(66266,377)*v-S.Rational(16165,754),'algebraic determinant reduction')
require(trace_alpha.subs(v,hi)<0 and det_alpha.subs(v,lo)>0,'strict Hurwitz signs throughout root bracket')

out={
 'status':'exact-candidate-replay-passed; independent review required',
 'scope':'only unqualified proportional-recruitment/standard-incidence extension; base constant-recruitment theorem untouched',
 'parameters':{'nu':'34/5','mu_S':'5','mu':['1','17'],'P':['1/10','9/10'],'beta':['24','1'],'M':[['1','0'],['0','17']]},
 'normalized_rhs':[str(S.factor(p)) for p in F_printed],
 'positive_saddle':{'s':'4/5','u':'1/10','v':'1/10','population_growth_rate':'1','J':[[str(z) for z in row] for row in J_saddle.tolist()], 'trace':str(S.trace(J_saddle)),'det':str(J_saddle.det())},
 'disease_free_sink':{'s':'1','u':'0','v':'0','J':[[str(z) for z in row] for row in J_dfe.tolist()],'trace':str(S.trace(J_dfe)),'det':str(J_dfe.det())},
 'positive_sink':{'v_definition':'unique root alpha in (13/100,7/50) of 41464*alpha^2-72863*alpha+9180','u':'125/754-37*alpha/377','s':'1-u-alpha','trace_reduction':str(trace_alpha),'det_reduction':str(det_alpha),'trace_upper_bound':str(trace_alpha.subs(v,hi)),'det_lower_bound':str(det_alpha.subs(v,lo))},
 'limits':['mu_1=1 < mu_S=5 is within literal source assumptions but is not a calibrated disease-mortality mechanism','stable states refer to fractions, not equilibria in absolute population numbers','no full finite-parameter classification, periodic-orbit claim, or historical-first claim']}
print(json.dumps(out,indent=2,sort_keys=True))
