"""Exact algebra and finite reaction checks; not a proof of ODE entrainment."""
from pathlib import Path
from fractions import Fraction as F
import hashlib, json, itertools
import sympy as s
R=Path(__file__).resolve().parent
a,b,c,d,u,x,y=s.symbols('a b c d u x y',positive=True)
h=(a+c)/c; q=a*d/c; r=a*d/(a+c); y0=b*c/(a*d)
dx=(a+c)*u-b*x-d*x*y; dy=c*u-d*x*y
identities={}
def zero(name,expression):
 value=s.factor(expression)
 if value!=0:raise RuntimeError(name+': '+str(value))
 identities[name]=str(value)
zero('input_cancels_in_w',dy-dx/h-r*x*(y0-y))
eta=y-y0; z=x-h*eta
zero('z_derivative',dx-h*dy-q*x*eta)
zero('log_identity_numerator',q*x*eta-q*z*eta-q*h*eta**2)
zero('constant_equilibrium_x',dx.subs({x:a*u/b,y:y0}))
zero('constant_equilibrium_y',dy.subs({x:a*u/b,y:y0}))
u0,omega=s.symbols('u0 omega',positive=True);x0=a*u0/b
J=s.Matrix([dx,dy]).jacobian([x,y]).subs({x:x0,y:y0})
B=b+d*(y0+x0); K=b*d*x0
zero('jacobian_trace',s.trace(J)+B)
zero('jacobian_determinant',J.det()-K)
lam=s.symbols('s')
H=(lam*s.eye(2)-J).adjugate()*s.Matrix([a+c,c])
zero('y_transfer_numerator',H[1]-c*lam)
zero('y_transfer_denominator',(lam*s.eye(2)-J).det()-(lam**2+B*lam+K))
v=s.symbols('v');den=(K-v)**2+B**2*v
zero('leading_deficit_peak_derivative',den-v*s.diff(den,v)-(K**2-v**2))
M,L,Umax,slack,ell=s.symbols('M L Umax slack ell',nonnegative=True)
zero('upper_x_face_decomposition',dx.subs(x,M)-((a+c)*(u-Umax)-b*slack-d*M*y)-(a+c)*Umax+b*(M-slack))
zero('slanted_face',r*x*(y0-y).subs(y,y0+ell+x/h)+r*x*(ell+x/h))
p,vv,x2,y1=s.symbols('p vv x2 y1')
zero('incremental_product',(x2+p)*(y1)-(x2)*(y1+vv)-(y1*p-x2*vv))
epsilon,Y,m=s.symbols('epsilon Y m',positive=True)
zero('weighted_first_column',(-b-d*y+(1+epsilon)*d*y)-(-b+epsilon*d*y))
zero('weighted_second_column',(-d*x+d*x/(1+epsilon))+epsilon*d*x/(1+epsilon))
zero('first_column_at_Y',(-b+epsilon*d*y).subs({epsilon:b/(2*d*Y),y:Y})+b/2)

rates=list(itertools.product([F(1,2),F(1),F(2)],repeat=4))
states=list(itertools.product([F(0),F(1,3),F(1),F(2),F(5)],repeat=2))
arrows=[({'S':1},{'X':1},'a'),({'X':1},{'S':1},'b'),({'S':1},{'X':1,'Y':1},'c'),({'X':1,'Y':1},{'S':1},'d')]
def reaction_rhs(av,bv,cv,dv,uv,xv,yv,network=arrows):
 rate=dict(zip('abcd',(av,bv,cv,dv)));conc={'S':uv,'X':xv,'Y':yv};out={'X':F(0),'Y':F(0)}
 for inputs,outputs,k in network:
  flux=rate[k]
  for species,power in inputs.items():flux*=conc[species]**power
  for species in out:out[species]+=flux*(outputs.get(species,0)-inputs.get(species,0))
 return out['X'],out['Y']
cases=0;face_cases=0
for av,bv,cv,dv in rates:
 hv=(av+cv)/cv;rv=av*dv/(av+cv);qv=av*dv/cv;y0v=bv*cv/(av*dv)
 for xv,yv in states:
  uv=F(7,3)
  actual=reaction_rhs(av,bv,cv,dv,uv,xv,yv)
  expected=((av+cv)*uv-bv*xv-dv*xv*yv,cv*uv-dv*xv*yv)
  if actual!=expected:raise RuntimeError('reaction mismatch')
  if actual[1]-actual[0]/hv!=rv*xv*(y0v-yv):raise RuntimeError('w mismatch')
  if actual[0]-hv*actual[1]!=qv*xv*(yv-y0v):raise RuntimeError('z mismatch')
  cases+=1
 for Uv,slackv,ellv in itertools.product([F(1),F(3)],[F(0),F(2)],[F(0),F(1)]):
  Mv=(av+cv)*Uv/bv+slackv;Lv=y0v+ellv
  for uv,yv in itertools.product([Uv/2,Uv],[F(0),Lv+Mv/hv]):
   if reaction_rhs(av,bv,cv,dv,uv,Mv,yv)[0]>0:raise RuntimeError('x face outward')
  xv=Mv/2;yv=Lv+xv/hv;rx,ry=reaction_rhs(av,bv,cv,dv,Uv,xv,yv)
  if ry-rx/hv>0:raise RuntimeError('slant outward')
  face_cases+=1

mutations=[]
for name,changed in [
 ('lose_X_in_split',[arrows[0],arrows[1],({'S':1},{'Y':1},'c'),arrows[3]]),
 ('lose_Y_consumption',[*arrows[:3],({'X':1},{'S':1},'d')]),
 ('catalytic_reverse_X',[*arrows[:3],({'X':1,'Y':1},{'S':1,'X':1},'d')]),
]:
 if reaction_rhs(F(2),F(3),F(5),F(7),F(11),F(13),F(17),changed)==reaction_rhs(F(2),F(3),F(5),F(7),F(11),F(13),F(17)):
  raise RuntimeError('mutation not detected')
 mutations.append(name)
if s.factor(dy.subs({x:a*u/b,y:y0+1}))==0:raise RuntimeError('wrong baseline accepted')
mutations.append('wrong_baseline')
if s.diff(a*u/b,u)==0:raise RuntimeError('source_output_X_index_accepted')
mutations.append('wrong_adapted_output_X')
if F(1)/(1-F(1,4))-1<=0:raise RuntimeError('bad contraction weight not rejected')
mutations.append('decrease_second_weight_breaks_contraction')
result={'scope':'Exact algebra and complete declared finite grids; no probability/ODE/formal proof certification','sympy_version':s.__version__,'symbolic_zero_identities':identities,'reaction_grid':{'rate_values':['1/2','1','2'],'parameter_tuples':len(rates),'xy_values':['0','1/3','1','2','5'],'states_per_parameter':len(states),'cases':cases},'invariant_polygon_face_cases':face_cases,'mutation_controls_rejected':mutations,'source_output_mapping':{'S':1,'X':3,'Y':2},'proof_sha256':hashlib.sha256((R/'proof.md').read_bytes()).hexdigest(),'model_sha256':hashlib.sha256((R/'model.json').read_bytes()).hexdigest()}
(R/'evidence').mkdir(exist_ok=True)
(R/'evidence/certificate.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps({'symbolic_identities':len(identities),'exact_reaction_cases':cases,'face_cases':face_cases,'mutations':len(mutations)}))
