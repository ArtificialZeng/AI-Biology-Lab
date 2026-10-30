"""Exact witnesses for a source-specific response-data audit; no numerical fits."""
from pathlib import Path
import json,hashlib,itertools
import sympy as S
R=Path(__file__).resolve().parent
model=json.loads((R/'model.json').read_text())
rates=S.symbols('a b c d e f');edges=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
def generator(values):
 A=S.zeros(4)
 for (i,j),k in zip(edges,values):A[i,j]=A[j,i]=k
 for j in range(4):A[j,j]=-sum(A[i,j] for i in range(4) if i!=j)
 return A
def valid(A):return all(sum(A[i,j] for i in range(4))==0 for j in range(4)) and all(A[i,j]>0 for i in range(4) for j in range(4) if i!=j)
def require(x,label):
 if not x:raise RuntimeError(label)
A=generator(rates);e1=S.eye(4)[:,0];F=S.Matrix([(A**k)[j,0] for j in [1,2] for k in [1,2,3]])
DF=F.jacobian(rates);A0=generator([S.Integer(1)]*6);D0=DF.subs(dict.fromkeys(rates,1))
expected=S.zeros(6)
for col,start in [(0,0),(1,3)]:
 for i,vv in enumerate([1,-8,48]):expected[start+i,col]=vv
require(D0==expected and D0.rank()==2,'uniform symmetric derivative')
Ap=generator(list(map(S.Rational,model['positive_symmetric_rates'])))
Dp=DF.subs(dict(zip(rates,map(S.Rational,model['positive_symmetric_rates']))))
require(Dp.det()==800 and Dp.rank()==6,'positive local repair')
v=S.Matrix(list(map(S.Rational,model['krylov_v'])));b=S.Matrix(list(map(S.Rational,model['perturbation_b'])))
epsilon=S.Rational(model['epsilon']);B=b*v.T;At=Ap+epsilon*B
require(valid(Ap) and valid(At),'positive complete generators')
Krylov=S.Matrix.hstack(e1,Ap*e1,Ap**2*e1)
require(Krylov.rank()==3 and v.T*Krylov==S.zeros(1,3),'krylov annihilation')
require((v.T*Ap**3*e1)[0]==1112,'next Krylov moment')
same=[At**k*e1==Ap**k*e1 for k in range(4)];require(all(same),'first three jets')
fourth=(At**4-Ap**4)*e1;require(fourth==S.Matrix([0,S.Rational(139,125),-S.Rational(139,125),0]),'different fourth moment')
imbalance=S.factor(At[1,0]*At[2,1]*At[0,2]-At[0,1]*At[1,2]*At[2,0]);require(imbalance==S.Rational(1,100),'nonreversible cycle witness')
C=S.Matrix([[S.Rational(t) for t in row] for row in model['circulant']]);require(valid(C),'circulant contract')
s=S.symbols('s');den=s*(s+6)*(s*s+10*s+26)
require(S.expand(C.charpoly(s).as_expr()-den)==0,'characteristic polynomial')
response=(s*S.eye(4)-C).inv()*e1
numerators=[2*s*s+18*s+39,s*s+13*s+39];reduced=[]
for out,num in zip([1,2],numerators):
 require(S.cancel(response[out]-num/den)==0,'resolvent identity')
 require(S.gcd(num,den)==1,'no cancelled poles')
 reduced.append({'output_species':out+1,'numerator_coefficients':S.Poly(num,s).all_coeffs(),'denominator_coefficients':S.Poly(den,s).all_coeffs(),'gcd':str(S.gcd(num,den)),'quadratic_remainder':str(S.rem(num,s*s+10*s+26,s))})
require(S.discriminant(s*s+10*s+26,s)==-4,'nonreal pole pair')
# Check the derivative at A0 through independent matrix power directional sums.
directional=[]
for edge in edges:
 H=S.zeros(4);i,j=edge;H[i,j]=H[j,i]=1;H[i,i]=H[j,j]=-1
 require(S.ones(1,4)*H==S.zeros(1,4),'direction column conservation')
 require(A0*H==H*A0==-4*H,'uniform tangent relation')
 for k in [1,2,3]:
  dAk=sum((A0**l*H*A0**(k-1-l) for l in range(k)),S.zeros(4))
  require(dAk==k*(-4)**(k-1)*H,'directional derivative')
  directional.append(True)
mutations=[]
wrongH=S.zeros(4);wrongH[0,1]=wrongH[1,0]=1
require(S.ones(1,4)*wrongH!=S.zeros(1,4),'missing diagonal mutation detected');mutations.append('omit_Markov_diagonal_derivatives')
bad=C.copy();bad[0,1]=0;bad[1,1]=-(bad[0,1]+bad[2,1]+bad[3,1]);require(not valid(bad),'zero edge detected');mutations.append('zero_offdiagonal_breaks_complete_graph')
bad=Ap+epsilon*S.Matrix([0,1,0,0])*v.T;require(not valid(bad),'mass leakage detected');mutations.append('perturbation_column_sum_not_zero')
require(((v+e1).T*e1)[0]!=0,'wrong annihilator detected');mutations.append('wrong_Krylov_annihilator')
require(not valid(Ap+B),'large epsilon detected');mutations.append('large_perturbation_loses_positivity')
require(S.cancel(response[1]-(numerators[0]+1)/den)!=0,'wrong resolvent detected');mutations.append('changed_resolvent_numerator')
result={'scope':'Exact symbolic/fixed-witness algebra; no new general reversible spectral theorem or human/formal certification','sympy_version':S.__version__,'rate_order':['12','13','14','23','24','34'],'uniform_derivative':D0.tolist(),'uniform_rank':D0.rank(),'generic_derivative':Dp.tolist(),'generic_determinant':Dp.det(),'generic_rank':Dp.rank(),'positive_symmetric_base':Ap.tolist(),'perturbed_generator':At.tolist(),'epsilon':str(epsilon),'krylov_v':list(v),'perturbation_b':list(b),'all_full_state_moments_equal_orders0to3':same,'fourth_full_state_moment_difference':list(fourth),'cycle_imbalance':str(imbalance),'circulant':C.tolist(),'circulant_reduced_responses':reduced,'irreducible_quadratic_discriminant':-4,'directional_identity_cases':len(directional),'mutation_controls_rejected':mutations,'proof_sha256':hashlib.sha256((R/'proof.md').read_bytes()).hexdigest(),'model_sha256':hashlib.sha256((R/'model.json').read_bytes()).hexdigest()}
(R/'evidence').mkdir(exist_ok=True)
(R/'evidence/certificate.json').write_text(json.dumps(result,sort_keys=True,indent=2,default=str)+'\n')
print(json.dumps({'uniform_rank':2,'generic_determinant':800,'matching_derivative_orders':[0,1,2,3],'first_difference_order':4,'positive_generators':3,'resolvent_pairs_checked':2,'directional_cases':len(directional),'mutation_controls':len(mutations)}))
