"""Exact symbolic check of the first slow-manifold jet and implemented F1."""
from pathlib import Path
import re,json,hashlib,argparse
import sympy as sp
u,v,w,a,b,y,e=sp.symbols('u v w a b y e')
p=8-a-b;s=40-u-v-w-a-b;A=u+v
h0=u+100*b
full=sp.Matrix([y*(1-e*A)-2*u,u-v,v-100*w*p+a,100*w*p-sp.Rational(19,10)*a,3*p*(s-e*y)-101*b])
g=-y+e*A*y+u+100*b
F0=full.subs({y:h0,e:0})
lie=sp.expand((sp.Matrix([h0]).jacobian([u,v,w,a,b])*F0)[0])
h1=sp.expand(A*h0-lie)
h=h0+e*h1
invariance=e*(sp.Matrix([h]).jacobian([u,v,w,a,b])*full.subs(y,h))[0]-g.subs(y,h)
checks_invariance=[sp.expand(invariance).coeff(e,k)==0 for k in [0,1]]
F1=full.subs(y,h).applyfunc(lambda z:sp.expand(z).coeff(e,1))
cpp=(Path(__file__).resolve().parent/'verify_period_sensitivity.cpp').read_text()
expr=re.search(r'IMap vf\("par:r;var:u,v,w,a,b,e;fun:([^;]+);"\)',cpp).group(1)
locals_={str(q):q for q in [u,v,w,a,b,e]};locals_['r']=sp.Rational(19,10)
parsed=sp.Matrix([sp.sympify(z,locals=locals_) for z in expr.split(',')])
checks_cpp=[sp.expand(parsed[i]-(F0[i]+e*F1[i]))==0 for i in range(5)]+[parsed[5]==0]
if not all(checks_invariance) or not all(checks_cpp):
    raise RuntimeError('First-jet identity failed')
result={'status':'passed','arithmetic':'exact rational polynomial identities','invariance_coefficients_orders_0_1_zero':checks_invariance,'CPP_equal_F0_plus_eF1_and_eprime_zero':checks_cpp,'h0':str(h0),'h1':str(h1),'F1':[str(z) for z in F1],'CPP_sha256':hashlib.sha256(cpp.encode()).hexdigest(),'scope':'First-order jet identities only; smooth slow-manifold existence and periodic-branch linkage require analytic audit'}
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[4]/'artifacts/validated-cycle/first-jet-audit.json')
args=parser.parse_args()
args.output.parent.mkdir(parents=True,exist_ok=True)
args.output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
