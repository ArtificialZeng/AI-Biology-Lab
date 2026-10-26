"""Literal monomials and exact-rational rank diagnostics, not an infinite proof."""
import flint,json
from flint import fmpq_mat

def require(condition,reason):
    if not condition:raise RuntimeError(reason)

def build_case(r):
    n=1<<r
    nu=[4*(n*x.bit_count()+x) for x in range(n)]
    mu=[-1]+[1]*(n-1)
    weights=[[1]*n for _ in range(3)]+[[v+1 for v in nu],[2]*n,[2+v for v in mu]]
    require(all(v>0 for block in weights for v in block),'strictly positive actual edge weights')
    rows=[[0]*(n*n) for _ in range(6*n)]
    margins=[];choices=[];deltas=[]
    for g in range(n):
        for h in range(n):
            k=g^h
            monomials=[[(0,g),(1,h),(2,k),(3,k),(4,g)],[(0,g),(1,h),(2,k),(3,h),(5,g)]]
            actual=[sum(weights[edge][label] for edge,label in m) for m in monomials]
            delta=actual[0]-actual[1]
            require(delta==nu[k]-nu[h]-mu[g],'actual full monomial comparison')
            require(delta!=0,'strict tropical chamber')
            first=delta<0
            expected=g!=0 and (k.bit_count(),k)<(h.bit_count(),h)
            require(first==expected,'Hamming then integer branch identity')
            require(g==0 or (delta<0)!=((nu[h]-nu[k]-mu[g])<0),'exactly one selected first in each nonzero XOR pair')
            for edge,label in monomials[0 if first else 1]:rows[edge*n+label][g*n+h]+=1
            margins.append(abs(delta));choices.append(1 if first else 2);deltas.append(delta)
    A=fmpq_mat(rows)
    basis=[]
    for i in range(n):
        v=[0]*(6*n);v[i]=1;v[5*n+i]=-1
        if i:v[4*n+i]=-1
        basis.append(v)
    for edge in [1,2,3]:
        v=[0]*(6*n)
        for i in range(n):
            v[edge*n+i]=1;v[5*n+i]=-1
            if i:v[4*n+i]=-1
        basis.append(v)
    v=[0]*(6*n);v[4*n]=1;basis.append(v)
    B=fmpq_mat(basis)
    require(B*A==fmpq_mat(n+4,n*n),'all proposed kernel vectors annihilate')
    require(B.rank()==n+4,'n+4 proposed vectors independent')
    rank=A.rank();expected_rank=5*n-4 if r>=3 else (4 if r==1 else 14)
    require(rank==expected_rank,'exact held-out or boundary rank')
    return {'r':r,'order':n,'shape':[6*n,n*n],'rank_Q':rank,'left_nullity':6*n-rank,'target_5n_minus_4':5*n-4,'inside_theorem':r>=3,'held_out_from_triage':r in [5,6],'minimum_strict_margin':min(margins),'first_count':choices.count(1),'second_count':choices.count(2),'kernel_basis_rank':B.rank(),'actual_positive_weights':True,'literal_comparisons_checked':n*n,'one_first_one_second_each_nonzero_pair':True,'nu':nu,'mu':mu,'controls_note':None if r>=3 else 'Selected-matrix rank only; r2 known variety dimension15 differs from this chamber rank14.'}

def build_certificate():
    require(flint.__version__=='0.9.0','python-flint0.9.0 required')
    return {'method':'literal source monomials; actual integer positive weights; FLINT exact rational matrix rank; explicit integer kernel vectors','finite_tests_not_infinite_proof':True,'floating_point_used':False,'cases':[build_case(r) for r in [1,2,3,4,5,6]]}
from archive_tools import run
if __name__ == "__main__":
    raise SystemExit(run(build_certificate,"evidence/elementary-certificate.json","elementary-two literal full-matrix rational ranks"))
