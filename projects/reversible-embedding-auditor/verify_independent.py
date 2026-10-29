"""Independent literal-input/spectral-polynomial certificate; no constructed U imported."""
import sys
sys.dont_write_bytecode = True
import sympy as s
from archive_tools import run, require


def build_certificate():
    require(s.__version__ == '1.13.3', 'Use pinned SymPy1.13.3')
    numerator=s.Matrix([[67239937,67152554,67152554,66890411],
                        [67152554,67123429,67123427,67036046],
                        [67152554,67123427,67123429,67036046],
                        [66890411,67036046,67036046,67472953]])
    M=numerator/s.Integer(268435456)
    A=s.Matrix([[-63,33,33,-3],[33,-79,29,17],[33,29,-79,17],[-3,17,17,-31]])/36
    B=s.Matrix([[9,-21,9,3],[9,-5,-7,3],[-21,33,-5,-7],[3,-7,3,1]])/12
    one=s.ones(4,1);I=s.eye(4);zero=s.zeros(4)
    require(M==M.T and M*one==one and all(x>0 for x in M), 'Rational M is not positive symmetric stochastic')
    lam=s.Rational(1,512);mu=s.Rational(1,512**3)
    E=s.ones(4)/4
    R=s.simplify((M-mu*I-(1-mu)*E)/(lam-mu))
    P=I-E-R
    for name,F,rank in [('E',E,1),('R',R,1),('P',P,2)]:
        require(F==F.T and F*F==F and s.trace(F)==rank, 'Spectral projector failed: '+name)
    require(E*R==R*E==E*P==P*E==R*P==P*R==zero, 'Spectral subspaces not orthogonal')
    require(M==E+lam*R+mu*P, 'Spectral reconstruction failed')
    z=s.Symbol('z')
    expected=(z-1)*(z-lam)*(z-mu)**2
    require(s.expand(M.charpoly(z).as_expr()-expected)==0, 'Characteristic polynomial mismatch')
    require(A==-R-3*P and B*B==-4*P, 'Literal generator does not match spectral algebra')
    require(A*B==B*A and B*E==E*B==B*R==R*B==zero and B*P==P*B==B,
            'Commuting exponential decomposition failed')
    require(A*one==B*one==s.zeros(4,1) and one.T*A==one.T*B==s.zeros(1,4), 'Zero-sum generator condition failed')
    require(B!=B.T, 'Nonreversible generator distinction lost')
    bounds={}
    for i in range(4):
        for j in range(4):
            if i==j:continue
            ca,cb=A[i,j],B[i,j]
            lower=ca*(6 if ca>=0 else 9)+cb*(3 if cb>=0 else s.Rational(22,7))
            require(lower>=0 and (ca!=0 or cb!=0), 'Strict affine off-diagonal bound failed')
            bounds[f'{i+1},{j+1}']=str(lower)
    require(A[0,3]==-s.Rational(1,12) and 9*A[0,3]==-s.Rational(3,4), 'Negative principal-log entry missing')
    return {'schema':'independent-reversible-embedding-v1',
            'arithmetic':'exact rational matrices; no numerical logarithm, exponential or eigensolver',
            'input_denominator':268435456,'spectrum':['1','2^-9','2^-27','2^-27'],
            'projector_ranks':{'E':1,'R':1,'P':2},
            'characteristic_polynomial':str(s.factor(expected)),
            'A_from_input':[[str(x) for x in row] for row in A.tolist()],
            'B_literal':[[str(x) for x in row] for row in B.tolist()],
            'generator':'Q=9*log(2)*A+pi*B; off-diagonals strictly positive',
            'offdiagonal_interval_infima':bounds,
            'principal_log_14':'-3*log(2)/4',
            'checked_exponential_prerequisites':'B^2=-4P; BE=BR=0; BP=B; A=-R-3P; AB=BA',
            'analytic_step':'exp(pi*B)=I and exp(9log(2)*A)=E+2^-9R+2^-27P; see proof',
            'scope':'one repeated-positive-spectrum reversible matrix; specified Theorem3.16 implication only',
            'limitations':['General nonprincipal embeddings are prior work.','No failure claim for algorithms restricted to reversible generators.','No open-set, generic-spectrum, or complete-classification claim.']}


if __name__=='__main__':
    raise SystemExit(run(build_certificate,'evidence/independent-certificate.json','literal-input spectral-polynomial check'))
