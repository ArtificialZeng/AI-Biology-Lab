"""Exact finite diagnostics of the literal cyclic-inversion 3-sunlet map.

The uniform theorem is proved in inversion_proof.md, not by this finite list.
Run with CPython3.13 and python-flint0.9.0; stdout is deterministic JSON.
"""
import hashlib
import json
import random
import flint
from flint import fmpq_mat


def require(condition, reason):
    if not condition:
        raise RuntimeError(reason)


def distance(x, q):
    return min(x % q, (-x) % q)


def literal_terms(g, h, q):
    k=q//2+1
    a,b,c=distance(g,q),distance(h,q),distance(g+h,q)
    return ([(0,a),(1,b),(2,c),(3,c),(4,a)],
            [(0,a),(1,b),(2,c),(3,b),(5,a)])


def kernel_vectors(k):
    basis=[]
    for a in range(k):
        row=[0]*(6*k)
        row[a]=1;row[5*k+a]=-1
        if a:row[4*k+a]=-1
        basis.append(row)
    for block in (1,2,3):
        row=[0]*(6*k)
        for a in range(k):
            row[block*k+a]=1;row[5*k+a]=-1
            if a:row[4*k+a]=-1
        basis.append(row)
    row=[0]*(6*k);row[4*k]=1;basis.append(row)
    return basis


def build_case(q):
    k=q//2+1
    mu=[-1]+[1]*(k-1)
    weights=[[1]*k for _ in range(3)]+[[4*a+1 for a in range(k)],[2]*k,[2+x for x in mu]]
    rows=[[0]*(q*q) for _ in range(6*k)]
    serialized=[];choices=[];margins=[];triples=set()
    for g in range(q):
        for h in range(q):
            first,second=literal_terms(g,h,q)
            a,b,c=distance(g,q),distance(h,q),distance(g+h,q)
            values=[sum(weights[edge][label] for edge,label in term) for term in (first,second)]
            difference=values[0]-values[1]
            require(difference==4*(c-b)-mu[a],'literal monomial weight difference')
            require(difference!=0,'strict weights')
            selected=0 if difference<0 else 1
            require((selected==1) if g==0 else ((selected==0)==(c<=b)),'exact branch rule')
            for edge,label in (first,second)[selected]:rows[edge*k+label][g*q+h]+=1
            serialized.append([g,h,first,second,values,selected])
            choices.append(selected);margins.append(abs(difference));triples.add((a,b,c))
    require(all(v>0 for block in weights for v in block),'actual weights positive')
    for a in range(1,k):
        require(choices[a*q]==1 and choices[a*q+(-a)%q]==0,'every nonzero orbit has both branches')
    matrix=fmpq_mat(rows)
    basis=fmpq_mat(kernel_vectors(k))
    require(basis*matrix==fmpq_mat(k+4,q*q),'explicit integer kernel annihilates')
    require(basis.rank()==k+4,'explicit kernel independent')
    rank=matrix.rank()
    expected=5*k-4 if q>=5 else (5 if q==3 else 10)
    require(rank==expected,'exact selected rank')
    # Separate ordinary Jacobian evaluations of the sum of BOTH literal terms.
    rng=random.Random(13013+q)
    jacobian_points=[];jacobian_ranks=[]
    for test in range(3):
        point=[rng.randint(1,9) for _ in range(6*k)]
        jac=[[0]*(6*k) for _ in range(q*q)]
        for g in range(q):
            for h in range(q):
                for term in literal_terms(g,h,q):
                    indices=[edge*k+label for edge,label in term]
                    for i,index in enumerate(indices):
                        product=1
                        for j,other in enumerate(indices):
                            if j!=i:product*=point[other]
                        jac[g*q+h][index]+=product
        jacobian_points.append(point);jacobian_ranks.append(fmpq_mat(jac).rank())
    require(max(jacobian_ranks)==expected,'ordinary Jacobian reaches the finite-case dimension')
    return {'q':q,'orbit_count_k':k,'raw_consistent_label_coordinates':q*q,
            'distinct_orbit_triples':len(triples),'shape':[6*k,q*q],
            'actual_positive_weights':weights,'minimum_strict_margin':min(margins),
            'first_columns':choices.count(0),'second_columns':choices.count(1),
            'rank_Q':rank,'left_nullity':6*k-rank,'kernel_basis_rank_Q':basis.rank(),
            'upper_bound_5k_minus4':5*k-4,'ambient_after_duplicate_deletion':len(triples),
            'inside_theorem':q>=5,'held_out_from_triage':q in (11,12,17,18),
            'literal_inputs_sha256':hashlib.sha256(json.dumps(serialized,separators=(',',':')).encode()).hexdigest(),
            'jacobian_integer_points':jacobian_points,'jacobian_rank_Q':jacobian_ranks,
            'control_note':None if q>=5 else 'Outside theorem; ambient cap below 5k-4. No new small-order claim.'}


def build_certificate():
    require(flint.__version__=='0.9.0','python-flint0.9.0 required')
    return {'method':'literal raw q^2 monomials, actual positive orbit-constant weights, FLINT exact rational ranks, explicit integer kernel and ordinary Jacobians',
            'finite_diagnostics_not_uniform_proof':True,'floating_point_used':False,
            'cases':[build_case(q) for q in (3,4,5,6,7,8,9,10,11,12,17,18)]}


from archive_tools import run
if __name__ == "__main__":
    raise SystemExit(run(build_certificate,'evidence/inversion-certificate.json','cyclic-inversion literal matrices, exact ranks and integer kernel'))
