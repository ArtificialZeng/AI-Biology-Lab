"""Standard-library exact finite check of the scalar lemma's constraints."""
import json,math

def rank_Q(input_rows):
    rows=[list(x) for x in input_rows]
    if not rows:return 0
    pivot_row=0
    for col in range(len(rows[0])):
        found=next((k for k in range(pivot_row,len(rows)) if rows[k][col]),None)
        if found is None:continue
        rows[pivot_row],rows[found]=rows[found],rows[pivot_row]
        pivot=rows[pivot_row][col]
        for k in range(pivot_row+1,len(rows)):
            factor=rows[k][col]
            if factor:
                changed=[pivot*x-factor*y for x,y in zip(rows[k],rows[pivot_row])]
                divisor=math.gcd(*changed)
                rows[k]=[x//divisor for x in changed] if divisor else changed
        pivot_row+=1
        if pivot_row==len(rows):break
    return pivot_row

def case(r,hamming=True):
    n=1<<r
    order=sorted(range(n),key=lambda x:(x.bit_count(),x)) if hamming else list(range(n))
    equations=[]
    for i,x in enumerate(order):
        for y in order[i+1:]:
            # t_g=F(g)-F(0), from the pair(0,g); eliminate it.
            v=[0]*n;v[y]+=1;v[x]-=1;v[x^y]-=1;v[0]+=1
            equations.append(v)
    rank=rank_Q(equations)
    expected=n-1 if hamming and r>=3 else n-2
    if rank!=expected:raise RuntimeError('Scalar constraint rank mismatch')
    constant_annihilates=all(sum(row)==0 for row in equations)
    topbit=[int(x&(1<<(r-1))!=0) for x in range(n)]
    topbit_annihilates=all(sum(a*b for a,b in zip(row,topbit))==0 for row in equations)
    if not constant_annihilates or topbit_annihilates!=(not hamming or r<3):raise RuntimeError('Boundary scalar witness mismatch')
    return {'r':r,'order':n,'ordering':'hamming_then_integer' if hamming else 'integer_order','constraint_count':len(equations),'constraint_rank_Q':rank,'scalar_solution_dimension':n-rank,'constant_annihilates':constant_annihilates,'highest_bit_annihilates':topbit_annihilates}

def build_certificate():
    return {'method':'exact integer elimination of all oriented-XOR scalar equations with t_g eliminated','finite_diagnostics_only':True,'cases':[case(r) for r in [1,2,3,4,5,6]]+[case(r,False) for r in [3,4]]}
from archive_tools import run
if __name__ == "__main__":
    raise SystemExit(run(build_certificate,"evidence/scalar-certificate.json","elementary-two oriented-XOR scalar constraints"))
