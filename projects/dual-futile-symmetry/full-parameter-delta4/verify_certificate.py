"""Solver/cache independent exact verification of det>=0 => Delta4>0.

Run: python -B verify_certificate.py
Only python-flint and the standard library are needed. Coefficients c_k are
rebuilt as principal-minor sums, and Delta4 as a 4x4 determinant, independently
of multiplier.py's 6x6 characteristic determinant and closed Delta4 formula.
"""
from pathlib import Path
from itertools import permutations,combinations
from hashlib import sha256
import json,sys
sys.dont_write_bytecode = True
from flint import fmpz_mpoly_ctx
from integrity import verify_manifest
H=Path(__file__).resolve().parent
verify_manifest(H)
from matrix import matrix,NAMES
def require(truth,message):
    if not truth:
        raise RuntimeError(message)
ctx=fmpz_mpoly_ctx.get(NAMES);zero=ctx.constant(0);one=ctx.constant(1)
def det(m):
    n=len(m);answer=zero
    for p in permutations(range(n)):
        term=one
        for i,j in enumerate(p):
            term*=m[i][j]
            if term.is_zero():break
        if term.is_zero():continue
        sign=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
        answer+=sign*term
    return answer
g=matrix(ctx.gens());cs=[one]
for k in range(1,7):
    cs.append(sum((det([[-g[i][j] for j in inds] for i in inds])
                   for inds in combinations(range(6),k)),zero))
one,c1,c2,c3,c4,c5,c6=cs
h4=[[cs[2*j-i+1] if 0<=2*j-i+1<=6 else zero
     for j in range(4)] for i in range(4)]
d4=det(h4)
data=json.loads((H/'certificate.json').read_text())
require(data['names']==list(NAMES) and data['denominator']==1,'Certificate convention mismatch')
require(all(int(v)>0 and sum(e)==4 and len(e)==16 and all(isinstance(x,int) and x>=0 for x in e) for e,v in data['Q']),'Malformed multiplier')
require(len({tuple(e) for e,v in data['Q']})==len(data['Q']),'Duplicate multiplier monomial')
q=ctx.from_dict({tuple(e):int(v) for e,v in data['Q']})
a,b4,b10,c,e1,e4,f7,f10,u2,u3,v5,v6,x8,x9,y11,y12=ctx.gens()
A=e4+v5+v6;B=f10+y11+y12;L=a+c+u2+u3+x8+x9
compact=(L+f7)*A**3+(e1*(v5+v6)+b4*v6)*A**2+\
        (L+e1)*B**3+(f7*(y11+y12)+b10*y12)*B**2+A*B*(A+B)**2
require(q==compact,'Multiplier differs from compact formula')
core=a*c*e4*f10*v6*y12
negative=ctx.from_dict({tuple(e):-int(v) for e,v in d4.terms() if v<0})
require(negative==core*q,'Negative-part identity fails')
r=d4-q*c6
require(all(v>0 for v in r.coeffs()),'Nonpositive remainder coefficient')
require(len(r)>0,'A nonzero positive remainder is needed for strictness')
witness=min(((tuple(map(int,e)),int(v)) for e,v in r.terms()),
            key=lambda ev:(sum(x!=0 for x in ev[0]),ev[0]))
require(r[witness[0]]==witness[1]>0,'Strict coefficient witness fails')
def stats(p):
    return {'terms':len(p),'degree':int(p.total_degree()),
            'negative_coefficients':sum(v<0 for v in p.coeffs()),
            'coefficient_sum':str(sum(p.coeffs()))}
receipt={'arithmetic':'FLINT arbitrary precision integers; no floating arithmetic',
 'construction':'c_k via principal-minor sums of -G; Delta4 via 4x4 Leibniz determinant',
 'names':NAMES,'characteristic_coefficients':[stats(p) for p in cs],
 'Delta4':stats(d4),'c6':stats(c6),'Q':stats(q),'R':stats(r),
 'identities':{'Q_equals_compact_formula':True,'negative_part_Delta4_equals_core_times_Q':True,
               'Delta4_equals_Q_times_c6_plus_R':True},
 'positive_remainder_witness':{'powers':list(witness[0]),'coefficient':witness[1]},
 'certificate_sha256':sha256((H/'certificate.json').read_bytes()).hexdigest(),
 'integrity_verified':True,'optimization_mode_preserves_checks':True,'verified':True}
print(json.dumps(receipt,indent=2),flush=True)
