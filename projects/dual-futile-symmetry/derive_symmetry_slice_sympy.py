#!/usr/bin/env python3
"""Second algebra implementation: exact SymPy derivation of the same slice."""
import json
import sympy as s

a,b,e,f,u,w,v,z = s.symbols("a b e f u w v z", positive=True)
h,j,k = a+e+u+w, 2*b+f+v+z, f+v+z
G = s.Matrix([[0,0,-w,0,z,0],[0,0,0,z,0,-w],[a,0,-h,-e,0,0],
              [-b,-b,-f,-b-f-v-z,-b,0],[-b,-b,0,-b,-b-f-v-z,-f],
              [0,a,0,0,-e,-h]])
P = s.Matrix([[1,0,0,1,0,0],[1,0,0,-1,0,0],[0,1,0,0,1,0],
              [0,0,1,0,0,1],[0,0,1,0,0,-1],[0,1,0,0,-1,0]])
A = s.Matrix([[0,-w,z],[a,-h,-e],[-2*b,-f,-j]])
B = s.Matrix([[0,-w,-z],[a,-h,-e],[0,-f,-k]])
if s.simplify(P.inv()*G*P-s.diag(A,B)) != s.zeros(6):
    raise ValueError("block identity")
result = {"sympy_version": s.__version__, "basis_determinant": str(P.det()), "blocks": {}}
for name, block in [("even",A),("odd",B)]:
    coeff = block.charpoly().all_coeffs()
    H = s.Poly(s.expand(coeff[1]*coeff[2]-coeff[3]), a,b,e,f,u,w,v,z)
    if any(c <= 0 for c in H.coeffs()):
        raise ValueError("nonpositive Routh coefficient")
    result["blocks"][name] = {
        "matrix": str(block), "coefficients": [str(s.factor(c)) for c in coeff],
        "Routh_monomials": len(H.terms()), "min_Routh_coefficient": str(min(H.coeffs()))}
if s.expand(B.charpoly().all_coeffs()[-1]-a*(w*k-f*z)) != 0:
    raise ValueError("odd threshold")
print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))
