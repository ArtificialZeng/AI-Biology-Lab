#!/usr/bin/env python3
"""Fresh exact check from serialized reactions, not an owner verifier.

Requires Python 3 and SymPy >= 1.13.3. All identities are checked over rational function fields.
Analytic conclusions use the explicit positivity certificates reported below;
the checker does not certify scientific novelty or runner provenance.
"""
import json

import sympy as S


def zero(e):
    return S.cancel(e) == 0


def main():
    a, b, c, d, p, q = S.symbols('a b c d p q', positive=True)
    x, y, z, t = S.symbols('x y z t')
    rates = (a, b, c, d, p, q)
    # Each pair is (reactant stoichiometry, product stoichiometry).
    reactions = [((0, 0), (1, 0)), ((0, 0), (0, 1)),
                 ((2, 0), (1, 0)), ((0, 2), (0, 1)),
                 ((2, 0), (1, 1)), ((0, 2), (1, 1))]
    columns = [S.Matrix(v) - S.Matrix(u) for u, v in reactions]
    N = S.Matrix.hstack(*columns)
    flux = S.Matrix([k * x**u[0] * y**u[1]
                     for k, (u, _) in zip(rates, reactions)])
    f = N * flux
    A, B = c + p, d + q
    Delta = S.expand(A * B - p * q)
    u = (B * a + q * b) / Delta
    v = (p * a + A * b) / Delta
    normalized = S.Matrix([x*x-u, y*y-v])
    M = S.Matrix([[-A, q], [p, -B]])
    assert all(zero(e) for e in f - M * normalized)
    assert N.rank() == 2 and not N.T.nullspace()
    assert all(any(e != 0 for e in col) for col in columns)
    partition_dimensions = []
    # Fix reaction 0 in the first subset: exactly 31 unordered nonempty splits.
    for mask in range(1, (1 << 6) - 1):
        if not mask & 1:
            continue
        left = [j for j in range(6) if (mask >> j) & 1]
        right = [j for j in range(6) if not ((mask >> j) & 1)]
        intersection = N[:, left].rank() + N[:, right].rank() - N.rank()
        partition_dimensions.append(intersection)
    assert len(partition_dimensions) == 31
    assert min(partition_dimensions) >= 1
    assert max(sum(k) for pair in reactions for k in pair) <= 2

    field = S.QQ.frac_field(*rates)
    bases = {}
    for order in [(x, y), (y, x)]:
        gb = S.groebner(list(f), *order, order='lex', domain=field)
        expected = [normalized[0 if var == x else 1] for var in order]
        assert all(zero(h.as_expr() - e) for h, e in zip(gb.polys, expected))
        bases[str(order)] = [str(h.as_expr()) for h in gb.polys]
    J = f.jacobian([x, y])
    assert zero(J.det() - 4*Delta*x*y)
    assert zero(M.det() - Delta)
    assert Delta == c*d + c*q + d*p
    assert S.simplify(u).is_positive and S.simplify(v).is_positive

    # Independent saturation elimination at generic nonzero U,V.
    U, V = S.symbols('U V', nonzero=True)
    abstract_field = S.QQ.frac_field(U, V)
    sat_gb = S.groebner([x*x-U, y*y-V, 1-z*x*y], z, x, y,
                       order='lex', domain=abstract_field)
    eliminated = [h.as_expr() for h in sat_gb.polys if not h.as_expr().has(z)]
    assert len(eliminated) == 2
    assert all(zero(h-e) for h, e in zip(eliminated, [x*x-U, y*y-V]))
    abstract_gb = S.groebner([x*x-U, y*y-V], x, y,
                            order='lex', domain=abstract_field)
    assert abstract_gb.reduce((x*y)**2 - U*V)[1] == 0

    # Radicality: both univariate factors are square-free; their tensor product
    # splits into four distinct points over the algebraic closure in char 0.
    assert S.gcd(S.Poly(x*x-U, x), S.Poly(2*x, x)).degree() == 0
    assert S.gcd(S.Poly(y*y-V, y), S.Poly(2*y, y)).degree() == 0

    # Polynomial identities underlying the positive-quadrant Lyapunov proof.
    s, h = S.symbols('s h', positive=True)
    W = p*(x-s)**2*(x+2*s)/3 + q*(y-h)**2*(y+2*h)/3
    rr = S.Matrix([x*x-s*s, y*y-h*h])
    assert zero(S.diff(W, x) - p*rr[0])
    assert zero(S.diff(W, y) - q*rr[1])
    DM = S.diag(p, q)*M
    assert DM == DM.T
    assert zero(DM.det() - p*q*Delta)
    dotW = (S.Matrix([S.diff(W, x), S.diff(W, y)]).T * M * rr)[0]
    negative_squares = p*A*(rr[0]-q/A*rr[1])**2 + q*Delta/A*rr[1]**2
    assert zero(dotW + negative_squares)
    assert zero(f[0].subs(x,0) - a-q*y*y)
    assert zero(f[1].subs(y,0) - b-p*x*x)

    # Eliminate t=x+ell*y from the abstract grid and verify both directions of
    # the displayed shape representation, with no use of a square-root branch.
    ell = S.symbols('ell', nonzero=True)
    D0 = U-ell*ell*V
    g = t**4 - 2*(U+ell*ell*V)*t*t + D0*D0
    X = ((3*U+ell*ell*V)*t-t**3)/(2*D0)
    Y = (t**3-(U+3*ell*ell*V)*t)/(2*ell*D0)
    mix_field = S.QQ.frac_field(U,V,ell)
    grid_gb = S.groebner([x*x-U,y*y-V],x,y,domain=mix_field)
    for expr in [g.subs(t,x+ell*y), X.subs(t,x+ell*y)-x,
                 Y.subs(t,x+ell*y)-y]:
        assert grid_gb.reduce(S.cancel(expr))[1] == 0
    for expr in [X*X-U,Y*Y-V,X+ell*Y-t]:
        assert zero(S.rem(S.cancel(expr),g,t))
    discriminant = S.factor(S.discriminant(g,t))
    assert zero(discriminant - 4096*ell**4*U**2*V**2*D0**2)

    result = {
        'schema':'fresh-exact-computational-evidence-v1',
        'ok':True,
        'checker_origin':'reconstructed from reactions; owner checker not read',
        'sympy_version':S.__version__,
        'serialized_reactions':[[list(u0),list(v0)] for u0,v0 in reactions],
        'rate_order':[str(k) for k in rates],
        'stoichiometric_columns':[list(map(int,col)) for col in columns],
        'stoichiometric_rank':2,
        'nontrivial_unordered_partitions_checked':len(partition_dimensions),
        'partition_intersection_dimension_histogram':{
            str(k):partition_dimensions.count(k) for k in sorted(set(partition_dimensions))},
        'ideal_identity_verified':True,
        'native_lex_bases':bases,
        'quotient_standard_monomials':['1','x','y','x*y'],
        'radicality_certificate':'square-free x^2-u and y^2-v; u,v>0',
        'jacobian_determinant':str(S.factor(J.det())),
        'saturation_elimination_verified':True,
        'positivity_certificates':{
            'Delta':str(Delta),'u_numerator':str(S.expand(B*a+q*b)),
            'v_numerator':str(S.expand(p*a+A*b))},
        'lyapunov_factorization_verified':True,
        'lyapunov_expression':str(W),
        'negative_derivative_sum_of_squares':str(negative_squares),
        'global_asymptotic_stability':'analytic compact-sublevel/inward-boundary argument in math.md',
        'mixed_reconstruction_both_directions_verified':True,
        'mixed_eliminant':str(g),
        'mixed_eliminant_discriminant':str(discriminant),
        'mixed_bad_slopes':'ell=0 or ell^2*v=u (u,v>0)',
        'scope':'Exact symbolic identities plus stated analytic argument; no novelty or formal gate certification.'}
    print(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2))


if __name__ == '__main__':
    main()
