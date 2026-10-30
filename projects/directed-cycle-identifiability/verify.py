#!/usr/bin/env python3
"""Exact candidate certificates; no floating-point arithmetic or external solvers."""
import argparse
import itertools
import json
from pathlib import Path
import sympy as sp

def check(condition, message="exact identity failed"):
    """Execute checks even under python -O."""
    if not condition:
        raise ValueError(message)

def coefficients(n, p, leaks):
    z = sp.Symbol('s')
    k = sp.symbols('k1:' + str(n + 1))
    e = {ell: sp.Symbol('e' + str(ell)) for ell in leaks}
    d = [k[i] + e.get(i + 1, 0) for i in range(n)]
    denominator = sp.expand(sp.prod(z + di for di in d) - sp.prod(k))
    numerator = sp.expand(sp.prod(k[:p - 1]) * sp.prod(z + d[i] for i in range(p, n)))
    c = sp.Matrix(sp.Poly(denominator, z).all_coeffs()[1:] + sp.Poly(numerator, z).all_coeffs())
    return z, k, e, d, denominator, numerator, c


def verify(corruption=None):
    z, k, e, d, denominator, numerator, c = coefficients(3, 2, (2,))
    parameters = list(k) + [e[2]]
    k1, k2, k3, e2 = parameters
    A = sp.Matrix([[-k1, 0, k3], [k1, -k2 - e2, 0], [0, k2, -k3]])
    resolvent_matrix = z * sp.eye(3) - A
    check(sp.expand(resolvent_matrix.det() - denominator) == 0)
    check(sp.expand(resolvent_matrix.adjugate()[1, 0] - numerator) == 0)
    # Independently eliminate x1 and x3 from the three ODEs symbolically:
    # x1=(s+k2+e2)*y/k1, x3=k2*y/(s+k3), (s+k1)*x1-k3*x3=u.
    eliminated_denominator = sp.expand((z+k1)*(z+k2+e2)*(z+k3) - k1*k2*k3)
    eliminated_numerator = k1*(z+k3)
    check(eliminated_denominator == denominator)
    check(sp.expand(eliminated_numerator - numerator) == 0)
    if corruption == "numerator_sign":
        check(sp.expand(-numerator - eliminated_numerator) == 0, "corrupted numerator")
    J = c.jacobian(parameters)
    values = (2, 1, 3, 2) if corruption == "off_hyperplane" else (2, 1, 3, 1)
    point = dict(zip(parameters, map(sp.Integer, values)))
    check(all(value > 0 for value in point.values()))
    check(sp.expand(d[1] - k1).subs(point) == 0)
    selected_rows = [0, 2, 3, 4]
    minor = sp.factor(J[selected_rows, :].det())
    if corruption == "minor_sign":
        minor = -minor
    check(minor == -k1**2*k3)
    Jp = J.subs(point)
    check(Jp.rank() == 4)
    check(minor.subs(point) == -12)
    check(sp.gcd(denominator.subs(point), numerator.subs(point)) == 1)
    # Rational inverse from a subvector: (c2,c0,b1,b0).
    C2, C0, B1, B0 = sp.symbols('C2 C0 B1 B0')
    inverse = [B1, C2-B1-B0/B1-C0/B0, B0/B1, C0/B0]
    inverse_substitution = {C2:c[0], C0:c[2], B1:c[3], B0:c[4]}
    for expression, parameter in zip(inverse, parameters):
        check(sp.cancel(expression.subs(inverse_substitution)-parameter) == 0)
    minors = []
    for rows in itertools.combinations(range(5), 4):
        determinant = sp.factor(J[list(rows), :].det())
        minors.append({'coefficient_rows_zero_based': list(rows), 'polynomial': str(determinant), 'point_value': int(determinant.subs(point))})
    # Corrected interior-leak claim. Test all leak sets of size <=2 for n=3,4,
    # every input1/outputp, and every eligible (a,ell). Not a proof for all n.
    cases = []
    for n in (3, 4):
        for leak_count in (1, 2):
            for leaks in itertools.combinations(range(1, n+1), leak_count):
                for p in range(2, n+1):
                    for ell in leaks:
                        if ell >= p:
                            continue
                        for a in range(1, p):
                            if a in leaks:
                                continue
                            _, ks, es, ds, den, num, coeffs = coefficients(n, p, leaks)
                            pars = list(ks) + list(es.values())
                            v = {ks[a-1]:ks[a-1], ks[ell-1]:-ks[ell-1], es[ell]:ks[ell-1]-ks[a-1]}
                            restriction = {es[ell]:ks[a-1]-ks[ell-1]}
                            for expression in list(coeffs):
                                directional = sum(sp.diff(expression, q)*v.get(q, 0) for q in pars)
                                check(sp.expand(directional.subs(restriction)) == 0)
                            cases.append({'n':n,'p':p,'leaks':list(leaks),'a':a,'ell':ell})
    result = {
        'status':'exact-algebraic-checks-passed',
        'source':'Ahmed et al., SIAM J Appl Dyn Syst 25(1):304-350 (2026), Conjecture4.6',
        'model':{'n':3,'input':[1],'output':[2],'leaks':[2]},
        'parameter_order':[str(q) for q in parameters],
        'positive_rational_point':{str(q):int(v) for q,v in point.items()},
        'hyperplane':'k2 + e2 = k1',
        'coefficient_order':['C2','C1','C0','B1','B0'],
        'coefficient_polynomials':[str(q) for q in c],
        'coefficient_values':[int(q) for q in c.subs(point)],
        'jacobian':[[str(q) for q in row] for row in J.tolist()],
        'jacobian_at_point':[[int(q) for q in row] for row in Jp.tolist()],
        'jacobian_rank_at_point':int(Jp.rank()),
        'selected_rows_zero_based':selected_rows,
        'selected_minor':str(minor),
        'selected_minor_point_value':int(minor.subs(point)),
        'all_maximal_minors':minors,
        'denominator':str(denominator), 'numerator':str(numerator),
        'denominator_numerator_gcd_at_point':'1',
        'rational_inverse_order':['k1','k2','k3','e2'],
        'rational_inverse':[str(expr) for expr in inverse],
        'inverse_domain':'B1 * B0 != 0; holds for every positive k1,k3',
        'construction_checks':['state-space determinant','adjugate numerator','direct ODE elimination','rational inverse composition'],
        'corrected_interior_kernel_checks':{'exact_case_count':len(cases),'scope':'n=3,4; leak sets size1,2; all p,a,ell with a not in Leak and a,ell<p','cases':cases},
        'limitations':['Exact script checks do not establish external mathematical review.','No researcher or journal has confirmed it.','No absolute priority claim: source-version and follow-up search can miss corrections.','Interior-leak general proof is written separately; finite checks do not prove it.']
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='write receipt to a separate file')
    parser.add_argument('--self-test', action='store_true', help='reject three deliberately invalid claims')
    args = parser.parse_args()
    result = verify()
    if args.self_test:
        rejected = []
        for corruption in ('off_hyperplane', 'minor_sign', 'numerator_sign'):
            try:
                verify(corruption)
            except ValueError:
                rejected.append(corruption)
            else:
                raise SystemExit('negative control accepted: ' + corruption)
        result['rejected_negative_controls'] = rejected
    rendered = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(rendered)
    else:
        print(rendered, end='')


if __name__ == '__main__':
    main()
