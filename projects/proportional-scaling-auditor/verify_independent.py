"""Separate construction from unnormalized population equations and elimination."""
import json
import sys
sys.dont_write_bytecode = True
import sympy as sp
from archive_tools import run, require


def build_certificate():
    require(sp.__version__ == '1.13.3', 'Use pinned SymPy1.13.3')
    S, X, Y, n, u, v = sp.symbols('S X Y n u v')
    total = S + X + Y
    force = (24*X + Y)/total
    raw = sp.Matrix([sp.Rational(34,5)*total - 5*S - S*force,
                     S*force/10 - X, 9*S*force/10 - 17*Y])
    change = {S:n*(1-u-v), X:n*u, Y:n*v}
    quotient = sp.Matrix([sp.cancel(((raw[1]*total-X*sum(raw))/total**2).subs(change)),
                          sp.cancel(((raw[2]*total-Y*sum(raw))/total**2).subs(change))])
    f, h = [sp.Poly(sp.expand(10*z),u,v).as_expr() for z in quotient]
    require(not any(z.has(n) for z in quotient), 'Fractions still depend on population size')
    resultant = sp.factor(sp.resultant(f,h,u))
    q = 41464*v**2-72863*v+9180
    require(sp.expand(resultant + 800*v*(10*v-1)*q) == 0, 'Independent elimination identity failed')
    J = quotient.jacobian([u,v])
    tests = {}
    for label, point, det, trace in [
        ('disease_free',{u:0,v:0},sp.Integer(5),-sp.Rational(183,10)),
        ('positive_saddle',{u:sp.Rational(1,10),v:sp.Rational(1,10)},-sp.Rational(211,125),-sp.Rational(953,50))]:
        require(quotient.subs(point) == sp.zeros(2,1), 'Specified point is not stationary')
        M = J.subs(point)
        require(M.det()==det and sp.trace(M)==trace, 'Independent Jacobian test failed')
        tests[label] = {'Jacobian':[[str(x) for x in row] for row in M.tolist()],
                        'determinant':str(det),'trace':str(trace)}
    lo, hi = sp.Rational(13,100), sp.Rational(7,50)
    # Polynomial derivative is negative throughout this bracket: its maximum is at hi.
    require(q.subs(v,lo)>0 and q.subs(v,hi)<0 and sp.diff(q,v).subs(v,hi)<0,
            'Algebraic root is not uniquely isolated')
    ua = sp.Rational(125,754)-sp.Rational(37,377)*v
    require(all(sp.rem(sp.together(z.subs(u,ua)),q,v)==0 for z in quotient),
            'Eliminated positive equilibrium does not solve both equations')
    trace = sp.rem(sp.together(sp.trace(J).subs(u,ua)),q,v)
    det = sp.rem(sp.together(J.det().subs(u,ua)),q,v)
    require(sp.degree(trace,v)<=1 and sp.diff(trace,v)>0 and trace.subs(v,hi)<0,
            'Positive equilibrium trace is not strictly negative on bracket')
    require(sp.degree(det,v)<=1 and sp.diff(det,v)>0 and det.subs(v,lo)>0,
            'Positive equilibrium determinant is not strictly positive on bracket')
    require(ua.subs(v,hi)>0 and (1-ua-v).subs(v,hi)>0,
            'Isolated equilibrium leaves physical triangle')
    growth = sp.cancel((sum(raw)/total).subs(change))
    require(growth.subs({u:sp.Rational(1,10),v:sp.Rational(1,10)})==1,
            'Balanced-growth interpretation failed')
    return {'schema':'independent-proportional-scaling-v1','arithmetic':'exact QQ polynomial identities and rational interval bounds',
            'normalized_rhs':[str(z) for z in quotient], 'elimination_resultant':str(resultant),
            'growth_rate_expression':str(growth), 'rational_tests':tests,
            'positive_sink':{'isolating_interval':[str(lo),str(hi)],'polynomial':str(q),
                             'u_expression':str(ua),'trace':str(trace),'determinant':str(det),
                             'trace_upper_bound':str(trace.subs(v,hi)),
                             'determinant_lower_bound':str(det.subs(v,lo))},
            'scope':'counterexample to the named unqualified proportional-scaling extension only; all stability is in fractions',
            'mortality_limitation':'mu1<muS; no claim for all mu_j>=muS or a calibrated disease mechanism'}


if __name__ == '__main__':
    raise SystemExit(run(build_certificate, 'evidence/independent-certificate.json', 'independent numbers-equation elimination'))
