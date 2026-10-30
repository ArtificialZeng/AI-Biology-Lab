"""Exact original-number Ito reconstruction and stationary-density identities."""
import sys
sys.dont_write_bytecode = True
import sympy as s
from archive_tools import run, require

def zero(value, message):
    reduced = s.factor(s.cancel(s.together(value)))
    require(reduced == 0, message + ': ' + str(reduced))
    return str(reduced)

def build_certificate():
    require(s.__version__ == '1.13.3', 'Use pinned SymPy1.13.3')
    u,v = s.symbols('R1 R2', positive=True)
    x,q,w = s.symbols('x q omega', positive=True)
    m1,m2,t1,t2,rho = s.symbols('mu1 mu2 sigma1 sigma2 rho', real=True)
    N = u+v
    cov = s.Matrix([[t1**2*u**2,rho*t1*t2*u*v], [rho*t1*t2*u*v,t2**2*v**2]])
    drift = s.Matrix([(m1-w*q)*u+w*(1-q)*v, w*q*u+(m2-w*(1-q))*v])
    coords = (u,v)
    def ito_drift(function):
        return sum(s.diff(function,coords[i])*drift[i] for i in range(2)) + sum(
            s.diff(function,coords[i],coords[j])*cov[i,j]/2 for i in range(2) for j in range(2))
    def variance(function):
        gradient = s.Matrix([s.diff(function,c) for c in coords])
        return (gradient.T*cov*gradient)[0]
    chart = {u:1-x,v:x}
    ratio_drift = s.factor(ito_drift(v/N).subs(chart))
    ratio_variance = s.factor(variance(v/N).subs(chart))
    log_drift = s.factor(ito_drift(s.log(N)).subs(chart))
    V = t1**2-2*rho*t1*t2+t2**2
    A = m2-m1+t1**2-rho*t1*t2
    h = m1*(1-x)+m2*x-(t1**2*(1-x)**2+2*rho*t1*t2*x*(1-x)+t2**2*x**2)/2
    out = {'method':'Ito generator applied to R2/(R1+R2) and log(R1+R2), full covariance',
           'ratio_drift':str(ratio_drift),'ratio_variance':str(ratio_variance),'log_population_drift':str(log_drift)}
    out['checks'] = {
     'general_ratio_drift':zero(ratio_drift-(w*(q-x)+x*(1-x)*(A-V*x)),'ratio drift'),
     'general_ratio_variance':zero(ratio_variance-V*x**2*(1-x)**2,'ratio variance'),
     'instantaneous_log_growth':zero(log_drift-h,'log growth')}
    # Only now specialize the target q to the interior fast-switching optimum A/V.
    fixed = {m2:m1-t1**2+rho*t1*t2+V*q}
    out['checks']['fixed_optimum_ratio_drift'] = zero(ratio_drift.subs(fixed)-(q-x)*(w+V*x*(1-x)), 'fixed-optimum drift')
    out['checks']['quadratic_growth_penalty'] = zero((h-h.subs(x,q)+V*(x-q)**2/2).subs(fixed),'quadratic penalty')
    phi_gap=(m1-t1**2/2)-(m2-t2**2/2)
    out['checks']['ordered_fitness_branch'] = zero(phi_gap.subs(fixed)-V*(s.Rational(1,2)-q),'fitness ordering')
    # Log-density copied as a formula from the complete unmodified author stationary function.
    a,b=w*(1-q),w*q
    alpha=-2*((a-b+m1-m2-rho*t1*t2+t2**2)*s.log(x)+(b-a-m1+m2-rho*t1*t2+t1**2)*s.log(1-x)+b/x-a/(x-1))/V
    zero_current_rhs=(2*ratio_drift-s.diff(ratio_variance,x))/ratio_variance
    out['checks']['author_density_zero_current'] = zero(s.diff(alpha,x)-zero_current_rhs,'author stationary density')
    lam=s.symbols('lambda',positive=True)
    D=(1-2*q)*s.log(x/(1-x))+q/x+(1-q)/(1-x)
    log_base=-2*(1-q)*s.log(x)-2*q*s.log(1-x)
    # Derivative avoids irrelevant additive normalizing constants in symbolic logarithm conventions.
    out['checks']['dimensionless_log_density'] = zero(s.diff(alpha.subs(fixed).subs(w,lam*V)-log_base+2*lam*D,x),'dimensionless density')
    out['checks']['potential_derivative'] = zero(s.diff(D,x)-(x-q)/(x**2*(1-x)**2),'potential derivative')
    # Directly simplify the exact author fitness expression, without evaluating any Julia code.
    author_h=m2-s.Rational(1,2)*(1-x)*(-2*m1+2*m2+(1-x)*t1**2)-rho*x*(1-x)*t1*t2-s.Rational(1,2)*x**2*t2**2
    out['checks']['author_fitness_integrand'] = zero(author_h-h,'author fitness integrand')
    out['dimensionless_stationary_density']='p_lambda(x) proportional to x^(-2(1-q))*(1-x)^(-2q)*exp(-2lambda*D_q(x))'
    out['potential']='D_q(x)=(1-2q)log(x/(1-x))+q/x+(1-q)/(1-x)'
    out['objective']='F(omega)=h(q)-V/2*E_lambda[(X-q)^2], lambda=omega/V'
    out['derivative_discriminator']='d/dlambda E_lambda[(X-q)^2]=-2 Cov_lambda((X-q)^2,D_q(X))'
    out['limits']='q in (0,1/2), lambda>0; covariance sign not established by this script; no finite population extinction or Julia execution'
    return out


if __name__ == '__main__':
    raise SystemExit(run(build_certificate, 'evidence/model-certificate.json',
                        'original-number symbolic Ito reconstruction'))
