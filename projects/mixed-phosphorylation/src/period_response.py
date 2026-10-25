"""Numerical periodic-orbit shooting and first-order slow-flow response.

All results are floating-point evidence, not interval certificates.
"""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from model import limit_rhs,limit_jac

OUT=Path(__file__).resolve().parents[1]/'evidence'
base=json.loads((OUT/'baseline_Kinf_LSODA.json').read_text())
last=None

def flow(y,T,forcing=False):
    n=5
    def aug(t,z):
        u=z[:n];J=limit_jac(t,u)
        out=[limit_rhs(t,u), (J@z[n:n+n*n].reshape(n,n)).ravel()]
        if forcing:
            a,b,w,c,d=u;p=8-c-d;f=limit_rhs(t,u)
            h0=a+100*d
            f1=np.array([-f[0]-100*f[4],0,0,0,-3*h0*p])
            out.append(J@z[n+n*n:]+f1)
        return np.concatenate(out)
    z=np.r_[y,np.eye(n).ravel(),np.zeros(n)] if forcing else np.r_[y,np.eye(n).ravel()]
    return solve_ivp(aug,(0,T),z,method='LSODA',rtol=2e-12,atol=2e-14,dense_output=True)

def evaluate(z):
    global last
    if last is not None and np.array_equal(z,last[0]):return last[1:]
    y=z[:5];T=z[5];sol=flow(y,T);end=sol.y[:,-1]
    f=np.r_[end[:5]-y,limit_rhs(0,y)[2]]
    jac=np.zeros((6,6));jac[:5,:5]=end[5:].reshape(5,5)-np.eye(5);jac[:5,5]=limit_rhs(T,end[:5]);jac[5,:5]=limit_jac(0,y)[2]
    last=(z.copy(),f,jac,sol)
    return f,jac,sol

z0=np.r_[base['last_peak_state'],base['period_mean']]
rr=root(lambda z:evaluate(z)[0],z0,jac=lambda z:evaluate(z)[1],tol=1e-10)
yy=rr.x[:5];T=rr.x[5];s=flow(yy,T,True);end=s.y[:,-1];M=end[5:30].reshape(5,5);v=end[30:]
B=np.zeros((6,6));B[:5,:5]=M-np.eye(5);B[:5,5]=limit_rhs(0,yy);B[5,:5]=limit_jac(0,yy)[2]
q=np.linalg.solve(B,np.r_[-v,0.])
# Obtain the first-order change at the minimum; event-time shift vanishes
# because the observable derivative is zero at an extremum.
from scipy.optimize import minimize_scalar
mn=minimize_scalar(lambda t:s.sol(t)[2],bounds=(T*.05,T*.95),method='bounded',options={'xatol':1e-13})
tm=mn.x;za=s.sol(tm);qtm=za[5:30].reshape(5,5)@q[:5]+za[30:]
eig=np.linalg.eigvals(M)
result={'status':'floating-point numerical evidence, not rigorous enclosure','shoot_success':bool(rr.success),'shoot_residual':float(np.max(np.abs(evaluate(rr.x)[0]))),'T_limit':float(T),'peak_state':yy.tolist(),'w_max':float(yy[2]),'w_min':float(mn.fun),'minimum_time':float(tm),'monodromy_eigenvalues':[[float(x.real),float(x.imag)] for x in eig],'period_first_order_coefficient':float(q[5]),'maximum_first_order_coefficient':float(q[2]),'minimum_first_order_coefficient':float(qtm[2]),'response_BVP_residual':float(np.max(np.abs(B@q+np.r_[v,0.]))),'phase_condition':'d(S2)/dt=0 at maximum','slow_flow_first_order':'F1=(-F0_u-100F0_b,0,0,0,-3*(u+100b)*(8-a-b))'}
(OUT/'period_response.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
np.savez(OUT/'limit_shooting.npz',T=T,y0=yy,t=s.t,z=s.y,q0=q[:5],T1=q[5])
