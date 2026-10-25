"""Independent floating-point check built from the published reaction network.

No imports from the primary implementation. This is numerical evidence, not proof.
Source: arXiv:1809.02886v2, reactions (1), species Table 3, rates Table 1.
"""
from __future__ import annotations
import argparse
import json
import time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root, minimize_scalar

RATES = np.array([1, 1, 1, 1, 100, 1, .9, 3, 1, 100.])
NAMES = ('S0','K','S0K','S1K','S2','P','S2P','S1','S1P')
# Each column is independently entered from the chemical arrows, not copied ODEs.
ARROWS = [({'S0':1,'K':1},{'S0K':1}),
          ({'S0K':1},{'S0':1,'K':1}),
          ({'S0K':1},{'S1K':1}),
          ({'S1K':1},{'S2':1,'K':1}),
          ({'S2':1,'P':1},{'S2P':1}),
          ({'S2P':1},{'S2':1,'P':1}),
          ({'S2P':1},{'S1':1,'P':1}),
          ({'S1':1,'P':1},{'S1P':1}),
          ({'S1P':1},{'S1':1,'P':1}),
          ({'S1P':1},{'S0':1,'P':1})]
N = np.zeros((9,10))
POW = np.zeros((10,9), dtype=int)
for j,(before,after) in enumerate(ARROWS):
    for name,v in before.items():
        N[NAMES.index(name),j] -= v
        POW[j,NAMES.index(name)] += v
    for name,v in after.items(): N[NAMES.index(name),j] += v
CONS = np.array([[0,1,1,1,0,0,0,0,0],
                 [0,0,0,0,0,1,1,0,1],
                 [1,0,1,1,1,0,1,1,1]],dtype=float)
assert np.array_equal(CONS @ N, np.zeros((3,10)))

def full_rhs(t,x):
    return N @ (RATES * np.prod(x[None,:] ** POW, axis=1))

def cs_jac(fun,y):
    out=np.empty((len(y),len(y)))
    for j in range(len(y)):
        z=np.asarray(y,dtype=complex); z[j]+=1e-28j
        out[:,j]=np.imag(fun(0,z))/1e-28
    return out

def full_jac(t,x): return cs_jac(full_rhs,x)

# Independent full model reduction used only for shooting/Floquet, not the main
# 9D simulations. Selected coordinates S0K,S1K,S2,S2P,S1,S1P.
SEL=np.array([2,3,4,6,7,8])
def reconstruct(z,K):
    x=np.zeros(9,dtype=np.result_type(z));x[SEL]=z
    x[0]=40-np.sum(z);x[1]=K-z[0]-z[1];x[5]=8-z[3]-z[5]
    return x
def reduced_rhs(K):
    return lambda t,z: full_rhs(t,reconstruct(z,K))[SEL]

# Infinite kinase: free S0 instantly rebinds after either dissociation or final
# dephosphorylation. Slow states (S0K,S1K,S2,S2P,S1,S1P).
def slow_six(t,z):
    a,b,c,d,e,f=z;p=8-d-f
    bind2=100*c*p;bind1=3*e*p
    return np.array([100*f-a,a-b,b-bind2+d,
                     bind2-1.9*d,.9*d-bind1+f,bind1-101*f])
def slow_reconstruct(y):
    # Independent 5D coordinates (S0K,S1K,S2P,S1,S1P).
    a,b,d,e,f=y
    return np.array([a,b,40-a-b-d-e-f,d,e,f])
def slow_five(t,y): return slow_six(t,slow_reconstruct(y))[[0,1,3,4,5]]

def integrate_and_measure(K,method='Radau',alt=False,duration=800.,tol=1e-10):
    if K is None:
        fun=slow_six
        y0=np.array([4.,3.,17.,2.,13.,1.]) if alt else np.array([40.,0.,0.,0.,0.,0.])
        idx=2
        expected=np.array([40.])
    else:
        fun=full_rhs
        y0=np.array([2.,K-5.,3.,2.,12.,5.,2.,18.,1.]) if alt else np.array([40.,K,0.,0.,0.,8.,0.,0.,0.])
        idx=4
        expected=np.array([K,8.,40.])
    def maxima(t,y):return fun(t,y)[idx]
    maxima.direction=-1
    def minima(t,y):return fun(t,y)[idx]
    minima.direction=1
    jac=lambda t,y:cs_jac(fun,y)
    started=time.monotonic()
    sol=solve_ivp(fun,(0.,duration),y0,method=method,rtol=tol,atol=tol*.01,
                  jac=jac if method in ('Radau','BDF','LSODA') else None,
                  events=(maxima,minima),dense_output=True)
    peaks=sol.t_events[0];peaks=peaks[peaks>duration*.7]
    troughs=sol.t_events[1];troughs=troughs[troughs>duration*.7]
    yy=sol.sol(np.linspace(0,duration,20001))
    cons=(CONS@yy) if K is not None else np.sum(yy,axis=0,keepdims=True)
    info={'K':K,'method':method,'initial_condition':'interior' if alt else ('all_S0K' if K is None else 'all_S0'),
          'duration':duration,'rtol':tol,'atol':tol*.01,'success':bool(sol.success),
          'message':sol.message,'wall_seconds':time.monotonic()-started,
          'function_evaluations':int(sol.nfev),'jacobian_evaluations':int(sol.njev),
          'max_conservation_error':float(np.max(np.abs(cons-expected[:,None]))),
          'minimum_sampled_species':float(np.min(yy)),
          'late_peak_count':len(peaks),'late_peak_times':peaks.tolist(),
          'last_period':float(peaks[-1]-peaks[-2]) if len(peaks)>1 else None,
          'late_period_spread':float(np.ptp(np.diff(peaks))) if len(peaks)>2 else None,
          'last_S2_max':float(sol.sol(peaks[-1])[idx]) if len(peaks) else None,
          'last_S2_min':float(sol.sol(troughs[-1])[idx]) if len(troughs) else None}
    if info['last_S2_max'] is not None and info['last_S2_min'] is not None:
        info['last_S2_peak_to_peak']=info['last_S2_max']-info['last_S2_min']
    if K is None:
        pfree=8-yy[3]-yy[5]
        info['minimum_sampled_free_phosphatase']=float(np.min(pfree))
    if len(peaks)>1:
        point=sol.sol(peaks[-2]+.23*(peaks[-1]-peaks[-2]))
        info['_shoot_guess']=(point[SEL] if K is not None else point[[0,1,3,4,5]]).tolist()
    return info

def shooting(K,guess,period,tol=2e-11):
    fun=reduced_rhs(K) if K is not None else slow_five
    y0=np.asarray(guess);n=len(y0);phase=fun(0,y0);phase/=np.linalg.norm(phase)
    method='Radau' if K is not None else 'DOP853'
    def flow(y,T):
        return solve_ivp(fun,(0,T),y,method=method,rtol=tol,atol=tol*.01,
                         jac=(lambda t,v:cs_jac(fun,v)) if method=='Radau' else None).y[:,-1]
    def residual(v):
        y=v[:-1];T=v[-1]
        return np.r_[flow(y,T)-y, np.dot(y-y0,phase)]
    fit=root(residual,np.r_[y0,period],method='hybr',options={'xtol':2e-9})
    yy=fit.x[:-1];T=fit.x[-1]
    def augmented(t,v):
        y=v[:n];m=v[n:].reshape(n,n)
        return np.r_[fun(t,y),(cs_jac(fun,y)@m).ravel()]
    v0=np.r_[yy,np.eye(n).ravel()]
    sol=solve_ivp(augmented,(0,T),v0,method=method,rtol=tol,atol=tol*.01,
                  dense_output=True)
    mon=sol.y[n:,-1].reshape(n,n);eig=np.linalg.eigvals(mon)
    samples=sol.sol(np.linspace(0,T,5001))[:n]
    sx=samples[2] if K is not None else np.array([slow_reconstruct(samples[:,j])[2] for j in range(samples.shape[1])])
    def s2(t):
        y=sol.sol(t)[:n]
        return y[2] if K is not None else slow_reconstruct(y)[2]
    grid=np.linspace(0,T,5001);imax=int(np.argmax(sx));imin=int(np.argmin(sx))
    def extremum(i,sgn):
        lo=grid[max(0,i-2)];hi=grid[min(len(grid)-1,i+2)]
        z=minimize_scalar(lambda t:sgn*s2(t),bounds=(lo,hi),method='bounded',options={'xatol':1e-12})
        return float(s2(z.x)),float(z.x)
    mx,tmx=extremum(imax,-1);mn,tmn=extremum(imin,1)
    neutral=int(np.argmin(abs(eig-1)))
    info={'K':K,'kind':'floating_point_periodic_shooting','method':method,
          'root_success':bool(fit.success),'root_message':str(fit.message),
          'period':float(T),'shooting_residual_infinity_norm':float(np.linalg.norm(residual(fit.x),ord=np.inf)),
          'independently_integrated_closure_norm':float(np.linalg.norm(sol.y[:n,-1]-yy,ord=np.inf)),
          'reduced_coordinate_initial_state':yy.tolist(),'S2_max':mx,'S2_min':mn,'S2_peak_to_peak':mx-mn,
          'S2_max_time':tmx,'S2_min_time':tmn,'floquet_multipliers':[{'real':float(z.real),'imag':float(z.imag),'abs':float(abs(z))} for z in eig],
          'neutral_multiplier_error':float(abs(eig[neutral]-1)),
          'largest_nontrivial_multiplier_abs':float(max(abs(np.delete(eig,neutral)))),
          'monodromy_determinant':float(np.linalg.det(mon)),
          'rtol':tol,'atol':tol*.01,'proof_status':'floating_point_only'}
    return info

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--quick',action='store_true');ap.add_argument('--no-shoot',action='store_true');args=ap.parse_args()
    results=[];shoots=[]
    cases=[(K,'Radau',False,1e-10) for K in [100.,1000.,10000.,None]]
    if not args.quick:
        cases += [(K,'BDF',True,2e-11) for K in [100.,1000.,10000.,None]]
    for K,method,alt,tol in cases:
        info=integrate_and_measure(K,method,alt,duration=400.,tol=tol)
        guess=info.pop('_shoot_guess',None);results.append(info)
        print(json.dumps(info),flush=True)
        if guess is not None and not args.no_shoot and method=='Radau':
            shot=shooting(K,guess,info['last_period']);shoots.append(shot);print(json.dumps(shot),flush=True)
        out={'scope':'Independent implementation from reaction arrows; all numerical claims are floating point, not proof.',
             'rates':RATES.tolist(),'totals':{'P':8.,'S':40.},'stoichiometric_matrix':N.tolist(),
             'conservation_matrix':CONS.tolist(),'trajectory_checks':results,'periodic_shooting':shoots}
        Path('evidence/independent_numerics.json').write_text(json.dumps(out,indent=2))

if __name__=='__main__':main()
