import json,sys,time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from model import finite_rhs, finite_jac, limit_rhs, limit_jac, physical_species

OUT=Path(__file__).resolve().parents[1]/'evidence'

def run(K,method='LSODA',rtol=1e-9,atol=1e-11,t_end=1000):
    if np.isinf(K):
        f=limit_rhs;jac=limit_jac; y0=np.array([10.,10.,1.,1.,.1]);wi=2
    else:
        f=lambda t,y:finite_rhs(t,y,K)
        jac=lambda t,y:finite_jac(t,y,K)
        y0=np.array([39.,0.,0.,1.,0.,0.]);wi=3
    def peak(t,y):return f(t,y)[wi]
    peak.direction=-1
    def trough(t,y):return f(t,y)[wi]
    trough.direction=1
    start=time.time()
    sol=solve_ivp(f,(0,t_end),y0,method=method,jac=jac,rtol=rtol,atol=atol,events=[peak,trough],dense_output=True,max_step=.5)
    tm=sol.t_events[0];ym=sol.y_events[0];tn=sol.t_events[1];yn=sol.y_events[1]
    keep=tm>t_end*.6;tm=tm[keep];ym=ym[keep]
    kn=tn>t_end*.6;tn=tn[kn];yn=yn[kn]
    periods=np.diff(tm)
    result={'K':str(K),'method':method,'rtol':rtol,'atol':atol,'success':sol.success,'message':sol.message,
      'seconds':time.time()-start,'peak_times':tm[-10:].tolist(),'periods':periods[-10:].tolist(),
      'period_mean':float(np.mean(periods[-10:])) if len(periods) else None,
      'peak_values':ym[-10:,wi].tolist() if len(ym) else [],'trough_values':yn[-10:,wi].tolist() if len(yn) else [],
      'last_peak_state':ym[-1].tolist() if len(ym) else [],'last_state':sol.y[:,-1].tolist(),
      'max_rhs_at_end':float(np.max(np.abs(f(t_end,sol.y[:,-1]))))}
    if not np.isinf(K):
        xx=np.stack([physical_species(y,K) for y in sol.y.T]);result['min_species_during_integration']=float(xx.min())
    np.savez(OUT/f'trajectory_K{K}_{method}.npz',t=sol.t,y=sol.y)
    (OUT/f'baseline_K{K}_{method}.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result),flush=True)
    return result

if __name__=='__main__':
    vals=[float(x) for x in sys.argv[1:]] or [100.0,1000.0,10000.0,np.inf]
    all_results=[run(k) for k in vals]
    (OUT/'baseline.json').write_text(json.dumps(all_results,indent=2))
