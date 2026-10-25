"""Scientific figures from separately recorded numerical results."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from model import finite_rhs,finite_jac,limit_rhs,limit_jac

ROOT=Path(__file__).resolve().parents[1]
rows=json.loads((ROOT/'evidence/independent_numerics.json').read_text())['periodic_shooting']
response=json.loads((ROOT/'evidence/period_response.json').read_text())
fig,axes=plt.subplots(1,3,figsize=(12,3.6),layout='constrained')
colors=['#ce5c37','#3c8b71','#4075bc','#222222']
for idx,(r,color) in enumerate(zip(rows,colors)):
    inf=not isinstance(r['K'],(float,int))
    K=np.inf if inf else r['K']
    if inf:
        y0=np.array(response['peak_state']);f=limit_rhs;jac=limit_jac;wi=2
    else:
        d=json.loads((ROOT/f'evidence/baseline_K{float(K)}_LSODA.json').read_text());y0=np.array(d['last_peak_state']);f=lambda t,y:finite_rhs(t,y,K);jac=lambda t,y:finite_jac(t,y,K);wi=3
    sol=solve_ivp(f,(0,10),y0,method='Radau',jac=jac,rtol=1e-10,atol=1e-12,dense_output=True)
    ts=np.linspace(0,10,1500)
    axes[0].plot(ts,sol.sol(ts)[wi],color=color,lw=1.3,ls='--' if inf else '-',label=r'$K=\infty$' if inf else f'$K={int(K)}$')
inv=np.array([0 if not isinstance(r['K'],(float,int)) else 1/r['K'] for r in rows])
period=np.array([r['period'] for r in rows]);amplitude=np.array([r['S2_peak_to_peak'] for r in rows])
for ax,vals,key,label in [(axes[1],period,'period_first_order_coefficient','Period'),(axes[2],amplitude,None,'Peak-to-peak amplitude')]:
    ax.scatter(inv,vals,c=colors,s=32,zorder=3,label='Independent shooting')
    xx=np.linspace(0,.0102,100)
    coeff=response[key] if key else response['maximum_first_order_coefficient']-response['minimum_first_order_coefficient']
    ax.plot(xx,vals[-1]+coeff*xx,color='#555555',ls='--',lw=1.2,label='First-order response (numerical)')
    ax.set(xlabel=r'$1/K$',ylabel=label)
    ax.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
    ax.legend(fontsize=7,frameon=False)
axes[0].set(xlabel='Time (model units)',ylabel=r'$S_2$',title='Finite concentration effects')
axes[0].legend(fontsize=7,frameon=False)
for a in axes:
    a.spines[['top','right']].set_visible(False)
    a.grid(alpha=.15)
axes[1].set_title('Convergence of period');axes[2].set_title('Convergence of amplitude')
out=ROOT/'evidence'
fig.savefig(out/'concentration_effects.pdf');fig.savefig(out/'concentration_effects.png',dpi=180)
print(out/'concentration_effects.png')
