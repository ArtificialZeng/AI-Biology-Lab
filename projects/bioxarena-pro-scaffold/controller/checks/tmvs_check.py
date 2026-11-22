# Master-controller independent dev-only check (Replogle K562<->RPE1): Transfer_Mean_VarScaled vs Transfer_Mean residual signal.
import numpy as np, sys
B='/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2/shared/lines/'
def load(l):
    c=np.load(B+l+'/ctrl_cells.npz'); d=np.load(B+l+'/pseudobulk_delta.npz')
    g=[x.strip() for x in open(B+l+'/genes.txt')]
    X=c['X_log1p']; var=X.var(axis=0).astype(np.float64)
    return dict(g=g,var=var,mu=c['ctrl_mean_libnorm'].astype(np.float64),D=d['delta_log1p'],A=d['half_a_delta_log1p'],Bh=d['half_b_delta_log1p'],t=list(d['targets']),rep=d['reproducible'])
L={l:load(l) for l in ['replogle_k562_essential','replogle_rpe1']}
common=sorted(set(L['replogle_k562_essential']['g'])&set(L['replogle_rpe1']['g'])); print('common genes',len(common))
for l in L:
    ix={g:i for i,g in enumerate(L[l]['g'])}; gi=np.array([ix[g] for g in common])
    for k in ['var','mu']: L[l][k]=L[l][k][gi]
    for k in ['D','A','Bh']: L[l][k]=L[l][k][:,gi]
cidx={g:i for i,g in enumerate(common)}
rng=np.random.default_rng(0)
for tgt,src in [('replogle_k562_essential','replogle_rpe1'),('replogle_rpe1','replogle_k562_essential')]:
    T,S=L[tgt],L[src]
    st={p:i for i,p in enumerate(S['t']) if S['rep'][i]}
    rows=[(i,st[p]) for i,p in enumerate(T['t']) if T['rep'][i] and p in st and p in cidx]
    ratio=np.sqrt(np.maximum(T['var'],1e-8)/np.maximum(S['var'],1e-8))
    # fold-change transfer through target libnorm means
    out={'TMVS':[], 'FCT':[]}; obs=[]; ha=[]; hb=[]; per=[]
    for ti,si in rows:
        g=cidx[T['t'][ti]]; d=T['D'][ti].copy(); d[g]=0
        m=np.argsort(-np.abs(d)); m=m[m!=g][:50]
        ds=S['D'][si][m]; ro=T['D'][ti][m]-ds
        tm=ds*ratio[m]-ds
        # FCT: approximate source fold change in mean space via log1p deltas, re-applied on target mean
        mus=S['mu'][m]; mut=T['mu'][m]
        fc=np.maximum(np.expm1(np.log1p(mus)+ds),0)/np.maximum(mus,1e-6)
        fct=(np.log1p(mut*fc)-np.log1p(mut))-ds
        obs.append(ro); out['TMVS'].append(tm); out['FCT'].append(fct)
        ha.append(T['A'][ti][m]-ds); hb.append(T['Bh'][ti][m]-ds); per.append(len(m))
    obs=np.array(obs); ha=np.array(ha); hb=np.array(hb)
    def stat(idx):
        o=obs[idx].ravel(); rab=np.corrcoef(ha[idx].ravel(),hb[idx].ravel())[0,1]
        ceil=np.sqrt(2*rab/(1+rab)) if rab>0 else np.nan
        res={}
        for k,v in out.items():
            p=np.array(v)[idx].ravel(); r=np.corrcoef(o,p)[0,1]
            slope=np.cov(p,o)[0,1]/o.var()
            res[k]=(r,r/ceil,slope)
        return res,rab
    full,rab=stat(np.arange(len(rows)))
    boots={k:[] for k in out}
    for b in range(400):
        idx=rng.integers(0,len(rows),len(rows)); res,_=stat(idx)
        for k in out: boots[k].append(res[k][1])
    print(f'{tgt} <- {src}: n_targets={len(rows)} rho_AB={rab:.3f}')
    for k in out:
        lo,hi=np.nanpercentile(boots[k],[2.5,97.5])
        print(f'  {k}: r={full[k][0]:+.4f} r_cal={full[k][1]:+.4f} [{lo:+.4f},{hi:+.4f}] slope={full[k][2]:+.4f}')
