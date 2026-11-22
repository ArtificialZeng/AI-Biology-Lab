# Master-controller dev-only check: is TMVS signal gene-specific or a global amplitude (regression-to-mean) effect?
import numpy as np
exec(open(__import__('os').path.join(__import__('os').path.dirname(__file__),'tmvs_check.py')).read().split("rng=np.random.default_rng(0)")[0])
rng=np.random.default_rng(1)
def orth(r,base):  # remove per-perturbation component along base (transfer mean) on the masked genes
    b=base/ max(np.linalg.norm(base),1e-12); return r-(r@b)*b
for tgt,src in [('replogle_k562_essential','replogle_rpe1'),('replogle_rpe1','replogle_k562_essential')]:
    T,S=L[tgt],L[src]
    st={p:i for i,p in enumerate(S['t']) if S['rep'][i]}
    rows=[(i,st[p]) for i,p in enumerate(T['t']) if T['rep'][i] and p in st and p in cidx]
    ratio=np.sqrt(np.maximum(T['var'],1e-8)/np.maximum(S['var'],1e-8))
    cglob=np.sqrt(T['var'].mean()/S['var'].mean())          # zero-shot gene-agnostic amplitude
    perm=rng.permutation(len(ratio))
    P={'TMVS':[], 'GlobalScale':[], 'ShuffledVar':[]}; O=[]; PO={k:[] for k in P}; OO=[]
    for ti,si in rows:
        g=cidx[T['t'][ti]]; d=T['D'][ti].copy(); d[g]=0
        m=np.argsort(-np.abs(d)); m=m[m!=g][:50]
        ds=S['D'][si][m]; ro=T['D'][ti][m]-ds
        preds={'TMVS':ds*ratio[m]-ds,'GlobalScale':ds*cglob-ds,'ShuffledVar':ds*ratio[perm][m]-ds}
        O.append(ro); OO.append(orth(ro,ds))
        for k,v in preds.items(): P[k].append(v); PO[k].append(orth(v,ds))
    O=np.concatenate(O); OO=np.concatenate(OO)
    beta=np.polyfit(np.concatenate([S['D'][si][np.argsort(-np.abs(T['D'][ti]))[:51]] for ti,si in rows]), np.concatenate([T['D'][ti][np.argsort(-np.abs(T['D'][ti]))[:51]] for ti,si in rows]),1)[0]
    print(f'{tgt}<-{src}: global var ratio c={cglob:.3f}; oracle amplitude beta(Dt~Ds, top genes)={beta:.3f}')
    for k in P:
        p=np.concatenate(P[k]); po=np.concatenate(PO[k])
        r=np.corrcoef(O,p)[0,1] if p.std()>0 else np.nan
        ro=np.corrcoef(OO,po)[0,1] if po.std()>1e-12 else np.nan
        print(f'  {k:12s} r_resid={r:+.4f}   r_orth(amplitude removed)={ro:+.4f}')
