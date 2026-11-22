# Master-controller independent held-out check (post-lock, PREREG v4 4cc89eb, CODE_FREEZE 1ca9f60): TMVS r_orth_cal per confirmatory line.
import numpy as np
B='/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2/shared/'
CG=[g.strip() for g in open(B+'COMMON_GENES.txt') if g.strip()]; G=len(CG)
LINES=['nadig_hepg2','nadig_jurkat','nourreddine_kolf21j','replogle_k562_essential','replogle_rpe1','vcc2025_h1']
L={}
for l in LINES:
    genes=[x.strip() for x in open(B+'lines/'+l+'/genes.txt')]; ix={g:i for i,g in enumerate(genes)}; gi=np.array([ix[g] for g in CG])
    c=np.load(B+'lines/'+l+'/ctrl_cells.npz'); X=c['X_log1p'][:,gi].astype(np.float64); var=X.var(axis=0); del X
    d=np.load(B+'lines/'+l+'/pseudobulk_delta.npz')
    L[l]=dict(var=var,D=d['delta_log1p'][:,gi],A=d['half_a_delta_log1p'][:,gi],Bh=d['half_b_delta_log1p'][:,gi],t=list(d['targets']),rep=d['reproducible'])
    print('loaded',l,flush=True)
cidx={g:i for i,g in enumerate(CG)}
def orth(r,b): b=b/max(np.linalg.norm(b),1e-12); return r-(r@b)*b
rngb=np.random.default_rng(20261121)
for t in ['nadig_hepg2','nadig_jurkat','vcc2025_h1']:
    T=L[t]; S=[s for s in LINES if s!=t]
    srcmap={s:{p:i for i,p in enumerate(L[s]['t']) if L[s]['rep'][i]} for s in S}
    Cbar=np.mean([L[s]['var'] for s in S],axis=0)
    ratio=np.sqrt(np.maximum(T['var'],1e-8)/np.maximum(Cbar,1e-8))
    O=[];A=[];Bb=[];P=[]; rows=[]
    for ti,p in enumerate(T['t']):
        if not T['rep'][ti] or p not in cidx: continue
        Sp=[s for s in S if p in srcmap[s]]
        if not Sp: continue
        dsrc=np.mean([L[s]['D'][srcmap[s][p]] for s in Sp],axis=0)
        g=cidx[p]; d=T['D'][ti].copy(); d[g]=0
        m=np.argsort(-np.abs(d)); m=m[m!=g][:50]; ds=dsrc[m]
        O.append(orth(T['D'][ti][m]-ds,ds)); A.append(orth(T['A'][ti][m]-ds,ds)); Bb.append(orth(T['Bh'][ti][m]-ds,ds)); P.append(orth(ds*ratio[m]-ds,ds))
    O,A,Bb,P=map(np.array,(O,A,Bb,P)); n=len(O)
    def st(idx):
        rab=np.corrcoef(A[idx].ravel(),Bb[idx].ravel())[0,1]; c=np.sqrt(2*rab/(1+rab)) if rab>0 else np.nan
        return np.corrcoef(O[idx].ravel(),P[idx].ravel())[0,1]/c, rab
    full,rab=st(np.arange(n)); bs=[st(rngb.integers(0,n,n))[0] for _ in range(1000)]
    lo,hi=np.nanpercentile(bs,[1.25,98.75])
    print(f'{t}: n_targets={n} rho_AB_orth={rab:.4f} TMVS r_orth_cal={full:+.4f} [97.5% CI {lo:+.4f},{hi:+.4f}]',flush=True)
