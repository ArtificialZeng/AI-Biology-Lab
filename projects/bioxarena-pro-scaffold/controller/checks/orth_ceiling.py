# Master-controller dev-only: split-half ceiling and bootstrap CI for amplitude-orthogonalized residual correlation (TMVS, FCT).
import numpy as np
exec(open(__import__('os').path.join(__import__('os').path.dirname(__file__),'tmvs_check.py')).read().split("rng=np.random.default_rng(0)")[0])
rng=np.random.default_rng(2)
def orth(r,b): b=b/max(np.linalg.norm(b),1e-12); return r-(r@b)*b
for tgt,src in [('replogle_k562_essential','replogle_rpe1'),('replogle_rpe1','replogle_k562_essential')]:
    T,S=L[tgt],L[src]
    st={p:i for i,p in enumerate(S['t']) if S['rep'][i]}
    rows=[(i,st[p]) for i,p in enumerate(T['t']) if T['rep'][i] and p in st and p in cidx]
    ratio=np.sqrt(np.maximum(T['var'],1e-8)/np.maximum(S['var'],1e-8))
    O=[];A=[];Bb=[];Pt=[];Pf=[];Om=[];Am=[];Bm=[];Ptm=[]
    for ti,si in rows:
        g=cidx[T['t'][ti]]; d=T['D'][ti].copy(); d[g]=0
        m=np.argsort(-np.abs(d)); m=m[m!=g][:50]
        ds=S['D'][si][m]
        O.append(orth(T['D'][ti][m]-ds,ds)); A.append(orth(T['A'][ti][m]-ds,ds)); Bb.append(orth(T['Bh'][ti][m]-ds,ds))
        Pt.append(orth(ds*ratio[m]-ds,ds))
        # all-nontarget mask version
        ma=np.arange(len(common)); ma=ma[ma!=g]; dsa=S['D'][si][ma]
        Om.append(orth(T['D'][ti][ma]-dsa,dsa)); Am.append(orth(T['A'][ti][ma]-dsa,dsa)); Bm.append(orth(T['Bh'][ti][ma]-dsa,dsa)); Ptm.append(orth(dsa*ratio[ma]-dsa,dsa))
    O,A,Bb,Pt=map(np.array,(O,A,Bb,Pt)); Om,Am,Bm,Ptm=map(np.array,(Om,Am,Bm,Ptm))
    def st_(idx,O,A,B,P):
        rab=np.corrcoef(A[idx].ravel(),B[idx].ravel())[0,1]
        c=np.sqrt(2*rab/(1+rab)) if rab>0 else np.nan
        return np.corrcoef(O[idx].ravel(),P[idx].ravel())[0,1]/c, rab
    for lab,arrs in [('top50',(O,A,Bb,Pt)),('all_nontarget',(Om,Am,Bm,Ptm))]:
        full,rab=st_(np.arange(len(rows)),*arrs)
        bs=[st_(rng.integers(0,len(rows),len(rows)),*arrs)[0] for _ in range(300)]
        lo,hi=np.nanpercentile(bs,[1.25,98.75])
        print(f'{tgt}<-{src} [{lab}] rho_AB_orth={rab:.3f}  TMVS r_orth_cal={full:+.4f} [97.5% CI {lo:+.4f},{hi:+.4f}]')
