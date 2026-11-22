# Copy of r11_fct_check.py (frozen in R11 manifest) extended with paired FCT - Mis_s for all 6 sources (r_orth_cal), for PREREG_E review.
# Written by Claude from the prereg text and the controller's own R10 check, not from A-T2/A-T3 code.
# Usage:
#   python3 -I r11_fct_check.py count   <r11_line>             # counts only (allowed before the code freeze)
#   python3 -I r11_fct_check.py holdout <r11_line> [n_boot]    # only after round11/status/CODE_FREEZE_R11.txt exists
#   python3 -I r11_fct_check.py dev     <r10_line> [n_boot]    # R10 lines (already unblinded) for parity checks
import sys, os, json, numpy as np
ROOT = '/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/'
B9 = ROOT + 'round9/p2/shared/'; B9L = B9 + 'lines/'
BASE = {'count': ROOT + 'round11/shared/lines/', 'holdout': ROOT + 'round11/shared/lines/', 'dev': ROOT + 'round10/shared/lines/'}
R9 = ['nadig_hepg2', 'nadig_jurkat', 'nourreddine_kolf21j', 'replogle_k562_essential', 'replogle_rpe1', 'vcc2025_h1']
EPS = 1e-8
LO, HI = 1.25, 98.75  # 97.5% two-sided (Bonferroni over 2 claims, PREREG_D v2 §4.1)

mode, tline = sys.argv[1], sys.argv[2]
n_boot = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
if mode == 'holdout':
    fz = ROOT + 'round11/status/CODE_FREEZE_R11.txt'
    if not (os.path.exists(fz) and os.path.getsize(fz) > 0):
        sys.exit('BLOCKED: CODE_FREEZE_R11.txt missing')

CG = [g.strip() for g in open(B9 + 'COMMON_GENES.txt') if g.strip()]
_appr, _prev, _alias = set(), {}, {}
for i, ln in enumerate(open('/Users/mac/.cache/r9_p2_raw/hgnc_complete_set.txt', encoding='utf-8')):
    f = ln.rstrip('\n').split('\t')
    if i == 0: H = {k: j for j, k in enumerate(f)}; continue
    if f[H['status']] != 'Approved': continue
    sym = f[H['symbol']]; _appr.add(sym)
    for col, dct in (('prev_symbol', _prev), ('alias_symbol', _alias)):
        for a in f[H[col]].strip('"').split('|'):
            if a: dct.setdefault(a, set()).add(sym)
def canon(x):
    if x in _appr: return x
    for dct in (_prev, _alias):
        if x in dct and len(dct[x]) == 1: return next(iter(dct[x]))
    return x

def load_meta(base, line):
    genes = [x.strip() for x in open(base + line + '/genes.txt')]
    d = np.load(base + line + '/pseudobulk_delta.npz')
    return genes, d, dict(t=[str(x) for x in d['targets']], rep=d['reproducible'], hv=d['half_valid'], n=d['n_cells'].astype(np.float64))

genes_t, dT, T = load_meta(BASE[mode], tline)
S = list(R9)
tix = {}
for i, g in enumerate(genes_t): tix.setdefault(canon(g), i)
GE, _seen = [], set()   # G_eval(t) = canon(R9 COMMON) ∩ genes(t), keeping the original R9 symbol
for g in CG:
    c = canon(g)
    if c in tix and c not in _seen: GE.append((c, g)); _seen.add(c)
gi_t = np.array([tix[c] for c, _ in GE]); gpos = {c: i for i, (c, _) in enumerate(GE)}; G = len(GE)

src = {}
for s in S:
    gs, dS, Ls = load_meta(B9L, s); ix = {g: i for i, g in enumerate(gs)}; gi = np.array([ix[g] for _, g in GE])
    mu_s = np.load(B9L + s + '/ctrl_cells.npz')['ctrl_mean_libnorm'].astype(np.float64)
    src[s] = dict(mu=np.maximum(mu_s[gi], EPS), D=dS['delta_log1p'][:, gi].astype(np.float64),
                  rows={canon(p): i for i, p in enumerate(Ls['t']) if Ls['rep'][i]})

rows = []
for ti, p in enumerate(T['t']):
    p = canon(p)
    if not T['rep'][ti] or p not in gpos: continue
    Sp = [s for s in S if p in src[s]['rows']]
    if Sp: rows.append((ti, p, Sp))
N = len(rows); NHV = int(sum(bool(T['hv'][ti]) for ti, _, _ in rows))
print(f'{tline}: eligible targets N={N}, half_valid={NHV}, |G_eval|={G}, sources={len(S)}', flush=True)
if mode == 'count':
    print(json.dumps(dict(line=tline, N=N, half_valid=NHV, G_eval=G, n_reproducible_on_t=int(np.sum(T['rep'])))))
    sys.exit(0)

tid = np.array([ti for ti, _, _ in rows])
Dfull = dT['delta_log1p'][tid][:, gi_t].astype(np.float64)
Afull = dT['half_a_delta_log1p'][tid][:, gi_t].astype(np.float64)
Bfull = dT['half_b_delta_log1p'][tid][:, gi_t].astype(np.float64)
HVr = np.array([bool(T['hv'][ti]) for ti in tid])
mu_t = np.maximum(np.load(BASE[mode] + tline + '/ctrl_cells.npz')['ctrl_mean_libnorm'].astype(np.float64)[gi_t], EPS)
TMf = np.zeros((N, G)); LFCf = np.zeros((N, G)); selfg = np.zeros(N, int)
for r, (ti, p, Sp) in enumerate(rows):
    TMf[r] = np.mean([src[s]['D'][src[s]['rows'][p]] for s in Sp], axis=0)
    LFCf[r] = np.mean([np.log(np.maximum(np.expm1(np.log1p(src[s]['mu']) + src[s]['D'][src[s]['rows'][p]]), EPS) / src[s]['mu']) for s in Sp], axis=0)
    selfg[r] = gpos[p]

def masks(kind):
    out = []
    for r in range(N):
        if kind == 'all_nontarget':
            m = np.ones(G, bool); m[selfg[r]] = False; m &= np.isfinite(Dfull[r]) & np.isfinite(TMf[r])
            out.append(np.flatnonzero(m))
        else:  # resp_top50: top 50 by the target line's observed |Δ|, self gene excluded
            dd = np.abs(Dfull[r]).copy(); dd[selfg[r]] = -1; dd[~np.isfinite(dd)] = -1
            out.append(np.argsort(-dd, kind='stable')[:50])
    return out

def orth_row(x, u): return x - np.dot(x, u) * u

class Mask:
    def __init__(self, kind):
        self.kind = kind; self.M = masks(kind)
        self.U, self.O = [], []
        self.cs = np.zeros((N, 6))   # ceiling sums per target: n, Σa, Σb, Σab, Σa², Σb²
        for r, m in enumerate(self.M):
            tm = TMf[r, m]; u = tm / max(np.linalg.norm(tm), 1e-12); self.U.append(u)
            self.O.append(orth_row(Dfull[r, m] - tm, u))
            if HVr[r]:
                a = Afull[r, m]; b = Bfull[r, m]; ok = np.isfinite(a) & np.isfinite(b)
                if ok.all():
                    a = orth_row(a - tm, u); b = orth_row(b - tm, u)
                else:  # pairwise-complete fallback (reported)
                    a = orth_row(np.where(ok, a, tm) - tm, u)[ok]; b = orth_row(np.where(ok, b, tm) - tm, u)[ok]
                self.cs[r] = [a.size, a.sum(), b.sum(), (a * b).sum(), (a * a).sum(), (b * b).sum()]
    def stats(self, mu_vec=None, P=None):
        st = np.zeros((N, 6))  # n, Σx, Σy, Σxy, Σx², Σy²  (x = observed residual, y = prediction residual)
        for r, m in enumerate(self.M):
            tm = TMf[r, m]
            if P is None:
                mu = mu_vec[m]; pr = np.log1p(mu * np.exp(LFCf[r, m])) - np.log1p(mu)
            else:
                pr = P[r, m]
            y = orth_row(pr - tm, self.U[r]); x = self.O[r]
            st[r] = [x.size, x.sum(), y.sum(), (x * y).sum(), (x * x).sum(), (y * y).sum()]
        return st

def pooled_r(st, w):
    n, sx, sy, sxy, sxx, syy = (w @ st)
    cov = sxy - sx * sy / n; vx = sxx - sx * sx / n; vy = syy - sy * sy / n
    return cov / np.sqrt(vx * vy) if vx > 0 and vy > 0 else np.nan

def cal(st, cs, w):
    rab = pooled_r(cs, w * HVr)
    if not rab > 0: return np.nan, rab
    return float(np.clip(pooled_r(st, w) / np.sqrt(2 * rab / (1 + rab)), -1, 1)), rab

preds_mu = {'FCT': mu_t, 'E1': np.full(G, np.median(mu_t))}
for s in S: preds_mu['Mis_' + s] = src[s]['mu']
res = {}
rng = np.random.default_rng(777)
W = [np.bincount(rng.integers(0, N, N), minlength=N).astype(float) for _ in range(n_boot)]
ones = np.ones(N)
for kind in ('all_nontarget', 'resp_top50'):
    mk = Mask(kind)
    ST = {k: mk.stats(mu_vec=v) for k, v in preds_mu.items()}
    pt = {k: cal(v, mk.cs, ones)[0] for k, v in ST.items()}
    rab_pt = cal(ST['FCT'], mk.cs, ones)[1]
    perm = [cal(mk.stats(mu_vec=mu_t[np.random.default_rng(k).permutation(G)]), mk.cs, ones)[0] for k in range(100)]
    BS = {k: np.array([cal(v, mk.cs, w)[0] for w in W]) for k, v in ST.items()}
    ci = {k: tuple(np.nanpercentile(b, [LO, HI])) for k, b in BS.items()}
    valid = {k: float(np.mean(np.isfinite(b))) for k, b in BS.items()}
    d3 = BS['FCT'] - BS['E1']
    mis_hi = max(ci[k][1] for k in ci if k.startswith('Mis_'))
    res[kind] = dict(n_genes_median=float(np.median([len(m) for m in mk.M])), rho_ab_orth=rab_pt,
                     point=pt, ci=ci, valid_frac=valid, perm_p975=float(np.nanpercentile(perm, 97.5)),
                     perm_median=float(np.nanmedian(perm)), max_mis_ci_hi=mis_hi,
                     FCT_minus_E1=(pt['FCT'] - pt['E1'], tuple(np.nanpercentile(d3, [LO, HI]))),
                     FCT_minus_Mis_kolf=(pt['FCT'] - pt['Mis_nourreddine_kolf21j'],
                                         tuple(np.nanpercentile(BS['FCT'] - BS['Mis_nourreddine_kolf21j'], [LO, HI]))),
                     paired={k[4:]: (pt['FCT'] - pt[k], tuple(np.nanpercentile(BS['FCT'] - BS[k], [LO, HI]))) for k in BS if k.startswith('Mis_')},
                     _BS_E1=BS['E1'])
a, t5 = res['all_nontarget'], res['resp_top50']
d4b = t5['_BS_E1'] - a['_BS_E1']
for k in res: res[k].pop('_BS_E1')
out = dict(line=tline, mode=mode, N=N, half_valid=NHV, G_eval=G, n_boot=n_boot, seed=777, results=res,
           H_D1=bool(a['ci']['FCT'][0] > 0 and a['valid_frac']['FCT'] >= 0.95 and a['rho_ab_orth'] > 0),
           H_D2=bool(a['ci']['FCT'][0] > max(a['max_mis_ci_hi'], a['perm_p975'])),
           H_D3=bool(a['FCT_minus_E1'][1][0] > 0),
           H_D4a=bool(t5['ci']['E1'][0] > 0),
           H_D4b=bool(np.nanpercentile(d4b, LO) > 0), E1_top50_minus_all=(t5['point']['E1'] - a['point']['E1'], tuple(np.nanpercentile(d4b, [LO, HI]))))
print(json.dumps(out, indent=1, default=float))
