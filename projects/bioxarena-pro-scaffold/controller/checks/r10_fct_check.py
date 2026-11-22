# Master-controller independent check for PREREG_C (dc4b4a5): FCT r_orth_cal, 6 mismatched control means,
# 100 gene-label permutations, joint bootstrap. Written by Claude from the prereg text, not from A-T2/A-T3 code.
# Usage:
#   python3 -I r10_fct_check.py dev <r9_target_line>            # R9 leave-one-line-out dev test (allowed any time)
#   python3 -I r10_fct_check.py holdout <r10_line> [n_boot]     # only after round10/status/CODE_FREEZE_R10.txt exists
import sys, os, json, numpy as np
ROOT = '/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/'
B9 = ROOT + 'round9/p2/shared/'; B9L = B9 + 'lines/'
B10 = ROOT + 'round10/shared/lines/'
R9 = ['nadig_hepg2', 'nadig_jurkat', 'nourreddine_kolf21j', 'replogle_k562_essential', 'replogle_rpe1', 'vcc2025_h1']
EPS = 1e-8

mode, tline = sys.argv[1], sys.argv[2]
n_boot = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
if mode == 'holdout' and not os.path.exists(ROOT + 'round10/status/CODE_FREEZE_R10.txt'):
    sys.exit('BLOCKED: CODE_FREEZE_R10.txt missing')

CG = [g.strip() for g in open(B9 + 'COMMON_GENES.txt') if g.strip()]
# HGNC canonicalisation (approved symbol; else unambiguous previous symbol; else unambiguous alias)
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

def load(base, line):
    genes = [x.strip() for x in open(base + line + '/genes.txt')]
    c = np.load(base + line + '/ctrl_cells.npz'); d = np.load(base + line + '/pseudobulk_delta.npz')
    return dict(genes=genes, mu=c['ctrl_mean_libnorm'].astype(np.float64), D=d['delta_log1p'], A=d['half_a_delta_log1p'],
                Bh=d['half_b_delta_log1p'], t=[str(x) for x in d['targets']], rep=d['reproducible'],
                hv=d['half_valid'], n=d['n_cells'].astype(np.float64))

T = load(B10 if mode == 'holdout' else B9L, tline)
S = [s for s in R9 if s != tline]
tix = {}
for i, g in enumerate(T['genes']): tix.setdefault(canon(g), i)
GE, _seen = [], set()                                   # G_eval(t) = COMMON_GENES ∩ target genes (canonical, original R9 symbol)
for g in CG:
    c = canon(g)
    if c in tix and c not in _seen: GE.append((c, g)); _seen.add(c)
gi_t = np.array([tix[c] for c, _ in GE]); gpos = {c: i for i, (c, _) in enumerate(GE)}
src = {}
for s in S:
    L = load(B9L, s); ix = {g: i for i, g in enumerate(L['genes'])}; gi = np.array([ix[g] for _, g in GE])
    src[s] = dict(mu=np.maximum(L['mu'][gi], EPS), D=L['D'][:, gi].astype(np.float64),
                  rows={canon(p): i for i, p in enumerate(L['t']) if L['rep'][i]})
mu_t = np.maximum(T['mu'][gi_t], EPS)

rows = []
for ti, p in enumerate(T['t']):
    p = canon(p)
    if not T['rep'][ti] or not T['hv'][ti] or p not in gpos: continue
    Sp = [s for s in S if p in src[s]['rows']]
    if not Sp: continue
    rows.append((ti, p, Sp))
N = len(rows); print(f'{tline}: eligible targets N={N}, |G_eval|={len(GE)}, sources={len(S)}', flush=True)

# per-target 50-gene masks and source quantities
M, TM, LFC, Dobs, Da, Db, ncell = [], [], [], [], [], [], []
for ti, p, Sp in rows:
    d = T['D'][ti][gi_t].astype(np.float64); g = gpos[p]
    dd = np.abs(d); dd[g] = -1; m = np.argsort(-dd, kind='stable')[:50]
    tm = np.mean([src[s]['D'][src[s]['rows'][p]] for s in Sp], axis=0)
    lfc = np.mean([np.log(np.maximum(np.expm1(np.log1p(src[s]['mu']) + src[s]['D'][src[s]['rows'][p]]), EPS) / src[s]['mu']) for s in Sp], axis=0)
    M.append(m); TM.append(tm[m]); LFC.append(lfc[m]); Dobs.append(d[m])
    Da.append(T['A'][ti][gi_t][m].astype(np.float64)); Db.append(T['Bh'][ti][gi_t][m].astype(np.float64)); ncell.append(T['n'][ti])
M, TM, LFC, Dobs, Da, Db = map(np.array, (M, TM, LFC, Dobs, Da, Db))

def orth(R):  # project each row off its Transfer_Mean vector
    u = TM / np.maximum(np.linalg.norm(TM, axis=1, keepdims=True), 1e-12)
    return R - np.sum(R * u, axis=1, keepdims=True) * u
def fct(mu_full):
    mu = mu_full[M]; return np.log1p(mu * np.exp(LFC)) - np.log1p(mu)

Oo, Ao, Bo = orth(Dobs - TM), orth(Da - TM), orth(Db - TM)
def score(P, idx):
    rab = np.corrcoef(Ao[idx].ravel(), Bo[idx].ravel())[0, 1]
    if not rab > 0: return np.nan
    return float(np.clip(np.corrcoef(Oo[idx].ravel(), P[idx].ravel())[0, 1] / np.sqrt(2 * rab / (1 + rab)), -1, 1))

preds = {'FCT': orth(fct(mu_t) - TM)}
for s in S: preds['Mis_' + s] = orth(fct(src[s]['mu']) - TM)
allidx = np.arange(N)
point = {k: score(v, allidx) for k, v in preds.items()}
perm = [score(orth(fct(mu_t[np.random.default_rng(k).permutation(len(GE))]) - TM), allidx) for k in range(100)]
rng = np.random.default_rng(777)  # controller's own seed (agents use 20261122)
B = {k: [] for k in preds}
for _ in range(n_boot):
    idx = rng.integers(0, N, N)
    for k, v in preds.items(): B[k].append(score(v, idx))
ci = {k: tuple(np.nanpercentile(v, [1.25, 98.75])) for k, v in B.items()}
mis_hi = max(ci[k][1] for k in preds if k.startswith('Mis_')); p975 = float(np.percentile(perm, 97.5))
out = dict(line=tline, mode=mode, N=N, FCT=point['FCT'], FCT_ci=ci['FCT'], max_mis_ci_hi=mis_hi, perm_p975=p975,
           mis={k: (point[k], ci[k]) for k in preds if k.startswith('Mis_')},
           H_C1=bool(ci['FCT'][0] > 0), H_C2=bool(ci['FCT'][0] > max(mis_hi, p975)))
print(json.dumps(out, indent=1, default=float))
