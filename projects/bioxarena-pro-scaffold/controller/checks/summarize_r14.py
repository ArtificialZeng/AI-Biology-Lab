import json, sys
t = open('controller/checks/r14_holdout_a549_myllia.out').read()
m = json.loads(t[t.index('{'):])
a = m['results']['all_nontarget']
r = lambda x: round(x, 4)
print('N', m['N'], 'G_eval', m['G_eval'], 'rho_orth', r(a['rho_ab_orth']))
print('FCT', r(a['point']['FCT']), [r(x) for x in a['ci']['FCT']], 'valid', a['valid_frac']['FCT'], 'perm975', r(a['perm_p975']))
for s, (pt, (lo, hi)) in a['paired'].items():
    print('  paired', s, r(pt), [r(lo), r(hi)])
print('max Mis hi', r(a['max_mis_ci_hi']), 'E1 all', r(a['point']['E1']))
print('E1 split', r(m['split']['E1']['point']), [r(x) for x in m['split']['E1']['ci']], 'H_E2b', r(m['H_E2b'][0]), [r(x) for x in m['H_E2b'][1]])
print('flags', {k: m[k] for k in ['H_E1a', 'H_E1b', 'H_E1c', 'H_E2a', 'H_E2b_pass', 'H_D2_unpaired_descriptive']})
print('SGM / CIPHER on all_nontarget not computed by this checker')
man = dict((l.split()[1], l.split()[0]) for l in open('round14/status/FREEZE_MANIFEST_R14.txt') if not l.startswith('#'))
import hashlib
h = hashlib.sha256(open('controller/checks/r14_check.py', 'rb').read()).hexdigest()
print('r14_check.py matches manifest:', man.get('controller/checks/r14_check.py') == h)
