import json, glob
adj = sorted(glob.glob('round14/results/adjudication_holdout_r14*.json'))
print('A-T2 files:', adj)
a = json.load(open(adj[0]))
print('k_conf', a.get('k_conf'), a.get('confirmatory_status'))
D = a.get('all_lines_descriptive', a)
t = D['H-E1']['per_line']['a549_myllia']; t2 = D['H-E2']['per_line']['a549_myllia']
s = open('controller/checks/r14_holdout_a549_myllia.out').read(); m = json.loads(s[s.index('{'):]); an = m['results']['all_nontarget']
pairs = [(t['fct_r_orth_cal'], an['point']['FCT']), (t['fct_ci975_lo'], an['ci']['FCT'][0]), (t['fct_ci975_hi'], an['ci']['FCT'][1]),
         (t['perm100_all_nontarget_q975'], an['perm_p975']), (t['descriptive_unpaired_max_mis_ci975_hi'], an['max_mis_ci_hi']),
         (t2['e1_split_r_orth_cal'], m['split']['E1']['point']), (t2['e1_split_ci975_lo'], m['split']['E1']['ci'][0]), (t2['e1_split_ci975_hi'], m['split']['E1']['ci'][1]),
         (t2['e1_split_minus_all_est'], m['H_E2b'][0]), (t2['e1_split_minus_all_ci975_lo'], m['H_E2b'][1][0]), (t2['e1_split_minus_all_ci975_hi'], m['H_E2b'][1][1])]
for sname, v in t['H-E1b_per_source'].items():
    mine = an['paired'][sname]
    est = v.get('est', v.get('paired_diff', v.get('estimate', v.get('diff_est'))))
    lo = v.get('ci975_lo', v.get('diff_ci975_lo')); hi = v.get('ci975_hi', v.get('diff_ci975_hi'))
    pairs += [(est, mine[0]), (lo, mine[1][0]), (hi, mine[1][1])]
flags = [t['H-E1a_pass'] == m['H_E1a'], t['H-E1b_pass'] == m['H_E1b'], t['H-E1c_pass'] == m['H_E1c'], t2['H-E2a_pass'] == m['H_E2a'], t2['H-E2b_pass'] == m['H_E2b_pass']]
print('A-T2 vs controller: %d numbers, max|diff| %.2e, flags agree %s (%d)' % (len(pairs), max(abs(x - y) for x, y in pairs), all(flags), len(flags)))
summ = sorted(glob.glob('round14/results/summary_holdout_r14*.json'))
S = json.load(open(summ[0]))
def find(o, key, path=''):
    if isinstance(o, dict):
        for k, v in o.items():
            if key in str(k).lower(): yield path + '/' + str(k), v
            yield from find(v, key, path + '/' + str(k))
for p, v in find(S, 'subgroup'):
    print(p); print(json.dumps(v, indent=1)[:2500]); break
