import json, glob
S = json.load(open(sorted(glob.glob('round14/results/summary_holdout_r14*.json'))[0]))
def find(o):
    if isinstance(o, dict):
        for k, v in o.items():
            if k == 'overlap_subgroups': return v
            r = find(v)
            if r is not None: return r
sg = find(S)['a549_myllia']
for name, g in sg.items():
    if not isinstance(g, dict) or 'fct' not in g: continue
    f = g['fct']
    print(f"{name}: n={g['n_targets']} FCT {f['estimate']:+.3f} [{f['ci975_lo']:+.3f}, {f['ci975_hi']:+.3f}]")
    for s, v in g['paired_fct_minus_mismatched'].items():
        print(f"   {s}: {v['estimate']:+.3f} [{v['ci975_lo']:+.3f}, {v['ci975_hi']:+.3f}]")
