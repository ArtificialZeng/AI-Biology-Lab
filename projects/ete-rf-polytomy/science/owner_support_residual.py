from pathlib import Path
import sys,json
R=Path(__file__).resolve().parent;runtime=Path(sys.argv[1]).resolve()
sys.path.insert(0,str(runtime));sys.path.insert(1,str((R.parent/'compbio-primary-triage/evidence/runtime').resolve()))
from ete4 import Tree
cases=[]
for a,b in [('((A,B),(C,D));','((A,C),(B,D));'),('((A,B,C),D);','(((A,B),C),D);'),('((A,B,C),D);','((A,C,B),D);')]:
 for unrooted in [False,True]:
  for correction in [False,True]:
   for which in [1,2]:
    kw={'unrooted_trees':unrooted,'correct_by_polytomy_size':correction,'min_support_t'+str(which):0.5}
    try:res=Tree(a).robinson_foulds(Tree(b),**kw);result={'rf':res[0],'max_rf':res[1]}
    except Exception as e:result={'error_type':type(e).__name__,'error':str(e)}
    cases.append({'a':a,'b':b,'kwargs':kw,**result})
data={'runtime':str(runtime),'cases':cases,'claim':'Observed inherited support-filter failures, not passes or input prohibitions.'}
(R/'evidence'/('support-residual-'+runtime.name+'.json')).write_text(json.dumps(data,indent=2)+'\n');print(json.dumps({'runtime':runtime.name,'cases':len(cases),'error_types':sorted(set(c.get('error_type','none') for c in cases))}))
