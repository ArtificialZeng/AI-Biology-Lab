"""Fresh finite event-count oracle; no owner helper and no candidate timing helper."""
import ast, fractions, hashlib, importlib, json, math, os, pathlib, random, subprocess, sys, urllib.request
P=pathlib.Path(__file__).resolve().parent
R=P.parent/'runtime'
F=fractions.Fraction
rng=random.Random(1402026)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def integrity():
 s=json.loads((R/'snapshot.json').read_text());bad=[];actual_files={}
 for name,h in s['files'].items():
  q=R/name
  actual_files[name]=sha(q) if q.is_file() else None
  if actual_files[name]!=h:bad.append(name)
 digest=hashlib.sha256(json.dumps(s['files'],sort_keys=True,separators=(',',':')).encode()).hexdigest()
 if digest!=s['snapshot_id']:bad.append('canonical digest')
 return {'snapshot_id':digest,'files_checked':len(s['files']),'failures':bad,'actual_files':actual_files}
def inventory(root):return {str(q.relative_to(root)):sha(q) for q in root.rglob('*') if q.is_file() and '.git' not in q.parts and '__pycache__' not in q.parts and q.suffix!='.pyc'}
def load(root):
 for key in list(sys.modules):
  if key=='treeswift' or key.startswith('treeswift.'):del sys.modules[key]
 sys.path.insert(0,str(root));module=importlib.import_module('treeswift');sys.path.pop(0);return module
# Each event merges d extant lineages at age a (backwards from present).
# The oracle integrates the lineage count forward and sums the integrated
# lineage time at each of the n-2 noncrown branching events, counting a d-way
# event d-1 times. It never reads the serialized tree or resolves polytomies.
def expectation(n,events):
 H=max(a for a,d in events);jumps={}
 for a,d in events:jumps[H-a]=jumps.get(H-a,0)+d-1
 count=1;integral=F(0);weighted=F(0);previous=F(0);first=True
 for time,jump in sorted(jumps.items()):
  integral+=count*(time-previous)
  weighted+=(jump-1 if first else jump)*integral
  count+=jump;previous=time;first=False
 integral+=count*(H-previous)
 if count!=n:raise AssertionError('bad lineage endpoint')
 if not integral:return None
 coefficient=(weighted/F(n-2)-integral/2)/integral
 return float(coefficient)*math.sqrt(12*(n-2))
def serialize(node,parent_age=None,labelmode=0,scale=F(1),subdivide=False):
 a,kids,label=node
 text=('('+','.join(serialize(c,a,labelmode,scale,subdivide) for c in kids)+')') if kids else ''
 text+= label if labelmode==0 else ('X' if labelmode==1 else '')
 if parent_age is not None:
  length=(parent_age-a)*scale
  if subdivide and length:
   text='('+text+':'+str(float(length/2))+')U:'+str(float(length/2))
  else:text+=':'+str(float(length))
 return text
cases=[]
for n in range(3,9):
 for z in range(60):
  active=[(F(0),[],chr(65+i)) for i in range(n)];events=[];age=F(0)
  while len(active)>1:
   age+=rng.choice([F(0),F(1,8),F(1,4),F(1),F(9,4),F(1,3),F(2,7)])
   d=rng.randint(2,len(active)) if z%3==0 else 2
   kids=[]
   for _ in range(d):kids.append(active.pop(rng.randrange(len(active))))
   active.append((age,kids,'I'));events.append((age,d))
  expected=expectation(n,events)
  if expected is None:continue
  for mode,scale,subdivide in [(0,F(1),False),(1,F(1,8),False),(2,F(8),True),(0,F(1,10**100),False),(2,F(10**100),True)]:
   text=serialize(active[0],labelmode=mode,scale=scale,subdivide=subdivide)
   if z%2:text='(('+text+':123)U:1e100)ROOT'
   else:text+=':1e100'
   cases.append({'name':f'n{n}-history{z}-mode{mode}','newick':text+';','expected':expected,'events':[[str(a),d] for a,d in events]})
cases.extend([{'name':'decisive','newick':'((A:1,B:1)I:9,C:10)R;','expected':5*math.sqrt(3)/7}, {'name':'opposite','newick':'((A:9,B:9)I:1,C:10)R;','expected':-25*math.sqrt(3)/29}, {'name':'zero','newick':'((A:2,B:2)I:3,C:5)R;','expected':0}, {'name':'star','newick':'(A:3,B:3,C:3,D:3);','expected':-math.sqrt(6)}])
pre=integrity();(P/'prehashes.json').write_text(json.dumps(pre,indent=2));
if pre['failures']:raise RuntimeError(pre)
a=inventory(R/'upstream');b=inventory(R/'candidate')
delta={'baseline_count':len(a),'candidate_count':len(b),'changed':[k for k in a.keys()&b.keys() if a[k]!=b[k]],'added':sorted(b.keys()-a.keys()),'removed':sorted(a.keys()-b.keys()),'baseline_files':a,'candidate_files':b}
def strip_gamma(path):
 tree=ast.parse(path.read_text())
 for cls in tree.body:
  if isinstance(cls,ast.ClassDef) and cls.name=='Tree':cls.body=[x for x in cls.body if not isinstance(x,(ast.FunctionDef,ast.AsyncFunctionDef)) or x.name!='gamma_statistic']
 return ast.dump(tree,include_attributes=False)
delta['ast_outside_gamma_equal']=strip_gamma(R/'upstream/treeswift/Tree.py')==strip_gamma(R/'candidate/treeswift/Tree.py')
(P/'inventory.json').write_text(json.dumps(delta,indent=2));
if delta['changed']!=['treeswift/Tree.py'] or delta['added']!=['.tests/test_gamma_statistic.py'] or not delta['ast_outside_gamma_equal']:raise RuntimeError('scope mismatch')
results={}
for source in ['upstream','candidate']:
 mod=load(R/source);fails=[];mutations=[];matches=0
 for c in cases:
  tree=mod.read_tree_newick(c['newick']);before=tree.newick()
  try:
   value=tree.gamma_statistic()
   if math.isclose(value,c['expected'],rel_tol=2e-12,abs_tol=2e-12):matches+=1
   else:fails.append(dict(c,actual=value))
  except Exception as e:fails.append(dict(c,error=type(e).__name__,message=str(e)))
  if tree.newick()!=before:mutations.append(c['name'])
 boundaries=[]
 for text in ['A;','(A:1,B:1);','(A:1,B:1)R:17;','((A:0,B:0)I:0,C:0)R;','((A:1,B:1)I,C:2)R;','((A:1,B:1)I:1,C:3)R;','((A:-1,B:-1)I:2,C:1)R;']:
  try:boundaries.append({'tree':text,'result':mod.read_tree_newick(text).gamma_statistic()})
  except Exception as e:boundaries.append({'tree':text,'exception':type(e).__name__,'message':str(e)})
 results[source]={'tested':len(cases),'matches':matches,'failures':fails,'input_mutations':mutations,'boundaries':boundaries}
(P/('results-optimized.json' if sys.flags.optimize else 'results.json')).write_text(json.dumps(results,indent=2))

if results['candidate']['failures'] or results['candidate']['input_mutations']:raise RuntimeError('candidate failure')
post=integrity();(P/'posthashes.json').write_text(json.dumps(post,indent=2));
if post['failures']:raise RuntimeError(post)
print(json.dumps({'optimized':sys.flags.optimize,'snapshot':{k:v for k,v in post.items() if k!='actual_files'},'cases':len(cases),'matches':{k:v['matches'] for k,v in results.items()},'scope':{k:v for k,v in delta.items() if not k.endswith('files')}}))
