"""Fresh direct canonical-column reference; no candidate helper used."""
import itertools,json,hashlib,os,pathlib,random,sys,warnings
from fractions import Fraction
import numpy as np
import cogent3
from cogent3.evolve import fast_distance as fd
warnings.filterwarnings('ignore', category=RuntimeWarning)
random.seed(27015)
methods=('pdist','hamming','jc69','tn93','logdet','paralinear')
checks=[]; failures=[]
def same(a,b):
 if a is None or b is None:return a is b
 return bool(np.isclose(a,b,rtol=1e-10,atol=1e-10,equal_nan=True))
def expected_rows(seqs,canon,method):
 result={}
 for a,b in itertools.combinations(seqs,2):
  pairs=[(x,y) for x,y in zip(seqs[a],seqs[b],strict=True) if x in canon and y in canon]
  n=len(pairs); k=sum(x!=y for x,y in pairs)
  if not n: d=float('nan')
  elif method=='pdist': d=float(Fraction(k,n))
  elif method=='hamming': d=k
  elif not k:d=0.
  elif method=='jc69': d=-.75*np.log(1-4*k/(3*n)) if Fraction(k,n)<Fraction(3,4) else float('nan')
  else:continue
  result[(a,b)]=(n,k,d)
 return result

def run_case(seqs,mt,canon,label,permutations=True,app=False):
 for names in itertools.permutations(seqs) if permutations else [tuple(seqs)]:
  aln=cogent3.make_aligned_seqs({n:seqs[n] for n in names},moltype=mt)
  for method in methods:
   if mt=='protein' and method in ('jc69','tn93'):continue
   calc=fd.get_distance_calculator(method,moltype=mt)
   try:
    with np.errstate(all='ignore'):calc(aln)
    for include in (True,False):
     got=calc.get_pairwise_distances(include_duplicates=include)
     for (a,b),(n,k,e) in expected_rows(seqs,canon,method).items():
      if a not in got.names or b not in got.names:continue
      checks.append(1)
      if not same(got[a,b],e):failures.append({'label':label,'order':names,'method':method,'include':include,'seqs':seqs,'pair':[a,b],'expected':e,'actual':float(got[a,b])})
    # independent pair-alone model comparison, including invalid tuples
    for a,b in itertools.combinations(seqs,2):
     if method in ('pdist','hamming','jc69'):continue
     pair=cogent3.make_aligned_seqs({a:seqs[a],b:seqs[b]},moltype=mt)
     solo=fd.get_distance_calculator(method,moltype=mt)
     with np.errstate(all='ignore'):solo(pair)
     checks.append(1)
     if not same(calc.dists[a,b],solo.dists[a,b]):failures.append({'label':label,'method':method,'order':names,'pair':[a,b],'actual':float(calc.dists[a,b]),'pair_alone':float(solo.dists[a,b])})
    if app and method in ('pdist','hamming','jc69'):
     key='fast_calc' if method=='jc69' else 'distance'
     got=cogent3.get_app('fast_slow_dist',moltype=mt,**{key:method})(aln)
     for (a,b),(_,_,e) in expected_rows(seqs,canon,method).items():
      checks.append(1)
      if not same(got[a,b],e):failures.append({'label':label,'app':method,'order':names,'pair':[a,b],'expected':e,'actual':float(got[a,b])})
   except Exception as exc:
    failures.append({'label':label,'order':names,'method':method,'exception':repr(exc),'seqs':seqs})

for mt,canon,profiles in [('dna','ACGT',[('A-','AA','AC'),('AN','AA','AC'),('AR','AG','AT'),('A?','AA','AC'),('--','NN','AC'),('A-','-A','NN'),('ACGTN','ACGT-','ACGA?'),('ACGTN','ACGT-','ACGT?')]),('rna','ACGU',[('U?','UU','UC'),('UN','UU','UC'),('A-','AA','AU'),('--','NN','AC')]),('protein','ACDEFGHIKLMNPQRSTVWY',[('MX','MM','MP'),('M?','MM','MP'),('M-','MM','MP'),('--','XX','AC'),('ACDX','ACD-','ACD?')])]:
 for i,p in enumerate(profiles):run_case(dict(zip('abc',p)),mt,canon,mt+str(i),app=i==0)
 # random wider sequences with distinct full rows; pairwise model invariance vs solo
 symbols=canon+('-N?' if mt!='protein' else '-X?')
 for j in range(18):
  p=[''.join(random.choices(symbols,k=12)) for _ in range(3)]
  run_case(dict(zip('abc',p)),mt,canon,mt+'random'+str(j),permutations=False)
# exact aliases plus other zero-mask rows, and multiple duplicate groups
for mt,canon,seqs in [('dna','ACGT',{'a':'ACGTN','b':'ACGT-','c':'ACGA?','d':'ACGA-'}),('rna','ACGU',{'a':'ACGU?','b':'ACGU-','c':'ACGA?'}),('protein','ACDEFGHIKLMNPQRSTVWY',{'a':'ACDX','b':'ACD-','c':'ACD?'})]:
 run_case(seqs,mt,canon,'aliases',permutations=True)
# Drop invalid semantics on no overlap plus valid finite subset, every applicable metric
for mt,canon in [('dna','ACGT'),('rna','ACGU'),('protein','ACDEFGHIKLMNPQRSTVWY')]:
 for method in methods:
  if mt=='protein' and method in ('jc69','tn93'):continue
  seqs={'a':'ACGTACGT' if mt=='dna' else 'ACGUACGU' if mt=='rna' else 'ACDEACDE','b':'ACGTACGT' if mt=='dna' else 'ACGUACGU' if mt=='rna' else 'ACDEACDE','c':'--------'}
  calc=fd.get_distance_calculator(method,moltype=mt)
  with np.errstate(all='ignore'):calc(cogent3.make_aligned_seqs(seqs,moltype=mt))
  got=calc.dists
  for a,b in [('a','c'),('b','c')]:
   checks.append(1)
   if not np.isnan(got[a,b]):failures.append({'dropinvalid':method,'mt':mt,'pair':[a,b],'value':float(got[a,b])})
  dropped=got.drop_invalid();checks.append(1)
  if dropped is not None and any(not np.isfinite(x) for x in dropped.array.flat): failures.append({'dropinvalid':method,'mt':mt})
result={'imported':str(pathlib.Path(fd.__file__).resolve()),'source_sha256':hashlib.sha256(pathlib.Path(fd.__file__).read_bytes()).hexdigest(),'cogent3':cogent3.__version__,'python':sys.executable,'checks':len(checks),'failures_count':len(failures),'failures':failures,'seed':27015,'scope':'17 explicit profiles x all 6 row orders; 54 independent random length12 triples; 3 multialias profiles incl all24 orders for four rows; public fast_slow_dist 3 profiles x all6 orders; all applicable six metrics and moltypes; drop_invalid no overlap. Exact column oracle pdist/hamming/JC69; model pair-alone invariance secondary oracle.'}
print(json.dumps(result,indent=2))
