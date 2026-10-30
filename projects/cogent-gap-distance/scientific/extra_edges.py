import json,pathlib,hashlib,sys,itertools,warnings
from fractions import Fraction
import numpy as np
import cogent3
from cogent3.evolve import fast_distance as fd
warnings.filterwarnings('ignore',category=RuntimeWarning)
rows=[]
for mt,canonical,seqs in [('dna','ACGT',{'a':'-N','b':'A-','c':'-A'}),('rna','ACGU',{'a':'-N','b':'U-','c':'-U'}),('protein','ACDEFGHIKLMNPQRSTVWY',{'a':'-X','b':'M-','c':'-M'})]:
 for method in ['pdist','hamming','jc69','tn93','logdet','paralinear']:
  if mt=='protein' and method in ['jc69','tn93']:continue
  for order in itertools.permutations(seqs):
   aln=cogent3.make_aligned_seqs({n:seqs[n] for n in order},moltype=mt)
   c=fd.get_distance_calculator(method,moltype=mt,invalid_raises=True)
   try:
    with np.errstate(all='ignore'):c(aln)
    rows.append({'case':'no_overlap_raise','moltype':mt,'method':method,'order':order,'ok':False})
   except ArithmeticError:rows.append({'case':'no_overlap_raise','moltype':mt,'method':method,'order':order,'ok':True})
# supported wider base moltypes must not be silently rejected or changed
for mt in ['text','bytes']:
 for method in ['pdist','hamming']:
  for order in itertools.permutations('abc'):
   seqs={'a':'ABCD','b':'ABCD','c':'ABCE'}
   try:
    aln=cogent3.make_aligned_seqs({n:seqs[n] for n in order},moltype=mt)
    c=fd.get_distance_calculator(method,moltype=mt);c(aln)
    expected=.25 if method=='pdist' else 1
    rows.append({'case':'wide_moltype','moltype':mt,'method':method,'order':order,'ok':bool(c.dists['a','b']==0 and c.dists['a','c']==expected and c.dists['b','c']==expected)})
   except Exception as e:rows.append({'case':'wide_moltype','moltype':mt,'method':method,'order':order,'exception':repr(e)})
# logdet adjustment=False is also a public option. Count direct mismatches for a near-complete DNA profile.
seqs={'a':'ACGT-ACGT','b':'ACGTAACGT','c':'ACGTGACGT'}
for order in itertools.permutations(seqs):
 aln=cogent3.make_aligned_seqs({n:seqs[n] for n in order},moltype='dna')
 c=fd.LogDetPair(moltype='dna',use_tk_adjustment=False);c(aln)
 values={}
 for a,b in itertools.combinations(seqs,2):
  pair=cogent3.make_aligned_seqs({a:seqs[a],b:seqs[b]},moltype='dna');x=fd.LogDetPair(moltype='dna',use_tk_adjustment=False);x(pair)
  values[a+b]=float(c.dists[a,b])
  rows.append({'case':'logdet_no_adjust','order':order,'pair':[a,b],'ok':bool(np.isclose(c.dists[a,b],x.dists[a,b],equal_nan=True))})
# Explicit public API witness receipt in both important orders, source path pinned.
witness=[]
for order in [('a','b','c'),('b','c','a')]:
 seqs={'a':'A-','b':'AA','c':'AC'};aln=cogent3.make_aligned_seqs({n:seqs[n] for n in order},moltype='dna')
 for method in ['pdist','hamming']:
  c=fd.get_distance_calculator(method,moltype='dna');c(aln)
  app=cogent3.get_app('fast_slow_dist',distance=method,moltype='dna')(aln)
  witness.append({'order':order,'method':method,'legacy':{a+b:float(c.dists[a,b]) for a,b in itertools.combinations(seqs,2)},'app':{a+b:float(app[a,b]) for a,b in itertools.combinations(seqs,2)}})
print(json.dumps({'imported':str(pathlib.Path(fd.__file__).resolve()),'source_sha256':hashlib.sha256(pathlib.Path(fd.__file__).read_bytes()).hexdigest(),'checks':len(rows),'checks_failed':sum(x.get('ok') is False for x in rows),'exceptions':sum('exception' in x for x in rows),'rows':rows,'witness':witness},indent=2))
