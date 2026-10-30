import itertools, json, pathlib, hashlib
import cogent3.evolve.fast_distance as distance
from fractions import Fraction
import cogent3
from cogent3.evolve.fast_distance import get_distance_calculator
seqs={"a":"A-","b":"AA","c":"AC"}
def exact(s,t):
    v=[(a,b) for a,b in zip(s,t) if a in "ACGT" and b in "ACGT"]
    return str(Fraction(sum(a!=b for a,b in v),len(v)))
rows=[]
for order in [("a","b","c"),("b","c","a")]:
    aln=cogent3.make_aligned_seqs({k:seqs[k] for k in order},moltype="dna")
    calc=get_distance_calculator("pdist",moltype="dna");calc(aln)
    legacy=calc.get_pairwise_distances()
    app=cogent3.get_app("fast_slow_dist",distance="pdist",moltype="dna")(aln)
    modern=aln.distance_matrix(calc="pdist")
    rows.append({"order":order,"pairs":[{"pair":[a,b],"exact":exact(seqs[a],seqs[b]),"legacy":float(legacy[a,b]),"app":float(app[a,b]),"modern":float(modern[a,b])} for a,b in itertools.combinations(seqs,2)]})
print(json.dumps({"package_version":cogent3.__version__,"source_path":str(pathlib.Path(distance.__file__)),"source_sha256":hashlib.sha256(pathlib.Path(distance.__file__).read_bytes()).hexdigest(),"rows":rows},indent=2))
