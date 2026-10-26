"""Independent degree-count oracle for actual compiled ETE RF calls.

No patched helper is imported. Trees are generated as leaf-labelled hierarchies.
The oracle prunes graphs explicitly and computes unresolved degree excess.
"""
from pathlib import Path
from functools import lru_cache
import sys,itertools,json,collections,hashlib,time
R=Path(__file__).resolve().parent
runtime=Path(sys.argv[1]).resolve()
sys.path.insert(0,str(runtime))
sys.path.insert(1,str((R.parent/'compbio-primary-triage/evidence/runtime').resolve()))
from ete4 import Tree
from ete4.core.tree import TreeError
import ete4
if not Path(ete4.__file__).resolve().is_relative_to(runtime):
    raise RuntimeError('wrong compiled package imported')

def partitions(items):
    if not items:
        yield ();return
    first,*tail=items
    for partition in partitions(tuple(tail)):
        yield ((first,),)+partition
        for i,block in enumerate(partition):
            yield partition[:i]+((first,)+block,)+partition[i+1:]

@lru_cache(None)
def hierarchies(labels):
    if len(labels)==1:return (labels[0],)
    found=set()
    for part in partitions(labels):
        if len(part)<2:continue
        for choices in itertools.product(*(hierarchies(tuple(sorted(block))) for block in part)):
            found.add(tuple(sorted(choices,key=repr)))
    return tuple(sorted(found,key=repr))

def newick(t):
    return t if isinstance(t,str) else '('+','.join(map(newick,t))+')'

def leafset(t):
    return {t} if isinstance(t,str) else set().union(*(leafset(c) for c in t))

def rooted_induce(t,common):
    if isinstance(t,str):return t if t in common else None
    children=tuple(x for c in t if (x:=rooted_induce(c,common)) is not None)
    if not children:return None
    return children[0] if len(children)==1 else children

def rooted_oracle(t,common):
    t=rooted_induce(t,common)
    if t is None or isinstance(t,str):return set(),0
    clusters=set();excess=0
    def visit(node):
        nonlocal excess
        if isinstance(node,str):return {node}
        excess+=len(node)-2
        leaves=set().union(*(visit(c) for c in node))
        if 1<len(leaves)<len(common):clusters.add(frozenset(leaves))
        return leaves
    visit(t)
    return clusters,excess

def unrooted_oracle(t,common):
    graph={};names={};counter=iter(range(10000))
    def build(node):
        v=next(counter);graph[v]=set()
        if isinstance(node,str):names[v]=node
        else:
            for child in node:
                c=build(child);graph[v].add(c);graph[c].add(v)
        return v
    build(t)
    retained={v for v,name in names.items() if name in common}
    def remove(v):
        for n in list(graph[v]):graph[n].remove(v)
        del graph[v]
    changed=True
    while changed:
        changed=False
        for v in list(graph):
            if v not in retained and len(graph[v])<2:
                remove(v);changed=True
    # Suppress degree-two vertices of the actual induced unrooted graph.
    changed=True
    while changed:
        changed=False
        for v in list(graph):
            if v not in retained and len(graph[v])==2:
                a,b=graph[v];remove(v);graph[a].add(b);graph[b].add(a);changed=True
    splits=set()
    for a in graph:
        for b in graph[a]:
            seen={b};todo=[a];side=set()
            while todo:
                v=todo.pop()
                if v in seen:continue
                seen.add(v)
                if v in retained:side.add(names[v])
                todo.extend(graph[v]-seen)
            if 1<len(side)<len(common)-1:
                splits.add(frozenset((frozenset(side),frozenset(common-side))))
    excess=sum(max(0,len(graph[v])-3) for v in graph if v not in retained)
    return splits,excess

counts=collections.Counter();failures=[];examples=[];start=time.monotonic()
def check(a,b,unrooted,label,mutability=False):
    sa,sb=newick(a)+';',newick(b)+';';common=leafset(a)&leafset(b)
    oracle=unrooted_oracle if unrooted else rooted_oracle
    p1,d1=oracle(a,common);p2,d2=oracle(b,common)
    expected_error=d1>0 and d2>0
    expected=len(p1^p2)-max(d1,d2)
    first,second=Tree(sa),Tree(sb)
    before=(first.write(),second.write()) if mutability else None
    try:
        actual=first.robinson_foulds(second,unrooted_trees=unrooted,
                                    correct_by_polytomy_size=True)
        ok=not expected_error and actual[0]==expected and actual[1]==len(p1)+len(p2)
        detail={'expected':expected,'rf':actual[0],'max_expected':len(p1)+len(p2),'max_rf':actual[1]}
    except TreeError as ex:
        ok=expected_error and 'Both trees contain polytomies' in str(ex)
        detail={'error':str(ex)}
    except Exception as ex:
        ok=False;detail={'unexpected_error':type(ex).__name__+':'+str(ex)}
    if mutability and before!=(first.write(),second.write()):ok=False;detail['mutated']=True
    counts[label]+=1
    if expected_error:counts['expected_both_induced_polytomy_errors']+=1
    if not ok:
        counts['failures']+=1
        if len(failures)<12:failures.append(dict(a=sa,b=sb,unrooted=unrooted,d1=d1,d2=d2,**detail))
    if len(examples)<3 and not ok:examples.append(detail)

# Every unordered rooted hierarchy on n distinct fixed labels, with no unary
# vertices; every pair. Viewed unrooted, these include all master-root placements
# at vertices and edges (degree-two roots), including all genuine polytomies.
family_counts={}
for n in range(2,6):
    alltrees=hierarchies(tuple('ABCDE'[:n]));rooted=[t for t in alltrees if len(t)==2]
    family_counts[str(n)]={'all_hierarchies':len(alltrees),'rooted_API_accepted_root_degree2':len(rooted)}
    for a in alltrees:
        for b in alltrees:check(a,b,True,'same_labels_unrooted_n'+str(n),n<=3)
    for a in rooted:
        for b in rooted:check(a,b,False,'same_labels_rooted_n'+str(n),n<=3)
# Explicit unequal taxa: all 26 rooted hierarchies of ABCD against each of the
# 236 hierarchies of ABCDE, in both orientations. This covers all placements of
# the removed leaf E, not only one hand-picked ghost-polytomy case.
small=hierarchies(tuple('ABCD'));large=hierarchies(tuple('ABCDE'))
for a in small:
    for b in large:
        check(a,b,True,'different_labels_unrooted')
        check(b,a,True,'different_labels_unrooted')
        if len(a)==2 and len(b)==2:
            check(a,b,False,'different_labels_rooted')
            check(b,a,False,'different_labels_rooted')
# Non-nested unequal leaf sets: all placements of the noncommon E and F.
def rename(t):
    return ('F' if t=='E' else t) if isinstance(t,str) else tuple(rename(c) for c in t)
for a in large:
    for b0 in large:
        b=rename(b0)
        check(a,b,True,'nonnested_labels_unrooted')
        if len(a)==2 and len(b)==2:check(a,b,False,'nonnested_labels_rooted')
# Full label permutations of the three presentations of one unrooted quartet.
forms=[(('A','B'),('C','D')),('A','B',('C','D')),(('A','B'),'C','D')]
for labels in itertools.permutations('ABCD'):
    mapping=dict(zip('ABCD',labels))
    def perm(t):return mapping[t] if isinstance(t,str) else tuple(perm(c) for c in t)
    for a in forms:
        for b in forms:check(perm(a),perm(b),True,'quartet_label_permutations',True)
# Unary paths/custom-attribute cases are exercised in upstream unit tests.
result={'runtime':str(runtime),'actual_package':ete4.__file__,'oracle':'independent explicit graph pruning/degree-two suppression and degree-excess counts; rooted hierarchy induction','domain_counts':family_counts,'cases':dict(counts),'failures':failures,'seconds':time.monotonic()-start,'finite_only':True,'nonzero_support_filters_not_verified_here':True}
(R/'evidence'/('independent-'+runtime.name+'.json')).write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result),flush=True)
if failures:raise SystemExit(1)
