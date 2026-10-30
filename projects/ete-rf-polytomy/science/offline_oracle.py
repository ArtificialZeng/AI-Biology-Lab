"""Offline graph identities only; this does not execute or simulate ETE."""
from pathlib import Path
from functools import lru_cache
import ast,itertools,json
source=Path(__file__).with_name('owner_native_audit.py')
selected={'partitions','hierarchies','newick','leafset','rooted_induce','rooted_oracle','unrooted_oracle'}
tree=ast.parse(source.read_text())
body=[n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef)) and n.name in selected]
if {n.name for n in body}!=selected:raise RuntimeError('owner function inventory changed')
ns={'lru_cache':lru_cache,'itertools':itertools}
exec(compile(ast.Module(body=body,type_ignores=[]),str(source),'exec'),ns)
counts={};identities=0
for n in range(2,6):
    labels=tuple('ABCDE'[:n]);trees=ns['hierarchies'](labels);counts[str(n)]=len(trees)
    for t in trees:
        for unrooted in [False,True]:
            parts,degree_excess=ns['unrooted_oracle' if unrooted else 'rooted_oracle'](t,set(labels))
            if degree_excess!=max(0,n-(3 if unrooted else 2))-len(parts):raise RuntimeError('degree identity failed')
            identities+=1
forms=[(('A','B'),('C','D')),('A','B',('C','D'))]
parts=[ns['unrooted_oracle'](t,set('ABCD'))[0] for t in forms]
if parts[0]!=parts[1]:raise RuntimeError('quartet topology mismatch')
print(json.dumps({'kind':'offline-exact-graph-identities','hierarchies':counts,'degree_identities':identities,'quartet_expected_RF':0,'native_ETE_API_executed':False},sort_keys=True))
