"""Portable replay of the frozen 714 exact-rational fixtures."""
from pathlib import Path
import importlib,json,math,os,subprocess,sys,tempfile
P=Path(__file__).resolve().parent
ROOT=P.parent
from make_cases import cases
from independent_oracle import oracle_newick
expected=json.loads((P/'exact-cases.json').read_text())
assert cases()==expected
output={'generated_case_count':len(expected),'tolerance':{'absolute':2e-13,'relative':2e-13},'sources':{}}
for source in ['upstream','candidate']:
    for key in list(sys.modules):
        if key=='treeswift' or key.startswith('treeswift.'):del sys.modules[key]
    sys.path.insert(0,str(ROOT/'runtime'/source));mod=importlib.import_module('treeswift');sys.path.pop(0)
    rows=[]
    for case in expected:
        assert oracle_newick(case['newick'])==case['oracle']
        tree=mod.read_tree_newick(case['newick']);before=tree.newick()
        row={'id':case['id'],'exact_expectation':case['oracle']}
        try:
            value=tree.gamma_statistic();row['actual']=value
            row['matches']=math.isclose(value,case['oracle']['value'],rel_tol=2e-13,abs_tol=2e-13)
        except Exception as exc:
            row.update(exception=type(exc).__name__,message=str(exc),matches=False)
        row['input_unchanged']=before==tree.newick();rows.append(row)
    summary={'matches':sum(x['matches'] for x in rows),'numeric_mismatches':sum(not x['matches'] and 'exception' not in x for x in rows),'exceptions':sum('exception' in x for x in rows),'mutations':sum(not x['input_unchanged'] for x in rows)}
    if source=='candidate':assert summary=={'matches':714,'numeric_mismatches':0,'exceptions':0,'mutations':0}
    else:assert summary=={'matches':24,'numeric_mismatches':622,'exceptions':68,'mutations':0}
    output['sources'][source]={'summary':summary,'cases':rows}
# Execute author-facing regressions and existing smoke suite in external CWD.
with tempfile.TemporaryDirectory(prefix='treeswift-tests-') as temp:
    for source in ['upstream','candidate']:
        env=os.environ.copy();env['PYTHONPATH']=str(ROOT/'runtime'/source);env['PYTHONOPTIMIZE']='0';env['PYTHONDONTWRITEBYTECODE']='1'
        focused=subprocess.run([sys.executable,'-B',str(ROOT/'runtime/candidate/.tests/test_gamma_statistic.py')],cwd=temp,env=env,capture_output=True,text=True)
        existing=subprocess.run([sys.executable,'-B',str(ROOT/'runtime'/source/'.tests/tests.py')],cwd=temp,env=env,capture_output=True,text=True)
        assert focused.returncode==(1 if source=='upstream' else 0)
        assert existing.returncode==0
(P/'owner-certificate.json').write_text(json.dumps(output,sort_keys=True,indent=2)+'\n')
