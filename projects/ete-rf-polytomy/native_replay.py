"""Compile pinned official source twice and execute the original audit workers."""
from pathlib import Path
import argparse,sys,os,json,hashlib,urllib.request,tarfile,shutil,subprocess
ROOT=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output-dir',required=True);ap.add_argument('--source-archive');a=ap.parse_args()
    raw=Path(a.output_dir).expanduser()
    if not raw.is_absolute() or any(p.is_symlink() for p in [raw,*raw.parents]):raise RuntimeError('invalid output')
    work=raw.resolve()
    if work==ROOT or ROOT in work.parents or work.exists():raise RuntimeError('output must be new and external')
    work.mkdir();(work/'evidence').mkdir()
    spec=json.loads((ROOT/'provenance/source-mapping.json').read_text());archive=work/'source.tar.gz'
    if a.source_archive:shutil.copyfile(Path(a.source_archive),archive)
    else:
        req=urllib.request.Request(spec['source_url'],headers={'User-Agent':'ETE-RF-reproducibility/1.0'})
        with urllib.request.urlopen(req,timeout=60) as r,archive.open('wb') as f:shutil.copyfileobj(r,f)
    if sha(archive)!=spec['source_archive_sha256']:raise RuntimeError('official archive checksum mismatch')
    extract=work/'extract';extract.mkdir()
    with tarfile.open(archive) as tar:tar.extractall(extract,filter='data')
    children=list(extract.iterdir())
    if len(children)!=1 or not children[0].is_dir():raise RuntimeError('unexpected source archive')
    shutil.move(str(children[0]),work/'baseline');extract.rmdir();shutil.copytree(work/'baseline',work/'candidate')
    if sha(work/'baseline/ete4/core/tree.pyx')!=spec['production_baseline_sha256']:raise RuntimeError('baseline source mismatch')
    env=os.environ.copy();env.pop('PYTHONOPTIMIZE',None);env['PYTHONDONTWRITEBYTECODE']='1';env['PIP_DISABLE_PIP_VERSION_CHECK']='1'
    steps=[]
    def run(name,cmd,cwd=work,expected=0):
        q=subprocess.run(cmd,cwd=cwd,env=env,capture_output=True,text=True,timeout=300)
        (work/'evidence'/(name+'.log')).write_text(q.stdout+'\n'+q.stderr)
        steps.append({'name':name,'returncode':q.returncode})
        if q.returncode!=expected:raise RuntimeError(name+' failed; inspect external log')
        return q
    run('apply-patch',['patch','-p1','--batch','--forward','-i',str(ROOT/'upstream.patch')],work/'candidate')
    for rel,expected in spec['patched_files'].items():
        if sha(work/'candidate'/rel)!=expected:raise RuntimeError('patched source mismatch: '+rel)
    for variant in ['baseline','candidate']:
        run('build-'+variant,[sys.executable,'-B','-m','pip','install','--no-deps','--no-build-isolation','--target',str(work/(variant+'-runtime')),str(work/variant)])
    for public,original in spec['owner_worker_mapping'].items():shutil.copyfile(ROOT/'science'/public,work/original)
    run('baseline-regressions',[sys.executable,'-B',str(work/'run_upstream_tests.py'),str(work/'baseline-runtime')],expected=1)
    baseline=json.loads((work/'evidence/upstream-tests-baseline-runtime.json').read_text())
    if baseline['failure_count']!=10 or baseline['error_count']!=8:raise RuntimeError('baseline failure signature changed')
    run('candidate-regressions',[sys.executable,'-B',str(work/'run_upstream_tests.py'),str(work/'candidate-runtime')])
    candidate=json.loads((work/'evidence/upstream-tests-candidate-runtime.json').read_text())
    if not candidate['successful'] or candidate['tests_run']!=6:raise RuntimeError('candidate regressions failed')
    run('independent-native-audit',[sys.executable,'-B',str(work/'independent_audit.py'),str(work/'candidate-runtime')])
    exact=json.loads((work/'evidence/independent-candidate-runtime.json').read_text())
    count=sum(v for k,v in exact['cases'].items() if k!='expected_both_induced_polytomy_errors')
    if count!=189204 or exact['failures']:raise RuntimeError('independent native audit mismatch')
    for variant in ['baseline','candidate']:run('support-'+variant,[sys.executable,'-B',str(work/'support_residual.py'),str(work/(variant+'-runtime'))])
    support=[json.loads((work/'evidence'/('support-residual-'+v+'-runtime.json')).read_text())['cases'] for v in ['baseline','candidate']]
    if support[0]!=support[1] or len(support[0])!=24:raise RuntimeError('support residual changed')
    result={'kind':'real-current-source-native-replay','commit':spec['current_commit'],'baseline_signature':{'failures':10,'errors':8},'candidate_focused_methods_passed':6,'independent_API_cases':count,'expected_both_induced_poly_exceptions':exact['cases']['expected_both_induced_polytomy_errors'],'support_residual_cases_unchanged':24,'support_filter_correctness_claimed':False,'whole_upstream_suite_claimed':False,'steps':steps}
    (work/'native-result.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print(json.dumps(result,sort_keys=True))
if __name__=='__main__':main()
