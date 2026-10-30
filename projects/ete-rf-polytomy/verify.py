"""Read-only packet verification; all scientific outputs are external."""
from pathlib import Path
import sys,os,json,hashlib,subprocess,tempfile,shutil,argparse
ROOT=Path(__file__).resolve().parent
MANIFEST='archive_manifest.json'
def fail(message):raise RuntimeError(message)
def inventory(root):
    found={}
    for p in root.rglob('*'):
        if p.is_symlink():fail('symlink in packet')
        if p.is_file():found[p.relative_to(root).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
        elif not p.is_dir():fail('special file in packet')
    return found

def load_unique(text):
    def pairs(values):
        d={}
        for k,v in values:
            if k in d:fail('duplicate JSON key')
            d[k]=v
        return d
    return json.loads(text,object_pairs_hook=pairs)
def verify_inputs(root):
    actual=inventory(root);doc=load_unique((root/MANIFEST).read_text());expected=doc['files']
    if not isinstance(expected,dict):fail('invalid inventory')
    for n,h in expected.items():
        q=Path(n)
        if q.is_absolute() or '..' in q.parts or n!=q.as_posix() or n==MANIFEST or not isinstance(h,str) or len(h)!=64:fail('invalid manifest entry')
    if {k:v for k,v in actual.items() if k!=MANIFEST}!=expected:fail('inventory or byte hash mismatch')
    allowed_dirs={str(Path(n).parent) for n in expected}
    allowed_dirs={str(p) for n in expected for p in Path(n).parents if str(p)!='.'}
    if {p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_dir()}!=allowed_dirs:fail('extra directory')
    return actual

def external(raw):
    p=Path(raw).expanduser()
    if not p.is_absolute():fail('output must be absolute')
    if any(x.is_symlink() for x in [p,*p.parents]):fail('symlink output component')
    p=p.resolve()
    if p==ROOT or ROOT in p.parents:fail('output inside packet')
    if p.exists():fail('output exists')
    if not p.parent.is_dir():fail('output parent missing')
    return p

def worker(root):
    env=os.environ.copy();env.pop('PYTHONOPTIMIZE',None);env['PYTHONDONTWRITEBYTECODE']='1'
    with tempfile.TemporaryDirectory(prefix='ete-offline-oracle-') as work:
        q=subprocess.run([sys.executable,'-B','-W','ignore',str(root/'science/offline_oracle.py')],cwd=work,env=env,capture_output=True,text=True)
    if q.returncode:fail(q.stderr or q.stdout)
    return load_unique(q.stdout)
def main():
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['verify','regenerate','native']);parser.add_argument('--output-dir');parser.add_argument('--source-archive');args=parser.parse_args()
    before=verify_inputs(ROOT)
    if args.command=='verify':
        actual=worker(ROOT);expected=load_unique((ROOT/'evidence/offline-certificate.json').read_text())
        if actual!=expected:fail('offline certificate mismatch')
        if verify_inputs(ROOT)!=before:fail('packet mutated')
        print(json.dumps(actual,sort_keys=True));return
    if not args.output_dir:fail('output directory required')
    out=external(args.output_dir)
    if args.command=='regenerate':
        shutil.copytree(ROOT,out)
        # The archive may be read-only. Only its exclusive external copy becomes writable.
        out.chmod(out.stat().st_mode | 0o700)
        for item in out.rglob('*'):
            item.chmod(item.stat().st_mode | (0o700 if item.is_dir() else 0o600))
        doc=worker(out);(out/'evidence/offline-certificate.json').write_text(json.dumps(doc,sort_keys=True,indent=2)+'\n');hashes=inventory(out);hashes.pop(MANIFEST)
        (out/MANIFEST).write_text(json.dumps({'algorithm':'sha256','files':hashes},sort_keys=True,indent=2)+'\n');verify_inputs(out)
        if verify_inputs(ROOT)!=before:fail('original packet mutated')
        print(json.dumps({'regenerated_complete_packet':True,'files':len(hashes)+1},sort_keys=True));return
    cmd=[sys.executable,'-B',str(ROOT/'native_replay.py'),'--output-dir',str(out)]
    if args.source_archive:cmd+=['--source-archive',args.source_archive]
    env=os.environ.copy();env.pop('PYTHONOPTIMIZE',None);env['PYTHONDONTWRITEBYTECODE']='1'
    q=subprocess.run(cmd,env=env)
    if q.returncode:fail('native replay failed')
    if verify_inputs(ROOT)!=before:fail('packet mutated')
if __name__=='__main__':main()
