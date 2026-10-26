#!/usr/bin/env python3
"""Portable read-only verification; regeneration is explicit and external."""
from pathlib import Path
import argparse, hashlib, json, os, shutil, subprocess, sys, tempfile
ROOT=Path(__file__).resolve().parent
MANIFEST='archive_manifest.json'
CERTS=('evidence/owner-certificate.json','evidence/independent-certificate.json')
class VerificationError(Exception): pass
def require(condition,message):
    if not condition: raise VerificationError(message)
def sha(data): return hashlib.sha256(data).hexdigest()
def encode(obj): return (json.dumps(obj,sort_keys=True,indent=2)+'\n').encode()
def inventory():
    out=[]
    for path in ROOT.rglob('*'):
        rel=path.relative_to(ROOT).as_posix()
        require(not path.is_symlink(),'Symlink forbidden: '+rel)
        if path.is_file() and rel!=MANIFEST: out.append(rel)
    return sorted(out)
def check_archive():
    doc=json.loads((ROOT/MANIFEST).read_bytes())
    require(doc.get('schema')=='two-graph-moran-matching-archive-v1','Unknown manifest schema')
    require(doc.get('algorithm')=='sha256','Unknown hash algorithm')
    require(sorted(doc.get('files',{}))==inventory(),'Unexpected or missing archive input')
    for name,expected in doc['files'].items():
        require(sha((ROOT/name).read_bytes())==expected,'Input hash mismatch: '+name)
    return doc
def compute():
    # Original assertions and writes run only in disposable workers, with normal
    # Python explicitly forced even if this public wrapper was invoked with -O.
    import sympy
    require(sympy.__version__=='1.13.3','Use pinned SymPy 1.13.3')
    with tempfile.TemporaryDirectory(prefix='moran-matching-') as folder:
        work=Path(folder); (work/'owner').mkdir(); (work/'fresh').mkdir()
        shutil.copyfile(ROOT/'science/owner_checker.py',work/'owner/verify_certificate.py')
        shutil.copyfile(ROOT/'science/independent_checker.py',work/'fresh/independent_checker.py')
        env=os.environ.copy(); env['PYTHONOPTIMIZE']='0'; env['PYTHONDONTWRITEBYTECODE']='1'
        for relative in ('owner/verify_certificate.py','fresh/independent_checker.py'):
            run=subprocess.run([sys.executable,'-B',str(work/relative)],cwd=work,env=env,capture_output=True,text=True)
            require(run.returncode==0,'Worker failed '+relative+': '+run.stderr+run.stdout)
        return {CERTS[0]:(work/'owner/evidence/six-polynomial-certificates.json').read_bytes(),CERTS[1]:(work/'fresh/results.json').read_bytes()}
def verify():
    manifest=check_archive(); computed=compute()
    for name,data in computed.items(): require(data==(ROOT/name).read_bytes(),'Certificate mismatch: '+name)
    check_archive()
    return {'ok':True,'mode':'verify','archive_modified':False,'public_files_checked':len(manifest['files']),'certificates':{name:sha(data) for name,data in computed.items()},'worker_python_optimization':0}
def regenerate(target):
    require(target.is_absolute(),'Output must be absolute')
    require(not target.is_symlink(),'Output symlink forbidden')
    target=target.resolve()
    require(target!=ROOT and ROOT not in target.parents,'Output must be outside project')
    require(not target.exists(),'Output must be a new exclusive directory')
    check_archive(); computed=compute()
    files={name:sha(computed[name] if name in computed else (ROOT/name).read_bytes()) for name in inventory()}
    manifest={'schema':'two-graph-moran-matching-archive-v1','algorithm':'sha256','files':files}
    target.mkdir(parents=False,exist_ok=False)
    for name,data in computed.items():
        path=target/name; path.parent.mkdir(parents=True,exist_ok=True); path.write_bytes(data)
    (target/MANIFEST).write_bytes(encode(manifest))
    return {'ok':True,'mode':'regenerate','archive_modified':False,'certificates':{name:sha(data) for name,data in computed.items()}}
def main():
    parser=argparse.ArgumentParser(description=__doc__); sub=parser.add_subparsers(dest='mode',required=True)
    sub.add_parser('verify'); r=sub.add_parser('regenerate');r.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    try: result=verify() if args.mode=='verify' else regenerate(args.output_dir)
    except Exception as exc:
        print(json.dumps({'ok':False,'error':str(exc)},sort_keys=True));return 1
    print(json.dumps(result,sort_keys=True));return 0
if __name__=='__main__': sys.exit(main())
