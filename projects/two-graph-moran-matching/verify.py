#!/usr/bin/env python3
"""Read-only verification of historical m=2 and all-matching certificates."""
from pathlib import Path
import argparse, hashlib, json, os, shutil, subprocess, sys, tempfile
ROOT=Path(__file__).resolve().parent
MANIFEST='archive_manifest.json'
CERTS=('evidence/owner-certificate.json','evidence/independent-certificate.json',
       'evidence/all-matching-owner-certificate.json','evidence/all-matching-labelled-certificate.json',
       'evidence/all-matching-independent-certificate.json')
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
    import sympy, mpmath
    require(sympy.__version__=='1.13.3','Use pinned SymPy 1.13.3')
    require(mpmath.__version__=='1.3.0','Use pinned mpmath 1.3.0')
    temporary_root=Path(tempfile.gettempdir()).resolve()
    require(temporary_root!=ROOT and ROOT not in temporary_root.parents,'Temporary workers must be outside project')
    with tempfile.TemporaryDirectory(prefix='moran-matching-',dir=temporary_root) as folder:
        work=Path(folder)
        env=os.environ.copy();env['PYTHONOPTIMIZE']='0';env['PYTHONDONTWRITEBYTECODE']='1'
        env.pop('MATH_JOB_ID',None)
        computed={}
        jobs=(('owner','science/owner_checker.py','verify_certificate.py','evidence/six-polynomial-certificates.json',CERTS[0]),
              ('fresh','science/independent_checker.py','independent_checker.py','results.json',CERTS[1]),
              ('all-owner','science/all_matching_owner.py','certificate.py','evidence/uniform-certificate.json',CERTS[2]),
              ('all-labelled','science/all_matching_labelled.py','labelled_check.py','evidence/labelled-check.json',CERTS[3]),
              ('all-fresh','science/all_matching_independent.py','independent_check.py','results.json',CERTS[4]))
        for label,source,script,output,destination in jobs:
            job=work/label;job.mkdir();(job/'evidence').mkdir()
            shutil.copyfile(ROOT/source,job/script)
            if label in ('all-labelled','all-fresh'):
                (job/'evidence/uniform-certificate.json').write_bytes(computed[CERTS[2]])
            run=subprocess.run([sys.executable,'-B',str(job/script)],cwd=job,env=env,capture_output=True,text=True)
            require(run.returncode==0,'Worker failed '+label+': '+run.stderr+run.stdout)
            require((job/output).is_file(),'Missing worker output: '+label)
            computed[destination]=(job/output).read_bytes()
        require(set(computed)==set(CERTS),'Incomplete certificate computation')
        return computed
def verify():
    manifest=check_archive();computed=compute()
    for name,data in computed.items():require(data==(ROOT/name).read_bytes(),'Certificate mismatch: '+name)
    check_archive()
    return {'ok':True,'mode':'verify','archive_modified':False,'public_files_checked':len(manifest['files']),
            'certificates':{name:sha(data) for name,data in computed.items()},'worker_python_optimization':0}
def regenerate(target):
    require(target.is_absolute(),'Output must be absolute')
    require(not target.is_symlink(),'Output symlink forbidden')
    target=target.resolve()
    require(target!=ROOT and ROOT not in target.parents,'Output must be outside project')
    require(not target.exists(),'Output must be a new exclusive directory')
    check_archive();computed=compute()
    files={name:sha(computed[name] if name in computed else (ROOT/name).read_bytes()) for name in inventory()}
    manifest={'schema':'two-graph-moran-matching-archive-v1','algorithm':'sha256','files':files}
    target.mkdir(parents=False,exist_ok=False)
    for name,data in computed.items():
        path=target/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
    (target/MANIFEST).write_bytes(encode(manifest))
    return {'ok':True,'mode':'regenerate','archive_modified':False,'certificates':{name:sha(data) for name,data in computed.items()}}
def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='mode',required=True)
    sub.add_parser('verify');r=sub.add_parser('regenerate');r.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    try:result=verify() if args.mode=='verify' else regenerate(args.output_dir)
    except Exception as exc:
        print(json.dumps({'ok':False,'error':str(exc)},sort_keys=True));return 1
    print(json.dumps(result,sort_keys=True));return 0
if __name__=='__main__':sys.exit(main())
