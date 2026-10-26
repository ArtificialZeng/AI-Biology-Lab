#!/usr/bin/env python3
"""Read-only relocation-safe manifest and actual-library certificate replay."""
import hashlib,json,os,subprocess,sys,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def check(ok,message):
    if not ok:raise RuntimeError(message)

def replay(script,variant,args=()):
    env=os.environ.copy()
    env["PYTHONPATH"]=str(ROOT/"vendor"/variant)
    env["PYTHONDONTWRITEBYTECODE"]="1"
    # All optional runtime caches go outside the archived project.
    with tempfile.TemporaryDirectory(prefix="mantis-audit-") as temp:
        env.update(MPLCONFIGDIR=temp,NUMBA_CACHE_DIR=temp,XDG_CACHE_HOME=temp)
        flags=["-B"]+(["-O"] if sys.flags.optimize else [])
        command=[sys.executable,*flags,str(ROOT/"scripts"/script),*args]
        out=subprocess.run(command,cwd=temp,env=env,text=True,capture_output=True,timeout=120)
        check(out.returncode==0,"audit failed: "+script+"\n"+out.stderr)
        return json.loads(out.stdout)

def verify():
    manifest=json.loads((ROOT/"archive_manifest.json").read_text())
    actual={f.relative_to(ROOT).as_posix() for f in ROOT.rglob("*") if f.is_file() and f.name!="archive_manifest.json"}
    check(actual==set(manifest["sha256"]),"archive file inventory changed")
    for name,want in manifest["sha256"].items():
        check(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==want,"archive digest mismatch: "+name)
    tasks=[("api_audit.py","baseline",["baseline"],"baseline-api.json"),("api_audit.py","candidate",["candidate"],"candidate-api.json"),("family_audit.py","candidate",[],"family-audit.json"),("reviewer_audit.py","candidate",[],"reviewer-audit.json")]
    for script,variant,args,certificate in tasks:
        got=replay(script,variant,args)
        expected=json.loads((ROOT/"certificates"/certificate).read_text())
        check(got==expected,"certificate replay mismatch: "+certificate)
    print(json.dumps({"verified":True,"source_commit":manifest["source_commit"],"actual_library_api_variants":["baseline","candidate"],"owner_helper_cases":352,"owner_public_api_cases":5,"independent_review_cases":1077,"archive_files":len(actual),"new_math_theorem":False,"external_human_review":False},indent=2,sort_keys=True))

if __name__=="__main__":
    if len(sys.argv)>2 or (len(sys.argv)==2 and sys.argv[1]!="verify"):raise SystemExit("usage: python -B verify.py [verify]")
    verify()
