#!/usr/bin/env python3
"""Read-only exact witness replay with assertion-enabled scientific children."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
REQUIRED = {"README.md", "LICENSE", "requirements.txt", "sources.md", "status.json",
            "proof.md", "check.py", "witness.json", "reproduce.py"}

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    work = args.work_dir.expanduser().resolve()
    require(work != ROOT and ROOT not in work.parents,
            "Work directory must be outside the project.")
    manifest = json.loads((ROOT / "manifest.json").read_text())
    require(set(manifest["files"]) == REQUIRED, "Unexpected manifest file set.")
    for name, expected in manifest["files"].items():
        target = ROOT / name
        require(not target.is_symlink() and target.resolve().parent == ROOT,
                "Manifest paths must be regular local files.")
        require(target.is_file(), "Missing manifest file: " + name)
        require(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
                "Manifest mismatch: " + name)
    work.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.pop("PYTHONOPTIMIZE", None)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    dependency = subprocess.run(
        [sys.executable, "-I", "-B", "-c",
         'import sys,sympy,json; assert sys.flags.optimize==0; '
         'assert sympy.__version__=="1.13.3"; '
         'print(json.dumps({"python":sys.version.split()[0],'
         '"sympy":sympy.__version__,"optimize":sys.flags.optimize}))'],
        cwd=work, env=env, capture_output=True, text=True)
    (work / "dependency.stdout").write_text(dependency.stdout)
    (work / "dependency.stderr").write_text(dependency.stderr)
    require(dependency.returncode == 0,
            "Dependency/assertion check failed; see external logs.")
    bootstrap = ('import sys,runpy; '
                 'assert sys.flags.optimize==0,"Scientific assertions disabled"; '
                 'target=sys.argv[1]; sys.argv=sys.argv[1:]; '
                 'runpy.run_path(target,run_name="__main__")')
    child = subprocess.run(
        [sys.executable, "-I", "-B", "-c", bootstrap, str(ROOT / "check.py")],
        cwd=work, env=env, capture_output=True, text=True)
    (work / "scientific.stdout").write_text(child.stdout)
    (work / "scientific.stderr").write_text(child.stderr)
    require(child.returncode == 0, "Scientific witness check failed; see external logs.")
    scientific = json.loads(child.stdout)
    require(scientific.get("status") == "exact-witness-checks-passed" and
            scientific.get("directed_edges") == 90 and
            scientific.get("source_and_target_degree_max") == 3 and
            all(scientific.get(k) is True for k in ["all_rates_positive",
                "complete_reversible_graph", "Darboux_identities_exact",
                "tangent_no_real_parameter_roots",
                "fourquarter_Floquet_integrand_strictly_negative",
                "finite_witness_not_global_classification"]),
            "Unexpected scientific checker output.")
    result = {"status": "exact-witness-checks-passed", "manifest_integrity": True,
              "dependency": json.loads(dependency.stdout), "scientific": scientific,
              "scope": "Exact finite witness checks; analytic argument in proof.md; no peer-review or priority claim."}
    (work / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
