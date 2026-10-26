"""Offline verification. Explicit pip setup is separate from verification."""
import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent
BASE = "6088541cdcb01c9c351651244534037260dfae987cc5995c9ba3b41cab6dff4a"
CAND = "643e54b3cad43d4b31f99fb6ffcd9ebbf3761a8ae582a631569988c661fbdce7"

def canonical(obj):
    return (json.dumps(obj, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def inventory(root):
    rows = {}
    for path in root.rglob("*"):
        if path.is_symlink():
            raise ValueError("Symlink in public package")
        if path.is_file() and path != root / "archive_manifest.json":
            rows[path.relative_to(root).as_posix()] = sha(path)
    return rows

def integrity():
    manifest = json.loads((ROOT / "archive_manifest.json").read_text())
    if inventory(ROOT) != manifest["files"]:
        raise ValueError("Public input manifest mismatch (missing, extra or changed file)")
    return manifest

def dependencies():
    if sys.version_info[:2] != (3, 13):
        raise RuntimeError("This package pins the tested Python 3.13 runtime")
    for row in (ROOT / "requirements.txt").read_text().splitlines():
        name, expected = row.split("==")
        if importlib.metadata.version(name) != expected:
            raise RuntimeError("Dependency version mismatch: " + name)
    spec = importlib.util.find_spec("cogent3")
    installed = Path(spec.origin).resolve().parent
    if sha(installed / "evolve/fast_distance.py") != BASE:
        raise RuntimeError("Pinned installed package source does not match baseline")
    return installed

def run_actual(overlay, digest, script, temp, mode="script"):
    env = os.environ.copy()
    env.update(PYTHONOPTIMIZE="0", PYTHONDONTWRITEBYTECODE="1", NUMBA_CACHE_DIR=str(temp / "numba"), PYTEST_DISABLE_PLUGIN_AUTOLOAD="1")
    env.pop("PYTHONPATH", None)
    args = [sys.executable, "-B", str(ROOT / "worker.py"), str(overlay), digest, mode, str(ROOT / script)]
    if mode == "pytest":
        args += [os.devnull, str(temp / "pytest"), str(temp / "pytest.xml")]
    done = subprocess.run(args, cwd=temp, env=env, text=True, capture_output=True)
    if done.returncode:
        raise RuntimeError("Scientific worker failed: " + script + "\n" + done.stdout + done.stderr)
    if mode == "pytest":
        suite = ET.parse(temp / "pytest.xml").getroot().find("testsuite")
        vals = {k: int(suite.get(k)) for k in ("tests", "errors", "failures", "skipped")}
        vals["passed"] = vals["tests"] - vals["errors"] - vals["failures"] - vals["skipped"]
        return vals
    value = json.loads(done.stdout)
    if value.get("source_sha256") != digest:
        raise RuntimeError("Scientific receipt source mismatch")
    return value

def calculate():
    installed = dependencies()
    tempbase = Path(tempfile.gettempdir()).resolve()
    if tempbase == ROOT or tempbase.is_relative_to(ROOT):
        raise RuntimeError("Temporary workers must be outside public source")
    with tempfile.TemporaryDirectory(prefix="cogent-audit-", dir=tempbase) as td:
        tmp = Path(td)
        outputs = {}
        for label, digest in (("baseline", BASE), ("candidate", CAND)):
            overlay = tmp / label
            shutil.copytree(installed, overlay / "cogent3", ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.nbc", "*.nbi"))
            shutil.copyfile(ROOT / "vendor" / (label + "_fast_distance.py"), overlay / "cogent3/evolve/fast_distance.py")
            outputs[label] = {name: run_actual(overlay, digest, "scientific/" + filename, tmp) for name, filename in (("exhaustive", "exhaustive_oracle.py"), ("public", "replay_public_api.py"), ("fresh", "independent_oracle.py"), ("edges", "extra_edges.py"))}
            if label == "candidate":
                tests = run_actual(overlay, digest, "tests/test_distance_missing.py", tmp, "pytest")
        b, c = outputs["baseline"], outputs["candidate"]
        for receipt in c["public"]["rows"]:
            for pair in receipt["pairs"]:
                numerator, _, denominator = pair["exact"].partition("/")
                exact = int(numerator) / (int(denominator) if denominator else 1)
                if any(pair[k] != exact for k in ("legacy", "app", "modern")):
                    raise RuntimeError("Actual public candidate value disagrees with exact count")
        if not any(pair["legacy"] != pair["modern"] for row in b["public"]["rows"] for pair in row["pairs"]):
            raise RuntimeError("Baseline counterexample did not reproduce")
        frozen = {json.dumps(x, sort_keys=True) for x in b["fresh"]["failures"]}
        inherited = c["fresh"]["failures"]
        if any(json.dumps(x, sort_keys=True) not in frozen or x.get("exception") != "StopIteration()" for x in inherited):
            raise RuntimeError("New candidate discrepancy or unpaired exception")
        if c["edges"]["checks_failed"] or c["edges"]["exceptions"]:
            raise RuntimeError("Candidate failed edge control")
        if any(x["mismatch_count"] for x in c["exhaustive"]["checks"]):
            raise RuntimeError("Candidate failed exact column oracle")
        return {"schema": 1, "actual_package": "cogent3==2026.10.22", "loaded_source_sha256": {"baseline": BASE, "candidate": CAND}, "public_witness": "baseline discrepancy reproduced; candidate exact in both orders", "exhaustive": {label: [{k: x[k] for k in ("method", "ordered_triples", "checked_pair_values", "mismatch_count")} for x in outputs[label]["exhaustive"]["checks"]] for label in ("baseline", "candidate")}, "fresh": {"baseline_checks": b["fresh"]["checks"], "baseline_discrepancies": b["fresh"]["failures_count"], "candidate_checks": c["fresh"]["checks"], "candidate_inherited_StopIteration": len(inherited), "candidate_new_discrepancies": 0}, "edges": {label: {k: outputs[label]["edges"][k] for k in ("checks", "checks_failed", "exceptions")} for label in ("baseline", "candidate")}, "candidate_regression_tests": tests}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("verify", "regenerate"))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    manifest = integrity()
    if args.mode == "regenerate":
        if args.output is None:
            raise ValueError("Regeneration requires an exclusive external output directory")
        out = args.output.resolve()
        if out == ROOT or out.is_relative_to(ROOT):
            raise ValueError("Regeneration output must be external")
        out.mkdir(parents=False, exist_ok=False)
    actual = canonical(calculate())
    expected = (ROOT / "certificate.json").read_bytes()
    if actual != expected:
        raise ValueError("Actual scientific certificate differs from frozen expectation")
    integrity()
    if args.mode == "regenerate":
        (out / "certificate.json").write_bytes(actual)
        (out / "archive_manifest.json").write_bytes((ROOT / "archive_manifest.json").read_bytes())
    print("Cogent3 audit verified: exact column, public API, paired inherited exceptions and model controls")

if __name__ == "__main__":
    main()
