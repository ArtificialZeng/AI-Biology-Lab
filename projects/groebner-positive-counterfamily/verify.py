#!/usr/bin/env python3
"""Read-only manifest and exact-certificate replay."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def replay(script, certificate):
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PYTHONOPTIMIZE="0")
    run = subprocess.run([sys.executable, "-B", str(ROOT/script)],
                         capture_output=True, text=True, check=True,
                         timeout=180, env=env)
    value = json.loads(run.stdout)
    require(value == json.loads((ROOT/certificate).read_text()), f"Certificate mismatch: {script}")
    return value


def verify():
    manifest = json.loads((ROOT/"archive_manifest.json").read_text())
    require(manifest.get("schema") == "groebner-positive-archive-v1", "Manifest schema differs")
    actual = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob("*")
              if p.is_file() and p != ROOT/"archive_manifest.json"}
    require(actual == set(manifest["files"]), "Archive file coverage differs")
    for name, expected in manifest["files"].items():
        path = ROOT/name
        require(not path.is_symlink() and path.resolve().is_relative_to(ROOT.resolve()), "Unsafe path")
        require(hashlib.sha256(path.read_bytes()).hexdigest() == expected, f"Hash mismatch: {name}")
    owner = replay("science/owner_checker.py", "evidence/owner-certificate.json")
    independent = replay("science/independent_checker.py", "evidence/independent-certificate.json")
    require(bool(independent), "Empty independent certificate")
    return {"schema": "groebner-positive-replay-v1", "files_verified": len(actual),
            "owner_exact_positive_cases": owner["exact_positive_cases"],
            "complete_unordered_partitions": owner["unordered_partitions_checked"],
            "independent_certificate_reproduced": True,
            "scope": "Original analytic proof with exact supporting checks; no human review, first-disproof or priority certification."}


if __name__ == "__main__":
    require(len(sys.argv) == 2 and sys.argv[1] == "verify", "Usage: verify.py verify")
    print(json.dumps(verify(), sort_keys=True, indent=2))
