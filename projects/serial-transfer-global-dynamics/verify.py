#!/usr/bin/env python3
"""Read-only archive/hash and two-implementation certificate replay."""
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


def replay(script, archived):
    # Deliberately keep scientific child assertions enabled when testing the
    # wrapper under -O or in an environment that sets PYTHONOPTIMIZE.
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PYTHONOPTIMIZE="0")
    completed = subprocess.run([sys.executable, "-B", str(ROOT/script)],
                               capture_output=True, text=True, check=True,
                               timeout=180, env=env)
    actual = json.loads(completed.stdout)
    expected = json.loads((ROOT/archived).read_text())
    require(actual == expected, f"Certificate mismatch: {script}")
    return actual


def verify():
    manifest = json.loads((ROOT/"archive_manifest.json").read_text())
    require(manifest.get("schema") == "serial-transfer-archive-v1", "Wrong manifest schema")
    expected = manifest["files"]
    actual = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob("*")
              if p.is_file() and p != ROOT/"archive_manifest.json"}
    require(actual == set(expected), "Archive file coverage differs")
    for name, digest in expected.items():
        path = ROOT/name
        require(not path.is_symlink() and path.resolve().is_relative_to(ROOT.resolve()),
                "Unsafe archive path")
        require(hashlib.sha256(path.read_bytes()).hexdigest() == digest,
                f"Hash mismatch: {name}")
    owner = replay("science/owner_checker.py", "evidence/owner-certificate.json")
    independent = replay("science/independent_checker.py", "evidence/independent-certificate.json")
    return {"schema": "serial-transfer-replay-v1", "files_verified": len(expected),
            "owner_exact_transfers": owner["exact"]["exact_constructed_transfers"],
            "owner_floating_transfers": owner["numerical"]["transfers_replayed"],
            "independent_certificate_reproduced": bool(independent),
            "scope": "Original analytic proof plus two finite computational checks; no human review or priority certification."}


if __name__ == "__main__":
    require(len(sys.argv) == 2 and sys.argv[1] == "verify", "Usage: verify.py verify")
    print(json.dumps(verify(), sort_keys=True, indent=2))
