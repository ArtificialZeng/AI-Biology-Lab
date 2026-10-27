#!/usr/bin/env python3
"""Read-only archive verification; regeneration targets a new external directory."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

if not __debug__ or os.environ.get("PYTHONOPTIMIZE") != "0":
    environment = dict(os.environ, PYTHONOPTIMIZE="0", PYTHONDONTWRITEBYTECODE="1")
    os.execve(sys.executable, [sys.executable, "-B", str(Path(__file__).resolve()), *sys.argv[1:]], environment)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("command", choices=("verify", "regenerate"))
parser.add_argument("--output", type=Path)
args = parser.parse_args()

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def check_manifest():
    manifest = json.loads((HERE / "archive_manifest.json").read_text())
    for relative, expected in manifest["files"].items():
        file = HERE / relative
        if not file.is_file() or digest(file) != expected:
            raise RuntimeError("Archive hash mismatch: " + relative)
    workflow = ROOT / manifest["workflow"]["path"]
    if digest(workflow) != manifest["workflow"]["sha256"]:
        raise RuntimeError("Workflow hash mismatch")

def recompute(destination):
    for name in ("local_algebra.py", "independent_checker.py"):
        shutil.copyfile(HERE / name, destination / name)
    environment = dict(os.environ, PYTHONOPTIMIZE="0", PYTHONDONTWRITEBYTECODE="1")
    for name in ("local_algebra.py", "independent_checker.py"):
        result = subprocess.run([sys.executable, "-B", str(destination / name)],
                                cwd=destination, env=environment, check=True,
                                capture_output=True, text=True)
        if result.returncode:
            raise RuntimeError("Child computation failed: " + name)

check_manifest()
if args.command == "regenerate":
    if args.output is None:
        parser.error("regenerate requires --output pointing to a new external directory")
    destination = args.output.expanduser().resolve()
    if destination.is_relative_to(ROOT):
        parser.error("regenerate output must be outside this checkout")
    if destination.exists():
        parser.error("regenerate output must not already exist")
    with tempfile.TemporaryDirectory(prefix="six-clock-generate-") as folder:
        temp = Path(folder)
        recompute(temp)
        destination.mkdir(parents=True, exist_ok=False)
        for name in ("local-algebra.json", "independent-certificate.json"):
            shutil.copyfile(temp / name, destination / name)
    check_manifest()
    print(json.dumps({"ok": True, "command": "regenerate", "output": str(destination),
                      "archived_evidence_modified": False, "assertions_enabled": __debug__}))
else:
    if args.output is not None:
        parser.error("verify does not accept --output")
    with tempfile.TemporaryDirectory(prefix="six-clock-verify-") as folder:
        temp = Path(folder)
        recompute(temp)
        for name in ("local-algebra.json", "independent-certificate.json"):
            actual = json.loads((temp / name).read_text())
            expected = json.loads((HERE / name).read_text())
            if actual != expected:
                raise RuntimeError("Recomputed certificate differs: " + name)
    check_manifest()
    print(json.dumps({"ok": True, "command": "verify", "certificates_compared": 2,
                      "archived_evidence_modified": False, "assertions_enabled": __debug__}))
