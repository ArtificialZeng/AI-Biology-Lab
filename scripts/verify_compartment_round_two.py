#!/usr/bin/env python3
"""Recompute round-two exact certificates while preserving archived evidence."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
ENV = {k: v for k, v in os.environ.items()
       if k in {"PATH", "HOME", "LANG", "LC_ALL", "TMPDIR", "TEMP", "TMP", "SYSTEMROOT"}}


def verify_manifest(project):
    manifest = json.loads((project / "archive_manifest.json").read_text())
    for filename, expected in manifest["sha256"].items():
        actual = hashlib.sha256((project / filename).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError("Archived mathematical input changed: " + str(project / filename))


def replay(project_name, script, receipt, extra=()):
    project = ROOT / "projects" / project_name
    expected = json.loads((project / "evidence" / receipt).read_text())
    for optimized in (False, True):
        args = [sys.executable] + (["-O"] if optimized else []) + [str(project / script), *extra]
        run = subprocess.run(args, cwd=ROOT, env=ENV, capture_output=True, text=True, check=True)
        if json.loads(run.stdout) != expected:
            raise ValueError("Exact receipt differs from archive: " + project_name + "/" + script)
    print("PASS ordinary/optimized; archive unchanged:", project_name + "/" + script)


def main():
    for project in ("cycle-global-fibers", "tree-two-leak-sensors"):
        verify_manifest(ROOT / "projects" / project)
    replay("cycle-global-fibers", "verify.py", "certificate.json", ("--self-test",))
    replay("cycle-global-fibers", "independent_verify.py", "independent_certificate.json")
    replay("tree-two-leak-sensors", "verify.py", "certificate.json")
    replay("tree-two-leak-sensors", "independent_verify.py", "independent_certificate.json")


if __name__ == "__main__":
    main()
