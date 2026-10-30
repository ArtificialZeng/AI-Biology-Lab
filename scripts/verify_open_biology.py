#!/usr/bin/env python3
"""Recompute both research packages without modifying their archived receipts."""
import json
import os
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
ENV = {key: value for key, value in os.environ.items()
       if key in {"PATH", "HOME", "LANG", "LC_ALL", "TMPDIR", "TEMP", "TMP",
                  "SYSTEMROOT", "VIRTUAL_ENV"}}


def replay(script, archived, self_test=False):
    expected = json.loads((ROOT / archived).read_text())
    for optimized in (False, True):
        command = [sys.executable]
        if optimized:
            command.append("-O")
        command.append(str(ROOT / script))
        if self_test:
            command.append("--self-test")
        result = subprocess.run(command, cwd=ROOT, env=ENV, check=True,
                                capture_output=True, text=True)
        if json.loads(result.stdout) != expected:
            raise ValueError("Recomputed receipt differs from archive: " + script)
    print("PASS ordinary/optimized and unchanged archive:", script)


def main():
    replay("projects/directed-cycle-identifiability/verify.py",
           "projects/directed-cycle-identifiability/evidence/certificate.json", True)
    replay("projects/directed-cycle-identifiability/independent_reconstruction.py",
           "projects/directed-cycle-identifiability/evidence/independent_certificate.json")
    replay("projects/dual-futile-symmetry/verify_symmetry_slice.py",
           "projects/dual-futile-symmetry/evidence/polynomial_certificate.json", True)
    replay("projects/dual-futile-symmetry/derive_symmetry_slice_sympy.py",
           "projects/dual-futile-symmetry/evidence/sympy_certificate.json")
    replay("projects/dual-futile-symmetry/independent_reconstruction.py",
           "projects/dual-futile-symmetry/evidence/independent_certificate.json")


if __name__ == "__main__":
    main()
