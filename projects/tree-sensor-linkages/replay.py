#!/usr/bin/env python3
"""Replay exact tree-linkage checks without modifying the archived evidence."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

PROJECT = Path(__file__).resolve().parent
ROOT = PROJECT
ENV = {k: v for k, v in os.environ.items() if k in {
    'PATH', 'HOME', 'LANG', 'LC_ALL', 'TMPDIR', 'TEMP', 'TMP', 'SYSTEMROOT'}}


def main():
    manifest = json.loads((PROJECT / 'archive_manifest.json').read_text())
    for filename, expected in manifest['sha256'].items():
        actual = hashlib.sha256((PROJECT / filename).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError('Archived mathematical input changed: ' + filename)
    for script, receipt in [('verify.py', 'certificate.json'),
                            ('independent_verify.py', 'independent_certificate.json')]:
        expected = json.loads((PROJECT / 'evidence' / receipt).read_text())
        for optimized in (False, True):
            args = [sys.executable] + (['-O'] if optimized else []) + [str(PROJECT / script)]
            result = subprocess.run(args, cwd=ROOT, env=ENV, capture_output=True,
                                    text=True, check=True)
            if json.loads(result.stdout) != expected:
                raise ValueError('Recomputed certificate differs: ' + script)
        print('PASS ordinary/optimized; archive unchanged:', script)


if __name__ == '__main__':
    main()
