#!/usr/bin/env python3
"""Archive I/O wrapper for the byte-preserved generic root-count algorithm."""
from pathlib import Path
import importlib.util
import sys
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from archive_tools import run, require


def build_certificate():
    path = Path(__file__).resolve().parent / 'verify_root_count.py'
    spec = importlib.util.spec_from_file_location('_frozen_root_count_algorithm', path)
    require(spec is not None and spec.loader is not None, 'Cannot load root-count algorithm')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.build_certificate()


if __name__ == '__main__':
    raise SystemExit(run(build_certificate, 'root-count/evidence/exact-root-certificate.json', 'generic root-count identities'))
