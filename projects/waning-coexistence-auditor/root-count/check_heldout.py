#!/usr/bin/env python3
"""Archive I/O wrapper for the byte-preserved independent held-out algorithm."""
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import importlib.util
import json
import sys
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from archive_tools import run, require


def build_certificate():
    path = Path(__file__).resolve().parent / 'evidence/heldout_check.py'
    spec = importlib.util.spec_from_file_location('_frozen_heldout_algorithm', path)
    require(spec is not None and spec.loader is not None, 'Cannot load held-out algorithm')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    output = StringIO()
    with redirect_stdout(output):
        module.main()
    value = json.loads(output.getvalue())
    require(value.get('ok') is True, 'Held-out algorithm reported failure')
    return value


if __name__ == '__main__':
    raise SystemExit(run(build_certificate, 'root-count/evidence/heldout-certificate.json', 'independent exact held-out root isolation'))
