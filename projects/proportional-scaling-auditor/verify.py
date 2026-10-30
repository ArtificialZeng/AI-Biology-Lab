"""Read-only archive wrapper for the preserved exact counterexample calculation."""
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import json
import runpy
import sys
sys.dont_write_bytecode = True
import sympy
from archive_tools import run, require, encode


def build_certificate():
    require(sympy.__version__ == '1.13.3', 'Use pinned SymPy1.13.3')
    captured = StringIO()
    with redirect_stdout(captured):
        namespace = runpy.run_path(str(Path(__file__).resolve().parent / 'check_candidate.py'))
    value = namespace['out']
    require(json.loads(captured.getvalue()) == value, 'Calculation output differs from returned evidence')
    return value


if __name__ == '__main__':
    raise SystemExit(run(build_certificate, 'evidence/exact-certificate.json', 'preserved quotient-rule calculation'))
