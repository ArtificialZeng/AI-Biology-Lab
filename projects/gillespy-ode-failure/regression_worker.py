"""Execute preserved upstream tests against one selected actual library tree."""
import importlib.util
import io
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parent
tree, scope = sys.argv[1:]
source = ROOT / 'vendor' / tree
sys.path[:0] = [str(source), str(source / 'test')]
spec = importlib.util.spec_from_file_location(
    'candidate_ode_regressions', ROOT / 'vendor/candidate/test/test_ode_solver.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
suite = (unittest.defaultTestLoader.loadTestsFromModule(module) if scope == 'whole'
         else unittest.defaultTestLoader.loadTestsFromTestCase(module.TestODEIntegrationFailure))
stream = io.StringIO()
result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
print(json.dumps({'worker_result': {'tests': result.testsRun,
                                  'failures': len(result.failures),
                                  'errors': len(result.errors),
                                  'skipped': len(result.skipped)}}))
if result.errors:
    print(stream.getvalue(), file=sys.stderr)
