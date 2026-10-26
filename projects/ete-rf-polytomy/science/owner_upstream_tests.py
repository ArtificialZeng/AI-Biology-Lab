from pathlib import Path
import sys,unittest,json,io,importlib
R=Path(__file__).resolve().parent
runtime=Path(sys.argv[1]).resolve();source=R/'candidate'
sys.path.insert(0,str(runtime));sys.path.insert(1,str((R.parent/'compbio-primary-triage/evidence/runtime').resolve()));sys.path.insert(2,str(source))
import ete4
if not Path(ete4.__file__).resolve().is_relative_to(runtime):raise RuntimeError('wrong runtime')
suite=unittest.TestSuite()
suite.addTests(unittest.defaultTestLoader.loadTestsFromName('tests.test_rf_polytomy_correction'))
suite.addTests(unittest.defaultTestLoader.loadTestsFromName('tests.test_tree.Test_Core_Tree.test_robinson_foulds_and_more'))
stream=io.StringIO();r=unittest.TextTestRunner(stream=stream,verbosity=2).run(suite)
data={'runtime':str(runtime),'actual_package':ete4.__file__,'tests_run':r.testsRun,'successful':r.wasSuccessful(),'failure_count':len(r.failures),'error_count':len(r.errors),'output':stream.getvalue()}
(R/'evidence'/('upstream-tests-'+runtime.name+'.json')).write_text(json.dumps(data,indent=2)+'\n');print(json.dumps(data),flush=True)
if not r.wasSuccessful():raise SystemExit(1)
