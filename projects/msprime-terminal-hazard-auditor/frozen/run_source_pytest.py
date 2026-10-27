import argparse,json,sys,hashlib
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--test-path',type=Path);p.add_argument('pytest_args',nargs='*');a=p.parse_args()
s=a.source.resolve();sys.path.insert(0,str(s))
import msprime,pytest
expected=json.loads((Path(__file__).parent/'expected_sources.json').read_text())[s.name]
def guard():
    assert Path(msprime.__file__).resolve().is_relative_to(s)
    assert Path(msprime._msprime.__file__).resolve().is_relative_to(s)
    assert Path(msprime.demography.__file__).resolve()==s/'msprime/demography.py'
    assert hashlib.sha256((s/'msprime/demography.py').read_bytes()).hexdigest()==expected['demography_sha256']
    assert hashlib.sha256(Path(msprime._msprime.__file__).read_bytes()).hexdigest()==expected['native_sha256']
guard()
print(json.dumps({'module':msprime.__file__,'native':msprime._msprime.__file__,'demography_sha256':hashlib.sha256((s/'msprime/demography.py').read_bytes()).hexdigest()}))
class Guard:
    def pytest_sessionstart(self,session):guard()
    def pytest_runtest_setup(self,item):guard()
    def pytest_sessionfinish(self,session,exitstatus):guard()
test=(a.test_path or s/'tests/test_demography.py').resolve()
result=pytest.main(['-o','addopts=','-p','no:cacheprovider','-q',str(test)]+a.pytest_args,plugins=[Guard()])
guard()
raise SystemExit(result)
