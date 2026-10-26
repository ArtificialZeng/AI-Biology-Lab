"""Run unchanged scientific scripts against a real, externally overlaid package."""
import hashlib
import pathlib
import runpy
import sys

if sys.flags.optimize != 0 or not sys.dont_write_bytecode:
    raise RuntimeError("Scientific workers require normal assertions and -B")
overlay, expected, mode, target, *extra = sys.argv[1:]
sys.path.insert(0, overlay)
import cogent3
from cogent3.evolve import fast_distance
source = pathlib.Path(fast_distance.__file__).resolve()
if cogent3.__version__ != "2026.10.22" or not source.is_relative_to(pathlib.Path(overlay).resolve()):
    raise RuntimeError("Wrong actual package identity")
if hashlib.sha256(source.read_bytes()).hexdigest() != expected:
    raise RuntimeError("Wrong actual loaded source")
if mode == "script":
    sys.argv = [target]
    runpy.run_path(target, run_name="__main__")
elif mode == "pytest":
    import pytest
    raise SystemExit(pytest.main([target, "-q", "-c", extra[0], "-p", "no:cacheprovider", "--basetemp", extra[1], "--junitxml", extra[2]]))
else:
    raise RuntimeError("Unknown scientific worker mode")
