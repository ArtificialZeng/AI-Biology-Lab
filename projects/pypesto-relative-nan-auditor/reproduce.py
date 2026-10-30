"""Verify the archive, run scoped utility checks, or replay real native tests."""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
HEAD = "053079aa71710a62c1db17f3d0342e96499a3768"
REMOTE = "https://github.com/ICB-DCM/pyPESTO.git"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_archive():
    manifest = json.loads((ROOT / "archive_manifest.json").read_text())
    for name, expected in manifest["files"].items():
        assert digest(ROOT / name) == expected, name
    assert manifest["head"] == HEAD
    assert digest(ROOT / "upstream.patch") == "6ca3c00f863e9d6b9cf57b6c88850b64e35b5381e337890914a467052a81a920"
    print(f"Archive verified: {len(manifest['files'])} bound files.")


def gaussian_oracle(kind):
    """Closed-form observed Gaussian profile, independent of pyPESTO code."""
    x0, y0, rate = {"finite": (1., 1., 1.), "zero": (0., 1., 0.),
                     "unmeasured": (0., 0., 1.), "invalid": (0., 0., 1.)}[kind]
    observations = [(0., .6), (1., .4), (2., .7)]
    if kind == "unmeasured":
        observations = observations[1:]
    residuals, derivatives = [], []
    for t, datum in observations:
        denominator = x0 + y0 + (rate + 1) * t
        if denominator == 0:
            raise ValueError("Undefined output at a measured point")
        prediction = (x0 + rate * t) / denominator
        residuals.append(prediction - datum)
        derivatives.append(t * (y0 + t) / denominator**2)
    variance = math.fsum(value**2 for value in residuals) / len(residuals)
    nllh = len(residuals) / 2 * (1 + math.log(2 * math.pi * variance))
    nllh += .5 * math.log(2 * math.pi) + .5 * .1**2
    gradient = math.fsum(a*b for a, b in zip(residuals, derivatives, strict=True)) / variance
    fim = math.fsum(value**2 for value in derivatives) / variance
    return {"fval": nllh, "gradient": gradient, "conditional_fim": fim}


def utilities(side):
    """Extract two exact archived function ASTs; this is not a package import."""
    import numpy as np
    source = ROOT / "archive" / side / "util.py"
    tree = ast.parse(source.read_text())
    names = {"compute_optimal_sigma", "compute_nllh"}
    nodes = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in names]
    assert {node.name for node in nodes} == names
    namespace = {"np": np}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(source), "exec"), namespace)
    return namespace


def light():
    """Check exact utility source plus an independently evaluated Gaussian oracle."""
    import numpy as np
    verify_archive()
    for kind in ["unmeasured", "finite", "zero"]:
        expected = gaussian_oracle(kind)
        for side in ["baseline", "candidate"]:
            rows = json.loads((ROOT / "archive" / f"probe-output-{side}.json").read_text())
            row = next(r for r in rows if r["kind"] == kind and r["orders"] == [0, 1, 2])
            for actual, oracle in [(row["fval"], expected["fval"]),
                                   (row["gradient"][0], expected["gradient"]),
                                   (row["hessian"][0][0], expected["conditional_fim"])]:
                assert math.isclose(actual, oracle, rel_tol=1e-6, abs_tol=1e-8)
    try:
        gaussian_oracle("invalid")
    except ValueError:
        pass
    else:
        raise AssertionError("The observed undefined ratio cannot have a finite likelihood")
    for side in ["baseline", "candidate"]:
        functions = utilities(side)
        data = [np.array([[.6], [.4], [.7]])]
        simulated = [np.array([[np.nan], [.5], [.5]])]
        full = [np.ones((3, 1), dtype=bool)]
        sigma = functions["compute_optimal_sigma"](data, simulated, full)
        nllh = functions["compute_nllh"](data, simulated, [np.full((3, 1), sigma)], full)
        assert bool(np.isfinite(nllh)) == (side == "baseline")
        missing = [np.array([[False], [True], [True]])]
        sigma = functions["compute_optimal_sigma"](data, simulated, missing)
        nllh = functions["compute_nllh"](data, simulated, [np.full((3, 1), sigma)], missing)
        expected = gaussian_oracle("unmeasured")["fval"] - (.5*math.log(2*math.pi)+.005)
        assert math.isclose(nllh, expected, rel_tol=1e-7, abs_tol=1e-9)
    print("Utility/analytic checks passed; archived native observations matched by tolerance.")
    print("This command does not execute AMICI or the public collector.")


def run(command, **kwargs):
    subprocess.run(command, check=True, **kwargs)


def bootstrap(work, local_source):
    verify_archive()
    work = work.resolve()
    assert not work.is_relative_to(ROOT), "Keep fetched full sources outside the public project"
    work.mkdir(parents=True, exist_ok=True)
    baseline, candidate = work / "baseline", work / "candidate"
    assert not baseline.exists() and not candidate.exists(), "Use a new work directory"
    if local_source:
        run(["git", "clone", "--no-hardlinks", str(local_source.resolve()), str(baseline)])
        run(["git", "checkout", "--detach", HEAD], cwd=baseline)
    else:
        run(["git", "init", str(baseline)])
        run(["git", "remote", "add", "origin", REMOTE], cwd=baseline)
        run(["git", "fetch", "--depth=1", "origin", HEAD], cwd=baseline)
        run(["git", "checkout", "--detach", "FETCH_HEAD"], cwd=baseline)
    assert subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=baseline, text=True).strip() == HEAD
    run(["git", "clone", "--no-hardlinks", str(baseline), str(candidate)])
    run(["git", "apply", "--check", str(ROOT / "upstream.patch")], cwd=candidate)
    run(["git", "apply", str(ROOT / "upstream.patch")], cwd=candidate)
    for side in ["baseline", "candidate"]:
        for name in ["util.py", "calculator.py"]:
            assert digest(work/side/"pypesto/hierarchical/relative"/name) == digest(ROOT/"archive"/side/name)
    assert digest(candidate/"test/hierarchical/test_relative_nonfinite.py") == digest(ROOT/"test_native_regression.py")
    print("Pinned source bootstrap and exact patch application verified.")


def native_one(work, selection):
    for side in ["baseline", "candidate"]:
        actual_head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=work/side, text=True).strip()
        assert actual_head == HEAD
        for name in ["util.py", "calculator.py"]:
            assert digest(work/side/"pypesto/hierarchical/relative"/name) == digest(ROOT/"archive"/side/name)
    assert digest(work/"candidate/test/hierarchical/test_relative_nonfinite.py") == digest(ROOT/"test_native_regression.py")
    sys.path.insert(0, str(work / selection))
    import pypesto
    import pytest
    assert Path(pypesto.__file__).is_relative_to(work / selection)
    counts = {"passed": 0, "failed": 0, "skipped": 0, "errors": 0}
    class CountReports:
        def pytest_runtest_logreport(self, report):
            if report.skipped:
                counts["skipped"] += 1
            elif report.when == "call":
                counts[report.outcome] += 1
            elif report.failed:
                counts["errors"] += 1
    config = work / "replay-pytest.ini"
    config.write_text("[pytest]\n")
    code = pytest.main(["-c", str(config), "-p", "no:cacheprovider", "-q",
                        str(ROOT/"test_native_regression.py")], plugins=[CountReports()])
    result = {"selection": selection, "head": actual_head, "pytest_exit_code": int(code), "counts": counts}
    (work/f"native-{selection}.json").write_text(json.dumps(result, indent=2)+"\n")
    expected = {"passed": 33, "failed": 6, "skipped": 0, "errors": 0} if selection == "baseline" else {"passed": 39, "failed": 0, "skipped": 0, "errors": 0}
    assert counts == expected, result
    assert code == (1 if selection == "baseline" else 0), result
    print(json.dumps(result))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["verify", "light", "bootstrap", "native", "_native_one"])
    parser.add_argument("--work", type=Path)
    parser.add_argument("--from-local", type=Path, help="Optional local Git source; commit/content are still checked")
    parser.add_argument("--selection", choices=["baseline", "candidate"])
    args = parser.parse_args()
    if args.command == "verify": verify_archive()
    elif args.command == "light": light()
    elif args.command == "bootstrap": bootstrap(args.work, args.from_local)
    elif args.command == "native":
        verify_archive()
        for selection in ["baseline", "candidate"]:
            run([sys.executable, str(Path(__file__).resolve()), "_native_one", "--work", str(args.work.resolve()), "--selection", selection])
    else: native_one(args.work.resolve(), args.selection)


if __name__ == "__main__":
    main()
