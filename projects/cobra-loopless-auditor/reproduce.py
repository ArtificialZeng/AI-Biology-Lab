#!/usr/bin/env python3
"""Reproduce the defect and repair in temporary copies of a local source checkout."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile


HERE = Path(__file__).resolve().parent
PIN = "a9c5ced29670bed5035fb19e87e7b0b588adb62d"
ERROR = "IndexError: arrays used as indices must be of integer (or boolean) type"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def run_tests(source, temporary):
    env = {"PATH": "/usr/bin:/bin", "HOME": str(Path.home()), "LANG": "en_US.UTF-8",
           "PYTHONPATH": str(source / "src"), "PYTEST_DISABLE_PLUGIN_AUTOLOAD": "1",
           "SKIP_MP": "1"}
    # Prevent an installed package from masking the source under examination.
    script = (
        "import cobra,pytest; from pathlib import Path; "
        "cobra.Configuration().processes=1; "
        "actual=Path(cobra.__file__).resolve(); "
        f"expected=Path({str(source / 'src')!r}).resolve(); "
        "assert actual.is_relative_to(expected), (actual,expected); "
        f"raise SystemExit(pytest.main({[str(HERE / 'test_independent.py'), '-q', '--tb=short', '-p', 'no:cacheprovider']!r}))"
    )
    result = subprocess.run([sys.executable, "-c", script], cwd=temporary, env=env,
                            capture_output=True, text=True)
    log = result.stdout + result.stderr
    count = lambda word: int(re.findall(r"(\d+) " + word, log)[-1]) if re.findall(r"(\d+) " + word, log) else 0
    return result.returncode, {"failed": count("failed"), "passed": count("passed"),
                               "indexing_exceptions": log.count(ERROR)}, log


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkout", type=Path, required=True,
                        help="unmodified local source at the documented upstream commit")
    parser.add_argument("--output", type=Path, help="write a separate receipt; default stdout")
    parser.add_argument("--logs", type=Path, help="optional separate directory for test logs")
    args = parser.parse_args()
    original = args.checkout.resolve()
    expected = json.loads((HERE / "source_hashes.json").read_text())
    for filename, sha in expected.items():
        require(digest(original / filename) == sha, "Unexpected source bytes: " + filename)
    with tempfile.TemporaryDirectory(prefix="cobra-loopless-") as directory:
        temporary = Path(directory)
        candidate = temporary / "candidate"
        shutil.copytree(original, candidate,
                        ignore=shutil.ignore_patterns(".git", "__pycache__", ".pytest_cache", ".venv"))
        subprocess.run(["git", "apply", "--check", str(HERE / "upstream.patch")],
                       cwd=candidate, check=True, capture_output=True, text=True)
        subprocess.run(["git", "apply", str(HERE / "upstream.patch")],
                       cwd=candidate, check=True, capture_output=True, text=True)
        red_code, red, red_log = run_tests(original, temporary)
        green_code, green, green_log = run_tests(candidate, temporary)
        if args.logs:
            args.logs.mkdir(parents=True, exist_ok=True)
            (args.logs / "baseline.log").write_text(red_log)
            (args.logs / "candidate.log").write_text(green_log)
        require(red_code == 1 and red == {"failed": 24, "passed": 23, "indexing_exceptions": 24},
                "Baseline did not reproduce the expected defect; inspect separate logs")
        require(green_code == 0 and green == {"failed": 0, "passed": 47, "indexing_exceptions": 0},
                "Candidate failed the independent contract checks; inspect separate logs")
        result = {"pinned_upstream": PIN, "source_sha256": expected,
                  "patch_sha256": digest(HERE / "upstream.patch"),
                  "independent_tests_sha256": digest(HERE / "test_independent.py"),
                  "baseline": red, "candidate": green,
                  "arithmetic_and_solver": "analytic flux oracles, GLPK numerical LP/MILP",
                  "scope": "empty internal selections, cyclic finder, and downstream FVA; GLPK only",
                  "limitation": "Software regression evidence; not a new biological mechanism or proof about all metabolic models."}
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
