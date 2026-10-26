#!/usr/bin/env python3
"""Reproduce signed FASTCC checks using read-only input and temporary source copies."""
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
PIN = "a9c5ced29670bed5035fb19e87e7b0b588adb62d"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def outside_project(path):
    resolved = path.resolve()
    require(resolved != HERE and HERE not in resolved.parents,
            "Write run outputs outside the project to preserve archived evidence")
    return resolved


def sanitize(text, temporary, checkout):
    for path, label in ((temporary, "<temporary>"), (checkout, "<checkout>"),
                        (HERE, "<project>"), (Path(sys.prefix), "<environment>")):
        text = text.replace(str(path), label)
    return re.sub(r"(?:\.\./)+(?=<)", "", text)


def verify_archive():
    manifest = json.loads((HERE / "archive_manifest.json").read_text())
    for relative, sha in manifest["sha256"].items():
        require(digest(HERE / relative) == sha, "Archived file changed: " + relative)


def run_tests(source, temporary, checkout):
    # No inherited credentials or plugin configuration; only these seven keys.
    env = {"PATH": os.defpath, "HOME": str(temporary / "home"), "LANG": "C.UTF-8",
           "TMPDIR": str(temporary / "tmp"), "PYTHONDONTWRITEBYTECODE": "1",
           "PYTEST_DISABLE_PLUGIN_AUTOLOAD": "1", "PYTHONPATH": str(source / "src")}
    script = '''
import importlib
from pathlib import Path
import cobra
import pytest
import sys
source = Path(sys.argv[1]).resolve()
module = importlib.import_module("cobra.flux_analysis.fastcc")
for actual, relative in ((cobra.__file__, "cobra/__init__.py"),
                         (module.__file__, "cobra/flux_analysis/fastcc.py")):
    if Path(actual).resolve() != source / "src" / relative:
        raise RuntimeError("Imported an installed package instead of the requested source")
cobra.Configuration().processes = 1
raise SystemExit(pytest.main([sys.argv[2], "--tb=short", "-q", "--noconftest",
                             "-c", sys.argv[3], "-p", "no:cacheprovider"]))
'''
    result = subprocess.run([sys.executable, "-c", script, str(source),
                             str(HERE / "test_independent.py"), str(temporary / "pytest.ini")],
                            cwd=temporary, env=env, capture_output=True, text=True)
    log = sanitize(result.stdout + result.stderr, temporary, checkout)
    def count(word):
        matches = re.findall(r"(\d+) " + word + r"\b", log)
        return int(matches[-1]) if matches else 0
    return {"returncode": result.returncode, "failed": count("failed"),
            "passed": count("passed"), "warnings": count("warnings?")}, log


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkout", type=Path, required=True,
                        help="unmodified checkout or source archive at the pinned commit")
    parser.add_argument("--output", type=Path, help="write a separate JSON receipt; default stdout")
    parser.add_argument("--logs", type=Path, help="optional separate directory for sanitized pytest logs")
    args = parser.parse_args()
    verify_archive()
    original = args.checkout.resolve()
    require(original.is_dir(), "The supplied checkout directory does not exist")
    output = outside_project(args.output) if args.output else None
    logs = outside_project(args.logs) if args.logs else None
    for destination in (output, logs):
        if destination:
            require(destination != original and original not in destination.parents,
                    "Write run outputs outside the input checkout")
    expected = json.loads((HERE / "source_hashes.json").read_text())
    for relative, sha in expected.items():
        require((original / relative).is_file(), "Missing frozen source file: " + relative)
        require(digest(original / relative) == sha, "Unexpected source bytes: " + relative)
    expected_python = {name for name in expected if name.endswith(".py")}
    actual_python = {path.relative_to(original).as_posix() for path in (original / "src/cobra").rglob("*.py")}
    require(actual_python == expected_python, "Unexpected Python source inventory")
    if (original / ".git").exists():
        revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=original,
                                  capture_output=True, text=True, check=True).stdout.strip()
        require(revision == PIN, "Git checkout is not at the pinned upstream commit")
        clean = subprocess.run(["git", "diff", "--quiet", "HEAD", "--"], cwd=original)
        require(clean.returncode == 0, "Use an unmodified tracked checkout")
    with tempfile.TemporaryDirectory(prefix="fastcc-direction-") as directory:
        temporary = Path(directory)
        (temporary / "tmp").mkdir()
        (temporary / "home").mkdir()
        (temporary / "pytest.ini").write_text("[pytest]\n")
        for label in ("baseline", "candidate"):
            shutil.copytree(original, temporary / label,
                            ignore=shutil.ignore_patterns(".git", ".venv", "__pycache__", ".pytest_cache", "*.pyc"))
        candidate = temporary / "candidate"
        subprocess.run(["git", "apply", "--check", str(HERE / "upstream.patch")], cwd=candidate,
                       capture_output=True, text=True, check=True)
        subprocess.run(["git", "apply", str(HERE / "upstream.patch")], cwd=candidate,
                       capture_output=True, text=True, check=True)
        baseline, baseline_log = run_tests(temporary / "baseline", temporary, original)
        repaired, repaired_log = run_tests(candidate, temporary, original)
        candidate_source_sha256 = digest(candidate / "src/cobra/flux_analysis/fastcc.py")
        if logs:
            logs.mkdir(parents=True, exist_ok=True)
            (logs / "baseline.txt").write_text(baseline_log)
            (logs / "candidate.txt").write_text(repaired_log)
        require(baseline == {"returncode": 1, "failed": 1, "passed": 4, "warnings": 0},
                "Expected baseline 1 failed / 4 passed; inspect separate logs")
        require(repaired == {"returncode": 0, "failed": 0, "passed": 5, "warnings": 0},
                "Expected candidate 5 passed without warnings; inspect separate logs")
    result = {"pinned_upstream": PIN, "source_hashes_sha256": digest(HERE / "source_hashes.json"),
              "baseline_fastcc_sha256": expected["src/cobra/flux_analysis/fastcc.py"],
              "candidate_fastcc_sha256": candidate_source_sha256,
              "patch_sha256": digest(HERE / "upstream.patch"),
              "tests_sha256": digest(HERE / "test_independent.py"),
              "baseline": baseline, "candidate": repaired,
              "source_import_checks": "cobra and FASTCC matched each temporary source copy",
              "environment_allowlist": ["HOME", "LANG", "PATH", "PYTEST_DISABLE_PLUGIN_AUTOLOAD",
                                        "PYTHONDONTWRITEBYTECODE", "PYTHONPATH", "TMPDIR"],
              "versions": {name: importlib.metadata.version(name) for name in
                           ("numpy", "optlang", "swiglpk", "pytest", "sympy")},
              "scope": "Sv=0 and reaction interval bounds; GLPK regression coverage",
              "evidence_level": "Independent AI review and deterministic software checks; no external human peer review"}
    rendered = json.dumps(result, indent=2) + "\n"
    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
