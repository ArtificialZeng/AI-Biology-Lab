#!/usr/bin/env python3
"""Check a pinned Optlang tree and its local GLPK wrapper patch in temporary copies."""
import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
COMMIT = "6427c6c10220ebff56b1e72591bc4ee43549d81b"
BASE_SHA = "7932237d98ee07c344cc6cf12d067aa8b0496132a46d83fce5c45893e500b8e8"
DEPENDENCIES = {
    "swiglpk": "5.0.13", "sympy": "1.14.0", "mpmath": "1.3.0",
    "pytest": "9.1.1", "iniconfig": "2.3.1", "packaging": "26.3",
    "pluggy": "1.7.0", "Pygments": "2.21.0",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def outside(path):
    resolved = path.resolve()
    require(resolved != HERE and HERE not in resolved.parents, "Keep outputs outside the archived project")
    return resolved


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkout", type=Path, required=True, help="local Optlang clone containing the pinned commit")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--logs", type=Path)
    args = parser.parse_args()
    require(sys.version_info >= (3, 11), "Use Python 3.11 or newer")
    checkout = args.checkout.resolve()
    require(checkout.is_dir(), "Missing local Optlang checkout")
    output = outside(args.output) if args.output else None
    logs = outside(args.logs) if args.logs else None
    if output:
        require(not output.exists(), "Refuse to overwrite an existing receipt")
    if logs:
        require(not logs.exists(), "Refuse to overwrite existing logs")
        logs.mkdir(parents=True)
    manifest = json.loads((HERE / "archive_manifest.json").read_text())
    for name, digest in manifest["sha256"].items():
        require(sha(HERE / name) == digest, "Archived input changed: " + name)
    result = {"optlang_commit": COMMIT, "baseline_glpk_interface_sha256": BASE_SHA,
              "wrapper_patch_sha256": sha(HERE / "wrapper/optlang-mip-validation.patch"),
              "test_sha256": sha(HERE / "wrapper/test_glpk_mip_validation.py"),
              "source_method": "git archive of pinned commit; worktree changes ignored",
              "python": ".".join(map(str, sys.version_info[:2])),
              "scope": "GLPK-only primal residual validation and one disabled-presolve recovery per _optimize",
              "native_library_modified": False, "checks": {}}
    with tempfile.TemporaryDirectory(prefix="optlang-glpk-wrapper-audit-") as folder:
        temp = Path(folder)
        (temp / "home").mkdir()
        (temp / "tmp").mkdir()

        def environment(tree):
            return {"PATH": os.defpath, "HOME": str(temp / "home"), "TMPDIR": str(temp / "tmp"),
                    "LANG": "C", "PYTHONDONTWRITEBYTECODE": "1", "PYTHONPATH": str(tree / "src"),
                    "PYTEST_DISABLE_PLUGIN_AUTOLOAD": "1"}

        result["environment_allowlist"] = sorted(environment(temp))

        def run(argv, tree, label, expected=0, binary=False):
            process = subprocess.run([str(x) for x in argv], cwd=tree, env=environment(tree),
                                     text=not binary, capture_output=True, timeout=120)
            if not binary and logs:
                clean = process.stdout + process.stderr
                for path, replacement in ((temp, "<temporary>"), (HERE, "<project>"), (checkout, "<checkout>")):
                    clean = clean.replace(str(path), replacement)
                (logs / (label + ".txt")).write_text(clean)
            require(process.returncode == expected, "Unexpected command result: " + label + "; inspect separate logs")
            return process.stdout

        archive = run(["git", "-C", checkout, "archive", "--format=tar", COMMIT], temp, "source-archive", binary=True)
        source = temp / "source"
        source.mkdir()
        with tarfile.open(fileobj=io.BytesIO(archive), mode="r:") as stream:
            for member in stream.getmembers():
                name = Path(member.name)
                require(not name.is_absolute() and ".." not in name.parts, "Unsafe source member")
                require(member.isfile() or member.isdir(), "Unexpected linked/special source member")
            filtering = {"filter": "data"} if hasattr(tarfile, "data_filter") else {}
            stream.extractall(source, **filtering)
        require(sha(source / "src/optlang/glpk_interface.py") == BASE_SHA, "Pinned source bytes differ")
        test = temp / "test_glpk_mip_validation.py"
        shutil.copyfile(HERE / "wrapper/test_glpk_mip_validation.py", test)
        import_code = """
import importlib.metadata as metadata
import json
import sys
from pathlib import Path
from optlang import glpk_interface as interface
import swiglpk
if Path(sys.argv[1]).resolve() not in Path(interface.__file__).resolve().parents:
    raise RuntimeError("Imported the wrong source tree")
expected = json.loads(sys.argv[2])
actual = {name: metadata.version(name) for name in expected}
if actual != expected or swiglpk.glp_version() != "5.0":
    raise RuntimeError("Unexpected dependency versions")
print(json.dumps({"versions": actual, "glpk": swiglpk.glp_version(), "source_import_verified": True}))
"""
        for label in ("baseline", "candidate"):
            tree = temp / label
            shutil.copytree(source, tree)
            if label == "candidate":
                patch = HERE / "wrapper/optlang-mip-validation.patch"
                run(["git", "apply", "--check", patch], tree, "patch-check")
                run(["git", "apply", patch], tree, "patch-apply")
            imported = json.loads(run([sys.executable, "-c", import_code, tree, json.dumps(DEPENDENCIES)], tree,
                                      label + "-import-check"))
            report_path = temp / (label + "-junit.xml")
            argv = [sys.executable, "-m", "pytest", "--tb=short", "-q", "-o", "addopts=",
                    "-o", "cache_dir=" + str(temp / (label + "-cache")), "--junitxml=" + str(report_path), test]
            if label == "baseline":
                argv.extend(["-k", "not test_primal_check_uses_absolute_or_relative_configured_tolerance"])
            stdout = run(argv, tree, label + "-pytest", expected=1 if label == "baseline" else 0)
            suite = ET.parse(report_path).getroot()
            summaries = list(suite.iter("testsuite"))
            counts = {key: sum(int(s.get(key, 0)) for s in summaries)
                      for key in ("tests", "failures", "errors", "skipped")}
            require(counts["errors"] == counts["skipped"] == 0, "Unexpected test errors or skips")
            counts["passed"] = counts["tests"] - counts["failures"]
            counts["deselected"] = 8 if label == "baseline" else 0
            require(counts["tests"] == (54 if label == "baseline" else 62), "Coverage changed")
            require(counts["failures"] == (27 if label == "baseline" else 0), "Regression outcome changed")
            if label == "baseline":
                require("8 deselected" in stdout, "Baseline helper exclusion changed")
            result["checks"][label] = {"counts": counts, "import": imported}
    rendered = json.dumps(result, indent=2) + "\n"
    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
