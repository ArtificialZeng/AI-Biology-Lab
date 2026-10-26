"""Verify exact witnesses and real API behavior in disposable package copies.

Archived evidence is read-only. Default output is JSON on stdout; --output may
write a new report outside evidence/. Installed packages are never edited.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE / "scripts"))
from exact_oracle import exact_oracle


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def clean_environment(root, source):
    home = root / "home"
    home.mkdir(exist_ok=True)
    return {
        "PATH": "/usr/bin:/bin",
        "HOME": str(home),
        "TMPDIR": str(root),
        "LANG": "C",
        "PYTHONDONTWRITEBYTECODE": "1",
        "PYTHONPATH": str(source),
        "PYTEST_DISABLE_PLUGIN_AUTOLOAD": "1",
    }


def sanitize(text):
    text = re.sub(r"/(?:Users|home)/[^\s\"']+", "<local-path>", text)
    return text[-4000:]


def verify_manifest():
    manifest = BASE / "archive_manifest.json"
    if not manifest.exists():
        raise RuntimeError("Missing archived-evidence manifest")
    for entry in json.loads(manifest.read_text())["files"]:
        path = BASE / entry["path"]
        if digest(path) != entry["sha256"]:
            raise RuntimeError("Archived evidence mismatch: " + entry["path"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--child", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.child:
        from probe import probe

        print(json.dumps(probe(args.child), sort_keys=True))
        return
    if sys.version_info[:2] != (3, 11):
        raise RuntimeError("This pinned audit requires Python 3.11")
    verify_manifest()
    source_manifest = json.loads((BASE / "source-manifest.json").read_text())
    spec = importlib.util.find_spec("memote")
    if spec is None or spec.origin is None:
        raise RuntimeError("Install the pinned requirements first")
    installed = Path(spec.origin).parent
    for relative, expected in source_manifest["files"].items():
        if digest(installed / relative) != expected:
            raise RuntimeError("Installed source differs from pin: " + relative)
    report = {
        "source_commit": source_manifest["commit"],
        "exact_oracle": exact_oracle(),
        "source_hashes_verified": True,
    }
    with tempfile.TemporaryDirectory(prefix="memote-cache-audit-") as directory:
        temporary = Path(directory)
        original = temporary / "original" / "src" / "memote"
        candidate = temporary / "candidate" / "src" / "memote"
        shutil.copytree(
            installed, original, ignore=shutil.ignore_patterns("__pycache__", "*.pyc")
        )
        shutil.copytree(original, candidate)
        candidate_root = candidate.parents[1]
        patch = BASE / "patches" / "clear-caches-between-api-runs.patch"
        patched = subprocess.run(
            ["/usr/bin/git", "apply", "--check", str(patch)],
            cwd=candidate_root,
            env=clean_environment(temporary, candidate.parent),
            capture_output=True,
            text=True,
        )
        if patched.returncode:
            raise RuntimeError(
                "Production patch does not apply: " + sanitize(patched.stderr)
            )
        subprocess.run(
            ["/usr/bin/git", "apply", str(patch)],
            cwd=candidate_root,
            env=clean_environment(temporary, candidate.parent),
            check=True,
            capture_output=True,
            text=True,
        )
        for label, package in (("baseline", original), ("candidate", candidate)):
            env = clean_environment(temporary, package.parent)
            child = subprocess.run(
                [
                    sys.executable,
                    "-B",
                    str(Path(__file__).resolve()),
                    "--child",
                    str(package),
                ],
                cwd=temporary,
                env=env,
                capture_output=True,
                text=True,
            )
            if child.returncode:
                raise RuntimeError(
                    label + " subprocess failed: " + sanitize(child.stderr)
                )
            report[label] = json.loads(child.stdout)
            report[label]["environment_allowlist"] = sorted(env)
    expected_failures = {
        "add_refresh",
        "remove_refresh",
        "annotation_refresh",
        "collection_refresh",
        "exception_retry_refresh",
    }
    failures = {
        name for name, passed in report["baseline"]["checks"].items() if not passed
    }
    if failures != expected_failures:
        raise RuntimeError("Baseline did not reproduce the expected stale results")
    if not all(report["candidate"]["checks"].values()):
        raise RuntimeError("Candidate regression failed")
    if report["exact_oracle"] != json.loads(
        (BASE / "evidence" / "exact-oracle.json").read_text()
    ):
        raise RuntimeError("Exact oracle differs from archived evidence")
    report["baseline_expected_failures"] = sorted(failures)
    report["result"] = "PASS"
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        output = args.output.resolve()
        if output == BASE or output.is_relative_to(BASE):
            raise RuntimeError(
                "Archived project is read-only; choose a different output path"
            )
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("x", encoding="utf-8") as stream:
            stream.write(text)
    print(text, end="")


if __name__ == "__main__":
    main()
