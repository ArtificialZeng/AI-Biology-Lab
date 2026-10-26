"""Read-only replay of the frozen Example 8 exact checks, including Python -O."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replay() -> dict:
    root = Path(__file__).resolve().parent
    manifest_path = root / "archive_manifest.json"
    manifest_digest = digest(manifest_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != "example8-exact-replay-v1":
        raise RuntimeError("Unsupported archive manifest schema")
    files = manifest["files"]
    expected_runs = {
        (script, optimized)
        for script in ("verify.py", "independent_verify.py")
        for optimized in (False, True)
    }
    actual_runs = {(run["script"], run["optimized"]) for run in manifest["runs"]}
    if len(manifest["runs"]) != 4 or actual_runs != expected_runs:
        raise RuntimeError("Manifest must contain exactly the four ordinary/optimized runs")

    def check_hashes() -> None:
        if digest(manifest_path) != manifest_digest:
            raise RuntimeError("Archive manifest changed during replay")
        for relative, expected in files.items():
            path = (root / relative).resolve()
            if not path.is_relative_to(root) or digest(path) != expected:
                raise RuntimeError(f"Archive/input digest mismatch: {relative}")

    check_hashes()
    env = {k: v for k, v in os.environ.items() if k in {"PATH", "HOME", "LANG", "LC_ALL", "TMPDIR", "TEMP", "TMP", "SYSTEMROOT"}}
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    checked = []
    for run in manifest["runs"]:
        arguments = [sys.executable, "-B"]
        if run["optimized"]:
            arguments.append("-O")
        arguments.append(str(root / run["script"]))
        completed = subprocess.run(arguments, env=env, capture_output=True, check=False)
        if completed.returncode != 0 or completed.stderr:
            raise RuntimeError(f"Checker failed: {run['script']} optimized={run['optimized']}")
        archived = (root / run["archive"]).read_bytes()
        if completed.stdout != archived:
            raise RuntimeError(f"Replay output differs: {run['archive']}")
        report = json.loads(completed.stdout)
        if report.get("status") != "passed" or report.get("sympy_version") != "1.13.3":
            raise RuntimeError("Unexpected checker status or SymPy version")
        checked.append(run["archive"])
    check_hashes()
    return {
        "status": "passed",
        "runs": checked,
        "bound_files": len(files),
        "archive_unchanged": True,
        "writes_requested": False,
        "scope": "Exact algebra replay; analytic proof and scope are documented separately",
    }


if __name__ == "__main__":
    try:
        sys.stdout.write(json.dumps(replay(), indent=2) + "\n")
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        sys.stderr.write(f"replay failed: {error}\n")
        raise SystemExit(1)
