"""Read-only deterministic archive replay; Python standard library only."""
from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
JOBS = (
    ("finite_states", "verify_states.py", "evidence/states.json"),
    ("constructive_designs", "verify_designs.py", "evidence/designs.json"),
    ("independent_oracle", "independent_verify.py", "evidence/independent.json"),
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")


def main() -> dict:
    manifest = json.loads((HERE / "archive_manifest.json").read_text(encoding="utf-8"))
    require(manifest.get("schema") == "rna-short-helix-archive-v1", "Unknown archive schema")
    bindings = manifest.get("sha256")
    require(isinstance(bindings, dict) and bool(bindings), "Missing archive file bindings")
    for name, expected in sorted(bindings.items()):
        path = (HERE / name).resolve()
        require(path.is_relative_to(HERE), f"Archive path leaves the project: {name}")
        require(path.is_file(), f"Missing archived file: {name}")
        require(digest(path.read_bytes()) == expected, f"Archive hash mismatch: {name}")
    result = []
    for label, program, archive in JOBS:
        require(program in bindings and archive in bindings, f"Unbound replay input: {label}")
        syntax = ast.parse((HERE / program).read_text(encoding="utf-8"))
        require(not any(isinstance(node, ast.Assert) for node in ast.walk(syntax)),
                f"Checks removable by Python -O are forbidden: {program}")
        expected = json.loads((HERE / archive).read_text(encoding="utf-8"))
        outputs = []
        for optimized in (False, True):
            command = [sys.executable, "-I", "-B"]
            if optimized:
                command.append("-O")
            command.append(str(HERE / program))
            run = subprocess.run(command, cwd=HERE, capture_output=True, check=True, timeout=120)
            observed = json.loads(run.stdout)
            require(observed == expected, f"Result differs from archive: {program}, optimized={optimized}")
            outputs.append(run.stdout)
        require(outputs[0] == outputs[1], f"Ordinary/-O stdout differs: {program}")
        result.append({"check": label, "archive": archive,
                       "canonical_result_sha256": digest(canonical(expected)),
                       "stdout_sha256": digest(outputs[0]),
                       "ordinary_matches_archive": True, "optimized_matches_archive": True,
                       "mode_stdout_identical": True})
    return {"schema": "rna-archive-replay-v1", "archive_files_verified": len(bindings),
            "manifest_sha256": digest((HERE / "archive_manifest.json").read_bytes()),
            "results": result}


if __name__ == "__main__":
    try:
        print(json.dumps(main(), indent=2))
    except (ValueError, OSError, subprocess.SubprocessError) as error:
        print(f"RNA archive replay failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error
