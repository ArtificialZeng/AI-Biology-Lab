#!/usr/bin/env python3
"""Build pinned baseline/candidate tskit and replay one exact reference audit."""
from __future__ import annotations

import argparse
import ast
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.request

ROOT = Path(__file__).resolve().parent
COMMIT = "f154d1190e0142c0499d2698d76eed0fe6429876"
URL = f"https://codeload.github.com/tskit-dev/tskit/tar.gz/{COMMIT}"
BASELINE_C_SHA256 = "2c5e9f77ad5901d6c05801abf7d46d2e8782187fe68551b7556138d27f21a616"
CASE_COUNT = 396


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def verify_packet():
    manifest = json.loads((ROOT / "manifest.json").read_text())
    for name, expected in manifest["sha256"].items():
        path = ROOT / name
        if not path.is_file() or digest(path) != expected:
            raise RuntimeError(f"Packet integrity mismatch: {name}")
    for name in ("audit.py", "reproduce.py"):
        ast.parse((ROOT / name).read_text(), filename=name)
    patch = (ROOT / "native.patch").read_text()
    changed = re.findall(r"^diff --git a/(\S+) b/(\S+)$", patch, re.MULTILINE)
    wanted = {"c/tskit/trees.c", "python/tests/test_tree_stats.py",
              "c/CHANGELOG.rst", "python/CHANGELOG.rst"}
    if {a for a, _ in changed} != wanted or any(a != b for a, b in changed):
        raise RuntimeError("Patch must touch exactly the four documented files")
    return {"packet_integrity": "verified", "files": len(manifest["sha256"]),
            "scientific_audit_executed": False}


def unpack(archive, destination):
    """Reject traversal/special files; create safe relative symlinks last."""
    destination.mkdir()
    with tarfile.open(archive, "r:gz") as tar:
        members = tar.getmembers()
        roots = {PurePosixPath(m.name).parts[0] for m in members if m.name}
        if len(roots) != 1:
            raise RuntimeError("Expected one archive root")
        entries = []
        for member in members:
            path = PurePosixPath(member.name)
            if path.is_absolute() or ".." in path.parts:
                raise RuntimeError("Unsafe archive member")
            relative = Path(*path.parts[1:])
            target = destination / relative
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
            elif member.isfile():
                entries.append((member, target))
            elif member.issym():
                link = PurePosixPath(member.linkname)
                resolved = (target.parent / member.linkname).resolve()
                if link.is_absolute() or not resolved.is_relative_to(destination.resolve()):
                    raise RuntimeError("Unsafe archive symlink")
                entries.append((member, target))
            else:
                raise RuntimeError("Unsupported archive member type")
        for member, target in entries:
            if member.isfile():
                target.parent.mkdir(parents=True, exist_ok=True)
                stream = tar.extractfile(member)
                if stream is None:
                    raise RuntimeError("Archive member has no content")
                with stream, target.open("wb") as output:
                    shutil.copyfileobj(stream, output)
                target.chmod(member.mode & 0o777)
        for member, target in entries:
            if member.issym():
                target.parent.mkdir(parents=True, exist_ok=True)
                target.symlink_to(member.linkname)
    lib = destination / "python/lib"
    if not lib.is_symlink() or lib.resolve() != (destination / "c").resolve():
        raise RuntimeError("Archive must preserve python/lib -> ../c")
    if digest(destination / "c/tskit/trees.c") != BASELINE_C_SHA256:
        raise RuntimeError("Pinned baseline C SHA-256 mismatch")


def run_logged(command, cwd, logfile):
    environment = os.environ.copy()
    environment.pop("PYTHONOPTIMIZE", None)
    environment.pop("PYTHONPATH", None)
    with Path(logfile).open("w") as output:
        process = subprocess.run(command, cwd=cwd, env=environment,
                                 stdout=output, stderr=subprocess.STDOUT)
    if process.returncode:
        raise RuntimeError(f"Command failed ({process.returncode}); see {logfile}")


def native_library(build_base):
    matches = [p for p in Path(build_base).glob("lib.*")
               if (p / "tskit/__init__.py").is_file() and list(p.glob("_tskit*"))]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one native build library in {build_base}")
    return matches[0].resolve()


def replay(label, library, workdir):
    # -E neutralises inherited PYTHONOPTIMIZE/PYTHONPATH. Insert the requested
    # native library explicitly; -B prevents bytecode writes, including reuse.
    bootstrap = ("import runpy,sys; native,script=sys.argv[1:]; "
                 "sys.path.insert(0,native); sys.argv=[script]; "
                 "runpy.run_path(script,run_name='__main__')")
    receipt_path = workdir / f"{label}-audit.json"
    run_logged([sys.executable, "-E", "-B", "-c", bootstrap, str(library),
                str(ROOT / "audit.py")], workdir, receipt_path)
    receipt = json.loads(receipt_path.read_text())
    for key in ("tskit_module", "extension"):
        if not Path(receipt[key]).resolve().is_relative_to(library):
            raise RuntimeError(f"{label} imported outside requested native library")
    if receipt["cases"] != CASE_COUNT:
        raise RuntimeError(f"Expected {CASE_COUNT} audit cases")
    if receipt["passes"] + receipt["failures"] != CASE_COUNT:
        raise RuntimeError("Invalid audit accounting")
    return {k: receipt[k] for k in ("cases", "passes", "failures", "extension_sha256",
                                     "tskit_version", "tskit_module", "extension")}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-packet-only", action="store_true")
    parser.add_argument("--workdir", type=Path, default=Path(tempfile.gettempdir()) / "tskit-trait-native-replay")
    parser.add_argument("--existing-baseline-build", type=Path)
    parser.add_argument("--existing-candidate-build", type=Path)
    args = parser.parse_args()
    integrity = verify_packet()
    if args.verify_packet_only:
        print(json.dumps(integrity, indent=2))
        return
    if bool(args.existing_baseline_build) != bool(args.existing_candidate_build):
        parser.error("Both existing build directories must be supplied together")
    workdir = args.workdir.resolve()
    if workdir.exists() and any(workdir.iterdir()):
        parser.error("Choose an empty or new workdir; previous evidence is preserved")
    workdir.mkdir(parents=True, exist_ok=True)
    libraries = {}
    source_receipt = None
    if args.existing_baseline_build:
        libraries = {"baseline": native_library(args.existing_baseline_build),
                     "candidate": native_library(args.existing_candidate_build)}
        mode = "existing_native_builds"
    else:
        mode = "fresh_pinned_source_builds"
        archive = workdir / "upstream.tar.gz"
        with urllib.request.urlopen(URL, timeout=120) as source, archive.open("wb") as output:
            total = 0
            while chunk := source.read(1024 * 1024):
                total += len(chunk)
                if total > 250 * 1024 * 1024:
                    raise RuntimeError("Source archive exceeds 250 MiB limit")
                output.write(chunk)
        source_receipt = {"url": URL, "commit": COMMIT,
                          "archive_sha256": digest(archive),
                          "baseline_c_sha256": BASELINE_C_SHA256}
        for label in ("baseline", "candidate"):
            tree = workdir / f"source-{label}"
            unpack(archive, tree)
            if label == "candidate":
                run_logged(["patch", "-p1", "--batch", "-i", str(ROOT / "native.patch")],
                           tree, workdir / "patch.log")
            build = workdir / f"build-{label}"
            run_logged([sys.executable, "-E", "-B", "setup.py", "build", "--build-base",
                        str(build)], tree / "python", workdir / f"{label}-build.log")
            libraries[label] = native_library(build)
    results = {label: replay(label, lib, workdir) for label, lib in libraries.items()}
    summary = {"checked_at_utc": datetime.now(timezone.utc).isoformat(),
               "mode": mode, "python": sys.version, "python_executable": sys.executable,
               "dependencies": {name: importlib.metadata.version(name)
                                for name in ("numpy", "jsonschema", "setuptools", "wheel")},
               "source": source_receipt, "native_patch_sha256": digest(ROOT / "native.patch"),
               "results": results, "scientific_audit_executed": True,
               "bounded_reproduction_pass": results["candidate"]["failures"] == 0
                                             and results["baseline"]["failures"] >= 1}
    write_json(workdir / "summary.json", summary)
    if results["candidate"]["failures"] != 0:
        raise RuntimeError("Candidate failed exact-reference audit")
    if results["baseline"]["failures"] < 1:
        raise RuntimeError("Baseline defect did not reproduce on this platform")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
