#!/usr/bin/env python3
"""Build pinned upstream sources in a separate work directory; never patch research trees."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.request

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(command, *, cwd, env):
    print(json.dumps({"command": [str(x) for x in command], "cwd": str(cwd)}), flush=True)
    subprocess.run([str(x) for x in command], cwd=cwd, env=env, check=True)


def extract(archive, destination):
    with tarfile.open(archive) as stream:
        members = stream.getmembers()
        roots = {Path(m.name).parts[0] for m in members if Path(m.name).parts}
        if len(roots) != 1:
            raise ValueError("Archive must contain exactly one root")
        # data rejects absolute paths, traversal, devices and escaping links.
        stream.extractall(destination, filter="data")
    return destination / roots.pop()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--variant", choices=["baseline", "candidate"], default="baseline")
    parser.add_argument("--reuse-python", type=Path, help="Use an existing interpreter without installing packages")
    parser.add_argument("--gsl-prefix", type=Path, help="Installed GSL include/ and lib/ prefix")
    args = parser.parse_args()
    manifest = json.loads((HERE / "manifest.json").read_text())
    work = args.work.resolve() if args.work else Path(tempfile.mkdtemp(prefix="msprime-replay-"))
    if work == HERE or work.is_relative_to(HERE):
        parser.error("--work must be outside the public package")
    work.mkdir(parents=True, exist_ok=True)
    sources = work / "sources"
    if sources.exists():
        parser.error("Existing sources/ refused: choose a fresh work directory")
    sources.mkdir()
    source_dirs = {}
    for item in manifest["sources"]:
        archive = work / (item["name"] + ".tar.gz")
        if not archive.exists():
            request = urllib.request.Request(item["url"], headers={"User-Agent": "msprime-replay-bootstrap"})
            with urllib.request.urlopen(request, timeout=60) as response:
                archive.write_bytes(response.read())
        if digest(archive) != item["sha256"]:
            raise ValueError("Archive checksum mismatch: " + item["name"])
        target = sources / item["name"]
        target.mkdir()
        source_dirs[item["name"]] = extract(archive, target)
    source = source_dirs["msprime"]
    gitlink = source / "git-submodules" / "tskit"
    if gitlink.exists():
        if list(gitlink.iterdir()):
            raise ValueError("Unexpected populated gitlink in upstream archive")
        gitlink.rmdir()
    # Preserve upstream internal links rather than recursively following test fixtures.
    shutil.copytree(source_dirs["tskit"], gitlink, symlinks=True)
    if digest(source / "msprime" / "demography.py") != manifest["baseline_demography_sha256"]:
        raise ValueError("Baseline source hash differs from frozen source")
    if not (gitlink / "c/subprojects/kastore/kastore.c").is_file():
        raise ValueError("Pinned tskit archive is missing the native kastore source")
    if args.variant == "candidate":
        frozen = manifest.get("frozen_patch")
        if not frozen:
            raise ValueError("Candidate build refused: root has not frozen a patch")
        patch = HERE / frozen["path"]
        if digest(patch) != frozen["sha256"]:
            raise ValueError("Frozen patch hash mismatch")
        # Avoid discovering an unrelated ancestor Git checkout around --work.
        run(["git", "init", "--quiet"], cwd=source, env=os.environ.copy())
        run(["git", "apply", "--check", patch], cwd=source, env=os.environ.copy())
        run(["git", "apply", patch], cwd=source, env=os.environ.copy())
        if digest(source / "msprime/demography.py") != frozen["demography_sha256"]:
            raise ValueError("Patched demography source hash mismatch")
    payload = {"status": "prepared-only", "source": str(source), "work": str(work),
               "variant": args.variant, "native_build_executed": False,
               "base_commit": manifest["base_commit"],
               "tskit_submodule_commit": manifest["tskit_submodule_commit"]}
    (work / "bootstrap-result.json").write_text(json.dumps(payload, indent=2) + "\n")
    if args.prepare_only:
        print(json.dumps(payload, indent=2))
        return
    env = os.environ.copy()
    if args.gsl_prefix:
        prefix = args.gsl_prefix.resolve()
        for child in ["include/gsl/gsl_version.h", "lib"]:
            if not (prefix / child).exists():
                raise ValueError("Invalid --gsl-prefix: missing " + child)
        env["CFLAGS"] = env.get("CFLAGS", "") + " " + shlex.quote("-I" + str(prefix / "include"))
        env["LDFLAGS"] = env.get("LDFLAGS", "") + " " + shlex.quote("-L" + str(prefix / "lib")) + " " + shlex.quote("-Wl,-rpath," + str(prefix / "lib"))
    elif not shutil.which("gsl-config"):
        raise RuntimeError("Install GSL development files or provide --gsl-prefix")
    if args.reuse_python:
        # Keep the venv entry path: resolving its symlink loses pyvenv.cfg detection.
        python = Path(os.path.abspath(args.reuse_python))
    else:
        run([sys.executable, "-m", "venv", work / "venv"], cwd=work, env=env)
        python = work / "venv" / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
        run([python, "-m", "pip", "install", "-r", HERE / "requirements-observed.txt"], cwd=work, env=env)
    env["SETUPTOOLS_SCM_PRETEND_VERSION"] = "1.4.5.dev2"
    run([python, "setup.py", "build_ext", "--inplace"], cwd=source, env=env)
    payload.update(status="native-build-completed", python=str(python), native_build_executed=True,
                   demography_sha256=digest(source / "msprime/demography.py"))
    (work / "bootstrap-result.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
