"""Read-only complete-file inventory/hash guard; checks remain active under -O."""
from pathlib import Path
from hashlib import sha256
import json
import os
import shutil
import sys
import tempfile

sys.dont_write_bytecode = True

def require(truth, message):
    if not truth:
        raise RuntimeError(message)

def inventory(base):
    base = Path(base)
    require(not base.is_symlink(), "Symlink package root")
    files = {}
    for directory, dirs, names in os.walk(base, followlinks=False):
        here = Path(directory)
        for name in dirs + names:
            require(not (here/name).is_symlink(), "Symlink forbidden: "+str((here/name).relative_to(base)))
        for name in names:
            path = here/name
            rel = path.relative_to(base).as_posix()
            if rel == "manifest.json" or path.suffix in (".pyc", ".pyo"):
                continue
            require(path.is_file(), "Nonregular package entry: "+rel)
            files[rel] = sha256(path.read_bytes()).hexdigest()
    return dict(sorted(files.items()))

def verify_manifest(base):
    base = Path(base)
    require(not (base/"manifest.json").is_symlink(), "Symlink manifest forbidden")
    data = json.loads((base/"manifest.json").read_text())
    expected = data.get("files")
    require(isinstance(expected, dict) and bool(expected), "Empty or invalid manifest inventory")
    actual = inventory(base)
    require(set(actual) == set(expected), "Manifest inventory mismatch; missing="+str(sorted(set(expected)-set(actual)))+" unexpected="+str(sorted(set(actual)-set(expected))))
    for rel, value in actual.items():
        require(expected[rel] == value, "Package hash mismatch: "+rel)
    digest = sha256(json.dumps(expected, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    require(data.get("digest") == digest, "Manifest mapping digest mismatch")
    return True

def self_test(base):
    verify_manifest(base)
    probes = []
    mutations = ["SOURCE.md", "PROOF.md", "evidence/exact-verification.json", "unexpected-file", "missing-file", "symlink"]
    for kind in mutations:
        with tempfile.TemporaryDirectory(prefix="delta4-integrity-") as td:
            target = Path(td)/"bundle"
            shutil.copytree(base, target, ignore=shutil.ignore_patterns("*.pyc", "*.pyo"))
            if kind == "unexpected-file":
                (target/"extra.txt").write_text("unexpected")
            elif kind == "missing-file":
                (target/"SOURCE.md").unlink()
            elif kind == "symlink":
                (target/"SOURCE.md").unlink()
                (target/"SOURCE.md").symlink_to("PROOF.md")
            else:
                path = target/kind
                path.write_bytes(path.read_bytes()+b"\nchanged\n")
            try:
                verify_manifest(target)
            except RuntimeError:
                probes.append({"probe": kind, "rejected": True})
            else:
                raise RuntimeError("Integrity guard accepted mutation: "+kind)
    return {"integrity_verified": True, "all_archived_files_in_inventory": True,
            "source_proof_archive_and_inventory_mutations_rejected": probes,
            "optimization_mode_preserves_checks": True}

if __name__ == "__main__":
    base = Path(__file__).resolve().parent
    result = self_test(base) if "--self-test" in sys.argv else {"integrity_verified": verify_manifest(base)}
    print(json.dumps(result, indent=2))
