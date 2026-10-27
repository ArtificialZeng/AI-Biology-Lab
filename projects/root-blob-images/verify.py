"""Read-only, relocation-safe replay of the original exact evidence packet."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent

def require(test, message):
    if not test:
        raise ValueError(message)

def hashes():
    files = {}
    for path in sorted(HERE.rglob("*")):
        require(not path.is_symlink(), "Symlink is outside archive contract")
        if path.is_file() and path != HERE/"archive_manifest.json":
            files[path.relative_to(HERE).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    return files

def run_child(filename, *args):
    env = dict(os.environ, PYTHONOPTIMIZE="0", PYTHONDONTWRITEBYTECODE="1")
    with tempfile.TemporaryDirectory(prefix="root-blob-replay-") as work:
        result = subprocess.run([sys.executable, "-B", str(HERE/filename), *args],
                                cwd=work, env=env, text=True, capture_output=True,
                                timeout=120, check=False)
    require(result.returncode == 0, f"{filename} failed: {result.stderr}")
    return json.loads(result.stdout)

def main():
    require(sys.argv[1:] in ([], ["verify"]), "Use verify")
    expected = json.loads((HERE/"archive_manifest.json").read_text())["files"]
    before = hashes()
    require(before == expected, "Archive file list or digest changed")
    exact = run_child("check_exact.py", "verify")
    require(exact["exact_verification"] is True, "Exact verification did not succeed")
    graph = run_child("check_graph.py")
    require(graph == json.loads((HERE/"graph-certificate.json").read_text()),
            "Graph certificate changed")
    require(hashes() == before, "Replay modified archive inputs")
    print(json.dumps({"verified": True, "manifest_files": len(before),
                      "target_leaf_patterns": exact["target_leaf_patterns"],
                      "constructed_cases": exact["constructed_cases"],
                      "simple_rootings_checked": graph["rooting_count"],
                      "displayed_paths_checked": graph["displayed_path_count"],
                      "scientific_child_optimization": 0,
                      "scope": "Finite exact evidence; general claims require proof.md"},
                     sort_keys=True))

if __name__ == "__main__":
    main()
