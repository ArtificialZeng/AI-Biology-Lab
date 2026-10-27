"""Read-only integrity and exact-algebra replay; no third-party dependencies."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def main():
    root = Path(__file__).resolve().parent
    manifest = json.loads((root / "manifest.json").read_text())
    expected = {"README.md", "proof.md", "sources.md", "witness.json", "check.py", "reproduce.py", "status.json", "LICENSE"}
    if set(manifest["files"]) != expected:
        raise ValueError("manifest file-set mismatch")
    for name, digest in manifest["files"].items():
        if hashlib.sha256((root / name).read_bytes()).hexdigest() != digest:
            raise ValueError("integrity mismatch: " + name)
    subprocess.run([sys.executable, "-I", "-B", str(root / "check.py")], check=True)
    print("Packet integrity and exact finite algebra replay passed; analytic proofs require review.")


if __name__ == "__main__":
    main()
