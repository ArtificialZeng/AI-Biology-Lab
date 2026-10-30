#!/usr/bin/env python3
"""Run selected real pytest cases with native import provenance, without an outcome oracle."""
import argparse
import datetime
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-demography-sha256")
    parser.add_argument("--test", action="append", help="pytest node ID; default existing trajectory class")
    args = parser.parse_args()
    source = args.source.resolve()
    sys.path.insert(0, str(source))
    import msprime
    import msprime._msprime
    import pytest

    init = Path(msprime.__file__).resolve()
    native = Path(msprime._msprime.__file__).resolve()
    demography = Path(msprime.demography.__file__).resolve()
    if not init.is_relative_to(source) or not native.is_relative_to(source):
        raise RuntimeError("Imported modules are not from requested native source tree")
    if demography != source / "msprime/demography.py":
        raise RuntimeError("Demography import provenance mismatch")
    if native.suffix not in [".so", ".pyd"]:
        raise RuntimeError("Native extension suffix required")
    source_hash = digest(demography)
    if args.expected_demography_sha256 and source_hash != args.expected_demography_sha256:
        raise RuntimeError("Frozen source hash mismatch")
    tests = args.test or [str(source / "tests/test_demography.py") + "::TestDemographyTrajectories"]
    payload = {"started_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "python": sys.version, "module": str(init), "native": str(native),
              "native_sha256": digest(native), "demography_sha256": source_hash,
              "tests": tests, "collected": 0, "reports": [],
              "packages": {name: importlib.metadata.version(name) for name in
                           ["numpy", "scipy", "tskit", "demes", "pytest"]}}

    class Recorder:
        def pytest_collection_finish(self, session):
            payload["collected"] = len(session.items)

        def pytest_runtest_logreport(self, report):
            if report.when == "call" or report.failed:
                payload["reports"].append({"nodeid": report.nodeid, "when": report.when,
                                           "outcome": report.outcome})

    code = int(pytest.main(["-o", "addopts=", "-p", "no:cacheprovider", "-q", *tests], plugins=[Recorder()]))
    payload.update(pytest_exit_code=code,
                  completed_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  scope="Only the collected named cases; no whole-suite or CI claim")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    raise SystemExit(code)


if __name__ == "__main__":
    main()
