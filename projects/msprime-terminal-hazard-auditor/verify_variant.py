#!/usr/bin/env python3
"""Replay frozen scientific inputs; precisely classify expected baseline failures."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--variant", required=True, choices=["baseline", "candidate"])
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    source, output = args.source.resolve(), args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((HERE / "manifest.json").read_text())
    expected = manifest["expected_outcomes"][args.variant]
    source_hash = manifest["baseline_demography_sha256"] if args.variant == "baseline" else manifest["frozen_patch"]["demography_sha256"]
    if digest(source / "msprime/demography.py") != source_hash:
        raise RuntimeError("Frozen source does not match selected variant")
    for item in manifest["frozen_tests"].values():
        if digest(HERE / item["path"]) != item["sha256"]:
            raise RuntimeError("Frozen scientific input hash mismatch")
    native_files = list((source / "msprime").glob("_msprime*.so")) + list((source / "msprime").glob("_msprime*.pyd"))
    if len(native_files) != 1:
        raise RuntimeError("Exactly one built native extension required")

    def run(label, command):
        with (output / (label + ".stdout")).open("w") as out, (output / (label + ".stderr")).open("w") as err:
            result = subprocess.run(command, stdout=out, stderr=err)
        return result.returncode

    def pytest_run(label, test):
        report = output / (label + ".json")
        code = run(label, [sys.executable, str(HERE / "neutral_replay.py"), "--source", str(source), "--output", str(report), "--expected-demography-sha256", source_hash, "--test", test])
        if not report.is_file():
            raise RuntimeError(label + " did not write a report")
        return code, json.loads(report.read_text())

    code, new = pytest_run("new-regressions", str(HERE / "frozen/test_terminal_contract.py"))
    counts = Counter(item["outcome"] for item in new["reports"] if item["when"] == "call")
    if (code != expected["new_exit"] or new["collected"] != expected["new_cases"]
            or len(new["reports"]) != expected["new_cases"]
            or counts != Counter(passed=expected["new_passed"], failed=expected["new_failed"])):
        raise RuntimeError("Unexpected new-regression result: " + repr((code, counts)))
    code, old = pytest_run("existing-regressions", str(source / "tests/test_demography.py") + "::TestDemographyTrajectories")
    if (code != 0 or old["collected"] != expected["existing_passed"]
            or len(old["reports"]) != expected["existing_passed"]
            or any(x["outcome"] != "passed" or x["when"] != "call" for x in old["reports"])):
        raise RuntimeError("Existing cases did not all execute and pass")

    # Preserve the frozen oracle verbatim. Bind only its platform-specific provenance
    # data beside a work copy, after validating immutable variant source hashes above.
    oracle_dir = output / "oracle-inputs"
    oracle_dir.mkdir(exist_ok=True)
    oracle = oracle_dir / "check_public_api.py"
    shutil.copyfile(HERE / "frozen/check_public_api.py", oracle)
    (oracle_dir / "expected_sources.json").write_text(json.dumps({source.name: {
        "demography_sha256": source_hash, "native_sha256": digest(native_files[0])}}, indent=2) + "\n")
    code = run("analytic-oracle", [sys.executable, str(oracle), "--source", str(source), "--output", str(output / "analytic-oracle.json")])
    analytic = json.loads((output / "analytic-oracle.json").read_text())
    failures = sum(not (item["rate_ok"] and item["survival_ok"]) for item in analytic["records"])
    if (code != 0 or analytic["cases"] != expected["analytic_cases"]
            or len(analytic["records"]) != expected["analytic_cases"]
            or analytic["failures"] != failures or failures != expected["analytic_failed"]):
        raise RuntimeError("Unexpected analytic oracle counts")
    summary = {"variant": args.variant, "matched_expected_outcome": True,
               "source": str(source), "demography_sha256": source_hash,
               "native_sha256": digest(native_files[0]), "new_cases": new["collected"],
               "new_counts": dict(counts), "new_pytest_exit": new["pytest_exit_code"],
               "existing_passed": old["collected"], "analytic_cases": analytic["cases"],
               "analytic_failed": failures, "skips_allowed": False,
               "fresh_dependency_environment_claimed": False}
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
