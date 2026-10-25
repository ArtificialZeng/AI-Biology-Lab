#!/usr/bin/env python3
"""Portable reproduction in an isolated copy; archived evidence is read-only."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT / "projects" / "mixed-phosphorylation"


def check_archive() -> int:
    manifest = json.loads((PROJECT / "archive_manifest.json").read_text())
    for name, expected in manifest["sha256"].items():
        observed = hashlib.sha256((PROJECT / name).read_bytes()).hexdigest()
        if observed != expected:
            raise RuntimeError(f"Archive integrity mismatch: {name}")
    print(f"Archive integrity: {len(manifest['sha256'])} files match", flush=True)
    return len(manifest["sha256"])


def validate_numerics(work: Path) -> dict:
    """Check numerical reproduction against archived observables, not a proof."""
    data = json.loads((work / "evidence/independent_numerics.json").read_text())
    reference = json.loads((PROJECT / "evidence/independent_numerics.json").read_text())
    trajectories = data["trajectory_checks"]
    shots = data["periodic_shooting"]
    if len(trajectories) != 8 or len(shots) != 4:
        raise RuntimeError("Numerical reproduction did not complete all trajectories and shooting cases")
    for row in trajectories:
        if (not row["success"] or row["max_conservation_error"] > 1e-6
                or row["minimum_sampled_species"] < -1e-7):
            raise RuntimeError("Numerical integration failed a conservation/positivity check")
    errors = []
    for row, ref in zip(shots, reference["periodic_shooting"]):
        if row["K"] != ref["K"] or not row["root_success"]:
            raise RuntimeError("Periodic shooting did not converge for the expected case")
        for name in ["period", "S2_min", "S2_max"]:
            error = abs(row[name] - ref[name])
            if not math.isfinite(error) or error > 1e-6:
                raise RuntimeError(f"Numerical reference mismatch for {name}, K={row['K']}")
            errors.append(error)
        if (row["shooting_residual_infinity_norm"] > 1e-7
                or row["neutral_multiplier_error"] > 1e-6
                or not 0 <= row["largest_nontrivial_multiplier_abs"] < 1):
            raise RuntimeError("Numerical shooting residual or stability check failed")
    return {"scope": "floating-point reproduction only; not a periodic-orbit proof",
            "trajectories_checked": len(trajectories), "shooting_cases_checked": len(shots),
            "observable_tolerance": 1e-6, "max_observable_reference_error": max(errors)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["verify", "regenerate", "numerics", "figure"])
    args = parser.parse_args()
    if sys.flags.optimize or os.environ.get("PYTHONOPTIMIZE"):
        parser.error("Verification requires Python assertions; unset PYTHONOPTIMIZE and do not use -O.")
    output = ROOT / "artifacts" / "mixed-phosphorylation" / args.mode
    output.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    result = {"mode": args.mode, "status": "running", "commands": [],
              "python_version": sys.version.split()[0]}
    verify = ["manuscript/verify_certificate.py", "src/verify_finite_kinase.py",
              "audit/verify_reduction.py"]
    plans = {
        "verify": [[s] for s in verify],
        "regenerate": [["src/finite_kinase_exact.py"], ["src/finite_kinase_signs.py"]]
                      + [[s] for s in verify],
        "numerics": [["src/baseline.py", "100", "1000", "10000", "inf"],
                     ["src/independent_numeric.py"], ["src/period_response.py"],
                     ["src/plot_results.py"]],
        "figure": [["src/plot_results.py"]],
    }
    try:
        result["archive_files_checked"] = check_archive()
        with tempfile.TemporaryDirectory(prefix="ai-biology-") as temporary:
            work = Path(temporary) / "project"
            shutil.copytree(PROJECT, work, ignore=shutil.ignore_patterns("__pycache__", "build", "build-*", "CAPD"))
            for index, command in enumerate(plans[args.mode], start=1):
                command_started = time.monotonic()
                print(f"[{index}/{len(plans[args.mode])}] python {' '.join(command)}", flush=True)
                log = output / f"{index:02d}-{Path(command[0]).stem}.log"
                environment = dict(os.environ, MPLBACKEND="Agg", PYTHONUNBUFFERED="1")
                with log.open("w") as stream:
                    process = subprocess.Popen([sys.executable, *command], cwd=work,
                                               env=environment, stdout=subprocess.PIPE,
                                               stderr=subprocess.STDOUT, text=True)
                    for line in process.stdout:
                        stream.write(line)
                        stream.flush()
                        # Numerical scripts emit large JSON payloads; keep them in the log.
                        if args.mode in {"verify", "regenerate", "figure"}:
                            print(line, end="", flush=True)
                    code = process.wait()
                result["commands"].append({"command": ["python", *command], "exit_code": code,
                                            "seconds": time.monotonic() - command_started,
                                            "log": log.name})
                if code:
                    raise RuntimeError(f"Command failed ({code}); inspect {log.relative_to(ROOT)}")
                print(f"Completed in {time.monotonic() - command_started:.1f}s", flush=True)
            if args.mode == "regenerate":
                for name in ["finite_kinase_algebra.json", "finite_kinase_sign_certificate.json"]:
                    fresh = json.loads((work / "evidence" / name).read_text())
                    archived = json.loads((PROJECT / "evidence" / name).read_text())
                    if fresh != archived:
                        raise RuntimeError(f"Regenerated certificate differs: {name}")
                result["regenerated_certificates_match"] = True
            if args.mode == "numerics":
                result["numerical_checks"] = validate_numerics(work)
            shutil.copytree(work / "evidence", output / "evidence", dirs_exist_ok=True)
            for name in ["manuscript/certificate_check.json", "audit/reduction_verification.json"]:
                source = work / name
                if source.exists():
                    shutil.copy2(source, output / source.name)
        result["status"] = "passed"
    except Exception as error:
        result["status"] = "failed"
        result["error"] = str(error)
    result["seconds"] = time.monotonic() - started
    (output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(f"{args.mode}: {result['status']} ({result['seconds']:.1f}s)", flush=True)
    print(f"Results: {output.relative_to(ROOT) / 'result.json'}", flush=True)
    if result["status"] != "passed":
        print(result["error"], file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
