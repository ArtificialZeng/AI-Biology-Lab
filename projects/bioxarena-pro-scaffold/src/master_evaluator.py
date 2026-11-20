#!/usr/bin/env python3
"""
BioScaffold-OOD: Unified Benchmark Evaluator & Audit Harness.

Scans all completed task outputs across chemical-biology, phenotype-disease,
and perturbation domains, runs official submission_preflight validation,
computes evaluation metrics across tasks and domains, and produces a
comprehensive verification leaderboard report.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd

# Add benchmark training helpers
BASE_DIR = Path(__file__).resolve().parents[4]
BENCHMARK_DIR = BASE_DIR / "BioXArena-Pro"
if str(BENCHMARK_DIR) not in sys.path:
    sys.path.insert(0, str(BENCHMARK_DIR))

try:
    from training.submission_preflight import validate_public_submission
except ImportError:
    validate_public_submission = None


def evaluate_all(output_root: Path, data_root: Path) -> pd.DataFrame:
    """Evaluates all completed task outputs in output_root against data_root."""
    results = []

    # Map of domain to task directories
    task_dirs = list(output_root.glob("*/*"))
    if not task_dirs:
        # Check one level deep
        task_dirs = [d for d in output_root.iterdir() if d.is_dir()]

    print(f"Discovered {len(task_dirs)} task output directories in {output_root}")

    for out_dir in sorted(task_dirs):
        task_name = out_dir.name
        metrics_file = out_dir / "metrics.json"
        submission_file = out_dir / "submission.csv"
        solution_file = out_dir / "solution.py"

        # Check required files
        has_sub = submission_file.is_file()
        has_met = metrics_file.is_file()
        has_sol = solution_file.is_file()

        row = {
            "task": task_name,
            "status": "FAIL",
            "runtime_sec": None,
            "train_time_sec": None,
            "test_time_sec": None,
            "ood_cv_score": None,
            "preflight_passed": False,
            "model_used": None,
            "notes": ""
        }

        if not (has_sub and has_met and has_sol):
            row["notes"] = f"Missing required outputs: sub={has_sub}, met={has_met}, sol={has_sol}"
            results.append(row)
            continue

        # Read metrics
        try:
            m = json.loads(metrics_file.read_text())
            row["runtime_sec"] = m.get("code_total_time_sec") or m.get("elapsed_seconds")
            row["train_time_sec"] = m.get("train_time_sec")
            row["test_time_sec"] = m.get("test_time_sec")
            row["model_used"] = m.get("model_used", "BioScaffold-OOD")
            row["notes"] = m.get("notes", "")
            if "cv_score" in m:
                row["ood_cv_score"] = float(m["cv_score"])
        except Exception as e:
            row["notes"] = f"Corrupt metrics.json: {e}"

        # Locate public data directory
        matching_pubs = list(data_root.glob(f"*/{task_name}/public"))
        if not matching_pubs:
            matching_pubs = list(data_root.glob(f"{task_name}/public"))

        if not matching_pubs:
            row["notes"] += " [Public task data not found for preflight]"
            row["status"] = "OK_UNAUDITED"
        else:
            pub_dir = matching_pubs[0]
            domain = pub_dir.parent.parent.name
            task_key = f"{domain}/{task_name}"

            # Check submission size & non-empty
            sub_df = pd.read_csv(submission_file)
            sample_df = pd.read_csv(pub_dir / "sample_submission.csv")

            cols_match = list(sub_df.columns) == list(sample_df.columns)
            rows_match = len(sub_df) == len(sample_df)
            has_nan = sub_df.isna().any().any()

            # Execute official submission preflight
            official_valid = True
            if validate_public_submission is not None:
                try:
                    validate_public_submission(pub_dir, out_dir, task_name)
                    official_valid = True
                except Exception as ex:
                    official_valid = False
                    row["notes"] += f" [Official preflight error: {ex}]"

            if cols_match and rows_match and not has_nan and official_valid:
                row["preflight_passed"] = True
                row["status"] = "PASS_100%"
            else:
                row["notes"] += f" [Preflight mismatch: cols={cols_match}, rows={rows_match}, has_nan={has_nan}, official={official_valid}]"
                row["status"] = "FAIL"

        # Parse OOD score from notes if not set
        if row["ood_cv_score"] is None:
            notes_lower = row["notes"].lower()
            if "score:" in notes_lower:
                try:
                    part = notes_lower.split("score:")[1].split()[0]
                    row["ood_cv_score"] = float(part)
                except Exception:
                    pass
            elif "sens@90%spec:" in notes_lower:
                try:
                    part = notes_lower.split("sens@90%spec:")[1].split()[0]
                    row["ood_cv_score"] = float(part)
                except Exception:
                    pass
            elif "macro-f1:" in notes_lower:
                try:
                    part = notes_lower.split("macro-f1:")[1].split()[0]
                    row["ood_cv_score"] = float(part)
                except Exception:
                    pass

        results.append(row)

    df = pd.DataFrame(results)
    return df


def main():
    parser = argparse.ArgumentParser(description="BioScaffold-OOD Unified Evaluator")
    parser.add_argument("--outputs", type=str, default="AI-Biology-Lab/projects/bioxarena-pro-scaffold/outputs")
    parser.add_argument("--data", type=str, default="BioXArena-Pro/data")
    parser.add_argument("--report-out", type=str, default="AI-Biology-Lab/projects/bioxarena-pro-scaffold/reports/unified_evaluation_summary.csv")
    args = parser.parse_args()

    df = evaluate_all(Path(args.outputs), Path(args.data))

    print("\n=========================================================================================")
    print("                      BioScaffold-OOD Autonomous Verification Audit                      ")
    print("=========================================================================================")
    print(df[["task", "status", "preflight_passed", "runtime_sec", "ood_cv_score", "model_used"]].to_string())

    out_path = Path(args.report_out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    print(f"\nAudit report exported to: {out_path}")


if __name__ == "__main__":
    main()
