"""Offline real-score witness; outputs go to repository-root artifacts/."""
import sys
sys.dont_write_bytecode = True
import argparse
import platform
import hashlib
import json
import numpy as np
import pandas as pd
from metric_utils import PROJECT, METRICS, evaluate, normalize_truth, official_function, output_directory, tie_statistics


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--strict-reference',action='store_true',help='Also require exact archived ordinal outputs')
    args=parser.parse_args()
    metadata = json.loads((PROJECT / "fixtures/manifest.json").read_text())
    path = PROJECT / "fixtures/GAL4_numeric.csv"
    if hashlib.sha256(path.read_bytes()).hexdigest() != metadata["fixture_sha256"]:
        raise ValueError("Fixture hash mismatch")
    frame = pd.read_csv(path)
    if frame.columns.tolist() != metadata["columns"] or len(frame) != metadata["rows"]:
        raise ValueError("Fixture schema mismatch")
    models = set(metadata["columns"]) - {"DMS_score"}
    official = official_function()
    records = evaluate(frame, metadata["assay"], models)
    gain, weights, ideal = normalize_truth(frame.DMS_score.to_numpy())
    permutations = [np.arange(len(frame)), np.arange(len(frame))[::-1], np.random.default_rng(1).permutation(len(frame))]
    for row in records:
        if row["status"] != "ok":
            raise ValueError("Sample contains nonfinite predictions")
        scores = frame[row["model"]].to_numpy()
        for label, permutation in zip(METRICS[:3], permutations):
            actual = float(official(frame.DMS_score.to_numpy()[permutation], scores[permutation]))
            if not np.isclose(actual, row[label], atol=2e-12, rtol=0):
                raise ValueError("Independent ordinal implementation differs from pinned function")
            stat = tie_statistics(gain[permutation], scores[permutation], weights, ideal)
            for field in ["tie_average", "tie_min", "tie_max"]:
                if not np.isclose(stat[field], row[field], atol=2e-12, rtol=0):
                    raise ValueError("Tie statistics changed under row permutation")
    expected = pd.read_csv(PROJECT / "results/cell_metrics.csv")
    expected = expected[expected.assay == metadata["assay"]].set_index("model")
    ordinal_comparison=[]
    for row in records:
        for field in METRICS:
            saved=float(expected.loc[row["model"],field])
            same=bool(np.isclose(row[field],saved,atol=2e-12,rtol=0))
            if field in ["tie_average","tie_min","tie_max"]:
                if not same: raise ValueError("Order-free sample statistic differs from archived census")
            else:
                ordinal_comparison.append({"model":row["model"],"metric":field,"archived":saved,"local":row[field],"difference":row[field]-saved,"matches_tolerance":same})
                if args.strict_reference and not same:
                    raise ValueError("Ordinal output differs under strict archived-reference mode")
    # Guard checks: do not silently normalize constant, NaN or overflowing truth.
    for bad in [np.ones(12), np.r_[np.nan, np.arange(11)], np.r_[-1e308, np.zeros(10), 1e308]]:
        try:
            normalize_truth(bad)
        except ValueError:
            pass
        else:
            raise ValueError("Invalid experimental values were accepted")
    destination = output_directory("sample")
    result = {"status": "passed", "numpy": np.__version__, "pandas": pd.__version__,
              "scope": "Numeric-only GAL4 fixture; primary Spearman conclusions not tested", "models": records,
              "environment":{"system":platform.system(),"machine":platform.machine(),"python":platform.python_version()},
              "ordinal_reference_comparison":ordinal_comparison}
    (destination / "sample_result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": "passed", "assay": metadata["assay"], "rows": len(frame),
                      "models": len(records), "ordinal_reference_differences":sum(not x["matches_tolerance"] for x in ordinal_comparison),
                      "max_ordinal_reference_absolute_difference":max(abs(x["difference"]) for x in ordinal_comparison),
                      "output": str(destination / "sample_result.json")}))


if __name__ == "__main__":
    main()
