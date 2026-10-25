"""Verify shipped result integrity and reconstruct hierarchical aggregates."""
import sys
sys.dont_write_bytecode = True
import hashlib
import json
import numpy as np
import pandas as pd
from metric_utils import PROJECT, METRICS, aggregate, output_directory


def main():
    manifest = json.loads((PROJECT / "provenance/file_hashes.json").read_text())
    for name, expected in manifest.items():
        if hashlib.sha256((PROJECT / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"File hash mismatch: {name}")
    cells = pd.read_csv(PROJECT / "results/cell_metrics.csv")
    summary = json.loads((PROJECT / "results/census_summary.json").read_text())
    if cells.duplicated(["assay", "model"]).any():
        raise ValueError("Duplicate assay/model cells")
    valid = cells[cells.status == "ok"]
    if not np.isfinite(valid[METRICS].to_numpy()).all():
        raise ValueError("Nonfinite metrics")
    if (len(valid), cells.assay.nunique(), cells.model.nunique()) != (20615, 217, 95):
        raise ValueError("Census coverage differs")
    if int(valid.groupby("assay").n.first().sum()) != summary["mutation_rows_in_distinct_assays"]:
        raise ValueError("Assay row count differs")
    metadata = pd.read_csv(PROJECT / "provenance/assay_groups.csv")
    reproduced = aggregate(cells, metadata)
    saved = pd.read_csv(PROJECT / "results/aggregate_metrics.csv", index_col="model")
    if not np.allclose(reproduced.loc[saved.index, METRICS], saved[METRICS], atol=2e-12, rtol=0):
        raise ValueError("Hierarchical aggregates differ")
    for metric in ["reversed", "shuffle_seed_1", "tie_average"]:
        ranks = reproduced[metric].rank(method="min", ascending=False).astype(int)
        if not ranks.loc[saved.index].equals(saved[metric + "_rank"]):
            raise ValueError(f"Ranks differ for {metric}")
    top = reproduced.original.nlargest(5).index.tolist()
    if any(reproduced[field].nlargest(5).index.tolist() != top for field in ["reversed", "shuffle_seed_1", "tie_average"]):
        raise ValueError("Top-five statement differs")
    result = {"status": "passed", "assays": 217, "models": 95, "finite_cells": len(valid),
              "scope": "Recomputes aggregates from shipped cells; does not regenerate all model scores"}
    (output_directory("verification") / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
