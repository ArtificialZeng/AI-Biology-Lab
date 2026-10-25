"""Optional numeric census replay from a local official ZIP or numeric cache.

No downloads are performed. A complete ZIP contains sequences, but this script
retains only DMS_score and explicitly configured numeric prediction columns.
Neither sequences nor mutations are written to the output directory.
"""
import sys
sys.dont_write_bytecode = True
import argparse
import hashlib
import io
import json
from pathlib import Path
import zipfile
import numpy as np
import pandas as pd
from metric_utils import PROJECT, aggregate, evaluate, output_directory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--score-zip", type=Path)
    source.add_argument("--numeric-cache", type=Path)
    parser.add_argument("--assay", help="Optional single assay; otherwise replay all 217 tables")
    args = parser.parse_args()
    manifest = json.loads((PROJECT / "provenance/download_manifest.json").read_text())
    config = json.loads((PROJECT / "provenance/model_config.json").read_text())
    models = set(config["model_list_zero_shot_substitutions_DMS"])
    members = manifest["members"]
    if args.assay:
        members = [item for item in members if Path(item["member"]).stem == args.assay]
        if len(members) != 1:
            raise ValueError("Unknown or ambiguous assay")
    archive = zipfile.ZipFile(args.score_zip) if args.score_zip else None
    records = []
    try:
        for item in members:
            assay = Path(item["member"]).stem
            if archive:
                payload = archive.read(item["member"])
                if hashlib.sha256(payload).hexdigest() != item["sha256"]:
                    raise ValueError(f"Official member hash mismatch: {assay}")
                frame = pd.read_csv(io.BytesIO(payload), usecols=lambda name: name == "DMS_score" or name in models)
                del payload
            else:
                path = args.numeric_cache / (assay + ".npz")
                if hashlib.sha256(path.read_bytes()).hexdigest() != item["numeric_cache_sha256"]:
                    raise ValueError(f"Numeric cache hash mismatch: {assay}")
                with np.load(path, allow_pickle=False) as saved:
                    frame = pd.DataFrame(saved["data"], columns=saved["columns"])
                frame = frame[[name for name in frame.columns if name == "DMS_score" or name in models]]
            records.extend(evaluate(frame, assay, models))
    finally:
        if archive:
            archive.close()
    cells = pd.DataFrame(records)
    expected = pd.read_csv(PROJECT / "results/cell_metrics.csv")
    expected = expected[expected.assay.isin(cells.assay.unique())]
    compare = cells.merge(expected, on=["assay", "model"], suffixes=("", "_saved"), validate="one_to_one")
    for field in ["original", "reversed", "shuffle_seed_1", "tie_average", "tie_min", "tie_max"]:
        if not np.allclose(compare[field], compare[field + "_saved"], atol=2e-12, rtol=0):
            raise ValueError(f"Replay differs in {field}; check recorded NumPy version")
    destination = output_directory("census-replay" if not args.assay else "assay-replay")
    cells.to_csv(destination / "cell_metrics.csv", index=False)
    aggregate(cells, pd.read_csv(PROJECT / "provenance/assay_groups.csv")).to_csv(destination / "aggregate_metrics.csv")
    result = {"status": "passed", "assays": cells.assay.nunique(), "cells": len(cells),
              "numpy": np.__version__, "pandas": pd.__version__,
              "scope": "complete frozen archive" if not args.assay else "one assay only"}
    (destination / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
