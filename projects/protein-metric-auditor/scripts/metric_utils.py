"""Numeric-only NDCG utilities. Tie expectation/extrema are standard methods."""
from pathlib import Path
import ast
import hashlib
import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
REPOSITORY = PROJECT.parent.parent if PROJECT.parent.name == "projects" else PROJECT
OFFICIAL_SHA256 = "0002bc1b031c02dc2d67fde53da092c8fd301415645e8e83b425b4183570eba2"
METRICS = ["original", "reversed", "shuffle_seed_1", "tie_average", "tie_min", "tie_max"]


def output_directory(name):
    path = REPOSITORY / "artifacts" / "protein-metric-auditor" / name
    path.mkdir(parents=True, exist_ok=True)
    return path


def official_function():
    path = PROJECT / "vendor/proteingym/performance_DMS_benchmarks.py"
    content = path.read_bytes()
    if hashlib.sha256(content).hexdigest() != OFFICIAL_SHA256:
        raise ValueError("Pinned official evaluator hash mismatch")
    nodes = [node for node in ast.parse(content).body
             if isinstance(node, ast.FunctionDef) and node.name in {"minmax", "calc_ndcg"}]
    if {node.name for node in nodes} != {"minmax", "calc_ndcg"}:
        raise ValueError("Pinned evaluator functions missing")
    scope = {"np": np, "pd": pd}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), "exec"), scope)
    return scope["calc_ndcg"]


def normalize_truth(y):
    y = np.asarray(y, dtype=float)
    if y.ndim != 1 or len(y) < 10 or not np.isfinite(y).all():
        raise ValueError("Need at least ten finite scalar experimental values")
    with np.errstate(over="ignore", invalid="ignore"):
        span = np.ptp(y)
    if not np.isfinite(span) or span <= 0:
        raise ValueError("Experimental span must be finite and positive")
    gain = (y - y.min()) / span
    k = int(np.floor(len(y) * .1))
    weights = 1 / np.log2(np.arange(2, k + 2))
    ideal = float(np.dot(np.sort(gain)[::-1][:k], weights))
    if not np.isfinite(gain).all() or not np.isfinite(ideal) or ideal <= 0:
        raise ValueError("Normalized gain/ideal DCG is invalid")
    return gain, weights, ideal


def ordinal_ndcg(gain, scores, weights, ideal):
    order = np.argsort(-scores)
    return float(np.dot(gain[order[:len(weights)]], weights) / ideal)


def tie_statistics(gain, scores, weights, ideal):
    """Expectation and bounds over permutations inside exact-score tie groups.

    A group crossing k is averaged over the whole group, including its rows
    outside k. Each model's extremal order can differ from other models' orders.
    """
    n, k = len(scores), len(weights)
    order = np.argsort(-scores)
    sorted_scores, ranked_gain = scores[order], gain[order]
    starts = np.r_[0, np.flatnonzero(sorted_scores[1:] != sorted_scores[:-1]) + 1]
    ends = np.r_[starts[1:], n]
    tied = ends - starts > 1
    dcg = float(np.dot(ranked_gain[:k], weights))
    lo = hi = average = dcg
    relevant = 0
    for start, end in zip(starts[tied & (starts < k)], ends[tied & (starts < k)]):
        stop = min(end, k)
        weight = weights[start:stop]
        old = float(np.dot(ranked_gain[start:stop], weight))
        values = np.sort(ranked_gain[start:end])
        width = stop - start
        lo += float(np.dot(values[:width], weight)) - old
        hi += float(np.dot(values[::-1][:width], weight)) - old
        average += float(values.mean() * weight.sum()) - old
        relevant += 1
    return {"original": dcg / ideal, "tie_average": average / ideal,
            "tie_min": lo / ideal, "tie_max": hi / ideal,
            "tie_range": (hi - lo) / ideal,
            "tied_observations": int(np.sum((ends - starts)[tied])),
            "relevant_tie_groups": relevant,
            "unique_prediction_values": len(starts)}


def evaluate(frame, assay, models):
    gain, weights, ideal = normalize_truth(frame["DMS_score"].to_numpy(dtype=float))
    n = len(gain)
    permutations = {"reversed": np.arange(n)[::-1],
                    "shuffle_seed_1": np.random.default_rng(1).permutation(n)}
    records = []
    for model in frame.columns:
        if model not in models:
            continue
        scores = frame[model].to_numpy(dtype=float)
        if not np.isfinite(scores).all():
            records.append({"assay": assay, "model": model, "n": n,
                            "status": "nonfinite_prediction"})
            continue
        result = tie_statistics(gain, scores, weights, ideal)
        for label, permutation in permutations.items():
            result[label] = ordinal_ndcg(gain[permutation], scores[permutation], weights, ideal)
        result["common_permutation_range"] = max(result[key] for key in METRICS[:3]) - min(result[key] for key in METRICS[:3])
        if not np.isfinite([result[key] for key in METRICS]).all():
            raise ValueError("Nonfinite metric output")
        records.append({"assay": assay, "model": model, "n": n, "status": "ok", **result})
    return records


def aggregate(cells, metadata):
    valid = cells[cells.status == "ok"].copy()
    merged = valid.merge(metadata, left_on="assay", right_on="DMS_id", how="left", validate="many_to_one")
    if merged.UniProt_ID.isna().any():
        raise ValueError("Assay metadata missing")
    return (merged.groupby(["model", "UniProt_ID", "coarse_selection_type"])[METRICS].mean()
            .groupby(["model", "coarse_selection_type"]).mean().groupby("model").mean())
