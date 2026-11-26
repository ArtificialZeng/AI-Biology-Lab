"""Controller: staged verification step 2 (Addendum A §11 item 5): A-T3 arrays vs main-line arrays, one line at a time.

A-T3 file: round16/audit/results/{model}_{ckpt}__{line}.npz (one per model and target line).
Main dir:  round16/results/arrays/<line>/ (observed.npz, preds_<arm>_{matched,mis_<s>,resample}.npz).

Checks:
- Target and gene order: exact string equality; target_indices exact.
- Observed-side arrays (delta_obs, halves, half_valid, delta_transfer, source_grand_mean, FCT): |diff| <= 1e-5.
- Predicted Δ of the primary arm (STATE_SE fixed_index_zero; STATE_HVG official_median_log1p): matched, every
  mismatched source and the within-line resample, |diff| <= 1e-4 (Addendum A tolerance for predicted Δ).
Usage: python3 -I r16_compare_arrays.py <at3_npz> <main_line_dir>
"""
import json
import sys
from pathlib import Path

import numpy as np

TOL_OBS, TOL_PRED = 1e-5, 1e-4
A = np.load(sys.argv[1], allow_pickle=False)
D = Path(sys.argv[2])
M = np.load(D / "observed.npz", allow_pickle=False)
model = str(A["model_key"])
prefix = {"STATE_SE": "preds_STATE_SE_fixed_index_zero_", "STATE_HVG": "preds_STATE_HVG_official_median_log1p_"}[model]
rows, bad = [], []


def same(tag, a, b):
    ok = [str(x) for x in a] == [str(x) for x in b]
    rows.append((tag, "exact", ok, None))
    ok or bad.append(tag)


def num(tag, a, b, tol):
    a, b = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    if a.shape != b.shape:
        rows.append((tag, f"shape {a.shape} vs {b.shape}", False, None)); bad.append(tag); return
    fa, fb = np.isfinite(a), np.isfinite(b)
    if not np.array_equal(fa, fb):
        rows.append((tag, "finite-mask differs", False, None)); bad.append(tag); return
    d = float(np.max(np.abs(a[fa] - b[fb]))) if fa.any() else 0.0
    ok = d <= tol
    rows.append((tag, f"max|diff|={d:.3g} (tol {tol:g}, max|x|={float(np.max(np.abs(b[fb]))) if fb.any() else 0:.3g})", ok, d))
    ok or bad.append(tag)


same("target_names", A["eval_perts"], M["target_names"])
same("gene_names", A["eval_genes"], M["gene_names"])
same("target_indices", A["target_gene_idx"], M["target_indices"])
for ka, km in [("delta_obs", "delta_obs"), ("delta_half_a", "delta_half_a"), ("delta_half_b", "delta_half_b"),
               ("half_valid", "half_valid"), ("delta_transfer", "delta_transfer"),
               ("source_grand_mean", "source_grand_mean"), ("delta_fct", "fct_delta_pred")]:
    num(ka, A[ka], M[km], TOL_OBS)
num(f"{model}.matched", A["delta_pred_matched"], np.load(D / f"{prefix}matched.npz")["delta_pred"], TOL_PRED)
for s in [str(x) for x in A["seen_sources"]] + [str(x) for x in A["unseen_sources"]]:
    num(f"{model}.mis_{s}", A[f"delta_pred_mismatched__{s}"], np.load(D / f"{prefix}mis_{s}.npz")["delta_pred"], TOL_PRED)
if "delta_pred_resample" in A.files:
    num(f"{model}.resample", A["delta_pred_resample"], np.load(D / f"{prefix}resample.npz")["delta_pred"], TOL_PRED)
else:
    rows.append((f"{model}.resample", "absent in A-T3 file", None, None))

line = str(A["target_line"])
for tag, info, ok, _ in rows:
    print(f"{line} {tag}: {'OK' if ok else ('n/a' if ok is None else 'FAIL')} {info}")
worst_pred = max([d for t, _, ok, d in rows if d is not None and t.startswith(model + ".")] or [0.0])
print(json.dumps({"line": line, "model": model, "n_checks": len(rows), "n_fail": len(bad), "worst_pred_diff": worst_pred}))
print(f"R16_ARRAYS_MATCH {line} {model}" if not bad else f"R16_ARRAYS_DIFF {line} {model} ({', '.join(bad)})")
