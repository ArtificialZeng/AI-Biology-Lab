"""Controller: compare main-line vs A-T3 SE-600M control embeddings (Addendum A: barcodes exact and ordered,
embeddings max |diff| <= 1e-5).
Usage: python3 -I r16_compare_embeddings.py <main_embeddings_dir> <at3_se600m_embeddings_r16.npz>
"""
import sys
from pathlib import Path

import numpy as np

TOL = 1e-5
md, A = Path(sys.argv[1]), np.load(sys.argv[2], allow_pickle=False)
ok_all = True
for ln in ["a549_ko", "cd4t_rest", "dld1", "feng_ipsc", "nadig_hepg2", "nadig_jurkat"]:
    for suf, mf in [("", f"x_state_{ln}.npz"), ("__resample", f"x_state_{ln}_resample.npz")]:
        ka = f"{ln}{suf}__X_state"
        has_m, has_a = (md / mf).exists(), ka in A.files
        if not has_m and not has_a:
            continue
        if has_m != has_a:
            print(f"{ln}{suf}: PRESENCE DIFF main={has_m} at3={has_a}"); ok_all = False; continue
        m = np.load(md / mf, allow_pickle=False)
        bc = m["control_barcodes"].shape == A[f"{ln}{suf}__barcodes"].shape and bool(np.all(m["control_barcodes"] == A[f"{ln}{suf}__barcodes"]))
        X, Y = m["X_state"].astype(np.float64), A[ka].astype(np.float64)
        d = float(np.max(np.abs(X - Y))) if X.shape == Y.shape else float("inf")
        ok = bc and d <= TOL
        ok_all &= ok
        print(f"{ln}{suf}: {'MATCH' if ok else 'DIFF'} barcodes_equal={bc} shape={X.shape}/{Y.shape} max_abs_diff={d:.3g} "
              f"max|X|={np.max(np.abs(X)):.3g}")
print("R16_EMBEDDINGS_MATCH" if ok_all else "R16_EMBEDDINGS_DIFF")
