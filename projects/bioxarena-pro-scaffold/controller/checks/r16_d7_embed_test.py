"""Controller diagnosis for DEVIATION_R16_D7 (SE-600M gene mapping on duplicate raw symbols).

Addendum A §11 item 3 locks: adata.var_names = raw symbols, then adata.var_names_make_unique() (first occurrence
kept, so duplicates never map to the same SE-600M global ID). arc-state's FilteredGenesCounts maps by var_names.
Main line calls var_names_make_unique(); A-T3's encode_with_loaded_se_model does not.

This test reruns A-T3's own frozen encoder (CleanRoomSE600MEncoder, MPS, seed 20261601/20261602) on the controller
dry-run control matrices (identical for both parties) and checks:
  (a) raw symbols (A-T3 as frozen)        -> reproduces A-T3's exported X_state
  (b) symbols made unique by anndata first -> reproduces the main line's exported X_state (Addendum A tolerance 1e-5)
Usage: round16/.venv/bin/python controller/checks/r16_d7_embed_test.py <dryrun_dir> [lines...]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from anndata.utils import make_index_unique

P = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(P / "round16" / "audit"))
from audit_r16_state_adapter import CleanRoomSE600MEncoder  # noqa: E402

dd = Path(sys.argv[1])
lines = sys.argv[2:] or ["nadig_hepg2", "nadig_jurkat", "dld1", "feng_ipsc"]
AT3 = np.load(P / "round16/audit/results/se600m_embeddings_r16.npz", allow_pickle=False)
enc = CleanRoomSE600MEncoder(se600m_dir=P / "round16/shared/models/SE-600M", device="mps", batch_size=16)
enc.load()
pe = set(enc._pe_keys)
ok_all = True
for ln in lines:
    z = np.load(dd / f"main__{ln}.npz", allow_pickle=False)
    genes = [str(g) for g in z["genes"]]
    s = pd.Series(genes)
    dups = sorted(set(s[s.duplicated()]))
    uniq = list(make_index_unique(pd.Index(genes)))
    print(f"{ln}: {len(dups)} duplicated symbols, {sum(g in pe for g in dups)} of them in SE-600M PE: {[g for g in dups if g in pe]}")
    for tag, raw, seed, mkey, akey in [("primary", z["raw"], 20261601, f"x_state_{ln}.npz", f"{ln}__X_state"),
                                       ("resample", z["r_raw"] if "r_raw" in z.files else None, 20261602, f"x_state_{ln}_resample.npz", f"{ln}__resample__X_state")]:
        mf = P / "round16/results/embeddings" / mkey
        if raw is None or not mf.exists() or akey not in AT3.files:
            continue
        M = np.load(mf, allow_pickle=False)["X_state"].astype(np.float64)
        A = AT3[akey].astype(np.float64)
        xa, _ = enc.encode_cells(ctrl_counts=raw, source_genes=genes, seed=seed, dataset_name=ln)
        xu, _ = enc.encode_cells(ctrl_counts=raw, source_genes=uniq, seed=seed, dataset_name=ln)
        d_a_at3 = float(np.max(np.abs(xa - A))); d_u_main = float(np.max(np.abs(xu - M))); d_a_main = float(np.max(np.abs(xa - M)))
        ok = d_a_at3 <= 1e-5 and d_u_main <= 1e-5
        ok_all &= ok
        print(f"  {tag}: raw-symbols vs A-T3 export {d_a_at3:.3g} | raw-symbols vs main {d_a_main:.3g} | unique-symbols vs main export {d_u_main:.3g} -> {'CONFIRMED' if ok else 'NOT CONFIRMED'}")
print("R16_D7_DIAGNOSIS_CONFIRMED" if ok_all else "R16_D7_DIAGNOSIS_NOT_CONFIRMED")
