"""Controller dry run: raw NTC control extraction only (no model, no prediction) for one party and one line.

Usage: round16/.venv/bin/python controller/checks/r16_dryrun_controls.py <main|at3> <line> <out_dir>
Writes <out_dir>/<party>__<line>.npz (barcodes, batch labels, gene names, raw counts, resample) or <party>__<line>.err.
"""
import sys
import time
import traceback
from pathlib import Path

import numpy as np

P = Path(__file__).resolve().parents[2]
NADIG = Path("/Users/mac/.cache/r9_p2_raw/nadig2025")
party, line, out = sys.argv[1], sys.argv[2], Path(sys.argv[3])
out.mkdir(parents=True, exist_ok=True)
t0 = time.time()
try:
    if party == "main":
        sys.path.insert(0, str(P / "round16" / "src"))
        import run_eval_r16 as M
        b = M.extract_local_ntc_bundle_r16(line, include_resample=True, raw_nadig_dir=NADIG)
        rec = dict(barcodes=b["barcodes"], batches=b["batch_labels"], genes=b["gene_names"], raw=b["raw_counts"],
                   r_barcodes=b.get("resample_barcodes"), r_raw=b.get("resample_raw_counts"))
    else:
        sys.path.insert(0, str(P / "round16" / "audit"))
        import audit_r16_raw_controls as A
        pri, res = A.extract_raw_control_pool_r16(line, custom_nadig_dir=NADIG, return_resample=True)
        rec = dict(barcodes=list(pri.cell_barcodes), batches=list(pri.batch_labels), genes=list(pri.gene_names), raw=pri.raw_counts,
                   r_barcodes=None if res is None else list(res.cell_barcodes), r_raw=None if res is None else res.raw_counts)
    arr = {k: np.asarray(v, dtype=str) for k, v in rec.items() if k in ("barcodes", "batches", "genes", "r_barcodes") and v is not None}
    arr.update({k: np.asarray(v, dtype=np.float32) for k, v in rec.items() if k in ("raw", "r_raw") and v is not None})
    np.savez(out / f"{party}__{line}.npz", **arr)
    print(f"OK {party} {line} n={len(rec['barcodes'])} G={len(rec['genes'])} resample={'none' if rec['r_barcodes'] is None else len(rec['r_barcodes'])} {time.time() - t0:.0f}s")
except Exception:
    (out / f"{party}__{line}.err").write_text(traceback.format_exc())
    print(f"ERR {party} {line} {time.time() - t0:.0f}s: {traceback.format_exc().strip().splitlines()[-1]}")
