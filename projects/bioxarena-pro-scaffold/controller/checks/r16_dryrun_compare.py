"""Controller: compare main-line vs A-T3 dry-run control extraction (r16_dryrun_controls.py outputs).
Addendum A tolerances: barcodes exact (ordered); raw counts must match exactly (integers) after name alignment.
Usage: python3 -I r16_dryrun_compare.py <dryrun_dir> [lines...]
"""
import sys
from pathlib import Path

import numpy as np

d = Path(sys.argv[1])
lines = sys.argv[2:] or ["a549_ko", "cd4t_rest", "dld1", "feng_ipsc", "nadig_hepg2", "nadig_jurkat"]
allok = True
for ln in lines:
    fm, fa = d / f"main__{ln}.npz", d / f"at3__{ln}.npz"
    if not (fm.exists() and fa.exists()):
        print(f"{ln}: MISSING main={fm.exists()} at3={fa.exists()}"); allok = False; continue
    m, a = np.load(fm), np.load(fa)
    r = {}
    r["barcodes_equal"] = m["barcodes"].shape == a["barcodes"].shape and bool(np.all(m["barcodes"] == a["barcodes"]))
    r["batches_equal"] = m["batches"].shape == a["batches"].shape and bool(np.all(m["batches"] == a["batches"]))
    gm, ga = list(m["genes"]), list(a["genes"])
    r["gene_order_equal"] = gm == ga
    common = sorted(set(gm) & set(ga)); r["genes"] = (len(gm), len(ga), len(common))
    if r["barcodes_equal"]:
        im = {g: i for i, g in enumerate(gm)}; ia = {g: i for i, g in enumerate(ga)}
        dup = len(im) != len(gm) or len(ia) != len(ga)
        X = m["raw"][:, [im[g] for g in common]]; Y = a["raw"][:, [ia[g] for g in common]]
        r["raw_max_abs_diff_common"] = float(np.max(np.abs(X - Y))) if X.size else None
        r["dup_gene_names"] = dup
    has_rm, has_ra = "r_barcodes" in m.files, "r_barcodes" in a.files
    r["resample"] = (has_rm, has_ra, bool(has_rm and has_ra and m["r_barcodes"].shape == a["r_barcodes"].shape and np.all(m["r_barcodes"] == a["r_barcodes"])))
    ok = r["barcodes_equal"] and r["batches_equal"] and r.get("raw_max_abs_diff_common") == 0.0 and (r["resample"][0] == r["resample"][1]) and (not r["resample"][0] or r["resample"][2])
    allok &= bool(ok)
    print(f"{ln}: {'MATCH' if ok else 'DIFF'} {r}")
print("DRYRUN_ALL_MATCH" if allok else "DRYRUN_HAS_DIFFERENCES")
