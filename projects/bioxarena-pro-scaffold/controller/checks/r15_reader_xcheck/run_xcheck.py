"""Controller cross-check: run A-T2 (round15/src) and A-T3 (round15/audit) readers on controller fixtures,
then compare every array in their line directories exactly (no tolerance).
Usage: python3 -I run_xcheck.py <project_root> <fixtures_dir> <out_dir>
"""
import json
import sys
import traceback
from pathlib import Path

import numpy as np

ROOT = Path(sys.argv[1]).resolve()
FIX = Path(sys.argv[2]).resolve()
OUT = Path(sys.argv[3]).resolve()
sys.path.insert(0, str(ROOT))

from round15.src import dataset_readers_r15 as at2  # noqa: E402
from round15.audit import independent_r15_loader as at3  # noqa: E402

CASES = {
    "case1_mixed_explicit": {"lines": ["cal33", "hela"], "ctx": None, "man": None, "expect": "OK"},
    "case2_context_only": {"lines": ["cal33", "hela"], "ctx": "context_map.json", "man": None, "expect": "OK"},
    "case2_context_nomap": {"dir": "case2_context_only", "lines": ["cal33"], "ctx": None, "man": None, "expect": "EXCLUDED"},
    "case3_per_line_files": {"lines": ["cal33", "hela"], "ctx": None, "man": "line_manifest.json", "expect": "OK"},
    "case3_no_manifest": {"dir": "case3_per_line_files", "lines": ["cal33"], "ctx": None, "man": None, "expect": "EXCLUDED"},
    "case4_nonint": {"lines": ["cal33"], "ctx": None, "man": None, "expect": "EXCLUDED"},
}


def h5_files(d):
    return sorted(p for p in d.glob("*.h5ad"))


def run_at2(case, d, line, ctx, man):
    out = OUT / "at2" / case
    out.mkdir(parents=True, exist_ok=True)
    kw = {}
    if ctx:
        kw["context_map"] = json.loads((d / ctx).read_text())
    if man:
        kw["line_manifest"] = json.loads((d / man).read_text())
    od = at2.prepare_r15_line(line=line, raw_dir=d, out_root=out, workers=2, **kw)
    return Path(od)


def run_at3(case, d, line, ctx, man):
    out = OUT / "at3" / case
    out.mkdir(parents=True, exist_ok=True)
    kw = {}
    if ctx:
        kw["context_to_line_map"] = json.loads((d / ctx).read_text())
    if man:
        kw["line_manifest"] = json.loads((d / man).read_text())
    at3.standardize_vcc_line_to_audit_line(line_name=line, h5_paths=h5_files(d), cache_out_root=out, n_workers=2, **kw)
    cands = [p for p in out.rglob("pseudobulk_delta.npz") if line in str(p)]
    if len(cands) != 1:
        raise RuntimeError(f"cannot locate A-T3 output for {line}: {cands}")
    return cands[0].parent


def exact_compare(a, b):
    res = {}
    for name in ("genes.txt", "targets.txt"):
        ta, tb = (a / name).read_text().split(), (b / name).read_text().split()
        res[name] = {"equal_order": ta == tb, "equal_set": set(ta) == set(tb), "n": [len(ta), len(tb)]}
    for npz in ("pseudobulk_delta.npz", "ctrl_cells.npz"):
        with np.load(a / npz, allow_pickle=True) as fa, np.load(b / npz, allow_pickle=True) as fb:
            ka, kb = set(fa.files), set(fb.files)
            res[npz + ":keys_only_in_at2"] = sorted(ka - kb)
            res[npz + ":keys_only_in_at3"] = sorted(kb - ka)
            for k in sorted(ka & kb):
                xa, xb = fa[k], fb[k]
                if xa.shape != xb.shape:
                    res[f"{npz}:{k}"] = {"pass": False, "shape": [list(xa.shape), list(xb.shape)]}
                    continue
                if xa.dtype.kind in "fc":
                    same_nan = bool(np.array_equal(np.isnan(xa), np.isnan(xb)))
                    fin = np.isfinite(xa) & np.isfinite(xb)
                    d = float(np.max(np.abs(xa[fin] - xb[fin]))) if fin.any() else 0.0
                    res[f"{npz}:{k}"] = {"pass": same_nan and d == 0.0, "max_abs_diff": d, "same_nan": same_nan}
                else:
                    res[f"{npz}:{k}"] = {"pass": bool(np.array_equal(xa, xb))}
    meta = {}
    for p, tag in ((a, "at2"), (b, "at3")):
        m = json.loads((p / "line_metadata.json").read_text()) if (p / "line_metadata.json").exists() else {}
        meta[tag] = {k: m.get(k) for k in ("n_cells_total", "n_ntc", "n_ntc_cells", "n_targets", "n_genes", "split_tier", "batches") if k in m}
    res["metadata_subset"] = meta
    arr_ok = all(v.get("pass", True) for k, v in res.items() if isinstance(v, dict) and "pass" in v)
    txt_ok = res["genes.txt"]["equal_order"] and res["targets.txt"]["equal_order"]
    res["ALL_EXACT"] = bool(arr_ok and txt_ok)
    return res


def main():
    report = {}
    for case, spec in CASES.items():
        d = FIX / spec.get("dir", case)
        for line in spec["lines"]:
            key = f"{case}/{line}"
            r = {"expect": spec["expect"]}
            outs = {}
            for tag, fn in (("at2", run_at2), ("at3", run_at3)):
                try:
                    outs[tag] = fn(case, d, line, spec["ctx"], spec["man"])
                    r[tag] = "OK"
                except Exception as e:  # record, never hide
                    r[tag] = "RAISED: " + (str(e).splitlines() or [repr(e)])[0][:300]
                    r[tag + "_trace_tail"] = traceback.format_exc().splitlines()[-3:]
            if "at2" in outs and "at3" in outs:
                r["compare"] = exact_compare(outs["at2"], outs["at3"])
            r["verdict_match_expectation"] = (
                (spec["expect"] == "OK" and r.get("compare", {}).get("ALL_EXACT", False))
                or (spec["expect"] == "EXCLUDED" and r["at2"].startswith("RAISED") and r["at3"].startswith("RAISED"))
            )
            report[key] = r
            print(key, "at2:", r["at2"][:90], "| at3:", r["at3"][:90], "| exact:", r.get("compare", {}).get("ALL_EXACT"), "| ok:", r["verdict_match_expectation"], flush=True)
    (OUT / "XCHECK_REPORT.json").write_text(json.dumps(report, indent=1, default=str))


if __name__ == "__main__":
    main()
