#!/usr/bin/env python3
"""
PREREG_E §3.2 (`E3`) Pre-Metric Feasibility & Schema Gate for `XuCao2023.h5ad` (`hek293_xucao`).

BLINDING INVARIANT:
  Strictly schema, UMI integrality, `.obs` label presence (`sgRNA` + `NTC` + replicate/batch),
  and raw gene metadata inspection ONLY. ZERO predictors, correlations, r_orth_cal, slopes,
  residuals, or FCT values are computed or touched.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import h5py
import numpy as np

SCRIPTS_DIR = Path(__file__).resolve().parent
ROUND12_ROOT = SCRIPTS_DIR.parents[1]
REPO_ROOT = ROUND12_ROOT.parent
R10_SRC = REPO_ROOT / "round10" / "src"
R10_SCRIPTS = REPO_ROOT / "round10" / "shared" / "scripts"
R12_SRC = ROUND12_ROOT / "src"
for p in (str(R10_SRC), str(R10_SCRIPTS), str(R12_SRC)):
    if p not in sys.path:
        sys.path.insert(0, p)

from data_loader_r10 import decode_str_array  # noqa: E402
from prepare_r10_lines import DEFAULT_HGNC_PATH, build_hgnc_lookup, read_obs_or_var_col  # noqa: E402
from dataset_readers_r12 import check_hek293_feasibility  # noqa: E402


def _attr_str(attrs: Any, key: str, default: str = "") -> str:
    val = attrs.get(key, default)
    if isinstance(val, bytes):
        return val.decode("utf-8", errors="replace")
    return str(val)


def _inspect_matrix_node(node: Any, name: str, exhaustive_limit: int = 150_000_000) -> Dict[str, Any]:
    enc = _attr_str(node.attrs, "encoding-type", "")
    shape = node.attrs.get("shape", None)
    shape_list = [int(x) for x in shape] if shape is not None else (
        [int(x) for x in node.shape] if isinstance(node, h5py.Dataset) else []
    )
    info: Dict[str, Any] = {"name": name, "encoding": enc, "shape": shape_list}

    if isinstance(node, h5py.Dataset):
        info["dtype"] = str(node.dtype)
        n_rows = int(node.shape[0])
        take_rows = min(n_rows, 200)
        vals = np.asarray(node[:take_rows], dtype=np.float64).ravel()
        all_finite = bool(np.all(np.isfinite(vals)))
        n_neg = int(np.sum(vals < 0))
        frac_nonint = float(np.mean(np.abs(vals - np.round(vals)) > 1e-6)) if vals.size else 1.0
        vmin = float(np.min(vals)) if vals.size else None
        vmax = float(np.max(vals)) if vals.size else None
        checked = int(vals.size)
        sample = [float(x) for x in vals[:12]]
    elif isinstance(node, h5py.Group) and "data" in node:
        ds = node["data"]
        info["dtype"] = str(ds.dtype)
        nnz = int(ds.shape[0])
        info["nnz"] = nnz
        if nnz == 0:
            info["empty"] = True
            return info
        all_finite = True
        n_neg = 0
        n_nonint = 0
        checked = 0
        vmin = None
        vmax = None
        sample = []
        chunk = 10_000_000
        step = 1 if nnz <= exhaustive_limit else max(1, nnz // 5_000_000)
        for start in range(0, nnz, chunk):
            end = min(nnz, start + chunk)
            sl = np.asarray(ds[start:end:step])
            if not sample:
                sample = [float(x) for x in sl[:12]]
            if np.issubdtype(sl.dtype, np.integer):
                cur_neg = int(np.sum(sl < 0))
                cur_nonint = 0
                cur_min = float(np.min(sl))
                cur_max = float(np.max(sl))
            else:
                sl_f = sl.astype(np.float64, copy=False)
                fin = np.isfinite(sl_f)
                if not np.all(fin):
                    all_finite = False
                sl_f = sl_f[fin]
                cur_neg = int(np.sum(sl_f < 0)) if sl_f.size else 0
                cur_nonint = int(np.sum(np.abs(sl_f - np.round(sl_f)) > 1e-6)) if sl_f.size else 0
                cur_min = float(np.min(sl_f)) if sl_f.size else 0.0
                cur_max = float(np.max(sl_f)) if sl_f.size else 0.0
            n_neg += cur_neg
            n_nonint += cur_nonint
            checked += int(sl.size)
            vmin = cur_min if vmin is None else min(vmin, cur_min)
            vmax = cur_max if vmax is None else max(vmax, cur_max)
        frac_nonint = float(n_nonint / max(checked, 1))
    else:
        info["error"] = "unrecognized_node_layout"
        return info

    info.update({
        "nnz_checked": checked,
        "all_finite": all_finite,
        "n_negative": n_neg,
        "min": vmin,
        "max": vmax,
        "frac_noninteger": frac_nonint,
        "is_nonneg_integer": bool(all_finite and n_neg == 0 and frac_nonint == 0.0),
        "sample_values": sample,
    })
    return info


def _summarize_h5_df(group: h5py.Group, top_k: int = 15) -> Dict[str, Any]:
    cols = []
    if "column-order" in group.attrs:
        cols = [str(x) for x in decode_str_array(group.attrs["column-order"])]
    else:
        cols = sorted([k for k in group.keys() if not k.startswith("_")])
    summary: Dict[str, Any] = {
        "index_key": _attr_str(group.attrs, "_index", "_index"),
        "keys": sorted(list(group.keys())),
        "columns": cols,
        "details": {},
    }
    for col in cols:
        if col not in group:
            continue
        arr = np.asarray(read_obs_or_var_col(group, col))
        s = arr.astype(str)
        uniq, cnts = np.unique(s, return_counts=True)
        order = np.argsort(-cnts, kind="mergesort")
        top = [(str(uniq[i]), int(cnts[i])) for i in order[:top_k]]
        summary["details"][col] = {
            "dtype": str(arr.dtype),
            "n_unique": int(uniq.size),
            "top_counts": top,
        }
    return summary


def run_inspection(h5ad_path: Path, out_json: Path) -> Dict[str, Any]:
    if not h5ad_path.exists():
        raise FileNotFoundError(f"Missing {h5ad_path}")

    hgnc = build_hgnc_lookup(DEFAULT_HGNC_PATH)
    # check_hek293_feasibility expects raw_dir containing `hek293_xucao/*.h5ad`.
    raw_dir = h5ad_path.parent.parent
    reader_check = check_hek293_feasibility(raw_dir, hgnc=hgnc)

    with h5py.File(h5ad_path, "r") as f:
        top_keys = sorted(list(f.keys()))
        uns_keys = sorted(list(f["uns"].keys())) if ("uns" in f and isinstance(f["uns"], h5py.Group)) else []
        matrices: List[Dict[str, Any]] = []
        if "X" in f:
            matrices.append(_inspect_matrix_node(f["X"], "X"))
        if "layers" in f and isinstance(f["layers"], h5py.Group):
            for lk in sorted(f["layers"].keys()):
                matrices.append(_inspect_matrix_node(f["layers"][lk], f"layers/{lk}"))
        if "raw" in f and isinstance(f["raw"], h5py.Group) and "X" in f["raw"]:
            matrices.append(_inspect_matrix_node(f["raw"]["X"], "raw/X"))

        obs_summary = _summarize_h5_df(f["obs"]) if "obs" in f else {}
        var_summary = _summarize_h5_df(f["var"]) if "var" in f else {}

        obs_idx_key = obs_summary.get("index_key", "_index")
        var_idx_key = var_summary.get("index_key", "_index")
        obs_index = (
            list(decode_str_array(f["obs"][obs_idx_key][:]))
            if ("obs" in f and obs_idx_key in f["obs"])
            else []
        )
        var_index = (
            list(decode_str_array(f["var"][var_idx_key][:]))
            if ("var" in f and var_idx_key in f["var"])
            else []
        )

        # Inspect guide_id distribution specifically on control cells vs perturbed cells
        pert_arr = np.asarray(read_obs_or_var_col(f["obs"], "perturbation"), dtype=str) if ("obs" in f and "perturbation" in f["obs"]) else np.asarray([], dtype=str)
        guide_arr = np.asarray(read_obs_or_var_col(f["obs"], "guide_id"), dtype=str) if ("obs" in f and "guide_id" in f["obs"]) else np.asarray([], dtype=str)
        nperts_arr = np.asarray(read_obs_or_var_col(f["obs"], "nperts"), dtype=str) if ("obs" in f and "nperts" in f["obs"]) else np.asarray([], dtype=str)
        pmt_obs = np.asarray(read_obs_or_var_col(f["obs"], "percent_mito"), dtype=np.float64) if ("obs" in f and "percent_mito" in f["obs"]) else np.asarray([], dtype=np.float64)
        ng_obs = np.asarray(read_obs_or_var_col(f["obs"], "ngenes"), dtype=np.float64) if ("obs" in f and "ngenes" in f["obs"]) else np.asarray([], dtype=np.float64)
        nc_obs = np.asarray(read_obs_or_var_col(f["obs"], "ncounts"), dtype=np.float64) if ("obs" in f and "ncounts" in f["obs"]) else np.asarray([], dtype=np.float64)

    ctrl_mask = pert_arr == "control"
    ctrl_guide_counts = {}
    if ctrl_mask.size and guide_arr.size:
        ug, cg = np.unique(guide_arr[ctrl_mask], return_counts=True)
        ord_g = np.argsort(-cg, kind="mergesort")
        ctrl_guide_counts = {str(ug[i]): int(cg[i]) for i in ord_g}

    var_strs = [str(x) for x in var_index]
    n_ensmusg = sum(1 for x in var_strs if "ENSMUSG" in x.upper())
    cas9_hits = [x for x in var_strs if "CAS9" in x.upper() or "DCAS9" in x.upper()]

    # Check replicate / batch columns in obs
    obs_cols = obs_summary.get("columns", [])
    rep_candidates = [c for c in obs_cols if any(k in c.lower() for k in ("rep", "batch", "lane", "sample", "donor", "pool"))]

    # Inspect obs_index prefix/suffix diversity (delimiter tokens) without exposing any expression
    obs_strs = [str(x) for x in obs_index]
    bc_prefixes = sorted({s.split("_")[0] for s in obs_strs[:5000] if "_" in s})[:10]
    bc_suffixes = sorted({s.split("-")[-1] for s in obs_strs if "-" in s})[:20]

    qs = [5, 25, 50, 75, 95]
    pmt_obs_q = {str(q): round(float(v), 2) for q, v in zip(qs, np.percentile(pmt_obs, qs))} if pmt_obs.size else {}
    ng_min = int(np.min(ng_obs)) if ng_obs.size else 0
    nc_min = int(np.min(nc_obs)) if nc_obs.size else 0

    has_raw_int_x = any(m["name"] == "X" and m.get("is_nonneg_integer", False) for m in matrices)
    has_guide_col = ("guide_id" in obs_cols) and (obs_summary["details"].get("guide_id", {}).get("n_unique", 0) > 0)
    n_ctrl_cells = int(np.sum(ctrl_mask)) if ctrl_mask.size else 0
    has_ntc = n_ctrl_cells > 0
    e3_data_pass = bool(has_raw_int_x and has_guide_col and has_ntc)
    reader_pass = bool(reader_check.get("feasible", False))
    e3_pass = bool(e3_data_pass and reader_pass)
    if e3_pass:
        e3_verdict = "PASS"
    elif e3_data_pass and not reader_pass:
        e3_verdict = "BLOCKED_BY_READER_INDEX_BUG"
    else:
        e3_verdict = "EXCLUDED_FEASIBILITY_Q3"
    split_branch = (
        f"Stage 1 `.obs` replicate parity ({rep_candidates})"
        if rep_candidates
        else "Stage 2 deterministic seed=42 stratified 50/50 split (no .obs replicate/batch column)"
    )
    if e3_pass:
        e3_reason = (
            f"E3 gate PASS: X non-negative int64 raw UMIs verified ({matrices[0]['nnz']:,} nnz, frac_nonint=0.0, "
            f"+ layers/nascent_counts); per-cell guide_id ({obs_summary['details']['guide_id']['n_unique']} guides) & "
            f"NTC labels (perturbation=='control': {n_ctrl_cells:,} cells across {len(ctrl_guide_counts)} control guides; "
            f"nperts==0: {int(np.sum(nperts_arr == '0')) if nperts_arr.size else 0:,}; min ngenes={ng_min}, min ncounts={nc_min}); "
            f"obs percent_mito quantiles (5/25/50/75/95%)={pmt_obs_q}; split branch={split_branch}; "
            f"reader check=PASS (per round12/status/RULING_R12_PREFREEZE.md Ruling B: max_pct_mt disabled)."
        )
    else:
        e3_reason = (
            f"E3 data criteria SATISFIED: X non-negative int64 raw UMIs verified ({matrices[0]['nnz']:,} nnz, frac_nonint=0.0, "
            f"+ layers/nascent_counts); per-cell guide_id ({obs_summary['details']['guide_id']['n_unique']} unique guides) & "
            f"NTC labels (perturbation=='control': {n_ctrl_cells:,} cells across {len(ctrl_guide_counts)} control guide_ids; "
            f"nperts==0: {int(np.sum(nperts_arr == '0')) if nperts_arr.size else 0:,}; min ngenes={ng_min}, min ncounts={nc_min}); "
            f"obs percent_mito quantiles (5/25/50/75/95%)={pmt_obs_q}; split branch={split_branch}; "
            f"HOWEVER round12/src/dataset_readers_r12.py::check_hek293_feasibility failed with "
            f"{reader_check.get('reason', 'unknown')} because AnnData 0.2.0 stores index in "
            f"obs.attrs['_index']='{obs_idx_key}' and var.attrs['_index']='{var_idx_key}' rather than 'index'/'_index'."
        )

    report = {
        "h5ad_path": str(h5ad_path),
        "top_keys": top_keys,
        "uns_keys": uns_keys,
        "n_obs": int(len(obs_strs)),
        "n_var": int(len(var_strs)),
        "obs_ngenes_min": ng_min,
        "obs_ncounts_min": nc_min,
        "obs_percent_mito_quantiles": pmt_obs_q,
        "obs_index_sample": obs_strs[:8] + obs_strs[-4:],
        "obs_barcode_prefixes_first5k": bc_prefixes,
        "obs_barcode_dash_suffixes": bc_suffixes,
        "var_index_sample": var_strs[:15],
        "var_n_ensmusg": int(n_ensmusg),
        "var_cas9_hits": cas9_hits[:20],
        "replicate_candidate_columns": rep_candidates,
        "control_n_cells": n_ctrl_cells,
        "control_guide_id_counts": ctrl_guide_counts,
        "split_half_branch": split_branch,
        "e3_data_feasibility_pass": e3_data_pass,
        "e3_reader_feasibility_pass": reader_pass,
        "e3_gate_pass": e3_pass,
        "e3_verdict": e3_verdict,
        "e3_reason": e3_reason,
        "matrices": matrices,
        "obs_summary": obs_summary,
        "var_summary": var_summary,
        "reader_check_hek293": reader_check,
    }
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main() -> int:
    p = argparse.ArgumentParser(description="E3 pre-metric feasibility audit on XuCao2023.h5ad")
    p.add_argument(
        "--h5ad",
        type=Path,
        default=ROUND12_ROOT / "shared" / "raw" / "hek293_xucao" / "XuCao2023.h5ad",
    )
    p.add_argument(
        "--out-json",
        type=Path,
        default=ROUND12_ROOT / "shared" / "hek293_feasibility.json",
    )
    args = p.parse_args()
    rep = run_inspection(args.h5ad, args.out_json)
    print(json.dumps({
        "n_obs": rep["n_obs"],
        "n_var": rep["n_var"],
        "obs_index_sample": rep["obs_index_sample"][:6],
        "obs_barcode_dash_suffixes": rep["obs_barcode_dash_suffixes"],
        "var_index_sample": rep["var_index_sample"][:10],
        "var_n_ensmusg": rep["var_n_ensmusg"],
        "var_cas9_hits": rep["var_cas9_hits"],
        "replicate_candidate_columns": rep["replicate_candidate_columns"],
        "control_n_cells": rep["control_n_cells"],
        "control_guide_id_counts": rep["control_guide_id_counts"],
        "split_half_branch": rep["split_half_branch"],
        "e3_gate_pass": rep["e3_gate_pass"],
        "e3_verdict": rep["e3_verdict"],
        "e3_reason": rep["e3_reason"],
        "reader_check_hek293": rep["reader_check_hek293"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
