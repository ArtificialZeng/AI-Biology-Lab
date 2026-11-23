#!/usr/bin/env python3
"""
Clean-Room Round 14 Raw Rebuild Cross-Check Runner (A-T3).

Strictly per `主控 10:33` & `PREREG_G.md`:
  - Runs independent `independent_r14_loader` on raw `GSE337804` `.h5ad` files into
    a temporary folder under `round14/audit/_xcheck_tmp`.
  - Never modifies `round14/shared/lines/`.
  - Never imports `round14/src/` or `round14/results/`.
  - Never computes any predictions, correlations, or statistical metrics; compares ONLY
    counts, slot selection, gene/target sets, and numeric array differences (`<= 1e-6`).
"""

from __future__ import annotations

import json
import shutil
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Set

import numpy as np

AUDIT_R14_DIR = Path(__file__).resolve().parent
REPO_ROOT = AUDIT_R14_DIR.parents[1]
R10_SCRIPTS_DIR = REPO_ROOT / "round10" / "shared" / "scripts"
R11_AUDIT_DIR = REPO_ROOT / "round11" / "audit"
for _p in (str(AUDIT_R14_DIR), str(R10_SCRIPTS_DIR), str(R11_AUDIT_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from independent_r11_loader import (  # noqa: E402
    DEFAULT_HGNC_TSV,
    R9_LINE_NAMES,
    R9_LINES_DIR,
    CanonicalHGNC,
    load_line_dataset,
    load_r9_common_genes,
)
from independent_r14_loader import (  # noqa: E402
    R14_LINES_DIR,
    R14_RAW_DIR,
    R14_STATUS_DIR,
    R14_TARGET_KEY,
    build_a549_myllia_audit_parts,
    enforce_r14_gate,
)
from prepare_r10_lines import build_hgnc_lookup, process_r10_line  # noqa: E402
from run_r14_audit import NOVEL_NOT_IN_R12_A549, OVERLAP_R12_A549  # noqa: E402


def _read_text_lines(path: Path) -> List[str]:
    return [ln.strip() for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]


def _array_diff(
    a: np.ndarray,
    b: np.ndarray,
    row_names: List[str] | None = None,
    col_names: List[str] | None = None,
) -> Dict[str, Any]:
    aa = np.asarray(a, dtype=np.float64)
    bb = np.asarray(b, dtype=np.float64)
    if aa.shape != bb.shape:
        return {"shape_match": False, "shape_rebuilt": list(aa.shape), "shape_shared": list(bb.shape), "max_abs_diff": float("inf")}
    nan_match = bool(np.array_equal(np.isnan(aa), np.isnan(bb)))
    fin = np.isfinite(aa) & np.isfinite(bb)
    abs_d = np.where(fin, np.abs(aa - bb), 0.0)
    max_abs = float(np.max(abs_d)) if np.any(fin) else 0.0
    flat_idx = int(np.argmax(abs_d)) if np.any(fin) else -1
    worst: Dict[str, Any] = {}
    if flat_idx >= 0 and max_abs > 0.0:
        coord = np.unravel_index(flat_idx, aa.shape)
        worst = {
            "coord": [int(c) for c in coord],
            "rebuilt_val": float(aa[coord]),
            "shared_val": float(bb[coord]),
        }
        if row_names is not None and len(coord) >= 1:
            worst["row_name"] = row_names[int(coord[0])]
        if col_names is not None and len(coord) >= 2:
            worst["col_name"] = col_names[int(coord[1])]
        elif col_names is not None and len(coord) == 1:
            worst["col_name"] = col_names[int(coord[0])]
    max_rel = float(np.max(abs_d[fin] / np.maximum(np.abs(bb[fin]), 1e-12))) if np.any(fin) else 0.0
    n_gt_1e6 = int(np.sum(abs_d > 1e-6)) if np.any(fin) else 0
    return {
        "shape_match": True,
        "shape": list(aa.shape),
        "nan_mask_match": nan_match,
        "max_abs_diff": max_abs,
        "max_rel_diff": max_rel,
        "n_entries_gt_1e6": n_gt_1e6,
        "worst_entry": worst,
        "within_1e6": bool(nan_match and max_abs <= 1e-6),
    }


def _compute_eligible_sets(
    line_dir: Path,
    common_genes: List[str],
    hgnc_canon: CanonicalHGNC,
    src_gene_sets: List[Set[str]],
    src_avail_targets: Set[str],
) -> Dict[str, Any]:
    """Compute aligned `|G_eval|`, reproducible target set, and evaluation target set (counts/sets ONLY)."""
    ld = load_line_dataset(
        line_dir=line_dir,
        representation="log1p",
        aligned_genes=common_genes,
        hgnc=hgnc_canon,
        min_gene_overlap=3000,
    )
    eval_genes = [g for g in ld.genes if all(g in ss for ss in src_gene_sets)]
    eval_gene_set = set(eval_genes)

    rep_targets = sorted(
        ld.target_names[i]
        for i in range(len(ld.target_names))
        if (ld.repr_mask is None or bool(ld.repr_mask[i]))
    )
    eval_targets = sorted(
        t for t in rep_targets if (t in eval_gene_set and t in src_avail_targets)
    )
    ov_set = set(OVERLAP_R12_A549)
    nov_set = set(NOVEL_NOT_IN_R12_A549)
    eval_ov = sorted(t for t in eval_targets if t in ov_set)
    eval_nov = sorted(t for t in eval_targets if t in nov_set)
    eval_other = sorted(t for t in eval_targets if t not in ov_set and t not in nov_set)
    return {
        "n_genes_aligned_r9c": len(ld.genes),
        "n_genes_eval": len(eval_genes),
        "n_targets_dedup": len(ld.target_names),
        "n_rep_targets": len(rep_targets),
        "rep_targets": rep_targets,
        "n_eval_targets": len(eval_targets),
        "eval_targets": eval_targets,
        "overlap_r12_a549": eval_ov,
        "novel_not_in_r12_a549": eval_nov,
        "other_targets": eval_other,
    }


def main() -> None:
    t0 = time.time()
    gate_meta = enforce_r14_gate(
        phase="reader_xcheck_r14",
        target_lines=(R14_TARGET_KEY,),
        status_dir=R14_STATUS_DIR,
        require_confirmatory_file=False,
    )

    raw_dir = R14_RAW_DIR / "gse337804"
    shared_dir = R14_LINES_DIR / R14_TARGET_KEY
    tmp_root = AUDIT_R14_DIR / "_xcheck_tmp"
    if tmp_root.exists():
        shutil.rmtree(tmp_root)
    tmp_root.mkdir(parents=True, exist_ok=True)

    try:
        hgnc_prep = build_hgnc_lookup(DEFAULT_HGNC_TSV)
        parts, audit_meta = build_a549_myllia_audit_parts(
            raw_dir=raw_dir,
            hgnc=hgnc_prep,
        )
        rebuilt_dir = process_r10_line(
            line_name=R14_TARGET_KEY,
            parts=parts,
            hgnc=hgnc_prep,
            chunk_size=8192,
            out_root=tmp_root,
            extra_meta={"audit_reader": audit_meta},
            n_workers=4,
        )

        # Load shared & rebuilt metadata
        shared_meta = json.loads((shared_dir / "line_metadata.json").read_text(encoding="utf-8"))
        rebuilt_meta = json.loads((rebuilt_dir / "line_metadata.json").read_text(encoding="utf-8"))

        shared_genes = _read_text_lines(shared_dir / "genes.txt")
        rebuilt_genes = _read_text_lines(rebuilt_dir / "genes.txt")
        shared_targets = _read_text_lines(shared_dir / "targets.txt")
        rebuilt_targets = _read_text_lines(rebuilt_dir / "targets.txt")

        # Load source R9_6 ONLY for G_eval gene intersection and reproducible target union (zero prediction/metric)
        hgnc_canon = CanonicalHGNC.from_tsv(DEFAULT_HGNC_TSV)
        common_genes = load_r9_common_genes(hgnc=hgnc_canon)
        src_gene_sets: List[Set[str]] = []
        src_avail_targets: Set[str] = set()
        for sn in R9_LINE_NAMES:
            sdata = load_line_dataset(
                R9_LINES_DIR / sn,
                representation="log1p",
                aligned_genes=common_genes,
                hgnc=hgnc_canon,
                min_gene_overlap=3000,
            )
            src_gene_sets.append(set(sdata.genes))
            for i, tn in enumerate(sdata.target_names):
                if sdata.repr_mask is None or bool(sdata.repr_mask[i]):
                    src_avail_targets.add(tn)

        elig_shared = _compute_eligible_sets(shared_dir, common_genes, hgnc_canon, src_gene_sets, src_avail_targets)
        elig_rebuilt = _compute_eligible_sets(rebuilt_dir, common_genes, hgnc_canon, src_gene_sets, src_avail_targets)

        # Exact array diffs aligned by gene order & target order
        with np.load(shared_dir / "ctrl_cells.npz", allow_pickle=False) as cs, np.load(
            rebuilt_dir / "ctrl_cells.npz", allow_pickle=False
        ) as cr:
            shared_ctrl_keys = sorted(cs.files)
            rebuilt_ctrl_keys = sorted(cr.files)
            if rebuilt_genes == shared_genes:
                g_perm = np.arange(len(rebuilt_genes), dtype=np.int64)
            else:
                g_map = {g: idx for idx, g in enumerate(rebuilt_genes)}
                g_perm = np.asarray([g_map[g] for g in shared_genes], dtype=np.int64)

            ctrl_diffs: Dict[str, Any] = {}
            for k in ("ctrl_mean_libnorm", "ctrl_mean_libnorm_full", "ctrl_sd_libnorm", "ctrl_mean_log1p", "X_log1p"):
                if k in cs.files and k in cr.files:
                    arr_r = cr[k][g_perm] if cr[k].ndim == 1 else cr[k][:, g_perm]
                    ctrl_diffs[k] = _array_diff(arr_r, cs[k])

        with np.load(shared_dir / "pseudobulk_delta.npz", allow_pickle=False) as ps, np.load(
            rebuilt_dir / "pseudobulk_delta.npz", allow_pickle=False
        ) as pr:
            shared_pb_keys = sorted(ps.files)
            rebuilt_pb_keys = sorted(pr.files)
            s_t_list = [str(x) for x in ps["targets"].tolist()] if "targets" in ps.files else shared_targets
            r_t_list = [str(x) for x in pr["targets"].tolist()] if "targets" in pr.files else rebuilt_targets
            if r_t_list == s_t_list:
                t_perm = np.arange(len(r_t_list), dtype=np.int64)
            else:
                t_map = {t: idx for idx, t in enumerate(r_t_list)}
                t_perm = np.asarray([t_map[t] for t in s_t_list], dtype=np.int64)

            pb_diffs: Dict[str, Any] = {}
            for k in (
                "delta_log1p",
                "half_a_delta_log1p",
                "half_b_delta_log1p",
                "knockdown_depth_log1p",
                "n_cells",
                "ad_counts",
                "energy_pvals",
            ):
                if k in ps.files and k in pr.files:
                    if pr[k].ndim == 2:
                        arr_r = pr[k][np.ix_(t_perm, g_perm)]
                        pb_diffs[k] = _array_diff(arr_r, ps[k], row_names=s_t_list, col_names=shared_genes)
                    else:
                        arr_r = pr[k][t_perm]
                        pb_diffs[k] = _array_diff(arr_r, ps[k], row_names=None, col_names=s_t_list)

            for bk in ("reproducible", "half_valid"):
                if bk in ps.files and bk in pr.files:
                    arr_r = np.asarray(pr[bk], dtype=bool)[t_perm]
                    arr_s = np.asarray(ps[bk], dtype=bool)
                    mismatch_idx = np.flatnonzero(arr_r != arr_s)
                    pb_diffs[bk] = {
                        "exact_match": bool(mismatch_idx.size == 0),
                        "n_mismatches": int(mismatch_idx.size),
                        "mismatch_targets": [s_t_list[int(i)] for i in mismatch_idx[:20]],
                    }

        out = {
            "elapsed_sec": round(time.time() - t0, 2),
            "gate_meta": gate_meta,
            "feasibility": audit_meta["feasibility"],
            "audit_meta": {k: v for k, v in audit_meta.items() if k != "feasibility"},
            "rebuilt_line_meta": rebuilt_meta,
            "shared_line_meta": shared_meta,
            "genes_exact_order_match": bool(rebuilt_genes == shared_genes),
            "genes_set_match": bool(set(rebuilt_genes) == set(shared_genes)),
            "targets_exact_order_match": bool(rebuilt_targets == shared_targets),
            "targets_set_match": bool(set(rebuilt_targets) == set(shared_targets)),
            "ctrl_keys_match": bool(rebuilt_ctrl_keys == shared_ctrl_keys),
            "pb_keys_match": bool(rebuilt_pb_keys == shared_pb_keys),
            "eligible_rebuilt": elig_rebuilt,
            "eligible_shared": elig_shared,
            "ctrl_diffs": ctrl_diffs,
            "pb_diffs": pb_diffs,
        }
        out_json = AUDIT_R14_DIR / "_xcheck_r14_summary.json"
        out_json.write_text(json.dumps(out, indent=2), encoding="utf-8")
        print(json.dumps({k: v for k, v in out.items() if k not in ("rebuilt_line_meta", "shared_line_meta")}, indent=2))
    finally:
        if tmp_root.exists():
            shutil.rmtree(tmp_root, ignore_errors=True)


if __name__ == "__main__":
    main()
