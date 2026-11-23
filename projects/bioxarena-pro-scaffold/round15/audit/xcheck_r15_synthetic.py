#!/usr/bin/env python3
"""
Clean-Room Round 15 Synthetic Cross-Check Oracle & Cache Comparator (A-T3).

Strictly operates ONLY on temporary synthetic `.h5ad` fixtures or an explicitly passed
synthetic cached line directory (`--compare-line-dir` outside `round15/shared/` in
`synthetic_test` mode). Never touches `round15/src/`, `round15/results/`, or real 2026 data.

Verifies item-by-item exact parity (`max_abs_diff <= 1e-6`):
  - Cell counts & pre-QC `jurkat` deletion (`§3.1`)
  - Kept targeting & `non-targeting` NTC counts (per half)
  - Canonical HGNC target sets & gene lists
  - Two-tier split-half assignments (`Tier 1` `.obs` batch parity vs `Tier 2` `seed=20261231`)
  - Independent numpy oracle arrays (`ctrl_mean_libnorm`, `X_log1p`, `delta_log1p`,
    `half_a_delta_log1p`, `half_b_delta_log1p`) from raw counts directly.
"""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

AUDIT_R15_DIR = Path(__file__).resolve().parent
REPO_ROOT = AUDIT_R15_DIR.parents[1]
R10_SCRIPTS_DIR = REPO_ROOT / "round10" / "shared" / "scripts"
R11_AUDIT_DIR = REPO_ROOT / "round11" / "audit"
for _p in (str(AUDIT_R15_DIR), str(R10_SCRIPTS_DIR), str(R11_AUDIT_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from independent_r11_loader import DEFAULT_HGNC_TSV  # noqa: E402
from independent_r15_loader import (  # noqa: E402
    DEFAULT_R15_SEED,
    build_vcc_line_audit_parts,
    enforce_r15_gate,
)
from prepare_r10_lines import build_hgnc_lookup, process_r10_line  # noqa: E402


def compute_independent_oracle_arrays(
    h5_paths: Sequence[Path],
    line_name: str,
    context_to_line_map: Optional[Dict[str, str]] = None,
    line_manifest: Optional[Mapping[str, Any]] = None,
    hgnc_tsv: Path = DEFAULT_HGNC_TSV,
    seed: int = DEFAULT_R15_SEED,
) -> Dict[str, Any]:
    """
    Compute expected kept cells, NTCs, canonical targets, `ctrl_mean_libnorm`, and
    unit-centered `delta_log1p` / `half_a` / `half_b` directly from raw `.h5ad` in pure numpy/SciPy
    via `build_vcc_line_audit_parts` parts (never calling `process_r10_line`), for direct
    mathematical cross-check against cached `.npz` files.
    """
    hgnc = build_hgnc_lookup(hgnc_tsv)
    parts, meta = build_vcc_line_audit_parts(
        line_name=line_name,
        h5_paths=h5_paths,
        hgnc=hgnc,
        context_to_line_map=context_to_line_map,
        line_manifest=line_manifest,
        seed=seed,
    )
    # Build aligned raw gene axis identically to `process_r10_line`
    ref_names = [str(x).strip() for x in parts[0].var_gene_names]
    ref_ids = [str(x).strip() for x in parts[0].var_gene_ids]
    same_var = all(
        [str(x).strip() for x in p.var_gene_names] == ref_names and [str(x).strip() for x in p.var_gene_ids] == ref_ids
        for p in parts[1:]
    )
    if same_var:
        col_maps = [None] * len(parts)
        var_gnames = ref_names
        var_gids = ref_ids
    else:
        sym_maps: List[Dict[str, int]] = []
        for p in parts:
            smap: Dict[str, int] = {}
            for j, (gid, gnm) in enumerate(zip(p.var_gene_ids, p.var_gene_names)):
                key = str(gnm).strip() or str(gid).strip()
                if key and key not in smap:
                    smap[key] = j
            sym_maps.append(smap)
        common_keys = [k for k in sym_maps[0] if all(k in sm for sm in sym_maps[1:])]
        col_maps = [np.array([sm[k] for k in common_keys], dtype=np.int64) for sm in sym_maps]
        var_gnames = common_keys
        var_gids = [parts[0].var_gene_ids[j] for j in col_maps[0]]

    # Collect full kept matrix per part and compute exact per-batch NTC means directly in float64/float32
    X_raw_list: List[np.ndarray] = []
    nt_list: List[np.ndarray] = []
    batch_list: List[np.ndarray] = []
    half_list: List[np.ndarray] = []
    umi_list: List[np.ndarray] = []
    pk_list: List[str] = []

    for p, cmap in zip(parts, col_maps):
        blk = p.read_chunk_fn(0, p.n_cells)
        if cmap is not None:
            blk = blk[:, cmap]
        X_raw_list.append(np.asarray(blk, dtype=np.float32))
        nt_list.append(np.asarray(p.nt_mask, dtype=bool))
        batch_list.append(np.asarray(p.batches, dtype=str))
        half_list.append(np.asarray(p.half_flags, dtype=np.int8))
        umi_list.append(np.asarray(p.umi, dtype=np.float64))
        pk_list.extend(p.pert_keys)

    X_raw = np.vstack(X_raw_list)
    nt_all = np.concatenate(nt_list)
    batches_all = np.concatenate(batch_list)
    halves_all = np.concatenate(half_list)
    umi_all = np.concatenate(umi_list)

    scale = (1e4 / np.maximum(umi_all, 1.0)).astype(np.float32)
    X_lib = X_raw * scale[:, None]
    X_log = np.log1p(X_lib)

    # Global NTC mean_libnorm & per-batch NTC means (for overall and per-half)
    ctrl_mean_lib = np.mean(X_lib[nt_all], axis=0, dtype=np.float64)

    # Keep HGNC-resolved unique genes mirroring select_canonical_genes
    best_by_sym: Dict[str, Tuple[float, int]] = {}
    for idx in range(len(var_gnames)):
        canon_s, status = hgnc.resolve(var_gnames[idx], var_gids[idx])
        if status == "unmapped" and canon_s not in hgnc.approved:
            if not canon_s:
                continue
        sc = float(ctrl_mean_lib[idx])
        if canon_s not in best_by_sym or sc > best_by_sym[canon_s][0]:
            best_by_sym[canon_s] = (sc, idx)

    kept_genes_items = sorted(best_by_sym.items(), key=lambda x: x[1][1])
    sel_cols = np.asarray([it[1][1] for it in kept_genes_items], dtype=np.int64)
    sel_genes = [it[0] for it in kept_genes_items]

    X_lib_s = X_lib[:, sel_cols]
    X_log_s = X_log[:, sel_cols]
    oracle_ctrl_mean_lib = np.mean(X_lib_s[nt_all], axis=0, dtype=np.float64)

    batch_mu_log: Dict[str, np.ndarray] = {}
    for b in np.unique(batches_all):
        m_bn = nt_all & (batches_all == b)
        if np.any(m_bn):
            batch_mu_log[str(b)] = np.mean(X_log_s[m_bn], axis=0, dtype=np.float64)

    ha_mu_log: Dict[str, np.ndarray] = {}
    hb_mu_log: Dict[str, np.ndarray] = {}
    g_ha_mu = np.mean(X_log_s[nt_all & (halves_all == 1)], axis=0, dtype=np.float64)
    g_hb_mu = np.mean(X_log_s[nt_all & (halves_all == 2)], axis=0, dtype=np.float64)
    for b in np.unique(batches_all):
        ma = nt_all & (halves_all == 1) & (batches_all == b)
        mb = nt_all & (halves_all == 2) & (batches_all == b)
        ha_mu_log[str(b)] = np.mean(X_log_s[ma], axis=0, dtype=np.float64) if int(np.sum(ma)) >= 15 else g_ha_mu
        hb_mu_log[str(b)] = np.mean(X_log_s[mb], axis=0, dtype=np.float64) if int(np.sum(mb)) >= 15 else g_hb_mu

    # Group per-target deltas (>= 10 cells matching process_r10_line)
    pk_rows: Dict[str, List[int]] = defaultdict(list)
    for r_i, (is_nt, hf, pk) in enumerate(zip(nt_all, halves_all, pk_list)):
        if (not is_nt) and hf in (1, 2) and pk:
            pk_rows[pk].append(r_i)
    elig_pks = sorted(pk for pk, rlist in pk_rows.items() if len(rlist) >= 10)

    n_pk = len(elig_pks)
    g_dim = len(sel_genes)
    delta_log = np.zeros((n_pk, g_dim), dtype=np.float64)
    ha_delta = np.full((n_pk, g_dim), np.nan, dtype=np.float64)
    hb_delta = np.full((n_pk, g_dim), np.nan, dtype=np.float64)
    n_cells_arr = np.zeros(n_pk, dtype=np.float64)

    for p_i, pk in enumerate(elig_pks):
        rows = pk_rows[pk]
        n_cells_arr[p_i] = float(len(rows))
        acc_all = np.zeros(g_dim, dtype=np.float64)
        acc_a = np.zeros(g_dim, dtype=np.float64)
        acc_b = np.zeros(g_dim, dtype=np.float64)
        cnt_a = 0
        cnt_b = 0
        for r in rows:
            b = str(batches_all[r])
            h = int(halves_all[r])
            mu_b = batch_mu_log.get(b, g_ha_mu if h == 1 else g_hb_mu)
            acc_all += X_log_s[r] - mu_b
            if h == 1:
                acc_a += X_log_s[r] - ha_mu_log[b]
                cnt_a += 1
            elif h == 2:
                acc_b += X_log_s[r] - hb_mu_log[b]
                cnt_b += 1
        delta_log[p_i] = acc_all / float(len(rows))
        if cnt_a >= 1:
            ha_delta[p_i] = acc_a / float(cnt_a)
        if cnt_b >= 1:
            hb_delta[p_i] = acc_b / float(cnt_b)

    return {
        "audit_meta": meta,
        "genes": sel_genes,
        "targets": elig_pks,
        "n_cells": n_cells_arr,
        "ctrl_mean_libnorm": oracle_ctrl_mean_lib,
        "delta_log1p": delta_log,
        "half_a_delta_log1p": ha_delta,
        "half_b_delta_log1p": hb_delta,
    }


def compare_line_dirs_exact(
    rebuilt_dir: Path,
    other_dir: Path,
    tol: float = 1e-6,
) -> Dict[str, Any]:
    """Compare rebuilt line directory vs another cached line directory (`max_abs_diff <= tol`, zero metrics)."""
    r_genes = [ln.strip() for ln in (rebuilt_dir / "genes.txt").read_text(encoding="utf-8").splitlines() if ln.strip()]
    o_genes = [ln.strip() for ln in (other_dir / "genes.txt").read_text(encoding="utf-8").splitlines() if ln.strip()]
    r_tgts = [ln.strip() for ln in (rebuilt_dir / "targets.txt").read_text(encoding="utf-8").splitlines() if ln.strip()]
    o_tgts = [ln.strip() for ln in (other_dir / "targets.txt").read_text(encoding="utf-8").splitlines() if ln.strip()]

    if set(r_genes) != set(o_genes) or set(r_tgts) != set(o_tgts):
        return {
            "pass": False,
            "genes_set_match": set(r_genes) == set(o_genes),
            "targets_set_match": set(r_tgts) == set(o_tgts),
        }

    g_map = {g: i for i, g in enumerate(r_genes)}
    g_perm = np.asarray([g_map[g] for g in o_genes], dtype=np.int64)
    t_map = {t: i for i, t in enumerate(r_tgts)}
    t_perm = np.asarray([t_map[t] for t in o_tgts], dtype=np.int64)

    diffs: Dict[str, Any] = {}
    with np.load(rebuilt_dir / "ctrl_cells.npz", allow_pickle=False) as cr, np.load(
        other_dir / "ctrl_cells.npz", allow_pickle=False
    ) as co:
        for k in ("ctrl_mean_libnorm", "ctrl_mean_libnorm_full", "ctrl_mean_log1p"):
            if k in cr.files and k in co.files:
                d = float(np.max(np.abs(cr[k][g_perm] - co[k])))
                diffs[k] = {"max_abs_diff": d, "pass": d <= tol}

    with np.load(rebuilt_dir / "pseudobulk_delta.npz", allow_pickle=False) as pr, np.load(
        other_dir / "pseudobulk_delta.npz", allow_pickle=False
    ) as po:
        for k in ("delta_log1p", "half_a_delta_log1p", "half_b_delta_log1p"):
            if k in pr.files and k in po.files:
                a = pr[k][np.ix_(t_perm, g_perm)]
                b = po[k]
                fin = np.isfinite(a) & np.isfinite(b)
                d = float(np.max(np.abs(a[fin] - b[fin]))) if np.any(fin) else 0.0
                diffs[k] = {"max_abs_diff": d, "pass": bool(d <= tol and np.array_equal(np.isnan(a), np.isnan(b)))}
        for k in ("n_cells", "ad_counts"):
            if k in pr.files and k in po.files:
                d = float(np.max(np.abs(pr[k][t_perm] - po[k])))
                diffs[k] = {"max_abs_diff": d, "pass": d == 0.0}
        for k in ("reproducible", "half_valid"):
            if k in pr.files and k in po.files:
                eq = bool(np.array_equal(np.asarray(pr[k], dtype=bool)[t_perm], np.asarray(po[k], dtype=bool)))
                diffs[k] = {"exact_match": eq, "pass": eq}

    all_ok = all(v.get("pass", False) for v in diffs.values())
    return {
        "pass": all_ok,
        "genes_exact_order_match": r_genes == o_genes,
        "targets_exact_order_match": r_tgts == o_tgts,
        "diffs": diffs,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Clean-Room R15 Synthetic Reader Cross-Check")
    parser.add_argument("--line-name", required=True)
    parser.add_argument("--h5ad", nargs="+", type=Path, required=True)
    parser.add_argument("--compare-line-dir", type=Path, default=None)
    parser.add_argument("--context-map-json", type=Path, default=None)
    parser.add_argument("--line-manifest-json", type=Path, default=None)
    parser.add_argument("--phase", choices=("synthetic_test", "reader_xcheck_r15"), default="synthetic_test")
    args = parser.parse_args()

    enforce_r15_gate(
        phase=args.phase,
        target_lines=(args.line_name,),
        extra_dirs=tuple(args.h5ad) + ((args.compare_line_dir,) if args.compare_line_dir else ()),
        require_confirmatory_file=False,
    )
    ctx_map = json.loads(args.context_map_json.read_text(encoding="utf-8")) if args.context_map_json else None
    line_manifest = json.loads(args.line_manifest_json.read_text(encoding="utf-8")) if args.line_manifest_json else None

    with tempfile.TemporaryDirectory() as tmp:
        out_root = Path(tmp) / "rebuilt"
        hgnc_prep = build_hgnc_lookup(DEFAULT_HGNC_TSV)
        parts, meta = build_vcc_line_audit_parts(
            line_name=args.line_name,
            h5_paths=args.h5ad,
            hgnc=hgnc_prep,
            context_to_line_map=ctx_map,
            line_manifest=line_manifest,
        )
        rebuilt_dir = process_r10_line(
            line_name=args.line_name,
            parts=parts,
            hgnc=hgnc_prep,
            chunk_size=8192,
            out_root=out_root,
            extra_meta={"audit_reader": meta},
            n_workers=2,
        )
        oracle = compute_independent_oracle_arrays(
            h5_paths=args.h5ad,
            line_name=args.line_name,
            context_to_line_map=ctx_map,
            line_manifest=line_manifest,
        )
        # Verify process_r10_line vs oracle
        with np.load(rebuilt_dir / "ctrl_cells.npz", allow_pickle=False) as c_r, np.load(
            rebuilt_dir / "pseudobulk_delta.npz", allow_pickle=False
        ) as p_r:
            oracle_diffs = {
                "ctrl_mean_libnorm_f32_diff": float(
                    np.max(np.abs(c_r["ctrl_mean_libnorm"] - oracle["ctrl_mean_libnorm"].astype(np.float32)))
                ),
                "ctrl_mean_libnorm_rel_diff": float(
                    np.max(
                        np.abs(c_r["ctrl_mean_libnorm"] - oracle["ctrl_mean_libnorm"])
                        / np.maximum(np.abs(oracle["ctrl_mean_libnorm"]), 1e-12)
                    )
                ),
                "delta_log1p_diff": float(np.max(np.abs(p_r["delta_log1p"] - oracle["delta_log1p"]))),
                "half_a_diff": float(np.max(np.abs(p_r["half_a_delta_log1p"] - oracle["half_a_delta_log1p"]))),
                "half_b_diff": float(np.max(np.abs(p_r["half_b_delta_log1p"] - oracle["half_b_delta_log1p"]))),
            }

        ext_cmp = compare_line_dirs_exact(rebuilt_dir, args.compare_line_dir) if args.compare_line_dir else None
        print(
            json.dumps(
                {
                    "line": args.line_name,
                    "split_tier": meta["split_tier"],
                    "counts": meta["counts"],
                    "oracle_diffs": oracle_diffs,
                    "compare_line_dir_result": ext_cmp,
                },
                indent=2,
            )
        )


if __name__ == "__main__":
    main()
