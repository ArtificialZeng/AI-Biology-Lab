#!/usr/bin/env python3
"""
Clean-Room Round 13 Reader Cross-Check (`mac_ctrl`, `mac_lps`) — Counts & Arrays Only (A-T3).

Strictly follows `PREREG_F.md`:
  - Requires `PREREG_F_COMMIT.txt` + `R13_LINES_AUTHORIZED.txt` / `R13_CODE_PATCHED`
    via `enforce_r13_gate("reader_xcheck_r13")`.
  - Never reads `round13/src/` and never computes any predictor, correlation, or metric.
  - Builds independent standardized caches under `round13/audit/reader_xcheck/` and compares
    counts, per-channel/half stats, guide thresholds, reproducible/eval target sets, `G_eval`,
    control means, and pseudobulk delta arrays item-by-item against `round13/shared/lines/`.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, Set

import numpy as np

AUDIT_R13_DIR = Path(__file__).resolve().parent
REPO_ROOT = AUDIT_R13_DIR.parents[1]
R10_SCRIPTS_DIR = REPO_ROOT / "round10" / "shared" / "scripts"
R11_AUDIT_DIR = REPO_ROOT / "round11" / "audit"
for _p in (str(AUDIT_R13_DIR), str(R10_SCRIPTS_DIR), str(R11_AUDIT_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from independent_r11_loader import (  # noqa: E402
    DEFAULT_HGNC_TSV,
    R9_LINE_NAMES,
    R9_LINES_DIR,
    CanonicalHGNC,
    load_r9_common_genes,
)
from independent_r13_loader import (  # noqa: E402
    R13_LINES_DIR,
    R13_RAW_DIR,
    R13_STATUS_DIR,
    build_gse327644_audit_parts,
    enforce_r13_gate,
)
from prepare_r10_lines import build_hgnc_lookup, process_r10_line  # noqa: E402


def _compute_u6_rep_hgnc(hgnc_canon: CanonicalHGNC) -> Set[str]:
    u6_rep: Set[str] = set()
    for sn in R9_LINE_NAMES:
        sdir = R9_LINES_DIR / sn
        tgts = [x.strip() for x in (sdir / "targets.txt").read_text(encoding="utf-8").splitlines() if x.strip()]
        with np.load(sdir / "pseudobulk_delta.npz", allow_pickle=False) as npz:
            rep = np.asarray(npz["reproducible"], dtype=bool)
            ens_arr = [str(x) for x in npz["target_ensembl_ids"]] if "target_ensembl_ids" in npz.files else [""] * len(tgts)
        for t_s, e_s, r_b in zip(tgts, ens_arr, rep):
            if r_b:
                u6_rep.add(hgnc_canon.resolve(t_s, e_s))
    return u6_rep


def _eval_sets(line_dir: Path, hgnc_canon: CanonicalHGNC, common_r9: Set[str], u6_rep: Set[str]) -> Dict[str, Any]:
    raw_genes = [x.strip() for x in (line_dir / "genes.txt").read_text(encoding="utf-8").splitlines() if x.strip()]
    canon_genes = [hgnc_canon.resolve(g) for g in raw_genes]
    g_eval_set = set(canon_genes) & common_r9

    raw_targets = [x.strip() for x in (line_dir / "targets.txt").read_text(encoding="utf-8").splitlines() if x.strip()]
    with np.load(line_dir / "pseudobulk_delta.npz", allow_pickle=False) as pb:
        rep = np.asarray(pb["reproducible"], dtype=bool)
        hv = np.asarray(pb["half_valid"], dtype=bool)
        ens_arr = [str(x) for x in pb["target_ensembl_ids"]] if "target_ensembl_ids" in pb.files else [""] * len(raw_targets)

    canon_targets = [hgnc_canon.resolve(t, e) for t, e in zip(raw_targets, ens_arr)]
    rep_set = sorted({ct for ct, r in zip(canon_targets, rep) if r})
    eval_set = sorted({
        ct for ct, r, h in zip(canon_targets, rep, hv)
        if r and h and (ct in u6_rep) and (ct in g_eval_set)
    })
    return {
        "n_genes": len(raw_genes),
        "raw_genes": raw_genes,
        "g_eval_count": len(g_eval_set),
        "raw_targets": raw_targets,
        "canon_targets": canon_targets,
        "rep_set": rep_set,
        "eval_set": eval_set,
    }


def _compare_lines(
    line_name: str,
    xdir: Path,
    sdir: Path,
    xmeta: Dict[str, Any],
    hgnc_canon: CanonicalHGNC,
    common_r9: Set[str],
    u6_rep: Set[str],
) -> Dict[str, Any]:
    xm_file = json.loads((xdir / "line_metadata.json").read_text(encoding="utf-8"))
    sm_file = json.loads((sdir / "line_metadata.json").read_text(encoding="utf-8"))

    xs = _eval_sets(xdir, hgnc_canon, common_r9, u6_rep)
    ss = _eval_sets(sdir, hgnc_canon, common_r9, u6_rep)

    # Compare ctrl_cells.npz arrays
    ctrl_diffs: Dict[str, float] = {}
    with np.load(xdir / "ctrl_cells.npz", allow_pickle=False) as xc, np.load(sdir / "ctrl_cells.npz", allow_pickle=False) as sc:
        n_ctrl_x = int(xc["n_ctrl_total"])
        n_ctrl_s = int(sc["n_ctrl_total"])
        x_batches = list(xc["batch_labels"].astype(str))
        s_batches = list(sc["batch_labels"].astype(str))
        common_g_eq = xs["raw_genes"] == ss["raw_genes"]
        if common_g_eq:
            for k in (
                "ctrl_mean_libnorm",
                "ctrl_mean_libnorm_full",
                "ctrl_mean_log1p",
                "ctrl_mean_log1p_full",
                "ctrl_mean_libnorm_sub20k",
                "ctrl_mean_log1p_sub20k",
                "half_a_ctrl_mean_libnorm",
                "half_b_ctrl_mean_libnorm",
            ):
                if k in xc.files and k in sc.files:
                    ctrl_diffs[k] = float(np.max(np.abs(np.asarray(xc[k], dtype=np.float64) - np.asarray(sc[k], dtype=np.float64))))
        else:
            # Align common genes if gene lists differ
            sg_pos = {g: i for i, g in enumerate(ss["raw_genes"])}
            shared_pairs = [(i, sg_pos[g]) for i, g in enumerate(xs["raw_genes"]) if g in sg_pos]
            xi = np.array([a for a, _ in shared_pairs], dtype=np.int64)
            si = np.array([b for _, b in shared_pairs], dtype=np.int64)
            for k in (
                "ctrl_mean_libnorm",
                "ctrl_mean_libnorm_full",
                "ctrl_mean_log1p",
                "ctrl_mean_log1p_full",
                "half_a_ctrl_mean_libnorm",
                "half_b_ctrl_mean_libnorm",
            ):
                if k in xc.files and k in sc.files:
                    ctrl_diffs[k] = float(np.max(np.abs(np.asarray(xc[k], dtype=np.float64)[xi] - np.asarray(sc[k], dtype=np.float64)[si])))

    # Compare pseudobulk_delta.npz arrays aligned on shared (target, gene)
    pb_diffs: Dict[str, float] = {}
    with np.load(xdir / "pseudobulk_delta.npz", allow_pickle=False) as xp, np.load(sdir / "pseudobulk_delta.npz", allow_pickle=False) as spb:
        st_pos = {t: i for i, t in enumerate(ss["raw_targets"])}
        tpairs = [(i, st_pos[t]) for i, t in enumerate(xs["raw_targets"]) if t in st_pos]
        xt_idx = np.array([a for a, _ in tpairs], dtype=np.int64)
        st_idx = np.array([b for _, b in tpairs], dtype=np.int64)

        sg_pos = {g: i for i, g in enumerate(ss["raw_genes"])}
        gpairs = [(i, sg_pos[g]) for i, g in enumerate(xs["raw_genes"]) if g in sg_pos]
        xg_idx = np.array([a for a, _ in gpairs], dtype=np.int64)
        sg_idx = np.array([b for _, b in gpairs], dtype=np.int64)

        for k in ("n_cells", "half_a_n_cells", "half_b_n_cells", "ad_counts", "energy_test_p_value", "knockdown_depth_log1p"):
            if k in xp.files and k in spb.files:
                a = np.asarray(xp[k], dtype=np.float64)[xt_idx]
                b = np.asarray(spb[k], dtype=np.float64)[st_idx]
                valid = np.isfinite(a) & np.isfinite(b)
                pb_diffs[k] = float(np.max(np.abs(a[valid] - b[valid]))) if np.any(valid) else 0.0

        for k in ("delta_log1p", "half_a_delta_log1p", "half_b_delta_log1p", "delta_libnorm", "delta_zscore"):
            if k in xp.files and k in spb.files:
                a = np.asarray(xp[k], dtype=np.float64)[np.ix_(xt_idx, xg_idx)]
                b = np.asarray(spb[k], dtype=np.float64)[np.ix_(st_idx, sg_idx)]
                pb_diffs[k] = float(np.max(np.abs(a - b)))

    return {
        "line": line_name,
        "xcheck_audit_meta": xmeta,
        "xcheck_line_meta": xm_file,
        "shared_line_meta": sm_file,
        "shared_raw_genes_count": ss["n_genes"],
        "xcheck_raw_genes_count": xs["n_genes"],
        "genes_exact_equal": xs["raw_genes"] == ss["raw_genes"],
        "shared_g_eval": ss["g_eval_count"],
        "xcheck_g_eval": xs["g_eval_count"],
        "shared_targets_count": len(ss["raw_targets"]),
        "xcheck_targets_count": len(xs["raw_targets"]),
        "targets_only_in_xcheck": sorted(set(xs["raw_targets"]) - set(ss["raw_targets"])),
        "targets_only_in_shared": sorted(set(ss["raw_targets"]) - set(xs["raw_targets"])),
        "targets_exact_equal": xs["raw_targets"] == ss["raw_targets"],
        "shared_rep_targets": ss["rep_set"],
        "xcheck_rep_targets": xs["rep_set"],
        "rep_targets_exact_equal": xs["rep_set"] == ss["rep_set"],
        "shared_eval_targets": ss["eval_set"],
        "xcheck_eval_targets": xs["eval_set"],
        "eval_targets_exact_equal": xs["eval_set"] == ss["eval_set"],
        "shared_n_ctrl_total": n_ctrl_s,
        "xcheck_n_ctrl_total": n_ctrl_x,
        "ctrl_batch_labels_equal": x_batches == s_batches,
        "n_shared_targets_compared": len(tpairs),
        "n_shared_genes_compared": len(gpairs),
        "ctrl_max_abs_diff": ctrl_diffs,
        "pb_max_abs_diff": pb_diffs,
    }


def main() -> None:
    gate = enforce_r13_gate("reader_xcheck_r13", target_lines=("mac_ctrl", "mac_lps"), status_dir=R13_STATUS_DIR)
    print("GATE_OK:", gate, flush=True)

    out_root = AUDIT_R13_DIR / "reader_xcheck"
    out_root.mkdir(parents=True, exist_ok=True)
    hgnc_prep = build_hgnc_lookup(DEFAULT_HGNC_TSV)
    hgnc_canon = CanonicalHGNC.from_tsv(DEFAULT_HGNC_TSV)
    common_r9 = set(load_r9_common_genes(hgnc=hgnc_canon))
    u6_rep = _compute_u6_rep_hgnc(hgnc_canon)

    results: Dict[str, Any] = {"gate": gate, "lines": {}}
    t0 = time.time()
    for ln in ("mac_ctrl", "mac_lps"):
        print(f"\n=== Independent clean-room rebuild: {ln} ===", flush=True)
        parts, xmeta = build_gse327644_audit_parts(line_name=ln, raw_dir=R13_RAW_DIR, hgnc=hgnc_prep)
        xdir = process_r10_line(
            line_name=ln,
            parts=parts,
            hgnc=hgnc_prep,
            chunk_size=4096,
            out_root=out_root,
            extra_meta={"audit_reader": xmeta},
            n_workers=4,
        )
        cmp_res = _compare_lines(
            line_name=ln,
            xdir=xdir,
            sdir=R13_LINES_DIR / ln,
            xmeta=xmeta,
            hgnc_canon=hgnc_canon,
            common_r9=common_r9,
            u6_rep=u6_rep,
        )
        results["lines"][ln] = cmp_res

    summary_p = out_root / "reader_xcheck_comparison.json"
    summary_p.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n=== XCHECK FINISHED in {time.time()-t0:.1f}s -> {summary_p} ===", flush=True)
    for ln, r in results["lines"].items():
        print(
            f"[{ln}] kept={r['xcheck_audit_meta']['n_kept_cells']} (shared={r['shared_line_meta'].get('audit_reader', {}).get('n_qc_kept_cells', r['shared_line_meta'].get('n_cells_raw'))}), "
            f"ntc_A/B={r['xcheck_audit_meta']['n_ntc_half_a']}/{r['xcheck_audit_meta']['n_ntc_half_b']}, "
            f"targets_eq={r['targets_exact_equal']} (only_x={r['targets_only_in_xcheck']}), "
            f"rep_eq={r['rep_targets_exact_equal']} ({len(r['xcheck_rep_targets'])} vs {len(r['shared_rep_targets'])}), "
            f"eval_eq={r['eval_targets_exact_equal']} ({len(r['xcheck_eval_targets'])} vs {len(r['shared_eval_targets'])}), "
            f"ctrl_log1p_diff={r['ctrl_max_abs_diff'].get('ctrl_mean_log1p')}, "
            f"delta_log1p_diff={r['pb_max_abs_diff'].get('delta_log1p')}",
            flush=True,
        )


if __name__ == "__main__":
    main()
