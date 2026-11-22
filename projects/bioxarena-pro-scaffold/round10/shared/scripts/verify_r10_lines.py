#!/usr/bin/env python3
"""
Counts-Only Verification & Confirmatory Eligibility Auditor for R10 Lines (A-T1).

BLINDING INVARIANT:
  Computes and reports ONLY dataset counts, finite checks, and target/gene set
  overlaps under PREREG_C §2.2 (P1-P5). Never loads or computes any model
  predictions, correlations, r_orth_cal, slopes, or FCT metrics.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from prepare_r10_lines import (
    DEFAULT_HGNC_PATH,
    HGNCLookup,
    R9_COMMON_GENES_PATH,
    R10_LINES_DIR,
    build_hgnc_lookup,
)

PROJECT_ROOT = Path(__file__).resolve().parents[3]
R9_LINES_DIR = PROJECT_ROOT / "round9" / "p2" / "shared" / "lines"
CONFIRMATORY_OUT = PROJECT_ROOT / "round10" / "status" / "CONFIRMATORY_LINES.txt"

CANDIDATE_CONFIRMATORY = ("dld1", "cd4t_rest")
SECONDARY_ONLY_LINES = ("wtc11_i3neuron", "calu3")
ALL_R10_LINES = CANDIDATE_CONFIRMATORY + SECONDARY_ONLY_LINES


def compute_r9_reference_sets(hgnc: HGNCLookup, r9_lines_dir: Optional[Path] = None) -> dict:
    """Collect R9 COMMON_GENES and R9 reproducible targets (both raw & canonical HGNC)."""
    r9_dir = r9_lines_dir or R9_LINES_DIR
    common_raw = [x.strip() for x in R9_COMMON_GENES_PATH.read_text(encoding="utf-8").splitlines() if x.strip()]
    common_canon = {hgnc.resolve(g, "")[0] for g in common_raw}

    rep_raw_union: set = set()
    rep_canon_union: set = set()
    per_r9_canon: Dict[str, set] = {}

    if r9_dir.exists():
        for ld in sorted(r9_dir.iterdir()):
            if not ld.is_dir():
                continue
            pb_p = ld / "pseudobulk_delta.npz"
            if not pb_p.exists():
                continue
            with np.load(pb_p, allow_pickle=False) as pb:
                tgs = [str(x).strip() for x in pb["targets"]]
                rep = np.asarray(pb["reproducible"], dtype=bool)
            line_canon = set()
            for tg, is_rep in zip(tgs, rep):
                if is_rep and tg:
                    rep_raw_union.add(tg)
                    ct = hgnc.resolve(tg, "")[0]
                    rep_canon_union.add(ct)
                    line_canon.add(ct)
            per_r9_canon[ld.name] = line_canon

    return {
        "r9_common_raw_set": set(common_raw),
        "r9_common_canon_set": common_canon,
        "r9_rep_raw_union": rep_raw_union,
        "r9_rep_canon_union": rep_canon_union,
        "per_r9_canon": per_r9_canon,
    }


def verify_and_count_line(line_dir: Path, hgnc: HGNCLookup, r9_ref: dict) -> dict:
    """Validate HDF5/NPZ contracts and return counts-only summary for one prepared R10 line."""
    req_files = ["genes.txt", "targets.txt", "ctrl_cells.npz", "pseudobulk_delta.npz", "line_metadata.json"]
    for fn in req_files:
        if not (line_dir / fn).exists():
            raise FileNotFoundError(f"Missing {fn} in {line_dir}")

    genes = [x.strip() for x in (line_dir / "genes.txt").read_text(encoding="utf-8").splitlines() if x.strip()]
    targets = [x.strip() for x in (line_dir / "targets.txt").read_text(encoding="utf-8").splitlines() if x.strip()]
    if len(genes) != len(set(genes)):
        raise ValueError(f"{line_dir.name}: duplicate gene symbols in genes.txt")
    if len(targets) != len(set(targets)):
        raise ValueError(f"{line_dir.name}: duplicate target symbols in targets.txt")

    G = len(genes)
    P = len(targets)
    gene_set = set(genes)
    meta = json.loads((line_dir / "line_metadata.json").read_text(encoding="utf-8"))

    with np.load(line_dir / "ctrl_cells.npz", allow_pickle=False) as ctrl:
        for k in (
            "X_log1p", "X_zscore", "X_libnorm",
            "ctrl_mean_log1p_full", "ctrl_mean_libnorm_full",
            "ctrl_mean_log1p_sub20k", "ctrl_mean_libnorm_sub20k",
            "half_a_ctrl_mean_libnorm", "half_b_ctrl_mean_libnorm",
            "n_ctrl_total", "batch_labels",
        ):
            if k not in ctrl:
                raise KeyError(f"{line_dir.name} ctrl_cells.npz missing {k}")
        for k in ("X_log1p", "X_zscore", "X_libnorm"):
            arr = ctrl[k]
            if arr.ndim != 2 or arr.shape[1] != G or not np.all(np.isfinite(arr)):
                raise ValueError(f"{line_dir.name} invalid ctrl_cells[{k}] shape={arr.shape} or non-finite")
        n_ctrl_sampled = int(ctrl["X_log1p"].shape[0])
        n_ctrl_total = int(ctrl["n_ctrl_total"])

    with np.load(line_dir / "pseudobulk_delta.npz", allow_pickle=False) as pb:
        for k in (
            "delta_log1p", "delta_zscore", "delta_libnorm",
            "half_a_delta_log1p", "half_b_delta_log1p",
            "half_a_delta_zscore", "half_b_delta_zscore",
            "half_a_delta_libnorm", "half_b_delta_libnorm",
            "half_a_n_cells", "half_b_n_cells", "half_valid",
            "targets", "target_raw_symbols", "target_ensembl_ids", "target_hgnc_status",
            "reproducible", "n_cells", "knockdown_depth_log1p", "knockdown_fraction",
            "ad_counts", "energy_test_p_value",
        ):
            if k not in pb:
                raise KeyError(f"{line_dir.name} pseudobulk_delta.npz missing {k}")
        if [str(x) for x in pb["targets"]] != targets:
            raise ValueError(f"{line_dir.name}: targets.txt mismatch with pseudobulk_delta.npz")
        for k in ("delta_log1p", "delta_zscore", "delta_libnorm", "half_a_delta_libnorm", "half_b_delta_libnorm"):
            arr = pb[k]
            if arr.shape != (P, G) or not np.all(np.isfinite(arr)):
                raise ValueError(f"{line_dir.name} invalid {k} shape={arr.shape} or non-finite")

        ha_n = pb["half_a_n_cells"]
        hb_n = pb["half_b_n_cells"]
        half_valid = np.asarray(pb["half_valid"], dtype=bool)
        if not np.array_equal(half_valid, (ha_n >= 2) & (hb_n >= 2)):
            raise ValueError(f"{line_dir.name}: half_valid contract violation")

        n_cells = pb["n_cells"]
        ad_c = pb["ad_counts"]
        ep = pb["energy_test_p_value"]
        rep = np.asarray(pb["reproducible"], dtype=bool)
        if not np.array_equal(rep, (n_cells >= 25) & (ad_c >= 50.0) & (ep < 0.001)):
            raise ValueError(f"{line_dir.name}: reproducible contract violation")

        kd_log = pb["knockdown_depth_log1p"]
        for i, t in enumerate(targets):
            if t in gene_set and not np.isfinite(kd_log[i]):
                raise ValueError(f"{line_dir.name}: missing knockdown_depth_log1p for self-gene {t}")

    # Aligned gene sets & PREREG_C §2.2 P5 shared reproducible overlap counts
    g_eval_canon = r9_ref["r9_common_canon_set"] & gene_set
    g_eval_raw = r9_ref["r9_common_raw_set"] & gene_set

    rep_targets = [t for t, r in zip(targets, rep) if r]
    rep_hv_targets = [t for t, r, hv in zip(targets, rep, half_valid) if r and hv]

    shared_rep_r9_all_canon = [t for t in rep_targets if t in r9_ref["r9_rep_canon_union"]]
    shared_rep_r9_eval_canon = [t for t in shared_rep_r9_all_canon if t in g_eval_canon]
    shared_rep_r9_eval_hv_canon = [t for t in rep_hv_targets if t in r9_ref["r9_rep_canon_union"] and t in g_eval_canon]
    shared_rep_r9_eval_raw = [t for t in rep_targets if t in r9_ref["r9_rep_raw_union"] and t in g_eval_raw]

    return {
        "line": line_dir.name,
        "raw_cells_after_qc": int(meta.get("raw_shape", [0, 0])[0]),
        "n_nt_cells": int(n_ctrl_total),
        "n_nt_sampled": int(n_ctrl_sampled),
        "n_genes_finite": int(meta.get("n_genes_finite_all_batches", G)),
        "n_genes_final": int(G),
        "g_eval_canon": int(len(g_eval_canon)),
        "g_eval_raw": int(len(g_eval_raw)),
        "n_targets_dedup": int(P),
        "n_targets_reproducible": int(len(rep_targets)),
        "n_targets_reproducible_half_valid": int(len(rep_hv_targets)),
        "n_shared_reproducible_with_r9_all_canon": int(len(shared_rep_r9_all_canon)),
        "n_shared_reproducible_eval_canon": int(len(shared_rep_r9_eval_canon)),
        "n_shared_reproducible_eval_half_valid_canon": int(len(shared_rep_r9_eval_hv_canon)),
        "n_shared_reproducible_eval_raw": int(len(shared_rep_r9_eval_raw)),
        "n_half_a_batches": int(meta.get("n_half_a_batches", 0)),
        "n_half_b_batches": int(meta.get("n_half_b_batches", 0)),
    }


def evaluate_confirmatory_lines(counts_by_line: Dict[str, dict]) -> Tuple[List[str], str]:
    """
    Mechanically determine confirmatory lines per PREREG_C §2.2:
      Candidate confirmatory lines: `dld1`, `cd4t_rest`.
      Criteria: `g_eval_canon >= 3000` AND `n_shared_reproducible_eval_canon >= 50`.
    `wtc11_i3neuron` and `calu3` remain secondary-only regardless of counts.
    """
    passing: List[str] = []
    lines_txt = ["# Round 10 Mechanical Confirmatory Line Gate (PREREG_C §2.2 / P2-P5)", ""]
    for ln in CANDIDATE_CONFIRMATORY:
        c = counts_by_line.get(ln)
        if c is None:
            lines_txt.append(f"{ln}: NOT_PROCESSED")
            continue
        ok = (c["g_eval_canon"] >= 3000) and (c["n_shared_reproducible_eval_canon"] >= 50)
        if ok:
            passing.append(ln)
        lines_txt.append(
            f"{ln}: {'PASS_CONFIRMATORY' if ok else 'FAIL_CONFIRMATORY'} "
            f"(eval_shared_rep={c['n_shared_reproducible_eval_canon']} >= 50, |G_eval|={c['g_eval_canon']} >= 3000)"
        )
    for ln in SECONDARY_ONLY_LINES:
        c = counts_by_line.get(ln)
        if c is not None:
            lines_txt.append(
                f"{ln}: SECONDARY_ONLY (eval_shared_rep={c['n_shared_reproducible_eval_canon']}, |G_eval|={c['g_eval_canon']})"
            )
    lines_txt.append("h9: EXCLUDED_PRE_DOWNLOAD (P2: same lineage as R9 H1 hESC)")
    lines_txt.append("")
    lines_txt.append(f"CONFIRMATORY_LINES={','.join(passing) if passing else 'NONE'}")
    lines_txt.append(f"K_CONF={len(passing)}")
    lines_txt.append(f"GATE_STATUS={'CONFIRMATORY_OK' if len(passing) >= 2 else 'STOP_BELOW_2_CONFIRMATORY'}")
    return passing, "\n".join(lines_txt) + "\n"


def main() -> None:
    p = argparse.ArgumentParser(description="Verify prepared R10 lines and count confirmatory eligibility (counts only).")
    p.add_argument("--lines-dir", type=Path, default=R10_LINES_DIR)
    p.add_argument("--write-confirmatory", action="store_true")
    p.add_argument("--confirmatory-out", type=Path, default=CONFIRMATORY_OUT)
    p.add_argument("--json-out", type=Path, default=None)
    args = p.parse_args()

    hgnc = build_hgnc_lookup(DEFAULT_HGNC_PATH)
    r9_ref = compute_r9_reference_sets(hgnc)
    results: Dict[str, dict] = {}
    for ln in ALL_R10_LINES:
        ld = args.lines_dir / ln
        if (ld / "line_metadata.json").exists():
            results[ln] = verify_and_count_line(ld, hgnc, r9_ref)
            print(json.dumps(results[ln], indent=2))
        else:
            print(f"[SKIP] {ln}: not yet prepared at {ld}")

    if args.json_out is not None:
        args.json_out.write_text(json.dumps(results, indent=2), encoding="utf-8")

    passing, conf_text = evaluate_confirmatory_lines(results)
    print("\n" + conf_text)
    if args.write_confirmatory:
        args.confirmatory_out.parent.mkdir(parents=True, exist_ok=True)
        args.confirmatory_out.write_text(conf_text, encoding="utf-8")
        print(f"Wrote {args.confirmatory_out} (passing={passing})")


if __name__ == "__main__":
    main()
