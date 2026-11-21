#!/usr/bin/env python3
"""
Unified data loader for `p2/shared/lines/<line>/` with strict pre-registration gate.

Enforces Rule 3 of P2_RULES.md:
  Before `p2/status/PREREG_COMMIT.txt` exists and is non-empty, evaluation is strictly
  restricted to the two Replogle development cell lines (K562, RPE1).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

P2_ROOT = Path(__file__).resolve().parents[2]
SHARED_DIR = P2_ROOT / "shared"
LINES_DIR = SHARED_DIR / "lines"
STATUS_DIR = P2_ROOT / "status"
PREREG_COMMIT_FILE = STATUS_DIR / "PREREG_COMMIT.txt"

REPLOGLE_DEV_PREFIXES = ("replogle", "k562", "rpe1")


@dataclass
class LineDataset:
    """Loaded single-line dataset aligned onto common genes."""

    name: str
    genes: List[str]
    gene_to_idx: Dict[str, int]
    ctrl_cells: np.ndarray                 # (N_ctrl, G) float64
    deltas: np.ndarray                     # (N_pert, G) float64
    target_names: Tuple[str, ...]          # length N_pert
    target_indices: np.ndarray             # (N_pert,) int64 in 0..G-1 or -1
    repr_mask: Optional[np.ndarray]        # (N_pert,) bool if available
    representation: str
    half_a_deltas: Optional[np.ndarray] = None  # (N_pert, G) float64
    half_b_deltas: Optional[np.ndarray] = None  # (N_pert, G) float64
    half_valid: Optional[np.ndarray] = None     # (N_pert,) bool
    ctrl_mean_libnorm: Optional[np.ndarray] = None  # (G,) float64 unperturbed libnorm mean (for FCT)


def is_replogle_dev_line(line_name: str) -> bool:
    low = line_name.lower()
    return any(tok in low for tok in REPLOGLE_DEV_PREFIXES)


def check_prereg_gate(line_names: Sequence[str]) -> None:
    """
    Raise RuntimeError if any non-Replogle line is requested before PREREG_COMMIT.txt exists.
    """
    non_dev = [ln for ln in line_names if not is_replogle_dev_line(ln)]
    if not non_dev:
        return
    if not PREREG_COMMIT_FILE.exists():
        raise RuntimeError(
            f"PREREG GATE BLOCKED: {PREREG_COMMIT_FILE} does not exist yet. "
            f"Cannot load/evaluate non-Replogle lines: {non_dev}"
        )
    content = PREREG_COMMIT_FILE.read_text(encoding="utf-8").strip()
    if not content:
        raise RuntimeError(
            f"PREREG GATE BLOCKED: {PREREG_COMMIT_FILE} is empty. "
            f"Cannot load/evaluate non-Replogle lines: {non_dev}"
        )


def list_available_lines() -> List[str]:
    if not LINES_DIR.exists():
        return []
    return sorted(p.name for p in LINES_DIR.iterdir() if p.is_dir())


def _read_text_lines(path: Path) -> List[str]:
    return [ln.strip() for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]


def _pick_npz_array(npz: np.lib.npyio.NpzFile, candidates: Sequence[str]) -> np.ndarray:
    for key in candidates:
        if key in npz.files:
            return np.asarray(npz[key])
    raise KeyError(f"None of keys {candidates} found in npz (available: {npz.files})")


def load_line_dataset(
    line_name: str,
    representation: str = "log1p",
    common_genes: Optional[Sequence[str]] = None,
) -> LineDataset:
    """
    Load one line from `p2/shared/lines/<line_name>/`, optionally subset/reorder onto `common_genes`.
    Checks PREREG_COMMIT gate first so non-Replogle files are never even read prior to commit.
    """
    check_prereg_gate([line_name])
    line_dir = LINES_DIR / line_name
    if not line_dir.exists():
        raise FileNotFoundError(f"Line directory not found: {line_dir}")

    genes = _read_text_lines(line_dir / "genes.txt")
    targets = _read_text_lines(line_dir / "targets.txt")

    rep_candidates_ctrl = {
        "log1p": ["X_log1p", "log1p", "ctrl_log1p", "X", "ctrl_cells", "data"],
        "zscore": ["X_zscore", "zscore", "ctrl_zscore", "Z", "X", "data"],
        "libnorm": ["X_libnorm", "libnorm", "ctrl_libnorm", "X_norm", "X"],
    }.get(representation, [representation, "X", "data"])

    rep_candidates_delta = {
        "log1p": ["delta_log1p", "log1p", "deltas_log1p", "delta", "deltas", "X"],
        "zscore": ["delta_zscore", "zscore", "deltas_zscore", "delta", "deltas", "X"],
        "libnorm": ["delta_libnorm", "libnorm", "deltas_libnorm", "delta", "deltas", "X"],
    }.get(representation, [representation, "delta", "deltas", "X"])

    with np.load(line_dir / "ctrl_cells.npz", allow_pickle=True) as z_ctrl:
        ctrl_raw = _pick_npz_array(z_ctrl, rep_candidates_ctrl).astype(np.float32, copy=False)
        ctrl_mean_libnorm = (
            np.asarray(z_ctrl["ctrl_mean_libnorm"], dtype=np.float64)
            if "ctrl_mean_libnorm" in z_ctrl.files
            else None
        )

    with np.load(line_dir / "pseudobulk_delta.npz", allow_pickle=True) as z_del:
        deltas_raw = _pick_npz_array(z_del, rep_candidates_delta).astype(np.float64)
        if "targets" in z_del.files:
            del_targets = [str(x) for x in z_del["targets"].tolist()]
            if len(del_targets) == deltas_raw.shape[0]:
                targets = del_targets
        repr_mask = None
        for rk in ("reproducible", "repr_mask", "is_reproducible"):
            if rk in z_del.files:
                repr_mask = np.asarray(z_del[rk]).astype(bool)
                break
        ha_key = f"half_a_delta_{representation}"
        hb_key = f"half_b_delta_{representation}"
        half_a_raw = (
            np.asarray(z_del[ha_key], dtype=np.float64) if ha_key in z_del.files else None
        )
        half_b_raw = (
            np.asarray(z_del[hb_key], dtype=np.float64) if hb_key in z_del.files else None
        )
        half_valid = (
            np.asarray(z_del["half_valid"]).astype(bool) if "half_valid" in z_del.files else None
        )

    if ctrl_raw.shape[1] != len(genes) or deltas_raw.shape[1] != len(genes):
        raise ValueError(
            f"{line_name}: gene dimension mismatch: genes={len(genes)}, "
            f"ctrl={ctrl_raw.shape}, deltas={deltas_raw.shape}"
        )
    if deltas_raw.shape[0] != len(targets):
        raise ValueError(
            f"{line_name}: target count mismatch: targets={len(targets)}, deltas={deltas_raw.shape[0]}"
        )

    if common_genes is not None:
        g_pos = {g: i for i, g in enumerate(genes)}
        missing = [g for g in common_genes if g not in g_pos]
        if missing:
            raise ValueError(f"{line_name}: missing {len(missing)} common genes (first 5: {missing[:5]})")
        col_idx = np.array([g_pos[g] for g in common_genes], dtype=np.int64)
        ctrl_raw = ctrl_raw[:, col_idx]
        deltas_raw = deltas_raw[:, col_idx]
        if half_a_raw is not None:
            half_a_raw = half_a_raw[:, col_idx]
        if half_b_raw is not None:
            half_b_raw = half_b_raw[:, col_idx]
        if ctrl_mean_libnorm is not None:
            ctrl_mean_libnorm = ctrl_mean_libnorm[col_idx]
        genes = list(common_genes)

    g2i = {g: i for i, g in enumerate(genes)}
    t_idx = np.array([g2i.get(t, -1) for t in targets], dtype=np.int64)

    return LineDataset(
        name=line_name,
        genes=genes,
        gene_to_idx=g2i,
        ctrl_cells=ctrl_raw,
        deltas=deltas_raw,
        target_names=tuple(targets),
        target_indices=t_idx,
        repr_mask=repr_mask,
        representation=representation,
        half_a_deltas=half_a_raw,
        half_b_deltas=half_b_raw,
        half_valid=half_valid,
        ctrl_mean_libnorm=ctrl_mean_libnorm,
    )


def load_common_genes() -> Optional[List[str]]:
    cg_file = SHARED_DIR / "COMMON_GENES.txt"
    if not cg_file.exists():
        return None
    return _read_text_lines(cg_file)
