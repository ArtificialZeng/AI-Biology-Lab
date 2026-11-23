#!/usr/bin/env python3
"""
Round 14 (`PREREG_G` / `GSE337804`) Data Loader, Hard Blinding Gate, & Overlap Registry.

Strictly reuses frozen `round10` / `round11` / `round12` primitives without altering any
file under `round10/`, `round11/`, `round12/`, or `round13/`:
  - Enforces hard blinding gate (`check_prereg_g_gate`): refuses `--phase holdout_r14`
    unless BOTH `PREREG_G_COMMIT.txt` and `CODE_FREEZE_R14.txt` exist and are non-empty.
  - Exposes `DEFAULT_R14_SEED = 20261224` and `CONFIRMATORY_LINES_R14.txt` parser.
  - Exposes `load_r12_a549_tested_targets` to partition evaluated targets cleanly into
    Round 12 `a549_ko` overlap vs novel sets for `PREREG_G` Item 7 descriptive reporting.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

ROUND14_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = ROUND14_ROOT.parent
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
R11_SRC_DIR = REPO_ROOT / "round11" / "src"
R12_SRC_DIR = REPO_ROOT / "round12" / "src"
for p in (str(R10_SRC_DIR), str(R11_SRC_DIR), str(R12_SRC_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from data_loader_r12 import (  # noqa: E402
    DEFAULT_HGNC_FILE,
    HGNCBridge,
    LineDatasetR10,
    R9_COMMON_GENES_FILE,
    R9_LINES_DIR,
    R10_LINES_DIR,
    R11_LINES_DIR,
    R12_LINES_DIR,
    list_standardized_lines,
    load_r9_common_genes,
    load_source_cohort_r12,
    load_standardized_line_r10,
    parse_confirmatory_lines_r12,
    read_gate_hash,
)

DEFAULT_R14_SEED = 20261224
R14_SHARED_DIR = ROUND14_ROOT / "shared"
R14_LINES_DIR = R14_SHARED_DIR / "lines"
R14_STATUS_DIR = ROUND14_ROOT / "status"
R12_A549_LINE_DIR = R12_LINES_DIR / "a549_ko"

PREREG_G_COMMIT_FILE = R14_STATUS_DIR / "PREREG_G_COMMIT.txt"
CODE_FREEZE_R14_FILE = R14_STATUS_DIR / "CODE_FREEZE_R14.txt"
R14_CONFIRMATORY_LINES_FILE = R14_STATUS_DIR / "CONFIRMATORY_LINES_R14.txt"


def check_prereg_g_gate(status_dir: Path = R14_STATUS_DIR) -> Tuple[str, str]:
    """
    Hard blinding lock: refuse Round 14 held-out evaluation unless BOTH
    `PREREG_G_COMMIT.txt` and `CODE_FREEZE_R14.txt` exist and are non-empty.
    """
    prereg_file = Path(status_dir) / "PREREG_G_COMMIT.txt"
    freeze_file = Path(status_dir) / "CODE_FREEZE_R14.txt"
    prereg_hash = read_gate_hash(prereg_file)
    freeze_hash = read_gate_hash(freeze_file)
    if not prereg_hash:
        raise RuntimeError(
            f"PREREG_G hard gate blocked: missing/empty {prereg_file}. "
            "Round 14 held-out lines cannot be evaluated before PREREG_G lock."
        )
    if not freeze_hash:
        raise RuntimeError(
            f"CODE_FREEZE_R14 hard gate blocked: missing/empty {freeze_file}. "
            "Round 14 held-out lines cannot be evaluated before Round 14 code freeze."
        )
    return prereg_hash, freeze_hash


def parse_confirmatory_lines_r14(conf_file: Path = R14_CONFIRMATORY_LINES_FILE) -> List[str]:
    """Parse non-comment confirmatory target lines from `CONFIRMATORY_LINES_R14.txt`."""
    return parse_confirmatory_lines_r12(conf_file)


def load_source_cohort_r14(
    source_mode: str,
    representation: str,
    common_genes: Sequence[str],
    hgnc_bridge: Optional[HGNCBridge],
) -> Dict[str, LineDatasetR10]:
    """Delegate identically to frozen `load_source_cohort_r12` (`r9_6` or `r9_r10_8`)."""
    return load_source_cohort_r12(
        source_mode=source_mode,
        representation=representation,
        common_genes=common_genes,
        hgnc_bridge=hgnc_bridge,
    )


def load_r12_a549_tested_targets(
    a549_dir: Path = R12_A549_LINE_DIR,
    hgnc_bridge: Optional[HGNCBridge] = None,
) -> Set[str]:
    """
    Load the canonical set of perturbation targets tested in frozen Round 12 `a549_ko`
    from `<a549_dir>/targets.txt` (canonicalized via `HGNCBridge` if provided).
    """
    tgt_file = Path(a549_dir) / "targets.txt"
    if not tgt_file.exists():
        return set()
    raw_targets = [
        line.strip()
        for line in tgt_file.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if hgnc_bridge is None:
        return set(raw_targets)
    out: Set[str] = set()
    for sym in raw_targets:
        canon = hgnc_bridge.canonicalize(sym)
        if canon:
            out.add(canon)
    return out


__all__ = [
    "CODE_FREEZE_R14_FILE",
    "DEFAULT_HGNC_FILE",
    "DEFAULT_R14_SEED",
    "HGNCBridge",
    "LineDatasetR10",
    "PREREG_G_COMMIT_FILE",
    "R9_COMMON_GENES_FILE",
    "R9_LINES_DIR",
    "R10_LINES_DIR",
    "R11_LINES_DIR",
    "R12_A549_LINE_DIR",
    "R12_LINES_DIR",
    "R14_CONFIRMATORY_LINES_FILE",
    "R14_LINES_DIR",
    "R14_SHARED_DIR",
    "R14_STATUS_DIR",
    "check_prereg_g_gate",
    "list_standardized_lines",
    "load_r9_common_genes",
    "load_r12_a549_tested_targets",
    "load_source_cohort_r14",
    "load_standardized_line_r10",
    "parse_confirmatory_lines_r14",
]
