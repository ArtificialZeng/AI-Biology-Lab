#!/usr/bin/env python3
"""
Round 13 (`PREREG_F`) Data Loader & Hard Blinding Gate.

Strictly reuses frozen `round10` / `round11` / `round12` primitives without altering any
file under `round10/`, `round11/`, or `round12/`.

  - Enforces hard blinding gate (`check_prereg_f_gate`): refuses `--phase holdout_r13`
    unless BOTH `PREREG_F_COMMIT.txt` and `CODE_FREEZE_R13.txt` exist and are non-empty.
  - Exposes `DEFAULT_R13_SEED = 20261223` and `CONFIRMATORY_LINES_R13.txt` parser.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

ROUND13_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = ROUND13_ROOT.parent
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
    list_standardized_lines,
    load_r9_common_genes,
    load_source_cohort_r12,
    load_standardized_line_r10,
    parse_confirmatory_lines_r12,
    read_gate_hash,
)

DEFAULT_R13_SEED = 20261223
R13_SHARED_DIR = ROUND13_ROOT / "shared"
R13_LINES_DIR = R13_SHARED_DIR / "lines"
R13_STATUS_DIR = ROUND13_ROOT / "status"

PREREG_F_COMMIT_FILE = R13_STATUS_DIR / "PREREG_F_COMMIT.txt"
CODE_FREEZE_R13_FILE = R13_STATUS_DIR / "CODE_FREEZE_R13.txt"
R13_CONFIRMATORY_LINES_FILE = R13_STATUS_DIR / "CONFIRMATORY_LINES_R13.txt"


def check_prereg_f_gate(status_dir: Path = R13_STATUS_DIR) -> Tuple[str, str]:
    """
    Hard blinding lock: refuse Round 13 held-out evaluation unless BOTH
    `PREREG_F_COMMIT.txt` and `CODE_FREEZE_R13.txt` exist and are non-empty.
    """
    prereg_file = Path(status_dir) / "PREREG_F_COMMIT.txt"
    freeze_file = Path(status_dir) / "CODE_FREEZE_R13.txt"
    prereg_hash = read_gate_hash(prereg_file)
    freeze_hash = read_gate_hash(freeze_file)
    if not prereg_hash:
        raise RuntimeError(
            f"PREREG_F hard gate blocked: missing/empty {prereg_file}. "
            "Round 13 held-out lines cannot be evaluated before PREREG_F lock."
        )
    if not freeze_hash:
        raise RuntimeError(
            f"CODE_FREEZE_R13 hard gate blocked: missing/empty {freeze_file}. "
            "Round 13 held-out lines cannot be evaluated before Round 13 code freeze."
        )
    return prereg_hash, freeze_hash


def parse_confirmatory_lines_r13(conf_file: Path = R13_CONFIRMATORY_LINES_FILE) -> List[str]:
    """Parse non-comment confirmatory target lines from `CONFIRMATORY_LINES_R13.txt`."""
    return parse_confirmatory_lines_r12(conf_file)


def load_source_cohort_r13(
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


__all__ = [
    "CODE_FREEZE_R13_FILE",
    "DEFAULT_HGNC_FILE",
    "DEFAULT_R13_SEED",
    "HGNCBridge",
    "LineDatasetR10",
    "PREREG_F_COMMIT_FILE",
    "R9_COMMON_GENES_FILE",
    "R9_LINES_DIR",
    "R10_LINES_DIR",
    "R11_LINES_DIR",
    "R13_CONFIRMATORY_LINES_FILE",
    "R13_LINES_DIR",
    "R13_SHARED_DIR",
    "R13_STATUS_DIR",
    "check_prereg_f_gate",
    "list_standardized_lines",
    "load_r9_common_genes",
    "load_source_cohort_r13",
    "load_standardized_line_r10",
    "parse_confirmatory_lines_r13",
]
