#!/usr/bin/env python3
"""
Round 12 Data Loader & Blinding Gate (`PREREG_E`).

Reuses frozen Round 10 and Round 11 loaders (`round10/src/data_loader_r10.py` and
`round11/src/data_loader_r11.py`) without modifying any Round 10 or Round 11 file.

Supports:
  1. Hard runtime blinding gate (`check_prereg_e_gate`): refuses `--phase holdout_r12`
     unless BOTH `PREREG_E_COMMIT.txt` and `CODE_FREEZE_R12.txt` exist and are non-empty.
  2. Development evaluation strictly on already-unblinded Round 9 (`dev_r9`), Round 10
     (`dev_r10`), and Round 11 (`dev_r11`) standardized lines; zero access to any Round 12
     held-out expression data prior to `PREREG_E` lock & code freeze.
  3. Primary 6-line Round 9 source pool (`r9_6`) and 8-line Round 9 + Round 10
     confirmatory (`dld1`, `cd4t_rest`) sensitivity source pool (`r9_r10_8`).
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

ROUND12_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = ROUND12_ROOT.parent
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
R11_SRC_DIR = REPO_ROOT / "round11" / "src"
for p in (str(R10_SRC_DIR), str(R11_SRC_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from data_loader_r10 import (  # noqa: E402
    DEFAULT_HGNC_FILE,
    HGNCBridge,
    LineDatasetR10,
    R9_COMMON_GENES_FILE,
    R9_LINES_DIR,
    R10_LINES_DIR,
    list_standardized_lines,
    load_r9_common_genes,
    load_standardized_line_r10,
    read_gate_hash,
)
from data_loader_r11 import (  # noqa: E402
    R11_LINES_DIR,
    load_source_cohort_r11,
    parse_confirmatory_lines_file,
)

R12_SHARED_DIR = ROUND12_ROOT / "shared"
R12_LINES_DIR = R12_SHARED_DIR / "lines"
R12_STATUS_DIR = ROUND12_ROOT / "status"

PREREG_E_COMMIT_FILE = R12_STATUS_DIR / "PREREG_E_COMMIT.txt"
CODE_FREEZE_R12_FILE = R12_STATUS_DIR / "CODE_FREEZE_R12.txt"
R12_CONFIRMATORY_LINES_FILE = R12_STATUS_DIR / "CONFIRMATORY_LINES_R12.txt"


def check_prereg_e_gate(status_dir: Path = R12_STATUS_DIR) -> Tuple[str, str]:
    """
    Hard blinding lock: refuse any Round 12 held-out evaluation unless BOTH
    `PREREG_E_COMMIT.txt` and `CODE_FREEZE_R12.txt` exist and are non-empty.
    """
    prereg_file = Path(status_dir) / "PREREG_E_COMMIT.txt"
    freeze_file = Path(status_dir) / "CODE_FREEZE_R12.txt"
    prereg_hash = read_gate_hash(prereg_file)
    freeze_hash = read_gate_hash(freeze_file)
    if not prereg_hash:
        raise RuntimeError(
            f"PREREG_E hard gate blocked: missing/empty {prereg_file}. "
            "Round 12 held-out lines cannot be evaluated before PREREG_E lock."
        )
    if not freeze_hash:
        raise RuntimeError(
            f"CODE_FREEZE_R12 hard gate blocked: missing/empty {freeze_file}. "
            "Round 12 held-out lines cannot be evaluated before Round 12 code freeze."
        )
    return prereg_hash, freeze_hash


def parse_confirmatory_lines_r12(conf_file: Path = R12_CONFIRMATORY_LINES_FILE) -> List[str]:
    """Parse non-comment lines from `CONFIRMATORY_LINES_R12.txt`."""
    return parse_confirmatory_lines_file(conf_file)


def load_source_cohort_r12(
    source_mode: str,
    representation: str,
    common_genes: Sequence[str],
    hgnc_bridge: Optional[HGNCBridge],
) -> Dict[str, LineDatasetR10]:
    """Delegate identically to `load_source_cohort_r11` (`r9_6` or `r9_r10_8`)."""
    return load_source_cohort_r11(
        source_mode=source_mode,
        representation=representation,
        common_genes=common_genes,
        hgnc_bridge=hgnc_bridge,
    )


__all__ = [
    "CODE_FREEZE_R12_FILE",
    "DEFAULT_HGNC_FILE",
    "HGNCBridge",
    "LineDatasetR10",
    "PREREG_E_COMMIT_FILE",
    "R9_COMMON_GENES_FILE",
    "R9_LINES_DIR",
    "R10_LINES_DIR",
    "R11_LINES_DIR",
    "R12_CONFIRMATORY_LINES_FILE",
    "R12_LINES_DIR",
    "R12_SHARED_DIR",
    "R12_STATUS_DIR",
    "check_prereg_e_gate",
    "list_standardized_lines",
    "load_r9_common_genes",
    "load_source_cohort_r12",
    "load_standardized_line_r10",
    "parse_confirmatory_lines_r12",
]
