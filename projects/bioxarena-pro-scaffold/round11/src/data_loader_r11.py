#!/usr/bin/env python3
"""
Round 11 Data Loader & Blinding Gate (`PREREG_D`).

Imports frozen Round 10 loader and canonical HGNC bridge (`round10/src/data_loader_r10.py`)
without modifying any frozen Round 10 file.

Supports:
  1. Hard runtime blinding gate (`check_prereg_d_gate`): blocks `--phase holdout_r11`
     unless BOTH `PREREG_D_COMMIT.txt` and `CODE_FREEZE_R11.txt` exist and are non-empty.
  2. Development testing strictly on already-unblinded Round 9 (`dev_r9`) and Round 10
     (`dev_r10`) lines; zero access to any Round 11 held-out line prior to freeze.
  3. Primary 6-line Round 9 source pool (`r9_6`) and sensitivity 8-line Round 9 +
     Round 10 confirmatory (`dld1`, `cd4t_rest`) source pool (`r9_r10_8`).
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

ROUND11_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = ROUND11_ROOT.parent
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
if str(R10_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(R10_SRC_DIR))

from data_loader_r10 import (  # noqa: E402
    DEFAULT_HGNC_FILE,
    HGNCBridge,
    LineDatasetR10,
    R9_COMMON_GENES_FILE,
    R9_LINES_DIR,
    R9_SOURCE_LINE_ORDER,
    R10_LINES_DIR,
    list_standardized_lines,
    load_r9_common_genes,
    load_standardized_line_r10,
    read_gate_hash,
)

R11_SHARED_DIR = ROUND11_ROOT / "shared"
R11_LINES_DIR = R11_SHARED_DIR / "lines"
R11_STATUS_DIR = ROUND11_ROOT / "status"

PREREG_D_COMMIT_FILE = R11_STATUS_DIR / "PREREG_D_COMMIT.txt"
CODE_FREEZE_R11_FILE = R11_STATUS_DIR / "CODE_FREEZE_R11.txt"
R11_CONFIRMATORY_LINES_FILE = R11_STATUS_DIR / "CONFIRMATORY_LINES_R11.txt"

R10_EXTRA_SOURCE_LINES: Tuple[str, ...] = ("dld1", "cd4t_rest")
R11_ALL_8_SOURCE_LINES: Tuple[str, ...] = R9_SOURCE_LINE_ORDER + R10_EXTRA_SOURCE_LINES


def check_prereg_d_gate(status_dir: Path = R11_STATUS_DIR) -> Tuple[str, str]:
    """
    Hard blinding lock: refuse any Round 11 held-out evaluation unless BOTH
    `PREREG_D_COMMIT.txt` and `CODE_FREEZE_R11.txt` exist and are non-empty.
    """
    prereg_file = Path(status_dir) / "PREREG_D_COMMIT.txt"
    freeze_file = Path(status_dir) / "CODE_FREEZE_R11.txt"
    prereg_hash = read_gate_hash(prereg_file)
    freeze_hash = read_gate_hash(freeze_file)
    if not prereg_hash:
        raise RuntimeError(
            f"PREREG_D hard gate blocked: missing/empty {prereg_file}. "
            "Round 11 held-out lines cannot be evaluated before PREREG_D lock."
        )
    if not freeze_hash:
        raise RuntimeError(
            f"CODE_FREEZE_R11 hard gate blocked: missing/empty {freeze_file}. "
            "Round 11 held-out lines cannot be evaluated before Round 11 code freeze."
        )
    return prereg_hash, freeze_hash


def parse_confirmatory_lines_file(conf_file: Path = R11_CONFIRMATORY_LINES_FILE) -> List[str]:
    target = conf_file
    if not target.exists() and target.name == "CONFIRMATORY_LINES_R11.txt":
        fallback = target.parent / "CONFIRMATORY_LINES.txt"
        if fallback.exists():
            target = fallback
    if not target.exists():
        return []
    txt = target.read_text(encoding="utf-8").strip()
    if not txt:
        return []
    return [x.strip() for x in txt.splitlines() if x.strip() and not x.strip().startswith("#")]


def load_source_cohort_r11(
    source_mode: str,
    representation: str,
    common_genes: Sequence[str],
    hgnc_bridge: Optional[HGNCBridge],
) -> Dict[str, LineDatasetR10]:
    """
    Load fixed source datasets strictly from unblinded Round 9 (`r9_6`) or
    Round 9 + Round 10 confirmatory `dld1`, `cd4t_rest` (`r9_r10_8`).
    """
    if source_mode not in ("r9_6", "r9_r10_8"):
        raise ValueError(f"Unsupported source_mode: {source_mode}")

    avail_r9 = list_standardized_lines(R9_LINES_DIR)
    if not avail_r9:
        raise RuntimeError(f"Missing Round 9 standardized lines under {R9_LINES_DIR}.")

    sources: Dict[str, LineDatasetR10] = {
        ln: load_standardized_line_r10(
            R9_LINES_DIR / ln,
            representation=representation,
            common_genes=common_genes,
            hgnc_bridge=hgnc_bridge,
        )
        for ln in avail_r9
    }

    if source_mode == "r9_r10_8":
        for ln in R10_EXTRA_SOURCE_LINES:
            line_dir = R10_LINES_DIR / ln
            if not line_dir.exists():
                raise RuntimeError(f"Missing Round 10 confirmatory source line directory: {line_dir}")
            sources[ln] = load_standardized_line_r10(
                line_dir,
                representation=representation,
                common_genes=common_genes,
                hgnc_bridge=hgnc_bridge,
            )
    return sources
