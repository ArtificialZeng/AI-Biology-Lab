#!/usr/bin/env python3
"""
Independent Round 11 Data Loader, Canonical HGNC Mapper, & Strict Blinding Gate (A-T3).

Written independently per PREREG_D v2 (§1, §2) and PREREG_D_ADDENDUM_A.md:
  - Strictly NEVER imports or accesses `round11/src/` or `round11/results/`.
  - Supports `--mode dev_r9r10` strictly on unblinded Round 9 (`round9/p2/shared/lines/`)
    and Round 10 (`round10/shared/lines/`) datasets, refusing any access to Round 11
    data paths before freeze.
  - Under `--mode holdout_r11`, requires both non-empty `PREREG_D_COMMIT.txt` AND
    `CODE_FREEZE_R11.txt` before reading any Round 11 line matrix or
    `CONFIRMATORY_LINES_R11.txt`.
  - Reads confirmatory candidate set strictly from {telohaec, hcasmc, cd4t_stim48},
    keeping `feng_ipsc` secondary-only (`H-D5`) with relaxed targeted-panel overlap floor.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd

ROUND11_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = ROUND11_ROOT.parent
R9_LINES_DIR = REPO_ROOT / "round9" / "p2" / "shared" / "lines"
R9_COMMON_GENES_PATH = REPO_ROOT / "round9" / "p2" / "shared" / "COMMON_GENES.txt"
R10_LINES_DIR = REPO_ROOT / "round10" / "shared" / "lines"
R11_LINES_DIR = ROUND11_ROOT / "shared" / "lines"
R11_STATUS_DIR = ROUND11_ROOT / "status"
DEFAULT_HGNC_TSV = Path("/Users/mac/.cache/r9_p2_raw/hgnc_complete_set.txt")

R9_LINE_NAMES: Tuple[str, ...] = (
    "nadig_hepg2",
    "nadig_jurkat",
    "nourreddine_kolf21j",
    "replogle_k562_essential",
    "replogle_rpe1",
    "vcc2025_h1",
)
R9_LINE_SET: Set[str] = set(R9_LINE_NAMES)

R10_SOURCE_EXT_NAMES: Tuple[str, ...] = ("dld1", "cd4t_rest")
R10_UNBLINDED_LINES: Tuple[str, ...] = R9_LINE_NAMES + ("dld1", "cd4t_rest", "wtc11_i3neuron", "calu3")
R10_UNBLINDED_SET: Set[str] = set(R10_UNBLINDED_LINES)

R11_CANDIDATE_CONFIRMATORY: Tuple[str, ...] = ("telohaec", "hcasmc", "cd4t_stim48")
R11_SECONDARY_ONLY: Tuple[str, ...] = ("feng_ipsc",)


@dataclass
class AuditLineData:
    """Standardized cell-line container aligned for Round 11 FCT / E1 / H-D5 evaluation."""

    name: str
    genes: List[str]
    target_names: List[str]
    target_indices: np.ndarray
    deltas: np.ndarray
    ctrl_cells: np.ndarray
    ctrl_mean_libnorm: np.ndarray
    ctrl_mean_libnorm_full: Optional[np.ndarray]
    repr_mask: Optional[np.ndarray]
    half_a_deltas: Optional[np.ndarray]
    half_b_deltas: Optional[np.ndarray]
    half_valid: Optional[np.ndarray]
    n_cells: Optional[np.ndarray]
    kd_depth_log1p: np.ndarray
    metadata: Dict[str, Any]


def _first_token(path: Path) -> str:
    if not path.exists():
        return ""
    text = path.read_text(encoding="utf-8").strip()
    if not text:
        return ""
    return text.splitlines()[0].strip().split()[0]


def enforce_r11_gate(
    target_lines: Sequence[str],
    dev_only: bool = False,
    status_dir: Path = R11_STATUS_DIR,
    extra_lines_dirs: Optional[Sequence[Path]] = None,
) -> Tuple[str, str]:
    """
    Strictly enforce PREREG_D §1 Step 5 blinding lock:
      - If `dev_only=True`: permit ONLY unblinded Round 9/10 lines and directories
        (`round9/p2/shared/lines` or `round10/shared/lines`), refusing any Round 11 target/dir.
      - If `dev_only=False` (`holdout_r11`): require BOTH `PREREG_D_COMMIT.txt` AND
        `CODE_FREEZE_R11.txt` to exist and be non-empty.
    """
    if dev_only:
        non_dev = [t for t in target_lines if t not in R10_UNBLINDED_SET]
        if non_dev:
            raise RuntimeError(
                f"PREREG_D dev blinding lock: refused Round 11 / unknown target lines {non_dev}."
            )
        if extra_lines_dirs is not None:
            r11_shared = (ROUND11_ROOT / "shared").resolve()
            for d in extra_lines_dirs:
                rp = Path(d).resolve()
                if rp == r11_shared or r11_shared in rp.parents:
                    raise RuntimeError(
                        f"PREREG_D dev blinding lock: refused accessing {d} inside round11/shared/."
                    )
        return _first_token(status_dir / "PREREG_D_COMMIT.txt"), ""

    prereg_hash = _first_token(status_dir / "PREREG_D_COMMIT.txt")
    freeze_hash = _first_token(status_dir / "CODE_FREEZE_R11.txt")
    if not prereg_hash:
        raise RuntimeError(
            f"PREREG_D hard lock: missing or empty {status_dir / 'PREREG_D_COMMIT.txt'}."
        )
    if not freeze_hash:
        raise RuntimeError(
            f"CODE_FREEZE_R11 hard lock: missing or empty {status_dir / 'CODE_FREEZE_R11.txt'}. "
            "Round 11 expression data cannot be accessed before code freeze."
        )
    return prereg_hash, freeze_hash


class CanonicalHGNC:
    """Read-only HGNC symbol canonicalizer using hgnc_complete_set.txt."""

    def __init__(
        self,
        approved: Set[str],
        ens_to_sym: Dict[str, str],
        prev_to_sym: Dict[str, str],
        alias_to_sym: Dict[str, str],
    ) -> None:
        self.approved = approved
        self.ens_to_sym = ens_to_sym
        self.prev_to_sym = prev_to_sym
        self.alias_to_sym = alias_to_sym

    @classmethod
    def from_tsv(cls, path: Path = DEFAULT_HGNC_TSV) -> Optional["CanonicalHGNC"]:
        if not Path(path).exists():
            return None
        df = pd.read_csv(
            path,
            sep="\t",
            usecols=["symbol", "status", "alias_symbol", "prev_symbol", "ensembl_gene_id"],
            dtype=str,
            keep_default_na=False,
        )
        approved: Set[str] = set()
        ens_map: Dict[str, Set[str]] = {}
        prev_map: Dict[str, Set[str]] = {}
        alias_map: Dict[str, Set[str]] = {}
        for row in df.itertuples(index=False):
            if row.status != "Approved":
                continue
            sym = row.symbol.strip()
            if not sym:
                continue
            approved.add(sym)
            ens = row.ensembl_gene_id.strip().split(".")[0]
            if ens:
                ens_map.setdefault(ens, set()).add(sym)
            if row.prev_symbol:
                for p in row.prev_symbol.split("|"):
                    pt = p.strip().strip('"')
                    if pt:
                        prev_map.setdefault(pt, set()).add(sym)
            if row.alias_symbol:
                for a in row.alias_symbol.split("|"):
                    at = a.strip().strip('"')
                    if at:
                        alias_map.setdefault(at, set()).add(sym)
        return cls(
            approved=approved,
            ens_to_sym={k: next(iter(v)) for k, v in ens_map.items() if len(v) == 1},
            prev_to_sym={k: next(iter(v)) for k, v in prev_map.items() if len(v) == 1},
            alias_to_sym={k: next(iter(v)) for k, v in alias_map.items() if len(v) == 1},
        )

    def resolve(self, raw_symbol: Optional[str], ensembl_id: Optional[str] = None) -> Optional[str]:
        sym = (raw_symbol or "").strip()
        ens = (ensembl_id or "").strip().split(".")[0]
        if sym and sym in self.approved:
            return sym
        if ens and ens in self.ens_to_sym:
            return self.ens_to_sym[ens]
        if sym and sym in self.prev_to_sym:
            return self.prev_to_sym[sym]
        if sym and sym in self.alias_to_sym:
            return self.alias_to_sym[sym]
        return sym if sym else None

    def dedup_indices(
        self,
        raw_symbols: Sequence[str],
        ensembl_ids: Optional[Sequence[str]] = None,
        priority: Optional[Sequence[float]] = None,
    ) -> Tuple[List[str], np.ndarray]:
        best: Dict[str, Tuple[float, int]] = {}
        for idx, rs in enumerate(raw_symbols):
            ens = ensembl_ids[idx] if ensembl_ids is not None else None
            canon = self.resolve(rs, ensembl_id=ens)
            if canon is None:
                continue
            sc = float(priority[idx]) if priority is not None else 0.0
            if not np.isfinite(sc):
                sc = -1e18
            prev = best.get(canon)
            if prev is None or sc > prev[0] or (sc == prev[0] and idx < prev[1]):
                best[canon] = (sc, idx)
        ordered = sorted(((pos, sym) for sym, (_, pos) in best.items()), key=lambda x: x[0])
        return [s for _, s in ordered], np.asarray([p for p, _ in ordered], dtype=np.int64)


def load_r9_common_genes(
    path: Path = R9_COMMON_GENES_PATH,
    hgnc: Optional[CanonicalHGNC] = None,
) -> List[str]:
    raw = [x.strip() for x in path.read_text(encoding="utf-8").splitlines() if x.strip()]
    if hgnc is None:
        return raw
    seen: Set[str] = set()
    out: List[str] = []
    for g in raw:
        cg = hgnc.resolve(g)
        if cg and cg not in seen:
            seen.add(cg)
            out.append(cg)
    return out


def _extract_array(npz_obj: Any, keys: Sequence[str]) -> np.ndarray:
    avail = set(npz_obj.files)
    for k in keys:
        if k in avail:
            return np.asarray(npz_obj[k], dtype=np.float64)
    raise KeyError(f"None of keys {keys} found; available={sorted(avail)}")


def load_line_dataset(
    line_dir: Path,
    representation: str = "log1p",
    aligned_genes: Optional[Sequence[str]] = None,
    hgnc: Optional[CanonicalHGNC] = None,
    min_gene_overlap: int = 3000,
) -> AuditLineData:
    """Load one standardized cell line (`round9/p2/shared/lines`, `round10/shared/lines`, or `round11/shared/lines`)."""
    if representation not in ("log1p", "zscore", "libnorm"):
        raise ValueError(f"Unsupported representation: {representation}")

    raw_genes = [x.strip() for x in (line_dir / "genes.txt").read_text(encoding="utf-8").splitlines() if x.strip()]
    raw_targets = [x.strip() for x in (line_dir / "targets.txt").read_text(encoding="utf-8").splitlines() if x.strip()]
    meta_file = line_dir / "line_metadata.json"
    meta = json.loads(meta_file.read_text(encoding="utf-8")) if meta_file.exists() else {}

    with np.load(line_dir / "pseudobulk_delta.npz", allow_pickle=False) as pb_npz:
        pb_targets = (
            [str(x) for x in pb_npz["targets"].tolist()]
            if "targets" in pb_npz.files
            else raw_targets
        )
        target_ens = (
            [str(x) for x in pb_npz["target_ensembl_ids"].tolist()]
            if "target_ensembl_ids" in pb_npz.files
            else None
        )
        deltas_raw = _extract_array(pb_npz, (f"delta_{representation}", representation, "delta"))
        repr_mask = np.asarray(pb_npz["reproducible"], dtype=bool) if "reproducible" in pb_npz.files else None
        ad_counts = np.asarray(pb_npz["ad_counts"], dtype=np.float64) if "ad_counts" in pb_npz.files else None
        n_cells = np.asarray(pb_npz["n_cells"], dtype=np.float64) if "n_cells" in pb_npz.files else None
        ha_raw = (
            np.asarray(pb_npz[f"half_a_delta_{representation}"], dtype=np.float64)
            if f"half_a_delta_{representation}" in pb_npz.files
            else None
        )
        hb_raw = (
            np.asarray(pb_npz[f"half_b_delta_{representation}"], dtype=np.float64)
            if f"half_b_delta_{representation}" in pb_npz.files
            else None
        )
        hv_raw = np.asarray(pb_npz["half_valid"], dtype=bool) if "half_valid" in pb_npz.files else None
        kd_raw = (
            np.asarray(pb_npz["knockdown_depth_log1p"], dtype=np.float64)
            if "knockdown_depth_log1p" in pb_npz.files
            else None
        )

    with np.load(line_dir / "ctrl_cells.npz", allow_pickle=False) as ctrl_npz:
        ctrl_raw = _extract_array(ctrl_npz, (f"X_{representation}", f"ctrl_{representation}", representation, "X"))
        ctrl_mean_lib = (
            np.asarray(ctrl_npz["ctrl_mean_libnorm"], dtype=np.float64)
            if "ctrl_mean_libnorm" in ctrl_npz.files
            else np.asarray(ctrl_npz["X_libnorm"], dtype=np.float64).mean(axis=0)
        )
        ctrl_mean_lib_full = (
            np.asarray(ctrl_npz["ctrl_mean_libnorm_full"], dtype=np.float64)
            if "ctrl_mean_libnorm_full" in ctrl_npz.files
            else None
        )

    if hgnc is not None:
        canon_genes, g_idx = hgnc.dedup_indices(raw_genes, priority=ctrl_mean_lib)
        deltas_raw = deltas_raw[:, g_idx]
        ctrl_raw = ctrl_raw[:, g_idx]
        ctrl_mean_lib = ctrl_mean_lib[g_idx]
        if ctrl_mean_lib_full is not None:
            ctrl_mean_lib_full = ctrl_mean_lib_full[g_idx]
        if ha_raw is not None:
            ha_raw = ha_raw[:, g_idx]
        if hb_raw is not None:
            hb_raw = hb_raw[:, g_idx]

        t_prio = np.zeros(len(pb_targets), dtype=np.float64)
        if repr_mask is not None:
            t_prio += repr_mask.astype(np.float64) * 1e9
        if ad_counts is not None:
            t_prio += np.nan_to_num(ad_counts, nan=0.0) * 1e4
        if n_cells is not None:
            t_prio += np.nan_to_num(n_cells, nan=0.0)
        canon_targets, t_idx = hgnc.dedup_indices(pb_targets, ensembl_ids=target_ens, priority=t_prio)
        deltas_raw = deltas_raw[t_idx]
        if repr_mask is not None:
            repr_mask = repr_mask[t_idx]
        if n_cells is not None:
            n_cells = n_cells[t_idx]
        if ha_raw is not None:
            ha_raw = ha_raw[t_idx]
        if hb_raw is not None:
            hb_raw = hb_raw[t_idx]
        if hv_raw is not None:
            hv_raw = hv_raw[t_idx]
        if kd_raw is not None:
            kd_raw = kd_raw[t_idx]
    else:
        canon_genes = list(raw_genes)
        canon_targets = list(pb_targets)

    if aligned_genes is not None:
        g_map = {g: i for i, g in enumerate(canon_genes)}
        kept_genes = [g for g in aligned_genes if g in g_map]
        if len(kept_genes) < min_gene_overlap:
            raise RuntimeError(
                f"Line {line_dir.name}: overlap {len(kept_genes)} < min_gene_overlap {min_gene_overlap}."
            )
        cols = np.asarray([g_map[g] for g in kept_genes], dtype=np.int64)
        genes_out = kept_genes
        deltas = deltas_raw[:, cols]
        ctrl_cells = ctrl_raw[:, cols]
        ctrl_mean_lib = ctrl_mean_lib[cols]
        if ctrl_mean_lib_full is not None:
            ctrl_mean_lib_full = ctrl_mean_lib_full[cols]
        ha = ha_raw[:, cols] if ha_raw is not None else None
        hb = hb_raw[:, cols] if hb_raw is not None else None
    else:
        genes_out = canon_genes
        deltas = deltas_raw
        ctrl_cells = ctrl_raw
        ha = ha_raw
        hb = hb_raw

    g_to_col = {g: i for i, g in enumerate(genes_out)}
    target_indices = np.asarray([g_to_col.get(t, -1) for t in canon_targets], dtype=np.int64)

    n_t = len(canon_targets)
    if kd_raw is not None and len(kd_raw) == n_t:
        kd_depth = np.where(target_indices >= 0, np.abs(kd_raw), np.nan)
    else:
        kd_depth = np.full(n_t, np.nan, dtype=np.float64)
        for i, gi in enumerate(target_indices):
            if int(gi) >= 0:
                kd_depth[i] = abs(float(deltas[i, int(gi)]))

    return AuditLineData(
        name=line_dir.name,
        genes=genes_out,
        target_names=list(canon_targets),
        target_indices=target_indices,
        deltas=np.asarray(deltas, dtype=np.float64),
        ctrl_cells=np.asarray(ctrl_cells, dtype=np.float64),
        ctrl_mean_libnorm=np.asarray(ctrl_mean_lib, dtype=np.float64),
        ctrl_mean_libnorm_full=np.asarray(ctrl_mean_lib_full, dtype=np.float64) if ctrl_mean_lib_full is not None else None,
        repr_mask=repr_mask,
        half_a_deltas=ha,
        half_b_deltas=hb,
        half_valid=hv_raw,
        n_cells=n_cells,
        kd_depth_log1p=kd_depth,
        metadata=meta,
    )


def parse_confirmatory_lines_r11(path: Path = R11_STATUS_DIR / "CONFIRMATORY_LINES_R11.txt") -> List[str]:
    """
    Deterministically parse `round11/status/CONFIRMATORY_LINES_R11.txt`.
    Only lines in R11_CANDIDATE_CONFIRMATORY ({telohaec, hcasmc, cd4t_stim48}) are accepted;
    `feng_ipsc` is permanently secondary-only (`PREREG_D §2.2.4`).
    """
    if not path.exists():
        return []
    text = path.read_text(encoding="utf-8")
    found: List[str] = []
    valid_set = set(R11_CANDIDATE_CONFIRMATORY)
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line.upper().startswith("CONFIRMATORY_LINES"):
            rhs = line.split("=", 1)[-1] if "=" in line else line.split(":", 1)[-1]
            for tok in rhs.replace(",", " ").split():
                t = tok.strip()
                if t in valid_set and t not in found:
                    found.append(t)
        elif line in valid_set and line not in found:
            found.append(line)
    return found
