#!/usr/bin/env python3
"""
Round 10 Data Loader, Canonical HGNC Bridge, HDF5 Streaming Readers, and Blinding Gate.

Supports:
  1. Strict `PREREG_C` execution gates (`PREREG_C_COMMIT.txt` + `CODE_FREEZE_R10.txt`).
  2. Canonical HGNC resolution via frozen `/Users/mac/.cache/r9_p2_raw/hgnc_complete_set.txt`
     bridging Round 9 frozen lines/targets/COMMON_GENES to Round 10 lines without mutating R9.
  3. Standardized `.npz`/`.json` loading (R9 schema + R10 `S7` `n_cells` / knockdown fields).
  4. Pure-Python `h5py` chunked loaders (<= 4 GB RAM footprint) matching `R10_DOWNLOAD_PLAN.md`:
     - DLD-1 (`GSE337988`) 10x `.filtered_feature_bc_matrix.h5` + `.crispr.umi_correct.h5` + `.hashing.umi_correct.h5`
     - CD4+ T (`GSE314342` Rest) 10x Flex `.h5` + S3 `sgrna_library_metadata` CSV
     - H9 Neurons (`IGVFDS9825FOBP`) RNA `.h5` + lightweight guide CSV metadata pre-check
     - Compact `.h5ad` CSR/dense chunked reader (`wtc11_i3neuron`, `calu3`)
  5. Deterministic lexicographical unit split-half assignment (channels / donors / HTO libs / lanes).
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROUND10_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = ROUND10_ROOT.parent
R9_SHARED_DIR = REPO_ROOT / "round9" / "p2" / "shared"
R9_LINES_DIR = R9_SHARED_DIR / "lines"
R9_COMMON_GENES_FILE = R9_SHARED_DIR / "COMMON_GENES.txt"

R10_SHARED_DIR = ROUND10_ROOT / "shared"
R10_LINES_DIR = R10_SHARED_DIR / "lines"
R10_STATUS_DIR = ROUND10_ROOT / "status"

PREREG_C_COMMIT_FILE = R10_STATUS_DIR / "PREREG_C_COMMIT.txt"
CODE_FREEZE_R10_FILE = R10_STATUS_DIR / "CODE_FREEZE_R10.txt"
CONFIRMATORY_LINES_FILE = R10_STATUS_DIR / "CONFIRMATORY_LINES.txt"
DEFAULT_HGNC_FILE = Path("/Users/mac/.cache/r9_p2_raw/hgnc_complete_set.txt")

R9_SOURCE_LINE_ORDER: Tuple[str, ...] = (
    "nadig_hepg2",
    "nadig_jurkat",
    "nourreddine_kolf21j",
    "replogle_k562_essential",
    "replogle_rpe1",
    "vcc2025_h1",
)


def decode_str_array(arr: Iterable[Any]) -> List[str]:
    """Decode bytes/scalar HDF5 arrays cleanly to Python strings."""
    return [
        x.decode("utf-8", errors="replace") if isinstance(x, (bytes, np.bytes_)) else str(x)
        for x in arr
    ]


class HGNCBridge:
    """Deterministic 4-tier canonical HGNC mapper built from `hgnc_complete_set.txt`."""

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
    def from_tsv(cls, path: Path = DEFAULT_HGNC_FILE) -> "HGNCBridge":
        if not Path(path).exists():
            raise FileNotFoundError(f"Missing HGNC table at {path}")
        df = pd.read_csv(
            path,
            sep="\t",
            usecols=["symbol", "status", "alias_symbol", "prev_symbol", "ensembl_gene_id"],
            dtype=str,
            keep_default_na=False,
        )
        approved: Set[str] = set()
        ens_counts: Dict[str, Set[str]] = {}
        prev_counts: Dict[str, Set[str]] = {}
        alias_counts: Dict[str, Set[str]] = {}

        for row in df.itertuples(index=False):
            if row.status != "Approved":
                continue
            sym = row.symbol.strip()
            if not sym:
                continue
            approved.add(sym)
            ens = row.ensembl_gene_id.strip().split(".")[0]
            if ens:
                ens_counts.setdefault(ens, set()).add(sym)
            if row.prev_symbol:
                for p in row.prev_symbol.split("|"):
                    pt = p.strip().strip('"')
                    if pt:
                        prev_counts.setdefault(pt, set()).add(sym)
            if row.alias_symbol:
                for a in row.alias_symbol.split("|"):
                    at = a.strip().strip('"')
                    if at:
                        alias_counts.setdefault(at, set()).add(sym)

        ens_to_sym = {k: next(iter(v)) for k, v in ens_counts.items() if len(v) == 1}
        prev_to_sym = {k: next(iter(v)) for k, v in prev_counts.items() if len(v) == 1}
        alias_to_sym = {k: next(iter(v)) for k, v in alias_counts.items() if len(v) == 1}
        return cls(
            approved=approved,
            ens_to_sym=ens_to_sym,
            prev_to_sym=prev_to_sym,
            alias_to_sym=alias_to_sym,
        )

    def canonicalize_with_source(
        self, raw_symbol: Optional[str], ensembl_id: Optional[str] = None
    ) -> Tuple[Optional[str], str]:
        """
        Resolve raw symbol + optional ENSG ID using locked precedence:
          1) exact Approved symbol
          2) unambiguous Ensembl gene ID
          3) unambiguous prev_symbol
          4) unambiguous alias_symbol
          5) unresolved -> original symbol (status='unmapped') or None if empty
        """
        sym = (raw_symbol or "").strip()
        ens = (ensembl_id or "").strip().split(".")[0]
        if sym and sym in self.approved:
            return sym, "approved"
        if ens and ens in self.ens_to_sym:
            return self.ens_to_sym[ens], "ensg"
        if sym and sym in self.prev_to_sym:
            return self.prev_to_sym[sym], "prev_symbol"
        if sym and sym in self.alias_to_sym:
            return self.alias_to_sym[sym], "alias_symbol"
        if not sym:
            return None, "empty"
        return sym, "unmapped"

    def canonicalize(self, raw_symbol: Optional[str], ensembl_id: Optional[str] = None) -> Optional[str]:
        canon, _ = self.canonicalize_with_source(raw_symbol, ensembl_id=ensembl_id)
        return canon

    def dedup_map_indices(
        self,
        raw_symbols: Sequence[str],
        ensembl_ids: Optional[Sequence[str]] = None,
        priority_scores: Optional[Sequence[float]] = None,
        require_approved: bool = False,
    ) -> Tuple[List[str], np.ndarray, List[str]]:
        """
        Canonicalize `raw_symbols` and deduplicate many-to-one collisions deterministically
        (keeping highest `priority_scores`, tie-broken by earliest index).
        """
        best: Dict[str, Tuple[float, int, str]] = {}
        for idx, rs in enumerate(raw_symbols):
            ens = ensembl_ids[idx] if ensembl_ids is not None else None
            canon, status = self.canonicalize_with_source(rs, ensembl_id=ens)
            if canon is None:
                continue
            if require_approved and canon not in self.approved:
                continue
            score = float(priority_scores[idx]) if priority_scores is not None else 0.0
            if not np.isfinite(score):
                score = -1e18
            prev = best.get(canon)
            if prev is None or (score > prev[0]) or (score == prev[0] and idx < prev[1]):
                best[canon] = (score, idx, status)
        # Preserve deterministic first-seen/kept order sorted by retained original index
        items = sorted(((info[1], sym, info[2]) for sym, info in best.items()), key=lambda x: x[0])
        kept_syms = [x[1] for x in items]
        kept_idx = np.asarray([x[0] for x in items], dtype=np.int64)
        kept_status = [x[2] for x in items]
        return kept_syms, kept_idx, kept_status


@dataclass
class LineDatasetR10:
    """Standardized expression/control container supporting both R9 and R10 lines."""

    name: str
    genes: List[str]
    target_names: List[str]
    target_indices: np.ndarray
    deltas: np.ndarray                  # (N_t, G) float64
    ctrl_cells: np.ndarray              # (N_c, G) float32/float64
    ctrl_mean_log1p: Optional[np.ndarray]
    ctrl_mean_libnorm: Optional[np.ndarray]
    ctrl_mean_libnorm_full: Optional[np.ndarray]
    repr_mask: Optional[np.ndarray]
    half_a_deltas: Optional[np.ndarray]
    half_b_deltas: Optional[np.ndarray]
    half_valid: Optional[np.ndarray]
    n_cells: Optional[np.ndarray]
    knockdown_depth_log1p: Optional[np.ndarray]
    knockdown_fraction: Optional[np.ndarray]
    target_ensembl_ids: Optional[List[str]]
    representation: str
    metadata: Dict[str, Any]


def read_gate_hash(path: Path) -> str:
    if not path.exists():
        return ""
    txt = path.read_text(encoding="utf-8").strip()
    if not txt:
        return ""
    return txt.splitlines()[0].strip().split()[0]


def check_prereg_c_gate(status_dir: Path = R10_STATUS_DIR) -> Tuple[str, str]:
    """
    Hard lock (`PREREG_C §1 Step 5`): refuse any Round 10 held-out evaluation unless BOTH
    `PREREG_C_COMMIT.txt` and `CODE_FREEZE_R10.txt` exist and are non-empty.
    """
    prereg_file = status_dir / "PREREG_C_COMMIT.txt"
    freeze_file = status_dir / "CODE_FREEZE_R10.txt"
    prereg_hash = read_gate_hash(prereg_file)
    freeze_hash = read_gate_hash(freeze_file)
    if not prereg_hash:
        raise RuntimeError(
            f"PREREG_C hard gate blocked: missing/empty {prereg_file}. "
            "Round 10 held-out lines cannot be evaluated before PREREG_C lock."
        )
    if not freeze_hash:
        raise RuntimeError(
            f"CODE_FREEZE_R10 hard gate blocked: missing/empty {freeze_file}. "
            "Round 10 held-out lines cannot be evaluated before Round 10 code freeze."
        )
    return prereg_hash, freeze_hash


def load_r9_common_genes(
    hgnc_bridge: Optional[HGNCBridge] = None,
    common_genes_file: Path = R9_COMMON_GENES_FILE,
) -> List[str]:
    """Load Round 9 fixed 6,148 common genes, optionally projecting through HGNCBridge."""
    raw = [x.strip() for x in common_genes_file.read_text(encoding="utf-8").splitlines() if x.strip()]
    if hgnc_bridge is None:
        return raw
    # Preserve unique canonicalized symbols deterministically
    seen: Set[str] = set()
    out: List[str] = []
    for s in raw:
        c = hgnc_bridge.canonicalize(s)
        if c and c not in seen:
            seen.add(c)
            out.append(c)
    return out


def list_standardized_lines(lines_dir: Path) -> List[str]:
    if not lines_dir.exists():
        return []
    out: List[str] = []
    for d in sorted(lines_dir.iterdir()):
        if d.is_dir() and (d / "pseudobulk_delta.npz").exists() and (d / "ctrl_cells.npz").exists():
            out.append(d.name)
    return out


def compute_knockdown_depth_array(
    deltas: np.ndarray,
    target_indices: np.ndarray,
    stored_kd_log1p: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Return `|k_p|` vector (length N_t); NaN where target self-gene is absent."""
    N = deltas.shape[0]
    if stored_kd_log1p is not None and len(stored_kd_log1p) == N:
        kd = np.abs(np.asarray(stored_kd_log1p, dtype=np.float64))
        # Ensure negative target indices stay NaN
        kd = np.where(np.asarray(target_indices) >= 0, kd, np.nan)
        return kd
    kd = np.full(N, np.nan, dtype=np.float64)
    for i, g_idx in enumerate(target_indices):
        if int(g_idx) >= 0:
            kd[i] = abs(float(deltas[i, int(g_idx)]))
    return kd


def load_standardized_line_r10(
    line_dir: Path,
    representation: str = "log1p",
    common_genes: Optional[Sequence[str]] = None,
    hgnc_bridge: Optional[HGNCBridge] = None,
    min_gene_overlap: int = 3000,
) -> LineDatasetR10:
    """
    Load a standardized perturbation dataset from `line_dir` (compatible with both
    Round 9 `round9/p2/shared/lines/<line>` and Round 10 `round10/shared/lines/<line>`).
    Applies read-only canonical HGNC symbol resolution when `hgnc_bridge` is provided.
    """
    if representation not in ("log1p", "zscore", "libnorm"):
        raise ValueError(f"Unsupported representation: {representation}")

    raw_genes = [x.strip() for x in (line_dir / "genes.txt").read_text(encoding="utf-8").splitlines() if x.strip()]
    raw_targets = [x.strip() for x in (line_dir / "targets.txt").read_text(encoding="utf-8").splitlines() if x.strip()]

    pb_npz = np.load(line_dir / "pseudobulk_delta.npz", allow_pickle=False)
    ctrl_npz = np.load(line_dir / "ctrl_cells.npz", allow_pickle=False)
    meta_path = line_dir / "line_metadata.json"
    meta = json.loads(meta_path.read_text(encoding="utf-8")) if meta_path.exists() else {}

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
    repr_mask = np.asarray(pb_npz["reproducible"], dtype=bool) if "reproducible" in pb_npz.files else None
    ad_counts = np.asarray(pb_npz["ad_counts"], dtype=np.float64) if "ad_counts" in pb_npz.files else None
    n_cells_raw = np.asarray(pb_npz["n_cells"], dtype=np.float64) if "n_cells" in pb_npz.files else None

    delta_key = f"delta_{representation}"
    x_key = f"X_{representation}"
    deltas_raw = np.asarray(pb_npz[delta_key], dtype=np.float64)
    ctrl_raw = np.asarray(ctrl_npz[x_key], dtype=np.float32)

    ctrl_mean_log1p = (
        np.asarray(ctrl_npz["ctrl_mean_log1p"], dtype=np.float64)
        if "ctrl_mean_log1p" in ctrl_npz.files
        else np.asarray(ctrl_npz["X_log1p"], dtype=np.float64).mean(axis=0)
    )
    ctrl_mean_libnorm = (
        np.asarray(ctrl_npz["ctrl_mean_libnorm"], dtype=np.float64)
        if "ctrl_mean_libnorm" in ctrl_npz.files
        else np.asarray(ctrl_npz["X_libnorm"], dtype=np.float64).mean(axis=0)
    )
    ctrl_mean_libnorm_full = (
        np.asarray(ctrl_npz["ctrl_mean_libnorm_full"], dtype=np.float64)
        if "ctrl_mean_libnorm_full" in ctrl_npz.files
        else None
    )

    ha_key = f"half_a_delta_{representation}"
    hb_key = f"half_b_delta_{representation}"
    ha_raw = np.asarray(pb_npz[ha_key], dtype=np.float64) if ha_key in pb_npz.files else None
    hb_raw = np.asarray(pb_npz[hb_key], dtype=np.float64) if hb_key in pb_npz.files else None
    hv_raw = np.asarray(pb_npz["half_valid"], dtype=bool) if "half_valid" in pb_npz.files else None
    kd_log1p_raw = (
        np.asarray(pb_npz["knockdown_depth_log1p"], dtype=np.float64)
        if "knockdown_depth_log1p" in pb_npz.files
        else None
    )
    kd_frac_raw = (
        np.asarray(pb_npz["knockdown_fraction"], dtype=np.float64)
        if "knockdown_fraction" in pb_npz.files
        else None
    )

    # Canonical HGNC mapping (read-only over both genes and targets)
    if hgnc_bridge is not None:
        canon_genes, g_keep_idx, _ = hgnc_bridge.dedup_map_indices(
            raw_genes,
            ensembl_ids=None,
            priority_scores=ctrl_mean_libnorm,
            require_approved=False,
        )
        deltas_raw = deltas_raw[:, g_keep_idx]
        ctrl_raw = ctrl_raw[:, g_keep_idx]
        ctrl_mean_log1p = ctrl_mean_log1p[g_keep_idx]
        ctrl_mean_libnorm = ctrl_mean_libnorm[g_keep_idx]
        if ctrl_mean_libnorm_full is not None:
            ctrl_mean_libnorm_full = ctrl_mean_libnorm_full[g_keep_idx]
        if ha_raw is not None:
            ha_raw = ha_raw[:, g_keep_idx]
        if hb_raw is not None:
            hb_raw = hb_raw[:, g_keep_idx]

        # Prioritize targets by reproducible -> ad_counts -> n_cells matching Addendum-B dedup
        t_priority = np.zeros(len(pb_targets), dtype=np.float64)
        if repr_mask is not None:
            t_priority += repr_mask.astype(np.float64) * 1e9
        if ad_counts is not None:
            t_priority += np.nan_to_num(ad_counts, nan=0.0) * 1e4
        if n_cells_raw is not None:
            t_priority += np.nan_to_num(n_cells_raw, nan=0.0)
        canon_targets, t_keep_idx, _ = hgnc_bridge.dedup_map_indices(
            pb_targets,
            ensembl_ids=target_ens,
            priority_scores=t_priority,
            require_approved=False,
        )
        deltas_raw = deltas_raw[t_keep_idx]
        if repr_mask is not None:
            repr_mask = repr_mask[t_keep_idx]
        if n_cells_raw is not None:
            n_cells_raw = n_cells_raw[t_keep_idx]
        if ha_raw is not None:
            ha_raw = ha_raw[t_keep_idx]
        if hb_raw is not None:
            hb_raw = hb_raw[t_keep_idx]
        if hv_raw is not None:
            hv_raw = hv_raw[t_keep_idx]
        if kd_log1p_raw is not None:
            kd_log1p_raw = kd_log1p_raw[t_keep_idx]
        if kd_frac_raw is not None:
            kd_frac_raw = kd_frac_raw[t_keep_idx]
        if target_ens is not None:
            target_ens = [target_ens[int(i)] for i in t_keep_idx]
    else:
        canon_genes = raw_genes
        canon_targets = pb_targets

    if common_genes is not None:
        g_pos = {g: i for i, g in enumerate(canon_genes)}
        aligned_genes = [g for g in common_genes if g in g_pos]
        if len(aligned_genes) < min_gene_overlap:
            raise RuntimeError(
                f"Line {line_dir.name} aligned gene count {len(aligned_genes)} < {min_gene_overlap}."
            )
        col_idx = np.asarray([g_pos[g] for g in aligned_genes], dtype=np.int64)
        genes_final = aligned_genes
        deltas = deltas_raw[:, col_idx]
        ctrl_cells = ctrl_raw[:, col_idx]
        ctrl_mean_log1p = ctrl_mean_log1p[col_idx]
        ctrl_mean_libnorm = ctrl_mean_libnorm[col_idx]
        if ctrl_mean_libnorm_full is not None:
            ctrl_mean_libnorm_full = ctrl_mean_libnorm_full[col_idx]
        ha = ha_raw[:, col_idx] if ha_raw is not None else None
        hb = hb_raw[:, col_idx] if hb_raw is not None else None
    else:
        genes_final = list(canon_genes)
        deltas = deltas_raw
        ctrl_cells = ctrl_raw
        ha = ha_raw
        hb = hb_raw

    gene_to_idx = {g: i for i, g in enumerate(genes_final)}
    target_indices = np.asarray([gene_to_idx.get(t, -1) for t in canon_targets], dtype=np.int64)
    kd_depth = compute_knockdown_depth_array(deltas, target_indices, stored_kd_log1p=kd_log1p_raw)

    return LineDatasetR10(
        name=line_dir.name,
        genes=genes_final,
        target_names=list(canon_targets),
        target_indices=target_indices,
        deltas=deltas,
        ctrl_cells=ctrl_cells,
        ctrl_mean_log1p=ctrl_mean_log1p,
        ctrl_mean_libnorm=ctrl_mean_libnorm,
        ctrl_mean_libnorm_full=ctrl_mean_libnorm_full,
        repr_mask=repr_mask,
        half_a_deltas=ha,
        half_b_deltas=hb,
        half_valid=hv_raw,
        n_cells=n_cells_raw,
        knockdown_depth_log1p=kd_depth,
        knockdown_fraction=kd_frac_raw,
        target_ensembl_ids=target_ens,
        representation=representation,
        metadata=meta,
    )


# ---------------------------------------------------------------------------
# Pure-Python HDF5 / AnnData Stream Loaders & Split-Half Helpers (PREREG_C)
# ---------------------------------------------------------------------------

_DLD1_PREFIX_RE = re.compile(r"^(.+?)_(?:P\d+|ENST\d*|TSS|\d+).*$")


def parse_dld1_guide_gene(guide_name: str) -> Tuple[bool, Optional[str]]:
    """
    Parse DLD-1 (`GSE337988`) sgRNA identifier per `R10_DOWNLOAD_PLAN.md §3.1`.
    Returns `(is_ntc, raw_gene_symbol)`.
    """
    g = guide_name.strip().replace("(W)", "")
    if not g:
        return False, None
    up = g.upper()
    if up.startswith("NTC") or "NON-TARGET" in up or "NONTARGET" in up:
        return True, None
    m = _DLD1_PREFIX_RE.match(g)
    if m:
        return False, m.group(1)
    if "_" in g:
        return False, g.split("_", 1)[0]
    return False, g


def split_units_even_odd(unit_ids: Sequence[str]) -> Tuple[Set[str], Set[str]]:
    """
    Pre-locked `PREREG_C §2.3` deterministic unit split:
    sort distinct biological/experimental unit IDs lexicographically:
      even indices [0, 2, ...] -> Half A
      odd  indices [1, 3, ...] -> Half B
    """
    unique_sorted = sorted(set(str(u) for u in unit_ids if str(u).strip()))
    half_a = {u for idx, u in enumerate(unique_sorted) if idx % 2 == 0}
    half_b = {u for idx, u in enumerate(unique_sorted) if idx % 2 == 1}
    return half_a, half_b


@dataclass
class TenxMatrixBundle:
    """In-memory sparse Cell x Feature CSR matrix read directly from 10x HDF5."""

    barcodes: List[str]
    feature_names: List[str]
    feature_ids: List[str]
    feature_types: List[str]
    X_csr: sp.csr_matrix  # shape (n_cells, n_features)


def read_10x_feature_bc_matrix_h5(
    h5_path: Path,
    feature_type_filter: Optional[str] = None,
) -> TenxMatrixBundle:
    """
    Pure-Python `h5py` reader for 10x CellRanger `.h5` matrices (`matrix/` CSC group).
    Returns cell x feature `scipy.sparse.csr_matrix` without densifying.
    """
    with h5py.File(h5_path, "r") as f:
        if "matrix" not in f:
            raise KeyError(f"{h5_path} does not contain CellRanger `matrix` HDF5 group.")
        m = f["matrix"]
        barcodes = decode_str_array(m["barcodes"][:])
        shape_raw = tuple(int(x) for x in m["shape"][:])  # (n_features, n_cells)
        data = m["data"][:]
        indices = m["indices"][:]
        indptr = m["indptr"][:]
        feat = m["features"]
        f_names = decode_str_array(feat["name"][:])
        f_ids = decode_str_array(feat["id"][:]) if "id" in feat else [""] * len(f_names)
        f_types = (
            decode_str_array(feat["feature_type"][:])
            if "feature_type" in feat
            else ["Gene Expression"] * len(f_names)
        )

    # CellRanger stores CSC (n_features, n_cells); transpose cleanly yields CSR (n_cells, n_features)
    csc = sp.csc_matrix((data, indices, indptr), shape=shape_raw)
    csr = csc.T.tocsr()
    if feature_type_filter is not None:
        keep = [i for i, ft in enumerate(f_types) if ft == feature_type_filter]
        keep_arr = np.asarray(keep, dtype=np.int64)
        csr = csr[:, keep_arr].tocsr()
        f_names = [f_names[i] for i in keep]
        f_ids = [f_ids[i] for i in keep]
        f_types = [f_types[i] for i in keep]

    return TenxMatrixBundle(
        barcodes=barcodes,
        feature_names=f_names,
        feature_ids=f_ids,
        feature_types=f_types,
        X_csr=csr,
    )


def call_single_features_from_umi(
    bundle: TenxMatrixBundle,
    min_umi: int = 1,
) -> Dict[str, Optional[str]]:
    """
    Identify barcodes with strictly 1 assigned feature above `min_umi` (0 or >=2 -> None).
    Used for DLD-1 `.crispr.umi_correct.h5` (single sgRNA) & `.hashing.umi_correct.h5` (singlet).
    """
    csr = bundle.X_csr
    out: Dict[str, Optional[str]] = {}
    for row_i, bc in enumerate(bundle.barcodes):
        start, end = int(csr.indptr[row_i]), int(csr.indptr[row_i + 1])
        cols = csr.indices[start:end]
        vals = csr.data[start:end]
        pos = np.where(vals >= min_umi)[0]
        if len(pos) == 1:
            out[bc] = bundle.feature_names[int(cols[int(pos[0])])]
        else:
            out[bc] = None
    return out


def load_dld1_gsm_channel(
    gex_h5: Path,
    crispr_h5: Path,
    hashing_h5: Path,
    channel_id: str,
    min_crispr_umi: int = 1,
    min_hash_umi: int = 1,
) -> Tuple[TenxMatrixBundle, List[bool], List[Optional[str]], List[str]]:
    """
    Load one DLD-1 (`GSE337988`) Low-MOI 10x channel (`SAM24482910_k`), filtering strictly
    to MULTI-seq singlets AND cells uniquely assigned to exactly 1 sgRNA (`P4`).
    """
    gex = read_10x_feature_bc_matrix_h5(gex_h5, feature_type_filter="Gene Expression")
    crispr_b = read_10x_feature_bc_matrix_h5(crispr_h5, feature_type_filter=None)
    hash_b = read_10x_feature_bc_matrix_h5(hashing_h5, feature_type_filter=None)

    crispr_calls = call_single_features_from_umi(crispr_b, min_umi=min_crispr_umi)
    hash_calls = call_single_features_from_umi(hash_b, min_umi=min_hash_umi)

    keep_rows: List[int] = []
    is_ntc_list: List[bool] = []
    raw_tgt_list: List[Optional[str]] = []
    kept_barcodes: List[str] = []

    for r_i, bc in enumerate(gex.barcodes):
        h_call = hash_calls.get(bc)
        c_call = crispr_calls.get(bc)
        if h_call is None or c_call is None:
            continue
        is_ntc, tgt = parse_dld1_guide_gene(c_call)
        if not is_ntc and tgt is None:
            continue
        keep_rows.append(r_i)
        is_ntc_list.append(is_ntc)
        raw_tgt_list.append(tgt)
        kept_barcodes.append(f"{channel_id}:{bc}")

    keep_arr = np.asarray(keep_rows, dtype=np.int64)
    filtered = TenxMatrixBundle(
        barcodes=kept_barcodes,
        feature_names=gex.feature_names,
        feature_ids=gex.feature_ids,
        feature_types=gex.feature_types,
        X_csr=gex.X_csr[keep_arr].tocsr() if len(keep_arr) > 0 else sp.csr_matrix((0, len(gex.feature_names))),
    )
    return filtered, is_ntc_list, raw_tgt_list, [channel_id] * len(kept_barcodes)


def read_cd4t_flex_h5(
    h5_path: Path,
    sgrna_meta_df: Optional[pd.DataFrame] = None,
    min_guide_umi: int = 3,
) -> Tuple[TenxMatrixBundle, np.ndarray, List[Optional[str]], List[Optional[str]]]:
    """
    Pure-Python `h5py` reader for `GSE314342` CD4+ T 10x Flex single-cell `.h5`.
    Extracts RNA GEX sparse matrix and enforces single-sgRNA cell filtering (`P3`, `P4`),
    mapping guides via `sgrna_library_metadata.suppl_table.csv` when provided.
    """
    with h5py.File(h5_path, "r") as f:
        if "matrix" not in f:
            raise KeyError(f"Unsupported CD4+ T Flex HDF5 layout in {h5_path}: expected `matrix/`.")
        feat_types = decode_str_array(f["matrix"]["features"]["feature_type"][:])

    gex = read_10x_feature_bc_matrix_h5(h5_path, feature_type_filter="Gene Expression")
    guide_types = {ft for ft in feat_types if "crispr" in ft.lower() or "guide" in ft.lower()}
    if guide_types:
        g_bundle = read_10x_feature_bc_matrix_h5(h5_path, feature_type_filter=next(iter(sorted(guide_types))))
    else:
        raise KeyError(f"No CRISPR Guide Capture features found in Flex HDF5 {h5_path}.")

    guide_to_info: Dict[str, Tuple[str, str, str]] = {}
    if sgrna_meta_df is not None:
        id_col = "sgrna_id" if "sgrna_id" in sgrna_meta_df.columns else sgrna_meta_df.columns[0]
        gname_col = "target_gene_name" if "target_gene_name" in sgrna_meta_df.columns else "gene_name"
        gid_col = "target_gene_id" if "target_gene_id" in sgrna_meta_df.columns else "gene_id"
        gtype_col = "guide_type" if "guide_type" in sgrna_meta_df.columns else None
        for r in sgrna_meta_df.itertuples(index=False):
            sd = str(getattr(r, id_col))
            gn = str(getattr(r, gname_col))
            ge = str(getattr(r, gid_col))
            gt = str(getattr(r, gtype_col)) if gtype_col else ("non-targeting" if "NTC" in gn.upper() else "targeting")
            guide_to_info[sd] = (gt.lower(), gn, ge)

    calls = call_single_features_from_umi(g_bundle, min_umi=min_guide_umi)
    keep_rows: List[int] = []
    is_ntc: List[bool] = []
    raw_targets: List[Optional[str]] = []
    raw_ens: List[Optional[str]] = []

    for r_i, bc in enumerate(gex.barcodes):
        gd = calls.get(bc)
        if gd is None:
            continue
        if gd in guide_to_info:
            gt, gn, ge = guide_to_info[gd]
            if gt not in ("targeting", "non-targeting"):
                continue
            nt_flag = gt == "non-targeting"
            keep_rows.append(r_i)
            is_ntc.append(nt_flag)
            raw_targets.append(None if nt_flag else gn)
            raw_ens.append(None if nt_flag else ge)
        else:
            nt_flag = gd.upper().startswith("NTC") or "NON-TARGET" in gd.upper()
            keep_rows.append(r_i)
            is_ntc.append(nt_flag)
            raw_targets.append(None if nt_flag else gd.split("_")[0])
            raw_ens.append(None)

    keep_arr = np.asarray(keep_rows, dtype=np.int64)
    filtered = TenxMatrixBundle(
        barcodes=[gex.barcodes[i] for i in keep_rows],
        feature_names=gex.feature_names,
        feature_ids=gex.feature_ids,
        feature_types=gex.feature_types,
        X_csr=gex.X_csr[keep_arr].tocsr() if len(keep_arr) > 0 else sp.csr_matrix((0, len(gex.feature_names))),
    )
    return filtered, np.asarray(is_ntc, dtype=bool), raw_targets, raw_ens


def iter_h5ad_csr_chunks_r10(
    h5_path: Path,
    row_chunk: int = 4096,
    layer_name: Optional[str] = None,
) -> Iterator[Tuple[int, int, np.ndarray]]:
    """
    Pure-Python `h5py` chunked iterator over `.h5ad` CSR/dense RNA counts (`<= 4 GB` RAM).
    Supports standard AnnData `X` or `layers/<layer_name>`.
    """
    with h5py.File(h5_path, "r") as f:
        node = f["layers"][layer_name] if layer_name else f["X"]
        if isinstance(node, h5py.Dataset):
            n_rows = int(node.shape[0])
            for s in range(0, n_rows, row_chunk):
                e = min(s + row_chunk, n_rows)
                yield s, e, np.asarray(node[s:e, :], dtype=np.float64)
            return

        enc = str(node.attrs.get("encoding-type", "")).lower()
        shape = tuple(int(x) for x in node.attrs["shape"])
        n_rows, n_cols = shape
        if "csc" in enc:
            # Transposed or CSC in rare `.h5ad`: materialize sparse slices once via CSC columns if small
            mat = sp.csc_matrix((node["data"][:], node["indices"][:], node["indptr"][:]), shape=shape).tocsr()
            for s in range(0, n_rows, row_chunk):
                e = min(s + row_chunk, n_rows)
                yield s, e, mat[s:e].toarray().astype(np.float64, copy=False)
            return
        for s in range(0, n_rows, row_chunk):
            e = min(s + row_chunk, n_rows)
            indptr = node["indptr"][s : e + 1]
            p0, p1 = int(indptr[0]), int(indptr[-1])
            data = node["data"][p0:p1]
            indices = node["indices"][p0:p1]
            rel_ptr = indptr - p0
            block = sp.csr_matrix((data, indices, rel_ptr), shape=(e - s, n_cols)).toarray()
            yield s, e, block.astype(np.float64, copy=False)


def count_igvf_single_sgrna_from_metadata(
    guide_meta_csv: Path,
    cell_single_guide_series: Sequence[str],
    r9_source_targets: Set[str],
    hgnc_bridge: HGNCBridge,
    min_cells: int = 25,
) -> Dict[str, Any]:
    """
    Implement `PREREG_C §2.2 P2` pre-download metadata check for `IGVFDS9825FOBP` H9 neurons:
    count shared canonical targets achieving `>= min_cells` single-sgRNA cells without
    touching any expression values (`n_shared_ge25 >= 50` gate).
    """
    meta_df = pd.read_csv(guide_meta_csv)
    g_col = "guide_id" if "guide_id" in meta_df.columns else ("feature_id" if "feature_id" in meta_df.columns else meta_df.columns[0])
    t_col = "intended_target_name" if "intended_target_name" in meta_df.columns else meta_df.columns[1]
    type_col = "type" if "type" in meta_df.columns else None

    guide_to_canon: Dict[str, Optional[str]] = {}
    for row in meta_df.itertuples(index=False):
        gid = str(getattr(row, g_col))
        raw_t = str(getattr(row, t_col))
        gt = str(getattr(row, type_col)).lower() if type_col else "targeting"
        if "non-targeting" in gt or "negative" in gt or raw_t.upper().startswith("NTC"):
            guide_to_canon[gid] = None
            continue
        ens = raw_t if raw_t.startswith("ENSG") else None
        sym = None if raw_t.startswith("ENSG") else raw_t
        guide_to_canon[gid] = hgnc_bridge.canonicalize(sym, ensembl_id=ens)

    counts: Dict[str, int] = {}
    n_ntc = 0
    for gid in cell_single_guide_series:
        if not gid:
            continue
        tgt = guide_to_canon.get(str(gid))
        if tgt is None:
            if str(gid) in guide_to_canon:
                n_ntc += 1
        else:
            counts[tgt] = counts.get(tgt, 0) + 1

    canon_r9 = {c for s in r9_source_targets if (c := hgnc_bridge.canonicalize(s)) is not None}
    shared_ge_min = sorted(t for t, c in counts.items() if c >= min_cells and t in canon_r9)
    return {
        "n_single_sgrna_ntc": n_ntc,
        "n_targets_total": len(counts),
        "n_targets_ge_min_cells": sum(1 for c in counts.values() if c >= min_cells),
        "n_shared_ge_min_cells": len(shared_ge_min),
        "passes_p2_50_target_gate": len(shared_ge_min) >= 50,
    }
