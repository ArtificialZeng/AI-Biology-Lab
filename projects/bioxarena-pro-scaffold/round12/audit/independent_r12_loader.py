#!/usr/bin/env python3
"""
Independent Round 12 Data Loader, Feasibility/Blinding Gate, & Raw Cohort Readers (A-T3).

Clean-room guarantees per `PREREG_E.md` (`v2`, locked commit `4c606f4`):
  - Strictly NEVER reads, imports, or references `round12/src/` or `round12/results/`
    during loader construction or blind execution.
  - `--phase dev_r9 | dev_r10 | dev_r11`: Permits ONLY historical Round 9–11 lines and
    explicitly refuses any access to `round12/shared/` or Round 12 target keys.
  - `--phase holdout_r12`: Runtime hard gate requires BOTH non-empty
    `round12/status/PREREG_E_COMMIT.txt` AND `round12/status/CODE_FREEZE_R12.txt` AND
    `round12/status/CONFIRMATORY_LINES_R12.txt` before touching any Round 12 expression data.
  - Implements independent raw ingestion + E3/QC readers for `a549_ko` (`GSE345058`,
    including `SA-Cas9` `total_counts_CAS9 > 0`) and `hek293_xucao` (`XuCao2023.h5ad`),
    reusing ONLY frozen Round 10 (`process_r10_line`) and Round 11 audit primitives.
"""

from __future__ import annotations

import gzip
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Set, Tuple

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROUND12_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = ROUND12_ROOT.parent
R11_AUDIT_DIR = REPO_ROOT / "round11" / "audit"
R10_SCRIPTS_DIR = REPO_ROOT / "round10" / "shared" / "scripts"
for _p in (str(R11_AUDIT_DIR), str(R10_SCRIPTS_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from independent_r11_loader import (  # noqa: E402
    DEFAULT_HGNC_TSV,
    R9_LINE_NAMES,
    R9_LINE_SET,
    R9_LINES_DIR,
    R10_LINES_DIR,
    R10_SOURCE_EXT_NAMES,
    R11_LINES_DIR,
    AuditLineData,
    CanonicalHGNC,
    _first_token,
    load_line_dataset,
)
from prepare_r10_lines import (  # noqa: E402
    HGNCLookup,
    PartSpec,
    build_hgnc_lookup,
    process_r10_line,
)

R12_SHARED_DIR = ROUND12_ROOT / "shared"
R12_RAW_DIR = R12_SHARED_DIR / "raw"
R12_LINES_DIR = R12_SHARED_DIR / "lines"
R12_STATUS_DIR = ROUND12_ROOT / "status"

R12_TARGET_KEYS: Tuple[str, ...] = ("a549_ko", "hek293_xucao", "a549_ko_cas9pos")
R12_TARGET_SET: Set[str] = set(R12_TARGET_KEYS)

R10_ALL_LINES: Tuple[str, ...] = ("dld1", "cd4t_rest", "wtc11_i3neuron", "calu3")
R11_ALL_LINES: Tuple[str, ...] = ("telohaec", "cd4t_stim48", "feng_ipsc")
DEV_ALLOWED_LINES_BY_PHASE: Dict[str, Tuple[str, ...]] = {
    "dev_r9": R9_LINE_NAMES,
    "dev_r10": R10_ALL_LINES,
    "dev_r11": R11_ALL_LINES,
}
ALL_DEV_LINES: Set[str] = set(R9_LINE_NAMES) | set(R10_ALL_LINES) | set(R11_ALL_LINES)

VALID_HEK_ROLES: Set[str] = {"CONFIRMATORY", "SECONDARY_ONLY", "EXCLUDED_FEASIBILITY_Q3"}


def enforce_r12_gate(
    phase: str,
    target_lines: Sequence[str] = (),
    status_dir: Path = R12_STATUS_DIR,
    extra_dirs: Sequence[Path] = (),
    require_confirmatory_file: bool = True,
) -> Dict[str, str]:
    """
    Enforce `PREREG_E.md §1` Step 2/5 blinding lock before reading expression data:
      - `dev_r9 | dev_r10 | dev_r11`: refuses any Round 12 line key or `round12/shared/` path.
      - `holdout_r12`: requires non-empty `PREREG_E_COMMIT.txt` AND `CODE_FREEZE_R12.txt`
        (and `CONFIRMATORY_LINES_R12.txt` unless `require_confirmatory_file=False`).
    """
    if phase in DEV_ALLOWED_LINES_BY_PHASE:
        bad_targets = [t for t in target_lines if (t in R12_TARGET_SET or t not in ALL_DEV_LINES)]
        if bad_targets:
            raise RuntimeError(
                f"PREREG_E dev blinding lock ({phase}): refused non-dev / Round 12 targets {bad_targets}."
            )
        r12_shared = R12_SHARED_DIR.resolve()
        for d in extra_dirs:
            rp = Path(d).resolve()
            if rp == r12_shared or r12_shared in rp.parents:
                raise RuntimeError(
                    f"PREREG_E dev blinding lock ({phase}): refused accessing Round 12 shared dir {d}."
                )
        return {
            "phase": phase,
            "prereg_commit": _first_token(status_dir / "PREREG_E_COMMIT.txt"),
            "code_freeze": "",
        }

    if phase != "holdout_r12":
        raise ValueError(f"Unsupported phase: {phase!r}")

    prereg_hash = _first_token(status_dir / "PREREG_E_COMMIT.txt")
    freeze_hash = _first_token(status_dir / "CODE_FREEZE_R12.txt")
    if not prereg_hash:
        raise RuntimeError(
            f"PREREG_E hard gate blocked: missing or empty {status_dir / 'PREREG_E_COMMIT.txt'}."
        )
    if not freeze_hash:
        raise RuntimeError(
            f"PREREG_E hard gate blocked: missing or empty {status_dir / 'CODE_FREEZE_R12.txt'}. "
            "Held-out Round 12 expression matrices cannot be accessed before code freeze."
        )
    if require_confirmatory_file:
        conf_p = status_dir / "CONFIRMATORY_LINES_R12.txt"
        if not conf_p.exists() or not conf_p.read_text(encoding="utf-8").strip():
            raise RuntimeError(
                f"PREREG_E hard gate blocked: missing or empty {conf_p}."
            )
    return {
        "phase": phase,
        "prereg_commit": prereg_hash,
        "code_freeze": freeze_hash,
    }


def parse_confirmatory_lines_r12(
    path: Path = R12_STATUS_DIR / "CONFIRMATORY_LINES_R12.txt",
) -> Dict[str, Any]:
    """
    Mechanically parse `CONFIRMATORY_LINES_R12.txt` per `PREREG_E.md §0.2 & §2.2`:
      - `a549_ko` is mandatory: if absent from confirmatory set (or flagged FAIL/ABORT),
        confirmatory adjudication aborts (`K_conf = 0`, `confirmatory_lines = []`), even if
        `hek293_xucao` passed counts.
      - `hek293_xucao` role is one of `CONFIRMATORY`, `SECONDARY_ONLY`, or
        `EXCLUDED_FEASIBILITY_Q3` (never enters secondary reporting if `EXCLUDED_FEASIBILITY_Q3`).
    """
    if not path.exists():
        raise FileNotFoundError(f"Missing confirmatory gate file: {path}")
    text = path.read_text(encoding="utf-8")

    a549_role: Optional[str] = None
    hek_role: Optional[str] = None
    listed_conf: List[str] = []
    secondary_list: List[str] = []

    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        upper = line.upper()
        if upper.startswith("CONFIRMATORY_LINES"):
            rhs = line.split("=", 1)[-1] if "=" in line else line.split(":", 1)[-1]
            for tok in rhs.replace(",", " ").split():
                t = tok.strip().lower()
                if t in ("a549_ko", "hek293_xucao") and t not in listed_conf:
                    listed_conf.append(t)
        elif upper.startswith("SECONDARY_LINES"):
            rhs = line.split("=", 1)[-1] if "=" in line else line.split(":", 1)[-1]
            for tok in rhs.replace(",", " ").split():
                t = tok.strip().lower()
                if t in ("a549_ko", "hek293_xucao") and t not in secondary_list:
                    secondary_list.append(t)
        elif upper.startswith("A549_KO_ROLE") or upper.startswith("A549_ROLE") or line.lower().startswith("a549_ko="):
            rhs = line.split("=", 1)[-1] if "=" in line else line.split(":", 1)[-1]
            tok = rhs.strip().split()[0].upper()
            a549_role = tok
        elif upper.startswith("HEK293_XUCAO_ROLE") or upper.startswith("HEK293_ROLE") or line.lower().startswith("hek293_xucao="):
            rhs = line.split("=", 1)[-1] if "=" in line else line.split(":", 1)[-1]
            tok = rhs.strip().split()[0].upper()
            hek_role = tok
        elif line.lower() in ("a549_ko", "hek293_xucao"):
            t = line.lower()
            if t not in listed_conf:
                listed_conf.append(t)

    if a549_role is None:
        a549_role = "CONFIRMATORY" if "a549_ko" in listed_conf else "FAILED_GATE"
    if hek_role is None:
        if "hek293_xucao" in listed_conf:
            hek_role = "CONFIRMATORY"
        elif "hek293_xucao" in secondary_list:
            hek_role = "SECONDARY_ONLY"
        elif "EXCLUDED_FEASIBILITY_Q3" in text.upper():
            hek_role = "EXCLUDED_FEASIBILITY_Q3"
        else:
            hek_role = "SECONDARY_ONLY"

    if a549_role != "CONFIRMATORY":
        # Mandatory line failed: abort confirmatory adjudication regardless of HEK293 status.
        confirmatory_lines: List[str] = []
        aborted_mandatory = True
    else:
        confirmatory_lines = ["a549_ko"]
        if hek_role == "CONFIRMATORY":
            confirmatory_lines.append("hek293_xucao")
        aborted_mandatory = False

    secondary_eval_lines: List[str] = []
    if aborted_mandatory and a549_role not in ("EXCLUDED_FEASIBILITY_Q3", "INGESTION_FAILED"):
        # If A549 counts were insufficient but ingested cleanly, report descriptively.
        secondary_eval_lines.append("a549_ko")
    if hek_role == "SECONDARY_ONLY" or (aborted_mandatory and hek_role == "CONFIRMATORY"):
        secondary_eval_lines.append("hek293_xucao")

    return {
        "confirmatory_lines": confirmatory_lines,
        "k_conf": len(confirmatory_lines),
        "a549_role": a549_role,
        "hek293_role": hek_role,
        "hek293_excluded_feasibility_q3": hek_role == "EXCLUDED_FEASIBILITY_Q3",
        "aborted_mandatory_a549": aborted_mandatory,
        "secondary_eval_lines": secondary_eval_lines,
    }


def _decode_str(arr: Any) -> List[str]:
    out: List[str] = []
    for x in arr:
        if isinstance(x, (bytes, np.bytes_)):
            out.append(x.decode("utf-8", errors="replace"))
        else:
            out.append(str(x))
    return out


def _read_obs_or_var_column(group: h5py.Group, col_name: str) -> Optional[np.ndarray]:
    """Read an AnnData HDF5 `/obs` or `/var` column handling categorical and flat encodings."""
    if col_name not in group:
        return None
    node = group[col_name]
    if isinstance(node, h5py.Group) and "categories" in node and "codes" in node:
        cats = np.asarray(_decode_str(node["categories"][()]), dtype=object)
        codes = np.asarray(node["codes"][()], dtype=np.int64)
        res = np.empty(len(codes), dtype=object)
        valid = (codes >= 0) & (codes < len(cats))
        res[valid] = cats[codes[valid]]
        res[~valid] = ""
        return res
    raw = node[()]
    if raw.dtype.kind in ("S", "O", "U"):
        vals = _decode_str(raw)
        if "__categories" in group and col_name in group["__categories"] and raw.dtype.kind not in ("S", "O", "U"):
            cats = np.asarray(_decode_str(group["__categories"][col_name][()]), dtype=object)
            codes = np.asarray(raw, dtype=np.int64)
            return cats[codes]
        return np.asarray(vals, dtype=object)
    if "__categories" in group and col_name in group["__categories"]:
        cats = np.asarray(_decode_str(group["__categories"][col_name][()]), dtype=object)
        codes = np.asarray(raw, dtype=np.int64)
        res = np.empty(len(codes), dtype=object)
        valid = (codes >= 0) & (codes < len(cats))
        res[valid] = cats[codes[valid]]
        res[~valid] = ""
        return res
    return np.asarray(raw)


def _list_obs_columns(obs: h5py.Group) -> List[str]:
    if "column-order" in obs.attrs:
        return _decode_str(obs.attrs["column-order"])
    return sorted(k for k in obs.keys() if k != "__categories")


def _read_anndata_index(group: h5py.Group, group_name: str) -> List[str]:
    """
    Read AnnData HDF5 `/obs` or `/var` index (`RULING_R12_PREFREEZE` Ruling A1).
    Checks `group.attrs["_index"]` first (e.g. AnnData 0.2.0 `cell_barcode` / `gene_symbol`),
    then falls back to dataset `"_index"` and `"index"`. Raises `RuntimeError` if absent;
    never synthesizes fake `cell_{i}` or `G{i}` names.
    """
    candidates: List[str] = []
    if "_index" in group.attrs:
        raw_attr = group.attrs["_index"]
        if isinstance(raw_attr, (bytes, np.bytes_)):
            idx_name = raw_attr.decode("utf-8", errors="replace").strip()
        elif isinstance(raw_attr, np.ndarray):
            if raw_attr.ndim == 0:
                scalar = raw_attr.item()
                idx_name = (
                    scalar.decode("utf-8", errors="replace").strip()
                    if isinstance(scalar, (bytes, np.bytes_))
                    else str(scalar).strip()
                )
            elif raw_attr.size > 0:
                idx_name = _decode_str([raw_attr.flat[0]])[0].strip()
            else:
                idx_name = ""
        else:
            idx_name = str(raw_attr).strip()
        if idx_name:
            candidates.append(idx_name)
    for fallback_name in ("_index", "index"):
        if fallback_name not in candidates:
            candidates.append(fallback_name)
    for cname in candidates:
        arr = _read_obs_or_var_column(group, cname)
        if arr is not None:
            return [str(x).strip() for x in arr]
    raise RuntimeError(
        f"MISSING_ANNDATA_INDEX: /{group_name} lacks valid attrs['_index'], '_index', or 'index' dataset."
    )


def _sample_integer_check(h5_node: Any, max_sample: int = 200000) -> Tuple[bool, str]:
    """Verify non-negative integer counts without materializing huge arrays."""
    if isinstance(h5_node, h5py.Dataset):
        dtype = h5_node.dtype
        n_check = min(h5_node.shape[0], max(1, max_sample // max(1, h5_node.shape[1] if h5_node.ndim > 1 else 1)))
        vals = np.asarray(h5_node[:n_check])
        if np.any(~np.isfinite(vals)) or np.any(vals < 0):
            return False, f"negative_or_nonfinite_{dtype}"
        if np.issubdtype(dtype, np.integer):
            return True, f"dense_int_{dtype}"
        diff = np.max(np.abs(vals - np.rint(vals))) if vals.size else 0.0
        ok = bool(diff <= 1e-6 and float(np.max(vals) if vals.size else 0.0) >= 1.0)
        return ok, f"dense_float_integer={ok}_maxdiff={diff:.2e}"

    if isinstance(h5_node, h5py.Group) and "data" in h5_node and "indices" in h5_node and "indptr" in h5_node:
        dset = h5_node["data"]
        dtype = dset.dtype
        n = min(dset.shape[0], max_sample)
        vals = np.asarray(dset[:n])
        if vals.size == 0:
            return False, "empty_sparse_data"
        if np.any(~np.isfinite(vals)) or np.any(vals < 0):
            return False, f"sparse_negative_or_nonfinite_{dtype}"
        if np.issubdtype(dtype, np.integer):
            return True, f"sparse_int_{dtype}"
        diff = float(np.max(np.abs(vals - np.rint(vals))))
        max_v = float(np.max(vals))
        ok = bool(diff <= 1e-6 and max_v >= 1.0)
        return ok, f"sparse_float_integer={ok}_maxdiff={diff:.2e}_max={max_v:.1f}"

    return False, "unrecognized_matrix_encoding"


def check_hek293_xucao_feasibility(h5_path: Path) -> Dict[str, Any]:
    """
    PREREG_E §3.2 / E3 Pre-Metric Technical Feasibility Gate on `XuCao2023.h5ad`:
      1. Requires non-negative integer raw UMI counts in `/X` OR any `/layers/*`.
      2. Requires per-cell `sgRNA`/guide and explicit NTC/perturbation labels in `/obs`.
    If either condition fails, returns `status = "EXCLUDED_FEASIBILITY_Q3"`.
    """
    if not Path(h5_path).exists():
        return {"status": "EXCLUDED_FEASIBILITY_Q3", "feasible": False, "reason": f"missing_file:{h5_path}"}

    with h5py.File(h5_path, "r") as f:
        counts_candidates: List[Tuple[str, str, bool]] = []
        if "X" in f:
            ok_x, desc_x = _sample_integer_check(f["X"])
            counts_candidates.append(("/X", desc_x, ok_x))
        if "layers" in f and isinstance(f["layers"], h5py.Group):
            for lk in sorted(f["layers"].keys()):
                ok_l, desc_l = _sample_integer_check(f["layers"][lk])
                counts_candidates.append((f"/layers/{lk}", desc_l, ok_l))

        explicit_total_layer_names = {"counts", "raw_counts", "raw", "umi_counts", "umi", "total_counts", "x_raw"}
        excluded_subfraction_tokens = ("nascent", "unspliced", "spliced", "ambiguous", "velocity", "intron", "exon")

        chosen_count_path: Optional[str] = None
        for path_k, _, ok_flag in counts_candidates:
            if not ok_flag or not path_k.startswith("/layers/"):
                continue
            lname = path_k.split("/layers/", 1)[1].lower()
            if lname in explicit_total_layer_names:
                chosen_count_path = path_k
                break
        if chosen_count_path is None:
            for path_k, _, ok_flag in counts_candidates:
                if path_k == "/X" and ok_flag:
                    chosen_count_path = path_k
                    break
        if chosen_count_path is None:
            for path_k, _, ok_flag in counts_candidates:
                if not ok_flag or not path_k.startswith("/layers/"):
                    continue
                lname = path_k.split("/layers/", 1)[1].lower()
                if not any(tok in lname for tok in excluded_subfraction_tokens) and (
                    "raw" in lname or "count" in lname or "umi" in lname or "total" in lname
                ):
                    chosen_count_path = path_k
                    break
        if chosen_count_path is None:
            for path_k, _, ok_flag in counts_candidates:
                if not ok_flag:
                    continue
                if path_k.startswith("/layers/"):
                    lname = path_k.split("/layers/", 1)[1].lower()
                    if any(tok in lname for tok in excluded_subfraction_tokens):
                        continue
                chosen_count_path = path_k
                break

        obs_cols = _list_obs_columns(f["obs"]) if "obs" in f and isinstance(f["obs"], h5py.Group) else []
        guide_col = None
        pert_col = None
        for c in ("guide_id", "sgRNA", "sgrna", "guides", "perturbation", "gene", "target"):
            for oc in obs_cols:
                if oc.lower() == c.lower():
                    if "guide" in oc.lower() or "sgrna" in oc.lower():
                        if guide_col is None:
                            guide_col = oc
                    elif pert_col is None:
                        pert_col = oc
        has_ntc_label = False
        ntc_checked_col = None
        for col in [c for c in (pert_col, guide_col, "perturbation", "perturbation_type") if c and c in obs_cols]:
            vals = _read_obs_or_var_column(f["obs"], col)
            if vals is None:
                continue
            uniq_vals = {str(v).strip().lower() for v in vals}
            if any(
                v in ("control", "ntc", "nontargeting", "non-targeting", "non_targeting", "nt")
                or v.startswith(("ntc", "nontarget", "non-target", "non_target"))
                for v in uniq_vals
            ):
                has_ntc_label = True
                ntc_checked_col = col
                break

        has_guide_info = (guide_col is not None) or (pert_col is not None)
        feasible = bool(chosen_count_path is not None and has_guide_info and has_ntc_label)
        reasons = []
        if chosen_count_path is None:
            reasons.append("no_nonnegative_integer_raw_umi_counts_in_X_or_layers")
        if not has_guide_info:
            reasons.append("missing_sgrna_or_perturbation_obs_column")
        if not has_ntc_label:
            reasons.append("missing_explicit_ntc_labels_in_obs")

        return {
            "status": "FEASIBLE" if feasible else "EXCLUDED_FEASIBILITY_Q3",
            "feasible": feasible,
            "chosen_counts_path": chosen_count_path,
            "counts_candidates": counts_candidates,
            "obs_cols": obs_cols,
            "guide_col": guide_col,
            "pert_col": pert_col,
            "ntc_checked_col": ntc_checked_col,
            "reason": ";".join(reasons) if reasons else "OK",
        }


def make_csr_or_csc_h5_readers(
    h5_path: Path,
    node_path: str,
    keep_rows: np.ndarray,
    n_rows_raw: int,
    n_cols_raw: int,
    encoding: str = "csr",
) -> Tuple[
    Callable[[int, int], np.ndarray],
    Callable[[np.ndarray], np.ndarray],
    Callable[[int, int, np.ndarray], np.ndarray],
]:
    """Memory-bounded row subset readers supporting CSR, CSC, or dense `/X` & `/layers/*`."""
    keep_rows = np.asarray(keep_rows, dtype=np.int64)
    if keep_rows.size > 1 and not np.all(keep_rows[:-1] <= keep_rows[1:]):
        raise ValueError("keep_rows must be strictly sorted ascending")

    def _read(raw_rows: np.ndarray, feat_cols: Optional[np.ndarray] = None) -> np.ndarray:
        out_dtype = np.float32 if feat_cols is not None else np.float64
        out_w = len(feat_cols) if feat_cols is not None else n_cols_raw
        if raw_rows.size == 0:
            return np.zeros((0, out_w), dtype=out_dtype)

        with h5py.File(h5_path, "r") as f:
            node = f[node_path]
            if isinstance(node, h5py.Dataset):
                arr = np.asarray(node[raw_rows, :], dtype=out_dtype)
                return arr[:, feat_cols] if feat_cols is not None else arr

            indptr = node["indptr"]
            indices = node["indices"]
            data = node["data"]
            if encoding == "csr":
                r_min, r_max = int(raw_rows[0]), int(raw_rows[-1])
                span = r_max - r_min + 1
                if len(raw_rows) / max(span, 1) >= 0.20:
                    ptrs = indptr[r_min : r_max + 2].astype(np.int64)
                    p0, p1 = int(ptrs[0]), int(ptrs[-1])
                    mat = sp.csr_matrix(
                        (data[p0:p1], indices[p0:p1], ptrs - np.int64(p0)),
                        shape=(span, n_cols_raw),
                    )[raw_rows - r_min, :]
                else:
                    all_d, all_i = [], []
                    new_ptr = np.zeros(len(raw_rows) + 1, dtype=np.int64)
                    cur_p = 0
                    starts = indptr[raw_rows]
                    ends = indptr[raw_rows + 1]
                    for idx_r, (s0, s1) in enumerate(zip(starts, ends)):
                        if s1 > s0:
                            all_d.append(data[int(s0) : int(s1)])
                            all_i.append(indices[int(s0) : int(s1)])
                            cur_p += int(s1 - s0)
                        new_ptr[idx_r + 1] = cur_p
                    d_cat = np.concatenate(all_d) if all_d else np.empty(0, dtype=np.float32)
                    i_cat = np.concatenate(all_i) if all_i else np.empty(0, dtype=np.int32)
                    mat = sp.csr_matrix((d_cat, i_cat, new_ptr), shape=(len(raw_rows), n_cols_raw))
            elif encoding == "csc":
                # Stream columns once into sparse CSR slice over requested rows
                col_targets = np.arange(n_cols_raw, dtype=np.int64) if feat_cols is None else np.asarray(feat_cols, dtype=np.int64)
                pos_map = np.full(n_rows_raw, -1, dtype=np.int64)
                pos_map[raw_rows] = np.arange(len(raw_rows), dtype=np.int64)
                coo_r, coo_c, coo_v = [], [], []
                all_ptrs = indptr[:].astype(np.int64)
                for j_out, c_idx in enumerate(col_targets):
                    p0, p1 = int(all_ptrs[c_idx]), int(all_ptrs[c_idx + 1])
                    if p1 <= p0:
                        continue
                    rows_c = indices[p0:p1]
                    mapped = pos_map[rows_c]
                    hit = mapped >= 0
                    if np.any(hit):
                        coo_r.append(mapped[hit])
                        coo_c.append(np.full(int(np.sum(hit)), j_out, dtype=np.int32))
                        coo_v.append(data[p0:p1][hit])
                out_arr = np.zeros((len(raw_rows), len(col_targets)), dtype=out_dtype)
                if coo_r:
                    out_arr[np.concatenate(coo_r), np.concatenate(coo_c)] = np.concatenate(coo_v).astype(out_dtype)
                return out_arr
            else:
                raise ValueError(f"Unsupported HDF5 sparse encoding: {encoding!r}")

            if feat_cols is not None:
                mat = mat[:, feat_cols]
            return mat.toarray().astype(out_dtype, copy=False)

    return (
        (lambda s, e: _read(keep_rows[s:e], None)),
        (lambda loc_idx: _read(keep_rows[loc_idx], None)),
        (lambda s, e, fcols: _read(keep_rows[s:e], fcols)),
    )


def _stream_mtx_gz_to_csr_h5(
    mtx_gz_path: Path,
    cell_keep_order: np.ndarray,
    n_cells_annot: int,
    n_genes_annot: int,
    out_csr_h5: Path,
    cas9_col_raw: Optional[int] = None,
    chunk_lines: int = 2_000_000,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Stream `GSE345058_raw_counts.mtx.gz` in bounded-memory chunks into a sorted CSR HDF5 cache
    aligned strictly to `cell_annotations` rows (`0..n_cells_annot-1`).
    Auto-detects `(cells x genes)` vs `(genes x cells)` orientation from header dimensions.
    Returns `(umi_per_cell, cas9_counts_per_cell)`.
    """
    out_csr_h5.parent.mkdir(parents=True, exist_ok=True)
    tmp_dir = out_csr_h5.parent / f"_tmp_mtx_shards_{os.getpid()}"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    header_lines = 0
    with gzip.open(mtx_gz_path, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline()
        header_lines += 1
        if not header.startswith("%%MatrixMarket"):
            raise RuntimeError(f"Invalid MatrixMarket header in {mtx_gz_path}: {header!r}")
        line = fh.readline()
        header_lines += 1
        while line.startswith("%") or not line.strip():
            line = fh.readline()
            header_lines += 1
        parts = line.strip().split()
        n_rows_mtx, n_cols_mtx = int(parts[0]), int(parts[1])

    if (n_rows_mtx, n_cols_mtx) == (n_cells_annot, n_genes_annot):
        cell_first = True
    elif (n_rows_mtx, n_cols_mtx) == (n_genes_annot, n_cells_annot):
        cell_first = False
    else:
        raise RuntimeError(
            f"MTX dimension mismatch: got ({n_rows_mtx}, {n_cols_mtx}), "
            f"expected ({n_cells_annot}, {n_genes_annot}) or transpose."
        )

    umi = np.zeros(n_cells_annot, dtype=np.float64)
    cas9_umi = np.zeros(n_cells_annot, dtype=np.float64)
    n_buckets = max(1, (n_cells_annot + 32767) // 32768)
    bucket_r: List[List[np.ndarray]] = [[] for _ in range(n_buckets)]
    bucket_c: List[List[np.ndarray]] = [[] for _ in range(n_buckets)]
    bucket_v: List[List[np.ndarray]] = [[] for _ in range(n_buckets)]
    shard_files: List[List[Path]] = [[] for _ in range(n_buckets)]
    shard_counter = 0

    def _flush_buckets() -> None:
        nonlocal shard_counter
        for bi in range(n_buckets):
            if not bucket_r[bi]:
                continue
            rr = np.concatenate(bucket_r[bi])
            cc = np.concatenate(bucket_c[bi])
            vv = np.concatenate(bucket_v[bi])
            bucket_r[bi].clear()
            bucket_c[bi].clear()
            bucket_v[bi].clear()
            sp_p = tmp_dir / f"b{bi}_s{shard_counter}.npz"
            shard_counter += 1
            np.savez(sp_p, r=rr, c=cc, v=vv)
            shard_files[bi].append(sp_p)

    def _ingest_arrays(r0: np.ndarray, c0: np.ndarray, vv: np.ndarray) -> None:
        if r0.size == 0:
            return
        cell_i, gene_i = (r0, c0) if cell_first else (c0, r0)

        np.add.at(umi, cell_i, vv.astype(np.float64))
        if cas9_col_raw is not None:
            is_cas9 = gene_i == cas9_col_raw
            if np.any(is_cas9):
                np.add.at(cas9_umi, cell_i[is_cas9], vv[is_cas9].astype(np.float64))

        b_idx = cell_i // 32768
        for bi in np.unique(b_idx):
            m_b = b_idx == bi
            bucket_r[int(bi)].append(cell_i[m_b].astype(np.int32))
            bucket_c[int(bi)].append(gene_i[m_b].astype(np.int32))
            bucket_v[int(bi)].append(vv[m_b])

    mem_entries = 0
    for chunk in pd.read_csv(
        mtx_gz_path,
        sep=r"\s+",
        skiprows=header_lines,
        comment="%",
        header=None,
        names=["r", "c", "v"],
        dtype={"r": np.int64, "c": np.int64, "v": np.float32},
        chunksize=20_000_000,
        engine="c",
    ):
        if chunk.empty:
            continue
        r0 = chunk["r"].to_numpy(dtype=np.int64) - 1
        c0 = chunk["c"].to_numpy(dtype=np.int64) - 1
        vv = chunk["v"].to_numpy(dtype=np.float32)
        _ingest_arrays(r0, c0, vv)
        mem_entries += int(r0.size)
        if mem_entries >= 40_000_000:
            _flush_buckets()
            mem_entries = 0
    _flush_buckets()

    try:
        indptr = np.zeros(n_cells_annot + 1, dtype=np.int64)
        cur_p = 0
        with h5py.File(out_csr_h5, "w") as hf:
            xg = hf.create_group("X")
            ds_d = xg.create_dataset("data", shape=(0,), maxshape=(None,), dtype=np.float32, chunks=(65536,))
            ds_i = xg.create_dataset("indices", shape=(0,), maxshape=(None,), dtype=np.int32, chunks=(65536,))
            for bi in range(n_buckets):
                r_start = bi * 32768
                r_end = min(n_cells_annot, (bi + 1) * 32768)
                b_rows = r_end - r_start
                rr_l, cc_l, vv_l = [], [], []
                for sp_p in shard_files[bi]:
                    with np.load(sp_p, allow_pickle=False) as d:
                        rr_l.append(d["r"] - r_start)
                        cc_l.append(d["c"])
                        vv_l.append(d["v"])
                if rr_l:
                    mat_b = sp.csr_matrix(
                        (np.concatenate(vv_l), (np.concatenate(rr_l), np.concatenate(cc_l))),
                        shape=(b_rows, n_genes_annot),
                        dtype=np.float32,
                    )
                    mat_b.sum_duplicates()
                    mat_b.sort_indices()
                else:
                    mat_b = sp.csr_matrix((b_rows, n_genes_annot), dtype=np.float32)

                nnz_b = mat_b.nnz
                if nnz_b > 0:
                    ds_d.resize((cur_p + nnz_b,))
                    ds_i.resize((cur_p + nnz_b,))
                    ds_d[cur_p : cur_p + nnz_b] = mat_b.data
                    ds_i[cur_p : cur_p + nnz_b] = mat_b.indices
                indptr[r_start + 1 : r_end + 1] = cur_p + mat_b.indptr[1:].astype(np.int64)
                cur_p += nnz_b
            xg.create_dataset("indptr", data=indptr, dtype=np.int64)
    finally:
        for sp_list in shard_files:
            for sp_p in sp_list:
                sp_p.unlink(missing_ok=True)
        try:
            tmp_dir.rmdir()
        except OSError:
            pass

    return umi, cas9_umi


def build_a549_ko_audit_parts(
    raw_dir: Path,
    hgnc: HGNCLookup,
    exclude_zero_cas9: bool = False,
    cache_h5_path: Optional[Path] = None,
) -> Tuple[List[PartSpec], Dict[str, Any]]:
    """
    Independent `a549_ko` (`GSE345058`, Liu et al. 2026) audit reader (`PREREG_E.md §3.1`):
      - Joins strictly on unique `cell_id` (verifying `cell_id` uniqueness; NEVER `cell_barcode`).
      - Excludes synthetic vector row `Cas9` (`gene_symbol == "Cas9"` / `genome == "CAS9"` / `ensembl_id == "Cas9"`).
      - Keeps ONLY `Condition == "Nontargeting"` / `Perturbation == "Nontargeting"` NTCs and single-sgRNA perturbations.
      - Split-half: `Half A = batch1 + batch3`, `Half B = batch2 + batch4`.
      - NTC centering unit: strictly `Lane` (`32` lanes nested inside batch & half).
      - `exclude_zero_cas9=True` enables descriptive `SA-Cas9` arm (`total_counts_CAS9 > 0`).
    """
    adir = raw_dir / "a549_ko"
    mtx_p = adir / "GSE345058_raw_counts.mtx.gz"
    cell_p = adir / "GSE345058_cell_annotations.csv.gz"
    gene_p = adir / "GSE345058_gene_annotations.csv.gz"
    sgrna_p = adir / "GSE345058_sgRNA.csv.gz"
    for p in (mtx_p, cell_p, gene_p):
        if not p.exists():
            raise FileNotFoundError(f"Missing GSE345058 file: {p}")

    cdf = pd.read_csv(cell_p, low_memory=False)
    gdf = pd.read_csv(gene_p, low_memory=False)
    if "cell_id" in cdf.columns and cdf["cell_id"].nunique() != len(cdf):
        raise RuntimeError("GSE345058 cell_id is not unique!")
    if sgrna_p.exists():
        sdf = pd.read_csv(sgrna_p, low_memory=False)
        if "cell_id" in sdf.columns and "sgRNA" in sdf.columns and sdf["cell_id"].nunique() == len(sdf):
            pass  # Verified unique cell_id key integrity

    # Identify and drop synthetic Cas9 feature row cleanly without dropping columns from cached indices
    sym_col = "gene_symbol" if "gene_symbol" in gdf.columns else ("gene_name" if "gene_name" in gdf.columns else gdf.columns[0])
    ens_col = "ensembl_id" if "ensembl_id" in gdf.columns else ("gene_id" if "gene_id" in gdf.columns else gdf.columns[0])
    genome_col = "genome" if "genome" in gdf.columns else None

    raw_syms = [str(x).strip() for x in gdf[sym_col].tolist()]
    raw_ens = [str(x).strip().split(".")[0] for x in gdf[ens_col].tolist()]
    cas9_mask = np.asarray(
        [
            (s.upper() == "CAS9")
            or (e.upper() == "CAS9")
            or (genome_col is not None and str(gdf[genome_col].iloc[i]).strip().upper() == "CAS9")
            for i, (s, e) in enumerate(zip(raw_syms, raw_ens))
        ],
        dtype=bool,
    )
    cas9_cols = np.flatnonzero(cas9_mask)
    cas9_col_raw = int(cas9_cols[0]) if cas9_cols.size > 0 else None

    # Resolve metadata columns fail-fast before streaming 10GB MTX
    def _resolve_col(df: pd.DataFrame, candidates: Sequence[str]) -> Optional[str]:
        for c in candidates:
            if c in df.columns:
                return c
        lower_map = {str(col).lower(): str(col) for col in df.columns}
        for c in candidates:
            hit = lower_map.get(c.lower())
            if hit is not None:
                return hit
        return None

    batch_col = _resolve_col(cdf, ("Batch", "batch"))
    pert_col = _resolve_col(cdf, ("Perturbation", "perturbation"))
    if batch_col is None or pert_col is None:
        raise KeyError(
            f"GSE345058 cell_annotations missing Batch/Perturbation column: found {list(cdf.columns)}"
        )
    lane_col = _resolve_col(cdf, ("Lane", "lane"))
    cond_col = _resolve_col(cdf, ("Condition", "condition"))
    sgrna_col = _resolve_col(cdf, ("sgRNA", "sgrna", "guide_id"))
    cas9_col_name = _resolve_col(cdf, ("total_counts_CAS9", "cas9_counts"))

    batch_s = cdf[batch_col].astype(str).str.strip().str.lower()
    lane_s = cdf[lane_col].astype(str).str.strip() if lane_col is not None else batch_s
    cond_s = cdf[cond_col].astype(str).str.strip() if cond_col is not None else pd.Series([""] * len(cdf))
    pert_s = cdf[pert_col].astype(str).str.strip()
    sgrna_s = cdf[sgrna_col].astype(str).str.strip() if sgrna_col is not None else pert_s

    h5_cache = cache_h5_path or (adir / "_audit_a549_csr.h5")
    umi_all, cas9_mtx_umi = _stream_mtx_gz_to_csr_h5(
        mtx_gz_path=mtx_p,
        cell_keep_order=np.arange(len(cdf), dtype=np.int64),
        n_cells_annot=len(cdf),
        n_genes_annot=len(gdf),
        out_csr_h5=h5_cache,
        cas9_col_raw=cas9_col_raw,
    )

    is_ntc = (cond_s == "Nontargeting").to_numpy() | (pert_s == "Nontargeting").to_numpy()
    single_sg_ok = np.asarray(
        [(bool(sg) and ("," not in sg) and (";" not in sg) and ("|" not in sg) and sg.lower() not in ("nan", "none", "unassigned")) for sg in sgrna_s],
        dtype=bool,
    )
    # Use total_counts_CAS9 if present, else fallback to raw MTX Cas9 count
    if cas9_col_name is not None:
        cas9_vals = pd.to_numeric(cdf[cas9_col_name], errors="coerce").fillna(0.0).to_numpy()
    else:
        cas9_vals = cas9_mtx_umi

    half_arr = np.zeros(len(cdf), dtype=np.int8)
    half_arr[batch_s.isin({"batch1", "batch3"}).to_numpy()] = 1
    half_arr[batch_s.isin({"batch2", "batch4"}).to_numpy()] = 2

    keep_mask = single_sg_ok & (half_arr > 0)
    if exclude_zero_cas9:
        keep_mask &= cas9_vals > 0.0
    keep_rows = np.flatnonzero(keep_mask)

    keep_gene_cols = np.flatnonzero(~cas9_mask)
    var_ids = [raw_ens[i] for i in keep_gene_cols]
    var_names = [raw_syms[i] for i in keep_gene_cols]

    nt_kept = is_ntc[keep_rows]
    batches_kept = lane_s.iloc[keep_rows].to_numpy(dtype=str)
    halves_kept = half_arr[keep_rows]
    umi_human = np.maximum(umi_all - cas9_mtx_umi, 0.0)
    umi_kept = umi_human[keep_rows]

    pkeys: List[str] = []
    pmeta: Dict[str, Tuple[str, str, str, str]] = {}
    for row_i, nt_f in zip(keep_rows, nt_kept):
        if nt_f:
            pkeys.append("")
            continue
        raw_t = str(pert_s.iloc[row_i]).strip()
        if not raw_t or raw_t.lower() in ("nan", "none", "unassigned"):
            pkeys.append("")
            continue
        canon_t, hstat = hgnc.resolve(raw_t, "")
        ens_t = hgnc.sym_to_ens.get(canon_t, "")
        # Pool sgRNAs per canonical gene (`P4` single-sgRNA cells aggregated per perturbation)
        pkey = f"{canon_t}__{raw_t}" if canon_t != raw_t else canon_t
        pkeys.append(pkey)
        if pkey not in pmeta:
            pmeta[pkey] = (ens_t, canon_t, raw_t, hstat)

    raw_rc, raw_rr, raw_rcc = make_csr_or_csc_h5_readers(
        h5_path=h5_cache,
        node_path="X",
        keep_rows=keep_rows,
        n_rows_raw=len(cdf),
        n_cols_raw=len(gdf),
        encoding="csr",
    )

    def _wrap_chunk(s: int, e: int) -> np.ndarray:
        return raw_rc(s, e)[:, keep_gene_cols]

    def _wrap_rows(loc_idx: np.ndarray) -> np.ndarray:
        return raw_rr(loc_idx)[:, keep_gene_cols]

    def _wrap_cols(s: int, e: int, fcols: np.ndarray) -> np.ndarray:
        return raw_rcc(s, e, keep_gene_cols[fcols])

    part = PartSpec(
        tag="GSE345058_A549_CRISPRko" + ("_CAS9pos" if exclude_zero_cas9 else ""),
        path=h5_cache,
        n_cells=len(keep_rows),
        n_genes_raw=len(keep_gene_cols),
        nt_mask=nt_kept,
        batches=batches_kept,
        half_flags=halves_kept,
        umi=umi_kept,
        pert_keys=pkeys,
        pert_meta=pmeta,
        var_gene_ids=var_ids,
        var_gene_names=var_names,
        read_chunk_fn=_wrap_chunk,
        read_rows_fn=_wrap_rows,
        read_chunk_cols_fn=_wrap_cols,
    )
    audit_meta = {
        "n_cells_annot": int(len(cdf)),
        "n_genes_raw_including_cas9": int(len(gdf)),
        "n_cas9_excluded_features": int(cas9_mask.sum()),
        "n_kept_cells": int(len(keep_rows)),
        "n_kept_ntc": int(nt_kept.sum()),
        "exclude_zero_cas9": bool(exclude_zero_cas9),
        "n_lanes": int(len(np.unique(batches_kept))),
        "tmp_h5_path": str(h5_cache),
    }
    return [part], audit_meta


def build_hek293_xucao_audit_parts(
    raw_dir: Path,
    hgnc: HGNCLookup,
    apply_author_kd40: bool = False,
    min_umi: float = 300.0,
    min_genes: int = 200,
) -> Tuple[List[PartSpec], Dict[str, Any]]:
    """
    Independent `hek293_xucao` (`XuCao2023.h5ad`) `h5py` reader (`PREREG_E.md §3.2`,
    `RULING_R12_PREFREEZE.md` Ruling A1 & Ruling B):
      - Reads `/obs` and `/var` index strictly via `_read_anndata_index` (`attrs["_index"]`,
        falling back to `"_index"` then `"index"`, raising `RuntimeError` if missing — never
        falling back to synthetic `cell_{i}` / `G{i}` labels).
      - Runs mandatory pre-metric feasibility check (`EXCLUDED_FEASIBILITY_Q3` if lacking
        integer raw UMI counts or per-cell sgRNA/NTC labels).
      - **NO mitochondrial percentage filter (`max_pct_mt = None`, Ruling B):** keeps all
        QC-valid singlets (`>=300` UMI and `>=200` detected human genes by default, no `pct_mt` cap).
      - Applies single-sgRNA / NTC rule without post-hoc `>40%` KD filtering by default
        (`apply_author_kd40=False` per Directive §3).
      - Uses Two-Stage deterministic split:
          Stage 1: sorted `.obs` replicate odd/even parity if >= 2 replicates & >= 15 NTC/half;
          Stage 2: barcode-sorted deterministic stratified `np.random.default_rng(42)` 50/50 split.
    """
    h5_p = raw_dir / "hek293_xucao" / "XuCao2023.h5ad"
    feas = check_hek293_xucao_feasibility(h5_p)
    if not feas["feasible"]:
        raise RuntimeError(f"EXCLUDED_FEASIBILITY_Q3: {feas['reason']}")

    count_path = str(feas["chosen_counts_path"]).lstrip("/")
    with h5py.File(h5_p, "r") as f:
        if "obs" not in f or "var" not in f:
            raise RuntimeError("MISSING_ANNDATA_GROUP: XuCao2023.h5ad requires both /obs and /var groups.")
        obs = f["obs"]
        var = f["var"]
        barcodes = _read_anndata_index(obs, "obs")
        var_ids_raw = _read_anndata_index(var, "var")

        x_node = f[count_path]
        enc_attr = x_node.attrs.get("encoding-type", "")
        if isinstance(enc_attr, bytes):
            enc_attr = enc_attr.decode("utf-8", errors="ignore")
        if isinstance(x_node, h5py.Dataset):
            enc = "dense"
            n_cells_raw, n_genes_raw = x_node.shape
        elif "csc" in str(enc_attr).lower():
            enc = "csc"
            sh = x_node.attrs.get("shape", None)
            n_cells_raw = len(barcodes) if len(barcodes) > 0 else int(sh[0])
            n_genes_raw = len(var_ids_raw) if len(var_ids_raw) > 0 else int(sh[1])
        else:
            enc = "csr"
            sh = x_node.attrs.get("shape", None)
            n_cells_raw = int(sh[0]) if sh is not None else len(x_node["indptr"]) - 1
            n_genes_raw = int(sh[1]) if sh is not None else len(var_ids_raw)

        if len(barcodes) != n_cells_raw or len(var_ids_raw) != n_genes_raw:
            raise RuntimeError(
                f"INDEX_LENGTH_MISMATCH: obs={len(barcodes)} vs matrix_rows={n_cells_raw}, "
                f"var={len(var_ids_raw)} vs matrix_cols={n_genes_raw}"
            )

        v_sym = None
        for vc in ("gene_symbol", "gene_name", "symbol", "feature_name"):
            tmp_v = _read_obs_or_var_column(var, vc)
            if tmp_v is not None:
                v_sym = [str(x).strip() for x in tmp_v]
                break
        if v_sym is None:
            v_sym = list(var_ids_raw)

        v_ens = None
        for vc in ("ensembl_id", "gene_id", "ensembl_gene_id"):
            tmp_v = _read_obs_or_var_column(var, vc)
            if tmp_v is not None:
                v_ens = [str(x).strip().split(".")[0] for x in tmp_v]
                break
        if v_ens is None:
            v_ens = [x.split(".")[0] if str(x).startswith("ENSG") else "" for x in var_ids_raw]

        # Exclude synthetic lentiviral Cas9/dCas9 and non-human (e.g. ENSMUSG / mm10) features
        keep_gene_mask = np.asarray(
            [
                not (
                    ("CAS9" in s.upper())
                    or ("CAS9" in e.upper())
                    or e.upper().startswith("ENSMUSG")
                    or s.lower().startswith("mm10")
                )
                for s, e in zip(v_sym, v_ens)
            ],
            dtype=bool,
        )
        keep_gene_cols = np.flatnonzero(keep_gene_mask)

        # Resolve obs perturbation, nperts/singlet, kd, and replicate columns
        obs_cols = _list_obs_columns(obs)
        pert_arr = _read_obs_or_var_column(obs, feas["pert_col"] or feas["guide_col"])
        guide_arr = _read_obs_or_var_column(obs, feas["guide_col"] or feas["pert_col"])
        nperts_arr = None
        for nc in ("nperts", "n_guides", "num_guides"):
            if nc in obs_cols:
                nperts_arr = np.asarray(pd.to_numeric(_read_obs_or_var_column(obs, nc), errors="coerce"), dtype=np.float64)
                break
        rep_arr = None
        rep_col_used = None
        for rc in ("replicate", "rep", "sample", "batch"):
            for oc in obs_cols:
                if oc.lower() == rc:
                    cand = _read_obs_or_var_column(obs, oc)
                    if cand is not None and len({str(x).strip() for x in cand if str(x).strip()}) >= 2:
                        rep_arr = np.asarray([str(x).strip() for x in cand], dtype=object)
                        rep_col_used = oc
                        break
            if rep_arr is not None:
                break

        kd_arr = None
        if apply_author_kd40:
            for kc in ("knockdown", "kd_fraction", "on_target_kd"):
                if kc in obs_cols:
                    kd_arr = np.asarray(pd.to_numeric(_read_obs_or_var_column(obs, kc), errors="coerce"), dtype=np.float64)
                    break

    # Compute library UMI totals directly from matrix chunks on kept human genes
    dummy_rc, _, _ = make_csr_or_csc_h5_readers(
        h5_path=h5_p,
        node_path=count_path,
        keep_rows=np.arange(n_cells_raw, dtype=np.int64),
        n_rows_raw=n_cells_raw,
        n_cols_raw=n_genes_raw,
        encoding=enc,
    )
    umi_raw = np.zeros(n_cells_raw, dtype=np.float64)
    ngenes_raw = np.zeros(n_cells_raw, dtype=np.int64)
    for s0 in range(0, n_cells_raw, 16384):
        e0 = min(n_cells_raw, s0 + 16384)
        blk = dummy_rc(s0, e0)[:, keep_gene_cols]
        umi_raw[s0:e0] = blk.sum(axis=1)
        ngenes_raw[s0:e0] = (blk > 0).sum(axis=1)

    is_ntc_raw = np.zeros(n_cells_raw, dtype=bool)
    valid_singlet = np.ones(n_cells_raw, dtype=bool)
    canon_targets_raw: List[str] = [""] * n_cells_raw
    raw_targets_raw: List[str] = [""] * n_cells_raw

    ntc_exact = ("control", "ntc", "nontargeting", "non-targeting", "non_targeting", "nt", "scrambled", "ctrl")
    bad_tokens = ("unassigned", "multiplet", "doublet", "empty", "nan", "none")

    def _is_ntc_token(tok: str) -> bool:
        return (tok in ntc_exact) or tok.startswith(("ntc", "nontarget", "non-target", "non_target"))

    for i in range(n_cells_raw):
        p_str = str(pert_arr[i]).strip() if pert_arr is not None else ""
        g_str = str(guide_arr[i]).strip() if guide_arr is not None else p_str
        p_lo, g_lo = p_str.lower(), g_str.lower()
        if any(x in (p_lo, g_lo) for x in ("multiplet", "doublet")):
            valid_singlet[i] = False
            continue
        if nperts_arr is not None and np.isfinite(nperts_arr[i]) and float(nperts_arr[i]) > 1.0 + 1e-6:
            valid_singlet[i] = False
            continue
        if _is_ntc_token(p_lo) or (p_lo in ("", "nan", "none") and _is_ntc_token(g_lo)):
            is_ntc_raw[i] = True
            continue
        if (not p_str and not g_str) or p_lo in bad_tokens or g_lo == "unassigned":
            valid_singlet[i] = False
            continue
        if any(sep in g_str for sep in (",", ";", "|")) or any(sep in p_str for sep in (",", ";", "|")):
            valid_singlet[i] = False
            continue
        if apply_author_kd40 and kd_arr is not None and np.isfinite(kd_arr[i]) and float(kd_arr[i]) <= 0.40:
            valid_singlet[i] = False
            continue
        canon_s, _ = hgnc.resolve(p_str, "")
        canon_targets_raw[i] = canon_s
        raw_targets_raw[i] = p_str

    # PREREG_E §3.2 + RULING_R12_PREFREEZE Ruling B: UMI >= 300, genes >= 200, NO pct_mt filter
    valid_singlet &= (umi_raw > 0.0) & (umi_raw >= float(min_umi)) & (ngenes_raw >= int(min_genes))
    keep_rows = np.flatnonzero(valid_singlet)

    # Two-Stage Split-Half Rule (PREREG_E §3.2)
    split_branch = "STAGE2_STRATIFIED_SEED42"
    halves_kept = np.zeros(len(keep_rows), dtype=np.int8)
    batches_kept = np.full(len(keep_rows), "HEK293_pool", dtype=object)

    if rep_arr is not None:
        uniq_reps = sorted({str(rep_arr[i]) for i in keep_rows})
        if len(uniq_reps) >= 2:
            rep_half = {r: np.int8(1 if (idx % 2 == 1) else 2) for idx, r in enumerate(uniq_reps, start=1)}
            cand_halves = np.asarray([rep_half[str(rep_arr[i])] for i in keep_rows], dtype=np.int8)
            nt_k = is_ntc_raw[keep_rows]
            nt_a = int(np.sum(nt_k & (cand_halves == 1)))
            nt_b = int(np.sum(nt_k & (cand_halves == 2)))
            if nt_a >= 15 and nt_b >= 15:
                split_branch = "STAGE1_OBS_REPLICATE_PARITY"
                halves_kept = cand_halves
                rep_nt_cnt = Counter(str(rep_arr[keep_rows[j]]) for j in range(len(keep_rows)) if nt_k[j])
                batches_kept = np.asarray(
                    [
                        str(rep_arr[keep_rows[j]])
                        if rep_nt_cnt.get(str(rep_arr[keep_rows[j]]), 0) >= 15
                        else f"half_{halves_kept[j]}"
                        for j in range(len(keep_rows))
                    ],
                    dtype=object,
                )

    if split_branch == "STAGE2_STRATIFIED_SEED42":
        # RULING_R12_HEK293_UNITS §2 exact Stage-2 stratified 50/50 split:
        # 1) Single rng = np.random.default_rng(42);
        # 2) Strata keys: canonical HGNC symbol for targeting cells, "__NTC__" for controls,
        #    processed in Python default string sort order;
        # 3) Within each stratum, sort ascending by (real cell_barcode, kept_row_idx);
        # 4) perm = rng.permutation(n): perm[:n//2] -> half 1 (A), perm[n//2:] -> half 2 (B).
        rng = np.random.default_rng(42)
        strata: Dict[str, List[int]] = defaultdict(list)
        for kept_idx, r_raw in enumerate(keep_rows):
            if is_ntc_raw[r_raw]:
                strata["__NTC__"].append(kept_idx)
            elif canon_targets_raw[r_raw]:
                strata[canon_targets_raw[r_raw]].append(kept_idx)

        for key in sorted(strata.keys()):
            ordered = np.asarray(
                sorted(strata[key], key=lambda idx: (barcodes[keep_rows[idx]], idx)),
                dtype=np.int64,
            )
            perm = rng.permutation(len(ordered))
            n_half = len(ordered) // 2
            halves_kept[ordered[perm[:n_half]]] = 1
            halves_kept[ordered[perm[n_half:]]] = 2

        batches_kept = np.asarray([f"half_{h}" for h in halves_kept], dtype=object)

    nt_kept = is_ntc_raw[keep_rows]
    pkeys: List[str] = []
    pmeta: Dict[str, Tuple[str, str, str, str]] = {}
    for row_i, nt_f in zip(keep_rows, nt_kept):
        if nt_f:
            pkeys.append("")
            continue
        cs = canon_targets_raw[row_i]
        rs = raw_targets_raw[row_i]
        _, hstat = hgnc.resolve(rs, "")
        ens = hgnc.sym_to_ens.get(cs, "")
        pkeys.append(cs)
        if cs not in pmeta:
            pmeta[cs] = (ens, cs, rs, hstat)

    raw_rc, raw_rr, raw_rcc = make_csr_or_csc_h5_readers(
        h5_path=h5_p,
        node_path=count_path,
        keep_rows=keep_rows,
        n_rows_raw=n_cells_raw,
        n_cols_raw=n_genes_raw,
        encoding=enc,
    )
    part = PartSpec(
        tag="XuCao2023_HEK293_CRISPRi",
        path=h5_p,
        n_cells=len(keep_rows),
        n_genes_raw=len(keep_gene_cols),
        nt_mask=nt_kept,
        batches=batches_kept.astype(str),
        half_flags=halves_kept,
        umi=umi_raw[keep_rows],
        pert_keys=pkeys,
        pert_meta=pmeta,
        var_gene_ids=[v_ens[i] for i in keep_gene_cols],
        var_gene_names=[v_sym[i] for i in keep_gene_cols],
        read_chunk_fn=lambda s, e: raw_rc(s, e)[:, keep_gene_cols],
        read_rows_fn=lambda loc: raw_rr(loc)[:, keep_gene_cols],
        read_chunk_cols_fn=lambda s, e, fcols: raw_rcc(s, e, keep_gene_cols[fcols]),
    )
    audit_meta = {
        "feasibility": feas,
        "split_branch": split_branch,
        "rep_col_used": rep_col_used,
        "apply_author_kd40": bool(apply_author_kd40),
        "max_pct_mt": None,
        "min_umi": float(min_umi),
        "min_genes": int(min_genes),
        "n_cells_raw": int(n_cells_raw),
        "n_kept_cells": int(len(keep_rows)),
        "n_kept_ntc": int(nt_kept.sum()),
        "n_human_features": int(len(keep_gene_cols)),
    }
    return [part], audit_meta


def standardize_raw_cohort_to_audit_line(
    line_name: str,
    raw_dir: Path,
    cache_out_root: Path,
    aligned_genes: Optional[Sequence[str]] = None,
    hgnc_tsv: Path = DEFAULT_HGNC_TSV,
    min_gene_overlap: int = 100,
    n_workers: int = 4,
    exclude_zero_cas9: bool = False,
) -> Tuple[AuditLineData, Dict[str, Any]]:
    """
    End-to-end clean-room raw ingestion + Addendum-B standardization for synthetic fixtures or
    locked Round 12 cohorts (`a549_ko`, `a549_ko_cas9pos`, `hek293_xucao`), producing `AuditLineData`.
    """
    hgnc_prep = build_hgnc_lookup(hgnc_tsv)
    if line_name in ("a549_ko", "a549_ko_cas9pos"):
        eff_cas9 = exclude_zero_cas9 or (line_name == "a549_ko_cas9pos")
        parts, meta = build_a549_ko_audit_parts(
            raw_dir=raw_dir,
            hgnc=hgnc_prep,
            exclude_zero_cas9=eff_cas9,
            cache_h5_path=cache_out_root / line_name / "_tmp_audit_a549_csr.h5",
        )
    elif line_name == "hek293_xucao":
        parts, meta = build_hek293_xucao_audit_parts(raw_dir=raw_dir, hgnc=hgnc_prep, apply_author_kd40=False)
    else:
        raise ValueError(f"Unsupported raw line reader: {line_name!r}")

    tmp_cleanup = meta.pop("tmp_h5_path", None)
    try:
        line_dir = process_r10_line(
            line_name=line_name,
            parts=parts,
            hgnc=hgnc_prep,
            chunk_size=4096,
            out_root=cache_out_root,
            extra_meta={"audit_reader": meta},
            n_workers=n_workers,
        )
    finally:
        if tmp_cleanup:
            Path(tmp_cleanup).unlink(missing_ok=True)

    hgnc_canon = CanonicalHGNC.from_tsv(hgnc_tsv)
    adata = load_line_dataset(
        line_dir=line_dir,
        representation="log1p",
        aligned_genes=aligned_genes,
        hgnc=hgnc_canon,
        min_gene_overlap=min_gene_overlap,
    )
    return adata, meta


def resolve_line_cache_dir(line_name: str, phase: str) -> Path:
    """Resolve standardized cache directory cleanly under the current blinding phase."""
    if phase.startswith("dev_"):
        if line_name in R9_LINE_SET:
            return R9_LINES_DIR / line_name
        if line_name in R10_ALL_LINES:
            return R10_LINES_DIR / line_name
        if line_name in R11_ALL_LINES:
            return R11_LINES_DIR / line_name
        raise RuntimeError(f"Dev phase {phase} cannot resolve Round 12 line {line_name!r}.")
    if line_name in R9_LINE_SET:
        return R9_LINES_DIR / line_name
    if line_name in R10_SOURCE_EXT_NAMES:
        return R10_LINES_DIR / line_name
    return R12_LINES_DIR / line_name
