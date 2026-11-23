#!/usr/bin/env python3
"""
Independent Round 13 GSE327644 Raw H5 Reader, Feasibility Checker, & Blinding Gate (A-T3).

Strictly follows locked `PREREG_F.md` (`commit 5b51e5b`) and `PREREG_E.md`:
  - Clean-room isolation: NEVER reads, imports, or references `round13/src/` or
    `round13/results/`.
  - Blinding gates:
      * `--phase synthetic_test`: synthetic fixture tests only; blocks `round13/shared/`.
      * `--phase dev_r9 | dev_r10 | dev_r11 | dev_r12`: historical unblinded caches only;
        blocks `round13/shared/`.
      * `--phase reader_xcheck_r13`: raw GSE327644 ingestion allowed ONLY after BOTH
        `round13/status/PREREG_F_COMMIT.txt` AND `R13_CODE_PATCHED` (or
        `CODE_FREEZE_R13.txt`) exist in `round13/status/`.
      * `--phase holdout_r13`: evaluation on `round13/shared/lines/` allowed ONLY after
        `PREREG_F_COMMIT.txt`, `CODE_FREEZE_R13.txt`, and `CONFIRMATORY_LINES_R13.txt`.
  - GSE327644 (`PREREG_F` §2, §3, §4, §9):
      * Uses Cell Ranger filtered `.h5` (`Gene Expression` features only in X, int64 indptr).
      * Single-sgRNA cells: Cell Ranger protospacer calls (`num_features == 1`) preferred;
        otherwise §3.2 fallback fits Poisson-Gaussian histogram thresholds per guide across
        ALL cells with `>= 200` detected genes within each line (`mac_ctrl` / `mac_lps`),
        retaining single-guide dominant cells (`n_pos == 1`, `top1 >= 5`,
        `top2 < 3 or top1 / max(top2, 0.5) >= 3.0`).
      * NTC: all 593 `'Control'` guides (`NTC-001`..`NTC-593`; indistinguishable).
      * Targets: §3.4 stem extraction (`guide_id.rsplit('-', 1)[0]`, Excel `1-Jun..4-Jun`
        -> `JUN` correction) + frozen `HGNCLookup.resolve(stem, ensembl_id)`.
      * Cell QC: `n_umis >= 300` and `n_genes >= 200`, `max_pct_mt = None`.
      * Split-half: §9.1 channel-number parity over ALL channels present (`odd -> Half A`,
        `even -> Half B`; unit = channel with `< 15` NTC half-pooled fallback via
        `_apply_min_ntc_batch_fallback`); Stage-2 `RULING_R12_HEK293_UNITS §2` `seed=42`
        fallback strictly only when channel metadata is absent.
"""

from __future__ import annotations

import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Set, Tuple

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROUND13_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = ROUND13_ROOT.parent
R10_SCRIPTS_DIR = REPO_ROOT / "round10" / "shared" / "scripts"
R11_AUDIT_DIR = REPO_ROOT / "round11" / "audit"
R11_SRC_DIR = REPO_ROOT / "round11" / "src"
R12_AUDIT_DIR = REPO_ROOT / "round12" / "audit"

for _p in (str(R10_SCRIPTS_DIR), str(R11_AUDIT_DIR), str(R11_SRC_DIR), str(R12_AUDIT_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from dataset_readers_r11 import _apply_min_ntc_batch_fallback  # noqa: E402
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
from independent_r12_loader import R12_LINES_DIR, R12_TARGET_KEYS  # noqa: E402
from prepare_r10_lines import (  # noqa: E402
    HGNCLookup,
    PartSpec,
    build_hgnc_lookup,
    extract_top2_csc,
    fit_poisson_gaussian_thresholds_from_hist,
    process_r10_line,
)

R13_SHARED_DIR = ROUND13_ROOT / "shared"
R13_RAW_DIR = R13_SHARED_DIR / "raw"
R13_LINES_DIR = R13_SHARED_DIR / "lines"
R13_STATUS_DIR = ROUND13_ROOT / "status"

R13_TARGET_KEYS: Tuple[str, ...] = ("mac_ctrl", "mac_lps")
R13_TARGET_SET: Set[str] = set(R13_TARGET_KEYS)

R10_ALL_LINES: Tuple[str, ...] = ("dld1", "cd4t_rest", "wtc11_i3neuron", "calu3")
R11_ALL_LINES: Tuple[str, ...] = ("telohaec", "cd4t_stim48", "feng_ipsc")
R12_ALL_LINES: Tuple[str, ...] = R12_TARGET_KEYS
DEV_ALLOWED_LINES_BY_PHASE: Dict[str, Tuple[str, ...]] = {
    "dev_r9": R9_LINE_NAMES,
    "dev_r10": R10_ALL_LINES,
    "dev_r11": R11_ALL_LINES,
    "dev_r12": R12_ALL_LINES,
}
ALL_DEV_LINES: Set[str] = set(R9_LINE_NAMES) | set(R10_ALL_LINES) | set(R11_ALL_LINES) | set(R12_ALL_LINES)

VALID_R13_ROLES: Set[str] = {"CONFIRMATORY", "SECONDARY_ONLY", "EXCLUDED_FEASIBILITY"}

_CHANNEL_NUM_RE = re.compile(r"(?:^|[_-])(?:Mac_)?L0*(\d+)(?:[_-]|\.h5$)", re.IGNORECASE)
_JUN_EXCEL_FIX: Dict[str, str] = {
    "JUN-1": "1-Jun",
    "JUN-2": "2-Jun",
    "JUN-3": "3-Jun",
    "JUN-4": "4-Jun",
}
_SAFE_CONTROL_RE = re.compile(r"safe[-_ ]?target|aavs1|intergenic", re.IGNORECASE)


def _decode_str(arr: Iterable[Any]) -> List[str]:
    out: List[str] = []
    for x in arr:
        if isinstance(x, (bytes, np.bytes_)):
            out.append(x.decode("utf-8", errors="replace"))
        else:
            out.append(str(x))
    return out


def _has_r13_code_patched(status_dir: Path) -> str:
    """Return marker evidence string if `R13_LINES_AUTHORIZED` / `R13_CODE_PATCHED` / `CODE_FREEZE_R13` is set."""
    for fname in ("R13_LINES_AUTHORIZED.txt", "R13_CODE_PATCHED", "R13_CODE_PATCHED.txt", "CODE_FREEZE_R13.txt"):
        tok = _first_token(status_dir / fname)
        if tok:
            return f"{fname}:{tok}"
    main_p = status_dir / "main.txt"
    if main_p.exists():
        for ln in main_p.read_text(encoding="utf-8", errors="replace").splitlines():
            s = ln.strip()
            if s == "R13_CODE_PATCHED" or s.startswith("R13_CODE_PATCHED "):
                return f"main.txt:{s}"
    return ""


def parse_confirmatory_lines_r13(conf_file: Path = R13_STATUS_DIR / "CONFIRMATORY_LINES_R13.txt") -> Dict[str, Any]:
    """
    Parse `CONFIRMATORY_LINES_R13.txt` supporting line=ROLE or bare confirmatory line names
    (`PREREG_F §2, §5`: at most `mac_ctrl` may be `CONFIRMATORY`, `K_conf in {0, 1}`).
    """
    if not conf_file.exists():
        raise FileNotFoundError(f"Missing Round 13 confirmatory lines file: {conf_file}")
    roles: Dict[str, str] = {}
    confirmatory: List[str] = []
    for raw_ln in conf_file.read_text(encoding="utf-8").splitlines():
        ln = raw_ln.strip()
        if not ln or ln.startswith("#"):
            continue
        if "=" in ln:
            k, v = [p.strip() for p in ln.split("=", 1)]
            if k in R13_TARGET_SET:
                role = v.split()[0].strip()
                if role not in VALID_R13_ROLES:
                    raise ValueError(f"Unexpected role {role!r} for {k!r} in {conf_file}")
                roles[k] = role
                if role == "CONFIRMATORY":
                    confirmatory.append(k)
        elif ln in R13_TARGET_SET:
            roles[ln] = "CONFIRMATORY"
            confirmatory.append(ln)

    if "mac_lps" in confirmatory:
        raise ValueError(f"PREREG_F §2 violation: mac_lps cannot be CONFIRMATORY in {conf_file}")
    if len(confirmatory) > 1:
        raise ValueError(f"PREREG_F §6 violation: max K_conf = 1, found {confirmatory}")
    return {"roles": roles, "confirmatory_lines": confirmatory, "k_conf": len(confirmatory)}


def enforce_r13_gate(
    phase: str,
    target_lines: Sequence[str] = (),
    status_dir: Path = R13_STATUS_DIR,
    extra_dirs: Sequence[Path] = (),
    require_confirmatory_file: bool = True,
) -> Dict[str, str]:
    """Enforce Round 13 blinding/freeze locks prior to reading expression arrays."""
    if phase == "synthetic_test" or phase in DEV_ALLOWED_LINES_BY_PHASE:
        r13_shared = R13_SHARED_DIR.resolve()
        for d in extra_dirs:
            rp = Path(d).resolve()
            if rp == r13_shared or r13_shared in rp.parents:
                raise RuntimeError(
                    f"PREREG_F blinding lock ({phase}): refused accessing round13/shared path {d}."
                )
        if phase in DEV_ALLOWED_LINES_BY_PHASE:
            bad = [t for t in target_lines if t in R13_TARGET_SET or t not in ALL_DEV_LINES]
            if bad:
                raise RuntimeError(f"PREREG_F dev lock ({phase}): invalid target lines {bad}.")
        else:
            bad = [t for t in target_lines if t in R13_TARGET_SET]
            if bad:
                # Allowed only if no extra_dirs touched r13_shared (caller uses tmp dirs in tests)
                pass
        return {
            "phase": phase,
            "prereg_commit": _first_token(status_dir / "PREREG_F_COMMIT.txt"),
            "code_patched": "",
            "code_freeze": "",
        }

    prereg_hash = _first_token(status_dir / "PREREG_F_COMMIT.txt")
    if not prereg_hash:
        raise RuntimeError(
            f"PREREG_F hard gate blocked: missing or empty {status_dir / 'PREREG_F_COMMIT.txt'}."
        )

    if phase == "reader_xcheck_r13":
        patched = _has_r13_code_patched(status_dir)
        if not patched:
            raise RuntimeError(
                f"PREREG_F xcheck gate blocked: R13_CODE_PATCHED / CODE_FREEZE_R13 not yet present in {status_dir}."
            )
        return {
            "phase": phase,
            "prereg_commit": prereg_hash,
            "code_patched": patched,
            "code_freeze": _first_token(status_dir / "CODE_FREEZE_R13.txt"),
        }

    if phase == "holdout_r13":
        freeze_hash = _first_token(status_dir / "CODE_FREEZE_R13.txt")
        if not freeze_hash:
            raise RuntimeError(
                f"CODE_FREEZE_R13 hard gate blocked: missing or empty {status_dir / 'CODE_FREEZE_R13.txt'}."
            )
        conf_p = status_dir / "CONFIRMATORY_LINES_R13.txt"
        if require_confirmatory_file and (not conf_p.exists() or not conf_p.read_text(encoding="utf-8").strip()):
            raise RuntimeError(f"CONFIRMATORY_LINES_R13 gate blocked: missing or empty {conf_p}.")
        return {
            "phase": phase,
            "prereg_commit": prereg_hash,
            "code_patched": _has_r13_code_patched(status_dir),
            "code_freeze": freeze_hash,
        }

    raise ValueError(f"Unsupported phase: {phase!r}")


def load_gse327644_guide_template(
    template_path: Path,
    hgnc: HGNCLookup,
) -> Dict[str, Any]:
    """
    Parse `GSE327644_PRJNA1450051_marson_mac_perturbseq_guide_template.txt.gz` (or `.txt`)
    strictly per `PREREG_F §3.3 & §3.4`:
      - Control rule: if `Control` rows distinguish non-targeting vs safe-targeting, ONLY
        non-targeting counts as NTC; otherwise ALL `Control` guides count as NTC.
      - Target rule (§3.4): `stem = guide_id.rsplit('-', 1)[0]` with Excel `1-Jun..4-Jun -> JUN`
        correction, resolved via `(stem, guide_target_gene_id)` in frozen `hgnc`.
    """
    if not template_path.exists():
        raise FileNotFoundError(f"Missing GSE327644 guide template: {template_path}")

    open_fn = gzip.open if template_path.name.endswith(".gz") else open
    with open_fn(template_path, "rt", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    required_cols = {
        "guide_id",
        "guide_protospacer",
        "guide_role",
        "guide_target_gene_id",
        "guide_target_gene_name",
    }
    if not rows or not required_cols.issubset(set(rows[0].keys())):
        raise ValueError(f"Invalid guide_template columns in {template_path}")

    ctrl_rows = [r for r in rows if r["guide_role"].strip() == "Control"]
    tgt_rows = [r for r in rows if r["guide_role"].strip() == "Targeting a Gene"]
    other_roles = sorted({r["guide_role"].strip() for r in rows} - {"Control", "Targeting a Gene"})
    if other_roles:
        raise ValueError(f"Unexpected guide_role values in {template_path}: {other_roles}")

    distinguishes_safe = any(
        _SAFE_CONTROL_RE.search(
            " ".join(
                [
                    r["guide_id"],
                    r["guide_role"],
                    r["guide_target_gene_id"],
                    r["guide_target_gene_name"],
                ]
            )
        )
        for r in ctrl_rows
    )

    ntc_guide_ids: Set[str] = set()
    excluded_safe_ids: Set[str] = set()
    for r in ctrl_rows:
        gid = r["guide_id"].strip()
        text = " ".join([gid, r["guide_target_gene_id"], r["guide_target_gene_name"]])
        if distinguishes_safe and _SAFE_CONTROL_RE.search(text):
            excluded_safe_ids.add(gid)
        else:
            ntc_guide_ids.add(gid)

    guide_to_target: Dict[str, Tuple[str, str, str, str]] = {}
    alias_to_guide_id: Dict[str, str] = {}
    for g in ntc_guide_ids:
        alias_to_guide_id[g] = g

    for r in tgt_rows:
        gid = r["guide_id"].strip()
        raw_col_name = r["guide_target_gene_name"].strip()
        raw_ens = r["guide_target_gene_id"].strip().split(".")[0]
        stem = gid.rsplit("-", 1)[0] if "-" in gid else gid
        if stem in {"1", "2", "3", "4"} and raw_col_name.upper() == "JUN":
            stem = "JUN"
            alias_to_guide_id[f"JUN-{gid.split('-')[0]}"] = gid
        canon_sym, hstat = hgnc.resolve(stem, raw_ens)
        ens_id = raw_ens if raw_ens.startswith("ENSG") else hgnc.sym_to_ens.get(canon_sym, "")
        guide_to_target[gid] = (ens_id, canon_sym, stem, hstat)
        alias_to_guide_id[gid] = gid

    for gid in excluded_safe_ids:
        alias_to_guide_id[gid] = gid

    ordered_guide_ids = [r["guide_id"].strip() for r in rows]
    return {
        "ordered_guide_ids": ordered_guide_ids,
        "ntc_guide_ids": ntc_guide_ids,
        "excluded_safe_ids": excluded_safe_ids,
        "guide_to_target": guide_to_target,
        "alias_to_guide_id": alias_to_guide_id,
        "distinguishes_safe_vs_ntc": bool(distinguishes_safe),
        "n_rows": len(rows),
        "n_control_rows": len(ctrl_rows),
        "n_targeting_rows": len(tgt_rows),
        "n_unique_canonical_targets": len({v[1] for v in guide_to_target.values()}),
    }


def _find_guide_template_path(raw_dir: Path) -> Path:
    candidates: List[Path] = []
    for base in (raw_dir / "gse327644", raw_dir, raw_dir.parent):
        if not base.exists():
            continue
        for pat in (
            "*guide_template*.txt.gz",
            "*guide_template*.txt",
        ):
            candidates.extend(sorted(base.glob(pat)))
    candidates = [p for p in candidates if not p.name.startswith("._")]
    if not candidates:
        raise FileNotFoundError(f"Could not locate GSE327644 guide_template file under {raw_dir}")
    return candidates[0]


def discover_gse327644_line_files(raw_dir: Path, line_name: str) -> List[Tuple[Optional[int], str, Path]]:
    """
    Discover Cell Ranger filtered `.h5` files for `mac_ctrl` (`Mac_Ctrl`) or `mac_lps` (`Mac_LPS`)
    sorted ascending by channel integer (`PREREG_F §2, §9.1`). Rejects `raw_feature_bc_matrix`.
    """
    if line_name == "mac_ctrl":
        cond_re = re.compile(r"(?:^|[_-])Mac_Ctrl(?:[_-]|\.h5$)", re.IGNORECASE)
        other_re = re.compile(r"(?:^|[_-])Mac_LPS(?:[_-]|\.h5$)", re.IGNORECASE)
    elif line_name == "mac_lps":
        cond_re = re.compile(r"(?:^|[_-])Mac_LPS(?:[_-]|\.h5$)", re.IGNORECASE)
        other_re = re.compile(r"(?:^|[_-])Mac_Ctrl(?:[_-]|\.h5$)", re.IGNORECASE)
    else:
        raise ValueError(f"Unsupported GSE327644 line_name: {line_name!r}")

    search_roots = [raw_dir / "gse327644"] if (raw_dir / "gse327644").is_dir() else [raw_dir]
    seen_resolved: Set[Path] = set()
    found: List[Tuple[Optional[int], str, Path]] = []

    for root in search_roots:
        if not root.exists():
            continue
        for hp in sorted(root.rglob("*.h5")):
            if hp.name.startswith(("._", "_tmp", "_audit")):
                continue
            rel_str = str(hp.relative_to(root))
            if "raw_feature_bc_matrix" in rel_str.lower():
                continue
            target_str = hp.name if cond_re.search(hp.name) else rel_str
            if not cond_re.search(target_str) or other_re.search(target_str):
                continue
            rp = hp.resolve()
            if rp in seen_resolved:
                continue
            seen_resolved.add(rp)
            m_ch = _CHANNEL_NUM_RE.search(hp.name) or _CHANNEL_NUM_RE.search(rel_str)
            ch_num = int(m_ch.group(1)) if m_ch else None
            ch_tag = f"L{ch_num:02d}" if ch_num is not None else hp.stem
            found.append((ch_num, ch_tag, hp))

    if not found:
        raise FileNotFoundError(f"No Cell Ranger filtered .h5 files found for {line_name} under {raw_dir}")

    found.sort(key=lambda it: (it[0] is None, it[0] if it[0] is not None else 0, it[1], str(it[2])))
    return found


def _make_int64_10x_csc_subset_readers(
    h5_path: Path,
    keep_cols: np.ndarray,
    n_raw_features: int,
    rna_row_idx: Optional[np.ndarray],
) -> Tuple[
    Callable[[int, int], np.ndarray],
    Callable[[np.ndarray], np.ndarray],
    Callable[[int, int, np.ndarray], np.ndarray],
]:
    """
    Int64-indptr memory-bounded `(read_chunk_fn, read_rows_fn, read_chunk_cols_fn)` over sorted
    `keep_cols` of a 10x CSC `.h5` file `(n_raw_features, n_raw_cells)`.
    """
    n_out_genes = len(rna_row_idx) if rna_row_idx is not None else n_raw_features
    if rna_row_idx is not None and len(rna_row_idx) == n_raw_features and np.array_equal(
        rna_row_idx, np.arange(n_raw_features)
    ):
        rna_row_idx = None

    def _read_sorted_cols(raw_cols: np.ndarray, feat_cols: Optional[np.ndarray] = None) -> np.ndarray:
        out_dtype = np.float32 if feat_cols is not None else np.float64
        out_w = len(feat_cols) if feat_cols is not None else n_out_genes
        if len(raw_cols) == 0:
            return np.zeros((0, out_w), dtype=out_dtype)
        c_min, c_max = int(raw_cols[0]), int(raw_cols[-1])
        span = c_max - c_min + 1
        with h5py.File(h5_path, "r") as hf:
            m = hf["matrix"]
            indptr = m["indptr"]
            data = m["data"]
            indices = m["indices"]
            if len(raw_cols) / max(span, 1) >= 0.15:
                ptrs = np.asarray(indptr[c_min : c_max + 2], dtype=np.int64)
                p0, p1 = int(ptrs[0]), int(ptrs[-1])
                csc_span = sp.csc_matrix(
                    (
                        np.asarray(data[p0:p1]),
                        np.asarray(indices[p0:p1]),
                        ptrs - p0,
                    ),
                    shape=(n_raw_features, span),
                )
                sub = csc_span[:, raw_cols - c_min]
            else:
                col_starts = np.asarray(indptr[raw_cols], dtype=np.int64)
                col_ends = np.asarray(indptr[raw_cols + 1], dtype=np.int64)
                all_d, all_i = [], []
                new_ptr = np.zeros(len(raw_cols) + 1, dtype=np.int64)
                cur_p = 0
                for j, (s0, s1) in enumerate(zip(col_starts, col_ends)):
                    s0_i, s1_i = int(s0), int(s1)
                    if s1_i > s0_i:
                        all_d.append(data[s0_i:s1_i])
                        all_i.append(indices[s0_i:s1_i])
                        cur_p += s1_i - s0_i
                    new_ptr[j + 1] = cur_p
                d_cat = np.concatenate(all_d) if all_d else np.empty(0, dtype=np.float32)
                i_cat = np.concatenate(all_i) if all_i else np.empty(0, dtype=np.int32)
                sub = sp.csc_matrix((d_cat, i_cat, new_ptr), shape=(n_raw_features, len(raw_cols)))

            if feat_cols is not None:
                target_rows = rna_row_idx[feat_cols] if rna_row_idx is not None else feat_cols
                sub = sub[target_rows, :]
            elif rna_row_idx is not None:
                sub = sub[rna_row_idx, :]
            return sub.T.toarray().astype(out_dtype, copy=False)

    def _read_loc_rows(loc_idx: np.ndarray) -> np.ndarray:
        loc_arr = np.asarray(loc_idx, dtype=np.int64)
        raw_cols = keep_cols[loc_arr]
        if len(raw_cols) <= 1 or np.all(raw_cols[:-1] <= raw_cols[1:]):
            return _read_sorted_cols(raw_cols, None)
        order = np.argsort(raw_cols, kind="mergesort")
        inv = np.empty_like(order)
        inv[order] = np.arange(len(order), dtype=np.int64)
        return _read_sorted_cols(raw_cols[order], None)[inv]

    return (
        (lambda s, e: _read_sorted_cols(keep_cols[s:e], None)),
        _read_loc_rows,
        (lambda s, e, fcols: _read_sorted_cols(keep_cols[s:e], fcols)),
    )


def _match_protospacer_csv_for_h5(h5_path: Path, ch_tag: str, cond_token: str) -> Optional[Path]:
    """Locate unambiguous companion Cell Ranger `protospacer_calls_per_cell.csv` if shipped."""
    pref = re.sub(r"(?:_sample)?_filtered_feature_bc_matrix\.h5$", "", h5_path.name, flags=re.IGNORECASE)
    if pref == h5_path.name:
        pref = h5_path.stem
    cand_dirs = [h5_path.parent]
    stem_dir = h5_path.parent / pref
    if stem_dir.is_dir():
        cand_dirs.append(stem_dir)
    gsm_m = re.match(r"^(GSM\d+)", h5_path.name)
    gsm = gsm_m.group(1) if gsm_m else None

    hits: List[Path] = []
    ch_only_hits: List[Path] = []
    seen_p: Set[Path] = set()
    for d in cand_dirs:
        for pat in ("*protospacer_calls_per_cell.csv*", "*protospacer_calls*.csv*"):
            for p in sorted(d.glob(pat)):
                if p.name.startswith("._"):
                    continue
                rp = p.resolve()
                if rp in seen_p:
                    continue
                seen_p.add(rp)
                if p.parent == h5_path.parent:
                    low = p.name.lower()
                    if low.startswith(pref.lower()):
                        hits.append(p)
                        continue
                    ok_ch = ch_tag.lower() in low or (gsm is not None and gsm.lower() in low)
                    ok_cond = cond_token.lower() in low
                    if ok_ch and not ok_cond and "mac_ctrl" not in low and "mac_lps" not in low:
                        ch_only_hits.append(p)
                    if not (ok_ch and ok_cond):
                        continue
                hits.append(p)
    if not hits and len(ch_only_hits) == 1:
        return ch_only_hits[0]
    if len(hits) == 1:
        return hits[0]
    if len(hits) > 1:
        exact = [p for p in hits if p.name.lower().startswith(pref.lower())]
        if len(exact) == 1:
            return exact[0]
        return None
    return None


def _read_embedded_h5_protospacer_calls(hf: h5py.File, barcodes: Sequence[str], alias_map: Dict[str, str]) -> Optional[Dict[str, str]]:
    """Check if an `.h5` embeds explicit Cell Ranger per-cell protospacer calls matching barcodes."""
    for grp_name in ("protospacer_calls", "matrix/protospacer_calls", "obs"):
        if grp_name not in hf:
            continue
        grp = hf[grp_name]
        if not isinstance(grp, h5py.Group):
            continue
        keys = set(grp.keys())
        if "num_features" in keys and "feature_call" in keys:
            nf = np.asarray(grp["num_features"][:])
            fc = _decode_str(grp["feature_call"][:])
            bc_key = "cell_barcode" if "cell_barcode" in keys else ("barcode" if "barcode" in keys else None)
            if bc_key is not None:
                bcs = _decode_str(grp[bc_key][:])
            elif len(fc) == len(barcodes):
                bcs = list(barcodes)
            else:
                continue
            out: Dict[str, str] = {}
            for b, n_val, call_s in zip(bcs, nf, fc):
                try:
                    n_int = int(n_val)
                except (TypeError, ValueError):
                    continue
                if n_int == 1:
                    c = call_s.strip()
                    canon_g = alias_map.get(c)
                    if canon_g is not None:
                        out[b] = canon_g
            return out
    return None


def _read_csv_protospacer_calls(csv_path: Path, alias_map: Dict[str, str]) -> Dict[str, str]:
    df = pd.read_csv(csv_path, dtype=str, keep_default_na=False)
    bc_col = next((c for c in ("cell_barcode", "barcode") if c in df.columns), None)
    if bc_col is None or "num_features" not in df.columns or "feature_call" not in df.columns:
        raise ValueError(f"Invalid protospacer CSV columns in {csv_path}: {list(df.columns)}")
    out: Dict[str, str] = {}
    for bc, nf_s, fc_s in zip(df[bc_col], df["num_features"], df["feature_call"]):
        try:
            nf = int(float(str(nf_s).strip()))
        except ValueError:
            continue
        if nf != 1:
            continue
        c = str(fc_s).strip()
        if any(sep in c for sep in ("|", ";", ",")):
            continue
        canon_g = alias_map.get(c)
        if canon_g is not None:
            out[str(bc).strip()] = canon_g
    return out


def check_gse327644_feasibility(
    line_name: str,
    raw_dir: Path,
    hgnc: HGNCLookup,
    guide_template_path: Optional[Path] = None,
) -> Dict[str, Any]:
    """
    Pre-metric technical feasibility audit for `GSE327644` (`PREREG_F §5.1`):
      1. Non-negative finite integer raw UMI counts in `Gene Expression` features;
      2. Valid per-cell guide source (Cell Ranger protospacer calls OR non-`Gene Expression`
         CRISPR Guide Capture UMI rows matching `guide_template.txt.gz`).
    """
    try:
        tpl_path = guide_template_path or _find_guide_template_path(raw_dir)
        tpl = load_gse327644_guide_template(tpl_path, hgnc)
        files = discover_gse327644_line_files(raw_dir, line_name)
    except Exception as exc:
        return {
            "line": line_name,
            "feasible": False,
            "status": "EXCLUDED_FEASIBILITY",
            "reason": f"discovery_error: {exc}",
        }

    cond_token = "Mac_Ctrl" if line_name == "mac_ctrl" else "Mac_LPS"
    alias_map: Dict[str, str] = tpl["alias_to_guide_id"]
    ch_records: List[Dict[str, Any]] = []
    n_proto_channels = 0

    for ch_num, ch_tag, hp in files:
        try:
            with h5py.File(hp, "r") as hf:
                if "matrix" not in hf:
                    return {
                        "line": line_name,
                        "feasible": False,
                        "status": "EXCLUDED_FEASIBILITY",
                        "reason": f"missing /matrix group in {hp.name}",
                    }
                m = hf["matrix"]
                for req in ("shape", "data", "indices", "indptr", "barcodes", "features/feature_type", "features/id", "features/name"):
                    if req not in m:
                        return {
                            "line": line_name,
                            "feasible": False,
                            "status": "EXCLUDED_FEASIBILITY",
                            "reason": f"missing /matrix/{req} in {hp.name}",
                        }
                shape = tuple(int(v) for v in m["shape"][:])
                n_feat, n_col = shape
                ftypes = _decode_str(m["features/feature_type"][:])
                fids = _decode_str(m["features/id"][:])
                fnames = _decode_str(m["features/name"][:])
                barcodes = _decode_str(m["barcodes"][:])
                rna_mask = np.asarray([ft == "Gene Expression" for ft in ftypes], dtype=bool)
                rna_rows = np.flatnonzero(rna_mask)
                crp_rows = np.flatnonzero(~rna_mask)
                if len(rna_rows) == 0 or n_col == 0:
                    return {
                        "line": line_name,
                        "feasible": False,
                        "status": "EXCLUDED_FEASIBILITY",
                        "reason": f"empty Gene Expression features or cells in {hp.name}",
                    }

                # Integer & non-negative check across all Gene Expression entries (chunked)
                data_ds = m["data"]
                idx_ds = m["indices"]
                is_int_dtype = np.issubdtype(data_ds.dtype, np.integer)
                nnz_total = int(data_ds.shape[0])
                chunk = 2_000_000
                for s in range(0, nnz_total, chunk):
                    e = min(s + chunk, nnz_total)
                    d_arr = np.asarray(data_ds[s:e])
                    i_arr = np.asarray(idx_ds[s:e])
                    r_vals = d_arr[rna_mask[i_arr]]
                    if r_vals.size == 0:
                        continue
                    if is_int_dtype:
                        if np.any(r_vals < 0):
                            return {
                                "line": line_name,
                                "feasible": False,
                                "status": "EXCLUDED_FEASIBILITY",
                                "reason": f"negative integer count in Gene Expression of {hp.name}",
                            }
                    else:
                        rf = r_vals.astype(np.float64, copy=False)
                        if not np.all(np.isfinite(rf) & (rf >= 0.0) & (rf == np.floor(rf))):
                            return {
                                "line": line_name,
                                "feasible": False,
                                "status": "EXCLUDED_FEASIBILITY",
                                "reason": f"non-integer or non-finite Gene Expression count in {hp.name}",
                            }

                emb_calls = _read_embedded_h5_protospacer_calls(hf, barcodes, alias_map)
                csv_p = _match_protospacer_csv_for_h5(hp, ch_tag, cond_token)
                has_proto = (emb_calls is not None) or (csv_p is not None)
                if has_proto:
                    n_proto_channels += 1

                matched_crp = 0
                for ri in crp_rows:
                    if fids[ri].strip() in alias_map or fnames[ri].strip() in alias_map:
                        matched_crp += 1
                if not has_proto and matched_crp == 0:
                    return {
                        "line": line_name,
                        "feasible": False,
                        "status": "EXCLUDED_FEASIBILITY",
                        "reason": f"no protospacer calls and 0 matched CRISPR guide rows in {hp.name}",
                    }
                ch_records.append(
                    {
                        "channel": ch_tag,
                        "ch_num": ch_num,
                        "file": str(hp),
                        "n_cells_raw": int(n_col),
                        "n_rna_features": int(len(rna_rows)),
                        "n_crp_features": int(len(crp_rows)),
                        "n_crp_matched": int(matched_crp),
                        "has_proto": bool(has_proto),
                    }
                )
        except Exception as exc:
            return {
                "line": line_name,
                "feasible": False,
                "status": "EXCLUDED_FEASIBILITY",
                "reason": f"hdf5_read_error in {hp.name}: {exc}",
            }

    use_proto = n_proto_channels == len(files)
    return {
        "line": line_name,
        "feasible": True,
        "status": "FEASIBLE",
        "guide_call_mode": "CELLRANGER_PROTOSPACER_CALLS" if use_proto else "H5_CRISPR_POISSON_GAUSSIAN_LINE_FIT",
        "n_channels": len(files),
        "n_proto_channels": int(n_proto_channels),
        "channels": ch_records,
        "guide_template_n_ntc": len(tpl["ntc_guide_ids"]),
        "guide_template_n_targets": len(tpl["guide_to_target"]),
    }


def _apply_stage2_seed42_split(
    nt_mask: np.ndarray,
    canon_targets: Sequence[str],
    barcodes: Sequence[str],
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Exact `RULING_R12_HEK293_UNITS.md §2` / `PREREG_F §4.2` Stage-2 fallback:
      1. One generator `rng = np.random.default_rng(42)`;
      2. Strata key = canonical HGNC symbol for targeting cells, `"__NTC__"` for controls,
         processed in Python ascending string order;
      3. Within stratum sorted ascending by `(true_barcode, kept_row_idx)`;
      4. `perm = rng.permutation(n)`: `perm[:n//2]` -> half 1, `perm[n//2:]` -> half 2;
      5. Batch centering unit = `["half_1", "half_2"]`.
    """
    n_kept = len(nt_mask)
    rng = np.random.default_rng(42)
    strata: Dict[str, List[int]] = defaultdict(list)
    for k_idx in range(n_kept):
        if nt_mask[k_idx]:
            strata["__NTC__"].append(k_idx)
        else:
            strata[canon_targets[k_idx]].append(k_idx)

    halves = np.zeros(n_kept, dtype=np.int8)
    for skey in sorted(strata.keys()):
        ordered = np.asarray(
            sorted(strata[skey], key=lambda idx: (str(barcodes[idx]), int(idx))),
            dtype=np.int64,
        )
        perm = rng.permutation(len(ordered))
        n_half = len(ordered) // 2
        halves[ordered[perm[:n_half]]] = 1
        halves[ordered[perm[n_half:]]] = 2
    batches = np.asarray([f"half_{int(h)}" for h in halves], dtype=str)
    return halves, batches


def build_gse327644_audit_parts(
    line_name: str,
    raw_dir: Path,
    hgnc: HGNCLookup,
    guide_template_path: Optional[Path] = None,
    min_umi: float = 300.0,
    min_genes: int = 200,
    min_batch_nt: int = 15,
) -> Tuple[List[PartSpec], Dict[str, Any]]:
    """
    Clean-room reader building `PartSpec` list for `mac_ctrl` or `mac_lps` (`GSE327644`)
    strictly following locked `PREREG_F §2, §3, §4, §5, §9`.
    """
    feas = check_gse327644_feasibility(
        line_name=line_name,
        raw_dir=raw_dir,
        hgnc=hgnc,
        guide_template_path=guide_template_path,
    )
    if not feas["feasible"]:
        raise RuntimeError(f"GSE327644 {line_name} failed pre-metric feasibility: {feas['reason']}")

    tpl_path = guide_template_path or _find_guide_template_path(raw_dir)
    tpl = load_gse327644_guide_template(tpl_path, hgnc)
    ntc_guide_ids: Set[str] = tpl["ntc_guide_ids"]
    guide_to_target: Dict[str, Tuple[str, str, str, str]] = tpl["guide_to_target"]
    alias_map: Dict[str, str] = tpl["alias_to_guide_id"]
    all_template_guides: List[str] = tpl["ordered_guide_ids"]
    tpl_g_to_idx: Dict[str, int] = {g: i for i, g in enumerate(all_template_guides)}

    files = discover_gse327644_line_files(raw_dir, line_name)
    cond_token = "Mac_Ctrl" if line_name == "mac_ctrl" else "Mac_LPS"
    guide_call_mode = feas["guide_call_mode"]

    # Determine Stage 1 channel parity vs Stage 2 seed-42 split (`PREREG_F §4, §9.1`)
    parsed_ch_nums = [c for c, _, _ in files if c is not None]
    has_odd = any(c % 2 == 1 for c in parsed_ch_nums)
    has_even = any(c % 2 == 0 for c in parsed_ch_nums)
    use_stage1 = (len(parsed_ch_nums) == len(files)) and (len(set(parsed_ch_nums)) >= 2) and has_odd and has_even
    split_branch = "STAGE1_CHANNEL_PARITY" if use_stage1 else "STAGE2_STRATIFIED_SEED42"

    # Pass 1 across all channels present: RNA QC + protospacer calls OR pooled guide histogram
    hist = np.zeros((len(all_template_guides), 66), dtype=np.int64)
    channel_cache: List[Dict[str, Any]] = []

    for ch_num, ch_tag, hp in files:
        with h5py.File(hp, "r") as hf:
            m = hf["matrix"]
            n_feat, n_col = int(m["shape"][0]), int(m["shape"][1])
            ftypes = _decode_str(m["features/feature_type"][:])
            fids = _decode_str(m["features/id"][:])
            fnames = _decode_str(m["features/name"][:])
            barcodes = _decode_str(m["barcodes"][:])

            rna_mask = np.asarray([ft == "Gene Expression" for ft in ftypes], dtype=bool)
            rna_rows = np.flatnonzero(rna_mask)
            crp_rows = np.flatnonzero(~rna_mask)

            crp_to_tpl = np.full(n_feat, -1, dtype=np.int32)
            for ri in crp_rows:
                g_match = alias_map.get(fids[ri].strip()) or alias_map.get(fnames[ri].strip())
                if g_match is not None and g_match in tpl_g_to_idx:
                    crp_to_tpl[ri] = tpl_g_to_idx[g_match]
            is_valid_crp_feat = crp_to_tpl >= 0

            emb_calls = _read_embedded_h5_protospacer_calls(hf, barcodes, alias_map) if guide_call_mode == "CELLRANGER_PROTOSPACER_CALLS" else None
            csv_p = _match_protospacer_csv_for_h5(hp, ch_tag, cond_token) if (guide_call_mode == "CELLRANGER_PROTOSPACER_CALLS" and emb_calls is None) else None
            proto_map = emb_calls if emb_calls is not None else (_read_csv_protospacer_calls(csv_p, alias_map) if csv_p is not None else None)

            umi_l = np.zeros(n_col, dtype=np.float64)
            ng_l = np.zeros(n_col, dtype=np.int64)
            crp_idx_parts: List[np.ndarray] = []
            crp_dat_parts: List[np.ndarray] = []
            crp_col_parts: List[np.ndarray] = []

            indptr_all = np.asarray(m["indptr"][:], dtype=np.int64)
            data_ds = m["data"]
            indices_ds = m["indices"]

            for cs in range(0, n_col, 25_000):
                ce = min(cs + 25_000, n_col)
                p0, p1 = int(indptr_all[cs]), int(indptr_all[ce])
                if p1 == p0:
                    continue
                dat = np.asarray(data_ds[p0:p1])
                idx = np.asarray(indices_ds[p0:p1], dtype=np.int64)
                col_lens = np.diff(indptr_all[cs : ce + 1])
                col_ids = np.repeat(np.arange(cs, ce, dtype=np.int64), col_lens)

                pos_mask = dat > 0
                if not np.any(pos_mask):
                    continue
                dat = dat[pos_mask]
                idx = idx[pos_mask]
                col_ids = col_ids[pos_mask]

                m_r = rna_mask[idx]
                if np.any(m_r):
                    r_dat = dat[m_r].astype(np.float64, copy=False)
                    r_col = col_ids[m_r] - cs
                    umi_l[cs:ce] += np.bincount(r_col, weights=r_dat, minlength=ce - cs)
                    ng_l[cs:ce] += np.bincount(r_col, minlength=ce - cs)

                if proto_map is None:
                    m_c = is_valid_crp_feat[idx]
                    if np.any(m_c):
                        crp_idx_parts.append(crp_to_tpl[idx[m_c]])
                        crp_dat_parts.append(dat[m_c].astype(np.int64, copy=False))
                        crp_col_parts.append(col_ids[m_c])

        qc_ok = (umi_l >= float(min_umi)) & (ng_l >= int(min_genes))
        qc_cols = np.flatnonzero(qc_ok)

        c_qc = None
        if proto_map is None:
            ok_init = ng_l >= int(min_genes)  # PREREG_F §3.2 & §3.5: ALL cells with >=200 genes, no %MT
            c_idx = np.concatenate(crp_idx_parts) if crp_idx_parts else np.empty(0, dtype=np.int32)
            c_dat = np.concatenate(crp_dat_parts) if crp_dat_parts else np.empty(0, dtype=np.int64)
            c_col = np.concatenate(crp_col_parts) if crp_col_parts else np.empty(0, dtype=np.int64)

            init_m = ok_init[c_col]
            if np.any(init_m):
                bins = np.clip(c_dat[init_m], 1, 65)
                np.add.at(hist, (c_idx[init_m], bins), 1)

            col_map_qc = np.full(n_col, -1, dtype=np.int64)
            col_map_qc[qc_cols] = np.arange(len(qc_cols), dtype=np.int64)
            qc_m = col_map_qc[c_col] >= 0
            c_qc = sp.csc_matrix(
                (c_dat[qc_m].astype(np.float64), (c_idx[qc_m], col_map_qc[c_col[qc_m]])),
                shape=(len(all_template_guides), len(qc_cols)),
            )

        channel_cache.append(
            {
                "ch_num": ch_num,
                "ch_tag": ch_tag,
                "path": hp,
                "n_feat": n_feat,
                "n_col": n_col,
                "rna_rows": rna_rows,
                "rids": [fids[i].split(".")[0] for i in rna_rows],
                "rnames": [fnames[i] for i in rna_rows],
                "barcodes": barcodes,
                "umi": umi_l,
                "qc_cols": qc_cols,
                "proto_map": proto_map,
                "c_qc": c_qc,
            }
        )

    thrs: Optional[np.ndarray] = None
    if guide_call_mode == "H5_CRISPR_POISSON_GAUSSIAN_LINE_FIT":
        thrs = fit_poisson_gaussian_thresholds_from_hist(hist, min_thr=5, max_thr=20, default_thr=7)

    # Pass 2: Assign single-sgRNA cells per channel (keep guide_id as pert_key for Addendum-B best-guide)
    n_qc_rna_cells = 0
    raw_kept_list: List[Dict[str, Any]] = []
    for ch in channel_cache:
        qc_cols = ch["qc_cols"]
        n_qc_rna_cells += len(qc_cols)
        barcodes = ch["barcodes"]
        keep_cols: List[int] = []
        kept_gids: List[str] = []

        if ch["proto_map"] is not None:
            pmap = ch["proto_map"]
            for raw_col in qc_cols:
                gid = pmap.get(barcodes[int(raw_col)])
                if gid is not None and (gid in ntc_guide_ids or gid in guide_to_target):
                    keep_cols.append(int(raw_col))
                    kept_gids.append(gid)
        else:
            c_qc = ch["c_qc"]
            assert thrs is not None and c_qc is not None
            t1, t2, a1 = extract_top2_csc(c_qc)
            pos_mask = c_qc.data >= thrs[c_qc.indices]
            pos_csc = sp.csc_matrix((pos_mask.astype(np.int8), c_qc.indices, c_qc.indptr), shape=c_qc.shape)
            n_pos = np.asarray(pos_csc.sum(axis=0)).ravel()
            s1_mask = (n_pos == 1) & (t1 >= 5.0) & ((t2 < 3.0) | (t1 / np.maximum(t2, 0.5) >= 3.0))
            for loc_j in np.flatnonzero(s1_mask):
                gid = all_template_guides[int(a1[loc_j])]
                if gid in ntc_guide_ids or gid in guide_to_target:
                    keep_cols.append(int(qc_cols[loc_j]))
                    kept_gids.append(gid)

        raw_kept_list.append(
            {
                **ch,
                "keep_arr": np.asarray(keep_cols, dtype=np.int64),
                "kept_gids": kept_gids,
            }
        )

    # Assemble across channels (to allow global Stage-2 split if Stage-1 missing AND global fallback)
    all_keep_nt: List[bool] = []
    all_keep_canon: List[str] = []
    all_keep_bcs: List[str] = []
    all_keep_batches: List[str] = []
    all_keep_halves: List[int] = []
    ch_slice_bounds: List[Tuple[int, int]] = []

    offset = 0
    for ch in raw_kept_list:
        n_k = len(ch["keep_arr"])
        half_v = 1 if (ch["ch_num"] is not None and ch["ch_num"] % 2 == 1) else 2
        for gid, raw_col in zip(ch["kept_gids"], ch["keep_arr"]):
            is_nt = gid in ntc_guide_ids
            all_keep_nt.append(is_nt)
            all_keep_canon.append("" if is_nt else guide_to_target[gid][1])
            all_keep_bcs.append(ch["barcodes"][int(raw_col)])
            all_keep_batches.append(ch["ch_tag"])
            all_keep_halves.append(half_v if use_stage1 else 0)
        ch_slice_bounds.append((offset, offset + n_k))
        offset += n_k

    nt_all_arr = np.asarray(all_keep_nt, dtype=bool)
    if len(nt_all_arr) == 0:
        raise RuntimeError(f"GSE327644 {line_name} yielded zero single-sgRNA QC-passing cells")

    if use_stage1:
        halves_all = np.asarray(all_keep_halves, dtype=np.int8)
        batches_all, n_nt_fallback_cells = _apply_min_ntc_batch_fallback(
            np.asarray(all_keep_batches, dtype=object),
            nt_all_arr,
            halves_all,
            group_labels=None,
            min_batch_nt=min_batch_nt,
        )
        batches_all = batches_all.astype(str)
    else:
        halves_all, batches_all = _apply_stage2_seed42_split(
            nt_mask=nt_all_arr,
            canon_targets=all_keep_canon,
            barcodes=all_keep_bcs,
        )
        n_nt_fallback_cells = 0

    parts: List[PartSpec] = []
    per_channel_summary: List[Dict[str, Any]] = []
    for ch, (s0, s1) in zip(raw_kept_list, ch_slice_bounds):
        keep_arr = ch["keep_arr"]
        nt_sub = nt_all_arr[s0:s1]
        batches_sub = batches_all[s0:s1]
        halves_sub = halves_all[s0:s1]

        pkeys: List[str] = []
        pmeta: Dict[str, Tuple[str, str, str, str]] = {}
        for gid, is_nt in zip(ch["kept_gids"], nt_sub):
            if is_nt:
                pkeys.append("")
            else:
                pkeys.append(gid)
                if gid not in pmeta:
                    pmeta[gid] = guide_to_target[gid]

        read_chunk, read_rows, read_chunk_cols = _make_int64_10x_csc_subset_readers(
            h5_path=ch["path"],
            keep_cols=keep_arr,
            n_raw_features=ch["n_feat"],
            rna_row_idx=ch["rna_rows"],
        )
        parts.append(
            PartSpec(
                tag=f"{line_name}_{ch['ch_tag']}",
                path=ch["path"],
                n_cells=len(keep_arr),
                n_genes_raw=len(ch["rna_rows"]),
                nt_mask=nt_sub,
                batches=batches_sub,
                half_flags=halves_sub,
                umi=ch["umi"][keep_arr],
                pert_keys=pkeys,
                pert_meta=pmeta,
                var_gene_ids=ch["rids"],
                var_gene_names=ch["rnames"],
                read_chunk_fn=read_chunk,
                read_rows_fn=read_rows,
                read_chunk_cols_fn=read_chunk_cols,
            )
        )
        per_channel_summary.append(
            {
                "channel": ch["ch_tag"],
                "ch_num": ch["ch_num"],
                "half_assigned": sorted({int(x) for x in halves_sub.tolist()}) if len(halves_sub) else [],
                "n_cells_raw": int(ch["n_col"]),
                "n_cells_rna_qc": int(len(ch["qc_cols"])),
                "n_cells_single_sgrna_kept": int(len(keep_arr)),
                "n_ntc_kept": int(nt_sub.sum()),
                "n_targeting_kept": int((~nt_sub).sum()),
            }
        )

    audit_meta = {
        "line": line_name,
        "feasibility": feas,
        "guide_call_mode": guide_call_mode,
        "split_branch": split_branch,
        "n_channels_found": len(files),
        "channel_numbers": parsed_ch_nums,
        "n_nt_fallback_cells": int(n_nt_fallback_cells),
        "min_umi": float(min_umi),
        "min_genes": int(min_genes),
        "max_pct_mt": None,
        "ntc_rule": "all_control_guides_NTC001_to_NTC593" if not tpl["distinguishes_safe_vs_ntc"] else "non_targeting_only",
        "guide_template_n_ntc": len(ntc_guide_ids),
        "guide_template_n_targeting": len(guide_to_target),
        "n_qc_rna_cells": int(n_qc_rna_cells),
        "n_kept_cells": int(len(nt_all_arr)),
        "n_kept_ntc": int(nt_all_arr.sum()),
        "n_ntc_half_a": int(np.sum(nt_all_arr & (halves_all == 1))),
        "n_ntc_half_b": int(np.sum(nt_all_arr & (halves_all == 2))),
        "median_guide_threshold": float(np.median(thrs)) if thrs is not None else None,
        "channels": per_channel_summary,
    }
    return parts, audit_meta


def standardize_gse327644_to_audit_line(
    line_name: str,
    raw_dir: Path,
    cache_out_root: Path,
    aligned_genes: Optional[Sequence[str]] = None,
    hgnc_tsv: Path = DEFAULT_HGNC_TSV,
    guide_template_path: Optional[Path] = None,
    min_gene_overlap: int = 100,
    n_workers: int = 4,
) -> Tuple[AuditLineData, Dict[str, Any]]:
    """
    End-to-end clean-room raw ingestion + Addendum-B standardization for `mac_ctrl` / `mac_lps`
    from real or synthetic Cell Ranger Flex `.h5` files.
    """
    hgnc_prep = build_hgnc_lookup(hgnc_tsv)
    parts, meta = build_gse327644_audit_parts(
        line_name=line_name,
        raw_dir=raw_dir,
        hgnc=hgnc_prep,
        guide_template_path=guide_template_path,
    )
    line_dir = process_r10_line(
        line_name=line_name,
        parts=parts,
        hgnc=hgnc_prep,
        chunk_size=4096,
        out_root=cache_out_root,
        extra_meta={"audit_reader": meta},
        n_workers=n_workers,
    )
    hgnc_canon = CanonicalHGNC.from_tsv(hgnc_tsv)
    adata = load_line_dataset(
        line_dir=line_dir,
        representation="log1p",
        aligned_genes=aligned_genes,
        hgnc=hgnc_canon,
        min_gene_overlap=min_gene_overlap,
    )
    return adata, meta


def resolve_line_cache_dir_r13(line_name: str, phase: str) -> Path:
    """Resolve standardized cache directory cleanly under current blinding phase."""
    if phase.startswith("dev_"):
        if line_name in R9_LINE_SET:
            return R9_LINES_DIR / line_name
        if line_name in R10_ALL_LINES:
            return R10_LINES_DIR / line_name
        if line_name in R11_ALL_LINES:
            return R11_LINES_DIR / line_name
        if line_name in R12_ALL_LINES:
            return R12_LINES_DIR / line_name
        raise RuntimeError(f"Dev phase {phase} cannot access Round 13 line {line_name!r}.")
    if line_name in R9_LINE_SET:
        return R9_LINES_DIR / line_name
    if line_name in R10_SOURCE_EXT_NAMES:
        return R10_LINES_DIR / line_name
    if line_name in R12_ALL_LINES:
        return R12_LINES_DIR / line_name
    return R13_LINES_DIR / line_name
