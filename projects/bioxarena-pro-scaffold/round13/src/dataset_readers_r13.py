#!/usr/bin/env python3
"""
Round 13 Dataset Ingestion Adapter for `GSE327644` (`PREREG_F` Stage 2).

Reads Cell Ranger Flex `.h5` matrices (`Gene Expression` + `CRISPR Guide Capture`),
reusing frozen Round 10 / Round 11 preprocessing primitives cleanly without modifying any
file under `round10/`, `round11/`, or `round12/`:

  1. Lines:
     - `mac_ctrl`: unstimulated control macrophage Flex channels (`L01`–`L04`; confirmatory candidate).
     - `mac_lps`: LPS-stimulated macrophage Flex channels (`L01`–`L04`; secondary descriptive only).
  2. Strict Pre-Metric Feasibility Gate (`check_gse327644_feasibility` / `EXCLUDED_FEASIBILITY`):
     - Verifies non-negative integer raw counts in `/matrix/data` and `/matrix` CSC layout.
     - Verifies valid single-cell protospacer calling (`protospacer_calls_per_cell.csv` or
       `CRISPR Guide Capture` integer UMI submatrix) plus non-empty NTC & HGNC target matches.
  3. Filtering & NTC Definition (`PREREG_F` Rule 3):
     - Strictly Cell Ranger `filtered` `.h5` matrices (when present alongside `raw`).
     - Single-sgRNA cells: companion Cell Ranger `protospacer_calls_per_cell` `num_features == 1`
       when present, else Cell Ranger Flex Poisson-Gaussian bimodal EM (`fit_poisson_gaussian_thresholds_from_hist`)
       with `n_pos == 1` and 3x dominance.
     - Deterministic NTC rule: if guide template distinguishes `non-targeting` from `safe-targeting`,
       ONLY `non-targeting` counts as NTC (excluding `safe-targeting`); otherwise ALL control
       guides count as NTC.
     - Cell RNA QC: `umi >= 300` AND `ngenes >= 200`, with NO mitochondrial percentage cap
       (`max_pct_mt=None`; Flex panels may omit mitochondrial transcripts).
  4. Half Split & Unit Centering (`PREREG_F` Rule 4):
     - Stage 1: Channel parity (`L01`, `L03` -> Half A `1`; `L02`, `L04` -> Half B `2`; odd->1,
       even->2); centering unit/batch = channel (`L01`..`L04`), falling back to `half_{h}_pooled`
       if a channel has `< 15` NTCs.
     - Stage 2 fallback (only when channel labels are absent or single-half): locked Ruling 2
       `np.random.default_rng(42)` stratified 50/50 split on barcode-sorted indices (`half_1` / `half_2`).
  5. Integration: outputs `PartSpec` list cleanly into frozen `process_r10_line`.
"""

from __future__ import annotations

import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROUND13_SRC_DIR = Path(__file__).resolve().parent
REPO_ROOT = ROUND13_SRC_DIR.parents[1]
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
R10_SCRIPTS_DIR = REPO_ROOT / "round10" / "shared" / "scripts"
R11_SRC_DIR = REPO_ROOT / "round11" / "src"
for p in (str(R10_SRC_DIR), str(R10_SCRIPTS_DIR), str(R11_SRC_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from data_loader_r10 import decode_str_array  # noqa: E402
from dataset_readers_r11 import (  # noqa: E402
    _apply_min_ntc_batch_fallback,
    _is_plain_single_gene_symbol,
    _pick_column,
)
from prepare_r10_lines import (  # noqa: E402
    HGNCLookup,
    PartSpec,
    build_hgnc_lookup,
    compute_10x_csc_rna_qc,
    extract_top2_csc,
    fit_poisson_gaussian_thresholds_from_hist,
    make_10x_csc_subset_readers,
    process_r10_line,
)

_LANE_RE = re.compile(r"(?:^|[_\-/])(?:Mac[_-]?)?L0*(\d+)(?:[_\-.]|$)", re.IGNORECASE)
_LANE_FALLBACK_RE = re.compile(r"(?:lane|channel|ch)[_-]?0*(\d+)(?:[_\-.]|$)", re.IGNORECASE)
_MULTI_GUIDE_RE = re.compile(r"[,;|+]")
_UNASSIGNED_TOKENS = frozenset({"", "none", "nan", "na", "null", "unassigned", "multiplet", "doublet"})
_EXPLICIT_NONTARGET_TOKENS = frozenset({"non-targeting", "non_targeting", "nontargeting", "non-target", "ntc"})
_GENERIC_CONTROL_TOKENS = frozenset(
    {"control", "ctrl", "ntc", "nt", "neg", "negative", "negative_control", "neg_ctrl", "scramble", "scrambled", "false"}
)
_SAFE_RE = re.compile(r"(?:safe[-_ ]?target|^0?safe(?:[-_]|$))", re.IGNORECASE)
_JUN_EXCEL_RE = re.compile(r"^([1-4])-jun$", re.IGNORECASE)
_JUN_GUIDE_RE = re.compile(r"^JUN-([1-4])$", re.IGNORECASE)


def _resolve_raw_dir(raw_dir: Path, line: str) -> Path:
    base = Path(raw_dir)
    for cand in (base / "gse327644", base / line, base):
        if cand.exists():
            return cand
    return base / "gse327644"


def _extract_gse327644_lane_tag(path: Path) -> Tuple[Optional[str], Optional[int]]:
    """Extract channel ID (`L01`..`L08`) and 1-based channel integer from path stem/parents."""
    for text in (path.stem, path.parent.name, path.name):
        for pat in (_LANE_RE, _LANE_FALLBACK_RE):
            m = pat.search(text)
            if m:
                num = int(m.group(1))
                return f"L{num:02d}", num
    return None, None


def _discover_gse327644_h5_files(search_dir: Path, line: str) -> Tuple[List[Path], Optional[Path]]:
    """
    Locate strictly filtered Cell Ranger Flex `.h5` files matching `line` (`mac_ctrl` vs `mac_lps`)
    and optional guide template metadata file under `search_dir`.
    """
    if line not in ("mac_ctrl", "mac_lps"):
        raise ValueError(f"Unsupported Round 13 GSE327644 line: {line!r} (expected 'mac_ctrl' or 'mac_lps')")

    all_h5 = sorted(
        p for p in search_dir.rglob("*.h5")
        if p.is_file() and not p.name.startswith("._") and not p.name.startswith(".")
    )
    filtered_h5 = [p for p in all_h5 if "raw_feature_bc_matrix" not in p.name.lower()]
    pool_h5 = filtered_h5 if filtered_h5 else all_h5

    want_tok = "ctrl" if line == "mac_ctrl" else "lps"
    other_tok = "lps" if line == "mac_ctrl" else "ctrl"

    cond_h5 = [
        p for p in pool_h5
        if want_tok in f"{p.parent.name}/{p.name}".lower()
        and other_tok not in f"{p.parent.name}/{p.name}".lower()
    ]
    if not cond_h5:
        if (search_dir / line).exists():
            sub_h5 = sorted(
                p for p in (search_dir / line).rglob("*.h5")
                if p.is_file() and not p.name.startswith("._") and "raw_feature_bc_matrix" not in p.name.lower()
            )
            cond_h5 = sub_h5
        elif search_dir.name == line and not any(other_tok in p.name.lower() for p in pool_h5):
            cond_h5 = pool_h5

    if not cond_h5:
        raise FileNotFoundError(f"No Cell Ranger Flex .h5 files found for {line!r} under {search_dir}")

    tmpl_cands: List[Path] = []
    for pat in (
        "*guide_template*",
        "*guide*meta*",
        "*feature_ref*",
        "*sgrna*table*",
        "*sgrna*.csv*",
        "*sgrna*.tsv*",
        "*sgrna*.txt*",
    ):
        for p in sorted(search_dir.rglob(pat)):
            if p.is_file() and not p.name.startswith("._") and p not in tmpl_cands:
                tmpl_cands.append(p)
    return cond_h5, (tmpl_cands[0] if tmpl_cands else None)


def _find_companion_protospacer_csv(h5_path: Path, search_dir: Path) -> Optional[Path]:
    """Locate optional Cell Ranger `protospacer_calls_per_cell.csv(.gz)` matching `h5_path`."""
    candidates: List[Path] = []
    for base in (h5_path.parent, h5_path.parent / "crispr_analysis", h5_path.parent.parent / "crispr_analysis"):
        if base.exists():
            for pat in ("protospacer_calls_per_cell.csv*", "*protospacer_calls*.csv*"):
                candidates.extend(sorted(p for p in base.glob(pat) if p.is_file() and not p.name.startswith("._")))
    if not candidates:
        return None
    if len(candidates) == 1 and candidates[0].parent in (h5_path.parent, h5_path.parent / "crispr_analysis"):
        return candidates[0]

    stem_low = h5_path.stem.replace("_sample_filtered_feature_bc_matrix", "").replace("_filtered_feature_bc_matrix", "").lower()
    for p in candidates:
        if stem_low and stem_low in p.name.lower():
            return p
    lane_tag, _ = _extract_gse327644_lane_tag(h5_path)
    want_cond = "ctrl" if "ctrl" in h5_path.name.lower() else ("lps" if "lps" in h5_path.name.lower() else "")
    for p in candidates:
        pl = p.name.lower()
        if (not lane_tag or lane_tag.lower() in pl) and (not want_cond or want_cond in pl):
            return p
    return candidates[0] if candidates[0].name.startswith("protospacer_calls_per_cell") else None


def _is_safe_targeting(type_str: str, gid: str, tgt: str) -> bool:
    for val in (type_str, gid, tgt):
        if val and _SAFE_RE.search(val):
            return True
    return type_str.strip().lower() in {"safe", "safe_targeting", "safe-targeting", "safe_ctrl", "safe-ctrl"}


def _read_tabular_auto(path: Path) -> pd.DataFrame:
    opener = gzip.open if str(path).endswith(".gz") else open
    first_data = ""
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if s and not s.startswith("#"):
                first_data = line
                break
    sep = "\t" if "\t" in first_data else ","
    return pd.read_csv(path, sep=sep, comment="#", dtype=str, keep_default_na=False)


def _clean_meta_token(val: str) -> str:
    s = str(val).strip()
    if s.lower() in {"", "na", "n/a", "nan", "none", "null"}:
        return ""
    return s


def _extract_target_stem_r13(gid: str, gn: str, raw_tgt: str) -> str:
    """
    Derive target symbol stem strictly per `PREREG_F §3.4`:
      - `stem = guide_id.rsplit('-', 1)[0]` (prevents hyphen truncation such as `HLA-A-1 -> HLA`
        or `PPAN-P2RY11-1 -> PPAN`).
      - Remaps Excel-converted `1-Jun..4-Jun` (where `guide_target_gene_name == 'JUN'` with
        `stem in {'1','2','3','4'}` OR `.h5` `1-Jun..4-Jun` feature names) to `'JUN'`.
    """
    gid_s = gid.strip()
    gn_s = gn.strip()
    raw_tgt_s = _clean_meta_token(raw_tgt)

    if _JUN_EXCEL_RE.match(gid_s) or _JUN_EXCEL_RE.match(gn_s):
        return "JUN"

    guide_src = gid_s if ("-" in gid_s or not gn_s or "-" not in gn_s) else gn_s
    if "-" in guide_src:
        stem, tail = guide_src.rsplit("-", 1)
        stem = stem.strip()
        tail = tail.strip()
        if raw_tgt_s.upper() == "JUN" and stem in {"1", "2", "3", "4"}:
            return "JUN"
        if "_" not in tail:
            return stem
        return raw_tgt_s if raw_tgt_s else guide_src.rsplit("_", 1)[0].strip()

    if raw_tgt_s:
        return raw_tgt_s
    prefix = guide_src.split("_")[0].strip()
    return prefix if _is_plain_single_gene_symbol(prefix) else ""


def _classify_guide_rows(
    rows: Sequence[Tuple[str, str, str, str, str]],
    hgnc: HGNCLookup,
) -> Dict[str, Any]:
    """
    Classify `(gid, gname, gtype, raw_tgt, raw_ens)` tuples under `PREREG_F` Rule 3 & §3.4:
      - If ANY row distinguishes `safe-targeting`, ONLY `non-targeting` rows qualify as NTC;
        `safe-targeting` rows are excluded strictly.
      - Otherwise ALL control rows qualify as NTC.
      - Targeting guides derive `stem = guide_id.rsplit('-', 1)[0]` (with `1-Jun..4-Jun -> JUN`)
        and resolve `(stem, target_gene_id)` via `HGNCLookup.resolve`.
    """
    has_safe = any(_is_safe_targeting(gt, gid, rt) or _is_safe_targeting(gt, gn, rt) for gid, gn, gt, rt, _ in rows)
    ntc_keys: Set[str] = set()
    safe_keys: Set[str] = set()
    guide_to_target: Dict[str, Tuple[str, str, str, str]] = {}

    for gid, gn, gt, raw_tgt, raw_ens in rows:
        keys = [k for k in (gid, gn) if k]
        if not keys:
            continue
        if _is_safe_targeting(gt, gid, raw_tgt) or _is_safe_targeting(gt, gn, raw_tgt):
            safe_keys.update(keys)
            continue

        gt_l = gt.strip().lower()
        rt_l = raw_tgt.strip().lower()
        gid_l = gid.strip().lower()
        gn_l = gn.strip().lower()

        is_explicit_nt = (
            gt_l in _EXPLICIT_NONTARGET_TOKENS
            or "non-target" in gt_l
            or "nontarget" in gt_l
            or "non_target" in gt_l
            or "non-target" in rt_l
            or "nontarget" in rt_l
            or "non-target" in gid_l
            or "nontarget" in gid_l
            or "non-target" in gn_l
            or "nontarget" in gn_l
        )
        is_explicit_role_targeting = (
            "target" in gt_l and "non" not in gt_l and "safe" not in gt_l and "untarget" not in gt_l
        )
        is_generic_ctrl = (
            not is_explicit_role_targeting
            and (
                gt_l in _GENERIC_CONTROL_TOKENS
                or rt_l in _GENERIC_CONTROL_TOKENS
                or gid_l.startswith(("ntc", "ctrl", "control", "neg"))
                or gn_l.startswith(("ntc", "ctrl", "control", "neg"))
            )
        )
        if is_explicit_nt or (not has_safe and is_generic_ctrl):
            ntc_keys.update(keys)
            continue

        raw_ens_s = _clean_meta_token(raw_ens).split(".")[0]
        cand_sym = _extract_target_stem_r13(gid, gn, raw_tgt)

        if cand_sym and _is_plain_single_gene_symbol(cand_sym):
            canon_sym, hstat = hgnc.resolve(cand_sym, raw_ens_s)
            if canon_sym in hgnc.approved:
                ens = raw_ens_s or hgnc.sym_to_ens.get(canon_sym, "")
                info = (ens, canon_sym, cand_sym, hstat)
                for k in keys:
                    guide_to_target[k] = info
                    if canon_sym == "JUN":
                        m_ex = _JUN_EXCEL_RE.match(k)
                        m_ju = _JUN_GUIDE_RE.match(k)
                        if m_ex:
                            guide_to_target.setdefault(f"JUN-{m_ex.group(1)}", info)
                        elif m_ju:
                            guide_to_target.setdefault(f"{m_ju.group(1)}-Jun", info)

    rule_tag = (
        "only_nontargeting_excluding_safetargeting"
        if has_safe
        else "all_controls_no_safetargeting_distinction"
    )
    return {
        "ntc_keys": ntc_keys,
        "safe_keys": safe_keys,
        "guide_to_target": guide_to_target,
        "has_safe_targeting": has_safe,
        "ntc_rule_used": rule_tag,
    }


def load_gse327644_guide_metadata(
    template_path: Optional[Path],
    sample_h5_path: Optional[Path],
    hgnc: HGNCLookup,
) -> Dict[str, Any]:
    """Load guide-to-target and NTC/safe-targeting mapping from template file or `.h5` features."""
    raw_rows: List[Tuple[str, str, str, str, str]] = []
    source_desc = ""

    if template_path is not None and template_path.exists():
        df = _read_tabular_auto(template_path)
        c_id = _pick_column(df.columns, ["guide_id", "sgRNA", "sgrna", "id", "feature_id", "guide_name", "name"])
        c_nm = _pick_column(df.columns, ["name", "guide_name", "sgRNA", "sgrna", "guide_id", "id"])
        c_tp = _pick_column(
            df.columns,
            ["guide_type", "guide_role", "role", "type", "category", "class", "feature_category", "control_type", "is_control", "targeting"],
        )
        c_tg = _pick_column(
            df.columns,
            ["guide_target_gene_name", "target_gene_name", "target_gene", "gene_symbol", "gene_name", "target", "gene", "symbol"],
        )
        c_en = _pick_column(df.columns, ["guide_target_gene_id", "target_gene_id", "ensembl_id", "gene_id"])
        if c_id is None and c_nm is None:
            raise KeyError(f"Guide template {template_path} missing guide identifier column: {list(df.columns)}")
        for i in range(len(df)):
            gid = str(df[c_id].iloc[i]).strip() if c_id else ""
            gn = str(df[c_nm].iloc[i]).strip() if c_nm else gid
            gt = str(df[c_tp].iloc[i]).strip() if c_tp else ""
            rt = str(df[c_tg].iloc[i]).strip() if c_tg else ""
            re_id = str(df[c_en].iloc[i]).strip().split(".")[0] if c_en else ""
            raw_rows.append((gid or gn, gn or gid, gt, rt, re_id))
        source_desc = f"template:{template_path.name}"
    elif sample_h5_path is not None and sample_h5_path.exists():
        with h5py.File(sample_h5_path, "r") as f:
            feat = f["matrix"]["features"]
            ftypes = np.array(decode_str_array(feat["feature_type"][:]))
            fids = decode_str_array(feat["id"][:])
            fnames = decode_str_array(feat["name"][:])
            f_tg = decode_str_array(feat["target_gene_name"][:]) if "target_gene_name" in feat else [""] * len(fids)
            f_te = decode_str_array(feat["target_gene_id"][:]) if "target_gene_id" in feat else [""] * len(fids)
            f_gt = decode_str_array(feat["guide_type"][:]) if "guide_type" in feat else [""] * len(fids)
            for idx, ft in enumerate(ftypes):
                if "crispr" in ft.lower() or "guide" in ft.lower():
                    raw_rows.append((fids[idx], fnames[idx], f_gt[idx], f_tg[idx], f_te[idx].split(".")[0]))
        source_desc = f"h5_features:{sample_h5_path.name}"
    else:
        raise FileNotFoundError("Missing both guide_template file and sample .h5 for GSE327644 guide metadata")

    meta = _classify_guide_rows(raw_rows, hgnc)
    meta["source"] = source_desc
    return meta


def _sample_h5_nonneg_int(ds: h5py.Dataset, max_elems: int = 2000) -> bool:
    n = min(int(ds.shape[0]), max_elems)
    if n == 0:
        return True
    arr = np.asarray(ds[:n], dtype=np.float64)
    return bool(np.all(np.isfinite(arr) & (arr >= 0.0) & (np.abs(arr - np.round(arr)) <= 1e-6)))


def _is_nonhuman_or_vector_var(sym: str, gid: str, genome: str = "") -> bool:
    su = sym.upper()
    gu = gid.upper()
    geu = genome.upper()
    if gu.startswith("ENSMUSG") or geu in ("CAS9", "MM10", "GRCM38", "GRCM39"):
        return True
    if su in ("CAS9", "DCAS9", "DCAS9-KRAB-MECP2", "KRAB", "MECP2_LENTIVIRUS") or "DCAS9" in su:
        return True
    return False


def _read_protospacer_calls_csv(csv_path: Path) -> Dict[str, Optional[str]]:
    """
    Parse Cell Ranger `protospacer_calls_per_cell.csv(.gz)`:
    keep strictly single-guide calls (`num_features == 1` & non-delimited `feature_call`).
    """
    opener = gzip.open if str(csv_path).endswith(".gz") else open
    out: Dict[str, Optional[str]] = {}
    with opener(csv_path, "rt", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames:
            return out
        bc_c = _pick_column(reader.fieldnames, ["cell_barcode", "barcode", "cell_id"])
        nf_c = _pick_column(reader.fieldnames, ["num_features", "n_features", "num_guides"])
        fc_c = _pick_column(reader.fieldnames, ["feature_call", "call", "guide_call", "protospacer_call"])
        if bc_c is None or fc_c is None:
            raise KeyError(f"protospacer_calls CSV {csv_path} missing cell_barcode/feature_call columns: {reader.fieldnames}")
        for row in reader:
            bc = str(row[bc_c]).strip()
            call_raw = str(row[fc_c]).strip()
            if nf_c is not None:
                try:
                    if int(float(row[nf_c])) != 1:
                        out[bc] = None
                        continue
                except Exception:
                    out[bc] = None
                    continue
            if not call_raw or call_raw.lower() in _UNASSIGNED_TOKENS or _MULTI_GUIDE_RE.search(call_raw):
                out[bc] = None
            else:
                out[bc] = call_raw
    return out


def check_gse327644_feasibility(
    raw_dir: Path,
    line: str = "mac_ctrl",
    hgnc: Optional[HGNCLookup] = None,
) -> Dict[str, Any]:
    """
    `PREREG_F` Rule 5 pre-metric feasibility gate (metadata & counts ONLY; never computes any
    expression correlation or perturbation response):
      1) Checks filtered `.h5` files exist for `line` (`mac_ctrl` / `mac_lps`);
      2) Verifies `/matrix` CSC layout & non-negative integer counts;
      3) Verifies `Gene Expression` features & per-cell guide calling source (`protospacer_calls`
         CSV OR non-empty `CRISPR Guide Capture` UMIs);
      4) Verifies NTC & canonical HGNC targeting guide overlap and NTC rule classification.
    """
    try:
        hgnc_l = hgnc or build_hgnc_lookup()
        search_dir = _resolve_raw_dir(raw_dir, line)
        h5_files, tmpl_path = _discover_gse327644_h5_files(search_dir, line)
        guide_meta = load_gse327644_guide_metadata(tmpl_path, h5_files[0], hgnc_l)
        ntc_keys: Set[str] = guide_meta["ntc_keys"]
        guide_to_target: Dict[str, Tuple[str, str, str, str]] = guide_meta["guide_to_target"]

        total_cells = 0
        n_rna_ref = n_crp_ref = 0
        n_csv_matched_files = 0
        lanes_detected: Set[str] = set()
        crp_nnz_sample = 0

        for lp in h5_files:
            ltag, _ = _extract_gse327644_lane_tag(lp)
            if ltag:
                lanes_detected.add(ltag)
            with h5py.File(lp, "r") as f:
                if "matrix" not in f:
                    return {"feasible": False, "reason": f"Missing /matrix group in {lp.name}"}
                m = f["matrix"]
                for req in ("shape", "barcodes", "data", "indices", "indptr", "features"):
                    if req not in m:
                        return {"feasible": False, "reason": f"Missing /matrix/{req} in {lp.name}"}
                if not _sample_h5_nonneg_int(m["data"]):
                    return {"feasible": False, "reason": f"Non-integer or negative counts in {lp.name}/matrix/data"}
                feat = m["features"]
                for req_f in ("feature_type", "id", "name"):
                    if req_f not in feat:
                        return {"feasible": False, "reason": f"Missing /matrix/features/{req_f} in {lp.name}"}
                n_feat, n_col = int(m["shape"][0]), int(m["shape"][1])
                total_cells += n_col
                ftypes = np.array(decode_str_array(feat["feature_type"][:]))
                fids = np.array(decode_str_array(feat["id"][:]))
                fnames = np.array(decode_str_array(feat["name"][:]))
                fgenomes = np.array(decode_str_array(feat["genome"][:])) if "genome" in feat else np.full(n_feat, "")

                rna_idx = [
                    i for i, ft in enumerate(ftypes)
                    if ft == "Gene Expression" and not _is_nonhuman_or_vector_var(fnames[i], fids[i], fgenomes[i])
                ]
                crp_idx = [i for i, ft in enumerate(ftypes) if "crispr" in ft.lower() or "guide" in ft.lower()]
                if len(rna_idx) == 0:
                    return {"feasible": False, "reason": f"Zero human Gene Expression features in {lp.name}"}
                n_rna_ref = len(rna_idx)
                n_crp_ref = len(crp_idx)

                csv_p = _find_companion_protospacer_csv(lp, search_dir)
                if csv_p is not None:
                    n_csv_matched_files += 1
                if crp_idx:
                    crp_set = set(crp_idx)
                    idx_s = np.asarray(m["indices"][: min(int(m["indices"].shape[0]), 200000)], dtype=np.int64)
                    crp_nnz_sample += sum(1 for v in idx_s if int(v) in crp_set)

                overlap_ntc = sum(1 for i in crp_idx if fids[i] in ntc_keys or fnames[i] in ntc_keys)
                overlap_tgt = sum(1 for i in crp_idx if fids[i] in guide_to_target or fnames[i] in guide_to_target)

        has_guide_calling = (n_csv_matched_files == len(h5_files)) or (n_crp_ref > 0 and crp_nnz_sample > 0)
        if not has_guide_calling:
            return {
                "feasible": False,
                "reason": "Missing per-cell protospacer calls CSV and zero CRISPR Guide Capture UMIs",
            }
        ok = bool(n_rna_ref > 0 and (overlap_ntc > 0 or len(ntc_keys) > 0) and (overlap_tgt > 0 or len(guide_to_target) > 0))
        return {
            "feasible": ok,
            "reason": "OK" if ok else "No resolvable NTC or targeting guides matched features",
            "line": line,
            "n_h5_files": len(h5_files),
            "total_cells_raw": int(total_cells),
            "n_rna_features": int(n_rna_ref),
            "n_crispr_features": int(n_crp_ref),
            "guide_call_source": (
                "protospacer_calls_per_cell_csv" if n_csv_matched_files == len(h5_files) else "h5_crispr_umi_em"
            ),
            "guide_meta_source": guide_meta["source"],
            "has_safe_targeting": bool(guide_meta["has_safe_targeting"]),
            "ntc_rule_used": guide_meta["ntc_rule_used"],
            "n_ntc_guides_meta": len(ntc_keys),
            "n_safe_guides_excluded": len(guide_meta["safe_keys"]),
            "n_target_guides_meta": len(guide_to_target),
            "lanes_detected": sorted(lanes_detected),
        }
    except Exception as exc:
        return {"feasible": False, "reason": f"{type(exc).__name__}: {exc}"}


def _assign_single_guide_from_h5_em(
    lane_records: Sequence[Dict[str, Any]],
    crp_names_ref: Sequence[str],
    min_guide_thr: int,
    max_guide_thr: int,
    default_guide_thr: int,
) -> np.ndarray:
    """Fit bimodal mixture thresholds on pooled `ngenes >= 200` CRISPR counts per `PREREG_F §3.2`."""
    n_crp = len(crp_names_ref)
    hist = np.zeros((n_crp, 66), dtype=np.int64)
    for rec in lane_records:
        if rec.get("fit_hist") is not None:
            hist += rec["fit_hist"]
        else:
            C_qc = rec["C_qc"]
            if C_qc.nnz > 0:
                r_idx, _ = C_qc.nonzero()
                np.add.at(hist, (r_idx, np.clip(C_qc.data.astype(np.int64, copy=False), 1, 65)), 1)
    return fit_poisson_gaussian_thresholds_from_hist(
        hist,
        min_thr=min_guide_thr,
        max_thr=max_guide_thr,
        default_thr=default_guide_thr,
    )


def build_gse327644_parts(
    raw_dir: Path,
    hgnc: HGNCLookup,
    line: str = "mac_ctrl",
    min_umi: int = 300,
    min_ngenes: int = 200,
    # PREREG_F Rule 3: NO mitochondrial% cutoff (Flex probe panels may omit mitochondrial transcripts).
    max_pct_mt: Optional[float] = None,
    min_batch_nt: int = 15,
    min_guide_thr: int = 5,
    default_guide_thr: int = 7,
    batch_cols: int = 20000,
    force_seed42_split: bool = False,
) -> Tuple[List[PartSpec], dict]:
    """
    Build `PartSpec` shards for `GSE327644` (`mac_ctrl` or `mac_lps`) Cell Ranger Flex `.h5`:
      - Feasibility gate enforced upfront (`EXCLUDED_FEASIBILITY` on failure).
      - RNA QC: `umi >= 300`, `ngenes >= 200`, `max_pct_mt=None`.
      - Single-sgRNA selection: Cell Ranger `protospacer_calls_per_cell.csv` (`num_features == 1`)
        when present alongside `.h5`, else pooled bimodal mixture (`n_pos == 1`, `dominance >= 3`).
      - Deterministic NTC rule (`non-targeting` only if `safe-targeting` exists, else all controls).
      - Stage 1: channel parity (`L01/L03 -> 1`, `L02/L04 -> 2`) and channel NTC subtraction
        (`< 15` NTC falls back to `half_{h}_pooled`); Stage 2 fallback: Ruling 2 `seed=42`.
    """
    feas = check_gse327644_feasibility(raw_dir=raw_dir, line=line, hgnc=hgnc)
    if not feas.get("feasible", False):
        raise RuntimeError(f"EXCLUDED_FEASIBILITY: {line} failed GSE327644 feasibility gate: {feas.get('reason')}")

    search_dir = _resolve_raw_dir(raw_dir, line)
    h5_files, tmpl_path = _discover_gse327644_h5_files(search_dir, line)
    guide_meta = load_gse327644_guide_metadata(tmpl_path, h5_files[0], hgnc)
    ntc_keys: Set[str] = guide_meta["ntc_keys"]
    safe_keys: Set[str] = guide_meta["safe_keys"]
    guide_to_target: Dict[str, Tuple[str, str, str, str]] = guide_meta["guide_to_target"]

    csv_paths = [_find_companion_protospacer_csv(lp, search_dir) for lp in h5_files]
    using_csv_calls = all(p is not None for p in csv_paths)
    lane_records: List[Dict[str, Any]] = []
    crp_ref_keys: Optional[List[str]] = None
    n_rna_qc_excl = 0

    for lp, csv_p in zip(h5_files, csv_paths):
        ltag, lnum = _extract_gse327644_lane_tag(lp)

        with h5py.File(lp, "r") as f:
            m = f["matrix"]
            n_feat, n_col = int(m["shape"][0]), int(m["shape"][1])
            barcodes = decode_str_array(m["barcodes"][:])
            ftypes = np.array(decode_str_array(m["features/feature_type"][:]))
            fids = np.array(decode_str_array(m["features/id"][:]))
            fnames = np.array(decode_str_array(m["features/name"][:]))
            fgenomes = np.array(decode_str_array(m["features/genome"][:])) if "genome" in m["features"] else np.full(n_feat, "")

            rna_rows = np.array(
                [
                    i for i, ft in enumerate(ftypes)
                    if ft == "Gene Expression" and not _is_nonhuman_or_vector_var(fnames[i], fids[i], fgenomes[i])
                ],
                dtype=np.int64,
            )
            crp_rows = np.array([i for i, ft in enumerate(ftypes) if "crispr" in ft.lower() or "guide" in ft.lower()], dtype=np.int64)
            rids = [fids[i].split(".")[0] or fnames[i] for i in rna_rows]
            rnames = [fnames[i] for i in rna_rows]
            mt_raw = np.zeros(n_feat, dtype=bool)
            mt_raw[rna_rows] = [str(x).startswith("MT-") for x in rnames]

        crp_keys = [
            fids[i] if (fids[i] in ntc_keys or fids[i] in safe_keys or fids[i] in guide_to_target) else fnames[i]
            for i in crp_rows
        ]
        if crp_ref_keys is None:
            crp_ref_keys = crp_keys

        umi_all, ng_all, pct_mt_all = compute_10x_csc_rna_qc(lp, rna_rows, mt_raw, batch_cols=batch_cols)
        # PREREG_F §3.2: fallback guide threshold EM fits on ALL cells with `ngenes >= 200` within each line.
        ok_fit = ng_all >= int(min_ngenes)
        # PREREG_F §3.5: retained single-sgRNA cells require BOTH `umi >= 300` AND `ngenes >= 200`.
        ok_qc = (umi_all >= float(min_umi)) & (ng_all >= int(min_ngenes))
        if max_pct_mt is not None:
            ok_qc &= pct_mt_all <= float(max_pct_mt)
        n_rna_qc_excl += int(np.sum(~ok_qc))
        qc_cols = np.where(ok_qc)[0]

        C_qc = sp.csc_matrix((len(crp_rows), len(qc_cols)), dtype=np.float64)
        fit_hist = np.zeros((len(crp_rows), 66), dtype=np.int64)
        if not using_csv_calls and len(crp_rows) > 0 and (len(qc_cols) > 0 or bool(np.any(ok_fit))):
            crp_remap = np.full(n_feat, -1, dtype=np.int32)
            crp_remap[crp_rows] = np.arange(len(crp_rows), dtype=np.int32)
            is_crp = crp_remap >= 0
            c_idx_l, c_dat_l, c_col_l = [], [], []
            col_map = np.full(n_col, -1, dtype=np.int32)
            col_map[qc_cols] = np.arange(len(qc_cols), dtype=np.int32)
            with h5py.File(lp, "r") as f:
                m = f["matrix"]
                indptr_all = np.asarray(m["indptr"][:], dtype=np.int64)
                for cs in range(0, n_col, batch_cols):
                    ce = min(cs + batch_cols, n_col)
                    p0, p1 = int(indptr_all[cs]), int(indptr_all[ce])
                    if p1 == p0:
                        continue
                    dat = m["data"][p0:p1]
                    idx = m["indices"][p0:p1]
                    clens = np.diff(indptr_all[cs : ce + 1])
                    cids = np.repeat(np.arange(cs, ce, dtype=np.int32), clens)
                    fit_m = is_crp[idx] & ok_fit[cids] & (dat > 0)
                    if np.any(fit_m):
                        np.add.at(
                            fit_hist,
                            (crp_remap[idx[fit_m]], np.clip(dat[fit_m].astype(np.int64, copy=False), 1, 65)),
                            1,
                        )
                    mapped = col_map[cids]
                    mc = is_crp[idx] & (mapped >= 0) & (dat > 0)
                    if np.any(mc):
                        c_idx_l.append(crp_remap[idx[mc]])
                        c_dat_l.append(dat[mc].astype(np.float64, copy=False))
                        c_col_l.append(mapped[mc])
            if c_idx_l:
                C_qc = sp.csr_matrix(
                    (np.concatenate(c_dat_l), (np.concatenate(c_idx_l), np.concatenate(c_col_l))),
                    shape=(len(crp_rows), len(qc_cols)),
                ).tocsc()

        lane_records.append(
            {
                "path": lp,
                "lane_tag": ltag,
                "lane_num": lnum,
                "n_feat": n_feat,
                "n_col": n_col,
                "barcodes": barcodes,
                "rna_rows": rna_rows,
                "rids": rids,
                "rnames": rnames,
                "crp_keys": crp_keys,
                "umi_all": umi_all,
                "qc_cols": qc_cols,
                "C_qc": C_qc,
                "fit_hist": fit_hist,
                "csv_path": csv_p,
            }
        )

    thrs: Optional[np.ndarray] = None
    if not using_csv_calls and crp_ref_keys:
        thrs = _assign_single_guide_from_h5_em(
            lane_records,
            crp_ref_keys,
            min_guide_thr=min_guide_thr,
            max_guide_thr=20,
            default_guide_thr=default_guide_thr,
        )

    all_kept_parts: List[Dict[str, Any]] = []
    n_safe_excl = n_unmapped_excl = n_multi_or_zero_excl = 0

    for rec in lane_records:
        qc_cols = rec["qc_cols"]
        barcodes = rec["barcodes"]
        keep_cols: List[int] = []
        nt_flags: List[bool] = []
        pkeys: List[str] = []
        kept_bcs: List[str] = []
        pmeta: Dict[str, Tuple[str, str, str, str]] = {}

        if using_csv_calls and rec["csv_path"] is not None:
            calls_map = _read_protospacer_calls_csv(rec["csv_path"])
            pairs = [(int(col_i), calls_map.get(barcodes[int(col_i)])) for col_i in qc_cols]
        else:
            C_qc = rec["C_qc"]
            crp_keys = rec["crp_keys"]
            t1, t2, a1 = extract_top2_csc(C_qc)
            assert thrs is not None
            pos_mask = C_qc.data >= thrs[C_qc.indices]
            pos_csc = sp.csc_matrix((pos_mask.astype(np.int8), C_qc.indices, C_qc.indptr), shape=C_qc.shape)
            n_pos = np.asarray(pos_csc.sum(axis=0)).ravel()
            s1 = (n_pos == 1) & (t1 >= float(min_guide_thr)) & ((t2 < 3.0) | (t1 / np.maximum(t2, 0.5) >= 3.0))
            pairs = [
                (int(qc_cols[loc_j]), crp_keys[int(a1[loc_j])] if bool(s1[loc_j]) else None)
                for loc_j in range(len(qc_cols))
            ]

        for raw_col, gcall in pairs:
            if not gcall:
                n_multi_or_zero_excl += 1
                continue
            if gcall in safe_keys:
                n_safe_excl += 1
                continue
            if gcall in ntc_keys:
                keep_cols.append(raw_col)
                nt_flags.append(True)
                pkeys.append("")
                kept_bcs.append(str(barcodes[raw_col]))
            elif gcall in guide_to_target:
                info = guide_to_target[gcall]
                keep_cols.append(raw_col)
                nt_flags.append(False)
                pkeys.append(info[1])
                pmeta.setdefault(info[1], info)
                kept_bcs.append(str(barcodes[raw_col]))
            else:
                n_unmapped_excl += 1

        if not keep_cols:
            continue

        order = np.argsort(np.asarray(keep_cols, dtype=np.int64), kind="mergesort")
        keep_arr = np.asarray(keep_cols, dtype=np.int64)[order]
        nt_arr = np.asarray(nt_flags, dtype=bool)[order]
        pkeys_ord = [pkeys[int(i)] for i in order]
        bcs_ord = [kept_bcs[int(i)] for i in order]

        all_kept_parts.append(
            {
                "rec": rec,
                "keep_arr": keep_arr,
                "nt_arr": nt_arr,
                "pkeys": pkeys_ord,
                "bcs": bcs_ord,
                "pmeta": pmeta,
            }
        )

    if not all_kept_parts:
        raise RuntimeError(f"EXCLUDED_FEASIBILITY: {line} yielded zero valid single-sgRNA cells.")

    # Concatenate global flags for Stage 1 / Stage 2 half + batch assignment across shards
    tot_kept = sum(len(x["keep_arr"]) for x in all_kept_parts)
    global_nt = np.concatenate([x["nt_arr"] for x in all_kept_parts])
    global_pkeys: List[str] = [pk for x in all_kept_parts for pk in x["pkeys"]]
    global_bcs: List[str] = [bc for x in all_kept_parts for bc in x["bcs"]]
    lane_nums = [x["rec"]["lane_num"] for x in all_kept_parts]
    distinct_lanes = sorted({x["rec"]["lane_tag"] for x in all_kept_parts if x["rec"]["lane_tag"] is not None})

    used_seed42_fallback = False
    half_rule_used = "channel_parity_L01L03_vs_L02L04"
    if (not force_seed42_split) and all(n is not None for n in lane_nums) and len(distinct_lanes) >= 2:
        global_half = np.concatenate(
            [
                np.full(len(x["keep_arr"]), 1 if int(x["rec"]["lane_num"]) % 2 == 1 else 2, dtype=np.int8)
                for x in all_kept_parts
            ]
        )
        raw_batches = np.concatenate(
            [np.full(len(x["keep_arr"]), str(x["rec"]["lane_tag"]), dtype=object) for x in all_kept_parts]
        )
        if np.sum(global_nt & (global_half == 1)) == 0 or np.sum(global_nt & (global_half == 2)) == 0:
            used_seed42_fallback = True
        else:
            global_batches, n_fb = _apply_min_ntc_batch_fallback(
                raw_batches, global_nt, global_half, min_batch_nt=min_batch_nt
            )
    else:
        used_seed42_fallback = True

    if used_seed42_fallback:
        half_rule_used = "seed42_stratified_50_50"
        rng42 = np.random.default_rng(42)
        global_half = np.ones(tot_kept, dtype=np.int8)
        grp_to_idx: Dict[str, List[int]] = defaultdict(list)
        for i, (is_nt, pk) in enumerate(zip(global_nt, global_pkeys)):
            grp_to_idx["__NTC__" if is_nt else pk].append(i)
        for g_key in sorted(grp_to_idx.keys()):
            idxs_sorted = sorted(grp_to_idx[g_key], key=lambda idx: (global_bcs[idx], idx))
            idxs = np.array(idxs_sorted, dtype=np.int64)
            perm = rng42.permutation(len(idxs))
            global_half[idxs[perm[len(idxs) // 2 :]]] = 2
        global_batches = np.where(global_half == 1, "half_1", "half_2").astype(str)
        n_fb = 0

    parts: List[PartSpec] = []
    offset = 0
    for pi, item in enumerate(all_kept_parts):
        rec = item["rec"]
        klen = len(item["keep_arr"])
        sl = slice(offset, offset + klen)
        offset += klen
        read_chunk, read_rows, read_chunk_cols = make_10x_csc_subset_readers(
            rec["path"], item["keep_arr"], rec["n_feat"], rec["rna_rows"]
        )
        tag = rec["lane_tag"] or f"part_{pi + 1}"
        parts.append(
            PartSpec(
                tag=f"{line}_{tag}",
                path=rec["path"],
                n_cells=klen,
                n_genes_raw=len(rec["rna_rows"]),
                nt_mask=item["nt_arr"],
                batches=global_batches[sl],
                half_flags=global_half[sl],
                umi=rec["umi_all"][item["keep_arr"]],
                pert_keys=item["pkeys"],
                pert_meta=item["pmeta"],
                var_gene_ids=rec["rids"],
                var_gene_names=rec["rnames"],
                read_chunk_fn=read_chunk,
                read_rows_fn=read_rows,
                read_chunk_cols_fn=read_chunk_cols,
            )
        )

    audit = {
        "line": line,
        "n_h5_files": len(h5_files),
        "guide_call_mode": "protospacer_calls_per_cell_csv" if using_csv_calls else "h5_crispr_umi_em",
        "guide_meta_source": guide_meta["source"],
        "has_safe_targeting": bool(guide_meta["has_safe_targeting"]),
        "ntc_rule_used": guide_meta["ntc_rule_used"],
        "min_umi": int(min_umi),
        "min_ngenes": int(min_ngenes),
        "max_pct_mt": (None if max_pct_mt is None else float(max_pct_mt)),
        "n_rna_qc_excluded": int(n_rna_qc_excl),
        "n_multi_or_zero_guide_excluded": int(n_multi_or_zero_excl),
        "n_safe_targeting_excluded": int(n_safe_excl),
        "n_unmapped_guide_excluded": int(n_unmapped_excl),
        "n_qc_kept_cells": int(tot_kept),
        "n_kept_ntc": int(global_nt.sum()),
        "n_kept_half_a_ntc": int(np.sum(global_nt & (global_half == 1))),
        "n_kept_half_b_ntc": int(np.sum(global_nt & (global_half == 2))),
        "n_low_ntc_channel_fallback_cells": int(n_fb),
        "used_seed42_fallback": bool(used_seed42_fallback),
        "half_rule_used": half_rule_used,
        "channels_used": sorted({str(b) for b in global_batches}),
        "median_guide_threshold": (None if thrs is None else float(np.median(thrs))),
    }
    return parts, audit


def prepare_r13_line(
    line: str,
    raw_dir: Path,
    out_root: Path,
    hgnc: Optional[HGNCLookup] = None,
    workers: int = 8,
    chunk_size: int = 8192,
) -> Path:
    """Standardize `GSE327644` (`mac_ctrl` or `mac_lps`) cleanly via frozen `process_r10_line`."""
    hgnc_l = hgnc or build_hgnc_lookup()
    parts, audit = build_gse327644_parts(raw_dir=raw_dir, hgnc=hgnc_l, line=line)
    return process_r10_line(
        line_name=line,
        parts=parts,
        hgnc=hgnc_l,
        chunk_size=chunk_size,
        out_root=out_root,
        extra_meta={"adapter_audit": audit},
        n_workers=workers,
    )


__all__ = [
    "build_gse327644_parts",
    "check_gse327644_feasibility",
    "load_gse327644_guide_metadata",
    "prepare_r13_line",
]
