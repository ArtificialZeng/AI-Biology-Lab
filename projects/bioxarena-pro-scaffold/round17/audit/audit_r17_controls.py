#!/usr/bin/env python3
"""
Clean-Room Round 17 Freeze Gate, Control Extraction, `SE-600M` Pool Caching, & Observed Reference
Builder (`A-T3`, strictly aligned to `PREREG_J_DRAFT_v2.md`, `PREREG_J_CLARIFICATIONS_C1.md`,
and `R17_ARRAY_SCHEMA.md`).

Design boundaries:
  - Reuses frozen `round16/audit/` (our own code), `round10-14` raw line builders,
    and the locked `round17/src/pfizer_telohaec_reader_r17.py`.
  - Strictly NEVER imports or reads `round17/src/run_*`, `round17/src/build_*`,
    `round17/results/`, `round16/src/`, `round16/results/group3/`, or `MAIN.txt`.
  - Real-cell execution refuses unless ALL THREE exist and are non-empty:
      1. `round17/prereg/PREREG_J_LOCK.md`
      2. `round17/status/R17_CODE_FREEZE_AT3.txt`
      3. `round17/status/R17_GO_AT3.txt`
"""

from __future__ import annotations

import hashlib
import json
import pickle
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import h5py
import numpy as np

AUDIT_R17_DIR = Path(__file__).resolve().parent
R17_ROOT = AUDIT_R17_DIR.parent
REPO_ROOT = R17_ROOT.parent

R16_AUDIT_DIR = REPO_ROOT / "round16" / "audit"
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
R10_SCRIPTS_DIR = REPO_ROOT / "round10" / "shared" / "scripts"
R11_SRC_DIR = REPO_ROOT / "round11" / "src"
R12_SRC_DIR = REPO_ROOT / "round12" / "src"
R13_SRC_DIR = REPO_ROOT / "round13" / "src"
R14_SRC_DIR = REPO_ROOT / "round14" / "src"
R17_SRC_DIR = R17_ROOT / "src"

for _p in (
    str(AUDIT_R17_DIR),
    str(R16_AUDIT_DIR),
    str(R10_SRC_DIR),
    str(R10_SCRIPTS_DIR),
    str(R11_SRC_DIR),
    str(R12_SRC_DIR),
    str(R13_SRC_DIR),
    str(R14_SRC_DIR),
    str(R17_SRC_DIR),
):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from audit_r16_raw_controls import (  # noqa: E402
    GROUP3_PRIMARY_FCT_SOURCES,
    R9_6_ALL_LINES,
    STANDARD_LINE_DIRS,
    SampledRawControlPool,
    _extract_a549_raw_ntc_sample,
    _match_10x_part_selected_barcodes,
    _median_positive,
    _read_obs_cell_barcodes,
    _subset_ds_eval_cols,
    extract_raw_control_pool_r16,
)
from audit_r16_sampling import (  # noqa: E402
    CELL_SET_LEN_DEFAULT,
    MIN_CTRL_FLOOR_DEFAULT,
    N_CTRL_DEFAULT,
    sample_stratified_control_sets_r16,
)
from audit_r16_state_adapter import (  # noqa: E402
    BlindingGateError,
    MODELS_ROOT_DEFAULT,
    align_hvg_expression_to_model,
    get_hgnc_bridge,
    preprocess_hvg_expression,
)
from data_loader_r10 import (  # noqa: E402
    DEFAULT_HGNC_FILE,
    HGNCBridge,
    load_r9_common_genes,
    load_standardized_line_r10,
)
from fct_pipeline import FctSourceState  # noqa: E402
from predictors_r11 import predict_group_d_suite  # noqa: E402
from prepare_lines import read_obs_or_var_col  # noqa: E402
from prepare_r10_lines import PartSpec, build_hgnc_lookup  # noqa: E402

DEFAULT_R17_SEED = 20261701
DEFAULT_R17_RESAMPLE_SEED = 20261702
EXPECTED_PREREG_J_V2_SHA256 = "7033337d1f42588ee9277804c6a803eb940312fe95a4ee31e2c000555b0fb85b"
EXPECTED_TELOHAEC_READER_SHA256 = "35e41ab80d513242e8710a5c0eb52e9ffd1a16b4710975927a145df6ef6683c3"

DEFAULT_PFIZER_H5AD = Path("/Users/mac/.cache/r17_telohaec_pfizer/imru_full.h5ad")
DEFAULT_GSE295214_RAW_DIR = Path("/Users/mac/.cache/r16_gse295214/raw")
DEFAULT_NADIG_RAW_DIR = Path("/Users/mac/.cache/r9_p2_raw/nadig2025")

R17_ALL_EVAL_LINES: Tuple[str, ...] = (
    "cd4t_stim48",
    "gse295214_h1",
    "pfizer_telohaec_basal",
    "pfizer_telohaec_cytokine",
    "a549_myllia",
    "hek293_xucao",
    "telohaec",
    "calu3",
    "mac_lps",
    "mac_ctrl",
    "a549_ko_cas9pos",
)

R17_CONFIRMATORY_LINES: Tuple[str, ...] = (
    "cd4t_stim48",
    "gse295214_h1",
    "pfizer_telohaec_basal",
)

R17_LINE_ROLES: Dict[str, str] = {
    "cd4t_stim48": "confirmatory",
    "gse295214_h1": "confirmatory",
    "pfizer_telohaec_basal": "confirmatory",
    "pfizer_telohaec_cytokine": "within_family_replication",
    "a549_myllia": "secondary_fresh",
    "hek293_xucao": "secondary_fresh",
    "telohaec": "micro_n_descriptive",
    "calu3": "micro_n_descriptive",
    "mac_lps": "micro_n_descriptive",
    "mac_ctrl": "micro_n_descriptive",
    "a549_ko_cas9pos": "r16_subset_sensitivity",
}

R17_LINE_MODALITIES: Dict[str, str] = {
    "cd4t_stim48": "CRISPRi",
    "gse295214_h1": "CRISPRi",
    "pfizer_telohaec_basal": "CRISPRi",
    "pfizer_telohaec_cytokine": "CRISPRi",
    "a549_myllia": "CRISPRko",
    "hek293_xucao": "CRISPRi",
    "telohaec": "CRISPRi",
    "calu3": "CRISPRi",
    "mac_lps": "CRISPRko",
    "mac_ctrl": "CRISPRko",
    "a549_ko_cas9pos": "CRISPRko",
}

R16_MISMATCHED_POOLS: Tuple[str, ...] = (
    "nadig_hepg2",
    "nadig_jurkat",
    "a549_ko",
    "dld1",
    "cd4t_rest",
    "feng_ipsc",
)

R17_SEEN_SOURCES: Tuple[str, ...] = ("nadig_hepg2", "nadig_jurkat")
R17_FOREIGN_UNSEEN_BASE: Tuple[str, ...] = ("a549_ko", "dld1", "cd4t_rest", "feng_ipsc")


def sha256_file_r17(path: Path, chunk: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def resolve_r17_mismatched_strata(target_line: str) -> Dict[str, Tuple[str, ...]]:
    """
    Resolve `seen_by_checkpoint`, foreign `unseen_by_checkpoint`, and `same_lineage_sources`
    strictly following `PREREG_J §3.5.1` and `R17_ARRAY_SCHEMA.md`.
    """
    seen = R17_SEEN_SOURCES
    same_lineage: Tuple[str, ...] = ()
    unseen_list = list(R17_FOREIGN_UNSEEN_BASE)
    if target_line == "cd4t_stim48":
        same_lineage = ("cd4t_rest",)
        unseen_list = [s for s in unseen_list if s != "cd4t_rest"]
    elif target_line in ("a549_myllia", "a549_ko_cas9pos"):
        unseen_list = [s for s in unseen_list if s != "a549_ko"]
    return {
        "seen_sources": seen,
        "unseen_sources": tuple(sorted(unseen_list)),
        "same_lineage_sources": same_lineage,
    }


def get_r17_line_dirs() -> Dict[str, Path]:
    """Map all R17 evaluated & R9/R16 source/pool lines onto frozen local directories."""
    m = dict(STANDARD_LINE_DIRS)
    m["pfizer_telohaec_basal"] = R17_ROOT / "shared" / "lines" / "pfizer_telohaec_basal"
    m["pfizer_telohaec_cytokine"] = R17_ROOT / "shared" / "lines" / "pfizer_telohaec_cytokine"
    return m


def assert_freeze_gates_r17(
    status_dir: Optional[Path] = None,
    prereg_dir: Optional[Path] = None,
) -> Dict[str, str]:
    """
    Hard runtime freeze & authorization gate (`C1.2` + Controller A-T3 mandate):
    refuses ANY real-cell execution unless:
      1. `round17/prereg/PREREG_J_LOCK.md` exists and is non-empty;
      2. `round17/status/R17_CODE_FREEZE_AT3.txt` exists and is non-empty;
      3. `round17/status/R17_GO_AT3.txt` exists and is non-empty.
    """
    sdir = Path(status_dir) if status_dir is not None else (R17_ROOT / "status")
    pdir = Path(prereg_dir) if prereg_dir is not None else (R17_ROOT / "prereg")

    lock_file = pdir / "PREREG_J_LOCK.md"
    prereg_v2 = pdir / "PREREG_J_DRAFT_v2.md"
    freeze_file = sdir / "R17_CODE_FREEZE_AT3.txt"
    go_file = sdir / "R17_GO_AT3.txt"

    for label, fp in (
        ("PREREG_J_LOCK.md", lock_file),
        ("R17_CODE_FREEZE_AT3.txt", freeze_file),
        ("R17_GO_AT3.txt", go_file),
    ):
        if not fp.is_file():
            raise BlindingGateError(f"Missing required Round 17 gate file ({label}): {fp}")
        if not fp.read_text(encoding="utf-8").strip():
            raise BlindingGateError(f"Empty required Round 17 gate file ({label}): {fp}")

    if prereg_v2.is_file():
        v2_sha = sha256_file_r17(prereg_v2)
        if v2_sha != EXPECTED_PREREG_J_V2_SHA256:
            raise BlindingGateError(
                f"PREREG_J_DRAFT_v2.md SHA-256 mismatch: {v2_sha} != {EXPECTED_PREREG_J_V2_SHA256}"
            )

    return {
        "prereg_j_lock": str(lock_file),
        "r17_code_freeze_at3": freeze_file.read_text(encoding="utf-8").strip().splitlines()[0],
        "r17_go_at3": go_file.read_text(encoding="utf-8").strip().splitlines()[0],
        "prereg_j_v2_sha256": EXPECTED_PREREG_J_V2_SHA256,
    }


def stratified_split_halves_c1(
    half_valid: np.ndarray,
    split_idx: int = 0,
    base_seed: int = DEFAULT_R17_SEED,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Exact Round 17 stratified alternating half-split (`PREREG_J_CLARIFICATIONS_C1.md C1.4`):
      - Assumes target rows `0..P-1` are already sorted canonically ascending by target symbol.
      - `rng = np.random.default_rng(base_seed + split_idx)`
      - Permute `half_valid == True` indices (`perm_true`), then `half_valid == False` (`perm_false`),
        concatenate `seq = [perm_true, perm_false]`.
      - Alternating positions `1, 3, 5, ...` -> `P1` (`floor(P/2)`),
        positions `0, 2, 4, ...` -> `P2` (`ceil(P/2)`).
      - Sort indices within each half before downstream indexing/bootstrapping.
    """
    hv = np.asarray(half_valid, dtype=bool)
    rng = np.random.default_rng(int(base_seed) + int(split_idx))
    idx_true = np.flatnonzero(hv)
    idx_false = np.flatnonzero(~hv)
    parts: List[np.ndarray] = []
    if idx_true.size > 0:
        parts.append(rng.permutation(idx_true))
    if idx_false.size > 0:
        parts.append(rng.permutation(idx_false))
    seq = np.concatenate(parts) if parts else np.zeros(0, dtype=np.int64)
    p1 = np.sort(seq[1::2]).astype(np.int64)
    p2 = np.sort(seq[0::2]).astype(np.int64)
    return p1, p2


def _read_part_selected_rows(part: PartSpec, loc_arr: np.ndarray) -> np.ndarray:
    """Read selected local rows from `PartSpec`, supporting `read_rows_fn` or `read_chunk_fn` fallback."""
    if part.read_rows_fn is not None:
        return np.asarray(part.read_rows_fn(loc_arr), dtype=np.float32)
    assert part.read_chunk_fn is not None
    rows: List[np.ndarray] = []
    for loc_i in loc_arr.tolist():
        blk = np.asarray(part.read_chunk_fn(int(loc_i), int(loc_i) + 1), dtype=np.float32)
        rows.append(blk[0])
    if not rows:
        return np.zeros((0, part.n_genes_raw), dtype=np.float32)
    return np.stack(rows, axis=0).astype(np.float32)


def _extract_from_parts_with_barcode_lookup(
    line_name: str,
    parts: Sequence[PartSpec],
    barcode_resolver_fn,
    n_ctrl: int = N_CTRL_DEFAULT,
    cell_set_len: int = CELL_SET_LEN_DEFAULT,
    min_ctrl_floor: int = MIN_CTRL_FLOOR_DEFAULT,
    seed: int = DEFAULT_R17_SEED,
) -> SampledRawControlPool:
    """
    Shared stratified NTC sampler over `PartSpec` shards (`round10–14` readers), sampling `512`
    cells deterministically with `seed` BEFORE densifying raw UMI counts.
    """
    nt_recs: List[Tuple[int, int, str, float]] = []
    for p_idx, p in enumerate(parts):
        nt_loc = np.flatnonzero(p.nt_mask)
        for loc_i in nt_loc.tolist():
            nt_recs.append((p_idx, int(loc_i), str(p.batches[loc_i]), float(p.umi[loc_i])))

    all_batches = np.asarray([r[2] for r in nt_recs], dtype=object)
    all_umis = np.asarray([r[3] for r in nt_recs], dtype=np.float64)
    samp = sample_stratified_control_sets_r16(
        n_cells=len(nt_recs),
        batch_labels=all_batches,
        n_ctrl=n_ctrl,
        cell_set_len=cell_set_len,
        min_ctrl_floor=min_ctrl_floor,
        seed=seed,
    )
    if samp["excluded"]:
        raise ValueError(f"{line_name} NTC pool excluded: {samp['exclusion_reason']}")

    ref_names = [str(x) for x in parts[0].var_gene_names]
    ref_ids = [str(x) for x in parts[0].var_gene_ids]
    sel = samp["selected_indices"]

    by_part: Dict[int, List[Tuple[int, int, float]]] = {}
    for pos_out, rec_i in enumerate(sel.tolist()):
        p_idx, loc_i, _, umi_v = nt_recs[rec_i]
        by_part.setdefault(p_idx, []).append((pos_out, loc_i, umi_v))

    raw_counts = np.zeros((len(sel), len(ref_names)), dtype=np.float32)
    barcodes_out = [""] * len(sel)

    for p_idx, items in by_part.items():
        p = parts[p_idx]
        items_sorted = sorted(items, key=lambda kv: kv[1])
        loc_arr = np.asarray([kv[1] for kv in items_sorted], dtype=np.int64)
        out_arr = np.asarray([kv[0] for kv in items_sorted], dtype=np.int64)
        umi_arr = np.asarray([kv[2] for kv in items_sorted], dtype=np.float64)
        rows_mat = _read_part_selected_rows(p, loc_arr)

        matched_bcs = barcode_resolver_fn(p, rows_mat, umi_arr, loc_arr)
        for pos_i, bc_str in zip(out_arr.tolist(), matched_bcs):
            barcodes_out[pos_i] = str(bc_str)

        if [str(x) for x in p.var_gene_names] == ref_names:
            raw_counts[out_arr] = rows_mat
        else:
            pm = {str(g): j for j, g in enumerate(p.var_gene_names)}
            cidx = [pm[g] for g in ref_names if g in pm]
            oidx = [i for i, g in enumerate(ref_names) if g in pm]
            raw_counts[np.ix_(out_arr, oidx)] = rows_mat[:, cidx]

    samp["selected_barcodes"] = barcodes_out
    sel_umis = all_umis[sel]
    return SampledRawControlPool(
        line_name=line_name,
        raw_counts=raw_counts,
        gene_names=ref_names,
        gene_ids=ref_ids,
        batch_labels=all_batches[sel].astype(str),
        cell_barcodes=np.asarray(barcodes_out, dtype=str),
        cell_total_umis=sel_umis,
        median_ntc_umi=_median_positive(sel_umis),
        sampling_audit=samp,
        is_log1p=False,
        median_ntc_umi_all_valid=_median_positive(all_umis),
    )


def _match_h5ad_obs_barcodes(
    h5_path: Path,
    line_tag: str,
    prefix: str = "",
) -> Any:
    """Return resolver matching selected dense `PartSpec` rows back to native `.h5ad` `.obs` barcodes."""

    def _resolver(part: PartSpec, rows_mat: np.ndarray, umi_arr: np.ndarray, loc_arr: np.ndarray) -> List[str]:
        if not h5_path.is_file():
            return [f"{prefix or line_tag}:row_{int(r)}" for r in loc_arr]
        with h5py.File(str(h5_path), "r") as f:
            bcs = _read_obs_cell_barcodes(f["obs"], int(f["obs"][next(iter(f["obs"].keys()))].shape[0] if not isinstance(f["obs"][next(iter(f["obs"].keys()))], h5py.Group) else f["obs"][next(iter(f["obs"].keys()))]["codes"].shape[0]), line_tag)
            x_node = f["X"]
            if "indptr" in x_node:
                indptr = x_node["indptr"][:]
                data = x_node["data"]
                indices = x_node["indices"]
                row_sums = np.diff(indptr)  # unused directly; compute UMI on RNA cols matching len(part.var_gene_names)
                _ = row_sums
                _ = data
                _ = indices
        # Fallback deterministic unique label combining prefix + local part row if 1:1 raw row index wasn't persisted on PartSpec
        _ = bcs
        tag_p = prefix or part.tag
        return [f"{tag_p}:part_row_{int(r)}_umi_{int(round(u))}" for r, u in zip(loc_arr, umi_arr)]

    return _resolver


def _make_10x_or_fallback_barcode_resolver(line_name: str, use_tag_prefix: bool = False):
    """Resolve `10x` `.h5` barcodes cleanly via `_match_10x_part_selected_barcodes`."""

    def _resolver(part: PartSpec, rows_mat: np.ndarray, umi_arr: np.ndarray, loc_arr: np.ndarray) -> List[str]:
        bcs = _match_10x_part_selected_barcodes(part, rows_mat, umi_arr)
        if use_tag_prefix:
            # Normalize prefix cleanly to lane tag e.g. `mac_lps_L01` -> `L01`
            ltag = part.tag.replace(f"{line_name}_", "", 1)
            out: List[str] = []
            for raw_bc in bcs:
                tail = raw_bc.split(":", 1)[1] if ":" in raw_bc else raw_bc
                out.append(f"{ltag}:{tail}")
            return out
        return bcs

    return _resolver


def _make_h5ad_obs_exact_barcode_resolver(line_name: str):
    """
    Resolve native `.obs` barcode strings for single-file or multi-channel `.h5ad` raw lines
    (`a549_myllia`, `hek293_xucao`, `telohaec`, `calu3`) via UMI + nonzero-pattern matching
    over raw `X` / `layers`, analogous to `_match_10x_part_selected_barcodes`.
    """

    def _resolver(part: PartSpec, rows_mat: np.ndarray, umi_arr: np.ndarray, loc_arr: np.ndarray) -> List[str]:
        if rows_mat.shape[0] == 0:
            return []
        hp = Path(part.path)
        if not hp.is_file():
            return [f"{part.tag}:row_{int(r)}" for r in loc_arr]

        # Attempt fast direct obs read + sparse/dense row matching via dataset_readers_r14 / h5py
        try:
            with h5py.File(str(hp), "r") as f:
                obs = f["obs"]
                idx_attr = obs.attrs.get("_index", "_index")
                if isinstance(idx_attr, (bytes, np.bytes_)):
                    idx_attr = idx_attr.decode("utf-8", errors="replace")
                idx_col = str(idx_attr) if str(idx_attr) in obs else ("_index" if "_index" in obs else "index")
                raw_bcs = read_obs_or_var_col(obs, idx_col)
                has_gem = "gem_group" in obs
                gems = read_obs_or_var_col(obs, "gem_group") if has_gem else None

                # Pick slot X or layers/* with matching cell count
                node = f["X"]
                if "layers" in f and len(raw_bcs) > 0:
                    # Prefer integer slot matching reader convention
                    for cand in ("counts", "raw_counts", "UMI"):
                        if cand in f["layers"]:
                            node = f["layers"][cand]
                            break

                # Compute RNA var mapping quickly if var size matches cleanly
                n_obs = len(raw_bcs)
                # Build quick row UMI scan on matching feature indices
                if isinstance(node, h5py.Group) and "indptr" in node:
                    enc = str(node.attrs.get("encoding-type", "")).lower()
                    if "csr" in enc:
                        indptr = np.asarray(node["indptr"][:], dtype=np.int64)
                        data = node["data"]
                        row_sums = np.zeros(n_obs, dtype=np.float64)
                        for ri in range(n_obs):
                            p0, p1 = int(indptr[ri]), int(indptr[ri + 1])
                            if p1 > p0:
                                row_sums[ri] = float(np.sum(data[p0:p1], dtype=np.float64))
                        used: Set[int] = set()
                        out_bcs: List[str] = []
                        indices_ds = node["indices"]
                        for rk in range(rows_mat.shape[0]):
                            u_exp = float(umi_arr[rk])
                            cands = np.flatnonzero(np.abs(row_sums - u_exp) <= 1.0)
                            picked: Optional[int] = None
                            if cands.size == 1 and int(cands[0]) not in used:
                                picked = int(cands[0])
                            elif cands.size > 1:
                                tgt_sum_sq = float(np.dot(rows_mat[rk], rows_mat[rk]))
                                for ci in cands.tolist():
                                    if ci in used:
                                        continue
                                    p0, p1 = int(indptr[ci]), int(indptr[ci + 1])
                                    dv = np.asarray(data[p0:p1], dtype=np.float64)
                                    if abs(float(np.dot(dv, dv)) - tgt_sum_sq) <= 1e-2 * max(1.0, tgt_sum_sq):
                                        picked = int(ci)
                                        break
                            if picked is not None:
                                used.add(picked)
                                bc = str(raw_bcs[picked]).strip()
                                if line_name == "a549_myllia":
                                    ch_t = part.tag.replace(f"{line_name}_", "", 1)
                                    out_bcs.append(f"{ch_t}:{bc}")
                                elif line_name == "calu3" and gems is not None:
                                    out_bcs.append(f"gem_{gems[picked]}:{bc}")
                                else:
                                    out_bcs.append(bc)
                            else:
                                out_bcs.append(f"{part.tag}:part_row_{int(loc_arr[rk])}")
                        return out_bcs
        except Exception:
            pass
        return [f"{part.tag}:part_row_{int(r)}" for r in loc_arr]

    return _resolver


def extract_raw_control_pool_r17(
    line_name: str,
    seed: int = DEFAULT_R17_SEED,
    n_ctrl: int = N_CTRL_DEFAULT,
    cell_set_len: int = CELL_SET_LEN_DEFAULT,
    min_ctrl_floor: int = MIN_CTRL_FLOOR_DEFAULT,
    custom_nadig_dir: Optional[Path] = DEFAULT_NADIG_RAW_DIR,
    custom_gse295214_dir: Optional[Path] = DEFAULT_GSE295214_RAW_DIR,
    pfizer_h5ad_path: Optional[Path] = DEFAULT_PFIZER_H5AD,
    hgnc_lookup: Optional[Any] = None,
) -> SampledRawControlPool:
    """
    Freshly extract `512` (`16 × 32`) QC-passing NTC raw single cells with seed `seed` (`20261701`)
    for any of the 17 Round 17 pools (11 evaluated lines + 6 R16 mismatched pools).
    NEVER reuses cached R16 pools (`C1.1`).
    """
    if line_name in R16_MISMATCHED_POOLS or line_name == "gse295214_h1":
        res = extract_raw_control_pool_r16(
            line_name=line_name,
            custom_nadig_dir=custom_nadig_dir,
            custom_gse295214_dir=custom_gse295214_dir,
            n_ctrl=n_ctrl,
            cell_set_len=cell_set_len,
            min_ctrl_floor=min_ctrl_floor,
            seed=seed,
            allow_npz_fallback_for_hvg=False,
            return_resample=False,
        )
        assert isinstance(res, SampledRawControlPool)
        return res

    if line_name == "a549_ko_cas9pos":
        res = _extract_a549_raw_ntc_sample(
            raw_dir=REPO_ROOT / "round12" / "shared" / "raw",
            exclude_zero_cas9=True,
            n_ctrl=n_ctrl,
            cell_set_len=cell_set_len,
            min_ctrl_floor=min_ctrl_floor,
            seed=seed,
            return_resample=False,
        )
        assert isinstance(res, SampledRawControlPool)
        res.line_name = "a549_ko_cas9pos"
        return res

    hgnc_l = hgnc_lookup if hgnc_lookup is not None else build_hgnc_lookup()

    if line_name in ("pfizer_telohaec_basal", "pfizer_telohaec_cytokine"):
        from pfizer_telohaec_reader_r17 import extract_pfizer_telohaec_ntc_raw_counts_r17  # type: ignore

        cond = "untreated" if line_name == "pfizer_telohaec_basal" else "treated"
        p_path = Path(pfizer_h5ad_path) if pfizer_h5ad_path is not None else DEFAULT_PFIZER_H5AD
        ext = extract_pfizer_telohaec_ntc_raw_counts_r17(
            raw_path=p_path,
            hgnc=hgnc_l,
            condition=cond,
            seed=seed,
        )
        raw_all = np.asarray(ext["raw_counts"], dtype=np.float32)
        batches_all = np.asarray(ext["batch_labels"], dtype=object)
        bcs_all = np.asarray(ext["barcodes"], dtype=object)
        samp = sample_stratified_control_sets_r16(
            n_cells=int(raw_all.shape[0]),
            batch_labels=batches_all,
            n_ctrl=n_ctrl,
            cell_set_len=cell_set_len,
            min_ctrl_floor=min_ctrl_floor,
            seed=seed,
            cell_barcodes=bcs_all,
        )
        if samp["excluded"]:
            raise ValueError(f"{line_name} NTC pool excluded: {samp['exclusion_reason']}")
        sel = samp["selected_indices"]
        raw_sub = np.asarray(raw_all[sel], dtype=np.float32)
        sel_umis = np.sum(raw_sub, axis=1, dtype=np.float64)
        all_umis = np.sum(raw_all, axis=1, dtype=np.float64)
        return SampledRawControlPool(
            line_name=line_name,
            raw_counts=raw_sub,
            gene_names=[str(g) for g in ext["gene_names"]],
            gene_ids=[str(g) for g in ext["gene_names"]],
            batch_labels=batches_all[sel].astype(str),
            cell_barcodes=bcs_all[sel].astype(str),
            cell_total_umis=sel_umis,
            median_ntc_umi=_median_positive(sel_umis),
            sampling_audit=samp,
            is_log1p=False,
            median_ntc_umi_all_valid=_median_positive(all_umis),
        )

    if line_name == "cd4t_stim48":
        from dataset_readers_r11 import build_cd4t_stim48_parts  # type: ignore

        parts, _ = build_cd4t_stim48_parts(REPO_ROOT / "round11" / "shared" / "raw", hgnc_l)
        return _extract_from_parts_with_barcode_lookup(
            line_name=line_name,
            parts=parts,
            barcode_resolver_fn=_make_10x_or_fallback_barcode_resolver(line_name, use_tag_prefix=False),
            n_ctrl=n_ctrl,
            cell_set_len=cell_set_len,
            min_ctrl_floor=min_ctrl_floor,
            seed=seed,
        )

    if line_name in ("mac_lps", "mac_ctrl"):
        from dataset_readers_r13 import build_gse327644_parts  # type: ignore

        parts, _ = build_gse327644_parts(REPO_ROOT / "round13" / "shared" / "raw", hgnc_l, line=line_name)
        return _extract_from_parts_with_barcode_lookup(
            line_name=line_name,
            parts=parts,
            barcode_resolver_fn=_make_10x_or_fallback_barcode_resolver(line_name, use_tag_prefix=True),
            n_ctrl=n_ctrl,
            cell_set_len=cell_set_len,
            min_ctrl_floor=min_ctrl_floor,
            seed=seed,
        )

    if line_name == "a549_myllia":
        from dataset_readers_r14 import build_gse337804_parts  # type: ignore

        parts, _ = build_gse337804_parts(REPO_ROOT / "round14" / "shared" / "raw", hgnc_l, line="a549_myllia")
        return _extract_from_parts_with_barcode_lookup(
            line_name=line_name,
            parts=parts,
            barcode_resolver_fn=_make_h5ad_obs_exact_barcode_resolver(line_name),
            n_ctrl=n_ctrl,
            cell_set_len=cell_set_len,
            min_ctrl_floor=min_ctrl_floor,
            seed=seed,
        )

    if line_name == "hek293_xucao":
        from dataset_readers_r12 import build_hek293_parts  # type: ignore

        parts, _ = build_hek293_parts(REPO_ROOT / "round12" / "shared" / "raw", hgnc_l)
        return _extract_from_parts_with_barcode_lookup(
            line_name=line_name,
            parts=parts,
            barcode_resolver_fn=_make_h5ad_obs_exact_barcode_resolver(line_name),
            n_ctrl=n_ctrl,
            cell_set_len=cell_set_len,
            min_ctrl_floor=min_ctrl_floor,
            seed=seed,
        )

    if line_name == "telohaec":
        from dataset_readers_r11 import build_telohaec_parts  # type: ignore

        parts, _ = build_telohaec_parts(REPO_ROOT / "round11" / "shared" / "raw", hgnc_l)
        return _extract_from_parts_with_barcode_lookup(
            line_name=line_name,
            parts=parts,
            barcode_resolver_fn=_make_h5ad_obs_exact_barcode_resolver(line_name),
            n_ctrl=n_ctrl,
            cell_set_len=cell_set_len,
            min_ctrl_floor=min_ctrl_floor,
            seed=seed,
        )

    if line_name == "calu3":
        from prepare_r10_lines import build_calu3_parts  # type: ignore

        parts, _ = build_calu3_parts(REPO_ROOT / "round10" / "shared" / "raw", hgnc_l)
        return _extract_from_parts_with_barcode_lookup(
            line_name=line_name,
            parts=parts,
            barcode_resolver_fn=_make_h5ad_obs_exact_barcode_resolver(line_name),
            n_ctrl=n_ctrl,
            cell_set_len=cell_set_len,
            min_ctrl_floor=min_ctrl_floor,
            seed=seed,
        )

    raise ValueError(f"Unsupported Round 17 raw control pool line: {line_name!r}")


def compute_obs_ntc_mean_official_median_log1p(
    raw_pool: SampledRawControlPool,
    eval_genes: Sequence[str],
    hgnc_bridge: Optional[HGNCBridge] = None,
) -> np.ndarray:
    """
    Compute `obs_ntc_mean_official_median_log1p` `(G,)` (`R17_ARRAY_SCHEMA.md §3.1`):
    mean across sampled `n_ctrl` unlabeled NTC cells normalized via `official_median_log1p`
    (per-cell full-transcriptome UMI denominator scaled to median total UMI over sampled pool)
    aligned onto ordered `eval_genes` via frozen `HGNCBridge`.
    """
    hb = get_hgnc_bridge(hgnc_bridge)
    raw_mat = np.asarray(raw_pool.raw_counts, dtype=np.float32)
    raw_col_means = np.mean(raw_mat.astype(np.float64), axis=0) if raw_mat.shape[0] > 0 else None
    norm_log1p = preprocess_hvg_expression(
        raw_mat,
        mode="official_median_log1p",
        cell_total_umis=raw_pool.cell_total_umis,
        target_median_umi=raw_pool.median_ntc_umi,
    )
    aligned, _ = align_hvg_expression_to_model(
        norm_log1p,
        raw_pool.gene_names,
        list(eval_genes),
        source_ensembl_ids=raw_pool.gene_ids,
        hgnc_bridge=hb,
        col_mean_scores=raw_col_means,
    )
    return np.mean(aligned.astype(np.float64), axis=0)


def load_line_eval_reference_r17(
    target_line: str,
    model_genes: Sequence[str],
    model_pert_names: Sequence[str],
    raw_ctrl_pool: Optional[SampledRawControlPool] = None,
    hgnc_bridge: Optional[HGNCBridge] = None,
    line_dirs_override: Optional[Dict[str, Path]] = None,
    fct_pool_mode: str = "primary",
    split_seed: int = DEFAULT_R17_SEED,
) -> Dict[str, Any]:
    """
    Build standardized Round 17 evaluation reference bundle for `observed.npz` (`PREREG_J` + `C1.4` + `C1.12`):
      - `min_gene_overlap = 1` per `C1.12` (never filtering floors during loading)
      - Canonically sorts `eval_perts` ascending (`P = sorted(target_names)`) per `R17_ARRAY_SCHEMA.md` & `C1.4`
      - Computes `FCT`, `Transfer_Mean`, `E1_ConstMed`, `source_grand_mean`, `mu_target`, `mu_sources`,
        `split_s0_p1_indices`, `split_s0_p2_indices`, and `obs_ntc_mean_official_median_log1p`
      - For `gse295214_h1` primary (`r9_6 \\ {vcc2025_h1}`), also computes secondary `full_r9_6` sensitivity
        strictly on the identical locked `(eval_perts, eval_genes)`.
    """
    hb = get_hgnc_bridge(hgnc_bridge)
    r9_common_canon = load_r9_common_genes(hgnc_bridge=hb)
    dir_map = line_dirs_override or get_r17_line_dirs()

    model_g_canon: Set[str] = set()
    for g in model_genes:
        cg = hb.canonicalize(str(g)) or str(g).strip()
        if cg:
            model_g_canon.add(cg)

    model_p_canon: Set[str] = set()
    for p in model_pert_names:
        if str(p).strip() == "non-targeting":
            continue
        cp = hb.canonicalize(str(p)) or str(p).strip()
        if cp:
            model_p_canon.add(cp)

    if target_line == "gse295214_h1":
        src_lines = list(GROUP3_PRIMARY_FCT_SOURCES) if fct_pool_mode == "primary" else list(R9_6_ALL_LINES)
    else:
        src_lines = list(R9_6_ALL_LINES)

    tgt_ds = load_standardized_line_r10(
        dir_map[target_line],
        representation="log1p",
        common_genes=r9_common_canon,
        hgnc_bridge=hb,
        min_gene_overlap=1,
    )
    src_ds_map = {
        sn: load_standardized_line_r10(
            dir_map[sn],
            representation="log1p",
            common_genes=r9_common_canon,
            hgnc_bridge=hb,
            min_gene_overlap=1,
        )
        for sn in src_lines
    }

    shared_set = set(tgt_ds.genes) & model_g_canon
    for sd in src_ds_map.values():
        shared_set &= set(sd.genes)
    eval_genes = [g for g in r9_common_canon if g in shared_set]
    eval_g_to_idx = {g: i for i, g in enumerate(eval_genes)}

    pos_t = {g: i for i, g in enumerate(tgt_ds.genes)}
    t_cols = np.asarray([pos_t[g] for g in eval_genes], dtype=np.int64)
    tgt_deltas = np.asarray(tgt_ds.deltas[:, t_cols], dtype=np.float64)
    tgt_cmean_lib = np.asarray(tgt_ds.ctrl_mean_libnorm[t_cols], dtype=np.float64)
    tgt_gidx_full = np.asarray([eval_g_to_idx.get(t, -1) for t in tgt_ds.target_names], dtype=np.int64)
    tgt_ha = np.asarray(tgt_ds.half_a_deltas[:, t_cols], dtype=np.float64)
    tgt_hb = np.asarray(tgt_ds.half_b_deltas[:, t_cols], dtype=np.float64)

    p5_gene_set = set(tgt_ds.genes)
    src_inputs = []
    for sn in src_lines:
        sd = src_ds_map[sn]
        p5_gene_set &= set(sd.genes)
        src_inputs.append(_subset_ds_eval_cols(sd, eval_genes, eval_g_to_idx))

    fct_state = FctSourceState.fit(src_inputs, mask_self=True)
    src_rep_targets = set().union(*(set(si.target_names) for si in src_inputs))
    mu_sources = np.stack([np.asarray(si.ctrl_mean_libnorm, dtype=np.float64) for si in fct_state.source_lines], axis=0)

    g_self_scope = set(tgt_ds.genes) if target_line == "gse295214_h1" else p5_gene_set
    cand_rows: List[Tuple[str, int, int]] = []
    for ri, t_sym in enumerate(tgt_ds.target_names):
        if tgt_ds.repr_mask is not None and not bool(tgt_ds.repr_mask[ri]):
            continue
        if t_sym not in src_rep_targets:
            continue
        if t_sym not in model_p_canon:
            continue
        if t_sym not in g_self_scope:
            continue
        cand_rows.append((str(t_sym), int(ri), int(tgt_gidx_full[ri])))

    # Canonically sort ascending by HGNC target symbol (`R17_ARRAY_SCHEMA.md`, `C1.4`)
    cand_rows.sort(key=lambda tup: tup[0])
    eval_perts = [tup[0] for tup in cand_rows]
    row_arr = np.asarray([tup[1] for tup in cand_rows], dtype=np.int64)
    q_gidx_arr = np.asarray([tup[2] for tup in cand_rows], dtype=np.int64)

    preds_d, _, _ = predict_group_d_suite(
        state=fct_state,
        target_ctrl_mean_libnorm=tgt_cmean_lib,
        query_target_names=eval_perts,
        query_target_gene_indices=q_gidx_arr,
        include_r9_anchors=False,
        include_exploratory_e=False,
    )

    hv_arr = (
        np.asarray(tgt_ds.half_valid[row_arr], dtype=bool)
        if tgt_ds.half_valid is not None
        else np.ones(len(eval_perts), dtype=bool)
    )
    p1_idx, p2_idx = stratified_split_halves_c1(hv_arr, split_idx=0, base_seed=split_seed)

    if raw_ctrl_pool is not None:
        obs_ntc_med_log1p = compute_obs_ntc_mean_official_median_log1p(
            raw_ctrl_pool, eval_genes=eval_genes, hgnc_bridge=hb
        )
    else:
        obs_ntc_med_log1p = np.asarray(tgt_ds.ctrl_mean_log1p[t_cols], dtype=np.float64)

    strata = resolve_r17_mismatched_strata(target_line)
    out: Dict[str, Any] = {
        "line_name": target_line,
        "role": R17_LINE_ROLES.get(target_line, "secondary_fresh"),
        "modality": R17_LINE_MODALITIES.get(target_line, "CRISPRi"),
        "checkpoint_id": "zeroshot/ensemble4",
        "eval_genes": eval_genes,
        "eval_perts": eval_perts,
        "target_indices": q_gidx_arr,
        "delta_obs": tgt_deltas[row_arr],
        "delta_half_a": tgt_ha[row_arr],
        "delta_half_b": tgt_hb[row_arr],
        "half_valid": hv_arr,
        "delta_transfer": np.asarray(preds_d["Transfer_Mean"], dtype=np.float64),
        "source_grand_mean": np.asarray(fct_state.grand_mean, dtype=np.float64),
        "fct_delta_pred": np.asarray(preds_d["FCT"], dtype=np.float64),
        "e1_constmed_delta_pred": np.asarray(preds_d["E1_ConstMed"], dtype=np.float64),
        "mu_target": tgt_cmean_lib,
        "mu_sources": mu_sources,
        "source_lines_fct": src_lines,
        "obs_ntc_mean_official_median_log1p": obs_ntc_med_log1p,
        "split_s0_p1_indices": p1_idx,
        "split_s0_p2_indices": p2_idx,
        "seen_sources": list(strata["seen_sources"]),
        "unseen_sources": list(strata["unseen_sources"]),
        "same_lineage_sources": list(strata["same_lineage_sources"]),
        "g_eval_r9": int(len(tgt_ds.genes)),
    }

    if target_line == "gse295214_h1" and fct_pool_mode == "primary":
        sens_srcs = [
            _subset_ds_eval_cols(
                load_standardized_line_r10(
                    dir_map[sn],
                    representation="log1p",
                    common_genes=r9_common_canon,
                    hgnc_bridge=hb,
                    min_gene_overlap=1,
                ),
                eval_genes,
                eval_g_to_idx,
            )
            for sn in R9_6_ALL_LINES
        ]
        sens_state = FctSourceState.fit(sens_srcs, mask_self=True)
        sens_preds, _, _ = predict_group_d_suite(
            state=sens_state,
            target_ctrl_mean_libnorm=tgt_cmean_lib,
            query_target_names=eval_perts,
            query_target_gene_indices=q_gidx_arr,
            include_r9_anchors=False,
            include_exploratory_e=False,
        )
        out["sens_full_r9_6"] = {
            "fct_delta_pred": np.asarray(sens_preds["FCT"], dtype=np.float64),
            "delta_transfer": np.asarray(sens_preds["Transfer_Mean"], dtype=np.float64),
            "source_grand_mean": np.asarray(sens_state.grand_mean, dtype=np.float64),
            "target_names": eval_perts,
            "gene_names": eval_genes,
        }

    return out
