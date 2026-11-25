#!/usr/bin/env python3
"""
Clean-Room Round 16 Arc STATE Inference Adapter (`A-T3`, `PREREG_I.md`, `R16_GATE2_SCHEMA.md`).

Gate-2 Locked Structural Behavior:
  - Hard freeze gate (`assert_freeze_gates_r16`): blocks any real-cell execution unless
    `PREREG_I_COMMIT.txt`, `CODE_FREEZE_R16.txt`, AND `R16_GO.txt` (issued >= `R16_GATE2_SCHEMA.md`)
    all exist and are non-empty.
  - `SE-600M`: loads directly from `model.safetensors` + `config.yaml` + `protein_embeddings.pt`,
    runs strictly in `float32` on `mps` (or `cpu`), never `fp16`/`bf16`.
  - `ST-SE`: decodes `2,058 -> 2,000` via trained `gene_decoder(preds)` and clips `[0.0, 14.0]`.
  - `ST-HVG`: reads `clip(preds, 0.0, 14.0)` directly (`gene_decoder` untrained across all 4 folds).
  - Grouping & batch covariates: `16 × 32` cell sets (`padded=True`), fixed `batch = 0` (`index_zero`)
    identical across matched and mismatched control arms (`bypass_zero`/`mean_batch_emb`/
    `multi_batch_ens8`/custom `batch` index available strictly for `§4.4` secondary checks).
  - Predictions: `delta_pred(p) = mean_cells(clip(Y(x, p, b))) - mean_cells(clip(Y(x, "non-targeting", b)))`
    evaluated on identical sampled control cell sets.
  - Gene & perturbation symbols: aligned via frozen R10 `HGNCBridge`.

Never reads or imports `round16/src/`, and never opens `adata_pred.h5ad`.
"""

from __future__ import annotations

import os
import pickle
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import torch
import yaml

os.environ.setdefault("OMP_NUM_THREADS", "1")

AUDIT_DIR = Path(__file__).resolve().parent
REPO_ROOT = AUDIT_DIR.parents[1]
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
for _p in (str(AUDIT_DIR), str(R10_SRC_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from audit_r16_sampling import (  # noqa: E402
    CELL_SET_LEN_DEFAULT,
    DEFAULT_R16_SEED,
    MIN_CTRL_FLOOR_DEFAULT,
    N_CTRL_DEFAULT,
    REPLOGLE4_ALL_LINES,
    libnorm_log1p,
    sample_stratified_control_sets_r16,
)
from data_loader_r10 import DEFAULT_HGNC_FILE, HGNCBridge  # noqa: E402

MODELS_ROOT_DEFAULT = REPO_ROOT / "round16" / "shared" / "models"
SE600M_DEFAULT_DIR = MODELS_ROOT_DEFAULT / "SE-600M"
RAW_COUNT_HEURISTIC_THRESHOLD = 35.0
EXPONENTIATED_UMIS_LIMIT = 5_000_000.0
CLIP_MAX_EXPR = 14.0

_CACHED_HGNC_BRIDGE: Optional[HGNCBridge] = None


class BlindingGateError(RuntimeError):
    """Raised if real-cell inference is attempted before PREREG_I lock, R16 freeze, & R16_GO."""


class TrainingLeakageError(RuntimeError):
    """Raised if target line overlaps checkpoint training set."""


def get_hgnc_bridge(
    hgnc_bridge: Optional[HGNCBridge] = None,
    hgnc_path: Path = DEFAULT_HGNC_FILE,
) -> HGNCBridge:
    """Return provided or lazily loaded frozen R10 `HGNCBridge`."""
    global _CACHED_HGNC_BRIDGE
    if hgnc_bridge is not None:
        return hgnc_bridge
    if _CACHED_HGNC_BRIDGE is None:
        _CACHED_HGNC_BRIDGE = HGNCBridge.from_tsv(hgnc_path)
    return _CACHED_HGNC_BRIDGE


def resolve_compute_device(device: str = "auto") -> str:
    """Resolve `'auto'` -> `'mps'` when available on Apple Silicon, else `'cpu'`."""
    dev_s = str(device).strip().lower()
    if dev_s == "auto":
        return "mps" if torch.backends.mps.is_available() else "cpu"
    if dev_s == "mps" and not torch.backends.mps.is_available():
        raise RuntimeError("Requested device='mps', but PyTorch MPS backend is unavailable.")
    if dev_s not in ("cpu", "mps"):
        raise ValueError(f"Unsupported device={device!r}; expected 'auto', 'mps', or 'cpu'.")
    return dev_s


def _ckpt_short_name(checkpoint_id: str) -> str:
    ck = str(checkpoint_id).strip().lower()
    for s in ("hepg2", "jurkat", "k562", "rpe1", "ensemble4"):
        if ck == s or ck.endswith(f"/{s}"):
            return s
    raise ValueError(f"Unsupported checkpoint_id={checkpoint_id!r}")


def assert_zero_training_leakage(
    model_key: str,
    target_line: str,
    checkpoint_id: str,
) -> None:
    """Enforce zero training leakage between `checkpoint_id` and `target_line` (`PREREG_I §2`)."""
    short = _ckpt_short_name(checkpoint_id)
    if short == "ensemble4":
        trained = set(REPLOGLE4_ALL_LINES)
    else:
        heldout = {
            "hepg2": "nadig_hepg2",
            "jurkat": "nadig_jurkat",
            "k562": "replogle_k562_essential",
            "rpe1": "replogle_rpe1",
        }[short]
        trained = set(REPLOGLE4_ALL_LINES) - {heldout}
    aliases = {
        "nadig_hepg2": {"nadig_hepg2", "hepg2"},
        "nadig_jurkat": {"nadig_jurkat", "jurkat"},
        "replogle_k562_essential": {"replogle_k562_essential", "k562"},
        "replogle_rpe1": {"replogle_rpe1", "rpe1"},
    }
    tgt = str(target_line).lower()
    for tr in trained:
        if tgt in aliases[tr]:
            raise TrainingLeakageError(
                f"{model_key} ({checkpoint_id}) was trained on {tr}; cannot evaluate {target_line!r}."
            )


def _read_nonempty_gate_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8").strip()


def assert_freeze_gates_r16(status_dir: Path) -> Dict[str, str]:
    """
    Hard runtime blinding gate (`R16_CONTROLLER_NOTES.md`, `R16_GATE2_SCHEMA.md`):
    refuse any real-cell run unless ALL three exist and are non-empty:
      1. `PREREG_I_COMMIT.txt`
      2. `CODE_FREEZE_R16.txt`
      3. `R16_GO.txt` (with mtime >= `R16_GATE2_SCHEMA.md` when `R16_GATE2_SCHEMA.md` is present).
    """
    sdir = Path(status_dir)
    prereg_lock = sdir / "PREREG_I_COMMIT.txt"
    freeze_lock = sdir / "CODE_FREEZE_R16.txt"
    go_lock = sdir / "R16_GO.txt"
    gate2_schema = sdir / "R16_GATE2_SCHEMA.md"

    prereg_txt = _read_nonempty_gate_text(prereg_lock)
    if not prereg_txt:
        raise BlindingGateError(f"Missing or empty PREREG_I lock file: {prereg_lock}")

    freeze_txt = _read_nonempty_gate_text(freeze_lock)
    if not freeze_txt:
        raise BlindingGateError(f"Missing or empty Round 16 code freeze file: {freeze_lock}")

    go_txt = _read_nonempty_gate_text(go_lock)
    if not go_txt:
        raise BlindingGateError(f"Missing or empty controller GO authorization file: {go_lock}")

    if gate2_schema.exists():
        if go_lock.stat().st_mtime < gate2_schema.stat().st_mtime:
            raise BlindingGateError(
                f"Stale controller GO file ({go_lock}) predates Gate 2 schema ({gate2_schema})."
            )

    return {
        "prereg_commit": prereg_txt.splitlines()[0].strip(),
        "code_freeze": freeze_txt.splitlines()[0].strip(),
        "r16_go": go_txt.splitlines()[0].strip(),
    }


def preprocess_hvg_expression(
    counts: np.ndarray,
    mode: str = "libnorm1e4_log1p",
    cell_total_umis: Optional[np.ndarray] = None,
    target_median_umi: Optional[float] = None,
) -> np.ndarray:
    """
    Convert non-negative count matrix `(N, G)` into `log1p`-normalized gene expression:
      - `'libnorm1e4_log1p'`: `log1p(1e4 * counts / cell_total_umi)`
      - `'official_median_log1p'`: matches `scanpy.pp.normalize_total(adata)` (`target_sum = median`)
        + `scanpy.pp.log1p(adata)` (`state/_cli/_tx/_preprocess_train.py:46-49`).
    """
    x = np.asarray(counts, dtype=np.float64)
    if x.ndim != 2:
        raise ValueError(f"Expected 2D matrix, got shape {x.shape}")
    if cell_total_umis is not None:
        totals = np.asarray(cell_total_umis, dtype=np.float64)
        if totals.shape != (x.shape[0],):
            raise ValueError(f"cell_total_umis shape {totals.shape} != ({x.shape[0]},)")
    else:
        totals = np.sum(x, axis=1)

    if mode == "libnorm1e4_log1p":
        target_sum = 1e4
    elif mode == "official_median_log1p":
        if target_median_umi is not None:
            target_sum = float(target_median_umi)
        else:
            pos = totals[totals > 0.0]
            target_sum = float(np.median(pos)) if pos.size > 0 else 1.0
    else:
        raise ValueError(f"Unsupported HVG preprocessing mode={mode!r}")

    safe_tot = np.where(totals > 0.0, totals, 1.0)[:, None]
    normed = np.where(totals[:, None] > 0.0, (x / safe_tot) * target_sum, 0.0)
    return np.log1p(np.maximum(normed, 0.0)).astype(np.float32)


def load_state_checkpoint_bundle(
    repo_name: str,
    ckpt_short: str,
    models_root: Path = MODELS_ROOT_DEFAULT,
    device: str = "cpu",
) -> Dict[str, Any]:
    """Load `StateTransitionPerturbationModel` and structural vocabularies cleanly from `models_root`."""
    from state.tx.models.state_transition import StateTransitionPerturbationModel

    ckpt_dir = Path(models_root) / repo_name / "zeroshot" / ckpt_short
    cfg_path = ckpt_dir / "config.yaml"
    vd_path = ckpt_dir / "var_dims.pkl"
    pert_path = ckpt_dir / "pert_onehot_map.pt"
    ckpt_path = ckpt_dir / "checkpoints" / "final.ckpt"

    for p in (cfg_path, vd_path, pert_path, ckpt_path):
        if not p.exists():
            raise FileNotFoundError(f"Missing STATE artifact: {p}")

    with open(cfg_path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    with open(vd_path, "rb") as f:
        var_dims = pickle.load(f)
    pert_map = torch.load(pert_path, map_location="cpu", weights_only=False)
    pert_lookup = {str(k): v.float() for k, v in pert_map.items()}

    model = StateTransitionPerturbationModel.load_from_checkpoint(
        str(ckpt_path), map_location="cpu", weights_only=False
    )
    model.eval()
    dev = torch.device(resolve_compute_device(device))
    model.to(dev, dtype=torch.float32)

    return {
        "repo": repo_name,
        "ckpt_short": ckpt_short,
        "cfg": cfg,
        "var_dims": var_dims,
        "gene_names": [str(g) for g in var_dims["gene_names"]],
        "pert_lookup": pert_lookup,
        "control_pert": str(cfg.get("data", {}).get("kwargs", {}).get("control_pert", "non-targeting")),
        "cell_set_len": int(getattr(model, "cell_sentence_len", 32)),
        "model": model,
        "device": dev,
    }


def _build_source_symbol_index(
    source_genes: Sequence[str],
    source_ensembl_ids: Optional[Sequence[str]] = None,
    hgnc_bridge: Optional[HGNCBridge] = None,
    mean_scores: Optional[np.ndarray] = None,
) -> Tuple[Dict[str, int], Dict[str, int]]:
    """
    Build deterministic `(raw_sym -> best_col, canon_sym -> best_col)` maps for `source_genes`.
    Collisions resolve by `(is_exact_canonical, mean_score, -col_idx)`.
    """
    raw_best: Dict[str, Tuple[Tuple[float, int], int]] = {}
    canon_best: Dict[str, Tuple[Tuple[int, float, int], int]] = {}

    for idx, sg in enumerate(source_genes):
        s = str(sg).strip()
        if not s:
            continue
        sc = float(mean_scores[idx]) if mean_scores is not None else 0.0
        r_prev = raw_best.get(s)
        r_key = (sc, -idx)
        if r_prev is None or r_key > r_prev[0]:
            raw_best[s] = (r_key, idx)

        if hgnc_bridge is not None:
            ens = source_ensembl_ids[idx] if source_ensembl_ids is not None else None
            canon = hgnc_bridge.canonicalize(s, ensembl_id=ens)
        else:
            canon = s
        if canon:
            c_key = (1 if s == canon else 0, sc, -idx)
            c_prev = canon_best.get(canon)
            if c_prev is None or c_key > c_prev[0]:
                canon_best[canon] = (c_key, idx)

    return {k: v[1] for k, v in raw_best.items()}, {k: v[1] for k, v in canon_best.items()}


def align_hvg_expression_to_model(
    ctrl_expr_log1p: np.ndarray,
    source_genes: Sequence[str],
    model_genes: Sequence[str],
    source_ensembl_ids: Optional[Sequence[str]] = None,
    hgnc_bridge: Optional[HGNCBridge] = None,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    Align control `log1p` expression onto `model_genes` (`2,000` HVG table), matching exact
    raw symbols first and rescuing HGNC legacy/alias symbols via frozen R10 `HGNCBridge`.
    """
    col_means = np.mean(np.asarray(ctrl_expr_log1p, dtype=np.float64), axis=0) if ctrl_expr_log1p.shape[0] > 0 else None
    raw_pos, canon_pos = _build_source_symbol_index(
        source_genes,
        source_ensembl_ids=source_ensembl_ids,
        hgnc_bridge=hgnc_bridge,
        mean_scores=col_means,
    )
    n_cells = ctrl_expr_log1p.shape[0]
    aligned = np.zeros((n_cells, len(model_genes)), dtype=np.float32)

    matched_idx: List[Tuple[int, int]] = []
    n_raw_exact = 0
    n_canon_rescued = 0
    for j, mg in enumerate(model_genes):
        mg_s = str(mg).strip()
        si = raw_pos.get(mg_s)
        if si is not None:
            matched_idx.append((j, si))
            n_raw_exact += 1
            continue
        if hgnc_bridge is not None:
            mg_c = hgnc_bridge.canonicalize(mg_s)
            if mg_c and mg_c in canon_pos:
                matched_idx.append((j, canon_pos[mg_c]))
                n_canon_rescued += 1

    if matched_idx:
        dst_cols = np.asarray([p[0] for p in matched_idx], dtype=np.int64)
        src_cols = np.asarray([p[1] for p in matched_idx], dtype=np.int64)
        aligned[:, dst_cols] = np.asarray(ctrl_expr_log1p[:, src_cols], dtype=np.float32)

    stats = {
        "n_model_genes": len(model_genes),
        "n_mapped_genes": len(matched_idx),
        "n_raw_exact": n_raw_exact,
        "n_hgnc_canon_rescued": n_canon_rescued,
        "mapping_coverage": float(len(matched_idx)) / float(max(len(model_genes), 1)),
    }
    return aligned, stats


def resolve_pert_lookup_keys(
    eval_perts: Sequence[str],
    pert_lookup: Dict[str, torch.Tensor],
    control_pert: str = "non-targeting",
    hgnc_bridge: Optional[HGNCBridge] = None,
) -> Dict[str, str]:
    """
    Map each requested target perturbation name to its key inside `pert_onehot_map.pt`,
    harmonizing legacy Replogle-4 symbols (`AARS`, `H2AFZ`, `TAZ`) with canonical HGNC
    symbols (`AARS1`, `H2AZ1`, `TAFAZZIN`) via `HGNCBridge`.
    """
    canon_to_raw: Dict[str, str] = {}
    if hgnc_bridge is not None:
        # Prefer exact match when key already equals its canonical symbol
        for rk in sorted(pert_lookup.keys()):
            if rk == control_pert:
                continue
            ck = hgnc_bridge.canonicalize(rk)
            if ck and (ck not in canon_to_raw or rk == ck):
                canon_to_raw[ck] = rk

    resolved: Dict[str, str] = {}
    missing: List[str] = []
    for p in eval_perts:
        ps = str(p).strip()
        if ps in pert_lookup:
            resolved[ps] = ps
            continue
        if hgnc_bridge is not None:
            cp = hgnc_bridge.canonicalize(ps)
            if cp and cp in canon_to_raw:
                resolved[ps] = canon_to_raw[cp]
                continue
        missing.append(ps)
    if missing:
        raise KeyError(f"Perturbations absent from pert_onehot_map.pt (even after HGNCBridge): {missing[:5]}")
    return resolved


def map_genes_to_se600m_pe(
    source_genes: Sequence[str],
    pe_keys: Sequence[str],
    source_ensembl_ids: Optional[Sequence[str]] = None,
    mean_counts: Optional[np.ndarray] = None,
    hgnc_bridge: Optional[HGNCBridge] = None,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, int]]:
    """
    Align source transcriptome genes onto `SE-600M` `protein_embeddings.pt` (`19,790` genes)
    using exact raw match first and `HGNCBridge` canonical fallback second, deduplicating
    any source-column collisions onto unique global PE indices cleanly.
    """
    pe_raw_to_idx: Dict[str, int] = {str(k): i for i, k in enumerate(pe_keys)}
    pe_canon_to_idx: Dict[str, int] = {}
    if hgnc_bridge is not None:
        for i, k in enumerate(pe_keys):
            ks = str(k)
            ck = hgnc_bridge.canonicalize(ks)
            if ck and (ck not in pe_canon_to_idx or ks == ck):
                pe_canon_to_idx[ck] = i

    # For each PE global index, pick at most one source column (exact > canon > higher mean count > earlier idx)
    best_by_pe: Dict[int, Tuple[Tuple[int, float, int], int, str]] = {}
    for col_i, sg in enumerate(source_genes):
        s = str(sg).strip()
        if not s:
            continue
        sc = float(mean_counts[col_i]) if mean_counts is not None else 0.0
        pe_idx = pe_raw_to_idx.get(s)
        match_kind = "raw_exact"
        pri = 2
        if pe_idx is None and hgnc_bridge is not None:
            ens = source_ensembl_ids[col_i] if source_ensembl_ids is not None else None
            ck = hgnc_bridge.canonicalize(s, ensembl_id=ens)
            if ck and ck in pe_canon_to_idx:
                pe_idx = pe_canon_to_idx[ck]
                match_kind = "hgnc_canon"
                pri = 1
        if pe_idx is None:
            continue
        cand_key = (pri, sc, -col_i)
        prev = best_by_pe.get(pe_idx)
        if prev is None or cand_key > prev[0]:
            best_by_pe[pe_idx] = (cand_key, col_i, match_kind)

    chosen = sorted(best_by_pe.items(), key=lambda kv: kv[1][1])
    src_cols = np.asarray([kv[1][1] for kv in chosen], dtype=np.int64)
    pe_global = np.asarray([kv[0] for kv in chosen], dtype=np.int64)
    n_raw = sum(1 for kv in chosen if kv[1][2] == "raw_exact")
    n_canon = sum(1 for kv in chosen if kv[1][2] == "hgnc_canon")
    stats = {
        "n_source_genes": len(source_genes),
        "n_pe_total": len(pe_keys),
        "n_mapped_genes": len(chosen),
        "n_raw_exact": n_raw,
        "n_hgnc_canon_rescued": n_canon,
    }
    return src_cols, pe_global, stats


def _is_raw_umi_row(counts_row: torch.Tensor) -> bool:
    """Exact replica of `state/emb/data/loader.py:431-453` `is_raw_integer_counts()`."""
    if float(torch.max(counts_row).item()) > RAW_COUNT_HEURISTIC_THRESHOLD:
        return True
    if int(torch.expm1(counts_row).sum().item()) > EXPONENTIATED_UMIS_LIMIT:
        return True
    return False


def encode_with_loaded_se_model(
    model: torch.nn.Module,
    pe_keys: Sequence[str],
    ctrl_counts: np.ndarray,
    source_genes: Sequence[str],
    source_ensembl_ids: Optional[Sequence[str]] = None,
    hgnc_bridge: Optional[HGNCBridge] = None,
    pad_length: int = 2048,
    cls_token_idx: int = 3,
    batch_size: int = 16,
    seed: int = DEFAULT_R16_SEED,
    device: Optional[torch.device] = None,
) -> Tuple[np.ndarray, Dict[str, int]]:
    """
    Encode control cells through a loaded `StateEmbeddingModel` strictly matching
    `arc-state 0.11.1` (`state/emb/data/loader.py:457-629`, `state/emb/nn/model.py:209-246, 288-332`,
    `state/emb/inference.py:172-185`) in `float32` using deterministic `Generator` tie-breaking.
    Returns `(X_state [N, output_dim + z_dim_ds], mapping_stats)`.
    """
    x_arr = np.asarray(ctrl_counts, dtype=np.float32)
    if x_arr.ndim != 2:
        raise ValueError(f"Expected 2D ctrl_counts, got shape {x_arr.shape}")
    if x_arr.shape[1] != len(source_genes):
        raise ValueError(f"ctrl_counts cols ({x_arr.shape[1]}) != len(source_genes) ({len(source_genes)})")

    mean_c = np.mean(x_arr, axis=0, dtype=np.float64) if x_arr.shape[0] > 0 else None
    src_cols, pe_global_np, map_stats = map_genes_to_se600m_pe(
        source_genes=source_genes,
        pe_keys=pe_keys,
        source_ensembl_ids=source_ensembl_ids,
        mean_counts=mean_c,
        hgnc_bridge=hgnc_bridge,
    )
    if src_cols.size < 2:
        raise ValueError(f"Too few source genes mapped to SE-600M protein embeddings: {map_stats}")

    src_cols_t = torch.as_tensor(src_cols, dtype=torch.long)
    pe_global_t = torch.as_tensor(pe_global_np, dtype=torch.long)
    n_mapped = int(src_cols_t.numel())
    n_cells = x_arr.shape[0]
    dev = device or next(model.parameters()).device

    if dev.type == "mps":
        for p in model.parameters():
            if p.dtype != torch.float32:
                raise ValueError(f"SE-600M on MPS requires strict float32, found {p.dtype}")

    gen = torch.Generator(device="cpu")
    gen.manual_seed(int(seed))

    # Build tokenized cell sentences identically to VCIDatasetSentenceCollator.sample_cell_sentences
    cell_sentences = torch.zeros((n_cells, pad_length), dtype=torch.long)
    cell_counts_pct = torch.zeros((n_cells, pad_length), dtype=torch.float32)

    for c in range(n_cells):
        raw_full = torch.from_numpy(x_arr[c : c + 1]).to(dtype=torch.float32)
        if torch.isnan(raw_full).any():
            raise ValueError(f"NaN detected in control cell row {c}")
        if torch.any(raw_full < 0.0):
            raw_full = torch.relu(raw_full)
        if _is_raw_umi_row(raw_full):
            raw_full = torch.log1p(raw_full)

        cell = raw_full[0, src_cols_t]
        csum = float(torch.sum(cell).item())
        expr_w = (cell / csum) if csum > 0.0 else torch.softmax(cell, dim=0)

        num_pos = int(torch.sum(cell > 0.0).item())
        perm = torch.randperm(n_mapped, generator=gen)
        shuffled = cell[perm]
        ranked_local = perm[torch.argsort(shuffled, descending=True)]

        local_tokens = torch.zeros(pad_length, dtype=torch.long)
        safe_cls = int(cls_token_idx) if 0 <= int(cls_token_idx) < n_mapped else 0
        local_tokens[0] = safe_cls
        if n_mapped >= pad_length - 1:
            local_tokens[1:] = ranked_local[: pad_length - 1]
        else:
            n_nz = min(num_pos, pad_length - 1)
            if n_nz > 0:
                local_tokens[1 : n_nz + 1] = ranked_local[:n_nz]
            rem = pad_length - 1 - n_nz
            if rem > 0:
                unexp = ranked_local[n_nz:]
                fill_pool = unexp if len(unexp) > 0 else ranked_local
                pick = torch.randint(len(fill_pool), (rem,), generator=gen)
                local_tokens[n_nz + 1 :] = fill_pool[pick]

        cell_counts_pct[c, :] = 100.0 * expr_w[local_tokens]
        cell_sentences[c, :] = pe_global_t[local_tokens]

    outs: List[np.ndarray] = []
    use_ds_tok = getattr(model, "dataset_token", None) is not None
    bs_step = max(1, int(batch_size))

    model.eval()
    with torch.no_grad():
        for i0 in range(0, n_cells, bs_step):
            i1 = min(i0 + bs_step, n_cells)
            bs = i1 - i0
            tok_idx = cell_sentences[i0:i1].to(dev)
            cnt_blk = cell_counts_pct[i0:i1].to(device=dev, dtype=torch.float32)

            src_emb = model.pe_embedding(tok_idx).to(dtype=torch.float32)
            src_emb = torch.nn.functional.normalize(src_emb, dim=2)
            src_emb[:, 0, :] = model.cls_token.expand(bs, -1)

            mask = torch.zeros((bs, pad_length), dtype=torch.bool, device=dev)
            if use_ds_tok:
                ds_tok = model.dataset_token.expand(bs, -1).unsqueeze(1)
                src_emb = torch.cat((src_emb, ds_tok), dim=1)
                mask = torch.cat((mask, torch.zeros((bs, 1), dtype=torch.bool, device=dev)), dim=1)

            _, emb, ds_emb = model.forward(src_emb, mask=mask, counts=cnt_blk, dataset_nums=None)
            emb_np = emb.detach().cpu().float().numpy()
            if ds_emb is not None and getattr(model, "dataset_embedder", None) is not None:
                ds_proj = model.dataset_embedder(ds_emb).detach().cpu().float().numpy()
                block = np.concatenate([emb_np, ds_proj], axis=-1).astype(np.float32)
            else:
                block = emb_np.astype(np.float32)
            outs.append(block)

    return np.concatenate(outs, axis=0), map_stats


class CleanRoomSE600MEncoder:
    """
    Direct `safetensors` loader + Apple Silicon `MPS float32` encoder for `arcinstitute/SE-600M`
    (`R16_GATE2_SCHEMA.md §4.1`). Never invokes `torch.load` on `model.safetensors` or `fp16`/`bf16` on MPS.
    """

    def __init__(
        self,
        se600m_dir: Path = SE600M_DEFAULT_DIR,
        device: str = "auto",
        batch_size: int = 16,
        hgnc_bridge: Optional[HGNCBridge] = None,
    ) -> None:
        self.se600m_dir = Path(se600m_dir)
        self.device_str = resolve_compute_device(device)
        self.device = torch.device(self.device_str)
        self.batch_size = int(batch_size)
        self.hgnc_bridge = hgnc_bridge
        self._model: Optional[torch.nn.Module] = None
        self._cfg: Any = None
        self._pe_keys: List[str] = []

    def load(self) -> None:
        if self._model is not None:
            return
        from omegaconf import OmegaConf
        from safetensors.torch import load_file as safetensors_load_file
        from state.emb.nn.model import StateEmbeddingModel

        cfg_path = self.se600m_dir / "config.yaml"
        pe_path = self.se600m_dir / "protein_embeddings.pt"
        st_path = self.se600m_dir / "model.safetensors"
        for p in (cfg_path, pe_path, st_path):
            if not p.exists():
                raise FileNotFoundError(f"Missing SE-600M artifact: {p}")

        cfg = OmegaConf.load(str(cfg_path))
        pe_dict = torch.load(pe_path, map_location="cpu", weights_only=False)
        pe_keys = [str(k) for k in pe_dict.keys()]
        pe_mat = torch.vstack([v.float() for v in pe_dict.values()])

        model = StateEmbeddingModel(
            token_dim=int(cfg.tokenizer.token_dim),
            d_model=int(cfg.model.emsize),
            nhead=int(cfg.model.nhead),
            d_hid=int(cfg.model.d_hid),
            nlayers=int(cfg.model.nlayers),
            output_dim=int(cfg.model.output_dim),
            dropout=0.0,
            cfg=cfg,
        )
        sd = safetensors_load_file(str(st_path), device="cpu")
        missing, _ = model.load_state_dict(sd, strict=False)
        if missing:
            raise RuntimeError(f"SE-600M safetensors load missing keys: {missing[:10]}")

        model.pe_embedding = torch.nn.Embedding.from_pretrained(pe_mat)
        model.eval()
        model.to(self.device, dtype=torch.float32)
        model.pe_embedding.to(self.device, dtype=torch.float32)

        self._model = model
        self._cfg = cfg
        self._pe_keys = pe_keys

    def encode_cells(
        self,
        ctrl_counts: np.ndarray,
        source_genes: Sequence[str],
        source_ensembl_ids: Optional[Sequence[str]] = None,
        seed: int = DEFAULT_R16_SEED,
        hgnc_bridge: Optional[HGNCBridge] = None,
    ) -> Tuple[np.ndarray, Dict[str, int]]:
        self.load()
        assert self._model is not None and self._cfg is not None
        hb = get_hgnc_bridge(hgnc_bridge or self.hgnc_bridge)
        return encode_with_loaded_se_model(
            model=self._model,
            pe_keys=self._pe_keys,
            ctrl_counts=ctrl_counts,
            source_genes=source_genes,
            source_ensembl_ids=source_ensembl_ids,
            hgnc_bridge=hb,
            pad_length=int(self._cfg.dataset.pad_length),
            cls_token_idx=int(self._cfg.dataset.cls_token_idx),
            batch_size=self.batch_size,
            seed=seed,
            device=self.device,
        )


class CleanRoomArcStateAdapter:
    """
    Independent adapter running `ST-SE-Replogle` (with `SE-600M`) and `ST-HVG-Replogle`
    across matched vs. mismatched sampled `16 × 32` control sets under Gate-2 locked rules.
    """

    def __init__(
        self,
        model_key: str = "STATE_SE",
        checkpoint_id: str = "zeroshot/hepg2",
        models_root: Path = MODELS_ROOT_DEFAULT,
        status_dir: Optional[Path] = None,
        device: str = "auto",
        neutral_batch_mode: str = "index_zero",
        delta_mode: str = "minus_pred_ntc",
        hvg_preproc_mode: str = "libnorm1e4_log1p",
        allow_synthetic_smoke: bool = False,
        batch_size_sets: int = 16,
        se600m_encoder: Optional[CleanRoomSE600MEncoder] = None,
        hgnc_bridge: Optional[HGNCBridge] = None,
        custom_batch_idx: Optional[int] = None,
    ) -> None:
        if model_key not in ("STATE_HVG", "STATE_SE"):
            raise ValueError(f"Unsupported model_key={model_key!r}")
        if neutral_batch_mode not in ("index_zero", "bypass_zero", "mean_batch_emb", "multi_batch_ens8"):
            raise ValueError(f"Unsupported neutral_batch_mode={neutral_batch_mode!r}")
        if delta_mode not in ("minus_pred_ntc", "minus_input_ctrl"):
            raise ValueError(f"Unsupported delta_mode={delta_mode!r}")
        if hvg_preproc_mode not in ("libnorm1e4_log1p", "official_median_log1p"):
            raise ValueError(f"Unsupported hvg_preproc_mode={hvg_preproc_mode!r}")

        self.model_key = model_key
        self.repo_name = "ST-HVG-Replogle" if model_key == "STATE_HVG" else "ST-SE-Replogle"
        self.checkpoint_id = checkpoint_id
        self.ckpt_short = _ckpt_short_name(checkpoint_id)
        self.models_root = Path(models_root)
        self.status_dir = Path(status_dir) if status_dir is not None else (REPO_ROOT / "round16" / "status")
        self.device = resolve_compute_device(device)
        self.neutral_batch_mode = neutral_batch_mode
        self.delta_mode = delta_mode
        self.hvg_preproc_mode = hvg_preproc_mode
        self.allow_synthetic_smoke = bool(allow_synthetic_smoke)
        self.batch_size_sets = int(batch_size_sets)
        self.hgnc_bridge = hgnc_bridge
        self.custom_batch_idx = custom_batch_idx
        self._se600m_encoder = se600m_encoder
        self._bundles: Dict[str, Dict[str, Any]] = {}

    def _ckpt_list(self) -> List[str]:
        if self.ckpt_short == "ensemble4":
            return ["hepg2", "jurkat", "k562", "rpe1"]
        return [self.ckpt_short]

    def _get_se_encoder(self) -> CleanRoomSE600MEncoder:
        if self._se600m_encoder is None:
            self._se600m_encoder = CleanRoomSE600MEncoder(
                se600m_dir=self.models_root / "SE-600M",
                device=self.device,
                batch_size=16,
                hgnc_bridge=self.hgnc_bridge,
            )
        return self._se600m_encoder

    def load_weights(self) -> None:
        if self.model_key == "STATE_SE":
            se_safetensors = self.models_root / "SE-600M" / "model.safetensors"
            se_prot = self.models_root / "SE-600M" / "protein_embeddings.pt"
            if not (se_safetensors.exists() and se_prot.exists()):
                raise FileNotFoundError(
                    "SE-600M weights incomplete (`model.safetensors` & `protein_embeddings.pt` required): "
                    f"{self.models_root / 'SE-600M'}"
                )
        for ck in self._ckpt_list():
            if ck not in self._bundles:
                self._bundles[ck] = load_state_checkpoint_bundle(
                    self.repo_name, ck, models_root=self.models_root, device=self.device
                )

    def _predict_single_batch_index(
        self,
        model: torch.nn.Module,
        dev: torch.device,
        pert_lookup: Dict[str, torch.Tensor],
        query_raw_keys: Sequence[str],
        ctrl_tensor: torch.Tensor,
        batch_mode: str,
        fixed_batch_idx: int,
        n_genes: int,
    ) -> np.ndarray:
        """Run all `query_raw_keys` across `ctrl_tensor` `(S, 32, E)` under one neutral batch scheme."""
        n_sets, set_len, e_in = ctrl_tensor.shape
        n_cells = n_sets * set_len
        out = np.empty((len(query_raw_keys), n_cells, n_genes), dtype=np.float32)

        orig_batch_enc = getattr(model, "batch_encoder", None)
        if batch_mode == "bypass_zero":
            model.batch_encoder = None
        elif batch_mode == "mean_batch_emb" and orig_batch_enc is not None:
            mean_w = orig_batch_enc.weight.detach().mean(dim=0, keepdim=True)

            class _MeanBatchEncoder(torch.nn.Module):
                def __init__(self, vec: torch.Tensor) -> None:
                    super().__init__()
                    self.register_buffer("vec", vec)

                def forward(self, b_idx: torch.Tensor) -> torch.Tensor:
                    return self.vec.view(1, 1, -1).expand(b_idx.shape[0], b_idx.shape[1], -1)

            model.batch_encoder = _MeanBatchEncoder(mean_w.to(dev))

        try:
            with torch.no_grad():
                chunk = max(1, self.batch_size_sets)
                for q_i, pkey in enumerate(query_raw_keys):
                    pvec = pert_lookup[pkey].to(device=dev, dtype=torch.float32)
                    outs: List[np.ndarray] = []
                    for s0 in range(0, n_sets, chunk):
                        s1 = min(s0 + chunk, n_sets)
                        bs = s1 - s0
                        c_blk = ctrl_tensor[s0:s1].reshape(bs * set_len, e_in)
                        p_blk = pvec.unsqueeze(0).expand(bs * set_len, -1)
                        batch: Dict[str, torch.Tensor] = {
                            "ctrl_cell_emb": c_blk,
                            "pert_emb": p_blk,
                        }
                        if batch_mode in ("index_zero", "mean_batch_emb", "multi_batch_ens8"):
                            bval = int(fixed_batch_idx)
                            batch["batch"] = torch.full((bs * set_len,), bval, dtype=torch.long, device=dev)

                        step_out = model.predict_step(batch, batch_idx=0, padded=True)
                        if self.model_key == "STATE_SE":
                            decoded = step_out.get("pert_cell_counts_preds")
                            if decoded is None:
                                decoded = model.gene_decoder(step_out["preds"])
                            arr = decoded.detach().cpu().numpy()
                        else:
                            # ST-HVG gene_decoder is untrained random init; ALWAYS read direct preds!
                            arr = step_out["preds"].detach().cpu().numpy()

                        if arr.shape[1] != n_genes:
                            raise RuntimeError(
                                f"Output gene dim mismatch for {self.model_key}: got {arr.shape[1]}, expected {n_genes}"
                            )
                        outs.append(np.clip(arr.astype(np.float32, copy=False), 0.0, CLIP_MAX_EXPR))
                    out[q_i] = np.concatenate(outs, axis=0)
        finally:
            model.batch_encoder = orig_batch_enc

        return out

    def _run_forward_sets_single_ckpt(
        self,
        bundle: Dict[str, Any],
        ctrl_sets_embed: np.ndarray,
        pert_names: Sequence[str],
        hgnc_bridge: Optional[HGNCBridge] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Evaluate a single checkpoint on `ctrl_sets_embed` of shape `(n_sets, 32, E_in)`
        for each perturbation in `pert_names` plus `'non-targeting'`.
        """
        model = bundle["model"]
        dev = bundle["device"]
        pert_lookup: Dict[str, torch.Tensor] = bundle["pert_lookup"]
        ntc_key = bundle["control_pert"]
        n_sets, set_len, e_in = ctrl_sets_embed.shape
        expected_set_len = int(bundle["cell_set_len"])
        if set_len != expected_set_len:
            raise ValueError(f"Control set length {set_len} != model cell_sentence_len {expected_set_len}")

        n_genes = len(bundle["gene_names"])
        ctrl_tensor = torch.tensor(ctrl_sets_embed, dtype=torch.float32, device=dev)
        if ntc_key not in pert_lookup:
            raise KeyError(f"Control key {ntc_key!r} absent from pert_onehot_map.pt")

        p_resolved = resolve_pert_lookup_keys(
            pert_names, pert_lookup, control_pert=ntc_key, hgnc_bridge=hgnc_bridge
        )
        query_raw_keys = [p_resolved[str(p).strip()] for p in pert_names] + [ntc_key]

        if self.neutral_batch_mode == "multi_batch_ens8":
            indices = (0, 7, 14, 21, 28, 35, 42, 49)
            acc = None
            for bi in indices:
                res_b = self._predict_single_batch_index(
                    model, dev, pert_lookup, query_raw_keys, ctrl_tensor, "multi_batch_ens8", bi, n_genes
                )
                acc = res_b.astype(np.float64) if acc is None else (acc + res_b.astype(np.float64))
            pred_gene_rows = (acc / float(len(indices))).astype(np.float32)
        else:
            bidx = int(self.custom_batch_idx) if self.custom_batch_idx is not None else 0
            pred_gene_rows = self._predict_single_batch_index(
                model, dev, pert_lookup, query_raw_keys, ctrl_tensor, self.neutral_batch_mode, bidx, n_genes
            )

        return pred_gene_rows[:-1], pred_gene_rows[-1]

    def predict_arm_deltas(
        self,
        target_line: str,
        control_source_line: str,
        ctrl_matrix: Optional[np.ndarray],
        source_genes: Sequence[str],
        eval_perts: Sequence[str],
        batch_labels: Optional[Sequence[Any]] = None,
        qc_mask: Optional[Sequence[bool]] = None,
        input_is_log1p: bool = False,
        seed: int = DEFAULT_R16_SEED,
        n_ctrl: int = N_CTRL_DEFAULT,
        cell_set_len: int = CELL_SET_LEN_DEFAULT,
        is_synthetic_data: bool = False,
        source_ensembl_ids: Optional[Sequence[str]] = None,
        cell_total_umis: Optional[np.ndarray] = None,
        target_median_umi: Optional[float] = None,
        precomputed_ctrl_embed: Optional[np.ndarray] = None,
        precomputed_sampling_audit: Optional[Dict[str, Any]] = None,
        precomputed_gene_mapping: Optional[Dict[str, Any]] = None,
        hgnc_bridge: Optional[HGNCBridge] = None,
        pre_sampled: bool = False,
    ) -> Dict[str, Any]:
        """
        Run one arm (matched or mismatched) for `target_line` using controls from
        `control_source_line`, returning pseudobulk predicted deltas `(P, G_model)`.
        """
        if not (is_synthetic_data and self.allow_synthetic_smoke):
            assert_freeze_gates_r16(self.status_dir)
        assert_zero_training_leakage(self.model_key, target_line, self.checkpoint_id)

        self.load_weights()
        hb = get_hgnc_bridge(hgnc_bridge or self.hgnc_bridge) if (
            hgnc_bridge is not None or self.hgnc_bridge is not None or DEFAULT_HGNC_FILE.exists()
        ) else None

        ref_bundle = next(iter(self._bundles.values()))
        model_genes = ref_bundle["gene_names"]

        if precomputed_ctrl_embed is not None:
            emb_arr = np.asarray(precomputed_ctrl_embed, dtype=np.float32)
            if emb_arr.ndim == 2:
                if emb_arr.shape[0] % cell_set_len != 0:
                    raise ValueError(f"precomputed_ctrl_embed rows {emb_arr.shape[0]} not divisible by {cell_set_len}")
                n_sets = emb_arr.shape[0] // cell_set_len
                ctrl_sets_embed = emb_arr.reshape(n_sets, cell_set_len, emb_arr.shape[1])
            elif emb_arr.ndim == 3:
                ctrl_sets_embed = emb_arr
            else:
                raise ValueError(f"Invalid precomputed_ctrl_embed shape {emb_arr.shape}")
            samp = precomputed_sampling_audit or {"n_sampled": int(ctrl_sets_embed.shape[0] * cell_set_len), "n_sets": int(ctrl_sets_embed.shape[0]), "sampling_mode": "precomputed"}
            map_stats = precomputed_gene_mapping or {"source": "precomputed"}
            ctrl_mean_gene = None
        else:
            if ctrl_matrix is None:
                raise ValueError("Provide either `ctrl_matrix` or `precomputed_ctrl_embed`.")
            if pre_sampled:
                if ctrl_matrix.shape[0] % cell_set_len != 0 or ctrl_matrix.shape[0] == 0:
                    raise ValueError(f"pre_sampled ctrl_matrix rows {ctrl_matrix.shape[0]} invalid for cell_set_len={cell_set_len}")
                n_sets_ps = ctrl_matrix.shape[0] // cell_set_len
                samp = precomputed_sampling_audit or {
                    "n_sampled": int(ctrl_matrix.shape[0]),
                    "n_sets": int(n_sets_ps),
                    "sampling_mode": "pre_sampled",
                }
                ctrl_sub = np.asarray(ctrl_matrix, dtype=np.float32)
                umi_sub = np.asarray(cell_total_umis, dtype=np.float64) if cell_total_umis is not None else None
            else:
                samp = sample_stratified_control_sets_r16(
                    n_cells=ctrl_matrix.shape[0],
                    qc_mask=qc_mask,
                    batch_labels=batch_labels,
                    n_ctrl=n_ctrl,
                    cell_set_len=cell_set_len,
                    min_ctrl_floor=MIN_CTRL_FLOOR_DEFAULT,
                    seed=seed,
                )
                if samp["excluded"]:
                    raise ValueError(
                        f"Control source {control_source_line!r} excluded under floor: {samp['exclusion_reason']}"
                    )
                sel_idx = samp["selected_indices"]
                ctrl_sub = np.asarray(ctrl_matrix[sel_idx], dtype=np.float32)
                umi_sub = np.asarray(cell_total_umis[sel_idx], dtype=np.float64) if cell_total_umis is not None else None

            if self.model_key == "STATE_HVG":
                if input_is_log1p:
                    ctrl_log1p = ctrl_sub
                else:
                    ctrl_log1p = preprocess_hvg_expression(
                        ctrl_sub,
                        mode=self.hvg_preproc_mode,
                        cell_total_umis=umi_sub,
                        target_median_umi=target_median_umi,
                    )
                ctrl_aligned, map_stats = align_hvg_expression_to_model(
                    ctrl_log1p,
                    source_genes,
                    model_genes,
                    source_ensembl_ids=source_ensembl_ids,
                    hgnc_bridge=hb,
                )
                ctrl_sets_embed = ctrl_aligned.reshape(samp["n_sets"], cell_set_len, len(model_genes))
                ctrl_mean_gene = np.mean(ctrl_aligned, axis=0, dtype=np.float64)
            else:
                if input_is_log1p and float(np.max(ctrl_sub)) <= RAW_COUNT_HEURISTIC_THRESHOLD:
                    raise ValueError(
                        "STATE_SE requires raw UMI counts (or precomputed SE-600M embeddings), got normalized log1p matrix."
                    )
                enc = self._get_se_encoder()
                ctrl_emb, map_stats = enc.encode_cells(
                    ctrl_counts=ctrl_sub,
                    source_genes=source_genes,
                    source_ensembl_ids=source_ensembl_ids,
                    seed=seed,
                    hgnc_bridge=hb,
                )
                ctrl_sets_embed = ctrl_emb.reshape(samp["n_sets"], cell_set_len, ctrl_emb.shape[1])
                ctrl_mean_gene = None

        per_ckpt_deltas: Dict[str, np.ndarray] = {}
        delta_sum = np.zeros((len(eval_perts), len(model_genes)), dtype=np.float64)
        for ck, bundle in self._bundles.items():
            pred_perts_cells, pred_ntc_cells = self._run_forward_sets_single_ckpt(
                bundle=bundle,
                ctrl_sets_embed=ctrl_sets_embed,
                pert_names=eval_perts,
                hgnc_bridge=hb,
            )
            mean_pert = np.mean(pred_perts_cells, axis=1, dtype=np.float64)
            if self.delta_mode == "minus_pred_ntc":
                baseline = np.mean(pred_ntc_cells, axis=0, dtype=np.float64)
            else:
                if ctrl_mean_gene is None:
                    raise ValueError("delta_mode='minus_input_ctrl' is undefined for latent STATE_SE embeddings.")
                baseline = ctrl_mean_gene
            d_ck = mean_pert - baseline[None, :]
            per_ckpt_deltas[ck] = d_ck
            delta_sum += d_ck

        ensemble_delta = delta_sum / float(len(self._bundles))
        return {
            "model_key": self.model_key,
            "checkpoint_id": self.checkpoint_id,
            "target_line": target_line,
            "control_source_line": control_source_line,
            "delta_mode": self.delta_mode,
            "neutral_batch_mode": self.neutral_batch_mode,
            "hvg_preproc_mode": self.hvg_preproc_mode,
            "model_genes": model_genes,
            "eval_perts": [str(p) for p in eval_perts],
            "delta_pred": ensemble_delta,
            "per_ckpt_deltas": per_ckpt_deltas,
            "ctrl_embed_centroid": np.mean(ctrl_sets_embed.reshape(-1, ctrl_sets_embed.shape[-1]), axis=0),
            "ctrl_sets_embed": ctrl_sets_embed.reshape(-1, ctrl_sets_embed.shape[-1]),
            "sampling_audit": samp,
            "gene_mapping": map_stats,
        }

    def run_matched_and_mismatched_arms(
        self,
        target_line: str,
        matched_ctrl_matrix: Optional[np.ndarray],
        matched_genes: Sequence[str],
        mismatched_ctrls: Dict[str, Tuple[Optional[np.ndarray], Sequence[str], Optional[Sequence[Any]]]],
        eval_perts: Sequence[str],
        matched_batch_labels: Optional[Sequence[Any]] = None,
        matched_qc_mask: Optional[Sequence[bool]] = None,
        input_is_log1p: bool = False,
        seed: int = DEFAULT_R16_SEED,
        n_ctrl: int = N_CTRL_DEFAULT,
        cell_set_len: int = CELL_SET_LEN_DEFAULT,
        is_synthetic_data: bool = False,
        source_total_umis: Optional[Dict[str, Tuple[Optional[np.ndarray], Optional[float]]]] = None,
        precomputed_embeddings: Optional[Dict[str, Tuple[np.ndarray, Dict[str, Any], Dict[str, Any]]]] = None,
        hgnc_bridge: Optional[HGNCBridge] = None,
        pre_sampled: bool = False,
        pre_sampled_audits: Optional[Dict[str, Dict[str, Any]]] = None,
    ) -> Dict[str, Any]:
        """Run matched arm + all stratified mismatched control arms using identical `eval_perts`."""
        umi_map = source_total_umis or {}
        emb_map = precomputed_embeddings or {}
        ps_audits = pre_sampled_audits or {}

        m_umi, m_med = umi_map.get(target_line, (None, None))
        m_pre = emb_map.get(target_line)
        matched_res = self.predict_arm_deltas(
            target_line=target_line,
            control_source_line=target_line,
            ctrl_matrix=matched_ctrl_matrix,
            source_genes=matched_genes,
            eval_perts=eval_perts,
            batch_labels=matched_batch_labels,
            qc_mask=matched_qc_mask,
            input_is_log1p=input_is_log1p,
            seed=seed,
            n_ctrl=n_ctrl,
            cell_set_len=cell_set_len,
            is_synthetic_data=is_synthetic_data,
            cell_total_umis=m_umi,
            target_median_umi=m_med,
            precomputed_ctrl_embed=m_pre[0] if m_pre else None,
            precomputed_sampling_audit=(m_pre[1] if m_pre else ps_audits.get(target_line)),
            precomputed_gene_mapping=m_pre[2] if m_pre else None,
            hgnc_bridge=hgnc_bridge,
            pre_sampled=pre_sampled,
        )

        mismatched_res: Dict[str, Dict[str, Any]] = {}
        centroid_diag: Dict[str, Dict[str, float]] = {}
        c_m = matched_res["ctrl_embed_centroid"].astype(np.float64)
        nrm_m = float(np.linalg.norm(c_m))

        for src_name, (mat_s, genes_s, batch_s) in sorted(mismatched_ctrls.items()):
            s_umi, s_med = umi_map.get(src_name, (None, None))
            s_pre = emb_map.get(src_name)
            res_s = self.predict_arm_deltas(
                target_line=target_line,
                control_source_line=src_name,
                ctrl_matrix=mat_s,
                source_genes=genes_s,
                eval_perts=eval_perts,
                batch_labels=batch_s,
                qc_mask=None,
                input_is_log1p=input_is_log1p,
                seed=seed,
                n_ctrl=n_ctrl,
                cell_set_len=cell_set_len,
                is_synthetic_data=is_synthetic_data,
                cell_total_umis=s_umi,
                target_median_umi=s_med,
                precomputed_ctrl_embed=s_pre[0] if s_pre else None,
                precomputed_sampling_audit=(s_pre[1] if s_pre else ps_audits.get(src_name)),
                precomputed_gene_mapping=s_pre[2] if s_pre else None,
                hgnc_bridge=hgnc_bridge,
                pre_sampled=pre_sampled,
            )
            mismatched_res[src_name] = res_s
            c_s = res_s["ctrl_embed_centroid"].astype(np.float64)
            nrm_s = float(np.linalg.norm(c_s))
            c_mc = c_m - np.mean(c_m)
            c_sc = c_s - np.mean(c_s)
            den = float(np.linalg.norm(c_mc) * np.linalg.norm(c_sc))
            cos_cent = float(np.dot(c_mc, c_sc) / den) if den > 1e-12 else float("nan")
            l2_cent = float(np.linalg.norm(c_m - c_s))
            centroid_diag[src_name] = {
                "centered_cosine": cos_cent,
                "l2_distance": l2_cent,
                "matched_norm": nrm_m,
                "mismatched_norm": nrm_s,
            }

        return {
            "matched": matched_res,
            "mismatched": mismatched_res,
            "control_centroid_diagnostics": centroid_diag,
        }
