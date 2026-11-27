#!/usr/bin/env python3
"""
Model-dependent Occam inference adapter for Round 16 Addendum B
(`PREREG_I_ADDENDUM_B_DRAFT_v2.md §3.3–§3.5, §5.1–§5.2, §6.1–§6.2`).

Strict compliance:
  - Imports Occam's unchanged upstream classes strictly via `sys.path` / `PYTHONPATH`
    from the hash-locked extracted source tree
    `round16_addB/shared/models/occam_src_9ddb890e/src/Occam-9ddb890e65ca1ec2f81882c78311002a75b52f67/src`.
  - ZERO patches to Occam source files; ZERO mocking or stubbing of missing dependencies.
  - ZERO reads of any removed `data/` files or unapproved assets.
  - CPU-first (`device="cpu"`, `torch.set_num_threads(6)` by default; MPS explicitly blocked
    unless `allow_mps=True` after controller confirms `A-T3` finished Round 16).
"""

from __future__ import annotations

import hashlib
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import torch

ADDB_SRC_DIR = Path(__file__).resolve().parent
ADDB_ROOT = ADDB_SRC_DIR.parent
REPO_ROOT = ADDB_ROOT.parent

if str(ADDB_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(ADDB_SRC_DIR))

from export_addb import sha256_file  # noqa: E402
from preprocess_addb import (  # noqa: E402
    OCCAM_N_ENCODER_GENES,
    OCCAM_TARGET_SUM,
    RHJ_GENE_AXIS_SHA256,
    RHJ_N_AXIS_GENES,
    RhjGeneAxis,
    load_rhj_gene_axis,
    load_rhj_perturbations_table,
)

DEFAULT_OCCAM_SRC_DIR = (
    ADDB_ROOT
    / "shared"
    / "models"
    / "occam_src_9ddb890e"
    / "src"
    / "Occam-9ddb890e65ca1ec2f81882c78311002a75b52f67"
    / "src"
)
DEFAULT_OCCAM_HF_DIR = ADDB_ROOT / "shared" / "models" / "occam_hf_acd080b7"
DEFAULT_OCCAM_CKPT_PATH = (
    DEFAULT_OCCAM_HF_DIR
    / "checkpoints"
    / "occam_seeds"
    / "rhj_fewshot"
    / "all_sources"
    / "seed_42"
    / "best_model.pth"
)
DEFAULT_OCCAM_GENE_AXIS_CSV = DEFAULT_OCCAM_HF_DIR / "benchmarks" / "rhj_fewshot" / "gene_axis.csv"
DEFAULT_OCCAM_PERTURBATIONS_CSV = DEFAULT_OCCAM_HF_DIR / "benchmarks" / "rhj_fewshot" / "perturbations.csv"
DEFAULT_OCCAM_SPLITS_MANIFEST_JSON = (
    DEFAULT_OCCAM_HF_DIR / "benchmarks" / "rhj_fewshot" / "splithalf" / "splits_manifest.json"
)
DEFAULT_OCCAM_CHECKSUMS_FILE = ADDB_ROOT / "R16_ADDB_CHECKSUMS.sha256"

LOCKED_BEST_MODEL_SHA256 = "11d5ba5515e79710b3a952532c60aa29515046117292528c59d452a0937685ff"
LOCKED_N_PERT2ID = 18216
LOCKED_FILM_HIDDEN_DIM = 256


def configure_occam_runtime(
    *,
    device: str = "cpu",
    num_threads: int = 6,
    seed: int = 20261601,
    allow_mps: bool = False,
) -> str:
    """
    Lock thread count, deterministic CPU execution, and device guardrails.
    Refuses MPS/CUDA unless `allow_mps=True` (strictly CPU while A-T3 runs R16).
    """
    dev = str(device).strip().lower()
    if dev != "cpu":
        if dev == "mps" and not allow_mps:
            raise RuntimeError(
                "MPS execution is prohibited until Master Controller confirms A-T3 has finished Round 16."
            )
        if dev not in ("mps",):
            raise ValueError(f"Unsupported Occam runtime device {device!r} (expected 'cpu' or gated 'mps').")
    if int(num_threads) > 0:
        torch.set_num_threads(int(num_threads))
    torch.manual_seed(int(seed))
    np.random.seed(int(seed) % (2**32))
    if dev == "cpu":
        try:
            torch.use_deterministic_algorithms(True)
        except Exception:
            pass
    return dev


def ensure_occam_src_on_path(occam_src_dir: Path = DEFAULT_OCCAM_SRC_DIR) -> Path:
    """Prepend the unchanged extracted Occam `src/` tree to `sys.path`."""
    src_path = Path(occam_src_dir).resolve()
    if not (src_path / "occam" / "vcell" / "pipeline.py").is_file():
        raise FileNotFoundError(
            f"Unchanged Occam source tree not found at {src_path} "
            f"(expected occam/vcell/pipeline.py)."
        )
    s = str(src_path)
    if s not in sys.path:
        sys.path.insert(0, s)
    return src_path


def verify_checksums_file_hf_entries(
    checksums_path: Path = DEFAULT_OCCAM_CHECKSUMS_FILE,
    models_root: Path = ADDB_ROOT / "shared" / "models",
) -> Dict[str, Any]:
    """Verify every `occam_hf_acd080b7/...` file hash recorded in `R16_ADDB_CHECKSUMS.sha256`."""
    cpath = Path(checksums_path)
    if not cpath.is_file():
        return {"checksums_file": str(cpath), "present": False, "n_checked": 0, "ok": False, "mismatches": ["MISSING_FILE"]}

    checked: Dict[str, str] = {}
    mismatches: List[str] = []
    for raw_line in cpath.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(None, 1)
        if len(parts) != 2:
            continue
        want_sha, rel = parts[0].strip(), parts[1].strip()
        if not rel.startswith("occam_hf_acd080b7/"):
            continue
        fpath = Path(models_root) / rel
        if not fpath.is_file():
            mismatches.append(f"MISSING:{rel}")
            continue
        got_sha = sha256_file(fpath)
        checked[rel] = got_sha
        if got_sha != want_sha:
            mismatches.append(f"MISMATCH:{rel}:{got_sha[:12]}!={want_sha[:12]}")
    return {
        "checksums_file": str(cpath),
        "present": True,
        "n_checked": len(checked),
        "ok": len(mismatches) == 0 and len(checked) > 0,
        "checked": checked,
        "mismatches": mismatches,
    }


def _strip_encoder_prefix(enc_state_dict: Mapping[str, torch.Tensor]) -> Dict[str, torch.Tensor]:
    out: Dict[str, torch.Tensor] = {}
    for k, v in enc_state_dict.items():
        clean_k = k.replace("module.", "", 1).replace("_orig_mod.", "", 1)
        if clean_k.startswith("encoder."):
            clean_k = clean_k[len("encoder."):]
        out[clean_k] = v
    return out


@dataclass
class Gate5ACheckpointReport:
    """Structured result of Gate 5A Item 1 (checksums, config, & strict=True load assertions)."""

    passed: bool
    ckpt_path: str
    ckpt_sha256: str
    gene_axis_content_sha256: str
    n_gene_names: int
    n_pert2id: int
    target_sum: float
    output_activation: str
    interaction_type: str
    film_mode: str
    film_mean_free: bool
    film_hidden_dim: int
    source_emb_is_none: bool
    modality_emb_is_none: bool
    n_batch_encoder_params: int
    strict_load_ok: bool
    sidecar_checksums_ok: bool
    failures: List[str]
    details: Dict[str, Any]


class OccamInferenceAdapter:
    """
    Clean-room inference adapter wrapping Occam's unchanged `VirtualCellModel`
    loaded via `VirtualCellTrainer.from_vcell_checkpoint`.
    """

    def __init__(
        self,
        model: torch.nn.Module,
        ckpt: Mapping[str, Any],
        rhj_axis: RhjGeneAxis,
        rhj_eligible_keys: Sequence[str],
        *,
        device: str = "cpu",
        batch_size: int = 512,
        gate5a_item1: Optional[Gate5ACheckpointReport] = None,
    ) -> None:
        self.model = model
        self.model.eval()
        self.device = str(device)
        self.batch_size = int(batch_size)
        self.encoder_gene_names: List[str] = [str(g) for g in ckpt["gene_names"]]
        self.gene_to_index: Dict[str, int] = {g: i for i, g in enumerate(self.encoder_gene_names)}
        self.pert2id: Dict[str, int] = {str(k): int(v) for k, v in ckpt["pert2id"].items()}
        self.target_sum: float = float(ckpt["target_sum"])
        self.rhj_axis = rhj_axis
        self.rhj_eligible_keys: List[str] = list(rhj_eligible_keys)
        self.rhj_enc_cols: np.ndarray = np.asarray(
            [self.gene_to_index[g] for g in rhj_axis.ensembl_ids], dtype=np.int64
        )
        self.gate5a_item1 = gate5a_item1

    @classmethod
    def load_verified(
        cls,
        *,
        ckpt_path: Path = DEFAULT_OCCAM_CKPT_PATH,
        occam_src_dir: Path = DEFAULT_OCCAM_SRC_DIR,
        gene_axis_csv: Path = DEFAULT_OCCAM_GENE_AXIS_CSV,
        perturbations_csv: Path = DEFAULT_OCCAM_PERTURBATIONS_CSV,
        checksums_file: Path = DEFAULT_OCCAM_CHECKSUMS_FILE,
        device: str = "cpu",
        num_threads: int = 6,
        allow_mps: bool = False,
        batch_size: int = 512,
        verify_file_hashes: bool = True,
    ) -> Tuple["OccamInferenceAdapter", Gate5ACheckpointReport]:
        """
        Execute Gate 5A Item 1 checksum + structural + `strict=True` load assertions
        and return `(adapter, report)`. Raises `RuntimeError` immediately on any failure.
        """
        dev = configure_occam_runtime(device=device, num_threads=num_threads, allow_mps=allow_mps)
        ensure_occam_src_on_path(occam_src_dir)

        cpath = Path(ckpt_path).resolve()
        failures: List[str] = []

        # 1. File SHA-256 & sidecar checksum verification
        ckpt_sha = sha256_file(cpath) if verify_file_hashes else "SKIPPED"
        if verify_file_hashes and ckpt_sha != LOCKED_BEST_MODEL_SHA256:
            failures.append(f"best_model.pth sha256 {ckpt_sha} != locked {LOCKED_BEST_MODEL_SHA256}")

        sidecar_res = (
            verify_checksums_file_hf_entries(checksums_file)
            if verify_file_hashes
            else {"ok": True, "n_checked": 0, "mismatches": []}
        )
        if verify_file_hashes and not sidecar_res["ok"]:
            failures.append(f"sidecar checksum failures: {sidecar_res['mismatches']}")

        # 2. Gene axis CONTENT hash verification (splits_manifest.json convention)
        rhj_axis = load_rhj_gene_axis(gene_axis_csv, verify_sha256=True)
        ptab = load_rhj_perturbations_table(perturbations_csv)

        # 3. Load Occam's unchanged VirtualCellTrainer.from_vcell_checkpoint
        from occam.vcell.pipeline import VirtualCellTrainer  # imported strictly from extracted tree

        raw_meta = torch.load(str(cpath), map_location="cpu", weights_only=False)
        recorded_act = str((raw_meta.get("decoder_config") or {}).get("output_activation", ""))
        del raw_meta

        trainer, ckpt = VirtualCellTrainer.from_vcell_checkpoint(
            str(cpath), device=dev, output_activation=recorded_act or None
        )
        model = trainer.model
        model.eval()

        # 4. Exact Gate 5A Item 1 metadata & architecture assertions
        is_vc = bool(ckpt.get("is_virtual_cell_checkpoint") is True)
        if not is_vc:
            failures.append(f"is_virtual_cell_checkpoint={ckpt.get('is_virtual_cell_checkpoint')!r} != True")

        tsum = float(ckpt.get("target_sum", -1.0))
        if tsum != OCCAM_TARGET_SUM:
            failures.append(f"target_sum={tsum} != {OCCAM_TARGET_SUM}")

        out_act = str((ckpt.get("decoder_config") or {}).get("output_activation", ""))
        if out_act != "softplus":
            failures.append(f"output_activation={out_act!r} != 'softplus'")

        ic = ckpt.get("interaction_config") or {}
        itype = str(ic.get("interaction_type", ""))
        fmode = str(ic.get("film_mode", ""))
        fmean_free = bool(ic.get("film_mean_free") is True)
        fhidden = int(ic.get("film_hidden_dim", -1))
        if itype != "film":
            failures.append(f"interaction_type={itype!r} != 'film'")
        if fmode != "multiplicative_only":
            failures.append(f"film_mode={fmode!r} != 'multiplicative_only'")
        if not fmean_free:
            failures.append(f"film_mean_free={ic.get('film_mean_free')!r} != True")
        if fhidden != LOCKED_FILM_HIDDEN_DIM:
            failures.append(f"film_hidden_dim={fhidden} != {LOCKED_FILM_HIDDEN_DIM}")

        n_genes = len(ckpt.get("gene_names") or [])
        n_perts = len(ckpt.get("pert2id") or {})
        if n_genes != OCCAM_N_ENCODER_GENES:
            failures.append(f"len(gene_names)={n_genes} != {OCCAM_N_ENCODER_GENES}")
        if n_perts != LOCKED_N_PERT2ID:
            failures.append(f"len(pert2id)={n_perts} != {LOCKED_N_PERT2ID}")

        bench = ckpt.get("benchmark") or {}
        ckpt_axis_sha = str(bench.get("gene_axis_sha256", ""))
        if ckpt_axis_sha != RHJ_GENE_AXIS_SHA256:
            failures.append(f"ckpt benchmark.gene_axis_sha256={ckpt_axis_sha!r} != {RHJ_GENE_AXIS_SHA256!r}")
        if rhj_axis.axis_sha256 != RHJ_GENE_AXIS_SHA256:
            failures.append(f"csv gene_axis_sha256={rhj_axis.axis_sha256!r} != {RHJ_GENE_AXIS_SHA256!r}")
        if len(rhj_axis.ensembl_ids) != RHJ_N_AXIS_GENES:
            failures.append(f"len(rhj_axis)={len(rhj_axis.ensembl_ids)} != {RHJ_N_AXIS_GENES}")

        src_emb_none = getattr(model, "source_emb", "MISSING") is None
        mod_emb_none = getattr(model, "modality_emb", "MISSING") is None
        if not src_emb_none:
            failures.append("model.source_emb is not None")
        if not mod_emb_none:
            failures.append("model.modality_emb is not None")

        batch_enc_params = [n for n, _ in model.named_parameters() if "batch_encoder" in n.lower()]
        batch_enc_ckpt_keys = [
            f"{sd_name}:{k}"
            for sd_name in ("encoder_state_dict", "gnn_state_dict", "decoder_state_dict", "interaction_state_dict")
            for k in (ckpt.get(sd_name) or {}).keys()
            if "batch_encoder" in str(k).lower()
        ]
        if batch_enc_params or batch_enc_ckpt_keys:
            failures.append(f"Unexpected batch_encoder entries: params={batch_enc_params}, ckpt={batch_enc_ckpt_keys}")

        # 5. Explicit strict=True load verification across all 4 submodules
        strict_details: Dict[str, Dict[str, int]] = {}
        strict_ok = True
        try:
            enc_sd = _strip_encoder_prefix(ckpt["encoder_state_dict"])
            for sub_name, sub_mod, sd in (
                ("encoder", model.encoder, enc_sd),
                ("gnn", model.gnn, ckpt["gnn_state_dict"]),
                ("decoder", model.decoder, ckpt["decoder_state_dict"]),
                ("interaction", model.interaction, ckpt["interaction_state_dict"]),
            ):
                m_keys, u_keys = sub_mod.load_state_dict(sd, strict=True)
                strict_details[sub_name] = {
                    "n_tensors": len(sd),
                    "missing": len(m_keys),
                    "unexpected": len(u_keys),
                }
                if m_keys or u_keys:
                    strict_ok = False
                    failures.append(f"{sub_name} strict=True returned missing={m_keys}, unexpected={u_keys}")
        except Exception as exc:
            strict_ok = False
            failures.append(f"strict=True load failed: {type(exc).__name__}: {exc}")

        model.invalidate_gnn_cache()
        model.eval()

        report = Gate5ACheckpointReport(
            passed=(len(failures) == 0),
            ckpt_path=str(cpath),
            ckpt_sha256=ckpt_sha,
            gene_axis_content_sha256=rhj_axis.axis_sha256,
            n_gene_names=n_genes,
            n_pert2id=n_perts,
            target_sum=tsum,
            output_activation=out_act,
            interaction_type=itype,
            film_mode=fmode,
            film_mean_free=fmean_free,
            film_hidden_dim=fhidden,
            source_emb_is_none=src_emb_none,
            modality_emb_is_none=mod_emb_none,
            n_batch_encoder_params=len(batch_enc_params),
            strict_load_ok=strict_ok,
            sidecar_checksums_ok=bool(sidecar_res["ok"]),
            failures=failures,
            details={
                "epoch": int(ckpt.get("epoch", -1)),
                "best_metric": float(ckpt.get("best_metric", float("nan"))),
                "benchmark_id": str(bench.get("id", "")),
                "axis_scope": str(bench.get("axis_scope", "")),
                "strict_submodules": strict_details,
                "sidecar_n_checked": int(sidecar_res["n_checked"]),
            },
        )
        if not report.passed:
            raise RuntimeError("Gate 5A Item 1 FAILED: " + "; ".join(failures))

        adapter = cls(
            model=model,
            ckpt=ckpt,
            rhj_axis=rhj_axis,
            rhj_eligible_keys=ptab["eligible_keys"],
            device=dev,
            batch_size=batch_size,
            gate5a_item1=report,
        )
        return adapter, report

    def encode_controls(self, x_enc: np.ndarray, *, batch_size: Optional[int] = None) -> np.ndarray:
        """Encode `(n_cells, 44728)` preprocessed controls into `(n_cells, 512)` `float32` latents."""
        x = np.ascontiguousarray(x_enc, dtype=np.float32)
        if x.ndim != 2 or x.shape[1] != len(self.encoder_gene_names):
            raise ValueError(f"Expected shape (N, {len(self.encoder_gene_names)}), got {x.shape}.")
        bs = int(batch_size or self.batch_size)
        self.model.eval()
        parts: List[np.ndarray] = []
        with torch.no_grad():
            for s in range(0, x.shape[0], bs):
                xb = torch.from_numpy(x[s : s + bs]).to(self.device)
                zb = self.model.encoder(xb).float().cpu().numpy()
                parts.append(zb)
        z = np.ascontiguousarray(np.concatenate(parts, axis=0), dtype=np.float32)
        if not np.all(np.isfinite(z)):
            raise ValueError("Encoded z_ctrl contains non-finite values.")
        return z

    def encode_rhj_axis_matrix(self, x_rhj: np.ndarray, *, batch_size: Optional[int] = None) -> np.ndarray:
        """Scatter `(n_cells, 6983)` normalized RHJ-axis rows into `(B, 44728)` chunks and encode."""
        x = np.ascontiguousarray(x_rhj, dtype=np.float32)
        if x.ndim != 2 or x.shape[1] != len(self.rhj_enc_cols):
            raise ValueError(f"Expected shape (N, {len(self.rhj_enc_cols)}), got {x.shape}.")
        bs = int(batch_size or self.batch_size)
        n_enc = len(self.encoder_gene_names)
        buf = np.zeros((bs, n_enc), dtype=np.float32)
        self.model.eval()
        parts: List[np.ndarray] = []
        with torch.no_grad():
            for s in range(0, x.shape[0], bs):
                e = min(s + bs, x.shape[0])
                m = e - s
                buf[:m].fill(0.0)
                buf[:m, self.rhj_enc_cols] = x[s:e]
                xb = torch.from_numpy(buf[:m]).to(self.device)
                zb = self.model.encoder(xb).float().cpu().numpy()
                parts.append(zb)
        z = np.ascontiguousarray(np.concatenate(parts, axis=0), dtype=np.float32)
        if not np.all(np.isfinite(z)):
            raise ValueError("Encoded RHJ z_ctrl contains non-finite values.")
        return z

    def predict_zero_pert_from_z(
        self,
        z_ctrl: np.ndarray,
        output_enc_cols: np.ndarray,
        *,
        batch_size: Optional[int] = None,
    ) -> np.ndarray:
        """
        Compute exact unperturbed zero-perturbation forward `D(z + interaction(z, 0))`
        (`pert_indices = [[-1]] * B`) restricted to `output_enc_cols`.
        """
        z = np.ascontiguousarray(z_ctrl, dtype=np.float32)
        cols_t = torch.as_tensor(np.asarray(output_enc_cols, dtype=np.int64), dtype=torch.long)
        bs = int(batch_size or self.batch_size)
        self.model.eval()
        parts: List[np.ndarray] = []
        with torch.no_grad():
            for s in range(0, z.shape[0], bs):
                zb = torch.from_numpy(z[s : s + bs]).to(self.device)
                out = self.model.forward_composed(zb, [[-1]] * zb.shape[0], "sum")["prediction"]
                parts.append(out[:, cols_t.to(out.device)].float().cpu().numpy())
        arr = np.ascontiguousarray(np.concatenate(parts, axis=0), dtype=np.float32)
        if not np.all(np.isfinite(arr)):
            raise ValueError("Zero-perturbation reconstruction contains non-finite values.")
        return arr

    def predict_single_pert_mean_from_z(
        self,
        z_ctrl: np.ndarray,
        pert_id: int,
        output_enc_cols: np.ndarray,
        *,
        compose_mode: str = "decode_seq_mean",
        batch_size: Optional[int] = None,
    ) -> np.ndarray:
        """Return `(len(output_enc_cols),)` `float64` mean predicted profile across `z_ctrl`."""
        z = np.ascontiguousarray(z_ctrl, dtype=np.float32)
        cols_t = torch.as_tensor(np.asarray(output_enc_cols, dtype=np.int64), dtype=torch.long)
        bs = int(batch_size or self.batch_size)
        pid = int(pert_id)
        acc = np.zeros(len(output_enc_cols), dtype=np.float64)
        self.model.eval()
        with torch.no_grad():
            for s in range(0, z.shape[0], bs):
                zb = torch.from_numpy(z[s : s + bs]).to(self.device)
                out = self.model.forward_composed(zb, [[pid]] * zb.shape[0], compose_mode)["prediction"]
                sub = out[:, cols_t.to(out.device)].float().cpu().numpy()
                if not np.all(np.isfinite(sub)):
                    raise ValueError(f"Non-finite prediction for pert_id={pid} at rows {s}:{s + zb.shape[0]}.")
                acc += sub.sum(axis=0, dtype=np.float64)
        return acc / float(z.shape[0])

    def predict_deltas_from_controls(
        self,
        x_enc: np.ndarray,
        pert_ids: Sequence[int],
        output_enc_cols: np.ndarray,
        *,
        z_ctrl: Optional[np.ndarray] = None,
        compose_mode: str = "decode_seq_mean",
        batch_size: Optional[int] = None,
    ) -> Dict[str, Any]:
        """
        Compute BOTH locked delta conventions (`baseline="reconstructed"` primary AND
        `baseline="real"` sensitivity) on exact `output_enc_cols`, returning `(z_ctrl, deltas)`.
        """
        z = self.encode_controls(x_enc, batch_size=batch_size) if z_ctrl is None else np.ascontiguousarray(z_ctrl, dtype=np.float32)
        out_cols = np.asarray(output_enc_cols, dtype=np.int64)
        recon_mean = self.predict_zero_pert_from_z(z, out_cols, batch_size=batch_size).mean(axis=0, dtype=np.float64)
        real_mean = np.asarray(x_enc[:, out_cols], dtype=np.float32).mean(axis=0, dtype=np.float64)

        n_p = len(pert_ids)
        n_g = len(out_cols)
        pred_means = np.empty((n_p, n_g), dtype=np.float64)
        for i, pid in enumerate(pert_ids):
            pred_means[i] = self.predict_single_pert_mean_from_z(
                z, int(pid), out_cols, compose_mode=compose_mode, batch_size=batch_size
            )
        return {
            "z_ctrl": z,
            "pred_means": pred_means,
            "recon_ctrl_mean": recon_mean,
            "real_ctrl_mean": real_mean,
            "delta_reconstructed": pred_means - recon_mean[None, :],
            "delta_real": pred_means - real_mean[None, :],
        }

    def predict_cell_pairs_from_z(
        self,
        z_ctrl_pool: np.ndarray,
        ctrl_indices: np.ndarray,
        pert_ids: Sequence[int],
        output_enc_cols: np.ndarray,
        *,
        compose_mode: str = "decode_seq_mean",
        batch_size: Optional[int] = None,
        progress_label: str = "",
    ) -> np.ndarray:
        """
        Predict one `(len(output_enc_cols),)` `float32` cell per `(ctrl_indices[i], pert_ids[i])`
        pair (used for exact Gate 5B split-half `predict_panel` reproduction).
        """
        n = len(ctrl_indices)
        if len(pert_ids) != n:
            raise ValueError(f"len(pert_ids)={len(pert_ids)} != len(ctrl_indices)={n}.")
        c_idx = np.asarray(ctrl_indices, dtype=np.int64)
        out_cols = np.asarray(output_enc_cols, dtype=np.int64)
        cols_t = torch.as_tensor(out_cols, dtype=torch.long)
        bs = int(batch_size or self.batch_size)

        out_mat = np.empty((n, len(out_cols)), dtype=np.float32)
        self.model.eval()
        with torch.no_grad():
            for s in range(0, n, bs):
                e = min(s + bs, n)
                zb = torch.from_numpy(np.ascontiguousarray(z_ctrl_pool[c_idx[s:e]], dtype=np.float32)).to(self.device)
                pidx = [[int(p)] for p in pert_ids[s:e]]
                pred = self.model.forward_composed(zb, pidx, compose_mode)["prediction"]
                sub = pred[:, cols_t.to(pred.device)].float().cpu().numpy()
                if not np.all(np.isfinite(sub)):
                    raise ValueError(f"[{progress_label}] Non-finite prediction at rows {s}:{e}.")
                out_mat[s:e] = sub
        return out_mat


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()
