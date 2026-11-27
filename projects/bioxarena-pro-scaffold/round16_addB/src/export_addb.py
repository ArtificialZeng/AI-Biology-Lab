#!/usr/bin/env python3
"""
Model-free auditable export of Occam encoder latents (`z_ctrl_<pool>.npz`) and
predicted delta arrays (`preds_<arm>.npz`) with SHA-256 manifests and three-party
comparison utilities (`PREREG_I_ADDENDUM_B_DRAFT_v2.md §6.2`).

Audit tolerances locked in `PREREG_I_ADDENDUM_B_DRAFT_v2.md §6.2`:
  - Layer 1 (control barcodes): exact 100% ordered string equality;
  - Layer 2 (`z_ctrl` encoder latents, `float32`, shape `(N, 512)`): `max_abs <= 1e-5`;
  - Layer 3 (`delta_pred` arrays, `float64`, shape `(P, G)`): `max_abs <= 1e-4`.

Zero Occam imports or vendored Occam code.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence, Set, Tuple

import numpy as np

OCCAM_MODEL_KEY = "Occam"
OCCAM_CHECKPOINT_ID = "rhj_fewshot/all_sources/seed_42"
OCCAM_LATENT_DIM = 512

VALID_INPUT_AXIS_MODES = ("rhj_6983_axis", "pairwise_common_axis", "full_44728_api")
VALID_BASELINE_MODES = ("reconstructed", "real")

Z_CTRL_CROSS_PARTY_ATOL = 1e-5
DELTA_CROSS_PARTY_ATOL = 1e-4


def sanitize_artifact_token(name: str) -> str:
    """Convert an arm or pool identifier into a deterministic filesystem-safe token."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(name).strip()).strip("_")
    if not cleaned or not re.search(r"[A-Za-z0-9]", cleaned):
        raise ValueError(f"Token must contain alphanumeric characters, got {name!r}.")
    return cleaned


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_array_bytes(arr: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(arr).tobytes()).hexdigest()


def save_occam_z_ctrl_latents(
    latents_dir: Path,
    z_ctrl_entries: Mapping[str, Mapping[str, Any]],
) -> Dict[str, str]:
    """
    Export per-pool Occam encoder latents `z_ctrl_<safe_key>.npz` + `z_ctrl_manifest_sha256.json`
    for staged Layer-2 cross-party comparison (`max_abs <= 1e-5`) prior to predicted deltas.
    """
    ldir = Path(latents_dir)
    ldir.mkdir(parents=True, exist_ok=True)
    manifest: Dict[str, str] = {}
    seen_fnames: Set[str] = set()

    for pool_key, entry in sorted(z_ctrl_entries.items()):
        safe_key = sanitize_artifact_token(pool_key)
        z_arr = np.ascontiguousarray(entry["Z_ctrl"], dtype=np.float32)
        if z_arr.ndim != 2 or z_arr.shape[0] <= 0 or z_arr.shape[1] != OCCAM_LATENT_DIM:
            raise ValueError(
                f"z_ctrl {pool_key!r} must have shape (N > 0, {OCCAM_LATENT_DIM}), got {z_arr.shape}."
            )
        if not np.all(np.isfinite(z_arr)):
            raise ValueError(f"z_ctrl {pool_key!r} contains non-finite values.")

        barcodes = [str(b) for b in entry["control_barcodes"]]
        if len(barcodes) != z_arr.shape[0]:
            raise ValueError(
                f"z_ctrl {pool_key!r}: len(control_barcodes)={len(barcodes)} != rows={z_arr.shape[0]}."
            )
        if len(set(barcodes)) != len(barcodes):
            raise ValueError(f"z_ctrl {pool_key!r}: duplicate control_barcodes.")

        seed_raw = entry["seed"]
        if isinstance(seed_raw, bool) or not isinstance(seed_raw, (int, np.integer)):
            raise ValueError(f"z_ctrl {pool_key!r}: integer seed required, got {seed_raw!r}.")

        input_axis_mode = str(entry.get("input_axis_mode", "rhj_6983_axis"))
        if input_axis_mode not in VALID_INPUT_AXIS_MODES:
            raise ValueError(f"Invalid input_axis_mode {input_axis_mode!r} for {pool_key!r}.")

        fname = f"z_ctrl_{safe_key}.npz"
        if fname in seen_fnames:
            raise ValueError(f"Duplicate sanitized latent filename {fname!r} for key {pool_key!r}.")
        seen_fnames.add(fname)

        z_sha = sha256_array_bytes(z_arr)
        audit_obj = dict(entry.get("audit") or {})
        audit_obj["z_ctrl_sha256"] = z_sha
        fpath = ldir / fname
        np.savez_compressed(
            fpath,
            Z_ctrl=z_arr,
            control_barcodes=np.asarray(barcodes, dtype=str),
            seed=np.int64(int(seed_raw)),
            pool_name=np.asarray(str(entry.get("pool_name", pool_key))),
            input_axis_mode=np.asarray(input_axis_mode),
            n_active_rhj_genes=np.int64(int(entry.get("n_active_rhj_genes", 0))),
            axis_mask_sha256=np.asarray(str(entry.get("axis_mask_sha256", ""))),
            z_ctrl_sha256=np.asarray(z_sha),
            audit_json=np.asarray(json.dumps(audit_obj, sort_keys=True)),
        )
        manifest[fname] = sha256_file(fpath)

    manifest_path = ldir / "z_ctrl_manifest_sha256.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def build_occam_delta_payload(
    *,
    delta_pred: np.ndarray,
    target_names: Sequence[str],
    gene_names: Sequence[str],
    control_source: str,
    control_barcodes: Sequence[str],
    input_axis_mode: str = "rhj_6983_axis",
    baseline_mode: str = "reconstructed",
    resolved_pert_keys: Optional[Sequence[str]] = None,
    resolved_pert_ids: Optional[Sequence[int]] = None,
    pair_counterpart: str = "",
    n_active_rhj_genes: int = 6983,
    axis_mask_sha256: str = "",
    model_key: str = OCCAM_MODEL_KEY,
    checkpoint_id: str = OCCAM_CHECKPOINT_ID,
) -> Dict[str, Any]:
    """
    Validate and assemble a standardized Round 16 Addendum B predicted-delta `.npz` payload
    matching Round 16's `preds_*.npz` schema plus Occam-specific provenance keys.
    """
    if input_axis_mode not in VALID_INPUT_AXIS_MODES:
        raise ValueError(f"Invalid input_axis_mode: {input_axis_mode!r}")
    if baseline_mode not in VALID_BASELINE_MODES:
        raise ValueError(f"Invalid baseline_mode: {baseline_mode!r}")

    delta_arr = np.ascontiguousarray(delta_pred, dtype=np.float64)
    if delta_arr.ndim != 2 or delta_arr.shape[0] <= 0 or delta_arr.shape[1] <= 0:
        raise ValueError(f"delta_pred must be non-empty 2D array, got shape {delta_arr.shape}.")
    if not np.all(np.isfinite(delta_arr)):
        raise ValueError("delta_pred contains non-finite values.")

    t_list = [str(x) for x in target_names]
    g_list = [str(x) for x in gene_names]
    bc_list = [str(x) for x in control_barcodes]
    if len(t_list) != delta_arr.shape[0]:
        raise ValueError(f"len(target_names)={len(t_list)} != delta_pred rows={delta_arr.shape[0]}.")
    if len(g_list) != delta_arr.shape[1]:
        raise ValueError(f"len(gene_names)={len(g_list)} != delta_pred cols={delta_arr.shape[1]}.")
    if not bc_list or len(set(bc_list)) != len(bc_list):
        raise ValueError("control_barcodes must be non-empty and unique.")

    rkeys = [str(x) for x in (resolved_pert_keys if resolved_pert_keys is not None else t_list)]
    if len(rkeys) != len(t_list):
        raise ValueError(f"len(resolved_pert_keys)={len(rkeys)} != len(target_names)={len(t_list)}.")
    rids = (
        np.asarray(resolved_pert_ids, dtype=np.int64)
        if resolved_pert_ids is not None
        else np.full(len(t_list), -1, dtype=np.int64)
    )
    if rids.shape != (len(t_list),):
        raise ValueError(f"resolved_pert_ids shape {rids.shape} != ({len(t_list)},).")

    if input_axis_mode == "pairwise_common_axis" and not str(pair_counterpart).strip():
        raise ValueError("pairwise_common_axis payloads MUST specify non-empty pair_counterpart.")

    return {
        "delta_pred": delta_arr,
        "target_names": np.asarray(t_list, dtype=str),
        "gene_names": np.asarray(g_list, dtype=str),
        "resolved_pert_keys": np.asarray(rkeys, dtype=str),
        "resolved_pert_ids": rids,
        "model_key": np.asarray(str(model_key)),
        "checkpoint_id": np.asarray(str(checkpoint_id)),
        "input_axis_mode": np.asarray(str(input_axis_mode)),
        "baseline_mode": np.asarray(str(baseline_mode)),
        "control_source": np.asarray(str(control_source)),
        "control_barcodes": np.asarray(bc_list, dtype=str),
        "pair_counterpart": np.asarray(str(pair_counterpart)),
        "n_active_rhj_genes": np.int64(int(n_active_rhj_genes)),
        "axis_mask_sha256": np.asarray(str(axis_mask_sha256)),
        "delta_sha256": np.asarray(sha256_array_bytes(delta_arr)),
    }


def save_occam_line_deltas(
    arrays_dir: Path,
    line_name: str,
    pred_arms: Mapping[str, Mapping[str, Any]],
) -> Dict[str, str]:
    """
    Save per-arm `preds_<safe_arm>.npz` files under `arrays_dir / line_name` plus
    `array_manifest_sha256.json` (`{filename: sha256}`).
    """
    ldir = Path(arrays_dir) / sanitize_artifact_token(line_name)
    ldir.mkdir(parents=True, exist_ok=True)
    manifest: Dict[str, str] = {}
    seen_fnames: Set[str] = set()

    for arm_name, payload in sorted(pred_arms.items()):
        fname = f"preds_{sanitize_artifact_token(arm_name)}.npz"
        if fname in seen_fnames:
            raise ValueError(f"Duplicate sanitized arm filename {fname!r} for arm {arm_name!r}.")
        seen_fnames.add(fname)
        if "delta_pred" not in payload or "control_barcodes" not in payload:
            raise KeyError(f"Arm {arm_name!r} payload missing required keys ('delta_pred', 'control_barcodes').")
        fpath = ldir / fname
        np.savez_compressed(fpath, **dict(payload))
        manifest[fname] = sha256_file(fpath)

    manifest_path = ldir / "array_manifest_sha256.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def verify_exported_manifest(directory: Path, manifest_filename: str) -> Dict[str, Any]:
    """Verify every file listed in a saved `*_manifest_sha256.json` exists and matches its SHA-256."""
    ddir = Path(directory)
    mpath = ddir / manifest_filename
    manifest = json.loads(mpath.read_text(encoding="utf-8"))
    bad: Dict[str, str] = {}
    for fname, want_sha in sorted(manifest.items()):
        fpath = ddir / fname
        if not fpath.is_file():
            bad[fname] = "MISSING"
            continue
        got_sha = sha256_file(fpath)
        if got_sha != want_sha:
            bad[fname] = f"SHA_MISMATCH({got_sha[:12]}!={want_sha[:12]})"
    return {
        "manifest_path": str(mpath),
        "n_files": len(manifest),
        "n_bad": len(bad),
        "ok": len(bad) == 0,
        "bad_files": bad,
    }


def compare_exported_npz_arrays(
    path_a: Path,
    path_b: Path,
    *,
    array_key: str,
    atol: float,
    check_ordered_keys: Sequence[str] = ("control_barcodes",),
) -> Tuple[bool, Dict[str, Any]]:
    """
    Three-party tolerance checker (`§6.2`):
      - Exact ordered equality on `check_ordered_keys` (e.g., `control_barcodes`, `target_names`, `gene_names`);
      - Exact shape match and `max(abs(A - B)) <= atol` on `array_key`.
    """
    za = np.load(path_a, allow_pickle=False)
    zb = np.load(path_b, allow_pickle=False)

    order_checks: Dict[str, bool] = {}
    for k in check_ordered_keys:
        if k in za.files and k in zb.files:
            order_checks[k] = bool(np.array_equal(za[k], zb[k]))
        else:
            order_checks[k] = False

    a = np.asarray(za[array_key], dtype=np.float64)
    b = np.asarray(zb[array_key], dtype=np.float64)
    shape_ok = a.shape == b.shape
    finite_ok = bool(np.all(np.isfinite(a)) and np.all(np.isfinite(b)))
    max_abs = float(np.max(np.abs(a - b))) if (shape_ok and finite_ok and a.size > 0) else float("inf")
    passed = bool(shape_ok and finite_ok and all(order_checks.values()) and max_abs <= float(atol))
    return passed, {
        "array_key": array_key,
        "shape_a": list(a.shape),
        "shape_b": list(b.shape),
        "shape_ok": shape_ok,
        "finite_ok": finite_ok,
        "order_checks": order_checks,
        "max_abs_diff": max_abs,
        "atol": float(atol),
        "passed": passed,
    }
