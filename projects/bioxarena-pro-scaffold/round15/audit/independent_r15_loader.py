#!/usr/bin/env python3
"""
Independent Round 15 VCC 2026 Reader, Pre-Release Gate, & Feasibility Audit (A-T3).

Strictly implements locked `PREREG_H.md` (v2 + Addendum A, addendum to `PREREG_E.md` commit `4c606f4`):
  - Clean-room isolation: NEVER reads, imports, or references `round15/src/`,
    `round15/results/`, or `round15/status/A-T2.txt`.
  - Candidate lines (`max K_conf = 5`): `cal33`, `hela`, `ac16`, `molm13`, `u87mg`.
    **`jurkat` is STRICTLY EXCLUDED** (homologous to `r9_6` training line `nadig_jurkat`).
  - Line-identity resolution (`§3.1` + `Addendum A`):
      * `(a)` Explicit `.obs` column (`cell_line`, `cell_type`, `line`), else
      * `(b)` `.obs['context']` via official `context_to_line_map`, else
      * `(c)` Official release manifest `line_manifest` (`file -> {line, url, quote}`),
        requiring documented URL + exact statement/entry (filename alone is forbidden).
        If `(a)/(b)` is present and disagrees with `(c)` for any cell, all affected lines
        are marked `EXCLUDED_FEASIBILITY`. `jurkat` rows are dropped before any QC/metric.
  - Access, License & Sunset Locks (`§1, §2`):
      * `synthetic_test`: synthetic `.h5ad` fixtures ONLY; strictly blocks `round15/shared/`.
      * `dev_r9..dev_r14`: historical unblinded caches only; blocks `round15/shared/`.
      * `reader_xcheck_r15`: requires `PREREG_H_COMMIT.txt`, `R15_LINES_AUTHORIZED.txt`, AND
        verified open-license (`CC0`/`CC BY`) status before reading `round15/shared/`.
      * `holdout_r15`: requires `PREREG_H_COMMIT.txt`, `CODE_FREEZE_R15.txt`, AND
        `CONFIRMATORY_LINES_R15.txt` (sunset deadline `2027-06-30`).
  - Pre-Metric HDF5 Integer Slot Verification & Multi-Slot Conflict Stop Guard (`§5`):
      * Reuses frozen 100%-entry scanner (`_full_slot_integer_check`, `_compare_two_integer_slots`).
      * Stops with `EXCLUDED_FEASIBILITY` if zero integer raw UMI slots exist OR if `>= 2`
        non-excluded integer slots disagree in shape, `/var` axis, or matrix values.
  - Dual `.obs` Schema Support (`2025` `target_gene`/`guide_id`/`batch` & `2026`
    `target_gene`/`context`/`ntc_id` with/without `batch`):
      * `NTC`: strictly `target_gene == "non-targeting"`; rejects safe-targeting/gene_desert/multiplet.
      * Single-construct: rejects heterogeneous multi-target delimiters (`[,;|+whitespace]`) & `n_guides > 1`,
        while permitting valid single-construct cassette IDs (`A|A`, `non-targeting_A|non-targeting_B`)
        when `target_gene` resolves cleanly to one perturbation; `PREREG_F §3.2` fallback if absent.
      * Target mapping: `HGNCLookup.resolve(target_gene, "")` with `pert_key = canon_sym`.
      * Cell QC: `n_umis >= 300`, `n_genes >= 200`, `max_pct_mt = None`.
  - Two-Tier Split-Half & NTC Subtraction (`§6`):
      * Tier 1: `.obs` batch/channel parity (`sorted` labels, odd -> Half A, even -> Half B, unit
        subtraction with `min_batch_nt = 15` half-pooled fallback) when `>= 2` units & `>= 15` NTC/half.
      * Tier 2: stratified barcode-sorted `np.random.default_rng(20261231)` 50/50 split (`batches = half_1/half_2`).
"""

from __future__ import annotations

import datetime as _dt
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Set, Tuple

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROUND15_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = ROUND15_ROOT.parent
R10_SCRIPTS_DIR = REPO_ROOT / "round10" / "shared" / "scripts"
R11_AUDIT_DIR = REPO_ROOT / "round11" / "audit"
R11_SRC_DIR = REPO_ROOT / "round11" / "src"
R12_AUDIT_DIR = REPO_ROOT / "round12" / "audit"
R13_AUDIT_DIR = REPO_ROOT / "round13" / "audit"
R14_AUDIT_DIR = REPO_ROOT / "round14" / "audit"

for _p in (
    str(R10_SCRIPTS_DIR),
    str(R11_AUDIT_DIR),
    str(R11_SRC_DIR),
    str(R12_AUDIT_DIR),
    str(R13_AUDIT_DIR),
    str(R14_AUDIT_DIR),
):
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
from independent_r13_loader import R13_LINES_DIR, R13_TARGET_KEYS  # noqa: E402
from independent_r14_loader import (  # noqa: E402
    R14_LINES_DIR,
    R14_TARGET_KEYS,
    _compare_two_integer_slots,
    _decode_str,
    _extract_var_metadata,
    _full_slot_integer_check,
    _is_excluded_layer_name,
    _list_group_columns,
    _read_anndata_index,
    _read_obs_or_var_column,
    make_slot_subset_readers,
)
from prepare_r10_lines import (  # noqa: E402
    HGNCLookup,
    PartSpec,
    build_hgnc_lookup,
    extract_top2_csc,
    fit_poisson_gaussian_thresholds_from_hist,
    process_r10_line,
)

R15_SHARED_DIR = ROUND15_ROOT / "shared"
R15_RAW_DIR = R15_SHARED_DIR / "raw" / "vcc2026"
R15_LINES_DIR = R15_SHARED_DIR / "lines"
R15_STATUS_DIR = ROUND15_ROOT / "status"

DEFAULT_R15_SEED = 20261231
PREREG_H_SUNSET_DATE = _dt.date(2027, 6, 30)

R15_CANDIDATE_LINES: Tuple[str, ...] = ("cal33", "hela", "ac16", "molm13", "u87mg")
R15_CANDIDATE_SET: Set[str] = set(R15_CANDIDATE_LINES)
R15_EXCLUDED_LINES: Set[str] = {"jurkat"}

R10_ALL_LINES: Tuple[str, ...] = ("dld1", "cd4t_rest", "wtc11_i3neuron", "calu3")
R11_ALL_LINES: Tuple[str, ...] = ("telohaec", "cd4t_stim48", "feng_ipsc")

DEV_ALLOWED_LINES_BY_PHASE: Dict[str, Tuple[str, ...]] = {
    "dev_r9": R9_LINE_NAMES,
    "dev_r10": R10_ALL_LINES,
    "dev_r11": R11_ALL_LINES,
    "dev_r12": R12_TARGET_KEYS,
    "dev_r13": R13_TARGET_KEYS,
    "dev_r14": R14_TARGET_KEYS,
}
ALL_DEV_LINES: Set[str] = set().union(*DEV_ALLOWED_LINES_BY_PHASE.values())

_MULTI_TARGET_DELIM_RE = re.compile(r"[,;|+\s]")
_BAD_TARGET_TOKENS: Set[str] = {
    "",
    "na",
    "nan",
    "none",
    "null",
    "unassigned",
    "multiplet",
    "doublet",
    "ambiguous",
    "negative",
    "empty",
    "undetected",
    "safe-targeting",
    "safe_targeting",
    "safetargeting",
    "gene_desert",
    "gene-desert",
}
_VALID_OPEN_LICENSES: Set[str] = {"CC0", "CC0-1.0", "CC0_1.0", "CC-BY", "CC-BY-4.0", "CC_BY", "CC_BY_4.0"}


def _normalize_line_token(raw_val: str) -> str:
    return re.sub(r"[-_\s]+", "", str(raw_val).strip().lower())


ALL_VCC_KNOWN_LINES: Set[str] = R15_CANDIDATE_SET | R15_EXCLUDED_LINES


def _assert_not_inside_r15_shared(path: Path, label: str) -> None:
    r15_shared = R15_SHARED_DIR.resolve()
    resolved = path.resolve()
    if resolved == r15_shared or r15_shared in resolved.parents:
        raise RuntimeError(f"PREREG_H blinding lock ({label}): refused accessing {resolved} inside {r15_shared}.")


def verify_r15_open_release_and_sunset(
    license_id: Optional[str] = None,
    requires_auth_or_terms: bool = False,
    check_date: Optional[_dt.date] = None,
    open_released: bool = False,
) -> Dict[str, Any]:
    """Enforce `PREREG_H §2` open-license (`CC0`/`CC BY`, no login/terms) and `2027-06-30` sunset rule."""
    today = check_date or _dt.date.today()
    if requires_auth_or_terms:
        return {"status": "EXCLUDED_ACCESS", "allowed": False, "reason": "requires_account_terms_or_verification"}
    if not open_released:
        if today > PREREG_H_SUNSET_DATE:
            return {
                "status": "PREREG_H_EXPIRED",
                "allowed": False,
                "reason": f"past_sunset_{PREREG_H_SUNSET_DATE.isoformat()}_no_open_release",
            }
        return {"status": "WAITING_OPEN_RELEASE", "allowed": False, "reason": "public_bucket_empty_or_unreleased"}
    norm_lic = re.sub(r"\s+", "-", str(license_id or "").strip().upper())
    if norm_lic not in _VALID_OPEN_LICENSES and not (norm_lic.startswith("CC0") or norm_lic.startswith("CC-BY") or norm_lic.startswith("CC_BY")):
        return {"status": "EXCLUDED_ACCESS", "allowed": False, "reason": f"non_open_license:{license_id!r}"}
    return {"status": "OPEN_LICENSE_OK", "allowed": True, "license": norm_lic}


def enforce_r15_gate(
    phase: str,
    target_lines: Sequence[str] = (),
    status_dir: Path = R15_STATUS_DIR,
    extra_dirs: Sequence[Path] = (),
    require_confirmatory_file: bool = True,
) -> Dict[str, str]:
    """
    Enforce `PREREG_H §1–§3, §11` blinding, open-release, and `Jurkat` exclusion gates:
      - Unconditionally rejects `"jurkat"` (`nadig_jurkat` is in `r9_6`).
      - `synthetic_test`: synthetic fixtures ONLY; strictly blocks `round15/shared/`.
      - `dev_r9..dev_r14`: historical unblinded caches only; blocks `round15/shared/` and VCC 2026 keys.
      - `reader_xcheck_r15`: requires non-empty `PREREG_H_COMMIT.txt` AND `R15_LINES_AUTHORIZED.txt`.
      - `holdout_r15`: requires `PREREG_H_COMMIT.txt`, `CODE_FREEZE_R15.txt`, AND `CONFIRMATORY_LINES_R15.txt`.
    """
    norm_tgts = [_normalize_line_token(t) for t in target_lines]
    if any(t in R15_EXCLUDED_LINES for t in norm_tgts):
        raise RuntimeError(f"PREREG_H §3 hard exclusion: 'jurkat' is forbidden (in r9_6 source nadig_jurkat): {target_lines}")

    if phase == "synthetic_test":
        for ed in extra_dirs:
            _assert_not_inside_r15_shared(ed, phase)
        return {"phase": phase, "mode": "SYNTHETIC_FIXTURE_ONLY"}

    if phase in DEV_ALLOWED_LINES_BY_PHASE:
        bad_targets = [t for t in target_lines if (_normalize_line_token(t) in R15_CANDIDATE_SET or t not in ALL_DEV_LINES)]
        if bad_targets:
            raise RuntimeError(f"PREREG_H dev blinding lock ({phase}): refused non-dev/R15 targets {bad_targets}.")
        for ed in extra_dirs:
            _assert_not_inside_r15_shared(ed, phase)
        return {"phase": phase, "mode": "DEV_HISTORICAL_UNBLINDED"}

    prereg_hash = _first_token(status_dir / "PREREG_H_COMMIT.txt")
    if not prereg_hash:
        raise RuntimeError(
            f"PREREG_H_COMMIT hard gate blocked ({phase}): missing or empty {status_dir / 'PREREG_H_COMMIT.txt'}."
        )

    if phase == "reader_xcheck_r15":
        auth_tok = _first_token(status_dir / "R15_LINES_AUTHORIZED.txt")
        if not auth_tok:
            raise RuntimeError(
                f"R15_LINES_AUTHORIZED hard gate blocked ({phase}): missing or empty "
                f"{status_dir / 'R15_LINES_AUTHORIZED.txt'}."
            )
        return {
            "phase": phase,
            "prereg_commit": prereg_hash,
            "r15_lines_authorized": auth_tok,
            "code_freeze": _first_token(status_dir / "CODE_FREEZE_R15.txt"),
        }

    if phase == "holdout_r15":
        freeze_hash = _first_token(status_dir / "CODE_FREEZE_R15.txt")
        if not freeze_hash:
            raise RuntimeError(
                f"CODE_FREEZE_R15 hard gate blocked: missing or empty {status_dir / 'CODE_FREEZE_R15.txt'}."
            )
        conf_p = status_dir / "CONFIRMATORY_LINES_R15.txt"
        if require_confirmatory_file and (not conf_p.exists() or not conf_p.read_text(encoding="utf-8").strip()):
            raise RuntimeError(f"CONFIRMATORY_LINES_R15 gate blocked: missing or empty {conf_p}.")
        return {
            "phase": phase,
            "prereg_commit": prereg_hash,
            "code_freeze": freeze_hash,
            "r15_lines_authorized": _first_token(status_dir / "R15_LINES_AUTHORIZED.txt"),
        }

    raise ValueError(f"Unsupported Round 15 phase: {phase!r}")


def parse_confirmatory_lines_r15(conf_path: Path) -> Dict[str, Any]:
    """
    Parse `CONFIRMATORY_LINES_R15.txt` deterministically (`PREREG_H §7–§8`, `max K_conf = 5`).
    Supports `cal33`, `hela`, `ac16`, `molm13`, `u87mg`; strictly rejects `jurkat` as confirmatory.
    """
    if not conf_path.exists():
        raise FileNotFoundError(f"Missing confirmatory gate file: {conf_path}")
    text = conf_path.read_text(encoding="utf-8")
    if not text.strip():
        raise RuntimeError(f"Empty confirmatory gate file: {conf_path}")

    role_by_line: Dict[str, str] = {}
    listed_conf: List[str] = []
    listed_sec: List[str] = []

    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        upper = line.upper()
        if upper.startswith("CONFIRMATORY_LINES"):
            rhs = line.split("=", 1)[-1] if "=" in line else line.split(":", 1)[-1]
            for tok in rhs.replace(",", " ").split():
                t = _normalize_line_token(tok)
                if t in R15_EXCLUDED_LINES:
                    raise RuntimeError("PREREG_H §3 violation: jurkat cannot appear in CONFIRMATORY_LINES.")
                if t in R15_CANDIDATE_SET and t not in listed_conf:
                    listed_conf.append(t)
        elif upper.startswith("SECONDARY_LINES"):
            rhs = line.split("=", 1)[-1] if "=" in line else line.split(":", 1)[-1]
            for tok in rhs.replace(",", " ").split():
                t = _normalize_line_token(tok)
                if t in R15_CANDIDATE_SET and t not in listed_sec:
                    listed_sec.append(t)
        else:
            for cand in R15_CANDIDATE_LINES:
                if upper.startswith(f"{cand.upper()}_ROLE") or _normalize_line_token(line.split("=")[0]) == cand:
                    rhs = line.split("=", 1)[-1] if "=" in line else line.split(":", 1)[-1]
                    role_by_line[cand] = rhs.strip().split()[0].upper()
                    break
            else:
                nt = _normalize_line_token(line)
                if nt in R15_CANDIDATE_SET and nt not in listed_conf:
                    listed_conf.append(nt)

    confirmatory_lines: List[str] = []
    secondary_eval_lines: List[str] = []
    excluded_lines: List[str] = ["jurkat"]

    for cand in R15_CANDIDATE_LINES:
        role = role_by_line.get(cand)
        if role is None:
            if cand in listed_conf:
                role = "CONFIRMATORY"
            elif cand in listed_sec:
                role = "SECONDARY_ONLY"
            else:
                role = "NOT_RELEASED"
        role_by_line[cand] = role
        if role == "CONFIRMATORY":
            confirmatory_lines.append(cand)
        elif role == "SECONDARY_ONLY":
            secondary_eval_lines.append(cand)
        elif role.startswith("EXCLUDED"):
            excluded_lines.append(cand)

    return {
        "confirmatory_lines": confirmatory_lines,
        "k_conf": len(confirmatory_lines),
        "secondary_eval_lines": secondary_eval_lines,
        "excluded_lines": excluded_lines,
        "role_by_line": role_by_line,
        "aborted_confirmatory": len(confirmatory_lines) == 0,
    }


def _extract_vcc_var_metadata(var: h5py.Group, group_label: str) -> Dict[str, Any]:
    """Support `var['_index']` gene names (2025 style) AND `gene_symbol`/`gene_names`/`feature_name`."""
    vm = _extract_var_metadata(var, group_label)
    if vm["var_syms"] == vm["var_index"]:
        for alt_col in ("gene_names", "gene_short_name", "name"):
            tmp = _read_obs_or_var_column(var, alt_col)
            if tmp is not None:
                vm["var_syms"] = [str(x).strip() for x in tmp]
                break
    return vm


def _inspect_vcc_obs_columns(obs: h5py.Group) -> Dict[str, Any]:
    obs_cols = _list_group_columns(obs)
    low_map = {c.lower(): c for c in obs_cols}

    target_col = next((low_map[k] for k in ("target_gene", "perturbation", "gene", "target") if k in low_map), None)
    guide_col = next((low_map[k] for k in ("guide_id", "sgrna_id", "sgrna", "guide", "construct_id") if k in low_map), None)
    ntc_id_col = next((low_map[k] for k in ("ntc_id", "control_id") if k in low_map), None)
    explicit_line_col = next((low_map[k] for k in ("cell_line", "cell_type", "line") if k in low_map), None)
    context_col = low_map.get("context")
    n_guides_col = next((low_map[k] for k in ("n_guides", "nperts", "num_guides", "n_sgrna") if k in low_map), None)

    batch_col = None
    for bc in ("batch", "flex_channel", "channel", "lane", "replicate", "rep", "sample"):
        if bc in low_map:
            batch_col = low_map[bc]
            break

    return {
        "obs_cols": obs_cols,
        "target_col": target_col,
        "guide_col": guide_col,
        "ntc_id_col": ntc_id_col,
        "explicit_line_col": explicit_line_col,
        "context_col": context_col,
        "n_guides_col": n_guides_col,
        "batch_col": batch_col,
    }


def _resolve_manifest_entry_for_path(
    hp: Path,
    line_manifest: Optional[Mapping[str, Any]],
) -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
    """
    Validate and resolve an official rule `(c)` single-line file entry from `line_manifest`.

    Each entry MUST be a mapping specifying:
      - `line` (or `cell_line`) in `cal33, hela, ac16, molm13, u87mg, jurkat`
      - `url` (or `source_url` / `manifest_url` / `documentation_url`) with URL scheme (`https://`, `gs://`, etc.)
      - `quote` (or `sentence` / `statement` / `manifest_entry` / `verbatim`) — exact documentation statement
    A filename alone or bare string without provenance URL + exact statement is strictly rejected.
    """
    if not line_manifest:
        return None, None

    hp_resolved = str(hp.resolve())
    hp_str = str(hp)
    hp_name = hp.name
    full_hits = []
    base_hits = []
    for k, v in line_manifest.items():
        ks = str(k).strip()
        if ks in (hp_resolved, hp_str):
            full_hits.append((ks, v))
        elif ks == hp_name:
            base_hits.append((ks, v))

    all_hits = full_hits + base_hits
    if not all_hits:
        return None, None

    parsed_hits: List[Dict[str, Any]] = []
    for ks, raw_entry in all_hits:
        if not isinstance(raw_entry, Mapping):
            return None, (
                f"invalid_line_manifest_entry_for_{ks!r}: expected dict with line, url, and exact quote "
                "(filename/bare string alone is forbidden by PREREG_H Addendum A)"
            )
        raw_line = raw_entry.get("line") if raw_entry.get("line") is not None else raw_entry.get("cell_line")
        norm_line = _normalize_line_token(str(raw_line or ""))
        if not raw_line or norm_line not in ALL_VCC_KNOWN_LINES:
            return None, f"invalid_line_in_line_manifest_for_{ks!r}:{raw_line!r}"

        url_val = ""
        for uk in ("url", "source_url", "manifest_url", "documentation_url"):
            if raw_entry.get(uk) is not None and str(raw_entry.get(uk)).strip():
                url_val = str(raw_entry.get(uk)).strip()
                break
        if not url_val or ("://" not in url_val and not url_val.lower().startswith("doi:")):
            return None, (
                f"missing_or_invalid_source_url_in_line_manifest_for_{ks!r}: {url_val!r} "
                "(PREREG_H Addendum A requires valid documentation URL e.g. https:// or gs://)"
            )

        quote_val = ""
        for qk in ("quote", "sentence", "statement", "manifest_entry", "verbatim", "exact_sentence"):
            if raw_entry.get(qk) is not None and str(raw_entry.get(qk)).strip():
                quote_val = str(raw_entry.get(qk)).strip()
                break
        if not quote_val:
            return None, (
                f"missing_verbatim_quote_in_line_manifest_for_{ks!r} "
                "(PREREG_H Addendum A requires exact sentence or manifest entry)"
            )

        parsed_hits.append(
            {
                "matched_key": ks,
                "file": hp_str,
                "file_name": hp_name,
                "line": norm_line,
                "url": url_val,
                "quote": quote_val,
            }
        )

    lines_seen = {p["line"] for p in parsed_hits}
    if len(lines_seen) > 1:
        return None, f"conflicting_manifest_entries_for_{hp_name}: {sorted(lines_seen)}"
    return parsed_hits[0], None


def _resolve_vcc_cell_line_identities(
    obs: h5py.Group,
    obs_info: Dict[str, Any],
    obs_len: int,
    hp: Optional[Path] = None,
    context_to_line_map: Optional[Dict[str, str]] = None,
    line_manifest: Optional[Mapping[str, Any]] = None,
) -> Dict[str, Any]:
    """
    `PREREG_H v2 §3.1 + Addendum A` Line-Identity Resolution (`(a) -> (b) -> (c)`):
      (a) Explicit column (`cell_line`, `cell_type`, `line`) normalized via removing `-`, `_`, spaces
          and lowercasing -> must match `cal33, hela, jurkat, ac16, molm13, u87mg` unambiguously.
      (b) Else `context` column mapped ONLY through published `context_to_line_map`.
      (c) Else official `line_manifest` (`file -> {line, url, quote}` from Arc release docs/README/manifest);
          filename alone without manifest URL + verbatim statement is strictly rejected.
      Conflict rule (`Addendum A`): whenever `(a)` or mapped `(b)` applies AND `(c)` also covers the file,
      `(a)/(b)` take precedence and every cell must agree with `(c)`; any disagreement marks all affected
      lines `EXCLUDED_FEASIBILITY`. Never infer cell line from expression.
    """
    manifest_entry: Optional[Dict[str, Any]] = None
    manifest_err: Optional[str] = None
    if hp is not None and line_manifest is not None:
        manifest_entry, manifest_err = _resolve_manifest_entry_for_path(hp, line_manifest)
    if manifest_err is not None:
        return {
            "unambiguous": False,
            "status": "EXCLUDED_FEASIBILITY",
            "reason": manifest_err,
        }

    exp_col = obs_info.get("explicit_line_col")
    ctx_col = obs_info.get("context_col")

    if exp_col is not None:
        raw_vals = _read_obs_or_var_column(obs, exp_col)
        if raw_vals is None or len(raw_vals) == 0:
            return {"unambiguous": False, "status": "EXCLUDED_FEASIBILITY", "reason": f"empty_explicit_line_column:{exp_col}"}
        resolved = [_normalize_line_token(str(v)) for v in raw_vals]
        bad = sorted({str(raw_vals[i]) for i, lb in enumerate(resolved) if lb not in ALL_VCC_KNOWN_LINES})
        if bad:
            return {
                "unambiguous": False,
                "status": "EXCLUDED_FEASIBILITY",
                "reason": f"ambiguous_or_unknown_line_in_{exp_col}: {bad[:10]}",
            }
        mode = f"EXPLICIT_LINE_COL:{exp_col}"
        if manifest_entry is not None:
            m_line = manifest_entry["line"]
            disagree = sorted({lb for lb in resolved if lb != m_line})
            if disagree:
                affected = sorted({m_line} | set(disagree))
                return {
                    "unambiguous": False,
                    "status": "EXCLUDED_FEASIBILITY",
                    "affected_lines": affected,
                    "manifest_provenance": manifest_entry,
                    "reason": (
                        f"LINE_IDENTITY_CONFLICT_PREREG_H_3_1_C: explicit column {exp_col!r} "
                        f"(observed={sorted(set(resolved))}) disagrees with official line_manifest "
                        f"({manifest_entry['file_name']} -> {m_line!r}); affected_lines={affected}"
                    ),
                }
            mode = f"{mode}+MANIFEST_VERIFIED"
        n_jurkat = sum(1 for lb in resolved if lb == "jurkat")
        return {
            "unambiguous": True,
            "mode": mode,
            "resolved_lines": resolved,
            "present_lines": sorted(set(resolved) - R15_EXCLUDED_LINES),
            "n_jurkat_pre_qc": int(n_jurkat),
            "manifest_provenance": manifest_entry,
        }

    if ctx_col is not None and context_to_line_map:
        raw_vals = _read_obs_or_var_column(obs, ctx_col)
        if raw_vals is None or len(raw_vals) == 0:
            return {"unambiguous": False, "status": "EXCLUDED_FEASIBILITY", "reason": f"empty_context_column:{ctx_col}"}
        # Match keys cleanly: exact stripped lowercase OR normalized (removing - _ spaces)
        exact_map: Dict[str, str] = {}
        norm_key_map: Dict[str, str] = {}
        for k, v in context_to_line_map.items():
            nv = _normalize_line_token(v)
            if nv not in ALL_VCC_KNOWN_LINES:
                return {"unambiguous": False, "status": "EXCLUDED_FEASIBILITY", "reason": f"invalid_target_line_in_context_map:{k!r}->{v!r}"}
            exact_map[str(k).strip().lower()] = nv
            norm_key_map[_normalize_line_token(k)] = nv

        resolved = []
        bad_ctx: Set[str] = set()
        for rv in raw_vals:
            s_raw = str(rv).strip()
            mapped = exact_map.get(s_raw.lower())
            if mapped is None:
                mapped = norm_key_map.get(_normalize_line_token(s_raw))
            if mapped is None:
                bad_ctx.add(s_raw)
                resolved.append("")
            else:
                resolved.append(mapped)
        if bad_ctx:
            return {
                "unambiguous": False,
                "status": "EXCLUDED_FEASIBILITY",
                "reason": f"unmapped_context_values_in_{ctx_col}: {sorted(bad_ctx)[:10]}",
            }
        mode = f"OFFICIAL_CONTEXT_MAP:{ctx_col}"
        if manifest_entry is not None:
            m_line = manifest_entry["line"]
            disagree = sorted({lb for lb in resolved if lb != m_line})
            if disagree:
                affected = sorted({m_line} | set(disagree))
                return {
                    "unambiguous": False,
                    "status": "EXCLUDED_FEASIBILITY",
                    "affected_lines": affected,
                    "manifest_provenance": manifest_entry,
                    "reason": (
                        f"LINE_IDENTITY_CONFLICT_PREREG_H_3_1_C: mapped context column {ctx_col!r} "
                        f"(observed={sorted(set(resolved))}) disagrees with official line_manifest "
                        f"({manifest_entry['file_name']} -> {m_line!r}); affected_lines={affected}"
                    ),
                }
            mode = f"{mode}+MANIFEST_VERIFIED"
        n_jurkat = sum(1 for lb in resolved if lb == "jurkat")
        return {
            "unambiguous": True,
            "mode": mode,
            "resolved_lines": resolved,
            "present_lines": sorted(set(resolved) - R15_EXCLUDED_LINES),
            "n_jurkat_pre_qc": int(n_jurkat),
            "manifest_provenance": manifest_entry,
        }

    if manifest_entry is not None:
        m_line = manifest_entry["line"]
        resolved = [m_line] * int(obs_len)
        n_jurkat = int(obs_len) if m_line == "jurkat" else 0
        return {
            "unambiguous": True,
            "mode": "OFFICIAL_LINE_MANIFEST",
            "resolved_lines": resolved,
            "present_lines": [] if m_line in R15_EXCLUDED_LINES else [m_line],
            "n_jurkat_pre_qc": n_jurkat,
            "manifest_provenance": manifest_entry,
        }

    if ctx_col is not None:
        return {
            "unambiguous": False,
            "status": "EXCLUDED_FEASIBILITY",
            "reason": (
                f"context_only_obs ({ctx_col}) without published Arc context_to_line_map or official "
                "line_manifest (PREREG_H v2 §3.1 + Addendum A: never inferred from expression or filename)"
            ),
        }

    return {
        "unambiguous": False,
        "status": "EXCLUDED_FEASIBILITY",
        "reason": "missing_explicit_line_col_context_map_and_official_line_manifest_PREREG_H_3_1_C",
    }


def _detect_vcc_obsm_guides(f: h5py.File, obs_len: int) -> Optional[Dict[str, Any]]:
    """Detect raw guide matrix in `/obsm/*` + names in `/uns/*` for `PREREG_F §3.2` fallback."""
    if "obsm" not in f or not isinstance(f["obsm"], h5py.Group):
        return None
    for ok in sorted(f["obsm"].keys()):
        if not any(tok in ok.lower() for tok in ("sgrna", "guide", "crispr")):
            continue
        node = f["obsm"][ok]
        g_names: Optional[List[str]] = None
        if "uns" in f and isinstance(f["uns"], h5py.Group):
            for uk in (f"{ok}_names", "sgrna_ids", "guide_ids", "crispr_guide_ids"):
                if uk in f["uns"]:
                    g_names = _decode_str(f["uns"][uk][()])
                    break
        if g_names is None:
            continue
        if isinstance(node, h5py.Dataset) and node.ndim == 2 and node.shape == (obs_len, len(g_names)):
            return {"source": f"obsm/{ok}", "encoding": "dense", "guide_names": g_names}
        if isinstance(node, h5py.Group) and {"data", "indices", "indptr"}.issubset(set(node.keys())):
            enc = "csc" if "csc" in str(node.attrs.get("encoding-type", "")).lower() else "csr"
            return {"source": f"obsm/{ok}", "encoding": enc, "guide_names": g_names}
    return None


def _map_guide_name_to_target(guide_id: str, hgnc: HGNCLookup) -> Tuple[str, str, str]:
    g_str = str(guide_id).strip()
    g_low = g_str.lower()
    if any(tok in g_low for tok in ("gene_desert", "safe-targeting", "safe_targeting")):
        return "EXCLUDED_SAFE_OR_BAD", "", ""
    if g_low.startswith("non-targeting") or g_low.startswith("ntc"):
        return "NTC", "", ""
    stem = g_str.rsplit("_", 1)[0] if "_" in g_str else (g_str.rsplit("-", 1)[0] if "-" in g_str else g_str)
    canon_s, hstat = hgnc.resolve(stem, "")
    if canon_s and (hstat != "unmapped" or canon_s in hgnc.approved):
        return "TARGET", canon_s, stem
    return "EXCLUDED_UNMAPPED", "", stem


def check_vcc_h5ad_feasibility(
    line_name: str,
    h5_paths: Sequence[Path],
    hgnc: Optional[HGNCLookup] = None,
    context_to_line_map: Optional[Dict[str, str]] = None,
    line_manifest: Optional[Mapping[str, Any]] = None,
) -> Dict[str, Any]:
    """
    `PREREG_H v2 §3, §3.1 + Addendum A, §5` Pre-Metric Technical Feasibility Gate on one candidate line (`cal33..u87mg`):
      1. Unconditionally blocks `jurkat`.
      2. Resolves cell-line identity strictly via `§3.1` (`(a) cell_line/cell_type/line` ->
         `(b) context` via published `context_to_line_map` -> `(c)` official `line_manifest`
         with documented URL + verbatim statement -> else `EXCLUDED_FEASIBILITY`), enforcing
         `(a)/(b)` vs `(c)` conflict exclusion and dropping all `jurkat` rows pre-QC.
      3. Checks `/raw/X`, `/X`, and non-excluded `/layers/*` across 100% of stored non-zero entries
         for finite non-negative exact integers (`>= 0`, `x == floor(x)`).
      4. Stops with `EXCLUDED_FEASIBILITY` (`MULTIPLE_INCONSISTENT_INTEGER_SLOTS`) if `>= 2`
         non-excluded integer slots disagree in shape, `/var` genes, or matrix values.
      5. Verifies explicit NTC (`target_gene == "non-targeting"`) within `line_name` cells
         (or fallback raw guide counts in `/obsm/*`).
    """
    norm_line = _normalize_line_token(line_name)
    if norm_line in R15_EXCLUDED_LINES:
        return {
            "line": line_name,
            "feasible": False,
            "status": "EXCLUDED_LINE_JURKAT",
            "reason": "jurkat_excluded_homologous_to_r9_source_nadig_jurkat",
        }
    if norm_line not in R15_CANDIDATE_SET:
        return {
            "line": line_name,
            "feasible": False,
            "status": "EXCLUDED_UNKNOWN_LINE",
            "reason": f"unknown_vcc_candidate_line:{line_name!r}",
        }
    if not h5_paths:
        return {"line": line_name, "feasible": False, "status": "EXCLUDED_FEASIBILITY", "reason": "no_h5ad_files"}

    _ = hgnc or build_hgnc_lookup()
    file_reports: List[Dict[str, Any]] = []
    other_line_files: List[Dict[str, Any]] = []
    file_audit_log: List[Dict[str, Any]] = []
    manifest_provenance: List[Dict[str, Any]] = []
    modes: List[str] = []

    for hp in h5_paths:
        if not hp.exists():
            return {"line": line_name, "feasible": False, "status": "EXCLUDED_FEASIBILITY", "reason": f"missing:{hp}"}
        try:
            with h5py.File(hp, "r") as f:
                if "obs" not in f or "var" not in f:
                    return {
                        "line": line_name,
                        "feasible": False,
                        "status": "EXCLUDED_FEASIBILITY",
                        "reason": f"MISSING_ANNDATA_GROUP in {hp.name}",
                    }
                barcodes = _read_anndata_index(f["obs"], "obs")
                obs_len = len(barcodes)

                obs_info = _inspect_vcc_obs_columns(f["obs"])
                line_res = _resolve_vcc_cell_line_identities(
                    f["obs"],
                    obs_info,
                    obs_len=obs_len,
                    hp=hp,
                    context_to_line_map=context_to_line_map,
                    line_manifest=line_manifest,
                )
                if not line_res["unambiguous"]:
                    affected = line_res.get("affected_lines", [])
                    if affected and norm_line not in affected:
                        rec_other = {
                            "file": str(hp),
                            "obs_len": obs_len,
                            "contributes_cells": False,
                            "n_line_raw_cells": 0,
                            "status": "OTHER_LINE_CONFLICT_ONLY",
                            "affected_lines": affected,
                            "reason": line_res["reason"],
                        }
                        other_line_files.append(rec_other)
                        file_audit_log.append(rec_other)
                        continue
                    return {
                        "line": line_name,
                        "feasible": False,
                        "status": line_res.get("status", "EXCLUDED_FEASIBILITY"),
                        "reason": line_res["reason"],
                        "affected_lines": affected,
                        "manifest_provenance": line_res.get("manifest_provenance"),
                    }
                resolved_lines: List[str] = line_res["resolved_lines"]
                # Drop ALL jurkat cells immediately at ingestion before any QC/NTC check (`§3.1`)
                n_jurkat_pre = int(line_res["n_jurkat_pre_qc"])
                line_raw_rows = [r for r, lb in enumerate(resolved_lines) if lb == norm_line]
                if not line_raw_rows:
                    rec_other = {
                        "file": str(hp),
                        "obs_len": obs_len,
                        "contributes_cells": False,
                        "n_line_raw_cells": 0,
                        "status": "ZERO_CELLS_OTHER_LINE_ONLY",
                        "line_identity_mode": line_res["mode"],
                        "present_lines": line_res["present_lines"],
                        "n_jurkat_dropped_pre_qc": n_jurkat_pre,
                        "manifest_provenance": line_res.get("manifest_provenance"),
                    }
                    other_line_files.append(rec_other)
                    file_audit_log.append(rec_other)
                    continue

                candidate_paths: List[str] = []
                excluded_layers: List[str] = []
                if "raw" in f and isinstance(f["raw"], h5py.Group) and "X" in f["raw"]:
                    candidate_paths.append("/raw/X")
                if "X" in f:
                    candidate_paths.append("/X")
                if "layers" in f and isinstance(f["layers"], h5py.Group):
                    for lk in sorted(f["layers"].keys()):
                        if _is_excluded_layer_name(lk):
                            excluded_layers.append(f"/layers/{lk}")
                        else:
                            candidate_paths.append(f"/layers/{lk}")

                slot_checks: List[Dict[str, Any]] = []
                integer_slots: List[Dict[str, Any]] = []
                for spath in candidate_paths:
                    chk = _full_slot_integer_check(f, spath, obs_len=obs_len)
                    slot_checks.append(
                        {
                            "slot_path": spath,
                            "is_integer": chk["is_integer"],
                            "reason": chk["reason"],
                            "shape": chk.get("shape"),
                            "nnz_pos": chk.get("nnz_pos"),
                            "sum_total": chk.get("sum_total"),
                        }
                    )
                    if chk["is_integer"]:
                        integer_slots.append(chk)

                if not integer_slots:
                    return {
                        "line": line_name,
                        "feasible": False,
                        "status": "EXCLUDED_FEASIBILITY",
                        "reason": f"no_integer_raw_umi_slot in {hp.name}: {slot_checks}",
                        "slot_checks": slot_checks,
                    }

                chosen = integer_slots[0]
                if len(integer_slots) >= 2:
                    for other in integer_slots[1:]:
                        eq_ok, eq_msg = _compare_two_integer_slots(f, chosen, other)
                        if not eq_ok:
                            return {
                                "line": line_name,
                                "feasible": False,
                                "status": "EXCLUDED_FEASIBILITY",
                                "reason": f"MULTIPLE_INCONSISTENT_INTEGER_SLOTS in {hp.name}: {eq_msg}",
                                "slot_checks": slot_checks,
                            }
                    slot_consistency = f"MULTIPLE_IDENTICAL_INTEGER_SLOTS_VERIFIED ({[s['slot_path'] for s in integer_slots]})"
                else:
                    slot_consistency = f"UNIQUE_INTEGER_SLOT ({chosen['slot_path']})"

                prov = line_res.get("manifest_provenance")
                if prov is not None:
                    prov_rec = dict(prov)
                    prov_rec["identity_mode"] = line_res["mode"]
                    manifest_provenance.append(prov_rec)
                else:
                    prov_rec = None

                if obs_info["target_col"] is None:
                    obsm_g = _detect_vcc_obsm_guides(f, obs_len)
                    if obsm_g is None:
                        return {
                            "line": line_name,
                            "feasible": False,
                            "status": "EXCLUDED_FEASIBILITY",
                            "reason": f"missing_target_gene_column_and_raw_guide_obsm in {hp.name}: obs_cols={obs_info['obs_cols']}",
                        }
                    modes.append("PREREG_F_OBSM_GUIDE_FALLBACK")
                    rep = {
                        "file": str(hp),
                        "obs_len": obs_len,
                        "contributes_cells": True,
                        "n_line_raw_cells": len(line_raw_rows),
                        "chosen_slot": chosen["slot_path"],
                        "encoding": chosen["encoding"],
                        "shape": chosen["shape"],
                        "slot_consistency": slot_consistency,
                        "slot_checks": slot_checks,
                        "excluded_layers": excluded_layers,
                        "obs_info": obs_info,
                        "line_identity_mode": line_res["mode"],
                        "manifest_provenance": prov_rec,
                        "n_jurkat_dropped_pre_qc": n_jurkat_pre,
                        "resolved_lines": resolved_lines,
                        "obsm_guides": obsm_g,
                    }
                    file_reports.append(rep)
                    file_audit_log.append(
                        {
                            "file": str(hp),
                            "obs_len": obs_len,
                            "contributes_cells": True,
                            "n_line_raw_cells": len(line_raw_rows),
                            "status": "CONTRIBUTING",
                            "line_identity_mode": line_res["mode"],
                            "n_jurkat_dropped_pre_qc": n_jurkat_pre,
                        }
                    )
                    continue

                tgt_vals = _read_obs_or_var_column(f["obs"], obs_info["target_col"])
                if tgt_vals is None or not any(str(tgt_vals[r]).strip().lower() == "non-targeting" for r in line_raw_rows):
                    return {
                        "line": line_name,
                        "feasible": False,
                        "status": "EXCLUDED_FEASIBILITY",
                        "reason": f"missing_target_gene_non_targeting_controls in {hp.name} for {line_name}",
                    }
                modes.append("OBS_TARGET_GENE_CALLS")

                rep = {
                    "file": str(hp),
                    "obs_len": obs_len,
                    "contributes_cells": True,
                    "n_line_raw_cells": len(line_raw_rows),
                    "chosen_slot": chosen["slot_path"],
                    "encoding": chosen["encoding"],
                    "shape": chosen["shape"],
                    "slot_consistency": slot_consistency,
                    "slot_checks": slot_checks,
                    "excluded_layers": excluded_layers,
                    "obs_info": obs_info,
                    "line_identity_mode": line_res["mode"],
                    "manifest_provenance": prov_rec,
                    "n_jurkat_dropped_pre_qc": n_jurkat_pre,
                    "resolved_lines": resolved_lines,
                }
                file_reports.append(rep)
                file_audit_log.append(
                    {
                        "file": str(hp),
                        "obs_len": obs_len,
                        "contributes_cells": True,
                        "n_line_raw_cells": len(line_raw_rows),
                        "status": "CONTRIBUTING",
                        "line_identity_mode": line_res["mode"],
                        "n_jurkat_dropped_pre_qc": n_jurkat_pre,
                    }
                )
        except Exception as exc:
            return {
                "line": line_name,
                "feasible": False,
                "status": "EXCLUDED_FEASIBILITY",
                "reason": f"hdf5_read_error in {hp.name}: {exc}",
            }

    if not file_reports:
        return {
            "line": line_name,
            "feasible": False,
            "status": "EXCLUDED_FEASIBILITY",
            "reason": (
                f"line {line_name!r} contributed 0 cells across all {len(h5_paths)} input files: "
                f"other_lines_seen={sorted({ln for r in other_line_files for ln in r.get('present_lines', [])})}"
            ),
            "other_line_files": other_line_files,
            "file_audit_log": file_audit_log,
        }

    if len(set(modes)) != 1:
        return {
            "line": line_name,
            "feasible": False,
            "status": "EXCLUDED_FEASIBILITY",
            "reason": f"INCONSISTENT_GUIDE_MODES: {modes}",
            "files": file_reports,
            "other_line_files": other_line_files,
            "file_audit_log": file_audit_log,
        }

    return {
        "line": line_name,
        "feasible": True,
        "status": "FEASIBLE",
        "guide_call_mode": modes[0],
        "manifest_provenance": manifest_provenance,
        "files": file_reports,
        "other_line_files": other_line_files,
        "file_audit_log": file_audit_log,
    }


def _classify_vcc_obs_cell(
    target_raw: str,
    guide_raw: str,
    ntc_id_raw: str,
    n_guides_val: float,
    hgnc: HGNCLookup,
) -> Tuple[str, str, str]:
    """
    Classify one cell across 2025 (`target_gene`, `guide_id`) or 2026 (`target_gene`, `ntc_id`) `.obs` schemas.
    Returns `(role, canon_sym, raw_sym)` where `role` in:
      `NTC`, `TARGET`, `EXCLUDED_MULTIGUIDE`, `EXCLUDED_SAFE_OR_BAD`, `EXCLUDED_UNMAPPED`.
    """
    if np.isfinite(n_guides_val) and n_guides_val > 1.0 + 1e-6:
        return "EXCLUDED_MULTIGUIDE", "", ""

    t_str = target_raw.strip()
    g_str = guide_raw.strip()
    n_str = ntc_id_raw.strip()
    t_low = t_str.lower()
    g_low = g_str.lower()
    n_low = n_str.lower()

    if any(tok in (t_low, g_low, n_low) for tok in ("multiplet", "doublet")):
        return "EXCLUDED_MULTIGUIDE", "", ""
    if t_low in _BAD_TARGET_TOKENS:
        return "EXCLUDED_SAFE_OR_BAD", "", ""
    if any(tok in (g_low, n_low) for tok in ("safe-targeting", "safe_targeting", "gene_desert")):
        return "EXCLUDED_SAFE_OR_BAD", "", ""

    # Exact `PREREG_H §5` NTC check: `target_gene == "non-targeting"`
    if t_low == "non-targeting":
        if g_low in ("unassigned", "multiplet", "doublet"):
            return "EXCLUDED_MULTIGUIDE", "", ""
        if n_low in ("multiplet", "doublet"):
            return "EXCLUDED_MULTIGUIDE", "", ""
        if g_str and _MULTI_TARGET_DELIM_RE.search(g_str):
            toks = [x.strip().lower() for x in _MULTI_TARGET_DELIM_RE.split(g_str) if x.strip()]
            if any(not ("non-targeting" in tk or tk.startswith("ntc")) for tk in toks):
                return "EXCLUDED_MULTIGUIDE", "", ""
        return "NTC", "", ""

    # Targeting gene must be single-gene token (reject heterogeneous `GENE1,GENE2` or spaces)
    if _MULTI_TARGET_DELIM_RE.search(t_str):
        return "EXCLUDED_MULTIGUIDE", "", ""

    canon_sym, hstat = hgnc.resolve(t_str, "")
    is_mapped = bool(canon_sym) and (hstat != "unmapped" or canon_sym in hgnc.approved)
    if not is_mapped:
        return "EXCLUDED_UNMAPPED", "", t_str

    # Validate optional `guide_id`: permit 2025 dual-guide single-construct cassettes (`GENE_1|GENE_2`)
    # ONLY if all component tokens refer to the same construct/gene; reject true mixed-target guides.
    if g_str and g_low not in ("na", "nan", "none", "null", "-", "."):
        if g_low == "unassigned":
            return "EXCLUDED_MULTIGUIDE", "", ""
        if "non-targeting" in g_low:
            return "EXCLUDED_MULTIGUIDE", "", ""
        if _MULTI_TARGET_DELIM_RE.search(g_str):
            toks = [x.strip() for x in _MULTI_TARGET_DELIM_RE.split(g_str) if x.strip()]
            stems: Set[str] = set()
            for tk in toks:
                stem = tk.split("_")[0].split("-")[0].strip()
                c_tk, s_tk = hgnc.resolve(stem, "")
                if s_tk != "unmapped" or c_tk in hgnc.approved:
                    stems.add(c_tk)
                elif stem:
                    stems.add(stem.upper())
            if len(stems) > 1 or (stems and canon_sym.upper() not in {x.upper() for x in stems} and t_str.upper() not in {x.upper() for x in stems}):
                return "EXCLUDED_MULTIGUIDE", "", ""

    return "TARGET", canon_sym, t_str


def _apply_tier2_seed20261231_split(
    nt_mask: np.ndarray,
    canon_targets: Sequence[str],
    barcodes: Sequence[str],
    seed: int = DEFAULT_R15_SEED,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    `PREREG_H §6` Tier-2 Fallback:
      Within each canonical target stratum (`canon_sym`) and `"__NTC__"`:
      sort kept cell indices ascending by `(barcode, idx)` and apply one generator
      `rng = np.random.default_rng(20261231)` iterating sorted strata keys (`50/50` A/B split).
    """
    n_kept = len(nt_mask)
    rng = np.random.default_rng(seed)
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


def build_vcc_line_audit_parts(
    line_name: str,
    h5_paths: Sequence[Path],
    hgnc: HGNCLookup,
    context_to_line_map: Optional[Dict[str, str]] = None,
    line_manifest: Optional[Mapping[str, Any]] = None,
    min_umi: float = 300.0,
    min_genes: int = 200,
    min_batch_nt: int = 15,
    seed: int = DEFAULT_R15_SEED,
) -> Tuple[List[PartSpec], Dict[str, Any]]:
    """
    Independent clean-room `VCC 2025/2026` `.h5ad` reader (`PREREG_H v2 §3, §3.1 + Addendum A, §5, §6`):
      - Resolves cell lines via `§3.1` (`(a) explicit cell_line/cell_type/line` -> `(b) context` via
        published `context_to_line_map` -> `(c)` official `line_manifest` with URL + exact sentence ->
        else `EXCLUDED_FEASIBILITY`).
      - Drops ALL `jurkat` rows upfront at ingestion BEFORE any UMI/QC count, centering, or metric.
      - Enforces integer raw UMIs + multi-slot guard, excludes `Cas9`/`ENSMUSG` before RNA QC,
        keeps `>= 300 UMI`, `>= 200 genes`, `max_pct_mt=None`.
      - Two-tier split-half: Tier 1 `.obs` channel/batch parity (`>=2` units & `>=15` NTC/half)
        with unit-level subtraction, else Tier 2 barcode-sorted stratified `seed=20261231`.
    """
    feas = check_vcc_h5ad_feasibility(
        line_name=line_name,
        h5_paths=h5_paths,
        hgnc=hgnc,
        context_to_line_map=context_to_line_map,
        line_manifest=line_manifest,
    )
    if not feas["feasible"]:
        raise RuntimeError(f"EXCLUDED_FEASIBILITY ({feas['status']}): {feas['reason']}")

    norm_line = _normalize_line_token(line_name)
    counts_audit: Dict[str, int] = defaultdict(int)
    file_records: List[Dict[str, Any]] = []

    for frep in feas["files"]:
        hp = Path(frep["file"])
        slot_path = str(frep["chosen_slot"])
        encoding = str(frep["encoding"])
        n_rows_raw, n_cols_raw = int(frep["shape"][0]), int(frep["shape"][1])
        obs_info = frep["obs_info"]
        counts_audit["n_jurkat_dropped_pre_qc"] += int(frep["n_jurkat_dropped_pre_qc"])

        # Pre-QC filter: keep ONLY indices for target `norm_line` (jurkat and other lines excluded upfront)
        resolved_lines: List[str] = frep["resolved_lines"]
        line_raw_rows = np.asarray([r for r, lb in enumerate(resolved_lines) if lb == norm_line], dtype=np.int64)
        n_line_raw = int(line_raw_rows.size)
        counts_audit["n_line_raw_cells"] += n_line_raw

        with h5py.File(hp, "r") as f:
            barcodes = _read_anndata_index(f["obs"], "obs")
            var_key = "raw/var" if slot_path == "/raw/X" else "var"
            vm = _extract_vcc_var_metadata(f[var_key], var_key)
            v_syms = vm["var_syms"]
            v_ens = vm["var_ens"]

            keep_gene_mask = np.asarray(
                [
                    not (
                        ("CAS9" in s.upper())
                        or ("CAS9" in e.upper())
                        or e.upper().startswith("ENSMUSG")
                        or s.lower().startswith("mm10")
                    )
                    for s, e in zip(v_syms, v_ens)
                ],
                dtype=bool,
            )
            keep_gene_cols = np.flatnonzero(keep_gene_mask)
            if keep_gene_cols.size == 0:
                raise RuntimeError(f"Zero human RNA features remaining in {hp.name}")

            tgt_arr = _read_obs_or_var_column(f["obs"], obs_info["target_col"]) if obs_info["target_col"] else None
            guide_arr = _read_obs_or_var_column(f["obs"], obs_info["guide_col"]) if obs_info["guide_col"] else None
            ntc_id_arr = _read_obs_or_var_column(f["obs"], obs_info["ntc_id_col"]) if obs_info["ntc_id_col"] else None
            n_g_arr = (
                np.asarray(pd.to_numeric(_read_obs_or_var_column(f["obs"], obs_info["n_guides_col"]), errors="coerce"), dtype=np.float64)
                if obs_info["n_guides_col"]
                else None
            )
            batch_col = obs_info["batch_col"]
            batch_arr = _read_obs_or_var_column(f["obs"], batch_col) if batch_col else None

        rc_line, _, _ = make_slot_subset_readers(hp, slot_path, line_raw_rows, n_rows_raw, n_cols_raw, encoding)
        umi_line = np.zeros(n_line_raw, dtype=np.float64)
        ngenes_line = np.zeros(n_line_raw, dtype=np.int64)
        for s0 in range(0, n_line_raw, 16384):
            e0 = min(n_line_raw, s0 + 16384)
            blk = rc_line(s0, e0)[:, keep_gene_cols]
            umi_line[s0:e0] = blk.sum(axis=1)
            ngenes_line[s0:e0] = (blk > 0).sum(axis=1)

        fallback_roles: Optional[List[Tuple[str, str, str]]] = None
        if tgt_arr is None and "obsm_guides" in frep:
            og = frep["obsm_guides"]
            g_names = og["guide_names"]
            g_rc, _, _ = make_slot_subset_readers(hp, og["source"], line_raw_rows, n_rows_raw, len(g_names), og["encoding"])
            g_mat = g_rc(0, n_line_raw)
            init_mask = ngenes_line >= int(min_genes)
            hist = np.zeros((len(g_names), 66), dtype=np.int64)
            rr_i, cc_i = np.nonzero(g_mat[init_mask] > 0)
            if rr_i.size > 0:
                vals_i = np.clip(g_mat[init_mask][rr_i, cc_i].astype(np.int64), 1, 65)
                np.add.at(hist, (cc_i.astype(np.int32), vals_i), 1)
            thrs = fit_poisson_gaussian_thresholds_from_hist(hist, min_thr=5, max_thr=20, default_thr=7)
            csc = sp.csc_matrix(g_mat.T, dtype=np.float64)
            t1, t2, a1 = extract_top2_csc(csc)
            pos_m = csc.data >= thrs[csc.indices]
            pos_csc = sp.csc_matrix((pos_m.astype(np.int8), csc.indices, csc.indptr), shape=csc.shape)
            n_pos = np.asarray(pos_csc.sum(axis=0)).ravel()
            s1_mask = (n_pos == 1) & (t1 >= 5.0) & ((t2 < 3.0) | (t1 / np.maximum(t2, 0.5) >= 3.0))
            fallback_roles = [
                _map_guide_name_to_target(g_names[int(a1[r])], hgnc) if s1_mask[r] else ("EXCLUDED_MULTIGUIDE", "", "")
                for r in range(n_line_raw)
            ]

        keep_rows_f: List[int] = []
        umi_kept_f: List[float] = []
        nt_mask_f: List[bool] = []
        pkeys_f: List[str] = []
        bcs_f: List[str] = []
        batches_f: List[str] = []
        pmeta_f: Dict[str, Tuple[str, str, str, str]] = {}

        for loc_r in range(n_line_raw):
            raw_r = int(line_raw_rows[loc_r])
            if not (umi_line[loc_r] >= float(min_umi) and ngenes_line[loc_r] >= int(min_genes)):
                counts_audit["n_failed_qc_umi_or_genes"] += 1
                continue

            if fallback_roles is not None:
                role, canon_s, raw_s = fallback_roles[loc_r]
            else:
                t_raw = str(tgt_arr[raw_r]) if tgt_arr is not None else ""
                g_raw = str(guide_arr[raw_r]) if guide_arr is not None else ""
                n_raw = str(ntc_id_arr[raw_r]) if ntc_id_arr is not None else ""
                ng_v = float(n_g_arr[raw_r]) if n_g_arr is not None else float("nan")
                role, canon_s, raw_s = _classify_vcc_obs_cell(t_raw, g_raw, n_raw, ng_v, hgnc)

            if role == "NTC":
                keep_rows_f.append(raw_r)
                umi_kept_f.append(float(umi_line[loc_r]))
                nt_mask_f.append(True)
                pkeys_f.append("")
                bcs_f.append(barcodes[raw_r])
                batches_f.append(str(batch_arr[raw_r]).strip() if batch_arr is not None else "")
                counts_audit["n_kept_ntc"] += 1
            elif role == "TARGET":
                keep_rows_f.append(raw_r)
                umi_kept_f.append(float(umi_line[loc_r]))
                nt_mask_f.append(False)
                pkeys_f.append(canon_s)
                bcs_f.append(barcodes[raw_r])
                batches_f.append(str(batch_arr[raw_r]).strip() if batch_arr is not None else "")
                if canon_s not in pmeta_f:
                    _, hstat = hgnc.resolve(raw_s, "")
                    ens = hgnc.sym_to_ens.get(canon_s, "")
                    pmeta_f[canon_s] = (ens, canon_s, raw_s, hstat)
                counts_audit["n_kept_targeting_cells"] += 1
            else:
                counts_audit[f"n_{role.lower()}"] += 1

        if keep_rows_f:
            file_records.append(
                {
                    "file": hp,
                    "slot_path": slot_path,
                    "encoding": encoding,
                    "n_rows_raw": n_rows_raw,
                    "n_cols_raw": n_cols_raw,
                    "keep_gene_cols": keep_gene_cols,
                    "v_ens": [v_ens[j] for j in keep_gene_cols],
                    "v_sym": [v_syms[j] for j in keep_gene_cols],
                    "umi_kept": np.asarray(umi_kept_f, dtype=np.float64),
                    "keep_arr": np.asarray(keep_rows_f, dtype=np.int64),
                    "nt_arr": np.asarray(nt_mask_f, dtype=bool),
                    "pkeys": pkeys_f,
                    "bcs": bcs_f,
                    "batch_col": batch_col,
                    "batch_labels": batches_f,
                    "pmeta": pmeta_f,
                    "line_identity_mode": frep["line_identity_mode"],
                }
            )

    if not file_records:
        raise RuntimeError(f"Line {line_name!r} yielded zero QC-valid cells")

    # Assemble combined cell list to evaluate PREREG_H §6 Tier 1 vs Tier 2 cleanly across files
    all_nt = np.concatenate([rec["nt_arr"] for rec in file_records])
    all_pk: List[str] = []
    all_bcs: List[str] = []
    all_batch_labels: List[str] = []
    for rec in file_records:
        all_pk.extend(rec["pkeys"])
        all_bcs.extend(rec["bcs"])
        all_batch_labels.extend(rec["batch_labels"])

    total_nt = int(np.sum(all_nt))
    if total_nt == 0:
        raise RuntimeError(f"Line {line_name!r} yielded zero NTC cells")

    uniq_batches = sorted({b for b in all_batch_labels if b and b.lower() not in ("na", "nan", "none")})
    split_tier = "TIER2_STRATIFIED_SEED20261231"
    halves_all = np.zeros(len(all_nt), dtype=np.int8)
    batches_final = np.empty(len(all_nt), dtype=object)
    n_nt_fallback = 0
    tier1_used_col = None

    if len(uniq_batches) >= 2 and len(uniq_batches) == len({b for b in all_batch_labels if b}):
        # Strict lexicographical 1-based parity: odd (1,3,...) -> Half A (1), even (2,4,...) -> Half B (2)
        b_to_half = {b: np.int8(1 if (idx % 2 == 1) else 2) for idx, b in enumerate(uniq_batches, start=1)}
        cand_h = np.asarray([b_to_half.get(b, 0) for b in all_batch_labels], dtype=np.int8)
        if np.all(cand_h > 0):
            nt_a = int(np.sum(all_nt & (cand_h == 1)))
            nt_b = int(np.sum(all_nt & (cand_h == 2)))
            if nt_a >= min_batch_nt and nt_b >= min_batch_nt:
                split_tier = "TIER1_OBS_BATCH_PARITY"
                halves_all = cand_h
                tier1_used_col = next((r["batch_col"] for r in file_records if r["batch_col"]), "batch")
                batches_final, n_nt_fallback = _apply_min_ntc_batch_fallback(
                    np.asarray(all_batch_labels, dtype=object),
                    all_nt,
                    halves_all,
                    group_labels=None,
                    min_batch_nt=min_batch_nt,
                )

    if split_tier != "TIER1_OBS_BATCH_PARITY":
        halves_all, batches_final = _apply_tier2_seed20261231_split(
            nt_mask=all_nt,
            canon_targets=all_pk,
            barcodes=all_bcs,
            seed=seed,
        )

    parts: List[PartSpec] = []
    offset = 0
    unique_kept_targets: Set[str] = set()
    for p_idx, rec in enumerate(file_records):
        n_k = len(rec["keep_arr"])
        sub_h = halves_all[offset : offset + n_k]
        sub_b = np.asarray(batches_final[offset : offset + n_k], dtype=str)
        offset += n_k
        for pk in rec["pkeys"]:
            if pk:
                unique_kept_targets.add(pk)

        keep_gcols = rec["keep_gene_cols"]
        raw_rc, raw_rr, raw_rcc = make_slot_subset_readers(
            h5_path=rec["file"],
            slot_path=rec["slot_path"],
            keep_rows=rec["keep_arr"],
            n_rows_raw=rec["n_rows_raw"],
            n_cols_raw=rec["n_cols_raw"],
            encoding=rec["encoding"],
        )

        def _make_fns(
            _rc: Callable[[int, int], np.ndarray] = raw_rc,
            _rr: Callable[[np.ndarray], np.ndarray] = raw_rr,
            _rcc: Callable[[int, int, np.ndarray], np.ndarray] = raw_rcc,
            _gcols: np.ndarray = keep_gcols,
        ) -> Tuple[
            Callable[[int, int], np.ndarray],
            Callable[[np.ndarray], np.ndarray],
            Callable[[int, int, np.ndarray], np.ndarray],
        ]:
            return (
                (lambda s, e: _rc(s, e)[:, _gcols]),
                (lambda loc: _rr(loc)[:, _gcols]),
                (lambda s, e, fcols: _rcc(s, e, _gcols[fcols])),
            )

        fn_rc, fn_rr, fn_rcc = _make_fns()
        parts.append(
            PartSpec(
                tag=f"VCC_{norm_line}_part{p_idx}",
                path=rec["file"],
                n_cells=n_k,
                n_genes_raw=len(keep_gcols),
                nt_mask=rec["nt_arr"],
                batches=sub_b,
                half_flags=sub_h,
                umi=rec["umi_kept"],
                pert_keys=rec["pkeys"],
                pert_meta=rec["pmeta"],
                var_gene_ids=rec["v_ens"],
                var_gene_names=rec["v_sym"],
                read_chunk_fn=fn_rc,
                read_rows_fn=fn_rr,
                read_chunk_cols_fn=fn_rcc,
            )
        )

    nt_a_final = int(np.sum(all_nt & (halves_all == 1)))
    nt_b_final = int(np.sum(all_nt & (halves_all == 2)))
    other_files = feas.get("other_line_files", [])
    counts_audit["n_files_contributing"] = len(feas["files"])
    counts_audit["n_files_zero_cells_other_lines"] = len(other_files)
    counts_audit["n_jurkat_dropped_in_other_files"] = sum(
        int(r.get("n_jurkat_dropped_pre_qc", 0)) for r in other_files
    )
    manifest_provs = []
    for frep in feas["files"]:
        if frep.get("manifest_provenance"):
            manifest_provs.append(dict(frep["manifest_provenance"]))
    audit_meta = {
        "line": norm_line,
        "feasibility": feas,
        "line_identity_modes": [frep["line_identity_mode"] for frep in feas["files"]],
        "line_manifest_provenance": manifest_provs,
        "other_line_files": other_files,
        "file_audit_log": feas.get("file_audit_log", []),
        "split_tier": split_tier,
        "tier1_batch_col": tier1_used_col,
        "n_tier1_batches": len(uniq_batches),
        "n_nt_fallback_cells": int(n_nt_fallback),
        "n_nt_half_a": nt_a_final,
        "n_nt_half_b": nt_b_final,
        "max_pct_mt": None,
        "min_umi": float(min_umi),
        "min_genes": int(min_genes),
        "n_unique_targets_kept": int(len(unique_kept_targets)),
        "counts": dict(counts_audit),
    }
    return parts, audit_meta


def partition_vcc_files_by_line(
    h5_paths: Sequence[Path],
    context_to_line_map: Optional[Dict[str, str]] = None,
    line_manifest: Optional[Mapping[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Scan a release batch of `.h5ad` files under `PREREG_H v2 §3.1 + Addendum A`:
      - Resolves cell-line identity per file using `(a) explicit column -> (b) mapped context ->
        (c) official line_manifest` (`file -> {line, url, quote}`).
      - Enforces pre-QC exclusion of `jurkat` (`EXCLUDED_LINE_JURKAT`).
      - Order-independently propagates `(a)/(b)` vs `(c)` disagreements: if ANY file disagrees
        between `(a)/(b)` and `(c)`, EVERY affected cell line across the batch is marked
        `EXCLUDED_FEASIBILITY` so corrupted metadata cannot leak into another clean file's verdict.
    """
    files_by_line: Dict[str, List[Path]] = defaultdict(list)
    provenance_by_line: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    excluded_lines: Dict[str, str] = {"jurkat": "EXCLUDED_LINE_JURKAT: homologous to r9_6 nadig_jurkat"}
    unresolved_files: List[Dict[str, Any]] = []
    n_jurkat_dropped_pre_qc = 0

    for hp in h5_paths:
        if not hp.exists():
            unresolved_files.append({"file": str(hp), "status": "EXCLUDED_FEASIBILITY", "reason": f"missing:{hp}"})
            continue
        m_entry, _ = _resolve_manifest_entry_for_path(hp, line_manifest)
        try:
            with h5py.File(hp, "r") as f:
                if "obs" not in f or "var" not in f:
                    unresolved_files.append(
                        {"file": str(hp), "status": "EXCLUDED_FEASIBILITY", "reason": f"MISSING_ANNDATA_GROUP in {hp.name}"}
                    )
                    continue
                obs_len = len(_read_anndata_index(f["obs"], "obs"))
                obs_info = _inspect_vcc_obs_columns(f["obs"])
                line_res = _resolve_vcc_cell_line_identities(
                    f["obs"],
                    obs_info,
                    obs_len=obs_len,
                    hp=hp,
                    context_to_line_map=context_to_line_map,
                    line_manifest=line_manifest,
                )
        except Exception as exc:
            unresolved_files.append(
                {"file": str(hp), "status": "EXCLUDED_FEASIBILITY", "reason": f"hdf5_read_error in {hp.name}: {exc}"}
            )
            continue

        if not line_res["unambiguous"]:
            affected = set(line_res.get("affected_lines", []))
            if m_entry is not None:
                affected.add(m_entry["line"])
            for ln in sorted(affected):
                if ln in R15_CANDIDATE_SET:
                    excluded_lines[ln] = f"EXCLUDED_FEASIBILITY: {line_res['reason']}"
                    files_by_line.pop(ln, None)
                    provenance_by_line.pop(ln, None)
            unresolved_files.append(
                {
                    "file": str(hp),
                    "status": line_res.get("status", "EXCLUDED_FEASIBILITY"),
                    "reason": line_res["reason"],
                    "affected_lines": sorted(affected),
                }
            )
            continue

        n_jurkat_dropped_pre_qc += int(line_res["n_jurkat_pre_qc"])
        prov = line_res.get("manifest_provenance")
        for ln in line_res["present_lines"]:
            if ln in excluded_lines or ln not in R15_CANDIDATE_SET:
                continue
            files_by_line[ln].append(hp)
            if prov is not None:
                rec = dict(prov)
                rec["identity_mode"] = line_res["mode"]
                provenance_by_line[ln].append(rec)

    clean_files_by_line = {ln: files for ln, files in sorted(files_by_line.items()) if ln not in excluded_lines}
    clean_prov_by_line = {ln: provs for ln, provs in sorted(provenance_by_line.items()) if ln not in excluded_lines}
    all_prov: List[Dict[str, Any]] = []
    for provs in clean_prov_by_line.values():
        all_prov.extend(provs)

    return {
        "files_by_line": clean_files_by_line,
        "excluded_lines": excluded_lines,
        "provenance_by_line": clean_prov_by_line,
        "all_manifest_provenance": all_prov,
        "unresolved_files": unresolved_files,
        "n_jurkat_dropped_pre_qc": int(n_jurkat_dropped_pre_qc),
    }


def format_r15_manifest_provenance_md(provenance_records: Sequence[Mapping[str, Any]]) -> str:
    """Format rule `(c)` URL + verbatim statement records cleanly for `R15_QC_COUNTS.md`."""
    if not provenance_records:
        return "- None (`line_identity` resolved solely via `.obs` column `(a)` or `(b)`)."
    lines = ["| File | Line | Identity Mode | Source URL | Exact Statement / Manifest Entry |", "| --- | --- | --- | --- | --- |"]
    for rec in provenance_records:
        fn = str(rec.get("file_name") or Path(str(rec.get("file", ""))).name)
        ln = str(rec.get("line", ""))
        mode = str(rec.get("identity_mode", "OFFICIAL_LINE_MANIFEST"))
        url = str(rec.get("url", "")).replace("|", "\\|")
        quote = str(rec.get("quote", "")).replace("|", "\\|").replace("\n", " ")
        lines.append(f"| `{fn}` | `{ln}` | `{mode}` | {url} | {quote} |")
    return "\n".join(lines)


def standardize_vcc_line_to_audit_line(
    line_name: str,
    h5_paths: Sequence[Path],
    cache_out_root: Path,
    aligned_genes: Optional[Sequence[str]] = None,
    hgnc_tsv: Path = DEFAULT_HGNC_TSV,
    context_to_line_map: Optional[Dict[str, str]] = None,
    line_manifest: Optional[Mapping[str, Any]] = None,
    min_gene_overlap: int = 100,
    n_workers: int = 2,
    phase: str = "synthetic_test",
    status_dir: Path = R15_STATUS_DIR,
    seed: int = DEFAULT_R15_SEED,
) -> Tuple[AuditLineData, Dict[str, Any]]:
    """End-to-end clean-room ingestion + `process_r10_line` (`chunk_size=8192`) for VCC lines."""
    enforce_r15_gate(
        phase=phase,
        target_lines=(line_name,),
        status_dir=status_dir,
        extra_dirs=tuple(h5_paths) + (cache_out_root,),
        require_confirmatory_file=False,
    )
    hgnc_prep = build_hgnc_lookup(hgnc_tsv)
    parts, meta = build_vcc_line_audit_parts(
        line_name=line_name,
        h5_paths=h5_paths,
        hgnc=hgnc_prep,
        context_to_line_map=context_to_line_map,
        line_manifest=line_manifest,
        seed=seed,
    )
    line_dir = process_r10_line(
        line_name=_normalize_line_token(line_name),
        parts=parts,
        hgnc=hgnc_prep,
        chunk_size=8192,
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


def resolve_r15_line_cache_dir(line_name: str, phase: str) -> Path:
    """Resolve standardized cache directory cleanly under the active Round 15 blinding phase."""
    norm = _normalize_line_token(line_name)
    if norm in R15_EXCLUDED_LINES:
        raise RuntimeError("PREREG_H §3: jurkat is strictly excluded.")
    if phase == "synthetic_test" or phase.startswith("dev_"):
        if line_name in R9_LINE_SET:
            return R9_LINES_DIR / line_name
        if line_name in R10_ALL_LINES:
            return R10_LINES_DIR / line_name
        if line_name in R11_ALL_LINES:
            return R11_LINES_DIR / line_name
        if line_name in R12_TARGET_KEYS:
            return R12_LINES_DIR / line_name
        if line_name in R13_TARGET_KEYS:
            return R13_LINES_DIR / line_name
        if line_name in R14_TARGET_KEYS:
            return R14_LINES_DIR / line_name
        raise RuntimeError(f"Phase {phase} cannot resolve Round 15 line {line_name!r}.")
    if line_name in R9_LINE_SET:
        return R9_LINES_DIR / line_name
    if line_name in R10_SOURCE_EXT_NAMES:
        return R10_LINES_DIR / line_name
    if norm in R15_CANDIDATE_SET:
        return R15_LINES_DIR / norm
    raise ValueError(f"Unknown line for Round 15 holdout evaluation: {line_name!r}")
