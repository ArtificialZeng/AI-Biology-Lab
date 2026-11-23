#!/usr/bin/env python3
"""
Unit & End-to-End Integration Tests for Round 15 (`PREREG_H` v2) Modules (`round15/src/`).

Runs EXCLUSIVELY on in-memory / temporary synthetic `.h5ad` & planted cohort fixtures:
  - Zero real VCC 2026 files downloaded or touched
  - Zero frozen files under `round9/`–`round14/` modified
Covers both 2025 & 2026 `.obs` schemas (each WITH & WITHOUT `batch` -> Tier 1 / Tier 2),
`§3.1` multi-line + `jurkat` pre-QC exclusion + context-only enforcement, integer raw UMI
slot checks (`SLOT_AMBIGUOUS`, derived-layer exclusion, `EXCLUDED_FEASIBILITY`), `PREREG_F §3.2`
EM fallback, `SA-H1`, `SA-Panel`, `§8.4` cross-wave IUT aggregation, and `PREREG_H` gate locks.
"""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

R15_SRC_DIR = Path(__file__).resolve().parent
REPO_ROOT = R15_SRC_DIR.parents[1]
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
R10_SCRIPTS_DIR = REPO_ROOT / "round10" / "shared" / "scripts"
R11_SRC_DIR = REPO_ROOT / "round11" / "src"
R12_SRC_DIR = REPO_ROOT / "round12" / "src"
for _p in (str(R15_SRC_DIR), str(R10_SRC_DIR), str(R10_SCRIPTS_DIR), str(R11_SRC_DIR), str(R12_SRC_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from data_loader_r10 import LineDatasetR10, load_standardized_line_r10  # noqa: E402
from prepare_r10_lines import HGNCLookup, build_hgnc_lookup  # noqa: E402
from test_unit_r11 import _make_planted_synthetic_cohort  # noqa: E402

from adjudicate_r15 import (  # noqa: E402
    adjudicate_r15_results,
    aggregate_r15_wave_adjudications,
    render_r15_markdown_report,
    render_r15_overall_markdown_report,
)
from data_loader_r15 import (  # noqa: E402
    ALLOWED_R15_LINES,
    DEFAULT_R15_SEED,
    check_prereg_h_gate,
    load_vcc2026_panel_overlap_targets,
    parse_confirmatory_lines_r15,
    validate_r15_candidate_line,
)
from dataset_readers_r15 import (  # noqa: E402
    build_vcc_parts,
    check_vcc_feasibility,
    prepare_r15_line,
)
from run_eval_r15 import evaluate_r15_cohort  # noqa: E402

_CACHED_HGNC: Optional[HGNCLookup] = None


def _get_hgnc() -> HGNCLookup:
    global _CACHED_HGNC
    if _CACHED_HGNC is None:
        _CACHED_HGNC = build_hgnc_lookup()
    return _CACHED_HGNC


def _pick_hgnc_genes(n_genes: int) -> List[str]:
    hgnc = _get_hgnc()
    picked = [g for g in sorted(hgnc.approved) if g and not g.startswith("MT-") and g != "CAS9"][:n_genes]
    if len(picked) < n_genes:
        raise RuntimeError(f"Expected >= {n_genes} HGNC genes, got {len(picked)}")
    return picked


def _write_sparse_group(grp: h5py.Group, mat: sp.spmatrix, layout: str) -> None:
    m = mat.tocsc() if layout == "csc" else mat.tocsr()
    grp.attrs["encoding-type"] = "csc_matrix" if layout == "csc" else "csr_matrix"
    grp.attrs["shape"] = np.array(m.shape, dtype=np.int64)
    dt = np.int32 if np.issubdtype(mat.dtype, np.integer) else np.float32
    grp.create_dataset("data", data=m.data.astype(dt))
    grp.create_dataset("indices", data=m.indices.astype(np.int32))
    grp.create_dataset("indptr", data=m.indptr.astype(np.int64))


def _write_str_dataset(grp: h5py.Group, name: str, vals: Sequence[str], categorical: bool = False) -> None:
    if not categorical:
        grp.create_dataset(name, data=np.array([str(x).encode("utf-8") for x in vals]))
        return
    uniq: List[str] = []
    pos: Dict[str, int] = {}
    codes = np.empty(len(vals), dtype=np.int32)
    for i, v in enumerate(vals):
        sv = str(v)
        if sv not in pos:
            pos[sv] = len(uniq)
            uniq.append(sv)
        codes[i] = pos[sv]
    cg = grp.create_group(name)
    cg.attrs["encoding-type"] = "categorical"
    cg.create_dataset("categories", data=np.array([u.encode("utf-8") for u in uniq]))
    cg.create_dataset("codes", data=codes)


def _write_synthetic_h5ad(
    path: Path,
    X_mat,
    obs_df: Optional[pd.DataFrame],
    var_names: Sequence[str],
    barcodes: Optional[Sequence[str]] = None,
    feature_types: Optional[Sequence[str]] = None,
    x_layout: str = "csr",
    layers: Optional[Dict[str, tuple]] = None,
    raw_tuple: Optional[tuple] = None,
    categorical_obs: bool = False,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    n_obs = X_mat.shape[0]
    with h5py.File(path, "w") as f:
        if x_layout == "dense":
            arr = X_mat.toarray() if sp.issparse(X_mat) else np.asarray(X_mat)
            f.create_dataset("X", data=arr)
        else:
            xg = f.create_group("X")
            _write_sparse_group(xg, sp.csc_matrix(X_mat) if x_layout == "csc" else sp.csr_matrix(X_mat), x_layout)

        if layers:
            lg = f.create_group("layers")
            for lk, (lm, llay) in layers.items():
                if llay == "dense":
                    lg.create_dataset(lk, data=(lm.toarray() if sp.issparse(lm) else np.asarray(lm)))
                else:
                    lsub = lg.create_group(lk)
                    _write_sparse_group(lsub, sp.csc_matrix(lm) if llay == "csc" else sp.csr_matrix(lm), llay)

        og = f.create_group("obs")
        og.attrs["_index"] = "_index"
        bcs = list(barcodes) if barcodes is not None else [f"{path.stem}_C{i:04d}" for i in range(n_obs)]
        _write_str_dataset(og, "_index", bcs, categorical=False)
        if obs_df is not None:
            for col in obs_df.columns:
                _write_str_dataset(og, col, obs_df[col].tolist(), categorical=categorical_obs)

        vg = f.create_group("var")
        vg.attrs["_index"] = "_index"
        _write_str_dataset(vg, "_index", var_names, categorical=False)
        if feature_types is not None:
            _write_str_dataset(vg, "feature_types", feature_types, categorical=False)

        if raw_tuple is not None:
            rmat, rlay, r_names, r_ftypes = raw_tuple
            rg = f.create_group("raw")
            if rlay == "dense":
                rg.create_dataset("X", data=(rmat.toarray() if sp.issparse(rmat) else np.asarray(rmat)))
            else:
                rxg = rg.create_group("X")
                _write_sparse_group(rxg, sp.csc_matrix(rmat) if rlay == "csc" else sp.csr_matrix(rmat), rlay)
            rvg = rg.create_group("var")
            rvg.attrs["_index"] = "_index"
            _write_str_dataset(rvg, "_index", r_names, categorical=False)
            if r_ftypes is not None:
                _write_str_dataset(rvg, "feature_types", r_ftypes, categorical=False)


class Round15UnitTests(unittest.TestCase):
    def test_1_gate_line_validator_and_panel_overlap_loader(self) -> None:
        self.assertEqual(DEFAULT_R15_SEED, 20261231)
        self.assertEqual(ALLOWED_R15_LINES, ("cal33", "hela", "ac16", "molm13", "u87mg"))

        # Candidate line validator & hard Jurkat rejection
        self.assertEqual(validate_r15_candidate_line("Cal-33"), "cal33")
        self.assertEqual(validate_r15_candidate_line("HeLa"), "hela")
        self.assertEqual(validate_r15_candidate_line("AC16"), "ac16")
        self.assertEqual(validate_r15_candidate_line("MOLM-13"), "molm13")
        self.assertEqual(validate_r15_candidate_line("U-87 MG"), "u87mg")
        for bad_j in ("jurkat", "Jurkat", "jurkat_ctx"):
            if bad_j == "jurkat_ctx":
                with self.assertRaises(ValueError):
                    validate_r15_candidate_line(bad_j)
            else:
                with self.assertRaisesRegex(ValueError, "EXCLUDED_LINE"):
                    validate_r15_candidate_line(bad_j)

        # PREREG_H + CODE_FREEZE_R15 hard gate
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)
            with self.assertRaises(RuntimeError):
                check_prereg_h_gate(p)
            (p / "PREREG_H_COMMIT.txt").write_text("prereg_h_abc123\n", encoding="utf-8")
            with self.assertRaises(RuntimeError):
                check_prereg_h_gate(p)
            (p / "CODE_FREEZE_R15.txt").write_text("freeze_r15_def456\n", encoding="utf-8")
            pr, fr = check_prereg_h_gate(p)
            self.assertEqual(pr, "prereg_h_abc123")
            self.assertEqual(fr, "freeze_r15_def456")

            # Confirmatory manifest parsing (0..5 lines, NONE -> K=0, Jurkat -> ValueError)
            cf = p / "CONF.txt"
            cf.write_text("# header\nCal-33\nHeLa\n", encoding="utf-8")
            self.assertEqual(parse_confirmatory_lines_r15(cf), ["cal33", "hela"])
            cf.write_text("NONE\n", encoding="utf-8")
            self.assertEqual(parse_confirmatory_lines_r15(cf), [])
            cf.write_text("cal33\nJurkat\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "EXCLUDED_LINE"):
                parse_confirmatory_lines_r15(cf)

            # Panel overlap file loader
            miss_res, miss_meta = load_vcc2026_panel_overlap_targets(overlap_file=None)
            self.assertIsNone(miss_res)
            self.assertFalse(miss_meta["available"])
            self.assertEqual(miss_meta["status"], "PANEL_OVERLAP_LIST_NOT_AVAILABLE")
            with self.assertRaises(FileNotFoundError):
                load_vcc2026_panel_overlap_targets(overlap_file=p / "nonexistent.txt")

            genes = _pick_hgnc_genes(5)
            pf = p / "overlap_129.txt"
            pf.write_text("\n".join(genes[:3]) + "\n", encoding="utf-8")
            ov, ov_meta = load_vcc2026_panel_overlap_targets(pf, hgnc_lookup=_get_hgnc())
            self.assertTrue(ov_meta["available"])
            self.assertEqual(ov, set(genes[:3]))

    def test_2_multiline_mixed_with_jurkat_and_context_only_rules(self) -> None:
        rng = np.random.default_rng(1502)
        hgnc = _get_hgnc()
        genes = _pick_hgnc_genes(240)
        g0, g1 = genes[0], genes[1]

        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            # (A) Multi-line mixed .h5ad containing Cal-33, HeLa, and Jurkat (with anomalous Jurkat counts/NTCs!)
            mix_dir = tdp / "mixed"
            rows: List[tuple] = []
            # 16 cal33 NTC, 6 cal33 g0, 6 cal33 g1 (plus 2 low-UMI cal33)
            for _ in range(16):
                rows.append(("Cal-33", "non-targeting", "non-targeting_1|non-targeting_2", "b1", False))
            for _ in range(6):
                rows.append(("cal_33", g0, f"{g0}_g1|{g0}_g2", "b1", False))
            for _ in range(6):
                rows.append(("cal33", g1, f"{g1}_g1|{g1}_g2", "b2", False))
            rows.append(("cal33", g0, f"{g0}_g1", "b1", True))
            # 18 HeLa NTC, 6 HeLa g0 (also batch b1/b2)
            for _ in range(18):
                rows.append(("HeLa", "non-targeting", "non-targeting_1", "b2", False))
            for _ in range(6):
                rows.append(("HeLa", g0, f"{g0}_g1", "b2", False))
            # 25 Jurkat cells with bogus target/NTC -> MUST be dropped PRIOR to ANY QC count!
            for _ in range(25):
                rows.append(("Jurkat", "non-targeting", "non-targeting_1", "b1", False))

            X = rng.integers(2, 5, size=(len(rows), len(genes)), dtype=np.int32)
            for i, (_, _, _, _, is_low) in enumerate(rows):
                if is_low:
                    X[i, :] = 0
                    X[i, :30] = 1  # fails QC for cal33
            obs_df = pd.DataFrame(
                {
                    "cell_line": [r[0] for r in rows],
                    "target_gene": [r[1] for r in rows],
                    "guide_id": [r[2] for r in rows],
                    "batch": [r[3] for r in rows],
                }
            )
            _write_synthetic_h5ad(mix_dir / "mixed_val.h5ad", X, obs_df, genes, x_layout="csr", categorical_obs=True)

            # Jurkat directly rejected
            feas_j = check_vcc_feasibility(mix_dir, line="jurkat", hgnc=hgnc)
            self.assertFalse(feas_j["feasible"])
            self.assertEqual(feas_j["status"], "EXCLUDED_LINE")

            # cal33: Jurkat dropped (25), HeLa excluded (24), only 29 cal33 raw rows checked (16 NTC + 12 tgt, 1 low QC)!
            feas_cal = check_vcc_feasibility(mix_dir, line="cal33", hgnc=hgnc)
            self.assertTrue(feas_cal["feasible"], feas_cal.get("reason"))
            self.assertEqual(feas_cal["n_jurkat_dropped"], 25)
            self.assertEqual(feas_cal["n_other_lines_dropped"], 24)
            self.assertEqual(feas_cal["n_target_line_raw_cells"], 29)
            self.assertEqual(feas_cal["n_obs_ntc_meta"], 16)

            parts_cal, audit_cal = build_vcc_parts(mix_dir, hgnc=hgnc, line="cal33", min_half_nt=8, min_batch_nt=8)
            self.assertEqual(audit_cal["n_jurkat_dropped_pre_qc"], 25)
            self.assertEqual(audit_cal["n_rna_qc_excluded"], 1)  # Jurkat NEVER touches RNA QC
            self.assertEqual(audit_cal["n_qc_kept_cells"], 28)
            self.assertEqual(audit_cal["n_kept_ntc"], 16)
            self.assertEqual(audit_cal["n_kept_unique_targets"], 2)

            # (B) Context-only `.obs` file (NO cell_line/cell_type/line) -> MUST fail EXCLUDED_FEASIBILITY without official map!
            ctx_dir = tdp / "ctx_only"
            ctx_obs = pd.DataFrame(
                {
                    "context": ["Context_Alpha"] * 24 + ["Context_Beta"] * 24 + ["Context_J"] * 12,
                    "target_gene": (["non-targeting"] * 12 + [g0] * 12) * 2 + ["non-targeting"] * 12,
                    "ntc_id": (["NTC_01"] * 12 + [""] * 12) * 2 + ["NTC_99"] * 12,
                }
            )
            X_ctx = rng.integers(2, 5, size=(len(ctx_obs), len(genes)), dtype=np.int32)
            _write_synthetic_h5ad(ctx_dir / "vcc2026_ctrl.h5ad", X_ctx, ctx_obs, genes, x_layout="csr")

            feas_no_map = check_vcc_feasibility(ctx_dir, line="cal33", hgnc=hgnc, context_map=None)
            self.assertFalse(feas_no_map["feasible"])
            self.assertEqual(feas_no_map["status"], "EXCLUDED_FEASIBILITY")
            with self.assertRaisesRegex(RuntimeError, "EXCLUDED_FEASIBILITY"):
                build_vcc_parts(ctx_dir, hgnc=hgnc, line="cal33", context_map=None)

            # With official Arc context_map JSON -> Context_Alpha -> cal33, Context_Beta -> hela, Context_J -> jurkat!
            map_json = tdp / "arc_ctx_map.json"
            map_json.write_text(
                json.dumps(
                    {
                        "source_url": "https://arcinstitute.org/example-official-context-table",
                        "mapping": {"Context_Alpha": "Cal-33", "Context_Beta": "HeLa", "Context_J": "Jurkat"},
                    }
                ),
                encoding="utf-8",
            )
            feas_with_map = check_vcc_feasibility(ctx_dir, line="cal33", hgnc=hgnc, context_map_file=map_json)
            self.assertTrue(feas_with_map["feasible"], feas_with_map.get("reason"))
            self.assertEqual(feas_with_map["n_jurkat_dropped"], 12)
            self.assertEqual(feas_with_map["n_other_lines_dropped"], 24)
            self.assertEqual(feas_with_map["n_target_line_raw_cells"], 24)
            parts_ctx, audit_ctx = build_vcc_parts(
                ctx_dir, hgnc=hgnc, line="cal33", context_map_file=map_json, min_half_nt=5, min_batch_nt=5
            )
            self.assertEqual(audit_ctx["split_tier"], "tier2_stratified_seed")
            self.assertEqual(audit_ctx["n_jurkat_dropped_pre_qc"], 12)
            self.assertEqual(audit_ctx["n_kept_ntc"], 12)
            self.assertEqual(audit_ctx["n_ntc_half_a"], 6)
            self.assertEqual(audit_ctx["n_ntc_half_b"], 6)

            # (C) File lacking BOTH explicit line col AND context col -> EXCLUDED_FEASIBILITY
            no_id_dir = tdp / "no_line_id"
            _write_synthetic_h5ad(
                no_id_dir / "nolabel.h5ad",
                X_ctx[:24],
                pd.DataFrame({"target_gene": ["non-targeting"] * 12 + [g0] * 12, "guide_id": ["ntc_1"] * 12 + [f"{g0}_1"] * 12}),
                genes,
            )
            self.assertEqual(check_vcc_feasibility(no_id_dir, line="cal33", hgnc=hgnc)["status"], "EXCLUDED_FEASIBILITY")

    def test_2b_addendum_a_single_line_files_and_manifest_conflicts(self) -> None:
        rng = np.random.default_rng(1522)
        hgnc = _get_hgnc()
        genes = _pick_hgnc_genes(240)
        g0 = genes[0]

        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            # Single-line files WITHOUT any explicit line column and WITHOUT context column
            per_line_dir = tdp / "per_line_release"
            obs_plain = pd.DataFrame(
                {
                    "target_gene": ["non-targeting"] * 16 + [g0] * 12,
                    "guide_id": ["non-targeting_1"] * 16 + [f"{g0}_1"] * 12,
                }
            )
            X_line = rng.integers(2, 5, size=(len(obs_plain), len(genes)), dtype=np.int32)
            _write_synthetic_h5ad(per_line_dir / "cal33.h5ad", X_line, obs_plain, genes, x_layout="csr")
            _write_synthetic_h5ad(per_line_dir / "ac16_panel.h5ad", X_line, obs_plain, genes, x_layout="csr")
            _write_synthetic_h5ad(per_line_dir / "jurkat_ctx.h5ad", X_line, obs_plain, genes, x_layout="csr")

            # 1) WITHOUT manifest -> EXCLUDED_FEASIBILITY even though filename literally says `cal33.h5ad`!
            f_no_man = check_vcc_feasibility(per_line_dir, line="cal33", hgnc=hgnc, line_manifest=None)
            self.assertFalse(f_no_man["feasible"])
            self.assertEqual(f_no_man["status"], "EXCLUDED_FEASIBILITY")
            with self.assertRaisesRegex(RuntimeError, "EXCLUDED_FEASIBILITY"):
                build_vcc_parts(per_line_dir, hgnc=hgnc, line="cal33", line_manifest=None)

            # 2) WITH official Arc manifest JSON (attaching URL + exact quote) -> ingests cleanly & drops Jurkat!
            man_path = tdp / "arc_manifest.json"
            man_path.write_text(
                json.dumps(
                    {
                        "source_url": "https://arcinstitute.org/vcc2026-release-readme",
                        "files": {
                            "cal33.h5ad": {
                                "line": "Cal-33",
                                "url": "https://arcinstitute.org/vcc2026-release-readme#cal33",
                                "quote": "cal33.h5ad contains strictly Cal-33 perturbation single-cell profiles.",
                            },
                            "ac16_panel.h5ad": {
                                "line": "AC16",
                                "url": "https://arcinstitute.org/vcc2026-release-readme#ac16",
                                "quote": "ac16_panel.h5ad contains strictly AC16 cardiomyocyte profiles.",
                            },
                            "jurkat_ctx.h5ad": {
                                "line": "Jurkat",
                                "url": "https://arcinstitute.org/vcc2026-release-readme#jurkat",
                                "quote": "jurkat_ctx.h5ad contains Jurkat validation cells.",
                            },
                        },
                    }
                ),
                encoding="utf-8",
            )
            f_man = check_vcc_feasibility(per_line_dir, line="cal33", hgnc=hgnc, line_manifest_file=man_path)
            self.assertTrue(f_man["feasible"], f_man.get("reason"))
            self.assertEqual(f_man["n_jurkat_dropped"], 28)
            self.assertEqual(f_man["n_other_lines_dropped"], 28)
            self.assertEqual(f_man["n_target_line_raw_cells"], 28)
            id_meta = f_man["file_slots"]["file_0"]["line_identity"]
            self.assertEqual(id_meta["rule"], "line_manifest:cal33.h5ad")
            self.assertEqual(id_meta["manifest_url"], "https://arcinstitute.org/vcc2026-release-readme#cal33")
            self.assertIn("Cal-33", id_meta["manifest_quote"])

            _, aud_man = build_vcc_parts(per_line_dir, hgnc=hgnc, line="ac16", line_manifest_file=man_path, min_half_nt=6, min_batch_nt=6)
            self.assertEqual(aud_man["n_jurkat_dropped_pre_qc"], 28)
            self.assertEqual(aud_man["n_qc_kept_cells"], 28)

            # 3) Manifest vs explicit cell_line conflict -> EXCLUDED_FEASIBILITY
            conf_exp_dir = tdp / "conflict_explicit"
            obs_conflict = obs_plain.copy()
            obs_conflict["cell_line"] = "HeLa"
            _write_synthetic_h5ad(conf_exp_dir / "cal33.h5ad", X_line, obs_conflict, genes, x_layout="csr")
            f_conf = check_vcc_feasibility(
                conf_exp_dir,
                line="cal33",
                hgnc=hgnc,
                line_manifest={"cal33.h5ad": {"line": "cal33", "url": "https://x", "quote": "Cal-33 only"}},
            )
            self.assertFalse(f_conf["feasible"])
            self.assertEqual(f_conf["status"], "EXCLUDED_FEASIBILITY")
            self.assertIn("Addendum A conflict", f_conf["reason"])

            # 4) Manifest vs mapped context conflict -> EXCLUDED_FEASIBILITY
            conf_ctx_dir = tdp / "conflict_context"
            obs_ctx_conf = obs_plain.copy()
            obs_ctx_conf["context"] = "CTX_B"
            _write_synthetic_h5ad(conf_ctx_dir / "cal33.h5ad", X_line, obs_ctx_conf, genes, x_layout="csr")
            f_ctx_conf = check_vcc_feasibility(
                conf_ctx_dir,
                line="cal33",
                hgnc=hgnc,
                context_map={"CTX_B": "hela"},
                line_manifest={"cal33.h5ad": {"line": "cal33", "url": "https://x", "quote": "Cal-33 only"}},
            )
            self.assertFalse(f_ctx_conf["feasible"])
            self.assertEqual(f_ctx_conf["status"], "EXCLUDED_FEASIBILITY")
            self.assertIn("Addendum A conflict", f_ctx_conf["reason"])

    def test_3_2025_schema_with_and_without_batch_plus_dual_cassettes_and_mt(self) -> None:
        rng = np.random.default_rng(1503)
        hgnc = _get_hgnc()
        genes_base = _pick_hgnc_genes(250)
        var_names = genes_base + ["MT-CO1", "CAS9"]
        g0, g1 = genes_base[0], genes_base[1]

        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            # 1) 2025 schema WITH batch (4 batches: each has 10 NTC -> 20 NTC in odd {b1,b3}, 20 in even {b2,b4})
            #    Tests Tier 1 + within-half min_batch_nt=15 fallback + dual-guide cassettes + MT retention + CAS9 drop
            rdir_w = tdp / "vcc2025_with_batch"
            rows_meta: List[tuple] = []
            for bn in ("ch_1", "ch_2", "ch_3", "ch_4"):
                for k in range(10):
                    rows_meta.append(("Cal-33", "non-targeting", f"non-targeting_{k}_a|non-targeting_{k}_b", bn, True))
                for _ in range(4):
                    rows_meta.append(("Cal-33", g0, f"{g0}_g1|{g0}_g2", bn, False))
                for _ in range(4):
                    rows_meta.append(("Cal-33", g1, f"{g1}_g1|{g1}_g2", bn, False))
            # Bad rows that must be excluded cleanly:
            rows_meta.append(("Cal-33", "safe-targeting", "safe_1", "ch_1", False))
            rows_meta.append(("Cal-33", "non-targeting", "gene_desert_1", "ch_1", False))
            rows_meta.append(("Cal-33", g0, f"{g0}_g1|{g1}_g1", "ch_1", False))  # conflicting dual cassette!
            rows_meta.append(("Cal-33", f"{g0}|{g1}", f"{g0}_g1|{g1}_g1", "ch_2", False))  # multi-target gene col!

            n_obs = len(rows_meta)
            X = np.zeros((n_obs, len(var_names)), dtype=np.int32)
            for r, (_, tg, _, _, is_mt) in enumerate(rows_meta):
                X[r, : len(genes_base)] = rng.integers(2, 5, size=len(genes_base), dtype=np.int32)
                if is_mt:
                    X[r, len(genes_base)] = 650  # ~45% MT-CO1 (must NOT drop!)
                if tg == g0:
                    X[r, :15] += 5
                elif tg == g1:
                    X[r, 15:30] += 5
                X[r, -1] = 20  # CAS9 vector must be stripped

            obs_df = pd.DataFrame(
                {
                    "cell_line": [x[0] for x in rows_meta],
                    "target_gene": [x[1] for x in rows_meta],
                    "guide_id": [x[2] for x in rows_meta],
                    "batch": [x[3] for x in rows_meta],
                }
            )
            _write_synthetic_h5ad(rdir_w / "vcc2025_cal33.h5ad", X, obs_df, var_names, x_layout="csr")

            out_root = tdp / "lines_out"
            out_dir_w = prepare_r15_line("cal33", raw_dir=rdir_w, out_root=out_root, hgnc=hgnc, workers=1)
            audit_w = json.loads((out_dir_w / "r15_prep_audit.json").read_text(encoding="utf-8"))
            self.assertEqual(audit_w["split_tier"], "tier1_batch_channel_parity")
            self.assertEqual(audit_w["odd_half_a_units"], ["ch_1", "ch_3"])
            self.assertEqual(audit_w["even_half_b_units"], ["ch_2", "ch_4"])
            self.assertEqual(audit_w["n_ntc_half_a"], 20)
            self.assertEqual(audit_w["n_ntc_half_b"], 20)
            self.assertEqual(set(audit_w["batches_used"]), {"half_1_pooled", "half_2_pooled"})
            self.assertEqual(audit_w["n_safe_targeting_excluded"], 2)
            self.assertEqual(audit_w["n_multi_or_zero_guide_excluded"], 2)
            self.assertEqual(audit_w["n_kept_ntc"], 40)
            self.assertEqual(audit_w["n_kept_unique_targets"], 2)
            ds_w = load_standardized_line_r10(out_dir_w, min_gene_overlap=50)
            self.assertEqual(len(ds_w.target_names), 2)
            self.assertIn("MT-CO1", ds_w.genes)
            self.assertNotIn("CAS9", ds_w.genes)
            self.assertEqual(len(ds_w.target_ensembl_ids), 2)
            for eid in ds_w.target_ensembl_ids:
                self.assertTrue(eid == "" or eid.startswith("ENSG"), f"Unexpected target_ensembl_id: {eid!r}")
            with np.load(out_dir_w / "pseudobulk_delta.npz", allow_pickle=True) as pb_w:
                self.assertEqual([str(x) for x in pb_w["target_ensembl_ids"]], list(ds_w.target_ensembl_ids))
                self.assertEqual([str(x) for x in pb_w["target_raw_symbols"]], sorted([g0, g1]))
                statuses_w = [str(x) for x in pb_w["target_hgnc_status"]]
                self.assertTrue(all(st in {"approved", "alias", "previous"} for st in statuses_w))
                self.assertFalse(any("_g" in st for st in statuses_w))

            # 2) 2025 schema WITHOUT batch -> triggers Tier 2 seed=20261231 stratified split!
            rdir_nb = tdp / "vcc2025_no_batch"
            obs_nb = obs_df.drop(columns=["batch"])
            bcs = [f"BC_{i:04d}" for i in range(n_obs)]
            _write_synthetic_h5ad(rdir_nb / "vcc2025_cal33_nobatch.h5ad", X, obs_nb, var_names, barcodes=bcs, x_layout="csc")

            out_root_nb = tdp / "lines_out_nb"
            out_dir_nb = prepare_r15_line("cal33", raw_dir=rdir_nb, out_root=out_root_nb, hgnc=hgnc, workers=1)
            audit_nb = json.loads((out_dir_nb / "r15_prep_audit.json").read_text(encoding="utf-8"))
            self.assertEqual(audit_nb["split_tier"], "tier2_stratified_seed")
            self.assertEqual(audit_nb["split_seed"], DEFAULT_R15_SEED)
            self.assertEqual(audit_nb["n_ntc_half_a"], 20)
            self.assertEqual(audit_nb["n_ntc_half_b"], 20)
            self.assertEqual(set(audit_nb["batches_used"]), {"half_1", "half_2"})
            ds_nb = load_standardized_line_r10(out_dir_nb, min_gene_overlap=50)
            self.assertEqual(len(ds_nb.target_names), 2)

    def test_4_2026_schema_with_and_without_batch_and_parity_fallback(self) -> None:
        rng = np.random.default_rng(1504)
        hgnc = _get_hgnc()
        genes = _pick_hgnc_genes(240)
        g0, g1 = genes[0], genes[1]
        ctx_map = {"CTX_M13": "MOLM13", "CTX_AC16": "AC16"}

        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            # Build 2026 Control-Pack `.obs`: target_gene, context, ntc_id (+ optional batch)
            rows_m: List[tuple] = []
            for bi, bname in enumerate(("unit_1", "unit_2")):
                for _ in range(16):
                    rows_m.append(("CTX_M13", "non-targeting", "NTC_ctrl_1", bname))
                for _ in range(6):
                    rows_m.append(("CTX_M13", g0, "", bname))
                for _ in range(6):
                    rows_m.append(("CTX_M13", g1, "", bname))
            # Conflict: target_gene=g0 while ntc_id claims NTC -> MUST be excluded!
            rows_m.append(("CTX_M13", g0, "NTC_ctrl_1", "unit_1"))
            # Second context AC16 rows cleanly ignored when preparing molm13
            for _ in range(10):
                rows_m.append(("CTX_AC16", "non-targeting", "NTC_ctrl_1", "unit_1"))

            X = rng.integers(2, 5, size=(len(rows_m), len(genes)), dtype=np.int32)
            # 4A) 2026 WITH batch (16 + 16 NTCs) -> Tier 1
            dir_wb = tdp / "vcc2026_with_batch"
            obs_wb = pd.DataFrame(
                {
                    "context": [r[0] for r in rows_m],
                    "target_gene": [r[1] for r in rows_m],
                    "ntc_id": [r[2] for r in rows_m],
                    "batch": [r[3] for r in rows_m],
                }
            )
            _write_synthetic_h5ad(dir_wb / "ctrl_2026.h5ad", X, obs_wb, genes)
            out_wb = prepare_r15_line("molm13", raw_dir=dir_wb, out_root=tdp / "out_wb", hgnc=hgnc, context_map=ctx_map, workers=1)
            aud_wb = json.loads((out_wb / "r15_prep_audit.json").read_text(encoding="utf-8"))
            self.assertEqual(aud_wb["split_tier"], "tier1_batch_channel_parity")
            self.assertEqual(aud_wb["n_multi_or_zero_guide_excluded"], 1)  # conflicting ntc_id on g0 rejected
            self.assertEqual(aud_wb["n_kept_ntc"], 32)
            ds_wb = load_standardized_line_r10(out_wb, min_gene_overlap=50)
            self.assertEqual(list(ds_wb.target_names), sorted([g0, g1]))
            for eid in ds_wb.target_ensembl_ids:
                self.assertTrue(eid == "" or eid.startswith("ENSG"), f"Unexpected target_ensembl_id: {eid!r}")
            with np.load(out_wb / "pseudobulk_delta.npz", allow_pickle=True) as pb_wb:
                self.assertEqual([str(x) for x in pb_wb["target_ensembl_ids"]], list(ds_wb.target_ensembl_ids))
                self.assertEqual([str(x) for x in pb_wb["target_raw_symbols"]], sorted([g0, g1]))
                statuses_wb = [str(x) for x in pb_wb["target_hgnc_status"]]
                self.assertTrue(all(st in {"approved", "alias", "previous"} for st in statuses_wb))
                self.assertFalse(any("_g" in st for st in statuses_wb))

            # 4B) 2026 WITH unbalanced batch (28 in unit_1, 4 in unit_2 < 15!) -> auto-falls back to Tier 2 (16 & 16)!
            dir_unbal = tdp / "vcc2026_unbal"
            obs_unbal = obs_wb.copy()
            # Reassign first 12 unit_2 NTCs to unit_1 -> 28 vs 4 NTCs
            ntc_u2 = [i for i, r in enumerate(rows_m) if r[0] == "CTX_M13" and r[1] == "non-targeting" and r[3] == "unit_2"][:12]
            for i in ntc_u2:
                obs_unbal.at[i, "batch"] = "unit_1"
            _write_synthetic_h5ad(dir_unbal / "ctrl_2026_unbal.h5ad", X, obs_unbal, genes)
            _, aud_unbal = build_vcc_parts(dir_unbal, hgnc=hgnc, line="molm13", context_map=ctx_map)
            self.assertEqual(aud_unbal["split_tier"], "tier2_stratified_seed")
            self.assertIn("insufficient_ntc_under_parity", aud_unbal["tier1_fallback_reason"])
            self.assertEqual(aud_unbal["n_ntc_half_a"], 16)
            self.assertEqual(aud_unbal["n_ntc_half_b"], 16)

            # 4C) 2026 WITHOUT batch (`target_gene`, `context`, `ntc_id` ONLY) -> Tier 2 seed=20261231
            dir_nb = tdp / "vcc2026_no_batch"
            obs_nb = obs_wb.drop(columns=["batch"])
            _write_synthetic_h5ad(dir_nb / "ctrl_2026_nb.h5ad", X, obs_nb, genes)
            out_nb = prepare_r15_line("molm13", raw_dir=dir_nb, out_root=tdp / "out_nb", hgnc=hgnc, context_map=ctx_map, workers=1)
            aud_nb = json.loads((out_nb / "r15_prep_audit.json").read_text(encoding="utf-8"))
            self.assertEqual(aud_nb["split_tier"], "tier2_stratified_seed")
            self.assertEqual(aud_nb["tier1_fallback_reason"], "missing_batch_column")
            self.assertEqual(aud_nb["split_seed"], DEFAULT_R15_SEED)
            self.assertEqual(aud_nb["n_ntc_half_a"], 16)
            self.assertEqual(aud_nb["n_ntc_half_b"], 16)
            self.assertEqual(load_standardized_line_r10(out_nb, min_gene_overlap=50).ctrl_cells.shape[0], 32)

    def test_5_integer_slots_slot_ambiguous_derived_exclusion_and_em_fallback(self) -> None:
        rng = np.random.default_rng(1505)
        hgnc = _get_hgnc()
        genes = _pick_hgnc_genes(240)
        g0 = genes[0]

        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            n_obs = 40
            X_int = rng.integers(2, 5, size=(n_obs, len(genes)), dtype=np.int32)
            X_flt = np.log1p(X_int.astype(np.float32))
            obs = pd.DataFrame(
                {
                    "cell_line": ["u87mg"] * n_obs,
                    "target_gene": ["non-targeting"] * 20 + [g0] * 20,
                    "guide_id": ["non-targeting_1"] * 20 + [f"{g0}_1"] * 20,
                }
            )

            # 1) Float /X + integer /layers/counts + excluded /layers/spliced & /layers/lognorm -> selects layers/counts
            d_layers = tdp / "layers_ok"
            _write_synthetic_h5ad(
                d_layers / "u87mg.h5ad",
                X_flt,
                obs,
                genes,
                x_layout="csr",
                layers={
                    "spliced": (X_int + 1, "csr"),
                    "lognorm": (X_flt, "csr"),
                    "counts": (X_int, "csr"),
                },
            )
            f_ok = check_vcc_feasibility(d_layers, line="u87mg", hgnc=hgnc)
            self.assertTrue(f_ok["feasible"], f_ok.get("reason"))
            self.assertEqual(f_ok["file_slots"]["file_0"]["x_slot_used"], "layers/counts")

            # 2) Conflicting integer /X and /raw/X -> SLOT_AMBIGUOUS
            d_amb = tdp / "slot_amb"
            _write_synthetic_h5ad(
                d_amb / "u87mg_amb.h5ad",
                X_int,
                obs,
                genes,
                raw_tuple=(X_int + 3, "csr", genes, None),
            )
            f_amb = check_vcc_feasibility(d_amb, line="u87mg", hgnc=hgnc)
            self.assertFalse(f_amb["feasible"])
            self.assertEqual(f_amb["status"], "SLOT_AMBIGUOUS")
            with self.assertRaisesRegex(RuntimeError, "SLOT_AMBIGUOUS"):
                build_vcc_parts(d_amb, hgnc=hgnc, line="u87mg")

            # 3) Deep fractional entry (at row 35) -> EXCLUDED_FEASIBILITY
            d_flt = tdp / "float_deep"
            X_bad = X_int.astype(np.float32)
            X_bad[35, 10] = 2.5
            _write_synthetic_h5ad(d_flt / "u87mg_bad.h5ad", X_bad, obs, genes)
            f_bad = check_vcc_feasibility(d_flt, line="u87mg", hgnc=hgnc)
            self.assertFalse(f_bad["feasible"])
            self.assertEqual(f_bad["status"], "EXCLUDED_FEASIBILITY")

            # 4) PREREG_F §3.2 EM fallback in mixed Jurkat + U87MG file (no .obs target/guide cols)
            d_em = tdp / "em_fallback"
            guides = ["non-targeting_1", f"{g0}_1"]
            var_em = genes + guides
            ftypes = ["Gene Expression"] * len(genes) + ["CRISPR Guide Capture"] * len(guides)
            # 20 Jurkat rows (ignored) + 16 QC U87MG + 55 ngenes>=200/low-UMI U87MG fit rows (>=50 non-zero/guide)
            n_em = 20 + 16 + 55
            mat = np.zeros((n_em, len(var_em)), dtype=np.int32)
            # Jurkat: strong guide counts that must NOT corrupt U87MG EM fit
            for r in range(20):
                mat[r, : len(genes)] = 3
                mat[r, len(genes)] = 50
            # U87MG QC cells (8 NTC, 8 g0)
            for j in range(16):
                r = 20 + j
                mat[r, : len(genes)] = rng.integers(2, 5, size=len(genes), dtype=np.int32)
                mat[r, len(genes) + (0 if j < 8 else 1)] = 25
            # U87MG fit-only cells (ngenes = 205 >= 200, UMI = 205 < 300)
            for j in range(55):
                r = 36 + j
                mat[r, :205] = 1
                mat[r, len(genes)] = 1 if j < 38 else 26
                mat[r, len(genes) + 1] = 1 if j < 38 else 28
            obs_em = pd.DataFrame({"cell_line": ["Jurkat"] * 20 + ["U-87MG"] * 71})
            _write_synthetic_h5ad(d_em / "u87mg_em.h5ad", mat, obs_em, var_em, feature_types=ftypes, x_layout="csr")
            parts_em, aud_em = build_vcc_parts(d_em, hgnc=hgnc, line="u87mg", min_half_nt=3, min_batch_nt=3)
            self.assertEqual(aud_em["guide_call_mode"], "h5ad_crispr_em_fallback")
            self.assertEqual(aud_em["n_jurkat_dropped_pre_qc"], 20)
            self.assertEqual(aud_em["n_kept_ntc"], 8)
            self.assertEqual(aud_em["n_qc_kept_cells"], 16)
            self.assertNotEqual(aud_em["median_guide_threshold"], 7.0)
            ens_g0, can_g0, raw_g0, st_g0 = parts_em[0].pert_meta[g0]
            self.assertTrue(ens_g0 == "" or ens_g0.startswith("ENSG"))
            self.assertEqual(can_g0, g0)
            self.assertEqual(raw_g0, g0)
            self.assertEqual(st_g0, "approved")

            # Also verify optional guide_lib_file fallback preserves exact (ensg, canon, raw_sym, status) tuple
            glib_csv = tdp / "guide_lib.csv"
            glib_csv.write_text(
                f"id,symbol,class\nnon-targeting_1,non-targeting,ctrl_NT\n{g0}_1,{g0},target\n",
                encoding="utf-8",
            )
            parts_lib, _ = build_vcc_parts(
                d_em, hgnc=hgnc, line="u87mg", guide_lib_file=glib_csv, min_half_nt=3, min_batch_nt=3
            )
            ens_l, can_l, raw_l, st_l = parts_lib[0].pert_meta[g0]
            self.assertTrue(ens_l == "" or ens_l.startswith("ENSG"))
            self.assertEqual((can_l, raw_l, st_l), (g0, g0, "approved"))

    def test_6_end_to_end_eval_sa_h1_sa_panel_and_multiwave_iut(self) -> None:
        # Build synthetic 6-source cohort (explicitly naming one source `vcc2025_h1` for SA-H1 test!)
        src6, tgt_tpl = _make_planted_synthetic_cohort(seed=1506, n_genes=150, n_targets=55, n_sources=6)
        old_keys = sorted(src6.keys())
        src_map = {k: src6[k] for k in old_keys[:-1]}
        h1_ds = src6[old_keys[-1]]
        src_map["vcc2025_h1"] = LineDatasetR10(**{**h1_ds.__dict__, "name": "vcc2025_h1"})

        # Build planted passing lines (cal33, hela, ac16) + 1 failing line (molm13 with scrambled deltas)
        t_cal = LineDatasetR10(**{**tgt_tpl.__dict__, "name": "cal33"})
        t_hela = LineDatasetR10(**{**tgt_tpl.__dict__, "name": "hela"})
        t_ac16 = LineDatasetR10(**{**tgt_tpl.__dict__, "name": "ac16"})

        rng_fail = np.random.default_rng(999)
        t_fail = LineDatasetR10(
            **{
                **tgt_tpl.__dict__,
                "name": "molm13",
                "deltas": rng_fail.normal(0.0, 0.25, size=tgt_tpl.deltas.shape),
                "half_a_deltas": rng_fail.normal(0.0, 0.25, size=tgt_tpl.half_a_deltas.shape),
                "half_b_deltas": rng_fail.normal(0.0, 0.25, size=tgt_tpl.half_b_deltas.shape),
            }
        )

        # Panel overlap split of eval targets (129-panel mock: first 28 targets in overlap, rest remainder)
        ov_set = set(tgt_tpl.target_names[:28])

        # Wave 1: cal33 + hela (K_conf_w1 = 2, PASS, pending={ac16, molm13, u87mg} -> "interim, more lines pending")
        w1 = evaluate_r15_cohort(
            target_ds_map={"cal33": t_cal, "hela": t_hela},
            source_datasets=src_map,
            confirmatory_lines=["cal33", "hela"],
            panel_overlap_targets=ov_set,
            panel_overlap_meta={"available": True, "status": "OK", "path": "syn_129.txt", "n_targets": len(ov_set)},
            phase="holdout_r15",
            wave_id="wave1_val",
            min_gene_overlap=100,
            n_boot=120,
            n_perm=10,
            seed=DEFAULT_R15_SEED,
        )
        adj_w1 = w1["adj"]
        self.assertEqual(adj_w1["k_conf"], 2)
        self.assertEqual(adj_w1["statement_label"], "interim, more lines pending")
        self.assertEqual(adj_w1["pending_lines"], ["ac16", "molm13", "u87mg"])
        self.assertEqual(adj_w1["Group_E_Claim1_Verdict"], "CONFIRMED")
        self.assertTrue(adj_w1["single_lab_platform_notice"])
        # Verify SA-H1 ran on 5 remaining sources excluding vcc2025_h1
        self.assertIn("cal33", w1["summary"]["sa_h1"])
        self.assertEqual(len(w1["summary"]["sa_h1"]["cal33"]["remaining_sources"]), 5)
        self.assertNotIn("vcc2025_h1", w1["summary"]["sa_h1"]["cal33"]["remaining_sources"])
        # Verify SA-Panel partitioned targets into overlap vs remainder
        self.assertEqual(w1["summary"]["sa_panel"]["cal33"]["in_val_test_overlap_panel"]["n_targets"], 28)
        self.assertEqual(w1["summary"]["sa_panel"]["cal33"]["panel_unique_remainder"]["n_targets"], 27)
        md_w1 = render_r15_markdown_report(w1["df_lines"], w1["df_boot"], w1["summary"], adj_w1, "syn_w1")
        self.assertIn("interim, more lines pending", md_w1)
        self.assertIn("SA-H1", md_w1)
        self.assertIn("SA-Panel", md_w1)
        self.assertIn("10x Flex CRISPRi", md_w1)

        # Wave 2A: ac16 PASS + u87mg excluded -> all 4 non-molm + let's make molm13 PASS for full pass vs FAIL for revocation!
        w2_pass = evaluate_r15_cohort(
            target_ds_map={"ac16": t_ac16},
            source_datasets=src_map,
            confirmatory_lines=["ac16"],
            excluded_lines=["molm13", "u87mg"],
            phase="holdout_r15",
            wave_id="wave2_test_pass",
            min_gene_overlap=100,
            n_boot=120,
            n_perm=10,
            seed=DEFAULT_R15_SEED,
        )
        overall_ok = aggregate_r15_wave_adjudications([adj_w1, w2_pass["adj"]])
        self.assertEqual(overall_ok["k_conf_total"], 3)
        self.assertEqual(overall_ok["Overall_Claim1_Verdict"], "CONFIRMED")
        self.assertEqual(overall_ok["statement_label"], "final")
        self.assertEqual(overall_ok["pending_lines"], [])
        md_ov_ok = render_r15_overall_markdown_report(overall_ok, "syn_overall")
        self.assertIn("10x Flex CRISPRi", md_ov_ok)

        # Wave 2B: molm13 FAILS -> overall Claim 1 MUST be REVOKED to NOT_CONFIRMED across all waves!
        w2_fail = evaluate_r15_cohort(
            target_ds_map={"molm13": t_fail},
            source_datasets=src_map,
            confirmatory_lines=["molm13"],
            excluded_lines=["ac16", "u87mg"],
            phase="holdout_r15",
            wave_id="wave2_test_fail",
            min_gene_overlap=100,
            n_boot=120,
            n_perm=10,
            seed=DEFAULT_R15_SEED,
        )
        self.assertEqual(w2_fail["adj"]["Group_E_Claim1_Verdict"], "NOT_CONFIRMED")
        overall_revoked = aggregate_r15_wave_adjudications([adj_w1, w2_fail["adj"]])
        self.assertEqual(overall_revoked["Overall_Claim1_Verdict"], "NOT_CONFIRMED")
        self.assertEqual(overall_revoked["H-E1"]["pass_count"], 2)
        self.assertEqual(overall_revoked["H-E1"]["required"], 3)
        # Reversed order (Wave 1 fails, Wave 2 passes) -> Wave 2 CANNOT rescue Wave 1!
        overall_no_rescue = aggregate_r15_wave_adjudications([w2_fail["adj"], adj_w1])
        self.assertEqual(overall_no_rescue["Overall_Claim1_Verdict"], "NOT_CONFIRMED")
        self.assertEqual(overall_no_rescue["wave_trajectory"][0]["cumulative_Claim1"], "NOT_CONFIRMED")
        self.assertEqual(overall_no_rescue["wave_trajectory"][1]["cumulative_Claim1"], "NOT_CONFIRMED")

        # K_conf = 0 and K_conf = 1 checks + sunset clause
        w0 = evaluate_r15_cohort(
            target_ds_map={"cal33": t_cal},
            source_datasets=src_map,
            confirmatory_lines=[],
            phase="holdout_r15",
            sunset_reached=True,
            min_gene_overlap=100,
            n_boot=120,
            n_perm=10,
            seed=DEFAULT_R15_SEED,
        )
        self.assertEqual(w0["adj"]["Group_E_Claim1_Verdict"], "CONFIRMATORY_ABORTED_K0")
        self.assertEqual(w0["adj"]["statement_label"], "final")

        w1_single = adjudicate_r15_results(
            w1["df_lines"], w1["df_boot"], w1["summary"], confirmatory_lines=["cal33"], wave_id="single_line"
        )
        self.assertEqual(w1_single["confirmatory_status"], "CONFIRMATORY_SINGLE_LINE_K1")
        self.assertTrue(w1_single["single_line_limited_generalizability_notice"])

        # CLI gate block check on `--phase holdout_r15` without gate files
        with tempfile.TemporaryDirectory() as td_cli:
            proc = subprocess.run(
                [
                    sys.executable,
                    "-I",
                    str(R15_SRC_DIR / "run_eval_r15.py"),
                    "--phase",
                    "holdout_r15",
                    "--gate-dir",
                    td_cli,
                ],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("PREREG_H hard gate blocked", proc.stderr + proc.stdout)

            # CLI cross-wave aggregator check (`adjudicate_r15.py --wave-jsons ...`)
            w1_j = Path(td_cli) / "w1.json"
            w2_j = Path(td_cli) / "w2.json"
            ov_j = Path(td_cli) / "overall.json"
            ov_md = Path(td_cli) / "overall.md"
            w1_j.write_text(json.dumps(adj_w1, indent=2), encoding="utf-8")
            w2_j.write_text(json.dumps(w2_pass["adj"], indent=2), encoding="utf-8")
            agg_proc = subprocess.run(
                [
                    sys.executable,
                    "-I",
                    str(R15_SRC_DIR / "adjudicate_r15.py"),
                    "--wave-jsons",
                    str(w1_j),
                    str(w2_j),
                    "--out-json",
                    str(ov_j),
                    "--out-md",
                    str(ov_md),
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(agg_proc.returncode, 0, agg_proc.stderr)
            ov_parsed = json.loads(ov_j.read_text(encoding="utf-8"))
            self.assertEqual(ov_parsed["Overall_Claim1_Verdict"], "CONFIRMED")
            self.assertEqual(ov_parsed["statement_label"], "final")
            self.assertIn("10x Flex CRISPRi", ov_md.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
