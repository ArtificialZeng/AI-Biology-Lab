#!/usr/bin/env python3
"""
Synthetic Fixture Unit & Cross-Check Suite for Round 15 (`PREREG_H v2 + Addendum A`, A-T3).

Strictly uses temporary synthetic `.h5ad` fixtures outside `round15/shared/`:
  - Verifies `§1–§2` open-license/sunset (`2027-06-30`) & `§3` hard `jurkat` exclusion.
  - Verifies `§3.1` multi-line `.h5ad` cell-line identity resolution:
      * Explicit `cell_line`/`cell_type`/`line` column (`Cal-33`, `HE_LA`, `Jurkat`) with
        pre-QC `Jurkat` removal (never enters any QC denominator or matrix depth calculation).
      * `context`-only `.h5ad` WITHOUT `cell_line`: strictly fails `EXCLUDED_FEASIBILITY`
        unless paired with an explicit `context_to_line_map` or official `line_manifest`.
  - Verifies `Addendum A` (`§3.1 rule (c)` one-file-per-line without `cell_line` column):
      * Valid official `line_manifest` (`line`, `url`, exact `quote`) -> ingests cleanly,
        records URL + quote provenance, and matches independent numpy oracle (`<= 1e-6`).
      * Missing `line_manifest` (or bare string / missing URL / missing quote) -> `EXCLUDED_FEASIBILITY`.
      * Disagreement between `(a)` explicit column or `(b)` mapped `context` and `(c)` `line_manifest` ->
        ALL affected lines marked `EXCLUDED_FEASIBILITY` across batch partition and single-line check.
  - Verifies `§5` multi-integer slot inconsistency stop guard, dual `2025`/`2026` `.obs` schemas,
    single-construct dual-guide cassettes vs multi-target rejection, `PREREG_F §3.2` `/obsm` fallback,
    and `max_pct_mt = None`.
  - Verifies `§6` Tier 1 (`.obs` batch parity) & Tier 2 (`seed = 20261231`), oracle array parity
    (`max_abs_diff <= 1e-6`), and `§8.4` release-wave cross-IUT aggregation + `SA-H1` & `SA-Panel`.
"""

from __future__ import annotations

import datetime as _dt
import sys
import tempfile
import unittest
from pathlib import Path
from typing import List

import numpy as np

AUDIT_R15_DIR = Path(__file__).resolve().parent
REPO_ROOT = AUDIT_R15_DIR.parents[1]
R14_AUDIT_DIR = REPO_ROOT / "round14" / "audit"
for _p in (
    str(AUDIT_R15_DIR),
    str(R14_AUDIT_DIR),
    str(REPO_ROOT / "round10" / "shared" / "scripts"),
    str(REPO_ROOT / "round11" / "audit"),
    str(REPO_ROOT / "round12" / "audit"),
):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from independent_r11_loader import AuditLineData, DEFAULT_HGNC_TSV  # noqa: E402
from independent_r12_metrics import evaluate_target_line_r12  # noqa: E402
from independent_r15_loader import (  # noqa: E402
    DEFAULT_R15_SEED,
    R15_SHARED_DIR,
    build_vcc_line_audit_parts,
    check_vcc_h5ad_feasibility,
    enforce_r15_gate,
    format_r15_manifest_provenance_md,
    parse_confirmatory_lines_r15,
    partition_vcc_files_by_line,
    standardize_vcc_line_to_audit_line,
    verify_r15_open_release_and_sunset,
)
from prepare_r10_lines import build_hgnc_lookup  # noqa: E402
from run_r15_audit import (  # noqa: E402
    adjudicate_overall_round15,
    adjudicate_release_wave,
    check_r15_line_thresholds,
    evaluate_sa_h1_leave_vcc2025_out,
    evaluate_sa_panel_overlap,
)
from test_r14_audit import _get_real_gene_pool, _write_synthetic_h5ad  # noqa: E402
from xcheck_r15_synthetic import compare_line_dirs_exact, compute_independent_oracle_arrays  # noqa: E402


class TestRound15Audit(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.hgnc = build_hgnc_lookup(DEFAULT_HGNC_TSV)
        cls.gene_pool = _get_real_gene_pool(240)

    def test_gates_open_license_sunset_and_jurkat_exclusion(self) -> None:
        # License & sunset gate (`2027-06-30`)
        self.assertEqual(
            verify_r15_open_release_and_sunset("CC0-1.0", requires_auth_or_terms=True, open_released=True)["status"],
            "EXCLUDED_ACCESS",
        )
        self.assertEqual(
            verify_r15_open_release_and_sunset(None, open_released=False, check_date=_dt.date(2027, 7, 1))["status"],
            "PREREG_H_EXPIRED",
        )
        self.assertTrue(
            verify_r15_open_release_and_sunset("CC-BY-4.0", requires_auth_or_terms=False, open_released=True)["allowed"]
        )

        # Jurkat blocked in every phase; round15/shared blocked in synthetic_test
        with self.assertRaises(RuntimeError):
            enforce_r15_gate("synthetic_test", target_lines=["jurkat"])
        with self.assertRaises(RuntimeError):
            enforce_r15_gate("synthetic_test", extra_dirs=[R15_SHARED_DIR / "raw"])
        with self.assertRaises(RuntimeError):
            enforce_r15_gate("holdout_r15", target_lines=["ac16"])

        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            conf_bad = tdir / "CONF_BAD.txt"
            conf_bad.write_text("CONFIRMATORY_LINES=cal33,jurkat\n", encoding="utf-8")
            with self.assertRaises(RuntimeError):
                parse_confirmatory_lines_r15(conf_bad)

            conf_ok = tdir / "CONF_OK.txt"
            conf_ok.write_text("CONFIRMATORY_LINES=cal33,ac16\nSECONDARY_LINES=hela\nU87MG_ROLE=EXCLUDED_FEASIBILITY\n", encoding="utf-8")
            p = parse_confirmatory_lines_r15(conf_ok)
            self.assertEqual(p["confirmatory_lines"], ["cal33", "ac16"])
            self.assertEqual(p["secondary_eval_lines"], ["hela"])
            self.assertEqual(p["role_by_line"]["u87mg"], "EXCLUDED_FEASIBILITY")

    def test_v2_3_1_multiline_jurkat_preqc_and_context_only_files(self) -> None:
        genes = self.gene_pool[:215] + [("ENSG00000198804", "MT-CO1"), ("CAS9", "Cas9")]
        v_ids = [g[0] for g in genes]
        v_syms = [g[1] for g in genes]
        rng = np.random.default_rng(1501)

        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            h5_mix = tdir / "mixed_multiline_with_jurkat.h5ad"

            # Build mixed file: Cal-33 (36 NTC, 20 HDAC3, 1 low-UMI, 1 high-%MT HDAC3), HE_LA (10), Jurkat (30 + 5 ultra-low UMI)
            obs_line: List[str] = []
            obs_tgt: List[str] = []
            obs_gid: List[str] = []
            obs_batch: List[str] = []
            rows_mat: List[np.ndarray] = []

            def _add(cl: str, tgt: str, gid: str, b: str, n_expr: int = 210, extra: int = 3, mt_u: int = 5) -> None:
                v = np.zeros(len(genes), dtype=np.int32)
                v[:n_expr] = extra
                v[215] = mt_u
                v[216] = 20  # Cas9 stripped before QC
                obs_line.append(cl)
                obs_tgt.append(tgt)
                obs_gid.append(gid)
                obs_batch.append(b)
                rows_mat.append(v)

            for i in range(36):
                _add("Cal-33", "non-targeting", f"non-targeting_A{(i%4)+1}|non-targeting_B{(i%4)+1}", f"Flex_0{(i%4)+1}")
            for i in range(20):
                _add("cal_33", "HDAC3", "HDAC3_1|HDAC3_2", f"Flex_0{(i%4)+1}")
            # High %MT Cal-33 cell (>70% MT, >=300 human UMI, >=200 genes -> MUST BE KEPT!)
            _add("CAL33", "HDAC3", "HDAC3_1", "Flex_01", n_expr=205, extra=2, mt_u=1200)
            # Low-UMI Cal-33 cell (120 genes -> <200 genes & <300 UMI -> MUST fail QC)
            _add("Cal-33", "HDAC3", "HDAC3_1", "Flex_02", n_expr=120, extra=1, mt_u=0)
            # Multi-target doublet on Cal-33 -> MUST exclude
            _add("Cal-33", "HDAC3,PTEN", "HDAC3_1|PTEN_1", "Flex_01")
            # HeLa cells
            for i in range(10):
                _add("HE_LA", "non-targeting" if i < 6 else "PTEN", "g", f"Flex_0{(i%2)+1}")
            # Jurkat cells (valid + invalid ultra-low UMI -> MUST BE DROPPED PRE-QC before ANY count)
            for i in range(25):
                _add("Jurkat", "non-targeting" if i < 15 else "CREBBP", "gJ", f"Flex_0{(i%2)+1}")
            for _ in range(7):
                _add("jurkat", "CREBBP", "gJ", "Flex_01", n_expr=50, extra=1, mt_u=0)

            _write_synthetic_h5ad(
                h5_mix,
                np.vstack(rows_mat),
                [f"M_{i:04d}" for i in range(len(obs_line))],
                v_ids,
                v_syms,
                obs_cols={
                    "cell_line": obs_line,
                    "target_gene": obs_tgt,
                    "guide_id": obs_gid,
                    "batch": obs_batch,
                },
                categorical_cols=("cell_line", "batch"),
            )

            # Jurkat directly blocked even on feasibility
            fj = check_vcc_h5ad_feasibility("jurkat", [h5_mix], hgnc=self.hgnc)
            self.assertFalse(fj["feasible"])
            self.assertEqual(fj["status"], "EXCLUDED_LINE_JURKAT")

            # Cal-33 extracts ONLY Cal-33 rows, purges all 32 Jurkat rows PRE-QC, counts ONLY the 1 Cal-33 QC failure!
            parts_c, meta_c = build_vcc_line_audit_parts("cal33", [h5_mix], hgnc=self.hgnc)
            self.assertEqual(meta_c["split_tier"], "TIER1_OBS_BATCH_PARITY")
            self.assertEqual(meta_c["counts"]["n_jurkat_dropped_pre_qc"], 32)
            # Cal-33 total raw rows = 36 NTC + 20 HDAC3 + 1 hiMT + 1 lowUMI + 1 doublet = 59 (NOT counting Jurkat or HeLa!)
            self.assertEqual(meta_c["counts"]["n_line_raw_cells"], 59)
            self.assertEqual(meta_c["counts"]["n_failed_qc_umi_or_genes"], 1)
            self.assertEqual(meta_c["counts"]["n_kept_ntc"], 36)
            self.assertEqual(meta_c["counts"]["n_kept_targeting_cells"], 21)
            self.assertEqual(meta_c["n_nt_half_a"], 18)
            self.assertEqual(meta_c["n_nt_half_b"], 18)

            # Ambiguous cell_line value -> EXCLUDED_FEASIBILITY
            h5_amb = tdir / "ambig_line.h5ad"
            _write_synthetic_h5ad(
                h5_amb,
                np.vstack(rows_mat[:10]),
                [f"A_{i}" for i in range(10)],
                v_ids,
                v_syms,
                obs_cols={"cell_line": ["Cal-33"] * 8 + ["MysteryLine"] * 2, "target_gene": ["non-targeting"] * 10},
                categorical_cols=(),
            )
            fa = check_vcc_h5ad_feasibility("cal33", [h5_amb], hgnc=self.hgnc)
            self.assertFalse(fa["feasible"])
            self.assertEqual(fa["status"], "EXCLUDED_FEASIBILITY")

            # Context-only file (2026 controls schema WITHOUT cell_line and WITHOUT batch)
            h5_ctx = tdir / "context_only_no_batch.h5ad"
            n_ctx = 60
            x_ctx = rng.poisson(6.0, size=(n_ctx, len(genes))).astype(np.int32)
            ctx_col = ["Context_D"] * 44 + ["Context_C"] * 16  # Context_C -> jurkat via official map
            tgt_col = ["non-targeting"] * 30 + ["HDAC3"] * 14 + ["non-targeting"] * 10 + ["HDAC3"] * 6
            ntc_col = [f"ntc_{(i%46)+1}" if i < 30 else "NA" for i in range(44)] + ["ntc_1"] * 10 + ["NA"] * 6
            _write_synthetic_h5ad(
                h5_ctx,
                x_ctx,
                [f"C26_{i:04d}" for i in range(n_ctx)],
                v_ids,
                v_syms,
                obs_cols={"context": ctx_col, "target_gene": tgt_col, "ntc_id": ntc_col},
                categorical_cols=("context",),
            )

            # Without context_to_line_map -> MUST fail EXCLUDED_FEASIBILITY!
            f_nomap = check_vcc_h5ad_feasibility("ac16", [h5_ctx], hgnc=self.hgnc, context_to_line_map=None)
            self.assertFalse(f_nomap["feasible"])
            self.assertEqual(f_nomap["status"], "EXCLUDED_FEASIBILITY")
            self.assertIn("context_only_obs", f_nomap["reason"])

            # With official context_to_line_map -> succeeds, drops Jurkat (16 rows) pre-QC, uses Tier 2 seed=20261231!
            cmap = {"Context_D": "AC16", "Context_C": "Jurkat"}
            parts_ctx, meta_ctx = build_vcc_line_audit_parts("ac16", [h5_ctx], hgnc=self.hgnc, context_to_line_map=cmap)
            self.assertEqual(meta_ctx["counts"]["n_jurkat_dropped_pre_qc"], 16)
            self.assertEqual(meta_ctx["counts"]["n_line_raw_cells"], 44)
            self.assertEqual(meta_ctx["split_tier"], "TIER2_STRATIFIED_SEED20261231")
            self.assertEqual(meta_ctx["n_nt_half_a"], 15)
            self.assertEqual(meta_ctx["n_nt_half_b"], 15)
            self.assertEqual(set(parts_ctx[0].batches), {"half_1", "half_2"})

    def test_slot_conflict_guard_dual_schema_tier2_and_oracle_parity(self) -> None:
        genes = self.gene_pool[:225]
        v_ids = [g[0] for g in genes]
        v_syms = [g[1] for g in genes]
        rng = np.random.default_rng(1502)

        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            n_c = 70
            x_base = rng.poisson(7.0, size=(n_c, len(genes))).astype(np.int32)
            # Conflict file: /X and /layers/counts differ by +5 in one cell
            h5_bad_slots = tdir / "slot_conflict.h5ad"
            x_diff = x_base.copy()
            x_diff[3, 10] += 5
            _write_synthetic_h5ad(
                h5_bad_slots,
                x_base,
                [f"S_{i}" for i in range(n_c)],
                v_ids,
                v_syms,
                obs_cols={"cell_line": ["U-87MG"] * n_c, "target_gene": ["non-targeting"] * 40 + ["PTEN"] * 30},
                layers={"counts": (x_diff, "csr", False)},
            )
            fc = check_vcc_h5ad_feasibility("u87mg", [h5_bad_slots], hgnc=self.hgnc)
            self.assertFalse(fc["feasible"])
            self.assertIn("MULTIPLE_INCONSISTENT_INTEGER_SLOTS", fc["reason"])

            # Valid 2026-style file (explicit cell_line, context, ntc_id, NO batch -> Tier 2) with /raw/X == /X
            h5_ok = tdir / "u87mg_tier2.h5ad"
            obs_u = {
                "cell_line": ["U-87MG"] * n_c,
                "context": ["D"] * n_c,
                "target_gene": ["non-targeting"] * 40 + ["HDAC3"] * 15 + ["PTEN"] * 15,
                "ntc_id": [f"ntc_{(i%10)+1}" for i in range(40)] + ["NA"] * 30,
            }
            _write_synthetic_h5ad(
                h5_ok,
                x_base,
                [f"U_{i:04d}" for i in range(n_c)],
                v_ids,
                v_syms,
                obs_cols=obs_u,
                raw_x_mat=x_base,
                raw_var_ids=v_ids,
                raw_var_syms=v_syms,
                categorical_cols=("cell_line", "context"),
            )
            out_cache = tdir / "cache"
            adata, meta = standardize_vcc_line_to_audit_line(
                "u87mg",
                [h5_ok],
                cache_out_root=out_cache,
                aligned_genes=v_syms,
                min_gene_overlap=100,
                n_workers=1,
                phase="synthetic_test",
                seed=DEFAULT_R15_SEED,
            )
            self.assertEqual(meta["split_tier"], "TIER2_STRATIFIED_SEED20261231")
            self.assertEqual(meta["n_nt_half_a"], 20)
            self.assertEqual(meta["n_nt_half_b"], 20)
            self.assertEqual(sorted(adata.target_names), ["HDAC3", "PTEN"])

            # Compare built cache against pure-numpy independent oracle AND self-cache comparator (`<= 1e-6`)
            oracle = compute_independent_oracle_arrays([h5_ok], "u87mg", seed=DEFAULT_R15_SEED)
            line_dir = out_cache / "u87mg"
            with np.load(line_dir / "ctrl_cells.npz", allow_pickle=False) as c_r, np.load(
                line_dir / "pseudobulk_delta.npz", allow_pickle=False
            ) as p_r:
                # On a 225-gene fixture ctrl_mean_libnorm ~ 44.4 (1 float32 ULP = 2^-18 = 3.81e-6): check <= 1 ULP + rel <= 1e-6
                ulp_bound = float(2.0 * np.spacing(np.float32(np.max(c_r["ctrl_mean_libnorm"]))))
                self.assertLessEqual(
                    float(np.max(np.abs(c_r["ctrl_mean_libnorm"] - oracle["ctrl_mean_libnorm"]))),
                    ulp_bound,
                )
                self.assertLessEqual(
                    float(
                        np.max(
                            np.abs(c_r["ctrl_mean_libnorm"] - oracle["ctrl_mean_libnorm"])
                            / np.maximum(np.abs(oracle["ctrl_mean_libnorm"]), 1e-12)
                        )
                    ),
                    1e-6,
                )
                self.assertLessEqual(float(np.max(np.abs(p_r["delta_log1p"] - oracle["delta_log1p"]))), 1e-6)
                self.assertLessEqual(float(np.max(np.abs(p_r["half_a_delta_log1p"] - oracle["half_a_delta_log1p"]))), 1e-6)
                self.assertLessEqual(float(np.max(np.abs(p_r["half_b_delta_log1p"] - oracle["half_b_delta_log1p"]))), 1e-6)
            cmp_self = compare_line_dirs_exact(line_dir, line_dir, tol=1e-6)
            self.assertTrue(cmp_self["pass"])

            # Cross-compare /raw/X CSC-only vs /X CSR-only builds of identical counts: strictly 0.0 diff!
            h5_csr = tdir / "u87mg_csr_only.h5ad"
            _write_synthetic_h5ad(
                h5_csr,
                x_base,
                [f"U_{i:04d}" for i in range(n_c)],
                v_ids,
                v_syms,
                obs_cols=obs_u,
                x_encoding="csr",
                categorical_cols=("cell_line", "context"),
            )
            out_cache2 = tdir / "cache2"
            standardize_vcc_line_to_audit_line(
                "u87mg",
                [h5_csr],
                cache_out_root=out_cache2,
                aligned_genes=v_syms,
                min_gene_overlap=100,
                n_workers=1,
                phase="synthetic_test",
                seed=DEFAULT_R15_SEED,
            )
            cmp_dual = compare_line_dirs_exact(line_dir, out_cache2 / "u87mg", tol=1e-6)
            self.assertTrue(cmp_dual["pass"])
            for k_arr, info_arr in cmp_dual["diffs"].items():
                if "max_abs_diff" in info_arr:
                    self.assertEqual(info_arr["max_abs_diff"], 0.0, msg=f"{k_arr} nonzero diff")

            # Also verify PREREG_F /obsm raw guide fallback works cleanly when target_gene is absent
            h5_obsm = tdir / "molm13_obsm_fallback.h5ad"
            all_g = ["non-targeting_1", "non-targeting_2", "HDAC3_1", "PTEN_1"]
            g_mat = np.zeros((n_c, len(all_g)), dtype=np.int32)
            for i in range(40):
                g_mat[i, i % 2] = 15
            for i in range(40, 55):
                g_mat[i, 2] = 15
            for i in range(55, 70):
                g_mat[i, 3] = 15
            _write_synthetic_h5ad(
                h5_obsm,
                x_base,
                [f"M_{i:04d}" for i in range(n_c)],
                v_ids,
                v_syms,
                obs_cols={"cell_line": ["MOLM13"] * n_c},
                obsm={"sgrna_counts": (g_mat, "csr")},
                uns_guide_names=all_g,
                categorical_cols=("cell_line",),
            )
            parts_m, meta_m = build_vcc_line_audit_parts("molm13", [h5_obsm], hgnc=self.hgnc)
            self.assertEqual(meta_m["counts"]["n_kept_ntc"], 40)
            self.assertEqual(meta_m["counts"]["n_kept_targeting_cells"], 30)

    def test_v2_8_4_release_waves_overall_iut_and_sa_h1_sa_panel(self) -> None:
        # 1) Release-wave cross-IUT logic (per §8.4 & §8.2)
        w1 = adjudicate_release_wave(
            wave_id="wave_1_val",
            per_line_summary={
                "cal33": {"n_eval_targets": 60, "n_genes_eval": 4000, "he1_pass": True, "he2_pass": True},
                "hela": {"n_eval_targets": 40, "n_genes_eval": 4000, "he1_pass": False, "he2_pass": True},
            },
            gate_info={
                "confirmatory_lines": ["cal33", "hela"],
                "secondary_eval_lines": [],
                "role_by_line": {"cal33": "CONFIRMATORY", "hela": "CONFIRMATORY", "ac16": "NOT_RELEASED", "molm13": "NOT_RELEASED", "u87mg": "NOT_RELEASED"},
            },
            check_date=_dt.date(2026, 12, 1),
        )
        self.assertEqual(w1["k_conf_wave"], 1)  # hela auto-demoted (<50)
        self.assertEqual(w1["interim_label"], "interim, more lines pending")
        self.assertFalse(w1["is_final_study_endpoint"])
        self.assertTrue(w1["claim1_confirmed_wave"])

        w2 = adjudicate_release_wave(
            wave_id="wave_2_test",
            per_line_summary={
                "ac16": {"n_eval_targets": 75, "n_genes_eval": 4500, "he1_pass": True, "he2_pass": True},
                "molm13": {"n_eval_targets": 70, "n_genes_eval": 4500, "he1_pass": False, "he2_pass": True},
            },
            gate_info={
                "confirmatory_lines": ["ac16", "molm13"],
                "secondary_eval_lines": [],
                "role_by_line": {"u87mg": "EXCLUDED_FEASIBILITY", "ac16": "CONFIRMATORY", "molm13": "CONFIRMATORY"},
            },
            all_resolved_lines_so_far=["cal33", "hela"],
            check_date=_dt.date(2027, 1, 15),
        )
        overall = adjudicate_overall_round15([w1, w2], check_date=_dt.date(2027, 1, 15))
        self.assertTrue(overall["is_final_study_endpoint"])
        self.assertEqual(overall["interim_label"], "FINAL")
        self.assertEqual(overall["k_conf_total"], 3)
        self.assertEqual(overall["confirmatory_lines"], ["cal33", "ac16", "molm13"])
        self.assertFalse(overall["claim1_confirmed"])  # molm13 revokes wave_1 confirmation
        self.assertTrue(overall["claim2_confirmed"])   # 3/3 across waves
        self.assertIn("Arc Institute VCC 2026", overall["mandatory_scope_disclosure"])

        # Wave 2 CANNOT rescue Wave 1 failure
        w1_fail = adjudicate_release_wave(
            "w1",
            {"cal33": {"n_eval_targets": 60, "n_genes_eval": 4000, "he1_pass": False, "he2_pass": True}},
            {"confirmatory_lines": ["cal33"], "secondary_eval_lines": [], "role_by_line": {"cal33": "CONFIRMATORY"}},
            check_date=_dt.date(2026, 12, 1),
        )
        ov_k1 = adjudicate_overall_round15([w1_fail], check_date=_dt.date(2026, 12, 1))
        self.assertEqual(ov_k1["k_conf_total"], 1)
        self.assertIn("单线，推广性有限", ov_k1["mandatory_scope_disclosure"])
        ov_norescue = adjudicate_overall_round15([w1_fail, w2], check_date=_dt.date(2027, 1, 15))
        self.assertFalse(ov_norescue["claim1_confirmed"])

        # Past sunset with 0 waves -> expired
        ov_exp = adjudicate_overall_round15([], check_date=_dt.date(2027, 7, 1))
        self.assertEqual(ov_exp["status"], "CONFIRMATORY_EXPIRED_UNRELEASED")
        self.assertTrue(ov_exp["is_final_study_endpoint"])

        # 2) SA-H1 (leave vcc2025_h1 out) & SA-Panel (129 overlap vs remaining) on synthetic AuditLineData
        g_names = [g[1] for g in self.gene_pool[:220]]
        t_names = g_names[:8]
        rng = np.random.default_rng(DEFAULT_R15_SEED)
        mu = rng.uniform(2.0, 20.0, size=len(g_names))

        def _fake_line(name: str, scale: float) -> AuditLineData:
            deltas = rng.normal(0.0, 0.3, size=(len(t_names), len(g_names)))
            ha = deltas + rng.normal(0.0, 0.05, size=deltas.shape)
            hb = deltas + rng.normal(0.0, 0.05, size=deltas.shape)
            return AuditLineData(
                name=name,
                genes=list(g_names),
                target_names=list(t_names),
                target_indices=np.arange(len(t_names), dtype=np.int64),
                deltas=deltas,
                ctrl_cells=rng.normal(0.0, 0.25, size=(40, len(g_names))),
                ctrl_mean_libnorm=mu * scale,
                ctrl_mean_libnorm_full=mu * scale,
                repr_mask=np.ones(len(t_names), dtype=bool),
                half_a_deltas=ha,
                half_b_deltas=hb,
                half_valid=np.ones(len(t_names), dtype=bool),
                n_cells=np.full(len(t_names), 50.0),
                kd_depth_log1p=np.full(len(t_names), -1.0),
                metadata={"synthetic": True},
            )

        tgt = _fake_line("cal33", 1.0)
        srcs = [_fake_line("vcc2025_h1", 0.9), _fake_line("nadig_hepg2", 1.3), _fake_line("replogle_rpe1", 0.7)]
        thr = check_r15_line_thresholds(tgt, srcs, n_nt_half_a=20, n_nt_half_b=20, min_targets=5, min_genes=100)
        self.assertTrue(thr["passes_thresholds"])

        ev_main = evaluate_target_line_r12(tgt, srcs, only_reproducible=False, n_boot=30, n_perm=10, seed=DEFAULT_R15_SEED)
        self.assertEqual(ev_main["n_sources"], 3)

        sa_h1 = evaluate_sa_h1_leave_vcc2025_out(tgt, srcs, only_reproducible=False, n_boot=30, n_perm=10, seed=DEFAULT_R15_SEED)
        self.assertEqual(sa_h1["n_sources"], 2)
        self.assertNotIn("vcc2025_h1", sa_h1["source_names"])

        sa_panel_pending = evaluate_sa_panel_overlap(tgt, srcs, shared_129_targets=None)
        self.assertEqual(sa_panel_pending["status"], "PENDING_OPEN_SHARED_129_TARGET_LIST")

        sa_panel_ok = evaluate_sa_panel_overlap(
            tgt,
            srcs,
            shared_129_targets=set(t_names[:4]),
            only_reproducible=False,
            n_boot=30,
            n_perm=10,
            seed=DEFAULT_R15_SEED,
        )
        self.assertEqual(sa_panel_ok["subgroups"]["shared_val_test_129"]["n_eval_targets"], 4)
        self.assertEqual(sa_panel_ok["subgroups"]["panel_unique_remaining"]["n_eval_targets"], 4)

    def test_addendum_a_single_line_file_manifest_and_conflicts(self) -> None:
        genes = self.gene_pool[:225]
        v_ids = [g[0] for g in genes]
        v_syms = [g[1] for g in genes]
        rng = np.random.default_rng(1503)

        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            n_cal = 70
            n_hela = 60
            n_jur = 30

            # One file per line, ZERO cell_line/cell_type/line columns (filename alone MUST NOT suffice!)
            h5_cal = tdir / "cal33.h5ad"
            h5_hela = tdir / "hela_screen.h5ad"
            h5_jur = tdir / "jurkat_ctrl.h5ad"

            x_cal = rng.poisson(6.5, size=(n_cal, len(genes))).astype(np.int32)
            x_hela = rng.poisson(6.5, size=(n_hela, len(genes))).astype(np.int32)
            x_jur = rng.poisson(6.5, size=(n_jur, len(genes))).astype(np.int32)

            _write_synthetic_h5ad(
                h5_cal,
                x_cal,
                [f"CAL_{i:04d}" for i in range(n_cal)],
                v_ids,
                v_syms,
                obs_cols={
                    "target_gene": ["non-targeting"] * 40 + ["HDAC3"] * 15 + ["PTEN"] * 15,
                    "ntc_id": [f"ntc_{(i%8)+1}" for i in range(40)] + ["NA"] * 30,
                },
                categorical_cols=(),
            )
            _write_synthetic_h5ad(
                h5_hela,
                x_hela,
                [f"HELA_{i:04d}" for i in range(n_hela)],
                v_ids,
                v_syms,
                obs_cols={
                    "target_gene": ["non-targeting"] * 36 + ["HDAC3"] * 12 + ["PTEN"] * 12,
                    "context": ["baseline"] * 30 + ["stimulated"] * 30,  # unmapped context (no context_to_line_map)
                },
                categorical_cols=("context",),
            )
            _write_synthetic_h5ad(
                h5_jur,
                x_jur,
                [f"JUR_{i:04d}" for i in range(n_jur)],
                v_ids,
                v_syms,
                obs_cols={"target_gene": ["non-targeting"] * n_jur},
                categorical_cols=(),
            )

            # 1) NO manifest (even with filename `cal33.h5ad`!) -> EXCLUDED_FEASIBILITY
            f_no_man = check_vcc_h5ad_feasibility("cal33", [h5_cal], hgnc=self.hgnc, line_manifest=None)
            self.assertFalse(f_no_man["feasible"])
            self.assertEqual(f_no_man["status"], "EXCLUDED_FEASIBILITY")
            f_hela_no_man = check_vcc_h5ad_feasibility("hela", [h5_hela], hgnc=self.hgnc, line_manifest=None)
            self.assertFalse(f_hela_no_man["feasible"])
            self.assertEqual(f_hela_no_man["status"], "EXCLUDED_FEASIBILITY")

            # 2) Bare string / missing URL / missing verbatim quote -> EXCLUDED_FEASIBILITY
            for bad_manifest in (
                {"cal33.h5ad": "cal33"},
                {"cal33.h5ad": {"line": "cal33", "url": "https://virtualcellchallenge.org/2026"}},
                {"cal33.h5ad": {"line": "cal33", "url": "README.md", "quote": "cal33.h5ad contains Cal-33 cells."}},
                {"cal33.h5ad": {"line": "cal33", "quote": "cal33.h5ad contains Cal-33 cells."}},
            ):
                fb = check_vcc_h5ad_feasibility("cal33", [h5_cal], hgnc=self.hgnc, line_manifest=bad_manifest)
                self.assertFalse(fb["feasible"], msg=f"bad manifest unexpectedly passed: {bad_manifest}")
                self.assertEqual(fb["status"], "EXCLUDED_FEASIBILITY")

            # 3) Official manifest (with URL + exact statement) -> ingests cleanly, purges jurkat pre-QC, records URL & quote!
            official_manifest = {
                "cal33.h5ad": {
                    "line": "Cal-33",
                    "url": "https://virtualcellchallenge.org/2026/release",
                    "quote": "File cal33.h5ad contains CRISPRi Perturb-seq cells from the Cal-33 cell line only.",
                },
                "hela_screen.h5ad": {
                    "line": "hela",
                    "url": "https://virtualcellchallenge.org/2026/release",
                    "quote": "File hela_screen.h5ad contains HeLa cells across baseline/stimulated contexts.",
                },
                "jurkat_ctrl.h5ad": {
                    "line": "jurkat",
                    "url": "https://virtualcellchallenge.org/2026/release",
                    "quote": "File jurkat_ctrl.h5ad contains Jurkat reference cells.",
                },
            }

            part = partition_vcc_files_by_line(
                [h5_cal, h5_hela, h5_jur],
                context_to_line_map=None,
                line_manifest=official_manifest,
            )
            self.assertEqual(part["files_by_line"], {"cal33": [h5_cal], "hela": [h5_hela]})
            self.assertIn("jurkat", part["excluded_lines"])
            self.assertEqual(part["n_jurkat_dropped_pre_qc"], n_jur)
            md_prov = format_r15_manifest_provenance_md(part["all_manifest_provenance"])
            self.assertIn("https://virtualcellchallenge.org/2026/release", md_prov)
            self.assertIn("Cal-33 cell line only", md_prov)

            out_cache = tdir / "cache_manifest"
            adata_cal, meta_cal = standardize_vcc_line_to_audit_line(
                "cal33",
                part["files_by_line"]["cal33"],
                cache_out_root=out_cache,
                aligned_genes=v_syms,
                min_gene_overlap=100,
                n_workers=1,
                phase="synthetic_test",
                line_manifest=official_manifest,
                seed=DEFAULT_R15_SEED,
            )
            self.assertEqual(meta_cal["line_identity_modes"], ["OFFICIAL_LINE_MANIFEST"])
            self.assertEqual(len(meta_cal["line_manifest_provenance"]), 1)
            self.assertEqual(
                meta_cal["line_manifest_provenance"][0]["quote"],
                official_manifest["cal33.h5ad"]["quote"],
            )
            self.assertEqual(sorted(adata_cal.target_names), ["HDAC3", "PTEN"])

            # Cross-check built cache vs pure-numpy oracle using line_manifest (`<= 1e-6`)
            oracle = compute_independent_oracle_arrays(
                [h5_cal],
                "cal33",
                line_manifest=official_manifest,
                seed=DEFAULT_R15_SEED,
            )
            line_dir = out_cache / "cal33"
            with np.load(line_dir / "ctrl_cells.npz", allow_pickle=False) as c_r, np.load(
                line_dir / "pseudobulk_delta.npz", allow_pickle=False
            ) as p_r:
                self.assertLessEqual(
                    float(
                        np.max(
                            np.abs(c_r["ctrl_mean_libnorm"] - oracle["ctrl_mean_libnorm"])
                            / np.maximum(np.abs(oracle["ctrl_mean_libnorm"]), 1e-12)
                        )
                    ),
                    1e-6,
                )
                self.assertLessEqual(float(np.max(np.abs(p_r["delta_log1p"] - oracle["delta_log1p"]))), 1e-6)
                self.assertLessEqual(float(np.max(np.abs(p_r["half_a_delta_log1p"] - oracle["half_a_delta_log1p"]))), 1e-6)
                self.assertLessEqual(float(np.max(np.abs(p_r["half_b_delta_log1p"] - oracle["half_b_delta_log1p"]))), 1e-6)

            # Asking for hela on cal33.h5ad via manifest -> EXCLUDED_FEASIBILITY
            f_wrong_req = check_vcc_h5ad_feasibility("hela", [h5_cal], hgnc=self.hgnc, line_manifest=official_manifest)
            self.assertFalse(f_wrong_req["feasible"])
            self.assertEqual(f_wrong_req["status"], "EXCLUDED_FEASIBILITY")

            # 4) Unmapped context column + valid manifest -> (c) succeeds (hela reads cleanly)
            parts_h, meta_h = build_vcc_line_audit_parts(
                "hela",
                [h5_hela],
                hgnc=self.hgnc,
                context_to_line_map=None,
                line_manifest=official_manifest,
            )
            self.assertEqual(meta_h["counts"]["n_kept_ntc"], 36)
            self.assertEqual(meta_h["counts"]["n_kept_targeting_cells"], 24)

            # 5) Conflict between (a) explicit cell_line column and (c) manifest -> affected lines EXCLUDED_FEASIBILITY!
            h5_conf_a = tdir / "conflict_explicit.h5ad"
            _write_synthetic_h5ad(
                h5_conf_a,
                x_cal[:50],
                [f"CA_{i:04d}" for i in range(50)],
                v_ids,
                v_syms,
                obs_cols={
                    "cell_line": ["Cal-33"] * 45 + ["HE_LA"] * 5,
                    "target_gene": ["non-targeting"] * 30 + ["HDAC3"] * 20,
                },
                categorical_cols=("cell_line",),
            )
            conf_manifest = dict(official_manifest)
            conf_manifest["conflict_explicit.h5ad"] = {
                "line": "cal33",
                "url": "https://virtualcellchallenge.org/2026/release",
                "quote": "conflict_explicit.h5ad contains only Cal-33 cells.",
            }
            fc_cal = check_vcc_h5ad_feasibility("cal33", [h5_conf_a], hgnc=self.hgnc, line_manifest=conf_manifest)
            self.assertFalse(fc_cal["feasible"])
            self.assertEqual(fc_cal["status"], "EXCLUDED_FEASIBILITY")
            self.assertEqual(fc_cal["affected_lines"], ["cal33", "hela"])
            fc_hela = check_vcc_h5ad_feasibility("hela", [h5_conf_a], hgnc=self.hgnc, line_manifest=conf_manifest)
            self.assertFalse(fc_hela["feasible"])
            self.assertEqual(fc_hela["status"], "EXCLUDED_FEASIBILITY")
            self.assertEqual(fc_hela["affected_lines"], ["cal33", "hela"])

            # Release-wide batch partition also propagates conflict to taint `hela` even though h5_hela is clean!
            part_conf = partition_vcc_files_by_line([h5_hela, h5_conf_a], line_manifest=conf_manifest)
            self.assertNotIn("cal33", part_conf["files_by_line"])
            self.assertNotIn("hela", part_conf["files_by_line"])
            self.assertIn("cal33", part_conf["excluded_lines"])
            self.assertIn("hela", part_conf["excluded_lines"])

            # 6) Conflict between (b) mapped context column and (c) manifest -> affected lines EXCLUDED_FEASIBILITY!
            h5_conf_b = tdir / "conflict_context.h5ad"
            _write_synthetic_h5ad(
                h5_conf_b,
                x_cal[:50],
                [f"CB_{i:04d}" for i in range(50)],
                v_ids,
                v_syms,
                obs_cols={
                    "context": ["Ctx_AC16"] * 50,
                    "target_gene": ["non-targeting"] * 35 + ["PTEN"] * 15,
                },
                categorical_cols=("context",),
            )
            conf_b_manifest = {
                "conflict_context.h5ad": {
                    "line": "molm13",
                    "url": "https://virtualcellchallenge.org/2026/release",
                    "quote": "conflict_context.h5ad contains MOLM-13 cells.",
                }
            }
            fcb = check_vcc_h5ad_feasibility(
                "ac16",
                [h5_conf_b],
                hgnc=self.hgnc,
                context_to_line_map={"Ctx_AC16": "AC16"},
                line_manifest=conf_b_manifest,
            )
            self.assertFalse(fcb["feasible"])
            self.assertEqual(fcb["status"], "EXCLUDED_FEASIBILITY")
            self.assertEqual(fcb["affected_lines"], ["ac16", "molm13"])

            # 7) Consistent (a) explicit column + (c) manifest -> (a) takes precedence + MANIFEST_VERIFIED succeeds!
            h5_agree = tdir / "agree_explicit_and_manifest.h5ad"
            _write_synthetic_h5ad(
                h5_agree,
                x_cal[:50],
                [f"AG_{i:04d}" for i in range(50)],
                v_ids,
                v_syms,
                obs_cols={
                    "cell_line": ["Cal-33"] * 50,
                    "target_gene": ["non-targeting"] * 35 + ["HDAC3"] * 15,
                },
                categorical_cols=("cell_line",),
            )
            agree_manifest = {
                "agree_explicit_and_manifest.h5ad": {
                    "line": "cal33",
                    "url": "https://virtualcellchallenge.org/2026/release",
                    "quote": "agree_explicit_and_manifest.h5ad contains Cal-33 cells.",
                }
            }
            f_ag = check_vcc_h5ad_feasibility("cal33", [h5_agree], hgnc=self.hgnc, line_manifest=agree_manifest)
            self.assertTrue(f_ag["feasible"])
            self.assertEqual(f_ag["files"][0]["line_identity_mode"], "EXPLICIT_LINE_COL:cell_line+MANIFEST_VERIFIED")

            # 8) Passing ALL release directory files (`[h5_cal, h5_hela, h5_jur]`) to `standardize_vcc_line_to_audit_line`:
            #    other-line files contribute 0 cells (logged cleanly, NO error) and output arrays match single-file load 0.0!
            all_dir_files = [h5_cal, h5_hela, h5_jur]
            out_dir_cache = tdir / "cache_all_dir"
            adata_cal_dir, meta_cal_dir = standardize_vcc_line_to_audit_line(
                "cal33",
                all_dir_files,
                cache_out_root=out_dir_cache,
                aligned_genes=v_syms,
                min_gene_overlap=100,
                n_workers=1,
                phase="synthetic_test",
                line_manifest=official_manifest,
                seed=DEFAULT_R15_SEED,
            )
            self.assertEqual(meta_cal_dir["counts"]["n_files_contributing"], 1)
            self.assertEqual(meta_cal_dir["counts"]["n_files_zero_cells_other_lines"], 2)
            self.assertEqual(meta_cal_dir["counts"]["n_jurkat_dropped_in_other_files"], n_jur)
            self.assertEqual(len(meta_cal_dir["other_line_files"]), 2)
            self.assertTrue(all(r["n_line_raw_cells"] == 0 and not r["contributes_cells"] for r in meta_cal_dir["other_line_files"]))
            self.assertEqual(adata_cal_dir.target_names, adata_cal.target_names)
            cmp_dir_vs_single = compare_line_dirs_exact(out_cache / "cal33", out_dir_cache / "cal33", tol=0.0)
            self.assertTrue(cmp_dir_vs_single["pass"])

            # Same release file list also ingests `hela` cleanly!
            adata_hela_dir, meta_hela_dir = standardize_vcc_line_to_audit_line(
                "hela",
                all_dir_files,
                cache_out_root=out_dir_cache,
                aligned_genes=v_syms,
                min_gene_overlap=100,
                n_workers=1,
                phase="synthetic_test",
                line_manifest=official_manifest,
                seed=DEFAULT_R15_SEED,
            )
            self.assertEqual(meta_hela_dir["counts"]["n_files_contributing"], 1)
            self.assertEqual(meta_hela_dir["counts"]["n_files_zero_cells_other_lines"], 2)
            self.assertEqual(sorted(adata_hela_dir.target_names), ["HDAC3", "PTEN"])

            # A line absent from ALL files in the release directory (`ac16`: 0 cells in all 3 files) -> EXCLUDED_FEASIBILITY!
            f_ac16_zero = check_vcc_h5ad_feasibility("ac16", all_dir_files, hgnc=self.hgnc, line_manifest=official_manifest)
            self.assertFalse(f_ac16_zero["feasible"])
            self.assertEqual(f_ac16_zero["status"], "EXCLUDED_FEASIBILITY")
            self.assertIn("contributed 0 cells across all 3 input files", f_ac16_zero["reason"])
            with self.assertRaises(RuntimeError):
                standardize_vcc_line_to_audit_line(
                    "ac16",
                    all_dir_files,
                    cache_out_root=out_dir_cache,
                    aligned_genes=v_syms,
                    min_gene_overlap=100,
                    n_workers=1,
                    phase="synthetic_test",
                    line_manifest=official_manifest,
                    seed=DEFAULT_R15_SEED,
                )

            # Even if a separate line's file (`hela`) has invalid normalized floats, passing the whole dir to `cal33` succeeds
            # while `hela` itself is blocked!
            h5_hela_float = tdir / "hela_float_bad.h5ad"
            _write_synthetic_h5ad(
                h5_hela_float,
                np.log1p(x_hela.astype(np.float32)),
                [f"HF_{i:04d}" for i in range(n_hela)],
                v_ids,
                v_syms,
                obs_cols={"target_gene": ["non-targeting"] * 36 + ["HDAC3"] * 24},
                categorical_cols=(),
            )
            float_manifest = {
                "cal33.h5ad": official_manifest["cal33.h5ad"],
                "hela_float_bad.h5ad": {
                    "line": "hela",
                    "url": "https://virtualcellchallenge.org/2026/release",
                    "quote": "hela_float_bad.h5ad contains HeLa cells.",
                },
            }
            f_cal_ok = check_vcc_h5ad_feasibility("cal33", [h5_cal, h5_hela_float], hgnc=self.hgnc, line_manifest=float_manifest)
            self.assertTrue(f_cal_ok["feasible"])
            f_hela_bad = check_vcc_h5ad_feasibility("hela", [h5_cal, h5_hela_float], hgnc=self.hgnc, line_manifest=float_manifest)
            self.assertFalse(f_hela_bad["feasible"])
            self.assertEqual(f_hela_bad["status"], "EXCLUDED_FEASIBILITY")


if __name__ == "__main__":
    unittest.main(verbosity=2)
