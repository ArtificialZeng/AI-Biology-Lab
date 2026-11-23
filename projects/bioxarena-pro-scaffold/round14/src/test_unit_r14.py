#!/usr/bin/env python3
"""
Unit & Integration Tests for Round 14 (`PREREG_G`, `GSE337804`) Modules (`round14/src/`).

Runs strictly on synthetic `.h5ad` / `sgRNA_lib.csv.gz` fixtures without touching any real
Round 14 expression data or modifying any frozen file under `round10/`–`round13/`.
"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from typing import Dict, Optional, Sequence

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

R14_SRC_DIR = Path(__file__).resolve().parent
REPO_ROOT = R14_SRC_DIR.parents[1]
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
R10_SCRIPTS_DIR = REPO_ROOT / "round10" / "shared" / "scripts"
R11_SRC_DIR = REPO_ROOT / "round11" / "src"
R12_SRC_DIR = REPO_ROOT / "round12" / "src"
for p in (str(R10_SRC_DIR), str(R10_SCRIPTS_DIR), str(R11_SRC_DIR), str(R12_SRC_DIR), str(R14_SRC_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from data_loader_r10 import LineDatasetR10, load_standardized_line_r10  # noqa: E402
from prepare_r10_lines import HGNCLookup, build_hgnc_lookup  # noqa: E402
from test_unit_r11 import _make_planted_synthetic_cohort  # noqa: E402
from run_eval_r12 import evaluate_target_against_sources_r12  # noqa: E402

from adjudicate_r14 import adjudicate_r14_results, render_r14_markdown_report  # noqa: E402
from data_loader_r14 import (  # noqa: E402
    DEFAULT_R14_SEED,
    HGNCBridge,
    check_prereg_g_gate,
    load_r12_a549_tested_targets,
)
from dataset_readers_r14 import (  # noqa: E402
    build_gse337804_parts,
    check_gse337804_feasibility,
    prepare_r14_line,
)
from run_eval_r14 import compute_overlap_subgroup_descriptive_r14  # noqa: E402

_CACHED_HGNC: Optional[HGNCLookup] = None


def _get_hgnc() -> HGNCLookup:
    global _CACHED_HGNC
    if _CACHED_HGNC is None:
        _CACHED_HGNC = build_hgnc_lookup()
    return _CACHED_HGNC


def _pick_hgnc_genes(n_genes: int) -> list:
    hgnc = _get_hgnc()
    picked = [g for g in sorted(hgnc.approved) if g and not g.startswith("MT-") and g not in {"CAS9"}][:n_genes]
    if len(picked) < n_genes:
        raise RuntimeError(f"Requested {n_genes} HGNC genes, got {len(picked)}")
    return picked


def _write_sparse_group(grp: h5py.Group, mat: sp.spmatrix, layout: str) -> None:
    if layout == "csc":
        m = mat.tocsc()
        enc = "csc_matrix"
    else:
        m = mat.tocsr()
        enc = "csr_matrix"
    grp.attrs["encoding-type"] = enc
    grp.attrs["shape"] = np.array(m.shape, dtype=np.int64)
    dt = np.int32 if np.issubdtype(mat.dtype, np.integer) else np.float32
    grp.create_dataset("data", data=m.data.astype(dt))
    grp.create_dataset("indices", data=m.indices.astype(np.int32))
    grp.create_dataset("indptr", data=m.indptr.astype(np.int64))


def _write_str_dataset(grp: h5py.Group, name: str, vals: Sequence[str], categorical: bool = False) -> None:
    if not categorical:
        grp.create_dataset(name, data=np.array([str(x).encode("utf-8") for x in vals]))
        return
    uniq: list = []
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
    var_ids: Optional[Sequence[str]] = None,
    feature_types: Optional[Sequence[str]] = None,
    x_layout: str = "csr",
    layers: Optional[Dict[str, tuple]] = None,  # name -> (mat, layout)
    raw_tuple: Optional[tuple] = None,  # (mat, layout, raw_var_names, raw_var_ids, raw_ftypes)
    categorical_obs: bool = False,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    n_obs = X_mat.shape[0]
    n_var = len(var_names)
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
        bcs = [f"{path.stem}_C{i:04d}" for i in range(n_obs)]
        _write_str_dataset(og, "_index", bcs, categorical=False)
        if obs_df is not None:
            for col in obs_df.columns:
                _write_str_dataset(og, col, obs_df[col].tolist(), categorical=categorical_obs)

        vg = f.create_group("var")
        vg.attrs["_index"] = "_index"
        _write_str_dataset(vg, "_index", var_names, categorical=False)
        if var_ids is not None:
            _write_str_dataset(vg, "ensembl_id", var_ids, categorical=False)
        if feature_types is not None:
            _write_str_dataset(vg, "feature_types", feature_types, categorical=False)
        assert n_var == len(var_names)

        if raw_tuple is not None:
            rmat, rlay, r_names, r_ids, r_ftypes = raw_tuple
            rg = f.create_group("raw")
            if rlay == "dense":
                rg.create_dataset("X", data=(rmat.toarray() if sp.issparse(rmat) else np.asarray(rmat)))
            else:
                rxg = rg.create_group("X")
                _write_sparse_group(rxg, sp.csc_matrix(rmat) if rlay == "csc" else sp.csr_matrix(rmat), rlay)
            rvg = rg.create_group("var")
            rvg.attrs["_index"] = "_index"
            _write_str_dataset(rvg, "_index", r_names, categorical=False)
            if r_ids is not None:
                _write_str_dataset(rvg, "ensembl_id", r_ids, categorical=False)
            if r_ftypes is not None:
                _write_str_dataset(rvg, "feature_types", r_ftypes, categorical=False)


class Round14UnitTests(unittest.TestCase):
    def test_1_blinding_gate_adjudicator_and_r12_overlap_loader(self) -> None:
        self.assertEqual(DEFAULT_R14_SEED, 20261224)
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)
            with self.assertRaises(RuntimeError):
                check_prereg_g_gate(p)
            (p / "PREREG_G_COMMIT.txt").write_text("prereg_g_111\n", encoding="utf-8")
            with self.assertRaises(RuntimeError):
                check_prereg_g_gate(p)
            (p / "CODE_FREEZE_R14.txt").write_text("freeze_g_222\n", encoding="utf-8")
            pr, fr = check_prereg_g_gate(p)
            self.assertEqual(pr, "prereg_g_111")
            self.assertEqual(fr, "freeze_g_222")

        hgnc_b = HGNCBridge.from_tsv()
        r12_tgts = load_r12_a549_tested_targets(hgnc_bridge=hgnc_b)
        self.assertGreater(len(r12_tgts), 400)

        src6, tgt_syn = _make_planted_synthetic_cohort(seed=141, n_genes=150, n_targets=55, n_sources=6)
        ds_myllia = LineDatasetR10(**{**tgt_syn.__dict__, "name": "a549_myllia"})
        l1, b1, ex1 = evaluate_target_against_sources_r12(
            ds_myllia, src6, include_r9_anchors=False, min_gene_overlap=100, n_boot=120, seed=DEFAULT_R14_SEED, n_perm=10
        )
        df_l = pd.DataFrame(l1)
        df_b = pd.DataFrame(b1)
        summary = {
            "phase": "holdout_r14",
            "source_mode": "r9_6",
            "source_lines": sorted(src6.keys()),
            "prereg_commit_hash": "prereg_g",
            "code_freeze_hash": "freeze_r14",
            "bootstrap_seed": DEFAULT_R14_SEED,
            "n_boot": 120,
            "overlap_subgroups": {},
            "fct_perm100": {"a549_myllia": ex1["perm100"]},
            "n5_cone_metrics": {"a549_myllia": ex1["n5_cone_metrics"]},
            "n5_cone_split_metrics": {"a549_myllia": ex1["n5_cone_split_metrics"]},
            "n5_perm100": {"a549_myllia": ex1["n5_perm100"]},
        }

        with self.assertRaises(ValueError):
            adjudicate_r14_results(df_l, df_b, summary, confirmatory_lines=["a549_ko"])
        adj0 = adjudicate_r14_results(df_l, df_b, summary, confirmatory_lines=[])
        self.assertEqual(adj0["confirmatory_status"], "NO_CONFIRMATORY_LINES_DESCRIPTIVE_ONLY")
        adj1 = adjudicate_r14_results(df_l, df_b, summary, confirmatory_lines=["a549_myllia"])
        self.assertEqual(adj1["confirmatory_status"], "CONFIRMATORY_SINGLE_LINE_K1")
        self.assertTrue(adj1["same_cell_line_independent_dataset_only_notice"])
        self.assertTrue(adj1["no_cross_round_pooling_notice"])
        md = render_r14_markdown_report(df_l, df_b, summary, adj1, "unit_r14")
        self.assertIn("PREREG_G", md)
        self.assertIn("20261224", md)
        self.assertIn("Same-Cell-Line", md)

    def test_2_h5ad_integer_slot_autodetect_layers_raw_and_excluded(self) -> None:
        rng = np.random.default_rng(1402)
        hgnc = _get_hgnc()
        genes = _pick_hgnc_genes(240)
        g0, g1 = genes[0], genes[1]
        raw_genes = genes + ["A2M"]  # Extra gene only in raw/var

        with tempfile.TemporaryDirectory() as td:
            raw_dir = Path(td) / "raw"
            gdir = raw_dir / "gse337804"
            gdir.mkdir(parents=True)

            # Build float-normalized X AND integer counts in layers/counts (ch_1) vs raw/X + raw/var (ch_2)
            for ch_n, use_raw in ((1, False), (2, True)):
                rows_meta = (
                    [("non-targeting_1", "non-targeting", -1)] * 16
                    + [(f"{g0}_1", g0, 0)] * 8
                    + [(f"{g1}_1", g1, 1)] * 8
                )
                n_obs = len(rows_meta)
                X_int = rng.integers(2, 5, size=(n_obs, len(genes)), dtype=np.int32)
                for r, (_, _, tidx) in enumerate(rows_meta):
                    if tidx >= 0:
                        X_int[r, :10] += 5
                        X_int[r, tidx] = 1
                X_float = np.log1p(X_int.astype(np.float32))
                obs_df = pd.DataFrame(
                    {
                        "sgrna_id": [x[0] for x in rows_meta],
                        "sgrna_symbol": [x[1] for x in rows_meta],
                    }
                )
                hp = gdir / f"GSE337804_A549_ch_{ch_n}.RDS.h5ad"
                if not use_raw:
                    _write_synthetic_h5ad(
                        hp,
                        X_float,
                        obs_df,
                        genes,
                        x_layout="csr",
                        layers={"counts": (X_int, "csr")},
                    )
                else:
                    X_raw_extra = np.hstack([X_int, rng.integers(1, 3, size=(n_obs, 1), dtype=np.int32)])
                    _write_synthetic_h5ad(
                        hp,
                        X_float,
                        obs_df,
                        genes,
                        x_layout="csr",
                        layers={"scaled": (X_float * 0.5, "dense")},
                        raw_tuple=(X_raw_extra, "csc", raw_genes, None, None),
                    )

            feas = check_gse337804_feasibility(raw_dir, hgnc=hgnc)
            self.assertTrue(feas["feasible"], feas.get("reason"))
            self.assertEqual(feas["channel_slots"]["ch_1"]["x_slot_used"], "layers/counts")
            self.assertEqual(feas["channel_slots"]["ch_2"]["x_slot_used"], "raw/X")

            parts, audit = build_gse337804_parts(raw_dir, hgnc=hgnc)
            self.assertEqual(audit["channel_slots"]["ch_1"]["x_slot_used"], "layers/counts")
            self.assertEqual(audit["channel_slots"]["ch_2"]["x_slot_used"], "raw/X")
            self.assertEqual(audit["n_qc_kept_cells"], 64)
            self.assertEqual(parts[0].n_genes_raw, len(genes))
            self.assertEqual(parts[1].n_genes_raw, len(raw_genes))

            # If NO slot contains raw integers -> EXCLUDED_FEASIBILITY
            bad_dir = Path(td) / "raw_bad" / "gse337804"
            bad_dir.mkdir(parents=True)
            for ch_n in (1, 2):
                obs_b = pd.DataFrame({"sgrna_id": ["non-targeting_1"] * 4 + [f"{g0}_1"] * 4, "sgrna_symbol": ["non-targeting"] * 4 + [g0] * 4})
                Xf = rng.uniform(0.1, 3.5, size=(8, len(genes))).astype(np.float32)
                _write_synthetic_h5ad(bad_dir / f"GSE337804_A549_ch_{ch_n}.RDS.h5ad", Xf, obs_b, genes, x_layout="csr")
            feas_bad = check_gse337804_feasibility(Path(td) / "raw_bad", hgnc=hgnc)
            self.assertFalse(feas_bad["feasible"])
            with self.assertRaisesRegex(RuntimeError, "EXCLUDED_FEASIBILITY"):
                build_gse337804_parts(Path(td) / "raw_bad", hgnc=hgnc)

    def test_3_gse337804_ctrlNT_vs_ctrlT_depletex_parity_and_process_r10(self) -> None:
        rng = np.random.default_rng(1403)
        hgnc = _get_hgnc()
        genes_base = _pick_hgnc_genes(260)
        # Include vector CAS9 feature (must be dropped) and MT-CO1 with ~45% MT reads (MUST be retained since max_pct_mt=None!)
        var_names = genes_base + ["MT-CO1", "CAS9"]
        g0, g1 = genes_base[0], genes_base[1]

        with tempfile.TemporaryDirectory() as td:
            raw_dir = Path(td) / "raw"
            gdir = raw_dir / "gse337804"
            gdir.mkdir(parents=True)

            # Write authentic GSE337804_sgRNA_lib.csv.gz schema distinguishing 20 ctrl_NT vs 20 ctrl_T
            lib_rows = []
            for k in range(1, 21):
                lib_rows.append({"id": f"non-targeting_{k}", "symbol": "non-targeting", "rank": k, "class": "ctrl_NT", "entrez": "", "guide": "AAAA", "count": 2000})
                lib_rows.append({"id": f"gene_desert_{k}", "symbol": "gene_desert", "rank": k, "class": "ctrl_T", "entrez": "", "guide": "CCCC", "count": 3000})
            for gi, gsym in enumerate((g0, g1)):
                for rk in range(1, 4):
                    lib_rows.append({"id": f"{gsym}_{rk}", "symbol": gsym, "rank": rk, "class": "target", "entrez": str(100 + gi), "guide": "GGGG", "count": 4500})
            pd.DataFrame(lib_rows).to_csv(gdir / "GSE337804_sgRNA_lib.csv.gz", index=False, compression="gzip")

            # ch_1: 20 NTC, 4 gene_desert safe-targeting, 10 g0, 10 g1, 2 multi-guide '|', 2 low-UMI (categorical obs)
            # ch_2: 12 NTC (< 15 -> falls back to half_2_pooled!), 4 gene_desert, 10 g0, 10 g1, 1 unassigned, 1 low-ngenes
            configs = [
                (1, "csr", True, 20),
                (2, "csc", False, 12),
            ]
            for ch_n, lay, cat_obs, n_ntc in configs:
                cells = []
                for k in range(n_ntc):
                    cells.append((f"non-targeting_{(k % 20) + 1}", "non-targeting", -1, False, True))
                for k in range(4):
                    cells.append((f"gene_desert_{k + 1}", "gene_desert", -1, False, False))
                for k in range(10):
                    cells.append((f"{g0}_{(k % 3) + 1}", g0, 0, False, False))
                for k in range(10):
                    cells.append((f"{g1}_{(k % 3) + 1}", g1, 1, False, False))
                cells.append((f"{g0}_1|{g1}_1", f"{g0}|{g1}", -1, False, False))  # multi-guide excluded
                cells.append(("None", "None", -1, False, False))  # unassigned excluded
                cells.append((f"{g0}_1", g0, 0, True, False))  # low RNA QC excluded

                mat = np.zeros((len(cells), len(var_names)), dtype=np.int32)
                for r, (_, _, tidx, is_low, is_high_mt) in enumerate(cells):
                    if is_low:
                        mat[r, :40] = 1  # < 300 UMI & < 200 genes
                    else:
                        mat[r, : len(genes_base)] = rng.integers(2, 5, size=len(genes_base), dtype=np.int32)
                        if is_high_mt:
                            mat[r, len(genes_base)] = 600  # ~45% MT-CO1 (must NOT be dropped!)
                        if tidx >= 0:
                            mat[r, :15] += 6
                            mat[r, tidx] = 1
                    mat[r, -1] = 10  # CAS9 vector gene (must be excluded from RNA features)

                obs_df = pd.DataFrame(
                    {
                        "sgrna_id": [c[0] for c in cells],
                        "sgrna_symbol": [c[1] for c in cells],
                    }
                )
                _write_synthetic_h5ad(
                    gdir / f"GSE337804_A549_ch_{ch_n}.RDS.h5ad",
                    mat,
                    obs_df,
                    var_names,
                    x_layout=lay,
                    categorical_obs=cat_obs,
                )

            feas = check_gse337804_feasibility(raw_dir, hgnc=hgnc)
            self.assertTrue(feas["feasible"], feas.get("reason"))
            self.assertEqual(feas["channels"], ["ch_1", "ch_2"])
            self.assertEqual(feas["n_obs_safe_meta"], 8)

            parts, audit = build_gse337804_parts(raw_dir, hgnc=hgnc)
            self.assertIsNone(audit["max_pct_mt"])
            self.assertEqual(audit["half_rule_used"], "ch_1_HalfA_vs_ch_2_HalfB")
            self.assertEqual(audit["n_safe_targeting_excluded"], 8)
            self.assertEqual(audit["n_multi_or_zero_guide_excluded"], 4)
            self.assertEqual(audit["n_rna_qc_excluded"], 2)
            self.assertEqual(audit["n_kept_half_a_ntc"], 20)
            self.assertEqual(audit["n_kept_half_b_ntc"], 12)
            self.assertEqual(set(audit["channels_used"]), {"ch_1", "half_2_pooled"})
            self.assertEqual(parts[0].half_flags[0], 1)
            self.assertEqual(parts[1].half_flags[0], 2)
            self.assertNotIn("CAS9", parts[0].var_gene_names)

            out_root = Path(td) / "lines"
            prepare_r14_line("a549_myllia", raw_dir, out_root, hgnc=hgnc, workers=1)
            ds = load_standardized_line_r10(out_root / "a549_myllia", min_gene_overlap=50)
            self.assertEqual(len(ds.target_names), 2)
            self.assertEqual(len(ds.genes), len(genes_base) + 1)  # base + MT-CO1 (CAS9 dropped)

    def test_4_gse337804_crispr_em_fallback_on_ngenes200(self) -> None:
        rng = np.random.default_rng(1404)
        hgnc = _get_hgnc()
        genes = _pick_hgnc_genes(240)
        g0 = genes[0]
        guides = ["non-targeting_1", "gene_desert_1", f"{g0}_1"]
        var_names = genes + guides
        ftypes = ["Gene Expression"] * len(genes) + ["CRISPR Guide Capture"] * len(guides)

        with tempfile.TemporaryDirectory() as td:
            raw_dir = Path(td) / "raw"
            gdir = raw_dir / "gse337804"
            gdir.mkdir(parents=True)
            pd.DataFrame(
                [
                    {"id": "non-targeting_1", "symbol": "non-targeting", "class": "ctrl_NT"},
                    {"id": "gene_desert_1", "symbol": "gene_desert", "class": "ctrl_T"},
                    {"id": f"{g0}_1", "symbol": g0, "class": "target"},
                ]
            ).to_csv(gdir / "GSE337804_sgRNA_lib.csv.gz", index=False, compression="gzip")

            # Omit .obs guide columns -> trigger PREREG_F §3.2 bimodal EM fitted on ngenes >= 200
            for ch_n in (1, 2):
                n_qc = 16  # 10 NTC, 2 safe, 4 g0 (< 50 entries in ok_qc)
                n_fit_only = 55  # ngenes=210 >= 200, umi=210 < 300 -> excluded from ok_qc but feeds EM fit!
                n_tot = n_qc + n_fit_only
                mat = np.zeros((n_tot, len(var_names)), dtype=np.int32)
                for r in range(n_qc):
                    mat[r, : len(genes)] = rng.integers(2, 5, size=len(genes), dtype=np.int32)
                    gr = len(genes) + (0 if r < 10 else (1 if r < 12 else 2))
                    mat[r, gr] = 24
                for j in range(n_fit_only):
                    r = n_qc + j
                    mat[r, :210] = 1  # ngenes = 210 >= 200, UMI = 210 < 300
                    for g_col in range(len(genes), len(var_names)):
                        mat[r, g_col] = 1 if j < 38 else 28
                _write_synthetic_h5ad(
                    gdir / f"GSE337804_A549_ch_{ch_n}.RDS.h5ad",
                    mat,
                    obs_df=None,
                    var_names=var_names,
                    feature_types=ftypes,
                    x_layout="csr",
                )

            feas = check_gse337804_feasibility(raw_dir, hgnc=hgnc)
            self.assertTrue(feas["feasible"], feas.get("reason"))
            self.assertEqual(feas["guide_call_mode"], "h5ad_crispr_em_fallback")

            parts, audit = build_gse337804_parts(raw_dir, hgnc=hgnc)
            self.assertEqual(audit["guide_call_mode"], "h5ad_crispr_em_fallback")
            self.assertNotEqual(audit["median_guide_threshold"], 7.0)  # EM fitted on ngenes>=200, not default 7.0!
            self.assertEqual(audit["n_safe_targeting_excluded"], 4)
            self.assertEqual(audit["n_kept_ntc"], 20)
            self.assertEqual(audit["n_qc_kept_cells"], 28)

    def test_5_overlap_subgroup_descriptive_reporting_r14(self) -> None:
        src6, tgt_syn = _make_planted_synthetic_cohort(seed=145, n_genes=160, n_targets=56, n_sources=6)
        ds_myllia = LineDatasetR10(**{**tgt_syn.__dict__, "name": "a549_myllia"})
        overlap_set = set(tgt_syn.target_names[:20])

        sg_dict, flat_rows = compute_overlap_subgroup_descriptive_r14(
            target_ds=ds_myllia,
            source_datasets=src6,
            r12_overlap_targets=overlap_set,
            include_r9_anchors=False,
            n_boot=120,
            seed=DEFAULT_R14_SEED,
        )
        self.assertEqual(sg_dict["in_r12_a549_overlap"]["n_targets"], 20)
        self.assertEqual(sg_dict["novel_not_in_r12_a549"]["n_targets"], 36)
        self.assertEqual(len(sg_dict["in_r12_a549_overlap"]["paired_fct_minus_mismatched"]), 6)
        self.assertEqual(len(sg_dict["novel_not_in_r12_a549"]["paired_fct_minus_mismatched"]), 6)
        self.assertEqual(len(flat_rows), 2 * (1 + 6))
        self.assertGreater(sg_dict["in_r12_a549_overlap"]["fct"]["estimate"], 0.1)
        self.assertGreater(sg_dict["novel_not_in_r12_a549"]["fct"]["estimate"], 0.1)

        l1, b1, ex1 = evaluate_target_against_sources_r12(
            ds_myllia, src6, include_r9_anchors=False, min_gene_overlap=100, n_boot=120, seed=DEFAULT_R14_SEED, n_perm=10
        )
        df_l, df_b = pd.DataFrame(l1), pd.DataFrame(b1)
        summary = {
            "phase": "holdout_r14",
            "source_mode": "r9_6",
            "source_lines": sorted(src6.keys()),
            "prereg_commit_hash": "prereg_g",
            "code_freeze_hash": "freeze_r14",
            "bootstrap_seed": DEFAULT_R14_SEED,
            "n_boot": 120,
            "overlap_subgroups": {"a549_myllia": sg_dict},
            "fct_perm100": {"a549_myllia": ex1["perm100"]},
            "n5_cone_metrics": {"a549_myllia": ex1["n5_cone_metrics"]},
            "n5_cone_split_metrics": {"a549_myllia": ex1["n5_cone_split_metrics"]},
            "n5_perm100": {"a549_myllia": ex1["n5_perm100"]},
        }
        adj = adjudicate_r14_results(df_l, df_b, summary, confirmatory_lines=["a549_myllia"])
        md = render_r14_markdown_report(df_l, df_b, summary, adj, "unit_subgroups_r14")
        self.assertIn("PREREG_G` Item 7", md)
        self.assertIn("In Round 12 `a549_ko` Overlap", md)
        self.assertIn("Novel (Not Tested in Round 12 `a549_ko`)", md)

    def test_6_slot_ambiguous_excluded_layers_full_int_check_and_explicit_x_layer(self) -> None:
        rng = np.random.default_rng(1406)
        hgnc = _get_hgnc()
        genes = _pick_hgnc_genes(240)
        g0 = genes[0]

        with tempfile.TemporaryDirectory() as td:
            # 1) X is integer (SCT-like different counts) AND raw/X is integer -> SLOT_AMBIGUOUS
            dir_amb = Path(td) / "amb" / "gse337804"
            dir_amb.mkdir(parents=True)
            for ch_n in (1, 2):
                n_obs = 40
                X_raw = rng.integers(2, 5, size=(n_obs, len(genes)), dtype=np.int32)
                X_sct = np.maximum(X_raw - 1, 0).astype(np.int32)
                obs_df = pd.DataFrame(
                    {
                        "sgrna_id": ["non-targeting_1"] * 20 + [f"{g0}_1"] * 20,
                        "sgrna_symbol": ["non-targeting"] * 20 + [g0] * 20,
                    }
                )
                _write_synthetic_h5ad(
                    dir_amb / f"GSE337804_A549_ch_{ch_n}.RDS.h5ad",
                    X_sct,
                    obs_df,
                    genes,
                    x_layout="csr",
                    raw_tuple=(X_raw, "csr", genes, None, None),
                )

            feas_amb = check_gse337804_feasibility(Path(td) / "amb", hgnc=hgnc, x_layer="auto")
            self.assertFalse(feas_amb["feasible"])
            self.assertEqual(feas_amb["status"], "SLOT_AMBIGUOUS")
            self.assertIn("median_total_umi", feas_amb["reason"])
            with self.assertRaisesRegex(RuntimeError, "SLOT_AMBIGUOUS"):
                build_gse337804_parts(Path(td) / "amb", hgnc=hgnc, x_layer="auto")

            # Controller explicit override x_layer="raw/X" resolves cleanly & logs candidate medians
            parts_exp, audit_exp = build_gse337804_parts(Path(td) / "amb", hgnc=hgnc, x_layer="raw/X")
            self.assertEqual(audit_exp["channel_slots"]["ch_1"]["x_slot_used"], "raw/X")
            self.assertGreater(audit_exp["channel_slots"]["ch_1"]["selected_median_total_umi"], 600.0)
            self.assertEqual(len(parts_exp), 2)

            # 2) layers/spliced (and unspliced/nascent/velocity/lognorm/scale) with integers MUST be excluded,
            #    while element-wise identical X & gene-permuted raw/X resolves cleanly to X!
            #    Also inject one decimal value 0.5 at row index 29 (> 24 head/mid sample) inside X of another test.
            dir_spliced = Path(td) / "spliced" / "gse337804"
            dir_spliced.mkdir(parents=True)
            perm_idx = np.arange(len(genes))[::-1]
            genes_rev = list(reversed(genes))
            for ch_n in (1, 2):
                n_obs = 40
                X_raw = rng.integers(2, 5, size=(n_obs, len(genes)), dtype=np.int32)
                X_spliced = rng.integers(1, 2, size=(n_obs, len(genes)), dtype=np.int32)
                obs_df = pd.DataFrame(
                    {
                        "sgrna_id": ["non-targeting_1"] * 20 + [f"{g0}_1"] * 20,
                        "sgrna_symbol": ["non-targeting"] * 20 + [g0] * 20,
                    }
                )
                _write_synthetic_h5ad(
                    dir_spliced / f"GSE337804_A549_ch_{ch_n}.RDS.h5ad",
                    X_raw,
                    obs_df,
                    genes,
                    x_layout="csr",
                    layers={
                        "spliced": (X_spliced, "csr"),
                        "unspliced_counts": (X_spliced, "csr"),
                        "nascent": (X_spliced, "csr"),
                        "velocity": (X_spliced, "csr"),
                        "lognorm": (X_spliced, "csr"),
                        "scaled": (X_spliced, "csr"),
                    },
                    raw_tuple=(X_raw[:, perm_idx], "csr", genes_rev, None, None),
                )

            feas_sp = check_gse337804_feasibility(Path(td) / "spliced", hgnc=hgnc, x_layer="auto")
            self.assertTrue(feas_sp["feasible"], feas_sp.get("reason"))
            ch1_diag = feas_sp["channel_slots"]["ch_1"]["slot_diagnostics"]
            self.assertEqual(feas_sp["channel_slots"]["ch_1"]["x_slot_used"], "X")
            self.assertEqual(ch1_diag["slots"]["layers/spliced"]["status"], "excluded_subfraction_or_derived_layer")
            self.assertEqual(ch1_diag["slots"]["layers/scaled"]["status"], "excluded_subfraction_or_derived_layer")
            self.assertEqual(ch1_diag["selection_mode"], "identical_multi_slot_prefer_X")
            self.assertEqual(ch1_diag["comparisons"]["X_vs_raw/X"]["identical"], True)

            # 3) Single decimal non-zero entry deep inside CSR data MUST fail full non-zero validation!
            dir_tail = Path(td) / "tail_float" / "gse337804"
            dir_tail.mkdir(parents=True)
            for ch_n in (1, 2):
                n_obs = 60
                X_int = rng.integers(2, 4, size=(n_obs, len(genes)), dtype=np.int32)
                X_flt = X_int.astype(np.float32)
                X_flt[35, 100] = 2.5  # single decimal buried mid-matrix
                obs_df = pd.DataFrame(
                    {
                        "sgrna_id": ["non-targeting_1"] * 30 + [f"{g0}_1"] * 30,
                        "sgrna_symbol": ["non-targeting"] * 30 + [g0] * 30,
                    }
                )
                _write_synthetic_h5ad(
                    dir_tail / f"GSE337804_A549_ch_{ch_n}.RDS.h5ad",
                    X_flt,
                    obs_df,
                    genes,
                    x_layout="csr",
                )
            feas_tail = check_gse337804_feasibility(Path(td) / "tail_float", hgnc=hgnc)
            self.assertFalse(feas_tail["feasible"])
            with self.assertRaisesRegex(RuntimeError, "EXCLUDED_FEASIBILITY"):
                build_gse337804_parts(Path(td) / "tail_float", hgnc=hgnc)


if __name__ == "__main__":
    unittest.main(verbosity=2)
