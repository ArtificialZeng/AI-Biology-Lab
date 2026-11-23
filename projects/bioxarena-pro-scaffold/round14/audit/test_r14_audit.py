#!/usr/bin/env python3
"""
Synthetic Fixture Unit & End-to-End Tests for Round 14 Clean-Room Audit (`PREREG_G`).

Strictly uses ONLY temporary synthetic `.h5ad` and `.csv.gz` fixtures.
Never touches `round14/shared/`, `round14/src/`, or real `GSE337804` data.
"""

from __future__ import annotations

import csv
import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import h5py
import numpy as np
import scipy.sparse as sp

AUDIT_R14_DIR = Path(__file__).resolve().parent
REPO_ROOT = AUDIT_R14_DIR.parents[1]
for _p in (
    str(AUDIT_R14_DIR),
    str(REPO_ROOT / "round10" / "shared" / "scripts"),
    str(REPO_ROOT / "round11" / "audit"),
    str(REPO_ROOT / "round12" / "audit"),
):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from independent_r11_loader import AuditLineData, DEFAULT_HGNC_TSV  # noqa: E402
from independent_r14_loader import (  # noqa: E402
    R14_SHARED_DIR,
    R14_TARGET_KEY,
    build_a549_myllia_audit_parts,
    check_gse337804_h5ad_feasibility,
    enforce_r14_gate,
    load_gse337804_sgrna_lib,
    parse_confirmatory_lines_r14,
    standardize_a549_myllia_to_audit_line,
)
from prepare_r10_lines import build_hgnc_lookup  # noqa: E402
from run_r14_audit import (  # noqa: E402
    DEFAULT_R14_SEED,
    NOVEL_NOT_IN_R12_A549,
    OVERLAP_R12_A549,
    adjudicate_round14,
    check_r14_target_thresholds,
    evaluate_r14_subgroups,
    _write_r14_outputs,
)
from independent_r12_metrics import evaluate_target_line_r12  # noqa: E402


def _write_sparse_group(
    parent: h5py.Group,
    name: str,
    mat: sp.spmatrix,
    encoding: str = "csr",
    float_cast: bool = False,
) -> None:
    grp = parent.create_group(name)
    if encoding == "csr":
        sp_mat = mat.tocsr()
        grp.attrs["encoding-type"] = "csr_matrix"
    elif encoding == "csc":
        sp_mat = mat.tocsc()
        grp.attrs["encoding-type"] = "csc_matrix"
    else:
        raise ValueError(encoding)
    grp.attrs["encoding-version"] = "0.1.0"
    grp.attrs["shape"] = np.asarray(sp_mat.shape, dtype=np.int64)
    dt = np.float32 if float_cast else sp_mat.data.dtype
    grp.create_dataset("data", data=np.asarray(sp_mat.data, dtype=dt))
    grp.create_dataset("indices", data=np.asarray(sp_mat.indices, dtype=np.int32))
    grp.create_dataset("indptr", data=np.asarray(sp_mat.indptr, dtype=np.int64))


def _write_categorical(group: h5py.Group, col: str, values: Sequence[str], legacy: bool = False) -> None:
    cats = sorted(set(values))
    cat_to_code = {c: i for i, c in enumerate(cats)}
    codes = np.asarray([cat_to_code[v] for v in values], dtype=np.int32)
    cats_arr = np.asarray(cats, dtype=object)
    dt_str = h5py.string_dtype(encoding="utf-8")
    if legacy:
        cc = group.require_group("__categories")
        cc.create_dataset(col, data=cats_arr, dtype=dt_str)
        group.create_dataset(col, data=codes)
    else:
        cg = group.create_group(col)
        cg.attrs["encoding-type"] = "categorical"
        cg.attrs["encoding-version"] = "0.2.0"
        cg.attrs["ordered"] = False
        cg.create_dataset("categories", data=cats_arr, dtype=dt_str)
        cg.create_dataset("codes", data=codes)


def _write_synthetic_h5ad(
    path: Path,
    mat: np.ndarray,
    barcodes: Sequence[str],
    var_ids: Sequence[str],
    var_syms: Sequence[str],
    obs_cols: Dict[str, Sequence[object]],
    var_extra: Optional[Dict[str, Sequence[str]]] = None,
    x_encoding: str = "csr",
    x_is_float: bool = False,
    raw_x_mat: Optional[np.ndarray] = None,
    raw_var_ids: Optional[Sequence[str]] = None,
    raw_var_syms: Optional[Sequence[str]] = None,
    layers: Optional[Dict[str, Tuple[np.ndarray, str, bool]]] = None,
    obsm: Optional[Dict[str, Tuple[np.ndarray, str]]] = None,
    uns_guide_names: Optional[Sequence[str]] = None,
    categorical_cols: Sequence[str] = ("sgrna_symbol",),
    legacy_cat: bool = False,
) -> None:
    dt_str = h5py.string_dtype(encoding="utf-8")
    with h5py.File(path, "w") as f:
        if x_encoding == "dense":
            dt = np.float32 if x_is_float else mat.dtype
            f.create_dataset("X", data=np.asarray(mat, dtype=dt))
        else:
            _write_sparse_group(f, "X", sp.csr_matrix(mat), encoding=x_encoding, float_cast=x_is_float)

        obs = f.create_group("obs")
        obs.attrs["_index"] = "_index"
        obs.create_dataset("_index", data=np.asarray(barcodes, dtype=object), dtype=dt_str)
        col_order = list(obs_cols.keys())
        obs.attrs["column-order"] = np.asarray(col_order, dtype=object)
        for ck, vals in obs_cols.items():
            if ck in categorical_cols:
                _write_categorical(obs, ck, [str(x) for x in vals], legacy=legacy_cat)
            elif all(isinstance(v, (int, float, np.integer, np.floating)) for v in vals):
                obs.create_dataset(ck, data=np.asarray(vals, dtype=np.float64))
            else:
                obs.create_dataset(ck, data=np.asarray([str(x) for x in vals], dtype=object), dtype=dt_str)

        var = f.create_group("var")
        var.attrs["_index"] = "_index"
        var.create_dataset("_index", data=np.asarray(var_ids, dtype=object), dtype=dt_str)
        var.create_dataset("gene_symbol", data=np.asarray(var_syms, dtype=object), dtype=dt_str)
        if var_extra:
            for vk, vv in var_extra.items():
                var.create_dataset(vk, data=np.asarray(vv, dtype=object), dtype=dt_str)

        if raw_x_mat is not None:
            raw_g = f.create_group("raw")
            _write_sparse_group(raw_g, "X", sp.csc_matrix(raw_x_mat), encoding="csc", float_cast=False)
            rvar = raw_g.create_group("var")
            rv_ids = list(raw_var_ids or var_ids)
            rv_syms = list(raw_var_syms or var_syms)
            rvar.attrs["_index"] = "_index"
            rvar.create_dataset("_index", data=np.asarray(rv_ids, dtype=object), dtype=dt_str)
            rvar.create_dataset("gene_symbol", data=np.asarray(rv_syms, dtype=object), dtype=dt_str)

        if layers:
            lg = f.create_group("layers")
            for lk, (lmat, lenc, lfloat) in layers.items():
                if lenc == "dense":
                    dt = np.float32 if lfloat else lmat.dtype
                    lg.create_dataset(lk, data=np.asarray(lmat, dtype=dt))
                else:
                    _write_sparse_group(lg, lk, sp.csr_matrix(lmat), encoding=lenc, float_cast=lfloat)

        if obsm:
            mg = f.create_group("obsm")
            for mk, (marr, menc) in obsm.items():
                if menc == "dense":
                    mg.create_dataset(mk, data=marr)
                else:
                    _write_sparse_group(mg, mk, sp.csr_matrix(marr), encoding=menc, float_cast=False)

        if uns_guide_names is not None:
            ug = f.create_group("uns")
            ug.create_dataset("sgrna_ids", data=np.asarray(list(uns_guide_names), dtype=object), dtype=dt_str)


def _get_real_gene_pool(n_needed: int = 260) -> List[Tuple[str, str]]:
    """Pull genuine HGNC approved symbols (excluding overlap/novel target names cleanly except extras) from DEFAULT_HGNC_TSV."""
    hgnc = build_hgnc_lookup(DEFAULT_HGNC_TSV)
    reserved = set(OVERLAP_R12_A549) | set(NOVEL_NOT_IN_R12_A549)
    picked: List[Tuple[str, str]] = []
    # Include the 10 synthetic perturbation target genes first so perturbations exist in expression matrix
    for g in list(OVERLAP_R12_A549[:5]) + list(NOVEL_NOT_IN_R12_A549[:5]):
        eid = hgnc.sym_to_ens.get(g, f"ENSG_FAKE_{g}")
        picked.append((eid or f"ENSG_{g}", g))
    for sym in sorted(hgnc.approved):
        if sym in reserved or sym.startswith(("MT-", "RPS", "RPL")) or "-" in sym:
            continue
        eid = hgnc.sym_to_ens.get(sym, "")
        if not eid.startswith("ENSG"):
            continue
        picked.append((eid, sym))
        if len(picked) >= n_needed:
            break
    return picked


def _write_synthetic_sgrna_lib(
    path: Path,
    target_genes: Sequence[Tuple[str, str]],
) -> List[Dict[str, str]]:
    """Write GSE337804-format sgRNA_lib.csv.gz with 20 ctrl_NT, 20 ctrl_T, and 3 guides/gene."""
    rows: List[Dict[str, str]] = []
    for k in range(1, 21):
        rows.append({
            "id": f"non-targeting_{k}",
            "symbol": "non-targeting",
            "rank": str(k),
            "class": "ctrl_NT",
            "entrez": "0",
            "guide": f"GNT{k:017d}",
            "count": "2000",
        })
    for k in range(1, 21):
        rows.append({
            "id": f"gene_desert_{k}",
            "symbol": "gene_desert",
            "rank": str(k),
            "class": "ctrl_T",
            "entrez": "0",
            "guide": f"GGD{k:017d}",
            "count": "3000",
        })
    for raw_sym, _ in target_genes:
        for rk in (1, 2, 3):
            rows.append({
                "id": f"{raw_sym}_{rk}",
                "symbol": raw_sym,
                "rank": str(rk),
                "class": "gene",
                "entrez": "1",
                "guide": f"G{raw_sym[:3].upper()}{rk:015d}",
                "count": "5000",
            })

    with gzip.open(path, "wt", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["id", "symbol", "rank", "class", "entrez", "guide", "count"])
        w.writeheader()
        w.writerows(rows)
    return rows


class TestRound14Audit(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.hgnc = build_hgnc_lookup(DEFAULT_HGNC_TSV)
        cls.gene_pool = _get_real_gene_pool(240)

    def test_blinding_gates_and_conf_parser(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            # synthetic_test must block round14/shared paths
            with self.assertRaises(RuntimeError):
                enforce_r14_gate("synthetic_test", extra_dirs=[R14_SHARED_DIR / "raw"])

            # dev phase must refuse a549_myllia
            with self.assertRaises(RuntimeError):
                enforce_r14_gate("dev_r12", target_lines=[R14_TARGET_KEY])

            # reader_xcheck_r14 requires BOTH PREREG_G_COMMIT.txt AND R14_LINES_AUTHORIZED.txt
            (tdir / "PREREG_G_COMMIT.txt").write_text("d34560d\n", encoding="utf-8")
            with self.assertRaises(RuntimeError):
                enforce_r14_gate("reader_xcheck_r14", target_lines=[R14_TARGET_KEY], status_dir=tdir)
            (tdir / "R14_LINES_AUTHORIZED.txt").write_text("AUTHORIZED_2026\n", encoding="utf-8")
            meta_rx = enforce_r14_gate("reader_xcheck_r14", target_lines=[R14_TARGET_KEY], status_dir=tdir)
            self.assertEqual(meta_rx["prereg_commit"], "d34560d")
            self.assertEqual(meta_rx["r14_lines_authorized"], "AUTHORIZED_2026")

            # holdout_r14 requires CODE_FREEZE_R14.txt and CONFIRMATORY_LINES_R14.txt
            with self.assertRaises(RuntimeError):
                enforce_r14_gate("holdout_r14", target_lines=[R14_TARGET_KEY], status_dir=tdir)
            (tdir / "CODE_FREEZE_R14.txt").write_text("freeze_sha\n", encoding="utf-8")
            with self.assertRaises(RuntimeError):
                enforce_r14_gate("holdout_r14", target_lines=[R14_TARGET_KEY], status_dir=tdir)
            conf_p = tdir / "CONFIRMATORY_LINES_R14.txt"
            conf_p.write_text("CONFIRMATORY_LINES=a549_myllia\n", encoding="utf-8")
            meta_h = enforce_r14_gate("holdout_r14", target_lines=[R14_TARGET_KEY], status_dir=tdir)
            self.assertEqual(meta_h["code_freeze"], "freeze_sha")

            parsed_ok = parse_confirmatory_lines_r14(conf_p)
            self.assertEqual(parsed_ok["k_conf"], 1)
            self.assertEqual(parsed_ok["confirmatory_lines"], [R14_TARGET_KEY])

            conf_fail = tdir / "CONFIRMATORY_FAIL.txt"
            conf_fail.write_text("A549_MYLLIA_ROLE=EXCLUDED_FEASIBILITY\n", encoding="utf-8")
            parsed_bad = parse_confirmatory_lines_r14(conf_fail)
            self.assertEqual(parsed_bad["k_conf"], 0)
            self.assertTrue(parsed_bad["excluded_feasibility"])

    def test_slot_verification_and_multi_integer_conflict_stops(self) -> None:
        v_ids = [g[0] for g in self.gene_pool[:20]]
        v_syms = [g[1] for g in self.gene_pool[:20]]
        obs = {
            "sgrna_id": ["non-targeting_1"] * 4 + ["HDAC3_1"] * 4,
            "sgrna_symbol": ["non-targeting"] * 4 + ["HDAC3"] * 4,
        }
        bcs = [f"BC_{i}" for i in range(8)]
        int_mat = np.full((8, 20), 5, dtype=np.int32)
        float_mat = np.log1p(int_mat).astype(np.float32)

        with tempfile.TemporaryDirectory() as tmp:
            raw_dir = Path(tmp) / "raw"
            adir = raw_dir / "a549_myllia"
            adir.mkdir(parents=True)
            lib_p = adir / "GSE337804_sgRNA_lib.csv.gz"
            _write_synthetic_sgrna_lib(lib_p, [("HDAC3", "HDAC3"), ("NOTAREALGENE_XYZ", "NOTAREALGENE_XYZ")])
            parsed_lib = load_gse337804_sgrna_lib(lib_p, self.hgnc)
            self.assertEqual(len(parsed_lib["ntc_guide_ids"]), 20)
            self.assertEqual(len(parsed_lib["excluded_safe_ids"]), 20)
            self.assertIn("NOTAREALGENE_XYZ_1", parsed_lib["unmapped_guide_ids"])
            self.assertNotIn("NOTAREALGENE_XYZ_1", parsed_lib["guide_to_target"])

            ch1 = adir / "GSE337804_A549_ch_1.RDS.h5ad"
            ch2 = adir / "GSE337804_A549_ch_2.RDS.h5ad"

            # Case 1: Float-only matrix -> EXCLUDED_FEASIBILITY
            _write_synthetic_h5ad(ch1, float_mat, bcs, v_ids, v_syms, obs, x_is_float=True)
            _write_synthetic_h5ad(ch2, float_mat, bcs, v_ids, v_syms, obs, x_is_float=True)
            f1 = check_gse337804_h5ad_feasibility(raw_dir, hgnc=self.hgnc, lib_path=lib_p)
            self.assertFalse(f1["feasible"])
            self.assertEqual(f1["status"], "EXCLUDED_FEASIBILITY")

            # Case 2: Multiple inconsistent integer slots (/X vs /layers/counts) -> MUST STOP!
            diff_mat = int_mat.copy()
            diff_mat[6, 3] += 10
            _write_synthetic_h5ad(
                ch1,
                int_mat,
                bcs,
                v_ids,
                v_syms,
                obs,
                layers={"counts": (diff_mat, "csr", False)},
            )
            _write_synthetic_h5ad(ch2, int_mat, bcs, v_ids, v_syms, obs)
            f2 = check_gse337804_h5ad_feasibility(raw_dir, hgnc=self.hgnc, lib_path=lib_p)
            self.assertFalse(f2["feasible"])
            self.assertEqual(f2["status"], "EXCLUDED_FEASIBILITY")
            self.assertIn("MULTIPLE_INCONSISTENT_INTEGER_SLOTS", f2["reason"])
            with self.assertRaises(RuntimeError):
                build_a549_myllia_audit_parts(raw_dir, hgnc=self.hgnc, lib_path=lib_p)

            # Case 3: Inconsistent integer shape/var between /raw/X and /X -> MUST STOP!
            _write_synthetic_h5ad(
                ch1,
                int_mat,
                bcs,
                v_ids,
                v_syms,
                obs,
                raw_x_mat=int_mat[:, :18],
                raw_var_ids=v_ids[:18],
                raw_var_syms=v_syms[:18],
            )
            f3 = check_gse337804_h5ad_feasibility(raw_dir, hgnc=self.hgnc, lib_path=lib_p)
            self.assertFalse(f3["feasible"])
            self.assertIn("MULTIPLE_INCONSISTENT_INTEGER_SLOTS", f3["reason"])

            # Case 4: Exact duplicate integer slots (/raw/X CSC == /X CSR) + ignored nascent_counts & log1p -> PASS!
            nascent_sub = np.maximum(int_mat - 2, 0).astype(np.int32)
            _write_synthetic_h5ad(
                ch1,
                int_mat,
                bcs,
                v_ids,
                v_syms,
                obs,
                raw_x_mat=int_mat,
                raw_var_ids=v_ids,
                raw_var_syms=v_syms,
                layers={
                    "nascent_counts": (nascent_sub, "csr", False),
                    "log1p": (float_mat, "dense", True),
                },
            )
            _write_synthetic_h5ad(ch2, int_mat, bcs, v_ids, v_syms, obs, x_encoding="dense")
            f4 = check_gse337804_h5ad_feasibility(raw_dir, hgnc=self.hgnc, lib_path=lib_p)
            self.assertTrue(f4["feasible"], msg=f4.get("reason"))
            self.assertIn("MULTIPLE_IDENTICAL_INTEGER_SLOTS_VERIFIED", f4["channels"][0]["slot_consistency"])

            # Case 5: Non-integer hidden beyond first 200k entries (last entry is float 1.5 -> full-scan catches it!)
            edge_mat = int_mat.astype(np.float32)
            edge_mat[-1, -1] = 1.5
            _write_synthetic_h5ad(ch1, edge_mat, bcs, v_ids, v_syms, obs, x_is_float=True)
            f5 = check_gse337804_h5ad_feasibility(raw_dir, hgnc=self.hgnc, lib_path=lib_p)
            self.assertFalse(f5["feasible"])

    def test_obs_filtering_ctrl_t_exclusion_qc_no_pctmt_and_pooling(self) -> None:
        # 220 human genes + MT-CO1 + Cas9 + ENSMUSG00000001
        genes = self.gene_pool[:220] + [("ENSG00000198804", "MT-CO1"), ("CAS9", "Cas9"), ("ENSMUSG00000001", "Gapdh")]
        v_ids = [g[0] for g in genes]
        v_syms = [g[1] for g in genes]
        mt_col = 220

        # Resolve legacy prev/alias test: pull one real previous symbol in hgnc mapping to e.g. PTEN or HDAC3
        alias_raw = next(
            (k for k, v in self.hgnc.prev_to_sym.items() if v in {"HDAC3", "PTEN"}),
            None,
        )
        tlist = [("HDAC3", "HDAC3"), ("PTEN", "PTEN")]
        if alias_raw:
            tlist.append((alias_raw, self.hgnc.prev_to_sym[alias_raw]))

        rng = np.random.default_rng(123)
        with tempfile.TemporaryDirectory() as tmp:
            raw_dir = Path(tmp) / "raw"
            adir = raw_dir / "a549_myllia"
            adir.mkdir(parents=True)
            lib_p = adir / "GSE337804_sgRNA_lib.csv.gz"
            _write_synthetic_sgrna_lib(lib_p, tlist)

            for ch_num, fname, shift_base in (
                (1, "GSE337804_A549_ch_1.RDS.h5ad", 0),
                (2, "GSE337804_A549_ch_2.RDS.h5ad", 5),
            ):
                rows_gid: List[str] = []
                rows_sym: List[str] = []
                rows_np: List[int] = []
                rows_expr: List[np.ndarray] = []

                def _add_cell(gid: str, sym: str, nperts: int, total_h: int = 400, mt_u: int = 10, n_expr_genes: int = 215, boost_col: Optional[int] = None) -> None:
                    vec = np.zeros(len(genes), dtype=np.int32)
                    vec[:n_expr_genes] = 1
                    if boost_col is not None:
                        vec[boost_col] += 25
                    vec[mt_col] = mt_u
                    vec[221] = 15  # Cas9 (must be stripped before RNA UMI / feature count)
                    vec[222] = 15  # ENSMUSG (must be stripped)
                    # Add remainder UMIs on random expressed human genes
                    rem = max(0, total_h - int(vec[:220].sum()))
                    if rem > 0:
                        idx_r = rng.integers(0, n_expr_genes, size=rem)
                        np.add.at(vec, idx_r, 1)
                    rows_gid.append(gid)
                    rows_sym.append(sym)
                    rows_np.append(nperts)
                    rows_expr.append(vec)

                # 18 ctrl_NT valid cells
                for k in range(18):
                    _add_cell(f"non-targeting_{(k % 20) + 1}", "non-targeting", 1, total_h=450)
                # 6 ctrl_T safe-targeting cells (strong shift on col 20 -> MUST BE EXCLUDED from both NTC & targets)
                for k in range(6):
                    _add_cell(f"gene_desert_{(k % 20) + 1}", "gene_desert", 1, total_h=450, boost_col=20)
                # High %MT singlet (>75% MT UMIs, >=300 human UMIs & >=200 genes -> MUST BE RETAINED, max_pct_mt=None)
                _add_cell("HDAC3_1", "HDAC3", 1, total_h=450, mt_u=2000, n_expr_genes=215, boost_col=0)
                # 9 valid cells across 3 guides of HDAC3 (3 per guide)
                for g_i in (1, 2, 3):
                    for _ in range(3):
                        _add_cell(f"HDAC3_{g_i}", "HDAC3", 1, total_h=450 + shift_base * 10, boost_col=15)
                # 6 valid cells of PTEN (+ 1 legacy alias guide if available)
                for g_i in (1, 2, 3):
                    for _ in range(2):
                        _add_cell(f"PTEN_{g_i}", "PTEN", 1, total_h=450, boost_col=16)
                if alias_raw:
                    _add_cell(f"{alias_raw}_1", alias_raw, 1, total_h=450, boost_col=16)
                # Excluded low-UMI (<300), low-gene (<200), multi-guide, nperts=2, unassigned
                _add_cell("HDAC3_1", "HDAC3", 1, total_h=260, mt_u=0, n_expr_genes=205)    # fails >=300 human UMI
                _add_cell("HDAC3_1", "HDAC3", 1, total_h=500, mt_u=0, n_expr_genes=150)    # fails >=200 human genes
                _add_cell("HDAC3_1,PTEN_1", "HDAC3,PTEN", 1, total_h=500)
                _add_cell("HDAC3_1", "HDAC3", 2, total_h=500)
                _add_cell("unassigned", "unassigned", 0, total_h=500)

                mat = np.vstack(rows_expr)
                bcs = [f"CH{ch_num}_BC_{i:04d}" for i in range(len(rows_gid))]
                _write_synthetic_h5ad(
                    adir / fname,
                    mat,
                    bcs,
                    v_ids,
                    v_syms,
                    obs_cols={"sgrna_id": rows_gid, "sgrna_symbol": rows_sym, "nperts": rows_np},
                    x_encoding="csr" if ch_num == 1 else "csc",
                    legacy_cat=(ch_num == 2),
                )

            parts, meta = build_a549_myllia_audit_parts(raw_dir, hgnc=self.hgnc, lib_path=lib_p)
            self.assertEqual(len(parts), 2)
            self.assertEqual(parts[0].batches[0], "ch_1")
            self.assertEqual(parts[1].batches[0], "ch_2")
            self.assertTrue(np.all(parts[0].half_flags == 1))
            self.assertTrue(np.all(parts[1].half_flags == 2))
            self.assertIsNone(meta["max_pct_mt"])
            self.assertEqual(meta["counts"]["n_kept_ntc_ctrl_NT"], 36)
            self.assertEqual(meta["counts"]["n_excluded_safe_targeting"], 12)
            # HDAC3 & PTEN exact channel counts (verifying alias resolution + low-UMI/low-gene exclusion)
            exp_hdac3 = 10 + (1 if alias_raw and self.hgnc.prev_to_sym[alias_raw] == "HDAC3" else 0)
            exp_pten = 6 + (1 if alias_raw and self.hgnc.prev_to_sym[alias_raw] == "PTEN" else 0)
            self.assertEqual(sum(1 for pk in parts[0].pert_keys if pk == "HDAC3"), exp_hdac3)
            self.assertEqual(sum(1 for pk in parts[1].pert_keys if pk == "HDAC3"), exp_hdac3)
            self.assertEqual(sum(1 for pk in parts[0].pert_keys if pk == "PTEN"), exp_pten)
            self.assertEqual(sum(1 for pk in parts[1].pert_keys if pk == "PTEN"), exp_pten)
            self.assertEqual(meta["counts"]["n_failed_qc_umi_or_genes"], 4)
            self.assertNotIn("gene_desert", set(parts[0].pert_keys) | set(parts[1].pert_keys))

    def test_prereg_f_raw_guide_fallback_and_end_to_end_s7_subgroups(self) -> None:
        # 1) Test PREREG_F raw guide fallback when .obs guide calls are absent
        genes = self.gene_pool[:215]
        v_ids = [g[0] for g in genes]
        v_syms = [g[1] for g in genes]
        rng = np.random.default_rng(77)

        with tempfile.TemporaryDirectory() as tmp:
            troot = Path(tmp)
            raw_fb = troot / "raw_fb" / "a549_myllia"
            raw_fb.mkdir(parents=True)
            lib_fb = raw_fb / "GSE337804_sgRNA_lib.csv.gz"
            lib_rows = _write_synthetic_sgrna_lib(lib_fb, [("HDAC3", "HDAC3"), ("CTCF", "CTCF")])
            all_gids = [r["id"] for r in lib_rows]
            gid_to_idx = {g: i for i, g in enumerate(all_gids)}

            for ch_num, fname in ((1, "GSE337804_A549_ch_1.RDS.h5ad"), (2, "GSE337804_A549_ch_2.RDS.h5ad")):
                n_c = 50
                x_mat = rng.poisson(6.0, size=(n_c, len(genes))).astype(np.int32)
                g_mat = np.zeros((n_c, len(all_gids)), dtype=np.int32)
                # 20 ctrl_NT, 20 HDAC3 (high guide UMI = 15), 5 ctrl_T (excluded!), 5 doublets (both guides=12 -> excluded!)
                for r in range(20):
                    g_mat[r, gid_to_idx[f"non-targeting_{(r % 20) + 1}"]] = 15
                for r in range(20, 40):
                    g_mat[r, gid_to_idx[f"HDAC3_{(r % 3) + 1}"]] = 15
                for r in range(40, 45):
                    g_mat[r, gid_to_idx[f"gene_desert_{(r % 5) + 1}"]] = 15
                for r in range(45, 50):
                    g_mat[r, gid_to_idx["HDAC3_1"]] = 12
                    g_mat[r, gid_to_idx["CTCF_1"]] = 12

                _write_synthetic_h5ad(
                    raw_fb / fname,
                    x_mat,
                    [f"FB{ch_num}_{i}" for i in range(n_c)],
                    v_ids,
                    v_syms,
                    obs_cols={"channel": [f"ch_{ch_num}"] * n_c},
                    obsm={"sgrna_counts": (g_mat, "csr")},
                    uns_guide_names=all_gids,
                    categorical_cols=(),
                )

            feas_fb = check_gse337804_h5ad_feasibility(raw_fb.parent, hgnc=self.hgnc, lib_path=lib_fb)
            self.assertTrue(feas_fb["feasible"], msg=feas_fb.get("reason"))
            self.assertEqual(feas_fb["guide_call_mode"], "PREREG_F_RAW_GUIDE_FALLBACK")
            parts_fb, meta_fb = build_a549_myllia_audit_parts(raw_fb.parent, hgnc=self.hgnc, lib_path=lib_fb)
            self.assertEqual(meta_fb["counts"]["n_kept_ntc_ctrl_NT"], 40)
            self.assertEqual(meta_fb["counts"]["n_excluded_safe_targeting"], 10)
            self.assertEqual(sum(len([pk for pk in p.pert_keys if pk == "HDAC3"]) for p in parts_fb), 40)

            # 2) End-to-end synthetic standardize -> evaluate_target_line_r12 -> §7 subgroups -> adjudicate_round14
            raw_e2e = troot / "raw_e2e" / "a549_myllia"
            raw_e2e.mkdir(parents=True)
            t_ov = list(OVERLAP_R12_A549[:5])
            t_nov = list(NOVEL_NOT_IN_R12_A549[:5])
            t_all = t_ov + t_nov
            lib_e2e = raw_e2e / "GSE337804_sgRNA_lib.csv.gz"
            _write_synthetic_sgrna_lib(lib_e2e, [(t, t) for t in t_all])

            e2e_genes = self.gene_pool[:240]
            e_ids = [g[0] for g in e2e_genes]
            e_syms = [g[1] for g in e2e_genes]
            rng_e = np.random.default_rng(2026)

            # Deterministic latent perturbation profiles over genes 10..210 so Addendum-B (>=50 AD, Energy p<0.001) passes!
            true_eff = rng_e.normal(0.0, 1.5, size=(len(t_all), len(e2e_genes)))
            true_eff[:, :10] = 0.0

            for ch_num, fname in ((1, "GSE337804_A549_ch_1.RDS.h5ad"), (2, "GSE337804_A549_ch_2.RDS.h5ad")):
                base_mu = np.full(len(e2e_genes), 8.0 + (1.5 if ch_num == 2 else 0.0))
                mat_rows: List[np.ndarray] = []
                gids: List[str] = []
                gsyms: List[str] = []
                for k in range(24):
                    mat_rows.append(rng_e.poisson(base_mu).astype(np.int32))
                    gids.append(f"non-targeting_{(k % 20) + 1}")
                    gsyms.append("non-targeting")
                for tidx, tname in enumerate(t_all):
                    # 18 cells/ch = 36 total pooled across 3 guides
                    lam = np.clip(base_mu * np.exp(true_eff[tidx]), 1.0, 80.0)
                    for r in range(18):
                        mat_rows.append(rng_e.poisson(lam).astype(np.int32))
                        gids.append(f"{tname}_{(r % 3) + 1}")
                        gsyms.append(tname)
                _write_synthetic_h5ad(
                    raw_e2e / fname,
                    np.vstack(mat_rows),
                    [f"E{ch_num}_{i:04d}" for i in range(len(gids))],
                    e_ids,
                    e_syms,
                    obs_cols={"sgrna_id": gids, "sgrna_symbol": gsyms},
                )

            cache_root = troot / "lines_out"
            status_fake = troot / "status"
            status_fake.mkdir()
            tdata, tmeta = standardize_a549_myllia_to_audit_line(
                raw_dir=raw_e2e.parent,
                cache_out_root=cache_root,
                aligned_genes=e_syms,
                lib_path=lib_e2e,
                min_gene_overlap=100,
                n_workers=1,
                phase="synthetic_test",
                status_dir=status_fake,
            )
            self.assertEqual(len(tdata.target_names), 10)
            self.assertEqual(tmeta["pooling_unit"], "canonical_hgnc_gene_symbol")

            # Construct 2 synthetic source AuditLineData sharing exact genes & targets
            def _make_src_line(sname: str, ctrl_scale: float) -> AuditLineData:
                ctrl_c = rng_e.normal(0.0, 0.25, size=(60, len(tdata.genes)))
                mu_lib = np.maximum(np.asarray(tdata.ctrl_mean_libnorm) * ctrl_scale, 0.1)
                return AuditLineData(
                    name=sname,
                    genes=list(tdata.genes),
                    target_names=list(t_all),
                    target_indices=np.arange(len(t_all), dtype=np.int64),
                    deltas=np.asarray(tdata.deltas, dtype=np.float64) * 0.8 + rng_e.normal(0.0, 0.05, size=tdata.deltas.shape),
                    ctrl_cells=ctrl_c,
                    ctrl_mean_libnorm=mu_lib,
                    ctrl_mean_libnorm_full=mu_lib,
                    repr_mask=np.ones(len(t_all), dtype=bool),
                    half_a_deltas=None,
                    half_b_deltas=None,
                    half_valid=None,
                    n_cells=np.full(len(t_all), 40.0),
                    kd_depth_log1p=np.full(len(t_all), -0.5),
                    metadata={"synthetic": True},
                )

            src_list = [_make_src_line("src_a", 0.6), _make_src_line("src_b", 1.4)]

            # Default PREREG_G threshold check (>=50 targets, >=3000 genes) MUST fail on 10 targets / 240 genes!
            thr_strict = check_r14_target_thresholds(tdata, src_list, only_reproducible=False)
            self.assertFalse(thr_strict["passes_thresholds"])

            ev_full = evaluate_target_line_r12(
                target_data=tdata,
                source_list=src_list,
                only_reproducible=False,
                n_boot=40,
                n_perm=10,
                seed=DEFAULT_R14_SEED,
            )
            sub_res = evaluate_r14_subgroups(
                target_data=tdata,
                source_list=src_list,
                only_reproducible=False,
                n_boot=40,
                n_perm=10,
                seed=DEFAULT_R14_SEED,
            )
            self.assertEqual(sub_res["subgroups"]["overlap_r12_a549"]["status"], "OK")
            self.assertEqual(sub_res["subgroups"]["novel_not_in_r12_a549"]["status"], "OK")
            self.assertEqual(sub_res["subgroups"]["overlap_r12_a549"]["n_eval_targets"], 5)
            self.assertEqual(sub_res["subgroups"]["novel_not_in_r12_a549"]["n_eval_targets"], 5)
            self.assertTrue(
                np.isfinite(sub_res["subgroups"]["overlap_r12_a549"]["method_eval"]["fct"]["all_nontarget"]["r_orth_cal_est"])
            )

            # Adjudication aborts cleanly at default thresholds (>=50 targets, >=3000 genes)
            adj_abort = adjudicate_round14(
                {R14_TARGET_KEY: ev_full},
                {"confirmatory_lines": [R14_TARGET_KEY], "a549_myllia_role": "CONFIRMATORY"},
            )
            self.assertEqual(adj_abort["k_conf"], 0)
            self.assertEqual(adj_abort["status"], "CONFIRMATORY_ABORTED_K0")

            # And runs K_conf=1 with mandatory scope notice when thresholds met (simulated via min_targets=8, min_genes=100)
            adj_k1 = adjudicate_round14(
                {R14_TARGET_KEY: ev_full},
                {"confirmatory_lines": [R14_TARGET_KEY], "a549_myllia_role": "CONFIRMATORY"},
                min_targets=8,
                min_genes=100,
            )
            self.assertEqual(adj_k1["k_conf"], 1)
            self.assertEqual(adj_k1["status"], "CONFIRMATORY_SINGLE_LINE_K1")
            self.assertTrue(adj_k1["same_line_replication_only"])
            self.assertIn("A549", adj_k1["mandatory_disclosure"])

            out_d = troot / "out_reports"
            _write_r14_outputs(
                {
                    "per_line": {R14_TARGET_KEY: ev_full},
                    "subgroups_s7": {R14_TARGET_KEY: sub_res},
                    "adjudication": adj_k1,
                },
                out_d,
                "synthetic_smoke",
            )
            self.assertTrue((out_d / "audit_r14_synthetic_smoke.json").exists())
            self.assertTrue((out_d / "audit_r14_lines_synthetic_smoke.csv").exists())
            self.assertTrue((out_d / "audit_r14_paired_he1b_synthetic_smoke.csv").exists())
            self.assertTrue((out_d / "audit_r14_subgroups_synthetic_smoke.csv").exists())
            loaded_j = json.loads((out_d / "audit_r14_synthetic_smoke.json").read_text(encoding="utf-8"))
            self.assertEqual(loaded_j["adjudication"]["k_conf"], 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
