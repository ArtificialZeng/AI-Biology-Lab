#!/usr/bin/env python3
"""
Pre-Download Synthetic-Fixture Unit Tests for Pfizer TeloHAEC Reader
(`round17/src/pfizer_telohaec_reader_r17.py`).

Tests strictly on temporary in-memory/tempdir synthetic `.h5ad` fixtures (zero real data,
zero downloads, zero `torch`/`state` imports):
  1. CSR, CSC, and Dense storage parity (`max_abs == 0.0`);
  2. Modern categorical (`codes/categories`) AND legacy AnnData categorical (`__categories`) decoding;
  3. Strict condition separation (`Untreated` basal primary vs `Treated` IL-1b+TNF-a secondary);
  4. Strict `NO-TARGET` NTC vs `SAFE_TARGET` / multi-guide / unassigned exclusion;
  5. Fail-closed non-integer (`EXCLUDED_FEASIBILITY`) and multi-slot conflict (`SLOT_AMBIGUOUS`) guards;
  6. Deterministic hardware-unit split-half (`flask` odd/even parity, `welltag` fallback, `seed` fallback);
  7. End-to-end synthetic `process_r10_line` standardization + NTC raw-count extraction parity.
"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from typing import Dict, List, Optional

import h5py
import numpy as np
import scipy.sparse as sp

ROUND17_SRC_DIR = Path(__file__).resolve().parent
REPO_ROOT = ROUND17_SRC_DIR.parents[1]
R10_SCRIPTS_DIR = REPO_ROOT / "round10" / "shared" / "scripts"
for _p in (str(ROUND17_SRC_DIR), str(R10_SCRIPTS_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from pfizer_telohaec_reader_r17 import (  # noqa: E402
    build_pfizer_telohaec_parts_r17,
    classify_pfizer_condition_r17,
    classify_pfizer_guide_r17,
    extract_pfizer_telohaec_ntc_raw_counts_r17,
)
from prepare_r10_lines import HGNCLookup, process_r10_line  # noqa: E402


def _make_mock_hgnc(genes: List[str]) -> HGNCLookup:
    approved = set(genes)
    sym_to_ens = {g: f"ENSG{100000 + i:011d}" for i, g in enumerate(genes)}
    ens_to_sym = {v: k for k, v in sym_to_ens.items()}
    return HGNCLookup(
        approved=approved,
        ens_to_sym=ens_to_sym,
        prev_to_sym={"AARS": "AARS1"},
        alias_to_sym={},
        sym_to_ens=sym_to_ens,
    )


def _write_str_dataset(grp: h5py.Group, name: str, values: List[str]) -> None:
    dt = h5py.string_dtype(encoding="utf-8")
    grp.create_dataset(name, data=np.asarray(values, dtype=object), dtype=dt)


def _write_categorical(
    grp: h5py.Group,
    name: str,
    values: List[str],
    legacy: bool = False,
) -> None:
    cats = sorted(set(values))
    c_map = {c: i for i, c in enumerate(cats)}
    codes = np.asarray([c_map[v] for v in values], dtype=np.int32)
    if legacy:
        cat_grp = grp.require_group("__categories")
        _write_str_dataset(cat_grp, name, cats)
        grp.create_dataset(name, data=codes)
    else:
        sub = grp.create_group(name)
        sub.attrs["encoding-type"] = "categorical"
        sub.create_dataset("codes", data=codes)
        _write_str_dataset(sub, "categories", cats)


def _write_synthetic_pfizer_h5ad(
    path: Path,
    X_dense: np.ndarray,
    obs_cols: Dict[str, List[str]],
    var_genes: List[str],
    layout: str = "csr",
    legacy_cat: bool = False,
    extra_layer: Optional[Dict[str, np.ndarray]] = None,
) -> None:
    n_obs, n_var = X_dense.shape
    with h5py.File(path, "w") as f:
        if layout == "dense":
            f.create_dataset("X", data=X_dense)
        elif layout == "csr":
            mat = sp.csr_matrix(X_dense)
            xg = f.create_group("X")
            xg.attrs["encoding-type"] = "csr_matrix"
            xg.attrs["shape"] = np.asarray([n_obs, n_var], dtype=np.int64)
            xg.create_dataset("data", data=mat.data)
            xg.create_dataset("indices", data=mat.indices.astype(np.int32))
            xg.create_dataset("indptr", data=mat.indptr.astype(np.int64))
        elif layout == "csc":
            mat = sp.csc_matrix(X_dense)
            xg = f.create_group("X")
            xg.attrs["encoding-type"] = "csc_matrix"
            xg.attrs["shape"] = np.asarray([n_obs, n_var], dtype=np.int64)
            xg.create_dataset("data", data=mat.data)
            xg.create_dataset("indices", data=mat.indices.astype(np.int32))
            xg.create_dataset("indptr", data=mat.indptr.astype(np.int64))
        else:
            raise ValueError(layout)

        if extra_layer:
            lg = f.create_group("layers")
            for lname, lmat in extra_layer.items():
                lg.create_dataset(lname, data=lmat)

        obs = f.create_group("obs")
        obs.attrs["_index"] = "_index"
        _write_str_dataset(obs, "_index", [f"cell_{i:05d}" for i in range(n_obs)])
        for k, vals in obs_cols.items():
            if k in ("condition", "gene_target", "flask", "welltag"):
                _write_categorical(obs, k, vals, legacy=legacy_cat)
            else:
                _write_str_dataset(obs, k, vals)

        var = f.create_group("var")
        var.attrs["_index"] = "_index"
        _write_str_dataset(var, "_index", var_genes)
        _write_str_dataset(var, "gene_name", var_genes)
        _write_str_dataset(var, "ensembl_id", [f"ENSG{100000 + i:011d}" for i in range(n_var)])


class TestPfizerTeloHAECReaderR17(unittest.TestCase):
    def setUp(self) -> None:
        self.genes = [f"GENE{i:04d}" for i in range(260)]
        self.genes[0] = "TP53"
        self.genes[1] = "MYC"
        self.genes[2] = "AARS1"
        self.hgnc = _make_mock_hgnc(self.genes)

    def _build_fixture_arrays(self) -> tuple:
        rng = np.random.default_rng(20261701)
        # 240 cells:
        #   - 120 Untreated: 48 NO-TARGET (12/flask? -> 4 flasks x 20 NTC = 80 NTC + 40 targets + 8 SAFE_TARGET + 6 multi/NA + 6 low-UMI)
        #   - 80 Treated: 40 NO-TARGET + 32 targets + 8 SAFE_TARGET
        rows_meta: List[Dict[str, str]] = []
        flasks = ["F1", "F2", "F3", "F4"]

        # Untreated: 4 flasks x 20 NO-TARGET = 80 NTC
        for fi, fl in enumerate(flasks):
            for j in range(20):
                rows_meta.append(
                    {
                        "condition": "Untreated",
                        "cellID": f"U_NT_{fl}_{j:02d}",
                        "gene_target": "NO-TARGET",
                        "guide": f"NO-TARGET_g{(j % 4) + 1}",
                        "flask": fl,
                        "welltag": f"W{fi + 1}",
                        "low_qc": "0",
                    }
                )
        # Untreated: 4 flasks x 10 TP53 + 4 flasks x 10 AARS (legacy alias -> AARS1)
        for fi, fl in enumerate(flasks):
            for j in range(10):
                rows_meta.append(
                    {
                        "condition": "Untreated",
                        "cellID": f"U_TP53_{fl}_{j:02d}",
                        "gene_target": "TP53",
                        "guide": "TP53_sg1",
                        "flask": fl,
                        "welltag": f"W{fi + 1}",
                        "low_qc": "0",
                    }
                )
            for j in range(10):
                rows_meta.append(
                    {
                        "condition": "Untreated",
                        "cellID": f"U_AARS_{fl}_{j:02d}",
                        "gene_target": "AARS",
                        "guide": "AARS_sg1",
                        "flask": fl,
                        "welltag": f"W{fi + 1}",
                        "low_qc": "0",
                    }
                )
        # Untreated exclusions: 8 SAFE_TARGET, 4 multi-guide, 2 unassigned, 4 low-UMI RNA QC fail
        for j in range(8):
            rows_meta.append(
                {
                    "condition": "Untreated",
                    "cellID": f"U_SAFE_{j:02d}",
                    "gene_target": "SAFE_TARGET",
                    "guide": f"SAFE_TARGET_g{j}",
                    "flask": flasks[j % 4],
                    "welltag": f"W{(j % 4) + 1}",
                    "low_qc": "0",
                }
            )
        for j in range(4):
            rows_meta.append(
                {
                    "condition": "Untreated",
                    "cellID": f"U_MULTI_{j:02d}",
                    "gene_target": "TP53,MYC",
                    "guide": "TP53_sg1;MYC_sg1",
                    "flask": flasks[j % 4],
                    "welltag": f"W{(j % 4) + 1}",
                    "low_qc": "0",
                }
            )
        for j in range(2):
            rows_meta.append(
                {
                    "condition": "Untreated",
                    "cellID": f"U_NA_{j:02d}",
                    "gene_target": "unassigned",
                    "guide": "NA",
                    "flask": flasks[j % 4],
                    "welltag": f"W{(j % 4) + 1}",
                    "low_qc": "0",
                }
            )
        for j in range(4):
            rows_meta.append(
                {
                    "condition": "Untreated",
                    "cellID": f"U_LOWQC_{j:02d}",
                    "gene_target": "NO-TARGET",
                    "guide": "NO-TARGET_g1",
                    "flask": flasks[j % 4],
                    "welltag": f"W{(j % 4) + 1}",
                    "low_qc": "1",
                }
            )

        # Treated (IL-1b + TNF-a): 4 flasks x 16 NO-TARGET + 4 flasks x 8 MYC + 4 SAFE_TARGET
        for fi, fl in enumerate(flasks):
            for j in range(16):
                rows_meta.append(
                    {
                        "condition": "Treated",
                        "cellID": f"T_NT_{fl}_{j:02d}",
                        "gene_target": "NO-TARGET",
                        "guide": "NO-TARGET_g1",
                        "flask": fl,
                        "welltag": f"W{fi + 1}",
                        "low_qc": "0",
                    }
                )
            for j in range(8):
                rows_meta.append(
                    {
                        "condition": "Treated",
                        "cellID": f"T_MYC_{fl}_{j:02d}",
                        "gene_target": "MYC",
                        "guide": "MYC_sg1",
                        "flask": fl,
                        "welltag": f"W{fi + 1}",
                        "low_qc": "0",
                    }
                )
        for j in range(4):
            rows_meta.append(
                {
                    "condition": "Treated",
                    "cellID": f"T_SAFE_{j:02d}",
                    "gene_target": "SAFE_TARGET",
                    "guide": "SAFE_TARGET_g1",
                    "flask": flasks[j % 4],
                    "welltag": f"W{(j % 4) + 1}",
                    "low_qc": "0",
                }
            )

        n_obs = len(rows_meta)
        n_var = len(self.genes)
        X = rng.integers(1, 5, size=(n_obs, n_var), dtype=np.int32)
        for i, rec in enumerate(rows_meta):
            if rec["low_qc"] == "1":
                X[i, :] = 0
                X[i, :10] = 2  # umi=20 < 300, ngenes=10 < 200
            elif rec["gene_target"] == "TP53":
                X[i, 0] = 0
                X[i, 5:60] += 4
            elif rec["gene_target"] == "MYC":
                X[i, 1] = 0
                X[i, 60:115] += 4

        obs_cols = {
            k: [r[k] for r in rows_meta]
            for k in ("condition", "cellID", "gene_target", "guide", "flask", "welltag")
        }
        return X, obs_cols

    def test_01_condition_and_guide_classification_rules(self) -> None:
        self.assertEqual(classify_pfizer_condition_r17("Untreated"), "untreated")
        self.assertEqual(classify_pfizer_condition_r17("Basal"), "untreated")
        self.assertEqual(classify_pfizer_condition_r17("Treated"), "treated")
        self.assertEqual(classify_pfizer_condition_r17("IL-1b+TNF-a"), "treated")
        self.assertEqual(classify_pfizer_condition_r17("mystery"), "other")

        cls_nt, info_nt, _ = classify_pfizer_guide_r17("NO-TARGET", "NO-TARGET_1", self.hgnc)
        self.assertEqual(cls_nt, "ntc")
        self.assertIsNone(info_nt)

        for safe_tok in ("SAFE_TARGET", "safe-targeting", "gene_desert_1", "intergenic_ctrl"):
            cls_s, _, _ = classify_pfizer_guide_r17(safe_tok, "g1", self.hgnc)
            self.assertEqual(cls_s, "safe")

        cls_mg, _, _ = classify_pfizer_guide_r17("TP53,MYC", "g1", self.hgnc)
        self.assertEqual(cls_mg, "exclude")

        cls_aars, info_aars, gk = classify_pfizer_guide_r17("AARS", "AARS_sg1", self.hgnc)
        self.assertEqual(cls_aars, "target")
        assert info_aars is not None
        self.assertEqual(info_aars[1], "AARS1")  # legacy prev_symbol resolved to approved HGNC
        self.assertEqual(gk, "AARS_sg1")

    def test_02_csr_csc_dense_and_legacy_categorical_parity(self) -> None:
        X, obs_cols = self._build_fixture_arrays()
        with tempfile.TemporaryDirectory() as tmpd:
            td = Path(tmpd)
            outputs = {}
            for layout, leg in [("csr", False), ("csc", False), ("dense", False), ("csr", True)]:
                p = td / f"imru_{layout}_leg{int(leg)}.h5ad"
                _write_synthetic_pfizer_h5ad(p, X, obs_cols, self.genes, layout=layout, legacy_cat=leg)
                parts, audit = build_pfizer_telohaec_parts_r17(p, self.hgnc, condition="untreated")
                self.assertEqual(len(parts), 1)
                part = parts[0]
                self.assertEqual(audit["n_safe_targeting_excluded"], 8)
                self.assertEqual(audit["n_multiguide_or_unassigned_excluded"], 6)
                self.assertEqual(audit["n_rna_qc_excluded"], 4)
                self.assertEqual(audit["n_kept_ntc"], 80)
                self.assertEqual(audit["n_kept_targeting_singlets"], 80)
                self.assertEqual(audit["n_kept_half_a_ntc"], 40)
                self.assertEqual(audit["n_kept_half_b_ntc"], 40)
                self.assertTrue(audit["per_half_ntc_ge15"])
                self.assertEqual(audit["half_rule_used"], "flask_lexicographic_odd_even")
                self.assertEqual(audit["half_a_units"], ["F1", "F3"])
                self.assertEqual(audit["half_b_units"], ["F2", "F4"])
                self.assertEqual(audit["n_unique_canon_target_genes"], 2)  # TP53 & AARS1
                mat = part.read_chunk_fn(0, part.n_cells)
                outputs[(layout, leg)] = (mat, part.nt_mask, part.half_flags, part.batches, audit["kept_barcodes_sha256"])

            ref_mat, ref_nt, ref_hf, ref_b, ref_bc = outputs[("csr", False)]
            for key, (m, nt, hf, b, bc) in outputs.items():
                self.assertTrue(np.array_equal(ref_mat, m), f"matrix mismatch for {key}")
                self.assertTrue(np.array_equal(ref_nt, nt), f"nt_mask mismatch for {key}")
                self.assertTrue(np.array_equal(ref_hf, hf), f"half_flags mismatch for {key}")
                self.assertTrue(np.array_equal(ref_b, b), f"batches mismatch for {key}")
                self.assertEqual(ref_bc, bc, f"barcodes hash mismatch for {key}")

    def test_03_condition_separation_and_ntc_extractor(self) -> None:
        X, obs_cols = self._build_fixture_arrays()
        with tempfile.TemporaryDirectory() as tmpd:
            p = Path(tmpd) / "imru_full.h5ad"
            _write_synthetic_pfizer_h5ad(p, X, obs_cols, self.genes, layout="csr")

            parts_u, aud_u = build_pfizer_telohaec_parts_r17(p, self.hgnc, condition="untreated")
            parts_t, aud_t = build_pfizer_telohaec_parts_r17(p, self.hgnc, condition="treated")

            self.assertEqual(aud_u["n_kept_ntc"], 80)
            self.assertEqual(aud_t["n_kept_ntc"], 64)
            self.assertEqual(aud_t["n_kept_targeting_singlets"], 32)
            self.assertEqual(aud_t["n_safe_targeting_excluded"], 4)
            self.assertEqual({v[1] for v in parts_u[0].pert_meta.values()}, {"TP53", "AARS1"})
            self.assertEqual({v[1] for v in parts_t[0].pert_meta.values()}, {"MYC"})

            # NTC raw extractor returns exact integer matrix & condition-prefixed barcodes
            ntc_u = extract_pfizer_telohaec_ntc_raw_counts_r17(p, self.hgnc, condition="untreated")
            self.assertEqual(ntc_u["raw_counts"].shape, (80, len(self.genes)))
            self.assertEqual(len(ntc_u["barcodes"]), 80)
            self.assertTrue(all(bc.startswith("untreated:U_NT_") for bc in ntc_u["barcodes"]))
            self.assertTrue(np.all(ntc_u["raw_counts"] == np.floor(ntc_u["raw_counts"])))

            # End-to-end frozen process_r10_line smoke test inside tempdir
            out_dir = process_r10_line(
                line_name="pfizer_telohaec_basal",
                parts=parts_u,
                hgnc=self.hgnc,
                out_root=Path(tmpd) / "std_lines",
                extra_meta={"adapter_audit": aud_u},
                n_workers=1,
            )
            self.assertTrue((out_dir / "pseudobulk_delta.npz").is_file())
            self.assertTrue((out_dir / "line_metadata.json").is_file())

    def test_04_fail_closed_non_integer_and_slot_ambiguity(self) -> None:
        X, obs_cols = self._build_fixture_arrays()
        with tempfile.TemporaryDirectory() as tmpd:
            td = Path(tmpd)
            # Non-integer log-normalized X -> EXCLUDED_FEASIBILITY
            p_float = td / "bad_float.h5ad"
            _write_synthetic_pfizer_h5ad(p_float, np.log1p(X.astype(np.float32)), obs_cols, self.genes, layout="dense")
            with self.assertRaises(RuntimeError) as ctx1:
                build_pfizer_telohaec_parts_r17(p_float, self.hgnc, condition="untreated")
            self.assertIn("EXCLUDED_FEASIBILITY", str(ctx1.exception))

            # Conflicting integer slots (X vs layers/counts) -> SLOT_AMBIGUOUS
            p_ambig = td / "ambig_slots.h5ad"
            _write_synthetic_pfizer_h5ad(
                p_ambig,
                X,
                obs_cols,
                self.genes,
                layout="csr",
                extra_layer={"counts": X + 1},
            )
            with self.assertRaises(RuntimeError) as ctx2:
                build_pfizer_telohaec_parts_r17(p_ambig, self.hgnc, condition="untreated")
            self.assertIn("SLOT_AMBIGUOUS", str(ctx2.exception))

            # Explicit x_layer resolves ambiguity cleanly
            parts_ok, aud_ok = build_pfizer_telohaec_parts_r17(
                p_ambig, self.hgnc, condition="untreated", x_layer="layers/counts"
            )
            self.assertEqual(aud_ok["x_slot_used"], "layers/counts")
            self.assertEqual(parts_ok[0].n_cells, 160)


if __name__ == "__main__":
    unittest.main()
