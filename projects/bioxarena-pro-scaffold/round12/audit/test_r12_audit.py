#!/usr/bin/env python3
"""
Clean-Room Unit & Synthetic Fixture Test Suite for `round12/audit/` (A-T3).

Strictly NEVER reads `round12/src/`, `round12/results/`, or any held-out Round 12 data.
Run strictly via:
    `python3 -I round12/audit/test_r12_audit.py -v`
"""

from __future__ import annotations

import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

AUDIT_R12_DIR = Path(__file__).resolve().parent
REPO_ROOT = AUDIT_R12_DIR.parents[1]
R11_AUDIT_DIR = REPO_ROOT / "round11" / "audit"
R11_EXP_DIR = R11_AUDIT_DIR / "exploratory_check"
R10_SCRIPTS_DIR = REPO_ROOT / "round10" / "shared" / "scripts"
for _p in (str(AUDIT_R12_DIR), str(R11_AUDIT_DIR), str(R11_EXP_DIR), str(R10_SCRIPTS_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from exploratory_core import evaluate_cross_mask_split_minus_all  # noqa: E402
from independent_r11_loader import AuditLineData  # noqa: E402
from independent_r11_predictors import align_source_and_target, build_predictions, fct_from_fc_bar  # noqa: E402
from independent_r12_loader import (  # noqa: E402
    build_a549_ko_audit_parts,
    build_hek293_xucao_audit_parts,
    check_hek293_xucao_feasibility,
    enforce_r12_gate,
    parse_confirmatory_lines_r12,
)
from independent_r12_metrics import (  # noqa: E402
    adjudicate_round12,
    align_source_and_target_r12,
    evaluate_target_line_r12,
)
from prepare_r10_lines import HGNCLookup  # noqa: E402
from run_n5_audit import _eval_n5_from_moments, _pred_n5_moments, prepare_n5_per_pert_data  # noqa: E402


def _make_tiny_hgnc(genes: list[str]) -> HGNCLookup:
    return HGNCLookup(
        approved=set(genes),
        ens_to_sym={f"ENSG{i:011d}": g for i, g in enumerate(genes, start=1)},
        sym_to_ens={g: f"ENSG{i:011d}" for i, g in enumerate(genes, start=1)},
        prev_to_sym={},
        alias_to_sym={},
    )


def _write_synthetic_h5ad(
    path: Path,
    mat: np.ndarray,
    encoding: str,
    use_layers_counts: bool,
    obs_df: pd.DataFrame,
    var_names: list[str],
    obs_index: list[str],
    float_x_only: bool = False,
) -> None:
    n_c, n_g = mat.shape
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as hf:
        def _write_mat(grp_or_ds_parent: h5py.Group, name: str, arr: np.ndarray) -> None:
            if encoding == "dense":
                grp_or_ds_parent.create_dataset(name, data=arr)
                return
            sp_m = sp.csr_matrix(arr) if encoding == "csr" else sp.csc_matrix(arr)
            g = grp_or_ds_parent.create_group(name)
            g.attrs["encoding-type"] = f"{encoding}_matrix"
            g.attrs["shape"] = np.asarray([n_c, n_g], dtype=np.int64)
            g.create_dataset("data", data=sp_m.data)
            g.create_dataset("indices", data=sp_m.indices.astype(np.int32))
            g.create_dataset("indptr", data=sp_m.indptr.astype(np.int64))

        if use_layers_counts:
            x_log = np.log1p(mat.astype(np.float32))
            _write_mat(hf, "X", x_log)
            lg = hf.create_group("layers")
            _write_mat(lg, "counts", mat.astype(np.int32))
        elif float_x_only:
            _write_mat(hf, "X", np.log1p(mat.astype(np.float32)) + 0.12345)
        else:
            _write_mat(hf, "X", mat.astype(np.int32))

        obs = hf.create_group("obs")
        obs.create_dataset("_index", data=np.asarray(obs_index, dtype="S"))
        obs.attrs["column-order"] = np.asarray(list(obs_df.columns), dtype="S")
        for col in obs_df.columns:
            s = obs_df[col]
            if np.issubdtype(s.dtype, np.number):
                obs.create_dataset(col, data=s.to_numpy())
            else:
                obs.create_dataset(col, data=np.asarray([str(x).encode("utf-8") for x in s], dtype="S"))

        var = hf.create_group("var")
        var.create_dataset("_index", data=np.asarray(var_names, dtype="S"))
        var.create_dataset("gene_symbol", data=np.asarray([v.encode("utf-8") for v in var_names], dtype="S"))


class TestRound12CleanRoomAudit(unittest.TestCase):
    def test_1_blinding_gates_and_confirmatory_parser(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            sdir = Path(td) / "status"
            sdir.mkdir(parents=True)
            (sdir / "PREREG_E_COMMIT.txt").write_text("4c606f4 LOCKED\n", encoding="utf-8")

            # holdout_r12 must refuse before CODE_FREEZE_R12.txt exists
            with self.assertRaises(RuntimeError):
                enforce_r12_gate("holdout_r12", status_dir=sdir)

            # dev_r10 must refuse Round 12 keys and Round 12 shared paths
            with self.assertRaises(RuntimeError):
                enforce_r12_gate("dev_r10", target_lines=["a549_ko"], status_dir=sdir)
            with self.assertRaises(RuntimeError):
                enforce_r12_gate(
                    "dev_r10",
                    target_lines=["dld1"],
                    status_dir=sdir,
                    extra_dirs=[REPO_ROOT / "round12" / "shared" / "raw"],
                )

            # Also verify live repository CLI hard blocks `--phase holdout_r12` cleanly right now
            proc = subprocess.run(
                [sys.executable, "-I", str(AUDIT_R12_DIR / "run_r12_audit.py"), "--phase", "holdout_r12"],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("PREREG_E hard gate blocked", proc.stderr)

            # Test CONFIRMATORY_LINES_R12 parser across K_conf = 0, 1 (EXCLUDED_FEASIBILITY_Q3 / SECONDARY_ONLY), 2
            p_conf = sdir / "CONFIRMATORY_LINES_R12.txt"
            p_conf.write_text(
                "CONFIRMATORY_LINES = a549_ko\nHEK293_XUCAO_ROLE = EXCLUDED_FEASIBILITY_Q3\n",
                encoding="utf-8",
            )
            g1 = parse_confirmatory_lines_r12(p_conf)
            self.assertEqual(g1["k_conf"], 1)
            self.assertEqual(g1["confirmatory_lines"], ["a549_ko"])
            self.assertTrue(g1["hek293_excluded_feasibility_q3"])
            self.assertEqual(g1["secondary_eval_lines"], [])

            p_conf.write_text(
                "A549_ROLE = FAILED_GATE\nHEK293_ROLE = CONFIRMATORY\nCONFIRMATORY_LINES = hek293_xucao\n",
                encoding="utf-8",
            )
            g0 = parse_confirmatory_lines_r12(p_conf)
            self.assertEqual(g0["k_conf"], 0)
            self.assertEqual(g0["confirmatory_lines"], [])
            self.assertTrue(g0["aborted_mandatory_a549"])

    def test_2_a549_reader_cell_id_cas9_lane_and_sacas9(self) -> None:
        rng = np.random.default_rng(101)
        genes = [f"G{i:03d}" for i in range(40)] + ["Cas9"]
        hgnc = _make_tiny_hgnc(genes[:-1])

        with tempfile.TemporaryDirectory() as td:
            raw_root = Path(td) / "raw"
            adir = raw_root / "a549_ko"
            adir.mkdir(parents=True)

            # Build 48 cells across batch1..batch4 with duplicate cell_barcodes across lanes (cell_id unique!)
            batches = ["batch1"] * 12 + ["batch2"] * 12 + ["batch3"] * 12 + ["batch4"] * 12
            lanes = [f"Lane{(i % 2) + 1}-{batches[i][-1]}" for i in range(48)]
            perts = (["Nontargeting"] * 4 + ["G001"] * 4 + ["G002"] * 4) * 4
            conds = ["Nontargeting" if p == "Nontargeting" else "Targeting" for p in perts]
            cas9_ct = [0 if (i % 4 == 0) else 5 for i in range(48)]

            cdf = pd.DataFrame({
                "cell_id": [f"cid_{i:04d}" for i in range(48)],
                "cell_barcode": [f"BC_{i % 12:02d}" for i in range(48)],  # intentionally collides across batches!
                "batch": batches,
                "Lane": lanes,
                "Condition": conds,
                "Perturbation": perts,
                "sgRNA": [f"sg_{p}_{i % 2}" for i, p in enumerate(perts)],
                "total_counts_CAS9": cas9_ct,
            })
            cdf.to_csv(adir / "GSE345058_cell_annotations.csv.gz", index=False, compression="gzip")

            gdf = pd.DataFrame({
                "ensembl_id": [f"ENSG{i:011d}.3" for i in range(1, 41)] + ["Cas9"],
                "gene_symbol": genes,
                "genome": ["GRCh38"] * 40 + ["CAS9"],
            })
            gdf.to_csv(adir / "GSE345058_gene_annotations.csv.gz", index=False, compression="gzip")

            # Write (genes x cells) tab-delimited MTX to verify orientation & Cas9 column exclusion
            counts = rng.integers(1, 10, size=(48, 41), dtype=np.int32)
            with gzip.open(adir / "GSE345058_raw_counts.mtx.gz", "wt", encoding="utf-8") as mf:
                mf.write("%%MatrixMarket matrix coordinate integer general\n")
                mf.write(f"41 48 {41 * 48}\n")
                for ci in range(48):
                    for gi in range(41):
                        mf.write(f"{gi + 1}\t{ci + 1}\t{int(counts[ci, gi])}\n")

            parts, meta = build_a549_ko_audit_parts(
                raw_dir=raw_root,
                hgnc=hgnc,
                exclude_zero_cas9=False,
                cache_h5_path=Path(td) / "tmp_all.h5",
            )
            p = parts[0]
            self.assertEqual(meta["n_cas9_excluded_features"], 1)
            self.assertEqual(p.n_genes_raw, 40)
            self.assertNotIn("Cas9", p.var_gene_names)
            self.assertEqual(p.n_cells, 48)
            self.assertEqual(int(p.nt_mask.sum()), 16)
            self.assertTrue(all("Lane" in b for b in p.batches))
            # Half A = batch1+batch3 (24 cells), Half B = batch2+batch4 (24 cells)
            self.assertEqual(int((p.half_flags == 1).sum()), 24)
            self.assertEqual(int((p.half_flags == 2).sum()), 24)
            chunk = p.read_chunk_fn(0, 10)
            self.assertEqual(chunk.shape, (10, 40))

            # SA-Cas9 subset must exclude the 12 cells with total_counts_CAS9 == 0
            parts_c9, _ = build_a549_ko_audit_parts(
                raw_dir=raw_root,
                hgnc=hgnc,
                exclude_zero_cas9=True,
                cache_h5_path=Path(td) / "tmp_c9.h5",
            )
            self.assertEqual(parts_c9[0].n_cells, 36)

    def test_3_hek293_e3_feasibility_and_stage1_stage2_splits(self) -> None:
        rng = np.random.default_rng(202)
        genes = [f"G{i:03d}" for i in range(30)] + ["dCas9-KRAB-MeCP2"]
        hgnc = _make_tiny_hgnc(genes[:-1])
        mat = rng.integers(1, 8, size=(120, 31), dtype=np.int32)
        barcodes = [f"CELL_{120 - i:03d}" for i in range(120)]
        obs_base = pd.DataFrame({
            "perturbation": ["control"] * 40 + ["G001"] * 40 + ["G002"] * 40,
            "sgRNA": [f"sg_{i}" for i in range(120)],
            "nperts": np.ones(120, dtype=np.int32),
            "knockdown": np.concatenate([np.zeros(40), np.full(40, 0.15), np.full(40, 0.85)]),  # 0.15 would be dropped by KD>40%!
        })

        with tempfile.TemporaryDirectory() as td:
            raw_root = Path(td) / "raw"
            hdir = raw_root / "hek293_xucao"
            h5_p = hdir / "XuCao2023.h5ad"

            # Case A: Float-only normalized X & no integer layers -> E3 EXCLUDED_FEASIBILITY_Q3
            _write_synthetic_h5ad(h5_p, mat, "csr", False, obs_base, genes, barcodes, float_x_only=True)
            f_bad = check_hek293_xucao_feasibility(h5_p)
            self.assertEqual(f_bad["status"], "EXCLUDED_FEASIBILITY_Q3")
            with self.assertRaises(RuntimeError):
                build_hek293_xucao_audit_parts(raw_root, hgnc)

            # Case B: Integer CSC `/layers/counts` + no replicate col -> Stage-2 stratified seed=42, NO >40% KD filter
            obs_s2 = obs_base.copy()
            _write_synthetic_h5ad(h5_p, mat, "csc", True, obs_s2, genes, barcodes)
            f_ok = check_hek293_xucao_feasibility(h5_p)
            self.assertEqual(f_ok["status"], "FEASIBLE")
            parts_s2, meta_s2 = build_hek293_xucao_audit_parts(raw_root, hgnc, apply_author_kd40=False)
            self.assertEqual(meta_s2["split_branch"], "STAGE2_STRATIFIED_SEED42")
            self.assertEqual(parts_s2[0].n_cells, 120)  # all 120 singlets kept (no KD>40% drop!)
            self.assertEqual(parts_s2[0].n_genes_raw, 30)  # dCas9 excluded

            # Verify determinism of Stage 2 regardless of cell row order
            perm_rows = rng.permutation(120)
            _write_synthetic_h5ad(
                h5_p,
                mat[perm_rows],
                "csc",
                True,
                obs_s2.iloc[perm_rows].reset_index(drop=True),
                genes,
                [barcodes[i] for i in perm_rows],
            )
            parts_s2_perm, _ = build_hek293_xucao_audit_parts(raw_root, hgnc, apply_author_kd40=False)
            self.assertEqual(int((parts_s2_perm[0].half_flags == 1).sum()), int((parts_s2[0].half_flags == 1).sum()))

            # Case C: Stage-1 `.obs` replicate parity (>=2 replicates, >=15 NTC/half)
            obs_s1 = obs_base.copy()
            obs_s1["replicate"] = [f"rep{(i % 4) + 1}" for i in range(120)]
            _write_synthetic_h5ad(h5_p, mat, "csr", False, obs_s1, genes, barcodes)
            parts_s1, meta_s1 = build_hek293_xucao_audit_parts(raw_root, hgnc, apply_author_kd40=False)
            self.assertEqual(meta_s1["split_branch"], "STAGE1_OBS_REPLICATE_PARITY")
            self.assertEqual(parts_s1[0].n_cells, 120)

    def test_4_adjudicator_and_exact_parity_against_frozen_r11(self) -> None:
        # Part 1: Verify adjudication logic across K_conf = 0, 1, 2 and partial failure
        pass_line = {"he1_pass": True, "he2_pass": True}
        fail_c1 = {"he1_pass": False, "he2_pass": True}

        adj_k0 = adjudicate_round12({}, {"confirmatory_lines": [], "a549_role": "FAILED_GATE", "hek293_role": "EXCLUDED_FEASIBILITY_Q3"})
        self.assertEqual(adj_k0["status"], "CONFIRMATORY_ABORTED_K0")

        adj_k1 = adjudicate_round12({"a549_ko": pass_line}, {"confirmatory_lines": ["a549_ko"], "a549_role": "CONFIRMATORY", "hek293_role": "SECONDARY_ONLY"})
        self.assertEqual(adj_k1["status"], "CONFIRMATORY_SINGLE_LINE_K1")
        self.assertTrue(adj_k1["limited_generalizability_notice"])
        self.assertIn("GENERALIZATION", adj_k1["modality_interpretation"])

        adj_k2_fail = adjudicate_round12(
            {"a549_ko": fail_c1, "hek293_xucao": pass_line},
            {"confirmatory_lines": ["a549_ko", "hek293_xucao"], "a549_role": "CONFIRMATORY", "hek293_role": "CONFIRMATORY"},
        )
        self.assertFalse(adj_k2_fail["claim1_confirmed"])
        self.assertTrue(adj_k2_fail["claim2_confirmed"])
        self.assertIn("CONFOUNDED", adj_k2_fail["modality_interpretation"])

        # Part 2: Exact <= 1e-12 numerical parity of N_5 and H-E2b vs frozen R11 exploratory_check
        rng = np.random.default_rng(202612)
        n_g, n_t = 120, 24
        genes = [f"G{i:03d}" for i in range(n_g)]
        targets = genes[:n_t]

        def _mk_line(name: str, scale: float) -> AuditLineData:
            ctrl_m = rng.uniform(0.2, 5.0, size=n_g)
            ctrl_cells = rng.normal(loc=ctrl_m, scale=0.4, size=(40, n_g))
            deltas = rng.normal(0.0, scale, size=(n_t, n_g))
            ha = deltas + rng.normal(0.0, 0.08, size=(n_t, n_g))
            hb = deltas + rng.normal(0.0, 0.08, size=(n_t, n_g))
            return AuditLineData(
                name=name,
                genes=genes,
                target_names=targets,
                target_indices=np.arange(n_t, dtype=np.int64),
                deltas=deltas,
                ctrl_cells=ctrl_cells,
                ctrl_mean_libnorm=ctrl_m,
                ctrl_mean_libnorm_full=ctrl_m,
                repr_mask=np.ones(n_t, dtype=bool),
                half_a_deltas=ha,
                half_b_deltas=hb,
                half_valid=np.ones(n_t, dtype=bool),
                n_cells=np.full(n_t, 50.0),
                kd_depth_log1p=np.full(n_t, 0.6),
                metadata={},
            )

        s1 = _mk_line("src1", 0.35)
        s2 = _mk_line("src2", 0.35)
        tgt = _mk_line("tgt_syn", 0.35)
        # Inject true multiplicative signal so FCT has positive signal
        fc_true = np.exp(0.5 * (s1.deltas + s2.deltas))
        sig = fct_from_fc_bar(fc_true, tgt.ctrl_mean_libnorm) + rng.normal(0.0, 0.05, size=(n_t, n_g))
        tgt.half_a_deltas = sig + rng.normal(0.0, 0.04, size=(n_t, n_g))
        tgt.half_b_deltas = sig + rng.normal(0.0, 0.04, size=(n_t, n_g))
        tgt.deltas = 0.5 * (tgt.half_a_deltas + tgt.half_b_deltas)

        cohort = align_source_and_target([s1, s2], tgt, only_reproducible=True)
        preds, fc_bar, grand_mean, _ = build_predictions(cohort, include_sa3_anchors=True, include_exploratory_e=True)
        t_mean = preds["transferred_mean"]

        out_r12 = evaluate_target_line_r12(
            target_data=tgt,
            source_list=[s1, s2],
            n_boot=200,
            n_perm=20,
            seed=20261201,
            include_n5=True,
            include_n5_split=True,
            include_r19=False,
            include_s7=False,
        )

        # Reference check against exploratory_core & run_n5_audit
        ref_sa8 = evaluate_cross_mask_split_minus_all(
            cohort.delta_obs,
            preds["e1_const_med"],
            t_mean,
            cohort.target_cols,
            cohort.half_a_obs,
            cohort.half_b_obs,
            cohort.half_valid,
            top_k=50,
            n_boot=200,
            seed=20261201,
        )
        self.assertAlmostEqual(out_r12["sa8_e1_split_minus_all"]["diff_est"], ref_sa8["diff_est"], places=12)
        self.assertAlmostEqual(out_r12["sa8_e1_split_minus_all"]["diff_ci975_lo"], ref_sa8["diff_ci975_lo"], places=12)

        per_p, xx_o, haa_o, hbb_o, hab_o, _ = prepare_n5_per_pert_data(
            cohort.delta_obs, t_mean, grand_mean, fc_bar, cohort.target_cols, cohort.half_a_obs, cohort.half_b_obs, cohort.half_valid
        )
        draws = np.random.default_rng(20261201).integers(0, len(per_p), size=(200, len(per_p)))
        for mn in ("fct", "transfer_mean_global_scale", "source_grand_mean", "e1_const_med"):
            xy, yy_o, yy_c = _pred_n5_moments(per_p, cohort.delta_obs, preds[mn], t_mean)
            ref_n5 = _eval_n5_from_moments(xy, yy_o, yy_c, xx_o, haa_o, hbb_o, hab_o, draws, 200)
            got_n5 = out_r12["n5_all"]["methods"][mn]
            self.assertAlmostEqual(got_n5["r_cone_sp_est"], ref_n5["r_cone_sp_est"], places=12)
            self.assertAlmostEqual(got_n5["eta_est"], ref_n5["eta_est"], places=12)
            if mn in ("transfer_mean_global_scale", "source_grand_mean"):
                self.assertLess(got_n5["eta_est"], 1e-12)
                self.assertAlmostEqual(got_n5["r_cone_sp_est"], 0.0, places=12)

        # Verify joint gene-universe intersection when a source misses target genes (`r9_r10_8`-style)
        s2_trimmed = AuditLineData(
            name=s2.name,
            genes=s2.genes[:-5],
            target_names=s2.target_names,
            target_indices=s2.target_indices,
            deltas=s2.deltas[:, :-5],
            ctrl_cells=s2.ctrl_cells[:, :-5],
            ctrl_mean_libnorm=s2.ctrl_mean_libnorm[:-5],
            ctrl_mean_libnorm_full=s2.ctrl_mean_libnorm_full[:-5] if s2.ctrl_mean_libnorm_full is not None else None,
            repr_mask=s2.repr_mask,
            half_a_deltas=s2.half_a_deltas[:, :-5] if s2.half_a_deltas is not None else None,
            half_b_deltas=s2.half_b_deltas[:, :-5] if s2.half_b_deltas is not None else None,
            half_valid=s2.half_valid,
            n_cells=s2.n_cells,
            kd_depth_log1p=s2.kd_depth_log1p,
            metadata={},
        )
        aligned_trim = align_source_and_target_r12([s1, s2_trimmed], tgt, only_reproducible=True)
        self.assertEqual(len(aligned_trim.eval_genes), n_g - 5)


if __name__ == "__main__":
    unittest.main()
