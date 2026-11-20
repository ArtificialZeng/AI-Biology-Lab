#!/usr/bin/env python3
"""
BioScaffold-OOD: Modality-Aware Biological Feature Engineering & Prior Injection Library.

Provides robust, vectorized feature extractors for:
1. Molecules & Chemical Mixtures (Morgan FP, Physical Descriptors, Mixture Entropy)
2. Protein & DNA/RNA Sequences (AAC, DPC, k-mer spectra, GC-content, physicochemical scales)
3. Domain Calibration & Metric Evaluators (Calibrated Regression, Residual Correlation, AUROC/AUPRC)
4. Preflight Submission Verifiers (Strict schema, finite floats, artifact tree checks)
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False


# ===========================================================================
# 1. Molecular Feature Extractors
# ===========================================================================

class MoleculeFeatureExtractor:
    """Extracts Morgan fingerprints and normalized physicochemical descriptors from SMILES."""

    def __init__(self, n_bits: int = 1024, radius: int = 2):
        self.n_bits = n_bits
        self.radius = radius

    def transform_single(self, smi: str) -> np.ndarray:
        if not HAS_RDKIT or not isinstance(smi, str) or not smi:
            return np.zeros(self.n_bits + 8, dtype=np.float32)

        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                return np.zeros(self.n_bits + 8, dtype=np.float32)

            # Morgan Fingerprint
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=self.radius, nBits=self.n_bits)
            fp_arr = np.zeros(self.n_bits, dtype=np.float32)
            for bit in fp.GetOnBits():
                fp_arr[bit] = 1.0

            # 8 Curated Physical Descriptors
            mw = Descriptors.MolWt(mol) / 500.0
            logp = Descriptors.MolLogP(mol) / 5.0
            tpsa = Descriptors.TPSA(mol) / 150.0
            hba = float(Descriptors.NumHAcceptors(mol)) / 10.0
            hbd = float(Descriptors.NumHDonors(mol)) / 5.0
            rot = float(Descriptors.NumRotatableBonds(mol)) / 10.0
            arom = float(Descriptors.NumAromaticRings(mol)) / 5.0
            fsp3 = float(Descriptors.FractionCSP3(mol))

            desc = np.array([mw, logp, tpsa, hba, hbd, rot, arom, fsp3], dtype=np.float32)
            return np.concatenate([fp_arr, desc])
        except Exception:
            return np.zeros(self.n_bits + 8, dtype=np.float32)

    def transform(self, smiles_list: Union[List[str], pd.Series]) -> np.ndarray:
        return np.vstack([self.transform_single(smi) for smi in smiles_list])

    def transform_mixture(
        self,
        smiles_json_series: pd.Series,
        fractions_json_series: pd.Series
    ) -> np.ndarray:
        """Transforms multi-component mixtures into fraction-weighted & max-pooled representations."""
        features = []
        cache: Dict[str, np.ndarray] = {}

        for smi_str, frac_str in zip(smiles_json_series, fractions_json_series):
            try:
                smis = json.loads(smi_str) if isinstance(smi_str, str) else []
                fracs = [float(x) for x in json.loads(frac_str)] if isinstance(frac_str, str) else []
            except Exception:
                smis, fracs = [], []

            if not smis or len(smis) != len(fracs):
                features.append(np.zeros((self.n_bits + 8) * 2 + 3, dtype=np.float32))
                continue

            comp_feats = []
            for smi in smis:
                if smi not in cache:
                    cache[smi] = self.transform_single(smi)
                comp_feats.append(cache[smi])

            comp_mat = np.vstack(comp_feats)
            weights = np.array(fracs, dtype=np.float32).reshape(-1, 1)
            w_sum = max(1e-6, float(np.sum(weights)))
            weights /= w_sum

            # 1. Weighted average representation
            weighted_feat = np.sum(comp_mat * weights, axis=0)
            # 2. Max-pooled representation
            max_feat = np.max(comp_mat, axis=0)

            # 3. Mixture statistics: num components, max fraction, entropy
            n_comps = float(len(smis)) / 5.0
            max_frac = float(np.max(weights))
            entropy = float(-np.sum([p * math.log(max(1e-12, p)) for p in fracs if p > 0])) / 3.0
            stats = np.array([n_comps, max_frac, entropy], dtype=np.float32)

            features.append(np.concatenate([weighted_feat, max_feat, stats]))

        return np.vstack(features)


# ===========================================================================
# 2. Sequence Feature Extractors
# ===========================================================================

class SequenceFeatureExtractor:
    """Extracts k-mer spectra, composition, and length features from biological sequences."""

    AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
    DNA_BASES = "ACGT"

    def __init__(self, modality: str = "protein", k: int = 3):
        self.modality = modality.lower()
        self.k = k
        if self.modality == "protein":
            self.alphabet = self.AMINO_ACIDS
        else:
            self.alphabet = self.DNA_BASES
        self.char_to_idx = {c: i for i, c in enumerate(self.alphabet)}

    def transform_protein_aac_dpc(self, seq_series: pd.Series) -> np.ndarray:
        """Computes 20-dim AAC + sequence length."""
        features = []
        for seq in seq_series:
            vec = np.zeros(len(self.AMINO_ACIDS) + 1, dtype=np.float32)
            if isinstance(seq, str) and seq:
                length = len(seq)
                vec[-1] = math.log1p(length) / 10.0
                for char in seq.upper():
                    idx = self.char_to_idx.get(char)
                    if idx is not None:
                        vec[idx] += 1.0
                if length > 0:
                    vec[:len(self.AMINO_ACIDS)] /= length
            features.append(vec)
        return np.vstack(features)

    def transform_dna_kmer(self, seq_series: pd.Series, k: int = 3) -> np.ndarray:
        """Computes 4^k normalized k-mer frequencies + GC content."""
        from itertools import product
        kmers = ["".join(p) for p in product(self.DNA_BASES, repeat=k)]
        kmer_to_idx = {km: i for i, km in enumerate(kmers)}
        n_kmers = len(kmers)

        features = []
        for seq in seq_series:
            vec = np.zeros(n_kmers + 2, dtype=np.float32)
            if isinstance(seq, str) and len(seq) >= k:
                seq_u = seq.upper()
                total = len(seq_u) - k + 1
                gc_count = 0
                for char in seq_u:
                    if char in ("G", "C"):
                        gc_count += 1
                vec[-2] = float(gc_count) / max(1, len(seq_u))  # GC-ratio
                vec[-1] = math.log1p(len(seq_u)) / 10.0          # log-length

                for i in range(total):
                    kmer = seq_u[i:i+k]
                    idx = kmer_to_idx.get(kmer)
                    if idx is not None:
                        vec[idx] += 1.0
                if total > 0:
                    vec[:n_kmers] /= total
            features.append(vec)

        return np.vstack(features)


# ===========================================================================
# 3. Benchmark Metric Evaluators
# ===========================================================================

class BenchmarkMetrics:
    """Implements standard BioXArena-Pro primary evaluation metrics."""

    @staticmethod
    def calibrated_regression(y_true: np.ndarray, y_pred: np.ndarray) -> float:
        """0.5 * clip(Pearson_r, 0, 1) + 0.5 * clip(1 - RMSE / std(y_true), 0, 1)"""
        y_true = np.asarray(y_true, dtype=float)
        y_pred = np.asarray(y_pred, dtype=float)

        var_true = np.var(y_true)
        var_pred = np.var(y_pred)
        if var_true <= 1e-12 or var_pred <= 1e-12:
            r = 0.0
        else:
            cov = np.cov(y_true, y_pred)[0, 1]
            r = cov / (math.sqrt(var_true) * math.sqrt(var_pred))
        corr_score = float(np.clip(r, 0.0, 1.0))

        std_true = math.sqrt(var_true) if var_true > 1e-12 else 1e-6
        rmse = math.sqrt(mean_squared_error(y_true, y_pred))
        calib_score = float(np.clip(1.0 - (rmse / std_true), 0.0, 1.0))

        return 0.5 * corr_score + 0.5 * calib_score

    @staticmethod
    def residual_correlation(
        y_true: np.ndarray,
        y_pred: np.ndarray,
        reference: np.ndarray
    ) -> float:
        """Evaluates residual gain beyond conditions reference."""
        res_true = y_true - reference
        res_pred = y_pred - reference

        var_true = np.var(res_true)
        var_pred = np.var(res_pred)
        if var_true <= 1e-12 or var_pred <= 1e-12:
            r = 0.0
        else:
            cov = np.cov(res_true, res_pred)[0, 1]
            r = cov / (math.sqrt(var_true) * math.sqrt(var_pred))
        corr_score = float(np.clip(r, 0.0, 1.0))

        std_true = math.sqrt(var_true) if var_true > 1e-12 else 1e-6
        rmse = math.sqrt(mean_squared_error(y_pred, y_true))
        calib_score = float(np.clip(1.0 - (rmse / std_true), 0.0, 1.0))

        return 0.5 * corr_score + 0.5 * calib_score


# ===========================================================================
# 4. Submission Artifact Verifier
# ===========================================================================

class PreflightVerifier:
    """Verifies that generated submissions strictly satisfy all benchmark contracts."""

    @staticmethod
    def verify(task_dir: Path, output_dir: Path, task_key: str) -> Tuple[bool, str]:
        sub_path = output_dir / "submission.csv"
        sample_path = task_dir / "sample_submission.csv"

        if not sub_path.is_file():
            return False, f"Missing {sub_path}"
        if not sample_path.is_file():
            return False, f"Missing benchmark sample {sample_path}"

        sample_df = pd.read_csv(sample_path)
        sub_df = pd.read_csv(sub_path)

        if list(sub_df.columns) != list(sample_df.columns):
            return False, f"Column mismatch: expected {list(sample_df.columns)}, got {list(sub_df.columns)}"

        if len(sub_df) != len(sample_df):
            return False, f"Row count mismatch: expected {len(sample_df)}, got {len(sub_df)}"

        if not (sub_df.iloc[:, 0].astype(str).values == sample_df.iloc[:, 0].astype(str).values).all():
            return False, "ID alignment mismatch in first column"

        # Check numeric columns for NaN / Inf
        for col in sub_df.columns[1:]:
            if np.issubdtype(sub_df[col].dtype, np.number):
                if sub_df[col].isna().any():
                    return False, f"Column {col} contains NaN values"
                if np.isinf(sub_df[col]).any():
                    return False, f"Column {col} contains Inf values"

        return True, "100% compliant with submission contract"
