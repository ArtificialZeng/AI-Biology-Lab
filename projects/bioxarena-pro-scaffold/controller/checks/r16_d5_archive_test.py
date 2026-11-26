"""Controller test for DEVIATION_R16_D5 (positional gene index in Arc's archived nadig_hepg2 adata_pred.h5ad).

Feeds Arc's own pseudobulk deltas (computed here independently, positional G_model mapping) back into the patched
verify_arc_archived_adata_pred_r16 as the "adapter" predictions:
  1. correct gene names  -> must pass with median and pooled Pearson ~1 on all supplied genes
  2. gene names shuffled -> must fail (the check is not vacuous)
Uses Arc's archive only; no Round 16 prediction is involved.
Usage: round16/.venv/bin/python controller/checks/r16_d5_archive_test.py
"""
import pickle
import sys
from pathlib import Path

import h5py
import numpy as np

P = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(P / "round16" / "src"))
import run_eval_r16 as M  # noqa: E402

ARC = Path("/Users/mac/.cache/r16_arc_archive/ST-SE-Replogle/zeroshot/hepg2/eval_final.ckpt/adata_pred.h5ad")
gm = [str(g) for g in pickle.load(open(P / "round16/shared/models/ST-SE-Replogle/zeroshot/hepg2/var_dims.pkl", "rb"))["gene_names"]]
with h5py.File(ARC, "r") as f:
    pert = np.array(M._read_h5ad_obs_col_r16(f["obs"], "gene"))
    perts = sorted(set(pert) - {"non-targeting"})[:60]
    keys = ["non-targeting"] + perts
    kidx = {k: i for i, k in enumerate(keys)}
    b = np.array([kidx.get(p, -1) for p in pert])
    X = f["X"]; S = np.zeros((len(keys), X.shape[1])); C = np.zeros(len(keys))
    for s in range(0, X.shape[0], 8192):
        blk = np.asarray(X[s:s + 8192], dtype=np.float64); bb = b[s:s + 8192]; ok = bb >= 0
        np.add.at(S, bb[ok], blk[ok]); np.add.at(C, bb[ok], 1)
D = S[1:] / C[1:, None] - S[0] / C[0]
cols = list(range(0, 2000, 4))  # 500 genes
genes = [gm[j] for j in cols]
deltas = D[:, cols]

res = M.verify_arc_archived_adata_pred_r16(ARC, deltas, perts, genes, hgnc_bridge=None)
print("correct names:", {k: res[k] for k in ("passed", "n_common_perts", "n_common_genes", "median_per_pert_pearson", "pooled_pearson") if k in res})
assert res["passed"] and res["n_common_genes"] == len(genes) and res["n_common_perts"] == len(perts)
assert res["median_per_pert_pearson"] > 0.999 and res["pooled_pearson"] > 0.999

shuf = list(np.random.default_rng(1).permutation(genes))
res2 = M.verify_arc_archived_adata_pred_r16(ARC, deltas, perts, shuf, hgnc_bridge=None)
print("shuffled names:", {k: res2[k] for k in ("passed", "n_common_genes", "median_per_pert_pearson", "pooled_pearson") if k in res2})
assert not res2["passed"]
print("R16_D5_TEST_PASS")
