"""Controller (no model): gene -> SE-600M global-index map built by arc-state 0.11.1 FilteredGenesCounts for
(a) main line: var['gene_name'] = raw, var_names = raw, var_names_make_unique(), gene_column='gene_name'
(b) A-T3:      var_names = raw (no make_unique), gene_column=None
on the controller dry-run control gene lists. Addendum A §11 item 3 requires no duplicate global IDs.
Usage: round16/.venv/bin/python controller/checks/r16_d7_mapping_test.py <dryrun_dir>
"""
import sys
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import torch
from omegaconf import OmegaConf
from state.emb.data.loader import FilteredGenesCounts

P = Path(__file__).resolve().parents[2]
SE = P / "round16/shared/models/SE-600M"
cfg = OmegaConf.load(SE / "config.yaml")
pe = torch.load(SE / "protein_embeddings.pt", weights_only=False)
dd = Path(sys.argv[1])
for ln in ["a549_ko", "cd4t_rest", "dld1", "feng_ipsc", "nadig_hepg2", "nadig_jurkat"]:
    genes = [str(g) for g in np.load(dd / f"main__{ln}.npz")["genes"]]
    X = np.zeros((2, len(genes)), np.float32)
    a = ad.AnnData(X=X.copy()); a.var["gene_name"] = genes; a.var_names = genes; a.var_names_make_unique()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = ad.AnnData(X=X.copy()); b.var_names = genes
    mm = FilteredGenesCounts(cfg, adata=a, adata_name="m", shape_dict={"m": a.shape}, protein_embeds=pe, gene_column="gene_name").ds_emb_map["m"]
    ma = FilteredGenesCounts(cfg, adata=b, adata_name="t", shape_dict={"t": b.shape}, protein_embeds=pe, gene_column=None).ds_emb_map["t"]
    def dupg(m):
        v = m[m != -1]; u, c = np.unique(v, return_counts=True); return int((c > 1).sum())
    diff_cols = np.flatnonzero(mm != ma)
    print(f"{ln}: mapped main={int((mm != -1).sum())} at3={int((ma != -1).sum())} | duplicate global IDs main={dupg(mm)} at3={dupg(ma)} | "
          f"columns mapped differently={len(diff_cols)} {[genes[j] for j in diff_cols][:12]}")
