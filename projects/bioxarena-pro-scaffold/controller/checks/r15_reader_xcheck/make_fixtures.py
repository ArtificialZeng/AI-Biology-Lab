"""Controller-made synthetic VCC-style h5ad fixtures for the R15 reader cross-check.

Written with h5py in AnnData on-disk format (no anndata dependency). Neither A-T2 nor A-T3 has seen these files.
Usage: python3 -I make_fixtures.py <common_genes.txt> <targets_pool.txt> <out_dir>
"""
import json
import sys
from pathlib import Path

import h5py
import numpy as np

STR = h5py.string_dtype("utf-8")


def w_str(grp, name, values):
    ds = grp.create_dataset(name, data=np.asarray(list(values), dtype=object), dtype=STR)
    ds.attrs["encoding-type"] = "string-array"
    ds.attrs["encoding-version"] = "0.2.0"


def w_cat(grp, name, values):
    cats = sorted(set(values))
    idx = {c: i for i, c in enumerate(cats)}
    g = grp.create_group(name)
    g.attrs["encoding-type"] = "categorical"
    g.attrs["encoding-version"] = "0.2.0"
    g.attrs["ordered"] = False
    g.create_dataset("codes", data=np.asarray([idx[v] for v in values], dtype=np.int16))
    w_str(g, "categories", cats)


def w_df(f, name, index, cols):
    g = f.create_group(name)
    g.attrs["encoding-type"] = "dataframe"
    g.attrs["encoding-version"] = "0.2.0"
    g.attrs["_index"] = "_index"
    g.attrs.create("column-order", np.asarray(list(cols.keys()), dtype=object), dtype=STR)
    w_str(g, "_index", index)
    for k, v in cols.items():
        w_cat(g, k, v)


def w_dict(f, name):
    g = f.create_group(name)
    g.attrs["encoding-type"] = "dict"
    g.attrs["encoding-version"] = "0.1.0"


def write_h5ad(path, X, obs_index, obs_cols, genes):
    X = X.tocsr()
    with h5py.File(path, "w") as f:
        f.attrs["encoding-type"] = "anndata"
        f.attrs["encoding-version"] = "0.1.0"
        g = f.create_group("X")
        g.attrs["encoding-type"] = "csr_matrix"
        g.attrs["encoding-version"] = "0.1.0"
        g.attrs["shape"] = np.asarray(X.shape, dtype=np.int64)
        g.create_dataset("data", data=X.data.astype(np.float32))
        g.create_dataset("indices", data=X.indices.astype(np.int32))
        g.create_dataset("indptr", data=X.indptr.astype(np.int64))
        w_df(f, "obs", obs_index, obs_cols)
        vg = f.create_group("var")
        vg.attrs["encoding-type"] = "dataframe"
        vg.attrs["encoding-version"] = "0.2.0"
        vg.attrs["_index"] = "_index"
        vg.attrs.create("column-order", np.asarray(["gene_id"], dtype=object), dtype=STR)
        w_str(vg, "_index", genes)
        w_str(vg, "gene_id", [f"ENSG{9000000000 + i:011d}" for i in range(len(genes))])
        for k in ("layers", "obsm", "varm", "obsp", "varp", "uns"):
            w_dict(f, k)


def simulate_line(rng, genes, targets, line_seed, n_per_target=60, n_ntc=400, n_batches=4):
    """Raw integer counts for one line: NTC + targets (self-knockdown + random program) + junk cells."""
    import scipy.sparse as sp

    lr = np.random.default_rng(line_seed)
    G = len(genes)
    gidx = {g: i for i, g in enumerate(genes)}
    base = np.exp(lr.normal(-1.0, 1.4, G))
    base = base / base.sum()
    line_mult = np.exp(lr.normal(0.0, 0.4, G))
    mu0 = base * line_mult
    mu0 = mu0 / mu0.sum()
    rows, labels, guides, batches = [], [], [], []

    def draw(mu, n, depth=4000):
        size = np.exp(rng.normal(0.0, 0.3, n)) * depth
        return rng.poisson(np.outer(size, mu))

    def add(mat, lab, guide):
        for r in mat:
            rows.append(sp.csr_matrix(r))
            labels.append(lab)
            guides.append(guide)
            batches.append(f"Flex_1_{int(rng.integers(1, n_batches + 1)):02d}")

    add(draw(mu0, n_ntc), "non-targeting", "non-targeting_00001|non-targeting_00002")
    for t in targets:
        mu = mu0.copy()
        mu[gidx[t]] *= 0.15
        prog = lr.choice(G, 120, replace=False)
        mu[prog] *= np.exp(lr.normal(0.0, 0.8, 120))
        mu = mu / mu.sum()
        add(draw(mu, n_per_target), t, f"{t}_1|{t}_2")
    add(draw(mu0, 12, depth=150), "non-targeting", "non-targeting_00001|non-targeting_00002")  # low-UMI junk
    add(draw(mu0, 6), f"{targets[0]},{targets[1]}", "multi")  # multi-target junk
    X = sp.vstack(rows).tocsr()
    return X, labels, guides, batches


def main():
    common = [g for g in Path(sys.argv[1]).read_text().split() if g]
    pool = [g for g in Path(sys.argv[2]).read_text().split() if g]
    out = Path(sys.argv[3])
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(777)
    targets = pool[:30]
    rest = sorted(set(common) - set(targets))
    genes = sorted(set(targets) | set(rng.choice(rest, 3800 - len(targets), replace=False).tolist()))

    lines = {"Cal-33": 101, "HeLa": 202, "Jurkat": 303}
    sims = {name: simulate_line(rng, genes, targets, s) for name, s in lines.items()}

    import scipy.sparse as sp

    def stack(names, with_line_col, with_ctx, with_batch, prefix):
        Xs, lab, gd, bt, ln, ctx = [], [], [], [], [], []
        ctx_code = {"Cal-33": "A", "HeLa": "B", "Jurkat": "C"}
        for nm in names:
            X, l, g, b = sims[nm]
            Xs.append(X)
            lab += l
            gd += g
            bt += b
            ln += [nm] * len(l)
            ctx += [ctx_code[nm]] * len(l)
        X = sp.vstack(Xs).tocsr()
        idx = [f"{prefix}{i:06d}-1" for i in range(X.shape[0])]
        cols = {"target_gene": lab, "guide_id": gd}
        if with_batch:
            cols["batch"] = bt
        if with_line_col:
            cols["cell_line"] = ln
        if with_ctx:
            cols["context"] = ctx
            cols["ntc_id"] = [g if l == "non-targeting" else "" for l, g in zip(lab, gd)]
        return X, idx, cols

    # case1: one mixed file, explicit cell_line column, batch present (Tier 1), includes Jurkat
    d = out / "case1_mixed_explicit"
    d.mkdir(exist_ok=True)
    X, idx, cols = stack(["Cal-33", "HeLa", "Jurkat"], True, False, True, "AAAC")
    write_h5ad(d / "adata_Validation.h5ad", X, idx, cols, genes)

    # case2: one mixed file, context only, no batch (Tier 2), official map given separately
    d = out / "case2_context_only"
    d.mkdir(exist_ok=True)
    X, idx, cols = stack(["Cal-33", "HeLa", "Jurkat"], False, True, False, "CCGT")
    write_h5ad(d / "adata_Validation.h5ad", X, idx, cols, genes)
    (d / "context_map.json").write_text(json.dumps({"A": "Cal-33", "B": "HeLa", "C": "Jurkat"}, indent=1))

    # case3: one file per line, no line column, batch present; manifest per Addendum A rule (c)
    d = out / "case3_per_line_files"
    d.mkdir(exist_ok=True)
    man = {}
    for nm, fn in (("Cal-33", "cal33_challenge.h5ad"), ("HeLa", "hela_challenge.h5ad")):
        X, idx, cols = stack([nm], False, False, True, "GGTA")
        write_h5ad(d / fn, X, idx, cols, genes)
        man[fn] = {"line": nm, "url": "https://example.invalid/synthetic-manifest", "quote": f"{fn} contains {nm} cells only (synthetic)"}
    (d / "line_manifest.json").write_text(json.dumps(man, indent=1))

    # case4: non-integer X (log-normalized) -> must be EXCLUDED_FEASIBILITY
    d = out / "case4_nonint"
    d.mkdir(exist_ok=True)
    X, idx, cols = stack(["Cal-33"], True, False, True, "TTTA")
    X = X.astype(np.float64)
    X.data = np.log1p(X.data * 1.37)
    write_h5ad(d / "adata_Validation.h5ad", X, idx, cols, genes)

    print(json.dumps({"n_genes": len(genes), "targets": targets, "cells_per_line": int(sims["Cal-33"][0].shape[0])}))


if __name__ == "__main__":
    main()
