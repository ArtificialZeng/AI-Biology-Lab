"""Controller structure-only check of imru_full.h5ad against the R17 TeloHAEC reader spec.
Reads: HDF5 keys, shapes, dtypes, encoding attrs, obs/var column names, condition-column categories,
feature-type categories, number of categories per relevant column, and control-label flags in the target column
(only labels that are NOT plain gene symbols are printed). Never reads X / layer / raw values, never counts cells
per perturbation.
"""
import json
import re
import sys

import h5py

path = sys.argv[1]
NTC = ("no-target", "no_target", "notarget", "non-targeting", "nontargeting", "non_targeting", "non-target", "ntc", "ctrl_nt")
SAFE = ("safe_target", "safe-target", "safetarget", "safe_targeting", "safe-targeting", "safetargeting", "safe",
        "intergenic", "gene_desert", "ctrl_t")
UNASS = {"", "na", "n/a", "<na>", "nan", "none", "null", "unassigned", "multiplet", "doublet", "ambiguous", "unknown", "multiple"}
PLAIN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9.\-]*$")
COND = ["condition", "cell_treat", "treatment", "stim", "stimulation"]
TGT = ["gene_target", "target_gene", "perturbation", "gene", "target"]
GUIDE = ["guide", "sgRNA", "sgrna", "guide_id", "sgrna_id"]
FLASK = ["flask", "replicate", "batch", "gem_group", "channel", "lane"]
WELL = ["welltag", "well", "plate_well"]
NG = ["num_features", "n_guides", "guide_count", "n_sgrna"]


def pick(keys, pref):
    ks = set(keys)
    for c in pref:
        if c in ks:
            return c
    low = {k.lower(): k for k in sorted(keys)}
    for c in pref:
        if c.lower() in low:
            return low[c.lower()]
    return None


def cats(g, col):
    """Return category labels of a categorical column (modern or legacy), or None for non-categoricals."""
    if col in g and isinstance(g[col], h5py.Group) and "categories" in g[col]:
        return [x.decode() if isinstance(x, bytes) else str(x) for x in g[col]["categories"][()]]
    if "__categories" in g and col in g["__categories"]:
        return [x.decode() if isinstance(x, bytes) else str(x) for x in g["__categories"][col][()]]
    return None


def desc(node):
    if isinstance(node, h5py.Dataset):
        return {"kind": "dataset", "shape": list(node.shape), "dtype": str(node.dtype)}
    d = {"kind": "group", "encoding": node.attrs.get("encoding-type", b"").decode() if isinstance(node.attrs.get("encoding-type", ""), bytes) else str(node.attrs.get("encoding-type", "")),
         "keys": sorted(node.keys())[:12]}
    if "shape" in node.attrs:
        d["shape"] = [int(x) for x in node.attrs["shape"]]
    if "data" in node and isinstance(node["data"], h5py.Dataset):
        d["data_dtype"] = str(node["data"].dtype)
    return d


out = {}
with h5py.File(path, "r") as f:
    out["top_keys"] = sorted(f.keys())
    out["X"] = desc(f["X"])
    out["layers"] = {k: desc(f["layers"][k]) for k in f["layers"]} if "layers" in f else {}
    out["raw"] = {k: desc(f["raw"][k]) for k in f["raw"]} if "raw" in f else None
    obs, var = f["obs"], f["var"]
    obs_cols = [k for k in obs.keys() if k != "__categories"]
    var_cols = [k for k in var.keys() if k != "__categories"]
    out["obs_columns"] = sorted(obs_cols)
    out["var_columns"] = sorted(var_cols)
    res = {"condition": pick(obs_cols, COND), "gene_target": pick(obs_cols, TGT), "guide": pick(obs_cols, GUIDE),
           "flask": pick(obs_cols, FLASK), "welltag": pick(obs_cols, WELL), "n_guides": pick(obs_cols, NG)}
    out["reader_resolved_obs_columns"] = res
    c = res["condition"]
    out["condition_categories"] = cats(obs, c) if c else None
    for k in ("flask", "welltag", "guide"):
        col = res[k]
        cc = cats(obs, col) if col else None
        out[f"{k}_n_categories"] = None if cc is None else len(cc)
    t = res["gene_target"]
    tc = cats(obs, t) if t else None
    if tc is None:
        out["gene_target"] = "NOT_CATEGORICAL_OR_MISSING"
    else:
        low = [x.strip().lower() for x in tc]
        ntc = [x for x, l in zip(tc, low) if l in NTC or l.startswith(("no-target", "no_target", "non-targeting", "nontargeting", "non_targeting"))]
        safe = [x for x, l in zip(tc, low) if l in SAFE or l.startswith(("safe_target", "safe-target", "safetarget", "safe_")) or "gene_desert" in l or "intergenic" in l]
        unass = [x for x, l in zip(tc, low) if l in UNASS]
        multi = sum(1 for x in tc if re.search(r"[,;|+]", x))
        nonplain = [x for x in tc if not PLAIN.match(x.strip())]
        out["gene_target"] = {"n_categories": len(tc), "ntc_labels": ntc, "safe_labels": safe, "unassigned_labels": unass,
                              "n_multi_guide_labels": multi, "non_plain_symbol_labels_first30": nonplain[:30],
                              "n_non_plain_symbol_labels": len(nonplain)}
    for fc in ("feature_types", "feature_type", "modality", "type"):
        cc = cats(var, fc) if fc in var_cols else None
        if cc is not None:
            out["var_feature_type_col"] = {fc: cc}
            break
print(json.dumps(out, indent=1))
