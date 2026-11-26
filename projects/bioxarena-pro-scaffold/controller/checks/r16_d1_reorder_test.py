"""Controller test for DEVIATION_R16_D1 (cd4t_rest lane gene-order alignment in round16/src/run_eval_r16.py).

Compares the patched extract_ntc_bundle_from_parts_r16 with the frozen 848bdac version on synthetic parts:
  A. all parts share one gene order        -> patched output is bitwise identical to frozen output
  B. some parts list the same genes permuted -> patched output equals frozen output on the unpermuted parts
  C. a part with a different gene set       -> still raises ValueError
  D. a permuted part with a duplicated name -> still raises ValueError
Usage: round16/.venv/bin/python controller/checks/r16_d1_reorder_test.py
"""
import importlib.util
import subprocess
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np

P = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(P / "round16" / "src"))
import run_eval_r16 as NEW  # noqa: E402

frozen_src = subprocess.run(["git", "-C", str(P), "show", "848bdac:projects/bioxarena-pro-scaffold/round16/src/run_eval_r16.py"],
                            capture_output=True, text=True, check=True).stdout
tmpd = Path(tempfile.mkdtemp(prefix="r16_d1_"))
(tmpd / "run_eval_r16_frozen.py").write_text(frozen_src)
spec = importlib.util.spec_from_file_location("run_eval_r16_frozen", tmpd / "run_eval_r16_frozen.py")
OLD = importlib.util.module_from_spec(spec); sys.modules["run_eval_r16_frozen"] = OLD; spec.loader.exec_module(OLD)
assert "DEVIATION_R16_D1" not in frozen_src and "DEVIATION_R16_D1" in Path(NEW.__file__).read_text()

rng = np.random.default_rng(7)
G = 60
genes = [f"G{i:03d}" for i in range(G)]


def make_part(tag, X, names, batch):
    X = np.asarray(X, dtype=np.float32)
    return SimpleNamespace(tag=tag, var_gene_names=list(names), nt_mask=np.ones(X.shape[0], bool),
                           batches=np.full(X.shape[0], batch), cell_barcodes=[f"{tag}:{i}" for i in range(X.shape[0])],
                           read_rows_fn=lambda rows, X=X: X[np.asarray(rows)], read_chunk_fn=None)


mats = [rng.poisson(2.0, size=(260, G)) for _ in range(4)]
truth = [make_part(f"L{k}", mats[k], genes, f"b{k % 2}") for k in range(4)]
perm = rng.permutation(G)
mixed = [truth[0], truth[1]] + [make_part(f"L{k}", mats[k][:, perm], [genes[j] for j in perm], f"b{k % 2}") for k in (2, 3)]


def bundle(mod, parts):
    return mod.extract_ntc_bundle_from_parts_r16(parts, seed=20261601, resample_seed=20261602, include_resample=True)


def same(a, b):
    for k in ["raw_counts", "resample_raw_counts"]:
        if (a[k] is None) != (b[k] is None) or (a[k] is not None and not np.array_equal(a[k], b[k])):
            return False
    return all(a[k] == b[k] for k in ["gene_names", "batch_labels", "barcodes", "resample_barcodes", "resample_batch_labels"])


old_a, new_a = bundle(OLD, truth), bundle(NEW, truth)
assert same(old_a, new_a), "A: patched differs from frozen when orders already match"
assert new_a["ctrl_audit"]["deviation_r16_d1_reordered_part_tags"] == []
print("A ok: identical-order parts -> bitwise identical to frozen 848bdac")

try:
    bundle(OLD, mixed); raise AssertionError("frozen code should raise on mixed order")
except ValueError:
    pass
new_b = bundle(NEW, mixed)
assert same(old_a, new_b), "B: reordered parts do not reproduce the unpermuted reference"
assert new_b["ctrl_audit"]["deviation_r16_d1_reordered_part_tags"] == ["L2", "L3"]
print("B ok: permuted parts reordered by name -> equals frozen output on unpermuted parts; tags", new_b["ctrl_audit"]["deviation_r16_d1_reordered_part_tags"])

bad_set = [truth[0], make_part("Lx", mats[1], genes[:-1] + ["OTHER"], "b1")] + truth[2:]
try:
    bundle(NEW, bad_set); raise AssertionError("C: different gene set must raise")
except ValueError:
    print("C ok: different gene set still raises")

dup = [genes[j] for j in perm]; dup[5] = dup[6]
bad_dup = [truth[0], truth[1], make_part("Ld", mats[2][:, perm], dup, "b0"), truth[3]]
try:
    bundle(NEW, bad_dup); raise AssertionError("D: duplicated name must raise")
except ValueError:
    print("D ok: duplicated gene name in a reordered part still raises")
print("R16_D1_TEST_PASS")
