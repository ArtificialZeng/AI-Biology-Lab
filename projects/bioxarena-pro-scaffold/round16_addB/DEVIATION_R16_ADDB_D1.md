# DEVIATION_R16_ADDB_D1 — Occam inference needs torch-geometric, which Addendum B v2 §6.1 wrongly called unnecessary

Written by Claude (master controller), 2026-11-27 05:15 +0800.

## What happened
- The main line built `round16_addB/.venv` on the primary target: Python 3.8.19 with the §6.1 pins, frozen in `R16_ADDB_VENV_FREEZE.txt`.
- Importing Occam's unchanged `occam.vcell.pipeline` from the hash-pinned source tree fails with `ModuleNotFoundError: No module named 'torch_geometric'` (main line `R16_ADDB_DEVIATION`, 05:10).
- `torch_geometric` is imported unconditionally on the inference path:
  - `src/occam/vcell/gnn.py:13` (`GATv2Conv, GCNConv, GINConv, GATConv`), reached from `vcell/model.py:21` and `vcell/pipeline.py:49`;
  - `src/occam/vcell/adjacency_map.py:12` (`to_undirected`, `add_self_loops`), called at runtime in `ConcatKnowledgeModel.empty_like_state`.
- Addendum B v2 §6.1 says "`torch-geometric` is unnecessary at inference for `adjacency_map`". That is wrong. Occam's `pyproject.toml:48` pins the extra `graph = ["torch-geometric==2.6.1", "pyarrow>=14.0"]` for exactly these code paths.
- The main line's AST scan finds no other missing third-party import on the inference path (`pyarrow` is not imported there; `wandb` is lazy inside `train()`), and no import of the 7 data members the controller removed.
- No workaround, patch or mock was attempted.

## Ruling
- Install `torch-geometric==2.6.1` into `round16_addB/.venv` only. Wheel `torch_geometric-2.6.1-py3-none-any.whl`: 1,135,632 B, PyPI SHA-256 `8faeb353f9655f7dbec44c5e0b44c721773bdfb279994da96b9b8b12fd30f427`, requires Python >= 3.8.
- Install it hash-pinned (requirements line with `--hash=sha256:8faeb353…`). Its declared runtime dependencies (aiohttp, fsspec, jinja2, numpy, psutil, pyparsing, requests, tqdm) may be resolved normally.
  - They must not change any §6.1 pin. If the resolver would change numpy, scipy, torch, pandas, anndata, h5py or scikit-learn, stop and report.
- Do not install `pyarrow`, and do not install the optional compiled extensions (torch-scatter, torch-sparse, pyg-lib).
- Update `R16_ADDB_VENV_FREEZE.txt` with the full resulting list.
- Scope: the total added download is a few MB of PyPI packages into the isolated venv. It is a declared dependency of the human-approved Occam code (approval item (c), interpreter and packages for `round16_addB/.venv`), and is far below the per-item approval threshold. The controller reports it to the Human project lead.
- Same rule for both parties: A-T3's environment for Addendum B uses the identical frozen venv specification.

## What does not change
- Every Addendum B hypothesis, metric, gate, threshold and arm.
- The Occam source tree, which stays unchanged.
- The MPS ban until A-T3 has finished R16.
