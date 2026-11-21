# `ArcInstitute/cell-eval2` Local Reproduction & VCC 2025/2026 Verification Guide

This document describes how to reproduce [`ArcInstitute/cell-eval2`](https://github.com/ArcInstitute/cell-eval2) locally on macOS (CPU) and records the exact verification runs executed against public VCC 2025 data (`H1-VCC-2025-training.h5ad`).

---

## 1. Local Environment & Installation (CPU / macOS)

### Verified Versions
- Repository / Version: `ArcInstitute/cell-eval2` **`v0.22.0`** (MIT License; Python `>=3.11`)
- Local source snapshot tested: `/tmp/round9_S4/repos/cell-eval2/cell-eval2-main`
- Isolated venv tested: `/tmp/round9_A_T4/venv` (Python `3.12.13`)
- Core installed package versions observed in venv:
  - `cell-eval2==0.22.0`
  - `pdex==0.3.0` (CPU parallel differential-expression backend)
  - `anndata==0.13.4`
  - `polars==2.1.0`
  - `scanpy==1.12.6`
  - `scipy==1.18.1`
  - `numpy==2.5.5`
  - `zstandard==0.25.0`
  - `pytest==9.2.0`
  - `torch`: not installed (base `cell-eval2` install is CPU-only; no CUDA or PyTorch needed when pinning `device=cpu` and `de.backend=pdex`)

### Reproducing the Clean Virtual Environment
```bash
# Option A: From local source checkout (editable install)
uv venv /tmp/round9_A_T4/venv --python 3.12
uv pip install --python /tmp/round9_A_T4/venv/bin/python \
  -e /tmp/round9_S4/repos/cell-eval2/cell-eval2-main \
  pdex==0.3.0 zstandard pytest

# Option B: Clean install directly from PyPI (no local clone needed)
uv venv /tmp/round9_A_T4/venv --python 3.12
uv pip install --python /tmp/round9_A_T4/venv/bin/python \
  "cell-eval2==0.22.0" "pdex==0.3.0" zstandard pytest

# Verify CLI
/tmp/round9_A_T4/venv/bin/cell-eval2 --help
```

> **macOS Gatekeeper / XProtect note**: The very first `import scipy` / `import anndata` / `import cell_eval2` on freshly unpacked wheels may take ~30–50 seconds while macOS scans new shared libraries; subsequent imports finish in ~1 second.

### Reproducing Upstream Preset & Golden-Score Tests
```bash
cd /tmp/round9_S4/repos/cell-eval2/cell-eval2-main && \
/tmp/round9_A_T4/venv/bin/pytest tests/test_vcc2026_preset.py tests/test_golden_score.py -v
```
- Observed result: **13 passed in 9.93s** (including exact golden `avg_score = 1.5659362305778723`).

---

## 2. Preset & Metric Differences (`vcc2026` vs `v1` / `v2` / VCC 2025)

Measured live from `EvalConfig.from_preset(...)` in `cell-eval2==0.22.0`:

| Preset / Profile | `metrics` | `version` | `control_source` | `input_type` | `target_sum` | `discrimination.distance` | `de.fdr_scope` |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **`vcc2026`** | `vcc2026` | `v2` | `real` | `counts` | `1000000.0` (CPM) | `cosine` | `per_pert` |
| `v2` | `full` | `v2` | `real` | `counts` | `1000000.0` | `cosine` | `per_pert` |
| `v1` *(deprecated)* | `full` | `v1` | `pred` | `lognorm` | `None` | `l1` | `global` |
| `cell-eval-0.7.6` *(deprecated)* | `full` | `v1` | `pred` | `lognorm` | `10000.0` | `l1` | `per_pert` |
| VCC 2025 (`--profile vcc`) | `vcc` (`de_wilcoxon_overlap`, `expr_mae`, `pds_l1`) | `v2` | `real` *(by default)* | passes `--input-type counts` | `1000000.0` | `cosine` | `per_pert` |

### Key Practical Requirements for `vcc2026` Local Runs
1. **Pin CPU flags on macOS**: Always pass `--set device=cpu --set de.backend=pdex` to avoid GPU/CUDA auto-detection ambiguities.
2. **Inputs must be sparse non-negative integer counts** (`--input-type counts` implied by `--preset vcc2026`, `allow_fractional_counts=False`, `allow_discrete=False`, `autodetect_input_type=False`).
3. **6 Competition-Scored Metrics in `vcc2026`**:
   - `de_wilcoxon_lfc_nmae`
   - `expr_mse_unbiased_capped_norm`
   - `de_wilcoxon_direction_fidelity_yield_raw`
   - `de_wilcoxon_direction_reach_raw`
   - `de_wilcoxon_sig_jaccard`
   - `pds_cosine`
   - Combined into `avg_score` when scored against a `prep-real-bundle` (`from_replicate` column scales baseline = 0, 5-split replicate anchor = 1). In addition, `agg_results.csv` emits 4 expression diagnostics (`expr_distance_unbiased`, `expr_mse_unbiased`, `expr_mse_unbiased_capped`, `expr_real_mass_ratio`).

---

## 3. Crucial Format Bridge: Upload `.h5ad` vs `cell-eval2` Scoring `.h5ad`

There is a strict two-stage distinction between the **challenge upload matrix** and **`cell-eval2` scoring input**:
1. **Upload `.h5ad` (validated by `vcc prep`)**:
   - **MUST NOT** include `non-targeting` control cells (`reject_controls=True`).
   - Must contain `obs["context"]`, `obs["target_gene"]`, exact gene ordering (`18,533` genes for real VCC 2026), exact cell counts (`400` cells/pert/context), and non-negative integer CSR counts.
2. **Scoring `.h5ad` passed to `cell-eval2 run -ap ... -ar ...`**:
   - `cell_eval2.io.validate_pair` requires `pred` and `real` to share the **exact same set of perturbation labels** — **including `non-targeting`**.
   - Even though `vcc2026` uses `control_source="real"` (replacing predicted control rows with `real` control cells during scoring), `pred` must still carry at least one valid `non-targeting` control block per scored context.
3. **Bridge helper**:
   - `vcc2026_format.attach_real_controls_for_celleval2()` extracts one context (when `context` is present) from an upload `.h5ad`, splices the real `non-targeting` cells from the reference `.h5ad`, aligns `var_names`, and writes a CSR count `.h5ad` ready for `cell-eval2 run`.

---

## 4. Exact `cell-eval2` End-to-End Scoring CLI Sequence

For competition-style scoring (`from_replicate` + `avg_score`) on a local reference dataset (`REAL_H5AD`), the exact verified flags (in `cell-eval2==0.22.0`) are:

```bash
CE2=/tmp/round9_A_T4/venv/bin/cell-eval2
FLAGS="--preset vcc2026 --pert-col target_gene --control non-targeting --set device=cpu --set de.backend=pdex"

# 1. Build 0-end baseline AND save its dispersed count matrix (--save-pred is required for prep-real-bundle!)
$CE2 baseline -ar "$REAL_H5AD" -o out/baseline --save-pred out/baseline/baseline_pred.h5ad $FLAGS

# 2. Run submission prediction (--adata-pred must already be bridged with non-targeting cells)
$CE2 run -ap "$PRED_FOR_CE2" -ar "$REAL_H5AD" -o out/run_pred $FLAGS

# 3. Prep real bundle (--baseline expects the .h5ad from --save-pred, and bundle name flag is --id)
$CE2 prep-real-bundle \
  --real "$REAL_H5AD" \
  --baseline out/baseline/baseline_pred.h5ad \
  -o out/real_bundle \
  --id my-verify-bundle \
  --force \
  $FLAGS

# 4. Score submission against the real bundle (emits scored.csv with from_replicate & avg_score)
$CE2 score \
  --user-agg out/run_pred/agg_results.csv \
  --real-bundle out/real_bundle \
  -o out/run_pred/scored.csv
```

> **Pitfalls discovered & verified**:
> - `prep-real-bundle` uses `--id` (NOT `--bundle-id`).
> - `prep-real-bundle --baseline` expects the `.h5ad` prediction matrix saved via `baseline --save-pred <path>` (NOT the `baseline/` output directory).

---

## 5. Observed Verification Results on Public VCC 2025 Data

Script: `run_celleval2_vcc2025_verify.py`  
Public input: `/tmp/round9_S4/repos/cell-eval2/cell-eval2-main/docs/data/H1-VCC-2025-training.h5ad` (`600` cells × `1,000` genes; 5 perturbations `MED12`, `SRC`, `STAT1`, `TET1`, `TMSB4X` at 100 cells each + `100` `non-targeting` controls).

Command executed:
```bash
/tmp/round9_A_T4/venv/bin/python \
  /Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2/vcc_linespecific/src/prospective/run_celleval2_vcc2025_verify.py
```

### Verification Arms
- **Signal arm** (`pred_signal_upload.h5ad`): Perturbed-to-control mean fold ratios estimated from `H1-VCC-2025-training.h5ad` and emitted as dispersed sparse integer counts via `emit_dispersed_counts_for_context()`, validated with upload validator (`passed=True`, `500` cells × `1,000` genes, `nnz=499,799`, `min_cell_total=9,198`, `max_cell_total=44,565`), then bridged with real `non-targeting` controls.
- **Null arm** (`pred_null_upload.h5ad`): Unit fold ratios (`1.0`, predicting no perturbation shift), validated (`passed=True`, `500` cells × `1,000` genes, `nnz=499,790`, `min_cell_total=9,197`, `max_cell_total=44,394`), then bridged.

### Observed `vcc2026` Scores (`--real-bundle`, `vcc2025-verify-bundle`)

| Metric | Signal (`from_baseline`) | Signal (`from_replicate`) | Null (`from_baseline`) | Null (`from_replicate`) |
| --- | ---: | ---: | ---: | ---: |
| `de_wilcoxon_lfc_nmae` | `0.8407` | `1.6825` | `-0.3774` | `-0.7553` |
| `expr_mse_unbiased_capped_norm` | `1.0000` | `1.0000` | `0.0000` | `0.0000` |
| `de_wilcoxon_direction_fidelity_yield_raw` | `1.0000` | `1.2298` | `0.0000` | `0.0000` |
| `de_wilcoxon_direction_reach_raw` | `1.0000` | `1.0000` | `-0.2969` | `-0.2969` |
| `de_wilcoxon_sig_jaccard` | `-0.0473` | `-0.1129` | `0.0000` | `0.0000` |
| `pds_cosine` | `1.0000` | `1.4063` | `-0.3333` | `-0.4688` |
| **`avg_score`** | `NaN` | **`1.0343`** | `NaN` | **`-0.2535`** |

*(Exact values in `/tmp/round9_A_T4/celleval2_vcc2025_verify/celleval2_vcc2025_verify.json`: `signal_avg_score_from_replicate = 1.034265584954988`, `null_avg_score_from_replicate = -0.2534912253149939`. Note: `from_baseline` for `avg_score` is `NaN` by design when scoring against a `--real-bundle`, because `from_replicate` carries the competition `avg_score`. Note also that on this small `1,000`-gene toy fixture, `de_wilcoxon_lfc_nmae` passed the `>=10` gated DE gene threshold on 1 of 5 perturbations and dispersed control resampling changed per-gene variance relative to the real replicate anchor on `de_wilcoxon_sig_jaccard`, while every other metric and the overall `avg_score` clearly separated signal from null.)*

### Observed VCC 2025 Legacy Profile (`--profile vcc --input-type counts`) on Signal Arm
- `de_wilcoxon_overlap`: mean `0.2952` (median `0.3333`, range `0.0000`–`0.6429`)
- `expr_mae`: mean `0.02568` (median `0.02552`, range `0.02407`–`0.02714`)
- `pds_l1`: mean `1.0000` (median `1.0000`, all 5 perturbations)

### Prospective Unit & Integration Tests (`test_prospective.py`)
```bash
/tmp/round9_A_T4/venv/bin/pytest -v \
  /Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2/vcc_linespecific/src/prospective/test_prospective.py
```
- Observed result: **6 passed in 6.25s** (format emission & per-cell cap enforcement, CSV parsers, 8 validation rejection modes, `.vcc` tar+zstd archive structure, full 3-context bundle + `SHA256SUMS.txt` + `cell-eval2` scoring bridge, and `p2/shared/lines/` loader + `PREREG_COMMIT.txt` gate enforcement).
