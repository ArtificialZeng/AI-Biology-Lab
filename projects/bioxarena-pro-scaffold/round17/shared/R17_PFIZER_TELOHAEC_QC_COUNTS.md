# Round 17 Pfizer TeloHAEC (`Zenodo 18792213`, `imru_full.h5ad`) Counts-Only QC Gate

Model-free post-lock counts-only QC (`PREREG_J_LOCK.md` commit `9789b0b`; zero model predictions, zero correlation metrics).

## Summary Table

| Condition | Standardized Line | Role | Passes `PREREG_I` Floors? | `N_eval^R17` (`>=50`) | `N_eval_P5_R9` | `|G_eval^R9|` (`>=3000`) | `|G_eval^STATE|` | Half A NTC (`>=15`) | Half B NTC (`>=15`) | Split-Half Rule |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---|
| `Untreated` (basal, primary unit) | `pfizer_telohaec_basal` | `PASS_CONFIRMATORY_PRIMARY` | `True` | `114` | `131` | `6141` | `1843` | `2414` | `3087` | `flask_lexicographic_odd_even` |
| `Treated` (`IL-1b+TNF-a`, secondary within-family) | `pfizer_telohaec_cytokine` | `SECONDARY_WITHIN_FAMILY_REPLICATION` | `True` | `116` | `134` | `6141` | `1843` | `2740` | `3120` | `flask_lexicographic_odd_even` |

## Full Counts-Only JSON Audit

```json
{
  "untreated": {
    "line": "pfizer_telohaec_basal",
    "condition": "untreated",
    "role": "PASS_CONFIRMATORY_PRIMARY",
    "passes_prereg_i_floors": true,
    "n_eval_p5_r9": 131,
    "n_eval_r17": 114,
    "n_genes_g_eval_r9": 6141,
    "n_genes_g_eval_state": 1843,
    "n_kept_half_a_ntc": 2414,
    "n_kept_half_b_ntc": 3087,
    "per_half_ntc_ge15": true,
    "adapter_audit": {
      "source_zenodo_id": "18792213",
      "source_file": "imru_full.h5ad",
      "condition_selected": "untreated",
      "condition_cell_counts_raw": {
        "untreated": 428417,
        "treated": 435698,
        "other": 0
      },
      "x_slot_used": "X",
      "layout": "csr",
      "int_slots_found": [
        "X"
      ],
      "slot_selection_mode": "unique_candidate",
      "obs_columns_resolved": {
        "condition": "condition",
        "gene_target": "gene_target",
        "guide": "guide",
        "cellID": "cellID",
        "flask": "flask",
        "welltag": "welltag",
        "n_guides": null
      },
      "n_raw_obs_total": 864115,
      "n_var_total": 38606,
      "n_var_rna_kept": 38606,
      "n_var_non_rna_excluded": 0,
      "n_safe_targeting_excluded": 6073,
      "n_multiguide_or_unassigned_excluded": 0,
      "n_unapproved_symbol_excluded": 0,
      "n_rna_qc_excluded": 0,
      "min_umi": 300.0,
      "min_ngenes": 200,
      "max_pct_mt": null,
      "n_qc_kept_cells": 422344,
      "n_kept_ntc": 5501,
      "n_kept_targeting_singlets": 416843,
      "n_kept_half_a_ntc": 2414,
      "n_kept_half_b_ntc": 3087,
      "per_half_ntc_ge15": true,
      "half_rule_used": "flask_lexicographic_odd_even",
      "half_a_units": [
        "flask1"
      ],
      "half_b_units": [
        "flask2"
      ],
      "n_low_ntc_batch_fallback_cells": 0,
      "n_unique_pert_keys": 4626,
      "n_unique_canon_target_genes": 1733,
      "kept_barcodes_sha256": "d499cc8067811b7e2b66d21be6f4300a33c4cc9569460e9ee2e326a58c02625d"
    },
    "elapsed_s": 816.72
  },
  "treated": {
    "line": "pfizer_telohaec_cytokine",
    "condition": "treated",
    "role": "SECONDARY_WITHIN_FAMILY_REPLICATION",
    "passes_prereg_i_floors": true,
    "n_eval_p5_r9": 134,
    "n_eval_r17": 116,
    "n_genes_g_eval_r9": 6141,
    "n_genes_g_eval_state": 1843,
    "n_kept_half_a_ntc": 2740,
    "n_kept_half_b_ntc": 3120,
    "per_half_ntc_ge15": true,
    "adapter_audit": {
      "source_zenodo_id": "18792213",
      "source_file": "imru_full.h5ad",
      "condition_selected": "treated",
      "condition_cell_counts_raw": {
        "untreated": 428417,
        "treated": 435698,
        "other": 0
      },
      "x_slot_used": "X",
      "layout": "csr",
      "int_slots_found": [
        "X"
      ],
      "slot_selection_mode": "unique_candidate",
      "obs_columns_resolved": {
        "condition": "condition",
        "gene_target": "gene_target",
        "guide": "guide",
        "cellID": "cellID",
        "flask": "flask",
        "welltag": "welltag",
        "n_guides": null
      },
      "n_raw_obs_total": 864115,
      "n_var_total": 38606,
      "n_var_rna_kept": 38606,
      "n_var_non_rna_excluded": 0,
      "n_safe_targeting_excluded": 6154,
      "n_multiguide_or_unassigned_excluded": 0,
      "n_unapproved_symbol_excluded": 0,
      "n_rna_qc_excluded": 0,
      "min_umi": 300.0,
      "min_ngenes": 200,
      "max_pct_mt": null,
      "n_qc_kept_cells": 429544,
      "n_kept_ntc": 5860,
      "n_kept_targeting_singlets": 423684,
      "n_kept_half_a_ntc": 2740,
      "n_kept_half_b_ntc": 3120,
      "per_half_ntc_ge15": true,
      "half_rule_used": "flask_lexicographic_odd_even",
      "half_a_units": [
        "flask1"
      ],
      "half_b_units": [
        "flask2"
      ],
      "n_low_ntc_batch_fallback_cells": 0,
      "n_unique_pert_keys": 4487,
      "n_unique_canon_target_genes": 1732,
      "kept_barcodes_sha256": "238c3f95f172c70c9977641753a013cb903c2a01e150fa80c65dda83d88becb5"
    },
    "elapsed_s": 847.86
  }
}
```
