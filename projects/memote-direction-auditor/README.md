# MEMOTE reaction-direction audit

**A two-line candidate repair restores reaction-orientation invariance in a published metabolic-model quality check.** This is a software result; no defect in a living organism or an actual genome-scale model is asserted.

For a stoichiometric column S, replacing it by -S while changing bounds [l,u] to [-u,-l] describes the same feasible fluxes. MEMOTE 0.17.0's `close_boundaries_sensibly` instead changes every non-reversible internal reaction to [0,1], including negative-only [-1000,0]. Thus equivalent representations can receive opposite ATP energy-cycle diagnoses.

Using real COBRApy models and the complete unmodified MEMOTE ATP detector, with five elemental/charge-balanced metabolites and no external input, we reproduce both false-positive 0→1 and false-negative 1→0 outcomes. The candidate preserves the allowed negative direction and repairs both cases. Positive/reversible/boundary normalization and the documented [0,0] reopening policy are retained. Integer flux witnesses and the equation `sigma*v+d=0` provide a direct exact oracle; GLPK is a separate computational check.

| Allowed reaction | Native maximum ATP dissipation | Original detector | Candidate |
| --- | ---: | ---: | ---: |
| Reverse-written hydrolysis, negative-only | 0 | 1 | 0 |
| Reverse-written charging, negative-only | 1 | 0 | 1 |

## Reproduce

Use Python 3.11 in an isolated environment, install `requirements.txt`, then run from the repository root:

```sh
python projects/memote-direction-auditor/verify.py
```

This runs original and candidate packages in separate temporary processes, checks pinned source bytes, and exercises the full ATP entry point, equivalent orientations, fixed negative fluxes and boundary controls. It does not change installed packages. Outputs go to ignored `artifacts/`; archived evidence remains unchanged.

[Candidate patch](patches/preserve-negative-direction.patch) · [Recorded integration evidence](evidence/direction_patch_verification.json) · [Exact toy evidence](evidence/exact_toy_witnesses.json) · [Licenses](THIRD_PARTY.md)

## Limits and prior work

The existing helper intentionally removes fixed magnitude constraints; reopening [0,0] is not presented as a newly discovered bug. The narrow issue is reversing an allowed negative direction. No full MEMOTE test-suite/CLI report, genome-scale impact or maintenance acceptance is claimed. The tested latest Yeast-GEM closed-boundary ATP LP gave zero; it is not an example of a new energy cycle. Earlier MEMOTE issues #218/#224 concern a different reversible dead-end issue; #570 already covers the compartment-coverage idea. Searches were not exhaustive. No upstream issue/PR or other message has been sent for this candidate.

Sources: [fixed helper](https://github.com/opencobra/memote/blob/bf76477b4087feefc8ad429ebd492e30ef10eebc/src/memote/support/helpers.py#L676-L703), [complete ATP detector](https://github.com/opencobra/memote/blob/bf76477b4087feefc8ad429ebd492e30ef10eebc/src/memote/support/consistency.py), [original MEMOTE paper](https://doi.org/10.1038/s41587-020-0446-y). Prepared with AI assistance; no external human review is claimed.
