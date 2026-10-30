# LWL85 twofold-site weighting audit

This is a bounded scientific-software correction, not a new biological discovery. Both Biopython's modern `calculate_dn_ds` and legacy `cal_dn_ds` use the fourfold count in two numerator terms that require the twofold count. It reproduces in released 1.88 and source commit `893f1543f87981927104fb7af9de43c4579dfe60`.

The independent oracle counts standard-code sites and single-nucleotide substitutions without Biopython helpers. It deliberately rejects threefold sites, cross-degeneracy changes and multi-difference codons. Official MEGA documentation and PAML's explicitly labelled original LWL85 formula agree on the weighting; PAML distinguishes LWL85m and LPB93. The original 1985 full equation pages were unavailable; no PAML or MEGA executable was run.

For a 25-codon synthetic alignment with L0=50, L2=20, L4=5, two twofold transitions and one zero-fold transversion, old dN/dS is 0.335276186644 versus corrected 0.083819046661. This fourfold ratio is specific to this fixture. A second fixture exercises the dN weight; an equal-L2/L4 control conceals the error. No effect on actual genes, selection findings or published papers has been established.

## Reproduce

Use Python 3.12 or 3.13 and an environment outside this project:

```sh
python -m pip install -r requirements.txt
python -B verify.py verify
python -O -B verify.py verify
```

Verification hashes all packet inputs, tests the installed pinned release, copies its package into a temporary directory, applies only the two coefficient corrections in each API, and compares both API directions with independent counts/formulas. It does not edit the installed package or archived evidence. Optimized wrappers explicitly run scientific children with optimization disabled. `compute` prints a newly recomputed certificate; it does not overwrite any file. Floating values are compared within stated tolerances and serialized to 12 decimal places; site counts and substitution fractions use exact rational arithmetic.

The upstream patch also adds three regression tests and refreshes existing dependent LWL85 expectations, including one tutorial output exposed by the full offline suite. See [statement](statement.md), [sources](sources.md), [review status](review-summary.md) and [third-party rights](THIRD_PARTY.md). AI assistance was substantial. External human and maintainer confirmation remain pending.
