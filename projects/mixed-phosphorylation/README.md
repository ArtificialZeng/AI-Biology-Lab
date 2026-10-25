# Mixed-mechanism phosphorylation: an explicit finite-kinase instability bound

**Status:** exact partial mathematical result; numerical oscillation evidence; novelty and external review pending.

[Paper](manuscript/main.pdf) · [Proof](evidence/finite_kinase_proof.md) · [Certificate](evidence/finite_kinase_sign_certificate.json) · [Numerical checks](verification.md) · [Prior-work comparison](literature/finite_threshold_check.md)

## Statement

Use the mixed-mechanism mass-action model of Conradi, Mincheva, and Shiu, equations (1)–(3), with Table 1 rates `(1,1,1,1,100,1,9/10,3,1,100)`, total phosphatase `8`, and total substrate `40`.

For every **finite total kinase K ≥ 4660166/209569**, the unique positive equilibrium has exactly two eigenvalues with positive real part and four with negative real part on its six-dimensional compatibility class, counting algebraic multiplicity. In particular, this holds for **every K ≥ 23**.

The sufficient bound is not claimed sharp. This supplies explicit witnesses for the first assertion of the source's Conjecture 6.2. It does not prove the periodic-orbit assertion or the full conjecture.

## What is established and what is not

The exact proof uses a rational equilibrium parameter, monotonic coverage of an unbounded kinase ray, and positivity certificates for six Hurwitz determinants. Equilibrium uniqueness, the general parametrization/Hurwitz machinery, and the existence of oscillations at some parameters are prior work.

Numerical trajectories, periodic shooting, Floquet calculations, and first-order responses support attracting oscillations and concentration-dependent corrections for selected parameters. They do not prove periodic-orbit existence, global attraction, uniqueness, or invariance of a coordinate projection. A period change at one witness parameter pair does not disprove the source's full existential statement.

Related 2023/2025 conference work directly overlaps this topic. A bounded literature comparison found no equivalent complete proof in the inspected materials, but unavailable poster and follow-up texts prevent a priority claim. No external human referee or experimental validation is reported.

## Reproduction

From the repository root:

```sh
python -m pip install -r requirements-exact.txt
python scripts/reproduce.py verify
```

For rebuilding evidence, numerical runs, and the figure, see [the full instructions](../../docs/REPRODUCIBILITY.md).

The archive contains portable adaptations of the research code. Mathematical certificate files and the manuscript were preserved; operational provenance was replaced by [public integrity metadata](archive_manifest.json). Historical derivation notes in `theory.md` may describe intermediate or conditional results. This README and `claim.json` define the current public claim.
