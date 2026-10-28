# Cyclic-pair necessity for three-strain full-PSD communities

A proposed proof/code note for the mixed-inhibition model of Athreya, Gokhale and Verma (2025), DOI10.1371/journal.pcbi.1013763. Prepared by Shawn Zeng (Zijian Zeng), AI Biology Lab, UCSI University.

For three strains, suppose every antibiotic channel has exactly one producer P, sensitive strain S and degrader D. With any finite number of channels, common positive inhibition/degradation strengths and positive constant effective growth rates, a locally asymptotically stable interior fixed point requires two distinct channel columns related by a cyclic rotation.

The proof reduces all possible column supports to a fixed-role classification. Without the required pair, two immune strains have a constant/incompatible abundance ratio, or the exact log-ratio derivative forces a multiplier greater than1. This is a necessary condition, not sufficiency.

The source conjectures a broader cyclicity necessity after numerical topology/parameter searches. This note supplies an analytic obstruction only for the full-PSD three-strain class. The general conjecture, other phenotype columns, intrinsic resistance, more strains, cycling, permanence and spatial settings remain outside the result.

## Read and reproduce

See proof.md for the self-contained argument and sources.md for the source contract and nearest prior comparison.

~~~sh
python -m pip install -r requirements.txt
python verify.py
python -O verify.py
~~~

The verifier reads this archive without rewriting it. It checks:

* All63 nonempty supports of the six PSD column types, including the15 no-cyclic-pair supports. The analytic orbit argument, not a finite prefix of channel counts, proves completeness for arbitrary finite channel multiplicities.
* Exact symbolic native-survival derivative identities for common P/D and the exact simplex-chart determinant.
* Two exact ordered-cost equilibrium constructions using rational exponential values, plus their floating native-map Jacobians.
* Sixty frozen floating-point native-map fixtures across all15 supports. Their positive effective rates are chosen to place a fixed point; they are not all claimed to arise from the original ordered-cost parameterization.
* An exact pairwise-invasion comparison showing why the inspected2016necessary bistable-pair condition alone does not exclude this class.
* Archive membership and hashes. These protect reproducibility, not mathematical truth.

The infinite-parameter theorem is the written analytic proof. Symbolic identities and finite numerical checks corroborate it; CI is not a proof of the theorem. Floating Jacobians use complex-step differentiation of a direct normalized native map and tolerance-based comparisons.

## Status and contribution

This is a proposed modest proof/code contribution. A separate fresh AI context reconstructed the native equations, derivatives, invariant and coordinate argument, and independently chose supporting fixtures. AI assistance was substantial in development, source discovery and checking. This is not external human peer review or formal proof-assistant certification. Author validation and priority confirmation are pending. No experimental, clinical, published-biological-result or major-discovery claim is made.

Original prose and independently written code in this directory are MIT-licensed (LICENSE). The primary paper and original author repository retain their own licenses. No third-party source code, paper PDF or private correspondence is included.
