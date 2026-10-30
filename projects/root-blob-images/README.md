# Exact images of a hidden phylogenetic root subgraph

An original proof-and-code packet for technical review, using the uniform
stationary JC, K2P and K3P displayed-tree models. External human confirmation
and historical-priority confirmation are pending.

The proposed main result classifies the joint distributions at the two sink
exits of any finite hidden ancestral two-sink DAG: they are exactly the
strictly positive-entry kernels of the selected model with all nontrivial
Fourier eigenvalues in (0,1). A constructive compensating kernel realizes each
target with strictly positive inheritance probabilities and no identity edges.
Consequently replacing that ancestral subgraph by a new root and two edges
preserves its full descendant distribution image in this specified setting.

A companion rational quartet example shows why compressing the subgraph to a
single edge is a different operation on a domain that allows zero transition
entries. In the literal strict-spectrum/nonnegative-entry K2P/K3P domain of
Brits et al. v3, an identifiable internal edge can have a zero entry which the
original root-trapping diamond cannot reproduce, even after every original
edge parameter is retuned. A support-cardinality argument also covers the
broader positive-definite nonidentity stochastic reading. The leaf distribution
itself is strictly positive.

These statements address different source formulations. The example does not
refute Sullivant's splittable-model hypotheses, does not use continuously
embeddable Kimura parameters, and does not settle arbitrary root laws, general
transition families, or the arbitrary-graph JC single-edge question.

[Proof](proof.md) · [exact scope](SCOPE.md) · [primary sources](SOURCE.md) ·
[graph argument](graph-proof.md) · [verification scope](VERIFICATION.md).

## Reproduce without external dependencies

Python 3.12 or 3.13, standard library only:

```sh
python -B projects/root-blob-images/verify.py verify
python -B -O projects/root-blob-images/verify.py verify
```

The replay checks the complete file manifest, all 256 rational target leaf
probabilities, observable internal-edge reconstruction, nine explicit
constructions against direct ancestral-variable summation, and the separate
graph certificate. It writes only to stdout; temporary process directories
are outside this project. Optimization of the wrapper does not remove its
checks. Scientific child processes explicitly run with optimization disabled.

The general graph theorem and global nonmembership depend on the analytical
proof, not the finite certificate. Graph enumeration covers five simple root
positions and twenty displayed paths; its stated domain does not include all
possible hidden expansions. The cut/direction proof handles root trapping.

AI assisted the derivation, implementation and separate-context review. These
are not external human peer review, proof-assistant formalization, journal
acceptance or trusted-runner certification. This packet makes no empirical
claim about errors in inferred biological histories. All code and notes in
this folder are original; primary papers are linked rather than redistributed.
