# Sources and formulation differences

- Jari Brits, Niels Holtgrefe, Leo van Iersel and Samuel Martin, *On Tree-Network
  Distinguishability and Full Identifiability of Phylogenetic Networks*,
  [arXiv:2607.12919v3](https://arxiv.org/html/2607.12919v3), 25 August 2026.
  Section 2.3 gives the Fourier domain and uniform-root parametrization.
  Section 5.1 defines 2-sub-blob suppression to a single edge; the paragraph
  after Lemma 5.1 asks for exact equality in the root-trapping case. Its stated
  domain permits nonnegative transition entries and explicitly places all
  nontrivial Fourier eigenvalues in (0,1). Our example uses that explicit
  spectral reading for the rooted theorem. The quartet obstruction is proved
  even for the broader nonidentity positive-definite stochastic definition.
  Section 4 also uses continuous-time wording; that additional interpretation
  is not refuted by the example.
- Seth Sullivant, *Phylogenetic Network Models as Graphical Models*,
  [Bulletin of Mathematical Biology (2026)](https://doi.org/10.1007/s11538-026-01695-3),
  also [arXiv:2507.23056v2](https://arxiv.org/abs/2507.23056v2), 14 July 2026.
  Definition 3.10 specifies splittability. Proposition 3.12 uses strictly
  positive-entry open equivariant models. Propositions 5.6–5.7 retain two hidden
  sinks and replace their ancestral subgraph by a new root with two edges;
  Section 5 ends with the exact-equality question under a suitable transition
  model restriction. Root-law families are separate assumptions there.

The current primary versions were checked on 27 October 2026. They retain the
above open questions. The contribution proposed here is an explicit uniform-root
JC/K2P/K3P sink-joint classification with an exact realization, and a boundary
distinction for the different single-edge target. Existing inclusion, closure,
displayed-tree and equivariant-model machinery belongs to the cited sources.
No absolute-priority claim is made. The zero-entry example does not refute the
original splittable-model rooted statement.

This folder contains original arguments and checking code. No source paper,
third-party source package, private correspondence or author endorsement is
redistributed or implied.
