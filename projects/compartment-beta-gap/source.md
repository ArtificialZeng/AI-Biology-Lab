# Sources, exact model and prior overlap

Checked on 2026-10-26. Primary sources:

1. David F. Anderson, Aidan S. Howells and Diego Rojas La Luz, *Stochastic Reaction Networks Within Interacting Compartments with Content-Dependent Fragmentation*, Bulletin of Mathematical Biology 88:106 (2026), [DOI 10.1007/s11538-026-01659-7](https://doi.org/10.1007/s11538-026-01659-7), current equation (7), Section 2.2, Open Problem 3.8, and Propositions 3.9–3.10. [Authors' arXiv v4](https://arxiv.org/abs/2511.10223v4). Proposition 3.10 already uses an inverse-log function for an auxiliary aggregate CRN; the candidate contribution is the actual outside-carrier correction, occupation bound and transfer to the intermediate mixture, not inverse-log functions in general. Use the current version: the aggregate-CRN explosivity question left open in older versions is already resolved in Proposition 3.10.
2. Byungwoong Yoo, *Fragmentation perturbations and explosion in interacting compartments*, public preprint v1.1.0, 2026-09-13, [DOI 10.5281/zenodo.22728357](https://doi.org/10.5281/zenodo.22728357). The official latest-version endpoint checked on 2026-10-26 still reported v1.1.0. Section 9 defines the exact polynomial-decay mixture and explicitly leaves `0 < beta < 1` unclassified by its arguments. The precise intermediate-gap attribution is to the primary TeX in the Zenodo v1.1.0 package, lines 789–813. Submission acknowledgement is distinct from journal acceptance.

The original state is a finite multiset of pairs `(e,s)`. Every compartment gains an enzyme at rate 1 and gains substrate at rate `es`; it fragments at rate `s`; every unordered pair merges at rate 2. There is no external compartment arrival/deletion. Fragmentation preserves each species. The original model permits arbitrary admissible allocation kernels, so its diagram alone does not specify one lifetime outcome.

For fixed `c,p` in `(0,1)`, the follow-up defines `psi_beta=(1-c/(1+s)^beta)psi_b+[c/(1+s)^beta]kappa_p`, where `psi_b` produces unordered daughters `{0,x}` and `kappa_p` allocates every molecule independently with bias `p`.

Credited prior outcomes are nonexplosion at beta=0, almost-sure explosion at beta>=1 for positive substrate, and zero-substrate nonexplosion. The candidate here concerns only positive-probability explosion for `10/11 < beta < 1`. No conclusion for `0 < beta <= 10/11`, almost-sure upgrade, necessary-and-sufficient kernel condition or absolute novelty claim is made.

No verified public source-code GitHub upstream was found for the cited source papers; the directly overlapping follow-up package is available on Zenodo. An upstream code PR target is therefore unverified. Third-party TeX, PDFs and code are linked rather than redistributed in this package.
