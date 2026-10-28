# Sources and relation to prior work

- Currie, Englander, Esparza-Lozano, Gross, Hill, Long, Olds, O'Connor,
  Ranasinghe and Sum, *Semialgebraic Conditions for Identifying Triangles in
  Phylogenetic Networks*, [arXiv:2606.26673v1](https://arxiv.org/html/2606.26673v1)
  (2026). Eq. (6), Fig. 3, Theorem 3.1 and its leaf permutations supply the
  native JC parameterization and complete strict homogeneous image criterion.
  Appendix B supplies the Fourier interpretation. These are used as published
  mathematical dependencies, not claimed as new. The source theorem is not
  challenged. Author code:
  [max-hill/semialgebraic-conditions-for-triangle-networks](https://github.com/max-hill/semialgebraic-conditions-for-triangle-networks),
  inspected at `8bcd201e2a5efe072789c3039960ea61b32661b4`.
- Matsen, Mossel and Steel, *Mixed-up Trees: the Structure of Phylogenetic
  Mixtures*, [arXiv:0705.4328v2](https://arxiv.org/html/0705.4328v2).
  Corollary 29, Theorem 30 and the following even-state/JC remark give important
  scalar-rate no-mimicking precedents for tree topology. They are not a triangle-
  orientation image criterion. No new general covariance or mixture method is
  claimed here.
- Frankel and Ané, *Summary Tests of Introgression Are Highly Sensitive to Rate
  Variation Across Lineages*, [doi:10.1093/sysbio/syad056](https://doi.org/10.1093/sysbio/syad056).
  This work includes scalar rates across genes and gamma site variation as well
  as lineage effects. Its inference-error results are important context; the
  present exact compatibility statement does not establish robustness of those
  tests or of a misspecified likelihood estimator.

The inspected primary sources support a modest source-specific extension for
fixed-inheritance JC triangles. The bounded overlap check is not proof of
priority: complete text clearance of Englander et al. (2025) was unavailable;
available primary descriptions place its main network results in a triangle-
free class. No theorem from that unavailable text is used in the proof.

All proof prose and Fraction verification code in this directory are original
contributions with substantial AI assistance. No third-party paper, figure,
private correspondence or author-repository source code is redistributed.
Separate AI contexts checked the proof and reconstructed the original-state
probability calculations. These checks are not external human peer review or
formal proof-assistant certification. Maintainer and independent human review
remain pending.
