# Sources and contract

Somya Mehra and Peter G. Taylor, *A multitype Markovian branching process with one-type population size dependence*, [arXiv:2609.12567v1](https://arxiv.org/html/2609.12567v1), 11 September 2026. The [version history](https://arxiv.org/abs/2609.12567) listed v1 on 27 October 2026.

The source defines the jump law in §2 and assumptions A/B/C. Theorem 3.1 asserts almost-sure extinction and all polynomial moments under A/B. Remark 3.1 asks about an exponential extinction-time moment. Its general model permits type-1 reproduction; the immune example in §7 imposes additional structure.

Our two examples target the literal §2 assumptions: the nonsingularity assertion following (7) and the transience assertion before (12) require extra conditions. They do not contradict a formulation that explicitly imposes those conditions, and are not counterexamples to the specific §7 biological model. The corrected theorem in proof.md explicitly states all its added hypotheses.

Classical tools: Down, Meyn and Tweedie (1995), *Exponential and uniform ergodicity of Markov processes*, Annals of Probability 23, 1671–1691, [DOI](https://doi.org/10.1214/aop/1176987798), [author-hosted PDF](https://faculty.eng.ufl.edu/meyn/wp-content/uploads/sites/671/archive/spm_files/Papers_pdf/dmt.pdf). General drift/return-time theory does not verify the model-specific low-controller strip. The proof here supplies that drift directly, then gives its finite-set absorption argument explicitly. Birth–death crossing identities and Hurwitz–Metzler Lyapunov vectors are standard methods.

Limited exact-title/identifier and related primary-source searches found no equivalent source-specific resolution on this date; this is not a priority certification. No third-party paper or private correspondence is included.

## First-moment extension

The current primary version/history and Remark 3.1 were checked again on 27 October 2026; v1 remains the listed version and the exponential-moment direction remains stated there. The updated theorem in first-moment-proof.md explicitly keeps a Hurwitz limiting matrix and zero accessibility, replaces bounded support by a uniform first offspring moment, and proves a corrected A/B corollary without C. It does not establish the literal A/B formulation. The earlier cutoff function is retained; expectation bounds, nonexplosion/localization and renewal are the strengthening and classical justification, rather than new general Markov-process theory.

The inspected primary sources did not establish an equivalent source-specific first-moment result. This limited comparison does not certify priority or justify a separate publication. The author's attributable public repositories did not include code for this paper; no unrelated upstream repository is targeted. Private correspondence is not included.
