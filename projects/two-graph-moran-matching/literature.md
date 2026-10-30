# Primary-source comparison and provenance

Checked 2026-10-26 UTC (local date 2026-10-27, Asia/Shanghai), using read-only browsing. The model and nearest theorem are taken from the 2023 journal version. No external outreach occurred.

## Main source and nearest theorem

Tkadlec, Kaveh, Chatterjee and Nowak (2023), Evolutionary dynamics of mutants that modify population structure, Journal of the Royal Society Interface 20:20230355, DOI https://doi.org/10.1098/rsif.2023.0355. Full author-institution PDF: https://research-explorer.ista.ac.at/download/14657/14673/2023_RoyalInterface_Tkadlec.pdf . Frozen in evidence/source-2023.pdf.

Inspected §§2.1–2.3 for type-specific birth–death dynamics, equal fitness, connected undirected unweighted graphs and the random singleton initial law; §§3.1–3.5 for the finite N=6 enumeration and known graph classes; the first Discussion question for the original open inequality; and Appendix A, Theorem A.1, for the one-deleted-edge formula and hitting function. The source's N=6 numerical enumeration includes our graph at that one population size. The all-N two-disjoint-edge statement is not asserted in A.1, which deletes exactly one edge.

The original arXiv record https://arxiv.org/abs/2111.10890 has only v1, dated 21 November 2021, under the title Natural selection of mutants that modify population structure. Its version page was inspected; the 2023 journal text controls this project.

## Essential predecessor

Melissourgos, Nikoletseas, Raptopoulos and Spirakis, An extension of the Moran process using type-specific connection graphs, Journal of Computer and System Sciences 124 (2022), 77–96, DOI https://doi.org/10.1016/j.jcss.2021.07.007. Current arXiv PDF: https://arxiv.org/pdf/1710.07365v3 , dated 26 July 2021.

Inspected definitions on pp.3–4, the §3 result overview, Theorem 16 on p.18, Corollary 18 on p.22 and Lemma 12 on pp.15–16. This paper uses the reversed argument order f(resident,mutant;r). Theorem 16 supplies a general asymptotic lower bound, useful in the regime r>2c. It does not state our finite-N strict bound at neutral r=1. Corollary 18 concerns circulants with constant r>2. Lemma 12 is a weak graph-swap inequality, rather than a neutral 1/N bound. The generalized isothermal theorem cannot be applied merely because two different graphs are regular; the 2023 source explicitly demonstrates failure of that intuition.

## Current bounded overlap search

The abstract and version record of the author's 2026 paper Fixation location in structured populations, https://arxiv.org/abs/2605.26411v3 (28 June 2026), were inspected. That paper studies a different observable and is not evidence that every later paper has been cleared. The abstract of the vacant-space ecological model, DOI 10.64898/2026.10.16.759554, was inspected; its full text was unavailable.

Targeted queries included: “Moran two disjoint edges”, “fixation matching type-specific”, “fixation complete graph two edges”, the exact 2023 title with 2026/proof/conjecture, “two-graph complete fixation 2025/2026”, and “two regular graphs fixation”. Search hits served as leads; substantive support comes from the inspected primary texts above. No equivalent all-N two-disjoint-edge bound was found in those sources. This bounded negative finding does not establish universal priority.

## Author software relationship and license

The verified official author code/data deposit is https://doi.org/10.6084/m9.figshare.16910170 , with API metadata https://api.figshare.com/v2/articles/16910170 . It is version 1, modified 2023-10-30T15:14:44Z. The metadata, README and numerics-small-pairs.py are frozen in evidence/. The deposit metadata states CC BY 4.0: https://creativecommons.org/licenses/by/4.0/ . Credit for that archived code belongs to the paper's authors. Our mathematical checkers are original implementations.

The README advertises numerical pairwise tables for N=4 and N=5 and complete-versus-all-graphs values for N=6. In numerics-small-pairs.py, gg2mat builds the transition equations and gg2fp averages singleton fixation probabilities uniformly. Our author_reproduce.py executes that function prefix after mechanical Python 2 compatibility changes only: removing the unused graphs_iso import and replacing time.clock, xrange and <>. Its author SciPy outputs at N=4,6,8 with m=2 match our exact rational values within 1e-12. These floating-point reproductions are supplementary checks, not part of the all-N proof.

No author GitHub repository or appropriate issue/PR endpoint was verified through repository searches and the author's homepage https://sites.google.com/view/pepa-tkadlec/research . Figshare is a provenance deposit, not a pull-request host. Any public result still requires fresh review. No upstream target may be invented, and no author contact is authorized by this research pass.
