# Primary sources and exact contribution boundary

Bounded screen checked28October2026UTC; it does not certify exhaustive priority.

* Brewster etal, The Mixed Birth-death/death-Birth Moran Process, arXiv:2511.18252v3 (13 January2026), ITCS2026, DOI10.4230/LIPIcs.ITCS.2026.29. https://arxiv.org/html/2511.18252v3
  Text after Lemma2.4 expresses a belief in general fitness/initial-set monotonicity; RemarkA.1 explains the normalization/clock difficulty. This is a historical source statement; no current universal open status is asserted.
* Brewster etal, Mixed updating in structured populations, PLOS Computational Biology22(9):e1014829 (30 September2026). https://doi.org/10.1371/journal.pcbi.1014829
  This is a distinct later manuscript, not the same arXiv file's revision. The native fixed-mixture model is in main sections2.1–2.2. Its dB probability delta equals1-lambda here. Final S1Appendix section8 Lemma4 compares fitness1 with fitness>=1 from the same initial set, not arbitrary ordered positive fitnesses. Star/cycle formulas and neutral bidegreed results are known source content.
  Final supplement: https://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1014829.s001&type=supplementary
* Author code: https://github.com/davidb2/mixed-updating/tree/01d9266a44f0bcf568671b92ff639ab032460983
  The original full-state implementation supplies a source-alignment reference. This directory's code is an independent rational reimplementation, with no source-file import/copy.
* Monk and van Schaik(2021), Martingales and the characteristic functions of absorption time on bipartite graphs, Royal Society Open Science8:210657. https://doi.org/10.1098/rsos.210657
  Complete-bipartite count reduction for pure Bd and its product-martingale analysis are prior work. Its Discussion notes that uniform-death-first updating breaks that product-martingale symmetry. Count lumping itself is not claimed as new.
* Massey(1987), Stochastic Orderings for Markov Processes on Partially Ordered Spaces, Mathematics of Operations Research12(2):350–367. https://doi.org/10.1287/moor.12.2.350
  Stochastic comparison/basic coupling is established machinery. Only metadata/abstract were inspected in this bounded screen; no numbered theorem is imported. The proof gives the required pair-generator construction directly.

The proposed delta combines all-positive-fitness/all-fixed-mixture fixation comparison on complete bipartite and balanced multipartite graphs with the sharp all-fitness classification for the native coordinatewise count generator. It is not a general-graph resolution, a new pure-Bd formula, or a claim of nonmonotone fixation on excluded graphs.

The analytic model has positive integer part sizes, at least two parts, no loops/weights/directed edges/mutation, constant resident fitness1 and mutant fitness r>0, and an independent constant choice of update rule. The exact oracle accepts rational parameters; its state guard is computational only. Count-process comparisons concern native L=P-I in common physical time, not an imposed shared attempt schedule.
