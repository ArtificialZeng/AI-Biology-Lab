# Third-party material and attribution

This directory contains original explanatory text, Python verification code, and computed JSON receipts prepared for AI Biology Lab. Original repository material is covered by the repository's [MIT license](../../LICENSE). No third-party paper, source archive, font, dataset, or copied upstream program is bundled here.

The mathematical model and established machinery are credited as follows:

- Haleš, Héliou, Maňuch, Ponty, and Stacho, *Combinatorial RNA design: designability and structure-approximating algorithm in Watson–Crick and Nussinov–Jacobson energy models*, Algorithmica 79 (2017), 835–856, [DOI](https://doi.org/10.1007/s00453-016-0196-x): separated coloring implies unique design; saturated degree criterion; obstruction motifs.
- Boury, Bulteau, and Ponty, *RNA inverse folding can be solved in linear time for structures without isolated stacks or base pairs*, [WABI 2024](https://doi.org/10.4230/LIPIcs.WABI.2024.19), [journal extension](https://doi.org/10.1186/s13015-025-00278-6): modulo-m machinery, fixed-modulus algorithm, all-long guarantee, and local helix transfers.
- Ashutosh S. Jogalekar, *Designability of RNA Targets with Up to Two Length-2 Helices*, [arXiv:2608.25194v1](https://arxiv.org/abs/2608.25194v1): F/Q local construction and at-most-two guarantee. The unlimited comparable-short statement here is a short corollary of that construction. The known three-root-branch failure is credited in that paper to the earlier Haleš example family. Its public [source repository](https://github.com/ajogalekar/rna-at-most-two-short-helices) licenses Lean source under Apache-2.0 and documentation/manuscript under CC BY 4.0; none of those source files is copied here.
- Jedwab, Petrie, and Simon, *An infinite class of unsaturated rooted trees corresponding to designable RNA secondary structures*, TCS 833 (2020), 147–163, [DOI](https://doi.org/10.1016/j.tcs.2020.05.046): prior unsaturated floral-tree class discussed in the literature comparison.
- Nussinov and Jacobson, [PNAS 77 (1980), 6309–6313](https://doi.org/10.1073/pnas.77.11.6309): the established secondary-structure dynamic-programming framework. The optimum/count code here is an original implementation using exact integers.

A short source-question quotation in `SOURCE.md` identifies the problem being refined. External links retain their publishers' and authors' terms. Python is used through its standard library; no vendored runtime or package is included.

AI assisted the derivation, implementation, literature search, and separate review. These activities do not establish publication priority or replace independent human expert review. The new statements and verification code have not been Lean formalized.
