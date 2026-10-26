# Frozen source record — 2026-10-26

Primary work: Saber Ahmed et al., *Identifiability of Directed-Cycle and Catenary Linear Compartmental Models*, SIAM Journal on Applied Dynamical Systems **25**(1), 304–350 (2026), DOI [10.1137/25M1728636](https://doi.org/10.1137/25M1728636).

Inspected primary versions:

- [Published full text](https://epubs.siam.org/eprint/IXIEMIECP3TISXXACWXU/full), Conjecture 4.6; the endpoint is present in the published version.
- [arXiv v3 HTML](https://arxiv.org/html/2412.05283v3#S4.SS1), dated 30 August 2025.
- [arXiv v3 PDF](https://arxiv.org/pdf/2412.05283v3), PDF page 24. A local research-input copy and pdftotext extraction are retained under ignored artifacts; these input files are not intended for redistribution.

Conjecture 4.6 has the following frozen mathematical reading. Write k_i=k_{i+1,i}, e_i=k_{0i}, d_i=k_i+e_i if i is a leak and d_i=k_i otherwise (indices cyclic). For every n>=3 directed cycle with input {1}, output {p}, 1<=p<=n, and generic local identifiability, for every

    a in {1,...,p-1} minus Leak,
    ell in Leak intersect {1,...,p},

it asserts

    {d_ell = k_a} is a subset of the coefficient-map singular locus.

The endpoint ell=p is part of the original quantified assertion. The proposed repair replacing its range by ell<p is a separate assertion, never silently substituted for the original.

Definition 4.1 defines the singular locus through failure of full rank of the coefficient-map Jacobian. Definition 3.1's exceptional family requires exactly two leaks. Definition 3.3 admits every model with at most one leak. Thus the model n=3, input{1}, output{2}, Leak={2} is within the conjecture's stated domain. Proposition 2.22(1) independently ensures its generic local identifiability.

Current-source search on 2026-10-26 used the DOI with correction/erratum terms, arXiv ID with conjecture number, and author/conjecture/singular-locus combinations. No directly matching correction or prior counterexample was found in those searches. This is a bounded search report, not a guarantee of priority or confirmation by the authors.
