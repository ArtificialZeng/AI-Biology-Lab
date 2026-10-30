# Frozen source question and contribution

**Source question.** Schreiber, *Coexistence and extinction in flow-kick systems: An invasion growth rate approach*, Journal of Mathematical Biology 91:50 (2025), DOI [10.1007/s00285-025-02283-0](https://doi.org/10.1007/s00285-025-02283-0), Sections 2.2, 3.3 and 5.2, equation (26).

The model has two species, four strictly positive constant growth rates rho_i(y) for y=0,1, common dilution 0<delta<1, exponential growth until total density one and alternating environments at each dilution. Immediately after dilution, x_1+x_2=delta. Starting from any x_1,x_2>0 on that state space and either phase, does

    [rho_1(0)/rho_2(0)+rho_1(1)/rho_2(1)]/2<1

force x_1(n)->0, x_2(n)->delta at successive post-transfer times?

The source computes invasion exponents and establishes a locally attracting resident-only cycle. Its current journal text explicitly leaves the full positive basin unresolved. The proposed proof answers this precise strict-inequality question affirmatively, with an explicit odds bound. This is not a challenge to the source's invasion or permanence theorems.

**Coherent extension.** The same one-transfer argument supplies complete positive-composition dynamics for any finite prescribed periodic sequence of strictly positive growth pairs, including neutral/equality cases, and a cumulative envelope for arbitrary prescribed sequences. The extension does not replace or narrow the original two-phase target.

**Evidence categories.** `proof.md` contains the quantified analytic argument. Symbolic identities and rational witness checks are exact algebra. The numerical certificate contains finite floating-point observations only. Neither simulation nor CI establishes historical priority, independent human review or the universal theorem.
