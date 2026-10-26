# Source model and attribution

Floyd, Dinner, Murugan and Vaikuntanathan (2025), *Nature Communications*
16:7184, DOI [10.1038/s41467-025-61873-0](https://doi.org/10.1038/s41467-025-61873-0),
published Supplementary Information §X, equations (153)–(158), Supplementary
Figure 14B, defines a three-species model in which C catalyzes A→B and is not
consumed. With x=exp(F_a/2)>0 and all rates positive, it is

    dA/dt = a B + b C - c A - d A C,
    dB/dt = d A C + e x C - a B - f B/x,
    dC/dt = c A + f B/x - b C - e x C,
    A+B+C=1.

The rate map is a=W_AB, b=W_AC, c=W_CA, d=W_BA, e=W_BC, f=W_CB.
The displayed equations control this map; the surrounding narrative contains
a replacement-rate label inconsistent with them. The published plotted rates
are (a,b,c,d,e,f)=(1,1,0.1,10,0.05,4).

The published supplement already demonstrates a nonmonotone B response for this
example. The authors' [figure-code repository](https://github.com/csfloyd/NonEqExpressivity)
also computes the general symbolic self-consistent quadratic branches before
substituting the plotted rates. Neither the known nonmonotone example nor the
existence of an explicit quadratic solution is claimed as new here.

The present result classifies uniqueness and B criticality for every positive
rate vector in this exact model. The broad nonlinear-network critical-point
question and the separate linear structural-compatibility question are not
settled. The consumed-C sibling in Supplementary Figure 14C is outside scope.
No stochastic finite-copy, experimental, or global dynamical-attractor claim is
made.

The source was checked against the published supplement and the complete
author notebook, using a bounded current-primary-source comparison dated
2026-10-26. See literature.md for the searched scope and its limitations.
